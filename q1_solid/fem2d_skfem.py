"""Q1：实心圆铜导线 2D 涡流场 FEM —— gmsh 网格 + meshio 桥 + scikit-fem 装配。

（NGSolve wheel 的 gmsh 导入损坏、netgen CSG2D 同心圆死循环 → 本路线为生产路线。）

模型（MQS 相量，RMS，e^{jωt}，复对称弱式——不共轭试探函数）：
    E_z = -(jωA_z + v)，J_z = σE_z
    (1) (1/μ)S·a - jωσ·M·a - σ·b·v = 0        （A 方程）
    (2)  jωσ·bᵀ·a + σ·A_cu·v   = -I           （总流约束）
    分块复对称线性方程组直接解（scipy spsolve）。
    R_ac = P/I²，P = σ∫_cu|jωa+v|²；交叉校验 Z = -v/I。
DC 极限自检：v = -I/(σA_cu) → E=+I/A 均匀、Z→R_dc（精确）。
"""
import math
import sys
import time

import numpy as np
import meshio
from scipy.sparse import bmat, csr_matrix
from scipy.sparse.linalg import spsolve
from skfem import MeshTri, Basis, BilinearForm, LinearForm
from skfem.element import ElementTriP2
from skfem.helpers import grad, dot

sys.path.insert(0, ".")
from sim.constants import SIGMA, RHO, MU0, MU_R, F0, I_RMS, skin_depth
from sim.gmsh_mesh import make_wire_mesh


def load_mesh(path, a=1e-3):
    """meshio → skfem.MeshTri；材料按几何分类（质心在铜盘内=1）。

    不用 gmsh physical 标签的原因：skfem 构网可能重排单元，标签会错位
    （DC 自检实测：总面积对、空间分布错 → 虚部爆炸）。保角网格下
    「质心在圆内」分类是精确的，且直接推广到 Q2/Q3 的多丝束。"""
    mm = meshio.read(path)
    tris = np.vstack([cb.data for cb in mm.cells if cb.type == "triangle"])
    mesh = MeshTri(mm.points[:, :2].T, tris.T)
    centers = mesh.p[:, mesh.t].mean(axis=1)
    ind_cu_elem = (np.hypot(centers[0], centers[1]) <= a).astype(float)
    return mesh, ind_cu_elem


def qp_elem_array(basis, elem_vals):
    """单元数组 (ne,) → 积分点数组 (ne, nqp)（skfem 数组参数约定）。"""
    nqp = basis.W.size
    return np.broadcast_to(elem_vals[:, None], (basis.nelems, nqp))


@BilinearForm
def _lap(u, v, w):
    return dot(grad(u), grad(v))


@BilinearForm
def _mass(u, v, w):
    return w.sig * u * v


@LinearForm
def _load(v, w):
    return w.sig * v


def solve(mesh, ind_cu, f=F0, i_rms=I_RMS, a=1e-3):
    w_omega = 2 * math.pi * f
    basis = Basis(mesh, ElementTriP2(), intorder=4)

    S = BilinearForm(lambda u, v, w: dot(grad(u), grad(v))).assemble(basis)
    sig_qp = qp_elem_array(basis, ind_cu)
    M_cu = BilinearForm(lambda u, v, w: w.sig * u * v).assemble(basis, sig=sig_qp)
    b_cu = LinearForm(lambda v, w: w.sig * v).assemble(basis, sig=sig_qp)
    area_cu = b_cu.sum()

    r_air = np.hypot(*mesh.p).max()
    n = basis.N
    K11 = S / (MU0 * MU_R) - 1j * w_omega * SIGMA * M_cu
    K12 = (-SIGMA * b_cu).reshape(-1, 1)
    K21 = (1j * w_omega * SIGMA * b_cu).reshape(1, -1)
    K22 = csr_matrix(SIGMA * area_cu).reshape(1, 1)
    K = bmat([[K11.tocsr(), csr_matrix(K12)], [csr_matrix(K21), K22]], format="csr")

    rhs = np.zeros(n + 1, dtype=complex)
    rhs[-1] = -i_rms

    # 显式 Dirichlet 集：外边界节点 dof + 外边界边 dof
    # （get_dofs(callable) 会绑到 2D 不存在的 edges 参数上 → 拿到 0 个 dof 的坑）
    bn = mesh.boundary_nodes()
    bf = mesh.boundary_facets()
    D = np.concatenate([basis.nodal_dofs[:, bn].flatten(),
                        basis.facet_dofs[:, bf].flatten()]).astype(int)
    keep = np.ones(n + 1, dtype=bool)
    keep[D] = False
    keep[-1] = True  # v 未知量不固定

    Kff = K[np.ix_(keep, keep)]
    x = np.zeros(n + 1, dtype=complex)
    x[keep] = spsolve(Kff.tocsr(), rhs[keep])

    a_vec, v = x[:-1], x[-1]
    # P = σ∫|jωa+v|² = σ[ω² aᴴMa + 2Re(jω·conj(v)·bᵀa) + |v|²A]
    p = SIGMA * (w_omega**2 * np.real(a_vec.conj() @ (M_cu @ a_vec))
                 + 2 * np.real(1j * w_omega * np.conj(v) * (b_cu @ a_vec))
                 + abs(v) ** 2 * area_cu)
    z = -v / i_rms  # 端口每米阻抗
    rdc = RHO / (math.pi * a * a)
    return dict(basis=basis, a=a_vec, v=v, p=p, rac=p / i_rms**2, z=z,
                rdc=rdc, ratio=p / i_rms**2 / rdc, area_cu=area_cu, dof=n)


def sample_jr(basis, a_vec, v, f=F0, a=1e-3, npts=400):
    """沿 +x 半径采样 |J(r)|（P2 插值，skfem 12.x 单参数调用；失败显式报错）。"""
    w_omega = 2 * math.pi * f
    rs = np.linspace(0.0, a, npts)
    pts = np.vstack([rs, np.zeros_like(rs)])  # skfem 要 (dim, npts)
    vals = np.asarray(basis.interpolator(a_vec)(pts), dtype=complex)
    jr = np.abs(SIGMA * (-(1j * w_omega * vals + v)))
    assert jr.max() > 1e3, f"J 采样异常 max={jr.max():.3e}"
    jr[0] = jr[1]
    return rs, jr

def exact_ratio(a=1e-3, f=F0):
    from scipy.special import jv
    k = (1 - 1j) / skin_depth(f)
    x = k * a
    z = k * RHO * jv(0, x) / (2 * math.pi * a * jv(1, x))
    return z.real / (RHO / (math.pi * a * a))


if __name__ == "__main__":
    a, f = 1e-3, F0
    tgt = exact_ratio(a, f)
    print(f"δ={skin_depth(f)*1e6:.2f}µm  精确 Bessel 目标 Rac/Rdc = {tgt:.5f}\n")
    print(f"{'refine':>7} {'ndof':>7} {'R_ac/R_dc':>11} {'误差%':>8} {'P[W/m]':>8} {'Z口径':>9} {'耗时s':>6}")
    for refine in (1, 2, 4):
        t0 = time.time()
        path = make_wire_mesh(a=a, f=f, refine=refine)
        mesh, ind_cu = load_mesh(path, a=a)
        res = solve(mesh, ind_cu, f=f, a=a)
        dt = time.time() - t0
        err = (res["ratio"] / tgt - 1) * 100
        zrac = res["z"].real
        print(f"{refine:7d} {res['dof']:7d} {res['ratio']:11.4f} {err:8.2f} "
              f"{res['p']:8.3f} {zrac*1e3:7.3f}m {dt:6.1f}")
        if refine == 2:
            rs, jr = sample_jr(res["basis"], res["a"], res["v"])
            np.savetxt("data/q1_fem_J_r.csv", np.column_stack([rs, jr]),
                       delimiter=",", header="r_m,J_amplitude", comments="")
            print("  -> data/q1_fem_J_r.csv（refine=2）")
    print("\n两口径校验：P/I² 与 Re(Z) 应一致（误差即数值残差）")
