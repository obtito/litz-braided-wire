"""Q2：N=127 丝 Litz 束的 2D 涡流仿真（并联物理）。

关键公式结构（与 Q1 同构，仅铜区域变为丝的并集）：
    并联丝束两端短接 → 所有丝共享同一每米压降 v（单全局未知量）；
    E_z = -(jωA_z + v)；J_z = σE_z；约束 Σ I_i = 20 A。
    → 逐丝电流 I_i = -σ(jω∫_{s_i}A dA + v·A_si) 自然涌现，
      外层丝偏载（邻近效应）无需任何人为设定。

输出：R_ac/R_dc（束）、逐丝电流表（η_I、按环统计）、与等截面实心基准对比。
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
from geo.strand_packing import hex_packing


def load_bundle(path, centers, strand_r):
    """meshio → skfem MeshTri + 铜指示 + 逐丝指示矩阵 (nelems, N)。"""
    mm = meshio.read(path)
    tris = np.vstack([cb.data for cb in mm.cells if cb.type == "triangle"])
    p_arr = np.asarray(mm.points)[:, :2]
    mesh = MeshTri(p_arr.T, tris.T)
    cent = p_arr[tris].mean(axis=1)          # (m, 2)
    d = np.hypot(cent[:, None, 0] - centers[None, :, 0],
                 cent[:, None, 1] - centers[None, :, 1])
    inside = d <= strand_r                    # (m, N)
    ind_cu = inside.any(axis=1).astype(float)
    strand_ind = inside.astype(float)
    return mesh, ind_cu, strand_ind


def solve_bundle(mesh, ind_cu, f=F0, i_rms=I_RMS):
    w = 2 * math.pi * f
    basis = Basis(mesh, ElementTriP2(), intorder=4)
    nqp = basis.W.size

    def qp(v):
        return np.broadcast_to(v[:, None], (basis.nelems, nqp))

    S = BilinearForm(lambda u, v, w_: dot(grad(u), grad(v))).assemble(basis)
    M_cu = BilinearForm(lambda u, v, w_: w_.sig * u * v).assemble(basis, sig=qp(ind_cu))
    b_cu = LinearForm(lambda v, w_: w_.sig * v).assemble(basis, sig=qp(ind_cu))
    area_cu = b_cu.sum()

    n = basis.N
    K11 = S / (MU0 * MU_R) - 1j * w * SIGMA * M_cu
    K12 = csr_matrix((-SIGMA * b_cu).reshape(-1, 1))
    K21 = csr_matrix((1j * w * SIGMA * b_cu).reshape(1, -1))
    K22 = csr_matrix(SIGMA * area_cu).reshape(1, 1)
    K = bmat([[K11.tocsr(), K12], [K21, K22]], format="csr")
    rhs = np.zeros(n + 1, dtype=complex)
    rhs[-1] = -i_rms

    bn = mesh.boundary_nodes()
    bf = mesh.boundary_facets()
    D = np.concatenate([basis.nodal_dofs[:, bn].flatten(),
                        basis.facet_dofs[:, bf].flatten()]).astype(int)
    keep = np.ones(n + 1, dtype=bool)
    keep[D] = False
    keep[-1] = True
    x = np.zeros(n + 1, dtype=complex)
    x[keep] = spsolve(K[np.ix_(keep, keep)].tocsr(), rhs[keep])

    a_vec, v = x[:-1], x[-1]
    p = SIGMA * (w**2 * np.real(a_vec.conj() @ (M_cu @ a_vec))
                 + 2 * np.real(1j * w * np.conj(v) * (b_cu @ a_vec))
                 + abs(v) ** 2 * area_cu)
    return dict(basis=basis, a=a_vec, v=v, p=p, area_cu=area_cu, dof=n)


def strand_currents(basis, strand_ind, a_vec, v, f=F0):
    """逐丝复电流 I_i（RMS 相量）。"""
    w = 2 * math.pi * f
    nqp = basis.W.size
    lf = LinearForm(lambda vv, w_: w_.sig * vv)
    currents, areas = [], []
    for i in range(strand_ind.shape[1]):
        b_i = lf.assemble(basis, sig=np.broadcast_to(strand_ind[:, i][:, None],
                                                     (basis.nelems, nqp)))
        A_i = b_i.sum()
        currents.append(-SIGMA * (1j * w * (b_i @ a_vec) + v * A_i))
        areas.append(A_i)
    return np.array(currents), np.array(areas)


if __name__ == "__main__":
    from sim.gmsh_mesh import make_bundle_mesh

    N = 127
    r_s = math.sqrt(5.0e-6 / (N * math.pi))
    centers, R_bundle = hex_packing(N, r_s)
    rdc = RHO / (N * math.pi * r_s**2)

    t0 = time.time()
    path = make_bundle_mesh(centers, r_s)
    mesh, ind_cu, strand_ind = load_bundle(path, centers, r_s)
    res = solve_bundle(mesh, ind_cu)
    Ii, Ai = strand_currents(res["basis"], strand_ind, res["a"], res["v"])
    dt = time.time() - t0
    print(f"[{dt:.1f}s] dof={res['dof']}  ΣI = {Ii.sum().real:.3f} A（校验=20）")
    print(f"A_cu = {res['area_cu']*1e6:.4f} mm²   R_dc = {rdc*1e3:.3f} mΩ/m")
    rac = res["p"] / I_RMS**2
    print(f"R_ac = {rac*1e3:.3f} mΩ/m   R_ac/R_dc = {rac/rdc:.4f}")
    z = -res["v"] / I_RMS
    print(f"Re[v]/I 口径 = {z.real*1e3:.3f} mΩ/m（应与 P/I² 一致）")
    imag = np.abs(Ii)
    imean = I_RMS / N
    eta = (imag.max() - imean) / imean
    print(f"逐丝 |I|：mean={imean:.4f}  max={imag.max():.4f}  min={imag.min():.4f}"
          f"  η_I = {eta*100:.1f}%")
    # 按环统计（环 k 丝数 6k，环半径 ≈ k·spacing）
    spacing = 2 * r_s * 1.02
    ring = np.round(np.hypot(centers[:, 0], centers[:, 1]) / spacing).astype(int)
    print("\n环 | 丝数 | 平均|I| [A] | 相对均值")
    for k in sorted(set(ring)):
        m = ring == k
        print(f"{k:2d} | {m.sum():4d} | {imag[m].mean():.4f} | {imag[m].mean()/imean:.3f}")
    np.savetxt("data/q2_strand_currents.csv",
               np.column_stack([np.arange(N), ring, centers[:, 0], centers[:, 1],
                                imag, np.angle(Ii, deg=True)]),
               delimiter=",",
               header="strand,ring,x_m,y_m,I_amplitude_A,I_phase_deg", comments="")
    import json
    json.dump(dict(N=N, d_mm=2*r_s*1e3, bundle_od_mm=2*R_bundle*1e3,
                   a_cu_mm2=res["area_cu"]*1e6, rdc_mohm=rdc*1e3,
                   rac_mohm=rac*1e3, ratio=rac/rdc,
                   eta_I=eta, i_max=imag.max(), i_min=imag.min(),
                   i_mean=imean),
              open("data/q2_results.json", "w"), ensure_ascii=False, indent=1)
    print("\n-> data/q2_strand_currents.csv, data/q2_results.json")
