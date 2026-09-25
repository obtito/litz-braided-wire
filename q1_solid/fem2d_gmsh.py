"""Q1：实心圆铜导线 2D 涡流场 FEM —— gmsh 网格 + NGSolve 求解（A-v 公式）。

（netgen CSG2D 对同心圆相减会死循环，故网格生成全部交给 gmsh 的距离场加密，
NGSolve 只负责装配与求解——这也正是 Q2/Q3 丝束横截面要走的生产路线。）

模型：MQS 相量（RMS，e^{jωt}）
    未知量：A_z（H¹，外边界 Dirichlet=0）+ v（全局标量=每米电压降）
    E_z = -(jωA_z + v)；J_z = σE_z
    方程 A：(1/μ)∇A·∇w = ∫_cu σ(jωA+v) w
    方程 v：∫_cu σ(jωA+v) = -I（总流约束，已按铜面积归一匹配 NumberFES 尺度）
    输出：R_ac = P/I²，P = ∫|J|²/σ；交叉校验 Re[v]/I；J(r) 曲线。
"""
import math
import sys

import numpy as np
import gmsh
import ngsolve as ng

sys.path.insert(0, ".")
from sim.constants import SIGMA, RHO, MU0, MU_R, F0, I_RMS, skin_depth


def make_mesh(a=1e-3, f=F0, refine=1, order=2, r_air_factor=30, path="data/q1_mesh.msh"):
    """gmsh 生成同心圆网格：铜盘 + 空气域，距离场加密（表面 0.45δ）。"""
    delta = skin_depth(f)
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("q1")
    r_air = r_air_factor * a
    occ = gmsh.model.occ
    disk = occ.addDisk(0, 0, 0, a, a)
    air = occ.addDisk(0, 0, 0, r_air, r_air)
    occ.cut([(2, air)], [(2, disk)], removeObject=True, removeTool=False)  # 保留铜盘！
    occ.synchronize()
    # 材料分组
    surfaces = gmsh.model.getEntities(2)
    cu_tag = [s[1] for s in surfaces
              if abs(gmsh.model.occ.getCenterOfMass(2, s[1])[0]) < 1e-12
              and _area(s[1]) < math.pi * a * a * 1.1]
    gmsh.model.addPhysicalGroup(2, cu_tag, name="copper")
    air_tags = [s[1] for s in surfaces if s[1] not in cu_tag]
    gmsh.model.addPhysicalGroup(2, air_tags, name="air")
    gmsh.model.addPhysicalGroup(1, [e[1] for e in gmsh.model.getEntities(1) if _on_outer(e[1], r_air)],
                                name="outer")
    # 距离场加密
    f_dist = gmsh.model.mesh.field.add("Distance")
    gmsh.model.mesh.field.setNumbers(f_dist, "EdgesList", [e[1] for e in gmsh.model.getEntities(1)
                                                           if not _on_outer(e[1], r_air)])
    f_th = gmsh.model.mesh.field.add("Threshold")
    gmsh.model.mesh.field.setNumber(f_th, "InField", f_dist)
    gmsh.model.mesh.field.setNumber(f_th, "SizeMin", 0.45 * delta / refine)
    gmsh.model.mesh.field.setNumber(f_th, "SizeMax", 3e-3)
    gmsh.model.mesh.field.setNumber(f_th, "DistMin", 0.0)
    gmsh.model.mesh.field.setNumber(f_th, "DistMax", 15 * delta)
    gmsh.model.mesh.field.setAsBackgroundMesh(f_th)
    gmsh.option.setNumber("Mesh.Algorithm", 6)  # Frontal-Delaunay
    gmsh.model.mesh.generate(2)
    # netgen 的 gmsh 导入器只认旧版 MSH 2.2 ASCII（4.x 会读出空网格）
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.option.setNumber("Mesh.Binary", 0)
    gmsh.write(path)
    gmsh.finalize()
    return path


def _area(tag):
    return gmsh.model.occ.getMass(2, tag)


def _on_outer(tag, r_air):
    bb = gmsh.model.getBoundingBox(1, tag)
    return abs(max(abs(bb[0]), abs(bb[3])) - r_air) < 1e-6 * r_air


def solve(mesh_path="data/q1_mesh.msh", a=1e-3, f=F0, i_rms=I_RMS, order=2):
    w = 2 * math.pi * f
    mesh = ng.Mesh(mesh_path)
    mesh.Curve(order)

    fesA = ng.H1(mesh, order=order, complex=True, dirichlet="outer")
    fesV = ng.FESpace("number", mesh, complex=True)
    fes = fesA * fesV
    (A, v), (Aw, vw) = fes.TnT()
    area_cu = math.pi * a * a

    bf = ng.BilinearForm(fes)
    bf += (1 / (MU0 * MU_R)) * ng.grad(A) * ng.grad(Aw) * ng.dx
    bf += -SIGMA * (1j * w * A + v) * Aw * ng.dx("copper")
    bf += -SIGMA * (1j * w * A + v) * vw * ng.dx("copper")
    lf = ng.LinearForm(fes)
    lf += (-i_rms / area_cu) * vw * ng.dx("copper")
    bf.Assemble()
    lf.Assemble()

    gfu = ng.GridFunction(fes)
    gfu.vec.data = bf.mat.Inverse(freedofs=fes.FreeDofs(), inverse="umfpack") * lf.vec
    return mesh, gfu


def evaluate(mesh, gfu, a=1e-3, f=F0, i_rms=I_RMS):
    w = 2 * math.pi * f
    A, vcf = gfu.components
    v = complex(vcf.vec[0])
    E = -(1j * w * A + v)
    Jz = SIGMA * E
    p = ng.Integrate(ng.Conj(Jz) * Jz / SIGMA, mesh, definedon=mesh.Materials("copper"))
    p = p.real if isinstance(p, complex) else p
    rdc = RHO / (math.pi * a * a)
    rs = np.linspace(0.0, a, 400)
    jr = np.empty_like(rs)
    for i, x in enumerate(rs):
        xq = min(max(x, 0.0), a * (1 - 1e-9))  # 界面/边界点夹回内部
        jr[i] = abs(SIGMA * (-(1j * w * A(mesh(float(xq), 0.0)) + v)))
    jr[0] = jr[1]
    return dict(p=p, rac=p / i_rms**2, rdc=rdc, ratio=p / i_rms**2 / rdc, r=rs, jr=jr, v=v)


def exact_ratio(a=1e-3, f=F0):
    from scipy.special import jv
    k = (1 - 1j) / skin_depth(f)
    x = k * a
    z = k * RHO * jv(0, x) / (2 * math.pi * a * jv(1, x))
    return z.real / (RHO / (math.pi * a * a))


if __name__ == "__main__":
    import time

    def fes_ndof(gfu):
        return len(gfu.vec)

    a, f = 1e-3, F0
    tgt = exact_ratio(a, f)
    print(f"δ={skin_depth(f)*1e6:.2f}µm  精确解目标 Rac/Rdc = {tgt:.5f}\n")
    print(f"{'refine':>7} {'ndof':>8} {'R_ac/R_dc':>11} {'误差%':>8} {'P[W/m]':>8} {'耗时s':>7}")
    for refine in (1, 2, 4):
        t0 = time.time()
        path = make_mesh(a=a, f=f, refine=refine)
        mesh, gfu = solve(path, a=a, f=f)
        res = evaluate(mesh, gfu, a=a, f=f)
        dt = time.time() - t0
        err = (res["ratio"] / tgt - 1) * 100
        print(f"{refine:7d} {fes_ndof(gfu):8d} {res['ratio']:11.4f} {err:8.2f} {res['p']:8.3f} {dt:7.1f}")
        if refine == 2:
            np.savetxt("data/q1_fem_J_r.csv", np.column_stack([res["r"], res["jr"]]),
                       delimiter=",", header="r_m,J_amplitude", comments="")
            print("  -> data/q1_fem_J_r.csv")
    print("\n校验：Re[v]/I 与 P/I² 应一致（同一电阻的两种提取）")
