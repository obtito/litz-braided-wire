"""Q3 切片引擎：换位方案 → R_ac/R_dc 的 2D 多站切片仿真（L2 主路线）。

方法（PAPERS_NOTES/BATTLE_PLAN 定稿口径）：
  损耗相关周期 = P₁（cage 约定下 r(z) 只依赖微束自转——刚性旋转定理的工程红利）。
  沿 P₁ 取 S 站：每站 positions(z_s) → gmsh 网格 → scikit-fem A-v 多导体解
  → P_s（W/m 轴长）、逐丝电流 I_i(z_s)。
  R_ac 分解：
    ① R_2d = mean_s(P_s)·m̄ / I²   （2D 可见损耗：丝内趋肤 + 丝间横向邻近；
       m̄=全束平均长度因子——丝内与丝间损耗同按丝长标度，一阶近似）
    ② P_h∥ = Σ_i ⟨Gs∥·|hs∥,i(z)|²⟩（2D 盲区解析修正：绞合自生轴向场驱动的
       平行场丝邻近损耗，Umetani (8)(22)；hs∥=I·r·sinθ/(π·α_L²)，逐站由轨迹求）
    ③ R_ac = (①+②)/I² ；R_dc = ρ·m̄/A_cu（决策 5 口径）→ R_ac/R_dc。
  输出 η_I 的周期 RMS（电流失衡的公平度量）与逐丝 RMS 电流（喂 D↔Rac 标定）。
"""
import json
import math
import sys
import time

import numpy as np
from scipy.special import jv

sys.path.insert(0, ".")
from sim.constants import SIGMA, RHO, MU0, F0, I_RMS, skin_depth
from sim.gmsh_mesh import make_bundle_mesh
from q2_litz.solve_bundle import load_bundle, solve_bundle, strand_currents


def kelvin_K(x):
    """Umetani (11)：K(x) = −x·[ber₂·ber′+bei₂·bei′]/(ber²+bei²)，ber₂+i·bei₂=J₂(x·e^{3iπ/4})。"""
    rot = x * np.exp(3j * np.pi / 4)
    b2r, b2i = jv(2, rot).real, jv(2, rot).imag
    ber, bei = jv(0, rot).real, jv(0, rot).imag
    h = 1e-6
    berp = (jv(0, (x + h) * np.exp(3j * np.pi / 4)).real - jv(0, (x - h) * np.exp(3j * np.pi / 4)).real) / (2 * h)
    beip = (jv(0, (x + h) * np.exp(3j * np.pi / 4)).imag - jv(0, (x - h) * np.exp(3j * np.pi / 4)).imag) / (2 * h)
    return -x * (b2r * berp + b2i * beip) / (ber ** 2 + bei ** 2)


def station_solve(positions, strand_r, refine=1, r_air_factor=10, path="data/q3_station.msh"):
    """一站：网格 + 解 + 逐丝电流。"""
    from geo.strand_packing import hex_packing  # noqa: F401（保持导入路径一致）
    bundle_r = np.hypot(*positions.T).max() + strand_r
    make_bundle_mesh(positions, strand_r, refine=refine,
                     r_air=r_air_factor * bundle_r, path=path)
    mesh, ind_cu, strand_ind = load_bundle(path, positions, strand_r)
    res = solve_bundle(mesh, ind_cu)
    Ii, _ = strand_currents(res["basis"], strand_ind, res["a"], res["v"])
    return res["p"], Ii, res["dof"], res["area_cu"]


def h_parallel_power(traj, station_data, strand_r, f=F0, i_tot=I_RMS):
    """P_h∥（W/m）：逐站逐丝 hs∥ = i_tot·r·sinθ/(π·α_L²)，θ=该丝该处的绞合切角。
    Gs∥ = 2πρ·K(γ_s)（Umetani (8)；γ_s=√2·a_s/δ 半径制）。"""
    w = 2 * math.pi * f
    gamma = strand_r * math.sqrt(w * MU0 * SIGMA)
    Gs = 2 * math.pi * RHO * kelvin_K(gamma)
    alpha_L = traj.bundle_radius()
    acc = 0.0
    for z, Ii in station_data:
        r, th = traj.polar(z)
        # 切角 tanθ = r·|dθ/dz|（数值微分，步长 1e-6·P₁）
        eps = traj.period_loss() * 1e-6 if traj.period_loss() > 0 else 1e-8
        _, th2 = traj.polar(z + eps)
        dth = np.angle(np.exp(1j * (th2 - th))) / eps
        tan = r * np.abs(dth)
        sin = tan / np.sqrt(1 + tan ** 2)
        hs = i_tot * (alpha_L - r) * tan / (math.pi * alpha_L ** 2)  # Umetani(22) 修正：(α_L−r) 结构（Biot-Savart 验证）
        acc += float(np.sum(Gs * hs ** 2))
    return acc / len(station_data)


def run_scheme(traj, strand_r, S=16, refine=1, tag="", verbose=True):
    """轨迹对象（positions/polar/period_loss/strand_length_factor/bundle_radius）→ 指标。"""
    t0 = time.time()
    period = traj.period_loss()
    if period > 0:
        zs = np.linspace(0, period, S, endpoint=False)
    else:
        zs = np.array([0.0])          # 刚性/直束：一站即足
    ps, currents, dofs, areas = [], [], [], []
    for k, z in enumerate(zs):
        p, Ii, dof, a_cu = station_solve(traj.positions(z), strand_r, refine=refine,
                                         path=f"data/q3_st_{tag}_{k}.msh")
        ps.append(p); currents.append(Ii); dofs.append(dof); areas.append(a_cu)
    mbar_raw = traj.strand_length_factor()
    mbar = float(np.mean(mbar_raw)) if np.ndim(mbar_raw) > 0 else float(mbar_raw)
    P2d = float(np.mean(ps)) * mbar
    station_data = list(zip(zs, currents))
    Php = h_parallel_power(traj, station_data, strand_r)
    I_abs = np.abs(np.array(currents))                 # (S,N)
    I_rms = np.sqrt((I_abs ** 2).mean(axis=0))         # 逐丝周期 RMS
    imean = I_RMS / traj.N
    eta = float(np.sqrt(((I_abs.mean(axis=0) - imean) ** 2).mean()) / imean)  # RMS 偏差
    a_cu_exact = traj.N * math.pi * strand_r ** 2
    rdc = RHO * mbar / a_cu_exact
    rac = (P2d + Php) / I_RMS ** 2
    out = dict(tag=tag, N=traj.N, S=len(zs), dof=int(np.mean(dofs)),
               A_mesh=float(np.mean(areas)) * 1e6, A_exact=a_cu_exact * 1e6,
               mbar=float(mbar), P_2d=P2d, P_hpar=Php, rac=rac, rdc=rdc,
               ratio=rac / rdc, eta_I_rms=eta,
               i_rms_min=float(I_rms.min()), i_rms_max=float(I_rms.max()),
               stations=len(zs), runtime=time.time() - t0)
    if verbose:
        print(f"[{tag}] N={traj.N} S={len(zs)} dof={out['dof']} | "
              f"R_ac/R_dc = {out['ratio']:.4f}（2D {P2d/I_RMS**2*1e3:.2f} + h∥ {Php/I_RMS**2*1e3:.2f} mΩ/m；"
              f"m̄={mbar:.4f}）| η_I(rms)={eta*100:.1f}% | {out['runtime']:.0f}s")
    return out


if __name__ == "__main__":
    from geo.trajectories import TwoLevelCounterTwist, RigidBundle

    results = {}
    # 方案 A：未绞合对照（同 N=126/d 同几何，严格控制变量）
    t2 = TwoLevelCounterTwist(a_cu=5.0e-6)                 # P1=12.5/P2=25
    ctrl = RigidBundle(t2.positions(0.0), t2.r_s)
    results["A_untwisted"] = run_scheme(ctrl, t2.r_s, S=2, tag="A")
    # 方案 C：两级反向绞
    results["C_twolevel"] = run_scheme(t2, t2.r_s, S=16, tag="C")
    # C 的 S 收敛性（S=8 对照）
    results["C_twolevel_S8"] = run_scheme(t2, t2.r_s, S=8, tag="C8")
    json.dump(results, open("data/q3_first_results.json", "w"),
              ensure_ascii=False, indent=1, default=float)
    print("-> data/q3_first_results.json")
