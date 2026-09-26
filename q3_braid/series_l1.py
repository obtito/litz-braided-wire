"""Q3 串联合成求解器（series-composed L1）——换位收益的正确计算。

定理背景（findings H+17）：纯置换调度不改变每站位置多重集 → 逐站独立并联求解的
均值对换位不敏感（multiset invariance，已数值验证至 0.03%）。真实换位收益来自
串联约束：各丝电流复振幅沿整根线恒定、终端均压（每丝总压降相等）——
迫使周期 RMS 电流趋于均分，消除失衡超额损耗（Cauchy–Schwarz 严格不等）。

公式（Guillod COMPEL2017 思路，对数核实现）：
    Z_series[i,j] = Σ_s [ R_str·F(γ_s)·δ_ij + jω·M(p_i(z_s), p_j(z_s)) ]
    M(a,b) = (μ₀/2π)·ln(R_air/d)，d=|a−b|（对角 d=r_s + 内电感 μ₀/8π）
    解 V·1 = Z_series·I，ΣI=20A → I_i（沿线恒定复振幅）
    损耗 = Σ_i R_str·F·|I_i|²·L_i ＋ 逐站横向邻近（Bessel K，由 ΣI 包络场）
           ＋ h∥（轴向场修正，同 slice_engine）
已知口径：对数核 L1 对未绞合偏低 ~10%（无邻近反应）；本引擎用于【方案间相对收益】，
未绞合绝对值以 FEM 4.653 为准做锚定偏置校正。
"""
import json
import math
import sys

import numpy as np

sys.path.insert(0, ".")
from sim.constants import SIGMA, RHO, MU0, F0, I_RMS
from q2_litz.impedance_model import strand_skin_factor
from q3_braid.slice_engine import kelvin_K


def build_M_positions(pos, r_s, r_air):
    d = np.hypot(pos[:, None, 0] - pos[None, :, 0], pos[:, None, 1] - pos[None, :, 1])
    np.fill_diagonal(d, r_s)
    M = (MU0 / (2 * math.pi)) * np.log(r_air / d)
    M[np.diag_indices_from(M)] += MU0 / (8 * math.pi)
    return M


def series_solve(traj, strand_r, S=16, f=F0, i_tot=I_RMS, r_air_factor=12):
    """轨迹 → 串联电流、分站损耗、指标。返回 dict。"""
    w = 2 * math.pi * f
    N = traj.N
    period = traj.period_loss() or 1.0
    _st = getattr(traj, "stations", getattr(traj, "xy", None))
    if _st is not None:
        # 离散换位：只在站位求和（插值会让换环丝在段中相互穿越 → 对数核奇异污染 Z）
        zs = np.arange(len(_st)) * (period / len(_st))
    else:
        zs = np.linspace(0, period, S, endpoint=False) if period > 0 else np.array([0.0])
    alpha_L = traj.bundle_radius()
    r_air = r_air_factor * alpha_L
    F = strand_skin_factor(strand_r, f)
    R_str = RHO / (math.pi * strand_r ** 2) * F

    Z = np.zeros((N, N), dtype=complex)
    stations_pos = [traj.positions(z) for z in zs]
    for pos in stations_pos:
        Z += np.diag(np.full(N, R_str)) + 1j * w * build_M_positions(pos, strand_r, r_air)
    I1 = np.linalg.solve(Z, np.ones(N))
    V = i_tot / I1.sum()
    I = V * I1
    imean = i_tot / N
    # 逐站损耗：欧姆（沿线恒定电流）+ 横向邻近（包络场近似，Bessel K）+ h∥
    gamma = strand_r * math.sqrt(w * MU0 * SIGMA)
    Gs_perp = 4 * math.pi * RHO * kelvin_K(gamma)
    Gs_par = 2 * math.pi * RHO * kelvin_K(gamma)
    mbar = float(np.mean(traj.strand_length_factor()))
    p_ohm = float(np.sum(R_str * np.abs(I) ** 2) * mbar)
    p_perp = p_hpar = 0.0
    eps = (period if period > 0 else 1.0) * 1e-6
    _st2 = getattr(traj, "stations", getattr(traj, "xy", None))
    dz_step = (period / max(len(_st2), 1)) if _st2 is not None else eps
    for z, pos in zip(zs, stations_pos):
        r = np.hypot(pos[:, 0], pos[:, 1])
        # 包络场（同壳半权稳定口径）：H(r_i) = [Σ_{r_j<r_i}|I_j| + 0.5·Σ_{r_j≈r_i,j≠i}|I_j|] / (2π r_i)
        # ——修复记录：纯排序前缀和会把同半径组的包围流按排序劈开（刚性外环被低估 ~2×，
        #    换位对比失真）；半壳=壳层中点的方位平均场，物理一致。
        order = np.argsort(r)
        rs, Is = r[order], np.abs(I)[order]
        cum = np.concatenate([[0.0], np.cumsum(Is)[:-1]])          # 严格内层之和
        half = np.empty(N)
        lo = 0
        for k in range(N):                                          # 同半径组半权
            if k == 0 or rs[k] - rs[k - 1] > 1e-12:
                lo = k
            hi = k
            while hi + 1 < N and rs[hi + 1] - rs[k] <= 1e-12:
                hi += 1
            grp = Is[lo:hi + 1].sum()
            half[order[k]] = 0.5 * (grp - Is[k])
        h_perp = (cum + half) / (2 * math.pi * np.maximum(r, 1e-12))
        p_perp += float(np.sum(Gs_perp * h_perp ** 2))
        # h∥
        _, th = traj.polar(z)
        _, th2 = traj.polar((z + dz_step) % period if period > 0 else z + eps)
        dth = np.angle(np.exp(1j * (th2 - th))) / dz_step
        tan = r * np.abs(dth)
        sin = tan / np.sqrt(1 + tan ** 2)
        hs_par = i_tot * r * sin / (math.pi * alpha_L ** 2)
        p_hpar += float(np.sum(Gs_par * hs_par ** 2))
    p_perp /= len(zs); p_hpar /= len(zs)
    a_cu = N * math.pi * strand_r ** 2
    rdc = RHO * mbar / a_cu
    iabs = np.abs(I)
    eta = float(np.sqrt(((iabs - imean) ** 2).mean()) / imean)
    return dict(N=N, S=len(zs), mbar=mbar, V=V, currents=I,
                eta_I=eta, p_ohm=p_ohm, p_perp=p_perp, p_hpar=p_hpar,
                p_total=p_ohm + p_perp + p_hpar,
                rac=(p_ohm + p_perp + p_hpar) / i_tot ** 2, rdc=rdc,
                ratio=(p_ohm + p_perp + p_hpar) / i_tot ** 2 / rdc)


def calibrated_eval(traj, strand_r, f=F0, i_tot=I_RMS, S=16, refine=1, tag="cal"):
    """FEM 锚定的串联合成：ohm 精确 + prox 由单站 FEM 标定 + 解析 ΣH² 只做比值。
    绝对值锚 FEM（我们的黄金基准），解析仅承担电流再分配的相对形状——
    规避解析 prox 系数/离散环 ΣH² 的绝对误差（findings H+17 教训）。"""
    from q3_braid.slice_engine import station_solve
    # 1) 参考站 FEM：总损 + 逐丝电流（该布局的瞬时并联解）
    p_fem, I_fem, dof, a_mesh = station_solve(traj.positions(0.0), strand_r,
                                              refine=refine, path=f"data/q3_cal_{tag}.msh")
    F = strand_skin_factor(strand_r, f)
    R_str = RHO / (math.pi * strand_r ** 2) * F
    ohm_fem = float(np.sum(R_str * np.abs(I_fem) ** 2))
    prox_ref = p_fem - ohm_fem                       # FEM 精确 prox（该布局/该分配）
    # 2) 串联电流（周期轨迹的 Z 合成）
    ser = series_solve(traj, strand_r, S=S, f=f, i_tot=i_tot)
    I_ser = ser["currents"]
    # 3) ΣH² 比值（解析形状）：I_fem（站 0 分配）vs I_ser（串联分配）
    def sum_h2(pos, currents):
        r = np.hypot(pos[:, 0], pos[:, 1])
        order = np.argsort(r)
        rs, Is = r[order], np.abs(currents)[order]
        cum = np.concatenate([[0.0], np.cumsum(Is)[:-1]])
        half = np.empty(len(r))
        lo = 0
        for k in range(len(r)):
            if k == 0 or rs[k] - rs[k - 1] > 1e-12:
                lo = k
            hi = k
            while hi + 1 < len(r) and rs[hi + 1] - rs[k] <= 1e-12:
                hi += 1
            half[order[k]] = 0.5 * (Is[lo:hi + 1].sum() - Is[k])
        h = (cum + half) / (2 * math.pi * np.maximum(r, 1e-12))
        return float(np.sum(h ** 2))
    pos0 = traj.positions(0.0)
    ratio_h2 = sum_h2(pos0, I_ser) / max(sum_h2(pos0, I_fem), 1e-12)
    prox_ser = prox_ref * ratio_h2
    # 4) 总成（ohm 用串联分配；h∥ 用 series_solve 的解析值——其系数小、影响弱）
    ohm_ser = float(np.sum(R_str * np.abs(I_ser) ** 2)) * ser["mbar"]
    p_total = ohm_ser + prox_ser + ser["p_hpar"]
    rdc = ser["rdc"]
    return dict(tag=tag, N=traj.N, dof=dof, mbar=ser["mbar"], eta_I=ser["eta_I"],
                ohm=ohm_ser, prox=prox_ser, prox_ref=prox_ref, h2_ratio=ratio_h2,
                p_hpar=ser["p_hpar"], p_total=p_total, rac=p_total / i_tot ** 2,
                rdc=rdc, ratio=p_total / i_tot ** 2 / rdc)


if __name__ == "__main__":
    import sys as _sys
    if "--calibrated" in _sys.argv:
        from geo.trajectories import TwoLevelCounterTwist, RigidBundle
        from q3_braid.rotoflip_traj import RotoflipTrajectory
        from q3_braid.multiring import ThreeRingRoto, RigidFrom
        t2 = TwoLevelCounterTwist(a_cu=5.0e-6)
        r126 = RotoflipTrajectory(which=0)
        t504 = ThreeRingRoto()
        out = {}
        for tag, traj, rs in [
                ("A_untw126", RigidBundle(t2.positions(0.0), t2.r_s), t2.r_s),
                ("R_layout_rigid126", RigidBundle(r126.positions(0.0), r126.r_s), r126.r_s),
                ("C_twolevel", t2, t2.r_s),
                ("R_rotoflip126", r126, r126.r_s),
                ("A504_rigid", RigidFrom(t504), t504.r_s),
                ("R504_threering", t504, t504.r_s)]:
            res = calibrated_eval(traj, rs, S=12, tag=tag[:10])
            out[tag] = {k: v for k, v in res.items()}
            print(f"[{tag:20s}] ratio={res['ratio']:.4f} η_I={res['eta_I']*100:5.1f}% "
                  f"(ohm {res['ohm']/400*1e3:5.2f} + prox {res['prox']/400*1e3:5.2f} "
                  f"[ref {res['prox_ref']/400*1e3:.2f} × H2比 {res['h2_ratio']:.2f}] "
                  f"+ h∥ {res['p_hpar']/400*1e3:5.2f} mΩ/m, m̄={res['mbar']:.4f})")
        json.dump(out, open("data/q3_calibrated.json", "w"), ensure_ascii=False, indent=1, default=float)
        print("-> data/q3_calibrated.json")
        raise SystemExit(0)
    from geo.trajectories import TwoLevelCounterTwist, RigidBundle
    from q3_braid.rotoflip_traj import RotoflipTrajectory
    from q3_braid.multiring import ThreeRingRoto, RigidFrom

    out = {}
    schemes = []
    t2 = TwoLevelCounterTwist(a_cu=5.0e-6)
    r126 = RotoflipTrajectory(which=0)
    schemes += [("A_untw126", RigidBundle(t2.positions(0.0), t2.r_s), t2.r_s),
                ("C_twolevel", t2, t2.r_s),
                ("R_rotoflip126", r126, r126.r_s),
                ("R_layout_rigid126", RigidBundle(r126.positions(0.0), r126.r_s), r126.r_s)]
    t504 = ThreeRingRoto()
    schemes += [("A504_rigid", RigidFrom(t504), t504.r_s),
                ("R504_threering", t504, t504.r_s)]
    for tag, traj, rs in schemes:
        S_traj = min(16, getattr(traj, "S", 16) or 16)
        res = series_solve(traj, rs, S=max(S_traj, 4))
        out[tag] = {k: v for k, v in res.items() if k != "currents"}
        out[tag]["i_max_ratio"] = float(np.abs(res["currents"]).max() / (I_RMS / traj.N))
        print(f"[{tag:18s}] ratio={res['ratio']:.4f} η_I={res['eta_I']*100:5.1f}% "
              f"(ohm {res['p_ohm']/400*1e3:5.2f} + prox {res['p_perp']/400*1e3:5.2f} "
              f"+ h∥ {res['p_hpar']/400*1e3:5.2f} mΩ/m, m̄={res['mbar']:.4f})")
    json.dump(out, open("data/q3_series_l1.json", "w"), ensure_ascii=False, indent=1, default=float)
    print("-> data/q3_series_l1.json")
