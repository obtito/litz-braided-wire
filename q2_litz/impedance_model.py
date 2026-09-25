"""Q2/Q3 共用 L1 引擎：股线阻抗矩阵模型（解析，秒级）。

物理（Guillod COMPEL2017 / Sullivan&Zhang 2014 的思路，本实现为对数核版）：
    每丝每米压降 V_i = R_s·F(β)·I_i + jω Σ_j M_ij I_j，并联约束 V_i ≡ V，Σ I_i = I_tot
    M_ij = (μ0/2π)·ln(R_air/d_ij)（d_ij=丝心距；对角 d_ii=r_s + 内电感 μ0/8π）
    —— 外层丝链通少 → L 小 → 抢流：未绞合束的电流分配由该方程完全决定（电感性再分配）。

规则绞合（单绞向）：每丝以固定半径 ρ 匀速换方位 → 方位角平均
    ⟨ln d⟩(点ρ_i ↔ 圆ρ_j) = ln max(ρ_i, ρ_j)  ⇒  M̂ 仅依赖环半径，
    同环丝严格等价 → 问题降为「环电流方程」（7 个未知量），径向失衡保留（绞合换不了层）。
损耗 = Σ R_s F |I_i|²（欧姆，含单丝趋肤）+ Σ π d⁴ ω²μ0²|H_i|²/(64ρ)（邻近涡流，H=包络电流场）。
绞合长度因子：每股实际长 1/cosα_k，α_k=arctan(2πρ_k/P)，计入 R_s 与 Rdc 口径。
"""
import math
import sys

import numpy as np
from scipy.special import jv

sys.path.insert(0, ".")
from sim.constants import SIGMA, RHO, MU0, F0, I_RMS, skin_depth


def strand_skin_factor(r_s, f):
    """孤立圆丝趋肤因子 F(β)=Rac/Rdc（精确 Bessel）。"""
    k = (1 - 1j) / skin_depth(f)
    x = k * r_s
    z = (x / 2) * jv(0, x) / jv(1, x)
    return z.real


def build_M(centers, r_s, r_air):
    """未绞合互感矩阵（含内电感对角修正）。"""
    d = np.hypot(centers[:, None, 0] - centers[None, :, 0],
                 centers[:, None, 1] - centers[None, :, 1])
    np.fill_diagonal(d, r_s)
    M = (MU0 / (2 * math.pi)) * np.log(r_air / d)
    M[np.diag_indices_from(M)] += MU0 / (8 * math.pi)
    return M


def solve_common_V(R_diag, M, i_tot):
    """解 V·1 = R_diag·I + jω M I，ΣI = i_tot。返回 (I[n], V, Z_total=V/i_tot)。"""
    w = 2 * math.pi * F0
    Z = np.diag(R_diag) + 1j * w * M
    I1 = np.linalg.solve(Z, np.ones(len(R_diag)))   # V=1 时的电流
    V = i_tot / I1.sum()
    I = V * I1
    return I, V, V / i_tot


def untwisted(centers, r_s, r_air, f=F0, i_tot=I_RMS):
    """未绞合：逐丝 Z 矩阵直解（供 FEM 交叉验证）。"""
    F = strand_skin_factor(r_s, f)
    R_diag = np.full(len(centers), RHO / (math.pi * r_s**2) * F)
    M = build_M(centers, r_s, r_air)
    return solve_common_V(R_diag, M, i_tot)


def ring_radii(centers, spacing):
    """六角环号：axial 坐标六角距离（构造意义下的真环——按半径聚类必错，
    因为六角环本身不是等半径圆：角点与边中点半径差可达 ~14%）。"""
    # 本晶格基矢 e1=(s,0), e2=(s/2, s√3/2)：先解格坐标 (k,j)，环号=max(|k|,|j|,|k+j|)
    j_l = np.round(2 * centers[:, 1] / (math.sqrt(3) * spacing))
    k_l = np.round(centers[:, 0] / spacing - j_l / 2)
    labels = np.maximum(np.maximum(np.abs(k_l), np.abs(j_l)),
                        np.abs(k_l + j_l)).astype(int)
    rad = np.hypot(centers[:, 0], centers[:, 1])
    ring_r = np.array([rad[labels == k].mean() for k in range(labels.max() + 1)])
    return labels, ring_r


def prox_loss_per_strand(d_s, h_rms, f=F0):
    """均匀横向场中圆丝的邻近涡流损耗（每米，Sullivan P=πℓd⁴ω²μ²|H|²/(64ρ)，RMS 口径）。"""
    w = 2 * math.pi * f
    return math.pi * d_s**4 * w**2 * MU0**2 * h_rms**2 / (64 * RHO)


def azimuthal_average_model(centers, r_s, pitch, r_air, f=F0, i_tot=I_RMS, with_prox=True):
    """方位角平均环模型。pitch=绞距[m]。返回指标字典。

    ⚠️ 适用域（2026-09-25 物理勘误）：单绞向绞合是刚性旋转——所有丝共转、相对位置
    不变，内部互感网络与未绞合完全相同，本模型对方位角取平均是【错误】用法；
    该平均仅适用于 (a) 外部横向场中的束（外场在共动系中旋转）或 (b) 多级换向结构
    （子束反向旋转→丝的束轴半径真实变化）。保留本函数供 Q3 多级/编织方案的快速评估。"""
    spacing = 2 * r_s * 1.02
    labels, ring_r = ring_radii(centers, spacing)
    n_rings = len(ring_r)
    counts = np.array([(labels == k).sum() for k in range(n_rings)])
    w = 2 * math.pi * f
    F = strand_skin_factor(r_s, f)

    # 环系统（未知量=环内单丝电流 i_k，同环严格等价）：
    #   V = R_sF·i_k + jω·Σ_l A_kl·i_l
    #   A_kl = n_l·(μ0/2π)(ln R_air − ln max(ρk,ρl))   （k≠l，方位角平均 ⟨ln d⟩=ln max）
    #   A_kk = (n_k−1)(μ0/2π)(ln R_air − ln ρ_k) + (μ0/2π)ln(R_air/r_s) + μ0/8π
    c0 = MU0 / (2 * math.pi)
    lnmax = np.log(np.maximum(ring_r[:, None], ring_r[None, :]))
    A = c0 * (math.log(r_air) - lnmax) * counts[None, :]
    intra = np.array([(counts[k] - 1) * c0 * (math.log(r_air) - math.log(ring_r[k]))
                      if counts[k] > 1 else 0.0
                      for k in range(n_rings)])
    self_l = c0 * math.log(r_air / r_s) + MU0 / (8 * math.pi)
    for k in range(n_rings):
        A[k, k] = intra[k] + self_l
    R_str = RHO / (math.pi * r_s**2) * F
    Z = np.diag(np.full(n_rings, R_str)) + 1j * w * A
    i1 = np.linalg.solve(Z, np.ones(n_rings))
    V = i_tot / (i1 * counts).sum()
    i_k = V * i1          # 环内单丝电流
    Ik = i_k * counts     # 环总电流

    # 绞合长度因子（每股，逐环）
    len_fac = np.sqrt(1.0 + (2 * math.pi * ring_r / pitch) ** 2)  # 1/cos α
    rdc_str = RHO / (math.pi * r_s**2)
    ohm_loss = float(np.sum(rdc_str * F * len_fac * np.abs(i_k) ** 2 * counts))  # 同环丝 |I_i| 相等
    if with_prox:
        # 环 k 处场：包络电流近似（内环之和 + 本环一半）
        enc = np.array([Ik[:k].sum() + Ik[k] / 2 for k in range(n_rings)])
        h_ring = np.where(ring_r > 0, enc / (2 * np.pi * np.maximum(ring_r, 1e-12)), 0.0)
        prox = float(np.sum([counts[k] * prox_loss_per_strand(2 * r_s, abs(h_ring[k]), f)
                             for k in range(n_rings)]))
    else:
        prox = 0.0
    p_tot = ohm_loss + prox
    rdc_bundle = RHO / (len(centers) * math.pi * r_s**2) * np.average(len_fac, weights=counts)
    rac = p_tot / i_tot**2
    per_strand_i = i_k
    return dict(ring_labels=labels, ring_r=ring_r, counts=counts,
                ring_currents=Ik, per_strand_i=i_k, v=V,
                p_ohm=ohm_loss, p_prox=prox, p=p_tot,
                rac=rac, rdc=rdc_bundle, ratio=rac / rdc_bundle,
                eta_I=(np.abs(per_strand_i).max() - i_tot / len(centers)) / (i_tot / len(centers)),
                ring_phase_deg=np.angle(i_k, deg=True))


if __name__ == "__main__":
    import json
    from geo.strand_packing import hex_packing

    N = 127
    r_s = math.sqrt(5.0e-6 / (N * math.pi))
    centers, Rb = hex_packing(N, r_s)
    r_air = 10 * Rb  # 与 FEM 截断一致

    labels, ring_r = ring_radii(centers, 2 * r_s * 1.02)
    print("环聚类校验：", [(k, int((labels == k).sum())) for k in range(labels.max() + 1)])
    print("环半径 [mm]：", np.round(ring_r * 1e3, 3).tolist())

    # --- 未绞合 L1 vs FEM 交叉验证 ---
    I_l1, V, Z = untwisted(centers, r_s, r_air)
    rdc_b = RHO / (len(centers) * math.pi * r_s ** 2)
    print(f"未绞合 L1 校准：Re(Z)/Rdc = {Z.real / rdc_b:.4f}（FEM 4.66，L1 偏低因未含邻近反应）")
    imag = np.abs(I_l1)
    imean = I_RMS / N
    fem = np.genfromtxt("data/q2_strand_currents.csv", delimiter=",", names=True)
    print(f"\n未绞合 L1：|I| max/min = {imag.max():.4f}/{imag.min():.4f}  η_I={(imag.max()-imean)/imean*100:.1f}%")
    print("环 | L1平均|I| | FEM平均|I| | 偏差%")
    for k in range(labels.max() + 1):
        m = labels == k
        a, b = imag[m].mean(), fem["I_amplitude_A"][m].mean()
        print(f"{k:2d} | {a:.4f} | {b:.4f} | {(a/b-1)*100:+6.1f}")

    # --- 方位角平均模型（仅供 Q3 多级/编织参考，不代表单绞向绞合！） ---
    tw = azimuthal_average_model(centers, r_s, pitch=30e-3, r_air=r_air)
    print(f"\n[参考·非单绞向] 方位角平均模型（P=30mm）：R_ac/R_dc = {tw['ratio']:.4f}"
          f"（欧姆 {tw['p_ohm']:.2f} W + 邻近 {tw['p_prox']:.2f} W）"
          f"  η_I = {tw['eta_I']*100:.1f}%")
    print("环电流相对均值：", np.round((tw['per_strand_i'] / imean), 3).tolist())
    def _cvt(v):
        import numbers
        if isinstance(v, np.ndarray):
            return [_cvt(x) for x in v]
        if isinstance(v, complex):
            return {"re": v.real, "im": v.imag}
        if isinstance(v, (np.floating, np.integer)):
            return float(v)
        return v
    lenf = np.average(np.sqrt(1.0 + (2 * math.pi * tw["ring_r"] / 30e-3) ** 2),
                      weights=tw["counts"])
    print(f"\n单绞向绞合（P=30mm）解析结论：长度因子 1/cosα 平均 = {lenf:.4f}，"
          f"欧姆与 Rdc 同乘 → R_ac/R_dc 相对未绞合 4.66 变化 <±2%（刚性旋转：电流分配不变）。")
    json.dump({k: _cvt(v) for k, v in tw.items() if k != "ring_labels"},
              open("data/q2_azimuthal_avg_l1.json", "w"), ensure_ascii=False, indent=1)
    # 环电流摘要另存实值表
    np.savetxt("data/q2_ring_currents.csv",
               np.column_stack([np.arange(len(tw['ring_r'])), tw['counts'],
                                tw['ring_r'] * 1e3, np.abs(tw['per_strand_i']),
                                np.angle(tw['per_strand_i'], deg=True)]),
               delimiter=",", header="ring,count,rho_mm,I_amplitude_A,phase_deg", comments="")
    print("-> data/q2_azimuthal_avg_l1.json, data/q2_ring_currents.csv")
