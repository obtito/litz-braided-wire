"""TASK A —— 赛题 2 任务 1：丝径 d–δ 权衡扫描（"为什么选这个丝径"）。

固定 A_cu = 5.0 mm²（J=4 A/mm² @ 20 A），扫丝径 d ∈ [0.06, 0.40] mm：
  N_exact = A_cu/(π(d/2)²) → 取最近完整六角数 N_hex = 3k(k+1)+1 → 反解 d_actual。
每档计算：
  β = d_actual/δ, F_skin(精确 Bessel), 未绞合 L1 逐丝矩阵解 ratio_L1 = Re(Z)/Rdc,
  铜填充率, 束外径, 以及 L1_avg_ref（方位角平均环模型参考列，见下）。

⚠️ 口径声明（写入 JSON notes，引用时必须带上）：
  1) L1 未绞合模型是丝状（filament）近似：设计点 N=127 处 Re(Z)/Rdc = 4.171，
     而 FEM = 4.66/4.68 —— 系统性偏低 ~10%（无邻近反应项）。本扫描只用于
     【形状/趋势与设计空间论证】，不用于报绝对值。
  2) L1_avg_ref 列来自 impedance_model.azimuthal_average_model（pitch→∞ 即未绞合、
     长度因子=1）。该模型对互感做方位角平均——按 2026-09-25 物理勘误，这对单绞向
     绞合是【错误】用法（刚性旋转不改内部分配）；此处仅作环平均 H 邻近损耗的
     下界参考（averaging model, NOT valid for single-lay）。其邻近项用 Sullivan
     P∝d⁴|H|² 小丝径公式，d≳δ 时偏高，大丝径端只看量级。
  3) 完整六角数按定义 3k(k+1)+1 生成。任务书列表给到 k=20 (N=1261)，但
     d=0.06 mm 档 N_exact≈1768 > 1261，故按公式延拓到 k=26（列表本身不足以覆盖
     扫描下端，延拓不改变定义）。
  4) 物理规则：单绞向绞合不改变内部电流分配（刚性旋转）——本扫描的未绞合结论
     对单绞向直接适用（仅差长度因子 1/cosα ≈ 1.00x）。

输出：data/q2_tradeoff.csv（全表）、data/q2_tradeoff.json（关键数）。
不生成图（主会话负责）。
"""
import json
import math
import sys

import numpy as np

sys.path.insert(0, ".")

from sim.constants import RHO, F0, I_RMS, skin_depth
from geo.strand_packing import hex_packing, copper_fill
from q2_litz.impedance_model import (
    strand_skin_factor,
    untwisted,
    azimuthal_average_model,
)

A_CU = 5.0e-6          # m²，固定铜截面
DELTA = skin_depth(F0)  # 147.77 um @200 kHz
R_DC_BUNDLE = RHO / A_CU  # 3.4483 mΩ/m —— 所有档位相同（N·πr²≡A_cu）
PITCH_INF = 1e9        # m，pitch→∞ ⇒ 长度因子 1/cosα = 1（未绞合口径）

# 完整六角数（中心+6 环）：k=0..26。任务书列出 k≤20，此处按同一公式延拓覆盖
# d=0.06 mm 档（N_exact≈1768 → N_hex=1801=k24）。
HEX_NUMBERS = [3 * k * (k + 1) + 1 for k in range(27)]


def nearest_hex(n):
    return min(HEX_NUMBERS, key=lambda h: (abs(h - n), h))


def scan_point(d_req, is_design=False):
    """单个丝径档的全部指标。d_req: 名义丝径 [m]。"""
    n_exact = A_CU / (math.pi * (d_req / 2) ** 2)
    n_hex = nearest_hex(n_exact)
    d_act = math.sqrt(4 * A_CU / (math.pi * n_hex))
    r_s = d_act / 2

    centers, bundle_r = hex_packing(n_hex, r_s)
    r_air = 10 * bundle_r  # 与 FEM 截断口径一致

    F = strand_skin_factor(r_s, F0)
    I, V, Z = untwisted(centers, r_s, r_air)
    # 任务书口径：ratio_L1 = Re(Z)/RHO · (N·π·r_s²) ≡ Re(Z)/R_dc_bundle
    ratio_l1 = Z.real / RHO * (n_hex * math.pi * r_s ** 2)
    imb = float(np.sum(np.abs(I) ** 2) / (n_hex * (I_RMS / n_hex) ** 2))  # Σ|I|²/(N·Ī²)

    # 参考列：方位角平均环模型（未绞合口径 pitch→∞），⚠️ 仅下界参考，见模块 docstring
    # （np.errstate：环 0 半径为 0 → ln0=-inf 警告；该对角元随后被覆盖，无害）
    with np.errstate(divide="ignore"):
        avg = azimuthal_average_model(centers, r_s, pitch=PITCH_INF, r_air=r_air)

    return {
        "d_req_mm": d_req * 1e3,
        "design_point": is_design,
        "N_exact": n_exact,
        "N_hex": n_hex,
        "d_actual_mm": d_act * 1e3,
        "beta": d_act / DELTA,
        "F_skin": float(F),
        "ratio_L1": float(ratio_l1),
        "imbalance_L1": float(imb),              # = ratio_L1/F_skin（欧姆再分配因子）
        "L1_avg_ref": float(avg["ratio"]),       # averaging model, NOT valid for single-lay
        "L1_avg_p_ohm_W_per_m": float(avg["p_ohm"]),
        "L1_avg_p_prox_W_per_m": float(avg["p_prox"]),  # d⁴ 公式，d>δ 端偏高，仅量级
        "fill": float(copper_fill(n_hex, r_s, bundle_r)),
        "OD_mm": 2 * bundle_r * 1e3,
        "R_ac_L1_ohm_per_m": float(Z.real),
        "P_L1_W_per_m": float(Z.real * I_RMS ** 2),
    }


def main():
    print(f"δ(200 kHz) = {DELTA*1e6:.2f} um   R_dc(bundle) = {R_DC_BUNDLE*1e3:.4f} mΩ/m "
          f"(恒定, N·πr²≡A_cu)\n六角数表 k=0..26: {HEX_NUMBERS[:8]}...{HEX_NUMBERS[-3:]}")

    d_grid = list(np.linspace(0.06e-3, 0.40e-3, 35)) + [0.2239e-3]
    rows = []
    for i, d in enumerate(d_grid):
        r = scan_point(d, is_design=(i == 35))
        rows.append(r)
        print(f"  d_req={r['d_req_mm']:.4f} mm  N={r['N_hex']:5d}  d_act={r['d_actual_mm']:.4f} mm  "
              f"β={r['beta']:.3f}  F={r['F_skin']:.4f}  ratio_L1={r['ratio_L1']:.4f}")

    # ---- 锚点校验（不过关即失败）----
    des = rows[35]
    assert abs(des["ratio_L1"] - 4.171) < 0.01, f"设计点 ratio_L1={des['ratio_L1']} ≠ 4.171"
    assert abs(des["F_skin"] - 1.00682) < 2e-4, f"设计点 F_skin={des['F_skin']} ≠ 1.00682"
    # 边界锚：N=1（等面积实心导线）—— L1 退化为单丝精确 Bessel
    solid = strand_skin_factor(math.sqrt(A_CU / math.pi), F0)
    assert abs(solid - 4.5295) < 2e-3, f"N=1 实心锚 {solid} ≠ 4.5295"
    print(f"\n锚点校验通过：N=127 设计点 ratio_L1={des['ratio_L1']:.4f} (锚 4.171), "
          f"F_skin={des['F_skin']:.5f} (锚 1.00682), β={des['beta']:.4f}; "
          f"N=1 实心边界 F_skin={solid:.4f} (锚 4.5295)")

    # ---- CSV ----
    cols = ["d_req_mm", "design_point", "N_exact", "N_hex", "d_actual_mm", "beta",
            "F_skin", "ratio_L1", "imbalance_L1", "L1_avg_ref",
            "L1_avg_p_ohm_W_per_m", "L1_avg_p_prox_W_per_m", "fill", "OD_mm",
            "R_ac_L1_ohm_per_m", "P_L1_W_per_m"]
    with open("data/q2_tradeoff.csv", "w") as f:
        f.write(",".join(cols) + "\n")
        for r in rows:
            f.write(",".join(
                "1" if r[c] is True else "0" if r[c] is False else f"{r[c]:.6g}"
                for c in cols) + "\n")
    print("-> data/q2_tradeoff.csv")

    # ---- 代表性行（按唯一 N 去重后取代表档）----
    uniq = {}
    for r in rows:
        uniq.setdefault(r["N_hex"], r)
    show = sorted(uniq.values(), key=lambda r: -r["d_actual_mm"])
    print("\n  d_act[mm] |   N  |  β=d/δ | F_skin | ratio_L1 | imb_L1 | L1_avg_ref | fill | OD[mm]")
    for r in show:
        tag = "  <- 设计点" if r["design_point"] else ""
        print(f"  {r['d_actual_mm']:.4f}   | {r['N_hex']:4d} | {r['beta']:6.3f} | "
              f"{r['F_skin']:6.4f} | {r['ratio_L1']:8.4f} | {r['imbalance_L1']:6.3f} | "
              f"{r['L1_avg_ref']:10.4f} | {r['fill']:.3f} | {r['OD_mm']:.3f}{tag}")

    # ---- 趋势分析（按唯一 N 去重，d 升序）----
    curve = sorted(uniq.values(), key=lambda r: r["d_actual_mm"])  # d 小 → 大
    rmin = min(curve, key=lambda r: r["ratio_L1"])
    rmax = max(curve, key=lambda r: r["ratio_L1"])

    n_coarse = max(uniq, key=lambda n: uniq[n]["d_actual_mm"])       # 最粗档 (N=37)
    n_fine_01 = nearest_hex(A_CU / (math.pi * (0.10e-3 / 2) ** 2))   # d≈0.10 mm 档
    r_coarse, r_fine = uniq[n_coarse]["ratio_L1"], uniq[n_fine_01]["ratio_L1"]

    # 平坦区（细径端）：取【最大】的 d，使 d'≤d 全部档位 ratio_L1 极差 < 3%×最细档值
    r_at_min_d = curve[0]["ratio_L1"]
    d_flat = curve[0]["d_actual_mm"]
    for r in curve:  # d 升序，一旦极差超阈值即停
        below = [x["ratio_L1"] for x in curve if x["d_actual_mm"] <= r["d_actual_mm"]]
        if max(below) - min(below) >= 0.03 * r_at_min_d:
            break
        d_flat = r["d_actual_mm"]

    out = {
        "meta": {
            "f_Hz": F0, "delta_um": DELTA * 1e6, "A_cu_mm2": A_CU * 1e6,
            "I_rms_A": I_RMS, "R_dc_bundle_mohm_per_m": R_DC_BUNDLE * 1e3,
            "grid": "linspace(0.06,0.40,35) mm + d=0.2239 mm (design)",
            "hex_numbers": "3k(k+1)+1, k=0..26 (task list k<=20 extended per definition "
                           "to cover N_exact(0.06mm)=1768 -> N=1801)",
        },
        "design_point": {
            "d_mm": des["d_actual_mm"], "N_hex": des["N_hex"], "beta": des["beta"],
            "F_skin": des["F_skin"], "ratio_L1": des["ratio_L1"],
            "imbalance_L1": des["imbalance_L1"], "L1_avg_ref": des["L1_avg_ref"],
            "fill": des["fill"], "OD_mm": des["OD_mm"],
            "R_ac_L1_ohm_per_m": des["R_ac_L1_ohm_per_m"],
            "P_L1_W_per_m": des["P_L1_W_per_m"],
        },
        "argmin_ratio_L1": {
            "d_actual_mm": rmin["d_actual_mm"], "d_req_mm": rmin["d_req_mm"],
            "N_hex": rmin["N_hex"], "ratio_L1": rmin["ratio_L1"],
            "caveat": "argmin sits at the COARSE end where L1 filament bias is largest "
                      "(F_skin=1.076, beta=2.8) — do NOT read as 'thicker is better'; "
                      "see headline_finding",
        },
        "argmax_ratio_L1": {
            "d_actual_mm": rmax["d_actual_mm"], "N_hex": rmax["N_hex"],
            "ratio_L1": rmax["ratio_L1"],
        },
        "trend": {
            "coarse_end": {"N_hex": n_coarse, "d_mm": uniq[n_coarse]["d_actual_mm"],
                           "ratio_L1": r_coarse},
            "d_0p10mm_region": {"N_hex": n_fine_01,
                                "d_mm": uniq[n_fine_01]["d_actual_mm"],
                                "ratio_L1": r_fine},
            "improvement_0p40_to_0p10": {
                "abs": r_coarse - r_fine,
                "rel_pct": (r_coarse - r_fine) / r_coarse * 100,
                "verdict": "NEGATIVE — untwisted ratio gets WORSE toward fine strands; "
                           "'improvement' only exists for F_skin, not for the untwisted bundle"},
            "finest_end": {"N_hex": curve[0]["N_hex"], "d_mm": curve[0]["d_actual_mm"],
                           "ratio_L1": r_at_min_d,
                           "note": "N=1801 filament limit 4.5326 ~= N=1 solid equal-area "
                                   "exact 4.5295 — untwisted fine bundle converges to "
                                   "solid-conductor behavior"},
            "flat_below_mm": d_flat,
            "flat_criterion": "largest d such that max-min of ratio_L1 over all d'<=d "
                              "< 3% of the finest-d value (fine end flattens ~4.53-4.56)",
            "curve_unique_N": [[r["N_hex"], r["d_actual_mm"], r["beta"],
                                r["F_skin"], r["ratio_L1"], r["imbalance_L1"],
                                r["L1_avg_ref"]] for r in curve],
        },
        "headline_finding": (
            "UNTWISTED scan: finer strands do NOT reduce Rac/Rdc. ratio_L1 rises "
            "3.43 (N=37, d=0.415mm) -> 4.17 (design N=127) -> ~4.53-4.56 flat for "
            "d<=~0.15mm, converging onto the solid-conductor value of the same copper "
            "area (exact 4.5295). Cause: inductive current crowding to the bundle "
            "periphery is set by bundle OD (~2.8-3.0mm, D/delta ~ 19 >> 1), not by "
            "strand size — imbalance_L1 stays 3.2-4.6 (current effectively on 1/3 of "
            "strands). Strand-diameter selection therefore CANNOT be justified by the "
            "untwisted curve: d is chosen so (a) per-strand skin F_skin(d) ~ 1 "
            "(0.68% penalty at the design point beta=1.515), and (b) once "
            "transposition/braiding equalizes positions (kills the 3-4.6x imbalance), "
            "the residual AC loss is F_skin + proximity budget, both monotone in d. "
            "Untwisted-Litz ~= solid is the classic result this scan quantifies; it is "
            "the reason Problem-2 needs the braided/transposed structure of Q3."
        ),
        "notes": [
            "L1 untwisted is a filament model with documented ~10% UNDERestimate vs FEM "
            "(design point N=127: L1 4.171 vs FEM 4.66-4.68, no proximity-reaction term); "
            "the bias GROWS with beta, so coarse-end rows (N=37/61, beta>2) are even more "
            "underestimated. Scan is for SHAPE/trend and design-space reasoning only.",
            "L1_avg_ref column: azimuthal_average_model with pitch->inf (len factor=1). "
            "Averaging model, NOT valid for single-lay twist (rigid rotation keeps internal "
            "current sharing unchanged); used ONLY as a ring-averaged-H lower-bound reference. "
            "Its d^4 prox term is a small-beta formula and OVERestimates for d>delta: at "
            "N=37 the prox part exceeds ohmic (see L1_avg_p_prox column) — coarse-end "
            "L1_avg_ref is magnitude-only, do not quote.",
            "Single-lay twist does not change internal current sharing (rigid rotation); "
            "untwisted conclusions carry over up to length factor 1/cos(alpha) ~ 1.00x.",
            "R_dc bundle = rho/A_cu = 3.4483 mohm/m identical for all rows by construction "
            "(N*pi*r_s^2 == A_cu exactly via d_actual).",
            "Boundary anchors passed: design point L1 4.1706 vs 4.171, F_skin 1.00682; "
            "N=1 solid equal-area strand_skin_factor 4.5295 (exact Bessel).",
        ],
    }
    with open("data/q2_tradeoff.json", "w") as f:
        json.dump(out, f, ensure_ascii=False, indent=1)
    print("-> data/q2_tradeoff.json")

    print(f"\n趋势摘要：argmin ratio_L1 = {rmin['ratio_L1']:.4f} @ d={rmin['d_actual_mm']:.4f} mm "
          f"(N={rmin['N_hex']}, 粗端·L1 偏差最大处，勿读作'越粗越好')")
    print(f"  粗端 N={n_coarse} (d={uniq[n_coarse]['d_actual_mm']:.3f} mm): {r_coarse:.4f}  →  "
          f"d≈0.10 mm (N={n_fine_01}): {r_fine:.4f}   "
          f"Δ = {(r_coarse-r_fine)/r_coarse*100:+.1f}%（未绞合：细径不降反升）")
    print(f"  最细端 N={curve[0]['N_hex']} (d={curve[0]['d_actual_mm']:.4f} mm): {r_at_min_d:.4f} "
          f"≈ N=1 实心等铜面积锚 4.5295 —— 未绞合细丝束 → 实心行为")
    print(f"  细径端平坦区：d ≤ {d_flat:.3f} mm（极差 <3% 判据）")


if __name__ == "__main__":
    main()
