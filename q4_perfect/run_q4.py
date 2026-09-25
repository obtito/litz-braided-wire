"""Q4 总装：判据自检 → 三基线对照（判据-物理一致性验证）→ rotoflip 扫描 → 汇总。

基线对照（任务 4，关键）
------------------------
① RigidBundle(TwoLevelCounterTwist.positions(0), twist_pitch=None)：未绞合直束；
②a RigidBundle(同 centers, twist_pitch=25mm)：单绞向的「分配不变」建模
   （twist 只进长度因子，positions 静止）→ D 应与 ① 完全相同；
②b RigidRotator(同 centers, pitch=25mm)：单绞向的【字面】轨迹——整束刚体
   旋转 pos(z)=R(2πz/P)·centers（实验室系轨迹完全不同！）。判据去旋后
   D(②b) 必须逐位回到 D(①)——「单绞向绞合=刚性旋转、不产生任何换位收益」
   在判据层的定理式验证；同时报告 lab 框架的 D_phi 会被虚假改善的数值，
   说明为什么 lab 框架不能作 headline。
③ TwoLevelCounterTwist 完整几何周期（25mm = lcm(P₁,P₂)）：两级反向绞的
   部分换位（微束内自转 ρ≤0.229mm 的径向迁移；微束圆心 |C| 不交换、
   18 根中心丝 ρ=0 完全不动）→ D 应低于 ①②，且降幅被上述两个结构性
   地板限制——与 Q3 计划的「未绞 5.16 → P1 等化 4.77 → 地板 4.43」
   （PAPERS_NOTES）小增益带一致。

输出：data/q4_baselines.json；与 rotoflip 扫描（data/q4_topopt.json，
rotoflip.run_sweep 生成）合并成汇总表 + Umetani 理论地板衔接说明。

地板正交性：R_ideal/Rdc(d=0.2248)=4.43、d=0.1124=1.92（Umetani 解析地板，
PAPERS_NOTES）——拓扑优化只压 D（换位充分性），丝径 d 才动地板；两轴正交，
最优设计 = 低 D × 细 d 的乘积叙事（commit 241c068）。
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from geo.trajectories import RigidBundle, TwoLevelCounterTwist
from q4_perfect import braid_check, criterion, rotoflip
from q4_perfect.criterion import criterion as crit, sample_positions

ROOT = Path(__file__).resolve().parents[1]


class RigidRotator:
    """字面单绞向：整束刚体旋转 pos(z)=R(2πz/P)·centers（相对构型恒定）。"""

    def __init__(self, centers, pitch):
        self.centers = np.asarray(centers, dtype=float)
        self.pitch = pitch
        self.N = len(self.centers)

    def positions(self, z):
        a = 2 * math.pi * z / self.pitch
        c, s = math.cos(a), math.sin(a)
        return self.centers @ np.array([[c, -s], [s, c]]).T

    __call__ = positions

    def polar(self, z):
        p = self.positions(z)
        return np.hypot(*p.T), np.arctan2(p[:, 1], p[:, 0])


def _run_case(name, traj, period, n_samples=256, n_shells=8):
    pos = sample_positions(traj, period, n_samples)
    res = crit(pos=pos, n_shells=n_shells)                      # derotated headline
    lab = crit(pos=pos, n_shells=n_shells, frame="lab")
    return {"name": name, "period_mm": period * 1e3, "n_samples": n_samples,
            "D": res["D"], "D_dwell": res["D_dwell"], "D_dwell_max": res["max_tv"],
            "fairness": res["fairness"], "area_bias": res["area_bias"],
            "D_phi": res["D_phi"], "phi_per_basis": res["phi_per_basis"],
            "D_phi_lab": lab["D_phi"], "D_lab": lab["D"],
            "phi_per_basis_lab": lab["phi_per_basis"]}, pos


def run_baselines(json_path=None, verbose=True):
    t = TwoLevelCounterTwist()                    # N=126, P1=12.5, P2=25mm
    centers = t.positions(0.0)
    P = t.period_geometric()                      # 25mm

    b1 = RigidBundle(centers, t.r_s, twist_pitch=None)           # ① 未绞合
    b2a = RigidBundle(centers, t.r_s, twist_pitch=25e-3)         # ②a 单绞向（分配不变）
    b2b = RigidRotator(centers, 25e-3)                           # ②b 字面刚体旋转
    cases = [
        _run_case("① untwisted (RigidBundle, pitch=None)", b1, P),
        _run_case("②a single-lay (RigidBundle, pitch=25mm)", b2a, P),
        _run_case("②b single-lay literal (rigid rotator 25mm)", b2b, 25e-3),
        _run_case("③ two-level counter-twist (full period)", t, P),
    ]
    results = [c[0] for c in cases]

    # ---- 判据-物理一致性断言 ----
    d1, d2a, d2b, d3 = (r["D"] for r in results)
    assert d2a == d1, "RigidBundle with/without pitch share positions: D identical"
    assert abs(d2b - d1) < 1e-9, \
        "literal rigid rotation must derotate back to D(untwisted) exactly"
    assert results[2]["D_phi_lab"] < results[0]["D_phi_lab"] - 0.1, \
        "lab frame would credit single-lay twist (spurious) -- keep derotated headline"
    assert d3 < d1, "two-level counter-twist must beat untwisted on D"

    if verbose:
        print("== Q4 baselines (criterion-physics consistency) ==")
        for r in results:
            print(f"  {r['name']:48s} D={r['D']:.4f} "
                  f"(dwell {r['D_dwell']:.4f}, phi {r['D_phi']:.4f} | "
                  f"lab phi {r['D_phi_lab']:.4f})")
        print(f"  ✓ ②a==① exactly (allocation unchanged); "
              f"②b derotates back to ① (|ΔD|={abs(d2b-d1):.1e}) while lab frame "
              f"would show D_phi {results[0]['D_phi_lab']:.3f}->"
              f"{results[2]['D_phi_lab']:.3f} (spurious credit)")
        print(f"  ✓ ③ D={d3:.4f} < ①② D={d1:.4f} "
              f"({(1-d3/d1)*100:.0f}% relative reduction; partial transposition "
              f"-- floored by |C| ring split + 18 dead centre strands)")

    payload = {"criterion": "D=0.5*D_dwell(mean TV to area-weighted null, 8 shells)"
                            "+0.5*D_phi(6-basis, derotated frame)",
               "setup": {"N": t.N, "P1_mm": t.p1 * 1e3, "P2_mm": t.p2 * 1e3,
                         "period_mm": P * 1e3, "samples": 256,
                         "strand_d_mm": 2 * t.r_s * 1e3,
                         "bundle_OD_mm": 2 * t.bundle_radius() * 1e3},
               "baselines": results,
               "assertions": {"d2a_equals_d1": d2a == d1,
                              "d2b_minus_d1": d2b - d1,
                              "d3_lt_d1": d3 < d1}}
    if json_path:
        with open(json_path, "w") as f:
            json.dump(payload, f, ensure_ascii=False, indent=1)
        if verbose:
            print(f"  -> wrote {json_path}")
    return results


def main():
    print("=============== Q4 perfect-transposition criterion & topology ===============")
    criterion.run_selfcheck()
    print()
    braid_check.run_selfcheck()
    print()
    rotoflip.run_selfcheck()
    print()
    baselines = run_baselines(json_path=ROOT / "data" / "q4_baselines.json")
    print()
    print("== rotoflip sweep (c, m_flip, variant) — D minimization ==")
    rows, top = rotoflip.run_sweep(
        json_path=ROOT / "data" / "q4_topopt.json", verbose=False)

    # ---- 汇总表 ----
    print("== summary: D across designs (derotated criterion, 8 shells) ==")
    for r in baselines:
        print(f"  {r['name']:48s} D={r['D']:.4f}")
    flip_rows = [r for r in rows if r["n_flips"] > 0
                 and r["topology"]["group_transitive"]]
    best = flip_rows[0]
    long_cycle = min((r for r in flip_rows
                      if r["S_stations"] >= 18),
                     key=lambda r: r["D_strand"], default=None)
    print(f"  rotoflip best (D-optimal, transitive): c={best['params']['c']}, "
          f"m_flip={best['params']['m_flip']}, {best['params']['variant']}: "
          f"D={best['D_strand']:.4f} (dwell {best['D_dwell_strand']:.4f}, "
          f"phi {best['D_phi_strand']:.4f}), S={best['S_stations']} stations, "
          f"closure p={best['topology']['p_min_closure']}")
    if long_cycle:
        p = long_cycle["params"]
        print(f"  rotoflip manufacturable pick (long cycle, few flips): c={p['c']}, "
              f"m_flip={p['m_flip']}, {p['variant']}: "
              f"D={long_cycle['D_strand']:.4f}, S={long_cycle['S_stations']}, "
              f"flips={long_cycle['n_flips']} (D within "
              f"{(long_cycle['D_strand']/best['D_strand']-1)*100:.1f}% of best)")
    noflip = [r for r in rows if r["n_flips"] == 0]
    if noflip:
        print(f"  no-flip control (rotation only): D={noflip[0]['D_strand']:.4f} "
              f"— rotation alone never transposes (topology: grp_T=0), "
              f"mirrors the rigid-rotation physics")
    print()
    print("== Umetani theory floor (orthogonal axis) ==")
    print("  R_ideal/Rdc floor vs strand diameter d (PAPERS_NOTES, "
          "Umetani 2021 analytic):")
    print("    d=0.2248mm -> 4.43 | d=0.1124mm -> 1.92 | d=0.08mm -> 1.49")
    print("  Topology optimization moves D only (equalization quality); strand")
    print("  diameter moves the floor (m̄·F(γs) + 3.369·(d/0.2248)²). The two are")
    print("  orthogonal design axes; the optimum is their product: rotoflip-class")
    print("  D at d=0.1124 would approach 1.92, which no amount of topology at")
    print("  d=0.2248 can reach (its floor is 4.43).")


if __name__ == "__main__":
    main()
