"""Q4 拓扑优化：「轮转+翻面」(rotoflip) 拓扑族——RESEARCH.md §7.2 构造性方案。

设计
----
M = 2·m 个载纱器槽位排成【双层同心环】（内环 m 个 @r_in，外环 m 个 @r_out，
两环角度对齐——翻面的「方位保持」映射因此是精确等角互换）。每个载纱器携带
一束 7 丝微束（与 geo.trajectories.TwoLevelCounterTwist 的微束同构），故
载纱器级 N=M=18、丝级 N=7M=126 与 Q2/Q3 基线可比。

调度（一站 = 一次转移）：
  · 轮转：每站内、外环各自循环移位 c 个槽位（角步 Δ=2πc/m；两环等角步 →
    对整缆而言是【刚性旋转成分】——判据层是规范自由度（criterion 去旋），
    物理上单靠轮转不产生径向换位，恰是「rigid rotation」结论的调度层镜像）；
  · 翻面：每 m_flip 站将内外环【整环互换】（成对纯交换）：
      variant="preserve"：内槽 j ↔ 外槽 j（等角，方位保持）；
      variant="reverse" ：内槽 j ↔ 外槽 (−j mod m)（角度反射，方位反转）。
    翻面是唯一的径向换位机制——「轮转管方位均匀（外场通道），翻面管径向
    遍历（内部损耗通道）」两层正交（§7.2）。
  闭合：迭代至槽位占用回到恒等（闭式周期缆）；braid word 取整周期逐站转移
  的相邻对换分解（braid_check.word_from_transition）。

拓扑判定（braid_check）：单元字 u = 一个翻面周期（m_flip 站）的转移字；
  · p_min = ord(π(u))：闭包所需单元重复数（w^p∈P_N 纯辫条件）；
  · ⟨π(u)⟩ 循环传递 ⟺ π(u) 为 N-轮换（§7.1 判据）；
  · 调度层传递 transitive_group(逐站转移像)：每丝在调度中到达每槽
    （无翻面时两环永不连通 → 必不传递：轮转无法换位的群论表述）。

几何：r_in 由载纱器（微束）相切条件 2πr_in/m ≥ 2·rb_eff·gap 反解，
      r_out = r_in + 2.12·rb_eff（径向间隙留隙）；丝级坐标 = 载纱器中心 +
      随整缆共转（累计角 sΔ）的微束内 7 丝偏移。
输出：逐站坐标喂 criterion（载纱器级 + 丝级）；参数扫描 (c, m_flip,
      variant) 最小化 D；D 最优前 3 组的 braid word + 逐站坐标 +
      拓扑判定 存 data/q4_topopt.json。
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from geo.trajectories import TwoLevelCounterTwist, strand_offsets_in_bundle
from q4_perfect.braid_check import (cycle_decomposition, perm_order, slot_map,
                                    summarize_word, tietze_table,
                                    transitive_group, word_from_transition)
from q4_perfect.criterion import criterion

GAP = 1.03          # 载纱器（微束）间最小间隙系数（与 trajectories 的 bundle_gap 一致量级）
RADIAL_MARGIN = 2.12  # r_out−r_in = 2.12·rb_eff（≥2·rb_eff 相切 + 裕量）


class RotoFlip:
    """「轮转+翻面」拓扑族的一个具体设计。

    参数：m=每环载纱器数（默认 9 → M=18，微束级与基线 18 束对齐）；
    c=每站轮转步数；m_flip=每几站翻面一次（0=永不翻面，对照用）；
    variant="preserve"|"reverse"；r_in/r_out 可显式给，否则由 packing 反解。
    """

    def __init__(self, m=9, c=1, m_flip=3, variant="preserve",
                 r_in=None, r_out=None, ref=None):
        assert variant in ("preserve", "reverse")
        assert 1 <= c <= m - 1 or c == 0
        if ref is None:
            ref = TwoLevelCounterTwist()          # 只借用其微束几何（rb_eff、7 丝偏移）
        self.m, self.c, self.m_flip, self.variant = m, c % m, m_flip, variant
        self.M = 2 * m
        self.r_s = ref.r_s
        self.offsets = strand_offsets_in_bundle(ref.r_s, margin=1.02)   # (7,2)
        self.rb_eff = float(np.hypot(*self.offsets.T).max() + self.r_s)
        self.r_in = float(r_in) if r_in else self.m * self.rb_eff * GAP / math.pi * 1.04
        self.r_out = float(r_out) if r_out else self.r_in + RADIAL_MARGIN * self.rb_eff
        self._assert_layout()
        # 槽位坐标：内槽 j (0..m−1) @角度 2πj/m、半径 r_in；外槽 m+j @同角、r_out
        ang = 2 * math.pi * np.arange(m) / m
        self.slots = np.vstack([
            np.stack([self.r_in * np.cos(ang), self.r_in * np.sin(ang)], axis=1),
            np.stack([self.r_out * np.cos(ang), self.r_out * np.sin(ang)], axis=1)])

    # ---- 几何安全 ----
    def _assert_layout(self):
        m, rb = self.m, self.rb_eff
        for r, tag in ((self.r_in, "inner"), (self.r_out, "outer")):
            spacing = 2 * r * math.sin(math.pi / m)
            assert spacing >= 2 * rb * GAP, \
                f"{tag} ring spacing {spacing*1e3:.4f}mm < {2*rb*GAP*1e3:.4f}mm"
        assert self.r_out - self.r_in >= 2 * rb * GAP, "radial gap too small"

    # ---- 调度 ----
    def _rotate(self, holder):
        """holder: 位置表（槽→丝）。内/外环各循环移位 c 位（丝随槽走）。"""
        m, c = self.m, self.c
        out = list(holder)
        for j in range(m):
            out[(j + c) % m] = holder[j]                    # 内环
            out[m + (j + c) % m] = holder[m + j]            # 外环
        return out

    def _flip(self, holder):
        """内外环整环互换（成对交换丝）。preserve: 等角 j↔j；reverse: j↔−j。"""
        m = self.m
        out = list(holder)
        for j in range(m):
            jp = m + (j if self.variant == "preserve" else (m - j) % m)
            out[j], out[jp] = holder[jp], holder[j]
        return out

    def schedule(self, max_stations=4000):
        """生成闭合调度：返回 (S, holder_list, transitions)。

        holder_list[s] = 站 s 的位置表（长度 M，槽→丝）；holder_list[0]=id，
        holder_list[S]=id（闭合，S=闭合周期站数，不含末尾重复站）。
        transitions[t] = 第 t 站转移的类型标签（"rot"/"rot+flip"）。
        """
        holder = list(range(self.M))
        holders = [tuple(holder)]
        transitions = []
        for t in range(1, max_stations + 1):
            holder = self._rotate(holder)
            if self.m_flip and t % self.m_flip == 0:
                holder = self._flip(holder)
                transitions.append("rot+flip")
            else:
                transitions.append("rot")
            if holder == list(range(self.M)):
                return len(holders), holders, transitions   # 闭合：不含重复末站
            holders.append(tuple(holder))
        raise RuntimeError(f"no closure within {max_stations} stations "
                           f"(m={self.m}, c={self.c}, m_flip={self.m_flip}, "
                           f"{self.variant})")

    # ---- 坐标与轨迹 ----
    def carrier_coords(self, holders):
        """载纱器级逐站坐标 (S, M, 2)。"""
        return np.array([[self.slots[p] for p in h] for h in holders])

    def strand_coords(self, holders):
        """丝级逐站坐标 (S, 7M, 2)：载纱器中心 + 随整缆共转的 7 丝偏移。

        共转角 = s·Δ（Δ=2πc/m，翻面不改变整缆取向）；去旋判据会把它作为
        规范自由度移除，残余即微束内 7 位的相对构型。
        """
        S = len(holders)
        out = np.empty((S, 7 * self.M, 2))
        delta = 2 * math.pi * self.c / self.m
        for s, h in enumerate(holders):
            psi = s * delta
            c, sn = math.cos(psi), math.sin(psi)
            rot_off = np.stack([c * self.offsets[:, 0] - sn * self.offsets[:, 1],
                                sn * self.offsets[:, 0] + c * self.offsets[:, 1]], axis=1)
            pos = np.empty((self.M, 2))
            for p, strand in enumerate(h):
                pos[strand] = self.slots[p]
            out[s] = (pos[:, None, :] + rot_off[None, :, :]).reshape(7 * self.M, 2)
        return out

    def braid(self, holders, transitions, unit_len=None):
        """整周期 braid word（相邻对换分解）+ 单元字（一个翻面周期）。

        站序列补上闭合末站（恒等）：整周期字含 S 个转移（末转移回到恒等，
        π(word_full)=id ⇒ 整周期缆自动纯辫）。unit_len=None 时单元 = 整周期；
        否则取前 unit_len 站转移（§7.1 的「最小重复单元」——周期缆 = 单元
        重复 p 次，纯辫闭包/传递性对单元字判定）。
        返回 (word_full, word_unit, station_gen_images, transition_words)。
        """
        seq = list(holders) + [tuple(range(self.M))]     # 闭合末站
        words, gen_images = [], []
        for t in range(len(seq) - 1):
            w = word_from_transition(seq[t], seq[t + 1])
            words.append(w)
            gen_images.append(slot_map(list(seq[t + 1]), self.M))
        word_full = [g for w in words for g in w]
        k = unit_len if unit_len is not None else len(words)
        word_unit = [g for w in words[:k] for g in w]
        return word_full, word_unit, gen_images, words

    # ---- 一站式评估 ----
    def evaluate(self, n_shells=8):
        """调度 + 坐标 + 判据 + 拓扑判定，全部打包。"""
        S, holders, transitions = self.schedule()
        cc, sc = self.carrier_coords(holders), self.strand_coords(holders)
        d_car = criterion(pos=cc, n_shells=n_shells)
        d_str = criterion(pos=sc, n_shells=n_shells)
        unit_len = self.m_flip if self.m_flip else 1
        word_full, word_unit, gens, _tw = self.braid(holders, transitions, unit_len)
        topo = summarize_word(word_unit, self.M, unit_gens=gens)
        n_flips = sum(1 for tr in transitions if tr == "rot+flip")
        return {
            "params": {"m": self.m, "M": self.M, "c": self.c,
                       "m_flip": self.m_flip, "variant": self.variant,
                       "r_in_mm": self.r_in * 1e3, "r_out_mm": self.r_out * 1e3,
                       "strand_d_mm": 2 * self.r_s * 1e3},
            "S_stations": S, "n_flips": n_flips,
            "D_strand": d_str["D"], "D_dwell_strand": d_str["D_dwell"],
            "D_phi_strand": d_str["D_phi"], "phi_per_basis": d_str["phi_per_basis"],
            "D_carrier": d_car["D"], "D_dwell_carrier": d_car["D_dwell"],
            "D_phi_carrier": d_car["D_phi"],
            "area_bias_strand": d_str["area_bias"],
            "fairness_strand": d_str["fairness"],
            "topology": {"unit_stations": unit_len,
                         "p_min_closure": topo["p_min"],
                         "cyclic_transitive": topo["cyclic_transitive"],
                         "group_transitive": topo["group_transitive"],
                         "word_len_unit": topo["word_length"],
                         "word_len_full": len(word_full),
                         "unit_cycles": [list(c) for c in topo["cycles"]]},
            "_private": {"holders": holders, "transitions": transitions,
                         "carrier_coords": cc, "strand_coords": sc,
                         "word_full": word_full, "word_unit": word_unit,
                         "gen_images": gens},
        }


# ---------------------------------------------------------------- 扫描
def run_sweep(m=9, c_list=None, m_flip_list=(0, 2, 3, 4, 5, 6, 9),
              variants=("preserve", "reverse"), topk=3, json_path=None,
              verbose=True):
    """参数扫描 (c, m_flip, variant)，按丝级 D 排序取前 topk 存 JSON。"""
    c_list = list(c_list) if c_list else list(range(1, m))
    rows = []
    for variant in variants:
        for mf in m_flip_list:
            for c in c_list:
                try:
                    ev = RotoFlip(m=m, c=c, m_flip=mf, variant=variant).evaluate()
                except RuntimeError as e:            # 未闭合：记录跳过
                    if verbose:
                        print(f"  SKIP c={c} m_flip={mf} {variant}: {e}")
                    continue
                rows.append(ev)
                if verbose:
                    t = ev["topology"]
                    print(f"  c={c} m_flip={str(mf):>2} {variant:8s} "
                          f"S={ev['S_stations']:>4} flips={ev['n_flips']:>3} "
                          f"D={ev['D_strand']:.4f} "
                          f"(dwell {ev['D_dwell_strand']:.4f} phi {ev['D_phi_strand']:.4f}) "
                          f"p_min={t['p_min_closure']:>3} "
                          f"cyc_T={int(t['cyclic_transitive'])} "
                          f"grp_T={int(t['group_transitive'])}")
    rows.sort(key=lambda r: (r["D_strand"], r["topology"]["word_len_full"]))
    top = rows[:topk]
    if json_path:
        payload = {
            "family": "rotoflip: rotation (c slots/station, both rings) + "
                      "flip (inner<->outer ring swap every m_flip stations)",
            "layout": {"m_per_ring": m, "M_carriers": 2 * m,
                       "n_strands": 7 * 2 * m,
                       "note": "slot radii auto-derived from micro-bundle "
                               "packing of geo.trajectories.TwoLevelCounterTwist"},
            "criterion": "q4_perfect.criterion: D=0.5*D_dwell+0.5*D_phi, "
                         "derotated frame, 8 shells, area-weighted null",
            "n_configs": len(rows),
            "top3": [_pack(ev) for ev in top],
            "all_configs": [{k: v for k, v in r.items() if k != "_private"}
                            for r in rows],
        }
        Path(json_path).parent.mkdir(parents=True, exist_ok=True)
        with open(json_path, "w") as f:
            json.dump(payload, f, ensure_ascii=False, indent=1)
        if verbose:
            print(f"  -> wrote {json_path} ({len(rows)} configs, top{topk} packed)")
    return rows, top


def _pack(ev):
    """top-3 打包：braid word（单元 + 整期，超长仅存单元）+ 逐站坐标 + 槽位表。"""
    p = ev["_private"]
    holders = p["holders"]
    word_full_flat = [x for g in p["word_full"] for x in g]
    keep_full = len(word_full_flat) <= 60000
    return {
        "params": ev["params"], "S_stations": ev["S_stations"],
        "n_flips": ev["n_flips"],
        "D_strand": ev["D_strand"], "D_dwell_strand": ev["D_dwell_strand"],
        "D_phi_strand": ev["D_phi_strand"], "phi_per_basis": ev["phi_per_basis"],
        "D_carrier": ev["D_carrier"],
        "area_bias_strand": ev["area_bias_strand"],
        "fairness_strand": ev["fairness_strand"],
        "topology": ev["topology"],
        "braid_word_unit_flat": [x for g in p["word_unit"] for x in g],
        "unit_repeats_to_close": ev["topology"]["p_min_closure"],
        "braid_word_full_flat": word_full_flat if keep_full else None,
        "word_convention": "flat [i1,e1,i2,e2,...]; sigma_i (1-based) swaps "
                           "adjacent slots i-1,i (0-based); applied left to "
                           "right; all crossings +1 (S_N image sign-independent)",
        "slot_table": [list(h) for h in holders],      # 站 s 槽 p 持丝（位置表）
        "station_coords_carrier": np.round(p["carrier_coords"] * 1e3, 6).tolist(),
        "coords_unit": "mm",
    }


# ---------------------------------------------------------------- 自检
def run_selfcheck():
    print("== rotoflip.py self-check ==")
    rf = RotoFlip(m=9, c=1, m_flip=3, variant="preserve")
    S, holders, transitions = rf.schedule()

    # 1) 每站占用是置换（无重槽/空槽）
    for h in holders:
        assert sorted(h) == list(range(rf.M)), "occupancy must be a permutation"
    # 2) 闭合
    assert holders[-1] != tuple(range(rf.M)) or S == 1
    nxt = rf._rotate(holders[-1])
    if rf.m_flip and S % rf.m_flip == 0:
        nxt = rf._flip(nxt)
    assert tuple(nxt) == tuple(range(rf.M)), "schedule must close after S stations"
    # 3) braid word 回放：整周期 Tietze 末像 = 恒等（纯辫闭合）；逐站复现槽位表
    word_full, word_unit, gens, twords = rf.braid(holders, transitions,
                                                  unit_len=rf.m_flip)
    tab = tietze_table(word_full, rf.M)
    assert tab[-1] == tuple(range(rf.M)), "full-cycle word must close to identity"
    seq = list(holders) + [tuple(range(rf.M))]
    img = list(range(rf.M))
    for t, w in enumerate(twords):
        for i, _ in w:
            img[i - 1], img[i] = img[i], img[i - 1]
        assert img == list(seq[t + 1]), f"station {t+1} replay mismatch"
    # 4) 翻面确实交换内外环（径向历史交替）＋轮转不换环
    radii = np.array([[rf.slots[p][0] ** 2 + rf.slots[p][1] ** 2 for p in h]
                      for h in holders])
    inner_r, outer_r = rf.r_in, rf.r_out
    for strand in (0, rf.M - 1):
        rr = np.sqrt(radii[:, strand])
        assert rr.min() < inner_r + 1e-12 and rr.max() > outer_r - 1e-12, \
            "every strand must visit both rings when flips are on"
    # 5) 无翻面对照：两环不连通（群论层）→ 不传递；有翻面：调度层传递
    rf0 = RotoFlip(m=9, c=1, m_flip=0, variant="preserve")
    S0, h0, tr0 = rf0.schedule()
    _, _, gens0, _ = rf0.braid(h0, tr0, unit_len=1)
    assert not transitive_group(gens0, rf0.M), \
        "rotation alone must never connect the two rings (no radial transposition)"
    assert transitive_group(gens, rf.M), \
        "flip-enabled schedule must be transitive at group level"
    # 6) 布局安全已在构造断言；打印包络
    print(f"  layout: m=9/ring, r_in={rf.r_in*1e3:.3f}mm r_out={rf.r_out*1e3:.3f}mm "
          f"rb_eff={rf.rb_eff*1e3:.4f}mm bundle_OD={2*(rf.r_out+rf.rb_eff)*1e3:.3f}mm")
    print(f"  preserve c=1 m_flip=3: S={S} stations, {sum(1 for t in transitions if 'flip' in t)} flips, "
          f"word_unit {len(word_unit)} gens, "
          f"p_min={perm_order(slot_map(list(holders[rf.m_flip]), rf.M))}")
    print(f"  D(carrier)={criterion(pos=rf.carrier_coords(holders))['D']:.4f}  "
          f"D(strand)={criterion(pos=rf.strand_coords(holders))['D']:.4f}")
    print("  PASS: permutation occupancy / closure / word replay / ring exchange / "
          "rotation-alone not transitive / flip transitive")
    return rf


if __name__ == "__main__":
    run_selfcheck()
    print()
    print("== rotoflip sweep (c, m_flip, variant) ==")
    run_sweep(json_path=Path(__file__).resolve().parents[1] / "data" / "q4_topopt.json")
