"""Q4 拓扑层判据：辫群 word → S_N 像（Tietze 跟踪表）、纯辫闭包、传递性、逐站槽位表。

理论（RESEARCH.md §7.1 拓扑层）
--------------------------------
交叉序列 = braid word w = σ_{i1}^{±1}σ_{i2}^{±1}… ∈ B_N（关系 σᵢσᵢ₊₁σᵢ=σᵢ₊₁σᵢσᵢ₊₁，
|i−j|≥2 可交换）。满同态 π: B_N → S_N（σᵢ ↦ (i, i+1)）。周期缆「完美」的拓扑必要条件：

  ① 纯辫闭包：w^p ∈ P_N = ker π（各丝回原位）⟺ π(w)^p = id；
     最小 p = ord(π(w)) = π(w) 轮换长度的 lcm（恒存在，但可能很大——
     「闭合所需周期数」本身就是设计指标）。
  ② 传递性：⟨π(w)⟩ ≤ S_N 传递作用（每丝到达每槽）。⟨π(w)⟩ 由 w 的像生成；
     单个 word 的像生成循环群 → 传递 ⟺ π(w) 是 N-轮换（单轨道）。
     多生成元版本 transitive_group 接受逐站转移像列表（调度层传递性）。

Tietze 跟踪表（本模块的「Tietze 表」）：逐前缀计算 S_N 像——位置表 holder[p]=
当前位于位置 p 的丝编号，每读一个 σᵢ 交换 holder[i−1]↔holder[i]（置换像与
交叉符号无关，故 ±1 像相同；符号只影响辫的同痕类，不影响 ①② 判定）。

word 约定：生成元列表 [(i, ±1), …]，1 ≤ i ≤ N−1；σ_i 交换相邻位置 (i−1, i)
（0 基位置），word 从左到右依次作用。word_from_transition 把「相邻两站的
槽位占用表」分解为相邻对换（选择排序式冒泡，确定性、每步 ≤N−1 个生成元，
全取 +1 代表——S_N 像与分解无关，判定不变）。

自检（任务书：轮转字 (σ1σ2…σ_{N−1})^N 的性质，N=5）
  ρ = σ1σ2…σ4：π(ρ) = 5-轮换（⟨π(ρ)⟩ 传递）、ord = 5；
  ρ^5：π = id ⇒ ρ^5 ∈ P_5（纯辫，闭包幂 p=1；即轮转字需 5 个周期闭合）；
  负例：σ1（ord=2 纯但 ⟨π⟩ 不传递，N≥3）、σ1σ3（两个对换，不传递）；
  辫关系 σ1σ2σ1 = σ2σ1σ2 的 S_N 像相等（3-轮换）；
  ρ 的整周期幂 ρ^k：k<5 时 π≠id（未闭合）。
"""
from __future__ import annotations

import math

import numpy as np


# ---------------------------------------------------------------- Tietze 表
def tietze_table(word, n):
    """word: [(i, ±1), …] → 逐前缀 S_N 像表 [π_0, π_1, …, π_L]。

    每个 π 是长度 n 的元组：π[p] = 站 p（该前缀后）位于位置 p 的丝编号。
    π_0 = 恒等。置换像与生成元符号无关（σ 与 σ^{-1} 同像）。
    """
    if not (1 <= n):
        raise ValueError("n >= 1 required")
    holder = list(range(n))
    table = [tuple(holder)]
    for g in word:
        i, sgn = (g if isinstance(g, (tuple, list)) else (g, +1))
        if not (1 <= i <= n - 1):
            raise ValueError(f"generator index {i} out of range 1..{n-1}")
        holder[i - 1], holder[i] = holder[i], holder[i - 1]
        table.append(tuple(holder))
    return table


def perm_image(word, n):
    """整条 word 的 S_N 像（位置表，π[p]=位于 p 的丝）。"""
    return np.array(tietze_table(word, n)[-1])


def slot_history(word, n):
    """逐站槽位表：strand→slot 视角（Tietze 表的转置解读）。

    返回 (L+1, n) 数组，行 k = 前 k 个生成元后每根丝所在的位置（槽位）。
    """
    table = np.array(tietze_table(word, n))          # (L+1, n) 位置→丝
    out = np.empty_like(table)
    for k, row in enumerate(table):
        out[k, row] = np.arange(n)                   # 丝 s 在位置 row[s]
    return out


# ---------------------------------------------------------------- 群论判定
def cycle_decomposition(perm):
    """perm: 位置表（π[p]=位于 p 的丝）或丝→槽位表（二者轮换结构相同）。

    返回轮换列表（每个是位置/丝的元组，长度≥2，不动点略去）。
    """
    perm = list(perm)
    n = len(perm)
    seen = [False] * n
    cycles = []
    for start in range(n):
        if seen[start]:
            continue
        cyc, p = [], start
        while not seen[p]:
            seen[p] = True
            cyc.append(p)
            p = perm[p]
        if len(cyc) > 1:
            cycles.append(tuple(cyc))
    return cycles


def perm_order(perm):
    """ord(π) = 轮换长度 lcm = 纯辫闭包最小幂 p（π(w)^p = id）。"""
    cycles = cycle_decomposition(perm)
    if not cycles:
        return 1
    return math.lcm(*(len(c) for c in cycles))


def orbits(group_perms, n):
    """⟨生成元集⟩ 作用下的轨道（BFS 闭包）。group_perms: 丝→槽位映射列表。"""
    gens = [np.asarray(p) for p in group_perms]
    seen = {False: set(), True: set()}
    orb = set()
    frontier = [0]
    while frontier:
        nxt = []
        for x in frontier:
            if x in orb:
                continue
            orb.add(x)
            for g in gens:
                nxt.append(int(g[x]))
        frontier = [y for y in nxt if y not in orb]
    return orb


def transitive_group(group_perms, n):
    """⟨π(w₁),π(w₂),…⟩ 是否在 {0..n−1} 上传递（单轨道 = 每丝可达每槽）。"""
    return len(orbits(group_perms, n)) == n


def summarize_word(word, n, unit_gens=None):
    """一条 word（或逐站转移像列表 unit_gens）的拓扑层体检表。

    word 模式：p_min = ord(π(w))（w^p∈P_N 的最小 p）；cyclic_transitive =
    ⟨π(w)⟩ 传递 ⟺ π(w) 为 N-轮换。
    unit_gens 模式（调度层）：把每个站转移的 S_N 像作为生成元，报
    group_transitive = ⟨全部转移⟩ 传递性（「每丝在调度中到达每槽」的
    操作性判据——闭式周期缆 w_full 的 π=id，必须用逐站生成元才见传递性）。
    """
    out = {"n": n, "word_length": len(word)}
    table = tietze_table(word, n)
    out["perm_image"] = table[-1]
    out["cycles"] = cycle_decomposition(table[-1])
    out["p_min"] = perm_order(table[-1])
    out["pure_at_p"] = True            # π(w)^{p_min}=id 恒成立；p_min 即闭合周期数
    out["cyclic_transitive"] = transitive_group([slot_map(table[-1], n)], n)
    if unit_gens is not None:
        out["group_transitive"] = transitive_group(unit_gens, n)
        out["n_unit_gens"] = len(unit_gens)
    return out


def slot_map(position_holder, n):
    """位置表 π[p]=丝 → 丝→槽位映射（orbits 用的生成元方向）。"""
    position_holder = np.asarray(position_holder, dtype=int)
    m = np.empty(n, dtype=int)
    m[position_holder] = np.arange(n)
    return m


# ---------------------------------------------------------------- word 构造
def word_from_transition(occ_prev, occ_next):
    """相邻两站槽位占用 → 相邻对换分解（确定性，选择排序式冒泡）。

    occ[p] = 站上位于位置 p 的丝编号（位置表）。产出 [(i,+1), …] 使
    依次作用后 occ_prev → occ_next；长度 ≤ N(N−1)/2。
    """
    cur = list(occ_prev)
    tgt = list(occ_next)
    n = len(cur)
    word = []
    for p in range(n):
        q = cur.index(tgt[p], p)          # 目标丝必在 ≥p 处（p 前已就位）
        while q > p:
            cur[q - 1], cur[q] = cur[q], cur[q - 1]
            word.append((q, +1))          # σ_q 交换 0 基位置 q−1, q
            q -= 1
    assert cur == tgt
    return word


def word_pretty(word):
    """[(i,±1),…] → 'σ2 σ3^{-1} …' 字符串。"""
    return " ".join(f"σ{i}" if s > 0 else f"σ{i}⁻¹" for i, s in word)


def rotation_word(n):
    """轮转字 ρ = σ1σ2…σ_{n−1}（自检与文档锚点）。"""
    return [(i, +1) for i in range(1, n)]


# ---------------------------------------------------------------- 自检
def run_selfcheck(n=5):
    print("== braid_check.py self-check ==")
    rho = rotation_word(n)
    w = rho * n                                        # (σ1σ2…σ_{n-1})^n

    pi_rho = perm_image(rho, n)
    cyc = cycle_decomposition(pi_rho)
    assert len(cyc) == 1 and len(cyc[0]) == n, "rotation word must be an n-cycle"
    assert perm_order(pi_rho) == n, "n-cycle has order n"

    tab = tietze_table(w, n)
    assert tab[-1] == tuple(range(n)), "(σ1…σ_{n-1})^n maps to identity in S_N"
    assert perm_image(w, n).tolist() == list(range(n))
    assert transitive_group([slot_map(pi_rho, n)], n), \
        "⟨π(ρ)⟩ must be transitive (n-cycle)"
    # k < n 的幂未闭合
    for k in range(1, n):
        assert perm_image(rho * k, n).tolist() != list(range(n)), \
            f"ρ^{k} must not close for k < n"
    # 纯辫闭包语义：w=ρ^n 本身 π=id ⇒ p_min=1（w∈P_N）；ρ 需 p=n 个周期
    assert perm_order(perm_image(w, n)) == 1
    # 负例
    assert not transitive_group([slot_map(perm_image([(1, +1)], n), n)], n), \
        "single σ1 not transitive for n>=3"
    assert not transitive_group([slot_map(perm_image([(1, +1), (3, +1)], n), n)], n), \
        "σ1σ3 not transitive"
    # 辫关系 σ1σ2σ1 = σ2σ1σ2（S_N 像相等）
    a = perm_image([(1, 1), (2, 1), (1, 1)], 3)
    b = perm_image([(2, 1), (1, 1), (2, 1)], 3)
    assert a.tolist() == b.tolist(), "braid relation must hold in S_N image"
    # word_from_transition 回放一致性
    occ0 = tuple(range(n))
    occ1 = tuple(tab[7])                               # 某中间前缀
    word_rt = word_from_transition(occ0, occ1)
    assert perm_image(word_rt, n).tolist() == list(occ1), "roundtrip failed"
    # 逐站槽位表
    hist = slot_history(w, n)
    assert hist.shape == (len(w) + 1, n)
    assert hist[0].tolist() == list(range(n)) and hist[-1].tolist() == list(range(n))

    print(f"  ρ=σ1…σ{n-1}: π(ρ) is {n}-cycle (transitive), ord={n} "
          f"⇒ ρ needs {n} periods to close; ρ^{n}: π=id ⇒ pure braid (p_min=1), "
          f"word length {len(w)}")
    print(f"  negatives: σ1 / σ1σ3 not transitive ✓; braid relation image ✓; "
          f"transition roundtrip ✓")
    print("  PASS")
    return {"rho_image": pi_rho.tolist(), "word_image": tab[-1]}


if __name__ == "__main__":
    run_selfcheck()
