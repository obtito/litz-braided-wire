"""Q4 完美换位判据计算器（任务 A1）。

物理含义
---------
litz/换位导线抑制趋肤+邻近损耗的机理是「均等化」：若每根丝在一个结构周期内
对泄漏场的暴露历史完全相同，则各支路电动势相等、环流为零（Pyrhönen 教科书
"perfect transposition"；所有经典 litz 损耗模型都建立在此假设上）。
本模块把「暴露历史相同」量化为两个离散指标（RESEARCH.md §7.1）：

D_dwell  径向停留分布偏差：把束截面切成 L 个等宽壳层（默认 8 层），统计每丝
         在各壳层的停留比例 p_i(ℓ)。完美 ⇔ 所有丝的 p_i 相同（跨丝公平）。
         零假设另报一个「面积均匀偏差」：完美设计的平均停留 p̄(ℓ) 还应接近
         壳层面积权重 q(ℓ)（= 该壳层面积占比）——衡量槽位集本身对截面的
         覆盖均匀性（设计属性，不进合成 D，单独报告）。

D_phi    φ-矩离散度：取与圆柱多极展开匹配的基函数
         Φ = {r, r·cosθ, r·sinθ, r², r²cos2θ, r²sin2θ}，
         逐丝计算周期平均 F_b(k)=⟨φ_b(pos_k)⟩_s，基内离散度
         (max_k−min_k)F_b / (2·max|φ_b|)，再对基取加权平均。
         物理对应：一阶矩 ↔ 外场/邻近耦合的偶极分量，二阶矩 ↔ 趋肤场的
         径向权重（r² 主导）。归一化因子 2max|φ_b| ≈ 该基的「对径静态」
         上界（半数丝冻结在最正、半数冻结在最负 → D_b=1）。

合成     D = w_dwell·D_dwell + w_phi·D_phi（默认各 0.5，可配）。

自检（__main__）
  1) 所有丝同轨迹 → D=0（离散度定义保证）；
  2) 静止未绞合（rigid）→ D 大（接近对径上界）；
  3) 完全随机轨迹 → D 很小（大数定律：随机混合同样均等化！）——
     这不是 bug，而是判据的正确行为：连续层指标只度量「均等化」，
     随机辫的失败在拓扑层（w^p 非纯辫、无周期闭合、不可制造），
     由 braid_check.py 负责。两层互补，见报告。
  4) 循环赛轮转 → D≈0。

接口说明（trajectories.py 缺口）
  geo.trajectories.radial_history 只返回 r(z)，丢弃了 φ-矩所需的 θ。
  本模块自带 sample_positions(traj, ...) 采全 (S,N,2)，不改 trajectories.py。
"""
from __future__ import annotations

import math

import numpy as np

BASIS_NAMES = ("r", "r cosθ", "r sinθ", "r²",
               "r²cos2θ", "r²sin2θ")
DEFAULT_PHI_WEIGHTS = {name: 1.0 / 6 for name in BASIS_NAMES}


def basis_values(pos):
    """pos: (S,N,2) 截面坐标序列 → dict 基名 -> (S,N) 基函数值。"""
    pos = np.asarray(pos, dtype=float)
    x, y = pos[..., 0], pos[..., 1]
    r = np.hypot(x, y)
    th = np.arctan2(y, x)
    return {
        "r": r,
        "r cosθ": x,
        "r sinθ": y,
        "r²": r * r,
        "r²cos2θ": r * r * np.cos(2 * th),
        "r²sin2θ": r * r * np.sin(2 * th),
    }


def phi_dispersion(pos, weights=None):
    """D_phi：各基的跨丝周期平均离散度 dict + 加权平均。

    每基归一化：range_k F_b(k) / (2·max_{s,k}|φ_b|)，对径静态=1，全同=0。
    """
    B = basis_values(pos)
    if weights is None:
        weights = DEFAULT_PHI_WEIGHTS
    per_basis = {}
    for name, v in B.items():
        F = v.mean(axis=0)                       # (N,) 每丝周期平均
        scale = 2.0 * float(np.abs(v).max())
        per_basis[name] = 0.0 if scale <= 0 else float((F.max() - F.min()) / scale)
    total_w = sum(weights.get(k, 0.0) for k in per_basis)
    d_phi = (sum(weights.get(k, 0.0) * v for k, v in per_basis.items())
             / total_w) if total_w > 0 else 0.0
    return {"per_basis": per_basis, "D_phi": float(d_phi)}


def dwell_deviation(pos=None, radii=None, n_shells=8, bundle_r=None):
    """D_dwell：壳层停留分布的跨丝最大偏差 + 平均剖面的面积均匀偏差。

    输入：pos (S,N,2) 或 radii (S,N)（radial_history 的输出，仅够本指标）。
    返回 dict：D_dwell（headline，完美=0）、area_bias（p̄ vs 面积权 q 的
    总变差，衡量槽位集覆盖）、shell_edges、p_mean。
    """
    if radii is None:
        if pos is None:
            raise ValueError("need pos or radii")
        radii = np.hypot(np.asarray(pos)[..., 0], np.asarray(pos)[..., 1])
    radii = np.asarray(radii, dtype=float)
    S, N = radii.shape
    if bundle_r is None:
        bundle_r = float(radii.max())
    edges = np.linspace(0.0, bundle_r, n_shells + 1)
    idx = np.clip(np.digitize(radii, edges[1:-1], right=True), 0, n_shells - 1)
    p = np.zeros((N, n_shells))
    for l in range(n_shells):
        p[:, l] = (idx == l).mean(axis=0)        # 每丝壳层停留比例
    p_mean = p.mean(axis=0)
    d_dwell = float(np.abs(p - p_mean[None, :]).max())
    q = (edges[1:] ** 2 - edges[:-1] ** 2) / bundle_r ** 2   # 面积权重
    area_bias = float(0.5 * np.abs(p_mean - q).sum())
    return {"D_dwell": d_dwell, "area_bias": area_bias,
            "shell_edges": edges.tolist(), "p_mean": p_mean.tolist()}


def criterion(pos=None, radii=None, n_shells=8,
              w_dwell=0.5, w_phi=0.5, phi_weights=None):
    """合成判据 D 及全部分解量。pos 与 radii 至少给一个（phi 需要 pos）。"""
    out = {"n_shells": n_shells, "w_dwell": w_dwell, "w_phi": w_phi}
    dw = dwell_deviation(pos=pos, radii=radii, n_shells=n_shells)
    out.update({k: dw[k] for k in ("D_dwell", "area_bias", "p_mean")})
    if pos is not None:
        ph = phi_dispersion(pos, weights=phi_weights)
        out["phi_per_basis"] = ph["per_basis"]
        out["D_phi"] = ph["D_phi"]
    else:
        out["D_phi"] = None
        w_phi = 0.0
    out["D"] = w_dwell * out["D_dwell"] + w_phi * (out["D_phi"] or 0.0)
    return out


def sample_positions(traj, period_len, n_samples=240):
    """替代 radial_history 的全信息采样：返回 (S,N,2)。"""
    zs = np.linspace(0.0, period_len, n_samples)
    return np.array([np.asarray(traj(z), dtype=float) for z in zs])


# ---------------------------------------------------------------- self-check
def _static_two_ring(m_in=3, m_out=6, r_in=1.0, r_out=2.0):
    slots = [(r_in * math.cos(2 * math.pi * j / m_in),
              r_in * math.sin(2 * math.pi * j / m_in)) for j in range(m_in)]
    slots += [(r_out * math.cos(2 * math.pi * j / m_out),
               r_out * math.sin(2 * math.pi * j / m_out)) for j in range(m_out)]
    return np.array(slots)[None, :, :]           # (1,N,2) 静止


def _random_visits(slots, S=240, seed=7):
    rng = np.random.default_rng(seed)
    N = len(slots)
    return slots[rng.integers(0, N, size=(S, N))]


def _round_robin(slots, S=240):
    N = len(slots)
    return np.array([slots[(np.arange(N) + s) % N] for s in range(S)])


def _fmt(res):
    phi = res.get("D_phi")
    phi_s = f"{phi:.4f}" if phi is not None else "  n/a "
    return (f"D={res['D']:.4f}  D_dwell={res['D_dwell']:.4f}  "
            f"D_phi={phi_s}  area_bias={res['area_bias']:.4f}")


def run_selfcheck():
    print("== criterion.py self-check ==")
    slots = _static_two_ring()[0]
    # identical：4 根丝走同一条轨迹（重叠采样）⇒ 离散度定义为 0
    same = np.repeat(slots[0][None][None], 240, axis=0).repeat(4, axis=1)
    cases = {
        "identical-traj (4 strands)": same,
        "rigid static two-ring": _static_two_ring(),
        "random iid slot visits": _random_visits(slots),
        "round-robin (perfect)": _round_robin(slots),
    }

    res = {k: criterion(pos=v, n_shells=8) for k, v in cases.items()}
    for k, v in res.items():
        print(f"  {k:28s} {_fmt(v)}")

    d_id = res["identical-traj (4 strands)"]["D"]
    d_rg = res["rigid static two-ring"]["D"]
    d_rd = res["random iid slot visits"]["D"]
    d_rr = res["round-robin (perfect)"]["D"]
    assert d_id < 1e-12, "identical trajectories must give D=0"
    assert d_rg > 0.3, "rigid bundle should be far from perfect"
    assert d_rd < d_rg, "random mixing equalizes (LLN): D_random < D_rigid"
    assert d_rr < d_rd, "round-robin should beat random"
    print(f"  PASS: identical=0; rigid={d_rg:.3f} (near upper bound); "
          f"random={d_rd:.3f} (LLN floor, not upper bound -- topological "
          f"layer needed); round-robin={d_rr:.4f}")
    print("  NOTE: random scores LOW because D measures equalization; random "
          "braid word fails purity/periodicity -> see braid_check.py.")
    return res


if __name__ == "__main__":
    run_selfcheck()
