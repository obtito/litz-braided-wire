"""Q4 完美换位判据计算器（定稿版，2026-09-25）。

物理含义（RESEARCH.md §7.1 三层判据的「离散+连续混合」层）
---------------------------------------------------------
litz/换位导线抑制趋肤+邻近损耗的机理是「均等化」：若每根丝在一个结构周期内对
泄漏场的暴露历史完全相同，则各支路电动势相等、环流为零（Pyrhönen 教科书
"perfect transposition"；经典 litz 损耗模型均建立在此假设上）。本模块把
「暴露历史相同」量化为两个指标 + 合成 D：

D_dwell  径向壳层停留分布对【面积加权零假设】的偏差（任务书口径）：截面切
         L 个等宽壳层（默认 8），壳层 ℓ 面积权重 q(ℓ)=(r²_{ℓ+1}−r²_ℓ)/R²。
         逐丝统计停留比例 p_i(ℓ)（各丝 Σ_ℓ p_i=1），headline 取
         mean_i TV(p_i, q)，TV=½Σ_ℓ|·| 为总变差距离。
         D_dwell=0 ⇔ 每丝停留分布=面积均匀（既跨丝公平、又真正覆盖截面）。
         子指标：max_i TV(p_i,q)（最差丝）、fairness=mean_i TV(p_i,p̄)
         （纯跨丝公平，槽位集覆盖无关）、area_bias=TV(p̄,q)（槽位集设计属性，
         衡量槽位布局本身对截面的面积均匀覆盖——双层环等离散布局的地板来源）。

D_phi    φ-矩离散度：基 Φ={r, r·cosθ, r·sinθ, r², r²cos2θ, r²sin2θ}（与圆柱
         多极展开匹配），逐丝周期平均 F_b(k)=⟨φ_b(pos_k)⟩_s，基内离散度
         (max_k−min_k)F_b / (2·max|φ_b|)（「对径静态」上界=1，全同=0），
         对基取算术平均（权重可配）。物理对应：r/r² ↔ 趋肤场的径向权重
         （内部损耗主通道），一次/二次角矩 ↔ 外场横置与多极邻近分量。

  frame 框架参数（关键）：
    "derotated"（默认）—— 计算矩之前逐站去除【整体刚性旋转】（对站 0 的
        Kabsch 单自由度最优旋转角）。理由：孤立束内部损耗只依赖相对构型
        （自场轴对称；单绞向绞合=整体旋转，内部互感网络与未绞合完全相同——
        本项目 Q2/Q3 的核心物理结论「rigid rotation」），整体旋转是规范自由度，
        判据必须对其免疫 ⇒ 刚性旋转不变性成为可断言的定理（自检 3）。
    "lab" —— 实验室系原样计算。横向外场通道会奖励任何 coherent 旋转（包括
        单绞向）——作为「外部场通道」补充指标单独报告，不作 headline，否则
        会给单绞向记上它并不具备的（内部损耗意义上的）换位功劳。

合成     D = w_dwell·D_dwell + w_phi·D_phi（默认各 0.5，可配；两指标均∈[0,1]）。

自检（run_selfcheck，任务书：同轨迹 D=0；随机轨迹接近上界）
  1) 同轨迹 → D=0：所有丝共走同一条【面积匹配】路径（各壳层停留占比=面积
     权重 q）⇒ D_dwell=0 且各丝 F_b 全同 ⇒ D_phi=0；
  2) 随机轨迹接近上界：随机槽位分配后【冻结】（未绞合随机束）⇒ p_i=δ_{ℓ(i)}，
     mean_i TV(δ,q)=1−Σq²≈0.875，D_phi≈对径上界 ⇒ D≈0.85（上界口径）。
     注意辨析：若让各丝独立随机跳槽（iid 混合），大数定律同样均等化、D 会
     很小——这不是判据失效，而是连续层只度量均等化；iid 随机辫的失败在
     拓扑层（无周期闭合、w^p 非纯辫），由 braid_check.py 负责，两层互补；
  3) 刚性旋转不变性：对任意逐站整体旋转 ψ_s（可>2π、非常速），derotated
     框架 D 严格不变（与静态逐位相等）；lab 框架则被「虚假改善」——两者
     对照即「单绞向绞合无换位收益」在判据层的验证；
  4) 轮转 round-robin（每丝均匀遍历全部槽位）→ D_phi=0 精确成立，
     D_dwell=布局地板 TV(槽位均匀分布, q)——离散槽位集的可达下界。

接口：pos (S,N,2) 为逐站截面坐标 (x,y)。sample_positions(traj, period, S)
     从轨迹对象采样（geo.trajectories 的 positions(z)->(N,2) 约定）；
     geo.trajectories.radial_history 的 (S,N,2)=[r,θ] 可取 [...,0] 喂 radii
     （仅够 D_dwell）。
"""
from __future__ import annotations

import math

import numpy as np

BASIS_NAMES = ("r", "r cosθ", "r sinθ", "r²", "r²cos2θ", "r²sin2θ")
DEFAULT_PHI_WEIGHTS = {name: 1.0 / 6 for name in BASIS_NAMES}


# ---------------------------------------------------------------- 基函数与框架
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


def derotate_stations(pos, return_angles=False, coh_threshold=0.7):
    """逐站去除整体刚性旋转（Kabsch 单自由度，对齐到站 0）。

    ψ_s = arg Σ_k z_s,k·conj(z_0,k)，去旋后坐标 = e^{−iψ_s}·z_s。
    刚性旋转轨迹 z_s = e^{iφ_s}z_0 ⇒ ψ_s=φ_s 精确还原静态构型（机器精度）。
    旋转不变量（r、r² 基、dwell）不受影响；只有角矩被规范化。

    相干门 coh_threshold：|Σ z_s conj(z_0)| / Σ|z_s||z_0| < 阈值时视该站
    「无相干整体旋转」（如纯重标记/轮转调度——同一槽位集换标签，相干度
    ~1/√N 的噪声），不施加去旋，避免把噪声角当作旋转去除（污染 D_phi）。
    刚性旋转、双层环整体轮转（等半径 e^{iΔ}）、内外互换（同角度实数倍
    半径互换）的相干度均为 1，正常去除。
    """
    pos = np.asarray(pos, dtype=float)
    z_ref = pos[0, :, 0] + 1j * pos[0, :, 1]
    S = pos.shape[0]
    out = np.empty_like(pos)
    angles = np.zeros(S)
    denom_ref = np.abs(z_ref)
    for s in range(S):
        z_s = pos[s, :, 0] + 1j * pos[s, :, 1]
        coh = np.sum(z_s * np.conj(z_ref))
        denom = float(np.sum(np.abs(z_s) * denom_ref))
        psi = 0.0
        if denom > 0 and abs(coh) >= coh_threshold * denom:
            psi = float(np.angle(coh))
        rot = np.exp(-1j * psi) * z_s
        out[s, :, 0], out[s, :, 1] = rot.real, rot.imag
        angles[s] = psi
    if return_angles:
        return out, angles
    return out


def phi_dispersion(pos, weights=None, frame="derotated"):
    """D_phi：各基跨丝周期平均离散度 dict + 加权平均。

    每基归一化 range_k F_b(k) / (2·max_{s,k}|φ_b|)：对径静态=1、全同=0。
    frame="derotated" 先去整体旋转（内部损耗通道，默认）；"lab" 实验室系。
    """
    pos = np.asarray(pos, dtype=float)
    if frame == "derotated":
        pos = derotate_stations(pos)
    elif frame != "lab":
        raise ValueError(f"unknown frame: {frame}")
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


# ---------------------------------------------------------------- D_dwell
def _tv(p, q):
    return float(0.5 * np.abs(np.asarray(p) - np.asarray(q)).sum())


def dwell_deviation(pos=None, radii=None, n_shells=8, bundle_r=None):
    """D_dwell：停留分布 vs 面积加权零假设 q。

    输入：pos (S,N,2) 或 radii (S,N)（radial_history 输出取 [...,0]）。
    返回 dict：D_dwell（headline=mean_i TV(p_i,q)，完美=0）、
    max_tv（最差丝）、fairness（mean_i TV(p_i,p̄)，纯跨丝公平）、
    area_bias（TV(p̄,q)，槽位集覆盖属性）、shell_edges、q、p_mean。
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
    q = (edges[1:] ** 2 - edges[:-1] ** 2) / bundle_r ** 2   # 面积权重零假设
    return {
        "D_dwell": float(np.mean([_tv(p[i], q) for i in range(N)])),
        "max_tv": float(np.max([_tv(p[i], q) for i in range(N)])),
        "fairness": float(np.mean([_tv(p[i], p_mean) for i in range(N)])),
        "area_bias": _tv(p_mean, q),
        "shell_edges": edges.tolist(), "q": q.tolist(),
        "p_mean": p_mean.tolist(),
    }


# ---------------------------------------------------------------- 合成
def criterion(pos=None, radii=None, n_shells=8, w_dwell=0.5, w_phi=0.5,
              phi_weights=None, frame="derotated"):
    """合成判据 D 及全部分解量。pos 与 radii 至少给一个（D_phi 需要 pos）。"""
    out = {"n_shells": n_shells, "w_dwell": w_dwell, "w_phi": w_phi,
           "frame": frame}
    dw = dwell_deviation(pos=pos, radii=radii, n_shells=n_shells)
    out.update({k: dw[k] for k in
                ("D_dwell", "max_tv", "fairness", "area_bias", "q", "p_mean")})
    if pos is not None:
        ph = phi_dispersion(pos, weights=phi_weights, frame=frame)
        out["phi_per_basis"] = ph["per_basis"]
        out["D_phi"] = ph["D_phi"]
    else:
        out["phi_per_basis"] = None
        out["D_phi"] = None
        w_phi = 0.0
    out["D"] = w_dwell * out["D_dwell"] + w_phi * (out["D_phi"] or 0.0)
    return out


def sample_positions(traj, period_len, n_samples=240):
    """全信息采样：轨迹对象 → (S,N,2) 逐站截面坐标（喂 criterion 的 pos）。"""
    zs = np.linspace(0.0, period_len, n_samples)
    return np.array([np.asarray(traj(z), dtype=float) for z in zs])


# ---------------------------------------------------------------- 自检用例
def _uniform_disc(n, r=1.0, seed=11):
    rng = np.random.default_rng(seed)
    rr = r * np.sqrt(rng.random(n))
    tt = 2 * math.pi * rng.random(n)
    return np.stack([rr * np.cos(tt), rr * np.sin(tt)], axis=1)


def _frozen_random(slots, S=64):
    """随机轨迹（任务书口径）：随机槽位分配后冻结 = 未绞合随机束。"""
    return np.repeat(slots[None, :, :], S, axis=0)


def _iid_random(slots, S=240, seed=7):
    """各丝独立随机跳槽（iid 混合）——LLN 会均等化，作对照辨析用。"""
    rng = np.random.default_rng(seed)
    n = len(slots)
    return slots[rng.integers(0, n, size=(S, n))]


def _round_robin(slots, S=None):
    n = len(slots)
    S = S or n
    return np.array([slots[(np.arange(n) + s) % n] for s in range(S)])


def _covering_slots(n_sh=8, r=1.0):
    """每个壳层都有槽位的确定性槽位集（同轨迹 D=0 自检用）：壳层 ℓ 放
    2(ℓ+1) 个槽位于壳层中径、角度均布 → 任意壳层非空，面积匹配路径可行。"""
    slots = []
    for l in range(n_sh):
        rr = r * (l + 0.5) / n_sh
        for j in range(2 * (l + 1)):
            a = 2 * math.pi * j / (2 * (l + 1)) + l * 0.37
            slots.append((rr * math.cos(a), rr * math.sin(a)))
    return np.array(slots)


def _area_matched_path(slots, n_sh=8, r=1.0, S=64):
    """构造「面积匹配」路径：各壳层停留站数 = 面积权重×S（同轨迹 D=0 用）。

    等宽壳层 q(ℓ)=(2ℓ+1)/n_sh²，取 S=n_sh² ⇒ 计数 (2ℓ+1) 为精确整数，
    每壳层在其槽位间轮转 → 该路径的停留分布 p≡q 精确成立。
    """
    edges = np.linspace(0.0, r, n_sh + 1)
    radii = np.hypot(slots[:, 0], slots[:, 1])
    idx = np.clip(np.digitize(radii, edges[1:-1], right=True), 0, n_sh - 1)
    groups = [np.where(idx == l)[0] for l in range(n_sh)]
    counts = np.array([2 * l + 1 for l in range(n_sh)])
    assert sum(counts) == S and all(len(g) > 0 for g in groups)
    order = []
    ptr = [0] * n_sh
    # 轮流填充：按壳层次序每轮消耗计数，保证任意前缀无病态偏置
    while len(order) < S:
        for l in range(n_sh):
            if counts[l] > 0:
                g = groups[l]
                order.append(slots[g[ptr[l] % len(g)]])
                ptr[l] += 1
                counts[l] -= 1
                if len(order) == S:
                    break
    return np.array(order)


def _rigid_rotation_traj(slots, psis):
    """整体刚性旋转轨迹：站 s 全体坐标旋转 ψ_s（可>2π、非常速）。"""
    out = np.empty((len(psis), len(slots), 2))
    for s, psi in enumerate(psis):
        c, sn = math.cos(psi), math.sin(psi)
        out[s, :, 0] = c * slots[:, 0] - sn * slots[:, 1]
        out[s, :, 1] = sn * slots[:, 0] + c * slots[:, 1]
    return out


def _fmt(res):
    phi = res.get("D_phi")
    phi_s = f"{phi:.4f}" if phi is not None else "  n/a "
    return (f"D={res['D']:.4f}  D_dwell={res['D_dwell']:.4f} "
            f"(max {res['max_tv']:.3f}, fair {res['fairness']:.3f}, "
            f"area {res['area_bias']:.3f})  D_phi={phi_s}")


def run_selfcheck():
    print("== criterion.py self-check (定稿版) ==")
    slots = _uniform_disc(42)
    slots_cov = _covering_slots()

    # 1) 同轨迹（面积匹配路径）→ D=0
    path = _area_matched_path(slots_cov, n_sh=8, r=1.0, S=64)
    same = np.repeat(path[:, None, :], 6, axis=1)      # 6 丝共走同一条路径
    res_same = criterion(pos=same, n_shells=8)

    # 2) 随机轨迹（随机分配冻结）→ 接近上界
    res_rand = criterion(pos=_frozen_random(slots))

    # 2b) iid 混合对照（辨析：LLN 均等化，失败在拓扑层）
    res_iid = criterion(pos=_iid_random(slots))

    # 3) 刚性旋转不变性（derotated 定理 vs lab 虚假改善）
    psis = 2.71 * np.arange(64) + 0.3 * np.sin(np.arange(64))
    rot = _rigid_rotation_traj(slots, psis)
    res_rot = criterion(pos=rot)
    res_rot_lab = criterion(pos=rot, frame="lab")
    res_stat_lab = criterion(pos=_frozen_random(slots), frame="lab")

    # 4) 轮转 round-robin → D_phi=0、D_dwell=布局地板
    rr = _round_robin(slots, S=210)                    # 210=42×5 整周期
    res_rr = criterion(pos=rr)

    cases = [("same-traj (area-matched path)", res_same),
             ("random frozen (untransposed)", res_rand),
             ("random iid mixing (LLN, topo layer!)", res_iid),
             ("rigid rotation, derotated", res_rot),
             ("rigid rotation, lab frame", res_rot_lab),
             ("frozen static, lab frame", res_stat_lab),
             ("round-robin (perfect schedule)", res_rr)]
    for k, v in cases:
        print(f"  {k:36s} {_fmt(v)}")

    # ---- 断言 ----
    assert res_same["D"] < 1e-9, "identical area-matched trajectories must give D=0"
    assert res_rand["D"] > 0.7, "random frozen allocation should be near upper bound"
    assert res_rand["D_dwell"] > 0.8 and res_rand["D_phi"] > 0.5
    assert res_iid["D"] < res_rand["D"] - 0.3, "iid mixing equalizes (LLN); " \
        "its failure is topological -> braid_check.py"
    assert abs(res_rot["D"] - res_rand["D"]) < 1e-9, \
        "derotated D must be invariant under arbitrary per-station rigid rotation"
    assert res_rot_lab["D_phi"] < res_stat_lab["D_phi"] - 0.2, \
        "lab frame rewards mere rotation (spurious) -- why derotation is default"
    assert res_rr["D_phi"] < 1e-12, "round-robin equalizes phi-moments exactly"
    assert res_rr["D_dwell"] < res_rand["D_dwell"] - 0.2
    print("  PASS: same-traj D=0; random-frozen near upper bound "
          f"(D={res_rand['D']:.3f}); rigid rotation leaves derotated D unchanged "
          f"({res_rot['D']:.6f} vs {res_rand['D']:.6f}) while lab frame would "
          f"falsely credit it (D_phi {res_stat_lab['D_phi']:.3f}->"
          f"{res_rot_lab['D_phi']:.3f}); round-robin D_phi=0 with layout floor "
          f"D_dwell={res_rr['D_dwell']:.4f}.")
    print("  NOTE: iid random mixing scores LOW (LLN) -- not a criterion bug; "
          "random braid words fail periodicity/purity at the topological layer "
          "(braid_check.py).")
    return {k: v for k, v in cases}


if __name__ == "__main__":
    run_selfcheck()
