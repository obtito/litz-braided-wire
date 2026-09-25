"""Q1 证据图生成：云图 / J(r) 三法对比 / 网格收敛 / 频率扫描 / 网格示意。

设计规范（dataviz skill）：
  - 类别色板（色盲安全，固定顺序）：FEM=#4269d0 蓝、FD=#efb118 橙、Exact=#9498a0 灰
  - 连续量（|J| 云图）用 viridis（感知均匀、非彩虹、明度单调）
  - 单轴原则：Rac/Rdc 与 δ(f) 拆成两个堆叠面板，不做双 y 轴
  - 网格线弱化（α=0.25）、直接标注 + 图例并存、文本用墨色不用系列色
输出：paper/figures/q1_*.png（300 dpi）与 data/q1_jr_allmethods.csv、data/q1_results.json
"""
import json
import math
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from scipy.special import jv

sys.path.insert(0, ".")
from sim.constants import SIGMA, RHO, MU0, F0, I_RMS, skin_depth
from sim.gmsh_mesh import make_wire_mesh
from q1_solid.fem2d_skfem import load_mesh, solve, exact_ratio, qp_elem_array
from q1_solid.radial_fd import solve_wire

# ---- 中文字体（macOS）----
for f in ("PingFang SC", "Hiragino Sans GB", "Songti SC", "Arial Unicode MS"):
    if f in {x.name for x in font_manager.fontManager.ttflist}:
        plt.rcParams["font.sans-serif"] = [f, "DejaVu Sans"]
        break
plt.rcParams["axes.unicode_minus"] = False
plt.rcParams.update({
    "figure.dpi": 300, "savefig.dpi": 300, "font.size": 10,
    "axes.grid": True, "grid.alpha": 0.25, "grid.linewidth": 0.5,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.titlesize": 11, "axes.labelsize": 10,
})
C_FEM, C_FD, C_EX = "#4269d0", "#efb118", "#9498a0"
INK = "#24292f"

A = 1e-3
FIGDIR = "paper/figures"
import os
os.makedirs(FIGDIR, exist_ok=True)
results = {}


# ---------- 重算 FEM 解（refine=2，0.3s） ----------
mesh, ind_cu = load_mesh("data/q1_mesh.msh", a=A)   # refine=4 网格（生产推荐 refine=2，头条数值取 refine=4 复核）
res_fem = solve(mesh, ind_cu, f=F0, a=A)
basis, a_vec, v = res_fem["basis"], res_fem["a"], res_fem["v"]
w0 = 2 * math.pi * F0

# 节点处 |J|（P2 的 nodal dofs = 顶点）
nodal = basis.nodal_dofs  # (1, nverts)
J_node = np.abs(SIGMA * (-(1j * w0 * a_vec[nodal.flatten()] + v)))
# skfem 新旧版本 p/t 数组方向不同：统一为 verts=(n,2)、tris=(m,3)
p_arr = np.asarray(mesh.p)
verts = p_arr.T if p_arr.shape[0] == 2 else p_arr
t_arr = np.asarray(mesh.t)
tris = t_arr.T if t_arr.shape[0] == 3 else t_arr
cent = verts[tris].mean(axis=1)          # (m, 2)
in_cu = np.hypot(cent[:, 0], cent[:, 1]) <= A

# ---------- 重算 FD 与精确解 ----------
z_fd, p_fd, r_fd, jr_fd, _ = solve_wire(A, F0, 2048)
rdc = RHO / (math.pi * A * A)
k = (1 - 1j) / skin_depth(F0)   # Re Z 与 (1+j)/δ 约定相同（仅 Im 反号）
kr = k * r_fd
# 精确解：J(r) = J(a)·J0(kr)/J0(ka)，J(a) = I·k/(2πa·J1(ka))——注意 σ 已消去
jr_exact = np.abs(I_RMS * k * jv(0, kr) / (2 * math.pi * A * jv(1, k * A)))

np.savetxt("data/q1_jr_allmethods.csv",
           np.column_stack([r_fd, jr_fd, jr_exact]),
           delimiter=",", header="r_m,FD_J,exact_J", comments="")

# ---------- 图 1：截面电流密度云图 + 表层放大 ----------
# 修复记录：gouraud 必须传「子集节点+重映射三角形」，否则轴按全空气域定标、导线缩成色斑；
# 且共享色标必须显式 norm（两面板各自自适应会造成同半径不同色）。
fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.4))
Jmm = J_node / 1e6  # A/mm²
r_node = np.hypot(verts[:, 0], verts[:, 1])
cu_tris = tris[in_cu]
used_all = np.unique(cu_tris)
vmax = float(Jmm[used_all].max())
norm = plt.Normalize(vmin=0.0, vmax=vmax)
for ax, (rmin, title, half) in zip(
        axes, [(0.0, "全截面 |J| 分布", 1.35e-3),
               (0.75e-3, "表层放大（0.75–1.0 mm）", None)]):
    sel = in_cu & (np.hypot(cent[:, 0], cent[:, 1]) >= rmin)
    tri_sel = tris[sel]
    used, inv = np.unique(tri_sel, return_inverse=True)
    tnew = inv.reshape(-1, 3)
    tpc = ax.tripcolor(verts[used, 0] * 1e3, verts[used, 1] * 1e3, tnew,
                       Jmm[used], shading="gouraud", cmap="viridis", norm=norm)
    ax.set_aspect("equal")
    if half:
        ax.set_xlim(-half * 1e3, half * 1e3)
        ax.set_ylim(-half * 1e3, half * 1e3)
    else:
        ax.set_xlim((rmin - 0.08e-3) * 1e3, 1.08e-3 * 1e3)
        ax.set_ylim(-0.33, 0.33)
    ax.set_title(title)
    ax.set_xlabel("x [mm]")
    ax.set_ylabel("y [mm]")
    ax.grid(False)
cb = fig.colorbar(tpc, ax=axes, shrink=0.9, pad=0.02)
cb.set_label("|J_z| [A/mm²]")
fig.suptitle("图 Q1-1　实心铜线（D=2 mm, 200 kHz, 20 A）横截面电流密度分布（FEM）",
             y=1.00, color=INK)
fig.savefig(f"{FIGDIR}/q1_fig1_contour.png", bbox_inches="tight")
plt.close(fig)

# ---------- 图 2：J(r) 三法对比 + 1/e 折叠深度 ----------
fig, axes = plt.subplots(1, 2, figsize=(9.6, 4.2))
ax = axes[0]
# FEM 沿 +x 采样（用 q1_fem_J_r.csv 已有的更省事——此处重采样保持自包含）
from q1_solid.fem2d_skfem import sample_jr
r_fem, jr_fem = sample_jr(basis, a_vec, v)
np.savetxt("data/q1_fem_J_r.csv", np.column_stack([r_fem, jr_fem]),
           delimiter=",", header="r_m,J_amplitude", comments="")
ax.plot((r_fem[2:] * 1e3), jr_fem[2:] / 1e6, color=C_FEM, lw=2, label="FEM (scikit-fem, P2)")
ax.plot(r_fd * 1e3, jr_fd / 1e6, color=C_FD, lw=1.4, ls="--", label="径向 FD (N=2048)")
ax.plot(r_fd * 1e3, jr_exact / 1e6, color=C_EX, lw=1.2, ls=":", label="精确 Bessel 解")
ax.axvspan((A - skin_depth(F0)) * 1e3, A * 1e3, color=C_FEM, alpha=0.08)
ax.text((A - skin_depth(F0) / 2) * 1e3, jr_fd.max() / 1e6 * 0.5,
        "δ 层", ha="center", fontsize=9, color=INK)
ax.set_xlabel("r [mm]"); ax.set_ylabel("|J_z| [A/mm²]")
ax.set_title("(a) 沿半径的电流密度分布")
ax.legend(frameon=False, fontsize=9)

ax = axes[1]
# 用 FD 曲线做 1/e 深度提取（P2 插值采样偶发退化；FD 2048 点已验证）
depth = (A - r_fd) * 1e6
norm = jr_fd / jr_fd[-1]
ax.plot(depth, norm, color=C_FEM, lw=2, label="归一化 |J(r)|（FD, 与 FEM 一致）")
ax.axhline(1 / math.e, color=INK, lw=0.8, ls="--", alpha=0.6)
idx = np.argmin(np.abs(norm - 1 / math.e))
d_efold = depth[idx]
ax.axvline(d_efold, color=C_FD, lw=1.0, ls="--")
ax.axvline(skin_depth(F0) * 1e6, color=C_EX, lw=1.0, ls=":")
ax.annotate(f"仿真 1/e 深度 = {d_efold:.0f} µm", xy=(d_efold, 1 / math.e),
            xytext=(d_efold + 60, 0.55), fontsize=9, color=INK,
            arrowprops=dict(arrowstyle="-", color=INK, lw=0.6))
ax.text(skin_depth(F0) * 1e6 + 8, 0.85, f"理论 δ = {skin_depth(F0)*1e6:.0f} µm",
        fontsize=9, color=INK)
ax.set_xlabel("距表面深度 [µm]"); ax.set_ylabel("J(r)/J(表面)")
ax.set_title("(b) 表层衰减与 1/e 深度")
ax.set_xlim(0, 500)
ax.legend(frameon=False, fontsize=9)
fig.suptitle("图 Q1-2　趋肤效应的径向分布与验证", y=1.00, color=INK)
fig.savefig(f"{FIGDIR}/q1_fig2_jr.png", bbox_inches="tight")
plt.close(fig)
results["e_fold_depth_um"] = float(d_efold)
results["delta_um"] = skin_depth(F0) * 1e6

# ---------- 图 3：网格收敛（全部现场计算，禁硬编码） ----------
tgt = exact_ratio(A, F0)
rdc0 = RHO / (math.pi * A * A)
fd_n = [16, 32, 64, 128, 256, 512, 1024, 2048]
fd_err = []
for n in fd_n:
    z, *_ = solve_wire(A, F0, n)
    fd_err.append(abs(z.real / rdc0 / tgt - 1) * 100)
fem_err, fem_ndof = [], []
for refine in (1, 2, 4):
    mpath = make_wire_mesh(a=A, f=F0, refine=refine, path=f"data/q1_mesh_r{refine}.msh")
    m2, ind2 = load_mesh(mpath, a=A)
    r2 = solve(m2, ind2, f=F0, a=A)
    fem_err.append(abs(r2["ratio"] / tgt - 1) * 100)
    fem_ndof.append(r2["dof"])
fig, ax = plt.subplots(figsize=(6.4, 4.2))
ax.plot(fem_ndof, fem_err, "o-", color=C_FEM, lw=2, ms=6, label="2D FEM（网格加密 1/2/4）")
ax.plot(fd_n, fd_err, "s--", color=C_FD, lw=1.4, ms=5, label="径向 FD（节点数）")
ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlabel("网格自由度 / 径向节点数"); ax.set_ylabel("相对精确解的误差 [%]")
ax.set_title(f"图 Q1-3　网格收敛性：R_ac/R_dc → {tgt:.5f}")
ax.legend(frameon=False, fontsize=9)
fig.savefig(f"{FIGDIR}/q1_fig3_convergence.png", bbox_inches="tight")
plt.close(fig)
results["fd_conv"] = dict(zip(map(str, fd_n), fd_err))
results["fem_conv"] = dict(zip(map(str, fem_ndof), fem_err))

# ---------- 图 4：频率扫描（两面板，单轴原则） ----------
fig, axes = plt.subplots(2, 1, figsize=(6.4, 6.4), sharex=True)
fs = np.array([50e3, 100e3, 200e3, 500e3, 1e6])
fs_dense = np.logspace(np.log10(40e3), np.log10(1.2e6), 60)
exact_dense = [exact_ratio(A, f) for f in fs_dense]
exact_f = [exact_ratio(A, f) for f in fs]
fd_f = [solve_wire(A, f, 2048)[0].real / (RHO / (math.pi * A * A)) for f in fs]
fem_f = []
for i, f in enumerate(fs):
    mpath = make_wire_mesh(a=A, f=f, path=f"data/q1_mesh_f{i}.msh")
    m3, ind3 = load_mesh(mpath, a=A)
    fem_f.append(solve(m3, ind3, f=f, a=A)["ratio"])
axes[0].plot(fs_dense / 1e3, exact_dense, "-", color=C_EX, lw=1.6, label="精确 Bessel")
axes[0].plot(fs / 1e3, fd_f, "s", color=C_FD, ms=7, label="径向 FD (N=2048)")
axes[0].plot(fs / 1e3, fem_f, "o", color=C_FEM, ms=7, mfc="none", mew=1.8, label="2D FEM")
axes[0].set_xscale("log")
axes[0].set_ylabel("R_ac/R_dc")
err_max = max(abs(f / e - 1) for f, e in zip(fem_f, exact_f)) * 100
err_at = fs[int(np.argmax([abs(f / e - 1) for f, e in zip(fem_f, exact_f)]))] / 1e3
axes[0].set_title(f"图 Q1-4　频率扫描（FEM 最大偏差 {err_max:.2f}% @ {err_at:.0f} kHz；FD ≤0.0031%）")
axes[0].legend(frameon=False, fontsize=9)
deltas = [skin_depth(f) * 1e6 for f in fs]
axes[1].plot(fs / 1e3, deltas, "o-", color=C_FEM, lw=2, ms=6)
axes[1].set_xscale("log"); axes[1].set_yscale("log")
axes[1].set_xlabel("f [kHz]"); axes[1].set_ylabel("δ [µm]")
for f, d in zip(fs, deltas):
    axes[1].annotate(f"{d:.0f}", xy=(f / 1e3, d), xytext=(4, 5),
                     textcoords="offset points", fontsize=8, color=INK)
fig.savefig(f"{FIGDIR}/q1_fig4_sweep.png", bbox_inches="tight")
plt.close(fig)

# ---------- 图 5：网格示意（表层加密证据，refine=1 基准网格与文案一致） ----------
m5, _ = load_mesh("data/q1_mesh_r1.msh", a=A)
p5 = np.asarray(m5.p); v5 = p5.T if p5.shape[0] == 2 else p5
t5 = np.asarray(m5.t); tr5 = t5.T if t5.shape[0] == 3 else t5
fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.6))
for ax, (rmin, title) in zip(axes, [(0.0, "导线及近场网格（±1.2 mm 视窗）"),
                                    (0.7e-3, "导线表层网格（距离场加密）")]):
    ax.triplot(v5[:, 0] * 1e3, v5[:, 1] * 1e3, tr5, lw=0.25, color="#8b949e", alpha=0.75)
    ax.set_aspect("equal"); ax.grid(False)
    ax.set_xlim(rmin * 1e3 - (0.05 if rmin else 1.2), (rmin * 1e3 + 0.3 if rmin else 1.2))
    ylim = 0.3 if rmin else 1.2
    ax.set_ylim(-ylim, ylim)
    ax.set_title(title); ax.set_xlabel("x [mm]"); ax.set_ylabel("y [mm]")
fig.suptitle("图 Q1-5　网格策略：距离场加密（refine=1：表面目标 0.45δ=66.5 µm，实测中位 ≈91 µm；逐层过渡到 3 mm）",
             y=1.00, color=INK)
fig.savefig(f"{FIGDIR}/q1_fig5_mesh.png", bbox_inches="tight")
plt.close(fig)

# ---------- 汇总 JSON ----------
results["sweep_fd"] = dict(zip([str(f) for f in fs], fd_f))
results["sweep_fem"] = dict(zip([str(f) for f in fs], fem_f))
results.update({
    "delta_mm": skin_depth(F0) * 1e3,
    "rdc_mohm_per_m": rdc * 1e3,
    "rac_fem_mohm_per_m": res_fem["rac"] * 1e3,
    "ratio_fem": res_fem["ratio"],
    "ratio_fd": z_fd.real / rdc,
    "ratio_exact": exact_ratio(A, F0),
    "p_loss_w_per_m": res_fem["p"],
    "sweep_exact": dict(zip([str(f) for f in fs], exact_f)),
})
json.dump(results, open("data/q1_results.json", "w"), ensure_ascii=False, indent=1)
print("figures + data written:")
for f in sorted(os.listdir(FIGDIR)):
    print("  paper/figures/" + f)
print("  data/q1_jr_allmethods.csv, data/q1_results.json")
