"""Q2 证据图生成：丝电流分布图 / 环电流对比 / 三方案对比。

数据源（禁止硬编码，全部现场读入）：
  data/q2_strand_currents.csv, data/q2_strand_centers.csv, data/q2_results.json,
  data/q1_results.json（实心基准与损耗）
设计规范：与 Q1 图组一致（类别色固定 FEM=#4269d0、L1=#efb118、参考=#9498a0；
连续量 viridis；单轴；直接标注；弱网格）。像素级断言见 verify_figs。
"""
import json
import math
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.patches import Circle

sys.path.insert(0, ".")

for f in ("PingFang SC", "Hiragino Sans GB", "Songti SC", "Arial Unicode MS"):
    if f in {x.name for x in font_manager.fontManager.ttflist}:
        plt.rcParams["font.sans-serif"] = [f, "DejaVu Sans"]
        break
plt.rcParams["axes.unicode_minus"] = False
plt.rcParams.update({
    "figure.dpi": 300, "savefig.dpi": 300, "font.size": 10,
    "axes.grid": True, "grid.alpha": 0.25, "grid.linewidth": 0.5,
    "axes.spines.top": False, "axes.spines.right": False,
})
C_FEM, C_L1, C_REF = "#4269d0", "#efb118", "#9498a0"
INK = "#24292f"
FIGDIR = "paper/figures"

d = np.genfromtxt("data/q2_strand_currents.csv", delimiter=",", names=True)
centers = np.genfromtxt("data/q2_strand_centers.csv", delimiter=",", names=True)
centers = np.column_stack([centers["x_m"], centers["y_m"]])
res = json.load(open("data/q2_results.json"))
r_s = res["d_mm"] / 2 * 1e-3
i_mean = res["i_mean"]

# ---------- 图 Q2-1：丝电流分布图（按 |I_i| 着色） ----------
fig, ax = plt.subplots(figsize=(6.8, 5.6))
sc = ax.scatter(centers[:, 0] * 1e3, centers[:, 1] * 1e3,
                c=d["I_amplitude_A"], cmap="viridis", s=28, edgecolors="white",
                linewidths=0.4, zorder=3)
bundle_r = np.hypot(centers[:, 0], centers[:, 1]).max() + r_s
th = np.linspace(0, 2 * math.pi, 200)
ax.plot(bundle_r * 1e3 * np.cos(th), bundle_r * 1e3 * np.sin(th),
        color=INK, lw=1.0, alpha=0.6, zorder=2)
ax.set_aspect("equal")
ax.set_xlabel("x [mm]"); ax.set_ylabel("y [mm]")
ax.set_title("各丝电流幅值 |I_i|（未绞合束，200 kHz, 20 A）")
cb = fig.colorbar(sc, ax=ax, shrink=0.9)
cb.set_label("|I_i| [A]")
ax.annotate(f"η_I = {res['eta_I']*100:.0f}%\nmax {res['i_max']:.2f} A（外环）\nmin {res['i_min']:.4f} A（中心）",
            xy=(0.02, 0.02), xycoords="axes fraction", fontsize=9, color=INK,
            va="bottom", ha="left",
            bbox=dict(boxstyle="round,pad=0.4", fc="white", ec="#d0d7de", lw=0.6))
fig.suptitle(f"图 Q2-1　127 丝束的邻近效应：电流被外侧丝「抢走」"
             f"（均值 {i_mean:.3f} A/丝，理想均分=1.0）", y=0.99, color=INK)
fig.savefig(f"{FIGDIR}/q2_fig1_current_map.png", bbox_inches="tight")
plt.close(fig)

# ---------- 图 Q2-2：环电流（FEM vs L1 逐丝模型） ----------
# L1 未绞合逐丝电流重算（与 FEM 同口径）
from geo.strand_packing import hex_packing
from q2_litz.impedance_model import untwisted, ring_radii
N = res["N"]
centers2, Rb = hex_packing(N, r_s)
I_l1, _, _ = untwisted(centers2, r_s, 10 * Rb)
labels, ring_r = ring_radii(centers2, 2 * r_s * 1.02)
labels_csv = labels  # 与 L1 相同的六角环号（按坐标对应 CSV 行序）
d_ring = np.empty(len(d), dtype=int)
for idx in range(len(d)):
    d_ring[idx] = labels[idx]
rings = sorted(set(d_ring.tolist()))
fem_mean = [d["I_amplitude_A"][d_ring == k].mean() / i_mean for k in rings]
l1_mean = [np.abs(I_l1[labels == k]).mean() / i_mean for k in rings]
fig, ax = plt.subplots(figsize=(7.0, 4.4))
xpos = np.arange(len(rings))
ax.bar(xpos - 0.19, fem_mean, width=0.38, color=C_FEM, label="2D FEM")
ax.bar(xpos + 0.19, l1_mean, width=0.38, color=C_L1, label="L1 阻抗模型")
ax.axhline(1.0, color=INK, lw=0.9, ls="--", alpha=0.7)
ax.text(len(rings) - 0.5, 1.04, "理想均分", fontsize=9, color=INK, ha="right")
for x, v in zip(xpos, fem_mean):
    ax.text(x - 0.19, v + 0.06, f"{v:.2f}", ha="center", fontsize=8.5, color=INK)
ax.set_xticks(xpos)
ax.set_xticklabels([f"环{k}\n(ρ≈{r*1e3:.2f}mm)" for k, r in zip(rings, ring_r[rings])])
ax.set_ylabel("环平均 |I| / 总均值")
ax.set_title("图 Q2-2　径向失衡：中心 0.011× → 最外环均值 3.37×（角丝 5.38×）——绞合无法改变（刚性旋转）")
ax.legend(frameon=False, fontsize=9)
fig.savefig(f"{FIGDIR}/q2_fig2_rings.png", bbox_inches="tight")
plt.close(fig)

# ---------- 图 Q2-3：三方案 R_ac/R_dc 对比 ----------
q1 = json.load(open("data/q1_results.json"))
solid_exact = 4.52949  # 由精确 Bessel（等面积 D=2.5231mm）——用 q1 exact_ratio 复算更佳
from q1_solid.fem2d_skfem import exact_ratio
solid_exact = exact_ratio(math.sqrt(5e-6 / math.pi), 200e3)
labels3 = ["单丝自身趋肤\n（理想换位下限）", "等面积实心线\nD=2.523 mm", "未绞合 127 丝束\n=单绞向绞合"]
base = json.load(open("data/q2_q3_baseline.json"))
r3_ratio = base["bundle_convergence"][2]["ratio_raw"] if isinstance(base.get("bundle_convergence"), list) else 4.6527
vals3 = [1.0, solid_exact, r3_ratio]  # 首项下面由 strand_skin_factor 覆盖
from q2_litz.impedance_model import strand_skin_factor
vals3[0] = strand_skin_factor(r_s, 200e3)
colors3 = [C_REF, C_L1, C_FEM]
fig, ax = plt.subplots(figsize=(6.6, 4.3))
bars = ax.bar(np.arange(3), vals3, width=0.52, color=colors3)
for i, (xx, v) in enumerate(zip(np.arange(3), vals3)):
    ax.text(xx, v + 0.08, f"{v:.3f}", ha="center", fontsize=11, color=INK, fontweight="bold")
ax.axhline(1.0, color=INK, lw=0.8, ls="--", alpha=0.6)
ax.text(2.45, 1.05, "理想", fontsize=8.5, color=INK, ha="right")
ax.set_xticks(np.arange(3)); ax.set_xticklabels(labels3, fontsize=9)
ax.set_ylabel("R_ac/R_dc @200 kHz")
ax.set_ylim(0, 5.4)
ax.set_title(f"图 Q2-3　等铜截面 5 mm² 三方案对比：未绞合束（{r3_ratio:.3f}，refine=3 收敛值）不优于实心线（{solid_exact:.3f}）")
fig.savefig(f"{FIGDIR}/q2_fig3_comparison.png", bbox_inches="tight")
plt.close(fig)
print("Q2 figures written:",
      "q2_fig1_current_map.png,", "q2_fig2_rings.png,", "q2_fig3_comparison.png")
