#!/usr/bin/env python
"""生成 Q4 展示网站的 demo 轨迹 JSON（写入 web/public/demo/）。

用法（在仓库任意位置）:
    ./.venv/bin/python web/tools/gen_demo.py

demo1 单绞向（对照）: r(z) 逐丝恒定、截面刚体旋转——RESEARCH.md §三/§七的
    「简单成绞 = 刚体旋转」，无任何径向换位，φ-矩判据 D 很大。
demo2 轮转+翻面（换位示意）: θ 均匀轮转（方位均匀）+ 每丝相位错开 2π/N 的
    周期性内外壳层交换（径向遍历）。r(z) 用 tanh 展平的正弦实现
    「壳层停留 + 平滑过渡」（RESEARCH.md §7.2 的轨迹观点）。
    这是示意性连续模型——q4_perfect 产出角齿轮事件序列的真实轨迹后，
    直接替换同名 JSON 即可，查看器零改动。

输出格式（查看器契约，见 web/src/braid.js）:
    { name, title, note, N, units, period_len, z_total, n_samples,
      strand_r, r_min, r_max,
      strands: [ { r: [...], theta: [...] } ],   # theta 为【展开】弧度
      metrics: { basis, mean_r, D } }

坐标单位 mm；z 均匀采样含端点；theta 单调不取模，供 Catmull-Rom 插值。
metrics.D = φ-矩判据（RESEARCH.md §7.1 离散层）：基 Φ={r, r²}，
D = max_{i,j}|ΣΦ_i − ΣΦ_j| / max_k|ΣΦ_k|，在一个周期上等权采样计算。
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np

OUT_DIR = Path(__file__).resolve().parents[1] / "public" / "demo"

# ---------- 几何参数（mm；与 Q2 设计同量级：束径 ~3mm、丝径 ~0.2mm） ----------
P = 24.0          # 换位周期（demo1 一个整旋转节距 / demo2 一个完整换位周期）
N_PERIODS = 2     # 展示两个周期，首尾截面严格相同（周期性）
SAMPLES_PER_PERIOD = 240


def sample_z():
    n_total = SAMPLES_PER_PERIOD * N_PERIODS + 1
    return np.linspace(0.0, P * N_PERIODS, n_total)


# ---------- demo1：单绞向（恒半径螺旋，刚体旋转） ----------
def demo1_helix() -> dict:
    N, strand_r = 19, 0.14
    radii = np.array([0.0] + [0.70] * 6 + [1.40] * 12)
    theta0 = np.array(
        [0.0]
        + [j * math.pi / 3 for j in range(6)]
        + [j * math.pi / 6 + math.pi / 12 for j in range(12)]
    )  # 1 心丝 + 内六角 + 外十二角（microbundle 双圈排布，见 geo/trajectories.py）
    z = sample_z()
    omega = 2 * math.pi / P
    strands = [
        {"r": np.full_like(z, radii[k]), "theta": theta0[k] + omega * z}
        for k in range(N)
    ]
    return {
        "name": "demo1-helix",
        "title": "单绞向（对照）",
        "note": "r(z) 恒定、截面刚体旋转：丝的径向位置终生不变——与未绞合束电磁等价 "
                "（findings.md H+8 物理勘误）。R_ac/R_dc ≈ 未绞合 4.65，比等截面实心还差 2.75%。",
        "N": N,
        "units": "mm",
        "period_len": P,
        "z_total": P * N_PERIODS,
        "n_samples": len(z),
        "strand_r": strand_r,
        "r_min": float(radii.min()),
        "r_max": float(radii.max()),
        "strands": [pack(s) for s in strands],
        **metrics(strands, z),
    }


# ---------- demo2：轮转+翻面（径向换位示意） ----------
def demo2_transpose() -> dict:
    N, strand_r = 12, 0.13
    r_in, r_out = 0.60, 1.45                 # 内外壳层半径
    r_mid, amp = 0.5 * (r_in + r_out), 0.5 * (r_out - r_in)
    flat = 2.0                               # tanh 展平：壳层停留 + 平滑过渡
    z = sample_z()
    omega = 2 * math.pi / P
    psi = 2 * math.pi * np.arange(N) / N    # 每丝径向相位错开 2π/N：轮流进出壳层
    strands = []
    for k in range(N):
        u = psi[k] - omega * z               # 径向迁移相位（与轮转反向 → 编织交叉可见）
        r = r_mid + amp * np.tanh(flat * np.sin(u)) / math.tanh(flat)
        theta = psi[k] + omega * z           # 整体轮转：方位均匀
        strands.append({"r": r, "theta": theta})
    return {
        "name": "demo2-transpose",
        "title": "轮转+翻面（换位示意）",
        "note": "θ 均匀轮转（方位均匀）+ 相位错开的周期壳层交换（径向遍历）：每根丝轮流"
                "占据所有径向位置 → φ-矩判据 D→0（完美换位，RESEARCH.md §7.1/§7.2）。"
                "示意性连续模型；q4_perfect 的角齿轮事件序列轨迹产出后替换本文件。",
        "N": N,
        "units": "mm",
        "period_len": P,
        "z_total": P * N_PERIODS,
        "n_samples": len(z),
        "strand_r": strand_r,
        "r_min": r_in,
        "r_max": r_out,
        "strands": [pack(s) for s in strands],
        **metrics(strands, z),
    }


# ---------- 工具 ----------
def pack(s: dict) -> dict:
    return {
        "r": [round(float(v), 5) for v in s["r"]],
        "theta": [round(float(v), 5) for v in s["theta"]],
    }


def metrics(strands: list[dict], z: np.ndarray) -> dict:
    """φ-矩判据 D（离散层）：基 Φ={r, r²}，一个周期等权采样。"""
    period = z <= z[-1] / N_PERIODS + 1e-9
    r_all = np.array([s["r"][period] for s in strands])        # (N, S)
    mean_r = r_all.mean(axis=1)
    out, basis = {}, ["r", "r^2"]
    for key, phi in zip(basis, [r_all, r_all ** 2]):
        sums = phi.sum(axis=1)                                  # Σ_s Φ(pos_k(s))
        scale = np.abs(sums).max()
        d = np.abs(sums[:, None] - sums[None, :]).max() / scale
        out[key] = round(float(d), 6)
    return {
        "metrics": {
            "basis": basis,
            "mean_r": [round(float(m), 5) for m in mean_r],
            "D": out,
            "D_max": round(float(max(out.values())), 6),
        }
    }


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    for demo in (demo1_helix, demo2_transpose):
        data = demo()
        path = OUT_DIR / f"{data['name']}.json"
        path.write_text(json.dumps(data, ensure_ascii=False), encoding="utf-8")
        print(
            f"wrote {path.relative_to(OUT_DIR.parents[1])}  "
            f"N={data['N']}  samples={data['n_samples']}  "
            f"D_max={data['metrics']['D_max']}"
        )


if __name__ == "__main__":
    main()
