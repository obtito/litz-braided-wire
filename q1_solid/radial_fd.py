"""Q1：实心圆铜导线趋肤效应——径向一维有限差分（FD）求解器。

物理模型（MQS 相量，e^{jωt}，全部有效值）：
    导体内 E_z(r) 满足  (1/r) d/dr ( r dE_z/dr ) = jωμσ E_z
    边界（安培环路，总电流 I 全部包围）：dE_z/dr |_{r=a} = jωμ I / (2πa)
    对称轴：dE_z/dr |_{r=0} = 0
    内阻抗（每米）：Z = E_z(a) / I ；R_ac = Re Z；L_int = Im Z / ω
    损耗（每米）：P = ∫ σ|E_z|² · 2πr dr ；交叉校验 R_ac = P / I²

FD：均匀径向网格 + 中心差分 -> 复数三对角方程组（numpy 直接求解）。
N → ∞ 时收敛于精确 Bessel 解 Z = kρ J0(ka)/(2πa J1(ka))，k = (1-j)/δ。
"""
import math
import sys
import numpy as np

sys.path.insert(0, ".")
from sim.constants import SIGMA, RHO, MU0, MU_R, F0, OMEGA, I_RMS, skin_depth


def solve_wire(a: float, f: float, n: int, i_rms: float = I_RMS):
    """返回 (Z[Ω/m], P[W/m], r[], |J(r)|[], J_phase[])。"""
    w = 2 * math.pi * f
    k2 = 1j * w * MU0 * MU_R * SIGMA
    h = a / (n - 1)
    r = np.linspace(0.0, a, n)

    # 三对角组装（有限体积 + 统一符号约定：行乘 -1 后保持与边界行同向）。
    # 内部节点 FV over [r-h/2, r+h/2]:
    #   rm·E_{i-1} - (rp+rm)·E_i + rp·E_{i+1} = k²·E_i·(r_i h²)
    # 轴心节点（镜像对称，精确到 O(h²)）：E_1 = E_0·(1 + k²h²/4)
    # 表面节点（Neumann，单侧差分）：(E_{n-1} - E_{n-2})/h = g
    main = np.zeros(n, dtype=complex)
    lower = np.zeros(n - 1, dtype=complex)
    upper = np.zeros(n - 1, dtype=complex)
    rhs = np.zeros(n, dtype=complex)

    main[0] = 4.0 / h**2 + k2
    upper[0] = -4.0 / h**2
    for i in range(1, n - 1):
        rp, rm = r[i] + h / 2, r[i] - h / 2
        c = 1.0 / (r[i] * h**2)
        lower[i - 1] += c * rm
        main[i] += -c * (rp + rm) - k2
        upper[i] += c * rp
    g = 1j * w * MU0 * MU_R * i_rms / (2 * math.pi * a)
    main[n - 1] = 1.0 / h
    lower[n - 2] = -1.0 / h
    rhs[n - 1] = g
    # 组装为稀疏求解
    from scipy.sparse import diags
    from scipy.sparse.linalg import spsolve
    A = diags([lower, main, upper], [-1, 0, 1], format="csc")
    E = spsolve(A, rhs)

    z = E[-1] / i_rms
    jr = SIGMA * np.abs(E)
    jp = np.angle(E, deg=True)
    # 损耗积分（梯形，含 2πr 面积元）
    p = np.trapezoid(SIGMA * np.abs(E) ** 2 * 2 * math.pi * r, r)
    return z, p, r, jr, jp


def exact_bessel(a: float, f: float):
    """精确 Bessel 解的内阻抗（每米），k=(1-j)/δ。"""
    from scipy.special import jv
    d = skin_depth(f)
    k = (1 - 1j) / d
    x = k * a
    return k * RHO * jv(0, x) / (2 * math.pi * a * jv(1, x))


if __name__ == "__main__":
    a = 1e-3  # Q1: D = 2 mm
    rdc = RHO / (math.pi * a * a)
    print(f"δ(200kHz) = {skin_depth(F0)*1e6:.2f} µm   a/δ = {a/skin_depth(F0):.3f}")
    print(f"R_dc = {rdc*1e3:.3f} mΩ/m\n")

    print("=== 网格收敛（N 径向节点数）@200 kHz ===")
    print(f"{'N':>6} {'R_ac/R_dc':>12} {'误差%':>8} {'P [W/m]':>9}")
    ze = exact_bessel(a, F0)
    target = ze.real / rdc
    print(f"{'exact':>6} {target:12.5f} {'—':>8} {I_RMS**2*ze.real:9.3f}")
    for n in (16, 32, 64, 128, 256, 512, 1024, 2048):
        z, p, *_ = solve_wire(a, F0, n)
        ratio = z.real / rdc
        err = (ratio / target - 1) * 100
        print(f"{n:6d} {ratio:12.5f} {err:8.2f} {p:9.3f}")

    print("\n=== 频率扫描（N=2048）vs 精确解 ===")
    print(f"{'f [kHz]':>8} {'δ[µm]':>8} {'FD R_ac/R_dc':>13} {'exact':>10} {'L_int[nH/m]':>12}")
    rows = []
    for f in (50e3, 100e3, 200e3, 500e3, 1e6):
        z, p, r, jr, jp = solve_wire(a, f, 2048)
        zx = exact_bessel(a, f)
        lint = z.imag / (2 * math.pi * f) * 1e9
        print(f"{f/1e3:8.0f} {skin_depth(f)*1e6:8.2f} {z.real/rdc:13.4f} {zx.real/rdc:10.4f} {lint:12.2f}")
        rows.append((f, skin_depth(f), z.real / rdc, zx.real / rdc, lint))
    np.savetxt("data/q1_sweep.csv", np.array(rows), delimiter=",",
               header="f_Hz,delta_m,Frac_over_Rdc,exact_over_Rdc,L_int_nH_per_m", comments="")
    z, p, r, jr, jp = solve_wire(a, F0, 2048)
    np.savetxt("data/q1_J_r.csv", np.column_stack([r, jr, jp]), delimiter=",",
               header="r_m,J_amplitude_A_per_m2,J_phase_deg", comments="")
    print("\n已写出 data/q1_sweep.csv 与 data/q1_J_r.csv")
    # J 的 e 折叠深度检查：J 降到表面值的 1/e 处
    js = jr[-1]
    idx = np.where(jr >= js / math.e)[0]
    d_eff = a - r[idx[0]] if len(idx) else float("nan")
    print(f"J(r) 从表面衰减到 1/e 的深度 ≈ {d_eff*1e6:.1f} µm（理论 δ = {skin_depth(F0)*1e6:.1f} µm）")
