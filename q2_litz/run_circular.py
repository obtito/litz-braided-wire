"""Q2 圆形排布复现实验：用我方 skfem 管线独立复现队友 COMSOL 的「圆形简单成束」结果。

对照数据（队友 COMSOL 6.2，2026-09-26，Q2_Circular 包）：
  圆形排布 Rac/Rdc = 4.40331，六角排布 = 4.65833，等面积实心(2.547mm) = 4.56990
  设计点：N=331，d=0.14mm，环心距 0.15mm，铜面积 5.0953 mm²，Jdc=3.925 A/mm²
物理口径与队友一致：每丝独立导体 + 共享电压（并联）+ 总流 20 A 约束。

附带网格边界收敛阶梯：丝边界段数 n≈πd/size_min，多边形面积亏损 ~1/n²，
size 逐级减半 → 面积亏损 2.3%→0.6%→0.15%，比值应收敛到队友 COMSOL 值。
产出：data/q2c_results.json, data/q2c_strand_currents.csv（取最细一档）
"""
import json
import math
import sys
import time

import numpy as np

sys.path.insert(0, ".")
from sim.constants import RHO, F0, I_RMS                     # noqa: E402
from sim.gmsh_mesh import make_bundle_mesh                   # noqa: E402
from q2_litz.solve_bundle import load_bundle, solve_bundle, strand_currents  # noqa: E402

N_RINGS = 10          # 环 k 含 6k 丝（k=0 为中心 1 丝）→ 总数 331
SPACING = 0.15e-3     # 环心距（= 丝径 0.14 + 铜间隙 0.01 mm）
R_S = 0.07e-3         # 丝半径 0.07 mm
REF = dict(comsol_circular=4.40330542877846, comsol_hex=4.65832849766564,
           comsol_solid=4.56990309812434, ideal_strand=1.0010481467460408)
SIZE_LADDER = [1.0, 0.5, 0.25]   # × 0.36·r_s 基准（= 现网 q2 默认 → 逐级减半）


def circular_centers():
    """同心圆环排布：环 k 半径 k·spacing，丝数 6k，环内均匀。与队友 layout_parameters.json 一致。"""
    pts = [(0.0, 0.0)]
    for k in range(1, N_RINGS + 1):
        for i in range(6 * k):
            th = 2 * math.pi * i / (6 * k)
            pts.append((k * SPACING * math.cos(th), k * SPACING * math.sin(th)))
    return np.array(pts)


def run_once(centers, size_scale, tag):
    a_strand = math.pi * R_S**2
    rdc = RHO / (len(centers) * a_strand)
    r_strand_dc = RHO / a_strand
    size_min = 0.36 * R_S * size_scale
    t0 = time.time()
    path = make_bundle_mesh(centers, R_S, size_min=size_min,
                            path=f"data/q2c_mesh_{tag}.msh")
    mesh, ind_cu, strand_ind = load_bundle(path, centers, R_S)
    res = solve_bundle(mesh, ind_cu)
    Ii, _ = strand_currents(res["basis"], strand_ind, res["a"], res["v"])
    dt = time.time() - t0
    rac = res["p"] / I_RMS**2
    ratio = rac / rdc
    a_analytic = len(centers) * a_strand
    deficit = res["area_cu"] / a_analytic - 1
    print(f"size×{size_scale:<4} ({size_min*1e6:.1f}μm) | dof={res['dof']:>6} [{dt:5.1f}s]"
          f" | 铜面积 {res['area_cu']*1e6:.4f} mm² (亏 {deficit*100:+.2f}%)"
          f" | Rac/Rdc = {ratio:.5f} (vs COMSOL {(ratio/REF['comsol_circular']-1)*100:+.2f}%)")
    return dict(size_scale=size_scale, ratio=ratio, deficit=deficit, res=res,
                Ii=Ii, rdc=rdc, r_strand_dc=r_strand_dc)


def main():
    centers = circular_centers()
    N = len(centers)
    assert N == 331, N
    print(f"N={N}  d={2*R_S*1e3}mm  环距 {SPACING*1e3}mm  "
          f"解析铜面积 {N*math.pi*R_S**2*1e6:.4f} mm²（队友 5.0953）\n")

    results = [run_once(centers, s, f"{int(s*100)}") for s in SIZE_LADDER]
    best = results[-1]                                        # 最细一档入库
    ratio, deficit, res = best["ratio"], best["deficit"], best["res"]
    Ii, rdc, r_strand_dc = best["Ii"], best["rdc"], best["r_strand_dc"]

    imag = np.abs(Ii)
    avg = I_RMS / N
    eta_mag = np.abs(imag - avg).max() / avg
    eta_cplx = np.abs(Ii - avg).max() / avg
    p_transport = float(np.sum(r_strand_dc * imag**2))
    p_total = float(res["p"])
    p_within = p_total - p_transport
    p_dc = I_RMS**2 * rdc

    ring = np.round(np.hypot(centers[:, 0], centers[:, 1]) / SPACING).astype(int)
    print("\n环 | 丝数 | 平均|I| [mA] | 环内极差 % | 损耗占比 %   （队友环均值: "
          "0.54/0.56/0.81/1.49/2.94/5.97/12.4/26.1/55.7/119.6/258.6）")
    for k in sorted(set(ring)):
        m = ring == k
        spread = (imag[m].max() - imag[m].min()) / imag[m].mean() * 100
        pl = float(np.sum(r_strand_dc * imag[m]**2)) / p_transport * 100
        print(f"{k:2d} | {m.sum():4d} | {imag[m].mean()*1e3:8.3f} | {spread:9.3f} | {pl:8.3f}")
    print(f"\nη_complex = {eta_cplx*100:.1f}%（队友 338.2）  "
          f"η_magnitude = {eta_mag*100:.1f}%（队友 328.0）")
    print(f"损耗 W/m：dc {p_dc:.4f} + 分布失衡 {p_transport-p_dc:.4f} + 丝内超额 {p_within:.4f}")

    np.savetxt("data/q2c_strand_currents.csv",
               np.column_stack([np.arange(N), ring, centers[:, 0], centers[:, 1],
                                imag, np.angle(Ii, deg=True)]),
               delimiter=",",
               header="strand,ring,x_m,y_m,I_amplitude_A,I_phase_deg", comments="")
    json.dump(dict(N=N, d_mm=2*R_S*1e3, layout="circular_concentric",
                   spacing_mm=SPACING*1e3, a_cu_mm2=res["area_cu"]*1e6,
                   area_deficit_percent=deficit*100,
                   size_ladder=[dict(size_scale=x["size_scale"], ratio=x["ratio"],
                                     area_deficit_percent=x["deficit"]*100)
                                for x in results],
                   rdc_mohm=rdc*1e3, rac_mohm=res["p"]/I_RMS**2*1e3, ratio=ratio,
                   eta_complex=eta_cplx, eta_magnitude=eta_mag,
                   p_dc_W=p_dc, p_transport_W=p_transport, p_within_W=p_within,
                   dof=int(res["dof"]),
                   ref_comsol_circular=REF["comsol_circular"],
                   delta_vs_comsol_percent=(ratio/REF["comsol_circular"]-1)*100),
              open("data/q2c_results.json", "w"), ensure_ascii=False, indent=1)
    print("\n-> data/q2c_results.json, data/q2c_strand_currents.csv")


if __name__ == "__main__":
    main()
