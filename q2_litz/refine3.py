"""TASK B：refine=3 第三收敛点 + Q2 全量数字固化。

产出：
  data/q2_convergence.csv   —— 束 refine 1/2/3 收敛表（raw + 面积修正口径）
  data/q2_q3_baseline.json  —— Q2/Q3 基线汇总（束收敛、实心基准 exact+FEM、
                                 单丝趋肤因子、逐 refine η_I、L1 交叉校验）

面积修正口径（消网格圆化误差）：ratio_areacorr = rac_raw·(A_exact/A_mesh)/rdc_exact
其中 A_exact = N·π·r_s² = 5 mm²（设计值），rdc_exact = ρ/(N·π·r_s²)。
若 refine=3 直接解内存失败 → 回退 refine=2.5（size_min = refine3 的 size_min×1.26）。
"""
import json
import math
import sys
import time

import numpy as np

sys.path.insert(0, ".")

from sim.constants import RHO, F0, I_RMS, skin_depth
from sim.gmsh_mesh import make_bundle_mesh, make_wire_mesh
from geo.strand_packing import hex_packing
from q2_litz.impedance_model import strand_skin_factor, untwisted
import q2_litz.solve_bundle as sb
import q1_solid.fem2d_skfem as q1

N = 127
A_CU_EXACT = 5.0e-6                    # m²（J=4 A/mm² @ 20 A）
r_s = math.sqrt(A_CU_EXACT / (N * math.pi))
RDC_EXACT = RHO / (N * math.pi * r_s**2)
DELTA = skin_depth(F0)
S0 = min(0.45 * DELTA, 0.36 * r_s)     # refine=1 的 size_min 基准

centers, R_bundle = hex_packing(N, r_s)


def run_bundle(refine, size_min=None):
    """单个 refine 级：建网格 → 解 → 逐丝电流。返回标量字典，大对象即用即删。"""
    t0 = time.time()
    path = f"data/q2_bundle_r{refine}.msh"
    try:
        make_bundle_mesh(centers, r_s, refine=refine, size_min=size_min, path=path)
    except Exception:
        import gmsh
        try:
            gmsh.finalize()            # 崩在 generate 时防止占用句柄
        except Exception:
            pass
        raise
    t_mesh = time.time() - t0

    t1 = time.time()
    mesh, ind_cu, strand_ind = sb.load_bundle(path, centers, r_s)
    res = sb.solve_bundle(mesh, ind_cu)
    t_solve = time.time() - t1

    t2 = time.time()
    Ii, Ai = sb.strand_currents(res["basis"], strand_ind, res["a"], res["v"])
    t_cur = time.time() - t2

    rac = res["p"] / I_RMS**2
    a_mesh = res["area_cu"]
    imag = np.abs(Ii)
    imean = I_RMS / N
    out = dict(refine=refine, dof=int(res["dof"]), a_mesh_mm2=float(a_mesh) * 1e6,
               rac_mohm_raw=float(rac) * 1e3, ratio_raw=float(rac / RDC_EXACT),
               ratio_areacorr=float(rac * (A_CU_EXACT / a_mesh) / RDC_EXACT),
               eta_I=float((imag.max() - imean) / imean),
               i_max=float(imag.max()), i_min=float(imag.min()),
               sum_I=float(Ii.sum().real),
               t_mesh_s=float(t_mesh), t_solve_s=float(t_solve),
               t_currents_s=float(t_cur))
    del res, mesh, strand_ind, Ii, Ai, imag
    print(f"[refine={refine}] dof={out['dof']}  A_mesh={out['a_mesh_mm2']:.4f} mm²  "
          f"ratio_raw={out['ratio_raw']:.4f}  ratio_areacorr={out['ratio_areacorr']:.4f}  "
          f"η_I={out['eta_I']*100:.1f}%  ΣI={out['sum_I']:.4f} A  "
          f"(mesh {t_mesh:.0f}s + solve {t_solve:.0f}s + currents {t_cur:.0f}s)")
    return out


def run_solid(refine):
    """等面积实心基准（与束同一套 refine 口径）。"""
    a_eq = math.sqrt(A_CU_EXACT / math.pi)
    t0 = time.time()
    path = f"data/q2_solid_eq_r{refine}.msh"
    make_wire_mesh(a=a_eq, refine=refine, path=path)
    mesh, ind_cu = q1.load_mesh(path, a=a_eq)
    res = q1.solve(mesh, ind_cu, a=a_eq)
    dt = time.time() - t0
    a_exact = math.pi * a_eq**2
    out = dict(refine=refine, dof=int(res["dof"]),
               a_mesh_mm2=float(res["area_cu"]) * 1e6,
               rac_mohm_raw=float(res["rac"]) * 1e3,
               ratio_raw=float(res["ratio"]),
               ratio_areacorr=float(res["rac"] * (a_exact / res["area_cu"])
                                    / (RHO / a_exact)),
               z_re_mohm=float(res["z"].real) * 1e3, t_total_s=float(dt))
    del res, mesh, ind_cu
    print(f"[solid refine={refine}] dof={out['dof']}  ratio_raw={out['ratio_raw']:.4f}  "
          f"ratio_areacorr={out['ratio_areacorr']:.4f}  ({dt:.0f}s)")
    return out


if __name__ == "__main__":
    t_all = time.time()
    print(f"N={N}  d={2*r_s*1e3:.4f} mm  δ={DELTA*1e6:.2f} µm  "
          f"R_dc={RDC_EXACT*1e3:.4f} mΩ/m  A_exact={A_CU_EXACT*1e6:.4f} mm²")

    rows = [run_bundle(1), run_bundle(2)]
    try:
        rows.append(run_bundle(3))
    except (MemoryError, RuntimeError) as e:
        print(f"!! refine=3 失败（{type(e).__name__}: {e}）→ 回退 refine=2.5"
              f"（size_min = refine3·1.26 = {S0/3*1.26*1e6:.2f} µm）")
        rows.append(run_bundle(2.5, size_min=S0 / 3 * 1.26))

    solid = [run_solid(1), run_solid(2)]
    a_eq = math.sqrt(A_CU_EXACT / math.pi)
    solid_exact_ratio = q1.exact_ratio(a=a_eq)

    # L1 丝阻抗矩阵交叉校验（秒级，127×127）
    I_l1, V_l1, Z_l1 = untwisted(centers, r_s, 10 * R_bundle)
    i_l1 = np.abs(I_l1)
    l1 = dict(reZ_over_rdc=float(Z_l1.real / RDC_EXACT),
              eta_I=float((i_l1.max() - I_RMS / N) / (I_RMS / N)),
              note="灯丝模型不含邻近反应 → 已知偏低 ~10%")

    # ---- 收敛表 CSV ----
    with open("data/q2_convergence.csv", "w") as f:
        f.write("refine,dof,A_mesh_mm2,ratio_raw,ratio_areacorr\n")
        for r in rows:
            f.write(f"{r['refine']},{r['dof']},{r['a_mesh_mm2']:.6f},"
                    f"{r['ratio_raw']:.6f},{r['ratio_areacorr']:.6f}\n")

    # ---- 全量基线 JSON ----
    baseline = dict(
        generated=time.strftime("%Y-%m-%d %H:%M:%S"),
        design=dict(N=N, d_mm=2 * r_s * 1e3, a_cu_exact_mm2=A_CU_EXACT * 1e6,
                    bundle_od_mm=2 * R_bundle * 1e3,
                    rdc_exact_mohm=RDC_EXACT * 1e3,
                    f0_kHz=F0 / 1e3, i_rms=I_RMS),
        skin=dict(delta_um=DELTA * 1e6, d_over_delta=2 * r_s / DELTA,
                  strand_skin_factor_F=strand_skin_factor(r_s, F0)),
        bundle_convergence=rows,
        solid_equal_area=dict(D_mm=2 * a_eq * 1e3, exact_ratio=solid_exact_ratio,
                              exact_rac_mohm=(RHO / (math.pi * a_eq**2))
                              * solid_exact_ratio * 1e3,
                              fem=solid),
        L1_filament_untwisted=l1,
        notes=[
            "面积修正 ratio_areacorr = rac_raw·(A_exact/A_mesh)/rdc_exact，消除网格圆化亏面积",
            "单层绞合=刚性旋转：内部电流分配不变（勿用方位角平均建模单绞向）",
            "束 R_ac/R_dc 收敛值 vs 实心 4.53：未绞合 litz 比实心差 ~3% —— Q3 换位动机",
        ],
    )
    with open("data/q2_q3_baseline.json", "w") as f:
        json.dump(baseline, f, ensure_ascii=False, indent=1)

    print(f"\n收敛表（束，N=127）：")
    print(f"{'refine':>7} {'dof':>8} {'A_mesh':>9} {'raw':>9} {'areacorr':>9}")
    for r in rows:
        print(f"{r['refine']:>7} {r['dof']:>8} {r['a_mesh_mm2']:>8.4f} "
              f"{r['ratio_raw']:>9.4f} {r['ratio_areacorr']:>9.4f}")
    print(f"实心等面积：exact={solid_exact_ratio:.4f}  FEM(r2)={solid[1]['ratio_raw']:.4f}"
          f"/{solid[1]['ratio_areacorr']:.4f}(修正)")
    print(f"单丝趋肤因子 F={baseline['skin']['strand_skin_factor_F']:.5f}"
          f"  L1 Re(Z)/Rdc={l1['reZ_over_rdc']:.4f}")
    print(f"\n-> data/q2_convergence.csv, data/q2_q3_baseline.json"
          f"  （总耗时 {time.time()-t_all:.0f}s）")
