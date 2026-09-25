"""Q2 截断复核 + Richardson 变体：r_air 10R→20R 与 raw 序列外推，写 data/q2_richardson_truncation.json。"""
import sys, math, json
sys.path.insert(0, ".")
import numpy as np
from scipy.optimize import curve_fit
from geo.strand_packing import hex_packing
from sim.gmsh_mesh import make_bundle_mesh
from q2_litz.solve_bundle import load_bundle, solve_bundle
from sim.constants import RHO, I_RMS

N, r_s = 127, math.sqrt(5.0e-6 / (127 * math.pi))
centers, Rb = hex_packing(N, r_s)
rdc = RHO / (N * math.pi * r_s ** 2)
out = {}
for raf in (10, 20):
    p = make_bundle_mesh(centers, r_s, refine=2, r_air=raf * Rb,
                         path=f"data/q2_bundle_air{raf}R.msh")
    m, ind, _ = load_bundle(p, centers, r_s)
    out[raf] = solve_bundle(m, ind)["p"] / I_RMS ** 2 / rdc
    print(f"r_air={raf:2d}R: ratio={out[raf]:.5f}")
trunc = (out[20] / out[10] - 1) * 100
y1, y2, y3 = 4.6833, 4.6578, 4.6527          # refine 1/2/3 raw（q2_convergence.csv）
h1, h2, h3 = 1.0, 1 / 2, 1 / 3
aitken = y3 - (y3 - y2) ** 2 / ((y3 - y2) - (y2 - y1))
popt, _ = curve_fit(lambda h, yinf, C, p: yinf + C * h ** p, [h1, h2, h3], [y1, y2, y3],
                    p0=[4.65, 0.03, 2.0])
yinf_r2 = (y3 * y1 - y2 ** 2) / (y1 - 2 * y2 + y3)
rich = {"aitken": round(aitken, 4), "lsq_p_fit": round(float(popt[0]), 4),
        "observed_p": round(float(popt[2]), 2), "assumed_r2_p2": round(yinf_r2, 4)}
print(f"Δtrunc = {trunc:+.4f}%  Richardson(raw) = {rich}")
json.dump({"truncation": {"air10R": out[10], "air20R": out[20], "change_pct": trunc, "refine": 2},
           "richardson_raw": rich}, open("data/q2_richardson_truncation.json", "w"), indent=1)
