"""Q1 截断复核：空气域 30a → 60a（refine=2）。产出 data/q1_mesh_air30/60.msh 与 Δ%。"""
import sys
sys.path.insert(0, ".")
from sim.gmsh_mesh import make_wire_mesh
from q1_solid.fem2d_skfem import load_mesh, solve

out = {}
for raf in (30, 60):
    p = make_wire_mesh(a=1e-3, refine=2, r_air_factor=raf, path=f"data/q1_mesh_air{raf}.msh")
    m, ind = load_mesh(p, a=1e-3)
    out[raf] = solve(m, ind)["ratio"]
    print(f"r_air={raf:2d}a: Rac/Rdc = {out[raf]:.6f}")
print(f"Δ = {(out[60]/out[30]-1)*100:+.5f}%  (<1e-4% 即可忽略)")
