"""【已废弃——勿用】NGSolve 直连路线的死胡同记录（2026-09-25）。

两个环境级坑（macOS pip wheel ngsolve 6.2.2607）：
  1. netgen CSG2D 对同心圆相减的网格生成会无限卡死（>10 min 无输出）；
  2. netgen 的 gmsh .msh 导入器读出空网格（ne=0，即使 MSH 2.2 ASCII 格式正确）。
结论：NGSolve 仅保留作后续 3D 验证的备选，2D 生产路线见 fem2d_skfem.py
（gmsh 距离场网格 + meshio 桥 + scikit-fem 复对称装配，0.2–0.8 s/解）。
"""
