"""gmsh 网格生成（距离场加密）——Q1/Q2/Q3 共用。

注意两个 gmsh↔其它工具互操作坑（都已踩过）：
  1. occ.cut 默认 removeTool=True 会把工具实体删掉 → 必须 removeTool=False；
  2. netgen 的 gmsh 导入器在本 wheel 上损坏（读出空网格），一律走 meshio 桥。
"""
import math

import gmsh

from sim.constants import skin_depth


def _finish(path):
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.option.setNumber("Mesh.Binary", 0)
    gmsh.write(path)
    gmsh.finalize()
    return path


def make_wire_mesh(a=1e-3, f=200e3, refine=1, r_air_factor=30,
                   size_min_factor=0.45, path="data/q1_mesh.msh", quiet=True):
    """实心导线：铜盘 + 空气环，导线表面距离场加密（size_min_factor·δ）。"""
    delta = skin_depth(f)
    gmsh.initialize()
    if quiet:
        gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("wire")
    occ = gmsh.model.occ
    disk = occ.addDisk(0, 0, 0, a, a)
    air = occ.addDisk(0, 0, 0, r_air_factor * a, r_air_factor * a)
    occ.cut([(2, air)], [(2, disk)], removeObject=True, removeTool=False)
    occ.synchronize()

    surfaces = gmsh.model.getEntities(2)
    cu, ar = [], []
    for _, tag in surfaces:
        m = occ.getMass(2, tag)
        if m < math.pi * a * a * 1.1:
            cu.append(tag)
        else:
            ar.append(tag)
    gmsh.model.addPhysicalGroup(2, cu, name="copper")
    gmsh.model.addPhysicalGroup(2, ar, name="air")

    r_air = r_air_factor * a
    outer_curves = [e[1] for e in gmsh.model.getEntities(1)
                    if _on_circle(e[1], r_air)]
    inner_curves = [e[1] for e in gmsh.model.getEntities(1)
                    if not _on_circle(e[1], r_air)]
    gmsh.model.addPhysicalGroup(1, outer_curves, name="outer")

    f_d = gmsh.model.mesh.field.add("Distance")
    gmsh.model.mesh.field.setNumbers(f_d, "EdgesList", inner_curves)
    f_t = gmsh.model.mesh.field.add("Threshold")
    gmsh.model.mesh.field.setNumber(f_t, "InField", f_d)
    gmsh.model.mesh.field.setNumber(f_t, "SizeMin", size_min_factor * delta / refine)
    gmsh.model.mesh.field.setNumber(f_t, "SizeMax", 3e-3)
    gmsh.model.mesh.field.setNumber(f_t, "DistMin", 0.0)
    gmsh.model.mesh.field.setNumber(f_t, "DistMax", 15 * delta)
    gmsh.model.mesh.field.setAsBackgroundMesh(f_t)
    gmsh.option.setNumber("Mesh.Algorithm", 6)
    gmsh.model.mesh.generate(2)
    return _finish(path)


def _on_circle(tag, r, tol=1e-3):
    x0, y0, _, x1, y1, _ = gmsh.model.getBoundingBox(1, tag)
    return abs(max(abs(x0), abs(y0), abs(x1), abs(y1)) - r) < tol * r


def make_bundle_mesh(centers, strand_r, f=200e3, refine=1, r_air=None,
                     size_min=None, path="data/q2_bundle.msh", quiet=True):
    """N 丝束横截面网格：丝圆（互不重叠）+ 空气域；全部丝边界距离场加密。

    size_min：丝内/表面目标单元尺寸（默认 min(0.45δ, 0.36·r_s)/refine——
    既要分辨趋肤层又要保证每丝直径 ≥5 个单元）。
    """
    import numpy as _np
    delta = skin_depth(f)
    if size_min is None:
        size_min = min(0.45 * delta, 0.36 * strand_r) / refine
    bundle_r = _np.hypot(centers[:, 0], centers[:, 1]).max() + strand_r * 1.1
    if r_air is None:
        r_air = 10 * bundle_r
    gmsh.initialize()
    if quiet:
        gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("bundle")
    occ = gmsh.model.occ
    strand_tags = [occ.addDisk(float(x), float(y), 0, strand_r, strand_r)
                   for x, y in centers]
    air = occ.addDisk(0, 0, 0, r_air, r_air)
    occ.cut([(2, air)], [(2, t) for t in strand_tags],
            removeObject=True, removeTool=False)
    occ.synchronize()
    # 材料分组：丝=保留的工具盘；空气=切割余量
    strand_area = math.pi * strand_r**2
    cu, ar = [], []
    for _, tag in gmsh.model.getEntities(2):
        m = occ.getMass(2, tag)
        (cu if m < strand_area * 1.1 else ar).append(tag)
    gmsh.model.addPhysicalGroup(2, cu, name="copper")
    gmsh.model.addPhysicalGroup(2, ar, name="air")
    outer = [e[1] for e in gmsh.model.getEntities(1) if _on_circle(e[1], r_air)]
    gmsh.model.addPhysicalGroup(1, outer, name="outer")
    # 距离场：所有丝边界
    inner = [e[1] for e in gmsh.model.getEntities(1) if e[1] not in outer]
    f_d = gmsh.model.mesh.field.add("Distance")
    gmsh.model.mesh.field.setNumbers(f_d, "EdgesList", inner)
    f_t = gmsh.model.mesh.field.add("Threshold")
    gmsh.model.mesh.field.setNumber(f_t, "InField", f_d)
    gmsh.model.mesh.field.setNumber(f_t, "SizeMin", size_min)
    gmsh.model.mesh.field.setNumber(f_t, "SizeMax", 3e-3)
    gmsh.model.mesh.field.setNumber(f_t, "DistMin", 0.0)
    gmsh.model.mesh.field.setNumber(f_t, "DistMax", 20 * delta)
    gmsh.model.mesh.field.setAsBackgroundMesh(f_t)
    gmsh.option.setNumber("Mesh.Algorithm", 6)
    gmsh.model.mesh.generate(2)
    return _finish(path)
