"""三环轮转+环循环拓扑（N=504=3环×24载纱器×7丝，d=0.1124mm）——
「d² 地板 × 换位」乘积效应的验证载体。

节拍（6 站一周期）：每站全部环轮转 c=5 槽；每 2 站执行环循环（载纱器 1→2→3→1 环，
环内槽号保持）。载纱器 = 换位单元（7 丝六角微束刚性随行——载纱器内不换位，
与 126 方案 R 同构，诚实口径）。
"""
import math
import sys

import numpy as np

sys.path.insert(0, ".")
from geo.trajectories import strand_offsets_in_bundle


class ThreeRingRoto:
    def __init__(self, a_cu=5.0e-6, n_rings=3, per_ring=24, strands_per=7,
                 c=5, flip_every=2, margin=1.02, gap=1.05, braid_pitch=40e-3):
        self.N = n_rings * per_ring * strands_per
        self.r_s = math.sqrt(a_cu / (self.N * math.pi))
        self.per_ring, self.n_rings, self.spp = per_ring, n_rings, strands_per
        self.c, self.flip_every = c, flip_every
        off = strand_offsets_in_bundle(self.r_s, margin)
        self.rb_eff = np.hypot(*off.T).max() + self.r_s
        # 环半径：环内相邻距 2R sin(π/24) ≥ 2·rb_eff·gap；环间径向差 ≥ 2.03·rb_eff
        r_min = 2 * self.rb_eff * gap / (2 * math.sin(math.pi / per_ring))
        self.radii = [r_min * (1 + i * (2.03 * self.rb_eff) / r_min) for i in range(n_rings)]
        self.offsets = [0.0, math.pi / per_ring, 0.0][:n_rings]
        self.S = self.flip_every * n_rings            # 6 站一周期
        self.P = braid_pitch
        self.off = off
        # 载纱器状态演化 -> 站点坐标
        self.stations = self._build_stations()

    def _carrier_pos(self, ring, slot):
        R = self.radii[ring]
        a = 2 * math.pi * slot / self.per_ring + self.offsets[ring]
        return np.array([R * math.cos(a), R * math.sin(a)])

    def _build_stations(self):
        # 每载纱器状态 (ring, slot)；初始 ring i 的槽 j
        state = [(i, j) for i in range(self.n_rings) for j in range(self.per_ring)]
        stations = []
        for k in range(self.S):
            new_state = []
            for ring, slot in state:
                s2 = (slot + self.c * k) % self.per_ring          # 每站轮转（累计）
                r2 = ring
                if k > 0 and k % self.flip_every == 0:            # 每 2 站环循环
                    r2 = (ring + 1) % self.n_rings
                new_state.append((r2, s2))
            state = new_state
            xy = np.array([self._carrier_pos(r, s) for r, s in state])   # (M,2)
            full = (xy[:, None, :] + self.off[None, :, :]).reshape(-1, 2)
            stations.append(full)
        return np.array(stations)                                  # (S, N, 2)

    # ---- 轨迹接口（与 RotoflipTrajectory 同构）----
    def _seg(self, z):
        u = (z / self.P * self.S) % self.S
        k = int(u)
        t = u - k
        return self.stations[k], self.stations[(k + 1) % self.S], t

    def positions(self, z):
        if np.isscalar(z):
            a, b, t = self._seg(z)
            return a + t * (b - a)
        return self.stations[int(z) % self.S]

    __call__ = positions

    def polar(self, z):
        p = self.positions(z)
        return np.hypot(p[:, 0], p[:, 1]), np.arctan2(p[:, 1], p[:, 0])

    def period_loss(self):
        return self.P                    # 物理节距（h∥/m̄ 量纲正确）

    def period_geometric(self):
        return float(self.S)

    def bundle_radius(self):
        return float(np.hypot(*self.stations[..., :].reshape(-1, 2).T).max() + self.r_s)

    def strand_length_factor(self, strand_idx=None):
        dz = self.P / self.S
        dxy = np.roll(self.stations, -1, axis=0) - self.stations
        seg = np.sqrt(1.0 + (np.hypot(dxy[..., 0], dxy[..., 1]) / dz) ** 2)
        per = seg.mean(axis=0)
        return float(per.mean()) if strand_idx is None else float(per[strand_idx])

    def check_overlap(self):
        dmin = np.inf
        for s in range(self.S):
            p = self.stations[s]
            dd = np.hypot(p[:, None, 0] - p[None, :, 0], p[:, None, 1] - p[None, :, 1])
            np.fill_diagonal(dd, np.inf)
            dmin = min(dmin, dd.min())
        return dmin


class RigidFrom:
    """把某站的截面冻结为刚性对照（同几何、零换位）。"""

    def __init__(self, base, station=0):
        self.centers = base.stations[station].copy()
        self.N = base.N
        self.r_s = base.r_s

    def positions(self, z):
        return self.centers

    __call__ = positions

    def polar(self, z):
        return np.hypot(*self.centers.T), np.arctan2(self.centers[:, 1], self.centers[:, 0])

    def period_loss(self):
        return 0.0

    def period_geometric(self):
        return 0.0

    def bundle_radius(self):
        return float(np.hypot(*self.centers.T).max() + self.r_s)

    def strand_length_factor(self, strand_idx=None):
        return 1.0


if __name__ == "__main__":
    import json
    import time
    from q3_braid.slice_engine import run_scheme
    from q4_perfect.criterion import criterion

    t0 = time.time()
    tri = ThreeRingRoto()
    dmin = tri.check_overlap()
    print(f"三环拓扑：N={tri.N} d={2*tri.r_s*1e3:.4f}mm 环半径={[round(r*1e3,2) for r in tri.radii]}mm "
          f"OD={2*tri.bundle_radius()*1e3:.2f}mm S={tri.S} 最小丝距={dmin*1e3:.4f}mm "
          f"(丝径{2*tri.r_s*1e3:.4f}, 安全={dmin>2*tri.r_s})")
    assert dmin > 2 * tri.r_s, "丝重叠"
    D = criterion(pos=tri.stations)
    D_ctl = criterion(pos=np.repeat(tri.stations[:1], tri.S, axis=0))
    print(f"判据 D：换位 {D['D']:.4f}（dwell {D['D_dwell']:.3f} / phi {D['D_phi']:.3f}） | 刚性对照 {D_ctl['D']:.4f}")

    res = {}
    res["A504_rigid"] = run_scheme(RigidFrom(tri), tri.r_s, S=1, tag="A504")
    res["R504_threering"] = run_scheme(tri, tri.r_s, S=tri.S, tag="R504")
    res["D_criterion"] = {k: D[k] for k in ("D", "D_dwell", "D_phi")}
    res["D_rigid"] = D_ctl["D"]
    json.dump(res, open("data/q3_d_sweep_504.json", "w"), ensure_ascii=False, indent=1, default=float)
    print(f"-> data/q3_d_sweep_504.json  [{time.time()-t0:.0f}s]")
