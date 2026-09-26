"""轮转+翻面拓扑的轨迹适配器：把 q4_topopt.json 的离散站点坐标变为
分段线性连续轨迹（positions/polar/period/长度因子），直接喂 slice_engine。"""
import json
import math
import sys

import numpy as np

sys.path.insert(0, ".")


class RotoflipTrajectory:
    """站点坐标（mm in JSON）→ 周期 P_braid 内的分段线性轨迹。"""

    def __init__(self, topopt_json="data/q4_topopt.json", which=0, braid_pitch=40e-3,
                 a_cu=5.0e-6):
        d = json.load(open(topopt_json))
        top = d["top3"][which] if "top3" in d else d["results"][which]
        carriers_mm = np.array(top["station_coords_carrier"])   # (S, M, 2) mm
        S, M, _ = carriers_mm.shape
        from geo.trajectories import strand_offsets_in_bundle
        off = strand_offsets_in_bundle(self._r_s_from(top))     # (7,2) m
        self.xy = (carriers_mm * 1e-3)[:, :, None, :] + off[None, None, :, :]  # (S, M*7, 2) m
        self.xy = self.xy.reshape(S, -1, 2)
        self.S, self.N, _ = self.xy.shape
        self.P = braid_pitch
        self.r_s = math.sqrt(a_cu / (self.N * math.pi))
        self.rho_max = np.hypot(self.xy[..., 0], self.xy[..., 1]).max()
        # 重叠安全自检
        dmin = np.inf
        for s in range(self.S):
            p = self.xy[s]
            dd = np.hypot(p[:, None, 0] - p[None, :, 0], p[:, None, 1] - p[None, :, 1])
            np.fill_diagonal(dd, np.inf)
            dmin = min(dmin, dd.min())
        assert dmin > 2 * self.r_s, f"丝重叠：{dmin*1e3:.4f} ≤ {2*self.r_s*1e3:.4f} mm"
        self._dmin = dmin

    @staticmethod
    def _r_s_from(top):
        return top["params"]["strand_d_mm"] * 1e-3 / 2

    def _seg(self, z):
        u = (z / self.P * self.S) % self.S
        k = int(u)
        t = u - k
        return self.xy[k], self.xy[(k + 1) % self.S], t

    def positions(self, z):
        a, b, t = self._seg(z)
        return a + t * (b - a)

    __call__ = positions

    def polar(self, z):
        p = self.positions(z)
        return np.hypot(p[:, 0], p[:, 1]), np.arctan2(p[:, 1], p[:, 0])

    def period_loss(self):
        return self.P

    def period_geometric(self):
        return self.P

    def bundle_radius(self):
        return float(self.rho_max + self.r_s)

    def strand_length_factor(self, strand_idx=None):
        """逐段 √(1+|Δxy/Δz|²) 的周期平均（分段线性：段内常数）。"""
        dxy = np.roll(self.xy, -1, axis=0) - self.xy          # (S,N,2)
        seg_fac = np.sqrt(1.0 + (np.hypot(dxy[..., 0], dxy[..., 1]) /
                                 (self.P / self.S)) ** 2)      # (S,N)
        per_strand = seg_fac.mean(axis=0)
        if strand_idx is None:
            return float(per_strand.mean())
        return float(per_strand[strand_idx])


if __name__ == "__main__":
    from q3_braid.slice_engine import run_scheme

    for which, name in [(0, "R_rotoflip_best"), (1, "R_rotoflip_2")]:
        traj = RotoflipTrajectory(which=which)
        print(f"{name}: S={traj.S} N={traj.N} m̄={traj.strand_length_factor():.4f} "
              f"OD≤{2*traj.bundle_radius()*1e3:.2f}mm 最小丝距={traj._dmin*1e3:.4f}mm")
        run_scheme(traj, traj.r_s, S=min(traj.S, 8), tag=name[:8])
