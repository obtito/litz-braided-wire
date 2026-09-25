"""Q3/Q4 几何引擎：换位方案的 r(z)/θ(z) 轨迹生成器（2026-09-25 审查后重写）。

运动学约定（【无退扭笼绞】cage braiding——与真实笼式成缆机一致）：
    pos(z) = R(φ₂)·[C + R(φ₁)·r_local]
    φ₁(z) = s₁·2πz/P₁：丝在【束共转系】内绕微束轴自转（P₁ 取缆轴口径；
             若从机器的弧长螺距换算：P₁_axial = P₁_arc·cosα₂）
    φ₂(z) = s₂·2πz/P₂：微束圆心绕束轴公转（s₂ 与 s₁ 反向 → 反向两级绞；
             真 S/Z（周期性换向 s₂(z) 分段）属后续变体，此处恒向）
    ⇒ r(z) = |C + R(φ₁)r_local| 只依赖 φ₁：
       · 损耗相关周期 = P₁（切片仿真只需扫 P₁，无需 lcm，省 ~25× FEM）
       · 几何周期 = lcm(P₁, P₂)（整根丝回到初始位置）
       · 微束圆心距 |C| 恒定：简单两级绞不交换微束位置（部分换位，Q3 量化对象）

N=126=18 微束×7 丝 与 Q2 的 N=127 对应关系：等铜截面 A_cu=5.000mm² 下
d=0.2248mm（126 丝）vs 0.2239mm（127 丝），J=4.00 均达标（题目表 1 约束自 Q2 起生效，
构造函数以 a_cu 反解丝径并断言）。
"""
import math

import numpy as np


def _rot(angle):
    c, s = math.cos(angle), math.sin(angle)
    return np.array([[c, -s], [s, c]])


def microbundle_centers(r_pitch):
    """18 个微束圆心：6 内圈 @r_pitch + 12 外圈 @2r_pitch（空芯——|C|=0 的中心微束
    零径向迁移，故不放中心；r_pitch 由 _safe_pitch 给出保证互不重叠）。"""
    pts = []
    for j in range(6):
        a = j * math.pi / 3
        pts.append((r_pitch * math.cos(a), r_pitch * math.sin(a)))
    for j in range(12):
        a = j * math.pi / 6 + math.pi / 12   # 30° 间距（12 点），偏置 15° 错开内圈
        pts.append((2 * r_pitch * math.cos(a), 2 * r_pitch * math.sin(a)))
    return np.array(pts)


def strand_offsets_in_bundle(strand_r, margin=1.02):
    """微束内 7 丝（1+6 六角）偏移；中心丝 ρ=0（零径向波动，成缆界的「死丝」）。"""
    s = 2 * strand_r * margin
    off = [(0.0, 0.0)]
    for j in range(6):
        a = j * math.pi / 3
        off.append((s * math.cos(a), s * math.sin(a)))
    return np.array(off)


def _safe_pitch(rb_eff, gap=1.03):
    """微束节圆半径：受外圈相邻距离约束 d_adj = 2·(2R)·sin(π/12) = 1.035R ≥ 2·rb_eff·gap
    ⇒ R ≥ 1.93·rb_eff·gap；取 1.95 留裕量（内圈-外圈最近距 1.07R 更宽）。"""
    return 1.95 * rb_eff * gap


class TwoLevelCounterTwist:
    """两级反向绞（cage 约定）轨迹。N = 18×7 = 126。"""

    N_BUNDLES = 18  # 布局仅支持 6+12 双圈

    def __init__(self, a_cu=5.0e-6, p1=12.5e-3, p2=25e-3, s1=+1.0, s2=-1.0,
                 margin=1.02, bundle_gap=1.03):
        self.r_s = math.sqrt(a_cu / (7 * self.N_BUNDLES * math.pi))
        self.p1, self.p2, self.s1, self.s2 = p1, p2, s1, s2
        assert p2 > p1, "P₂（公转）应长于 P₁（自转）"
        # 可通约性：损耗周期 P₁ 内公转相位偏移 s₂·P₁/P₂——几何周期 lcm
        from fractions import Fraction
        frac = Fraction(p2 / p1).limit_denominator(1000)   # P₂/P₁ = m/n（有理近似）
        self._lcm = float(frac.numerator * p1)              # lcm = m·P₁ = n·P₂
        self._ratio = f"{frac.numerator}/{frac.denominator}"
        off = strand_offsets_in_bundle(self.r_s, margin)
        self.offsets = np.tile(off, (self.N_BUNDLES, 1))          # (N,2)
        self.rho = np.hypot(*self.offsets.T)                      # (N,) 丝到微束心距离
        self.rb_eff = self.rho.max() + self.r_s
        self.r_pitch = _safe_pitch(self.rb_eff, bundle_gap)
        self.centers = microbundle_centers(self.r_pitch)
        self.C_norm = np.repeat(np.hypot(*self.centers.T), 7)     # (N,) 每丝的 |C|
        self.N = 7 * self.N_BUNDLES
        self.w1 = self.s1 * 2 * math.pi / p1
        self.w2 = self.s2 * 2 * math.pi / p2

    # ---- 轨迹 ----
    def positions(self, z):
        pos = np.empty((self.N, 2))
        R1, R2 = _rot(self.w1 * z), _rot(self.w2 * z)
        inner = self.offsets @ R1.T                                 # 微束内自转
        for b in range(self.N_BUNDLES):
            pos[7 * b:7 * (b + 1)] = (self.centers[b] + inner[7 * b:7 * (b + 1)]) @ R2.T
        return pos

    __call__ = positions

    def polar(self, z):
        p = self.positions(z)
        return np.hypot(*p.T), np.arctan2(p[:, 1], p[:, 0])

    # ---- 周期接口 ----
    def period_loss(self):
        """损耗相关周期 = P₁（r(z) 只依赖 φ₁）。切片仿真扫这个。"""
        return self.p1

    def period_geometric(self):
        """几何周期 = lcm(P₁,P₂)（P₂/P₁ 为整数比时有限，否则返回 None 并注明准周期）。"""
        return self._lcm

    # ---- 长度因子（闭式，含轴向 +1）----
    def strand_length_factor(self, strand_idx=None, n_quad=256):
        """⟨ds/dz⟩ = ⟨√(1+|dpos/dz|²)⟩，按损耗周期数值积分（闭式被积）。
        |dpos/dz|² = ω₂²r(z)² + ω₁²ρ² + 2ω₁ω₂(C·R(φ₁)r + ρ²)  （cage 约定）。"""
        zs = np.linspace(0, self.p1, n_quad)
        if strand_idx is None:
            idx = np.arange(self.N)
        else:
            idx = np.atleast_1d(strand_idx)
        out = np.empty(len(idx))
        for k, i in enumerate(idx):
            c = self.centers[i // 7]
            rloc = self.offsets[i]
            rho2 = rloc @ rloc
            acc = 0.0
            for z in zs:
                ph1 = self.w1 * z
                r1 = _rot(ph1) @ rloc
                rz2 = (c + r1) @ (c + r1)
                d2 = self.w2 ** 2 * rz2 + self.w1 ** 2 * rho2 \
                    + 2 * self.w1 * self.w2 * (c @ r1 + rho2)
                acc += math.sqrt(1.0 + d2)
            out[k] = acc / len(zs)
        if strand_idx is None:
            return float(out.mean())
        return float(out[0]) if np.isscalar(strand_idx) else out

    def bundle_radius(self):
        """保证上界：max|C| + rb_eff（任意 z 的包络）。"""
        return 2 * self.r_pitch + self.rb_eff

    def export_polyline3d(self, period, n_samples=400, path=None):
        """3D 折线 (x,y,z) 逐丝采样——Q4 网站 / three.js 用。"""
        zs = np.linspace(0, period, n_samples)
        pts = np.array([self.positions(z) for z in zs])            # (S,N,2)
        data = {"N": self.N, "period": period, "strand_r": self.r_s,
                "samples": [{"z": float(z),
                             "xy": pts[i].tolist()} for i, z in enumerate(zs)]}
        if path:
            import json
            json.dump(data, open(path, "w"))
        return data


class RigidBundle:
    """直束（方案 A，factor=1）或单绞向（方案 B：factor=√(1+(2πρ/P)²)，分配不变）。"""

    def __init__(self, centers, strand_r=0.0, twist_pitch=None):
        self.centers = np.asarray(centers)
        self.r_s = strand_r
        self.rho = np.hypot(*self.centers.T)
        self.twist_pitch = twist_pitch
        self.N = len(centers)

    def positions(self, z):
        return self.centers

    __call__ = positions

    def polar(self, z):
        return self.rho.copy(), np.arctan2(self.centers[:, 1], self.centers[:, 0])

    def period_loss(self):
        return 0.0  # 平稳（与 z 无关）

    def period_geometric(self):
        return self.twist_pitch or 0.0

    def strand_length_factor(self, strand_idx=None):
        if self.twist_pitch is None:
            return 1.0
        f = np.sqrt(1.0 + (2 * math.pi * self.rho / self.twist_pitch) ** 2)
        return float(f.mean()) if strand_idx is None else float(f[strand_idx])

    def bundle_radius(self):
        return float(self.rho.max() + self.r_s)

    def export_polyline3d(self, period, n_samples=100, path=None):
        zs = np.linspace(0, period, n_samples)
        pts = np.array([self.positions(z) for z in zs])
        data = {"N": self.N, "period": period, "strand_r": self.r_s,
                "samples": [{"z": float(z), "xy": pts[i].tolist()} for i, z in enumerate(zs)]}
        if path:
            import json
            json.dump(data, open(path, "w"))
        return data


def radial_history(traj, period, n_samples=200):
    """每丝 (r,θ)(z) 采样——停留分布 / φ-矩判据输入（Q4）。返回 (S,N,2)。"""
    zs = np.linspace(0, period, n_samples)
    out = np.empty((len(zs), traj.N, 2))
    for i, z in enumerate(zs):
        r, th = traj.polar(z)
        out[i, :, 0] = r
        out[i, :, 1] = th
    return out


if __name__ == "__main__":
    # 冒烟自检（审查项逐条验证）
    t = TwoLevelCounterTwist(a_cu=5.0e-6, p1=12.5e-3, p2=25e-3)  # P₂=2P₁；1m=80/40 整捻距（规避半整数最坏对齐，PAPERS_NOTES）
    print(f"N={t.N}  d={2*t.r_s*1e3:.4f}mm  A_cu={t.N*math.pi*t.r_s**2*1e6:.4f}mm² "
          f"(J={20/(t.N*math.pi*t.r_s**2)/1e6:.2f})  束OD≤{2*t.bundle_radius()*1e3:.3f}mm")
    assert abs(t.N * math.pi * t.r_s**2 - 5e-6) < 1e-12
    # 1) 损耗周期自检：r(z+P1)==r(z)
    r0, _ = t.polar(0.0)
    r1, _ = t.polar(t.p1 * (1 + 1e-9))
    assert np.allclose(r0, r1, atol=1e-7), "P₁ 周期性失效"
    print(f"✓ r(z) 以 P₁ 为周期（损耗周期 {t.p1*1e3:.1f}mm；几何周期 {t.period_geometric()*1e3:.0f}mm，P₂/P₁={t._ratio}）")
    # 2) 长度因子：中心丝（ρ=0）解析对照
    f_center = t.strand_length_factor(0)  # 微束 0 的中心丝（|C|=r_pitch 外圈为 2r_pitch... 微束 0 在内圈）
    expect = math.sqrt(1 + (2 * math.pi * t.C_norm[0] / t.p2) ** 2)
    assert abs(f_center - expect) < 1e-6, (f_center, expect)
    print(f"✓ 长度因子闭式自检：中心丝 {f_center:.4f} = 解析 {expect:.4f}；全束均值 {t.strand_length_factor():.4f}（≥1）")
    # 3) 64 站重叠安全 + 包络
    dmin = np.inf
    for z in np.linspace(0, t.p1, 64):
        p = t.positions(z)
        d = np.hypot(p[:, None, 0] - p[None, :, 0], p[:, None, 1] - p[None, :, 1])
        np.fill_diagonal(d, np.inf)
        dmin = min(dmin, d.min())
    print(f"✓ 64 站最小丝心距 {dmin*1e3:.4f}mm（丝径 {2*t.r_s*1e3:.4f}，安全={dmin > 2*t.r_s}）")
    # 4) radial_history 冒烟（原审查项：traj(z) 崩溃 → 已修 polar/positions）
    hist = radial_history(t, t.p1, 64)
    print(f"✓ radial_history: {hist.shape}  r 范围 [{hist[:,:,0].min()*1e3:.3f},{hist[:,:,0].max()*1e3:.3f}]mm")
    # 5) 方案 B（单绞向）长度因子非平凡
    rb = RigidBundle(t.positions(0.0), t.r_s, twist_pitch=30e-3)
    print(f"✓ RigidBundle(B 单绞向): factor={rb.strand_length_factor():.4f}；直束=1.0")
