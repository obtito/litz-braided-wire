"""Q2/Q3 几何基元：N 等径丝在圆束内的堆积。

方案：六方晶格布点 + 圆域裁剪 + 间距二分，使丝数恰为 N。
确定性（无随机数）、局部密排 ≈ π/2√3、与真实紧压束一致。
输出每丝圆心坐标 (x_i, y_i)，保证互不重叠且在束圆内。
"""
import math

import numpy as np


def hex_packing(n_strands: int, strand_r: float, margin: float = 1.02):
    """返回 (centers[n,2], bundle_r)。margin：丝心间距 = 2r·margin（留漆膜几何间隙）。"""
    spacing = 2 * strand_r * margin

    def count_at(bundle_r):
        """半径 bundle_r 内的六方晶格点数（含原点）。"""
        if bundle_r <= 0:
            return 1, np.zeros((1, 2))
        rows = int(bundle_r / (spacing * math.sqrt(3) / 2)) + 2
        pts = [(0.0, 0.0)]
        for j in range(-rows, rows + 1):
            y = j * spacing * math.sqrt(3) / 2
            x_off = 0.0 if j % 2 == 0 else spacing / 2
            kmax = int((bundle_r - abs(x_off)) / spacing) + 1
            for k in range(-kmax, kmax + 1):
                x = x_off + k * spacing
                if math.hypot(x, y) <= bundle_r and (x, y) != (0.0, 0.0):
                    pts.append((x, y))
        return len(pts), np.array(pts)

    # 二分束半径，使丝数恰为 n_strands（hex 数 1,7,19,37,61,91,127,169 之外取最近可解 N）
    lo, hi = spacing, spacing * (math.sqrt(n_strands) + 2)
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        cnt, _ = count_at(mid)
        if cnt >= n_strands:
            hi = mid
        else:
            lo = mid
    # hi 是「恰好容纳 ≥n 的最小半径」；取该半径下的全部点并裁剪到 n
    cnt, pts = count_at(hi)
    if cnt < n_strands:
        raise RuntimeError(f"二分失败：count={cnt}")
    # 按半径排序取最内 n 个（若恰为完整六角数则全部保留）
    order = np.argsort(np.hypot(pts[:, 0], pts[:, 1]))
    chosen = pts[order[:n_strands]]
    used_r = float(np.hypot(chosen[:, 0], chosen[:, 1]).max())
    bundle_r = used_r + strand_r * margin  # 含丝半径与间隙
    return chosen, bundle_r


def copper_fill(n_strands, strand_r, bundle_r):
    return n_strands * math.pi * strand_r**2 / (math.pi * bundle_r**2)


if __name__ == "__main__":
    import sys
    sys.path.insert(0, ".")
    # 完整六角环数（1+6·3k(k+1)）：127 = 0..6 环；由目标铜截面反解丝径
    N = 127
    A_CU = 5.0e-6  # mm² -> m²（J=4 A/mm² @ 20A）
    r_s = math.sqrt(A_CU / (N * math.pi))
    centers, R = hex_packing(N, r_s)
    a_cu = N * math.pi * r_s**2
    print(f"N = {len(centers)}")
    print(f"束外径 D_bundle = {2*R*1e3:.3f} mm")
    print(f"铜截面 A_cu = {a_cu*1e6:.3f} mm²（目标 5 mm²，J = {20/a_cu/1e6:.2f} A/mm²）")
    print(f"束内铜填充率 = {copper_fill(N, r_s, R)*100:.1f}%（六方上限 90.7%）")
    # 最近丝心距（不得 < 2r）
    from itertools import combinations
    dmin = min(math.hypot(*(centers[i] - centers[j])) for i, j in combinations(range(len(centers)), 2))
    print(f"丝径 d = {2*r_s*1e3:.4f} mm；最近丝心距 = {dmin*1e3:.4f} mm（>d，安全）")
    np.savetxt("data/q2_strand_centers.csv", centers, delimiter=",",
               header="x_m,y_m", comments="")
    print("-> data/q2_strand_centers.csv")
