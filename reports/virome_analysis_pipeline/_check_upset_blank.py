# -*- coding: utf-8 -*-
"""量：柱条 x 区间 / 点阵真实列中心 / 两者间隙（**排除柱缘抗锯齿假点**）。

⚠⚠ 坑：柱条右缘的抗锯齿是**暖色**的，会被"棕点"掩膜命中 → 必须用「点在**列内纵向分布**
   是否像一列点」来筛，最稳的办法是要求该列段在**多个行带**都有墨迹（真点列 7 行都有点）。
   这里用更简单的判据：段宽 20..45 px **且** 该段在 y 方向跨越 >= 5 个点行。
"""
import sys
import numpy as np
from PIL import Image

DOT_ROWS_Y = (1150, 1800)


def dot_cols(a):
    R, G, B = a[..., 0], a[..., 1], a[..., 2]
    brown = (R > 150) & (G > 40) & (G < 150) & (B < 95) & ((R - B) > 80)
    m = np.zeros_like(brown); m[DOT_ROWS_Y[0]:DOT_ROWS_Y[1], 900:] = True
    bm = brown & m
    cb = np.where(bm.sum(axis=0) > 0)[0]
    if len(cb) == 0:
        return []
    segs, s = [], cb[0]
    for i in range(1, len(cb)):
        if cb[i] - cb[i - 1] > 8:
            segs.append((s, cb[i - 1])); s = cb[i]
    segs.append((s, cb[-1]))
    out = []
    for a0, a1 in segs:
        if not (20 <= a1 - a0 <= 45):
            continue
        col = bm[:, a0:a1 + 1].sum(axis=1)
        idx = np.where(col > 0)[0]
        # 真点列：纵向至少出现 5 段（7 行点）
        seg2, s2 = [], idx[0]
        for i in range(1, len(idx)):
            if idx[i] - idx[i - 1] > 10:
                seg2.append((s2, idx[i - 1])); s2 = idx[i]
        seg2.append((s2, idx[-1]))
        if len(seg2) >= 5:
            out.append((int(a0), int(a1)))
    return out


def measure(p, tag):
    a = np.asarray(Image.open(p).convert("RGB")).astype(np.int16)
    R, G, B = a[..., 0], a[..., 1], a[..., 2]
    ink = (a.min(axis=2) < 245)
    grey = (abs(R - G) < 10) & (abs(G - B) < 10) & (R > 55) & (R < 110)
    band = grey[1170:1610, :]
    cols = np.where(band.sum(axis=0) > 20)[0]
    cols = cols[cols > 250]
    bar_r = int(cols.max())
    dots = dot_cols(a)
    if not dots:
        print(f"  {tag}: 未检出点列"); return
    d0 = dots[0]
    gap = d0[0] - bar_r - 1
    z = ink[1170:1610, bar_r + 1:d0[0]]
    free = (z.sum(axis=0) == 0)
    best = cur = 0
    for v in free:
        cur = cur + 1 if v else 0
        best = max(best, cur)
    dcen = [(x + y) // 2 for x, y in dots]
    sp = [dcen[i + 1] - dcen[i] for i in range(len(dcen) - 1)]
    print(f"  {tag:<16} 柱条右缘 {bar_r:>5} | 点列 {len(dots):>2} 个 首列 {d0[0]:>5} "
          f"| 间隙 {gap:>4} px（无墨迹带 {best:>4} px）| 列间距 {min(sp) if sp else 0}..{max(sp) if sp else 0}")


if __name__ == "__main__":
    for arg in sys.argv[1:]:
        r, p = arg.split("=", 1)
        measure(p, r)
