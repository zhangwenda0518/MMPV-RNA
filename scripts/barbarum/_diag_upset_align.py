#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""_diag_upset_align.py — 量 UpSet 图里「交集柱」与「点阵列」的 x 位置是否对齐。

判据：两组 x 中心应逐一重合（TOL 2 px @ 原图分辨率）。
⚠ 背景条纹是 grey92(235) —— 阈值必须 < 215 才不会被当成墨迹。
"""
import sys
from PIL import Image

path = sys.argv[1]
INK = 215
im = Image.open(path).convert("RGB")
W, H = im.size
px = im.load()
print(f"图像 {W} x {H}")


def is_ink(x, y, thr=INK):
    r, g, b = px[x, y]
    return r < thr or g < thr or b < thr


# ── 找柱区基线：柱区里"暗像素"最多的那一行 ──
y_lo, y_hi = int(H * 0.14), int(H * 0.62)
best_y, best_c = None, -1
for y in range(y_lo, y_hi):
    c = sum(1 for x in range(W) if is_ink(x, y, 150))
    if c > best_c:
        best_c, best_y = c, y
print(f"柱区基线 y={best_y}（该行暗像素 {best_c}）")

# 基线上方 8 px 扫一横带：每个柱（哪怕很小）都有墨
band = range(max(0, best_y - 8), best_y - 1)


def runs(prof, thr=0):
    out, x = [], 0
    while x < len(prof):
        if prof[x] > thr:
            s = x
            while x < len(prof) and prof[x] > thr:
                x += 1
            out.append(((s + x - 1) / 2, s, x - 1))
        else:
            x += 1
    return out


bars = [0] * W
for x in range(W):
    bars[x] = sum(1 for y in band if is_ink(x, y))
rb = runs(bars)
print(f"\n柱 {len(rb)} 根（中心）:")
print("  " + " ".join(f"{c:.0f}" for c, _, _ in rb))

# ── 点阵：扫矩阵区，按列统计圆点墨量 ──
y_lo2, y_hi2 = int(H * 0.60), int(H * 0.97)
dots = [0] * W
for x in range(W):
    dots[x] = sum(1 for y in range(y_lo2, y_hi2) if is_ink(x, y))
rd = runs(dots, thr=max(4, 0.02 * (y_hi2 - y_lo2)))
print(f"\n点阵列 {len(rd)} 个（中心）:")
print("  " + " ".join(f"{c:.0f}" for c, _, _ in rd))

print("\n对齐残差（点阵中心 - 最近柱中心）:")
for c, _, _ in rd:
    if rb:
        near = min(rb, key=lambda r: abs(r[0] - c))
        flag = "" if abs(c - near[0]) <= 2 else "   <-- 未对齐"
        print(f"  dot@{c:7.1f}  ->  bar@{near[0]:7.1f}   d={c - near[0]:+6.1f}{flag}")
