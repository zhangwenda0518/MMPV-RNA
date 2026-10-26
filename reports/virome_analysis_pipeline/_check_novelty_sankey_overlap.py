#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""fig_novelty_sankey.py 的版式体检（判据与 `_check_discovery_overlap.py` 一致）。

⚠ 为什么不能直接用 `_check_discovery_overlap.py`：
   该脚本会剔除 `ax.axison == False` 的轴，而桑基图**全部轴都是 `axis("off")`**
   → 直接跑会得到「0 对重叠」的**假通过**。这里改成遍历**所有** Text 图元。

判据：**0 对重叠**（TOL 0.5 px，同 `_check_discovery_overlap.py`）+ **无 TEXT-OFF-CANVAS**。
注意：文字压在彩带上不算碰撞（本图有意为之，靠白色描边保证可读），只判 text-text。

用法：
  "C:/Python312/python.exe" reports/virome_analysis_pipeline/_check_novelty_sankey_overlap.py
"""
import os
import sys
import tempfile
import itertools

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

TOL = 0.5
INRANGE_PX = 0.0

import fig_novelty_sankey as M   # noqa: E402

_captured = []
_orig = matplotlib.figure.Figure.savefig


def _spy(self, fname, *a, **kw):
    _captured.append(self)
    return None


matplotlib.figure.Figure.savefig = _spy
sys.argv = [sys.argv[0], "--outdir", tempfile.mkdtemp(prefix="sankey_qa_")]
try:
    M.main()
finally:
    matplotlib.figure.Figure.savefig = _orig

assert _captured, "没有捕获到 figure —— main() 未调用 savefig？"
fig = _captured[0]
# main() 里 plt.close(fig) 会把 canvas 摘掉 → 重新挂一个 Agg canvas 才能取 renderer
from matplotlib.backends.backend_agg import FigureCanvasAgg   # noqa: E402
FigureCanvasAgg(fig)
fig.canvas.draw()
ren = fig.canvas.get_renderer()
W, H = fig.canvas.get_width_height()

items = []
for ax in fig.axes:
    for t in ax.texts:
        items.append((ax, t))
for t in fig.texts:
    items.append((None, t))
for lg in fig.legends:
    for t in lg.get_texts():
        items.append((None, t))

def bb(t):
    try:
        return t.get_window_extent(renderer=ren)
    except Exception:
        return None

boxes = []
for ax, t in items:
    s = (t.get_text() or "").strip()
    if not s:
        continue
    b = bb(t)
    if b is None or b.width <= 0 or b.height <= 0:
        continue
    boxes.append((s, b))

# 1) 越出画布
off = [(s, b) for s, b in boxes
       if b.x0 < -INRANGE_PX or b.y0 < -INRANGE_PX or b.x1 > W + INRANGE_PX or b.y1 > H + INRANGE_PX]
print(f"画布 {W} x {H} px；Text 图元 {len(boxes)} 个")
if off:
    print(f"❌ TEXT-OFF-CANVAS {len(off)} 个：")
    for s, b in off[:12]:
        print(f"   {s[:44]!r:48s} bbox=({b.x0:.1f},{b.y0:.1f},{b.x1:.1f},{b.y1:.1f})")
else:
    print("✅ 无 TEXT-OFF-CANVAS")

# 2) 两两重叠
pairs = []
for (s1, b1), (s2, b2) in itertools.combinations(boxes, 2):
    ox = min(b1.x1, b2.x1) - max(b1.x0, b2.x0)
    oy = min(b1.y1, b2.y1) - max(b1.y0, b2.y0)
    if ox > TOL and oy > TOL:
        pairs.append((s1, s2, ox, oy))
print(f"{'✅' if not pairs else '❌'} 重叠对 {len(pairs)}（TOL {TOL} px）")
for s1, s2, ox, oy in pairs[:15]:
    print(f"   {s1[:36]!r}  ×  {s2[:36]!r}   overlap {ox:.1f} x {oy:.1f} px")

sys.exit(1 if (pairs or off) else 0)
