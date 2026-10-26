#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""`_VariationRich` 版式自动碰撞检测。

做法（沿用项目既有验证纪律）：渲染后把每个 Text 图元的**显示坐标 bbox** 取出来，
两两求交，报出所有重叠对与越界图元。比肉眼盯缩略图可靠。

用法: python _check_var_rich_overlap.py [--tol 0.5] [--data-dir A:B:...] [--out DIR]
"""
import os
import sys
import itertools
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import make_sci_composite_figs as M  # noqa: E402


def _arg(name, default=None):
    return sys.argv[sys.argv.index(name) + 1] if name in sys.argv else default


TOL = float(_arg("--tol", 0.5))
_dd = _arg("--data-dir")
DATA = ([d for d in _dd.split(os.pathsep) if d] if _dd else [
    os.path.abspath(os.path.join(HERE, "_archive", "fig_revisions",
                                 "sci_composite", "_diag", "data"))])
OUT = _arg("--out", os.path.abspath(os.path.join(
    HERE, "_archive", "fig_revisions", "sci_composite", "_try_dbg")))
os.makedirs(OUT, exist_ok=True)
print(f"[data] {DATA}")
print(f"[out ] {OUT}")

_fig = {}


def _collect(fig):
    items = []
    for ai, ax in enumerate(fig.axes):
        if not ax.get_visible():
            continue
        cand = []
        ab = ax.get_window_extent()
        if ax.axison:                      # axis("off") 的轴刻度不显示 -> 别当图元
            if ax.xaxis.get_visible():
                cand += [("xtick", t) for t in ax.get_xticklabels()]
                cand += [("xtick", t) for t in ax.get_xticklabels(minor=True)]
            if ax.yaxis.get_visible():
                cand += [("ytick", t) for t in ax.get_yticklabels()]
                cand += [("ytick", t) for t in ax.get_yticklabels(minor=True)]
            cand += [("xlabel", ax.xaxis.label), ("ylabel", ax.yaxis.label)]
            cand += [("title", ax.title), ("title", ax._left_title),
                     ("title", ax._right_title)]
        cand += [("text", t) for t in ax.texts]
        for lg in ([ax.get_legend()] if ax.get_legend() else []):
            cand += [("legend", t) for t in lg.get_texts()]
            if lg.get_title():
                cand += [("legend", lg.get_title())]
        for ch in ax.child_axes:
            for lg in ([ch.get_legend()] if ch.get_legend() else []):
                cand += [("legend", t) for t in lg.get_texts()]
                if lg.get_title():
                    cand += [("legend", lg.get_title())]
        for kind, t in cand:
            s = t.get_text().strip()
            if not s or not t.get_visible():
                continue
            items.append((f"ax{ai}({ax.get_xlabel()[:18]}|{ax.get_title()[:16]}):{kind}",
                          t, kind, ab))
    return items


def _in_range(kind, b, ab):
    """刻度标签超出 view limit 时 matplotlib 不绘制 -> 从检测里剔除。"""
    if kind == "xtick":
        cx = 0.5 * (b.x0 + b.x1)
        return ab.x0 - 24 <= cx <= ab.x1 + 24
    if kind == "ytick":
        cy = 0.5 * (b.y0 + b.y1)
        return ab.y0 - 24 <= cy <= ab.y1 + 24
    return True


def _bb(t, rend):
    try:
        return t.get_window_extent(renderer=rend)
    except Exception:
        return None


def check(fig, tag=""):
    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
    fw, fh = fig.canvas.get_width_height()
    items = []
    for k, t, kind, ab in _collect(fig):
        b = _bb(t, rend)
        if b is None or b.width < 0.1 or b.height < 0.1:
            continue
        if not _in_range(kind, b, ab):
            continue
        items.append((k, t, b))
    print(f"\n===== {tag}  canvas {fw}x{fh}px  texts={len(items)} =====")
    bad = 0
    for (k1, t1, b1), (k2, t2, b2) in itertools.combinations(items, 2):
        ix = min(b1.x1, b2.x1) - max(b1.x0, b2.x0)
        iy = min(b1.y1, b2.y1) - max(b1.y0, b2.y0)
        if ix > TOL and iy > TOL:
            bad += 1
            print(f"  [OVERLAP {ix:6.1f}x{iy:6.1f}px] {k1} {t1.get_text()[:34]!r}"
                  f"  <->  {k2} {t2.get_text()[:34]!r}")
    for k, t, b in items:
        if ":text" in k and (b.x0 < -TOL or b.x1 > fw + TOL):
            print(f"  [LETTER-OFF-CANVAS] {k} {t.get_text()[:12]!r} bbox="
                  f"({b.x0:.0f},{b.y0:.0f},{b.x1:.0f},{b.y1:.0f}) canvas=({fw},{fh})")
    print(f"  -> {bad} overlapping pairs")
    return bad


_orig = plt.Figure.savefig


def spy(self, fname, *a, **kw):
    _fig["f"] = self
    check(self, os.path.basename(str(fname)))
    return _orig(self, fname, *a, **kw)


plt.Figure.savefig = spy
M.build_var_rich(DATA, OUT, "Fig_KnownVirome_Composite")
