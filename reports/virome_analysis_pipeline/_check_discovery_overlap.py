#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""discovery / known-virus 静态图版式自动碰撞检测（通用）。

判据与 `_check_main_overlap.py` **完全一致**（同一套 `_collect` / `_in_range` /
两两求交），只是把被检测对象换成 `--module` 指定的作图模块的 `main()`：

    劫持 `plt.Figure.savefig` → `fig.canvas.draw()` → 遍历每个 Text 图元
    （轴刻度 / xlabel / ylabel / title / ax.texts / legend）的**显示坐标 bbox**
    两两求交；剔除 `ax.axison == False` 的轴与超出 view limit 的刻度；
    判据 **0 overlapping pairs**，另报越出画布的字母（面板字母是外挂 text）。

用法（`--tol` / `--inrange-px` / `--module` 之外的参数原样转给 `main()`）：
    python _check_discovery_overlap.py --tol 0.5 --inrange-px 0 -- \
        --root <RNA-Lycium_barbarum_out> --outdir <DIR>
    python _check_discovery_overlap.py --module redraw_barbarum_detection_figs -- \
        -v ... -m ... -d ... --assembly-root ...
"""
import importlib
import os
import sys
import itertools

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

_argv = list(sys.argv)
TOL = 0.5
if "--tol" in _argv:
    i = _argv.index("--tol")
    TOL = float(_argv[i + 1])
    del _argv[i:i + 2]
INRANGE_PX = 0.0
if "--inrange-px" in _argv:
    i = _argv.index("--inrange-px")
    INRANGE_PX = float(_argv[i + 1])
    del _argv[i:i + 2]
MODULE = "make_discovery_sci_figs"
if "--module" in _argv:
    i = _argv.index("--module")
    MODULE = _argv[i + 1]
    del _argv[i:i + 2]
if "--" in _argv:
    i = _argv.index("--")
    main_argv = _argv[i + 1:]
else:
    main_argv = _argv[1:]
sys.argv = [sys.argv[0]] + main_argv
print(f"[module] {MODULE}   [tol ] {TOL}   [inrange-px] {INRANGE_PX}")
print(f"[argv] {' '.join(main_argv)}")

M = importlib.import_module(MODULE)  # noqa: E402


def _collect(fig):
    items = []
    for ai, ax in enumerate(fig.axes):
        if not ax.get_visible():
            continue
        cand = []
        ab = ax.get_window_extent()
        if ax.axison:
            if ax.xaxis.get_visible():
                cand += [("xtick", t) for t in ax.get_xticklabels()]
                cand += [("xtick", t) for t in ax.get_xticklabels(minor=True)]
            if ax.yaxis.get_visible():
                cand += [("ytick", t) for t in ax.get_yticklabels()]
                cand += [("ytick", t) for t in ax.get_yticklabels(minor=True)]
        # axis("off") 只隐藏轴线/刻度/spine，**不**隐藏 title 与轴标签
        cand += [("xlabel", ax.xaxis.label), ("ylabel", ax.yaxis.label)]
        cand += [("title", ax.title), ("title", ax._left_title),
                 ("title", ax._right_title)]
        cand += [("text", t) for t in ax.texts]
        for lg in ([ax.get_legend()] if ax.get_legend() else []):
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
    if kind == "xtick":
        cx = 0.5 * (b.x0 + b.x1)
        return ab.x0 - INRANGE_PX <= cx <= ab.x1 + INRANGE_PX
    if kind == "ytick":
        cy = 0.5 * (b.y0 + b.y1)
        return ab.y0 - INRANGE_PX <= cy <= ab.y1 + INRANGE_PX
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
            print(f"  [TEXT-OFF-CANVAS] {k} {t.get_text()[:12]!r} bbox="
                  f"({b.x0:.0f},{b.y0:.0f},{b.x1:.0f},{b.y1:.0f}) canvas=({fw},{fh})")
    print(f"  -> {bad} overlapping pairs")
    return bad


_orig = plt.Figure.savefig
_TOTAL = {"bad": 0, "figs": 0}


def spy(self, fname, *a, **kw):
    _TOTAL["figs"] += 1
    _TOTAL["bad"] += check(self, os.path.basename(str(fname)))
    return _orig(self, fname, *a, **kw)


plt.Figure.savefig = spy
M.main()
print(f"\n[summary] {_TOTAL['figs']} figures checked, "
      f"{_TOTAL['bad']} overlapping pairs total")
sys.exit(1 if _TOTAL["bad"] else 0)
