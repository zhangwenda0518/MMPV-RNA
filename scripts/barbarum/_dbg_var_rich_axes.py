#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""列出 `_VariationRich` 各轴的 figure 坐标，用来精确判断标题/xlabel/字母是否重叠。"""
import os
import sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import make_sci_composite_figs as M  # noqa: E402

DATA = os.path.abspath(os.path.join(HERE, "..", "..", "reports",
                                    "sci_composite", "_diag", "data"))
OUT = os.path.abspath(os.path.join(HERE, "..", "..", "reports",
                                   "sci_composite", "_try_dbg"))
os.makedirs(OUT, exist_ok=True)

_orig = plt.Figure.savefig


def spy(self, fname, *a, **kw):
    fig_h = self.get_figheight()
    print(f"--- fig {self.get_figwidth():.2f} x {fig_h:.2f} in ---")
    rows = []
    for ax in self.axes:
        if not ax.get_visible():
            continue
        p = ax.get_position()
        title = ax.get_title()
        xl = ax.get_xlabel()
        rows.append((p.y1, p.y0, p.x0, p.x1, title[:46], xl[:34]))
    for y1, y0, x0, x1, t, xl in sorted(rows, reverse=True):
        print(f"y1={y1:.4f}({y1 * fig_h:6.2f}in) y0={y0:.4f}({y0 * fig_h:6.2f}in) "
              f"x=[{x0:.3f},{x1:.3f}] T='{t}' XL='{xl}'")
    return _orig(self, fname, *a, **kw)


plt.Figure.savefig = spy
M.build_var_rich(DATA, OUT, "Fig_KnownVirome_Composite")
