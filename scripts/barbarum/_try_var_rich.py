#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""本机快速迭代 `_VariationRich` 版式：只跑 build_var_rich，不碰主图其余面板。

用法: python _try_var_rich.py [out_dir]
"""
import os
import sys
import matplotlib
matplotlib.use("Agg")

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import make_sci_composite_figs as M  # noqa: E402

DATA = os.path.abspath(os.path.join(
    HERE, "..", "..", "reports", "sci_composite", "_diag", "data"))
OUT = sys.argv[1] if len(sys.argv) > 1 else os.path.abspath(os.path.join(
    HERE, "..", "..", "reports", "sci_composite", "_try"))
os.makedirs(OUT, exist_ok=True)

if not os.path.isdir(DATA):
    raise SystemExit(f"[err] data dir not found: {DATA}")

print(f"[data] {DATA}")
print(f"[out ] {OUT}")
M.build_var_rich(DATA, OUT, "Fig_KnownVirome_Composite")
