#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""把 /tmp/chk_chimera.py 的成品路径改指向 v66f count 产物，生成对照脚本。"""
import io
import os

SRC = "/tmp/chk_chimera.py"
DST = "/tmp/chk_chimera_v66f.py"
OLD = 'INTEG = os.path.join(BASE, "Votus.integrated", "final_integrated_classification.tsv")'
NEW = 'INTEG = "/tmp/v66f_count/final_integrated_classification.tsv"'

src = io.open(SRC, encoding="utf-8").read()
if OLD not in src:
    print("ERROR: 目标行未找到")
    raise SystemExit(1)
io.open(DST, "w", encoding="utf-8", newline="").write(src.replace(OLD, NEW))
print("written:", DST, os.path.getsize(DST), "B")
