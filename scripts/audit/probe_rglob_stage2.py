#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""诊断 generate_pipeline_report.py 的 _stage_files 对 02_filtering 的收集：
Path.rglob 尾随 ** 是否只返回目录、不返回文件。"""
import sys
from pathlib import Path

ds = Path(sys.argv[1] if len(sys.argv) > 1
          else "/home/zhangwenda/MMPV-paper/goji-virome/03_known_virus/known_virus_all_v2")
sp = ds / "02_filtering"
print("spath:", sp, sp.is_dir())

for pat in ("**/metadata_association/**", "**/metadata_association/*",
            "**/metadata_association/**/*", "**/metadata_association*"):
    try:
        r = list(sp.rglob(pat))
    except Exception as e:
        print("%-32s ERROR %s" % (pat, e))
        continue
    nf = sum(1 for f in r if f.is_file())
    nd = sum(1 for f in r if f.is_dir())
    print("%-32s total=%-6d files=%-6d dirs=%d" % (pat, len(r), nf, nd))

md = sp / "metadata_association"
print("metadata_association 直连 rglob('*'):", sum(1 for f in md.rglob("*") if f.is_file()))
print("其中 png:", sum(1 for f in md.rglob("*.png") if f.is_file()))
print("其中 tsv:", [str(f.name) for f in md.rglob("*.tsv") if f.is_file()])
