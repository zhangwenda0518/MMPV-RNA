#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""改进 Stage 10 报告 stage-level 文件的图/表标题: 用相对 stage 目录的路径替代裸文件名。

背景: 元数据区块 135 张图里有 8 个 Profile_* 子目录，裸文件名（如 vs_BioProject.png）
会重复 8 次，无法区分。改为 nested 文件显示 "04_Virus_Specific_Profiles/Profile_X/vs_BioProject.png"。
"""
import hashlib
import io
import os
import shutil
import sys

P = "/home/zhangwenda/MMPV-RNA/virome_analysis_pipeline/generate_pipeline_report.py"
BAK = P + ".bak_caption_20260915"
OLD = """            fn = f.name
            if fn.endswith('.tsv') or fn.endswith('.csv'):"""
NEW = """            try:
                fn = str(f.relative_to(spath)) if f.parent != spath else f.name
            except ValueError:
                fn = f.name
            if fn.endswith('.tsv') or fn.endswith('.csv'):"""


def md5(path):
    return hashlib.md5(open(path, "rb").read()).hexdigest()


src = io.open(P, encoding="utf-8").read()
print("before: md5=%s size=%d" % (md5(P), os.path.getsize(P)))
if NEW in src:
    print("ALREADY PATCHED")
    sys.exit(0)
if src.count(OLD) != 1:
    print("ERROR: 目标片段出现 %d 次，未改动" % src.count(OLD))
    sys.exit(1)
if not os.path.exists(BAK):
    shutil.copy2(P, BAK)
    print("backup ->", BAK)
io.open(P, "w", encoding="utf-8", newline="").write(src.replace(OLD, NEW))
print("after : md5=%s size=%d" % (md5(P), os.path.getsize(P)))
print("CAPTION PATCH OK")
