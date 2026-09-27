#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""修 generate_pipeline_report.py 的 stage 2 元数据收集模式。

根因（实机实测）: Path.rglob("**/metadata_association/**") 在 Python 3.10 下
只返回目录、不返回文件（实测 total=13 / files=0 / dirs=13），
于是 _stage_files 的 `if f.is_file()` 把 135 张元数据图 + 汇总表全部过滤掉，
[2b] meta 区块在 Stage 10 报告里从不出现。
修复: 模式改为 "**/metadata_association/**/*"（实测 files=272，含 135 png + 1 tsv）。
"""
import hashlib
import io
import os
import shutil
import sys

P = "/home/zhangwenda/MMPV-RNA/virome_analysis_pipeline/generate_pipeline_report.py"
BAK = P + ".bak_metaembed_20260915"
OLD = '2:["**/metadata_association/**"],'
NEW = '2:["**/metadata_association/**/*"],'


def md5(path):
    return hashlib.md5(open(path, "rb").read()).hexdigest()


src = io.open(P, encoding="utf-8").read()
print("before: md5=%s size=%d" % (md5(P), os.path.getsize(P)))
print("crlf count: %d" % src.count("\r\n"))

if OLD not in src:
    if NEW in src:
        print("ALREADY PATCHED, nothing to do")
        sys.exit(0)
    print("ERROR: 目标字符串未找到，未改动")
    sys.exit(1)

if src.count(OLD) != 1:
    print("ERROR: 目标字符串出现 %d 次，未改动" % src.count(OLD))
    sys.exit(1)

if not os.path.exists(BAK):
    shutil.copy2(P, BAK)
    print("backup -> %s (md5=%s)" % (BAK, md5(BAK)))
else:
    print("backup 已存在，不覆盖: %s (md5=%s)" % (BAK, md5(BAK)))

out = src.replace(OLD, NEW)
io.open(P, "w", encoding="utf-8", newline="").write(out)
print("after : md5=%s size=%d" % (md5(P), os.path.getsize(P)))
assert NEW in io.open(P, encoding="utf-8").read()
print("PATCH OK")
