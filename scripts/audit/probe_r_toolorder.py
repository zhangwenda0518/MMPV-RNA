#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""把 pipeline 里所有 virus_classifier_analysis.R* 版本的 TOOL_BIAS / RANK_DEPTH_WEIGHTS 抓出来对比。"""
import glob, os, hashlib, time

PIPE = "/home/zhangwenda/MMPV-RNA/virome_discovery_pipeline"
KEYS = ("TOOL_BIAS", "RANK_DEPTH_WEIGHTS", "TOOL_PRIORITY", "PRIMARY", "fallback")

files = sorted(glob.glob(os.path.join(PIPE, "virus_classifier_analysis.R*")))
print("=== 匹配到的 R 脚本版本 %d 个 ===" % len(files))
for p in files:
    st = os.stat(p)
    h = hashlib.md5(open(p, "rb").read()).hexdigest()[:8]
    print("  %-62s %8d B  md5=%s  mtime=%s" % (
        os.path.basename(p), st.st_size, h,
        time.strftime("%Y-%m-%d %H:%M", time.localtime(st.st_mtime))))

for p in files:
    print("\n---- %s ----" % os.path.basename(p))
    with open(p, "r", encoding="utf-8", errors="replace") as f:
        for i, line in enumerate(f, 1):
            if any(k in line for k in KEYS):
                print("  %5d | %s" % (i, line.rstrip()[:200]))
