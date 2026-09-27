#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""核对 Stage 10 报告 HTML 是否真的嵌入了元数据关联（[2b] meta 区块）的图与表。

用法: python3 check_report_embed.py <dataset_dir>
判定方式: 取源文件 base64 的中段子串（确定性），在 HTML 里做子串匹配。
"""
import base64
import glob
import os
import sys


def main():
    ds = sys.argv[1].rstrip("/")
    html = os.path.join(ds, "10_report", "Pipeline_Summary_Report.html")
    with open(html, "r", encoding="utf-8", errors="replace") as fh:
        body = fh.read()
    print("HTML: %s  (%.1f MB, %d 字符)" % (html, os.path.getsize(html) / 1e6, len(body)))

    md = os.path.join(ds, "02_filtering", "metadata_association")
    files = []
    for ext in ("*.png", "*.jpg", "*.svg"):
        files += glob.glob(os.path.join(md, "**", ext), recursive=True)
    tsvs = glob.glob(os.path.join(md, "**", "*.tsv"), recursive=True)
    files.sort()
    tsvs.sort()

    hit = miss = 0
    for f in files:
        raw = open(f, "rb").read()
        if len(raw) < 400:
            continue
        b = base64.b64encode(raw).decode()
        probe = b[len(b) // 3: len(b) // 3 + 80]
        if probe in body:
            hit += 1
        else:
            miss += 1
            print("  MISS 图: %s (%d B)" % (os.path.relpath(f, ds), len(raw)))
    print("图嵌入: 命中 %d / 未命中 %d（共 %d 张）" % (hit, miss, hit + miss))

    thit = tmiss = 0
    for f in tsvs:
        txt = open(f, "r", encoding="utf-8", errors="replace").read()
        lines = [ln for ln in txt.splitlines() if ln.strip()]
        if len(lines) < 2:
            continue
        probe = lines[1].split("\t")[0][:40]
        if probe and probe in body:
            thit += 1
        else:
            tmiss += 1
            print("  MISS 表: %s (首数据行首字段=%r)" % (os.path.relpath(f, ds), probe))
    print("表嵌入: 命中 %d / 未命中 %d（共 %d 个）" % (thit, tmiss, thit + tmiss))

    for key in ("metadata", "Metadata", "关联", "Global_Summary"):
        print("  HTML 关键词 %-14s 出现 %d 次" % (key, body.count(key)))


if __name__ == "__main__":
    main()
