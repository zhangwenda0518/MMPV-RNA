#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""下游影响面扫描：哪些产物表带 Genus/Species，且会被 1,246 行置空波及。

输入：
  - 校准表 calibration_20260914/final_integrated_classification.calibrated.tsv
    (只取 calib_action == blank 的 contig)
  - 数据集目录树
输出：命中文件清单 + 受影响行数
"""
import csv
import os
import sys

BASE = "/home/zhangwenda/MMPV-paper/goji-virome/02_novel_virus/RNA-Lycium_barbarum_out"
CAL = os.path.join(BASE, "05_Taxonomy", "Votus.integrated", "calibration_20260914",
                   "final_integrated_classification.calibrated.tsv")

PRUNE = {"01_Clean", "02_Deplete", "02b_Filter", "03_Assembly", "04_CLUSTER",
         "logs", "__pycache__", "tmp", "temp", "intermediate", "bam", "fastq", "fasta"}
EXTS = (".tsv", ".csv", ".txt", ".md", ".html", ".json", ".xlsx")
MAXSZ = 40 * 1024 * 1024
IDCOLS = ["contig_id", "Nucleotide", "seq_name", "query", "contig", "qseqid",
          "seq_id", "sequence_id", "member", "id", "Representative", "ref_id",
          "centroid_id", "Contig_id", "contigID", "ID"]


def p(*a):
    print(*a, flush=True)


blank = set()
with open(CAL, encoding="utf-8", errors="replace") as fh:
    for d in csv.DictReader(fh, delimiter="\t"):
        if (d.get("calib_action") or "").strip().strip('"') == "blank":
            blank.add((d.get("contig_id") or "").strip().strip('"'))
p("置空 contig 数 = %d" % len(blank))

hits = []
scanned = 0
for root, dirs, files in os.walk(BASE):
    dirs[:] = [d for d in dirs if d not in PRUNE]
    for fn in files:
        if not fn.lower().endswith(EXTS):
            continue
        fp = os.path.join(root, fn)
        try:
            if os.path.getsize(fp) > MAXSZ:
                continue
        except OSError:
            continue
        scanned += 1
        try:
            with open(fp, encoding="utf-8", errors="replace") as fh:
                head = fh.readline().rstrip("\n")
                if "Genus" not in head:
                    continue
                sep = "\t" if head.count("\t") >= head.count(",") else ","
                cols = [c.strip().strip('"') for c in head.split(sep)]
                idc = next((c for c in IDCOLS if c in cols), None)
                if idc is None:
                    hits.append((fp, len(cols), "无识别 ID 列", -1))
                    continue
                gi = cols.index(idc)
                nrow = 0
                napp = 0
                for line in fh:
                    nrow += 1
                    f2 = line.split(sep)
                    if gi < len(f2):
                        k = f2[gi].strip().strip('"')
                        if k in blank:
                            napp += 1
                hits.append((fp, nrow, idc, napp))
        except Exception as e:
            p("   [err] %s %s" % (fp, e))

p("扫描文件 %d 个，带 Genus 列 %d 个\n" % (scanned, len(hits)))
hits.sort(key=lambda x: -x[3])
tot_aff = 0
for fp, nrow, idc, napp in hits:
    rel = fp.replace(BASE + os.sep, "")
    tot_aff += max(napp, 0)
    p("   %-6d / %-8d  [%s]  %s" % (max(napp, 0), nrow if nrow >= 0 else -1, idc, rel))
p("\n受影响行合计 = %d" % tot_aff)
