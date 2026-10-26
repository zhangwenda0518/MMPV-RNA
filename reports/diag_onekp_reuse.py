# -*- coding: utf-8 -*-
"""核实 42 个缺失样本在 onekp-virome 项目里是否已有可复用的 01_detection 结果。

发现: /home/zhangwenda/MMPV-paper/onekp-virome/onekp-analysis/01_detection/bam/
      下存在 ERR3487367.sorted.bam （bam 本体，v2 里只剩 .bai 死索引）

若这批样本在 onekp-analysis 里已有完整的 01_detection 产物，
那么"补齐 01_detection"可能不需要重跑 —— 可以直接复用其summary 表。
只读核实。
"""
import os
import re

# 本文件支持的环境变量：MMPV_SERVER_ROOT、MMPV_PAPER_ROOT（未设时保持原服务器字面量，路径逐字不变）
_SERVER = os.environ.get("MMPV_SERVER_ROOT", "/home/zhangwenda")
_PAPER  = os.environ.get("MMPV_PAPER_ROOT",  f"{_SERVER}/MMPV-paper")
_GOJI   = f"{_PAPER}/goji-virome"
_ONEKP  = f"{_PAPER}/onekp-virome"


D = f"{_GOJI}/02_novel_virus/RNA-Lycium_barbarum_out/01_Assembly"
B = f"{_GOJI}/03_known_virus/known_virus_barbarum"
OK = f"{_ONEKP}/onekp-analysis"

RUN = re.compile(r"^((?:CRR|SRR|ERR|DRR)[0-9]+)")
RS = re.compile(r"(?:CRR|SRR|ERR|DRR)[0-9]+")

asm = {}
for d in os.listdir(D):
    p = os.path.join(D, d)
    if os.path.isdir(p) and RUN.match(d):
        c = [f for f in os.listdir(p) if "contig.fasta" in f]
        asm[RUN.match(d).group(0)] = (
            d, os.path.getsize(os.path.join(p, c[0])) if c else 0)

mirror = set()
for root, dirs, files in os.walk(os.path.join(B, "01_detection")):
    for n in list(dirs) + list(files):
        for m in RS.finditer(n):
            mirror.add(m.group(0))
missing = sorted(set(asm) - mirror)

print("onekp-analysis 是否存在:", os.path.isdir(OK))
if os.path.isdir(OK):
    print("顶层:", ", ".join(sorted(os.listdir(OK))[:12]))
    print()

    # 1) bam 本体覆盖
    bamdir = os.path.join(OK, "01_detection", "bam")
    have_bam = {}
    if os.path.isdir(bamdir):
        for f in os.listdir(bamdir):
            if f.endswith(".sorted.bam"):
                m = RS.search(f)
                if m:
                    have_bam[m.group(0)] = os.path.getsize(os.path.join(bamdir, f))

    print("=" * 76)
    print("onekp 01_detection/bam 里 .bam 本体数 = %d" % len(have_bam))
    print("  42 个缺失目标中, onekp 有 bam 的 = %d" % len(set(missing) & set(have_bam)))
    print("=" * 76)
    for r in missing:
        if r in have_bam:
            print("   %-13s bam=%8.1f MB   (权威源 contig=%d B)" % (
                r, have_bam[r] / 1e6, asm[r][1]))

    # 2) summary 表是否覆盖这 42 个
    sd = os.path.join(OK, "01_detection", "summary")
    print()
    print("=" * 76)
    print("onekp 01_detection/summary 内容")
    print("=" * 76)
    if os.path.isdir(sd):
        for f in sorted(os.listdir(sd)):
            p = os.path.join(sd, f)
            if not os.path.isfile(p):
                continue
            if os.path.splitext(f)[1] not in (".tsv", ".txt", ".csv"):
                print("   %-44s (非表文件, %.1f KB)" % (f, os.path.getsize(p) / 1e3))
                continue
            runs = set()
            nrows = 0
            delim = "\t" if f.endswith((".tsv", ".txt")) else ","
            with open(p, encoding="utf-8", errors="replace") as fh:
                fh.readline()
                for line in fh:
                    if not line.strip():
                        continue
                    nrows += 1
                    for m in RS.finditer(line):
                        runs.add(m.group(0))
            print("   %-44s 行=%-7d 含run=%-5d 缺失run命中=%d" % (
                f, nrows, len(runs), len(runs & set(missing))))
    else:
        print("   (无 summary 目录)")

    # 3) 参考数据库是否可得（重跑定量需要）
    print()
    print("=" * 76)
    print("参考基因组(重跑定量必需)")
    print("=" * 76)
    for root, dirs, files in os.walk(OK):
        for f in files:
            if re.search(r"reference|ref_genome|virus_ref", f, re.I) and \
               f.endswith((".fa", ".fasta", ".fna")):
                print("   %s  (%.1f MB)" % (
                    os.path.join(root, f).replace(OK + "/", ""),
                    os.path.getsize(os.path.join(root, f)) / 1e6))