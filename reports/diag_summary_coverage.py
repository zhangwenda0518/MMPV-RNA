# -*- coding: utf-8 -*-
"""最后一环：42 个缺失样本能否进入已知病毒结果表。

判据= summary 表里的样本列是否含这 42 个 run。
若不含 -> 补不了结果，只能补bam/stat/logs（属于中间产物）。
只读。
"""
import csv
import os
import re
import sys

# 本文件支持的环境变量：MMPV_SERVER_ROOT、MMPV_PAPER_ROOT（未设时保持原服务器字面量，路径逐字不变）
_SERVER = os.environ.get("MMPV_SERVER_ROOT", "/home/zhangwenda")
_PAPER  = os.environ.get("MMPV_PAPER_ROOT",  f"{_SERVER}/MMPV-paper")
_GOJI   = f"{_PAPER}/goji-virome"


V2 = f"{_GOJI}/03_known_virus/known_virus_all_v2"
BAR = f"{_GOJI}/03_known_virus/known_virus_barbarum"
D = f"{_GOJI}/02_novel_virus/RNA-Lycium_barbarum_out/01_Assembly"

RUN = re.compile(r"^((?:CRR|SRR|ERR|DRR)[0-9]+)")
RS = re.compile(r"(?:CRR|SRR|ERR|DRR)[0-9]+")

asm = set()
for d in os.listdir(D):
    p = os.path.join(D, d)
    if os.path.isdir(p) and RUN.match(d):
        asm.add(RUN.match(d).group(0))

mirror = set()
for root, dirs, files in os.walk(os.path.join(BAR, "01_detection")):
    for n in list(dirs) + list(files):
        for m in RS.finditer(n):
            mirror.add(m.group(0))
missing = sorted(asm - mirror)
miss_set = set(missing)
print("权威源 %d / 镜像(detection) %d / 缺失 %d\n" % (len(asm), len(mirror), len(missing)))


def scan_summary(root, tag):
    sd = os.path.join(root, "01_detection", "summary")
    if not os.path.isdir(sd):
        print("[%s] 无 summary 目录" % tag)
        return
    print("=" * 74)
    print("[%s] 01_detection/summary 内容" % tag)
    print("=" * 74)
    for f in sorted(os.listdir(sd)):
        p = os.path.join(sd, f)
        if not os.path.isfile(p):
            continue
        sz = os.path.getsize(p)
        nrun = set()
        nrows = 0
        cols = []
        try:
            with open(p, encoding="utf-8", errors="replace") as fh:
                delim = "\t" if f.endswith((".tsv", ".txt")) else ","
                rdr = csv.reader(fh, delimiter=delim)
                head = next(rdr, None)
                cols = head or []
                for row in rdr:
                    nrows += 1
                    for c in row:
                        for m in RS.finditer(c or ""):
                            nrun.add(m.group(0))
        except Exception as e:
            print("   %-46s <读取失败: %s>" % (f, e))
            continue
        hit = len(nrun & miss_set)
        print("   %-46s %8d B  行数=%-7d 含run=%-5d 缺失run命中=%d" % (
            f, sz, nrows, len(nrun), hit))
        # 找出样本列
        if nrows and cols:
            samp = [i for i, c in enumerate(cols)
                    if re.search(r"sample|run|id", c, re.I)]
            if samp:
                print("        样本相关列: %s" % [cols[i] for i in samp][:5])


scan_summary(V2, "v2")
print()
scan_summary(BAR, "barbarum")

# 02_filtering 的 high_conf 里有没有这 42 个
print()
print("=" * 74)
print("02_filtering/high_conf.summary.tsv 覆盖")
print("=" * 74)
for tag, root in (("v2", V2), ("barbarum", BAR)):
    p = os.path.join(root, "02_filtering", "high_conf.summary.tsv")
    if not os.path.isfile(p):
        print("   [%s] 不存在" % tag)
        continue
    runs = set()
    with open(p, encoding="utf-8", errors="replace") as fh:
        head = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            for m in RS.finditer(line):
                runs.add(m.group(0))
    print("   [%s] 行=%-7d 含run=%-5d 缺失run命中=%d  列=%s" % (
        tag, sum(1 for _ in open(p, encoding='utf-8', errors='replace')) - 1,
        len(runs), len(runs & miss_set), head[:6]))