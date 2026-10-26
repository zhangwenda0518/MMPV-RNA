# -*- coding: utf-8 -*-
"""核查：权威源 445 个 run 里，那42 个未进镜像的 run 在 tab 里的物种标签是什么。

只读诊断，不写任何数据。
"""
import os
import re
import sys
from collections import Counter

# 本文件支持的环境变量：MMPV_SERVER_ROOT、MMPV_PAPER_ROOT（未设时保持原服务器字面量，路径逐字不变）
_SERVER = os.environ.get("MMPV_SERVER_ROOT", "/home/zhangwenda")
_PAPER  = os.environ.get("MMPV_PAPER_ROOT",  f"{_SERVER}/MMPV-paper")
_GOJI   = f"{_PAPER}/goji-virome"


sys.path.insert(0, f"{_GOJI}/00_scripts")

D = f"{_GOJI}/02_novel_virus/RNA-Lycium_barbarum_out/01_Assembly"
B = f"{_GOJI}/03_known_virus/known_virus_barbarum"
DATA = f"{_GOJI}/01_data"

RUN = re.compile(r"^(CRR|SRR|ERR|DRR)[0-9]+")
RS = re.compile(r"(?:CRR|SRR|ERR|DRR)[0-9]+")

# 1) 权威源 run -> 目录名，并统计 contig 大小
asm = {}
for d in sorted(os.listdir(D)):
    p = os.path.join(D, d)
    if not os.path.isdir(p):
        continue
    m = RUN.match(d)
    if not m:
        continue
    contig = [f for f in os.listdir(p) if "contig.fasta" in f]
    sz = os.path.getsize(os.path.join(p, contig[0])) if contig else 0
    asm[m.group(0)] = (d, sz)

# 2) 镜像里的 run
mirror = set()
for root, dirs, files in os.walk(B):
    for n in dirs + files:
        for m in RS.finditer(n):
            mirror.add(m.group(0))

missing = sorted(set(asm) - mirror)

# 3) 从 tab 找这些 run 的物种标签
TAB_FILES = ["gsa.tab", "ncbi.DNA.tab", "ncbi.tab", "ncbi.RNA.tab"]
run2sp = {}
for tf in TAB_FILES:
    p = os.path.join(DATA, tf)
    if not os.path.isfile(p):
        continue
    with open(p, encoding="utf-8", errors="replace") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        # 找 run 列与物种列
        ridx = None
        sidx = None
        for i, h in enumerate(header):
            hl = h.strip().lower()
            if ridx is None and hl in ("run", "run_acc", "sra_run", "run_accession"):
                ridx = i
            if sidx is None and hl in ("species", "organism", "taxid_name", "species_name"):
                sidx = i
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if ridx is None or ridx >= len(parts):
                continue
            r = parts[ridx].strip()
            if not RS.fullmatch(r):
                m = RS.search(r)
                if not m:
                    continue
                r = m.group(0)
            sp = parts[sidx].strip() if sidx is not None and sidx < len(parts) else "(无物种列)"
            run2sp.setdefault(r, set()).add(sp)

print("=" * 78)
print("权威源 %d run / 镜像 %d run / 缺失 %d" % (len(asm), len(mirror), len(missing)))
print("=" * 78)

# 按 contig 大小分层
big = [r for r in missing if asm[r][1] >= 100_000]
small = [r for r in missing if 0 < asm[r][1] < 100_000]
empty = [r for r in missing if asm[r][1] == 0]
print("\n【A】有实质组装 (>=100KB): %d 个" % len(big))
for r in big:
    print("   %-14s %-22s %9d B   物种标签=%s" % (
        r, asm[r][0], asm[r][1], sorted(run2sp.get(r, {"(tab 中无此 run)"}))))

print("\n【B】有小组装 (0 < x < 100KB): %d 个" % len(small))
for r in small:
    print("   %-14s %-22s %9d B   物种标签=%s" % (
        r, asm[r][0], asm[r][1], sorted(run2sp.get(r, {"(tab 中无此 run)"}))))

print("\n【C】空组装 (0 B): %d 个 —— 本就无组装结果，不构成遗漏" % len(empty))
print("   ", ", ".join(empty))

print("\n【D】缺失 run 的物种标签分布（全部 42 个）:")
cnt = Counter()
for r in missing:
    for sp in run2sp.get(r, {"(tab 中无此 run)"}):
        cnt[sp] += 1
for sp, n in cnt.most_common():
    print("   %-45s %d" % (sp[:45], n))