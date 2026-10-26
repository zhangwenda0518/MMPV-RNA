# -*- coding: utf-8 -*-
"""精确核实：42 个缺失样本在 v2 与 v3 各自的存在情况（逐阶段）。

用户问：v2 / v3 是不是也没有这批样本？
必须区分三种"没有"：
  (a) 目录里完全没有该run 的任何文件
  (b) 有中间产物(bam索引/log/stat) 但无结果表行
  (c) 有结果表行
只读。
"""
import os
import re

# 本文件支持的环境变量：MMPV_SERVER_ROOT、MMPV_PAPER_ROOT（未设时保持原服务器字面量，路径逐字不变）
_SERVER = os.environ.get("MMPV_SERVER_ROOT", "/home/zhangwenda")
_PAPER  = os.environ.get("MMPV_PAPER_ROOT",  f"{_SERVER}/MMPV-paper")
_GOJI   = f"{_PAPER}/goji-virome"


D = f"{_GOJI}/02_novel_virus/RNA-Lycium_barbarum_out/01_Assembly"
BASE = f"{_GOJI}/03_known_virus"
BAR = os.path.join(BASE, "known_virus_barbarum")
STAGES = ["01_detection", "02_filtering", "03_variants", "04_post_analysis",
          "05_assembly", "06_extraction", "07_capheine", "08_similarity",
          "09_dvg", "10_report"]

RUN = re.compile(r"^((?:CRR|SRR|ERR|DRR)[0-9]+)")
RS = re.compile(r"(?:CRR|SRR|ERR|DRR)[0-9]+")

asm = {}
for d in os.listdir(D):
    p = os.path.join(D, d)
    if os.path.isdir(p) and RUN.match(d):
        c = [f for f in os.listdir(p) if "contig.fasta" in f]
        asm[RUN.match(d).group(0)] = (
            d, os.path.getsize(os.path.join(p, c[0])) if c else 0)

bar_runs = set()
for root, dirs, files in os.walk(BAR):
    for n in list(dirs) + list(files):
        for m in RS.finditer(n):
            bar_runs.add(m.group(0))

targets = sorted(set(asm) - bar_runs)
tset = set(targets)
print("权威源 %d run / barbarum 覆盖 %d run / 待查 %d 个"
      % (len(asm), len(bar_runs), len(targets)))
print()

for tag in ("known_virus_all_v2", "known_virus_all_v3"):
    root = os.path.join(BASE, tag)
    print("=" * 78)
    print("[%s]" % tag)
    print("=" * 78)
    if not os.path.isdir(root):
        print("   目录不存在")
        continue
    present = set()
    for root_d, dirs, files in os.walk(root):
        for n in list(dirs) + list(files):
            for m in RS.finditer(n):
                present.add(m.group(0))
    print("  该目录出现的 run 总数 = %d" % len(present))
    print("  42 个目标中出现 = %d ; 完全不存在 = %d"
          % (len(present & tset), len(tset - present)))

    # 逐阶段
    print("\n  逐阶段命中（仅列出有命中的阶段）:")
    for st in STAGES:
        sp = os.path.join(root, st)
        if not os.path.isdir(sp):
            continue
        got = set()
        nfile = 0
        for r, dirs, files in os.walk(sp):
            for n in list(dirs) + list(files):
                nfile += 1
                for m in RS.finditer(n):
                    got.add(m.group(0))
        hit = got & tset
        if hit:
            print("    %-18s 命中 %2d/%d  (该阶段文件总数 %d)"
                  % (st, len(hit), len(tset), nfile))
    print()

# 结果表层面的最终确认
print("=" * 78)
print("结果表层面（这才是'分析过'的判据）")
print("=" * 78)
for tag in ("known_virus_all_v2", "known_virus_all_v3", "known_virus_barbarum"):
    p = os.path.join(BASE, tag, "01_detection", "summary",
                     "all_viruses.best.summary.tsv")
    if not os.path.isfile(p):
        print("  [%s] 无 all_viruses.best.summary.tsv" % tag)
        continue
    runs = set()
    rows = 0
    with open(p, encoding="utf-8", errors="replace") as fh:
        fh.readline()
        for line in fh:
            if line.strip():
                rows += 1
                for m in RS.finditer(line):
                    runs.add(m.group(0))
    print("  [%s] best.summary 行=%-5d 含run=%-4d  42目标命中=%d"
          % (tag, rows, len(runs), len(runs & tset)))