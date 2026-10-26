# -*- coding: utf-8 -*-
"""逐 run 精确核实：v2 / v3 到底有没有这些样本。

用户问「这 46 个全都缺失嘛」——需要区分三种状态并给逐条证据：
  (a) 目录里完全没有任何该 run 的文件      => 真缺失
  (b) 只有中间产物(bai/log/stat), 无结果行  => 跑过但没检出
  (c) 结果表 best.summary 里有行           => 真分析过
只读，不改任何数据。
"""
import csv
import os
import re

# 本文件支持的环境变量：MMPV_SERVER_ROOT、MMPV_PAPER_ROOT（未设时保持原服务器字面量，路径逐字不变）
_SERVER = os.environ.get("MMPV_SERVER_ROOT", "/home/zhangwenda")
_PAPER  = os.environ.get("MMPV_PAPER_ROOT",  f"{_SERVER}/MMPV-paper")
_GOJI   = f"{_PAPER}/goji-virome"


D = f"{_GOJI}/02_novel_virus/RNA-Lycium_barbarum_out/01_Assembly"
BASE = f"{_GOJI}/03_known_virus"
BAR = os.path.join(BASE, "known_virus_barbarum")

RUN = re.compile(r"^((?:CRR|SRR|ERR|DRR)[0-9]+)")
RS = re.compile(r"(?:CRR|SRR|ERR|DRR)[0-9]+")

# 1) 权威源
asm = {}
for d in os.listdir(D):
    p = os.path.join(D, d)
    if os.path.isdir(p) and RUN.match(d):
        c = [f for f in os.listdir(p) if "contig.fasta" in f]
        asm[RUN.match(d).group(0)] = (
            d, os.path.getsize(os.path.join(p, c[0])) if c else 0)

# 2) barbarum 覆盖
bar = set()
for root, dirs, files in os.walk(BAR):
    for n in list(dirs) + list(files):
        for m in RS.finditer(n):
            bar.add(m.group(0))

targets = sorted(set(asm) - bar)
tset = set(targets)

print("=" * 80)
print("数字核对")
print("=" * 80)
print("  权威源 run           = %d" % len(asm))
print("  barbarum 覆盖 run= %d" % len(bar))
print("  差值(权威源 - barbarum) = %d" % (len(asm) - len(bar)))
print("  其中在权威源却不在 barbarum = %d" % len(targets))
print("  barbarum 有但权威源无        = %d" % len(bar - set(asm)))
print()
print("  >>> 你说的 46 与我算的 %d 不一致。差异来源有两种可能：" % len(targets))
print("      (1) 权威源 446 个目录里有 1 个不是 run 形态 -> 实际 run 数 %d" % len(asm))
print("      (2) 之前某次统计口径不同（含 scaffolds 而非 _clean 目录等）")
print("  >>> 权威源 446 个目录的形态核查：")
odd = [d for d in os.listdir(D)
       if os.path.isdir(os.path.join(D, d)) and not RUN.match(d)]
for d in odd:
    print("        非 run 形态目录: %s" % d)
print("        共 %d 个" % len(odd))
print()

# 3) 逐 run 在 v2 / v3 的状态
def scan(root):
    """返回 {run: set(相对路径第一层)}"""
    out = {}
    for r, dirs, files in os.walk(root):
        for n in list(dirs) + list(files):
            for m in RS.finditer(n):
                rel = os.path.relpath(os.path.join(r, n), root)
                out.setdefault(m.group(0), set()).add(rel.split(os.sep)[0])
    return out


def best_runs(root):
    p = os.path.join(root, "01_detection", "summary",
                     "all_viruses.best.summary.tsv")
    if not os.path.isfile(p):
        return set(), 0
    runs, rows = set(), 0
    with open(p, encoding="utf-8", errors="replace") as fh:
        fh.readline()
        for line in fh:
            if line.strip():
                rows += 1
                for m in RS.finditer(line):
                    runs.add(m.group(0))
    return runs, rows


v2 = scan(os.path.join(BASE, "known_virus_all_v2"))
v3 = scan(os.path.join(BASE, "known_virus_all_v3"))
v2b, v2rows = best_runs(os.path.join(BASE, "known_virus_all_v2"))
v3b, v3rows = best_runs(os.path.join(BASE, "known_virus_all_v3"))

print("=" * 80)
print("v2 状态: 出现 run=%d | best.summary 行=%d 含run=%d"
      % (len(v2), v2rows, len(v2b)))
print("v3 状态: 出现 run=%d | best.summary 行=%d 含run=%d"
      % (len(v3), v3rows, len(v3b)))
print("=" * 80)

print("\n%-14s %10s | %-28s | %-28s | %s"
      % ("run", "contig", "v2", "v3", "best.summary(v2/v3)"))
print("-" * 108)
cnt = {"v2_none": 0, "v2_mid": 0, "v2_best": 0,
       "v3_none": 0, "v3_mid": 0, "v3_best": 0}
for r in targets:
    sz = asm[r][1]
    s2 = sorted(v2.get(r, [])) or ["— 完全缺失 —"]
    s3 = sorted(v3.get(r, [])) or ["— 完全缺失 —"]
    in2 = r in v2b
    in3 = r in v3b
    if r not in v2:
        cnt["v2_none"] += 1
    elif in2:
        cnt["v2_best"] += 1
    else:
        cnt["v2_mid"] += 1
    if r not in v3:
        cnt["v3_none"] += 1
    elif in3:
        cnt["v3_best"] += 1
    else:
        cnt["v3_mid"] += 1
    print("%-14s %10d | %-28s | %-28s | %s / %s"
          % (r, sz, ",".join(s2)[:27], ",".join(s3)[:27],
             "有" if in2 else "无", "有" if in3 else "无"))

print("-" * 108)
print("\n汇总（%d 个目标）:" % len(targets))
print("  v2: 完全缺失 %d | 仅有中间产物 %d | best.summary 有行 %d"
      % (cnt["v2_none"], cnt["v2_mid"], cnt["v2_best"]))
print("  v3: 完全缺失 %d | 仅有中间产物 %d | best.summary 有行 %d"
      % (cnt["v3_none"], cnt["v3_mid"], cnt["v3_best"]))