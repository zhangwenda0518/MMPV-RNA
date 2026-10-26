# -*- coding: utf-8 -*-
"""复核：权威组装源 vs 已知病毒侧 barbarum 的样本集差异（只读，不写）。"""
import os
import re

# 本文件支持的环境变量：MMPV_SERVER_ROOT、MMPV_PAPER_ROOT（未设时保持原服务器字面量，路径逐字不变）
_SERVER = os.environ.get("MMPV_SERVER_ROOT", "/home/zhangwenda")
_PAPER  = os.environ.get("MMPV_PAPER_ROOT",  f"{_SERVER}/MMPV-paper")
_GOJI   = f"{_PAPER}/goji-virome"


D = f"{_GOJI}/02_novel_virus/RNA-Lycium_barbarum_out/01_Assembly"
RUN_RE = re.compile(r"^(CRR|SRR|ERR|DRR)[0-9]+")

asm = {}
for d in sorted(os.listdir(D)):
    p = os.path.join(D, d)
    if not os.path.isdir(p):
        continue
    m = RUN_RE.match(d)
    if m:
        asm[m.group(0)] = d
print("权威源: 样本目录 %d 个, run %d 个" % (len(os.listdir(D)), len(asm)))

TABS = [
    f"{_GOJI}/01_data/gsa.tab",
    f"{_GOJI}/01_data/ncbi.DNA.tab",
    f"{_GOJI}/01_data/ncbi.tab",
    f"{_GOJI}/01_data/ncbi.RNA.tab",
]
runs = set()
for f in TABS:
    if not os.path.isfile(f):
        continue
    with open(f, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            for m in re.finditer(r"(?:CRR|SRR|ERR|DRR)[0-9]+", line):
                runs.add(m.group(0))
print("tab 文件 run 总数 =", len(runs))

only_asm = sorted(set(asm) - runs)
only_tab = sorted(runs - set(asm))
print()
print("权威源有 / tab 无 (潜在新增) =", len(only_asm))
print("tab 有 / 权威源无 (需核实)   =", len(only_tab))

# 已知病毒侧镜像里实际存在的 run
B = f"{_GOJI}/03_known_virus/known_virus_barbarum"
mirror = set()
RUN_SHAPE = re.compile(r"(?:CRR|SRR|ERR|DRR)[0-9]+")
for root, dirs, files in os.walk(B):
    for name in dirs + files:
        for m in RUN_SHAPE.finditer(name):
            mirror.add(m.group(0))
print()
print("barbarum 镜像中出现的 run =", len(mirror))
print("  权威源 ∩ 镜像 =", len(set(asm) & mirror))
print("  权威源 有 / 镜像 无 =", len(set(asm) - mirror))
print("  镜像 有 / 权威源 无 =", len(mirror - set(asm)))