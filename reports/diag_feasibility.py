# -*- coding: utf-8 -*-
"""可行性核查：补齐 42 个样本的 01_detection 需要什么输入、是否具备。

结论依据(已实测):
  - batch_virus_depth.py -i要的是 FASTQ/FASTA reads 目录（不是组装 contig）
  - v2/01_detection/bam 里 .bam 本体 = 0，只有 1156 个 .bai 死索引
  - 所有 00a_CleanData / 01_clean 目录都是空的（已清理）
  -> 唯一可行输入是 clean reads 库 /home/zhangwenda/sra_public_orders/fq_flat

本脚本只读核实：42 个目标样本的 reads 是否在 fq_flat 里、各自多大。
"""
import os
import re

# 本文件支持的环境变量：MMPV_SERVER_ROOT、MMPV_PAPER_ROOT（未设时保持原服务器字面量，路径逐字不变）
_SERVER = os.environ.get("MMPV_SERVER_ROOT", "/home/zhangwenda")
_PAPER  = os.environ.get("MMPV_PAPER_ROOT",  f"{_SERVER}/MMPV-paper")
_GOJI   = f"{_PAPER}/goji-virome"


FQ = f"{_SERVER}/sra_public_orders/fq_flat"
D = f"{_GOJI}/02_novel_virus/RNA-Lycium_barbarum_out/01_Assembly"
B = f"{_GOJI}/03_known_virus/known_virus_barbarum"

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

# 建 fq_flat 索引：run -> [文件名...]
idx = {}
if os.path.isdir(FQ):
    for f in os.listdir(FQ):
        for m in RS.finditer(f):
            idx.setdefault(m.group(0), []).append(f)

print("=" * 76)
print("fq_flat 索引: %d 个文件,覆盖 %d 个 run" % (len(os.listdir(FQ)), len(idx)))
print("缺失目标:%d 个" % len(missing))
print("=" * 76)

has_fq, no_fq = [], []
for r in missing:
    (has_fq if r in idx else no_fq).append(r)

print("\n【A】clean reads 在 fq_flat 里可用: %d 个" % len(has_fq))
tot = 0
for r in has_fq:
    files = idx[r]
    sz = sum(os.path.getsize(os.path.join(FQ, f)) for f in files)
    tot += sz
    print("   %-13s contig=%9d B  fq=%d 文件 %7.1f MB" % (
        r, asm[r][1], len(files), sz / 1e6))
print("   合计 %.1f GB" % (tot / 1e9))

print("\n【B】clean reads 缺失(需重新从 SRA/ENA 下载): %d 个" % len(no_fq))
for r in no_fq[:50]:
    print("   %-13s contig=%9d B" % (r, asm[r][1]))

# 已有结果里的样本是否都在 fq_flat（验证该库是否就是v2 当初的输入）
ok = sum(1 for r in mirror if r in idx)
print("\n【交叉验证】现有 %d 个已覆盖 run 中，%d 个在 fq_flat 里有 reads" % (
    len(mirror), ok))