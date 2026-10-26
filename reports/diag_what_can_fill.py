# -*- coding: utf-8 -*-
"""复核结论落地前的最后一步：确认 42 个缺失样本能补什么。

已知: 42 个 run 在 v2 里只出现于 01_detection。
需确认:
  1. 它们在 v2/01_detection 下的具体文件形态(bam/stat/log/quant)
  2. barbarum 的 01_detection 里对应目录是什么形态(软链 or 真实文件)
  3. 是否有 quantification 结果 —— 有定量才能进 02_filtering 与已知病毒表
只读。
"""
import os
import re
from collections import Counter

# 本文件支持的环境变量：MMPV_SERVER_ROOT、MMPV_PAPER_ROOT（未设时保持原服务器字面量，路径逐字不变）
_SERVER = os.environ.get("MMPV_SERVER_ROOT", "/home/zhangwenda")
_PAPER  = os.environ.get("MMPV_PAPER_ROOT",  f"{_SERVER}/MMPV-paper")
_GOJI   = f"{_PAPER}/goji-virome"


D = f"{_GOJI}/02_novel_virus/RNA-Lycium_barbarum_out/01_Assembly"
V2 = f"{_GOJI}/03_known_virus/known_virus_all_v2/01_detection"
BAR = f"{_GOJI}/03_known_virus/known_virus_barbarum/01_detection"
RUN = re.compile(r"^((?:CRR|SRR|ERR|DRR)[0-9]+)")
RS = re.compile(r"(?:CRR|SRR|ERR|DRR)[0-9]+")

asm = {}
for d in os.listdir(D):
    p = os.path.join(D, d)
    if os.path.isdir(p) and RUN.match(d):
        contig = [f for f in os.listdir(p) if "contig.fasta" in f]
        asm[RUN.match(d).group(0)] = (
            d, os.path.getsize(os.path.join(p, contig[0])) if contig else 0)

mirror = set()
for root, dirs, files in os.walk(BAR):
    for n in list(dirs) + list(files):
        for m in RS.finditer(n):
            mirror.add(m.group(0))

missing = sorted(set(asm) - mirror)
print("权威源 %d / 镜像 %d / 缺失 %d\n" % (len(asm), len(mirror), len(missing)))

# v2 的 01_detection 结构
print("=" * 74)
print("v2/01_detection 顶层结构")
print("=" * 74)
for e in sorted(os.listdir(V2))[:20]:
    p = os.path.join(V2, e)
    kind = "dir " if os.path.isdir(p) else "file"
    n = len(os.listdir(p)) if os.path.isdir(p) else os.path.getsize(p)
    print("   [%s] %-34s %s" % (kind, e, n))

# 缺失 run 在 v2 01_detection 里的文件
print()
print("=" * 74)
print("42 个缺失 run 在 v2/01_detection 下的文件（按后缀统计）")
print("=" * 74)
suffix = Counter()
sample_files = None
for r in missing:
    hits = []
    for root, dirs, files in os.walk(V2):
        for n in list(dirs) + list(files):
            if r in n:
                hits.append(os.path.join(root, n))
    for h in hits:
        suffix[os.path.splitext(h)[1] or "(dir)"] += 1
    if sample_files is None:
        sample_files = sorted(hits)
for s, n in suffix.most_common():
    print("   %-12s %d" % (s, n))
print("\n示例（首个缺失 run 的全部文件）:")
for h in (sample_files or [])[:10]:
    print("   ", os.path.relpath(h, V2))

# barbarum 01_detection 结构对照
print()
print("=" * 74)
print("barbarum/01_detection 顶层结构（对照）")
print("=" * 74)
for e in sorted(os.listdir(BAR))[:20]:
    p = os.path.join(BAR, e)
    kind = "dir " if os.path.isdir(p) else "file"
    real = sum(1 for _, _, fs in os.walk(p) for _ in fs) if os.path.isdir(p) else 1
    link = sum(1 for _, _, fs in os.walk(p) for f in fs
               if os.path.islink(os.path.join(p, f))) if os.path.isdir(p) else 0
    print("   [%s] %-30s 条目数=%d" % (kind, e, len(os.listdir(p)) if os.path.isdir(p) else real))

# 定量文件是否存在（决定能否补进已知病毒表）
print()
print("=" * 74)
print("定量产物核查（quant / counts / tpms）")
print("=" * 74)
for sub in ("quant", "stats", "stat"):
    p = os.path.join(V2, sub)
    if os.path.isdir(p):
        files = os.listdir(p)
        hit = sum(1 for f in files if any(r in f for r in missing))
        print("   v2/%-8s 存在, 缺失run 相关文件 %d/%d" % (sub, hit, len(files)))