# -*- coding: utf-8 -*-
"""复核：权威源 vs v2/v3/barbarum 三方样本集（只读）。

要回答的核心问题：42 个缺失样本是「被镜像过滤器排除了」，
还是「v2/v3 源头就没有」—— 两者的补救方式完全不同。
"""
import os
import re

# 本文件支持的环境变量：MMPV_SERVER_ROOT、MMPV_PAPER_ROOT（未设时保持原服务器字面量，路径逐字不变）
_SERVER = os.environ.get("MMPV_SERVER_ROOT", "/home/zhangwenda")
_PAPER  = os.environ.get("MMPV_PAPER_ROOT",  f"{_SERVER}/MMPV-paper")
_GOJI   = f"{_PAPER}/goji-virome"


D = f"{_GOJI}/02_novel_virus/RNA-Lycium_barbarum_out/01_Assembly"
BASE = f"{_GOJI}/03_known_virus"
RUN = re.compile(r"^(CRR|SRR|ERR|DRR)[0-9]+")
RS = re.compile(r"(?:CRR|SRR|ERR|DRR)[0-9]+")


def scan_runs(root):
    out = set()
    for r, dirs, files in os.walk(root):
        for n in dirs + files:
            for m in RS.finditer(n):
                out.add(m.group(0))
    return out


asm = {}
for d in sorted(os.listdir(D)):
    p = os.path.join(D, d)
    if os.path.isdir(p):
        m = RUN.match(d)
        if m:
            contig = [f for f in os.listdir(p) if "contig.fasta" in f]
            asm[m.group(0)] = (d, os.path.getsize(os.path.join(p, contig[0])) if contig else 0)

v2 = scan_runs(os.path.join(BASE, "known_virus_all_v2"))
v3 = scan_runs(os.path.join(BASE, "known_virus_all_v3"))
bar = scan_runs(os.path.join(BASE, "known_virus_barbarum"))

print("权威源 run           = %d" % len(asm))
print("v2run           = %d" % len(v2))
print("v3    run           = %d" % len(v3))
print("barbarum run          = %d" % len(bar))
print()
print("权威源∩v2            = %d" % len(set(asm) & v2))
print("权威源 ∩ v2 ∩ barbarum   = %d" % len(set(asm) & v2 & bar))
print("v2 有 / 权威源无         = %d" % len(v2 - set(asm)))
print()

miss = sorted(set(asm) - bar)
# 关键判定：这42 个在 v2 里有没有？
in_v2 = [r for r in miss if r in v2]
not_v2 = [r for r in miss if r not in v2]
print("=" * 74)
print("【判定】镜像缺失的 %d 个 run，源头 v2 是否存在：" % len(miss))
print("  v2 里有（=> 被镜像过滤器排除了，需修过滤器）: %d" % len(in_v2))
print("  v2 里也没有（=> 源头就没跑，补救=重跑分析）: %d" % len(not_v2))
print("=" * 74)

def show(title, lst):
    if not lst:
        print("\n%s: 无" % title)
        return
    print("\n%s (%d 个):" % (title, len(lst)))
    for r in lst:
        d, sz = asm[r]
        print("   %-13s %-22s contig=%9d B" % (r, d, sz))

show("A) v2 中存在却被镜像漏掉", in_v2)
show("B) v2 中本就不存在", not_v2)