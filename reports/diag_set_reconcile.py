# -*- coding: utf-8 -*-
"""复核数字对不上的部分：权威源∩v2=404 但 barbarum=403，缺失却是 42。

三个集合的口径不同：
  - scan_runs() 在整个目录树上抓**任何**含 run 串的文件/目录名（含日志名 sample_XXX.log）
  - 权威源集合只认<Run>_clean 目录
  - barbarum 镜像里的 run 分布在多个阶段，可能同名不同层级

需要逐 run 定位：它在哪一层被计入？为何有的进镜像有的没进？
只读诊断。
"""
import os
import re

# 本文件支持的环境变量：MMPV_SERVER_ROOT、MMPV_PAPER_ROOT（未设时保持原服务器字面量，路径逐字不变）
_SERVER = os.environ.get("MMPV_SERVER_ROOT", "/home/zhangwenda")
_PAPER  = os.environ.get("MMPV_PAPER_ROOT",  f"{_SERVER}/MMPV-paper")
_GOJI   = f"{_PAPER}/goji-virome"


D = f"{_GOJI}/02_novel_virus/RNA-Lycium_barbarum_out/01_Assembly"
B = f"{_GOJI}/03_known_virus/known_virus_barbarum"
RUN = re.compile(r"^((?:CRR|SRR|ERR|DRR)[0-9]+)")
RS = re.compile(r"(?:CRR|SRR|ERR|DRR)[0-9]+")

asm = set()
for d in os.listdir(D):
    p = os.path.join(D, d)
    if os.path.isdir(p) and RUN.match(d):
        asm.add(RUN.match(d).group(0))

# 记录每个 run 在镜像里"出现于哪些相对路径的前两层"
where = {}
for root, dirs, files in os.walk(B):
    for n in list(dirs) + list(files):
        for m in RS.finditer(n):
            r = m.group(0)
            rel = os.path.relpath(os.path.join(root, n), B)
            parts = rel.split(os.sep)
            stage = parts[0] if parts else "?"
            where.setdefault(r, set()).add(stage)

mirror = set(where)
print("权威源 =", len(asm))
print("镜像中出现的 run =", len(mirror))
print("交集 =", len(asm & mirror))
print("镜像有/权威源无 =", len(mirror - asm))
print()

only_asm = sorted(asm - mirror)
# 对每个缺失 run，看它在 v2 里位于哪个阶段 —— 判断是否"阶段性缺失"
print("=" * 76)
print("42 个缺失 run 在 v2 中的分布（判断是全阶段缺还是个别阶段缺）")
print("=" * 76)
V2 = f"{_GOJI}/03_known_virus/known_virus_all_v2"
RS2 = re.compile(r"(?:CRR|SRR|ERR|DRR)[0-9]+")
v2where = {}
for root, dirs, files in os.walk(V2):
    for n in list(dirs) + list(files):
        for m in RS2.finditer(n):
            rel = os.path.relpath(os.path.join(root, n), V2)
            v2where.setdefault(m.group(0), set()).add(rel.split(os.sep)[0])

from collections import Counter
cnt = Counter()
for r in only_asm:
    stages = v2where.get(r, set())
    for s in stages:
        cnt[s] += 1
print("\n缺失 run 在 v2 里出现过的阶段（计数, 一个 run 可出现在多个阶段）:")
for s, n in cnt.most_common():
    print("   %-22s %d" % (s, n))
print("\n其中 v2 里完全找不到的:", sum(1 for r in only_asm if r not in v2where))

print("\n明细（前 45）: run | v2 阶段")
for r in only_asm[:45]:
    print("   %-13s %s" % (r, sorted(v2where.get(r, {"(v2 无)"}))))