# -*- coding: utf-8 -*-
"""干跑 + 执行：把 v2 里 42 个缺失样本的 01_detection 中间产物整合进 barbarum 镜像。

映射关系（v2 与 barbarum 的 01_detection 同构，扁平结构）：
    bam/<run>.sorted.bam.bai
    logs/sample_<run>.log
    stat/<run>.SiteDepth.gz
    stat/<run>.chr.stat.gz
-> 全部建为**绝对路径 symlink**，与 barbarum 现有风格一致
   （实测现有 link 指向 v2 的绝对路径）。

安全约束：
  1. 只在 barbarum 侧新建 symlink，**绝不写 v2**（避免之前的写穿污染）
  2. 已存在的条目一律跳过（幂等），不覆盖
  3. 干跑模式只打印计划，不落盘
"""
import os
import sys

# 本文件支持的环境变量：MMPV_SERVER_ROOT、MMPV_PAPER_ROOT（未设时保持原服务器字面量，路径逐字不变）
_SERVER = os.environ.get("MMPV_SERVER_ROOT", "/home/zhangwenda")
_PAPER  = os.environ.get("MMPV_PAPER_ROOT",  f"{_SERVER}/MMPV-paper")
_GOJI   = f"{_PAPER}/goji-virome"


DRY = "--apply" not in sys.argv

D = f"{_GOJI}/02_novel_virus/RNA-Lycium_barbarum_out/01_Assembly"
V2_DET = f"{_GOJI}/03_known_virus/known_virus_all_v2/01_detection"
BAR_DET = f"{_GOJI}/03_known_virus/known_virus_barbarum/01_detection"

import re
RUN = re.compile(r"^((?:CRR|SRR|ERR|DRR)[0-9]+)")
RS = re.compile(r"(?:CRR|SRR|ERR|DRR)[0-9]+")

asm = set()
for d in os.listdir(D):
    p = os.path.join(D, d)
    if os.path.isdir(p) and RUN.match(d):
        asm.add(RUN.match(d).group(0))

bar = set()
for root, dirs, files in os.walk(BAR_DET):
    for n in list(dirs) + list(files):
        for m in RS.finditer(n):
            bar.add(m.group(0))

targets = sorted(asm - bar)
print("权威源 %d run | barbarum %d run | 待补 %d" % (len(asm), len(bar), len(targets)))
print("模式: %s" % ("执行" if not DRY else "干跑(不写盘)"))
print()

# 建立 42 个目标 run 在 v2 里的实际文件清单（以 v2 为准，不靠猜文件名）
SUBDIRS = ["bam", "logs", "stat"]
plan = []          # (src, dst)
missing_in_v2 = [] # run 在 v2 里一个文件都没有

for r in targets:
    found = False
    for sd in SUBDIRS:
        sdpath = os.path.join(V2_DET, sd)
        if not os.path.isdir(sdpath):
            continue
        for fn in os.listdir(sdpath):
            if r in fn:
                plan.append((os.path.join(sdpath, fn),
                             os.path.join(BAR_DET, sd, fn)))
                found = True
    if not found:
        missing_in_v2.append(r)

print("=" * 74)
print("计划整合条目: %d 个（来自 %d 个 run）" % (len(plan), len(targets) - len(missing_in_v2)))
print("=" * 74)

# 按子目录分类展示 + 检测冲突
from collections import Counter
per_sd = Counter()
todo, skip_exist, broken_src = [], [], []
for src, dst in plan:
    sd = os.path.relpath(src, V2_DET).split(os.sep)[0]
    per_sd[sd] += 1
    if not os.path.isfile(src):
        broken_src.append(src)
    elif os.path.lexists(dst):
        skip_exist.append(dst)
    else:
        todo.append((src, dst))

for sd, n in per_sd.most_common():
    print("  %-8s %d 个" % (sd, n))
print("\n  将新建   : %d" % len(todo))
print("  已存在跳过: %d" % len(skip_exist))
print("  源缺失   : %d" % len(broken_src))
if missing_in_v2:
    print("  v2 中无任何文件的 run: %d %s" % (len(missing_in_v2), missing_in_v2[:5]))

if broken_src:
    print("\n[错误] 以下源不存在，无法建链:")
    for s in broken_src[:5]:
        print("   ", s)
    sys.exit(2)

if DRY:
    print("\n[干跑] 样例（前 8 条）:")
    for src, dst in todo[:8]:
        print("   %s\n     -> %s" % (src.replace(V2_DET + "/", "v2/"),
                                     dst.replace(BAR_DET + "/", "barbarum/")))
    print("\n[干跑] 未落盘。确认无误后加 --apply 执行。")
    sys.exit(0)

# ── 执行：先全部校验通过再写 ──
for src, dst in todo:
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    os.symlink(src, dst)

print("\n[完成] 新建 symlink %d 个" % len(todo))

# ── 验证：链接是否可用（不能是断链）──
ok = broken = 0
for src, dst in todo:
    if os.path.isfile(dst):
        ok += 1
    else:
        broken += 1
print("[验证] 可解析 %d / 断链 %d" % (ok, broken))