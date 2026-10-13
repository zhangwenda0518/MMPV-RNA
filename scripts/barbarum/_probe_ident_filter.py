#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""_probe_ident_filter.py — 探针：确认 02a_Identification / 02b_Filter 的样本集合与逐工具产物形态。

只读，不写任何产物。
"""
import collections
import glob
import os
import sys
from pathlib import Path

BASE = Path(os.path.expanduser(
    "~/MMPV-paper/goji-virome/02_novel_virus/RNA-Lycium_barbarum_out"))
A = BASE / "02a_Identification"
B = BASE / "02b_Filter"
TOOLS = ["viralverify", "metabuli", "blast", "genomad", "virbot", "rdrpcatch"]

print("BASE =", BASE)
print("BASE.is_symlink() =", BASE.is_symlink(), "->", os.path.realpath(BASE))
print("BASE.parent =", BASE.parent)
print()

a_dirs = sorted(d.name for d in A.iterdir() if d.is_dir())
b_dirs = sorted(d.name for d in B.iterdir() if d.is_dir())
print(f"02a dirs = {len(a_dirs)}    02b dirs = {len(b_dirs)}")

a_id = [s for s in a_dirs if glob.glob(str(A / s / "*virus.all.result.id"))]
print(f"02a 有 *virus.all.result.id = {len(a_id)}")
# 空文件？
nonempty = [s for s in a_id
            if sum(1 for _ in open(A / s / glob.glob(str(A / s / '*virus.all.result.id'))[0])) > 0]
print(f"  其中非空 = {len(nonempty)}")

b_uni = [s for s in b_dirs
         if (B / s / "uniprot" / "all_candidates_uniprot.fasta").is_file()]
print(f"02b 有 uniprot/all_candidates_uniprot.fasta = {len(b_uni)}")
b_uni_nonempty = [s for s in b_uni
                  if (B / s / "uniprot" / "all_candidates_uniprot.fasta").stat().st_size > 2]
print(f"  其中 >2 字节（有保留序列）= {len(b_uni_nonempty)}")

sa, sb = set(a_id), set(b_uni)
print(f"a∩b = {len(sa & sb)}   a-b = {len(sa - sb)}   b-a = {len(sb - sa)}")
if sa - sb:
    print("  a-b 样例:", sorted(sa - sb)[:10])
if sb - sa:
    print("  b-a 样例:", sorted(sb - sa)[:10])
print()

# 逐工具产物存在性
cnt = collections.Counter()
for s in a_id:
    for t in TOOLS:
        if glob.glob(str(A / s / f"*virus.{t}.result.id")):
            cnt[t] += 1
    if glob.glob(str(A / s / "*viroids.result.id")):
        cnt["viroid_blastn"] += 1
print("逐工具 id 文件存在样本数（分母 =", len(a_id), "）:")
for t in TOOLS + ["viroid_blastn"]:
    print(f"  {t:16s} {cnt.get(t, 0)}")
print()

# 逐样本 union 与各工具大小 + 保留数（抽 3 个样本）
def read_ids(p):
    with open(p) as f:
        return set(x.strip() for x in f if x.strip())

print("抽样 3 个样本:")
for s in a_id[:3]:
    ad = A / s
    per = {}
    for t in TOOLS:
        g = glob.glob(str(ad / f"*virus.{t}.result.id"))
        per[t] = len(read_ids(g[0])) if g else None
    gv = glob.glob(str(ad / "*viroids.result.id"))
    viroid = read_ids(gv[0]) if gv else set()
    union = read_ids(glob.glob(str(ad / "*virus.all.result.id"))[0])
    toolsets = {t: read_ids(glob.glob(str(ad / f"*virus.{t}.result.id"))[0])
                for t in TOOLS if glob.glob(str(ad / f"*virus.{t}.result.id"))}
    recomputed = set().union(*toolsets.values()) if toolsets else set()
    up = B / s / "uniprot" / "all_candidates_uniprot.fasta"
    retained = set()
    if up.is_file():
        with open(up) as f:
            for line in f:
                if line.startswith(">"):
                    retained.add(line[1:].split()[0])
    cdd = B / s / "cdd" / "all_candidates_cdd_filter.fasta"
    cdd_ids = set()
    if cdd.is_file():
        with open(cdd) as f:
            for line in f:
                if line.startswith(">"):
                    cdd_ids.add(line[1:].split()[0])
    print(f"  [{s}] 各工具={per}  viroid={len(viroid)}")
    print(f"       union(id文件)={len(union)}  重算并集={len(recomputed)}  相等={union == recomputed}")
    print(f"       viroid∩union = {len(viroid & union)}")
    print(f"       cdd_kept={len(cdd_ids)}  union-cdd = {len(union - cdd_ids)}  cdd-union = {len(cdd_ids - union)}")
    print(f"       retained={len(retained)}  retained⊆cdd = {retained <= cdd_ids}  retained⊆union = {retained <= union}")
    print(f"       cdd - retained = {len(cdd_ids - retained)}  (union - retained) = {len(union - retained)}")
print()

# 全域汇总
tot_union = 0
tot_viroid = 0
tot_retained = 0
tot_cdd = 0
bad = 0
for s in a_id:
    ad = A / s
    g = glob.glob(str(ad / "*virus.all.result.id"))
    tot_union += len(read_ids(g[0]))
    gv = glob.glob(str(ad / "*viroids.result.id"))
    if gv:
        tot_viroid += len(read_ids(gv[0]))
    up = B / s / "uniprot" / "all_candidates_uniprot.fasta"
    if up.is_file():
        with open(up) as f:
            tot_retained += sum(1 for line in f if line.startswith(">"))
    cdd = B / s / "cdd" / "all_candidates_cdd_filter.fasta"
    if cdd.is_file():
        with open(cdd) as f:
            tot_cdd += sum(1 for line in f if line.startswith(">"))
print(f"全域汇总: Σunion={tot_union}  Σviroid={tot_viroid}  Σretained={tot_retained}  Σcdd_kept={tot_cdd}")
print(f"  （对照 ident_summary TOTAL=556773, viroid=604, filter_per_sample Uniprot_Filtered=94826）")
