#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""诊断：cascade 逐级淘汰之后，为什么行内 Family 与 Genus 仍会互不配套。

输入（闸门关闭的 cascade 产物）：
    /tmp/pregate_cascade/final_integrated_classification.tsv   裸共识结果
    /tmp/pregate_cascade/standardized_<TOOL>.tsv               7 个工具各自的 8 阶元取值
    参照表 ~/database/taxonomy/genus_family_ref.tsv

判据：
    「自洽」= 存在某个工具同时报出该行的 Family 与该行的 Genus。
    不自洽行按成因分类，回答「科是这一家定的、属是那一家定的」怎么发生。
    参照表侧按小写属名匹配（R 层同口径），区分「双参照真矛盾」「单侧」「参照不可判」。
"""
import csv
import glob
import os
from collections import Counter

OUT = os.environ.get("DIAG_OUT", "/tmp/pregate_cascade")
REF = os.path.expanduser("~/database/taxonomy/genus_family_ref.tsv")
EMPTY = {"", "na", "n/a", "nan", "-", "none", "null"}


def norm(x):
    x = (x or "").strip().strip('"').strip()
    return "" if x.lower() in EMPTY else x


def fam_set(x):
    """VMR 侧用 ;科; 包裹以承载一属多科。"""
    x = norm(x)
    if not x:
        return set()
    return {p.strip() for p in x.split(";") if p.strip()}


# ---------- 参照表：genus(小写) -> (科集合, 来源标记) ----------
ref = {}
with open(REF, newline="", encoding="utf-8", errors="replace") as f:
    rd = csv.DictReader(f, delimiter="\t")
    for r in rd:
        g = norm(r.get("Genus")).lower()
        if not g:
            continue
        ncbi = fam_set(r.get("NCBI_Family"))
        vmr = fam_set(r.get("VMR_Family"))
        if not ncbi and not vmr:
            continue
        prev = ref.get(g)
        if prev is None:
            ref[g] = (ncbi, vmr)
        else:
            ref[g] = (prev[0] | ncbi, prev[1] | vmr)
print("参照表属数（小写去重）%d" % len(ref))

# ---------- 逐工具表 ----------
tools = {}
for p in sorted(glob.glob(os.path.join(OUT, "standardized_*.tsv"))):
    name = os.path.basename(p)[len("standardized_"):-len(".tsv")]
    d = {}
    with open(p, newline="", encoding="utf-8", errors="replace") as f:
        rd = csv.DictReader(f, delimiter="\t")
        for r in rd:
            d[norm(r["contig_id"])] = (norm(r.get("Family")), norm(r.get("Genus")))
    tools[name] = d
names = sorted(tools)
print("工具数 %d: %s" % (len(names), ", ".join(names)))

# ---------- 裸共识结果 ----------
p = os.path.join(OUT, "final_integrated_classification.tsv")
rows = list(csv.DictReader(open(p, newline="", encoding="utf-8", errors="replace"), delimiter="\t"))
print("裸共识行数 %d" % len(rows))

stat = Counter()
ref_stat = Counter()
fam_ct = Counter()
examples = []
for r in rows:
    cid = norm(r["contig_id"])
    fam, gen = norm(r.get("Family")), norm(r.get("Genus"))
    if not fam and not gen:
        stat["科属全空"] += 1
        continue
    if not fam:
        stat["有属无科"] += 1
        continue
    if not gen:
        stat["有科无属"] += 1
        continue
    fam_tools = [t for t in names if tools[t].get(cid, ("", ""))[0] == fam]
    gen_tools = [t for t in names if tools[t].get(cid, ("", ""))[1] == gen]
    both = sorted(set(fam_tools) & set(gen_tools))
    if both:
        stat["自洽(存在单一工具同时报出这对)"] += 1
        rf = ref.get(gen.lower())
        if rf is None:
            ref_stat["自洽行 属不在参照表"] += 1
        elif fam.lower() in {x.lower() for x in rf[0]} and fam.lower() in {x.lower() for x in rf[1]}:
            ref_stat["自洽行 双参照都相容"] += 1
        else:
            ref_stat["自洽行 参照侧另有说法"] += 1
        continue

    stat["不配套"] += 1
    fam_ct[fam] += 1
    # 参照侧裁决力
    rf = ref.get(gen.lower())
    if rf is None:
        ref_stat["不配套行 属不在参照表(参照不可判)"] += 1
    else:
        fl = fam.lower()
        in_n = fl in {x.lower() for x in rf[0]}
        in_v = fl in {x.lower() for x in rf[1]}
        if in_n and in_v:
            ref_stat["不配套行 参照说行内科是对的(属另有定型科但不排斥)"] += 1
        elif in_n or in_v:
            ref_stat["不配套行 参照单侧相容"] += 1
        elif rf[0] and rf[1]:
            ref_stat["不配套行 双参照一致反驳行内科"] += 1
        else:
            ref_stat["不配套行 参照单侧且有反驳"] += 1
    # 成因
    if not fam_tools:
        cat = "科无任何工具支持(兜底填出)"
    elif not gen_tools:
        cat = "属无任何工具支持(填充或由种提属)"
    else:
        silent = [t for t in gen_tools if not tools[t].get(cid, ("", ""))[0]]
        if len(silent) == len(gen_tools):
            cat = "沉默逃逸: 报属的工具在科级全无值"
        elif silent:
            cat = "部分沉默: 报属的工具中部分在科级无值"
        else:
            cat = "各阶元赢家不同工具: 报属的工具全都报了别的科"
    stat[cat] += 1
    if len(examples) < 6:
        gen_fams = sorted({tools[t].get(cid, ("", ""))[0] for t in gen_tools} - {""})
        examples.append((cid, fam, fam_tools, gen, gen_tools, gen_fams, cat,
                         (sorted(rf[0]) if rf else []), (sorted(rf[1]) if rf else [])))

print()
for k, v in stat.most_common():
    print("  %-42s %6d" % (k, v))
print()
for k, v in ref_stat.most_common():
    print("  [参照侧] %-46s %6d" % (k, v))
print()
print("不配套行的共识科分布 top10:")
for k, v in fam_ct.most_common(10):
    print("  %-32s %5d" % (k, v))
print()
print("样例（contig | 共识科 <- 报科工具 | 共识属 <- 报属工具 | 报属者自报的科 | 参照 NCBI / VMR | 成因）")
for cid, fam, ft, gen, gt, gf, cat, rn, rv in examples:
    print("  %s\n    科=%s <- %s\n    属=%s <- %s (报属者自报科=%s)\n    参照 NCBI=%s VMR=%s | %s"
          % (cid, fam, ",".join(ft), gen, ",".join(gt) or "∅", ",".join(gf) or "∅",
             ",".join(rn) or "∅", ",".join(rv) or "∅", cat))
