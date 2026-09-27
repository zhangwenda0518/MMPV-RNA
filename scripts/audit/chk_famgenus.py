#!/usr/bin/env python3
"""科-属校准原型：用标准化的逐工具谱系量出「补这条校准会动多少行」

规则（与 harmonize_genus_species 同风格，纯内生一致性，不引外部库）：
  对每行取共识 Family = F*
  S = 自己报的 Family == F* 的那些工具
  候选 = S 里各工具报的 Genus（按 build_consensus 的权重 64*TOOL_BIAS 计票）
  若共识 Genus 不在候选里 -> 用候选里票最高的替换；候选为空则该行无可替换项
"""
import csv
import importlib.util
from collections import Counter, defaultdict
from pathlib import Path

BASE = Path("/home/zhangwenda/MMPV-paper/goji-virome/02_novel_virus/"
            "RNA-Lycium_barbarum_out/05_Taxonomy/Votus.integrated")
PIPE = Path("/home/zhangwenda/MMPV-RNA/virome_discovery_pipeline")
VMR = Path("/home/zhangwenda/plant_virus_db/1.virus-host_db/B-ictv/VMR_MSL41.tsv")

TOOLS = ["ACVirus", "CAT", "diamond_lca", "genomad", "metabuli", "mmseqs", "VITAP"]
BIAS = {"ACVirus": 1.2, "VITAP": 1.1, "mmseqs": 1.0, "metabuli": 1.0,
        "CAT": 0.9, "genomad": 0.9, "diamond_lca": 0.8}
GW = {t: 64.0 * BIAS.get(t, 0.8) for t in TOOLS}
RANKS = ["Realm", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]

spec = importlib.util.spec_from_file_location("ann", PIPE / "utils" / "annotate_nucleic_acid.py")
ann = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ann)
idx = ann.load_vmr(VMR)

# ---- 载入逐工具标准化表 ----
per_tool = defaultdict(dict)   # contig_id -> {tool: {rank: value}}
for t in TOOLS:
    p = BASE / ("standardized_%s.tsv" % t)
    if not p.is_file():
        print("缺 %s" % p.name)
        continue
    with open(p, newline="", encoding="utf-8", errors="replace") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            per_tool[r.get("contig_id")][t] = r
print("已载入 %d 条 contig 的逐工具谱系（%d 个工具）" % (len(per_tool), len(TOOLS)))

with open(BASE / "final_integrated_classification.tsv", newline="", encoding="utf-8",
          errors="replace") as f:
    rows = list(csv.DictReader(f, delimiter="\t"))
print("共识表 %d 行" % len(rows))

B = ann.is_blank
V = lambda v: not B(v)

# ============================================================
# Part A: 逐个工具摊开那条争议 contig
# ============================================================
tgt = next((r["contig_id"] for r in rows
            if r["contig_id"].startswith("CRR1126135_clean_NODE_68_length_2879")), None)
if tgt:
    print()
    print("=" * 78)
    print("Part A: %s" % tgt)
    print("=" * 78)
    print("  %-13s %-14s %-16s %-14s %-18s %s" % ("Tool", "Family", "Genus", "Order", "Class", "Species"))
    for t in TOOLS:
        d = per_tool.get(tgt, {}).get(t)
        if not d:
            print("  %-13s (该工具未输出此 contig)" % t)
            continue
        print("  %-13s %-14s %-16s %-14s %-18s %s" % (
            t, d.get("Family", "") or "-", d.get("Genus", "") or "-",
            d.get("Order", "") or "-", d.get("Class", "") or "-",
            (d.get("Species", "") or "-")[:28]))
    print("  工具报出的属: %s" % {
        t: (per_tool[tgt][t].get("Genus") or "-") for t in TOOLS if t in per_tool.get(tgt, {})})

# ============================================================
# Part B: 全校准影响
# ============================================================
st = Counter()
fix_rows, blank_rows, keep_rows = [], [], []
for r in rows:
    cid = r["contig_id"]
    F, G = r.get("Family", ""), r.get("Genus", "")
    if not V(F) or not V(G):
        st["无科或无属, 不动"] += 1
        continue
    tl = per_tool.get(cid, {})
    if not tl:
        st["逐工具表里没有这条 contig"] += 1
        continue
    S = [t for t in tl if (tl[t].get("Family") or "") == F]
    if not S:
        st["没有任何工具报这个科(异常)"] += 1
        continue
    cand = Counter()
    for t in S:
        g = tl[t].get("Genus")
        if V(g):
            cand[g] += GW[t]
    if G in cand:
        st["属已自洽"] += 1
        continue
    if cand:
        st["属与科不一致, 有替换候选"] += 1
        fix_rows.append((r, cand.most_common(1)[0][0]))
    else:
        st["属与科不一致, 无替换候选"] += 1
        blank_rows.append(r)

print()
print("=" * 78)
print("Part B: 校准统计（只针对科和属都非空的行）")
print("=" * 78)
for k, v in st.most_common():
    print("  %-30s %6d" % (k, v))

print()
print("  会改属的合计 : %d 行 (%5.2f%%)" % (len(fix_rows) + len(blank_rows),
                                       100.0 * (len(fix_rows) + len(blank_rows)) / len(rows)))
print("    ├ 换成家族内其它工具报的属 : %d" % len(fix_rows))
print("    └ 家族内无人报属 -> 置空   : %d" % len(blank_rows))

# ============================================================
# Part C: 对 Nucleic_acid 列的影响
# ============================================================
def rep(row, newg):
    d = dict(row)
    d["Genus"] = newg
    return d


def col(rows_iter, mutate=None):
    out = []
    for r in rows_iter:
        rr = rep(r, mutate(r)) if mutate else r
        b, _ = ann.decide(rr, idx)
        out.append(b)
    return out


fix_map = {r["contig_id"]: g for r, g in fix_rows}
blank_ids = {r["contig_id"] for r in blank_rows}

before = col(rows)
after_R1 = col(rows, lambda r: fix_map.get(r["contig_id"],
                                           ("NA" if r["contig_id"] in blank_ids else r.get("Genus", ""))))
after_R2 = col(rows, lambda r: fix_map.get(r["contig_id"], r.get("Genus", "")))

print()
print("=" * 78)
print("Part C: 对 Nucleic_acid 的影响")
print("=" * 78)
for nm, vec in [("校准前", before), ("R1 替换或置空", after_R1), ("R2 只替换不置空", after_R2)]:
    c = Counter(v or "NA" for v in vec)
    print("  %-16s DNA %6d  RNA %6d  NA %5d" % (nm, c["DNA"], c["RNA"], c["NA"]))
for nm, vec in [("R1 替换或置空", after_R1), ("R2 只替换不置空", after_R2)]:
    d = sum(1 for x, y in zip(before, vec) if x != y)
    dist = Counter("%s->%s" % (x or "NA", y or "NA") for x, y in zip(before, vec) if x != y)
    print("  %-16s 与校准前分歧 %4d 行  %s" % (nm, d, dict(dist)))

# 306 行打架行
clash = [r for r in rows
         if ann.resolve_rank(r["Family"], "Family", idx) and ann.resolve_rank(r["Genus"], "Genus", idx)
         and ann.resolve_rank(r["Family"], "Family", idx) != ann.resolve_rank(r["Genus"], "Genus", idx)]
print()
print("  306 行科属打架行里：")
fixed = sum(1 for r in clash if r["contig_id"] in fix_map)
blanked = sum(1 for r in clash if r["contig_id"] in blank_ids)
print("    被替换属的 %d 行, 被置空属的 %d 行, 未触及 %d 行" % (
    fixed, blanked, len(clash) - fixed - blanked))
