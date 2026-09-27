#!/usr/bin/env python3
"""摊出 final_integrated_classification.tsv 里「一行两个分类」的原始行

回答：这张表是不是每条 contig 只有一个分类？
"""
import csv
import importlib.util
from collections import Counter
from pathlib import Path

PIPE = Path("/home/zhangwenda/MMPV-RNA/virome_discovery_pipeline")
VMR = Path("/home/zhangwenda/plant_virus_db/1.virus-host_db/B-ictv/VMR_MSL41.tsv")
SRC = Path("/home/zhangwenda/MMPV-paper/goji-virome/02_novel_virus/"
           "RNA-Lycium_barbarum_out/05_Taxonomy/Votus.integrated/"
           "final_integrated_classification.tsv")

spec = importlib.util.spec_from_file_location("ann", PIPE / "utils" / "annotate_nucleic_acid.py")
ann = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ann)
idx = ann.load_vmr(VMR)
RANKS = ann.RANKS

with open(SRC, newline="", encoding="utf-8", errors="replace") as f:
    rd = csv.DictReader(f, delimiter="\t")
    HEAD = rd.fieldnames
    rows = list(rd)

print("表头 (%d 列): %s" % (len(HEAD), HEAD))
print("行数: %d" % len(rows))

# ---- 科属两层直接打架 ----
both = [r for r in rows
        if ann.resolve_rank(r["Family"], "Family", idx) and ann.resolve_rank(r["Genus"], "Genus", idx)]
clash = [r for r in both
         if ann.resolve_rank(r["Family"], "Family", idx) != ann.resolve_rank(r["Genus"], "Genus", idx)]
print()
print("=== 科与属都判得出大类的行 ===")
print("  两层都可判 : %d" % len(both))
print("  两层打架   : %d (%5.2f%% of all, %5.2f%% of 两层可判)" % (
    len(clash), 100.0 * len(clash) / len(rows), 100.0 * len(clash) / max(len(both), 1)))
print("  打架形态   : %s" % dict(Counter(
    "%s(科)->%s(属)" % (ann.resolve_rank(r["Family"], "Family", idx),
                        ann.resolve_rank(r["Genus"], "Genus", idx)) for r in clash)))
print("  涉及的科   : %s" % dict(Counter(r["Family"] for r in clash).most_common(6)))


def dump(r, tag):
    print()
    print("  ── [%s] ──" % tag)
    for k in HEAD:
        v = r[k]
        extra = ""
        if k in RANKS and v and v != "NA":
            b = ann.resolve_rank(v, k, idx)
            extra = "   -> VMR: %s" % (b if b else "?")
        print("     %-22s %-42s%s" % (k, (v[:42] if v else ""), extra))


print()
print("=" * 70)
print("=== 样例 A：Family=Partitiviridae 但 Genus=Biavirus ===")
for r in [x for x in rows if x["Family"] == "Partitiviridae" and x["Genus"] == "Biavirus"][:2]:
    dump(r, "Partitiviridae + Biavirus")

print()
print("=" * 70)
print("=== 样例 B：Family=Caulimoviridae 但被判 RNA ===")
n = 0
for r in rows:
    if r["Family"] != "Caulimoviridae":
        continue
    per = {k: ann.resolve_rank(r[k], k, idx) for k in RANKS}
    per = {k: v for k, v in per.items() if v}
    if len(set(per.values())) > 1 and per.get("Family") == "DNA":
        dump(r, "Caulimoviridae 但谱系有 RNA 层")
        n += 1
        if n >= 2:
            break

print()
print("=" * 70)
print("=== 参考库怎么说这两个名字 ===")
for rank, name in [("Family", "Partitiviridae"), ("Genus", "Biavirus"),
                   ("Family", "Caulimoviridae"), ("Genus", "Mimivirus"),
                   ("Genus", "Sirevirus"), ("Genus", "Geminivirus")]:
    gc = idx["lv"][rank].get(name)
    print("  %-8s %-20s -> %s" % (rank, name, dict(gc) if gc else "未收录"))
