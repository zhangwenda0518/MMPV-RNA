#!/usr/bin/env python3
"""Biavirus 是什么、为什么进了结果

1. VMR_MSL41 里 Biavirus 的全部字段（科/类/基因组/物种）
2. 种名里含 Guapo 的所有 VMR 记录（看 Guapo partitivirus 归在哪）
3. Schizomimiviridae 是个什么科（界/门/纲/目/基因组/属清单）
4. 那批 Genus=Biavirus 的植物病毒 contig，逐工具看谁投了 Biavirus、谁投了 Partitiviridae
"""
import csv
from collections import Counter, defaultdict
from pathlib import Path

DT = Path("/home/zhangwenda/MMPV-paper/goji-virome/02_novel_virus/RNA-Lycium_barbarum_out")
VMR = Path("/home/zhangwenda/plant_virus_db/1.virus-host_db/B-ictv/VMR_MSL41.tsv")
STD = DT / "05_Taxonomy/Votus.integrated"
TOOLS = ["ACVirus", "CAT", "VITAP", "diamond_lca", "genomad", "metabuli", "mmseqs"]


def clean(v):
    return (v or "").strip().strip('"')


def load(p, **kw):
    with open(p, newline="", encoding="utf-8", errors="replace") as f:
        return list(csv.DictReader(f, delimiter="\t", **kw))


vmr = load(VMR)
cols = list(vmr[0].keys())
print("VMR 列 (%d): %s\n" % (len(cols), " | ".join(cols)))

gen_col = next((c for c in cols if c.lower() == "genome"), None)
print("=" * 90)
print("【1】VMR 里 Genus = Biavirus 的记录")
print("=" * 90)
bia = [r for r in vmr if clean(r.get("Genus")).lower() == "biavirus"]
print("记录数: %d" % len(bia))
for k in ["Realm", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", gen_col, "Species", "ICTV_ID"]:
    if k:
        vals = Counter(clean(r.get(k)) or "(空)" for r in bia)
        print("  %-16s %s" % (k, dict(vals)))
for r in bia[:5]:
    print("    样本行: %s" % {k: clean(r.get(k)) for k in
                           ("Realm", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", gen_col, "Species")})

print("\n" + "=" * 90)
print("【2】种名含 Guapo 的 VMR 记录")
print("=" * 90)
gu = [r for r in vmr if "guapo" in clean(r.get("Species")).lower()]
print("记录数: %d" % len(gu))
for r in gu:
    print("  %s | Family=%s Genus=%s Species=%s Genome=%s" %
          (clean(r.get("Realm")), clean(r.get("Family")), clean(r.get("Genus")),
           clean(r.get("Species")), clean(r.get(gen_col)) if gen_col else "-"))

print("\n" + "=" * 90)
print("【3】Schizomimiviridae 是什么")
print("=" * 90)
sch = [r for r in vmr if clean(r.get("Family")) == "Schizomimiviridae"]
print("记录数: %d" % len(sch))
for k in ["Realm", "Kingdom", "Phylum", "Subphylum", "Class", "Order", "Genus", gen_col]:
    if k and k in cols:
        print("  %-16s %s" % (k, dict(Counter(clean(r.get(k)) or "(空)" for r in sch))))
print("  该科属清单: %s" % ", ".join(sorted({clean(r.get("Genus")) for r in sch if clean(r.get("Genus"))})))

print("\n" + "=" * 90)
print("【4】植物病毒子集里 Genus=Biavirus 的 contig，逐工具投票")
print("=" * 90)
ids = set()
with open(DT / "09_Virome_Analysis/All_plant.viruses.fasta", encoding="utf-8", errors="replace") as f:
    for ln in f:
        if ln.startswith(">"):
            ids.add(ln[1:].split()[0].strip())
base = {r["contig_id"]: r for r in load(DT / "05_Taxonomy/Votus.integrated/final_integrated_classification.tsv")}
targets = [c for c, r in base.items() if c in ids and clean(r.get("Genus")) == "Biavirus"]
print("命中 %d 个 contig\n" % len(targets))

votes = {t: {r.get("contig_id"): r for r in load(STD / ("standardized_%s.tsv" % t))} for t in TOOLS}
for cid in sorted(targets):
    r = base[cid]
    print("-" * 90)
    print("%s" % cid)
    print("  最终: Family=%s Genus=%s Species=%s  primary_tool=%s" %
          (clean(r.get("Family")), clean(r.get("Genus")), clean(r.get("Species")), clean(r.get("primary_tool"))))
    print("  %-12s %-22s %-22s %s" % ("工具", "Family", "Genus", "Species"))
    gv = Counter()
    for t in TOOLS:
        row = votes[t].get(cid)
        if row is None:
            print("  %-12s (该 contig 无此工具结果)" % t)
            continue
        f, g, s = clean(row.get("Family")), clean(row.get("Genus")), clean(row.get("Species"))
        print("  %-12s %-22s %-22s %s" % (t, f or "NA", g or "NA", s or "NA"))
        if g:
            gv[g] += 1
    print("  Genus 票: %s" % dict(gv))
