#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""09 层 DNA/RNA 补列：只读预演（不覆盖任何原有文件）

落点：<proj>_out/09_Virome_Analysis/all_plant_analysis/All_plant.viruses_info.tsv
      <proj>_out/09_Virome_Analysis/HQ_analysis/HQ_plant_viruses_info.tsv
做法：追加 3 列（不动任何现有列）
      Genome_Type      DNA / RNA
      Genome_Structure ds / ss / ds-RT / ss-RT（VMR Genome 原文拆解）
      Genome_Source    species / genus / tax / none（取值来源，便于追溯）
取值优先级：VMR Species 精确匹配 > VMR Genus 唯一取值兜底 > Realm/Kingdom 分类规则 > none

预演产物写 /tmp/dnarna_preview/ 下，不改 09 层。
"""
import csv, os, glob, json
from collections import Counter, defaultdict

VMR = os.path.expanduser("~/plant_virus_db/1.virus-host_db/B-ictv/VMR_MSL41.tsv")
ROOT = os.path.expanduser("~/MMPV-paper")
OUTDIR = "/tmp/dnarna_preview"

# ---------- 1. VMR 字典 ----------
sp2genome, genus2genomes = {}, defaultdict(set)
with open(VMR, newline="", encoding="utf-8", errors="replace") as f:
    rd = csv.reader(f, delimiter="\t")
    next(rd)
    for r in rd:
        if len(r) < 27:
            continue
        g = r[25].strip().strip('"')      # Genome
        genus = r[15].strip().strip('"')  # Genus
        sp = r[17].strip().strip('"')     # Species
        if sp and g:
            sp2genome.setdefault(sp, g)
        if genus and g:
            genus2genomes[genus].add(g)

def split_genome(g):
    """ssRNA(+) -> (RNA, ss); dsDNA-RT -> (DNA, ds-RT)"""
    if not g:
        return "", ""
    t = "DNA" if "DNA" in g.upper() else ("RNA" if "RNA" in g.upper() else "")
    struc = []
    gl = g.upper()
    if gl.startswith("DS"):
        struc.append("ds")
    elif gl.startswith("SS"):
        struc.append("ss")
    if "RT" in gl:
        struc.append("RT")
    if "+" in g:
        struc.append("(+)")
    elif "-" in g:
        struc.append("(-)")
    return t, "-".join(struc)

# Realm / Kingdom 分类规则（兜底）
REALM_RULE = {
    "Riboviria": "RNA",
    "Ribozyviria": "RNA",
    "Monodnaviria": "DNA",
    "Varidnaviria": "DNA",
    "Duplodnaviria": "DNA",
    "Adnaviria": "DNA",
}
KINGDOM_RULE = {
    "Orthornavirae": "RNA",
    "Pararnavirae": "DNA",
    "Loebvirae": "DNA",
    "Sangervirae": "DNA",
    "Shotokuvirae": "DNA",
    "Trapavirae": "DNA",
    "Bamfordvirae": "DNA",
    "Helvetiavirae": "DNA",
    "Riboviria": "RNA",
}

def tax_call(row, idx):
    realm = row[idx["Realm"]].strip().strip('"') if "Realm" in idx else ""
    king = row[idx["Kingdom"]].strip().strip('"') if "Kingdom" in idx else ""
    if realm in REALM_RULE:
        return REALM_RULE[realm]
    if king in KINGDOM_RULE:
        return KINGDOM_RULE[king]
    return ""

# ---------- 2. 目标文件 ----------
TARGETS = [
    "09_Virome_Analysis/all_plant_analysis/All_plant.viruses_info.tsv",
    "09_Virome_Analysis/HQ_analysis/HQ_plant_viruses_info.tsv",
]
files = []
for rel in TARGETS:
    files += sorted(glob.glob(os.path.join(ROOT, "*", "02_novel_virus", "*_out", rel)))
    files += sorted(glob.glob(os.path.join(ROOT, "*", "onekp-virus", rel)))

print("目标文件 %d 个\n" % len(files))
grand = Counter()
miss_species = Counter()
conflict_pairs = Counter()
conflict_examples = []
hdr_out = ["Genome_Type", "Genome_Structure", "Genome_Source"]

for p in files:
    with open(p, newline="", encoding="utf-8", errors="replace") as f:
        rd = csv.reader(f, delimiter="\t")
        hdr = [h.strip().strip('"') for h in next(rd)]
        idx = {h: i for i, h in enumerate(hdr)}
        rows = list(rd)
    stat = Counter()
    out_rows = []
    for r in rows:
        sp = r[idx["Species"]].strip().strip('"') if "Species" in idx else ""
        ge = r[idx["Genus"]].strip().strip('"') if "Genus" in idx else ""
        gtype, gstruc, src = "", "", "none"
        if sp and sp in sp2genome:
            gtype, gstruc = split_genome(sp2genome[sp]); src = "species"
        elif ge and len(genus2genomes.get(ge, ())) == 1:
            gtype, gstruc = split_genome(next(iter(genus2genomes[ge]))); src = "genus"
        elif ge and len(genus2genomes.get(ge, ())) > 1:
            src = "genus_ambiguous"
        if not gtype:
            t = tax_call(r, idx)
            if t:
                gtype, src = t, "tax" if src == "none" else src + "+tax"
        tax_t = tax_call(r, idx)
        if gtype and tax_t and gtype != tax_t:
            stat["conflict_tax_vs_value"] += 1
            conflict_pairs[(gtype, tax_t)] += 1
            if len(conflict_examples) < 12:
                conflict_examples.append((p.split("MMPV-paper/")[-1], sp or ge, gtype, tax_t,
                                          r[idx["Family"]].strip().strip('"') if "Family" in idx else ""))
        stat[src] += 1
        if src == "none" or src == "genus_ambiguous":
            miss_species[sp or ge or "(空)"] += 1
        out_rows.append(r + [gtype, gstruc, src])

    n = len(rows)
    hit = stat["species"] + stat["genus"] + stat["tax"] + stat["species+tax"] + stat["genus+tax"] + stat["tax+tax"]
    print("%s" % p.split("MMPV-paper/")[-1])
    print("   行=%d | species=%d genus=%d tax兜底=%d 属级有歧义=%d 完全未命中=%d | 与分类冲突=%d" % (
        n, stat["species"], stat["genus"], stat["tax"] + stat["species+tax"] + stat["genus+tax"],
        stat["genus_ambiguous"], stat["none"], stat["conflict_tax_vs_value"]))
    print("   DNA/RNA 分布: %s" % dict(Counter(x[-3] for x in out_rows if x[-3])))
    grand.update(stat)

    dst = os.path.join(OUTDIR, p.split("MMPV-paper/")[-1].replace("/", "__"))
    os.makedirs(OUTDIR, exist_ok=True)
    with open(dst, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter="\t", quoting=csv.QUOTE_MINIMAL)
        w.writerow(hdr + hdr_out)
        w.writerows(out_rows)

print("\n===== 合计 =====")
tot = sum(v for k, v in grand.items() if k not in ("conflict_tax_vs_value",))
for k in ("species", "genus", "tax", "species+tax", "genus+tax", "tax+tax", "genus_ambiguous", "none"):
    if grand[k]:
        print("  %-16s %d (%.2f%%)" % (k, grand[k], 100.0 * grand[k] / tot))
print("  冲突(tax vs VMR) %d  %s" % (grand["conflict_tax_vs_value"], dict(conflict_pairs)))

print("\n===== 未完全命中的物种/属名（top 40，用于判断是否命名口径差异） =====")
for k, v in miss_species.most_common(40):
    print("  %-56s x%d" % (k[:56], v))
print("\n===== 冲突样例 =====")
for c in conflict_examples:
    print("  %s | %s | 取值判=%s 分类判=%s | Family=%s" % c)
