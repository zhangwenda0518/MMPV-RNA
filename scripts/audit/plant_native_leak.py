"""原生病毒为什么躺在植物表里: 级别优先规则的漏网路径 + 全库宿主分离能力盘点。

只读。
"""
import csv
import os
from collections import Counter

CA = "/home/zhangwenda/MMPV-RNA/database/cross_analysis"
RUN = "/home/zhangwenda/MMPV-paper/goji-virome/02_novel_virus/RNA-Lycium_barbarum_out"
PLANT = f"{RUN}/06_HostPrediction/C9_ICTV_result/Plant.classified.tsv"
FULL = f"{RUN}/05_Taxonomy/Votus.integrated/final_integrated_classification.tsv"
FULLHOST = f"{RUN}/06_HostPrediction/C9_ICTV_result/classification_result.tsv"

CATS = ["Algae", "Animal_other", "Arachnida", "Archaea", "Aves", "Bacteria",
        "Fungi", "Human", "Insecta", "Mammalia", "Oomycetes", "Plant", "Protist"]
NATIVE = {"Protist", "Algae"}


def rows_of(path):
    with open(path, errors="replace") as f:
        for d in csv.DictReader(f, delimiter="\t"):
            yield {(k.strip() if k else k): (v.strip().strip('"') if isinstance(v, str) else v)
                   for k, v in d.items()}


def load_prob(level):
    out = {}
    p = f"{CA}/{level.lower()}_host_probability.tsv"
    for r in rows_of(p):
        n = r.get(level, "")
        if n:
            out[n] = r
    return out


def fl(r, k, d=0.0):
    try:
        return float(r.get(k) or d)
    except Exception:
        return d


fam_p, gen_p, sp_p = load_prob("Family"), load_prob("Genus"), load_prob("Species")
nat_gen = {k for k, v in gen_p.items() if v["Predicted_Host"] in NATIVE}
nat_fam = {k for k, v in fam_p.items() if v["Predicted_Host"] in NATIVE}

print("=" * 96)
print("[A] 植物表里带着原生科/属名字的行, 它们靠什么进的植物表")
print("=" * 96)
hdr = f"{'contig_id':<42} {'Family':<16} {'Genus':<18} {'Species':<26} {'det':<20} {'conf':<18} {'IC':>6}"
print(hdr)
sus = []
for r in rows_of(PLANT):
    if r.get("Family") in nat_fam or r.get("Genus") in nat_gen:
        sus.append(r)
for r in sorted(sus, key=lambda x: (x.get("Family", ""), x.get("Genus", "")))[:40]:
    print(f"{r.get('contig_id','')[:42]:<42} {r.get('Family','')[:16]:<16} {r.get('Genus','')[:18]:<18} "
          f"{r.get('Species','')[:26]:<26} {r.get('Determination_Level','')[:20]:<20} "
          f"{r.get('Confidence_Level','')[:18]:<18} {fl(r,'Integrated_Confidence'):>6.3f}")
print(f"\n  合计 {len(sus)} 行")
print("  Determination_Level 构成:", dict(Counter(r.get("Determination_Level", "") for r in sus)))
print("  这些行命中的属级条目(原生):", dict(Counter(r.get("Genus") for r in sus if r.get("Genus") in nat_gen)))
print("  这些行命中的科级条目(原生):", dict(Counter(r.get("Family") for r in sus if r.get("Family") in nat_fam)))

# 这些行的 Species 级条目证据有多弱
print("\n  这些行的 Species 级条目证据(用于解释为何被压过):")
for r in sus:
    sp = r.get("Species", "")
    if sp in sp_p:
        v = sp_p[sp]
        print(f"    {sp[:40]:<40} species_host={v['Predicted_Host']:<12} records={v['Total_Records']:>5} "
              f"lvl={v['Confidence_Level']:<20} IC={v['Integrated_Confidence']}")
    else:
        print(f"    {sp[:40]:<40} (无 Species 级条目)  该行判定={r.get('Determination_Level')}")

print("\n" + "=" * 96)
print("[B] 关键科在概率表中的状态")
print("=" * 96)
for name, tb, tag in [("Schizomimiviridae", fam_p, "Family"), ("Mimiviridae", fam_p, "Family"),
                      ("Marseilleviridae", fam_p, "Family"), ("Phycodnaviridae", fam_p, "Family"),
                      ("Partitiviridae", fam_p, "Family"), ("Biavirus", gen_p, "Genus"),
                      ("Pandoravirus", gen_p, "Genus"), ("Chlorovirus", gen_p, "Genus")]:
    v = tb.get(name)
    if v:
        print(f"  [{tag}] {name:<20} host={v['Predicted_Host']:<12} P_Max={v['P_Max']:<8} "
              f"records={v['Total_Records']:>6} Shannon={v['Shannon_Entropy']:<8} "
              f"P(Plant)={v.get('P(Plant)','-'):<8} P(Protist)={v.get('P(Protist)','-')}")
    else:
        print(f"  [{tag}] {name:<20} 表中缺失")

print("\n" + "=" * 96)
print("[C] 全库宿主判定能力: 按 Predicted_Host 与 Determination_Level 交叉")
print("=" * 96)
src = FULLHOST if os.path.exists(FULLHOST) else PLANT
print(f"  源: {src}")
cross = Counter()
for r in rows_of(src):
    cross[(r.get("Predicted_Host", "?"), r.get("Determination_Level", "?"))] += 1
tot = sum(cross.values())
print(f"  总行数 {tot}")
print(f"  {'Predicted_Host':<14} {'Species':>8} {'Genus':>8} {'Family':>8} {'Order':>8} {'None':>8} {'合计':>8}")
for cat in CATS + ["Unknown"]:
    sl = {lvl: sum(v for (c, lvl), v in cross.items() if c == cat and lvl.startswith(prefix))
          for prefix, lvl in [("Species", "Species"), ("Genus", "Genus"), ("Family", "Family"), ("Order", "Order")]}
    nl = sum(v for (c, lvl), v in cross.items() if c == cat and not lvl.startswith(("Species", "Genus", "Family", "Order")))
    row_tot = sum(v for (c, l), v in cross.items() if c == cat)
    if row_tot == 0:
        continue
    print(f"  {cat:<14} {sl['Species']:>8} {sl['Genus']:>8} {sl['Family']:>8} {sl['Order']:>8} {nl:>8} {row_tot:>8}")

print("\n" + "=" * 96)
print("[D] 全库中判为原生/藻类的行, 其科属构成")
print("=" * 96)
nat_rows = [r for r in rows_of(src) if r.get("Predicted_Host") in NATIVE]
print(f"  行数 {len(nat_rows)}")
print("  科分布 top15:", dict(Counter(r.get("Family", "") for r in nat_rows).most_common(15)))
print("  属分布 top15:", dict(Counter(r.get("Genus", "") for r in nat_rows).most_common(15)))
print("  Determination_Level:", dict(Counter(r.get("Determination_Level", "") for r in nat_rows)))
