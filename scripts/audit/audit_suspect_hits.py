"""交叉: 20 个可疑属在当前 Plant.classified.fasta (20054 条) 里的实际命中量。

只有真正出现在植物结果里的属才有清理意义。
"""
import csv
from collections import Counter

OUT = "/home/zhangwenda/data-test/out"
PLANT_FA = f"{OUT}/06_HostPrediction/host_classified_fasta/Plant.classified.fasta"
CLS = f"{OUT}/06_HostPrediction/C9_ICTV_result/classification_result.tsv"

SUSPECT = {
    "Partitiviridae": ["Betapartitivirus", "Cryspovirus", "Gammapartitivirus"],
    "Alphaflexiviridae": ["Botrexvirus", "Sclerodarnavirus"],
    "Alphasatellitidae": ["Draflysatellite"],
    "Amalgaviridae": ["Unirnavirus", "Zybavirus"],
    "Endornaviridae": ["Betaendornavirus"],
    "Geplanaviridae": ["Acciovirus", "Aguamentivirus", "Aparevirus", "Avadavirus",
                       "Densaugvirus", "Engorgivirus", "Expellivirus", "Impedivirus",
                       "Oblivivirus", "Protegovirus", "Wingardivirus"],
}

# contig -> (Family, Genus)
info = {}
with open(CLS, errors="replace") as f:
    for row in csv.DictReader(f, delimiter="\t"):
        info[row["contig_id"]] = ((row.get("Family", "") or "").strip(),
                                  (row.get("Genus", "") or "").strip())

# 遍历 Plant.classified.fasta
hits = Counter()
fam_hits = Counter()
total = 0
with open(PLANT_FA) as f:
    for line in f:
        if not line.startswith(">"):
            continue
        total += 1
        cid = line[1:].split()[0]
        fm, gen = info.get(cid, ("", ""))
        fam_hits[fm] += 1
        for fam, gens in SUSPECT.items():
            if gen in gens:
                hits[(fam, gen)] += 1
                break

print(f"Plant.classified.fasta 当前: {total} 条\n")
print("=== 20 个可疑属的实际命中 ===")
print(f"{'科':<22} {'属':<22} {'命中':>6}")
print("-" * 54)
tot = 0
for (fam, gen), n in hits.most_common():
    print(f"{fam:<22} {gen:<22} {n:>6}")
    tot += n
print("-" * 54)
print(f"{'合计':<46} {tot:>6}")

print("\n=== 这 20 属所属科的总体命中 (对照) ===")
for fam in SUSPECT:
    print(f"  {fam:<24} {fam_hits.get(fam,0):>6}")
