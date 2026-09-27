"""两件事: 1) taxa.txt 是否全量 MSL41 子集; 2) Ourmiaviridae 那 3 行为何没被闸门置空."""
import csv
import os
from collections import Counter

MSL = "/home/zhangwenda/database/virus-db/acvirus_db/taxa.txt"
REF = os.path.expanduser("~/database/taxonomy/genus_family_ref.tsv")
BASE = ("/home/zhangwenda/MMPV-paper/goji-virome/02_novel_virus/RNA-Lycium_barbarum_out"
        "/05_Taxonomy/Votus.integrated")
CAL = f"{BASE}/calibration_20260914/final_integrated_classification.calibrated.tsv"

# 1. taxa.txt 的界/纲分布
realm = Counter()
cls = Counter()
has_auto = 0
fams = set()
with open(MSL, errors="replace") as f:
    for r in csv.DictReader(f):
        realm[(r.get("Realm") or "").strip()] += 1
        cls[(r.get("Class") or "").strip()] += 1
        fam = (r.get("Family") or "").strip()
        fams.add(fam.lower())
print("taxa.txt Realm 分布:", dict(realm.most_common(8)))
print("\ntaxa.txt 是否含噬菌体纲 Caudoviricetes:", "Caudoviricetes" in cls)
print("含 Autographiviridae:", "autographiviridae" in fams)
print("含 Hepaciviridae:", "hepaciviridae" in fams, "| 含 Pestiviridae:", "pestiviridae" in fams,
      "| 含 Ourmiaviridae:", "ourmiaviridae" in fams, "| 含 Botourmiaviridae:", "botourmiaviridae" in fams,
      "| 含 Parahypoviridae:", "parahypoviridae" in fams, "| 含 Zimmerviridae:", "zimmerviridae" in fams,
      "| 含 Ambiguiviridae:", "ambiguiviridae" in fams)
print("Class 前 8:", dict(cls.most_common(8)))

# 2. 参照表里 Ourmiavirus / Iotaourmiavirus 的定型科
print("\n=== 参照表里这两个属 ===")
with open(REF, errors="replace") as f:
    for r in csv.DictReader(f, delimiter="\t"):
        if (r.get("Genus") or "").strip().lower() in ("ourmiavirus", "iotaourmiavirus", "alphambiguivirus", "orthohepacivirus", "jouyvirus"):
            print("  ", dict(r))

# 3. Ourmiaviridae 的 3 行全列
print("\n=== 校准表里 Family=Ourmiaviridae 的行 ===")
with open(CAL, errors="replace") as f:
    rd = csv.DictReader(f, delimiter="\t")
    cols = rd.fieldnames
    calib_cols = [c for c in cols if c.startswith("calib")]
    for r in rd:
        if (r.get("Family") or "").strip() == "Ourmiaviridae":
            print("  contig:", r.get("contig_id"))
            print("   ", {k: r.get(k) for k in ("Realm", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "confidence", "primary_tool")})
            print("   ", {k: r.get(k) for k in calib_cols})
