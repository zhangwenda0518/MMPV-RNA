import os
import csv

M = "/home/zhangwenda/plant_virus_db/1.virus-host_db/C-host_classify/classified_clean/miss"

targets = ["Mimiviridae", "Metaviridae", "Biavirus", "Partitiviridae", "Schizomimiviridae"]

files = [
    "Plant_L1_kingdom.tsv",
    "Plant_L2b_non_plant_taxon.tsv",
    "Plant_L2_collision.tsv",
    "Plant_L2d_protist.tsv",
    "Plant_L3_phage.tsv",
    "Plant_L4a_blacklist_family.tsv",
    "Plant_L4b_fungal_genus.tsv",
    "Plant_L4c_pathogen_name.tsv",
    "Plant_L4e_env_CRESS.tsv",
    "Plant_all_removed.tsv",
    "Plant_S_conflict.tsv",
]

for fn in files:
    p = os.path.join(M, fn)
    if not os.path.isfile(p):
        continue
    n = 0
    hits = []
    with open(p, errors="replace") as f:
        head = f.readline().rstrip("\n").split("\t")
        print(f"\n=== {fn} ===")
        print(f"  列: {head[:8]}")
        f.seek(0)
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            n += 1
            blob = "\t".join(str(v) for v in row.values())
            for t in targets:
                if t in blob:
                    hits.append((t, row))
        print(f"  行数: {n}")
    if hits:
        print(f"  *** 命中 {len(hits)} 条 ***")
        seen = set()
        for t, row in hits[:6]:
            key = (t,) + tuple(list(row.values())[:3])
            if key in seen:
                continue
            seen.add(key)
            vals = list(row.values())
            print(f"    [{t}] " + " | ".join(str(v)[:45] for v in vals[:6]))
