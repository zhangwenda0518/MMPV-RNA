import csv
import collections

F = "/home/zhangwenda/plant_virus_db/1.virus-host_db/C-host_classify/host_extract/Final_Virus_Host_Lineage.tsv"

fams = ["Mimiviridae", "Metaviridae", "Biavirus"]
hits = {k: [] for k in fams}

with open(F, errors="replace") as f:
    r = csv.DictReader(f, delimiter="\t")
    for row in r:
        lin = row.get("Virus_lineage", "") or ""
        for k in fams:
            if k in lin or (k == "Biavirus" and "Biavirus" in lin):
                if len(hits[k]) < 12:
                    hits[k].append(row)
                break

for k in fams:
    print(f"\n{'='*70}\n=== {k}: 样例 ===")
    for row in hits[k][:8]:
        lin = row.get("Virus_lineage", "")
        parts = [p for p in lin.split(";") if p]
        fam = next((p for p in parts if p.endswith("viridae")), "?")
        gen = next((p for p in parts if p.endswith("virus") and p != fam), "?")
        print(f"  {row.get('Accession','')[:12]:12} | {row.get('Virus_Name','')[:38]:38} | fam={fam[:22]:22} | host={row.get('Host_Name','')[:25]:25} | {row.get('Host_lineage','')[:50]}")
