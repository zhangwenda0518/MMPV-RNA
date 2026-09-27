"""同步更新派生 TSV: 从各元数据表中剔除黑名单 contig。

覆盖:
  - 09_Virome_Analysis/HQ_analysis/HQ_plant_viruses_info.tsv  (9415 -> 8916 行含表头)
  - 10_Reports/HQ_plant_viruses_info.tsv
  - 10_Reports/All_plant.viruses_info.tsv
  - 10_Reports/plant_virus_summary.tsv
  - 10_Reports/rescue_evidence_scored.tsv
  - 10_Reports/cdd_evidence_report.tsv
  - 10_Reports/plant_virus_cluster_summary.tsv
"""
import os

OUT = "/home/zhangwenda/data-test/out"
kill = set()
with open("/tmp/blacklist_kill_ids.txt") as f:
    for l in f:
        l = l.strip()
        if l:
            kill.add(l)
print(f"剔除集: {len(kill)}\n")

# (路径, 首列是否 contig_id)
TARGETS = [
    (f"{OUT}/09_Virome_Analysis/HQ_analysis/HQ_plant_viruses_info.tsv", 0),
    (f"{OUT}/10_Reports/HQ_plant_viruses_info.tsv", 0),
    (f"{OUT}/10_Reports/All_plant.viruses_info.tsv", 0),
    (f"{OUT}/10_Reports/plant_virus_summary.tsv", 0),
    (f"{OUT}/10_Reports/rescue_evidence_scored.tsv", 0),
    (f"{OUT}/10_Reports/cdd_evidence_report.tsv", 0),
    (f"{OUT}/10_Reports/plant_virus_cluster_summary.tsv", 0),
    (f"{OUT}/10_Reports/All_plant.viruses_info.tsv", 0),
]

print(f"{'文件':<62} {'前':>7} {'后':>7} {'删':>6}")
print("-" * 88)
for path, col in TARGETS:
    if not os.path.isfile(path):
        print(f"{path.replace(OUT+'/',''):<62} 缺失")
        continue
    tmp = path + ".tmp_kill"
    nb = na = 0
    with open(path, errors="replace") as fin, open(tmp, "w") as fout:
        for i, line in enumerate(fin):
            if i == 0:
                fout.write(line)
                continue
            parts = line.split("\t")
            cid = parts[col].strip().strip('"') if len(parts) > col else ""
            nb += 1
            if cid in kill:
                continue
            na += 1
            fout.write(line)
    os.replace(tmp, path)
    rel = path.replace(OUT + "/", "")
    print(f"{rel:<62} {nb:>7} {na:>7} {nb-na:>6}")
