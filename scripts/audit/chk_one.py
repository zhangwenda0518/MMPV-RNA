#!/usr/bin/env python3
"""单 contig 状态追踪: 各版本终表行 + 七工具原始谱系 + standardized 行"""
import csv
from pathlib import Path

PREFIX = "CRR1126135_clean_NODE_68_length_2879_cov_8"
PIPE = Path("/home/zhangwenda/MMPV-RNA/virome_discovery_pipeline")
RAW = Path("/home/zhangwenda/MMPV-paper/goji-virome/02_novel_virus/RNA-Lycium_barbarum_out"
           "/05_Taxonomy/Votus.classed/Votus_combined_taxonomy.tsv")
FNAME = "final_integrated_classification.tsv"
TABLES = [
    ("Aug12表(旧)", Path("/home/zhangwenda/MMPV-paper/goji-virome/02_novel_virus/RNA-Lycium_barbarum_out"
                         "/05_Taxonomy/Votus.integrated") / FNAME),
    ("base(现行)", Path("/tmp/rc_base") / FNAME),
    ("R1-last", Path("/tmp/rc_patch") / FNAME),
    ("R2-last", Path("/tmp/rc_v_last_noblank") / FNAME),
    ("R1-first", Path("/tmp/rc_v_first_blank") / FNAME),
    ("R2-first", Path("/tmp/rc_v_first_noblank") / FNAME),
]
TOOLS = ["ACVirus", "CAT", "diamond_lca", "genomad", "metabuli", "mmseqs", "VITAP"]
RANKS = ["Realm", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]


def rd(p):
    with open(p, newline="", encoding="utf-8", errors="replace") as f:
        return list(csv.DictReader(f, delimiter="\t"))


print("=== 1. 各版本终表中的该 contig ===")
for name, p in TABLES:
    if not p.exists():
        print("\n[%s] 缺失" % name)
        continue
    hit = [r for r in rd(p) if r["contig_id"].startswith(PREFIX)]
    if not hit:
        print("\n[%s] 未命中" % name)
        continue
    r = hit[0]
    print("\n[%s]  contig_id=%s" % (name, r["contig_id"]))
    print("  primary_tool=%s  completeness=%s  confidence=%s"
          % (r.get("primary_tool"), r.get("completeness"), r.get("confidence")))
    for k in RANKS:
        print("    %-8s %-28s | %s" % (k, r.get(k, ""), r.get(k + "_agree", "")))

print("\n\n=== 2. 原始 combined 表里七个工具的谱系 ===")
rows = [r for r in rd(RAW) if r["seq_name"].startswith(PREFIX)]
print("命中 %d 行" % len(rows))
for r in rows:
    lin = " > ".join((r.get(k) or "NA") for k in RANKS)
    print("  [%-11s] %s" % (r.get("tool"), lin))

print("\n\n=== 3. standardized_*.tsv 里的行 ===")
for t in TOOLS:
    p = Path("/tmp/rc_base") / ("standardized_%s.tsv" % t)
    if not p.exists():
        continue
    hit = [r for r in rd(p) if r["contig_id"].startswith(PREFIX)]
    for r in hit:
        lin = " > ".join((r.get(k) or "NA") for k in RANKS)
        print("  [%-11s] %s" % (t, lin))
