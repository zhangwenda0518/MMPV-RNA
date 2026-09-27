#!/usr/bin/env python3
"""打印该 contig 在当前真实产物与校准后产物里的完整行（全部列）"""
import csv
from pathlib import Path

PREFIX = "CRR1126135_clean_NODE_68_length_2879_cov_8"
FILES = [
    ("现行真实产物 (05_Taxonomy/Votus.integrated)",
     Path("/home/zhangwenda/MMPV-paper/goji-virome/02_novel_virus/RNA-Lycium_barbarum_out"
          "/05_Taxonomy/Votus.integrated/final_integrated_classification.tsv")),
    ("校准后 (/tmp/rc_live, 已部署代码跑出)",
     Path("/tmp/rc_live/final_integrated_classification.tsv")),
]

for label, p in FILES:
    if not p.exists():
        print("[%s] 文件不存在: %s\n" % (label, p))
        continue
    with open(p, newline="", encoding="utf-8", errors="replace") as f:
        rd = csv.DictReader(f, delimiter="\t")
        cols = rd.fieldnames
        hit = [r for r in rd if r["contig_id"].startswith(PREFIX)]
    print("=== %s ===" % label)
    print("    文件: %s" % p)
    print("    修改时间: %s   列数: %d" % (
        __import__("datetime").datetime.fromtimestamp(p.stat().st_mtime).strftime("%Y-%m-%d %H:%M"), len(cols)))
    if not hit:
        print("    未命中\n")
        continue
    r = hit[0]
    print("    命中 %d 行; contig_id = %s" % (len(hit), r["contig_id"]))
    for c in cols:
        print("      %-28s = %s" % (c, r.get(c, "")))
    print()

# 顺带看 9a 是否留下了带 Nucleic_acid 的副本
base = Path("/home/zhangwenda/MMPV-paper/goji-virome/02_novel_virus/RNA-Lycium_barbarum_out/05_Taxonomy")
print("=== 05_Taxonomy 下与 final/annotated 相关的文件 ===")
for p in sorted(base.rglob("*classification*.tsv")) + sorted(base.rglob("*annotated*.tsv")):
    print("  %s  (%.0f B, %s)" % (p, p.stat().st_size,
          __import__("datetime").datetime.fromtimestamp(p.stat().st_mtime).strftime("%m-%d %H:%M")))
