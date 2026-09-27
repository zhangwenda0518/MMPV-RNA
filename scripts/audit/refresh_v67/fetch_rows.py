#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Print full rows (quoted TSV) for a few contig_ids, using csv so quotes are stripped."""
import csv
import sys

jobs = [
    ("onekp NEW",
     "/home/zhangwenda/MMPV-paper/onekp-virome/onekp-virus/05_Taxonomy/Votus.integrated/final_integrated_classification.tsv",
     {"ERR2040124_clean_NODE_2483_length_2313_cov_14.886344",
      "ERR2040122_clean_NODE_23705_length_598_cov_3.039640"}),
    ("barbarum NEW",
     "/home/zhangwenda/MMPV-paper/goji-virome/02_novel_virus/RNA-Lycium_barbarum_out/05_Taxonomy/Votus.integrated/final_integrated_classification.tsv",
     {"CRR1440122_clean_NODE_16_length_2430_cov_3.085702",
      "CRR732652_clean_NODE_613_length_1528_cov_1.500343"}),
    ("Fusarium NEW",
     "/home/zhangwenda/MMPV-paper/goji-virome/02_novel_virus/RNA-Fusarium_nematophilum_out/05_Taxonomy/Votus.integrated/final_integrated_classification.tsv",
     {"SRR29516379_clean_NODE_7132_length_1749_cov_26.553699", "contig_869"}),
    ("ruthenicum NEW",
     "/home/zhangwenda/MMPV-paper/goji-virome/02_novel_virus/RNA-Lycium_ruthenicum_out/05_Taxonomy/Votus.integrated/final_integrated_classification.tsv",
     {"CRR527056_clean_NODE_4696_length_1872_cov_28.481934",
      "CRR126177_clean_NODE_3084_length_672_cov_1.757930"}),
]

for label, path, ids in jobs:
    print("### " + label)
    with open(path, encoding="utf-8", errors="ignore", newline="") as f:
        rd = csv.DictReader(f, delimiter="\t", quotechar='"')
        cols = rd.fieldnames
        for row in rd:
            if row.get("contig_id") in ids:
                print("  " + row["contig_id"])
                print("    " + " | ".join("%s=%s" % (c, row.get(c)) for c in cols))
    print()
