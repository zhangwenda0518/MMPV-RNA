#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Find a few NEW-product onekp rows for a given (Genus, Species) pair (quoted TSV)."""
import csv

path = "/home/zhangwenda/MMPV-paper/onekp-virome/onekp-virus/05_Taxonomy/Votus.integrated/final_integrated_classification.tsv"
want = [("Varicellovirus", "Pandoravirus salinus"),
        ("Medusavirus", "Pandoravirus inopinatum"),
        ("Alphaorpheovirus", "Pandoravirus quercus")]
got = {k: 0 for k in want}
with open(path, encoding="utf-8", errors="ignore", newline="") as f:
    rd = csv.DictReader(f, delimiter="\t", quotechar='"')
    for row in rd:
        g, s = row.get("Genus"), row.get("Species")
        for k in want:
            if got[k] < 1 and g == k[0] and s == k[1]:
                print("%s | Family=%s Genus=%s Species=%s" % (row["contig_id"], row.get("Family"), g, s))
                got[k] += 1
        if all(v >= 1 for v in got.values()):
            break
