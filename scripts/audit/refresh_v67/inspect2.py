#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Read-only check of rankedlineage.dmp column semantics under the assumed mapping."""
import sys

COLS = ["tax_id", "tax_name", "species", "genus", "family", "order",
        "class", "phylum", "kingdom", "superkingdom"]
F = "/home/zhangwenda/database/taxonomy/rankedlineage.dmp"
targets = {"homo sapiens", "narnaviridae", "escherichia coli", "root",
           "riboviria", "orthornavirae", "wolframvirales", "sars-cov-2"}

n = 0
lens = {}
with open(F, encoding="utf-8", errors="ignore") as f:
    for i, line in enumerate(f, 1):
        if i <= 5:
            lens[len(line.split("\t|\t"))] = lens.get(len(line.split("\t|\t")), 0) + 1
        parts = [p.strip().strip("|").strip() for p in line.split("\t|\t")]
        if len(parts) < 10:
            continue
        rec = dict(zip(COLS, parts))
        nm = rec["tax_name"].lower()
        if nm in targets:
            print("%-16s id=%-10s sp=[%s] ge=[%s] fa=[%s] or=[%s] cl=[%s] ph=[%s] ki=[%s] sk=[%s]" % (
                rec["tax_name"], rec["tax_id"], rec["species"], rec["genus"], rec["family"],
                rec["order"], rec["class"], rec["phylum"], rec["kingdom"], rec["superkingdom"]))
            targets.discard(nm)
        if not targets or i > 3000000:
            break
print("parts_len_first5:", lens)
# count how many rows have species column nonempty over full file (approx via streaming)
cnt = {"species": 0, "family": 0, "superkingdom": 0}
with open(F, encoding="utf-8", errors="ignore") as f:
    for line in f:
        parts = [p.strip().strip("|").strip() for p in line.split("\t|\t")]
        if len(parts) < 10:
            continue
        rec = dict(zip(COLS, parts))
        for k in cnt:
            if rec[k]:
                cnt[k] += 1
print("nonempty counts:", cnt)
