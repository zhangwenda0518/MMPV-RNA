#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""prep_eval_sample.py — 把 step3_identification_eval 评测集整理成管线样本结构"""
import os

D = "/home/zhangwenda/snap/simulator-data/step3_identification_eval"
OUT = "/home/zhangwenda/gp_eval"
os.makedirs(OUT + "/asm/eval", exist_ok=True)
os.makedirs(OUT + "/reads", exist_ok=True)

labels = {}
with open(D + "/sequence_labels_category.tsv") as fh:
    fh.readline()
    for line in fh:
        f = line.rstrip("\n").split("\t")
        if len(f) >= 3:
            labels[f[0]] = (f[1], f[3] if len(f) > 3 else "")
with open(D + "/sequence_labels.tsv") as fh:
    fh.readline()
    for line in fh:
        f = line.rstrip("\n").split("\t")
        if len(f) >= 2 and f[0] not in labels:
            labels[f[0]] = (f[1], "")

n, unmatched = 0, 0
rows = ["eval_id\torig_id\tlabel\tcategory\tlength"]
with open(D + "/evaluation_sequences.fasta") as fin, \
     open(OUT + "/asm/eval/eval_megahit.contig.fasta", "w") as fout:
    w = False
    for line in fin:
        if line.startswith(">"):
            oid = line[1:].rstrip("\n").split()[0]
            n += 1
            eid = "eval_%05d" % n
            lbl, cat = labels.get(oid, ("unmatched", ""))
            if lbl == "unmatched":
                unmatched += 1
            rows.append("%s\t%s\t%s\t%s" % (eid, oid, lbl, cat))
            fout.write(">%s\n" % eid)
            w = True
        else:
            seq = line.strip()
            rows[-1] += "\t%d" % len(seq)
            fout.write(seq + "\n")
            w = False

with open(OUT + "/id_map.tsv", "w") as fh:
    fh.write("\n".join(rows) + "\n")

import collections
cnt = collections.Counter(r.split("\t")[2] for r in rows[1:])
print("sequences:", n, "unmatched:", unmatched)
print("label 分布:", dict(cnt))
