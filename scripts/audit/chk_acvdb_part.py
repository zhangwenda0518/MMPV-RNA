#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
chk_acvdb_part.py — ACVirus 库里到底有没有 Partitiviridae，以及它们在这批 ORF 的命中里排第几
"""
import os, re, subprocess, collections

A = "/home/zhangwenda/database/virus-db/acvirus_db"
ACC = re.compile(r"\b([A-Z]{1,4}[0-9]{5,}(?:\.[0-9]+)?)\b")

for f in ["taxa.txt", "vmr.tsv", "taxon_min_coverage.csv"]:
    p = os.path.join(A, f)
    print("=" * 110)
    print(p)
    if not os.path.exists(p):
        print("  缺")
        continue
    lines = open(p, encoding="utf-8", errors="replace").read().splitlines()
    print("  行数 %d ; 前 3 行:" % len(lines))
    for l in lines[:3]:
        print("    " + l[:150])
    hit = [l for l in lines if re.search(r"partiti|Partitivir", l)]
    print("  含 partitivir 的行 %d ; 前 10:" % len(hit))
    for l in hit[:10]:
        print("    " + l[:150])

# 收集 Partitiviridae 相关 accession
pacc = set()
for f in ["taxa.txt", "vmr.tsv"]:
    p = os.path.join(A, f)
    if not os.path.exists(p):
        continue
    for l in open(p, encoding="utf-8", errors="replace"):
        if re.search(r"partiti", l, re.I):
            pacc.update(ACC.findall(l))
print("\n从 taxon 表抽到 partitivirus 相关 accession %d 个: %s" % (len(pacc), sorted(pacc)[:20]))

# 全库 top200 命中
o = "/tmp/acvdb_top200.tsv"
cmd = ["/home/zhangwenda/biosoft/binary/diamond", "blastp", "-q", "/tmp/orf13.faa",
       "-d", os.path.join(A, "database.dmnd"), "--outfmt", "6", "qseqid", "sseqid",
       "pident", "length", "evalue", "bitscore", "--max-target-seqs", "200", "--quiet", "-o", o]
subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
rows = [l.split("\t") for l in open(o).read().splitlines()]
print("\ntop200 命中 %d 条" % len(rows))

# subject accession（去掉 _ORFindex）
sub_acc = set()
for r in rows:
    m = ACC.search(r[1])
    sub_acc.add(m.group(1) if m else r[1].split("_")[0])
inter = pacc & sub_acc
print("top200 里的 subject accession 数 %d ; 与 partitivirus accession 交集 %d : %s"
      % (len(sub_acc), len(inter), sorted(inter)[:20]))

# 也看 subject 是否命中已知的 partitivirus 基因组
known = {"WNO13897", "WNO13952", "WZH58866", "WZH58921", "QDH91349", "WZH58629", "WZH60815", "WZH58367"}
print("已知近亲 accession 前缀出现在 top200 里: %s" % (sorted(known & sub_acc) or "无"))
