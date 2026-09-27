#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
chk_acvdb_rank.py — 同一批 ORF 打 ACVirus 自己的库（acvirus_db/database.dmnd）
看：① 有没有 partitivirus 命中；② top 命中落在哪些基因组上
"""
import os, subprocess, re
from collections import defaultdict

ADB = "/home/zhangwenda/database/virus-db/acvirus_db"
DIA = "/home/zhangwenda/biosoft/binary/diamond"
FAA = os.path.join(ADB, "all_virus.faa")

for p in [os.path.join(ADB, "database.dmnd"), FAA]:
    print("  %s  %s bytes" % ("OK " if os.path.exists(p) else "缺 ", p))

cmd = [DIA, "blastp", "-q", "/tmp/orf13.faa", "-d", os.path.join(ADB, "database.dmnd"),
       "--outfmt", "6", "qseqid", "sseqid", "pident", "length", "evalue", "bitscore",
       "--max-target-seqs", "20", "--quiet", "-o", "/tmp/acvdb_rank.tsv"]
print("  $ " + " ".join(cmd))
r = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
print("  exit=%d %s" % (r.returncode, r.stdout.decode()[:300]))

rows = [l.rstrip("\n").split("\t") for l in open("/tmp/acvdb_rank.tsv")]
print("  命中 %d 条" % len(rows))
per = defaultdict(list)
for x in rows:
    per[x[0]].append(x)

# 把 subject ID 映射回 all_virus.faa 的注释
want = {x[1] for x in rows}
hdr = {}
if os.path.exists(FAA):
    with open(FAA, "r", encoding="utf-8", errors="replace") as fi:
        for line in fi:
            if line.startswith(">"):
                key = line[1:].split()[0]
                if key in want:
                    hdr[key] = line.rstrip("\n")[1:160]
                    if len(hdr) == len(want):
                        break
print("  subject 头信息拿到 %d/%d" % (len(hdr), len(want)))

PART = re.compile(r"partiti|Guapo|Berere|Chronawyn|Chronabint|Partitiviridae", re.I)
BIG = re.compile(r"mimivir|phycodnavir|tupanvirus|prymnesium|schizomimivir|iridovir|"
                 r"Prymnesium|Phaeocystis|Ostreococcus|Bathycoccus|cotonvirus|pandoravirus", re.I)

n_part = n_big = 0
print("\n" + "=" * 134)
for q in sorted(per):
    hits = per[q]
    print("\n  [%s]  top %d" % (q, min(len(hits), 6)))
    for i, h in enumerate(hits[:6]):
        txt = hdr.get(h[1], "")
        tag = "P" if PART.search(h[1] + txt) else ("B" if BIG.search(h[1] + txt) else "?")
        if tag == "P":
            n_part += 1
        if tag == "B":
            n_big += 1
        print("    #%-2d %-18s pid=%-6s len=%-5s ev=%-11s [%s] %s" % (i + 1, h[1], h[2], h[3], h[4], tag, txt[:78]))

allt = " ".join(hdr.values())
print("\n" + "=" * 134)
print("ACVirus 库内 top20 命中：partitivirus 标记 %d 条，巨病毒标记 %d 条" % (n_part, n_big))
print("抽 3 个 subject 头看命名法：")
for k in list(hdr)[:3]:
    print("   %s" % hdr[k])
