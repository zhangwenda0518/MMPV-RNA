#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
chk_rank.py — 同库排名对照：ACVirus 的 ORF 打同一个 RVDB 库，
看 partitivirus 与巨病毒（含 Biavirus/HG999358）谁排在前面、差多少 e-value
"""
import os, json, subprocess, urllib.request
from collections import defaultdict

DB = "/home/zhangwenda/database/virus-db/RVDB-v31/RVDB_viroids.diamond_db/U-RVDBv31.0-prot_unique.dmnd"
DIA = "/home/zhangwenda/biosoft/binary/diamond"

BIG = ["mimivir", "phycodnavir", "iridovir", "poxvir", "marseillevir", "pithovir",
       "tupanvirus", "pandoravir", "cotonvirus", "chlorovirus", "prymnesium",
       "schizomimivir", "megaviricetes", "mollivirus", "faustovirus", "kaumoebavirus",
       "medusavirus", "pacmanvirus", "cedratvirus", "orpheovirus", "yokavirus? "]
PART = ["partitivir", "chronawyn", "chronabint", "berere", "guapo"]

cmd = [DIA, "blastp", "-q", "/tmp/orf13.faa", "-d", DB, "--outfmt", "6",
       "qseqid", "sseqid", "pident", "length", "evalue", "bitscore",
       "--max-target-seqs", "60", "--quiet", "-o", "/tmp/rank60.tsv"]
print("$ " + " ".join(cmd))
r = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
print("exit=%d %s" % (r.returncode, r.stdout.decode()[:300]))

rows = []
subjects = set()
for line in open("/tmp/rank60.tsv"):
    p = line.rstrip("\n").split("\t")
    rows.append(p)
    subjects.add(p[1])
print("命中 %d 条，去重 subject %d 个" % (len(rows), len(subjects)))

# eutils 批量标注
info = {}
ids = sorted(subjects)
for i in range(0, len(ids), 150):
    chunk = ids[i:i + 150]
    url = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=protein&id="
           + ",".join(chunk) + "&retmode=json")
    try:
        d = json.load(urllib.request.urlopen(url, timeout=90))["result"]
        for k, v in d.items():
            if k == "uids":
                continue
            info[k] = (v.get("organism", "?") or "?", v.get("taxid", "?"), v.get("title", "") or "")
    except Exception as e:
        print("[eutils 失败 %d] %s" % (i, e))
print("标注到手 %d/%d 个" % (len(info), len(ids)))

per = defaultdict(list)
for p in rows:
    per[p[0]].append(p)

print("\n" + "=" * 132)
for q in sorted(per):
    hits = per[q]
    print("\n  [%s]" % q)
    best_part = best_big = None
    for h in hits:
        org, txid, title = info.get(h[1], ("?", "?", ""))
        blob = (org + " " + title).lower()
        isb = any(k in blob for k in BIG)
        isp = any(k in blob for k in PART)
        if isp and best_part is None:
            best_part = (h, org, title)
        if isb and best_big is None:
            best_big = (h, org, title)
    if best_part:
        h, org, title = best_part
        print("   最佳 partitivirus : rank %-3d %-14s pid=%-6s len=%-5s ev=%-11s bits=%-7s %s"
              % (hits.index(h) + 1, h[1], h[2], h[3], h[4], h[5], org[:40]))
    else:
        print("   最佳 partitivirus : 无")
    if best_big:
        h, org, title = best_big
        print("   最佳 巨病毒       : rank %-3d %-14s pid=%-6s len=%-5s ev=%-11s bits=%-7s %s"
              % (hits.index(h) + 1, h[1], h[2], h[3], h[4], h[5], org[:40]))
    else:
        print("   最佳 巨病毒       : 无（top60 内不出现）")
    print("   前 3 名：")
    for h in hits[:3]:
        org, txid, title = info.get(h[1], ("?", "?", ""))
        print("      %-14s pid=%-6s len=%-5s ev=%-11s bits=%-7s %s" % (h[1], h[2], h[3], h[4], h[5], org[:45]))
