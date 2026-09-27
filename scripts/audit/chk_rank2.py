#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
chk_rank2.py — 只做标注：读 /tmp/rank60.tsv，逐 ORF 给出
「最佳 partitivirus 命中」vs「最佳巨病毒命中」的排名与 e-value
"""
import json, urllib.request, time
from collections import defaultdict

BIG = ["mimivir", "phycodnavir", "iridovir", "poxvir", "marseillevir", "pithovir",
       "tupanvirus", "pandoravir", "cotonvirus", "chlorovirus", "prymnesium",
       "schizomimivir", "megaviricetes", "mollivirus", "faustovirus", "kaumoebavirus",
       "medusavirus", "pacmanvirus", "cedratvirus", "orpheovirus", "virophage",
       "bialevir", "biavirus", "aureococcus", "chrysochromulina", "heterosigma",
       "phaeocystis", "eclovirus", "coccolitho"]
PART = ["partitivir", "chronawyn", "chronabint", "berere", "guapo", "amalgavir",
        "curvulari", "gammapartiti", "alphapartiti", "betapartiti", "deltapartiti"]

rows = [l.rstrip("\n").split("\t") for l in open("/tmp/rank60.tsv")]
subj = sorted({r[1] for r in rows})
print("命中 %d 条，去重 subject %d 个" % (len(rows), len(subj)))

info = {}
for i in range(0, len(subj), 40):
    chunk = subj[i:i + 40]
    url = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=protein&id="
           + ",".join(chunk) + "&retmode=json")
    ok = False
    for attempt in range(3):
        try:
            d = json.load(urllib.request.urlopen(url, timeout=90))["result"]
            for k, v in d.items():
                if k == "uids":
                    continue
                acc = v.get("accessionversion") or k
                info[acc] = (v.get("organism") or v.get("title", "?"), v.get("taxid", "?"),
                             v.get("title", "") or "")
            ok = True
            break
        except Exception as e:
            print("[重试 %d/3 chunk %d] %s" % (attempt + 1, i, e))
            time.sleep(3)
    if not ok:
        print("[失败 chunk %d] %s ..." % (i, ",".join(chunk[:3])))
print("标注到手 %d/%d" % (len(info), len(subj)))
if not info:
    raise SystemExit("标注全空，无法分类")

per = defaultdict(list)
for r in rows:
    per[r[0]].append(r)


def cls(acc):
    org, txid, title = info.get(acc, ("?", "?", ""))
    blob = (org + " " + title).lower()
    if any(k in blob for k in PART):
        return "P", org, title
    if any(k in blob for k in BIG):
        return "B", org, title
    return "?", org, title


print("\n" + "=" * 134)
for q in sorted(per):
    hits = per[q]
    bp = bb = None
    for i, h in enumerate(hits):
        c, org, title = cls(h[1])
        if c == "P" and bp is None:
            bp = (i + 1, h, org)
        if c == "B" and bb is None:
            bb = (i + 1, h, org)
    print("\n  [%s]  命中 %d 条" % (q, len(hits)))
    if bp:
        rk, h, org = bp
        print("   最佳 partitivirus  rank=%-3d %-14s pid=%-6s len=%-5s ev=%-11s  %s"
              % (rk, h[1], h[2], h[3], h[4], org[:52]))
    else:
        print("   最佳 partitivirus  无")
    if bb:
        rk, h, org = bb
        print("   最佳 巨病毒        rank=%-3d %-14s pid=%-6s len=%-5s ev=%-11s  %s"
              % (rk, h[1], h[2], h[3], h[4], org[:52]))
    else:
        print("   最佳 巨病毒        无（top60 内不出现）")
    for i, h in enumerate(hits[:5]):
        c, org, title = cls(h[1])
        print("     #%-2d %-14s pid=%-6s len=%-5s ev=%-11s [%s] %s" % (i + 1, h[1], h[2], h[3], h[4], c, org[:44]))
