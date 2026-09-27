#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""打印 VMR_MSL41.tsv 表头及索引，并给出 Potyvirus rapae 一行的真实取值"""
import csv, os
VMR = os.path.expanduser("~/plant_virus_db/1.virus-host_db/B-ictv/VMR_MSL41.tsv")
with open(VMR, newline="", encoding="utf-8", errors="replace") as f:
    rd = csv.reader(f, delimiter="\t")
    hdr = next(rd)
    print("列数=%d" % len(hdr))
    for i, h in enumerate(hdr):
        print("  [%2d] %s" % (i, h))
    n = 0
    for r in rd:
        if len(r) > 17 and r[17].strip().strip('"') in ("Potyvirus rapae", "Rimosavirus zeae", "Potexvirus cymbidii"):
            print("\n样例:")
            for i, (h, v) in enumerate(zip(hdr, r)):
                print("  [%2d] %-28s = %s" % (i, h, v))
            n += 1
            if n >= 2:
                break
