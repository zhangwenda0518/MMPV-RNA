#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Cell-by-cell cross-check: my independent R1 vs the existing chimera.txt.

Reads /tmp/refresh_v67/<label>/chimera.txt (online=NEW product, v66f=OLD product)
and /tmp/indep_verify/logs/<label>.{new,old}.log, compares the 7 pair counts,
row counts and conflict rates. Prints OK/DIFF per cell. Read-only.
"""
import os
import re

LABELS = ["Alternaria", "amarum", "Aphis", "barbarum", "chinense", "Fusarium", "onekp", "ruthenicum"]
PAIRS = ["Realm-Kingdom", "Kingdom-Phylum", "Phylum-Class", "Class-Order",
         "Order-Family", "Family-Genus", "Genus-Species"]


def parse_chimera(path):
    out = {}
    with open(path, encoding="utf-8", errors="replace") as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < 11 or p[0] == "数据集":
                continue
            out[p[1]] = {
                "rows": int(p[2]),
                "pairs": [int(x) for x in p[3:10]],
                "conf": int(p[10]),
                "rate": p[11].strip(),
            }
    return out


def parse_log(path):
    rows = None
    r1 = {}
    conf = None
    in_r1 = False
    with open(path, encoding="utf-8", errors="replace") as f:
        for line in f:
            m = re.match(r"数据行数: (\d+)", line)
            if m:
                rows = int(m.group(1))
            if line.startswith("--- 口径 R1"):
                in_r1 = True
                m2 = re.search(r"矛盾行 (\d+)", line)
                conf = int(m2.group(1)) if m2 else None
                continue
            if line.startswith("--- 口径 R2"):
                in_r1 = False
            if in_r1:
                m3 = re.match(r"\s+(\S+)\s+(\d+)\s*$", line)
                if m3:
                    r1[m3.group(1)] = int(m3.group(2))
    return rows, r1, conf


def main():
    allok = True
    for lab in LABELS:
        ch = parse_chimera("/tmp/refresh_v67/%s/chimera.txt" % lab)
        for tag, logf in (("online", "new"), ("v66f", "old")):
            c = ch.get(tag)
            rows, r1, conf = parse_log("/tmp/indep_verify/logs/%s.%s.log" % (lab, logf))
            ok = True
            msgs = []
            if c is None:
                msgs.append("chimera.txt missing %s" % tag)
                ok = False
            else:
                if c["rows"] != rows:
                    msgs.append("rows %d!=%d" % (c["rows"], rows)); ok = False
                if c["conf"] != conf:
                    msgs.append("conf %d!=%d" % (c["conf"], conf)); ok = False
                for i, pr in enumerate(PAIRS):
                    if c["pairs"][i] != r1.get(pr):
                        msgs.append("%s %d!=%s" % (pr, c["pairs"][i], r1.get(pr))); ok = False
            allok = allok and ok
            print("%-11s %-7s %s %s" % (lab, tag, "OK" if ok else "DIFF", ";".join(msgs)))
    print("ALL_CELLS_MATCH" if allok else "HAS_DIFFS")


if __name__ == "__main__":
    main()
