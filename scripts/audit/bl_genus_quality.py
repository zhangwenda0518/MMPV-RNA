#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""属名质量核查：统计表里 Genus 的占位/垃圾值，并交叉 diff 表看替换是否引入垃圾。"""
import argparse
import csv
import sys
from collections import Counter

JUNK = set(x.lower() for x in [
    "virus", "viruses", "unclassified", "unidentified", "unknown", "uncultured",
    "na", "n/a", "none", "null", "other", "environmental", "samples", "s",
    "viridae", "virinae", "virales", "viricetes", "viricota", "riboviria",
])


def norm(v):
    if v is None:
        return ""
    v = str(v).strip()
    while v.startswith('"') and v.endswith('"') and len(v) >= 2:
        v = v[1:-1]
    return v.strip()


def read_col(path, col="Genus"):
    out = []
    with open(path, newline="", encoding="utf-8", errors="replace") as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        fn = [c.strip() for c in (rd.fieldnames or [])]
        if col not in fn:
            return out
        for r in rd:
            out.append((norm(r.get("contig_id")), norm(r.get(col))))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--table", required=True)
    ap.add_argument("--label", default="")
    ap.add_argument("--diff", default=None)
    args = ap.parse_args()

    rows = read_col(args.table)
    filled = [(k, g) for k, g in rows if g]
    junk = [(k, g) for k, g in filled if g.lower() in JUNK]
    print("[%s] 总行 %d | Genus 非空 %d | 占位垃圾值 %d" % (args.label, len(rows), len(filled), len(junk)))
    if junk:
        c = Counter(g for _, g in junk)
        print("  占位值分布: %s" % ", ".join("%s=%d" % (a, b) for a, b in c.most_common(12)))
        for k, g in junk[:6]:
            print("    %-50s %s" % (k[:50], g))

    # 疑似非属名：含空格、含大写开头但以 virus 结尾之外的形式、纯数字
    weird = [(k, g) for k, g in filled if " " in g or g.isupper() or g.isdigit()]
    print("  疑似非属名(含空格/全大写/纯数字): %d" % len(weird))
    for k, g in weird[:6]:
        print("    %-50s %s" % (k[:50], g))

    if args.diff:
        changed = []
        with open(args.diff, newline="", encoding="utf-8", errors="replace") as fh:
            rd = csv.DictReader(fh, delimiter="\t")
            for r in rd:
                if norm(r.get("column")) == "Genus":
                    changed.append((norm(r.get("contig_id")), norm(r.get("old")), norm(r.get("new"))))
        to_junk = [(k, a, b) for k, a, b in changed if b and b.lower() in JUNK]
        to_blank = [(k, a, b) for k, a, b in changed if not b]
        print("")
        print("[diff] Genus 变动 %d 行: 改为占位垃圾 %d | 置空 %d | 改为真实属名 %d"
              % (len(changed), len(to_junk), len(to_blank), len(changed) - len(to_junk) - len(to_blank)))
        for k, a, b in to_junk[:8]:
            print("    %-50s %s -> %s" % (k[:50], a, b))
    return 0


if __name__ == "__main__":
    sys.exit(main())
