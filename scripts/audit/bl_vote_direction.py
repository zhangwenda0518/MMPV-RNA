#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""方向性判据：科层 vs 属层谁的工具支持更强。

对每个"科属打架"行，从 combined 长表数出两边各有多少工具支持：
  n_old = 报 old_family 且自带 Genus 的工具数（完整谱系票）
  n_new = 报 new_family 且自带 Genus 的工具数（即属所隐含的科）
输出 n_old vs n_new 的胜负分布，判断"以属修科"是否会二次污染。
"""
import argparse
import csv
import sys
from collections import Counter, defaultdict


def norm(v):
    if v is None:
        return None
    v = str(v).strip().strip('"').strip()
    while v.endswith("*"):
        v = v[:-1].strip()
    if not v or v.upper() in ("NA", "N/A", "-"):
        return None
    return v


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--log", required=True)
    ap.add_argument("--combined", required=True)
    args = ap.parse_args()

    rows = {}
    with open(args.log, newline="", encoding="utf-8", errors="replace") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            cid = norm(r.get("contig_id"))
            if cid:
                rows[cid] = (norm(r.get("old_family")), norm(r.get("new_family")))
    per = defaultdict(list)
    with open(args.combined, newline="", encoding="utf-8", errors="replace") as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        for r in rd:
            cid = norm(r.get("seq_name"))
            if cid in rows:
                per[cid].append((norm(r.get("tool")), norm(r.get("Family")), norm(r.get("Genus"))))

    win = Counter()
    for cid, (old, new) in rows.items():
        tools = per.get(cid) or []
        n_old = sum(1 for t, f, g in tools if f and old and f.lower() == old.lower() and g)
        n_old_bare = sum(1 for t, f, g in tools if f and old and f.lower() == old.lower() and not g)
        n_new = sum(1 for t, f, g in tools if f and new and f.lower() == new.lower() and g)
        if n_old > n_new:
            win["old_family_more_lineage_votes"] += 1
        elif n_new > n_old:
            win["new_family_more_lineage_votes"] += 1
        else:
            win["tie"] += 1
        if n_old == 0:
            win["old_family_zero_lineage_vote"] += 1
        if n_new == 0:
            win["new_family_zero_lineage_vote"] += 1
        if n_old_bare and n_old == 0:
            win["old_family_only_bare_votes"] += 1

    print("被改写行数: %d (combined 命中 %d)" % (len(rows), len(per)))
    print("")
    print("=== 工具级完整谱系票对比（old_family vs 属所隐含科）===")
    for k, v in win.most_common():
        print("  %-34s %6d" % (k, v))
    tot = len(rows)
    print("")
    print("  以科为准的比例 = (old>=new) / all = %d/%d = %.1f%%"
          % (win["old_family_more_lineage_votes"] + win["tie"], tot,
             100.0 * (win["old_family_more_lineage_votes"] + win["tie"]) / max(tot, 1)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
