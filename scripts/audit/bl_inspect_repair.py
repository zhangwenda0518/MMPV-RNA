#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""核验 repair_family_from_genus 的方向性：被改写行的"科"主张是裸科还是自带属的完整谱系。

对被改写行，从 combined 长表取回每个工具原始主张：
  - old_family 被哪些工具报了；其中多少工具同时给了 Genus（完整谱系 = 强主张）
  - new_family 被哪些工具报了；其中多少工具同时给了 Genus
  - 有没有工具报了 Genus 却把 Family 归到 old_family 一侧（真冲突）
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
    ap.add_argument("--sample", type=int, default=12)
    args = ap.parse_args()

    rows = {}
    with open(args.log, newline="", encoding="utf-8", errors="replace") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            cid = norm(r.get("contig_id"))
            if cid:
                rows[cid] = {"old": norm(r.get("old_family")), "new": norm(r.get("new_family")),
                             "genus": norm(r.get("genus")), "action": norm(r.get("action"))}
    print("被改写行数: %d" % len(rows))

    per = defaultdict(dict)
    with open(args.combined, newline="", encoding="utf-8", errors="replace") as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        cols = [c.strip().strip('"') for c in (rd.fieldnames or [])]
        cmap = {c: c for c in cols}
        cn = cmap.get("seq_name")
        ck = cmap.get("tool")
        cf = cmap.get("Family")
        cg = cmap.get("Genus")
        cs = cmap.get("Species")
        for r in rd:
            cid = norm(r.get(cn))
            if cid in rows:
                per[cid][norm(r.get(ck))] = {
                    "Family": norm(r.get(cf)), "Genus": norm(r.get(cg)), "Species": norm(r.get(cs)),
                }
    print("combined 命中 %d / %d" % (len(per), len(rows)))

    stat = Counter()
    detail = []
    for cid, meta in rows.items():
        tools = per.get(cid) or {}
        old, new, gen = meta["old"], meta["new"], meta["genus"]
        fam_old_with_genus = [t for t, v in tools.items()
                              if v["Family"] and old and v["Family"].lower() == old.lower() and v["Genus"]]
        fam_old_bare = [t for t, v in tools.items()
                        if v["Family"] and old and v["Family"].lower() == old.lower() and not v["Genus"]]
        fam_new_with_genus = [t for t, v in tools.items()
                              if v["Family"] and new and v["Family"].lower() == new.lower() and v["Genus"]]
        gen_reporters = [t for t, v in tools.items() if v["Genus"] and gen and v["Genus"].lower() == gen.lower()]
        gen_fam_matches_new = [t for t in gen_reporters
                               if tools[t]["Family"] and new and tools[t]["Family"].lower() == new.lower()]
        gen_fam_matches_old = [t for t in gen_reporters
                               if tools[t]["Family"] and old and tools[t]["Family"].lower() == old.lower()]
        if fam_old_with_genus:
            stat["old_family_also_has_genus"] += 1
        elif fam_old_bare:
            stat["old_family_only_bare"] += 1
        else:
            stat["old_family_no_tool_found"] += 1
        if fam_new_with_genus:
            stat["new_family_also_has_genus"] += 1
        if gen_reporters and len(gen_fam_matches_new) == len(gen_reporters):
            stat["genus_tool_family_all_match_new"] += 1
        elif gen_reporters:
            stat["genus_tool_family_partly_match_new"] += 1
        if gen_fam_matches_old:
            stat["genus_tool_family_matches_old"] += 1
        detail.append((cid, old, new, gen, fam_old_bare, fam_old_with_genus,
                       gen_reporters, gen_fam_matches_new, gen_fam_matches_old, tools))

    print("")
    print("=== 方向性统计 ===")
    for k, v in stat.most_common():
        print("  %-34s %6d" % (k, v))

    print("")
    print("=== 样本（%d 行）===" % args.sample)
    # 优先展示 old_family 自带属的行（潜在反向风险）
    detail.sort(key=lambda x: (len(x[5]) == 0, -len(x[5])))
    for d in detail[:args.sample]:
        cid, old, new, gen = d[0], d[1], d[2], d[3]
        print("")
        print("%s  len/cov=%s" % (cid, cid.split("length_")[-1] if "length_" in cid else "?"))
        print("  改写: Family %s -> %s | Genus %s" % (old, new, gen))
        for t in sorted(d[9]):
            v = d[9][t]
            print("    %-14s F=%-22s G=%-22s S=%s" % (t, v["Family"] or "-", v["Genus"] or "-", v["Species"] or "-"))
    print("")
    return 0


if __name__ == "__main__":
    sys.exit(main())
