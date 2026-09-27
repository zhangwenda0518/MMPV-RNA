#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""票源级一致性检查：相邻等级的投票工具集合是否相交。

不依赖任何参照库，纯读 *_agree 列。
动机：参照库相容性判据抓不到"内部自洽但整条谱系由单一工具主导"的行。
只读，不改任何文件。
"""
import csv
import os
from collections import Counter, defaultdict

BASE = "/home/zhangwenda/MMPV-paper/goji-virome/02_novel_virus/RNA-Lycium_barbarum_out/05_Taxonomy"
INTEG = os.path.join(BASE, "Votus.integrated", "final_integrated_classification.tsv")

RANKS = ["Realm", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
PAIRS = [(RANKS[i], RANKS[i + 1]) for i in range(len(RANKS) - 1)]  # (高, 低)


def norm(v):
    if v is None:
        return None
    v = v.strip().strip('"').strip()
    if not v or v.upper() == "NA":
        return None
    return v


def parse_agree(v):
    """'3/5: CAT,diamond_lca,mmseqs' -> set of tools"""
    v = norm(v)
    if not v or ":" not in v:
        return None
    tools = v.split(":", 1)[1]
    s = {t.strip() for t in tools.split(",") if t.strip()}
    return s or None


def main():
    rows = []
    with open(INTEG, newline="", encoding="utf-8", errors="replace") as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        for r in rd:
            rows.append({k.strip().strip('"'): norm(v) for k, v in r.items()})
    print("行数: %d" % len(rows))

    pair_counts = Counter()
    disjoint_rows = 0
    rows_detail = []
    for d in rows:
        if not d.get("contig_id"):
            continue
        bad = []
        for hi, lo in PAIRS:
            sh = parse_agree(d.get(hi + "_agree"))
            sl = parse_agree(d.get(lo + "_agree"))
            if not sh or not sl:
                continue
            if sh.isdisjoint(sl):
                bad.append((hi, lo, ",".join(sorted(sh)), ",".join(sorted(sl))))
                pair_counts["%s <-> %s" % (hi, lo)] += 1
        if bad:
            disjoint_rows += 1
        rows_detail.append({"row": d, "bad": bad})

    print("")
    print("--- 票源不相交（相邻等级无共同投票工具）---")
    print("  受影响行数: %d (%.2f%%)" % (disjoint_rows, 100.0 * disjoint_rows / len(rows)))
    for k, v in pair_counts.most_common():
        print("  %-22s %6d" % (k, v))

    # 单源属级
    single_gen = 0
    single_gen_disjoint = 0
    for r in rows_detail:
        d = r["row"]
        sg = parse_agree(d.get("Genus_agree"))
        sf = parse_agree(d.get("Family_agree"))
        if sg and len(sg) == 1:
            single_gen += 1
            if sf and sg.isdisjoint(sf):
                single_gen_disjoint += 1
    print("")
    print("--- 属级单源（Genus_agree 只有 1 个工具）---")
    print("  行数: %d" % single_gen)
    print("  其中科级票源与属级票源完全不相交: %d" % single_gen_disjoint)

    # 全等级单源频谱
    print("")
    print("--- 各等级'唯一表态工具数=1'的行数 ---")
    for L in RANKS:
        c = 0
        for r in rows_detail:
            s = parse_agree(r["row"].get(L + "_agree"))
            if s and len(s) == 1:
                c += 1
        print("  %-9s %6d" % (L, c))

    # Biavirus 专项：自洽那批的票源
    print("")
    print("--- Biavirus 行（Genus=Biavirus）票源分布 ---")
    bv = [r for r in rows_detail if (r["row"].get("Genus") or "").lower() == "biavirus"]
    print("  总行数: %d" % len(bv))
    print("  Family_agree 分布 top 10:")
    for k, v in Counter((r["row"].get("Family_agree") or "NA") for r in bv).most_common(10):
        print("     %-32s %5d" % (k, v))
    print("  Genus_agree 分布 top 10:")
    for k, v in Counter((r["row"].get("Genus_agree") or "NA") for r in bv).most_common(10):
        print("     %-32s %5d" % (k, v))
    print("  Species_agree 分布 top 10:")
    for k, v in Counter((r["row"].get("Species_agree") or "NA") for r in bv).most_common(10):
        print("     %-32s %5d" % (k, v))
    print("  primary_tool 分布:")
    for k, v in Counter((r["row"].get("primary_tool") or "NA") for r in bv).most_common(10):
        print("     %-32s %5d" % (k, v))

    # 单源自洽族：Family/Genus/Species 三层票源全为同一单一工具
    print("")
    print("--- 三层票源全为同一单一工具（单源自洽）的行数 ---")
    same_source = 0
    same_source_acvirus = 0
    for r in rows_detail:
        d = r["row"]
        sf = parse_agree(d.get("Family_agree"))
        sg = parse_agree(d.get("Genus_agree"))
        ss = parse_agree(d.get("Species_agree"))
        sets = [s for s in (sf, sg, ss) if s]
        if len(sets) >= 2 and all(len(s) == 1 for s in sets) and len(set.union(*sets)) == 1:
            same_source += 1
            if list(sets[0])[0] == "ACVirus":
                same_source_acvirus += 1
    print("  行数: %d（其中工具=ACVirus: %d）" % (same_source, same_source_acvirus))

    # 两种判据的交集/并集
    fam_gen_ref = [r for r in rows_detail if any(b[0] == "Family" and b[1] == "Genus" for b in r["bad"])]
    print("")
    print("--- 参照库判据(科-属跨树) 与 票源判据(不相交) 的关系 ---")
    print("  科-属跨树行数（参照库）: %d" % len(fam_gen_ref))


if __name__ == "__main__":
    main()
