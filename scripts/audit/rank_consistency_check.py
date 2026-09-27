#!/usr/bin/env python3
"""rank_consistency_check.py -- 跨阶元一致性残余错误量化

对 final_integrated_classification.tsv 逐行检查「Species 名蕴含的属」与
「Genus 列取值」是否自洽，以及 Genus 列是否被塞进物种级名字（含空格）。

只读，不改任何产物。用于 v6.7 刷新后残余分类错误的分层量化。

用法: python3 rank_consistency_check.py <final_integrated_classification.tsv>
"""
import csv
import sys
from collections import Counter


def clean(v):
    v = (v or "").strip()
    if v.startswith('"') and v.endswith('"') and len(v) >= 2:
        v = v[1:-1]
    return v.strip()


def main(path):
    with open(path, newline="") as fh:
        r = csv.reader(fh, delimiter="\t")
        header = next(r)
        idx = {name.strip().strip('"'): i for i, name in enumerate(header)}
        c_fam = idx.get("Family")
        c_gen = idx.get("Genus")
        c_spe = idx.get("Species")
        total = 0
        both = 0
        mismatch = 0
        multiword_genus = 0
        pairs = Counter()
        mw = Counter()
        for row in r:
            total += 1
            fam = clean(row[c_fam]) if c_fam is not None else ""
            gen = clean(row[c_gen]) if c_gen is not None else ""
            spe = clean(row[c_spe]) if c_spe is not None else ""
            if gen and gen != "NA" and " " in gen:
                multiword_genus += 1
                mw[gen] += 1
            if not gen or gen == "NA" or not spe or spe == "NA":
                continue
            both += 1
            first = spe.split()[0]
            gfirst = gen.split()[0]
            if first.lower() != gfirst.lower():
                mismatch += 1
                pairs[(gfirst, first, fam)] += 1
    rate = (100.0 * mismatch / both) if both else 0.0
    print("%s\ttotal=%d\tgen_spe_both=%d\tcross_rank_mismatch=%d\t%.2f%%\tmultiword_genus=%d"
          % (path.split("/")[-3] if "/" in path else path, total, both, mismatch, rate, multiword_genus))
    for (g, s, fam), n in pairs.most_common(5):
        print("    MISMATCH\t%d\tGenus=%s\tSpecies_prefix=%s\tFamily=%s" % (n, g, s, fam))
    for g, n in mw.most_common(5):
        print("    MULTIWORD_GENUS\t%d\t%s" % (n, g))


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit("usage: rank_consistency_check.py <final_integrated_classification.tsv>")
    main(sys.argv[1])
