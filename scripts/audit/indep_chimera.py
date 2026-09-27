#!/usr/bin/env python3
"""独立实现的分类矛盾率计算器（用于交叉验证，不复用管线代码）

判据：相邻阶元对 (H 高, L 低)，Species<Genus<Family<Order<Class<Phylum<Kingdom<Realm(superkingdom)
  当且仅当 ①行内 H 有值 ②行内 L 有值 ③参照库里 L 有定型父级（对应 H 层级）非空
  才比较，两者不等即为该行该对上的矛盾。

为量化方法学敏感性，同时给出两种参照索引口径：
  R1 = 按 tax_name 在 rankedlineage.dmp 首次出现行取整条谱系（chimera_multi.py 口径）
  R2 = 按 L 层级列值聚合所有出现行的 H 层级取值集合（多父级集合口径）

尾随星号（软标注）与 NA/N/A 的处理：星号剥离后比较（口径与 chimera_multi.py 一致），
并单独统计含星号标注的行数。

用法: python3 indep_chimera.py <final_integrated_classification.tsv> [rankedlineage.dmp]
"""
import csv
import sys
from collections import Counter

TSV = sys.argv[1] if len(sys.argv) > 1 else ""
DMP = sys.argv[2] if len(sys.argv) > 2 else "/home/zhangwenda/database/taxonomy/rankedlineage.dmp"

LEVELS = ["Species", "Genus", "Family", "Order", "Class", "Phylum", "Kingdom", "Realm"]
DMP_COL = {"Realm": "superkingdom", "Kingdom": "kingdom", "Phylum": "phylum",
           "Class": "class", "Order": "order", "Family": "family",
           "Genus": "genus", "Species": "species"}
COLS = ["tax_id", "tax_name", "species", "genus", "family", "order",
        "class", "phylum", "kingdom", "superkingdom"]


def norm(v):
    """返回 (值或 None, 是否带尾随星号)"""
    if v is None:
        return None, False
    v = v.strip().strip('"').strip()
    star = False
    while v.endswith("*"):
        star = True
        v = v[:-1].strip()
    if not v or v.upper() in ("NA", "N/A"):
        return None, star
    return v, star


def parse_dmp(path):
    r1 = {}
    r2 = {}
    with open(path, encoding="utf-8", errors="ignore") as f:
        for line in f:
            parts = [p.strip().strip("|").strip() for p in line.split("\t|\t")]
            if len(parts) < 10:
                continue
            rec = dict(zip(COLS, parts))
            nm = rec["tax_name"].lower()
            if nm and nm not in r1:
                r1[nm] = {L: rec[DMP_COL[L]] for L in LEVELS}
            for i in range(len(LEVELS) - 1):
                low, high = LEVELS[i], LEVELS[i + 1]
                lv, hv = rec[DMP_COL[low]], rec[DMP_COL[high]]
                if lv and hv:
                    r2.setdefault((low, lv.lower()), set()).add(hv)
    return r1, r2


def main():
    r1, r2 = parse_dmp(DMP)
    print("参照库: R1(首次出现谱系) name=%d | R2(多父级集合) 键=%d" % (len(r1), len(r2)))

    pair_order = [("%s-%s" % (LEVELS[i + 1], LEVELS[i]), LEVELS[i + 1], LEVELS[i])
                  for i in range(len(LEVELS) - 1)][::-1]

    per_pair = {"R1": Counter(), "R2": Counter()}
    bad_rows = {"R1": set(), "R2": set()}
    total = 0
    starred_rows = 0

    with open(TSV, encoding="utf-8", errors="ignore") as f:
        rd = csv.DictReader(f, delimiter="\t", quotechar='"')
        for idx, raw in enumerate(rd):
            total += 1
            row = {}
            star_here = False
            for k, v in raw.items():
                val, star = norm(v)
                row[k.strip()] = val
                star_here = star_here or star
            if star_here:
                starred_rows += 1

            for name, high, low in pair_order:
                hv, lv = row.get(high), row.get(low)
                if not hv or not lv:
                    continue

                e = r1.get(lv.lower())
                if e is not None:
                    pv = (e.get(high) or "").strip()
                    if pv and pv != hv:
                        per_pair["R1"][name] += 1
                        bad_rows["R1"].add(idx)

                ps = r2.get((low, lv.lower()))
                if ps and hv not in ps:
                    per_pair["R2"][name] += 1
                    bad_rows["R2"].add(idx)

    print("输入: " + TSV)
    print("数据行数: %d | 含星号/软标注行: %d" % (total, starred_rows))
    for tag in ("R1", "R2"):
        n = len(bad_rows[tag])
        print("--- 口径 %s：矛盾行 %d (%.2f%%)" % (tag, n, n / total * 100 if total else 0.0))
        for name, _, _ in pair_order:
            print("    %-15s %6d" % (name, per_pair[tag][name]))
    print("SUMMARY\t%d\t%d\t%.4f\t%d\t%.4f" % (
        total, len(bad_rows["R1"]), len(bad_rows["R1"]) / total * 100 if total else 0.0,
        len(bad_rows["R2"]), len(bad_rows["R2"]) / total * 100 if total else 0.0))


if __name__ == "__main__":
    main()
