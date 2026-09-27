#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
chk_agree.py —— 独立复算 v6.4 的 *_agree 列（不依赖 R，纯从原始输入重建）

口径（严格对齐 virus_classifier_analysis.R.cascade_v64）：
  long      : melt(standardized_*.tsv, id=c(contig_id,Tool), measure=TAX_LEVELS) 后按
              is_valid_value_vec 过滤（invalid_strs 见下），bad 值不进入 long
  n_tot     : long 中 (contig_id, Rank) 的行数 = 该阶元有合法取值的工具数
  cons_final: 产物 final_integrated_classification.tsv 的最终值，同样过 is_valid_value_vec
  n_ag      : 取值与最终值 tolower 后相等的工具数（最终值为空则恒 0）
  agree_str : n_ag>0 -> "n_ag/n_tot: tool1,tool2"（工具名按 ASCII 升序，C locale 同 R setorder）
              否则         -> "0/n_tot"

用法：python3 chk_agree.py <run_dir> <out_txt>
"""
import csv
import os
import sys
from collections import defaultdict

TAX_LEVELS = ["Realm", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
INVALID = {"-", "NA", "na", "N/A", "no rank", "undefined", "unknown", "null", "default", "Unclassified"}


def valid(v):
    if v is None:
        return False
    if v == "":
        return False
    return v not in INVALID


def read_tsv(path):
    with open(path, "r", encoding="utf-8", newline="") as fh:
        rd = csv.reader(fh, delimiter="\t", quotechar='"')
        head = next(rd)
        rows = [r for r in rd if r]
    return head, rows


def main():
    run_dir = sys.argv[1]
    out_txt = sys.argv[2]
    out = open(out_txt, "w", encoding="utf-8")

    def emit(s=""):
        out.write(s + "\n")

    # ---- 1. 读 7 个工具的标准化表，重建 long ----
    tools = []
    long = defaultdict(list)          # (contig, rank) -> [(tool, taxon)]
    contigs_by_tool = {}
    for fn in sorted(os.listdir(run_dir)):
        if not (fn.startswith("standardized_") and fn.endswith(".tsv")):
            continue
        tool = fn[len("standardized_"):-len(".tsv")]
        tools.append(tool)
        head, rows = read_tsv(os.path.join(run_dir, fn))
        idx = {name: i for i, name in enumerate(head)}
        missing = [l for l in TAX_LEVELS if l not in idx]
        if missing:
            emit("FATAL: %s 缺列 %s" % (fn, missing))
            return 1
        cset = set()
        for r in rows:
            c = r[idx["contig_id"]]
            cset.add(c)
            for lv in TAX_LEVELS:
                v = r[idx[lv]]
                if valid(v):
                    long[(c, lv)].append((tool, v))
        contigs_by_tool[tool] = cset
        emit("输入: %-12s 行 %6d  contig %6d" % (tool, len(rows), len(cset)))

    tools_sorted = sorted(tools)
    emit("工具 %d 个（ASCII 序）: %s" % (len(tools), ",".join(tools_sorted)))

    # ---- 2. 读产物最终值 ----
    head, rows = read_tsv(os.path.join(run_dir, "final_integrated_classification.tsv"))
    idx = {name: i for i, name in enumerate(head)}
    for lv in TAX_LEVELS:
        if lv not in idx or (lv + "_agree") not in idx:
            emit("FATAL: 产物缺列 %s / %s_agree" % (lv, lv))
            return 1

    n_cell = 0
    n_eq = 0
    by_rank = {lv: [0, 0] for lv in TAX_LEVELS}   # [比对格, 不一致格]
    samples = []
    agree_all_zero = {lv: 0 for lv in TAX_LEVELS}

    for r in rows:
        c = r[idx["contig_id"]]
        for lv in TAX_LEVELS:
            n_cell += 1
            by_rank[lv][0] += 1
            prod_val = r[idx[lv]]
            prod_agree = r[idx[lv + "_agree"]]
            pool = long.get((c, lv), [])
            n_tot = len(pool)
            if valid(prod_val):
                lo = prod_val.lower()
                hit = sorted(t for t, v in pool if v.lower() == lo)
                n_ag = len(hit)
                exp = "%d/%d: %s" % (n_ag, n_tot, ",".join(hit)) if n_ag > 0 else "0/%d" % n_tot
            else:
                n_ag = 0
                exp = "0/%d" % n_tot
                agree_all_zero[lv] += 1
            if exp == prod_agree:
                n_eq += 1
            else:
                by_rank[lv][1] += 1
                if len(samples) < 25:
                    samples.append("%s | %-8s | 产物值=%s | 值=%s | 复算=%s | long池=%s"
                                   % (c[:44], lv, prod_val, prod_agree, exp,
                                      ",".join("%s:%s" % (t, v) for t, v in pool)[:150]))

    emit("")
    emit("==== 逐格复算 *_agree：比对 %d 格，一致 %d 格，不一致 %d 格 ====" % (n_cell, n_eq, n_cell - n_eq))
    emit("%-9s %8s %8s %8s" % ("Rank", "格数", "不一致", "最终值空"))
    for lv in TAX_LEVELS:
        emit("%-9s %8d %8d %8d" % (lv, by_rank[lv][0], by_rank[lv][1], agree_all_zero[lv]))
    if samples:
        emit("")
        emit("---- 不一致样例（最多 25 条）----")
        for s in samples:
            emit(s)

    # ---- 3. 附带：一致率格式自检（分子分母关系） ----
    fmt_bad = 0
    fmt_sample = []
    for r in rows:
        for lv in TAX_LEVELS:
            s = r[idx[lv + "_agree"]]
            if s.count("/") < 1:
                fmt_bad += 1
                if len(fmt_sample) < 5:
                    fmt_sample.append("%s | %s -> %s" % (r[idx["contig_id"]][:44], lv, s))
    emit("")
    emit("格式自检: 非 a/b 结构的格子 %d 个%s" % (fmt_bad, ("  例: " + "; ".join(fmt_sample)) if fmt_sample else ""))

    # ---- 4. contig 覆盖 ----
    prod_contigs = set(r[idx["contig_id"]] for r in rows)
    emit("")
    emit("产物 contig %d 个" % len(prod_contigs))
    for t in tools_sorted:
        only_prod = len(prod_contigs - contigs_by_tool[t])
        emit("  %-12s 覆盖产物 contig 缺口中 %d 个（该工具无此行）" % (t, only_prod))
    out.close()
    return 0


if __name__ == "__main__":
    sys.exit(main())
