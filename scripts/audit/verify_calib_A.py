#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""校验 apply_calib_A 的产物：校准表相对原表只应改动指定格子。

检查项：
 1. 行数一致；列数 = 原 20 + 6；列序正确
 2. 逐行逐列比较原始 20 列，统计变化格子，并核对变化位置只落在
    calib_action=blank 行的 Genus / Species 上
 3. 打印第一个意外变化的样例（若有）
"""
import csv
import os
import sys

BASE = "/home/zhangwenda/MMPV-paper/goji-virome/02_novel_virus/RNA-Lycium_barbarum_out/05_Taxonomy/Votus.integrated"
SRC = os.path.join(BASE, "final_integrated_classification.tsv")
OUT = os.path.join(BASE, "calibration_20260914", "final_integrated_classification.calibrated.tsv")

CALIB_COLS = ["calib_action", "calib_prev_genus", "calib_prev_species",
              "calib_evidence", "calib_flag_src", "calib_flag_single_ref"]


def load(path):
    with open(path, newline="", encoding="utf-8", errors="replace") as fh:
        rd = csv.reader(fh, delimiter="\t")
        head = [c.strip().strip('"') for c in next(rd)]
        body = []
        for r in rd:
            body.append(r)
    return head, body


def main():
    h1, b1 = load(SRC)
    h2, b2 = load(OUT)
    print("原表 列数 %d / 行数 %d" % (len(h1), len(b1)))
    print("产物 列数 %d / 行数 %d" % (len(h2), len(b2)))
    ok = True
    if len(b2) != len(b1):
        print("!! 行数不一致")
        ok = False
    if h2 != h1 + CALIB_COLS:
        print("!! 列名/列序不符")
        print("   产物头: %s" % h2)
        ok = False

    n = min(len(b1), len(b2))
    changed = {}
    unexpected = []
    blank_rows = 0
    blank_species_filled = 0
    for i in range(n):
        a, b = b1[i], b2[i]
        act = b[len(h1) + 0].strip()
        if act == "blank":
            blank_rows += 1
        for j, col in enumerate(h1):
            va = a[j] if j < len(a) else ""
            vb = b[j] if j < len(b) else ""
            if va != vb:
                changed[col] = changed.get(col, 0) + 1
                legal = (act == "blank" and col in ("Genus", "Species") and vb == "")
                if not legal:
                    unexpected.append((i + 2, col, act, va, vb))
        if act == "blank" and a[11] != b[11] and b[11] == "":
            blank_species_filled += 1

    print("")
    print("=== 相对原表发生变化的格子 ===")
    for col, c in sorted(changed.items(), key=lambda x: -x[1]):
        print("  %-12s %6d" % (col, c))
    print("")
    print("calib_action=blank 行数: %d" % blank_rows)
    print("其中 Species 由非空置空的格数: %d" % blank_species_filled)
    print("非法变化（非 blank 行的 Gen/species 或落空以外）: %d" % len(unexpected))
    for row in unexpected[:10]:
        print("   行 %d 列 %s action=%s: %r -> %r" % row)
    if unexpected:
        ok = False

    # 追加列自检
    print("")
    print("=== 追加列自检 ===")
    bad_prev = 0
    bad_act = 0
    ev_blank = 0
    for i in range(n):
        a, b = b1[i], b2[i]
        act = b[len(h1) + 0].strip()
        pg = b[len(h1) + 1].strip()
        ps = b[len(h1) + 2].strip()
        ev = b[len(h1) + 3].strip()
        if pg != (a[10] or "").strip():
            bad_prev += 1
        if ps != (a[11] or "").strip():
            bad_prev += 1
        if act not in ("", "blank", "review", "conflict_one_side"):
            bad_act += 1
        if act == "blank" and not ev:
            ev_blank += 1
    print("calib_prev_* 与原值不符: %d" % bad_prev)
    print("calib_action 取值非法: %d" % bad_act)
    print("blank 行缺 evidence: %d" % ev_blank)
    if bad_prev or bad_act or ev_blank:
        ok = False

    print("")
    print("结论: %s" % ("PASS 除指定格外逐字一致" if ok else "FAIL 见上"))
    return 0 if ok else 2


if __name__ == "__main__":
    sys.exit(main())
