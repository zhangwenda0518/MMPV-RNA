#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""桥接：闸门前后的逐行差异，把「裸共识残留多少行」与「闸门动了多少行」对齐。

输入：
    /tmp/pregate_cascade/final_integrated_classification.tsv   闸门关闭
    /tmp/rc64b_cascade/final_integrated_classification.tsv     闸门开启（正式）
输出：Genus 变化行的动作分布（替换 / 置空 / 占位清理），以及与「工具层不自洽」判据的重合度。
"""
import csv
import os
from collections import Counter

PRE = os.environ.get("DIAG_PRE", "/tmp/pregate_cascade/final_integrated_classification.tsv")
POST = os.environ.get("DIAG_POST", "/tmp/rc64b_cascade/final_integrated_classification.tsv")
EMPTY = {"", "na", "n/a", "nan", "-", "none", "null"}


def norm(x):
    x = (x or "").strip().strip('"').strip()
    return "" if x.lower() in EMPTY else x


def load(p):
    d = {}
    with open(p, newline="", encoding="utf-8", errors="replace") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            d[norm(r["contig_id"])] = (norm(r.get("Family")), norm(r.get("Genus")), norm(r.get("Species")))
    return d


pre, post = load(PRE), load(POST)
print("闸门关闭行数 %d | 闸门开启行数 %d | 交集 %d"
      % (len(pre), len(post), len(set(pre) & set(post))))

act = Counter()
fam_changed = 0
for cid in set(pre) & set(post):
    pf, pg, ps = pre[cid]
    qf, qg, qs = post[cid]
    if pf != qf:
        fam_changed += 1
    if pg == qg:
        if ps != qs:
            act["仅 Species 变化"] += 1
        continue
    if not qg:
        act["Genus 被置空"] += 1
    elif not pg:
        act["Genus 被补上(原为空)"] += 1
    else:
        act["Genus 被替换"] += 1
    if pg and qg and pg.lower() == qg.lower():
        act["  其中仅大小写变化"] += 1
        continue
    if pf and qf and pf.lower() != qf.lower():
        act["  其中 Family 同时被改写"] += 1
        act["  其中 Family 被改写: %s -> %s" % (pf, qf)] += 1

print("Family 发生变化的行数 %d" % fam_changed)
print()
for k, v in act.most_common(20):
    print("  %-44s %6d" % (k, v))
