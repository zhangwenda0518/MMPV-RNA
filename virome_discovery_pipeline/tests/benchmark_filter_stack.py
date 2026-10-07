#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""benchmark_filter_stack.py — 过滤方法基准 (CDD → UniProt → graph 逐层归因)
在 step3_identification_eval 设计集 (2,344 contigs) 上评测, 246 运行。"""
import collections
import os
import sys

ROOT = "/home/zhangwenda/gp_eval"
N_TOTAL_POS = 1344
CLASSES = ["positive", "negative_A", "negative_B", "negative_C"]
CLS_CN = {"positive": "正类(病毒)", "negative_A": "negA随机",
          "negative_B": "negB保守域陷阱", "negative_C": "negC内源EVE"}

LAYERS = [
    ("L0 identification(候选池)", ROOT + "/ident/eval/eval_virus.all.candidate.fasta"),
    ("L1 +CDD域过滤", ROOT + "/filter/eval/cdd/eval_virus.all.candidate_cdd_filter.fasta"),
    ("L2 +UniProt同源过滤", ROOT + "/filter/eval/uniprot/eval_virus.all.candidate_cdd_filter_uniprot.fasta"),
    ("L3 +graph图裁决(最终)", ROOT + "/graph/eval/eval.virus.candidate.final.fasta"),
]

# ── 真值映射 ──
idmap, covmap = {}, {}
for line in open(ROOT + "/id_map.tsv"):
    f = line.rstrip("\n").split("\t")
    if f[0] == "eval_id":
        continue
    idmap[f[0]] = f[2]
    cov = ""
    for p in f[1].split("|"):
        if p.startswith("cov"):
            cov = p[3:]
    covmap[f[0]] = cov or "?"

report = []


def P(s=""):
    print(s)
    report.append(s)


def ids_of(path):
    if not os.path.isfile(path):
        return None
    return [l[1:].split()[0] for l in open(path) if l.startswith(">")]


P("=" * 78)
P("过滤方法基准 — step3_identification_eval 设计集 (2,344 contigs)")
P("  正类 1,344 (植物病毒, 覆盖度 40-100% 七档) | "
  "negA 随机 500 | negB 保守域陷阱 300 | negC 内源EVE 200")
P("  过滤栈: L1 CDD域(mmseqs) → L2 UniProt同源(diamond blastx) → L3 graph图裁决")
P("=" * 78)

P()
P("【表1】逐层归因: 每层之后各类剩余条数 (累计移除 = 上一层-本层)")
P("%-28s %7s %9s %9s %9s %9s %8s %8s" % (
    "层", "总数", "正类", "negA", "negB", "negC", "精确率", "召回"))
prev = None
for name, fa in LAYERS:
    ids = ids_of(fa)
    if ids is None:
        P("%-28s  <产物缺失>" % name)
        continue
    cnt = collections.Counter(idmap.get(i, "?") for i in ids)
    tp = cnt["positive"]
    fp = len(ids) - tp
    prec = tp / len(ids) if ids else 0.0
    rec = tp / N_TOTAL_POS
    P("%-28s %7d %9d %9d %9d %9d %7.1f%% %7.1f%%" % (
        name, len(ids), cnt["positive"], cnt["negative_A"],
        cnt["negative_B"], cnt["negative_C"], 100 * prec, 100 * rec))
    if prev is not None:
        removed = collections.Counter()
        for i in prev:
            c = idmap.get(i, "?")
            removed[c] += 1
        for i in ids:
            removed[idmap.get(i, "?")] -= 1
        det = ", ".join("%s %d" % (CLS_CN[c], -removed[c]) for c in CLASSES)
        P("%-28s       └ 移除: %s" % ("", det))
    prev = set(ids)

final_ids = prev
P()
P("【表2】最终过滤栈的二分类指标 (正类=病毒, 负类=全部非病毒)")
tp = sum(1 for i in final_ids if idmap.get(i) == "positive")
fp = len(final_ids) - tp
fn = N_TOTAL_POS - tp
tn = (2344 - N_TOTAL_POS) - fn
prec = tp / (tp + fp) if tp + fp else 0
rec = tp / N_TOTAL_POS
f1 = 2 * prec * rec / (prec + rec) if prec + rec else 0
P("  TP=%d  FP=%d  FN=%d  TN=%d" % (tp, fp, fn, tn))
P("  精确率 %.1f%% | 召回率 %.1f%% | F1 %.1f%%" % (100 * prec, 100 * rec, 100 * f1))

P()
P("【表3】正类召回 vs 覆盖度梯度 (graph 最终层)")
cov_stat = collections.defaultdict(lambda: [0, 0])
for i in idmap:
    if idmap[i] != "positive":
        continue
    cov = covmap.get(i, "?")
    cov_stat[cov][1] += 1
    if i in final_ids:
        cov_stat[cov][0] += 1
for cov in sorted(cov_stat, key=lambda x: (len(x), x)):
    kept, tot = cov_stat[cov]
    P("  cov%-6s: %d/%d = %5.1f%%" % (cov, kept, tot, 100 * kept / tot if tot else 0))

P()
P("【表4】三档负类的拒识率 (1 - 保留/总数, graph 最终层)")
for c in CLASSES[1:]:
    tot = sum(1 for v in idmap.values() if v == c)
    kept = sum(1 for i in final_ids if idmap.get(i) == c)
    P("  %s: 拒识 %d/%d = %5.1f%%" % (CLS_CN[c], tot - kept, tot, 100 * (tot - kept) / tot))

P()
P("【表5】graph 裁决分布 (逐类, verdict TSV)")
vd = collections.Counter()
vtsv = ROOT + "/graph/eval/eval.graph_verdict.tsv"
if os.path.isfile(vtsv):
    fh = open(vtsv, encoding="utf-8")
    head = fh.readline().rstrip("\n").split("\t")
    vi = head.index("verdict") if "verdict" in head else None
    for line in fh:
        f = line.rstrip("\n").split("\t")
        if vi is not None and len(f) > vi and f[0] in idmap:
            vd[(idmap[f[0]], f[vi])] += 1
    fh.close()
for cls in CLASSES:
    row = {k[1]: v for k, v in vd.items() if k[0] == cls}
    P("  %s: %s" % (CLS_CN[cls], dict(row)))

with open(ROOT + "/filter_benchmark_report.txt", "w", encoding="utf-8") as fh:
    fh.write("\n".join(report) + "\n")
print("\n报告已存: %s/filter_benchmark_report.txt" % ROOT)
