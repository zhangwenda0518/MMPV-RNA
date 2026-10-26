#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
virus_spot_dist.py — 按病毒种统计每个 run 的 STAT spot 分布
===========================================================

输入
----
  --best-summary  all_viruses.best.summary.tsv（goji 已知病毒 best 表）
                  只取 Sample + taxid + Species_NCBI 三列作为病毒清单与分母
  --stat-raw      plant_virus_spot_dist.py 抓下来的 STAT JSON 缓存目录
                  （文件名 = {acc}.json，内容 = [payload]）

spot 口径
---------
  对每个目标病毒 taxid，在 STAT tax_table 里取：
    self_spots  = 该 taxid 那一行的 self_count（直接命中，不含后代）
    sub_spots   = 该 taxid 那一行的 total_count（含全部后代）
                  ⚠ 只有当该 taxid 出现在表里才有；表里没这行则为 0

  self_spots 是「这个种自己的 k-mer」，sub_spots 是「这个种及其下游株系」。
  类病毒（PSTVd 等）通常没有下游株系节点，两者相等。

输出
----
  virus_spot_summary.tsv   每病毒一行：检出 run 数 / 分位数 / 阈值通过率
  virus_spot_per_sample.tsv 每 run × 病毒 长表（spot > 0 才有行）
  VIRUS_SPOT_REPORT.md     人读报告

用法
----
    python virus_spot_dist.py \
        --best-summary ~/MMPV-paper/goji-virome/03_known_virus/known_virus_all_v2/01_detection/summary/all_viruses.best.summary.tsv \
        --stat-raw ~/MMPV-paper/pv_spot_work/final/stat_raw \
        -o ~/MMPV-paper/pv_spot_work/virus_spot
"""

from __future__ import annotations

import argparse
import csv
import glob
import json
import os
import statistics
import sys
from collections import defaultdict


def quantiles(vals: list[float]) -> dict[str, float]:
    if not vals:
        return {}
    v = sorted(vals)
    n = len(v)

    def q(p: float) -> float:
        if n == 1:
            return v[0]
        i = (n - 1) * p
        lo, hi = int(i), min(int(i) + 1, n - 1)
        return v[lo] + (v[hi] - v[lo]) * (i - lo)

    return {"n": n, "min": v[0], "p10": q(0.10), "p25": q(0.25), "median": q(0.5),
            "p75": q(0.75), "p90": q(0.90), "p95": q(0.95), "max": v[-1],
            "mean": statistics.fmean(v), "sum": sum(v)}


THRESHOLDS = [1, 10, 50, 100, 500, 1000, 5000, 10000]


def main() -> int:
    ap = argparse.ArgumentParser(description="按病毒种统计每 run 的 STAT spot 分布")
    ap.add_argument("--best-summary", required=True)
    ap.add_argument("--stat-raw", required=True)
    ap.add_argument("-o", "--outdir", required=True)
    ap.add_argument("--top", type=int, default=0,
                    help="只报告检出 run 数最多的前 N 个病毒（0 = 全部）")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    # ---- 1) 病毒清单与分母：来自 best.summary.tsv ----
    wanted: dict[int, str] = {}          # taxid -> 名称
    best_samples: dict[int, set[str]] = defaultdict(set)
    all_samples: set[str] = set()
    with open(args.best_summary, "r", encoding="utf-8-sig", errors="replace") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            s = (r.get("Sample") or "").strip()
            if not s:
                continue
            all_samples.add(s)
            try:
                tid = int((r.get("taxid") or "").strip())
            except ValueError:
                continue
            nm = (r.get("Species_NCBI") or r.get("Adjusted_Species") or "").strip()
            wanted.setdefault(tid, nm or f"taxid:{tid}")
            best_samples[tid].add(s)
    print(f"[in] best.summary: {len(wanted)} 个病毒种, {len(all_samples)} 个样本")

    # ---- 2) 扫 STAT 缓存 ----
    files = sorted(glob.glob(os.path.join(args.stat_raw, "*.json")))
    print(f"[in] STAT 缓存 {len(files)} 个文件")

    # taxid -> {acc: (self, sub)}
    per_virus: dict[int, dict[str, tuple[float, float]]] = defaultdict(dict)
    acc_meta: dict[str, dict] = {}
    n_ok = 0
    for fp in files:
        acc = os.path.basename(fp)[:-5]
        try:
            obj = json.load(open(fp, "r", encoding="utf-8"))
        except Exception:
            continue
        payload = obj[0] if isinstance(obj, list) else obj
        if not isinstance(payload, dict):
            continue
        n_ok += 1
        tot = payload.get("tax_totals", {}) or {}
        acc_meta[acc] = {
            "analysed": int(tot.get("analysed") or tot.get("total") or 0),
            "identified": int(tot.get("identified") or 0),
            # 该 run 的 STAT 用的是哪一版分类库 —— 决定新 taxid 能否出现
            "taxver": (tot.get("last_update") or "")[:4],
        }
        for row in payload.get("tax_table", []) or []:
            try:
                tid = int(row.get("tax_id"))
            except (TypeError, ValueError):
                continue
            if tid not in wanted:
                continue
            try:
                sc = float(row.get("self_count") or 0)
                tc = float(row.get("total_count") or 0)
            except (TypeError, ValueError):
                continue
            per_virus[tid][acc] = (sc, tc)
    print(f"[in] 可用 STAT payload {n_ok} 个")

    # ---- 3) 汇总 ----
    stat_samples = set(acc_meta)
    rows_summary = []
    long_rows = []
    for tid, nm in wanted.items():
        hits = per_virus.get(tid, {})
        self_vals = [v[0] for v in hits.values() if v[0] > 0]
        sub_vals = [v[1] for v in hits.values() if v[1] > 0]
        q = quantiles(self_vals)
        in_stat = len(stat_samples & best_samples.get(tid, set()))
        rec = {
            "taxid": tid, "virus": nm,
            "n_best_summary": len(best_samples.get(tid, set())),
            "n_in_stat_pool": in_stat,
            "n_spot_pos": len(self_vals),
            "sum_spot": q.get("sum", 0),
            "min": q.get("min", 0), "p10": q.get("p10", 0), "p25": q.get("p25", 0),
            "median": q.get("median", 0), "p75": q.get("p75", 0),
            "p90": q.get("p90", 0), "p95": q.get("p95", 0),
            "max": q.get("max", 0), "mean": q.get("mean", 0),
            "median_sub": quantiles(sub_vals).get("median", 0),
            "max_sub": quantiles(sub_vals).get("max", 0),
            # 该 taxid 是否出现在 STAT 表里（区分「没检出」与「分类库太旧没这个节点」）
            "n_taxid_present": sum(1 for acc in acc_meta if acc in hits),
        }
        # 「组装确认但 STAT 里连节点都没有」的样本，其分类库年份分布
        absent = sorted((best_samples.get(tid, set()) & stat_samples) - set(hits))
        vers = defaultdict(int)
        for acc in absent:
            vers[acc_meta[acc]["taxver"] or "?"] += 1
        rec["absent_taxver"] = ",".join(f"{k}:{v}" for k, v in sorted(vers.items()))
        rec["n_absent_no_node"] = len(absent)
        for t in THRESHOLDS:
            rec[f"ge{t}"] = sum(1 for x in self_vals if x >= t)
        rows_summary.append(rec)
        for acc, (sc, tc) in hits.items():
            if sc <= 0 and tc <= 0:
                continue
            long_rows.append({"taxid": tid, "virus": nm, "acc": acc,
                              "self_spots": sc, "sub_spots": tc,
                              "analysed": acc_meta.get(acc, {}).get("analysed", 0)})

    rows_summary.sort(key=lambda r: (-r["n_spot_pos"], -r["sum_spot"]))
    if args.top:
        rows_summary = rows_summary[: args.top]
        keep = {r["taxid"] for r in rows_summary}
        long_rows = [r for r in long_rows if r["taxid"] in keep]

    # ---- 4) 写表 ----
    scols = ["taxid", "virus", "n_best_summary", "n_in_stat_pool", "n_taxid_present",
             "n_spot_pos", "n_absent_no_node", "absent_taxver",
             "sum_spot", "min", "p10", "p25", "median", "p75", "p90", "p95",
             "max", "mean", "median_sub", "max_sub"] + [f"ge{t}" for t in THRESHOLDS]
    sp = os.path.join(args.outdir, "virus_spot_summary.tsv")
    with open(sp, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=scols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in rows_summary:
            w.writerow({k: (f"{v:,.2f}" if isinstance(v, float) else v)
                        for k, v in r.items()})

    lp = os.path.join(args.outdir, "virus_spot_per_sample.tsv")
    with open(lp, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["taxid", "virus", "acc", "self_spots",
                                          "sub_spots", "analysed"],
                           delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in sorted(long_rows, key=lambda x: (x["taxid"], -x["self_spots"])):
            w.writerow(r)

    # ---- 5) 报告 ----
    rp = os.path.join(args.outdir, "VIRUS_SPOT_REPORT.md")
    with open(rp, "w", encoding="utf-8") as f:
        f.write("# 按病毒种的 STAT spot 分布（goji）\n\n")
        f.write(f"- 病毒清单来源：`{os.path.basename(args.best_summary)}`"
                f"（{len(wanted)} 种 / {len(all_samples)} 样本）\n")
        f.write(f"- STAT 缓存：{n_ok} 个 run（SRR/ERR/DRR；CRR 无 STAT，不在此表）\n")
        f.write("- spot = STAT `self_count`（该种节点直接命中，不含后代株系）\n")
        f.write("- `sub_spots` = `total_count`（含该种下游株系）；"
                "无株系节点时与 self 相等\n")
        f.write("- ⚠ `n_best_summary` 是组装流程确认的样本数；"
                "`n_spot_pos` 是 STAT 里 spot>0 的 run 数，两者口径不同\n")
        f.write("- ⚠ **`n_taxid_present` 是「STAT 表里存在这个 taxid 节点」的 run 数**。"
                "NCBI STAT 是各 run 提交时用自己的分类库版本算的，"
                "新建立的 taxid 在旧 run 的 STAT 里根本不存在 —— "
                "此时 spot=0 是「分类库没这个节点」，不是「没检出病毒」。\n\n")

        # 分类库版本导致的漏检，单列出来
        stale = [r for r in rows_summary
                 if r["n_best_summary"] > 0 and r["n_taxid_present"] < r["n_best_summary"] * 0.5]
        if stale:
            f.write("## ⚠ 分类库版本导致的系统性缺测\n\n")
            f.write("下表病毒的组装确认样本数远多于「STAT 里存在该 taxid」的 run 数，"
                    "说明多数样本的 STAT 是用旧分类库跑的，本表 spot 分布**严重低估**。\n\n")
            f.write("| 病毒 | taxid | best表样本 | STAT中存在该节点 | 有spot | 无节点样本的分类库年份 |\n")
            f.write("|---|---:|---:|---:|---:|---|\n")
            for r in sorted(stale, key=lambda x: -x["n_best_summary"])[:15]:
                f.write(f"| {r['virus']} | {r['taxid']} | {r['n_best_summary']} | "
                        f"{r['n_taxid_present']} | {r['n_spot_pos']} | "
                        f"{r['absent_taxver'] or '-'} |\n")
            f.write("\n")
            f.write("年份是该 run 的 STAT payload 里 `tax_totals.last_update`。"
                    "同一病毒在不同年份的库里有/没有该节点，说明缺测纯粹由建库时间造成。\n\n")

        f.write("## 概览（按 STAT 检出 run 数排序）\n\n")
        f.write("| 病毒 | taxid | best表样本 | STAT检出 | spot合计 | 中位数 | p90 | max |\n")
        f.write("|---|---:|---:|---:|---:|---:|---:|---:|\n")
        for r in rows_summary[:40]:
            f.write(f"| {r['virus']} | {r['taxid']} | {r['n_best_summary']} | "
                    f"{r['n_spot_pos']} | {r['sum_spot']:,.0f} | "
                    f"{r['median']:,.1f} | {r['p90']:,.1f} | {r['max']:,.0f} |\n")
        if len(rows_summary) > 40:
            f.write(f"\n_（共 {len(rows_summary)} 种，此处列前 40）_\n")
        f.write("\n")

        # 重点病毒逐个展开
        focus = [r for r in rows_summary
                 if r["n_spot_pos"] >= 10][:12]
        if focus:
            f.write("## 重点病毒：spot 阈值通过率\n\n")
            f.write("分母 = 该病毒在 STAT 里 spot>0 的 run 数。\n\n")
            f.write("| 病毒 | 检出run | " + " | ".join(f"≥{t}" for t in THRESHOLDS) + " |\n")
            f.write("|---|---:|" + "---:|" * len(THRESHOLDS) + "\n")
            for r in focus:
                cells = " | ".join(f"{r[f'ge{t}']} ({100.0*r[f'ge{t}']/r['n_spot_pos']:.0f}%)"
                                   for t in THRESHOLDS)
                f.write(f"| {r['virus']} | {r['n_spot_pos']} | {cells} |\n")
            f.write("\n")

            f.write("## 重点病毒：分位数明细\n\n")
            f.write("| 病毒 | n | min | p10 | p25 | median | p75 | p90 | p95 | max |\n")
            f.write("|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|\n")
            for r in focus:
                f.write(f"| {r['virus']} | {r['n_spot_pos']} | {r['min']:,.0f} | "
                        f"{r['p10']:,.1f} | {r['p25']:,.1f} | {r['median']:,.1f} | "
                        f"{r['p75']:,.1f} | {r['p90']:,.1f} | {r['p95']:,.1f} | "
                        f"{r['max']:,.0f} |\n")
            f.write("\n")

        # PSTVd / GCVA 专项
        for tid, label in ((12892, "PSTVd"), (3112997, "GCVA")):
            r = next((x for x in rows_summary if x["taxid"] == tid), None)
            if not r:
                continue
            f.write(f"## {label}（taxid {tid}）专项\n\n")
            f.write(f"- best.summary 确认 {r['n_best_summary']} 个样本；"
                    f"STAT 表里存在该 taxid 节点的 run {r['n_taxid_present']} 个；"
                    f"其中 spot>0 的 **{r['n_spot_pos']}** 个\n")
            if r["n_absent_no_node"]:
                f.write(f"- ⚠ 另有 **{r['n_absent_no_node']}** 个组装确认样本的 STAT 里"
                        f"**完全没有这个 taxid 节点**（分类库年份：{r['absent_taxver']}），"
                        f"它们的 spot=0 是缺测，不是阴性\n")
            f.write(f"- spot 分位数：min={r['min']:,.0f} / p25={r['p25']:,.1f} / "
                    f"median={r['median']:,.1f} / p75={r['p75']:,.1f} / "
                    f"p90={r['p90']:,.1f} / max={r['max']:,.0f}\n")
            f.write(f"- spot 合计 {r['sum_spot']:,.0f}\n")
            f.write(f"- 阈值通过：")
            f.write(" / ".join(f"≥{t}: {r[f'ge{t}']} ({100.0*r[f'ge{t}']/r['n_spot_pos']:.0f}%)"
                               for t in THRESHOLDS))
            f.write("\n\n")
            sub = sorted([x for x in long_rows if x["taxid"] == tid],
                         key=lambda x: -x["self_spots"])[:10]
            if sub:
                f.write("| run | self spot | sub spot | analysed |\n|---|---:|---:|---:|\n")
                for x in sub:
                    f.write(f"| {x['acc']} | {x['self_spots']:,.0f} | "
                            f"{x['sub_spots']:,.0f} | {int(x['analysed']):,} |\n")
                f.write("\n")

    print(f"[done] {sp}\n[done] {lp}\n[done] {rp}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
