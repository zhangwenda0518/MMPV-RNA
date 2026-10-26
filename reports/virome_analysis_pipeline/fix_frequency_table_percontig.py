#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""fix_frequency_table_percontig.py — 把 `frequency_table.tsv` 的 `n_samples` 从
**属/代表序列级** 修正为 **逐 contig 级**。

为什么需要这个脚本（2026-10-13，用户：「你这发生率有点问题」）
--------------------------------------------------------------
`virome_pipeline.py::_write_frequency_table()` 旧实现按 `Rep_Accession` 分组：

    freq[acc].add(sample)          # acc = Rep_Accession（属/物种级代表序列）

但 `Rep_Accession` 在 `rescue_detection_summary.tsv` 里**只有 31 个不同值**
（`batch_virus_depth21.py::asm_df` 是 `group_by(["Sample", "taxid"])` 聚合的），
而写表时却**逐 contig_id 输出** → 同一代表序列下的 19–23 条 contig
**全部拿到同一个 n_samples**。实测后果：

    SRR15037478_clean_NODE_23  (cov 0.70×)  → 现行 272 (87.18%) / 真值 178
    CRR527049_clean_NODE_8     (cov 45.87×) → 现行 272 (87.18%) / 真值 186
    SRR23515444_clean_NODE_1   (cov 324.6×) → 现行 295 (94.55%) / 真值 198

真值来源 = `Segment_Accessions` 列（该样本下同 taxid 的**全部** contig，
`Rep_Accession` 恒为其成员 → 代表序列自身也仍被正确统计）。

本脚本对**已生成的结果树**做同等修正（不必重跑整条 pipeline），
上游 `virome_pipeline.py` 已同步修正（真正重跑时会直接产出正确表）。

⚠ 分母：本脚本**不改分母**（沿用原表的 `prevalence_% / n_samples` 反推值），
  只重算 `n_samples` 与 `prevalence_%` 的**分子**。分母必须由调用方确认。
⚠ 自检：修正后 n_samples 必须 **≤ 分母**，且逐 contig 值必须 **≤ 原属级值**
  （属级是并集，逐 contig 是子集，不可能更大）—— 违反即报错，不静默写出。

用法：
    "C:/Python312/python.exe" reports/virome_analysis_pipeline/fix_frequency_table_percontig.py \
        --summary <root>/10_Reports/rescue_detection_summary.tsv \
        --table   <root>/10_Reports/frequency_table.tsv \
        [--denom 312] [--dry-run]
"""
from __future__ import annotations

import argparse
import csv
import sys
from collections import defaultdict
from pathlib import Path


def build_per_contig(summary: Path):
    """从 summary 展开 `Segment_Accessions` → {contig_id: set(sample)}。"""
    freq = defaultdict(set)
    n_rep_only = 0
    with summary.open(encoding="utf-8", errors="ignore") as f:
        hdr = f.readline().rstrip("\n").split("\t")
        col = {h: i for i, h in enumerate(hdr)}
        if "Rep_Accession" not in col:
            sys.exit("❌ summary 缺 `Rep_Accession` 列")
        ac, sc, gc = col["Rep_Accession"], col.get("Sample", 0), col.get("Segment_Accessions")
        if gc is None:
            sys.exit("❌ summary 缺 `Segment_Accessions` 列（无法做逐 contig 展开）")
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) <= max(ac, sc, gc):
                continue
            s = p[sc].strip().strip('"')
            if not s.startswith(("SRR", "CRR", "ERR")):
                continue
            segs = [a.strip().strip('"') for a in (p[gc] or "").split(",") if a.strip()]
            if segs:
                for a in segs:
                    freq[a].add(s)
            else:
                freq[p[ac].strip().strip('"')].add(s)
                n_rep_only += 1
    if n_rep_only:
        print(f"   ⚠ {n_rep_only} 行缺 Segment_Accessions → 回退 Rep_Accession 口径")
    return freq


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--summary", required=True)
    ap.add_argument("--table", required=True)
    ap.add_argument("--denom", type=int, default=None,
                    help="发生率分母（队列样本总数）。不传则从原表 n_samples/prevalence_% 反推。")
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()

    summary, table = Path(a.summary), Path(a.table)
    for p in (summary, table):
        if not p.is_file():
            sys.exit(f"❌ 文件不存在：{p}")

    freq = build_per_contig(summary)
    rows = list(csv.DictReader(table.open(encoding="utf-8"), delimiter="\t"))
    if not rows:
        sys.exit("❌ 表为空")
    print(f"summary 展开 → {len(freq)} 条 contig；表 {len(rows)} 行")

    # ── 分母：优先 --denom，否则从原表反推（n_samples / prevalence_% 的比值）──
    if a.denom:
        DENOM = a.denom
        src = "--denom"
    else:
        cands = []
        for r in rows:
            try:
                ns, pv = int(r["n_samples"]), float(r["prevalence_%"])
            except (ValueError, KeyError):
                continue
            if ns > 0 and pv > 0:
                cands.append(ns / (pv / 100.0))
        if not cands:
            sys.exit("❌ 无法从原表反推分母，请用 --denom 显式给出")
        DENOM = round(sorted(cands)[len(cands) // 2])
        src = "反推自原表 n_samples/prevalence_%"
    print(f"分母 = {DENOM} 样本  [{src}]")

    # ── 重算 ──
    n_changed = 0
    out_rows = []
    for r in rows:
        cid = r["contig_id"]
        old = int(r["n_samples"]) if r.get("n_samples") not in (None, "") else 0
        new = len(freq.get(cid, set()))
        if new != old:
            n_changed += 1
        # ⚠ 逐 contig 值是属级（并集）的子集 → 不可能比原值大
        if new > old:
            sys.exit(f"❌ 逐 contig 值 {new} > 原属级值 {old}（{cid}）"
                     "—— 属级应为并集，逻辑不成立")
        if new > DENOM:
            sys.exit(f"❌ n_samples {new} > 分母 {DENOM}（{cid}）")
        r2 = dict(r)
        r2["n_samples"] = str(new)
        r2["prevalence_%"] = f"{new / DENOM * 100:.2f}"
        out_rows.append(r2)

    print(f"n_samples 变更 {n_changed} / {len(out_rows)} 行")
    if a.dry_run:
        print("(--dry-run: 未写盘)")
        for r in out_rows:
            if r["contig_id"] in freq:
                pass
        return

    bak = table.with_suffix(table.suffix + ".bak_percontig")
    if not bak.exists():
        bak.write_bytes(table.read_bytes())
        print(f"已备份 → {bak}")
    with table.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), delimiter="\t",
                           lineterminator="\n")
        w.writeheader()
        w.writerows(out_rows)
    print(f"OK  {table}")


if __name__ == "__main__":
    main()
