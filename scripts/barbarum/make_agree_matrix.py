#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""make_agree_matrix.py — 从 final_integrated_classification.tsv 生成「分类工具赞成票」逐 contig 矩阵。

口径（与 fig_upset_classification.R 的 Species 版严格同构）：
  输入列 `<Rank>_agree` 形如 `"3/4: ACVirus,CAT,diamond_lca"`：
    · 分子 = 赞成票数，分母 = 在该 rank 上有投票的工具数；
    · 冒号后 = **投了赞成票**的工具集合（taxon 一致）。
  → 本脚本取「冒号后的工具集合」映射成 7 列 0/1 矩阵。

⚠ 与图①（UpsetTools）**不是同一口径**：
    图① = 「谁检出了这条 contig」（21,630 全集）；
    本图 = 「谁在该 rank 上投了赞成票」（子集，故行数更少）。

用法：
  python make_agree_matrix.py --rank Family
  python make_agree_matrix.py --rank Species --out upset_data/agree_matrix_Species.tsv
"""
import argparse
import csv
import sys
from pathlib import Path

# 7 个分类工具的固定顺序（与 upset_common.R 的 TOOLS7 严格一致，勿改）
TOOLS7 = ["diamond_lca", "mmseqs", "CAT", "metabuli", "genomad", "VITAP", "ACVirus"]

# 默认输入：本地镜像的整合表
DEFAULT_TSV = ("reports/_discovery_sci_20261012/data/10_Reports/"
               "final_integrated_classification.tsv")


def parse_agree(cell):
    """`"3/4: ACVirus,CAT,diamond_lca"` -> set of tools。无冒号 / 空 -> 空集。"""
    if not cell or ":" not in cell:
        return set()
    tail = cell.split(":", 1)[1].strip()
    if not tail:
        return set()
    return {t.strip() for t in tail.split(",") if t.strip()}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rank", default="Family",
                    help="阶元名（Realm/Kingdom/Phylum/Class/Order/Family/Genus/Species）")
    ap.add_argument("--input", default=DEFAULT_TSV)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    rank = args.rank
    col = f"{rank}_agree"
    inp = Path(args.input)
    if not inp.exists():
        sys.exit(f"输入不存在：{inp}")
    out = Path(args.out) if args.out else Path("upset_data") / f"agree_matrix_{rank}.tsv"

    rows = []
    n_zero = 0          # 该 rank 无任何赞成票（zero_no_vote）
    with inp.open(encoding="utf-8") as f:
        r = csv.DictReader(f, delimiter="\t")
        if col not in (r.fieldnames or []):
            sys.exit(f"列 {col} 不存在；可用列：{r.fieldnames}")
        for row in r:
            tools = parse_agree(row.get(col, ""))
            if not tools:
                n_zero += 1
                continue
            # 未知工具名要暴露，不能静默丢
            unknown = tools - set(TOOLS7)
            if unknown:
                sys.exit(f"出现未知工具名 {sorted(unknown)}（不在 TOOLS7）")
            rec = {"contig_id": row["contig_id"]}
            for t in TOOLS7:
                rec[t] = 1 if t in tools else 0
            rec["n_agree"] = sum(rec[t] for t in TOOLS7)
            rec["primary_tool"] = row.get("primary_tool", "")
            rows.append(rec)

    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, delimiter="\t",
                           fieldnames=["contig_id"] + TOOLS7 + ["n_agree", "primary_tool"])
        w.writeheader()
        w.writerows(rows)

    dist = {}
    for r in rows:
        dist[r["n_agree"]] = dist.get(r["n_agree"], 0) + 1
    print(f"rank = {rank}")
    print(f"  输入 {inp.name}：{len(rows) + n_zero} 行；其中 {n_zero} 行该 rank 无赞成票（不入图）")
    print(f"  入图 = {len(rows)} 条 vOTU x {len(TOOLS7)} 工具")
    print(f"  赞成工具数分布 = {dict(sorted(dist.items()))}")
    print(f"  写出 -> {out}")


if __name__ == "__main__":
    main()
