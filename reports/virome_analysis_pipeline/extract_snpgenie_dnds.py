#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
extract_snpgenie_dnds.py
========================
把 `virus-SNPGenie` 下所有 GCVA 样本的 `product_results.txt` 汇总成**一张长表**，
供 `make_sci_composite_figs.py --ef self` 的自绘面板 E 使用。

为什么需要它
------------
`04_post_analysis/` 里 snpgenie 已经出的成品图（`Fig05_dNdS_per_gene` 之类）**不能直接嵌进主图**：
它们是白底灰网格 + 自带字体的位图，与主图 A–D/G 的 Okabe-Ito 语言不统一（见 rev.h 决策）。
但**数值都在** `product_results.txt` 里（15 列），所以重新汇总成表 → 自绘。

product -> 基因名的映射**不是猜的**，来自权威 GTF：
    <Final_Virus_Results2>/virus-annotations/OR489165.1.gtf
      CDS 436-1950   transcript_id "1" gene_name "N"
      CDS 2499-4046  transcript_id "2" gene_name "P"
      CDS 4305-5024  transcript_id "3" gene_name "P4"
      CDS 5480-6340  transcript_id "4" gene_name "M"
      CDS 6742-7368  transcript_id "5" gene_name "G"
      CDS 7562-14344 transcript_id "6" gene_name "L"
（CDS 长度 1515/1548/720/861/627/6783 nt 与 `N_sites + S_sites` 反推的
 1512/1545/717/858/624/6780 互相印证 —— 映射正确，不是"看着像"。）

用法
----
    export LD_LIBRARY_PATH=/home/zhangwenda/mambaforge/lib:$LD_LIBRARY_PATH
    python extract_snpgenie_dnds.py \
        --base /home/zhangwenda/MMPV-paper/goji-virome/01_data/sra_rna.data4/\
test/Final_Virus_Results2/virus-SNPGenie \
        --out  <composite 输出目录>/_data/gcva_dnds.tsv

    实测（2026-10-10）：`files=161 rows=966` = **161 样本 × 6 基因**。
输出列
------
    sample gene product piN piS N_diffs_vs_ref S_diffs_vs_ref N_sites S_sites
    mean_dN_vs_ref mean_dS_vs_ref mean_gdiv_polymorphic mean_N_gdiv mean_S_gdiv
"""
from __future__ import annotations

import argparse
import csv
import glob
import os
import sys

# 默认值 = barbarum 实测路径；换物种/换批次时用 --base / --pattern 覆盖
DEFAULT_BASE = ("/home/zhangwenda/MMPV-paper/goji-virome/01_data/sra_rna.data4/"
                "test/Final_Virus_Results2/virus-SNPGenie")
DEFAULT_PATTERN = "*Cytorhabdovirus*"

GENE_OF_PRODUCT = {"1": "N", "2": "P", "3": "P4", "4": "M", "5": "G", "6": "L"}

COLS = ["piN", "piS", "N_diffs_vs_ref", "S_diffs_vs_ref", "N_sites", "S_sites",
        "mean_dN_vs_ref", "mean_dS_vs_ref", "mean_gdiv_polymorphic",
        "mean_N_gdiv", "mean_S_gdiv"]


def extract(base: str, pattern: str, out_path: str, log=print) -> int:
    """返回写出的数据行数。"""
    files = sorted(glob.glob(os.path.join(
        base, pattern, "SNPGenie_Results", "product_results.txt")))
    if not files:
        raise SystemExit("[err] 没找到 product_results.txt：%s"
                         % os.path.join(base, pattern, "SNPGenie_Results"))
    os.makedirs(os.path.dirname(os.path.abspath(out_path)), exist_ok=True)
    n_files = n_rows = 0
    unknown = {}
    with open(out_path, "w", newline="") as fo:
        fo.write("sample\tgene\tproduct\t" + "\t".join(COLS) + "\n")
        for p in files:
            # <sample>/SNPGenie_Results/product_results.txt  ->  <sample>
            sample = os.path.basename(p.split("/SNPGenie_Results/")[0]).split(".")[0]
            n_files += 1
            with open(p, newline="") as fh:
                for row in csv.DictReader(fh, delimiter="\t"):
                    prod = (row.get("product") or "").strip()
                    gene = GENE_OF_PRODUCT.get(prod)
                    if gene is None:                      # 映射表没覆盖 -> 记账，不静默丢
                        unknown[prod] = unknown.get(prod, 0) + 1
                        gene = "?"
                    vals = [(row.get(c) or "").strip() for c in COLS]
                    fo.write("\t".join([sample, gene, prod] + vals) + "\n")
                    n_rows += 1
    if unknown:
        log("  [warn] 未映射的 product 编号：%s" % unknown)
    log("files=%d rows=%d -> %s" % (n_files, n_rows, out_path))
    return n_rows


def main() -> int:
    ap = argparse.ArgumentParser(
        description="汇总 SNPGenie 的 product_results.txt 成长表（供自绘 E 使用）")
    ap.add_argument("--base", default=DEFAULT_BASE, help="virus-SNPGenie 目录")
    ap.add_argument("--pattern", default=DEFAULT_PATTERN,
                    help="样本子目录的 glob，默认 %s" % DEFAULT_PATTERN)
    ap.add_argument("--out", required=True, help="输出 tsv，建议 <composite>/_data/gcva_dnds.tsv")
    args = ap.parse_args()
    extract(args.base, args.pattern, args.out, log=lambda s: print(s, file=sys.stderr))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
