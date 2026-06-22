#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
FastViromeExplorer Pro 结果二次过滤脚本 (增强版)
功能：关键词、覆盖度、Reads、TPM、ANI、Poisson 等多维度双向过滤 + 数据概览

用法:
  # 基础过滤
  python filter_summary.py -i all_viruses.best.summary.tsv -o filtered.tsv

  # 严格过滤 + 统计
  python filter_summary.py -i summary.tsv -o filtered.tsv \\
      --min_cov 20 --min_depth 1 --min_reads 10 --min_tpm 1 --summary report.tsv

  # 仅看统计，不写文件
  python filter_summary.py -i summary.tsv -o /dev/null --dry-run --summary report.tsv
"""

import argparse
import sys
import os

try:
    import polars as pl
except ImportError:
    print("ERROR: polars not found. Run: pip install polars")
    sys.exit(1)

# 支持的列名映射 (pipeline 输出的标准列名)
KNOWN_COLS = {
    "coverage": "Rep_Coverage(%)",
    "depth": "Rep_MeanDepth",
    "reads": "Asm_EM_Reads",
    "tpm": "Asm_TPM",
    "ani": "Avg_Read_ANI",
    "abund": "Asm_Rel_Abund(%)",
    "support": "Predicted_Support",
    "poisson": "Poisson_Ratio",
    "unique_pct": "Unique(%)",
}

# keyword 搜索的列
KEYWORD_COLS = ["Adjusted_Species", "Species_NCBI", "Species_ICTV",
                "Sample", "Rep_Accession"]


def print_summary(df, args):
    """打印数据概览"""
    n = len(df)
    print("=" * 55)
    print(f"  Data Summary  |  {args.input}")
    print("=" * 55)
    print(f"  Total records:           {n}")
    print(f"  Unique samples:          {df['Sample'].n_unique() if 'Sample' in df.columns else 'N/A'}")
    print(f"  Unique viruses (species):{df['Adjusted_Species'].n_unique() if 'Adjusted_Species' in df.columns else 'N/A'}")

    # 覆盖度/深度/Reads 分布
    for key, label in [("coverage", "Coverage(%)"), ("depth", "MeanDepth"),
                       ("reads", "EM_Reads"), ("tpm", "TPM")]:
        col = KNOWN_COLS.get(key)
        if col and col in df.columns:
            s = df[col]
            print(f"  {label:>12} — min:{s.min():.1f}  Q1:{s.quantile(0.25):.1f}  "
                  f"median:{s.median():.1f}  Q3:{s.quantile(0.75):.1f}  max:{s.max():.1f}")

    # 每样本病毒数分布
    if "Sample" in df.columns and "Adjusted_Species" in df.columns:
        smp = df.group_by("Sample").agg(pl.col("Adjusted_Species").n_unique().alias("n_virus"))
        print(f"  Viruses per sample — min:{smp['n_virus'].min()}  "
              f"median:{smp['n_virus'].median()}  max:{smp['n_virus'].max()}")

    # 每种病毒出现样本数 (Top 10)
    if "Adjusted_Species" in df.columns and "Sample" in df.columns:
        vsp = df.group_by("Adjusted_Species").agg(
            pl.col("Sample").n_unique().alias("n_samples"),
            pl.col("Rep_Coverage(%)").mean().alias("avg_cov"),
            pl.col("Rep_MeanDepth").mean().alias("avg_depth"),
            pl.col("Asm_EM_Reads").sum().alias("total_reads"),
        ).sort("n_samples", "total_reads", descending=[True, True])
        print(f"\n  Top viruses by sample count:")
        for row in vsp.head(10).iter_rows(named=True):
            print(f"    {row['Adjusted_Species'][:55]:<55} "
                  f"samples:{int(row['n_samples']):>2}  cov:{row['avg_cov']:.0f}%  "
                  f"dp:{row['avg_depth']:.0f}x  reads:{int(row['total_reads'])}")
    print("=" * 55 + "\n")


def write_summary(df, out_path):
    """输出每样本和每病毒的统计 TSV"""
    rows = []
    if "Sample" in df.columns and "Adjusted_Species" in df.columns:
        smp = df.group_by("Sample").agg([
            pl.col("Adjusted_Species").n_unique().alias("n_viruses"),
            pl.col("Asm_EM_Reads").sum().alias("total_EM_reads"),
            pl.col("Rep_Coverage(%)").mean().alias("mean_coverage"),
        ]).sort("Sample")
        smp.write_csv(str(out_path).rsplit('.', 1)[0] + ".per_sample.tsv", separator='\t')
        print(f"  Per-sample summary -> {str(out_path).rsplit('.', 1)[0]}.per_sample.tsv")

        vsp = df.group_by("Adjusted_Species").agg([
            pl.col("Sample").n_unique().alias("n_samples"),
            pl.col("Rep_Coverage(%)").mean().alias("avg_cov_%"),
            pl.col("Rep_MeanDepth").mean().alias("avg_depth"),
            pl.col("Asm_EM_Reads").sum().alias("total_EM_reads"),
            pl.col("Asm_TPM").sum().alias("total_TPM"),
        ]).sort("n_samples", "total_EM_reads", descending=[True, True])
        vsp.write_csv(str(out_path).rsplit('.', 1)[0] + ".per_virus.tsv", separator='\t')
        print(f"  Per-virus summary  -> {str(out_path).rsplit('.', 1)[0]}.per_virus.tsv")


def main():
    parser = argparse.ArgumentParser(
        description="summary.tsv 二次过滤 + 数据概览 (增强版)")

    parser.add_argument("-i", "--input", required=True, help="输入 TSV")
    parser.add_argument("-o", "--output", required=True, help="通过过滤的输出 TSV")
    parser.add_argument("-x", "--discarded", help="被过滤掉的输出 TSV (默认自动命名)")

    # keyword
    parser.add_argument("-k", "--keyword", help="关键词过滤 (跨 5 列搜索)")

    # 下限过滤
    parser.add_argument("-c", "--min_cov", type=float, default=0.0, help="最低覆盖率 %%")
    parser.add_argument("-d", "--min_depth", type=float, default=0.0, help="最低平均深度")
    parser.add_argument("-r", "--min_reads", type=float, default=0.0, help="最低 EM_Reads")
    parser.add_argument("--min_tpm", type=float, default=0.0, help="最低 TPM")
    parser.add_argument("--min_ani", type=float, default=0.0, help="最低 Avg_Read_ANI")
    parser.add_argument("--min_abund", type=float, default=0.0, help="最低 Rel_Abund(%%)")
    parser.add_argument("--min_support", type=float, default=0.0, help="最低 Predicted_Support")
    parser.add_argument("--min_poisson", type=float, default=0.0, help="最低 Poisson_Ratio")
    parser.add_argument("--min_unique", type=float, default=0.0, help="最低 Unique(%%)")

    # control
    parser.add_argument("--dry-run", action="store_true", help="仅打印统计不写文件")
    parser.add_argument("--summary", help="输出每样本/每病毒统计 TSV 前缀")

    args = parser.parse_args()

    if not os.path.exists(args.input):
        print(f"ERROR: input not found: {args.input}")
        sys.exit(1)

    discarded_file = args.discarded
    if not discarded_file:
        base, ext = os.path.splitext(args.output)
        discarded_file = f"{base}.discarded{ext}"

    print(f"Reading: {args.input}")
    try:
        df = pl.read_csv(args.input, separator='\t', ignore_errors=True)
    except Exception as e:
        print(f"ERROR reading file: {e}")
        sys.exit(1)

    n_initial = len(df)
    print(f"  Initial records: {n_initial}")

    # 数据概览 (过滤前)
    print_summary(df, args)

    # 输出统计 TSV
    if args.summary:
        write_summary(df, args.summary)

    # ---- 构建过滤掩码 ----
    mask = pl.lit(True)

    filters_applied = []

    def _col_ok(col_name):
        if col_name not in df.columns:
            print(f"  WARN: column '{col_name}' not found, skipping filter")
            return False
        return True

    # keyword (5 列搜索)
    if args.keyword:
        pattern = f"(?i){args.keyword}"
        kw_mask = pl.lit(False)
        for col in KEYWORD_COLS:
            if col in df.columns:
                kw_mask = kw_mask | pl.col(col).fill_null("").str.contains(pattern)
        if kw_mask.is_not_null().any():
            mask = mask & kw_mask
            filters_applied.append(f"keyword='{args.keyword}' (across {len([c for c in KEYWORD_COLS if c in df.columns])} cols)")

    # 数值下限过滤 (统一逻辑)
    thresholds = {
        "min_cov": ("Rep_Coverage(%)", "coverage >= {v}%"),
        "min_depth": ("Rep_MeanDepth", "depth >= {v}x"),
        "min_reads": ("Asm_EM_Reads", "reads >= {v}"),
        "min_tpm": ("Asm_TPM", "TPM >= {v}"),
        "min_ani": ("Avg_Read_ANI", "ANI >= {v}"),
        "min_abund": ("Asm_Rel_Abund(%)", "Rel_Abund >= {v}%"),
        "min_support": ("Predicted_Support", "Support >= {v}"),
        "min_poisson": ("Poisson_Ratio", "Poisson >= {v}"),
        "min_unique": ("Unique(%)", "Unique >= {v}%"),
    }

    for arg_name, (col, desc) in thresholds.items():
        val = getattr(args, arg_name)
        if val > 0:
            if _col_ok(col):
                mask = mask & (pl.col(col) >= val)
                filters_applied.append(desc.format(v=val))

    # 打印应用的过滤条件
    if filters_applied:
        print(f"\n  Filters applied:")
        for f in filters_applied:
            print(f"    - {f}")
    else:
        print("\n  No filters applied (all records pass)")

    # 拆分
    passed = df.filter(mask)
    failed = df.filter(~mask)

    print(f"\n  Result: {len(passed)} passed / {len(failed)} discarded "
          f"({n_initial - len(passed) - len(failed)} dropped as invalid)")

    if args.dry_run:
        print("\n  DRY-RUN mode: no files written.")
    else:
        if len(passed) > 0:
            passed.write_csv(args.output, separator='\t')
            print(f"  Passed -> {args.output}")
        else:
            print(f"  No records passed, output not created.")
        if len(failed) > 0:
            failed.write_csv(discarded_file, separator='\t')
            print(f"  Discarded -> {discarded_file}")


if __name__ == "__main__":
    main()
