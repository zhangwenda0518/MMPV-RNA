#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
📋 SRA Run 列表提取器 v1.0
===========================
从 GSA/SRA 合并 CSV 中快速提取 Run ID 列表，
支持按数据库、组织、策略等字段筛选。

用法:
  python extract_sra_list.py --input SRA_GSA_Merged_Final.csv --output sra.list
  python extract_sra_list.py --input merged.csv --output sra.list --filter-db SRA
  python extract_sra_list.py --input merged.csv --output sra.list --filter-col Tissue --filter-val leaf
"""

import os
import sys
import argparse
import pandas as pd

if sys.platform == 'win32':
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')


def main():
    parser = argparse.ArgumentParser(
        description="📋 从 merged CSV 提取 SRA Run 列表"
    )
    parser.add_argument('--input', '-i', required=True,
                        help='输入的 SRA_GSA_Merged_Final.csv 路径')
    parser.add_argument('--output', '-o', required=True,
                        help='输出的 Run 列表文件 (.list 或 .txt)')
    parser.add_argument('--filter-db', choices=['SRA', 'GSA'],
                        help='仅提取指定数据库的 Run')
    parser.add_argument('--filter-col', help='按列名过滤')
    parser.add_argument('--filter-val', help='按列值过滤 (与 --filter-col 配合使用)')
    parser.add_argument('--head', type=int, help='仅输出前 N 个 Run')

    args = parser.parse_args()

    if not os.path.exists(args.input):
        print(f"❌ 找不到输入文件: {args.input}")
        return 1

    df = pd.read_csv(args.input)

    # 查找 Run 列
    run_col = None
    for col in df.columns:
        if col.strip().lower() == 'run':
            run_col = col
            break

    if run_col is None:
        print(f"❌ 未找到 'Run' 列！可用列: {list(df.columns)}")
        return 1

    # 过滤
    if args.filter_db and 'Database' in df.columns:
        df = df[df['Database'] == args.filter_db]
        print(f"  过滤数据库: {args.filter_db} → {len(df)} 条")

    if args.filter_col and args.filter_val:
        if args.filter_col in df.columns:
            mask = df[args.filter_col].astype(str).str.contains(
                args.filter_val, case=False, na=False
            )
            df = df[mask]
            print(f"  过滤 {args.filter_col}={args.filter_val} → {len(df)} 条")
        else:
            print(f"⚠ 列 '{args.filter_col}' 不存在，跳过过滤")

    # 提取 Run 列表
    runs = df[run_col].dropna().astype(str).str.strip()
    runs = runs[runs != ''].tolist()

    if args.head:
        runs = runs[:args.head]

    # 写入文件
    os.makedirs(os.path.dirname(args.output) or '.', exist_ok=True)
    with open(args.output, 'w', encoding='utf-8') as f:
        f.write('\n'.join(runs))

    print(f"✅ 成功提取 {len(runs)} 个 Run ID → {args.output}")
    return 0


if __name__ == '__main__':
    exit(main())
