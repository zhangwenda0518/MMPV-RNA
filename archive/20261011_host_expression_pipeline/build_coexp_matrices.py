#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""build_coexp_matrices.py — coexp 阶段第一步: 按病毒分组切分 TPM 矩阵
=========================================================================
把 02_Quant 的 TPM 矩阵按 04_Labels 的分组列切成"阳性/阴性"两个子矩阵,
供 TEA-GCN 分别构建网络。这是 TEA-GCN 差异共表达 (Δz) 的数据准备步。

前置闸门 (源项目 ruthenicum 教训: 两组网络样本量差 >5× 会使 zScore 基准
不可比, 产生系统性单向伪差):
  - 每组样本数 >= --min-group (默认 30)
  - 两组样本量比 <= --max-imbalance (默认 3.0)

闸门不过时:
  --on-gate-fail error (默认)  直接退出 (显式跑 coexp 时的行为)
  --on-gate-fail skip          写 group_summary.tsv 标记 SKIPPED 后退出 0
                               (编排器 --stage all 下用, 让后续阶段继续)

输出 (05_Coexp/matrices/):
  positive.tsv  negative.tsv  group_summary.tsv

用法:
  python build_coexp_matrices.py --tpm .../TPM_matrix.tsv \
      --labels .../sample_labels.tsv --output_dir .../05_Coexp/matrices \
      [--samples .../qc_passed_samples.txt] [--force]
"""

import argparse
import sys
from pathlib import Path

import pandas as pd


def gate_fail(msg: str, args, out: Path) -> None:
    """闸门失败: 按模式处理。"""
    print(f"[GATE] {msg}", file=sys.stderr)
    summary = pd.DataFrame([{
        "status": "SKIPPED", "reason": msg,
        "n_positive": "", "n_negative": "", "n_excluded": "",
        "min_group": args.min_group, "max_imbalance": args.max_imbalance,
    }])
    if args.on_gate_fail == "skip":
        out.mkdir(parents=True, exist_ok=True)
        summary.to_csv(out / "group_summary.tsv", sep="\t", index=False)
        print("[GATE] skip 模式: 已写 group_summary.tsv, 退出 0 (后续阶段继续)")
        sys.exit(0)
    sys.exit(2)


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description="按病毒分组切分 TPM 矩阵 (TEA-GCN 输入)")
    ap.add_argument("--tpm", required=True, help="TPM 矩阵 (02_Quant/matrix/TPM_matrix.tsv)")
    ap.add_argument("--labels", required=True, help="sample_labels.tsv (04_Labels)")
    ap.add_argument("--output_dir", required=True, help="输出目录 (05_Coexp/matrices)")
    ap.add_argument("--samples", help="QC 通过样本名单 (03_QC/qc_passed_samples.txt)")
    ap.add_argument("--group-col", default="Virus_any",
                    help="分组列 (默认 Virus_any; 也可用 <病毒>_positive)")
    ap.add_argument("--pos-level", default="+", help="阳性组水平 (默认 '+')")
    ap.add_argument("--neg-level", default="-", help="阴性组水平 (默认 '-')")
    ap.add_argument("--min-group", type=int, default=30,
                    help="每组样本数下限 (默认 30, goji barbarum 为 171+)")
    ap.add_argument("--max-imbalance", type=float, default=3.0,
                    help="两组样本量比上限 (默认 3.0; ruthenicum 40/331=8.3 教训)")
    ap.add_argument("--on-gate-fail", choices=["error", "skip"], default="error",
                    help="闸门不过时 error=退出 / skip=标记后正常退出")
    ap.add_argument("--force", action="store_true", help="强制重切")
    args = ap.parse_args(argv)

    out = Path(args.output_dir)
    pos_path, neg_path, sum_path = out / "positive.tsv", out / "negative.tsv", out / "group_summary.tsv"
    if not args.force and pos_path.is_file() and neg_path.is_file() and sum_path.is_file():
        print(f"切分产物已存在, 跳过: {out}")
        return

    tpm = pd.read_csv(args.tpm, sep="\t", index_col=0)
    labels = pd.read_csv(args.labels, sep="\t", dtype=str)
    if "Sample" not in labels.columns or args.group_col not in labels.columns:
        gate_fail(f"labels 缺 Sample 或 {args.group_col} 列", args, out)
    # 2026-10-10 修复: 元数据 Sample 重复时 set_index 产生重复索引, 后续
    # lab[s] 会返回 Series 而非标量, `lab[s] == level` 触发 "truth value of a
    # Series is ambiguous" 崩溃。改为显式报错并提示去重 (重复标签本身就是数据
    # 错误, 不能静默取首值)。
    dup_mask = labels["Sample"].duplicated(keep=False)
    if dup_mask.any():
        dups = sorted(labels.loc[dup_mask, "Sample"].unique())[:5]
        print(f"[ERROR] labels 的 Sample 列存在重复 "
              f"({int(dup_mask.sum())} 行, 例: {dups} ...); "
              f"重复标签会导致分组取值歧义, 请先对 labels 去重", file=sys.stderr)
        sys.exit(2)
    lab = labels.set_index("Sample")[args.group_col]

    # 样本交集: TPM 列 ∩ 标签行 ∩ QC 名单
    samples = [s for s in tpm.columns if s in lab.index]
    n_no_label = tpm.shape[1] - len(samples)
    if args.samples:
        with Path(args.samples).open(encoding="utf-8") as fh:
            passed = {l.strip() for l in fh if l.strip()}
        n_before = len(samples)
        samples = [s for s in samples if s in passed]
        print(f"QC 过滤: {n_before} → {len(samples)} (剔除 {n_before - len(samples)})")

    pos = [s for s in samples if lab[s] == args.pos_level]
    neg = [s for s in samples if lab[s] == args.neg_level]
    n_other = len(samples) - len(pos) - len(neg)
    print(f"分组 [{args.group_col}]: 阳性({args.pos_level}) {len(pos)} / "
          f"阴性({args.neg_level}) {len(neg)} / 其他 {n_other} (无标签 {n_no_label})")

    # ── 闸门 ──
    if len(pos) < args.min_group or len(neg) < args.min_group:
        gate_fail(f"组样本数不足 (阳 {len(pos)} / 阴 {len(neg)}, 需各 >= {args.min_group}); "
                  f"小样本建议改用 06_DE + 07_LoadScan, 或调低 --min-group 明知风险重跑",
                  args, out)
    ratio = max(len(pos), len(neg)) / max(1, min(len(pos), len(neg)))
    if ratio > args.max_imbalance:
        gate_fail(f"两组失衡 {ratio:.1f}× (>{args.max_imbalance}×): zScore 基准不可比, "
                  f"Δz 会系统性偏向小样本组方向 (ruthenicum 40/331 教训); "
                  f"请平衡采样或调高 --max-imbalance 明知风险重跑", args, out)

    # ── 切分输出 (保持基因行序不变) ──
    out.mkdir(parents=True, exist_ok=True)
    tpm[pos].to_csv(pos_path, sep="\t")
    tpm[neg].to_csv(neg_path, sep="\t")
    pd.DataFrame([{
        "status": "OK", "reason": "",
        "n_positive": len(pos), "n_negative": len(neg),
        "n_excluded_no_label": n_no_label, "n_excluded_other_level": n_other,
        "imbalance_ratio": round(ratio, 3),
        "min_group": args.min_group, "max_imbalance": args.max_imbalance,
        "group_col": args.group_col,
    }]).to_csv(sum_path, sep="\t", index=False)

    est_gb = 2 * 70 * tpm.shape[0] ** 2 / 1e9   # 2 网络 × (~30B RAW + ~40B Completed)/边·基因行
    print(f"切分完成: {tpm.shape[0]} 基因 × ({len(pos)}+{len(neg)}) 样本")
    print(f"  {pos_path}")
    print(f"  {neg_path}")
    if est_gb > 50:
        print(f"[WARN] 预计磁盘占用 ~{est_gb:.0f} GB (RAW+Completed 双网络), 请确认空间")


if __name__ == "__main__":
    main()


def main_testable(argv):
    """测试入口: 显式传 argv 列表, 不动 sys.argv。"""
    main(argv)
