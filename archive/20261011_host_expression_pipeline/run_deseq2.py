#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""run_deseq2.py — 差异表达分析 (⑧ 宿主表达模块, 05_DE 阶段)
=========================================================================
DESeq2 (pydeseq2) 做"病毒阳性 vs 阴性"（或任意二分组）宿主基因差异表达。

设计要点 (沿自 goji-virome 与 LSTrAP 精简工作流两版实现):
  - 输入用 estimated counts (salmon NumReads_matrix.tsv), 不用 TPM
  - 计数取整: pydeseq2 要求整数; salmon 的 estimated counts 按四舍五入转 int
  - 低表达预过滤: 全样本计数和 < --min-total-count (默认 10) 的基因剔除,
    既提速也避免 DESeq2 独立过滤的警告
  - 分组严重失衡 (>5×, 如 goji ruthenicum 的 40 vs 331) 时打印告警:
    DESeq2 本身可跑, 但结论应谨慎 (与小样本侧共表达网络分析的已知偏置呼应)

输入:
  --counts   counts 矩阵 (基因×样本; 02_Quant/matrix/NumReads_matrix.tsv)
  --labels   sample_labels.tsv (04_Labels 产物; 含分组列)
  --samples  [可选] 只用名单内样本 (如 03_QC/qc_passed_samples.txt)

输出 (默认落在 05_DE/):
  deseq2_all_results.tsv   全基因: baseMean, log2FoldChange, lfcSE, stat, pvalue, padj
  deseq2_degs.tsv          显著 DEG (padj < --padj 且 |log2FC| >= --lfc)
  deseq2_summary.tsv       对比/样本数/上下调计数
  volcano_plot.png         火山图 (matplotlib 可用时)

用法:
  python run_deseq2.py \
      --counts <08_HostExpression>/02_Quant/matrix/NumReads_matrix.tsv \
      --labels <08_HostExpression>/04_Labels/sample_labels.tsv \
      --output_dir <08_HostExpression>/05_DE
"""

import argparse
import inspect
import sys
from pathlib import Path

import numpy as np
import pandas as pd


def die(msg: str) -> None:
    print(f"[ERROR] {msg}", file=sys.stderr)
    sys.exit(1)


def load_inputs(counts_path: Path, labels_path: Path,
                samples_file: Path = None):
    """读取 counts 与样本标签, 取交集样本并校验分组。"""
    counts = pd.read_csv(counts_path, sep="\t", index_col=0)
    labels = pd.read_csv(labels_path, sep="\t", dtype=str)
    if "Sample" not in labels.columns:
        die(f"{labels_path} 缺 Sample 列; 请用 04_Labels 阶段产物")

    keep = None
    if samples_file is not None:
        with samples_file.open(encoding="utf-8") as fh:
            keep = {l.strip() for l in fh if l.strip()}
        missing = keep - set(counts.columns)
        if missing:
            print(f"[WARN] 名单中 {len(missing)} 个样本无定量结果, 已忽略: "
                  f"{sorted(missing)[:5]}{'...' if len(missing) > 5 else ''}",
                  file=sys.stderr)

    lab = labels.set_index("Sample")
    sample_ids = [c for c in counts.columns
                  if c in lab.index and (keep is None or c in keep)]
    dropped = len(counts.columns) - len(sample_ids)
    if dropped:
        print(f"[WARN] {dropped} 个定量样本无标签或不在名单内, 已剔除",
              file=sys.stderr)
    if not sample_ids:
        die("counts 与 labels 无交集样本")
    return counts[sample_ids], lab.loc[sample_ids]


def prefilter_counts(counts: pd.DataFrame, min_total: float) -> pd.DataFrame:
    """剔除低表达基因 (全样本计数和 < min_total)。"""
    totals = counts.sum(axis=1)
    kept = counts.loc[totals >= min_total]
    print(f"预过滤: 基因 {counts.shape[0]} → {kept.shape[0]} "
          f"(全样本计数和 >= {min_total:g})")
    return kept


def run_deseq(counts: pd.DataFrame, meta: pd.DataFrame,
              group_col: str, ref_level: str, alt_level: str,
              threads: int):
    """跑 pydeseq2 (兼容 0.4 的 design_factors 与 0.5+ 的 design 两代 API)。

    2026-10-10 修复: 旧实现用 `try: design_factors=… except TypeError: design=…`
    兜底 —— 任何参数错误 (拼写/缺参/类型不符) 都会被误判为 API 版本差异而静默
    重试, 掩盖根因。现改为**检查 DeseqDataSet 构造签名**决定用哪一代参数, 只在
    签名确实不含目标形参时才回退, 真实 TypeError 不再被吞。
    """
    try:
        from pydeseq2.dds import DeseqDataSet
        from pydeseq2.ds import DeseqStats
    except ImportError:
        die("未安装 pydeseq2。安装: pip install pydeseq2\n"
            "  (或 conda install -c conda-forge pydeseq2)")

    counts_int = counts.round().astype(int).T          # pydeseq2: 样本×基因
    metadata = meta[[group_col]].copy()
    metadata[group_col] = metadata[group_col].astype("category")

    kwargs = dict(counts=counts_int, metadata=metadata,
                  ref_level=[group_col, ref_level], n_cpus=threads)
    try:
        dds_params = inspect.signature(DeseqDataSet.__init__).parameters
    except (TypeError, ValueError):                    # C 扩展/无法内省 → 兜底
        dds_params = {}
    if "design_factors" in dds_params:                 # pydeseq2 <=0.4
        dds = DeseqDataSet(design_factors=group_col, **kwargs)
    elif "design" in dds_params:                       # pydeseq2 >=0.5
        dds = DeseqDataSet(design=f"~{group_col}", **kwargs)
    else:
        die("pydeseq2 的 DeseqDataSet 既不接受 design_factors 也不接受 design; "
            f"实际形参: {sorted(dds_params) or '无法内省'} —— 版本不兼容")
    dds.deseq2()

    stat = DeseqStats(dds, contrast=[group_col, alt_level, ref_level],
                      n_cpus=threads)
    stat.summary()
    return stat.results_df


def plot_volcano(res: pd.DataFrame, padj_cut: float, lfc_cut: float,
                 out_png: Path, title: str) -> bool:
    """火山图 (matplotlib 缺失时跳过, 不影响主产物)。"""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("[WARN] 未安装 matplotlib, 跳过火山图")
        return False

    sub = res[["log2FoldChange", "padj"]].dropna()
    neglogp = -np.log10(sub["padj"].clip(lower=1e-300))
    sig = (sub["padj"] < padj_cut) & (sub["log2FoldChange"].abs() >= lfc_cut)
    fig, ax = plt.subplots(figsize=(8, 7))
    ax.scatter(sub.loc[~sig, "log2FoldChange"], neglogp[~sig],
               c="#9e9e9e", s=6, alpha=0.5, label=f"NS ({int((~sig).sum())})")
    ax.scatter(sub.loc[sig, "log2FoldChange"], neglogp[sig],
               c="#d62728", s=7, alpha=0.7, label=f"DEG ({int(sig.sum())})")
    ax.axvline(lfc_cut, color="grey", ls="--", lw=0.8)
    ax.axvline(-lfc_cut, color="grey", ls="--", lw=0.8)
    ax.set_xlabel("log2 fold change")
    ax.set_ylabel("-log10(padj)")
    ax.set_title(title)
    ax.legend(loc="upper right", frameon=False)
    fig.tight_layout()
    fig.savefig(out_png, dpi=300)
    plt.close(fig)
    return True


def main() -> None:
    ap = argparse.ArgumentParser(description="pydeseq2 差异表达 (病毒阳性 vs 阴性)")
    ap.add_argument("--counts", required=True,
                    help="counts 矩阵 (基因×样本; NumReads_matrix.tsv)")
    ap.add_argument("--labels", required=True,
                    help="sample_labels.tsv (04_Labels 产物)")
    ap.add_argument("--output_dir", required=True, help="输出目录 (05_DE)")
    ap.add_argument("--group-col", default="Virus_any",
                    help="分组列 (默认 Virus_any; 也可用 <病毒>_positive)")
    ap.add_argument("--ref-level", default="-", help="对照组水平 (默认 '-')")
    ap.add_argument("--alt-level", default="+", help="处理组水平 (默认 '+')")
    ap.add_argument("--samples", help="只用名单内样本 (如 03_QC/qc_passed_samples.txt)")
    ap.add_argument("--padj", type=float, default=0.05, help="padj 阈值 (默认 0.05)")
    ap.add_argument("--lfc", type=float, default=1.0, help="|log2FC| 阈值 (默认 1.0)")
    ap.add_argument("--min-total-count", type=float, default=10.0,
                    help="低表达预过滤 (全样本计数和阈值, 默认 10)")
    ap.add_argument("--strata", default=None,
                    help="分层列 (如 Tissue): 每层内各做一次对比, 结果写 "
                         "deseq2_strata/<层>.tsv; 公共队列建议开启 "
                         "(组织构成混杂防护, goji ERAD 教训)")
    ap.add_argument("--threads", "-t", type=int, default=8)
    args = ap.parse_args()

    counts, lab = load_inputs(Path(args.counts), Path(args.labels),
                              Path(args.samples) if args.samples else None)
    if args.group_col not in lab.columns:
        die(f"分组列 {args.group_col!r} 不在 labels 中; 可选: "
            f"{[c for c in lab.columns if c not in ('Sample',)][:10]}...")
    groups = lab[args.group_col].value_counts()
    for lvl in (args.ref_level, args.alt_level):
        if groups.get(lvl, 0) < 2:
            die(f"分组 {lvl!r} 只有 {groups.get(lvl, 0)} 个样本, DESeq2 至少需 2 个\n"
                f"当前分组分布: {dict(groups)}")
    ratio = groups.max() / groups.min()
    if ratio > 5:
        print(f"[WARN] 分组严重失衡 ({groups.max()} vs {groups.min()}, "
              f"{ratio:.1f}×) —— DESeq2 仍可运行, 但结论谨慎解释",
              file=sys.stderr)
    print(f"分组: {dict(groups)}")

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    contrast = f"{args.group_col}: {args.alt_level} vs {args.ref_level}"
    summary_rows = []

    def emit(tag: str, sub_counts: pd.DataFrame, sub_lab: pd.DataFrame,
             is_global: bool = False) -> None:
        """跑一次对比并落盘 (tag = all 或层名)。

        2026-10-10 修复: 旧实现用 `if tag == "all"` 区分全局/分层 —— 若分层列
        恰好有个层名叫 "all", 该层会覆盖全局结果文件 deseq2_all_results.tsv。
        现改用显式 is_global 标志, 层名与全局解耦。
        """
        sub_counts = prefilter_counts(sub_counts, args.min_total_count)
        res = run_deseq(sub_counts, sub_lab, args.group_col,
                        args.ref_level, args.alt_level, args.threads)
        res.index.name = "gene"
        sig = res[(res["padj"].notna()) & (res["padj"] < args.padj)
                  & (res["log2FoldChange"].abs() >= args.lfc)]
        degs = sig.sort_values("padj")
        if is_global:
            res.to_csv(out_dir / "deseq2_all_results.tsv", sep="\t")
            degs.to_csv(out_dir / "deseq2_degs.tsv", sep="\t")
            plot_volcano(res, args.padj, args.lfc, out_dir / "volcano_plot.png",
                         f"{contrast}  (DEG: {len(degs)})" if len(degs) else
                         f"{contrast}  (no DEG at thresholds)")
        else:
            strata_dir = out_dir / "deseq2_strata"
            strata_dir.mkdir(exist_ok=True)
            res.to_csv(strata_dir / f"{tag}.tsv", sep="\t")
        sub_groups = sub_lab[args.group_col].value_counts()
        n_up = int((degs["log2FoldChange"] > 0).sum())
        n_down = int((degs["log2FoldChange"] < 0).sum())
        summary_rows.append({
            "contrast": contrast, "stratum": tag,
            "n_samples": int(len(sub_lab)),
            f"n_{args.ref_level}": int(sub_groups.get(args.ref_level, 0)),
            f"n_{args.alt_level}": int(sub_groups.get(args.alt_level, 0)),
            "n_genes_tested": int(res["padj"].notna().sum()),
            "padj_cutoff": args.padj, "lfc_cutoff": args.lfc,
            "n_deg": int(len(degs)), "n_up": n_up, "n_down": n_down,
        })
        print(f"[{tag}] DEG {len(degs)} (上调 {n_up} / 下调 {n_down})")

    # 汇总 (全部有标签样本)
    emit("all", counts, lab, is_global=True)
    # 分层 (每层内独立对比; 至少两组各 2 样本)
    if args.strata:
        if args.strata not in lab.columns:
            print(f"[WARN] 分层列 {args.strata!r} 不在 labels 中, 跳过分层",
                  file=sys.stderr)
        else:
            strata = lab[args.strata].fillna("unknown").astype(str).str.strip()
            for tag, idx in sorted(strata.groupby(strata).groups.items()):
                sub_lab = lab.loc[idx]
                sub_groups = sub_lab[args.group_col].value_counts()
                if min(sub_groups.get(args.ref_level, 0),
                       sub_groups.get(args.alt_level, 0)) < 2:
                    print(f"[strata {tag}] 两组各需 >=2 样本 "
                          f"({dict(sub_groups)}), 跳过")
                    continue
                print(f"── 分层 {tag}: {dict(sub_groups)} ──")
                emit(str(tag), counts[idx], sub_lab)

    pd.DataFrame(summary_rows).to_csv(out_dir / "deseq2_summary.tsv",
                                      sep="\t", index=False)
    all_row = summary_rows[0]
    print(f"全基因结果: {out_dir / 'deseq2_all_results.tsv'} "
          f"({all_row['n_genes_tested']} 基因受检)")
    print(f"显著 DEG  : {out_dir / 'deseq2_degs.tsv'} "
          f"({all_row['n_deg']} 个: 上调 {all_row['n_up']} / 下调 {all_row['n_down']}; "
          f"padj<{args.padj}, |log2FC|>={args.lfc})")
    if args.strata:
        print("注意: 汇总(all)结果在异质公共队列上受组织构成混杂影响 "
              "(goji ERAD 教训), 结论以分层结果与 05_Coexp/07_LoadScan 交叉验证为准")


if __name__ == "__main__":
    main()
