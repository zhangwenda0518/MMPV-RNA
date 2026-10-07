#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""diff_coexp.py — TEA-GCN 双网络差异共表达 ΔzScore(MR) (coexp 阶段最终步)
=========================================================================
比较"阳性网络"与"阴性网络"两个 Completed_GCN 的 zScore(Co-exp_Str_MR):
Δz = z(阳性) − z(阴性)。每基因取 |Δz| 最大的 Top-N 邻居为候选对, 用候选对
Δz 采样分布的 0.5%/99.5% 分位阈值判显著。

重写自 goji 的 diff_coexp_strict.py + _post.py, 改进:
  - Target 行集按名对齐 (源版按位置对齐; 两网络 selected k 不同导致 NaN
    过滤模式不同时, 行集天然不一致 — 取交集, 单侧缺失边不计入 Δz)
  - 自边按 Target==基因名 精确剔除 (源版固定 dz[0]=0, 只对第 1 个基因正确)
  - Target 输出真实基因名 (源版 Target_idx 0 基索引, 曾致 7,547/3,908 口径误用)
  - Top-N 与采样步长参数化 (源版硬编码 20 与 ::392)

输入: 两个 Completed_GCN 目录 (每基因一文件, 4 列; 第 0/3 列 = Target/zScore)
输出 (05_Coexp/):
  diff_coexp_candidates.tsv    候选对 (每基因 Top-N, 有向行)
  diff_coexp_significant.tsv   显著对 (分位阈值过滤)
  diff_coexp_gene_aggregate.tsv 每源基因的显著对计数/平均 Δz (hub 排序)
  diff_coexp_summary.tsv       阈值与计数 (有向行 + 唯一对两种口径)
"""

import argparse
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd


def list_gene_files(gcn_dir: Path):
    return sorted((p for p in gcn_dir.iterdir() if p.is_file()),
                  key=lambda p: p.name)


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description="TEA-GCN 双网络 ΔzScore(MR) 差异共表达")
    ap.add_argument("--pos-gcn", required=True, help="阳性组 Completed_GCN 目录")
    ap.add_argument("--neg-gcn", required=True, help="阴性组 Completed_GCN 目录")
    ap.add_argument("--output_dir", required=True, help="输出目录 (05_Coexp)")
    ap.add_argument("--top-n", type=int, default=20, help="每基因取 |Δz| Top-N 邻居 (默认 20)")
    ap.add_argument("--q-low", type=float, default=0.5, help="显著阈值下分位 (默认 0.5)")
    ap.add_argument("--q-high", type=float, default=99.5, help="显著阈值上分位 (默认 99.5)")
    args = ap.parse_args(argv)

    pos_dir, neg_dir = Path(args.pos_gcn), Path(args.neg_gcn)
    if not pos_dir.is_dir() or not neg_dir.is_dir():
        sys.exit(f"[ERROR] Completed_GCN 目录不存在: {pos_dir} / {neg_dir}")
    pos_files = list_gene_files(pos_dir)
    if not pos_files:
        sys.exit(f"[ERROR] 阳性网络目录为空: {pos_dir}")
    neg_names = {p.name for p in list_gene_files(neg_dir)}
    pos_names = {p.name for p in pos_files}
    if pos_names != neg_names:
        only_pos = sorted(pos_names - neg_names)[:3]
        only_neg = sorted(neg_names - pos_names)[:3]
        sys.exit(f"[ERROR] 两网络基因文件集不一致 "
                 f"(pos {len(pos_names)} vs neg {len(neg_names)}); "
                 f"仅阳有 {only_pos}... 仅阴有 {only_neg}... — "
                 f"两矩阵必须同基因集 (切分只动列)")
    genes = [p.name for p in pos_files]
    n = len(genes)
    print(f"双网络: {n} 基因; 每基因 Top-{args.top_n} |Δz| 候选")

    step = max(1, n // 100)          # 每基因采样 ~100 个 Δz 构建分布
    rows, deltas = [], []
    align_warned = False
    t0 = time.time()
    for i, g in enumerate(genes):
        zp = pd.read_csv(pos_dir / g, sep="\t", usecols=[0, 3])
        zn = pd.read_csv(neg_dir / g, sep="\t", usecols=[0, 3])
        tp = zp.iloc[:, 0].astype(str).to_numpy()
        tn = zn.iloc[:, 0].astype(str).to_numpy()
        zv = zp.iloc[:, 1].to_numpy(dtype=float)
        znv = zn.iloc[:, 1].to_numpy(dtype=float)
        if len(tp) == len(tn) and bool((tp == tn).all()):
            targets, za, zb = tp, zv, znv          # 快路径: 行序逐行一致
        else:
            # 慢路径: 两网络 (常因 selected k 不同 → NaN 过滤模式不同) 行集
            # 不一致 — 按 Target 名取交集对齐, 单侧缺失的边不计入 Δz
            if not align_warned:
                print(f"[对齐] 两网络行集不同 (pos {len(tp)} vs neg {len(tn)}); "
                      f"按 Target 名交集对齐, 单侧缺失边不计入", flush=True)
                align_warned = True
            pidx = {t: k for k, t in enumerate(tp)}
            keep = np.fromiter((t in pidx for t in tn), dtype=bool,
                               count=len(tn))
            targets = tn[keep]
            za = np.fromiter((zv[pidx[t]] for t in targets), dtype=float,
                             count=len(targets))
            zb = znv[keep]
        dz = za - zb
        dz[targets == g] = np.nan          # 自边精确剔除
        valid = dz[~np.isnan(dz)]
        if valid.size:
            deltas.extend(valid[::step].tolist())
        score = np.abs(np.nan_to_num(dz, nan=0.0))
        for j in np.argsort(-score)[:args.top_n]:
            if np.isnan(dz[j]) or abs(dz[j]) < 1e-9:
                continue
            rows.append((g, targets[j], za[j], zb[j], dz[j]))
        if (i + 1) % 2000 == 0:
            print(f"  {i + 1}/{n}  {time.time() - t0:.0f}s", flush=True)

    cand = pd.DataFrame(rows, columns=["Gene", "Target", "z_pos", "z_neg", "Delta_z"])
    deltas = np.asarray(deltas, dtype=float)
    deltas = deltas[~np.isnan(deltas)]
    if deltas.size == 0:
        # 网络退化 (两网络 zScore 全 NaN / 全被自边剔除) → 分位阈值无从计算。
        # 2026-10-03: 旧行为在 np.percentile(空数组) 处抛 IndexError, 信息量低;
        # 现给出明确错误并以退出码 3 结束 (与 coexp_network_qc 的 QC-FAIL 约定一致)。
        print("[ERROR] Δz 采样分布为空: 两网络 zScore 全为 NaN (网络退化) 或所有 "
              "Δz 均被剔除, 无法计算分位阈值 → 退出码 3", file=sys.stderr)
        sys.exit(3)
    lo, hi = (np.percentile(deltas, args.q_low),
              np.percentile(deltas, args.q_high))
    print(f"Δz 分布 (采样 {len(deltas)}): mean={deltas.mean():.4f} sd={deltas.std():.4f} "
          f"p{args.q_low:g}={lo:.4f} p{args.q_high:g}={hi:.4f}")

    sig = cand[(cand["Delta_z"] > hi) | (cand["Delta_z"] < lo)].copy()
    sig = sig.dropna(subset=["Delta_z"]).sort_values("Delta_z", ascending=False)

    pair = sig.apply(lambda r: tuple(sorted((r["Gene"], r["Target"]))), axis=1)
    n_uniq = int(pair.nunique())
    n_enh = int((sig["Delta_z"] > 0).sum())
    n_wea = int((sig["Delta_z"] < 0).sum())
    print(f"显著对: {len(sig)} 有向行 = {n_uniq} 唯一对 "
          f"(增强 {n_enh} / 减弱 {n_wea}; 涉及基因 "
          f"{pd.unique(pd.concat([sig['Gene'], sig['Target']])).size})")

    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    cand.to_csv(out / "diff_coexp_candidates.tsv", sep="\t", index=False)
    sig.to_csv(out / "diff_coexp_significant.tsv", sep="\t", index=False)

    agg = (sig.groupby("Gene", as_index=False)
              .agg(n_sig_pairs=("Delta_z", "size"),
                   mean_delta_z=("Delta_z", "mean"),
                   max_abs_delta_z=("Delta_z", lambda x: float(np.abs(x).max())))
              .sort_values(["n_sig_pairs", "max_abs_delta_z"], ascending=False))
    agg.to_csv(out / "diff_coexp_gene_aggregate.tsv", sep="\t", index=False)

    pd.DataFrame([{
        "n_genes": n, "top_n": args.top_n,
        "q_low": args.q_low, "q_high": args.q_high,
        "threshold_low": round(float(lo), 4), "threshold_high": round(float(hi), 4),
        "n_significant_rows": len(sig), "n_unique_pairs": n_uniq,
        "n_enhanced": n_enh, "n_weakened": n_wea,
        "n_genes_involved": int(pd.unique(pd.concat([sig["Gene"], sig["Target"]])).size),
    }]).to_csv(out / "diff_coexp_summary.tsv", sep="\t", index=False)

    for f in ("diff_coexp_candidates.tsv", "diff_coexp_significant.tsv",
              "diff_coexp_gene_aggregate.tsv", "diff_coexp_summary.tsv"):
        print(f"输出: {out / f}")


if __name__ == "__main__":
    main()


def main_testable(argv):
    """测试入口: 显式传 argv 列表, 不动 sys.argv。"""
    main(argv)
