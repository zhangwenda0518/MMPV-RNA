#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""diff_coexp.py — TEA-GCN 双网络差异共表达 ΔzScore(MR) (coexp 阶段最终步)
=========================================================================
比较"阳性网络"与"阴性网络"两个 Completed_GCN 的 zScore(Co-exp_Str_MR):
Δz = z(阳性) − z(阴性)。每基因取 |Δz| 最大的 Top-N 邻居为候选对, 用分位阈值
筛出尾部候选对。

阈值口径 (2026-10-10 修复): 候选对是"每基因 |Δz| Top-N"的**极值子集**, 而
全总体 Δz 采样近似"随机一条边"的零分布。把总体分位线施加于极值子集天然偏松
(显著对数随阈值口径漂移)。本脚本据此提供 --threshold-basis:
  population (默认)  总体分位线 (与历史结果一致)
  candidates         候选对自洽切 (等价符号翻转置换零模型: 零假设下 Δz 关于 0
                     对称, 每基因 Top-N 选择不变、仅符号翻转, 故候选集分位即
                     零分布分位 → 固定为候选的前/后 q%)
无论哪种口径, 本脚本输出的都是**描述性尾部候选对**, 不是校准统计检验; 真正
p 值须在样本标签层面置换并整网重建 (本脚本不做)。脚本会打印"零假设下预期命中
数"(候选数 × 尾部概率) 以量化选择偏倚。

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
    ap.add_argument("--threshold-basis", choices=["population", "candidates"],
                    default="population",
                    help="阈值分位取自哪一分布: population=全总体 Δz 采样 "
                         "(默认, 与历史一致) / candidates=候选对自洽切")
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
    if cand.empty or deltas.size == 0:
        # 网络退化 (两网络 zScore 全 NaN / 全被自边剔除) → 无候选对或 Δz 采样
        # 分布为空, 分位阈值无从计算。2026-10-03: 旧行为在 np.percentile(空
        # 数组) 处抛 IndexError, 信息量低; 现给出明确错误并以退出码 3 结束
        # (与 coexp_network_qc 的 QC-FAIL 约定一致)。
        print("[ERROR] 候选对/Δz 采样分布为空: 两网络 zScore 全为 NaN (网络退化) "
              "或所有 Δz 均被剔除, 无法计算分位阈值 → 退出码 3", file=sys.stderr)
        sys.exit(3)

    # ── 显著阈值 (口径见文件头) ──
    if args.threshold_basis == "candidates":
        ref = cand["Delta_z"].to_numpy(dtype=float)
        lo, hi = (float(np.percentile(ref, args.q_low)),
                  float(np.percentile(ref, args.q_high)))
        # 自洽切: 含端点, 否则极值本身落在分位点上会被 > / < 排除而全空
        sel = (cand["Delta_z"] >= hi) | (cand["Delta_z"] <= lo)
        print(f"阈值基准=candidates (候选 n={ref.size}): "
              f"p{args.q_low:g}={lo:.4f} p{args.q_high:g}={hi:.4f}")
    else:
        lo, hi = (float(np.percentile(deltas, args.q_low)),
                  float(np.percentile(deltas, args.q_high)))
        sel = (cand["Delta_z"] > hi) | (cand["Delta_z"] < lo)
        print(f"Δz 总体采样 (n={deltas.size}): mean={deltas.mean():.4f} "
              f"sd={deltas.std():.4f}  p{args.q_low:g}={lo:.4f} "
              f"p{args.q_high:g}={hi:.4f}")

    # 选择偏倚诊断: 零假设下, 候选对按尾部概率随机命中多少行
    tail_prob = (args.q_low + (100.0 - args.q_high)) / 100.0
    n_cand = int(len(cand))
    exp_chance = n_cand * tail_prob

    sig = cand[sel].copy()
    sig = sig.dropna(subset=["Delta_z"]).sort_values("Delta_z", ascending=False)
    print(f"零假设下预期命中 ~{exp_chance:.1f} 行 "
          f"(候选 {n_cand} × 尾部概率 {tail_prob:.3g}); 实测 {len(sig)} 行")

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
        "threshold_basis": args.threshold_basis,
        "q_low": args.q_low, "q_high": args.q_high,
        "threshold_low": round(float(lo), 4), "threshold_high": round(float(hi), 4),
        "n_candidates": n_cand,
        "expected_sig_by_chance": round(float(exp_chance), 2),
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
