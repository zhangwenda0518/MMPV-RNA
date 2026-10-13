#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""viral_load_scan.py — 病毒载量×全基因组相关扫描 (⑧ 宿主表达模块, 08_LoadScan 阶段)
=========================================================================
把"病毒阳性/阴性"二分组升级为连续剂量变量: 逐病毒, 对每个宿主基因做
载量 × 基因表达 (TPM) 的 Spearman 相关 + BH-FDR 多重校正。

范式参考: Torres-Rodriguez et al. 2024 Plant J (群体表达-性状关联);
goji-virome 用同一设计在果实组织内找到 12,475 个 GCVA 载量响应基因。

实现说明: Spearman 用"秩变换 + 向量化 Pearson"精确实现 (与
scipy.stats.spearmanr 的 rho/p 逐值一致, 但整矩阵一次算完), 无新依赖。

输入:
  --tpm       TPM 矩阵 (基因×样本; 02_Quant/matrix/TPM_matrix.tsv)
  --labels    virus_labels.tsv (04_Labels 长表; 载量来自 Load 列)
  --meta      [可选] sample_labels.tsv —— 提供 Tissue 列供 --tissue 过滤

输出 (默认落在 06_LoadScan/):
  <病毒>_genome_scan.tsv   逐病毒: gene, virus, n_samples, rho, pvalue, padj
  load_scan_summary.tsv    汇总: virus, n_samples, n_positives, n_genes_tested,
                           n_sig_padj05, n_sig_positive_direction

用法:
  python viral_load_scan.py \
      --tpm <08_HostExpression>/02_Quant/matrix/TPM_matrix.tsv \
      --labels <08_HostExpression>/04_Labels/virus_labels.tsv \
      --meta <08_HostExpression>/04_Labels/sample_labels.tsv \
      --tissue fruit,berry \
      --output_dir <08_HostExpression>/06_LoadScan
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import rankdata

from derive_virus_labels import sanitize_virus


def bh_fdr(pvals) -> np.ndarray:
    """Benjamini-Hochberg step-up FDR (NaN 安全; 与 R p.adjust('BH') 一致)。"""
    p = np.asarray(pvals, dtype=float)
    out = np.full(p.shape, np.nan)
    idx = np.where(~np.isnan(p))[0]
    if idx.size == 0:
        return out
    order = idx[np.argsort(p[idx], kind="mergesort")]
    ranked = p[order]
    m = order.size
    q = ranked * m / np.arange(1, m + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]          # step-up 单调化
    out[order] = np.clip(q, 0.0, 1.0)
    return out


def spearman_scan(load: np.ndarray, X: np.ndarray):
    """基因×样本矩阵 X 对载量 load 的 Spearman rho 与 p (向量化)。

    返回 (rho, pvalue); 零方差基因 (秩全等) 记 NaN。
    p 用 t 近似: t = rho*sqrt((n-2)/(1-rho^2)), 与 scipy.stats.spearmanr 同式。
    """
    n = len(load)
    lr = rankdata(load).astype(float)
    lr -= lr.mean()
    Xr = rankdata(X, axis=1).astype(float)
    Xr -= Xr.mean(axis=1, keepdims=True)

    num = Xr @ lr
    den = np.sqrt((Xr ** 2).sum(axis=1) * (lr ** 2).sum())
    with np.errstate(invalid="ignore", divide="ignore"):
        rho = np.where(den > 0, num / den, np.nan)
    rho = np.clip(rho, -1.0, 1.0)

    with np.errstate(invalid="ignore", divide="ignore"):
        t = rho * np.sqrt((n - 2) / (1 - rho ** 2))
    p = 2 * stats.t.sf(np.abs(t), df=n - 2)
    return rho, p


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description="病毒载量×全基因组 Spearman 扫描 + BH-FDR")
    ap.add_argument("--tpm", required=True, help="TPM 矩阵 (基因×样本)")
    ap.add_argument("--labels", required=True, help="virus_labels.tsv (长表)")
    ap.add_argument("--output_dir", required=True, help="输出目录 (06_LoadScan)")
    ap.add_argument("--meta", help="sample_labels.tsv (提供 Tissue 供 --tissue 过滤)")
    ap.add_argument("--tissue", help="范围过滤 (按 --strata-col 列子串匹配, 逗号分隔)")
    ap.add_argument("--strata-col", default="Tissue",
                    help="组内扫描的分层列名 (默认 Tissue; 可换任意元数据列)")
    ap.add_argument("--min-positives", type=int, default=5,
                    help="该病毒阳性样本数下限, 低于则跳过 (默认 5)")
    ap.add_argument("--samples",
                    help="QC 通过样本名单 (03_QC/qc_passed_samples.txt); "
                         "给定则扫描仅限名单内样本")
    ap.add_argument("--min-samples", type=int, default=8,
                    help="参与相关分析的样本数下限 (默认 8)")
    args = ap.parse_args(argv)

    tpm = pd.read_csv(args.tpm, sep="\t", index_col=0)
    vlabels = pd.read_csv(args.labels, sep="\t")

    sample_ids = list(tpm.columns)
    # QC 过滤 (剔除 03_QC 判 fail 的低比对/低 reads 样本)
    if args.samples:
        with Path(args.samples).open(encoding="utf-8") as fh:
            passed = {l.strip() for l in fh if l.strip()}
        kept = [s for s in sample_ids if s in passed]
        print(f"QC 过滤: {len(sample_ids)} → {len(kept)} "
              f"(剔除 {len(sample_ids) - len(kept)})")
        sample_ids = kept
    # 范围过滤 (基于 sample_labels 的分层列, --strata-col)
    if args.tissue:
        if not args.meta:
            print("[WARN] --tissue 需配合 --meta (sample_labels.tsv), 已忽略过滤",
                  file=sys.stderr)
        else:
            meta = pd.read_csv(args.meta, sep="\t", dtype=str)
            if args.strata_col not in meta.columns:
                print(f"[WARN] sample_labels 无 {args.strata_col} 列, 已忽略 --tissue",
                      file=sys.stderr)
            else:
                pats = [p.strip().lower() for p in args.tissue.split(",") if p.strip()]
                meta["_t"] = meta[args.strata_col].fillna("").str.lower()
                keep = set(meta.loc[meta["_t"].str.contains("|".join(pats),
                                                           regex=True, na=False),
                                    "Sample"])
                sample_ids = [s for s in sample_ids if s in keep]
                print(f"组织过滤 [{args.tissue}]: 样本 {len(tpm.columns)} → "
                      f"{len(sample_ids)}")

    if len(sample_ids) < args.min_samples:
        print(f"[ERROR] 参与分析样本数 {len(sample_ids)} < {args.min_samples}; "
              f"放宽过滤或调低 --min-samples", file=sys.stderr)
        sys.exit(1)

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    genes = tpm.index.to_numpy()

    def scan_scope(stratum: str, scope_ids: list) -> list:
        """对一个样本范围做逐病毒扫描; 追加带 stratum 列的 summary 行。"""
        rows = []

        def _load_of(load_map: dict, s: str) -> float:
            """取样本载量; 缺行/None/NaN 一律记 0.0。

            2026-10-10 修复: 旧式 `float(load_map.get(s, 0.0) or 0.0)` 对 NaN
            失效 —— NaN 是真值, `nan or 0.0` 返回 nan, 会混入 Spearman 并被
            rankdata 排到最大秩, 污染 rho/p。
            """
            v = load_map.get(s, 0.0)
            if v is None:
                return 0.0
            try:
                v = float(v)
            except (TypeError, ValueError):
                return 0.0
            return v if np.isfinite(v) else 0.0

        def _pos_of(pos_map: dict, s: str) -> int:
            v = pos_map.get(s, 0)
            if v is None or (isinstance(v, float) and not np.isfinite(v)):
                return 0
            try:
                return int(v)
            except (TypeError, ValueError):
                return 0

        X = tpm[scope_ids].to_numpy(dtype=float)
        for virus, g in vlabels.groupby("Virus"):
            load_map = dict(zip(g["Sample"], g["Load"]))
            pos_map = dict(zip(g["Sample"], g["Positive"]))
            load = np.array([_load_of(load_map, s) for s in scope_ids],
                            dtype=float)
            n_pos = int(sum(_pos_of(pos_map, s) for s in scope_ids))

            if n_pos < args.min_positives:
                print(f"[跳过] [{stratum}] {virus}: 阳性 {n_pos} < {args.min_positives}")
                continue
            if np.unique(load).size < 2:
                print(f"[跳过] [{stratum}] {virus}: 载量无变异")
                continue

            rho, p = spearman_scan(load, X)
            padj = bh_fdr(p)
            n_sig = int(np.nansum(padj < 0.05))
            n_pos_dir = int(np.nansum((padj < 0.05) & (rho > 0)))

            scan = pd.DataFrame({
                "gene": genes, "virus": virus, "stratum": stratum,
                "n_samples": len(scope_ids),
                "rho": rho, "pvalue": p, "padj": padj,
            }).sort_values("pvalue")
            suffix = f"__{sanitize_virus(stratum)}" if stratum != "all" else ""
            fname = out_dir / f"{sanitize_virus(virus)}{suffix}_genome_scan.tsv"
            scan.to_csv(fname, sep="\t", index=False)

            rows.append({
                "virus": virus, "stratum": stratum,
                "n_samples": len(scope_ids), "n_positives": n_pos,
                "n_genes_tested": int(np.sum(~np.isnan(rho))),
                "n_sig_padj05": n_sig, "n_sig_positive_direction": n_pos_dir,
                "scan_file": fname.name,
            })
            print(f"[{stratum}] {virus}: 阳性 {n_pos}/{len(scope_ids)}  "
                  f"显著基因 {n_sig} (正向 {n_pos_dir})  → {fname.name}")
        return rows

    summary_rows = scan_scope("all", sample_ids)
    if not summary_rows:
        print("[ERROR] 没有任何病毒通过阳性数/变异门槛; "
              "检查标签表或调低 --min-positives", file=sys.stderr)
        sys.exit(1)

    # 组织内自动扫描 (剂量×组织双防混杂; goji: 组织平衡与剂量建模互补)
    if args.meta and not args.tissue:
        meta = pd.read_csv(args.meta, sep="\t", dtype=str)
        if args.strata_col in meta.columns:
            meta["_t"] = meta[args.strata_col].fillna("").str.strip().str.lower()
            meta = meta[meta["_t"] != ""]
            by_tissue = {}
            for s in sample_ids:
                t = meta.loc[meta["Sample"] == s, "_t"]
                if not t.empty:
                    by_tissue.setdefault(t.iloc[0], []).append(s)
            for tissue in sorted(by_tissue):
                ids = by_tissue[tissue]
                if len(ids) < args.min_samples:
                    print(f"[跳过] 组织 {tissue}: 样本 {len(ids)} < {args.min_samples}")
                    continue
                print(f"── 组织内扫描 [{tissue}] ({len(ids)} 样本) ──")
                summary_rows.extend(scan_scope(tissue, ids))
            if any(r["stratum"] != "all" for r in summary_rows):
                print("注意: 组织内结果与全样本结果互为校验 "
                      "(全样本受组织构成混杂影响, 组织内受功效限制)")
        else:
            print(f"[INFO] sample_labels 无 {args.strata_col} 列, 跳过组内自动扫描")

    summary = pd.DataFrame(summary_rows)
    summary.to_csv(out_dir / "load_scan_summary.tsv", sep="\t", index=False)
    print(f"汇总: {out_dir / 'load_scan_summary.tsv'} ({len(summary)} 行, "
          f"stratum 列区分 全样本(all)/组织内)")


if __name__ == "__main__":
    main()


def main_testable(argv):
    """测试入口: 显式传 argv 列表, 不动 sys.argv。"""
    main(argv)
