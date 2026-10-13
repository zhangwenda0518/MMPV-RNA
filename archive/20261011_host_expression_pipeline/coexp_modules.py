#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""coexp_modules.py — TEA-GCN 网络 → DEG 限定共表达模块检测 (coexp 子步 "modules")
=========================================================================
方法学依据: KingdomStress 预印本 (bioRxiv 10.64898/2026.04.17.719161) Fig 5 —
"对每个物种×胁迫组合, 网络过滤到 z >= 1.7 且两端均为差异表达基因的边,
再做 Louvain 社区发现 (resolution = 1.0) 得到胁迫响应共表达模块"。

流程 (每组网络独立):
  1. 只读 DEG 基因的逐基因文件 (文件名 = 源基因; 其余文件直接跳过, 省一遍全网络 IO)
  2. 取 z >= --z-threshold 且 Target != 源基因 且两端均在 DEG 表内的边 → networkx Graph
  3. Louvain (weight=z, resolution=--resolution, seed=42; networkx<2.8 回退
     greedy_modularity_communities 并打印说明)
  4. 模块×病毒载量 Spearman 关联 (可选 --tpm + --virus-labels): 模块成员基因
     TPM 横向 z-score 取均值 → 模块表达谱, 与每种病毒的连续载量相关;
     全部 (模块, 病毒) 行做 BH-FDR。通用分组模式 (无 virus_labels) 自动跳过。

前置: coexp_network_qc.py 的 network_qc.tsv (若有): Max_z < --z-threshold 的
组按 Fig S8 规则跳过 (分布压缩时模块检测必一无所获)。

输入: --gcn-dir 组名=Completed_GCN 目录 (可重复), --deg deseq2_degs.tsv
输出 (--output_dir, 建议 05_Coexp/modules/):
  module_table_<组>.tsv   Gene, Module, Degree (长表)
  module_edges_<组>.tsv   Gene1, Gene2, Z, Module1, Module2 (网络边表, 可视化用)
  module_summary.tsv      Group, Module, N_Genes, Mean_Degree, Max_Degree, Hub_Genes
  module_load_assoc.tsv   Group, Module, N_Genes, Virus, N_Samples, Rho, Pvalue, Padj
"""

import argparse
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

try:
    from networkx.algorithms.community import louvain_communities
    _LOUVAIN = True
except ImportError:                                    # networkx < 2.8
    from networkx.algorithms.community import greedy_modularity_communities
    _LOUVAIN = False
import networkx as nx


def read_deg(path: Path) -> set:
    """deseq2_degs.tsv: 首列 = 基因名 (run_deseq2 写出时 index 名为 gene)。"""
    df = pd.read_csv(path, sep="\t", index_col=0)
    return set(df.index.astype(str))


def read_qc_maxz(path: Path) -> dict:
    """network_qc.tsv: Group → Max_z (缺失文件返回空 dict)。"""
    if not path.is_file():
        return {}
    df = pd.read_csv(path, sep="\t")
    if "Group" not in df.columns or "Max_z" not in df.columns:
        return {}
    return dict(zip(df["Group"].astype(str), pd.to_numeric(df["Max_z"],
                                                           errors="coerce")))


def build_graph(gcn_dir: Path, deg: set, z_threshold: float):
    """只扫 DEG 基因文件, 收集 z≥阈值 且两端在 DEG 内的非自边 → 无向图。"""
    G = nx.Graph()
    files = [f for f in sorted(gcn_dir.iterdir(), key=lambda p: p.name)
             if f.is_file() and f.name in deg]
    if not files:
        return G, 0
    t0 = time.time()
    for i, f in enumerate(files):
        df = pd.read_csv(f, sep="\t", usecols=[0, 3])
        targets = df.iloc[:, 0].astype(str).to_numpy()
        z = pd.to_numeric(df.iloc[:, 1], errors="coerce").to_numpy(dtype=float)
        for t, zv in zip(targets, z):
            if t == f.name or t not in deg or zv != zv or zv < z_threshold:
                continue
            if G.has_edge(f.name, t):
                if zv > G[f.name][t]["weight"]:
                    G[f.name][t]["weight"] = float(zv)
            else:
                G.add_edge(f.name, t, weight=float(zv))
        if (i + 1) % 2000 == 0:
            print(f"  {gcn_dir.name}: {i + 1}/{len(files)} DEG 文件  "
                  f"{time.time() - t0:.0f}s", flush=True)
    return G, len(files)


def detect_modules(G, resolution: float, seed: int):
    if _LOUVAIN:
        return louvain_communities(G, weight="weight", resolution=resolution,
                                   seed=seed)
    print("[WARN] networkx 无 louvain_communities (<2.8), 回退 "
          "greedy_modularity_communities (分辨率参数不可用)")
    return greedy_modularity_communities(G, weight="weight")


def module_tables(G, communities, group: str):
    """长表 (Gene, Module, Degree) + 汇总行 (Group, Module, N_Genes, ...)。"""
    rows, sums = [], []
    for mid, comm in enumerate(sorted(communities, key=len, reverse=True)):
        members = sorted(comm)
        degs = sorted(((G.degree[n], n) for n in members), reverse=True)
        for d, g in degs:
            rows.append((g, mid, int(d)))
        hub = ";".join(g for _, g in degs[:3])
        dl = [d for d, _ in degs]
        sums.append({"Group": group, "Module": mid, "N_Genes": len(members),
                     "Mean_Degree": round(float(np.mean(dl)), 2) if dl else 0.0,
                     "Max_Degree": int(max(dl)) if dl else 0,
                     "Hub_Genes": hub})
    table = pd.DataFrame(rows, columns=["Gene", "Module", "Degree"])
    summary = pd.DataFrame(sums)
    return table, summary


def _bh_fdr(pvals: np.ndarray) -> np.ndarray:
    p = np.asarray(pvals, dtype=float)
    out = np.full(p.shape, np.nan)
    ok = ~np.isnan(p)
    pv = p[ok]
    if pv.size == 0:
        return out
    order = np.argsort(pv)
    ranked = pv[order] * pv.size / (np.arange(pv.size) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    res = np.empty(pv.size)
    res[order] = np.clip(ranked, 0, 1)
    out[ok] = res
    return out


def load_association(mod_tables: dict, tpm_path: Path, vlabels_path: Path,
                     min_samples: int = 8, samples_filter=None):
    """模块表达谱 × 每种病毒连续载量 Spearman; 全部行 BH-FDR。

    samples_filter: QC 通过样本集合 (None = 不过滤)。"""
    from scipy.stats import spearmanr
    tpm = pd.read_csv(tpm_path, sep="\t", index_col=0)
    vlab = pd.read_csv(vlabels_path, sep="\t")
    need = {"Sample", "Virus", "Load"}
    if not need.issubset(vlab.columns):
        print(f"[WARN] virus_labels 缺列 {need - set(vlab.columns)}; 跳过载量关联")
        return None
    allowed = set(samples_filter) if samples_filter is not None else None
    samples = [s for s in tpm.columns
               if s in set(vlab["Sample"].astype(str))
               and (allowed is None or s in allowed)]
    if len(samples) < min_samples:
        print(f"[WARN] 可用样本 {len(samples)} < {min_samples}; 跳过载量关联")
        return None
    tpm = tpm[samples]
    rows = []
    for group, table in mod_tables.items():
        if table.empty:
            continue
        for mid, grp in table.groupby("Module"):
            genes = [g for g in grp["Gene"] if g in tpm.index]
            if not genes:
                continue
            X = tpm.loc[genes].to_numpy(dtype=float)
            mu = np.nanmean(X, axis=1, keepdims=True)
            sd = np.nanstd(X, axis=1, keepdims=True)
            sd[sd == 0] = np.nan
            Z = (X - mu) / sd
            with np.errstate(invalid="ignore"):
                profile = np.nanmean(Z, axis=0)
            for virus, vg in vlab.groupby("Virus"):
                load_map = dict(zip(vg["Sample"].astype(str),
                                    pd.to_numeric(vg["Load"], errors="coerce")))
                xs, ys = [], []
                for s, pv in zip(samples, profile):
                    lv = load_map.get(s)
                    if lv is None or lv != lv or pv != pv:
                        continue
                    xs.append(pv)
                    ys.append(lv)
                if len(xs) < min_samples or np.std(ys) == 0:
                    continue
                rho, p = spearmanr(xs, ys)
                rows.append({"Group": group, "Module": int(mid),
                             "N_Genes": len(genes), "Virus": virus,
                             "N_Samples": len(xs),
                             "Rho": round(float(rho), 4),
                             "Pvalue": float(p)})
    if not rows:
        return None
    res = pd.DataFrame(rows)
    res["Padj"] = _bh_fdr(res["Pvalue"].to_numpy())
    return res.sort_values("Padj")


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(
        description="TEA-GCN 网络 DEG 限定共表达模块检测 (KingdomStress Fig5 流程)")
    ap.add_argument("--gcn-dir", action="append", required=True,
                    help="组名=Completed_GCN 目录 (可重复)")
    ap.add_argument("--deg", required=True, help="DEG 表 (deseq2_degs.tsv; 首列基因)")
    ap.add_argument("--network-qc", help="coexp_network_qc 的 network_qc.tsv "
                                         "(Max_z<阈值 的组自动跳过)")
    ap.add_argument("--z-threshold", type=float, default=1.7,
                    help="模块边阈值 (论文 Fig5: 1.7)")
    ap.add_argument("--resolution", type=float, default=1.0,
                    help="Louvain resolution (论文 Fig5: 1.0)")
    ap.add_argument("--seed", type=int, default=42, help="Louvain 随机种子")
    ap.add_argument("--tpm", help="TPM 矩阵 (载量关联用; 与 --virus-labels 成对)")
    ap.add_argument("--virus-labels", help="virus_labels.tsv (Sample/Virus/Load)")
    ap.add_argument("--min-assoc-samples", type=int, default=8,
                    help="载量关联最少样本数 (默认 8)")
    ap.add_argument("--samples",
                    help="QC 通过样本名单 (03_QC/qc_passed_samples.txt); "
                         "限定载量关联的样本范围")
    ap.add_argument("--output_dir", required=True, help="输出目录 (05_Coexp/modules)")
    args = ap.parse_args(argv)

    deg = read_deg(Path(args.deg))
    print(f"DEG: {len(deg)} 个 (z>={args.z_threshold:g}, Louvain "
          f"resolution={args.resolution:g}, seed={args.seed})")
    if not deg:
        sys.exit("[ERROR] DEG 表为空 — 06_DE 未产出显著基因; 模块检测无从限定")
    maxz = read_qc_maxz(Path(args.network_qc)) if args.network_qc else {}

    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    mod_tables, sum_rows = {}, []
    for spec in args.gcn_dir:
        group, _, path = spec.partition("=")
        if not path:
            sys.exit(f"[ERROR] --gcn-dir 需要 组名=路径 形式, 得到: {spec}")
        if group in maxz and maxz[group] == maxz[group] \
                and maxz[group] < args.z_threshold:
            print(f"[跳过] [{group}] 网络 Max_z={maxz[group]:.3f} < "
                  f"{args.z_threshold:g} (Fig S8: 分布压缩, 模块检测必空)")
            mod_tables[group] = pd.DataFrame(columns=["Gene", "Module", "Degree"])
            continue
        print(f"构建 DEG 限定网络 [{group}]: {path}")
        G, n_read = build_graph(Path(path), deg, args.z_threshold)
        print(f"  DEG 文件 {n_read}, 图: {G.number_of_nodes()} 节点 / "
              f"{G.number_of_edges()} 边")
        if G.number_of_nodes() == 0:
            print(f"[WARN] [{group}] 无 z≥{args.z_threshold:g} 且两端均为 DEG 的边; "
                  f"该组无模块")
            mod_tables[group] = pd.DataFrame(columns=["Gene", "Module", "Degree"])
            continue
        comms = detect_modules(G, args.resolution, args.seed)
        table, summary = module_tables(G, comms, group)
        mod_tables[group] = table
        sum_rows.append(summary)
        print(f"  模块数: {len(comms)} (最大 {max(len(c) for c in comms)} 基因)")
        table.to_csv(out / f"module_table_{group}.tsv", sep="\t", index=False)

        # 边表导出 (网络可视化原料: 模块内布局图按此绘边)
        gmap = dict(zip(table["Gene"], table["Module"]))
        edges = [(u, v, float(d.get("weight", 0.0)),
                  gmap.get(u, -1), gmap.get(v, -1))
                 for u, v, d in G.edges(data=True)]
        pd.DataFrame(edges, columns=["Gene1", "Gene2", "Z",
                                     "Module1", "Module2"]).to_csv(
            out / f"module_edges_{group}.tsv", sep="\t", index=False)

    if not sum_rows:
        print("[WARN] 所有组均无模块产出; 仅写空 summary")
        pd.DataFrame(columns=["Group", "Module", "N_Genes", "Mean_Degree",
                              "Max_Degree", "Hub_Genes"]).to_csv(
            out / "module_summary.tsv", sep="\t", index=False)
        return
    summary = pd.concat(sum_rows, ignore_index=True)
    summary.to_csv(out / "module_summary.tsv", sep="\t", index=False)
    print(f"输出: {out / 'module_summary.tsv'} ({len(summary)} 模块)")

    if args.tpm and args.virus_labels:
        qc_ids = None
        if args.samples:
            with open(args.samples, encoding="utf-8") as fh:
                qc_ids = {l.strip() for l in fh if l.strip()}
        assoc = load_association(mod_tables, Path(args.tpm),
                                 Path(args.virus_labels), args.min_assoc_samples,
                                 samples_filter=qc_ids)
        if assoc is not None and len(assoc):
            assoc.to_csv(out / "module_load_assoc.tsv", sep="\t", index=False)
            n_sig = int((assoc["Padj"] < 0.05).sum())
            print(f"输出: {out / 'module_load_assoc.tsv'} "
                  f"({len(assoc)} 行, Padj<0.05 {n_sig} 行)")
        else:
            print("载量关联: 无可计算行 (样本不足/载量无变异), 未写文件")


if __name__ == "__main__":
    main()


def main_testable(argv):
    """测试入口: 显式传 argv 列表, 不动 sys.argv。"""
    main(argv)
