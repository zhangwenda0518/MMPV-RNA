#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""coexp_network_qc.py — TEA-GCN 网络 z 分布 QC 闸门 (coexp 子步 "qc")
=========================================================================
方法学依据: KingdomStress 预印本 (bioRxiv 10.64898/2026.04.17.719161) Fig S8 —
样本量不足会整体压低 Mutual Rank z-score 分布, 该文 36 物种中 5 个因
max z < 1.7 被剔除模块分析。本脚本对每个 Completed_GCN 网络流式计算:

  Max_z / Mean_z / N_z_ge_1.5 / N_z_ge_1.7 (份额) → Status 判级
    PASS  max_z >= 1.7   (模块检测可行, Δz 可信)
    LOW   1.5 <= max_z < 1.7  (分布压缩, 勉强可用; 模块检测将一无所获)
    FAIL  max_z < 1.5    (分布严重截断, Δz 跨网络比较不可信)

Δz 是两个网络各自 z 标准化的差值: 任一侧分布被样本量压缩, Δz 就被
量纲差主导而非真实 Rewiring — 这是 Δz 的守门员检查 (论文 5/36 物种的教训)。

流式实现: 逐基因文件只取 col0(Target)/col3(zScore), 增量累计 count/sum/max/
阈值计数, 不驻留全量 z (O(N²) 量级不可入内存); 中位数等需全量数据的统计量
不输出。自边 (Target==文件名) 计入 NaN 排除口径。

输入: --gcn-dir <组名>=<Completed_GCN 目录> (可重复, 如 positive=/ negative=/)
输出 (05_Coexp/):
  network_qc.tsv   每网络一行的 z 分布统计 + Status

用法:
  python coexp_network_qc.py --gcn-dir positive=teagcn_positive/Completed_GCN/<net> \
      --gcn-dir negative=teagcn_negative/Completed_GCN/<net> \
      --group-summary matrices/group_summary.tsv --output_dir 05_Coexp
"""

import argparse
import csv
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

Z_LOW = 1.5        # 论文初筛阈值: 低于此 = 分布截断 (FAIL)
Z_MODULE = 1.7     # 论文模块检测阈值: 达不到则模块分析必空 (LOW)


def list_gene_files(gcn_dir: Path):
    return sorted((p for p in gcn_dir.iterdir() if p.is_file()),
                  key=lambda p: p.name)


def scan_network(gcn_dir: Path):
    """流式扫描一个 Completed_GCN 目录, 返回 z 分布统计 dict。"""
    files = list_gene_files(gcn_dir)
    if not files:
        sys.exit(f"[ERROR] 网络目录为空: {gcn_dir}")
    n_edges = n_nan = n_ge15 = n_ge17 = 0
    z_sum = 0.0
    z_max = -np.inf
    targets = set()
    t0 = time.time()
    for i, f in enumerate(files):
        df = pd.read_csv(f, sep="\t", usecols=[0, 3])
        targets.update(df.iloc[:, 0].astype(str))
        z = pd.to_numeric(df.iloc[:, 1], errors="coerce").to_numpy(dtype=float)
        # 自边: Target==源基因名 (文件名) — 与 diff_coexp 同口径剔除
        src_self = (df.iloc[:, 0].astype(str).to_numpy() == f.name)
        keep = ~np.isnan(z) & ~src_self
        valid = z[keep]
        n_nan += int(np.isnan(z).sum())
        n_edges += int(keep.sum())
        if valid.size:
            z_sum += float(valid.sum())
            z_max = max(z_max, float(valid.max()))
            n_ge15 += int((valid >= Z_LOW).sum())
            n_ge17 += int((valid >= Z_MODULE).sum())
        if (i + 1) % 2000 == 0:
            print(f"  {gcn_dir.name}: {i + 1}/{len(files)}  {time.time() - t0:.0f}s",
                  flush=True)
    mean_z = z_sum / n_edges if n_edges else float("nan")
    return {
        "N_Genes": len(files), "N_Edges": n_edges, "N_SelfEdge_NaN": n_nan,
        "Max_z": round(z_max, 4) if np.isfinite(z_max) else float("nan"),
        "Mean_z": round(mean_z, 4) if n_edges else float("nan"),
        "N_z_ge_1.5": n_ge15,
        "Share_z_ge_1.5": round(n_ge15 / n_edges, 6) if n_edges else float("nan"),
        "N_z_ge_1.7": n_ge17,
        "Share_z_ge_1.7": round(n_ge17 / n_edges, 6) if n_edges else float("nan"),
        "N_Targets": len(targets),
    }


def judge(stats: dict) -> str:
    mz = stats["Max_z"]
    if mz != mz or mz < Z_LOW:          # NaN 或低于初筛线
        return "FAIL"
    if mz < Z_MODULE:
        return "LOW"
    return "PASS"


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(
        description="TEA-GCN 网络 z 分布 QC (KingdomStress Fig S8 物种剔除规则的工程化)")
    ap.add_argument("--gcn-dir", action="append", required=True,
                    help="组名=Completed_GCN 目录 (可重复, 如 positive=... negative=...)")
    ap.add_argument("--group-summary",
                    help="build_coexp_matrices 的 group_summary.tsv (并入每组样本数)")
    ap.add_argument("--output_dir", required=True, help="输出目录 (05_Coexp)")
    ap.add_argument("--on-fail", choices=["warn", "error"], default="warn",
                    help="存在 FAIL 网络时: warn=打印警告照常退出 / error=exit 3 "
                         "(编排器默认自行读 network_qc.tsv 决定闸门, 此参数供单跑)")
    args = ap.parse_args(argv)

    n_samples = {}
    if args.group_summary and Path(args.group_summary).is_file():
        with open(args.group_summary, encoding="utf-8") as fh:
            rows = list(csv.DictReader(fh, delimiter="\t"))
        if rows:
            n_samples = {"positive": rows[0].get("n_positive", ""),
                         "negative": rows[0].get("n_negative", "")}

    records = []
    for spec in args.gcn_dir:
        group, _, path = spec.partition("=")
        if not path:
            sys.exit(f"[ERROR] --gcn-dir 需要 组名=路径 形式, 得到: {spec}")
        print(f"扫描网络 [{group}]: {path}")
        stats = scan_network(Path(path))
        status = judge(stats)
        stats.update({"Group": group, "N_Samples": n_samples.get(group, ""),
                      "Status": status})
        records.append(stats)
        print(f"  → Max_z={stats['Max_z']}  z≥1.5: {stats['N_z_ge_1.5']} "
              f"({stats['Share_z_ge_1.5']:.2%})  z≥1.7: {stats['N_z_ge_1.7']} "
              f"→ {status}")

    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    cols = ["Group", "N_Genes", "N_Edges", "N_SelfEdge_NaN", "N_Targets",
            "Max_z", "Mean_z", "N_z_ge_1.5", "Share_z_ge_1.5",
            "N_z_ge_1.7", "Share_z_ge_1.7", "N_Samples", "Status"]
    pd.DataFrame(records, columns=cols).to_csv(
        out / "network_qc.tsv", sep="\t", index=False)
    print(f"输出: {out / 'network_qc.tsv'}")

    fails = [r["Group"] for r in records if r["Status"] == "FAIL"]
    lows = [r["Group"] for r in records if r["Status"] == "LOW"]
    if fails:
        msg = (f"网络 z 分布截断 (max_z<{Z_LOW}): {fails} — Δz 跨网络比较不可信 "
               f"(参考 KingdomStress Fig S8: 样本量不足压低 MR z 分布)")
        if args.on_fail == "error":
            print(f"[QC-FAIL] {msg}", file=sys.stderr)
            sys.exit(3)
        print(f"[QC-WARN] {msg}", file=sys.stderr)
    if lows:
        print(f"[QC-WARN] 网络分布偏低 (1.5≤max_z<1.7, 模块检测将无产出): {lows}",
              file=sys.stderr)


if __name__ == "__main__":
    main()


def main_testable(argv):
    """测试入口: 显式传 argv 列表, 不动 sys.argv。"""
    main(argv)
