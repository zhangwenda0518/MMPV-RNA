#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""rank_transform_mp.py — TEA-GCN Rank_transform 的 fork 分片并行版
=========================================================================
复刻 teagcn/main/Rank_transform.py 的 full_attributes=False 路径
(extract_ranks → MR=√(upper·lower) → z-score → write_MR_attributes_only),
把两处 O(N²) 的逐行 Python 循环 (读 RAW 与写 Completed) 按源基因分片到
fork 多进程; 计算语义与串行版完全一致, 输出目标为**逐字节一致**
(同一标准: 当年 mp 网络构建 39,224/39,224 md5 一致)。

并行正确性 (免锁, 每数组单写者):
  边 (a,b), a<b 按 gene_dict 下标: upper 只由 source a 写, lower 只由
  source b 写 — 天然单写者, 无需锁。cor 列不进共享数组 (write 阶段从
  RAW 文件原样重读, 与串行版语义一致)。

内存: 2 个 float64 共享数组 ≈ 8B×N(N-1)/2×2 (39,224 基因 ≈ 12.3 GB),
worker 经 fork 继承映射; 主进程额外驻留 MR/zMR 两个同规模数组。

用法 (与 Rank_transform.py 同参 + -w):
  python rank_transform_mp.py -o <work_dir> [-n net_name] [-w workers]

独立实现 (不 import 上游模块); Windows 无 fork 时自动退化为单进程执行
(同一代码路径, 结果不变)。
"""

import argparse
import os
import sys
import time

import numpy as np


def calculate_edge_index(s, t, num_nodes):
    """与上游 data_processing.network.calculate_edge_index 逐字一致。"""
    if s > t:
        s, t = t, s  # Ensure s < t
        return False, num_nodes * s - ((s + 1) * s) // 2 + t - s - 1
    else:
        return True, num_nodes * s - ((s + 1) * s) // 2 + t - s - 1


def get_zscore(np_array, flip_scores=True):
    """与上游 Rank_transform.get_zscore 逐字一致。"""
    std_deviation = np.nanstd(np_array)
    ret_array = np_array - np.nanmean(np_array)
    if flip_scores:
        return -1 * (ret_array / std_deviation)
    else:
        return (ret_array / std_deviation)


def list_network_paths(output_dir):
    """与上游相同: 收集 RAW_GCN 下含 'k' 的网络目录 (PCC/SCC/bicor/TEA)。"""
    bgcn = os.path.join(output_dir, "RAW_GCN")
    paths = []
    for cc in ["PCC", "SCC", "bicor", "TEA"]:
        d = os.path.join(bgcn, cc)
        if os.path.exists(d):
            for netdir in os.listdir(d):
                if "k" in netdir:
                    paths.append(os.path.join(d, netdir))
    return paths


def _extract_worker(shard, genes, gene_dict, num_nodes, net_dir,
                    upper_arr, lower_arr):
    """分片解析源基因文件。每数组仅本 worker 可写的元素被写入 (见模块头)。"""
    for source in shard:
        s_idx = gene_dict[source]
        with open(os.path.join(net_dir, source), "r") as f:
            for line_no, line in enumerate(f):
                if line_no != 0 and line != "":
                    target, cor, rank = line.split("\n")[0].split("\t")
                    if target == source:
                        continue
                    t_idx = gene_dict[target]
                    upper, ei = calculate_edge_index(s_idx, t_idx, num_nodes)
                    if upper:
                        upper_arr[ei] = rank
                    else:
                        lower_arr[ei] = rank
    return len(shard)


def _write_worker(shard, genes, gene_dict, num_nodes, net_dir, out_dir,
                  MR_arr, zMR_arr, max_zscore_MR):
    """分片写 Completed 文件 — write_MR_attributes_only 循环体的逐字复刻。"""
    for source in shard:
        with open(os.path.join(net_dir, source), "r") as fin:
            with open(os.path.join(out_dir, source), "w") as fout:
                for line_no, line in enumerate(fin):
                    if line_no == 0:
                        fout.write("Target\tCo-exp_Str\tCo-exp_Str_MR\t"
                                   "zScore(Co-exp_Str_MR)\n")
                    elif line != "":
                        target, cor, rank = line.split("\n")[0].split("\t")
                        cor = np.round(float(cor), 5)
                        if target == source:
                            mr, zmr = 1.0, max_zscore_MR
                        else:
                            s_idx = gene_dict[source]
                            t_idx = gene_dict[target]
                            _, ei = calculate_edge_index(s_idx, t_idx, num_nodes)
                            mr, zmr = MR_arr[ei], zMR_arr[ei]
                            mr, zmr = np.round(mr, 3), np.round(zmr, 7)
                        fout.write(f"{target}\t{cor}\t{mr}\t{zmr}\n")
    return len(shard)


def _shard_list(items, n):
    """把列表切成 n 个尽量均匀的连续片段 (含空片段)。"""
    n = max(1, n)
    k, m = divmod(len(items), n)
    out, start = [], 0
    for i in range(n):
        end = start + k + (1 if i < m else 0)
        out.append(items[start:end])
        start = end
    return out


def process_network(net_dir, out_dir, workers):
    """对单个网络目录执行 并行 extract → 串行 MR/z → 并行 write。"""
    t0 = time.time()
    genes = os.listdir(net_dir)                    # 与串行版同序 (listdir)
    gene_dict = {g: i for i, g in enumerate(genes)}
    num_nodes = len(genes)
    num_edges = num_nodes * (num_nodes - 1) // 2
    print(f"网络: {net_dir}\n基因 {num_nodes}, 边 {num_edges}, "
          f"workers {workers}")

    use_mp = workers > 1 and hasattr(os, "fork")
    if workers > 1 and not use_mp:
        print("[WARN] 本平台无 fork, 退化为单进程 (结果不变)")

    # ── 并行 extract (共享数组: fork 继承, 每元素单写者免锁) ──
    if use_mp:
        import multiprocessing as mp
        arrays = [mp.Array("d", num_edges, lock=False) for _ in range(2)]
        upper_arr, lower_arr = (np.frombuffer(a, dtype=np.float64)
                                for a in arrays)
        upper_arr[:] = np.nan
        lower_arr[:] = np.nan
        shards = _shard_list(genes, workers)
        procs = []
        for shard in shards:
            if not shard:
                continue
            p = mp.get_context("fork").Process(
                target=_extract_worker,
                args=(shard, genes, gene_dict, num_nodes, net_dir,
                      upper_arr, lower_arr))
            p.start()
            procs.append(p)
        done = 0
        for p in procs:
            p.join()
            done += 1
            if p.exitcode != 0:
                sys.exit(f"[ERROR] extract worker 异常退出 "
                         f"(exitcode={p.exitcode}); 结果不完整, 请重跑")
            print(f"extract shard {done}/{len(procs)} joined  "
                  f"{time.time() - t0:.0f}s")
    else:
        upper_arr = np.full(num_edges, np.nan)
        lower_arr = np.full(num_edges, np.nan)
        _extract_worker(genes, genes, gene_dict, num_nodes, net_dir,
                        upper_arr, lower_arr)
        print(f"extract (单进程) 完成  {time.time() - t0:.0f}s")

    # ── MR / z-score (串行, numpy 向量化, 与串行版同式) ──
    MR_arr = np.sqrt(upper_arr * lower_arr)
    del upper_arr, lower_arr
    zMR_arr = get_zscore(MR_arr, flip_scores=True)
    max_zscore_MR = np.nanmax(zMR_arr)
    print(f"MR/z-score 完成  {time.time() - t0:.0f}s")

    # ── 并行 write ──
    os.makedirs(out_dir, exist_ok=True)
    if use_mp:
        import multiprocessing as mp
        shards = [s for s in _shard_list(genes, workers) if s]
        procs = []
        for shard in shards:
            p = mp.get_context("fork").Process(
                target=_write_worker,
                args=(shard, genes, gene_dict, num_nodes, net_dir, out_dir,
                      MR_arr, zMR_arr, max_zscore_MR))
            p.start()
            procs.append(p)
        for i, p in enumerate(procs, 1):
            p.join()
            if p.exitcode != 0:
                sys.exit(f"[ERROR] write worker 异常退出 "
                         f"(exitcode={p.exitcode}); 结果不完整, 请重跑")
            print(f"write shard {i}/{len(procs)} joined  "
                  f"{time.time() - t0:.0f}s")
    else:
        _write_worker(genes, genes, gene_dict, num_nodes, net_dir, out_dir,
                      MR_arr, zMR_arr, max_zscore_MR)
        print(f"write (单进程) 完成  {time.time() - t0:.0f}s")
    print(f"网络完成: {out_dir}  总耗时 {time.time() - t0:.0f}s")


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(
        description="Rank_transform 并行版 (输出与串行逐字节一致)")
    ap.add_argument("-o", "--output_dir", required=True,
                    help="Run_TEA-GCN 的工作目录 (含 RAW_GCN/)")
    ap.add_argument("-n", "--net_name", default="All",
                    help="网络名过滤 (默认 All)")
    ap.add_argument("-w", "--workers", type=int, default=8,
                    help="fork 进程数 (extract 与 write 各一轮)")
    args = ap.parse_args(argv)

    paths = list_network_paths(args.output_dir)
    if args.net_name != "All":
        paths = [p for p in paths if args.net_name.lower() in p.lower()]
    if not paths:
        sys.exit(f"[ERROR] RAW_GCN 下未找到匹配的网络: {args.output_dir}")
    print("待处理网络:")
    for p in paths:
        print(f"  {p}")

    for net_dir in paths:
        out_dir = os.path.join(
            args.output_dir, "Completed_GCN",
            os.path.basename(net_dir.replace("\\", "/").rstrip("/")))
        process_network(net_dir, out_dir, args.workers)
    print("rank_transform_mp 全部完成")


if __name__ == "__main__":
    main()
