#!/usr/bin/env python3
"""
virus_distribution.py — 病毒分布分析 v1.0
==========================================
输入: cluster 结果目录 (含 vclust_clusters.tsv + centroids + rmdup)
输出: 09_Virome_Analysis/ 下的分布表格

产出:
  sample_virus_counts.tsv      每个样本的病毒 contig 数 + centroids 数
  cluster_sample_matrix.tsv    每个 cluster 在哪些样本中出现
  cluster_summary_enriched.tsv 每个 cluster 的样本数、contig 数、centroids 长度
  rmdup_stats.tsv             rmdup 去冗余统计
  flye_contribution.tsv       Flye 共组装贡献统计 (如果存在)
"""

import argparse, sys, os
from pathlib import Path
from collections import defaultdict, Counter

try:
    import polars as pl
    HAS_POLARS = True
except ImportError:
    HAS_POLARS = False

try:
    import pandas as pd
    HAS_PANDAS = True
except ImportError:
    HAS_PANDAS = False

from Bio import SeqIO


def load_clusters(ctsv_path):
    """解析 vclust_clusters.tsv → {cluster_id: [member_ids]}"""
    clusters = defaultdict(list)
    with open(ctsv_path) as f:
        f.readline()  # skip header
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 2:
                clusters[parts[1].strip()].append(parts[0].strip())
    # 去重成员
    return {cid: list(set(members)) for cid, members in clusters.items()}


def extract_sample(contig_id):
    """从 contig ID 提取样本名"""
    if contig_id.startswith("flye_"):
        return "__FLYE__"  # Flye 共组装, 不归属单样本
    if "_all_" in contig_id:
        return contig_id.split("_all_")[0]
    return contig_id.split("_")[0]


def load_host_map(host_tsv):
    """加载宿主预测 → {contig_id: host_category}"""
    host_map = {}
    if not host_tsv or not os.path.isfile(host_tsv):
        return host_map
    with open(host_tsv) as f:
        header = f.readline().rstrip('\n').split('\t')
        try:
            cid_col = header.index("contig_id")
        except ValueError:
            cid_col = 0
        try:
            host_col = header.index("Final_Host")
        except ValueError:
            host_col = 1
        for line in f:
            cols = line.rstrip('\n').split('\t')
            if len(cols) > max(cid_col, host_col):
                cid = cols[cid_col]
                host = cols[host_col] if cols[host_col] and cols[host_col] not in ("NA", "N/A", "") else "Unknown"
                host_map[cid] = host
    return host_map


def load_centroids(fasta_path):
    """加载 centroids → {id: length}"""
    centroids = {}
    for rec in SeqIO.parse(fasta_path, "fasta"):
        centroids[rec.id] = len(rec.seq)
    return centroids


def categorize_host(host_name):
    """将宿主名归类到大类"""
    host_lower = host_name.lower()
    if "plant" in host_lower or "viridiplantae" in host_lower or "streptophyta" in host_lower:
        return "Plant"
    if "bacteria" in host_lower or "proteobacteria" in host_lower or "firmicutes" in host_lower \
       or "actinobacteria" in host_lower or "cyanobacteria" in host_lower:
        return "Bacteria"
    if "fungi" in host_lower or "ascomycota" in host_lower or "basidiomycota" in host_lower:
        return "Fungi"
    if "invertebrate" in host_lower or "insect" in host_lower or "arthropod" in host_lower \
       or "nematod" in host_lower or "mollus" in host_lower or "annelid" in host_lower:
        return "Invertebrate"
    if "vertebrate" in host_lower or "mammal" in host_lower or "aves" in host_lower \
       or "fish" in host_lower or "reptil" in host_lower or "amphib" in host_lower \
       or "human" in host_lower or "homo" in host_lower:
        return "Vertebrate"
    if "archae" in host_lower:
        return "Archaea"
    if host_name == "Unknown" or host_name == "":
        return "Unknown"
    return "Other"
    """加载 centroids → {id: length}"""
    centroids = {}
    for rec in SeqIO.parse(fasta_path, "fasta"):
        centroids[rec.id] = len(rec.seq)
    return centroids


def main():
    p = argparse.ArgumentParser(description="病毒分布分析")
    p.add_argument("--cluster-dir", required=True, help="04_CLUSTER 目录")
    p.add_argument("--output-dir", required=True, help="输出目录")
    p.add_argument("--host-tsv", default=None, help="宿主预测 TSV (ensemble_host_summary.tsv)")
    p.add_argument("--host-filter", default="Plant", help="目标宿主 (默认: Plant)")
    args = p.parse_args()

    cluster_dir = Path(args.cluster_dir)
    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # ── 加载数据 ──
    ctsv = cluster_dir / "3_vclust" / "vclust_clusters.tsv"
    centroids_fa = cluster_dir / "4_centroids" / "final_centroids.fasta"
    if not centroids_fa.is_file():
        centroids_fa = cluster_dir / "04_centroids" / "final_centroids.fasta"  # 旧回退
    if not centroids_fa.is_file():
        centroids_fa = cluster_dir / "centroids" / "final_centroids.fasta"  # 更旧回退
    rmdup_fa = cluster_dir / "3_vclust" / "all.cluster.ref.rmdup.fasta"

    if not ctsv.is_file():
        print(f"错误: clusters 文件不存在 {ctsv}", file=sys.stderr)
        sys.exit(1)
    if not centroids_fa.is_file():
        print(f"错误: centroids 不存在 {centroids_fa}", file=sys.stderr)
        sys.exit(1)

    print(f"加载 clusters: {ctsv}")
    clusters = load_clusters(ctsv)
    print(f"  {len(clusters)} 个 cluster")

    print(f"加载 centroids: {centroids_fa}")
    centroids = load_centroids(centroids_fa)
    print(f"  {len(centroids)} 条 centroids")

    # 加载宿主预测
    host_map = load_host_map(args.host_tsv)
    if host_map:
        print(f"加载宿主预测: {args.host_tsv} ({len(host_map)} 条)")

    # ── 1. 每样本病毒分布 ──
    print("\n── 1. 样本级病毒分布 ──")
    sample_contigs = Counter()
    sample_centroids = Counter()
    sample_host_cats = defaultdict(Counter)  # sample → {host_cat: centroids_count}
    flye_contig_count = 0
    flye_centroid_count = 0

    for cid, members in clusters.items():
        c_samples = set()
        for m in members:
            s = extract_sample(m)
            if s == "__FLYE__":
                flye_contig_count += 1
            else:
                sample_contigs[s] += 1
                c_samples.add(s)
        # centroids 归属: 代表序列的样本
        rep = min(members, key=lambda x: len(x)) if members else ""
        rep_sample = extract_sample(rep)
        rep_host_cat = categorize_host(host_map.get(rep, "Unknown")) if host_map else ""
        if rep_sample == "__FLYE__":
            flye_centroid_count += 1
        else:
            for s in c_samples:
                sample_centroids[s] += 1
                if rep_host_cat:
                    sample_host_cats[s][rep_host_cat] += 1

    rows = []
    all_samples = sorted(set(list(sample_contigs.keys()) + list(sample_centroids.keys())))
    # 收集所有宿主大类作为列
    all_host_cats = set()
    for cats in sample_host_cats.values():
        all_host_cats.update(cats.keys())
    host_cat_columns = sorted(all_host_cats)

    for s in all_samples:
        row = {
            "sample": s,
            "contigs_in_clusters": sample_contigs.get(s, 0),
            "centroids": sample_centroids.get(s, 0),
        }
        for cat in host_cat_columns:
            row[cat] = sample_host_cats.get(s, {}).get(cat, 0)
        rows.append(row)

    if flye_contig_count > 0:
        rows.append({
            "sample": "__FLYE__ (共组装)",
            "contigs_in_clusters": flye_contig_count,
            "centroids": flye_centroid_count,
        })

    _write_table(rows, out_dir / "sample_virus_counts.tsv", "sample")
    print(f"  → {out_dir / 'sample_virus_counts.tsv'} ({len(rows)} 行)")

    # ── 2. Cluster × Sample 矩阵 ──
    print("\n── 2. Cluster 样本分布矩阵 ──")
    cluster_sample_map = {}  # cluster_id → {sample: contig_count}
    cluster_info = []        # 每 cluster 汇总

    for cid, members in sorted(clusters.items()):
        sample_counts = Counter()
        flye_members = 0
        for m in members:
            s = extract_sample(m)
            if s == "__FLYE__":
                flye_members += 1
            else:
                sample_counts[s] += 1

        n_samples = len(sample_counts)
        n_contigs = len(members)
        rep = min(members, key=lambda x: len(x)) if members else ""
        rep_len = centroids.get(rep, 0)
        is_flye = rep.startswith("flye_")

        cluster_sample_map[cid] = dict(sample_counts)
        cluster_info.append({
            "cluster_id": cid,
            "num_contigs": n_contigs,
            "num_samples": n_samples,
            "centroid_length": rep_len,
            "centroid_id": rep[:80],
            "is_flye_centroid": "Yes" if is_flye else "",
            "flye_members": flye_members,
        })

    # 写入矩阵
    all_samples_set = set()
    for smap in cluster_sample_map.values():
        all_samples_set.update(smap.keys())
    all_samples_list = sorted(all_samples_set)

    matrix_rows = []
    for cid, smap in sorted(cluster_sample_map.items()):
        row = {"cluster_id": cid}
        for s in all_samples_list:
            row[s] = smap.get(s, 0)
        row["total"] = sum(smap.values())
        matrix_rows.append(row)

    _write_table(matrix_rows, out_dir / "cluster_sample_matrix.tsv", "cluster_id")
    print(f"  → {out_dir / 'cluster_sample_matrix.tsv'} ({len(matrix_rows)} 行 × {len(all_samples_list)+2} 列)")

    # 写入 cluster 汇总
    _write_table(cluster_info, out_dir / "cluster_summary_enriched.tsv", "cluster_id")
    print(f"  → {out_dir / 'cluster_summary_enriched.tsv'} ({len(cluster_info)} 行)")

    # ── 3. rmdup 统计 ──
    print("\n── 3. rmdup 去冗余统计 ──")
    rmdup_rows = []
    n_before = len(clusters)
    n_after = n_before
    n_removed = 0
    removed_bp = 0

    if rmdup_fa.is_file():
        rmdup_ids = set()
        rmdup_lens = []
        for rec in SeqIO.parse(str(rmdup_fa), "fasta"):
            rmdup_ids.add(rec.id)
            rmdup_lens.append(len(rec.seq))
        n_after = len(rmdup_ids)
        n_removed = n_before - n_after

        # 统计被移除的代表序列
        removed_lens = []
        for cid, members in clusters.items():
            rep = min(members, key=lambda x: len(x)) if members else ""
            if rep not in rmdup_ids and rep in centroids:
                removed_lens.append(centroids[rep])
        removed_bp = sum(removed_lens)

        rmdup_rows.append({
            "stage": "vclust 聚类后",
            "num_centroids": n_before,
            "total_bp": sum(centroids.values()),
            "n50": _n50(list(centroids.values())),
            "n90": _n90(list(centroids.values())),
        })
        kept_lens = [centroids[rid] for rid in rmdup_ids if rid in centroids]
        rmdup_rows.append({
            "stage": "rmdup 去冗余后",
            "num_centroids": n_after,
            "total_bp": sum(kept_lens),
            "n50": _n50(kept_lens),
            "n90": _n90(kept_lens),
        })
        rmdup_rows.append({
            "stage": "移除",
            "num_centroids": n_removed,
            "total_bp": removed_bp,
            "n50": _n50(removed_lens) if removed_lens else 0,
            "n90": _n90(removed_lens) if removed_lens else 0,
        })
        _write_table(rmdup_rows, out_dir / "rmdup_stats.tsv", "stage")
        print(f"  → {out_dir / 'rmdup_stats.tsv'} (移除 {n_removed} 条, {removed_bp:,} bp)")
    else:
        rmdup_rows.append({
            "stage": "vclust 聚类后 (未 rmdup)",
            "num_centroids": n_before,
            "total_bp": sum(centroids.values()),
            "n50": _n50(list(centroids.values())),
            "n90": _n90(list(centroids.values())),
        })
        _write_table(rmdup_rows, out_dir / "rmdup_stats.tsv", "stage")
        print(f"  → {out_dir / 'rmdup_stats.tsv'} (未运行 rmdup)")

    # ── 4. Flye 贡献统计 ──
    flye_centroids = [cid for cid in centroids if cid.startswith("flye_")]
    if flye_centroids:
        print("\n── 4. Flye 共组装贡献 ──")
        flye_rows = []
        for cid in flye_centroids:
            flye_rows.append({
                "centroid_id": cid,
                "length": centroids[cid],
            })

        flye_rows.append({})
        flye_rows.append({
            "centroid_id": f"Flye centroids: {len(flye_centroids)} 条",
            "length": sum(centroids[c] for c in flye_centroids),
        })
        normal_count = len(centroids) - len(flye_centroids)
        flye_rows.append({
            "centroid_id": f"正常 centroids: {normal_count} 条",
            "length": sum(centroids[c] for c in centroids if not c.startswith("flye_")),
        })
        flye_rows.append({
            "centroid_id": f"总计: {len(centroids)} 条",
            "length": sum(centroids.values()),
        })

        _write_table(flye_rows, out_dir / "flye_contribution.tsv", "centroid_id")
        print(f"  → {out_dir / 'flye_contribution.tsv'} (Flye {len(flye_centroids)}/{len(centroids)} 条, "
              f"{sum(centroids[c] for c in flye_centroids)/max(sum(centroids.values()),1)*100:.1f}%)")

    # ── 汇总 ──
    print(f"\n{'='*60}")
    print(f"分析完成 → {out_dir}")
    print(f"  sample_virus_counts.tsv      {len(rows)} 个样本/组")
    print(f"  cluster_sample_matrix.tsv    {len(matrix_rows)} 个 cluster")
    print(f"  cluster_summary_enriched.tsv {len(cluster_info)} 个 cluster")
    print(f"  rmdup_stats.tsv              rmdup 前后对比")
    if flye_centroids:
        print(f"  flye_contribution.tsv        Flye 贡献 {len(flye_centroids)} 条")


def _n50(lengths):
    if not lengths:
        return 0
    lengths = sorted(lengths, reverse=True)
    half = sum(lengths) / 2
    cum = 0
    for l in lengths:
        cum += l
        if cum >= half:
            return l
    return lengths[-1]


def _n90(lengths):
    if not lengths:
        return 0
    lengths = sorted(lengths, reverse=True)
    target = sum(lengths) * 0.9
    cum = 0
    for l in lengths:
        cum += l
        if cum >= target:
            return l
    return lengths[-1]


def _write_table(rows, out_path, sort_key=None):
    """写入 TSV，使用 polars/pandas/纯 Python 回退"""
    if HAS_POLARS:
        df = pl.DataFrame(rows)
        if sort_key and sort_key in df.columns:
            df = df.sort(sort_key)
        df.write_csv(str(out_path), separator="\t")
    elif HAS_PANDAS:
        df = pd.DataFrame(rows)
        if sort_key and sort_key in df.columns:
            df = df.sort_values(sort_key)
        df.to_csv(str(out_path), sep="\t", index=False)
    else:
        if not rows:
            with open(out_path, "w") as f:
                f.write("\n")
            return
        keys = list(rows[0].keys())
        with open(out_path, "w") as f:
            f.write("\t".join(keys) + "\n")
            for row in rows:
                f.write("\t".join(str(row.get(k, "")) for k in keys) + "\n")


if __name__ == "__main__":
    main()
