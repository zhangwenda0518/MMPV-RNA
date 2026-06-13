#!/usr/bin/env python3
"""
cluster_pipeline.py -- 病毒基因组聚类管道 v3.0
==============================================

流程:
  1. seqkit 长度过滤
  2. CD-HIT 参考引导预聚类 (可选, --ref-genomes)
     合并 contig + 参考 → vclust cd-hit → 拆分 known/novel
  3. vclust Leiden 聚类 (仅 novel 部分)
  4. 输出 centroids + per-cluster 拆分

后续: taxonomy → host → rescue_pipeline.py (三支路拯救)

依赖: vclust, seqkit
"""

import argparse, subprocess, sys, os, shutil
from pathlib import Path
from datetime import datetime
from collections import defaultdict

import pandas as pd
from Bio import SeqIO

try:
    import polars as pl
    HAS_POLARS = True
except ImportError:
    HAS_POLARS = False


# ══════════════════════════════════════════════════════════════
# 工具
# ══════════════════════════════════════════════════════════════

def run(cmd, desc="", timeout=None, check=True):
    label = desc or " ".join(str(c) for c in cmd)
    print(f"  [RUN] {label}", flush=True)
    try:
        return subprocess.run(cmd, capture_output=True, text=True, timeout=timeout, check=check)
    except subprocess.CalledProcessError as e:
        print(f"  [FAIL] 退出码 {e.returncode}")
        if e.stderr:
            print(f"        {e.stderr.strip()[:500]}")
        if check:
            sys.exit(1)
        return None


# ══════════════════════════════════════════════════════════════
# vclust
# ══════════════════════════════════════════════════════════════

def run_vclust_prefilter_align_cluster(fasta_in, out_dir, ani=0.95, qcov=0.85, threads=4, algorithm="leiden"):
    """vclust 三件套: prefilter → align → cluster。返回 clusters.tsv 路径"""
    d = Path(out_dir); d.mkdir(parents=True, exist_ok=True)
    pre, atsv, ctsv, ids = d / "vclust_prefilter.txt", d / "vclust_ani.tsv", d / "vclust_clusters.tsv", d / "vclust_ani.ids.tsv"
    run(["vclust", "prefilter", "-i", str(fasta_in), "-o", str(pre), "--min-ident", str(ani), "--threads", str(threads)], "vclust prefilter")
    run(["vclust", "align", "--filter", str(pre), "-i", str(fasta_in), "--out-ani", str(ani), "--out-qcov", str(qcov), "-o", str(atsv), "--out-aln", str(d / "vclust_ani.aln.tsv"), "--threads", str(threads)], "vclust align")
    run(["vclust", "cluster", "-i", str(atsv), "-o", str(ctsv), "--ids", str(ids), "--algorithm", algorithm, "--metric", "ani", "--ani", str(ani), "--qcov", str(qcov), "--out-repr"], "vclust cluster")
    return str(ctsv)


# ══════════════════════════════════════════════════════════════
# CD-HIT 参考引导预聚类 (vclust --algorithm cd-hit)
# ══════════════════════════════════════════════════════════════

TAG_REF = "ref|"
TAG_OUR = "our|"


def tag_and_merge_fasta(contig_fa, ref_fa, out_dir):
    """去重参考 → 给 ID 加前缀 → 合并: ref → ref|xxx, contig → our|xxx。返回 (merged_fa, ref_id_set)"""
    d = Path(out_dir); d.mkdir(parents=True, exist_ok=True)
    merged = d / "cdhit_combined.fasta"

    # 0. vclust deduplicate 去重参考基因组 (处理 ICTV/NCBI 重复)
    dedup_ref = d / "ref_dedup.fasta"
    if ref_fa and os.path.isfile(ref_fa):
        run(["vclust", "deduplicate", "-i", ref_fa, "-o", str(dedup_ref)], "vclust dedup 参考")
        ref_fa = str(dedup_ref)

    # 1. 标记 + 合并
    ref_ids = set()
    with open(merged, "w") as out:
        if ref_fa and os.path.isfile(ref_fa):
            for rec in SeqIO.parse(ref_fa, "fasta"):
                ref_ids.add(rec.id)
                out.write(f">{TAG_REF}{rec.id}\n{str(rec.seq)}\n")
        for rec in SeqIO.parse(contig_fa, "fasta"):
            out.write(f">{TAG_OUR}{rec.id}\n{str(rec.seq)}\n")

    n_ref = len(ref_ids)
    n_our = sum(1 for _ in SeqIO.parse(contig_fa, "fasta"))
    print(f"  CD-HIT 合并: {n_ref} 条参考 (去重后) + {n_our} 条 contig → {merged}")
    return str(merged), ref_ids


def strip_tags(seq_id):
    """去除 ref| 或 our| 前缀"""
    if seq_id.startswith(TAG_REF):
        return seq_id[len(TAG_REF):]
    if seq_id.startswith(TAG_OUR):
        return seq_id[len(TAG_OUR):]
    return seq_id


def split_cdhit_clusters(clusters_dict, ref_id_set):
    """
    根据代表序列前缀拆分 cd-hit 聚类结果。
    返回 (known_clusters, novel_clusters, association_map)

    known_clusters: 代表以 ref| 开头的簇
    novel_clusters: 代表以 our| 开头的簇
    association_map: {our_contig_id: ref_accession} 已知关联映射
    """
    known = {}
    novel = {}
    association = {}

    for cname, info in clusters_dict.items():
        rep = info["ref"]  # cd-hit 中最长序列为代表
        if rep.startswith(TAG_REF):
            known[cname] = info
            # 记录 contig→参考 关联
            ref_acc = strip_tags(rep)
            for member in info["members"]:
                if member.startswith(TAG_OUR):
                    association[strip_tags(member)] = ref_acc
        else:
            novel[cname] = info

    return known, novel, association


def run_cdhit_reference_clustering(input_fa, ref_fa, out_dir, ani=0.95, qcov=0.85, threads=4, min_length=500):
    """
    CD-HIT 参考引导预聚类:
      1. 合并 contig + 参考 (加 ref|/our| 前缀)
      2. vclust cd-hit 聚类
      3. 拆分 known / novel
      4. 写入 known 簇和 novel contig 子集
      返回 (novel_fa, known_centroids_fa, n_known_clusters, n_novel_contigs, association_map)
    """
    d = Path(out_dir) / "2_cdhit"; d.mkdir(parents=True, exist_ok=True)

    # 1. 合并
    merged_fa, ref_ids = tag_and_merge_fasta(input_fa, ref_fa, d)

    # 2. vclust cd-hit 聚类
    ctsv = run_vclust_prefilter_align_cluster(merged_fa, d / "vclust_cdhit", ani, qcov, threads, algorithm="cd-hit")
    all_clusters = parse_vclust_clusters(ctsv)
    print(f"  CD-HIT 聚类: {len(all_clusters)} 个簇")

    # 3. 拆分
    known, novel, association = split_cdhit_clusters(all_clusters, ref_ids)
    n_known = len(known)
    n_novel_contigs = sum(len(info["members"]) for info in novel.values())
    print(f"  已知关联簇: {n_known} (含参考基因组)")
    print(f"  新颖簇:     {len(novel)} ({n_novel_contigs} 条 contig)")

    # 4. 写入 known 簇
    known_dir = d / "known_clusters"; known_dir.mkdir(exist_ok=True)
    known_centroids = []
    for cname, info in known.items():
        # 收集簇内所有成员 (去前缀)
        records = []
        for member_id in info["members"]:
            clean_id = strip_tags(member_id)
            for rec in SeqIO.parse(merged_fa, "fasta"):
                if rec.id == member_id:
                    rec.id = clean_id; rec.description = ""
                    records.append(rec)
                    if member_id == info["ref"]:
                        known_centroids.append(rec)
                    break
        if records:
            cluster_fa = known_dir / f"{cname}.all.fasta"
            SeqIO.write(records, cluster_fa, "fasta")

    # known centroids
    known_centroids_fa = str(d / "known_centroids.fasta")
    if known_centroids:
        SeqIO.write(known_centroids, known_centroids_fa, "fasta")
        print(f"  known centroids: {len(known_centroids)} 条 → {known_centroids_fa}")

    # 写入 association 表
    assoc_tsv = d / "known_association.tsv"
    with open(assoc_tsv, "w") as af:
        af.write("contig_id\tref_accession\n")
        for contig_id, ref_acc in association.items():
            af.write(f"{contig_id}\t{ref_acc}\n")
    print(f"  association: {len(association)} 条 → {assoc_tsv}")

    # 5. 写入 novel contig 子集 (去前缀)
    novel_fa = str(d / "novel_contigs.fasta")
    novel_records = []
    novel_ids = set()
    for _, info in novel.items():
        for member_id in info["members"]:
            novel_ids.add(member_id)
    for rec in SeqIO.parse(merged_fa, "fasta"):
        if rec.id in novel_ids:
            rec.id = strip_tags(rec.id); rec.description = ""
            if len(rec.seq) >= min_length:
                novel_records.append(rec)
    SeqIO.write(novel_records, novel_fa, "fasta")
    n_written = sum(1 for _ in novel_records)
    print(f"  novel contig: {n_written} 条 (≥{min_length}bp) → {novel_fa}")

    return novel_fa, known_centroids_fa, n_known, n_written, association


# ══════════════════════════════════════════════════════════════

def parse_vclust_clusters(tsv):
    raw = defaultdict(list)
    with open(tsv) as f:
        f.readline()
        for line in f:
            p = line.strip().split()
            if len(p) >= 2:
                raw[p[1].strip()].append(p[0].strip())
    return {f"cluster_{cid}": {"ref": m[0], "members": list(set(m))} for cid, m in raw.items()}


# ══════════════════════════════════════════════════════════════
# 聚类后处理: 统计 + 物理拆分 (参照 run_cluster_pipeline.py)
# ══════════════════════════════════════════════════════════════

def load_fasta_info(fasta_path):
    """加载 FASTA, 返回 {id: {"header": str, "seq": str, "length": int}}"""
    info = {}
    cur_id, cur_seq, cur_header = None, [], ""
    with open(fasta_path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if cur_id:
                    s = "".join(cur_seq)
                    info[cur_id] = {"header": cur_header, "seq": s, "length": len(s)}
                cur_header, cur_id, cur_seq = line, line[1:].split()[0], []
            else:
                cur_seq.append(line)
        if cur_id:
            s = "".join(cur_seq)
            info[cur_id] = {"header": cur_header, "seq": s, "length": len(s)}
    return info


def extract_cluster_files(clusters, fasta_info, out_dir):
    """物理拆分: 全局代表序列 + per-cluster ref.fasta + per-cluster all.fasta"""
    d = Path(out_dir) / "split_fastas"; d.mkdir(parents=True, exist_ok=True)
    global_ref = Path(out_dir) / "all.cluster.ref.fasta"

    with open(global_ref, "w") as grf:
        for cname, info in clusters.items():
            ref_id = info["ref"]
            members = info["members"]

            # 全局代表合集
            if ref_id in fasta_info:
                grf.write(f"{fasta_info[ref_id]['header']}\n{fasta_info[ref_id]['seq']}\n")

            # cluster_X.ref.fasta
            with open(d / f"{cname}.ref.fasta", "w") as rf:
                if ref_id in fasta_info:
                    rf.write(f"{fasta_info[ref_id]['header']}\n{fasta_info[ref_id]['seq']}\n")

            # cluster_X.all.fasta
            with open(d / f"{cname}.all.fasta", "w") as af:
                for mid in members:
                    if mid in fasta_info:
                        af.write(f"{fasta_info[mid]['header']}\n{fasta_info[mid]['seq']}\n")

    print(f"  拆分完成: {len(clusters)} 簇 → {d}/")
    return str(global_ref), str(d)


def compute_cluster_stats(clusters, fasta_info, out_dir):
    """Polars 统计 + 全局摘要"""
    basic_stats = []
    for cname, info in clusters.items():
        ref_id = info["ref"]
        members = info["members"]
        ref_len = fasta_info.get(ref_id, {}).get("length", 0)
        other_members = [m for m in members if m != ref_id]
        other_lens = [fasta_info.get(m, {}).get("length", 0) for m in other_members if m in fasta_info]
        basic_stats.append({
            "Cluster_ID": cname, "Ref_ID": ref_id, "Ref_Length": ref_len,
            "Total_Size": len(members), "Other_Members_Count": len(other_members),
            "Other_Max_Len": max(other_lens) if other_lens else 0,
            "Other_Avg_Len": round(sum(other_lens) / len(other_lens), 2) if other_lens else 0.0,
            "Other_Min_Len": min(other_lens) if other_lens else 0,
        })

    if HAS_POLARS:
        df = pl.DataFrame(basic_stats, schema_overrides={"Other_Max_Len": pl.Int64, "Other_Min_Len": pl.Int64})
        df = df.with_columns(
            pl.col("Cluster_ID").str.extract(r"(\d+)").cast(pl.Int64).alias("sort_key")
        ).sort("sort_key").drop("sort_key")
        df.write_csv(os.path.join(out_dir, "cluster_summary.tsv"), separator="\t")
    else:
        pd.DataFrame(basic_stats).to_csv(os.path.join(out_dir, "cluster_summary.tsv"), sep="\t", index=False)

    total_clusters = len(basic_stats)
    total_input = len(fasta_info)
    sizes = [s["Total_Size"] for s in basic_stats]
    singletons = sum(1 for s in sizes if s == 1)
    summary = (
        "===========================================================\n"
        f"VCLUST 聚类统计概览\n"
        "===========================================================\n"
        f"  输入序列总数:    {total_input:,}\n"
        f"  总聚类簇数:      {total_clusters:,}\n"
        f"  整体去冗余率:    {(1 - total_clusters/total_input)*100:.2f}%\n"
        "  ---------------------------------------------------------\n"
        f"  单例簇数:        {singletons:,} ({singletons/total_clusters*100:.1f}%)\n"
        f"  多序列簇数:      {total_clusters - singletons:,}\n"
        f"  最大簇大小:      {max(sizes):,}\n"
        f"  平均簇大小:      {sum(sizes)/len(sizes):.2f}\n"
        "===========================================================\n"
    )
    print("\n" + summary)
    with open(os.path.join(out_dir, "global_summary.txt"), "w") as f:
        f.write(summary)
    return total_clusters, total_input, singletons


# ══════════════════════════════════════════════════════════════

# ══════════════════════════════════════════════════════════════
# 主流程
# ══════════════════════════════════════════════════════════════

def main():
    p = argparse.ArgumentParser(description="WVDB 病毒基因组三支路拯救管道 v2.0")
    p.add_argument("--input-fasta", "-i", required=True, help="输入 FASTA (规范命名后)")
    p.add_argument("--output-dir", "-o", required=True)
    p.add_argument("--fastq-dir", "-fq", default=".", help="原始 reads 目录 (聚类不需要, 保留兼容)")
    p.add_argument("--threads", "-t", type=int, default=64)
    p.add_argument("--min-length", type=int, default=500, help="病毒最小长度 bp (默认: 500)")
    p.add_argument("--ani", type=float, default=0.95)
    p.add_argument("--qcov", type=float, default=0.85)
    p.add_argument("--skip-vclust", action="store_true")
    p.add_argument("--vclust-cluster-file", default=None)
    p.add_argument("--stop-after-vclust", action="store_true",
                   help="仅运行到 vclust + 统计 + 拆分, 跳过三支路级联拯救 (供外部按宿主过滤后分批拯救)")
    p.add_argument("--ref-genomes", default=None,
                   help="ICTV/NCBI 参考基因组 FASTA (启用 CD-HIT 参考引导预聚类)")
    p.add_argument("--cdhit-ani", type=float, default=0.95,
                   help="CD-HIT 预聚类 ANI 阈值 (默认 0.95)")
    p.add_argument("--cdhit-qcov", type=float, default=0.85,
                   help="CD-HIT 预聚类 QCOV 阈值 (默认 0.85)")
    p.add_argument("--resume", action="store_true", help="断点续传")
    args = p.parse_args()

    start = datetime.now()
    print(f"开始: {start.strftime('%Y-%m-%d %H:%M:%S')}")

    out = Path(args.output_dir); out.mkdir(parents=True, exist_ok=True)

    # ── Step 1: seqkit 最小长度过滤 (快速, 仅此一步) ──
    # 注: genome_rmDuplicates 和 ID 截断已移除 — CD-HIT 参考引导预聚类天然处理冗余和碎片
    d1 = out / "1_seqkit"; d1.mkdir(exist_ok=True)
    virus_fa = d1 / "virus.candidate.fasta"
    id_map = {}  # 不再截断 ID, id_map 为空则 skip ID 还原

    if not args.resume or not virus_fa.is_file() or virus_fa.stat().st_size == 0:
        print("\n── Step 1: seqkit 长度过滤 ──")
        run(["seqkit", "seq", "-g", "-m", str(args.min_length), "-j", str(args.threads),
             args.input_fasta, "-o", str(virus_fa)],
            f"seqkit: ≥{args.min_length}bp")
        n_virus = sum(1 for _ in SeqIO.parse(virus_fa, "fasta"))
        print(f"  病毒基因组 (≥{args.min_length}bp): {n_virus:,}")
        input_fa = str(virus_fa)
    else:
        print(f"  [RESUME] Step 1 已有结果 → {virus_fa}")
        input_fa = str(virus_fa)

    # ── Step 2a: CD-HIT 参考引导预聚类 (可选) ──
    known_centroids_fa = None
    n_known_clusters = 0
    association_map = {}

    if args.ref_genomes and os.path.isfile(args.ref_genomes):
        print("\n── Step 2a: CD-HIT 参考引导预聚类 ──")
        novel_fa, known_centroids_fa, n_known_clusters, n_novel, association_map = \
            run_cdhit_reference_clustering(input_fa, args.ref_genomes, out,
                                           args.cdhit_ani, args.cdhit_qcov,
                                           args.threads, args.min_length)
        print(f"  CD-HIT 结果: {n_known_clusters} 已知簇 → CheckV, {n_novel} 新颖 contig → vclust Leiden")
        # 保存 association 表到 centroids
        assoc_out = out / "centroids" / "known_association.tsv"
        assoc_out.parent.mkdir(parents=True, exist_ok=True)
        src_assoc = Path(out) / "2_cdhit" / "known_association.tsv"
        if src_assoc.is_file():
            shutil.copy(src_assoc, assoc_out)
        # 此后 vclust 仅处理 novel 部分
        input_fa = novel_fa
    elif args.ref_genomes:
        print(f"\n  [WARN] --ref-genomes 指定的文件不存在: {args.ref_genomes}")

    # ── Step 2b: vclust Leiden (novel 部分或全部) ──
    d2 = out / "3_vclust"; log_dir_2 = d2 / "logs"; log_dir_2.mkdir(parents=True, exist_ok=True)
    ctsv = d2 / "vclust_clusters.tsv"
    global_ref_fa = d2 / "all.cluster.ref.fasta"

    if args.resume and ctsv.is_file() and global_ref_fa.is_file() and global_ref_fa.stat().st_size > 0:
        print(f"\n  [RESUME] Step 2 已有结果, 跳过 vclust")
        clusters = parse_vclust_clusters(ctsv)
        fasta_info = load_fasta_info(input_fa)
        tc, ts, sgl = compute_cluster_stats(clusters, fasta_info, d2)
    elif args.skip_vclust and args.vclust_cluster_file:
        ctsv = Path(args.vclust_cluster_file)
        print(f"\n── 复用 vclust: {ctsv}")
        clusters = parse_vclust_clusters(ctsv)
        fasta_info = load_fasta_info(input_fa)
        tc, ts, sgl = compute_cluster_stats(clusters, fasta_info, d2)
    else:
        print("\n── Step 2: vclust 初步聚类 ──")
        ctsv = run_vclust_prefilter_align_cluster(input_fa, d2, args.ani, args.qcov, args.threads)
        clusters = parse_vclust_clusters(ctsv)
        fasta_info = load_fasta_info(input_fa)
        tc, ts, sgl = compute_cluster_stats(clusters, fasta_info, d2)
        global_ref_fa = d2 / "all.cluster.ref.fasta"
        with open(global_ref_fa, "w") as grf:
            for cname, info in clusters.items():
                rid = info["ref"]
                if rid in fasta_info:
                    grf.write(f"{fasta_info[rid]['header']}\n{fasta_info[rid]['seq']}\n")
    print(f"  代表序列合集: {global_ref_fa}")

    if args.stop_after_vclust:
        # 产出 centroids (供 taxonomy/host 阶段读取)
        # CD-HIT known + vclust novel 合并, 统一走 taxonomy → host → rescue
        final_dir = out / "centroids"; final_dir.mkdir(parents=True, exist_ok=True)
        centroids_fa = final_dir / "final_centroids.fasta"
        known_id_file = final_dir / "known_ids.txt"

        n_known = 0
        known_ids = set()

        # 先写入 CD-HIT known centroids (已有完整参考基因组)
        if known_centroids_fa and os.path.isfile(known_centroids_fa):
            with open(centroids_fa, "w") as cf:
                for rec in SeqIO.parse(known_centroids_fa, "fasta"):
                    cf.write(f">{rec.id}\n{str(rec.seq)}\n")
                    known_ids.add(rec.id)
            n_known = len(known_ids)
            # 记录哪些 centroids 来自 CD-HIT known (供 rescue 阶段识别)
            with open(known_id_file, "w") as kf:
                for kid in sorted(known_ids):
                    kf.write(f"{kid}\n")

        # 追加 vclust Leiden centroids (novel contig 簇代表)
        n_novel = 0
        with open(centroids_fa, "a") as cf:
            for cname, info in clusters.items():
                rid = info["ref"]
                if rid in fasta_info:
                    cf.write(f"{fasta_info[rid]['header']}\n{fasta_info[rid]['seq']}\n")
                    n_novel += 1

        n_total = n_known + n_novel
        src = f" ({n_known} CD-HIT known + {n_novel} vclust novel)" if n_known else ""
        print(f"  代表序列合集: {centroids_fa} ({n_total} 条{src})")

        # 物理拆分: CD-HIT known 子文件也加入 split_fastas (供 rescue 识别)
        extracts_dir = d2 / "split_fastas"; extracts_dir.mkdir(parents=True, exist_ok=True)
        if n_known:
            cdhit_known_dir = Path(out) / "2_cdhit" / "known_clusters"
            if cdhit_known_dir.is_dir():
                for all_fa in cdhit_known_dir.glob("cluster_*.all.fasta"):
                    dest = extracts_dir / f"cdhit_{all_fa.name}"
                    shutil.copy(all_fa, dest)
            print(f"  CD-HIT known 拆分: {cdhit_known_dir} → {extracts_dir}")

        # 物理拆分 vclust novel
        extract_cluster_files(clusters, fasta_info, d2)
        print(f"  vclust novel 拆分: {d2}/split_fastas/")

        elapsed = (datetime.now() - start).total_seconds()
        print(f"\n{'=' * 60}")
        print(f"  [STOP] --stop-after-vclust: 已产出 clusters + centroids + 拆分文件, 跳过三支路拯救")
        if n_known:
            print(f"  CD-HIT known: {n_known} 条 (免 rescue, 已有完整参考)")
            print(f"  known IDs:    {known_id_file}")
        print(f"  vclust novel: {n_novel} 条 (后续进 rescue)")
        print(f"  clusters:     {ctsv}")
        print(f"  centroids:    {centroids_fa}")
        print(f"  耗时:         {elapsed / 60:.1f} min")
        print(f"{'=' * 60}")
        return

