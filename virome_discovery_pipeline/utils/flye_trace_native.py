#!/usr/bin/env python3
"""
Flye --subassemblies 输入→输出映射生成器
基于 read_alignment_dump + edge_mapping.tsv + assembly_info.txt 精确三级关联
(依赖修改后的 Flye 源码在 30-contigger 输出 edge_mapping.tsv)
"""

import argparse, os, re, sys
from collections import defaultdict

def n50(lengths):
    if not lengths: return 0
    lengths = sorted(lengths, reverse=True)
    half = sum(lengths) / 2
    cum = 0
    for l in lengths:
        cum += l
        if cum >= half: return l
    return lengths[-1]

def fmt(n):
    if n >= 1_000_000: return f"{n/1_000_000:.1f}M"
    elif n >= 1_000: return f"{n/1_000:.1f}K"
    return str(n)

def load_edge_mapping(edge_mapping_file):
    """edge_mapping.tsv: 重排edgeID -> [原始edgeID] 和反向 [原始edgeID] -> 重排edgeID"""
    edge_to_orig = {}   # "edge_458" -> {543, ...}
    orig_to_edge = {}   # "543" -> "edge_458"
    with open(edge_mapping_file) as f:
        for line in f:
            line = line.strip()
            if not line: continue
            parts = line.split("\t")
            if len(parts) < 2: continue
            edge_name = parts[0]           # "edge_458"
            orig_ids = parts[1].split(",") # ["543", ...]
            edge_to_orig[edge_name] = set(orig_ids)
            for oid in orig_ids:
                orig_to_edge[oid] = edge_name
    return edge_to_orig, orig_to_edge

def load_contig_paths(info_file):
    """assembly_info.txt -> contig -> [重排edgeID]"""
    contig_edges = {}
    with open(info_file) as f:
        for line in f:
            if line.startswith("contig_"):
                parts = line.strip().split("\t")
                cid, length = parts[0], parts[1]
                path = parts[7] if len(parts) > 7 else ""
                edges = []
                for token in path.split(","):
                    token = token.strip("*").lstrip("-")
                    if token.isdigit():
                        edges.append(token)
                contig_edges[cid] = (edges, int(length))
    return contig_edges

def main():
    parser = argparse.ArgumentParser(description="Flye --subassemblies 输入→输出映射生成器")
    parser.add_argument("-i", "--input", required=True, help="Flye 输出目录路径")
    parser.add_argument("-o", "--output", default="-", help="输出 TSV 文件路径")
    parser.add_argument("-s", "--summary", default="", help="统计摘要文件")
    parser.add_argument("--full", action="store_true", help="显示全部输入序列名称")
    parser.add_argument("-p", "--plot", default="", help="统计图表")

    args = parser.parse_args()
    flye_dir = args.input

    dump_file = os.path.join(flye_dir, "20-repeat", "read_alignment_dump")
    info_file = os.path.join(flye_dir, "assembly_info.txt")
    edge_mapping_file = os.path.join(flye_dir, "30-contigger", "edge_mapping.tsv")

    if not os.path.isfile(edge_mapping_file):
        print("错误: 找不到 edge_mapping.tsv，请用修改后的 Flye 重新运行 contigger 阶段", file=sys.stderr)
        sys.exit(1)

    # === 1. read_alignment_dump -> 输入contig -> {原始edge _id} ===
    # dump 第2列是 edgeId 的内部 _id (偶数=正向, 奇数=反向)
    read_to_orig_edges = defaultdict(set)
    read_lengths = {}
    with open(dump_file) as f:
        for line in f:
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 4 or cols[1] != "Aln":
                continue
            try:
                edge_id = int(cols[2])   # 内部 _id
            except ValueError:
                continue
            # 统一成正向 _id (偶数)
            fwd_id = edge_id - (edge_id % 2)
            # 提取 read 名称 (第3列的第一个 token, 去掉 +/- 前缀)
            read_tokens = cols[3].split()
            if not read_tokens:
                continue
            read_name = read_tokens[0]
            if read_name[0] in "+-":
                read_name = read_name[1:]
            # read 长度 (第3列的第4个 token, 即 read_len)
            read_len = 0
            if len(read_tokens) >= 4:
                try:
                    read_len = int(read_tokens[3])
                except ValueError:
                    pass
            read_to_orig_edges[read_name].add(str(fwd_id))
            if read_name not in read_lengths or read_len > read_lengths[read_name]:
                read_lengths[read_name] = read_len

    # === 2. edge_mapping: 原始edgeID -> 重排edgeID ===
    edge_to_orig, orig_to_edge = load_edge_mapping(edge_mapping_file)

    # === 3. assembly_info: contig -> [重排edgeID] ===
    contig_edges = load_contig_paths(info_file)

    # === 4. 三级关联: 输入contig -> 原始edge -> 重排edge -> contig ===
    # 反向: 重排edge -> {输入contig}
    edge_to_reads = defaultdict(set)
    unmapped = 0
    for read_name, orig_edges in read_to_orig_edges.items():
        for oid in orig_edges:
            if oid in orig_to_edge:
                renamed_edge = orig_to_edge[oid]
                edge_to_reads[renamed_edge].add(read_name)
            else:
                unmapped += 1

    # === 5. 输出映射表 ===
    out = open(args.output, "w") if args.output != "-" else sys.stdout
    header = ["output_contig", "length", "num_input_reads", "num_graph_edges", "input_reads"]
    out.write("\t".join(header) + "\n")

    output_lengths = []
    input_counts = []
    node_counts = []
    total_input_len = sum(read_lengths.values())
    input_lens = list(read_lengths.values())

    for cid, (edges, length) in sorted(contig_edges.items(), key=lambda x: int(x[0].split("_")[1])):
        all_reads = set()
        for e in edges:
            # e 是重排 edge ID（数字），对应 edge_name = "edge_" + e
            edge_name = f"edge_{e}"
            all_reads |= edge_to_reads.get(edge_name, set())
        
        output_lengths.append(length)
        node_counts.append(len(edges))
        input_counts.append(len(all_reads))
        
        if args.full or len(all_reads) <= 20:
            reads_str = "; ".join(sorted(all_reads))
        else:
            reads_str = "; ".join(sorted(all_reads)[:20]) + f" ... (+{len(all_reads)-20} more)"
        if not reads_str:
            reads_str = "(no reads mapped)"
        
        out.write("\t".join([cid, str(length), str(len(all_reads)), str(len(edges)), reads_str]) + "\n")

    if out is not sys.stdout:
        out.close()

    n_output = len(output_lengths)
    total_output_len = sum(output_lengths)

    # === 6. 摘要 ===
    summary = f"""============================================================
  Flye --subassemblies 组装统计摘要 (edge_mapping 精确关联)
============================================================
[输出]  Flye contigs: {sum(1 for _ in open(os.path.join(flye_dir, 'assembly.fasta')) if _.startswith('>'))}
  有映射 contigs: {n_output}  总碱基: {total_output_len:,} bp
[输入]  唯一源 contigs: {len(read_lengths)}
  总碱基: {total_input_len:,} bp
  未关联的原始edge引用: {unmapped}
"""
    print(summary, file=sys.stderr)
    if args.summary:
        with open(args.summary, "w") as sf:
            sf.write(summary)

if __name__ == "__main__":
    main()
