#!/usr/bin/env python3
"""
Flye --subassemblies 输入→输出映射生成器
从 Flye 内部文件 (read_alignment_dump + assembly_info.txt) 直接提取映射关系，无需 minimap2
"""

import argparse, os, re, sys
from collections import defaultdict

def n50(lengths):
    """计算 N50"""
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

def n90(lengths):
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

def fmt(n):
    """格式化数字"""
    if n >= 1_000_000:
        return f"{n/1_000_000:.1f}M"
    elif n >= 1_000:
        return f"{n/1_000:.1f}K"
    return str(n)

def _plot(input_lens, output_lengths, node_counts, input_counts,
          n_input, n_output, total_in, total_out,
          extended, merged, singleton, outfile):
    """生成统计图表"""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.ticker as ticker
    import numpy as np

    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.size": 9,
        "axes.titlesize": 11,
        "axes.labelsize": 9,
        "figure.dpi": 150,
    })

    fig, axes = plt.subplots(2, 3, figsize=(16, 10))
    ((ax1, ax2, ax3), (ax4, ax5, ax6)) = axes

    # --- ax1: 长度分布对比 ---
    bins = np.logspace(np.log10(200), np.log10(max(max(input_lens), max(output_lengths)) + 1000), 60)
    ax1.hist(input_lens, bins=bins, alpha=0.6, label=f"Input (n={n_input})", color="steelblue", edgecolor="white", linewidth=0.3)
    ax1.hist(output_lengths, bins=bins, alpha=0.8, label=f"Output (n={n_output})", color="darkorange", edgecolor="white", linewidth=0.3)
    ax1.set_xscale("log")
    ax1.set_xlabel("Length (bp)")
    ax1.set_ylabel("Count")
    ax1.set_title("Contig Length Distribution")
    ax1.legend(fontsize=8)
    ax1.set_xlim(200, None)

    # --- ax2: 输入参与度饼图 ---
    bins_input = [(1,1), (2,5), (6,20), (21,100), (101,99999)]
    labels_input = ["1 read", "2-5", "6-20", "21-100", ">100"]
    sizes = [sum(1 for n in input_counts if lo <= n <= hi) for lo, hi in bins_input]
    colors = ["#e8e8e8", "#b3cde3", "#6497b1", "#ffb347", "#e76f51"]
    wedges, texts, autotexts = ax2.pie(sizes, labels=labels_input, colors=colors,
                                         autopct="%1.1f%%", startangle=90,
                                         textprops={"fontsize": 8})
    ax2.set_title("Output Contigs by # Input Reads")

    # --- ax3: 散点图 - 输出长度 vs 输入数量 ---
    sizes_pt = [min(max(n, 1) ** 0.5 * 3, 80) for n in input_counts]
    colors_scatter = ["#e76f51" if n > 1 else "#6497b1" for n in node_counts]
    ax3.scatter(input_counts, output_lengths, s=sizes_pt, c=colors_scatter,
                alpha=0.5, edgecolors="none")
    ax3.set_xscale("log")
    ax3.set_yscale("log")
    ax3.set_xlabel("# Input Reads Merged")
    ax3.set_ylabel("Output Length (bp)")
    ax3.set_title("Output Length vs Merge Count")
    ax3.axhline(y=1000, color="gray", linestyle="--", linewidth=0.5, alpha=0.4)
    # legend
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], marker="o", color="w", markerfacecolor="#e76f51", markersize=8, label="Multi-node (extended)"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor="#6497b1", markersize=8, label="Single node"),
    ]
    ax3.legend(handles=legend_elements, fontsize=7, loc="lower right")

    # --- ax4: 图节点复杂度分布 ---
    bins_node = [(1,1), (2,2), (3,5), (6,10), (11,50), (51,999)]
    labels_node = ["1", "2", "3-5", "6-10", "11-50", ">50"]
    counts_node = [sum(1 for n in node_counts if lo <= n <= hi) for lo, hi in bins_node]
    colors_node = ["#e8e8e8", "#b3cde3", "#6497b1", "#ffb347", "#e76f51", "#c0392b"]
    bars = ax4.bar(range(len(labels_node)), counts_node, color=colors_node, edgecolor="white", linewidth=0.5)
    ax4.set_xticks(range(len(labels_node)))
    ax4.set_xticklabels(labels_node)
    ax4.set_xlabel("Graph Nodes per Contig")
    ax4.set_ylabel("Count")
    ax4.set_title("Graph Complexity Distribution")
    for bar, cnt in zip(bars, counts_node):
        if cnt > 0:
            ax4.text(bar.get_x() + bar.get_width()/2, bar.get_height() + max(counts_node)*0.02,
                     str(cnt), ha="center", fontsize=8)

    # --- ax5: 压缩效果对比 ---
    metrics = ["Contigs", "Bases (Mb)", "N50 (Kb)"]
    inp_vals = [n_input, total_in / 1e6, n50(input_lens) / 1000 if input_lens else 0]
    out_vals = [n_output, total_out / 1e6, n50(output_lengths) / 1000 if output_lengths else 0]
    x = np.arange(len(metrics))
    w = 0.35
    bars1 = ax5.bar(x - w/2, inp_vals, w, label="Input", color="steelblue", edgecolor="white", linewidth=0.5)
    bars2 = ax5.bar(x + w/2, out_vals, w, label="Output", color="darkorange", edgecolor="white", linewidth=0.5)
    ax5.set_xticks(x)
    ax5.set_xticklabels(metrics)
    ax5.set_title("Assembly Compression")
    ax5.legend(fontsize=8)
    for bar, val in zip(bars1, inp_vals):
        ax5.text(bar.get_x() + bar.get_width()/2, bar.get_height() + max(inp_vals + out_vals)*0.02,
                 f"{val:.1f}" if val < 100 else str(int(val)), ha="center", fontsize=8)
    for bar, val in zip(bars2, out_vals):
        ax5.text(bar.get_x() + bar.get_width()/2, bar.get_height() + max(inp_vals + out_vals)*0.02,
                 f"{val:.1f}" if val < 100 else str(int(val)), ha="center", fontsize=8)

    # --- ax6: 延伸 vs 去冗余 分类统计 ---
    categories = ["Extended\n(multi-node)", "Merged\n(multi-input)", "Singleton\n(1→1 pass)"]
    values = [extended, merged, singleton]
    colors_cat = ["#e76f51", "#ffb347", "#6497b1"]
    bars = ax6.bar(range(len(categories)), values, color=colors_cat, edgecolor="white", linewidth=0.5)
    ax6.set_xticks(range(len(categories)))
    ax6.set_xticklabels(categories, fontsize=8)
    ax6.set_ylabel("Count")
    ax6.set_title(f"Extension vs Dedup (total={n_output})")
    for bar, val in zip(bars, values):
        ax6.text(bar.get_x() + bar.get_width()/2, bar.get_height() + max(values)*0.02,
                 f"{val}\n({val/n_output*100:.1f}%)", ha="center", fontsize=8)

    fig.suptitle(f"Flye --subassemblies Assembly Report\nIn: {n_input:,} reads ({total_in/1e6:.1f} Mb)  →  Out: {n_output:,} contigs ({total_out/1e6:.1f} Mb)",
                 fontsize=13, fontweight="bold", y=1.01)
    plt.tight_layout()
    fig.savefig(outfile, dpi=150, bbox_inches="tight", facecolor="white")
    plt.close()
    print(f"图表已保存: {outfile}", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(
        description="从 Flye --subassemblies 输出目录提取 输入序列→输出contig 的映射关系",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例:
  python flye_trace_native.py -i fly.fasta2/ -o mapping.tsv
  python flye_trace_native.py -i fly.fasta2/ -o mapping.tsv --full --summary summary.tsv
  python flye_trace_native.py -i fly.fasta2/              # 直接打印到终端
        """
    )
    parser.add_argument("-i", "--input", required=True, metavar="DIR",
                        help="Flye 输出目录路径")
    parser.add_argument("-o", "--output", default="-", metavar="FILE",
                        help="输出 TSV 文件路径 (默认: 打印到终端)")
    parser.add_argument("-s", "--summary", default="", metavar="FILE",
                        help="额外输出统计摘要文件")
    parser.add_argument("--full", action="store_true",
                        help="显示全部输入序列名称 (默认只显示前 20 条)")
    parser.add_argument("--no-reads", action="store_true",
                        help="不显示输入序列名称列 (精简输出)")
    parser.add_argument("-p", "--plot", default="", metavar="FILE",
                        help="生成统计图表 (png/pdf/svg)")
    args = parser.parse_args()

    flye_dir = args.input
    out = open(args.output, "w", encoding="utf-8") if args.output != "-" else sys.stdout

    # === 检查必要文件 ===
    dump_file = os.path.join(flye_dir, "20-repeat", "read_alignment_dump")
    info_file = os.path.join(flye_dir, "assembly_info.txt")

    for fname, fpath in [("read_alignment_dump", dump_file), ("assembly_info.txt", info_file)]:
        if not os.path.exists(fpath):
            print(f"错误: 找不到 {fname}，请确认 -i 指向的是 Flye 输出目录", file=sys.stderr)
            sys.exit(1)

    # === 1. read_alignment_dump → input_read → {edge_id}, 同时收集read长度 ===
    read_to_edges = defaultdict(set)
    read_lengths = {}  # read_name → max_length (取最大值，因为一条read可能有多条alignment)

    with open(dump_file) as f:
        for line in f:
            m = re.match(
                r'\tAln\s+\d+\s+([+-])(\S+)\s+\d+\s+\d+\s+(\d+)\s+[+-]edge_(\d+)_\d+_[+-]disjointig_\d+',
                line
            )
            if m:
                read_name = m.group(2)
                read_len = int(m.group(3))
                edge_id = m.group(4)
                read_to_edges[read_name].add(edge_id)
                if read_name not in read_lengths or read_len > read_lengths[read_name]:
                    read_lengths[read_name] = read_len

    # === 2. assembly_info.txt → final_contig → {edge_ids} ===
    contig_edges = {}
    contig_meta = {}  # cid → (circ, repeat, mult) from assembly_info

    with open(info_file) as f:
        for line in f:
            if line.startswith("contig_"):
                parts = line.strip().split("\t")
                cid, length, cov, circ, repeat, mult, alt, path = parts
                edges = set()
                for token in path.split(","):
                    token = token.strip("*").lstrip("-")
                    if token.isdigit():
                        edges.add(token)
                contig_edges[cid] = (edges, int(length))
                contig_meta[cid] = (circ, repeat, int(mult))

    # === 3. 反向索引 ===
    edge_to_reads = defaultdict(set)
    for read_name, edges in read_to_edges.items():
        for e in edges:
            edge_to_reads[e].add(read_name)

    # === 4. 收集统计数据 ===
    output_lengths = []
    node_counts = []
    input_counts = []
    extended_contigs = 0    # 多节点 = 确实延长
    merged_contigs = 0      # 多输入 = 多个碎片合并
    singleton_contigs = 0   # 1输入, 1节点
    circular_contigs = 0
    repeat_contigs = 0

    total_input_len = sum(read_lengths.values())
    input_lens = list(read_lengths.values())

    for cid, (edges, length) in contig_edges.items():
        all_reads = set()
        for e in edges:
            all_reads |= edge_to_reads.get(e, set())

        if not all_reads and not edges:
            continue

        output_lengths.append(length)
        n_nodes = len(edges)
        n_inputs = len(all_reads)
        node_counts.append(n_nodes)
        input_counts.append(n_inputs)

        if n_nodes > 1:
            extended_contigs += 1
        if n_inputs > 1:
            merged_contigs += 1
        if n_inputs == 1 and n_nodes == 1:
            singleton_contigs += 1

        circ, repeat, _ = contig_meta.get(cid, ("N", "N", 1))
        if circ == "Y":
            circular_contigs += 1
        if repeat == "Y":
            repeat_contigs += 1

    n_output = len(output_lengths)
    n_mapped_input = len(read_lengths)
    n_total_reads_mapped = len(read_to_edges)
    total_output_len = sum(output_lengths)

    # === 5. 输出映射表 ===
    header = ["output_contig", "length", "num_input_reads", "num_graph_edges"]
    if not args.no_reads:
        header.append("input_reads")
    header.append("graph_edges")
    out.write("\t".join(header) + "\n")

    for cid, (edges, length) in sorted(contig_edges.items(), key=lambda x: int(x[0].split("_")[1])):
        all_reads = set()
        for e in edges:
            all_reads |= edge_to_reads.get(e, set())

        if not all_reads and not edges:
            continue

        if not args.no_reads:
            if args.full or len(all_reads) <= 20:
                reads_str = "; ".join(sorted(all_reads))
            else:
                reads_str = "; ".join(sorted(all_reads)[:20])
                reads_str += f" ... (+{len(all_reads)-20} more)"
            if not reads_str:
                reads_str = "(no reads mapped)"

        edges_str = ",".join(sorted(edges, key=int)) if edges else "(none)"
        row = [cid, str(length), str(len(all_reads)), str(len(edges))]
        if not args.no_reads:
            row.append(reads_str)
        row.append(edges_str)
        out.write("\t".join(row) + "\n")

    if out is not sys.stdout:
        out.close()

    # === 6. 绘图 ===
    if args.plot:
        _plot(input_lens, output_lengths, node_counts, input_counts,
              n_mapped_input, n_output, total_input_len, total_output_len,
              extended_contigs, merged_contigs, singleton_contigs,
              args.plot)

    # === 7. 输出统计摘要 ===
    summary_lines = [
        "=" * 60,
        "  Flye --subassemblies 组装统计摘要",
        "=" * 60,
        "",
        "[输入序列]",
        f"  总数:        {n_mapped_input}",
        f"  总碱基:      {fmt(total_input_len)} bp",
        f"  最大长度:    {max(input_lens)} bp",
        f"  N50 / N90:   {n50(input_lens)} / {n90(input_lens)} bp",
        f"  平均长度:    {int(sum(input_lens)/len(input_lens))} bp",
        "",
        "[输出 contig]",
        f"  总数:        {n_output}",
        f"  总碱基:      {fmt(total_output_len)} bp",
        f"  最大长度:    {max(output_lengths)} bp",
        f"  N50 / N90:   {n50(output_lengths)} / {n90(output_lengths)} bp",
        f"  平均长度:    {int(total_output_len/n_output)} bp",
        "",
        "[压缩效果]",
        f"  contig 数:   {n_mapped_input} → {n_output}  (压缩 {n_mapped_input/n_output:.1f}x)",
        f"  碱基数:      {fmt(total_input_len)} → {fmt(total_output_len)} bp  ({total_output_len/total_input_len*100:.1f}%)",
        f"  N50 提升:    {n50(input_lens)} → {n50(output_lengths)} bp  ({n50(output_lengths)/max(1,n50(input_lens)):.1f}x)",
        "",
        "[延伸 vs 去冗余]",
        f"  确实延伸:    {extended_contigs} 条 contig (多节点串联, 占 {extended_contigs/n_output*100:.1f}%)",
        f"  碎片合并:    {merged_contigs} 条 contig (多输入合并, 占 {merged_contigs/n_output*100:.1f}%)",
        f"  单纯保留:    {singleton_contigs} 条 contig (1→1 直接通过)",
        f"  环状:        {circular_contigs} 条",
        f"  含重复:      {repeat_contigs} 条",
        "",
        "[图复杂度分布]",
    ]

    # 节点数分布
    bins = [(1, 1), (2, 2), (3, 5), (6, 10), (11, 50), (51, 9999)]
    labels = ["1 节点 (未延伸)", "2 节点", "3-5 节点", "6-10 节点", "11-50 节点", ">50 节点"]
    for (lo, hi), label in zip(bins, labels):
        cnt = sum(1 for n in node_counts if lo <= n <= hi)
        if cnt > 0:
            summary_lines.append(f"  {label:<18} {cnt:>5} 条")

    summary_lines += [
        "",
        "[输入参与度分布]",
    ]
    bins = [(1, 1), (2, 5), (6, 20), (21, 100), (101, 9999)]
    labels = ["1 条输入", "2-5 条", "6-20 条", "21-100 条", ">100 条"]
    for (lo, hi), label in zip(bins, labels):
        cnt = sum(1 for n in input_counts if lo <= n <= hi)
        if cnt > 0:
            summary_lines.append(f"  {label:<18} {cnt:>5} 条 contig")

    summary_lines += [
        "",
        "[Top 10 合并最多的 contig]",
    ]

    top_merged = sorted(contig_edges.items(),
                        key=lambda x: len(set().union(*[edge_to_reads.get(e, set()) for e in x[1][0]])) if x[1][0] else 0,
                        reverse=True)[:10]
    for cid, (edges, length) in top_merged:
        n_inputs = len(set().union(*[edge_to_reads.get(e, set()) for e in edges])) if edges else 0
        summary_lines.append(f"  {cid:<16} {length:>6} bp  ← {n_inputs:>4} 条输入碎片")

    summary_lines += [
        "",
        "=" * 60,
    ]

    summary_text = "\n".join(summary_lines)

    if args.summary:
        with open(args.summary, "w", encoding="utf-8") as sf:
            sf.write(summary_text + "\n")
        print(f"统计摘要已写入: {args.summary}", file=sys.stderr)

    print(summary_text, file=sys.stderr)
    print(f"\n输出文件: {args.output if args.output != '-' else '(stdout)'}", file=sys.stderr)

if __name__ == "__main__":
    main()
