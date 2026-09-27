#!/usr/bin/env python3
"""
regenerate_coverage_plot.py — 独立重绘 Step 12 Coverage Visualization
====================================================================
从已有组装中间文件重新生成覆盖度可视化图，不重跑组装流程。

用法:
  python regenerate_coverage_plot.py \
    -r <参考 FASTA> \
    -d <组装样本目录 (含各 Step .fasta 文件)> \
    -o <输出目录>

示例:
  python regenerate_coverage_plot.py \
    -r ref_MW648525.1.ref.fasta \
    -d CRR527041.MW648525.1/ \
    -o CRR527041.MW648525.1/12.Coverage_Visualization/
"""

import argparse
import os
import sys
import tempfile
import subprocess
from pathlib import Path
import re

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
from Bio import SeqIO

# ── 从 virus-full.py 复制的核心函数 ──

def draw_alignment_block(ax, aln, rect_y, color):
    start, map_span = aln['start'], aln['map_span']
    ax.add_patch(patches.Rectangle((start, rect_y), map_span, 0.4,
        facecolor='#ff9999', edgecolor='none', alpha=0.9, zorder=2))
    for seg in aln['match_segments']:
        ax.add_patch(patches.Rectangle((seg['start'], rect_y), seg['len'], 0.4,
            facecolor=color, edgecolor='none', alpha=0.9, zorder=3))
    ax.add_patch(patches.Rectangle((start, rect_y), map_span, 0.4,
        facecolor='none', edgecolor='black', lw=1, zorder=4))


def find_n_regions(fasta_path, min_n_run=3):
    """从 FASTA 中提取连续 N 区域 [(start, end), ...]"""
    regions = []
    try:
        for rec in SeqIO.parse(fasta_path, "fasta"):
            seq = str(rec.seq).upper()
            in_n = False; n_start = 0
            for i, base in enumerate(seq):
                if base == 'N':
                    if not in_n:
                        n_start = i; in_n = True
                else:
                    if in_n and (i - n_start) >= min_n_run:
                        regions.append((n_start, i))
                    in_n = False
            if in_n and (len(seq) - n_start) >= min_n_run:
                regions.append((n_start, len(seq)))
            break
    except Exception:
        pass
    return regions


def draw_n_overlay(ax, n_regions, y_top, y_bottom, color='#222222', alpha=0.55):
    for n_start, n_end in n_regions:
        ax.add_patch(patches.Rectangle(
            (n_start, y_bottom), n_end - n_start, y_top - y_bottom,
            facecolor=color, edgecolor='none', alpha=alpha, zorder=6, hatch='////'))
    return len(n_regions)


# ── PAF 解析 ──

def run_minimap2(ref_fasta, query_fasta, paf_output, threads=4):
    cmd = f"minimap2 -t {threads} -c -x asm10 {ref_fasta} {query_fasta} > {paf_output} 2>/dev/null"
    subprocess.run(cmd, shell=True, check=True)


def parse_paf(paf_path):
    """与 virus-full.py parse_paf_with_cigar 逻辑一致"""
    alignments, ref_length = [], 0
    if not os.path.exists(paf_path): return alignments, ref_length
    with open(paf_path) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 11: continue
            q_name, q_len, strand = parts[0], int(parts[1]), parts[4]
            r_name, r_len, r_start, r_end = parts[5], int(parts[6]), int(parts[7]), int(parts[8])
            ref_length = max(ref_length, r_len)
            cigar_str = next((tag[5:] for tag in parts[12:] if tag.startswith("cg:Z:")), "")
            match_segments = []
            if cigar_str:
                curr_r = r_start
                for length_str, op in re.findall(r'(\d+)([MIDNSHP=X])', cigar_str):
                    length = int(length_str)
                    if op in ['M', '=', 'X']:
                        match_segments.append({'start': curr_r, 'len': length}); curr_r += length
                    elif op in ['D', 'N']: curr_r += length
            else:
                match_segments = [{'start': r_start, 'len': r_end - r_start}]
            alignments.append({
                'q_name': q_name, 'q_len': q_len, 'strand': strand,
                'start': r_start, 'end': r_end,
                'map_span': r_end - r_start, 'match_segments': match_segments,
            })
    alignments.sort(key=lambda x: x['start'])
    return alignments, ref_length


# ── 主绘图函数 ──

def plot_coverage_compact(tool_data, ref_length, out_prefix, n_regions_map=None):
    """n_regions_map: {step_name: [(start,end), ...]} 每个步骤各自的 N 区域"""
    colors = plt.cm.Set2.colors
    num_tools = len([t for t in tool_data.values() if t])
    fig, ax = plt.subplots(figsize=(16, max(5, num_tools * 1.8)))
    # 参考灰条 + 长度标注
    ax.add_patch(patches.Rectangle((0, 0), ref_length, 0.4,
        facecolor='lightgray', edgecolor='black', lw=1, zorder=2))
    ax.text(ref_length + ref_length * 0.005, 0.2, f'{ref_length} bp',
            va='center', fontsize=10, color='gray', fontstyle='italic')
    # 参考条上的 N（来自最终结果）
    ref_n = n_regions_map.get('__ref__', []) if n_regions_map else []
    if ref_n:
        draw_n_overlay(ax, ref_n, 0.4, 0, color='#222222', alpha=0.55)
    y_pos, yticks_pos, yticks_labels, tool_idx = 1.0, [0.2], ['Reference'], 0
    for tool_name, data in tool_data.items():
        alignments = data.get('alignments', data) if isinstance(data, dict) else data
        if not alignments: continue
        color = colors[tool_idx % len(colors)]; tool_idx += 1
        # 该步骤的背景区
        ax.add_patch(patches.Rectangle((-ref_length*0.05, y_pos - 0.1),
            ref_length*1.1, 0.6, facecolor='whitesmoke', edgecolor='none', alpha=0.3, zorder=1))
        for aln in alignments:
            draw_alignment_block(ax, aln, y_pos, color)
            if aln['map_span'] > ref_length * 0.08:
                short_name = aln['q_name'][:14] + ".." if len(aln['q_name']) > 16 else aln['q_name']
                ax.text(aln['start'] + aln['map_span']/2, y_pos + 0.2,
                        f"{short_name}\n{aln['q_len']}bp",
                        ha='center', va='center', fontsize=8,
                        color='black', fontweight='bold', zorder=5)
        # 该步骤自己的 N 区域（如果有）
        step_n = n_regions_map.get(tool_name, []) if n_regions_map else []
        if step_n:
            draw_n_overlay(ax, step_n, y_pos + 0.4, y_pos, color='#333333', alpha=0.5)
        yticks_pos.append(y_pos + 0.2); yticks_labels.append(tool_name); y_pos += 1.2
    ax.set_xlim(-ref_length * 0.15, ref_length * 1.08); ax.set_ylim(-0.5, y_pos)
    for s in ['top', 'right', 'left']: ax.spines[s].set_visible(False)
    ax.set_yticks(yticks_pos); ax.set_yticklabels(yticks_labels, fontsize=12, fontweight='bold')
    ax.set_xlabel('Genomic Coordinates (bp)', fontsize=14, fontweight='bold')
    ax.set_title('Genome Assembly Evolution (Compact View)', fontsize=18, fontweight='bold', pad=20)
    ax.grid(axis='x', linestyle='--', alpha=0.5, zorder=0)
    has_any_n = any(n_regions_map.values()) if n_regions_map else False
    legend_elements = [patches.Patch(color='#ff9999', label='Gap (Deletion)')]
    if has_any_n:
        legend_elements.append(patches.Patch(facecolor='#333333', alpha=0.5,
            hatch='////', label='N-Regions'))
    ax.legend(handles=legend_elements, loc='upper right', frameon=True)
    plt.tight_layout()
    plt.savefig(f"{out_prefix}_Compact.png", dpi=300)
    plt.savefig(f"{out_prefix}_Compact.pdf", dpi=300)
    plt.savefig(f"{out_prefix}_Compact.png", dpi=300)
    plt.close()


def plot_coverage_stacked(tool_data, ref_length, out_prefix, n_regions_map=None):
    colors = plt.cm.Set2.colors
    draw_info, y_cursor, yticks_pos, yticks_labels = {}, 1.0, [0.2], ['Reference']
    for tool_name, data in tool_data.items():
        alignments = data.get('alignments', data) if isinstance(data, dict) else data
        if not alignments: continue
        levels = []
        for aln in alignments:
            placed = False
            for i, end_pos in enumerate(levels):
                if aln['start'] > end_pos + (ref_length * 0.01):
                    aln['level'] = i; levels[i] = aln['end']; placed = True; break
            if not placed: aln['level'] = len(levels); levels.append(aln['end'])
        max_level = len(levels); tool_height = max_level * 0.6
        draw_info[tool_name] = {'alignments': alignments, 'base_y': y_cursor, 'max_level': max_level,
                                'height': tool_height}
        yticks_pos.append(y_cursor + (tool_height / 2) - 0.1)
        yticks_labels.append(tool_name)
        y_cursor += tool_height + 0.6
    fig, ax = plt.subplots(figsize=(16, max(5, y_cursor * 0.8)))
    ax.add_patch(patches.Rectangle((0, 0), ref_length, 0.4,
        facecolor='lightgray', edgecolor='black', lw=1, zorder=2))
    ax.text(ref_length + ref_length * 0.005, 0.2, f'{ref_length} bp',
            va='center', fontsize=10, color='gray', fontstyle='italic')
    ref_n = n_regions_map.get('__ref__', []) if n_regions_map else []
    if ref_n:
        draw_n_overlay(ax, ref_n, 0.4, 0, color='#222222', alpha=0.55)
    tool_idx = 0
    for tool_name, info in draw_info.items():
        color = colors[tool_idx % len(colors)]; tool_idx += 1; base_y = info['base_y']
        ax.add_patch(patches.Rectangle((-ref_length*0.05, base_y - 0.1),
            ref_length*1.1, info['max_level']*0.6, facecolor='whitesmoke',
            edgecolor='none', alpha=0.5, zorder=1))
        for aln in info['alignments']:
            rect_y = base_y + aln['level'] * 0.6
            draw_alignment_block(ax, aln, rect_y, color)
            if aln['map_span'] > ref_length * 0.04:
                short_name = aln['q_name'][:15] + ".." if len(aln['q_name']) > 17 else aln['q_name']
                ax.text(aln['start'] + aln['map_span']/2, rect_y + 0.2,
                        f"{short_name}\n{aln['q_len']}bp",
                        ha='center', va='center', fontsize=8,
                        color='black', fontweight='bold', zorder=5)
        # 每步骤自己的 N 区域
        step_n = n_regions_map.get(tool_name, []) if n_regions_map else []
        if step_n:
            draw_n_overlay(ax, step_n, base_y + info['height'], base_y - 0.1,
                           color='#333333', alpha=0.5)
    ax.set_xlim(-ref_length * 0.15, ref_length * 1.08); ax.set_ylim(-0.5, y_cursor)
    for s in ['top', 'right', 'left']: ax.spines[s].set_visible(False)
    ax.set_yticks(yticks_pos); ax.set_yticklabels(yticks_labels, fontsize=12, fontweight='bold')
    ax.set_xlabel('Genomic Coordinates (bp)', fontsize=14, fontweight='bold')
    ax.set_title('Genome Assembly Evolution (Stacked View)', fontsize=18, fontweight='bold', pad=20)
    ax.grid(axis='x', linestyle='--', alpha=0.5, zorder=0)
    has_any_n = any(n_regions_map.values()) if n_regions_map else False
    legend_elements = [patches.Patch(color='#ff9999', label='Gap (Deletion)')]
    if has_any_n:
        legend_elements.append(patches.Patch(facecolor='#333333', alpha=0.5,
            hatch='////', label='N-Regions'))
    ax.legend(handles=legend_elements, loc='upper right', frameon=True)
    plt.tight_layout()
    plt.savefig(f"{out_prefix}_Stacked.png", dpi=300)
    plt.savefig(f"{out_prefix}_Stacked.pdf", dpi=300)
    plt.savefig(f"{out_prefix}_Stacked.png", dpi=300)
    plt.close()


# ── 主函数 ──

def main():
    parser = argparse.ArgumentParser(description="独立重绘 Step 12 Coverage Visualization")
    parser.add_argument('-r', '--reference', required=True, help='参考基因组 FASTA')
    parser.add_argument('-d', '--sample-dir', required=True, help='组装样本目录 (含各 Step .fasta)')
    parser.add_argument('-o', '--outdir', required=True, help='输出目录')
    parser.add_argument('-t', '--threads', type=int, default=4, help='minimap2 线程数')
    args = parser.parse_args()

    sample_dir = Path(args.sample_dir)
    out_dir = Path(args.outdir)
    out_dir.mkdir(parents=True, exist_ok=True)
    ref_fasta = Path(args.reference)

    # 查找所有中间步骤 FASTA (与 virus-full.py valid_plots 对齐)
    step_files = {}
    step_patterns = [
        ('1.DeNovo_Cleaned',        sample_dir / '2.RefineC_Merge_Raw' / '2.refinec_merged_raw.fasta'),
        ('2.Ref_Merged_1',          sample_dir / '4.Ref_Merged_1' / '4.ref_merged_1.fasta'),
        ('3.Shiver_Cleanup',        sample_dir / '3.Shiver_Cleanup' / '3.cleaned_viral_contigs.fasta'),
        ('4.PVGA_Extension',        sample_dir / '5.PVGA_Extension' / '5.pvga_evaluated_cut.fasta'),
        ('5.Pre_Fusion_Merge',      sample_dir / '6.Pre_Fusion_Merge' / '6.pre_fusion_merged.fasta'),
        ('6.rmDup_Purification',    sample_dir / '7.Final_rmDup' / '7.rmDup_cleaned.fasta'),
        ('7.Ref_Merged_2',          sample_dir / '8.Ref_Merged_2' / '8.ref_merged_2.fasta'),
        ('8.Iterative_Consensus',   sample_dir / '9.Consensus_Polish' / '9.final_consensus.fasta'),
        ('9.Gap_Filled',            sample_dir / '10.Gap_Filling' / '10.gap_filled_final.fasta'),
        ('10.Ultimate_Result',      sample_dir / '11.Ultimate_Circular_Result.fasta'),
    ]

    for name, path in step_patterns:
        if path.exists() and path.stat().st_size > 0:
            step_files[name] = path

    if not step_files:
        print(f"错误: 在 {sample_dir} 中未找到任何 FASTA 文件")
        sys.exit(1)

    print(f"找到 {len(step_files)} 个中间步骤 FASTA:")
    for n, p in step_files.items():
        print(f"  {n}: {p.name}")

    # 计算每步的 N 区域
    n_regions_map = {}
    for step_name, fa_path in step_files.items():
        n_regs = find_n_regions(str(fa_path))
        if n_regs:
            n_regions_map[step_name] = n_regs
            print(f"  N 区域 ({step_name}): {sum(e-s for s,e in n_regs)} bp")

    # 参考条上的 N 取自最终结果
    for final_key in ['10.Ultimate_Result', '9.Gap_Filled', '8.Iterative_Consensus']:
        if final_key in n_regions_map:
            n_regions_map['__ref__'] = n_regions_map[final_key]
            break

    # minimap2 比对各步骤 vs 参考
    tool_data = {}
    global_ref_len = 0
    # 先用参考自身长度初始化
    ref_seq_len = len(str(next(SeqIO.parse(str(ref_fasta), "fasta")).seq))
    global_ref_len = ref_seq_len

    with tempfile.TemporaryDirectory() as tmpdir:
        for step_name, fa_path in step_files.items():
            paf_file = os.path.join(tmpdir, f"{step_name.replace('/', '_')}.paf")
            print(f"\n比对: {step_name} ...")
            run_minimap2(str(ref_fasta), str(fa_path), paf_file, args.threads)
            aln, r_len = parse_paf(paf_file)
            global_ref_len = max(global_ref_len, r_len)
            # 如果比对未覆盖全参考(末端N被软剪切)，扩展最后的比对块到参考全长
            if aln and aln[-1]['end'] < global_ref_len:
                last = aln[-1]
                extra = global_ref_len - last['end']
                last['end'] = global_ref_len
                last['map_span'] += extra
                if last['match_segments'] is None:
                    last['match_segments'] = [{'start': last['start'], 'len': last['map_span']}]
            tool_data[step_name] = {'alignments': aln, 'q_len': sum(1 for _ in SeqIO.parse(str(fa_path), "fasta"))}
            print(f"  → {len(aln)} 条比对, 参考长度={global_ref_len}")

    if global_ref_len == 0:
        print("错误: 无法确定参考长度")
        sys.exit(1)

    # 生成图
    sample_name = sample_dir.name
    plot_prefix = str(out_dir / f"{sample_name}_Coverage")
    print(f"\n生成 Coverage Visualization...")
    plot_coverage_compact(tool_data, global_ref_len, plot_prefix, n_regions_map=n_regions_map)
    plot_coverage_stacked(tool_data, global_ref_len, plot_prefix, n_regions_map=n_regions_map)
    print(f"完成 → {out_dir}/")
    print(f"  {plot_prefix}_Compact.png")
    print(f"  {plot_prefix}_Stacked.png")


if __name__ == '__main__':
    main()
