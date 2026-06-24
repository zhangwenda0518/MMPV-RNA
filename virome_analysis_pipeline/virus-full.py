#!/usr/bin/env python3
"""
OmniVirusAssembler - 生产级最终版本 (Production Edition V9.0 - Ultimate)
架构升级：
- [逻辑修正] 升级“提前交卷”为“抛光快车道”：即使早期骨架完美，也强制执行 Step 9 迭代抛光与 Step 11 环化检测，确保单碱基精度！
- [V9.0] 融合 Step 9 Reads级别迭代抛光出图功能，自动解析 GenBank (.gb) 绘制基因轨道与测序深度双轨图。
- [V9.0] 彻底解除 Shiver-like 过滤的基础长度限制，保留任何有效的同源碎片。
- [V9.0] 引入强大的自动化中间临时文件清理机制 (Cleanup Temp Files) 释放磁盘空间。
- [V8.7] 新增 Step 10：双引擎 Gap Filling (gmcloser 结构填补 + abyss-sealer Reads布隆填补)。
- [V8.7] 修复了 argparse 参数组丢失导致 mafft_args 报错的问题。
- [V8.6] 新增 PVGA 假连接熔断机制 (评估并斩断 >100bp 的错误延伸 Deletion)。
- [V8.6] 废弃 RagTag，启用全新 12 步精炼流水线。
- [V8.6] 架构重排：提前 rmDup 净化骨架物料 -> 纯净 Divine Fusion -> 后置 Consensus 全局迭代抛光。
- [V8.5] Divine Fusion: 100% 内存复现 SHIVER 算法，打造绝对实心骨架。
- [V8.5] 环化闭合检测: 精准切割，附带 [Circular=True] 标记。
"""

import sys
import os
import argparse
import subprocess
import shutil
import re
import math
import gzip
import shlex
import logging
import tempfile
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed

SCRIPT_DIR = Path(__file__).resolve().parent

try:
    from tqdm import tqdm
except ImportError:
    sys.exit("[致命错误] 未安装 tqdm。请运行 'pip install tqdm'")

try:
    import pandas as pd
    import numpy as np
    from Bio import SeqIO
    ADVANCED_STATS = True
except ImportError:
    ADVANCED_STATS = False

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import matplotlib.patches as patches
    import matplotlib.ticker as ticker
    CAN_PLOT = True
except ImportError:
    CAN_PLOT = False

# =====================================================================
# 日志系统配置
# =====================================================================
class TqdmLoggingHandler(logging.Handler):
    def __init__(self, level=logging.NOTSET):
        super().__init__(level)
    def emit(self, record):
        try:
            msg = self.format(record)
            tqdm.write(msg)
            self.flush()
        except Exception:
            self.handleError(record)

def setup_logger(out_dir, dry_run=False):
    logger = logging.getLogger("OmniVirusAssembler")
    logger.setLevel(logging.DEBUG)
    logger.handlers.clear() 
  
    formatter_console = logging.Formatter('%(asctime)s | %(levelname)-7s | %(message)s', datefmt='%H:%M:%S')
    formatter_file = logging.Formatter('%(asctime)s | %(levelname)-7s | %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
  
    ch = TqdmLoggingHandler()
    ch.setLevel(logging.INFO) 
    ch.setFormatter(formatter_console)
    logger.addHandler(ch)
  
    if not dry_run and out_dir:
        log_dir = os.path.join(out_dir, "logs")
        os.makedirs(log_dir, exist_ok=True)
        log_file = os.path.join(log_dir, "pipeline_detailed.log")
        fh = logging.FileHandler(log_file, mode='a')
        fh.setLevel(logging.DEBUG) 
        fh.setFormatter(formatter_file)
        logger.addHandler(fh)
    return logger

# =====================================================================
# 全局变量与底层安全辅助函数
# =====================================================================
VALID_EXTENSIONS = ['.fastq', '.fq', '.fastq.gz', '.fq.gz', '.fasta', '.fa', '.fasta.gz', '.fa.gz']

def is_valid(fp):
    if not fp: return False
    p = Path(fp)
    return p.exists() and p.stat().st_size > 0

def check_dependencies(tools):
    return [tool for tool in tools if shutil.which(tool) is None]

def is_fastq(filepath):
    return any(str(filepath).lower().endswith(ext) for ext in ['.fq', '.fastq', '.fq.gz', '.fastq.gz'])

def smart_open(filename, mode='rt'):
    if str(filename).lower().endswith('.gz'): return gzip.open(filename, mode, encoding='utf-8')
    return open(filename, mode, encoding='utf-8')

def safe_concat_to_fasta(input_files, output_file):
    with open(output_file, 'wt', encoding='utf-8') as out_f:
        for f in input_files:
            f_str = str(f)
            if f_str.lower().endswith('.gz'):
                with gzip.open(f_str, 'rt', encoding='utf-8') as in_f: shutil.copyfileobj(in_f, out_f)
            else:
                with open(f_str, 'rt', encoding='utf-8') as in_f: shutil.copyfileobj(in_f, out_f)

def extract_and_move_fasta(source_file: Path, target_file: Path):
    try:
        if str(source_file).lower().endswith('.gz'):
            with gzip.open(source_file, 'rt', encoding='utf-8') as f_in, open(target_file, 'wt', encoding='utf-8') as f_out: shutil.copyfileobj(f_in, f_out)
        else: shutil.copy2(source_file, target_file)
        return True
    except Exception: return False

def run_cmd(cmd, cwd=None, log_file=None, logger=None, sample_name="SYS", ignore_error=False, use_bash=False):
    is_str_cmd = isinstance(cmd, str) or use_bash
    if use_bash and isinstance(cmd, list): cmd = " ".join([str(x) for x in cmd])
    cmd_for_log = cmd if is_str_cmd else ' '.join(cmd)
  
    if logger and not ignore_error: logger.info(f"[{sample_name}] [CMD] {cmd_for_log}")
    try:
        kwargs = {'cwd': cwd, 'text': True, 'shell': is_str_cmd}
        if use_bash: kwargs['executable'] = '/bin/bash'
      
        if log_file:
            with open(log_file, 'w') as f:
                kwargs['stdout'] = f
                kwargs['stderr'] = subprocess.STDOUT
                proc = subprocess.run(cmd, **kwargs)
                if proc.returncode != 0 and logger and not ignore_error: logger.debug(f"[{sample_name}] [CMD FAILED] 日志: {log_file}")
                return proc.returncode == 0
        else:
            kwargs['capture_output'] = True
            proc = subprocess.run(cmd, **kwargs)
            if proc.returncode != 0 and logger and not ignore_error: logger.debug(f"[{sample_name}] [CMD ERROR] {proc.stderr.strip()}")
            return proc.returncode == 0
    except Exception as e:
        if logger and not ignore_error: logger.error(f"[{sample_name}] [SYS ERROR] {str(e)}")
        return False

def reverse_complement(seq):
    trans = str.maketrans('ATCGNacgnt', 'TAGCNtgcna')
    return seq.translate(trans)[::-1]

def read_first_sequence(fasta_file):
    if not fasta_file or not os.path.exists(fasta_file): return "", ""
    seq_name, seq_data = "", []
    with smart_open(fasta_file, 'rt') as f:
        for line in f:
            line = line.strip()
            if not line: continue
            if line.startswith(">"):
                if seq_name: break 
                seq_name = line[1:].split()[0]
            else: seq_data.append(line)
    return seq_name, "".join(seq_data).upper()

def parse_fasta_string(fasta_str):
    seqs, name, seq_data = {}, "", []
    for line in fasta_str.splitlines():
        line = line.strip()
        if not line: continue
        if line.startswith('>'):
            if name: seqs[name] = "".join(seq_data).upper()
            name = line[1:].split()[0]
            seq_data = []
        else: seq_data.append(line)
    if name: seqs[name] = "".join(seq_data).upper()
    return seqs

def write_fasta(file_path, seq_name, sequence):
    sequence = sequence.upper()
    with open(file_path, 'w') as f:
        f.write(f">{seq_name}\n")
        for i in range(0, len(sequence), 80): f.write(sequence[i:i+80] + "\n")

# =====================================================================
# 核心算法模块 (PVGA评估, Divine Fusion, Circularity)
# =====================================================================

def evaluate_and_split_pvga(ref_fasta, pvga_raw, pvga_out, threads, gap_size=100, min_len=200, logger=None, sample="Sample"):
    if not is_valid(pvga_raw):
        open(pvga_out, 'a').close(); return False
    logger.info(f"[{sample}] 评估 PVGA 延伸结果，检测超大 Gap (>{gap_size}bp) 并强制斩断...")
    ref_name, ref_seq = read_first_sequence(ref_fasta)
    combined_fasta = f">{ref_name}\n{ref_seq}\n"
    pvga_seqs = parse_fasta_string(Path(pvga_raw).read_text())
    if not pvga_seqs:
        open(pvga_out, 'a').close(); return False
    for k, v in pvga_seqs.items(): combined_fasta += f">{k}\n{v}\n"

    mafft_cmd = ['mafft', '--thread', str(threads), '--quiet', '--auto', '-']
    res = subprocess.run(mafft_cmd, input=combined_fasta, capture_output=True, text=True)
    if res.returncode != 0:
        logger.warning(f"[{sample}] PVGA 评估比划失败，继承原序列。")
        shutil.copy(pvga_raw, pvga_out); return False

    aln_seqs = parse_fasta_string(res.stdout)
    gap_pat = re.compile(r'-{%d,}' % gap_size)
    split_count = 0
    with open(pvga_out, 'w') as out_f:
        for name, a_seq in aln_seqs.items():
            if name == ref_name: continue
            frag_idx = 1
            for frag in gap_pat.split(a_seq):
                clean_frag = frag.replace('-', '').replace('N', '').upper()
                if len(clean_frag) >= min_len:
                    out_f.write(f">{name}_cut{frag_idx}\n{clean_frag}\n"); frag_idx += 1
            if frag_idx > 2: split_count += 1
    if split_count > 0 and logger: logger.info(f"[{sample}] 成功打断 {split_count} 条跨越巨大 Deletion 的 PVGA 错误连接序列！")
    return True

def divine_fusion_shiver_style(ref_file, contig_file, pvga_file, out_file, threads, split_gap_size=100, min_contig_len=200, mafft_args="--auto", logger=None, sample="Sample"):
    _, ref_seq = read_first_sequence(ref_file)
    combined_fasta = f">REFERENCE_BASE\n{ref_seq}\n"
    if pvga_file and is_valid(pvga_file):
        _, pvga_seq = read_first_sequence(pvga_file)
        if pvga_seq: combined_fasta += f">PVGA_SEQ\n{pvga_seq}\n"
    if contig_file and is_valid(contig_file):
        for i, (_, cseq) in enumerate(parse_fasta_string(Path(contig_file).read_text()).items()): 
            combined_fasta += f">CONTIG_{i}\n{cseq}\n"

    if logger: logger.info(f"[{sample}] Divine Fusion (SHIVER Core): 执行 MAFFT 全局骨架对齐...")
    mafft_cmd = ['mafft', '--thread', str(threads), '--quiet'] + (shlex.split(mafft_args) if mafft_args else ['--auto']) + ['-']
    result = subprocess.run(mafft_cmd, input=combined_fasta, capture_output=True, text=True, check=True)
  
    aln = parse_fasta_string(result.stdout)
    aln_ref = aln.get('REFERENCE_BASE', '').upper()
    raw_contigs = {k: v.upper() for k, v in aln.items() if k != 'REFERENCE_BASE'}
    aln_len = len(aln_ref)

    split_contigs = []
    gap_pat = re.compile(r'-{%d,}' % split_gap_size)
    for cname, cseq in raw_contigs.items():
        last_end = 0; valid_blocks = []
        for match in gap_pat.finditer(cseq):
            valid_blocks.append((last_end, match.start())); last_end = match.end()
        valid_blocks.append((last_end, aln_len))
        for frag_idx, (s, e) in enumerate(valid_blocks):
            frag_seq = cseq[s:e]
            true_bases = len(frag_seq) - frag_seq.count('-') 
            if true_bases >= min_contig_len:
                split_contigs.append({'name': f"{cname}_frag{frag_idx}", 'seq': ('-' * s) + frag_seq + ('-' * (aln_len - e)), 'true_len': true_bases})

    split_contigs.sort(key=lambda x: x['true_len'], reverse=True)
    fused_seq_chars = []
    for i in range(aln_len):
        chosen_base = '-'
        for contig in split_contigs:
            b = contig['seq'][i]
            if b not in ['-', 'N', '?']: 
                chosen_base = b; break
        if chosen_base == '-':
            b_ref = aln_ref[i]
            chosen_base = b_ref if b_ref not in ['-', 'N', '?'] else 'N'
        fused_seq_chars.append(chosen_base)

    final_seq = "".join(fused_seq_chars).replace('-', '')
    write_fasta(out_file, f"{sample}_Divine_Fusion", final_seq)
    if logger: logger.info(f"[{sample}] Divine Fusion 完美闭环。骨架长度: {len(final_seq)}bp")
    return True

def check_and_trim_circularity(fasta_in, fasta_out, work_dir, avg_read_len=150, min_identity=95.0, logger=None, sample="Sample"):
    if logger: logger.info(f"[{sample}] 启动 Virseqimprover 环化检测引擎...")
    p_dir = Path(work_dir) / "circularity_check"; p_dir.mkdir(exist_ok=True, parents=True)
    seq_name, seq_str = read_first_sequence(fasta_in)
    scaffold_len = len(seq_str)
  
    extract_len = math.ceil(avg_read_len * 0.8) if scaffold_len < 500 else int(avg_read_len)
    if scaffold_len <= extract_len * 2:
        shutil.copy(fasta_in, fasta_out); return False

    subj_file, query_file = p_dir / "subject_head.fasta", p_dir / "query_tail.fasta"
    write_fasta(str(subj_file), seq_name, seq_str[ : scaffold_len - (extract_len * 2)])
    write_fasta(str(query_file), seq_name, seq_str[scaffold_len - (extract_len * 2) : scaffold_len - extract_len])
  
    blast_db = p_dir / "subj_db"
    run_cmd(['makeblastdb', '-in', str(subj_file), '-dbtype', 'nucl', '-out', str(blast_db)], ignore_error=True)
    blast_res = p_dir / "blast_res.txt"
    run_cmd(['blastn', '-query', str(query_file), '-db', str(blast_db), '-outfmt', '6 qseqid sseqid pident length qstart qend sstart send', '-out', str(blast_res), '-num_threads', '1'], ignore_error=True)
  
    is_circular, subject_start, min_aln_len = False, 0, int(extract_len * 0.90)
    if blast_res.exists():
        with open(blast_res, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) < 8: continue
                pident, aln_len, sstart = float(parts[2]), int(parts[3]), int(parts[6])
                if pident >= min_identity and aln_len >= min_aln_len:
                    is_circular, subject_start = True, sstart; break
  
    if is_circular:
        if logger: logger.info(f"[{sample}] 🌟 成功检测到环化特征！切除首尾重叠冗余区...")
        write_fasta(fasta_out, f"{seq_name}_[Circular=True]", seq_str[ : subject_start - 1] + seq_str[scaffold_len - (extract_len * 2) :])
        return True
    else: shutil.copy(fasta_in, fasta_out); return False

# =====================================================================
# V9.0 新增: Coverage & Polishing 可视化模块
# =====================================================================
def run_minimap2_plot(ref_fasta, query_fasta, paf_out):
    subprocess.run(f"minimap2 -x asm10 -c {ref_fasta} {query_fasta} > {paf_out} 2>/dev/null", shell=True)

def parse_paf_with_cigar(paf_file):
    alignments, ref_length = [], 0
    if not os.path.exists(paf_file): return alignments, ref_length
    with open(paf_file, 'r') as f:
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
            else: match_segments = [{'start': r_start, 'len': r_end - r_start}]
            alignments.append({'q_name': q_name, 'q_len': q_len, 'strand': strand, 'start': r_start, 'end': r_end, 'map_span': r_end - r_start, 'match_segments': match_segments})
    alignments.sort(key=lambda x: x['start'])
    return alignments, ref_length

def draw_alignment_block(ax, aln, rect_y, color):
    start, map_span = aln['start'], aln['map_span']
    ax.add_patch(patches.Rectangle((start, rect_y), map_span, 0.4, facecolor='#ff9999', edgecolor='none', alpha=0.9, zorder=2))
    for seg in aln['match_segments']:
        ax.add_patch(patches.Rectangle((seg['start'], rect_y), seg['len'], 0.4, facecolor=color, edgecolor='none', alpha=0.9, zorder=3))
    ax.add_patch(patches.Rectangle((start, rect_y), map_span, 0.4, facecolor='none', edgecolor='black', lw=1, zorder=4))

def find_n_regions(fasta_path, min_n_run=3):
    """从 FASTA 中提取连续 N 区域 [(start, end), ...]"""
    regions = []
    try:
        from Bio import SeqIO
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
            break  # 只处理第一条序列
    except Exception:
        pass
    return regions

def draw_n_overlay(ax, n_regions, y_top, y_bottom, color='#333333', alpha=0.45):
    """在指定 y 范围叠加 N 区域标记"""
    for n_start, n_end in n_regions:
        ax.add_patch(patches.Rectangle(
            (n_start, y_bottom), n_end - n_start, y_top - y_bottom,
            facecolor=color, edgecolor='none', alpha=alpha, zorder=6, hatch='////'))
    return len(n_regions)

def plot_coverage_compact(tool_data, ref_length, out_prefix, n_regions=None):
    colors = plt.cm.Set2.colors
    num_tools = len([t for t in tool_data.values() if t])
    fig, ax = plt.subplots(figsize=(16, max(5, num_tools * 1.8)))

    # 参考灰条 + 长度标注
    ax.add_patch(patches.Rectangle((0, 0), ref_length, 0.4, facecolor='lightgray', edgecolor='black', lw=1, zorder=2))
    ax.text(ref_length + ref_length * 0.005, 0.2, f'{ref_length} bp', va='center', fontsize=10, color='gray', fontstyle='italic')

    # 叠加 N 区域
    if n_regions:
        draw_n_overlay(ax, n_regions, 0.4, 0, color='#222222', alpha=0.55)

    y_pos, yticks_pos, yticks_labels, tool_idx = 1.0, [0.2], ['Reference'], 0
    for tool_name, alignments in tool_data.items():
        if not alignments: continue
        color = colors[tool_idx % len(colors)]; tool_idx += 1
        for aln in alignments:
            draw_alignment_block(ax, aln, y_pos, color)
            if aln['map_span'] > ref_length * 0.08:
                short_name = aln['q_name'][:14] + ".." if len(aln['q_name']) > 16 else aln['q_name']
                ax.text(aln['start'] + aln['map_span']/2, y_pos + 0.2, f"{short_name}\n{aln['q_len']}bp", ha='center', va='center', fontsize=8, color='black', fontweight='bold', zorder=5)
        yticks_pos.append(y_pos + 0.2); yticks_labels.append(tool_name); y_pos += 1.2
    ax.set_xlim(-ref_length * 0.15, ref_length * 1.05); ax.set_ylim(-0.5, y_pos)
    for s in ['top', 'right', 'left']: ax.spines[s].set_visible(False)
    ax.set_yticks(yticks_pos); ax.set_yticklabels(yticks_labels, fontsize=12, fontweight='bold')
    ax.set_xlabel('Reference Genomic Coordinates (bp)', fontsize=14, fontweight='bold')
    ax.set_title('Genome Assembly Evolution (Compact View)', fontsize=18, fontweight='bold', pad=20)
    ax.grid(axis='x', linestyle='--', alpha=0.5, zorder=0)
    legend_elements = [patches.Patch(color='#ff9999', label='Gap (Deletion)')]
    if n_regions:
        legend_elements.append(patches.Patch(facecolor='#222222', alpha=0.55, hatch='////', label=f'N-Regions ({len(n_regions)} sites)'))
    ax.legend(handles=legend_elements, loc='upper right', frameon=True)
    plt.tight_layout()
    plt.savefig(f"{out_prefix}_Compact.png", dpi=300); plt.savefig(f"{out_prefix}_Compact.pdf", dpi=300); plt.close()

def plot_coverage_stacked(tool_data, ref_length, out_prefix, n_regions=None):
    colors = plt.cm.Set2.colors
    draw_info, y_cursor, yticks_pos, yticks_labels = {}, 1.0, [0.2], ['Reference']
    for tool_name, alignments in tool_data.items():
        if not alignments: continue
        levels = [] 
        for aln in alignments:
            placed = False
            for i, end_pos in enumerate(levels):
                if aln['start'] > end_pos + (ref_length * 0.01):
                    aln['level'] = i; levels[i] = aln['end']; placed = True; break
            if not placed: aln['level'] = len(levels); levels.append(aln['end'])
        max_level = len(levels); tool_height = max_level * 0.6
        draw_info[tool_name] = {'alignments': alignments, 'base_y': y_cursor, 'max_level': max_level}
        yticks_pos.append(y_cursor + (tool_height / 2) - 0.1); yticks_labels.append(tool_name); y_cursor += tool_height + 0.6 
  
    fig, ax = plt.subplots(figsize=(16, max(5, y_cursor * 0.8)))
    # 参考灰条 + 长度标注
    ax.add_patch(patches.Rectangle((0, 0), ref_length, 0.4, facecolor='lightgray', edgecolor='black', lw=1, zorder=2))
    ax.text(ref_length + ref_length * 0.005, 0.2, f'{ref_length} bp', va='center', fontsize=10, color='gray', fontstyle='italic')
    # 叠加 N 区域
    if n_regions:
        draw_n_overlay(ax, n_regions, 0.4, 0, color='#222222', alpha=0.55)

    tool_idx = 0
    for tool_name, info in draw_info.items():
        color = colors[tool_idx % len(colors)]; tool_idx += 1; base_y = info['base_y']
        ax.add_patch(patches.Rectangle((-ref_length*0.05, base_y - 0.1), ref_length*1.1, info['max_level']*0.6, facecolor='whitesmoke', edgecolor='none', alpha=0.5, zorder=1))
        for aln in info['alignments']:
            rect_y = base_y + aln['level'] * 0.6; draw_alignment_block(ax, aln, rect_y, color)
            if aln['map_span'] > ref_length * 0.04:
                short_name = aln['q_name'][:15] + ".." if len(aln['q_name']) > 17 else aln['q_name']
                ax.text(aln['start'] + aln['map_span']/2, rect_y + 0.2, f"{short_name}\n{aln['q_len']}bp", ha='center', va='center', fontsize=8, color='black', fontweight='bold', zorder=5)

    ax.set_xlim(-ref_length * 0.15, ref_length * 1.05); ax.set_ylim(-0.5, y_cursor)
    for s in ['top', 'right', 'left']: ax.spines[s].set_visible(False)
    ax.set_yticks(yticks_pos); ax.set_yticklabels(yticks_labels, fontsize=12, fontweight='bold')
    ax.set_xlabel('Reference Genomic Coordinates (bp)', fontsize=14, fontweight='bold')
    ax.set_title('Genome Assembly Evolution (Stacked View)', fontsize=18, fontweight='bold', pad=20)
    ax.grid(axis='x', linestyle='--', alpha=0.5, zorder=0)
    legend_elements = [patches.Patch(color='#ff9999', label='Gap (Deletion)')]
    if n_regions:
        legend_elements.append(patches.Patch(facecolor='#222222', alpha=0.55, hatch='////', label=f'N-Regions ({len(n_regions)} sites)'))
    ax.legend(handles=legend_elements, loc='upper right', frameon=True)
    plt.tight_layout()
    plt.savefig(f"{out_prefix}_Stacked.png", dpi=300); plt.savefig(f"{out_prefix}_Stacked.pdf", dpi=300); plt.close()

# --- GB 解析与测序深度绘图 ---
def parse_genbank(gb_file):
    genes_data = []
    try:
        record = SeqIO.read(gb_file, "genbank")
        seq_length = len(record.seq)
        for feature in record.features:
            if feature.type == "gene":
                start, end = int(feature.location.start), int(feature.location.end)
                gene_name = feature.qualifiers.get("gene", ["Unknown"])[0]
                genes_data.append({"name": gene_name, "start": start, "end": end})
        return genes_data, seq_length
    except Exception as e:
        import logging
        logging.getLogger("virus-full").warning(f"GenBank 文件解析失败 ({gb_file}): {e}, 基因轨道将为空")
        return [], 0

def draw_gene_track(ax_gene, genes_data, max_pos):
    ax_gene.set_facecolor('#f8f9fa')
    gene_colors = ['#8dd3c7', '#ffffb3', '#bebada', '#fb8072', '#80b1d3', '#fdb462', '#b3de69']
    for i, gene in enumerate(genes_data):
        start, end, name = gene['start'], gene['end'], gene['name']
        width = end - start
        color = gene_colors[i % len(gene_colors)]
        rect = patches.Rectangle((start, 0.2), width, 0.6, facecolor=color, edgecolor='black', linewidth=0.8)
        ax_gene.add_patch(rect)
        ax_gene.text(start + width / 2, 0.5, name, ha='center', va='center', fontweight='bold', fontsize=11)

    xticks = np.arange(0, max_pos + 512, 512)
    ax_gene.set_xlim(0, max_pos)
    ax_gene.set_xticks(xticks)
    ax_gene.set_xlabel('Genome Position (bp)', fontweight='bold', fontsize=13)
    ax_gene.set_ylim(0, 1)
    ax_gene.set_yticks([])
    for spine in ['top', 'right', 'left']: ax_gene.spines[spine].set_visible(False)
    ax_gene.spines['bottom'].set_color('gray')
    ax_gene.spines['bottom'].set_linewidth(1.5)

def generate_polishing_plots(tsv_file, gb_file, out_dir, logger, sample):
    if not CAN_PLOT: return
    logger.info(f"[{sample}] Step 9 (Plot): 开始生成测序深度与碱基组成图表...")
    out_dir = Path(out_dir)
    df = pd.read_csv(tsv_file, sep='\t').rename(columns={'-': 'Other'})
    x = df['Pos']
    has_gb, max_pos = False, int(df['Pos'].max())
    if gb_file and Path(gb_file).exists():
        genes_data, gb_seq_len = parse_genbank(gb_file)
        if genes_data: has_gb, max_pos = True, max(max_pos, gb_seq_len)

    # 图1: 碱基组成
    fig1, axes1 = plt.subplots(2, 1, figsize=(16, 9), facecolor='#f8f9fa', gridspec_kw={'height_ratios': [6, 1]}, sharex=True) if has_gb else plt.subplots(1, 1, figsize=(16, 7), facecolor='#f8f9fa')
    ax1_cov = axes1[0] if has_gb else axes1
    fig1.subplots_adjust(hspace=0.05)
    ax1_cov.set_facecolor('#f8f9fa')
    ax1_cov.set_title(f'[{sample}] Viral Genome Base Composition', fontweight='bold', fontsize=16, pad=15)
    colors, labels = ['#66e146', '#ffb944', '#ed4242', '#418aed', '#000000'], ['A', 'C', 'G', 'T', 'Other']
    ax1_cov.stackplot(x, df['A'], df['C'], df['G'], df['T'], df['Other'], labels=labels, colors=colors, linewidth=0)
    ax1_cov.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, p: f'{int(y):,}'))
    ax1_cov.set_ylabel('Count (Linear Scale)', fontweight='bold', fontsize=13)
    ax1_cov.legend(handles=[patches.Patch(color=c, label=l) for c, l in zip(colors, labels)], title='$\\bf{Legend}$', loc='upper left', bbox_to_anchor=(1.01, 1), frameon=False, fontsize=11)
    for spine in ['top', 'right', 'bottom']: ax1_cov.spines[spine].set_visible(False)
    ax1_cov.spines['left'].set_color('gray')
    ax1_cov.tick_params(axis='x', length=0)
    if not has_gb: ax1_cov.set_xlabel('Genome Position (bp)', fontweight='bold', fontsize=13)
    else: draw_gene_track(axes1[1], genes_data, max_pos)
    fig1.savefig(out_dir / f'{sample}_plot_1_base_composition.png', dpi=300, bbox_inches='tight')
    plt.close(fig1)

    # 图2: 测序总深度
    fig2, axes2 = plt.subplots(2, 1, figsize=(16, 9), facecolor='#f8f9fa', gridspec_kw={'height_ratios': [6, 1]}, sharex=True) if has_gb else plt.subplots(1, 1, figsize=(16, 7), facecolor='#f8f9fa')
    ax2_cov = axes2[0] if has_gb else axes2
    fig2.subplots_adjust(hspace=0.05)
    ax2_cov.set_facecolor('#f8f9fa')
    ax2_cov.set_title(f'[{sample}] Viral Genome Sequencing Read Depth', fontweight='bold', fontsize=16, pad=15)
    y_total = df['Total']
    ax2_cov.fill_between(x, y_total, color='#4A90E2', alpha=0.8, linewidth=0)
    mean_depth = y_total.mean()
    ax2_cov.axhline(y=mean_depth, color='#e74c3c', linestyle='--', linewidth=2, label=f'Average Depth: {mean_depth:,.0f}x')
    ax2_cov.legend(loc='upper right', frameon=False, fontsize=12)
    ax2_cov.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, p: f'{int(y):,}'))
    ax2_cov.set_ylabel('Total Read Depth', fontweight='bold', fontsize=13)
    for spine in ['top', 'right', 'bottom']: ax2_cov.spines[spine].set_visible(False)
    ax2_cov.spines['left'].set_color('gray')
    ax2_cov.tick_params(axis='x', length=0)
    if not has_gb: ax2_cov.set_xlabel('Genome Position (bp)', fontweight='bold', fontsize=13)
    else: draw_gene_track(axes2[1], genes_data, max_pos)
    fig2.savefig(out_dir / f'{sample}_plot_2_total_depth.png', dpi=300, bbox_inches='tight')
    plt.close(fig2)

# =====================================================================
# Temp File Cleanup
# =====================================================================
def cleanup_temp_files(work_dir, logger, sample):
    logger.info(f"[{sample}] ===== 开始清理中间临时文件，释放磁盘空间 =====")
    work_p = Path(work_dir)
  
    step1 = work_p / "1.DeNovo_Assembly"
    if step1.exists():
        if (step1 / "penguin").exists() and (step1 / "penguin").is_file(): shutil.move(str(step1/"penguin"), str(step1/"penguin.contigs.fasta"))
        if (step1 / "spades/contigs.fasta").exists(): shutil.copy(str(step1/"spades/contigs.fasta"), str(step1/"spades.contig.fasta"))
        if (step1 / "megahit/final.contigs.fa").exists(): shutil.copy(str(step1/"megahit/final.contigs.fa"), str(step1/"megahit.contigs.fasta"))
        for p in [step1/"megahit", step1/"spades"]:
            if p.exists() and p.is_dir(): shutil.rmtree(p)
        for p in step1.glob("*split_tmp*"):
            if p.is_dir(): shutil.rmtree(p)

    step2 = work_p / "2.RefineC_Merge_Raw"
    if step2.exists():
        for p in step2.glob("all_tools_rc_tmp.*"):
            if p.is_dir(): shutil.rmtree(p)
            else: p.unlink()

    step3 = work_p / "3.Shiver_Cleanup"
    if step3.exists() and (step3 / "blast_db").exists(): shutil.rmtree(step3 / "blast_db")

    step6 = work_p / "6.Pre_Fusion_Merge"
    if step6.exists() and (step6 / "pre_fusion_tmp").exists(): shutil.rmtree(step6 / "pre_fusion_tmp")

    step7 = work_p / "7.Final_rmDup"
    if step7.exists() and (step7 / "tmp_rmdup").exists(): shutil.rmtree(step7 / "tmp_rmdup")

    step10 = work_p / "10.Gap_Filling"
    if step10.exists():
        if (step10 / "temp").exists(): shutil.rmtree(step10 / "temp")
        for p in step10.glob("closed_assembly*"):
            if p.is_file() and not p.name.endswith("gap_filled_final.fasta"): p.unlink()
    logger.info(f"[{sample}] 清理工作完成，工作目录已瘦身！")

# =====================================================================
# RefineC & Shiver Cleanup Tools
# =====================================================================
def run_refinec_split(input_fasta, out_dir, prefix, threads, frag_len, logger, sample):
    if not is_valid(input_fasta): return None
    split_dir = Path(out_dir) / f"{prefix}_split_tmp"; split_dir.mkdir(exist_ok=True, parents=True)
    run_cmd(['refineC', 'split', '--threads', str(threads), '--contigs', str(input_fasta), '--prefix', prefix, '--output', str(split_dir), '--frag-min-len', str(frag_len)], cwd=None, log_file=split_dir / f"{prefix}_split.log", logger=logger, sample_name=sample, ignore_error=True)
    final_out = Path(out_dir) / f"{prefix}.split.fasta"
    found_file = next((f for f in split_dir.rglob("*.split.fasta*") if f.is_file() and f.stat().st_size > 0), None) or next((f for f in Path(out_dir).glob(f"*{prefix}*.split.fasta*") if f.is_file() and f.stat().st_size > 0), None)
    extract_and_move_fasta(found_file if found_file else input_fasta, final_out)
    return final_out

def run_refinec_merge(input_fasta, out_dir, prefix, threads, logger, sample):
    if not is_valid(input_fasta): return None
    merge_dir = Path(out_dir) / f"{prefix}_tmp"; merge_dir.mkdir(exist_ok=True, parents=True)
    logger.info(f"[{sample}] 执行 refineC merge ({prefix})...")
    run_cmd(['refineC', 'merge', '--threads', str(threads), '--contigs', str(input_fasta), '--prefix', prefix, '--output', str(merge_dir)], cwd=None, log_file=merge_dir / f"{prefix}_merge.log", logger=logger, sample_name=sample)
    final_out = Path(out_dir) / f"{prefix}.merged.fasta"
    found_file = next((f for f in merge_dir.rglob("*.merged.fasta*") if f.is_file() and f.stat().st_size > 0), None) or next((f for f in Path(out_dir).glob(f"*{prefix}*.merged.fasta*") if f.is_file() and f.stat().st_size > 0), None)
    if found_file: extract_and_move_fasta(found_file, final_out); return final_out
    return None

def run_shiver_cleanup(raw_contigs, ref_fasta, out_fasta, work_dir, args, logger, sample):
    logger.info(f"[{sample}] 执行 Shiver-like Contig 剪切与过滤 (已解禁长度限制)...")
    work_p = Path(work_dir); blast_db_dir = work_p / "blast_db"; blast_db_dir.mkdir(exist_ok=True)
    local_ref = blast_db_dir / "target_ref.fasta"; shutil.copy(ref_fasta, local_ref)
    ref_name, _ = read_first_sequence(local_ref)
    run_cmd(['makeblastdb', '-in', str(local_ref), '-dbtype', 'nucl', '-out', str(blast_db_dir/"refdb")], logger=logger, sample_name=sample, ignore_error=True)
  
    raw_seqs = parse_fasta_string(Path(raw_contigs).read_text())
    filtered_fasta = work_p / "filtered_raw.fasta"
    with open(filtered_fasta, 'w') as f:
        for name, seq in raw_seqs.items():
            if len(seq) > 0: f.write(f">{name}\n{seq.upper()}\n")
    if not is_valid(filtered_fasta): return False

    blast_out = work_p / "blast_hits.tsv"
    run_cmd(['blastn', '-query', str(filtered_fasta), '-db', str(blast_db_dir/"refdb"), '-outfmt', '6 qseqid sseqid evalue pident qlen qstart qend sstart send', '-out', str(blast_out), '-evalue', '1e-5'], logger=logger, sample_name=sample, ignore_error=True)
    if not blast_out.exists() or blast_out.stat().st_size == 0: return False
      
    hit_dict = {}
    with open(blast_out, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) < 9: continue
            qseqid, qstart, qend, sstart, send = parts[0], int(parts[5]), int(parts[6]), int(parts[7]), int(parts[8])
            if qstart > qend: qstart, qend = qend, qstart
            hit_dict.setdefault(qseqid, []).append((qstart, qend, sstart, send))
          
    corrected_fasta = work_p / "corrected_contigs.fasta"; corrected_count = 0
    with open(corrected_fasta, 'w') as out_f:
        for qseqid, hits in hit_dict.items():
            if qseqid not in raw_seqs: continue
            original_seq = raw_seqs[qseqid]; hits.sort(key=lambda x: x[0])
            for i, (qs, qe, ss, se) in enumerate(hits):
                frag_seq = reverse_complement(original_seq[qs-1:qe]) if ss > se else original_seq[qs-1:qe]
                if len(frag_seq) > 0: out_f.write(f">{qseqid}_frag{i+1}\n{frag_seq.upper()}\n"); corrected_count += 1
    if corrected_count == 0: return False

    aln_fasta = work_p / "aligned_for_cut.fasta"
    run_cmd(['mafft', '--quiet', '--add', str(corrected_fasta), str(local_ref)], log_file=aln_fasta, logger=logger, sample_name=sample, ignore_error=True)
  
    aln_seqs = parse_fasta_string(Path(aln_fasta).read_text())
    gap_regex = re.compile(r'-{%d,}' % args.split_gap_size)
    with open(out_fasta, 'w') as final_f:
        for name, aln_seq in aln_seqs.items():
            if name == ref_name: continue
            frag_idx = 1
            for frag in gap_regex.split(aln_seq.strip('-').strip('?')):
                clean_frag = frag.replace('-', '').replace('?', '').upper()
                if len(clean_frag) > 0: final_f.write(f">{name}_split{frag_idx}\n{clean_frag}\n"); frag_idx += 1
    return True

def calculate_sequence_stats(fasta_file):
    if not is_valid(fasta_file): return 0, 0, 0
    total_len, non_n_count, contig_lengths = 0, 0, []
    with smart_open(fasta_file, 'rt') as f:
        current_seq = []
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_seq:
                    seq_str = "".join(current_seq).upper().replace('?', 'N').replace('-', 'N')
                    total_len += len(seq_str); non_n_count += len(seq_str) - seq_str.count('N')
                    contig_lengths.extend([len(c) for c in seq_str.split('N') if len(c) > 0])
                current_seq = []
            else: current_seq.append(line)
        if current_seq:
            seq_str = "".join(current_seq).upper().replace('?', 'N').replace('-', 'N')
            total_len += len(seq_str); non_n_count += len(seq_str) - seq_str.count('N')
            contig_lengths.extend([len(c) for c in seq_str.split('N') if len(c) > 0])
          
    if not contig_lengths: return total_len, non_n_count, 0
    contig_lengths.sort(reverse=True)
    target = sum(contig_lengths) / 2.0; cumsum = 0
    for l in contig_lengths:
        cumsum += l
        if cumsum >= target: return total_len, non_n_count, l
    return total_len, non_n_count, 0

def check_assembly_quality(fasta_path, ref_len, chk_len, chk_n50):
    if not is_valid(fasta_path): return False, 0, 0, 0, 0
    tot_len, non_n, n50 = calculate_sequence_stats(fasta_path)
    return (tot_len >= chk_len * ref_len) and (n50 >= chk_n50 * ref_len), tot_len, non_n, (tot_len - non_n), n50

def parse_reads_input(inputs, logger):
    if not inputs: return {}
    if len(inputs) == 1 and os.path.isdir(inputs[0]):
        raw_samples = {}
        for f in Path(inputs[0]).glob("*"):
            if not f.is_file() or not any(f.name.lower().endswith(ext) for ext in VALID_EXTENSIONS): continue
            name, abs_f = f.name, str(f.absolute())
            r1_m = re.search(r'_(R?1|1)(?:_clean)?\.(fast[aq]|f[aq])(?:\.gz)?$', name, re.IGNORECASE)
            r2_m = re.search(r'_(R?2|2)(?:_clean)?\.(fast[aq]|f[aq])(?:\.gz)?$', name, re.IGNORECASE)
            if r1_m: raw_samples.setdefault(name[:r1_m.start()], {})['r1'] = abs_f
            elif r2_m: raw_samples.setdefault(name[:r2_m.start()], {})['r2'] = abs_f
            else:
                p = name
                for ext in VALID_EXTENSIONS:
                    if p.lower().endswith(ext): p = p[:-len(ext)]; break
                raw_samples.setdefault(p, {})['r1'] = abs_f; raw_samples[p]['single'] = True
        res = {k: v for k, v in raw_samples.items() if 'r1' in v}
        for k, v in res.items(): v['single'] = v.get('single', 'r2' not in v)
        return res
      
    files = [os.path.abspath(f) for f in inputs if os.path.isfile(f)]
    if not files: return {}
    sample_name = os.path.basename(files[0])
    for ext in VALID_EXTENSIONS:
        if sample_name.lower().endswith(ext): sample_name = sample_name[:-len(ext)]; break
    sample_name = re.sub(r'_(R?1|1)(?:_clean)?$', '', sample_name, flags=re.IGNORECASE)
    res = {sample_name: {'r1': files[0], 'single': True}}
    if len(files) >= 2: res[sample_name]['r2'] = files[1]; res[sample_name]['single'] = False
    return res

# =====================================================================
# V9.0 核心流水线引擎 (12 步终极版 - 支持快车道模式)
# =====================================================================
def run_pipeline(sample, sample_data, orig_ref, ref_len, out_root, active_tools, args, logger):
    sample_dir = Path(out_root) / sample
    sample_dir.mkdir(parents=True, exist_ok=True)
  
    # --- 架构目录定义 ---
    d1_asm      = sample_dir / "1.DeNovo_Assembly"
    d2_rc_merge = sample_dir / "2.RefineC_Merge_Raw"
    d3_cleanup  = sample_dir / "3.Shiver_Cleanup"
    d4_ref_m1   = sample_dir / "4.Ref_Merged_1"
    d5_pvga     = sample_dir / "5.PVGA_Extension"
    d6_pre_m    = sample_dir / "6.Pre_Fusion_Merge"
    d7_rmdup    = sample_dir / "7.Final_rmDup"
    d8_ref_m2   = sample_dir / "8.Ref_Merged_2"
    d9_cons     = sample_dir / "9.Consensus_Polish"
    d10_gapfill = sample_dir / "10.Gap_Filling"
    d11_circ    = sample_dir / "11.Circularity_Check"
    d12_plot    = sample_dir / "12.Coverage_Visualization"
  
    # --- 文件路径定义 ---
    f2_asm      = d2_rc_merge / "2.refinec_merged_raw.fasta"
    f3_clean    = d3_cleanup  / "3.cleaned_viral_contigs.fasta"
    f4_ref_m1   = d4_ref_m1   / "4.ref_merged_1.fasta"
    f5_pvga_raw = d5_pvga     / "5.pvga_raw.fasta"
    f5_pvga_cut = d5_pvga     / "5.pvga_evaluated_cut.fasta" 
    f6_pre_m    = d6_pre_m    / "6.pre_fusion_merged.fasta"
    f7_rmdup    = d7_rmdup    / "7.rmDup_cleaned.fasta"
    f8_ref_m2   = d8_ref_m2   / "8.ref_merged_2.fasta"
    f9_cons     = d9_cons     / "9.final_consensus.fasta"
    f10_gapf    = d10_gapfill / "10.gap_filled_final.fasta"
    f11_circ    = sample_dir  / "11.Ultimate_Circular_Result.fasta"
    f13_stats   = sample_dir  / "13.Evolution_Stats.tsv"
    f_ok        = sample_dir  / "Fully-assembled.ok"

    asm_r1, asm_r2, asm_is_single = sample_data['assembly']['r1'], sample_data['assembly'].get('r2'), sample_data['assembly']['single']
    pvga_r1, pvga_r2, pvga_is_single = sample_data['pvga']['r1'], sample_data['pvga'].get('r2'), sample_data['pvga']['single']
    cons_r1, cons_r2, cons_is_single = sample_data['consensus']['r1'], sample_data['consensus'].get('r2'), sample_data['consensus']['single']
    cons_reads_args = [cons_r1] if cons_is_single else [cons_r1, cons_r2]
  
    stats_history = []
    def log_stat(step_name, file_path):
        is_perf, t_len, non_n, n_c, n50 = check_assembly_quality(str(file_path), ref_len, args.chk_len, args.chk_n50)
        stats_history.append({'Step': step_name, 'Total_Length': t_len, 'Non_N_Bases': non_n, 'N_Count': n_c, 'Contig_N50': n50})
        return is_perf

    def finish_pipeline(msg, is_perfect_early=False, final_eval_file=None):
        pd.DataFrame(stats_history).to_csv(f13_stats, sep='\t', index=False)
      
        if CAN_PLOT:
            d12_plot.mkdir(exist_ok=True)
            plot_inputs_ordered = [
                ("2.RefineC_Raw", f2_asm),
                ("3.Shiver_Cleaned", f3_clean),
                ("4.Ref_Merged_1", f4_ref_m1),
                ("5.PVGA_Split", f5_pvga_cut),
                ("6.Pre_Fusion_Merge", f6_pre_m),
                ("7.rmDup_Cleaned", f7_rmdup),
                ("8.Ref_Merged_2", f8_ref_m2),
                ("9.Consensus_Polished", f9_cons),
                ("10.Gap_Filled", f10_gapf),
                ("11.Ultimate_Result", f11_circ)
            ]
            # 自动过滤掉因为“快车道”而未生成的中间步骤文件，不会报错
            valid_plots = [item for item in plot_inputs_ordered if is_valid(item[1])]
          
            if valid_plots:
                logger.info(f"[{sample}] Step 12: 绘制全生命周期演化图谱...")
                tool_data, global_ref_len = {}, 0
                with tempfile.TemporaryDirectory(dir=os.environ.get("TMPDIR")) as tmpdir:
                    for t_name, f_path in valid_plots:
                        paf_f = os.path.join(tmpdir, f"{t_name}.paf")
                        run_minimap2_plot(orig_ref, f_path, paf_f)
                        aln, r_len = parse_paf_with_cigar(paf_f)
                        global_ref_len = max(global_ref_len, r_len)
                        tool_data[t_name] = aln
              
                if global_ref_len > 0:
                    # 从最终结果提取 N 区域用于高亮
                    n_regs = []
                    for final_fa in [f11_circ, f10_gapf, f9_cons]:
                        if final_fa.exists():
                            n_regs = find_n_regions(str(final_fa))
                            if n_regs: break
                    plot_prefix = str(d12_plot / f"{sample}_Coverage")
                    plot_coverage_compact(tool_data, global_ref_len, plot_prefix, n_regions=n_regs)
                    plot_coverage_stacked(tool_data, global_ref_len, plot_prefix, n_regions=n_regs)
      
        if not getattr(args, 'keep_tmp', False):
            cleanup_temp_files(sample_dir, logger, sample)

        is_ok = is_perfect_early
        if final_eval_file and not is_ok:
            is_ok, _, _, _, _ = check_assembly_quality(str(final_eval_file), ref_len, args.chk_len, args.chk_n50)
      
        if is_ok:
            open(f_ok, 'a').close()
            return True, f"{msg} (达到完美阈值 Fully-assembled.ok)", sample, stats_history
        else:
            return True, f"{msg} (流程结束，但未能达到设定完美阈值)", sample, stats_history

    # ==========================
    # 快车道状态控制变量
    # ==========================
    fast_track_source = None
    is_perfect_early = False

    # --- Step 1 & 2: Assembly ---
    d1_asm.mkdir(exist_ok=True); d2_rc_merge.mkdir(exist_ok=True)
    if not is_valid(f2_asm):
        logger.info(f"[{sample}] Step 1&2: De novo 组装与初步合并...")
        contig_files = []
        if 'megahit' in active_tools:
            mh_out = d1_asm / "megahit"
            cmd_mh = ['megahit', '-t', str(args.threads), '--out-dir', str(mh_out), '--out-prefix', 'megahit', '--min-contig-len', str(args.min_contig_len)] + (['-r', asm_r1] if asm_is_single else ['-1', asm_r1, '-2', asm_r2])
            if run_cmd(cmd_mh, cwd=d1_asm, log_file=d1_asm/"mh.log", logger=logger, sample_name=sample) and (mh_out/"megahit.contigs.fa").exists():
                split_fa = run_refinec_split(mh_out/"megahit.contigs.fa", d1_asm, "megahit", args.threads, args.refineC_frag, logger, sample)
                if split_fa: contig_files.append(split_fa)
        if 'spades' in active_tools:
            sp_out = d1_asm / "spades"
            cmd_sp = ['rnaviralspades.py', '-t', str(args.threads), '-o', str(sp_out)] + (['-s', asm_r1] if asm_is_single else ['-1', asm_r1, '-2', asm_r2])
            if run_cmd(cmd_sp, cwd=d1_asm, log_file=d1_asm/"sp.log", logger=logger, sample_name=sample) and (sp_out/"contigs.fasta").exists():
                split_fa = run_refinec_split(sp_out/"contigs.fasta", d1_asm, "spades", args.threads, args.refineC_frag, logger, sample)
                if split_fa: contig_files.append(split_fa)
        if 'penguin' in active_tools:
            pg_out, pg_tmp = d1_asm / "penguin", d1_asm / "pg_tmp"
            cmd_pg = ['penguin', 'guided_nuclassemble'] + ([asm_r1] if asm_is_single else [asm_r1, asm_r2]) + ['--threads', str(args.threads), str(pg_out), str(pg_tmp), '--min-contig-len', str(args.min_contig_len)]
            if run_cmd(cmd_pg, cwd=d1_asm, log_file=d1_asm/"pg.log", logger=logger, sample_name=sample) and pg_out.exists():
                pg_fasta_list = list(pg_out.glob("*.fa*"))
                if pg_fasta_list: 
                    split_fa = run_refinec_split(max(pg_fasta_list, key=lambda x: x.stat().st_size), d1_asm, "penguin", args.threads, args.refineC_frag, logger, sample)
                    if split_fa: contig_files.append(split_fa)

        if contig_files:
            cat_fa = d2_rc_merge / "all_tools_split.fasta"
            with open(cat_fa, 'w') as outf:
                for c in contig_files: outf.write(Path(c).read_text() + "\n")
            merged_fa = run_refinec_merge(cat_fa, d2_rc_merge, "all_tools_rc", args.threads, logger, sample)
            extract_and_move_fasta(merged_fa if merged_fa else max(contig_files, key=lambda x: x.stat().st_size), f2_asm)
        else: open(f2_asm, 'a').close()

    # --- Step 3: Shiver Cleanup ---
    d3_cleanup.mkdir(exist_ok=True)
    if not is_valid(f3_clean):
        if is_valid(f2_asm) and not run_shiver_cleanup(f2_asm, orig_ref, f3_clean, d3_cleanup, args, logger, sample): open(f3_clean, 'a').close()
        elif not is_valid(f2_asm): open(f3_clean, 'a').close()
    
    # 检查是否触发抛光快车道
    if log_stat('1.DeNovo_Cleaned', f3_clean):
        fast_track_source = f3_clean
        is_perfect_early = True
        logger.info(f"[{sample}] 🌟 触发快车道：初步净化骨架质量已达完美标准，跳过后续延伸融合，直接进入 Step 9 迭代抛光！")

    # 如果未触发快车道，执行标准构建流程
    if not fast_track_source:
        # --- Step 4: Ref Merge 1 ---
        d4_ref_m1.mkdir(exist_ok=True)
        if not is_valid(f4_ref_m1):
            if is_valid(f3_clean): divine_fusion_shiver_style(orig_ref, str(f3_clean), None, str(f4_ref_m1), args.threads, args.split_gap_size, args.min_contig_len, args.mafft_args, logger, sample)
            else: shutil.copy(orig_ref, f4_ref_m1)

        # --- Step 5: PVGA 延伸 & 评估打断 ---
        d5_pvga.mkdir(exist_ok=True)
        if not is_valid(f5_pvga_cut):
            logger.info(f"[{sample}] Step 5: 运行 PVGA 延伸...")
            pvga_input = d5_pvga / f"{sample}_pvga_reads.fasta"
            if not pvga_is_single and pvga_r2 and is_fastq(pvga_r1):
                run_cmd(['bbmerge.sh', f'in1={pvga_r1}', f'in2={pvga_r2}', f'out={d5_pvga}/m.fastq', f'outu1={d5_pvga}/u1.fastq', f'outu2={d5_pvga}/u2.fastq', 'overwrite=t'], cwd=d5_pvga, log_file=d5_pvga/"bbmerge.log", logger=logger, sample_name=sample, ignore_error=True)
                safe_concat_to_fasta([d5_pvga/'m.fastq', d5_pvga/'u1.fastq', d5_pvga/'u2.fastq'] if (d5_pvga/'m.fastq').exists() else [pvga_r1, pvga_r2], d5_pvga/'all.fastq')
                run_cmd(['reformat.sh', f'in={d5_pvga}/all.fastq', f'out={pvga_input}', 'overwrite=t'], cwd=d5_pvga, logger=logger, sample_name=sample)
            else: safe_concat_to_fasta([pvga_r1] + ([pvga_r2] if not pvga_is_single and pvga_r2 else []), pvga_input)

            pvga_sandbox = d5_pvga / "pvga_sandbox"; pvga_sandbox.mkdir(exist_ok=True)
            pvga_out = pvga_sandbox / "pvga_res"
            cmd_pvga = ['pvga', '-r', str(pvga_input.absolute()), '-b', str(f4_ref_m1.absolute()), '-o', str(pvga_out.absolute())] + (shlex.split(args.pvga_args) if args.pvga_args else [])
            run_cmd(cmd_pvga, cwd=pvga_sandbox, log_file=d5_pvga/"pvga.log", logger=logger, sample_name=sample, ignore_error=True)
          
            collected_pvga = [f for f in pvga_out.glob("*.fa*") if "scaffold" not in f.name.lower() and f.stat().st_size > 0] if pvga_out.exists() else []
            if collected_pvga:
                with open(f5_pvga_raw, 'w') as outf:
                    for idx, pvga_f in enumerate(collected_pvga):
                        for head, seq in parse_fasta_string(Path(pvga_f).read_text()).items(): outf.write(f">PVGA_SRC{idx}_{head}\n{seq}\n")
            else: open(f5_pvga_raw, 'a').close()
            try: shutil.rmtree(pvga_sandbox)
            except Exception: pass

            evaluate_and_split_pvga(orig_ref, f5_pvga_raw, f5_pvga_cut, args.threads, gap_size=args.split_gap_size, logger=logger, sample=sample)
        
        # 再次检查是否触发快车道
        if log_stat('2.PVGA_Extension', f5_pvga_cut):
            fast_track_source = f5_pvga_cut
            is_perfect_early = True
            logger.info(f"[{sample}] 🌟 触发快车道：PVGA延伸骨架已达完美标准，跳过后续繁琐融合，直接进入 Step 9 迭代抛光！")

    if not fast_track_source:
        # --- Step 6: 纯净聚合 (Pre-Fusion Merge) ---
        d6_pre_m.mkdir(exist_ok=True)
        if not is_valid(f6_pre_m):
            logger.info(f"[{sample}] Step 6: 原始纯净 Contig 与打断后的 PVGA 序列融合...")
            all_sources = d6_pre_m / "sources.fasta"
            with open(all_sources, 'w') as outf:
                for f in [f3_clean, f5_pvga_cut]: 
                    if is_valid(f): outf.write(Path(f).read_text() + "\n")
            merged_fa = run_refinec_merge(all_sources, d6_pre_m, "pre_fusion", args.threads, logger, sample)
            extract_and_move_fasta(merged_fa if merged_fa else f3_clean, f6_pre_m)
        log_stat('3.Pre_Fusion_Merge', f6_pre_m)

        # --- Step 7: Final rmDup (前置去冗余) ---
        d7_rmdup.mkdir(exist_ok=True)
        if not is_valid(f7_rmdup):
            logger.info(f"[{sample}] Step 7: 执行严格去冗余，净化融合骨架物料...")
            tmp_pl = d7_rmdup / "tmp_rmdup"; tmp_pl.mkdir(exist_ok=True)
            cmd_rmdup = (f"perl {args.rmdup_script} --length {args.rmdup_len} --coverage {args.rmdup_cov} "
                         f"--identity {args.rmdup_iden} --evalue {args.rmdup_evalue} --CPU {args.threads} --tmp {tmp_pl} {f6_pre_m} > {f7_rmdup}")
            run_cmd(cmd_rmdup, cwd=d7_rmdup, log_file=d7_rmdup/"rmdup.log", logger=logger, sample_name=sample)
            if not is_valid(f7_rmdup): shutil.copy(f6_pre_m, f7_rmdup)
        log_stat('4.rmDup_Purification', f7_rmdup)

        # --- Step 8: Ref Merged 2 (终极 Divine Fusion) ---
        d8_ref_m2.mkdir(exist_ok=True)
        if not is_valid(f8_ref_m2):
            logger.info(f"[{sample}] Step 8: 执行 Divine Fusion 构建最终无缝实心骨架...")
            divine_fusion_shiver_style(orig_ref, str(f7_rmdup), None, str(f8_ref_m2), args.threads, args.split_gap_size, args.min_contig_len, args.mafft_args, logger, sample)
        log_stat('5.Fusion_Skeleton_Solid', f8_ref_m2)

    # --- 快车道去冗余: 即使骨架完美，也要执行 rmDup 净化 ---
    if fast_track_source and is_valid(fast_track_source):
        d7_rmdup.mkdir(exist_ok=True)
        if not is_valid(f7_rmdup):
            logger.info(f"[{sample}] 快车道去冗余: 对完美骨架执行 rmDup 净化...")
            tmp_pl = d7_rmdup / "tmp_rmdup"; tmp_pl.mkdir(exist_ok=True)
            cmd_rmdup = (f"perl {args.rmdup_script} --length {args.rmdup_len} --coverage {args.rmdup_cov} "
                         f"--identity {args.rmdup_iden} --evalue {args.rmdup_evalue} --CPU {args.threads} --tmp {tmp_pl} {fast_track_source} > {f7_rmdup}")
            run_cmd(cmd_rmdup, cwd=d7_rmdup, log_file=d7_rmdup/"rmdup.log", logger=logger, sample_name=sample)
            if not is_valid(f7_rmdup):
                shutil.copy(fast_track_source, f7_rmdup)
            else:
                fast_track_source = f7_rmdup  # 替换为去冗余后的骨架
        log_stat('4.rmDup_Purification', f7_rmdup if is_valid(f7_rmdup) else fast_track_source)

    # --- Step 9: Iterative Consensus (后置迭代抛光) 【V9.0 引入快车道输入支持】 ---
    d9_cons.mkdir(exist_ok=True); iter_process_dir = d9_cons / "Iteration_Process"; iter_process_dir.mkdir(exist_ok=True)

    # 根据是否快车道决定抛光的输入文件
    pre_polish_fasta = fast_track_source if fast_track_source else f8_ref_m2
    
    if not is_valid(f9_cons):
        logger.info(f"[{sample}] Step 9: 对目标骨架进行 Reads 级别迭代抛光纠错 (消除早期的测序偏差)...")
        active_ref = str(pre_polish_fasta)
        reads_str = " ".join(shlex.quote(str(p)) for p in cons_reads_args)
        last_tsv = None
      
        for i in range(1, args.iter + 1):
            raw_cons, filled_cons = iter_process_dir / f"iter{i}_raw.fasta", iter_process_dir / f"iter{i}_filled.fasta"
            position_tsv = iter_process_dir / f"iter{i}_position.tsv"
            reads_bam = iter_process_dir / f"iter{i}_mapped.bam"
          
            cmd_compound = (
                f"minimap2 -t {args.threads} -a -x {args.mm2_preset} {active_ref} {reads_str} | "
                f"tee >(viral_consensus -i - -r {active_ref} -o {raw_cons} --min_qual {args.vc_min_qual} --min_depth {args.vc_min_depth} -op {position_tsv}) | "
                f"samtools view -b -@ {args.threads} > {reads_bam}"
            )
            run_cmd(cmd_compound, cwd=iter_process_dir, logger=logger, sample_name=sample, use_bash=True)
          
            if reads_bam.exists() and reads_bam.stat().st_size > 0:
                sorted_bam = iter_process_dir / f"iter{i}_mapped.sorted.bam"
                run_cmd(f"samtools sort -@ {args.threads} {reads_bam} -o {sorted_bam} && samtools index {sorted_bam}", cwd=iter_process_dir, logger=logger, sample_name=sample, use_bash=True)
                reads_bam.unlink()

            if is_valid(raw_cons):
                shutil.copy(raw_cons, f9_cons)
                last_tsv = position_tsv
                if i == args.iter: break
              
                _, cons_seq = read_first_sequence(raw_cons)
                _, ref_seq = read_first_sequence(active_ref)
                try:
                    res = subprocess.run(['mafft', '--quiet'] + (shlex.split(args.mafft_args) if args.mafft_args else ['--auto']) + ['-'], input=f">C\n{cons_seq}\n>R\n{ref_seq}\n", capture_output=True, text=True, check=True)
                    aln = parse_fasta_string(res.stdout)
                    new_cons = [CB if CB in ['A','C','G','T'] else RB for CB, RB in zip(aln['C'].upper(), aln['R'].upper())]
                    write_fasta(str(filled_cons), f"{sample}_Iter{i}", "".join([b for b in new_cons if b != '-']))
                    active_ref = str(filled_cons)
                except Exception: break
            else: break
          
        if not is_valid(f9_cons): shutil.copy(pre_polish_fasta, f9_cons)
      
        # 迭代完成后，利用最终生成的 position_tsv 进行测序深度和基因轨道出图
        if last_tsv and last_tsv.exists() and CAN_PLOT:
            generate_polishing_plots(last_tsv, args.gb, d9_cons, logger, sample)
          
    log_stat('6.Iterative_Consensus', f9_cons)

    # --- Step 10: Gap Filling (gmcloser + abyss-sealer) ---
    d10_gapfill.mkdir(exist_ok=True)
    if fast_track_source:
        # 快车道模式下，骨架的 N50 已经达标，通常无需再进行繁琐的 Gap 填补
        if not is_valid(f10_gapf):
            logger.info(f"[{sample}] 快车道状态：跳过 Step 10 (Gap Filling)，直接进入环化检测...")
            shutil.copy(f9_cons, f10_gapf)
    else:
        if not is_valid(f10_gapf):
            logger.info(f"[{sample}] Step 10: 启动双引擎 Gap 深度填补 (gmcloser + abyss-sealer)...")
          
            # 1. 准备 gmcloser 物料
            gm_queries = d10_gapfill / "gmcloser_queries.fasta"
            with open(gm_queries, 'w') as outf:
                for f in [f3_clean, f5_pvga_cut, f7_rmdup]:
                    if is_valid(f): outf.write(Path(f).read_text() + "\n")
          
            # 2. 运行 gmcloser (结构性回填)
            gm_prefix = "closed_assembly"
            gm_cmd = [
                'gmcloser', '-t', str(f9_cons), '-q', str(gm_queries), '-p', gm_prefix,
                '-l', str(args.gm_l), '-i', str(args.gm_i), '-d', '50', '-mm', '50', 
                '-ms', '50', '-mi', '95', '-et', '-n', str(args.threads)
            ]
            run_cmd(gm_cmd, cwd=d10_gapfill, log_file=d10_gapfill/"gmcloser.log", logger=logger, sample_name=sample, ignore_error=True)
          
            gm_res = d10_gapfill / f"{gm_prefix}.closed.fa"
            active_gap_scaffold = gm_res if is_valid(gm_res) else f9_cons
            if not is_valid(gm_res):
                logger.warning(f"[{sample}] gmcloser 填充未产生有效结果，直接使用 Consensus 进行下级 abyss-sealer 填补...")

            # 3. 运行 abyss-sealer (Reads 布隆填补)
            sealer_prefix = "sealer_closed_viral"
            sealer_cmd = ['abyss-sealer', '-S', str(active_gap_scaffold), '-b', args.sealer_b]
            for k in args.sealer_k.split(','): sealer_cmd.extend(['-k', k.strip()])
            sealer_cmd.extend(['-F', str(args.sealer_F), '-L', str(args.sealer_L), '-P', str(args.threads), '-o', sealer_prefix])
            sealer_cmd.extend(cons_reads_args)
          
            run_cmd(sealer_cmd, cwd=d10_gapfill, log_file=d10_gapfill/"sealer.log", logger=logger, sample_name=sample, ignore_error=True)
          
            sealer_res = d10_gapfill / f"{sealer_prefix}_scaffold.fa"
            if is_valid(sealer_res):
                logger.info(f"[{sample}] Gap 填补成功！")
                shutil.copy(sealer_res, f10_gapf)
            else:
                logger.warning(f"[{sample}] abyss-sealer 填补未改变序列或失败，继承上游序列...")
                shutil.copy(active_gap_scaffold, f10_gapf)
              
    log_stat('7.Gap_Filled', f10_gapf)

    # --- Step 11: Circularity Check & Trim ---
    if not is_valid(f11_circ):
        check_and_trim_circularity(f10_gapf, f11_circ, d11_circ, avg_read_len=args.avg_read_len, min_identity=args.circ_identity, logger=logger, sample=sample)
    log_stat('8.Ultimate_Result', f11_circ)
  
    final_msg = "快车道极速抛光组装完成" if is_perfect_early else "全流程双引擎 Gap 填补深度组装执行完毕"
    return finish_pipeline(final_msg, is_perfect_early=is_perfect_early, final_eval_file=f11_circ)

# =====================================================================
# Main 启动器
# =====================================================================
def main():
    parser = argparse.ArgumentParser(description="OmniVirusAssembler (Production Edition V9.0)")
    g_io = parser.add_argument_group('数据输入与运行模式')
    g_io.add_argument('-r', '--reference', required=True, help="原始参考基因组文件")
    g_io.add_argument('-o', '--output', default="Ultimate_Result", help="输出路径")
    g_io.add_argument('-t', '--threads', type=int, default=8, help="单样本线程数")
    g_io.add_argument('-j', '--jobs', type=int, default=2, help="并行样本数")
    g_io.add_argument('--dry-run', action='store_true', help="预演模式")
    g_io.add_argument('--keep-tmp', action='store_true', help="保留所有中间临时文件 (默认开启自动瘦身)")
  
    g_reads = parser.add_argument_group('解耦 Reads 挂载')
    g_reads.add_argument('--assembly_reads', nargs='+', required=True, help="De novo 组装数据")
    g_reads.add_argument('--pvga_reads', nargs='+', help="PVGA 延伸数据 (缺省继承)")
    g_reads.add_argument('--consensus_reads', nargs='+', help="Consensus 抛光数据 (缺省继承)")

    g_chk = parser.add_argument_group('智能熔断机制')
    g_chk.add_argument('--chk-len', type=float, default=0.98, help="全长>=参考比例")
    g_chk.add_argument('--chk-n50', type=float, default=0.95, help="N50>=参考比例")

    g_denovo = parser.add_argument_group('De novo 参数与净化')
    g_denovo.add_argument('--assembly_tools', default="all", help="启用工具 (megahit,spades,penguin,all)")
    g_denovo.add_argument('--refineC-frag', type=int, default=300, help="refineC 最小拆分长度")
    g_denovo.add_argument('--min-contig-len', type=int, default=200, help="最短过滤长度")
    g_denovo.add_argument('--split-gap-size', type=int, default=100, help="多大 Gap 引发打断 (包含PVGA打断)")

    g_align = parser.add_argument_group('融合与延伸')
    g_align.add_argument('--pvga-args', default="", help="PVGA 附加参数 (如 -n 10)")
    g_align.add_argument('--mafft-args', default="--auto", help="MAFFT 附加参数")

    g_cons = parser.add_argument_group('抛光迭代与出图 (V9.0 升级)')
    g_cons.add_argument('-n', '--iter', type=int, default=3, help="抛光循环次数")
    g_cons.add_argument('--mm2-preset', default="sr", help="Minimap2 预设")
    g_cons.add_argument('--vc-min-qual', type=int, default=10, help="VC 最小质量")
    g_cons.add_argument('--vc-min-depth', type=int, default=1, help="VC 最小深度")
    g_cons.add_argument('--gb', default=None, help="参考基因组 GenBank (.gb) 文件，用于绘制带基因结构的测序深度图")

    g_gapfill = parser.add_argument_group('Gap 深度填补 (V8.7 新增)')
    g_gapfill.add_argument('--gm-l', type=int, default=150, help="gmcloser reads 长度 (-l)")
    g_gapfill.add_argument('--gm-i', type=int, default=350, help="gmcloser 插入片段长度 (-i)")
    g_gapfill.add_argument('--sealer-b', default="500M", help="abyss-sealer 布隆过滤器大小 (-b)")
    g_gapfill.add_argument('--sealer-k', default="120,96,80,64,48", help="abyss-sealer K-mer 梯度列表")
    g_gapfill.add_argument('--sealer-F', type=int, default=150, help="abyss-sealer Flank 长度 (-F)")
    g_gapfill.add_argument('--sealer-L', type=int, default=1000, help="abyss-sealer 探索上限 (-L)")

    g_rmdup = parser.add_argument_group('清理与环化验证')
    g_rmdup.add_argument('--rmdup-script', default=str(SCRIPT_DIR / "utils/genome_rmDuplicates.pl"), help="去重脚本")
    g_rmdup.add_argument('--rmdup-len', type=int, default=100000)
    g_rmdup.add_argument('--rmdup-cov', type=float, default=0.90)
    g_rmdup.add_argument('--rmdup-iden', type=float, default=0.95)
    g_rmdup.add_argument('--rmdup-evalue', default="1e-10")
    g_rmdup.add_argument('--avg-read-len', type=int, default=150, help="测序读长 (用于环化头尾提取)")
    g_rmdup.add_argument('--circ-identity', type=float, default=95.0, help="判定环化的最低相似度阈值")

    args = parser.parse_args()
    logger = setup_logger(args.output, dry_run=args.dry_run)

    required_tools = ['mafft', 'minimap2', 'viral_consensus', 'refineC', 'blastn', 'makeblastdb', 'pvga', 'gmcloser', 'abyss-sealer', 'samtools']
    if shutil.which('bbmerge.sh') is None and shutil.which('reformat.sh') is None: required_tools.append('bbmap(缺失)')
    missing = check_dependencies(required_tools)
    if missing: sys.exit(logger.error(f"[致命错误] 缺失依赖: {', '.join(missing)}"))
    if not CAN_PLOT: logger.warning("[提示] 未安装 matplotlib，将跳过绘图模块。")

    active_tools = {'megahit', 'spades', 'penguin'} if 'all' in args.assembly_tools.lower() else set(t.strip().lower() for t in args.assembly_tools.split(','))
    orig_ref = os.path.abspath(args.reference); _, ref_seq_str = read_first_sequence(orig_ref)
    if not ref_seq_str:
        logger.error(f"[致命错误] 参考序列为空或无法读取: {orig_ref}")
        sys.exit(1)
    ref_len = len(ref_seq_str)
  
    asm_dict = parse_reads_input(args.assembly_reads, logger)
    pvga_dict = parse_reads_input(args.pvga_reads, logger) if args.pvga_reads else {}
    cons_dict = parse_reads_input(args.consensus_reads, logger) if args.consensus_reads else {}
    master_samples = {s: {'assembly': a, 'pvga': pvga_dict.get(s, a), 'consensus': cons_dict.get(s, a)} for s, a in asm_dict.items()}

    if args.dry_run: sys.exit(0)

    out_dir = Path(os.path.abspath(args.output))
    logger.info("="*60)
    logger.info("OmniVirusAssembler 生产线已启动 (V9.0)".center(58))
    logger.info("="*60)

    all_stats = []
    with ThreadPoolExecutor(max_workers=args.jobs) as executor:
        futures = {executor.submit(run_pipeline, s, data, orig_ref, ref_len, out_dir, active_tools, args, logger): s for s, data in master_samples.items()}
        for future in tqdm(as_completed(futures), total=len(master_samples), desc="批处理进度", dynamic_ncols=True):
            s_name = futures[future]
            try:
                success, msg, _, stats = future.result()
                if success:
                    logger.info(f"✓ [{s_name}] {msg}")
                    for st in stats: st['Sample'] = s_name
                    all_stats.extend(stats)
            except Exception as exc: logger.error(f"✗ [{s_name}] 运行时崩溃: {str(exc)}")

    if all_stats:
        df = pd.DataFrame(all_stats)
        cols = ['Sample'] + [c for c in df.columns if c != 'Sample']
        df[cols].to_csv(out_dir / ("Global_Evolution_Stats.tsv" if ADVANCED_STATS else "Global_Evolution_Stats.csv"), sep='\t' if ADVANCED_STATS else ',', index=False)
        logger.info(f"\n[大功告成] 整体指标已汇总保存。")

if __name__ == "__main__":
    main()
