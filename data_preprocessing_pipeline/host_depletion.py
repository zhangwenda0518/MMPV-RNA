#!/usr/bin/env python3
"""
终极独立版去宿主 + 去rRNA混合管道 (Kraken2 初筛 + 精筛 + Ribodetector)
【数据可视化版】集成 seqkit stats 与绘图模块，智能识别 rna-short，自适应报表目录结构。
*支持 FASTA/FASTQ 及其 gz 压缩格式自适应处理*
*新增：每样本资源使用统计（时间/内存/CPU）并生成汇总报表*
"""

import os
import re
import sys
import time
import glob
import shutil
import atexit
import signal
import logging
import argparse
import tempfile
import warnings
import subprocess
import concurrent.futures
from threading import Lock

# 尝试导入绘图依赖
try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import pandas as pd
    PLOT_AVAILABLE = True
except ImportError:
    PLOT_AVAILABLE = False

# ==========================================
# 0. 终端颜色 UI 与动态进度条
# ==========================================
class UI:
    CYAN = '\033[96m'
    GREEN = '\033[92m'
    YELLOW = '\033[93m'
    RED = '\033[91m'
    PURPLE = '\033[95m'
    GRAY = '\033[90m'
    BOLD = '\033[1m'
    RESET = '\033[0m'

    print_lock = Lock()
    total_tasks = 0
    completed_tasks = 0
    success_tasks = 0
    fail_tasks = 0

    @staticmethod
    def update_progress(success=True):
        with UI.print_lock:
            UI.completed_tasks += 1
            if success: UI.success_tasks += 1
            else: UI.fail_tasks += 1
            
            percent = int((UI.completed_tasks / UI.total_tasks) * 100) if UI.total_tasks > 0 else 100
            bar_len = 30
            filled_len = int(bar_len * UI.completed_tasks // UI.total_tasks) if UI.total_tasks > 0 else bar_len
            bar = '█' * filled_len + '░' * (bar_len - filled_len)
            
            color = UI.GREEN if UI.fail_tasks == 0 else UI.YELLOW
            sys.stdout.write(f"\r\033[K{color}{UI.BOLD}进度: [{bar}] {percent}% ({UI.completed_tasks}/{UI.total_tasks}) | 成功: {UI.success_tasks} | 失败: {UI.fail_tasks}{UI.RESET}")
            sys.stdout.flush()

# ==========================================
# 1. 状态机、断点管理、Seqkit与绘图模块
# ==========================================
class CheckpointManager:
    def __init__(self, outdir, force=False):
        os.makedirs(outdir, exist_ok=True)
        self.filepath = os.path.join(outdir, '.checkpoints')
        self.completed = set()
        if force and os.path.exists(self.filepath): os.remove(self.filepath)
        if os.path.exists(self.filepath):
            with open(self.filepath, 'r') as f:
                self.completed = set(line.strip() for line in f if line.strip())

    def is_done(self, sample_id): return sample_id in self.completed
    def mark_done(self, sample_id):
        with UI.print_lock:
            self.completed.add(sample_id)
            with open(self.filepath, 'a') as f: f.write(sample_id + '\n')

class SeqkitStats:
    report_lock = Lock()

    @staticmethod
    def run(sample_name, stage, file_paths, threads, summary_tsv):
        valid_files = [f for f in file_paths if f and os.path.exists(f)]
        if not valid_files: return
        cmd = ['seqkit', 'stats', '-T', '-j', str(threads)] + valid_files
        try:
            res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True, check=True)
            lines = res.stdout.strip().split('\n')
            if len(lines) < 2: return
            with SeqkitStats.report_lock:
                write_header = not os.path.exists(summary_tsv)
                with open(summary_tsv, 'a', encoding='utf-8') as f:
                    header = lines[0].split('\t')
                    if write_header: f.write("Sample\tStage\t" + "\t".join(header) + "\n")
                    for line in lines[1:]: f.write(f"{sample_name}\t{stage}\t{line}\n")
        except Exception as e:
            logging.error(f"[{sample_name}] Seqkit 统计在 {stage} 阶段失败: {e}")

def get_avg_length(file_path, threads=1):
    cmd = f"seqkit stats -T -j {threads} {file_path}"
    try:
        res = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        lines = res.stdout.strip().split('\n')
        if len(lines) > 1:
            header = lines[0].split('\t')
            data = lines[1].split('\t')
            if 'avg_len' in header:
                idx = header.index('avg_len')
                return float(data[idx])
    except Exception: pass
    return 100.0

def remove_ansi_escape(text):
    ansi_escape = re.compile(r'\x1B(?:[@-Z\\-_]|\[[0-?]*[ -/]*[@-~])')
    return ansi_escape.sub('', text)

def generate_rrna_report(tool_logs_dir, output_path):
    results = []
    for fname in os.listdir(tool_logs_dir):
        if fname.endswith('_ribodetector.log'):
            sample = fname.replace('_ribodetector.log', '')
            total, non_rrna, rrna = 0, 0, 0
            with open(os.path.join(tool_logs_dir, fname), 'r') as f:
                for line in f:
                    clean_line = remove_ansi_escape(line)
                    if 'Processed' in clean_line and 'sequences in total' in clean_line:
                        m = re.search(r'Processed (\d+) sequences in total', clean_line)
                        if m: total = int(m.group(1))
                    elif 'Detected' in clean_line and 'non-rRNA sequences' in clean_line:
                        m = re.search(r'Detected (\d+) non-rRNA sequences', clean_line)
                        if m: non_rrna = int(m.group(1))
                    elif 'Detected' in clean_line and 'rRNA sequences' in clean_line:
                        m = re.search(r'Detected (\d+) rRNA sequences', clean_line)
                        if m: rrna = int(m.group(1))
            if total > 0:
                results.append({
                    'sample': sample, 'total': total, 'non_rRNA': non_rrna, 'rRNA': rrna,
                    'non_rRNA_percent': (non_rrna / total * 100), 'rRNA_percent': (rrna / total * 100)
                })
    if results:
        with open(output_path, 'w') as f:
            f.write("Sample\tTotal_sequences\tnon_rRNA\trRNA\tnon_rRNA(%)\trRNA(%)\n")
            for res in sorted(results, key=lambda x: x['sample']):
                f.write(f"{res['sample']}\t{res['total']}\t{res['non_rRNA']}\t{res['rRNA']}\t{res['non_rRNA_percent']:.2f}%\t{res['rRNA_percent']:.2f}%\n")
        logger.info(f"📊 rRNA 去除统计报告已生成: {output_path}")

# --- 绘图函数整合 ---
def plot_data(pivot_df, plot_type, output_file, figsize=(10, 6), dpi=100, max_labels=50, million=False):
    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
    ylabel = 'Read Count' + (' (millions)' if million else '')
    
    num_samples = len(pivot_df)
    too_many_samples = num_samples > max_labels

    if plot_type == 'line':
        pivot_df.T.plot(kind='line', marker='o', ax=ax, legend=not too_many_samples)
        ax.set_xlabel('Stage')
        ax.set_ylabel(ylabel)
        ax.set_title('Read count changes across clean stages')
        if not too_many_samples:
            ax.legend(title='Sample', bbox_to_anchor=(1.05, 1), loc='upper left')
        ax.grid(True, linestyle='--', alpha=0.6)
        plt.xticks(rotation=45, ha='right')
    else:  # bar
        pivot_df.plot(kind='bar', ax=ax)
        ax.set_xlabel('Sample')
        ax.set_ylabel(ylabel)
        ax.set_title('Read count at each clean stage by sample')
        ax.legend(title='Stage', bbox_to_anchor=(1.05, 1), loc='upper left')
        ax.grid(axis='y', linestyle='--', alpha=0.6)
        if too_many_samples: ax.set_xticks([])
        else: plt.xticks(rotation=45, ha='right')

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        plt.tight_layout()
        
    plt.savefig(output_file, dpi=dpi, bbox_inches='tight')
    plt.close(fig)

def generate_plots(tsv_file, output_prefix):
    try:
        df = pd.read_csv(tsv_file, sep='\t', header=0)
        required_cols = ['Sample', 'Stage', 'num_seqs']
        if not all(col in df.columns for col in required_cols):
            logger.warning(f"缺少必要的列 {required_cols}，跳过绘图。")
            return
        
        df = df[required_cols].copy()
        df_grouped = df.groupby(['Sample', 'Stage'], as_index=False)['num_seqs'].mean()
        pivot_df = df_grouped.pivot(index='Sample', columns='Stage', values='num_seqs')
        stages = sorted(pivot_df.columns.tolist())
        pivot_df = pivot_df[stages]

        if pivot_df.empty: return

        line_out = output_prefix + '_line.png'
        bar_out = output_prefix + '_bar.png'
        
        plot_data(pivot_df, 'line', line_out)
        plot_data(pivot_df, 'bar', bar_out)
        logger.info(f"📈 可视化图表已生成: {line_out} 和 {bar_out}")
    except Exception as e:
        logger.error(f"绘图失败: {e}")

# ==========================================
# 2. 日志与基础配置
# ==========================================
def configure_logger(name, debug=False, quiet=False, filename=None):
    logger = logging.getLogger(name)
    if logger.hasHandlers(): logger.handlers.clear()
    format_style = '%(asctime)s [%(levelname)s] %(message)s'
    log_level = logging.DEBUG if debug else (logging.WARNING if quiet else logging.INFO)
    if filename: logging.basicConfig(format=format_style, level=log_level, filename=filename, filemode='a')
    else: logging.basicConfig(format=format_style, level=log_level)
    return logger

logger = logging.getLogger(__name__)

# ==========================================
# 3. 智能自动配对引擎 
# ==========================================
def auto_pair_files(file_list, outdir):
    os.makedirs(outdir, exist_ok=True)
    unpaired = set([os.path.abspath(f) for f in file_list])
    tasks = []
    pe_pattern = re.compile(r'^(.*?)([._-])(R1|r1|1)([._](?:fastq|fq|fasta|fa)(?:\.gz)?)$', re.IGNORECASE)
    
    for f in list(unpaired):
        if f not in unpaired: continue
        basename = os.path.basename(f)
        dirname = os.path.dirname(f)
        match = pe_pattern.match(basename)
        if match:
            prefix, sep, r1_tag, suffix = match.groups()
            r2_tag = 'R2' if r1_tag.upper() == 'R1' else '2'
            if r1_tag.islower(): r2_tag = r2_tag.lower()
            r2_basename = f"{prefix}{sep}{r2_tag}{suffix}"
            r2_file = os.path.join(dirname, r2_basename)
            if r2_file in unpaired:
                if not suffix.lower().endswith('.gz'): suffix += '.gz'
                clean_r1 = os.path.join(outdir, f"{prefix}_clean{sep}{r1_tag}{suffix}")
                clean_r2 = os.path.join(outdir, f"{prefix}_clean{sep}{r2_tag}{suffix}")
                tasks.append( ([f, r2_file], [clean_r1, clean_r2]) )
                unpaired.remove(f)
                unpaired.remove(r2_file)
                continue
    
    for f in sorted(unpaired):
        basename = os.path.basename(f)
        match_se = re.match(r'^(.*?)([._](?:fastq|fq|fasta|fa)(?:\.gz)?)$', basename, re.IGNORECASE)
        if match_se:
            prefix, suffix = match_se.groups()
            if not suffix.lower().endswith('.gz'): suffix += '.gz'
            clean_out = os.path.join(outdir, f"{prefix}_clean{suffix}")
        else:
            if not basename.lower().endswith('.gz'): basename += '.gz'
            clean_out = os.path.join(outdir, f"{basename}_clean")
        tasks.append( ([f], [clean_out]) )
    return tasks

def extract_sample_name(filename):
    basename = os.path.basename(filename)
    match = re.match(r'^(.*?)([._-])(R1|r1|1)([._](?:fastq|fq|fasta|fa)(?:\.gz)?)$', basename, re.IGNORECASE)
    if match: return match.group(1)
    else: return re.sub(r'[._](?:fastq|fq|fasta|fa)(?:\.gz)?$', '', basename, flags=re.IGNORECASE)

def check_dependencies(tools):
    missing = [tool for tool in tools if shutil.which(tool) is None]
    if missing:
        logger.error(f"环境依赖缺失！找不到以下命令: {', '.join(missing)}")
        sys.exit(1)

def verify_indices(kr2_idx, step2_tool, step2_idx):
    for f in ['hash.k2d', 'opts.k2d', 'taxo.k2d']:
        if not os.path.isfile(os.path.join(kr2_idx, f)):
            logger.error(f"Kraken2 数据库缺失: {f}")
            sys.exit(1)
    if not step2_idx:
        return  # --steps 不含 align，无需验证第二步索引
    if step2_tool == 'bowtie2' and not (glob.glob(f"{step2_idx}.1.bt2") or glob.glob(f"{step2_idx}.1.bt2l")):
        logger.error(f"找不到 Bowtie2 索引: {step2_idx}")
        sys.exit(1)
    elif step2_tool == 'hisat2' and not (glob.glob(f"{step2_idx}.1.ht2") or glob.glob(f"{step2_idx}.1.ht2l")):
        logger.error(f"找不到 HISAT2 索引: {step2_idx}")
        sys.exit(1)
    elif step2_tool == 'minimap2' and not (os.path.isfile(step2_idx) and step2_idx.endswith(('.mmi', '.fa', '.fasta', '.fna', '.gz'))):
        logger.error(f"找不到 Minimap2 索引: {step2_idx}")
        sys.exit(1)

def validate_args(args):
    for arg in args:
        if type(arg) == str:
            chars = re.findall(r'[^a-zA-Z0-9\. !#=/_-]', str(arg))
            if not re.match(r'^[a-zA-Z0-9\. !#=/_-]+$|^(?![\s\S])', str(arg)): return False, arg, chars
    return True, '', []

# ==========================================
# 4. 管道执行引擎 (增强 CPU 统计)
# ==========================================
resource_lock = Lock()

def execute(cmds, pipe=False, log_prefix="", logs_dir="", sample_id="", step_name=""):
    logger.debug(f'{log_prefix} 执行命令流: {cmds}')
    returncodes, procs = [], []
    stdin = None
    log_file_handles = []

    os.makedirs(logs_dir, exist_ok=True)

    for i, cmd in enumerate(cmds):
        tool_name = os.path.basename(cmd[0])
        log_file_path = os.path.join(logs_dir, f"{sample_id}.{step_name}.{i}_{tool_name}.log")
        err_f = open(log_file_path, 'w', encoding='utf-8')
        log_file_handles.append((tool_name, err_f, log_file_path))

        err_f.write(f"=== 时间: {time.strftime('%Y-%m-%d %H:%M:%S')} ===\n")
        err_f.write(f"=== 命令: {' '.join(cmd)} ===\n\n")
        err_f.flush()

        if pipe and i < len(cmds) - 1: out_target = subprocess.PIPE
        else: out_target = subprocess.DEVNULL
            
        proc = subprocess.Popen(cmd, stdin=stdin, stdout=out_target, stderr=err_f)
        if stdin: stdin.close()
        stdin = proc.stdout if pipe else None
        procs.append(proc)

    peak_mem_mb = 0
    # CPU 统计：每个进程的最大 CPU jiffies 总和
    pid_max_cpu = {}   # pid -> max jiffies
    # 初始化
    for proc in procs:
        pid_max_cpu[proc.pid] = 0

    while True:
        all_done = True
        current_mem = 0
        for proc in procs:
            if proc.poll() is None:
                all_done = False
                # Linux 下通过 /proc 统计内存和 CPU
                if sys.platform == 'linux':
                    try:
                        with open(f"/proc/{proc.pid}/statm", 'r') as f:
                            current_mem += int(f.read().split()[1]) * 4 / 1024
                    except: pass
                    # 读取 CPU 时间
                    try:
                        with open(f"/proc/{proc.pid}/stat", 'r') as f:
                            fields = f.read().split()
                            # 字段索引 13: utime, 14: stime (0-based)
                            if len(fields) > 14:
                                utime = int(fields[13])
                                stime = int(fields[14])
                                jiffies = utime + stime
                                if jiffies > pid_max_cpu.get(proc.pid, 0):
                                    pid_max_cpu[proc.pid] = jiffies
                    except: pass
                else:
                    # 非 Linux 不做统计，CPU 时间返回 0
                    pass
        if current_mem > peak_mem_mb:
            peak_mem_mb = current_mem
        if all_done: break
        time.sleep(0.5)

    for proc in procs:
        proc.wait()
        returncodes.append(proc.returncode)

    # 计算总 CPU 时间（秒），假设 jiffies 时钟频率为 100 Hz
    total_jiffies = sum(pid_max_cpu.values())
    if sys.platform == 'linux':
        try:
            clk_tck = os.sysconf('SC_CLK_TCK')
        except:
            clk_tck = 100
        total_cpu_seconds = total_jiffies / clk_tck
    else:
        total_cpu_seconds = 0.0

    for i, c in enumerate(returncodes):
        tool_name, err_f, log_file_path = log_file_handles[i]
        err_f.close()
        if c != 0:
            logger.error(f"{log_prefix} ❌ {tool_name} 失败 (Exit: {c}) 日志: {log_file_path}")
            try:
                with open(log_file_path, 'r') as f:
                    lines = [line.strip() for line in f.readlines() if line.strip()]
                    if lines: logger.error(f"{log_prefix} 🔍 [探针] -> {' | '.join(lines[-3:])}")
            except: pass

    return returncodes, peak_mem_mb, total_cpu_seconds

class Aligners:
    @staticmethod
    def kraken2(index, seq1, classified_out=None, unclassified_out=None, seq2=None, threads=1, confidence=0.4, report_out=None, output=None, memory_mapping=True):
        # --memory-mapping: kraken2 以 mmap 方式访问只读索引, 多进程经 page cache 共享同一份 DB。
        # 这样 --jobs>1 并行时内存只占一份 (~DB 大小), 而不是每个样本各把 DB 读进堆内存 (223G/进程)。
        cmd = ['kraken2', '--threads', str(threads), '--db', index, '--confidence', str(confidence)]
        if memory_mapping:
            cmd += ['--memory-mapping']
        if output: cmd += ['--output', output]
        else: cmd += ['--output', '/dev/null']
        if report_out: cmd += ['--report', report_out]
        if classified_out: cmd += ['--classified-out', classified_out]
        if unclassified_out: cmd += ['--unclassified-out', unclassified_out]
        if seq2: cmd += ['--paired', seq1, seq2]
        else: cmd += [seq1]
        return cmd

    @staticmethod
    def bowtie2(index, seq1, seq2=None, threads=1, options=[], is_fasta=False):
        fmt_flag = '-f' if is_fasta else '-q'
        cmd = ['bowtie2', '-p', str(threads), '-x', index, fmt_flag] + options
        if seq2: cmd += ['-1', seq1, '-2', seq2]
        else: cmd += ['-U', seq1]
        return cmd

    @staticmethod
    def hisat2(index, seq1, seq2=None, threads=1, options=[], is_fasta=False):
        cmd = ['hisat2', '-p', str(threads), '-x', index]
        if is_fasta: cmd += ['-f']
        cmd += options
        if seq2: cmd += ['-1', seq1, '-2', seq2]
        else: cmd += ['-U', seq1]
        return cmd

    @staticmethod
    def minimap2(index, seq1, seq2=None, threads=1, seq_type='dna-short', options=[]):
        cmd = ['minimap2', '-t', str(threads), '-a']
        if not options:
            if seq_type == 'dna-short': cmd += ['-x', 'sr']
            elif seq_type == 'rna-short': cmd += ['-x', 'splice:sr']
            elif seq_type == 'nanopore': cmd += ['-x', 'map-ont']
            elif seq_type == 'pacbio': cmd += ['-x', 'map-pb']
        else: cmd += options
        cmd += [index, seq1]
        if seq2: cmd += [seq2]
        return cmd

class SAM:
    @staticmethod
    def view_filter(out_bam, threads=1, paired=True, mfilter=True):
        cmd = ['samtools', 'view', '-@', str(threads), '-b']
        if paired:
            if mfilter: cmd += ['-f', '12', '-F', '256']
            else: cmd += ['-F', '268']
        else:
            if mfilter: cmd += ['-f', '4', '-F', '256']
            else: cmd += ['-F', '260']
        cmd += ['-o', out_bam, '-']
        return cmd

    @staticmethod
    def sort_name(in_bam, out_bam, threads=1, tmp_prefix=None):
        cmd = ['samtools', 'sort', '-n', '-@', str(threads), '-o', out_bam]
        if tmp_prefix: cmd += ['-T', tmp_prefix]
        cmd += [in_bam]
        return cmd

    @staticmethod
    def seq_extract(in_bam, out1, out2=None, threads=1, is_fasta=False):
        extract_tool = 'fasta' if is_fasta else 'fastq'
        cmd = ['samtools', extract_tool, '-@', str(threads), '-n']
        if out1 and out2: cmd += ['-1', out1, '-2', out2, '-0', '/dev/null', '-s', '/dev/null']
        else: cmd += ['-0', out1]
        cmd += [in_bam]
        return cmd

def write_resource_usage(sample_id, elapsed, peak_mem_mb, cpu_time_seconds, avg_cpus, logs_dir):
    """写入每样本资源日志和全局资源汇总表"""
    resource_log_dir = os.path.join(logs_dir, "logs")
    os.makedirs(resource_log_dir, exist_ok=True)
    single_log = os.path.join(resource_log_dir, f"{sample_id}.resource.log")
    summary_tsv = os.path.join(logs_dir, "host_depletion_resource_usage.tsv")

    log_lines = [
        f"Sample: {sample_id}",
        f"Elapsed time (s): {elapsed:.2f}",
        f"Peak memory (MB): {peak_mem_mb:.2f}",
        f"Total CPU time (s): {cpu_time_seconds:.2f}",
        f"Average CPUs used: {avg_cpus:.2f}",
        f"Timestamp: {time.strftime('%Y-%m-%d %H:%M:%S')}"
    ]
    with open(single_log, 'w') as f:
        f.write('\n'.join(log_lines) + '\n')

    with resource_lock:
        write_header = not os.path.exists(summary_tsv)
        with open(summary_tsv, 'a', encoding='utf-8') as f:
            if write_header:
                f.write("Sample\tElapsed_Seconds\tPeak_Mem_MB\tCPU_Time_Seconds\tAvg_CPUs\n")
            f.write(f"{sample_id}\t{elapsed:.2f}\t{peak_mem_mb:.2f}\t{cpu_time_seconds:.2f}\t{avg_cpus:.2f}\n")

# ==========================================
# 4b. Kraken2 read extraction helpers (--keep-taxids)
# ==========================================
def process_kraken_output_line(line):
    """Parse a Kraken2 --output line, return (taxid_int, read_id_str)."""
    l_vals = line.split('\t')
    if len(l_vals) < 5:
        return -1, ''
    tax_id_str = l_vals[2]
    if "taxid" in tax_id_str:
        tax_id_str = tax_id_str.split("taxid ")[-1][:-1]
    read_id = l_vals[1]
    if tax_id_str == 'A':
        return 81077, read_id
    elif tax_id_str == 'U':
        return 0, read_id
    try:
        return int(tax_id_str), read_id
    except ValueError:
        return -1, read_id


def build_children_taxids(report_file, target_taxids):
    """Parse a Kraken2 report file, build taxonomy tree, return target_taxids + all descendant taxIDs."""
    class _Node:
        __slots__ = ('taxid', 'level_num', 'children', 'parent')
        def __init__(self, taxid, level_num, parent=None):
            self.taxid = taxid; self.level_num = level_num; self.children = []; self.parent = parent
        def add_child(self, node): self.children.append(node)

    target_set = set(target_taxids)
    if not os.path.exists(report_file):
        return target_set

    base_nodes = {}
    with open(report_file, 'r') as f:
        prev_node = None
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 5:
                continue
            try:
                taxid = int(parts[-3])
            except (ValueError, IndexError):
                try:
                    taxid = int(parts[-2])
                except (ValueError, IndexError):
                    continue
            if taxid == 0:
                continue
            spaces = len(parts[-1]) - len(parts[-1].lstrip(' '))
            level_num = spaces // 2

            if taxid == 1:
                root = _Node(taxid, level_num)
                prev_node = root
                if taxid in target_set:
                    base_nodes[taxid] = root
                continue

            while prev_node and level_num != (prev_node.level_num + 1):
                prev_node = prev_node.parent
            if prev_node:
                curr = _Node(taxid, level_num, parent=prev_node)
                prev_node.add_child(curr)
                prev_node = curr
                if taxid in target_set:
                    base_nodes[taxid] = curr

    for tid in base_nodes:
        stack = list(base_nodes[tid].children)
        while stack:
            node = stack.pop()
            target_set.add(node.taxid)
            stack.extend(node.children)

    if not base_nodes:
        logger.warning(f"Kraken2 report 中未找到目标 taxID {target_taxids}，将仅匹配 taxID 本身（可能漏掉子节点）")
    return target_set


def seqtk_extract_reads(fastq_in, read_ids, fastq_out, threads=2):
    """Use seqtk subseq to extract reads by ID from a (possibly gzipped) FASTQ."""
    if not read_ids:
        return False
    fd, id_file = tempfile.mkstemp(prefix="k2ids_", suffix=".txt")
    try:
        with os.fdopen(fd, 'w') as f:
            for rid in read_ids:
                f.write(f"{rid}\n{rid}/1\n{rid}/2\n")
        is_gz = fastq_in.endswith('.gz')
        if is_gz:
            # decompress on the fly, avoid bash process substitution
            cmd = ['seqtk', 'subseq', '/dev/stdin', id_file]
            with open(fastq_out, 'w') as out_f, \
                 subprocess.Popen(['pigz', '-dc', fastq_in], stdout=subprocess.PIPE, stderr=subprocess.DEVNULL) as pz:
                subprocess.run(cmd, stdin=pz.stdout, stdout=out_f, check=True, stderr=subprocess.DEVNULL)
        else:
            subprocess.run(['seqtk', 'subseq', fastq_in, id_file],
                          stdout=open(fastq_out, 'w'), check=True, stderr=subprocess.DEVNULL)
        return os.path.getsize(fastq_out) > 0
    except Exception:
        return False
    finally:
        if os.path.exists(id_file):
            os.remove(id_file)


class HybridPipeline:
    def run_single(self, task_id, step2_tool, seq_type, kr2_idx, step2_idx, seq1, out1, seq2=None, out2=None,
                   threads=1, mfilter=True, config='', confidence=0.4, tmp_dir=None, logs_dir=None,
                   do_rrna=False, keep_rrna=False, chunk_size=256,
                   rrna_tool='ribodetector', silva_index=None,
                   keep_taxids=None, memory_mapping=True):
        
        start_time = time.time()
        max_mem_used = 0.0
        total_cpu_time = 0.0  # 累计该样本所有步骤的子进程 CPU 时间
        
        sample_name = extract_sample_name(seq1)
        steps_str = os.environ.get('HOST_DEPLETION_STEPS', 'kraken2,align,rrna')
        steps_set = set(steps_str.split(','))
        do_kraken2 = 'kraken2' in steps_set
        do_align = 'align' in steps_set
        is_fasta = bool(re.search(r'[._](?:fasta|fa)(?:\.gz)?$', seq1, re.IGNORECASE))
        seq_ext_name = 'FastA' if is_fasta else 'FastQ'
        
        prefix = f"[任务 {task_id:02d} | {sample_name}]"
        logger.info(f"{prefix} 🚀 启动 ({'双端' if seq2 else '单端'} | {seq_ext_name})")
        
        valid, arg, chars = validate_args([seq1, seq2, out1, out2] + config.split())
        if not valid:
            logger.error(f"{prefix} ❌ 安全拦截：参数包含非法字符: '{arg}'")
            return 1

        extra_options = config.split() if config else []
        if step2_tool == 'bowtie2' and not extra_options: extra_options = ['--end-to-end']

        outdir_base = os.path.dirname(os.path.abspath(out1))
        
        tool_logs_dir = os.path.join(logs_dir, "logs")
        k2_report_dir = os.path.join(logs_dir, "kraken2_report")
        os.makedirs(tool_logs_dir, exist_ok=True)
        os.makedirs(k2_report_dir, exist_ok=True)

        report_file = os.path.join(k2_report_dir, f"{sample_name}_kraken2_report.txt")
        summary_tsv = os.path.join(logs_dir, "host_depletion_seqkit_summary.tsv")

        actual_out1 = out1[:-3] if out1.endswith('.gz') else out1
        actual_out2 = (out2[:-3] if out2.endswith('.gz') else out2) if out2 else None
        ext_tmp = '.fasta' if is_fasta else '.fastq'

        with tempfile.TemporaryDirectory(dir=tmp_dir) as temp_dir:
            # --- [Stage 1] Raw Data Stats ---
            logger.info(f"{prefix} -> [Stage 1] Seqkit 统计原始数据...")
            SeqkitStats.run(sample_name, "1_Raw", [seq1, seq2], threads, summary_tsv)

            # --- [Stage 2] Kraken2 (可选) ---
            if do_kraken2:
                keep_taxids_list = [int(t) for t in str(keep_taxids).split(',') if t.strip()] if keep_taxids else []
                do_keep_taxids = bool(keep_taxids_list)

                if do_keep_taxids:
                    logger.info(f"{prefix} -> [Stage 2] Kraken2 初筛 (Conf: {confidence}, 保留 taxIDs: {keep_taxids_list})...")
                    # PE: Kraken2 with --paired requires '#' placeholder → produces _1 / _2 files
                    # SE: plain filename, single file
                    if seq2:
                        kr2_unclass_out  = f"{temp_dir}/k2_unclass#.fastq"
                        kr2_class_out    = f"{temp_dir}/k2_class#.fastq"
                        # After Kraken2, actual files are k2_unclass_1.fastq / k2_unclass_2.fastq etc.
                        kr2_unclass_1 = f"{temp_dir}/k2_unclass_1.fastq"
                        kr2_unclass_2 = f"{temp_dir}/k2_unclass_2.fastq"
                        kr2_class_1   = f"{temp_dir}/k2_class_1.fastq"
                        kr2_class_2   = f"{temp_dir}/k2_class_2.fastq"
                    else:
                        kr2_unclass_out  = f"{temp_dir}/k2_unclass{ext_tmp}"
                        kr2_class_out    = f"{temp_dir}/k2_class{ext_tmp}"
                    kr2_output_file  = f"{temp_dir}/kraken2_output.txt"

                    kr2_cmd = Aligners.kraken2(kr2_idx, seq1,
                        classified_out=kr2_class_out,
                        unclassified_out=kr2_unclass_out,
                        seq2=seq2, threads=threads, confidence=confidence,
                        report_out=report_file, output=kr2_output_file,
                        memory_mapping=memory_mapping)
                    ret, mem, cpu_time = execute([kr2_cmd], pipe=False, log_prefix=prefix,
                        logs_dir=tool_logs_dir, sample_id=sample_name, step_name="step1_kraken2")
                    max_mem_used = max(max_mem_used, mem)
                    total_cpu_time += cpu_time
                    if any(c != 0 for c in ret): return 1

                    # --- Parse report to expand child taxIDs (e.g. 10239 → all virus descendants) ---
                    expanded_taxids = build_children_taxids(report_file, keep_taxids_list)
                    # 注意: 不把 taxID=0(unclassified) 加进 target。unclassified reads 不在 classified
                    # 池里, 加进来只会让 seqtk 为几百条病毒去空扫上百万个 id(极慢 + 吃内存)。
                    # unclassified reads 由后面 merge 直接整箱保留(--unclassified-out), 无需在此列出。
                    logger.info(f"{prefix}   展开后 target taxIDs: {len(expanded_taxids)} 个")

                    # --- Parse kraken2 output to collect read IDs for target taxIDs ---
                    target_read_ids = set()
                    with open(kr2_output_file, 'r') as f:
                        for line in f:
                            tid, rid = process_kraken_output_line(line)
                            if tid in expanded_taxids:
                                target_read_ids.add(rid)
                    n_viral = len(target_read_ids)
                    logger.info(f"{prefix}   Kraken2 classified pool 中检出目标 taxID reads: {n_viral}")

                    # --- Extract viral reads from classified FASTQ(s) with seqtk, then merge ---
                    if seq2:
                        if n_viral > 0:
                            kr2_viral_1 = f"{temp_dir}/k2_viral_1.fastq"
                            kr2_viral_2 = f"{temp_dir}/k2_viral_2.fastq"
                            # R1/R2 抽取并行 (各自单线程 seqtk, 互不依赖)
                            with concurrent.futures.ThreadPoolExecutor(max_workers=2) as _ex:
                                _f1 = _ex.submit(seqtk_extract_reads, kr2_class_1, target_read_ids, kr2_viral_1, threads)
                                _f2 = _ex.submit(seqtk_extract_reads, kr2_class_2, target_read_ids, kr2_viral_2, threads)
                                ok1, ok2 = _f1.result(), _f2.result()
                            if ok1 and ok2:
                                merged_1 = f"{temp_dir}/k2_merged_1.fastq"
                                merged_2 = f"{temp_dir}/k2_merged_2.fastq"
                                for uncl, viral, merged in [(kr2_unclass_1, kr2_viral_1, merged_1),
                                                           (kr2_unclass_2, kr2_viral_2, merged_2)]:
                                    with open(merged, 'wb') as mf:
                                        for i, src in enumerate([uncl, viral]):
                                            with open(src, 'rb') as sf:
                                                shutil.copyfileobj(sf, mf)
                                            if i == 0:
                                                mf.write(b'\n')  # ensure newline before appending second file
                                temp1, temp2 = merged_1, merged_2
                                logger.info(f"{prefix}   合并完成: unclassified + {n_viral} viral reads (PE)")
                            else:
                                temp1, temp2 = kr2_unclass_1, kr2_unclass_2
                                logger.warning(f"{prefix}   seqtk 提取失败，仅保留 unclassified reads")
                            # Clean up classified and viral temp files
                            for _f in [kr2_class_1, kr2_class_2, kr2_viral_1, kr2_viral_2, kr2_output_file]:
                                if os.path.exists(_f): os.remove(_f)
                        else:
                            temp1, temp2 = kr2_unclass_1, kr2_unclass_2
                            for _f in [kr2_class_1, kr2_class_2, kr2_output_file]:
                                if os.path.exists(_f): os.remove(_f)
                            logger.info(f"{prefix}   目标 taxID reads 为空，仅保留 unclassified")
                    else:
                        # SE path
                        if n_viral > 0:
                            kr2_viral_out = f"{temp_dir}/k2_viral{ext_tmp}"
                            ok = seqtk_extract_reads(kr2_class_out, target_read_ids, kr2_viral_out, threads=threads)
                            if ok:
                                merged_out = f"{temp_dir}/k2_merged{ext_tmp}"
                                with open(merged_out, 'wb') as mf:
                                    for src in [kr2_unclass_out, kr2_viral_out]:
                                        with open(src, 'rb') as sf:
                                            shutil.copyfileobj(sf, mf)
                                temp1, temp2 = merged_out, None
                                logger.info(f"{prefix}   合并完成: unclassified + {n_viral} viral reads")
                            else:
                                temp1, temp2 = kr2_unclass_out, None
                                logger.warning(f"{prefix}   seqtk 提取失败，仅保留 unclassified reads")
                            for _f in [kr2_class_out, kr2_viral_out, kr2_output_file]:
                                if os.path.exists(_f): os.remove(_f)
                        else:
                            temp1, temp2 = kr2_unclass_out, None
                            for _f in [kr2_class_out, kr2_output_file]:
                                if os.path.exists(_f): os.remove(_f)
                            logger.info(f"{prefix}   目标 taxID reads 为空，仅保留 unclassified")
                else:
                    logger.info(f"{prefix} -> [Stage 2] Kraken2 初筛 (Conf: {confidence})...")
                    kr2_out = f"{temp_dir}/out#{ext_tmp}" if seq2 and actual_out2 else f"{temp_dir}/out_1{ext_tmp}"
                    class_out, unclass_out = (None, kr2_out) if mfilter else (kr2_out, None)

                    kr2_cmd = Aligners.kraken2(kr2_idx, seq1, classified_out=class_out, unclassified_out=unclass_out, seq2=seq2, threads=threads, confidence=confidence, report_out=report_file, memory_mapping=memory_mapping)
                    ret, mem, cpu_time = execute([kr2_cmd], pipe=False, log_prefix=prefix, logs_dir=tool_logs_dir, sample_id=sample_name, step_name="step1_kraken2")
                    max_mem_used = max(max_mem_used, mem)
                    total_cpu_time += cpu_time
                    if any(c != 0 for c in ret): return 1

                    temp1, temp2 = f"{temp_dir}/out_1{ext_tmp}", f"{temp_dir}/out_2{ext_tmp}" if seq2 else None

                SeqkitStats.run(sample_name, "2_Kraken2_Filtered", [temp1, temp2], threads, summary_tsv)
            else:
                temp1, temp2 = seq1, seq2
                logger.info(f"{prefix} -> [Stage 2] Kraken2 跳过，原始数据直接进入下一步...")

            # --- [Stage 4] Aligner Filter (可选) ---
            if do_align:
                logger.info(f"{prefix} -> [Stage 4] {step2_tool} 精筛并内存直滤为小体积 BAM...")
                filtered_bam = os.path.join(temp_dir, f"{sample_name}_filtered.bam")
                if step2_tool == 'bowtie2': align_cmd = Aligners.bowtie2(step2_idx, temp1, seq2=temp2, threads=threads, options=extra_options, is_fasta=is_fasta)
                elif step2_tool == 'hisat2': align_cmd = Aligners.hisat2(step2_idx, temp1, seq2=temp2, threads=threads, options=extra_options, is_fasta=is_fasta)
                elif step2_tool == 'minimap2': align_cmd = Aligners.minimap2(step2_idx, temp1, seq2=temp2, threads=threads, seq_type=seq_type, options=extra_options)
                view_cmd = SAM.view_filter(out_bam=filtered_bam, threads=threads, paired=bool(seq2), mfilter=mfilter)

                ret, mem, cpu_time = execute([align_cmd, view_cmd], pipe=True, log_prefix=prefix, logs_dir=tool_logs_dir, sample_id=sample_name, step_name="step2_align_filter")
                max_mem_used = max(max_mem_used, mem)
                total_cpu_time += cpu_time
                if any(c != 0 for c in ret): return 1

                # --- [Stage 5] Sort ---
                target_bam = filtered_bam
                if seq2:
                    logger.info(f"{prefix} -> [Stage 5] BAM Name-Sort (确保双端精确提取)...")
                    sorted_bam = os.path.join(temp_dir, f"{sample_name}_sorted.bam")
                    sort_cmd = SAM.sort_name(in_bam=filtered_bam, out_bam=sorted_bam, threads=threads, tmp_prefix=os.path.join(temp_dir, "sort_tmp"))
                    ret, mem, cpu_time = execute([sort_cmd], pipe=False, log_prefix=prefix, logs_dir=tool_logs_dir, sample_id=sample_name, step_name="step3_sort")
                    max_mem_used = max(max_mem_used, mem)
                    total_cpu_time += cpu_time
                    if any(c != 0 for c in ret): return 1
                    target_bam = sorted_bam
                else:
                    logger.info(f"{prefix} -> [Stage 5] 单端数据，跳过 Name-Sort...")

                # --- [Stage 6] Sequence Extract ---
                extract_tgt1 = f"{temp_dir}/{sample_name}_hostclean.1{ext_tmp}" if do_rrna else actual_out1
                extract_tgt2 = (f"{temp_dir}/{sample_name}_hostclean.2{ext_tmp}" if do_rrna else actual_out2) if seq2 else None

                extract_cmd = SAM.seq_extract(in_bam=target_bam, out1=extract_tgt1, out2=extract_tgt2, threads=threads, is_fasta=is_fasta)
                ret, mem, cpu_time = execute([extract_cmd], pipe=False, log_prefix=prefix, logs_dir=tool_logs_dir, sample_id=sample_name, step_name="step4_extract")
                max_mem_used = max(max_mem_used, mem)
                total_cpu_time += cpu_time
                if any(c != 0 for c in ret): return 1

                SeqkitStats.run(sample_name, "3_Host_Filtered", [extract_tgt1, extract_tgt2], threads, summary_tsv)
            else:
                extract_tgt1 = f"{temp_dir}/{sample_name}_hostclean.1{ext_tmp}" if do_rrna else actual_out1
                extract_tgt2 = (f"{temp_dir}/{sample_name}_hostclean.2{ext_tmp}" if do_rrna else actual_out2) if seq2 else None
                shutil.copy2(temp1, extract_tgt1)
                if seq2 and extract_tgt2: shutil.copy2(temp2, extract_tgt2)
                SeqkitStats.run(sample_name, "3_Host_Filtered", [extract_tgt1, extract_tgt2], threads, summary_tsv)

            # --- [Stage 7] rRNA Depletion (可选) ---
            do_rrna_step = do_rrna and ('rrna' in steps_set)
            if do_rrna_step:
                if rrna_tool == 'silva' and silva_index:
                    # ── SILVA 模式: Bowtie2 比对 SILVA 库, --un 直接输出未比对 reads ──
                    logger.info(f"{prefix} -> [Stage 7] SILVA Bowtie2 剔除 rRNA...")

                    # --un-conc-gz 会在文件名后加 .1/.2 后缀, 需先写临时前缀再 rename
                    silva_tmp_prefix = os.path.join(temp_dir, f"{sample_name}_silva_unmapped")
                    silva_cmd = ["bowtie2", "--very-sensitive-local", "-p", str(threads),
                                 "-x", silva_index]
                    if seq2:
                        silva_cmd += ["-1", extract_tgt1, "-2", extract_tgt2,
                                      "--un-conc-gz", silva_tmp_prefix]
                    else:
                        silva_cmd += ["-U", extract_tgt1,
                                      "--un-gz", actual_out1]
                    ret, mem, cpu_time = execute([silva_cmd], pipe=True, log_prefix=prefix,
                                                  logs_dir=tool_logs_dir, sample_id=sample_name,
                                                  step_name="step5_silva_bowtie2")
                    max_mem_used = max(max_mem_used, mem)
                    total_cpu_time += cpu_time
                    if any(c != 0 for c in ret): return 1

                    # bowtie2 --un-conc-gz 产出 silva_tmp_prefix.1 / silva_tmp_prefix.2
                    # rename 到 actual_out1 / actual_out2
                    if seq2:
                        os.rename(silva_tmp_prefix + ".1", actual_out1)
                        os.rename(silva_tmp_prefix + ".2", actual_out2)

                    # keep_rrna: 重新比对并保留 rRNA reads
                    if keep_rrna:
                        rrna_dir = os.path.join(outdir_base, "rrna")
                        os.makedirs(rrna_dir, exist_ok=True)
                        rrna_out1 = os.path.join(rrna_dir, f"{sample_name}.silva_rrna.1{ext_tmp}") if seq2 else os.path.join(rrna_dir, f"{sample_name}.silva_rrna{ext_tmp}")
                        rrna_tmp_prefix = os.path.join(temp_dir, f"{sample_name}_silva_rrna_tmp")
                        rrna_cmd = ["bowtie2", "--very-sensitive-local", "-p", str(threads),
                                    "-x", silva_index]
                        if seq2:
                            rrna_cmd += ["-1", extract_tgt1, "-2", extract_tgt2,
                                         "--al-conc-gz", rrna_tmp_prefix]
                        else:
                            rrna_cmd += ["-U", extract_tgt1,
                                         "--al-gz", rrna_out1]
                        ret2, mem2, cpu2 = execute([rrna_cmd], pipe=True, log_prefix=prefix,
                                                    logs_dir=tool_logs_dir, sample_id=sample_name,
                                                    step_name="step5b_silva_rrna")
                        max_mem_used = max(max_mem_used, mem2)
                        total_cpu_time += cpu2
                        if seq2:
                            os.rename(rrna_tmp_prefix + ".1", rrna_out1)
                            os.rename(rrna_tmp_prefix + ".2", rrna_out2)

                else:
                    # ── Ribodetector 模式 (原有逻辑) ──
                    logger.info(f"{prefix} -> [Stage 7] Ribodetector 剔除 rRNA 序列...")
                    avg_len = get_avg_length(extract_tgt1, threads)
                    length_threshold = int(avg_len * 0.5) if avg_len > 0 else 50
                    logger.debug(f"{prefix} 推测平均读长 {avg_len}，设定过滤阈值 {length_threshold}")

                    if keep_rrna:
                        rrna_dir = os.path.join(outdir_base, "rrna")
                        os.makedirs(rrna_dir, exist_ok=True)
                        rrna_out1 = os.path.join(rrna_dir, f"{sample_name}.rrna.1{ext_tmp}") if seq2 else os.path.join(rrna_dir, f"{sample_name}.rrna{ext_tmp}")
                        rrna_out2 = os.path.join(rrna_dir, f"{sample_name}.rrna.2{ext_tmp}") if seq2 else None
                    else:
                        rrna_out1 = os.path.join(temp_dir, f"{sample_name}.rrna.1{ext_tmp}") if seq2 else os.path.join(temp_dir, f"{sample_name}.rrna{ext_tmp}")
                        rrna_out2 = os.path.join(temp_dir, f"{sample_name}.rrna.2{ext_tmp}") if seq2 else None

                    rrna_log = os.path.join(tool_logs_dir, f"{sample_name}_ribodetector.log")
                    ribo_cmd = ["ribodetector_cpu", "-t", str(threads), "-l", str(length_threshold), "-i", extract_tgt1]
                    if seq2: ribo_cmd.append(extract_tgt2)
                    ribo_cmd.extend(["-e", "rrna", "--chunk_size", str(chunk_size), "--log", rrna_log, "-o", actual_out1])
                    if seq2: ribo_cmd.append(actual_out2)
                    ribo_cmd.extend(["-r", rrna_out1])
                    if seq2: ribo_cmd.append(rrna_out2)

                    ret, mem, cpu_time = execute([ribo_cmd], pipe=False, log_prefix=prefix, logs_dir=tool_logs_dir, sample_id=sample_name, step_name="step5_ribodetector")
                    max_mem_used = max(max_mem_used, mem)
                    total_cpu_time += cpu_time
                    if any(c != 0 for c in ret): return 1

                    if keep_rrna and out1.endswith('.gz'):
                        zip_tool = 'pigz' if shutil.which('pigz') else 'gzip'
                        files_to_zip = [rrna_out1]
                        if seq2: files_to_zip.append(rrna_out2)
                        ret, mem, cpu_time = execute([[zip_tool, '-f'] + files_to_zip], pipe=False, log_prefix=prefix, logs_dir=tool_logs_dir, sample_id=sample_name, step_name="step6_compress_rrna")
                        max_mem_used = max(max_mem_used, mem)
                        total_cpu_time += cpu_time

                SeqkitStats.run(sample_name, "4_rRNA_Filtered", [actual_out1, actual_out2], threads, summary_tsv)

            # --- [Stage 8] Compress Main output ---
            if out1.endswith('.gz'):
                logger.info(f"{prefix} -> GZ 压缩最终序列...")
                files_to_compress = [actual_out1]
                if actual_out2: files_to_compress.append(actual_out2)
                zip_tool = 'pigz' if shutil.which('pigz') else 'gzip'
                cmd_zip = [zip_tool, '-f']
                if zip_tool == 'pigz': cmd_zip += ['-p', str(threads)]
                cmd_zip += files_to_compress
                ret, mem, cpu_time = execute([cmd_zip], pipe=False, log_prefix=prefix, logs_dir=tool_logs_dir, sample_id=sample_name, step_name="step7_compress_main")
                max_mem_used = max(max_mem_used, mem)
                total_cpu_time += cpu_time
                if any(c != 0 for c in ret): return 1

        elapsed = time.time() - start_time
        avg_cpus = total_cpu_time / elapsed if elapsed > 0 else 0.0

        # 写入资源使用记录
        write_resource_usage(sample_name, elapsed, max_mem_used, total_cpu_time, avg_cpus, logs_dir)

        logger.info(f"{prefix} ✅ 成功完成 (耗时: {elapsed:.2f} 秒 | 峰值内存: {max_mem_used:.2f} MB | CPU时间: {total_cpu_time:.2f}秒 [均{avg_cpus:.2f}核])")
        return 0

# ==========================================
# 5. 任务调度与解析
# ==========================================
def main():
    parser = argparse.ArgumentParser(description="智能序列去宿主 + 去 rRNA 混合管道")
    
    parser.add_argument('--tool', choices=['bowtie2', 'hisat2', 'minimap2'], help="步骤2精细比对工具 (--steps 含 align 时必需)")
    parser.add_argument('--seq-type', choices=['dna-short', 'rna-short', 'nanopore', 'pacbio'], help="测序数据类型 (--steps 含 align 或 rrna 时必需)")
    parser.add_argument('-k', '--kraken2_index', required=True, help="Kraken2 宿主数据库")
    parser.add_argument('-x', '--step2_index', help="第二步工具索引前缀 (--steps 含 align 时必需)")
    
    group_in = parser.add_mutually_exclusive_group(required=True)
    group_in.add_argument('-I', '--input-dir', type=str, help="输入目录自动扫描配对")
    group_in.add_argument('-i', '--input', action='append', nargs='+', help="指定一个或多个输入文件")
    
    parser.add_argument('-O', '--outdir', type=str, help="统一输出存放目录")
    parser.add_argument('-o', '--output', action='append', nargs='+', help="手动指定输出路径")
    
    parser.add_argument('--jobs', type=int, default=1, help="并发处理样本数")
    parser.add_argument('-t', '--threads', type=int, default=4, help="每样本分配线程数")
    parser.add_argument('-T', '--tmp', type=str, help="临时文件存放位置")
    parser.add_argument('--logs_dir', type=str, default='host_depletion_logs', help="底层运行日志及报告保存目录")
    parser.add_argument('--force', action='store_true', help="强制重新运行所有样本")
    
    parser.add_argument('--confidence', type=float, default=0.4, help="Kraken2 分类置信度阈值")
    parser.add_argument('-f', '--filter', choices=['true', 'false'], default='true', help="true: 去除宿主; false: 提取宿主")
    
    parser.add_argument('--rrna', action='store_true', help="开启 rRNA 剔除")
    parser.add_argument('--rrna_tool', default='ribodetector', choices=['ribodetector', 'silva'],
                        help="rRNA 剔除工具: ribodetector (默认, 仅 rna-short) / silva (Bowtie2 比对 SILVA 库)")
    parser.add_argument('--silva_index', help="SILVA Bowtie2 索引前缀 (--rrna_tool silva 时必需)")
    parser.add_argument('--keep_rrna', action='store_true', help="将分离出来的rRNA序列保存到 rrna/ 文件夹中 (默认丢弃)")
    parser.add_argument('--keep-taxids', default='10239', help="Kraken2 中保留的 taxID (逗号分隔), 默认 10239=Viruses; 设为空则不保留任何分类reads")
    parser.add_argument('--vmtouch', action='store_true', help="用 vmtouch 预热 DB 到 page cache (配合 --memory-mapping, 多样本并行时推荐, 避免首样本冷启动)")
    parser.add_argument('--no-memory-mapping', action='store_true', help="禁用 --memory-mapping (每样本把 DB 读进堆内存; 仅单样本/调试, 多样本会爆内存)")
    parser.add_argument('--chunk_size', type=int, default=256, help="ribodetector_cpu chunk_size (默认: 256)")
    parser.add_argument('--rrna_report', type=str, default='ribodetector.report.txt', help="rRNA统计表名")
    
    parser.add_argument('--steps', default='kraken2,align,rrna',
                        choices=['kraken2', 'align', 'rrna', 'kraken2,align', 'kraken2,rrna',
                                 'align,rrna', 'kraken2,align,rrna'],
                        help="消融实验：指定执行哪些步骤 (default: kraken2,align,rrna)")
    parser.add_argument('-c', '--config', default='', help="透传额外参数")
    parser.add_argument('-d', '--debug', action='store_true', help="开启详细日志")

    args = parser.parse_args()

    # --- memory-mapping 模式: 默认开启, kraken2 mmap 共享只读索引, 支持多样本并行且内存仅占一份 ---
    memory_mapping = not args.no_memory_mapping

    # --- 智能步骤调整: --keep-taxids 启用时自动跳过 align ---
    # k2_pluspfp 等通用库已通过 Kraken2 + keep-taxids 完成宿主过滤和病毒保留
    if args.keep_taxids:
        steps_list = args.steps.split(',')
        if 'align' in steps_list:
            steps_list = [s for s in steps_list if s != 'align']
            args.steps = ','.join(steps_list) if steps_list else ''
            sys.stderr.write("💡 --keep-taxids 已启用，自动跳过 align 阶段 (通用库 Kraken2 已覆盖宿主过滤)\n")

    os.environ['HOST_DEPLETION_STEPS'] = args.steps
    global logger
    logger = configure_logger(__name__, debug=args.debug)

    # 进程守护: 停止主命令时连带杀掉所有子孙 (kraken2/seqtk/worker), 避免后台残留
    try:
        import process_guard
        process_guard.install(logger)
    except Exception:
        pass

    if args.tmp: os.makedirs(args.tmp, exist_ok=True)
    os.makedirs(args.logs_dir, exist_ok=True)

    tasks = []
    checkpoint_dir = "" 

    if args.input_dir:
        if not args.outdir:
            logger.error("使用目录扫描 (-I/--input-dir) 时，必须指定统一输出目录 (-O/--outdir)！")
            sys.exit(1)
        checkpoint_dir = args.outdir
        files = []
        for ext in ('*.fastq', '*.fq', '*.fastq.gz', '*.fq.gz',
                     '*.fasta', '*.fa', '*.fasta.gz', '*.fa.gz',
                     '*_fastq.gz', '*_fq.gz', '*_fasta.gz', '*_fa.gz'):
            files.extend(glob.glob(os.path.join(args.input_dir, ext)))
        if not files:
            logger.error(f"目录 {args.input_dir} 中未找到序列文件！")
            sys.exit(1)
        tasks = auto_pair_files(files, args.outdir)

    elif args.input:
        if args.outdir:
            checkpoint_dir = args.outdir
            flat_files = [f for sublist in args.input for f in sublist]
            tasks = auto_pair_files(flat_files, args.outdir)
        elif args.output:
            if len(args.input) != len(args.output): sys.exit(1)
            checkpoint_dir = os.path.dirname(os.path.abspath(args.output[0][0]))
            for ins, outs in zip(args.input, args.output): tasks.append((ins, outs))

    if not tasks: sys.exit(1)

    # 默认临时目录 = 输出目录，避免 /tmp 空间不足 (Kraken2 classified/unclassified 文件很大)
    if not args.tmp:
        args.tmp = checkpoint_dir
        os.makedirs(args.tmp, exist_ok=True)

    # rRNA 智能条件限制
    actual_do_rrna = False
    rrna_tool = 'ribodetector'
    silva_index = None
    if args.rrna:
        rrna_tool = args.rrna_tool
        if rrna_tool == 'silva':
            # SILVA 模式: Bowtie2 比对, 不限 seq_type
            if not args.silva_index:
                logger.error("--rrna_tool silva 需要 --silva_index (SILVA Bowtie2 索引)")
                sys.exit(1)
            actual_do_rrna = True
            silva_index = args.silva_index
        elif args.seq_type == 'rna-short':
            actual_do_rrna = True
        else:
            logger.warning(f"检测到 '--seq-type {args.seq_type}'，ribodetector 仅支持 rna-short. 已自动禁用. 可改用 --rrna_tool silva")
            actual_do_rrna = False

    deps = ['kraken2', 'samtools', 'seqkit']
    if 'align' in args.steps.split(','):
        if not args.tool:
            logger.error("--steps 含 align 时必须提供 --tool"); sys.exit(1)
        deps.append(args.tool)
        if not args.step2_index:
            logger.error("--steps 含 align 时必须提供 -x/--step2_index"); sys.exit(1)
    if actual_do_rrna and rrna_tool == 'ribodetector':
        if not args.seq_type:
            logger.error("--steps 含 rrna (ribodetector) 时必须提供 --seq-type"); sys.exit(1)
        deps.append('ribodetector_cpu')
    elif actual_do_rrna and rrna_tool == 'silva':
        deps.append('bowtie2')
    if args.keep_taxids:
        deps.append('seqtk')
    check_dependencies(deps)
    if 'align' in args.steps.split(',') and not args.step2_index:
        logger.error("--steps 含 align 时必须提供 -x/--step2_index")
        sys.exit(1)
    verify_indices(args.kraken2_index, args.tool, args.step2_index)

    cm = CheckpointManager(checkpoint_dir, force=args.force)

    pe_count = sum(1 for t in tasks if len(t[0]) == 2)
    se_count = sum(1 for t in tasks if len(t[0]) == 1)
    mfilter = True if args.filter == 'true' else False

    print("\n")
    logger.info("============== 智能序列清理管道启动 ==============")
    logger.info(f"数据类型: {args.seq_type.upper() if args.seq_type else 'N/A'} | 第2步工具: {args.tool.upper() if args.tool else 'N/A (已跳过)'}")
    tool_label = "SILVA/Bowtie2" if rrna_tool == 'silva' else "ribodetector_cpu"
    if actual_do_rrna: logger.info(f"rRNA过滤: 已开启 (使用 {tool_label})")
    logger.info(f"报表与图表存储: {args.logs_dir}/")
    logger.info(f"解析任务: 包含 {pe_count} 个双端，{se_count} 个单端。")
    logger.info(f"并发配置: {args.jobs} 并发任务 × {args.threads} 线程/任务")
    logger.info("="*60)

    # --- memory-mapping: kraken2 以 mmap 访问只读索引, --jobs>1 时全体样本共享 page cache 中的同一份 DB ---
    # 可选 vmtouch 预热: 提前把 .k2d 页读进 page cache, 避免前几个样本冷启动时 page-fault 拖慢。
    # 注意: 不再复制到 /dev/shm (那会让每进程另占一份堆内存), mmap 直接共享磁盘 DB 的页缓存。
    do_k2 = 'kraken2' in args.steps.split(',')
    if do_k2 and memory_mapping and args.vmtouch and shutil.which('vmtouch'):
        db_dir = os.path.expanduser(args.kraken2_index)
        k2d_files = sorted(glob.glob(os.path.join(db_dir, '*.k2d')))
        if k2d_files:
            db_size_gb = sum(os.path.getsize(f) for f in k2d_files) / (1024**3)
            logger.info(f"🔥 vmtouch 预热 DB 到 page cache ({len(k2d_files)} 文件 ~{db_size_gb:.0f}GB, 首次约 2-4 分钟)...")
            for f in k2d_files:
                subprocess.run(['vmtouch', '-t', f], check=False,
                              stdin=subprocess.DEVNULL, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            logger.info("✅ 预热完成 (kraken2 --memory-mapping 多进程共享此缓存, 内存仅占一份)")
        else:
            logger.warning(f"vmtouch: 在 {db_dir} 未找到 .k2d 文件")
    elif do_k2 and memory_mapping:
        logger.info("🧠 memory-mapping 模式 (kraken2 mmap 共享 DB; 未预热, 首样本会自动 fault-in)")
    elif do_k2 and args.vmtouch:
        logger.warning("vmtouch 未安装, 跳过预热 (kraken2 --memory-mapping 首样本会自动 fault-in)")

    pipeline = HybridPipeline()
    UI.total_tasks = len(tasks)
    start_time_global = time.time()

    with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as executor:
        futures = {}
        for idx, (ins, outs) in enumerate(tasks, 1):
            seq1, seq2 = ins[0], ins[1] if len(ins) > 1 else None
            out1, out2 = outs[0], outs[1] if len(outs) > 1 else None
            sample_id = extract_sample_name(seq1)
            
            if cm.is_done(sample_id):
                logger.info(f"{UI.GRAY}⏭ [跳过] 任务 {idx:02d} | {sample_id} (历史成功){UI.RESET}")
                UI.update_progress(success=True)
                continue

            future = executor.submit(
                pipeline.run_single,
                task_id=idx, step2_tool=args.tool, seq_type=args.seq_type,
                kr2_idx=args.kraken2_index, step2_idx=args.step2_index,
                seq1=seq1, seq2=seq2, out1=out1, out2=out2,
                threads=args.threads, mfilter=mfilter, config=args.config,
                confidence=args.confidence, tmp_dir=args.tmp, logs_dir=args.logs_dir,
                do_rrna=actual_do_rrna, keep_rrna=args.keep_rrna, chunk_size=args.chunk_size,
                rrna_tool=rrna_tool, silva_index=silva_index,
                keep_taxids=args.keep_taxids, memory_mapping=memory_mapping
            )
            futures[future] = sample_id

        for future in concurrent.futures.as_completed(futures):
            sample_id = futures[future]
            try:
                if future.result() == 0: 
                    cm.mark_done(sample_id) 
                    UI.update_progress(success=True)
                else: 
                    UI.update_progress(success=False)
            except Exception as e:
                logger.error(f"[{sample_id}] 任务发生严重异常: {e}")
                UI.update_progress(success=False)

    # === 后处理：报表解析与绘图 ===

    print("\n")
    logger.info("============== 开始生成数据分析报告 ==============")
    if actual_do_rrna:
        tool_logs_dir = os.path.join(args.logs_dir, "logs")
        rrna_report_path = os.path.join(args.logs_dir, args.rrna_report)
        generate_rrna_report(tool_logs_dir, rrna_report_path)

    summary_tsv_path = os.path.join(args.logs_dir, "host_depletion_seqkit_summary.tsv")
    if os.path.exists(summary_tsv_path):
        if PLOT_AVAILABLE:
            plot_prefix = os.path.join(args.logs_dir, "host_depletion_plot")
            generate_plots(summary_tsv_path, plot_prefix)
        else:
            logger.warning("未检测到 matplotlib 或 pandas，已跳过可视化绘图步骤。")

    # 资源汇总表位置提示
    resource_summary = os.path.join(args.logs_dir, "host_depletion_resource_usage.tsv")
    if os.path.exists(resource_summary):
        logger.info(f"📋 资源使用汇总表已生成: {resource_summary}")
    else:
        logger.info("无样本完成，未生成资源汇总表。")

    end_time_global = time.time()
    logger.info("============== 批量任务执行完毕 ==============")
    logger.info(f"总耗时: {(end_time_global - start_time_global)/60:.2f} 分钟")
    logger.info(f"成功: {UI.success_tasks} | 失败: {UI.fail_tasks} | 总计: {UI.total_tasks}")
    sys.exit(1 if UI.fail_tasks > 0 else 0)

if __name__ == '__main__':
    main()
