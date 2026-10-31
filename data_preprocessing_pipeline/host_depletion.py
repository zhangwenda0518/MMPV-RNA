#!/usr/bin/env python3
"""
终极独立版去宿主 + 去rRNA混合管道 (分类初筛 + 精筛 + Ribodetector)
初筛分类器支持 Kraken2 与 Kun-peng 两种后端 (--classifier):
  · kraken2 : kraken2 --classified-out/--unclassified-out 直接分流 FASTQ
  · kunpeng : kun_peng classify 输出 Kraken 兼容 per-read 结果 + kreport2,
              再用 seqkit grep 按 read ID 反选 (Kun-peng 本身不写分流 FASTQ)
两者输出格式 (per-read 5 列 + kreport2) 与 downstream 完全一致。
【数据可视化版】集成 seqkit stats 与绘图模块，智能识别 rna-short，自适应报表目录结构。
*支持 FASTA/FASTQ 及其 gz 压缩格式自适应处理*
*新增：每样本资源使用统计（时间/内存/CPU）并生成汇总报表*
"""

# ── 图件四格式导出 (SVG / PNG / PDF / PNG@1200dpi) — 默认开启 ──────────
# 2026-10-20: 所有出图脚本统一口径; 用 MMPV_FIG_EXPORT=0 可整体关闭,
# MMPV_FIG_FORMATS=svg,png,pdf,hires 可裁剪, MMPV_FIG_HIRES_DPI 可改高清 dpi。
try:
    from mmpv_common.figure_export import install_savefig_hook as _install_fig_hook
except ImportError:  # 仓库根不在 sys.path (脚本被单独执行) → 向上找 mmpv_common/
    import os as _os
    import sys as _sys
    _d = _os.path.dirname(_os.path.abspath(__file__))
    while _d and not _os.path.isdir(_os.path.join(_d, "mmpv_common")):
        _parent = _os.path.dirname(_d)
        if _parent == _d:
            break
        _d = _parent
    if _d and _d not in _sys.path:
        _sys.path.insert(0, _d)
    from mmpv_common.figure_export import install_savefig_hook as _install_fig_hook
_install_fig_hook()
# ──────────────────────────────────────────────────────────────────────

import os
import re
import sys
import time
import glob
import gzip
import json
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
    # 中文样本名/标签: 微软雅黑 → 黑体 回退, 负号用 ASCII 形 (免 Glyph missing 方框)
    matplotlib.rcParams['font.sans-serif'] = ['Microsoft YaHei', 'SimHei', 'DejaVu Sans']
    matplotlib.rcParams['axes.unicode_minus'] = False
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
            with open(os.path.join(tool_logs_dir, fname), 'r', encoding='utf-8', errors='replace') as f:
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
        with open(output_path, 'w', encoding='utf-8') as f:   # 样本名可为中文: Windows 默认 GBK 会崩
            f.write("Sample\tTotal_sequences\tnon_rRNA\trRNA\tnon_rRNA(%)\trRNA(%)\n")
            for res in sorted(results, key=lambda x: x['sample']):
                f.write(f"{res['sample']}\t{res['total']}\t{res['non_rRNA']}\t{res['rRNA']}\t{res['non_rRNA_percent']:.2f}%\t{res['rRNA_percent']:.2f}%\n")
        logger.info(f"📊 rRNA 去除统计报告已生成: {output_path}")

# --- 绘图函数整合 ---
def plot_data(pivot_df, plot_type, output_file, figsize=(10, 6), dpi=300, max_labels=50, million=False):
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
    """读取 seqkit 汇总表绘制 line/bar PNG; 失败打印 traceback 并返回 False (不静默)"""
    try:
        df = pd.read_csv(tsv_file, sep='\t', header=0)
        required_cols = ['Sample', 'Stage', 'num_seqs']
        if not all(col in df.columns for col in required_cols):
            logger.warning(f"缺少必要的列 {required_cols}，跳过绘图。")
            return False

        df = df[required_cols].copy()
        # num_seqs 可能带千分位逗号/非数值 (seqkit 变体或手工表): 统一转数值, 丢弃不可解析行
        df['num_seqs'] = pd.to_numeric(
            df['num_seqs'].astype(str).str.replace(',', '', regex=False), errors='coerce')
        df = df.dropna(subset=['num_seqs'])
        if df.empty:
            logger.warning("num_seqs 列无可用数值，跳过绘图。")
            return False
        df_grouped = df.groupby(['Sample', 'Stage'], as_index=False)['num_seqs'].mean()
        pivot_df = df_grouped.pivot(index='Sample', columns='Stage', values='num_seqs')
        stages = sorted(pivot_df.columns.tolist())
        pivot_df = pivot_df[stages]

        if pivot_df.empty:
            logger.warning("汇总表无有效样本×阶段数据，跳过绘图。")
            return False

        line_out = output_prefix + '_line.png'
        bar_out = output_prefix + '_bar.png'
        
        plot_data(pivot_df, 'line', line_out)
        plot_data(pivot_df, 'bar', bar_out)
        logger.info(f"📈 可视化图表已生成: {line_out} 和 {bar_out}")
        return True
    except Exception as e:
        logger.error(f"绘图失败: {e}", exc_info=True)
        return False

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
    # 允许直接传绝对路径 (如 --kunpeng-bin /opt/bin/kun_peng), which() 对路径同样有效,
    # 但 Windows 上无扩展名的路径 which() 会补 PATHEXT; 这里额外兜底 isfile。
    missing = [tool for tool in tools
               if shutil.which(tool) is None
               and not (os.path.isfile(tool) and os.access(tool, os.X_OK))]
    if missing:
        logger.error(f"环境依赖缺失！找不到以下命令: {', '.join(missing)}")
        sys.exit(1)

# ── 分类后端 (kraken2 / kunpeng) ──────────────────────────────────────
# Kraken2 库布局: hash.k2d + opts.k2d + taxo.k2d
# Kun-peng 库布局: hash_*.k2d + hash_config.k2d + opts.k2d + taxo.k2d
# `kun_peng hashshard --db <k2db>` 会**原地追加**分片 (保留 hash.k2d), 因此
# 转换后的库同时可被 kraken2 与 kun_peng 使用, 无需复制整库。
KRAKEN2_DB_FILES = ('hash.k2d', 'opts.k2d', 'taxo.k2d')
KUNPENG_DB_FILES = ('hash_config.k2d', 'opts.k2d', 'taxo.k2d')
KUNPENG_SHARD_RE = re.compile(r'^hash_\d+\.k2d$')


def kunpeng_shards(db_dir):
    return sorted(p for p in glob.glob(os.path.join(db_dir, 'hash_*.k2d'))
                  if KUNPENG_SHARD_RE.match(os.path.basename(p)))


def kunpeng_db_ready(db_dir):
    """Kun-peng 库是否可直接 classify (hash 分片 + 元数据齐全)。"""
    return all(os.path.isfile(os.path.join(db_dir, f)) for f in KUNPENG_DB_FILES) \
        and bool(kunpeng_shards(db_dir))


def prepare_kunpeng_db(db_dir, kunpeng_bin='kun_peng', hash_capacity='1G'):
    """确保 --classifier kunpeng 的库可用。

    · 已是 Kun-peng 分片库 → 直接返回;
    · 只有 Kraken2 库 (hash.k2d) → 调 `kun_peng hashshard` 原地转换;
    · 两者都没有 → 返回原路径, 由 verify_indices 报错。
    """
    db_dir = os.path.expanduser(db_dir)
    if kunpeng_db_ready(db_dir):
        return db_dir
    if not os.path.isfile(os.path.join(db_dir, 'hash.k2d')):
        return db_dir
    if os.path.isfile(os.path.join(db_dir, 'hash_config.k2d')):
        # 半转换状态: hashshard 自身会因为 hash_config 已存在而 panic,
        # 明确提示而不是让用户看 Rust panic。
        logger.error(f"Kun-peng 库处于半转换状态 (hash_config.k2d 存在但缺少 hash_*.k2d): {db_dir}")
        logger.error("  请删除 hash_config.k2d 与 hash_*.k2d 后重跑, 或换用干净的 Kraken2 库。")
        sys.exit(1)
    logger.info(f"🔧 Kraken2 库缺少 Kun-peng 分片, 自动执行 hashshard (容量 {hash_capacity})...")
    try:
        extra_gb = os.path.getsize(os.path.join(db_dir, 'hash.k2d')) / 1e9
        logger.info(f"   hashshard 会把 hash 表复制成分片, 预计额外占用 ~{extra_gb:.0f} GB 磁盘")
    except OSError:
        pass
    logger.info(f"   {kunpeng_bin} hashshard --db {db_dir} --hash-capacity {hash_capacity}")
    try:
        subprocess.run([kunpeng_bin, 'hashshard', '--db', db_dir,
                        '--hash-capacity', str(hash_capacity)], check=True)
    except FileNotFoundError:
        logger.error(f"找不到 Kun-peng 可执行文件: {kunpeng_bin}")
        sys.exit(1)
    except subprocess.CalledProcessError as e:
        logger.error(f"hashshard 失败 (exit={e.returncode}); 若目标库已有 hash_config.k2d 请先清理")
        sys.exit(1)
    if not kunpeng_db_ready(db_dir):
        logger.error(f"hashshard 结束但库校验失败: {db_dir}")
        sys.exit(1)
    logger.info(f"✅ Kun-peng 分片就绪: {len(kunpeng_shards(db_dir))} 个分片")
    return db_dir


def verify_indices(kr2_idx, step2_tool, step2_idx, classifier='kraken2'):
    if classifier == 'kunpeng':
        for f in KUNPENG_DB_FILES:
            if not os.path.isfile(os.path.join(kr2_idx, f)):
                logger.error(f"Kun-peng 数据库缺失: {f}")
                sys.exit(1)
        if not kunpeng_shards(kr2_idx):
            logger.error(f"Kun-peng 数据库缺少 hash_*.k2d 分片: {kr2_idx}")
            sys.exit(1)
    else:
        for f in KRAKEN2_DB_FILES:
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
            # Windows 路径用反斜杠, 先归一成 / 再做白名单校验 (源路径不受影响)
            probe = str(arg).replace('\\', '/')
            chars = re.findall(r'[^a-zA-Z0-9\. !#=/:_-]', probe)
            if not re.match(r'^[a-zA-Z0-9\. !#=/:_-]+$|^(?![\s\S])', probe): return False, arg, chars
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
    def kunpeng(index, seq1, seq2=None, threads=1, confidence=0.4,
                output_dir=None, chunk_dir=None, min_hit_groups=2,
                bin_path='kun_peng'):
        """Kun-peng 分类命令 (Kraken2 兼容输出)。

        Kun-peng 与 kraken2 的关键差异:
          · 不写 --classified-out/--unclassified-out, 只有 output_*.txt + *.kreport2;
          · 必须提供干净的 --chunk-dir 与 --output-dir;
          · 双端用 -P 且把 R1/R2 作为两个位置参数传入, 输出是**每个 read 对一行**
            (不是每个 mate 一行), read ID 取 R1 去 /1 后的 base ID。
        """
        # seqkmer 的并行读入断言 n_threads > 2 (0.7.11 实测 -p 2 直接 panic);
        # 这里统一抬到下限, 避免用户 -t 2 时跑到一半崩掉。
        threads = max(3, int(threads))
        cmd = [bin_path, 'classify', '--db', index,
               '--chunk-dir', chunk_dir, '--output-dir', output_dir,
               '-p', str(threads), '-g', str(int(min_hit_groups))]
        if confidence and float(confidence) > 0:
            cmd += ['-T', str(confidence)]
        if seq2:
            cmd.append('-P')
        cmd.append(seq1)
        if seq2:
            cmd.append(seq2)
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
# 4b. Kraken2/Kun-peng read extraction helpers (--keep-taxids)
# ==========================================
def process_kraken_output_line(line):
    """Parse a Kraken2 / Kun-peng --output line, return (taxid_int, read_id_str)."""
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


def report_unclassified_count(report_file):
    """从 kreport(2) 的 U 行取未分类计数 (Kraken2/Kun-peng 通用); 缺失返回 0。"""
    try:
        with open(report_file, 'r', encoding='utf-8', errors='replace') as f:
            for line in f:
                parts = line.rstrip('\n').split('\t')
                if len(parts) >= 5 and parts[3].strip() == 'U' and parts[4].strip() == '0':
                    return int(float(parts[1]))
    except (OSError, ValueError, IndexError):
        pass
    return 0


def build_children_taxids(report_file, target_taxids):
    """Parse a Kraken2/Kun-peng kreport(2), build taxonomy tree, return target_taxids + descendants."""
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
        logger.warning(f"分类 report 中未找到目标 taxID {target_taxids}，将仅匹配 taxID 本身（可能漏掉子节点）")
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


# ==========================================
# 4c. Kun-peng read extraction helpers
# ==========================================
def strip_read_suffix(rid):
    """去掉 read ID 尾部的 /1 /2 (Kun-peng 输出的是 R1 base ID)。"""
    if rid.endswith('/1') or rid.endswith('/2'):
        return rid[:-2]
    return rid


def _open_text_auto(path):
    if str(path).lower().endswith('.gz'):
        return gzip.open(path, 'rt', encoding='utf-8', errors='replace')
    return open(path, 'r', encoding='utf-8', errors='replace')


def detect_read_suffix(path, sample_n=100):
    """探测 FASTQ/FASTA 首 token 的配对后缀 ('' / '/1' / '/2')。

    MGI/华大: @LH00278:...:26992/1 (后缀在首 token) → 需要把 /1 /2 写进
    seqkit 模式文件才能精确命中; Illumina/SRA 的后缀在第二 token,
    首 token 无后缀 → 返回 ''。
    """
    counts = {}
    first = ''
    try:
        with _open_text_auto(path) as fh:
            line = fh.readline()
            if not line:
                return ''
            step = 2 if line.startswith('>') else 4
            idx = heads = 0
            while line:
                s = line.rstrip('\r\n')
                if s:
                    if idx % step == 0:
                        tok = s.split(None, 1)[0]
                        suf = tok[-2:] if (tok.endswith('/1') or tok.endswith('/2')) else ''
                        if heads == 0:
                            first = suf
                        counts[suf] = counts.get(suf, 0) + 1
                        heads += 1
                        if heads >= sample_n:
                            break
                    idx += 1
                line = fh.readline()
    except OSError as e:
        logger.warning(f"配对后缀探测失败（按无后缀处理）: {os.path.basename(str(path))}: {e}")
        return ''
    if len(counts) > 1:
        dist = '，'.join(f'{k or "(无)"} ×{v}' for k, v in counts.items())
        logger.warning(f"配对后缀不一致 ({os.path.basename(str(path))}: {dist})，按首条取 {first or '(无)'}")
    return first


def write_suffixed_patterns(src_pattern, dst_pattern, suffix):
    """把 base-ID 模式文件复制一份并追加后缀 (MGI header 需要)。"""
    if not suffix:
        return src_pattern
    with open(src_pattern, 'r', encoding='utf-8') as fin, \
         open(dst_pattern, 'w', encoding='utf-8') as fout:
        for line in fin:
            rid = line.rstrip('\n')
            if rid:
                fout.write(rid + suffix + '\n')
    return dst_pattern


def seqkit_filter_by_ids(seqkit, src, dst, pattern_file, invert, threads=4):
    """seqkit grep [-v] -f <patterns>：按 read ID 精确分流 (hash 查表)。"""
    cmd = [seqkit, 'grep']
    if invert:
        cmd.append('-v')
    cmd += ['-f', str(pattern_file), '-j', str(threads), '-o', str(dst), str(src)]
    subprocess.run(cmd, check=True)


def kractor_filter_reads(kractor_bin, seq1, seq2, out1, out2, raw_txt, raw_report,
                         keep_taxids_list, mfilter, prefix=''):
    """用 kractor 直接按 taxid 分流，返回 total_output_records。

    mfilter=True  → 保留 unclassified(0) + keep_taxids(含子代)
    mfilter=False → --exclude 上述集合，只留宿主
    kractor 依赖完整的 per-read 输出；缺行时调用方应先回退 seqkit。
    """
    for p in (out1, out2):
        if p and os.path.exists(p):
            os.remove(p)
    cmd = [kractor_bin, '-i', str(seq1)]
    if seq2:
        cmd += ['-i', str(seq2)]
    cmd += ['-o', str(out1)]
    if seq2:
        cmd += ['-o', str(out2)]
    cmd += ['-k', str(raw_txt)]
    if keep_taxids_list:
        if not raw_report or not os.path.isfile(raw_report):
            raise RuntimeError("kractor --children 需要 kreport2，但文件不存在")
        cmd += ['-r', str(raw_report), '--children']
    if not mfilter:
        cmd.append('--exclude')
    cmd += ['-t', '0'] + [str(t) for t in keep_taxids_list]
    cmd.append('--summary')
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"kractor 失败 (exit={proc.returncode}): "
                           f"{(proc.stderr or '')[-400:]}")
    out = proc.stdout or ''
    lo, hi = out.find('{'), out.rfind('}')
    if lo < 0 or hi < lo:
        raise RuntimeError("kractor --summary 未返回 JSON")
    try:
        summary = json.loads(out[lo:hi + 1])
        return int(summary['total_output_records'])
    except (ValueError, KeyError) as e:
        raise RuntimeError(f"kractor summary 解析失败: {e}")


def parse_kunpeng_output(out_txt, report_file, host_pattern,
                         keep_taxids_list, prefix, sample_name):
    """解析 Kun-peng per-read 输出 → 宿主 read ID 模式文件 + 统计。

    out_txt: `output_*.txt` 路径列表 (单样本 1 个; batch 传入该样本对应的那个)。
    """
    expanded = build_children_taxids(report_file, keep_taxids_list) if keep_taxids_list else set()
    n_classified = n_host = n_rescued = n_u_lines = n_output_lines = 0
    with open(host_pattern, 'w', encoding='utf-8') as pf:
        for txt in out_txt:
            with open(txt, 'r', encoding='utf-8', errors='replace') as fh:
                for line in fh:
                    tid, rid = process_kraken_output_line(line)
                    if not rid:
                        continue
                    n_output_lines += 1
                    if tid <= 0:
                        n_u_lines += 1    # U 行: 非宿主 (部分版本会写, 部分不写)
                        continue
                    n_classified += 1
                    if expanded and tid in expanded:
                        n_rescued += 1    # 病毒等保留类群: 不算宿主
                        continue
                    pf.write(strip_read_suffix(rid) + '\n')
                    n_host += 1
    # 总数以 kreport 的 U 行为准: output 里的 U 行与它重叠, 不能相加
    # (0.7.11 实测 50k 对: output 47,075 C + 1,861 U, kreport U=2,925;
    #  真实总数 = 47,075 + 2,925 = 50,000)。
    n_unclassified = report_unclassified_count(report_file)
    if n_unclassified == 0 and n_u_lines:
        n_unclassified = n_u_lines      # 无 kreport 时的兜底
    total = n_classified + n_unclassified
    stats = {'records': total, 'classified': n_classified,
             'unclassified': n_unclassified, 'output_u_lines': n_u_lines,
             'output_lines': n_output_lines,
             'total': total, 'host': n_host, 'rescued': n_rescued,
             'outputs': len(out_txt)}
    logger.info(f"{prefix}   Kun-peng 分类 {n_classified:,} 条 + 未分类 {n_unclassified:,} 条"
                f" (合计 {total:,}); 宿主 {n_host:,} 条"
                f" ({(n_host / max(total, 1)) * 100:.2f}%); 救回保留 taxID {n_rescued:,} 条")
    if n_host > 5_000_000:
        # seqkit grep -f 会把全部模式装进内存哈希; 大样本并行时这是主要内存项
        logger.warning(f"{prefix}   宿主 ID {n_host:,} 个; seqkit grep 约需 "
                       f"{n_host * 80 / 1e9:.1f} GB 内存做模式匹配, 注意控制 --jobs")
    return stats


def kunpeng_split_reads(seq1, seq2, host_pattern, n_host, mfilter, temp_dir,
                        ext_tmp, prefix, sample_name, threads,
                        splitter='seqkit', kractor_bin='kractor',
                        raw_txt=None, raw_report=None, keep_taxids_list=None,
                        stats=None):
    """按 taxid/宿主 ID 分流 (mfilter=True 保留非宿主; False 只留宿主)。

    splitter='kractor' / 'auto': 直接消费 Kun-peng 的 per-read 输出 + kreport2，
    无需先构建宿主 ID 集合；输出条数必须与预期一致，否则回退 seqkit。
    """
    expected = None
    if stats:
        expected = (stats.get('total', 0) - n_host) if mfilter else n_host
    complete = (bool(stats)
                and stats.get('output_lines', -1) == stats.get('total', -2))
    if splitter in ('auto', 'kractor') and raw_txt and os.path.isfile(raw_txt):
        if not complete:
            logger.warning(f"{prefix}   Kun-peng 输出不完整 "
                           f"({stats.get('output_lines')}/{stats.get('total')})，"
                           f"跳过 kractor，回退 seqkit")
        else:
            dst1 = os.path.join(temp_dir, f'kp_filtered_1{ext_tmp}')
            dst2 = os.path.join(temp_dir, f'kp_filtered_2{ext_tmp}') if seq2 else None
            try:
                n = kractor_filter_reads(
                    kractor_bin, seq1, seq2, dst1, dst2, raw_txt, raw_report,
                    keep_taxids_list or [], mfilter, prefix)
                expected_records = (expected * (2 if seq2 else 1)
                                    if expected is not None else None)
                if expected_records is not None and n != expected_records:
                    raise RuntimeError(f"输出 {n} 条记录，预期 {expected_records} 条"
                                       f"（{expected} 对）")
                logger.info(f"{prefix}   kractor 分流完成: "
                            f"{(n // 2) if seq2 else n:,} 条"
                            f"{' (paired)' if seq2 else ''}")
                return dst1, dst2
            except Exception as e:
                logger.warning(f"{prefix}   kractor 不可用/计数不符，回退 seqkit: {e}")
                for p in (dst1, dst2):
                    if p and os.path.exists(p):
                        try:
                            os.remove(p)
                        except OSError:
                            pass

    seqkit = shutil.which('seqkit') or 'seqkit'
    inputs = [seq1] + ([seq2] if seq2 else [])
    outs = []
    for i, src in enumerate(inputs, 1):
        dst = os.path.join(temp_dir, f'kp_filtered_{i}{ext_tmp}')
        if n_host == 0:
            # 无宿主命中: invert 模式原样保留输入 (比空跑一遍 seqkit 省一大截 I/O)。
            # 注意 gz 输入必须解压后再写, 否则 dst 是 .fastq 却装着 gzip 字节流。
            if mfilter:
                if str(src).lower().endswith('.gz'):
                    with gzip.open(src, 'rb') as fin, open(dst, 'wb') as fout:
                        shutil.copyfileobj(fin, fout)
                else:
                    shutil.copyfile(src, dst)
            else:
                open(dst, 'w').close()
        else:
            suffix = detect_read_suffix(src)
            pat = write_suffixed_patterns(
                host_pattern, os.path.join(temp_dir, f'{sample_name}_kp_host_ids_{i}.txt'), suffix)
            try:
                seqkit_filter_by_ids(seqkit, src, dst, pat, invert=mfilter, threads=threads)
            except subprocess.CalledProcessError as e:
                logger.error(f"{prefix}   seqkit grep 失败 (exit={e.returncode}): {src}")
                return None
            except FileNotFoundError:
                logger.error(f"{prefix}   找不到 seqkit; --classifier kunpeng 需要 seqkit 做 read 分流")
                return None
        outs.append(dst)
    return outs[0], (outs[1] if len(outs) > 1 else None)


def kunpeng_classify_and_split(kr2_idx, seq1, seq2, threads, confidence,
                               keep_taxids_list, mfilter, temp_dir,
                               report_file, ext_tmp, prefix, sample_name,
                               tool_logs_dir, kunpeng_bin='kun_peng',
                               min_hit_groups=2, keep_classifier_output=False,
                               keep_dir=None, splitter='seqkit',
                               kractor_bin='kractor'):
    """单样本 Kun-peng 初筛 + 分流 (--kunpeng-batch-size 1 时使用)。

    批量模式请用 kunpeng_batch_classify(): 一次进程处理整批, 分片只加载一次。
    """
    keep_classifier_output = bool(keep_classifier_output and keep_dir)
    raw_base = keep_dir if keep_classifier_output else temp_dir
    kp_out = os.path.join(raw_base, 'kunpeng_out')
    kp_chunk = os.path.join(raw_base, 'kunpeng_chunk')
    # --chunk-dir 必须是干净目录 (上游硬性要求)
    for d in (kp_out, kp_chunk):
        if os.path.isdir(d):
            shutil.rmtree(d, ignore_errors=True)
        os.makedirs(d, exist_ok=True)

    kr2_cmd = Aligners.kunpeng(kr2_idx, seq1, seq2=seq2, threads=threads,
                               confidence=confidence, output_dir=kp_out,
                               chunk_dir=kp_chunk, min_hit_groups=min_hit_groups,
                               bin_path=kunpeng_bin)
    ret, mem, cpu_time = execute([kr2_cmd], pipe=False, log_prefix=prefix,
                                 logs_dir=tool_logs_dir, sample_id=sample_name,
                                 step_name="step1_kunpeng")
    stats = {'mem': mem, 'cpu_time': cpu_time, 'records': 0,
             'host': 0, 'rescued': 0, 'outputs': 0}
    if any(c != 0 for c in ret):
        return None, None, stats

    out_txt = sorted(glob.glob(os.path.join(kp_out, 'output_*.txt')))
    # 优先取 output_1.kreport2; 多文件时 resolve 还会写一个 output_min-max 合并
    # 报告, 字典序会排在前面, 但它不适合当单样本报告。
    preferred_report = os.path.join(kp_out, 'output_1.kreport2')
    reports = ([preferred_report] if os.path.isfile(preferred_report)
               else sorted(glob.glob(os.path.join(kp_out, '*.kreport2'))))
    stats['outputs'] = len(out_txt)
    if reports:
        # 与 kraken2 --report 同路径, 供 build_children_taxids / 报告复用
        shutil.copyfile(reports[0], report_file)
    elif keep_taxids_list:
        logger.warning(f"{prefix}   Kun-peng 未生成 kreport2, --keep-taxids 只能匹配 taxID 本身")

    host_pattern = os.path.join(temp_dir, f'{sample_name}_kp_host_ids.txt')
    stats.update(parse_kunpeng_output(out_txt, report_file, host_pattern,
                                      keep_taxids_list, prefix, sample_name))
    temp1, temp2 = kunpeng_split_reads(seq1, seq2, host_pattern, stats['host'],
                                       mfilter, temp_dir, ext_tmp, prefix,
                                       sample_name, threads,
                                       splitter=splitter, kractor_bin=kractor_bin,
                                       raw_txt=out_txt[0] if out_txt else None,
                                       raw_report=report_file,
                                       keep_taxids_list=keep_taxids_list,
                                       stats=stats)

    if not keep_classifier_output:
        shutil.rmtree(kp_out, ignore_errors=True)
        shutil.rmtree(kp_chunk, ignore_errors=True)
    else:
        logger.info(f"{prefix}   已保留 Kun-peng 中间输出: {kp_out}")
    if temp1 is None:
        return None, None, stats
    return temp1, temp2, stats


def kunpeng_batch_classify(batch, db, kunpeng_bin, threads, confidence,
                           keep_taxids_list, logs_dir, work_root,
                           mode='classify', min_hit_groups=2,
                           keep_output=False, splitter='seqkit',
                           kractor_bin='kractor', stage_shm=True):
    """批量 Kun-peng 分类: 一个进程处理整批, 索引/分片只加载一次。

    batch: [(task_idx, [seq1, seq2?], [out1, out2?], sample_id), ...]
    返回 ({task_idx: {'host_pattern': path, 'stats': stats}}, batch_dir)。

    mode='classify' (Kun-peng-M, 分片按需加载, 低内存) 支持 .list 输入;
    mode='direct'   (Kun-peng-F, 全量加载, 最快) 不支持 .list, 直接传所有路径。
    PE/SE 会分成两次调用 (同一个 -P 批次不能混单端)。
    """
    os.makedirs(work_root, exist_ok=True)
    batch_dir = tempfile.mkdtemp(prefix='kunpeng_batch_', dir=work_root)
    batch_id = os.path.basename(batch_dir)
    raw_dir = os.path.join(logs_dir, 'kunpeng_raw', batch_id) if keep_output else batch_dir
    report_dir = os.path.join(logs_dir, 'kraken2_report')
    os.makedirs(report_dir, exist_ok=True)
    out_root = os.path.join(raw_dir, 'out')
    chunk_root = os.path.join(raw_dir, 'chunk')
    os.makedirs(out_root, exist_ok=True)
    os.makedirs(chunk_root, exist_ok=True)

    result = {}
    keep_raw = keep_output or splitter in ('auto', 'kractor')

    # --- 输入预载 /dev/shm: 阵列 ~30MB/s 喂不饱批分类线程; 先把本批输入整批拷进
    # --- 内存盘 (有 cpx 用 cpx), splitr/direct/kractor 全改读内存盘, 批结束清理。
    # --- 空间不足或校验失败则跳过, 退回阵列直读 (原行为)。
    stage_dir = None
    path_map = {}
    if stage_shm:
        try:
            # 预载预算 = shm 剩余 × 0.85 - 5GB (给 kraken2 挂载/其他用户留余量);
            # 逐样本自适应: 按批内顺序尽量多装, 单样本装不下就跳过试下一个 (部分预载)
            _budget = int(shutil.disk_usage('/dev/shm').free * 0.85) - (5 << 30)
            _picked, _used = [], 0
            for idx, ins, _o, _s in batch:
                _paths = [os.path.abspath(p) for p in ins]
                _sz = sum(os.path.getsize(p) for p in _paths)
                if _sz > 0 and _used + int(_sz * 1.05) <= _budget:
                    _picked.append((idx, _paths))
                    _used += int(_sz * 1.05)
            all_files = sorted({p for _i, ps in _picked for p in ps})
            in_total = sum(os.path.getsize(p) for p in all_files)
            # 解压直读判定: kunpeng 的 gz 解压是单线程 (~40MB/s/文件), 100 线程喂不饱;
            # 预载时并行解压成明文 FASTA, 分类读明文 (扫描是并行的, 能吃满线程)。
            # 明文体积 ≈ gz × 4.2, 装不下就退回 gz 预载。
            _gz_files = [p for p in all_files if p.lower().endswith('.gz')]
            _decompress = bool(_gz_files) and shutil.disk_usage('/dev/shm').free > int(in_total * 4.2) + (10 << 30)
            if _decompress:
                logger.info(f"   解压直读: {len(_gz_files)}/{len(all_files)} 个 gz 输入将并行解压为明文 FASTA"
                            f" (绕过单线程解压瓶颈)")
            if in_total > 0:
                stage_dir = os.path.join('/dev/shm', f'kunpeng_stage_{os.getpid()}_{batch_id}')
                os.makedirs(stage_dir, exist_ok=True)
                if len(_picked) == len(batch):
                    logger.info(f"   📥 输入预载 /dev/shm: 整批 {len(_picked)} 样本 "
                                f"{in_total/(1024**3):.1f}GB → {stage_dir}")
                else:
                    logger.info(f"   📥 输入预载 /dev/shm: 自适应 {len(_picked)}/{len(batch)} 样本 "
                                f"{in_total/(1024**3):.1f}GB (预算 {_budget/(1024**3):.0f}GB), 其余阵列直读")
                import threading as _th
                _stop = _th.Event()

                def _stage_progress():
                    while not _stop.wait(30):
                        try:
                            _d = sum(os.path.getsize(os.path.join(r, x))
                                     for r, _dirs, fs in os.walk(stage_dir) for x in fs)
                            logger.info(f"   ⏳ 预载进度: {_d/(1024**3):.1f}/{in_total/(1024**3):.1f}GB "
                                        f"({100.0*_d/max(1,in_total):.1f}%)")
                        except Exception:
                            pass
                _pt = _th.Thread(target=_stage_progress, daemon=True)
                _pt.start()
                cpx_bin = shutil.which('cpx')
                _ok = False
                if _decompress:
                    import concurrent.futures as _cf
                    _gz_src = [p for p in all_files if p.lower().endswith('.gz')]
                    _plain_src = [p for p in all_files if not p.lower().endswith('.gz')]
                    _pigz = shutil.which('pigz') or shutil.which('gzip')
                    logger.info(f"   复制工具: {os.path.basename(_pigz)} -dc × 8 并行 (解压直读)")
                    def _dgz(s):
                        d = os.path.join(stage_dir,
                                         os.path.basename(s)[:-3] if s.lower().endswith('.gz')
                                         else os.path.basename(s))
                        try:
                            with open(d, 'wb') as _out:
                                _r = subprocess.run([_pigz, '-dc', s], stdout=_out,
                                                    stderr=subprocess.DEVNULL,
                                                    stdin=subprocess.DEVNULL)
                            return _r.returncode == 0 and os.path.getsize(d) > 0
                        except Exception:
                            return False
                    with _cf.ThreadPoolExecutor(max_workers=8) as _ex:
                        _res = list(_ex.map(_dgz, _gz_src))
                    for p in _plain_src:
                        _res.append(_dgz(p))
                    _ok = all(_res)
                elif cpx_bin:
                    logger.info("   复制工具: cpx (并行)")
                    _r = subprocess.run([cpx_bin, '-j', '8', '-t', stage_dir, *all_files],
                                        check=False, stdin=subprocess.DEVNULL,
                                        stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
                    _ok = (_r.returncode == 0)
                else:
                    logger.info("   复制工具: shutil.copyfile (未安装 cpx)")
                if not _ok:
                    for p in all_files:
                        shutil.copyfile(p, os.path.join(stage_dir, os.path.basename(p)))
                _stop.set(); _pt.join(timeout=2)
                if _decompress:
                    _ok = all(os.path.isfile(os.path.join(stage_dir, os.path.basename(p)[:-3] if p.lower().endswith('.gz') else os.path.basename(p))) and
                              os.path.getsize(os.path.join(stage_dir, os.path.basename(p)[:-3] if p.lower().endswith('.gz') else os.path.basename(p))) > 0
                              for p in all_files)
                else:
                    _ok = all(os.path.isfile(os.path.join(stage_dir, os.path.basename(p))) and
                              os.path.getsize(os.path.join(stage_dir, os.path.basename(p))) == os.path.getsize(p)
                              for p in all_files)
                if _ok:
                    for idx, _paths in _picked:
                        path_map[idx] = [os.path.join(stage_dir,
                                                      os.path.basename(p)[:-3] if p.lower().endswith('.gz')
                                                      else os.path.basename(p))
                                         for p in _paths]
                else:
                    logger.warning("   预载校验失败, 本批退回阵列直读")
                    shutil.rmtree(stage_dir, ignore_errors=True)
                    stage_dir = None
            else:
                logger.info("   /dev/shm 空间不足, 本批输入直接从阵列读取")
        except OSError as e:
            logger.warning(f"   输入预载跳过 ({e}), 退回阵列直读")
            stage_dir = None

    groups = [(True, [x for x in batch if len(x[1]) > 1]),
              (False, [x for x in batch if len(x[1]) == 1])]
    for paired, items in groups:
        if not items:
            continue
        tag = 'pe' if paired else 'se'
        out_dir = os.path.join(out_root, tag)
        chunk_dir = os.path.join(chunk_root, tag)
        os.makedirs(out_dir, exist_ok=True)
        os.makedirs(chunk_dir, exist_ok=True)

        cmd = [kunpeng_bin, mode, '--db', db, '--output-dir', out_dir,
               '-p', str(max(3, int(threads))), '-g', str(int(min_hit_groups))]
        if confidence and float(confidence) > 0:
            cmd += ['-T', str(confidence)]
        if paired:
            cmd.append('-P')
        if mode == 'direct':
            # direct 没有 .list 支持, 直接展开输入路径
            for _idx, ins, _outs, _sid in items:
                cmd.extend(os.path.abspath(p) for p in (path_map.get(_idx) or ins))
        else:
            cmd += ['--chunk-dir', chunk_dir]
            # Kun-peng 只在扩展名为 .txt 时把单文件当作输入清单展开
            list_file = os.path.join(batch_dir, f'{batch_id}_{tag}.txt')
            with open(list_file, 'w', encoding='utf-8') as lf:
                for _idx, ins, _outs, _sid in items:
                    lf.write('\n'.join(os.path.abspath(p) for p in (path_map.get(_idx) or ins)) + '\n')
            cmd.append(list_file)

        logger.info(f"   [batch {batch_id}] {mode} {'PE' if paired else 'SE'} "
                    f"× {len(items)} 样本 (索引只加载一次)...")
        ret, mem, cpu_time = execute([cmd], pipe=False,
                                     log_prefix=f"[batch {batch_id}]",
                                     logs_dir=os.path.join(logs_dir, 'logs'),
                                     sample_id=f'batch_{batch_id}_{tag}',
                                     step_name=f'step1_kunpeng_{mode}_{tag}')
        if any(c != 0 for c in ret):
            raise RuntimeError(f"Kun-peng 批量分类失败 (mode={mode}, {tag})")
        logger.info(f"   [batch {batch_id}] {tag} 完成: 峰值内存 {mem:.1f}MB | "
                    f"CPU {cpu_time:.1f}s | {len(items)} 样本共享一次索引加载")

        for j, (idx, ins, _outs, sample_id) in enumerate(items, 1):
            txt = os.path.join(out_dir, f'output_{j}.txt')
            kre = os.path.join(out_dir, f'output_{j}.kreport2')
            report_file = os.path.join(report_dir, f'{sample_id}_kraken2_report.txt')
            if os.path.isfile(kre):
                shutil.copyfile(kre, report_file)
            elif keep_taxids_list:
                logger.warning(f"[batch {batch_id}] {sample_id}: 未生成 kreport2")
            host_pattern = os.path.join(batch_dir, f'host_ids_{idx}.txt')
            prefix = f"[batch {batch_id} | {sample_id}]"
            if os.path.isfile(txt):
                stats = parse_kunpeng_output([txt], report_file, host_pattern,
                                             keep_taxids_list, prefix, sample_id)
            else:
                # 空输入/全部被质量过滤时上游可能不写输出; 按 0 reads 处理
                open(host_pattern, 'w').close()
                stats = {'records': 0, 'classified': 0, 'unclassified': 0,
                         'output_u_lines': 0, 'output_lines': 0, 'total': 0, 'host': 0,
                         'rescued': 0, 'outputs': 0}
                logger.warning(f"{prefix} 无 per-read 输出 (按 0 reads 处理)")
            stats['mem'] = 0.0        # 批量资源记在 batch 日志, 不重复计入每样本
            stats['cpu_time'] = 0.0
            result[idx] = {'host_pattern': host_pattern, 'stats': stats,
                           'raw_txt': txt if keep_raw else None,
                           'raw_report': kre if keep_raw else None,
                           'staged_ins': path_map.get(idx)}
            if not keep_output and not keep_raw:
                for f in (txt, kre):
                    try:
                        os.remove(f)
                    except OSError:
                        pass

    if not keep_output:
        # keep_raw: out_root 留给 run_single 用 kractor 后自行删除；
        # batch_dir 由调用方在整批下游结束后清理。
        if not keep_raw:
            shutil.rmtree(out_root, ignore_errors=True)
        shutil.rmtree(chunk_root, ignore_errors=True)
    else:
        logger.info(f"   已保留批量 Kun-peng 中间输出: {raw_dir}")
    return result, batch_dir, stage_dir

class HybridPipeline:
    def run_single(self, task_id, step2_tool, seq_type, kr2_idx, step2_idx, seq1, out1, seq2=None, out2=None,
                   threads=1, mfilter=True, config='', confidence=0.4, tmp_dir=None, logs_dir=None,
                   do_rrna=False, keep_rrna=False, chunk_size=256,
                   rrna_tool='ribodetector', silva_index=None,
                   keep_taxids=None, memory_mapping=True,
                   classifier='kraken2', kunpeng_bin='kun_peng',
                   min_hit_groups=2, keep_classifier_output=False,
                   kunpeng_precomputed=None,
                   rrna_impl='subprocess', rrna_threads=None,
                   splitter='auto', kractor_bin='kractor'):
        
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

            # --- [Stage 2] 分类初筛 (kraken2 / kunpeng) ---
            if do_kraken2:
                keep_taxids_list = [int(t) for t in str(keep_taxids).split(',') if t.strip()] if keep_taxids else []
                do_keep_taxids = bool(keep_taxids_list)

                if classifier == 'kunpeng':
                    if kunpeng_precomputed:
                        # 批量模式: 分类已在主线程一次完成 (分片只加载一次),
                        # 这里只做本样本的 taxid 分流 (kractor 或 seqkit)。
                        kp_stats = kunpeng_precomputed['stats']
                        logger.info(f"{prefix} -> [Stage 2] 使用批量 Kun-peng 结果 "
                                    f"(宿主 {kp_stats.get('host', 0):,} / "
                                    f"{kp_stats.get('total', 0):,})...")
                        temp1, temp2 = kunpeng_split_reads(
                            seq1, seq2, kunpeng_precomputed['host_pattern'],
                            kp_stats.get('host', 0), mfilter, temp_dir, ext_tmp,
                            prefix, sample_name, threads,
                            splitter=splitter, kractor_bin=kractor_bin,
                            raw_txt=kunpeng_precomputed.get('raw_txt'),
                            raw_report=kunpeng_precomputed.get('raw_report'),
                            keep_taxids_list=keep_taxids_list, stats=kp_stats)
                        if not keep_classifier_output:
                            for p in (kunpeng_precomputed.get('raw_txt'),
                                      kunpeng_precomputed.get('raw_report')):
                                if p and os.path.exists(p):
                                    try:
                                        os.remove(p)
                                    except OSError:
                                        pass
                        if temp1 is None:
                            return 1
                    else:
                        logger.info(f"{prefix} -> [Stage 2] Kun-peng 初筛 (Conf: {confidence}, "
                                    f"保留 taxIDs: {keep_taxids_list if do_keep_taxids else '无'}, "
                                    f"min-hit-groups={min_hit_groups}, splitter={splitter})...")
                        temp1, temp2, kp_stats = kunpeng_classify_and_split(
                            kr2_idx, seq1, seq2, threads, confidence, keep_taxids_list,
                            mfilter, temp_dir, report_file, ext_tmp, prefix, sample_name,
                            tool_logs_dir, kunpeng_bin=kunpeng_bin,
                            min_hit_groups=min_hit_groups,
                            keep_classifier_output=keep_classifier_output,
                            keep_dir=os.path.join(logs_dir, 'kunpeng_raw', sample_name),
                            splitter=splitter, kractor_bin=kractor_bin)
                    max_mem_used = max(max_mem_used, kp_stats.get('mem', 0.0))
                    total_cpu_time += kp_stats.get('cpu_time', 0.0)
                    if temp1 is None:
                        return 1
                    if kp_stats['records'] == 0:
                        logger.warning(f"{prefix}   Kun-peng 无分类记录 (输出为空)")
                elif do_keep_taxids:
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
                    # 输出未压缩, 统一交给 [Stage 8] gzip 压缩 (避免 gz(gz(fasta)) 双重压缩)
                    logger.info(f"{prefix} -> [Stage 7] SILVA Bowtie2 剔除 rRNA...")

                    # --un-conc 会在文件名后加 .1/.2 后缀, 需先写临时前缀再 rename
                    silva_tmp_prefix = os.path.join(temp_dir, f"{sample_name}_silva_unmapped")
                    silva_cmd = ["bowtie2", "--very-sensitive-local", "-p", str(threads),
                                 "-x", silva_index]
                    if seq2:
                        silva_cmd += ["-1", extract_tgt1, "-2", extract_tgt2,
                                      "--un-conc", silva_tmp_prefix]
                    else:
                        silva_cmd += ["-U", extract_tgt1,
                                      "--un", actual_out1]
                    ret, mem, cpu_time = execute([silva_cmd], pipe=True, log_prefix=prefix,
                                                  logs_dir=tool_logs_dir, sample_id=sample_name,
                                                  step_name="step5_silva_bowtie2")
                    max_mem_used = max(max_mem_used, mem)
                    total_cpu_time += cpu_time
                    if any(c != 0 for c in ret): return 1

                    # bowtie2 --un-conc 产出 silva_tmp_prefix.1 / silva_tmp_prefix.2
                    # rename 到 actual_out1 / actual_out2
                    if seq2:
                        os.rename(silva_tmp_prefix + ".1", actual_out1)
                        os.rename(silva_tmp_prefix + ".2", actual_out2)

                    # keep_rrna: 重新比对并保留 rRNA reads (gzip 内容 + .fq.gz 名, 见 doc/03-host_depletion.md)
                    if keep_rrna:
                        rrna_dir = os.path.join(outdir_base, "rrna")
                        os.makedirs(rrna_dir, exist_ok=True)
                        rrna_out1 = os.path.join(rrna_dir, f"{sample_name}.silva_rrna.1.fq.gz") if seq2 else os.path.join(rrna_dir, f"{sample_name}.silva_rrna.fq.gz")
                        rrna_out2 = os.path.join(rrna_dir, f"{sample_name}.silva_rrna.2.fq.gz") if seq2 else None
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
                        elif not os.path.exists(rrna_out1):
                            logger.warning(f"{prefix} ⚠️ SILVA rRNA 输出未生成: {rrna_out1} "
                                           f"(bowtie2 --al-gz 可能未比对到任何 rRNA reads)")

                else:
                    # ── Ribodetector 模式 (原有逻辑) ──
                    logger.info(f"{prefix} -> [Stage 7] Ribodetector 剔除 rRNA 序列...")
                    avg_len = get_avg_length(extract_tgt1, threads)
                    length_threshold = int(avg_len * 0.5) if avg_len > 0 else 50
                    logger.debug(f"{prefix} 推测平均读长 {avg_len}，设定过滤阈值 {length_threshold}")

                    if keep_rrna:
                        rrna_dir = os.path.join(outdir_base, "rrna")
                        os.makedirs(rrna_dir, exist_ok=True)
                        # 命名对齐 doc/03-host_depletion.md: {sample}.rrna.[1|2].fq.gz (Stage 6 压缩后 .fq → .fq.gz)
                        rrna_out1 = os.path.join(rrna_dir, f"{sample_name}.rrna.1.fq") if seq2 else os.path.join(rrna_dir, f"{sample_name}.rrna.fq")
                        rrna_out2 = os.path.join(rrna_dir, f"{sample_name}.rrna.2.fq") if seq2 else None
                    else:
                        rrna_out1 = os.path.join(temp_dir, f"{sample_name}.rrna.1{ext_tmp}") if seq2 else os.path.join(temp_dir, f"{sample_name}.rrna{ext_tmp}")
                        rrna_out2 = os.path.join(temp_dir, f"{sample_name}.rrna.2{ext_tmp}") if seq2 else None

                    rrna_log = os.path.join(tool_logs_dir, f"{sample_name}_ribodetector.log")
                    rrna_done = False
                    if rrna_impl == 'inprocess':
                        try:
                            from ribodetector_runner import get_runner
                            runner = get_runner(
                                seq_len=length_threshold,
                                threads=(rrna_threads or threads),
                                batch_size=1024,
                                chunk_size=chunk_size,
                                ensure='rrna')
                            t_rrna = time.time()
                            stats = runner.classify(
                                extract_tgt1, extract_tgt2,
                                actual_out1, actual_out2,
                                rrna_out1 if keep_rrna else None,
                                rrna_out2 if keep_rrna else None,
                                log_path=rrna_log)
                            logger.info(f"{prefix}   in-process Ribodetector: "
                                        f"total={stats['total']:,}, "
                                        f"non-rRNA={stats['non_rRNA']:,}, "
                                        f"rRNA={stats['rRNA']:,} "
                                        f"({time.time() - t_rrna:.1f}s)")
                            rrna_done = True
                        except Exception as e:
                            logger.warning(f"{prefix} ⚠️ in-process Ribodetector 失败 "
                                           f"({type(e).__name__}: {e})，回退 subprocess")

                    if not rrna_done:
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
    parser.add_argument('-k', '--kraken2_index', required=True, help="宿主分类数据库 (Kraken2 库或 Kun-peng 分片库)")
    parser.add_argument('-x', '--step2_index', help="第二步工具索引前缀 (--steps 含 align 时必需)")

    # 初筛分类后端: kraken2 / kunpeng。两者共用同一宿主库:
    # `kun_peng hashshard` 会把 Kraken2 库原地转成 Kun-peng 分片库 (保留 hash.k2d),
    # 因此同一个 --kraken2_index 可在两种后端间切换。
    parser.add_argument('--classifier', choices=['kraken2', 'kunpeng'], default='kraken2',
                        help="初筛分类器 (默认: kraken2; kunpeng 需 kun_peng 可执行文件, "
                             "会自动用 hashshard 把 Kraken2 库转成分片格式)")
    parser.add_argument('--kunpeng-bin', default=os.environ.get('MMPV_KUNPENG_BIN', 'kun_peng'),
                        help="Kun-peng 可执行文件 (默认: $MMPV_KUNPENG_BIN 或 PATH 中的 kun_peng)")
    parser.add_argument('--kunpeng-hash-capacity', default='1G',
                        help="自动 hashshard 转换时的分片容量 (默认: 1G ≈ 4GiB/片)")
    parser.add_argument('--min-hit-groups', type=int, default=2,
                        help="Kun-peng 最小 hit group 数 (默认: 2; kraken2 无此参数, 忽略)")
    parser.add_argument('--kunpeng-mode', choices=['classify', 'direct'], default='direct',
                        help="Kun-peng 分类模式: direct=全量加载(默认, 最快, 需 ≈分片总量的 RAM) / "
                             "classify=分片按需加载(低内存)")
    parser.add_argument('--kunpeng-threads', type=int, default=None,
                        help="批量分类进程的线程数 (默认: --threads × --jobs 的乘积, 如 10×10=100; "
                             "阶段②仍按 --threads/--jobs; 显式传值则覆盖乘积)")
    parser.add_argument('--kunpeng-batch-size', type=int, default=16,
                        help="每批一起分类的样本数 (默认: 16; 0=全部样本一批; 1=每样本单独启动)。"
                             "批量模式下一个 kun_peng 进程处理整批, 索引只加载一次 (热启动), "
                             "这是 Kun-peng 速度优势的前提。")
    parser.add_argument('--splitter', choices=['auto', 'seqkit', 'kractor'], default='auto',
                        help="Kun-peng 的 reads 分流实现: auto=优先 kractor(计数不符/缺失自动回退 seqkit) / "
                             "seqkit=始终用 seqkit grep / kractor=强制 kractor")
    parser.add_argument('--kractor-bin', default=os.environ.get('MMPV_KRACTOR_BIN', 'kractor'),
                        help="补丁版 kractor 可执行文件 (默认 $MMPV_KRACTOR_BIN 或 PATH 中的 kractor)")
    parser.add_argument('--keep-classifier-output', action='store_true',
                        help="保留 Kun-peng 的 per-read 输出与 chunk 中间文件 (体积大, 默认清理)")
    
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
    
    parser.add_argument('--confidence', type=float, default=0.4, help="分类置信度阈值 (kraken2 / kunpeng 通用)")
    parser.add_argument('-f', '--filter', choices=['true', 'false'], default='true', help="true: 去除宿主; false: 提取宿主")
    
    parser.add_argument('--rrna', action='store_true', help="开启 rRNA 剔除")
    parser.add_argument('--rrna_tool', default='ribodetector', choices=['ribodetector', 'silva'],
                        help="rRNA 剔除工具: ribodetector (默认, 仅 rna-short) / silva (Bowtie2 比对 SILVA 库)")
    parser.add_argument('--rrna_impl', default='inprocess', choices=['inprocess', 'subprocess'],
                        help="Ribodetector 实现: inprocess=复用 ONNX 模型(默认, 每线程只加载一次) / "
                             "subprocess=每次启动 ribodetector_cpu (旧行为)")
    parser.add_argument('--rrna_threads', type=int, default=None,
                        help="每个样本的 Ribodetector 推理线程数 (默认: max(1, threads//jobs))")
    parser.add_argument('--silva_index', help="SILVA Bowtie2 索引前缀 (--rrna_tool silva 时必需)")
    parser.add_argument('--keep_rrna', action='store_true', help="将分离出来的rRNA序列保存到 rrna/ 文件夹中 (默认丢弃)")
    parser.add_argument('--keep-taxids', default='10239', help="分类结果中保留的 taxID (逗号分隔), 默认 10239=Viruses; 设为空则不保留任何分类reads")
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

    # --- 智能步骤调整: --keep-taxids 启用且未提供比对索引时自动跳过 align ---
    # 自建宿主库(提供 -x/--step2_index) → Kraken2 + align 双层去除;
    # 通用库如 k2_pluspfp(无比对索引) → 仅 Kraken2 去除。
    if args.keep_taxids and not args.step2_index:
        steps_list = args.steps.split(',')
        if 'align' in steps_list:
            steps_list = [s for s in steps_list if s != 'align']
            args.steps = ','.join(steps_list) if steps_list else ''
            sys.stderr.write("💡 --keep-taxids 已启用且未提供比对索引 (-x/--step2_index)，align 跳过; 提供索引则 Kraken2+align 双层执行\n")
    elif not args.keep_taxids:
        # 显式传空串: 不保留任何分类 reads, 也不跳过 align
        sys.stderr.write("💡 --keep-taxids 为空字符串，已禁用分类 reads 保留 (align 不跳过)\n")

    os.environ['HOST_DEPLETION_STEPS'] = args.steps
    global logger
    logger = configure_logger(__name__, debug=args.debug)

    # 进程守护: 停止主命令时连带杀掉所有子孙 (kraken2/seqtk/worker), 避免后台残留
    try:
        try:
            import process_guard
        except ImportError:
            # process_guard 在 virome_analysis_pipeline/, 单独执行本脚本时不在 sys.path
            _pg_sib = os.path.abspath(os.path.join(
                os.path.dirname(os.path.abspath(__file__)), '..', 'virome_analysis_pipeline'))
            if os.path.isdir(_pg_sib) and _pg_sib not in sys.path:
                sys.path.insert(0, _pg_sib)
            import process_guard
        process_guard.install_global_guard(logger)
    except Exception:
        logger.warning("process_guard 未找到，进程守护未启用 (中断时子孙进程可能残留)")

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
        # --- 输入快照 (硬链接): 并发 clean/下载循环会改写 2.fasta, 分类读到中途消失/截断的
        # --- 文件会整批炸掉。同文件系统硬链接瞬时完成、零拷贝; clean 的删除/重写只影响它自己那份。
        # [P0] 快照刷新机制: 硬链接共享 inode 元数据 —— 源文件重写后 size/mtime 必然变化,
        # 与快照文件 stat 不一致即判定输入已更新: 重链 + 从 .checkpoints 剔除该样本 (下轮重算)。
        # --force 则整目录清空 (配合 --force 已清 .checkpoints), 全量重建。
        try:
            _snap = os.path.join(args.outdir, '.deplete_snapshot')
            if args.force and os.path.isdir(_snap):
                shutil.rmtree(_snap)
                logger.info("🧹 --force: 清空输入快照, 本轮全量重建")
            os.makedirs(_snap, exist_ok=True)
            _linked = _refreshed = 0
            _stale_ids = set()
            for t in tasks:
                _new_ins = []
                for _p in t[0]:
                    _tgt = os.path.join(_snap, os.path.basename(_p))
                    _src_st = os.stat(_p)
                    if os.path.exists(_tgt):
                        _snap_st = os.stat(_tgt)
                        if (_snap_st.st_size, _snap_st.st_mtime_ns) == (_src_st.st_size, _src_st.st_mtime_ns):
                            _new_ins.append(_tgt); _linked += 1
                            continue
                        # 源已重写 → 旧快照作废: 删旧链 + 重链, 并剔除样本 checkpoint
                        try: os.remove(_tgt)
                        except OSError: pass
                        _stale_ids.add(extract_sample_name(os.path.basename(_p)))
                        _refreshed += 1
                    try:
                        os.link(_p, _tgt)
                    except OSError:
                        shutil.copyfile(_p, _tgt)
                    _new_ins.append(_tgt)
                    _linked += 1
                t[0][:] = _new_ins
            if _refreshed:
                logger.info(f"🔄 输入快照刷新: {_refreshed} 个文件源已重写 → 重链; 剔除 {len(_stale_ids)} 个样本 checkpoint 待重算")
            if _stale_ids:
                _ckpt_file = os.path.join(args.outdir, '.checkpoints')
                if os.path.exists(_ckpt_file):
                    with open(_ckpt_file) as f:
                        _keep = [l for l in f if l.strip() not in _stale_ids]
                    with open(_ckpt_file, 'w') as f:
                        f.writelines(_keep)
            logger.info(f"📸 输入快照(硬链接): {_linked} 个文件 → {_snap}")
        except OSError as e:
            logger.warning(f"输入快照失败, 直接使用原路径 ({e})")

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

    classifier_bin = args.kunpeng_bin if args.classifier == 'kunpeng' else 'kraken2'
    rrna_impl = args.rrna_impl
    rrna_threads = args.rrna_threads or max(1, args.threads // max(1, args.jobs))
    if actual_do_rrna and rrna_tool == 'ribodetector' and rrna_impl == 'inprocess':
        try:
            import ribodetector  # noqa: F401
            import onnxruntime  # noqa: F401
        except ImportError as e:
            logger.warning(f"in-process Ribodetector 不可用 ({e})，回退 subprocess")
            rrna_impl = 'subprocess'
    # seqkit: 每步统计 + Kun-peng 的 read 分流; samtools 仅 align 阶段 (view/sort) 需要
    deps = [classifier_bin, 'seqkit']
    if 'align' in args.steps.split(','):
        if not args.tool:
            logger.error("--steps 含 align 时必须提供 --tool"); sys.exit(1)
        deps.append('samtools')
        deps.append(args.tool)
        if not args.step2_index:
            logger.error("--steps 含 align 时必须提供 -x/--step2_index"); sys.exit(1)
    if actual_do_rrna and rrna_tool == 'ribodetector':
        if not args.seq_type:
            logger.error("--steps 含 rrna (ribodetector) 时必须提供 --seq-type"); sys.exit(1)
        if rrna_impl == 'subprocess':
            deps.append('ribodetector_cpu')
    elif actual_do_rrna and rrna_tool == 'silva':
        deps.append('bowtie2')
    # kraken2 用 seqtk 救回 taxid reads; kunpeng 用 seqkit 分流 (seqkit 已在基础依赖中)
    if args.keep_taxids and args.classifier == 'kraken2':
        deps.append('seqtk')
    check_dependencies(deps)
    if 'align' in args.steps.split(',') and not args.step2_index:
        logger.error("--steps 含 align 时必须提供 -x/--step2_index")
        sys.exit(1)
    if args.classifier == 'kunpeng':
        prepare_kunpeng_db(args.kraken2_index, args.kunpeng_bin, args.kunpeng_hash_capacity)
    verify_indices(args.kraken2_index, args.tool, args.step2_index, args.classifier)

    cm = CheckpointManager(checkpoint_dir, force=args.force)

    pe_count = sum(1 for t in tasks if len(t[0]) == 2)
    se_count = sum(1 for t in tasks if len(t[0]) == 1)
    mfilter = True if args.filter == 'true' else False

    print("\n")
    logger.info("============== 智能序列清理管道启动 ==============")
    logger.info(f"数据类型: {args.seq_type.upper() if args.seq_type else 'N/A'} | 第2步工具: {args.tool.upper() if args.tool else 'N/A (已跳过)'}")
    logger.info(f"初筛分类器: {args.classifier.upper()}"
                + (f" ({args.kunpeng_bin}, min-hit-groups={args.min_hit_groups})"
                   if args.classifier == 'kunpeng' else ''))
    tool_label = "SILVA/Bowtie2" if rrna_tool == 'silva' else f"ribodetector ({rrna_impl})"
    if actual_do_rrna: logger.info(f"rRNA过滤: 已开启 (使用 {tool_label})")
    logger.info(f"报表与图表存储: {args.logs_dir}/")
    logger.info(f"解析任务: 包含 {pe_count} 个双端，{se_count} 个单端。")
    logger.info(f"并发配置: {args.jobs} 并发任务 × {args.threads} 线程/任务")
    logger.info("="*60)

    # --- DB 加速层 ---
    # kraken2: 先把库挂载到 /dev/shm 再使用 (标准流程, 不直接用阵列路径)。
    #   mmap 映射 tmpfs 文件时进程直接共享 tmpfs 页 → 零重复内存、免疫页缓存驱逐;
    #   副本跨运行复用 (三件套大小校验一致即跳过重拷), 空间不足时回退 vmtouch 预热/原路径。
    # kunpeng: 库留在原路径 (direct 无 chunk 写洪峰, 页缓存自愈); --vmtouch 只预热分片
    #   hash_*.k2d (双格式库连 hash.k2d 一起热 = 460G 互挤, 反而把分片挤出去)。
    do_k2 = 'kraken2' in args.steps.split(',')
    kr2_shm_dir, kr2_shm_created = None, False
    if do_k2 and args.classifier == 'kraken2':
        db_dir = os.path.expanduser(args.kraken2_index)
        if db_dir.startswith('/dev/shm'):
            logger.info(f"📦 kraken2 库已在 /dev/shm: {db_dir}")
        elif all(os.path.isfile(os.path.join(db_dir, f)) for f in ('hash.k2d', 'opts.k2d', 'taxo.k2d')):
            trio = ['hash.k2d', 'opts.k2d', 'taxo.k2d']
            srcs = [os.path.join(db_dir, f) for f in trio]
            total = sum(os.path.getsize(s) for s in srcs)
            shm_dir = os.path.join('/dev/shm', 'mmpv_kr2_' + os.path.basename(os.path.normpath(db_dir)))

            def _shm_valid():
                return all(os.path.isfile(os.path.join(shm_dir, f)) and
                           os.path.getsize(os.path.join(shm_dir, f)) == os.path.getsize(s)
                           for f, s in zip(trio, srcs))

            if _shm_valid():
                logger.info(f"📦 复用 /dev/shm 里的 kraken2 库副本: {shm_dir} (三件套大小校验一致)")
            else:
                _st = os.statvfs('/dev/shm')
                _free = _st.f_bavail * _st.f_frsize
                if _free < int(total * 1.02):
                    logger.warning(f"/dev/shm 剩余 {_free/(1024**3):.0f}GB < 需要 {total/(1024**3):.0f}GB, "
                                   f"kraken2 直接使用原路径 {db_dir} (可先清理旧副本后重跑)")
                else:
                    os.makedirs(shm_dir, exist_ok=True)
                    logger.info(f"📦 挂载 kraken2 库到 /dev/shm: {shm_dir} "
                                f"({total/(1024**3):.0f}GB, 一次性拷贝; 用完释放: rm -rf {shm_dir})...")
                    cpx_bin = shutil.which('cpx')
                    copied_ok = False
                    # 挂载进度条: 每 30s 报一次已复制量 (cpx/cp 都适用)
                    import threading as _th
                    _stop = _th.Event()

                    def _mount_progress():
                        while not _stop.wait(30):
                            try:
                                _done = sum(os.path.getsize(os.path.join(shm_dir, x))
                                            for x in os.listdir(shm_dir)
                                            if os.path.isfile(os.path.join(shm_dir, x)))
                                logger.info(f"   ⏳ 挂载进度: {_done/(1024**3):.1f}/{total/(1024**3):.1f}GB "
                                            f"({100.0*_done/max(1,total):.1f}%)")
                            except Exception:
                                pass
                    _pt = _th.Thread(target=_mount_progress, daemon=True)
                    _pt.start()
                    if cpx_bin:
                        logger.info("   复制工具: cpx (并行)")
                        _r = subprocess.run([cpx_bin, '-j', '8', '-t', shm_dir, *srcs],
                                            check=False, stdin=subprocess.DEVNULL,
                                            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
                        copied_ok = (_r.returncode == 0 and _shm_valid())
                        if not copied_ok:
                            logger.warning("   cpx 复制失败/校验不一致, 回退 shutil.copyfile")
                    if not cpx_bin:
                        logger.info("   复制工具: shutil.copyfile (未安装 cpx)")
                    if not copied_ok:
                        for f, s in zip(trio, srcs):
                            shutil.copyfile(s, os.path.join(shm_dir, f))
                    _stop.set()
                    _pt.join(timeout=2)
                    if _shm_valid():
                        logger.info("✅ kraken2 库 shm 副本就绪 (mmap 共享 tmpfs 页, 零重复内存)")
                        kr2_shm_dir, kr2_shm_created = shm_dir, True
                    else:
                        logger.warning("shm 副本校验不一致, 回退原路径")
                        shm_dir = None
            if shm_dir:
                args.kraken2_index = shm_dir
            elif args.vmtouch and shutil.which('vmtouch'):
                k2d_files = sorted(glob.glob(os.path.join(db_dir, '*.k2d')))
                if k2d_files:
                    db_size_gb = sum(os.path.getsize(f) for f in k2d_files) / (1024**3)
                    logger.info(f"🔥 vmtouch 预热(回退方案) ({len(k2d_files)} 文件 ~{db_size_gb:.0f}GB)...")
                    for f in k2d_files:
                        subprocess.run(['vmtouch', '-t', f], check=False,
                                      stdin=subprocess.DEVNULL, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        else:
            logger.info(f"kraken2 库 {db_dir} 缺少 hash/opts/taxo 三件套, 跳过 /dev/shm 挂载")
    elif do_k2 and args.vmtouch and shutil.which('vmtouch'):
        db_dir = os.path.expanduser(args.kraken2_index)
        k2d_files = sorted(glob.glob(os.path.join(db_dir, 'hash_[0-9]*.k2d')))
        if k2d_files:
            db_size_gb = sum(os.path.getsize(f) for f in k2d_files) / (1024**3)
            logger.info(f"🔥 vmtouch 预热分片到 page cache ({len(k2d_files)} 文件 ~{db_size_gb:.0f}GB)...")
            for f in k2d_files:
                subprocess.run(['vmtouch', '-t', f], check=False,
                              stdin=subprocess.DEVNULL, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            logger.info("✅ 预热完成 (多进程共享此缓存, 内存仅占一份)")
        else:
            logger.warning(f"vmtouch: 在 {db_dir} 未找到 hash_*.k2d 分片")
    if do_k2 and args.classifier == 'kraken2' and memory_mapping             and not os.path.expanduser(args.kraken2_index).startswith('/dev/shm'):
        logger.info("🧠 memory-mapping 模式 (kraken2 mmap 共享 DB; 首样本自动 fault-in)")
    if do_k2 and args.classifier == 'kunpeng':
        logger.info("🧠 Kun-peng shard 模式 (按需加载 hash_*.k2d; vmtouch 可选)")

    keep_taxids_list = [int(t) for t in str(args.keep_taxids).split(',') if t.strip()] if args.keep_taxids else []
    splitter = args.splitter
    if args.classifier == 'kunpeng' and splitter in ('auto', 'kractor'):
        found = (shutil.which(args.kractor_bin)
                 or (os.path.isfile(args.kractor_bin) and os.access(args.kractor_bin, os.X_OK)))
        if not found:
            if splitter == 'kractor':
                logger.error(f"找不到 kractor: {args.kractor_bin}; "
                             f"请安装补丁版或改用 --splitter seqkit")
                sys.exit(1)
            logger.info("kractor 未找到，splitter 自动回退 seqkit")
            splitter = 'seqkit'
        else:
            logger.info(f"Kun-peng reads 分流: {splitter} (kractor={args.kractor_bin})")
    pipeline = HybridPipeline()
    UI.total_tasks = len(tasks)
    start_time_global = time.time()

    # 预扫描: 已完成的样本不进入分类批
    pending = []
    for idx, (ins, outs) in enumerate(tasks, 1):
        sample_id = extract_sample_name(ins[0])
        if cm.is_done(sample_id):
            logger.info(f"{UI.GRAY}⏭ [跳过] 任务 {idx:02d} | {sample_id} (历史成功){UI.RESET}")
            UI.update_progress(success=True)
            continue
        pending.append((idx, ins, outs, sample_id))

    def _run_group(group, pre_map=None):
        with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as executor:
            futures = {}
            for idx, ins, outs, sample_id in group:
                seq1, seq2 = ins[0], ins[1] if len(ins) > 1 else None
                out1, out2 = outs[0], outs[1] if len(outs) > 1 else None
                # 批模式下用 /dev/shm 预载副本跑阶段② (kractor 分流也吃内存盘)
                _pre = pre_map.get(idx) if pre_map else None
                if _pre and _pre.get('staged_ins'):
                    seq1 = _pre['staged_ins'][0]
                    seq2 = _pre['staged_ins'][1] if len(_pre['staged_ins']) > 1 else None
                future = executor.submit(
                    pipeline.run_single,
                    task_id=idx, step2_tool=args.tool, seq_type=args.seq_type,
                    kr2_idx=args.kraken2_index, step2_idx=args.step2_index,
                    seq1=seq1, seq2=seq2, out1=out1, out2=out2,
                    threads=args.threads, mfilter=mfilter, config=args.config,
                    confidence=args.confidence, tmp_dir=args.tmp, logs_dir=args.logs_dir,
                    do_rrna=actual_do_rrna, keep_rrna=args.keep_rrna, chunk_size=args.chunk_size,
                    rrna_tool=rrna_tool, silva_index=silva_index,
                rrna_impl=rrna_impl, rrna_threads=rrna_threads,
                    keep_taxids=args.keep_taxids, memory_mapping=memory_mapping,
                    classifier=args.classifier, kunpeng_bin=args.kunpeng_bin,
                    min_hit_groups=args.min_hit_groups,
                    keep_classifier_output=args.keep_classifier_output,
                    splitter=splitter, kractor_bin=args.kractor_bin,
                    kunpeng_precomputed=(pre_map.get(idx) if pre_map else None)
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

    # 清扫历史运行残留的 /dev/shm 预载目录 (只清进程已死的; 存活运行的不碰)
    for _d in glob.glob('/dev/shm/kunpeng_stage_*'):
        try:
            _pid = int(os.path.basename(_d).split('_')[2])
            if not os.path.isdir(f'/proc/{_pid}'):
                shutil.rmtree(_d, ignore_errors=True)
                logger.info(f"🧹 清理已死运行残留的 shm 预载目录: {_d}")
        except (ValueError, IndexError):
            pass

    # Kun-peng 批量分类: 一个进程处理整批, 索引/分片只加载一次 (热启动)。
    # 这是论文中 Kun-peng-M/F 速度优势的前提; 每样本单独启动 (batch-size=1)
    # 会退化为冷启动, 大库上固定加载开销会重复 N 次。
    use_batch = (args.classifier == 'kunpeng' and do_k2
                 and args.kunpeng_batch_size != 1 and pending)
    if use_batch:
        size = (args.kunpeng_batch_size if args.kunpeng_batch_size > 0
                else len(pending))
        batches = [pending[i:i + size] for i in range(0, len(pending), size)]
        if args.kunpeng_mode == 'direct':
            db_gb = sum(os.path.getsize(f) for f in
                        glob.glob(os.path.join(args.kraken2_index, 'hash_*.k2d'))) / 1e9
            logger.warning(f"--kunpeng-mode direct: 整库加载，预计峰值内存 ≈ {db_gb:.0f} GB；"
                           f"请确认可用内存足够")
        logger.info(f"Kun-peng 批量分类: {len(pending)} 个样本 / {len(batches)} 批 "
                    f"(每批 ≤{size}, mode={args.kunpeng_mode})")
        for bi, group in enumerate(batches, 1):
            logger.info(f"▶ batch {bi}/{len(batches)}: {len(group)} 个样本，"
                        f"索引只加载一次...")
            try:
                pre_map, batch_dir, stage_dir = kunpeng_batch_classify(
                    group, args.kraken2_index, args.kunpeng_bin,
                    (args.kunpeng_threads or (args.threads * max(1, args.jobs))),
                    args.confidence, keep_taxids_list, args.logs_dir,
                    os.path.join(args.tmp, 'kunpeng_batch'),
                    mode=args.kunpeng_mode,
                    min_hit_groups=args.min_hit_groups,
                    keep_output=args.keep_classifier_output,
                    splitter=splitter, kractor_bin=args.kractor_bin)
            except Exception as e:
                # 单批分类失败不再中止整轮: 该批全部样本记失败 (与 kraken2 单样本失败
                # 同一套记账), 继续下一批; 批内已成功样本的 .checkpoints marker 不受影响,
                # 循环结束后由统一的 exit 1 收尾。
                logger.error(f"Kun-peng 批量分类失败 (batch {bi}/{len(batches)}, "
                             f"{len(group)} 个样本): {e}")
                for _idx, _ins, _outs, _sample_id in group:
                    logger.error(f"[{_sample_id}] 批量分类失败, 记为失败")
                    UI.update_progress(success=False)
                continue
            try:
                _run_group(group, pre_map)
            finally:
                if not args.keep_classifier_output:
                    shutil.rmtree(batch_dir, ignore_errors=True)
                if stage_dir:
                    shutil.rmtree(stage_dir, ignore_errors=True)
    else:
        if args.classifier == 'kunpeng' and do_k2 and args.kunpeng_batch_size == 1:
            logger.warning("--kunpeng-batch-size 1: 每样本单独启动 Kun-peng（冷启动），"
                           "大批量会明显变慢")
        _run_group(pending, None)

    # === 运行结束: 移除本运行创建的 /dev/shm 内容 (按需重建; 复用他人挂载的不动) ===
    if kr2_shm_created and kr2_shm_dir and os.path.isdir(kr2_shm_dir):
        shutil.rmtree(kr2_shm_dir, ignore_errors=True)
        logger.info(f"🧹 运行结束, 已移除本运行创建的 kraken2 shm 挂载: {kr2_shm_dir}")
    for _d in glob.glob(f'/dev/shm/kunpeng_stage_{os.getpid()}_*'):
        shutil.rmtree(_d, ignore_errors=True)

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
            if not generate_plots(summary_tsv_path, plot_prefix):
                logger.warning("⚠️ 可视化图表未生成 (详见上文错误/traceback); 报告与汇总表不受影响")
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
