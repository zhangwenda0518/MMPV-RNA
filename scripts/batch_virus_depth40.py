#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
FastViromeExplorer Pro + QuantPipe 终极统一管线 (V42.4 Resource Monitor Edition)
================================================================================
【核心特性】
1. [智能索引] 自动探测并复用 Reference 所在目录及历史目录的索引，完美跳过重复建库。
2. [泊松打假] 引入强大的 Poisson Ratio 建模，精准剔除Reads局部堆叠引起的假阳性。
3. [全引擎支持] 统一接入传统比对(Bowtie2/BWA等) 与 极速伪比对(Kallisto/Salmon)。
4. [双轨过滤] 完美继承 genes_cov 活跃转录区挽救机制。
5. [断点续传] Batch Parquet 断点保护 + 完善的 Tqdm/File 分流日志系统。
6. [Spawn并发] 强制启用 Spawn 多进程模式，彻底消灭 Polars/C 扩展带来的 Fork 死锁！
7. [资源监控] 自动测算并生成 样本级/全局级 运行时间与内存峰值(Max RSS)硬件报表。
================================================================================
"""

import argparse
import os
import sys
import subprocess
import gzip
import shutil
import time
import re
import datetime
import logging
import platform
import multiprocessing as mp
from pathlib import Path
from collections import defaultdict, Counter
from concurrent.futures import ProcessPoolExecutor, as_completed

try:
    from tqdm import tqdm
except ImportError:
    tqdm = None

try:
    import polars as pl
except ImportError:
    print("❌ 致命错误: 必须安装 polars. (pip install polars)")
    sys.exit(1)

try:
    import pysam
except ImportError:
    print("❌ 致命错误: 必须安装 pysam. (pip install pysam)")
    sys.exit(1)

try:
    import colorlog
except ImportError:
    colorlog = None

# 尝试导入 Unix 系统资源限制模块，用于精确抓取内存
try:
    import resource
except ImportError:
    resource = None

# ==================== R 绘图脚本 ====================
R_PLOT_SCRIPT = r"""#!/usr/bin/env Rscript
suppressWarnings({ suppressPackageStartupMessages({ library(ggplot2); library(dplyr); library(tidyr); library(optparse); library(viridis) }) })
option_list <- list(
  make_option(c("-i", "--input"), type = "character", default = "all_viruses.best.summary.tsv", help = "输入文件"),
  make_option(c("-o", "--output"), type = "character", default = "virus_analysis_plots", help = "输出前缀"),
  make_option(c("-w", "--width"), type = "numeric", default = 10),
  make_option(c("-e", "--height"), type = "numeric", default = 8),
  make_option(c("--log10-transform"), type = "logical", default = FALSE, action = "store_true")
)
opt <- parse_args(OptionParser(option_list = option_list))
if (!file.exists(opt$input)) { stop("输入文件不存在") }
data <- read.delim(opt$input, check.names = FALSE)
if(nrow(data) == 0) { cat("数据为空，跳过绘图\n"); q() }
sp_col <- if("Adjusted_Species" %in% colnames(data)) "Adjusted_Species" else "Species"
data$Display_Name <- paste0(data[[sp_col]], "\n(TaxID: ", data$taxid, ")")
data$Display_Name <- sapply(data$Display_Name, function(x) paste(strwrap(x, width = 45), collapse = "\n"))
metrics <- c("Asm_EM_Reads", "Asm_CPM", "Asm_FPKM", "Avg_Read_ANI", "Poisson_Ratio")
available_metrics <- intersect(metrics, colnames(data))
plot_data <- data %>% select(Display_Name, Sample, all_of(available_metrics)) %>% pivot_longer(cols = all_of(available_metrics), names_to = "Metric", values_to = "Value") %>% filter(!is.na(Value))
if (opt$`log10-transform`) {
  for (m in unique(plot_data$Metric)) { min_v <- min(plot_data$Value[plot_data$Metric==m & plot_data$Value>0], na.rm=TRUE); plot_data$Value[plot_data$Metric==m & plot_data$Value<=0] <- min_v/10 }
}
if ("Asm_EM_Reads" %in% colnames(data)) {
    medians <- aggregate(Value ~ Display_Name, data=plot_data[plot_data$Metric=="Asm_EM_Reads",], median)
    plot_data$Display_Name <- factor(plot_data$Display_Name, levels=medians[order(medians$Value), "Display_Name"])
}
p <- ggplot(plot_data, aes(x=Display_Name, y=Value)) + geom_boxplot(aes(fill=Display_Name), alpha=0.6, outlier.shape=NA) + geom_point(aes(color=Display_Name), position=position_jitter(width=0.2, height=0), alpha=0.6) + facet_wrap(~ Metric, scales="free_x", ncol=length(available_metrics)) + scale_fill_viridis_d(option="turbo") + scale_color_viridis_d(option="turbo") + theme_bw(base_size=13) + theme(legend.position="none", axis.text.y=element_text(size=9, face="italic")) + coord_flip()
if (opt$`log10-transform`) p <- p + scale_y_log10()
ggsave(sprintf("%s_multi_metrics.pdf", opt$output), plot=p, width=opt$width*1.5, height=max(opt$height, length(unique(data$Display_Name))*0.8), dpi=300)
"""

# ==================== 日志与基础工具 ====================
class TqdmLoggingHandler(logging.Handler):
    def emit(self, record):
        try:
            msg = self.format(record)
            if tqdm is not None: tqdm.write(msg)
            else: sys.stdout.write(msg + '\n')
            sys.stdout.flush()
        except Exception:
            self.handleError(record)

def setup_logging(out_root, verbose=False):
    logs_dir = out_root / "logs"
    logs_dir.mkdir(parents=True, exist_ok=True)
    log_file = logs_dir / f"pipeline_{datetime.datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
    logger = logging.getLogger("VPipe_Ultimate")
    logger.setLevel(logging.DEBUG if verbose else logging.INFO)
    logger.handlers = []
    fh = logging.FileHandler(log_file, mode='a', encoding='utf-8')
    fh.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
    logger.addHandler(fh)
    ch = TqdmLoggingHandler()
    if colorlog:
        ch.setFormatter(colorlog.ColoredFormatter('%(log_color)s[%(asctime)s] %(levelname)s - %(message)s', datefmt='%H:%M:%S', log_colors={'DEBUG': 'cyan', 'INFO': 'green', 'WARNING': 'yellow', 'ERROR': 'red'}))
    else:
        ch.setFormatter(logging.Formatter('[%(asctime)s] %(levelname)s - %(message)s', "%H:%M:%S"))
    logger.addHandler(ch)
    return logger, logs_dir

def format_time(seconds):
    m, s = divmod(seconds, 60); h, m = divmod(m, 60)
    if h > 0: return f"{int(h)}h {int(m)}m {int(s)}s"
    elif m > 0: return f"{int(m)}m {int(s)}s"
    return f"{s:.1f}s"

def get_memory_mb():
    if not resource: return 0.0
    # Linux (KB), MacOS (Bytes)
    div = 1024 * 1024 if platform.system() == 'Darwin' else 1024
    ru_self = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    ru_child = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
    return (max(ru_self, ru_child)) / div

def _smart_open(filename, mode='rt'):
    return gzip.open(filename, mode) if str(filename).endswith('.gz') else open(filename, mode)

def _get_ref_read_count_safe(bam, ref):
    try:
        for stat in bam.get_index_statistics():
            if stat.contig == ref: return stat.mapped
    except Exception: pass
    try: return bam.count(contig=ref, read_callback="all")
    except Exception: pass
    return sum(1 for _ in bam.fetch(ref) if not _.is_unmapped)

METRICS_SCHEMA = {
    'Sample': pl.String, 'Accession': pl.String, 'Length': pl.Int64,
    'taxid': pl.String, 'Species': pl.String, 'Segment': pl.String,
    'Coverage(%)': pl.Float64, 'MeanDepth': pl.Float64,
    'EM_Reads': pl.Float64, 'Uniq_Reads': pl.Int64, 'Multi_Reads': pl.Int64,
    'Sample_Total_Mapped': pl.Int64, 'Sample_Global_Reads': pl.Int64,
    'Avg_Read_Len': pl.Float64
}

# ==================== 主管线类 ====================
class UnifiedVirusPipeline:
    def __init__(self, args):
        self.args = args
        self.output_dir = Path(args.output_dir).resolve()
        global logger
        logger, self.logs_dir = setup_logging(self.output_dir, args.verbose)
        self.logger = logger
        
        self.samples = []
        self.index_path = None
        self.ref_length_dict = {}
        self.taxid_clusters = {}  
        self.tax_map = {}
        
        self.is_pseudo = self.args.tool in ['kallisto', 'salmon']
        
        self.check_tools()
        self._load_reference_lengths()
        self._load_taxid_clusters()
        self._load_ref_info_smart()
        
    def write_sample_log(self, sample_name, message, level="info"):
        log_file = self.logs_dir / f"sample_{sample_name}.log"
        try:
            with open(log_file, 'a', encoding='utf-8') as f: 
                f.write(f"[{datetime.datetime.now().strftime('%H:%M:%S')}] {message}\n")
        except Exception:
            pass

    def check_tools(self):
        required_tools = {'samtools': 'BAM处理', 'pandepth': '深度计算'}
        if self.args.use_coverm and not self.is_pseudo: required_tools['coverm'] = '序列清洗'
        tool_map = {'kallisto': 'kallisto', 'salmon': 'salmon', 'bowtie2': 'bowtie2-build', 'bwa': 'bwa', 'bwa-mem2': 'bwa-mem2', 'hisat2': 'hisat2-build', 'minimap2': 'minimap2', 'strobealign': 'strobealign'}
        if self.args.tool in tool_map: required_tools[tool_map[self.args.tool]] = f'{self.args.tool} 核心程序'
        missing = [f"{t} ({d})" for t, d in required_tools.items() if shutil.which(t) is None]
        if missing: 
            self.logger.error("❌ 缺少必要工具:"); [self.logger.error(f"  - {m}") for m in missing]; sys.exit(1)

    def _load_reference_lengths(self):
        with _smart_open(self.args.reference, 'rt') as f:
            curr_id, curr_len = None, 0
            for line in f:
                if line.startswith('>'):
                    if curr_id: self.ref_length_dict[curr_id] = curr_len
                    curr_id, curr_len = line.strip().split()[0][1:], 0
                else: curr_len += len(line.strip())
            if curr_id: self.ref_length_dict[curr_id] = curr_len

    def _load_taxid_clusters(self):
        if getattr(self.args, 'taxid_clusters', None) and os.path.exists(self.args.taxid_clusters):
            with open(self.args.taxid_clusters, 'r') as f:
                for line in f:
                    if line.startswith('object_taxid'): continue
                    parts = line.strip().split('\t')
                    if len(parts) >= 2: self.taxid_clusters[parts[0]] = parts[1]

    def _load_ref_info_smart(self):
        if not self.args.ref_info or not os.path.exists(self.args.ref_info):
            self.logger.error("❌ 必须提供 --ref_info (本地 TSV)")
            sys.exit(1)
            
        acc_synonyms = ['Accession', 'accession', 'Virus GENBANK accession', 'ID']
        tax_synonyms = ['Taxid', 'taxonomy_id', 'taxid']
        sp_synonyms = ['Species_NCBI', 'Species_ICTV', 'taxonomy', 'Species', 'description', 'Virus name(s)']
        seg_synonyms = ['Segment', 'segment']

        with open(self.args.ref_info, 'r', encoding='utf-8') as f:
            header, idx_acc, idx_tax, idx_sp, idx_seg = None, -1, -1, -1, -1
            for line in f:
                if line.startswith('####') or not line.strip(): continue
                parts = line.rstrip('\n').split('\t')
                if header is None:
                    header = [h.strip() for h in parts]
                    idx_acc = next((header.index(s) for s in acc_synonyms if s in header), -1)
                    idx_tax = next((header.index(s) for s in tax_synonyms if s in header), -1)
                    idx_sp = next((header.index(s) for s in sp_synonyms if s in header), -1)
                    idx_seg = next((header.index(s) for s in seg_synonyms if s in header), -1)
                    if idx_acc == -1: self.logger.error(f"❌ 找不到 Accession 列: {header[:5]}"); sys.exit(1)
                    continue
                if len(parts) <= idx_acc: continue
                acc = parts[idx_acc].strip()
                info = {
                    'taxid': self.taxid_clusters.get(parts[idx_tax].strip(), parts[idx_tax].strip()) if idx_tax != -1 and len(parts) > idx_tax else "Unannotated",
                    'species': parts[idx_sp].strip() if idx_sp != -1 and len(parts) > idx_sp else "Unannotated",
                    'segment': parts[idx_seg].strip() if idx_seg != -1 and len(parts) > idx_seg else ""
                }
                self.tax_map[acc] = info
                self.tax_map[acc.split('.')[0]] = info

    def setup_output_directory(self):
        subdirs = ['bam', 'stat', 'summary', 'index', 'plots', 'batches']
        if self.is_pseudo: subdirs.append('pseudo_quant')
        for subdir in subdirs: (self.output_dir / subdir).mkdir(parents=True, exist_ok=True)

    def build_index(self):
        tool = self.args.tool
        ref_path = Path(self.args.reference).absolute()
        
        if tool in ['strobealign', 'minimap2']:
            self.index_path = str(ref_path); return

        def _check_index_exists(prefix_str, ext_list, is_dir=False):
            if not prefix_str: return False
            p = Path(prefix_str)
            if is_dir: return p.exists() and p.is_dir()
            return any(Path(f"{prefix_str}{ext}").exists() for ext in ext_list)

        found_existing = False
        potential_prefixes = [str(ref_path), str(ref_path.with_suffix(''))]
        
        if self.is_pseudo:
            pseudo_idx_path = ref_path.with_suffix(f'.{tool}_idx')
            if tool == 'salmon' and pseudo_idx_path.exists() and pseudo_idx_path.is_dir(): self.index_path = str(pseudo_idx_path); found_existing = True
            elif tool == 'kallisto' and pseudo_idx_path.exists() and pseudo_idx_path.is_file(): self.index_path = str(pseudo_idx_path); found_existing = True
            if not found_existing:
                out_idx = self.output_dir / 'index' / f"{ref_path.name}.{tool}_idx"
                if tool == 'salmon' and out_idx.exists() and out_idx.is_dir(): self.index_path = str(out_idx); found_existing = True
                elif tool == 'kallisto' and out_idx.exists() and out_idx.is_file(): self.index_path = str(out_idx); found_existing = True
        else:
            tool_exts = {'bwa': ['.bwt'], 'bwa-mem2': ['.bwt.2bit.64', '.0123'], 'bowtie2': ['.1.bt2', '.1.bt2l'], 'hisat2': ['.1.ht2', '.1.ht2l']}
            exts_to_check = tool_exts.get(tool, [])
            for pref in potential_prefixes:
                if _check_index_exists(pref, exts_to_check): self.index_path = pref; found_existing = True; break
            if not found_existing:
                out_prefix = str(self.output_dir / 'index' / f"{ref_path.stem}_{tool}")
                if _check_index_exists(out_prefix, exts_to_check): self.index_path = out_prefix; found_existing = True

        if found_existing:
            self.logger.info(f"✅ [智能探测] 检测到已存在的 {tool} 索引，跳过构建 -> {self.index_path}")
            return

        self.logger.info(f"🏗️ [构建索引] 未检测到 {tool} 索引，准备新建 (这可能需要一些时间)...")
        if self.is_pseudo:
            new_idx_path = self.output_dir / 'index' / f"{ref_path.name}.{tool}_idx"
            self.index_path = str(new_idx_path)
            idx_threads = str(min(self.args.threads, 16))
            cmd = ['salmon', 'index', '-k', '31', '-i', self.index_path, '-t', str(ref_path), '--threads', idx_threads] if tool == 'salmon' else ['kallisto', 'index', '-k', '31', '-i', self.index_path, '--threads', idx_threads, str(ref_path)]
        else:
            prefix = self.output_dir / 'index' / f"{ref_path.stem}_{tool}"
            self.index_path = str(prefix)
            cmds = {'bwa': ['bwa', 'index', '-p', self.index_path, str(ref_path)], 'bwa-mem2': ['bwa-mem2', 'index', '-p', self.index_path, str(ref_path)], 'bowtie2': ['bowtie2-build', '--threads', str(self.args.threads), str(ref_path), self.index_path], 'hisat2': ['hisat2-build', '-p', str(self.args.threads), str(ref_path), self.index_path]}
            cmd = cmds.get(tool)
            
        index_log_file = self.logs_dir / f"index_build_{tool}.log"
        with open(index_log_file, "w") as flog: subprocess.run(cmd, stdout=flog, stderr=subprocess.STDOUT, check=True)
        self.logger.info(f"✅ [构建完毕] 新 {tool} 索引已成功生成至 -> {self.index_path}")

    def find_samples(self):
        input_dir = Path(self.args.input_dir) if self.args.input_dir else None
        valid_exts = ['.fastq.gz', '.fq.gz', '.fasta.gz', '.fa.gz', '.fastq', '.fq', '.fasta', '.fa']
        samples = []
        if input_dir:
            file_dict = defaultdict(dict)
            for f in input_dir.iterdir():
                if f.is_dir(): continue
                ext_found = next((ext for ext in valid_exts if f.name.lower().endswith(ext)), None)
                if not ext_found: continue
                base_name = f.name[:-len(ext_found)]
                
                # 修复1: 强制要求必须包含分隔符 [._-]+ 防止切断带数字的 SRR 编号 (例如 SRR15037502 变成 SRR1503750)
                m = re.search(r'([._\-]+)(R?[12])(?:_\d+)?$', base_name, re.IGNORECASE)
                if m and m.group(2): 
                    direction = 'R1' if '1' in m.group(2) else 'R2'
                    sname = base_name[:m.start()] or base_name
                else: 
                    direction = 'SE'
                    sname = base_name
                
                # 修复2: 加入 _clean 和 .clean 的匹配与清理
                sname = re.sub(r'(_unmapped|\.unmapped|_trimmed|\.trimmed|_clean|\.clean)$', '', sname, flags=re.IGNORECASE) or sname
                
                file_dict[sname][direction] = str(f)

            for sname, reads in file_dict.items():
                if self.args.single_end:
                    for d, path in reads.items(): samples.append({'name': sname if len(reads) == 1 else f"{sname}_{d}", 'r1': path, 'r2': None})
                else:
                    if 'R1' in reads and 'R2' in reads: samples.append({'name': sname, 'r1': reads['R1'], 'r2': reads['R2']})
                    elif 'R1' in reads: samples.append({'name': sname, 'r1': reads['R1'], 'r2': None})
                    elif 'R2' in reads: samples.append({'name': sname, 'r1': reads['R2'], 'r2': None})
                    elif 'SE' in reads: samples.append({'name': sname, 'r1': reads['SE'], 'r2': None})

        if self.args.sample_list:
            with open(self.args.sample_list, 'r') as f:
                for line in f:
                    parts = line.strip().split()
                    if len(parts) >= 2 and not line.startswith('#') and Path(parts[1]).exists():
                        samples.append({'name': parts[0], 'r1': parts[1], 'r2': parts[2] if len(parts)>2 and Path(parts[2]).exists() else None})
        
        self.samples = list({s['name']: s for s in samples}.values())
        if not self.samples: self.logger.error("❌ 未找到样本文件。"); sys.exit(1)
        
        # 🟢 [新增] 打印数据类型统计
        se_count = sum(1 for s in self.samples if not s.get('r2'))
        pe_count = len(self.samples) - se_count
        self.logger.info("📊 [输入数据统计]")
        self.logger.info(f"   ├─ 总计发现样本: {len(self.samples)} 个")
        self.logger.info(f"   ├─ 单端样本 (SE): {se_count} 个")
        self.logger.info(f"   └─ 双端样本 (PE): {pe_count} 个")

    def get_average_read_length(self, file_path):
        is_fasta = any(file_path.lower().endswith(ext) for ext in ['.fa', '.fasta', '.fa.gz', '.fasta.gz'])
        awk_script = "'/^>/ {if (seqlen > 0) { lenSum+=seqlen; readCount++; seqlen=0 }} !/^>/ { seqlen += length($0) } END { if (seqlen > 0) { lenSum+=seqlen; readCount++; }; if (readCount > 0) print lenSum/readCount; else print 0; }'" if is_fasta else "'NR%4 == 2 {lenSum+=length($0); readCount++;} END {if (readCount > 0) print lenSum/readCount; else print 0}'"
        cmd = f"gzip -dc '{file_path}' | awk {awk_script}" if file_path.endswith(".gz") else f"awk {awk_script} '{file_path}'"
        try:
            res = subprocess.run(cmd, shell=True, capture_output=True, text=True)
            return float(res.stdout.strip())
        except Exception: return 150.0

    def process_sample(self, sample):
        t0 = time.time()
        try:
            if self.is_pseudo:
                res = self.process_sample_pseudo(sample)
            else:
                res = self.process_sample_traditional(sample)
        except Exception as e:
            self.write_sample_log(sample['name'], f"❌ 处理失败: {str(e)}", "error")
            return []
            
        t1 = time.time()
        elapsed = t1 - t0
        
        # 从 /usr/bin/time -v 日志提取准确的峰值内存和CPU（比 getrusage 更可靠）
        worker_mem_mb = 0.0
        cpu_pct = 0.0
        slf = self.logs_dir / f"sample_{sample['name']}.log"
        if slf.exists():
            try:
                with open(slf) as f:
                    txt = f.read()
                    m_mem = re.search(r'Maximum resident set size \(kbytes\):\s+(\d+)', txt)
                    if m_mem: worker_mem_mb = float(m_mem.group(1)) / 1024.0
                    m_cpu = re.search(r'Percent of CPU this job got:\s+(\d+)', txt)
                    if m_cpu: cpu_pct = float(m_cpu.group(1))
            except: pass
        self.write_sample_log(sample['name'], f"✅ [性能监控] 耗时: {format_time(elapsed)} | 峰值内存(MAX_RSS): {worker_mem_mb:.2f} MB | CPU利用率: {cpu_pct:.0f}%")
        return res

    def process_sample_pseudo(self, sample):
        sname = sample['name']; bam_output = self.output_dir / 'bam' / f'{sname}.sorted.bam'
        quant_dir = self.output_dir / 'pseudo_quant' / sname; quant_dir.mkdir(exist_ok=True)
        threads = self.args.align_threads; avg_len = self.get_average_read_length(sample['r1'])
        samtools_pipe = f"samtools view -@ {min(2, threads)} -b -F 0x04 - | samtools sort -@ {min(2, threads)} -o '{bam_output}'"
        
        if self.args.tool == 'salmon':
            read_arg = f"-1 '{sample['r1']}' -2 '{sample['r2']}'" if sample['r2'] else f"-r '{sample['r1']}'"
            full_cmd = f"/usr/bin/time -v salmon quant -i '{self.index_path}' -p {threads} -l A {read_arg} -o '{quant_dir}' --writeMappings | {samtools_pipe}"
            sf_file = quant_dir / "quant.sf"
        else:
            quant_args = f"--pseudobam '{sample['r1']}' '{sample['r2']}'" if sample['r2'] else f"--single -l 200 -s 50 --pseudobam '{sample['r1']}'"
            full_cmd = f"/usr/bin/time -v kallisto quant -i '{self.index_path}' -t {threads} -o '{quant_dir}' {quant_args} | {samtools_pipe}"
            sf_file = quant_dir / "abundance.tsv"
            
        with open(self.logs_dir / f"sample_{sname}.log", "a") as flog:
            subprocess.run(full_cmd, shell=True, executable='/bin/bash', check=True, stderr=flog)
            pysam.index(str(bam_output))
            
        stat_prefix = self.output_dir / 'stat' / sname
        subprocess.run(['pandepth', '-a', '-i', str(bam_output), '-o', str(stat_prefix), '-t', str(threads) ], capture_output=True)
        
        df_quant = pl.read_csv(sf_file, separator='\t').select([pl.col("Name" if self.args.tool == 'salmon' else "target_id").alias("Accession"), pl.col("NumReads" if self.args.tool == 'salmon' else "est_counts").alias("Reads")]).filter(pl.col("Reads") > 0)
        
        depth_dict = {}
        stat_file = str(stat_prefix) + '.chr.stat.gz'
        if os.path.exists(stat_file):
            try:
                with gzip.open(stat_file, 'rt') as f:
                    for line in f:
                        if not line.startswith('#'):
                            p = line.split('\t')
                            if len(p) >= 6: depth_dict[p[0].strip()] = {'Cov': float(p[4]), 'Dep': float(p[5])}
            except Exception: pass
            
        total_mapped = df_quant["Reads"].sum(); global_reads = total_mapped 

        results = []
        for row in df_quant.iter_rows(named=True):
            acc = row['Accession']; reads = row['Reads']; d_info = depth_dict.get(acc, {'Cov': 0.0, 'Dep': 0.0})
            results.append({'Sample': sname, 'Accession': acc, 'Length': self.ref_length_dict.get(acc, 1), 'taxid': self.tax_map.get(acc, {}).get('taxid', 'Unannotated'), 'Species': self.tax_map.get(acc, {}).get('species', 'Unannotated'), 'Segment': self.tax_map.get(acc, {}).get('segment', ''), 'Coverage(%)': d_info['Cov'], 'MeanDepth': d_info['Dep'], 'EM_Reads': float(reads), 'Uniq_Reads': int(reads), 'Multi_Reads': 0, 'Sample_Total_Mapped': int(total_mapped), 'Sample_Global_Reads': int(global_reads), 'Avg_Read_Len': avg_len})
            
        if not self.args.keep_tmp:
            try: shutil.rmtree(quant_dir)
            except Exception: pass
        return results

    def process_sample_traditional(self, sample):
        sname = sample['name']; raw_bam = self.output_dir / 'bam' / f'{sname}.raw.bam'; filt_bam = self.output_dir / 'bam' / f'{sname}.sorted.bam'
        threads = self.args.align_threads; avg_len = self.get_average_read_length(sample['r1'])
        
        cmd_str = ""
        if self.args.tool == 'strobealign': cmd_str = f"/usr/bin/time -v strobealign -t {threads} '{self.index_path}' '{sample['r1']}' {' ' + repr(sample['r2']) if sample['r2'] else ''}"
        elif self.args.tool == 'minimap2': cmd_str = f"/usr/bin/time -v minimap2 -ax {'sr' if sample['r2'] else 'map-ont'} -t {threads} '{self.index_path}' '{sample['r1']}' {' ' + repr(sample['r2']) if sample['r2'] else ''}"
        elif self.args.tool in ['bwa', 'bwa-mem2']: cmd_str = f"/usr/bin/time -v {self.args.tool} mem -v 1 -t {threads} '{self.index_path}' '{sample['r1']}' {' ' + repr(sample['r2']) if sample['r2'] else ''}"
        elif self.args.tool == 'bowtie2':
            fmt = '-f' if any(sample['r1'].lower().endswith(ext) for ext in ['.fa', '.fasta', '.fa.gz', '.fasta.gz']) else ''
            io_args = f"-1 '{sample['r1']}' -2 '{sample['r2']}'" if sample['r2'] else f"-U '{sample['r1']}'"
            cmd_str = f"/usr/bin/time -v bowtie2 -p {threads} -x '{self.index_path}' {fmt} {io_args}"
        elif self.args.tool == 'hisat2':
            fmt = '-f' if any(sample['r1'].lower().endswith(ext) for ext in ['.fa', '.fasta', '.fa.gz', '.fasta.gz']) else ''
            io_args = f"-1 '{sample['r1']}' -2 '{sample['r2']}'" if sample['r2'] else f"-U '{sample['r1']}'"
            cmd_str = f"/usr/bin/time -v hisat2 -p {threads} -x '{self.index_path}' {fmt} {io_args}"
            
        with open(self.logs_dir / f"sample_{sname}.log", "a") as flog:
            subprocess.run(f"set -o pipefail; {cmd_str} | samtools view -b -o '{raw_bam}'", shell=True, executable='/bin/bash', check=True, stderr=flog)
            
        if self.args.use_coverm:
            filter_cmd = ['coverm', 'filter', '-b', str(raw_bam), '-o', str(filt_bam) + '.unsorted', '--min-read-aligned-length', str(self.args.min_aln_len), '--min-read-percent-identity', str(self.args.min_pid), '--min-read-aligned-percent', str(self.args.min_aln_prop), '--include-secondary', '-t', '1']
            if sample['r2']: filter_cmd.append('--proper-pairs-only')
            subprocess.run(filter_cmd, capture_output=True)
            subprocess.run(['samtools', 'sort', '-@', '1', '-o', str(filt_bam), str(filt_bam) + '.unsorted'])
        else:
            subprocess.run(f"samtools view -F 4 -b '{raw_bam}' | samtools sort -@ 1 -o '{filt_bam}'", shell=True, executable='/bin/bash')
            
        pysam.index(str(filt_bam))
        if not self.args.keep_tmp:
            if raw_bam.exists(): raw_bam.unlink()
            if Path(str(filt_bam) + '.unsorted').exists(): Path(str(filt_bam) + '.unsorted').unlink()
            
        qc_res = subprocess.run(['samtools', 'view', '-c', str(filt_bam)], capture_output=True, text=True)
        total_mapped = int(qc_res.stdout.strip() or 1); global_reads = total_mapped
        
        read_best_data = {}
        with pysam.AlignmentFile(filt_bam, "rb") as bam:
            for read in bam.fetch(until_eof=True):
                if read.is_unmapped: continue
                ref = bam.get_reference_name(read.reference_id)
                taxid = self.tax_map.get(ref, {}).get('taxid', ref)
                try: score = read.get_tag('AS')
                except KeyError: score = read.query_alignment_length - read.get_tag('NM', 0)
                
                qname = read.query_name
                if qname not in read_best_data: read_best_data[qname] = {'max_score': score, 'refs': {ref: score}, 'taxids': {taxid}}
                else:
                    curr_max = read_best_data[qname]['max_score']
                    if score > curr_max: read_best_data[qname] = {'max_score': score, 'refs': {ref: score}, 'taxids': {taxid}}
                    elif score == curr_max: read_best_data[qname]['refs'][ref] = score; read_best_data[qname]['taxids'].add(taxid)

        em_counts, uniq_counts, multi_counts = defaultdict(float), defaultdict(int), defaultdict(int)
        uniq_taxid_counts = Counter()
        for data in read_best_data.values():
            if len(data['taxids']) == 1: uniq_taxid_counts[list(data['taxids'])[0]] += 1

        for data in read_best_data.values():
            taxids, refs_dict = list(data['taxids']), data['refs']
            refs = list(refs_dict.keys())
            if len(taxids) == 1:
                total_score = sum(refs_dict.values())
                for r in refs:
                    weight = (refs_dict[r] / total_score) if total_score > 0 else (1.0 / len(refs))
                    em_counts[r] += weight; uniq_counts[r] += 1
            else:
                total_evidence = sum(uniq_taxid_counts[tid] for tid in taxids)
                if total_evidence > 0:
                    for r in refs:
                        tid = self.tax_map.get(r, {}).get('taxid', r)
                        ev = uniq_taxid_counts[tid]
                        if ev > 0:
                            ref_n = sum(1 for x in refs if self.tax_map.get(x, {}).get('taxid', x) == tid)
                            em_counts[r] += (ev / total_evidence) / ref_n; multi_counts[r] += 1

        stat_prefix = self.output_dir / 'stat' / sname
        subprocess.run(['pandepth', '-a', '-i', str(filt_bam), '-o', str(stat_prefix), '-t', str(threads) ], capture_output=True)
        depth_dict = {}
        stat_file = str(stat_prefix) + '.chr.stat.gz'
        if os.path.exists(stat_file):
            try:
                with gzip.open(stat_file, 'rt') as f:
                    for line in f:
                        if not line.startswith('#'):
                            p = line.split()
                            if len(p) >= 6: depth_dict[p[0]] = {'Cov': float(p[4]), 'Dep': float(p[5])}
            except Exception: pass
            
        results = []
        for virus, em_mapped in em_counts.items():
            if em_mapped <= 0: continue
            d_info = depth_dict.get(virus, {'Cov': 0, 'Dep': 0})
            results.append({'Sample': sname, 'Accession': virus, 'Length': self.ref_length_dict.get(virus, 1), 'taxid': self.tax_map.get(virus, {}).get('taxid', 'Unannotated'), 'Species': self.tax_map.get(virus, {}).get('species', 'Unannotated'), 'Segment': self.tax_map.get(virus, {}).get('segment', ''), 'Coverage(%)': float(d_info['Cov']), 'MeanDepth': float(d_info['Dep']), 'EM_Reads': float(em_mapped), 'Uniq_Reads': int(uniq_counts.get(virus, 0)), 'Multi_Reads': int(multi_counts.get(virus, 0)), 'Sample_Total_Mapped': int(total_mapped), 'Sample_Global_Reads': int(global_reads), 'Avg_Read_Len': avg_len})
        return results

    def _compute_ani_pi_worker(self, task):
        sname, bam_path, ref = task
        if not os.path.exists(bam_path): return {'Sample': sname, 'Rep_Accession': ref, 'Avg_Read_ANI': None, 'Avg_Pi': None}
        ani_sum, ani_cnt, base_counts = 0.0, 0, defaultdict(Counter)
        try:
            with pysam.AlignmentFile(bam_path, "rb") as bam:
                if ref not in bam.references: return {'Sample': sname, 'Rep_Accession': ref, 'Avg_Read_ANI': None, 'Avg_Pi': None}
                total_reads = _get_ref_read_count_safe(bam, ref)
                if total_reads == 0: return {'Sample': sname, 'Rep_Accession': ref, 'Avg_Read_ANI': None, 'Avg_Pi': None}
                step = max(1, total_reads // 10000)
                for idx, read in enumerate(bam.fetch(ref)):
                    if idx % step != 0 or read.is_unmapped: continue
                    aln_len = read.query_alignment_length
                    if aln_len > 0:
                        try: nm = read.get_tag('NM')
                        except KeyError: nm = 0
                        ani_sum += (aln_len - nm) / aln_len; ani_cnt += 1
                    for qpos, rpos, ref_base in read.get_aligned_pairs(matches_only=True, with_seq=True):
                        if rpos is not None and qpos is not None: base_counts[rpos][read.query_sequence[qpos].upper()] += 1
                    if ani_cnt >= 10000: break
        except Exception: pass
        avg_ani = (ani_sum / ani_cnt) * 100.0 if ani_cnt > 0 else None
        pi_sum, covered_pos = 0.0, 0
        for rpos, counts in base_counts.items():
            total = sum(counts.values())
            if total > 1: pi_sum += 1.0 - sum((c/total)**2 for c in counts.values()); covered_pos += 1
        avg_pi = (pi_sum / covered_pos) if covered_pos > 0 else None
        return {'Sample': sname, 'Rep_Accession': ref, 'Avg_Read_ANI': round(avg_ani, 2) if avg_ani else None, 'Avg_Pi': round(avg_pi, 5) if avg_pi else None}

    def summarize_results_polars(self, df):
        if len(df) == 0:
            self.logger.warning("⚠️ 没有产生任何有效结果。")
            return
            
        df = df.with_columns([
            (pl.col("EM_Reads") * 1000.0 / pl.col("Length")).alias("Seq_RPK"),
            pl.when(pl.col("Segment") != "").then(pl.col("Segment") + ":" + pl.col("Accession")).otherwise(pl.col("Accession")).alias("Seg_Acc_Str")
        ])
        df.write_csv(str(self.output_dir / "summary" / f"all_viruses.raw.tsv"), separator='\t')
        
        best_acc_df = df.sort(["Sample", "taxid", "EM_Reads", "Coverage(%)"], descending=[False, False, True, True]).group_by(["Sample", "taxid"]).first()
        rep_df = pl.DataFrame({"Sample": best_acc_df["Sample"], "taxid": best_acc_df["taxid"], "Rep_Accession": best_acc_df["Accession"], "Rep_Reads": best_acc_df["EM_Reads"], "Base_Parsed_Species": best_acc_df["Species"]})
        
        sample_meta_df = df.group_by("Sample").agg([pl.col("Seq_RPK").sum().alias("Sample_Total_RPK"), pl.col("Sample_Total_Mapped").first(), pl.col("Sample_Global_Reads").first(), pl.col("Avg_Read_Len").first()])
        
        asm_df = df.group_by(["Sample", "taxid"]).agg([
            pl.col("EM_Reads").sum().alias("Asm_EM_Reads"), pl.col("Seq_RPK").sum().alias("Asm_RPK_Sum"),
            pl.col("Uniq_Reads").sum(), pl.col("Multi_Reads").sum(),
            pl.col("Length").sum().alias("Asm_Length"), pl.col("Seg_Acc_Str").unique().str.join(",").alias("Segment_Accessions")
        ])
        
        asm_df = asm_df.join(rep_df, on=["Sample", "taxid"], how="left").join(sample_meta_df, on="Sample", how="left")
        
        asm_df = asm_df.with_columns([
            pl.when((pl.col("Uniq_Reads") + pl.col("Multi_Reads")) > 0).then((pl.col("Uniq_Reads") / (pl.col("Uniq_Reads") + pl.col("Multi_Reads")) * 100).round(2)).otherwise(0.0).alias("Unique(%)"),
            (pl.col("Asm_EM_Reads") * 1e6 / pl.col("Sample_Total_Mapped")).round(2).alias("Asm_CPM"),
            (pl.col("Asm_EM_Reads") * 1e6 / pl.col("Sample_Global_Reads")).round(2).alias("Asm_RPM"),
            pl.when(pl.col("Asm_Length") > 0).then((pl.col("Asm_EM_Reads") * 1e9 / (pl.col("Sample_Total_Mapped") * pl.col("Asm_Length")))).otherwise(0.0).round(2).alias("Asm_FPKM"),
            pl.when(pl.col("Sample_Total_RPK") > 0).then((pl.col("Asm_RPK_Sum") * 1e6 / pl.col("Sample_Total_RPK")).round(2)).otherwise(0.0).alias("Asm_TPM"),
            (pl.col("Asm_EM_Reads") / pl.col("Sample_Total_Mapped") * 100.0).round(4).alias("Asm_Rel_Abund(%)")
        ])
        
        rep_seq_stats = df.select(["Sample", "Accession", "Length", "Coverage(%)", "MeanDepth"]).rename({"Accession": "Rep_Accession", "Coverage(%)": "Rep_Coverage(%)", "MeanDepth": "Rep_MeanDepth", "Length": "Rep_Length"}).unique(subset=["Sample", "Rep_Accession"], keep="first", maintain_order=True)
        asm_df = asm_df.join(rep_seq_stats, on=["Sample", "Rep_Accession"], how="left")

        asm_df = asm_df.with_columns([(1.0 - (-(pl.col("Rep_Reads") * pl.col("Avg_Read_Len")) / pl.col("Rep_Length")).exp()).alias("Predicted_Support")])
        asm_df = asm_df.with_columns([pl.when(pl.col("Predicted_Support") > 0).then((pl.col("Rep_Coverage(%)") / 100.0) / pl.col("Predicted_Support")).otherwise(0.0).alias("Poisson_Ratio")])
        
        ref_meta_df = pl.read_csv(self.args.ref_info, separator='\t', ignore_errors=True, truncate_ragged_lines=True, quote_char=None)
        acc_col_name = next((col for col in ['Accession', 'accession', 'Virus GENBANK accession', 'ID'] if col in ref_meta_df.columns), None)
        
        if acc_col_name:
            ref_meta_df = ref_meta_df.with_columns(pl.col(acc_col_name).cast(pl.Utf8).str.strip_chars().str.replace(r"\.\d+$", "").alias("_safe_acc"))
            asm_df = asm_df.with_columns(pl.col("Rep_Accession").str.replace(r"\.\d+$", "").alias("_safe_acc"))
            asm_df = asm_df.join(ref_meta_df, on="_safe_acc", how="left").drop("_safe_acc")

        filter_base = (pl.col("Sample_Total_Mapped") > 0) & (pl.col("Uniq_Reads") >= self.args.min_uniq_reads) & (pl.col("Asm_TPM") >= self.args.min_tpm)
        track_a_pass = (pl.col("Rep_Coverage(%)") >= self.args.coverage) & (pl.col("Poisson_Ratio") >= self.args.ratio) & (pl.col("Rep_MeanDepth") >= self.args.meandepth)
        track_b_pass = pl.lit(False)
        if self.args.genes_cov and os.path.exists(self.args.genes_cov):
            try:
                genes_df = pl.read_csv(self.args.genes_cov, separator='\t', ignore_errors=True)
                if 'seqid' in genes_df.columns and 'gene_total_cov' in genes_df.columns and 'gene_avr_cov' in genes_df.columns:
                    asm_df = asm_df.join(genes_df.select(['seqid', 'gene_total_cov', 'gene_avr_cov']), left_on="Rep_Accession", right_on="seqid", how="left")
                    track_b_pass = (pl.col("gene_total_cov").fill_null(0) >= self.args.min_gene_total_cov) & (pl.col("gene_avr_cov").fill_null(0) >= self.args.min_gene_avr_cov)
            except Exception: pass
                
        final_df = asm_df.filter(filter_base & (track_a_pass | track_b_pass))
        final_df.write_csv(str(self.output_dir / "summary" / "all_viruses.summary.tsv"), separator='\t')
        pre_best_df = final_df.sort(["Sample", "taxid", "Asm_EM_Reads"], descending=[False, False, True]).unique(subset=["Sample", "taxid"], keep="first", maintain_order=True)
        
        if not self.is_pseudo:
            self.logger.info(f"🧬 [进化测算] 正在对存活的 {len(pre_best_df)} 个核心代表株进行深度进化测算 (ANI/Pi)...")
            tasks = [(row['Sample'], str(self.output_dir / 'bam' / f"{row['Sample']}.sorted.bam"), row['Rep_Accession']) for row in pre_best_df.iter_rows(named=True)]
            ani_pi_results = []
            ctx = mp.get_context('spawn')
            with ProcessPoolExecutor(max_workers=self.args.threads, mp_context=ctx) as executor:
                futures = {executor.submit(self._compute_ani_pi_worker, task): i for i, task in enumerate(tasks)}
                if tqdm:
                    with tqdm(total=len(futures), desc="ANI测算", dynamic_ncols=True, colour="blue") as pbar:
                        for future in as_completed(futures):
                            if (res := future.result()): ani_pi_results.append(res)
                            pbar.update(1)
                else:
                    for future in as_completed(futures):
                        if (res := future.result()): ani_pi_results.append(res)

            ap_df = pl.DataFrame(ani_pi_results) if ani_pi_results else pl.DataFrame(schema={"Sample": pl.String, "Rep_Accession": pl.String, "Avg_Read_ANI": pl.Float64, "Avg_Pi": pl.Float64})
            merged_df = pre_best_df.join(ap_df, on=["Sample", "Rep_Accession"], how="left")
        else:
            merged_df = pre_best_df.with_columns([pl.lit(None).cast(pl.Float64).alias("Avg_Read_ANI"), pl.lit(None).cast(pl.Float64).alias("Avg_Pi")])

        sp_candidates = ['Species_NCBI', 'Species_ICTV', 'Base_Parsed_Species']
        available_sp_cols = [c for c in sp_candidates if c in merged_df.columns]
        
        if available_sp_cols:
            merged_df = merged_df.with_columns(pl.coalesce([pl.col(c) for c in available_sp_cols]).fill_null("Unknown").alias("_Final_Species_Target"))
        else:
            merged_df = merged_df.with_columns(pl.lit("Unknown").alias("_Final_Species_Target"))

        sp_thresh = self.args.sp_thresh
        df_confirmed = merged_df.filter((pl.col("Avg_Read_ANI").is_not_null()) & (pl.col("Avg_Read_ANI") >= sp_thresh) | (pl.col("Avg_Read_ANI").is_null())).with_columns(pl.col("_Final_Species_Target").alias("Adjusted_Species"))
        df_novel = merged_df.filter((pl.col("Avg_Read_ANI").is_not_null()) & (pl.col("Avg_Read_ANI") < sp_thresh)).with_columns(pl.concat_str([pl.lit("s__unclassified_"), pl.col("_Final_Species_Target").str.replace_all(" ", "_")]).alias("Adjusted_Species"))
        
        base_cols = ["Sample", "taxid", "Adjusted_Species", "Species_NCBI", "Species_ICTV", "Rep_Accession", "Rep_Length", "Rep_Coverage(%)", "Rep_MeanDepth", "Asm_EM_Reads", "Uniq_Reads", "Multi_Reads", "Unique(%)", "Avg_Read_ANI", "Avg_Pi", "Asm_CPM", "Asm_RPM", "Asm_FPKM", "Asm_TPM", "Asm_Rel_Abund(%)", "Predicted_Support", "Poisson_Ratio", "Segment_Accessions", "Rep_Reads"]
        if 'gene_total_cov' in merged_df.columns:
            base_cols.extend(['gene_total_cov', 'gene_avr_cov'])
        final_cols = [c for c in base_cols if c in df_confirmed.columns]
        
        if len(df_novel) > 0: df_novel.select(final_cols).write_csv(str(self.output_dir / "summary" / "all_viruses.unclassified.tsv"), separator='\t')
        
        if len(df_confirmed) > 0:
            best_csv = self.output_dir / "summary" / "all_viruses.best.summary.tsv"
            df_confirmed.select(final_cols).write_csv(str(best_csv), separator='\t')
            plots_dir = self.output_dir / 'plots'
            with open(plots_dir / 'virus_frequency_plot.R', 'w', encoding='utf-8') as f: f.write(R_PLOT_SCRIPT)
            if shutil.which('Rscript'): subprocess.run(['Rscript', str(plots_dir / 'virus_frequency_plot.R'), '-i', str(best_csv), '-o', str(plots_dir / 'virus_analysis'), '--log10-transform'], capture_output=True)
            
        self.generate_report_txt(len(df), len(final_df), len(df_confirmed), len(df_novel), df_confirmed)

    def _export_sample_resource_report(self):
        """从 sample_*.log 提取每个样本的 Time/Mem/CPU，导出 TSV"""
        import csv, re as re_mod
        out_path = self.logs_dir / "sample_resource_usage.tsv"
        rows = []
        for fn in sorted(self.logs_dir.glob("sample_*.log")):
            sname = fn.name[7:-4]
            t_sec, mem_mb, cpu = 0.0, 0.0, 0.0
            with open(fn) as f:
                for line in f:
                    if '耗时:' in line and '峰值内存' in line:
                        m = re_mod.search(r'耗时:\s*(.+?)\s*\|.*?峰值内存.*?:\s*([\d.]+).*?CPU利用率:\s*([\d.]+)', line)
                        if m:
                            t_str, mem_mb, cpu = m.group(1), float(m.group(2)), float(m.group(3))
                            t_sec = 0.0
                            for part in t_str.strip().split():
                                if 'h' in part: t_sec += float(part.replace('h',''))*3600
                                elif 'm' in part: t_sec += float(part.replace('m',''))*60
                                elif 's' in part: t_sec += float(part.replace('s',''))
            rows.append([sname, self.args.tool.upper(), "Success", f"{t_sec:.2f}", f"{mem_mb:.2f}", f"{cpu:.0f}"])
        with open(out_path, 'w', newline='') as f:
            w = csv.writer(f, delimiter='\t')
            w.writerow(["Sample","Tool","Status","Time(s)","Peak_Memory(MB)","CPU(%)"])
            w.writerows(rows)
        self.logger.info(f"📊 样本资源报表已导出: {out_path} ({len(rows)} samples)")

    def generate_report_txt(self, raw_c, sum_c, best_c, unc_c, best_df):
        with open(self.output_dir / 'summary' / 'analysis_report.txt', 'w', encoding='utf-8') as f:
            f.write("=" * 80 + "\n宏病毒定量全景报告 (V42.4 Edition)\n" + "=" * 80 + "\n\n")
            f.write(f"  核心引擎: {self.args.tool.upper()} ({'极速伪比对' if self.is_pseudo else '传统比对'})\n")
            f.write(f"  处理样本: {len(self.samples)} | 留存(Raw): {raw_c} | 达标(Summary): {sum_c}\n")
            f.write(f"  确诊白名单(Best): {best_c} | 降级新种(Unclassified): {unc_c}\n")
            if not self.is_pseudo: f.write(f"  质控策略: {'CoverM 严格清洗' if self.args.use_coverm else '关闭 CoverM (基础提取)'}\n")
            f.write(f"  打假策略: 启动 Poisson Ratio 建模, 要求覆盖比值 >= {self.args.ratio}\n")
            f.write(f"  过滤模式: {'双轨制 (全基因组泊松打假 + 活跃转录区) 启用' if self.args.genes_cov else '单轨全基因组泊松过滤'}\n\n")
            
            if len(best_df) > 0:
                f.write("【确诊群落概览】\n" + "-" * 60 + "\n")
                asm_counts = best_df.get_column('Adjusted_Species').value_counts()
                for row in asm_counts.iter_rows(named=True):
                    asm = row['Adjusted_Species']
                    vd = best_df.filter(pl.col('Adjusted_Species') == asm)
                    count = row.get('count', row.get('counts', len(vd)))
                    f.write(f"🎯 {asm}: 检出 {count} 例 (群体检出率 {(count/len(self.samples))*100:.1f}%)\n")
                    f.write(f"   ├─ 平均 CPM: {vd['Asm_CPM'].mean():.2f} | 平均 FPKM: {vd['Asm_FPKM'].mean():.2f}\n")
                    if not self.is_pseudo: f.write(f"   ├─ 测定 ANI: {vd['Avg_Read_ANI'].mean() or 0:.2f}% | Pi: {vd['Avg_Pi'].mean() or 0:.5f}\n")
                    f.write(f"   └─ 覆盖度: {vd['Rep_Coverage(%)'].mean():.2f}% | 泊松打假得分(Ratio): {vd['Poisson_Ratio'].mean():.2f}\n\n")

    def run_pipeline(self):
        pipeline_start_time = time.time()
        self.logger.info("=" * 60)
        self.logger.info("🚀 V42.4 Smart-Index & Resource Monitor Pipeline 启动")
        self.logger.info("=" * 60)
        
        self.setup_output_directory(); self.build_index(); self.find_samples()
        batches = [self.samples[i:i + self.args.batch_size] for i in range(0, len(self.samples), self.args.batch_size)]
        
        for i, batch in enumerate(batches):
            batch_idx = i + 1; batch_file = self.output_dir / 'batches' / f'batch_{batch_idx:04d}.parquet'
            if self.args.resume and batch_file.exists(): 
                self.logger.info(f"⏭️ 批次断点触发: [批次 {batch_idx}/{len(batches)}] 已完成，整体跳过..."); continue
            
            self.logger.info(f"📦 开始计算 [批次 {batch_idx}/{len(batches)}] (内含 {len(batch)} 个样本)...")
            res_list = []
            ctx = mp.get_context('spawn')
            with ProcessPoolExecutor(max_workers=self.args.threads, mp_context=ctx) as executor:
                futures = {executor.submit(self.process_sample, s): s for s in batch}
                if tqdm:
                    with tqdm(total=len(batch), desc=f"批次 {batch_idx}", dynamic_ncols=True, colour="green") as pbar:
                        for future in as_completed(futures):
                            sample_info = futures[future]; pbar.set_postfix_str(f"最新完成: {sample_info['name']}")
                            if (result := future.result()): res_list.extend(result)
                            pbar.update(1)
                else:
                    for future in as_completed(futures):
                        if (result := future.result()): res_list.extend(result)
            pl.DataFrame(res_list, schema=METRICS_SCHEMA).write_parquet(batch_file) if res_list else pl.DataFrame(schema=METRICS_SCHEMA).write_parquet(batch_file)
            
        self.logger.info("🧩 所有批次处理完成，正在合并数据并执行分类学与全局丰度重整...")
        all_batch_files = list((self.output_dir / 'batches').glob('batch_*.parquet'))
        if all_batch_files: self.summarize_results_polars(pl.concat([pl.read_parquet(f) for f in all_batch_files]))
        
        # 🟢 导出每个样本的资源消耗 TSV
        self._export_sample_resource_report()
        
        pipeline_elapsed = time.time() - pipeline_start_time
        peak_mem = get_memory_mb()
        
        # 🟢 [新增] 最终输出资源总报表
        self.logger.info("=" * 60)
        self.logger.info("🎉 核心管线全部运算结束！")
        self.logger.info(f"⏱️ 整体总耗时: {format_time(pipeline_elapsed)}")
        self.logger.info(f"💾 全局峰值内存: {peak_mem:.2f} MB")
        self.logger.info("=" * 60)
        
        with open(self.logs_dir / 'pipeline_hardware_report.txt', 'w', encoding='utf-8') as f:
            f.write("【FastViromeExplorer Pro 运行资源监控报表】\n")
            f.write(f"时间: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write("-" * 40 + "\n")
            f.write(f"处理样本总数: {len(self.samples)}\n")
            f.write(f"并行策略: {self.args.threads} 并发样本 / {self.args.align_threads} 比对线程\n")
            f.write("-" * 40 + "\n")
            f.write(f"整体总耗时: {format_time(pipeline_elapsed)}\n")
            f.write(f"全局峰值内存: {peak_mem:.2f} MB\n")

def main():
    try:
        mp.set_start_method('spawn')
    except RuntimeError:
        pass

    parser = argparse.ArgumentParser(description='宏病毒鉴定与精确定量统一管线 (V42.4 Resource Monitor Edition)')
    io_group = parser.add_argument_group('📥 输入输出')
    io_group.add_argument('-i', '--input_dir', required=False, help='输入 FASTQ/FASTA 文件夹')
    io_group.add_argument('--sample_list', type=str, help='指定样本列表文件 (替代 -i)')
    io_group.add_argument('-o', '--output_dir', default='./virus_out', help='输出目录')
    io_group.add_argument('--single_end', action='store_true', help='强制单端模式')
    db_group = parser.add_argument_group('🗄️ 数据库参数')
    db_group.add_argument('-r', '--reference', required=True, help='全局参考基因组 FASTA')
    db_group.add_argument('--ref_info', type=str, required=True, help='本地参考信息 TSV 文件')
    db_group.add_argument('--taxid_clusters', type=str, help='同义 TaxID 映射文件')
    perf_group = parser.add_argument_group('🚀 核心引擎与并发')
    perf_group.add_argument('--tool', choices=['kallisto', 'salmon', 'bowtie2', 'bwa', 'minimap2', 'strobealign', 'bwa-mem2', 'hisat2'], default='bowtie2', help='比对工具 (伪比对/传统比对)')
    perf_group.add_argument('-t', '--threads', type=int, default=8, help='并发样本数')
    perf_group.add_argument('--align_threads', type=int, default=4, help='单样本内部线程数')
    perf_group.add_argument('--batch_size', type=int, default=20, help='批次刷盘保护数量')
    filt_group = parser.add_argument_group('🛡️ 过滤控制')
    filt_group.add_argument('--coverage', type=float, default=10.0, help='A轨: 绝对全长覆盖度下限(%%) (默认10.0)')
    filt_group.add_argument('--ratio', type=float, default=0.3, help='A轨: 泊松分布覆盖度比值下限 (默认 0.3)')
    filt_group.add_argument('--meandepth', type=float, default=0.0, help='A轨: 代表序列最小平均深度 (可选)')
    filt_group.add_argument('--min_tpm', type=float, default=0.0, help='最少 TPM 值')
    filt_group.add_argument('--min_uniq_reads', type=int, default=1, help='最少独特比对读数')
    filt_group.add_argument('--sp_thresh', type=float, default=95.0, help='物种ANI识别阈值%% (仅限传统比对)')
    filt_group.add_argument('--genes_cov', type=str, help='B轨: 导入转录覆盖率文件启动【双轨制】放行策略')
    filt_group.add_argument('--min_gene_total_cov', type=float, default=80.0, help='B轨: 最低转录区总覆盖')
    filt_group.add_argument('--min_gene_avr_cov', type=float, default=5.0, help='B轨: 最低转录区平均覆盖')
    misc_group = parser.add_argument_group('⚙️ 其他设置')
    misc_group.add_argument('--resume', action='store_true', help='开启双重断点续传')
    misc_group.add_argument('--keep_tmp', action='store_true', help='保留庞大的中间 BAM 文件')
    misc_group.add_argument('--use_coverm', action='store_true', help='【仅传统比对】启用 CoverM 严格清洗机制')
    misc_group.add_argument('--min_aln_len', type=int, default=80, help='CoverM: 最小比对长度')
    misc_group.add_argument('--min_aln_prop', type=float, default=0.85, help='CoverM: 最小比对比例')
    misc_group.add_argument('--min_pid', type=float, default=0.90, help='CoverM: 最小序列相似度')
    misc_group.add_argument('--verbose', action='store_true', help='输出详细底层日志')
    
    args = parser.parse_args()
    if not args.input_dir and not args.sample_list: parser.error("必须提供 -i (输入目录) 或 --sample_list (样本列表)")
    UnifiedVirusPipeline(args).run_pipeline()

if __name__ == '__main__':
    main()
