#!/usr/bin/env python3
"""
宏病毒鉴定与精确定量管线 (Virus Quantification - V34.3 Batch-Processing Master)
================================================================================
【V34.3 工业级批次架构与早期提取】
1. [批次架构] 引入 Batch 机制，按批次落盘临时结果(.parquet)，彻底杜绝超大队列内存泄漏。
2. [强悍断点] 实现 "批次级" 与 "样本级" 双重断点续跑，完美应对集群意外中断。
3. [早期提取] 重构 --extract 逻辑：在 CoverM 过滤前基于初始 raw.bam 提取 reads，保证信息完整度。
4. [智能寻址] 完美支持同一目录下单/双端混合、孤立 R2 等复杂情况。
5. [严谨降级] 设立 0.1% 真实突变下限，严格隔离测算失败(Failed)与真实新种(Novel)。
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
from pathlib import Path
from collections import defaultdict, Counter
from concurrent.futures import ProcessPoolExecutor, as_completed
import logging
import colorlog
from tqdm import tqdm
import polars as pl
import pysam
import numpy as np

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
data$Display_Name <- paste0(data$Adjusted_Species, "\n(TaxID: ", data$taxid, ")")
data$Display_Name <- sapply(data$Display_Name, function(x) paste(strwrap(x, width = 45), collapse = "\n"))
metrics <- c("Asm_EM_Reads", "Asm_CPM", "Asm_FPKM", "Avg_Read_ANI")
available_metrics <- intersect(metrics, colnames(data))
plot_data <- data %>% select(Display_Name, Sample, all_of(available_metrics)) %>% pivot_longer(cols = all_of(available_metrics), names_to = "Metric", values_to = "Value") %>% filter(!is.na(Value))
if (opt$`log10-transform`) {
  for (m in unique(plot_data$Metric)) { min_v <- min(plot_data$Value[plot_data$Metric==m & plot_data$Value>0], na.rm=TRUE); plot_data$Value[plot_data$Metric==m & plot_data$Value<=0] <- min_v/10 }
}
medians <- aggregate(Value ~ Display_Name, data=plot_data[plot_data$Metric=="Asm_EM_Reads",], median)
plot_data$Display_Name <- factor(plot_data$Display_Name, levels=medians[order(medians$Value), "Display_Name"])
p <- ggplot(plot_data, aes(x=Display_Name, y=Value)) + geom_boxplot(aes(fill=Display_Name), alpha=0.6, outlier.shape=NA) + geom_point(aes(color=Display_Name), position=position_jitter(width=0.2, height=0), alpha=0.6) + facet_wrap(~ Metric, scales="free_x", ncol=length(available_metrics)) + scale_fill_viridis_d(option="turbo") + scale_color_viridis_d(option="turbo") + theme_bw(base_size=13) + theme(legend.position="none", axis.text.y=element_text(size=9, face="italic")) + coord_flip()
if (opt$`log10-transform`) p <- p + scale_y_log10()
ggsave(sprintf("%s_multi_metrics.pdf", opt$output), plot=p, width=opt$width*1.5, height=max(opt$height, length(unique(data$Display_Name))*0.8), dpi=300)
"""

class TqdmLoggingHandler(logging.Handler):
    def __init__(self, level=logging.NOTSET):
        super().__init__(level)
        self.formatter = colorlog.ColoredFormatter(
            '%(log_color)s[%(asctime)s] %(levelname)s - %(message)s', 
            datefmt='%H:%M:%S', 
            log_colors={'DEBUG': 'cyan', 'INFO': 'green', 'WARNING': 'yellow', 'ERROR': 'red', 'CRITICAL': 'bold_red'}
        )
    def emit(self, record):
        try:
            msg = self.format(record)
            tqdm.write(msg)
            self.flush()
        except Exception:
            self.handleError(record)

def setup_logging(verbose=False):
    logger = logging.getLogger("QuantPipe_V34")
    logger.setLevel(logging.DEBUG if verbose else logging.INFO)
    handler = TqdmLoggingHandler()
    logger.addHandler(handler)
    return logger

logger = setup_logging()

class Timer:
    def __init__(self, name=""): self.name = name
    def __enter__(self): self.start_time = time.time(); return self
    def __exit__(self, exc_type, exc_val, exc_tb): logger.info(f"{'✅ 完成' if exc_type is None else '❌ 失败'}: {self.name} [耗时: {time.time() - self.start_time:.1f}秒]")

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

# Polars 数据结构 Schema (用于强制对齐批次结果类型)
METRICS_SCHEMA = {
    'Sample': pl.String, 'Accession': pl.String, 'Length': pl.Int64,
    'taxid': pl.String, 'Species': pl.String, 'Segment': pl.String,
    'Coverage(%)': pl.Float64, 'MeanDepth': pl.Float64,
    'EM_Reads': pl.Float64, 'Uniq_Reads': pl.Int64, 'Multi_Reads': pl.Int64,
    'Sample_Total_Mapped': pl.Int64, 'Sample_Global_Reads': pl.Int64
}

class VirusQuantificationPipeline:
    def __init__(self, args):
        self.args = args
        self.samples = []
        self.output_dir = Path(args.output_dir)
        self.index_path = None
        self.ref_length_dict = {}
        self.taxid_clusters = {}  
        self.tax_map = {}         
        
        global logger
        logger.setLevel(logging.DEBUG if self.args.verbose else logging.INFO)
        
        self.check_tools()
        self._load_reference_lengths()
        self._load_taxid_clusters()
        self._load_ref_info()
        
    def write_sample_log(self, sample_name, message, level="info"):
        log_file = self.output_dir / 'logs' / f"{sample_name}.pipeline.log"
        with open(log_file, 'a', encoding='utf-8') as f:
            f.write(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {message}\n")
        
        if self.args.verbose or level == "error":
            if level == "info": logger.info(f"[{sample_name}] {message}")
            elif level == "warning": logger.warning(f"[{sample_name}] {message}")
            elif level == "error": logger.error(f"[{sample_name}] {message}")

    def check_tools(self):
        required_tools = {'samtools': 'BAM处理', 'coverm': '序列清洗', 'pandepth': '深度计算'}
        tool_map = {
            'bowtie2': 'bowtie2-build', 'bwa': 'bwa', 'bwa-mem2': 'bwa-mem2',
            'hisat2': 'hisat2-build', 'minimap2': 'minimap2', 'strobealign': 'strobealign'
        }
        if self.args.tool in tool_map:
            required_tools[tool_map[self.args.tool]] = f'{self.args.tool} 建库/比对工具'
        
        if self.args.extract:
            required_tools['pigz'] = '多线程压缩工具'
            
        missing = [f"{t} ({d})" for t, d in required_tools.items() if shutil.which(t) is None]
        if missing: 
            logger.error("❌ 缺少必要工具:"); [logger.error(f"  - {m}") for m in missing]
            sys.exit(1)

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
        if self.args.taxid_clusters and os.path.exists(self.args.taxid_clusters):
            with open(self.args.taxid_clusters, 'r') as f:
                for line in f:
                    if line.startswith('object_taxid'): continue
                    parts = line.strip().split('\t')
                    if len(parts) >= 2: self.taxid_clusters[parts[0]] = parts[1]

    def _load_ref_info(self):
        if not self.args.ref_info or not os.path.exists(self.args.ref_info):
            logger.error("❌ 必须提供 --ref_info (包含 Accession, taxid, Species, Segment 等列的本地 TSV)")
            sys.exit(1)
        with open(self.args.ref_info, 'r', encoding='utf-8') as f:
            header = None; idx_acc, idx_tax, idx_sp, idx_seg = 0, 3, 4, 13
            for line in f:
                if line.startswith('####'): continue
                parts = line.rstrip('\n').split('\t')
                if header is None and 'Accession' in line:
                    header = parts
                    if 'Accession' in header: idx_acc = header.index('Accession')
                    if 'taxid' in header: idx_tax = header.index('taxid')
                    if 'Species' in header: idx_sp = header.index('Species')
                    if 'Segment' in header: idx_seg = header.index('Segment')
                    continue
                if len(parts) <= max(idx_acc, idx_tax): continue
                acc = parts[idx_acc].strip()
                raw_taxid = parts[idx_tax].strip() if len(parts) > idx_tax else "Unannotated"
                true_taxid = self.taxid_clusters.get(raw_taxid, raw_taxid)
                
                info = {
                    'taxid': true_taxid,
                    'species': parts[idx_sp].strip() if len(parts) > idx_sp else "Unannotated",
                    'segment': parts[idx_seg].strip() if len(parts) > idx_seg else ""
                }
                self.tax_map[acc] = info; self.tax_map[acc.split('.')[0]] = info

    def setup_output_directory(self):
        subdirs = ['bam', 'stat', 'summary', 'logs', 'index', 'plots', 'batches', 'reads', 'fasta']
        for subdir in subdirs: (self.output_dir / subdir).mkdir(parents=True, exist_ok=True)

    def build_index(self):
        tool = self.args.tool
        ref_path = Path(self.args.reference).resolve()
        
        if tool in ['strobealign', 'minimap2']:
            self.index_path = str(ref_path)
            return
            
        prefix = self.output_dir / 'index' / f"{ref_path.stem}_{tool}"
        self.index_path = str(prefix)
        cmd = []
        
        if tool == 'bwa':
            if Path(f"{prefix}.bwt").exists(): return
            cmd = ['bwa', 'index', '-p', str(prefix), str(ref_path)]
        elif tool == 'bwa-mem2':
            if Path(f"{prefix}.bwt.2bit.64").exists() or Path(f"{prefix}.0123").exists(): return
            cmd = ['bwa-mem2', 'index', '-p', str(prefix), str(ref_path)]
        elif tool == 'bowtie2':
            if Path(f"{prefix}.1.bt2").exists() or Path(f"{prefix}.1.bt2l").exists(): return
            cmd = ['bowtie2-build', '--threads', str(self.args.threads), str(ref_path), str(prefix)]
        elif tool == 'hisat2':
            if Path(f"{prefix}.1.ht2").exists() or Path(f"{prefix}.1.ht2l").exists(): return
            cmd = ['hisat2-build', '-p', str(self.args.threads), str(ref_path), str(prefix)]
            
        logger.info(f"🏗️ 构建全局 {tool} 索引...")
        with Timer(f"构建 {tool} 索引"):
            result = subprocess.run(cmd, capture_output=True, text=True)
            if result.returncode != 0:
                logger.error(f"❌ 索引构建失败: {result.stderr}")
                sys.exit(1)

    def find_samples(self):
        """【深度重构】强大的单/双端混合寻址逻辑"""
        input_dir = Path(self.args.input_dir) if self.args.input_dir else None
        valid_exts = ['.fastq.gz', '.fq.gz', '.fasta.gz', '.fa.gz', '.fastq', '.fq', '.fasta', '.fa']
        samples = []
        
        if input_dir:
            file_dict = defaultdict(dict)
            for f in input_dir.iterdir():
                if f.is_dir(): continue
                name = f.name
                
                ext_found = None
                for ext in valid_exts:
                    if name.lower().endswith(ext):
                        ext_found = ext; break
                if not ext_found: continue
                
                base_name = name[:-len(ext_found)]
                m = re.search(r'([._\-]?)(R?[12])(_\d+)?$', base_name, re.IGNORECASE)
                
                if m and m.group(2):
                    direction = 'R1' if '1' in m.group(2) else 'R2'
                    sname = base_name[:m.start()]
                    if not sname: sname = base_name
                else:
                    direction = 'SE'; sname = base_name
                
                sname_clean = re.sub(r'(_unmapped|\.unmapped|_trimmed|\.trimmed)$', '', sname, flags=re.IGNORECASE)
                sname = sname_clean if sname_clean else sname
                file_dict[sname][direction] = str(f)

            for sname, reads in file_dict.items():
                if self.args.single_end:
                    for d, path in reads.items():
                        sn = sname if len(reads) == 1 else f"{sname}_{d}"
                        samples.append({'name': sn, 'r1': path, 'r2': None})
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
        
        unique_samples, seen = [], set()
        for s in samples:
            if s['name'] not in seen:
                unique_samples.append(s); seen.add(s['name'])
        
        self.samples = unique_samples
        if not self.samples:
            logger.error("❌ 未找到任何样本文件，请检查目录或单双端命名规则。")
            sys.exit(1)
            
        pe_count = sum(1 for s in self.samples if s['r2'])
        logger.info(f"✅ 成功寻址 {len(self.samples)} 个测序样本 (自动识别出: {pe_count} PE, {len(self.samples)-pe_count} SE)")

    def _extract_reads_from_bam(self, bam_file, sname, is_fasta, is_paired):
        """【新增】在CoverM过滤前，基于原版 BAM 提取 Reads"""
        out_dir = self.output_dir / 'reads'
        ext_name = "fasta.gz" if is_fasta else "fastq.gz"
        cmd_tool = "fasta" if is_fasta else "fastq"

        r1_out = out_dir / f"{sname}_virus_1.{ext_name}"
        r2_out = out_dir / f"{sname}_virus_2.{ext_name}"
        rs_out = out_dir / f"{sname}_virus_single.{ext_name}"

        # 样本级断点续跑：提取过就不再提取
        if self.args.resume and (r1_out.exists() or rs_out.exists()):
            return True

        if not Path(bam_file).exists(): return False

        self.write_sample_log(sname, f"-> [早期提取] 正从 {Path(bam_file).name} 中无损提取原始 Mapped Reads...", level="info")
        try:
            if is_paired:
                cmd = (
                    f"samtools {cmd_tool} -F 4 "
                    f"-1 >(pigz -p 2 -c > '{r1_out}') "
                    f"-2 >(pigz -p 2 -c > '{r2_out}') "
                    f"-0 /dev/null -s /dev/null '{bam_file}'"
                )
            else:
                cmd = f"samtools {cmd_tool} -F 4 '{bam_file}' | pigz -p 2 -c > '{rs_out}'"

            subprocess.run(cmd, shell=True, executable='/bin/bash', stderr=subprocess.PIPE, check=True)
            
            for f in [r1_out, r2_out, rs_out]:
                if f.exists() and f.stat().st_size <= 50: f.unlink()
            return True
        except Exception as e:
            self.write_sample_log(sname, f"Reads 提取异常: {e}", level="warning")
            return False

    def get_global_bam_qc(self, bam_file, sample_name):
        qc_file = self.output_dir / 'stat' / f"{sample_name}.global_qc.txt"
        if self.args.resume and qc_file.exists() and qc_file.stat().st_size > 0:
            return qc_file
            
        cmd = f"samtools stats '{bam_file}' | grep ^SN | cut -f 2- > '{qc_file}'"
        try:
            subprocess.run(cmd, shell=True, check=True, executable='/bin/bash')
            with open(qc_file, 'r') as f:
                lines = f.readlines()
                total, mapped, mq0, unmapped = 0, 0, 0, 0
                for line in lines:
                    if 'raw total sequences' in line: total = line.split('\t')[1].strip()
                    elif line.startswith('reads mapped:'): mapped = line.split('\t')[1].strip()
                    elif line.startswith('reads unmapped:'): unmapped = line.split('\t')[1].strip()
                    elif line.startswith('reads MQ0:'): mq0 = line.split('\t')[1].strip()
                
                log_msg = f"-> [质控报告] 总序列: {total} | 成功比对: {mapped} | 未比对: {unmapped} | 多重比对(MQ0): {mq0}"
                self.write_sample_log(sample_name, log_msg, level="info")
            return qc_file
        except Exception as e:
            self.write_sample_log(sample_name, f"❌ BAM 统计失败: {e}", level="warning")
            return None

    def align_and_coverm(self, sample):
        sname = sample['name']
        raw_bam = self.output_dir / 'bam' / f'{sname}.raw.bam'
        filt_bam = self.output_dir / 'bam' / f'{sname}.sorted.bam'
        bai_path = str(filt_bam) + '.bai'
        qc_file = self.output_dir / 'stat' / f"{sname}.global_qc.txt"
        
        is_fasta = any(sample['r1'].lower().endswith(ext) for ext in ['.fa', '.fasta', '.fa.gz', '.fasta.gz'])
        is_paired = bool(sample['r2'])
        
        # 样本级断点续跑：如果 sorted.bam 存在，直接跳过计算
        if self.args.resume and filt_bam.exists() and Path(bai_path).exists() and filt_bam.stat().st_size > 1024:
            self.write_sample_log(sname, f"⏭️ 发现已完成的 BAM，跳过比对清洗步骤", level="info")
            mapped_res = subprocess.run(['samtools', 'view', '-c', str(filt_bam)], capture_output=True, text=True)
            total_mapped = int(mapped_res.stdout.strip() or 1)
            
            global_reads = total_mapped
            if qc_file.exists():
                with open(qc_file, 'r') as f:
                    for line in f:
                        if 'raw total sequences' in line:
                            global_reads = int(line.split('\t')[1].strip())
            
            # 补齐提取逻辑：如果之前由于意外中断未提取，现在补提
            if self.args.extract:
                self._extract_reads_from_bam(raw_bam if raw_bam.exists() else filt_bam, sname, is_fasta, is_paired)
                
            return filt_bam, total_mapped, global_reads

        self.write_sample_log(sname, f"▶️ [步骤 1/4] 开始序列比对 ({self.args.tool})", level="info")
        
        tool = self.args.tool
        align_cmd = []
        inner_threads = min(2, self.args.align_threads) 
        
        if tool == 'strobealign':
            align_cmd = ['strobealign', '-t', str(inner_threads), self.index_path, sample['r1']] + ([sample['r2']] if sample['r2'] else [])
        elif tool == 'minimap2':
            preset = 'sr' if sample['r2'] else 'map-ont'
            align_cmd = ['minimap2', '-ax', preset, '-t', str(inner_threads), self.index_path, sample['r1']] + ([sample['r2']] if sample['r2'] else [])
        elif tool == 'bwa':
            align_cmd = ['bwa', 'mem', '-v', '1', '-t', str(inner_threads), self.index_path, sample['r1']] + ([sample['r2']] if sample['r2'] else [])
        elif tool == 'bwa-mem2':
            align_cmd = ['bwa-mem2', 'mem', '-v', '1', '-t', str(inner_threads), self.index_path, sample['r1']] + ([sample['r2']] if sample['r2'] else [])
        elif tool == 'bowtie2':
            align_cmd = ['bowtie2', '-p', str(inner_threads), '-x', self.index_path]
            if is_fasta: align_cmd.append('-f')
            if sample['r2']: align_cmd.extend(['-1', sample['r1'], '-2', sample['r2']])
            else: align_cmd.extend(['-U', sample['r1']])
        elif tool == 'hisat2':
            align_cmd = ['hisat2', '-p', str(inner_threads), '-x', self.index_path]
            if is_fasta: align_cmd.append('-f')
            if sample['r2']: align_cmd.extend(['-1', sample['r1'], '-2', sample['r2']])
            else: align_cmd.extend(['-U', sample['r1']])

        cmd_str = ' '.join(shlex.quote(str(x)) for x in align_cmd)
        raw_bam_str = shlex.quote(str(raw_bam))
        
        align_proc = subprocess.run(
            f"set -o pipefail; {cmd_str} | samtools view -b -o {raw_bam_str}", 
            shell=True, executable='/bin/bash', capture_output=True, text=True
        )
        if align_proc.returncode != 0:
            err_msg = align_proc.stderr[-500:] if align_proc.stderr else "Unknown Fatal Error"
            self.write_sample_log(sname, f"❌ 比对器抛出异常: {err_msg}", level="error")
            raise RuntimeError(f"Alignment crashed for {sname}")
        
        self.get_global_bam_qc(raw_bam, sname)
        total_raw_res = subprocess.run(['samtools', 'view', '-c', str(raw_bam)], capture_output=True, text=True)
        global_reads = int(total_raw_res.stdout.strip() or 1)
        
        mapped_res = subprocess.run(['samtools', 'view', '-c', '-F', '4', str(raw_bam)], capture_output=True, text=True)
        total_mapped = int(mapped_res.stdout.strip() or 1)

        # 【重点】在这里执行早期 Reads 提取（CoverM 过滤前），保留最大信息量
        if self.args.extract:
            self._extract_reads_from_bam(raw_bam, sname, is_fasta, is_paired)

        self.write_sample_log(sname, f"▶️ [步骤 2/4] 开始 CoverM 质量清洗与过滤", level="info")
        
        filter_cmd = [
            'coverm', 'filter', '-b', str(raw_bam), '-o', str(filt_bam) + '.unsorted',
            '--min-read-aligned-length', str(self.args.min_aln_len), 
            '--min-read-percent-identity', str(self.args.min_pid),
            '--min-read-aligned-percent', str(self.args.min_aln_prop), 
            '--include-secondary', '-t', '1'
        ]
        
        if is_paired:
            filter_cmd.append('--proper-pairs-only')
        
        filter_res = subprocess.run(filter_cmd, capture_output=True, text=True)
        if filter_res.returncode != 0: 
            self.write_sample_log(sname, f"CoverM 警告: {filter_res.stderr[:200]}", level="warning")
            
        subprocess.run(['samtools', 'sort', '-@', '1', '-o', str(filt_bam), str(filt_bam) + '.unsorted'])
        
        if not os.path.exists(bai_path):
            try: pysam.index(str(filt_bam))
            except Exception: pass
        
        if raw_bam.exists(): raw_bam.unlink()
        unsorted_bam = Path(str(filt_bam) + '.unsorted')
        if unsorted_bam.exists(): unsorted_bam.unlink()
        
        return filt_bam, total_mapped, global_reads

    def run_em_taxid_allocation(self, filt_bam_path, sname):
        self.write_sample_log(sname, f"▶️ [步骤 3/4] 开始 EM 多重比对丰度分配", level="info")
        
        read_best_data = {}
        with pysam.AlignmentFile(filt_bam_path, "rb") as bam:
            for read in bam.fetch(until_eof=True):
                if read.is_unmapped: continue
                ref = bam.get_reference_name(read.reference_id)
                taxid = self.tax_map.get(ref, {}).get('taxid', ref)
                try: score = read.get_tag('AS')
                except KeyError: score = read.query_alignment_length - read.get_tag('NM', 0)
                
                qname = read.query_name
                if qname not in read_best_data:
                    read_best_data[qname] = {'max_score': score, 'refs': {ref: score}, 'taxids': {taxid}}
                else:
                    curr_max = read_best_data[qname]['max_score']
                    if score > curr_max:
                        read_best_data[qname] = {'max_score': score, 'refs': {ref: score}, 'taxids': {taxid}}
                    elif score == curr_max:
                        read_best_data[qname]['refs'][ref] = score
                        read_best_data[qname]['taxids'].add(taxid)

        uniq_taxid_counts = defaultdict(int)
        for data in read_best_data.values():
            if len(data['taxids']) == 1: uniq_taxid_counts[list(data['taxids'])[0]] += 1

        em_counts, uniq_counts, multi_counts = defaultdict(float), defaultdict(int), defaultdict(int)
        for data in read_best_data.values():
            taxids = list(data['taxids'])
            refs_dict = data['refs']
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

        return em_counts, uniq_counts, multi_counts

    def process_sample(self, sample):
        sname = sample['name']
        try:
            sorted_bam, total_mapped, global_reads = self.align_and_coverm(sample)
            if total_mapped <= 1: 
                self.write_sample_log(sname, "-> 过滤后 mapped reads ≤ 1，提前终止计算", level="info")
                return []
            
            em_counts, uniq_counts, multi_counts = self.run_em_taxid_allocation(sorted_bam, sname)
            
            self.write_sample_log(sname, f"▶️ [步骤 4/4] 开始 Pandepth 深度与覆盖度测算", level="info")
            stat_prefix = self.output_dir / 'stat' / sname
            stat_file = str(stat_prefix) + '.chr.stat.gz'
            
            pan_res = subprocess.run(['pandepth', '-a', '-i', str(sorted_bam), '-o', str(stat_prefix), '-t', '1'], capture_output=True, text=True)
            depth_dict = {}
            if pan_res.returncode != 0:
                self.write_sample_log(sname, f"Pandepth 失败/警告: {pan_res.stderr[:200]}", level="warning")
                
            if os.path.exists(stat_file):
                try:
                    with gzip.open(stat_file, 'rt') as f:
                        for line in f:
                            if not line.startswith('#'):
                                p = line.split()
                                if len(p) >= 6: depth_dict[p[0]] = {'Cov': float(p[4]), 'Dep': float(p[5])}
                except Exception as e: pass
                        
            results = []
            for virus, em_mapped in em_counts.items():
                if em_mapped <= 0: continue
                u_reads, m_reads = uniq_counts.get(virus, 0), multi_counts.get(virus, 0)
                length_val = self.ref_length_dict.get(virus, 1)
                d_info = depth_dict.get(virus, {'Cov': 0, 'Dep': 0})
                    
                results.append({
                    'Sample': sname, 'Accession': virus, 'Length': length_val,
                    'taxid': self.tax_map.get(virus, {}).get('taxid', 'Unannotated'),
                    'Species': self.tax_map.get(virus, {}).get('species', 'Unannotated'),
                    'Segment': self.tax_map.get(virus, {}).get('segment', ''),
                    'Coverage(%)': float(d_info['Cov']), 'MeanDepth': float(d_info['Dep']),
                    'EM_Reads': float(em_mapped), 'Uniq_Reads': int(u_reads), 'Multi_Reads': int(m_reads),
                    'Sample_Total_Mapped': int(total_mapped), 'Sample_Global_Reads': int(global_reads)
                })
            self.write_sample_log(sname, "✅ 定量全流程计算完毕！", level="info")
            return results
        except Exception as e:
            self.write_sample_log(sname, f"❌ 样本运行失败: {e}", level="error")
            return []

    def _compute_ani_pi_worker(self, task):
        import pysam
        sname, bam_path, ref = task
        if not os.path.exists(bam_path): 
            return {'Sample': sname, 'Rep_Accession': ref, 'Avg_Read_ANI': None, 'Avg_Pi': None}
        
        ani_sum, ani_cnt = 0.0, 0
        base_counts = defaultdict(Counter)
        
        try:
            with pysam.AlignmentFile(bam_path, "rb") as bam:
                if ref not in bam.references: 
                    return {'Sample': sname, 'Rep_Accession': ref, 'Avg_Read_ANI': None, 'Avg_Pi': None}
                
                total_reads = _get_ref_read_count_safe(bam, ref)
                if total_reads == 0: 
                    return {'Sample': sname, 'Rep_Accession': ref, 'Avg_Read_ANI': None, 'Avg_Pi': None}
                
                step = max(1, total_reads // 10000)
                
                for idx, read in enumerate(bam.fetch(ref)):
                    if idx % step != 0: continue
                    if read.is_unmapped: continue
                    
                    aln_len = read.query_alignment_length
                    if aln_len > 0:
                        try: nm = read.get_tag('NM')
                        except KeyError: nm = 0
                        ani_sum += (aln_len - nm) / aln_len
                        ani_cnt += 1
                        
                    for qpos, rpos, ref_base in read.get_aligned_pairs(matches_only=True, with_seq=True):
                        if rpos is not None and qpos is not None:
                            base_counts[rpos][read.query_sequence[qpos].upper()] += 1
                    
                    if ani_cnt >= 10000: break
                        
        except Exception: pass
            
        avg_ani = (ani_sum / ani_cnt) * 100.0 if ani_cnt > 0 else None
        pi_sum, covered_pos = 0.0, 0
        for rpos, counts in base_counts.items():
            total = sum(counts.values())
            if total > 1: 
                pi = 1.0 - sum((c/total)**2 for c in counts.values())
                pi_sum += pi
                covered_pos += 1
                
        avg_pi = (pi_sum / covered_pos) if covered_pos > 0 else None
        return {
            'Sample': sname, 
            'Rep_Accession': ref, 
            'Avg_Read_ANI': round(avg_ani, 2) if avg_ani is not None else None, 
            'Avg_Pi': round(avg_pi, 5) if avg_pi is not None else None
        }

    def extract_fastas_with_python(self, target_ids, out_dir):
        """只负责提取参考病毒的 FASTA 基因组序列 (Reads 提取已经在早期批次处理过了)"""
        target_set = set(target_ids)
        found = set()
        current_id, current_seq = None, []
        
        def save_record(vid, seq_list):
            if vid in target_set:
                vname = re.sub(r'[\\/*?:"<>| ]', "_", vid)
                tax = self.tax_map.get(vid, {}).get('species', 'Unannotated')
                safe_tax = re.sub(r'[\\/*?:"<>| ]', "_", tax)[:50]
                
                v_dir = out_dir / f"{safe_tax}_{vname}"
                v_dir.mkdir(parents=True, exist_ok=True)
                with open(v_dir / f"{safe_tax}_{vname}.ref.fasta", 'w') as outf:
                    outf.write(f">{vid}\n" + "".join(seq_list) + "\n")
                found.add(vid)
                
        with _smart_open(self.args.reference, 'rt') as f:
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    if current_id: save_record(current_id, current_seq)
                    current_id = line[1:].split()[0]
                    current_seq = []
                elif line: current_seq.append(line)
            if current_id: save_record(current_id, current_seq)

    def summarize_results_polars(self, df):
        """【批次汇总】在所有批次完成之后执行延迟合并测算"""
        if len(df) == 0:
            logger.warning("⚠️ 没有任何样本产生合格计算结果，流程退出。")
            return None
        
        ext, sep = "tsv", "\t"
        
        df = df.with_columns([
            (pl.col("EM_Reads") * 1000.0 / pl.col("Length")).alias("Seq_RPK"),
            pl.when(pl.col("Segment") != "").then(
                pl.col("Segment") + ":" + pl.col("Accession")
            ).otherwise(pl.col("Accession")).alias("Seg_Acc_Str")
        ])
        
        df.write_csv(str(self.output_dir / "summary" / f"all_viruses.raw.{ext}"), separator=sep)
        
        best_acc_df = df.sort(["Sample", "taxid", "EM_Reads", "Coverage(%)"], descending=[False, False, True, True]).group_by(["Sample", "taxid"]).first()
        rep_df = pl.DataFrame({
            "Sample": best_acc_df["Sample"], "taxid": best_acc_df["taxid"],
            "Rep_Accession": best_acc_df["Accession"], "Rep_Reads": best_acc_df["EM_Reads"]
        })
        
        sample_meta_df = df.group_by("Sample").agg([
            pl.col("Seq_RPK").sum().alias("Sample_Total_RPK"),
            pl.col("Sample_Total_Mapped").first(),
            pl.col("Sample_Global_Reads").first()
        ])
        
        asm_df = df.group_by(["Sample", "taxid"]).agg([
            pl.col("EM_Reads").sum().alias("Asm_EM_Reads"), 
            pl.col("Seq_RPK").sum().alias("Asm_RPK_Sum"),
            pl.col("Uniq_Reads").sum(), pl.col("Multi_Reads").sum(),
            pl.col("Length").sum().alias("Asm_Length"), 
            pl.col("Species").first(), 
            pl.col("Seg_Acc_Str").unique().str.join(",").alias("Segment_Accessions")
        ])
        
        asm_df = asm_df.join(rep_df, on=["Sample", "taxid"], how="left")
        asm_df = asm_df.join(sample_meta_df, on="Sample", how="left")
        
        asm_df = asm_df.with_columns([
            pl.when((pl.col("Uniq_Reads") + pl.col("Multi_Reads")) > 0)
              .then((pl.col("Uniq_Reads") / (pl.col("Uniq_Reads") + pl.col("Multi_Reads")) * 100).round(2))
              .otherwise(0.0).alias("Unique(%)"),
            
            (pl.col("Asm_EM_Reads") * 1e6 / pl.col("Sample_Total_Mapped")).round(2).alias("Asm_CPM"),
            (pl.col("Asm_EM_Reads") * 1e6 / pl.col("Sample_Global_Reads")).round(2).alias("Asm_RPM"),
            
            pl.when(pl.col("Asm_Length") > 0)
              .then((pl.col("Asm_EM_Reads") * 1e9 / (pl.col("Sample_Total_Mapped") * pl.col("Asm_Length"))))
              .otherwise(0.0).round(2).alias("Asm_FPKM"),
            
            pl.when(pl.col("Sample_Total_RPK") > 0)
              .then((pl.col("Asm_RPK_Sum") * 1e6 / pl.col("Sample_Total_RPK")).round(2))
              .otherwise(0.0).alias("Asm_TPM"),
            
            (pl.col("Asm_EM_Reads") / pl.col("Sample_Total_Mapped") * 100.0).round(4).alias("Asm_Rel_Abund(%)")
        ])
        
        rep_seq_stats = df.select(["Sample", "Accession", "Length", "Coverage(%)", "MeanDepth"]).rename({
            "Accession": "Rep_Accession", "Coverage(%)": "Rep_Coverage(%)",
            "MeanDepth": "Rep_MeanDepth", "Length": "Rep_Length"
        }).unique(subset=["Sample", "Rep_Accession"], keep="first", maintain_order=True)
        
        asm_df = asm_df.join(rep_seq_stats, on=["Sample", "Rep_Accession"], how="left")
        
        final_df = asm_df.filter(
            (pl.col("Rep_Coverage(%)") >= self.args.coverage) & 
            (pl.col("Rep_MeanDepth") >= self.args.meandepth) &   
            (pl.col("Sample_Total_Mapped") > 0) &
            (pl.col("Uniq_Reads") >= self.args.min_uniq_reads) & 
            (pl.col("Asm_TPM") >= self.args.min_tpm)
        )
        
        final_df.write_csv(str(self.output_dir / "summary" / f"all_viruses.summary.{ext}"), separator=sep)
        
        pre_best_df = final_df.sort(["Sample", "taxid", "Asm_EM_Reads"], descending=[False, False, True]).unique(subset=["Sample", "taxid"], keep="first", maintain_order=True)
        
        logger.info(f"\n🧬 正在对存活的 {len(pre_best_df)} 个核心代表株进行深度进化测算 (ANI/Pi)...")
        tasks = [(row['Sample'], str(self.output_dir / 'bam' / f"{row['Sample']}.sorted.bam"), row['Rep_Accession']) for row in pre_best_df.iter_rows(named=True)]
        
        ani_pi_results = []
        with ProcessPoolExecutor(max_workers=self.args.threads) as executor:
            futures = {executor.submit(self._compute_ani_pi_worker, task): i for i, task in enumerate(tasks)}
            with tqdm(total=len(futures), desc="[ANI/Pi进化测算]", unit="对", bar_format='{l_bar}{bar}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}]') as pbar:
                for future in as_completed(futures):
                    try:
                        if (res := future.result()): ani_pi_results.append(res)
                    except Exception: pass
                    pbar.update(1)

        if len(ani_pi_results) == 0:
            ap_df = pl.DataFrame(schema={"Sample": pl.String, "Rep_Accession": pl.String, "Avg_Read_ANI": pl.Float64, "Avg_Pi": pl.Float64})
        else:
            ap_df = pl.DataFrame(ani_pi_results)
            for col in ["Sample", "Rep_Accession", "Avg_Read_ANI", "Avg_Pi"]:
                if col not in ap_df.columns: ap_df = ap_df.with_columns(pl.lit(None, pl.Float64).alias(col))
        
        merged_df = pre_best_df.join(ap_df, on=["Sample", "Rep_Accession"], how="left")

        sp_thresh = self.args.sp_thresh
        df_confirmed = merged_df.filter((pl.col("Avg_Read_ANI").is_not_null()) & (pl.col("Avg_Read_ANI") >= sp_thresh)).with_columns(pl.col("Species").alias("Adjusted_Species"))
        df_novel = merged_df.filter((pl.col("Avg_Read_ANI").is_not_null()) & (pl.col("Avg_Read_ANI") > 0.1) & (pl.col("Avg_Read_ANI") < sp_thresh)).with_columns(pl.concat_str([pl.lit("s__unclassified_"), pl.col("Species").str.replace_all(" ", "_")]).alias("Adjusted_Species"))
        df_failed = merged_df.filter((pl.col("Avg_Read_ANI").is_null()) | (pl.col("Avg_Read_ANI") <= 0.1)).with_columns(pl.concat_str([pl.lit("[ANI_FAILED]_"), pl.col("Species")]).alias("Adjusted_Species"))
        
        final_cols = ["Sample", "taxid", "Adjusted_Species", "Species", "Rep_Accession", "Rep_Length", "Rep_Coverage(%)", "Rep_MeanDepth", "Asm_EM_Reads", "Uniq_Reads", "Multi_Reads", "Unique(%)", "Avg_Read_ANI", "Avg_Pi", "Asm_CPM", "Asm_RPM", "Asm_FPKM", "Asm_TPM", "Asm_Rel_Abund(%)", "Segment_Accessions", "Rep_Reads"]
        existing_cols = [c for c in final_cols if c in df_confirmed.columns]
        
        if len(df_novel) > 0:
            df_novel.select(existing_cols).write_csv(str(self.output_dir / "summary" / f"all_viruses.unclassified.{ext}"), separator=sep)
        
        df_white = pl.concat([df_confirmed.select(existing_cols), df_failed.select(existing_cols)])
        if len(df_white) > 0:
            best_csv = self.output_dir / "summary" / f"all_viruses.best.summary.{ext}"
            df_white.write_csv(str(best_csv), separator=sep)
            
            plots_dir = self.output_dir / 'plots'
            with open(plots_dir / 'virus_frequency_plot.R', 'w', encoding='utf-8') as f: f.write(R_PLOT_SCRIPT)
            if shutil.which('Rscript'): subprocess.run(['Rscript', str(plots_dir / 'virus_frequency_plot.R'), '-i', str(best_csv), '-o', str(plots_dir / 'virus_analysis'), '--log10-transform'], capture_output=True)
            
        self.generate_report_txt(len(df), len(final_df), len(df_white), len(df_novel), df_white)
        
        # 早期 Reads 已经提取完毕，只需要最后提取阳性参考基因组 FASTA
        if self.args.extract and len(df_white) > 0:
            logger.info("\n📦 正在提取阳性核心参考基因组 (FASTA)...")
            target_ids = df_white['Rep_Accession'].unique().to_list()
            self.extract_fastas_with_python(target_ids, self.output_dir / 'fasta')
            logger.info("✅ 病毒序列(Reads)与参考基因组(FASTA) 全部提取完毕！")

        time.sleep(0.5)
        if not self.args.keep_bam:
            logger.info("🧹 正在清理巨型 BAM 与中间计算文件...")
            for d in ['bam']:
                d_path = self.output_dir / d
                if d_path.exists(): shutil.rmtree(d_path)

    def generate_report_txt(self, raw_c, sum_c, best_c, unc_c, best_df):
        with open(self.output_dir / 'summary' / 'analysis_report.txt', 'w', encoding='utf-8') as f:
            f.write("=" * 80 + "\n宏病毒定量全景报告 (Virus Quantification V34.3 Batch-Processing Master)\n" + "=" * 80 + "\n\n")
            f.write(f"  处理样本: {len(self.samples)} | 留存(Raw): {raw_c} | 达标(Summary): {sum_c}\n")
            f.write(f"  确诊白名单(Best): {best_c} | 降级新种(Unclassified): {unc_c}\n\n")
            if len(best_df) > 0:
                f.write("【确诊白名单概览】\n" + "-" * 60 + "\n")
                asm_counts = best_df.get_column('Adjusted_Species').value_counts()
                for row in asm_counts.iter_rows(named=True):
                    asm = row['Adjusted_Species']
                    vd = best_df.filter(pl.col('Adjusted_Species') == asm)
                    count = row.get('count', row.get('counts', len(vd)))
                    if len(vd) > 0:
                        ani_mean, pi_mean = vd['Avg_Read_ANI'].mean(), vd['Avg_Pi'].mean()
                        ani_str = f"{ani_mean:.2f}%" if ani_mean is not None else "N/A"
                        pi_str = f"{pi_mean:.5f}" if pi_mean is not None else "N/A"
                        
                        f.write(f"🎯 {asm}: 检出 {count} 例 (群体检出率 {(count/len(self.samples))*100:.1f}%)\n")
                        f.write(f"   ├─ 平均 CPM: {vd['Asm_CPM'].mean():.2f} | 平均 RPM: {vd['Asm_RPM'].mean():.2f}\n")
                        f.write(f"   ├─ 平均 FPKM: {vd['Asm_FPKM'].mean():.2f} | 平均 TPM: {vd['Asm_TPM'].mean():.2f}\n")
                        f.write(f"   ├─ 测定 ANI: {ani_str} | 多态性 Pi: {pi_str}\n")
                        f.write(f"   ├─ 代表序列长度: {vd['Rep_Length'].first():.0f} bp\n")
                        f.write(f"   ├─ 平均覆盖度: {vd['Rep_Coverage(%)'].mean():.2f}% | 平均深度: {vd['Rep_MeanDepth'].mean():.2f}x\n")
                        f.write(f"   └─ 代表株(Rep): {vd['Rep_Accession'][0]}\n\n")

    def run_pipeline(self):
        logger.info("=" * 60); logger.info("💥 Virus Quantification Pipeline (V34.3 Batch-Processing Master)"); logger.info("=" * 60)
        with Timer("整体管线运行"):
            self.setup_output_directory()
            self.build_index()
            self.find_samples()
            
            batches = [self.samples[i:i + self.args.batch_size] for i in range(0, len(self.samples), self.args.batch_size)]
            
            for i, batch in enumerate(batches):
                batch_idx = i + 1
                batch_file = self.output_dir / 'batches' / f'batch_{batch_idx:04d}.parquet'

                if self.args.resume and batch_file.exists():
                    try:
                        pl.read_parquet(batch_file)
                        logger.info(f"⏭️ 批次断点触发: [批次 {batch_idx}/{len(batches)}] 已完成，整体跳过...")
                        continue
                    except Exception:
                        logger.warning(f"⚠️ [批次 {batch_idx}/{len(batches)}] 存档文件损坏，重新计算...")

                logger.info(f"📦 开始计算 [批次 {batch_idx}/{len(batches)}] (内含 {len(batch)} 个样本)...")
                res_list = []
                with ProcessPoolExecutor(max_workers=self.args.threads) as executor:
                    futures = {executor.submit(self.process_sample, s): s for s in batch}
                    with tqdm(total=len(batch), desc=f"[批次 {batch_idx}] 定量", unit="样本", bar_format='{l_bar}{bar}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}, {postfix}]') as pbar:
                        for future in as_completed(futures):
                            sample_info = futures[future]
                            pbar.set_postfix_str(f"最新完成: {sample_info['name']}")
                            try:
                                if (result := future.result()): res_list.extend(result)
                            except Exception as e:
                                logger.error(f"❌ 样本 {sample_info['name']} 抛出致命异常: {e}")
                            pbar.update(1)
                
                if res_list:
                    pl.DataFrame(res_list, schema=METRICS_SCHEMA).write_parquet(batch_file)
                else:
                    pl.DataFrame(schema=METRICS_SCHEMA).write_parquet(batch_file)

            logger.info("🧩 所有批次处理完成，正在合并所有批次数据并执行延迟计算...")
            all_batch_files = list((self.output_dir / 'batches').glob('batch_*.parquet'))
            if not all_batch_files:
                logger.warning("⚠️ 没有发现任何批次结果，流程提前退出。")
                return
            
            combined_df = pl.concat([pl.read_parquet(f) for f in all_batch_files])
            self.summarize_results_polars(combined_df)


def main():
    parser = argparse.ArgumentParser(
        description='宏病毒鉴定与精确定量管线 V34.3',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
【参数详细说明 - Detailed Parameter Guide】

【输入与运行参数】
  -i / --input_dir    输入序列文件夹 (自动识别 FASTQ/FASTA 及 gz 命名格式)
  -r / --reference    全局参考基因组 FASTA 文件
  --ref_info          包含 Accession, taxid, Species 等信息的 TSV 文件
  --single_end        强制将输入目录中的文件按单端处理
  --sample_list       直接提供样本列表(格式: Sample R1 [R2])
  --resume            开启双重断点续传（批次级 + BAM/Reads样本级跳过）

【比对与质控参数】
  --tool              序列比对工具 (默认: bowtie2)
                      可选: bwa, bowtie2, minimap2, strobealign, bwa-mem2, hisat2
  --min_aln_len       最小比对长度 (默认: 80 bp)
  --min_pid           最小序列相似度 (默认: 0.90 = 90%%)
  --min_aln_prop      最小比对比例 (默认: 0.85 = 85%%)
  ⚠️ 若检测到输入为双端序列，将自动激活 CoverM --proper-pairs-only 逻辑

【丰度与覆盖度过滤参数】
  --min_uniq_reads    最少独特比对读数 (默认: 1)
  --min_tpm           最少 TPM 值 (默认: 0.0)
  --coverage          代表序列最小覆盖度百分比 (默认: 90.0%%)
  --meandepth         代表序列最小平均深度 (默认: 10.0x)

【进化学指标参数】
  --sp_thresh         物种 ANI 识别阈值 (默认: 95.0%%)
                      >= 95%% : 确诊种
                      > 0.1%% and < 95%% : 新种/未分类
                      <= 0.1%% : 检测失败

【性能优化与输出控制参数】
  -t / --threads      总计同时处理的样本并发数 (默认: 8)
  --align_threads     每个样本内部分配的比对线程数 (默认: 4)
  --batch_size        批次处理大小，防爆内存 (默认: 20个样本一批)
  --extract           【早期截获】基于未清洗的原始 BAM，提取并无损落地 Fastq/Fasta 和 参考序列
  --keep_bam          保留所有中间产物与 BAM 文件 (默认: False，执行清理以节约磁盘)
  --verbose           关闭静默模式，将每个样本的【比对/清洗/计算】步骤细节输出至终端
        """
    )
    
    parser.add_argument('-i', '--input_dir', required=False, help='输入 FASTQ/FASTA 文件夹')
    parser.add_argument('-r', '--reference', required=True, help='全局参考基因组 FASTA')
    parser.add_argument('-o', '--output_dir', default='./virus_out', help='输出目录 (默认: ./virus_out)')
    parser.add_argument('--ref_info', type=str, required=True, help='本地参考信息 TSV 文件')
    parser.add_argument('--taxid_clusters', type=str, help='同义 TaxID 映射文件')
    
    parser.add_argument('--single_end', action='store_true', help='强制单端模式')
    parser.add_argument('--sample_list', type=str, help='指定样本列表文件 (替代 -i)')
    parser.add_argument('--resume', action='store_true', help='断点续传/跳过已有BAM')
    parser.add_argument('--extract', action='store_true', help='提取比对到病毒的原始序列文件和参考序列')
    
    parser.add_argument('--tool', choices=['bwa', 'bowtie2', 'strobealign', 'minimap2', 'bwa-mem2', 'hisat2'], 
                        default='bowtie2', help='序列比对工具 (默认: bowtie2)')
    parser.add_argument('-t', '--threads', type=int, default=8, help='样本并行处理数 (默认: 8)')
    parser.add_argument('--align_threads', type=int, default=4, help='单样本比对线程数 (默认: 4)')
    parser.add_argument('--batch_size', type=int, default=20, help='防爆内存的批次大小 (默认: 20)')
    
    parser.add_argument('--min_aln_len', type=int, default=80, help='最小比对长度 (默认: 80 bp)')
    parser.add_argument('--min_aln_prop', type=float, default=0.85, help='最小比对比例 (默认: 0.85)')
    parser.add_argument('--min_pid', type=float, default=0.90, help='最小序列相似度 (默认: 0.90)')
    parser.add_argument('--min_uniq_reads', type=int, default=1, help='最少独特比对读数 (默认: 1)')
    parser.add_argument('--sp_thresh', type=float, default=95.0, help='物种ANI识别阈值%% (默认: 95.0)')
    
    parser.add_argument('--coverage', type=float, default=90.0, help='代表序列最小覆盖度%% (默认: 90.0)')
    parser.add_argument('--meandepth', type=float, default=10.0, help='代表序列最小平均深度x (默认: 10.0)')
    parser.add_argument('--min_tpm', type=float, default=0.0, help='最少 TPM 值 (默认: 0.0)')
    
    parser.add_argument('--keep_bam', action='store_true', help='保留 BAM 文件')
    parser.add_argument('--verbose', action='store_true', help='强制将所有子进程详情打印到屏幕')
    
    args = parser.parse_args()
    if not args.input_dir and not args.sample_list:
        parser.error("必须提供 -i (输入目录) 或 --sample_list (样本列表)")
        
    VirusQuantificationPipeline(args).run_pipeline()

if __name__ == '__main__':
    main()
