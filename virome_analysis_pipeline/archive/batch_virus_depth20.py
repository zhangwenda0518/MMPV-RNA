#!/usr/bin/env python3
"""
工业级宏病毒全栈定性定量与组装解剖管线 (Zero-Compromise Ultimate Edition)
================================================================================
【V20 底层黑科技引擎 (能力提升)】
1. [CoverM 洗流] 双端配对感知(-f 2)及三大硬核阈值过滤。
2. [图论聚类] 迭代式 Reads 共享网络聚类，彻底消灭同源交叉多重比对难题。
3. [多进程安全] 采用 JSON 缓存 Micro-metrics 完美支持 Resume。
4. [动态分类学] 多节段聚合与动态分类学降级 (ANI < 90% 自动降级为 unclassified)。
5. [Pysam 防崩] 实时擦除巨型 BAM Header，彻底治愈 viral_consensus/变异检测 崩溃。

【V13 顶级工程体系 (全量保留)】
6. [三级结果输出] all_viruses.raw / summary / best.summary 三级分流与详尽全景报告。
7. [全自动R绘图] 内置 R_PLOT_SCRIPT，生成多维度 (TPM/FPKM/Depth) 病毒丰度箱线/散点图。
8. [多进程架构] Reads提取、突变检测、共识构建全部配备独立 tqdm 进度条与并发池。
9. [优美目录系统] 基于 Taxonomy_VirusName 自动生成嵌套归档目录。
================================================================================
"""

import argparse
import os
import sys
import subprocess
import glob
import gzip
import shutil
import time
import math
import resource
import platform
import json
from datetime import datetime
from pathlib import Path
from collections import defaultdict, Counter
from concurrent.futures import ProcessPoolExecutor, as_completed
import logging
import colorlog
from tqdm import tqdm
import polars as pl
import pandas as pd
import pysam
import numpy as np
import re

# ==================== 内置 R 绘图脚本 (完整保留 v13) ====================
R_PLOT_SCRIPT = r"""#!/usr/bin/env Rscript
suppressWarnings({
  suppressPackageStartupMessages({
    library(ggplot2)
    library(dplyr)
    library(tidyr)
    library(optparse)
    if (!require("viridis", quietly = TRUE)) install.packages("viridis", dependencies = TRUE, repos="http://cran.rstudio.com/")
    library(viridis)
  })
})

option_list <- list(
  make_option(c("-i", "--input"), type = "character", default = "all_viruses.best.summary.tsv", help = "输入文件"),
  make_option(c("-o", "--output"), type = "character", default = "virus_analysis_plots", help = "输出目录或文件前缀"),
  make_option(c("-w", "--width"), type = "numeric", default = 10, help = "图片宽度"),
  make_option(c("-e", "--height"), type = "numeric", default = 8, help = "图片高度"),
  make_option(c("-m", "--modes"), type = "character", default = "all", help = "绘图模式"),
  make_option(c("-p", "--point-size"), type = "numeric", default = 3, help = "散点大小"),
  make_option(c("-d", "--dpi"), type = "numeric", default = 300, help = "图片分辨率"),
  make_option(c("-t", "--theme"), type = "character", default = "classic", help = "ggplot主题"),
  make_option(c("--log10-transform"), type = "logical", default = FALSE, action = "store_true"),
  make_option(c("--multi-plot"), type = "logical", default = FALSE, action = "store_true"),
  make_option(c("--format"), type = "character", default = "pdf", help = "输出图片格式")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (!file.exists(opt$input)) { stop(sprintf("错误: 输入文件不存在: %s", opt$input)) }

data <- if (grepl("\\.csv$", opt$input)) read.csv(opt$input, check.names = FALSE) else read.delim(opt$input, check.names = FALSE)

# 适配多节段聚合后的列名
tax_col <- if("Adjusted_Species" %in% colnames(data)) "Adjusted_Species" else if("taxonomy" %in% colnames(data)) "taxonomy" else "Virus"
acc_col <- if("Assembly" %in% colnames(data)) "Assembly" else "Virus"

data$Display_Name <- ifelse(data[[tax_col]] != "Unannotated" & data[[tax_col]] != "-", paste0(data[[tax_col]], "\n(", data[[acc_col]], ")"), data[[acc_col]])
data$Display_Name <- sapply(data$Display_Name, function(x) paste(strwrap(x, width = 40), collapse = "\n"))

# 适配度量列名
depth_col <- if("MeanDepth" %in% colnames(data)) "MeanDepth" else "MeanDepth"
fpkm_col <- if("Asm_FPKM" %in% colnames(data)) "Asm_FPKM" else "FPKM"
rpm_col <- if("Asm_RPM" %in% colnames(data)) "Asm_RPM" else "RPM"
tpm_col <- if("Asm_TPM" %in% colnames(data)) "Asm_TPM" else "TPM"

required_columns <- list("MeanDepth" = c("Display_Name", depth_col), "FPKM" = c("Display_Name", fpkm_col), "RPM" = c("Display_Name", rpm_col), "TPM" = c("Display_Name", tpm_col))
modes <- if(opt$modes == "all") names(required_columns) else strsplit(trimws(opt$modes), ",")[[1]]
available_modes <- modes[sapply(modes, function(m) all(required_columns[[m]] %in% colnames(data)))]
if (length(available_modes) == 0) stop("错误: 没有可用的分析模式。")

output_dir <- dirname(opt$output)
if (output_dir != "." && !dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

unique_viruses <- length(unique(data$Display_Name))
theme_func <- switch(opt$theme, "classic"=theme_classic, "minimal"=theme_minimal, "bw"=theme_bw, theme_bw)

if (opt$`multi-plot` && length(available_modes) > 1) {
  # 映射真实列名用于数据透视
  actual_cols <- sapply(available_modes, function(m) required_columns[[m]][2])
  plot_data_long <- data %>% select(Display_Name, Sample, all_of(actual_cols)) %>%
    pivot_longer(cols = all_of(actual_cols), names_to = "Metric", values_to = "Value") %>% filter(!is.na(Value))
  
  if (opt$`log10-transform`) {
    for (metric in unique(plot_data_long$Metric)) {
      md <- plot_data_long[plot_data_long$Metric == metric, ]
      min_p <- min(md$Value[md$Value > 0], na.rm=TRUE)
      plot_data_long$Value[plot_data_long$Metric == metric & plot_data_long$Value <= 0] <- min_p / 10
    }
  }
  
  base_metric <- actual_cols[1]
  medians <- aggregate(Value ~ Display_Name, data=plot_data_long[plot_data_long$Metric==base_metric,], median)
  plot_data_long$Display_Name <- factor(plot_data_long$Display_Name, levels=medians[order(medians$Value), "Display_Name"])
  
  p_facet <- ggplot(plot_data_long, aes(x=Display_Name, y=Value)) +
    geom_boxplot(aes(fill=Display_Name), alpha=0.6, outlier.shape=NA) +
    geom_point(aes(color=Display_Name), position=position_jitter(width=0.2, height=0), alpha=0.6) +
    facet_wrap(~ Metric, scales="free_x", ncol=length(available_modes)) +
    scale_fill_viridis_d(option="turbo") + scale_color_viridis_d(option="turbo") +
    labs(title="Comprehensive Multi-metric Virus Abundance", x="Taxonomy (Accession)", y="Value") +
    theme_func(base_size=13) + theme(plot.title=element_text(hjust=0.5, face="bold"), axis.text.y=element_text(size=10, face="italic"), legend.position="none") + coord_flip()
  
  if (opt$`log10-transform`) p_facet <- p_facet + scale_y_log10()
  
  fn_f <- sprintf("%s_multi_metrics%s.%s", opt$output, ifelse(opt$`log10-transform`, "_log10", ""), opt$format)
  dynamic_height_facet <- max(opt$height, unique_viruses * 0.8)
  ggsave(fn_f, plot=p_facet, width=opt$width * 1.5, height=dynamic_height_facet, dpi=opt$dpi, bg="white")
}
cat("所有图表生成完成!\n")
"""

def setup_logging(verbose=False):
    handler = colorlog.StreamHandler()
    handler.setFormatter(colorlog.ColoredFormatter(
        '%(log_color)s[%(asctime)s] %(levelname)s - %(message)s',
        datefmt='%H:%M:%S',
        log_colors={'DEBUG': 'cyan', 'INFO': 'green', 'WARNING': 'yellow', 'ERROR': 'red', 'CRITICAL': 'bold_red'}
    ))
    logger = colorlog.getLogger("VirusPipe")
    logger.addHandler(handler)
    logger.setLevel(logging.DEBUG if verbose else logging.INFO)
    return logger

logger = setup_logging()

class Timer:
    def __init__(self, name=""): self.name = name
    def __enter__(self):
        self.start_time = time.time()
        logger.info(f"⏱️  开始: {self.name}")
        return self
    def __exit__(self, exc_type, exc_val, exc_tb):
        duration = time.time() - self.start_time
        if exc_type is None: logger.info(f"✅ 完成: {self.name} [耗时: {self.format_duration(duration)}]")
        else: logger.error(f"❌ 失败: {self.name} [耗时: {self.format_duration(duration)}]")
    @staticmethod
    def format_duration(seconds):
        if seconds < 60: return f"{seconds:.1f}秒"
        elif seconds < 3600: return f"{seconds/60:.1f}分钟"
        else: return f"{seconds/3600:.1f}小时"

class VirusMetagenomicsPipeline:
    def __init__(self, args):
        self.args = args
        self.samples = []
        self.output_dir = Path(args.output_dir)
        self.index_path = None
        self.ref_length_dict = {} 
        
        self.df_best = None
        self.total_raw = 0
        self.total_summary = 0

        self.setup_output_directory()
        self.check_env()
        self._load_reference_lengths()

    def check_env(self):
        required_tools = {'samtools': 'BAM处理', 'coverm': '序列清洗过滤', 'pandepth': '深度计算', self.args.tool: '序列比对'}
        if self.args.tool == 'bowtie2': required_tools['bowtie2-build'] = 'Bowtie2建库'
        elif self.args.tool == 'hisat2': required_tools['hisat2-build'] = 'HISAT2建库'
        if self.args.consensus: required_tools['viral_consensus'] = '共识构建'
        if getattr(self.args, 'call_variants', False): required_tools[self.args.variant_caller] = '变异检测'
        if getattr(self.args, 'extract_reads', False): required_tools['pigz'] = '多线程压缩'
            
        missing = [f"{t} ({d})" for t, d in required_tools.items() if shutil.which(t) is None]
        if missing:
            logger.error("❌ 缺少必要工具:")
            for m in missing: logger.error(f"  - {m}")
            sys.exit(1)
        logger.info("✅ 所有必要工具可用")

    def _load_reference_lengths(self):
        with open(self.args.reference, 'r') as f:
            curr_id = None
            curr_len = 0
            for line in f:
                if line.startswith('>'):
                    if curr_id: self.ref_length_dict[curr_id] = curr_len
                    curr_id = line.strip().split()[0][1:]
                    curr_len = 0
                else:
                    curr_len += len(line.strip())
            if curr_id: self.ref_length_dict[curr_id] = curr_len

    def setup_output_directory(self):
        subdirs = ['bam', 'stat', 'summary', 'logs', 'index', 'plots']
        if getattr(self.args, 'extract_reads', False): subdirs.append('reads')
        if self.args.consensus: subdirs.append('consensus')
        if getattr(self.args, 'call_variants', False): subdirs.append('variants')
            
        with Timer("创建输出目录结构"):
            for subdir in subdirs: 
                (self.output_dir / subdir).mkdir(parents=True, exist_ok=True)
            log_file = self.output_dir / 'logs' / f'pipeline_fusion_{datetime.now().strftime("%Y%m%d_%H%M%S")}.log'
            fh = logging.FileHandler(log_file, encoding='utf-8')
            fh.setFormatter(logging.Formatter('[%(asctime)s] %(levelname)s - %(message)s'))
            logger.addHandler(fh)

    def write_sample_log(self, sample_name, content, mode='a'):
        log_file = self.output_dir / 'logs' / f"{sample_name}.log"
        with open(log_file, mode, encoding='utf-8') as f: f.write(content + "\n")

    def build_index(self):
        tool = self.args.tool
        ref_path = Path(self.args.reference).resolve()
        if tool in ['strobealign', 'minimap2']:
            self.index_path = str(ref_path)
            return
            
        prefix = self.output_dir / 'index' / f"{ref_path.stem}_{tool}"
        self.index_path = str(prefix)
        
        cmd = []
        if tool == 'bwa': cmd = ['bwa', 'index', '-p', str(prefix), str(ref_path)]
        elif tool == 'bwa-mem2': cmd = ['bwa-mem2', 'index', '-p', str(prefix), str(ref_path)]
        elif tool == 'bowtie2' and not (Path(f"{prefix}.1.bt2").exists() or Path(f"{prefix}.1.bt2l").exists()):
            cmd = ['bowtie2-build', '--threads', str(self.args.threads), str(ref_path), str(prefix)]
        elif tool == 'hisat2': cmd = ['hisat2-build', '-p', str(self.args.threads), str(ref_path), str(prefix)]
            
        if cmd:
            logger.info(f"🏗️ 构建 {tool} 索引...")
            subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    def find_samples(self):
        input_dir = Path(self.args.input_dir) if self.args.input_dir else None
        fq_exts = ['.fq', '.fastq', '.fq.gz', '.fastq.gz']
        samples = []
        if input_dir and self.args.single_end:
            for ext in fq_exts:
                for fastq in input_dir.glob(f'*{ext}'):
                    sname = fastq.name.replace(ext, '').replace('.gz', '')
                    for suf in ['_1', '_R1', '.R1', '.1', '_unmapped', '.unmapped', '_trimmed']:
                        if suf in sname: sname = sname.split(suf)[0]
                    samples.append({'name': sname, 'r1': str(fastq), 'r2': None})
        elif input_dir:
            patterns = [('*_1.f*q*', '*_2.f*q*'), ('*_R1*.f*q*', '*_R2*.f*q*'), ('*.R1.*', '*.R2.*'), 
                        ('*.1.f*q*', '*.2.f*q*'), ('*_unmapped.R1.fq.gz', '*_unmapped.R2.fq.gz')]
            found_files = set()
            for p1, p2 in patterns:
                for r1_file in input_dir.glob(p1):
                    if r1_file in found_files: continue
                    r1_name = r1_file.name
                    r2_name = r1_name.replace('_1.', '_2.').replace('_R1', '_R2').replace('.R1.', '.R2.').replace('.1.', '.2.')
                    if (input_dir / r2_name).exists():
                        sname = r1_name.split('.')[0]
                        for suf in ['_1', '_R1', '_unmapped']: sname = sname.replace(suf, '')
                        samples.append({'name': sname, 'r1': str(r1_file), 'r2': str(input_dir / r2_name)})
                        found_files.update([r1_file, input_dir / r2_name])
                        
        if self.args.sample_list:
            with open(self.args.sample_list, 'r') as f:
                for line in f:
                    p = line.strip().split()
                    if len(p) >= 2: samples.append({'name': p[0], 'r1': p[1], 'r2': p[2] if len(p)>2 else None})
                    
        self.samples = list({v['name']: v for v in samples}.values())
        if not self.samples: sys.exit(1)
        logger.info(f"✅ 找到 {len(self.samples)} 个测序样本")

    # ==================== V20 核心洗流引擎 ====================
    def align_and_coverm_filter(self, sample):
        sname = sample['name']
        bam_dir = self.output_dir / 'bam'
        raw_bam = bam_dir / f"{sname}.raw.bam"
        filt_sec_bam = bam_dir / f"{sname}.filtered_with_sec.bam"
        metrics_cache = self.output_dir / 'stat' / f"{sname}_total_mapped.txt"
        
        if self.args.resume and filt_sec_bam.exists() and metrics_cache.exists():
            with open(metrics_cache, 'r') as f: total_mapped = int(f.read().strip())
            return filt_sec_bam, total_mapped

        self.write_sample_log(sname, f"\n--- 1. 比对与 CoverM 严格洗流 ({self.args.tool}) ---")
        
        align_cmd = []
        if self.args.tool == 'bowtie2':
            align_cmd = ['bowtie2', '-p', str(self.args.align_threads), '-x', self.index_path]
            align_cmd.extend(['-1', sample['r1'], '-2', sample['r2']] if sample['r2'] else ['-U', sample['r1']])
        elif self.args.tool == 'bwa':
            align_cmd = ['bwa', 'mem', '-v', '1', '-t', str(self.args.align_threads), self.index_path, sample['r1']]
            if sample['r2']: align_cmd.append(sample['r2'])

        view_cmd = ['samtools', 'view', '-b', '-o', str(raw_bam)]
        ap = subprocess.Popen(align_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        vp = subprocess.Popen(view_cmd, stdin=ap.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        ap.stdout.close()
        ap.communicate()
        vp.communicate()

        count_res = subprocess.run(['samtools', 'view', '-c', '-F', '4', str(raw_bam)], capture_output=True, text=True)
        total_mapped = int(count_res.stdout.strip() or 0)
        with open(metrics_cache, 'w') as f: f.write(str(total_mapped))

        pid_val = self.args.min_pid * 100 if self.args.min_pid <= 1.0 else self.args.min_pid
        prop_val = self.args.min_aln_prop * 100 if self.args.min_aln_prop <= 1.0 else self.args.min_aln_prop
        filter_cmd = [
            'coverm', 'filter', '-b', str(raw_bam), '-o', str(filt_sec_bam) + '.unsorted',
            '--min-read-aligned-length', str(self.args.min_aln_len),
            '--min-read-percent-identity', str(pid_val),
            '--min-read-aligned-percent', str(prop_val),
            '--include-secondary', '-t', str(self.args.align_threads)
        ]
        if self.args.strict_paired and sample['r2']: filter_cmd.append('--proper-pairs-only')
        
        subprocess.run(filter_cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        subprocess.run(['samtools', 'sort', '-@', str(min(4, self.args.threads)), '-o', str(filt_sec_bam), str(filt_sec_bam) + '.unsorted'], check=True)
        pysam.index(str(filt_sec_bam))
        
        raw_bam.unlink()
        Path(str(filt_sec_bam) + '.unsorted').unlink()
        return filt_sec_bam, total_mapped

    def build_read_sharing_network(self, sname, bam_path):
        self.write_sample_log(sname, f"--- 2. 图论网络聚类 (消除多重比对冗余) ---")
        read_to_contigs = defaultdict(set)
        contig_primary_reads = defaultdict(int)
        
        with pysam.AlignmentFile(bam_path, "rb") as bam:
            for read in bam.fetch(until_eof=True):
                if read.is_unmapped: continue
                contig = bam.get_reference_name(read.reference_id)
                read_to_contigs[read.query_name].add(contig)
                if not read.is_secondary and not read.is_supplementary:
                    contig_primary_reads[contig] += 1
                    
        shared_counts = defaultdict(int)
        for qname, contigs in read_to_contigs.items():
            clist = list(contigs)
            if len(clist) > 1:
                for i in range(len(clist)):
                    for j in range(i + 1, len(clist)):
                        pair = tuple(sorted([clist[i], clist[j]]))
                        shared_counts[pair] += 1
                        
        adjacency = defaultdict(set)
        for (c1, c2), shared_n in shared_counts.items():
            min_total = min(contig_primary_reads.get(c1, 0), contig_primary_reads.get(c2, 0))
            if min_total > 0 and (shared_n / min_total) >= self.args.cluster_thresh:
                adjacency[c1].add(c2)
                adjacency[c2].add(c1)
                
        visited = set()
        clusters = []
        for contig in contig_primary_reads.keys():
            if contig not in visited:
                cluster = set()
                queue = [contig]
                while queue:
                    curr = queue.pop(0)
                    if curr not in visited:
                        visited.add(curr)
                        cluster.add(curr)
                        queue.extend(list(adjacency[curr] - visited))
                clusters.append(cluster)
                
        exemplars = set()
        for cl in clusters:
            sorted_cl = sorted(list(cl), key=lambda x: -contig_primary_reads.get(x, 0))
            exemplars.add(sorted_cl[0])
            
        return exemplars

    def generate_primary_bam_and_metrics(self, sname, filt_sec_bam, exemplars):
        final_bam_path = self.output_dir / 'bam' / f"{sname}.final.bam"
        json_cache = self.output_dir / 'stat' / f"{sname}.micro_metrics.json"
        
        if self.args.resume and final_bam_path.exists() and json_cache.exists():
            with open(json_cache, 'r') as f: return final_bam_path, json.load(f)

        metrics = []
        with pysam.AlignmentFile(filt_sec_bam, "rb") as in_bam:
            out_bam = pysam.AlignmentFile(final_bam_path, "wb", header=in_bam.header)
            for contig in exemplars:
                reads = list(in_bam.fetch(contig))
                primary_reads = [r for r in reads if not r.is_secondary and not r.is_supplementary]
                if not primary_reads: continue
                
                unique_n, multi_n, ani_list = 0, 0, []
                contig_len = in_bam.get_reference_length(contig)
                base_counts = [Counter() for _ in range(contig_len)]
                
                for r in primary_reads:
                    out_bam.write(r)
                    if r.mapping_quality >= 10: unique_n += 1
                    else: multi_n += 1
                    
                    aln_len = r.query_alignment_length
                    if aln_len > 0:
                        try: nm = r.get_tag('NM')
                        except KeyError: nm = 0 
                        ani_list.append((aln_len - nm) / aln_len)
                    
                    for qpos, ref_pos, ref_base in r.get_aligned_pairs(matches_only=True, with_seq=True):
                        if ref_pos is not None and ref_pos < contig_len and qpos is not None:
                            base_counts[ref_pos][r.query_sequence[qpos].upper()] += 1

                avg_read_ani = float(np.mean(ani_list)) if ani_list else 0.0
                pi_list = []
                for c in base_counts:
                    total = sum(c.values())
                    if total > 0:
                        max_count = max(c.values())
                        pi_list.append((total - max_count) / total)
                avg_pi = float(np.mean(pi_list)) if pi_list else 0.0

                metrics.append({
                    'Sample': sname, 'Accession': contig, 'MappedReads': len(primary_reads),
                    'Unique_Reads': unique_n, 'Multi_Reads': multi_n,
                    'Read_ANI': round(avg_read_ani, 4), 'Diversity_Pi': round(avg_pi, 4)
                })
            out_bam.close()
        
        pysam.index(str(final_bam_path))
        with open(json_cache, 'w') as f: json.dump(metrics, f)
        if filt_sec_bam.exists(): filt_sec_bam.unlink()
        return final_bam_path, metrics

    def run_pandepth(self, final_bam, sname, metrics, total_mapped):
        stat_prefix = self.output_dir / 'stat' / sname
        if not (self.args.resume and Path(str(stat_prefix) + '.chr.stat.gz').exists()):
            subprocess.run(['pandepth', '-a', '-i', str(final_bam), '-o', str(stat_prefix), '-t', str(min(10, self.args.align_threads))], capture_output=True)
        
        stat_file = str(stat_prefix) + '.chr.stat.gz'
        depth_dict = {}
        with gzip.open(stat_file, 'rt') as f:
            for line in f:
                if line.startswith('#'): continue
                p = line.strip().split()
                if len(p) >= 6: depth_dict[p[0]] = {'CoveredBases': int(p[2]), 'TotalDepth': int(p[3]), 'Coverage(%)': float(p[4]), 'MeanDepth': float(p[5])}
        
        for m in metrics:
            d = depth_dict.get(m['Accession'], {})
            m['CoveredBases'] = d.get('CoveredBases', 0)
            m['TotalDepth'] = d.get('TotalDepth', 0)
            m['Coverage(%)'] = d.get('Coverage(%)', 0.0)
            m['MeanDepth'] = d.get('MeanDepth', 0.0)
            
            length = self.ref_length_dict.get(m['Accession'], 1)
            m['Length'] = length
            
            m['RPM'] = (m['MappedReads'] * 1e6) / total_mapped if total_mapped else 0
            m['Relative_Abundance(%)'] = round((m['MappedReads'] / total_mapped * 100), 4) if total_mapped else 0
            m['FPKM'] = (m['MappedReads'] * 1e9) / (total_mapped * length) if total_mapped and length else 0
            m['_RPK'] = (m['MappedReads'] * 1000) / length if length else 0

        sum_rpk = sum(m['_RPK'] for m in metrics)
        for m in metrics:
            m['TPM'] = (m['_RPK'] * 1e6) / sum_rpk if sum_rpk else 0
            del m['_RPK']
            
        return metrics

    def process_single_sample(self, sample):
        sname = sample['name']
        log_mode = 'a' if self.args.resume else 'w'
        self.write_sample_log(sname, f"\n{'='*50}\n[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] 开始处理: {sname}", mode=log_mode)
        try:
            filt_sec_bam, total_mapped = self.align_and_coverm_filter(sample)
            exemplars = self.build_read_sharing_network(sname, str(filt_sec_bam))
            if not exemplars: return []
            
            final_bam, micro_metrics = self.generate_primary_bam_and_metrics(sname, filt_sec_bam, exemplars)
            full_metrics = self.run_pandepth(final_bam, sname, micro_metrics, total_mapped)
            return full_metrics
        except Exception as e:
            logger.error(f"❌ 样本 {sname} 致命错误: {e}")
            return []

    # ==================== V13/V20 分类学与三级输出体系 ====================
    def fast_load_rvdb_taxonomy(self, found_viruses):
        if not getattr(self.args, 'rvdb_taxon', None) or not found_viruses: return {}
        target_exact = set(found_viruses)
        target_base = {v.split('.')[0] for v in found_viruses}
        taxon_dict = {}
        try:
            with open(self.args.rvdb_taxon, 'r', encoding='utf-8') as f:
                header = None
                for line in f:
                    if line.startswith('####'): continue
                    parts = line.rstrip('\n').split('\t')
                    if header is None:
                        header = parts
                        idx_acc, idx_desc = header.index('accession'), header.index('description')
                        idx_tax = header.index('taxonomy')
                        continue
                    if len(parts) <= max(idx_acc, idx_desc, idx_tax): continue
                    acc = parts[idx_acc].strip()
                    acc_base = acc.split('.')[0]
                    if acc in target_exact or acc_base in target_base:
                        info = {'description': parts[idx_desc], 'taxonomy': parts[idx_tax]}
                        taxon_dict[acc] = info
                        taxon_dict[acc_base] = info
            return taxon_dict
        except Exception: return {}

    def fetch_ncbi_taxonomy(self, accessions):
        try: from Bio import Entrez, SeqIO
        except ImportError: return {}
        ncbi_dict = {}
        cache_file = self.output_dir / 'summary' / 'ncbi_taxonomy_cache.tsv'
        
        if cache_file.exists() and cache_file.stat().st_size > 0:
            try:
                df_cache = pd.read_csv(cache_file, sep='\t', dtype=str)
                for _, row in df_cache.iterrows(): ncbi_dict[row['accession']] = {'description': str(row['description']), 'taxonomy': str(row['taxonomy'])}
            except Exception: cache_file.unlink()
                
        to_fetch = [acc for acc in accessions if acc not in ncbi_dict]
        if not to_fetch: return ncbi_dict

        Entrez.email = getattr(self.args, 'email', None)
        Entrez.api_key = getattr(self.args, 'api_key', None)
        
        write_header = not cache_file.exists()
        with open(cache_file, 'a', encoding='utf-8') as cf:
            if write_header: cf.write("accession\tdescription\ttaxonomy\n")
            for acc in to_fetch:
                try:
                    handle = Entrez.efetch(db="nuccore", id=acc, rettype="gb", retmode="text")
                    record = SeqIO.read(handle, "genbank")
                    handle.close()
                    desc, tax = record.description, record.annotations.get('organism', 'Unannotated')
                    ncbi_dict[acc] = {'description': desc, 'taxonomy': tax}
                    cf.write(f"{acc}\t{desc}\t{tax}\n")
                    cf.flush()
                except: pass
                time.sleep(0.3)
        return ncbi_dict

    def summarize_results(self, all_metrics):
        if not all_metrics: return None
        df = pl.DataFrame(all_metrics)
        
        found_viruses = df['Accession'].unique().to_list()
        tax_map = {}
        if getattr(self.args, 'rvdb_taxon', None): tax_map.update(self.fast_load_rvdb_taxonomy(found_viruses))
        unannotated = [v for v in found_viruses if v not in tax_map]
        if unannotated and getattr(self.args, 'email', None): tax_map.update(self.fetch_ncbi_taxonomy(unannotated))
        
        def get_tax(acc, key): return tax_map.get(acc, {}).get(key, "Unannotated")
        df = df.with_columns([
            pl.col("Accession").map_elements(lambda x: get_tax(x, 'taxonomy'), return_dtype=pl.Utf8).alias("Assembly"), 
            pl.col("Accession").map_elements(lambda x: get_tax(x, 'description'), return_dtype=pl.Utf8).alias("description")
        ])
        
        ext = self.args.format
        sep = "," if ext == "csv" else "\t"
        
        # 1. 原始输出 (Raw)
        df.write_csv(str(self.output_dir / "summary" / f"all_viruses.raw.{ext}"), separator=sep)
        self.total_raw = len(df)
        
        # 多节段物理聚合
        asm_df = df.group_by(["Sample", "Assembly"]).agg([
            pl.col("MappedReads").sum(), pl.col("Unique_Reads").sum(), pl.col("Multi_Reads").sum(),
            pl.col("CoveredBases").sum(), pl.col("Length").sum().alias("Asm_Length"),
            pl.col("Coverage(%)").mean().alias("Coverage(%)"), pl.col("MeanDepth").mean().alias("MeanDepth"),
            pl.col("Read_ANI").mean().alias("Avg_Read_ANI"), pl.col("Diversity_Pi").mean().alias("Avg_Pi"),
            pl.col("TPM").sum().alias("Asm_TPM"), pl.col("RPM").sum().alias("Asm_RPM"),
            pl.col("FPKM").sum().alias("Asm_FPKM"), pl.col("Relative_Abundance(%)").sum().alias("Asm_Rel_Abund(%)"),
            pl.col("description").first(), pl.col("Accession").alias("Segment_Accessions")
        ])
        asm_df = asm_df.with_columns((pl.col("Unique_Reads") / pl.col("MappedReads") * 100).round(2).alias("Unique(%)"))
        
        # 动态分类学降级
        sp_thresh = self.args.sp_thresh
        asm_df = asm_df.with_columns(
            pl.when(pl.col("Avg_Read_ANI") < sp_thresh)
            .then(pl.concat_str([pl.lit("s__unclassified_"), pl.col("Assembly").str.replace_all(" ", "_")]))
            .otherwise(pl.col("Assembly")).alias("Adjusted_Species")
        ).with_columns(pl.col("Segment_Accessions").list.join(","))
        
        # 2. 严格阈值过滤 (Summary)
        final_df = asm_df.filter(
            (pl.col("Coverage(%)") >= self.args.coverage) &
            (pl.col("MeanDepth") >= self.args.meandepth) &
            (pl.col("Asm_FPKM") >= self.args.min_fpkm) &
            (pl.col("Asm_TPM") >= self.args.min_tpm)
        )
        final_df.write_csv(str(self.output_dir / "summary" / f"all_viruses.summary.{ext}"), separator=sep)
        self.total_summary = len(final_df)
        
        # 3. 最优代表株筛选 (Best: 针对同类 Taxonomy，优先完整度，其次深度)
        best_df = final_df.sort(
            ["Sample", "Adjusted_Species", "Coverage(%)", "MeanDepth"], descending=[False, False, True, True]
        ).unique(subset=["Sample", "Adjusted_Species"], keep="first", maintain_order=True)
        
        best_csv = self.output_dir / "summary" / f"all_viruses.best.summary.{ext}"
        best_df.write_csv(str(best_csv), separator=sep)
        self.df_best = best_df
        
        logger.info(f"✅ 三级结果保存完成! (Raw: {self.total_raw} -> Summary: {self.total_summary} -> Best: {len(best_df)})")
        
        self.generate_report_txt()
        self.generate_plots(best_csv)
        
        # 返回 Best 结果字典，供下游流程无缝对接
        return best_df.to_dicts()

    # ==================== V13 下游工程体系 (全量恢复) ====================
    def generate_plots(self, best_summary_file):
        logger.info("\n📊 开始调用 R 脚本生成全自动可视化图表...")
        plots_dir = self.output_dir / 'plots'
        r_script_path = plots_dir / 'virus_frequency_plot.R'
        with open(r_script_path, 'w', encoding='utf-8') as f: f.write(R_PLOT_SCRIPT)
            
        if shutil.which('Rscript') is None:
            logger.warning("⚠️ 未检测到 Rscript 环境，跳过绘图步骤。")
            return
            
        cmd = ['Rscript', str(r_script_path), '-i', str(best_summary_file), '-o', str(plots_dir / 'virus_analysis'), '--multi-plot', '--log10-transform']
        try:
            subprocess.run(cmd, check=True, capture_output=True, text=True)
            logger.info(f"✅ 分析图表生成成功，已保存至: {plots_dir}/")
        except subprocess.CalledProcessError as e: logger.error(f"❌ R 脚本图表生成失败:\n{e.stderr}")

    def generate_report_txt(self):
        report_file = self.output_dir / 'summary' / 'analysis_report.txt'
        df_best = self.df_best
        total_samples = len(self.samples)
        
        with open(report_file, 'w') as f:
            f.write("=" * 80 + "\n")
            f.write(f"工业级宏病毒丰度与多节段全景分析报告 (Fusion Edition)\n")
            f.write("=" * 80 + "\n\n")
            
            f.write(f"【管线配置与阈值】\n")
            f.write(f"  比对工具: {self.args.tool}\n")
            f.write(f"  序列清洗 (CoverM): 长度>{self.args.min_aln_len}, 匹配度>{self.args.min_pid}, 覆盖比例>{self.args.min_aln_prop}\n")
            f.write(f"  图论网络聚类阈值: {self.args.cluster_thresh} (用于根除多重交叉比对)\n")
            f.write(f"  分类学降级 ANI 阈值: {self.args.sp_thresh} (低于此值将降级为 s__unclassified)\n")
            f.write(f"  宏观丰度硬核阈值: Coverage >= {self.args.coverage}%, MeanDepth >= {self.args.meandepth}, FPKM >= {self.args.min_fpkm}\n")
            f.write(f"  变异/共识生成: {'开启' if self.args.call_variants or self.args.consensus else '未开启'}\n\n")

            f.write(f"【核心数据流漏斗统计】\n")
            f.write(f"  1. 处理样本总数: {total_samples}\n")
            f.write(f"  2. 全量检测留存 (Raw): {self.total_raw} 条 (已通过图论洗流)\n")
            f.write(f"  3. 过滤达标的病毒 (Summary): {self.total_summary} 条\n")
            f.write(f"  4. 同分类下最优代表基因组 (Best): {len(df_best)} 条\n\n")

            if df_best is not None and len(df_best) > 0:
                f.write("【最终确诊核心群落全景 (基于 Best Representative)】\n")
                f.write("-" * 60 + "\n")
                asm_counts = df_best['Adjusted_Species'].value_counts()
                count_col_name = 'count' if 'count' in asm_counts.columns else 'counts'
                for row in asm_counts.iter_rows(named=True):
                    asm, count = row['Adjusted_Species'], row[count_col_name]
                    vd = df_best.filter(pl.col('Adjusted_Species') == asm)
                    mean_tpm, mean_fpkm, mean_rpm = vd['Asm_TPM'].mean(), vd['Asm_FPKM'].mean(), vd['Asm_RPM'].mean()
                    f.write(f"🎯 {asm}:\n")
                    f.write(f"   ├─ 样本检出数: {count}/{total_samples}\n")
                    f.write(f"   ├─ 核心平均 TPM:  {mean_tpm:.2f} (该毒株在整个群落中的相对占比)\n")
                    f.write(f"   ├─ 核心平均 FPKM: {mean_fpkm:.2f} (长度归一化后的绝对载量)\n")
                    f.write(f"   └─ 核心平均 RPM:  {mean_rpm:.2f} (无偏差的绝对测序文库载量)\n\n")

    def extract_fastas_with_python(self, target_ids, virus_to_tax, out_dir):
        """恢复 v13 优美的 Taxonomy 嵌套目录提取法"""
        target_set = set(target_ids)
        current_id, current_seq = None, []
        
        def save_record(vid, seq_list):
            if vid in target_set:
                vname = re.sub(r'[\\/*?:"<>| ]', "_", vid)
                tax = virus_to_tax.get(vid, 'Unannotated')
                safe_tax = re.sub(r'[\\/*?:"<>| ]', "_", tax)
                if len(safe_tax) > 50: safe_tax = safe_tax[:50]
                folder_name = f"{safe_tax}_{vname}"
                v_dir = out_dir / folder_name
                v_dir.mkdir(parents=True, exist_ok=True)
                with open(v_dir / f"{folder_name}.ref.fasta", 'w') as outf:
                    outf.write(f">{vid}\n" + "".join(seq_list) + "\n")
                    
        with open(self.args.reference, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    if current_id: save_record(current_id, current_seq)
                    current_id = line[1:].split()[0]
                    current_seq = []
                else: current_seq.append(line)
            if current_id: save_record(current_id, current_seq)

    def _extract_pure_bam_pysam(self, global_bam_path, acc, out_bam_path):
        """恢复 v20 的无敌 Pysam 纯化术，斩断 awk 的表头错误"""
        try:
            with pysam.AlignmentFile(global_bam_path, "rb") as in_bam:
                target_tid = in_bam.get_tid(acc)
                if target_tid < 0: return False
                new_header = {'HD': {'VN': '1.0', 'SO': 'coordinate'}, 'SQ': [{'LN': in_bam.lengths[target_tid], 'SN': acc}]}
                with pysam.AlignmentFile(out_bam_path, "wb", header=new_header) as out_bam:
                    for read in in_bam.fetch(acc):
                        read.reference_id = 0
                        if read.next_reference_id == target_tid: read.next_reference_id = 0
                        elif read.next_reference_id >= 0:
                            read.next_reference_id = -1
                            read.mate_is_unmapped = True
                        out_bam.write(read)
            return True
        except Exception: return False

    def _run_single_read_extraction(self, sample):
        bam_file = self.output_dir / 'bam' / f"{sample}.final.bam"
        reads_dir = self.output_dir / 'reads'
        r1 = reads_dir / f"{sample}_virus_1.fastq.gz"
        r2 = reads_dir / f"{sample}_virus_2.fastq.gz"
        rs = reads_dir / f"{sample}_virus_single.fastq.gz"
        
        if self.args.resume and (r1.exists() or r2.exists() or rs.exists()): return True
        if not bam_file.exists(): return False
        try:
            cmd = f"samtools fastq --threads 4 -F 4 -1 >(pigz -c > '{r1}') -2 >(pigz -c > '{r2}') '{bam_file}' | pigz -c > '{rs}'"
            subprocess.run(cmd, shell=True, check=True, executable='/bin/bash', stderr=subprocess.DEVNULL)
            for f in [r1, r2, rs]:
                if f.exists() and f.stat().st_size <= 50: f.unlink()
            return True
        except: return False

    def extract_mapped_reads(self, best_results):
        logger.info("\n📦 开始并行提取比对上病毒的序列 Reads (仅限阳性样本)...")
        reads_dir = self.output_dir / 'reads'
        reads_dir.mkdir(parents=True, exist_ok=True)
        unique_samples = list(set(r['Sample'] for r in best_results))
        success_count = 0
        with Timer("阳性病毒 Reads 提取"):
            with ProcessPoolExecutor(max_workers=min(len(unique_samples), self.args.threads)) as executor:
                futures = {executor.submit(self._run_single_read_extraction, s): s for s in unique_samples}
                with tqdm(total=len(unique_samples), desc="Reads提取", bar_format='{l_bar}{bar}| {n_fmt}/{total_fmt}[{elapsed}<{remaining}]') as pbar:
                    for future in as_completed(futures):
                        if future.result(): success_count += 1
                        pbar.update(1)
            logger.info(f"✅ 成功提取了 {success_count} 个样本的 Reads")

    def _run_single_variant_calling(self, task):
        sample, virus, taxonomy = task
        vname = re.sub(r'[\\/*?:"<>| ]', "_", virus)
        safe_tax = re.sub(r'[\\/*?:"<>| ]', "_", taxonomy)[:50]
        folder_name = f"{safe_tax}_{vname}"
        virus_dir = self.output_dir / 'variants' / folder_name
        
        bam_file = self.output_dir / 'bam' / f"{sample}.final.bam"
        ref_fasta = virus_dir / f"{folder_name}.ref.fasta"
        fixed_bam = virus_dir / f"{sample}.{folder_name}.fixed.bam"
        out_file = virus_dir / f"{sample}.{folder_name}.variants.vcf"
        
        if self.args.resume and out_file.exists(): return True
        if not ref_fasta.exists() or not bam_file.exists(): return False
        
        try:
            if not Path(str(ref_fasta) + '.fai').exists(): subprocess.run(f"samtools faidx '{ref_fasta}'", shell=True, check=True)
            if not self._extract_pure_bam_pysam(str(bam_file), virus, str(fixed_bam)): return False
            subprocess.run(f"samtools index '{fixed_bam}'", shell=True, check=True)
            
            if self.args.variant_caller == 'freebayes':
                subprocess.run(f"freebayes -p 1 -f '{ref_fasta}' '{fixed_bam}' > '{out_file}'", shell=True, stderr=subprocess.DEVNULL)
            elif self.args.variant_caller == 'ivar':
                prefix = str(out_file).replace('.vcf', '')
                subprocess.run(f"samtools mpileup -A -aa -d 0 -Q 0 --reference '{ref_fasta}' '{fixed_bam}' | ivar variants -r '{ref_fasta}' -p '{prefix}'", shell=True, stderr=subprocess.DEVNULL)
                
            fixed_bam.unlink(); Path(str(fixed_bam)+'.bai').unlink()
            return True
        except: return False

    def run_variants_calling(self, best_results):
        logger.info(f"\n🔬 开始并行执行变异检测 ({self.args.variant_caller}, 按 Taxonomy 归档)...")
        virus_to_tax = {}
        tasks = []
        for r in best_results:
            accs = r['Segment_Accessions'].split(',')
            for acc in accs:
                virus_to_tax[acc] = r.get('Adjusted_Species', 'Unannotated')
                tasks.append((r['Sample'], acc, r.get('Adjusted_Species', 'Unannotated')))
                
        with Timer("核心变异检测"):
            self.extract_fastas_with_python(list(virus_to_tax.keys()), virus_to_tax, self.output_dir / 'variants')
            success_count = 0
            with ProcessPoolExecutor(max_workers=self.args.threads) as executor:
                futures = {executor.submit(self._run_single_variant_calling, t): t for t in tasks}
                with tqdm(total=len(tasks), desc="变异检测", bar_format='{l_bar}{bar}| {n_fmt}/{total_fmt}[{elapsed}<{remaining}]') as pbar:
                    for future in as_completed(futures):
                        if future.result(): success_count += 1
                        pbar.update(1)
            logger.info(f"✅ 成功提取 {success_count} 个核心变异检测结果")

    def _run_single_consensus(self, task):
        sample, virus, mean_depth, taxonomy = task
        vname = re.sub(r'[\\/*?:"<>| ]', "_", virus)
        safe_tax = re.sub(r'[\\/*?:"<>| ]', "_", taxonomy)[:50]
        folder_name = f"{safe_tax}_{vname}"
        virus_dir = self.output_dir / 'consensus' / folder_name
        
        bam_file = self.output_dir / 'bam' / f"{sample}.final.bam"
        ref_fasta = virus_dir / f"{folder_name}.ref.fasta"
        fixed_bam = virus_dir / f"{sample}.{folder_name}.fixed.bam"
        out_fasta = virus_dir / f"{sample}.{folder_name}.consensus.fasta"
        
        if self.args.resume and out_fasta.exists(): return True
        if not ref_fasta.exists() or not bam_file.exists(): return False
        
        min_depth = int(math.floor(mean_depth))
        if min_depth > 10: min_depth = 10
        elif min_depth < 1: min_depth = 1
        
        try:
            if not self._extract_pure_bam_pysam(str(bam_file), virus, str(fixed_bam)): return False
            subprocess.run(f"samtools index '{fixed_bam}'", shell=True, check=True)
            subprocess.run(f"viral_consensus -i '{fixed_bam}' -r '{ref_fasta}' -o '{out_fasta}' --min_depth {min_depth}", shell=True, check=True, stderr=subprocess.DEVNULL)
            fixed_bam.unlink(); Path(str(fixed_bam)+'.bai').unlink()
            return True
        except: return False

    def build_consensus_sequences(self, best_results):
        logger.info("\n🧬 开始并行生成核心共识序列 (动态阈值 + Taxonomy 归档)...")
        virus_to_tax = {}
        tasks = []
        for r in best_results:
            accs = r['Segment_Accessions'].split(',')
            for acc in accs:
                virus_to_tax[acc] = r.get('Adjusted_Species', 'Unannotated')
                tasks.append((r['Sample'], acc, r['MeanDepth'], r.get('Adjusted_Species', 'Unannotated')))
                
        with Timer("核心共识构建"):
            self.extract_fastas_with_python(list(virus_to_tax.keys()), virus_to_tax, self.output_dir / 'consensus')
            success_count = 0
            with ProcessPoolExecutor(max_workers=self.args.threads) as executor:
                futures = {executor.submit(self._run_single_consensus, t): t for t in tasks}
                with tqdm(total=len(tasks), desc="共识生成", bar_format='{l_bar}{bar}| {n_fmt}/{total_fmt}[{elapsed}<{remaining}]') as pbar:
                    for future in as_completed(futures):
                        if future.result(): success_count += 1
                        pbar.update(1)
            logger.info(f"✅ 成功归档 {success_count} 个核心共识序列文件")

    # ==================== 管线执行主轴 ====================
    def run(self):
        logger.info("=" * 60)
        logger.info("💥 Ultimate Metagenomics Virus Pipeline (Zero-Compromise)")
        logger.info("=" * 60)
        
        with Timer("整体管线运行"):
            self.build_index()
            self.find_samples()
            
            max_parallel_jobs = self.args.threads
            if getattr(self.args, 'parallel', False):
                max_parallel_jobs = self.args.parallel_jobs if self.args.parallel_jobs else max(1, self.args.threads // self.args.align_threads)

            all_metrics = []
            with ProcessPoolExecutor(max_workers=max_parallel_jobs) as executor:
                futures = {executor.submit(self.process_single_sample, s): s for s in self.samples}
                with tqdm(total=len(self.samples), desc="主比对与微观提取") as pbar:
                    for future in as_completed(futures):
                        res = future.result()
                        if res: all_metrics.extend(res)
                        pbar.update(1)
                        
            best_results = self.summarize_results(all_metrics)
            
            if best_results:
                if getattr(self.args, 'extract_reads', False): self.extract_mapped_reads(best_results)
                if getattr(self.args, 'call_variants', False): self.run_variants_calling(best_results)
                if getattr(self.args, 'consensus', False): self.build_consensus_sequences(best_results)
            else:
                logger.warning("⚠️ 未找到任何满足阈值的病毒记录")

        logger.info("=" * 60)
        logger.info("✨ 分析管线全流程运行结束 (所有特性已满载)")

def main():
    parser = argparse.ArgumentParser(description='工业级宏病毒全栈管线 (Zero-Compromise Edition)', formatter_class=argparse.RawDescriptionHelpFormatter)
    
    io_g = parser.add_argument_group('输入输出参数')
    io_g.add_argument('-i', '--input_dir', help='输入文件夹路径')
    io_g.add_argument('-r', '--reference', required=True, help='参考基因组fasta文件')
    io_g.add_argument('-o', '--output_dir', default='./virus_out', help='输出目录')
    io_g.add_argument('-l', '--sample_list', help='样本列表文件')
    io_g.add_argument('--single_end', action='store_true', help='使用单端数据')

    align_g = parser.add_argument_group('高级流程与变异控制参数')
    align_g.add_argument('--tool', choices=['strobealign', 'minimap2', 'bwa', 'bwa-mem2', 'bowtie2', 'hisat2'], default='bowtie2', help='比对工具')
    align_g.add_argument('--align_threads', type=int, default=8, help='每样本比对线程')
    align_g.add_argument('--extract_reads', action='store_true', help='提取阳性reads并输出为fastq.gz')
    align_g.add_argument('--consensus', action='store_true', help='生成共识序列')
    align_g.add_argument('--call_variants', action='store_true', help='开启变异检测')
    align_g.add_argument('--variant_caller', choices=['freebayes', 'ivar', 'lofreq'], default='freebayes', help='变异检测工具')
    align_g.add_argument('--resume', action='store_true', help='断点续跑')

    flt_g = parser.add_argument_group('高级过滤与图论聚类 (V20 特性)')
    flt_g.add_argument('--strict_paired', action='store_true', help='强制要求 Proper Pair')
    flt_g.add_argument('--min_aln_len', type=int, default=80, help='CoverM: 最小比词长度')
    flt_g.add_argument('--min_aln_prop', type=float, default=0.85, help='CoverM: Read被覆盖比例')
    flt_g.add_argument('--min_pid', type=float, default=0.90, help='CoverM: 最小序列相似度')
    flt_g.add_argument('--cluster_thresh', type=float, default=0.30, help='图论: Read-Sharing 聚类阈值')
    flt_g.add_argument('--sp_thresh', type=float, default=0.90, help='分类学: 降级为 unclassified 的 ANI 阈值')
    
    thresh_g = parser.add_argument_group('丰度过滤参数 (V13 特性)')
    thresh_g.add_argument('--coverage', type=float, default=90.0, help='Coverage(%%)阈值')
    thresh_g.add_argument('--meandepth', type=float, default=10.0, help='MeanDepth阈值')
    thresh_g.add_argument('--min_fpkm', type=float, default=0.0, help='FPKM最小值过滤')
    thresh_g.add_argument('--min_tpm', type=float, default=0.0, help='TPM最小值过滤')

    para_g = parser.add_argument_group('并行处理参数')
    para_g.add_argument('--parallel', action='store_true', help='并行处理样本')
    para_g.add_argument('--parallel_jobs', type=int, help='并行任务数')
    para_g.add_argument('-t', '--threads', type=int, default=4, help='总线程数')

    ctrl_g = parser.add_argument_group('兼容与数据库参数')
    ctrl_g.add_argument('--format', choices=['csv', 'tsv'], default='tsv', help='输出表格格式')
    ctrl_g.add_argument('--rvdb_taxon', type=str, help='RVDB 分类注释文件')
    ctrl_g.add_argument('--email', type=str, help='NCBI 联网邮箱')
    ctrl_g.add_argument('--api_key', type=str, help='NCBI API Key')
    ctrl_g.add_argument('--verbose', action='store_true', help='输出 DEBUG 日志')
    
    args = parser.parse_args()
    if not args.input_dir and not args.sample_list: parser.error("必须提供 -i 或 -l")
    if args.verbose: logger.setLevel(logging.DEBUG)
    
    VirusMetagenomicsPipeline(args).run()

if __name__ == '__main__':
    main()
