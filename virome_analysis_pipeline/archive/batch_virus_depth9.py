#!/usr/bin/env python3
"""
病毒序列比对、深度统计、多维度定量分析、共识序列构建与可视化完整流程
1. 包含 FPKM, RPM, TPM 三大定量体系
2. 包含 RVDB 极速流式解析 + NCBI Entrez 智能联网补全 (针对 Unannotated)
3. 包含最优同源代表株筛选机制 (Best Representative by Taxonomy)
4. 包含智能深度截断与按子目录整理的共识序列构建
5. 包含断点续跑(--resume)、单样本全量日志追踪及全自动 R 绘图出图
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
from pathlib import Path
from collections import defaultdict
from datetime import datetime
from concurrent.futures import ProcessPoolExecutor, as_completed
import logging
import colorlog
from tqdm import tqdm
import pandas as pd
import re
import io

# ==================== 内置 R 绘图脚本 ====================
R_PLOT_SCRIPT = r"""#!/usr/bin/env Rscript
# 加载必要的包
suppressWarnings({
  suppressPackageStartupMessages({
    library(ggplot2)
    library(dplyr)
    library(tidyr)
    library(optparse)
    if (!require("viridis", quietly = TRUE)) install.packages("viridis", dependencies = TRUE, repos="http://cran.rstudio.com/")
    library(viridis)
    library(reshape2)
  })
})

option_list <- list(
  make_option(c("-i", "--input"), type = "character", default = "all_viruses.summary.tsv", help = "输入文件"),
  make_option(c("-o", "--output"), type = "character", default = "virus_analysis_plots", help = "输出目录或文件前缀"),
  make_option(c("-w", "--width"), type = "numeric", default = 10, help = "图片宽度"),
  make_option(c("-e", "--height"), type = "numeric", default = 8, help = "图片高度"),
  make_option(c("-m", "--modes"), type = "character", default = "all", help = "绘图模式"),
  make_option(c("-p", "--point-size"), type = "numeric", default = 3, help = "散点大小"),
  make_option(c("-d", "--dpi"), type = "numeric", default = 300, help = "图片分辨率"),
  make_option(c("-t", "--theme"), type = "character", default = "classic", help = "ggplot主题"),
  make_option(c("-c", "--color-palette"), type = "character", default = "Set2", help = "颜色调色板"),
  make_option(c("--title"), type = "character", default = "Virus Analysis", help = "图表标题前缀"),
  make_option(c("-s", "--show-stats"), type = "logical", default = FALSE, action = "store_true"),
  make_option(c("--show-names"), type = "logical", default = FALSE, action = "store_true"),
  make_option(c("--order-by-mean"), type = "logical", default = FALSE, action = "store_true"),
  make_option(c("--log10-transform"), type = "logical", default = FALSE, action = "store_true"),
  make_option(c("--multi-plot"), type = "logical", default = FALSE, action = "store_true"),
  make_option(c("--format"), type = "character", default = "png", help = "输出图片格式")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (!file.exists(opt$input)) { stop(sprintf("错误: 输入文件不存在: %s", opt$input)) }
cat(sprintf("正在读取数据文件: %s\n", opt$input))

# 智能兼容 CSV 和 TSV
if (grepl("\\.csv$", opt$input)) {
  data <- read.csv(opt$input, stringsAsFactors = FALSE, check.names = FALSE)
} else {
  data <- read.delim(opt$input, stringsAsFactors = FALSE, check.names = FALSE)
}

required_columns <- list("MeanDepth" = c("Virus", "MeanDepth"), "FPKM" = c("Virus", "FPKM"), "RPM" = c("Virus", "RPM"), "TPM" = c("Virus", "TPM"))
modes <- if(opt$modes == "all") names(required_columns) else strsplit(trimws(opt$modes), ",")[[1]]

available_modes <- c()
for (mode in modes) {
  if (all(required_columns[[mode]] %in% colnames(data))) { available_modes <- c(available_modes, mode) }
}
if (length(available_modes) == 0) stop("错误: 没有可用的分析模式。")

output_dir <- dirname(opt$output)
if (output_dir != "." && !dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

prepare_data <- function(data, value_column) {
  if (!"Sample" %in% colnames(data)) data$Sample <- paste0("Sample", 1:nrow(data))
  data <- data[!is.na(data[[value_column]]), ]
  if (opt$`order-by-mean`) {
    value_means <- aggregate(data[[value_column]] ~ Virus, data, mean)
    data$Virus <- factor(data$Virus, levels = value_means[order(-value_means[[value_column]]), "Virus"])
  } else {
    data$Virus <- factor(data$Virus, levels = unique(data$Virus))
  }
  data <- data %>% group_by(Virus) %>% mutate(
    Q1 = quantile(.data[[value_column]], 0.25, na.rm = TRUE),
    Q3 = quantile(.data[[value_column]], 0.75, na.rm = TRUE),
    IQR = Q3 - Q1,
    lower_bound = Q1 - 1.5 * IQR, upper_bound = Q3 + 1.5 * IQR,
    is_outlier = .data[[value_column]] < lower_bound | .data[[value_column]] > upper_bound
  ) %>% ungroup()
  return(data)
}

theme_func <- switch(opt$theme, "classic"=theme_classic, "minimal"=theme_minimal, "bw"=theme_bw, "light"=theme_light, theme_classic)
color_scale <- if(opt$`color-palette` == "viridis") scale_color_viridis_d() else scale_color_brewer(palette = opt$`color-palette`)
fill_scale <- if(opt$`color-palette` == "viridis") scale_fill_viridis_d() else scale_fill_brewer(palette = opt$`color-palette`)

create_boxplot <- function(data, value_column, title_suffix = "", y_label = value_column) {
  plot_data <- prepare_data(data, value_column)
  y_trans <- NULL
  if (opt$`log10-transform`) {
    if (any(plot_data[[value_column]] <= 0, na.rm = TRUE)) {
      min_pos <- min(plot_data[[value_column]][plot_data[[value_column]] > 0], na.rm = TRUE)
      plot_data[[value_column]][plot_data[[value_column]] <= 0] <- min_pos / 10
    }
    y_trans <- scale_y_log10()
    y_label <- paste0(y_label, " (log10)")
  }
  
  p <- ggplot(plot_data, aes(x = Virus, y = .data[[value_column]])) +
    geom_boxplot(aes(fill = Virus), alpha = 0.7, outlier.shape = NA, width = 0.6) +
    geom_point(aes(color = Virus), position = position_jitter(width = 0.2, height = 0), size = opt$`point-size`, alpha = 0.7) +
    fill_scale + color_scale +
    labs(title = paste(opt$title, title_suffix), x = "Virus", y = y_label) +
    theme_func(base_size = 14) +
    theme(plot.title = element_text(hjust=0.5, face="bold"), axis.text.x = element_text(angle=45, hjust=1), legend.position="none")
  
  if (!is.null(y_trans)) p <- p + y_trans + annotation_logticks(sides = "l")
  return(p)
}

for (mode in available_modes) {
  p <- create_boxplot(data, mode, paste("-", mode, "Distribution"))
  fn <- sprintf("%s_%s_boxplot%s.%s", opt$output, tolower(mode), ifelse(opt$`log10-transform`, "_log10", ""), opt$format)
  ggsave(fn, plot=p, width=opt$width, height=opt$height, dpi=opt$dpi, bg="white")
}

if (opt$`multi-plot` && length(available_modes) > 1) {
  plot_data_long <- data %>% select(Virus, Sample, all_of(available_modes)) %>%
    pivot_longer(cols = all_of(available_modes), names_to = "Metric", values_to = "Value") %>% filter(!is.na(Value))
  
  plot_data_long$Metric <- factor(plot_data_long$Metric, levels = available_modes)
  
  if (opt$`log10-transform`) {
    for (metric in unique(plot_data_long$Metric)) {
      md <- plot_data_long[plot_data_long$Metric == metric, ]
      if (any(md$Value <= 0, na.rm=TRUE)) {
        min_p <- min(md$Value[md$Value > 0], na.rm=TRUE)
        plot_data_long$Value[plot_data_long$Metric == metric & plot_data_long$Value <= 0] <- min_p / 10
      }
    }
  }
  
  p_facet <- ggplot(plot_data_long, aes(x=Virus, y=Value)) +
    geom_boxplot(aes(fill=Virus), alpha=0.7, outlier.shape=NA) +
    geom_point(aes(color=Virus), position=position_jitter(width=0.2, height=0), alpha=0.5) +
    facet_wrap(~ Metric, scales="free_y", ncol=2) +
    fill_scale + color_scale +
    labs(title=paste(opt$title, "- Multiple Metrics"), x="Virus", y="Value") +
    theme_func(base_size=12) +
    theme(axis.text.x=element_text(angle=45, hjust=1), legend.position="none")
  
  if (opt$`log10-transform`) p_facet <- p_facet + scale_y_log10()
  
  fn_f <- sprintf("%s_multi_metrics%s.%s", opt$output, ifelse(opt$`log10-transform`, "_log10", ""), opt$format)
  ggsave(fn_f, plot=p_facet, width=opt$width*1.5, height=opt$height*1.2, dpi=opt$dpi, bg="white")
}
cat("所有图表生成完成!\n")
"""

# 设置彩色全局日志
def setup_logging(verbose=False):
    handler = colorlog.StreamHandler()
    handler.setFormatter(colorlog.ColoredFormatter(
        '%(log_color)s%(asctime)s - %(levelname)s - %(message)s',
        datefmt='%H:%M:%S',
        log_colors={'DEBUG': 'cyan', 'INFO': 'green', 'WARNING': 'yellow', 'ERROR': 'red', 'CRITICAL': 'bold_red'}
    ))
    logger = colorlog.getLogger(__name__)
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
        else: logger.error(f"❌ 失败: {self.name}[耗时: {self.format_duration(duration)}]")
    @staticmethod
    def format_duration(seconds):
        if seconds < 60: return f"{seconds:.1f}秒"
        elif seconds < 3600: return f"{seconds/60:.1f}分钟"
        else: return f"{seconds/3600:.1f}小时"

class VirusAnalysisPipeline:
    def __init__(self, args):
        self.args = args
        self.samples =[]
        self.output_dir = Path(args.output_dir)
        self.index_path = None
        
        # 针对 --resume 的参数隐式设置
        if self.args.resume:
            self.args.skip_alignment = True
            self.args.skip_depth = True

        self.check_tools()
        
    def check_tools(self):
        required_tools = {'samtools': '测序数据处理', 'pandepth': '深度计算', self.args.tool: '序列比对'}
        if self.args.tool == 'bowtie2': required_tools['bowtie2-build'] = 'Bowtie2索引构建'
        elif self.args.tool == 'hisat2': required_tools['hisat2-build'] = 'HISAT2索引构建'
        if self.args.consensus:
            required_tools['viral_consensus'] = '共识序列构建'
            required_tools['awk'] = 'BAM头处理工具'
            
        missing_tools =[f"{t} ({d})" for t, d in required_tools.items() if shutil.which(t) is None]
        if missing_tools:
            logger.error("❌ 缺少必要的工具:")
            for m in missing_tools: logger.error(f"  - {m}")
            sys.exit(1)
        logger.info("✅ 所有必要工具可用")

    def write_sample_log(self, sample_name, content, mode='a'):
        log_file = self.output_dir / 'logs' / f"{sample_name}.log"
        with open(log_file, mode, encoding='utf-8') as f:
            f.write(content + "\n")

    def fast_load_rvdb_taxonomy(self, found_viruses):
        if not self.args.rvdb_taxon or not found_viruses: return None
        logger.info(f"📂 正在快速流式扫描本地 RVDB 注释文件...")
        target_exact = set(found_viruses)
        target_base = {v.split('.')[0] for v in found_viruses}
        taxon_dict = {}
        found_count = 0
        try:
            with open(self.args.rvdb_taxon, 'r', encoding='utf-8') as f:
                header = None
                idx_acc, idx_desc, idx_tax, idx_taxid = -1, -1, -1, -1
                for line in f:
                    if line.startswith('####'): continue
                    parts = line.rstrip('\n').split('\t')
                    if header is None:
                        header = parts
                        try:
                            idx_acc, idx_desc = header.index('accession'), header.index('description')
                            idx_tax, idx_taxid = header.index('taxonomy'), header.index('taxonomy_id')
                        except ValueError: return None
                        continue
                    if len(parts) <= max(idx_acc, idx_desc, idx_tax, idx_taxid): continue
                    acc = parts[idx_acc].strip()
                    acc_base = acc.split('.')[0]
                    if acc in target_exact or acc_base in target_base or acc in target_base:
                        info = {'description': parts[idx_desc] or '', 'taxonomy': parts[idx_tax] or '', 'taxonomy_id': parts[idx_taxid] or ''}
                        taxon_dict[acc], taxon_dict[acc_base] = info, info
                        found_count += 1
            logger.info(f"✅ 本地扫描完成，成功提取到 {found_count} 条分类信息")
            return taxon_dict
        except Exception as e:
            logger.error(f"❌ 读取 RVDB 失败: {e}")
            return None

    def fetch_ncbi_taxonomy(self, accessions):
        try:
            from Bio import Entrez, SeqIO
        except ImportError:
            logger.warning("⚠️ 找不到 biopython 模块，跳过 NCBI 在线解析。")
            return {}

        Entrez.email = self.args.email
        if self.args.api_key:
            Entrez.api_key = self.args.api_key
            sleep_time = 0.11
        else:
            sleep_time = 0.35

        ncbi_dict = {}
        logger.info(f"🌐 正在通过 NCBI 联网获取 {len(accessions)} 个本地未注释病毒的分类信息...")

        for acc in tqdm(accessions, desc="NCBI 解析进度", unit="seq"):
            max_retries = 3
            success = False
            for attempt in range(max_retries):
                try:
                    handle = Entrez.efetch(db="nuccore", id=acc, rettype="gb", retmode="text")
                    record = SeqIO.read(handle, "genbank")
                    handle.close()
                    tax_id = "-"
                    for feature in record.features:
                        if feature.type == "source":
                            for db_xref in feature.qualifiers.get('db_xref',[]):
                                if db_xref.startswith('taxon:'): tax_id = db_xref.split(':')[1]; break
                            break
                    ncbi_dict[acc] = {'description': record.description, 'taxonomy': record.annotations.get('organism', '-'), 'taxonomy_id': tax_id}
                    success = True
                    break
                except Exception as e:
                    time.sleep(1)
            if not success:
                ncbi_dict[acc] = {'description': "Unannotated (NCBI Error)", 'taxonomy': "Unannotated", 'taxonomy_id': "-"}
            time.sleep(sleep_time)

        return ncbi_dict
    
    def setup_output_directory(self):
        subdirs =['bam', 'stat', 'summary', 'logs', 'index', 'plots']
        if self.args.consensus: subdirs.append('consensus')
        with Timer("创建输出目录结构"):
            for subdir in subdirs: (self.output_dir / subdir).mkdir(parents=True, exist_ok=True)
            has_file_handler = any(isinstance(h, logging.FileHandler) for h in logger.handlers)
            if not has_file_handler:
                log_file = self.output_dir / 'logs' / f'analysis_global_{datetime.now().strftime("%Y%m%d_%H%M%S")}.log'
                file_handler = logging.FileHandler(log_file, encoding='utf-8')
                file_handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s', datefmt='%Y-%m-%d %H:%M:%S'))
                logger.addHandler(file_handler)
        return self.output_dir
        
    def build_index(self):
        tool = self.args.tool
        ref_path = Path(self.args.reference).resolve()
        if tool in ['strobealign', 'minimap2']:
            self.index_path = str(ref_path)
            return
        prefix = self.output_dir / 'index' / f"{ref_path.stem}_{tool}"
        self.index_path = str(prefix)
        
        if tool == 'bwa':
            if Path(f"{prefix}.bwt").exists(): return
            cmd = ['bwa', 'index', '-p', str(prefix), str(ref_path)]
        elif tool == 'bwa-mem2':
            if Path(f"{prefix}.bwt.2bit.64").exists() or Path(f"{prefix}.0123").exists(): return
            cmd =['bwa-mem2', 'index', '-p', str(prefix), str(ref_path)]
        elif tool == 'bowtie2':
            if Path(f"{prefix}.1.bt2").exists() or Path(f"{prefix}.1.bt2l").exists(): return
            cmd =['bowtie2-build', '--threads', str(self.args.threads), str(ref_path), str(prefix)]
        elif tool == 'hisat2':
            if Path(f"{prefix}.1.ht2").exists() or Path(f"{prefix}.1.ht2l").exists(): return
            cmd =['hisat2-build', '-p', str(self.args.threads), str(ref_path), str(prefix)]
            
        logger.info(f"🏗️ 构建 {tool} 索引...")
        with Timer(f"构建 {tool} 索引"):
            result = subprocess.run(cmd, capture_output=True, text=True)
            if result.returncode != 0:
                logger.error(f"❌ 索引构建失败: {result.stderr}")
                sys.exit(1)

    def find_samples(self):
        input_dir = Path(self.args.input_dir) if self.args.input_dir else None
        fq_exts =['.fq', '.fastq', '.fq.gz', '.fastq.gz']
        samples =[]
        if input_dir and self.args.single_end:
            for ext in fq_exts:
                for fastq in input_dir.glob(f'*{ext}'):
                    sname = fastq.name.replace(ext, '').replace('.gz', '')
                    for suf in['_1', '_R1', '.R1', "R1_10239",'.1', '_unmapped', '.unmapped', '_trimmed', '.trimmed']:
                        if suf in sname: sname = sname.split(suf)[0]
                    samples.append({'name': sname, 'r1': str(fastq), 'r2': None})
        elif input_dir:
            patterns =[('*_1.f*q*', '*_2.f*q*'), ('*_R1*.f*q*', '*_R2*.f*q*'), ('*.R1.*', '*.R2.*'), 
                        ('*.1.f*q*', '*.2.f*q*'), ('*_1_*.f*q*', '*_2_*.f*q*'), 
                        ('*_unmapped.R1.fq.gz', '*_unmapped.R2.fq.gz'), ('*unmapped.R1.fq.gz', '*unmapped.R2.fq.gz'), 
                        ('*_1_unmapped.*', '*_2_unmapped.*'), ('*.unmapped.R1_10239.*', '*.unmapped.R2_10239.*')]
            found_files = set()
            for p1, p2 in patterns:
                for r1_file in input_dir.glob(p1):
                    if r1_file in found_files: continue
                    r1_name, r2_name = r1_file.name, None
                    if '_1.' in r1_name: r2_name = r1_name.replace('_1.', '_2.')
                    elif '_R1' in r1_name: r2_name = r1_name.replace('_R1', '_R2')
                    elif '.R1.' in r1_name: r2_name = r1_name.replace('.R1.', '.R2.')
                    elif '.1.' in r1_name: r2_name = r1_name.replace('.1.', '.2.')
                    elif '_1_' in r1_name: r2_name = r1_name.replace('_1_', '_2_')
                    elif '.unmapped.R1.' in r1_name: r2_name = r1_name.replace('.unmapped.R1.', '.unmapped.R2.')
                    elif 'unmapped.R1.' in r1_name: r2_name = r1_name.replace('unmapped.R1.', 'unmapped.R2.')
                    elif '_1_unmapped.' in r1_name: r2_name = r1_name.replace('_1_unmapped.', '_2_unmapped.')
                    elif '.unmapped.R1_10239.' in r1_name: r2_name = r1_name.replace('.unmapped.R1_10239.', '.unmapped.R2_10239.')
                    
                    if r2_name and (input_dir / r2_name).exists():
                        sname = r1_name
                        for suf in['_1', '_R1', '.R1', '.1', '_unmapped', '.unmapped', '_1_unmapped']:
                            if suf in sname: sname = sname.split(suf)[0]
                        sname = sname.split('.')[0]
                        samples.append({'name': sname, 'r1': str(r1_file), 'r2': str(input_dir / r2_name)})
                        found_files.update([r1_file, input_dir / r2_name])
            if not samples:
                all_fqs =[]
                for ext in fq_exts: all_fqs.extend(input_dir.glob(f'*{ext}'))
                fg = defaultdict(list)
                for f in all_fqs:
                    prefix = re.sub(r'[._](R?[12]|unmapped\.R?[12])[._].*', '', f.name)
                    if '.' in prefix: prefix = prefix.split('.')[0]
                    fg[prefix].append(str(f))
                for prefix, files in fg.items():
                    if len(files) == 2:
                        r1, r2 = None, None
                        for f in files:
                            if any(p in f for p in['_1', '_R1', '.R1.', '.1.', 'unmapped.R1']): r1 = f
                            elif any(p in f for p in['_2', '_R2', '.R2.', '.2.', 'unmapped.R2']): r2 = f
                        if r1 and r2: samples.append({'name': prefix, 'r1': r1, 'r2': r2})
        if self.args.sample_list:
            with open(self.args.sample_list, 'r') as f:
                for line in f:
                    parts = line.strip().split()
                    if len(parts) >= 2 and not line.startswith('#') and Path(parts[1]).exists():
                        samples.append({'name': parts[0], 'r1': parts[1], 'r2': parts[2] if len(parts)>2 and Path(parts[2]).exists() else None})
        
        unique_samples, seen = [], set()
        for s in samples:
            if s['name'] not in seen: unique_samples.append(s); seen.add(s['name'])
        self.samples = unique_samples
        if not self.samples:
            logger.error("❌ 未找到任何样本文件。")
            sys.exit(1)
        logger.info(f"✅ 找到 {len(self.samples)} 个测序样本")
    
    def align_sample(self, sample):
        sample_name = sample['name']
        bam_dir = self.output_dir / 'bam'
        sorted_bam = bam_dir / f'{sample_name}.sorted.bam'
        
        if sorted_bam.exists() and self.args.skip_alignment:
            self.write_sample_log(sample_name, f"⏭️[断点续跑] 跳过比对，直接使用已有BAM: {sorted_bam}")
            return sorted_bam
            
        self.write_sample_log(sample_name, f"\n--- 1. 序列比对 (使用 {self.args.tool}) ---")
        tool = self.args.tool
        align_cmd = []
        if tool == 'strobealign':
            align_cmd =['strobealign', '-t', str(self.args.align_threads), self.index_path, sample['r1']]
            if sample['r2']: align_cmd.append(sample['r2'])
        elif tool == 'minimap2':
            preset = 'sr' if sample['r2'] else 'map-ont'
            align_cmd =['minimap2', '-ax', preset, '-t', str(self.args.align_threads), self.index_path, sample['r1']]
            if sample['r2']: align_cmd.append(sample['r2'])
        elif tool == 'bwa':
            align_cmd =['bwa', 'mem', '-v', '1', '-t', str(self.args.align_threads), self.index_path, sample['r1']]
            if sample['r2']: align_cmd.append(sample['r2'])
        elif tool == 'bwa-mem2':
            align_cmd =['bwa-mem2', 'mem', '-v', '1', '-t', str(self.args.align_threads), self.index_path, sample['r1']]
            if sample['r2']: align_cmd.append(sample['r2'])
        elif tool == 'bowtie2':
            align_cmd =['bowtie2', '-p', str(self.args.align_threads), '-x', self.index_path]
            if sample['r2']: align_cmd.extend(['-1', sample['r1'], '-2', sample['r2']])
            else: align_cmd.extend(['-U', sample['r1']])
        elif tool == 'hisat2':
            align_cmd =['hisat2', '-p', str(self.args.align_threads), '-x', self.index_path]
            if sample['r2']: align_cmd.extend(['-1', sample['r1'], '-2', sample['r2']])
            else: align_cmd.extend(['-U', sample['r1']])
        
        sort_cmd =['samtools', 'sort', '-@', str(min(4, self.args.threads)), '-o', str(sorted_bam)]
        self.write_sample_log(sample_name, f"Align CMD: {' '.join(align_cmd)}")
        self.write_sample_log(sample_name, f"Sort CMD: {' '.join(sort_cmd)}")
        
        align_proc = subprocess.Popen(align_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        sort_proc = subprocess.Popen(sort_cmd, stdin=align_proc.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        align_proc.stdout.close()
        _, align_stderr = align_proc.communicate()
        sort_stdout, sort_stderr = sort_proc.communicate()
        
        if align_proc.returncode != 0:
            self.write_sample_log(sample_name, f"❌ 比对失败:\n{align_stderr}")
            return None
        if sort_proc.returncode != 0:
            self.write_sample_log(sample_name, f"❌ 排序失败:\n{sort_stderr}")
            return None
            
        self.write_sample_log(sample_name, f"✅ 比对成功。\n[Stderr 截取输出]:\n{align_stderr[:1000]}")
        
        index_cmd =['samtools', 'index', str(sorted_bam)]
        subprocess.run(index_cmd, capture_output=True, text=True)
        return sorted_bam
            
    def get_idxstats(self, bam_file, sample_name):
        self.write_sample_log(sample_name, "\n--- 2. 提取测序绝对文库大小 (samtools idxstats) ---")
        cmd =['samtools', 'idxstats', str(bam_file)]
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            self.write_sample_log(sample_name, f"❌ idxstats 失败:\n{result.stderr}")
            return None
            
        stats, total_reads, total_rpk = {}, 0, 0.0
        for line in result.stdout.strip().split('\n'):
            parts = line.split('\t')
            if len(parts) >= 4:
                ref, length, mapped, unmapped = parts[0], int(parts[1]), int(parts[2]), int(parts[3])
                stats[ref] = mapped
                total_reads += (mapped + unmapped)
                if length > 0 and mapped > 0: total_rpk += (mapped * 1000.0) / length
                
        self.write_sample_log(sample_name, f"Total Library Reads: {total_reads} | Total RPK: {total_rpk:.2f}")
        return {'stats': stats, 'global_total_reads': total_reads, 'total_rpk': total_rpk}
    
    def run_pandepth(self, bam_file, sample_name):
        output_prefix = self.output_dir / 'stat' / sample_name
        stat_file = output_prefix.with_suffix('.chr.stat.gz')
        if stat_file.exists() and self.args.skip_depth:
            self.write_sample_log(sample_name, f"⏭️  [断点续跑] 跳过深度计算，使用已有文件: {stat_file}")
            return stat_file
            
        self.write_sample_log(sample_name, "\n--- 3. 深度计算 (pandepth) ---")
        pandepth_cmd =['pandepth', '-a', '-i', str(bam_file), '-o', str(output_prefix), '-t', str(min(10, self.args.threads))]
        result = subprocess.run(pandepth_cmd, capture_output=True, text=True)
        if result.returncode != 0:
            self.write_sample_log(sample_name, f"❌ pandepth 失败:\n{result.stderr}")
            return None
        self.write_sample_log(sample_name, "✅ pandepth 运行成功。")
        return stat_file
    
    def parse_stat_file(self, stat_file, sample_name, idx_data):
        if not idx_data: return[]
        global_total_reads = idx_data['global_total_reads']
        total_rpk = idx_data['total_rpk']
        idx_stats = idx_data['stats']
        results =[]
        try:
            with gzip.open(stat_file, 'rt') as f:
                for line in f:
                    if not line or line.startswith('#'): continue
                    parts = line.split()
                    if len(parts) < 6: continue
                    virus_name, length, covered_site, total_depth, coverage, mean_depth = parts
                    
                    length_val, covered_site_val, total_depth_val = int(length), int(covered_site), int(total_depth)
                    coverage_val, mean_depth_val = float(coverage), float(mean_depth)
                    mapped_reads = idx_stats.get(virus_name, 0)
                    
                    rpm = (mapped_reads * 1e6) / global_total_reads if global_total_reads > 0 else 0.0
                    fpkm = (mapped_reads * 1e9) / (global_total_reads * length_val) if global_total_reads > 0 and length_val > 0 else 0.0
                    tpm = ((mapped_reads * 1000.0) / length_val * 1e6) / total_rpk if total_rpk > 0 and length_val > 0 else 0.0
                        
                    if coverage_val > self.args.coverage and mean_depth_val > self.args.meandepth:
                        if fpkm >= self.args.min_fpkm and tpm >= self.args.min_tpm:
                            results.append({
                                'Sample': sample_name, 'Virus': virus_name, 'Length': length_val,
                                'CoveredSite': covered_site_val, 'TotalDepth': total_depth_val,
                                'Coverage(%)': coverage_val, 'MeanDepth': mean_depth_val,
                                'MappedReads': mapped_reads, 'RPM': round(rpm, 2),
                                'FPKM': round(fpkm, 2), 'TPM': round(tpm, 2)
                            })
            self.write_sample_log(sample_name, f"✅ 数据解析完毕, 筛选出 {len(results)} 条符合阈值的记录。")
            return results
        except Exception as e:
            self.write_sample_log(sample_name, f"❌ 统计文件解析异常: {str(e)}")
            return[]
    
    def process_sample(self, sample):
        sample_name = sample['name']
        start_time = datetime.now()
        
        log_mode = 'a' if self.args.resume else 'w'
        status_str = "[断点续跑]" if self.args.resume else "[全新运行]"
        self.write_sample_log(sample_name, f"\n{'='*50}\n[{start_time.strftime('%Y-%m-%d %H:%M:%S')}] {status_str} 开始处理测序样本: {sample_name}", mode=log_mode)
        
        try:
            bam_file = self.align_sample(sample)
            if not bam_file: return None
            
            idx_data = self.get_idxstats(bam_file, sample_name)
            if not idx_data: return None
            
            stat_file = self.run_pandepth(bam_file, sample_name)
            if not stat_file: return None
            
            res = self.parse_stat_file(stat_file, sample_name, idx_data)
            
            end_time = datetime.now()
            duration = end_time - start_time
            self.write_sample_log(sample_name, f"[{end_time.strftime('%Y-%m-%d %H:%M:%S')}] 🎉 样本整体处理流程结束 [耗时: {Timer.format_duration(duration.total_seconds())}]")
            return res
        except Exception as e:
            self.write_sample_log(sample_name, f"❌ 样本执行期间发生致命错误: {str(e)}")
            return None
    
    def run_pipeline(self):
        start_time = datetime.now()
        logger.info("=" * 60)
        logger.info(f"🚀 开始病毒分析管线 (工具: {self.args.tool} | 断点续跑: {self.args.resume})")
        logger.info("=" * 60)
        
        self.setup_output_directory()
        self.build_index()
        self.find_samples()
        
        all_results =[]
        if self.args.parallel and self.args.threads > 1:
            max_parallel_jobs = min(max(1, self.args.threads // self.args.align_threads), len(self.samples))
            if self.args.parallel_jobs: max_parallel_jobs = min(self.args.parallel_jobs, len(self.samples))
            logger.info(f"🔀 开始并行处理任务 (最大并行数: {max_parallel_jobs})")
            
            with ProcessPoolExecutor(max_workers=max_parallel_jobs) as executor:
                future_to_sample = {executor.submit(self.process_sample, s): s for s in self.samples}
                with tqdm(total=len(self.samples), desc="处理进度", unit="样本", bar_format='{l_bar}{bar}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}, {postfix}, {rate_fmt}]') as pbar:
                    for future in as_completed(future_to_sample):
                        sample = future_to_sample[future]
                        pbar.set_postfix_str(f"处理: {sample['name']} 中")
                        res = future.result()
                        if res is not None: all_results.extend(res)
                        pbar.update(1)
        else:
            with tqdm(self.samples, desc="处理进度", unit="样本", bar_format='{l_bar}{bar}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}, {postfix}, {rate_fmt}]') as pbar:
                for sample in pbar:
                    pbar.set_postfix_str(f"处理: {sample['name']} 中")
                    res = self.process_sample(sample)
                    if res is not None: all_results.extend(res)
        
        if all_results:
            best_records, best_summary_file = self.save_results(all_results)
            
            if self.args.consensus and best_records:
                self.build_consensus_sequences(best_records)
                
            # 全自动 R 语言制图 (依赖于最佳核心文件)
            if best_summary_file:
                self.generate_plots(best_summary_file)
        else:
            logger.warning("⚠️ 未找到任何满足过滤条件的病毒记录")
        
        total_duration = datetime.now() - start_time
        logger.info("=" * 60)
        logger.info("✨ 分析管线全流程运行结束")
        logger.info(f"⏱️  管线总耗时: {Timer.format_duration(total_duration.total_seconds())}")
        logger.info("=" * 60)
    
    def save_results(self, all_results):
        with Timer("汇总去重与保存结果"):
            df_all = pd.DataFrame(all_results)
            total_raw_records = len(df_all)
            
            if 'Virus' in df_all.columns:
                found_viruses = df_all['Virus'].unique().tolist()
                df_all['description'] = "Unannotated"
                df_all['taxonomy'] = "Unannotated"
                df_all['taxonomy_id'] = "-"
                
                # a. RVDB 快速匹配
                if self.args.rvdb_taxon:
                    taxon_db = self.fast_load_rvdb_taxonomy(found_viruses)
                    if taxon_db:
                        def get_tax_info(vid, key):
                            if vid in taxon_db: return taxon_db[vid][key]
                            base_id = vid.split('.')[0]
                            if base_id in taxon_db: return taxon_db[base_id][key]
                            return "Unannotated" if key != "taxonomy_id" else "-"
                            
                        df_all['description'] = df_all['Virus'].apply(lambda x: get_tax_info(x, 'description'))
                        df_all['taxonomy'] = df_all['Virus'].apply(lambda x: get_tax_info(x, 'taxonomy'))
                        df_all['taxonomy_id'] = df_all['Virus'].apply(lambda x: get_tax_info(x, 'taxonomy_id'))

                # b. NCBI Entrez 补全漏网之鱼
                unannotated_viruses = df_all[df_all['taxonomy'] == 'Unannotated']['Virus'].unique().tolist()
                if unannotated_viruses:
                    if self.args.email:
                        ncbi_db = self.fetch_ncbi_taxonomy(unannotated_viruses)
                        if ncbi_db:
                            for vid in unannotated_viruses:
                                if vid in ncbi_db:
                                    mask = df_all['Virus'] == vid
                                    df_all.loc[mask, 'description'] = ncbi_db[vid]['description']
                                    df_all.loc[mask, 'taxonomy'] = ncbi_db[vid]['taxonomy']
                                    df_all.loc[mask, 'taxonomy_id'] = ncbi_db[vid]['taxonomy_id']

                cols = list(df_all.columns)
                for col in ['taxonomy_id', 'taxonomy', 'description']:
                    if col in cols: cols.remove(col)
                    cols.insert(cols.index('Virus') + 1, col)
                df_all = df_all[cols]

            summary_dir = self.output_dir / 'summary'
            summary_file = summary_dir / f'all_viruses.summary.{self.args.format}'
            
            if self.args.format == 'csv': df_all.to_csv(summary_file, index=False)
            else: df_all.to_csv(summary_file, sep='\t', index=False)
            
            # --- 【核心优化】：同源去除，提取最优核心序列 ---
            if 'taxonomy' in df_all.columns:
                df_all['group_tax'] = df_all.apply(lambda x: x['taxonomy'] if x['taxonomy'] not in['Unannotated', '-'] else x['Virus'], axis=1)
            else:
                df_all['group_tax'] = df_all['Virus']
                
            df_sorted = df_all.sort_values(by=['Sample', 'group_tax', 'Coverage(%)', 'MeanDepth'], ascending=[True, True, False, False])
            df_best = df_sorted.drop_duplicates(subset=['Sample', 'group_tax'], keep='first').copy()
            
            df_best.drop(columns=['group_tax'], inplace=True)
            df_all.drop(columns=['group_tax'], inplace=True)

            best_summary_file = summary_dir / f'all_viruses.best.summary.{self.args.format}'
            if self.args.format == 'csv': df_best.to_csv(best_summary_file, index=False)
            else: df_best.to_csv(best_summary_file, sep='\t', index=False)
            
            logger.info(f"✅ 全量冗余结果已保存 ({total_raw_records} 条记录)")
            logger.info(f"✅ 最优代表株结果已保存至 all_viruses.best.summary ({len(df_best)} 条核心记录)")
            
            self.generate_report(df_best, total_raw_records)
            return df_best.to_dict('records'), best_summary_file

    def generate_plots(self, best_summary_file):
        """利用内置 R 脚本对最优结果进行自动可视化"""
        logger.info("\n📊 开始调用 R 脚本生成全自动可视化图表...")
        plots_dir = self.output_dir / 'plots'
        r_script_path = plots_dir / 'virus_frequency_plot.R'
        
        with open(r_script_path, 'w', encoding='utf-8') as f:
            f.write(R_PLOT_SCRIPT)
            
        if shutil.which('Rscript') is None:
            logger.warning("⚠️ 未检测到 Rscript 环境，跳过绘图步骤。如果您需要绘图，请确保已安装 R。")
            return
            
        # 运行制图命令 (利用 --multi-plot 开启多维度拼图展示，并启用 log10 适应极差)
        cmd =[
            'Rscript', str(r_script_path),
            '-i', str(best_summary_file),
            '-o', str(plots_dir / 'virus_analysis'),
            '--multi-plot',
            '--log10-transform'
        ]
        
        try:
            res = subprocess.run(cmd, check=True, capture_output=True, text=True)
            logger.info(f"✅ 分析图表生成成功，已保存至: {plots_dir}/")
        except subprocess.CalledProcessError as e:
            logger.error(f"❌ R 脚本图表生成失败，Rscript报错:\n{e.stderr}")

    # ================= 共识序列处理核心模块 =================
    def extract_fastas_with_python(self, target_ids, out_dir):
        target_set = set(target_ids)
        found = set()
        current_id, current_seq = None,[]
        
        def save_record(vid, seq_list):
            if vid in target_set:
                vname = re.sub(r'[\\/*?:"<>|]', "_", vid)
                v_dir = out_dir / vname
                v_dir.mkdir(parents=True, exist_ok=True)
                with open(v_dir / f"{vname}.ref.fasta", 'w') as outf:
                    outf.write(f">{vid}\n" + "".join(seq_list) + "\n")
                found.add(vid)
                
        with open(self.args.reference, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    if current_id: save_record(current_id, current_seq)
                    current_id = line[1:].split()[0]
                    current_seq =[]
                else:
                    current_seq.append(line)
            if current_id: save_record(current_id, current_seq)

    def _run_single_consensus(self, task):
        sample, virus, mean_depth = task
        consensus_dir = self.output_dir / 'consensus'
        
        vname = re.sub(r'[\\/*?:"<>|]', "_", virus)
        virus_dir = consensus_dir / vname
        
        bam_file = self.output_dir / 'bam' / f"{sample}.sorted.bam"
        ref_fasta = virus_dir / f"{vname}.ref.fasta"
        fixed_bam = virus_dir / f"{sample}.{vname}.fixed.bam"
        out_fasta = virus_dir / f"{sample}.{vname}.consensus.fasta"
        
        # --- 断点续跑 Checkpoint 拦截 ---
        if self.args.resume and out_fasta.exists() and out_fasta.stat().st_size > 0:
            self.write_sample_log(sample, f"\n⏭️ [断点续跑] 跳过共识序列生成，因目标文件已存在: {out_fasta.name}")
            return True
        
        if not ref_fasta.exists() or not bam_file.exists(): return False
            
        min_depth = int(math.floor(mean_depth))
        if min_depth > 10: min_depth = 10
        elif min_depth < 1: min_depth = 1
        
        log_block = [f"\n[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] === 共识序列生成: 病毒 {virus} ==="]
        log_block.append(f"提取依据: 最优覆盖代表株 | 动态 min_depth 设置为: {min_depth}")
        
        try:
            awk_cmd = f"awk -v v='{virus}' '/^@SQ/ && $2 != \"SN:\"v {{next}} {{print}}'"
            pipe_cmd = f"samtools view -h '{bam_file}' '{virus}' | {awk_cmd} | samtools view -b -o '{fixed_bam}'"
            log_block.append(f"过滤BAM头部 CMD: {pipe_cmd}")
            subprocess.run(pipe_cmd, shell=True, check=True, stderr=subprocess.PIPE, text=True)
            
            vc_cmd = f"viral_consensus -i '{fixed_bam}' -r '{ref_fasta}' -o '{out_fasta}' --min_depth {min_depth}"
            log_block.append(f"viral_consensus CMD: {vc_cmd}")
            res = subprocess.run(vc_cmd, shell=True, check=True, capture_output=True, text=True)
            log_block.append(f"✅ 构建成功! 细节输出:\n{res.stderr.strip()}")
            
            if fixed_bam.exists(): fixed_bam.unlink()
            self.write_sample_log(sample, "\n".join(log_block))
            return True
        except subprocess.CalledProcessError as e:
            log_block.append(f"❌ 共识序列生成失败:\n{e.stderr}")
            self.write_sample_log(sample, "\n".join(log_block))
            return False

    def build_consensus_sequences(self, best_results):
        logger.info("\n🧬 开始并行生成各核心病毒的共识序列...")
        consensus_dir = self.output_dir / 'consensus'
        
        with Timer("核心共识序列构建"):
            unique_viruses = list(set(r['Virus'] for r in best_results))
            self.extract_fastas_with_python(unique_viruses, consensus_dir)
            
            tasks = [(r['Sample'], r['Virus'], r['MeanDepth']) for r in best_results]
            success_count = 0
            
            with ProcessPoolExecutor(max_workers=min(len(tasks), self.args.threads)) as executor:
                futures = {executor.submit(self._run_single_consensus, t): t for t in tasks}
                with tqdm(total=len(tasks), desc="构建进度", unit="文件", bar_format='{l_bar}{bar}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}, {postfix}, {rate_fmt}]') as pbar:
                    for future in as_completed(futures):
                        task = futures[future]
                        pbar.set_postfix_str(f"提取: {task[0]} ({task[1][:10]}..)中")
                        if future.result(): success_count += 1
                        pbar.update(1)
            
            logger.info(f"✅ 成功归档了 {success_count} 个核心共识序列文件")

    def generate_report(self, df_best, total_raw_records):
        report_file = self.output_dir / 'summary' / 'analysis_report.txt'
        with open(report_file, 'w') as f:
            f.write("=" * 80 + "\n")
            f.write(f"病毒丰度核心分析报告 (基于全测序文库的绝对定量 + 分类去冗余)\n")
            f.write("=" * 80 + "\n\n")
            
            f.write(f"分析时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"比对工具: {self.args.tool}\n")
            f.write(f"参考基因组: {self.args.reference}\n")
            if self.args.rvdb_taxon: f.write(f"本地注释文件: {self.args.rvdb_taxon}\n")
            if self.args.email: f.write(f"开启 NCBI API 补全 (绑定邮箱: {self.args.email})\n")
            f.write(f"基础过滤条件: Coverage > {self.args.coverage}%, MeanDepth > {self.args.meandepth}\n")
            f.write(f"定量过滤要求: FPKM >= {self.args.min_fpkm}, TPM >= {self.args.min_tpm}\n")
            f.write(f"共识提取状态: {'已针对核心代表株生成 (存于 consensus 下)' if self.args.consensus else '否'}\n\n")
            
            f.write(f"总处理样本数: {len(self.samples)}\n")
            f.write(f"初筛包含冗余记录总数: {total_raw_records} (保存在 all_viruses.summary)\n")
            f.write(f"同源去重后核心记录数: {len(df_best)} (保存在 all_viruses.best.summary)\n\n")
            
            if 'Virus' in df_best.columns:
                f.write("同源去重后的最终确诊核心病毒群落概览:\n")
                f.write("-" * 60 + "\n")
                virus_counts = df_best['Virus'].value_counts()
                for virus, count in virus_counts.items():
                    vd = df_best[df_best['Virus'] == virus]
                    desc = vd['description'].iloc[0] if 'description' in vd.columns else ''
                    mean_tpm, mean_fpkm, mean_rpm = vd['TPM'].mean(), vd['FPKM'].mean(), vd['RPM'].mean()
                    f.write(f"{virus} ({desc[:40]}...):\n")
                    f.write(f"  样本检出数: {count}\n")
                    f.write(f"  核心平均 TPM:  {mean_tpm:.2f} (该毒株在整个群落中的相对占比)\n")
                    f.write(f"  核心平均 FPKM: {mean_fpkm:.2f} (长度归一化后的绝对载量)\n")
                    f.write(f"  核心平均 RPM:  {mean_rpm:.2f} (无偏差的绝对测序文库载量)\n\n")
        
def main():
    parser = argparse.ArgumentParser(
        description='病毒分析终极版 (带断点续跑、最优毒株筛选、共识自动提取与 R制图)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
增强功能说明:
  --resume   : 开启断点续跑。将自动跳过已完成比对、深度计算和共识提取的文件，节约重复调试时间。
  --consensus: 基于去重后的同源最优代表序列自动生成共识FASTA，按病毒名分目录保存并智能截断阈值。
  全自动绘图 : 流程跑完后会自动将 R 脚本写入 plots 目录并执行，生成高分辨率的数据可视化箱线图。
        """
    )
    
    input_group = parser.add_argument_group('输入输出参数')
    input_group.add_argument('--input_dir', type=str, help='输入文件夹路径（包含fastq文件）')
    input_group.add_argument('--reference', type=str, required=True, help='参考基因组fasta文件路径')
    input_group.add_argument('--rvdb_taxon', type=str, help='RVDB 分类注释文件路径 (如: RVDB_Taxon_Current.tab)')
    input_group.add_argument('--email', type=str, help='[可选] 提供有效邮箱，用于启用 NCBI 联网获取未知病毒注释')
    input_group.add_argument('--api_key', type=str, help='[可选] NCBI API Key，可以提高联网获取的访问频率')
    input_group.add_argument('--output_dir', type=str, default='./virus_analysis', help='输出目录路径 [默认: ./virus_analysis]')
    input_group.add_argument('--sample_list', type=str, help='样本列表文件（每行: 样本名 fastq1 fastq2）')
    input_group.add_argument('--single_end', action='store_true', help='使用单端测序数据')
    
    align_group = parser.add_argument_group('高级流程与控制参数')
    align_group.add_argument('--tool', type=str, choices=['strobealign', 'minimap2', 'bwa', 'bwa-mem2', 'bowtie2', 'hisat2'], default='strobealign', help='比对工具选择 [默认: strobealign]')
    align_group.add_argument('--align_threads', type=int, default=8, help='每个样本比对的线程数[默认: 8]')
    align_group.add_argument('--consensus', action='store_true', help='基于去重后的最佳代表株自动生成归档的共识序列')
    align_group.add_argument('--resume', action='store_true', help='开启断点续跑模式，自动跳过已存在的结果文件')
    
    filter_group = parser.add_argument_group('过滤参数')
    filter_group.add_argument('--coverage', type=float, default=90.0, help='Coverage(%%)阈值[默认: 90]')
    filter_group.add_argument('--meandepth', type=float, default=10.0, help='MeanDepth阈值 [默认: 10]')
    filter_group.add_argument('--min_fpkm', type=float, default=0.0, help='FPKM最小值过滤 [默认: 0]')
    filter_group.add_argument('--min_tpm', type=float, default=0.0, help='TPM最小值过滤 [默认: 0]')
    
    parallel_group = parser.add_argument_group('并行处理参数')
    parallel_group.add_argument('--parallel', action='store_true', help='并行处理样本')
    parallel_group.add_argument('--parallel_jobs', type=int, help='并行任务数（如果不指定，自动计算）')
    parallel_group.add_argument('--threads', type=int, default=4, help='总线程数[默认: 4]')
    
    control_group = parser.add_argument_group('兼容控制参数')
    control_group.add_argument('--skip_alignment', action='store_true', help='跳过比对步骤 (若使用 --resume 则自动包含此功能)')
    control_group.add_argument('--skip_depth', action='store_true', help='跳过深度计算步骤 (若使用 --resume 则自动包含此功能)')
    control_group.add_argument('--format', type=str, choices=['csv', 'tsv'], default='tsv', help='输出文件格式[默认: tsv]')
    
    other_group = parser.add_argument_group('其他参数')
    other_group.add_argument('--verbose', action='store_true', help='输出详细日志信息')
    
    args = parser.parse_args()
    if not args.input_dir and not args.sample_list: parser.error("必须提供 --input_dir 或 --sample_list 参数")
    if not os.path.exists(args.reference): parser.error(f"参考基因组文件不存在: {args.reference}")
    if args.input_dir and not os.path.exists(args.input_dir): parser.error(f"输入目录不存在: {args.input_dir}")
    if args.verbose: logger.setLevel(logging.DEBUG)
    
    pipeline = VirusAnalysisPipeline(args)
    pipeline.run_pipeline()

if __name__ == '__main__':
    main()
