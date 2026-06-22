#!/usr/bin/env python3
"""
病毒序列比对、深度统计、多维度定量分析、共识序列构建与可视化完整流程
1. 包含 FPKM, RPM, TPM 三大定量体系
2. 包含 RVDB 极速流式解析 + NCBI Entrez 智能联网补全 (含实时缓存断点续跑)
3. 包含最优同源代表株筛选机制 (Best Representative by Taxonomy)
4. 包含阳性reads提取、变异检测及断点续跑(--resume)、全自动 R 绘图出图
5. 三级输出逻辑: raw (全量比对) -> summary (阈值过滤) -> best (同源去重最优)
6. [新增] 完美兼容 5 大软件的 Unique_Reads / Multi_Reads 分离统计算法
7. [新增] 样本级全局 BAM 测序比对质量评估 (基于 samtools stats)
8. [新增] 完美兼容 .fa/.fasta 输入及单双端数据混合自动发现
9. [修复] 强制清理分类名中的单引号/括号等特殊字符，防止引起shell断点报错
10. [优化] 过滤前置：仅对通过阈值过滤的有效病毒进行分类学和联网注释，极大提升速度
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

if ("taxonomy" %in% colnames(data)) {
  data$Display_Name <- ifelse(data$taxonomy != "Unannotated" & data$taxonomy != "-", paste0(data$taxonomy, "\n(", data$Virus, ")"), data$Virus)
} else {
  data$Display_Name <- data$Virus
}
data$Display_Name <- sapply(data$Display_Name, function(x) paste(strwrap(x, width = 40), collapse = "\n"))

required_columns <- list("MeanDepth" = c("Display_Name", "MeanDepth"), "FPKM" = c("Display_Name", "FPKM"), "RPM" = c("Display_Name", "RPM"), "TPM" = c("Display_Name", "TPM"))
modes <- if(opt$modes == "all") names(required_columns) else strsplit(trimws(opt$modes), ",")[[1]]
available_modes <- modes[sapply(modes, function(m) all(required_columns[[m]] %in% colnames(data)))]
if (length(available_modes) == 0) stop("错误: 没有可用的分析模式。")

output_dir <- dirname(opt$output)
if (output_dir != "." && !dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

unique_viruses <- length(unique(data$Display_Name))
should_flip <- unique_viruses > 5
theme_func <- switch(opt$theme, "classic"=theme_classic, "minimal"=theme_minimal, "bw"=theme_bw, theme_bw)

if (opt$`multi-plot` && length(available_modes) > 1) {
  plot_data_long <- data %>% select(Display_Name, Sample, all_of(available_modes)) %>%
    pivot_longer(cols = all_of(available_modes), names_to = "Metric", values_to = "Value") %>% filter(!is.na(Value))
  plot_data_long$Metric <- factor(plot_data_long$Metric, levels = available_modes)
  if (opt$`log10-transform`) {
    for (metric in unique(plot_data_long$Metric)) {
      md <- plot_data_long[plot_data_long$Metric == metric, ]
      min_p <- min(md$Value[md$Value > 0], na.rm=TRUE)
      plot_data_long$Value[plot_data_long$Metric == metric & plot_data_long$Value <= 0] <- min_p / 10
    }
  }
  
  base_metric <- available_modes[1]
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
        if exc_type is None: logger.info(f"✅ 完成: {self.name}[耗时: {self.format_duration(duration)}]")
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
        
        if self.args.resume:
            self.args.skip_alignment = True
            self.args.skip_depth = True

        self.check_tools()
        
    def check_tools(self):
        required_tools = {'samtools': '测序数据处理', 'pandepth': '深度计算', self.args.tool: '序列比对'}
        if self.args.tool == 'bowtie2': required_tools['bowtie2-build'] = 'Bowtie2索引构建'
        elif self.args.tool == 'hisat2': required_tools['hisat2-build'] = 'HISAT2索引构建'
        
        if getattr(self.args, 'extract_reads', False):
            required_tools['pigz'] = '多线程压缩工具'
            
        if self.args.consensus:
            required_tools['viral_consensus'] = '共识序列构建'
            required_tools['awk'] = 'BAM头处理工具'
            
        if getattr(self.args, 'call_variants', False):
            required_tools[self.args.variant_caller] = '变异检测工具'
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

        ncbi_dict = {}
        cache_file = self.output_dir / 'summary' / 'ncbi_taxonomy_cache.tsv'
        
        if self.args.resume and cache_file.exists():
            try:
                df_cache = pd.read_csv(cache_file, sep='\t', dtype=str)
                for _, row in df_cache.iterrows():
                    ncbi_dict[row['accession']] = {
                        'description': str(row['description']),
                        'taxonomy': str(row['taxonomy']),
                        'taxonomy_id': str(row['taxonomy_id'])
                    }
                logger.info(f"🔄 [断点续跑] 从本地缓存加载了 {len(ncbi_dict)} 条历史 NCBI 获取记录")
            except Exception as e:
                logger.warning(f"⚠️ 读取 NCBI 缓存失败，将重新获取: {e}")

        to_fetch =[acc for acc in accessions if acc not in ncbi_dict]
        
        if not to_fetch:
            logger.info("✅ 所有需要注释的序列 NCBI 信息已在缓存中找到，跳过联网获取。")
            return ncbi_dict

        Entrez.email = self.args.email
        if self.args.api_key:
            Entrez.api_key = self.args.api_key
            sleep_time = 0.11
        else:
            sleep_time = 0.35

        logger.info(f"🌐 正在通过 NCBI 联网获取 {len(to_fetch)} 个存活病毒的分类信息...")

        write_header = not cache_file.exists()
        with open(cache_file, 'a', encoding='utf-8') as cf:
            if write_header:
                cf.write("accession\tdescription\ttaxonomy\ttaxonomy_id\n")
                
            for acc in tqdm(to_fetch, desc="NCBI 解析进度", unit="seq"):
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
                                
                        desc = record.description
                        tax = record.annotations.get('organism', '-')
                        
                        ncbi_dict[acc] = {'description': desc, 'taxonomy': tax, 'taxonomy_id': tax_id}
                        cf.write(f"{acc}\t{desc}\t{tax}\t{tax_id}\n")
                        cf.flush()
                        success = True
                        break
                    except Exception as e:
                        time.sleep(1)
                        
                if not success:
                    ncbi_dict[acc] = {'description': "Unannotated (NCBI Error)", 'taxonomy': "Unannotated", 'taxonomy_id': "-"}
                    cf.write(f"{acc}\tUnannotated (NCBI Error)\tUnannotated\t-\n")
                    cf.flush()
                    
                time.sleep(sleep_time)

        return ncbi_dict
    
    def setup_output_directory(self):
        subdirs =['bam', 'stat', 'summary', 'logs', 'index', 'plots']
        if getattr(self.args, 'extract_reads', False): subdirs.append('reads')
        if self.args.consensus: subdirs.append('consensus')
        if getattr(self.args, 'call_variants', False): subdirs.append('variants')
            
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
        if tool in['strobealign', 'minimap2']:
            self.index_path = str(ref_path)
            return
        prefix = self.output_dir / 'index' / f"{ref_path.stem}_{tool}"
        self.index_path = str(prefix)
        
        if tool == 'bwa':
            if Path(f"{prefix}.bwt").exists(): return
            cmd =['bwa', 'index', '-p', str(prefix), str(ref_path)]
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
        valid_exts = ['.fq', '.fastq', '.fq.gz', '.fastq.gz', '.fa', '.fasta', '.fa.gz', '.fasta.gz']
        samples = []
        
        if input_dir:
            all_files = []
            for ext in valid_exts:
                all_files.extend(input_dir.glob(f'*{ext}'))
                
            file_names = {f.name: f for f in all_files}
            processed_files = set()
            
            for f in all_files:
                if str(f) in processed_files: continue
                fname = f.name
                
                if self.args.single_end:
                    r2_name = None
                else:
                    r2_name = None
                    if '_1.' in fname: r2_name = fname.replace('_1.', '_2.')
                    elif '_R1' in fname: r2_name = fname.replace('_R1', '_R2')
                    elif '.R1.' in fname: r2_name = fname.replace('.R1.', '.R2.')
                    elif '.1.' in fname: r2_name = fname.replace('.1.', '.2.')
                    elif '_1_' in fname: r2_name = fname.replace('_1_', '_2_')
                    elif '.unmapped.R1.' in fname: r2_name = fname.replace('.unmapped.R1.', '.unmapped.R2.')
                    elif 'unmapped.R1.' in fname: r2_name = fname.replace('unmapped.R1.', 'unmapped.R2.')
                    elif '_1_unmapped.' in fname: r2_name = fname.replace('_1_unmapped.', '_2_unmapped.')
                
                # 跳过可能是由于R1衍生出的R2文件，避免其独立成单端
                if not self.args.single_end and any(p in fname for p in ['_2.', '_R2', '.R2.', '.2.', '_2_', '.unmapped.R2.', 'unmapped.R2.', '_2_unmapped.']):
                    r1_name_guess = None
                    if '_2.' in fname: r1_name_guess = fname.replace('_2.', '_1.')
                    elif '_R2' in fname: r1_name_guess = fname.replace('_R2', '_R1')
                    elif '.R2.' in fname: r1_name_guess = fname.replace('.R2.', '.R1.')
                    elif '.2.' in fname: r1_name_guess = fname.replace('.2.', '.1.')
                    elif '_2_' in fname: r1_name_guess = fname.replace('_2_', '_1_')
                    elif '.unmapped.R2.' in fname: r1_name_guess = fname.replace('.unmapped.R2.', '.unmapped.R1.')
                    elif 'unmapped.R2.' in fname: r1_name_guess = fname.replace('unmapped.R2.', 'unmapped.R1.')
                    elif '_2_unmapped.' in fname: r1_name_guess = fname.replace('_2_unmapped.', '_1_unmapped.')
                    
                    if r1_name_guess and r1_name_guess in file_names: 
                        continue # 交给匹配到 R1 的逻辑中一并处理
                
                # 获取干净的 sample name
                sname = fname
                for ext in valid_exts:
                    if sname.endswith(ext):
                        sname = sname[:-len(ext)]
                        break
                        
                for suf in ['_1', '_R1', '.R1', '.1', '_unmapped', '.unmapped', '_1_unmapped', '_2', '_R2', '.R2', '.2']:
                    if sname.endswith(suf): 
                        sname = sname[:-len(suf)]
                        break
                
                if r2_name and r2_name in file_names:
                    samples.append({'name': sname, 'r1': str(f), 'r2': str(file_names[r2_name])})
                    processed_files.add(str(f))
                    processed_files.add(str(file_names[r2_name]))
                else:
                    samples.append({'name': sname, 'r1': str(f), 'r2': None})
                    processed_files.add(str(f))

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
        align_cmd =[]

        # 智能识别 FASTA 输入，动态控制底层参数
        r1_is_fa = any(sample['r1'].lower().endswith(ext) for ext in ['.fa', '.fasta', '.fa.gz', '.fasta.gz'])

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
            if r1_is_fa: align_cmd.append('-f') # 支持FASTA
            if sample['r2']: align_cmd.extend(['-1', sample['r1'], '-2', sample['r2']])
            else: align_cmd.extend(['-U', sample['r1']])
        elif tool == 'hisat2':
            align_cmd =['hisat2', '-p', str(self.args.align_threads), '-x', self.index_path]
            if r1_is_fa: align_cmd.append('-f') # 支持FASTA
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

    # ================= [新增] 样本级 BAM 全局测序比对质量评估 =================
    def get_global_bam_qc(self, bam_file, sample_name):
        qc_file = self.output_dir / 'stat' / f"{sample_name}.global_qc.txt"
        
        if self.args.resume and qc_file.exists():
            self.write_sample_log(sample_name, f"⏭️  [断点续跑] 跳过全局 BAM 统计，直接使用已有文件: {qc_file.name}")
            return qc_file
            
        self.write_sample_log(sample_name, "\n--- 2.1 全局比对质量诊断 (samtools stats) ---")
        cmd = f"samtools stats '{bam_file}' | grep ^SN | cut -f 2- > '{qc_file}'"
        
        try:
            subprocess.run(cmd, shell=True, check=True, executable='/bin/bash')
            # 提取核心指标写入日志
            with open(qc_file, 'r') as f:
                lines = f.readlines()
                total, mapped, mq0, unmapped = 0, 0, 0, 0
                for line in lines:
                    if 'raw total sequences' in line: total = line.split('\t')[1].strip()
                    elif line.startswith('reads mapped:'): mapped = line.split('\t')[1].strip()
                    elif line.startswith('reads unmapped:'): unmapped = line.split('\t')[1].strip()
                    elif line.startswith('reads MQ0:'): mq0 = line.split('\t')[1].strip()
                
                log_msg = (
                    f"✅ 全局质控完成:\n"
                    f"  - 总测序 Reads  : {total}\n"
                    f"  - 未比对 (Unmapped) : {unmapped}\n"
                    f"  - 总比对 (Mapped)   : {mapped}\n"
                    f"  - 其中多重比对 (MQ0): {mq0} (占总比对的 {float(mq0)/float(mapped)*100:.1f}% 如果 Mapped>0)\n"
                ) if int(mapped) > 0 else f"✅ 全局质控完成: 0 条 reads 成功比对。"
                self.write_sample_log(sample_name, log_msg)
            return qc_file
        except Exception as e:
            self.write_sample_log(sample_name, f"❌ 全局 BAM 统计失败:\n{e}")
            return None

    def get_idxstats(self, bam_file, sample_name):
        self.write_sample_log(sample_name, "\n--- 2.2 提取测序绝对文库大小 (samtools idxstats) ---")
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

    # ================= [终极优化版] 提取 Unique / Multi Reads 统计 =================
    def get_unique_multi_stats(self, bam_file, sample_name):
        stat_out = self.output_dir / 'stat' / f"{sample_name}.um_stat.tsv"
        
        if self.args.resume and stat_out.exists():
            self.write_sample_log(sample_name, f"⏭️  [断点续跑] 跳过 Unique/Multi 解析，直接使用已有文件: {stat_out.name}")
            stats = {}
            with open(stat_out, 'r') as f:
                for line in f:
                    parts = line.strip().split('\t')
                    if len(parts) == 3:
                        stats[parts[0]] = {'unique': int(parts[1]), 'multi': int(parts[2])}
            return stats
            
        self.write_sample_log(sample_name, "\n--- 2.5 提取唯一/多重比对精细统计 (Unique/Multi-mapped) ---")
        
        awk_script = """
        {
            ref = $3;
            if (ref == "*") next; 
            
            # 使用行业金标准: MAPQ < 10 为多重比对
            if ($5 < 10) {
                multi[ref]++
            } else {
                uniq[ref]++
            }
        }
        END {
            for (r in uniq) {
                m = multi[r]+0;
                print r "\\t" uniq[r] "\\t" m;
                delete multi[r];
            }
            for (r in multi) {
                print r "\\t0\\t" multi[r];
            }
        }
        """
        
        cmd = f"samtools view -F 2308 '{bam_file}' | awk '{awk_script}' > '{stat_out}'"
        
        try:
            subprocess.run(cmd, shell=True, check=True, executable='/bin/bash')
            self.write_sample_log(sample_name, f"✅ Unique/Multi 统计完成 (基于 MAPQ 标准)。")
        except subprocess.CalledProcessError as e:
            self.write_sample_log(sample_name, f"❌ Unique/Multi 统计失败:\n{e}")
            return {}

        stats = {}
        with open(stat_out, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) == 3:
                    stats[parts[0]] = {'unique': int(parts[1]), 'multi': int(parts[2])}
        return stats

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
        um_stats = idx_data.get('um_stats', {}) 
        
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
                    
                    if mapped_reads == 0: 
                        continue 
                    
                    # 提取该病毒的 Unique 与 Multi reads 数
                    virus_um = um_stats.get(virus_name, {'unique': 0, 'multi': 0})
                    unique_reads = virus_um['unique']
                    multi_reads = virus_um['multi']
                    
                    rpm = (mapped_reads * 1e6) / global_total_reads if global_total_reads > 0 else 0.0
                    fpkm = (mapped_reads * 1e9) / (global_total_reads * length_val) if global_total_reads > 0 and length_val > 0 else 0.0
                    tpm = ((mapped_reads * 1000.0) / length_val * 1e6) / total_rpk if total_rpk > 0 and length_val > 0 else 0.0
                        
                    results.append({
                        'Sample': sample_name, 'Virus': virus_name, 'Length': length_val,
                        'CoveredSite': covered_site_val, 'TotalDepth': total_depth_val,
                        'Coverage(%)': coverage_val, 'MeanDepth': mean_depth_val,
                        'MappedReads': mapped_reads,
                        'Unique_Reads': unique_reads,  
                        'Multi_Reads': multi_reads,    
                        'RPM': round(rpm, 2),
                        'FPKM': round(fpkm, 2), 'TPM': round(tpm, 2)
                    })
            self.write_sample_log(sample_name, f"✅ 数据解析完毕, 提取出 {len(results)} 条含比稳记录的 Raw 数据。")
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
            
            # [新增] 样本级全局 BAM 质控
            self.get_global_bam_qc(bam_file, sample_name)
            
            idx_data = self.get_idxstats(bam_file, sample_name)
            if not idx_data: return None
            
            # 调用 Unique/Multi 解析模块
            um_stats = self.get_unique_multi_stats(bam_file, sample_name)
            if um_stats: idx_data['um_stats'] = um_stats
            
            stat_file = self.run_pandepth(bam_file, sample_name)
            if not stat_file: return None
            
            res = self.parse_stat_file(stat_file, sample_name, idx_data)
            
            end_time = datetime.now()
            duration = end_time - start_time
            self.write_sample_log(sample_name, f"[{end_time.strftime('%Y-%m-%d %H:%M:%S')}] 🎉 样本整体处理流程结束[耗时: {Timer.format_duration(duration.total_seconds())}]")
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
                with tqdm(total=len(self.samples), desc="处理进度", unit="样本", bar_format='{l_bar}{bar}| {n_fmt}/{total_fmt}[{elapsed}<{remaining}, {postfix}, {rate_fmt}]') as pbar:
                    for future in as_completed(future_to_sample):
                        sample = future_to_sample[future]
                        pbar.set_postfix_str(f"处理: {sample['name']} 中")
                        res = future.result()
                        if res is not None: all_results.extend(res)
                        pbar.update(1)
        else:
            with tqdm(self.samples, desc="处理进度", unit="样本", bar_format='{l_bar}{bar}| {n_fmt}/{total_fmt}[{elapsed}<{remaining}, {postfix}, {rate_fmt}]') as pbar:
                for sample in pbar:
                    pbar.set_postfix_str(f"处理: {sample['name']} 中")
                    res = self.process_sample(sample)
                    if res is not None: all_results.extend(res)
        
        best_summary_file = None
        if all_results:
            best_records, best_summary_file = self.save_results(all_results)
            
            # --- 分支0: Reads提取 ---
            if getattr(self.args, 'extract_reads', False) and best_records:
                self.extract_mapped_reads(best_records)
                
            # --- 分支1: 变异检测 ---
            if getattr(self.args, 'call_variants', False) and best_records:
                self.run_variants_calling(best_records)
                
            # --- 分支2: 共识提取 ---
            if getattr(self.args, 'consensus', False) and best_records:
                self.build_consensus_sequences(best_records)
                
            # --- 出图 ---
            if best_summary_file:
                self.generate_plots(best_summary_file)
        else:
            logger.warning("⚠️ 未找到任何有 MappedReads 的病毒记录")
        
        total_duration = datetime.now() - start_time
        logger.info("=" * 60)
        logger.info("✨ 分析管线全流程运行结束")
        logger.info(f"⏱️  管线总耗时: {Timer.format_duration(total_duration.total_seconds())}")
        
        if best_summary_file:
            self.print_plotting_instructions(best_summary_file)
        
        logger.info("=" * 60)

    def save_results(self, all_results):
        with Timer("汇总去重与保存结果"):
            df_all = pd.DataFrame(all_results)
            summary_dir = self.output_dir / 'summary'
            
            # 【优化修改点 1】: 将全量未注释数据直接保存为 raw 文件（不再进行几千条无关序列的冗余注释）
            raw_summary_file = summary_dir / f'all_viruses.raw.{self.args.format}'
            if self.args.format == 'csv': df_all.to_csv(raw_summary_file, index=False)
            else: df_all.to_csv(raw_summary_file, sep='\t', index=False)
            total_raw_records = len(df_all)
            
            # 【优化修改点 2】: 提前进行阈值过滤
            df_filtered = df_all[
                (df_all['Coverage(%)'] > self.args.coverage) &
                (df_all['MeanDepth'] > self.args.meandepth) &
                (df_all['FPKM'] >= self.args.min_fpkm) &
                (df_all['TPM'] >= self.args.min_tpm)
            ].copy()
            total_filtered_records = len(df_filtered)

            # 【优化修改点 3】: 仅对存活在 summary 中的病毒提取分类信息
            if not df_filtered.empty and 'Virus' in df_filtered.columns:
                # 核心逻辑变更：只提取过滤后剩下的目标病毒用于注释查询
                found_viruses = df_filtered['Virus'].unique().tolist()
                
                df_filtered['description'] = "Unannotated"
                df_filtered['taxonomy'] = "Unannotated"
                df_filtered['taxonomy_id'] = "-"
                
                if self.args.rvdb_taxon:
                    taxon_db = self.fast_load_rvdb_taxonomy(found_viruses)
                    if taxon_db:
                        def get_tax_info(vid, key):
                            if vid in taxon_db: return taxon_db[vid][key]
                            base_id = vid.split('.')[0]
                            if base_id in taxon_db: return taxon_db[base_id][key]
                            return "Unannotated" if key != "taxonomy_id" else "-"
                            
                        df_filtered['description'] = df_filtered['Virus'].apply(lambda x: get_tax_info(x, 'description'))
                        df_filtered['taxonomy'] = df_filtered['Virus'].apply(lambda x: get_tax_info(x, 'taxonomy'))
                        df_filtered['taxonomy_id'] = df_filtered['Virus'].apply(lambda x: get_tax_info(x, 'taxonomy_id'))

                unannotated_viruses = df_filtered[df_filtered['taxonomy'] == 'Unannotated']['Virus'].unique().tolist()
                if unannotated_viruses and self.args.email:
                    # 现在 NCBI API 只会收到过滤后少量需要查询的 accession id
                    ncbi_db = self.fetch_ncbi_taxonomy(unannotated_viruses)
                    if ncbi_db:
                        for vid in unannotated_viruses:
                            if vid in ncbi_db:
                                mask = df_filtered['Virus'] == vid
                                df_filtered.loc[mask, 'description'] = ncbi_db[vid]['description']
                                df_filtered.loc[mask, 'taxonomy'] = ncbi_db[vid]['taxonomy']
                                df_filtered.loc[mask, 'taxonomy_id'] = ncbi_db[vid]['taxonomy_id']

                cols = list(df_filtered.columns)
                for col in ['taxonomy_id', 'taxonomy', 'description']:
                    if col in cols: cols.remove(col)
                    cols.insert(cols.index('Virus') + 1, col)
                df_filtered = df_filtered[cols]
                
            # 保存注入注释信息后的 Summary 文件
            filtered_summary_file = summary_dir / f'all_viruses.summary.{self.args.format}'
            if self.args.format == 'csv': df_filtered.to_csv(filtered_summary_file, index=False)
            else: df_filtered.to_csv(filtered_summary_file, sep='\t', index=False)

            # --- 生成 Best 数据集 ---
            if not df_filtered.empty and 'taxonomy' in df_filtered.columns:
                df_filtered['group_tax'] = df_filtered.apply(lambda x: x['taxonomy'] if x['taxonomy'] not in ['Unannotated', '-'] else x['Virus'], axis=1)
            else:
                df_filtered['group_tax'] = df_filtered['Virus'] if not df_filtered.empty else []
                
            if not df_filtered.empty:
                df_sorted = df_filtered.sort_values(by=['Sample', 'group_tax', 'Coverage(%)', 'MeanDepth'], ascending=[True, True, False, False])
                df_best = df_sorted.drop_duplicates(subset=['Sample', 'group_tax'], keep='first').copy()
                df_best.drop(columns=['group_tax'], inplace=True)
            else:
                df_best = pd.DataFrame()

            best_summary_file = summary_dir / f'all_viruses.best.summary.{self.args.format}'
            if self.args.format == 'csv': df_best.to_csv(best_summary_file, index=False)
            else: df_best.to_csv(best_summary_file, sep='\t', index=False)
            
            logger.info(f"✅ 全量 Mapped 记录已保存至 all_viruses.raw ({total_raw_records} 条, 不含冗余注释)")
            logger.info(f"✅ 阈值过滤记录已保存至 all_viruses.summary ({total_filtered_records} 条, 已注释)")
            logger.info(f"✅ 最优代表株结果已保存至 all_viruses.best.summary ({len(df_best)} 条核心记录)")
            
            self.generate_report(df_best, total_raw_records, total_filtered_records)
            return df_best.to_dict('records'), best_summary_file

    def generate_plots(self, best_summary_file):
        logger.info("\n📊 开始调用 R 脚本生成全自动可视化图表...")
        plots_dir = self.output_dir / 'plots'
        r_script_path = plots_dir / 'virus_frequency_plot.R'
        
        with open(r_script_path, 'w', encoding='utf-8') as f:
            f.write(R_PLOT_SCRIPT)
            
        if shutil.which('Rscript') is None:
            logger.warning("⚠️ 未检测到 Rscript 环境，跳过绘图步骤。")
            return
            
        cmd =[
            'Rscript', str(r_script_path),
            '-i', str(best_summary_file),
            '-o', str(plots_dir / 'virus_analysis'),
            '--multi-plot',
            '--log10-transform'
        ]
        
        try:
            subprocess.run(cmd, check=True, capture_output=True, text=True)
            logger.info(f"✅ 分析图表生成成功，已保存至: {plots_dir}/")
        except subprocess.CalledProcessError as e:
            logger.error(f"❌ R 脚本图表生成失败:\n{e.stderr}")

    def extract_fastas_with_python(self, target_ids, virus_to_tax, out_dir):
        target_set = set(target_ids)
        found = set()
        current_id, current_seq = None,[]
        
        def save_record(vid, seq_list):
            if vid in target_set:
                # 核心修复: 强力清洗一切引起shell断裂的符号
                vname = re.sub(r"[\\/*?:\"<>| '\(\)\[\]{}]", "_", vid)
                tax = virus_to_tax.get(vid, 'Unannotated')
                safe_tax = re.sub(r"[\\/*?:\"<>| '\(\)\[\]{}]", "_", tax)
                if len(safe_tax) > 50: safe_tax = safe_tax[:50]
                
                folder_name = f"{safe_tax}_{vname}"
                v_dir = out_dir / folder_name
                v_dir.mkdir(parents=True, exist_ok=True)
                
                with open(v_dir / f"{folder_name}.ref.fasta", 'w') as outf:
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

    def _run_single_read_extraction(self, sample):
        bam_file = self.output_dir / 'bam' / f"{sample}.sorted.bam"
        reads_dir = self.output_dir / 'reads'
        
        r1 = reads_dir / f"{sample}_virus_1.fastq.gz"
        r2 = reads_dir / f"{sample}_virus_2.fastq.gz"
        rs = reads_dir / f"{sample}_virus_single.fastq.gz"
        
        if self.args.resume and (r1.exists() or r2.exists() or rs.exists()):
            self.write_sample_log(sample, f"\n⏭️ [断点续跑] 跳过 reads 提取，目标文件已存在: {reads_dir.name}/")
            return True
            
        if not bam_file.exists(): return False

        log_block = [f"\n[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] === 提取比对上病毒的序列 Reads ==="]
        try:
            cmd = (
                f"samtools fastq -F 4 "
                f"-1 >(pigz -c > '{r1}') "
                f"-2 >(pigz -c > '{r2}') "
                f"'{bam_file}' | pigz -c > '{rs}'"
            )
            subprocess.run(cmd, shell=True, check=True, executable='/bin/bash', stderr=subprocess.PIPE)
            
            for f in[r1, r2, rs]:
                if f.exists() and f.stat().st_size <= 50:
                    f.unlink()
                    
            log_block.append(f"✅ {sample} 病毒 reads 提取成功! 存入目录: reads/")
            self.write_sample_log(sample, "\n".join(log_block))
            return True
        except subprocess.CalledProcessError as e:
            err_msg = e.stderr.decode() if isinstance(e.stderr, bytes) else (e.stderr or str(e))
            log_block.append(f"❌ reads 提取失败:\n{err_msg}")
            self.write_sample_log(sample, "\n".join(log_block))
            return False

    def extract_mapped_reads(self, best_results):
        logger.info("\n📦 开始并行提取比对上病毒的序列 Reads (仅限阳性样本)...")
        reads_dir = self.output_dir / 'reads'
        reads_dir.mkdir(parents=True, exist_ok=True)
        
        unique_samples = list(set(r['Sample'] for r in best_results))
        success_count = 0
        
        with Timer("阳性病毒 Reads 提取"):
            with ProcessPoolExecutor(max_workers=min(len(unique_samples), self.args.threads)) as executor:
                futures = {executor.submit(self._run_single_read_extraction, s): s for s in unique_samples}
                with tqdm(total=len(unique_samples), desc="Reads提取", unit="样本", bar_format='{l_bar}{bar}| {n_fmt}/{total_fmt}[{elapsed}<{remaining}, {postfix}, {rate_fmt}]') as pbar:
                    for future in as_completed(futures):
                        sample = futures[future]
                        pbar.set_postfix_str(f"{sample} 中")
                        if future.result(): success_count += 1
                        pbar.update(1)
            
            logger.info(f"✅ 成功提取了 {success_count} 个阳性样本的病毒序列 Reads")

    def _run_single_variant_calling(self, task):
        sample, virus, taxonomy = task
        variants_dir = self.output_dir / 'variants'
        
        # 核心修复: 强力清洗一切引起shell断裂的符号
        vname = re.sub(r"[\\/*?:\"<>| '\(\)\[\]{}]", "_", virus)
        safe_tax = re.sub(r"[\\/*?:\"<>| '\(\)\[\]{}]", "_", taxonomy)
        if len(safe_tax) > 50: safe_tax = safe_tax[:50]
        
        folder_name = f"{safe_tax}_{vname}"
        virus_dir = variants_dir / folder_name
        
        bam_file = self.output_dir / 'bam' / f"{sample}.sorted.bam"
        ref_fasta = virus_dir / f"{folder_name}.ref.fasta"
        fixed_bam = virus_dir / f"{sample}.{folder_name}.fixed.bam"
        
        caller = self.args.variant_caller
        out_ext = 'tsv' if caller == 'ivar' else 'vcf'
        out_file = virus_dir / f"{sample}.{folder_name}.variants.{out_ext}"
        
        if self.args.resume and out_file.exists() and out_file.stat().st_size > 0:
            self.write_sample_log(sample, f"\n⏭️ [断点续跑] 跳过变异检测，目标文件已存在: {out_file.name}")
            return True
            
        if not ref_fasta.exists() or not bam_file.exists(): return False
            
        log_block = [f"\n[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] === 变异检测 ({caller}): {taxonomy} ({virus}) ==="]
        
        try:
            if not Path(str(ref_fasta) + '.fai').exists():
                subprocess.run(f"samtools faidx '{ref_fasta}'", shell=True, check=True, stderr=subprocess.PIPE)
                
            awk_cmd = f"awk -v v='{virus}' '/^@SQ/ && $2 != \"SN:\"v {{next}} {{print}}'"
            pipe_cmd = f"samtools view -h '{bam_file}' '{virus}' | {awk_cmd} | samtools view -b -o '{fixed_bam}'"
            subprocess.run(pipe_cmd, shell=True, check=True, stderr=subprocess.PIPE)
            subprocess.run(f"samtools index '{fixed_bam}'", shell=True, check=True, stderr=subprocess.PIPE)
            
            if caller == 'freebayes':
                cmd = f"freebayes -p 1 -f '{ref_fasta}' '{fixed_bam}' > '{out_file}'"
                subprocess.run(cmd, shell=True, check=True, stderr=subprocess.PIPE)
                
            elif caller == 'ivar':
                pileup = virus_dir / f"{sample}.pileup.txt"
                prefix = str(out_file).replace('.tsv', '')
                gff_param = f"-g '{self.args.ref_gff}'" if self.args.ref_gff else ""
                cmd1 = f"samtools mpileup -A -aa -d 0 -Q 0 --reference '{ref_fasta}' '{fixed_bam}' > '{pileup}'"
                cmd2 = f"cat '{pileup}' | ivar variants -r '{ref_fasta}' {gff_param} -p '{prefix}'"
                subprocess.run(cmd1, shell=True, check=True, stderr=subprocess.PIPE)
                subprocess.run(cmd2, shell=True, check=True, stderr=subprocess.PIPE)
                if pileup.exists(): pileup.unlink()
                
            elif caller == 'lofreq':
                tmp_bam = virus_dir / f"{sample}.tmp.lofreq.bam"
                cmd1 = f"lofreq indelqual -f '{ref_fasta}' --dindel -o '{tmp_bam}' '{fixed_bam}'"
                cmd2 = f"lofreq call -f '{ref_fasta}' --call-indels -o '{out_file}' '{tmp_bam}'"
                subprocess.run(cmd1, shell=True, check=True, stderr=subprocess.PIPE)
                subprocess.run(f"samtools index '{tmp_bam}'", shell=True, check=True, stderr=subprocess.PIPE)
                subprocess.run(cmd2, shell=True, check=True, stderr=subprocess.PIPE)
                if tmp_bam.exists(): tmp_bam.unlink()
                if Path(str(tmp_bam) + '.bai').exists(): Path(str(tmp_bam) + '.bai').unlink()
                
            log_block.append(f"✅ {sample} 变异检测 ({caller}) 成功! 存入目录: {folder_name}/")
            
            if fixed_bam.exists(): fixed_bam.unlink()
            if Path(str(fixed_bam) + '.bai').exists(): Path(str(fixed_bam) + '.bai').unlink()
            
            self.write_sample_log(sample, "\n".join(log_block))
            return True
        except subprocess.CalledProcessError as e:
            err_msg = e.stderr if hasattr(e, 'stderr') and e.stderr else str(e)
            log_block.append(f"❌ 变异检测失败:\n{err_msg}")
            self.write_sample_log(sample, "\n".join(log_block))
            return False

    def run_variants_calling(self, best_results):
        logger.info(f"\n🔬 开始并行执行变异检测 (工具: {self.args.variant_caller}, 按 Taxonomy 归档)...")
        variants_dir = self.output_dir / 'variants'
        
        virus_to_tax = {r['Virus']: r.get('taxonomy', 'Unannotated') for r in best_results}
        
        with Timer("核心代表株变异检测"):
            unique_viruses = list(set(r['Virus'] for r in best_results))
            self.extract_fastas_with_python(unique_viruses, virus_to_tax, variants_dir)
            
            tasks = [(r['Sample'], r['Virus'], r.get('taxonomy', 'Unannotated')) for r in best_results]
            success_count = 0
            
            with ProcessPoolExecutor(max_workers=min(len(tasks), self.args.threads)) as executor:
                futures = {executor.submit(self._run_single_variant_calling, t): t for t in tasks}
                with tqdm(total=len(tasks), desc=f"变异检测", unit="任务", bar_format='{l_bar}{bar}| {n_fmt}/{total_fmt}[{elapsed}<{remaining}, {postfix}, {rate_fmt}]') as pbar:
                    for future in as_completed(futures):
                        task = futures[future]
                        short_tax = str(task[2])[:10]
                        pbar.set_postfix_str(f"{task[0]} ({short_tax}..)中")
                        if future.result(): success_count += 1
                        pbar.update(1)
            
            logger.info(f"✅ 成功提取了 {success_count} 个按分类学命名的核心病毒变异检测结果 (VCF/TSV)")

    def _run_single_consensus(self, task):
        sample, virus, mean_depth, taxonomy = task
        consensus_dir = self.output_dir / 'consensus'
        
        # 核心修复: 强力清洗一切引起shell断裂的符号
        vname = re.sub(r"[\\/*?:\"<>| '\(\)\[\]{}]", "_", virus)
        safe_tax = re.sub(r"[\\/*?:\"<>| '\(\)\[\]{}]", "_", taxonomy)
        if len(safe_tax) > 50: safe_tax = safe_tax[:50]
        
        folder_name = f"{safe_tax}_{vname}"
        virus_dir = consensus_dir / folder_name
        
        bam_file = self.output_dir / 'bam' / f"{sample}.sorted.bam"
        ref_fasta = virus_dir / f"{folder_name}.ref.fasta"
        fixed_bam = virus_dir / f"{sample}.{folder_name}.fixed.bam"
        out_fasta = virus_dir / f"{sample}.{folder_name}.consensus.fasta"
        
        if self.args.resume and out_fasta.exists() and out_fasta.stat().st_size > 0:
            self.write_sample_log(sample, f"\n⏭️ [断点续跑] 跳过共识提取，目标文件已存在: {out_fasta.name}")
            return True
        
        if not ref_fasta.exists() or not bam_file.exists(): return False
            
        min_depth = int(math.floor(mean_depth))
        if min_depth > 10: min_depth = 10
        elif min_depth < 1: min_depth = 1
        
        log_block = [f"\n[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] === 共识生成: {taxonomy} ({virus}) ==="]
        
        try:
            awk_cmd = f"awk -v v='{virus}' '/^@SQ/ && $2 != \"SN:\"v {{next}} {{print}}'"
            pipe_cmd = f"samtools view -h '{bam_file}' '{virus}' | {awk_cmd} | samtools view -b -o '{fixed_bam}'"
            subprocess.run(pipe_cmd, shell=True, check=True, stderr=subprocess.PIPE, text=True)
            
            vc_cmd = f"viral_consensus -i '{fixed_bam}' -r '{ref_fasta}' -o '{out_fasta}' --min_depth 6"
            subprocess.run(vc_cmd, shell=True, check=True, capture_output=True, text=True)
            log_block.append(f"✅ {sample} 共识序列构建成功! 存入目录: {folder_name}/")
            
            if fixed_bam.exists(): fixed_bam.unlink()
            self.write_sample_log(sample, "\n".join(log_block))
            return True
        except subprocess.CalledProcessError as e:
            log_block.append(f"❌ 共识序列生成失败:\n{e.stderr}")
            self.write_sample_log(sample, "\n".join(log_block))
            return False

    def build_consensus_sequences(self, best_results):
        logger.info("\n🧬 开始并行生成各核心病毒的共识序列 (按 Taxonomy 归档)...")
        consensus_dir = self.output_dir / 'consensus'
        
        virus_to_tax = {r['Virus']: r.get('taxonomy', 'Unannotated') for r in best_results}
        
        with Timer("核心共识序列构建"):
            unique_viruses = list(set(r['Virus'] for r in best_results))
            self.extract_fastas_with_python(unique_viruses, virus_to_tax, consensus_dir)
            
            tasks = [(r['Sample'], r['Virus'], r['MeanDepth'], r.get('taxonomy', 'Unannotated')) for r in best_results]
            success_count = 0
            
            with ProcessPoolExecutor(max_workers=min(len(tasks), self.args.threads)) as executor:
                futures = {executor.submit(self._run_single_consensus, t): t for t in tasks}
                with tqdm(total=len(tasks), desc="构建进度", unit="文件", bar_format='{l_bar}{bar}| {n_fmt}/{total_fmt}[{elapsed}<{remaining}, {postfix}, {rate_fmt}]') as pbar:
                    for future in as_completed(futures):
                        task = futures[future]
                        short_tax = str(task[3])[:10]
                        pbar.set_postfix_str(f"{task[0]} ({short_tax}..)中")
                        if future.result(): success_count += 1
                        pbar.update(1)
            
            logger.info(f"✅ 成功归档了 {success_count} 个按分类学命名的核心共识序列文件")

    def print_plotting_instructions(self, best_summary_file):
        plots_out_prefix = self.output_dir / 'plots' / 'virus_analysis'
        logger.info("\n" + "=" * 60)
        logger.info("📊[后续绘图及可视化提示]")
        logger.info("多样本分析的 R 语言出图脚本已包含并执行完毕。若需绘制全基因组的单碱基深度折线图：")
        logger.info("请在命令行中使用 plot_virus_depth.py 工具：")
        logger.info(f"\033[1;36mpython ~/bin/plot_virus_depth.py -d {self.output_dir}/stat -m {best_summary_file} -o {self.output_dir}/plots/\033[0m")
        logger.info("=" * 60)

    def generate_report(self, df_best, total_raw_records, total_filtered_records):
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
            f.write(f"Reads提取状态: {'已提取阳性样本病毒序列 (存于 reads 下)' if getattr(self.args, 'extract_reads', False) else '否'}\n")
            f.write(f"共识提取状态: {'已针对核心代表株生成 (存于 consensus 下)' if self.args.consensus else '否'}\n")
            f.write(f"变异检测状态: {'已使用 ' + self.args.variant_caller + ' 进行检测 (存于 variants 下)' if getattr(self.args, 'call_variants', False) else '否'}\n\n")
            
            f.write(f"总处理样本数: {len(self.samples)}\n")
            f.write(f"[第一级] Raw 包含全量未过滤记录总数: {total_raw_records} (保存在 all_viruses.raw.{self.args.format})\n")
            f.write(f"[第二级] Summary 阈值过滤后记录总数: {total_filtered_records} (保存在 all_viruses.summary.{self.args.format})\n")
            f.write(f"[第三级] Best 同源去重后核心记录数: {len(df_best)} (保存在 all_viruses.best.summary.{self.args.format})\n\n")
            
            if 'Virus' in df_best.columns:
                f.write("同源去重后的最终确诊核心病毒群落概览:\n")
                f.write("-" * 60 + "\n")
                virus_counts = df_best['Virus'].value_counts()
                for virus, count in virus_counts.items():
                    vd = df_best[df_best['Virus'] == virus]
                    desc = vd['description'].iloc[0] if 'description' in vd.columns else ''
                    mean_tpm, mean_fpkm, mean_rpm = vd['TPM'].mean(), vd['FPKM'].mean(), vd['RPM'].mean()
                    
                    mean_unique = vd['Unique_Reads'].mean() if 'Unique_Reads' in vd.columns else 0
                    mean_multi = vd['Multi_Reads'].mean() if 'Multi_Reads' in vd.columns else 0
                    
                    f.write(f"{virus} ({desc[:40]}...):\n")
                    f.write(f"  样本检出数: {count}\n")
                    f.write(f"  核心平均 Unique_Reads: {mean_unique:.1f} (特异性高置信度比对)\n")
                    f.write(f"  核心平均 Multi_Reads:  {mean_multi:.1f} (跨物种多重同源比对)\n")
                    f.write(f"  核心平均 TPM:  {mean_tpm:.2f} (该毒株在整个群落中的相对占比)\n")
                    f.write(f"  核心平均 FPKM: {mean_fpkm:.2f} (长度归一化后的绝对载量)\n")
                    f.write(f"  核心平均 RPM:  {mean_rpm:.2f} (无偏差的绝对测序文库载量)\n\n")
        
def main():
    parser = argparse.ArgumentParser(
        description='病毒分析终极版 (带断点续跑、三级过滤、最优毒株筛选、Reads提取、共识及突变检测)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
增强功能说明:
  --resume       : 开启断点续跑。将自动跳过已完成的比对、深度计算、序列提取、共识提取和变异检测。
  --extract_reads: 将阳性样本中比对到病毒的Reads单独提取出来，生成 sample_virus_1/2/single.fastq.gz。
  --consensus    : 基于去重最优代表序列生成共识。
  --call_variants: 基于去重最优代表序列进行变异检测，支持 freebayes, ivar, lofreq。
        """
    )
    
    input_group = parser.add_argument_group('输入输出参数')
    input_group.add_argument('--input_dir', type=str, help='输入文件夹路径（包含fastq或fasta文件）')
    input_group.add_argument('--reference', type=str, required=True, help='参考基因组fasta文件路径')
    input_group.add_argument('--rvdb_taxon', type=str, help='RVDB 分类注释文件路径 (如: RVDB_Taxon_Current.tab)')
    input_group.add_argument('--email', type=str, help='[可选] 提供有效邮箱，用于启用 NCBI 联网获取未知病毒注释')
    input_group.add_argument('--api_key', type=str, help='[可选] NCBI API Key，可以提高联网获取的访问频率')
    input_group.add_argument('--output_dir', type=str, default='./virus_analysis', help='输出目录路径 [默认: ./virus_analysis]')
    input_group.add_argument('--sample_list', type=str, help='样本列表文件（每行: 样本名 fq/fa1 fq/fa2）')
    input_group.add_argument('--single_end', action='store_true', help='强制以单端模式处理所有测序数据')
    
    align_group = parser.add_argument_group('高级流程与变异控制参数')
    align_group.add_argument('--tool', type=str, choices=['strobealign', 'minimap2', 'bwa', 'bwa-mem2', 'bowtie2', 'hisat2'], default='strobealign', help='比对工具选择[默认: strobealign]')
    align_group.add_argument('--align_threads', type=int, default=8, help='每个样本比对的线程数[默认: 8]')
    align_group.add_argument('--extract_reads', action='store_true', help='提取比对上病毒的 reads，自动分离并输出为单双端 fastq.gz')
    align_group.add_argument('--consensus', action='store_true', help='基于去重后的最佳代表株自动生成带分类学归档的共识序列')
    align_group.add_argument('--call_variants', action='store_true', help='基于去重后的最佳代表株开启变异检测')
    align_group.add_argument('--variant_caller', type=str, choices=['freebayes', 'ivar', 'lofreq'], default='freebayes', help='变异检测工具选择 [默认: freebayes]')
    align_group.add_argument('--ref_gff', type=str, help='[可选] 参考基因组GFF注释文件路径 (仅使用iVar时适用)')
    align_group.add_argument('--resume', action='store_true', help='开启断点续跑模式，自动跳过已存在的耗时结果及请求')
    
    filter_group = parser.add_argument_group('过滤参数')
    filter_group.add_argument('--coverage', type=float, default=90.0, help='Coverage(%%)阈值[默认: 90]')
    filter_group.add_argument('--meandepth', type=float, default=10.0, help='MeanDepth阈值[默认: 10]')
    filter_group.add_argument('--min_fpkm', type=float, default=0.0, help='FPKM最小值过滤[默认: 0]')
    filter_group.add_argument('--min_tpm', type=float, default=0.0, help='TPM最小值过滤[默认: 0]')
    
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
    if args.ref_gff and not os.path.exists(args.ref_gff): parser.error(f"指定的GFF文件不存在: {args.ref_gff}")
    if args.input_dir and not os.path.exists(args.input_dir): parser.error(f"输入目录不存在: {args.input_dir}")
    if args.verbose: logger.setLevel(logging.DEBUG)
    
    pipeline = VirusAnalysisPipeline(args)
    pipeline.run_pipeline()

if __name__ == '__main__':
    main()
