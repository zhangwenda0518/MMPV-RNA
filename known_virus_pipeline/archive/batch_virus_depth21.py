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

⭐【全新突破: 生物学 Unique 引擎】
打破常规 MAPQ 限制！引入 TaxID 感知追踪，只要 Read 命中同物种所有变异株，即算 Unique！
⭐【全新突破: 多节段与 TaxID 智能解析】
完全适配包含 Accession/taxid/Species/Segment 的信息表，精准实现物理聚合！

【V13 顶级工程体系 (全量保留)】
6. [三级结果输出] all_viruses.raw / summary / best.summary 三级分流与详尽全景报告。
7. [全自动R绘图] 内置 R_PLOT_SCRIPT，生成多维度 (TPM/FPKM/Depth) 病毒丰度图。
8. [多进程架构] Reads提取、突变检测、共识构建全部配备独立 tqdm 进度条与并发池。
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

tax_col <- if("Adjusted_Species" %in% colnames(data)) "Adjusted_Species" else if("Species" %in% colnames(data)) "Species" else "Virus"
acc_col <- if("taxid" %in% colnames(data)) "taxid" else "Assembly"

data$Display_Name <- ifelse(data[[tax_col]] != "Unannotated" & data[[tax_col]] != "-", paste0(data[[tax_col]], "\n(TaxID: ", data[[acc_col]], ")"), data[[acc_col]])
data$Display_Name <- sapply(data$Display_Name, function(x) paste(strwrap(x, width = 45), collapse = "\n"))

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
    labs(title="Comprehensive Multi-metric Virus Abundance", x="Taxonomy (TaxID)", y="Value") +
    theme_func(base_size=13) + theme(plot.title=element_text(hjust=0.5, face="bold"), axis.text.y=element_text(size=9, face="italic"), legend.position="none") + coord_flip()
  
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
        self.total_raw, self.total_summary = 0, 0
        self.setup_output_directory()
        self.check_env()
        self._load_reference_lengths()

    def check_env(self):
        required_tools = {'samtools': 'BAM处理', 'coverm': '序列清洗过滤', 'pandepth': '深度计算', self.args.tool: '序列比对'}
        if self.args.tool == 'bowtie2': required_tools['bowtie2-build'] = 'Bowtie2建库'
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
                else: curr_len += len(line.strip())
            if curr_id: self.ref_length_dict[curr_id] = curr_len

    def setup_output_directory(self):
        subdirs = ['bam', 'stat', 'summary', 'logs', 'index', 'plots']
        if getattr(self.args, 'extract_reads', False): subdirs.append('reads')
        if self.args.consensus: subdirs.append('consensus')
        if getattr(self.args, 'call_variants', False): subdirs.append('variants')
            
        with Timer("创建输出目录结构"):
            for subdir in subdirs: (self.output_dir / subdir).mkdir(parents=True, exist_ok=True)
            log_file = self.output_dir / 'logs' / f'pipeline_fusion_{datetime.now().strftime("%Y%m%d_%H%M%S")}.log'
            fh = logging.FileHandler(log_file, encoding='utf-8')
            fh.setFormatter(logging.Formatter('[%(asctime)s] %(levelname)s - %(message)s'))
            logger.addHandler(fh)

    def write_sample_log(self, sample_name, content, mode='a'):
        with open(self.output_dir / 'logs' / f"{sample_name}.log", mode, encoding='utf-8') as f: f.write(content + "\n")

    def build_index(self):
        tool, ref_path = self.args.tool, Path(self.args.reference).resolve()
        if tool in ['strobealign', 'minimap2']:
            self.index_path = str(ref_path)
            return
        prefix = self.output_dir / 'index' / f"{ref_path.stem}_{tool}"
        self.index_path = str(prefix)
        cmd = []
        if tool == 'bwa': cmd = ['bwa', 'index', '-p', str(prefix), str(ref_path)]
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
                    for suf in ['_1', '_R1', '.R1', '.1', '_unmapped', '.unmapped', '_trimmed']: sname = sname.replace(suf, '')
                    samples.append({'name': sname, 'r1': str(fastq), 'r2': None})
        elif input_dir:
            patterns = [('*_1.f*q*', '*_2.f*q*'), ('*_R1*.f*q*', '*_R2*.f*q*'), ('*.R1.*', '*.R2.*'), ('*.1.f*q*', '*.2.f*q*')]
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
            with open(metrics_cache, 'r') as f: return filt_sec_bam, int(f.read().strip())

        self.write_sample_log(sname, f"\n--- 1. 比对与 CoverM 严格洗流 ({self.args.tool}) ---")
        align_cmd = []
        if self.args.tool == 'bowtie2':
            align_cmd = ['bowtie2', '-p', str(self.args.align_threads), '-x', self.index_path]
            align_cmd.extend(['-1', sample['r1'], '-2', sample['r2']] if sample['r2'] else ['-U', sample['r1']])
        elif self.args.tool == 'bwa':
            align_cmd = ['bwa', 'mem', '-v', '1', '-t', str(self.args.align_threads), self.index_path, sample['r1']]
            if sample['r2']: align_cmd.append(sample['r2'])

        ap = subprocess.Popen(align_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        vp = subprocess.Popen(['samtools', 'view', '-b', '-o', str(raw_bam)], stdin=ap.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        ap.stdout.close(); ap.communicate(); vp.communicate()

        count_res = subprocess.run(['samtools', 'view', '-c', '-F', '4', str(raw_bam)], capture_output=True, text=True)
        total_mapped = int(count_res.stdout.strip() or 0)
        with open(metrics_cache, 'w') as f: f.write(str(total_mapped))

        pid_val = self.args.min_pid * 100 if self.args.min_pid <= 1.0 else self.args.min_pid
        prop_val = self.args.min_aln_prop * 100 if self.args.min_aln_prop <= 1.0 else self.args.min_aln_prop
        filter_cmd = [
            'coverm', 'filter', '-b', str(raw_bam), '-o', str(filt_sec_bam) + '.unsorted',
            '--min-read-aligned-length', str(self.args.min_aln_len), '--min-read-percent-identity', str(pid_val),
            '--min-read-aligned-percent', str(prop_val), '--include-secondary', '-t', str(self.args.align_threads)
        ]
        if self.args.strict_paired and sample['r2']: filter_cmd.append('--proper-pairs-only')
        
        subprocess.run(filter_cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        subprocess.run(['samtools', 'sort', '-@', str(min(4, self.args.threads)), '-o', str(filt_sec_bam), str(filt_sec_bam) + '.unsorted'], check=True)
        pysam.index(str(filt_sec_bam))
        raw_bam.unlink(); Path(str(filt_sec_bam) + '.unsorted').unlink()
        return filt_sec_bam, total_mapped

    def build_read_sharing_network(self, sname, bam_path):
        self.write_sample_log(sname, f"--- 2. 图论网络聚类 (消除多重比对冗余) ---")
        read_to_contigs = defaultdict(set)
        contig_primary_reads = defaultdict(int)
        hit_contigs = set()
        
        with pysam.AlignmentFile(bam_path, "rb") as bam:
            for read in bam.fetch(until_eof=True):
                if read.is_unmapped: continue
                contig = bam.get_reference_name(read.reference_id)
                read_to_contigs[read.query_name].add(contig)
                hit_contigs.add(contig)
                if not read.is_secondary and not read.is_supplementary: contig_primary_reads[contig] += 1
                    
        shared_counts = defaultdict(int)
        for qname, contigs in read_to_contigs.items():
            clist = list(contigs)
            if len(clist) > 1:
                for i in range(len(clist)):
                    for j in range(i + 1, len(clist)): shared_counts[tuple(sorted([clist[i], clist[j]]))] += 1
                        
        adjacency = defaultdict(set)
        for (c1, c2), shared_n in shared_counts.items():
            min_total = min(contig_primary_reads.get(c1, 0), contig_primary_reads.get(c2, 0))
            if min_total > 0 and (shared_n / min_total) >= self.args.cluster_thresh:
                adjacency[c1].add(c2); adjacency[c2].add(c1)
                
        visited, clusters = set(), []
        for contig in contig_primary_reads.keys():
            if contig not in visited:
                cluster = set()
                queue = [contig]
                while queue:
                    curr = queue.pop(0)
                    if curr not in visited:
                        visited.add(curr); cluster.add(curr); queue.extend(list(adjacency[curr] - visited))
                clusters.append(cluster)
                
        exemplars = set(sorted(list(cl), key=lambda x: -contig_primary_reads.get(x, 0))[0] for cl in clusters)
        return exemplars, read_to_contigs, list(hit_contigs)

    # ⭐ [全新突破: TaxID 生物学 Unique 判定引擎]
    def generate_primary_bam_and_metrics(self, sname, filt_sec_bam, exemplars, read_to_contigs, tax_map):
        final_bam_path = self.output_dir / 'bam' / f"{sname}.final.bam"
        json_cache = self.output_dir / 'stat' / f"{sname}.micro_metrics.json"
        
        if self.args.resume and final_bam_path.exists() and json_cache.exists():
            with open(json_cache, 'r') as f: return final_bam_path, json.load(f)

        self.write_sample_log(sname, f"--- 3. 提取物种级 Unique (TaxID Level) 与微观多态性 ---")
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
                    
                    # ⭐ 核心逻辑：只要 Read 追踪命中的所有 Contigs 都同属一个 TaxID，即判定为绝对 Unique！
                    mapped_contigs = read_to_contigs[r.query_name]
                    mapped_taxids = set()
                    for c in mapped_contigs:
                        taxid = tax_map.get(c, {}).get('taxid', 'Unannotated')
                        if taxid in ['Unannotated', '-']: mapped_taxids.add(c) # 无注释时以ID自身孤立防合并
                        else: mapped_taxids.add(str(taxid))
                        
                    if len(mapped_taxids) == 1: unique_n += 1
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
                pi_list = [ (sum(c.values()) - max(c.values())) / sum(c.values()) for c in base_counts if sum(c.values()) > 0 ]
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
        
        depth_dict = {}
        with gzip.open(str(stat_prefix) + '.chr.stat.gz', 'rt') as f:
            for line in f:
                if not line.startswith('#'):
                    p = line.strip().split()
                    if len(p) >= 6: depth_dict[p[0]] = {'CoveredBases': int(p[2]), 'TotalDepth': int(p[3]), 'Coverage(%)': float(p[4]), 'MeanDepth': float(p[5])}
        
        for m in metrics:
            d = depth_dict.get(m['Accession'], {})
            m.update({k: d.get(k, 0) for k in ['CoveredBases', 'TotalDepth', 'Coverage(%)', 'MeanDepth']})
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
            exemplars, read_to_contigs, hit_contigs = self.build_read_sharing_network(sname, str(filt_sec_bam))
            if not exemplars: return []
            
            tax_map = {}
            if getattr(self.args, 'ref_info', None): tax_map.update(self.fast_load_ref_info(hit_contigs))
            unannotated = [v for v in hit_contigs if v not in tax_map]
            if unannotated and getattr(self.args, 'email', None): tax_map.update(self.fetch_ncbi_taxonomy(unannotated))
            
            final_bam, micro_metrics = self.generate_primary_bam_and_metrics(sname, filt_sec_bam, exemplars, read_to_contigs, tax_map)
            return self.run_pandepth(final_bam, sname, micro_metrics, total_mapped)
        except Exception as e:
            logger.error(f"❌ 样本 {sname} 致命错误: {e}")
            return []

    # ⭐ [全新突破: Ref_Info 智能表单解析器]
    def fast_load_ref_info(self, found_viruses):
        if not getattr(self.args, 'ref_info', None) or not found_viruses: return {}
        target_exact, target_base = set(found_viruses), {v.split('.')[0] for v in found_viruses}
        info_dict = {}
        try:
            with open(self.args.ref_info, 'r', encoding='utf-8') as f:
                header = None
                idx_acc, idx_tax, idx_sp, idx_seg, idx_desc = 0, 3, 4, 13, 8  # 默认对应 1,4,5,14,9 列
                for line in f:
                    if line.startswith('####'): continue
                    parts = line.rstrip('\n').split('\t')
                    if header is None and 'Accession' in line:
                        header = parts
                        if 'Accession' in header: idx_acc = header.index('Accession')
                        if 'taxid' in header: idx_tax = header.index('taxid')
                        if 'Species' in header: idx_sp = header.index('Species')
                        if 'Segment' in header: idx_seg = header.index('Segment')
                        if 'GenBank_Title' in header: idx_desc = header.index('GenBank_Title')
                        continue
                    
                    if len(parts) <= max(idx_acc, idx_tax): continue
                    acc = parts[idx_acc].strip()
                    acc_base = acc.split('.')[0]
                    
                    if acc in target_exact or acc_base in target_base:
                        info_dict[acc] = info_dict[acc_base] = {
                            'taxid': parts[idx_tax].strip() if len(parts) > idx_tax else "Unannotated",
                            'species': parts[idx_sp].strip() if len(parts) > idx_sp else "Unannotated",
                            'segment': parts[idx_seg].strip() if len(parts) > idx_seg else "",
                            'description': parts[idx_desc].strip() if len(parts) > idx_desc else ""
                        }
            return info_dict
        except Exception as e:
            logger.warning(f"读取参考信息表失败: {e}")
            return {}

    def fetch_ncbi_taxonomy(self, accessions):
        try: from Bio import Entrez, SeqIO
        except ImportError: return {}
        ncbi_dict = {}
        cache_file = self.output_dir / 'summary' / 'ncbi_taxonomy_cache.tsv'
        if cache_file.exists() and cache_file.stat().st_size > 0:
            try:
                df_cache = pd.read_csv(cache_file, sep='\t', dtype=str).fillna("")
                for _, row in df_cache.iterrows(): 
                    ncbi_dict[row['accession']] = {'taxid': str(row['taxid']), 'species': str(row['species']), 'segment': '', 'description': str(row['description'])}
            except Exception: cache_file.unlink()
                
        to_fetch = [acc for acc in accessions if acc not in ncbi_dict]
        if not to_fetch: return ncbi_dict

        Entrez.email = getattr(self.args, 'email', None)
        Entrez.api_key = getattr(self.args, 'api_key', None)
        with open(cache_file, 'a', encoding='utf-8') as cf:
            if not cache_file.exists() or cache_file.stat().st_size == 0: cf.write("accession\ttaxid\tspecies\tdescription\n")
            for acc in to_fetch:
                try:
                    handle = Entrez.efetch(db="nuccore", id=acc, rettype="gb", retmode="text")
                    record = SeqIO.read(handle, "genbank")
                    handle.close()
                    desc = record.description
                    species, taxid = "Unannotated", "Unannotated"
                    for feat in record.features:
                        if feat.type == 'source':
                            species = feat.qualifiers.get('organism', ['Unannotated'])[0]
                            for xref in feat.qualifiers.get('db_xref', []):
                                if xref.startswith('taxon:'): taxid = xref.split(':')[1]
                    if taxid == "Unannotated": taxid = species
                    ncbi_dict[acc] = {'taxid': taxid, 'species': species, 'segment': '', 'description': desc}
                    cf.write(f"{acc}\t{taxid}\t{species}\t{desc}\n")
                    cf.flush()
                except: pass
                time.sleep(0.3)
        return ncbi_dict

    def summarize_results(self, all_metrics):
        if not all_metrics: return None
        df = pl.DataFrame(all_metrics)
        found_viruses = df['Accession'].unique().to_list()
        tax_map = {}
        if getattr(self.args, 'ref_info', None): tax_map.update(self.fast_load_ref_info(found_viruses))
        unannotated = [v for v in found_viruses if v not in tax_map]
        if unannotated and getattr(self.args, 'email', None): tax_map.update(self.fetch_ncbi_taxonomy(unannotated))
        
        def get_tax(acc, key): return tax_map.get(acc, {}).get(key, "" if key == 'segment' else "Unannotated")
        
        df = df.with_columns([
            pl.col("Accession").map_elements(lambda x: get_tax(x, 'taxid'), return_dtype=pl.Utf8).alias("taxid"),
            pl.col("Accession").map_elements(lambda x: get_tax(x, 'species'), return_dtype=pl.Utf8).alias("Species"),
            pl.col("Accession").map_elements(lambda x: get_tax(x, 'segment'), return_dtype=pl.Utf8).alias("Segment"),
            pl.col("Accession").map_elements(lambda x: get_tax(x, 'description'), return_dtype=pl.Utf8).alias("description")
        ])
        
        # 多节段格式化显示: DNA-A:HM007119.1
        df = df.with_columns(
            pl.when(pl.col("Segment") != "").then(pl.col("Segment") + ":" + pl.col("Accession"))
            .otherwise(pl.col("Accession")).alias("Seg_Acc_Str")
        )
        
        ext, sep = self.args.format, "," if self.args.format == "csv" else "\t"
        df.write_csv(str(self.output_dir / "summary" / f"all_viruses.raw.{ext}"), separator=sep)
        self.total_raw = len(df)
        
        # 基于 TaxID 物理聚合所有同物种的不同节段
        asm_df = df.group_by(["Sample", "taxid"]).agg([
            pl.col("MappedReads").sum(), pl.col("Unique_Reads").sum(), pl.col("Multi_Reads").sum(),
            pl.col("CoveredBases").sum(), pl.col("Length").sum().alias("Asm_Length"),
            pl.col("Coverage(%)").mean().alias("Coverage(%)"), pl.col("MeanDepth").mean().alias("MeanDepth"),
            pl.col("Read_ANI").mean().alias("Avg_Read_ANI"), pl.col("Diversity_Pi").mean().alias("Avg_Pi"),
            pl.col("TPM").sum().alias("Asm_TPM"), pl.col("RPM").sum().alias("Asm_RPM"),
            pl.col("FPKM").sum().alias("Asm_FPKM"), pl.col("Relative_Abundance(%)").sum().alias("Asm_Rel_Abund(%)"),
            pl.col("Species").first(), pl.col("description").first(), pl.col("Seg_Acc_Str").alias("Segment_Accessions")
        ])
        
        asm_df = asm_df.with_columns((pl.col("Unique_Reads") / pl.col("MappedReads") * 100).round(2).alias("Unique(%)"))
        
        sp_thresh = self.args.sp_thresh
        asm_df = asm_df.with_columns(
            pl.when(pl.col("Avg_Read_ANI") < sp_thresh)
            .then(pl.concat_str([pl.lit("s__unclassified_"), pl.col("Species").str.replace_all(" ", "_")]))
            .otherwise(pl.col("Species")).alias("Adjusted_Species")
        ).with_columns(pl.col("Segment_Accessions").list.join(","))
        
        final_df = asm_df.filter(
            (pl.col("Coverage(%)") >= self.args.coverage) & (pl.col("MeanDepth") >= self.args.meandepth) &
            (pl.col("Asm_FPKM") >= self.args.min_fpkm) & (pl.col("Asm_TPM") >= self.args.min_tpm)
        )
        final_df.write_csv(str(self.output_dir / "summary" / f"all_viruses.summary.{ext}"), separator=sep)
        self.total_summary = len(final_df)
        
        best_df = final_df.sort(["Sample", "Adjusted_Species", "Coverage(%)", "MeanDepth"], descending=[False, False, True, True]).unique(subset=["Sample", "Adjusted_Species"], keep="first", maintain_order=True)
        best_csv = self.output_dir / "summary" / f"all_viruses.best.summary.{ext}"
        best_df.write_csv(str(best_csv), separator=sep)
        self.df_best = best_df
        
        logger.info(f"✅ 三级结果保存完成! (Raw: {self.total_raw} -> Summary: {self.total_summary} -> Best: {len(best_df)})")
        self.generate_report_txt()
        self.generate_plots(best_csv)
        return best_df.to_dicts()

    def generate_plots(self, best_summary_file):
        plots_dir = self.output_dir / 'plots'
        with open(plots_dir / 'virus_frequency_plot.R', 'w', encoding='utf-8') as f: f.write(R_PLOT_SCRIPT)
        if shutil.which('Rscript') is None: return
        subprocess.run(['Rscript', str(plots_dir / 'virus_frequency_plot.R'), '-i', str(best_summary_file), '-o', str(plots_dir / 'virus_analysis'), '--multi-plot', '--log10-transform'], capture_output=True)

    def generate_report_txt(self):
        with open(self.output_dir / 'summary' / 'analysis_report.txt', 'w') as f:
            f.write("=" * 80 + "\n工业级宏病毒丰度与多节段全景分析报告 (Fusion Edition)\n" + "=" * 80 + "\n\n")
            f.write(f"【核心数据流漏斗统计】\n  1. 处理样本: {len(self.samples)}\n  2. 留存 (Raw): {self.total_raw}\n  3. 达标 (Summary): {self.total_summary}\n  4. 最优代表 (Best): {len(self.df_best)}\n\n")
            if self.df_best is not None and len(self.df_best) > 0:
                f.write("【确诊群落全景】\n" + "-" * 60 + "\n")
                asm_counts = self.df_best['Adjusted_Species'].value_counts()
                c_name = 'count' if 'count' in asm_counts.columns else 'counts'
                for row in asm_counts.iter_rows(named=True):
                    asm, count = row['Adjusted_Species'], row[c_name]
                    vd = self.df_best.filter(pl.col('Adjusted_Species') == asm)
                    f.write(f"🎯 {asm}:\n   ├─ 检出数: {count}/{len(self.samples)}\n   ├─ 生物特异 Unique: {vd['Unique(%)'].mean():.1f}%\n   ├─ 多态性(Pi): {vd['Avg_Pi'].mean():.4f}\n   ├─ 平均 TPM: {vd['Asm_TPM'].mean():.2f}\n   └─ 涵盖节段: {vd['Segment_Accessions'][0]}\n\n")

    def _extract_pure_bam_pysam(self, global_bam_path, acc, out_bam_path):
        try:
            with pysam.AlignmentFile(global_bam_path, "rb") as in_bam:
                target_tid = in_bam.get_tid(acc)
                if target_tid < 0: return False
                new_header = {'HD': {'VN': '1.0', 'SO': 'coordinate'}, 'SQ': [{'LN': in_bam.lengths[target_tid], 'SN': acc}]}
                with pysam.AlignmentFile(out_bam_path, "wb", header=new_header) as out_bam:
                    for read in in_bam.fetch(acc):
                        read.reference_id, read.next_reference_id = 0, (0 if read.next_reference_id == target_tid else -1)
                        if read.next_reference_id == -1: read.mate_is_unmapped = True
                        out_bam.write(read)
            return True
        except Exception: return False

    def extract_fastas_with_python(self, target_ids, virus_to_tax, out_dir):
        target_set = set(target_ids)
        current_id, current_seq = None, []
        def save_record(vid, seq_list):
            if vid in target_set:
                vname = re.sub(r'[\\/*?:"<>| ]', "_", vid)
                safe_tax = re.sub(r'[\\/*?:"<>| ]', "_", virus_to_tax.get(vid, 'Unannotated'))[:50]
                v_dir = out_dir / f"{safe_tax}_{vname}"
                v_dir.mkdir(parents=True, exist_ok=True)
                with open(v_dir / f"{safe_tax}_{vname}.ref.fasta", 'w') as outf: outf.write(f">{vid}\n" + "".join(seq_list) + "\n")
        with open(self.args.reference, 'r') as f:
            for line in f:
                if line.startswith('>'):
                    if current_id: save_record(current_id, current_seq)
                    current_id, current_seq = line[1:].strip().split()[0], []
                else: current_seq.append(line)
            if current_id: save_record(current_id, current_seq)

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

    def _get_pure_acc(self, segment_str):
        return [item.split(':')[-1] for item in segment_str.split(',')]

    # ⭐ 将闭包函数提出来变成可 Pickle 的类方法
    def _run_single_variant_calling(self, task):
        sample, virus, taxonomy = task
        safe_tax = re.sub(r'[\\/*?:\'<>| ]', '_', taxonomy)[:50]
        safe_virus = re.sub(r'[\\/*?:\'<>| ]', '_', virus)
        folder = f"{safe_tax}_{safe_virus}"
        
        v_dir = self.output_dir / 'variants' / folder
        ref, out = v_dir / f"{folder}.ref.fasta", v_dir / f"{sample}.{folder}.variants.vcf"
        bam, fixed = self.output_dir / 'bam' / f"{sample}.final.bam", v_dir / f"{sample}.{folder}.fixed.bam"
        
        if self.args.resume and out.exists(): return True
        try:
            if not Path(str(ref)+'.fai').exists(): subprocess.run(f"samtools faidx '{ref}'", shell=True)
            if not self._extract_pure_bam_pysam(str(bam), virus, str(fixed)): return False
            subprocess.run(f"samtools index '{fixed}'", shell=True)
            if self.args.variant_caller == 'freebayes': subprocess.run(f"freebayes -p 1 -f '{ref}' '{fixed}' > '{out}'", shell=True, stderr=subprocess.DEVNULL)
            else: subprocess.run(f"samtools mpileup -A -aa -d 0 -Q 0 --reference '{ref}' '{fixed}' | ivar variants -r '{ref}' -p '{str(out)[:-4]}'", shell=True, stderr=subprocess.DEVNULL)
            
            # 清理临时文件，安全删除
            if fixed.exists(): fixed.unlink()
            bai_file = Path(str(fixed)+'.bai')
            if bai_file.exists(): bai_file.unlink()
            return True
        except: return False

    def run_variants_calling(self, best_results):
        logger.info(f"\n🔬 开始执行变异检测 ({self.args.variant_caller})...")
        tasks, virus_to_tax = [], {}
        for r in best_results:
            for acc in self._get_pure_acc(r['Segment_Accessions']):
                virus_to_tax[acc] = r.get('Adjusted_Species', 'Unannotated')
                tasks.append((r['Sample'], acc, r.get('Adjusted_Species', 'Unannotated')))
        self.extract_fastas_with_python(list(virus_to_tax.keys()), virus_to_tax, self.output_dir / 'variants')
        
        # 使用提出来的类方法进行多进程映射
        with ProcessPoolExecutor(max_workers=self.args.threads) as executor:
            list(tqdm(executor.map(self._run_single_variant_calling, tasks), total=len(tasks), desc="变异检测"))

    # ⭐ 将闭包函数提出来变成可 Pickle 的类方法
    def _run_single_consensus(self, task):
        sample, virus, depth, taxonomy = task
        safe_tax = re.sub(r'[\\/*?:\'<>| ]', '_', taxonomy)[:50]
        safe_virus = re.sub(r'[\\/*?:\'<>| ]', '_', virus)
        folder = f"{safe_tax}_{safe_virus}"
        
        v_dir = self.output_dir / 'consensus' / folder
        fixed = v_dir / f"{sample}.{folder}.fixed.bam"
        ref, out = v_dir / f"{folder}.ref.fasta", v_dir / f"{sample}.{folder}.consensus.fasta"
        
        if self.args.resume and out.exists(): return True
        try:
            if not self._extract_pure_bam_pysam(str(self.output_dir / 'bam' / f"{sample}.final.bam"), virus, str(fixed)): return False
            subprocess.run(f"samtools index '{fixed}'", shell=True)
            subprocess.run(f"viral_consensus -i '{fixed}' -r '{ref}' -o '{out}' --min_depth {max(1, min(10, int(math.floor(depth))))}", shell=True, stderr=subprocess.DEVNULL)
            
            # 清理临时文件，安全删除
            if fixed.exists(): fixed.unlink()
            bai_file = Path(str(fixed)+'.bai')
            if bai_file.exists(): bai_file.unlink()
            return True
        except: return False

    def build_consensus_sequences(self, best_results):
        logger.info("\n🧬 开始生成核心共识序列...")
        tasks, virus_to_tax = [], {}
        for r in best_results:
            for acc in self._get_pure_acc(r['Segment_Accessions']):
                virus_to_tax[acc] = r.get('Adjusted_Species', 'Unannotated')
                tasks.append((r['Sample'], acc, r['MeanDepth'], r.get('Adjusted_Species', 'Unannotated')))
        self.extract_fastas_with_python(list(virus_to_tax.keys()), virus_to_tax, self.output_dir / 'consensus')
        
        # 使用提出来的类方法进行多进程映射
        with ProcessPoolExecutor(max_workers=self.args.threads) as executor:
            list(tqdm(executor.map(self._run_single_consensus, tasks), total=len(tasks), desc="共识生成"))

    def run(self):
        logger.info("=" * 60); logger.info("💥 Ultimate Virus Pipeline (Taxonomy & Multi-Seg Edition)"); logger.info("=" * 60)
        with Timer("整体管线运行"):
            self.build_index(); self.find_samples()
            max_p = self.args.parallel_jobs if getattr(self.args, 'parallel', False) and self.args.parallel_jobs else max(1, self.args.threads // self.args.align_threads)
            all_metrics = []
            with ProcessPoolExecutor(max_workers=max_p) as executor:
                for res in tqdm(executor.map(self.process_single_sample, self.samples), total=len(self.samples), desc="主比对与特征提取"):
                    if res: all_metrics.extend(res)
                        
            best = self.summarize_results(all_metrics)
            if best:
                if getattr(self.args, 'extract_reads', False): self.extract_mapped_reads(best)
                if getattr(self.args, 'call_variants', False): self.run_variants_calling(best)
                if getattr(self.args, 'consensus', False): self.build_consensus_sequences(best)
            else: logger.warning("⚠️ 未找到任何满足阈值的病毒记录")

def main():
    parser = argparse.ArgumentParser(description='工业级宏病毒全栈管线 (Zero-Compromise Edition)')
    io_g = parser.add_argument_group('输入输出参数')
    io_g.add_argument('-i', '--input_dir', help='输入文件夹路径')
    io_g.add_argument('-r', '--reference', required=True, help='参考基因组fasta文件')
    io_g.add_argument('-o', '--output_dir', default='./virus_out', help='输出目录')
    io_g.add_argument('-l', '--sample_list', help='样本列表文件')
    io_g.add_argument('--single_end', action='store_true', help='使用单端数据')
    
    # ⭐ 修复参数名：同时兼容 --ref_info 和 --rvdb_taxon，防止你以前的脚本命令跑空
    io_g.add_argument('--ref_info', '--rvdb_taxon', dest='ref_info', type=str, help='参考信息 TSV 文件 (支持 Accession, taxid, Species, Segment 等列)')

    align_g = parser.add_argument_group('高级流程参数')
    align_g.add_argument('--tool', choices=['bwa', 'bowtie2', 'hisat2'], default='bowtie2', help='比对工具')
    align_g.add_argument('--align_threads', type=int, default=8, help='每样本比对线程')
    align_g.add_argument('--extract_reads', action='store_true', help='提取阳性reads')
    align_g.add_argument('--consensus', action='store_true', help='生成共识序列')
    align_g.add_argument('--call_variants', action='store_true', help='开启变异检测')
    align_g.add_argument('--variant_caller', choices=['freebayes', 'ivar'], default='freebayes')
    align_g.add_argument('--resume', action='store_true', help='断点续跑')

    flt_g = parser.add_argument_group('高级过滤与图论聚类 (V20 特性)')
    flt_g.add_argument('--strict_paired', action='store_true')
    flt_g.add_argument('--min_aln_len', type=int, default=80)
    flt_g.add_argument('--min_aln_prop', type=float, default=0.85)
    flt_g.add_argument('--min_pid', type=float, default=0.90)
    flt_g.add_argument('--cluster_thresh', type=float, default=0.30)
    flt_g.add_argument('--sp_thresh', type=float, default=0.90)
    
    thresh_g = parser.add_argument_group('丰度过滤参数 (V13 特性)')
    thresh_g.add_argument('--coverage', type=float, default=90.0)
    thresh_g.add_argument('--meandepth', type=float, default=10.0)
    thresh_g.add_argument('--min_fpkm', type=float, default=0.0)
    thresh_g.add_argument('--min_tpm', type=float, default=0.0)

    para_g = parser.add_argument_group('并行与杂项')
    para_g.add_argument('--parallel', action='store_true')
    para_g.add_argument('--parallel_jobs', type=int)
    para_g.add_argument('-t', '--threads', type=int, default=4)
    para_g.add_argument('--format', choices=['csv', 'tsv'], default='tsv')
    para_g.add_argument('--email', type=str)
    para_g.add_argument('--api_key', type=str)
    para_g.add_argument('--verbose', action='store_true')
    
    args = parser.parse_args()
    if not args.input_dir and not args.sample_list: parser.error("必须提供 -i 或 -l")
    if args.verbose: logger.setLevel(logging.DEBUG)
    VirusMetagenomicsPipeline(args).run()

if __name__ == '__main__':
    main()
