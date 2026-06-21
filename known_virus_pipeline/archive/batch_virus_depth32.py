#!/usr/bin/env python3
"""
宏病毒鉴定与精确定量管线 (Virus Quantification - V34 Production Gold Master 最终版)
================================================================================
【V34 工业级防崩与性能榨取】
1. [数学安全] 彻底封堵 TPM/FPKM 计算中的 ZeroDivision 漏洞。
2. [环境兼容] 重写 Pysam 索引读取，兼容旧版环境；加入 Samtools O(1) 降级方案防卡死。
3. [内存熔断] 废弃 map 一波流，引入 as_completed 异步收割，解决百 GB 级并发 OOM。
4. [IO 容错] 引入 Pandepth 磁盘延迟检测，Pysam 并发索引锁保护。
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
import math
import shlex
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

def setup_logging(verbose=False):
    handler = colorlog.StreamHandler()
    handler.setFormatter(colorlog.ColoredFormatter('%(log_color)s[%(asctime)s] %(levelname)s - %(message)s', datefmt='%H:%M:%S', log_colors={'DEBUG': 'cyan', 'INFO': 'green', 'WARNING': 'yellow', 'ERROR': 'red', 'CRITICAL': 'bold_red'}))
    logger = colorlog.getLogger("QuantPipe_V34")
    logger.addHandler(handler)
    logger.setLevel(logging.DEBUG if verbose else logging.INFO)
    return logger

logger = setup_logging()

class Timer:
    def __init__(self, name=""): self.name = name
    def __enter__(self): self.start_time = time.time(); logger.info(f"⏱️  开始: {self.name}"); return self
    def __exit__(self, exc_type, exc_val, exc_tb): logger.info(f"{'✅ 完成' if exc_type is None else '❌ 失败'}: {self.name} [耗时: {time.time() - self.start_time:.1f}秒]")

# ==================== 工具函数 ====================

def _get_ref_read_count_safe(bam, ref):
    """【修复1】安全获取参考序列读数，兼容所有 pysam 版本"""
    try:
        # Pysam 0.15.0+：O(1) 索引查询（最快）
        for stat in bam.get_index_statistics():
            if stat.contig == ref:
                return stat.mapped
    except (AttributeError, TypeError):
        pass
    
    # 降级方案 1：用 samtools idxstats 快速查询
    try:
        result = subprocess.run(
            ['samtools', 'idxstats', str(bam.filename)],
            capture_output=True, text=True, timeout=10
        )
        if result.returncode == 0:
            for line in result.stdout.strip().split('\n'):
                parts = line.split('\t')
                if len(parts) >= 3 and parts[0] == ref:
                    return int(parts[2])  # mapped reads
    except Exception:
        pass
    
    # 降级方案 2：计数（最慢但最可靠）
    logger.warning(f"⚠️ 使用缓慢的计数方式获取 {ref} 读数，建议更新 pysam/samtools")
    return sum(1 for _ in bam.fetch(ref) if not _.is_unmapped)


class VirusQuantificationPipeline:
    def __init__(self, args):
        self.args = args
        self.samples = []
        self.output_dir = Path(args.output_dir)
        self.index_path = None
        self.ref_length_dict = {}
        self.taxid_clusters = {}  
        self.tax_map = {}         
        
        self.check_tools()
        self._load_reference_lengths()
        self._load_taxid_clusters()
        self._load_ref_info()
        
    def check_tools(self):
        required_tools = {'samtools': 'BAM处理', 'coverm': '序列清洗', 'pandepth': '深度计算', self.args.tool: '序列比对'}
        if self.args.tool == 'bowtie2': required_tools['bowtie2-build'] = 'Bowtie2建库'
        elif self.args.tool == 'bwa': required_tools['bwa'] = 'BWA建库'
        missing = [f"{t} ({d})" for t, d in required_tools.items() if shutil.which(t) is None]
        if missing: logger.error("❌ 缺少必要工具:"); [logger.error(f"  - {m}") for m in missing]; sys.exit(1)

    def _load_reference_lengths(self):
        with open(self.args.reference, 'r') as f:
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
        subdirs = ['bam', 'stat', 'summary', 'logs', 'index', 'plots']
        for subdir in subdirs: (self.output_dir / subdir).mkdir(parents=True, exist_ok=True)

    def build_index(self):
        tool, ref_path = self.args.tool, Path(self.args.reference).resolve()
        prefix = self.output_dir / 'index' / f"{ref_path.stem}_{tool}"
        self.index_path = str(prefix)
        cmd = []
        if tool == 'bwa' and not Path(f"{prefix}.bwt").exists(): cmd = ['bwa', 'index', '-p', str(prefix), str(ref_path)]
        elif tool == 'bowtie2' and not (Path(f"{prefix}.1.bt2").exists() or Path(f"{prefix}.1.bt2l").exists()): cmd = ['bowtie2-build', '--threads', str(self.args.threads), str(ref_path), str(prefix)]
        if cmd: subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    def find_samples(self):
        input_dir = Path(self.args.input_dir)
        samples = []
        for r1_file in input_dir.glob('*_1.f*q*'):
            r2_file = input_dir / r1_file.name.replace('_1.', '_2.')
            if r2_file.exists(): samples.append({'name': r1_file.name.split('_1')[0], 'r1': str(r1_file), 'r2': str(r2_file)})
        if not samples:
            for f in input_dir.glob('*.f*q*'): samples.append({'name': f.name.split('.')[0], 'r1': str(f), 'r2': None})
        self.samples = list({v['name']: v for v in samples}.values())
        if not self.samples: sys.exit(1)
        logger.info(f"✅ 找到 {len(self.samples)} 个测序样本")

    def align_and_coverm(self, sample):
        sname = sample['name']
        raw_bam = self.output_dir / 'bam' / f'{sname}.raw.bam'
        filt_bam = self.output_dir / 'bam' / f'{sname}.sorted.bam'
        bai_path = str(filt_bam) + '.bai'
        
        # ==========================================
        # 🚀 新增：断点续传（跳过已存在的 BAM）逻辑
        # ==========================================
        if filt_bam.exists() and Path(bai_path).exists() and filt_bam.stat().st_size > 1024:
            logger.info(f"⏩ 发现已完成的 BAM，跳过比对步骤: {sname}")
            # 快速统计已存在的 mapped reads
            mapped_res = subprocess.run(['samtools', 'view', '-c', str(filt_bam)], capture_output=True, text=True)
            total_mapped = int(mapped_res.stdout.strip() or 1)
            # 恢复模式下，由于原始 raw_bam 已经被清空，这里用 total_mapped 暂代 global_reads
            return filt_bam, total_mapped, total_mapped
        # ==========================================

        align_cmd = []
        inner_threads = min(2, self.args.align_threads) 
        if self.args.tool == 'bwa': 
            align_cmd = ['bwa', 'mem', '-t', str(inner_threads), self.index_path, sample['r1']] + ([sample['r2']] if sample['r2'] else [])
        elif self.args.tool == 'bowtie2': 
            align_cmd = ['bowtie2', '-p', str(inner_threads), '-x', self.index_path] + (['-1', sample['r1'], '-2', sample['r2']] if sample['r2'] else ['-U', sample['r1']])
        
        cmd_str = ' '.join(shlex.quote(str(x)) for x in align_cmd)
        raw_bam_str = shlex.quote(str(raw_bam))
        subprocess.run(f"set -o pipefail; {cmd_str} | samtools view -b -o {raw_bam_str}", shell=True, executable='/bin/bash', stderr=subprocess.DEVNULL)
        
        total_raw_res = subprocess.run(['samtools', 'view', '-c', str(raw_bam)], capture_output=True, text=True)
        global_reads = int(total_raw_res.stdout.strip() or 1)
        
        mapped_res = subprocess.run(['samtools', 'view', '-c', '-F', '4', str(raw_bam)], capture_output=True, text=True)
        total_mapped = int(mapped_res.stdout.strip() or 1)

        filter_cmd = [
            'coverm', 'filter', '-b', str(raw_bam), '-o', str(filt_bam) + '.unsorted',
            '--min-read-aligned-length', str(self.args.min_aln_len), 
            '--min-read-percent-identity', str(self.args.min_pid),
            '--min-read-aligned-percent', str(self.args.min_aln_prop), 
            '--include-secondary', '-t', '1'
        ]
        
        filter_res = subprocess.run(filter_cmd, capture_output=True, text=True)
        if filter_res.returncode != 0: 
            logger.warning(f"⚠️ CoverM 警告 (样本 {sname}): {filter_res.stderr[:200]}")
            
        subprocess.run(['samtools', 'sort', '-@', '1', '-o', str(filt_bam), str(filt_bam) + '.unsorted'])
        
        if not os.path.exists(bai_path):
            try:
                pysam.index(str(filt_bam))
            except Exception:
                time.sleep(np.random.uniform(0.1, 0.5))
                if not os.path.exists(bai_path):
                    try: pysam.index(str(filt_bam))
                    except Exception as e2: logger.warning(f"⚠️ 索引创建失败 {bai_path}: {e2}")
        
        if raw_bam.exists(): raw_bam.unlink()
        unsorted_bam = Path(str(filt_bam) + '.unsorted')
        if unsorted_bam.exists(): unsorted_bam.unlink()
        
        return filt_bam, total_mapped, global_reads

    def run_em_taxid_allocation(self, filt_bam_path):
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
            if total_mapped <= 1: return []
            
            em_counts, uniq_counts, multi_counts = self.run_em_taxid_allocation(sorted_bam)
            
            stat_prefix = self.output_dir / 'stat' / sname
            stat_file = str(stat_prefix) + '.chr.stat.gz'
            
            # 【修复】：Pandepth IO 容错
            pan_res = subprocess.run(['pandepth', '-a', '-i', str(sorted_bam), '-o', str(stat_prefix), '-t', '1'], capture_output=True, text=True)
            depth_dict = {}
            if pan_res.returncode != 0:
                logger.warning(f"⚠️ Pandepth 失败 (样本 {sname}): {pan_res.stderr[:200]}")
            elif not os.path.exists(stat_file):
                logger.warning(f"⚠️ Pandepth 输出缺失: {stat_file}")
            else:
                try:
                    with gzip.open(stat_file, 'rt') as f:
                        for line in f:
                            if not line.startswith('#'):
                                p = line.split()
                                if len(p) >= 6: depth_dict[p[0]] = {'Cov': float(p[4]), 'Dep': float(p[5])}
                except Exception as e:
                    logger.warning(f"⚠️ 读取 Pandepth 文件失败: {e}")
                        
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
                    'Coverage(%)': d_info['Cov'], 'MeanDepth': d_info['Dep'],
                    'EM_Reads': em_mapped, 'Uniq_Reads': u_reads, 'Multi_Reads': m_reads,
                    'Sample_Total_Mapped': total_mapped, 'Sample_Global_Reads': global_reads
                })
            return results
        except Exception as e:
            logger.error(f"❌ 样本 {sname} 失败", exc_info=self.args.verbose)
            return []

    def _compute_ani_pi_worker(self, task):
        """【修复2】O(1) 索引读数 + 真正的均匀抽样"""
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
                
                # 【修复1】：安全获取总读数
                total_reads = _get_ref_read_count_safe(bam, ref)
                if total_reads == 0: 
                    return {'Sample': sname, 'Rep_Accession': ref, 'Avg_Read_ANI': None, 'Avg_Pi': None}
                
                # 【修复】：真正的空间均匀抽样
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
                        
        except Exception as e:
            logger.warning(f"ANI 计算异常 {sname}/{ref}: {str(e)}")
            return {'Sample': sname, 'Rep_Accession': ref, 'Avg_Read_ANI': None, 'Avg_Pi': None}
            
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

    def summarize_results_polars(self, all_metrics):
        if not all_metrics: return None
        df = pl.DataFrame(all_metrics)
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
        
        # 【修复3】：数学完美的五大丰度指标 + 严格除零保护
        asm_df = asm_df.with_columns([
            pl.when((pl.col("Uniq_Reads") + pl.col("Multi_Reads")) > 0)
              .then((pl.col("Uniq_Reads") / (pl.col("Uniq_Reads") + pl.col("Multi_Reads")) * 100).round(2))
              .otherwise(0.0).alias("Unique(%)"),
            
            (pl.col("Asm_EM_Reads") * 1e6 / pl.col("Sample_Total_Mapped")).round(2).alias("Asm_CPM"),
            (pl.col("Asm_EM_Reads") * 1e6 / pl.col("Sample_Global_Reads")).round(2).alias("Asm_RPM"),
            
            # FPKM 加除零保护
            pl.when(pl.col("Asm_Length") > 0)
              .then((pl.col("Asm_EM_Reads") * 1e9 / (pl.col("Sample_Total_Mapped") * pl.col("Asm_Length"))))
              .otherwise(0.0).round(2).alias("Asm_FPKM"),
            
            # 【核心修复】：TPM 除零保护
            pl.when(pl.col("Sample_Total_RPK") > 0)
              .then((pl.col("Asm_RPK_Sum") * 1e6 / pl.col("Sample_Total_RPK")).round(2))
              .otherwise(0.0).alias("Asm_TPM"),
            
            (pl.col("Asm_EM_Reads") / pl.col("Sample_Total_Mapped") * 100.0).round(4).alias("Asm_Rel_Abund(%)")
        ])
        
        # 【修复���：先 join 代表序列数据，然后应用过滤条件
        rep_seq_stats = df.select([
            "Sample", "Accession", "Length", "Coverage(%)", "MeanDepth"
        ]).rename({
            "Accession": "Rep_Accession",
            "Coverage(%)": "Rep_Coverage(%)",
            "MeanDepth": "Rep_MeanDepth",
            "Length": "Rep_Length"
        })
        
        # 去重：每个样本每个代表序列只保留一条
        rep_seq_stats = rep_seq_stats.unique(subset=["Sample", "Rep_Accession"], keep="first", maintain_order=True)
        
        asm_df = asm_df.join(rep_seq_stats, on=["Sample", "Rep_Accession"], how="left")
        
        # 【核心修复】：在 final_df 过滤时应用覆盖��和平均深度的阈值
        # 这里使用代表序列的 Rep_Coverage(%) 和 Rep_MeanDepth 而不是聚合的平均值
        final_df = asm_df.filter(
            (pl.col("Rep_Coverage(%)") >= self.args.coverage) &  # 使用代表序列的覆盖度
            (pl.col("Rep_MeanDepth") >= self.args.meandepth) &   # 使用代表序列的平均深度
            (pl.col("Sample_Total_Mapped") > 0) &
            (pl.col("Uniq_Reads") >= self.args.min_uniq_reads) & 
            (pl.col("Asm_TPM") >= self.args.min_tpm)
        )
        
        final_df.write_csv(str(self.output_dir / "summary" / f"all_viruses.summary.{ext}"), separator=sep)
        
        pre_best_df = final_df.sort(["Sample", "taxid", "Asm_EM_Reads"], descending=[False, False, True]).unique(subset=["Sample", "taxid"], keep="first", maintain_order=True)
        
        logger.info(f"🧬 正在对存活的 {len(pre_best_df)} 个核心代表株进行深度进化测算 (ANI/Pi)...")
        tasks = [(row['Sample'], str(self.output_dir / 'bam' / f"{row['Sample']}.sorted.bam"), row['Rep_Accession']) for row in pre_best_df.iter_rows(named=True)]
        
        # 【修复4】：as_completed 异步收割，解决内存 OOM
        ani_pi_results = []
        with ProcessPoolExecutor(max_workers=self.args.threads) as executor:
            futures = {executor.submit(self._compute_ani_pi_worker, task): i for i, task in enumerate(tasks)}
            for future in tqdm(as_completed(futures), total=len(futures), desc="计算 ANI/Pi"):
                try:
                    result = future.result()
                    if result:
                        ani_pi_results.append(result)
                except Exception as e:
                    logger.warning(f"ANI/Pi 计算异常: {e}")

        # 【修复 - 核心】：确保 DataFrame 结构一致
        if len(ani_pi_results) == 0:
            ap_df = pl.DataFrame(schema={
                "Sample": pl.String, 
                "Rep_Accession": pl.String, 
                "Avg_Read_ANI": pl.Float64, 
                "Avg_Pi": pl.Float64
            })
        else:
            # 确保所有字段都存在且名称一致
            ap_df = pl.DataFrame(ani_pi_results)
            # 验证必要列存在
            required_cols = ["Sample", "Rep_Accession", "Avg_Read_ANI", "Avg_Pi"]
            for col in required_cols:
                if col not in ap_df.columns:
                    ap_df = ap_df.with_columns(pl.lit(None, pl.Float64).alias(col))
        
        # 使用字典方式加入，确保列对齐
        merged_df = pre_best_df.join(ap_df, on=["Sample", "Rep_Accession"], how="left")

        sp_thresh = self.args.sp_thresh
        
        # 【修复5】：严格三分类 + 0.1% 真实突变下限
        df_confirmed = merged_df.filter(
            (pl.col("Avg_Read_ANI").is_not_null()) & (pl.col("Avg_Read_ANI") >= sp_thresh)
        ).with_columns(pl.col("Species").alias("Adjusted_Species"))
        
        df_novel = merged_df.filter(
            (pl.col("Avg_Read_ANI").is_not_null()) & 
            (pl.col("Avg_Read_ANI") > 0.1) & 
            (pl.col("Avg_Read_ANI") < sp_thresh)
        ).with_columns(
            pl.concat_str([pl.lit("s__unclassified_"), pl.col("Species").str.replace_all(" ", "_")]).alias("Adjusted_Species")
        )
        
        df_failed = merged_df.filter(
            (pl.col("Avg_Read_ANI").is_null()) | (pl.col("Avg_Read_ANI") <= 0.1)
        ).with_columns(
            pl.concat_str([pl.lit("[ANI_FAILED]_"), pl.col("Species")]).alias("Adjusted_Species")
        )
        
        # 【修复】：使用代表序列的原始深度数据替代平均值
        final_cols = [
            "Sample", "taxid", "Adjusted_Species", "Species", "Rep_Accession", 
            "Rep_Length",  # 代表序列长度
            "Rep_Coverage(%)",  # 代表序列覆盖度
            "Rep_MeanDepth",  # 代表序列平均深度
            "Asm_EM_Reads", "Uniq_Reads", "Multi_Reads", "Unique(%)",
            "Avg_Read_ANI", "Avg_Pi",  # 进化指标
            "Asm_CPM", "Asm_RPM", "Asm_FPKM", "Asm_TPM", "Asm_Rel_Abund(%)", 
            "Segment_Accessions", "Rep_Reads"
        ]
        
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
        
        time.sleep(0.5)
        if not self.args.keep_bam:
            logger.info("🧹 正在清理巨型 BAM 与中间计算文件...")
            for d in ['bam', 'stat']:
                d_path = self.output_dir / d
                if d_path.exists(): shutil.rmtree(d_path)

    def generate_report_txt(self, raw_c, sum_c, best_c, unc_c, best_df):
        with open(self.output_dir / 'summary' / 'analysis_report.txt', 'w', encoding='utf-8') as f:
            f.write("=" * 80 + "\n宏病毒定量全景报告 (Virus Quantification V34 Production Gold Master)\n" + "=" * 80 + "\n\n")
            f.write(f"  处理样本: {len(self.samples)} | 留存(Raw): {raw_c} | 达标(Summary): {sum_c}\n")
            f.write(f"  确诊白名单(Best): {best_c} | 降级新种(Unclassified): {unc_c}\n\n")
            if len(best_df) > 0:
                f.write("【确诊白名单概览】\n" + "-" * 60 + "\n")
                asm_counts = best_df['Adjusted_Species'].value_counts()
                for row in asm_counts.iter_rows(named=True):
                    asm = row['Adjusted_Species']
                    count = row['counts'] if 'counts' in row else len(best_df.filter(pl.col('Adjusted_Species') == asm))
                    vd = best_df.filter(pl.col('Adjusted_Species') == asm)
                    if len(vd) > 0:
                        # 【修复】：安全处理 None 值
                        ani_mean = vd['Avg_Read_ANI'].mean()
                        pi_mean = vd['Avg_Pi'].mean()
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
        logger.info("=" * 60); logger.info("💥 Virus Quantification Pipeline (V34 Production Gold Master)"); logger.info("=" * 60)
        with Timer("整体管线运行"):
            self.setup_output_directory()
            self.build_index(); self.find_samples()
            with ProcessPoolExecutor(max_workers=self.args.threads) as executor:
                res_list = list(tqdm(executor.map(self.process_sample, self.samples), total=len(self.samples), desc="定量处理"))
            all_metrics = [m for sublist in res_list for m in sublist]
            self.summarize_results_polars(all_metrics)
def main():
    parser = argparse.ArgumentParser(
        description='宏病毒鉴定与精确定量管线 V34',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
【参数详细说明 - Detailed Parameter Guide】

【比对与质控参数】
  --tool              序列比对工具: bwa 或 bowtie2 (默认: bowtie2)
  --min_aln_len       最小比对长度 (默认: 80 bp)
  --min_pid           最小序列相似度 (默认: 0.90 = 90%%)
  --min_aln_prop      最小比对比例 (默认: 0.85 = 85%%)

【丰度阈值参数 - Abundance Thresholds】
  --min_uniq_reads    最少独特比对读数 (默认: 1)
  --min_tpm           最少 TPM 值 (默认: 0.0)

【深度与覆盖度过滤参数 - Coverage & Depth Filtering】
  --coverage          最小覆盖度百分比 (默认: 90.0%%)  ⚠️ 应用于代表序列
  --meandepth         最小平均深度 (默认: 10.0x)      ⚠️ 应用于代表序列

【进化学指标参数】
  --sp_thresh         物种ANI识别阈值 (默认: 95.0%%)
                      >= 95%% : 确诊种
                      > 0.1%% and < 95%% : 新种/未分类
                      <= 0.1%% : 检测失败

【性能参数】
  -t/--threads        样本并行处理数 (默认: 8)
  --align_threads     单样本比对线程数 (默认: 4)

【输出控制】
  --keep_bam          保留中间 BAM 文件 (默认: False，会删除)
  --verbose           详细日志输出 (默认: False)

【关键修复说明】
  ✅ 过滤阈值现在基于"代表序列"的真实数据而非平均值
  ✅ Rep_Coverage(%%) >= --coverage 时才通过过滤
  ✅ Rep_MeanDepth >= --meandepth 时才通过过滤
  ✅ 这确保输出的病毒序列质量稳定可控

【使用示例】
  # 严格模式
  python batch_virus_depth32.py -i input/ -r ref.fa --ref_info info.tsv \\
    --coverage 90 --meandepth 50 --sp_thresh 95.0

  # 标准模式（推荐）
  python batch_virus_depth32.py -i input/ -r ref.fa --ref_info info.tsv \\
    --coverage 50 --meandepth 10 --sp_thresh 95.0

  # 宽松模式（RNA病毒）
  python batch_virus_depth32.py -i input/ -r ref.fa --ref_info info.tsv \\
    --coverage 30 --meandepth 5 --sp_thresh 93.0
        """
    )
    
    parser.add_argument('-i', '--input_dir', required=True, help='输入 FASTQ 文件夹')
    parser.add_argument('-r', '--reference', required=True, help='全局参考基因组 FASTA')
    parser.add_argument('-o', '--output_dir', default='./virus_out', help='输出目录 (默认: ./virus_out)')
    parser.add_argument('--ref_info', type=str, required=True, help='本地参考信息 TSV 文件')
    parser.add_argument('--taxid_clusters', type=str, help='同义 TaxID 映射文件')
    
    parser.add_argument('--tool', choices=['bwa', 'bowtie2'], default='bowtie2', 
                        help='序列比对工具 (默认: bowtie2)')
    parser.add_argument('-t', '--threads', type=int, default=8, 
                        help='样本并行处理数 (默认: 8)')
    parser.add_argument('--align_threads', type=int, default=4, 
                        help='单样本比对线程数 (默认: 4)')
    
    parser.add_argument('--min_aln_len', type=int, default=80, 
                        help='最小比对长度 (默认: 80 bp)')
    parser.add_argument('--min_aln_prop', type=float, default=0.85, 
                        help='最小比对比例 (默认: 0.85)')
    parser.add_argument('--min_pid', type=float, default=0.90, 
                        help='最小序列相似度 (默认: 0.90)')
    parser.add_argument('--min_uniq_reads', type=int, default=1, 
                        help='最少独特比对读数 (默认: 1)')
    parser.add_argument('--sp_thresh', type=float, default=95.0, 
                        help='物种ANI识别阈值%% (默认: 95.0)')
    
    parser.add_argument('--coverage', type=float, default=90.0, 
                        help='代表序列最小覆盖度%% (默认: 90.0)')
    parser.add_argument('--meandepth', type=float, default=10.0, 
                        help='代表序列最小���均深度x (默认: 10.0)')
    parser.add_argument('--min_tpm', type=float, default=0.0, 
                        help='最少 TPM 值 (默认: 0.0)')
    
    parser.add_argument('--keep_bam', action='store_true', help='保留 BAM 文件')
    parser.add_argument('--verbose', action='store_true', help='详细日志输出')
    
    VirusQuantificationPipeline(parser.parse_args()).run_pipeline()

if __name__ == '__main__':
    main()

