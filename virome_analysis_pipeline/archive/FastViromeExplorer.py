#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
FastViromeExplorer Pro - Pseudo-Alignment V38 Edition (Enhanced Abundance)
终极极速版：流式排序 BAM + Pandepth 物理深度 + Polars 极速聚合
规范列序: Sample Accession Taxonomy Taxid Length CoveredSite TotalDepth Coverage(%) MeanDepth Reads TPM CPM RPM FPKM Rel_Abund(%)
新增：日志分离记录 + 样本动静态统计 + Tqdm无缝兼容 + Sorted BAM 提速
"""

import sys
import os
import subprocess
import argparse
import logging
import re
import time
import datetime
import gzip
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed

try:
    from tqdm import tqdm
except ImportError:
    tqdm = None

try:
    import polars as pl
except ImportError:
    print("❌ 致命错误: 必须安装 polars. 运行极速引擎请先执行: pip install polars")
    sys.exit(1)

# 可选的绘图库
try:
    import matplotlib
    matplotlib.use('Agg') 
    import matplotlib.pyplot as plt
    import seaborn as sns
except ImportError:
    plt = None
    sns = None

# ==========================================
# 兼容 Tqdm 的自定义日志处理器
# ==========================================
class TqdmLoggingHandler(logging.Handler):
    def emit(self, record):
        try:
            msg = self.format(record)
            if tqdm is not None:
                tqdm.write(msg)
            else:
                sys.stdout.write(msg + '\n')
            sys.stdout.flush()
        except Exception:
            self.handleError(record)

logger = logging.getLogger("Pseudo-V38")

def setup_logging(out_root):
    logs_dir = out_root / "logs"
    logs_dir.mkdir(parents=True, exist_ok=True)
    
    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    log_file = logs_dir / f"pipeline_run_{timestamp}.log"
    
    logger.setLevel(logging.INFO)
    logger.handlers = [] # 清除默认 handlers
    
    # 文件日志 (详细记录)
    fh = logging.FileHandler(log_file, mode='a', encoding='utf-8')
    fh.setLevel(logging.INFO)
    fh.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
    logger.addHandler(fh)
    
    # 屏幕控制台日志 (兼容 Tqdm)
    ch = TqdmLoggingHandler()
    ch.setLevel(logging.INFO)
    ch.setFormatter(logging.Formatter('%(asctime)s - %(message)s', "%H:%M:%S"))
    logger.addHandler(ch)
    
    return logs_dir

# ==========================================
# 自动建库与智能索引推断模块
# ==========================================
def build_index_if_needed(config):
    ref_path = Path(config.reference).resolve()
    
    if not ref_path.exists():
        logger.error(f"❌ 提供的参考 FASTA 不存在: {config.reference}")
        sys.exit(1)
        
    # 智能推断索引的存放路径 (与 FASTA 序列同级目录)
    if config.use_salmon:
        idx_path = ref_path.with_suffix('.salmon_idx')
    else:
        idx_path = ref_path.with_suffix('.kallisto_idx')
        
    config.index_file = str(idx_path)
    index_exists = False
    
    # 智能检查索引是否已完整存在
    if config.use_salmon:
        if idx_path.exists() and idx_path.is_dir() and any(idx_path.iterdir()):
            index_exists = True
    else:
        if idx_path.exists() and idx_path.is_file():
            index_exists = True
            
    if index_exists:
        logger.info(f"✅ [数据库检测] 寻址成功: 检测到可用索引库 -> {config.index_file}，跳过建库。")
        return
        
    # 组装并执行建库命令
    index_log_file = config.logs_dir / "index_build.log"
    logger.info(f"🔨 [数据库检测] 未检测到索引，开始使用 {ref_path.name} 全自动构建至: {idx_path.name}")
    logger.info(f"📁 [建库日志] 详细建库日志重定向至: {index_log_file}")
    
    total_threads = config.parallel * config.threads_per_sample
    
    if config.use_salmon:
        cmd = f"salmon index -k 31 -i {config.index_file} -t {config.reference} --threads {total_threads}"
    else:
        cmd = f"kallisto index -k 31 -i {config.index_file} --threads {total_threads} {config.reference}"
        
    try:
        with open(index_log_file, "w") as flog:
            flog.write(f"[{datetime.datetime.now()}] 执行建库命令: {cmd}\n{'='*60}\n")
            flog.flush()
            subprocess.run(cmd, shell=True, check=True, stdout=flog, stderr=subprocess.STDOUT)
        logger.info("✅ 索引自动构建成功！即将进入比对流程...")
    except subprocess.CalledProcessError as e:
        logger.error(f"❌ 索引构建失败，进程退出码: {e.returncode}。请查看 {index_log_file} 获取详情。")
        sys.exit(1)

# ==========================================
# 核心工作类：单样本计算流程
# ==========================================
class PseudoViromeWorker:
    def __init__(self, config, sample_info):
        self.config = config
        self.sample_id = sample_info['id']
        self.read1 = sample_info['r1']
        self.read2 = sample_info['r2']
        self.out_dir = Path(sample_info['out_dir']).resolve()
        self.index_file = config.index_file
        self.out_dir.mkdir(parents=True, exist_ok=True)
        # 单样本运行详细日志定向到统一的 logs 目录
        self.sample_log = self.config.logs_dir / f"sample_{self.sample_id}.log"

    def run_command(self, command):
        with open(self.sample_log, "a") as log_file:
            log_file.write(f"\n[{datetime.datetime.now()}] EXEC: {command}\n")
            log_file.flush()
            subprocess.run(command, shell=True, check=True, stderr=log_file, stdout=log_file)

    def get_average_read_length(self):
        is_fasta = any(self.read1.lower().endswith(ext) for ext in ['.fa', '.fasta', '.fa.gz', '.fasta.gz'])
        awk_script = "'/^>/ {if (seqlen > 0) { lenSum+=seqlen; readCount++; seqlen=0 }} !/^>/ { seqlen += length($0) } END { if (seqlen > 0) { lenSum+=seqlen; readCount++; }; if (readCount > 0) print lenSum/readCount; else print 0; }'" if is_fasta else "'NR%4 == 2 {lenSum+=length($0); readCount++;} END {if (readCount > 0) print lenSum/readCount; else print 0}'"
        cmd = f"gzip -dc {self.read1} | awk {awk_script}" if self.read1.endswith(".gz") else f"awk {awk_script} {self.read1}"
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        try:
            return float(result.stdout.strip())
        except ValueError:
            return 150.0

    def run_pseudo_alignment_and_depth(self):
        bam_output = self.out_dir / "reads_mapped.bam"
        threads = self.config.threads_per_sample
        
        # 🔥 优化：流式转 BAM 并排序，极大提升写入和 pandepth 读取速度
        samtools_pipe = f"samtools view -@ {threads} -b -F 0x04 - | samtools sort -@ {threads} -o {bam_output}"
        
        if self.config.use_salmon:
            read_arg = f"-1 {self.read1} -2 {self.read2}" if self.read2 else f"-r {self.read1}"
            quant_cmd = f"salmon quant -i {self.index_file} -p {threads} -l A {read_arg} -o {self.out_dir} --writeMappings"
            full_cmd = f"{quant_cmd} | {samtools_pipe}"
        else:
            quant_args = f"--pseudobam {self.read1} {self.read2}" if self.read2 else f"--single -l 200 -s 50 --pseudobam {self.read1}"
            full_cmd = f"kallisto quant -i {self.index_file} -t {threads} -o {self.out_dir} {quant_args} | {samtools_pipe}"
        
        self.run_command(full_cmd)

        # Pandepth 现读取排序后的 BAM，速度更快
        pandepth_prefix = self.out_dir / "pandepth"
        pandepth_cmd = f"pandepth -i {bam_output} -o {pandepth_prefix} -a -t {threads}"
        self.run_command(pandepth_cmd)

    def generate_parquet(self, avg_read_len):
        sf = "quant.sf" if self.config.use_salmon else "abundance.tsv"
        quant_path = self.out_dir / sf
        if not quant_path.exists(): return
        
        if self.config.use_salmon:
            df_quant = pl.read_csv(quant_path, separator='\t').select([
                pl.col("Name").alias("Accession"),
                pl.col("Length"),
                pl.col("TPM"),
                pl.col("NumReads").alias("Reads")
            ])
        else:
            df_quant = pl.read_csv(quant_path, separator='\t').select([
                pl.col("target_id").alias("Accession"),
                pl.col("length").alias("Length"),
                pl.col("tpm").alias("TPM"),
                pl.col("est_counts").alias("Reads")
            ])
            
        df_quant = df_quant.filter(pl.col("Reads") > 0)

        # 安全解析 Pandepth，过滤尾部脏数据
        cov_path = self.out_dir / "pandepth.chr.stat.gz"
        cov_data = []
        if cov_path.exists():
            try:
                with gzip.open(cov_path, 'rt') as f:
                    for line in f:
                        if line.startswith('#') or not line.strip(): 
                            continue
                        p = line.split('\t')
                        if len(p) >= 6:
                            try:
                                cov_data.append({
                                    "Accession": p[0].strip(),
                                    "CoveredSite": int(p[2]),
                                    "TotalDepth": float(p[3]),
                                    "Coverage(%)": float(p[4]),
                                    "MeanDepth": float(p[5])
                                })
                            except ValueError:
                                pass
            except Exception as e:
                logger.warning(f"样本 {self.sample_id} 深度解析异常: {e}")

        if cov_data:
            df_cov = pl.DataFrame(cov_data)
        else:
            df_cov = pl.DataFrame({"Accession": [], "CoveredSite": [], "TotalDepth": [], "Coverage(%)": [], "MeanDepth": []},
                                 schema={"Accession": pl.Utf8, "CoveredSite": pl.Int64, "TotalDepth": pl.Float64, "Coverage(%)": pl.Float64, "MeanDepth": pl.Float64})

        df_merged = df_quant.join(df_cov, on="Accession", how="left").fill_null(0.0)
        df_merged = df_merged.with_columns([
            pl.lit(self.sample_id).alias("Sample"),
            pl.lit(avg_read_len).alias("Avg_Read_Len")
        ])
        
        df_merged.write_parquet(self.out_dir / f"{self.sample_id}.batch.parquet")

    def cleanup(self):
        if self.config.keep_tmp: return
        # 优化：清理 .bam 替代了 .sam
        for ext in ["reads_mapped.bam", "pandepth.SiteDepth.gz", "quant.sf", "abundance.tsv", "abundance.h5"]:
            p = self.out_dir / ext
            if p.is_file(): p.unlink()

    def run(self):
        st = time.time()
        batch_file = self.out_dir / f"{self.sample_id}.batch.parquet"
        
        if self.config.resume and batch_file.exists():
            return (True, self.sample_id, "Skipped", "0s")

        try:
            # 清理之前的日志文件以便重新记录
            if self.sample_log.exists():
                self.sample_log.unlink()
                
            avg_len = self.get_average_read_length()
            self.run_pseudo_alignment_and_depth()
            self.generate_parquet(avg_len)
            self.cleanup()
            return (True, self.sample_id, "Success", str(datetime.timedelta(seconds=int(time.time()-st))))
        except Exception as e:
            return (False, self.sample_id, str(e), "0")

# ==========================================
# 多进程安全包装函数
# ==========================================
def run_pseudo_worker(config, sample_info):
    worker = PseudoViromeWorker(config, sample_info)
    return worker.run()

# ==========================================
# 目录解析逻辑
# ==========================================
def scan_directory_for_samples(dir_path):
    path_obj = Path(dir_path)
    if not path_obj.exists() or not path_obj.is_dir(): return []
    
    file_list = sorted(list(path_obj.glob("*.fastq*")) + list(path_obj.glob("*.fq*")) + list(path_obj.glob("*.fa*")) + list(path_obj.glob("*.fasta*")))
    samples = []
    processed_files = set()
    pattern = re.compile(r"^(.*?)(_R1|-R1|_1|\.1)([\._].*)$")
    
    for f in file_list:
        if str(f) in processed_files: continue
        match = pattern.match(f.name)
        if match:
            prefix, r1_tag, suffix = match.groups()
            r2_tag = r1_tag.replace('1', '2')
            r2_path = f.with_name(f"{prefix}{r2_tag}{suffix}")
            
            if r2_path.exists():
                samples.append({'id': prefix.rstrip('_-.='), 'r1': str(f), 'r2': str(r2_path)})
                processed_files.update([str(f), str(r2_path)])
                continue
        
        clean_id = f.name
        for ext in ['.fastq.gz', '.fq.gz', '.fasta.gz', '.fa.gz', '.fastq', '.fq', '.fasta', '.fa']:
            if clean_id.lower().endswith(ext):
                clean_id = clean_id[:-len(ext)]
                break
        samples.append({'id': clean_id, 'r1': str(f), 'r2': None})
        processed_files.add(str(f))
        
    return samples

# ==========================================
# 全局 V38 聚合引擎
# ==========================================
def generate_global_summary(config, out_root, total_samples_count):
    logger.info("📊 开始 Polars V38 聚合引擎 (节段物种聚合 + 双轨制过滤 + 规范排版)...")
    summary_dir = out_root / "summary"
    summary_dir.mkdir(parents=True, exist_ok=True)
    
    parquet_files = list(out_root.rglob("*.batch.parquet"))
    if not parquet_files:
        logger.warning("未发现分析结果，全局汇总取消。")
        return
        
    lazy_dfs = [pl.scan_parquet(f) for f in parquet_files]
    df = pl.concat(lazy_dfs)
    
    sample_totals = df.group_by("Sample").agg(
        pl.col("Reads").sum().alias("Sample_Total_Mapped")
    )
    df = df.join(sample_totals.lazy(), on="Sample", how="left")
    
    df = df.with_columns([
        (1.0 - (-(pl.col("Reads") * pl.col("Avg_Read_Len")) / pl.col("Length")).exp()).alias("Predicted_Support")
    ])
    
    df = df.with_columns(
        pl.when(pl.col("Predicted_Support") > 0)
          .then((pl.col("Coverage(%)") / 100.0) / pl.col("Predicted_Support"))
          .otherwise(0.0)
          .alias("Poisson_Ratio")
    )
    
    cond_a = (
        (pl.col("Reads") >= config.min_reads) & 
        (pl.col("Coverage(%)") / 100.0 >= config.coverage) & 
        (pl.col("Poisson_Ratio") >= config.ratio)
    )
    
    cond_b = pl.lit(False)
    if config.genes_cov_dict:
        df_genes = pl.DataFrame([
            {"Accession": k, "g_tot": v['gene_total_cov'], "g_avr": v['gene_avr_cov']} 
            for k, v in config.genes_cov_dict.items()
        ])
        df = df.join(df_genes.lazy(), on="Accession", how="left").fill_null(0.0)
        cond_b = (
            (pl.col("Reads") >= config.min_reads) &
            (pl.col("g_tot") >= config.min_gene_total_cov) &
            (pl.col("g_avr") >= config.min_gene_avr_cov)
        )
    
    df = df.with_columns(
        pl.when(cond_a & ~cond_b).then(pl.lit("Track_A"))
          .when(~cond_a & cond_b).then(pl.lit("Track_B"))
          .when(cond_a & cond_b).then(pl.lit("Track_A+B"))
          .otherwise(pl.lit("Filtered_Out"))
          .alias("Filter_Track")
    )
    
    df_collected = df.collect()
    
    logger.info(f"🔗 正在合并参考库分类学与节段元数据: {config.ref_info}")
    ref_meta = pl.read_csv(config.ref_info, separator='\t', ignore_errors=True)
    
    acc_col = next((c for c in ['Accession', 'accession', 'ID', 'seqid'] if c in ref_meta.columns), None)
    if acc_col:
        ref_meta = ref_meta.with_columns(pl.col(acc_col).cast(pl.Utf8).str.strip_chars())
        df_collected = df_collected.with_columns(pl.col("Accession").cast(pl.Utf8).str.strip_chars())
        
        tax_col = next((c for c in ['Species_NCBI', 'Species_ICTV', 'Species', 'Taxonomy'] if c in ref_meta.columns), None)
        taxid_col = next((c for c in ['Taxid', 'taxid', 'TaxID', 'taxID'] if c in ref_meta.columns), None)
        seg_col = next((c for c in ['Segment', 'segment'] if c in ref_meta.columns), None)
        
        keep_cols = [acc_col]
        if tax_col: keep_cols.append(tax_col)
        if taxid_col: keep_cols.append(taxid_col)
        if seg_col: keep_cols.append(seg_col)
        
        ref_meta_subset = ref_meta.select(keep_cols)
        
        rename_mapping = {}
        if tax_col and tax_col != "Taxonomy": rename_mapping[tax_col] = "Taxonomy"
        if taxid_col and taxid_col != "Taxid": rename_mapping[taxid_col] = "Taxid"
        if seg_col and seg_col != "Segment": rename_mapping[seg_col] = "Segment"
        if rename_mapping:
            ref_meta_subset = ref_meta_subset.rename(rename_mapping)
            
        df_collected = df_collected.join(ref_meta_subset, left_on="Accession", right_on=acc_col, how="left")
        
    if "Taxonomy" not in df_collected.columns: df_collected = df_collected.with_columns(pl.lit("NA").alias("Taxonomy"))
    else: df_collected = df_collected.with_columns(pl.col("Taxonomy").fill_null("NA"))
        
    if "Taxid" not in df_collected.columns: df_collected = df_collected.with_columns(pl.lit("NA").alias("Taxid"))
    else: df_collected = df_collected.with_columns(pl.col("Taxid").cast(pl.Utf8).fill_null("NA"))
        
    if "Segment" not in df_collected.columns: df_collected = df_collected.with_columns(pl.lit("").alias("Segment"))
    else: df_collected = df_collected.with_columns(pl.col("Segment").cast(pl.Utf8).fill_null(""))

    raw_out = summary_dir / "all_viruses.raw.tsv"
    df_collected.write_csv(raw_out, separator='\t')

    df_best = df_collected.filter(pl.col("Filter_Track") != "Filtered_Out")
    
    if len(df_best) == 0:
        logger.warning("⚠️ 所有数据均被双轨制过滤模型拦截，无高置信度阳性结果。")
        return

    # 节段聚合缝合逻辑
    df_best = df_best.with_columns([
        pl.when(pl.col("Segment") != "")
          .then(pl.col("Segment") + ":" + pl.col("Accession"))
          .otherwise(pl.col("Accession"))
          .alias("Seg_Acc_Str")
    ])

    df_agg = df_best.group_by(["Sample", "Taxid", "Taxonomy"]).agg([
        pl.col("Seg_Acc_Str").unique().str.join(" | ").alias("Segment_Accessions"),
        pl.col("Length").sum().alias("Total_Length"),
        pl.col("CoveredSite").sum().alias("Total_CoveredSite"),
        pl.col("TotalDepth").sum().alias("Sum_TotalDepth"),
        pl.col("Reads").sum().alias("Total_Reads"),
        pl.col("TPM").sum().alias("Total_TPM"),
        ((pl.col("MeanDepth") * pl.col("Length")).sum() / pl.col("Length").sum()).round(2).alias("Weighted_MeanDepth"),
        pl.col("Filter_Track").first().alias("Filter_Track"),
        pl.col("Poisson_Ratio").mean().round(4).alias("Avg_Poisson_Ratio"),
        pl.col("Sample_Total_Mapped").first().alias("Sample_Total_Mapped")
    ])

    df_agg = df_agg.with_columns([
        (pl.col("Total_CoveredSite") / pl.col("Total_Length") * 100.0).round(2).alias("Overall_Coverage(%)"),
        (pl.col("Total_Reads") * 1e6 / pl.col("Sample_Total_Mapped")).round(2).alias("Asm_CPM"),
        (pl.col("Total_Reads") * 1e6 / pl.col("Sample_Total_Mapped")).round(2).alias("Asm_RPM"),
        (pl.col("Total_Reads") * 1e9 / (pl.col("Total_Length") * pl.col("Sample_Total_Mapped"))).round(2).alias("Asm_FPKM"),
        (pl.col("Total_Reads") / pl.col("Sample_Total_Mapped") * 100.0).round(4).alias("Asm_Rel_Abund(%)")
    ])

    df_agg = df_agg.rename({
        "Segment_Accessions": "Accession",
        "Total_Length": "Length",
        "Total_CoveredSite": "CoveredSite",
        "Sum_TotalDepth": "TotalDepth",
        "Overall_Coverage(%)": "Coverage(%)",
        "Weighted_MeanDepth": "MeanDepth",
        "Total_Reads": "Reads",
        "Total_TPM": "TPM",
        "Avg_Poisson_Ratio": "Poisson_Ratio"
    })

    target_base_cols = [
        "Sample", "Accession", "Taxonomy", "Taxid", "Length", 
        "CoveredSite", "TotalDepth", "Coverage(%)", "MeanDepth", 
        "Reads", "TPM", "Asm_CPM", "Asm_RPM", "Asm_FPKM", "Asm_Rel_Abund(%)", 
        "Poisson_Ratio", "Filter_Track"
    ]
    
    base_cols = [c for c in target_base_cols if c in df_agg.columns]
    df_final = df_agg.select(base_cols)
    df_final = df_final.sort(by=["Sample", "TPM"], descending=[False, True])
    
    best_out = summary_dir / "all_viruses.best.summary.tsv"
    df_final.write_csv(best_out, separator='\t')
    
    try:
        matrix_tpm = df_final.pivot(index=["Taxonomy", "Taxid"], on="Sample", values="TPM", aggregate_function="sum").fill_null(0.0)
        matrix_tpm.write_csv(summary_dir / "Virus_TPM_Matrix.tsv", separator='\t')
        
        matrix_reads = df_final.pivot(index=["Taxonomy", "Taxid"], on="Sample", values="Reads", aggregate_function="sum").fill_null(0.0)
        matrix_reads.write_csv(summary_dir / "Virus_Reads_Matrix.tsv", separator='\t')
        
        pd_best = df_final.to_pandas()
        with open(summary_dir / "Virus_Sample_Distribution.txt", "w") as f:
            f.write("#Taxonomy\tTaxid\tSamples\n")
            for (name, tid), group in pd_best.groupby(['Taxonomy', 'Taxid']):
                f.write(f"{name}\t{tid}\t{' '.join(sorted(group['Sample'].unique()))}\n")
                
        with open(summary_dir / "analysis_report.txt", "w") as f:
            f.write(f"================================================================================\n")
            f.write(f"FastViromeExplorer Pro 病毒丰度核心分析报告 (V38 极速增强版)\n")
            f.write(f"================================================================================\n\n")
            f.write(f"分析时间: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"总处理样本数: {total_samples_count}\n\n")
            f.write(f"--- 核心确诊群落概览 (Top 流行病原体) ---\n")
            f.write(f"------------------------------------------------------------\n")
            stats = pd_best.groupby(['Taxonomy', 'Taxid']).size().reset_index(name='Freq').sort_values('Freq', ascending=False)
            for i, r in stats.head(30).iterrows():
                vd = pd_best[(pd_best['Taxonomy'] == r['Taxonomy']) & (pd_best['Taxid'] == r['Taxid'])]
                f.write(f"{r['Taxonomy']} (TaxID: {r['Taxid']}):\n")
                f.write(f"  样本检出频次: {r['Freq']} / {total_samples_count} (检出率 {(r['Freq']/total_samples_count)*100:.1f}%)\n")
                f.write(f"  阳性样本平均 Reads: {vd['Reads'].mean():.1f}\n")
                f.write(f"  阳性样本平均 TPM:   {vd['TPM'].mean():.2f}\n")
                if 'Asm_CPM' in vd.columns:
                    f.write(f"  阳性样本平均 CPM:   {vd['Asm_CPM'].mean():.2f}\n")
                f.write("\n")

        if plt and sns:
            top_plot = stats.head(30).copy()
            top_plot['Perc'] = (top_plot['Freq'] / total_samples_count) * 100
            plt.figure(figsize=(12, 8))
            sns.barplot(data=top_plot, x='Perc', y='Taxonomy', color='c')
            plt.title(f"Top Virus Prevalence (Total Samples: {total_samples_count})", fontsize=15)
            plt.xlabel("Prevalence (% of Samples)", fontsize=12)
            plt.tight_layout()
            plt.savefig(summary_dir / "Virus_Prevalence_Plot.png", dpi=300)

    except Exception as e:
        logger.error(f"矩阵或图表生成中出现部分失败: {e}")

    logger.info(f"✨ 汇总报表及矩阵生成完毕！核心报告位于: {summary_dir.name}/")


def main():
    parser = argparse.ArgumentParser(
        description="FastViromeExplorer Pro: 批量并发病毒检测与全自动定量汇总管线 (V38极速增强版)",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    io_group = parser.add_argument_group('📥 输入与输出')
    io_group.add_argument("--input_dirs", nargs='+', required=True, help="输入目录列表，自动扫描 fastq/fasta")
    io_group.add_argument("-o", "--out_dir", required=True, help="输出主目录")
    
    db_group = parser.add_argument_group('🗄️ 数据库与索引')
    db_group.add_argument("-r", "--reference", required=True, help="参考基因组 FASTA。脚本将自动在其同目录下查找或构建索引库")
    db_group.add_argument("--ref_info", required=True, help="包含分类学注释的数据库元信息表 (如 ref_info.tsv)")
    
    perf_group = parser.add_argument_group('🚀 并发与性能')
    perf_group.add_argument("-p", "--parallel", type=int, default=1, help="并发样本处理数 [默认: 1]")
    perf_group.add_argument("-t", "--threads", type=int, default=4, help="单样本内部分配线程数 [默认: 4]")
    
    track_a_group = parser.add_argument_group('🛡️ 轨道A (全长基因组物理过滤)')
    track_a_group.add_argument("-cr", "--ratio", type=float, default=0.3, help="泊松分布 Coverage Ratio 阈值 [默认: 0.3]")
    track_a_group.add_argument("-co", "--coverage", type=float, default=0.1, help="全长覆盖率阈值 [默认: 0.1]")
    track_a_group.add_argument("-cn", "--min_reads", type=int, default=10, help="最低有效比对 Reads 数 [默认: 10]")
    
    track_b_group = parser.add_argument_group('🧬 轨道B (活跃转录区基因级挽救)')
    track_b_group.add_argument("--genes_cov", help="传入 genes_cov.tsv 即可激活该双轨挽救系统")
    track_b_group.add_argument("--min_gene_total_cov", type=float, default=80.0, help="最低转录区总覆盖比例(%%) [默认: 80.0]")
    track_b_group.add_argument("--min_gene_avr_cov", type=float, default=5.0, help="最低平均基因覆盖深度 [默认: 5.0]")

    misc_group = parser.add_argument_group('⚙️ 高级控制')
    misc_group.add_argument("-salmon", action="store_true", help="使用 Salmon 替代 Kallisto 定量")
    misc_group.add_argument("--resume", action="store_true", help="断点续跑：极速跳过已存在结果的样本")
    misc_group.add_argument("--keep_tmp", action="store_true", help="保留底层产生的无序 BAM 等庞大临时文件")

    args = parser.parse_args()
    args.use_salmon = args.salmon
    args.threads_per_sample = args.threads

    # 1. 初始化输出主目录并挂载全局 Logging 系统
    out_root = Path(args.out_dir).resolve()
    args.logs_dir = setup_logging(out_root)
    
    logger.info("="*60)
    logger.info("🚀 FastViromeExplorer Pro - 极速比对定量管线启动")
    logger.info("="*60)

    # 2. 检查与构建数据库 (日志被导流到日志文件)
    build_index_if_needed(args)

    # 3. 扫描目录与生成统计信息
    tasks = []
    for d in args.input_dirs: tasks.extend(scan_directory_for_samples(d))
    if not tasks:
        logger.error("❌ 未在输入目录中找到任何 FASTQ/FASTA 文件。")
        sys.exit(1)

    total_samples = len(tasks)
    se_count = sum(1 for t in tasks if t['r2'] is None)
    pe_count = sum(1 for t in tasks if t['r2'] is not None)

    logger.info("📊 [样本扫描统计]")
    logger.info(f"   ▶ 总计样本数 : {total_samples}")
    logger.info(f"   ▶ 单端数据(SE): {se_count}")
    logger.info(f"   ▶ 双端数据(PE): {pe_count}")
    logger.info("="*60)

    # 载入轨道B的基因深度表(如提供)
    args.genes_cov_dict = {}
    if args.genes_cov and os.path.exists(args.genes_cov):
        import pandas as pd
        try:
            gdf = pd.read_csv(args.genes_cov, sep='\t')
            args.genes_cov_dict = gdf.set_index('seqid')[['gene_total_cov', 'gene_avr_cov']].to_dict('index')
            logger.info("✅ 轨道 B (活跃转录区挽救) 已激活。")
        except Exception as e:
            logger.warning(f"⚠️ 解析 {args.genes_cov} 失败: {e}")

    samples_dir = out_root / "samples"
    for t in tasks: t['out_dir'] = str(samples_dir / t['id'])

    # 4. 多进程并发比对 (Tqdm进度条与状态跟踪)
    with ProcessPoolExecutor(max_workers=min(args.parallel, total_samples)) as exe:
        futures = {exe.submit(run_pseudo_worker, args, t): t['id'] for t in tasks}
        
        if tqdm: pbar = tqdm(total=total_samples, desc="比对与深度计算", dynamic_ncols=True, colour="green")
        
        for future in as_completed(futures):
            try:
                success, sid, msg, t_str = future.result()
                if success: 
                    logger.info(f"✔️ 样本 [{sid}] 运行完毕 | 状态: {msg} | 耗时: {t_str}")
                else:
                    logger.error(f"❌ 样本 [{sid}] 运行崩溃 | 异常: {msg}")
                    
                if tqdm: 
                    pbar.set_postfix_str(f"最新完成: {sid}")
                    pbar.update(1)
            except Exception as e:
                logger.error(f"并发调度器严重异常: {e}")
                
        if tqdm: pbar.close()

    # 5. 全局报表输出
    logger.info("="*60)
    generate_global_summary(args, out_root, total_samples)
    logger.info("🎉 全部流程执行完毕！")

if __name__ == "__main__":
    main()
