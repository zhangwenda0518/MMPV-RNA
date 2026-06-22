#!/usr/bin/env python3
"""
batch_virus_downstream.py — 病毒确诊后处理、节段病毒优化与变异分析终极管线
【VAP专属表头 | 真实测序深度全指标修正 | 盲区侦测 | Shannon/AF & SNPGenie 稳健版 | 🚀全局命令日志100%捕获版】
"""

import argparse
import logging
import math
import os
import re
import shutil
import subprocess
import sys
import time
import urllib.request
import gzip
import traceback
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime
from pathlib import Path

import colorlog
import numpy as np
import pandas as pd
import pysam
from tqdm import tqdm
import matplotlib.pyplot as plt
import seaborn as sns

# ==========================================
# 1. 全局日志与增强版系统调用工具
# ==========================================
def setup_logging(verbose: bool = False) -> logging.Logger:
    handler = colorlog.StreamHandler()
    handler.setFormatter(colorlog.ColoredFormatter(
        "%(log_color)s%(asctime)s - %(levelname)s - %(message)s",
        datefmt="%H:%M:%S",
        log_colors={"DEBUG": "cyan", "INFO": "green", "WARNING": "yellow", "ERROR": "red", "CRITICAL": "bold_red"},
    ))
    log = colorlog.getLogger("vap")
    log.handlers.clear()
    log.addHandler(handler)
    log.setLevel(logging.DEBUG if verbose else logging.INFO)
    return log

logger = setup_logging()

class Timer:
    def __init__(self, name: str = ""): self.name = name
    def __enter__(self):
        self._t = time.time(); logger.info(f"⏱️  开始: {self.name}"); return self
    def __exit__(self, exc_type, *_):
        dur = time.time() - self._t
        if exc_type is None: logger.info(f"✅ 完成: {self.name} [{self._fmt(dur)}]")
        else: logger.error(f"❌ 失败: {self.name} [{self._fmt(dur)}]")

    @staticmethod
    def _fmt(s: float) -> str:
        if s < 60: return f"{s:.1f} 秒"
        if s < 3600: return f"{s / 60:.1f} 分钟"
        return f"{s / 3600:.1f} 小时"

def safe_name(s: str, max_len: int = 100) -> str:
    s = str(s)
    s = re.sub(r'[^A-Za-z0-9\-.]', '_', s)
    s = re.sub(r'_+', '_', s)
    return s.strip('_.')[:max_len]

def run_cmd(cmd: str, log_path: str = None, master_log: str = None, check: bool = True):
    """
    🔥 超级增强版命令执行工具：
    1. 强制开启 pipefail 捕获管道错误
    2. 将命令本身、退出码、标准输出(STDOUT)、标准错误(STDERR) 全部写入独立的日志文件
    3. 同时追加写入到全局的 master_log 中，确保绝对不漏掉任何一条命令的执行细节
    """
    full_cmd = f"set -o pipefail; {cmd}"
    start_time = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    
    result = subprocess.run(full_cmd, shell=True, executable="/bin/bash", capture_output=True, text=True)
    
    # 拼装详细的日志内容
    log_content = f"\n[{start_time}] CMD: {cmd}\nEXIT_CODE: {result.returncode}\n"
    if result.stdout: log_content += f"--- STDOUT ---\n{result.stdout.strip()}\n"
    if result.stderr: log_content += f"--- STDERR ---\n{result.stderr.strip()}\n"
    log_content += "-" * 80 + "\n"
    
    # 写入单独的文件日志（如 SRR123.variants.log）
    if log_path:
        lp = Path(log_path)
        lp.parent.mkdir(parents=True, exist_ok=True)
        with open(lp, "a", encoding="utf-8") as f: f.write(log_content)
            
    # 同步写入全局大日志（master_commands.log）
    if master_log:
        ml = Path(master_log)
        ml.parent.mkdir(parents=True, exist_ok=True)
        with open(ml, "a", encoding="utf-8") as f: f.write(log_content)
            
    if check and result.returncode != 0:
        raise subprocess.CalledProcessError(result.returncode, cmd, output=result.stdout, stderr=result.stderr)
    return result


# ==========================================
# 2. 变异等位基因频率 & 虚拟 GTF 工具函数
# ==========================================
def extract_allele_frequency(vcf_path, out_tsv):
    """智能解析 VCF，提取等位基因突变频率 (AF) 保存为 TSV 用于画图"""
    if not Path(vcf_path).exists() or Path(vcf_path).stat().st_size == 0: return
    records = []
    with open(vcf_path, 'r') as f:
        for line in f:
            if line.startswith('#'): continue
            cols = line.strip().split('\t')
            if len(cols) < 8: continue
            chrom, pos, ref, alt, info = cols[0], cols[1], cols[3], cols[4], cols[7]
            
            freq = 0.0
            af_match = re.search(r'\bAF=([\d\.]+)', info)
            if af_match: 
                freq = float(af_match.group(1))
            else:
                ao_match = re.search(r'\bAO=([\d,]+)', info)
                dp_match = re.search(r'\bDP=([\d]+)', info)
                if ao_match and dp_match:
                    try:
                        ao = max([int(x) for x in ao_match.group(1).split(',')])
                        dp = int(dp_match.group(1))
                        if dp > 0: freq = ao / dp
                    except: pass
                    
            if freq > 0:
                records.append(f"{chrom}\t{pos}\t{ref}\t{alt}\t{freq:.4f}")
                
    if records:
        with open(out_tsv, 'w') as f:
            f.write("CHROM\tPOS\tREF\tALT\tALT_FREQ\n")
            f.write("\n".join(records) + "\n")

def generate_dummy_gtf(fasta_path, gtf_path):
    """如果NCBI没有提供CDS注释，通过FASTA长度生成虚拟全基因组GTF，拯救SNPGenie"""
    try:
        with open(fasta_path, 'r') as f:
            header = f.readline().strip()[1:].split()[0]
            seq = "".join([line.strip() for line in f])
        length = len(seq)
        with open(gtf_path, 'w') as f:
            f.write(f"{header}\tVAP_Fallback\tCDS\t1\t{length}\t.\t+\t0\tgene_id \"{header}_CDS\"; transcript_id \"{header}_TX\";\n")
    except Exception as e:
        logger.warning(f"生成虚拟 GTF 失败: {e}")

class VODKA2_PostProcessor:
    @staticmethod
    def run_blast_validation(results_txt: Path, ref_fasta: Path, dvg_type: str, n_shift: int = 5, evalue: float = 0.001):
        if not results_txt.exists() or results_txt.stat().st_size == 0: return None
        df = pd.read_csv(results_txt, sep='\t')
        if df.empty: return None
        b_fa, r_fa = results_txt.with_suffix('.B.fa'), results_txt.with_suffix('.R.fa')
        with open(b_fa, 'w') as fb, open(r_fa, 'w') as fr:
            for _, row in df.iterrows():
                ac, read_id = row['A_C(A_C)'], row['READ_ID']
                seq = str(row['MAPPED_ONLY']).replace('<', '').replace('>', '')
                if ':::' in seq:
                    b_seq, r_seq = seq.split(':::')
                    fb.write(f">{ac}-{read_id}\n{b_seq}\n")
                    fr.write(f">{ac}-{read_id}\n{r_seq}\n")
        blast_out_b, blast_out_r = results_txt.with_suffix('.B.blast'), results_txt.with_suffix('.R.blast')
        blast_cmd_base = f"blastn -max_hsps 1 -db '{ref_fasta}' -outfmt '6 qseqid qstart qend sseqid sstart send sstrand' -word_size 11 -gapopen 5 -gapextend 2 -penalty -3 -reward 2 -evalue {evalue} -perc_identity 0.1"
        run_cmd(f"{blast_cmd_base} -query '{b_fa}' -out '{blast_out_b}'", check=False)
        run_cmd(f"{blast_cmd_base} -query '{r_fa}' -out '{blast_out_r}'", check=False)
        cols = ['qseqid', 'qstart', 'qend', 'sseqid', 'sstart', 'send', 'sstrand']
        try:
            db, dr = pd.read_csv(blast_out_b, sep='\t', names=cols), pd.read_csv(blast_out_r, sep='\t', names=cols)
        except: return None
        if db.empty or dr.empty: return None
        db['ac'], dr['ac'] = db['qseqid'].apply(lambda x: x.split('-')[0]), dr['qseqid'].apply(lambda x: x.split('-')[0])
        valid_ac = set()
        for ac in set(db['ac']).intersection(set(dr['ac'])):
            b_hits, r_hits = db[db['ac'] == ac], dr[dr['ac'] == ac]
            try: b_pos, r_pos = int(ac.split('_')[0]), int(ac.split('_')[1].split('(')[0])
            except: continue
            for _, bh in b_hits.iterrows():
                for _, rh in r_hits.iterrows():
                    if dvg_type == 'DEL' and bh['sstrand'] != rh['sstrand']: continue
                    if dvg_type == 'CB' and bh['sstrand'] == rh['sstrand']: continue
                    r1s, r1e = sorted([bh['sstart'], bh['send']])
                    r2s, r2e = sorted([rh['sstart'], rh['send']])
                    if dvg_type == 'DEL':
                        if (r1e < r2s and (b_pos - n_shift <= r1e <= b_pos + n_shift) and (r_pos - n_shift <= r2s <= r_pos + n_shift)) or \
                           (r2e < r1s and (b_pos - n_shift <= r2e <= b_pos + n_shift) and (r_pos - n_shift <= r1s <= r_pos + n_shift)): valid_ac.add(ac)
                    elif dvg_type == 'CB':
                        if (r1s < r2s and (b_pos - n_shift <= r1s <= b_pos + n_shift) and (r_pos - n_shift <= r2s <= r_pos + n_shift)) or \
                           (r2s < r1s and (b_pos - n_shift <= r2s <= b_pos + n_shift) and (r_pos - n_shift <= r1s <= r_pos + n_shift)): valid_ac.add(ac)
        for f in [b_fa, r_fa, blast_out_b, blast_out_r]: Path(f).unlink(missing_ok=True)
        return df[df['A_C(A_C)'].isin(valid_ac)].copy()

    @staticmethod
    def generate_report_and_plot(df: pd.DataFrame, genome_len: int, dvg_type: str, out_prefix: str, n_shift: int = 5):
        if df is None or df.empty: return
        df['BREAK_POSITION'] = df['A_C(A_C)'].apply(lambda x: int(x.split('_')[0]))
        df['REJOIN_POSITION'] = df['A_C(A_C)'].apply(lambda x: int(x.split('_')[1].split('(')[0]))
        if dvg_type == 'DEL':
            df['LENGTH'] = genome_len - df['REJOIN_POSITION'] + df['BREAK_POSITION'] + 1
            df['DELETION_SIZE'] = df['REJOIN_POSITION'] - df['BREAK_POSITION'] - 1
            df['PERC_STANDARD'] = (df['LENGTH'] / genome_len) * 100
        else:
            df['LENGTH'] = (genome_len - df['BREAK_POSITION'] + 1) + (genome_len - df['REJOIN_POSITION'] + 1)
            df['LOOP_SIZE'] = df['REJOIN_POSITION'] - df['BREAK_POSITION']
            df['STEM_SIZE'] = (genome_len - df['REJOIN_POSITION'] + 1) * 2
            df['PERC_STEM'] = (df['STEM_SIZE'] / df['LENGTH']) * 100
        df = df.sort_values(by=['LENGTH', 'BREAK_POSITION'])
        species_list, b_ref, r_ref, s_ref = [], 0, 0, 0
        for _, row in df.iterrows():
            b, r, s = row['BREAK_POSITION'], row['REJOIN_POSITION'], row['LENGTH']
            if s == s_ref and b <= b_ref + n_shift: species_list.append(f"{b_ref}_{r_ref}")
            else: b_ref, r_ref, s_ref = b, r, s; species_list.append(f"{b_ref}_{r_ref}")
        df['SPECIES'] = species_list
        df.to_csv(f"{out_prefix}.all-info_{dvg_type}.N{n_shift}.txt", sep='\t', index=False)
        
        mode_records = []
        for species, group in df.groupby('SPECIES'):
            mode_b = group['BREAK_POSITION'].mode()[0]
            mode_r = group[group['BREAK_POSITION'] == mode_b]['REJOIN_POSITION'].iloc[0]
            group_copy = group.copy()
            group_copy['mode'] = f"{mode_b}_{mode_r}"
            mode_records.append(group_copy)
        mode_df = pd.concat(mode_records)
        mode_df.to_csv(f"{out_prefix}.all-info_{dvg_type}.N{n_shift}_mode.txt", sep='\t', index=False)
        
        plot_data = mode_df['mode'].value_counts().reset_index()
        plot_data.columns = ['mode', 'N']
        plot_data['BREAK_POSITION'] = plot_data['mode'].apply(lambda x: int(x.split('_')[0]))
        plot_data['REJOIN_POSITION'] = plot_data['mode'].apply(lambda x: int(x.split('_')[1]))
        
        plt.figure(figsize=(6, 5))
        sns.set_theme(style="whitegrid")
        main_pts, gray_pts = plot_data[plot_data['N'] >= 2], plot_data[plot_data['N'] < 2]
        if not gray_pts.empty: sns.scatterplot(data=gray_pts, x='BREAK_POSITION', y='REJOIN_POSITION', size='N', color='lightgray', alpha=0.5, legend=False)
        if not main_pts.empty: sns.scatterplot(data=main_pts, x='BREAK_POSITION', y='REJOIN_POSITION', size='N', color='#b6429f' if dvg_type == 'CB' else '#4293b6', alpha=0.7, sizes=(20, 200))
        plt.plot([0, genome_len], [0, genome_len], 'k--', alpha=0.2)
        plt.xlim(0, genome_len); plt.ylim(0, genome_len)
        plt.xlabel('Break position', fontweight='bold'); plt.ylabel('Rejoin position', fontweight='bold')
        plt.title(f"{Path(out_prefix).name} {dvg_type}", fontweight='bold')
        plt.legend(title='Junction reads', bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.tight_layout()
        plt.savefig(f"{out_prefix}.all-info_{dvg_type}.N{n_shift}_mode_plot.png", dpi=300, bbox_inches='tight')
        plt.close()
# ==========================================
# 3. 多进程 Worker 函数 (带完善日志流)
# ==========================================
def worker_align(args):
    sname, fq, index_prefix, out_bam_str, threads, resume, log_file, master_log = args
    out_bam = Path(out_bam_str)
    t_io = min(4, threads) 
    try:
        if not (out_bam.exists() and resume):
            fq_arg = f"-1 '{fq['r1']}' -2 '{fq['r2']}'" if fq["r2"] else f"-U '{fq['r1']}'"
            format_flag = "-f" if fq["is_fasta"] else "-q"
            cmd = f"bowtie2 --local -p {threads} {format_flag} -x '{index_prefix}' {fq_arg} | samtools sort -@ {t_io} -o '{out_bam}'"
            run_cmd(cmd, log_path=log_file, master_log=master_log, check=False)
        
        bai_file = Path(str(out_bam) + ".bai")
        if out_bam.exists() and (not bai_file.exists() or bai_file.stat().st_mtime < out_bam.stat().st_mtime):
            run_cmd(f"samtools index -@ {t_io} '{out_bam}'", log_path=log_file, master_log=master_log, check=False)
            
        return sname, str(out_bam), fq["is_fasta"]
    except Exception as e:
        err = f"\n[Python Exception] worker_align: {traceback.format_exc()}\n"
        with open(log_file, "a") as lf: lf.write(err)
        with open(master_log, "a") as mlf: mlf.write(err)
        return sname, None, fq["is_fasta"]

def worker_calc_metrics(args):
    s, v, row_dict, bam_path = args
    stats = {'read_count': 0, 'mean_coverage': 0.0, 'covered_bases': 0, 'Sites_0X': 0, 'Sites_LowCov': 0, 'Pi': 0.0, 'Shannon': 0.0, 'ANI': 0.0, 'contig_length': 0, 'bam_total_reads': 0}
    
    if bam_path and Path(bam_path).exists():
        try:
            with pysam.AlignmentFile(bam_path, "rb") as bamfile:
                try: stats['bam_total_reads'] = bamfile.mapped + bamfile.unmapped
                except: pass
                
                if v in bamfile.references:
                    length = bamfile.get_reference_length(v)
                    if length > 0:
                        cov = np.zeros(length, dtype=np.uint32)
                        bc = [Counter() for _ in range(length)]
                        idents = []
                        reads = [r for r in bamfile.fetch(v) if not r.is_secondary and not r.is_unmapped]
                        
                        if reads:
                            for r in reads:
                                for blk in r.get_blocks(): cov[blk[0]:blk[1]] += 1
                                for q, ref, base in r.get_aligned_pairs(matches_only=True, with_seq=True):
                                    if q is not None and ref is not None and ref < length:
                                        b = r.query_sequence[q].upper()
                                        if b in 'ACGTN': bc[ref][b] += 1
                                            
                                align_len = r.query_alignment_length
                                nm = r.get_tag('NM') if r.has_tag('NM') else None
                                if nm is not None and align_len > 0: 
                                    idents.append((align_len - nm) / align_len)
                            
                            pi_list = []
                            shannon_list = []
                            for pos, c in enumerate(bc):
                                if cov[pos] > 0 and c and sum(c.values()) > 0:
                                    tot = sum(c.values())
                                    pi_list.append((tot - max(c.values())) / tot)
                                    h = sum(-(count/tot)*math.log(count/tot) for count in c.values() if count > 0)
                                    shannon_list.append(h)
                                       
                            stats.update({
                                'read_count': len(reads), 
                                'mean_coverage': float(np.mean(cov)), 
                                'covered_bases': int(np.count_nonzero(cov)), 
                                'Sites_0X': int(np.count_nonzero(cov == 0)),               
                                'Sites_LowCov': int(np.count_nonzero((cov > 0) & (cov < 10))), 
                                'Pi': float(np.mean(pi_list)) if pi_list else 0.0, 
                                'Shannon': float(np.mean(shannon_list)) if shannon_list else 0.0, 
                                'ANI': float(np.mean(idents)) if idents else 0.0, 
                                'contig_length': length
                            })
        except Exception: pass
            
    row_dict.update(stats)
    return row_dict
import subprocess
import tempfile
import traceback
from pathlib import Path

import subprocess # 确保导入了这个模块
import shlex
import tempfile
from pathlib import Path

import subprocess
import tempfile
import traceback
from pathlib import Path

def worker_extract_reads(args):
    sample, virus, bam_path, is_fasta, r1_str, r2_str, rs_str, threads, resume, log_file, master_log = args

    # 1. 彻底清理变量中的换行符和空格
    virus = str(virus).strip()
    bam_path = str(bam_path).strip()
    r1, r2, rs = Path(r1_str), Path(r2_str), Path(rs_str)

    # 2. 只有当文件存在 且 大小大于 50 bytes 时，才真正跳过（防止被空文件欺骗）
    if resume and r1.exists() and r1.stat().st_size > 50:
        return True

    t_io = min(4, threads)
    ext_tool = "samtools fasta" if is_fasta else "samtools fastq"

    with tempfile.TemporaryDirectory(prefix="extract_") as tmpdir:
        tmp_names = Path(tmpdir) / "reads.names"
        tmp_r1 = Path(tmpdir) / "R1.tmp"
        tmp_r2 = Path(tmpdir) / "R2.tmp"
        tmp_rs = Path(tmpdir) / "RS.tmp"
        tmp_0 = Path(tmpdir) / "0.tmp"

        try:
            # === 步骤 1: 提取名单 ===
            cmd_get_names = f"samtools view -@ {t_io} '{bam_path}' '{virus}' | cut -f1 | awk '!seen[$0]++' > '{tmp_names}'"
            res_names = subprocess.run(cmd_get_names, shell=True, executable='/bin/bash', capture_output=True, text=True)

            if res_names.returncode != 0:
                print(f"\n❌ [名单提取报错] {sample}: {res_names.stderr}")
                return False

            if not tmp_names.exists() or tmp_names.stat().st_size == 0:
                return True # 确实没有reads，正常退出

            # === 步骤 2: 提取序列 ===
            cmd_ext = (
                f"samtools view -@ {t_io} -h -N '{tmp_names}' '{bam_path}' | "
                f"samtools collate -O -u -@ {t_io} - | "
                f"{ext_tool} -@ {t_io} -1 '{tmp_r1}' -2 '{tmp_r2}' -s '{tmp_rs}' -0 '{tmp_0}' -n -"
            )
            res_ext = subprocess.run(cmd_ext, shell=True, executable='/bin/bash', capture_output=True, text=True)

            if res_ext.returncode != 0:
                print(f"\n❌ [序列提取报错] {sample}: {res_ext.stderr}")
                return False

            # 合并单端
            if tmp_0.exists() and tmp_0.stat().st_size > 0:
                subprocess.run(f"cat '{tmp_0}' >> '{tmp_rs}'", shell=True, executable='/bin/bash')

            # === 步骤 3: 压缩输出 ===
            if tmp_r1.exists() and tmp_r1.stat().st_size > 0:
                subprocess.run(f"pigz -p {t_io} -c '{tmp_r1}' > '{r1}'", shell=True, executable='/bin/bash', check=True)
            if tmp_r2.exists() and tmp_r2.stat().st_size > 0:
                subprocess.run(f"pigz -p {t_io} -c '{tmp_r2}' > '{r2}'", shell=True, executable='/bin/bash', check=True)
            if tmp_rs.exists() and tmp_rs.stat().st_size > 0:
                subprocess.run(f"pigz -p {t_io} -c '{tmp_rs}' > '{rs}'", shell=True, executable='/bin/bash', check=True)

            return True

        except Exception as e:
            print(f"\n❌ 样本 {sample} 发生 Python 异常:\n{traceback.format_exc()}")
            return False

def worker_consensus(args):
    sample, virus, bam_path, ref_fa, out_fa_str, fixed_bam_str, depth, qual, freq, ambig, threads, resume, log_file, master_log = args
    out_fa, fixed_bam = Path(out_fa_str), Path(fixed_bam_str)
    if resume and out_fa.exists() and out_fa.stat().st_size > 0: return True
    t_io = min(4, threads)
    awk_filter = f"awk -v v='{virus}' '/^@SQ/ && $2 != \"SN:\"v {{next}} {{print}}'"
    extract_cmd = f"samtools view -@ {t_io} -h '{bam_path}' '{virus}' | {awk_filter} | samtools view -@ {t_io} -b | samtools sort -@ {t_io} -o '{fixed_bam}'"
    vc_cmd = f"viral_consensus -i '{fixed_bam}' -r '{ref_fa}' -o '{out_fa}' -q {qual} -d {depth} -f {freq} -a {ambig}"
    try:
        run_cmd(extract_cmd, log_path=log_file, master_log=master_log)
        run_cmd(f"samtools index -@ {t_io} '{fixed_bam}'", log_path=log_file, master_log=master_log)
        run_cmd(vc_cmd, log_path=log_file, master_log=master_log)
        return True
    except Exception as e:
        err = f"\n[Python Exception] worker_consensus: {traceback.format_exc()}\n"
        with open(log_file, "a") as lf: lf.write(err)
        with open(master_log, "a") as mlf: mlf.write(err)
        return False
    finally:
        for f in (fixed_bam, Path(str(fixed_bam)+".bai")):
            if f.exists(): f.unlink()

def parse_ann_to_tsv(ann_vcf: str, out_tsv: str):
    opener = gzip.open if ann_vcf.endswith('.gz') else open
    try:
        with opener(ann_vcf, 'rt') as f_in, open(out_tsv, 'w') as f_out:
            f_out.write("CHROM\tPOS\tREF\tALT\tGENE\tEFFECT\tIMPACT\tDNA_CHANGE\tAA_CHANGE\n")
            for line in f_in:
                if line.startswith("#"): continue
                cols = line.strip().split("\t")
                if len(cols) < 8: continue
                ann_field = [x for x in cols[7].split(";") if x.startswith("ANN=")]
                if not ann_field: continue
                first_ann = ann_field[0][4:].split(",")[0].split("|")
                if len(first_ann) > 10: f_out.write("\t".join([cols[0], cols[1], cols[3], cols[4], first_ann[3], first_ann[1], first_ann[2], first_ann[9], first_ann[10]]) + "\n")
    except Exception as e: logger.error(f"解析 TSV 失败: {e}")

def worker_variants(args):
    (sample, virus, mean_depth, bam_path, ref_fa, raw_out_str, clean_vcf_str, fixed_bam_str, caller, db_ready, snpeff_jar, snpeff_config, snpeff_db_name, snpeff_mem, disable_dyn, se_dir_str, run_snpgenie, gtf_dir_str, sg_out_dir, threads, resume, log_file, master_log) = args
    raw_out, clean_vcf, fixed_bam = Path(raw_out_str), Path(clean_vcf_str), Path(fixed_bam_str)
    t_io = min(4, threads)
    
    if not (resume and clean_vcf.exists() and clean_vcf.stat().st_size > 0):
        awk_filter = f"awk -v v='{virus}' '/^@SQ/ && $2 != \"SN:\"v {{next}} {{print}}'"
        extract_cmd = f"samtools view -@ {t_io} -h '{bam_path}' '{virus}' | {awk_filter} | samtools view -@ {t_io} -b | samtools sort -@ {t_io} -o '{fixed_bam}'"
        try:
            run_cmd(extract_cmd, log_path=log_file, master_log=master_log)
            run_cmd(f"samtools index -@ {t_io} '{fixed_bam}'", log_path=log_file, master_log=master_log)
            if caller == "freebayes": 
                run_cmd(f"freebayes -p 1 -f '{ref_fa}' '{fixed_bam}' > '{raw_out}'", log_path=log_file, master_log=master_log)
            if caller != "ivar":
                dp = 100 if disable_dyn else (5 if mean_depth < 50 else (15 if mean_depth < 1000 else 100))
                frq = 0.5
                saf = 1 if dp == 5 else (2 if dp == 15 else 10)
                flt = f"QUAL>20 && INFO/DP>={dp} && INFO/SAF>={saf} && INFO/SAR>={saf} && (INFO/AO/INFO/DP)>{frq}"
                soft = Path(str(clean_vcf).replace(".filtered.", ".soft."))
                run_cmd(f"bcftools filter --threads {t_io} -s FAIL -i '{flt}' -Ov -o '{soft}' '{raw_out}'", log_path=log_file, master_log=master_log)
                run_cmd(f"bcftools filter --threads {t_io} -i 'FILTER==\"PASS\"' -Ov -o '{clean_vcf}' '{soft}'", log_path=log_file, master_log=master_log)
                if soft.exists(): soft.unlink()
        except Exception as e:
            err = f"\n[Python Exception] worker_variants extraction: {traceback.format_exc()}\n"
            with open(log_file, "a") as lf: lf.write(err)
            with open(master_log, "a") as mlf: mlf.write(err)
            return False
        finally:
            for tmp in (fixed_bam, Path(str(fixed_bam)+".bai")):
                if tmp.exists(): tmp.unlink()

    if clean_vcf.exists() and clean_vcf.stat().st_size > 0:
        
        # 🌟 1. 提取变异等位基因频率表 (AF)
        af_tsv = Path(raw_out_str).parent / f"{sample}.{virus}.allele_frequencies.tsv"
        extract_allele_frequency(str(clean_vcf), str(af_tsv))
        
        # 🌟 2. 运行 SnpEff 注释
        if db_ready:
            se_dir = Path(se_dir_str)
            ann_vcf, sum_tsv = se_dir / f"{clean_vcf.name.replace('.filtered.vcf', '.ann.vcf')}", se_dir / f"{clean_vcf.name.replace('.filtered.vcf', '.annotation_summary.tsv')}"
            # 注意：把 > '{ann_vcf}' 留在命令里，因为这是 SnpEff 正常的输出重定向，同时 run_cmd 会完美捕获它的 STDERR 运行日志
            cmd_ann = f"java -Xmx{snpeff_mem} -jar '{snpeff_jar}' ann -c '{snpeff_config}' -noStats {snpeff_db_name} '{clean_vcf}' > '{ann_vcf}'"
            try: 
                run_cmd(cmd_ann, log_path=log_file, master_log=master_log)
                parse_ann_to_tsv(str(ann_vcf), str(sum_tsv))
            except Exception as e:
                with open(master_log, "a") as mlf: mlf.write(f"\n[Python Error] SnpEff failed: {e}\n")
            
        # 🌟 3. 运行 SNPGenie
        gtf_file = Path(gtf_dir_str) / f"{virus}.gtf"
        if run_snpgenie:
            if not gtf_file.exists() or gtf_file.stat().st_size == 0:
                generate_dummy_gtf(str(ref_fa), str(gtf_file))
                
            sg_dir = Path(sg_out_dir) / f"{clean_vcf.name.replace('.filtered.vcf', '')}"
            sg_dir.mkdir(parents=True, exist_ok=True)
            try:
                vcf_work, fa_work, gtf_work = sg_dir/f"{sample}.vcf", sg_dir/f"{virus}.fasta", sg_dir/f"{virus}.gtf"
                shutil.copy(clean_vcf, vcf_work); shutil.copy(ref_fa, fa_work); shutil.copy(gtf_file, gtf_work)
                
                cmd_sg = f"cd '{sg_dir}' && snpgenie.pl --vcfformat=2 --snpreport='{vcf_work.name}' --fastafile='{fa_work.name}' --gtffile='{gtf_work.name}'"
                # 不用check=False，利用 run_cmd 全面捕获输出
                run_cmd(cmd_sg, log_path=log_file, master_log=master_log, check=False)
                
                for f in [vcf_work, fa_work, gtf_work]: f.unlink(missing_ok=True)
            except Exception as e:
                err = f"\n[Python Exception] SNPGenie failed: {traceback.format_exc()}\n"
                with open(log_file, "a") as lf: lf.write(err)
                with open(master_log, "a") as mlf: mlf.write(err)
                
    return True

def worker_build_vodka2_db(args):
    virus, ref_fa, bp, rl, db_dir_str, core_pl, threads, log_file, master_log = args
    db_dir = Path(db_dir_str)
    
    cb_fa = db_dir / f"{virus}.CB.{bp}.{rl}.fasta"
    del_fa = db_dir / f"{virus}.DEL.{bp}.{rl}.fasta"
    cb_idx = db_dir / f"{virus}.CB.{bp}.{rl}"
    del_idx = db_dir / f"{virus}.DEL.{bp}.{rl}"
    
    try:
        if not Path(f"{cb_idx}.1.bt2").exists() and not Path(f"{cb_idx}.1.bt2l").exists():
            run_cmd(f"perl '{core_pl}' build --type cb --fasta '{ref_fa}' --bases {bp} --size {rl} --out '{cb_fa}'", log_path=log_file, master_log=master_log)
            run_cmd(f"bowtie2-build --threads {threads} '{cb_fa}' '{cb_idx}'", log_path=log_file, master_log=master_log)
            
        if not Path(f"{del_idx}.1.bt2").exists() and not Path(f"{del_idx}.1.bt2l").exists():
            run_cmd(f"perl '{core_pl}' build --type del --fasta '{ref_fa}' --bases {bp} --size {rl} --out '{del_fa}'", log_path=log_file, master_log=master_log)
            run_cmd(f"bowtie2-build --threads {threads} '{del_fa}' '{del_idx}'", log_path=log_file, master_log=master_log)
            
        return virus, str(cb_idx), str(del_idx), bp
    except Exception as e:
        err = f"\n[Python Exception] worker_build_vodka2_db: {traceback.format_exc()}\n"
        with open(log_file, "a") as lf: lf.write(err)
        with open(master_log, "a") as mlf: mlf.write(err)
        return virus, None, None, bp

def worker_vodka2(args):
    # 🌟 第一层防弹衣：拦截参数解包错误
    try:
        sample, virus, bam_path, is_fasta, ref_fa, cb_idx, del_idx, readlen, bp, core_pl, target_dir_str, threads, log_file, master_log, resume = args
    except Exception as e:
        # 如果解包失败，尝试提取 master_log 路径并写入报错
        fallback_log = args[13] if len(args) > 13 else "VODKA2_CRITICAL_ERROR.log"
        with open(fallback_log, "a") as f:
            f.write(f"\n[CRITICAL ERROR] worker_vodka2 参数解包失败! 期望 15 个参数，实际收到 {len(args)} 个。\n错误信息: {e}\n")
        raise RuntimeError(f"参数解包失败: {e}")

    # 🌟 第二层防弹衣：拦截所有业务逻辑错误
    try:
        target_dir = Path(target_dir_str)

        cb_plot = target_dir / "CB_OUT" / f"{sample}_CB.all-info_CB.N5_mode_plot.png"
        del_plot = target_dir / "DEL_OUT" / f"{sample}_DEL.all-info_DEL.N5_mode_plot.png"
        if resume and cb_plot.exists() and del_plot.exists():
            return True

        genome_len = 0
        with open(ref_fa, 'r') as f:
            for line in f:
                if not line.startswith('>'): genome_len += len(line.strip())

        ext_tool = "samtools fasta" if is_fasta else "samtools fastq"
        format_flag = "-f" if is_fasta else "-q"
        unmapped_reads = target_dir / (f"{sample}_unmapped.fa" if is_fasta else f"{sample}_unmapped.fq")

        cmd_extract = f"samtools view -@ {threads} -b -f 4 '{bam_path}' | {ext_tool} -@ {threads} - > '{unmapped_reads}'"

        run_cmd(cmd_extract, log_path=log_file, master_log=master_log)
        if not unmapped_reads.exists() or unmapped_reads.stat().st_size < 100: return True

        if not Path(f"{ref_fa}.nhr").exists() and not Path(f"{ref_fa}.nsq").exists():
            run_cmd(f"makeblastdb -in '{ref_fa}' -dbtype nucl", log_path=log_file, master_log=master_log, check=False)

        for dvg_type, idx in [('CB', cb_idx), ('DEL', del_idx)]:
            out_dir = target_dir / f"{dvg_type}_OUT"
            out_dir.mkdir(exist_ok=True)
            idx_name = Path(idx).name

            cmd_run = (
                f"bowtie2 -p {threads} -x '{idx}' {format_flag} --very-fast-local -U '{unmapped_reads}' --no-sq --no-unal --mp 0,0 | "
                f"perl '{core_pl}' parse --sam - --readlen {readlen} --bp {bp} --sample {sample}_{dvg_type} --outdir '{out_dir}' --index '{idx_name}'"
            )
            run_cmd(cmd_run, log_path=log_file, master_log=master_log, check=False)

            results_txt = out_dir / "results" / f"{sample}_{dvg_type}_{idx_name}_RESULTS.txt"
            if results_txt.exists():
                confirmed_df = VODKA2_PostProcessor.run_blast_validation(results_txt, Path(ref_fa), dvg_type)
                if confirmed_df is not None and not confirmed_df.empty:
                    out_prefix = str(out_dir / f"{sample}_{dvg_type}")
                    VODKA2_PostProcessor.generate_report_and_plot(confirmed_df, genome_len, dvg_type, out_prefix)

    except Exception as e:
        err = f"\n[Python Exception] worker_vodka2: {traceback.format_exc()}\n"
        with open(log_file, "a") as lf: lf.write(err)
        with open(master_log, "a") as mlf: mlf.write(err)
        raise RuntimeError(f"VODKA2 业务逻辑崩溃: {e}")
    finally:
        if 'unmapped_reads' in locals() and unmapped_reads.exists(): unmapped_reads.unlink()

    return True
    
# ==========================================
# 4. 核心流水线 
# ==========================================
class PostProcessPipeline:
    def __init__(self, args):
        self.args = args
        self.out = Path(args.output_dir)

        self.d_fasta       = self.out / "virus-fasta"
        self.d_index       = self.out / "virus-index"
        self.d_bam         = self.out / "virus-bam"
        self.d_consensus   = self.out / "virus-consensus"
        self.d_variants    = self.out / "virus-variants"
        self.d_snpeff      = self.out / "virus-SnpEff"
        self.d_snpgenie    = self.out / "virus-SNPGenie" 
        self.d_reads       = self.out / "reads"
        self.d_summary     = self.out / "summary"
        self.d_ind_reports = self.out / "summary" / "individual_virus_reports"
        self.d_logs        = self.out / "logs"

        self.index_prefix: str = ""
        self.harmonized_fasta = self.d_index / "harmonized_references.fasta"
        self.snpeff_db_name = f"vap_db_{datetime.now():%Y%m%d_%H%M%S}"
        
        # 🚀 新增：全流程主日志文件路径
        self.master_log = str(self.d_logs / "master_commands.log")

        self._validate_args(); self._check_tools(); self._setup_dirs(); self._add_file_logger()

    def _validate_args(self):
        a = self.args
        for f in (a.summary, a.reference, a.info):
            if not Path(f).exists(): logger.error(f"❌ 文件不存在: {f}"); sys.exit(1)
        if a.fastq and not Path(a.fastq).is_dir(): logger.error(f"❌ --fastq 目录不存在: {a.fastq}"); sys.exit(1)

    def _check_tools(self):
        need = {"samtools": "BAM 操作"}
        if self.args.fastq: need.update({"bowtie2": "序列比对", "bowtie2-build": "索引构建"})
        if self.args.extract_reads: need["pigz"] = "多线程压缩"
        if self.args.consensus: need["viral_consensus"] = "共识构建"
        if self.args.call_variants:
            need[self.args.variant_caller] = "变异检测"
            if self.args.variant_caller in ("freebayes", "lofreq"): need["bcftools"] = "VCF 过滤"
        if self.args.snpeff: need["java"] = "SnpEff 注释引擎"
        if self.args.snpgenie: need["snpgenie.pl"] = "SNPGenie"
        missing = [f"{t} ({d})" for t, d in need.items() if not shutil.which(t)]
        if missing: logger.error("❌ 缺少必要工具:\n  " + "\n  ".join(f"- {m}" for m in missing)); sys.exit(1)

    def _setup_dirs(self):
        dirs = [self.d_fasta, self.d_summary, self.d_ind_reports, self.d_logs, self.d_bam]
        if self.args.fastq: dirs.append(self.d_index)
        if self.args.extract_reads: dirs.append(self.d_reads)
        if self.args.consensus: dirs.append(self.d_consensus)
        if self.args.call_variants: 
            dirs.append(self.d_variants)
            if self.args.snpeff: dirs.append(self.d_snpeff)
            if self.args.snpgenie: dirs.append(self.d_snpgenie)
        for d in dirs: d.mkdir(parents=True, exist_ok=True)

    def _add_file_logger(self):
        log_path = self.d_logs / f"pipeline_{datetime.now():%Y%m%d_%H%M%S}.log"
        if not any(isinstance(h, logging.FileHandler) for h in logger.handlers):
            fh = logging.FileHandler(log_path, encoding="utf-8")
            fh.setFormatter(logging.Formatter("%(asctime)s - %(levelname)s - %(message)s", datefmt="%Y-%m-%d %H:%M:%S"))
            logger.addHandler(fh)

    def load_data_and_harmonize(self) -> pd.DataFrame:
        sep = "," if Path(self.args.summary).suffix.lower() == ".csv" else "\t"
        df = pd.read_csv(self.args.summary, sep=sep, dtype=str)
        df.columns = df.columns.str.strip()
        df['Sample'] = df['Sample'].apply(lambda x: re.sub(r'(?i)(_clean|_trimmed|_filtered|_val|_fastp)+$', '', str(x).strip()))
        
        col_mappings = {'Virus': ['Rep_Accession', 'Accession', 'Virus GENBANK accession'], 'taxonomy': ['Adjusted_Species', 'Species_NCBI', 'Species_ICTV', 'Taxonomy'], 'MeanDepth': ['Rep_MeanDepth', 'MeanDepth', 'Depth'], 'Taxid': ['taxid', 'Taxid']}
        for target, fallbacks in col_mappings.items():
            if target not in df.columns:
                for col in fallbacks:
                    if col in df.columns: df[target] = df[col]; break

        if "taxonomy" not in df.columns: df["taxonomy"] = "Unannotated"
        if "Virus" not in df.columns: logger.error("❌ 无法在表中找到病毒 Accession 列！"); sys.exit(1)

        info_df = pd.read_csv(self.args.info, sep='\t', dtype=str)
        info_cols = [c for c in info_df.columns if c in ['Accession', 'Segment']]
        info_df = info_df[info_cols].drop_duplicates(subset=['Accession'])
        
        df = df.merge(info_df, left_on='Virus', right_on='Accession', how='left', suffixes=('', '_info'))
        if "Segment_info" in df.columns: df["Segment"] = df.get("Segment", pd.Series(dtype=str)).fillna(df["Segment_info"])
        df['Segment'] = df.get('Segment', pd.Series(dtype=str)).fillna('Unsegmented')

        logger.info(f"📋 读取汇总表: 识别出 {len(df)} 记录 | 包含 {df['Sample'].nunique()} 个独立样本")
        logger.info("\n🧬 开始跨样本参考基因组强制统一 (基于频率投票):")
        
        valid_tax = df[~df['taxonomy'].isin(['Unannotated', '-', ''])]
        tax_to_best_virus = {}
        for (tax, seg), group in valid_tax.groupby(['taxonomy', 'Segment']): tax_to_best_virus[(tax, seg)] = group['Virus'].value_counts().idxmax()
        df['Original_Virus'] = df['Virus']
        df['Virus'] = df.apply(lambda row: tax_to_best_virus.get((row['taxonomy'], row['Segment']), row['Virus']), axis=1)
        
        self.raw_input_df = df.copy() 
        return df.drop_duplicates(subset=['Sample', 'taxonomy', 'Segment', 'Virus'])

    def extract_virus_fastas(self, df: pd.DataFrame) -> dict:
        target_set = set(df["Virus"].unique().tolist())
        found_map, seq_buf, vid_cur = {}, [], None

        def _flush():
            if vid_cur not in target_set or not seq_buf: return
            folder = f"ref_{safe_name(vid_cur)}"
            vdir = self.d_fasta / folder; vdir.mkdir(parents=True, exist_ok=True)
            ref_fa = vdir / f"{folder}.ref.fasta"
            if not ref_fa.exists():
                with open(ref_fa, "w") as f: f.write(f">{vid_cur}\n" + "".join(seq_buf) + "\n")
            if not Path(str(ref_fa) + ".fai").exists():
                run_cmd(f"samtools faidx '{ref_fa}'", log_path=str(self.d_logs/"samtools_faidx.log"), master_log=self.master_log)
            found_map[vid_cur] = ref_fa

        logger.info(f"\n🔍 从库中提取 {len(target_set)} 条代表株序列 → virus-fasta/")
        with open(self.args.reference, "r") as fh:
            for line in fh:
                line = line.rstrip()
                if line.startswith(">"): _flush(); vid_cur = line[1:].split()[0]; seq_buf = []
                else: seq_buf.append(line)
        _flush()
        return found_map

    def resolve_bam_map(self, samples: list, ref_map: dict) -> dict:
        if self.args.fastq:
            self.index_prefix = str(self.d_index / "harmonized_bowtie2")
            if not Path(f"{self.index_prefix}.1.bt2").exists():
                with open(self.harmonized_fasta, "w") as out_f:
                    for fp in ref_map.values():
                        with open(fp, "r") as in_f: out_f.write(in_f.read() + "\n")
                run_cmd(f"bowtie2-build --threads {self.args.threads} '{self.harmonized_fasta}' '{self.index_prefix}'", log_path=str(self.d_logs/"bowtie2_build.log"), master_log=self.master_log)

            fq_dir = Path(self.args.fastq)
            logger.info(f"\n📁 正在扫描 {fq_dir} (包含所有子文件夹)...")
            all_files = [f for f in fq_dir.rglob("*") if f.is_file() and any(ext in f.name.lower() for ext in ['.fq', '.fastq', '.fa', '.fasta'])]
            sample_files, unmatched = {}, []
            
            for sname in samples:
                s_clean = sname.strip().lower()
                matched = []
                for f in all_files:
                    f_spaced = f.name.lower().replace('_', ' ').replace('.', ' ').replace('-', ' ')
                    s_spaced = s_clean.replace('_', ' ').replace('.', ' ').replace('-', ' ')
                    if re.search(r'\b' + re.escape(s_spaced) + r'\b', f_spaced) or f.name.lower().startswith(s_clean + "_") or f.name.lower().startswith(s_clean + "."): matched.append(f)
                matched = list(set(matched))
                if matched:
                    r1, r2, is_fasta = None, None, False
                    for f in matched:
                        nl = f.name.lower()
                        if any(x in nl for x in ['.fa', '.fasta']): is_fasta = True
                        if any(x in nl for x in ['_r2', '_2.', '.r2', '_2_']): r2 = f
                        elif any(x in nl for x in ['_r1', '_1.', '.r1', '_1_']): r1 = f
                        else:
                            if not r1: r1 = f 
                    if r1: sample_files[sname] = {"r1": str(r1), "r2": str(r2) if r2 else None, "is_fasta": is_fasta}
                else: unmatched.append(sname)
            
            if unmatched:
                logger.warning(f"⚠️ 警告: 有 {len(unmatched)} 个样本未能找到对应的原始序列文件！")
            
            targets = [s for s in samples if s in sample_files]
            logger.info(f"\n🚀 开始精准靶向重比对 - 成功挂载 {len(targets)} 个样本的数据包 (Jobs:{self.args.jobs}, Threads/Job:{self.args.threads}) ...")
            tasks = [(s, sample_files[s], self.index_prefix, str(self.d_bam / f"{s}.sorted.bam"), self.args.threads, self.args.resume, str(self.d_logs / f"{s}.align.log"), self.master_log) for s in targets]

            bam_map = {}
            with ProcessPoolExecutor(max_workers=min(len(tasks), self.args.jobs)) as ex:
                for fut in tqdm(as_completed({ex.submit(worker_align, t): t for t in tasks}), total=len(tasks), desc="并行比对"):
                    s, bp, isfa = fut.result(); bam_map[s] = {"bam": bp, "is_fasta": isfa}
            return bam_map
        else:
            bam_map = {}
            for f in sorted(Path(self.args.bam).glob("*.bam")):
                stem = f.name.replace(".sorted.bam", "").replace(".sort.bam", "").replace(".bam", "")
                if stem not in bam_map or "sorted" in f.name:
                    bam_map[stem] = {"bam": str(f), "is_fasta": False}
                    bai = Path(str(f) + ".bai")
                    if not bai.exists() or bai.stat().st_mtime < f.stat().st_mtime:
                        run_cmd(f"samtools index -@ {self.args.threads} '{f}'", log_path=str(self.d_logs/"samtools_index.log"), master_log=self.master_log)
            return {s: bam_map[s] for s in samples if s in bam_map}

    def generate_final_summaries(self, df: pd.DataFrame, bam_map: dict) -> pd.DataFrame:
        logger.info(f"\n📊 并行读取 BAM 重新独立计算各丰度及物理指标...")
        tasks = [(r["Sample"], r["Virus"], r.to_dict(), bam_map[r["Sample"]]['bam'] if r["Sample"] in bam_map else None) for _, r in df.iterrows()]
        raw_rows = []
        with ProcessPoolExecutor(max_workers=min(len(tasks), self.args.jobs)) as ex:
            for fut in tqdm(as_completed([ex.submit(worker_calc_metrics, t) for t in tasks]), total=len(tasks), desc="基础指标重算"): raw_rows.append(fut.result())
        
        raw_df = pd.DataFrame(raw_rows)
        
        sample_lib_size = {}
        for col in ['Total_Reads', 'Clean_Reads', 'Library_Size', 'TotalReads']:
            if col in self.raw_input_df.columns:
                for _, r in self.raw_input_df.iterrows():
                    try: 
                        val = float(r[col])
                        if val > 0: sample_lib_size[r['Sample']] = val
                    except: pass
                if sample_lib_size: break
                
        if not sample_lib_size:
            for s in raw_df['Sample'].unique():
                bt = raw_df[raw_df['Sample'] == s]['bam_total_reads'].max()
                if pd.notna(bt) and bt > 0: sample_lib_size[s] = float(bt)

        agg_rows = []
        for (sample, taxonomy), group in raw_df.groupby(['Sample', 'taxonomy']):
            reads, total_len = group['read_count'].sum(), group['contig_length'].sum()
            agg_rows.append({
                'Sample': sample, 'Accession': ",".join(group['Virus'].unique()), 'Taxonomy': taxonomy, 'Taxid': ",".join(group['Taxid'].unique()), 'Segment': ",".join(group['Segment'].unique()),
                'Length': int(total_len), 'CoveredBases': int(group['covered_bases'].sum()), 'Reads': int(reads),
                'Sites_0X': int(group['Sites_0X'].sum()), 'Sites_LowCov': int(group['Sites_LowCov'].sum()), 
                'Avg_Read_ANI(%)': np.average(group['ANI'], weights=group['read_count']) * 100 if reads > 0 else 0.0,
                'Pi_avr': np.average(group['Pi'], weights=group['read_count']) if reads > 0 else 0.0,
                'Shannon_avr': np.average(group['Shannon'], weights=group['read_count']) if reads > 0 else 0.0
            })
        agg_df = pd.DataFrame(agg_rows)
        
        mapped_sums = agg_df.groupby('Sample')['Reads'].sum().to_dict()
        agg_df['Total_Reads'] = agg_df['Sample'].apply(lambda x: sample_lib_size.get(x, 0) or mapped_sums.get(x, 1))
        
        agg_df['RPK'] = agg_df.apply(lambda x: x['Reads'] / (x['Length'] / 1000.0) if x['Length'] > 0 else 0.0, axis=1)
        agg_df = agg_df.merge(agg_df.groupby('Sample')['RPK'].sum().reset_index().rename(columns={'RPK': 'Total_RPK'}), on='Sample')
        
        def _calc_metrics(r):
            tr = r['Total_Reads'] 
            trpk = r['Total_RPK'] 
            
            rpm = (r['Reads'] / tr * 1e6) if tr > 0 else 0.0
            rpkm = (r['RPK'] / (tr / 1e6)) if tr > 0 else 0.0
            tpm = (r['RPK'] / trpk * 1e6) if trpk > 0 else 0.0
            cov = (r['CoveredBases'] / r['Length'] * 100) if r['Length'] > 0 else 0.0
            
            return pd.Series({
                'Covered%': round(cov, 2), 'Sites_0X': r['Sites_0X'], 'Sites_LowCov(1-9X)': r['Sites_LowCov'],
                'CPM': round(rpm, 2), 'RPM': round(rpm, 2), 'FPKM': round(rpkm, 2), 'RPKM': round(rpkm, 2), 'TPM': round(tpm, 2), 
                'Avg_Read_ANI(%)': round(r['Avg_Read_ANI(%)'], 2), 'Pi_avr': round(r['Pi_avr'], 4), 'Shannon_avr': round(r['Shannon_avr'], 4)
            })

        final_df = pd.concat([agg_df[['Sample', 'Accession', 'Taxonomy', 'Taxid', 'Length', 'Reads', 'Segment']], agg_df.apply(_calc_metrics, axis=1)], axis=1)
        
        final_df = final_df[['Sample', 'Accession', 'Taxonomy', 'Taxid', 'Length', 'Covered%', 'Sites_0X', 'Sites_LowCov(1-9X)', 'Reads', 'CPM', 'RPM', 'FPKM', 'RPKM', 'TPM', 'Avg_Read_ANI(%)', 'Pi_avr', 'Shannon_avr', 'Segment']]

        out_file = self.d_summary / "all_summary.tsv"
        final_df.to_csv(out_file, sep='\t', index=False)
        for _, r in final_df.iterrows(): pd.DataFrame([r.to_dict()]).to_csv(self.d_ind_reports / f"{r['Sample']}_{safe_name(r['Taxonomy'])}.summary.tsv", sep='\t', index=False)

        logger.info("\n🕸️ 正在生成共感染丰度矩阵 (Coinfection Matrix)...")
        coinf_df = final_df.pivot_table(index='Sample', columns='Taxonomy', values='Reads', aggfunc='sum').fillna(0).astype(int)
        coinf_file = self.d_summary / "Coinfection_Matrix_Reads.tsv"
        coinf_df.to_csv(coinf_file, sep='\t')
        
        raw_df.rename(columns={'read_count': 'Recalc_Reads', 'mean_coverage': 'Recalc_MeanDepth'}, inplace=True)
        return raw_df

    def run_extract_reads(self, df: pd.DataFrame, bam_map: dict):
        tasks = []
        for _, r in df.iterrows():
            if r["Sample"] in bam_map:
                s, v, isf = r["Sample"], r["Virus"], bam_map[r["Sample"]]["is_fasta"]
                t = r.get("taxonomy", "Unannotated")
                
                # 🌟 生成极其干爽的层级目录名称 (替换掉异常字符与点号)
                L1 = f"{safe_name(t).replace('.', '_')}_{safe_name(v).replace('.', '_')}"
                L2 = f"{safe_name(s).replace('.', '_')}_{safe_name(v).replace('.', '_')}"
                
                # 建立双层嵌套文件夹: virus-reads / Taxonomy_Acc / Sample_Acc /
                target_dir = self.d_reads / L1 / L2
                target_dir.mkdir(parents=True, exist_ok=True)
                
                ext = 'fasta.gz' if isf else 'fastq.gz'
                r1 = str(target_dir / f"{L2}_R1.{ext}")
                r2 = str(target_dir / f"{L2}_R2.{ext}")
                rs = str(target_dir / f"{L2}_single.{ext}")
                
                # 💡 精妙之举：把本来堆在统一日志夹的 log 降维下放到这个样本自己的文件夹下
                log_file = str(target_dir / f"{L2}_extract_reads.log")
                
                tasks.append((s, v, bam_map[s]["bam"], isf, r1, r2, rs, self.args.threads, self.args.resume, log_file, self.master_log))
        
        logger.info(f"\n📦 靶向提取特异性病毒 Reads (Jobs:{self.args.jobs}, Threads/Job:{min(4, self.args.threads)})...")
        with ProcessPoolExecutor(max_workers=min(len(tasks), self.args.jobs)) as ex:
            for _ in tqdm(as_completed([ex.submit(worker_extract_reads, t) for t in tasks]), total=len(tasks)): pass
    
    def run_consensus(self, df: pd.DataFrame, bam_map: dict, ref_map: dict):
        tasks = []
        for _, r in df.iterrows():
            if r["Sample"] in bam_map and r["Virus"] in ref_map:
                s, v, t = r["Sample"], r["Virus"], r.get("taxonomy", "Unannotated")
                
                # 🌟 统一层级目录命名
                L1 = f"{safe_name(t).replace('.', '_')}_{safe_name(v).replace('.', '_')}"
                L2 = f"{safe_name(s).replace('.', '_')}_{safe_name(v).replace('.', '_')}"
                
                # virus-consensus / Taxonomy_Acc / Sample_Acc /
                vdir = self.d_consensus / L1 / L2
                vdir.mkdir(parents=True, exist_ok=True)
                
                out_fa = str(vdir / f"{L2}.consensus.fasta")
                fixed_bam = str(vdir / f"{L2}.fixed.bam")
                log_file = str(vdir / f"{L2}_consensus.log")
                
                d = max(1, min(10, int(math.floor(float(r.get("Recalc_MeanDepth", 0.0)) / 2)))) if self.args.vc_depth == 0 else self.args.vc_depth
                tasks.append((s, v, bam_map[s]["bam"], str(ref_map[v]), out_fa, fixed_bam, d, self.args.vc_qual, self.args.vc_freq, self.args.vc_ambig, self.args.threads, self.args.resume, log_file, self.master_log))
        
        logger.info(f"\n🧬 共识序列构建 (已自动屏蔽低覆盖度为 N) (Jobs:{self.args.jobs}, Threads/Job:{min(4, self.args.threads)})...")
        with ProcessPoolExecutor(max_workers=min(len(tasks), self.args.jobs)) as ex:
            for _ in tqdm(as_completed([ex.submit(worker_consensus, t) for t in tasks]), total=len(tasks)): pass

    def build_snpeff_db(self, accessions: list):
        jar, cfg = os.path.expanduser(self.args.snpeff_jar), os.path.expanduser(self.args.snpeff_config)
        db_path = os.path.join(os.path.dirname(cfg), "data", self.snpeff_db_name)
        
        # SnpEff 强制要求的内部文件路径
        snpeff_gbk_file = os.path.join(db_path, "genes.gbk")
        
        # 定义你想保存的输出目录
        gtf_dir = self.out / "virus-annotations"
        gtf_dir.mkdir(parents=True, exist_ok=True)
        
        logger.info(f"\n📚 为 SnpEff 下载 GenBank, 构建数据库及提取单独的 GB/GTF ({self.snpeff_db_name})...")
        os.makedirs(db_path, exist_ok=True)
        
        # 1. 批量下载所有序列，直接保存给 SnpEff 做大底库
        with open(snpeff_gbk_file, "wb") as f_out:
            for i in range(0, len(accessions), 200):
                url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={','.join(accessions[i:i+200])}&rettype=gbwithparts&retmode=text"
                try:
                    with urllib.request.urlopen(url) as response: 
                        f_out.write(response.read())
                except Exception as e: 
                    logger.error(f"GenBank 文件下载失败: {e}")
                    return False
                
        # ==========================================
        # 2. 🌟 核心升级：使用 Biopython 拆分单独的 GB 文件并提取 GTF
        # ==========================================
        if getattr(self.args, 'snpgenie', False) and os.path.exists(snpeff_gbk_file):
            try:
                from Bio import SeqIO
                # 读取刚才下载好给 SnpEff 的大文件
                for gb in SeqIO.parse(snpeff_gbk_file, 'gb'):
                    acc = gb.id
                    
                    # 🌟 修改点：单独保存这个病毒的 .gb 文件
                    individual_gb_path = gtf_dir / f"{acc}.gb"
                    SeqIO.write(gb, individual_gb_path, "genbank")
                    
                    # 提取并保存这个病毒的 .gtf 文件
                    gtf_path = gtf_dir / f"{acc}.gtf"
                    tid = 0
                    with open(gtf_path, 'w') as f_gtf:
                        for f in gb.features:
                            if f.type == 'CDS':
                                tid += 1
                                comments = []
                                keys = f.qualifiers.keys()
                                
                                comments.append(f'transcript_id "{tid}"')
                                if 'gene' in keys:
                                    comments.append(f'gene_id "{tid}"')
                                    comments.append(f'gene_name "{f.qualifiers["gene"][0]}"')
                                elif 'label' in keys and 'gene' not in keys:
                                    comments.append(f'gene_id "{tid}"')
                                    comments.append(f'gene_name "{f.qualifiers["label"][0]}"')
                                else:
                                    comments.append(f'gene_id "{tid}"')
                                    
                                strand_val = f.location.strand
                                strand = '+' if strand_val == 1 else ('-' if strand_val == -1 else '.')
                                start = int(f.location.start) + 1
                                end = int(f.location.end)
                                f_gtf.write(f"{acc}\tgb2gtf\tCDS\t{start}\t{end}\t.\t{strand}\t0\t{' ; '.join(comments)}\n")
            except ImportError:
                logger.error("❌ 缺少 Biopython 库，无法提取准确的 CDS 坐标。请运行 `pip install biopython`")
            except Exception as e:
                logger.warning(f"⚠️ 解析 GenBank 提取 GTF/GB 时出现异常: {e}")

        # 3. 将数据库写入 SnpEff config 并启动 SnpEff 建库
        try:
            # 将文件读写也放入 try 块中，捕获路径或权限错误
            with open(cfg, "r") as f:
                if f"{self.snpeff_db_name}.genome" not in f.read():
                    with open(cfg, "a") as fw: 
                        fw.write(f"\n{self.snpeff_db_name}.genome : {self.snpeff_db_name}\n")
            
            # 执行建库命令，即使失败也记录日志
            run_cmd(f"java -Xmx{self.args.snpeff_mem} -jar '{jar}' build -genbank -noCheckCds -noCheckProtein  -v -c '{cfg}' {self.snpeff_db_name}", log_path=str(self.d_logs/"snpeff_build.log"), master_log=self.master_log)
            return True
            
        except Exception as e:
            # 终极捕获：打印出到底是哪一行 Python 代码崩溃了
            err_msg = f"❌ SnpEff 建库阶段发生严重 Python/系统错误:\n{traceback.format_exc()}"
            logger.error(err_msg)
            with open(self.master_log, "a") as mlf: 
                mlf.write(f"\n{err_msg}\n")
            return False

    def run_variants_and_snpeff(self, df: pd.DataFrame, bam_map: dict, ref_map: dict):
        db_ready = self.build_snpeff_db(df["Virus"].unique().tolist()) if self.args.snpeff and self.args.variant_caller != "ivar" else False
        tasks = []
        gtf_dir = str(self.out / "virus-annotations")
        
        for _, r in df.iterrows():
            if r["Sample"] in bam_map and r["Virus"] in ref_map:
                s, v, t = r["Sample"], r["Virus"], r.get("taxonomy", "Unannotated")
                
                # 🌟 统一层级目录命名
                L1 = f"{safe_name(t).replace('.', '_')}_{safe_name(v).replace('.', '_')}"
                L2 = f"{safe_name(s).replace('.', '_')}_{safe_name(v).replace('.', '_')}"
                
                # 分别为不同步骤创建双层嵌套目录
                vdir = self.d_variants / L1 / L2
                se_dir = self.d_snpeff / L1 / L2
                vdir.mkdir(parents=True, exist_ok=True)
                se_dir.mkdir(parents=True, exist_ok=True)
                
                # 动态把 SNPGenie 也塞进这样的层级大纲中
                sg_target_dir = str(self.d_snpgenie / L1)  # worker_variants 会自动在里面创建 L2(依据 VCF 名字)
                
                ext = 'tsv' if self.args.variant_caller=='ivar' else 'vcf'
                raw_out = str(vdir / f"{L2}.variants.{ext}")
                clean_vcf = str(vdir / f"{L2}.filtered.vcf")
                fixed_bam = str(vdir / f"{L2}.fixed.bam")
                log_file = str(vdir / f"{L2}_variants.log")
                
                tasks.append((s, v, float(r.get("Recalc_MeanDepth", 0.0)), bam_map[s]["bam"], str(ref_map[v]), raw_out, clean_vcf, fixed_bam, self.args.variant_caller, db_ready, os.path.expanduser(self.args.snpeff_jar), os.path.expanduser(self.args.snpeff_config), self.snpeff_db_name, self.args.snpeff_mem, self.args.disable_dynamic_vcf, str(se_dir), getattr(self.args, 'snpgenie', False), gtf_dir, sg_target_dir, self.args.threads, self.args.resume, log_file, self.master_log))

        logger.info(f"\n🔬 变异检测与 SnpEff/SNPGenie (Jobs:{self.args.jobs}, Threads/Job:{min(4, self.args.threads)})...")
        with ProcessPoolExecutor(max_workers=min(len(tasks), self.args.jobs)) as ex:
            for _ in tqdm(as_completed([ex.submit(worker_variants, t) for t in tasks]), total=len(tasks)): pass
    
    def prepare_vodka2_dbs(self, ref_map: dict):
        logger.info(f"\n🧬 准备 VODKA2 断点数据库 (共 {len(ref_map)} 个病毒)...")
        db_dir = self.out / "virus-VODKA2-DBs"
        db_dir.mkdir(exist_ok=True)
        self.vodka2_db_map = getattr(self, 'vodka2_db_map', {})

        tasks = []
        # 🌟 智能线程分配：总线程 16，并发任务 4，则每个 bowtie2-build 分配 4 个线程。防止 CPU 爆炸！
        bt2_threads = max(1, self.args.threads // self.args.jobs)

        for virus, ref_fa in ref_map.items():
            seq_len = sum(len(line.strip()) for line in open(ref_fa) if not line.startswith('>'))
            bp = seq_len if self.args.vodka2_bp <= 0 else min(self.args.vodka2_bp, seq_len)
            rl = self.args.vodka2_readlen
            log_f = str(self.d_logs / f"vodka2_build_{virus}.log")

            # 只有当索引不存在时，才加入任务队列 (断点续传逻辑)
            cb_idx = db_dir / f"{virus}.CB.{bp}.{rl}"
            del_idx = db_dir / f"{virus}.DEL.{bp}.{rl}"
            if not (Path(f"{cb_idx}.1.bt2").exists() or Path(f"{cb_idx}.1.bt2l").exists()) or \
               not (Path(f"{del_idx}.1.bt2").exists() or Path(f"{del_idx}.1.bt2l").exists()):
                tasks.append((virus, ref_fa, bp, rl, str(db_dir), self.args.vodka2_core, bt2_threads, log_f, self.master_log))
            else:
                # 如果已经建好库了，直接存入 map
                self.vodka2_db_map[virus] = {"cb": str(cb_idx), "del": str(del_idx), "bp": bp}

        if tasks:
            # 🌟 启动多进程并发建库，并显示进度条
            with ProcessPoolExecutor(max_workers=min(len(tasks), self.args.jobs)) as ex:
                for fut in tqdm(as_completed([ex.submit(worker_build_vodka2_db, t) for t in tasks]), total=len(tasks), desc="VODKA2 建库"):
                    v, cb_idx_res, del_idx_res, actual_bp = fut.result()
                    if cb_idx_res and del_idx_res:
                        self.vodka2_db_map[v] = {"cb": cb_idx_res, "del": del_idx_res, "bp": actual_bp}
    def run_vodka2_analysis(self, df: pd.DataFrame, bam_map: dict, ref_map: dict):
        self.prepare_vodka2_dbs(ref_map)
        self.d_vodka2 = self.out / "virus-VODKA2"
        tasks = []
        for _, r in df.iterrows():
            if r["Sample"] in bam_map and r["Virus"] in self.vodka2_db_map:
                s, v, t = r["Sample"], r["Virus"], r.get("taxonomy", "Unannotated")
                if r.get("Recalc_Reads", 0) < 50: continue
                L1, L2 = f"{safe_name(t).replace('.', '_')}_{safe_name(v).replace('.', '_')}", f"{safe_name(s).replace('.', '_')}_{safe_name(v).replace('.', '_')}"
                target_dir = self.d_vodka2 / L1 / L2
                target_dir.mkdir(parents=True, exist_ok=True)
                cb_idx, del_idx, bp = self.vodka2_db_map[v]["cb"], self.vodka2_db_map[v]["del"], self.vodka2_db_map[v]["bp"]
                is_fasta = bam_map[s]["is_fasta"]
                
                # 🌟 严格打包 15 个参数
                task_args = (s, v, bam_map[s]["bam"], is_fasta, str(ref_map[v]), cb_idx, del_idx, self.args.vodka2_readlen, bp, self.args.vodka2_core, str(target_dir), self.args.threads, str(target_dir / f"{L2}_vodka2.log"), self.master_log, self.args.resume)
                tasks.append(task_args)
                
        logger.info(f"\n🚀 启动 VODKA2 DVG 鉴定 (Jobs:{self.args.jobs}, Threads/Job:{self.args.threads})...")
        with ProcessPoolExecutor(max_workers=min(len(tasks), self.args.jobs)) as ex:
            for fut in tqdm(as_completed([ex.submit(worker_vodka2, t) for t in tasks]), total=len(tasks)):
                # 🌟 第三层防弹衣：强行获取子进程结果！如果子进程崩溃，立刻在屏幕上打印红色报错！
                try:
                    fut.result() 
                except Exception as e:
                    logger.error(f"❌ VODKA2 子进程发生致命崩溃: {e}")

    def run(self):
        start = datetime.now()
        logger.info("=" * 65)
        logger.info("🚀 病毒分析终极管线 (全指标修正+极速并行版 + 盲区侦测 + SNPGenie/AF稳健版)")
        logger.info(f"📁 核心命令日志将全部输出至: {self.master_log}")
        logger.info("=" * 65)
        df = self.load_data_and_harmonize()
        ref_map = self.extract_virus_fastas(df)
        bam_map = self.resolve_bam_map(df["Sample"].unique().tolist(), ref_map)
        if not bam_map: 
            logger.error("❌ 无可用 BAM！请检查样本名与文件名的对应关系。")
            sys.exit(1)
        df = self.generate_final_summaries(df, bam_map)
        if self.args.extract_reads: self.run_extract_reads(df, bam_map)
        if self.args.consensus: self.run_consensus(df, bam_map, ref_map)
        if self.args.call_variants: self.run_variants_and_snpeff(df, bam_map, ref_map)
        if getattr(self.args, 'vodka2', False): self.run_vodka2_analysis(df, bam_map, ref_map)
        logger.info("=" * 65)
        logger.info(f"✨ 流水线完成  总耗时: {Timer._fmt((datetime.now() - start).total_seconds())}")


def main():
    parser = argparse.ArgumentParser(description="病毒宏基因组精细化后处理分析 (全能自适应版)")
    req = parser.add_argument_group("必须参数")
    req.add_argument("--summary", required=True, help="上游丰度表")
    req.add_argument("--info", required=True, help="病毒信息库 (包含 Segment 分类)")
    req.add_argument("--reference", required=True, help="超级参考基因组 FASTA")
    src_ex = parser.add_argument_group("数据来源").add_mutually_exclusive_group(required=True)
    src_ex.add_argument("--fastq", help="FASTQ/FASTA 文件夹 (格式自适应)")
    src_ex.add_argument("--bam", help="已有 BAM 文件夹")
    mod = parser.add_argument_group("功能开关")
    mod.add_argument("--extract_reads", action="store_true")
    mod.add_argument("--consensus", action="store_true")
    mod.add_argument("--call_variants", action="store_true")
    mod.add_argument("--variant_caller", choices=["freebayes", "ivar", "lofreq"], default="freebayes")
    mod.add_argument("--disable_dynamic_vcf", action="store_true")
    se = parser.add_argument_group("SnpEff 参数")
    se.add_argument("--snpeff", action="store_true")
    se.add_argument("--snpeff_jar", default="~/snpEff/snpEff.jar")
    se.add_argument("--snpeff_config", default="~/snpEff/snpEff.config")
    se.add_argument("--snpeff_mem", default="4g")
    
    se.add_argument("--snpgenie", action="store_true", help="启用 SNPGenie 进化指标分析 (依赖 perl 及 snpgenie.pl 环境变量)")
    
    vd = parser.add_argument_group("VODKA2 DVG 参数")
    vd.add_argument("--vodka2", action="store_true", help="启用 VODKA2 缺陷病毒基因组(DVG)鉴定")
    vd.add_argument("--vodka2_core", default="./VODKA2_CoreEngine.pl", help="VODKA2_CoreEngine.pl 脚本路径")
    vd.add_argument("--vodka2_readlen", type=int, default=150, help="测序 Read 长度 (默认 150)")
    vd.add_argument("--vodka2_bp", type=int, default=0, help="断点序列截取长度 (默认 0，表示自动使用病毒全长进行无死角扫描)")

    vc = parser.add_argument_group("共识参数")
    vc.add_argument("-q", "--vc_qual", type=int, default=20)
    vc.add_argument("-d", "--vc_depth", type=int, default=0)
    vc.add_argument("-f", "--vc_freq", type=float, default=0.5)
    vc.add_argument("-a", "--vc_ambig", type=str, default="N")
    ctl = parser.add_argument_group("控制参数")
    ctl.add_argument("--output_dir", default="./post_analysis")
    ctl.add_argument("--threads", type=int, default=8)
    ctl.add_argument("--jobs", type=int, default=4)
    ctl.add_argument("--resume", action="store_true")
    args = parser.parse_args()
    PostProcessPipeline(args).run()

if __name__ == "__main__": main()
