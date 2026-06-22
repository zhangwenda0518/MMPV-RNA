#!/usr/bin/env python3
"""
batch_virus_downstream.py — 病毒确诊后处理、节段病毒优化与变异分析终极管线
【VAP专属表头 | 真实测序深度全指标修正 | 盲区侦测 | Shannon/AF & SNPGenie 稳健版 | 融合 DI-tector 终极版】
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
import tempfile
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime
from pathlib import Path

import colorlog
import numpy as np
import pandas as pd
import pysam
from tqdm import tqdm


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
    full_cmd = f"set -o pipefail; {cmd}"
    start_time = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    
    result = subprocess.run(full_cmd, shell=True, executable="/bin/bash", capture_output=True, text=True)
    
    log_content = f"\n[{start_time}] CMD: {cmd}\nEXIT_CODE: {result.returncode}\n"
    if result.stdout: log_content += f"--- STDOUT ---\n{result.stdout.strip()}\n"
    if result.stderr: log_content += f"--- STDERR ---\n{result.stderr.strip()}\n"
    log_content += "-" * 80 + "\n"
    
    if log_path:
        lp = Path(log_path)
        lp.parent.mkdir(parents=True, exist_ok=True)
        with open(lp, "a", encoding="utf-8") as f: f.write(log_content)
            
    if master_log:
        ml = Path(master_log)
        ml.parent.mkdir(parents=True, exist_ok=True)
        with open(ml, "a", encoding="utf-8") as f: f.write(log_content)
            
    if check and result.returncode != 0:
        raise subprocess.CalledProcessError(result.returncode, cmd, output=result.stdout, stderr=result.stderr)
    return result


# ==========================================
# 2. 功能辅助工具
# ==========================================
def extract_allele_frequency(vcf_path, out_tsv):
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
    try:
        with open(fasta_path, 'r') as f:
            header = f.readline().strip()[1:].split()[0]
            seq = "".join([line.strip() for line in f])
        length = len(seq)
        with open(gtf_path, 'w') as f:
            f.write(f"{header}\tVAP_Fallback\tCDS\t1\t{length}\t.\t+\t0\tgene_id \"{header}_CDS\"; transcript_id \"{header}_TX\";\n")
    except Exception as e:
        logger.warning(f"生成虚拟 GTF 失败: {e}")

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


# ==========================================
# 3. 多进程 Worker 函数群
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

def worker_extract_reads(args):
    sample, virus, bam_path, is_fasta, r1_str, r2_str, rs_str, threads, resume, log_file, master_log = args
    virus = str(virus).strip()
    bam_path = str(bam_path).strip()
    r1, r2, rs = Path(r1_str), Path(r2_str), Path(rs_str)

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
            cmd_get_names = f"samtools view -@ {t_io} '{bam_path}' '{virus}' | cut -f1 | awk '!seen[$0]++' > '{tmp_names}'"
            res_names = subprocess.run(cmd_get_names, shell=True, executable='/bin/bash', capture_output=True, text=True)

            if res_names.returncode != 0:
                print(f"\n❌ [名单提取报错] {sample}: {res_names.stderr}")
                return False

            if not tmp_names.exists() or tmp_names.stat().st_size == 0:
                return True 

            cmd_ext = (
                f"samtools view -@ {t_io} -h -N '{tmp_names}' '{bam_path}' | "
                f"samtools collate -O -u -@ {t_io} - | "
                f"{ext_tool} -@ {t_io} -1 '{tmp_r1}' -2 '{tmp_r2}' -s '{tmp_rs}' -0 '{tmp_0}' -n -"
            )
            res_ext = subprocess.run(cmd_ext, shell=True, executable='/bin/bash', capture_output=True, text=True)

            if res_ext.returncode != 0:
                print(f"\n❌ [序列提取报错] {sample}: {res_ext.stderr}")
                return False

            if tmp_0.exists() and tmp_0.stat().st_size > 0:
                subprocess.run(f"cat '{tmp_0}' >> '{tmp_rs}'", shell=True, executable='/bin/bash')

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
        af_tsv = Path(raw_out_str).parent / f"{sample}.{virus}.allele_frequencies.tsv"
        extract_allele_frequency(str(clean_vcf), str(af_tsv))
        
        if db_ready:
            se_dir = Path(se_dir_str)
            ann_vcf, sum_tsv = se_dir / f"{clean_vcf.name.replace('.filtered.vcf', '.ann.vcf')}", se_dir / f"{clean_vcf.name.replace('.filtered.vcf', '.annotation_summary.tsv')}"
            cmd_ann = f"java -Xmx{snpeff_mem} -jar '{snpeff_jar}' ann -c '{snpeff_config}' -noStats {snpeff_db_name} '{clean_vcf}' > '{ann_vcf}'"
            try: 
                run_cmd(cmd_ann, log_path=log_file, master_log=master_log)
                parse_ann_to_tsv(str(ann_vcf), str(sum_tsv))
            except Exception as e:
                with open(master_log, "a") as mlf: mlf.write(f"\n[Python Error] SnpEff failed: {e}\n")
            
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
                run_cmd(cmd_sg, log_path=log_file, master_log=master_log, check=False)
                for f in [vcf_work, fa_work, gtf_work]: f.unlink(missing_ok=True)
            except Exception as e:
                err = f"\n[Python Exception] SNPGenie failed: {traceback.format_exc()}\n"
                with open(log_file, "a") as lf: lf.write(err)
                with open(master_log, "a") as mlf: mlf.write(err)
                
    return True

def worker_ditector(args):
    (sample, virus, bam_path, ref_fa, target_dir_str, script_path, threads, resume, log_file, master_log) = args
    target_dir = Path(target_dir_str)

    counts_out = target_dir / f"{sample}_{virus}_counts.txt"
    if resume and counts_out.exists() and counts_out.stat().st_size > 0: return True

    t_io = min(4, threads)
    cov_txt = target_dir / f"{sample}_{virus}_true_coverage.txt"
    # 🌟 修改点 1：将扩展名改为 .fa
    subset_fa = target_dir / f"{sample}_{virus}_abnormal.fa"
    tmp_sam = target_dir / f"{sample}_{virus}_extract.tmp.sam"

    try:
        # 步骤 1: 精准提取目标病毒深度分布
        depth_cmd = f"samtools depth -a -r '{virus}' '{bam_path}' > '{cov_txt}'"
        run_cmd(depth_cmd, log_path=log_file, master_log=master_log)

        # 步骤 2: 极速提取策略 (异常比对reads + 完全没比对上的reads)
        awk_extract = f"awk -v v='{virus}' '$1~/^@/ || $6~/S|H/'"
        cmd_extract_1 = f"samtools view -@ {t_io} -h '{bam_path}' '{virus}' | {awk_extract} > '{tmp_sam}'"
        cmd_extract_2 = f"samtools view -@ {t_io} '{bam_path}' '*' >> '{tmp_sam}'"

        # 🌟 修改点 2：改用 samtools fasta，并使用终端重定向 '>'，彻底阻止 reads 流入 log
        cmd_fastx = f"samtools fasta -@ {t_io} '{tmp_sam}' > '{subset_fa}'"

        run_cmd(cmd_extract_1, log_path=log_file, master_log=master_log)
        run_cmd(cmd_extract_2, log_path=log_file, master_log=master_log)
        run_cmd(cmd_fastx, log_path=log_file, master_log=master_log)

        # 步骤 3: 运行 DI-tector_06_ProMax.py
        python_bin = sys.executable
        di_tag = f"{sample}_{virus}"
        # 🌟 修改点 3：输入改为 .fa，并加入核心开关 `-f` 告诉 DI-tector 这是 FASTA 格式
        di_cmd = f"{python_bin} '{script_path}' '{ref_fa}' '{subset_fa}' -f -o '{target_dir}' -t '{di_tag}' -x {threads} -d -n 3 -c '{cov_txt}'"
        run_cmd(di_cmd, log_path=log_file, master_log=master_log, check=False)

        return True
    except Exception as e:
        err = f"\n[Python Exception] worker_ditector: {traceback.format_exc()}\n"
        with open(log_file, "a") as lf: lf.write(err)
        with open(master_log, "a") as mlf: mlf.write(err)
        return False
    finally:
        if tmp_sam.exists(): tmp_sam.unlink()
        if subset_fa.exists(): subset_fa.unlink()

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
        self.d_ditector    = self.out / "virus-DI-tector" 
        self.d_summary     = self.out / "summary"
        self.d_ind_reports = self.out / "summary" / "individual_virus_reports"
        self.d_logs        = self.out / "logs"

        self.index_prefix: str = ""
        self.harmonized_fasta = self.d_index / "harmonized_references.fasta"
        self.snpeff_db_name = f"vap_db_{datetime.now():%Y%m%d_%H%M%S}"
        
        self.master_log = str(self.d_logs / "master_commands.log")

        self._validate_args(); self._check_tools(); self._setup_dirs(); self._add_file_logger()

    def _validate_args(self):
        a = self.args
        for f in (a.summary, a.reference, a.info):
            if not Path(f).exists(): logger.error(f"❌ 文件不存在: {f}"); sys.exit(1)
        if a.fastq and not Path(a.fastq).is_dir(): logger.error(f"❌ --fastq 目录不存在: {a.fastq}"); sys.exit(1)
        if a.ditector and not Path(a.ditector_script).exists():
            logger.error(f"❌ 找不到 DI-tector 脚本，请检查路径: {a.ditector_script}"); sys.exit(1)

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
        if self.args.ditector: need["bwa"] = "DI-tector 核心要求" 
            
        missing = [f"{t} ({d})" for t, d in need.items() if not shutil.which(t)]
        if missing: logger.error("❌ 缺少必要工具:\n  " + "\n  ".join(f"- {m}" for m in missing)); sys.exit(1)

    def _setup_dirs(self):
        dirs = [self.d_fasta, self.d_summary, self.d_ind_reports, self.d_logs, self.d_bam]
        if self.args.fastq: dirs.append(self.d_index)
        if self.args.extract_reads: dirs.append(self.d_reads)
        if self.args.consensus: dirs.append(self.d_consensus)
        if self.args.ditector: dirs.append(self.d_ditector) 
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

            # 🌟 新增：强制为单个参考基因组生成 BWA 索引，专供 DI-tector 使用
            if not Path(str(ref_fa) + ".bwt").exists():
                run_cmd(f"bwa index '{ref_fa}'", log_path=str(self.d_logs/"bwa_index.log"), master_log=self.master_log)

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
            logger.info(f"\n🚀 开始精准靶向重比对 (Jobs:{self.args.jobs}, Threads/Job:{self.args.threads}) ...")
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
                L1 = f"{safe_name(t).replace('.', '_')}_{safe_name(v).replace('.', '_')}"
                L2 = f"{safe_name(s).replace('.', '_')}_{safe_name(v).replace('.', '_')}"
                target_dir = self.d_reads / L1 / L2
                target_dir.mkdir(parents=True, exist_ok=True)
                ext = 'fasta.gz' if isf else 'fastq.gz'
                r1, r2, rs = str(target_dir / f"{L2}_R1.{ext}"), str(target_dir / f"{L2}_R2.{ext}"), str(target_dir / f"{L2}_single.{ext}")
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
                L1 = f"{safe_name(t).replace('.', '_')}_{safe_name(v).replace('.', '_')}"
                L2 = f"{safe_name(s).replace('.', '_')}_{safe_name(v).replace('.', '_')}"
                vdir = self.d_consensus / L1 / L2
                vdir.mkdir(parents=True, exist_ok=True)
                out_fa, fixed_bam = str(vdir / f"{L2}.consensus.fasta"), str(vdir / f"{L2}.fixed.bam")
                log_file = str(vdir / f"{L2}_consensus.log")
                d = max(1, min(10, int(math.floor(float(r.get("Recalc_MeanDepth", 0.0)) / 2)))) if self.args.vc_depth == 0 else self.args.vc_depth
                tasks.append((s, v, bam_map[s]["bam"], str(ref_map[v]), out_fa, fixed_bam, d, self.args.vc_qual, self.args.vc_freq, self.args.vc_ambig, self.args.threads, self.args.resume, log_file, self.master_log))
        
        logger.info(f"\n🧬 共识序列构建 (已自动屏蔽低覆盖度为 N) (Jobs:{self.args.jobs}, Threads/Job:{min(4, self.args.threads)})...")
        with ProcessPoolExecutor(max_workers=min(len(tasks), self.args.jobs)) as ex:
            for _ in tqdm(as_completed([ex.submit(worker_consensus, t) for t in tasks]), total=len(tasks)): pass

    def build_snpeff_db(self, accessions: list):
        jar, cfg = os.path.expanduser(self.args.snpeff_jar), os.path.expanduser(self.args.snpeff_config)
        db_path = os.path.join(os.path.dirname(cfg), "data", self.snpeff_db_name)
        snpeff_gbk_file = os.path.join(db_path, "genes.gbk")
        gtf_dir = self.out / "virus-annotations"
        gtf_dir.mkdir(parents=True, exist_ok=True)

        logger.info(f"\n📚 为 SnpEff 下载 GenBank, 构建数据库及提取单独的 GB/GTF ({self.snpeff_db_name})...")
        os.makedirs(db_path, exist_ok=True)

        with open(snpeff_gbk_file, "wb") as f_out:
            for i in range(0, len(accessions), 200):
                url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={','.join(accessions[i:i+200])}&rettype=gbwithparts&retmode=text"
                try:
                    with urllib.request.urlopen(url) as response:
                        f_out.write(response.read())
                except Exception as e:
                    logger.error(f"GenBank 文件下载失败: {e}")
                    return False

        if getattr(self.args, 'snpgenie', False) and os.path.exists(snpeff_gbk_file):
            try:
                from Bio import SeqIO
                for gb in SeqIO.parse(snpeff_gbk_file, 'gb'):
                    acc = gb.id
                    individual_gb_path = gtf_dir / f"{acc}.gb"
                    SeqIO.write(gb, individual_gb_path, "genbank")

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
                logger.error("❌ 缺少 Biopython 库，运行 `pip install biopython` 以获取精确 GTF。")
            except Exception as e:
                logger.warning(f"⚠️ 解析 GenBank 提取 GTF/GB 异常: {e}")

        try:
            with open(cfg, "r") as f:
                if f"{self.snpeff_db_name}.genome" not in f.read():
                    with open(cfg, "a") as fw:
                        fw.write(f"\n{self.snpeff_db_name}.genome : {self.snpeff_db_name}\n")

            run_cmd(f"java -Xmx{self.args.snpeff_mem} -jar '{jar}' build -genbank -noCheckCds -noCheckProtein  -v -c '{cfg}' {self.snpeff_db_name}", log_path=str(self.d_logs/"snpeff_build.log"), master_log=self.master_log)
            return True
        except Exception as e:
            err_msg = f"❌ SnpEff 建库发生异常:\n{traceback.format_exc()}"
            logger.error(err_msg)
            with open(self.master_log, "a") as mlf: mlf.write(f"\n{err_msg}\n")
            return False

    def run_variants_and_snpeff(self, df: pd.DataFrame, bam_map: dict, ref_map: dict):
        db_ready = self.build_snpeff_db(df["Virus"].unique().tolist()) if self.args.snpeff and self.args.variant_caller != "ivar" else False
        tasks = []
        gtf_dir = str(self.out / "virus-annotations")

        for _, r in df.iterrows():
            if r["Sample"] in bam_map and r["Virus"] in ref_map:
                s, v, t = r["Sample"], r["Virus"], r.get("taxonomy", "Unannotated")
                L1 = f"{safe_name(t).replace('.', '_')}_{safe_name(v).replace('.', '_')}"
                L2 = f"{safe_name(s).replace('.', '_')}_{safe_name(v).replace('.', '_')}"

                vdir = self.d_variants / L1 / L2
                se_dir = self.d_snpeff / L1 / L2
                vdir.mkdir(parents=True, exist_ok=True)
                se_dir.mkdir(parents=True, exist_ok=True)

                sg_target_dir = str(self.d_snpgenie / L1)

                ext = 'tsv' if self.args.variant_caller=='ivar' else 'vcf'
                raw_out = str(vdir / f"{L2}.variants.{ext}")
                clean_vcf = str(vdir / f"{L2}.filtered.vcf")
                fixed_bam = str(vdir / f"{L2}.fixed.bam")
                log_file = str(vdir / f"{L2}_variants.log")

                tasks.append((s, v, float(r.get("Recalc_MeanDepth", 0.0)), bam_map[s]["bam"], str(ref_map[v]), raw_out, clean_vcf, fixed_bam, self.args.variant_caller, db_ready, os.path.expanduser(self.args.snpeff_jar), os.path.expanduser(self.args.snpeff_config), self.snpeff_db_name, self.args.snpeff_mem, self.args.disable_dynamic_vcf, str(se_dir), getattr(self.args, 'snpgenie', False), gtf_dir, sg_target_dir, self.args.threads, self.args.resume, log_file, self.master_log))

        logger.info(f"\n🔬 变异检测与 SnpEff/SNPGenie (Jobs:{self.args.jobs}, Threads/Job:{min(4, self.args.threads)})...")
        with ProcessPoolExecutor(max_workers=min(len(tasks), self.args.jobs)) as ex:
            for _ in tqdm(as_completed([ex.submit(worker_variants, t) for t in tasks]), total=len(tasks)): pass

    # 🌟 新增：批量运行 DI-tector 的分发函数
    def run_ditector(self, df: pd.DataFrame, bam_map: dict, ref_map: dict):
        tasks = []
        for _, r in df.iterrows():
            if r["Sample"] in bam_map and r["Virus"] in ref_map:
                s, v, t = r["Sample"], r["Virus"], r.get("taxonomy", "Unannotated")

                # 遵循优雅的嵌套结构规则
                L1 = f"{safe_name(t).replace('.', '_')}_{safe_name(v).replace('.', '_')}"
                L2 = f"{safe_name(s).replace('.', '_')}_{safe_name(v).replace('.', '_')}"

                target_dir = self.d_ditector / L1 / L2
                target_dir.mkdir(parents=True, exist_ok=True)

                log_file = str(target_dir / f"{L2}_ditector.log")

                tasks.append((
                    s, v, bam_map[s]["bam"], str(ref_map[v]), str(target_dir),
                    os.path.abspath(self.args.ditector_script),
                    self.args.threads, self.args.resume, log_file, self.master_log
                ))

        logger.info(f"\n🧬 DI-tector DVG重组探测 + 非零中位数丰度精准校正 (Jobs:{self.args.jobs}, Threads/Job:{min(4, self.args.threads)})...")
        with ProcessPoolExecutor(max_workers=min(len(tasks), self.args.jobs)) as ex:
            for _ in tqdm(as_completed([ex.submit(worker_ditector, t) for t in tasks]), total=len(tasks)): pass

    def run(self):
        start = datetime.now()
        logger.info("=" * 65)
        logger.info("🚀 病毒宏基因组终极分析管线 (全指标修正 + SNPGenie + DI-tector_ProMax)")
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

        # 🌟 触发 DI-tector 分析流程
        if self.args.ditector: self.run_ditector(df, bam_map, ref_map)

        logger.info("=" * 65)
        logger.info(f"✨ 流水线绝密级分析完成  总耗时: {Timer._fmt((datetime.now() - start).total_seconds())}")


def main():
    parser = argparse.ArgumentParser(description="病毒宏基因组精细化后处理分析 (全能自适应超级融合版)")
    req = parser.add_argument_group("必须参数")
    req.add_argument("--summary", required=True, help="上游丰度表")
    req.add_argument("--info", required=True, help="病毒信息库 (包含 Segment 分类)")
    req.add_argument("--reference", required=True, help="超级参考基因组 FASTA")

    src_ex = parser.add_argument_group("数据来源").add_mutually_exclusive_group(required=True)
    src_ex.add_argument("--fastq", help="FASTQ/FASTA 文件夹 (格式自适应)")
    src_ex.add_argument("--bam", help="已有 BAM 文件夹")

    mod = parser.add_argument_group("基础功能开关")
    mod.add_argument("--extract_reads", action="store_true")
    mod.add_argument("--consensus", action="store_true")
    mod.add_argument("--call_variants", action="store_true")
    mod.add_argument("--variant_caller", choices=["freebayes", "ivar", "lofreq"], default="freebayes")
    mod.add_argument("--disable_dynamic_vcf", action="store_true")

    # 🌟 新增 DI-tector 面板
    dt = parser.add_argument_group("DI-tector 缺陷型重组病毒(DVG)挖掘面板")
    dt.add_argument("--ditector", action="store_true", help="启用 DI-tector 进行 DVG 检测与全深度精准丰度定量")
    dt.add_argument("--ditector_script", default="./DI-tector_06_ProMax.py", help="你的 DI-tector_06_ProMax.py 脚本路径 (默认当前目录下寻找)")

    se = parser.add_argument_group("SnpEff & 进化生信图谱参数")
    se.add_argument("--snpeff", action="store_true")
    se.add_argument("--snpeff_jar", default="~/snpEff/snpEff.jar")
    se.add_argument("--snpeff_config", default="~/snpEff/snpEff.config")
    se.add_argument("--snpeff_mem", default="4g")
    se.add_argument("--snpgenie", action="store_true", help="启用 SNPGenie (依赖 perl及环境)")

    vc = parser.add_argument_group("共识参数")
    vc.add_argument("-q", "--vc_qual", type=int, default=20)
    vc.add_argument("-d", "--vc_depth", type=int, default=0)
    vc.add_argument("-f", "--vc_freq", type=float, default=0.5)
    vc.add_argument("-a", "--vc_ambig", type=str, default="N")

    ctl = parser.add_argument_group("工程化与并行控制参数")
    ctl.add_argument("--output_dir", default="./post_analysis")
    ctl.add_argument("--threads", type=int, default=8, help="单任务调用的工具线程数 (如bowtie2, samtools)")
    ctl.add_argument("--jobs", type=int, default=4, help="同时并行处理的最大样本/病毒任务数")
    ctl.add_argument("--resume", action="store_true", help="断点续传机制，跳过已生成文件的步骤")

    args = parser.parse_args()
    PostProcessPipeline(args).run()

if __name__ == "__main__": main()
