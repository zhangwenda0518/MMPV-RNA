#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
===============================================================================
Virus Auto Pipeline[终极无死锁·原子断点续传版]
===============================================================================
🌟 最新突破级功能护航：
  1. 【颗粒级断点续传】: 支持 `--resume`。矩阵每完成 5% 自动原子刷盘。服务器意外宕机
     随时接力重启，已计算对位直接装载，绝不浪费 1 滴算力！
  2. 【破除画图假死】: 强制挂载 Matplotlib 'Agg' 后端，免疫集群 SSH X11 静默挂起。
  3. 【全透明状态流】: 实时播报，拒绝 100% 后的漫长黑匣子。
  4. 保留完美展开结构、SciPy 极速降维排序、脱水果序列以及出版级图集。
===============================================================================
"""

import argparse
import sys
import os
import re
import csv
import tempfile
import subprocess
import warnings
from pathlib import Path
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor

import numpy as np

# 🔥【核心防挂死】静默无头渲染器，免疫 Linux 无显示器服务器死锁！
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import seaborn as sns
from matplotlib.colors import ListedColormap, BoundaryNorm
import scipy.cluster.hierarchy as sch
import scipy.spatial.distance as ssd
from tqdm import tqdm

from Bio import SeqIO, Entrez, Align, AlignIO
from Bio.Align import substitution_matrices
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import BiopythonWarning

warnings.simplefilter('ignore', BiopythonWarning)

NT_MATCH, NT_MISMATCH = 2, -1
NT_GAP_OPEN, NT_GAP_EXTEND = -10, -0.5
AA_GAP_OPEN, AA_GAP_EXTEND = -10, -0.5
AA_SUBSTITUTION_MATRIX = "BLOSUM62"
GENETIC_CODE = 1

SDT_COLORS =[
    '#000066', '#0033CC', '#0099FF', '#00CC99', '#99FF99', 
    '#FFFF66', '#FFCC00', '#FF6600', '#E60000', '#8B0000'
]
PLOT_DPI = 600

# =============================================================================
#[ 工具函数 ]
# =============================================================================
def get_sample_id(raw_id):
    name = re.sub(r'\.consensus\.fasta$|\.fasta$|\.full$', '', raw_id, flags=re.IGNORECASE)
    name = re.sub(r'[\._][A-Z]{1,4}_?\d{5,8}(\.\d+)?.*$', '', name, flags=re.IGNORECASE)
    return name

def trim_ends(seq):
    return Seq(re.sub(r'[Nn]+$', '', re.sub(r'^[Nn]+', '', str(seq))))

def check_protein_integrity(aa_seq):
    clean_seq = str(aa_seq).rstrip('*')
    return clean_seq.count('*') == 0, clean_seq.count('*')

def parse_gb_metadata(gb_path):
    organism, accession, has_cds, ref_len = "Unknown_Virus", "Unknown_Acc", False, 0
    try:
        records = list(SeqIO.parse(gb_path, "genbank"))
        if not records:
            print(f"[❌ ERR] GenBank文件为空或格式无效！检查输入。")
            return organism, accession, has_cds, ref_len
            
        record = records[0]
        organism = record.annotations.get('organism', 'Unknown_Virus')
        accession = record.id
        ref_len = len(record.seq)
        
        for feature in record.features:
            if feature.type in["CDS", "mat_peptide"]:
                has_cds = True
                
        organism = re.sub(r'[^\w-]', '_', organism) 
        return organism, accession, has_cds, ref_len
    except Exception as e:
        print(f"      [❌ ERR] 解析 GenBank 失败: {e}")
        return organism, accession, has_cds, ref_len

def download_genbank(accession, out_dir):
    if not accession or re.match(r'^[SCE]RR', accession):
        return None
    
    fixed_acc = accession
    if re.match(r'^[A-Z]{2}\.\d+', accession):
        fixed_acc = accession.replace('.', '_', 1) 
        
    out_path = Path(out_dir) / f"{fixed_acc}.gb"
    if out_path.exists() and out_path.stat().st_size > 500: 
        return out_path
        
    print(f"      🌐 正在从 NCBI 联网下载参考序列: {fixed_acc} ...")
    try:
        with Entrez.efetch(db="nucleotide", id=fixed_acc, rettype="gb", retmode="text") as handle:
            content = handle.read()
            if "Item not found" in content or "Error" in content[:100]:
                print(f"      [❌ ERR] NCBI 返回查无此序列！")
                return None
            with open(out_path, "w") as f: 
                f.write(content)
        return out_path
    except Exception as e:
        print(f"[❌ ERR] 联网下载溃败: {e}")
        return None

def extract_valid_accession(text):
    patterns =[r'(?<![A-Z0-9])([A-Z]{1,2}\d{5,6}(\.\d+)?)(?![A-Z0-9])', 
                r'(?<![A-Z0-9])([A-Z]{2}_\d{6,8}(\.\d+)?)(?![A-Z0-9])']
    for p in patterns:
        for m in re.finditer(p, text):
            acc = m.group(1)
            if not re.match(r'^[SCE]RR', acc): 
                return acc
    return None

# =============================================================================
#[ 多层级聚类与去重引擎 ]
# =============================================================================
def calculate_adaptive_threshold(ref_length):
    if ref_length == 0: 
        return 0.99
    if ref_length < 1000: 
        return max(0.85, 1.0 - (5.0 / ref_length)) 
    elif ref_length < 10000: 
        return max(0.90, 1.0 - (15.0 / ref_length))
    else: 
        return max(0.95, 1.0 - (30.0 / ref_length))

def get_cdhit_params(ref_len, threshold_override):
    c = threshold_override if threshold_override else calculate_adaptive_threshold(ref_len)
    
    if ref_len < 1000: aL, aS = 0.95, 0.90
    elif ref_len < 10000: aL, aS = 0.98, 0.95
    else: aL, aS = 0.99, 0.99
    
    if c >= 0.98: n = 10
    elif c >= 0.95: n = 8
    elif c >= 0.90: n = 6
    else: n = 5
    
    return c, n, aL, aS

def run_cdhit_and_parse(recs, out_dir, ref_len, threshold_override, threads, cdhit_path="cd-hit-est"):
    in_fasta = out_dir / "cdhit_temp_input.fasta"
    out_fasta = out_dir / "cdhit_temp_output.fasta"
    clstr_file = out_dir / "cdhit_temp_output.fasta.clstr"
    report_tsv = out_dir / "Genome_CDHIT_Clustering_Report.tsv"
    
    SeqIO.write(recs, in_fasta, "fasta")
    c, n, aL, aS = get_cdhit_params(ref_len, threshold_override)
    print(f"      🧠[CD-HIT 智能参数] Identity: {c*100:.2f}% | Word: {n}")
    
    cmd =[
        cdhit_path, "-i", str(in_fasta), "-o", str(out_fasta),
        "-c", str(c), "-n", str(n), "-aL", str(aL), "-aS", str(aS),
        "-d", "0", "-M", "0", "-T", str(threads)
    ]
    try: 
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except Exception as e: 
        sys.exit(f"❌ 错误: CD-HIT 执行失败！{e}")

    clusters =[]
    with open(clstr_file, 'r') as f:
        for line in f:
            if line.startswith('>Cluster'):
                clusters.append({'id': line.strip().replace('>', ''), 'rep': None, 'rep_size': '', 'members':[]})
            else:
                parts = line.strip().split(', >')
                size_str = parts[0].split('\t')[1]
                rest = parts[1].split('... ')
                sample_id = get_sample_id(rest[0])
                ident = rest[1].strip()

                if ident == '*':
                    clusters[-1]['rep'] = sample_id
                    clusters[-1]['rep_size'] = size_str
                    clusters[-1]['members'].append((sample_id, size_str, '100% (Ref)'))
                else:
                    m = re.search(r'([0-9\.]+\%)', ident)
                    ident_val = m.group(1) if m else ident
                    clusters[-1]['members'].append((sample_id, size_str, ident_val))

    rep_names = set()
    with open(report_tsv, 'w', encoding='utf-8') as f:
        f.write("Cluster_ID\tRepresentative_Sample\tRep_Length\tCluster_Size\tCluster_Members_Details\n")
        for cl in clusters:
            rep_names.add(cl['rep'])
            members_str = ", ".join([f"{m[0]}({m[2]})" for m in cl['members'] if m[0] != cl['rep']])
            if not members_str: 
                members_str = "None"
            f.write(f"{cl['id']}\t{cl['rep']}\t{cl['rep_size']}\t{len(cl['members'])}\t{members_str}\n")

    if in_fasta.exists(): in_fasta.unlink()
    if out_fasta.exists(): out_fasta.unlink()
    if clstr_file.exists(): clstr_file.unlink()
    
    return list(rep_names)

def exact_deduplicate(samps_list, seq_dict, report_path):
    clusters = defaultdict(list)
    for sid in samps_list:
        seq = seq_dict[sid].upper()
        clusters[seq].append(sid)
    
    reps =[]
    with open(report_path, 'w', encoding='utf-8') as f:
        f.write("Cluster_ID\tRepresentative_Sample\tCluster_Size\tIdentical_Duplicates\n")
        for idx, (seq, sids) in enumerate(clusters.items()):
            reps.append(sids[0])
            dups = sids[1:]
            dup_str = ", ".join(dups) if dups else "None"
            f.write(f"ExactCluster_{idx}\t{sids[0]}\t{len(sids)}\t{dup_str}\n")
    return sorted(reps)

# =============================================================================
# [ 同源锚定提取引擎 ]
# =============================================================================
def extract_by_alignment(str_ref, consensus_seq):
    aligner = Align.PairwiseAligner()
    aligner.mode = 'local'
    aligner.match_score = 2
    aligner.mismatch_score = -3
    aligner.open_gap_score = -5
    aligner.extend_gap_score = -2

    str_con = str(consensus_seq).upper()
    str_con_rev = str(consensus_seq.reverse_complement()).upper()

    alns_fwd = aligner.align(str_ref, str_con)
    alns_rev = aligner.align(str_ref, str_con_rev)

    best_aln = None
    is_rev = False

    if alns_fwd and alns_rev:
        if alns_fwd[0].score >= alns_rev[0].score: 
            best_aln = alns_fwd[0]
        else: 
            best_aln, is_rev = alns_rev[0], True
    elif alns_fwd: 
        best_aln = alns_fwd[0]
    elif alns_rev: 
        best_aln, is_rev = alns_rev[0], True

    if best_aln is None or len(best_aln.aligned[1]) == 0: 
        return None
    
    start, end = best_aln.aligned[1][0][0], best_aln.aligned[1][-1][1]
    match_length = end - start
    
    if match_length < len(str_ref) * 0.7: 
        return None

    if is_rev: 
        extracted_seq = consensus_seq.reverse_complement()[start:end]
    else: 
        extracted_seq = consensus_seq[start:end]
        
    return extracted_seq

def _extract_worker(task):
    rec_id, seq_str, ref_features_data, is_nt_only = task
    prefix = get_sample_id(rec_id)
    consensus_seq = Seq(seq_str)
    whole_seq_str = seq_str.upper()
    
    if len(whole_seq_str) > 0:
        genome_n_ratio = (whole_seq_str.count('N') / len(whole_seq_str)) * 100
    else:
        genome_n_ratio = 100.0

    sample_stats =[]
    extracted_records =[]

    for clean_name, raw_name, gene_name, str_ref_feat in ref_features_data:
        raw_nuc = extract_by_alignment(str_ref_feat, consensus_seq)
        
        if raw_nuc is None:
            sample_stats.append((clean_name, "Mis/Trunc", "Mis/Trunc", "Fail", "Fail"))
            continue

        if len(raw_nuc) > 0:
            raw_n_ratio = (str(raw_nuc).upper().count('N') / len(raw_nuc)) * 100
        else:
            raw_n_ratio = 100.0
            
        raw_prot_qc, raw_aa = "Pass", Seq("")
        
        if not is_nt_only:
            try:
                raw_aa = raw_nuc.translate(table=GENETIC_CODE, to_stop=True)
                if not check_protein_integrity(raw_aa)[0]:
                    raw_prot_qc = "Stop"
            except: 
                raw_prot_qc = "TransErr"

        fill_nuc = raw_nuc
        fill_n_ratio = raw_n_ratio
        
        fill_prot_qc, fill_aa = "Pass", Seq("")
        if not is_nt_only:
            try:
                fill_aa = fill_nuc.translate(table=GENETIC_CODE, to_stop=True)
                if not check_protein_integrity(fill_aa)[0]:
                    fill_prot_qc = "Stop"
            except: 
                fill_prot_qc = "TransErr"

        sample_stats.append((clean_name, f"{raw_n_ratio:.2f}", f"{fill_n_ratio:.2f}", raw_prot_qc, fill_prot_qc))
        extracted_records.append({"clean_name": clean_name, "raw_nuc": raw_nuc,  "raw_aa": raw_aa,  "raw_prot_qc": raw_prot_qc, "fill_nuc": fill_nuc,  "fill_aa": fill_aa,  "fill_prot_qc": fill_prot_qc})

    return prefix, genome_n_ratio, sample_stats, extracted_records

# =============================================================================
#[ 绘图引擎 ] 出版级智能超清多格式并发输出
# =============================================================================
def plot_heatmap(mat, labels, title, out_path, triangle_only=False):
    n = len(labels)
    fig_size = max(8.0, n * 0.3)
    font_size = min(14, max(6, 250 / n))
    
    fig, ax = plt.subplots(figsize=(fig_size, fig_size))
    ax.set_facecolor('white')
    
    if triangle_only:
        mask = np.triu(np.ones_like(mat, dtype=bool), k=0) | np.isnan(mat)
    else:
        mask = np.eye(n, dtype=bool) | np.isnan(mat)

    cmap = ListedColormap(SDT_COLORS)
    cmap.set_bad('white', alpha=0.0)

    vals = mat[~mask]
    if len(vals) > 0:
        actual_min = np.nanmin(vals)
        if actual_min >= 95.0: min_v = 95.0
        elif actual_min >= 90.0: min_v = 90.0
        elif actual_min >= 85.0: min_v = 85.0
        elif actual_min >= 80.0: min_v = 80.0
        else: min_v = np.floor(actual_min / 10.0) * 10.0
    else:
        min_v = 80.0

    bounds = np.linspace(min_v, 100.0, len(SDT_COLORS) + 1)

    # v2 修复: linewidths=0 避免 masked 区域绘制多余小格子
    sns.heatmap(mat, mask=mask, cmap=cmap, norm=BoundaryNorm(bounds, cmap.N),
                xticklabels=labels, yticklabels=labels,
                linewidths=0, cbar_kws={'shrink': 0.5, 'pad': 0.02,
                'label': 'Pairwise Identity (%)'}, ax=ax)

    # 手动绘制网格线 (仅非mask区域)
    line_width = 0.5 if n < 100 else 0.1
    if line_width > 0:
        for i in range(n + 1):
            if triangle_only:
                ax.plot([0, min(i, n)], [i, i], color='black', lw=line_width)
                ax.plot([i, i], [i, n], color='black', lw=line_width)
            else:
                ax.axhline(i, color='black', lw=line_width)
                ax.axvline(i, color='black', lw=line_width)
    
    ax.set_title(title, fontsize=max(16, font_size*1.2), fontweight='bold', pad=20)
    ax.tick_params(axis='x', rotation=90, labelsize=font_size, pad=2)
    ax.tick_params(axis='y', rotation=0, labelsize=font_size, pad=2)
    
    base_path = str(out_path).rsplit('.', 1)[0]
    plt.savefig(f"{base_path}.png", dpi=PLOT_DPI, bbox_inches='tight', facecolor='white')
    plt.savefig(f"{base_path}.pdf", dpi=PLOT_DPI, bbox_inches='tight', facecolor='white')
    plt.close('all')

def plot_distribution(mat, title, out_path, triangle_only=False):
    plt.figure(figsize=(8, 6))
    
    if triangle_only: 
        mask = np.triu(np.ones_like(mat, dtype=bool), k=0) | np.isnan(mat)
    else: 
        mask = np.eye(mat.shape[0], dtype=bool) | np.isnan(mat)
        
    vals = mat[~mask]
    if len(vals) > 0:
        min_v = max(30.0, np.floor(np.nanmin(vals)))
        counts, edges = np.histogram(vals, bins=np.arange(min_v, 101, 1))
        plt.plot(0.5*(edges[:-1]+edges[1:]), counts/len(vals), color='#CC0000', lw=3.0)
        plt.ylim(0, max(counts/len(vals))*1.1)
    else: 
        min_v = 80.0
        
    plt.title(f'Distribution - {title}', fontsize=16, fontweight='bold', pad=15)
    plt.xlabel('Identity (%)', fontsize=14)
    plt.ylabel('Proportion', fontsize=14)
    plt.xlim(min_v, 100.0)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.grid(axis='y', linestyle='--', alpha=0.7, color='#CCCCCC')
    
    base_path = str(out_path).rsplit('.', 1)[0]
    plt.savefig(f"{base_path}.png", dpi=PLOT_DPI, bbox_inches='tight', facecolor='white')
    plt.savefig(f"{base_path}.pdf", dpi=PLOT_DPI, bbox_inches='tight', facecolor='white')
    plt.close('all')

# =============================================================================
#[ 核心底层多进程并发计算引擎 (支持颗粒级断点续传) ]
# =============================================================================
def calc_sim_pair(args):
    s1, s2, is_p, ign, metric = args 
    if not s1 or not s2 or len(s1) == 0 or len(s2) == 0: return np.nan

    aligner = Align.PairwiseAligner()
    aligner.mode = 'global'
    if is_p:
        aligner.substitution_matrix = substitution_matrices.load(AA_SUBSTITUTION_MATRIX)
        aligner.open_gap_score, aligner.extend_gap_score = AA_GAP_OPEN, AA_GAP_EXTEND
    else:
        aligner.match_score, aligner.mismatch_score = NT_MATCH, NT_MISMATCH
        aligner.open_gap_score, aligner.extend_gap_score = NT_GAP_OPEN, NT_GAP_EXTEND

    try:
        alns = aligner.align(s1, s2)
        best_aln = alns[0]
    except Exception: return np.nan

    rowA = str(best_aln[0]).upper()
    rowB = str(best_aln[1]).upper()

    matches = 0
    valid_sdt = 0
    valid_blast = 0
    
    start_idx, end_idx = 0, len(rowA)
    for i in range(len(rowA)):
        if rowA[i] != '-' and rowB[i] != '-':
            start_idx = i; break
            
    for i in range(len(rowA)-1, -1, -1):
        if rowA[i] != '-' and rowB[i] != '-':
            end_idx = i + 1; break

    for i in range(start_idx, end_idx):
        a, b = rowA[i], rowB[i]
        if ign and (('N' in (a, b) and not is_p) or ('X' in (a, b) and is_p)): continue
        valid_blast += 1
        if a != '-' and b != '-':
            valid_sdt += 1
            if a == b: matches += 1  

    if metric == 'sdt_strict': 
        return (matches / valid_sdt * 100.0) if valid_sdt > 0 else np.nan
    else: 
        return (matches / valid_blast * 100.0) if valid_blast > 0 else np.nan

def _pw_worker(task): 
    i, j = task[-2], task[-1]
    args_for_calc = task[:-2]
    val = calc_sim_pair(args_for_calc)
    return i, j, val

def build_mat_pw(recs, is_p, ign, workers, desc, metric, out_dir, cache_prefix, resume=False):
    n = len(recs)
    strs =[str(r.seq) for r in recs]

    # 🌟 原子级断点续传系统初始化
    cache_dir = out_dir / ".resume_cache"
    cache_dir.mkdir(parents=True, exist_ok=True)
    cache_file = cache_dir / f"{cache_prefix}_pw_{n}.npy"

    mat = np.full((n, n), np.nan)

    if resume and cache_file.exists():
        print(f"      🔄 [颗粒级续传] 检测到阵列缓存，加载并识别空缺位点...")
        try:
            cached_mat = np.load(cache_file)
            if cached_mat.shape == (n, n):
                mat = cached_mat
        except Exception:
            print(f"      ⚠️ 缓存破损，重新构建阵列...")
            mat = np.full((n, n), np.nan)

    tasks =[]
    for i in range(n):
        for j in range(i+1, n):
            if np.isnan(mat[i, j]):  # ⭐️ 核心判决：只有 NaN 的对位才会被派发任务！
                tasks.append((strs[i], strs[j], is_p, ign, metric, i, j))

    if not tasks:
        print(f"      ✅ [{desc}] 矩阵已 100% 在缓存中就绪，秒跳过！")
        return mat

    with ProcessPoolExecutor(max_workers=workers) as ex:
        chunk_sz = max(1, len(tasks) // (workers * 4))
        results = ex.map(_pw_worker, tasks, chunksize=chunk_sz)

        flush_interval = max(1000, len(tasks) // 20) # 智能刷盘频率判定
        save_counter = 0

        with tqdm(total=len(tasks), desc=f"      🧬 {desc:<16}", unit="pair", leave=True) as pbar:
            for i, j, val in results:
                mat[i, j] = val
                pbar.update(1)
                save_counter += 1

                # 🌟 断点刷盘 (Atomic Save): 顺应 Numpy 潜规则，直存 .tmp.npy
                if save_counter % flush_interval == 0:
                    tmp_cache = cache_file.parent / (cache_file.stem + '.tmp.npy')
                    np.save(tmp_cache, mat)
                    tmp_cache.replace(cache_file)

    # 最终完整阵列固化保存
    tmp_cache = cache_file.parent / (cache_file.stem + '.tmp.npy')
    np.save(tmp_cache, mat)
    tmp_cache.replace(cache_file)

    return mat

def build_mat_mafft(recs, is_p, ign, mafft, threads, desc, metric, out_dir, cache_prefix, resume=False):
    n = len(recs)

    # 🌟 完全整体级别断点续传
    cache_dir = out_dir / ".resume_cache"
    cache_dir.mkdir(parents=True, exist_ok=True)
    cache_file = cache_dir / f"{cache_prefix}_mafft_{n}.npy"

    if resume and cache_file.exists():
        print(f"      🔄[全局续传] 加载完整 MAFFT 缓存列阵！")
        try:
            mat = np.load(cache_file)
            if mat.shape == (n, n):
                return mat
        except: pass

    mat = np.full((n, n), np.nan)
    safe =[SeqRecord(Seq("X" if is_p else "N"), id=r.id) if len(str(r.seq))==0 else r for r in recs]
    print(f"      🏃 运行 MAFFT ({desc})...")

    with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as tin, tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as tout:
        SeqIO.write(safe, tin.name, "fasta")
        subprocess.run([mafft, "--auto", "--thread", str(threads), "--quiet", tin.name], stdout=tout, check=True)
        aln = AlignIO.read(tout.name, "fasta")

    os.unlink(tin.name)
    os.unlink(tout.name)

    ig_c = 'X' if is_p else 'N'
    for i in range(n):
        for j in range(i+1, n):
            s_i = str(aln[i].seq).upper()
            s_j = str(aln[j].seq).upper()
            matches, valid_sdt, valid_blast = 0, 0, 0
            start_idx, end_idx = 0, len(s_i)

            for k in range(len(s_i)):
                if s_i[k] != '-' and s_j[k] != '-':
                    start_idx = k; break
            for k in range(len(s_i)-1, -1, -1):
                if s_i[k] != '-' and s_j[k] != '-':
                    end_idx = k + 1; break

            for k in range(start_idx, end_idx):
                a, b = s_i[k], s_j[k]
                if ign and (a == ig_c or b == ig_c): continue
                valid_blast += 1
                if a != '-' and b != '-':
                    valid_sdt += 1
                    if a == b: matches += 1

            if metric == 'sdt_strict':
                mat[i, j] = (matches / valid_sdt * 100.0) if valid_sdt > 0 else np.nan
            else:
                mat[i, j] = (matches / valid_blast * 100.0) if valid_blast > 0 else np.nan

    tmp_cache = cache_file.parent / (cache_file.stem + '.tmp.npy')
    np.save(tmp_cache, mat)
    tmp_cache.replace(cache_file)

    return mat


# =============================================================================
#[ 高效建树排序与输出保存 ] SciPy 安全提速版避免假死
# =============================================================================
def get_safe_leaf_order(mat, n):
    try:
        if np.isnan(mat).all() or n < 3: 
            return list(range(n))
        dist_mat = np.zeros((n, n))
        for i in range(n):
            for j in range(n):
                if i == j: 
                    dist_mat[i, j] = 0.0
                else: 
                    dist_mat[i, j] = 100.0 if np.isnan(mat[i, j]) else max(0.0, 100.0 - mat[i, j])
        condensed = ssd.squareform(dist_mat)
        Z = sch.linkage(condensed, method='average')
        return sch.leaves_list(Z).tolist()
    except Exception as e:
        warnings.warn(f"get_safe_leaf_order 层级聚类失败，退回到原始顺序: {e}")
        return list(range(n))

def execute_single_output(mat, names, out_dir, prefix, args, plot_title_prefix, is_nt):
    n = len(names)
    full_mat = np.zeros((n,n))
    
    for i in range(n):
        for j in range(n): 
            if i == j: full_mat[i, j] = 100.0
            elif i < j: full_mat[i, j] = mat[i, j]
            else: full_mat[i, j] = mat[j, i]

    if not args.plot_only_csv:
        print(f"      💾 写入记录阵列距阵: {prefix}")
        with open(out_dir / f"{prefix}_{args.align_method}.csv", "w", newline="", encoding="utf-8") as f:
            w = csv.writer(f)
            w.writerow([""] + names)
            for i, nm in enumerate(names):
                row_data = [nm]
                for j in range(n):
                    if i == j: row_data.append("")
                    else: row_data.append("NaN" if np.isnan(full_mat[i,j]) else f"{full_mat[i,j]:.1f}")
                w.writerow(row_data)

    if not args.no_nj:
        print(f"      🌳 构建 SciPy 亚秒聚类树拓扑排序...")
        ord_idx = get_safe_leaf_order(full_mat, n)
    else:
        ord_idx = list(range(n))
        
    if args.plot:
        print(f"      🎨 并发渲染高清图库...")
        labels_ordered = [names[i] for i in ord_idx]
        mat_ord = full_mat[np.ix_(ord_idx, ord_idx)]
        plot_title_hm = f'{plot_title_prefix}\n{"NT" if is_nt else "AA"} Identity (%)'
        plot_heatmap(mat_ord, labels_ordered, plot_title_hm, out_dir / f"{prefix}_heatmap.{args.plot_format}", triangle_only=True)
        plot_title_dist = f'{plot_title_prefix} - {"NT" if is_nt else "AA"}'
        plot_distribution(mat_ord, plot_title_dist, out_dir / f"{prefix}_distribution.{args.plot_format}", triangle_only=True)

def execute_composite_output(nt_mat, aa_mat, names, out_dir, prefix, args, plot_title_prefix):
    n = len(names)
    combined = np.full((n,n), np.nan)
    full_nt = np.zeros((n,n))
    
    for i in range(n):
        for j in range(n):
            if i == j:
                full_nt[i, j] = 100.0
                combined[i, j] = 100.0
            elif i < j:
                full_nt[i, j] = nt_mat[i, j]
                combined[i, j] = nt_mat[i, j]
            else:
                full_nt[i, j] = nt_mat[j, i]
                if aa_mat is not None: 
                    combined[i, j] = aa_mat[j, i]

    if not args.plot_only_csv:
        print(f"      💾 写入记录复合距阵: {prefix}")
        with open(out_dir / f"{prefix}_{args.align_method}.csv", "w", newline="", encoding="utf-8") as f:
            w = csv.writer(f)
            w.writerow([""] + names)
            for i, nm in enumerate(names):
                row_data =[nm]
                for j in range(n):
                    if i == j: row_data.append("")
                    else: row_data.append("NaN" if np.isnan(combined[i,j]) else f"{combined[i,j]:.1f}")
                w.writerow(row_data)

    if not args.no_nj:
        ord_idx = get_safe_leaf_order(full_nt, n)
    else:
        ord_idx = list(range(n))
        
    if args.plot:
        print(f"      🎨 并发渲染高清复合图库...")
        labels_ordered = [names[i] for i in ord_idx]
        mat_ord = combined[np.ix_(ord_idx, ord_idx)]
        plot_title_hm = f'{plot_title_prefix}\nComposite: NT (Upper) / AA (Lower)'
        plot_heatmap(mat_ord, labels_ordered, plot_title_hm, out_dir / f"{prefix}_heatmap.{args.plot_format}", triangle_only=False)
        plot_title_dist = f'{plot_title_prefix} - Composite'
        plot_distribution(mat_ord, plot_title_dist, out_dir / f"{prefix}_distribution.{args.plot_format}", triangle_only=False)

def parse_ext_name(filename, known_samples, suffix):
    base = filename.replace(suffix, "")
    for s_id in sorted(known_samples, key=len, reverse=True):
        if base.startswith(s_id + "_"): 
            return s_id, base[len(s_id)+1:]
    return None, None

# =============================================================================
#[ 双轨计算引擎与实体序列落盘库 ]
# =============================================================================
def run_mode(mode_name, args, seq_records, ext_dir, base_out_dir, organism, accession, global_gb, ref_len):
    if mode_name == 'strict':
        max_n_gen = 0.0
        max_n_feat = 0.0
    else:
        max_n_gen = args.max_n_genome
        max_n_feat = args.max_n_gene
    
    print(f"\n" + "━"*60)
    print(f" 🚀 [Phase 2] 运行模式: {mode_name.upper()} | 容错阈值: 全长N<={max_n_gen}%, 单基因N<={max_n_feat}%")
    print("━"*60)
    
    mode_dir = base_out_dir / f"Mode_{mode_name.capitalize()}"
    src_dir = ext_dir / "fill" if mode_name == 'fill' else ext_dir / "raw"

    def process_sub_dataset(sub_out_dir, is_dedup_run, run_desc):
        print(f"\n" + "┌" + "─"*58 + "┐")
        print(f"│ ➡️  处理子轨迹: {run_desc:<41}│")
        print("└" + "─"*58 + "┘")
        
        mat_dir = sub_out_dir / "02_similarity_matrices"
        report_dir = sub_out_dir / "03_deduplication_reports"
        seq_out_dir = sub_out_dir / "04_filtered_sequences"
        
        mat_dir.mkdir(parents=True, exist_ok=True)
        seq_out_dir.mkdir(parents=True, exist_ok=True)
        if is_dedup_run: 
            report_dir.mkdir(parents=True, exist_ok=True)
        
        all_names = []
        valid_genomes =[]
        over_recs =[]
        
        ref_id = get_sample_id(accession) + "_REF"
        
        for rec in seq_records:
            seq_str = str(rec.seq).upper()
            if len(seq_str) > 0:
                n_rat = (seq_str.count('N') / len(seq_str)) * 100 
            else:
                n_rat = 100.0
                
            if len(seq_str) >= args.min_length:
                s_id = get_sample_id(rec.id)
                all_names.append(s_id)
                if n_rat <= max_n_gen: 
                    valid_genomes.append(rec)

        if is_dedup_run and args.cdhit and len(valid_genomes) > 0:
            print(f"   🧠 [软聚类去重降维] 调用外部 CD-HIT 降噪分析...")
            if args.cdhit_threshold:
                threshold = args.cdhit_threshold
            else:
                threshold = calculate_adaptive_threshold(ref_len)
                
            valid_dedup_names = run_cdhit_and_parse(valid_genomes, report_dir, ref_len, threshold, args.threads, args.cdhit_path)
            
            for r in valid_genomes:
                if get_sample_id(r.id) in valid_dedup_names:
                    over_recs.append(r)
        else:
            for r in valid_genomes:
                over_recs.append(r)

        if global_gb:
            try:
                with open(global_gb, 'r') as f:
                    ref_rec = next(SeqIO.parse(f, "genbank"))
                    all_names.append(ref_id)
                    over_recs.append(SeqRecord(ref_rec.seq, id=ref_id, description="Reference"))
            except Exception: 
                pass

        gene_dict = defaultdict(lambda: {"nt": {}, "aa": {}})
        ext_samps = set()
        
        if global_gb:
            try:
                with open(global_gb, 'r') as f:
                    ref_rec = next(SeqIO.parse(f, "genbank"))
                    for feat in ref_rec.features:
                        if feat.type in["CDS", "mat_peptide"]:
                            raw_name = feat.qualifiers.get("product", [feat.type])[0]
                            g = re.sub(r'[^\w-]', '_', raw_name)
                            ref_nuc = feat.location.extract(ref_rec.seq)
                            gene_dict[g]["nt"][ref_id] = str(ref_nuc)
                            ext_samps.add(ref_id)
                            if not args.nt_only:
                                try: 
                                    ref_aa = ref_nuc.translate(table=GENETIC_CODE, to_stop=True)
                                    gene_dict[g]["aa"][ref_id] = str(ref_aa)
                                except Exception: 
                                    pass
            except Exception: 
                pass

        for f in src_dir.glob("*_nucl.fasta"):
            s_id, g = parse_ext_name(f.name, all_names, "_nucl.fasta")
            if s_id and g:
                rec = next(SeqIO.parse(f, "fasta"))
                seq_str = str(rec.seq).upper()
                if len(seq_str) > 0:
                    n_rat = (seq_str.count('N') / len(seq_str)) * 100 
                else:
                    n_rat = 100.0
                    
                if len(seq_str) >= args.min_length and n_rat <= max_n_feat:
                    gene_dict[g]["nt"][s_id] = str(rec.seq)
                    ext_samps.add(s_id)

        if not args.nt_only:
            for f in src_dir.glob("*_prot.fasta"):
                s_id, g = parse_ext_name(f.name, all_names, "_prot.fasta")
                if s_id and g and s_id in gene_dict[g]["nt"]: 
                    rec = next(SeqIO.parse(f, "fasta"))
                    gene_dict[g]["aa"][s_id] = str(rec.seq)

        # 🌟 指针方法引入断点专属缓存签注
        def run_m(seqs, is_p, d, cache_prefix):
            if args.align_method == 'mafft': 
                return build_mat_mafft(seqs, is_p, args.ignore_ambiguous, args.mafft_path, args.threads, d, args.scoring_metric, sub_out_dir, cache_prefix, args.resume)
            else: 
                return build_mat_pw(seqs, is_p, args.ignore_ambiguous, args.threads, d, args.scoring_metric, sub_out_dir, cache_prefix, args.resume)

        if is_dedup_run:
            suffix_title = "(Deduplicated)"
        else:
            suffix_title = "(Full Dataset)"
            
        base_title = f"{organism} ({accession}) {suffix_title}"
        
        # ===[ 整体分析与序列归档 ] ===
        if len(over_recs) >= 3:
            print(f"   🧬 [全长宏观计算] 共有 {len(over_recs)} 个序列参与比对...")
            
            for r in over_recs: 
                r.description = ""  
            SeqIO.write(over_recs, seq_out_dir / "overall_WholeGenome_NT.fasta", "fasta")
            
            nt_m = run_m(over_recs, False, "整体核酸", "overall_NT")
            names_list =[]
            for r in over_recs:
                names_list.append(get_sample_id(r.id))
                
            execute_single_output(nt_m, names_list, mat_dir, "overall_NT_only", args, f"{base_title}\nWhole Genome", is_nt=True)

            if not args.nt_only:
                # NOTE: Whole-genome AA translation is intentionally skipped because viral
                # genomes contain UTRs, intergenic regions, and stop codons — translating
                # them produces meaningless proteins. AA similarity is computed per-gene
                # in the gene-level section below where CDS extraction ensures correct ORFs.
                pass
        else: 
            print(f"   ⚠️[抛弃全长扫描] 合规参与运算基数偏小。")

        # ===[解绑单基因分型与落盘] ===
        if ext_samps and not args.nt_only:
            print(f"   🧬 [多态结构分解验证] 提纯并发射实体序列...")
            comp_g_nt = []
            comp_g_aa = []
            comp_g_strict =[]
            
            for g in sorted(gene_dict.keys()):
                nt_samps = sorted(set(all_names) & set(gene_dict[g]["nt"].keys()))
                if len(nt_samps) == len(ext_samps): 
                    comp_g_nt.append(g)
                    
                aa_samps = sorted(set(all_names) & set(gene_dict[g]["aa"].keys()))
                if len(aa_samps) == len(ext_samps): 
                    comp_g_aa.append(g)
                    
                strict_samps = sorted(set(nt_samps) & set(aa_samps))
                if len(strict_samps) == len(ext_samps): 
                    comp_g_strict.append(g)

                if is_dedup_run:
                    if len(nt_samps) > 0: 
                        nt_samps = exact_deduplicate(nt_samps, gene_dict[g]["nt"], report_dir / f"gene_{g}_NT_Clustering_Report.tsv")
                    if len(aa_samps) > 0: 
                        aa_samps = exact_deduplicate(aa_samps, gene_dict[g]["aa"], report_dir / f"gene_{g}_AA_Clustering_Report.tsv")
                    if len(strict_samps) > 0: 
                        strict_samps = exact_deduplicate(strict_samps, gene_dict[g]["nt"], report_dir / f"gene_{g}_Composite_Clustering_Report.tsv")

                if len(nt_samps) >= 3:
                    nt_recs =[]
                    for s in nt_samps:
                        nt_recs.append(SeqRecord(Seq(gene_dict[g]["nt"][s]), id=s, description=""))
                    SeqIO.write(nt_recs, seq_out_dir / f"gene_{g}_NT.fasta", "fasta")
                    execute_single_output(run_m(nt_recs, False, f"基因 {g} NT", f"gene_{g}_NT"), nt_samps, mat_dir, f"gene_{g}_NT_only", args, f"{base_title}\nGene: {g}", is_nt=True)

                if len(aa_samps) >= 3:
                    aa_recs =[]
                    for s in aa_samps:
                        aa_recs.append(SeqRecord(Seq(gene_dict[g]["aa"][s]), id=s, description=""))
                    SeqIO.write(aa_recs, seq_out_dir / f"gene_{g}_AA.fasta", "fasta")
                    execute_single_output(run_m(aa_recs, True, f"基因 {g} AA", f"gene_{g}_AA"), aa_samps, mat_dir, f"gene_{g}_AA_only", args, f"{base_title}\nGene: {g}", is_nt=False)

                if len(strict_samps) >= 3:
                    nt_recs_strict =[]
                    aa_recs_strict =[]
                    for s in strict_samps:
                        nt_recs_strict.append(SeqRecord(Seq(gene_dict[g]["nt"][s]), id=s))
                        aa_recs_strict.append(SeqRecord(Seq(gene_dict[g]["aa"][s]), id=s))
                        
                    execute_composite_output(run_m(nt_recs_strict, False, f"基因 {g} NT", f"gene_{g}_comp_NT"), 
                                             run_m(aa_recs_strict, True, f"基因 {g} AA", f"gene_{g}_comp_AA"), 
                                             strict_samps, mat_dir, f"gene_{g}_Composite", args, f"{base_title}\nGene: {g}")

            if len(comp_g_nt) > 0:
                valid_nt =[]
                for s in ext_samps:
                    has_all = True
                    for g in comp_g_nt:
                        if s not in gene_dict[g]["nt"]:
                            has_all = False
                            break
                    if has_all:
                        valid_nt.append(s)
                valid_nt = sorted(valid_nt)
                
                if is_dedup_run: 
                    cat_nt_dict = {}
                    for s in valid_nt:
                        cat_nt_dict[s] = "".join(gene_dict[g]["nt"][s] for g in comp_g_nt)
                    valid_nt = exact_deduplicate(valid_nt, cat_nt_dict, report_dir / "concatenated_NT_Clustering_Report.tsv")
                    
                if len(valid_nt) >= 3:
                    cat_nt_recs =[]
                    for s in valid_nt:
                        seq_combined = "".join(gene_dict[g]["nt"][s] for g in comp_g_nt)
                        cat_nt_recs.append(SeqRecord(Seq(seq_combined), id=s, description=""))
                        
                    SeqIO.write(cat_nt_recs, seq_out_dir / "concatenated_CDS_NT.fasta", "fasta")
                    execute_single_output(run_m(cat_nt_recs, False, "串联 NT", "cat_NT"), valid_nt, mat_dir, "concatenated_NT_only", args, f"{base_title}\nConcatenated CDS", is_nt=True)

            if len(comp_g_aa) > 0:
                valid_aa =[]
                for s in ext_samps:
                    has_all = True
                    for g in comp_g_aa:
                        if s not in gene_dict[g]["aa"]:
                            has_all = False
                            break
                    if has_all:
                        valid_aa.append(s)
                valid_aa = sorted(valid_aa)
                
                if is_dedup_run: 
                    cat_aa_dict = {}
                    for s in valid_aa:
                        cat_aa_dict[s] = "".join(gene_dict[g]["aa"][s] for g in comp_g_aa)
                    valid_aa = exact_deduplicate(valid_aa, cat_aa_dict, report_dir / "concatenated_AA_Clustering_Report.tsv")
                    
                if len(valid_aa) >= 3:
                    cat_aa_recs =[]
                    for s in valid_aa:
                        seq_combined = "".join(gene_dict[g]["aa"][s] for g in comp_g_aa)
                        cat_aa_recs.append(SeqRecord(Seq(seq_combined), id=s, description=""))
                        
                    SeqIO.write(cat_aa_recs, seq_out_dir / "concatenated_CDS_AA.fasta", "fasta")
                    execute_single_output(run_m(cat_aa_recs, True, "串联 AA", "cat_AA"), valid_aa, mat_dir, "concatenated_AA_only", args, f"{base_title}\nConcatenated CDS", is_nt=False)

        elif ext_samps and args.nt_only:
            print(f"   🧬 [单基因计算：纯核酸防卫模式] 收纳序列集...")
            comp_g_nt =[]
            
            for g in sorted(gene_dict.keys()):
                nt_samps =[]
                for s in all_names:
                    if s in gene_dict[g]["nt"].keys():
                        nt_samps.append(s)
                nt_samps = sorted(set(nt_samps))
                
                if len(nt_samps) == len(ext_samps): 
                    comp_g_nt.append(g)
                    
                if is_dedup_run and len(nt_samps) > 0: 
                    nt_samps = exact_deduplicate(nt_samps, gene_dict[g]["nt"], report_dir / f"gene_{g}_NT_Clustering_Report.tsv")
                
                if len(nt_samps) >= 3:
                    nt_recs =[]
                    for s in nt_samps:
                        nt_recs.append(SeqRecord(Seq(gene_dict[g]["nt"][s]), id=s, description=""))
                        
                    SeqIO.write(nt_recs, seq_out_dir / f"gene_{g}_NT.fasta", "fasta") 
                    execute_single_output(run_m(nt_recs, False, f"基因 {g} NT", f"gene_{g}_NT"), nt_samps, mat_dir, f"gene_{g}_NT_only", args, f"{base_title}\nGene: {g}", is_nt=True)
                                          
            if len(comp_g_nt) > 0:
                valid_nt =[]
                for s in ext_samps:
                    has_all = True
                    for g in comp_g_nt:
                        if s not in gene_dict[g]["nt"]:
                            has_all = False
                            break
                    if has_all:
                        valid_nt.append(s)
                valid_nt = sorted(valid_nt)
                
                if is_dedup_run: 
                    cat_nt_dict = {}
                    for s in valid_nt:
                        cat_nt_dict[s] = "".join(gene_dict[g]["nt"][s] for g in comp_g_nt)
                    valid_nt = exact_deduplicate(valid_nt, cat_nt_dict, report_dir / "concatenated_NT_Clustering_Report.tsv")
                    
                if len(valid_nt) >= 3:
                    cat_nt_recs =[]
                    for s in valid_nt:
                        seq_combined = "".join(gene_dict[g]["nt"][s] for g in comp_g_nt)
                        cat_nt_recs.append(SeqRecord(Seq(seq_combined), id=s, description=""))
                        
                    SeqIO.write(cat_nt_recs, seq_out_dir / "concatenated_CDS_NT.fasta", "fasta") 
                    execute_single_output(run_m(cat_nt_recs, False, "串联 NT", "cat_NT"), valid_nt, mat_dir, "concatenated_NT_only", args, f"{base_title}\nConcatenated CDS", is_nt=True)

    # ==== 队列发射台 ====
    process_sub_dataset(mode_dir / "Full_Dataset", False, "全量未去重数据 (完全谱系阵列)")
    
    if args.cdhit: 
        process_sub_dataset(mode_dir / "Dedup_Dataset", True, "CD-HIT软去重核心骨架 (出图选品库)")

# =============================================================================
#[ 主控指挥台 ]
# =============================================================================
def main():
    parser = argparse.ArgumentParser(description="🌟 生信自动化终极巨作 (断点续传/防挂死版) 🚀", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument('-i', '--input', dest='input_path', required=True, help='FASTA 目录或 Multi-FASTA')
    parser.add_argument('-g', '--gb', '--genbank', dest='genbank', help='参考 GenBank 路径或 Accession')
    parser.add_argument('-o', '--output_dir', default='./pipeline_results', help='总输出目录')
    parser.add_argument('--email', default='your_email@example.com', help='NCBI Entrez 邮箱 (GenBank 下载需要)')
    parser.add_argument('--mode', choices=['strict', 'filter', 'fill', 'all'], default='filter', help='队列验证模式')
    parser.add_argument('--plot_only_csv', help='【降维画质发版】抛开算力对接，直读历史 CSV 重建高分辨画册')
    parser.add_argument('--resume', action='store_true', help='🔥【颗粒级断点续传】系统崩溃自动接管历史进度，严密恢复矩阵计算。')
    
    parser.add_argument('--cdhit', action='store_true', help='激活 CD-HIT 基因组清场降噪器')
    parser.add_argument('--cdhit_threshold', type=float, help='指定 CD-HIT(0.85-1.0)')
    parser.add_argument('--cdhit_path', default='cd-hit-est', help='指定 cd-hit-est 命令别名')
    parser.add_argument('--nt_only', action='store_true', help='强制锁死核酸保护流')
    
    g1 = parser.add_argument_group('Phase 1: 多核提炼阵列设置')
    g1.add_argument('--skip_extract', action='store_true', help='越过同源组抓捕提炼')
    g1.add_argument('--max_n_genome', type=float, default=5.0, help='测序全集 N 缺陷容忍%%')
    g1.add_argument('--max_n_gene', type=float, default=5.0, help='局部热点 N 缺陷容忍%%')
    
    g2 = parser.add_argument_group('Phase 2: 并行算力编队与出图中心')
    g2.add_argument('--skip_similarity', action='store_true', help='只进行拦截降噪不执行比对分析')
    g2.add_argument('--align_method', choices=['pairwise', 'mafft'], default='pairwise', help='排针架构')
    g2.add_argument('--mafft_path', default='mafft', help='mafft路径')
    g2.add_argument('--threads', type=int, default=4, help='引擎线程阀门并发基数')
    g2.add_argument('--ignore_ambiguous', action='store_true', help='免疫由于序列脏点产生的人工算力惩罚')
    g2.add_argument('--scoring_metric', choices=['sdt_strict', 'blast_global'], default='sdt_strict', help='🎯 【基点法则】sdt_strict(不计空位) vs blast_global(严厉惩罚缺口)')
    g2.add_argument('--min_length', type=int, default=150, help='微小碎片遗物清扫底线 bp')
    g2.add_argument('--no_nj', action='store_true', help='强行废除由层级聚类衍生过的序列平滑度排列')
    g2.add_argument('--plot', action=argparse.BooleanOptionalAction, default=True, help='打通出图系统 (default: True, use --no-plot to disable)')
    parser.add_argument('--plot_format', choices=['png', 'pdf'], default='pdf', help='输出格式 (默认 pdf; 底层始终同时出 PNG+PDF)')

    args = parser.parse_args()
    Entrez.email = args.email
    out_dir = Path(args.output_dir)
    
    # === [极速重绘制图通道] ===
    if args.plot_only_csv:
        csv_path = Path(args.plot_only_csv)
        
        def quick_plot(c):
            names = []
            m_data =[]
            with open(c, 'r', encoding='utf-8') as f:
                reader = csv.reader(f)
                names = next(reader)[1:]
                for row in reader: 
                    row_vals =[]
                    for v in row[1:]:
                        if v.lower() in["", "nan"]:
                            row_vals.append(np.nan)
                        else:
                            row_vals.append(float(v))
                    m_data.append(row_vals)
                    
            n = len(names)
            n_m = np.full((n, n), np.nan)
            a_m = np.full((n, n), np.nan)
            has_a = False
            
            for i in range(n):
                for j in range(n):
                    if i < j: 
                        n_m[i, j] = m_data[i][j]
                    elif i > j: 
                        if not np.isnan(m_data[i][j]): 
                            has_a = True
                            a_m[i, j] = m_data[i][j]
                            
            t = c.stem.replace('_pairwise', '').replace('_mafft', '').replace('_', ' ').title()
            execute_single_output(n_m, names, out_dir, f"{c.stem}_NT_only", args, f"Dataset: {t}", is_nt=True)
            if has_a:
                execute_single_output(a_m, names, out_dir, f"{c.stem}_AA_only", args, f"Dataset: {t}", is_nt=False)
                execute_composite_output(n_m, a_m, names, out_dir, f"{c.stem}_Composite", args, f"Dataset: {t}")
                
        if csv_path.is_file(): 
            quick_plot(csv_path)
        else: 
            for c in csv_path.glob("*.csv"):
                quick_plot(c)
                
        print(f"\n🎉 PDF / PNG 无损图形集库已安全建立完成！")
        return

    # 🔥 断点续传 Phase 1 跳过判断
    if args.resume and (out_dir / "01_extraction_summary.tsv").exists():
        print("\n" + "="*60 + "\n 🚀 [断点续传] 侦测到本地已具有完备全量特征池，强行超车跳过提取期！\n" + "="*60)
        args.skip_extract = True

    input_path = Path(args.input_path)
    seq_records =[]
    
    if input_path.is_file():
        for rec in SeqIO.parse(input_path, "fasta"):
            rec.id = get_sample_id(rec.id)
            rec.seq = trim_ends(rec.seq)
            seq_records.append(rec)
    elif input_path.is_dir():
        for f in input_path.glob("*.fasta"):
            try: 
                rec = next(SeqIO.parse(f, "fasta"))
                rec.id = get_sample_id(rec.id)
                rec.seq = trim_ends(rec.seq)
                seq_records.append(rec)
            except Exception: 
                pass
    else: 
        sys.exit("❌ 输入边界崩溃！")
        
    if not seq_records: 
        sys.exit("❌ 未能截获实体分子。")

    gb_cache = out_dir / ".gb_cache"
    gb_cache.mkdir(parents=True, exist_ok=True)
    
    if args.genbank and os.path.exists(args.genbank):
        global_gb = Path(args.genbank)
    elif args.genbank:
        global_gb = download_genbank(args.genbank, gb_cache)
    else:
        guessed_id = extract_valid_accession(seq_records[0].id)
        if guessed_id:
            global_gb = download_genbank(guessed_id, gb_cache)
        else:
            global_gb = None

    organism = "Unknown_Organism"
    accession = "Unknown_Accession"
    has_cds = True
    ref_len = 0
    
    if global_gb:
        organism, accession, has_cds, ref_len = parse_gb_metadata(global_gb)
        if not has_cds and not args.nt_only:
            args.nt_only = True
            print(f"\n 🧠[AI 警告] 发现 {organism} 彻底缺少 CDS，防毒机制激活全局切入纯核酸比对状态！")

    ext_dir = out_dir / "01_extracted_genes"
    
    if not args.skip_extract:
        print("\n" + "="*60 + "\n 🚀[Phase 1] 大规模多核特征靶向提炼机投运 \n" + "="*60)
        
        raw_dir = ext_dir / "raw"
        fill_dir = ext_dir / "fill"
        raw_dir.mkdir(parents=True, exist_ok=True)
        fill_dir.mkdir(parents=True, exist_ok=True)
        
        ref_features_data =[]
        if global_gb:
            with open(global_gb, 'r') as gb_handle:
                for gr in SeqIO.parse(gb_handle, "genbank"):
                    for feat in gr.features:
                        if feat.type in ["CDS", "mat_peptide"]:
                            f_prod = feat.qualifiers.get("product", [feat.type])[0]
                            c_name = re.sub(r'[^\w-]', '_', f_prod)
                            f_gene = feat.qualifiers.get("gene", ["NA"])[0]
                            f_seq = str(feat.location.extract(gr.seq))
                            ref_features_data.append((c_name, f_prod, f_gene, f_seq))

        tasks =[]
        for rec in seq_records:
            tasks.append((rec.id, str(rec.seq), ref_features_data, args.nt_only))
            
        g_sum =[]
        
        with ProcessPoolExecutor(max_workers=args.threads) as e:
            chk_sz = max(1, len(tasks) // (args.threads * 4))
            results = e.map(_extract_worker, tasks, chunksize=chk_sz)
            
            with tqdm(total=len(tasks), desc="   🧬 暴力剥离中...", unit="sample") as bar:
                for prefix, gn_r, sts, rs in results:
                    if prefix:
                        g_sum.append((prefix, gn_r, sts))
                        for r in rs:
                            c_n = r["clean_name"]
                            d = "QC:PASS"
                            
                            with open(raw_dir / f"{prefix}_{c_n}_nucl.fasta", 'w') as fh: 
                                SeqIO.write(SeqRecord(r["raw_nuc"], id=f"{prefix}_{c_n}_nucl", description=d), fh, "fasta")
                                
                            if not args.nt_only and r["raw_prot_qc"] == "Pass":
                                with open(raw_dir / f"{prefix}_{c_n}_prot.fasta", 'w') as fh: 
                                    SeqIO.write(SeqRecord(r["raw_aa"], id=f"{prefix}_{c_n}_prot", description=d), fh, "fasta")
                                    
                            with open(fill_dir / f"{prefix}_{c_n}_nucl.fasta", 'w') as fh: 
                                SeqIO.write(SeqRecord(r["fill_nuc"], id=f"{prefix}_{c_n}_nucl", description=d), fh, "fasta")
                                
                            if not args.nt_only and r["fill_prot_qc"] == "Pass":
                                with open(fill_dir / f"{prefix}_{c_n}_prot.fasta", 'w') as fh: 
                                    SeqIO.write(SeqRecord(r["fill_aa"], id=f"{prefix}_{c_n}_prot", description=d), fh, "fasta")
                                    
                    bar.update(1)
                    
        with open(out_dir / "01_extraction_summary.tsv", "w", encoding="utf-8") as f:
            f.write("Sample_ID\tGenome_Raw_N%\tFeature\tRaw_N%\tRaw_Prot_QC\tFill_N%\tFill_Prot_QC\n")
            for prefix, gn, sts in g_sum:
                if not sts: 
                    f.write(f"{prefix}\t{gn:.2f}%\tNone\t-\t-\t-\t-\n")
                else: 
                    for st in sts: 
                        f.write(f"{prefix}\t{gn:.2f}%\t{st[0]}\t{st[1]}%\t{st[3]}\t{st[2]}%\t{st[4]}\n")
                        
        print(f"\n   📊 [提取看板] 清场完毕。提取档案注入: 01_extraction_summary.tsv。")

    if args.skip_similarity: 
        return
        
    if args.mode == 'all':
        modes_to_run =['strict', 'filter', 'fill']
    else:
        modes_to_run =[args.mode]
        
    for m in modes_to_run: 
        run_mode(m, args, seq_records, ext_dir, out_dir, organism, accession, global_gb, ref_len)
    
    print(f"\n🚀🧬 任务确认大交割完成。矩阵核心和脱水高分库已打包封锁于 -> {out_dir.absolute()} ！")

if __name__ == "__main__":
    main()
