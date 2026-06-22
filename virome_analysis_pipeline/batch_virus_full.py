#!/usr/bin/env python3
"""
batch_virus_full.py — 桥接 VAP 产出与 OmniVirusAssembler 的自动化批量执行引擎
特性: 
1. 严格规范化 Taxonomy 命名，消除空格、符号及括号。
2. 按 [Taxonomy_Accession]/[Sample_Accession]/ 结构进行目录分发归档。
3. 严格的双端 Reads 识别正则匹配 (修复 SRR 样本名冲突问题)。
4. 支持基于 Covered% 过滤无效或极低丰度靶标。
"""

import argparse
import os
import re
import sys
import subprocess
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
import pandas as pd
from tqdm import tqdm

def safe_name(s: str, max_len: int = 100) -> str:
    s = str(s)
    s = re.sub(r'[^A-Za-z0-9\-.]', '_', s)
    s = re.sub(r'_+', '_', s)
    return s.strip('_.')[:max_len]

def normalize_taxonomy(tax_str: str) -> str:
    if pd.isna(tax_str) or not str(tax_str).strip():
        return "Unknown_Virus"
    s = str(tax_str).strip()
    s = re.sub(r'[^A-Za-z0-9]', '_', s)
    s = re.sub(r'_+', '_', s)
    return s.strip('_')

def find_raw_reads(sample_name, clean_data_dir):
    """严格正则匹配原始 Clean Data"""
    valid_exts = ['.fastq.gz', '.fq.gz', '.fastq', '.fq', '.fasta.gz', '.fa.gz', '.fasta', '.fa']
    matched_files = []
    
    s_clean = sample_name.strip().lower()
    for root, _, files in os.walk(clean_data_dir):
        for f in files:
            if not any(f.lower().endswith(ext) for ext in valid_exts): continue
            
            f_spaced = f.lower().replace('_', ' ').replace('.', ' ').replace('-', ' ')
            s_spaced = s_clean.replace('_', ' ').replace('.', ' ').replace('-', ' ')
            
            if re.search(r'\b' + re.escape(s_spaced) + r'\b', f_spaced) or f.lower().startswith(s_clean + "_") or f.lower().startswith(s_clean + "."):
                matched_files.append(os.path.join(root, f))
                
    r1, r2 = None, None
    for f in set(matched_files):
        fname = os.path.basename(f)
        # 强制包含分隔符，避免样本名冲突
        if re.search(r'(_R2|_2\.|_2_|\.R2|\.2\.)', fname, re.IGNORECASE):
            r2 = f
        elif re.search(r'(_R1|_1\.|_1_|\.R1|\.1\.)', fname, re.IGNORECASE):
            r1 = f
            
    reads_list = []
    if r1: reads_list.append(r1)
    if r2: reads_list.append(r2)
    # 兜底单端
    if not reads_list and matched_files:
        reads_list.append(list(set(matched_files))[0])
    return reads_list

def find_extracted_reads(sample, virus, reads_dir):
    """严格正则匹配靶向提取的 reads，屏蔽 SRR 字符干扰"""
    safe_v = safe_name(virus)
    prefix = f"{sample}.{safe_v}"
    
    matched = []
    for f in os.listdir(reads_dir):
        if f.startswith(prefix) and any(ext in f for ext in ['.fastq', '.fq', '.fasta', '.fa']):
            matched.append(os.path.join(reads_dir, f))
            
    r1, r2, rs = None, None, None
    for f in matched:
        fname = os.path.basename(f)
        # 严格限定前缀下划线或点，SRR2344 不会被误判为 R2
        if re.search(r'(_R1|_1\.|_1_|\.R1|\.1\.)', fname, re.IGNORECASE):
            r1 = f
        elif re.search(r'(_R2|_2\.|_2_|\.R2|\.2\.)', fname, re.IGNORECASE):
            r2 = f
        elif re.search(r'(single|se)', fname, re.IGNORECASE):
            rs = f
            
    res = []
    if r1: res.append(r1)
    if r2: res.append(r2)
    if rs and not (r1 and r2): res.append(rs)
    if not res and matched:
        # 如果格式太奇怪都没匹配上，原样返回防止崩溃
        res = sorted(matched)
    return res

def worker_run_virus_full(task):
    sample, virus, norm_tax, ref_fasta, ext_reads, raw_reads, out_dir, args = task
    
    log_file = out_dir.parent / f"{out_dir.name}_assembly.log"
    
    cmd = [
        "python3", args.virus_full_script,
        "-r", str(ref_fasta),
        "-o", str(out_dir),
        "-t", str(args.threads),
        "-j", "1",  
        "--assembly_tools", args.assembly_tools
    ]
    
    cmd.extend(["--assembly_reads"] + ext_reads)
    
    if raw_reads:
        cmd.extend(["--pvga_reads"] + raw_reads)
        cmd.extend(["--consensus_reads"] + raw_reads)
        
    if args.gb and os.path.exists(args.gb):
        cmd.extend(["--gb", args.gb])
        
    if args.extra_args:
        cmd.extend(args.extra_args.split())

    try:
        with open(log_file, "w") as log:
            log.write(f"CMD: {' '.join(cmd)}\n\n")
            proc = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, text=True)
            if proc.returncode == 0:
                return f"✅ [{sample}] {norm_tax} ({virus}) 组装完成"
            else:
                return f"❌ [{sample}] {norm_tax} ({virus}) 组装失败 (查看日志: {log_file})"
    except Exception as e:
        return f"❌ [{sample}] {norm_tax} 发生系统级错误: {e}"

def main():
    parser = argparse.ArgumentParser(description="批量调度 virus-full.py 进行全长组装打磨 (分类与样本重组归档版)")
    
    parser.add_argument("--downstream_dir", required=True, help="batch_virus_downstream.py 输出目录")
    parser.add_argument("--summary", required=True, help="确诊 summary 表格")
    parser.add_argument("--clean_data", required=True, help="原始 Clean Data 目录")
    parser.add_argument("--virus_full_script", default="virus-full.py", help="virus-full.py 路径")
    
    parser.add_argument("--outdir", default="./virus_assemblies_final", help="输出总目录")
    parser.add_argument("--min_covered", type=float, default=0.0, help="过滤 Covered%% 小于或等于该值的记录 (默认: 0.0)")
    parser.add_argument("--assembly_tools", default="all", help="启用的组装工具")
    parser.add_argument("--gb", help="全局 GenBank (.gb) 文件")
    
    parser.add_argument("--jobs", type=int, default=4, help="同时并行处理任务数")
    parser.add_argument("--threads", type=int, default=8, help="单任务线程数")
    parser.add_argument("--extra_args", type=str, default="", help="传递给 virus-full.py 的其它参数")
    
    args = parser.parse_args()
    
    downstream = Path(args.downstream_dir)
    if not downstream.exists():
        sys.exit(f"❌ 找不到 Downstream 目录: {downstream}")
        
    df = pd.read_csv(args.summary, sep='\t')
    
    if 'Sample' not in df.columns:
        sys.exit("❌ Summary 表格缺失 'Sample' 列。")
        
    acc_col = 'Accession' if 'Accession' in df.columns else 'Rep_Accession' if 'Rep_Accession' in df.columns else None
    if not acc_col:
        sys.exit("❌ Summary 表格缺失 'Accession' 列。")
        
    out_base = Path(args.outdir)
    out_base.mkdir(parents=True, exist_ok=True)
    
    tasks = []
    skipped_cov = 0
    print(f"🔍 正在解析任务并建立规范化目录结构...")
    
    for _, row in df.iterrows():
        if str(row['Sample']).strip() == 'Taxonomy':
            continue
            
        sample = str(row['Sample']).strip()
        virus = str(row[acc_col]).strip()
        
        cov_val = 100.0  
        if 'Covered%' in row:
            try:
                cov_val = float(str(row['Covered%']).replace('%', ''))
            except ValueError:
                cov_val = 0.0
                
        if cov_val <= args.min_covered:
            skipped_cov += 1
            continue

        tax_raw = 'Unknown_Virus'
        for col in ['Taxonomy', 'Adjusted_Species', 'Species_NCBI', 'Species']:
            if col in df.columns and pd.notna(row.get(col)):
                tax_raw = row[col]; break
        norm_tax = normalize_taxonomy(tax_raw)
        
        safe_v = safe_name(virus)
        safe_s = safe_name(sample)
        
        ref_fasta = downstream / "virus-fasta" / f"ref_{safe_v}" / f"ref_{safe_v}.ref.fasta"
        if not ref_fasta.exists():
            continue
            
        ext_reads = find_extracted_reads(sample, virus, str(downstream / "virus_reads"))
        if not ext_reads:
            continue
            
        raw_reads = find_raw_reads(sample, args.clean_data)
        
        dir_L1 = f"{norm_tax}_{safe_v}"
        dir_L2 = f"{safe_s}_{safe_v}"
        
        task_outdir = out_base / dir_L1 / dir_L2
        task_outdir.mkdir(parents=True, exist_ok=True)
        
        tasks.append((sample, virus, norm_tax, ref_fasta, ext_reads, raw_reads, task_outdir, args))

    print(f"ℹ️  根据 --min_covered {args.min_covered} 的设定，过滤了 {skipped_cov} 个低覆盖度任务。")
    
    if not tasks:
        sys.exit("❌ 没有生成任何有效的组装任务，请检查过滤阈值和文件路径。")

    print(f"🚀 构建了 {len(tasks)} 个病毒全长组装任务！开始并行执行 (并发数: {args.jobs})...")
    
    success_count = 0
    print(f"[Assembly] 0/{len(tasks)} done", flush=True)
    with ProcessPoolExecutor(max_workers=args.jobs) as executor:
        futures = [executor.submit(worker_run_virus_full, t) for t in tasks]
        done = 0
        for future in tqdm(as_completed(futures), total=len(tasks), desc="Assembly", file=sys.stdout, ncols=80, mininterval=5):
            res = future.result()
            done += 1
            if "✅" in res:
                success_count += 1
            if done % 10 == 0 or done == len(tasks):
                print(f"[Assembly] {done}/{len(tasks)} ({success_count} OK)", flush=True)
            
    print("-" * 60)
    print(f"🎉 批量组装完毕！成功: {success_count} / 总计: {len(tasks)}")
    print(f"📁 完美的分类结果目录已保存在: {out_base.absolute()}")

if __name__ == "__main__":
    main()
