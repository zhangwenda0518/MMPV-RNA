#!/usr/bin/env python3
"""
virus_batch_qc_extract.py (v9 TSV 汇总输出版)
更新：将汇总表输出格式改为 TSV (\t 分隔符)，更符合生信分析习惯。
"""

import argparse
import sys
import os
import re
import time
from pathlib import Path
from Bio import SeqIO, Entrez
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

Entrez.email = os.environ.get("NCBI_EMAIL", "your_email@example.com")

def calculate_assembly_stats(sequences_dict):
    lengths = [len(str(seq.seq)) for seq in sequences_dict.values()]
    total_len = sum(lengths)
    total_n = sum(str(seq.seq).upper().count('N') for seq in sequences_dict.values())
    total_gc = sum(str(seq.seq).upper().count('G') + str(seq.seq).upper().count('C') for seq in sequences_dict.values())

    lengths.sort(reverse=True)
    n50, running_sum = 0, 0
    for l in lengths:
        running_sum += l
        if running_sum >= total_len / 2.0:
            n50 = l
            break

    gc_ratio = (total_gc / total_len * 100) if total_len > 0 else 0
    n_ratio = (total_n / total_len * 100) if total_len > 0 else 0
    return len(lengths), total_len, n50, gc_ratio, n_ratio

def check_protein_integrity(aa_seq):
    clean_seq = str(aa_seq).rstrip('*')
    stop_count = clean_seq.count('*')
    return stop_count == 0, stop_count

def download_genbank(accession, out_dir):
    if re.match(r'^[SCE]RR', accession): return None
    out_path = Path(out_dir) / f"{accession}.gb"
    if out_path.exists(): return out_path
    
    print(f"\n 🌐 正在从 NCBI 联网下载参考序列: {accession} ...")
    try:
        with Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="text") as handle:
            with open(out_path, "w") as f: f.write(handle.read())
        return out_path
    except Exception as e:
        print(f" [❌ ERR] 下载失败: {e}")
        return None

def extract_valid_accession(text):
    patterns = [r'(?<![A-Z0-9])([A-Z]{1,2}\d{5,6}(\.\d+)?)(?![A-Z0-9])', r'(?<![A-Z0-9])([A-Z]{2}_\d{6,8}(\.\d+)?)(?![A-Z0-9])']
    for p in patterns:
        for m in re.finditer(p, text):
            acc = m.group(1)
            if not re.match(r'^[SCE]RR', acc): return acc
    return None

def process_single_sample(fasta_path, gb_path, out_dir, fill_n=False):
    prefix = fasta_path.name.replace(".consensus.fasta", "").replace(".fasta", "")
    
    print(f"\n" + "━"*60)
    print(f" 📂 当前处理样本: {prefix}")
    print("━"*60)

    try:
        consensus_seqs = SeqIO.to_dict(SeqIO.parse(fasta_path, "fasta"))
    except: return None

    # 初始化样本统计
    sample_stats = {"total": 0, "pass": 0, "fail": 0, "failed_genes": []}

    contigs, total_bp, n50, gc, n_pct = calculate_assembly_stats(consensus_seqs)
    print(f" 📊 全局 N 比例: {n_pct:.2f}% | 序列长度: {total_bp} bp")
    print(f" {'状态':<6} | {'特征名称':<20} | {'原N碱基%':<8} | {'填补数' if fill_n else '内部终止符'}")
    print("-" * 55)

    with open(gb_path, 'r') as gb_handle:
        for gb_record in SeqIO.parse(gb_handle, "genbank"):
            target_contig = list(consensus_seqs.values())[0] if len(consensus_seqs) == 1 else consensus_seqs.get(gb_record.id)
            if not target_contig: continue

            for feature in gb_record.features:
                if feature.type in ["CDS", "mat_peptide"]:
                    raw_name = feature.qualifiers.get("product", [feature.type])[0]
                    clean_name = re.sub(r'[^\w-]', '_', raw_name)
                    gene_name = feature.qualifiers.get("gene", ["NA"])[0]

                    try:
                        nuc_seq = feature.location.extract(target_contig).seq
                        ref_nuc_seq = feature.location.extract(gb_record.seq)
                    except: continue

                    orig_n_count = str(nuc_seq).upper().count('N')
                    orig_n_ratio = (orig_n_count / len(nuc_seq)) * 100
                    filled_count = 0

                    if fill_n and orig_n_count > 0:
                        seq_chars = list(str(nuc_seq).upper())
                        ref_chars = list(str(ref_nuc_seq).upper())
                        for i in range(len(seq_chars)):
                            if seq_chars[i] == 'N' and i < len(ref_chars):
                                seq_chars[i] = ref_chars[i]
                                filled_count += 1
                        nuc_seq = Seq("".join(seq_chars))

                    current_n_count = str(nuc_seq).upper().count('N')
                    current_n_ratio = (current_n_count / len(nuc_seq)) * 100

                    try:
                        aa_seq = nuc_seq.translate(table=1, cds=False)
                    except: 
                        sample_stats["fail"] += 1
                        sample_stats["failed_genes"].append(raw_name)
                        continue

                    is_intact, stop_count = check_protein_integrity(aa_seq)

                    status_icon = "✅ PASS"
                    if not is_intact or current_n_ratio >= 5.0:
                        status_icon = "❌ FAIL"
                        sample_stats["fail"] += 1
                        sample_stats["failed_genes"].append(raw_name)
                    else:
                        sample_stats["pass"] += 1
                    
                    sample_stats["total"] += 1

                    if fill_n:
                        print(f" [{status_icon:^6}] | {raw_name[:20]:<20} | {orig_n_ratio:>5.2f}%   | 🔧 +{filled_count} bp")
                    else:
                        print(f" [{status_icon:^6}] | {raw_name[:20]:<20} | {orig_n_ratio:>5.2f}%   | {stop_count}")

                    fill_info = f"Filled:{filled_count}bp" if fill_n else ""
                    desc = f"Gene:{gene_name} N_ratio:{current_n_ratio:.2f}% Stops:{stop_count} QC:{status_icon[-4:]} {fill_info}".strip()

                    prot_id = f"{prefix}_{clean_name}_prot"
                    with open(out_dir / f"{prot_id}.fasta", 'w') as f:
                        SeqIO.write(SeqRecord(aa_seq, id=prot_id, description=desc), f, "fasta")

                    nucl_id = f"{prefix}_{clean_name}_nucl"
                    with open(out_dir / f"{nucl_id}.fasta", 'w') as f:
                        SeqIO.write(SeqRecord(nuc_seq, id=nucl_id, description=desc), f, "fasta")

    return prefix, sample_stats

def main():
    parser = argparse.ArgumentParser(description='病毒序列提取工具 (支持独立指定 TSV 汇总报告路径)')
    parser.add_argument('-i', '--input', required=True, help='输入目录或单个 FASTA 文件')
    parser.add_argument('-g', '--genbank', help='指定 GB 文件或 Accession (不填自动解析)')
    parser.add_argument('-o', '--out', default='extracted_features', help='输出序列 FASTA 的目录')
    # 修改点 1: 帮助文档提示更改为 TSV
    parser.add_argument('-s', '--summary', help='指定 TSV 汇总表的完整输出路径和文件名 (例: reports/my_summary.tsv)')
    parser.add_argument('--fill_n', action='store_true', help='开启此选项，将使用参考基因组序列填补共识序列中的 N 碱基')
    args = parser.parse_args()

    input_path = Path(args.input)
    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)
    gb_cache = out_dir / "gb_cache"
    gb_cache.mkdir(exist_ok=True)

    files = [input_path] if input_path.is_file() else list(input_path.glob("*.fasta"))
    if not files: sys.exit("❌ 未找到任何 .fasta 文件")

    print(f"\n🚀 发现 {len(files)} 个序列文件，准备开始任务...")
    if args.fill_n: print("⚠️  注意: 已开启 [--fill_n] 模式，将使用参考序列修补未知区域。")

    global_gb = None
    if args.genbank:
        global_gb = Path(args.genbank) if os.path.exists(args.genbank) else download_genbank(args.genbank, gb_cache)

    global_summary = {}

    for f in files:
        current_gb = global_gb or download_genbank(extract_valid_accession(f.name) or extract_valid_accession(str(f.absolute())), gb_cache)
        if current_gb:
            res = process_single_sample(f, current_gb, out_dir, args.fill_n)
            if res:
                prefix, stats = res
                global_summary[prefix] = stats

    print("\n\n" + "█"*60)
    print(" 📑 全局任务执行汇总报告")
    print("█"*60)
    total_processed = len(global_summary)
    perfect_samples = sum(1 for s in global_summary.values() if s["fail"] == 0 and s["total"] > 0)
    
    print(f"处理样本总数 : {total_processed}")
    print(f"全基因通过数 : {perfect_samples} (无质控失败特征)")
    print("-" * 60)
    print(f"{'样本前缀':<25} | {'总提取':<6} | {'通过':<6} | {'失败及原因'}")
    print("-" * 60)
    
    for prefix, stats in global_summary.items():
        pass_rate = "🟢" if stats["fail"] == 0 else ("🟡" if stats["pass"] > 0 else "🔴")
        failed_str = ",".join(stats["failed_genes"])[:20] + ("..." if len(",".join(stats["failed_genes"])) > 20 else "")
        failed_info = failed_str if stats["fail"] > 0 else "-"
        short_prefix = prefix[:23] + ".." if len(prefix) > 25 else prefix
        print(f"{pass_rate} {short_prefix:<23} | {stats['total']:<6} | {stats['pass']:<6} | {failed_info}")
    
    print("█"*60 + "\n")

    # --- 修改点 2: 确定 TSV 文件的输出路径 ---
    if args.summary:
        tsv_file_path = Path(args.summary)
        tsv_file_path.parent.mkdir(parents=True, exist_ok=True)
    else:
        # 修改点 3: 默认后缀改为 .tsv
        tsv_file_path = out_dir / "extraction_summary.tsv"

    # --- 修改点 4: 使用 \t 分隔符写入文件 ---
    try:
        with open(tsv_file_path, 'w', encoding='utf-8') as tsv_file:
            # 写入表头 (\t 分隔)
            tsv_file.write("Sample_ID\tTotal_Features\tPassed_QC\tFailed_QC\tFailed_Genes_List\n")
            # 写入每一行数据 (\t 分隔)
            for prefix, stats in global_summary.items():
                failed_full_list = "; ".join(stats["failed_genes"]) if stats["fail"] > 0 else "None"
                tsv_file.write(f"{prefix}\t{stats['total']}\t{stats['pass']}\t{stats['fail']}\t{failed_full_list}\n")
        print(f" 🎉 汇总表格已成功保存至: {tsv_file_path.absolute()}")
    except Exception as e:
        print(f" [❌ ERR] 写入汇总文件失败: {e}")

if __name__ == "__main__":
    main()
