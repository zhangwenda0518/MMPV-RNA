#!/usr/bin/env python3
import os
import argparse
import gzip
import glob
import sys
import time
import re

# ---------------------------------------------------------
# 第三方模块防爆引入 
# ---------------------------------------------------------
try:
    from tqdm import tqdm
except ImportError:
    def tqdm(iterable, *args, **kwargs):
        return iterable

# ---------------------------------------------------------
# 模块：基因组版本侦测器 (Auto-Build Detection - Regex Edition)
# ---------------------------------------------------------
def detect_build_context(vcf_path):
    acc_pattern = re.compile(r'(?:_|^|-)([a-zA-Z]{1,2}_?\d{5,8}\.\d+)')
    
    base_name = os.path.basename(vcf_path)
    parent_dir = os.path.basename(os.path.dirname(os.path.abspath(vcf_path)))

    match = acc_pattern.search(base_name)
    if match: return match.group(1).upper()
        
    match_dir = acc_pattern.search(parent_dir)
    if match_dir: return match_dir.group(1).upper()

    try:
        f = gzip.open(vcf_path, 'rt') if vcf_path.endswith('.gz') else open(vcf_path, 'r')
        for line in f:
            if not line.startswith('#'): break 
            if line.startswith('##reference='):
                ref = line.strip().split('=', 1)[1]
                ref = ref.replace('file://', '')
                base_ref = os.path.basename(ref)
                for ext in['.fasta', '.fa', '.fna', '.mmi', '.gz']:
                    base_ref = base_ref.replace(ext, '')
                
                match_ref = acc_pattern.search(base_ref)
                f.close()
                return match_ref.group(1).upper() if match_ref else base_ref
        f.close()
    except Exception:
        pass

    cleaned_base = base_name.replace('.ann.vcf', '').replace('.vcf.gz', '').replace('.vcf', '')
    if '_' in cleaned_base:
        return cleaned_base.split('_', 1)[-1]
    return cleaned_base

def count_vcf_lines(vcf_path):
    total = 0
    header = 0
    try:
        f = gzip.open(vcf_path, 'rt', encoding='utf-8') if vcf_path.endswith('.gz') else open(vcf_path, 'r', encoding='utf-8')
        for line in f:
            total += 1
            if line.startswith('#'): header += 1
        f.close()
    except Exception:
        pass
    return total, header

# ---------------------------------------------------------
# 核心模块：VCF 解析器与漏斗统计 (多重 DP/AF 雷达引擎)
# ---------------------------------------------------------
def parse_vcf(snpeff_file, output_file=None, min_dp=0, min_af=0.0, build="auto", filter_pass=True, quiet=False):
    if output_file is None:
        outdir = os.path.dirname(os.path.abspath(snpeff_file))
        base = os.path.basename(snpeff_file).replace('.ann.vcf.gz', '').replace('.ann.vcf', '').replace('.vcf.gz', '').replace('.vcf', '')
        output_file = os.path.join(outdir, f"{base}.maf")

    tumor_sample_barcode = os.path.basename(snpeff_file).replace('.ann.vcf.gz', '').replace('.ann.vcf', '').replace('.vcf.gz', '').replace('.vcf', '')

    resolved_build = build
    if str(build).lower() == "auto":
        resolved_build = detect_build_context(snpeff_file)

    total_lines, header_lines = count_vcf_lines(snpeff_file)
    
    stats = {
        "Total_Input_Variants": 0,
        "Dropped_Empty_ALT": 0,
        "Dropped_Not_PASS": 0,
        "Dropped_Low_DP": 0,
        "Dropped_Low_AF": 0,
        "Final_Retained_MAF": 0
    }

    try:
        infile = gzip.open(snpeff_file, 'rt', encoding='utf-8') if snpeff_file.endswith('.gz') else open(snpeff_file, 'r', encoding='utf-8')
    except Exception as e:
        if not quiet: print(f"[ERROR] Failed to open {snpeff_file}: {e}", file=sys.stderr)
        return False, None, resolved_build

    start_time = time.time()
    
    try:
        with open(output_file, 'w', encoding='utf-8') as outfile:
            header =[
                "Hugo_Symbol", "Entrez_Gene_Id", "Center", "NCBI_Build", "Chromosome", 
                "Start_Position", "End_Position", "Strand", "Variant_Classification", 
                "Variant_Type", "Reference_Allele", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2", 
                "Tumor_Sample_Barcode", "Protein_Change", "i_TumorVAF_WU", "i_transcript_name"
            ]
            outfile.write("\t".join(header) + "\n")

            pbar = tqdm(infile, total=total_lines, desc=f"Mapping {tumor_sample_barcode[:15]}", leave=False, disable=quiet)

            for line in pbar:
                if line.startswith('#'): continue
                
                stats["Total_Input_Variants"] += 1
                fields = line.rstrip('\n').split('\t')
                
                if len(fields) < 8: continue  

                filter_status = fields[6]
                alt_allele = fields[4]

                if alt_allele == '.' or alt_allele == "":
                    stats["Dropped_Empty_ALT"] += 1
                    continue
                
                if filter_pass and filter_status not in ["PASS", "."]:
                    stats["Dropped_Not_PASS"] += 1
                    continue

                chromosome = fields[0].replace("chr", "")
                try: start_pos = int(fields[1])
                except ValueError: continue 
                
                ref_allele, tumor_seq_allele1 = fields[3], fields[3]
                end_pos = start_pos

                if len(ref_allele) == 1 and len(alt_allele) == 1: var_type = "SNP"
                elif len(ref_allele) > len(alt_allele): var_type = "DEL"
                elif len(ref_allele) < len(alt_allele): var_type = "INS"
                else: var_type = "Complex"

                variant_classification, hugo_symbol, transcript_name, protein_change = "NA", "NA", "NA", "NA"
                
                # ====== 1. 解析 INFO 字段 (提取 SnpEff ANN 以及 潜在的 DP/AF) ======
                info_field = fields[7]
                info_dict = {}
                for item in info_field.split(';'):
                    if '=' in item:
                        k, v = item.split('=', 1)
                        info_dict[k] = v
                    else:
                        info_dict[item] = True 
                
                # SnpEff ANN 解析
                if "ANN" in info_dict:
                    ann_fields = info_dict["ANN"].split(',')[0].split('|')
                    if len(ann_fields) > 10:
                        variant_classification = ann_fields[1] if ann_fields[1] else "NA"
                        hugo_symbol = ann_fields[3] if ann_fields[3] else "NA"
                        transcript_name = ann_fields[6] if ann_fields[6] else "NA"
                        protein_change = ann_fields[10] if ann_fields[10] else "NA"

                # ====== 2. 全域 DP/AF 搜捕引擎 ======
                dp, ad, tumor_vaf = "NA", "NA", "NA"
                
                # [触手 A]: 从 INFO 区强行挖掘 (LoFreq / iVar / FreeBayes 标准)
                if "DP" in info_dict and str(info_dict["DP"]).isdigit():
                    dp = int(info_dict["DP"])
                if "AF" in info_dict:
                    try:
                        # AF 可能是多等位基因逗号分隔的，取第一个
                        tumor_vaf = float(info_dict["AF"].split(',')[0])
                    except ValueError:
                        pass

                # [触手 B]: 试图从 FORMAT / SAMPLE 区提取 (GATK / Mutect2 标准，若存在则覆盖)
                if len(fields) >= 10:
                    format_fields = fields[8].split(':')
                    sample_values = fields[-1].split(':')
                    if 'DP' in format_fields:
                        dp_idx = format_fields.index('DP')
                        if dp_idx < len(sample_values) and sample_values[dp_idx].isdigit(): 
                            dp = int(sample_values[dp_idx])
                    if 'AD' in format_fields:
                        ad_idx = format_fields.index('AD')
                        if ad_idx < len(sample_values):
                            parts = sample_values[ad_idx].split(',')
                            if len(parts) >= 2 and parts[-1].isdigit(): 
                                ad = int(parts[-1])
                                # 如果能拿到特异的变异深度(AD)和总深度(DP)，亲自算一次保证极高准确性
                                if isinstance(dp, int) and dp > 0:
                                    tumor_vaf = float(ad) / float(dp)

                # ====== 3. 严格数据筛洗漏斗 ======
                if min_dp > 0:
                    # 如果此时 dp 还是 "NA" 拿不到，就证明真的没有深度信息，直接斩掉
                    if not isinstance(dp, int) or dp < min_dp:
                        stats["Dropped_Low_DP"] += 1
                        continue
                        
                if min_af > 0.0:
                    if not isinstance(tumor_vaf, float) or tumor_vaf < min_af:
                        stats["Dropped_Low_AF"] += 1
                        continue

                str_tumor_vaf = f"{tumor_vaf:.5f}" if isinstance(tumor_vaf, float) else "NA"
                
                maf_line =[
                    hugo_symbol, "NA", "NA", resolved_build, chromosome, str(start_pos), str(end_pos), "NA",
                    variant_classification, var_type, ref_allele, tumor_seq_allele1, alt_allele,
                    tumor_sample_barcode, protein_change, str_tumor_vaf, transcript_name
                ]
                outfile.write("\t".join(maf_line) + "\n")
                stats["Final_Retained_MAF"] += 1

            pbar.close() 

        if not quiet:
            time_elapsed = time.time() - start_time
            print("")
            print("┌" + "─" * 55 + "┐")
            print(f"│    [Statistical Breakdown] {tumor_sample_barcode[:25]:<25} │")
            print("├" + "─" * 55 + "┤")
            print(f"│  Total Input Mutational Events : {stats['Total_Input_Variants']:<19} │")
            print(f"│  > Removed (Filter != PASS)    : {stats['Dropped_Not_PASS']:<19} │")
            print(f"│  > Removed (Empty ALT Allele)  : {stats['Dropped_Empty_ALT']:<19} │")
            print(f"│  > Removed (DP Depth < {min_dp:<3})    : {stats['Dropped_Low_DP']:<19} │")
            print(f"│  > Removed (VAF/AF < {min_af:<4})    : {stats['Dropped_Low_AF']:<19} │")
            print("├" + "─" * 55 + "┤")
            print(f"│  FINAL MAF VARIANTS RETAINED   : {stats['Final_Retained_MAF']:<19} │")
            print("└" + "─" * 55 + "┘")
            print(f"   ✓ [SUCCESS] Metadata Genome Build = '{resolved_build}'")
            print(f"   ✓ [SUCCESS] Written to: {os.path.abspath(output_file)}")
            print(f"   ⏱ [TIME] Parsing completed in {time_elapsed:.2f} seconds.\n")

        return True, stats, resolved_build
        
    except Exception as e:
        if not quiet: print(f"[ERROR] An unexpected error occurred while parsing VCF: {e}", file=sys.stderr)
        return False, None, resolved_build
    finally:
        infile.close()

# ---------------------------------------------------------
# 主进程入口处理
# ---------------------------------------------------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Robustly convert viral SnpEff annotated VCF to MAF format."
    )
    parser.add_argument("input_path", nargs='?', help="Input SnpEff-annotated VCF file or directory")
    parser.add_argument("-i", "--input", type=str, help="Input VCF file or directory")
    parser.add_argument("-o", "--output", type=str, help="Output file/directory")
    parser.add_argument("-minDP", type=int, default=0, help="Minimum sequencing read depth (DP) filter. (default: 0)")
    parser.add_argument("-minAF", type=float, default=0.0, help="Minimum Allele Frequency (AF) filter. (default: 0.0)")
    parser.add_argument("--build", type=str, default="auto", help="NCBI Build metadata. (default: auto)")
    parser.add_argument("--filter-pass", action='store_true', dest="filterPASS", help="Filter strictly for PASS. (default: OFF)")

    args = parser.parse_args()
    input_path = os.path.abspath(args.input if args.input else args.input_path)
    if not input_path:
        parser.error("No input provided.")

    master_start_time = time.time()

    if os.path.isdir(input_path):
        vcf_files = glob.glob(os.path.join(input_path, "**", "*.vcf"), recursive=True) + \
                    glob.glob(os.path.join(input_path, "**", "*.vcf.gz"), recursive=True)
        
        if not vcf_files:
            print(f"[WARNING] No VCF files found deep inside directory:" )
            sys.exit(1)

        outdir = os.path.abspath(args.output) if args.output else None
        if outdir: os.makedirs(outdir, exist_ok=True)
        
        print(f"\n[SYSTEM] Aggregated Batch Mode Initiated. Processing Cohort (N={len(vcf_files)}) ...")

        cohort_stats = {
            "Total_Input_Variants": 0, "Dropped_Empty_ALT": 0, "Dropped_Not_PASS": 0,
            "Dropped_Low_DP": 0, "Dropped_Low_AF": 0, "Final_Retained_MAF": 0
        }
        
        for vcf_file in tqdm(vcf_files, desc="Converting VCFs", unit="file"):
            base = os.path.basename(vcf_file)
            if base.endswith('.gz'): base = base[:-3]
            if base.endswith('.vcf'): base = base[:-4]
            output_file_path = os.path.join(outdir, f"{base}.maf") if outdir else None
            
            success, stats, _ = parse_vcf(vcf_file, output_file_path, args.minDP, args.minAF, args.build, args.filterPASS, quiet=True)
            
            if success and stats:
                for k in cohort_stats: cohort_stats[k] += stats[k]
            
        time_elapsed = time.time() - master_start_time
        
        print("")
        print("┌" + "─" * 55 + "┐")
        print(f"│    [Cohort Aggregate Statistics] N = {len(vcf_files):<17} │")
        print("├" + "─" * 55 + "┤")
        print(f"│  Total Input Mutational Events : {cohort_stats['Total_Input_Variants']:<19} │")
        print(f"│  > Removed (Filter != PASS)    : {cohort_stats['Dropped_Not_PASS']:<19} │")
        print(f"│  > Removed (Empty ALT Allele)  : {cohort_stats['Dropped_Empty_ALT']:<19} │")
        print(f"│  > Removed (DP Depth < {args.minDP:<3})    : {cohort_stats['Dropped_Low_DP']:<19} │")
        print(f"│  > Removed (VAF/AF < {args.minAF:<4})    : {cohort_stats['Dropped_Low_AF']:<19} │")
        print("├" + "─" * 55 + "┤")
        print(f"│  FINAL MAF VARIANTS RETAINED   : {cohort_stats['Final_Retained_MAF']:<19} │")
        print("└" + "─" * 55 + "┘")
        print(f"   ✓[SUCCESS] Entire cohort smoothly parsed in {time_elapsed:.2f} s.\n")

    elif os.path.isfile(input_path):
        output_file_path = None
        if args.output:
            if os.path.isdir(args.output) or args.output.endswith(os.sep):
                outdir = os.path.abspath(args.output)
                os.makedirs(outdir, exist_ok=True)
                base = os.path.basename(input_path)
                if base.endswith('.gz'): base = base[:-3]
                if base.endswith('.vcf'): base = base[:-4]
                output_file_path = os.path.join(outdir, f"{base}.maf")
            else:
                output_file_path = os.path.abspath(args.output)
                
        print(f"\n[SYSTEM] Single-file Mode.")
        parse_vcf(input_path, output_file_path, args.minDP, args.minAF, args.build, args.filterPASS, quiet=False)
    else:
        print(f"[ERROR] Origin '{input_path}' is structurally invalid or does not exist.", file=sys.stderr)
        sys.exit(1)
