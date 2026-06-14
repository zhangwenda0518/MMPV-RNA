#!/usr/bin/env python3
import numpy as np
import pandas as pd
import argparse
import os
import textwrap
import re
import gc
import shutil
import subprocess
import multiprocessing
import tempfile
from concurrent.futures import ProcessPoolExecutor, as_completed

try:
    from tqdm import tqdm
except ImportError:
    print("⚠️ 未检测到 tqdm，建议运行 `conda install tqdm`。")
    def tqdm(iterable, **kwargs): return iterable

# ── 列名兼容映射 (新旧格式) ──
COL_ALIASES = {
    'Virus':       ['Rep_Accession', 'Virus', 'Accession'],
    'taxonomy':    ['Adjusted_Species', 'Species_NCBI', 'Species_ICTV', 'taxonomy'],
    'MeanDepth':   ['Rep_MeanDepth', 'MeanDepth'],
    'FPKM':        ['Asm_FPKM', 'FPKM'],
    'RPM':         ['Asm_RPM', 'RPM'],
    'TPM':         ['Asm_TPM', 'TPM'],
    'Coverage(%)': ['Rep_Coverage(%)', 'Coverage(%)'],
    'Length':      ['Rep_Length', 'Length'],
    'MappedReads': ['Asm_EM_Reads', 'MappedReads'],
}

def _normalize_df(df):
    """统一列名到标准格式"""
    for target, candidates in COL_ALIASES.items():
        for cand in candidates:
            if cand in df.columns:
                df = df.rename(columns={cand: target})
                break
    return df

def smooth(y, box_pts):
    if len(y) < box_pts: return y
    box = np.ones(box_pts) / box_pts
    return np.convolve(y, box, mode='same')

def create_info_text(v_stat):
    info =[
        f"Taxonomy: {textwrap.shorten(str(v_stat.get('taxonomy', 'Unknown')), width=45, placeholder='...')}",
        f"Accession: {v_stat['Virus']}",
        f"Length: {v_stat['Length']:,.0f} bp",
        f"Coverage: {v_stat['Coverage(%)']:.2f}%",
        f"Mean Depth: {v_stat['MeanDepth']:.2f}x",
        "-"*25,
        f"Mapped Reads: {v_stat['MappedReads']:,.0f}",
        f"RPM: {v_stat.get('RPM', 0):.2f}",
        f"FPKM: {v_stat.get('FPKM', 0):.2f}",
        f"TPM: {v_stat.get('TPM', 0):.2f}"
    ]
    return "\n".join(info)

def fetch_gb_via_efetch(accession, save_dir):
    os.makedirs(save_dir, exist_ok=True)
    out_file = os.path.join(save_dir, f"{accession}.gb")
    if os.path.exists(out_file) and os.path.getsize(out_file) > 0: return out_file
    cmd =["efetch", "-db", "nuccore", "-id", accession, "-format", "gb"]
    try:
        with open(out_file, "w") as f:
            subprocess.run(cmd, stdout=f, stderr=subprocess.PIPE, check=True)
        return out_file
    except subprocess.CalledProcessError:
        if os.path.exists(out_file): os.remove(out_file)
        return None

def parse_gb_with_biopython(gb_path):
    try:
        from Bio import SeqIO
    except ImportError: return[]
    
    features =[]
    try:
        for record in SeqIO.parse(gb_path, "genbank"):
            for feat in record.features:
                if feat.type in["CDS", "gene"]:
                    name = feat.qualifiers.get("gene", feat.qualifiers.get("product", ["Unknown"]))[0]
                    features.append({
                        'type': feat.type, 'start': int(feat.location.start), 
                        'end': int(feat.location.end), 'strand': feat.location.strand, 'name': name
                    })
    except Exception: return[]
    
    seen = {}
    for f in features:
        if f['name'] == 'Unknown': continue
        coord = (f['start'], f['end'])
        if coord not in seen: seen[coord] = f
        else:
            if f['type'] == 'CDS' and seen[coord]['type'] == 'gene': seen[coord] = f
            elif seen[coord]['name'] == 'Unknown' and f['name'] != 'Unknown': seen[coord]['name'] = f['name']
    final = list(seen.values())
    final.sort(key=lambda x: x['start'])
    return final

# ================= 核心重构：解除 rg 一切限制的终极提取 =================
def process_single_sample(task_args):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import matplotlib.patches as patches
    
    sample, v_stats_list, depth_file, sample_out_dir, window, fontsize, genes_dict, threads_per_job = task_args
    results_log =[]
    
    try:
        valid_viruses = [str(v['Virus']).strip() for v in v_stats_list]
        valid_set = set(valid_viruses)
        
        # 1. 严格写入匹配规则，确保带 Tab 符匹配第一列
        fd, tmp_pattern_file = tempfile.mkstemp(suffix='.txt', prefix='rg_patterns_')
        with os.fdopen(fd, 'w') as f:
            for vir in valid_set: f.write(f"{vir}\t\n")
        
        has_rg = shutil.which('rg') is not None
        has_rapidgzip = shutil.which('rapidgzip') is not None
        has_pigz = shutil.which('pigz') is not None

        # --- 【修复核心】：追加 -a, --no-ignore, -N, -I 彻底解除 rg 的束缚 ---
        rg_flags = f"-a --no-ignore -N -I -F -j {threads_per_job}"
        grep_flags = "-a -F"
        
        if depth_file.endswith('.gz'):
            if has_rapidgzip:
                cat_cmd = f"rapidgzip -d -c -k -P {threads_per_job} '{depth_file}'"
            elif has_pigz:
                cat_cmd = f"pigz -dc -p {threads_per_job} '{depth_file}'"
            else:
                cat_cmd = f"gzip -dc '{depth_file}'"
                
            if has_rg:
                grep_cmd = f"rg {rg_flags} -f '{tmp_pattern_file}'"
            else:
                grep_cmd = f"LC_ALL=C grep {grep_flags} -f '{tmp_pattern_file}'"
                
            full_cmd = f"{cat_cmd} | {grep_cmd}"
        else:
            if has_rg:
                full_cmd = f"rg {rg_flags} -f '{tmp_pattern_file}' '{depth_file}'"
            else:
                full_cmd = f"LC_ALL=C grep {grep_flags} -f '{tmp_pattern_file}' '{depth_file}'"
        
        # 2. 捕捉 stderr 以便后续诊断
        proc = subprocess.Popen(full_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        
        err_msg = ""
        try:
            depth_df = pd.read_csv(proc.stdout, sep='\t', header=None, names=['chr', 'position', 'depth'])
        except pd.errors.EmptyDataError:
            depth_df = pd.DataFrame(columns=['chr', 'position', 'depth'])
            # 只有提取失败时，才读取底层的报错信息
            err_msg = proc.stderr.read().strip()
        finally:
            proc.stdout.close()
            proc.stderr.close()
            proc.wait()
            if os.path.exists(tmp_pattern_file):
                os.remove(tmp_pattern_file) 

        # 如果真的提取为空，打印出底层的命令报错（如果有）
        if depth_df.empty:
            fail_reason = "提取为空 (文件内没有匹配的坐标行)"
            if err_msg:
                fail_reason += f" | 底层命令报错: {err_msg[:200]}"
            return sample, False, fail_reason,[]

        # 3. 开始正常绘图...
        for v_stat in v_stats_list:
            virus = str(v_stat['Virus']).strip()
            chrom_df = depth_df[depth_df['chr'] == virus]
            if chrom_df.empty:
                results_log.append((virus, False, "提取结果中未包含该株的坐标"))
                continue

            x, y = chrom_df['position'].values, chrom_df['depth'].values
            virus_genes = genes_dict.get(virus, genes_dict.get(virus.split('.')[0],[]))
            mean_depth = v_stat['MeanDepth']
            info_text = create_info_text(v_stat)
            
            tax = str(v_stat.get('taxonomy', 'Unknown'))
            safe_tax = re.sub(r'[^\w\-.]', '_', tax)
            safe_vname = re.sub(r'[^\w\-.]', '_', virus)
            title = f"Sample: {sample}\n{textwrap.shorten(tax, width=65, placeholder='...')} ({virus})"
            output_path = os.path.join(sample_out_dir, f"{sample}_{safe_tax}_{safe_vname}_depth.pdf")

            y_smooth = smooth(y, window)
            max_x = max(x) if len(x) > 0 else (virus_genes[-1]['end'] if virus_genes else 1000)

            if virus_genes:
                fig, (ax, ax_g) = plt.subplots(2, 1, figsize=(12, 7.5), sharex=True, gridspec_kw={'height_ratios':[5, 1], 'hspace': 0.08})
            else:
                fig, ax = plt.subplots(figsize=(12, 6))
                ax_g = None

            ax.fill_between(x, 0, y_smooth, alpha=0.3, color='#1f77b4')
            ax.plot(x, y_smooth, color='#1f77b4', linewidth=1.5, label=f'Smoothed Depth (Window={window})')
            ax.axhline(y=mean_depth, color='#d62728', linestyle='--', linewidth=2, label=f'Mean Depth: {mean_depth:.2f}x')
            ax.annotate(info_text, xy=(0.02, 0.96), xycoords='axes fraction', bbox=dict(boxstyle="round,pad=0.6", fc="#f8f9fa", ec="#ced4da", alpha=0.9), fontsize=fontsize, family='monospace', ha='left', va='top')

            ax.set_title(title, fontsize=14, fontweight='bold', pad=15)
            ax.set_ylabel('Sequencing Depth (x)', fontsize=12)
            ax.legend(loc='upper right', bbox_to_anchor=(0.98, 0.98), framealpha=0.9)
            ax.grid(True, linestyle=':', alpha=0.6)
            ax.set_ylim(bottom=0)
            ax.set_xlim(0, max_x)

            if ax_g is not None:
                colors =['#8dd3c7', '#ffffb3', '#bebada', '#fb8072', '#80b1d3', '#fdb462', '#b3de69', '#fccde5', '#d9d9d9']
                y_center, height = 0.5, 0.4
                for i, gene in enumerate(virus_genes):
                    start, end = gene['start'], gene['end']
                    head_length = min(max_x * 0.015, (end - start) * 0.4) 
                    y_top, y_bottom = y_center + height/2, y_center - height/2
                    if gene['strand'] >= 0: x_coords =[start, end - head_length, end, end - head_length, start]
                    else: x_coords =[end, start + head_length, start, start + head_length, end]
                    
                    poly = patches.Polygon(np.column_stack((x_coords,[y_bottom, y_bottom, y_center, y_top, y_top])), 
                                           closed=True, facecolor=colors[i % len(colors)], edgecolor='#444444', alpha=0.9)
                    ax_g.add_patch(poly)
                    short_name = gene['name'][:10] + '..' if len(gene['name']) > 12 else gene['name']
                    ax_g.text(start + (end - start)/2, y_center, short_name, ha='center', va='center', fontsize=8, fontweight='bold')
                
                ax_g.set_ylim(0, 1)
                ax_g.set_yticks([])
                for spine in['top', 'right', 'left']: ax_g.spines[spine].set_visible(False)
                ax_g.set_xlabel('Genome Position (bp)', fontsize=12)
                ax.tick_params(labelbottom=False)
            else:
                ax.set_xlabel('Genome Position (bp)', fontsize=12)

            plt.tight_layout()
            plt.savefig(output_path, bbox_inches='tight', dpi=300)
            plt.close(fig) 
            results_log.append((virus, True, output_path))
            
        del depth_df
        gc.collect()
        return sample, True, "极速解析与绘图成功", results_log

    except Exception as e:
        return sample, False, f"发生致命崩溃: {str(e)}",[]


def main():
    parser = argparse.ArgumentParser(description='生信终极测序深度绘图工具（强制绕过 rg 保护限制版）')
    parser.add_argument('-d', '--depth_dir', required=True)
    parser.add_argument('-m', '--summary', required=True)
    parser.add_argument('-o', '--output', required=True)
    parser.add_argument('-n', '--sample', required=False)
    parser.add_argument('-g', '--add_genes', action='store_true')
    parser.add_argument('--gbk_dir', required=False, default='gbk_files')
    parser.add_argument('-w', '--window', type=int, default=100)
    parser.add_argument('-f', '--fontsize', type=float, default=9)
    parser.add_argument('-t', '--threads', type=int, default=0, help="处理的样本并发数。")
    args = parser.parse_args()

    total_cores = multiprocessing.cpu_count()
    max_workers = args.threads if args.threads > 0 else min(16, max(1, total_cores - 1))
    threads_per_job = max(1, min(4, total_cores // max_workers))
    
    print("="*60)
    print(f"🚀 生信图表极速渲染引擎已启动 (样本并发数: {max_workers})")
    print(f"⚡ 检测底层加速器: rapidgzip [{'✅' if shutil.which('rapidgzip') else '❌'}] | ripgrep (rg) [{'✅' if shutil.which('rg') else '❌'}]")
    print("="*60)

    summary_df = pd.read_csv(args.summary, sep=',' if args.summary.endswith('.csv') else '\t')
    summary_df = _normalize_df(summary_df)
    samples_to_process =[args.sample] if args.sample else summary_df['Sample'].unique().tolist()
    
    genes_dict = {}
    if args.add_genes:
        print(f"\n🧬[双轨模式已开启] 准备 GenBank 注释...")
        unique_viruses = summary_df[summary_df['Sample'].isin(samples_to_process)]['Virus'].unique()
        os.makedirs(args.gbk_dir, exist_ok=True)
        for virus_id in tqdm(unique_viruses, desc="提取注释", leave=False):
            gb_path = fetch_gb_via_efetch(virus_id, args.gbk_dir) or os.path.join(args.gbk_dir, f"{virus_id}.gb")
            if os.path.exists(gb_path):
                features = parse_gb_with_biopython(gb_path)
                if features:
                    genes_dict[virus_id] = features
                    genes_dict[virus_id.split('.')[0]] = features 

    os.makedirs(args.output, exist_ok=True)
    executor = ProcessPoolExecutor(max_workers=max_workers)
    future_to_sample = {}

    print(f"\n▶ 正在下发智能解析任务...")
    
    for sample in samples_to_process:
        sample_summary = summary_df[summary_df['Sample'] == sample]
        if sample_summary.empty: continue

        possible_files =[
            f"{sample}.SiteDepth", f"{sample}.chr.SiteDepth",
            f"{sample}.SiteDepth.gz", f"{sample}.chr.SiteDepth.gz"
        ]
        depth_file = None
        for pf in possible_files:
            pf_path = os.path.join(args.depth_dir, pf)
            if os.path.exists(pf_path):
                depth_file = pf_path
                break
                
        if not depth_file: 
            print(f"⚠️ 找不到样本 {sample} 的深度文件，跳过...")
            continue

        sample_out_dir = os.path.join(args.output, str(sample))
        os.makedirs(sample_out_dir, exist_ok=True)

        task_args = (sample, sample_summary.to_dict('records'), depth_file, sample_out_dir, args.window, args.fontsize, genes_dict, threads_per_job)
        future = executor.submit(process_single_sample, task_args)
        future_to_sample[future] = sample

    success_plots, total_plots = 0, 0

    with tqdm(total=len(future_to_sample), desc="📊 总体进度", unit="sample") as pbar:
        for future in as_completed(future_to_sample):
            sample = future_to_sample[future]
            try:
                s_name, s_success, s_msg, results_log = future.result()
                if not s_success:
                    # 如果再次失败，这里会打印出 rg 底层的真实报错（比如内存不足、权限问题等）
                    tqdm.write(f"⚠️[{sample} 失败] {s_msg}")
                else:
                    for virus, v_success, info in results_log:
                        total_plots += 1
                        if v_success:
                            success_plots += 1
                            tqdm.write(f"  ✅ [{sample}] 绘制成功: {os.path.basename(info)}")
                        else:
                            tqdm.write(f"  ❌ [{sample}] {virus} 画图报错: {info}")
            except Exception as exc:
                tqdm.write(f"💥 [严重崩溃] {sample}: {exc}")
            pbar.update(1)

    executor.shutdown()
    print(f"\n🎉 运行完成！(生成图表: {success_plots} / {total_plots})")

if __name__ == "__main__":
    main()
