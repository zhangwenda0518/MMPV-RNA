#!/usr/bin/env python3
"""
batch_plot_virus_depth.py — 病毒组学三合一可视化引擎
=================================================================
模式:
  --mode all      一键全部生成 (默认)
  --mode depth    测序深度图 (per-virus coverage + 基因轨道)
  --mode freq     丰度分布箱线图 (MeanDepth/FPKM/RPM/TPM)
  --mode meta     病毒 vs 样本元数据关联图

深度模式 (--mode depth):
  python batch_plot_virus_depth.py --mode depth \\
    -d stat/ -m all_viruses.best.summary.tsv -o plots/ -t 8 -g

丰度模式 (--mode freq):
  python batch_plot_virus_depth.py --mode freq \\
    -m all_viruses.best.summary.tsv -o plots/ --log10

元数据模式 (--mode meta):
  python batch_plot_virus_depth.py --mode meta \\
    -m all_viruses.best.summary.tsv --meta metadata.tsv -o plots/
"""

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.ticker as ticker

import argparse
import io
import os
import sys
import textwrap
import re
import gc
import shlex
import shutil
import subprocess
import multiprocessing
import tempfile
import warnings
from concurrent.futures import ProcessPoolExecutor, as_completed

warnings.filterwarnings("ignore", category=UserWarning)

try:
    from tqdm import tqdm
except ImportError:
    def tqdm(iterable, **kwargs): return iterable

# =====================================================================
# 通用工具
# =====================================================================

def smooth(y, box_pts):
    if len(y) < box_pts: return y
    box = np.ones(box_pts) / box_pts
    return np.convolve(y, box, mode='same')

def safe_name(s):
    return re.sub(r'[^\w\-.]', '_', str(s)).strip('_')

def create_info_text(v_stat):
    info = [
        f"Taxonomy: {textwrap.shorten(str(v_stat.get('taxonomy', 'Unknown')), width=45, placeholder='...')}",
        f"Accession: {v_stat['Virus']}",
        f"Length: {v_stat['Length']:,.0f} bp",
        f"Coverage: {v_stat['Coverage(%)']:.2f}%",
        f"Mean Depth: {v_stat['MeanDepth']:.2f}x",
        "-" * 25,
        f"Mapped Reads: {v_stat['MappedReads']:,.0f}",
        f"RPM: {v_stat.get('RPM', 0):.2f}",
        f"FPKM: {v_stat.get('FPKM', 0):.2f}",
        f"TPM: {v_stat.get('TPM', 0):.2f}"
    ]
    return "\n".join(info)

def fetch_gb_via_efetch(accession, save_dir):
    os.makedirs(save_dir, exist_ok=True)
    out_file = os.path.join(save_dir, f"{accession}.gb")
    if os.path.exists(out_file) and os.path.getsize(out_file) > 0:
        return out_file
    cmd = ["efetch", "-db", "nuccore", "-id", accession, "-format", "gb"]
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
    except ImportError:
        return []
    features = []
    try:
        for record in SeqIO.parse(gb_path, "genbank"):
            for feat in record.features:
                if feat.type in ["CDS", "gene"]:
                    name = feat.qualifiers.get("gene", feat.qualifiers.get("product", ["Unknown"]))[0]
                    features.append({
                        'type': feat.type, 'start': int(feat.location.start),
                        'end': int(feat.location.end), 'strand': feat.location.strand, 'name': name
                    })
    except Exception:
        return []
    seen = {}
    for f in features:
        if f['name'] == 'Unknown': continue
        coord = (f['start'], f['end'])
        if coord not in seen:
            seen[coord] = f
        else:
            if f['type'] == 'CDS' and seen[coord]['type'] == 'gene':
                seen[coord] = f
            elif seen[coord]['name'] == 'Unknown' and f['name'] != 'Unknown':
                seen[coord]['name'] = f['name']
    final = list(seen.values())
    final.sort(key=lambda x: x['start'])
    return final


# =====================================================================
# 模块 1: 测序深度图 (原 batch_plot_virus_depth)
# =====================================================================

def process_single_sample(task_args):
    sample, v_stats_list, depth_file, sample_out_dir, window, fontsize, genes_dict, threads_per_job = task_args
    results_log = []

    try:
        valid_viruses = [str(v['Virus']).strip() for v in v_stats_list]
        valid_set = set(valid_viruses)
        if not valid_set:
            return sample, False, "无有效病毒可处理", []

        tmp_dir = Path(os.environ.get("TMPDIR", "/tmp"))
        fd, tmp_pattern_file = tempfile.mkstemp(suffix='.txt', prefix='rg_patterns_', dir=str(tmp_dir))
        with os.fdopen(fd, 'w') as f:
            for vir in valid_set: f.write(f"{vir}\n")

        has_rg = shutil.which('rg') is not None
        has_rapidgzip = shutil.which('rapidgzip') is not None
        has_pigz = shutil.which('pigz') is not None

        rg_flags = f"-a --no-ignore -N -I -F -j {threads_per_job}"
        grep_flags = "-a -F"

        if depth_file.endswith('.gz'):
            if has_rapidgzip:
                cat_cmd = f"rapidgzip -d -c -k -P {threads_per_job} {shlex.quote(depth_file)}"
            elif has_pigz:
                cat_cmd = f"pigz -dc -p {threads_per_job} {shlex.quote(depth_file)}"
            else:
                cat_cmd = f"gzip -dc {shlex.quote(depth_file)}"
            if has_rg:
                grep_cmd = f"rg {rg_flags} -f '{tmp_pattern_file}'"
            else:
                grep_cmd = f"LC_ALL=C grep {grep_flags} -f '{tmp_pattern_file}'"
            full_cmd = f"{cat_cmd} | {grep_cmd}"
        else:
            if has_rg:
                full_cmd = f"rg {rg_flags} -f {shlex.quote(tmp_pattern_file)} {shlex.quote(depth_file)}"
            else:
                full_cmd = f"LC_ALL=C grep {grep_flags} -f {shlex.quote(tmp_pattern_file)} {shlex.quote(depth_file)}"

        result = subprocess.run(full_cmd, shell=True, capture_output=True, text=True)
        err_msg = result.stderr.strip() if result.returncode != 0 else ""
        try:
            depth_df = pd.read_csv(io.StringIO(result.stdout), sep='\t', header=None, names=['chr', 'position', 'depth'])
        except (pd.errors.EmptyDataError, pd.errors.ParserError):
            depth_df = pd.DataFrame(columns=['chr', 'position', 'depth'])
            err_msg = err_msg or result.stderr.strip()
        finally:
            if os.path.exists(tmp_pattern_file):
                os.remove(tmp_pattern_file)

        if depth_df.empty:
            fail_reason = "提取为空 (文件内没有匹配的坐标行)"
            if err_msg:
                fail_reason += f" | 底层命令报错: {err_msg[:200]}"
            return sample, False, fail_reason, []

        for v_stat in v_stats_list:
            virus = str(v_stat['Virus']).strip()
            chrom_df = depth_df[depth_df['chr'] == virus]
            if chrom_df.empty:
                results_log.append((virus, False, "提取结果中未包含该株的坐标"))
                continue

            x, y = chrom_df['position'].values, chrom_df['depth'].values
            virus_genes = genes_dict.get(virus, genes_dict.get(virus.split('.')[0], []))
            mean_depth = v_stat['MeanDepth']
            info_text = create_info_text(v_stat)

            tax = str(v_stat.get('taxonomy', 'Unknown'))
            safe_tax = safe_name(tax)
            safe_vname = safe_name(virus)
            title = f"Sample: {sample}\n{textwrap.shorten(tax, width=65, placeholder='...')} ({virus})"
            output_path = os.path.join(sample_out_dir, f"{sample}_{safe_tax}_{safe_vname}_depth.pdf")

            y_smooth = smooth(y, window)
            max_x = max(x) if len(x) > 0 else (virus_genes[-1]['end'] if virus_genes else 1000)

            if virus_genes:
                fig, (ax, ax_g) = plt.subplots(2, 1, figsize=(12, 7.5), sharex=True,
                                               gridspec_kw={'height_ratios': [5, 1], 'hspace': 0.08})
            else:
                fig, ax = plt.subplots(figsize=(12, 6))
                ax_g = None

            ax.fill_between(x, 0, y_smooth, alpha=0.3, color='#1f77b4')
            ax.plot(x, y_smooth, color='#1f77b4', linewidth=1.5, label=f'Smoothed Depth (Window={window})')
            ax.axhline(y=mean_depth, color='#d62728', linestyle='--', linewidth=2, label=f'Mean Depth: {mean_depth:.2f}x')
            ax.annotate(info_text, xy=(0.02, 0.96), xycoords='axes fraction',
                        bbox=dict(boxstyle="round,pad=0.6", fc="#f8f9fa", ec="#ced4da", alpha=0.9),
                        fontsize=fontsize, family='monospace', ha='left', va='top')

            ax.set_title(title, fontsize=14, fontweight='bold', pad=15)
            ax.set_ylabel('Sequencing Depth (x)', fontsize=12)
            ax.legend(loc='upper right', bbox_to_anchor=(0.98, 0.98), framealpha=0.9)
            ax.grid(True, linestyle=':', alpha=0.6)
            ax.set_ylim(bottom=0)
            ax.set_xlim(0, max_x)

            if ax_g is not None:
                colors = ['#8dd3c7', '#ffffb3', '#bebada', '#fb8072', '#80b1d3', '#fdb462', '#b3de69', '#fccde5',
                          '#d9d9d9']
                y_center, height = 0.5, 0.4
                for i, gene in enumerate(virus_genes):
                    start, end = gene['start'], gene['end']
                    head_length = min(max_x * 0.015, (end - start) * 0.4)
                    y_top, y_bottom = y_center + height / 2, y_center - height / 2
                    if gene['strand'] >= 0:
                        x_coords = [start, end - head_length, end, end - head_length, start]
                    else:
                        x_coords = [end, start + head_length, start, start + head_length, end]
                    poly = patches.Polygon(
                        np.column_stack((x_coords, [y_bottom, y_bottom, y_center, y_top, y_top])),
                        closed=True, facecolor=colors[i % len(colors)], edgecolor='#444444', alpha=0.9)
                    ax_g.add_patch(poly)
                    short_name = gene['name'][:10] + '..' if len(gene['name']) > 12 else gene['name']
                    ax_g.text(start + (end - start) / 2, y_center, short_name, ha='center', va='center',
                              fontsize=8, fontweight='bold')
                ax_g.set_ylim(0, 1)
                ax_g.set_yticks([])
                for spine in ['top', 'right', 'left']: ax_g.spines[spine].set_visible(False)
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
        return sample, False, f"发生致命崩溃: {str(e)}", []


def run_depth_mode(args):
    total_cores = multiprocessing.cpu_count()
    max_workers = args.threads if args.threads > 0 else min(16, max(1, total_cores - 1))
    threads_per_job = max(1, min(4, total_cores // max_workers))

    print("=" * 60)
    print(f"🚀 模式: 测序深度图 (样本并发数: {max_workers})")
    print(f"⚡ 底层加速: rapidgzip [{'✅' if shutil.which('rapidgzip') else '❌'}]  "
          f"ripgrep [{'✅' if shutil.which('rg') else '❌'}]")
    print("=" * 60)

    summary_df = pd.read_csv(args.summary, sep=',' if args.summary.endswith('.csv') else '\t')
    samples_to_process = [args.sample] if args.sample else summary_df['Sample'].unique().tolist()

    genes_dict = {}
    if args.add_genes:
        print(f"\n🧬 准备 GenBank 注释...")
        unique_viruses = summary_df[summary_df['Sample'].isin(samples_to_process)]['Virus'].unique()
        os.makedirs(args.gbk_dir, exist_ok=True)
        for virus_id in tqdm(unique_viruses, desc="提取注释"):
            gb_path = fetch_gb_via_efetch(virus_id, args.gbk_dir) or os.path.join(args.gbk_dir, f"{virus_id}.gb")
            if os.path.exists(gb_path):
                feats = parse_gb_with_biopython(gb_path)
                if feats:
                    genes_dict[virus_id] = feats
                    genes_dict[virus_id.split('.')[0]] = feats

    os.makedirs(args.output, exist_ok=True)
    executor = ProcessPoolExecutor(max_workers=max_workers)
    future_to_sample = {}

    print(f"\n▶ 下发智能解析任务...")
    for sample in samples_to_process:
        sample_summary = summary_df[summary_df['Sample'] == sample]
        if sample_summary.empty: continue
        possible_files = [
            f"{sample}.SiteDepth", f"{sample}.chr.SiteDepth",
            f"{sample}.SiteDepth.gz", f"{sample}.chr.SiteDepth.gz"
        ]
        depth_file = None
        for pf in possible_files:
            pf_path = os.path.join(args.depth_dir, pf)
            if os.path.exists(pf_path):
                depth_file = pf_path; break
        if not depth_file:
            print(f"⚠️ 找不到样本 {sample} 的深度文件，跳过...")
            continue
        sample_out_dir = os.path.join(args.output, str(sample))
        os.makedirs(sample_out_dir, exist_ok=True)
        task = (sample, sample_summary.to_dict('records'), depth_file, sample_out_dir,
                args.window, args.fontsize, genes_dict, threads_per_job)
        future_to_sample[executor.submit(process_single_sample, task)] = sample

    success_plots, total_plots = 0, 0
    with tqdm(total=len(future_to_sample), desc="📊 深度图进度", unit="sample") as pbar:
        for future in as_completed(future_to_sample):
            sample = future_to_sample[future]
            try:
                s_name, s_success, s_msg, results_log = future.result()
                if not s_success:
                    tqdm.write(f"⚠️ [{sample} 失败] {s_msg}")
                else:
                    for virus, v_success, info in results_log:
                        total_plots += 1
                        if v_success:
                            success_plots += 1
                            tqdm.write(f"  ✅ [{sample}] {os.path.basename(info)}")
                        else:
                            tqdm.write(f"  ❌ [{sample}] {virus}: {info}")
            except Exception as exc:
                tqdm.write(f"💥 [{sample}] {exc}")
            pbar.update(1)
    executor.shutdown()
    print(f"\n🎉 深度图完成！(成功: {success_plots} / {total_plots})")


# =====================================================================
# 模块 2: 丰度分布箱线图 (R virus_frequency_plot.R → Python)
# =====================================================================

def run_frequency_mode(args):
    print("=" * 60)
    print(f"📊 模式: 病毒丰度分布箱线图")
    print("=" * 60)

    df = pd.read_csv(args.summary, sep=',' if args.summary.endswith('.csv') else '\t')

    # 智能 Display_Name
    if 'taxonomy' in df.columns:
        df['Display_Name'] = df.apply(
            lambda r: f"{r['taxonomy']}\n({r['Virus']})"
            if str(r.get('taxonomy', 'Unannotated')) not in ['Unannotated', '-', 'nan']
            else str(r['Virus']),
            axis=1
        )
    elif 'Adjusted_Species' in df.columns:
        df['Display_Name'] = df.apply(
            lambda r: f"{r['Adjusted_Species']}\n({r.get('Rep_Accession', r.get('Virus', ''))})", axis=1
        )
    else:
        df['Display_Name'] = df['Virus'] if 'Virus' in df.columns else df.iloc[:, 0]

    # 文本换行
    df['Display_Name'] = df['Display_Name'].apply(
        lambda x: '\n'.join(textwrap.wrap(str(x), width=40))
    )

    # 可用的丰度指标
    metric_col_map = {
        'MeanDepth': 'MeanDepth',
        'FPKM': 'FPKM',
        'RPM': 'RPM',
        'TPM': 'TPM',
        'CPM': 'CPM',
        'Coverage(%)': 'Coverage(%)',
    }

    if args.freq_metrics:
        metrics = [m.strip() for m in args.freq_metrics.split(',')]
    else:
        metrics = ['MeanDepth', 'FPKM', 'RPM', 'TPM']

    available = [m for m in metrics if m in df.columns]
    if not available:
        print(f"❌ 未找到指定的丰度列。可用列: {list(df.columns)}")
        sys.exit(1)

    os.makedirs(args.output, exist_ok=True)
    unique_count = df['Display_Name'].nunique()
    should_flip = unique_count > 5

    # 颜色映射
    viridis_colors = plt.cm.turbo(np.linspace(0, 1, max(unique_count, 1)))

    # Shared color map for consistent coloring across all subplots
    all_names = df['Display_Name'].unique()
    color_map = {name: viridis_colors[i % len(viridis_colors)] for i, name in enumerate(all_names)}

    def make_boxplot(data, value_col, out_prefix):
        plot_df = data.dropna(subset=[value_col]).copy()
        if plot_df.empty: return

        # 按中位数排序
        medians = plot_df.groupby('Display_Name')[value_col].median().sort_values()
        plot_df['Display_Name'] = pd.Categorical(plot_df['Display_Name'], categories=medians.index, ordered=True)

        fig, ax = plt.subplots(figsize=(args.width, max(6, unique_count * 0.5)))

        positions = list(range(len(medians)))

        bp = ax.boxplot(
            [plot_df[plot_df['Display_Name'] == n][value_col].values for n in medians.index],
            positions=positions, patch_artist=True, widths=0.6,
            showfliers=False, vert=not should_flip
        )

        for i, (patch, name) in enumerate(zip(bp['boxes'], medians.index)):
            patch.set_facecolor(color_map[name])
            patch.set_alpha(0.6)

        for i, name in enumerate(medians.index):
            vals = plot_df[plot_df['Display_Name'] == name][value_col].values
            jitter = np.random.uniform(-0.2, 0.2, len(vals))
            if should_flip:
                ax.scatter(vals, np.full(len(vals), i) + jitter, alpha=0.7, s=args.point_size,
                          color=color_map[name], edgecolors='white', linewidth=0.3)
            else:
                ax.scatter(np.full(len(vals), i) + jitter, vals, alpha=0.7, s=args.point_size,
                          color=color_map[name], edgecolors='white', linewidth=0.3)

        if should_flip:
            ax.set_yticks(positions)
            ax.set_yticklabels(medians.index, fontsize=9, style='italic')
            ax.set_xlabel(value_col, fontweight='bold', fontsize=12)
        else:
            ax.set_xticks(positions)
            ax.set_xticklabels(medians.index, fontsize=9, style='italic', rotation=45, ha='right')
            ax.set_ylabel(value_col, fontweight='bold', fontsize=12)

        if args.log10:
            ax.set_xscale('log') if should_flip else ax.set_yscale('log')

        ax.set_title(f'Virus Abundance — {value_col} Distribution', fontweight='bold', fontsize=14, pad=15)
        ax.grid(axis='y' if should_flip else 'x', linestyle=':', alpha=0.4)

        plt.tight_layout()
        fmt = args.format
        fn = os.path.join(args.output, f"{out_prefix}_{value_col.lower()}{'_log10' if args.log10 else ''}.{fmt}")
        plt.savefig(fn, dpi=args.dpi, bbox_inches='tight')
        plt.close()
        print(f"  📈 {os.path.basename(fn)}")

    # 逐一单指标图
    for m in available:
        make_boxplot(df, m, "freq")

    # 综合面板图
    if args.multi_plot and len(available) > 1:
        print(f"  📐 生成综合多指标面板图...")
        plot_data = df.melt(id_vars=['Display_Name'], value_vars=available,
                            var_name='Metric', value_name='Value').dropna()

        base_metric = available[0]
        medians = plot_data[plot_data['Metric'] == base_metric].groupby('Display_Name')['Value'].median().sort_values()
        plot_data['Display_Name'] = pd.Categorical(plot_data['Display_Name'], categories=medians.index, ordered=True)

        n_cols = len(available)
        fig, axes = plt.subplots(1, n_cols, figsize=(args.width * 1.5, max(args.height, unique_count * 0.5)))
        if n_cols == 1: axes = [axes]

        for ax, m in zip(axes, available):
            mdf = plot_data[plot_data['Metric'] == m]
            positions = list(range(len(medians)))
            bp = ax.boxplot(
                [mdf[mdf['Display_Name'] == n]['Value'].values for n in medians.index],
                positions=positions, patch_artist=True, widths=0.6,
                showfliers=False, vert=False
            )
            for patch, name in zip(bp['boxes'], medians.index):
                c = color_map.get(name, '#1f77b4')
                patch.set_facecolor(c); patch.set_alpha(0.6)
            for i, name in enumerate(medians.index):
                vals = mdf[mdf['Display_Name'] == name]['Value'].values
                jitter = np.random.uniform(-0.15, 0.15, len(vals))
                ax.scatter(vals, np.full(len(vals), i) + jitter, alpha=0.5, s=args.point_size * 0.7,
                          color=color_map.get(name, '#1f77b4'), edgecolors='white', linewidth=0.2)
            ax.set_yticks(positions)
            ax.set_yticklabels(medians.index, fontsize=8, style='italic')
            ax.set_title(m, fontweight='bold', fontsize=11)
            if args.log10: ax.set_xscale('log')

        fig.suptitle('Comprehensive Multi-metric Virus Abundance', fontweight='bold', fontsize=16, y=1.01)
        plt.tight_layout()
        fn = os.path.join(args.output, f"freq_multi_metrics{'_log10' if args.log10 else ''}.{args.format}")
        plt.savefig(fn, dpi=args.dpi, bbox_inches='tight')
        plt.close()
        print(f"  📈 {os.path.basename(fn)}")

    print("\n🎉 丰度分布图完成！")


# =====================================================================
# 模块 3: 病毒 vs 元数据关联图 (原 virus_metadata_plot.py)
# =====================================================================

def clean_ai_label(val):
    s = str(val).strip()
    s = re.sub(r'_AI$', '', s, flags=re.IGNORECASE)
    if s.lower() in ['not_provided', 'unknown', 'nan', 'none', '', '<na>']:
        return 'Unknown'
    return s

def run_metadata_mode(args):
    if not args.meta:
        # 尝试从默认路径找
        script_dir = os.path.dirname(os.path.abspath(__file__))
        gsa_script = os.path.join(os.path.dirname(script_dir), "metadata", "gsa_sra.info.py")
        default_outdir = args.meta_outdir or os.path.join(os.path.dirname(args.output), "sra_results")
        default_meta = os.path.join(default_outdir, "Global_Unified_Metadata_Core13.tsv")

        if os.path.exists(default_meta):
            meta_file = default_meta
        elif args.sra_list and os.path.exists(args.sra_list) and os.path.exists(gsa_script):
            print(f"⏳ 未找到元数据，自动运行 metadata/gsa_sra.info.py ...")
            cmd = [sys.executable, gsa_script, "-i", args.sra_list, "-o", default_outdir,
                   "-m", "both", "-t", str(args.meta_threads)]
            if args.deepseek_api: cmd += ["--deepseek-api", args.deepseek_api]
            if args.ncbi_api: cmd += ["--ncbi-api", args.ncbi_api]
            if args.email: cmd += ["--email", args.email]
            if subprocess.run(cmd).returncode != 0:
                print(f"❌ 元数据生成失败"); sys.exit(1)
            meta_file = default_meta
        else:
            print("❌ 请提供 --meta 或 --sra_list"); sys.exit(1)
    else:
        meta_file = args.meta

    if not os.path.exists(meta_file):
        print(f"❌ 元数据文件不存在: {meta_file}"); sys.exit(1)

    print("=" * 60)
    print(f"📊 模式: 病毒 vs 元数据关联分析")
    print(f"  Virus:  {args.summary}")
    print(f"  Meta:   {meta_file}")
    print("=" * 60)

    os.makedirs(args.output, exist_ok=True)

    # 读取
    df_v = pd.read_csv(args.summary, sep=',' if args.summary.endswith('.csv') else '\t')
    df_m = pd.read_csv(meta_file, sep=',' if meta_file.endswith('.csv') else '\t')

    # 列清洗
    df_v.columns = [c.strip().replace('\ufeff', '') for c in df_v.columns]
    df_m.columns = [c.strip().replace('\ufeff', '') for c in df_m.columns]

    # 主键对齐 — 病毒表用 Sample，元数据表已有 Run/query_id
    if 'Run' not in df_v.columns:
        for c in df_v.columns:
            if c.lower() in ['sample', 'run', 'run_accession', 'query_id', 'srr']:
                df_v.rename(columns={c: 'Run'}, inplace=True); break
        else:
            df_v.rename(columns={df_v.columns[0]: 'Run'}, inplace=True)

    if 'Run' not in df_m.columns:
        for c in df_m.columns:
            if c.lower() in ['run', 'run_accession', 'query_id', 'srr']:
                df_m.rename(columns={c: 'Run'}, inplace=True); break
        else:
            df_m.rename(columns={df_m.columns[0]: 'Run'}, inplace=True)

    meta_features = ['ScientificName', 'BioProject', 'CenterName', 'Tissue', 'Source',
                     'Location', 'Age_GrowthStage', 'CollectionDate', 'LibrarySource', 'TaxID']

    for col in meta_features:
        if col in df_m.columns:
            df_m[col] = df_m[col].apply(clean_ai_label)
        else:
            df_m[col] = 'Unknown'

    # 合并
    merged = pd.merge(df_v, df_m, on='Run', how='left')
    for col in meta_features:
        if col not in merged.columns: merged[col] = 'Unknown'
    merged['ScientificName'] = merged['ScientificName'].fillna('Unknown_Host')

    matched = (merged['ScientificName'] != 'Unknown_Host').sum()
    print(f"✅ 合并完成: {len(merged)} 条 (匹配元数据: {matched})")

    # --- 图1: 全局共感染概览 ---
    print("📊 模块 1/4: 全局共感染概览...")
    d1 = os.path.join(args.output, "01_Global_Summary"); os.makedirs(d1, exist_ok=True)

    # 获取 species/taxonomy 列
    sp_col = next((c for c in ['Taxonomy', 'taxonomy', 'Adjusted_Species', 'Species'] if c in merged.columns),
                  merged.columns[1])
    host_counts = merged.groupby(['Run', 'ScientificName'])[sp_col].nunique().reset_index()
    host_counts.columns = ['Run', 'ScientificName', 'Virus_Count']
    host_counts['Count_Label'] = host_counts['Virus_Count'].astype(str) + ' Virus(es)'

    ct = pd.crosstab(host_counts['ScientificName'], host_counts['Count_Label'])
    fig, ax = plt.subplots(figsize=(12, 7))
    ct.plot(kind='bar', stacked=True, colormap='Set2', ax=ax, edgecolor='black', linewidth=0.5)
    ax.set_title('Global Co-infection Overview across Hosts', fontsize=14, fontweight='bold')
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')
    ax.legend(title='Concurrent Infections', bbox_to_anchor=(1.05, 1))
    plt.tight_layout()
    plt.savefig(os.path.join(d1, "01_Coinfection_by_Host.pdf"), dpi=300, bbox_inches='tight'); plt.close()

    # --- 图2: Top 15 病毒组合 ---
    combos = merged.groupby('Run')[sp_col].apply(lambda x: ' + \n'.join(sorted(x))).reset_index()
    combo_counts = combos[sp_col].value_counts().head(15)
    fig, ax = plt.subplots(figsize=(12, 8))
    colors = plt.cm.Spectral(np.linspace(0, 1, len(combo_counts)))
    ax.barh(range(len(combo_counts)), combo_counts.values, color=colors)
    ax.set_yticks(range(len(combo_counts)))
    ax.set_yticklabels(combo_counts.index, fontsize=10)
    for i, v in enumerate(combo_counts.values):
        ax.text(v + combo_counts.max() * 0.01, i, str(v), va='center', fontweight='bold')
    ax.set_title('Top 15 Viral Infection Combinations', fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(os.path.join(d1, "02_Top_Combinations.pdf"), dpi=300, bbox_inches='tight'); plt.close()

    # --- 图3: 全局特征 vs 病毒谱 ---
    print("📊 模块 2/4: 全局特征 vs 病毒谱...")
    d2 = os.path.join(args.output, "02_Features_vs_Viruses"); os.makedirs(d2, exist_ok=True)
    top_viruses = merged[sp_col].value_counts().head(12).index
    df_filt = merged[merged[sp_col].isin(top_viruses)]

    for feat in meta_features:
        if feat not in df_filt.columns: continue
        top_feats = df_filt[feat].value_counts().head(8).index
        df_p = df_filt.copy()
        df_p[feat] = df_p[feat].apply(lambda x: x if x in top_feats else 'Other')
        ct = pd.crosstab(df_p[sp_col], df_p[feat])
        if ct.empty: continue
        ct['Total'] = ct.sum(1); ct = ct.sort_values('Total', ascending=True).drop(columns='Total')
        fig, ax = plt.subplots(figsize=(14, max(6, len(ct) * 0.5)))
        ct.plot(kind='barh', stacked=True, ax=ax, colormap='tab20', edgecolor='black', linewidth=0.5)
        ax.set_title(f'Viral Taxonomy across {feat}', fontsize=14, fontweight='bold')
        ax.set_xlabel('Detections'); ax.legend(title=feat, bbox_to_anchor=(1.05, 1))
        plt.tight_layout()
        plt.savefig(os.path.join(d2, f"Virus_vs_{safe_name(feat)}.pdf"), dpi=300, bbox_inches='tight'); plt.close()

    # --- 图4: 感染复杂度拆解 (1种/2种/3+) ---
    print("📊 模块 3/4: 感染复杂度拆解...")
    d3 = os.path.join(args.output, "03_Complexity_Breakdown"); os.makedirs(d3, exist_ok=True)

    sample_info = merged.groupby('Run').agg({
        sp_col: lambda x: sorted(set(x)),
        'ScientificName': 'first'
    }).reset_index()
    sample_info['Complexity'] = sample_info[sp_col].apply(len)
    sample_info['Combo'] = sample_info[sp_col].apply(lambda x: ' + \n'.join(x))
    sample_info['Level'] = sample_info['Complexity'].apply(lambda x: x if x <= 2 else '3+')

    for level in [1, 2, '3+']:
        df_lvl = sample_info[sample_info['Level'] == level]
        if df_lvl.empty: continue
        top_combos = df_lvl['Combo'].value_counts().head(12).index
        df_lt = df_lvl[df_lvl['Combo'].isin(top_combos)]
        ct = pd.crosstab(df_lt['Combo'], df_lt['ScientificName'])
        ct['Total'] = ct.sum(1); ct = ct.sort_values('Total', ascending=True).drop(columns='Total')
        if ct.empty: continue
        fig, ax = plt.subplots(figsize=(14, max(5, len(ct) * 0.6)))
        ct.plot(kind='barh', stacked=True, ax=ax, colormap='tab20', edgecolor='black', linewidth=0.5)
        title = "Single Infections" if level == 1 else f"{level}-Virus Co-infections"
        ax.set_title(f'{title} Breakdown by Host', fontsize=14, fontweight='bold')
        ax.legend(title='Host', bbox_to_anchor=(1.05, 1))
        plt.tight_layout()
        plt.savefig(os.path.join(d3, f"Level_{level}_Breakdown.pdf"), dpi=300, bbox_inches='tight'); plt.close()

    # --- 图5: Top 8 病毒专属画像 ---
    print("📊 模块 4/4: 病毒专属特写画像...")
    d4 = os.path.join(args.output, "04_Virus_Profiles"); os.makedirs(d4, exist_ok=True)
    top8 = merged[sp_col].value_counts().head(8).index

    for virus in top8:
        vdir = os.path.join(d4, f"Profile_{safe_name(virus)}"); os.makedirs(vdir, exist_ok=True)
        dv = merged[merged[sp_col] == virus]
        for feat in meta_features:
            if feat not in dv.columns or feat == 'ScientificName': continue
            top_f = dv[feat].value_counts().head(6).index
            dv_f = dv.copy(); dv_f[feat] = dv_f[feat].apply(lambda x: x if x in top_f else 'Other')
            ct = pd.crosstab(dv_f[feat], dv_f['ScientificName'])
            if ct.empty: continue
            ct['Total'] = ct.sum(1); ct = ct.sort_values('Total', ascending=False).drop(columns='Total')
            fig, ax = plt.subplots(figsize=(10, 6))
            ct.plot(kind='bar', stacked=True, ax=ax, colormap='Set1', edgecolor='black', linewidth=0.5)
            ax.set_title(f"'{virus[:40]}' across {feat}", fontsize=12, fontweight='bold')
            ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')
            ax.legend(title='Host', bbox_to_anchor=(1.05, 1))
            plt.tight_layout()
            plt.savefig(os.path.join(vdir, f"vs_{safe_name(feat)}.pdf"), dpi=300, bbox_inches='tight'); plt.close()

    # 保存合并表
    merged.to_csv(os.path.join(args.output, "Viral_Infection_Metadata_Merged.tsv"), sep='\t', index=False)
    print("\n🎉 元数据关联图全部完成！")


# =====================================================================
# 主入口
# =====================================================================

def main():
    parser = argparse.ArgumentParser(
        description='病毒组学三合一可视化引擎: depth / freq / meta / all',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    # ── 模式选择 ──
    g = parser.add_argument_group('模式选择')
    g.add_argument('--mode', default='all', choices=['depth', 'freq', 'meta', 'all'],
                   help='depth=测序深度图 | freq=丰度箱线图 | meta=元数据关联 | all=全部 (default: depth)')

    # ── 通用参数 ──
    g = parser.add_argument_group('通用参数')
    g.add_argument('-m', '--summary', required=True, help='all_viruses.best.summary.tsv 路径')
    g.add_argument('-o', '--output', required=True, help='输出目录')

    # ── 深度图参数 ──
    g = parser.add_argument_group('深度图 (--mode depth)')
    g.add_argument('-d', '--depth_dir', help='pandepth 深度文件目录 (.SiteDepth.gz)')
    g.add_argument('-n', '--sample', help='指定单个样本 (缺省全部)')
    g.add_argument('-g', '--add_genes', action='store_true', help='叠加 GenBank 基因轨道')
    g.add_argument('--gbk_dir', default='gbk_files', help='GenBank 缓存目录')
    g.add_argument('-w', '--window', type=int, default=100, help='平滑窗口 (default: 100)')
    g.add_argument('-f', '--fontsize', type=float, default=9, help='注释字体')
    g.add_argument('-t', '--threads', type=int, default=0, help='并发数 (default: auto)')

    # ── 丰度图参数 ──
    g = parser.add_argument_group('丰度图 (--mode freq)')
    g.add_argument('--freq_metrics', default=None,
                   help='指标列表逗号分隔 (default: MeanDepth,FPKM,RPM,TPM)')
    g.add_argument('--log10', action='store_true', help='log10 坐标')
    g.add_argument('--multi_plot', action='store_true', help='生成综合多指标面板图')
    g.add_argument('--width', type=float, default=12, help='图宽 (default: 12)')
    g.add_argument('--height', type=float, default=8, help='图高 (default: 8)')
    g.add_argument('--point_size', type=float, default=3, help='散点大小')
    g.add_argument('--dpi', type=int, default=300, help='分辨率')
    g.add_argument('--format', default='pdf', choices=['pdf', 'png'], help='输出格式')

    # ── 元数据图参数 ──
    g = parser.add_argument_group('元数据图 (--mode meta)')
    g.add_argument('--meta', default=None, help='元数据 TSV 路径')
    g.add_argument('--sra_list', default=None, help='SRA 列表 (自动生成元数据)')
    g.add_argument('--meta_outdir', default=None, help='元数据输出目录')
    g.add_argument('--deepseek_api', default=None, help='DeepSeek API Key')
    g.add_argument('--ncbi_api', default=None, help='NCBI API Key')
    g.add_argument('--email', default=None, help='NCBI 邮箱')
    g.add_argument('--meta_threads', type=int, default=4, help='元数据解析线程')

    args = parser.parse_args()

    modes = ['depth', 'freq', 'meta'] if args.mode == 'all' else [args.mode]

    for mode in modes:
        if mode == 'depth':
            if not args.depth_dir:
                print("⚠️ --depth_dir 未指定，跳过深度图模式")
                continue
            run_depth_mode(args)

        elif mode == 'freq':
            run_frequency_mode(args)

        elif mode == 'meta':
            run_metadata_mode(args)

    print("\n✅ 全部可视化任务完成！")


if __name__ == "__main__":
    main()
