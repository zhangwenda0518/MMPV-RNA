#!/usr/bin/env python3
"""
batch_plot_virus_depth.py — 病毒组学四合一可视化引擎
=================================================================
模式:
  --mode all      一键全部生成 (默认)
  --mode depth    测序深度图 (per-virus coverage + GenBank 基因轨道)
  --mode freq     丰度分布箱线图 (MeanDepth/FPKM/RPM/TPM)
  --mode sample   样本病毒数量分布 + 病毒发生率表
  --mode coabundance  病毒共丰度关联 (丰度热图 + 相关热图 + 网络)

深度模式 (--mode depth):
  python batch_plot_virus_depth.py --mode depth \\
    -d stat/ -m all_viruses.best.summary.tsv -o plots/ -t 8 -g

丰度模式 (--mode freq):
  python batch_plot_virus_depth.py --mode freq \\
    -m all_viruses.best.summary.tsv -o plots/ --log10

样本分布 (--mode sample):
  python batch_plot_virus_depth.py --mode sample \\
    -m all_viruses.best.summary.tsv -o plots/

共丰度分析 (--mode coabundance):
  python batch_plot_virus_depth.py --mode coabundance \\
    -m all_viruses.best.summary.tsv -o plots/

  产出: TPM丰度热图 + 相对丰度堆叠柱形 + Spearman相关热图 + 显著性网络
"""

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.ticker as ticker
from matplotlib.lines import Line2D

try:
    import seaborn as sns
    HAS_SEABORN = True
except ImportError:
    HAS_SEABORN = False

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
from pathlib import Path
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from itertools import combinations

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

def _resolve_taxonomy_col(df):
    """Return the best taxonomy/species column name from summary df."""
    for c in ['taxonomy', 'Adjusted_Species', 'Species_NCBI', 'Species_ICTV', 'Species']:
        if c in df.columns: return c
    return None

def _resolve_acc_col(df):
    """Return the best accession column name from summary df."""
    for c in ['Virus', 'Rep_Accession', 'Accession']:
        if c in df.columns: return c
    return None

def _make_display_name(row, tax_col, acc_col):
    """Build taxonomy\\n(accession) display label."""
    tax = str(row.get(tax_col, '')) if tax_col else ''
    acc = str(row.get(acc_col, '')) if acc_col else str(row.iloc[0])
    if tax and tax not in ['Unannotated', '-', 'nan', 'None', '']:
        return f"{tax}\n({acc})"
    return acc

def _load_summary(summary_path):
    """Load and normalize summary TSV/CSV."""
    df = pd.read_csv(summary_path, sep=',' if summary_path.endswith('.csv') else '\t')
    df.columns = [c.strip().replace('\ufeff', '') for c in df.columns]
    # Column alias normalization — only rename first matching source per target
    alias_rules = [
        (['Rep_Accession', 'Accession'], 'Virus'),
        (['Adjusted_Species', 'Species_NCBI', 'Species_ICTV', 'Species'], 'taxonomy'),
        (['Rep_MeanDepth'], 'MeanDepth'),
        (['Rep_Coverage(%)'], 'Coverage(%)'),
        (['Asm_TPM'], 'TPM'),
        (['Asm_FPKM'], 'FPKM'),
        (['Asm_RPM'], 'RPM'),
        (['Asm_EM_Reads'], 'MappedReads'),
        (['Rep_Length'], 'Length'),
    ]
    for sources, target in alias_rules:
        for src in sources:
            if src in df.columns and src != target:
                df = df.rename(columns={src: target})
                break
    return df


# =====================================================================
# 模块 1: 测序深度图
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
            depth_df = pd.read_csv(io.StringIO(result.stdout), sep='\t', header=None,
                                   names=['chr', 'position', 'depth'])
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
            png_path = os.path.join(sample_out_dir, f"{sample}_{safe_tax}_{safe_vname}_depth.png")

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
            ax.axhline(y=mean_depth, color='#d62728', linestyle='--', linewidth=2,
                       label=f'Mean Depth: {mean_depth:.2f}x')
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
                colors = ['#8dd3c7', '#ffffb3', '#bebada', '#fb8072', '#80b1d3', '#fdb462',
                          '#b3de69', '#fccde5', '#d9d9d9']
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
            plt.savefig(png_path, bbox_inches='tight', dpi=300)
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
    print(f" 模式: 测序深度图 (样本并发数: {max_workers})")
    print(f" 底层加速: rapidgzip [{'OK' if shutil.which('rapidgzip') else '--'}]  "
          f"ripgrep [{'OK' if shutil.which('rg') else '--'}]")
    print("=" * 60)

    summary_df = _load_summary(args.summary)
    samples_to_process = [args.sample] if args.sample else summary_df['Sample'].unique().tolist()

    genes_dict = {}
    if args.add_genes:
        print(f"\n 准备 GenBank 注释...")
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

    print(f"\n 下发智能解析任务...")
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
            print(f" 找不到样本 {sample} 的深度文件，跳过...")
            continue
        sample_out_dir = os.path.join(args.output, str(sample))
        os.makedirs(sample_out_dir, exist_ok=True)
        task = (sample, sample_summary.to_dict('records'), depth_file, sample_out_dir,
                args.window, args.fontsize, genes_dict, threads_per_job)
        future_to_sample[executor.submit(process_single_sample, task)] = sample

    success_plots, total_plots = 0, 0
    with tqdm(total=len(future_to_sample), desc=" 深度图进度", unit="sample") as pbar:
        for future in as_completed(future_to_sample):
            sample = future_to_sample[future]
            try:
                s_name, s_success, s_msg, results_log = future.result()
                if not s_success:
                    tqdm.write(f"  [{sample} 失败] {s_msg}")
                else:
                    for virus, v_success, info in results_log:
                        total_plots += 1
                        if v_success:
                            success_plots += 1
                            tqdm.write(f"    [{sample}] {os.path.basename(info)}")
                        else:
                            tqdm.write(f"    [{sample}] {virus}: {info}")
            except Exception as exc:
                tqdm.write(f"  [{sample}] {exc}")
            pbar.update(1)
    executor.shutdown()
    print(f"\n 深度图完成! (成功: {success_plots} / {total_plots})")


# =====================================================================
# 模块 2: 丰度分布箱线图
# =====================================================================

def run_frequency_mode(args):
    print("=" * 60)
    print(f" 模式: 病毒丰度分布箱线图")
    print("=" * 60)

    df = _load_summary(args.summary)
    tax_col = _resolve_taxonomy_col(df)
    acc_col = _resolve_acc_col(df)

    df['Display_Name'] = df.apply(lambda r: _make_display_name(r, tax_col, acc_col), axis=1)
    df['Display_Name'] = df['Display_Name'].apply(
        lambda x: '\n'.join(textwrap.wrap(str(x), width=40))
    )

    # Append virus prevalence: n=X/Y (Z%)
    total_samples = df['Sample'].nunique()
    prevalence = df.groupby('Display_Name')['Sample'].nunique()
    df['Display_Name'] = df['Display_Name'].apply(
        lambda x: f"{x}\nn={prevalence[x]}/{total_samples} ({prevalence[x]/total_samples*100:.0f}%)"
    )

    if args.freq_metrics:
        metrics = [m.strip() for m in args.freq_metrics.split(',')]
    else:
        metrics = ['MeanDepth', 'FPKM', 'RPM', 'TPM']

    available = [m for m in metrics if m in df.columns]
    if not available:
        print(f" 未找到指定的丰度列。可用列: {list(df.columns)}")
        sys.exit(1)

    os.makedirs(args.output, exist_ok=True)
    unique_count = df['Display_Name'].nunique()
    should_flip = unique_count > 5

    # -- seaborn style matching R theme_bw --
    if HAS_SEABORN:
        sns.set_style("ticks")
        sns.set_context("paper", font_scale=1.2)

    turbo_palette = (sns.color_palette("turbo", max(unique_count, 1))
                     if HAS_SEABORN else
                     plt.cm.turbo(np.linspace(0, 1, max(unique_count, 1))))

    def make_boxplot(data, value_col, out_prefix):
        plot_df = data.dropna(subset=[value_col]).copy()
        if plot_df.empty:
            return

        medians = plot_df.groupby('Display_Name')[value_col].median().sort_values()
        plot_df['Display_Name'] = pd.Categorical(
            plot_df['Display_Name'], categories=medians.index, ordered=True)

        fig, ax = plt.subplots(figsize=(args.width, max(8, unique_count * 0.6)))

        if HAS_SEABORN:
            if should_flip:
                sns.boxplot(data=plot_df, y='Display_Name', x=value_col, hue='Display_Name', ax=ax,
                           palette=turbo_palette, legend=False, width=0.6, linewidth=1.0,
                           fliersize=0, saturation=0.85)
                sns.stripplot(data=plot_df, y='Display_Name', x=value_col, hue='Display_Name', ax=ax,
                             palette=turbo_palette, legend=False, size=3.5, alpha=0.7,
                             jitter=0.2, edgecolor='white', linewidth=0.3)
            else:
                sns.boxplot(data=plot_df, x='Display_Name', y=value_col, hue='Display_Name', ax=ax,
                           palette=turbo_palette, legend=False, width=0.6, linewidth=1.0,
                           fliersize=0, saturation=0.85)
                sns.stripplot(data=plot_df, x='Display_Name', y=value_col, hue='Display_Name', ax=ax,
                             palette=turbo_palette, legend=False, size=3.5, alpha=0.7,
                             jitter=0.2, edgecolor='white', linewidth=0.3)
        else:
            positions = list(range(len(medians)))
            bp = ax.boxplot(
                [plot_df[plot_df['Display_Name'] == n][value_col].values
                 for n in medians.index],
                positions=positions, patch_artist=True, widths=0.6,
                showfliers=False, vert=not should_flip)
            for i, (patch, name) in enumerate(zip(bp['boxes'], medians.index)):
                patch.set_facecolor(turbo_palette[i % len(turbo_palette)])
                patch.set_alpha(0.6)
            for i, name in enumerate(medians.index):
                vals = plot_df[plot_df['Display_Name'] == name][value_col].values
                jitter = np.random.uniform(-0.2, 0.2, len(vals))
                if should_flip:
                    ax.scatter(vals, np.full(len(vals), i) + jitter, alpha=0.7,
                              s=args.point_size ** 2,
                              color=turbo_palette[i % len(turbo_palette)],
                              edgecolors='white', linewidth=0.3)
                else:
                    ax.scatter(np.full(len(vals), i) + jitter, vals, alpha=0.7,
                              s=args.point_size ** 2,
                              color=turbo_palette[i % len(turbo_palette)],
                              edgecolors='white', linewidth=0.3)

        if not should_flip:
            ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right',
                              fontsize=10, style='italic')

        if args.log10:
            if should_flip:
                ax.set_xscale('log')
            else:
                ax.set_yscale('log')

        ax.set_title(f'Virus Abundance — {value_col} Distribution',
                    fontweight='bold', fontsize=16, pad=15)
        if should_flip:
            ax.set_xlabel(value_col, fontweight='bold', fontsize=14)
            ax.set_ylabel('')
        else:
            ax.set_ylabel(value_col, fontweight='bold', fontsize=14)
            ax.set_xlabel('')
        ax.grid(axis='y', linestyle=':', alpha=0.3)
        if HAS_SEABORN:
            sns.despine(ax=ax)
        plt.tight_layout()
        fn = os.path.join(args.output,
                         f"{out_prefix}_{value_col.lower()}"
                         f"{'_log10' if args.log10 else ''}.{args.format}")
        plt.savefig(fn, dpi=args.dpi, bbox_inches='tight')
        plt.savefig(fn.replace(".pdf",".png"), dpi=args.dpi, bbox_inches='tight')
        plt.close()
        print(f"   {os.path.basename(fn)}")

    for m in available:
        make_boxplot(df, m, "freq")

    # -- multi-metric facet panel --
    if args.multi_plot and len(available) > 1:
        print(f"   生成综合多指标面板图...")
        plot_data = df.melt(id_vars=['Display_Name'], value_vars=available,
                            var_name='Metric', value_name='Value').dropna()
        base_metric = available[0]
        medians = (plot_data[plot_data['Metric'] == base_metric]
                   .groupby('Display_Name')['Value'].median().sort_values())
        plot_data['Display_Name'] = pd.Categorical(
            plot_data['Display_Name'], categories=medians.index, ordered=True)

        n_cols = len(available)
        fig, axes = plt.subplots(1, n_cols,
                                figsize=(args.width * 1.8,
                                         max(args.height, unique_count * 0.65)))
        if n_cols == 1:
            axes = [axes]

        for i, (ax, m) in enumerate(zip(axes, available)):
            mdf = plot_data[plot_data['Metric'] == m]
            if HAS_SEABORN:
                sns.boxplot(data=mdf, y='Display_Name', x='Value', hue='Display_Name', ax=ax,
                           palette=turbo_palette, legend=False, width=0.6, linewidth=1.0,
                           fliersize=0, saturation=0.85)
                sns.stripplot(data=mdf, y='Display_Name', x='Value', hue='Display_Name', ax=ax,
                             palette=turbo_palette, legend=False, size=2.5, alpha=0.55,
                             jitter=0.2, edgecolor='white', linewidth=0.2)
            else:
                positions = list(range(len(medians)))
                bp = ax.boxplot(
                    [mdf[mdf['Display_Name'] == n]['Value'].values
                     for n in medians.index],
                    positions=positions, patch_artist=True, widths=0.6,
                    showfliers=False, vert=False)
                for pi, (patch, _) in enumerate(zip(bp['boxes'], medians.index)):
                    patch.set_facecolor(turbo_palette[pi % len(turbo_palette)])
                    patch.set_alpha(0.6)
                for pi, name in enumerate(medians.index):
                    vals = mdf[mdf['Display_Name'] == name]['Value'].values
                    ax.scatter(vals,
                              np.full(len(vals), pi) +
                              np.random.uniform(-0.15, 0.15, len(vals)),
                              alpha=0.5, s=args.point_size,
                              color=turbo_palette[pi % len(turbo_palette)],
                              edgecolors='white', linewidth=0.2)
            ax.set_title(m, fontweight='bold', fontsize=13)
            if args.log10:
                ax.set_xscale('log')
            if HAS_SEABORN:
                sns.despine(ax=ax)

        # Hide Y-axis labels on panels 2..N, keep only leftmost
        for i in range(1, n_cols):
            axes[i].set_ylabel('')
            axes[i].tick_params(left=False, labelleft=False)

        fig.suptitle('Comprehensive Multi-metric Virus Abundance',
                    fontweight='bold', fontsize=18, y=1.01)
        plt.tight_layout()
        fn = os.path.join(args.output,
                         f"freq_multi_metrics"
                         f"{'_log10' if args.log10 else ''}.{args.format}")
        plt.savefig(fn, dpi=args.dpi, bbox_inches='tight')
        plt.savefig(fn.replace(".pdf",".png"), dpi=args.dpi, bbox_inches='tight')
        plt.close()
        print(f"   {os.path.basename(fn)}")

    print("\n 丰度分布图完成!")

# =====================================================================
# 模块 3: 样本病毒数量分布 + 病毒发生率
# =====================================================================

def run_sample_mode(args):
    print("=" * 60)
    print(" 模式: 样本病毒数量分布 + 病毒发生率")
    print("=" * 60)

    df = _load_summary(args.summary)
    tax_col = _resolve_taxonomy_col(df)
    acc_col = _resolve_acc_col(df)

    os.makedirs(args.output, exist_ok=True)

    # --- 每样本病毒数 ---
    sample_counts = df.groupby('Sample').apply(
        lambda g: g[[c for c in [tax_col, acc_col] if c]].drop_duplicates().shape[0]
    ).reset_index(name='Virus_Count')
    sample_counts = sample_counts.sort_values('Virus_Count', ascending=False)

    # 柱形图
    fig, ax = plt.subplots(figsize=(max(10, len(sample_counts) * 0.3), 6))
    colors = plt.cm.viridis(np.linspace(0.2, 0.9, len(sample_counts)))
    bars = ax.bar(range(len(sample_counts)), sample_counts['Virus_Count'], color=colors, edgecolor='#333', linewidth=0.5)
    for i, (_, row) in enumerate(sample_counts.iterrows()):
        ax.text(i, row['Virus_Count'] + 0.1, str(row['Virus_Count']), ha='center', fontsize=8, fontweight='bold')
    ax.set_xticks(range(len(sample_counts)))
    ax.set_xticklabels(sample_counts['Sample'], rotation=90, ha='center', fontsize=7)
    ax.set_ylabel('Virus Species Count', fontweight='bold')
    ax.set_title(f'Detected Virus Count per Sample (n={len(sample_counts)})', fontweight='bold', fontsize=14)
    ax.set_ylim(bottom=0)
    plt.tight_layout()
    fn_bar = os.path.join(args.output, "sample_virus_count_bar.pdf")
    fn_bar_png = os.path.join(args.output, "sample_virus_count_bar.png")
    plt.savefig(fn_bar, dpi=300, bbox_inches='tight')
    plt.savefig(fn_bar_png, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"   {os.path.basename(fn_bar)}")

    # 饼图 (1种 / 2种 / 3+种)
    def categorize(n):
        if n == 1: return '1 virus'
        if n == 2: return '2 viruses'
        return '3+ viruses'

    sample_counts['Category'] = sample_counts['Virus_Count'].apply(categorize)
    pie_data = sample_counts['Category'].value_counts()
    pie_order = ['1 virus', '2 viruses', '3+ viruses']
    pie_vals = [pie_data.get(k, 0) for k in pie_order]
    pie_colors = ['#66c2a5', '#fc8d62', '#8da0cb']

    fig, ax = plt.subplots(figsize=(7, 7))
    wedges, texts, autotexts = ax.pie(pie_vals, labels=pie_order, autopct='%1.1f%%',
                                       colors=pie_colors, startangle=90,
                                       textprops={'fontsize': 12, 'fontweight': 'bold'})
    for at in autotexts: at.set_fontsize(13)
    ax.set_title('Sample Co-infection Complexity', fontweight='bold', fontsize=14)
    fn_pie = os.path.join(args.output, "sample_virus_count_pie.pdf")
    fn_pie_png = os.path.join(args.output, "sample_virus_count_pie.png")
    plt.savefig(fn_pie, dpi=300, bbox_inches='tight')
    plt.savefig(fn_pie_png, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"   {os.path.basename(fn_pie)}")

    # 数据表
    sample_counts[['Sample', 'Virus_Count', 'Category']].to_csv(
        os.path.join(args.output, "sample_virus_count.csv"), index=False)

    # --- 每病毒发生率 ---
    virus_df = df[[c for c in [tax_col, acc_col, 'Sample'] if c]].drop_duplicates()
    virus_occ = virus_df.groupby([c for c in [tax_col, acc_col] if c]).agg(
        Sample_Count=('Sample', 'nunique')
    ).reset_index()
    total_samples = df['Sample'].nunique()
    virus_occ['Prevalence(%)'] = (virus_occ['Sample_Count'] / total_samples * 100).round(1)

    # 构建显示名
    tax_c = tax_col if tax_col in virus_occ.columns else None
    acc_c = acc_col if acc_col in virus_occ.columns else None
    virus_occ['Display'] = virus_occ.apply(lambda r: _make_display_name(r, tax_c, acc_c), axis=1)
    virus_occ = virus_occ.sort_values('Sample_Count', ascending=True)

    # 横向柱形图
    fig, ax = plt.subplots(figsize=(10, max(5, len(virus_occ) * 0.35)))
    colors = plt.cm.plasma(np.linspace(0.1, 0.9, len(virus_occ)))
    bars = ax.barh(range(len(virus_occ)), virus_occ['Sample_Count'], color=colors, edgecolor='#333', linewidth=0.5)
    for i, (_, row) in enumerate(virus_occ.iterrows()):
        ax.text(row['Sample_Count'] + max(virus_occ['Sample_Count']) * 0.01, i,
                f"{row['Sample_Count']} ({row['Prevalence(%)']}%)",
                va='center', fontsize=8, fontweight='bold')
    ax.set_yticks(range(len(virus_occ)))
    ax.set_yticklabels(virus_occ['Display'], fontsize=8, style='italic')
    ax.set_xlabel('Number of Samples', fontweight='bold')
    ax.set_title(f'Virus Occurrence across {total_samples} Samples', fontweight='bold', fontsize=14)
    plt.tight_layout()
    fn_occ = os.path.join(args.output, "virus_occurrence_bar.pdf")
    fn_occ_png = os.path.join(args.output, "virus_occurrence_bar.png")
    plt.savefig(fn_occ, dpi=300, bbox_inches='tight')
    plt.savefig(fn_occ_png, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"   {os.path.basename(fn_occ)}")

    # 表
    virus_occ_out = virus_occ.drop(columns=['Display'])
    virus_occ_out.to_csv(os.path.join(args.output, "virus_occurrence.csv"), index=False)

    print(f"\n 样本分布图完成! (样本={total_samples}, 病毒={len(virus_occ)})")


# =====================================================================
# 模块 4: 病毒共丰度关联 (丰度热图 + 相对丰度 + Spearman)
# =====================================================================

def run_coabundance_mode(args):
    print("=" * 60)
    print(" 模式: 病毒共丰度关联网络")
    print("=" * 60)

    df = _load_summary(args.summary)
    tax_col = _resolve_taxonomy_col(df)
    acc_col = _resolve_acc_col(df)

    # 构建样本×病毒表 — add Virus_Label to df
    label_key = list({c for c in [tax_col, acc_col] if c})
    label_df = df[label_key].drop_duplicates()
    label_df['Virus_Label'] = label_df.apply(lambda r: _make_display_name(r, tax_col, acc_col).replace('\n', ' '), axis=1)
    df = df.merge(label_df, on=label_key, how='left')

    # TPM 矩阵
    tpm_pivot = df.pivot_table(index='Sample', columns='Virus_Label', values='TPM', aggfunc='mean').fillna(0)
    # 二进制矩阵 (用于 prevalence 过滤)
    bin_pivot = (tpm_pivot > 0).astype(int)

    viruses = tpm_pivot.columns.tolist()
    n_viruses = len(viruses)

    if n_viruses < 2:
        print(f"  仅 {n_viruses} 种病毒, 无法做共现分析")
        return

    # 过滤: 至少出现在 min_occurrence 个样本中
    min_occ = getattr(args, 'coab_min_occurrence', 2)
    virus_prevalence = bin_pivot.sum(axis=0)
    keep = virus_prevalence[virus_prevalence >= min_occ].index.tolist()
    if len(keep) < 2:
        print(f"  过滤后 ({min_occ}+ 样本) 不足 2 种病毒, 跳过")
        return
    bin_pivot = bin_pivot[keep]
    tpm_pivot = tpm_pivot[keep]
    viruses = keep
    n_viruses = len(viruses)
    total_samples = bin_pivot.shape[0]
    print(f"  病毒数: {n_viruses} (出现 >= {min_occ} 样本), 样本数: {total_samples}")

    os.makedirs(args.output, exist_ok=True)

    # ── 丰度热图 (samples × viruses, log10 TPM+1) ──
    tpm_display = tpm_pivot.copy()
    tpm_display.columns = [c.replace('\n', ' ').split('(')[-1].rstrip(')') if '(' in c else c[:15]
                           for c in tpm_display.columns]
    log_tpm = np.log10(tpm_display + 1)
    fig, ax = plt.subplots(figsize=(max(8, n_viruses * 1.2), max(5, total_samples * 0.6)))
    sns.heatmap(log_tpm, annot=True, fmt=".1f", cmap="YlGnBu", ax=ax,
               cbar_kws={'label': 'log10(TPM + 1)'})
    ax.set_title(f'Virus Abundance Heatmap\n({n_viruses} viruses, {total_samples} samples)',
                fontweight='bold', fontsize=13)
    ax.set_ylabel('Sample')
    ax.set_xlabel('')
    plt.xticks(rotation=45, ha='right', fontsize=8)
    plt.yticks(rotation=0, fontsize=9)
    plt.tight_layout()
    fn_ah = os.path.join(args.output, f"coabundance_abundance_heatmap.{args.format}")
    plt.savefig(fn_ah, dpi=args.dpi, bbox_inches='tight')
    plt.savefig(fn_ah.replace(".pdf",".png"), dpi=args.dpi, bbox_inches='tight'); plt.close()
    print(f"   {os.path.basename(fn_ah)}")

    # ── 相对丰度堆叠柱形图 ──
    if 'TPM' in df.columns:
        rel_pivot = df.pivot_table(index='Sample', columns='Virus_Label', values='TPM',
                                   aggfunc='mean').fillna(0)
        rel_pivot = rel_pivot[keep]
        rel_pivot_pct = rel_pivot.div(rel_pivot.sum(axis=1), axis=0) * 100
        rel_pivot_pct.columns = [c.replace('\n', ' ').split('(')[-1].rstrip(')') if '(' in c else c[:15]
                                 for c in rel_pivot_pct.columns]
        fig, ax = plt.subplots(figsize=(max(10, total_samples * 0.8), 6))
        rel_pivot_pct.plot(kind='bar', stacked=True, cmap='Set2', ax=ax, edgecolor='#333', linewidth=0.5)
        ax.set_title(f'Virus Relative Abundance Composition per Sample\n({total_samples} samples)',
                    fontweight='bold', fontsize=13)
        ax.set_ylabel('Relative Abundance (%)', fontweight='bold')
        ax.set_xlabel('Sample')
        ax.legend(title='', bbox_to_anchor=(1.02, 1), loc='upper left', fontsize=8)
        plt.xticks(rotation=0)
        plt.tight_layout()
        fn_ra = os.path.join(args.output, f"coabundance_relabund_stacked.{args.format}")
        plt.savefig(fn_ra, dpi=args.dpi, bbox_inches='tight')
        plt.savefig(fn_ra.replace(".pdf",".png"), dpi=args.dpi, bbox_inches='tight'); plt.close()
        print(f"   {os.path.basename(fn_ra)}")

    # ── 描述性共现率矩阵 ──
    print("\n [Co-occurrence] 计算病毒两两共现率...")
    bin_mat = bin_pivot.values.T  # viruses × samples
    cooc_mat = np.full((n_viruses, n_viruses), np.nan)
    annot_mat = np.full((n_viruses, n_viruses), '', dtype=object)
    edges = []
    for i, j in combinations(range(n_viruses), 2):
        a = bin_mat[i]; b = bin_mat[j]
        n_both = int(np.sum(a * b))
        n_i_only = int(np.sum(a * (1 - b)))
        n_j_only = int(np.sum((1 - a) * b))
        n_neither = int(np.sum((1 - a) * (1 - b)))
        rate = n_both / total_samples * 100
        cooc_mat[i, j] = rate
        cooc_mat[j, i] = rate
        annot_mat[i, j] = f'{n_both}/{total_samples}\n{rate:.0f}%'
        annot_mat[j, i] = annot_mat[i, j]
        edges.append({
            'from': viruses[i], 'to': viruses[j],
            'both': n_both, 'only_A': n_i_only, 'only_B': n_j_only, 'neither': n_neither,
            'cooccur_rate(%)': round(rate, 1)
        })
    # 对角线 = 自身发生率
    for i in range(n_viruses):
        n_i = int(bin_mat[i].sum())
        rate_i = n_i / total_samples * 100
        cooc_mat[i, i] = rate_i
        annot_mat[i, i] = f'{n_i}/{total_samples}\n{rate_i:.0f}%'

    edges_df = pd.DataFrame(edges)
    edges_df.to_csv(os.path.join(args.output, "coabundance_cooccurrence.csv"), index=False)

    short_names = [v.replace('\n', ' ').split('(')[-1].rstrip(')') if '(' in v else v[:12]
                   for v in viruses]
    mask = np.triu(np.ones_like(cooc_mat, dtype=bool), k=1)
    fig, ax = plt.subplots(figsize=(max(7, n_viruses * 1.0), max(5, n_viruses * 0.8)))
    sns.heatmap(cooc_mat, annot=annot_mat, fmt='', cmap="YlOrRd", vmin=0, vmax=100,
               mask=mask, ax=ax, cbar_kws={'label': 'Co-occurrence Rate (%)', 'shrink': 0.8},
               xticklabels=short_names, yticklabels=short_names,
               linewidths=0.5, linecolor='white')
    ax.set_title(f'Virus Co-occurrence Rate\n({n_viruses} viruses, {total_samples} samples)',
                fontweight='bold', fontsize=12)
    plt.xticks(rotation=45, ha='right', fontsize=8)
    plt.yticks(rotation=0, fontsize=8)
    plt.tight_layout()
    fn_co = os.path.join(args.output, f"coabundance_cooccurrence_heatmap.{args.format}")
    plt.savefig(fn_co, dpi=args.dpi, bbox_inches='tight')
    plt.savefig(fn_co.replace(".pdf",".png"), dpi=args.dpi, bbox_inches='tight'); plt.close()
    print(f"   {os.path.basename(fn_co)}")

    # ── Spearman 丰度相关性热图 ──
    corr_mat = tpm_pivot.corr(method='spearman')
    sc_names = [v.replace('\n', ' ').split('(')[-1].rstrip(')') if '(' in v else v[:12]
                for v in corr_mat.columns]
    fig, ax = plt.subplots(figsize=(max(6, n_viruses * 0.8), max(5, n_viruses * 0.7)))
    smask = np.triu(np.ones_like(corr_mat, dtype=bool), k=1)
    sns.heatmap(corr_mat, annot=True, fmt=".2f", cmap="vlag", vmin=-1, vmax=1, center=0,
               mask=smask, ax=ax, cbar_kws={'label': 'Spearman ρ', 'shrink': 0.8},
               xticklabels=sc_names, yticklabels=sc_names)
    ax.set_title(f'Spearman Correlation between Virus Species\n({n_viruses} viruses, {total_samples} samples)',
                fontweight='bold', fontsize=12)
    plt.xticks(rotation=45, ha='right', fontsize=8)
    plt.yticks(rotation=0, fontsize=8)
    plt.tight_layout()
    fn_sc = os.path.join(args.output, f"coabundance_spearman_heatmap.{args.format}")
    plt.savefig(fn_sc, dpi=args.dpi, bbox_inches='tight')
    plt.savefig(fn_sc.replace(".pdf",".png"), dpi=args.dpi, bbox_inches='tight'); plt.close()
    print(f"   {os.path.basename(fn_sc)}")

    print("\n 共丰度分析完成!")

# =====================================================================
# 主入口
# =====================================================================

def main():
    parser = argparse.ArgumentParser(
        description='病毒组学四合一可视化引擎: depth / freq / sample / coabundance / all',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    # ── 模式选择 ──
    g = parser.add_argument_group('模式选择')
    g.add_argument('--mode', default='all',
                   choices=['depth', 'freq', 'sample', 'coabundance', 'all'],
                   help='depth=测序深度图 | freq=丰度箱线图 | sample=样本分布 | coabundance=共丰度网络 | all=全部 (default)')

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

    # ── 共丰度参数 ──
    g = parser.add_argument_group('共丰度 (--mode coabundance)')
    g.add_argument('--coab_min_occurrence', type=int, default=2,
                   help='病毒最少出现样本数 (default: 2)')

    args = parser.parse_args()

    all_modes = ['depth', 'freq', 'sample', 'coabundance']
    modes = all_modes if args.mode == 'all' else [args.mode]

    for mode in modes:
        if mode == 'depth':
            if not args.depth_dir:
                print("  --depth_dir 未指定，跳过深度图模式\n")
                continue
            run_depth_mode(args)

        elif mode == 'freq':
            run_frequency_mode(args)

        elif mode == 'sample':
            run_sample_mode(args)

        elif mode == 'coabundance':
            run_coabundance_mode(args)

    print("\n 全部可视化任务完成!")


if __name__ == "__main__":
    main()
