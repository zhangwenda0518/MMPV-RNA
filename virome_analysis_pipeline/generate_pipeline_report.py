#!/usr/bin/env python3
"""
generate_pipeline_report.py — Interactive HTML Report Generator
===============================================================
Scans pipeline output, generates an interactive HTML report with:
  - Left sidebar navigation
  - Embedded charts from post-hoc analysis
  - Summary data tables
  - AI interpretation prompts
"""

import argparse
import os
import sys
import base64
import shutil
from pathlib import Path
from datetime import datetime

# 跨管线统一 I/O 布局 (mmpv_common/, 仓库根): 目录名随 MMPV_IO_LAYOUT 解析
# (编排器已 normalize 环境变量, 子进程导入本模块时快照即正确布局)
_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)
from mmpv_common.io_layout import layout_dirs as _layout_dirs, locate_analysis_dir
_D = _layout_dirs(os.environ.get("MMPV_IO_LAYOUT"))

# --phylo-root (2026-10-15): ⑤ 系统发育管线输出根 (results_* 目录), main() 里赋值。
# capheine 正选择分析已迁 ⑤ (select 模块), 本管线新跑不会再产出 07_capheine;
# 提供此参数后 Stage 7 各触点额外扫描 <phylo_root>/**/select|05_select/capheine 产物。
PHYLO_ROOT = None

try:
    from mmpv_common.table_widgets import TW_BLOCK as _TW
except Exception:
    _TW = ""   # 共享表格控件不可达时静默降级为原生表格


def _esc(v):
    """HTML-escape a value to prevent XSS injection."""
    return (str(v).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
            .replace('"', "&quot;").replace("'", "&#39;"))


def _fnum(v, default=0.0):
    """容错 float: summary TSV 用 "-" 占位表示缺失, 直接 float() 会 ValueError 打断全报告。"""
    try:
        if v is None or (isinstance(v, str) and v.strip() in ("", "-", "NA", "N/A", "nan", ".")):
            return default
        return float(v)
    except (TypeError, ValueError):
        return default


_SKIP_IMAGES = False   # --no_images: 跳过图片字节内嵌(卡片仅留文件名占位)

# [P2-24b] 目录名尾部 accession 形态 (与 utils/extract_full_fasta.py 的 _ACC_RE 同步放宽):
#   NC_028268.1 / MW648525.1 (旧形态) + TA104.1 / U13908.2 / Y_U13908.2 / 数字段前缀 (新增)
_ACC_PAT = (r'(\d+[A-Z]{2}\d{5,6}\.\d+'
            r'|[A-Z]{1,2}_?\d{5,6}\.\d+'
            r'|[A-Z]{1,2}_\d+\.\d+'
            r'|[A-Z]{2,4}_?\d{5,9}\.\d{1,2})')


def _img_placeholder(title):
    """--no_images / 图片超限时的文件名占位卡 (不内嵌字节)。"""
    return (f'<div class="chart-card"><div class="chart-title">{_esc(title)}</div>'
            f'<div style="padding:10px;color:#888;font-size:11px">[image omitted]</div></div>')

try:
    import pandas as pd
    HAS_PANDAS = True
except ImportError:
    HAS_PANDAS = False

try:
    import numpy as _np
    _INT_TYPES = (int, _np.integer)   # Series.sum()/nunique() 返回 numpy.int64, isinstance(x, int) 为 False
except ImportError:
    _INT_TYPES = (int,)


def safe_read_csv(fp, sep="\t"):
    if not Path(fp).exists(): return None
    try:
        if HAS_PANDAS:
            df = pd.read_csv(fp, sep=sep)
            return df if len(df) > 0 else None
    except Exception:
        return None
    return None


def img_to_base64(path, max_kb=2000):
    """Convert image to base64 for embedding. Skip if > max_kb."""
    if not path or not Path(path).exists():
        return None
    if _SKIP_IMAGES:   # P3-19: --no_images 生效点 (所有内嵌都经此函数)
        return None
    size_kb = Path(path).stat().st_size / 1024
    if size_kb > max_kb:
        return None
    with open(path, 'rb') as f:
        return base64.b64encode(f.read()).decode('utf-8')


def collect_charts(post_dir, virus_acc):
    """Collect key chart paths for a virus (fuzzy match by accession in dir name)."""
    charts = {}
    vdir = Path(post_dir) / virus_acc
    if not Path(post_dir).is_dir():
        return charts
    if not vdir.exists():
        for d in Path(post_dir).iterdir():
            if d.is_dir() and virus_acc in d.name:
                vdir = d; break
        else:
            return charts

    chart_patterns = {
        'vcf_viz': [
            ('Figure1A_All_Variants_Landscape.png', 'Genome-wide Variant Landscape'),
            ('Figure1B_Top50_Variants_Landscape.png', 'Top50 Variant Landscape'),
            ('Figure2_TsTv_Pie.png', 'Ts/Tv Ratio'),
            ('Figure3_Functional_Pie.png', 'Functional Classification Pie'),
            ('Figure4_Clustermap.png', 'Variant Clustermap'),
            ('Figure5_AFS.png', 'Allele Frequency Spectrum'),
            ('Figure6_Variant_Density.png', 'Variant Density'),
            ('Figure7_AF_Violin.png', 'AF Violin'),
            ('Figure8_PopGen_Dynamics.png', 'PopGen Sliding Window'),
            ('Figure9_PCA_Lineages.png', 'PCA 2D Clustering'),
            ('Figure9_PCA_Lineages_3D.png', 'PCA 3D Clustering'),
        ],
        'snpeff_macro': [
            ('Figure_1_Manhattan_Mut_Landscape.png', 'Mutation Manhattan Plot'),
            ('Figure_2_Gene_Payload.png', 'Gene Mutation Payload'),
            ('Figure_3_IntraHost_Diversity.png', 'Quasispecies Diversity'),
            ('Figure_4_Lineage_Clustermap.png', 'Missense Lineage Clustermap'),
            ('Figure_5_Ultimate_Quasispecies_OncoPrint.png', 'Mutation OncoPrint'),
        ],
        'maftools': [
            ('mafSummary_TCGA.png', 'MAF Mutation Types'),
            ('Oncoplot.png', 'Mutation OncoPrint'),
            ('TiTv_Summary.png', 'Ti/Tv Summary'),
            ('03c_WholeGenome_Lollipop_Map.png', 'Whole-Genome Lollipop Map'),
            ('04_SomaticInteractions_Network.png', 'Somatic Interactions Network'),
            ('03b_Lollipop_All_Stitched.png', 'Stitched Gene Lollipop Chart'),
        ],
        'snpgenie': [
            ('Fig01_VAF_Spectrum.png', 'VAF Spectrum'),
            ('Fig02_Adjusted_iSNV_Density.png', 'Depth-adjusted iSNV Density'),
            ('Fig03_InterHost_dNdS.png', 'dN vs dS Joint Distribution'),
            ('Fig04_IntraHost_pi.png', 'πN vs πS Joint Distribution'),
            ('Fig05_Gene_dNdS_Stats.png', 'Per-gene dN/dS'),
            ('Fig06_Bootstrapped_dNdS.png', 'Bootstrap Significance'),
            ('Fig07_Mutational_Spectrum.png', 'Mutational Spectrum'),
            ('Fig08_Top_Hotspots.png', 'Top Hotspots'),
            ('Fig09_DualTrack_Window.png', 'Dual-Track Sliding Window'),
            ('Fig10_Trinity_Landscape.png', 'Trinity Landscape'),
            ('Fig11a_PCA_2D.png', '2D PCA Clustering'),
            ('Fig11b_PCA_3D.png', '3D PCA Clustering'),
        ],
        'vcf_merge/figs': [
            ('pca.png', 'Population PCA'),
            ('mantel_ibd.png', 'Isolation-by-Distance Mantel Test'),
            ('distance_clustermap.png', 'Sample Distance Clustermap'),
            ('afs.png', 'Merged Allele Frequency Spectrum'),
            ('dendrogram.png', 'Sample Dendrogram'),
        ],
    }

    for subdir, patterns in chart_patterns.items():
        sd = vdir / subdir
        if not sd.exists(): continue
        for fname, label in patterns:
            fp = sd / fname
            # Fallback: maftools/vcf_viz files may have prefix in filename
            if not fp.exists():
                candidates = list(sd.glob(f"*{fname}"))
                fp = candidates[0] if candidates else fp
            if fp.exists():
                b64 = img_to_base64(fp, max_kb=3500)
                charts[f"{subdir}_{fname}"] = {
                    'label': label,
                    'path': str(fp.relative_to(post_dir.parent)),
                    'base64': b64,
                    'ext': fp.suffix[1:],
                }
    maf_sd = vdir / "maftools"
    if maf_sd.is_dir():
        sk = maf_sd / "SKIPPED_EMPTY_INPUT.txt"
        if sk.exists():
            rsn = ""
            try:
                rsn = sk.read_text(errors='replace').splitlines()[0].strip()
            except Exception:
                pass
            short = rsn.replace("maftools skipped: ", "", 1) if rsn else "empty input / no non-synonymous variants"
            charts["maftools_SKIPPED"] = {"label": "Skipped \u2014 " + short, "path": "", "base64": None, "ext": "txt"}
        for lf in sorted(maf_sd.glob("maftools_03_Lollipop_*.png")):
            gname = lf.stem.replace("maftools_03_Lollipop_", "")
            b64 = img_to_base64(lf, max_kb=3500)
            charts[f"maftools_lollipop_{gname}"] = {
                "label": f"Lollipop {gname}",
                "path": str(lf.relative_to(post_dir.parent)),
                "base64": b64,
                "ext": "png",
            }
    return charts


def collect_dvg_charts(dvg_dir, virus_acc):
    """Collect DVG charts from Stage 9 output."""
    charts = {}
    plots_dir = Path(dvg_dir) / "Summary_Analysis_Report" / "Virus_Specific_Plots"
    if not plots_dir.exists(): return charts
    vdir = plots_dir / virus_acc
    if not vdir.exists():
        short = virus_acc.split(".")[0] if "." in virus_acc else virus_acc
        for d in plots_dir.iterdir():
            if d.is_dir() and (virus_acc in d.name or short in d.name):
                vdir = d; break
        else: return charts
    # 图类白名单 (与生产者 virema_summary_report.R 的 ggsave 输出对齐):
    # 前 4 类为 Circos/Arc/NT 频率/Top 事件, 后 4 类为单病毒散点与
    # 删除/重复覆盖度、donor 断点 NT 偏好 (R 侧 prefix_{name}.png)。
    for suffix, label in [("Circos_4Track.png","DVG Circos"),("Arc_Diagram.png","DVG Arc"),
            ("Acceptor_NT_Freq.png","DVG Acceptor NT"),("Top_Events.png","DVG Top Events"),
            ("Scatterplot.png","DVG Scatterplot"),("Deletions.png","DVG Deletions Coverage"),
            ("Duplications.png","DVG Duplications Coverage"),("Donor_NT_Freq.png","DVG Donor NT")]:
        for f in vdir.glob(f"*{suffix}"):
            b64 = img_to_base64(f, 3000)
            charts[f"dvg_{suffix}"] = {"label":label,"path":str(f.relative_to(dvg_dir)),"base64":b64,"ext":f.suffix[1:]}
    return charts


def collect_virus_data(out_dir, summary_in):
    """Collect per-virus statistics."""
    out = Path(out_dir)
    viruses = {}
    df = safe_read_csv(summary_in)
    if df is None or not HAS_PANDAS: return viruses

    acc_col = next((c for c in ["Rep_Accession", "Accession"] if c in df.columns), None)
    sp_col = next((c for c in ["Adjusted_Species", "Species_NCBI"] if c in df.columns), None)

    for _, row in df.iterrows():
        acc = str(row.get(acc_col, ""))
        if not acc: continue
        viruses[acc] = {
            "species": str(row.get(sp_col, acc)),
            "cpm": _fnum(row.get("Asm_CPM")),
            "coverage": _fnum(row.get("Rep_Coverage(%)")),
            "depth": _fnum(row.get("Rep_MeanDepth")),
            "poisson": _fnum(row.get("Poisson_Ratio")),
            "reads": _fnum(row.get("Asm_EM_Reads")),
        }

    # Count per-virus samples
    for acc in list(viruses.keys()):
        n = (df[acc_col].astype(str) == acc).sum() if acc_col and HAS_PANDAS else "?"
        viruses[acc]["n_samples"] = n

    return viruses


def _collect_variant_summary(out_dir):
    """Aggregate per-virus variant stats from S3 all_summary.tsv."""
    out = Path(out_dir)
    summary_tsv = out / os.path.join(_D['a_align'], 'summary/all_summary.tsv')
    df = safe_read_csv(summary_tsv)
    if df is None: return {}
    acc_col = next((c for c in ["Accession"] if c in df.columns), None)
    if not acc_col: return {}
    # "-" 占位列以字符串 dtype 读入 -> .mean() 抛 TypeError; 先整体数值化 (P0-2 同类)
    if HAS_PANDAS:
        for _c in ("Covered%", "Rep_MeanDepth", "Pi_avr", "Shannon_avr"):
            if _c in df.columns:
                df[_c] = pd.to_numeric(df[_c], errors="coerce")

    def _mean_or_none(grp, col):
        if col not in grp.columns: return None
        m = grp[col].mean()
        return float(m) if pd.notna(m) else None   # 全 "-" -> None

    result = {}
    for acc, grp in df.groupby(acc_col):
        samples = grp["Sample"].nunique() if "Sample" in df.columns else len(grp)
        cov = _mean_or_none(grp, "Covered%")
        depth = _mean_or_none(grp, "Rep_MeanDepth")
        pi = _mean_or_none(grp, "Pi_avr")
        shannon = _mean_or_none(grp, "Shannon_avr")
        try:
            length = _fnum(grp["Length"].iloc[0]) if "Length" in df.columns else None
        except Exception:
            length = None
        entry = {
            "samples": samples,
            "avg_cov": cov, "avg_depth": depth,
            "avg_pi": pi, "avg_shannon": shannon,
            "ref_length": length,
        }
        # [P2-24c] 节段病毒 Accession 为逗号拼接 (batch_virus_variants.py: ",".join(...)),
        # 报告侧按单 acc 查询必然 miss → 注册"拼接键 + 各单段键", 任一命中即取同一组统计
        for key in {str(acc)} | {k.strip() for k in str(acc).split(",") if k.strip()}:
            result[key] = entry
    return result


def _collect_assembly_summary(out_dir):
    """Scan S4 assembly dirs for per-virus per-sample assembly evolution stats.

    Returns dict keyed by accession, each with samples list containing:
      name, denovo_len, denovo_n50, final_len, final_n50, final_n, stages
    where stages is a list of {step, length, n50, n_count} for full evolution tracking.
    """
    out = Path(out_dir)
    asm_dir = _resolve_dir(out_dir, _D['a_assembly'])
    if not asm_dir.is_dir(): return {}
    result = {}
    for vdir in asm_dir.iterdir():
        if not vdir.is_dir() or vdir.name.startswith("run_"): continue
        samples = []
        for sdir in vdir.iterdir():
            if not sdir.is_dir(): continue
            stats_f = sdir / "Global_Evolution_Stats.tsv"
            if not stats_f.is_file():
                samples.append({"name": sdir.name})
                continue
            try:
                sdf = safe_read_csv(stats_f)
                if sdf is None or len(sdf) == 0:
                    samples.append({"name": sdir.name})
                    continue
                # First & last rows
                first = sdf.iloc[0]; last = sdf.iloc[-1]
                denovo_len = float(first.get("Total_Length", 0))
                denovo_n50 = float(first.get("Contig_N50", 0))
                final_len = float(last.get("Total_Length", 0))
                final_n50 = float(last.get("Contig_N50", 0))
                final_n = int(float(last.get("N_Count", 0)))
                # Full stage evolution
                stages = []
                for _, row in sdf.iterrows():
                    stages.append({
                        "step": str(row.get("Step", "")).split(".")[-1] if "." in str(row.get("Step", "")) else str(row.get("Step", "")),
                        "length": float(row.get("Total_Length", 0)),
                        "n50": float(row.get("Contig_N50", 0)),
                        "n_count": int(float(row.get("N_Count", 0))),
                    })
                samples.append({
                    "name": sdir.name,
                    "denovo_len": denovo_len, "denovo_n50": denovo_n50,
                    "final_len": final_len, "final_n50": final_n50, "final_n": final_n,
                    "stages": stages,
                })
            except Exception:
                samples.append({"name": sdir.name})
        # Extract accession from dir name
        import re as _re
        m = _re.search(_ACC_PAT, vdir.name)
        acc = m.group() if m else vdir.name
        result[acc] = {"dir_name": vdir.name, "samples": samples}
    return result


def _collect_detection_summary(out_dir):
    """Parse S1 best-summary TSV to aggregate per-virus detection metrics."""
    out = Path(out_dir)
    tsv = out / os.path.join(_D['a_detect'], 'summary/all_viruses.best.summary.tsv')
    if not tsv.is_file():
        tsv = out / os.path.join(_D['a_detect'], 'summary/all_viruses.summary.tsv')
    df = safe_read_csv(tsv)
    if df is None: return {}
    raw_tsv = out / os.path.join(_D['a_detect'], 'summary/all_viruses.raw.tsv')
    raw_df = safe_read_csv(raw_tsv)
    total_raw = len(raw_df) if raw_df is not None else None
    acc_col = next((c for c in ["Rep_Accession"] if c in df.columns), None)
    if not acc_col: return {}
    # "-" 占位列会以字符串 dtype 读入 -> .mean() 直接 TypeError; 先整体数值化
    if HAS_PANDAS:
        for _c in ("Asm_CPM", "Asm_FPKM", "Rep_Coverage(%)", "Rep_MeanDepth", "Poisson_Ratio"):
            if _c in df.columns:
                df[_c] = pd.to_numeric(df[_c], errors="coerce")
    viruses = {}
    for acc, grp in df.groupby(acc_col):
        sp_col = next((c for c in ["Adjusted_Species", "Species_NCBI"] if c in df.columns), None)
        species = str(grp[sp_col].iloc[0]) if sp_col else str(acc)

        def _gmean(col):
            if col not in grp.columns: return None
            m = grp[col].mean()
            return float(m) if pd.notna(m) else None   # 全 "-" 组 -> None, 不再崩

        viruses[str(acc)] = {
            "species": species,
            "n_samples": len(grp),
            "avg_cpm": _gmean("Asm_CPM"),
            "avg_fpkm": _gmean("Asm_FPKM"),
            "avg_cov": _gmean("Rep_Coverage(%)"),
            "avg_depth": _gmean("Rep_MeanDepth"),
            "avg_poisson": _gmean("Poisson_Ratio"),
        }
    return {"viruses": viruses, "total_best": len(df), "total_raw": total_raw}


def _collect_capheine_summary(out_dir, virus_acc):
    """Parse DRHIP combined_summary.csv for per-gene BUSTED selection results."""
    out = Path(out_dir)
    # capheine 已迁 ⑤ virome_phylo_pipeline (select 模块); 07_capheine 仅作 legacy 产物候选
    vdir = _find_virus_dir(out / "07_capheine", virus_acc)
    if not vdir and PHYLO_ROOT is not None:
        # 新链路: <phylo_root>/**/<acc>/.../select|05_select/capheine/drhip/combined_summary.csv
        pr = Path(PHYLO_ROOT)
        if pr.is_dir():
            for p in sorted(pr.rglob(os.path.join("drhip", "combined_summary.csv"))):
                if virus_acc in p.parts:
                    vdir = p.parent.parent
                    break
    if not vdir: return None
    cs = vdir / "drhip" / "combined_summary.csv"
    df = safe_read_csv(cs)
    if df is None: return None
    genes = []
    for _, row in df.iterrows():
        gene_raw = str(row.get("gene", ""))
        gene_name = gene_raw.split(".part_")[-1] if ".part_" in gene_raw else gene_raw
        genes.append({
            "gene": gene_name,
            "omega3": _fnum(row.get("BUSTED_omega3")),
            "pval": _fnum(row.get("BUSTED_pval"), 1.0),
            "positive": int(_fnum(row.get("positive_sites"))),
            "negative": int(_fnum(row.get("negative_sites"))),
            "total_sites": int(_fnum(row.get("sites"))),
            "n_seq": int(_fnum(row.get("N"))),
        })
    sig_count = sum(1 for g in genes if g["pval"] < 0.05 and g["positive"] > 0)
    return {"genes": genes, "sig_count": sig_count, "total_genes": len(genes), "dir_name": vdir.name}


def _collect_similarity_data(out_dir, virus_acc):
    """Find similarity heatmaps and pairwise stats for a virus."""
    out = Path(out_dir)
    # 布局感知: standard 07_Similarity / legacy 07_similarity / 旧 08_similarity 均兼容
    vdir = _find_virus_dir(locate_analysis_dir(out, "similarity"), virus_acc)
    if not vdir: return {"available": False, "reason": "no similarity directory"}
    # [P2-24d] 原实现硬编码 Mode_Filter → virus_auto_pipeline 以 --mode strict/fill/all 跑出的
    # 结果永远取不到; 改为按优先级查找 Mode_* 候选目录, 存在即用 (最后 glob 兜底任意 Mode_*)
    mat_dir = None
    for _mode in ("Mode_Filter", "Mode_Strict", "Mode_Fill", "Mode_All"):
        _cand = vdir / _mode / "Full_Dataset" / "02_similarity_matrices"
        if _cand.is_dir():
            mat_dir = _cand
            break
    if mat_dir is None:
        for _cand in sorted(vdir.glob("Mode_*/Full_Dataset/02_similarity_matrices")):
            if _cand.is_dir():
                mat_dir = _cand
                break
    if mat_dir is None:
        return {"available": False, "reason": "empty (single sample?)"}
    heatmap_png = mat_dir / "overall_NT_only_heatmap.png"
    dist_png = mat_dir / "overall_NT_only_distribution.png"
    # [P2-24e] 相似性 csv 文件名错配修复: 生产者 virus_auto_pipeline.py 写
    #   {prefix}_{align_method}.csv (execute_single_output/execute_composite_output,
    #   --align_method 默认 sdt_exact) → 旧硬编码 overall_NT_only_pairwise.csv
    #   在默认跑法下必然 miss。改为 glob overall_NT_only*.csv, 命中任一即用;
    #   多命中时优先级 pairwise > sdt_exact > 其他 (mafft/复合矩阵等)。
    _csv_prio = {"overall_NT_only_pairwise.csv": 0, "overall_NT_only_sdt_exact.csv": 1}
    _csv_cands = sorted(mat_dir.glob("overall_NT_only*.csv"))
    csv_f = (min(_csv_cands, key=lambda p: (_csv_prio.get(p.name, 2), p.name))
             if _csv_cands else mat_dir / "overall_NT_only_pairwise.csv")
    if not heatmap_png.exists():
        return {"available": False, "reason": "empty matrices"}
    result = {"available": True, "images": {}, "stats": {}}
    for f, label in [(heatmap_png, "NT Identity Heatmap"), (dist_png, "Identity Distribution")]:
        if f.is_file():
            b64 = img_to_base64(f, max_kb=3000)
            if b64:
                result["images"][label] = b64
    if csv_f.is_file() and HAS_PANDAS:
        try:
            pwd = pd.read_csv(csv_f, index_col=0)
            vals = []
            for i in range(len(pwd.columns)):
                for j in range(i + 1, len(pwd.columns)):
                    v = pwd.iloc[i, j]
                    if not pd.isna(v):
                        vals.append(float(v))
            if vals:
                result["stats"] = {"n_pairs": len(vals), "n_samples": len(pwd.columns),
                    "min_id": min(vals), "max_id": max(vals), "mean_id": sum(vals) / len(vals)}
        except Exception: pass
    return result


def _collect_extract_stats(out_dir):
    """Parse S5 extracted FASTA files for per-virus stats."""
    out = Path(out_dir)
    ext_dir = _resolve_dir(out_dir, _D['a_extract'])
    if not ext_dir.is_dir(): return {}
    result = {}
    for vdir in ext_dir.iterdir():
        if not vdir.is_dir() or vdir.name.startswith("run_"): continue
        lengths = []; n_counts = []; total_bp = 0; total_n = 0
        for fa in list(vdir.glob("*.fasta")) + list(vdir.glob("*.fa")):
            try:
                with open(fa) as fh:
                    seq = ""
                    for line in fh:
                        if line.startswith(">"):
                            if seq:
                                lengths.append(len(seq)); n_counts.append(seq.upper().count("N"))
                                total_bp += len(seq); total_n += seq.upper().count("N")
                            seq = ""
                        else:
                            seq += line.strip()
                    if seq:
                        lengths.append(len(seq)); n_counts.append(seq.upper().count("N"))
                        total_bp += len(seq); total_n += seq.upper().count("N")
            except Exception:
                pass
        if not lengths: continue
        # Extract accession from dir name
        import re as _re
        m = _re.search(_ACC_PAT, vdir.name)
        acc = m.group() if m else vdir.name
        result[acc] = {
            "dir_name": vdir.name,
            "n_contigs": len(lengths),
            "total_bp": total_bp,
            "avg_len": total_bp / len(lengths),
            "min_len": min(lengths),
            "max_len": max(lengths),
            "total_n": total_n,
            "avg_n_pct": (total_n / total_bp * 100) if total_bp > 0 else 0,
        }
    return result


def _find_virus_dir(stage_path, virus_acc):
    """Fuzzy match a virus directory by accession substring."""
    p = Path(stage_path)
    if not p.is_dir(): return None
    if (p / virus_acc).is_dir(): return p / virus_acc
    short = virus_acc.split(".")[0] if "." in virus_acc else virus_acc
    for d in p.iterdir():
        if d.is_dir() and (virus_acc in d.name or short in d.name):
            return d
    return None


def _build_virus_summary_table(viruses):
    """Build a compact per-virus metrics table for the Virus Summary section."""
    rows = ""
    for acc, data in sorted(viruses.items()):
        sp = data.get("species", acc)
        cpm = data.get("cpm", "?")
        cov = data.get("coverage", "?")
        depth = data.get("depth", "?")
        poisson = data.get("poisson", "?")
        reads = data.get("reads", "?")
        n = data.get("n_samples", "?")
        cpm_s = f"{cpm:.1f}" if isinstance(cpm, (int, float)) else str(cpm)
        cov_s = f"{cov:.1f}%" if isinstance(cov, (int, float)) else str(cov)
        depth_s = f"{depth:.1f}x" if isinstance(depth, (int, float)) else str(depth)
        poisson_s = f"{poisson:.2f}" if isinstance(poisson, (int, float)) else str(poisson)
        reads_s = f"{reads:.0f}" if isinstance(reads, (int, float)) else str(reads)
        rows += f'<tr><td><strong>{_esc(sp)}</strong><br><span style="font-size:10px;color:#888">{_esc(acc)}</span></td><td>{n}</td><td>{cov_s}</td><td>{depth_s}</td><td>{cpm_s}</td><td>{poisson_s}</td><td>{reads_s}</td></tr>'
    tbl = f'<div class="chart-card"><table style="width:100%;font-size:12px;border-collapse:collapse"><thead><tr><th>Species</th><th>Samples</th><th>Coverage</th><th>Depth</th><th>CPM</th><th>Poisson</th><th>Reads</th></tr></thead><tbody>{rows}</tbody></table></div>'
    return tbl


def _build_stage_data_json(viruses, variant_summary, assembly_summary, detection_summary, extract_stats, overview_metrics):
    """Build a JSON string with detailed per-stage data for browser-based AI."""
    import json as _json
    ds = detection_summary or {}
    # Build comprehensive per-virus data
    virus_detail = []
    for acc in sorted(viruses.keys()):
        sp = viruses[acc].get("species", acc)
        n = viruses[acc].get("n_samples", "?")
        cov = viruses[acc].get("coverage", "?")
        depth = viruses[acc].get("depth", "?")
        cpm = viruses[acc].get("cpm", "?")
        vs = variant_summary.get(acc, {})
        asm = assembly_summary.get(acc, {})
        es = extract_stats.get(acc, {})
        detail = f"{sp} ({acc}): {n} samples"
        if cov: detail += f", coverage={cov}%"
        if depth: detail += f", depth={depth}x"
        if cpm: detail += f", CPM={cpm}"
        if vs.get("avg_pi") is not None: detail += f", nucleotide diversity pi={vs['avg_pi']:.4f}"
        valid_asm = [s for s in asm.get("samples", []) if s.get("final_len")]
        if valid_asm:
            best_len = max(s.get("final_len", 0) or 0 for s in valid_asm)
            n_with_n = sum(1 for s in valid_asm if s.get("final_n", 0) > 0)
            detail += f", {len(valid_asm)} assemblies (best={best_len:.0f}bp"
            if n_with_n > 0: detail += f", {n_with_n} with gaps"
            detail += ")"
        if es.get("n_contigs"): detail += f", {es['n_contigs']} extracted contigs (avg {es['avg_len']:.0f}bp, N={es['avg_n_pct']:.1f}%)"
        virus_detail.append(detail)

    # Rich per-virus detail for AI
    virus_rich = []
    for acc in sorted(viruses.keys()):
        sp = viruses[acc].get("species", acc)
        n = viruses[acc].get("n_samples", "?")
        cov = viruses[acc].get("coverage", "?")
        depth = viruses[acc].get("depth", "?")
        cpm = viruses[acc].get("cpm", "?")
        vs = variant_summary.get(acc, {})
        asm = assembly_summary.get(acc, {})
        es = extract_stats.get(acc, {})
        entry = {"species": sp, "accession": acc, "samples": n, "coverage": cov, "depth": depth, "cpm": cpm}
        if vs.get("avg_pi") is not None: entry["pi"] = round(vs["avg_pi"], 4)
        if vs.get("avg_shannon") is not None: entry["shannon"] = round(vs["avg_shannon"], 4)
        valid = [s for s in asm.get("samples", []) if s.get("final_len")]
        if valid:
            best = max(s.get("final_len", 0) or 0 for s in valid)
            n_gap = sum(1 for s in valid if s.get("final_n", 0) > 0)
            entry["assembly"] = f"{len(valid)} assemblies, best={best:.0f}bp"
            if n_gap > 0: entry["assembly"] += f", {n_gap} with gaps(N)"
        if es.get("n_contigs"): entry["extract"] = f"{es['n_contigs']} contigs, avg={es['avg_len']:.0f}bp, N%={es['avg_n_pct']:.1f}%"
        virus_rich.append(entry)
    def _json_safe(o):
        import numpy as _np
        if isinstance(o, (_np.integer,)): return int(o)
        if isinstance(o, (_np.floating,)): return float(o)
        if isinstance(o, (_np.ndarray,)): return o.tolist()
        raise TypeError
    virus_detail_str = _json.dumps(virus_rich, ensure_ascii=False, default=_json_safe)

    data = {
        "samples": overview_metrics.get("n_samples","?"),
        "viruses": len(viruses),
        "_virus_detail": virus_rich,
        "_virus_detail_str": virus_detail_str,
        "s1": {"desc": "Salmon pseudo-alignment + Poisson filtering", "raw": overview_metrics.get("raw_detections","?"), "high_conf": overview_metrics.get("best_detections","?"), "virus_detail": virus_rich},
        "s2": {"desc": "Multi-dimensional filter: cov>=50%, depth>=5x, reads>=100", "pre": overview_metrics.get("raw_detections","?"), "post": overview_metrics.get("filtered_records","?")},
        "s3": {"desc": "FreeBayes + SnpEff + SNPGenie", "per_virus": {acc: {"samples": variant_summary.get(acc,{}).get("samples","?"), "cov": f"{variant_summary.get(acc,{}).get('avg_cov','?'):.1f}%" if variant_summary.get(acc,{}).get('avg_cov') is not None else "?", "depth": f"{variant_summary.get(acc,{}).get('avg_depth','?'):.1f}x" if variant_summary.get(acc,{}).get('avg_depth') is not None else "?", "pi": f"{variant_summary.get(acc,{}).get('avg_pi','?'):.4f}" if variant_summary.get(acc,{}).get('avg_pi') is not None else "?"} for acc in sorted(viruses.keys())}},
        "s4": {"desc": "12-step de novo assembly", "per_virus": {}},
        "s5": {"desc": "Contig extraction + N-fill", "per_virus": {}},
        "s7": {"desc": "HyPhy selection: FEL/MEME/BUSTED/PRIME", "coding_viruses": sum(1 for acc in viruses if variant_summary.get(acc, {}).get("ref_length", 0) > 1000)},
        "s9": {"desc": "ViReMa DVG detection", "viruses_with_dvg": overview_metrics.get("dvg_viruses","?")},
    }
    for acc in sorted(viruses.keys()):
        asm = assembly_summary.get(acc, {})
        valid = [s for s in asm.get("samples", []) if s.get("final_len")]
        if valid:
            best = max(s.get("final_len", 0) or 0 for s in valid)
            data["s4"]["per_virus"][acc] = f"{len(valid)} assemblies, best={best:.0f}bp"
        es = extract_stats.get(acc, {})
        if es:
            data["s5"]["per_virus"][acc] = f"{es['n_contigs']} contigs, avg={es['avg_len']:.0f}bp, N%={es['avg_n_pct']:.1f}%"
    return _json.dumps(data, ensure_ascii=False, default=_json_safe)


# ── Backward-compatible directory resolver (module-level, shared by all functions) ──
_DIR_FALLBACK = {
    _D['a_align']: "03_Variants",
    _D['a_variants']: "03_Variants",
    _D['a_post']: "04_PostAnalysis",
    _D['a_assembly']: "04_assembly",
    _D['a_extract']: "05_extraction",
}

def _resolve_dir(base_dir, path_segment):
    """Resolve directory path with new→old naming fallback."""
    d = Path(base_dir) / path_segment
    if d.is_dir():
        return d
    if path_segment in _DIR_FALLBACK:
        fb = Path(base_dir) / _DIR_FALLBACK[path_segment]
        if fb.is_dir():
            return fb
    return d


# ── Pipeline methods section (auto-generated from stage metadata) ──
# [P2-24f] 编号以 auto_known_virus.py 的 STAGES 注册表为准 (本管线 v3.0 起为 9 阶段):
#   1 detect / 2 filter / 3 variants / 4 post / 5 full / 6 extract / 7 similarity /
#   8 dvg / 9 report。Selection (capheine) 已迁 ⑤, 作遗留附录用字母后缀 "6b" 占位,
#   不占用 1-9 数字槽, 保证编号唯一递增 (原先 DVG 与 Report 都标 9)。
STAGE_METHODS = [
    (1, "Detection", "batch_virus_depth.py",
     "Virus detection was performed using Salmon pseudo-alignment in quantification mode "
     "against a curated plant virus reference database. Read counts were normalized to "
     "CPM and FPKM. A Poisson Ratio filter was applied to distinguish uniformly covered "
     "true viral reads from localized spurious alignments (threshold ≥ 0.3). Parameters: "
     "min coverage 10%, min mean depth 0.5×, min TPM 1.0, species ANI threshold 95%."),

    (2, "Filtering", "utils/filter_summary.py",
     "High-confidence detections were retained by requiring Rep_Coverage ≥ 50%, "
     "Rep_MeanDepth ≥ 5×, and ≥ 100 estimated mapped reads. Optional keyword-based "
     "taxonomic filtering available for targeted virus groups."),

    (3, "Variant Calling", "batch_virus_variants.py",
     "Variants were called using iVar in haploid mode with minimum quality 20, "
     "minimum depth 5×, and minimum allele frequency 5%. Functional annotation "
     "was performed with SnpEff using custom databases built from viral GenBank records. "
     "Population genetic parameters (π, dN/dS) were computed using SNPGenie."),

    (4, "Post-hoc Visualization", "virus_vcf_pipeline.py + virus_variants_analyzer.py",
     "Per-virus VCF files were merged with bcftools and subjected to sample-level QC "
     "filtering. Jaccard and Hamming distance matrices were computed from SNP genotypes, "
     "and NJ trees constructed via UPGMA. PCA was performed on kinship matrices. "
     "Nine publication-grade figures were generated per virus: genomic mutation landscape, "
     "Ts/Tv ratio, functional impact, hierarchical clustering, allele frequency spectrum, "
     "variant density with gene annotation, AF by mutation type, sliding-window π + Tajima's D, "
     "and PCA sub-lineage analysis with auto-optimized K-means clustering."),

    (5, "Full-length Assembly", "virus-full.py (OmniVirusAssembler V9)",
     "De novo genome assembly employed a 12-step pipeline: (1-2) MEGAHIT/SPAdes/PenguIN "
     "multi-tool assembly with refineC merging; (3) Shiver-like orientation correction; "
     "(4-8) Divine Fusion reference-guided gap resolution, PVGA read extension with "
     ">100bp gap splitting, rmDup deduplication at 95% identity; (9) 3-round minimap2 + "
     "viral_consensus iterative polishing; (10) Dual-engine gap filling via gmcloser + "
     "abyss-sealer; (11) Circularity detection via terminal BLASTN self-alignment; "
     "(12) Full-lifecycle coverage visualization. Fast-track mode bypassed intermediate "
     "fusion steps for high-quality initial skeletons. Assembly quality assessed at "
     "≥98% reference coverage and ≥95% N50 for perfect designation."),

    (6, "Contig Extraction", "utils/extract_full_fasta.py",
     "The longest assembled contig per virus was extracted from the 12-step pipeline output. "
     "N-filled regions were resolved using reference-guided gap filling. Contigs with "
     ">10% ambiguous bases were excluded."),

    ("6b", "Selection Analysis", "virome_phylo_pipeline (⑤ output; v3.0 migration)",
     "Codon-level positive selection is produced by ⑤ virome_phylo_pipeline (select module) since "
     "the v3.0 migration; the legacy capheine_pipeline.py / 07_capheine directory is retained only "
     "as an artifact candidate. The analysis integrates "
     "CAWLIGN codon-aware alignment, HyPhy CLN duplicate removal, IQ-TREE maximum likelihood "
     "phylogeny (GTR+I+G), and six HyPhy methods: FEL (Fixed Effects Likelihood), "
     "MEME (Mixed Effects Model of Evolution), PRIME (Property-Informed Models of Evolution), "
     "BUSTED (Branch-Site Unrestricted Statistical Test), CONTRAST-FEL, and RELAX. "
     "Results were aggregated via DRHIP. Non-coding viruses (viroids) were automatically excluded."),

    (7, "Similarity Panorama", "virus_auto_pipeline.py",
     "Pairwise nucleotide and amino acid similarity matrices were computed for each virus "
     "across all samples using Biopython pairwise alignment (SDT-strict mode, no gap penalty). "
     "CD-HIT soft deduplication at adaptive identity thresholds removed redundant sequences. "
     "SciPy hierarchical clustering produced optimally ordered triangular heatmaps at 600 DPI. "
     "Multi-gene concatenated CDS matrices were also generated when multiple CDS regions were available."),

    (8, "DVG & Recombination", "batch_virema_dvg.py",
     "Defective viral genomes and recombination events were detected using ViReMa "
     "with seed length 25, micro-indel threshold 15 bp. Recombination landscapes "
     "were visualized via Circos plots showing junction networks. Viruses with "
     "≥80% genome coverage were included in DVG analysis."),

    (9, "Report Generation", "generate_pipeline_report.py",
     "This interactive HTML report was auto-generated from pipeline outputs, embedding "
     "all figures, tables, and result files. Includes AI-assisted interpretation via "
     "LLM integration (OpenAI/DeepSeek/Ollama compatible)."),
]


def _build_methods_section(viruses, overview_metrics, total_records):
    """Generate dynamic methods section from STAGE_METHODS metadata."""
    n_viruses = len(viruses)

    # Pipeline summary table
    rows = ""
    for num, name, script, _ in STAGE_METHODS:
        rows += f"<tr><td>{num}. {name}</td><td>{_esc(script)}</td></tr>"

    # Stage details (collapsible)
    details = ""
    for num, name, script, desc in STAGE_METHODS:
        details += f'''<details><summary><strong>{name} (Stage {num})</strong> — <code>{_esc(script)}</code></summary>
        <div class="suggestion" style="padding:10px 16px;font-size:13px;line-height:1.7">{desc}</div></details>\n'''

    return f'''
    <h3 id="paper-reference">Pipeline Methods</h3>
    <div class="paper-block">
        <p>The Known Virus Analysis Pipeline processed <strong>{n_viruses} virus species</strong> across
        {overview_metrics.get("n_samples","?")} samples through 9 automated stages.</p>

        <h4>Stage Overview</h4>
        <table>
            <tr><th>Stage</th><th>Primary Script</th></tr>
            {rows}
        </table>

        <h4>Detailed Methods</h4>
        {details}

        <h4>Key Quantitative Results</h4>
        <ul>
            <li><strong>Detection</strong>: {total_records} virus-sample records identified across {n_viruses} species</li>
            <li><strong>Assembly</strong>: {overview_metrics.get("assemblies","?")} de novo assemblies generated</li>
            <li><strong>Variant Calling</strong>: SNVs called at minimum 5× depth with 5% allele frequency (iVar haploid mode)</li>
            <li><strong>Selection</strong>: HyPhy codon-level testing (FEL/MEME/BUSTED/PRIME) for coding viruses</li>
        </ul>
    </div>'''


def _build_funnel_section(out_dir):
    """五阶段漏斗 (原始检出→best→高置信→全长组装→提取留存): 调 stage_funnel_summary 生成表+图并渲染 HTML。

    任何失败都返回 '' (报告其余部分不受影响); 产物落在 09_Report/stage_funnel/。
    """
    out = Path(out_dir)
    try:
        _here = str(Path(__file__).parent)
        if _here not in sys.path:
            sys.path.insert(0, _here)
        import stage_funnel_summary as _sfs
        res = _sfs.build(out, out / os.path.join(_D['a_report'], 'stage_funnel'), make_plots=True)
    except Exception as e:
        print(f"  (stage funnel section skipped: {e})")
        return ""
    stages = res.get("stages")
    if stages is None or len(stages) == 0 or int(stages["Records"].fillna(0).sum()) == 0:
        return ""

    # 指标卡片: 每阶段记录数 + 环比留存
    cards = ""
    for _, r in stages.iterrows():
        ret = r.get("Retention_vs_Previous(%)", "")
        ret_s = f' <span style="color:#c0392b;font-size:11px">&#8600;{ret}%</span>' if ret != "" else ""
        cards += (f'<div class="card"><div class="value">{r["Records"]}</div>'
                  f'<div class="label">{_esc(str(r["Stage"]))}{ret_s}<br>'
                  f'<span style="font-size:10px">{r["Species"]} species &middot; {r["Samples"]} samples</span></div></div>')

    # 图: 漏斗图 + 物种留存热图
    charts = ""
    for png, title in [(res.get("funnel_png"), "stage_funnel.png &mdash; records &amp; species per stage"),
                       (res.get("heatmap_png"), "virus_stage_matrix.png &mdash; per-virus retention")]:
        if png and Path(png).is_file():
            b64 = img_to_base64(str(png), 2000)
            if b64:
                charts += (f'<div class="chart-card"><div class="chart-title">{title}</div>'
                           f'<img src="data:image/png;base64,{b64}" loading="lazy" alt="{title}"></div>')
            else:
                charts += _img_placeholder(title)
    chart_html = f'<div class="chart-gallery">{charts}</div>' if charts else ""

    # 漏斗摘要表 (5 行, 静态)
    th = ''.join(f'<th>{_esc(c)}</th>' for c in stages.columns)
    trs = ""
    for _, r in stages.iterrows():
        trs += '<tr>' + ''.join(f'<td>{_esc(str(v))}</td>' for v in r.values) + '</tr>'
    summary_tbl = (f'<div class="tb-scroll"><table style="font-size:12px">'
                   f'<thead><tr>{th}</tr></thead><tbody>{trs}</tbody></table></div>')

    # 物种 × 阶段矩阵表 (分页, 复用 pg_tbl)
    matrix_tbl = ""
    mpath = Path(res["outdir"]) / "virus_stage_matrix.tsv"
    if mpath.is_file():
        try:
            rows = [ln.rstrip("\n").split("\t") for ln in
                    mpath.read_text(encoding="utf-8", errors="replace").splitlines() if ln.strip()]
            if len(rows) > 1:
                tid = 'tbl_funnel_matrix'
                ar = rows[1:]
                pp = 10
                tp = max((len(ar) + pp - 1) // pp, 1)
                mh = ''.join(f'<th>{_esc(c)}</th>' for c in rows[0])
                matrix_tbl = f'<h3 style="color:var(--ink-secondary);font-size:14px;margin-top:12px">Per-Virus Stage Matrix ({len(ar)} species)</h3>'
                matrix_tbl += f'<div class="pg-nav" id="{tid}_nav" style="display:flex;gap:8px;align-items:center;padding:6px 10px;background:#f8f9fa;border-bottom:1px solid #e0e0e0;font-size:12px">'
                matrix_tbl += f'<button onclick="pg_tbl(\'{tid}\',-{pp})" style="padding:3px 10px;cursor:pointer;border:1px solid #ccc;border-radius:3px;background:#fff">&laquo; Prev</button>'
                matrix_tbl += f'<span id="{tid}_info" style="color:#666">Page 1 of {tp}</span>'
                matrix_tbl += f'<button onclick="pg_tbl(\'{tid}\',{pp})" style="padding:3px 10px;cursor:pointer;border:1px solid #ccc;border-radius:3px;background:#fff">Next &raquo;</button>'
                matrix_tbl += f'</div><div class="tb-scroll"><table id="{tid}"><thead><tr>{mh}</tr></thead><tbody>'
                for ri, r in enumerate(ar):
                    cls = ' style="display:none"' if ri >= pp else ''
                    matrix_tbl += f'<tr{cls}><td>' + '</td><td>'.join(_esc(c) for c in r) + '</td></tr>'
                matrix_tbl += '</tbody></table></div>'
        except Exception:
            pass

    outdir_note = (f'<p style="font-size:11px;color:var(--ink-secondary)">Data &amp; figures: '
                   f'<code>{_esc(str(res["outdir"]))}</code> '
                   f'(stage_funnel_summary.tsv / virus_stage_matrix.tsv / *.png *.pdf)</p>')

    return (f'<section id="funnel"><h2>Stage Funnel: Detection &rarr; Full-Length</h2>'
            f'<p style="color:var(--ink-secondary);font-size:13px;line-height:1.6">'
            f'Five-stage retention from raw detection to usable full-length genomes: '
            f'raw detection (pre-Poisson) &rarr; best (Poisson-filtered, per-species representative) &rarr; '
            f'high-confidence filter &rarr; de novo full-length assembly &rarr; extracted full-length (length/N filtered). '
            f'Stage names: raw detections &rarr; best detections &rarr; filtered detections &rarr; '
            f'full-length assemblies &rarr; filtered full-length (extraction retained).</p>'
            f'<div class="card-row" style="margin-top:8px">{cards}</div>'
            f'<h3 style="color:var(--ink-secondary);font-size:14px;margin-top:12px">Funnel Summary</h3>'
            f'{summary_tbl}{chart_html}{matrix_tbl}{outdir_note}</section>')


def generate_html(out_dir, out_html, viruses, ai_html="", stage_summaries=None):
    """Generate interactive HTML report with sidebar navigation."""
    now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    out = Path(out_dir)

    def _d(seg):
        return _resolve_dir(out, seg)

    post_dir = _resolve_dir(out, _D['a_post'])

    total_records = sum(_fnum(v.get("n_samples")) for v in viruses.values()) if isinstance(
        next(iter(viruses.values()), {}).get("n_samples", 0), _INT_TYPES) else "?"

    # Pre-compute global overview metrics
    overview_metrics = {}
    # Detection counts
    raw_tsv_ov = out / os.path.join(_D['a_detect'], 'summary/all_viruses.raw.tsv')
    raw_ov = safe_read_csv(raw_tsv_ov)
    best_tsv_ov = out / os.path.join(_D['a_detect'], 'summary/all_viruses.best.summary.tsv')
    best_ov = safe_read_csv(best_tsv_ov)
    overview_metrics["raw_detections"] = len(raw_ov) if raw_ov is not None else "?"
    overview_metrics["best_detections"] = len(best_ov) if best_ov is not None else "?"
    overview_metrics["n_samples"] = best_ov["Sample"].nunique() if best_ov is not None and "Sample" in best_ov.columns else "?"
    # Filter counts
    hc_tsv_ov = out / os.path.join(_D['a_filter'], 'high_conf.summary.tsv')
    hc_ov = safe_read_csv(hc_tsv_ov)
    overview_metrics["filtered_records"] = len(hc_ov) if hc_ov is not None else "?"
    # v21: metadata association coverage (编排器现写 04_Variants/metadata_association;
    #      02_Filtering 为 legacy 旧位置回退)
    meta_tsv_ov = post_dir / 'metadata_association/Viral_Infection_with_Metadata_Summary.tsv'
    if not meta_tsv_ov.is_file():
        meta_tsv_ov = out / os.path.join(_D['a_filter'], 'metadata_association/Viral_Infection_with_Metadata_Summary.tsv')
    meta_ov = safe_read_csv(meta_tsv_ov)
    overview_metrics["meta_records"] = len(meta_ov) if meta_ov is not None else 0
    # Assembly counts
    asm_dir_ov = _d(_D['a_assembly'])
    n_asm = 0
    if asm_dir_ov.is_dir():
        for vd in asm_dir_ov.iterdir():
            if vd.is_dir() and not vd.name.startswith("run_"):
                n_asm += sum(1 for d in vd.iterdir() if d.is_dir())
    overview_metrics["assemblies"] = n_asm
    # DVG counts
    dvg_dir_ov = out / os.path.join(_D['a_dvg'], 'Summary_Analysis_Report/Virus_Specific_Plots')
    overview_metrics["dvg_viruses"] = sum(1 for d in dvg_dir_ov.iterdir() if d.is_dir()) if dvg_dir_ov.is_dir() else 0

    # Pre-load stage summaries
    detection_summary = _collect_detection_summary(out_dir)
    variant_summary = _collect_variant_summary(out_dir)
    assembly_summary = _collect_assembly_summary(out_dir)
    extract_stats = _collect_extract_stats(out_dir)

    # Helper: collect stage-level files for embedding
    def _stage_files(outd, sn):
        """Return HTML for stage-level result files (plots/tables from stage directory)."""
        dirs = {1:_D['a_detect'], 2:_D['a_filter'],
                4:_D['a_assembly'],
                6:_D['a_post'], 7:"07_capheine"}
        patterns = {1:["summary/all_viruses.best.summary.tsv","plots/virus_analysis/freq_multi_metrics_log10.png"],
                    2:["**/metadata_association/**/*"],
                    # 4=组装: 每任务覆盖度双轨图在 12.Coverage_Visualization/ 下;
                    # 该目录可被 clean_intermediates.py 清理 → 由下方 fallback 补
                    # 06_Extraction 的 assembly_stats, 保证 Stage 5 区块不空
                    4:["**/12.Coverage_Visualization/*_Coverage_*", "**/*final_coverage*"],
                    6:[],
                    7:["**/selection_per_gene.*"]}
        if sn not in dirs: return "", ""
        spath = outd / dirs[sn]
        pats = patterns.get(sn, ["*.png"])
        if sn == 2:
            # metadata_association 写于 04_Variants (编排器); 02_Filtering 为 legacy 回退。
            # spath 直接指向 metadata_association 目录时, 模式改为全量收集
            ma_post = outd / _D['a_post'] / "metadata_association"
            if ma_post.is_dir():
                spath = ma_post
                pats = ["*", "**/*"]
        if not spath.exists(): return "", ""
        html_tables = ""; html_charts = ""
        collected = []
        if sn == 4:
            # 组装块兜底: 12.Coverage_Visualization 被清理后, 用 06_Extraction 的
            # 统计表/图补位 (assembly_stats.* 由 full 尾段恒定产出)
            _ext = _resolve_dir(outd, _D['a_extract'])
            for _n in ("assembly_stats.png", "assembly_stats.pdf", "assembly_stats.csv"):
                _f = _ext / _n
                if _f.is_file():
                    collected.append(_f)
        for pat in pats:
            for f in spath.rglob(pat) if "**" in pat else spath.glob(pat):
                if f.is_file(): collected.append(f)
        if sn == 7 and PHYLO_ROOT is not None:
            # --phylo-root: 追加 ⑤ 管线 capheine 位点统计图 (selection_per_gene.png/pdf)
            pr = Path(PHYLO_ROOT)
            if pr.is_dir():
                for f in sorted(pr.rglob("selection_per_gene.*")):
                    if f.is_file(): collected.append(f)
        # Dedup: skip PDF when PNG with same stem exists
        png_stems = {f.stem for f in collected if f.suffix.lower() in ('.png','.jpg','.jpeg','.svg')}
        for f in collected:
            if f.suffix.lower() == '.pdf' and f.stem in png_stems:
                continue
            try:
                fn = str(f.relative_to(spath)) if f.parent != spath else f.name
            except ValueError:
                # ⑤ phylo 产物不在 spath 下: 以 <病毒acc>/文件名 作标签
                try:
                    fn = f"{f.parts[f.parts.index('gcva_analysis') - 1]}/{f.name}"
                except ValueError:
                    fn = f.name
            if fn.endswith('.tsv') or fn.endswith('.csv'):
                try:
                    rows = []
                    with open(str(f),'r',encoding='utf-8',errors='replace') as fh:
                        for line in fh:
                            rows.append([c.strip() for c in line.rstrip().split('\t')])
                    if rows and len(rows) > 1:
                        tid = 'tbl_' + str(abs(hash(str(f))))[:8]
                        ar = rows[1:]; tr = len(ar); pp = 10; tp = (tr+pp-1)//pp
                        th = ''.join(f'<th>{_esc(c)}</th>' for c in rows[0])
                        html_tables += f'<div class="chart-card"><div class="chart-title">{_esc(fn)} ({tr} rows, {tp} pages)</div>'
                        html_tables += f'<div class="pg-nav" id="{tid}_nav" style="display:flex;gap:8px;align-items:center;padding:6px 10px;background:#f8f9fa;border-bottom:1px solid #e0e0e0;font-size:12px">'
                        html_tables += f'<button onclick="pg_tbl(\'{tid}\',-{pp})" style="padding:3px 10px;cursor:pointer;border:1px solid #ccc;border-radius:3px;background:#fff">&laquo; Prev</button>'
                        html_tables += f'<span id="{tid}_info" style="color:#666">Page 1 of {tp}</span>'
                        html_tables += f'<button onclick="pg_tbl(\'{tid}\',{pp})" style="padding:3px 10px;cursor:pointer;border:1px solid #ccc;border-radius:3px;background:#fff">Next &raquo;</button>'
                        html_tables += f'</div><div class="tb-scroll"><table id="{tid}"><thead><tr>{th}</tr></thead><tbody>'
                        for ri, r in enumerate(ar):
                            cls = ' style="display:none"' if ri >= pp else ''
                            html_tables += f'<tr{cls}><td>' + '</td><td>'.join(_esc(c) for c in r) + '</td></tr>'
                        html_tables += '</tbody></table></div></div>'
                except Exception: pass
            elif fn.endswith('.pdf'):
                b64 = img_to_base64(str(f), 1000)
                if b64:
                    html_charts += f'<div class="chart-card"><div class="chart-title">{_esc(fn)}</div><object data="data:application/pdf;base64,{b64}" type="application/pdf" width="100%" height="500px"><p>PDF</p></object></div>'
                else:
                    html_charts += _img_placeholder(fn)
            elif fn.endswith(('.png','.jpg','.jpeg','.svg')):
                b64 = img_to_base64(str(f), 2000)
                if b64:
                    mime = "image/png"
                    if fn.endswith('.jpg') or fn.endswith('.jpeg'): mime = "image/jpeg"
                    elif fn.endswith('.svg'): mime = "image/svg+xml"
                    html_charts += f'<div class="chart-card"><div class="chart-title">{_esc(fn)}</div><img src="data:{mime};base64,{b64}" loading="lazy" alt="{_esc(fn)}"></div>'
                else:
                    html_charts += _img_placeholder(fn)
        return html_tables, html_charts

    # Build sidebar nav items + stage-per-virus body sections
    stage_cfg = [
        ("S1","Detection","Salmon/Kallisto pseudo-alignment + Poisson ratio false-positive filtering.",False),
        ("S2","Filter","Multi-dimensional filter: coverage, depth, reads, TPM, Poisson ratio, ANI, keyword.",False),
        ("S3","Variants","FreeBayes/iVar/LoFreq variant calling. Dynamic VCF filtering. SnpEff + SNPGenie.",True),
        ("S4","Assembly","12-step de novo assembly: MEGAHIT/SPAdes to PVGA to polishing to gap-filling.",True),
        ("S5","Extract","Longest contig extraction with N-fill via reference global pairwise alignment.",True),
        ("S6","Post-hoc","Post-hoc per-virus: mutation landscape, OncoPrint, dN/dS, PCA.",True),
        ("S7","Capheine","HyPhy positive selection: FEL/MEME/PRIME/BUSTED/CONTRASTFEL/RELAX.",True),
        ("S8","Similarity","Pairwise genome similarity: SDT-style NT/AA identity matrices + clustering.",True),
        ("S9","DVG","ViReMa DVG detection: Circos 4-track recombination plots + arc diagrams.",True),
    ]

    nav_items = '<li class="nav-item"><a href="#overview" class="nav-link active">Global Overview</a></li>'
    nav_items += '<li class="nav-item"><a href="#funnel" class="nav-link">Stage Funnel</a></li>'
    nav_items += '<li class="nav-header" style="padding:8px 20px 4px;font-size:12px;color:var(--ink-secondary);font-weight:600;border-top:1px solid var(--border);margin-top:4px">Pipeline Stages</li>'
    nav_items += '<div class="nav-sub">'

    stage_sections = ""
    for i, (sid, sname, sdesc, has_sub) in enumerate(stage_cfg):
        sn = i + 1
        nav_items += f'<li class="nav-item"><a href="#stage-{sn}" class="nav-link stage-toggle" style="font-size:12px;padding-left:24px;font-weight:600">{sid}: {sname}</a></li>'
        if has_sub:
            nav_items += '<div class="nav-sub stage-sub">'
        # Stage overview section
        stage_sections += f'<section id="stage-{sn}"><h2>{sid}: {sname} <button onclick="runStageAI({sn})" class="ai-stage-btn" title="AI summarize this stage">AI</button></h2><p style="color:var(--ink-secondary);font-size:13px;line-height:1.6">{sdesc}</p>'
        if stage_summaries and sn in stage_summaries:
            stage_sections += f'<details class="ai-stage-detail"><summary>AI Summary</summary><p style="font-size:12px;color:var(--ink);line-height:1.6;padding:6px 10px;background:var(--accent-subtle);border-radius:4px;margin:4px 0">{_esc(stage_summaries[sn])}</p></details>'
        # Stage 1 & 2: cards go RIGHT HERE, before _stage_files
        if sn == 1 and detection_summary:
            ds = detection_summary
            total_raw = ds.get("total_raw", "?")
            total_best = ds.get("total_best", "?")
            n_virus_raw = "?"
            # Load raw TSV for raw-level stats
            raw_tsv_s1 = out / os.path.join(_D['a_detect'], 'summary/all_viruses.raw.tsv')
            raw_s1 = safe_read_csv(raw_tsv_s1)
            if raw_s1 is not None:
                n_samples = raw_s1["Sample"].nunique() if "Sample" in raw_s1.columns else "?"
                acc_raw = next((c for c in ["Accession","Rep_Accession"] if c in raw_s1.columns), None)
                n_virus_raw = raw_s1[acc_raw].nunique() if acc_raw else "?"
            else:
                n_samples = "?"
            # Best summary stats
            best_tsv_s1 = out / os.path.join(_D['a_detect'], 'summary/all_viruses.best.summary.tsv')
            best_s1 = safe_read_csv(best_tsv_s1)
            if best_s1 is not None:
                acc_best = next((c for c in ["Rep_Accession","Accession"] if c in best_s1.columns), None)
                n_virus_best = best_s1[acc_best].nunique() if acc_best else "?"
                n_records_best = len(best_s1)
            else:
                n_virus_best = "?"; n_records_best = "?"
            # Poisson filtered = raw - best
            poisson_filtered = total_raw - total_best if isinstance(total_raw,int) and isinstance(total_best,int) else "?"
            rate = f"{total_best}/{total_raw} ({total_best/total_raw*100:.0f}%)" if isinstance(total_raw,int) and isinstance(total_best,int) and total_raw>0 else "?"
            stage_sections += f'<div class="card-row" style="margin-top:8px"><div class="card"><div class="value">{total_raw or "?"}</div><div class="label">Raw Detections (raw.tsv)</div></div><div class="card"><div class="value">{poisson_filtered}</div><div class="label">Poisson Filtered</div></div><div class="card"><div class="value">{total_best}</div><div class="label">High-Conf (best.tsv)</div></div><div class="card"><div class="value">{n_samples}</div><div class="label">Samples</div></div><div class="card"><div class="value">{n_virus_best}</div><div class="label">Virus Species</div></div><div class="card"><div class="value">{rate}</div><div class="label">Retention</div></div></div>'
        if sn == 2:
            # Quick-load data for cards (must load before using)
            raw_tsv_s2 = out / os.path.join(_D['a_detect'], 'summary/all_viruses.raw.tsv')
            hc_tsv_s2 = out / os.path.join(_D['a_filter'], 'high_conf.summary.tsv')
            ps_tsv_s2 = out / os.path.join(_D['a_filter'], 'filter_stats.per_sample.tsv')
            raw_s2 = safe_read_csv(raw_tsv_s2); hc_s2 = safe_read_csv(hc_tsv_s2); ps_s2 = safe_read_csv(ps_tsv_s2)
            raw_n = len(raw_s2) if raw_s2 is not None else "?"
            hc_n = len(hc_s2) if hc_s2 is not None else "?"
            rate2 = f"{hc_n}/{raw_n} ({hc_n/raw_n*100:.0f}%)" if isinstance(raw_n,int) and isinstance(hc_n,int) and raw_n>0 else "?"
            n_samples2 = len(ps_s2) if ps_s2 is not None else "?"
            total_disc = raw_n - hc_n if isinstance(raw_n,int) and isinstance(hc_n,int) else "?"
            stage_sections += f'<div class="card-row" style="margin-top:8px"><div class="card"><div class="value">{raw_n}</div><div class="label">Raw Records</div></div><div class="card"><div class="value">{hc_n}</div><div class="label">High-Conf</div></div><div class="card"><div class="value">{total_disc}</div><div class="label">Filtered Out</div></div><div class="card"><div class="value">{n_samples2}</div><div class="label">Samples</div></div><div class="card"><div class="value">{rate2}</div><div class="label">Retention</div></div></div>'
        # Add stage-level result files
        sfiles_tables, sfiles_charts = _stage_files(out, sn)
        if sfiles_tables or sfiles_charts:
            stage_sections += '<h3 style="color:var(--ink-secondary);font-size:14px;margin-top:12px">Stage Results</h3>'
        if sfiles_tables:
            stage_sections += sfiles_tables
        if sfiles_charts:
            if sn == 2:
                stage_sections += sfiles_charts
            else:
                stage_sections += '<div class="chart-gallery">' + sfiles_charts + '</div>'
        # Stage 2: Custom filter stats with pre/post comparison
        if sn == 2:
            raw_tsv = out / os.path.join(_D['a_detect'], 'summary/all_viruses.raw.tsv')
            raw = safe_read_csv(raw_tsv)
            pre_sample = {}; pre_virus = {}
            if raw is not None:
                raw_acc = next((c for c in ["Accession","Rep_Accession"] if c in raw.columns), None)
                if raw_acc and "Sample" in raw.columns:
                    pre_sample = raw.groupby("Sample")[raw_acc].nunique().to_dict()
                    pre_virus = raw.groupby(raw_acc)["Sample"].nunique().to_dict()
            # Also build species->accession map from best TSV for per-virus lookup
            best_tsv = out / os.path.join(_D['a_detect'], 'summary/all_viruses.best.summary.tsv')
            best = safe_read_csv(best_tsv)
            sp2acc = {}
            if best is not None:
                bac = next((c for c in ["Rep_Accession","Accession"] if c in best.columns), None)
                spc = next((c for c in ["Adjusted_Species","Species_NCBI","Species"] if c in best.columns), None)
                if bac and spc:
                    for _, r in best.iterrows():
                        sp2acc[str(r[spc])] = str(r[bac])
            # Parse filter_stats
            ps_tsv = out / os.path.join(_D['a_filter'], 'filter_stats.per_sample.tsv')
            pv_tsv = out / os.path.join(_D['a_filter'], 'filter_stats.per_virus.tsv')
            ps_df = safe_read_csv(ps_tsv); pv_df = safe_read_csv(pv_tsv)
            # [P2-23a] Post-Filter 列须取过滤后统计: filter_stats.*.tsv 是过滤前写出 (旧版语义),
            # 过滤判定完成后追加的 .passed 表才是真 Post-Filter; 旧版无 .passed 时回退原表并注明
            ps_passed_df = safe_read_csv(out / os.path.join(_D['a_filter'], 'filter_stats.per_sample.passed.tsv'))
            post_map = None
            if (ps_passed_df is not None and "Sample" in ps_passed_df.columns
                    and "n_viruses" in ps_passed_df.columns):
                post_map = {str(r["Sample"]): int(_fnum(r["n_viruses"])) for _, r in ps_passed_df.iterrows()}
            # Per-sample table
            if ps_df is not None:
                post_hdr = "Post-Filter" if post_map is not None else "Post-Filter (pre-filter stats, legacy)"
                stage_sections += '<h3 style="color:var(--ink-secondary);font-size:14px;margin-top:12px">Per-Sample Filtering</h3>'
                stage_sections += f'<div class="tb-scroll"><table style="font-size:12px"><tr><th>Sample</th><th>Pre-Filter</th><th>{post_hdr}</th><th>Filtered</th><th>Retention</th><th>Total Reads</th><th>Mean Cov%</th></tr>'
                for _, r in ps_df.iterrows():
                    s = str(r["Sample"]); pre = pre_sample.get(s, 0)
                    if post_map is not None:
                        post = post_map.get(s, 0)   # .passed 表中无该样本 → 过滤后为 0 条
                    else:
                        post = int(r["n_viruses"])  # 回退: 原表实为过滤前统计 (表头已注明)
                    flt = pre - post
                    rate = f"{post/max(pre,1)*100:.0f}%" if isinstance(pre, (int,float)) else "?"
                    stage_sections += f'<tr><td>{s}</td><td>{pre}</td><td>{post}</td><td>{flt}</td><td>{rate}</td><td>{r["total_EM_reads"]}</td><td>{r["mean_coverage"]}%</td></tr>'
                stage_sections += '</table></div>'
            # Build per-virus max/avg/min from high_conf
            hc_tsv = out / os.path.join(_D['a_filter'], 'high_conf.summary.tsv')
            hc = safe_read_csv(hc_tsv)
            if hc is not None:
                # "-" 占位列字符串 dtype -> min()/max()/mean() 崩; 先数值化
                for _c in ("Rep_Coverage(%)", "Rep_MeanDepth", "Asm_EM_Reads", "Asm_TPM"):
                    if _c in hc.columns:
                        hc[_c] = pd.to_numeric(hc[_c], errors="coerce")
            virus_stats = {}
            if hc is not None:
                ac = next((c for c in ["Rep_Accession","Accession"] if c in hc.columns), None)
                sc = next((c for c in ["Adjusted_Species","Species_NCBI","Species"] if c in hc.columns), None)
                if ac:
                    def _stat(grp, col, how):
                        if col is None or col not in grp.columns: return 0
                        v = getattr(grp[col], how)()
                        return float(v) if pd.notna(v) else 0
                    for acc, grp in hc.groupby(ac):
                        sp = str(grp[sc].iloc[0])[:30] if sc else str(acc)
                        cov = "Rep_Coverage(%)" if "Rep_Coverage(%)" in hc.columns else None
                        depth = "Rep_MeanDepth" if "Rep_MeanDepth" in hc.columns else None
                        reads = "Asm_EM_Reads" if "Asm_EM_Reads" in hc.columns else None
                        tpm = "Asm_TPM" if "Asm_TPM" in hc.columns else None
                        virus_stats[acc] = {"species": sp,
                            "cov_min": _stat(grp, cov, "min"), "cov_max": _stat(grp, cov, "max"), "cov_avg": _stat(grp, cov, "mean"),
                            "depth_min": _stat(grp, depth, "min"), "depth_max": _stat(grp, depth, "max"), "depth_avg": _stat(grp, depth, "mean"),
                            "reads_min": _stat(grp, reads, "min"), "reads_max": _stat(grp, reads, "max"), "reads_avg": _stat(grp, reads, "mean"),
                            "tpm_min": _stat(grp, tpm, "min"), "tpm_max": _stat(grp, tpm, "max"), "tpm_avg": _stat(grp, tpm, "mean")}
            # Per-virus table
            if pv_df is not None:
                # [P2-23a] 同 Per-Sample: Post 列改读过滤后 .passed 表 (旧版无该表 → 回退原表并注明)
                pv_passed_df = safe_read_csv(out / os.path.join(_D['a_filter'], 'filter_stats.per_virus.passed.tsv'))
                post_virus_map = None
                if (pv_passed_df is not None and "Adjusted_Species" in pv_passed_df.columns
                        and "n_samples" in pv_passed_df.columns):
                    post_virus_map = {str(r["Adjusted_Species"]): int(_fnum(r["n_samples"])) for _, r in pv_passed_df.iterrows()}
                post_hdr_v = "Post" if post_virus_map is not None else "Post(pre-filter, legacy)"
                stage_sections += '<h3 style="color:var(--ink-secondary);font-size:14px;margin-top:16px">Per-Virus Filtering</h3>'
                stage_sections += f'<div class="tb-scroll"><table style="font-size:11px"><tr><th>Species</th><th>Pre</th><th>{post_hdr_v}</th><th>Flt</th><th>Ret</th>'
                for m in ["min","avg","max"]:
                    stage_sections += f'<th>Cov% {m}</th><th>Depth {m}</th><th>Reads {m}</th><th>TPM {m}</th>'
                stage_sections += '</tr>'
                for _, r in pv_df.iterrows():
                    sp = str(r["Adjusted_Species"])[:25]; acc = sp2acc.get(str(r["Adjusted_Species"]), "")
                    pre = pre_virus.get(acc, 0) if acc else 0
                    if post_virus_map is not None:
                        post = post_virus_map.get(str(r["Adjusted_Species"]), 0)  # .passed 表无该物种 → 过滤后 0
                    else:
                        post = int(r["n_samples"])  # 回退: 原表实为过滤前统计 (表头已注明)
                    flt = pre - post
                    rate = f"{post/max(pre,1)*100:.0f}%" if isinstance(pre,(int,float)) else "?"
                    vs = virus_stats.get(acc, {})
                    cells = f'<td>{sp}</td><td>{pre}</td><td>{post}</td><td>{flt}</td><td>{rate}</td>'
                    for m in ["min","avg","max"]:
                        c = f'{vs.get("cov_"+m,0):.1f}%' if vs else "-"
                        d = f'{vs.get("depth_"+m,0):.1f}' if vs else "-"
                        rd = f'{vs.get("reads_"+m,0):.0f}' if vs else "-"
                        tp = f'{vs.get("tpm_"+m,0):.0f}' if vs else "-"
                        cells += f'<td>{c}</td><td>{d}</td><td>{rd}</td><td>{tp}</td>'
                    stage_sections += f'<tr>{cells}</tr>'
                stage_sections += '</table></div>'
        # [P2-24g] Stage 2: filter_summary_plot 无条件内嵌 (原内嵌分支挂在 has_sub 循环内,
        # 而 S2 has_sub=False → 永不执行); png 优先, 缺失/超限时回退既有 pdf <object> 通道
        if sn == 2:
            fp_png = out / _D['a_filter'] / "filter_summary_plot.png"
            fp_pdf = out / _D['a_filter'] / "filter_summary_plot.pdf"
            s2_plot_html = ""
            if fp_png.is_file():
                b64_p = img_to_base64(fp_png, 3000)   # 超限 / --no_images → None
                if b64_p:
                    s2_plot_html = (f'<div class="chart-card"><div class="chart-title">Filter Summary</div>'
                                    f'<img src="data:image/png;base64,{b64_p}" loading="lazy" alt="Filter Summary"></div>')
            if not s2_plot_html and fp_pdf.is_file():
                b64_d = img_to_base64(fp_pdf, 1000)
                if b64_d:
                    s2_plot_html = (f'<div class="chart-card"><div class="chart-title">Filter Summary (PDF)</div>'
                                    f'<object data="data:application/pdf;base64,{b64_d}" type="application/pdf" width="100%" height="500px"><p>PDF</p></object></div>')
            if not s2_plot_html and fp_png.is_file():
                s2_plot_html = _img_placeholder("Filter Summary")
            if s2_plot_html:
                stage_sections += s2_plot_html
        # Stage 1: Detection table + co-infection (cards already inserted above before _stage_files)
        if sn == 1 and detection_summary:
            ds = detection_summary
            dv = ds.get("viruses", {})
            if dv:
                stage_sections += '<div class="tb-scroll"><table style="font-size:11px;margin:8px 0"><tr><th>Species</th><th>Samples</th><th>Avg Cov%</th><th>Avg Depth</th><th>Avg CPM</th><th>Avg FPKM</th><th>Poisson</th></tr>'
                for acc, v in sorted(dv.items(), key=lambda x: x[1].get("n_samples",0), reverse=True):
                    sp = _esc(v["species"][:30]); n = v["n_samples"]
                    cov = f'{v["avg_cov"]:.1f}%' if v.get("avg_cov") is not None else "-"
                    depth = f'{v["avg_depth"]:.1f}x' if v.get("avg_depth") is not None else "-"
                    cpm = f'{v["avg_cpm"]:.0f}' if v.get("avg_cpm") is not None else "-"
                    fpkm = f'{v["avg_fpkm"]:.0f}' if v.get("avg_fpkm") is not None else "-"
                    poisson = f'{v["avg_poisson"]:.2f}' if v.get("avg_poisson") is not None else "-"
                    stage_sections += f'<tr><td>{sp}</td><td>{n}</td><td>{cov}</td><td>{depth}</td><td>{cpm}</td><td>{fpkm}</td><td>{poisson}</td></tr>'
                stage_sections += '</table></div>'
            coi_tsv = out / os.path.join(_D['a_align'], 'summary/Coinfection_Matrix_Reads.tsv')
            if coi_tsv.is_file():
                try:
                    coi_rows = []
                    with open(str(coi_tsv),'r',encoding='utf-8',errors='replace') as fh:
                        for line in fh: coi_rows.append([c.strip() for c in line.rstrip().split('\t')])
                    if coi_rows and len(coi_rows) > 1:
                        stage_sections += '<details style="margin-top:8px"><summary style="font-weight:600;color:var(--accent);cursor:pointer;font-size:13px">Co-infection Matrix (reads per sample)</summary>'
                        stage_sections += '<div class="tb-scroll" style="margin-top:4px"><table style="font-size:11px"><tr>'
                        for c in coi_rows[0]: stage_sections += f'<th>{_esc(c[:15])}</th>'
                        stage_sections += '</tr>'
                        for r in coi_rows[1:]:
                            stage_sections += '<tr>'
                            for c in r: stage_sections += f'<td>{_esc(c[:12])}</td>'
                            stage_sections += '</tr>'
                        stage_sections += '</table></div></details>'
                except Exception: pass
        # Stage 1: Add batch_plot visuals (after cards/table/co-infection)
        if sn == 1:
            # viz 输出路径: 新编排器写 01_Detection/sample_distribution|coabundance,
            # batch_virus_depth 写 01_Detection/plots/ 下 —— 两个候选根都要查
            # [P1-12] 原实现只在"plots/ 存在时"二选一, 两处产物并存时会漏取顶层图; 改为逐图双根查找
            detect_dir = out / _D['a_detect']
            plots_roots = [detect_dir, detect_dir / 'plots']
            keep = {"sample_distribution": ["sample_virus_count_bar.png", "virus_occurrence_bar.png"], "coabundance": ["coabundance_cooccurrence_heatmap.png", "coabundance_spearman_heatmap.png"]}
            stage_sections += '<h3 style="color:var(--ink-secondary);font-size:14px;margin-top:12px">Sample Distribution &amp; Co-abundance</h3><div class="chart-gallery">'
            for bd, keep_files in keep.items():
                for fn in keep_files:
                    for root in plots_roots:   # 顶层优先, 命中即止 (防同一张图两处内嵌)
                        f = root / bd / fn
                        if not f.is_file(): continue
                        b64 = img_to_base64(str(f), 2000)
                        if b64: stage_sections += f'<div class="chart-card"><div class="chart-title">{_esc(fn)}</div><img src="data:image/png;base64,{b64}" loading="lazy" alt="{_esc(fn)}"></div>'
                        else: stage_sections += _img_placeholder(fn)
                        break
            stage_sections += '</div>'
            # [P2-24h] S1 depth_plots 零覆盖修复: batch_plot_virus_depth.py --mode depth 产物
            # (auto_known_virus.py 接线为 01_Detection/depth_plots/{sample}/*_depth.png);
            # 按样本分组折叠内嵌, 总上限 40 张, 超出仅列文件名 (避免巨型 HTML)
            depth_root = None
            for _root in (detect_dir / 'depth_plots', detect_dir / 'plots' / 'depth_plots'):
                if _root.is_dir():
                    depth_root = _root
                    break
            if depth_root is not None:
                by_sample = {}
                for sdir in sorted(depth_root.iterdir()):
                    if not sdir.is_dir(): continue
                    files = [f for f in sorted(sdir.glob("*_depth.png")) if f.is_file()]
                    if files: by_sample[sdir.name] = files
                if by_sample:
                    _DEPTH_CAP = 40
                    n_total = sum(len(v) for v in by_sample.values())
                    stage_sections += '<h3 style="color:var(--ink-secondary);font-size:14px;margin-top:12px">Sequencing Depth Plots</h3>'
                    stage_sections += (f'<p style="color:var(--ink-secondary);font-size:11px;margin:2px 0">'
                                       f'{n_total} per-virus depth plots across {len(by_sample)} samples '
                                       f'(max {_DEPTH_CAP} embedded; remaining listed by filename)</p>')
                    n_embedded = 0
                    overflow = []
                    for sname in sorted(by_sample):
                        stage_sections += f'<details style="margin:4px 0"><summary style="font-weight:600;color:var(--accent);cursor:pointer;font-size:13px">{_esc(sname)} ({len(by_sample[sname])} plots)</summary><div class="chart-gallery" style="margin-top:6px">'
                        for f in by_sample[sname]:
                            if n_embedded >= _DEPTH_CAP:
                                overflow.append(f"{sname}/{f.name}")
                                continue
                            b64 = img_to_base64(str(f), 2000)
                            if b64:
                                stage_sections += f'<div class="chart-card"><div class="chart-title">{_esc(f.name)}</div><img src="data:image/png;base64,{b64}" loading="lazy" alt="{_esc(f.name)}"></div>'
                            else:
                                stage_sections += _img_placeholder(f.name)
                            n_embedded += 1
                        stage_sections += '</div></details>'
                    if overflow:
                        stage_sections += f'<details style="margin:4px 0"><summary style="font-size:12px;color:var(--ink-secondary);cursor:pointer">Not embedded &mdash; {len(overflow)} plots over the {_DEPTH_CAP}-plot cap</summary>'
                        stage_sections += ''.join(f'<p style="font-size:11px;color:var(--ink-secondary);margin:2px 0">{_esc(o)}</p>' for o in overflow)
                        stage_sections += '</details>'
        # Stage 9: DVG per-virus statistics table
        if sn == 9:
            dvg_m = out / os.path.join(_D['a_dvg'], 'Summary_Analysis_Report/Matrix_Sample_Wise_Recombination_Statistics.csv')
            if dvg_m.is_file():
                try:
                    import re as _re2
                    mdf = safe_read_csv(dvg_m)
                    if mdf is not None and len(mdf) > 0:
                        a = {}; fn_col = mdf.columns[0]
                        for _, r in mdf.iterrows():
                            fn = str(r[fn_col])
                            mm = _re2.search(r'([A-Z]{2,4}_?\d{5,9}\.\d{1,2})', fn)
                            acc = mm.group(1) if mm else (fn.split('_')[1] if '_' in fn else '?')
                            if acc not in a: a[acc] = {"r":0,"j":0,"s":set()}
                            a[acc]["r"] += int(float(r.iloc[1])) if len(r)>1 else 0
                            a[acc]["j"] += int(float(r.iloc[2])) if len(r)>2 else 0
                            a[acc]["s"].add(fn.split('_')[0])
                        if a:
                            stage_sections += '<details style="margin-top:8px"><summary style="font-weight:600;color:var(--accent);cursor:pointer;font-size:13px">DVG Events per Virus</summary>'
                            stage_sections += '<div class="tb-scroll" style="margin-top:4px"><table style="font-size:11px"><tr><th>Virus</th><th>DVG Reads</th><th>Jump Events</th><th>Samples</th></tr>'
                            for acc in sorted(a.keys()):
                                v = a[acc]; stage_sections += f'<tr><td>{_esc(acc)}</td><td>{v["r"]}</td><td>{v["j"]}</td><td>{len(v["s"])}</td></tr>'
                            stage_sections += '</table></div></details>'
                except Exception: pass
        # Stage 5: Extract summary table with per-virus FASTA stats
        if sn == 5 and extract_stats:
            rows = ""
            for acc, data in sorted(viruses.items()):
                sp = data.get("species", acc)[:30]
                es = extract_stats.get(acc, {})
                if es:
                    # Show assembly source count for filtering context
                    asm_data = assembly_summary.get(acc, {})
                    n_asm = len(asm_data.get("samples", []))
                    filtered = n_asm - es["n_contigs"] if n_asm > es["n_contigs"] else 0
                    filter_str = f' (filtered {filtered})' if filtered > 0 else ""
                    rows += f'<tr><td>{_esc(sp)}</td><td>{_esc(acc)}</td><td>{es["n_contigs"]}{filter_str}</td><td>{es["avg_len"]:.0f}</td><td>{es["min_len"]:,}</td><td>{es["max_len"]:,}</td><td>{es["avg_n_pct"]:.1f}%</td></tr>'
                else:
                    asm_data = assembly_summary.get(acc, {})
                    n_asm = len(asm_data.get("samples", []))
                    note = f'0 (all {n_asm} filtered)' if n_asm > 0 else '-'
                    rows += f'<tr><td>{_esc(sp)}</td><td>{_esc(acc)}</td><td>{note}</td><td>-</td><td>-</td><td>-</td><td>-</td></tr>'
            stage_sections += '<div class="tb-scroll"><table style="font-size:12px;margin-top:8px"><tr><th>Species</th><th>Accession</th><th>Numbers</th><th>Avg bp</th><th>Min bp</th><th>Max bp</th><th>N%</th></tr>' + rows + '</table></div>'
        # [P2-24j] Stage 5: assembly_stats 图内嵌 (原先只由 copy_results_to_report 复制、不内嵌);
        # 阶段级单图 (非每病毒), png 优先, 仅 pdf 时走既有 <object> 通道 (阈值 1000KB),
        # 超限/缺 png 时保留原复制行为并给文件名占位卡
        if sn == 5:
            _ext_dir = _resolve_dir(out_dir, _D['a_extract'])
            asm_png = _ext_dir / "assembly_stats.png"
            asm_pdf = _ext_dir / "assembly_stats.pdf"
            s5_plot_html = ""
            if asm_png.is_file():
                b64_a = img_to_base64(str(asm_png), 2000)
                if b64_a:
                    s5_plot_html = (f'<div class="chart-card"><div class="chart-title">Assembly Stats</div>'
                                    f'<img src="data:image/png;base64,{b64_a}" loading="lazy" alt="assembly_stats.png"></div>')
            if not s5_plot_html and asm_pdf.is_file():
                b64_ap = img_to_base64(str(asm_pdf), 1000)
                if b64_ap:
                    s5_plot_html = (f'<div class="chart-card"><div class="chart-title">Assembly Stats (PDF)</div>'
                                    f'<object data="data:application/pdf;base64,{b64_ap}" type="application/pdf" width="100%" height="500px"><p>PDF</p></object></div>')
                else:
                    s5_plot_html = _img_placeholder("assembly_stats.pdf (copied to S5_Assembly_Stats/)")
            elif not s5_plot_html and asm_png.is_file():
                s5_plot_html = _img_placeholder("assembly_stats.png (copied to S5_Assembly_Stats/)")
            if s5_plot_html:
                stage_sections += '<h3 style="color:var(--ink-secondary);font-size:14px;margin-top:12px">Assembly Stats Plot</h3>'
                stage_sections += s5_plot_html
        if has_sub:
            stage_sections += '<div style="margin-left:20px">'
            for acc, data in sorted(viruses.items()):
                sp = data.get("species", acc)
                vsa = f"s{sn}-virus-{acc}"
                # Sidebar sub-item
                nav_items += f'<li class="nav-item"><a href="#{vsa}" class="nav-link" style="padding-left:40px;font-size:11px;color:var(--ink-secondary)">&bull; {_esc(sp[:20])}</a></li>'
                # Per-stage per-virus sub-section
                stage_sections += f'<section id="{vsa}" style="margin:8px 0">'
                stage_sections += f'<h3 style="color:#2980b9;font-size:16px">{_esc(sp)} <span style="color:var(--ink-secondary);font-size:12px">({_esc(acc)})</span> <button onclick="runVirusAI({sn},\'{_esc(acc)}\')" class="ai-stage-btn" title="AI summarize this virus">AI</button></h3>'
                # Stage-specific content
                if sn == 6:
                    if data.get("n_samples", 0) <= 1:
                        stage_sections += '<p style="color:var(--ink-secondary);font-size:12px">Single sample &mdash; post-hoc analysis skipped (requires &ge;2 samples).</p>'
                        stage_sections += '</section>'
                        continue
                    # Kruskal-Wallis dN/dS test (if available)
                    kw_csv = _d(_D['a_post']) / _find_virus_dir(_d(_D['a_post']), acc).name / "snpgenie/Stats_Kruskal_Wallis.csv" if _find_virus_dir(_d(_D['a_post']), acc) else None
                    if kw_csv and kw_csv.is_file():
                        try:
                            kw = safe_read_csv(kw_csv)
                            if kw is not None and len(kw) > 0:
                                p = float(kw.iloc[0].get("P_value", 1))
                                stat = float(kw.iloc[0].get("Statistic", 0))
                                genes = str(kw.iloc[0].get("Genes_Compared", ""))
                                sig = "significant" if p < 0.05 else "not significant"
                                stage_sections += f'<p style="font-size:11px;color:var(--ink-secondary);margin:4px 0">Kruskal-Wallis dN/dS across genes: H={stat:.1f}, p={p:.3f} ({sig})</p>'
                        except Exception: pass
                    # Variant calls summary table
                    vcf_viz_dir = _find_virus_dir(post_dir, acc)
                    if vcf_viz_dir:
                        var_csv = vcf_viz_dir / "vcf_viz/Table1_Unified_HighQuality_Variants.csv"
                        if var_csv.is_file():
                            try:
                                vdf = safe_read_csv(var_csv)
                                if vdf is not None and len(vdf) > 0:
                                    n_variants = len(vdf)
                                    ts = int((vdf["Molecular_Type"] == "Transition (Ts)").sum()) if "Molecular_Type" in vdf.columns else 0
                                    tv = int((vdf["Molecular_Type"] == "Transversion (Tv)").sum()) if "Molecular_Type" in vdf.columns else 0
                                    mean_af = vdf["ALT_FREQ"].mean() if "ALT_FREQ" in vdf.columns else 0
                                    stage_sections += f'<p style="font-size:11px;color:var(--ink-secondary);margin:4px 0">{n_variants} high-quality variants: {ts} transitions, {tv} transversions, mean AF={mean_af:.2f}</p>'
                                    # Show top 5 variants
                                    stage_sections += '<details style="margin:2px 0"><summary style="font-size:11px;color:var(--accent);cursor:pointer">Top variants</summary>'
                                    stage_sections += '<div class="tb-scroll"><table style="font-size:10px"><tr><th>Sample</th><th>Pos</th><th>Ref</th><th>Alt</th><th>AF</th><th>Type</th></tr>'
                                    cols = {"Sample_ID":"?", "POS":"?", "REF":"?", "ALT":"?", "ALT_FREQ":"?", "Molecular_Type":"?"}
                                    for _, r in vdf.head(8).iterrows():
                                        s = str(r.get("Sample_ID","?"))[:8]
                                        pos = str(r.get("POS","?"))
                                        ref = str(r.get("REF","?"))[:10]
                                        alt = str(r.get("ALT","?"))[:10]
                                        af = f'{float(r.get("ALT_FREQ",0))*100:.1f}%' if r.get("ALT_FREQ") else "?"
                                        mt = str(r.get("Molecular_Type","?"))[:20]
                                        stage_sections += f'<tr><td>{s}</td><td>{pos}</td><td>{ref}</td><td>{alt}</td><td>{af}</td><td>{mt}</td></tr>'
                                    stage_sections += '</table></div></details>'
                            except Exception: pass
                    charts = collect_charts(post_dir, acc)
                    if charts:
                        # Group by analysis type
                        s6_groups = [
                            ("vcf_viz", "Variant Landscape"),
                            ("snpeff_macro", "SnpEff Macro"),
                            ("maftools", "MAF Analysis"),
                            ("snpgenie", "SNPGenie"),
                            ("vcf_merge", "VCF Merge"),
                        ]
                        for prefix, group_label in s6_groups:
                            group_charts = {k: v for k, v in charts.items() if k.startswith(prefix)}
                            if not group_charts: continue
                            gid = f"s6-virus-{acc}-{prefix}"
                            nav_items += f'<li class="nav-item"><a href="#{gid}" class="nav-link" style="padding-left:52px;font-size:10px;color:var(--ink-secondary)">&ndash; {group_label}</a></li>'
                            stage_sections += f'<section id="{gid}" style="margin:4px 0;padding:6px 8px;border:1px solid var(--border-light);border-radius:4px">'
                            stage_sections += f'<h4 style="font-size:12px;color:var(--ink-secondary);margin:0 0 4px;font-weight:600">{group_label} <span style="font-weight:400;font-size:10px">({len(group_charts)} charts)</span></h4>'
                            stage_sections += '<div class="chart-gallery">'
                            for ci in sorted(group_charts.values(), key=lambda x: x.get('label','')):
                                if ci.get('ext') == 'txt':
                                    stage_sections += f'<div class="skip-chip">\u23f8 {_esc(ci.get("label",""))}</div>'
                                    continue
                                if ci.get('base64'):
                                    b64 = ci['base64']; ext = ci.get('ext','')
                                    if ext == 'pdf':
                                        stage_sections += f'<div class="chart-card"><div class="chart-title">{ci.get("label","")}</div><object data="data:application/pdf;base64,{b64}" type="application/pdf" width="100%" height="500px"><p>PDF</p></object></div>'
                                    else:
                                        stage_sections += f'<div class="chart-card"><div class="chart-title">{ci.get("label","")}</div><img src="data:image/png;base64,{b64}" loading="lazy" alt="{ci.get("label","")}"></div>'
                                else:
                                    stage_sections += _img_placeholder(ci.get("label",""))
                            stage_sections += '</div></section>'
                    else:
                        stage_sections += '<p style="color:var(--ink-secondary);font-size:12px">Post-hoc charts not available for this virus.</p>'
                elif sn == 9:
                    dcharts = collect_dvg_charts(out / _D['a_dvg'], acc)
                    if dcharts:
                        stage_sections += '<div class="chart-gallery">'
                        for ci in sorted(dcharts.values(), key=lambda x: x.get('label','')):
                            if ci.get('base64'):
                                b64 = ci['base64']; ext = ci.get('ext','')
                                if ext == 'pdf':
                                    stage_sections += f'<div class="chart-card"><div class="chart-title">{ci.get("label","")}</div><object data="data:application/pdf;base64,{b64}" type="application/pdf" width="100%" height="500px"><p>PDF</p></object></div>'
                                else:
                                    stage_sections += f'<div class="chart-card"><div class="chart-title">{ci.get("label","")}</div><img src="data:image/png;base64,{b64}" loading="lazy" alt="{ci.get("label","")}"></div>'
                            else:
                                stage_sections += _img_placeholder(ci.get("label",""))
                        stage_sections += '</div>'
                    else:
                        stage_sections += '<p style="color:var(--ink-secondary);font-size:12px">DVG charts not available for this virus.</p>'
                else:
                    # S3/S4/S7/S8: locate per-virus result files
                    if sn == 5:
                        # Count extracted contigs per virus
                        ext_dir = _resolve_dir(out_dir, _D['a_extract'])
                        vd = _find_virus_dir(ext_dir, acc)
                        if vd:
                            fasta_files = list(vd.glob("*.fasta")) + list(vd.glob("*.fa"))
                            n_contigs = len(fasta_files)
                            stage_sections += f'<p style="font-size:12px;color:#666">{n_contigs} cleaned contig(s) extracted (N-fill via reference alignment).</p>'
                        else:
                            stage_sections += '<p style="color:var(--ink-secondary);font-size:12px">No extracted contigs found for this virus.</p>'
                    elif sn == 3:
                        # Stage 3: Variants - aggregated metrics
                        vs = variant_summary.get(acc, {})
                        if vs:
                            pi_s = f"{vs['avg_pi']:.4f}" if vs.get('avg_pi') is not None else "-"
                            sh_s = f"{vs['avg_shannon']:.4f}" if vs.get('avg_shannon') is not None else "-"
                            cov_s = f"{vs['avg_cov']:.1f}%" if vs.get('avg_cov') is not None else "-"
                            depth_s = f"{vs['avg_depth']:.1f}x" if vs.get('avg_depth') is not None else "-"
                            ref_len = f"{vs['ref_length']:.0f} bp" if vs.get('ref_length') is not None else "-"
                            stage_sections += '<table style="font-size:12px;margin:8px 0"><tr><th>Samples</th><th>Ref Length</th><th>Avg Cov</th><th>Avg Depth</th><th>Avg Pi</th><th>Shannon</th></tr>'
                            stage_sections += f'<tr><td>{vs.get("samples","?")}</td><td>{ref_len}</td><td>{cov_s}</td><td>{depth_s}</td><td>{pi_s}</td><td>{sh_s}</td></tr></table>'
                        else:
                            stage_sections += '<p style="color:var(--ink-secondary);font-size:12px">No variant summary available.</p>'
                        # [P2-24k] 产物清单: S3 无可内嵌图 (VCF/consensus/SnpEff 注释均为文本产物),
                        # 列出本病毒的实际文件名, 消除"该阶段在报告里无实体产物痕迹"
                        # [阶段重排] virus-consensus 在 03_Extract (a_align); VCF/SnpEff 在 04_Variants
                        artifacts = []
                        _ART_CAP = 15
                        for base_key, sub_rel, globs in (
                            ('a_variants', "virus-variants", ("*.filtered.vcf", "*.variants.vcf")),
                            ('a_align',    "virus-consensus", ("*.consensus.fasta", "*.fasta", "*.fa")),
                            ('a_variants', "virus-SnpEff", ("*.ann.vcf", "*.annotation_summary.tsv")),
                        ):
                            if len(artifacts) >= _ART_CAP: break
                            subd = out / _D[base_key] / sub_rel
                            if not subd.is_dir(): continue
                            vd = _find_virus_dir(subd, acc)
                            if not vd: continue
                            _seen = set()
                            for g in globs:
                                for f in sorted(vd.rglob(g)):
                                    if len(artifacts) >= _ART_CAP: break
                                    if not f.is_file() or f in _seen: continue
                                    _seen.add(f)
                                    artifacts.append(f"{sub_rel}/{vd.name}/{f.relative_to(vd).as_posix()}")
                        if artifacts:
                            stage_sections += f'<details style="margin:6px 0"><summary style="font-weight:600;color:var(--accent);cursor:pointer;font-size:13px">Stage artifacts ({len(artifacts)}) &mdash; VCF / consensus / SnpEff</summary>'
                            for a in artifacts:
                                stage_sections += f'<p style="font-size:11px;color:var(--ink-secondary);margin:2px 0"><code>{_esc(a)}</code></p>'
                            stage_sections += '</details>'
                    elif sn == 4:
                        # Stage 4: Assembly - per-sample stage evolution
                        asm_data = assembly_summary.get(acc, {})
                        samples = asm_data.get("samples", [])
                        valid = [s for s in samples if s.get("final_len")]
                        if valid:
                            ref_len = variant_summary.get(acc, {}).get("ref_length")
                            stage_sections += '<table style="font-size:11px;margin:8px 0"><tr><th>Sample</th><th>DeNovo</th><th>Ultimate</th><th>Growth</th><th>vs Ref</th><th>N50</th><th>N</th></tr>'
                            for s in sorted(valid, key=lambda x: x["name"]):
                                dlen = s["denovo_len"]; flen = s["final_len"]
                                if dlen > 0:
                                    diff = flen - dlen
                                    pct = (flen/dlen - 1) * 100
                                    sign = "+" if diff >= 0 else ""
                                    growth = f"{sign}{diff:.0f} bp ({sign}{pct:.0f}%)"
                                else:
                                    growth = "-"
                                vs_ref = f"{flen/ref_len*100:.1f}%" if ref_len and ref_len > 0 else "-"
                                n50 = f"{s['final_n50']:.0f}" if s.get("final_n50") else "-"
                                n = s.get("final_n", 0)
                                n_cls = ' style="color:#e74c3c;font-weight:bold"' if n > 0 else ''
                                n_str = f'<span{n_cls}>{n}</span>' if n > 0 else '0'
                                stage_sections += f'<tr><td>{_esc(s["name"])}</td><td>{dlen:.0f}</td><td>{flen:.0f}</td><td>{growth}</td><td>{vs_ref}</td><td>{n50}</td><td>{n_str}</td></tr>'
                            stage_sections += '</table>'
                            # Show stage evolution detail (collapsed) if multi-stage
                            if len(valid) > 0 and len(valid[0].get("stages", [])) > 2:
                                # Collect all unique stage names in order
                                all_steps = []
                                seen = set()
                                for s in valid:
                                    for st in s.get("stages", []):
                                        if st["step"] not in seen:
                                            all_steps.append(st["step"])
                                            seen.add(st["step"])
                                if len(all_steps) >= 3:
                                    stage_sections += '<details style="margin:4px 0"><summary style="font-size:11px;color:#2980b9;cursor:pointer">Stage detail (click to expand)</summary>'
                                    stage_sections += '<div class="tb-scroll"><table style="font-size:10px;margin:4px 0"><tr><th>Sample</th>'
                                    for step_name in all_steps:
                                        stage_sections += f'<th>{step_name}</th>'
                                    stage_sections += '</tr>'
                                    for s in valid:
                                        stage_sections += f'<tr><td>{_esc(s["name"])}</td>'
                                        stage_map = {st["step"]: st["length"] for st in s.get("stages", [])}
                                        for step_name in all_steps:
                                            val = stage_map.get(step_name)
                                            stage_sections += f'<td>{val:.0f}</td>' if val is not None else '<td>-</td>'
                                        stage_sections += '</tr>'
                                    stage_sections += '</table></div></details>'
                        else:
                            stage_sections += '<p style="color:var(--ink-secondary);font-size:12px">No assembly data found for this virus.</p>'
                    elif sn == 7:
                        # Stage 7: Capheine - BUSTED per-gene selection summary
                        cap = _collect_capheine_summary(out_dir, acc)
                        if cap and cap.get("genes"):
                            genes = cap["genes"]
                            sig = cap["sig_count"]
                            stage_sections += f'<p style="font-size:12px;color:#666">{cap["total_genes"]} genes analyzed. '
                            if sig > 0:
                                stage_sections += f'<span style="color:#e74c3c">{sig} genes show positive selection (BUSTED p&lt;0.05).</span>'
                            else:
                                stage_sections += 'No significant positive selection detected (BUSTED p&gt;0.05).'
                            stage_sections += '</p>'
                            stage_sections += '<div class="tb-scroll"><table style="font-size:11px;margin:8px 0"><tr><th>Gene</th><th>N</th><th>&omega;3</th><th>p-value</th><th>+</th><th>&minus;</th><th>Sites</th></tr>'
                            for g in genes:
                                oc = ' style="color:#e74c3c;font-weight:bold"' if g["omega3"] > 1 else ''
                                pc = ' style="color:#27ae60;font-weight:bold"' if g["pval"] < 0.05 else ''
                                pv = f'{g["pval"]:.2e}' if g["pval"] < 0.01 else f'{g["pval"]:.3f}'
                                stage_sections += f'<tr><td><strong>{_esc(g["gene"])}</strong></td><td>{g["n_seq"]}</td><td{oc}>{g["omega3"]:.3f}</td><td{pc}>{pv}</td><td>{g["positive"]}</td><td>{g["negative"]}</td><td>{g["total_sites"]}</td></tr>'
                            stage_sections += '</table></div>'
                        else:
                            if (out / "07_capheine").is_dir():
                                stage_sections += '<p style="color:var(--ink-secondary);font-size:12px">No coding sequences available (non-coding RNA virus) or capheine not run.</p>'
                            else:
                                stage_sections += '<p style="color:var(--ink-secondary);font-size:12px">CAPHEINE has moved to virome_phylo_pipeline (pipeline 5, select module); positive selection results are no longer included in this report.</p>'
                    elif sn == 8:
                        # Stage 8: Similarity - pairwise identity heatmaps
                        sim = _collect_similarity_data(out_dir, acc)
                        if sim and sim.get("available"):
                            s = sim.get("stats", {})
                            if s:
                                stage_sections += f'<p style="font-size:12px;color:#666">Pairwise NT identity: mean {s["mean_id"]:.1f}%, range {s["min_id"]:.1f}%&ndash;{s["max_id"]:.1f}% ({s["n_pairs"]} pairs, {s["n_samples"]} samples)</p>'
                            stage_sections += '<div class="chart-gallery">'
                            for label, b64 in sim.get("images", {}).items():
                                stage_sections += f'<div class="chart-card"><div class="chart-title">{label}</div><img src="data:image/png;base64,{b64}" loading="lazy" alt="{label}"></div>'
                            stage_sections += '</div>'
                        else:
                            stage_sections += '<p style="color:var(--ink-secondary);font-size:12px">Insufficient samples for pairwise similarity (requires &ge;2 samples).</p>'
                stage_sections += '</section>'
            stage_sections += '</div>'
            nav_items += '</div>'  # close stage-sub
        stage_sections += '</section>'

    nav_items += '</div>'  # close Pipeline Stages nav-sub
    nav_items += '<li class="nav-header" style="padding:8px 20px 4px;font-size:12px;color:var(--ink-secondary);font-weight:600;border-top:1px solid var(--border);margin-top:4px">Virus Summary</li>'
    nav_items += '<div class="nav-sub">'
    nav_items += '<li class="nav-item"><a href="#virus-summary" class="nav-link" style="font-size:11px;padding-left:24px;color:var(--sidebar-ink)">Metrics</a></li>'
    nav_items += '<li class="nav-item"><a href="#paper-reference" class="nav-link" style="font-size:11px;padding-left:24px;color:var(--sidebar-ink)">Paper Reference</a></li>'
    nav_items += '</div>'  # close Virus Summary nav-sub

    funnel_section = _build_funnel_section(out)
    html = f"""<!DOCTYPE html><html lang="en"><head><meta charset="UTF-8"><meta name="viewport" content="width=device-width, initial-scale=1.0">
<!-- impeccable-disable numbered-section-markers,em-dash-overuse: stage IDs are pipeline identifiers; CJK scientific copy -->
<title>Virus Pipeline Report</title>
<style>
/* ═══════════════════════════════════════════
   Virus Pipeline Report — Product Register
   Palette: Data-dense dashboard (blue data + amber highlights)
   Type: Charter display / Fira Sans body / Fira Code data
   ═══════════════════════════════════════════ */
@import url('https://fonts.googleapis.com/css2?family=Fira+Code:wght@400;500;600&family=Fira+Sans:wght@400;500;600;700&display=swap');
:root {{
  --ink: #16233d;
  --ink-secondary: #475569;
  --accent: #1e40af;
  --accent-hover: #16368f;
  --accent-subtle: #dbeafe;
  --accent-soft: #eff4ff;
  --accent-2: #0e7f8c;
  --surface: #ffffff;
  --surface-alt: #F8FAFC;
  --surface-hover: #edf1f7;
  --border: #dce3ee;
  --border-light: #ecf0f5;
  --sidebar-bg: #f3f5f8;
  --sidebar-ink: #4a5568;
  --sidebar-active: #1e40af;
  --success: #0d9488;
  --warning: #d97706;
  --danger: #dc2626;
  --radius: 8px;
  --shadow-sm: 0 1px 2px rgba(22,35,61,.05);
  --shadow-md: 0 6px 18px rgba(22,35,61,.09);
  --font-display: Charter, 'Bitstream Charter', 'Source Serif Pro', Georgia, 'Times New Roman', serif;
  --font-mono: 'Fira Code', ui-monospace, 'Cascadia Mono', Consolas, monospace;
}}

* {{box-sizing:border-box;margin:0;padding:0}}

body {{
  font-family:'Fira Sans', -apple-system, BlinkMacSystemFont, 'Segoe UI', system-ui, sans-serif;
  color: var(--ink); background: var(--surface-alt); display:flex; min-height:100vh;
  -webkit-font-smoothing: antialiased;
}}

/* ── Sidebar ── */
.sidebar {{
  position:fixed; left:0; top:0; width:248px; height:100vh;
  background: var(--sidebar-bg); color: var(--sidebar-ink);
  overflow-y:auto; padding:20px 0; z-index:100;
  border-right: 1px solid var(--border);
}}
.sidebar h3 {{
  padding:0 20px 12px; color: var(--ink); font-size:16px; font-weight:700;
  border-bottom: 1px solid var(--border); margin-bottom:8px; letter-spacing: -.01em;
}}
.nav-item {{list-style:none}}
.nav-link {{
  display:block; padding:7px 20px; color: var(--sidebar-ink);
  text-decoration:none; font-size:13px; transition: background 0.15s, color 0.15s;
  border-radius: 0 20px 20px 0; margin-right:8px;
}}
.nav-link:hover {{background: var(--surface-hover); color: var(--ink)}}
.nav-link.active {{background: var(--accent-subtle); color: var(--accent); font-weight:600}}

.main {{
  margin-left:248px; padding:32px 40px; flex:1; max-width:1440px;
}}

/* ── Typography ── */
h1 {{
  font-family: var(--font-display);
  color: var(--ink); font-size:31px; font-weight:700; margin-bottom:4px;
  padding-bottom:12px; border-bottom: 2px solid var(--border); letter-spacing: -.015em;
}}
h2 {{
  color: var(--ink); font-size:19px; font-weight:650; margin:28px 0 8px; letter-spacing: -.01em;
  position:relative; padding-left:14px;
}}
h2::before {{content:'';position:absolute;left:0;top:4px;bottom:4px;width:4px;border-radius:2px;background:linear-gradient(180deg,var(--accent),#7aa5f8)}}
h3 {{font-size:15px; font-weight:600; color: var(--ink); margin:16px 0 6px}}

/* ── Metrics cards ── */
.metrics {{display:flex; flex-wrap:wrap; gap:10px; margin:12px 0}}
.metric {{
  background: var(--surface); border:1px solid var(--border);
  border-radius: var(--radius); padding:12px 18px; text-align:center; min-width:88px;
  transition: box-shadow 0.15s;
}}
.metric:hover {{box-shadow: var(--shadow-sm)}}
.metric .value {{display:block; font-size:22px; font-weight:700; color: var(--accent); font-variant-numeric: tabular-nums}}
.metric .unit {{display:block; font-size:10.5px; color: var(--ink-secondary); margin-top:2px}}

.card-row {{display:flex; flex-wrap:wrap; gap:10px; margin:12px 0}}
.card {{
  background: var(--surface); border:1px solid var(--border);
  border-radius: var(--radius); padding:14px 18px; text-align:center; min-width:100px; flex:1;
  transition: box-shadow 0.15s;
}}
.card:hover {{box-shadow: var(--shadow-sm)}}
.card .value {{display:block; font-size:24px; font-weight:700; color: var(--accent); word-break:break-all; font-variant-numeric: tabular-nums}}
.card .label {{display:block; font-size:10px; color: var(--ink-secondary); margin-top:2px; font-weight:500}}

/* ── v21 additions: hero flow strip, status chips, print ── */
.flow-strip {{display:flex;align-items:stretch;gap:0;margin:16px 0 4px;flex-wrap:wrap;row-gap:6px}}
.flow-node {{
  display:inline-flex;flex-direction:column;justify-content:center;gap:1px;
  padding:7px 13px;background:var(--surface);border:1px solid var(--border);
  border-radius:5px;font-size:10px;font-weight:600;color:var(--ink-secondary);
  box-shadow:var(--shadow-sm); transition:border-color .15s, transform .15s;
}}
.flow-node b {{font-family:var(--font-display);font-size:12.5px;color:var(--accent);letter-spacing:.03em;line-height:1.1}}
.flow-node small {{font-weight:400;font-size:9px;color:var(--ink-secondary);opacity:.9;line-height:1.2}}
.flow-node:hover {{border-color:var(--accent);transform:translateY(-1px)}}
.flow-arrow {{align-self:center;color:#b9c2cf;padding:0 6px;font-size:13px;user-select:none}}
.skip-chip {{
  display:flex;align-items:center;gap:8px;padding:11px 14px;margin:4px 0;
  background:#fffbeb;border:1px solid #fde68a;
  border-radius:var(--radius);font-size:12px;line-height:1.5;color:#92400e;
}}
@media print {{
  .sidebar,#back-to-top,#ai-panel,#ai-btn,.pg-nav,.chart-dl-btn,.ai-stage-btn{{display:none!important}}
  .main{{margin-left:0;padding:0;max-width:none}}
  body{{display:block;background:#fff}}
  section{{page-break-inside:avoid}}
}}

/* ── v22: a11y, rendering perf, data-typography ── */
:focus-visible {{outline:2px solid var(--accent); outline-offset:2px; border-radius:2px}}
::selection {{background:#bfd7ff; color:var(--ink)}}
h1, h2 {{text-wrap:balance}}
table {{font-variant-numeric: tabular-nums}}
.tb-scroll thead th {{position:sticky; top:0; z-index:2; background:#f1f5f9; box-shadow:inset 0 -1px 0 var(--border)}}
.main > section {{content-visibility:auto; contain-intrinsic-size:auto 640px}}
.sidebar::-webkit-scrollbar {{width:8px}}
.sidebar::-webkit-scrollbar-thumb {{background:#cdd6e2; border-radius:4px; border:2px solid var(--sidebar-bg)}}
.sidebar::-webkit-scrollbar-thumb:hover {{background:#b3c1d4}}
.flow-node b, .card .value {{font-feature-settings:'tnum' 1}}

/* ── Chart gallery ── */
.chart-gallery {{
  display:grid; grid-template-columns: 1fr;
  gap:12px; margin:12px 0;
}}
.chart-card {{
  background: var(--surface); border:1px solid var(--border);
  border-radius: var(--radius); overflow:hidden;
  transition: box-shadow 0.18s, transform 0.18s, border-color 0.18s;
}}
.chart-card:hover {{box-shadow: var(--shadow-md); transform:translateY(-1px); border-color:#c9d4e4}}
.chart-title {{
  background:#fafbfd; padding:7px 11px; font-size:11.5px;
  font-weight:600; color: var(--ink-secondary); border-bottom: 1px solid var(--border-light);
  display:flex; align-items:center; gap:6px;
}}
.chart-title::before {{content:'';width:6px;height:6px;border-radius:50%;background:linear-gradient(135deg,var(--accent),#7aa5f8);flex-shrink:0}}
.chart-card img {{
  width:100%; height:auto; max-height:80vh; object-fit:contain; display:block; background:#fff;
}}
.chart-placeholder {{padding:40px; text-align:center; background: var(--surface-alt)}}
.chart-placeholder a {{color: var(--accent); text-decoration:none; font-size:13px}}

/* ── Tables ── */
table {{border-collapse:collapse; width:100%; font-size:12px}}
th,td {{padding:7px 10px; text-align:left; border-bottom:1px solid var(--border-light)}}
th {{
  background: var(--surface-alt); font-size:11px; color: var(--ink-secondary);
  font-weight:600; position:sticky; top:0; z-index:1;
}}
td {{color: var(--ink)}}
tbody tr:nth-child(even) td {{background: #f9fafb}}
tbody tr:hover td {{background: var(--accent-subtle)}}
.tb-scroll {{
  overflow-x:auto; max-height:500px; overflow-y:auto;
  border:1px solid var(--border); border-radius: var(--radius);
}}

/* ── Per-virus sub-sections (no side-stripe!) ── */
section[id^="s"][id*="-virus-"] {{
  border: 1px solid var(--border);
  border-radius: var(--radius);
  padding: 10px 14px; margin: 8px 0 10px;
  background: var(--surface);
  transition: box-shadow 0.15s;
}}
section[id^="s"][id*="-virus-"]:hover {{box-shadow: var(--shadow-sm)}}
section[id^="s"][id*="-virus-"] h3 {{font-size:13.5px; margin:0 0 4px; color: var(--ink); font-weight:600}}

/* ── Paper Reference ── */
.paper-block {{
  background: var(--surface); border:1px solid var(--border);
  border-radius: var(--radius); padding:14px 18px; margin:10px 0;
}}
.paper-block table {{font-size:12px}}
.paper-block ul {{padding-left:20px; line-height:1.8; color: var(--ink-secondary)}}
.suggestion {{
  background: #fefdf5; border:1px solid #f0d060; border-radius: var(--radius);
  padding:10px 14px; margin:8px 0; font-size:12.5px; line-height:1.7; color: var(--ink-secondary);
}}

/* ── Collapsible & Misc ── */
details {{margin:6px 0}}
summary {{cursor:pointer; font-weight:600; color: var(--accent); padding:4px 0; font-size:12.5px}}
section {{margin-bottom:22px; padding-top:8px; scroll-margin-top:20px}}
.footer {{
  margin-top:40px; padding:20px; text-align:center;
  color: var(--ink-secondary); font-size:11.5px; border-top:1px solid var(--border);
}}

/* ── Back to top ── */
#back-to-top {{
  position:fixed; bottom:24px; right:24px; width:38px; height:38px;
  background: var(--accent); color:#fff; border:none; border-radius:50%;
  font-size:16px; cursor:pointer; opacity:0; transform:translateY(10px);
  transition: opacity 0.2s, transform 0.2s, background 0.15s;
  z-index:200; box-shadow: var(--shadow-md);
}}
#back-to-top.visible {{opacity:1; transform:translateY(0)}}
#back-to-top:hover {{background: var(--accent-hover)}}

/* ── Sidebar collapsible groups ── */
.nav-header {{
  cursor:pointer; user-select:none;
  display:flex; justify-content:space-between; align-items:center;
  font-weight:600; color: var(--ink-secondary);
}}
.nav-header::after {{
  content:'\\25B2'; font-size:8px; transition:transform 0.15s; margin-left:auto;
}}
.nav-header.collapsed::after {{transform:rotate(180deg)}}
.nav-sub {{overflow:hidden; max-height:2000px}}
.nav-sub.collapsed,.stage-sub {{max-height:0}}
.stage-toggle::after {{content:' \\25BC';font-size:7px;margin-left:4px;opacity:.4}}
.stage-toggle.expanded::after {{content:' \\25B2'}}

/* ── Lightbox ── */
.lightbox {{display:none;position:fixed;top:0;left:0;width:100vw;height:100vh;background:rgba(0,0,0,.88);z-index:9999;justify-content:center;align-items:center;flex-direction:column}}
.lightbox.show {{display:flex}}
.lightbox img {{max-width:92vw;max-height:82vh;object-fit:contain;border-radius:4px;box-shadow:0 8px 32px rgba(0,0,0,.4)}}
.lightbox-close {{position:absolute;top:16px;right:24px;color:#fff;font-size:32px;cursor:pointer;line-height:1;opacity:.7;transition:opacity .15s;z-index:10000;background:none;border:none}}
.lightbox-close:hover {{opacity:1}}
.lightbox-dl {{margin-top:12px;padding:8px 18px;background:rgba(255,255,255,.15);color:#fff;border:1px solid rgba(255,255,255,.25);border-radius:4px;font-size:13px;cursor:pointer;transition:background .15s;z-index:10000}}
.lightbox-dl:hover {{background:rgba(255,255,255,.25)}}
.lightbox-caption {{color:#ccc;font-size:12px;margin-top:8px;max-width:90vw;text-align:center;z-index:10000}}
/* Download button on chart cards */
.chart-card {{position:relative}}
.chart-dl-btn {{position:absolute;top:4px;right:8px;background:rgba(0,0,0,.06);border:none;border-radius:3px;padding:2px 8px;font-size:10px;cursor:pointer;color:#555;z-index:5;transition:background .15s;display:none}}
.chart-card:hover .chart-dl-btn {{display:block}}
.chart-dl-btn:hover {{background:rgba(0,0,0,.12);color:#222}}
/* ── AI Summary ── */
.ai-section {{margin:16px 0;border:1px solid var(--border);border-radius:var(--radius);overflow:hidden}}
.ai-header {{display:flex;align-items:center;gap:10px;padding:10px 14px;background:linear-gradient(135deg,#667eea22,#764ba222);border-bottom:1px solid var(--border)}}
.ai-title {{font-weight:700;font-size:14px;color:#5a4fcf}}
.ai-model-tag {{font-size:10px;color:var(--ink-secondary);background:var(--surface-alt);padding:2px 8px;border-radius:10px}}
.ai-brief {{padding:12px 14px;font-size:13px;line-height:1.7;color:var(--ink)}}
.ai-detail-toggle {{border-top:1px solid var(--border-light)}}
.ai-detail-toggle summary {{cursor:pointer;font-size:12px;font-weight:600;color:var(--accent);padding:8px 14px;background:var(--surface-alt);user-select:none}}
.ai-detail-toggle summary:hover {{background:var(--surface-hover)}}
.ai-detail-content {{padding:8px 14px 12px}}
.ai-block {{margin:8px 0}}
.ai-block-title {{font-weight:600;font-size:11px;color:var(--accent);margin-bottom:2px}}
.ai-block-text {{font-size:12px;line-height:1.7;color:var(--ink-secondary)}}
.ai-error {{padding:12px 14px;color:var(--ink-secondary);font-size:12px;font-style:italic;background:var(--surface-alt);border-radius:var(--radius)}}
/* ── In-page AI Panel ── */
#ai-btn {{position:fixed;top:16px;right:24px;width:36px;height:36px;background:linear-gradient(135deg,#667eea,#764ba2);color:#fff;border:none;border-radius:50%;font-size:11px;font-weight:700;cursor:pointer;z-index:300;box-shadow:0 2px 8px rgba(102,126,234,.3);transition:transform .15s}}
#ai-btn:hover {{transform:scale(1.1)}}
#ai-panel {{position:fixed;top:60px;right:24px;width:300px;max-height:80vh;overflow-y:auto;background:var(--surface);border:1px solid var(--border);border-radius:var(--radius);box-shadow:0 8px 24px rgba(0,0,0,.12);z-index:299}}
.ai-stage-summary {{font-size:11px;line-height:1.6;color:var(--ink);padding:6px 10px;margin:4px 0;background:var(--accent-subtle);border-radius:4px;display:none}}
.ai-copy-btn {{font-size:9px;padding:1px 6px;margin-left:8px;background:var(--surface);border:1px solid var(--border);border-radius:3px;color:var(--ink-secondary);cursor:pointer;vertical-align:middle}}
.ai-copy-btn:hover {{background:var(--accent-subtle);color:var(--accent)}}
.ai-stage-btn {{font-size:10px;font-weight:700;padding:1px 8px;margin-left:8px;background:linear-gradient(135deg,#667eea,#764ba2);color:#fff;border:none;border-radius:10px;cursor:pointer;vertical-align:middle;opacity:.8;transition:opacity .15s}}
.ai-stage-btn:hover {{opacity:1}}
/* ── Chat Q&A ── */
#chat-btn {{position:fixed;bottom:72px;right:24px;width:38px;height:38px;background:linear-gradient(135deg,#667eea,#764ba2);color:#fff;border:none;border-radius:50%;font-size:16px;font-weight:700;cursor:pointer;z-index:200;box-shadow:0 2px 8px rgba(102,126,234,.3);transition:transform .15s}}
#chat-btn:hover {{transform:scale(1.1)}}
#chat-panel {{position:fixed;bottom:120px;right:24px;width:360px;max-height:500px;background:var(--surface);border:1px solid var(--border);border-radius:var(--radius);box-shadow:0 8px 24px rgba(0,0,0,.12);z-index:199;display:flex;flex-direction:column}}
.chat-header {{display:flex;justify-content:space-between;align-items:center;padding:8px 12px;border-bottom:1px solid var(--border);font-size:13px;font-weight:600;color:var(--ink)}}
.chat-messages {{flex:1;overflow-y:auto;padding:8px 10px;max-height:360px;display:flex;flex-direction:column;gap:6px}}
.chat-msg {{font-size:12px;line-height:1.5;padding:6px 10px;border-radius:8px;max-width:85%}}
.chat-msg.user {{align-self:flex-end;background:var(--accent);color:#fff}}
.chat-msg.ai {{align-self:flex-start;background:var(--surface-alt);color:var(--ink)}}
.chat-input-row {{display:flex;gap:6px;padding:8px 10px;border-top:1px solid var(--border-light)}}
.chat-input-row input {{flex:1;padding:6px 8px;border:1px solid var(--border);border-radius:4px;font-size:12px}}
.chat-send-btn {{padding:6px 12px;background:var(--accent);color:#fff;border:none;border-radius:4px;font-size:12px;cursor:pointer}}
/* ── Reduced motion ── */
@media (prefers-reduced-motion: reduce) {{
  *, *::after {{transition-duration:0s!important; animation-duration:0s!important}}
}}

/* ── Responsive（2026-10-17 用户裁定：侧栏常驻，不再随宽度隐藏）── */
@media (max-width:560px) {{
  #back-to-top{{right:12px; bottom:12px}}
}}
@media print {{
  .sidebar{{display:none}} #back-to-top{{display:none}} .main{{margin-left:0}}
}}
</style>
<script>
// ── pg_tbl: TSV table pagination ──
function pg_tbl(tid,step){{
    var t=document.getElementById(tid);
    if(!t)return;
    var rows=t.getElementsByTagName('tbody')[0].getElementsByTagName('tr');
    var info=document.getElementById(tid+'_info');
    var pp=Math.abs(step);var cur=0;
    for(var i=0;i<rows.length;i++){{if(rows[i].style.display!=='none'){{cur=i;break;}}}}
    var page=Math.floor(cur/pp);
    var npage=page+(step>0?1:-1);
    var total=Math.ceil(rows.length/pp);
    if(npage<0)npage=0;if(npage>=total)npage=total-1;
    for(var i=0;i<rows.length;i++){{
        var p=Math.floor(i/pp);
        rows[i].style.display=(p===npage)?'':'none';
    }}
    if(info)info.textContent='Page '+(npage+1)+' of '+total;
}}
// ── Sidebar group collapse ──
document.addEventListener('DOMContentLoaded',function(){{
    var headers=document.querySelectorAll('.nav-header');
    headers.forEach(function(h){{
        h.addEventListener('click',function(){{
            this.classList.toggle('collapsed');
            var sub=this.nextElementSibling;
            while(sub&&sub.classList.contains('nav-sub')){{
                sub.classList.toggle('collapsed');
                sub=sub.nextElementSibling;
            }}
        }});
    }});
    // Stage sub-item toggles
    document.querySelectorAll('.stage-toggle').forEach(function(link){{
        link.addEventListener('click',function(e){{
            e.preventDefault();
            var sub=this.parentElement.nextElementSibling;
            if(sub&&sub.classList.contains('stage-sub')){{
                var isCollapsed=sub.style.maxHeight==='0px'||!sub.style.maxHeight;
                sub.style.maxHeight=isCollapsed?'2000px':'0px';
                this.classList.toggle('expanded',isCollapsed);
            }}
            // Navigate after short delay
            var href=this.getAttribute('href');
            if(href){{setTimeout(function(){{window.location=href;}},50);}}
        }});
    }});
}});
// ── Back to top ──
window.addEventListener('scroll',function(){{
    var btn=document.getElementById('back-to-top');
    if(!btn)return;
    if(window.scrollY>400)btn.classList.add('visible');
    else btn.classList.remove('visible');
}});
function scrollTop(){{window.scrollTo({{top:0,behavior:'smooth'}});}}
// ── Active nav highlight on scroll ──
window.addEventListener('scroll',function(){{
    var links=document.querySelectorAll('.nav-link');
    var fromTop=window.scrollY+100;
    var current='';
    links.forEach(function(l){{
        var id=l.getAttribute('href');
        if(!id||id[0]!=='#')return;
        var sec=document.getElementById(id.substring(1));
        if(sec&&sec.offsetTop<=fromTop)current=id;
    }});
    links.forEach(function(l){{
        l.classList.toggle('active',l.getAttribute('href')===current);
    }});
}});
</script></head><body>
<button id="back-to-top" onclick="scrollTop()" title="Back to top" aria-label="Back to top">&#9650;</button>

<button id="chat-btn" onclick="toggleChat()" title="Ask questions about this report" aria-label="Open report Q&amp;A">?</button>
<div id="chat-panel" style="display:none">
    <div class="chat-header">
        <span>Report Q&amp;A</span>
        <button onclick="toggleChat()" style="background:none;border:none;color:var(--ink-secondary);font-size:16px;cursor:pointer">&times;</button>
    </div>
    <div class="chat-messages" id="chat-msgs"></div>
    <div class="chat-input-row">
        <input id="chat-input" type="text" placeholder="Ask about this report…" aria-label="Ask a question about this report" onkeydown="if(event.key==='Enter')sendChat()">
        <button onclick="sendChat()" class="chat-send-btn">Send</button>
    </div>
</div>
<nav class="sidebar">
    <h3>Pipeline Report</h3>
    <ul>{nav_items}</ul>
</nav>
<main class="main">
    <h1>Known Virus Pipeline Report</h1>
    <p style="color:var(--ink-secondary);font-size:13px">Generated: {now} &nbsp;|&nbsp; Directory: <span style="font-family:var(--font-mono);font-size:11.5px">{out_dir}</span></p>
    <div class="flow-strip" aria-label="Pipeline stages">
        <div class="flow-node"><b>S1</b><small>Detection</small></div><span class="flow-arrow">&#8594;</span>
        <div class="flow-node"><b>S2</b><small>Filter&nbsp;+&nbsp;Meta</small></div><span class="flow-arrow">&#8594;</span>
        <div class="flow-node"><b>S3</b><small>Variants</small></div><span class="flow-arrow">&#8594;</span>
        <div class="flow-node"><b>S4</b><small>Assembly</small></div><span class="flow-arrow">&#8594;</span>
        <div class="flow-node"><b>S5</b><small>Extract</small></div><span class="flow-arrow">&#8594;</span>
        <div class="flow-node"><b>S6</b><small>Post-hoc</small></div><span class="flow-arrow">&#8594;</span>
        <div class="flow-node"><b>S7</b><small>Selection</small></div><span class="flow-arrow">&#8594;</span>
        <div class="flow-node"><b>S8</b><small>Similarity</small></div><span class="flow-arrow">&#8594;</span>
        <div class="flow-node"><b>S9</b><small>DVG</small></div>
    </div>

    {ai_html}

    <script id="stage-data" type="application/json">{_build_stage_data_json(viruses, variant_summary, assembly_summary, detection_summary, extract_stats, overview_metrics)}</script>

    <button id="ai-btn" onclick="toggleAIPanel()" title="AI Settings">AI</button>
    <div id="ai-panel" style="display:none">
        <div style="display:flex;justify-content:space-between;align-items:center;padding:10px 14px;border-bottom:1px solid var(--border)">
            <span style="font-weight:700;font-size:14px;color:#5a4fcf">AI Settings</span>
            <button onclick="toggleAIPanel()" style="background:none;border:none;font-size:18px;cursor:pointer;color:var(--ink-secondary)">&times;</button>
        </div>
        <div style="padding:10px 14px">
            <label style="font-size:11px;color:var(--ink-secondary)">API Key (saved in this browser)</label>
            <input id="ai-api-key" type="password" placeholder="sk-..." style="width:100%;padding:6px 8px;margin:4px 0 8px;border:1px solid var(--border);border-radius:4px;font-size:12px" autocomplete="off">
            <label style="font-size:11px;color:var(--ink-secondary)">Provider</label>
            <select id="ai-provider" onchange="onProviderChange()" style="width:100%;padding:6px 8px;margin:4px 0 8px;border:1px solid var(--border);border-radius:4px;font-size:12px">
                <option value="openai">OpenAI</option>
                <option value="deepseek">DeepSeek</option>
                <option value="moonshot">Kimi (Moonshot)</option>
                <option value="ollama">Ollama (localhost)</option>
                <option value="custom">Custom</option>
            </select>
            <label style="font-size:11px;color:var(--ink-secondary)">Model</label>
            <input id="ai-model" value="gpt-4o-mini" style="width:100%;padding:6px 8px;margin:4px 0 8px;border:1px solid var(--border);border-radius:4px;font-size:12px">
        </div>
    </div>

    <section id="overview">
        <h2>Global Overview</h2>
        <div class="card-row">
            <div class="card"><div class="value">{overview_metrics.get("n_samples","?")}</div><div class="label">Samples</div></div>
            <div class="card"><div class="value">{len(viruses)}</div><div class="label">Virus Species</div></div>
            <div class="card"><div class="value">{overview_metrics.get("meta_records",0)}</div><div class="label">Meta-linked Records</div></div>
            <div class="card"><div class="value">{overview_metrics.get("raw_detections","?")}</div><div class="label">Raw Detections</div></div>
            <div class="card"><div class="value">{overview_metrics.get("best_detections","?")}</div><div class="label">High-Confidence</div></div>
            <div class="card"><div class="value">{overview_metrics.get("filtered_records","?")}</div><div class="label">Post-Filter Records</div></div>
            <div class="card"><div class="value">{overview_metrics.get("assemblies","?")}</div><div class="label">Assemblies</div></div>
            <div class="card"><div class="value">{overview_metrics.get("dvg_viruses","?")}</div><div class="label">Viruses with DVG</div></div>
            <div class="card"><div class="value">{len(stage_cfg)}</div><div class="label">Pipeline Stages</div></div>
        </div>
    </section>

    {funnel_section}

    {stage_sections}

    <section id="virus-summary">
        <h2>Virus Summary</h2>
        {_build_virus_summary_table(viruses)}
        {_build_methods_section(viruses, overview_metrics, total_records)}
    </section>

    <div class="footer">Generated by known_virus_pipeline — {now}</div>
</main>
<div class="lightbox" id="lightbox" onclick="if(event.target===this)closeLightbox()">
    <button class="lightbox-close" onclick="closeLightbox()" aria-label="Close">&times;</button>
    <!-- impeccable-disable-next-line broken-image -- src is assigned by openLightbox JS -->
    <img id="lightbox-img" alt="">
    <div class="lightbox-caption" id="lightbox-caption"></div>
    <button class="lightbox-dl" id="lightbox-dl" onclick="downloadLightbox()">Download</button>
</div>
<script>
// ── Lightbox ──
var lbImg=null,lbDl=null,lbCap=null;
document.addEventListener('DOMContentLoaded',function(){{
    var lb=document.getElementById('lightbox');
    lbImg=document.getElementById('lightbox-img');
    lbDl=document.getElementById('lightbox-dl');
    lbCap=document.getElementById('lightbox-caption');
    // Dynamic grid columns: 1→1col, 2→2col, 3→3col, 4+→4col max
    document.querySelectorAll('.chart-gallery').forEach(function(gallery){{
        var n=gallery.querySelectorAll('.chart-card').length;
        var cols;
        if(n<=1) cols=1;
        else if(n===2) cols=2;
        else if(n===3) cols=3;
        else cols=4;
        gallery.style.gridTemplateColumns='repeat('+cols+', 1fr)';
    }});
    // Add click handlers to all chart images
    document.querySelectorAll('.chart-card img').forEach(function(img){{
        img.style.cursor='pointer';
        img.addEventListener('click',function(e){{
            e.stopPropagation();
            lbImg.src=this.src;
            lbCap.textContent=this.alt||'';
            lbDl.setAttribute('data-src',this.src);
            lbDl.setAttribute('data-name',(this.alt||'chart').replace(/[^a-zA-Z0-9_.-]/g,'_')+'.png');
            lb.classList.add('show');
            document.body.style.overflow='hidden';
        }});
    }});
    // Add download buttons to chart cards with images
    document.querySelectorAll('.chart-card').forEach(function(card){{
        var img=card.querySelector('img');
        if(!img)return;
        var btn=document.createElement('button');
        btn.className='chart-dl-btn';
        btn.textContent='Save';
        btn.setAttribute('aria-label','Download this chart as PNG');
        btn.title='Download image';
        btn.addEventListener('click',function(e){{
            e.stopPropagation();
            downloadImg(img.src,(img.alt||'chart').replace(/[^a-zA-Z0-9_.-]/g,'_')+'.png');
        }});
        card.appendChild(btn);
    }});
}});
function closeLightbox(){{
    var lb=document.getElementById('lightbox');
    lb.classList.remove('show');
    document.body.style.overflow='';
}}
function downloadLightbox(){{
    var dl=document.getElementById('lightbox-dl');
    downloadImg(dl.getAttribute('data-src'),dl.getAttribute('data-name'));
}}
function downloadImg(src,name){{
    var a=document.createElement('a');
    a.href=src;a.download=name;
    document.body.appendChild(a);a.click();
    document.body.removeChild(a);
}}
// ── Chat Q&A ──
function toggleChat(){{
    var p=document.getElementById('chat-panel');
    p.style.display=p.style.display==='none'?'flex':'none';
}}
function sendChat(){{
    var inp=document.getElementById('chat-input');
    var q=inp.value.trim();if(!q)return;
    var key=document.getElementById('ai-api-key').value||sessionStorage.getItem('ai_key');
    if(!key){{toggleAIPanel();alert('Please enter your API key first');return;}}
    if(!sessionStorage.getItem('ai_key')){{sessionStorage.setItem('ai_key',key);sessionStorage.setItem('ai_provider',document.getElementById('ai-provider').value);sessionStorage.setItem('ai_model',document.getElementById('ai-model').value);}}
    var msgs=document.getElementById('chat-msgs');
    msgs.innerHTML+='<div class=\"chat-msg user\">'+q.replace(/</g,'&lt;')+'</div>';
    inp.value='';msgs.scrollTop=msgs.scrollHeight;
    var thinking=document.createElement('div');thinking.className='chat-msg ai';thinking.textContent='...';
    msgs.appendChild(thinking);
    var url=getAIURL();
    var ctx=JSON.stringify(stageData);
    fetch(url,{{method:'POST',headers:{{'Content-Type':'application/json','Authorization':'Bearer '+key}},
        body:JSON.stringify({{model:document.getElementById('ai-model').value,
        messages:[{{role:'system',content:'You are a virology research assistant. Answer questions about this pipeline report concisely in English. Use these data: '+ctx}},{{role:'user',content:q}}],
        temperature:0.3,max_tokens:300}})
    }}).then(r=>r.json()).then(d=>{{
        thinking.textContent=d.choices[0].message.content.trim();
        msgs.scrollTop=msgs.scrollHeight;
    }}).catch(e=>{{thinking.textContent='Error: '+e;}});
}}
// Close lightbox on Escape
document.addEventListener('keydown',function(e){{if(e.key==='Escape')closeLightbox();}});
// ── In-page AI ──
var sysMsg='You are a research scientist in virology. When writing academic summaries, strictly follow: 1. Base descriptions ONLY on the provided data and charts; never fabricate any numbers or findings. 2. When citing specific values, use objective wording such as the data show, detected, observed. 3. When citing a chart, use wording such as shown in Figure X, Table Y shows. 4. Use academic hedging such as suggest, imply, possibly when uncertain. 5. When the data do not support it, state clearly that nothing was detected or no significant difference was found. 6. Output in English, academic written style.';
var stageData=JSON.parse(document.getElementById('stage-data').textContent);
var saved=sessionStorage.getItem('ai_key');
if(saved)document.getElementById('ai-api-key').value=saved;
document.getElementById('ai-provider').value=sessionStorage.getItem('ai_provider')||'openai';
document.getElementById('ai-model').value=sessionStorage.getItem('ai_model')||'gpt-4o-mini';
function onProviderChange(){{
    var p=document.getElementById('ai-provider').value;
    var models={{openai:'gpt-4o-mini',deepseek:'deepseek-v4-flash',moonshot:'kimi-k2.6',ollama:'qwen2.5:7b',custom:''}};
    document.getElementById('ai-model').value=models[p]||'';
}}
function toggleAIPanel(){{
    var p=document.getElementById('ai-panel');
    p.style.display=p.style.display==='none'?'block':'none';
}}
function getAIURL(){{
    var p=document.getElementById('ai-provider').value;
    if(p==='openai')return'https://api.openai.com/v1/chat/completions';
    if(p==='deepseek')return'https://api.deepseek.com/v1/chat/completions';
    if(p==='moonshot')return'https://api.moonshot.cn/v1/chat/completions';
    if(p==='ollama')return'http://localhost:11434/v1/chat/completions';
    return prompt('Enter API base URL:')||'';
}}
function runStageAI(sn){{
    var key=document.getElementById('ai-api-key').value||sessionStorage.getItem('ai_key');
    if(!key){{toggleAIPanel();alert('Please enter your API key first');return;}}
    sessionStorage.setItem('ai_key',key);
    sessionStorage.setItem('ai_provider',document.getElementById('ai-provider').value);
    sessionStorage.setItem('ai_model',document.getElementById('ai-model').value);
    var btn=document.querySelector('#stage-'+sn+' .ai-stage-btn');
    var target=document.querySelector('#stage-'+sn+' .ai-stage-summary');
    if(!target){{
        target=document.createElement('p');target.className='ai-stage-summary';
        var h2=document.querySelector('#stage-'+sn+' h2');
        h2.parentNode.insertBefore(target,h2.nextSibling);
    }}
    target.style.display='block';target.textContent='Generating...';
    btn.disabled=true;
    var sd=stageData;
    // Build compact per-virus summary for prompts
    var vd=sd._virus_detail||[];var vdBrief=[];
    for(var i=0;i<vd.length;i++){{var v=vd[i];vdBrief.push(v.species+' ('+v.accession+'): '+v.samples+'samples, cov='+v.coverage+'%, depth='+v.depth+'x, CPM='+v.cpm+(v.pi!=null?', pi='+v.pi:'')+(v.assembly?', '+v.assembly:''));}}
    var vdStr=vdBrief.join('\\n');
    var prompts={{
        1:'Stage 1 Salmon pseudo-alignment + Poisson-filtered viral quantification.\\nPer virus: '+vdStr+'\\nEmbedded charts: freq_multi_metrics_log10, sample_virus_count_bar, virus_occurrence_bar, cooccurrence/spearman heatmap, Co-infection_Matrix.\\n\\nWrite a 2-3 sentence academic summary in English, briefly reviewing for each virus the number of detected samples, abundance and coverage. Cite specific numbers. Do not fabricate.',
        2:'Stage 2 Multi-dimensional threshold filtering (coverage >=50%, depth >=5x, reads >=100).\\n'+sd.s2.pre+' records -> '+sd.s2.post+' retained.\\nPer-virus detail: '+vdStr+'\\nEmbedded tables: Per-Sample/Per-Virus Filtering.\\n\\nBriefly review, for each virus, the post-filter coverage and depth. Cite specific numbers. Do not fabricate.',
        3:'Stage 3 FreeBayes+SnpEff+SNPGenie.\\nPer-virus data: '+JSON.stringify(sd.s3.per_virus)+'\\n\\nFor each virus, review its nucleotide diversity (π). High π = rich quasispecies; low π = conserved. Cite specific values. Do not fabricate.',
        4:'Stage 4 12-step de novo assembly.\\nPer-virus data: '+JSON.stringify(sd.s4.per_virus)+'\\n\\nFor each virus, review assembly completeness (N%, N50). Cite specific values. Do not fabricate.',
        5:'Stage 5 Extraction + N gap filling.\\nPer-virus data: '+JSON.stringify(sd.s5.per_virus)+'\\n\\nReview contig completeness (N%=0 means complete). Cite specific values. Do not fabricate.',
        6:'Stage 6 Multi-dimensional variant visualization chart set, grouped by analysis type.\\nChart content (per virus): Variant Landscape (genome-wide mutation distribution), Ts/Tv Pie (transition/transversion ratio pie -- a ratio >1 suggests random mutation dominance), Functional Pie (variant functional classification), Clustermap (cross-sample variant clustering), AFS (allele frequency spectrum -- low-frequency dominant suggests rare variants, high-frequency dominant suggests fixed differences), dN/dS scatter plots and 3D PCA (three-dimensional clustering of population-genetic parameters).\\n'+sd.viruses+' virus species are included in the analysis (single-sample viruses skipped).\\n\\nBased on the chart types above, write a 2-3 sentence academic summary. Chart names may be cited. Do not fabricate patterns that were not observed.',
        7:'Stage 7 Codon-level natural selection testing with HyPhy (FEL/MEME/BUSTED/PRIME).\\nData: '+sd.s7.coding_viruses+'/'+sd.viruses+' virus species have coding sequences; the rest are non-coding RNA (viroids/satellite RNA, which do not encode proteins and therefore are not amenable to dN/dS analysis).\\nEmbedded table: per-gene BUSTED selection results (gene name, ω3 value, p-value, positive selection sites, negative selection sites, total sites).\\n\\nBased on the data above, write a 2-3 sentence academic summary. Note: ω3≈1 with p>0.05 indicates no significant positive selection signal; a prevalence of negative selection sites suggests purifying selection dominates.',
        8:'Stage 8 SDT-style pairwise NT identity analysis.\\nInter-sample sequences of '+sd.viruses+' virus species were compared.\\nEmbedded charts: NT Identity Heatmap (pairwise identity heatmap + hierarchical clustering tree), Identity Distribution (identity histogram), pairwise statistics (mean/range/number of pairs).\\n\\nBased on the data above, write a 2-3 sentence academic summary. Identity >95% usually indicates the same strain; 85-95% is intra-species variation. State only the observed identity range.',
        9:'Stage 9 DVG and recombination detection with ViReMa.\\nData: '+sd.s9.viruses_with_dvg+'/'+sd.viruses+' virus species have DVG events.\\nEmbedded charts: Circos 4-track (mutation density/gene annotation/deletion coverage/duplication coverage four-track circular plot), Arc Diagram (donor-acceptor junction network arcs), Acceptor NT Frequency (acceptor breakpoint nucleotide frequency), Top Events (most frequent recombination events).\\nEmbedded table: DVG Events per Virus (DVG reads, jump events, samples involved per virus).\\n\\nBased on the data above, write a 2-3 sentence academic summary. Cite charts and tables. Discuss the DIP hypothesis only when the data support it.'
    }};
    var prompt=prompts[sn]||('Pipeline overview: '+sd.samples+' samples, '+sd.viruses+' virus species, 10-stage analysis. Write 2-3 sentence academic summary based ONLY on provided data.');
    prompt+='\\n\\nOutput in English. Academic written style. Cite chart names. Do not fabricate data. Use cautious wording where uncertain.'
    var url=getAIURL();
    fetch(url,{{method:'POST',headers:{{'Content-Type':'application/json','Authorization':'Bearer '+key}},
        body:JSON.stringify({{model:document.getElementById('ai-model').value,
        messages:[{{role:'system',content:sysMsg}},{{role:'user',content:prompt}}],
        temperature:0.3,max_tokens:400}})
    }}).then(r=>r.json()).then(d=>{{
        target.textContent=d.choices[0].message.content.trim();
        btn.disabled=false;
    }}).catch(e=>{{target.textContent='Error: '+e;btn.disabled=false;}});
    // Add copy button after result display
    var cp=target.nextElementSibling;if(!cp||!cp.classList.contains('ai-copy-btn')){{
        cp=document.createElement('button');cp.className='ai-copy-btn';cp.textContent='Copy';cp.title='Copy to clipboard';
        cp.addEventListener('click',function(){{navigator.clipboard.writeText(target.textContent);this.textContent='Copied!';setTimeout(()=>this.textContent='Copy',1500);}});
        target.parentNode.insertBefore(cp,target.nextSibling);
    }}
}}
function runVirusAI(sn,acc){{
    var key=document.getElementById('ai-api-key').value||sessionStorage.getItem('ai_key');
    if(!key){{toggleAIPanel();alert('Please enter your API key first');return;}}
    sessionStorage.setItem('ai_key',key);sessionStorage.setItem('ai_provider',document.getElementById('ai-provider').value);sessionStorage.setItem('ai_model',document.getElementById('ai-model').value);
    var vsa='s'+sn+'-virus-'+acc;
    // Escape dots in accession for CSS selector
    var sel='[id=\"'+vsa+'\"]';
    var target=document.querySelector(sel+' .ai-stage-summary');
    if(!target){{
        target=document.createElement('p');target.className='ai-stage-summary';
        var h=document.querySelector(sel+' h3');if(h)h.parentNode.insertBefore(target,h.nextSibling);
    }}
    target.style.display='block';target.textContent='Generating...';
    var sd=stageData;var vd=sd._virus_detail||[];
    var vinfo='';for(var i=0;i<vd.length;i++){{if(vd[i].accession===acc){{vinfo=JSON.stringify(vd[i]);break;}}}}
    if(!vinfo)vinfo=acc+' (data not found)';
    var stageFocus={{3:'This virus is in the Stage 3 Variants section. Embedded chart: per-virus variant metrics table (Samples/Cov/Depth/Pi/Shannon).\\n\\nReview the nucleotide diversity (π) and Shannon index of this virus. High π = rich quasispecies diversity (actively replicating quasispecies population); low π = highly conserved (strong purifying selection or a recent bottleneck). Assess variant-calling reliability in light of coverage and depth. State only facts supported by the data.',4:'This virus is in the Stage 4 Assembly section. Embedded charts: per-sample assembly evolution table (DeNovo->Ultimate length/N50/N) and collapsible stage details (length change per step).\\n\\nReview the assembly completeness of this virus: number of assemblies, best length, N content (N>0 = unclosed gaps), and percent of the reference genome (>100% suggests multimer/tandem-repeat assembly, 100% = complete monomer, <100% = partial assembly). Use the stage evolution to discuss how assembly improved.',5:'This virus is in the Stage 5 Extract section. Embedded table: extraction statistics (Numbers/Avg/Min/Max bp/N%).\\n\\nReview the extraction results of this virus: contig count, mean/min/max length, N% (0 = complete, no gaps). N%>0 indicates unfilled regions. Compare with the reference genome to assess extraction completeness.',6:'This virus is in the Stage 6 Post-hoc section. Embedded charts (grouped): Variant Landscape (genome-wide mutation distribution), Ts/Tv Pie (transition/transversion ratio >1 = random mutation dominance), Functional Pie (functional classification), Clustermap (cross-sample variant lineage), AFS (allele frequency spectrum), dN/dS scatter plots and 3D PCA (population-genetic parameter clustering).\\n\\nReview with the charts: the distribution of mutations along the genome (hotspots or not), the evolutionary driver reflected by the Ts/Tv ratio, the biological meaning of dN/dS departing from neutrality (purifying vs positive selection vs neutral evolution), and whether PCA clustering reveals distinct evolutionary lineages.',7:'This virus is in the Stage 7 Capheine section. Embedded table: per-gene BUSTED selection results (gene/ω3/p-value/positive sites/negative sites/total sites).\\n\\nReview the selection pressure on the genes of this virus: ω3>1 suggests the gene may be under positive selection (requires p<0.05), while ω3≈1 with p>0.05 indicates no significant positive selection signal. A prevalence of negative selection sites suggests purifying selection dominates (strong functional constraint). For non-coding RNA viruses (viroids/satellite RNA), note that dN/dS analysis is not applicable -- their evolutionary constraint comes from RNA secondary structure rather than protein coding.',8:'This virus is in the Stage 8 Similarity section. Embedded charts: NT Identity Heatmap (pairwise identity heatmap + hierarchical clustering tree), Identity Distribution (identity histogram), pairwise statistics (mean/range/number of pairs).\\n\\nReview the inter-sample NT identity of this virus: >95% usually indicates the same strain (recent transmission), 85-95% is intra-species variation, <85% is a different species. Discuss the evolutionary relationships revealed by the clustering tree. High identity suggests a recent common ancestor or a single introduction event; low identity suggests long-term local evolution or multiple introductions.',9:'This virus is in the Stage 9 DVG section. Embedded charts: Circos 4-track (mutation density/gene annotation/deletion coverage/duplication coverage), Arc Diagram (donor-acceptor junction network), Acceptor NT Frequency (acceptor breakpoint nucleotide preference), Top Events (most frequent recombination events). Embedded table: DVG Events per Virus (DVG reads/jump events/samples involved).\\n\\nReview the DVG/recombination events of this virus: DVG reads and jump events reflect recombination activity. Biological significance of DVGs (defective viral genomes): defective interfering particles (DIPs) can competitively suppress standard virus replication and affect infection progression and pathogenicity. The recombination type (deletion/duplication) relates to the error-prone template switching of the viral RNA polymerase.'}};
    var focus=stageFocus[sn]||'Review the key metrics of this virus such as detection, coverage and abundance.';
    var prompt='You are a virology researcher. Based on the following virus data, write a 1-2 sentence academic summary in English. '+focus+'\\n\\nVirus data: '+vinfo+'\\n\\nState only facts supported by the data. Cite specific numbers. Do not fabricate. In English.';
    fetch(getAIURL(),{{method:'POST',headers:{{'Content-Type':'application/json','Authorization':'Bearer '+key}},
        body:JSON.stringify({{model:document.getElementById('ai-model').value,
        messages:[{{role:'system',content:sysMsg}},{{role:'user',content:prompt}}],
        temperature:0.3,max_tokens:300}})
    }}).then(r=>r.json()).then(d=>{{
        target.textContent=d.choices[0].message.content.trim();
        var cp=target.nextElementSibling;if(!cp||!cp.classList.contains('ai-copy-btn')){{
            cp=document.createElement('button');cp.className='ai-copy-btn';cp.textContent='Copy';cp.title='Copy to clipboard';
            cp.addEventListener('click',function(){{navigator.clipboard.writeText(target.textContent);this.textContent='Copied!';setTimeout(()=>this.textContent='Copy',1500);}});
            target.parentNode.insertBefore(cp,target.nextSibling);
        }}
    }}).catch(e=>{{target.textContent='Error: '+e;}});
}}
</script>
{_TW}</body></html>"""

    with open(out_html, "w", encoding="utf-8") as f:
        f.write(html)
    return str(out_html)


# Stage 2 过滤阈值: 与主流程 auto_known_virus.py 默认值 (-c 50 / -d 5 / -r 100,
# 即 STAGE_METHODS 记载的 Rep_Coverage≥50%, Rep_MeanDepth≥5x, reads≥100) 保持一致
_FILTER_MIN_COV, _FILTER_MIN_DEPTH, _FILTER_MIN_READS = 50.0, 5.0, 100.0


def ensure_summary_plots(out_dir):
    import subprocess as _sp
    out = Path(out_dir)
    script_dir = Path(__file__).resolve().parent
    filter_dir = out / _D['a_filter']
    fp_plot = filter_dir / "filter_summary_plot.pdf"
    if not fp_plot.exists():
        high_conf_tsv = filter_dir / "high_conf.summary.tsv"
        best_tsv = None
        for cand in (out / _D['a_detect'] / "summary" / "all_viruses.best.summary.tsv",
                     filter_dir / "all_viruses.best.summary.tsv"):
            if cand.is_file():
                best_tsv = cand
                break
        # 仅当输入确为 high_conf 层时才允许回写 high_conf.summary.tsv;
        # 用 best 层重跑时只写预览 TSV, 不得凭空生成"高置信"层 (阈值用主流程同款)。
        if high_conf_tsv.is_file():
            src, dst = high_conf_tsv, high_conf_tsv
        elif best_tsv is not None:
            src, dst = best_tsv, filter_dir / "filter_summary_preview.tsv"
        else:
            src = dst = None
        if src is not None:
            fs = script_dir / "utils" / "filter_summary.py"
            if fs.is_file():
                # [P1] 正式产物保护: 在临时目录跑过滤画图, 只搬回 plot 文件;
                # 不再用报告侧常量阈值原地改写正式 high_conf.summary.tsv (阈值不同会污染正式层)
                import tempfile as _tf, shutil as _shutil
                with _tf.TemporaryDirectory() as _td:
                    _tdp = Path(_td)
                    _tmp_in = _tdp / src.name
                    _shutil.copy2(src, _tmp_in)
                    _sp.run([sys.executable, str(fs), "-i", str(_tmp_in), "-o", str(_tmp_in),
                             "-c", str(_FILTER_MIN_COV), "-d", str(_FILTER_MIN_DEPTH),
                             "-r", str(_FILTER_MIN_READS), "--plot"], capture_output=True)
                    for _pn in ("filter_summary_plot.pdf", "filter_summary_plot.png"):
                        _pp = _tdp / _pn
                        if _pp.is_file():
                            _shutil.move(str(_pp), str(filter_dir / _pn))
    asm_plot = _resolve_dir(out, _D['a_extract']) / "assembly_stats.pdf"
    if not asm_plot.exists():
        asm_dir = _resolve_dir(out_dir, _D['a_assembly'])
        if asm_dir.is_dir():
            es = script_dir / "utils" / "extract_full_fasta.py"
            if es.is_file():
                # [P1] 正式产物保护: extract 画图跑在临时目录, 只搬回 plot/stats,
                # 不再用 --max_n 100 --min_len 100 (无 --fill) 覆盖正式 .full.fasta 产物
                import tempfile as _tf, shutil as _shutil
                with _tf.TemporaryDirectory() as _td:
                    _tdp = Path(_td)
                    _sp.run([sys.executable, str(es), "-d", str(asm_dir), "-o", str(_tdp),
                             "--plot", "--max_n", "100", "--min_len", "100"], capture_output=True)
                    _dst_extract = out / _D['a_extract']
                    for _pn in ("assembly_stats.pdf", "assembly_stats.png", "assembly_stats.csv"):
                        _pp = _tdp / _pn
                        if _pp.is_file():
                            _shutil.move(str(_pp), str(_dst_extract / _pn))


def copy_results_to_report(out_dir, report_dir):
    """Copy key result files from each stage into the report directory.

    Copies summary tables, key plots, and per-virus charts to the report
    stage dir (layout-aware: 09_Report in standard, 10_report in legacy)
    so users can easily share the complete report package.
    """
    out = Path(out_dir); rep = Path(report_dir)
    rep.mkdir(parents=True, exist_ok=True)
    copied = []; skipped = []

    def _cp(src_rel, dst_rel=None):
        src = Path(src_rel) if Path(src_rel).is_absolute() else out / src_rel
        if not src.exists(): skipped.append(str(src_rel)); return
        if dst_rel:
            dst = rep / dst_rel
        else:
            dst = rep / src.relative_to(out) if str(src).startswith(str(out)) else rep / src.name
        if src.is_file():
            dst.parent.mkdir(parents=True, exist_ok=True)
            if not dst.exists() or src.stat().st_mtime > dst.stat().st_mtime:
                shutil.copy2(src, dst); copied.append(str(src_rel))
        elif src.is_dir():
            for f in src.rglob("*"):
                if f.is_file():
                    d = dst / f.relative_to(src)
                    d.parent.mkdir(parents=True, exist_ok=True)
                    if not d.exists() or f.stat().st_mtime > d.stat().st_mtime:
                        shutil.copy2(f, d); copied.append(str(f.relative_to(out)))

    # Stage 1: Detection summary
    _cp(os.path.join(_D['a_detect'], 'summary/all_viruses.best.summary.tsv'), "S1_Detection/all_viruses.best.summary.tsv")

    # Stage 2: Filter results
    _cp(os.path.join(_D['a_filter'], 'high_conf.summary.tsv'), "S2_Filter/high_conf.summary.tsv")
    _cp(os.path.join(_D['a_filter'], 'filter_summary_plot.png'), "S2_Filter/filter_summary_plot.png")
    _cp(os.path.join(_D['a_filter'], 'filter_summary_plot.pdf'), "S2_Filter/filter_summary_plot.pdf")
    # Stage 2b: Metadata association results ([2b] meta block output;
    # 编排器写 04_Variants/metadata_association, 02_Filtering 为 legacy 回退)
    _ma = _resolve_dir(out, _D['a_post']) / 'metadata_association'
    if not _ma.is_dir():
        _ma = out / _D['a_filter'] / 'metadata_association'
    _cp(str(_ma), "S2_Filter/metadata_association")

    # Stage 3: Variants summary
    _cp(os.path.join(_D['a_align'], 'summary/all_summary.tsv'), "S3_Variants/all_summary.tsv")
    _cp(os.path.join(_D['a_align'], 'summary/Coinfection_Matrix_Reads.tsv'), "S3_Variants/Coinfection_Matrix_Reads.tsv")

    # Stage 5: Assembly stats
    _cp(os.path.join(_D['a_extract'], 'assembly_stats.png'), "S5_Assembly_Stats/assembly_stats.png")
    _cp(os.path.join(_D['a_extract'], 'assembly_stats.pdf'), "S5_Assembly_Stats/assembly_stats.pdf")

    # Stage 9b: Stage funnel (五阶段漏斗表+图, 由 _build_funnel_section 生成)
    _cp(os.path.join(_D['a_report'], 'stage_funnel'), "S9_StageFunnel")

    # Stage 6: Post-hoc charts (per-virus)
    post_dir = _resolve_dir(out, _D['a_post'])
    if post_dir.is_dir():
        for vdir in post_dir.iterdir():
            if not vdir.is_dir() or vdir.name.startswith("run_"): continue
            for subd in vdir.iterdir():
                if subd.is_dir():
                    _cp(str(subd.relative_to(out)), f"S6_Post_hoc/{vdir.name}/{subd.name}")

    # Stage 7: Capheine (已迁 ⑤; 07_capheine 为 legacy 产物候选)
    cap_dir = out / "07_capheine"
    if cap_dir.is_dir():
        for vdir in cap_dir.iterdir():
            if vdir.is_dir() and not vdir.name.startswith("run_"):
                _cp(str(vdir.relative_to(out)), f"S7_Capheine/{vdir.name}")
    # Stage 7 (phylo): --phylo-root 提供时额外打包 ⑤ capheine 位点统计图
    if PHYLO_ROOT is not None:
        pr = Path(PHYLO_ROOT)
        if pr.is_dir():
            for f in sorted(pr.rglob("selection_per_gene.*")):
                if f.is_file():
                    try:
                        acc = f.parts[f.parts.index("gcva_analysis") - 1]
                    except ValueError:
                        acc = f.parent.parent.name
                    _cp(str(f), f"S7_Capheine/{acc}/{f.name}")

    # Stage 8: Similarity (recursive: heatmaps are nested deep); 布局感知目录
    sim_dir = locate_analysis_dir(out, "similarity")
    if sim_dir.is_dir():
        for vdir in sim_dir.iterdir():
            if vdir.is_dir():
                for f in vdir.rglob("*heatmap*"):
                    _cp(str(f.relative_to(out)), f"S8_Similarity/{vdir.name}/{f.name}")
                for f in vdir.rglob("*composite*"):
                    _cp(str(f.relative_to(out)), f"S8_Similarity/{vdir.name}/{f.name}")
                for f in vdir.rglob("*similarity*"):
                    _cp(str(f.relative_to(out)), f"S8_Similarity/{vdir.name}/{f.name}")

    # Stage 9: DVG plots
    dvg_plots = out / os.path.join(_D['a_dvg'], 'Summary_Analysis_Report/Virus_Specific_Plots')
    if dvg_plots.is_dir():
        for vdir in dvg_plots.iterdir():
            if vdir.is_dir():
                _cp(str(vdir.relative_to(out)), f"S9_DVG/{vdir.name}")

    return copied, skipped


def generate_ai_summary(viruses, variant_summary, assembly_summary, detection_summary, extract_stats, overview_metrics):
    """Call LLM for an AI-powered research brief + IMRaD summary."""
    import json as _json, re, urllib.request
    provider = getattr(generate_ai_summary, '_provider', 'openai')
    model = getattr(generate_ai_summary, '_model', 'gpt-4o-mini')
    api_key = getattr(generate_ai_summary, '_api_key', '')
    base_url = getattr(generate_ai_summary, '_base_url', '')
    if not api_key:
        return '<div class="ai-error">AI summary requires --ai-key</div>'

    # Build pipeline stats for prompt
    n_samples = overview_metrics.get("n_samples", "?")
    n_virus = len(viruses)
    raw_d = overview_metrics.get("raw_detections", "?")
    best_d = overview_metrics.get("best_detections", "?")
    filtered = overview_metrics.get("filtered_records", "?")
    n_asm = overview_metrics.get("assemblies", "?")
    dvg_n = overview_metrics.get("dvg_viruses", "?")

    virus_lines = []
    for acc, data in sorted(viruses.items()):
        sp = data.get("species", acc)[:40]
        n = data.get("n_samples", "?")
        cov = data.get("coverage", "?")
        vs = variant_summary.get(acc, {})
        asm = assembly_summary.get(acc, {})
        es = extract_stats.get(acc, {})
        line = f"  {sp} ({acc}): {n} samples, cov {cov}%"
        if vs: line += f", pi={vs.get('avg_pi','?'):.4f}" if vs.get('avg_pi') else ""
        if asm: line += f", asm={len(asm.get('samples',[]))} assembled"
        if es: line += f", {es.get('n_contigs','?')} contigs"
        virus_lines.append(line)

    prompt = f"""You are a virology research scientist. Write an English research summary for a Known Virus Analysis Pipeline report.

## Pipeline Data
- Samples: {n_samples}
- Virus species detected: {n_virus}
- Raw detections: {raw_d} → High-confidence: {best_d} → Post-filter: {filtered}
- Assemblies: {n_asm}
- Viruses with DVG/recombination: {dvg_n}
- 10-stage pipeline: Detection → Filter → Variants → Assembly → Extract → Post-hoc → Capheine → Similarity → DVG → Report

## Per-Virus Results
{chr(10).join(virus_lines)}

## Output Format
[BRIEF]
150-200 word English research highlight. Cover: study purpose, sample scale, key findings (viruses detected, assembly quality, selection evidence), main conclusions. Style: journal Highlights.

[DETAILED]
IMRaD structure (300-400 words, English):
[BACKGROUND] 1-2 sentences. Known virus surveillance context.
[METHODS] 3-4 sentences. Summarize the 10-stage pipeline.
[RESULTS] 4-5 sentences. Report quantitative findings with numbers.
[DISCUSSION] 2-3 sentences. Significance, limitations, next steps.

Requirements: Professional academic tone, cite exact numbers from data, English output, only output the format requested."""

    try:
        url = base_url or "https://api.openai.com/v1/chat/completions"
        body = _json.dumps({
            "model": model, "messages": [
                {"role": "system", "content": "You are a virology research scientist. Output in English, objective, precise, cite data."},
                {"role": "user", "content": prompt}
            ], "temperature": 0.3, "max_tokens": 2000
        }).encode()
        req = urllib.request.Request(url, data=body, headers={
            "Content-Type": "application/json", "Authorization": f"Bearer {api_key}"
        })
        with urllib.request.urlopen(req, timeout=300) as resp:
            result = _json.loads(resp.read())
        raw = result["choices"][0]["message"]["content"].strip()
        print(f"  AI Summary generated ({len(raw)} chars)")

        brief_m = re.search(r'\[BRIEF\]\s*(.*?)(?=\[DETAILED\]|\Z)', raw, re.DOTALL | re.IGNORECASE)
        brief = brief_m.group(1).strip() if brief_m else raw[:500]
        det_m = re.search(r'\[DETAILED\]\s*(.*)', raw, re.DOTALL | re.IGNORECASE)
        detailed = det_m.group(1).strip() if det_m else ""

        def _sec(text, tag):
            m = re.search(rf'\[{tag}\]\s*(.*?)(?=\[(?:BACKGROUND|METHODS|RESULTS|DISCUSSION)\]|\Z)', text, re.DOTALL | re.IGNORECASE)
            return m.group(1).strip() if m else ""

        secs = {}
        if detailed:
            secs = {k: _sec(detailed, k) for k in ["BACKGROUND","METHODS","RESULTS","DISCUSSION"]}
        if not any(secs.values()):
            secs = {"BACKGROUND": detailed or raw}

        labels = {"BACKGROUND": ("Background", "ai-bg"), "METHODS": ("Methods", "ai-methods"),
                   "RESULTS": ("Results", "ai-results"), "DISCUSSION": ("Discussion", "ai-discussion")}
        detail_items = ""
        for k in ("BACKGROUND","METHODS","RESULTS","DISCUSSION"):
            t = secs.get(k, "")
            if not t: continue
            lb, cl = labels[k]
            detail_items += f'<div class="ai-block {cl}"><div class="ai-block-title">{lb}</div><div class="ai-block-text">{_esc(t)}</div></div>'

        return f'''<div class="ai-section" id="ai-summary">
<div class="ai-header"><span class="ai-title">AI Research Brief</span><span class="ai-model-tag">{_esc(model)}</span></div>
<div class="ai-brief">{_esc(brief)}</div>
<details class="ai-detail-toggle"><summary>IMRaD Detailed Report</summary><div class="ai-detail-content">{detail_items}</div></details>
</div>'''
    except Exception as e:
        print(f"  [WARN] AI Summary failed: {e}")
        return f'<div class="ai-error">AI summary failed: {e}</div>'


def generate_stage_ai_summaries(viruses, variant_summary, assembly_summary, detection_summary, extract_stats, overview_metrics):
    """Generate a 1-2 sentence AI summary for each pipeline stage."""
    import json as _json, re, urllib.request
    provider = getattr(generate_ai_summary, '_provider', 'openai')
    model = getattr(generate_ai_summary, '_model', 'gpt-4o-mini')
    api_key = getattr(generate_ai_summary, '_api_key', '')
    base_url = getattr(generate_ai_summary, '_base_url', '')
    if not api_key: return {}

    n_samples = overview_metrics.get("n_samples", "?")
    n_virus = len(viruses)
    ds = detection_summary
    dv = ds.get("viruses", {}) if ds else {}
    raw_d = overview_metrics.get("raw_detections", "?")
    best_d = overview_metrics.get("best_detections", "?")
    filtered = overview_metrics.get("filtered_records", "?")
    n_asm = overview_metrics.get("assemblies", "?")

    # Stage-specific data
    s1_info = f"Raw={raw_d}, High-conf={best_d}, {n_virus} species in {n_samples} samples"
    s2_info = f"Pre={raw_d}→Post={filtered}, retention={filtered}/{raw_d} records"
    s3_lines = []
    for acc, v in sorted(viruses.items()):
        vs = variant_summary.get(acc, {})
        if vs:
            s3_lines.append(f"{v.get('species',acc)[:25]}: {vs.get('samples','?')}samples, pi={vs.get('avg_pi','?'):.4f}" if vs.get('avg_pi') else f"{v.get('species',acc)[:25]}: {vs.get('samples','?')}samples")
    s3_info = "; ".join(s3_lines[:5])
    s4_lines = []
    for acc in sorted(viruses.keys()):
        asm = assembly_summary.get(acc, {})
        samples = asm.get("samples", [])
        valid = [s for s in samples if s.get("final_len")]
        if valid:
            best_len = max(s.get("final_len", 0) or 0 for s in valid)
            s4_lines.append(f"{acc}: {len(valid)}asm, best={best_len:.0f}bp")
    s4_info = "; ".join(s4_lines[:5])
    s5_lines = []
    for acc in sorted(viruses.keys()):
        es = extract_stats.get(acc, {})
        if es:
            s5_lines.append(f"{acc}: {es.get('n_contigs','?')}contigs, avg={es.get('avg_len','?'):.0f}bp")
    s5_info = "; ".join(s5_lines[:5])
    s7_coding = sum(1 for acc in viruses if variant_summary.get(acc, {}).get("ref_length", 0) > 1000)
    s7_info = f"{s7_coding}/{n_virus} viruses have coding sequences (capheine analysis)"
    s9_info = f"{overview_metrics.get('dvg_viruses','?')}/{n_virus} viruses with DVG events"

    prompt = f"""You are a virology research scientist. For each pipeline stage below, write ONE concise sentence (English) describing what was found. Be specific with numbers.

Pipeline: Known Virus Analysis (9 stages), {n_samples} samples, {n_virus} virus species.

Stage 1 Detection: {s1_info}
Stage 2 Filter: {s2_info}
Stage 3 Variants: {s3_info}
Stage 4 Assembly: {s4_info}
Stage 5 Extract: {s5_info}
Stage 6 Post-hoc: per-virus mutation landscapes, dN/dS, PCA (see charts)
Stage 7 Capheine: {s7_info}
Stage 8 Similarity: pairwise NT identity matrices
Stage 9 DVG: {s9_info}

Output format (EXACTLY):
S1: <one English sentence>
S2: <one English sentence>
...
S9: <one English sentence>

Requirements: One sentence each. Cite numbers. English. Professional tone. No extra text."""

    try:
        url = base_url or "https://api.openai.com/v1/chat/completions"
        body = _json.dumps({
            "model": model, "messages": [
                {"role": "system", "content": "You are a virology research scientist. Output ONLY in the requested format."},
                {"role": "user", "content": prompt}
            ], "temperature": 0.3, "max_tokens": 800
        }).encode()
        req = urllib.request.Request(url, data=body, headers={
            "Content-Type": "application/json", "Authorization": f"Bearer {api_key}"
        })
        with urllib.request.urlopen(req, timeout=300) as resp:
            raw = _json.loads(resp.read())["choices"][0]["message"]["content"].strip()
        print(f"  Per-stage AI summaries generated ({len(raw)} chars)")

        summaries = {}
        for m in re.finditer(r'S(\d):\s*(.+?)(?=\nS\d:|\Z)', raw, re.DOTALL):
            summaries[int(m.group(1))] = m.group(2).strip()
        return summaries
    except Exception as e:
        print(f"  [WARN] Stage AI summaries failed: {e}")
        return {}


def main():
    parser = argparse.ArgumentParser(description="Interactive HTML Pipeline Report Generator")
    parser.add_argument("-d", "--dir", required=True, help="Pipeline output root directory")
    parser.add_argument("-o", "--output", default=None, help="HTML report output path")
    parser.add_argument("--no_images", action="store_true", help="Skip image embedding (faster, smaller)")
    parser.add_argument("--phylo-root", default=None,
                        help="⑤ 系统发育管线输出根 (results_* 目录); 提供后 Stage 7 正选择节额外扫描 "
                             "capheine 产物 (selection_per_gene 位点图 + drhip combined_summary 汇总)")
    # --ai-api-key: 编排器 auto_known_virus.py 传此名; --ai-key 为历史名, 两者等价
    parser.add_argument("--ai-key", "--ai-api-key", dest="ai_key", default="",
                        help="OpenAI/DeepSeek API key for AI summary")
    parser.add_argument("--ai-provider", default="openai", choices=["openai","ollama","deepseek","moonshot","custom"])
    parser.add_argument("--ai-model", default="gpt-4o-mini", help="LLM model name")
    parser.add_argument("--ai-base-url", default="", help="Custom API base URL")

    args = parser.parse_args()
    # --no_images: 关闭所有图片字节内嵌 (卡片保留文件名占位)
    global _SKIP_IMAGES
    _SKIP_IMAGES = bool(args.no_images)
    # --phylo-root: ⑤ 管线输出根, Stage 7 正选择节据此扫描 capheine 产物
    global PHYLO_ROOT
    PHYLO_ROOT = Path(args.phylo_root) if getattr(args, "phylo_root", None) else None
    # Configure AI summary
    if args.ai_key:
        generate_ai_summary._provider = args.ai_provider
        generate_ai_summary._model = args.ai_model
        generate_ai_summary._api_key = args.ai_key
        if args.ai_base_url:
            generate_ai_summary._base_url = args.ai_base_url
        elif args.ai_provider == "ollama":
            generate_ai_summary._base_url = "http://localhost:11434/v1/chat/completions"
        elif args.ai_provider == "deepseek":
            generate_ai_summary._base_url = "https://api.deepseek.com/v1/chat/completions"
        elif args.ai_provider == "moonshot":
            generate_ai_summary._base_url = "https://api.moonshot.cn/v1/chat/completions"
        else:
            generate_ai_summary._base_url = ""
    out = Path(args.dir)
    if not out.exists():
        sys.exit(f"Directory not found: {args.dir}")

    # [P2-24i] legacy 存量项目: 分析根可能是 legacy 布局的 09a/09 两代分析目录
    # (standard 布局时分析根 = --dir 本身) → 三代候选根都查; 阶段目录名由 locate_analysis_dir 探测。
    # 2026-10-13: 候选名改由 io_layout 注册表解析 (原字面量会被 check_io_handoffs 判 FAIL)
    candidate_roots = [out]
    try:
        from mmpv_common.io_layout import build_discovery_dirs as _bdd
        candidate_roots += [d for d in _bdd(out, out, "legacy")['_analysis_dirs'] if d.is_dir()]
    except Exception:
        pass
    summary_in = None
    for p in ["high_conf.summary.tsv", "all_viruses.best.summary.tsv", "all_viruses.summary.tsv"]:
        for root in candidate_roots:
            for loc in [locate_analysis_dir(root, "filter") / p,
                        locate_analysis_dir(root, "detect") / "summary" / p]:
                if loc.is_file(): summary_in = loc; break
            if summary_in: break
        if summary_in: break
    if summary_in is None:
        sys.exit("No summary file found")

    ensure_summary_plots(args.dir)
    print("Collecting virus data...")
    viruses = collect_virus_data(args.dir, summary_in)
    print(f"  Found {len(viruses)} viruses")

    out_html = Path(args.output) if args.output else out / "Pipeline_Summary_Report.html"

    # Generate AI summary if API key provided
    ai_html = ""
    if getattr(generate_ai_summary, '_api_key', ''):
        print("Generating AI summary...")
        # Quick-load data needed for AI summary
        from pathlib import Path as _P
        vs = _collect_variant_summary(args.dir)
        asm = _collect_assembly_summary(args.dir)
        ds = _collect_detection_summary(args.dir)
        es = _collect_extract_stats(args.dir)
        # Build overview metrics (subset needed for AI)
        ov = {"n_samples": "?", "n_viruses": len(viruses)}
        raw_t = safe_read_csv(_P(args.dir) / os.path.join(_D['a_detect'], 'summary/all_viruses.raw.tsv'))
        best_t = safe_read_csv(_P(args.dir) / os.path.join(_D['a_detect'], 'summary/all_viruses.best.summary.tsv'))
        hc_t = safe_read_csv(_P(args.dir) / os.path.join(_D['a_filter'], 'high_conf.summary.tsv'))
        ov["raw_detections"] = len(raw_t) if raw_t is not None else "?"
        ov["best_detections"] = len(best_t) if best_t is not None else "?"
        ov["filtered_records"] = len(hc_t) if hc_t is not None else "?"
        n_as = 0; d_as = _P(args.dir) / _D['a_assembly']
        if d_as.is_dir():
            for vd in d_as.iterdir():
                if vd.is_dir() and not vd.name.startswith("run_"):
                    n_as += sum(1 for d in vd.iterdir() if d.is_dir())
        ov["assemblies"] = n_as
        d_dvg = _P(args.dir) / os.path.join(_D['a_dvg'], 'Summary_Analysis_Report/Virus_Specific_Plots')
        ov["dvg_viruses"] = sum(1 for d in d_dvg.iterdir() if d.is_dir()) if d_dvg.is_dir() else 0
        ai_html = generate_ai_summary(viruses, vs, asm, ds, es, ov)
        stage_sums = generate_stage_ai_summaries(viruses, vs, asm, ds, es, ov)
    else:
        stage_sums = {}

    print("Generating interactive HTML...")
    generate_html(args.dir, out_html, viruses, ai_html, stage_sums)
    print(f"  Report: {out_html} ({out_html.stat().st_size / 1024:.0f} KB)")

    # Copy result files to report directory
    report_dir = locate_analysis_dir(out, "report")
    print("Copying result files to report directory...")
    copied, skipped = copy_results_to_report(args.dir, report_dir)
    print(f"  Copied {len(copied)} files, skipped {len(skipped)} missing")
    for s in sorted(skipped):
        print(f"    (not found) {s}")


if __name__ == "__main__":
    main()
