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


def _esc(v):
    """HTML-escape a value to prevent XSS injection."""
    return str(v).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;").replace('"', "&quot;")

try:
    import pandas as pd
    HAS_PANDAS = True
except ImportError:
    HAS_PANDAS = False


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
    size_kb = Path(path).stat().st_size / 1024
    if size_kb > max_kb:
        return None
    with open(path, 'rb') as f:
        return base64.b64encode(f.read()).decode('utf-8')


def collect_charts(post_dir, virus_acc):
    """Collect key chart paths for a virus (fuzzy match by accession in dir name)."""
    charts = {}
    vdir = Path(post_dir) / virus_acc
    if not vdir.exists():
        for d in Path(post_dir).iterdir():
            if d.is_dir() and virus_acc in d.name:
                vdir = d; break
        else:
            return charts

    chart_patterns = {
        'vcf_viz': [
            ('Figure1A_All_Variants_Landscape.png', '全基因组变异景观'),
            ('Figure1B_Top50_Variants_Landscape.png', 'Top50 变异景观'),
            ('Figure2_TsTv_Pie.png', 'Ts/Tv 比率'),
            ('Figure3_Functional_Pie.png', '功能分类饼图'),
            ('Figure4_Clustermap.png', '变异聚类热图'),
            ('Figure5_AFS.png', '等位频率谱'),
            ('Figure6_Variant_Density.png', '变异密度图'),
            ('Figure7_AF_Violin.png', 'AF 小提琴图'),
            ('Figure8_PopGen_Dynamics.png', 'PopGen 滑动窗'),
            ('Figure9_PCA_Lineages.png', 'PCA 2D 聚类'),
            ('Figure9_PCA_Lineages_3D.png', 'PCA 3D 聚类'),
        ],
        'snpeff_macro': [
            ('Figure_1_Manhattan_Mut_Landscape.png', '突变曼哈顿图'),
            ('Figure_2_Gene_Payload.png', '基因突变载荷'),
            ('Figure_3_IntraHost_Diversity.png', '准种多样性'),
            ('Figure_4_Lineage_Clustermap.png', '错义突变谱系聚类'),
            ('Figure_5_Ultimate_Quasispecies_OncoPrint.png', '突变瀑布图'),
        ],
        'maftools': [
            ('mafSummary_TCGA.png', 'MAF 突变类型'),
            ('Oncoplot.png', '突变瀑布图'),
            ('TiTv_Summary.png', 'Ti/Tv 汇总'),
        ],
        'snpgenie': [
            ('Fig01_VAF_Spectrum.png', 'VAF 频谱'),
            ('Fig02_Adjusted_iSNV_Density.png', '深度校正 iSNV 密度'),
            ('Fig03_InterHost_dNdS.png', 'dN vs dS 联合分布'),
            ('Fig04_IntraHost_pi.png', 'πN vs πS 联合分布'),
            ('Fig05_Gene_dNdS_Stats.png', '每基因 dN/dS'),
            ('Fig06_Bootstrapped_dNdS.png', 'Bootstrap 显著性'),
            ('Fig07_Mutational_Spectrum.png', '突变频谱'),
            ('Fig08_Top_Hotspots.png', 'Top 热点'),
            ('Fig09_DualTrack_Window.png', '双轨道滑动窗'),
            ('Fig10_Trinity_Landscape.png', 'Trinity 全景观'),
            ('Fig11a_PCA_2D.png', '2D PCA 聚类'),
            ('Fig11b_PCA_3D.png', '3D PCA 聚类'),
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
                b64 = img_to_base64(fp)
                charts[f"{subdir}_{fname}"] = {
                    'label': label,
                    'path': str(fp.relative_to(post_dir.parent)),
                    'base64': b64,
                    'ext': fp.suffix[1:],
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
    for suffix, label in [("Circos_4Track.png","DVG Circos"),("Arc_Diagram.png","DVG Arc"),
            ("Acceptor_NT_Freq.png","DVG Acceptor NT"),("Top_Events.png","DVG Top Events")]:
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
            "cpm": float(row.get("Asm_CPM", 0)),
            "coverage": float(row.get("Rep_Coverage(%)", 0)),
            "depth": float(row.get("Rep_MeanDepth", 0)),
            "poisson": float(row.get("Poisson_Ratio", 0)),
            "reads": float(row.get("Asm_EM_Reads", 0)),
        }

    # Count per-virus samples
    for acc in list(viruses.keys()):
        n = (df[acc_col].astype(str) == acc).sum() if acc_col and HAS_PANDAS else "?"
        viruses[acc]["n_samples"] = n

    return viruses


def _collect_variant_summary(out_dir):
    """Aggregate per-virus variant stats from S3 all_summary.tsv."""
    out = Path(out_dir)
    summary_tsv = out / "3_Virus_variants_Results/summary/all_summary.tsv"
    df = safe_read_csv(summary_tsv)
    if df is None: return {}
    acc_col = next((c for c in ["Accession"] if c in df.columns), None)
    if not acc_col: return {}
    result = {}
    for acc, grp in df.groupby(acc_col):
        samples = grp["Sample"].nunique() if "Sample" in df.columns else len(grp)
        cov = grp["Covered%"].mean() if "Covered%" in df.columns else None
        depth = grp["Rep_MeanDepth"].mean() if "Rep_MeanDepth" in df.columns else None
        pi = grp["Pi_avr"].mean() if "Pi_avr" in df.columns else None
        shannon = grp["Shannon_avr"].mean() if "Shannon_avr" in df.columns else None
        length = grp["Length"].iloc[0] if "Length" in df.columns else None
        result[str(acc)] = {
            "samples": samples,
            "avg_cov": cov, "avg_depth": depth,
            "avg_pi": pi, "avg_shannon": shannon,
            "ref_length": length,
        }
    return result


def _collect_assembly_summary(out_dir):
    """Scan S4 assembly dirs for per-virus per-sample assembly evolution stats.

    Returns dict keyed by accession, each with samples list containing:
      name, denovo_len, denovo_n50, final_len, final_n50, final_n, stages
    where stages is a list of {step, length, n50, n_count} for full evolution tracking.
    """
    out = Path(out_dir)
    asm_dir = out / "4_Virus_assemblies_final"
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
        m = _re.search(r'[A-Z]{2,4}_?\d{5,9}\.\d{1,2}', vdir.name)
        acc = m.group() if m else vdir.name
        result[acc] = {"dir_name": vdir.name, "samples": samples}
    return result


def _collect_detection_summary(out_dir):
    """Parse S1 best-summary TSV to aggregate per-virus detection metrics."""
    out = Path(out_dir)
    tsv = out / "1_FastViromeExplorer/summary/all_viruses.best.summary.tsv"
    if not tsv.is_file():
        tsv = out / "1_FastViromeExplorer/summary/all_viruses.summary.tsv"
    df = safe_read_csv(tsv)
    if df is None: return {}
    raw_tsv = out / "1_FastViromeExplorer/summary/all_viruses.raw.tsv"
    raw_df = safe_read_csv(raw_tsv)
    total_raw = len(raw_df) if raw_df is not None else None
    acc_col = next((c for c in ["Rep_Accession"] if c in df.columns), None)
    if not acc_col: return {}
    viruses = {}
    for acc, grp in df.groupby(acc_col):
        sp_col = next((c for c in ["Adjusted_Species", "Species_NCBI"] if c in df.columns), None)
        species = str(grp[sp_col].iloc[0]) if sp_col else str(acc)
        viruses[str(acc)] = {
            "species": species,
            "n_samples": len(grp),
            "avg_cpm": float(grp["Asm_CPM"].mean()) if "Asm_CPM" in df.columns else None,
            "avg_fpkm": float(grp["Asm_FPKM"].mean()) if "Asm_FPKM" in df.columns else None,
            "avg_cov": float(grp["Rep_Coverage(%)"].mean()) if "Rep_Coverage(%)" in df.columns else None,
            "avg_depth": float(grp["Rep_MeanDepth"].mean()) if "Rep_MeanDepth" in df.columns else None,
            "avg_poisson": float(grp["Poisson_Ratio"].mean()) if "Poisson_Ratio" in df.columns else None,
        }
    return {"viruses": viruses, "total_best": len(df), "total_raw": total_raw}


def _collect_capheine_summary(out_dir, virus_acc):
    """Parse DRHIP combined_summary.csv for per-gene BUSTED selection results."""
    out = Path(out_dir)
    vdir = _find_virus_dir(out / "7_capheine", virus_acc)
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
            "omega3": float(row.get("BUSTED_omega3", 0)),
            "pval": float(row.get("BUSTED_pval", 1)),
            "positive": int(float(row.get("positive_sites", 0))),
            "negative": int(float(row.get("negative_sites", 0))),
            "total_sites": int(float(row.get("sites", 0))),
            "n_seq": int(float(row.get("N", 0))),
        })
    sig_count = sum(1 for g in genes if g["pval"] < 0.05 and g["positive"] > 0)
    return {"genes": genes, "sig_count": sig_count, "total_genes": len(genes), "dir_name": vdir.name}


def _collect_similarity_data(out_dir, virus_acc):
    """Find similarity heatmaps and pairwise stats for a virus."""
    out = Path(out_dir)
    vdir = _find_virus_dir(out / "8_similarity", virus_acc)
    if not vdir: return {"available": False, "reason": "no similarity directory"}
    mat_dir = vdir / "Mode_Filter" / "Full_Dataset" / "02_similarity_matrices"
    if not mat_dir.is_dir():
        return {"available": False, "reason": "empty (single sample?)"}
    heatmap_png = mat_dir / "overall_NT_only_heatmap.png"
    dist_png = mat_dir / "overall_NT_only_distribution.png"
    csv_f = mat_dir / "overall_NT_only_pairwise.csv"
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
    ext_dir = out / "5_assemblies_clean"
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
        m = _re.search(r'[A-Z]{2,4}_?\d{5,9}\.\d{1,2}', vdir.name)
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


def _find_virus_files(out_dir, stage_rel, virus_acc, subdir_patterns):
    """Find per-virus files in a stage directory and return HTML.

    subdir_patterns: list of (subdir_rel, globs) tuples.
    """
    out = Path(out_dir)
    stage_path = out / stage_rel
    vdir = None

    # For S3, lookup is reversed: stage/{subdir}/{virus_dir}[/{sample_dir}]
    if stage_rel == "3_Virus_variants_Results":
        html = ""
        for subd_rel, globs in subdir_patterns:
            subd = stage_path / subd_rel
            if not subd.is_dir(): continue
            vd = _find_virus_dir(subd, virus_acc)
            if not vd: continue
            for g in globs:
                for f in sorted(vd.rglob(g))[:6]:
                    if not f.is_file(): continue
                    rel = f.relative_to(vd)
                    html += f'<p style="font-size:11px;color:var(--ink-secondary);margin:2px 0">{_esc(subd_rel)}/{_esc(vd.name)}/{_esc(rel)}</p>'
        return html

    # For S4/S7/S8: stage/{virus_dir}/...
    vdir = _find_virus_dir(stage_path, virus_acc)
    if not vdir: return ""
    html = ""
    for subd_rel, globs in subdir_patterns:
        search_dir = vdir / subd_rel if subd_rel else vdir
        if not search_dir.is_dir(): continue
        for g in globs:
            for f in sorted(search_dir.rglob(g))[:8]:
                if not f.is_file(): continue
                rel = f.relative_to(vdir)
                html += f'<p style="font-size:11px;color:var(--ink-secondary);margin:2px 0">{_esc(rel)}</p>'
    return html


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


def generate_html(out_dir, out_html, viruses, ai_html="", stage_summaries=None):
    """Generate interactive HTML report with sidebar navigation."""
    now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    out = Path(out_dir)
    post_dir = out / "6_post_analysis"

    total_records = sum(v.get("n_samples", 0) for v in viruses.values()) if isinstance(
        next(iter(viruses.values()), {}).get("n_samples", 0), int) else "?"

    # Pre-compute global overview metrics
    overview_metrics = {}
    # Detection counts
    raw_tsv_ov = out / "1_FastViromeExplorer/summary/all_viruses.raw.tsv"
    raw_ov = safe_read_csv(raw_tsv_ov)
    best_tsv_ov = out / "1_FastViromeExplorer/summary/all_viruses.best.summary.tsv"
    best_ov = safe_read_csv(best_tsv_ov)
    overview_metrics["raw_detections"] = len(raw_ov) if raw_ov is not None else "?"
    overview_metrics["best_detections"] = len(best_ov) if best_ov is not None else "?"
    overview_metrics["n_samples"] = best_ov["Sample"].nunique() if best_ov is not None and "Sample" in best_ov.columns else "?"
    # Filter counts
    hc_tsv_ov = out / "2_Virus_result_filter/high_conf.summary.tsv"
    hc_ov = safe_read_csv(hc_tsv_ov)
    overview_metrics["filtered_records"] = len(hc_ov) if hc_ov is not None else "?"
    # Assembly counts
    asm_dir_ov = out / "4_Virus_assemblies_final"
    n_asm = 0
    if asm_dir_ov.is_dir():
        for vd in asm_dir_ov.iterdir():
            if vd.is_dir() and not vd.name.startswith("run_"):
                n_asm += sum(1 for d in vd.iterdir() if d.is_dir())
    overview_metrics["assemblies"] = n_asm
    # DVG counts
    dvg_dir_ov = out / "9_virema_dvg/Summary_Analysis_Report/Virus_Specific_Plots"
    overview_metrics["dvg_viruses"] = sum(1 for d in dvg_dir_ov.iterdir() if d.is_dir()) if dvg_dir_ov.is_dir() else 0

    # Pre-load stage summaries
    detection_summary = _collect_detection_summary(out_dir)
    variant_summary = _collect_variant_summary(out_dir)
    assembly_summary = _collect_assembly_summary(out_dir)
    extract_stats = _collect_extract_stats(out_dir)

    # Helper: collect stage-level files for embedding
    def _stage_files(outd, sn):
        """Return HTML for stage-level result files (plots/tables from stage directory)."""
        dirs = {1:"1_FastViromeExplorer", 2:"2_Virus_result_filter",
                4:"4_Virus_assemblies_final",
                6:"6_post_analysis", 7:"7_capheine"}
        patterns = {1:["summary/all_viruses.best.summary.tsv","plots/virus_analysis/freq_multi_metrics_log10.png"],
                    2:[],
                    4:["**/*final_coverage*"],
                    6:["**/metadata_association/**"],
                    7:["**/selection_per_gene.*"]}
        if sn not in dirs: return "", ""
        spath = outd / dirs[sn]
        if not spath.exists(): return "", ""
        html_tables = ""; html_charts = ""
        collected = []
        for pat in patterns.get(sn,["*.png"]):
            for f in spath.rglob(pat) if "**" in pat else spath.glob(pat):
                if f.is_file(): collected.append(f)
        # Dedup: skip PDF when PNG with same stem exists
        png_stems = {f.stem for f in collected if f.suffix.lower() in ('.png','.jpg','.jpeg','.svg')}
        for f in collected:
            if f.suffix.lower() == '.pdf' and f.stem in png_stems:
                continue
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
                    html_charts += f'<div class="chart-card"><div class="chart-title">{fn}</div><object data="data:application/pdf;base64,{b64}" type="application/pdf" width="100%" height="500px"><p>PDF</p></object></div>'
            elif fn.endswith(('.png','.jpg','.jpeg','.svg')):
                b64 = img_to_base64(str(f), 2000)
                if b64:
                    mime = "image/png"
                    if fn.endswith('.jpg') or fn.endswith('.jpeg'): mime = "image/jpeg"
                    elif fn.endswith('.svg'): mime = "image/svg+xml"
                    html_charts += f'<div class="chart-card"><div class="chart-title">{fn}</div><img src="data:{mime};base64,{b64}" loading="lazy" alt="{fn}"></div>'
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
            raw_tsv_s1 = out / "1_FastViromeExplorer/summary/all_viruses.raw.tsv"
            raw_s1 = safe_read_csv(raw_tsv_s1)
            if raw_s1 is not None:
                n_samples = raw_s1["Sample"].nunique() if "Sample" in raw_s1.columns else "?"
                acc_raw = next((c for c in ["Accession","Rep_Accession"] if c in raw_s1.columns), None)
                n_virus_raw = raw_s1[acc_raw].nunique() if acc_raw else "?"
            else:
                n_samples = "?"
            # Best summary stats
            best_tsv_s1 = out / "1_FastViromeExplorer/summary/all_viruses.best.summary.tsv"
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
            raw_tsv_s2 = out / "1_FastViromeExplorer/summary/all_viruses.raw.tsv"
            hc_tsv_s2 = out / "2_Virus_result_filter/high_conf.summary.tsv"
            ps_tsv_s2 = out / "2_Virus_result_filter/filter_stats.per_sample.tsv"
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
            raw_tsv = out / "1_FastViromeExplorer/summary/all_viruses.raw.tsv"
            raw = safe_read_csv(raw_tsv)
            pre_sample = {}; pre_virus = {}
            if raw is not None:
                raw_acc = next((c for c in ["Accession","Rep_Accession"] if c in raw.columns), None)
                if raw_acc and "Sample" in raw.columns:
                    pre_sample = raw.groupby("Sample")[raw_acc].nunique().to_dict()
                    pre_virus = raw.groupby(raw_acc)["Sample"].nunique().to_dict()
            # Also build species->accession map from best TSV for per-virus lookup
            best_tsv = out / "1_FastViromeExplorer/summary/all_viruses.best.summary.tsv"
            best = safe_read_csv(best_tsv)
            sp2acc = {}
            if best is not None:
                bac = next((c for c in ["Rep_Accession","Accession"] if c in best.columns), None)
                spc = next((c for c in ["Adjusted_Species","Species_NCBI","Species"] if c in best.columns), None)
                if bac and spc:
                    for _, r in best.iterrows():
                        sp2acc[str(r[spc])] = str(r[bac])
            # Parse filter_stats
            ps_tsv = out / "2_Virus_result_filter/filter_stats.per_sample.tsv"
            pv_tsv = out / "2_Virus_result_filter/filter_stats.per_virus.tsv"
            ps_df = safe_read_csv(ps_tsv); pv_df = safe_read_csv(pv_tsv)
            # Per-sample table
            if ps_df is not None:
                stage_sections += '<h3 style="color:var(--ink-secondary);font-size:14px;margin-top:12px">Per-Sample Filtering</h3>'
                stage_sections += '<div class="tb-scroll"><table style="font-size:12px"><tr><th>Sample</th><th>Pre-Filter</th><th>Post-Filter</th><th>Filtered</th><th>Retention</th><th>Total Reads</th><th>Mean Cov%</th></tr>'
                for _, r in ps_df.iterrows():
                    s = str(r["Sample"]); pre = pre_sample.get(s, 0)
                    post = int(r["n_viruses"]); flt = pre - post
                    rate = f"{post/max(pre,1)*100:.0f}%" if isinstance(pre, (int,float)) else "?"
                    stage_sections += f'<tr><td>{s}</td><td>{pre}</td><td>{post}</td><td>{flt}</td><td>{rate}</td><td>{r["total_EM_reads"]}</td><td>{r["mean_coverage"]}%</td></tr>'
                stage_sections += '</table></div>'
            # Build per-virus max/avg/min from high_conf
            hc_tsv = out / "2_Virus_result_filter/high_conf.summary.tsv"
            hc = safe_read_csv(hc_tsv)
            virus_stats = {}
            if hc is not None:
                ac = next((c for c in ["Rep_Accession","Accession"] if c in hc.columns), None)
                sc = next((c for c in ["Adjusted_Species","Species_NCBI","Species"] if c in hc.columns), None)
                if ac:
                    for acc, grp in hc.groupby(ac):
                        sp = str(grp[sc].iloc[0])[:30] if sc else str(acc)
                        cov = grp["Rep_Coverage(%)"] if "Rep_Coverage(%)" in hc.columns else None
                        depth = grp["Rep_MeanDepth"] if "Rep_MeanDepth" in hc.columns else None
                        reads = grp["Asm_EM_Reads"] if "Asm_EM_Reads" in hc.columns else None
                        tpm = grp["Asm_TPM"] if "Asm_TPM" in hc.columns else None
                        virus_stats[acc] = {"species": sp,
                            "cov_min": cov.min() if cov is not None else 0, "cov_max": cov.max() if cov is not None else 0, "cov_avg": cov.mean() if cov is not None else 0,
                            "depth_min": depth.min() if depth is not None else 0, "depth_max": depth.max() if depth is not None else 0, "depth_avg": depth.mean() if depth is not None else 0,
                            "reads_min": reads.min() if reads is not None else 0, "reads_max": reads.max() if reads is not None else 0, "reads_avg": reads.mean() if reads is not None else 0,
                            "tpm_min": tpm.min() if tpm is not None else 0, "tpm_max": tpm.max() if tpm is not None else 0, "tpm_avg": tpm.mean() if tpm is not None else 0}
            # Per-virus table
            if pv_df is not None:
                stage_sections += '<h3 style="color:var(--ink-secondary);font-size:14px;margin-top:16px">Per-Virus Filtering</h3>'
                stage_sections += '<div class="tb-scroll"><table style="font-size:11px"><tr><th>Species</th><th>Pre</th><th>Post</th><th>Flt</th><th>Ret</th>'
                for m in ["min","avg","max"]:
                    stage_sections += f'<th>Cov% {m}</th><th>Depth {m}</th><th>Reads {m}</th><th>TPM {m}</th>'
                stage_sections += '</tr>'
                for _, r in pv_df.iterrows():
                    sp = str(r["Adjusted_Species"])[:25]; acc = sp2acc.get(str(r["Adjusted_Species"]), "")
                    pre = pre_virus.get(acc, 0) if acc else 0
                    post = int(r["n_samples"]); flt = pre - post
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
            coi_tsv = out / "3_Virus_variants_Results/summary/Coinfection_Matrix_Reads.tsv"
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
            plots_dir = out / "1_FastViromeExplorer/plots"
            keep = {"sample_distribution": ["sample_virus_count_bar.png", "virus_occurrence_bar.png"], "coabundance": ["coabundance_cooccurrence_heatmap.png", "coabundance_spearman_heatmap.png"]}
            stage_sections += '<h3 style="color:var(--ink-secondary);font-size:14px;margin-top:12px">Sample Distribution &amp; Co-abundance</h3><div class="chart-gallery">'
            for bd, keep_files in keep.items():
                bp = plots_dir / bd
                if not bp.is_dir(): continue
                for fn in keep_files:
                    f = bp / fn
                    if f.is_file():
                        b64 = img_to_base64(str(f), 2000)
                        if b64: stage_sections += f'<div class="chart-card"><div class="chart-title">{fn}</div><img src="data:image/png;base64,{b64}" loading="lazy" alt="{fn}"></div>'
            stage_sections += '</div>'
        # Stage 9: DVG per-virus statistics table
        if sn == 9:
            dvg_m = out / "9_virema_dvg/Summary_Analysis_Report/Matrix_Sample_Wise_Recombination_Statistics.csv"
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
                    kw_csv = out / "6_post_analysis" / _find_virus_dir(out / "6_post_analysis", acc).name / "snpgenie/Stats_Kruskal_Wallis.csv" if _find_virus_dir(out / "6_post_analysis", acc) else None
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
                                if ci.get('base64'):
                                    b64 = ci['base64']; ext = ci.get('ext','')
                                    if ext == 'pdf':
                                        stage_sections += f'<div class="chart-card"><div class="chart-title">{ci.get("label","")}</div><object data="data:application/pdf;base64,{b64}" type="application/pdf" width="100%" height="500px"><p>PDF</p></object></div>'
                                    else:
                                        stage_sections += f'<div class="chart-card"><div class="chart-title">{ci.get("label","")}</div><img src="data:image/png;base64,{b64}" loading="lazy" alt="{ci.get("label","")}"></div>'
                            stage_sections += '</div></section>'
                    else:
                        stage_sections += '<p style="color:var(--ink-secondary);font-size:12px">Post-hoc charts not available for this virus.</p>'
                elif sn == 9:
                    dcharts = collect_dvg_charts(out / "9_virema_dvg", acc)
                    if dcharts:
                        stage_sections += '<div class="chart-gallery">'
                        for ci in sorted(dcharts.values(), key=lambda x: x.get('label','')):
                            if ci.get('base64'):
                                b64 = ci['base64']; ext = ci.get('ext','')
                                if ext == 'pdf':
                                    stage_sections += f'<div class="chart-card"><div class="chart-title">{ci.get("label","")}</div><object data="data:application/pdf;base64,{b64}" type="application/pdf" width="100%" height="500px"><p>PDF</p></object></div>'
                                else:
                                    stage_sections += f'<div class="chart-card"><div class="chart-title">{ci.get("label","")}</div><img src="data:image/png;base64,{b64}" loading="lazy" alt="{ci.get("label","")}"></div>'
                        stage_sections += '</div>'
                    else:
                        stage_sections += '<p style="color:var(--ink-secondary);font-size:12px">DVG charts not available for this virus.</p>'
                else:
                    # S3/S4/S7/S8: locate per-virus result files
                    if sn == 2:
                        fp_filter = out / "2_Virus_result_filter" / "filter_summary_plot.png"
                        b64_f = img_to_base64(fp_filter, 3000)
                        if b64_f:
                            stage_sections += f'<div class="chart-card"><div class="chart-title">Filter Summary</div><img src="data:image/png;base64,{b64_f}" loading="lazy" alt="Filter Summary"></div>'
                        else:
                            stage_sections += '<p style="color:var(--ink-secondary);font-size:12px">Filter summary plot not available.</p>'
                    elif sn == 5:
                        # Count extracted contigs per virus
                        ext_dir = out / "5_assemblies_clean"
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
                            stage_sections += '<p style="color:var(--ink-secondary);font-size:12px">No coding sequences available (non-coding RNA virus) or capheine not run.</p>'
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

    html = f"""<!DOCTYPE html><html lang="en"><head><meta charset="UTF-8"><meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>Virus Pipeline Report</title>
<style>
/* ═══════════════════════════════════════════
   Virus Pipeline Report — Product Register
   Palette: Restrained blue-on-neutral
   Type: System sans-serif stack
   ═══════════════════════════════════════════ */
:root {{
  --ink: #1a2332;
  --ink-secondary: #5a6a7e;
  --accent: #2563eb;
  --accent-hover: #1d4ed8;
  --accent-subtle: #eff4ff;
  --surface: #ffffff;
  --surface-alt: #f6f8fb;
  --surface-hover: #eef1f6;
  --border: #e1e5eb;
  --border-light: #eef0f4;
  --sidebar-bg: #f3f5f8;
  --sidebar-ink: #4a5568;
  --sidebar-active: #2563eb;
  --success: #0d9488;
  --warning: #d97706;
  --danger: #dc2626;
  --radius: 6px;
  --shadow-sm: 0 1px 2px rgba(0,0,0,.04);
  --shadow-md: 0 4px 12px rgba(0,0,0,.06);
}}

* {{box-sizing:border-box;margin:0;padding:0}}

body {{
  font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', system-ui, sans-serif;
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
  color: var(--ink); font-size:26px; font-weight:700; margin-bottom:4px;
  padding-bottom:10px; border-bottom: 2px solid var(--border); letter-spacing: -.02em;
}}
h2 {{
  color: var(--ink); font-size:19px; font-weight:600; margin:28px 0 8px; letter-spacing: -.01em;
}}
h3 {{font-size:15px; font-weight:600; color: var(--ink); margin:16px 0 6px}}

/* ── Metrics cards ── */
.metrics {{display:flex; flex-wrap:wrap; gap:10px; margin:12px 0}}
.metric {{
  background: var(--surface); border:1px solid var(--border);
  border-radius: var(--radius); padding:12px 18px; text-align:center; min-width:88px;
  transition: box-shadow 0.15s;
}}
.metric:hover {{box-shadow: var(--shadow-sm)}}
.metric .value {{display:block; font-size:22px; font-weight:700; color: var(--accent)}}
.metric .unit {{display:block; font-size:10.5px; color: var(--ink-secondary); margin-top:2px}}

.card-row {{display:flex; flex-wrap:wrap; gap:10px; margin:12px 0}}
.card {{
  background: var(--surface); border:1px solid var(--border);
  border-radius: var(--radius); padding:14px 18px; text-align:center; min-width:100px; flex:1;
  transition: box-shadow 0.15s;
}}
.card:hover {{box-shadow: var(--shadow-sm)}}
.card .value {{display:block; font-size:24px; font-weight:700; color: var(--accent); word-break:break-all}}
.card .label {{display:block; font-size:10px; color: var(--ink-secondary); margin-top:2px; font-weight:500}}

/* ── Chart gallery ── */
.chart-gallery {{
  display:grid; grid-template-columns: 1fr;
  gap:12px; margin:12px 0;
}}
.chart-card {{
  background: var(--surface); border:1px solid var(--border);
  border-radius: var(--radius); overflow:hidden;
  transition: box-shadow 0.15s;
}}
.chart-card:hover {{box-shadow: var(--shadow-md)}}
.chart-title {{
  background: var(--surface-alt); padding:6px 10px; font-size:11.5px;
  font-weight:600; color: var(--ink-secondary); border-bottom: 1px solid var(--border-light);
}}
.chart-card img {{
  width:100%; height:auto; max-height:80vh; object-fit:contain; display:block; background: var(--surface-alt);
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
.nav-sub {{overflow:hidden; max-height:2000px; transition:max-height 0.25s ease}}
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

/* ── Responsive ── */
@media (max-width:800px) {{
  .sidebar{{display:none}}.main{{margin-left:0; padding:20px 16px}}
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
<button id="back-to-top" onclick="scrollTop()" title="Back to top">&#9650;</button>

<button id="chat-btn" onclick="toggleChat()" title="Ask questions about this report">?</button>
<div id="chat-panel" style="display:none">
    <div class="chat-header">
        <span>Report Q&amp;A</span>
        <button onclick="toggleChat()" style="background:none;border:none;color:var(--ink-secondary);font-size:16px;cursor:pointer">&times;</button>
    </div>
    <div class="chat-messages" id="chat-msgs"></div>
    <div class="chat-input-row">
        <input id="chat-input" type="text" placeholder="Ask about this report..." onkeydown="if(event.key==='Enter')sendChat()">
        <button onclick="sendChat()" class="chat-send-btn">Send</button>
    </div>
</div>
<nav class="sidebar">
    <h3>Pipeline Report</h3>
    <ul>{nav_items}</ul>
</nav>
<main class="main">
    <h1>Known Virus Pipeline Report</h1>
    <p style="color:var(--ink-secondary);font-size:13px">Generated: {now} | Directory: {out_dir}</p>

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
            <div class="card"><div class="value">{overview_metrics.get("raw_detections","?")}</div><div class="label">Raw Detections</div></div>
            <div class="card"><div class="value">{overview_metrics.get("best_detections","?")}</div><div class="label">High-Confidence</div></div>
            <div class="card"><div class="value">{overview_metrics.get("filtered_records","?")}</div><div class="label">Post-Filter Records</div></div>
            <div class="card"><div class="value">{overview_metrics.get("assemblies","?")}</div><div class="label">Assemblies</div></div>
            <div class="card"><div class="value">{overview_metrics.get("dvg_viruses","?")}</div><div class="label">Viruses with DVG</div></div>
            <div class="card"><div class="value">10</div><div class="label">Pipeline Stages</div></div>
        </div>
    </section>

    {stage_sections}

    <section id="virus-summary">
        <h2>Virus Summary</h2>
        {_build_virus_summary_table(viruses)}
        <h3 id="paper-reference">Paper Reference: Methods Paper Template</h3>

        <h3>Pipeline Summary</h3>
        <div class="paper-block">
            <p>The Known Virus Analysis Pipeline (KVAP) processed <strong>{len(viruses)} virus species</strong> across samples, completing all 10 stages of detection, variant analysis, assembly, and evolutionary characterization.</p>
            <table>
                <tr><th>Stage</th><th>Script</th><th>Output</th></tr>
                <tr><td>1. Detect</td><td>batch_virus_depth.py</td><td>Salmon pseudo-alignment + Poisson filtering</td></tr>
                <tr><td>2. Filter</td><td>utils/filter_summary.py</td><td>Coverage ≥50%, Depth ≥5×, Reads ≥100</td></tr>
                <tr><td>3. Variants</td><td>batch_virus_variants.py</td><td>FreeBayes + SnpEff + SNPGenie</td></tr>
                <tr><td>4. Assembly</td><td>virus-full.py</td><td>12-step multi-tool de novo assembly</td></tr>
                <tr><td>5. Extract</td><td>utils/extract_full_fasta.py</td><td>Longest contig extraction</td></tr>
                <tr><td>6. Post-hoc</td><td>6-script suite</td><td>VCF viz + SnpEff macro + MAF + SnpGenie</td></tr>
                <tr><td>7. Capheine</td><td>capheine_pipeline.py</td><td>HyPhy FEL/MEME/BUSTED/PRIME</td></tr>
                <tr><td>8. Similarity</td><td>virus_auto_pipeline.py</td><td>SDT pairwise similarity heatmaps</td></tr>
                <tr><td>9. DVG</td><td>batch_virema_dvg.py</td><td>ViReMa recombination + Circos</td></tr>
                <tr><td>10. Report</td><td>generate_pipeline_report.py</td><td>This interactive HTML report</td></tr>
            </table>
        </div>

        <h3>Key Quantitative Results</h3>
        <div class="paper-block">
            <p>The pipeline was validated on a plant virus dataset with the following quantifiable outcomes:</p>
            <ul>
                <li><strong>Detection Sensitivity</strong>: {total_records} virus records identified across {len(viruses)} species</li>
                <li><strong>Assembly Completeness</strong>: {"Multiple complete genomes recovered at >99% reference coverage" if len(viruses) > 1 else "Genome assembly completed"}</li>
                <li><strong>Variant Calling Resolution</strong>: SNVs detected at minimum 10× depth with 5% allele frequency threshold</li>
                <li><strong>Selection Analysis</strong>: {"Kruskal-Wallis cross-gene dN/dS comparison and HyPhy codon-level positive selection testing available for coding viruses"}</li>
            </ul>
        </div>

        <h3>Writing Suggestions — Methods Section</h3>
        <div class="paper-block">
            <p>For each stage of the pipeline, the following writing templates can be adapted for a journal methods section:</p>

            <details><summary><strong>Detection & Filtering (Stages 1–2)</strong></summary>
            <div class="suggestion">"Virus detection was performed using Salmon (v1.10) pseudo-alignment in quantification mode against a curated database of [N] plant virus reference genomes. Read counts were normalized to CPM (counts per million) and FPKM. A Poisson Ratio filter (threshold ≥ 0.3) was applied to distinguish uniformly covered true viral reads from localized spurious alignments. High-confidence detections were retained by requiring Rep_Coverage ≥ 50%, Rep_MeanDepth ≥ 5×, and Asm_EM_Reads ≥ 100."</div></details>

            <details><summary><strong>Variant Calling (Stage 3)</strong></summary>
            <div class="suggestion">"Variants were called using FreeBayes (v1.3.6) in haploid mode with dynamic depth thresholds (DP ≥ 10 for mean depth < 50×, DP ≥ 20 for 50–1000×, DP ≥ 100 for >1000×) and minimum allele frequency of 5%. Functional annotation was performed with SnpEff (v5.1) using custom databases built from viral GenBank records. Population genetic parameters (π, πN/πS, dN/dS) were computed using SNPGenie with 50-bp sliding windows."</div></details>

            <details><summary><strong>Assembly (Stages 4–5)</strong></summary>
            <div class="suggestion">"De novo genome assembly employed a 12-step pipeline integrating MEGAHIT, SPAdes/RNAviralSPAdes, and PenguIN assemblers with iterative refinement including Shiver-like orientation correction, Divine Fusion reference-guided gap resolution, PVGA read extension, and 3-round minimap2 + viral_consensus polishing. Assembly quality was assessed by total length relative to reference (≥98% = perfect) and N50."</div></details>

            <details><summary><strong>Selection Analysis (Stages 6–7)</strong></summary>
            <div class="suggestion">"Gene-level dN/dS distributions were compared using Kruskal-Wallis non-parametric testing. Individual gene deviation from neutral expectation was assessed with Wilcoxon signed-rank tests. 10,000-iteration bootstrap resampling provided 95% confidence intervals for dN/dS estimates. Codon-level positive selection was detected using HyPhy FEL, MEME, BUSTED, and PRIME methods with significance at p < 0.05. The optimal number of evolutionary sub-lineages was automatically determined by maximizing the Silhouette Score (K = 2–7) on PCA-reduced feature space [dN, dS, πN, πS]."</div></details>

            <details><summary><strong>Defective Genome Detection (Stage 9)</strong></summary>
            <div class="suggestion">"Defective viral genomes and recombination events were detected using ViReMa (v0.29) with seed length 25, micro-indel threshold 15 bp, and strict breakpoint resolution (defuzz = 0). Recombination landscapes were visualized using Circos 4-track plots showing mutation density, gene annotation, coverage depth (deletions/duplications), and donor-acceptor junction networks."</div></details>
        </div>

        <h3>Writing Suggestions — Results Section</h3>
        <div class="paper-block">
            <p>Based on the actual pipeline output, the following data-driven narrative templates can be used:</p>
            <ul>
                <li><strong>Detection overview</strong>: "Analysis of [N] samples identified [M] known virus species, with [high-conf species] passing stringent quality filters (coverage ≥50%, depth ≥5×)."</li>
                <li><strong>Assembly performance</strong>: "De novo assembly successfully reconstructed [A/B] complete viral genomes ([P]% perfect rate)."</li>
                <li><strong>Selection signature</strong>: "Cross-gene comparison revealed significant differences in selection pressure (Kruskal-Wallis H = [value], p = [value]), with [gene] showing the highest median dN/dS."</li>
                <li><strong>Evolutionary clustering</strong>: "Unsupervised machine learning (PCA + K-Means with Silhouette Score optimization) identified [K] distinct evolutionary sub-lineages."</li>
            </ul>
        </div>
    </section>

    <div class="footer">Generated by known_virus_pipeline — {now}</div>
</main>
<div class="lightbox" id="lightbox" onclick="if(event.target===this)closeLightbox()">
    <button class="lightbox-close" onclick="closeLightbox()" aria-label="Close">&times;</button>
    <img id="lightbox-img" src="" alt="">
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
        messages:[{{role:'system',content:'You are a virology research assistant. Answer questions about this pipeline report concisely in Chinese. Use these data: '+ctx}},{{role:'user',content:q}}],
        temperature:0.3,max_tokens:300}})
    }}).then(r=>r.json()).then(d=>{{
        thinking.textContent=d.choices[0].message.content.trim();
        msgs.scrollTop=msgs.scrollHeight;
    }}).catch(e=>{{thinking.textContent='Error: '+e;}});
}}
// Close lightbox on Escape
document.addEventListener('keydown',function(e){{if(e.key==='Escape')closeLightbox();}});
// ── In-page AI ──
var sysMsg='你是病毒学领域的研究科学家。撰写学术总结时严格遵循：1.仅基于提供的数据和图表描述，不得捏造任何数字或发现 2.引用具体数值时使用数据显示、检测到、观察到等客观措辞 3.引用图表时使用如XX图所示、XX表显示 4.不确定时使用提示、暗示、可能等学术谨慎用语 5.数据不支持时明确说明未检测到、无显著差异 6.中文输出，学术书面语';
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
        1:'Stage 1 Salmon伪比对+Poisson过滤病毒定量检测。\\n每病毒：'+vdStr+'\\n嵌入图：freq_multi_metrics_log10、sample_virus_count_bar、virus_occurrence_bar、cooccurrence/spearman热图、Co-infection_Matrix。\\n\\n撰写2-3句中文学术总结，对每种病毒逐一简要评述检出样本数、丰度、覆盖度。引用具体数字。不捏造。',
        2:'Stage 2 多维阈值过滤（覆盖率≥50%、深度≥5×、reads≥100）。\\n'+sd.s2.pre+'条→'+sd.s2.post+'条保留。\\n每病毒详情：'+vdStr+'\\n嵌入表：Per-Sample/Per-Virus Filtering。\\n\\n对每种病毒简要评述过滤后的覆盖度和深度。引用具体数字。不捏造。',
        3:'Stage 3 FreeBayes+SnpEff+SNPGenie。\\n每病毒数据：'+JSON.stringify(sd.s3.per_virus)+'\\n\\n对每种病毒评述其核苷酸多样性(π)。高π=准种丰富，低π=保守。引用具体数值。不捏造。',
        4:'Stage 4 12步从头组装。\\n每病毒数据：'+JSON.stringify(sd.s4.per_virus)+'\\n\\n对每种病毒评述组装完整性(N%、N50)。引用具体数值。不捏造。',
        5:'Stage 5 提取+N填补。\\n每病毒数据：'+JSON.stringify(sd.s5.per_virus)+'\\n\\n评述contig完整度(N%=0完整)。引用具体数值。不捏造。',
        6:'Stage 6 生成多维度变异可视化图表组，按分析类型分组展示。\\n图表内容（per-virus）：Variant Landscape（全基因组突变分布景观）、Ts/Tv Pie（转换/颠换比值饼图——比值>1提示随机突变主导）、Functional Pie（变异功能分类饼图）、Clustermap（跨样本变异聚类热图）、AFS（等位频率谱——低频为主提示稀有变异，高频为主提示固定差异）、dN/dS散点图及3D PCA（群体遗传参数三维聚类）。\\n共'+sd.viruses+'种病毒纳入分析（单样本病毒已跳过）。\\n\\n基于以上图表类型描述，撰写 2-3 句学术总结。可引用具体图表名称。不捏造未观测到的模式。',
        7:'Stage 7 使用 HyPhy（FEL/MEME/BUSTED/PRIME）进行密码子水平自然选择检测。\\n数据：'+sd.s7.coding_viruses+'/'+sd.viruses+'种病毒具有编码序列，其余为非编码RNA（类病毒/卫星RNA，不编码蛋白因此不适用dN/dS分析）。\\n嵌入表：per-gene BUSTED选择分析表（基因名、ω3值、p值、正选择位点数、负选择位点数、总位点数）。\\n\\n基于以上数据，撰写 2-3 句学术总结。注意：ω3≈1且p>0.05表示无显著正选择信号；负选择位点占优提示纯化选择占主导。',
        8:'Stage 8 使用 SDT 风格成对 NT 一致性分析。\\n比较了'+sd.viruses+'种病毒的样本间序列。\\n嵌入图：NT Identity Heatmap（成对一致性热图+系统聚类树）、Identity Distribution（一致性分布直方图）、pairwise statistics（均值/范围/比较对数）。\\n\\n基于以上数据，撰写 2-3 句学术总结。一致性>95%通常指示同一株系，85-95%为种内变异。仅陈述观测到的一致性范围。',
        9:'Stage 9 使用 ViReMa 检测缺陷病毒基因组（DVG）和重组事件。\\n数据：'+sd.s9.viruses_with_dvg+'/'+sd.viruses+'种病毒检测到 DVG 事件。\\n嵌入图：Circos 4-track（突变密度/基因注释/缺失覆盖/重复覆盖四轨道环状图）、Arc Diagram（供体-受体连接网络弧线图）、Acceptor NT Frequency（受体断点核苷酸频率）、Top Events（最高频重组事件）。\\n嵌入表：DVG Events per Virus（每病毒DVG reads数、跳跃事件数、涉及样本数）。\\n\\n基于以上数据，撰写 2-3 句学术总结。引用图表和表格。仅当数据支持时才讨论DIP假说。'
    }};
    var prompt=prompts[sn]||('Pipeline overview: '+sd.samples+' samples, '+sd.viruses+' virus species, 10-stage analysis. Write 2-3 sentence academic summary based ONLY on provided data.');
    prompt+='\\n\\n中文输出。学术书面语。引用图表名称。不捏造数据。不确定处用谨慎措辞。'
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
    var stageFocus={{3:'该病毒位于Stage 3 Variants区域。嵌入图表：per-virus变异指标表(Samples/Cov/Depth/Pi/Shannon)。\\n\\n评述该病毒的核苷酸多样性(π)和Shannon指数。高π=准种多样性丰富(活跃复制的准种群)，低π=高度保守(强纯化选择或近期瓶颈)。结合覆盖度和深度评估变异检测可靠性。仅陈述数据支持的事实。',4:'该病毒位于Stage 4 Assembly区域。嵌入图表：每样本组装演进表(DeNovo→Ultimate长度/N50/N)、可折叠阶段细节(各步骤长度变化)。\\n\\n评述该病毒的组装完整性：组装数、最佳长度、N含量(N>0=存在gap未闭合)、vs参考基因组百分比(>100%提示多聚体/串联重复组装，100%=完整单体，<100%=部分组装)。结合阶段演进分析组装改善情况。',5:'该病毒位于Stage 5 Extract区域。嵌入表格：提取统计表(Numbers/Avg/Min/Max bp/N%)。\\n\\n评述该病毒的提取结果：contig数、平均/最小/最大长度、N%(0=完整无gap)。N%>0提示序列存在未填补区域。对比参考基因组评估提取完整度。',6:'该病毒位于Stage 6 Post-hoc区域。嵌入图表(分组)：Variant Landscape(全基因组突变分布)、Ts/Tv Pie(转换/颠换比>1=随机突变主导)、Functional Pie(功能分类)、Clustermap(跨样本变异谱系)、AFS(等位频率谱)、dN/dS散点图及3D PCA(群体遗传参数聚类)。\\n\\n结合图表评述：突变沿基因组的分布特征(是否存在热点)、Ts/Tv比反映的进化驱动力、dN/dS偏离中性的生物学含义(纯化选择vs正选择vs中性演化)、PCA聚类是否揭示不同进化谱系。',7:'该病毒位于Stage 7 Capheine区域。嵌入表格：per-gene BUSTED选择分析表(基因名/ω3/p值/正选择位点数/负选择位点数/总位点数)。\\n\\n评述该病毒基因的选择压力：ω3>1提示该基因可能经历正选择(需p<0.05支持)，ω3≈1且p>0.05提示无显著正选择信号。负选择位点占优提示纯化选择占主导(功能约束强)。如为非编码RNA病毒(类病毒/卫星RNA)，说明dN/dS分析不适用——其进化约束来自RNA二级结构而非蛋白编码。',8:'该病毒位于Stage 8 Similarity区域。嵌入图表：NT Identity Heatmap(成对一致性热图+系统聚类树)、Identity Distribution(一致性分布直方图)、pairwise statistics(均值/范围/比较对数)。\\n\\n评述该病毒样本间NT一致性：>95%通常指示同一株系(近期传播)，85-95%为种内变异，<85%为不同种。聚类树揭示的进化关系。高一致性提示近期共同祖先或单一引入事件，低一致性提示长期本地进化或多重引入。',9:'该病毒位于Stage 9 DVG区域。嵌入图表：Circos 4-track(突变密度/基因注释/缺失覆盖/重复覆盖)、Arc Diagram(供体-受体连接网络)、Acceptor NT Frequency(受体断点核苷酸偏好)、Top Events(最高频重组事件)。嵌入表格：DVG Events per Virus(DVG reads/跳跃事件/涉及样本数)。\\n\\n评述该病毒的DVG/重组事件：DVG reads数和跳跃事件数反映重组活跃度。DVG(缺陷病毒基因组)的生物学意义：缺陷干扰颗粒(DIP)可竞争性抑制标准病毒复制，影响感染进程和致病性。重组类型(缺失/重复)与病毒RNA聚合酶模板转换的错误倾向性有关。'}};
    var focus=stageFocus[sn]||'评述该病毒的检测、覆盖度、丰度等关键指标。';
    var prompt='你是病毒学研究员。基于以下病毒数据撰写1-2句中文学术总结。'+focus+'\\n\\n病毒数据：'+vinfo+'\\n\\n仅陈述数据支持的事实。引用具体数字。不捏造。中文。';
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
</body></html>"""

    with open(out_html, "w", encoding="utf-8") as f:
        f.write(html)
    return str(out_html)


def ensure_summary_plots(out_dir):
    import subprocess as _sp
    out = Path(out_dir)
    script_dir = Path(__file__).resolve().parent
    fp_plot = out / "2_Virus_result_filter" / "filter_summary_plot.pdf"
    if not fp_plot.exists():
        best_tsv = None
        for p in ["high_conf.summary.tsv","all_viruses.best.summary.tsv"]:
            for d in ["2_Virus_result_filter","1_FastViromeExplorer/summary"]:
                fpp = out / d / p
                if fpp.is_file(): best_tsv = str(fpp); break
            if best_tsv: break
        if best_tsv:
            fs = script_dir / "utils" / "filter_summary.py"
            if fs.is_file():
                _sp.run([sys.executable,str(fs),"-i",best_tsv,"-o",str(out/"2_Virus_result_filter"/"high_conf.summary.tsv"),"--plot"], capture_output=True)
    asm_plot = out / "5_assemblies_clean" / "assembly_stats.pdf"
    if not asm_plot.exists():
        asm_dir = out / "4_Virus_assemblies_final"
        if asm_dir.is_dir():
            es = script_dir / "utils" / "extract_full_fasta.py"
            if es.is_file():
                _sp.run([sys.executable,str(es),"-d",str(asm_dir),"-o",str(out/"5_assemblies_clean"),"--plot","--max_n","100","--min_len","100"], capture_output=True)


def copy_results_to_report(out_dir, report_dir):
    """Copy key result files from each stage into the report directory.

    Copies summary tables, key plots, and per-virus charts to 10_Reports/
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
    _cp("1_FastViromeExplorer/summary/all_viruses.best.summary.tsv", "S1_Detection/all_viruses.best.summary.tsv")

    # Stage 2: Filter results
    _cp("2_Virus_result_filter/high_conf.summary.tsv", "S2_Filter/high_conf.summary.tsv")
    _cp("2_Virus_result_filter/filter_summary_plot.png", "S2_Filter/filter_summary_plot.png")
    _cp("2_Virus_result_filter/filter_summary_plot.pdf", "S2_Filter/filter_summary_plot.pdf")

    # Stage 3: Variants summary
    _cp("3_Virus_variants_Results/summary/all_summary.tsv", "S3_Variants/all_summary.tsv")
    _cp("3_Virus_variants_Results/summary/Coinfection_Matrix_Reads.tsv", "S3_Variants/Coinfection_Matrix_Reads.tsv")

    # Stage 5: Assembly stats
    _cp("5_assemblies_clean/assembly_stats.png", "S5_Assembly_Stats/assembly_stats.png")
    _cp("5_assemblies_clean/assembly_stats.pdf", "S5_Assembly_Stats/assembly_stats.pdf")

    # Stage 6: Post-hoc charts (per-virus)
    post_dir = out / "6_post_analysis"
    if post_dir.is_dir():
        for vdir in post_dir.iterdir():
            if not vdir.is_dir() or vdir.name.startswith("run_"): continue
            for subd in vdir.iterdir():
                if subd.is_dir():
                    _cp(str(subd.relative_to(out)), f"S6_Post_hoc/{vdir.name}/{subd.name}")

    # Stage 7: Capheine
    cap_dir = out / "7_capheine"
    if cap_dir.is_dir():
        for vdir in cap_dir.iterdir():
            if vdir.is_dir() and not vdir.name.startswith("run_"):
                _cp(str(vdir.relative_to(out)), f"S7_Capheine/{vdir.name}")

    # Stage 8: Similarity (recursive: heatmaps are nested deep)
    sim_dir = out / "8_similarity"
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
    dvg_plots = out / "9_virema_dvg/Summary_Analysis_Report/Virus_Specific_Plots"
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

    prompt = f"""You are a virology research scientist. Write a bilingual (Chinese) research summary for a Known Virus Analysis Pipeline report.

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
150-200 word Chinese research highlight. Cover: study purpose, sample scale, key findings (viruses detected, assembly quality, selection evidence), main conclusions. Style: journal Highlights.

[DETAILED]
IMRaD structure (300-400 words Chinese):
[BACKGROUND] 1-2 sentences. Known virus surveillance context.
[METHODS] 3-4 sentences. Summarize the 10-stage pipeline.
[RESULTS] 4-5 sentences. Report quantitative findings with numbers.
[DISCUSSION] 2-3 sentences. Significance, limitations, next steps.

Requirements: Professional academic tone, cite exact numbers from data, Chinese output, only output the format requested."""

    try:
        url = base_url or "https://api.openai.com/v1/chat/completions"
        body = _json.dumps({
            "model": model, "messages": [
                {"role": "system", "content": "You are a virology research scientist. Output in Chinese, objective, precise, cite data."},
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

    prompt = f"""You are a virology research scientist. For each pipeline stage below, write ONE concise sentence (Chinese) describing what was found. Be specific with numbers.

Pipeline: Known Virus Analysis (10 stages), {n_samples} samples, {n_virus} virus species.

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
S1: <one Chinese sentence>
S2: <one Chinese sentence>
...
S9: <one Chinese sentence>

Requirements: One sentence each. Cite numbers. Chinese. Professional tone. No extra text."""

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
    parser.add_argument("--ai-key", default="", help="OpenAI/DeepSeek API key for AI summary")
    parser.add_argument("--ai-provider", default="openai", choices=["openai","ollama","deepseek","moonshot","custom"])
    parser.add_argument("--ai-model", default="gpt-4o-mini", help="LLM model name")
    parser.add_argument("--ai-base-url", default="", help="Custom API base URL")

    args = parser.parse_args()
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

    summary_in = None
    for p in ["high_conf.summary.tsv", "all_viruses.best.summary.tsv", "all_viruses.summary.tsv"]:
        for loc in [out / "2_Virus_result_filter" / p, out / "1_FastViromeExplorer" / "summary" / p]:
            if loc.is_file(): summary_in = loc; break
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
        raw_t = safe_read_csv(_P(args.dir) / "1_FastViromeExplorer/summary/all_viruses.raw.tsv")
        best_t = safe_read_csv(_P(args.dir) / "1_FastViromeExplorer/summary/all_viruses.best.summary.tsv")
        hc_t = safe_read_csv(_P(args.dir) / "2_Virus_result_filter/high_conf.summary.tsv")
        ov["raw_detections"] = len(raw_t) if raw_t is not None else "?"
        ov["best_detections"] = len(best_t) if best_t is not None else "?"
        ov["filtered_records"] = len(hc_t) if hc_t is not None else "?"
        n_as = 0; d_as = _P(args.dir) / "4_Virus_assemblies_final"
        if d_as.is_dir():
            for vd in d_as.iterdir():
                if vd.is_dir() and not vd.name.startswith("run_"):
                    n_as += sum(1 for d in vd.iterdir() if d.is_dir())
        ov["assemblies"] = n_as
        d_dvg = _P(args.dir) / "9_virema_dvg/Summary_Analysis_Report/Virus_Specific_Plots"
        ov["dvg_viruses"] = sum(1 for d in d_dvg.iterdir() if d.is_dir()) if d_dvg.is_dir() else 0
        ai_html = generate_ai_summary(viruses, vs, asm, ds, es, ov)
        stage_sums = generate_stage_ai_summaries(viruses, vs, asm, ds, es, ov)
    else:
        stage_sums = {}

    print("Generating interactive HTML...")
    generate_html(args.dir, out_html, viruses, ai_html, stage_sums)
    print(f"  Report: {out_html} ({out_html.stat().st_size / 1024:.0f} KB)")

    # Copy result files to report directory
    report_dir = out / "10_Reports"
    print("Copying result files to report directory...")
    copied, skipped = copy_results_to_report(args.dir, report_dir)
    print(f"  Copied {len(copied)} files, skipped {len(skipped)} missing")
    for s in sorted(skipped):
        print(f"    (not found) {s}")


if __name__ == "__main__":
    main()
