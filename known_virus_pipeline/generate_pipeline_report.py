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
from pathlib import Path
from datetime import datetime

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


def img_to_base64(path, max_kb=200):
    """Convert image to base64 for embedding. Skip if > max_kb."""
    if not path or not Path(path).exists():
        return None
    size_kb = Path(path).stat().st_size / 1024
    if size_kb > max_kb:
        return None
    with open(path, 'rb') as f:
        return base64.b64encode(f.read()).decode('utf-8')


def collect_charts(post_dir, virus_acc):
    """Collect key chart paths for a virus."""
    charts = {}
    vdir = Path(post_dir) / virus_acc
    if not vdir.exists():
        return charts

    chart_patterns = {
        'vcf_viz': [
            ('Figure1A_All_Variants_Landscape.png', '全基因组变异景观'),
            ('Figure2_TsTv_Pie.png', 'Ts/Tv 比率'),
            ('Figure5_AFS.png', '等位频率谱'),
            ('Figure8_PopGen_Dynamics.png', 'PopGen 滑动窗'),
        ],
        'snpeff_macro': [
            ('Figure_1_Manhattan_Mut_Landscape.pdf', '突变曼哈顿图'),
            ('Figure_2_Gene_Payload.pdf', '基因突变载荷'),
            ('Figure_3_IntraHost_Diversity.pdf', '准种多样性'),
        ],
        'maftools': [
            ('mafSummary_TCGA.pdf', 'MAF 突变类型'),
            ('Oncoplot.pdf', '突变瀑布图'),
        ],
        'snpgenie': [
            ('Fig03_InterHost_dNdS.png', 'dN vs dS 联合分布'),
            ('Fig05_Gene_dNdS_Stats.png', '每基因 dN/dS'),
            ('Fig06_Bootstrapped_dNdS.png', 'Bootstrap 显著性'),
            ('Fig11a_PCA_2D.png', '2D PCA 聚类'),
            ('Fig11b_PCA_3D.png', '3D PCA 聚类'),
        ],
    }

    for subdir, patterns in chart_patterns.items():
        sd = vdir / subdir
        if not sd.exists(): continue
        for fname, label in patterns:
            fp = sd / fname
            if fp.exists():
                b64 = img_to_base64(fp)
                charts[f"{subdir}_{fname}"] = {
                    'label': label,
                    'path': str(fp.relative_to(post_dir.parent)),
                    'base64': b64,
                    'ext': fp.suffix[1:],
                }
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


def generate_html(out_dir, out_html, viruses):
    """Generate interactive HTML report with sidebar navigation."""
    now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    out = Path(out_dir)
    post_dir = out / "5_post_analysis"

    total_records = sum(v.get("n_samples", 0) for v in viruses.values()) if isinstance(
        next(iter(viruses.values()), {}).get("n_samples", 0), int) else "?"

    # Build sidebar nav items
    nav_items = '<li class="nav-item"><a href="#overview" class="nav-link active">Global Overview</a></li>'
    for acc, data in sorted(viruses.items()):
        sp = data.get("species", acc)[:30]
        nav_items += f'<li class="nav-item"><a href="#virus-{acc}" class="nav-link">{sp}</a></li>'

    # Build per-virus sections with embedded charts
    virus_sections = ""
    for acc, data in sorted(viruses.items()):
        sp = data.get("species", acc)
        cpm = data.get("cpm", "?")
        cov = data.get("coverage", "?")
        depth = data.get("depth", "?")
        poisson = data.get("poisson", "?")
        reads = data.get("reads", "?")
        n = data.get("n_samples", "?")

        charts = collect_charts(post_dir, acc)

        # Build chart gallery
        chart_html = ""
        for chart_id, chart_info in charts.items():
            label = chart_info['label']
            if chart_info.get('base64'):
                ext = chart_info['ext']
                mime = 'image/png' if ext == 'png' else 'application/pdf'
                chart_html += f"""<div class="chart-card">
                    <div class="chart-title">{label}</div>
                    <img src="data:{mime};base64,{chart_info['base64']}" alt="{label}" loading="lazy" />
                </div>"""
            else:
                rel_path = chart_info['path']
                chart_html += f"""<div class="chart-card">
                    <div class="chart-title">{label}</div>
                    <div class="chart-placeholder">
                        <a href="../{rel_path}" target="_blank">Open {label} →</a>
                    </div>
                </div>"""

        virus_sections += f"""<section id="virus-{acc}">
            <h2>{sp}</h2>
            <p class="accession">{acc}</p>
            <div class="metrics">
                <div class="metric"><span class="value">{n}</span><span class="unit">Samples</span></div>
                <div class="metric"><span class="value">{cov:.1f}%</span><span class="unit">Coverage</span></div>
                <div class="metric"><span class="value">{cpm:.1f}</span><span class="unit">CPM</span></div>
                <div class="metric"><span class="value">{depth:.1f}x</span><span class="unit">Depth</span></div>
                <div class="metric"><span class="value">{poisson:.2f}</span><span class="unit">Poisson</span></div>
                <div class="metric"><span class="value">{reads:.0f}</span><span class="unit">Reads</span></div>
            </div>
            <div class="chart-gallery">
                {chart_html or "<p>No charts found in post-hoc output.</p>"}
            </div>
        </section>
        """

    html = f"""<!DOCTYPE html><html lang="en"><head><meta charset="UTF-8"><meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>Virus Pipeline Report</title>
<style>
* {{box-sizing:border-box;margin:0;padding:0}}
body {{font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',Arial,sans-serif;color:#333;display:flex;min-height:100vh}}
.sidebar {{position:fixed;left:0;top:0;width:260px;height:100vh;background:#1a1a2e;color:#e0e0e0;overflow-y:auto;padding:20px 0;z-index:100}}
.sidebar h3 {{padding:0 20px 15px;color:#3498db;font-size:16px;border-bottom:1px solid #333;margin-bottom:10px}}
.nav-item {{list-style:none}}
.nav-link {{display:block;padding:8px 20px;color:#bbb;text-decoration:none;font-size:13px;transition:all 0.2s}}
.nav-link:hover,.nav-link.active {{color:#fff;background:#16213e}}
.main {{margin-left:260px;padding:30px 40px;flex:1;max-width:1200px}}
h1 {{color:#2c3e50;font-size:28px;margin-bottom:5px;border-bottom:3px solid #3498db;padding-bottom:10px}}
h2 {{color:#2980b9;font-size:22px;margin:30px 0 10px}}
.accession {{color:#888;font-size:13px;margin-bottom:15px}}
.metrics {{display:flex;flex-wrap:wrap;gap:12px;margin:15px 0}}
.metric {{background:#f0f4f8;border-radius:8px;padding:12px 20px;text-align:center;min-width:80px}}
.metric .value {{display:block;font-size:22px;font-weight:bold;color:#2980b9}}
.metric .unit {{display:block;font-size:11px;color:#888;margin-top:2px}}
.chart-gallery {{display:grid;grid-template-columns:repeat(auto-fill,minmax(350px,1fr));gap:15px;margin:15px 0}}
.chart-card {{background:#fff;border:1px solid #e0e0e0;border-radius:8px;overflow:hidden;box-shadow:0 1px 3px rgba(0,0,0,.08)}}
.chart-title {{background:#f8f9fa;padding:8px 12px;font-size:13px;font-weight:600;color:#555;border-bottom:1px solid #e0e0e0}}
.chart-card img {{width:100%;height:auto;display:block}}
.chart-placeholder {{padding:40px;text-align:center;background:#fafafa}}
.chart-placeholder a {{color:#3498db;text-decoration:none;font-size:14px}}
.card-row {{display:flex;flex-wrap:wrap;gap:12px;margin:15px 0}}
.card {{background:#f0f4f8;border-radius:8px;padding:15px 20px;text-align:center;min-width:120px;flex:1}}
.card .value {{display:block;font-size:26px;font-weight:bold;color:#2980b9}}
.card .label {{display:block;font-size:11px;color:#888;margin-top:3px}}
.footer {{margin-top:40px;padding:20px;text-align:center;color:#aaa;font-size:12px;border-top:1px solid #eee}}
section {{margin-bottom:30px;padding-top:10px}}
@media (max-width:800px) {{.sidebar{{display:none}}.main{{margin-left:0}}}}
</style></head><body>
<nav class="sidebar">
    <h3>Pipeline Report</h3>
    <ul>{nav_items}</ul>
</nav>
<main class="main">
    <h1>Known Virus Pipeline Report</h1>
    <p style="color:#888;font-size:13px">Generated: {now} | Directory: {out_dir}</p>

    <section id="overview">
        <h2>Global Overview</h2>
        <div class="card-row">
            <div class="card"><div class="value">{len(viruses)}</div><div class="label">Virus Species</div></div>
            <div class="card"><div class="value">{total_records}</div><div class="label">Records</div></div>
            <div class="card"><div class="value">{out_dir}</div><div class="label">Output</div></div>
        </div>
    </section>

    {virus_sections}

    <div class="footer">Generated by known_virus_pipeline — {now}</div>
</main></body></html>"""

    with open(out_html, "w", encoding="utf-8") as f:
        f.write(html)
    return str(out_html)


def main():
    parser = argparse.ArgumentParser(description="Interactive HTML Pipeline Report Generator")
    parser.add_argument("-d", "--dir", required=True, help="Pipeline output root directory")
    parser.add_argument("-o", "--output", default=None, help="HTML report output path")
    parser.add_argument("--no_images", action="store_true", help="Skip image embedding (faster, smaller)")

    args = parser.parse_args()
    out = Path(args.dir)
    if not out.exists():
        sys.exit(f"Directory not found: {args.dir}")

    summary_in = None
    for p in ["high_conf.summary.tsv", "all_viruses.best.summary.tsv", "all_viruses.summary.tsv"]:
        fp = out / "1_FastViromeExplorer" / "summary" / p
        if fp.exists():
            summary_in = fp; break
    if summary_in is None:
        sys.exit("No summary file found")

    print("Collecting virus data...")
    viruses = collect_virus_data(args.dir, summary_in)
    print(f"  Found {len(viruses)} viruses")

    out_html = Path(args.output) if args.output else out / "Pipeline_Summary_Report.html"
    print("Generating interactive HTML...")
    generate_html(args.dir, out_html, viruses)
    print(f"  Report: {out_html} ({out_html.stat().st_size / 1024:.0f} KB)")


if __name__ == "__main__":
    main()
