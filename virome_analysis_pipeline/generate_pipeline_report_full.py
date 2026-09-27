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
from mmpv_common.io_layout import layout_dirs as _layout_dirs
_D = _layout_dirs(os.environ.get("MMPV_IO_LAYOUT", "legacy"))

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
    """Generate interactive HTML report with per-stage per-virus navigation."""
    now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    out = Path(out_dir)
    rpt_data = out / _D['d_reports'] / "report_data"
    copied = copy_stage_results(out_dir, rpt_data) if out else {}
    post_dir = out / "6_post_analysis"

    all_stages = [
        ("S1","Detection","stage_1","Salmon/Kallisto pseudo-alignment or traditional alignment. Poisson ratio false-positive filtering. Dual-track coverage filter.",False),
        ("S2","Filter","stage_2","Multi-dimensional filter: coverage, depth, reads, TPM, Poisson ratio, ANI, keyword search.",False),
        ("S3","Variants","stage_3","FreeBayes/iVar/LoFreq variant calling. Dynamic VCF filtering. SnpEff annotation + SNPGenie dN/dS.",True),
        ("S4","Post-hoc","stage_4","Per-virus post-hoc: mutation landscape, OncoPrint, dN/dS, PCA. Pipeline-level: sample distribution, co-abundance, depth, metadata.",True),
        ("S5","Assembly","stage_5","12-step de novo assembly: MEGAHIT/SPAdes to refineC to Divine Fusion to PVGA to polishing to gap-filling to circularization.",True),
        ("S6","Extract","stage_6","Longest contig extraction with N-fill via reference global pairwise alignment.",True),
        ("S7","Similarity","stage_7","Pairwise genome similarity: SDT-style NT/AA identity matrices, hierarchical clustering, CD-HIT deduplication.",True),
        ("S8","DVG","stage_8","ViReMa DVG detection: ViReMa to Circos 4-track recombination plots + arc diagrams.",True),
    ]
    virus_list = sorted(viruses.items())

    def _render_file(fp, fname):
        fn_low = fname.lower()
        if fn_low.endswith('.tsv') or fn_low.endswith('.csv'):
            try:
                rows = []
                with open(str(fp), 'r', encoding='utf-8', errors='replace') as tf:
                    for line in tf:
                        cols = [c.strip() for c in line.rstrip('\\n').split('\\t')]
                        rows.append(cols)
                if rows and len(rows) > 1:
                    tid = 'tbl_' + str(abs(hash(str(fp))))[:8]
                    all_rows = rows[1:]; total_rows = len(all_rows)
                    pp = 10; tp = (total_rows + pp - 1) // pp
                    th = ''.join('<th>' + c + '</th>' for c in rows[0])
                    h = '<div class="chart-card"><div class="chart-title">' + fname + ' (' + str(total_rows) + ' rows, ' + str(tp) + ' pages)</div>'
                    h += '<div class="pg-nav" id="' + tid + '_nav" style="display:flex;gap:8px;align-items:center;padding:6px 10px;background:#f8f9fa;border-bottom:1px solid #e0e0e0;font-size:12px">'
                    h += '<button onclick="pg_tbl(' + repr(tid) + ',-' + str(pp) + ')" style="padding:3px 10px;cursor:pointer;border:1px solid #ccc;border-radius:3px;background:#fff">&laquo; Prev</button>'
                    h += '<span id="' + tid + '_info" style="color:#666">Page 1 of ' + str(tp) + '</span>'
                    h += '<button onclick="pg_tbl(' + repr(tid) + ',' + str(pp) + ')" style="padding:3px 10px;cursor:pointer;border:1px solid #ccc;border-radius:3px;background:#fff">Next &raquo;</button>'
                    h += '</div><div class="tb-scroll"><table id="' + tid + '"><thead><tr>' + th + '</tr></thead><tbody>'
                    for ri, r in enumerate(all_rows):
                        cls = ' style="display:none"' if ri >= pp else ''
                        h += '<tr' + cls + '><td>' + '</td><td>'.join(r) + '</td></tr>'
                    h += '</tbody></table></div></div>'
                    return h
            except Exception: pass
            return '<div class="chart-card"><div class="chart-title">' + fname + '</div><div class="chart-placeholder">Cannot parse table</div></div>'
        elif fn_low.endswith('.pdf'):
            png_b64 = pdf_to_png_base64(str(fp), max_kb=1500)
            if png_b64:
                return '<div class="chart-card"><div class="chart-title">' + fname + '</div><img src="data:image/png;base64,' + png_b64 + '" loading="lazy"></div>'
            pdf_b64 = img_to_base64(str(fp), max_kb=800)
            if pdf_b64:
                return '<div class="chart-card"><div class="chart-title">' + fname + '</div><object data="data:application/pdf;base64,' + pdf_b64 + '" type="application/pdf" width="100%" height="500px"><p>PDF preview</p></object></div>'
            return '<div class="chart-card"><div class="chart-title">' + fname + '</div><div class="chart-placeholder">PDF too large</div></div>'
        elif fn_low.endswith(('.png','.jpg','.jpeg','.gif','.svg')):
            b64 = img_to_base64(str(fp), max_kb=1000)
            if b64:
                mime = "image/png"
                if fn_low.endswith('.jpg') or fn_low.endswith('.jpeg'): mime = "image/jpeg"
                elif fn_low.endswith('.svg'): mime = "image/svg+xml"
                elif fn_low.endswith('.gif'): mime = "image/gif"
                return '<div class="chart-card"><div class="chart-title">' + fname + '</div><img src="data:' + mime + ';base64,' + b64 + '" loading="lazy"></div>'
        return '<div class="chart-card"><div class="chart-title">' + fname + '</div><div class="chart-placeholder">File: ' + fname + '</div></div>'

    # ── Sidebar ──
    nav = '<li class="nav-item"><a href="#overview" class="nav-link active">Global Overview</a></li>'
    nav += '<li class="nav-header" style="padding:8px 20px 4px;font-size:11px;color:#3498db;font-weight:bold;border-top:1px solid #333;margin-top:4px">Pipeline Stages</li>'
    for sid, sname, sk, _, hs in all_stages:
        sa = "stage-" + sid.lower()
        nav += '<li class="nav-item"><a href="#' + sa + '" class="nav-link" style="font-size:11px;padding-left:24px;font-weight:bold">' + sid + ': ' + sname + '</a></li>'
        if hs:
            sn = int(sid[1:])
            for acc, data in virus_list:
                sp = data.get("species", acc)[:18]
                va = "s" + str(sn) + "-virus-" + acc
                nav += '<li class="nav-item"><a href="#' + va + '" class="nav-link" style="padding-left:40px;font-size:10px;color:#999">&bull; ' + sp + '</a></li>'
    nav += '<li class="nav-header" style="padding:8px 20px 4px;font-size:11px;color:#3498db;font-weight:bold;border-top:1px solid #333;margin-top:4px">Virus Summary</li>'
    for acc, data in virus_list:
        sp = data.get("species", acc)[:25]
        nav += '<li class="nav-item"><a href="#virus-' + acc + '" class="nav-link" style="padding-left:24px;font-size:10px;color:#999">  ' + sp + '</a></li>'

    # ── Body ──
    body = ""
    for sid, sname, sk, desc, hs in all_stages:
        sa = "stage-" + sid.lower()
        body += '<section id="' + sa + '"><h2>' + sid + ': ' + sname + '</h2>'
        body += '<p style="color:#555;font-size:13px;line-height:1.6;margin:6px 0">' + desc + '</p>'
        if sk in copied:
            body += '<h3 style="color:#666;font-size:14px;margin-top:12px">Stage-level Results</h3><div class="chart-gallery">'
            for fn in copied[sk]:
                fp = rpt_data / sk / fn
                if fp.is_dir(): continue
                body += _render_file(fp, fn)
            body += '</div>'
        if hs:
            sn = int(sid[1:])
            for acc, data in virus_list:
                sp = data.get("species", acc)
                va = "s" + str(sn) + "-virus-" + acc
                body += '<section id="' + va + '" style="margin-left:20px;border-left:3px solid #3498db;padding-left:15px">'
                body += '<h3 style="color:#2980b9;font-size:16px">' + sp + ' <span style="color:#888;font-size:12px">(' + acc + ')</span></h3>'
                ch = ""
                if sn == 4:
                    pc = collect_charts(post_dir, acc)
                    for ci in sorted(pc.values(), key=lambda x: x['label']):
                        if ci.get('base64'):
                            if ci['ext'] == 'pdf':
                                ch += '<div class="chart-card"><div class="chart-title">' + ci['label'] + '</div><object data="data:application/pdf;base64,' + ci['base64'] + '" type="application/pdf" width="100%" height="500px"><p>PDF</p></object></div>'
                            else:
                                ch += '<div class="chart-card"><div class="chart-title">' + ci['label'] + '</div><img src="data:image/png;base64,' + ci['base64'] + '" loading="lazy"></div>'
                elif sn == 8:
                    dc = collect_dvg_charts(out / "9_virema_dvg", acc)
                    for ci in sorted(dc.values(), key=lambda x: x['label']):
                        if ci.get('base64'):
                            if ci['ext'] == 'pdf':
                                ch += '<div class="chart-card"><div class="chart-title">' + ci['label'] + '</div><object data="data:application/pdf;base64,' + ci['base64'] + '" type="application/pdf" width="100%" height="500px"><p>PDF</p></object></div>'
                            else:
                                ch += '<div class="chart-card"><div class="chart-title">' + ci['label'] + '</div><img src="data:image/png;base64,' + ci['base64'] + '" loading="lazy"></div>'
                if ch:
                    body += '<div class="chart-gallery">' + ch + '</div>'
                else:
                    body += '<p style="color:#999;font-size:12px;margin:10px 0">No charts available for this virus at Stage ' + str(sn) + '.</p>'
                body += '</section>'
        body += '</section>'

    # Virus Summary
    vs = ""
    for acc, data in virus_list:
        sp = data.get("species", acc)
        cpm = str(data.get("cpm", "?")); cov = str(data.get("coverage", "?")); depth = str(data.get("depth", "?"))
        poisson = str(data.get("poisson", "?")); reads = str(data.get("reads", "?")); n = str(data.get("n_samples", "?"))
        pc = collect_charts(post_dir, acc); dc = collect_dvg_charts(out / "9_virema_dvg", acc); pc.update(dc)
        ch = ""
        for ci in sorted(pc.values(), key=lambda x: x['label']):
            if ci.get('base64'):
                if ci['ext'] == 'pdf':
                    ch += '<div class="chart-card"><div class="chart-title">' + ci['label'] + '</div><object data="data:application/pdf;base64,' + ci['base64'] + '" type="application/pdf" width="100%" height="500px"><p>PDF</p></object></div>'
                else:
                    ch += '<div class="chart-card"><div class="chart-title">' + ci['label'] + '</div><img src="data:image/png;base64,' + ci['base64'] + '" loading="lazy"></div>'
        vs += '<section id="virus-' + acc + '"><h2>' + sp + '</h2><p class="accession">' + acc + '</p>'
        vs += '<div class="metrics">'
        vs += '<div class="metric"><span class="value">' + n + '</span><span class="unit">Samples</span></div>'
        vs += '<div class="metric"><span class="value">' + cov + '%</span><span class="unit">Coverage</span></div>'
        vs += '<div class="metric"><span class="value">' + cpm + '</span><span class="unit">CPM</span></div>'
        vs += '<div class="metric"><span class="value">' + depth + 'x</span><span class="unit">Depth</span></div>'
        vs += '<div class="metric"><span class="value">' + poisson + '</span><span class="unit">Poisson</span></div>'
        vs += '<div class="metric"><span class="value">' + reads + '</span><span class="unit">Reads</span></div>'
        vs += '</div><div class="chart-gallery">' + (ch or "<p>No charts found.</p>") + '</div></section>'

    total_records = sum(v.get("n_samples", 0) for v in viruses.values())

    html = '''<!DOCTYPE html><html lang="en"><head><meta charset="UTF-8"><meta name="viewport" content="width=device-width, initial-scale=1.0"><title>Virus Pipeline Report</title>
<style>
*{{box-sizing:border-box;margin:0;padding:0}}
body{{font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',Arial,sans-serif;color:#333;display:flex;min-height:100vh}}
.sidebar{{position:fixed;left:0;top:0;width:260px;height:100vh;background:#1a1a2e;color:#e0e0e0;overflow-y:auto;padding:20px 0;z-index:100}}
.sidebar h3{{padding:0 20px 15px;color:#3498db;font-size:16px;border-bottom:1px solid #333;margin-bottom:10px}}
.nav-item{{list-style:none}}
.nav-link{{display:block;padding:6px 20px;color:#bbb;text-decoration:none;font-size:12px;transition:all 0.2s}}
.nav-link:hover,.nav-link.active{{color:#fff;background:#16213e}}
.main{{margin-left:260px;padding:30px 40px;flex:1;max-width:1400px}}
h1{{color:#2c3e50;font-size:28px;margin-bottom:5px;border-bottom:3px solid #3498db;padding-bottom:10px}}
h2{{color:#2980b9;font-size:22px;margin:30px 0 10px}}
h3{{color:#2c3e50;font-size:16px;margin:15px 0 8px}}
.accession{{color:#888;font-size:13px;margin-bottom:15px}}
.metrics{{display:flex;flex-wrap:wrap;gap:12px;margin:15px 0}}
.metric{{background:#f0f4f8;border-radius:8px;padding:12px 20px;text-align:center;min-width:80px}}
.metric .value{{display:block;font-size:22px;font-weight:bold;color:#2980b9}}
.metric .unit{{display:block;font-size:11px;color:#888;margin-top:2px}}
.chart-gallery{{display:grid;grid-template-columns:repeat(auto-fill,minmax(350px,1fr));gap:15px;margin:15px 0}}
.chart-card{{background:#fff;border:1px solid #e0e0e0;border-radius:8px;overflow:hidden;box-shadow:0 1px 3px rgba(0,0,0,.08)}}
.chart-title{{background:#f8f9fa;padding:8px 12px;font-size:13px;font-weight:600;color:#555;border-bottom:1px solid #e0e0e0}}
.chart-card img{{width:100%;height:auto;display:block}}
.chart-placeholder{{padding:40px;text-align:center;background:#fafafa}}
.chart-placeholder a{{color:#3498db;text-decoration:none;font-size:14px}}
.tb-scroll{{max-height:500px;overflow-y:auto;border:1px solid #eee;border-radius:4px}}
.tb-scroll table{{width:100%;border-collapse:collapse;font-size:11px}}
.tb-scroll th{{background:#1a1a2e;color:#fff;padding:6px 8px;text-align:left;position:sticky;top:0;z-index:1}}
.tb-scroll td{{padding:4px 8px;border-bottom:1px solid #eee;white-space:nowrap}}
.tb-scroll tr:nth-child(even){{background:#f8f9fa}}
.tb-scroll tr:hover{{background:#e8f0fe}}
.pg-nav button:hover{{background:#3498db;color:#fff}}
.card-row{{display:flex;flex-wrap:wrap;gap:12px;margin:15px 0}}
.card{{background:#f0f4f8;border-radius:8px;padding:15px 20px;text-align:center;min-width:120px;flex:1}}
.card .value{{display:block;font-size:26px;font-weight:bold;color:#2980b9}}
.card .label{{display:block;font-size:11px;color:#888;margin-top:3px}}
.footer{{margin-top:40px;padding:20px;text-align:center;color:#aaa;font-size:12px;border-top:1px solid #eee}}
section{{margin-bottom:30px;padding-top:10px}}
@media(max-width:800px){{.sidebar{{display:none}}.main{{margin-left:0}}}}
</style></head><body>
<nav class="sidebar"><h3>Pipeline Report</h3><ul>''' + nav + '''</ul></nav>
<main class="main">
<h1>Known Virus Pipeline Report</h1>
<p style="color:#888;font-size:13px">Generated: ''' + now + ''' | Directory: ''' + str(out_dir) + '''</p>
<section id="overview"><h2>Global Overview</h2>
<div class="card-row">
<div class="card"><div class="value">''' + str(len(viruses)) + '''</div><div class="label">Virus Species</div></div>
<div class="card"><div class="value">''' + str(total_records) + '''</div><div class="label">Records</div></div>
</div></section>
''' + body + '''
''' + vs + '''
<section id="paper-reference"><h2>Paper Reference</h2>
<div class="paper-block">
<p>The Known Virus Analysis Pipeline processed <strong>''' + str(len(viruses)) + ''' virus species</strong> across samples.</p>
<table><tr><th>Stage</th><th>Script</th><th>Output</th></tr>
<tr><td>1. Detect</td><td>batch_virus_depth.py</td><td>Salmon pseudo-alignment + Poisson filter</td></tr>
<tr><td>2. Filter</td><td>utils/filter_summary.py</td><td>Coverage &ge;50%, Depth &ge;5x, Reads &ge;100</td></tr>
<tr><td>3. Variants</td><td>batch_virus_variants.py</td><td>FreeBayes + SnpEff + SNPGenie</td></tr>
<tr><td>4. Post-hoc</td><td>6-script suite</td><td>VCF viz + SnpEff macro + MAF + SnpGenie</td></tr>
<tr><td>5. Assembly</td><td>virus-full.py</td><td>12-step de novo assembly</td></tr>
<tr><td>6. Extract</td><td>utils/extract_full_fasta.py</td><td>Longest contig extraction</td></tr>
<tr><td>7. Similarity</td><td>virus_auto_pipeline.py</td><td>SDT pairwise similarity</td></tr>
<tr><td>8. DVG</td><td>batch_virema_dvg.py</td><td>ViReMa recombination + Circos</td></tr>
<tr><td>9. Report</td><td>generate_pipeline_report.py</td><td>This interactive report</td></tr>
</table></div></section>
<p style="color:#888;font-size:11px;margin-top:10px">Report data: 10_Reports/report_data/</p>
<div class="footer">Generated by known_virus_pipeline &mdash; ''' + now + '''</div>
</main>
<script>
function pg_tbl(tid, offset) {{
    var tbl = document.getElementById(tid);
    var rows = tbl.getElementsByTagName('tbody')[0].getElementsByTagName('tr');
    var per = Math.abs(offset);
    var cur = 0;
    for (var i = 0; i < rows.length; i++) {{
        if (rows[i].style.display !== 'none') {{ cur = Math.floor(i / per) * per; break; }}
    }}
    var next = cur + offset;
    if (next < 0 || next >= rows.length) return;
    for (var i = 0; i < rows.length; i++) {{
        rows[i].style.display = (i >= next && i < next + per) ? '' : 'none';
    }}
    var page = Math.floor(next / per) + 1;
    var total = Math.ceil(rows.length / per);
    document.getElementById(tid + '_info').textContent = 'Page ' + page + ' of ' + total;
}}
</script>
</body></html>'''

    with open(out_html, "w", encoding="utf-8") as f:
        f.write(html)
    return str(out_html)

