#!/usr/bin/env python3
"""
generate_report.py — Public Metadata Pipeline HTML Report Generator
===================================================================
Reads search output + info output + plot output, generates interactive HTML.

Usage:
  python generate_report.py -d public_data_pipeline_output/
  python generate_report.py -d public_data_pipeline_output/ --ai-key sk-xxx --ai-provider deepseek
"""

import argparse, os, sys, base64, csv, json as _json, re, urllib.request
from pathlib import Path
from datetime import datetime
from collections import Counter

# 跨管线统一 I/O 布局 (mmpv_common/, 仓库根): 目录名随 MMPV_IO_LAYOUT 解析
# (编排器已 normalize 环境变量, 子进程导入本模块时快照即正确布局)
_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)
from mmpv_common.io_layout import layout_dirs as _layout_dirs
_D = _layout_dirs(os.environ.get("MMPV_IO_LAYOUT", "legacy"))
_HOSTREF_ROOT = str(Path(_D['h_genome']).parent)   # legacy=host_reference / standard=05_HostRef

try:
    import pandas as pd
    HAS_PANDAS = True
except ImportError:
    HAS_PANDAS = False

SCRIPT_DIR = Path(__file__).resolve().parent

if sys.platform == 'win32':
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')

def safe_read_csv(fp, sep=","):
    """Read CSV/TSV, auto-detect separator."""
    p = Path(fp)
    if not p.is_file(): return None
    try:
        if not HAS_PANDAS: return None
        # Try comma first, then tab
        try:
            df = pd.read_csv(p, sep=sep, nrows=5)
            if len(df.columns) <= 1 and sep == ",":
                df = pd.read_csv(p, sep="\t")
        except:
            if sep == ",":
                df = pd.read_csv(p, sep="\t")
            else:
                return None
        return pd.read_csv(p, sep="\t" if sep != "," else ",") if sep != "," else df
    except Exception:
        return None

def img_to_base64(path, max_kb=2000):
    if not path or not Path(path).exists(): return None
    if Path(path).stat().st_size / 1024 > max_kb: return None
    with open(path, 'rb') as f:
        return base64.b64encode(f.read()).decode('utf-8')

def _esc(s):
    return str(s).replace("&","&amp;").replace("<","&lt;").replace(">","&gt;").replace('"',"&quot;")

# ═══════════════════════════════════════════
# Data Collection
# ═══════════════════════════════════════════

def collect_search_summary(out_dir):
    """Read SRA_GSA_Merged_Final.csv and return summary stats."""
    out = Path(out_dir)
    csv_path = None
    for p in out.rglob("SRA_GSA_Merged_Final.csv"):
        csv_path = p; break
    if not csv_path:
        for p in out.rglob("*search*/*Merged*.csv"):
            csv_path = p; break
    if not csv_path: return None
    df = safe_read_csv(csv_path)
    if df is None: return None
    total = len(df)
    sra_n = int((df["Database"] == "SRA").sum()) if "Database" in df.columns else 0
    gsa_n = int((df["Database"] == "GSA").sum()) if "Database" in df.columns else 0
    n_bioprojects = df["BioProject"].nunique() if "BioProject" in df.columns else "?"
    n_species = df["ScientificName"].nunique() if "ScientificName" in df.columns else "?"
    # Temporal range
    years = "?"
    if "ReleaseDate" in df.columns:
        try:
            dates = pd.to_datetime(df["ReleaseDate"], errors="coerce")
            ymin = int(dates.dt.year.min()) if dates.notna().any() else "?"
            ymax = int(dates.dt.year.max()) if dates.notna().any() else "?"
            years = f"{ymin}-{ymax}"
        except: pass
    return {"total": total, "sra": sra_n, "gsa": gsa_n, "bioprojects": n_bioprojects,
            "species": n_species, "years": years, "path": str(csv_path.relative_to(out)) if csv_path else "?"}

def collect_info_summary(out_dir):
    """Read Core14 metadata and return field completeness + distributions."""
    out = Path(out_dir)
    csv_path = None
    for p in out.rglob("Global_Unified_Metadata_Core14.csv"):
        csv_path = p; break
    if not csv_path:
        for p in out.rglob("*info*/**/Core14*.csv"):
            csv_path = p; break
    if not csv_path: return None
    df = safe_read_csv(csv_path)
    if df is None: return None
    total = len(df)
    fields = {}
    for col in df.columns:
        non_null = int(df[col].notna().sum())
        non_empty = int((df[col].astype(str).str.strip() != "").sum())
        pct = f"{non_empty/total*100:.0f}%" if total > 0 else "?"
        fields[col] = {"complete": non_empty, "pct": pct}
    # Top distributions
    dist = {}
    for col in ["Source", "Tissue", "Location", "LibrarySource", "CenterName"]:
        if col not in df.columns: continue
        counts = df[col].astype(str).value_counts().head(8).to_dict()
        dist[col] = [(k, v) for k, v in counts.items() if k and k != "nan"]
    # PMID coverage
    pmid_n = int((df["PMID"].notna() & (df["PMID"].astype(str).str.strip() != "")).sum()) if "PMID" in df.columns else 0
    return {"total": total, "fields": fields, "distributions": dist, "pmid_coverage": pmid_n,
            "path": str(csv_path.relative_to(out)) if csv_path else "?"}

def collect_plot_images(out_dir):
    """Find SCI landscape plot images."""
    out = Path(out_dir)
    images = {}
    for name, pat in [("landscape", "Combined_Landscape_Full.png"), ("landscape_pdf", "Combined_Landscape_Full.pdf")]:
        for p in out.rglob(pat):
            images[name] = str(p)
            break
    return images

def collect_download_summary(out_dir):
    """download_report_*.csv 状态统计 + failed_sra_*.txt 行数 + down/ 目录文件盘点."""
    out = Path(out_dir)
    status, methods, n_rows = {}, {}, 0
    for p in out.rglob("download_report_*.csv"):
        try:
            with open(p, encoding="utf-8-sig") as f:
                for r in csv.DictReader(f):
                    n_rows += 1
                    st = (r.get("Status") or "?").strip()
                    me = (r.get("Method") or "?").strip()
                    status[st] = status.get(st, 0) + 1
                    methods[me] = methods.get(me, 0) + 1
        except Exception:
            continue
    failed, failed_files = 0, []
    for p in out.rglob("failed_sra_*.txt"):
        try:
            n = sum(1 for line in p.read_text(encoding="utf-8").splitlines() if line.strip())
        except Exception:
            n = 0
        failed += n
        failed_files.append(p.name)
    fq = [f for f in out.rglob("*.fastq.gz") if f.is_file()]
    fq_gb = sum(f.stat().st_size for f in fq) / 1024 ** 3
    sra_left = sum(1 for f in out.rglob("*.sra") if f.is_file())
    return {"rows": n_rows, "status": status, "methods": methods, "failed": failed,
            "failed_files": failed_files[:5], "fastq_count": len(fq),
            "fastq_gb": fq_gb, "sra_left": sra_left}

RUN_PATTERNS = [
    re.compile(r"^(?P<acc>.+?)_[12]\.(?:fastq|fq)\.gz$"),
    re.compile(r"^(?P<acc>.+)\.(?:fastq|fq)\.gz$"),
]

def write_handoff_manifest(out_dir):
    """扫描 down/ 的 FASTQ.GZ, 对齐 Core14 元数据, 产出下游交接清单 sample_handoff.csv.

    交棒给 data_preprocessing_pipeline: R1/R2_Path 即 clean-data.py --input 的文件级清单。
    """
    out = Path(out_dir)
    files = {}
    for f in sorted(out.rglob("*.fastq.gz")) + sorted(out.rglob("*.fq.gz")):
        if not f.is_file():
            continue
        acc = None
        for pat in RUN_PATTERNS:
            m = pat.match(f.name)
            if m:
                acc = m.group("acc")
                break
        if not acc:
            continue
        e = files.setdefault(acc, {})
        if re.search(r"_1\.(?:fastq|fq)\.gz$", f.name):
            e["r1"] = f
        elif re.search(r"_2\.(?:fastq|fq)\.gz$", f.name):
            e["r2"] = f
        else:
            e.setdefault("se", f)
    meta_rows = {}
    core_csv = next(iter(out.rglob("Global_Unified_Metadata_Core14.csv")), None)
    if core_csv:
        try:
            with open(core_csv, encoding="utf-8-sig") as f:
                for r in csv.DictReader(f):
                    run = (r.get("Run") or "").strip()
                    if run:
                        meta_rows[run] = r
        except Exception:
            pass
    accs = sorted(set(files) | set(meta_rows))
    handoff = out / "sample_handoff.csv"
    with open(handoff, "w", encoding="utf-8-sig", newline="") as f:
        w = csv.writer(f)
        w.writerow(["Run", "Layout", "R1_Path", "R2_Path", "SE_Path", "Size_MB",
                    "ScientificName", "Tissue", "Location", "BioProject", "Status"])
        n_ready = 0
        for acc in accs:
            e = files.get(acc, {})
            r1, r2, se = e.get("r1"), e.get("r2"), e.get("se")
            size = sum(x.stat().st_size for x in (r1, r2, se) if x is not None) / 1024 ** 2
            layout = "PE" if (r1 or r2) else ("SE" if se else "?")
            m = meta_rows.get(acc, {})
            ready = bool(r1 or r2 or se)
            if ready:
                n_ready += 1
            def _s(p):
                return str(p.resolve()) if p is not None else ""
            w.writerow([acc, layout, _s(r1), _s(r2), _s(se), f"{size:.1f}",
                        m.get("ScientificName", ""), m.get("Tissue", ""), m.get("Location", ""),
                        m.get("BioProject", ""), "ready" if ready else "missing"])
    return {"path": str(handoff), "runs": len(accs), "ready": n_ready,
            "missing": len(accs) - n_ready}

def collect_hostref_summary(out_dir):
    """扫描 hostref 阶段产物 (work_dir/host_reference/): 基因组 FASTA + 四索引就绪状态."""
    out = Path(out_dir)
    hr = None
    for cand in [out / _HOSTREF_ROOT, out / "build_host_pipeline_output"]:
        if cand.is_dir():
            hr = cand
            break
    if not hr:
        for p in out.rglob("hostdb"):
            if p.is_dir():
                hr = p.parent
                break
    if not hr:
        return None
    genome = None
    for pat in ["all.genome.uniq.fasta", "genome/*.fasta", "genome/*.fna", "*.fasta", "*.fna"]:
        genome = next(iter(hr.rglob(pat)), None)
        if genome:
            break
    hostdb = hr / "hostdb" if (hr / "hostdb").is_dir() else hr

    def _has(*pats):
        return any(hostdb.rglob(p) for p in pats)

    idx = {
        "Kraken2": _has("*.k2d"),
        "Bowtie2": _has("*.1.bt2", "*.1.bt2l"),
        "HISAT2": _has("*.1.ht2"),
        "Minimap2": _has("*.mmi"),
    }
    return {"root": str(hr), "genome": str(genome) if genome else None,
            "genome_mb": round(genome.stat().st_size / 1024 ** 2, 1) if genome else None,
            "indexes": idx, "built": sum(idx.values())}

# ═══════════════════════════════════════════
# AI Summary
# ═══════════════════════════════════════════

def generate_ai_summary(search_data, info_data, plot_images):
    provider = getattr(generate_ai_summary, '_provider', 'openai')
    model = getattr(generate_ai_summary, '_model', 'gpt-4o-mini')
    api_key = getattr(generate_ai_summary, '_api_key', '')
    base_url = getattr(generate_ai_summary, '_base_url', '')
    if not api_key: return ""

    s = search_data or {}
    i = info_data or {}
    prompt = f"""你是生物信息学研究员。基于以下公共数据检索管线结果撰写中文研究简报。

## 检索结果
- 检索到 {s.get('total','?')} 条记录: SRA={s.get('sra','?')}, GSA={s.get('gsa','?')}
- 涉及 {s.get('bioprojects','?')} 个 BioProject, {s.get('species','?')} 个物种
- 时间跨度: {s.get('years','?')}

## 元数据质量
- Core14 统一元数据: {i.get('total','?')} 条记录
- PMID 覆盖率: {i.get('pmid_coverage','?')} 条有文献链接

## 元数据字段完整度
{_json.dumps({k: v.get('pct','?') for k,v in (i.get('fields',{}) or {}).items()}, ensure_ascii=False)}

输出格式:
[BRIEF]
150-200字中文简报: 研究目的、检索规模、关键发现（数据库分布、元数据质量、文献覆盖）

[DETAILED]
按 [BACKGROUND][METHODS][RESULTS][DISCUSSION] 格式输出 300-400 字中文详报。

仅陈述数据支持的事实。中文。学术书面语。不捏造。"""

    try:
        url = base_url or "https://api.openai.com/v1/chat/completions"
        body = _json.dumps({
            "model": model, "messages": [
                {"role": "system", "content": "你是生物信息学研究员。仅基于数据陈述。中文输出。"},
                {"role": "user", "content": prompt}
            ], "temperature": 0.3, "max_tokens": 1500
        }).encode()
        req = urllib.request.Request(url, data=body, headers={
            "Content-Type": "application/json", "Authorization": f"Bearer {api_key}"
        })
        with urllib.request.urlopen(req, timeout=300) as resp:
            raw = _json.loads(resp.read())["choices"][0]["message"]["content"].strip()
        brief_m = re.search(r'\[BRIEF\]\s*(.*?)(?=\[DETAILED\]|\Z)', raw, re.DOTALL | re.IGNORECASE)
        brief = brief_m.group(1).strip() if brief_m else raw[:500]
        det_m = re.search(r'\[DETAILED\]\s*(.*)', raw, re.DOTALL | re.IGNORECASE)
        detailed = det_m.group(1).strip() if det_m else ""
        def _sec(t, tag):
            m = re.search(rf'\[{tag}\]\s*(.*?)(?=\[(?:BACKGROUND|METHODS|RESULTS|DISCUSSION)\]|\Z)', t, re.DOTALL | re.IGNORECASE)
            return m.group(1).strip() if m else ""
        secs = {}
        if detailed:
            for k in ["BACKGROUND","METHODS","RESULTS","DISCUSSION"]:
                secs[k] = _sec(detailed, k)
        if not any(secs.values()): secs["RESULTS"] = detailed or raw
        labels = {"BACKGROUND":"Background","METHODS":"Methods","RESULTS":"Results","DISCUSSION":"Discussion"}
        detail_items = ""
        for k in ["BACKGROUND","METHODS","RESULTS","DISCUSSION"]:
            t = secs.get(k, "")
            if not t: continue
            detail_items += f'<div class="ai-block"><div class="ai-block-title">{labels[k]}</div><div class="ai-block-text">{_esc(t)}</div></div>'
        return f'''<div class="ai-section" id="ai-summary">
<div class="ai-header"><span class="ai-title">AI Research Brief</span><span class="ai-model-tag">{_esc(model)}</span></div>
<div class="ai-brief">{_esc(brief)}</div>
<details class="ai-detail-toggle"><summary>IMRaD Detailed Report</summary><div class="ai-detail-content">{detail_items}</div></details>
</div>'''
    except Exception as e:
        print(f"  [WARN] AI Summary failed: {e}")
        return f'<div class="ai-error">AI summary failed: {e}</div>'

# ═══════════════════════════════════════════
# HTML Generation
# ═══════════════════════════════════════════

def generate_html(out_dir, search_data, info_data, plot_images, ai_html="", out_html=None,
                  download_data=None, handoff_data=None, hostref_data=None):
    out = Path(out_dir); now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    s = search_data or {}; i = info_data or {}; d = download_data or {}
    h = handoff_data or {}; hr = hostref_data or {}

    # Build stage data JSON for browser AI
    sd = {
        "_overview": {"records": s.get("total","?"), "sra": s.get("sra","?"), "gsa": s.get("gsa","?"),
                      "bioprojects": s.get("bioprojects","?"), "species": s.get("species","?"), "years": s.get("years","?")},
        "s1_search": {"total": s.get("total","?"), "sra": s.get("sra","?"), "gsa": s.get("gsa","?"),
                      "bioprojects": s.get("bioprojects","?"), "species": s.get("species","?"), "years": s.get("years","?")},
        "s2_info": {"total": i.get("total","?"), "pmid_coverage": i.get("pmid_coverage","?"),
                    "fields": {k: v.get("pct","?") for k,v in (i.get("fields",{}) or {}).items()}},
        "s3_down": {"report_rows": d.get("rows","?"), "status": d.get("status",{}),
                    "failed": d.get("failed","?"), "fastq_files": d.get("fastq_count","?"),
                    "fastq_gb": round(d.get("fastq_gb",0),1), "sra_left": d.get("sra_left","?")},
        "s4_plot": {"note": "6-panel SCI landscape figure: temporal trends, database composition, top organizations, tissues, locations, growth stages"},
        "s5_hostref": {"genome": hr.get("genome"), "genome_mb": hr.get("genome_mb"),
                       "indexes": hr.get("indexes", {}), "built": hr.get("built", "?")},
        "s6_handoff": {"runs": h.get("runs","?"), "ready": h.get("ready","?"),
                       "missing": h.get("missing","?"), "manifest": h.get("path","sample_handoff.csv")},
    }
    stage_data_json = _json.dumps(sd, ensure_ascii=False)

    # ── Sidebar ──
    nav_items = '<li class="nav-item"><a href="#overview" class="nav-link active">Overview</a></li>'
    nav_items += '<li class="nav-header">Pipeline Stages</li>'
    nav_items += '<div class="nav-sub">'
    for sid, label in [("stage-search","S1: Search"),("stage-info","S2: Info"),("stage-down","S3: Download"),("stage-plot","S4: Plot"),("stage-hostref","S5: HostRef"),("stage-handoff","S6: Handoff")]:
        nav_items += f'<li class="nav-item"><a href="#{sid}" class="nav-link stage-toggle" style="font-size:12px;padding-left:24px;font-weight:600">{label}</a></li>'
    nav_items += '</div>'
    nav_items += '<li class="nav-header">Report</li><div class="nav-sub">'
    nav_items += '<li class="nav-item"><a href="#ai-summary-area" class="nav-link" style="font-size:12px;padding-left:24px">AI Summary</a></li>'
    nav_items += '</div>'

    # ── Global Overview ──
    overview = f'''<section id="overview"><h2>Public Metadata Pipeline Report</h2>
    <p style="color:var(--ink-secondary);font-size:13px;line-height:1.6">Dual-engine (NCBI SRA + CNCB GSA) species-level metadata search, unification, download, and visualization.
    Generated: {now}</p>
    <div class="card-row"><div class="card"><div class="value">{s.get("total","?")}</div><div class="label">Total Records</div></div>
    <div class="card"><div class="value">{s.get("sra","?")}</div><div class="label">SRA</div></div>
    <div class="card"><div class="value">{s.get("gsa","?")}</div><div class="label">GSA</div></div>
    <div class="card"><div class="value">{s.get("bioprojects","?")}</div><div class="label">BioProjects</div></div>
    <div class="card"><div class="value">{s.get("species","?")}</div><div class="label">Species</div></div>
    <div class="card"><div class="value">{s.get("years","?")}</div><div class="label">Year Range</div></div></div></section>'''

    # ── Stage 1: Search ──
    s1 = f'''<section id="stage-search"><h2>S1: Search <button onclick="runStageAI('s1_search')" class="ai-stage-btn">AI</button></h2>
    <p style="color:var(--ink-secondary);font-size:13px">Dual-engine search across NCBI SRA (E-utilities API) and CNCB GSA (web scraping). Species Latin name + TaxID query.</p>'''
    if s:
        s1 += f'<div class="card-row"><div class="card"><div class="value">{s.get("total","?")}</div><div class="label">Merged Records</div></div>'
        s1 += f'<div class="card"><div class="value">{s.get("sra","?")}</div><div class="label">SRA</div></div>'
        s1 += f'<div class="card"><div class="value">{s.get("gsa","?")}</div><div class="label">GSA</div></div></div>'
        s1 += f'<p style="font-size:12px;color:var(--ink-secondary)">Source: {_esc(s.get("path","?"))}</p>'
    else:
        s1 += '<p style="color:var(--ink-secondary)">Search output not found. Run Stage 1 (search) first.</p>'
    s1 += '</section>'

    # ── Stage 2: Info ──
    s2 = f'''<section id="stage-info"><h2>S2: Info <button onclick="runStageAI('s2_info')" class="ai-stage-btn">AI</button></h2>
    <p style="color:var(--ink-secondary);font-size:13px">Deep metadata unification: SRA XML parsing + GSA web scraping + AI inference (DeepSeek/Kimi) for 13 core fields. BioProject → PubMed tracing.</p>'''
    if i:
        s2 += f'<div class="card-row"><div class="card"><div class="value">{i.get("total","?")}</div><div class="label">Unified Records</div></div>'
        s2 += f'<div class="card"><div class="value">{i.get("pmid_coverage","?")}</div><div class="label">With PubMed</div></div></div>'
        # Field completeness table
        fields = i.get("fields", {})
        if fields:
            s2 += '<details style="margin-top:8px"><summary style="font-weight:600;color:var(--accent);cursor:pointer;font-size:13px">Field Completeness</summary>'
            s2 += '<div class="tb-scroll"><table style="font-size:12px"><tr><th>Field</th><th>Complete</th><th>Rate</th></tr>'
            for fn, fv in sorted(fields.items()):
                s2 += f'<tr><td>{_esc(fn)}</td><td>{fv.get("complete","?")}</td><td>{fv.get("pct","?")}</td></tr>'
            s2 += '</table></div></details>'
        # Distribution tables
        dists = i.get("distributions", {})
        for col, items in dists.items():
            if not items: continue
            s2 += f'<details style="margin-top:6px"><summary style="font-weight:600;color:var(--accent);cursor:pointer;font-size:12px">Top {_esc(col)}</summary>'
            s2 += '<div class="tb-scroll"><table style="font-size:12px"><tr><th>Value</th><th>Count</th></tr>'
            for val, cnt in items[:10]:
                s2 += f'<tr><td>{_esc(val[:50])}</td><td>{cnt}</td></tr>'
            s2 += '</table></div></details>'
        s2 += f'<p style="font-size:12px;color:var(--ink-secondary);margin-top:4px">Source: {_esc(i.get("path","?"))}</p>'
    else:
        s2 += '<p style="color:var(--ink-secondary)">Info output not found. Run Stage 2 (info) first.</p>'
    s2 += '</section>'

    # ── Stage 3: Download ──
    s3 = f'''<section id="stage-down"><h2>S3: Download <button onclick="runStageAI('s3_down')" class="ai-stage-btn">AI</button></h2>
    <p style="color:var(--ink-secondary);font-size:13px">Smart SRA downloader: dual-protocol (aria2c NGDC + prefetch NCBI), progress tracking, failed retry list.</p>'''
    if d.get("rows") or d.get("fastq_count"):
        chips = "".join(f'<div class="card"><div class="value">{v}</div><div class="label">Status: {_esc(k)}</div></div>'
                        for k, v in sorted(d.get("status", {}).items()))
        s3 += f'<div class="card-row"><div class="card"><div class="value">{d.get("rows","?")}</div><div class="label">Report Rows</div></div>{chips}'
        s3 += f'<div class="card"><div class="value">{d.get("failed","?")}</div><div class="label">Failed (retry list)</div></div>'
        s3 += f'<div class="card"><div class="value">{d.get("fastq_count","?")}</div><div class="label">FASTQ.GZ Files</div></div>'
        s3 += f'<div class="card"><div class="value">{round(d.get("fastq_gb",0),1)}</div><div class="label">Data Size (GB)</div></div>'
        if d.get("sra_left"):
            s3 += f'<div class="card"><div class="value" style="color:#d97706">{d["sra_left"]}</div><div class="label">.SRA Not Converted</div></div>'
        s3 += '</div>'
        if d.get("methods"):
            s3 += f'<p style="font-size:12px;color:var(--ink-secondary)">Methods: {_esc(", ".join(f"{k}×{v}" for k, v in sorted(d["methods"].items())))}</p>'
        if d.get("failed_files"):
            s3 += f'<p style="font-size:12px;color:var(--ink-secondary)">Retry lists: {_esc(", ".join(d["failed_files"]))}</p>'
    else:
        s3 += '<p style="color:var(--ink-secondary)">Download report not found. Run Stage 3 (down) first.</p>'
    s3 += '</section>'

    # ── Stage 4: Plot ──
    s4 = f'''<section id="stage-plot"><h2>S4: Plot <button onclick="runStageAI('s4_plot')" class="ai-stage-btn">AI</button></h2>
    <p style="color:var(--ink-secondary);font-size:13px">SCI-grade 6-panel landscape figure: (A) Temporal Distribution, (B) Database Proportion, (C) Top Organizations, (D) Top Tissues, (E) Top Locations, (F) Top Growth Stages.</p>'''
    if plot_images.get("landscape"):
        b64 = img_to_base64(plot_images["landscape"], 5000)
        if b64:
            s4 += f'<div class="chart-card"><div class="chart-title">Combined Landscape (6-panel SCI Figure)</div><img src="data:image/png;base64,{b64}" loading="lazy" alt="SCI Landscape Figure"></div>'
    else:
        s4 += '<p style="color:var(--ink-secondary)">Plot not found. Run Stage 4 (plot) to generate the 6-panel figure.</p>'
    s4 += '</section>'

    # ── Stage 5: HostRef ──
    s_hr = f'''<section id="stage-hostref"><h2>S5: HostRef 宿主参考基因组</h2>
    <p style="color:var(--ink-secondary);font-size:13px">Host genome acquisition (datasets / ncbi-genome-download / FTP / gget, or user FASTA)
    and Kraken2 / Bowtie2 / HISAT2 / Minimap2 index building — consumed by the host-depletion stage downstream.</p>'''
    if hr:
        idx = hr.get("indexes", {})
        genome_mb = hr.get("genome_mb")
        genome_badge = f" ({genome_mb} MB)" if genome_mb else ""
        genome_color = "#16a34a" if hr.get("genome") else "#d97706"
        s_hr += f'<div class="card-row"><div class="card"><div class="value">{hr.get("built","?")}/4</div><div class="label">Indexes Built</div></div>'
        s_hr += f'<div class="card"><div class="value" style="font-size:14px;color:{genome_color}">{"OK" if hr.get("genome") else "N/A"}</div><div class="label">Genome FASTA{genome_badge}</div></div>'
        for name, ok in idx.items():
            color = "#16a34a" if ok else "#d97706"
            mark = "✓" if ok else "✗"
            s_hr += f'<div class="card"><div class="value" style="font-size:14px;color:{color}">{mark}</div><div class="label">{_esc(name)}</div></div>'
        s_hr += '</div>'
        s_hr += f'<p style="font-size:12px;color:var(--ink-secondary)">Root: <code>{_esc(hr.get("root",""))}</code></p>'
    else:
        s_hr += '<p style="color:var(--ink-secondary)">Host reference not found under <code>host_reference/</code>. Run hostref + hostdb stages first, or pass --host-fasta.</p>'
    s_hr += '</section>'

    # ── Stage 6: Handoff ──
    s5 = f'''<section id="stage-handoff"><h2>S6: Handoff 交接清单</h2>
    <p style="color:var(--ink-secondary);font-size:13px">Downstream handoff manifest: Run → FASTQ absolute paths + Core14 metadata.
    Consumed by <code>data_preprocessing_pipeline/clean-data.py --input</code> as the file-level sample list.</p>'''
    if h:
        s5 += f'<div class="card-row"><div class="card"><div class="value">{h.get("runs","?")}</div><div class="label">Runs Total</div></div>'
        s5 += f'<div class="card"><div class="value" style="color:#16a34a">{h.get("ready","?")}</div><div class="label">Ready (files found)</div></div>'
        s5 += f'<div class="card"><div class="value" style="color:#d97706">{h.get("missing","?")}</div><div class="label">Missing (metadata only)</div></div></div>'
        s5 += f'<p style="font-size:12px;color:var(--ink-secondary)">Manifest: <code>{_esc(h.get("path","sample_handoff.csv"))}</code> — 列: Run / Layout / R1_Path / R2_Path / SE_Path / Size_MB / ScientificName / Tissue / Location / BioProject / Status</p>'
    else:
        s5 += '<p style="color:var(--ink-secondary)">Handoff manifest not generated (no down/ files and no Core14 metadata found).</p>'
    s5 += '</section>'

    # ── AI Summary area ──
    ai_area = f'<section id="ai-summary-area">{ai_html}</section>' if ai_html else ''

    html = f'''<!DOCTYPE html><html lang="en"><head><meta charset="UTF-8"><meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>Public Metadata Pipeline Report</title>
<style>
:root {{
  --ink: #1a2332; --ink-secondary: #5a6a7e; --accent: #2563eb; --accent-hover: #1d4ed8;
  --accent-subtle: #eff4ff; --surface: #ffffff; --surface-alt: #f6f8fb; --surface-hover: #eef1f6;
  --border: #e1e5eb; --border-light: #eef0f4; --sidebar-bg: #f3f5f8; --sidebar-ink: #4a5568;
  --radius: 6px; --shadow-sm: 0 1px 2px rgba(0,0,0,.04); --shadow-md: 0 4px 12px rgba(0,0,0,.06);
}}
* {{box-sizing:border-box;margin:0;padding:0}}
body {{font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',system-ui,sans-serif;color:var(--ink);background:var(--surface-alt);display:flex;min-height:100vh;-webkit-font-smoothing:antialiased}}
.sidebar {{position:fixed;left:0;top:0;width:248px;height:100vh;background:var(--sidebar-bg);color:var(--sidebar-ink);overflow-y:auto;padding:20px 0;z-index:100;border-right:1px solid var(--border)}}
.sidebar h3 {{padding:0 20px 12px;color:var(--ink);font-size:16px;font-weight:700;border-bottom:1px solid var(--border);margin-bottom:8px}}
.nav-item {{list-style:none}}
.nav-link {{display:block;padding:7px 20px;color:var(--sidebar-ink);text-decoration:none;font-size:13px;transition:background 0.15s,color 0.15s;border-radius:0 20px 20px 0;margin-right:8px}}
.nav-link:hover {{background:var(--surface-hover);color:var(--ink)}}
.nav-link.active {{background:var(--accent-subtle);color:var(--accent);font-weight:600}}
.nav-header {{cursor:pointer;user-select:none;display:flex;justify-content:space-between;align-items:center;padding:8px 20px 4px;font-size:12px;color:var(--ink-secondary);font-weight:600;border-top:1px solid var(--border);margin-top:4px}}
.nav-header::after {{content:'\\25B2';font-size:8px;transition:transform 0.15s}}
.nav-header.collapsed::after {{transform:rotate(180deg)}}
.nav-sub {{overflow:hidden;max-height:2000px;transition:max-height 0.25s ease}}
.nav-sub.collapsed {{max-height:0}}
.main {{margin-left:248px;padding:32px 40px;flex:1;max-width:1440px}}
h1 {{color:var(--ink);font-size:26px;font-weight:700;margin-bottom:4px;padding-bottom:10px;border-bottom:2px solid var(--border)}}
h2 {{color:var(--ink);font-size:19px;font-weight:600;margin:28px 0 8px}}
section {{margin-bottom:22px;padding-top:8px;scroll-margin-top:20px}}
.card-row {{display:flex;flex-wrap:wrap;gap:10px;margin:12px 0}}
.card {{background:var(--surface);border:1px solid var(--border);border-radius:var(--radius);padding:14px 18px;text-align:center;min-width:100px;flex:1;transition:box-shadow 0.15s}}
.card:hover {{box-shadow:var(--shadow-sm)}}
.card .value {{display:block;font-size:24px;font-weight:700;color:var(--accent);word-break:break-all}}
.card .label {{display:block;font-size:10px;color:var(--ink-secondary);margin-top:2px;font-weight:500}}
table {{border-collapse:collapse;width:100%;font-size:12px}}
th,td {{padding:7px 10px;text-align:left;border-bottom:1px solid var(--border-light)}}
th {{background:var(--surface-alt);font-size:11px;color:var(--ink-secondary);font-weight:600;position:sticky;top:0;z-index:1}}
tbody tr:nth-child(even) td {{background:#f9fafb}}
tbody tr:hover td {{background:var(--accent-subtle)}}
.tb-scroll {{overflow-x:auto;max-height:500px;overflow-y:auto;border:1px solid var(--border);border-radius:var(--radius)}}
.chart-card {{background:var(--surface);border:1px solid var(--border);border-radius:var(--radius);overflow:hidden;transition:box-shadow 0.15s}}
.chart-card:hover {{box-shadow:var(--shadow-md)}}
.chart-title {{background:var(--surface-alt);padding:6px 10px;font-size:11.5px;font-weight:600;color:var(--ink-secondary);border-bottom:1px solid var(--border-light)}}
.chart-card img {{width:100%;max-height:90vh;height:auto;object-fit:contain;display:block;background:var(--surface-alt)}}
details {{margin:6px 0}}
summary {{cursor:pointer;font-weight:600;color:var(--accent);padding:4px 0;font-size:12.5px}}
.footer {{margin-top:40px;padding:20px;text-align:center;color:var(--ink-secondary);font-size:11.5px;border-top:1px solid var(--border)}}
.chart-gallery {{display:grid;grid-template-columns:1fr;gap:12px;margin:12px 0}}
/* AI */
#ai-btn{{position:fixed;top:16px;right:24px;width:36px;height:36px;background:linear-gradient(135deg,#667eea,#764ba2);color:#fff;border:none;border-radius:50%;font-size:11px;font-weight:700;cursor:pointer;z-index:300;box-shadow:0 2px 8px rgba(102,126,234,.3)}}
#ai-panel{{position:fixed;top:60px;right:24px;width:300px;max-height:80vh;overflow-y:auto;background:var(--surface);border:1px solid var(--border);border-radius:var(--radius);box-shadow:0 8px 24px rgba(0,0,0,.12);z-index:299}}
.ai-stage-summary{{font-size:11px;color:var(--ink);padding:6px 10px;margin:4px 0;background:var(--accent-subtle);border-radius:4px;display:none}}
.ai-stage-btn{{font-size:10px;font-weight:700;padding:1px 8px;margin-left:6px;background:linear-gradient(135deg,#667eea,#764ba2);color:#fff;border:none;border-radius:10px;cursor:pointer;vertical-align:middle;opacity:.8}}
.ai-stage-btn:hover{{opacity:1}}
.ai-section{{margin:16px 0;border:1px solid var(--border);border-radius:var(--radius);overflow:hidden}}
.ai-header{{display:flex;align-items:center;gap:10px;padding:10px 14px;background:linear-gradient(135deg,#667eea22,#764ba222);border-bottom:1px solid var(--border)}}
.ai-title{{font-weight:700;font-size:14px;color:#5a4fcf}}
.ai-model-tag{{font-size:10px;color:var(--ink-secondary);background:var(--surface-alt);padding:2px 8px;border-radius:10px}}
.ai-brief{{padding:12px 14px;font-size:13px;line-height:1.7;color:var(--ink)}}
.ai-detail-toggle{{border-top:1px solid var(--border-light)}}
.ai-detail-toggle summary{{cursor:pointer;font-size:12px;font-weight:600;color:var(--accent);padding:8px 14px;background:var(--surface-alt)}}
.ai-detail-content{{padding:8px 14px 12px}}
.ai-block{{margin:8px 0}}
.ai-block-title{{font-weight:600;font-size:11px;color:var(--accent);margin-bottom:2px}}
.ai-block-text{{font-size:12px;line-height:1.7;color:var(--ink-secondary)}}
.ai-error{{padding:12px 14px;color:var(--ink-secondary);font-size:12px;font-style:italic;background:var(--surface-alt);border-radius:var(--radius)}}
#chat-btn{{position:fixed;bottom:72px;right:24px;width:38px;height:38px;background:linear-gradient(135deg,#667eea,#764ba2);color:#fff;border:none;border-radius:50%;font-size:16px;font-weight:700;cursor:pointer;z-index:200;box-shadow:0 2px 8px rgba(102,126,234,.3)}}
#chat-panel{{position:fixed;bottom:120px;right:24px;width:360px;max-height:500px;background:var(--surface);border:1px solid var(--border);border-radius:var(--radius);box-shadow:0 8px 24px rgba(0,0,0,.12);z-index:199;display:flex;flex-direction:column}}
.chat-header{{display:flex;justify-content:space-between;align-items:center;padding:8px 12px;border-bottom:1px solid var(--border);font-size:13px;font-weight:600}}
.chat-messages{{flex:1;overflow-y:auto;padding:8px 10px;max-height:360px;display:flex;flex-direction:column;gap:6px}}
.chat-msg{{font-size:12px;line-height:1.5;padding:6px 10px;border-radius:8px;max-width:85%}}
.chat-msg.user{{align-self:flex-end;background:var(--accent);color:#fff}}
.chat-msg.ai{{align-self:flex-start;background:var(--surface-alt);color:var(--ink)}}
.chat-input-row{{display:flex;gap:6px;padding:8px 10px;border-top:1px solid var(--border-light)}}
.chat-input-row input{{flex:1;padding:6px 8px;border:1px solid var(--border);border-radius:4px;font-size:12px}}
.chat-send-btn{{padding:6px 12px;background:var(--accent);color:#fff;border:none;border-radius:4px;font-size:12px;cursor:pointer}}
@media (max-width:800px){{.sidebar{{display:none}}.main{{margin-left:0;padding:20px 16px}}}}
@media print{{.sidebar{{display:none}}#ai-btn,#chat-btn,#ai-panel,#chat-panel{{display:none}}.main{{margin-left:0}}}}
@media (prefers-reduced-motion:reduce){{*,*::after{{transition-duration:0s!important;animation-duration:0s!important}}}}
</style>
<script id="stage-data" type="application/json">{stage_data_json}</script>
</head><body>
<button id="ai-btn" onclick="toggleAIPanel()" title="AI Settings">AI</button>
<div id="ai-panel" style="display:none">
  <div style="display:flex;justify-content:space-between;align-items:center;padding:10px 14px;border-bottom:1px solid var(--border)">
    <span style="font-weight:700;font-size:14px;color:#5a4fcf">AI Settings</span>
    <button onclick="toggleAIPanel()" style="background:none;border:none;font-size:18px;cursor:pointer;color:var(--ink-secondary)">&times;</button>
  </div>
  <div style="padding:10px 14px">
    <label style="font-size:11px;color:var(--ink-secondary)">API Key (saved in browser)</label>
    <input id="ai-api-key" type="password" placeholder="sk-..." style="width:100%;padding:6px 8px;margin:4px 0 8px;border:1px solid var(--border);border-radius:4px;font-size:12px" autocomplete="off">
    <label style="font-size:11px;color:var(--ink-secondary)">Provider</label>
    <select id="ai-provider" onchange="onProviderChange()" style="width:100%;padding:6px 8px;margin:4px 0 8px;border:1px solid var(--border);border-radius:4px;font-size:12px">
      <option value="openai">OpenAI</option><option value="deepseek">DeepSeek</option><option value="moonshot">Kimi (Moonshot)</option><option value="ollama">Ollama</option><option value="custom">Custom</option>
    </select>
    <label style="font-size:11px;color:var(--ink-secondary)">Model</label>
    <input id="ai-model" value="gpt-4o-mini" style="width:100%;padding:6px 8px;margin:4px 0 8px;border:1px solid var(--border);border-radius:4px;font-size:12px">
  </div>
</div>
<button id="chat-btn" onclick="toggleChat()" title="Ask about this report">?</button>
<div id="chat-panel" style="display:none">
  <div class="chat-header"><span>Report Q&A</span><button onclick="toggleChat()" style="background:none;border:none;color:var(--ink-secondary);font-size:16px;cursor:pointer">&times;</button></div>
  <div class="chat-messages" id="chat-msgs"></div>
  <div class="chat-input-row">
    <input id="chat-input" type="text" placeholder="Ask about this report..." onkeydown="if(event.key==='Enter')sendChat()">
    <button onclick="sendChat()" class="chat-send-btn">Send</button>
  </div>
</div>
<nav class="sidebar"><h3>Metadata Pipeline</h3><ul>{nav_items}</ul></nav>
<main class="main">
  {overview}
  {s1}{s2}{s3}{s4}{s_hr}{s5}
  {ai_area}
  <div class="footer">Public Metadata Pipeline — {now}</div>
</main>
<script>
(function(){{
  document.querySelectorAll('.nav-header').forEach(function(h){{
    h.addEventListener('click',function(){{
      this.classList.toggle('collapsed');
      var s=this.nextElementSibling;
      while(s&&s.classList.contains('nav-sub')){{s.classList.toggle('collapsed');s=s.nextElementSibling;}}
    }});
  }});
  document.querySelectorAll('.chart-gallery').forEach(function(g){{
    var n=g.querySelectorAll('.chart-card').length;var c=n<=1?1:n===2?2:n===3?3:4;
    g.style.gridTemplateColumns='repeat('+c+', 1fr)';
  }});
}})();
var stageData=JSON.parse(document.getElementById('stage-data').textContent);
var saved=sessionStorage.getItem('ai_key');
if(saved)document.getElementById('ai-api-key').value=saved;
document.getElementById('ai-provider').value=sessionStorage.getItem('ai_provider')||'openai';
document.getElementById('ai-model').value=sessionStorage.getItem('ai_model')||'gpt-4o-mini';
function onProviderChange(){{
  var p=document.getElementById('ai-provider').value;
  var models={{openai:'gpt-4o-mini',deepseek:'deepseek-v4-pro',moonshot:'kimi-k2.6',ollama:'qwen2.5:7b',custom:''}};
  document.getElementById('ai-model').value=models[p]||'';
}}
function toggleAIPanel(){{var p=document.getElementById('ai-panel');p.style.display=p.style.display==='none'?'block':'none';}}
function getAIURL(){{
  var p=document.getElementById('ai-provider').value;
  if(p==='openai')return'https://api.openai.com/v1/chat/completions';
  if(p==='deepseek')return'https://api.deepseek.com/v1/chat/completions';
  if(p==='moonshot')return'https://api.moonshot.cn/v1/chat/completions';
  if(p==='ollama')return'http://localhost:11434/v1/chat/completions';
  return prompt('Enter API base URL:')||'';
}}
function runStageAI(sn){{
  var key=document.getElementById('ai-api-key').value;
  if(!key){{toggleAIPanel();alert('Please enter your API key first');return;}}
  sessionStorage.setItem('ai_key',key);sessionStorage.setItem('ai_provider',document.getElementById('ai-provider').value);sessionStorage.setItem('ai_model',document.getElementById('ai-model').value);
  var secId='stage-'+sn.replace('s1_','').replace('s2_','').replace('s3_','').replace('s4_','');
  if(sn.startsWith('s1'))secId='stage-search';
  else if(sn.startsWith('s2'))secId='stage-info';
  else if(sn.startsWith('s3'))secId='stage-down';
  else if(sn.startsWith('s4'))secId='stage-plot';
  var target=document.querySelector('#'+secId+' .ai-stage-summary');
  if(!target){{target=document.createElement('p');target.className='ai-stage-summary';var h=document.querySelector('#'+secId+' h2');if(h)h.parentNode.insertBefore(target,h.nextSibling);}}
  target.style.display='block';target.textContent='Generating...';
  var btn=document.querySelector('#'+secId+' .ai-stage-btn');if(btn)btn.disabled=true;
  var info=stageData[sn]||{{}};if(typeof info==='object')info=JSON.stringify(info);
  var sysMsg='你是生物信息学研究员。仅基于数据陈述。中文输出。学术书面语。不捏造。';
  fetch(getAIURL(),{{method:'POST',headers:{{'Content-Type':'application/json','Authorization':'Bearer '+key}},
    body:JSON.stringify({{model:document.getElementById('ai-model').value,
    messages:[{{role:'system',content:sysMsg}},{{role:'user',content:'基于以下管线数据撰写1-2句中文学术总结: '+info}}],
    temperature:0.3,max_tokens:300}})
  }}).then(r=>r.json()).then(d=>{{target.textContent=d.choices[0].message.content.trim();if(btn)btn.disabled=false;}}).catch(e=>{{target.textContent='Error: '+e;if(btn)btn.disabled=false;}});
}}
function toggleChat(){{var p=document.getElementById('chat-panel');p.style.display=p.style.display==='none'?'flex':'none';}}
function sendChat(){{
  var inp=document.getElementById('chat-input');var q=inp.value.trim();if(!q)return;
  var key=sessionStorage.getItem('ai_key');if(!key){{toggleAIPanel();alert('Please enter API key first');return;}}
  var msgs=document.getElementById('chat-msgs');
  msgs.innerHTML+='<div class=\"chat-msg user\">'+q.replace(/</g,'&lt;')+'</div>';inp.value='';msgs.scrollTop=msgs.scrollHeight;
  var thinking=document.createElement('div');thinking.className='chat-msg ai';thinking.textContent='...';msgs.appendChild(thinking);
  fetch(getAIURL(),{{method:'POST',headers:{{'Content-Type':'application/json','Authorization':'Bearer '+key}},
    body:JSON.stringify({{model:document.getElementById('ai-model').value,
    messages:[{{role:'system',content:'你是生物信息学助手。基于以下元数据管线数据回答。中文。仅陈述数据支持的事实: '+JSON.stringify(stageData)}},{{role:'user',content:q}}],
    temperature:0.3,max_tokens:300}})
  }}).then(r=>r.json()).then(d=>{{thinking.textContent=d.choices[0].message.content.trim();msgs.scrollTop=msgs.scrollHeight;}}).catch(e=>{{thinking.textContent='Error: '+e;}});
}}
</script>
</body></html>'''

    if out_html:
        with open(out_html, "w", encoding="utf-8") as f:
            f.write(html)
    return html

# ═══════════════════════════════════════════
# Main
# ═══════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(description="Public Metadata Pipeline — HTML Report Generator")
    parser.add_argument("-d", "--dir", required=True, help="Pipeline output root directory")
    parser.add_argument("-o", "--output", default=None, help="HTML output path")
    parser.add_argument("--ai-key", default="", help="API key for AI summary")
    parser.add_argument("--ai-provider", default="openai", choices=["openai","deepseek","moonshot","ollama","custom"])
    parser.add_argument("--ai-model", default="gpt-4o-mini", help="LLM model name")
    parser.add_argument("--ai-base-url", default="", help="Custom API base URL")
    args = parser.parse_args()

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

    out = Path(args.dir)
    if not out.exists():
        sys.exit(f"Directory not found: {args.dir}")

    print("Collecting search summary...")
    search_data = collect_search_summary(args.dir)
    print(f"  {'Found' if search_data else 'Not found'}")

    print("Collecting info summary...")
    info_data = collect_info_summary(args.dir)
    print(f"  {'Found' if info_data else 'Not found'}")

    print("Collecting plot images...")
    plot_images = collect_plot_images(args.dir)
    print(f"  {'Found' if plot_images else 'Not found'}")

    print("Collecting download summary...")
    download_data = collect_download_summary(args.dir)
    print(f"  {download_data.get('rows', 0)} report rows, {download_data.get('fastq_count', 0)} fastq.gz, "
          f"{download_data.get('failed', 0)} failed, {download_data.get('sra_left', 0)} .sra unconverted")

    print("Writing handoff manifest (sample_handoff.csv)...")
    handoff_data = write_handoff_manifest(args.dir)
    print(f"  {handoff_data['ready']}/{handoff_data['runs']} runs ready → {handoff_data['path']}")

    print("Collecting host reference summary...")
    hostref_data = collect_hostref_summary(args.dir)
    if hostref_data:
        print(f"  indexes built: {hostref_data['built']}/4, genome: {'yes' if hostref_data['genome'] else 'no'}")
    else:
        print("  Not found (run hostref stage first)")

    ai_html = ""
    if getattr(generate_ai_summary, '_api_key', ''):
        print("Generating AI summary...")
        ai_html = generate_ai_summary(search_data, info_data, plot_images)

    out_html = Path(args.output) if args.output else out / "Pipeline_Summary_Report.html"
    print("Generating HTML report...")
    generate_html(args.dir, search_data, info_data, plot_images, ai_html, out_html,
                  download_data, handoff_data, hostref_data)
    print(f"  Report: {out_html} ({out_html.stat().st_size / 1024:.0f} KB)")

if __name__ == "__main__":
    main()
