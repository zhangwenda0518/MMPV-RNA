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
_D = _layout_dirs(os.environ.get("MMPV_IO_LAYOUT"))
_HOSTREF_ROOT = str(Path(_D['h_genome']).parent)   # legacy=host_reference / standard=05_HostRef

try:
    import pandas as pd
    HAS_PANDAS = True
except ImportError:
    HAS_PANDAS = False

SCRIPT_DIR = Path(__file__).resolve().parent

if sys.platform == 'win32':
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')

# 字段完整度判定: 这些"看着有值"的占位串一律不算已填 (astype(str) 的产物含 nan/None)
NULL_TOKENS = {"", "nan", "none", "na", "n/a", "null"}

def safe_read_csv(fp, sep=","):
    """Read CSV/TSV, auto-detect separator; 返回**全表**.

    前几行只用来探测分隔符/列结构, 绝不把预览当结果返回 ——
    历史 bug: sep="," 分支直接 return 了 nrows=5 的预览, 报告里的统计
    (行数/字段完整度/分布) 全是按 5 行算的, 交付 HTML 的统计数字失实。
    """
    p = Path(fp)
    if not p.is_file(): return None
    try:
        if not HAS_PANDAS: return None
        use_sep = sep
        # Try comma first, then tab (预览只用于判定, 下面重读全表)
        try:
            df = pd.read_csv(p, sep=sep, nrows=5)
            if len(df.columns) <= 1 and sep == ",":
                use_sep = "\t"
        except:
            if sep == ",":
                use_sep = "\t"
            else:
                return None
        return pd.read_csv(p, sep=use_sep)
    except Exception:
        return None

def img_to_base64(path, max_kb=2000):
    if not path or not Path(path).exists(): return None
    if Path(path).stat().st_size / 1024 > max_kb: return None
    with open(path, 'rb') as f:
        return base64.b64encode(f.read()).decode('utf-8')

def _esc(s):
    return str(s).replace("&","&amp;").replace("<","&lt;").replace(">","&gt;").replace('"',"&quot;")

def _fmt_pct(v, nd=1):
    """把百分数格式化成给人看的样子 —— 避免表格里出现 4.6e-05 这种科学计数法。

    只用于展示; TSV/JSON 里保留原始精度给下游。
    """
    try:
        x = float(v)
    except (TypeError, ValueError):
        return str(v)
    if x == 0:
        return "0"
    if abs(x) < 0.001:
        return "<0.001"
    if abs(x) < 1:
        return f"{x:.4f}".rstrip("0").rstrip(".")
    return f"{x:.{nd}f}"

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
        non_null = df[col].notna()
        # astype(str) 会把 NaN 变 "nan"、NaT 变 "NaT" —— 不算已填, 否则完整度虚高
        vals = df[col].astype(str).str.strip().str.lower()
        non_empty = int((non_null & ~vals.isin(NULL_TOKENS)).sum())
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
    """Find SCI landscape plot images.

    同一次运行 plot 阶段会对 search / info 两个数据源各出一套 (from_search/,
    from_info/), 这里全部收集 (list), 不再只取首个匹配。
    """
    out = Path(out_dir)
    images = {}
    for name, pat in [("landscape", "Combined_Landscape_Full.png"), ("landscape_pdf", "Combined_Landscape_Full.pdf")]:
        hits = sorted(out.rglob(pat))
        if hits:
            images[name] = [str(p) for p in hits]
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

def collect_prescreen_summary(out_dir, top_n=15):
    """读取预筛产物 (sample_decision.tsv + prescreen_meta.json) 汇总成报告用数据。

    纯 csv/json 读取, 不依赖 pandas —— 预筛是独立阶段, 不应因缺 pandas 而丢失。
    返回 None 表示没跑过 prescreen 阶段。
    """
    out = Path(out_dir)
    dec_path = next(iter(out.rglob("sample_decision.tsv")), None)
    if not dec_path:
        return None
    try:
        with open(dec_path, "r", encoding="utf-8", errors="replace") as f:
            rows = list(csv.DictReader(f, delimiter="\t"))
    except Exception:
        return None
    if not rows:
        return None

    def _f(v, d=0.0):
        try:
            return float(v)
        except (TypeError, ValueError):
            return d

    decision = Counter((r.get("decision") or "").strip() for r in rows)
    status = Counter((r.get("status") or "").strip() for r in rows)
    omatch = Counter((r.get("organism_match") or "").strip() for r in rows)
    signal = Counter((r.get("virus_signal") or "").strip() for r in rows)

    ok_rows = [r for r in rows if (r.get("status") or "").strip() == "OK"]
    viral_rows = [r for r in ok_rows if _f(r.get("virus_pct_ident")) > 0]

    vt = Counter()
    for r in ok_rows:
        tv = (r.get("top_virus") or "").strip()
        if tv and tv != "Viruses (未细分)":
            vt[tv] += 1

    mism = sorted([r for r in ok_rows if (r.get("organism_match") or "").strip() == "MISMATCH"],
                  key=lambda x: -_f(x.get("top_organism_pct")))[:top_n]
    top_runs = sorted(viral_rows, key=lambda x: -_f(x.get("virus_pct_ident")))[:top_n]

    by_cra = Counter()
    for r in rows:
        if (r.get("status") or "").strip() != "OK":
            by_cra[(r.get("cra") or "(unresolved)").strip()] += 1

    meta = {}
    mp = dec_path.parent / "prescreen_meta.json"
    if mp.is_file():
        try:
            meta = _json.loads(mp.read_text(encoding="utf-8"))
        except Exception:
            meta = {}

    try:
        rel = str(dec_path.relative_to(out))
    except ValueError:
        rel = str(dec_path)

    return {
        "path": rel, "total": len(rows), "ok": len(ok_rows), "with_virus": len(viral_rows),
        "decision": dict(decision), "status": dict(status),
        "organism_match": dict(omatch), "virus_signal": dict(signal),
        "virus_taxa": vt.most_common(top_n), "mismatch": mism, "top_runs": top_runs,
        "unresolved_by_cra": by_cra.most_common(12), "meta": meta,
        "virus_profile": _collect_virus_profile(dec_path.parent, out, top_n),
    }


def _collect_virus_profile(prescreen_dir, out, top_n=15):
    """读 stat_virus_profile.py 的产物 (病毒层面画像), 供 HTML 展示。

    病毒本身也是结果 —— 有哪些病毒 / 各占多少 / 在多少样本里出现。
    返回 None 表示没跑过画像 (老版本产物, 或 --no-virus-profile)。
    纯 csv/json 读取, 不依赖 pandas。
    """
    vdir = Path(prescreen_dir) / "virus_profile"
    prev_tsv = vdir / "virus_prevalence.tsv"
    if not prev_tsv.is_file():
        return None
    try:
        with open(prev_tsv, "r", encoding="utf-8", errors="replace") as f:
            prow = list(csv.DictReader(f, delimiter="\t"))
    except Exception:
        return None
    if not prow:
        return None

    def _f(v, d=0.0):
        try:
            return float(v)
        except (TypeError, ValueError):
            return d

    vmeta = {}
    mj = vdir / "virus_profile_meta.json"
    if mj.is_file():
        try:
            vmeta = _json.loads(mj.read_text(encoding="utf-8"))
        except Exception:
            vmeta = {}

    by_leaf = vmeta.get("prevalence_by", "leaf") == "leaf"
    ncol = "n_leaf_samples" if by_leaf else "n_samples"
    pcol = "leaf_prevalence_pct" if by_leaf else "prevalence_pct"
    min_n = int(vmeta.get("min_samples", 2) or 2)

    ranked = sorted(prow, key=lambda r: -_f(r.get(ncol)))
    shown = [r for r in ranked if _f(r.get(ncol)) >= min_n][:top_n]

    figs = sorted((vdir / "figures").glob("*.png")) if (vdir / "figures").is_dir() else []
    # 优先展示：检出率条形图 → 组合图
    pick = ([f for f in figs if "prevalence_bar" in f.name]
            + [f for f in figs if "Combined" in f.name])

    try:
        rel = str(vdir.relative_to(out))
    except ValueError:
        rel = str(vdir)

    return {
        "dir": rel,
        "scope": vmeta.get("plot_scope", "leaf-node" if by_leaf else "any-rank"),
        "n_taxa": vmeta.get("n_viral_taxa_total", len(prow)),
        "n_leaf": sum(1 for r in prow if (r.get("is_leaf_node") or "") == "yes"),
        "n_samples": vmeta.get("n_samples_with_taxonomy", "?"),
        "n_viral": vmeta.get("n_samples_with_virus", "?"),
        "by_class": vmeta.get("by_class", {}) or {},
        "rows": [{
            "virus": r.get("virus", ""), "class": r.get("virus_class", ""),
            "rank": r.get("rank", ""), "n": r.get(ncol, ""), "pct": r.get(pcol, ""),
            "mean_pct": r.get("mean_pct", ""), "share": r.get("hits_share_pct", ""),
        } for r in shown],
        "figures": [str(f) for f in pick],
    }

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
    prompt = f"""You are a bioinformatics researcher. Write an English research brief based on the results of the public data search pipeline below.

## Search Results
- Records found: {s.get('total','?')} (SRA={s.get('sra','?')}, GSA={s.get('gsa','?')})
- BioProjects: {s.get('bioprojects','?')}, species: {s.get('species','?')}
- Temporal range: {s.get('years','?')}

## Metadata Quality
- Core14 unified metadata: {i.get('total','?')} records
- PMID coverage: {i.get('pmid_coverage','?')} records with literature links

## Metadata Field Completeness
{_json.dumps({k: v.get('pct','?') for k,v in (i.get('fields',{}) or {}).items()}, ensure_ascii=False)}

Output format:
[BRIEF]
150-200 word English brief: research objective, search scale, key findings (database distribution, metadata quality, literature coverage)

[DETAILED]
Output a 300-400 word English detailed report in [BACKGROUND][METHODS][RESULTS][DISCUSSION] format.

State only facts supported by the data. English. Academic prose. Do not fabricate."""

    try:
        url = base_url or "https://api.openai.com/v1/chat/completions"
        body = _json.dumps({
            "model": model, "messages": [
                {"role": "system", "content": "You are a bioinformatics researcher. State only facts supported by the data. Output in English."},
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
                  download_data=None, handoff_data=None, hostref_data=None, prescreen_data=None):
    out = Path(out_dir); now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    s = search_data or {}; i = info_data or {}; d = download_data or {}
    h = handoff_data or {}; hr = hostref_data or {}; ps = prescreen_data or {}

    # Build stage data JSON for browser AI
    sd = {
        "_overview": {"records": s.get("total","?"), "sra": s.get("sra","?"), "gsa": s.get("gsa","?"),
                      "bioprojects": s.get("bioprojects","?"), "species": s.get("species","?"), "years": s.get("years","?")},
        "s1_search": {"total": s.get("total","?"), "sra": s.get("sra","?"), "gsa": s.get("gsa","?"),
                      "bioprojects": s.get("bioprojects","?"), "species": s.get("species","?"), "years": s.get("years","?")},
        "s2_info": {"total": i.get("total","?"), "pmid_coverage": i.get("pmid_coverage","?"),
                    "fields": {k: v.get("pct","?") for k,v in (i.get("fields",{}) or {}).items()}},
        "s3_prescreen": {"total": ps.get("total","?"), "ok": ps.get("ok","?"),
                         "with_virus": ps.get("with_virus","?"),
                         "decision": ps.get("decision",{}),
                         "organism_match": ps.get("organism_match",{}),
                         "virus_taxa": [list(x) for x in (ps.get("virus_taxa") or [])],
                         "virus_profile": {
                             "scope": (ps.get("virus_profile") or {}).get("scope"),
                             "n_taxa": (ps.get("virus_profile") or {}).get("n_taxa"),
                             "n_leaf": (ps.get("virus_profile") or {}).get("n_leaf"),
                             "by_class": (ps.get("virus_profile") or {}).get("by_class"),
                             "top": [[r.get("virus"), r.get("n"), r.get("pct")]
                                     for r in ((ps.get("virus_profile") or {}).get("rows") or [])[:15]],
                         }},
        "s3_down": {"report_rows": d.get("rows","?"), "status": d.get("status",{}),
                    "failed": d.get("failed","?"), "fastq_files": d.get("fastq_count","?"),
                    "fastq_gb": round(d.get("fastq_gb",0),1), "sra_left": d.get("sra_left","?")},
        "s4_plot": {"note": "6-panel SCI landscape figure: temporal trends, database composition, top organizations, tissues, locations, growth stages"},
        "s5_hostref": {"genome": hr.get("genome"), "genome_mb": hr.get("genome_mb"),
                       "indexes": hr.get("indexes", {}), "built": hr.get("built", "?")},
        "s6_handoff": {"runs": h.get("runs","?"), "ready": h.get("ready","?"),
                       "missing": h.get("missing","?"), "manifest": h.get("path","sample_handoff.csv")},
    }
    # 数据里若含 "</script>" 会提前截断 <script> 脚本块 —— JSON 层把 "<" 转义为
    # \u003c (json.loads / JSON.parse 均能还原, 语义不变; 全文件仅此一处注入点)
    stage_data_json = _json.dumps(sd, ensure_ascii=False).replace("<", "\\u003c")

    # ── Sidebar ──
    nav_items = '<li class="nav-item"><a href="#overview" class="nav-link active">Overview</a></li>'
    nav_items += '<li class="nav-header">Pipeline Stages</li>'
    nav_items += '<div class="nav-sub">'
    for sid, label in [("stage-search","S1: Search"),("stage-info","S2: Info"),("stage-prescreen","S3: Prescreen"),("stage-down","S4: Download"),("stage-plot","S5: Plot"),("stage-hostref","S6: HostRef"),("stage-handoff","S7: Handoff")]:
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
            s2 += '<div class="tb-scroll"><table style="font-size:12px"><thead><tr><th>Field</th><th>Complete</th><th>Rate</th></tr></thead><tbody>'
            for fn, fv in sorted(fields.items()):
                s2 += f'<tr><td>{_esc(fn)}</td><td>{fv.get("complete","?")}</td><td>{fv.get("pct","?")}</td></tr>'
            s2 += '</tbody></table></div></details>'
        # Distribution tables
        dists = i.get("distributions", {})
        for col, items in dists.items():
            if not items: continue
            s2 += f'<details style="margin-top:6px"><summary style="font-weight:600;color:var(--accent);cursor:pointer;font-size:12px">Top {_esc(col)}</summary>'
            s2 += '<div class="tb-scroll"><table style="font-size:12px"><thead><tr><th>Value</th><th>Count</th></tr></thead><tbody>'
            for val, cnt in items[:10]:
                s2 += f'<tr><td>{_esc(val[:50])}</td><td>{cnt}</td></tr>'
            s2 += '</tbody></table></div></details>'
        s2 += f'<p style="font-size:12px;color:var(--ink-secondary);margin-top:4px">Source: {_esc(i.get("path","?"))}</p>'
    else:
        s2 += '<p style="color:var(--ink-secondary)">Info output not found. Run Stage 2 (info) first.</p>'
    s2 += '</section>'

    # ── Stage 3: Prescreen ──
    s_ps = f'''<section id="stage-prescreen"><h2>S3: Prescreen — Viral Pre-screening &amp; Species Identification <button onclick="runStageAI('s3_prescreen')" class="ai-stage-btn">AI</button></h2>
    <p style="color:var(--ink-secondary);font-size:13px">Based on <b>pre-computed taxonomic profiles</b> from NCBI STAT / CNCB GSA — answers two questions <b>without downloading any sequencing data</b>:
    (1) <b>Plant species identification</b>: cross-check the observed top organism against the metadata ScientificName to flag mislabeled samples;
    (2) <b>Viral pre-screening</b>: assign a decision from the viral fraction to filter the download stage and avoid wasting resources.</p>'''
    if ps:
        dec = ps.get("decision", {})
        s_ps += '<div class="card-row">'
        s_ps += f'<div class="card"><div class="value">{ps.get("total","?")}</div><div class="label">Runs</div></div>'
        s_ps += f'<div class="card"><div class="value">{ps.get("ok","?")}</div><div class="label">With Taxonomy</div></div>'
        s_ps += f'<div class="card"><div class="value">{ps.get("with_virus","?")}</div><div class="label">Virus Detected</div></div>'
        for k, label, color in (("DOWNLOAD", "Recommended", "#16a34a"),
                                ("REVIEW", "Needs Review", "#d97706"),
                                ("SKIP", "Skippable", "#6b7280"),
                                ("LOCAL_CHECK", "Local Check", "#2563eb")):
            s_ps += (f'<div class="card"><div class="value" style="color:{color}">{dec.get(k, 0)}</div>'
                     f'<div class="label">{label}</div></div>')
        s_ps += '</div>'

        s_ps += ('<p style="font-size:12px;color:var(--ink-secondary);margin-top:6px">'
                 'decision: <b>DOWNLOAD</b> virus detected and organism match consistent · '
                 '<b>REVIEW</b> virus detected but top organism conflicts with metadata (still downloaded, for manual review) · '
                 '<b>SKIP</b> taxonomy available but no viral signal · '
                 '<b>LOCAL_CHECK</b> no pre-computed taxonomy available, reads must be downloaded for local alignment</p>')

        # ── 物种鉴定 ──
        om = ps.get("organism_match", {})
        if any(om.values()):
            s_ps += '<details style="margin-top:8px" open><summary style="font-weight:600;color:var(--accent);cursor:pointer;font-size:13px">Plant Species Identification (metadata vs observed top organism)</summary>'
            s_ps += '<div class="tb-scroll"><table style="font-size:12px"><thead><tr><th>Call</th><th>Count</th><th>Meaning</th></tr></thead><tbody>'
            for k, dsc in (("MATCH", "consistent"), ("PARTIAL", "same genus, different species (related / hybrid)"),
                           ("MISMATCH", "different genus with sufficient fraction — metadata likely mislabeled"),
                           ("COARSE", "observed only above genus rank (insufficient resolution)"),
                           ("WEAK", "observed fraction too low (insufficient evidence)"),
                           ("UNKNOWN", "one side missing, no conclusion")):
                v = om.get(k, 0)
                hl = ' style="color:#b91c1c;font-weight:700"' if (k == "MISMATCH" and v) else ''
                s_ps += f'<tr><td{hl}>{k}</td><td{hl}>{v}</td><td>{dsc}</td></tr>'
            s_ps += '</tbody></table></div>'
            mism = ps.get("mismatch") or []
            if mism:
                s_ps += ('<div class="tb-scroll" style="margin-top:6px"><table style="font-size:12px">'
                         '<thead><tr><th>Run</th><th>Metadata Species</th><th>Observed Top Species</th><th>Observed %</th>'
                         '<th>Confidence</th><th>Top Virus</th></tr></thead><tbody>')
                for r in mism:
                    s_ps += (f'<tr><td>{_esc(r.get("acc",""))}</td>'
                             f'<td>{_esc(r.get("meta_species",""))}</td>'
                             f'<td style="color:#b91c1c">{_esc(r.get("top_organism",""))}</td>'
                             f'<td>{_esc(r.get("top_organism_pct",""))}</td>'
                             f'<td>{_esc(r.get("organism_conf",""))}</td>'
                             f'<td>{_esc(r.get("top_virus",""))}</td></tr>')
                s_ps += '</tbody></table></div>'
                s_ps += (f'<p style="font-size:11px;color:var(--ink-secondary)">'
                         f'Showing the top {len(mism)} samples with suspected metadata errors'
                         f' (full list: organism_match=MISMATCH in sample_decision.tsv)</p>')
            s_ps += '</details>'

        # ── 病毒类群分布 ──
        vt = ps.get("virus_taxa") or []
        if vt:
            s_ps += '<details style="margin-top:6px"><summary style="font-weight:600;color:var(--accent);cursor:pointer;font-size:12px">Top Detected Viral Taxa</summary>'
            s_ps += '<div class="tb-scroll"><table style="font-size:12px"><thead><tr><th>Viral Taxon (most specific node)</th><th>Samples</th></tr></thead><tbody>'
            for name, cnt in vt:
                s_ps += f'<tr><td>{_esc(name)}</td><td>{cnt}</td></tr>'
            s_ps += '</tbody></table></div></details>'

        # ── 病毒占比 Top ──
        tr = ps.get("top_runs") or []
        if tr:
            s_ps += '<details style="margin-top:6px"><summary style="font-weight:600;color:var(--accent);cursor:pointer;font-size:12px">Top Samples by Viral Fraction</summary>'
            s_ps += ('<div class="tb-scroll"><table style="font-size:12px"><thead><tr><th>Run</th><th>Source</th>'
                     '<th>Virus %</th><th>Virus % (identified)</th><th>Top Virus</th><th>Reads</th>'
                     '<th>Top Non-Viral Organism</th></tr></thead><tbody>')
            for r in tr:
                s_ps += (f'<tr><td>{_esc(r.get("acc",""))}</td><td>{_esc(r.get("source",""))}</td>'
                         f'<td>{_esc(r.get("virus_pct",""))}</td><td>{_esc(r.get("virus_pct_ident",""))}</td>'
                         f'<td>{_esc(r.get("top_virus",""))}</td><td>{_esc(r.get("analysed",""))}</td>'
                         f'<td>{_esc(r.get("top_organism",""))}</td></tr>')
            s_ps += '</tbody></table></div></details>'

        # ── 无分类数据 ──
        uc = ps.get("unresolved_by_cra") or []
        if uc:
            n_unres = sum(v for _, v in uc)
            s_ps += '<details style="margin-top:6px"><summary style="font-weight:600;color:var(--accent);cursor:pointer;font-size:12px">Samples without pre-computed taxonomy (need local alignment)</summary>'
            s_ps += '<div class="tb-scroll"><table style="font-size:12px"><thead><tr><th>Project (CRA)</th><th>Runs without data</th></tr></thead><tbody>'
            for k, v in uc:
                s_ps += f'<tr><td>{_esc(k)}</td><td>{v}</td></tr>'
            s_ps += '</tbody></table></div>'
            s_ps += (f'<p style="font-size:11px;color:var(--ink-secondary)">Runs without pre-computed taxonomy: {n_unres}. '
                     f'This is <b>not an analysis failure</b> — the upstream simply has no pre-computed data for these runs, so species can only be determined by downloading reads for local alignment.</p>')
            s_ps += '</details>'

        # ---- 病毒层面画像：病毒本身也是结果 ----
        vp = ps.get("virus_profile") or {}
        if vp:
            s_ps += ('<div style="margin-top:16px;padding-top:12px;border-top:1px dashed var(--border)">')
            s_ps += ('<h3 style="font-size:14px;margin:0 0 8px">Virus Profile '
                     '<span style="font-weight:400;font-size:11px;color:var(--ink-secondary)">'
                     '— which viruses / their share / in how many samples'
                     f' (scope: {_esc(vp.get("scope",""))})</span></h3>')
            s_ps += '<div class="card-row">'
            s_ps += (f'<div class="card"><div class="value">{vp.get("n_taxa","?")}</div>'
                     f'<div class="label">Viral Taxa</div></div>')
            s_ps += (f'<div class="card"><div class="value">{vp.get("n_leaf","?")}</div>'
                     f'<div class="label">Leaf Nodes (finest detected)</div></div>')
            s_ps += (f'<div class="card"><div class="value">{vp.get("n_viral","?")}</div>'
                     f'<div class="label">Samples with Virus</div></div>')
            for c, n in sorted((vp.get("by_class") or {}).items(), key=lambda x: -x[1])[:3]:
                s_ps += (f'<div class="card"><div class="value" style="font-size:20px">{n}</div>'
                         f'<div class="label">{_esc(c)}</div></div>')
            s_ps += '</div>'

            for fp in (vp.get("figures") or [])[:2]:
                b64 = img_to_base64(fp, 6000)
                if b64:
                    nm = Path(fp).name
                    s_ps += (f'<div class="chart-card"><div class="chart-title">{_esc(nm)}</div>'
                             f'<img src="data:image/png;base64,{b64}" loading="lazy" alt="{_esc(nm)}"></div>')
                elif Path(fp).is_file():
                    mb = Path(fp).stat().st_size / 1024 ** 2
                    s_ps += (f'<div class="chart-card"><div class="chart-title">{_esc(Path(fp).name)}</div>'
                             f'<p style="padding:12px;color:#d97706;font-size:13px">⚠ Image too large ({mb:.1f} MB) — not embedded, open the file directly: <code>{_esc(fp)}</code></p></div>')

            s_ps += ('<details style="margin-top:8px" open><summary style="font-weight:600;'
                     'color:var(--accent);cursor:pointer;font-size:13px">'
                     f'Top Cross-Sample Prevalence (scope: {_esc(vp.get("scope",""))})</summary>')
            s_ps += ('<div class="tb-scroll"><table style="font-size:12px">'
                     '<thead><tr><th>Virus</th><th>Class</th><th>Rank</th><th>Detected Samples</th>'
                     '<th>Prevalence %</th><th>Mean Abundance %</th><th>Hit Share %</th></tr></thead><tbody>')
            for r in vp.get("rows", []):
                s_ps += (f'<tr><td>{_esc(r.get("virus",""))}</td><td>{_esc(r.get("class",""))}</td>'
                         f'<td>{_esc(r.get("rank",""))}</td><td>{_esc(r.get("n",""))}</td>'
                         f'<td>{_esc(_fmt_pct(r.get("pct"), 1))}</td>'
                         f'<td>{_esc(_fmt_pct(r.get("mean_pct"), 3))}</td>'
                         f'<td>{_esc(_fmt_pct(r.get("share"), 2))}</td></tr>')
            s_ps += '</tbody></table></div></details>'
            s_ps += (f'<p style="font-size:11px;color:var(--ink-secondary);margin-top:4px">'
                     f'Source: <code>{_esc(vp.get("dir",""))}</code> — '
                     f'fast k-mer classification, intended as <b>pre-screening</b>: low-abundance viruses may be missed, and '
                     f'phages / human viruses (e.g. Human gammaherpesvirus 4) are common background. '
                     f'Confirm conclusions with downstream alignment.</p>')
            s_ps += '</div>'

        s_ps += f'<p style="font-size:12px;color:var(--ink-secondary);margin-top:6px">Source: <code>{_esc(ps.get("path","?"))}</code></p>'
    else:
        s_ps += '<p style="color:var(--ink-secondary)">Pre-screen output not found. Run Stage 3 (prescreen) first — this stage downloads no data, it only reads the upstream pre-computed taxonomy.</p>'
    s_ps += '</section>'

    # ── Stage 4: Download ──
    s3 = f'''<section id="stage-down"><h2>S4: Download <button onclick="runStageAI('s3_down')" class="ai-stage-btn">AI</button></h2>
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
        s3 += '<p style="color:var(--ink-secondary)">Download report not found. Run Stage 4 (down) first.</p>'
    s3 += '</section>'

    # ── Stage 4: Plot ──
    s4 = f'''<section id="stage-plot"><h2>S5: Plot <button onclick="runStageAI('s4_plot')" class="ai-stage-btn">AI</button></h2>
    <p style="color:var(--ink-secondary);font-size:13px">SCI-grade 6-panel landscape figure: (A) Temporal Distribution, (B) Database Proportion, (C) Top Organizations, (D) Top Tissues, (E) Top Locations, (F) Top Growth Stages.</p>'''
    if plot_images.get("landscape") or plot_images.get("landscape_pdf"):
        for fp in plot_images.get("landscape") or []:
            # 双数据源 (from_search/ from_info/) 各出一套图时全部嵌入, 子目录名作小标题
            src = Path(fp).parent.name or "plot"
            b64 = img_to_base64(fp, 5000)
            if b64:
                s4 += f'<div class="chart-card"><div class="chart-title">Combined Landscape (6-panel SCI Figure) — {_esc(src)}</div><img src="data:image/png;base64,{b64}" loading="lazy" alt="SCI Landscape Figure ({_esc(src)})"></div>'
            else:
                mb = Path(fp).stat().st_size / 1024 ** 2
                s4 += (f'<div class="chart-card"><div class="chart-title">Combined Landscape (6-panel SCI Figure) — {_esc(src)}</div>'
                       f'<p style="padding:12px;color:#d97706;font-size:13px">⚠ Image too large ({mb:.1f} MB) — not embedded, open the file directly: <code>{_esc(fp)}</code></p></div>')
        for fp in plot_images.get("landscape_pdf") or []:
            try:
                href = str(Path(fp).relative_to(out)).replace("\\", "/")
            except ValueError:
                href = str(fp).replace("\\", "/")
            s4 += (f'<p style="font-size:12px;margin-top:4px">PDF (vector): '
                   f'<a href="{_esc(href)}">{_esc(fp)}</a></p>')
    else:
        s4 += '<p style="color:var(--ink-secondary)">Plot not found. Run Stage 5 (plot) to generate the 6-panel figure.</p>'
    s4 += '</section>'

    # ── Stage 5: HostRef ──
    s_hr = f'''<section id="stage-hostref"><h2>S6: HostRef — Host Reference Genome</h2>
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
    s5 = f'''<section id="stage-handoff"><h2>S7: Handoff Manifest</h2>
    <p style="color:var(--ink-secondary);font-size:13px">Downstream handoff manifest: Run → FASTQ absolute paths + Core14 metadata.
    Consumed by <code>data_preprocessing_pipeline/clean-data.py --input</code> as the file-level sample list.</p>'''
    if h:
        s5 += f'<div class="card-row"><div class="card"><div class="value">{h.get("runs","?")}</div><div class="label">Runs Total</div></div>'
        s5 += f'<div class="card"><div class="value" style="color:#16a34a">{h.get("ready","?")}</div><div class="label">Ready (files found)</div></div>'
        s5 += f'<div class="card"><div class="value" style="color:#d97706">{h.get("missing","?")}</div><div class="label">Missing (metadata only)</div></div></div>'
        s5 += f'<p style="font-size:12px;color:var(--ink-secondary)">Manifest: <code>{_esc(h.get("path","sample_handoff.csv"))}</code> — columns: Run / Layout / R1_Path / R2_Path / SE_Path / Size_MB / ScientificName / Tissue / Location / BioProject / Status</p>'
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
  {s1}{s2}{s_ps}{s3}{s4}{s_hr}{s5}
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
  else if(sn.startsWith('s3_prescreen'))secId='stage-prescreen';
  else if(sn.startsWith('s3'))secId='stage-down';
  else if(sn.startsWith('s4'))secId='stage-plot';
  var target=document.querySelector('#'+secId+' .ai-stage-summary');
  if(!target){{target=document.createElement('p');target.className='ai-stage-summary';var h=document.querySelector('#'+secId+' h2');if(h)h.parentNode.insertBefore(target,h.nextSibling);}}
  target.style.display='block';target.textContent='Generating...';
  var btn=document.querySelector('#'+secId+' .ai-stage-btn');if(btn)btn.disabled=true;
  var info=stageData[sn]||{{}};if(typeof info==='object')info=JSON.stringify(info);
  var sysMsg='You are a bioinformatics researcher. State only facts supported by the data. Output in English. Academic prose. Do not fabricate.';
  fetch(getAIURL(),{{method:'POST',headers:{{'Content-Type':'application/json','Authorization':'Bearer '+key}},
    body:JSON.stringify({{model:document.getElementById('ai-model').value,
    messages:[{{role:'system',content:sysMsg}},{{role:'user',content:'Write a 1-2 sentence academic summary in English based on the following pipeline data: '+info}}],
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
    messages:[{{role:'system',content:'You are a bioinformatics assistant. Answer based on the metadata pipeline data below. English. State only facts supported by the data: '+JSON.stringify(stageData)}},{{role:'user',content:q}}],
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

    print("Collecting pre-screen summary...")
    prescreen_data = collect_prescreen_summary(args.dir)
    if prescreen_data:
        _dec = prescreen_data.get("decision", {})
        print(f"  {prescreen_data['total']} runs, {prescreen_data['with_virus']} with virus signal, "
              f"decision: " + ", ".join(f"{k}={v}" for k, v in sorted(_dec.items())))
    else:
        print("  Not found (run prescreen stage first)")

    ai_html = ""
    if getattr(generate_ai_summary, '_api_key', ''):
        print("Generating AI summary...")
        ai_html = generate_ai_summary(search_data, info_data, plot_images)

    out_html = Path(args.output) if args.output else out / "Pipeline_Summary_Report.html"
    print("Generating HTML report...")
    generate_html(args.dir, search_data, info_data, plot_images, ai_html, out_html,
                  download_data, handoff_data, hostref_data, prescreen_data)
    print(f"  Report: {out_html} ({out_html.stat().st_size / 1024:.0f} KB)")

if __name__ == "__main__":
    main()
