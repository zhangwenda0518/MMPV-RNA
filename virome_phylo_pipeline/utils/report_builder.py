#!/usr/bin/env python3
"""
report_builder.py — 模块化系统发育分析汇总 HTML 报告

按 8 大模块 (data → phylogeny → recomb → select → time → geography → popgen) 分节,
每节展示: stage 状态徽章 + 该模块 summary 表 (csv/tsv/json → HTML table) + plots (png base64 内嵌)。
单文件零依赖输出 phylo_report.html。
"""

import os
import json
import csv
import base64
import html as _html
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from datetime import datetime

# 2026-09-15 修复 (P0-A): HTML 报告与 phylo_summary.csv 共用同一取数入口。
# 修复前 HTML 读 metrics['beta'/'r_squared'/'p_value'] (TreeTime RTT),
# 而 CSV 读 td_mean_rate/pic_r/td_relaxed_clock_p (treedater/PIC) —— 同名列
# 不同来源, 同一病毒两份产物数字互相打脸。现两处都走 headline_metrics()。
from utils.summary_metrics import headline_metrics, fmt_metric

# 跨管线统一 I/O 布局: ⑤ per-virus 模块目录名随 MMPV_IO_LAYOUT 解析
try:
    from mmpv_common.io_layout import ph_dir
except ImportError:  # 脱离包环境直接运行时自举
    import os as _os, sys as _sys
    _REPO_ROOT = _os.path.dirname(_os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
    if _REPO_ROOT not in _sys.path:
        _sys.path.insert(0, _REPO_ROOT)
    from mmpv_common.io_layout import ph_dir

# ─────────────────────────────────────────────────────────────
# 模块定义: 每模块的目录、状态 stage 集、summary 文件、plots glob
# 相对每个病毒 work_dir 解析
# ─────────────────────────────────────────────────────────────
MODULES = [
    {
        "key": "data", "name": "Data (prep + online)",
        "dir": ph_dir("data"),
        "stages": ["prep", "online"],
        # data 模块本身无汇总表, 但 seqgrouper 的分组饼图/统计表落在这里
        "summaries": ["ncbi_harvest/seqgroup_*.csv"],
        "plots": ["ncbi_harvest/*.png", "ncbi_harvest/*.pdf", "ncbi_harvest/*.svg",
                  "**/*.pdf"],
        "notes": ["combined.fasta 样本+参考序列", "dates.csv / sample_metadata.csv", "ncbi_harvest 公共序列收集"],
    },
    {
        "key": "phylogeny", "name": "Phylogeny (align + saturation + tree + splitstree)",
        "dir": ph_dir("phylogeny"),
        "stages": ["align", "saturation", "tree", "splitstree"],
        "summaries": ["saturation_out/saturation_report.csv", "qc/align_qc_report.tsv"],
        "plots": ["saturation_out/*.png", "qc/*.png", "tree_circular.png", "splitstree/*.png",
                  "gene_analysis/*.png", "cross_validation/*.png",
                  "qc/*.pdf", "tree_circular.pdf", "**/*.pdf"],
        "notes": ["MAFFT 比对 + IQ-TREE ML 树", "Xia 替换饱和检验", "比对 QC"],
    },
    {
        "key": "popgen", "name": "Population genetics (popgen + host)",
        "dir": ph_dir("popgen"),
        "stages": ["popgen", "host"],
        "summaries": ["dnasp/results.tsv", "dnasp/ld_pairs.tsv", "location_pops.tsv"],
        "plots": ["dnasp/figures/*.png", "host_analysis/*.png", "distance/*.png",
                  "popgen_metrics.png", "**/*.pdf"],
        "notes": ["单倍型多样性 + Tajima's D / Fu's Fs", "宿主分化分析"],
    },
    {
        "key": "recomb", "name": "Recombination (rdp5)",
        "dir": ph_dir("recomb"),
        "stages": ["rdp5"],
        "summaries": ["rdp5/recombination_report.md", "rdp5/*.csv"],
        "plots": ["rdp5/*.png", "rdp5/figures/*.png", "**/*.pdf"],
        "notes": ["RDP5 七方法 + mask 重组区"],
    },
    {
        "key": "select", "name": "Selection (capheine)",
        "dir": ph_dir("select"),
        "stages": ["capheine"],
        "summaries": ["capheine/drhip/combined_sites.csv", "capheine/multiqc/*summary*.txt"],
        "plots": ["capheine/**/*.png", "capheine/**/*.pdf"],
        "notes": ["HyPhy FEL/MEME/BUSTED 等 + DRHIP 位点整合"],
    },
    {
        "key": "time", "name": "Time (rtt + temporal + beast + gene_dating)",
        "dir": ph_dir("time"),
        "stages": ["rtt", "temporal", "beast", "gene_dating"],
        "summaries": ["clock_summary.tsv", "temporal_signal/drt_summary.csv",
                      "beast/merged/parameter_summary.csv",
                      "treetime_rtt/robust/robust_report.json"],
        "plots": ["treetime_rtt/*.png", "treetime_rtt/robust/*.png", "treedater_ltt/*.png",
                  "temporal_signal/real/*.png", "beast/postprocess/*.png",
                  "beast/merged/*.png", "gene_dating/*.png", "**/*.pdf"],
        "notes": ["Root-to-tip + DRT 时间信号", "BEAST 贝叶斯定年", "clock_summary 汇总各方法速率"],
    },
    {
        "key": "geography", "name": "Geography (phylogeo + geo)",
        "dir": ph_dir("geography"),
        "stages": ["phylogeo", "geo", "geo_paths"],
        "summaries": ["phylogeography/merged/parameter_summary.json",
                      "tempmig/migration_yearly_summary.csv", "geo_analysis/spatiotemporal_counts.csv",
                      # geo_paths (2026-09-16): 事件级路径/episode/LTL/扩散统计
                      "geo_analysis/pathways/pathway_qc.json",
                      "geo_analysis/pathways/pathways.csv", "geo_analysis/pathways/episodes.csv",
                      "geo_analysis/pathways/ltl_summary.csv",
                      "geo_analysis/dispersion/dispersion_summary.json",
                      "geo_analysis/timescale_consistency.json",
                      "geo_analysis/markov_jumps.json"],
        # visualization/ 里既有 SpreaD3 地图 (html) 也有发表级多面板图 (pdf)
        "plots": ["visualization/*.png", "geo_analysis/*.png", "tempmig/*.png",
                  "visualization/*.gif",
                  "**/*.pdf", "visualization/*.html", "geo_analysis/*.html"],
        "notes": ["CTMC+BSSVS 祖先地理重建", "Mantel 遗传-地理相关", "时序迁移",
                  "geo_paths: 事件级传播路径/episode/三级置信 + 扩散统计 (2026-09-16)"],
    },
]

# ─────────────────────────────────────────────────────────────
# 可渲染产物白名单 (E1, 2026-09-16)
# ─────────────────────────────────────────────────────────────
# 历史坑: collect_plots 只认 `.suffix == ".png"`, 而管线地理/系统地理类产物
# 几乎全是 PDF (phylogeo_multipanel.pdf / mantel_test.pdf / tempmig.pdf /
# rrt_results.pdf / VirSpaceTime.pdf / sample_map.pdf / mcmc_diagnostics.pdf …)
# → 报告 Geography 段长期「有产物但空图」且不报任何错。
# 同一条 `.png` 过滤还把 MODULES['data'] 里写的 `ncbi_harvest/*.pdf` 变成死 glob。
#
# 现在按扩展名分派: 图片 → <img>, PDF → <embed>+下载链接 (浏览器原生 PDF 阅读器)。
PLOT_IMG_EXTS = (".png", ".jpg", ".jpeg", ".gif", ".svg", ".webp")
PLOT_PDF_EXTS = (".pdf",)
# 自包含交互图 (phylogeo_map.html / skyline.html …): 体积可达数 MB (内联 plotly.min.js),
# 不做内联, 报告里给**相对路径链接**——报告与产物同在 work_dir 下一棵树内, 链接有效。
PLOT_LINK_EXTS = (".html",)
PLOT_EXTS = PLOT_IMG_EXTS + PLOT_PDF_EXTS

# 同一 stem 多格式共存时选谁内联 (png 最通用; svg 次之; pdf 兜底, 体积最大; html 只做链接)
_EXT_PRIORITY = {".png": 0, ".svg": 1, ".jpg": 2, ".jpeg": 3,
                 ".gif": 4, ".webp": 5, ".pdf": 6, ".html": 9}
_MIME = {".png": "image/png", ".jpg": "image/jpeg", ".jpeg": "image/jpeg",
         ".gif": "image/gif", ".svg": "image/svg+xml", ".webp": "image/webp",
         ".pdf": "application/pdf"}

# 内联体积上限 (base64 膨胀 ~33%, 报告是单文件, 不能无限膨胀)
MAX_IMG_KB = 500
MAX_PDF_KB = 6000


# ─────────────────────────────────────────────────────────────
# 工具函数
# ─────────────────────────────────────────────────────────────

def embed_asset(path: str,
                max_img_kb: int = MAX_IMG_KB,
                max_pdf_kb: int = MAX_PDF_KB) -> Tuple[str, str, str]:
    """产物 → (kind, data_uri, reason)。

    kind ∈ {'img', 'pdf', ''}；不内联时 data_uri 为空且 **reason 非空** —— 
    调用方据此在报告里显式写出「未内联原因」，不静默丢产物。
    """
    if not path:
        return "", "", "路径为空"
    p = Path(path)
    if not p.is_file():
        return "", "", "文件不存在"
    ext = p.suffix.lower()
    if ext not in PLOT_EXTS:
        return "", "", f"扩展名 {ext or '(无)'} 不在可渲染白名单 {PLOT_EXTS}"
    kind = "pdf" if ext in PLOT_PDF_EXTS else "img"
    size_kb = p.stat().st_size / 1024
    cap = max_pdf_kb if kind == "pdf" else max_img_kb
    if size_kb > cap:
        return kind, "", f"{size_kb:.0f} KB > 内联上限 {cap} KB"
    try:
        with open(path, "rb") as f:
            b64 = base64.b64encode(f.read()).decode()
    except OSError as e:
        return kind, "", f"读取失败: {e}"
    return kind, f"data:{_MIME[ext]};base64,{b64}", ""


def _rel_href(path: str, output_dir: Optional[str]) -> str:
    if output_dir:
        try:
            return os.path.relpath(path, output_dir).replace("\\", "/")
        except ValueError:
            pass
    return Path(path).name


def asset_card(path: str, output_dir: Optional[str] = None) -> Tuple[str, str]:
    """产物 → (HTML 卡片, 跳过原因)。不内联时仍给相对路径下载链接。"""
    name = Path(path).name
    label = _html.escape(Path(path).stem.replace("_", " ").title())
    if Path(path).suffix.lower() in PLOT_LINK_EXTS:
        href = _html.escape(_rel_href(path, output_dir), quote=True)
        return (f'<div class="img-card"><div class="img-label">{label}'
                f'<span class="tag">HTML 交互图</span></div>'
                f'<p class="asset-note">↗ <a href="{href}">在新标签打开 {_html.escape(name, quote=True)}</a>'
                f'（单文件自包含, 零外链）</p></div>'), ""
    kind, uri, reason = embed_asset(path)
    if uri:
        if kind == "pdf":
            return (f'<div class="img-card"><div class="img-label">{label}'
                    f'<span class="tag">PDF</span></div>'
                    f'<embed class="pdf-embed" type="application/pdf" src="{uri}" /></div>'), ""
        return (f'<div class="img-card"><div class="img-label">{label}</div>'
                f'<img src="{uri}" /></div>'), ""
    href = _html.escape(_rel_href(path, output_dir), quote=True)
    return (f'<div class="img-card off"><div class="img-label">{label}'
            f'<span class="tag warn">未内联</span></div>'
            f'<p class="asset-note">⚠ {_html.escape(reason, quote=True)}<br>'
            f'<a href="{href}" download>下载 {_html.escape(name, quote=True)}</a></p></div>'), reason


def collect_plots(vdir: str, patterns: List[str], max_imgs: int = 20,
                  skipped: Optional[List[Tuple[str, str]]] = None) -> List[str]:
    """按 glob 收集模块可渲染产物 (png/svg/jpg/gif/webp/pdf/html), 同 stem 取优先级最高者。

    skipped: 传入 list 时把每个「被排除的候选 + 原因」追加进去 (上限之外、
            同 stem 低优先级、扩展名不在白名单), 由调用方写进报告 —— 排除不等于没发生。
    """
    by_stem: Dict[Tuple[str, str], Path] = {}
    for pat in patterns:
        for p in sorted(Path(vdir).glob(pat)):
            if not p.is_file():
                continue
            ext = p.suffix.lower()
            if ext not in PLOT_EXTS and ext not in PLOT_LINK_EXTS:
                if skipped is not None:
                    skipped.append((str(p), f"扩展名 {ext or '(无)'} 不在可渲染白名单"))
                continue
            key = (p.parent.as_posix(), p.stem.lower())
            cur = by_stem.get(key)
            if cur is None:
                by_stem[key] = p
            elif _EXT_PRIORITY[ext] < _EXT_PRIORITY[cur.suffix.lower()]:
                if skipped is not None:
                    skipped.append((str(cur), f"同名单格式 {ext} 优先级更高, 不重复嵌入"))
                by_stem[key] = p
            elif skipped is not None:
                skipped.append((str(p), f"同名单格式 {cur.suffix.lower()} 优先级更高, 不重复嵌入"))
    out = [str(p) for p in sorted(by_stem.values(),
                                  key=lambda x: (x.parent.as_posix(), x.stem.lower()))]
    if len(out) > max_imgs:
        if skipped is not None:
            for _p in out[max_imgs:]:
                skipped.append((_p, f"超出单模块上限 max_imgs={max_imgs}"))
        out = out[:max_imgs]
    return out


def missing_declared_summaries(vdir: str, m: Dict, stages_done) -> List[str]:
    """模块 stage 已完成、但声明产物在磁盘上找不到的条目 (E2/G3 复发护栏)。

    历史坑: report_builder 声明 `geo_analysis/spatiotemporal_counts.csv`, 产出端
    若写在别的目录, 报告只会「静默少一张表」。此处把「声明了但没找到」显式列出。
    """
    if not vdir or not stages_done:
        return []
    if not any(s in stages_done for s in m["stages"]):
        return []
    missing = []
    for sp in m["summaries"]:
        if any(ch in sp for ch in "*?["):
            hit = list(Path(vdir).glob(f"{m['dir']}/{sp}"))
        else:
            hit = [p for p in [os.path.join(vdir, m["dir"], sp)] if os.path.exists(p)]
        if not hit:
            missing.append(f"{m['dir']}/{sp}")
    return missing


def _fmt_cell(v) -> str:
    s = str(v)
    # 科学计数/长小数紧凑化
    if len(s) > 12:
        try:
            f = float(s)
            s = f"{f:.4g}"
        except ValueError:
            pass
    return s


def render_table(path: str, max_rows: int = 50) -> str:
    """csv/tsv → HTML 表 (最多 max_rows 行)"""
    if not path or not os.path.exists(path):
        return ""
    ext = Path(path).suffix.lower()
    try:
        if ext in (".tsv", ".txt"):
            rows = list(csv.reader(open(path, encoding="utf-8", errors="replace"), delimiter="\t"))
        else:
            rows = list(csv.reader(open(path, encoding="utf-8", errors="replace")))
    except Exception:
        return ""
    if not rows:
        return ""
    # 去全空行
    rows = [r for r in rows if any(c.strip() for c in r)]
    header, body = rows[0], rows[1:1 + max_rows]
    thead = "".join(f"<th>{_fmt_cell(c)}</th>" for c in header)
    trs = []
    for r in body:
        trs.append("<tr>" + "".join(f"<td>{_fmt_cell(c)}</td>" for c in r) + "</tr>")
    extra = f'<p class="tbl-more">… {len(rows)-1-max_rows} more rows (see source file)</p>' if len(rows) - 1 > max_rows else ""
    return (f'<table class="mod-table"><thead><tr>{thead}</tr></thead>'
            f'<tbody>{"".join(trs)}</tbody></table>{extra}')


def render_json(path: str, max_rows: int = 40) -> str:
    """json → key/value 两列表 (嵌套 dict 展开)"""
    if not path or not os.path.exists(path):
        return ""
    try:
        data = json.load(open(path, encoding="utf-8", errors="replace"))
    except Exception:
        return ""
    if isinstance(data, dict):
        items = list(data.items())[:max_rows]
        trs = "".join(f'<tr><td class="k">{_fmt_cell(k)}</td><td>{_fmt_cell(v)}</td></tr>'
                      for k, v in items)
        return f'<table class="mod-table kv"><thead><tr><th>Parameter</th><th>Value</th></tr></thead><tbody>{trs}</tbody></table>'
    if isinstance(data, list):
        return render_kv_list(data, max_rows)
    return ""


def render_kv_list(items: List, max_rows: int = 30) -> str:
    rows = []
    for it in items[:max_rows]:
        if isinstance(it, dict):
            rows.append("<tr>" + "".join(f"<td>{_fmt_cell(v)}</td>" for v in it.values()) + "</tr>")
    if not rows:
        return ""
    cols = list(items[0].keys())
    thead = "".join(f"<th>{_fmt_cell(c)}</th>" for c in cols)
    return f'<table class="mod-table"><thead><tr>{thead}</tr></thead><tbody>{"".join(rows)}</tbody></table>'


def render_md(path: str, max_chars: int = 3000) -> str:
    """md/txt → <pre> 简单展示"""
    if not path or not os.path.exists(path):
        return ""
    try:
        text = open(path, encoding="utf-8", errors="replace").read(max_chars)
    except Exception:
        return ""
    text = text.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
    return f'<pre class="md-pre">{text}</pre>'


def stage_badge(stage: str, r: Dict) -> str:
    vdir = r.get("virus_dir", "")
    done = stage in r.get("stages_done", [])
    skipped = stage in r.get("stages_skipped", [])
    failed = stage in r.get("stages_failed", [])
    # 本轮未涉及该 stage 时回退看 checkpoint (report 单独重跑场景)
    if not (done or skipped or failed) and vdir:
        if os.path.exists(os.path.join(vdir, ".checkpoints", f"{stage}.done")):
            done = True
    if done:
        return f'<span class="badge ok">{stage}</span>'
    if skipped:
        return f'<span class="badge skip">{stage} · skip</span>'
    if failed:
        return f'<span class="badge fail">{stage} · fail</span>'
    return f'<span class="badge off">{stage}</span>'


# ─────────────────────────────────────────────────────────────
# 主函数
# ─────────────────────────────────────────────────────────────

def build_html_report(
    results: List[Dict],
    output_dir: str,
    title: str = "PhyloPipeline Report",
) -> str:
    """
    从分析结果列表生成模块化 HTML 汇总报告。

    results 每个元素包含 virus, phase, success, outputs, metrics,
    stages_done, stages_skipped, stages_failed, 以及病毒 work_dir (virus_dir)。
    """
    summary_rows = []
    virus_sections = []

    for i, r in enumerate(results):
        vname = r.get("virus", f"Virus_{i}")
        metrics = r.get("metrics", {})
        success = r.get("success", False)
        vdir = r.get("virus_dir") or r.get("output_dir") or ""

        # ── 总表行: 加 Mantel / TMRCA / stage 计数 ──
        # 2026-09-15 (P0-A): 头条指标统一走 headline_metrics, 与 phylo_summary.csv
        # 完全同源; 每个值后面带来源标签, 读者不必猜是 treedater / PIC / RTT。
        hm = headline_metrics(metrics)
        mantel = metrics.get("mantel_r", "")
        tmrca_cell = fmt_metric(hm["tmrca"], hm["tmrca_source"])
        n_done = len(r.get("stages_done", []))
        n_skip = len(r.get("stages_skipped", []))
        n_fail = len(r.get("stages_failed", []))
        # report 单独重跑: 本轮 stages 空时以 checkpoint 为准
        if not (n_done or n_skip or n_fail) and vdir:
            cp_dir = os.path.join(vdir, ".checkpoints")
            if os.path.isdir(cp_dir):
                n_done = len([f for f in os.listdir(cp_dir) if f.endswith(".done")])
        overall = "PASS" if (n_fail == 0 and (n_done > 0 or success)) else ("FAIL" if n_fail else "PARTIAL")
        summary_rows.append(f"""
        <tr>
            <td>{vname}</td>
            <td><span class="badge {'ok' if overall == 'PASS' else ('fail' if overall == 'FAIL' else 'skip')}">{overall}</span></td>
            <td>{n_done} / {n_skip} / {n_fail}</td>
            <td>{fmt_metric(hm['rate'], hm['rate_source'])}</td>
            <td>{fmt_metric(hm['r2'], hm['r2_source'])}</td>
            <td>{tmrca_cell}</td>
            <td>{mantel if mantel != '' else '—'}</td>
            <td>{metrics.get("n_pathways", '') if metrics.get("n_pathways", '') != '' else '—'}</td>
            <td>{metrics.get("geo_mean_velocity_km_yr", '') if metrics.get("geo_mean_velocity_km_yr", '') != '' else '—'}</td>
        </tr>""")

        # ── 每模块一节 ──
        mod_html = []
        for m in MODULES:
            mdir = os.path.join(vdir, m["dir"]) if vdir else ""
            badges = "".join(stage_badge(s, r) for s in m["stages"])
            has_dir = mdir and os.path.isdir(mdir)

            # summary 内容
            tables = []
            for sp in m["summaries"]:
                matched = []
                if any(ch in sp for ch in "*?["):
                    if vdir:
                        matched = [str(p) for p in sorted(Path(vdir).glob(f"{m['dir']}/{sp}"))][:5]
                else:
                    sp_full = os.path.join(vdir, m["dir"], sp) if vdir else ""
                    if sp_full and os.path.exists(sp_full):
                        matched = [sp_full]
                for sp_full in matched:
                    ext = Path(sp_full).suffix.lower()
                    label = os.path.relpath(sp_full, os.path.join(vdir, m["dir"])) if vdir else Path(sp).name
                    if ext == ".json":
                        t = render_json(sp_full)
                    elif ext == ".md":
                        t = render_md(sp_full)
                    else:
                        t = render_table(sp_full)
                    if t:
                        tables.append(f'<div class="tbl-wrap"><div class="tbl-title">{label}</div>{t}</div>')

            # plots (E1: png/svg/jpg/gif/webp/pdf; 同 stem 取高优先级; 被排除的显式列出)
            plot_skips = []
            imgs = collect_plots(mdir, m["plots"], skipped=plot_skips) if has_dir else []
            img_tags, not_inlined = [], []
            for p in imgs:
                card, reason = asset_card(p, output_dir)
                img_tags.append(card)
                if reason:
                    not_inlined.append((p, reason))

            # E2: stage 已完成但声明的 summary 产物不在磁盘 → 显式列出 (G3 复发护栏)
            missing = missing_declared_summaries(vdir, m, r.get("stages_done", []))
            warn_html = ""
            if missing:
                warn_html += ('<p class="asset-note">⚠ 声明但未找到产物 (产出路径与报告 glob '
                              '可能不一致): ' + _html.escape("、".join(missing)) + '</p>')
            if not_inlined:
                warn_html += ('<p class="asset-note">⚠ 有产物未内联: ' + _html.escape(
                    "；".join(f"{Path(p).name} ({why})" for p, why in not_inlined)) + '</p>')
            if plot_skips:
                warn_html += ('<p class="asset-note">· 同模块未嵌入的候选 ('
                              f'{len(plot_skips)} 项): ' + _html.escape(
                    "；".join(f"{Path(p).name} ({why})" for p, why in plot_skips[:8])
                    + ("…" if len(plot_skips) > 8 else "")) + '</p>')

            # 未运行且无目录 → 整节灰显
            state_cls = "mod" if (has_dir or tables or img_tags) else "mod off"
            notes = "".join(f'<li>{n}</li>' for n in m.get("notes", []))
            body = ""
            if tables:
                body += "".join(tables)
            if img_tags:
                body += f'<div class="img-grid">{"".join(img_tags)}</div>'
            if warn_html:
                body += warn_html
            if not body:
                body = f'<p class="mod-empty">{"未运行 / 不适用 (skip)" if not has_dir else "该模块无表格与图片产物"}</p>'
            mod_html.append(f"""
            <div class="{state_cls}">
                <div class="mod-head">
                    <h4>{m['name']}</h4>
                    <div class="badges">{badges}</div>
                </div>
                {f'<ul class="mod-notes">{notes}</ul>' if notes and not tables else ''}
                {body}
            </div>""")

        virus_sections.append(f"""
        <div class="virus-card">
            <h3>{vname} <span class="badge {'ok' if overall == 'PASS' else ('fail' if overall == 'FAIL' else 'skip')}">{overall}</span></h3>
            <div class="metrics">
                <span>β={fmt_metric(hm['rate'], hm['rate_source'])}</span>
                <span>R²={fmt_metric(hm['r2'], hm['r2_source'])}</span>
                <span>P={fmt_metric(hm['p'], hm['p_source'])}</span>
                <span>TMRCA={fmt_metric(hm['tmrca'], hm['tmrca_source'])}</span>
            </div>
            {''.join(mod_html)}
        </div>""")

    # Full HTML
    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>{title}</title>
<style>
    * {{ margin: 0; padding: 0; box-sizing: border-box; }}
    body {{ font-family: -apple-system, 'Segoe UI', Arial, sans-serif; background: #f5f6fa; color: #2c3e50; padding: 20px; }}
    h1 {{ text-align: center; margin: 20px 0; font-size: 28px; color: #1a1a2e; }}
    .summary-table {{ width: 100%; max-width: 1100px; margin: 20px auto; border-collapse: collapse; background: white; border-radius: 8px; overflow: hidden; box-shadow: 0 2px 8px rgba(0,0,0,0.08); }}
    .summary-table th {{ background: #1a1a2e; color: white; padding: 12px 16px; text-align: left; font-size: 13px; }}
    .summary-table td {{ padding: 10px 16px; border-bottom: 1px solid #eee; font-size: 13px; }}
    .summary-table tr:hover {{ background: #f8f9ff; }}
    .virus-card {{ background: white; border-radius: 10px; padding: 24px; margin: 20px auto; max-width: 1200px; box-shadow: 0 2px 12px rgba(0,0,0,0.06); }}
    .virus-card h3 {{ font-size: 18px; margin-bottom: 12px; display: flex; align-items: center; gap: 10px; }}
    .badge {{ font-size: 11px; padding: 3px 10px; border-radius: 12px; font-weight: 600; white-space: nowrap; }}
    .badge.ok {{ background: #27ae6018; color: #27ae60; }}
    .badge.fail {{ background: #e74c3c18; color: #e74c3c; }}
    .badge.skip {{ background: #f39c1218; color: #b9770e; }}
    .badge.off {{ background: #95a5a622; color: #7f8c8d; }}
    .metrics {{ display: flex; gap: 16px; margin-bottom: 16px; font-size: 13px; color: #666; }}
    .metrics span {{ background: #f0f0f5; padding: 4px 10px; border-radius: 6px; font-family: monospace; }}
    .mod {{ background: #fafafc; border: 1px solid #e8e8ef; border-radius: 8px; padding: 16px; margin: 14px 0; }}
    .mod.off {{ opacity: 0.55; }}
    .mod-head {{ display: flex; justify-content: space-between; align-items: center; flex-wrap: wrap; gap: 8px; margin-bottom: 8px; }}
    .mod-head h4 {{ font-size: 15px; color: #1a1a2e; }}
    .badges {{ display: flex; gap: 6px; flex-wrap: wrap; }}
    .mod-notes {{ margin: 4px 0 8px 18px; font-size: 12px; color: #888; }}
    .mod-empty {{ font-size: 13px; color: #aaa; font-style: italic; padding: 6px 0; }}
    .tbl-wrap {{ margin: 10px 0; }}
    .tbl-title {{ font-size: 12px; font-weight: 600; color: #555; margin-bottom: 4px; font-family: monospace; }}
    .mod-table {{ border-collapse: collapse; width: 100%; max-height: 400px; display: block; overflow: auto; background: white; }}
    .mod-table th {{ background: #2c3e50; color: white; padding: 6px 10px; font-size: 12px; text-align: left; position: sticky; top: 0; }}
    .mod-table td {{ padding: 5px 10px; border-bottom: 1px solid #eee; font-size: 12px; font-family: monospace; }}
    .mod-table.kv td.k {{ font-weight: 600; width: 40%; }}
    .tbl-more {{ font-size: 11px; color: #999; padding: 4px; }}
    .md-pre {{ background: white; border: 1px solid #eee; border-radius: 6px; padding: 12px; font-size: 12px; max-height: 300px; overflow: auto; white-space: pre-wrap; }}
    .img-grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(320px, 1fr)); gap: 16px; margin-top: 12px; }}
    .img-card {{ background: white; border-radius: 8px; overflow: hidden; border: 1px solid #eee; }}
    .img-card.off {{ border-style: dashed; }}
    .img-label {{ font-size: 12px; padding: 8px 12px; background: #1a1a2e; color: white; font-weight: 500; display: flex; justify-content: space-between; align-items: center; gap: 8px; }}
    .img-card img {{ width: 100%; display: block; }}
    .pdf-embed {{ width: 100%; height: 460px; display: block; background: #f0f0f5; }}
    .tag {{ font-size: 10px; padding: 2px 7px; border-radius: 9px; background: #ffffff2e; font-weight: 600; }}
    .tag.warn {{ background: #e74c3c; }}
    .asset-note {{ font-size: 11px; color: #b9770e; background: #fdf6e3; border-left: 3px solid #f39c12; padding: 6px 10px; margin-top: 8px; line-height: 1.6; }}
    .asset-note a {{ color: #2c3e50; }}
    .footer {{ text-align: center; margin-top: 40px; padding: 20px; color: #999; font-size: 12px; }}
</style>
</head>
<body>
<h1>{title}</h1>
<p style="text-align:center;color:#888;margin-bottom:30px">
    Generated: {datetime.now().strftime('%Y-%m-%d %H:%M')} · MMPV-RNA PhyloPipeline · VirPhyKit
</p>

<h2 style="max-width:1100px;margin:0 auto 10px;font-size:16px;">Summary</h2>
<table class="summary-table">
    <thead>
        <tr><th>Virus</th><th>Status</th><th>done/skip/fail</th><th>β (rate)</th><th>R²</th><th>TMRCA</th><th>Mantel r</th><th>Geo paths</th><th>Geo velocity</th></tr>
    </thead>
    <tbody>
        {''.join(summary_rows)}
    </tbody>
</table>
<p style="max-width:1100px;margin:8px auto 0;color:#888;font-size:12px;">
    括号内为指标来源。优先级: 速率 td_mean_rate(treedater) → pic_rate(PIC)；
    R²  pic_r(PIC) → r_squared(RTT)；P  td_relaxed_clock_p(treedater) → p_value(RTT)；
    TMRCA  beast_tmrca(BEAST) → td_tmrca(treedater)。
    与 <code>phylo_summary.csv</code> 同源同值（该表另存 <code>*_source</code> 与 <code>alt_*</code> 备用口径列）。
</p>

<h2 style="max-width:1200px;margin:30px auto 10px;font-size:16px;">Per-virus Results by Module</h2>
{''.join(virus_sections) if virus_sections else '<p style="text-align:center;color:#999;">No results to display</p>'}

<div class="footer">MMPV-RNA PhyloPipeline · Powered by VirPhyKit (Yin &amp; Gao, 2025)</div>
</body>
</html>"""

    report_path = os.path.join(output_dir, "phylo_report.html")
    with open(report_path, 'w', encoding='utf-8') as f:
        f.write(html)

    return report_path
