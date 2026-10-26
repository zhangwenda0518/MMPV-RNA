#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""eve_html.py — 02/03 阶段 HTML 汇总的共用渲染与表型识别.

为什么单独抽一层: 02_classify 和 03_report 都要出"阶段汇总 HTML", 但它们的输入
表来源不同 (有基因组 / 无基因组 / 只有 ML 分)。把"表头 -> 表型"的识别和版式放在
这里, 两个阶段脚本就只剩"收集哪几张表"这点差异, 不会各写一套颜色和列宽。

约定:
  * 纯标准库 (无 pandas), 服务器/本机都能直接跑
  * 表型由**表头**识别, 不靠文件名 —— 阶段产物被改名/搬目录后仍能正确汇总
  * 数值缺失/非数字不抛异常, 只当 0 跳过
"""

import csv
import gzip
import sys
from collections import Counter
from datetime import datetime
from html import escape
from pathlib import Path

_FALLBACK_CSS = """
:root{--success:#2b7a3d;--warning:#b58900;--warning-light:#fdf6e3;
      --danger:#b91c1c;--surface-alt:#f3f4f6;--ink-secondary:#6b7280}
body{font-family:-apple-system,'Segoe UI','Microsoft YaHei',sans-serif;
     margin:24px;color:#1f2937;background:#fff;max-width:1180px}
h1,h2,h3{color:#111827}table{border-collapse:collapse;margin:8px 0}
th,td{border:1px solid #d1d5db;padding:4px 10px;font-size:14px}
th{background:#f3f4f6}code,pre{background:#f9fafb;font-size:13px}
"""


def _find_repo_root():
    for base in Path(__file__).resolve().parents:
        if (base / "mmpv_common" / "report_style.py").is_file():
            return base
    return None


try:
    from mmpv_common.report_style import get_base_css
except ImportError:
    _root = _find_repo_root()
    if _root and str(_root) not in sys.path:
        sys.path.insert(0, str(_root))
    try:
        from mmpv_common.report_style import get_base_css
    except ImportError:
        def get_base_css():          # 独立运行兜底: 无 mmpv_common 也能出报告
            return _FALLBACK_CSS


_CSS = get_base_css() + """
h1{font-size:22px;border-bottom:3px solid var(--success);padding-bottom:8px}
h2{font-size:16px;margin-top:30px;color:var(--success);
   border-left:4px solid var(--success);padding-left:8px}
h3{font-size:14px;margin-top:18px}
table{border-collapse:collapse;margin:10px 0;font-size:13px;background:#fff}
th,td{border:1px solid #ccc;padding:4px 10px;text-align:left;vertical-align:top}
th{background:var(--surface-alt)}
tr:nth-child(even){background:#fafafa}
td.num{text-align:right;font-variant-numeric:tabular-nums}
.meta{color:var(--ink-secondary);font-size:12px;margin:4px 0}
.note{background:var(--warning-light);border-left:4px solid var(--warning);
      padding:8px 12px;font-size:12px;margin:12px 0}
.cards{display:flex;gap:12px;flex-wrap:wrap;margin:10px 0}
.card{border:1px solid #e5e7eb;border-radius:6px;padding:8px 14px;min-width:110px}
.card .k{font-size:11px;color:var(--ink-secondary)}
.card .v{font-size:18px;font-weight:600;font-variant-numeric:tabular-nums}
.bars{margin:8px 0}
.bar{display:grid;grid-template-columns:minmax(120px,240px) 1fr 90px;
     align-items:center;gap:8px;margin:2px 0;font-size:12px}
.bar .track{background:var(--surface-alt);height:12px;border-radius:3px;
            overflow:hidden}
.bar .fill{height:12px;background:var(--success)}
.bar .val{text-align:right;font-variant-numeric:tabular-nums;
          color:var(--ink-secondary)}
.missing{color:var(--danger);font-weight:600}
.ok{color:var(--success);font-weight:600}
"""

# 表型 -> (章节标题, 说明). 顺序即识别优先级, 越具体的列名越靠前.
KIND_INFO = {
    "final_tiers": ("最终分层 (final_tiers)", "多路证据合并后的 tier 分布"),
    "dbitscore":   ("Δbit-score 双库竞争", "Δ = EXO − EVE; Δ<0 偏内源, Δ>0 偏真病毒"),
    "ref_classify": ("双库分类 (ref_class)", "候选对 EVE 元件库 vs 真实病毒库的命中比较"),
    "refgenome":   ("参考基因组整合证据", "候选回比宿主染色体后的侧翼/剪接判据"),
    "chimera":     ("转录本嵌合", "宿主-病毒-宿主结构: 无基因组时的整合旁证"),
    "locus":       ("位点架构", "同一宿主 scaffold 上病毒组件是否齐全"),
    "nanovir":     ("Nanoviridae 无监督降解打分", "经验零分布的 q 值; 不走监督模型"),
    "ml":          ("ML 序列打分", "16 科分科模型; call 由阈值判定"),
    "ml_legacy":   ("ML 序列打分 (旧列名)", "p_eve_cls: 早期部署表的概率列"),
    "verdict":     ("verdict 表", "逐候选判决, 具体口径由上游脚本决定"),
    "generic":     ("其他表", "未识别为已知表型, 仅列概览"),
}


# ---------------------------------------------------------------- 读表

def read_table(path):
    """读 TSV -> (headers, rows). 支持 .gz; 空表/单列不报错."""
    p = Path(path)
    if not p.is_file():
        return [], []
    opener = gzip.open if p.suffix == ".gz" else open
    with opener(p, "rt", encoding="utf-8", errors="replace", newline="") as fh:
        rd = csv.reader(fh, delimiter="\t")
        try:
            headers = next(rd)
        except StopIteration:
            return [], []
        headers = [h.strip() for h in headers]
        rows = []
        for cells in rd:
            if not cells or (len(cells) == 1 and not cells[0].strip()):
                continue
            rows.append({h: (cells[i] if i < len(cells) else "")
                         for i, h in enumerate(headers)})
    return headers, rows


def kind_of(headers):
    """按表头判断表型; 返回 KIND_INFO 的 key."""
    hs = set(headers)
    if "final_tier" in hs:
        return "final_tiers"
    if "delta_bitscore" in hs:
        return "dbitscore"
    if "ref_class" in hs:
        return "ref_classify"
    if "rg_call" in hs:
        return "refgenome"
    if "chimera_call" in hs:
        return "chimera"
    if "locus_arch" in hs:
        return "locus"
    if "q_value" in hs or ("p_final" in hs and "z_stop" in hs):
        return "nanovir"
    if "score" in hs and "call" in hs:
        return "ml"
    if "p_eve_cls" in hs:
        return "ml_legacy"
    if "verdict" in hs:
        return "verdict"
    return "generic"


# ---------------------------------------------------------------- 小工具

def _f(x, default=None):
    try:
        v = float(x)
    except (TypeError, ValueError):
        return default
    if v != v:                       # NaN
        return default
    return v


def _fmt_int(n):
    return f"{int(n):,}"


def _fmt_pct(a, b):
    return f"{(100.0 * a / b):.1f}%" if b else "-"


def _median(vals):
    vals = sorted(vals)
    if not vals:
        return None
    n = len(vals)
    return vals[n // 2] if n % 2 else (vals[n // 2 - 1] + vals[n // 2]) / 2.0


def cards(items):
    """items: [(label, value)]"""
    out = ['<div class="cards">']
    for k, v in items:
        out.append(f'<div class="card"><div class="k">{escape(str(k))}</div>'
                   f'<div class="v">{escape(str(v))}</div></div>')
    out.append("</div>")
    return "".join(out)


def bars(counts, total=None, max_bars=14, color=None):
    """{label: n} -> 横向条形图 (纯 div/CSS)."""
    if not counts:
        return '<div class="meta">无可统计字段</div>'
    total = total or sum(counts.values())
    items = sorted(counts.items(), key=lambda kv: (-kv[1], str(kv[0])))[:max_bars]
    peak = max((n for _l, n in items), default=0) or 1
    style = f' style="background:{color}"' if color else ""
    out = ['<div class="bars">']
    for label, n in items:
        pct = 100.0 * n / peak
        out.append(
            '<div class="bar">'
            f'<div title="{escape(str(label))}">{escape(str(label) or "(空)")}</div>'
            f'<div class="track"><div class="fill" style="width:{pct:.1f}%"{style}></div></div>'
            f'<div class="val">{_fmt_int(n)} · {_fmt_pct(n, total)}</div>'
            "</div>")
    out.append("</div>")
    return "".join(out)


def hist_svg(values, bins=24, width=680, height=150, color="#2b7a3d"):
    """数值直方图 (内联 SVG, 不依赖 JS)."""
    vals = [v for v in values if v is not None]
    if not vals:
        return ""
    lo, hi = min(vals), max(vals)
    if hi == lo:
        hi = lo + 1.0
    counts = [0] * bins
    step = (hi - lo) / bins
    for v in vals:
        i = min(int((v - lo) / step), bins - 1)
        if i < 0:
            i = 0
        counts[i] += 1
    peak = max(counts) or 1
    pad_l, pad_b = 34, 20
    plot_w = width - pad_l - 8
    plot_h = height - pad_b - 8
    bw = plot_w / bins
    parts = [f'<svg viewBox="0 0 {width} {height}" width="100%" '
             f'style="max-width:{width}px;font:10px sans-serif">']
    for i, c in enumerate(counts):
        h = plot_h * c / peak
        x = pad_l + i * bw
        y = 8 + plot_h - h
        parts.append(f'<rect x="{x:.1f}" y="{y:.1f}" width="{max(bw - 1, 1):.1f}" '
                     f'height="{h:.1f}" fill="{color}" opacity="0.85"/>')
    parts.append(f'<line x1="{pad_l}" y1="{8 + plot_h}" '
                 f'x2="{width - 8}" y2="{8 + plot_h}" stroke="#9ca3af"/>')
    parts.append(f'<text x="2" y="{8 + plot_h}" fill="#6b7280">{peak}</text>')
    parts.append(f'<text x="2" y="14" fill="#6b7280">0</text>')
    parts.append(f'<text x="{pad_l}" y="{height - 4}" fill="#6b7280">'
                 f'{lo:.3g}</text>')
    parts.append(f'<text x="{width - 8}" y="{height - 4}" fill="#6b7280" '
                 f'text-anchor="end">{hi:.3g}</text>')
    parts.append("</svg>")
    return "".join(parts)


def table_html(headers, rows, limit=25, num_cols=()):
    if not headers:
        return '<div class="meta">空表</div>'
    num_cols = set(num_cols)
    out = ["<table><thead><tr>"]
    for h in headers:
        out.append(f"<th>{escape(h)}</th>")
    out.append("</tr></thead><tbody>")
    for r in rows[:limit]:
        out.append("<tr>")
        for h in headers:
            cls = ' class="num"' if h in num_cols else ""
            v = r.get(h, "")
            if h in num_cols:
                fv = _f(v)
                v = "" if fv is None else f"{fv:g}"
            out.append(f"<td{cls}>{escape(str(v))}</td>")
        out.append("</tr>")
    out.append("</tbody></table>")
    if len(rows) > limit:
        out.append(f'<div class="meta">仅显示前 {limit} / {_fmt_int(len(rows))} 行</div>')
    return "".join(out)


# ---------------------------------------------------------------- 各表型汇总

def _section_ml(rows, top):
    key = "score" if "score" in rows[0] else "p_eve_cls"
    call_col = "call" if "call" in rows[0] else None
    vals = [_f(r.get(key)) for r in rows]
    vals = [v for v in vals if v is not None]
    calls = Counter(r.get(call_col, "") for r in rows) if call_col else Counter()
    fam = Counter(r.get("family", "") for r in rows)
    c = cards([("候选数", _fmt_int(len(rows))),
               ("有分值", _fmt_int(len(vals))),
               ("中位分值", f"{_median(vals):.3f}" if vals else "-"),
               ("EVE_like", _fmt_int(calls.get("EVE_like", 0)))])
    html = [c]
    if calls:
        html.append("<h3>call 分布</h3>" + bars(calls, total=len(rows)))
    html.append("<h3>分值分布</h3>" + hist_svg(vals))
    if fam:
        html.append("<h3>科分布</h3>" + bars(fam, total=len(rows), max_bars=18))
    # 每科 call 率
    if call_col:
        per = {}
        for r in rows:
            f = r.get("family", "")
            d = per.setdefault(f, [0, 0])
            d[0] += 1
            if r.get(call_col) == "EVE_like":
                d[1] += 1
        prof = [(f, n, e, 100.0 * e / n if n else 0.0)
                for f, (n, e) in per.items()]
        prof.sort(key=lambda x: -x[1])
        html.append("<h3>分科统计</h3>")
        html.append(table_html(
            ["family", "n", "EVE_like", "EVE_like_pct"],
            [{"family": f, "n": n, "EVE_like": e, "EVE_like_pct": f"{p:.1f}"}
             for f, n, e, p in prof],
            limit=30, num_cols=("n", "EVE_like", "EVE_like_pct")))
        html.append("<h3>分值最高的候选</h3>")
        top_rows = sorted(rows, key=lambda r: -(_f(r.get(key), -1e9)))
        html.append(table_html(list(rows[0].keys()), top_rows, limit=top))
    return "".join(html)


def _section_dbitscore(rows, top):
    vals = [_f(r.get("delta_bitscore")) for r in rows]
    vals = [v for v in vals if v is not None]
    neg = sum(1 for v in vals if v < 0)
    pos = sum(1 for v in vals if v > 0)
    zero = len(vals) - neg - pos
    c = cards([("候选数", _fmt_int(len(rows))),
               ("Δ 可计算", _fmt_int(len(vals))),
               ("Δ<0 偏内源", f"{_fmt_int(neg)} ({_fmt_pct(neg, len(vals))})"),
               ("Δ>0 偏病毒", f"{_fmt_int(pos)} ({_fmt_pct(pos, len(vals))})"),
               ("Δ=0", _fmt_int(zero))])
    html = [c, "<h3>Δ bit-score 分布 (负=内源侧)</h3>" + hist_svg(vals)]
    fam = Counter(r.get("family", "") for r in rows)
    if fam:
        html.append("<h3>科分布</h3>" + bars(fam, total=len(rows), max_bars=18))
    # 分科 Δ 方向
    per = {}
    for r in rows:
        v = _f(r.get("delta_bitscore"))
        if v is None:
            continue
        d = per.setdefault(r.get("family", ""), [0, 0])
        d[0] += 1
        if v < 0:
            d[1] += 1
    if per:
        prof = [(f, n, e, 100.0 * e / n) for f, (n, e) in per.items()]
        prof.sort(key=lambda x: -x[1])
        html.append("<h3>分科 Δ<0 比例</h3>")
        html.append(table_html(
            ["family", "n", "delta_negative", "delta_negative_pct"],
            [{"family": f, "n": n, "delta_negative": e,
              "delta_negative_pct": f"{p:.1f}"} for f, n, e, p in prof],
            limit=30, num_cols=("n", "delta_negative", "delta_negative_pct")))
    html.append("<h3>Δ 最偏内源的候选</h3>")
    top_rows = sorted([r for r in rows if _f(r.get("delta_bitscore")) is not None],
                      key=lambda r: _f(r.get("delta_bitscore"), 0.0))
    html.append(table_html(list(rows[0].keys()), top_rows, limit=top,
                           num_cols=("delta_bitscore", "exo_bitscore",
                                     "eve_bitscore", "max_bitscore")))
    return "".join(html)


def _section_call_table(rows, top, col, num_cols=()):
    calls = Counter(r.get(col, "") for r in rows)
    html = [cards([("候选数", _fmt_int(len(rows))),
                   ("类别数", _fmt_int(len(calls)))])]
    html.append(f"<h3>{escape(col)} 分布</h3>" + bars(calls, total=len(rows)))
    html.append("<h3>明细</h3>")
    headers = list(rows[0].keys())
    html.append(table_html(headers, rows, limit=top, num_cols=num_cols))
    return "".join(html)


def _section_nanovir(rows, top):
    qs = [_f(r.get("q_value")) for r in rows]
    qs = [v for v in qs if v is not None]
    calls = Counter(r.get("call", "") for r in rows)
    html = [cards([("候选数", _fmt_int(len(rows))),
                   ("EVE_like", _fmt_int(calls.get("EVE_like", 0))),
                   ("中位 q", f"{_median(qs):.4g}" if qs else "-"),
                   ("q<0.05", _fmt_int(sum(1 for v in qs if v < 0.05)))])]
    html.append("<h3>call 分布</h3>" + bars(calls, total=len(rows)))
    html.append("<h3>q 值分布</h3>" + hist_svg(qs))
    html.append("<h3>明细</h3>")
    html.append(table_html(list(rows[0].keys()), rows, limit=top,
                           num_cols=("len_nt", "stop_density", "orf_max_frac",
                                     "gc", "z_stop", "z_orf", "z_gc", "stat",
                                     "p_final", "q_value", "p_mosaic",
                                     "p_drop2")))
    return "".join(html)


def _section_final_tiers(rows, top):
    """final_tiers: tier 分布 + 证据组合 + 明细."""
    order = ["tier1_endogenous_strong", "tier2_endogenous_weak",
             "tier3_no_evidence", "tier4_exogenous"]
    tiers = Counter(r.get("final_tier", "") for r in rows)
    c = cards([(k.replace("tier", "T"), _fmt_int(tiers.get(k, 0)))
               for k in order] + [("候选数", _fmt_int(len(rows)))])
    html = [c, "<h3>tier 分布</h3>" + bars(tiers, total=len(rows))]
    # 证据组合: 直接用上游写好的 evidence_list, 只取证据类型名
    combos = Counter()
    for r in rows:
        ev = r.get("evidence_list", "") or ""
        kinds = []
        for tag, name in (("命中已知EVE元件", "同源EVE"),
                          ("转录组嵌合", "嵌合"),
                          ("位点组件齐全", "位点齐全"),
                          ("多拷贝家族", "多拷贝"),
                          ("命中真实病毒", "命中原病毒")):
            if tag in ev:
                kinds.append(name)
        combos["+".join(kinds) if kinds else "无证据"] += 1
    html.append("<h3>证据组合</h3>" + bars(combos, total=len(rows), max_bars=16))
    html.append("<h3>候选明细 (按 tier 排序)</h3>")
    ranked = sorted(rows, key=lambda r: (order.index(r.get("final_tier", ""))
                                         if r.get("final_tier", "") in order
                                         else 99))
    html.append(table_html(list(rows[0].keys()), ranked, limit=top,
                           num_cols=("n_endo_evidence", "n_exo_evidence",
                                     "eve_bits", "virus_bits", "locus_ncomp",
                                     "multicopy_size", "bam_mean_depth")))
    return "".join(html)


def _section_generic(rows, top, headers):
    html = [cards([("行数", _fmt_int(len(rows))),
                   ("列数", _fmt_int(len(headers)))])]
    html.append(table_html(headers, rows, limit=top))
    return "".join(html)


def section_for(kind, rows, top=25):
    """按表型生成章节 HTML (不含 <h2>)."""
    if not rows:
        return '<div class="meta">空表</div>'
    if kind in ("ml", "ml_legacy"):
        return _section_ml(rows, top)
    if kind == "final_tiers":
        return _section_final_tiers(rows, top)
    if kind == "dbitscore":
        return _section_dbitscore(rows, top)
    if kind == "nanovir":
        return _section_nanovir(rows, top)
    if kind == "ref_classify":
        return _section_call_table(rows, top, "ref_class",
                                   ("eve_bits", "virus_bits", "delta"))
    if kind == "refgenome":
        return _section_call_table(rows, top, "rg_call",
                                   ("flank_gap", "host_seg_cov", "best_host_pid"))
    if kind == "chimera":
        return _section_call_table(rows, top, "chimera_call",
                                   ("left_cov", "right_cov", "mid_gap",
                                    "host_identity"))
    if kind == "locus":
        return _section_call_table(rows, top, "locus_arch",
                                   ("locus_ncomp", "locus_members", "locus_span"))
    if kind == "verdict":
        return _section_call_table(rows, top, "verdict")
    return _section_generic(rows, top, list(rows[0].keys()))


def section(kind, title, rows, top=25, extra_note=""):
    body = section_for(kind, rows, top=top)
    note = f'<div class="note">{extra_note}</div>' if extra_note else ""
    return f"<h2>{escape(title)}</h2>{note}{body}"


# ---------------------------------------------------------------- 页面

def render_page(title, subtitle, sections, footer=""):
    stamp = datetime.now().strftime("%Y-%m-%d %H:%M")
    return f"""<!DOCTYPE html>
<html lang="zh-CN"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{escape(title)}</title><style>{_CSS}</style></head>
<body>
<h1>{escape(title)}</h1>
<div class="meta">{subtitle}</div>
{''.join(sections)}
<div class="meta" style="margin-top:36px;border-top:1px solid #e5e7eb;padding-top:8px">
生成时间 {stamp} · virome_discovery_pipeline/eve_discriminate{(' · ' + footer) if footer else ''}
</div>
</body></html>
"""


def write_page(path, html):
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(html, encoding="utf-8")
    return p


def overview_table(entries):
    """entries: [(name, path, kind_key, n_rows)] -> 概览表 HTML."""
    rows = []
    for name, path, kind, n in entries:
        label = KIND_INFO.get(kind, KIND_INFO["generic"])[0]
        rows.append({"table": name, "type": label, "rows": _fmt_int(n),
                     "path": str(path)})
    return table_html(["table", "type", "rows", "path"], rows,
                      limit=100, num_cols=("rows",))
