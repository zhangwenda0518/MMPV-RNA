#!/usr/bin/env python3
"""
panvirome_novel_known_report.py — 新病毒判定报告生成器
====================================================
读 panvirome_novel_known.py 输出的 novel_vs_known_final.tsv，生成一个
自包含的 HTML 报告（深色 SCI 风格 + Chart.js），含:
  1. 判定分布柱状图
  2. 新病毒候选清单(可排序表格)
  3. 判定方法说明

用法:
  python panvirome_novel_known_report.py --tsv novel_vs_known_final.tsv \
      --out novel_vs_known_report.html
"""

import argparse
import csv
import json
from collections import Counter
from pathlib import Path

# 判定类别 -> 颜色（深色 SCI 配色）
JUDGE_COLORS = {
    "已知病毒": "#66bb6a",
    "已知病毒(蛋白)": "#a5d6a7",
    "新种候选": "#ffa726",
    "新种候选(蛋白)": "#ffcc80",
    "新属/远缘候选": "#ef5350",
    "远缘新病毒候选": "#ab47bc",
    "非病毒(噬菌体)": "#90a4ae",
    "非病毒(宿主)": "#78909c",
    "远缘未定": "#8d6e63",
    "未定": "#9e9e9e",
}

NOVEL_CLASSES = {"新种候选", "新种候选(蛋白)", "新属/远缘候选", "远缘新病毒候选"}

# 图内文字英文化: 仅用于图表标签显示 (数据列 final_judge 原值不改, 表格仍用原值)
JUDGE_EN = {
    "已知病毒": "Known virus",
    "已知病毒(蛋白)": "Known virus (protein)",
    "新种候选": "Novel species candidate",
    "新种候选(蛋白)": "Novel species candidate (protein)",
    "新属/远缘候选": "Novel genus / distant candidate",
    "远缘新病毒候选": "Distant novel virus candidate",
    "非病毒(噬菌体)": "Non-viral (phage)",
    "非病毒(宿主)": "Non-viral (host)",
    "远缘未定": "Distant, undetermined",
    "未定": "Undetermined",
}


def parse_args():
    p = argparse.ArgumentParser(description="新病毒判定报告生成器")
    p.add_argument("--tsv", required=True, help="novel_vs_known_final.tsv 路径")
    p.add_argument("--out", default="novel_vs_known_report.html", help="输出 HTML")
    p.add_argument("--title", default="Novel vs Known Virus", help="报告标题")
    return p.parse_args()


def esc(v):
    return (v or "").replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def fnum(v, default=0.0):
    """容错 float: novel_vs_known_final.tsv 用 "-" 占位表示无命中, 直接 float() 会崩。"""
    try:
        if v is None or str(v).strip() in ("", "-", "NA", "N/A", "nan", "."):
            return default
        return float(v)
    except (TypeError, ValueError):
        return default


def load_chart_js_inline():
    """离线自包含: 复用仓库内已有 Chart.js (virome_discovery_pipeline/utils/chart.min.js)。

    report_pipeline.py 用同一份 js 内联; 找不到时返回 "" -> 报告降级为表格, 绝不引 CDN。
    """
    cand = Path(__file__).resolve().parent.parent / "virome_discovery_pipeline" / "utils" / "chart.min.js"
    if cand.is_file():
        try:
            return cand.read_text(encoding="utf-8")
        except Exception:
            return ""
    return ""


def main():
    args = parse_args()
    rows = list(csv.DictReader(open(args.tsv, encoding="utf-8"), delimiter="\t"))

    judge_cnt = Counter(r["final_judge"] for r in rows)
    novel = [r for r in rows if r["final_judge"] in NOVEL_CLASSES]
    novel.sort(key=lambda r: -fnum(r.get("blast_pident")))

    # 判定分布 (图表标签走 JUDGE_EN 英文化; 颜色/计数仍按原值取)
    judge_keys = list(judge_cnt.keys())
    labels = [JUDGE_EN.get(k, k) for k in judge_keys]
    data = [judge_cnt[k] for k in judge_keys]
    colors = [JUDGE_COLORS.get(k, "#9e9e9e") for k in judge_keys]

    # 新病毒候选表格
    table_rows = []
    for r in novel:
        cid = r["contig_id"]
        short = cid.split("_clean")[0] if "_clean" in cid else cid[:24]
        table_rows.append(
            f"<tr>"
            f"<td>{esc(short)}</td>"
            f"<td>{esc(r['final_judge'])}</td>"
            f"<td>{esc(r.get('aa_species','')[:40])}</td>"
            f"<td>{esc(r.get('cdd_top','')[:48])}</td>"
            f"<td>{esc(r.get('blast_pident','-'))}</td>"
            f"<td>{esc(r.get('blast_alnlen','-'))}</td>"
            f"<td>{esc(r.get('identify_tools','-'))}</td>"
            f"<td>{esc(r.get('classify_tool','-'))}</td>"
            f"</tr>"
        )

    # 离线自包含: 内联仓库自带 Chart.js; 缺失则降级为表格 (不引 CDN, 无网络时不再白图)
    chart_js = load_chart_js_inline()
    if chart_js:
        head_js = "<script>\n" + chart_js + "\n</script>"
        dist_html = '<div class="chartbox"><canvas id="dist"></canvas></div>'
        _opts = {
            "type": "bar",
            "data": {"labels": labels,
                     "datasets": [{"label": "vOTUs", "data": data, "backgroundColor": colors}]},
            "options": {"indexAxis": "y", "responsive": True, "plugins": {"legend": {"display": False}},
                        "scales": {"x": {"ticks": {"color": "#9aa0a6"}, "grid": {"color": "#262c35"}},
                                   "y": {"ticks": {"color": "#e8eaed"}, "grid": {"color": "#262c35"}}}},
        }
        chart_block = ("<script>\nnew Chart(document.getElementById('dist'), "
                       + json.dumps(_opts, ensure_ascii=False) + ");\n</script>")
    else:
        head_js = ""
        _rows_html = "".join(f"<tr><td>{esc(k)}</td><td>{judge_cnt[k]}</td></tr>" for k in judge_keys)
        dist_html = ('<div class="note">Local inline Chart.js not found '
                     '(<code>virome_discovery_pipeline/utils/chart.min.js</code>); '
                     'offline runs do not load the CDN, so the verdict distribution is shown as a table instead:</div>'
                     '<table style="max-width:520px"><thead><tr><th>Verdict</th><th>Count</th></tr></thead>'
                     f'<tbody>{_rows_html}</tbody></table>')
        chart_block = ""

    html = f"""<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8">
<title>{esc(args.title)}</title>
{head_js}
<style>
:root {{ --bg:#111418; --card:#1a1f26; --ink:#e8eaed; --ink2:#9aa0a6; --accent:#42a5f5; }}
body {{ background:var(--bg); color:var(--ink); font-family:-apple-system,'Segoe UI',Roboto,sans-serif; margin:0; padding:24px; }}
h1 {{ font-size:22px; font-weight:600; margin:0 0 4px; }}
.sub {{ color:var(--ink2); font-size:13px; margin-bottom:20px; }}
.card {{ background:var(--card); border:1px solid #262c35; border-radius:10px; padding:18px; margin-bottom:18px; }}
h2 {{ font-size:16px; margin:0 0 12px; color:var(--accent); }}
.chartbox {{ max-width:640px; }}
table {{ width:100%; border-collapse:collapse; font-size:12px; }}
th,td {{ text-align:left; padding:7px 8px; border-bottom:1px solid #262c35; }}
th {{ color:var(--ink2); font-weight:500; position:sticky; top:0; background:var(--card); }}
tr:hover td {{ background:#232a33; }}
.note {{ color:var(--ink2); font-size:12px; line-height:1.6; }}
</style></head><body>
<h1>{esc(args.title)}</h1>
<div class="sub">{len(rows)} rescue vOTUs &middot; {len(novel)} novel virus candidates</div>

<div class="card"><h2>Verdict Distribution</h2>
{dist_html}</div>

<div class="card"><h2>Novel Virus Candidate List ({len(novel)})</h2>
<div style="max-height:520px;overflow:auto">
<table>
<thead><tr><th>contig</th><th>Verdict</th><th>Protein hit</th><th>CDD domain</th><th>02blast%</th><th>alnlen</th><th>Identification tool</th><th>Classification tool</th></tr></thead>
<tbody>{''.join(table_rows)}</tbody>
</table></div></div>

<div class="card"><h2>Verdict Rules</h2>
<div class="note">
<b>Three-layer evidence</b>: nucleic acid identity (blastn vs nt virus DB) &rarr; protein identity (mmseqs vs ICTV protein DB) &rarr; CDD domains (distinguishing RdRp/RT from host domains).<br>
<b>Decision rules</b>: nt identity &ge;85% with coverage &ge;50% = known virus; 50-85% = novel species candidate; when there is no nt hit, use protein (&ge;70% known, 50-70% novel species); when protein is also weak, use CDD (RdRp/RT present = novel genus/distant candidate; host domains only = non-viral).<br>
<b>Provenance</b>: blast_* are the best hits of diamond vs the RVDB virus database at stage 02 (Identification); identify_tools are the tools that called it viral at stage 02; classify_tool is the primary_tool at stage 05 (Taxonomy).<br>
<b>Note</b>: CDD domains are hard evidence against taxonomy tools mislabelling a host gene as a virus (e.g. an extremely short blast alnlen plus a host enzyme in CDD = false hit).
</div></div>

{chart_block}
</body></html>"""

    with open(args.out, "w", encoding="utf-8") as f:
        f.write(html)
    print(f"输出: {args.out} ({len(rows)} 条, 新病毒候选 {len(novel)} 条)")


if __name__ == "__main__":
    main()
