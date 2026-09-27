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
from collections import Counter

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


def parse_args():
    p = argparse.ArgumentParser(description="新病毒判定报告生成器")
    p.add_argument("--tsv", required=True, help="novel_vs_known_final.tsv 路径")
    p.add_argument("--out", default="novel_vs_known_report.html", help="输出 HTML")
    p.add_argument("--title", default="Novel vs Known Virus", help="报告标题")
    return p.parse_args()


def esc(v):
    return (v or "").replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def main():
    args = parse_args()
    rows = list(csv.DictReader(open(args.tsv, encoding="utf-8"), delimiter="\t"))

    judge_cnt = Counter(r["final_judge"] for r in rows)
    novel = [r for r in rows if r["final_judge"] in NOVEL_CLASSES]
    novel.sort(key=lambda r: -(float(r["blast_pident"]) if r.get("blast_pident") else 0))

    # 判定分布
    labels = list(judge_cnt.keys())
    data = [judge_cnt[k] for k in labels]
    colors = [JUDGE_COLORS.get(k, "#9e9e9e") for k in labels]

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
            f"<td>{r.get('cdd_top','')[:48]}</td>"
            f"<td>{r.get('blast_pident','-')}</td>"
            f"<td>{r.get('blast_alnlen','-')}</td>"
            f"<td>{esc(r.get('identify_tools','-'))}</td>"
            f"<td>{esc(r.get('classify_tool','-'))}</td>"
            f"</tr>"
        )

    html = f"""<!DOCTYPE html>
<html lang="zh"><head><meta charset="utf-8">
<title>{esc(args.title)}</title>
<script src="https://cdn.jsdelivr.net/npm/chart.js"></script>
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
<div class="sub">共 {len(rows)} 条 rescue vOTU · 新病毒候选 {len(novel)} 条</div>

<div class="card"><h2>判定分布</h2>
<div class="chartbox"><canvas id="dist"></canvas></div></div>

<div class="card"><h2>新病毒候选清单（{len(novel)} 条）</h2>
<div style="max-height:520px;overflow:auto">
<table>
<thead><tr><th>contig</th><th>判定</th><th>蛋白命中</th><th>CDD 结构域</th><th>02blast%</th><th>alnlen</th><th>鉴定工具</th><th>分类工具</th></tr></thead>
<tbody>{''.join(table_rows)}</tbody>
</table></div></div>

<div class="card"><h2>判定方法</h2>
<div class="note">
<b>三层证据</b>：核酸 identity（blastn 对 nt 病毒库）→ 蛋白 identity（mmseqs 对 ICTV 蛋白库）→ CDD 结构域（区分 RdRp/RT vs 宿主结构域）。<br>
<b>判定规则</b>：核酸 identity ≥85% 且覆盖 ≥50% = 已知病毒；50~85% = 新种候选；核酸无命中时看蛋白（≥70% 已知，50~70% 新种）；蛋白也弱时看 CDD（有 RdRp/RT = 新属/远缘候选，仅宿主结构域 = 非病毒）。<br>
<b>溯源</b>：blast_* 为 02 鉴定阶段 diamond 对 RVDB 病毒库的 best-hit；identify_tools 为 02 阶段判成病毒的工具；classify_tool 为 05 分类阶段 primary_tool。<br>
<b>注意</b>：CDD 结构域是戳穿"分类工具误标宿主基因为病毒"的硬证据（例：blast alnlen 极短 + CDD 显示宿主酶 = 假命中）。
</div></div>

<script>
new Chart(document.getElementById('dist'), {{
  type:'bar',
  data:{{ labels:{labels}, datasets:[{{ label:'vOTUs', data:{data}, backgroundColor:{colors} }}] }},
  options:{{ indexAxis:'y', responsive:true, plugins:{{ legend:{{display:false}} }},
    scales:{{ x:{{ ticks:{{color:'#9aa0a6'}}, grid:{{color:'#262c35'}} }}, y:{{ ticks:{{color:'#e8eaed'}}, grid:{{color:'#262c35'}} }} }} }}
}});
</script>
</body></html>"""

    with open(args.out, "w", encoding="utf-8") as f:
        f.write(html)
    print(f"输出: {args.out} ({len(rows)} 条, 新病毒候选 {len(novel)} 条)")


if __name__ == "__main__":
    main()
