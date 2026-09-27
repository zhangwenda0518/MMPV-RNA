#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
eve_report.py — EVE 筛查结果 HTML 报告生成器
==============================================
在 eve_screen.py --merge 之后运行, 把汇总 TSV 渲染成自包含 HTML 报告,
补齐 ⑥ 管线与 ①-⑤ 一致的"末阶段出 HTML 报告"闭环。

输入 (全部为 --merge 产物, 均可缺省, 缺哪个跳过对应章节):
  <outdir>/kingdom_summary.tsv        genome, viral_supported, host_like, undetermined
  <outdir>/family_by_genome.tsv       family, <genome1>, <genome2>, ...
  <outdir>/kingdom_loci_all.tsv.gz    genome, locus, verdict, ref_family, ..., rvdb_bitscore (13 列)

输出:
  <outdir>/04_Summary/eve_report.html (默认; --out 可覆盖)

用法:
  python eve_report.py -o out/eve/
  python eve_report.py -o out/eve/ --top 30

纯标准库实现 (无 pandas/matplotlib 依赖), 服务器/本机均可直接跑。
"""

import argparse
import csv
import gzip
import os
import sys
from datetime import datetime
from html import escape
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

# 与筛查编排器共用同一目录注册表 (mmpv_common.io_layout, 带独立运行兜底)
try:
    from eve_scan_core import STAGE_DIRS  # noqa: E402
except Exception:  # pragma: no cover - 脱离仓库环境兜底
    STAGE_DIRS = {"loci": "01_Loci", "verdict": "02_Verdict",
                  "rvdb": "03_RVDB", "summary": "04_Summary"}

VERDICTS = ("viral_supported", "host_like", "undetermined")
# 展示名统一用短英文标识 (与 TSV 字段可对照); 中文语义统一放底部"判定口径"说明
VERDICT_LABEL = {"viral_supported": "Viral-supported (内源病毒化石)",
                 "host_like": "Host-like (宿主基因)",
                 "undetermined": "Undetermined (未判定)"}
VERDICT_COLOR = {"viral_supported": "#2b7a3d",
                 "host_like": "#b58900",
                 "undetermined": "#6b7280"}

_CSS = """
body{font-family:'Segoe UI','Microsoft YaHei',Helvetica,Arial,sans-serif;margin:24px;
     color:#222;background:#fafafa;line-height:1.45}
h1{font-size:22px;border-bottom:3px solid #2b7a3d;padding-bottom:8px}
h2{font-size:16px;margin-top:28px;color:#2b7a3d}
table{border-collapse:collapse;margin:10px 0;font-size:13px;background:#fff}
th,td{border:1px solid #ccc;padding:4px 10px;text-align:left}
th{background:#eef4ee;position:sticky;top:0}
tr:nth-child(even){background:#f6f8f6}
td.num{text-align:right;font-variant-numeric:tabular-nums}
.meta{color:#666;font-size:12px;margin:4px 0}
.cards{display:flex;gap:14px;margin:14px 0;flex-wrap:wrap}
.card{background:#fff;border:1px solid #ddd;border-radius:8px;padding:12px 18px;
      flex:1 1 0;min-width:150px;text-align:center}
.card .v{font-size:26px;font-weight:700}
.bar-track{background:#e8e8e8;border-radius:3px;width:100%;height:16px;overflow:hidden}
.bar-fill{height:16px;border-radius:3px 0 0 3px;min-width:2px}
.note{background:#fff8e6;border-left:4px solid #b58900;padding:8px 12px;
      font-size:12px;margin:10px 0}
"""


def _num(x):
    try:
        return int(x)
    except (TypeError, ValueError):
        return 0


def load_kingdom_summary(outdir: Path):
    """→ [(genome, viral_supported, host_like, undetermined), ...] 或 None。"""
    p = outdir / "kingdom_summary.tsv"
    if not p.is_file():
        return None
    rows = []
    with open(p, encoding="utf-8") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            rows.append((r.get("genome", ""), _num(r.get("viral_supported")),
                         _num(r.get("host_like")), _num(r.get("undetermined"))))
    return rows


def load_family_by_genome(outdir: Path):
    """→ [(family, {genome: count}, total), ...] 按 total 降序, 或 None。"""
    p = outdir / "family_by_genome.tsv"
    if not p.is_file():
        return None
    out = []
    with open(p, encoding="utf-8") as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader, None)
        genomes = header[1:] if header else []
        for row in reader:
            if not row:
                continue
            fam, counts = row[0], [_num(x) for x in row[1:]]
            out.append((fam, dict(zip(genomes, counts)), sum(counts)))
    out.sort(key=lambda x: -x[2])
    return out


def load_loci_all(outdir: Path):
    """→ verdict×family 计数, (family,bitscore,top) 摘要, 总行数; 文件缺省返回 None。"""
    p = outdir / "kingdom_loci_all.tsv.gz"
    if not p.is_file():
        return None
    vf = {}
    fam_tot = {}
    top = []
    n = 0
    with gzip.open(p, "rt", encoding="utf-8") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            n += 1
            v = r.get("verdict", "")
            fam = r.get("ref_family", "") or "(no_family)"
            vf[(v, fam)] = vf.get((v, fam), 0) + 1
            fam_tot[fam] = fam_tot.get(fam, 0) + 1
            try:
                bs = float(r.get("ref_bitscore") or 0)
            except ValueError:
                bs = 0.0
            top.append((bs, r.get("genome", ""), r.get("locus", ""),
                        v, fam, r.get("rvdb_stitle", "")))
    top.sort(key=lambda x: -x[0])
    return vf, fam_tot, top, n


def render(outdir: Path, top_n: int) -> str:
    """汇总 → 自包含 HTML 字符串。"""
    ks = load_kingdom_summary(outdir)
    fg = load_family_by_genome(outdir)
    la = load_loci_all(outdir)

    n_genomes = len(ks) if ks else 0
    tot = {v: sum(r[i + 1] for r in ks) for i, v in enumerate(VERDICTS)} if ks else {}
    n_loci = sum(tot.values())
    now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    h = ["<!DOCTYPE html><html><head><meta charset='utf-8'>",
         "<title>MMPV-EVE 筛查报告</title>",
         f"<style>{_CSS}</style></head><body>",
         "<h1>MMPV ⑥ EVE 内源性病毒筛查报告</h1>",
         f"<p class='meta'>生成时间: {escape(now)} | 输入根: {escape(str(outdir))} | "
         f"数据源: kingdom_summary.tsv / family_by_genome.tsv / kingdom_loci_all.tsv.gz</p>"]

    # ── 总览卡片 ──
    h.append("<div class='cards'>")
    h.append(f"<div class='card'><div class='v'>{n_genomes}</div>"
             f"<div class='meta'>基因组 (screened)</div></div>")
    h.append(f"<div class='card'><div class='v'>{n_loci}</div>"
             f"<div class='meta'>EVE 位点 (loci)</div></div>")
    for v in VERDICTS:
        pct = (100.0 * tot[v] / n_loci) if n_loci else 0.0
        h.append(f"<div class='card'><div class='v' style='color:{VERDICT_COLOR[v]}'>"
                 f"{tot.get(v, 0)}</div><div class='meta'>{VERDICT_LABEL[v]} "
                 f"({pct:.1f}%)</div></div>")
    h.append("</div>")

    if ks is None:
        h.append("<div class='note'>未找到 kingdom_summary.tsv —— 请先运行 "
                 "<code>eve_screen.py --merge</code> 生成汇总, 再重跑本报告。</div>")

    # ── 判定分布条 (固定宽 track 内按占比渲染, 避免 % 宽度在 auto 单元格中塌缩) ──
    if n_loci:
        h.append("<h2>位点判定分布</h2><table><tr><th>verdict</th><th>数量</th>"
                 "<th>占比</th><th style='width:300px'>占比条 (相对全部 loci)</th></tr>")
        for v in VERDICTS:
            pct = 100.0 * tot.get(v, 0) / n_loci
            h.append(f"<tr><td>{VERDICT_LABEL[v]}</td>"
                     f"<td class='num'>{tot.get(v, 0)}</td>"
                     f"<td class='num'>{pct:.1f}%</td>"
                     f"<td><div class='bar-track'>"
                     f"<div class='bar-fill' style='width:{max(pct, 0.8):.1f}%;"
                     f"background:{VERDICT_COLOR[v]}'></div></div></td></tr>")
        h.append("</table>")

    # ── 逐基因组表 ──
    if ks:
        h.append("<h2>逐基因组汇总 (按 viral_supported 降序)</h2>"
                 "<table><tr><th>genome</th>"
                 f"<th>{VERDICT_LABEL['viral_supported'].split(' (')[0]}</th>"
                 f"<th>{VERDICT_LABEL['host_like'].split(' (')[0]}</th>"
                 f"<th>{VERDICT_LABEL['undetermined'].split(' (')[0]}</th>"
                 "<th>loci 合计</th></tr>")
        for g, nv, nh, nu in sorted(ks, key=lambda r: (-r[1], r[0])):
            h.append(f"<tr><td>{escape(g)}</td><td class='num'>{nv}</td>"
                     f"<td class='num'>{nh}</td><td class='num'>{nu}</td>"
                     f"<td class='num'>{nv + nh + nu}</td></tr>")
        h.append("</table>")

    # ── 家族维度 ──
    if fg:
        h.append(f"<h2>病毒科 × 基因组 (Top {min(top_n, len(fg))} / 共 {len(fg)} 科)</h2>"
                 "<table><tr><th>family</th><th>loci 合计</th><th>各基因组分布</th></tr>")
        for fam, per_g, total in fg[:top_n]:
            dist = ", ".join(f"{escape(g)}:{c}" for g, c in sorted(per_g.items())
                             if c > 0) or "-"
            h.append(f"<tr><td>{escape(fam)}</td><td class='num'>{total}</td>"
                     f"<td>{dist}</td></tr>")
        h.append("</table>")

    # ── verdict × family 交叉 + Top 位点 ──
    if la:
        vf, fam_tot, top, n_all = la
        n_fam = len(fam_tot)
        h.append(f"<h2>ref_family × verdict 交叉 (Top {min(top_n, n_fam)} / 共 {n_fam} 科)</h2>"
                 "<table><tr><th>family</th>")
        for v in VERDICTS:
            h.append(f"<th>{VERDICT_LABEL[v].split(' (')[0]}</th>")
        h.append("<th>合计</th></tr>")
        for fam in sorted(fam_tot, key=lambda k: -fam_tot[k])[:top_n]:
            cells = "".join(f"<td class='num'>{vf.get((v, fam), 0)}</td>"
                            for v in VERDICTS)
            h.append(f"<tr><td>{escape(fam)}</td>{cells}"
                     f"<td class='num'>{fam_tot[fam]}</td></tr>")
        h.append("</table>")

        h.append(f"<h2>Top {min(top_n, len(top))} 位点 (按 ref_bitscore)</h2>"
                 "<table><tr><th>bitscore</th><th>genome</th><th>locus</th>"
                 "<th>verdict</th><th>family</th><th>RVDB 最佳命中</th></tr>")
        for bs, g, locus, v, fam, stitle in top[:top_n]:
            color = VERDICT_COLOR.get(v, "#000")
            h.append(f"<tr><td class='num'>{bs:.0f}</td><td>{escape(g)}</td>"
                     f"<td>{escape(locus)}</td>"
                     f"<td style='color:{color}'>{VERDICT_LABEL[v].split(' (')[0]}</td>"
                     f"<td>{escape(fam)}</td><td>{escape((stitle or '')[:80])}</td></tr>")
        h.append(f"</table><p class='meta'>loci_all 总行数: {n_all}</p>")

    h.append("<div class='note'>判定口径: Viral-supported=病毒侧 bitscore 占优"
             "(内源病毒化石); Host-like=植物/宿主侧占优; Undetermined=双侧均不占优。"
             "生成工具: endogenous_virus_pipeline/eve_report.py (MMPV ⑥ 管线)。</div>")
    h.append("</body></html>")
    return "\n".join(h)


def main():
    ap = argparse.ArgumentParser(
        description="EVE 筛查汇总 → 自包含 HTML 报告 (在 eve_screen --merge 之后运行)")
    ap.add_argument("-o", "--outdir", required=True,
                    help="EVE 筛查输出根目录 (含 kingdom_summary.tsv 等 --merge 产物)")
    ap.add_argument("--out", default=None,
                    help=f"输出 HTML 路径 (默认: <outdir>/{STAGE_DIRS['summary']}/eve_report.html)")
    ap.add_argument("--top", type=int, default=20,
                    help="家族/位点明细表展示条数 (默认: 20)")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    if not outdir.is_dir():
        sys.exit(f"[ERROR] 输出目录不存在: {outdir}")

    html = render(outdir, args.top)
    out_path = Path(args.out) if args.out else outdir / STAGE_DIRS["summary"] / "eve_report.html"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", encoding="utf-8") as f:
        f.write(html)
    print(f"[OK] EVE 报告 → {out_path} ({len(html) / 1024:.1f} KB)")


if __name__ == "__main__":
    main()
