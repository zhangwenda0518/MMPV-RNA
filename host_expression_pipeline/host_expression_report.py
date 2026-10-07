#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""host_expression_report.py — HTML 汇总报告 (⑧ 宿主表达模块, 09_Report 阶段)
=========================================================================
把各阶段产物汇总为一页 HTML (纯标准库, 无外部资源; 火山图 base64 内嵌,
文件可直接发送/归档)。任何一段产物缺失时该段自动跳过并在页面标注。

用法:
  python host_expression_report.py --root <08_HostExpression> \
      --output <08_HostExpression>/09_Report/host_expression_report.html
"""

import argparse
import base64
import html
import json
from pathlib import Path

CSS = """
body{font-family:-apple-system,'Segoe UI','Microsoft YaHei',sans-serif;
     margin:0;background:#f5f6f8;color:#222;line-height:1.55}
.wrap{max-width:1080px;margin:0 auto;padding:28px 20px 60px}
h1{font-size:22px;margin:0 0 4px} h2{font-size:16px;margin:28px 0 10px;
   border-left:4px solid #2f6fb3;padding-left:8px}
.sub{color:#667;font-size:13px;margin-bottom:18px}
.cards{display:flex;gap:12px;flex-wrap:wrap;margin:14px 0}
.card{background:#fff;border:1px solid #e2e5ea;border-radius:8px;
      padding:12px 18px;min-width:130px}
.card .v{font-size:22px;font-weight:600;color:#2f6fb3}
.card .k{font-size:12px;color:#667}
table{border-collapse:collapse;background:#fff;width:100%;font-size:13px}
th,td{border:1px solid #e2e5ea;padding:5px 8px;text-align:left}
th{background:#eef2f7} tr:nth-child(even) td{background:#fafbfc}
.fail{color:#c0392b;font-weight:600}
.miss{color:#b58b00}
img{max-width:100%;border:1px solid #e2e5ea;border-radius:6px;background:#fff}
.note{font-size:12px;color:#667;margin-top:6px}
"""


def _read_tsv(path: Path):
    if not path.is_file():
        return None
    with path.open(encoding="utf-8") as fh:
        rows = [l.rstrip("\n").split("\t") for l in fh if l.strip()]
    return rows if rows else None


def _fmt_cell(c: str) -> str:
    """长浮点数显示为 4 位有效数字 (样品名/基因名等非数值原样保留)。"""
    if len(c) < 9:
        return c
    try:
        v = float(c)
    except (TypeError, ValueError):
        return c
    return f"{v:.4g}"


def _table(rows, highlight_col=None, fail_pred=None) -> str:
    if not rows:
        return "<p class='miss'>（无此段产物）</p>"
    head, body = rows[0], rows[1:]
    out = ["<table><thead><tr>"]
    out += [f"<th>{html.escape(c)}</th>" for c in head]
    out.append("</tr></thead><tbody>")
    for r in body:
        out.append("<tr>")
        for i, c in enumerate(r):
            cls = ""
            if fail_pred and i == highlight_col and fail_pred(c):
                cls = " class='fail'"
            out.append(f"<td{cls}>{html.escape(_fmt_cell(c))}</td>")
        out.append("</tr>")
    out.append("</tbody></table>")
    return "".join(out)


def _img_b64(path: Path, alt: str) -> str:
    if not path.is_file():
        return ""
    b64 = base64.b64encode(path.read_bytes()).decode("ascii")
    return (f"<p><img alt='{html.escape(alt)}' "
            f"src='data:image/png;base64,{b64}'></p>")


def build_report(root: Path) -> str:
    qc_rows = _read_tsv(root / "03_QC" / "qc_summary.tsv")
    srows = _read_tsv(root / "04_Labels" / "sample_labels.tsv")
    vrows = _read_tsv(root / "04_Labels" / "virus_columns.tsv")
    vlabels = _read_tsv(root / "04_Labels" / "virus_labels.tsv")
    coexp_sum = _read_tsv(root / "05_Coexp" / "diff_coexp_summary.tsv")
    coexp_hub = _read_tsv(root / "05_Coexp" / "diff_coexp_gene_aggregate.tsv")
    netqc = _read_tsv(root / "05_Coexp" / "network_qc.tsv")
    mod_sum = _read_tsv(root / "05_Coexp" / "modules" / "module_summary.tsv")
    mod_assoc = _read_tsv(root / "05_Coexp" / "modules" / "module_load_assoc.tsv")
    de_sum = _read_tsv(root / "06_DE" / "deseq2_summary.tsv")
    degs = _read_tsv(root / "06_DE" / "deseq2_degs.tsv")
    scan = _read_tsv(root / "08_LoadScan" / "load_scan_summary.tsv")
    enrich_sum = _read_tsv(root / "07_Enrich" / "enrich_summary.tsv")
    enrich_files = sorted((root / "07_Enrich").glob("enrich_*.tsv")
                          if (root / "07_Enrich").is_dir() else [])

    run_json = root / "02_Quant" / "metadata" / "run_summary.json"
    run_meta = {}
    if run_json.is_file():
        try:
            run_meta = json.loads(run_json.read_text(encoding="utf-8"))
        except json.JSONDecodeError:
            pass

    n_samples = max(len(srows) - 1, 0) if srows else 0
    n_pos = n_neg = 0
    if srows:
        hdr = {c: i for i, c in enumerate(srows[0])}
        if "Virus_any" in hdr:
            j = hdr["Virus_any"]
            n_pos = sum(1 for r in srows[1:] if len(r) > j and r[j] == "+")
            n_neg = len(srows) - 1 - n_pos
    n_deg = (len(degs) - 1) if degs else 0
    n_pairs = "-"
    if coexp_sum and len(coexp_sum) > 1:
        try:
            j = coexp_sum[0].index("n_unique_pairs")
            n_pairs = coexp_sum[1][j]
        except (ValueError, IndexError):
            pass
    n_mod = (len(mod_sum) - 1) if mod_sum else 0
    n_quant = run_meta.get("samples_total")
    n_pe = run_meta.get("paired_end_samples")
    n_se = run_meta.get("single_end_samples")
    if n_quant is None:                      # 缺 run_summary.json: 回退数 quant 目录
        qdir = root / "02_Quant" / "quant"
        n_quant = (sum(1 for p in qdir.iterdir() if p.is_dir())
                   if qdir.is_dir() else "?")
        n_pe = n_se = "?"

    parts = [
        "<!DOCTYPE html><html lang='zh'><head><meta charset='utf-8'>",
        "<title>MMPV 宿主表达报告</title>",
        f"<style>{CSS}</style></head><body><div class='wrap'>",
        "<h1>⑧ 宿主表达分析报告（可选补充模块）</h1>",
        f"<div class='sub'>输出根: {html.escape(str(root))}</div>",
        "<div class='cards'>",
        f"<div class='card'><div class='v'>{n_quant}</div>"
        f"<div class='k'>定量样本 (PE {n_pe} / SE {n_se})</div></div>",
        f"<div class='card'><div class='v'>{n_samples}</div>"
        f"<div class='k'>标签样本</div></div>",
        f"<div class='card'><div class='v'>{n_pos} / {n_neg}</div>"
        f"<div class='k'>病毒阳性 / 阴性</div></div>",
        f"<div class='card'><div class='v'>{n_pairs}</div>"
        f"<div class='k'>Δz 显著对 (coexp)</div></div>",
        f"<div class='card'><div class='v'>{n_mod}</div>"
        f"<div class='k'>共表达模块 (z≥1.7∩DEG)</div></div>",
        f"<div class='card'><div class='v'>{n_deg}</div>"
        f"<div class='k'>显著 DEG (交叉参考)</div></div>",
        "</div>",
        "<h2>样本定量 QC (03_QC)</h2>",
        _table(qc_rows, highlight_col=4, fail_pred=lambda c: c == "0"),
        "<div class='note'>pass=0 的样本 (低比对率/低 reads) 已在 DESeq2 "
        "阶段剔除（若启用了 QC 过滤）。</div>",
        "<h2>病毒检出与分组 (04_Labels)</h2>",
        _table(vrows),
    ]

    if vlabels and len(vlabels) > 1:
        top = sorted((r for r in vlabels[1:] if len(r) >= 7),
                     key=lambda r: -_safe_float(r[6]))[:15]
        parts.append("<h3>检出明细 (按覆盖度前 15 行)</h3>")
        parts.append(_table([vlabels[0]] + top))

    parts.append("<h2>差异共表达 TEA-GCN Δz (05_Coexp, 主分析)</h2>")
    parts.append(_table(coexp_sum))
    if coexp_hub and len(coexp_hub) > 1:
        parts.append("<h3>Hub 基因 (显著对数 Top 15)</h3>")
        parts.append(_table([coexp_hub[0]] + coexp_hub[1:16]))
        parts.append(_img_b64(root / "09_Report" / "figures" /
                              "fig01_hub_ranking.png", "hub ranking"))
    parts.append("<div class='note'>阈值: 候选对 Δz 采样分布 0.5%/99.5% 分位; "
                 "有向行数与唯一对数两种口径见 summary。<b>解读警示</b>: 候选对"
                 "本身是每基因 |Δz| Top-N 极值, 相对全总体分位线天然偏松 (最大值"
                 "统计效应), 显著对数量是阈值装置的函数而非生物学计数; 请基于 "
                 "hub 排序、|Δz| 排名与下游模块/富集解读, 必要时以置换零模型校准"
                 "(组内打乱样本标签重建网络的经验分布)。</div>")

    # 网络 z 分布 QC (KingdomStress Fig S8 闸门: 分布截断 → Δz/模块不可信)
    parts.append("<h2>网络 z 分布 QC (05_Coexp, Fig S8 闸门)</h2>")
    qc_status_i = netqc[0].index("Status") if netqc and "Status" in netqc[0] \
        else None
    parts.append(_table(netqc, highlight_col=qc_status_i,
                        fail_pred=lambda c: c in ("FAIL", "LOW")))
    parts.append("<div class='note'>PASS: max_z≥1.7 (模块检测可行); LOW: "
                 "1.5≤max_z&lt;1.7 (分布压缩); FAIL: max_z&lt;1.5 (Δz 跨网络比较"
                 "不可信)。判级依据 KingdomStress Fig S8: 样本量不足会整体压低 "
                 "MR z-score 分布。</div>")

    # DEG 限定共表达模块 (KingdomStress Fig5 流程)
    parts.append("<h2>DEG 限定共表达模块 (05_Coexp/modules, "
                 "z≥1.7 → Louvain res=1.0)</h2>")
    if mod_sum and len(mod_sum) > 1:
        body = mod_sum[1:]
        hdr = mod_sum[0]
        ng_i = hdr.index("N_Genes") if "N_Genes" in hdr else 2
        body = sorted(body, key=lambda r: -_safe_float(r[ng_i]))[:15]
        parts.append(_table([hdr] + body))
        parts.append("<div class='note'>按模块规模 Top 15; Hub_Genes 为模块内"
                     "度最高的前 3 个基因。模块划分基于单一分辨率 (res=1.0, "
                     "seed=42), 未做 bootstrap 稳定性评估; 跨分辨率一致性可"
                     "用 --module-resolution 扫描核验。</div>")
    else:
        parts.append("<p class='miss'>（无模块产物 — 需 06_DE 的 DEG 表与 "
                     "z≥1.7 网络; 06_DE 完成后自动补跑）</p>")
    if mod_assoc and len(mod_assoc) > 1:
        body = mod_assoc[1:]
        hdr = mod_assoc[0]
        pj = hdr.index("Padj") if "Padj" in hdr else len(hdr) - 1
        body = sorted(body, key=lambda r: _safe_float(r[pj]))[:15]
        parts.append("<h3>模块 × 病毒载量 Spearman 关联 (Padj 前 15)</h3>")
        parts.append(_table([hdr] + body))
    parts.append(_img_b64(root / "09_Report" / "figures" /
                          "fig02_module_load.png", "module load assoc"))
    parts.append(_img_b64(root / "09_Report" / "figures" /
                          "fig05_circular_pathway.png",
                          "circular pathway network"))
    for side in ("positive", "negative"):
        parts.append(_img_b64(
            root / "09_Report" / "figures" /
            f"fig04_module_network_{side}.png", f"module network {side}"))

    parts.append("<h2>差异表达 (06_DE, 交叉参考)</h2>")
    parts.append(_table(de_sum))
    if degs and len(degs) > 1:
        body = degs[1:]
        hdr = degs[0]
        padj_i = hdr.index("padj") if "padj" in hdr else len(hdr) - 1
        body = sorted(body, key=lambda r: _safe_float(r[padj_i]))[:20]
        parts.append("<h3>Top 20 DEG (按 padj)</h3>")
        parts.append(_table([hdr] + body))
    parts.append(_img_b64(root / "06_DE" / "volcano_plot.png", "volcano"))

    parts.append("<h2>病毒载量×全基因组扫描 (08_LoadScan, 剂量验证)</h2>")
    parts.append(_table(scan))

    # 功能富集 (07_Enrich): 汇总 + 每基因集显著 term Top 8
    parts.append("<h2>功能富集 ORA (07_Enrich, 超几何+BH)</h2>")
    parts.append(_table(enrich_sum))
    parts.append(_img_b64(root / "09_Report" / "figures" /
                          "fig03_enrichment.png", "enrichment dotplot"))
    for f in enrich_files:
        if f.name == "enrich_summary.tsv":
            continue
        rows = _read_tsv(f)
        if rows and len(rows) > 1:
            parts.append(f"<h3>{html.escape(f.stem.replace('enrich_', ''))} "
                         f"显著 term (Top 8)</h3>")
            parts.append(_table([rows[0]] + rows[1:9]))

    parts.append("<div class='note'>MMPV ⑧ 宿主表达模块 · "
                 "数据来源见各阶段目录 · 本页无外部依赖</div>")
    parts.append("</div></body></html>")
    return "".join(parts)


def _safe_float(s: str) -> float:
    try:
        return float(s)
    except (TypeError, ValueError):
        return float("inf")


def main() -> None:
    ap = argparse.ArgumentParser(description="⑧ 宿主表达 HTML 汇总报告")
    ap.add_argument("--root", required=True, help="08_HostExpression 输出根")
    ap.add_argument("--output", help="HTML 路径 (默认 <root>/09_Report/host_expression_report.html)")
    args = ap.parse_args()

    root = Path(args.root)
    out = Path(args.output) if args.output else \
        root / "09_Report" / "host_expression_report.html"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(build_report(root), encoding="utf-8")
    print(f"报告: {out}")


if __name__ == "__main__":
    main()
