#!/usr/bin/env python3
"""
preprocess_report.py — 数据清洗管道 report 阶段 (汇总 + 交接清单)
==================================================================
汇聚 clean-data.py (fastp JSON) 与 host_depletion.py (seqkit 汇总表 /
ribodetector 报告 / 比对日志) 的既有统计产物, 生成:

  preprocessing_summary.tsv   每样本各步 reads 数 / 留存率 / 状态
  assembly_ready.list         最终 reads 文件清单 (下游 --input_reads 交接)
  preprocessing_report.html   交互式 HTML 报告 (Chart.js 内嵌, 离线可用)

纯解析不重算: 全部输入是上游已落盘的统计文件, 缺失项留空, 不报错。

用法:
  python preprocess_report.py --output-dir out/
  python preprocess_report.py --clean-dir out/00a_CleanData --deplete-dir out/00b_HostDepletion
  python data_preprocessing.py --stage all ...   # report 阶段自动收尾
"""

import argparse
import os
import csv
import json
import re
import sys
from pathlib import Path

# 跨管线统一 I/O 布局 (mmpv_common/, 仓库根): 目录名随 MMPV_IO_LAYOUT 解析
_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)
from mmpv_common.io_layout import layout_dirs as _layout_dirs
_D = _layout_dirs(os.environ.get("MMPV_IO_LAYOUT", "legacy"))

SCRIPT_DIR = Path(__file__).resolve().parent

if sys.platform == 'win32':
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')


def _num(s):
    if s in (None, ""):
        return None
    try:
        return int(str(s).replace(",", "").strip())
    except ValueError:
        try:
            return float(str(s).replace(",", "").strip())
        except ValueError:
            return None


# ═══════════════════════════════════════════
# 收集: 各上游统计文件 → {sample: {...}}
# ═══════════════════════════════════════════

def collect_fastp(clean_dir):
    """clean-data.py 落盘的每样本 fastp JSON: logs/{sample}_fastp_report.json"""
    out = {}
    if not clean_dir:
        return out
    for p in sorted(Path(clean_dir).rglob("*_fastp_report.json")):
        sample = p.name[: -len("_fastp_report.json")]
        rec = {"raw": None, "fastp": None, "q30": None, "dup": None, "path": str(p)}
        try:
            d = json.loads(p.read_text(encoding="utf-8"))
            summ = d.get("summary", {})
            rec["raw"] = _num((summ.get("before_filtering") or {}).get("total_reads"))
            rec["fastp"] = _num((summ.get("after_filtering") or {}).get("total_reads"))
            q30 = (summ.get("after_filtering") or {}).get("q30_rate")
            rec["q30"] = round(float(q30) * 100, 1) if q30 is not None else None
            dup = (d.get("duplication") or {}).get("rate")
            rec["dup"] = round(float(dup) * 100, 1) if dup is not None else None
        except Exception:
            pass
        out[sample] = rec
    return out


def collect_seqkit(deplete_dir):
    """host_depletion.py 的 {logs}/host_depletion_seqkit_summary.tsv (Sample/Stage/num_seqs...)"""
    out = {}
    if not deplete_dir:
        return out
    for p in Path(deplete_dir).rglob("host_depletion_seqkit_summary.tsv"):
        try:
            with open(p, encoding="utf-8") as f:
                rows = list(csv.DictReader(f, delimiter="\t"))
        except Exception:
            continue
        for r in rows:
            sample = (r.get("Sample") or "").strip()
            stage = (r.get("Stage") or "").strip()
            key = None
            if stage.startswith("1_"):
                key = "input"
            elif stage.startswith("2_"):
                key = "kraken2"
            elif stage.startswith("3_"):
                key = "host"
            elif stage.startswith("4_"):
                key = "rrna"
            if key and sample:
                n = _num(r.get("num_seqs"))
                if n is not None:
                    out.setdefault(sample, {})[key] = n
    return out


def collect_rrna(deplete_dir):
    """ribodetector 汇总报告 (--rrna_report, 默认 ribodetector.report.txt)"""
    out = {}
    if not deplete_dir:
        return out
    for p in Path(deplete_dir).rglob("*.txt"):
        name = p.name.lower()
        if "ribodetector" not in name and "rrna" not in name:
            continue
        if "kraken2" in name:
            continue
        try:
            with open(p, encoding="utf-8") as f:
                for r in csv.DictReader(f, delimiter="\t"):
                    sample = (r.get("Sample") or "").strip()
                    if not sample:
                        continue
                    out[sample] = {
                        "total": _num(r.get("Total_sequences")),
                        "non_rrna": _num(r.get("non_rRNA")),
                        "rrna": _num(r.get("rRNA")),
                    }
        except Exception:
            continue
    return out


def collect_aligner(deplete_dir):
    """bowtie2/hisat2 stderr 汇总行 (从未被解析过的盲点): 'overall alignment rate  95.5%'"""
    out = {}
    if not deplete_dir:
        return out
    pat = re.compile(r"overall alignment rate\D+([\d.]+)", re.IGNORECASE)
    for p in Path(deplete_dir).rglob("*.step2_align_filter.*.log"):
        sample = p.name.split(".")[0]
        try:
            m = pat.search(p.read_text(encoding="utf-8", errors="replace"))
        except Exception:
            continue
        if m:
            try:
                out[sample] = float(m.group(1))
            except ValueError:
                pass
    return out


FINAL_PATTERNS = [
    (re.compile(r"^(?P<s>.+)_clean_1\.(?:fastq|fq|fa)\.gz$"), "pe"),
    (re.compile(r"^(?P<s>.+)_clean_2\.(?:fastq|fq|fa)\.gz$"), "pe"),
    (re.compile(r"^(?P<s>.+)_clean\.(?:fastq|fq|fa)\.gz$"), "se"),
    (re.compile(r"^(?P<s>.+)_1\.(?:fastq|fq|fa)\.gz$"), "pe"),
    (re.compile(r"^(?P<s>.+)_2\.(?:fastq|fq|fa)\.gz$"), "pe"),
    (re.compile(r"^(?P<s>.+)\.(?:fastq|fq|fa)\.gz$"), "se"),
]


def collect_final_files(deplete_dir, clean_dir):
    """最终 reads 文件: 去宿主输出优先, 回退 clean 的 3.clumpify / 2.fasta"""
    files = {}
    roots = []
    if deplete_dir:
        roots.append(Path(deplete_dir))
    if clean_dir:
        for sub in ("3.clumpify", "2.fasta"):
            d = Path(clean_dir) / sub
            if d.is_dir():
                roots.append(d)
    for root in roots:
        for f in sorted(root.rglob("*.gz")):
            if not f.is_file():
                continue
            for pat, layout in FINAL_PATTERNS:
                m = pat.match(f.name)
                if m:
                    files.setdefault(m.group("s"), []).append(f)
                    break
    return files


# ═══════════════════════════════════════════
# 汇总
# ═══════════════════════════════════════════

def build_rows(clean_dir, deplete_dir, warn_retained):
    fastp = collect_fastp(clean_dir)
    seqkit = collect_seqkit(deplete_dir)
    rrna = collect_rrna(deplete_dir)
    aligner = collect_aligner(deplete_dir)
    files = collect_final_files(deplete_dir, clean_dir)

    samples = sorted(set(fastp) | set(seqkit) | set(rrna) | set(files))
    rows = []
    for s in samples:
        fp = fastp.get(s, {})
        sk = seqkit.get(s, {})
        rr = rrna.get(s, {})
        input_n = sk.get("input")
        after_k2 = sk.get("kraken2")
        after_host = sk.get("host")
        after_rrna = sk.get("rrna") if sk.get("rrna") is not None else rr.get("non_rrna")
        final = files.get(s, [])
        final_n = len(final)

        base = input_n if input_n is not None else fp.get("fastp")
        retained_n = next((v for v in (after_rrna, after_host, after_k2) if v is not None), None)
        retained = round(retained_n / base * 100, 1) if (retained_n is not None and base) else None

        if retained is None:
            status = "NO_DATA" if not final else "OK"
        elif retained < warn_retained:
            status = "LOW_RETAINED"
        else:
            status = "PASS"
        if not final:
            status = "NO_OUTPUT" if status != "NO_DATA" else status

        rows.append({
            "Sample": s,
            "Raw_Reads": fp.get("raw"),
            "Fastp_Clean_Reads": fp.get("fastp"),
            "Fastp_Q30(%)": fp.get("q30"),
            "Depletion_Input": input_n,
            "After_Kraken2": after_k2,
            "After_HostAlign": after_host,
            "After_rRNA": after_rrna,
            "Retained(%)": retained,
            "Aligner_Rate(%)": aligner.get(s),
            "Final_Files": final_n,
            "Status": status,
        })
    return rows


def write_summary_tsv(rows, out_path):
    cols = ["Sample", "Raw_Reads", "Fastp_Clean_Reads", "Fastp_Q30(%)", "Depletion_Input",
            "After_Kraken2", "After_HostAlign", "After_rRNA", "Retained(%)",
            "Aligner_Rate(%)", "Final_Files", "Status"]
    with open(out_path, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)
        # TOTAL 行 (数值列求和, 留存率按输入加权)
        tot = {c: 0 for c in cols[1:-2] if c not in ("Fastp_Q30(%)", "Retained(%)", "Aligner_Rate(%)")}
        wsum = wbase = 0
        for r in rows:
            for c in tot:
                v = r.get(c)
                if isinstance(v, (int, float)):
                    tot[c] += v
            if r.get("Retained(%)") is not None and r.get("Depletion_Input"):
                wsum += r["Retained(%)"] / 100 * r["Depletion_Input"]
                wbase += r["Depletion_Input"]
        total = {"Sample": "TOTAL", **tot,
                 "Retained(%)": round(wsum / wbase * 100, 1) if wbase else None,
                 "Status": ""}
        w.writerow(total)


def write_ready_list(files, out_path):
    """交接清单: 每行一个最终 reads 文件绝对路径 (按样本分组)"""
    with open(out_path, "w", encoding="utf-8", newline="") as f:
        for s in sorted(files):
            for pth in files[s]:
                f.write(str(pth.resolve()) + "\n")


# ═══════════════════════════════════════════
# HTML
# ═══════════════════════════════════════════

def _esc(s):
    return str(s).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def _chartjs_inline():
    p = SCRIPT_DIR.parent / "virome_discovery_pipeline" / "utils" / "chart.min.js"
    try:
        return "<script>" + p.read_text(encoding="utf-8") + "</script>"
    except Exception:
        return ""


def generate_html(rows, out_path, warn_retained):
    now = __import__("datetime").datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    n = len(rows)
    n_pass = sum(1 for r in rows if r["Status"] == "PASS")
    n_low = sum(1 for r in rows if r["Status"] == "LOW_RETAINED")
    n_bad = sum(1 for r in rows if r["Status"] in ("NO_OUTPUT", "NO_DATA"))
    rets = [r["Retained(%)"] for r in rows if r["Retained(%)"] is not None]
    median = sorted(rets)[len(rets) // 2] if rets else None

    chart_data = {"samples": [r["Sample"] for r in rows if r["Retained(%)"] is not None],
                  "retained": [r["Retained(%)"] for r in rows if r["Retained(%)"] is not None],
                  "warn": warn_retained}
    chartjs = _chartjs_inline()

    def fmt(v):
        return "" if v is None else (f"{v:,}" if isinstance(v, int) else v)

    body_rows = ""
    for r in rows:
        cls = {"PASS": "st-pass", "LOW_RETAINED": "st-low"}.get(r["Status"], "st-bad")
        body_rows += f'<tr class="{cls}">' + "".join(
            f"<td>{fmt(r[c])}</td>" for c in
            ["Sample", "Raw_Reads", "Fastp_Clean_Reads", "Fastp_Q30(%)", "Depletion_Input",
             "After_Kraken2", "After_HostAlign", "After_rRNA", "Retained(%)",
             "Aligner_Rate(%)", "Final_Files", "Status"]) + "</tr>"

    html = f'''<!DOCTYPE html><html lang="zh"><head><meta charset="UTF-8">
<title>Data Preprocessing Report</title>
<style>
:root{{--ink:#1a2332;--sub:#5a6a7e;--acc:#2563eb;--bg:#f6f8fb;--bd:#e1e5eb;--ok:#16a34a;--warn:#d97706;--bad:#dc2626}}
*{{box-sizing:border-box;margin:0;padding:0}}
body{{font-family:-apple-system,'Segoe UI',system-ui,sans-serif;color:var(--ink);background:var(--bg);padding:32px 40px}}
h1{{font-size:22px;margin-bottom:2px}} h2{{font-size:16px;margin:24px 0 8px}}
.meta{{color:var(--sub);font-size:12px;margin-bottom:16px}}
.cards{{display:flex;gap:10px;flex-wrap:wrap;margin:12px 0}}
.card{{background:#fff;border:1px solid var(--bd);border-radius:6px;padding:14px 20px;min-width:120px;flex:1;text-align:center}}
.card .v{{font-size:24px;font-weight:700;color:var(--acc)}} .card .l{{font-size:10px;color:var(--sub)}}
table{{border-collapse:collapse;width:100%;background:#fff;font-size:12px;border:1px solid var(--bd);border-radius:6px;overflow:hidden}}
th,td{{padding:6px 10px;text-align:left;border-bottom:1px solid #eef0f4;white-space:nowrap}}
th{{background:var(--bg);font-size:11px;color:var(--sub);position:sticky;top:0}}
tbody tr:hover td{{background:#eff4ff}}
.st-pass td:last-child{{color:var(--ok);font-weight:600}}
.st-low td:last-child{{color:var(--warn);font-weight:600}}
.st-bad td:last-child{{color:var(--bad);font-weight:600}}
.chart-card{{background:#fff;border:1px solid var(--bd);border-radius:6px;padding:14px;margin:12px 0;max-width:900px}}
.note{{background:#eff4ff;border:1px solid var(--bd);border-radius:6px;padding:10px 14px;font-size:12px;color:var(--sub);margin:10px 0}}
code{{background:var(--bg);padding:1px 6px;border-radius:4px;font-size:11px}}
</style></head><body>
<h1>Data Preprocessing Report — 数据清洗与去宿主</h1>
<div class="meta">Generated: {now} &nbsp;|&nbsp; 阈值: Retained &lt; {warn_retained}% → LOW_RETAINED</div>
<div class="cards">
<div class="card"><div class="v">{n}</div><div class="l">Samples</div></div>
<div class="card"><div class="v" style="color:var(--ok)">{n_pass}</div><div class="l">PASS</div></div>
<div class="card"><div class="v" style="color:var(--warn)">{n_low}</div><div class="l">LOW_RETAINED</div></div>
<div class="card"><div class="v" style="color:var(--bad)">{n_bad}</div><div class="l">NO_OUTPUT / NO_DATA</div></div>
<div class="card"><div class="v">{median if median is not None else "?"}</div><div class="l">Median Retained (%)</div></div>
</div>
<h2>交付物 / Handoff</h2>
<div class="note">下游交接清单: <code>assembly_ready.list</code> (最终 reads 绝对路径, 逐行一条) ·
<code>preprocessing_summary.tsv</code> (本表机器可读版)。<br>
下一步: <code>virome_discovery_pipeline/virome_pipeline.py --input_reads &lt;00b_HostDepletion 目录&gt;</code></div>
<div class="chart-card"><canvas id="ret" height="110"></canvas></div>
<h2>每样本明细 / Per-Sample Summary</h2>
<table><thead><tr>{"".join(f"<th>{c}</th>" for c in ["Sample","Raw_Reads","Fastp_Clean_Reads","Fastp_Q30(%)","Depletion_Input","After_Kraken2","After_HostAlign","After_rRNA","Retained(%)","Aligner_Rate(%)","Final_Files","Status"])}</tr></thead>
<tbody>{body_rows}</tbody></table>
<script>var CD={json.dumps(chart_data)};var hasChart={ "true" if chartjs else "false"};</script>
{chartjs}
<script>
if(hasChart&&typeof Chart!=='undefined'&&CD.samples.length){{
new Chart(document.getElementById('ret'),{{type:'bar',data:{{labels:CD.samples,datasets:[{{label:'Retained (%)',data:CD.retained,backgroundColor:'rgba(37,99,235,.6)'}}]}},
options:{{indexAxis:'y',plugins:{{legend:{{display:false}}}},scales:{{x:{{title:{{display:true,text:'Retained (%)'}}}}}}}}}});
}}else{{document.querySelector('.chart-card').style.display='none';}}
</script>
</body></html>'''
    out_path.write_text(html, encoding="utf-8")


# ═══════════════════════════════════════════
# Main
# ═══════════════════════════════════════════

def main():
    p = argparse.ArgumentParser(description="数据清洗管道 report 阶段: 汇总 + 交接清单")
    p.add_argument("--output-dir", "-o", help="输出根目录 (自动定位 00a_CleanData/00b_HostDepletion, 报告也写到这里)")
    p.add_argument("--clean-dir", help="clean-data.py 输出目录 (默认 <output-dir>/00a_CleanData)")
    p.add_argument("--deplete-dir", help="host_depletion.py 输出目录 (默认 <output-dir>/00b_HostDepletion)")
    p.add_argument("--warn-retained", type=float, default=20.0,
                   help="留存率低于该百分比标记 LOW_RETAINED (默认: 20)")
    args = p.parse_args()

    if not args.output_dir and not (args.clean_dir or args.deplete_dir):
        p.error("需要 --output-dir 或 --clean-dir/--deplete-dir")
    outdir = Path(args.output_dir).resolve() if args.output_dir else None
    clean_dir = Path(args.clean_dir).resolve() if args.clean_dir else (outdir / _D['d_clean'] if outdir else None)
    deplete_dir = Path(args.deplete_dir).resolve() if args.deplete_dir else (outdir / _D['d_hostdep'] if outdir else None)
    if not outdir:
        outdir = (deplete_dir or clean_dir).parent

    print(f"[report] clean:   {clean_dir or '(未提供)'}")
    print(f"[report] deplete: {deplete_dir or '(未提供)'}")

    rows = build_rows(clean_dir, deplete_dir, args.warn_retained)
    if not rows:
        print("[report] 未找到任何样本统计 (fastp JSON / seqkit 汇总 / 输出文件), 不生成报告")
        sys.exit(1)

    tsv = outdir / "preprocessing_summary.tsv"
    write_summary_tsv(rows, tsv)

    files = collect_final_files(deplete_dir, clean_dir)
    ready = outdir / "assembly_ready.list"
    write_ready_list(files, ready)

    html = outdir / "preprocessing_report.html"
    generate_html(rows, html, args.warn_retained)

    n_ready = sum(len(v) for v in files.values())
    print(f"[report] 样本: {len(rows)}  |  最终 reads 文件: {n_ready}")
    print(f"[report] 汇总表:   {tsv}")
    print(f"[report] 交接清单: {ready}")
    print(f"[report] HTML:     {html}")


if __name__ == "__main__":
    main()
