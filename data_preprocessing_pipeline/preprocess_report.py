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
import base64
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

_D = _layout_dirs(os.environ.get("MMPV_IO_LAYOUT"))

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


def _median(vals):
    """中位数: 偶数个取中间两数均值 (旧代码取上中位, 两样本时偏大)"""
    s = sorted(vals)
    if not s:
        return None
    mid = len(s) // 2
    return s[mid] if len(s) % 2 else (s[mid - 1] + s[mid]) / 2


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
    """host_depletion.py 的 {logs}/host_depletion_seqkit_summary.tsv (Sample/Stage/num_seqs...)

    PE 样本每 stage 落两行 (R1/R2 各一, seqkit stats 逐文件一行): 必须按端
    累加为合计, 旧代码后行覆盖前行 → 只剩 R2 的量。此处口径与 fastp JSON 的
    total_reads (本就两端合计) 一致, 故两者可直接对比/回退, 不再叠加。

    重跑语义 (按 (Sample, Stage) 分组, 只取"最后一批"):
      - 稳定路径 (1_Raw 输入 / 无 rrna 时的终产物): 重跑写同名文件 → 同名行
        视为同一批的覆盖写, 以后写为准, 不重复累加;
      - 临时随机路径 (2_Kraken2_Filtered / 3_Host_Filtered 的中间产物):
        重跑整批落在新临时目录 → 旧批与新批文件名完全不同, 按 file 列去重
        会双倍累加 (After_Kraken2/After_HostAlign 虚高), 故以「文件所在目录
        变化」标记新一批, 只累加最后一批。
    """
    out = {}
    if not deplete_dir:
        return out
    # (sample, stage键) → {批目录: {file: num_seqs}}; 批目录 = 该批文件所在目录
    # (同一批由一次 SeqkitStats.run 写出, 文件同目录; 重跑的临时批落在新目录)。
    # 跨多个汇总表 (历史运行残留) 共用一个批次流, 按路径排序保证 "最后一批" 稳定
    groups = {}
    last_dir = {}
    for p in sorted(Path(deplete_dir).rglob("host_depletion_seqkit_summary.tsv")):
        try:
            with open(p, encoding="utf-8") as f:
                rows = list(csv.DictReader(f, delimiter="\t"))
        except Exception:
            continue
        for i, r in enumerate(rows):
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
            if not (key and sample):
                continue
            n = _num(r.get("num_seqs"))
            if n is None:
                continue
            fname = (r.get("file") or "").strip() or f"__row{i}"
            fdir = os.path.dirname(fname)
            per_dir = groups.setdefault((sample, key), {})
            per_dir.setdefault(fdir, {})[fname] = n   # 同名 (稳定路径重跑) 后写覆盖
            last_dir[(sample, key)] = fdir            # 行序上最后出现的批即最后一批
    for (sample, key), per_dir in groups.items():
        out.setdefault(sample, {})[key] = sum(per_dir[last_dir[(sample, key)]].values())
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
    """bowtie2/hisat2 stderr 汇总行: '95.50% overall alignment rate' (百分比在前)

    真实输出顺序是「百分比 → overall alignment rate」, 旧正则按反向匹配, 恒不命中。
    两种排版都收: '<pct>% overall alignment rate' 与 'overall alignment rate: <pct>%'。
    日志在但无该行 (如 samtools / minimap2 的日志) 显式记 UNPARSED, 不留空。
    """
    out = {}
    if not deplete_dir:
        return out
    pat = re.compile(r"([\d.]+)%\s*overall alignment rate|overall alignment rate[:\s]+([\d.]+)",
                     re.IGNORECASE)
    for p in Path(deplete_dir).rglob("*.step2_align_filter.*.log"):
        # 日志命名 {sample}.step2_align_filter.{i}_{tool}.log; 样本名自带点
        # (如 ERR1.2) 时 split(".")[0] 会截断, 按 stage 标记 rsplit 取全名
        sample = p.name.rsplit(".step2_align_filter", 1)[0]
        try:
            text = p.read_text(encoding="utf-8", errors="replace")
        except Exception:
            continue
        m = pat.search(text)
        if not m:
            out.setdefault(sample, "UNPARSED")
            continue
        try:
            out[sample] = float(m.group(1) or m.group(2))
        except (TypeError, ValueError):
            out.setdefault(sample, "UNPARSED")
    return out


# 最终 reads 各级产物命名 → (样本, 端) 规则, 长后缀优先。
# 注意: host_depletion.py 的 auto_pair_files 保留输入的读端标签, 输入带 R 标签
# (X_R1.fa.gz) 时产物是 X_clean_R1/X_clean_R2; 不显式列出就会落到最后的通用兜底,
# 把 R1/R2 拆成两个幻影样本。
FINAL_PATTERNS = [
    (re.compile(r"^(?P<s>.+)_clean_1\.(?:fastq|fq|fa)\.gz$"), "pe1"),
    (re.compile(r"^(?P<s>.+)_clean_2\.(?:fastq|fq|fa)\.gz$"), "pe2"),
    (re.compile(r"^(?P<s>.+)_clean_R1\.(?:fastq|fq|fa)\.gz$"), "pe1"),
    (re.compile(r"^(?P<s>.+)_clean_R2\.(?:fastq|fq|fa)\.gz$"), "pe2"),
    (re.compile(r"^(?P<s>.+)_clean\.(?:fastq|fq|fa)\.gz$"), "se"),
    (re.compile(r"^(?P<s>.+)_1\.(?:fastq|fq|fa)\.gz$"), "pe1"),
    (re.compile(r"^(?P<s>.+)_2\.(?:fastq|fq|fa)\.gz$"), "pe2"),
    (re.compile(r"^(?P<s>.+)\.(?:fastq|fq|fa)\.gz$"), "se"),
]

# run_bbnorm.py 产物: {输入名去读端}_norm_{R1,R2,SE}.fq.gz; 去读端后样本名可能残留
# 上游 _clean / 扩展名 (SE 走 rsplit('.') 命名), 由可选组吃掉。
BBNORM_PATTERNS = [
    (re.compile(r"^(?P<s>.+?)(?:_clean)?(?:[._-][Rr]?[12])?_norm_R1\.(?:fastq|fq)\.gz$"), "pe1"),
    (re.compile(r"^(?P<s>.+?)(?:_clean)?(?:[._-][Rr]?[12])?_norm_R2\.(?:fastq|fq)\.gz$"), "pe2"),
    (re.compile(r"^(?P<s>.+?)(?:_clean)?(?:\.(?:fastq|fq|fasta|fa))?(?:[._-][Rr]?[12])?_norm_SE\.(?:fastq|fq)\.gz$"), "se"),
]


def _collect_source(root, patterns, source, found):
    """扫描一个来源目录, 归入 {样本: {来源: {端: Path}}}"""
    for f in sorted(root.rglob("*.gz")):
        if not f.is_file():
            continue
        if "rrna" in f.relative_to(root).parts[:-1]:  # --keep_rrna 的 rrna/ 子目录产物不是最终交付 reads
            continue
        for pat, role in patterns:
            m = pat.match(f.name)
            if m:
                found.setdefault(m.group("s"), {}).setdefault(source, {})[role] = f
                break


def collect_final_files(deplete_dir, clean_dir, bbnorm_dir=None):
    """最终 reads 文件: bbnorm 归一化 > 去宿主 > clean, 每样本只取一套完整 PE/SE

    返回 {样本: {"source": "normalized"|"depleted"|"clean", "files": [Path, ...]}}。
    旧实现把去宿主根与 clean 子目录的产物全量累积到同一列表, 同一样本既列
    deplete 又列 clean (下游会把同一份数据交两遍); 且 bbnorm 产物完全没进清单。
    """
    sources = []                                  # (来源标签, 根目录, 规则) 优先级从高到低
    if bbnorm_dir and Path(bbnorm_dir).is_dir():
        sources.append(("normalized", Path(bbnorm_dir), BBNORM_PATTERNS))
    if deplete_dir:
        sources.append(("depleted", Path(deplete_dir), FINAL_PATTERNS))
    if clean_dir:
        for sub in ("3.clumpify", "2.fasta"):
            d = Path(clean_dir) / sub
            if d.is_dir():
                sources.append(("clean", d, FINAL_PATTERNS))

    found = {}
    for label, root, pats in sources:
        _collect_source(root, pats, label, found)

    files = {}
    for s, per_src in found.items():
        for label, _root, _pats in sources:        # 取第一套完整 PE/SE
            roles = per_src.get(label)
            if not roles:
                continue
            if "pe1" in roles and "pe2" in roles:
                files[s] = {"source": label, "files": [roles["pe1"], roles["pe2"]]}
                break
            if "se" in roles:
                files[s] = {"source": label, "files": [roles["se"]]}
                break
    return files


# ═══════════════════════════════════════════
# 汇总
# ═══════════════════════════════════════════

def build_rows(clean_dir, deplete_dir, warn_retained, bbnorm_dir=None):
    fastp = collect_fastp(clean_dir)
    seqkit = collect_seqkit(deplete_dir)
    rrna = collect_rrna(deplete_dir)
    aligner = collect_aligner(deplete_dir)
    files = collect_final_files(deplete_dir, clean_dir, bbnorm_dir)

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
        finfo = files.get(s) or {}
        final = finfo.get("files", [])
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
            "Final_Source": finfo.get("source"),
            "Status": status,
        })
    return rows


SUMMARY_COLS = ["Sample", "Raw_Reads", "Fastp_Clean_Reads", "Fastp_Q30(%)", "Depletion_Input",
                "After_Kraken2", "After_HostAlign", "After_rRNA", "Retained(%)",
                "Aligner_Rate(%)", "Final_Files", "Final_Source", "Status"]


def write_summary_tsv(rows, out_path):
    cols = SUMMARY_COLS
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
    """交接清单: 每行一个最终 reads 文件绝对路径 (按样本分组, 每样本一套)"""
    with open(out_path, "w", encoding="utf-8", newline="") as f:
        for s in sorted(files):
            for pth in files[s]["files"]:
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


def _png_data_uri(path):
    """PNG → data URI (内嵌 HTML, 离线可看); 读不到返回 None"""
    try:
        return "data:image/png;base64," + base64.b64encode(Path(path).read_bytes()).decode("ascii")
    except Exception:
        return None


def find_deplete_plots(deplete_dir):
    """定位 host_depletion.py 的 line/bar PNG (logs_dir 默认 <deplete>/logs/)"""
    if not deplete_dir:
        return []
    root = Path(deplete_dir)
    out = []
    for d in (root / "logs", root):
        for kind in ("line", "bar"):
            p = d / f"host_depletion_plot_{kind}.png"
            if p.is_file() and p not in out:
                out.append(p)
    return out


def generate_html(rows, out_path, warn_retained, plot_paths=None):
    now = __import__("datetime").datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    n = len(rows)
    n_pass = sum(1 for r in rows if r["Status"] == "PASS")
    n_low = sum(1 for r in rows if r["Status"] == "LOW_RETAINED")
    n_bad = sum(1 for r in rows if r["Status"] in ("NO_OUTPUT", "NO_DATA"))
    rets = [r["Retained(%)"] for r in rows if r["Retained(%)"] is not None]
    median = _median(rets)

    chart_data = {"samples": [r["Sample"] for r in rows if r["Retained(%)"] is not None],
                  "retained": [r["Retained(%)"] for r in rows if r["Retained(%)"] is not None],
                  "warn": warn_retained}
    chartjs = _chartjs_inline()
    # 样本名进 <script> 前转义 '<' (json.dumps 不转义, </script> 会截断脚本块)
    chart_json = json.dumps(chart_data).replace("<", "\\u003c")

    def fmt(v):
        return "" if v is None else (f"{v:,}" if isinstance(v, int) else v)

    body_rows = ""
    for r in rows:
        cls = {"PASS": "st-pass", "LOW_RETAINED": "st-low"}.get(r["Status"], "st-bad")
        body_rows += f'<tr class="{cls}">' + "".join(
            f"<td>{_esc(fmt(r[c]))}</td>" for c in SUMMARY_COLS) + "</tr>"

    plot_imgs = ""
    for p in (plot_paths or []):
        uri = _png_data_uri(p)
        if uri:
            plot_imgs += (f'<div class="plot-card"><div class="plot-t">{_esc(Path(p).stem)}</div>'
                          f'<img src="{uri}" alt="{_esc(Path(p).name)}"></div>')
    plots_html = (f'<h2>Host Depletion Plots</h2><div class="plots">{plot_imgs}</div>'
                  if plot_imgs else "")

    html = f'''<!DOCTYPE html><html lang="en"><head><meta charset="UTF-8">
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
.plots{{display:flex;gap:12px;flex-wrap:wrap;margin:12px 0}}
.plot-card{{background:#fff;border:1px solid var(--bd);border-radius:6px;padding:14px;flex:1;min-width:360px}}
.plot-t{{font-size:12px;color:var(--sub);margin-bottom:8px}}
.plot-card img{{width:100%;height:auto}}
.note{{background:#eff4ff;border:1px solid var(--bd);border-radius:6px;padding:10px 14px;font-size:12px;color:var(--sub);margin:10px 0}}
code{{background:var(--bg);padding:1px 6px;border-radius:4px;font-size:11px}}
</style></head><body>
<h1>Data Preprocessing Report — Cleaning &amp; Host Depletion</h1>
<div class="meta">Generated: {now} &nbsp;|&nbsp; Threshold: Retained &lt; {warn_retained}% → LOW_RETAINED</div>
<div class="cards">
<div class="card"><div class="v">{n}</div><div class="l">Samples</div></div>
<div class="card"><div class="v" style="color:var(--ok)">{n_pass}</div><div class="l">PASS</div></div>
<div class="card"><div class="v" style="color:var(--warn)">{n_low}</div><div class="l">LOW_RETAINED</div></div>
<div class="card"><div class="v" style="color:var(--bad)">{n_bad}</div><div class="l">NO_OUTPUT / NO_DATA</div></div>
<div class="card"><div class="v">{median if median is not None else "?"}</div><div class="l">Median Retained (%)</div></div>
</div>
<h2>Handoff</h2>
<div class="note">Downstream handoff list: <code>assembly_ready.list</code> (absolute paths of final reads, one per line; one set per sample;
see the Final_Source column below: normalized=bbnorm &gt; depleted &gt; clean) ·
<code>preprocessing_summary.tsv</code> (machine-readable version of this table).<br>
Next step: <code>virome_discovery_pipeline/virome_pipeline.py --input_reads &lt;{_D['d_hostdep']} dir&gt;</code></div>
<div class="chart-card"><canvas id="ret" height="110"></canvas></div>
{plots_html}
<h2>Per-Sample Summary</h2>
<table><thead><tr>{"".join(f"<th>{c}</th>" for c in SUMMARY_COLS)}</tr></thead>
<tbody>{body_rows}</tbody></table>
<script>var CD={chart_json};var hasChart={ "true" if chartjs else "false"};</script>
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

    # bbnorm 产物 (--stage bbnorm, 仅在显式跑过时存在): 交接优先级 bbnorm > 去宿主 > clean
    bbnorm_dir = outdir / _D['d_bbnorm']
    bbnorm_dir = bbnorm_dir if bbnorm_dir.is_dir() else None

    rows = build_rows(clean_dir, deplete_dir, args.warn_retained, bbnorm_dir)
    if not rows:
        print("[report] 未找到任何样本统计 (fastp JSON / seqkit 汇总 / 输出文件), 不生成报告")
        sys.exit(1)

    tsv = outdir / "preprocessing_summary.tsv"
    write_summary_tsv(rows, tsv)

    files = collect_final_files(deplete_dir, clean_dir, bbnorm_dir)
    ready = outdir / "assembly_ready.list"
    write_ready_list(files, ready)

    html = outdir / "preprocessing_report.html"
    plots = find_deplete_plots(deplete_dir)
    generate_html(rows, html, args.warn_retained, plots)
    if not plots:
        print("[report] 未找到 host_depletion_plot_{line,bar}.png, HTML 不含可视化小节")

    # AI 解读 section (可选, 纯增量: 离线/未配 key 时自动跳过, 不改动报告)
    try:
        from mmpv_common.ai import append_ai_section_to_html
        from mmpv_common.run_digest import guess_state_dir
        n_pass = sum(1 for r in rows if r["Status"] == "PASS")
        n_low = sum(1 for r in rows if r["Status"] == "LOW_RETAINED")
        n_bad = sum(1 for r in rows if r["Status"] in ("NO_OUTPUT", "NO_DATA"))
        rets = [r["Retained(%)"] for r in rows
                if r["Retained(%)"] is not None]
        append_ai_section_to_html(html, {
            "pipeline": "preprocessing", "n_samples": len(rows),
            "pass": n_pass, "low_retained": n_low, "no_output_or_no_data": n_bad,
            "median_retained_pct": _median(rets),
        }, state_dir=guess_state_dir(outdir))
    except Exception:
        pass
    # 过程面板 (确定性, 无需 AI): 项目 90_Handoff/state 在时嵌入全过程总览
    try:
        from mmpv_common.process_report import append_process_panel_to_html
        append_process_panel_to_html(html, outdir)
    except Exception:
        pass

    n_ready = sum(len(v["files"]) for v in files.values())
    print(f"[report] 样本: {len(rows)}  |  最终 reads 文件: {n_ready}")
    print(f"[report] 汇总表:   {tsv}")
    print(f"[report] 交接清单: {ready}")
    print(f"[report] HTML:     {html}")


if __name__ == "__main__":
    main()
