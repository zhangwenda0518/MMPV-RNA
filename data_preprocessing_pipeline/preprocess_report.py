#!/usr/bin/env python3
"""
preprocess_report.py — 数据清洗管道 report 阶段 (汇总 + 交接清单)
==================================================================
汇聚 clean-data.py (fastp JSON) 与 host_depletion.py (seqkit 汇总表 /
ribodetector 报告 / 比对日志) 的既有统计产物, 生成:

  preprocessing_summary.tsv   每样本各步 reads 数 / 留存率 / 状态 (末行 TOTAL)
  assembly_ready.list         最终 reads 文件清单 (下游 --input_reads 交接)
  preprocessing_report.html   交互式 HTML 报告 (分节 + 侧栏导航, Chart.js 内嵌离线可用)

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

try:
    from mmpv_common.table_widgets import TW_BLOCK as _TW
except Exception:
    _TW = ""   # 共享表格控件不可达时静默降级为原生表格

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
            "Fastp_Dup(%)": fp.get("dup"),   # fastp duplication.rate, 缺失留空
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
                "Aligner_Rate(%)", "Final_Files", "Final_Source", "Status", "Fastp_Dup(%)"]

# 比率列: TOTAL 行不求和 (留存率按输入加权, 其余留空)
_RATE_COLS = ("Fastp_Q30(%)", "Fastp_Dup(%)", "Retained(%)", "Aligner_Rate(%)")


def _total_row(rows):
    """TOTAL 行 (TSV 与 HTML 表尾共用同一份计算, 两处保持一致)。

    计数列求和; 比率列不求和 — Retained(%) 按 Depletion_Input 加权,
    其余比率列与字符串列 (Final_Source / Status) 留空。
    """
    total = {"Sample": "TOTAL"}
    for c in SUMMARY_COLS[1:]:
        if c in _RATE_COLS or c in ("Final_Source", "Status"):
            total[c] = ""
        else:
            total[c] = sum(r.get(c) for r in rows if isinstance(r.get(c), (int, float)))
    wsum = wbase = 0
    for r in rows:
        if r.get("Retained(%)") is not None and r.get("Depletion_Input"):
            wsum += r["Retained(%)"] / 100 * r["Depletion_Input"]
            wbase += r["Depletion_Input"]
    total["Retained(%)"] = round(wsum / wbase * 100, 1) if wbase else ""
    return total


def write_summary_tsv(rows, out_path):
    cols = SUMMARY_COLS
    with open(out_path, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)
        w.writerow(_total_row(rows))     # TOTAL 行 (数值列求和, 留存率按输入加权)


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


# ═══════════════════════════════════════════
# 老版式 (③ virome_discovery_pipeline/report_pipeline.py): CSS/JS 原样拷贝
# ═══════════════════════════════════════════
# ② 本就是从老一体化 pipeline (现 ③ report_pipeline.py 的前身) 拆出来的, 报告
# 沿用同一套 "stage 卡" 视觉语言: 固定侧栏 (Pipeline Modules) + 00a/00b/[00c]
# 三张 stage 卡 + stage-tables/stage-charts。下列规则自 ③ 的 <style>/<script>
# 逐条拷贝 (仅保留本报告用到的部分; 单一花括号, 保持非 f-string 常量以便原样内联)。

_STYLE_RULES = """:root{
  --bg:#F8F4ED;--card-bg:#FCFAF5;--text:#3B3D3F;--muted:#6B6F73;
  --navy:#1B365D;--indigo:#1B365D;--blue:#537D96;--green:#4A6B4A;
  --red:#8B3A3A;--amber:#9D5F4D;--border:rgba(122,96,88,.18);
  --shadow:0 1px 3px rgba(60,48,40,.05);--shadow-lg:0 3px 10px rgba(60,48,40,.08);
  --radius:4px;--radius-sm:2px;
}
*{box-sizing:border-box;margin:0;padding:0}
body{font-family:'EB Garamond','Noto Serif SC','Source Han Serif SC','Songti SC','STSong',serif;background:var(--bg);color:var(--text);line-height:1.7}
.container{max-width:1200px;margin:0 auto;padding:0 20px 40px;padding-left:240px}
.hero{background:#1B365D;color:#FFFDF7;padding:40px 32px 32px;border-radius:0 0 8px 8px;margin-bottom:28px;position:relative;overflow:hidden}
.hero h1{font-size:26px;font-weight:500;margin-bottom:6px;position:relative;z-index:1}
.hero .subtitle{font-size:13px;opacity:.85;position:relative;z-index:1}
.hero .gen-time{font-size:12px;opacity:.7;margin-top:6px;position:relative;z-index:1}
.kpi-row{display:grid;grid-template-columns:repeat(auto-fit,minmax(180px,1fr));gap:14px;margin-bottom:20px}
.kpi-card{background:var(--card-bg);border-radius:var(--radius);padding:18px 20px;box-shadow:var(--shadow);text-align:center;border:1px solid var(--border)}
.kpi-icon{font-size:24px;margin-bottom:6px}
.kpi-value{font-size:13px;font-weight:600;color:var(--text);line-height:1.4}
.kpi-label{font-size:11px;color:var(--muted);margin-top:4px;text-transform:uppercase;letter-spacing:.5px}
.sidebar{position:fixed;top:0;left:0;width:224px;height:100vh;background:var(--card-bg);border-right:1px solid var(--border);box-shadow:2px 0 8px rgba(0,0,0,.04);z-index:200;display:flex;flex-direction:column;padding-top:12px;overflow-y:auto}
.sb-title{font-size:11px;font-weight:700;color:var(--muted);text-transform:uppercase;letter-spacing:.8px;padding:8px 16px 6px;margin:0}
.sb-item{display:flex;align-items:center;gap:10px;padding:9px 16px 9px 19px;text-decoration:none;color:var(--text);font-size:12.5px;font-weight:500;transition:all .15s}
.sb-item:hover{background:rgba(83,125,150,.08);color:#3F6179}
.sb-dot{width:10px;height:10px;border-radius:50%;flex-shrink:0}
.sb-label{overflow:hidden;text-overflow:ellipsis;white-space:nowrap}
.sb-pass{background:rgba(74,107,74,.05)}.sb-fail{background:rgba(139,58,58,.05)}.sb-skip{background:rgba(107,111,115,.035)}
.sb-active{background:rgba(83,125,150,.1);font-weight:500;color:#3F6179}
.stage{background:var(--card-bg);border-radius:var(--radius);box-shadow:var(--shadow);margin-bottom:20px;overflow:hidden}
.stage-pass{border:1px solid rgba(74,107,74,.38);background:rgba(74,107,74,.045)}.stage-fail{border:1px solid rgba(139,58,58,.38);background:rgba(139,58,58,.045)}.stage-skip{border:1px solid var(--border);background:rgba(107,111,115,.04)}
.stage-header{display:flex;align-items:center;gap:16px;padding:18px 22px;border-bottom:1px solid var(--border);background:var(--bg)}
.stage-icon{width:44px;height:44px;border-radius:var(--radius-sm);display:flex;align-items:center;justify-content:center;font-size:11px;font-weight:500;color:#FFFDF7;flex-shrink:0}
.stage-title{flex:1;min-width:0}
.stage-title h2{font-size:16px;font-weight:500;color:var(--text);margin-bottom:2px}
.stage-metric{font-size:12px;color:var(--muted);display:block;overflow:hidden;text-overflow:ellipsis;white-space:nowrap}
.stage-badge{font-size:11px;font-weight:500;padding:4px 12px;border-radius:var(--radius-sm);flex-shrink:0}
.s-pass{background:rgba(74,107,74,.12);color:var(--green)}.s-fail{background:rgba(139,58,58,.12);color:var(--red)}.s-skip{background:rgba(0,0,0,.05);color:var(--muted)}
.stage-charts{display:grid;gap:16px;padding:18px 22px;grid-template-columns:repeat(auto-fit,minmax(0,1fr))}
.chart-box{background:var(--bg);border-radius:var(--radius-sm);padding:12px;border:1px solid var(--border);position:relative;min-width:0}
.pct-btn{font-size:11px;padding:4px 14px;border:1px solid var(--blue);border-radius:var(--radius-sm);background:var(--blue);cursor:pointer;color:#FFFDF7;font-weight:500;transition:all .15s}
.pct-btn:hover{background:var(--accent-hover,#3F6179);border-color:#3F6179}
.stage-tables{padding:0 22px 14px;display:flex;flex-direction:column;gap:8px}
.stage-table-detail{border-top:1px solid var(--border);padding:8px 0 4px}
.stage-table-detail summary{cursor:pointer;font-size:12px;padding:2px 0;color:var(--muted)}
.stage-table-detail summary:hover{color:var(--blue)}
.table-wrap{background:var(--card-bg);border-radius:var(--radius);box-shadow:var(--shadow);overflow:hidden;margin-bottom:20px}
table{width:100%;border-collapse:collapse;font-size:13px}
th{background:var(--bg);text-align:left;padding:10px 16px;font-size:11px;color:var(--muted);text-transform:uppercase;letter-spacing:.5px;border-bottom:2px solid var(--border)}
td{padding:9px 16px;border-bottom:1px solid var(--border)}
.tr-pass{background:rgba(74,107,74,.07)}.tr-fail{background:rgba(139,58,58,.07)}.tr-skip{background:rgba(107,111,115,.05)}
.tb-badge{display:inline-block;padding:1px 8px;border-radius:var(--radius-sm);font-size:11px;font-weight:500}
.tb-pass{background:rgba(74,107,74,.12);color:var(--green)}.tb-fail{background:rgba(139,58,58,.12);color:var(--red)}.tb-skip{background:rgba(0,0,0,.05);color:var(--muted)}
.app-table{width:100%;border-collapse:collapse;font-size:11px;margin-bottom:8px}
.app-table th{background:#f5f5f5;padding:6px 8px;font-size:11px;text-align:left;border:1px solid #e0e0e0;white-space:nowrap}
.app-table td{padding:4px 8px;border:1px solid #f0f0f0;max-width:200px;overflow:hidden;text-overflow:ellipsis;white-space:nowrap;font-size:11px}
.footer{text-align:center;padding:24px;color:var(--muted);font-size:12px;line-height:1.8}
.footer a{color:var(--blue);text-decoration:none}
@media(prefers-reduced-motion:reduce){.kpi-card,.sb-item,.pct-btn,#chat-btn,#ai-btn,.stage-table-detail summary{transition:none!important}
.chart-box img{transition:none!important}}
/* 断点 560px（非 ③ 的 768px）：预览窗格/分屏常见 600-800px 宽，768 会把侧栏藏掉 */
@media(max-width:560px){
  .sidebar{display:none}
  .container{padding-left:20px}
  .hero{padding:24px 20px 20px}.hero h1{font-size:20px}
  .kpi-row{grid-template-columns:repeat(2,1fr)}
  .stage-header{flex-wrap:wrap;gap:10px}
  .stage-charts{grid-template-columns:1fr!important}
  .footer{font-size:12px}}
@media print{
  body{background:#fff;font-size:11px}
  .hero{background:#1B365D!important;-webkit-print-color-adjust:exact}
  .stage,.table-wrap,.kpi-card{box-shadow:none;border:1px solid var(--border);break-inside:avoid}
  .sidebar{display:none}
  .container{padding-left:20px}
}
"""

# ③ <script> 内实现逐行拷贝 (lightbox / scroll-spy / 堆叠图 % 切换 / 表格分页)
_LIGHTBOX_JS = """function openLightbox(src){document.getElementById('lightbox-img').src=src;document.getElementById('lightbox').style.display='flex';}
function closeLightbox(){document.getElementById('lightbox').style.display='none';}
"""

_SPY_JS = """// Sidebar scroll-spy
(function(){
  const items=document.querySelectorAll('.sb-item');
  const stages=document.querySelectorAll('.stage[id]');
  if(!items.length||!stages.length)return;
  function onScroll(){
    let current='';
    stages.forEach(s=>{if(s.getBoundingClientRect().top<=160)current=s.id});
    items.forEach(a=>{
      a.classList.toggle('sb-active',a.getAttribute('href')==='#'+current);
    });
  }
  window.addEventListener('scroll',onScroll,{passive:true});
  onScroll();
})();
"""

_PCT_JS = """// Stacked bar chart % toggle
(function(){
  var _orig={};
  window.toggleStackedPct=function(btn){
    var cid=btn.getAttribute('data-chart');
    var chart=Chart.getChart(cid);
    if(!chart)return;
    var ds=chart.data.datasets;
    var n=chart.data.labels.length;
    // Value axis: horizontal charts (indexAxis='y') have the value axis on x; vertical charts on y
    var valAxis=chart.options.indexAxis==='y'?'x':'y';
    if(btn.getAttribute('data-mode')==='abs'){
      if(!_orig[cid])_orig[cid]=ds.map(function(d){return d.data.slice()});
      for(var i=0;i<n;i++){
        var tot=0;
        for(var j=0;j<ds.length;j++)tot+=Number(_orig[cid][j][i])||0;
        for(var j=0;j<ds.length;j++)ds[j].data[i]=tot>0?((Number(_orig[cid][j][i])||0)/tot*100):0;
      }
      chart.options.scales[valAxis].title={display:true,text:'%'};
      chart.options.scales[valAxis].max=100;chart.options.scales[valAxis].min=0;
      btn.setAttribute('data-mode','pct');btn.textContent='Show Count';
    }else{
      if(_orig[cid])for(var j=0;j<ds.length;j++)ds[j].data=_orig[cid][j].slice();
      chart.options.scales[valAxis].title={display:true,text:'Reads'};
      chart.options.scales[valAxis].max=undefined;chart.options.scales[valAxis].min=0;
      btn.setAttribute('data-mode','abs');btn.textContent='Show %';
    }
    chart.update();
  };
})();
"""

# 表格分页 (③ 的 pag-row 变体逐行拷贝; __TBID__ = 表 id, 行数 >10 时才输出)
_PAGER_JS_TMPL = """(function(){
  var tbody=document.getElementById('__TBID__').querySelector('tbody');
  var rows=Array.from(tbody.querySelectorAll('.pag-row'));
  var perPage=10,totalPages=Math.ceil(rows.length/perPage);
  var pager=document.getElementById('pager___TBID__');
  var page=0;
  function show(p){
    page=Math.max(0,Math.min(p,totalPages-1));
    rows.forEach(function(r,i){r.style.display=(i>=page*perPage&&i<(page+1)*perPage)?'':'none'});
    pager.innerHTML='<button onclick="void(0)" '+(page===0?'disabled':'')+' style="border:1px solid #ccc;background:#fff;padding:3px 10px;border-radius:4px;cursor:pointer;font-size:11px">Prev</button>'+'<span>Page <b>'+(page+1)+'</b> of '+totalPages+'</span>'+'<button '+(page===totalPages-1?'disabled':'')+' style="border:1px solid #ccc;background:#fff;padding:3px 10px;border-radius:4px;cursor:pointer;font-size:11px">Next</button>';
    pager.querySelectorAll('button')[0].onclick=function(){show(page-1)};
    pager.querySelectorAll('button')[1].onclick=function(){show(page+1)};
  }
  show(0);
})();
"""


def _rel_href(base_dir, target):
    """报告输出目录 → 目标文件的相对链接 (posix 风格; 报告与产物可跨层级)"""
    rel = os.path.relpath(str(Path(target).resolve()), str(Path(base_dir).resolve()))
    return Path(rel).as_posix()


def _esc_attr(s):
    """HTML 属性值转义 (_esc 不处理引号, href 等属性值需额外转义)"""
    return _esc(s).replace('"', "&quot;")


def find_fastp_reports(clean_dir):
    """clean-data.py 的每样本 fastp HTML 报告: {clean}/logs/{sample}_fastp_report.html"""
    if not clean_dir:
        return []
    return [p for p in sorted(Path(clean_dir).rglob("*_fastp_report.html")) if p.is_file()]


def find_kraken2_reports(deplete_dir):
    """host_depletion.py 的每样本 kraken2 报告: {logs}/kraken2_report/{sample}_kraken2_report.txt"""
    if not deplete_dir:
        return []
    root = Path(deplete_dir)
    dirs = [d for d in [root / "logs" / "kraken2_report"]
            + sorted(root.rglob("kraken2_report")) if d.is_dir()]
    out, seen = [], set()
    for d in dirs:
        for f in sorted(d.glob("*_kraken2_report.txt")):
            if f.is_file() and f not in seen:
                seen.add(f)
                out.append(f)
    return out


def find_resource_usage(deplete_dir):
    """host_depletion.py 的每样本资源用量表: {logs}/host_depletion_resource_usage.tsv"""
    if not deplete_dir:
        return None
    cands = sorted(Path(deplete_dir).rglob("host_depletion_resource_usage.tsv"))
    return cands[0] if cands else None


def generate_html(rows, out_path, warn_retained, plot_paths=None,
                  clean_dir=None, deplete_dir=None):
    """老版式 stage 卡 HTML 报告 (③ report_pipeline 同款视觉语言)。

    结构: 固定侧栏 (Pipeline Modules: 00a/00b/[00c]) + 三张 stage 卡
    (#stage-CleanData / #stage-HostDep / [#stage-BBNorm]) + 页尾 Handoff
    details。bbnorm 卡仅在检测到归一化产物 (Final_Source=normalized 或
    03_BBnorm 目录) 时渲染, 侧栏同步隐藏。可选 clean_dir/deplete_dir 只用于
    给「每样本明细」表 (fastp HTML 报告 / kraken2 报告 / 资源用量表) 生成
    相对链接; 省略或文件缺失时不渲染对应块。统计数值全部来自 rows,
    原 3 参调用形式保持兼容。
    """
    out_path = Path(out_path)
    base_dir = out_path.parent                    # 相对链接基准 = 报告输出目录
    now = __import__("datetime").datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    n = len(rows)
    n_pass = sum(1 for r in rows if r["Status"] == "PASS")
    n_low = sum(1 for r in rows if r["Status"] == "LOW_RETAINED")
    n_bad = sum(1 for r in rows if r["Status"] in ("NO_OUTPUT", "NO_DATA"))
    rets = [r["Retained(%)"] for r in rows if r["Retained(%)"] is not None]
    median = _median(rets)

    chartjs = _chartjs_inline()

    def fmt(v):
        return "" if v is None else (f"{v:,}" if isinstance(v, int) else v)

    def dash(v):
        return v if v is not None else "?"

    # ── Per-Sample Summary 行 (老表格语言: tr-pass/tr-fail/tr-skip + tb-badge 状态) ──
    _TR_CLS = {"PASS": "tr-pass", "OK": "tr-pass", "LOW_RETAINED": "tr-fail",
               "NO_OUTPUT": "tr-fail", "NO_DATA": "tr-skip"}
    _TB_CLS = {"PASS": "tb-pass", "OK": "tb-pass", "LOW_RETAINED": "tb-fail",
               "NO_OUTPUT": "tb-fail", "NO_DATA": "tb-skip"}
    _TBID = "stagetbl_s00b_preprocessing_summary"
    use_pagination = len(rows) > 10           # >10 行走 .pag-row + Prev/Next (③ 同规则)
    summary_body = ""
    for r in rows:
        cls = _TR_CLS.get(r["Status"], "tr-skip") + (" pag-row" if use_pagination else "")
        row_attr = ' style="display:none"' if use_pagination else ""
        cells = ""
        for c in SUMMARY_COLS:
            if c == "Status":
                cells += (f'<td><span class="tb-badge {_TB_CLS.get(r["Status"], "tb-skip")}">'
                          f'{_esc(r["Status"])}</span></td>')
            else:
                cells += f"<td>{_esc(fmt(r.get(c)))}</td>"
        summary_body += f'<tr class="{cls}"{row_attr}>{cells}</tr>'
    # 表尾 TOTAL 行 (与 TSV 的 TOTAL 共用 _total_row; 比率列不求和)
    tfoot = ""
    if rows:
        _tot = _total_row(rows)
        tfoot = ('<tfoot><tr style="font-weight:600;background:#f5f5f5">'
                 + "".join(f"<td>{_esc(fmt(_tot.get(c)))}</td>" for c in SUMMARY_COLS)
                 + "</tr></tfoot>")
    # 表格分页 (③ 同款: >10 行时输出空 pager div, JS 按 .pag-row 逐页显示)
    pager_html = pager_js = ""
    if use_pagination:
        pager_html = (f'<div id="pager_{_TBID}" style="display:flex;justify-content:center;'
                      f'gap:8px;padding:8px;font-size:11px;color:var(--muted)"></div>')
        pager_js = _PAGER_JS_TMPL.replace("__TBID__", _TBID)

    # ── stage-charts 数据: 留存率堆叠条 (老 chart_s00b: Retained vs Removed) ──
    chart_rows = []
    for r in rows:
        base_n = r.get("Depletion_Input")
        if base_n is None:
            base_n = r.get("Fastp_Clean_Reads")
        if base_n and r.get("Retained(%)") is not None:
            ret_n = int(round(r["Retained(%)"] / 100 * base_n))
            chart_rows.append((r["Sample"], ret_n, max(0, int(base_n) - ret_n)))
    chart_data = {"labels": [c[0] for c in chart_rows],
                  "retained": [c[1] for c in chart_rows],
                  "removed": [c[2] for c in chart_rows],
                  "warn": warn_retained}
    # 样本名进 <script> 前转义 '<' (json.dumps 不转义, </script> 会截断脚本块)
    chart_json = json.dumps(chart_data).replace("<", "\\u003c")
    chart_boxes = ""
    if chart_rows:
        chart_boxes += ('<div class="chart-box">'
                        '<div style="display:flex;justify-content:flex-end;margin-bottom:6px">'
                        '<button class="pct-btn" data-chart="chart_s00b" data-mode="abs" '
                        'onclick="toggleStackedPct(this)">Show %</button></div>'
                        '<canvas id="chart_s00b" style="max-height:320px"></canvas></div>')
    for p in (plot_paths or []):
        uri = _png_data_uri(p)
        if uri:
            chart_boxes += (f'<div class="chart-box"><div class="chart-title">{_esc(Path(p).stem)}</div>'
                            f'<img src="{uri}" loading="lazy" onclick="openLightbox(this.src)" '
                            f'title="Click to zoom" alt="{_esc(Path(p).name)}" '
                            f'style="width:100%;height:auto;max-height:650px;object-fit:contain;'
                            f'display:block;cursor:zoom-in"></div>')
    charts_html = (f'<div class="stage-charts" style="grid-template-columns:1fr">{chart_boxes}</div>'
                   if chart_boxes else "")

    # ── 00a Data Preprocessing 聚合 (fastp QC) ──
    q30s = [r["Fastp_Q30(%)"] for r in rows if r.get("Fastp_Q30(%)") is not None]
    dups = [r["Fastp_Dup(%)"] for r in rows if r.get("Fastp_Dup(%)") is not None]
    n_clean = sum(1 for r in rows if isinstance(r.get("Fastp_Clean_Reads"), (int, float)))
    mean_q30 = round(sum(q30s) / len(q30s), 1) if q30s else None
    mean_dup = round(sum(dups) / len(dups), 1) if dups else None
    mean_ret = round(sum(rets) / len(rets), 1) if rets else None
    clean_metric = " · ".join([f"{n_clean} Samples"]
                              + ([f"mean Q30 {mean_q30}%"] if mean_q30 is not None else [])
                              + ([f"mean dup {mean_dup}%"] if mean_dup is not None else [])
                              + ([f"mean retained {mean_ret}%"] if mean_ret is not None else []))
    clean_pass = n_clean > 0
    # 每样本 fastp 报告 → 00a stage-tables 内的 details (文件不存在则不渲染)
    fp_tables = ""
    fp_reports = find_fastp_reports(clean_dir)
    if fp_reports:
        by_sample = {r["Sample"]: r for r in rows}
        fp_rows = ""
        for p in fp_reports:
            sname = p.name[:-len("_fastp_report.html")] or p.name
            rr = by_sample.get(sname) or {}
            fp_rows += (f'<tr><td>{_esc(sname)}</td>'
                        f'<td>{_esc(fmt(rr.get("Raw_Reads")))}</td>'
                        f'<td>{_esc(fmt(rr.get("Fastp_Clean_Reads")))}</td>'
                        f'<td>{_esc(fmt(rr.get("Fastp_Q30(%)")))}</td>'
                        f'<td>{_esc(fmt(rr.get("Fastp_Dup(%)")))}</td>'
                        f'<td><a href="{_esc_attr(_rel_href(base_dir, p))}">{_esc(p.name)}</a></td></tr>')
        fp_tables = f'''<div class="stage-tables">
<details class="stage-table-detail" open>
<summary>Per-sample fastp QC reports ({len(fp_reports)}) <span style="color:var(--muted);font-size:11px">({len(fp_reports)} rows)</span></summary>
<div style="overflow-x:auto;margin-top:6px">
<table class="app-table" id="stagetbl_s00a_fastp_qc"><thead><tr><th>Sample</th><th>Raw_Reads</th><th>Clean_Reads</th><th>Q30(%)</th><th>Fastp_Dup(%)</th><th>Report</th></tr></thead>
<tbody>{fp_rows}</tbody></table>
</div>
</details>
</div>'''

    # ── 00a stage 卡 (老版式: stage-header + stage-icon + stage-metric + stage-badge) ──
    stage_cleandata = f'''<section class="stage {"stage-pass" if clean_pass else "stage-skip"}" id="stage-CleanData">
  <div class="stage-header">
    <div class="stage-icon" style="background:#607d8b">CLN</div>
    <div class="stage-title">
<h2>00a Data Preprocessing</h2>
<span class="stage-metric">{clean_metric}</span>
</div>
    <span class="stage-badge {"s-pass" if clean_pass else "s-skip"}">{"✓ PASS" if clean_pass else "○ SKIP"}</span>
  </div>
{fp_tables}
</section>'''

    # ── 00b Host Depletion 聚合 + 徽章 (由该阶段样本状态推导) ──
    dep_in = [r["Depletion_Input"] for r in rows if isinstance(r.get("Depletion_Input"), (int, float))]
    align_rates = [r["Aligner_Rate(%)"] for r in rows
                   if isinstance(r.get("Aligner_Rate(%)"), (int, float))]
    dep_metric = " · ".join([f"{len(dep_in)} Samples"]
                            + ([f"median aligner rate {_median(align_rates)}%"]
                               if align_rates else []))
    if n_pass > 0:                                     # 有 PASS 样本 → ✓ PASS
        dep_stage_cls, dep_badge, dep_badge_txt = "stage-pass", "s-pass", "✓ PASS"
    elif rows and n_bad == len(rows):                  # 统计齐但全部异常 → ✗ FAIL
        dep_stage_cls, dep_badge, dep_badge_txt = "stage-fail", "s-fail", "✗ FAIL"
    else:                                              # 无 STATS 样本 → ○ SKIP
        dep_stage_cls, dep_badge, dep_badge_txt = "stage-skip", "s-skip", "○ SKIP"
    dep_sb_cls = {"stage-pass": "sb-pass", "stage-fail": "sb-fail",
                  "stage-skip": "sb-skip"}[dep_stage_cls]
    # kraken2 / 资源用量 → 00b stage-tables 内的 details (缺失则不渲染)
    dep_tables = ""
    k2_reports = find_kraken2_reports(deplete_dir)
    if k2_reports:
        k2_rows = "".join(
            f'<tr><td>{_esc(p.name[:-len("_kraken2_report.txt")] or p.name)}</td>'
            f'<td><a href="{_esc_attr(_rel_href(base_dir, p))}">{_esc(p.name)}</a></td></tr>'
            for p in k2_reports)
        dep_tables += f'''<details class="stage-table-detail" open>
<summary>{len(k2_reports)} per-sample Kraken2 reports <span style="color:var(--muted);font-size:11px">({len(k2_reports)} rows)</span></summary>
<div style="overflow-x:auto;margin-top:6px">
<table class="app-table" id="stagetbl_s00b_kraken2"><thead><tr><th>Sample</th><th>Report</th></tr></thead>
<tbody>{k2_rows}</tbody></table>
</div>
</details>'''
    res_tsv = find_resource_usage(deplete_dir)
    if res_tsv:
        dep_tables += (f'<details class="stage-table-detail">'
                       f'<summary>Host Depletion Resource Usage '
                       f'<span style="color:var(--muted);font-size:11px">({_esc(res_tsv.name)})</span></summary>'
                       f'<div style="padding:6px 0 2px;font-size:12px">Resource usage table: '
                       f'<a href="{_esc_attr(_rel_href(base_dir, res_tsv))}" '
                       f'style="color:var(--blue)">{_esc(res_tsv.name)}</a></div></details>')

    # ── 00b stage 卡: Per-Sample Summary 表 (app-table + tfoot TOTAL) + 明细表 + 图 ──
    stage_hostdep = f'''<section class="stage {dep_stage_cls}" id="stage-HostDep">
  <div class="stage-header">
    <div class="stage-icon" style="background:#546e7a">DEP</div>
    <div class="stage-title">
<h2>00b Host Depletion</h2>
<span class="stage-metric">{dep_metric}</span>
</div>
    <span class="stage-badge {dep_badge}">{dep_badge_txt}</span>
  </div>
  <div class="stage-tables">
<details class="stage-table-detail" open>
<summary>preprocessing_summary.tsv <span style="color:var(--muted);font-size:11px">({len(rows)} rows)</span></summary>
<div style="overflow-x:auto;margin-top:6px">
<table class="app-table" id="{_TBID}"><thead><tr>{"".join(f"<th>{c}</th>" for c in SUMMARY_COLS)}</tr></thead>
<tbody>{summary_body}</tbody>{tfoot}</table>
{pager_html}
</div>
</details>
{dep_tables}
  </div>
{charts_html}
</section>'''

    # ── 00c BBNorm: 仅检测到归一化产物/目录时渲染 (侧栏同步) ──
    n_norm = sum(1 for r in rows if r.get("Final_Source") == "normalized")
    n_norm_files = sum(r.get("Final_Files") or 0
                       for r in rows if r.get("Final_Source") == "normalized")
    bbnorm_on = n_norm > 0 or (base_dir / _D["d_bbnorm"]).is_dir()
    bb_metric = f"{n_norm} samples normalized · {n_norm_files} files"
    bb_pass = n_norm > 0
    stage_bbnorm = ""
    if bbnorm_on:
        stage_bbnorm = f'''<section class="stage {"stage-pass" if bb_pass else "stage-skip"}" id="stage-BBNorm">
  <div class="stage-header">
    <div class="stage-icon" style="background:#00838f">BBN</div>
    <div class="stage-title">
<h2>00c BBNorm</h2>
<span class="stage-metric">{bb_metric}</span>
</div>
    <span class="stage-badge {"s-pass" if bb_pass else "s-skip"}">{"✓ PASS" if bb_pass else "○ SKIP"}</span>
  </div>
</section>'''

    # ── 侧栏 (老报告同款: Pipeline Modules + 色点 CLN/DEP/BBN; bbnorm 关闭时不出现) ──
    sb_items = (f'<a href="#stage-CleanData" class="sb-item {"sb-pass" if clean_pass else "sb-skip"}" '
                f'title="{_esc_attr(clean_metric)}">\n'
                '<span class="sb-dot" style="background:#607d8b"></span>\n'
                '<span class="sb-label">00a Data Preprocessing</span>\n</a>\n')
    sb_items += (f'<a href="#stage-HostDep" class="sb-item {dep_sb_cls}" '
                 f'title="{_esc_attr(dep_metric)}">\n'
                 '<span class="sb-dot" style="background:#546e7a"></span>\n'
                 '<span class="sb-label">00b Host Depletion</span>\n</a>\n')
    if bbnorm_on:
        sb_items += (f'<a href="#stage-BBNorm" class="sb-item {"sb-pass" if bb_pass else "sb-skip"}" '
                     f'title="{_esc_attr(bb_metric)}">\n'
                     '<span class="sb-dot" style="background:#00838f"></span>\n'
                     '<span class="sb-label">00c BBNorm</span>\n</a>\n')

    # ── KPI 行 (老 kpi-card 版式) ──
    kpi_cards = (
        f'<div class="kpi-card"><div class="kpi-icon">🧬</div><div class="kpi-value">{n}</div>'
        '<div class="kpi-label">Samples</div></div>'
        f'<div class="kpi-card"><div class="kpi-icon">✅</div><div class="kpi-value">{n_pass}</div>'
        '<div class="kpi-label">PASS</div></div>'
        f'<div class="kpi-card"><div class="kpi-icon">⚠️</div><div class="kpi-value">{n_low}</div>'
        '<div class="kpi-label">LOW_RETAINED</div></div>'
        f'<div class="kpi-card"><div class="kpi-icon">🚫</div><div class="kpi-value">{n_bad}</div>'
        '<div class="kpi-label">NO_OUTPUT / NO_DATA</div></div>'
        f'<div class="kpi-card"><div class="kpi-icon">📊</div><div class="kpi-value">{dash(median)}</div>'
        '<div class="kpi-label">Median Retained (%)</div></div>')

    # ── 页尾 Handoff (紧凑 details, 不再是独立侧栏节) ──
    n_final_files = sum(r.get("Final_Files") or 0 for r in rows)
    handoff_html = f'''<details class="stage-table-detail" id="handoff" style="background:var(--card-bg);border-radius:var(--radius);box-shadow:var(--shadow);border-top:none;padding:0;margin-bottom:20px">
<summary style="padding:14px 22px;font-size:13px;font-weight:500;color:var(--text)">Handoff — assembly_ready.list ({n_final_files} files) + preprocessing_summary.tsv</summary>
<div style="padding:2px 22px 16px;font-size:12px;color:var(--muted);line-height:1.9">Downstream handoff list: <code>assembly_ready.list</code> (absolute paths of final reads, one per line; one set per sample;
see the Final_Source column in the per-sample summary: normalized=bbnorm &gt; depleted &gt; clean) ·
<code>preprocessing_summary.tsv</code> (machine-readable version of that table).<br>
Next step: <code>virome_discovery_pipeline/virome_pipeline.py --input_reads &lt;{_D['d_hostdep']} dir&gt;</code></div>
</details>'''

    html = f'''<!DOCTYPE html><html lang="en"><head><meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<!-- impeccable-disable numbered-section-markers,em-dash-overuse: stage IDs are pipeline identifiers; CJK scientific copy -->
<title>Data Preprocessing Report</title>
<style>
{_STYLE_RULES}</style>
</head>
<body>
<div id="lightbox" onclick="closeLightbox()" style="display:none;position:fixed;top:0;left:0;width:100%;height:100%;background:rgba(0,0,0,.85);z-index:999;cursor:zoom-out;align-items:center;justify-content:center">
  <img id="lightbox-img" src="data:image/gif;base64,R0lGODlhAQABAIAAAAAAAP///yH5BAEAAAAALAAAAAABAAEAAAIBRAA7" alt="" style="max-width:95%;max-height:95%;object-fit:contain" onclick="event.stopPropagation()">
</div>
<div class="sidebar">
  <div class="sb-title">Pipeline Modules</div>
  {sb_items}
</div>
<div class="container">
<div class="hero">
  <h1>Data Preprocessing Report</h1>
  <div class="subtitle">clean-data.py → host_depletion.py [→ run_bbnorm.py] — QC · Host Depletion · Normalization</div>
  <div class="gen-time">Generated: {now} &nbsp;|&nbsp; Retained &lt; {warn_retained}% → LOW_RETAINED &nbsp;|&nbsp; {n} samples</div>
</div>
<div class="kpi-row">{kpi_cards}</div>
{stage_cleandata}
{stage_hostdep}
{stage_bbnorm}
{handoff_html}
<div class="footer">
  <strong>MMPV-RNA</strong> — ② Data Preprocessing (data_preprocessing_pipeline)<br>
  Generated by <code>preprocess_report.py</code> &nbsp;|&nbsp; {now}
</div>
</div>
<script>var CD={chart_json};var hasChart={ "true" if chartjs else "false"};</script>
{chartjs}
<script>
{_LIGHTBOX_JS}
if(hasChart&&typeof Chart!=='undefined'&&CD.labels.length){{
new Chart(document.getElementById('chart_s00b'),{{type:'bar',data:{{labels:CD.labels,datasets:[{{label:'Retained (non-host)',data:CD.retained,backgroundColor:'#009E73'}},{{label:'Removed (host+rRNA)',data:CD.removed,backgroundColor:'#D55E00'}}]}},options:{{responsive:true,indexAxis:'y',plugins:{{title:{{display:true,text:'Host Depletion: Reads Retained vs Removed'}},legend:{{display:true}}}},scales:{{x:{{stacked:true,beginAtZero:true,title:{{text:'Reads'}}}},y:{{stacked:true}}}}}}}});
}}else{{var _cb=document.getElementById('chart_s00b');if(_cb&&_cb.parentNode)_cb.parentNode.style.display='none';}}
{_SPY_JS}{_PCT_JS}{pager_js}
if(hasChart&&typeof Chart!=='undefined')setTimeout(function(){{document.querySelectorAll('.pct-btn').forEach(function(b){{toggleStackedPct(b)}})}},100);
</script>
{_TW}</body></html>'''
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
    # 显式传 clean/deplete 目录: 供报告内「每样本明细」小节 (fastp/kraken2 报告
    # 与资源表) 计算相对链接, 不改变任何统计口径
    generate_html(rows, html, args.warn_retained, plots,
                  clean_dir=clean_dir, deplete_dir=deplete_dir)
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
