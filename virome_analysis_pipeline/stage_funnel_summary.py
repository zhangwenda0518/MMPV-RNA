#!/usr/bin/env python3
"""
stage_funnel_summary.py — 已知病毒五阶段检出→全长漏斗综合表与图
==================================================================
串联五个阶段的产物, 回答"病毒从原始检出到最后可用全长序列到底留存多少":

  1. 原始检出     01_Detection/summary/all_viruses.raw.tsv          (Poisson 去假阳性前)
  2. best 检出    01_Detection/summary/all_viruses.best.summary.tsv (Poisson 过滤 + 物种最优代表)
  3. 过滤检出     02_Filtering/high_conf.summary.tsv                (Stage2 多维高置信过滤)
  4. 全长组装     05_Assembly/{Virus}/{Sample}/Global_Evolution_Stats.tsv (de novo 组装终态)
  5. 组装过滤     06_Extraction/{Virus}/{Sample}.full.fasta         (最长 contig 提取 + 长度/N 过滤后留存)

输出 (默认 09_Report/stage_funnel/):
  stage_funnel_summary.tsv    每阶段一行: 记录数(病毒×样本)/物种数/样本数/留存率
  virus_stage_matrix.tsv      每物种一行: 五阶段各自的病毒×样本记录数 + 组装长度/N 含量
  stage_funnel.png/.pdf       连续梯形漏斗 (记录数 + 物种数双面板, 只标数量)
  virus_stage_matrix.png/.pdf 物种×阶段留存气泡点阵 (圆=记录数, 灰空格=淘汰, Top N 物种)

用法:
  python stage_funnel_summary.py -d <pipeline 根目录> [-o 输出目录] [--no_plots] [--top_n 40]

供 generate_pipeline_report.py 在 Stage 9 报告时 import 调用 (build 函数)。
与 virus_detection_summary.py 的关系: 那是 panvirome 风格的旧综合表 (无图/无报告集成),
本脚本面向已知病毒管线报告, 覆盖完整五阶段口径。
"""
import argparse
import os
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

# 跨管线统一 I/O 布局 (mmpv_common/, 仓库根): 目录名随 MMPV_IO_LAYOUT 解析
# (编排器已 normalize 环境变量, 子进程导入本模块时快照即正确布局)
_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)
try:
    from mmpv_common.io_layout import layout_dirs as _layout_dirs
    _D = _layout_dirs(os.environ.get("MMPV_IO_LAYOUT"))
except Exception:
    _D = {"a_detect": "01_Detection", "a_filter": "02_Filtering",
          "a_assembly": "05_Assembly", "a_extract": "06_Extraction",
          "a_report": "09_Report"}

ACC_RE = re.compile(r'([A-Z]{1,4}_?\d{5,9}\.\d{1,2})')

# (key, 英文图标签, 中文说明) — 顺序即漏斗顺序
STAGES = [
    ("raw",       "Raw Detection",            "原始检出 (Poisson 过滤前)"),
    ("best",      "Best (Poisson-filtered)",  "best 检出 (Poisson 过滤+物种最优代表)"),
    ("highconf",  "High-Conf Filter",         "过滤检出 (多维高置信过滤)"),
    ("assembly",  "Full-Length Assembly",     "全长组装 (de novo, 含终态统计)"),
    ("extracted", "Extracted Full-Length",    "全长组装过滤 (提取+长度/N 过滤留存)"),
]
STAGE_KEYS = [s[0] for s in STAGES]


def _extract_acc(text):
    m = ACC_RE.search(str(text))
    return m.group(1) if m else (str(text).strip() or None)


def _read_tsv(path):
    p = Path(path)
    if not p.is_file():
        return None
    try:
        df = pd.read_csv(p, sep="\t", dtype=str)
        return df if len(df) > 0 else None
    except Exception:
        return None


def _pairs(df, acc_cols, sample_col="Sample"):
    """从 summary 表提取 (acc, sample) 记录集合; 返回 DataFrame[Acc, Sample] 或 None"""
    if df is None:
        return None
    acc_col = next((c for c in acc_cols if c in df.columns), None)
    if not acc_col or sample_col not in df.columns:
        return None
    out = pd.DataFrame({"Acc": df[acc_col].astype(str), "Sample": df[sample_col].astype(str)})
    return out.drop_duplicates()


# ── 五阶段收集器 ─────────────────────────────────────────────

def collect_raw(base):
    """raw 表的 acc 列是 Accession (best 前的形态); 同时返回自身 Species 列的映射。"""
    df = _read_tsv(Path(base) / _D['a_detect'] / "summary" / "all_viruses.raw.tsv")
    pairs = _pairs(df, ["Accession", "Rep_Accesion", "Rep_Accession"])
    acc2sp = {}
    if df is not None:
        acc_col = next((c for c in ["Accession", "Rep_Accesion", "Rep_Accession"] if c in df.columns), None)
        sp_col = next((c for c in ["Species", "Adjusted_Species", "Species_NCBI"] if c in df.columns), None)
        if acc_col and sp_col:
            for acc, sp in df[[acc_col, sp_col]].drop_duplicates().values:
                sp_s = str(sp).strip()
                if sp_s and sp_s.lower() not in ("nan", "none", "-", "unannotated"):
                    acc2sp[str(acc)] = sp_s
    return pairs, acc2sp


def collect_best(base):
    return _read_best_like(Path(base) / _D['a_detect'] / "summary" / "all_viruses.best.summary.tsv")


def _read_best_like(path):
    df = _read_tsv(path)
    pairs = _pairs(df, ["Rep_Accession", "Accession"])
    if df is None:
        return None, {}
    sp_col = next((c for c in ["Adjusted_Species", "Species_NCBI", "Species"] if c in df.columns), None)
    acc_col = next((c for c in ["Rep_Accession", "Accession"] if c in df.columns), None)
    acc2sp = {}
    if acc_col and sp_col:
        for acc, sp in df[[acc_col, sp_col]].drop_duplicates().values:
            sp_s = str(sp).strip()
            acc2sp[str(acc)] = sp_s if sp_s and sp_s.lower() not in ("nan", "none", "-") else str(acc)
    return pairs, acc2sp


def collect_highconf(base):
    return _read_best_like(Path(base) / _D['a_filter'] / "high_conf.summary.tsv")


def collect_assembly(base):
    """扫描 05_Assembly/{Virus}/{Sample}/: 有 Global_Evolution_Stats.tsv 视为组装成功。

    返回 (records, attempted): records = DataFrame[Acc, Sample, Final_len, N_count, N50],
    attempted = 成功∪失败的全部 (acc, sample) 目录数。
    """
    asm_dir = Path(base) / _D['a_assembly']
    rows, attempted = [], 0
    if not asm_dir.is_dir():
        return None, 0
    for vdir in sorted(asm_dir.iterdir()):
        if not vdir.is_dir() or vdir.name.startswith("run_"):
            continue
        acc = _extract_acc(vdir.name)
        for sdir in sorted(vdir.iterdir()):
            if not sdir.is_dir():
                continue
            attempted += 1
            stats_f = sdir / "Global_Evolution_Stats.tsv"
            sdf = _read_tsv(stats_f)
            if sdf is None:
                continue
            last = sdf.iloc[-1]
            def _num(col):
                try:
                    return float(last.get(col, 0) or 0)
                except (TypeError, ValueError):
                    return 0.0
            rows.append({"Acc": acc, "Sample": sdir.name,
                         "Final_len": _num("Total_Length"),
                         "N_count": int(_num("N_Count")),
                         "N50": _num("Contig_N50")})
    return (pd.DataFrame(rows) if rows else None), attempted


def collect_extract(base):
    """扫描 06_Extraction/{Virus}/*.fasta: 留存文件即通过长度/N 过滤的全长序列。

    样本名从 FASTA 头恢复 (提取器写入的 header = 原样本目录名)。
    返回 DataFrame[Acc, Sample, Len, N_pct] 或 None。
    """
    ext_dir = Path(base) / _D['a_extract']
    rows = []
    if not ext_dir.is_dir():
        return None
    for vdir in sorted(ext_dir.iterdir()):
        if not vdir.is_dir() or vdir.name.startswith("run_"):
            continue
        acc = _extract_acc(vdir.name)
        fas = list(vdir.glob("*.full.fasta")) or list(vdir.glob("*.fasta")) + list(vdir.glob("*.fa"))
        for fa in sorted(fas):
            sample, seqlen, ncnt = None, 0, 0
            try:
                cur = []
                with open(fa, "r", encoding="utf-8", errors="replace") as f:
                    for line in f:
                        line = line.strip()
                        if line.startswith(">"):
                            # 提取产物为单记录 fasta; header 即原样本目录名
                            if sample is None:
                                parts = line[1:].split()
                                sample = parts[0] if parts else fa.stem
                        elif line:
                            cur.append(line)
                seq = "".join(cur)
                seqlen, ncnt = len(seq), seq.upper().count("N")
                if sample is None:
                    sample = fa.stem  # 无 header 的裸文件, 退回文件名
            except Exception:
                continue
            if seqlen == 0:
                continue
            rows.append({"Acc": acc, "Sample": sample, "Len": seqlen,
                         "N_pct": round(ncnt / seqlen * 100, 2)})
    return pd.DataFrame(rows) if rows else None


# ── 汇总 ─────────────────────────────────────────────────────

def build(base, outdir=None, make_plots=True, top_n=40, dpi=300):
    """主入口: 收集五阶段 → 写两张 TSV → 出两张图。

    返回 dict(stages=漏斗统计 DataFrame, matrix=物种矩阵 DataFrame,
              outdir=输出目录, funnel_png/heatmap_png=图路径或 None)。
    任何阶段缺失都不报错, 该阶段记 0/None, 报告端可据此降级展示。
    """
    base = Path(base)
    if outdir is None:
        outdir = base / _D['a_report'] / "stage_funnel"
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    raw, acc2sp_raw = collect_raw(base)
    best, acc2sp_best = collect_best(base)
    hc, acc2sp_hc = collect_highconf(base)
    asm, asm_attempted = collect_assembly(base)
    ext = collect_extract(base)

    # acc → 物种名: best 优先, high_conf 次之, raw 表 Species 列兜底
    acc2sp = dict(acc2sp_raw)
    acc2sp.update(acc2sp_hc)
    acc2sp.update(acc2sp_best)

    def _sp_of(acc):
        return acc2sp.get(acc, acc)

    def _stage_records(df):
        if df is None or len(df) == 0:
            return pd.DataFrame(columns=["Acc", "Sample"])
        return df

    stage_dfs = {
        "raw": _stage_records(raw),
        "best": _stage_records(best),
        "highconf": _stage_records(hc),
        "assembly": _stage_records(asm),
        "extracted": _stage_records(ext),
    }

    # ── 漏斗统计表 ──
    rows = []
    prev_records = None
    raw_records = None
    for key, en, cn in STAGES:
        df = stage_dfs[key]
        n_rec = len(df)
        n_sp = df["Acc"].map(_sp_of).nunique() if n_rec else 0
        n_smp = df["Sample"].nunique() if n_rec else 0
        # 上一级/原始级为 0 时留存率不可计算, 留空 (避免 nan)
        ret_prev = (n_rec / prev_records * 100) if prev_records else None
        ret_raw = (n_rec / raw_records * 100) if raw_records else None
        rows.append({"Stage": en, "Stage_CN": cn, "Records": n_rec,
                     "Species": int(n_sp), "Samples": int(n_smp),
                     "Retention_vs_Previous(%)": round(ret_prev, 1) if ret_prev is not None else "",
                     "Retention_vs_Raw(%)": round(ret_raw, 1) if ret_raw is not None else ""})
        if key == "raw":
            raw_records = n_rec
        prev_records = n_rec
    summary_df = pd.DataFrame(rows)
    summary_df.to_csv(outdir / "stage_funnel_summary.tsv", sep="\t", index=False)

    # ── 物种 × 阶段矩阵 ──
    all_pairs = pd.concat([stage_dfs[k].assign(_stage=k)[["Acc", "Sample", "_stage"]]
                           for k in STAGE_KEYS if len(stage_dfs[k]) > 0],
                          ignore_index=True) if any(len(stage_dfs[k]) > 0 for k in STAGE_KEYS) else None
    mat_rows = []
    if all_pairs is not None:
        all_pairs = all_pairs.drop_duplicates()
        all_pairs["Species"] = all_pairs["Acc"].map(_sp_of)
        cnt = (all_pairs.pivot_table(index="Species", columns="_stage", values="Sample",
                                     aggfunc="count", fill_value=0)
               .reindex(columns=STAGE_KEYS, fill_value=0).astype(int))
        # acc / 组装 / 提取度量
        acc_map = all_pairs.drop_duplicates("Species").set_index("Species")["Acc"].to_dict()
        for sp, r in cnt.iterrows():
            acc = acc_map.get(sp, sp)
            row = {"Species": sp, "Rep_Accession": acc}
            row.update({k: int(r[k]) for k in STAGE_KEYS})
            if asm is not None and len(asm) and acc in set(asm["Acc"]):
                g = asm[asm["Acc"] == acc]
                row["assembly_final_len_avg"] = int(g["Final_len"].mean())
                row["assembly_N50_max"] = int(g["N50"].max())
            else:
                row["assembly_final_len_avg"] = ""
                row["assembly_N50_max"] = ""
            if ext is not None and len(ext) and acc in set(ext["Acc"]):
                g = ext[ext["Acc"] == acc]
                row["extracted_len_avg"] = int(g["Len"].mean())
                row["extracted_N_pct_avg"] = round(g["N_pct"].mean(), 2)
                row["full_length"] = "Y"
            else:
                row["extracted_len_avg"] = ""
                row["extracted_N_pct_avg"] = ""
                row["full_length"] = "N"
            mat_rows.append(row)
    matrix_df = pd.DataFrame(mat_rows)
    if len(matrix_df) > 0:
        matrix_df = matrix_df.sort_values(["extracted", "raw"], ascending=False).reset_index(drop=True)
        cols = ["Species", "Rep_Accession"] + STAGE_KEYS + [
            "assembly_final_len_avg", "assembly_N50_max", "extracted_len_avg",
            "extracted_N_pct_avg", "full_length"]
        matrix_df = matrix_df[[c for c in cols if c in matrix_df.columns]]
        matrix_df.to_csv(outdir / "virus_stage_matrix.tsv", sep="\t", index=False)

    # ── 出图 ──
    funnel_png = heatmap_png = None
    if make_plots:
        try:
            funnel_png = plot_funnel(summary_df, outdir, dpi=dpi)
            heatmap_png = plot_matrix(matrix_df, outdir, dpi=dpi, top_n=top_n)
        except Exception as e:
            print(f"[warn]绘图失败 (表格已写出): {e}")

    # 控制台摘要
    print("=" * 62)
    print(" 五阶段漏斗: 原始检出 → best → 高置信过滤 → 全长组装 → 提取留存")
    print("=" * 62)
    for _, r in summary_df.iterrows():
        ret = r["Retention_vs_Previous(%)"]
        ret_s = f" ← {ret}%" if ret != "" else ""
        print(f"  {r['Stage']:<26} {r['Records']:>7} records | {r['Species']:>4} species{ret_s}")
    print(f"  输出 -> {outdir}")

    return {"stages": summary_df, "matrix": matrix_df, "outdir": str(outdir),
            "funnel_png": funnel_png, "heatmap_png": heatmap_png,
            "assembly_attempted": asm_attempted}


# ── 绘图 (matplotlib, 与管线出图风格一致: PDF+PNG 300dpi) ──────

_SHORT = {"raw": "Raw\nDetection", "best": "Best\n(Poisson)", "highconf": "High-Conf\nFilter",
          "assembly": "Full-Length\nAssembly", "extracted": "Extracted\nFull-Length"}


def plot_funnel(summary_df, outdir, dpi=300):
    """连续梯形漏斗: 层间连接带成整体漏斗形; 左面板病毒×样本记录数, 右面板物种数。

    层内只标数量, 不再标注留存率 (表里有, 图面保持干净)。
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Polygon

    stages = summary_df["Stage"].tolist()
    _LABELS = {"Raw Detection": "Raw\nDetection",
               "Best (Poisson-filtered)": "Best\n(Poisson)",
               "High-Conf Filter": "High-Conf\nFilter",
               "Full-Length Assembly": "Full-Length\nAssembly",
               "Extracted Full-Length": "Extracted\nFull-Length"}
    labels = [_LABELS.get(s, s) for s in stages]
    n = len(stages)
    blues = plt.cm.Blues(np.linspace(0.88, 0.42, n))

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    for ax, col, title in zip(axes, ["Records", "Species"],
                              ["Records (virus × sample)", "Virus species"]):
        vals = summary_df[col].astype(float).values
        vmax = max(vals.max(), 1)
        row_h, gap = 0.62, 0.55
        ys_top = np.zeros(n)

        def _digit_color(fc):
            r, g, b = fc[0], fc[1], fc[2]
            return "white" if 0.299 * r + 0.587 * g + 0.114 * b < 0.5 else "#17395c"

        for i in range(n):
            y_i = -i * (row_h + gap)
            ys_top[i] = y_i
            ax.add_patch(plt.Rectangle((-vals[i] / 2, y_i - row_h), vals[i], row_h,
                                       facecolor=blues[i], edgecolor="white",
                                       linewidth=1.2, zorder=3))
            ax.text(0, y_i - row_h / 2, f"{int(vals[i]):,}", ha="center", va="center",
                    fontsize=12, fontweight="bold", color=_digit_color(blues[i]), zorder=4)
            if i > 0:  # 层间梯形连接带 (淡色, 只造形不标数)
                poly = Polygon([
                    (-vals[i - 1] / 2, ys_top[i - 1] - row_h),
                    (vals[i - 1] / 2, ys_top[i - 1] - row_h),
                    (vals[i] / 2, y_i), (-vals[i] / 2, y_i)],
                    closed=True, facecolor=blues[i], alpha=0.28,
                    edgecolor="none", zorder=2)
                ax.add_patch(poly)
        # ylim 下界必须盖住最后一层条的底边, 否则末层被裁 (视觉上变矮)
        ax.set_ylim(-(n - 1) * (row_h + gap) - row_h - 0.18, row_h * 0.6)
        ax.set_xlim(-vmax * 0.62, vmax * 0.9)
        ax.set_yticks(ys_top - row_h / 2)
        ax.set_yticklabels(labels, fontsize=10)
        ax.xaxis.set_visible(False)
        for spine in ax.spines.values():
            spine.set_visible(False)
        ax.set_title(title, fontweight="bold", fontsize=13, pad=10)
    fig.suptitle("Virus Detection & Assembly Funnel — Raw to Full-Length",
                 fontweight="bold", fontsize=15, y=1.02)
    plt.tight_layout()
    fn = os.path.join(str(outdir), "stage_funnel")
    fig.savefig(fn + ".pdf", dpi=dpi, bbox_inches="tight")
    fig.savefig(fn + ".png", dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return fn + ".png"


def plot_matrix(matrix_df, outdir, dpi=300, top_n=40):
    """物种 × 阶段留存气泡点阵: 有检出=圆泡 (大小/深浅=记录数), 淘汰=灰空格。

    Top N 物种 (按 extracted/raw 降序); 物种多时点阵比色块热图更耐读。
    """
    if matrix_df is None or len(matrix_df) == 0:
        return None
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    vals = matrix_df[STAGE_KEYS].astype(int).values
    labels = matrix_df["Species"].astype(str)
    note = ""
    if len(matrix_df) > top_n:
        vals = vals[:top_n]
        labels = labels.iloc[:top_n]
        note = f" (top {top_n} of {len(matrix_df)} species by retention)"

    nr, nc = vals.shape
    fig, ax = plt.subplots(figsize=(8.2, max(3.2, 0.62 * nr + 1.6)))
    ax.set_facecolor("#f5f5f5")
    vmax = max(vals.max(), 1)
    for i in range(nr):
        for j in range(nc):
            v = vals[i, j]
            if v > 0:
                size = 340 + 1050 * (v / vmax)
                ax.scatter(j, i, s=size, c=[plt.cm.Blues(0.35 + 0.55 * v / vmax)],
                           edgecolors="white", linewidths=1.4, zorder=3)
                ax.text(j, i, str(v), ha="center", va="center", fontsize=9,
                        fontweight="bold",
                        color="white" if v / vmax > 0.55 else "#1a3a5c", zorder=4)
    ax.set_xticks(range(nc))
    ax.set_xticklabels([_SHORT[k] for k in STAGE_KEYS], fontsize=9.5)
    ax.set_yticks(range(nr))
    ax.set_yticklabels(labels, fontsize=8, style="italic")
    ax.set_ylim(nr - 0.5 + 0.45, -0.68)  # 行序自上而下 + 底部留气泡余量
    ax.set_xlim(-0.6, nc - 0.4)
    ax.xaxis.set_ticks_position("top")
    for j in np.arange(0.5, nc - 1, 1):
        ax.axvline(j, color="white", lw=2.2)
    for i in np.arange(0.5, nr - 1, 1):
        ax.axhline(i, color="white", lw=2.2)
    for spine in ax.spines.values():
        spine.set_color("#dddddd")
    ax.set_title(f"Per-Virus Retention Across Stages — bubble = records (virus × sample){note}",
                 fontweight="bold", fontsize=12, pad=30)
    plt.tight_layout()
    fn = os.path.join(str(outdir), "virus_stage_matrix")
    fig.savefig(fn + ".pdf", dpi=dpi, bbox_inches="tight")
    fig.savefig(fn + ".png", dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return fn + ".png"


def main():
    ap = argparse.ArgumentParser(description="已知病毒五阶段漏斗综合表与图")
    ap.add_argument("-d", "--dir", required=True, help="管线输出根目录 (含 01_Detection 等)")
    ap.add_argument("-o", "--output", default=None,
                    help=f"输出目录 (默认: <dir>/{_D['a_report']}/stage_funnel)")
    ap.add_argument("--no_plots", action="store_true", help="只出表格, 不绘图")
    ap.add_argument("--top_n", type=int, default=40, help="热图最多展示的物种数 (默认 40)")
    ap.add_argument("--dpi", type=int, default=300, help="图片分辨率 (默认 300)")
    args = ap.parse_args()

    build(args.dir, args.output, make_plots=not args.no_plots,
          top_n=args.top_n, dpi=args.dpi)


if __name__ == "__main__":
    main()
