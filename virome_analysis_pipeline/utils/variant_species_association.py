#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
variant_species_association.py — 变异位点 × 宿主物种 关联检验 (Fisher exact + BH)

设计目标
========
对单个病毒的 样本×位点 矩阵, 逐位点检验 "该位点的 ALT 携带是否在某一宿主物种中富集"。
这是 04_post_analysis 的组间统计补件, 与 metadata_association (病毒级) 平级互补。

输入
====
1) 位点矩阵 (二选一布局, 自动识别):
   - snpeff_macro/Matrix_05_Continuous_AF_Matrix.csv   行=Variant_ID, 列=Sample, 值=AF(0填充)
   - vcf_merge/matrix/snp_matrix.tsv                    行=Sample,    列=CHROM_POS_REF_ALT, 值=0/1/NA
2) 元数据 TSV/CSV: 至少含 Run 列与 Host 列 (列名可用参数覆盖)
   宿主按子串归入 --species 给出的物种桶 (默认枸杞属四物种), 未匹配样本剔除并计数。

方法 (算法均为自写标准实现, 无 GPL 代码)
========================================
- 携带定义: 该位点值 >= --af-min (默认 0.05); NA 视为未测, 不计入分母。
- 单侧 Fisher exact (greater): 物种内携带率 vs 其余样本携带率,
  超几何分布尾概率 P(X >= c_g), lgamma 稳定实现。
- BH (Benjamini-Hochberg) 校正: 在"同一物种比较"内部对所有受检位点做 step-up 校正。
- 效应量: Rate_Ratio = (c_g/n_g) / ((C-c_g)/(N-n_g)), 分母为 0 时记 Inf。

输出 (-o 目录)
==============
- Species_Variant_Association_Full.tsv   全部检验行
- SpeciesSpecific_Significant.tsv        BH < alpha 且 Rate_Ratio > 1 的位点
- Association_Summary.tsv                每物种一行汇总
- species_variant_association.png/.pdf   各物种显著位点数条形图 (无显著时跳过)

用法示例
========
python utils/variant_species_association.py \
    -m 04_post_analysis/<virus>/snpeff_macro/Matrix_05_Continuous_AF_Matrix.csv \
    -M sra_rna.data4/global_metadata/Global_Unified_Metadata_Core14.tsv \
    --run-col Run --host-col Host \
    -o 04_post_analysis/<virus>/species_association --label PSTVd
"""

import argparse
import math
import os
import re
import sys

try:
    import pandas as pd
except ImportError:
    sys.exit("[致命错] 需要 pandas (管线 post 阶段环境已自带)")

# Okabe-Ito 主色 (SCI 浅色白底规范)
COLOR_MAIN = "#0072B2"
COLOR_NEUTRAL = "#BFBFBF"


# ==========================================
# 输入校验 (io_utils 风格: 文件名 + 行号)
# ==========================================
class InputError(ValueError):
    pass


def _err(path, lineno, msg):
    raise InputError(f"{os.path.basename(path)}:{lineno}: {msg}")


def load_matrix(path):
    """读入 样本×位点 矩阵, 自动识别两种布局并统一为 行=样本, 列=位点。返回 (df, src_note)。"""
    if not os.path.isfile(path):
        raise InputError(f"矩阵文件不存在: {path}")
    sep = "," if path.lower().endswith(".csv") else "\t"
    df = pd.read_csv(path, sep=sep, index_col=0, dtype=str)

    first_header = str(df.columns[0]).strip()
    idx_name = str(df.index.name or "").strip().lower()

    if first_header.lower() == "sample" or idx_name == "sample":
        # snp_matrix.tsv 布局: 行已是样本; 但首列名是 Sample 时 index_col=0 吃掉的其实是 Sample 列
        # 此时 index.name == 'Sample', 列是位点 → 正确
        note = "layout=samples_rows (snp_matrix 风格)"
    else:
        # Matrix_05 布局: 行=Variant_ID, 列=样本 → 转置
        df = df.T
        df.index.name = "Sample"
        note = "layout=sites_rows (Matrix_05 风格, 已转置)"

    dup_samples = df.index[df.index.duplicated()].unique().tolist()
    if dup_samples:
        _err(path, 0, f"重复样本名: {dup_samples[:5]}")
    dup_sites = df.columns[df.columns.duplicated()].unique().tolist()
    if dup_sites:
        _err(path, 0, f"重复位点ID: {dup_sites[:5]}")

    df = df.apply(pd.to_numeric, errors="coerce")  # 非数值/NA → NaN
    if df.shape[1] == 0:
        _err(path, 1, "未解析到任何位点列")
    return df, note


def parse_site_position(label):
    """从位点标签提取基因组坐标。兼容 CHROM_POS_REF_ALT / GENE_POS_REF>ALT 等;
    取倒数第 3 个下划线字段优先, 失败则回退最后一个纯数字字段; 全失败返回 0。"""
    parts = str(label).strip().split("_")
    if len(parts) >= 3:
        cand = parts[-3]
        if cand.isdigit():
            return int(cand)
    for p in reversed(parts):
        if p.isdigit():
            return int(p)
    m = re.findall(r"\d+", str(label))
    return int(m[-1]) if m else 0


def _resolve_host_col(meta, host_col):
    """宿主列自动检测: 显式指定优先, 否则在常见列名中找。"""
    if host_col:
        return host_col
    for cand in ("Host", "ScientificName", "Species", "host", "Adjusted_Species"):
        if cand in meta.columns:
            return cand
    return ""


def load_metadata(path, run_col, host_col, species_list):
    """读元数据, 返回 (run→物种桶 dict, 统计信息 dict)。校验带文件名+行号。"""
    if not os.path.isfile(path):
        raise InputError(f"元数据文件不存在: {path}")
    sep = "," if path.lower().endswith(".csv") else "\t"
    header_line = open(path, encoding="utf-8-sig", errors="replace").readline()
    if sep not in header_line:
        _err(path, 1, f"分隔符异常: 期望 '{'TAB' if sep == chr(9) else 'comma'}', 请检查扩展名与实际格式")
    meta = pd.read_csv(path, sep=sep, dtype=str)
    if run_col not in meta.columns:
        _err(path, 1, f"缺少必需列 '{run_col}'; 实际列: {list(meta.columns)[:15]}")
    host_col = _resolve_host_col(meta, host_col)
    if not host_col:
        _err(path, 1, f"未找到宿主列 (尝试 Host/ScientificName/Species); 实际列: {list(meta.columns)[:15]}")

    sp_lower = [s.lower() for s in species_list]
    run_map, seen = {}, {}
    unmatched_hosts = 0
    n_rows = 0
    for i, row in meta.iterrows():
        lineno = i + 2
        run = str(row[run_col]).strip()
        host = str(row[host_col]).strip()
        if not run or run.lower() == "nan":
            continue
        n_rows += 1
        if run in seen:
            _err(path, lineno, f"Run 重复出现: '{run}' (首次在第 {seen[run]} 行)")
        seen[run] = lineno
        h_low = host.lower()
        bucket = None
        for name, key in zip(species_list, sp_lower):
            if key in h_low:
                bucket = name
                break
        if bucket is None:
            unmatched_hosts += 1
            continue
        run_map[run] = bucket
    stats = {"rows": n_rows, "mapped": len(run_map), "unmatched_host": unmatched_hosts}
    return run_map, stats


def load_groups(matrix_df, run_map, min_group):
    """把矩阵样本映射到物种桶; 返回 ({物种: [样本]}, 被剔除样本列表)。"""
    groups, dropped = {}, []
    for s in matrix_df.index.astype(str):
        sp = run_map.get(s.strip())
        if sp is None:
            dropped.append(s)
        else:
            groups.setdefault(sp, []).append(s)
    usable = {g: v for g, v in sorted(groups.items()) if len(v) >= min_group}
    too_small = {g: len(v) for g, v in groups.items() if len(v) < min_group}
    return usable, dropped, too_small


# ==========================================
# 统计内核 (纯 stdlib, 自写实现)
# ==========================================
def _logcomb(n, k):
    if k < 0 or k > n:
        return -math.inf
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)


def fisher_right(a, b, c, d):
    """单侧 (greater) Fisher exact。
    2x2:       | 携带ALT | 未携带 |
      物种 g   |    a    |   b    |
      其余样本 |    c    |   d    |
    返回 P(X >= a), 即物种 g 携带数不低于观测值的超几何尾概率。"""
    n_row, A, N = a + b, a + c, a + b + c + d
    if N == 0 or n_row == 0 or n_row == N or A == 0 or A == N:
        return 1.0
    logs = [_logcomb(n_row, k) + _logcomb(N - n_row, A - k) - _logcomb(N, A)
            for k in range(a, min(n_row, A) + 1)]
    m = max(logs)
    return min(1.0, math.exp(m) * sum(math.exp(x - m) for x in logs))


def bh_adjust(pvals):
    """Benjamini-Hochberg step-up。输入原 p 列表, 返回同序 q 值列表。"""
    n = len(pvals)
    if n == 0:
        return []
    order = sorted(range(n), key=lambda i: pvals[i])
    q = [0.0] * n
    running = 1.0
    for rank in range(n, 0, -1):
        idx = order[rank - 1]
        running = min(running, pvals[idx] * n / rank)
        q[idx] = running
    return q


# ==========================================
# 主检验流程
# ==========================================
def run_association(M, groups, af_min, alpha):
    samples = [str(s) for s in M.index]
    site_ids = list(M.columns)
    positions = [parse_site_position(s) for s in site_ids]
    sample_species = {}
    for sp, members in groups.items():
        for s in members:
            sample_species[str(s)] = sp
    g_arr = [sample_species.get(s) for s in samples]

    rows = []
    for j, site in enumerate(site_ids):
        col = M.iloc[:, j]
        tested_mask = col.notna()
        carrier_mask = tested_mask & (col >= af_min)
        C = int(carrier_mask.sum())
        N = int(tested_mask.sum())
        if C == 0 or C == N:
            continue  # 固定位点或零携带, 无信息
        t_arr = tested_mask.to_numpy()
        c_arr = carrier_mask.to_numpy()
        for sp in sorted(groups):
            a = b = c = d = 0
            for i in range(len(samples)):
                if not t_arr[i]:
                    continue
                car = bool(c_arr[i])
                if g_arr[i] == sp:
                    a += car
                    b += not car
                else:
                    c += car
                    d += not car
            n_g, rest = a + b, c + d
            rate_g = a / n_g if n_g else float("nan")
            rate_r = c / rest if rest else float("nan")
            ratio = (rate_g / rate_r) if rate_r > 0 else float("inf")
            p = fisher_right(a, b, c, d)
            rows.append({
                "Site_ID": site, "POS": positions[j], "Species": sp,
                "Carriers_In_Species": a, "Tested_In_Species": n_g,
                "Carrier_Rate_Species_pct": round(rate_g * 100, 2),
                "Carriers_Elsewhere": c, "Tested_Elsewhere": rest,
                "Carrier_Rate_Elsewhere_pct": round(rate_r * 100, 2),
                "Rate_Ratio": ratio,
                "P_Fisher_OneSided": p,
            })

    full = pd.DataFrame(rows)
    if full.empty:
        return full
    for sp, sub in full.groupby("Species", sort=True):
        q = bh_adjust(sub["P_Fisher_OneSided"].tolist())
        full.loc[sub.index, "P_BH"] = [round(x, 6) for x in q]
    full["Rate_Ratio"] = full["Rate_Ratio"].replace(float("inf"), 1e6).round(3)
    full = full.sort_values(["Species", "P_BH", "POS"],
                            na_position="last").reset_index(drop=True)
    return full


def fmt_ratio(v):
    return "Inf" if v >= 1e6 else v


def plot_summary(full, sig, out_png, out_pdf, alpha):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    sp_order = sorted(full["Species"].unique())
    total = [int((full["Species"] == s).sum()) for s in sp_order]
    hits = [int((sig["Species"] == s).sum()) if not sig.empty else 0 for s in sp_order]

    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.bar(sp_order, total, color=COLOR_NEUTRAL, edgecolor="black",
           linewidth=0.6, label="Tested sites")
    ax.bar(sp_order, hits, color=COLOR_MAIN, edgecolor="black",
           linewidth=0.6, label=f"Significant (BH < {alpha})")
    for x, h in zip(range(len(sp_order)), hits):
        if h:
            ax.text(x, h, str(h), ha="center", va="bottom",
                    fontsize=10, fontweight="bold", color="#111111")
    ax.set_ylabel("Number of variant sites", fontsize=11)
    ax.set_xlabel("Host species", fontsize=11)
    ax.set_title("Host-species-associated variants per comparison",
                 fontsize=12, fontweight="bold")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend(frameon=False, fontsize=9)
    plt.tight_layout()
    fig.savefig(out_png, dpi=300, facecolor="white")
    fig.savefig(out_pdf, facecolor="white")
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(
        description="变异位点 × 宿主物种 Fisher exact + BH 关联检验 (单病毒)",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("-m", "--matrix", required=True,
                    help="样本×位点矩阵 (Matrix_05_Continuous_AF_Matrix.csv 或 snp_matrix.tsv)")
    ap.add_argument("-M", "--metadata", required=True, help="元数据 TSV/CSV (需含 Run 与 Host 列)")
    ap.add_argument("-o", "--outdir", required=True, help="输出目录")
    ap.add_argument("--run-col", default="Run", help="元数据中 Run ID 列名")
    ap.add_argument("--host-col", default="",
                    help="元数据中宿主列名 (默认自动检测: Host/ScientificName/Species)")
    ap.add_argument("--species", nargs="+",
                    default=["Lycium barbarum", "Lycium ruthenicum",
                             "Nicotiana chinense", "Nicotiana amarum"],
                    help="宿主物种桶 (按子串大小写不敏感匹配 host-col 值)")
    ap.add_argument("--af-min", type=float, default=0.05, help="ALT 携带判定阈值")
    ap.add_argument("--min-group", type=int, default=3, help="物种参与检验的最小样本数")
    ap.add_argument("--alpha", type=float, default=0.05, help="BH 显著性阈值")
    ap.add_argument("--label", default="", help="病毒名标签 (用于标题与日志)")
    ns = ap.parse_args()

    os.makedirs(ns.outdir, exist_ok=True)
    tag = ns.label or os.path.basename(os.path.dirname(os.path.dirname(os.path.abspath(ns.matrix))))
    print(f"\n{'=' * 64}\n 变异×宿主物种关联检验 | {tag}\n{'=' * 64}")

    try:
        M, layout = load_matrix(ns.matrix)
        print(f"[矩阵] {M.shape[0]} 样本 × {M.shape[1]} 位点 ({layout})")

        run_map, mstat = load_metadata(ns.metadata, ns.run_col, ns.host_col, ns.species)
        print(f"[元数据] 读入 {mstat['rows']} 行, 映射成功 {mstat['mapped']}, "
              f"宿主无法归类剔除 {mstat['unmatched_host']}")

        groups, dropped, too_small = load_groups(M, run_map, ns.min_group)
        if dropped:
            print(f"[警告] {len(dropped)} 个矩阵样本不在元数据中, 剔除: {dropped[:8]}"
                  f"{' ...' if len(dropped) > 8 else ''}")
        if too_small:
            print(f"[警告] 以下物种样本数 < {ns.min_group}, 不参与检验: {too_small}")
        if not groups:
            sys.exit("[致命错] 无任何物种满足最小样本数, 请检查 --run-col/--host-col/--species 或矩阵样本命名")
        print(f"[分组] 检验对象: " + ", ".join(f"{g}(n={len(v)})" for g, v in groups.items()))

        full = run_association(M, groups, ns.af_min, ns.alpha)
        if full.empty:
            print("[结果] 所有位点均固定或零携带, 无可检验变异。")
            return
        full.insert(0, "Label", tag)
        full_out = full.copy()
        full_path = os.path.join(ns.outdir, "Species_Variant_Association_Full.tsv")
        full_out.assign(Rate_Ratio=full_out["Rate_Ratio"].map(fmt_ratio)).to_csv(
            full_path, sep="\t", index=False)

        sig = full_out[(full_out["P_BH"] < ns.alpha) & (full_out["Rate_Ratio"] > 1)]
        sig_path = os.path.join(ns.outdir, "SpeciesSpecific_Significant.tsv")
        sig.assign(Rate_Ratio=sig["Rate_Ratio"].map(fmt_ratio)).to_csv(
            sig_path, sep="\t", index=False)

        summ_rows = []
        for sp in sorted(groups):
            tot = int((full_out["Species"] == sp).sum())
            h = int((sig["Species"] == sp).sum()) if not sig.empty else 0
            top = sig[sig["Species"] == sp].head(3)["Site_ID"].tolist() if h else []
            summ_rows.append({"Species": sp, "N_Samples": len(groups[sp]),
                              "Sites_Tested": tot, "Sites_Significant": h,
                              "Top_Hits": "; ".join(top)})
        pd.DataFrame(summ_rows).to_csv(
            os.path.join(ns.outdir, "Association_Summary.tsv"), sep="\t", index=False)

        png = os.path.join(ns.outdir, "species_variant_association.png")
        pdf = os.path.join(ns.outdir, "species_variant_association.pdf")
        plot_summary(full_out, sig, png, pdf, ns.alpha)

        print(f"[结果] 受检组合 {len(full_out)} 行; 显著位点 "
              f"{len(sig)} 个 (BH<{ns.alpha}, Rate_Ratio>1)")
        for _, r in pd.DataFrame(summ_rows).iterrows():
            print(f"   - {r['Species']}: {r['Sites_Significant']}/{r['Sites_Tested']}"
                  + (f"  top: {r['Top_Hits']}" if r["Top_Hits"] else ""))
        print(f"[输出] {full_path}\n[输出] {sig_path}\n[图]   {png}")
    except InputError as e:
        sys.exit(f"[输入错误] {e}")


if __name__ == "__main__":
    main()
