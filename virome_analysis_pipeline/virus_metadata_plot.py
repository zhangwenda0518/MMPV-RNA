#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
virus_metadata_plot.py — 病毒组 vs 宿主元数据 联合分析与绘图引擎 (发表级重构 v2)

v2 核心改进 (相对 v1 的 135 张原始计数交叉表):
  1. 用整个筛选队列 (denominator) 计算「真实流行率 / 检出率」, 而非原始计数。
     因为公共 SRA 数据里宿主/组织的测序量严重不均衡 (L. barbarum 被测 400+ 次,
     L. chinense 只测 41 次), 原始计数会系统性误导。
  2. 每个 prevalence 带 Wilson 95% 置信区间。
  3. 用 Fisher exact + Benjamini-Hochberg FDR 做「富集/缺失」显著性检验, 标记 * / ** / ***。
  4. 病毒两两共检做「零模型」: 观测 vs 独立假设下的期望, 证明是把协同/拮抗而非巧合共现。
  5. 剔除纯数据采集痕迹 (BioProject / CenterName / TaxID / LibrarySource), 不把它们当生物学信号。
  6. 统一发表级样式: Okabe-Ito 色盲友好离散色 / viridis 热图、白底、300dpi、全英文、
     无文字重叠、顶-N 截断、单双栏友好。
  7. 输出收敛为 ~7 张主图 + 全部 CSV 矩阵, 取代原来 135 张。

用法示例:
  python virus_metadata_plot.py \
      -v 02_filtering/high_conf.meta_adapter.tsv \
      -m sra_rna.data4/global_metadata/Global_Unified_Metadata_Core14.tsv \
      -o metadata_association_pub

兼容旧接口 (-v / -m / -o 不变)。旧行为用 --legacy 可复现。
"""

import os
import sys
import argparse
import warnings
import re
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import fisher_exact, norm

warnings.filterwarnings("ignore", category=UserWarning, module="scipy")
warnings.filterwarnings("ignore", category=FutureWarning, module="seaborn")
try:
    warnings.filterwarnings("ignore", category=pd.errors.SettingWithCopyWarning)
except AttributeError:
    pass

# ----------------------------------------------------------------------------- 样式规范 (发表级)
OKABE_ITO = ["#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
             "#D55E00", "#CC79A7", "#000000", "#999999", "#E69F00",
             "#56B4E9", "#009E73"]
MATPLOTLIB_RCPARAMS = {
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "axes.unicode_minus": False,
    "figure.facecolor": "white",
    "axes.facecolor": "white",
    "savefig.facecolor": "white",
    "axes.edgecolor": "#333333",
    "axes.labelcolor": "#000000",
    "text.color": "#000000",
    "xtick.color": "#000000",
    "ytick.color": "#000000",
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "axes.titlesize": 12,
    "axes.labelsize": 11,
    "axes.linewidth": 0.8,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "legend.frameon": False,
    "figure.dpi": 120,
    "savefig.dpi": 300,
}
plt.rcParams.update(MATPLOTLIB_RCPARAMS)


def safe_filename(name):
    return re.sub(r"[^a-zA-Z0-9_\-]", "_", str(name)).strip("_")


def clean_label(val):
    """去除 _AI 后缀, 统一缺失值 -> Unknown, 规范化中文/英文逗号。"""
    s = str(val).strip()
    s = re.sub(r"_AI$", "", s, flags=re.IGNORECASE)
    s = s.replace("\uff0c", ", ")
    if s.lower() in ["not_provided", "unknown", "nan", "none", "", "<na>", "na"]:
        return "Unknown"
    return s


def simplify_location(loc):
    """把 'China, Ningxia, Yinchuan' 收敛为省份 'Ningxia' 以降低类别数; 国外保留国家。"""
    s = clean_label(loc)
    if s == "Unknown":
        return "Unknown"
    parts = [p.strip() for p in s.split(",")]
    # 中国: 取省级 (第二段); 非中国: 取国家 (第一段)
    if parts and parts[0].lower().startswith("china"):
        return parts[1] if len(parts) > 1 and parts[1] else "China"
    return parts[0] if parts else "Unknown"


def short_host(name):
    """宿主名缩写, 让热图 x 轴标签不重叠。"""
    s = clean_label(name)
    s = s.replace("Lycium chinense var. potaninii", "L. c. var. potaninii")
    s = s.replace("Lycium sp. 'barbarum/ruthenicum'", "L. sp. barb./ruthen.")
    s = s.replace("Alternaria alternata", "A. alternata")
    s = s.replace("Aphis gossypii", "A. gossypii")
    s = s.replace("Fusarium nematophilum", "F. nematophilum")
    s = s.replace("Lycium barbarum", "L. barbarum")
    s = s.replace("Lycium ruthenicum", "L. ruthenicum")
    s = s.replace("Lycium chinense", "L. chinense")
    s = s.replace("Lycium amarum", "L. amarum")
    s = s.replace("Lycium", "L.")
    return s


def wilson_ci(k, n, z=1.96):
    """Wilson score 95% CI for binomial proportion k/n."""
    if n == 0:
        return (0.0, 0.0)
    p = k / n
    denom = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / denom
    margin = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / denom
    return (max(0.0, centre - margin), min(1.0, centre + margin))


def bh_fdr(pvals):
    """Benjamini-Hochberg FDR; 返回独立排序后的 q 值序列 (按输入顺序)。"""
    p = np.asarray(pvals, dtype=float)
    n = p.size
    q = np.ones(n)
    if n == 0:
        return q
    order = np.argsort(p)
    ranks = np.empty(n, dtype=int)
    ranks[order] = np.arange(n)
    ranked = p[order]
    cummin = np.minimum.accumulate(ranked * n / np.arange(1, n + 1))
    cummin = np.minimum.accumulate(cummin[::-1])[::-1]
    min_ = np.minimum(1.0, cummin)
    q[order] = min_
    return q


def significance_stars(q):
    if q < 0.001:
        return "***"
    if q < 0.01:
        return "**"
    if q < 0.05:
        return "*"
    return ""


class VirusMetadataAnalyzer:
    def __init__(self, virus_file, meta_file, out_dir, top_n=12, min_prev=0.02):
        self.virus_file = virus_file
        self.meta_file = meta_file
        self.out_dir = out_dir
        self.top_n = top_n
        self.min_prev = min_prev
        os.makedirs(os.path.join(out_dir, "figures"), exist_ok=True)
        os.makedirs(os.path.join(out_dir, "tables"), exist_ok=True)

        # 生物学相关的组织/来源/地理/宿主特征 (剔除纯数据采集痕迹)
        self.bio_features = ["ScientificName", "Location", "Source", "Tissue",
                             "Age_GrowthStage", "CollectionDate"]
        # 数据出身 (非生物学), 单独导出, 不画进主图
        self.provenance_features = ["BioProject", "CenterName", "TaxID", "LibrarySource"]

    # ---------------- 数据载入与清洗 ----------------
    def _clean_df(self, df):
        df.columns = [str(c).replace("\ufeff", "").strip() for c in df.columns]
        for col in df.columns:
            if df[col].dtype == "object":
                df[col] = df[col].map(clean_label)
        return df

    def load(self):
        # 队列 = 所有被筛选的样本 (分母)
        self.cohort = self._clean_df(pd.read_csv(self.meta_file, sep=None, engine="python",
                                                 skipinitialspace=True))
        self.cohort = self.cohort.rename(columns={c: "Run" for c in self.cohort.columns
                                                  if c.lower() in ["run", "run_accession",
                                                                   "sample", "query_id"]})
        if "Run" not in self.cohort.columns:
            self.cohort = self.cohort.rename(columns={self.cohort.columns[0]: "Run"})
        self.cohort["Run"] = self.cohort["Run"].astype(str).str.strip()
        self.cohort = self.cohort.drop_duplicates(subset=["Run"], keep="last")
        self.n_cohort = len(self.cohort)
        print(f"📥 队列 (分母) 样本数: {self.n_cohort}")

        # 检出 = Run -> 病毒列表 (可一个样本多个病毒)
        inf = self._clean_df(pd.read_csv(self.virus_file, sep=None, engine="python",
                                         skipinitialspace=True))
        inf = inf.rename(columns={c: "Run" for c in inf.columns if c.lower() in
                                  ["sample", "sample_id", "query", "query_id", "id",
                                   "srr", "run_accession"]})
        if "Run" not in inf.columns:
            inf = inf.rename(columns={inf.columns[0]: "Run"})
        virus_col = "Taxonomy" if "Taxonomy" in inf.columns else inf.columns[-1]
        self.infect = inf[["Run", virus_col]].rename(columns={virus_col: "Virus"})
        self.infect["Run"] = self.infect["Run"].astype(str).str.strip()
        self.infect["Virus"] = self.infect["Virus"].map(clean_label)
        self.infect = self.infect[self.infect["Virus"] != "Unknown"]
        print(f"📥 检出记录: {len(self.infect)} (涉及 "
              f"{self.infect['Run'].nunique()} 个样本)")

        # 只有队列内的检出才有效 (已核实全部在内)
        self.infect = self.infect[self.infect["Run"].isin(self.cohort["Run"])]
        self.infect = self.infect.drop_duplicates(subset=["Run", "Virus"])
        # 每样本病毒数 (共感染)
        self.infect["Virus"] = self.infect["Virus"].astype(str)

    # ---------------- 流行率矩阵 ----------------
    def prevalence_matrix(self, feature_col, normalize_province=False):
        """返回 (piv, mat): piv = virus x category 的 prevalence 矩阵; mat 为含计数/CI 的明细表."""
        if feature_col not in self.cohort.columns:
            return None, None

        cohort = self.cohort.copy()
        if normalize_province and feature_col == "Location":
            cohort[feature_col] = cohort[feature_col].map(simplify_location)
        else:
            cohort[feature_col] = cohort[feature_col].map(clean_label)

        # 感染表带上该特征列 (从队列对齐)
        mapping = cohort.set_index("Run")[feature_col]
        inf = self.infect.copy()
        inf[feature_col] = inf["Run"].map(mapping)

        # 病毒排序 + 截断 (先按 >= min_prev 过滤, 再截断到 top_n)
        overall = self.infect["Virus"].value_counts()
        overall = overall[overall.index != "Unknown"]
        keep_viruses = []
        for v, cnt in overall.items():
            if cnt / self.n_cohort >= self.min_prev or len(keep_viruses) < self.top_n:
                keep_viruses.append(v)
            if len(keep_viruses) >= self.top_n:
                break
        inf = inf[inf["Virus"].isin(keep_viruses)]

        cats = cohort[feature_col].value_counts().index.tolist()
        # 类别太碎就只留前 N + Other
        if len(cats) > 12:
            top_cats = cohort[feature_col].value_counts().head(11).index.tolist()
            cohort[feature_col] = cohort[feature_col].where(
                cohort[feature_col].isin(top_cats), "Other")
            inf[feature_col] = inf[feature_col].where(
                inf[feature_col].isin(top_cats), "Other")
            cats = cohort[feature_col].value_counts().index.tolist()

        rows = []
        for v in keep_viruses:
            inf_v = inf[inf["Virus"] == v]
            for c in cats:
                total_c = int((cohort[feature_col] == c).sum())
                inf_c = int((inf_v[feature_col] == c).sum())
                if total_c == 0:
                    continue
                prev = inf_c / total_c
                lo, hi = wilson_ci(inf_c, total_c)
                rows.append({"Virus": v, feature_col: c, "inf": inf_c, "total": total_c,
                             "prevalence": prev, "lo": lo, "hi": hi})
        mat = pd.DataFrame(rows)
        if mat.empty:
            return None, None
        piv = mat.pivot(index="Virus", columns=feature_col, values="prevalence")
        return piv, mat

    # ---------------- 富集 Fisher 检验 ----------------
    def enrichment(self, mat, feature_col):
        """对每 (virus, category) 做 Fisher exact + FDR。返回含 q 的 meta 表。"""
        cohort = self.cohort.copy()
        if feature_col not in cohort.columns:
            return
        cohort[feature_col] = cohort[feature_col].map(clean_label)
        records = []
        for _, row in mat.iterrows():
            v, cat = row["Virus"], row[feature_col]
            inf_c = row["inf"]
            total_c = row["total"]
            inf_v_all = int((self.infect["Virus"] == v).sum())
            inf_not_c = inf_v_all - inf_c
            total_not_c = self.n_cohort - total_c
            uninf_c = total_c - inf_c
            uninf_not_c = total_not_c - inf_not_c
            table = [[inf_c, uninf_c], [inf_not_c, uninf_not_c]]
            if min(inf_c, uninf_c, inf_not_c, uninf_not_c) < 0:
                continue
            # 用 odds ratio 符号判断方向 (OR>1 富集, OR<1 缺失)
            odds, p = fisher_exact(table)
            records.append({"Virus": v, feature_col: cat, "fc": inf_c, "total": total_c,
                            "OR": odds, "p": p})
        en = pd.DataFrame(records)
        if not en.empty:
            en["q"] = bh_fdr(en["p"].values)
        return en

    # ---------------- 共检零模型 ----------------
    def cooccurrence(self):
        top = self.infect["Virus"].value_counts()
        top = top[top.index != "Unknown"].index[:self.top_n]
        n = self.n_cohort
        records = []
        for i, v1 in enumerate(top):
            for v2 in top[i + 1:]:
                runs_v1 = set(self.infect[self.infect["Virus"] == v1]["Run"])
                runs_v2 = set(self.infect[self.infect["Virus"] == v2]["Run"])
                both = len(runs_v1 & runs_v2)
                only_v1 = len(runs_v1 - runs_v2)
                only_v2 = len(runs_v2 - runs_v1)
                neither = n - len(runs_v1 | runs_v2)
                exp = n * (len(runs_v1) / n) * (len(runs_v2) / n)
                odds, p = fisher_exact([[both, only_v1], [only_v2, neither]])
                enr = (both / exp) if exp > 0 else np.nan
                records.append({"VirusA": v1, "VirusB": v2, "observed": both, "expected": exp,
                                "OR": odds, "p": p, "log2_enrich": np.log2(enr) if enr > 0 else np.nan})
        co = pd.DataFrame(records)
        if not co.empty:
            co["q"] = bh_fdr(co["p"].values)
        return co

    # ---------------- 绘图 ----------------
    def _prevalence_heatmap(self, piv, mat, feature_col, fname, title):
        """prevalence 热图: PowerNorm 色彩 + 仅标注非零细胞 + 列头带 n=队列规模。"""
        if piv is None or piv.empty:
            return
        en = self.enrichment(mat, feature_col)
        # 只保留有阳性细胞的行/列, 避免全零噪声
        mat = mat.copy()
        piv = piv.copy()
        # 列按队列规模重排
        sizes = mat.groupby(feature_col)["total"].first()
        cols = sizes.sort_values(ascending=False).index.tolist()
        piv = piv.reindex(columns=[c for c in cols if c in piv.columns])
        # 行按整体流行率重排
        row_sizes = mat.groupby("Virus")["total"].sum()
        row_prev = mat.groupby("Virus")["inf"].sum() / row_sizes
        rows = row_prev.sort_values(ascending=False).index.tolist()
        piv = piv.reindex(index=[r for r in rows if r in piv.index])

        data = piv.values * 100.0
        counts = {}
        for _, r in mat.iterrows():
            counts[(r["Virus"], r[feature_col])] = int(r["inf"])
        # 掩蔽全零细胞 (浅灰), 让信号突出
        masked = np.ma.masked_where(data <= 0.0 + 1e-9, data)
        maxnz = np.nanmax(np.where(data > 0, data, np.nan)) if np.any(data > 0) else 1.0
        vmax = max(5.0, np.ceil(maxnz))

        fig, ax = plt.subplots(figsize=(max(8, 0.7 * piv.shape[1] + 2), 0.42 * piv.shape[0] + 2))
        cmap = plt.get_cmap("viridis").copy()
        cmap.set_bad("#f2f2f2")  # 零/缺失
        norm = matplotlib.colors.PowerNorm(gamma=0.5, vmin=0, vmax=vmax)
        im = ax.imshow(masked, cmap=cmap, aspect="auto", norm=norm)
        ax.set_xticks(np.arange(piv.shape[1]))
        ax.set_yticks(np.arange(piv.shape[0]))
        col_labels = [f"{short_host(c)}\n(n={int(sizes.get(c, 0))})" for c in piv.columns]
        ax.set_xticklabels(col_labels, rotation=45, ha="right")
        ax.set_yticklabels(piv.index)
        # 仅标注非零细胞: value% + (显著性星)
        for r in range(piv.shape[0]):
            for c in range(piv.shape[1]):
                key = (piv.index[r], piv.columns[c])
                cnt = counts.get(key, 0)
                if cnt == 0:
                    continue
                v = data[r, c]
                star = ""
                if en is not None:
                    row_en = en[(en["Virus"] == piv.index[r]) & (en[feature_col] == piv.columns[c])]
                    if not row_en.empty:
                        star = significance_stars(row_en["q"].iloc[0])
                txt = f"{v:.0f}%{star}" if v >= 10 else f"{v:.1f}%{star}"
                ax.text(c, r, txt, ha="center", va="center", fontsize=8,
                        color="black" if v / vmax < 0.55 else "white")
        ax.set_title(title, pad=14)
        cb = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.03)
        cb.set_label("Detection rate (%, sqrt-scaled)")
        plt.tight_layout()
        plt.savefig(os.path.join(self.out_dir, "figures", fname), bbox_inches="tight")
        plt.close()

    def _host_bar(self, piv, mat, fname):
        """主宿主 (Lycium, n>=15) × 核心病毒的检出率图 + 显著性星。
        误差棒在宿主样本量小时会宽到 0-100%, 反而误导, 改用显著性星。"""
        en = self.enrichment(mat, "ScientificName") if "ScientificName" in mat.columns else None
        sizes = mat.groupby("ScientificName")["total"].first()
        # 只保留主 Lycium 宿主 (样本量 >= 15)
        hosts = [h for h in piv.columns if sizes.get(h, 0) >= 15]
        # 只保留 Top 6 病毒
        virs = piv.index[:6]
        x = np.arange(len(hosts))
        width = 0.8 / max(1, len(virs))
        fig, ax = plt.subplots(figsize=(max(8, 0.75 * len(hosts) + 1.5), 5))
        for i, v in enumerate(virs):
            vals, stars = [], []
            vrow = mat[mat["Virus"] == v]
            for h in hosts:
                row = vrow[vrow.ScientificName == h]
                if not row.empty:
                    vals.append(row["prevalence"].iloc[0] * 100)
                    st = ""
                    if en is not None:
                        r_en = en[(en["Virus"] == v) & (en.ScientificName == h)]
                        if not r_en.empty:
                            st = significance_stars(r_en["q"].iloc[0])
                    stars.append(st)
                else:
                    vals.append(0.0); stars.append("")
            xpos = x + i * width
            ax.bar(xpos, vals, width, label=v, color=OKABE_ITO[i % len(OKABE_ITO)],
                   edgecolor="white", linewidth=0.6)
            for xi, val, st in zip(xpos, vals, stars):
                if val > 0.05:
                    ax.text(xi, val + 1.5, st, ha="center", va="bottom", fontsize=9, weight="bold")
        ax.set_xticks(x + width * (len(virs) - 1) / 2)
        ax.set_xticklabels([f"{h}\n(n={int(sizes.get(h,0))})" for h in hosts], rotation=20, ha="right")
        ax.set_xlabel("")  # 组标签已在 xtick
        ax.set_ylabel("Detection rate (%)")
        maxval = 0.0
        for v in virs:
            vrow = mat[mat["Virus"] == v]
            for h in hosts:
                row = vrow[vrow.ScientificName == h]
                if not row.empty:
                    maxval = max(maxval, row["prevalence"].iloc[0] * 100)
        ax.set_ylim(0, max(45, maxval * 1.2))
        ax.set_title("Core virus detection rate across hosts (*q<0.05, **q<0.01, ***q<0.001)", pad=12)
        ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=8, ncol=1)
        plt.tight_layout()
        plt.savefig(os.path.join(self.out_dir, "figures", fname), bbox_inches="tight")
        plt.close()

    def _coinfection(self):
        # 每样本病毒数
        per = self.infect.groupby("Run")["Virus"].nunique()
        per = per[per > 0]
        bins = per.value_counts().sort_index()
        fig, ax = plt.subplots(figsize=(6.5, 4.5))
        ax.bar([f"{k}" for k in bins.index], bins.values, color=OKABE_ITO[0])
        ax.set_xlabel("Number of viruses per sample")
        ax.set_ylabel("Number of samples")
        ax.set_title("Co-infection structure (infected samples)", pad=12)
        for i, v in enumerate(bins.values):
            ax.text(i, v + bins.max() * 0.01, str(int(v)), ha="center", fontsize=9)
        plt.tight_layout()
        plt.savefig(os.path.join(self.out_dir, "figures", "fig05_coinfection_structure.png"),
                    bbox_inches="tight")
        plt.close()
        # 概览率
        total = int(per.sum())
        frac = bins / total * 100
        return bins, frac

    def _cooccurrence_heatmap(self, co):
        if co is None or co.empty:
            return
        virs = sorted(set(co["VirusA"]) | set(co["VirusB"]))
        mat = pd.DataFrame(np.nan, index=virs, columns=virs)
        sig = pd.DataFrame("", index=virs, columns=virs)
        for _, r in co.iterrows():
            mat.loc[r["VirusA"], r["VirusB"]] = r["log2_enrich"]
            mat.loc[r["VirusB"], r["VirusA"]] = r["log2_enrich"]
            st = significance_stars(r["q"])
            sig.loc[r["VirusA"], r["VirusB"]] = st
            sig.loc[r["VirusB"], r["VirusA"]] = st
        fig, ax = plt.subplots(figsize=(0.5 * len(virs) + 2, 0.5 * len(virs) + 2))
        data = mat.values
        im = ax.imshow(data, cmap="RdBu_r", aspect="auto", vmin=-np.nanmax(np.abs(data)),
                       vmax=np.nanmax(np.abs(data)))
        ax.set_xticks(np.arange(len(virs)))
        ax.set_yticks(np.arange(len(virs)))
        ax.set_xticklabels(virs, rotation=45, ha="right", fontsize=8)
        ax.set_yticklabels(virs, fontsize=8)
        for r in range(len(virs)):
            for c in range(len(virs)):
                if not np.isnan(data[r, c]):
                    lbl = ("+" if data[r, c] > 0 else "") + f"{data[r, c]:.1f}\n{sig.iloc[r, c]}"
                    ax.text(c, r, lbl, ha="center", va="center", fontsize=7,
                            color="white" if abs(data[r, c]) > np.nanmax(np.abs(data)) * 0.6 else "black")
        ax.set_title("log2(observed/expected) virus co-occurrence\n(→ enrichment →)",
                     pad=12)
        fig.colorbar(im, ax=ax, fraction=0.035, pad=0.03, label="log2 enrichment")
        plt.tight_layout()
        plt.savefig(os.path.join(self.out_dir, "figures", "fig06_virus_cooccurrence.png"),
                    bbox_inches="tight")
        plt.close()
        return mat

    def _overall_prevalence_dot(self):
        overall = self.infect["Virus"].value_counts()
        overall = overall[overall.index != "Unknown"].head(self.top_n) / self.n_cohort * 100
        fig, ax = plt.subplots(figsize=(6.5, max(4, 0.42 * len(overall) + 1)))
        order = overall.sort_values().index
        vals = [overall[v] for v in order]
        # Wilson CI
        cnt = self.infect["Virus"].value_counts()
        lo = [wilson_ci(cnt.get(v, 0), self.n_cohort)[0] * 100 for v in order]
        hi = [wilson_ci(cnt.get(v, 0), self.n_cohort)[1] * 100 for v in order]
        y = np.arange(len(order))
        ax.errorbar(vals, y, xerr=[np.abs(np.array([v - l for v, l in zip(vals, lo)])),
                                   np.abs(np.array([h - v for v, h in zip(vals, hi)]))],
                    fmt="o", color=OKABE_ITO[3], ecolor="#555", capsize=3)
        ax.set_yticks(y)
        ax.set_yticklabels(order, fontsize=8)
        ax.set_xlabel("Detection rate (% of screened cohort)")
        ax.set_title("Core virome prevalence (top viruses)", pad=12)
        plt.tight_layout()
        plt.savefig(os.path.join(self.out_dir, "figures", "fig07_overall_prevalence.png"),
                    bbox_inches="tight")
        plt.close()

    # ---------------- 主流程 ----------------
    def run(self):
        self.load()

        print("\n📊 模块 1: 核心病毒组 prevalence 热图 ...")
        for feature, fname, title in [
            ("ScientificName", "fig01_prevalence_by_host.png",
             "Virus detection rate by host (across screened cohort)"),
            ("Tissue", "fig02_prevalence_by_tissue.png",
             "Virus detection rate by tissue"),
            ("Location", "fig03_prevalence_by_location.png",
             "Virus detection rate by province"),
        ]:
            piv, mat = self.prevalence_matrix(feature, normalize_province=(feature == "Location"))
            if piv is not None:
                self._prevalence_heatmap(piv, mat, feature, fname, title)
                mat.to_csv(os.path.join(self.out_dir, "tables",
                                        f"prevalence_{feature}.csv"), index=False, encoding="utf-8-sig")

        print("\n📊 模块 2: 宿主检出率柱状图 (带 CI + 显著性) ...")
        piv, mat = self.prevalence_matrix("ScientificName")
        if piv is not None:
            self._host_bar(piv, mat, "fig04_host_range.png")
            en = self.enrichment(mat, "ScientificName")
            if en is not None:
                en.to_csv(os.path.join(self.out_dir, "tables", "enrichment_host.csv"),
                          index=False, encoding="utf-8-sig")

        print("\n📊 模块 3: 共感染结构与病毒两两共检零模型 ...")
        bins, frac = self._coinfection()
        co = self.cooccurrence()
        if co is not None:
            self._cooccurrence_heatmap(co)
            co.to_csv(os.path.join(self.out_dir, "tables", "cooccurrence.csv"),
                      index=False, encoding="utf-8-sig")

        print("\n📊 模块 4: 总体流行率 dot 图 ...")
        self._overall_prevalence_dot()

        # 输出表: 每样本病毒列表 + 队列对齐后的完整表
        out = self.cohort.merge(self.infect.groupby("Run")["Virus"]
                                .apply(lambda x: " + ".join(sorted(set(x)))).reset_index(),
                                on="Run", how="left").fillna({"Virus": "No virus detected"})
        out.to_csv(os.path.join(self.out_dir, "tables", "Sample_Virus_Cohort_Summary.tsv"),
                   sep="\t", index=False, encoding="utf-8-sig")
        print(f"\n🎉 完成! 图 → {os.path.join(self.out_dir, 'figures')} | "
              f"表 → {os.path.join(self.out_dir, 'tables')}")


def main():
    parser = argparse.ArgumentParser(description="病毒组 vs 宿主元数据 发表级分析引擎")
    parser.add_argument("-v", "--virus", required=True, help="病毒检出汇总 (Run->Taxonomy)")
    parser.add_argument("-m", "--meta", default=None, help="队列元数据表 (含全部分母)")
    parser.add_argument("-o", "--outdir", default="./virus_metadata_pub", help="输出目录")
    parser.add_argument("--top-n", type=int, default=12, help="保留核心病毒数 (默认12)")
    parser.add_argument("--min-prev", type=float, default=0.02, help="最低整体流行率 (默认0.02)")
    parser.add_argument("--legacy", action="store_true", help="保留 v1 行为 (原始计数交叉表, 已弃用)")
    args = parser.parse_args()

    if args.legacy:
        print("⚠️ --legacy 指定的旧行为将在下个版本移除。推荐直接使用新增的发表级分析。")
        from importlib import import_module
        legacy = Path(__file__).with_name("virus_metadata_plot_legacy.py")
        if legacy.exists():
            rc = subprocess.call([sys.executable, str(legacy),
                                  "-v", args.virus, "-m", args.meta, "-o", args.outdir])
            sys.exit(rc)
        else:
            print("❌ 未找到 legacy 脚本。")
            sys.exit(1)

    meta_file = args.meta
    if not meta_file or not os.path.exists(meta_file):
        # 自动尝试定位全局元数据
        script_dir = Path(__file__).resolve().parent
        for cand in [
            Path(meta_file).expanduser() if meta_file else None,
            script_dir.parent.parent / "sra_rna.data4" / "global_metadata" / "Global_Unified_Metadata_Core14.tsv",
            Path.home() / "virus" / "data-2026" / "sra_rna.data4" / "global_metadata" / "Global_Unified_Metadata_Core14.tsv",
        ]:
            if cand and cand.exists():
                meta_file = str(cand)
                break
        else:
            print("❌ 未找到队列元数据。请用 -m 指定 852 样本队列元数据表。")
            sys.exit(1)

    analyzer = VirusMetadataAnalyzer(args.virus, meta_file, args.outdir,
                                     top_n=args.top_n, min_prev=args.min_prev)
    analyzer.run()


if __name__ == "__main__":
    main()
