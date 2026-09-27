# -*- coding: utf-8 -*-
"""AF 兜底路径冒烟测试: build_matrix_from_af_tables / qc_stats_from_af_tables /
parse_site_pos / 距离矩阵非零性。"""
import os
import shutil
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)

import virus_vcf_pipeline as vp

TMP = os.path.join(HERE, "_tmp_af")
shutil.rmtree(TMP, ignore_errors=True)

failures = []


def check(name, cond, detail=""):
    print(f"[{'PASS' if cond else 'FAIL'}] {name}" + (f" | {detail}" if detail else ""))
    if not cond:
        failures.append(name)


# ---- 构造伪输入目录: <Run>_<Contig>.filtered.vcf + <Run>.<Contig>.allele_frequencies.tsv ----
def make_sample(run_contig, run, contig, rows):
    d = os.path.join(TMP, f"{run}_{contig}")
    os.makedirs(d)
    open(os.path.join(d, f"{run}_{contig}.filtered.vcf"), "w").close()
    with open(os.path.join(d, f"{run}.{contig}.allele_frequencies.tsv"), "w") as f:
        f.write("CHROM\tPOS\tREF\tALT\tALT_FREQ\n")
        for r in rows:
            f.write("\t".join(str(x) for x in [contig] + list(r)) + "\n")


make_sample("S_A_C1.1", "S_A", "C1.1", [
    ("10", "A", "G", 0.40),          # 高 AF → 1
    ("20", "C", "T", 0.02),          # 有记录低 AF → 显式 0
    ("30", "G", "A", 0.60),          # 多等位之一
    ("30", "G", "GA", 0.25),         # 多等位之二 (插入)
    ("40", "CN", "C", 0.90),         # 缺失
])
make_sample("S_B_C1.1", "S_B", "C1.1", [
    ("20", "C", "T", 0.80),          # S_A 低此处, S_B 高
    ("50", "T", "C", 0.10),          # 仅 S_B 有
])
os.makedirs(os.path.join(TMP, "S_C_C1.1"))
open(os.path.join(TMP, "S_C_C1.1", "S_C_C1.1.filtered.vcf"), "w").close()  # 无 AF 表 → 跳过

vcfs = []
for base, dirs, files in os.walk(TMP):
    for fn in files:
        if fn.endswith(".filtered.vcf"):
            vcfs.append(os.path.join(base, fn))

M, samples, labels = vp.build_matrix_from_af_tables(
    vcfs, af_min=0.05,
    out_bin_tsv=os.path.join(TMP, "snp_matrix.tsv"),
    out_af_tsv=os.path.join(TMP, "af_matrix.tsv"))

check("无表样本被跳过", samples == ["S_A", "S_B"], str(samples))
check("位点数 = 6 (含多等位拆分)", len(labels) == 6, str(labels))
iA, iB = samples.index("S_A"), samples.index("S_B")
lab = {l: k for k, l in enumerate(labels)}
check("高 AF 记 1", M[iA, lab["C1.1_10_A_G"]] == 1.0)
check("低 AF 记录显式 0 (修复 NA 灾难)", M[iA, lab["C1.1_20_C_T"]] == 0.0)
check("未记录位点显式 0", M[iB, lab["C1.1_10_A_G"]] == 0.0)
check("多等位 ALT 各占一列",
      "C1.1_30_G_A" in lab and "C1.1_30_G_GA" in lab)
check("插入/缺失入列", "C1.1_40_CN_C" in lab and "C1.1_50_T_C" in lab)
check("连续 AF 矩阵保真", abs(np.loadtxt(os.path.join(TMP, "af_matrix.tsv"),
      delimiter="\t", skiprows=1, usecols=1 + lab["C1.1_10_A_G"])[0] - 0.40) < 1e-9)

# ---- 距离矩阵不再全零 ----
D = vp.pairwise_distance_matrix(M, metric="both")
check("Hamming 非零", np.isfinite(D["hamming"][iA, iB]) and D["hamming"][iA, iB] > 0,
      f"d={D['hamming'][iA, iB]:.3f}")
check("Jaccard 非零", D["jaccard"][iA, iB] > 0, f"d={D['jaccard'][iA, iB]:.3f}")

# ---- QC 统计兼容 ----
stats = vp.qc_stats_from_af_tables(vcfs, af_min=0.05)
check("QC: S_A 携带数 = 4", stats["S_A"]["n_variants"] == 4,
      str({k: stats['S_A'][k] for k in ('n_variants', 'ts_count', 'tv_count')}))
check("QC: S_A ts/tv 合理 (3 ts, indel 排除)",
      stats["S_A"]["ts_count"] == 3 and stats["S_A"]["tv_count"] == 0, str(stats["S_A"]))

# ---- parse_site_pos 对带下划线 contig 的修复 ----
check("parse_site_pos RefSeq contig", vp.parse_site_pos("NC_002030.1_49_G_A") == 49)
check("parse_site_pos 普通标签", vp.parse_site_pos("NODE_1_100_A_T") == 100)
check("parse_site_pos 兜底 0", vp.parse_site_pos("weird_label") == 0)

print("\n" + ("ALL PASS ✔" if not failures else f"FAILURES: {failures}"))
sys.exit(1 if failures else 0)
