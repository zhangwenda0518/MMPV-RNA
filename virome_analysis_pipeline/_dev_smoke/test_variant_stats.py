# -*- coding: utf-8 -*-
"""一次性冒烟测试: variant_species_association.py + virus_haplotype_qst.py
构造合成矩阵与元数据, 端到端运行两个 CLI 并校验输出; 统计内核对拍 scipy。"""
import os
import shutil
import subprocess
import sys

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
UTILS = os.path.join(ROOT, "utils")
TMP = os.path.join(HERE, "_tmp")
shutil.rmtree(TMP, ignore_errors=True)
os.makedirs(TMP)

failures = []


def check(name, cond, detail=""):
    tag = "PASS" if cond else "FAIL"
    print(f"[{tag}] {name}" + (f" | {detail}" if detail else ""))
    if not cond:
        failures.append(name)


# ---------- 合成数据 ----------
# 8 样本: S1-S4 = Lycium barbarum, S5-S8 = Nicotiana amarum
samples = [f"S{i}" for i in range(1, 9)]
rows = {
    "ORF1_100_A>G":   [0.42, 0.51, 0.38, 0.60, 0.0, 0.0, 0.01, 0.02],  # barbarum 特异
    "ORF1_240_C>T":   [0.03, 0.01, 0.02, 0.00, 0.55, 0.48, 0.62, 0.40],  # amarum 特异
    "ORF2_310_G>A":   [0.30, 0.28, 0.35, 0.22, 0.31, 0.29, 0.02, 0.26],  # 两组均携带, 无关联
    "ORF2_500_T>C":   [0.0] * 8,                                        # 固定参考, 应被跳过
}
mat_csv = os.path.join(TMP, "Matrix_05_like.csv")
with open(mat_csv, "w", encoding="utf-8") as f:
    f.write("," + ",".join(samples) + "\n")
    for vid, vals in rows.items():
        f.write(vid + "," + ",".join(str(v) for v in vals) + "\n")

meta_tsv = os.path.join(TMP, "meta.tsv")
with open(meta_tsv, "w", encoding="utf-8") as f:
    f.write("Run\tHost\n")
    for i, s in enumerate(samples):
        host = "Lycium barbarum" if i < 4 else "Nicotiana amarum"
        f.write(f"{s}\t{host}\n")

# ---------- 单元: Fisher / BH ----------
sys.path.insert(0, UTILS)
import variant_species_association as vsa

try:
    from scipy.stats import fisher_exact
    ok = True
    for tab in [([5, 0], [0, 5]), ([3, 1], [1, 5]), ([0, 6], [4, 2]), ([7, 2], [1, 9])]:
        a, b = tab[0][0], tab[0][1]
        c, d = tab[1][0], tab[1][1]
        _, p_ref = fisher_exact([[a, b], [c, d]], alternative="greater")
        p_mine = vsa.fisher_right(a, b, c, d)
        if abs(p_ref - p_mine) > 1e-9:
            ok = False
            print(f"    mismatch {tab}: scipy={p_ref} mine={p_mine}")
    check("fisher_right 对拍 scipy (4 组表)", ok)
except ImportError:
    print("[跳过] 本机无 scipy, Fisher 仅自检边界")
check("fisher 边界 [[5,0],[0,5]]=1/252",
      abs(vsa.fisher_right(5, 0, 0, 5) - 1 / 252) < 1e-9)

q = vsa.bh_adjust([0.01, 0.04, 0.03])
check("BH 对拍 R p.adjust 期望", all(abs(x - y) < 1e-12 for x, y in zip(q, [0.03, 0.04, 0.04])), str(q))

# ---------- 端到端: 关联脚本 ----------
out_assoc = os.path.join(TMP, "assoc_out")
r1 = subprocess.run([sys.executable, os.path.join(UTILS, "variant_species_association.py"),
                     "-m", mat_csv, "-M", meta_tsv, "-o", out_assoc,
                     "--label", "SMOKE", "--alpha", "0.2"],
                    capture_output=True, text=True)
print(r1.stdout[-1200:] if r1.stdout else "", r1.stderr[-400:] if r1.returncode else "", sep="")
check("关联脚本退出码 0", r1.returncode == 0)
full_p = os.path.join(out_assoc, "Species_Variant_Association_Full.tsv")
full_df = pd.read_csv(full_p, sep="\t")
check("受检组合 6 行 (3 可检位点 x 2 物种)", len(full_df) == 6, f"n={len(full_df)}")
sig_p = os.path.join(out_assoc, "SpeciesSpecific_Significant.tsv")
check("Full 表存在", os.path.isfile(full_p))
if os.path.isfile(sig_p):
    import pandas as pd
    sig = pd.read_csv(sig_p, sep="\t")
    hits = set(zip(sig["Site_ID"], sig["Species"]))
    check("barbarum 富集位点检出", ("ORF1_100_A>G", "Lycium barbarum") in hits, str(hits))
    check("amarum 富集位点检出", ("ORF1_240_C>T", "Nicotiana amarum") in hits)
    check("无关联位点未入选", ("ORF2_310_G>A", "Lycium barbarum") not in hits)
    check("热图已生成", os.path.isfile(os.path.join(out_assoc, "species_variant_association.png")))

# ---------- 端到端: QST 脚本 ----------
out_qst = os.path.join(TMP, "qst_out")
r2 = subprocess.run([sys.executable, os.path.join(UTILS, "virus_haplotype_qst.py"),
                     "-m", mat_csv, "-M", meta_tsv, "-o", out_qst,
                     "--label", "SMOKE", "--hamming-threshold", "0.15"],
                    capture_output=True, text=True)
print(r2.stdout[-1500:] if r2.stdout else "", r2.stderr[-500:] if r2.returncode else "", sep="")
check("QST 脚本退出码 0", r2.returncode == 0)
for fn in ["QST_Definitions.tsv", "Sample_QST_Assignment.tsv",
           "QST_by_Host_Species.tsv", "QST_Diversity_Stats.tsv",
           "Private_QST_List.tsv", "Variant_Token_Catalog.tsv"]:
    check(f"QST 输出 {fn}", os.path.isfile(os.path.join(out_qst, fn)))
check("QST 热图 png", os.path.isfile(os.path.join(out_qst, "qst_token_heatmap.png")))
if r2.returncode == 0:
    import pandas as pd
    div = pd.read_csv(os.path.join(out_qst, "QST_Diversity_Stats.tsv"), sep="\t")
    print(div.to_string(index=False))
    defs_n = len(pd.read_csv(os.path.join(out_qst, "QST_Definitions.tsv"), sep="\t"))
    check("3 个 QST 簇", defs_n == 3, f"got {defs_n}")
    am = div[div["Species"] == "Nicotiana amarum"].iloc[0]
    check("amarum Hd≈0.5 且有私有簇缺失时至少检出多簇",
          abs(float(am["Haplotype_Diversity_Hd"]) - 0.5) < 1e-6,
          str(dict(am)))
    ba = div[div["Species"] == "Lycium barbarum"].iloc[0]
    check("barbarum 私有 QST >= 1", int(ba["Private_QST_Count"]) >= 1)
    check("Hd 数值在 [0,1]", all(0 <= float(v) <= 1 for v in div["Haplotype_Diversity_Hd"]))

# ---------- 布局识别: snp_matrix 风格 ----------
mat_tsv = os.path.join(TMP, "snp_matrix_like.tsv")
with open(mat_tsv, "w", encoding="utf-8") as f:
    f.write("Sample\tNODE_1_50_A_G\tNODE_1_90_C_T\n")
    f.write("S1\t1\t0\nS2\t1\tNA\nS3\t0\t1\nS4\tNA\t1\n")
out_qst2 = os.path.join(TMP, "qst_out2")
r3 = subprocess.run([sys.executable, os.path.join(UTILS, "virus_haplotype_qst.py"),
                     "-m", mat_tsv, "-o", out_qst2, "--label", "LAY", "--no-plot"],
                    capture_output=True, text=True)
check("snp_matrix 布局识别 (samples_rows)", r3.returncode == 0 and
      "samples_rows" in r3.stdout, r3.stdout[-200:] + r3.stderr[-200:])

# ---------- 校验报错路径 ----------
bad_meta = os.path.join(TMP, "bad_meta.tsv")
with open(bad_meta, "w", encoding="utf-8") as f:
    f.write("Run\tHost\nS1\tLycium barbarum\nS1\tNicotiana amarum\n")
r4 = subprocess.run([sys.executable, os.path.join(UTILS, "variant_species_association.py"),
                     "-m", mat_csv, "-M", bad_meta, "-o", os.path.join(TMP, "x")],
                    capture_output=True, text=True)
check("重复 Run 报错带行号", r4.returncode != 0 and "bad_meta.tsv:3" in r4.stderr,
      r4.stderr.strip()[-160:])

print("\n" + ("ALL PASS ✔" if not failures else f"FAILURES: {failures}"))
sys.exit(1 if failures else 0)
