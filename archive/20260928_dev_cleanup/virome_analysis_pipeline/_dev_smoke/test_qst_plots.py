# -*- coding: utf-8 -*-
"""QST 绘图增强冒烟测试: 分组bar/pie、性状箱线、网络图、动态字号。"""
import os
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
UTILS = os.path.join(ROOT, "utils")
TMP = os.path.join(HERE, "_tmp_plots")
shutil.rmtree(TMP, ignore_errors=True)
os.makedirs(TMP)

failures = []


def check(name, cond, detail=""):
    print(f"[{'PASS' if cond else 'FAIL'}] {name}" + (f" | {detail}" if detail else ""))
    if not cond:
        failures.append(name)


# ---- 合成数据: 12 样本, 6 位点, 三种物种 + Location + 数值 trait ----
samples = [f"S{i}" for i in range(1, 13)]
rows = {
    "ORF1_100_A_G": [0.42, 0.51, 0.38, 0.60, 0.0, 0.0, 0.01, 0.02, 0.0, 0.0, 0.0, 0.02],
    "ORF1_240_C_T": [0.03, 0.01, 0.02, 0.00, 0.55, 0.48, 0.62, 0.40, 0.10, 0.12, 0.08, 0.09],
    "ORF2_310_G_A": [0.30, 0.28, 0.35, 0.22, 0.31, 0.29, 0.02, 0.26, 0.90, 0.85, 0.88, 0.95],
    "ORF2_500_T_C": [0.0] * 12,
    "ORF3_600_G_T": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.70, 0.65, 0.60, 0.55, 0.50, 0.45],
    "ORF3_700_A_C": [0.2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
}
mat_csv = os.path.join(TMP, "mat.csv")
with open(mat_csv, "w") as f:
    f.write("," + ",".join(samples) + "\n")
    for vid, vals in rows.items():
        f.write(vid + "," + ",".join(str(v) for v in vals) + "\n")

meta_tsv = os.path.join(TMP, "meta.tsv")
with open(meta_tsv, "w") as f:
    f.write("Run\tHost\tLocation\tDepth\n")
    for i, s in enumerate(samples):
        sp = ["Lycium barbarum", "Lycium ruthenicum", "Nicotiana amarum"][i % 3]
        loc = ["China", "Tibet", "Italy"][i % 3]
        depth = 100 + i * 7 + (50 if i % 2 else 0)
        f.write(f"{s}\t{sp}\t{loc}\t{depth}\n")

out = os.path.join(TMP, "out")
r = subprocess.run(
    [sys.executable, os.path.join(UTILS, "virus_haplotype_qst.py"),
     "-m", mat_csv, "-M", meta_tsv, "-o", out, "--label", "PLOTTEST",
     "--trait-cols", "Depth", "--network-min-count", "1"],
    capture_output=True, text=True)
print(r.stdout[-1400:] if r.stdout else "", r.stderr[-600:] if r.returncode else "", sep="")
check("QST 绘图脚本退出码 0", r.returncode == 0)

for fn in ["qst_token_heatmap.png", "qst_group_stacked.png", "qst_group_pies.png",
           "qst_haplotype_network.png", "qst_Depth_trait_boxplot.png",
           "QST_Network_Edges.tsv"]:
    check(f"输出 {fn}", os.path.isfile(os.path.join(out, fn)))

# 网络图边表内容
if os.path.isfile(os.path.join(out, "QST_Network_Edges.tsv")):
    with open(os.path.join(out, "QST_Network_Edges.tsv")) as f:
        lines = f.read().strip().split("\n")
    check("网络图边表有 MST 边", len(lines) >= 2, f"{len(lines) - 1} 条边")

# 分组 bar/pie 存在即验证绘制成功 (进程退出码 0 且文件存在)

# ---- --group-col Location 覆盖分组 ----
out2 = os.path.join(TMP, "out2")
r2 = subprocess.run(
    [sys.executable, os.path.join(UTILS, "virus_haplotype_qst.py"),
     "-m", mat_csv, "-M", meta_tsv, "-o", out2, "--label", "LOC",
     "--group-col", "Location", "--no-plot"],
    capture_output=True, text=True)
check("--group-col 分组退出码 0", r2.returncode == 0)
if r2.returncode == 0:
    with open(os.path.join(out2, "QST_by_Host_Species.tsv")) as f:
        header = f.readline()
    check("Location 分组生效", "China" in header and "Tibet" in header, header.strip())

# ---- 无元数据 (物种层跳过, 网络图节点全 NA) ----
out3 = os.path.join(TMP, "out3")
r3 = subprocess.run(
    [sys.executable, os.path.join(UTILS, "virus_haplotype_qst.py"),
     "-m", mat_csv, "-o", out3, "--label", "NOMETA"],
    capture_output=True, text=True)
check("无元数据退出码 0", r3.returncode == 0)
check("无元数据网络图仍出", os.path.isfile(os.path.join(out3, "qst_haplotype_network.png")))

print("\n" + ("ALL PASS ✔" if not failures else f"FAILURES: {failures}"))
sys.exit(1 if failures else 0)
