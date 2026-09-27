#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""用修复后的 _write_rescue_report 重新生成报告, 不重跑流程"""
import sys, os
from pathlib import Path
from Bio import SeqIO

sys.path.insert(0, "/home/zhangwenda/MMPV-RNA/virome_discovery_pipeline")
import rescue_pipeline as rp

D = Path("/home/zhangwenda/data-test/out/08_Rescue/Plant")

# 复现 main() 里的路径推导
centroids_fa = os.environ.get("CENTROIDS")
if not centroids_fa:
    cand0 = "/home/zhangwenda/data-test/out/04_CLUSTER/rescue_centroids.fasta"
    centroids_fa = cand0 if os.path.isfile(cand0) else None
print("centroids:", centroids_fa)

centroids_records = list(SeqIO.parse(centroids_fa, "fasta"))
print("centroids 条数:", len(centroids_records))

# 各分支产物
fa_a = D / "branch_a" / "branchA_pass.fasta"
fa_b = D / "branch_b" / "branchB_pass.fasta"
fa_c = D / "branch_c" / "branchC_pass.fasta"
fa_d = D / "branch_d" / "branchD_pass.fasta"

def cnt(p):
    return sum(1 for _ in SeqIO.parse(str(p), "fasta")) if p.is_file() else 0

cnt_a, cnt_b, cnt_c, cnt_d = cnt(fa_a), cnt(fa_b), cnt(fa_c), cnt(fa_d)
print(f"分支产物条数: A={cnt_a} B={cnt_b} C={cnt_c} D={cnt_d}")

final_fa = D / "centroids" / "final_centroids.fasta"
n_final = sum(1 for _ in SeqIO.parse(str(final_fa), "fasta")) if final_fa.is_file() else 0
print("最终 vOTU:", n_final)

# 找 taxonomy / genus_len (从原始启动命令中提取的真实路径)
tax_tsv = "/home/zhangwenda/data-test/out/05_Taxonomy/Votus.integrated/final_integrated_classification.tsv"
if not os.path.isfile(tax_tsv):
    tax_tsv = None

genus_len = "/home/zhangwenda/database/virus-db/db/genus_lens"
if not os.path.isfile(genus_len):
    genus_len = None

# 分支 A: prepass 计数 (不在 centroids 中, 免拯救)
pre_pass_file = D / "branch_a_prepass.txt"
n_pre = sum(1 for _ in open(pre_pass_file)) if pre_pass_file.is_file() else 0
cnt_a += n_pre
print(f"分支 A prepass 追加: +{n_pre}")

print("taxonomy:", tax_tsv)
print("genus_len:", genus_len)
print("")

clusters = {}
try:
    ctsv = "/home/zhangwenda/data-test/out/04_CLUSTER/vclust_clusters.tsv"
    if os.path.isfile(ctsv):
        clusters = rp.parse_vclust_clusters(ctsv)
except Exception as e:
    print("clusters 加载跳过:", e)

rp._write_rescue_report(
    D, centroids_records, clusters,
    str(fa_a) if fa_a.is_file() else None,
    str(fa_b) if fa_b.is_file() else None,
    str(fa_c) if fa_c.is_file() else None,
    str(fa_d) if fa_d.is_file() else None,
    cnt_a, cnt_b, cnt_c, cnt_d, n_final,
    taxonomy_tsv=tax_tsv, genus_lens_path=genus_len, min_vsi_len=2000)
print("\n报告已重新生成")
