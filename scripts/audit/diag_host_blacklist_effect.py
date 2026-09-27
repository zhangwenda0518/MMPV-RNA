#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""端到端验证 run_host_prediction.run_ensemble 的字段碰撞对黑名单的影响。

结论口径:
  A. 现网代码 (Family/Genus 被 suffix 成 _x/_y) → is_blacklisted 拿到 None
  B. 修正后 (把 _x 还原成 Family/Genus)      → is_blacklisted 能真正生效
  比较两版 Final_Host 分布差异, 并列出翻转行。
"""
import importlib.util
import os
from collections import Counter

import pandas as pd

RHP = "/home/zhangwenda/MMPV-RNA/virome_discovery_pipeline/run_host_prediction.py"
spec = importlib.util.spec_from_file_location("rhp", RHP)
rhp = importlib.util.module_from_spec(spec)
spec.loader.exec_module(rhp)

TREE = "/home/zhangwenda/MMPV-paper/goji-virome/02_novel_virus/RNA-Lycium_barbarum_out"
TAX = os.path.join(TREE, "05_Taxonomy/Votus.integrated/final_integrated_classification.tsv")
OUT = os.path.join(TREE, "06_HostPrediction")
C9 = os.path.join(OUT, "C9_ICTV_result/classification_result.tsv")
RVH = os.path.join(OUT, "RVH_result/result.csv")
PB2 = os.path.join(OUT, "phabox2_output/final_prediction/cherry_prediction.tsv")

print("名单规模: 科兜底 %d / 属黑名单 %d" % (
    len(rhp.NON_PLANT_FAMILIES_FALLBACK), len(rhp.NON_PLANT_GENERA)))

# ── 复刻 run_ensemble 的合并逻辑 ──
df = pd.read_csv(TAX, sep="\t")
print("tax 列数 %d; 含 Family=%s Genus=%s" % (len(df.columns), "Family" in df.columns,
                                            "Genus" in df.columns))
_c9_cols = ["contig_id", "Predicted_Host"]
_c9_extra = [c for c in ("Family", "Genus", "Determination_Level")
             if c in pd.read_csv(C9, sep="\t", nrows=0).columns]
_c9_cols += _c9_extra
c9 = pd.read_csv(C9, sep="\t")[_c9_cols].rename(columns={"Predicted_Host": "Host_ICTV"})
merged = df.merge(c9, on="contig_id", how="left")
if os.path.isfile(RVH):
    rvh = pd.read_csv(RVH)
    id_col = rvh.columns[0]
    if "Unnamed" in id_col or "y|virus order" in rvh.columns:
        rvh = rvh.rename(columns={id_col: "contig_id"})
    use = ["contig_id"] + [c for c in ["pred|L1", "pred|L2", "evidence"] if c in rvh.columns]
    merged = merged.merge(rvh[use], on="contig_id", how="left")
if os.path.isfile(PB2):
    pb2 = pd.read_csv(PB2, sep="\t").rename(columns={"Accession": "contig_id"})
    use = [c for c in ["contig_id", "Host", "Host_NCBI_lineage", "Host_GTDB_lineage"]
           if c in pb2.columns]
    merged = merged.merge(pb2[use], on="contig_id", how="left")

print("\n[A] 现网字段视角: merged 有 'Family' 列 = %s → row.get('Family') = %r"
      % ("Family" in merged.columns, merged.iloc[0].get("Family")))
res_a = merged.apply(rhp.decision_tree_cascade, axis=1)
ca = Counter(h for h, _ in res_a)
print("    现网 Final_Host 分布: %s" % dict(ca))
n_bl_a = sum(1 for _, row in merged.iterrows()
             if rhp.is_blacklisted(row.get("Family"), row.get("Genus"),
                                   row.get("Determination_Level")))
print("    现网 is_blacklisted 命中行数: %d" % n_bl_a)

# ── 修正视角: 用 C9 侧的 _y 列 (原意) 与 tax 侧 _x 列分别算 ──
for tag, fam_col, gen_col in [("B1 用 tax 侧 (.x)", "Family_x", "Genus_x"),
                              ("B2 用 C9 侧 (.y)", "Family_y", "Genus_y")]:
    n_bl = sum(1 for _, row in merged.iterrows()
               if rhp.is_blacklisted(row.get(fam_col), row.get(gen_col),
                                     row.get("Determination_Level")))
    print("    %s is_blacklisted 命中行数: %d" % (tag, n_bl))

# ── C/F/G 一致性: C9 与 05 表是否同源 ──
if "Family_x" in merged.columns and "Family_y" in merged.columns:
    fx = merged["Family_x"].fillna("NA").astype(str).str.strip()
    fy = merged["Family_y"].fillna("NA").astype(str).str.strip()
    gx = merged["Genus_x"].fillna("NA").astype(str).str.strip()
    gy = merged["Genus_y"].fillna("NA").astype(str).str.strip()
    print("\n[C] Family 不一致行数: %d ; Genus 不一致行数: %d (共 %d)"
          % (int((fx != fy).sum()), int((gx != gy).sum()), len(merged)))
    a = merged.assign(r=res_a).groupby(merged["Family_x"].fillna("NA").astype(str), dropna=False)
    print("以上不一致样例: %s" % list(zip(fx[fx != fy][:5], fy[fx != fy][:5])))

# ── 用修正字段重跑决策树, 看翻转 ──
merged_b = merged.rename(columns={"Family_x": "Family", "Genus_x": "Genus"})
res_b = merged_b.apply(rhp.decision_tree_cascade, axis=1)
cb = Counter(h for h, _ in res_b)
print("\n[B] 修正后 Final_Host 分布: %s" % dict(cb))
diff = [(merged_b.iloc[i]["contig_id"], res_a.iloc[i][0], res_b.iloc[i][0],
         merged.iloc[i].get("Family_x"), merged.iloc[i].get("Genus_x"))
        for i in range(len(merged)) if res_a.iloc[i][0] != res_b.iloc[i][0]]
print("    翻转行数: %d" % len(diff))
for d in diff[:25]:
    print("      %s | %s -> %s | %s/%s" % d)

# ── 落盘对照 ──
disk = pd.read_csv(os.path.join(OUT, "ensemble_host_summary.tsv"), sep="\t",
                   usecols=["contig_id", "Final_Host"])
m2 = merged[["contig_id"]].merge(disk, on="contig_id", how="left")
print("\n[D] 落盘 Final_Host 分布: %s"
      % dict(Counter(m2["Final_Host"].fillna("MISSING"))))
print("    现网代码 vs 落盘 不一致行数: %d"
      % int((res_a.reset_index(drop=True) != m2["Final_Host"].reset_index(drop=True)).sum()))
