#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""诊断: 05 表与 C9 表 merge 后 Family/Genus 是否被 suffix 掉 (决定黑名单是否死代码)。"""
import hashlib
import os
import pandas as pd

A = "/home/zhangwenda/MMPV-paper/goji-virome/02_novel_virus/RNA-Lycium_barbarum_out"
B = "/home/zhangwenda/virus/data-2026/data-test/RNA-Lycium_barbarum_out"
C = "/home/zhangwenda/MMPV-paper/goji-virome/02_novel_virus/RNA-Lycium_barbarum_out"

print("== 树存在性 ==")
for p in (A, B):
    tax = os.path.join(p, "05_Taxonomy/Votus.integrated/final_integrated_classification.tsv")
    print("%s  05表存在=%s" % (p, os.path.isfile(tax)))
    if os.path.isfile(tax):
        with open(tax, "rb") as f:
            print("     md5=%s  size=%d" % (hashlib.md5(f.read()).hexdigest(), os.path.getsize(tax)))

print("\n== merge 测试 ==")
tax_path = os.path.join(A, "05_Taxonomy/Votus.integrated/final_integrated_classification.tsv")
c9_path = os.path.join(A, "06_HostPrediction/C9_ICTV_result/classification_result.tsv")
tax = pd.read_csv(tax_path, sep="\t", nrows=20)
print("tax 列: %s" % list(tax.columns))
_c9_cols = ["contig_id", "Predicted_Host"]
c9_head = pd.read_csv(c9_path, sep="\t", nrows=0)
_c9_cols += [c for c in ("Family", "Genus", "Determination_Level") if c in c9_head.columns]
c9 = pd.read_csv(c9_path, sep="\t", nrows=200)[_c9_cols].rename(
    columns={"Predicted_Host": "Host_ICTV"})
m = tax.merge(c9, on="contig_id", how="left")
print("merge 后列: %s" % list(m.columns))
print("Family in merged? %s ; Genus in merged? %s" % ("Family" in m.columns, "Genus" in m.columns))
if "Family" not in m.columns:
    print(">>> row.get('Family') = %r  → is_blacklisted 拿到空科, 名单失效" % m.iloc[0].get("Family"))
else:
    print(">>> row.get('Family') = %r" % m.iloc[0].get("Family"))

print("\n== 现网 ensemble_host_summary 头 ==")
s = os.path.join(A, "06_HostPrediction/ensemble_host_summary.tsv")
h = pd.read_csv(s, sep="\t", nrows=0)
print("%s" % list(h.columns))
print("Family_x/Family_y 存在? %s / %s" % ("Family_x" in h.columns, "Family_y" in h.columns))
