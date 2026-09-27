#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""回归验证: 加白名单前(bak) vs 加白名单后(现装) 是否 *逐行* 完全等价。

做法 (不依赖任何聚合口径, 直接比 decision_tree_cascade 的逐行输出):
  1) 两个模块分别用 importlib 加载, 互不干扰;
  2) 两棵树的 05 表各构造一次 merged frame (同一份输入);
  3) 分别用两版 decision_tree_cascade 判完所有行, 逐行比对;
  4) 另外把两棵树里出现过的全部 (Family, Genus, Determination_Level) 组合
     喂给两版 is_blacklisted, 逐三元组比对 (确保否决函数也没变)。
任何一处不一致都会打印出来, 全等则打印 PASS。
"""
import importlib.util
import os
import sys

import pandas as pd

PIPE = "/home/zhangwenda/MMPV-RNA/virome_discovery_pipeline"
CUR = os.path.join(PIPE, "run_host_prediction.py")
BAK = os.path.join(PIPE, "run_host_prediction.py.bak_plantwl_20260915")

TREES = {
    "goji": "/home/zhangwenda/MMPV-paper/goji-virome/02_novel_virus/RNA-Lycium_barbarum_out",
    "onekp": "/home/zhangwenda/MMPV-paper/onekp-virome/onekp-virus",
}


def load_mod(path, name):
    # 备份文件后缀不是 .py, spec_from_file_location 推不出 loader, 显式给 SourceFileLoader
    loader = importlib.machinery.SourceFileLoader(name, path)
    spec = importlib.util.spec_from_loader(name, loader)
    m = importlib.util.module_from_spec(spec)
    loader.exec_module(m)
    return m


def load_frame(tree):
    tax = os.path.join(tree, "05_Taxonomy/Votus.integrated/final_integrated_classification.tsv")
    out = os.path.join(tree, "06_HostPrediction")
    df = pd.read_csv(tax, sep="\t")
    c9 = pd.read_csv(os.path.join(out, "C9_ICTV_result/classification_result.tsv"), sep="\t")
    c9 = c9[[c for c in ["contig_id", "Predicted_Host", "Determination_Level"] if c in c9.columns]]
    df = df.merge(c9.rename(columns={"Predicted_Host": "Host_ICTV"}), on="contig_id", how="left")
    rvp = os.path.join(out, "RVH_result/result.csv")
    if os.path.isfile(rvp):
        rvh = pd.read_csv(rvp)
        idc = rvh.columns[0]
        if "Unnamed" in idc or "y|virus order" in rvh.columns:
            rvh = rvh.rename(columns={idc: "contig_id"})
        df = df.merge(rvh[[c for c in ["contig_id", "pred|L1", "pred|L2", "evidence"]
                           if c in rvh.columns]], on="contig_id", how="left")
    pbp = os.path.join(out, "phabox2_output/final_prediction/cherry_prediction.tsv")
    if os.path.isfile(pbp):
        pb = pd.read_csv(pbp, sep="\t").rename(columns={"Accession": "contig_id"})
        df = df.merge(pb[[c for c in ["contig_id", "Host", "Host_NCBI_lineage"]
                          if c in pb.columns]], on="contig_id", how="left")
    return df


cur = load_mod(CUR, "rhp_cur")
bak = load_mod(BAK, "rhp_bak")
print("现装 md5 侧常量: 科黑 %d 属黑 %d 科白 %d 属白 %d"
      % (len(cur.NON_PLANT_FAMILIES_FALLBACK), len(cur.NON_PLANT_GENERA),
         len(cur.PLANT_FAMILIES_WHITELIST), len(cur.PLANT_GENERA_WHITELIST)))
print("备份侧常量:     科黑 %d 属黑 %d"
      % (len(bak.NON_PLANT_FAMILIES_FALLBACK), len(bak.NON_PLANT_GENERA)))
same_f = set(cur.NON_PLANT_FAMILIES_FALLBACK) == set(bak.NON_PLANT_FAMILIES_FALLBACK)
same_g = set(cur.NON_PLANT_GENERA) == set(bak.NON_PLANT_GENERA)
print("黑名单是否完全相同: 科=%s 属=%s" % (same_f, same_g))
assert same_f and same_g, "黑名单被动过, 必须人工检查"

bad = 0
for tag, tree in TREES.items():
    df = load_frame(tree)
    h_cur = [r[0] for r in df.apply(cur.decision_tree_cascade, axis=1)]
    h_bak = [r[0] for r in df.apply(bak.decision_tree_cascade, axis=1)]
    diff = [i for i in range(len(df)) if h_cur[i] != h_bak[i]]
    print("\n[%s] 行 %d | 现装 Plant=%d | 备份 Plant=%d | 逐行差异 %d"
          % (tag, len(df), h_cur.count("Plant"), h_bak.count("Plant"), len(diff)))
    for i in diff[:10]:
        r = df.iloc[i]
        print("   差异 %s F=%s G=%s det=%s 现装=%s 备份=%s"
              % (r.get("contig_id"), r.get("Family"), r.get("Genus"),
                 r.get("Determination_Level"), h_cur[i], h_bak[i]))
    bad += len(diff)

    triples = set()
    for c, g, d in zip(df.get("Family"), df.get("Genus"), df.get("Determination_Level")):
        triples.add((c, g, d))
    d2 = [t for t in triples if bool(cur.is_blacklisted(*t)) != bool(bak.is_blacklisted(*t))]
    print("  is_blacklisted 三元组 %d 个, 两版判定不一致 %d 个" % (len(triples), len(d2)))
    for t in d2[:10]:
        print("   三元组差异 %s 现装=%s 备份=%s" % (t, cur.is_blacklisted(*t), bak.is_blacklisted(*t)))
    bad += len(d2)

print("\n总差异: %d  %s" % (bad, "PASS 逐行等价" if bad == 0 else "FAIL 存在不一致"))
sys.exit(0 if bad == 0 else 1)
