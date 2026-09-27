#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""独立复算: 用 pandas + 字符串搜索两种不同实现对 Plant.tsv 重算科/属集合,
与代码里的白名单常量逐项比对 (防生成脚本的解析漏洞被固化进代码)。"""
import importlib.util

import pandas as pd

P = "/home/zhangwenda/MMPV-RNA/virome_discovery_pipeline/run_host_prediction.py"
PLANT = "/home/zhangwenda/plant_virus_db/1.virus-host_db/C-host_classify/classified_clean/Plant.tsv"

spec = importlib.util.spec_from_file_location("rhp", P)
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)

df = pd.read_csv(PLANT, sep="\t", dtype=str, low_memory=False)
print("Plant.tsv 行 %d 列 %d: %s" % (len(df), df.shape[1], list(df.columns)[:8]))
seg = df.iloc[:, 3].fillna("").str.split(";")
nseg = seg.apply(len)
print("Virus_lineage 段数分布: %s" % nseg.value_counts().to_dict())

f_vec = seg.apply(lambda p: p[5].strip() if len(p) > 5 else "")
g_vec = seg.apply(lambda p: p[6].strip() if len(p) > 6 else "")
f_set = set(f_vec[f_vec != ""])
g_set = set(g_vec[g_vec != ""])
print("pandas 复算: 科 %d / 属 %d" % (len(f_set), len(g_set)))
print("代码常量:   科白 %d / 属白 %d (含 ICTV 零覆盖 9 科)"
      % (len(m.PLANT_FAMILIES_WHITELIST), len(m.PLANT_GENERA_WHITELIST)))

miss_f = sorted(f_set - set(m.PLANT_FAMILIES_WHITELIST))
miss_g = sorted(g_set - set(m.PLANT_GENERA_WHITELIST))
extra_f = sorted(set(m.PLANT_FAMILIES_WHITELIST) - f_set)
extra_g = sorted(set(m.PLANT_GENERA_WHITELIST) - g_set)
print("Plant.tsv 有而白名单缺: 科 %s 属 %s" % (miss_f or "无", miss_g or "无"))
print("白名单有而 Plant.tsv 无: 科 %s" % (extra_f or "无"))
print("白名单有而 Plant.tsv 无的属: %s" % (extra_g or "无"))

# 第二种实现: 纯文本逐行 index 切分 (不复用 pandas 的 split)
f2, g2 = set(), set()
with open(PLANT, encoding="utf-8", errors="replace") as fh:
    fh.readline()
    for line in fh:
        c = line.rstrip("\n").split("\t")
        if len(c) < 4:
            continue
        p = c[3].split(";")
        if len(p) > 5 and p[5].strip():
            f2.add(p[5].strip())
        if len(p) > 6 and p[6].strip():
            g2.add(p[6].strip())
print("纯文本复算: 科 %d / 属 %d | 与 pandas 结果一致: 科=%s 属=%s"
      % (len(f2), len(g2), f_set == f2, g_set == g2))

ok = not miss_f and not miss_g and f_set == f2 and g_set == g2
print("\n%s" % ("PASS 白名单与 Plant.tsv 独立复算一致" if ok else "FAIL 白名单与复算不一致"))
