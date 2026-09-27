#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""解析 run_host_prediction.py 的非植物科名单, 区分原始 38 与 0915 追加, 并做植物科交叉核对。"""
import ast
import re
import sys

SRC = r"D:\桌面\延伸基因组\MMPV-RNA\virome_discovery_pipeline\run_host_prediction.py"
# ICTV 现行宿主含植物 (或植物侵染证据) 的科: 不得登记为非植物科
PLANT_HOST = {
    "Ourmiaviridae", "Mitoviridae", "Botourmiaviridae", "Partitiviridae", "Endornaviridae",
    "Kanorauviridae", "Metaviridae", "Tombusviridae", "Virgaviridae", "Potyviridae",
    "Betaflexiviridae", "Secoviridae", "Bromoviridae", "Rhabdoviridae", "Caulimoviridae",
    "Alphaflexiviridae", "Tymoviridae", "Closteroviridae", "Solemoviridae", "Tospoviridae",
    "Amalgaviridae", "Fimoviridae", "Phenuiviridae", "Aspiviridae", "Reoviridae",
    "Geminiviridae", "Nanoviridae", "Luteoviridae", "Sobemoviridae", "Virga-like",
}
src = open(SRC, encoding="utf-8").read()
tree = ast.parse(src)
orig, appended = None, []
for node in tree.body:
    if isinstance(node, ast.Assign) and getattr(node.targets[0], "id", "") == "NON_PLANT_FAMILIES_FALLBACK":
        orig = [e.value for e in node.value.elts]
    if isinstance(node, ast.AugAssign) and getattr(node.target, "id", "") == "NON_PLANT_FAMILIES_FALLBACK":
        appended = [e.value for e in node.value.elts]
print("原始科数 %d ; 追加科数 %d ; 合计 %d" % (len(orig or []), len(appended), len(orig or []) + len(appended)))
print("\n原始 38 科全表:\n  %s" % ", ".join(orig or []))
hit_o = sorted(set(orig or []) & PLANT_HOST)
hit_a = sorted(set(appended) & PLANT_HOST)
print("\n原始表中 ICTV 宿主含植物的科: %s" % (hit_o or "无"))
print("追加表中 ICTV 宿主含植物的科: %s" % (hit_a or "无"))
dup = sorted(set(orig or []) & set(appended))
print("重复登记: %s" % (dup or "无"))
print("\n追加名单字母序检查: %s" % ("有序" if appended == sorted(appended) else "非字母序"))
print("追加名单全表 (%d):" % len(appended))
for i in range(0, len(appended), 8):
    print("  " + " ".join(appended[i:i + 8]))
