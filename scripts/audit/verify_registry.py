#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""注册表一致性自检: 导入 run_host_prediction.py 并断言白名单/黑名单关系。"""
import importlib.util
import sys

P = "/home/zhangwenda/MMPV-RNA/virome_discovery_pipeline/run_host_prediction.py"
spec = importlib.util.spec_from_file_location("rhp", P)
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)

fams, gens = set(m.NON_PLANT_FAMILIES_FALLBACK), set(m.NON_PLANT_GENERA)
wf, wg = set(m.PLANT_FAMILIES_WHITELIST), set(m.PLANT_GENERA_WHITELIST)
ictv = set(m.ICTV_PLANT_INCLUDING_FAMILIES)
print("导入 OK | 科黑 %d / 属黑 %d / 科白 %d / 属白 %d / ICTV 含植物科 %d"
      % (len(fams), len(gens), len(wf), len(wg), len(ictv)))
print("FAMILY_FIRST_VETO=%s | WHITELIST_OVERRIDES_FAMILY_VETO=%s"
      % (m.FAMILY_FIRST_VETO, m.WHITELIST_OVERRIDES_FAMILY_VETO))
bad_f, bad_g = sorted(wf & fams), sorted(wg & gens)
print("白黑交集: 科=%s 属=%s" % (bad_f or "空", bad_g or "空"))
print("ICTV 含植物 12 科是否被拉黑: %s" % sorted(ictv & fams))
missing = sorted(ictv - wf)
print("ICTV 含植物科未进白名单的: %s" % (missing or "无"))
assert not bad_f and not bad_g, "白黑交集非空"
assert not (ictv & fams), "ICTV 含植物科被拉黑"
assert not missing, "ICTV 含植物科未登记白名单"
print("全部断言通过")
