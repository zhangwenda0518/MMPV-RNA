#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""打印属黑名单全表 + 核查特定科/属是否在名单中 (定位 Aspiviridae/Miraophiovirus 这类否决来源)。"""
import importlib.util

P = "/home/zhangwenda/MMPV-RNA/virome_discovery_pipeline/run_host_prediction.py"
spec = importlib.util.spec_from_file_location("rhp", P)
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)

G = sorted(m.NON_PLANT_GENERA)
print("属黑名单 %d 条:" % len(G))
for i in range(0, len(G), 6):
    print("  " + " ".join("%-24s" % g for g in G[i:i + 6]))
print("\n指定属在名单中?")
for g in ["Miraophiovirus", "Ophiovirus", "Biavirus", "Rimosavirus", "Tethysvirus",
          "Alphambiguivirus", "Betambiguivirus", "Ourmiavirus", "Potyvirus", "Potexvirus",
          "Fadolivirus", "Cotonvirus", "Metavirus", "Betanucleorhabdovirus", "Iflavirus"]:
    print("   %-24s %s" % (g, g in m.NON_PLANT_GENERA))
print("\n指定科在名单中?")
for f in ["Aspiviridae", "Ophioviridae", "Tombusviridae", "Partitiviridae", "Ambiguiviridae",
          "Astroviridae", "Kyanoviridae", "Schizomimiviridae", "Nimaviridae", "Papillomaviridae"]:
    print("   %-20s %s" % (f, f in m.NON_PLANT_FAMILIES_FALLBACK))
