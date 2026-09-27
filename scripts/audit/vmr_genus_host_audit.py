#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""系统性复核: 用 ICTV VMR MSL41 的 "Host source" 列, 逐个核对宿主判定脚本里的
64 属黑名单与 38 科兜底名单, 把 "Plant.tsv 零记录" 这一代理判据换成 ICTV 权威宿主字段。
只读, 不修改管线文件。用法: python3 /tmp/vmr_genus_host_audit.py
"""
import importlib.machinery
import importlib.util
import os
import re
import sys
from collections import Counter, defaultdict

PIPE = "/home/zhangwenda/MMPV-RNA/virome_discovery_pipeline"
CUR = os.path.join(PIPE, "run_host_prediction.py")
VMR = "/home/zhangwenda/plant_virus_db/1.virus-host_db/B-ictv/VMR_MSL41.tsv"
PLANT_PAT = re.compile(r"plant|viridiplant|angiosperm|embryophyta|streptophyta", re.I)


def load_mod(path, name):
    loader = importlib.machinery.SourceFileLoader(name, path)
    spec = importlib.util.spec_from_loader(name, loader)
    m = importlib.util.module_from_spec(spec)
    loader.exec_module(m)
    return m


m = load_mod(CUR, "rhp_vmr")

# VMR 列定位: 用关键字模糊定位 (表头编码可能不标准)
lines = open(VMR, encoding="utf-8", errors="replace").read().splitlines()
hdr = lines[0].split("\t")
print("VMR 表头 %d 列" % len(hdr))


def find_col(*keywords):
    for i, c in enumerate(hdr):
        low = c.lower()
        if all(k in low for k in keywords):
            return i
    return None


c_id = find_col("isolate", "id")
c_fam = find_col("family")
c_gen = find_col("genus")
c_spe = find_col("specie") or None
c_host = find_col("hos")
print("定位: isolate=%s family=%s genus=%s species=%s host=%s"
      % (c_id, c_fam, c_gen, c_spe, c_host))
if c_host is None:
    for i, c in enumerate(hdr):
        print("  col%d=%r" % (i, c))
    sys.exit("找不到宿主列")

host_vals = Counter()
gen_host = defaultdict(Counter)
gen_fam = defaultdict(Counter)
fam_host = defaultdict(Counter)
for line in lines[1:]:
    p = line.split("\t")
    if len(p) <= max(c_fam, c_gen, c_host):
        continue
    fam = p[c_fam].strip()
    gen = p[c_gen].strip()
    host = p[c_host].strip()
    host_vals[host] += 1
    if gen:
        gen_host[gen][host] += 1
        gen_fam[gen][fam] += 1
    if fam:
        fam_host[fam][host] += 1

print()
print("Host source 取值 top20:")
for h, n in host_vals.most_common(20):
    print("   %-40s %d" % (h[:40], n))


def host_has_plant(text):
    return bool(PLANT_PAT.search(text or ""))


def report_plantish(kind, entries, lookup):
    print()
    print("=" * 78)
    print("VMR 宿主列含 plant 的 %s (共 %d 个待核)" % (kind, len(entries)))
    print("=" * 78)
    hits = []
    missing = []
    for e in entries:
        hv = lookup.get(e)
        if not hv:
            missing.append(e)
            continue
        joined = "; ".join(sorted(hv))
        if host_has_plant(joined):
            hits.append((e, hv))
    print("  命中 plant: %d / VMR 无记录: %d" % (len(hits), len(missing)))
    for e, hv in sorted(hits, key=lambda x: -sum(x[1].values())):
        print("   ★ %-22s %s" % (e, ", ".join("%s(%d)" % (k[:40], v) for k, v in hv.most_common(4))))
    if missing:
        print("  VMR 无记录(可能改名/非现行/拼写错):")
        for e in missing:
            print("   ? %s" % e)
    return hits


g_hits = report_plantish("属黑名单", m.NON_PLANT_GENERA, gen_host)
print()
# 同名字面近似 (改名/拼写) 检查: 黑名单条目在 VMR 里查不到, 但存在大小写/近似名
for e in list(m.NON_PLANT_GENERA):
    if e in gen_host:
        continue
    near = [g for g in gen_host if g.lower() == e.lower()]
    if near:
        print("  [名称大小写差异] %s ~ %s" % (e, near))
print()
f_hits = report_plantish("科兜底黑名单(撤回后生效集)", m.NON_PLANT_FAMILIES_FALLBACK, fam_host)

print()
print("=" * 78)
print("白名单 12 科 (ICTV 含植物) 在 VMR 的实况")
print("=" * 78)
for f in sorted(m.ICTV_PLANT_INCLUDING_FAMILIES) if hasattr(m, "ICTV_PLANT_INCLUDING_FAMILIES") else []:
    hv = fam_host.get(f)
    print("  %-20s %s" % (f, ", ".join("%s(%d)" % (k[:30], v) for k, v in hv.most_common(4)) if hv else "VMR 无记录"))
print()
print("DONE")
