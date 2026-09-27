#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""量各工具在「属->科」「属->域」两个可参照判据上的实测一致率，用来给平票兜底排序。

判据来源：genus_family_ref.tsv（Genus, NCBI_Family, NCBI_n, VMR_Family, VMR_n, Domain）
  A. 属可查时，工具报的 Family 是否等于参照科（NCBI 或 VMR 任一命中即算对）
  B. 属可查时，工具报的 Realm 是否等于参照 Domain
  C. 属可查时，工具报的 Genus 是否属于参照库已知属
用法: python3 measure_tool_consistency.py <combined.tsv> <ref.tsv>
"""
import csv
import sys
from collections import defaultdict

LV = ["Realm", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
EMPTY = {"", "na", "n/a", "nan", "-", "none", "null", "unclassified", "unknown"}


def nz(x):
    x = (x or "").strip()
    return "" if x.lower() in EMPTY else x


combined, ref_path = sys.argv[1], sys.argv[2]

# ── 参照表 ──
ref = {}
with open(ref_path, newline="", encoding="utf-8", errors="replace") as f:
    for r in csv.DictReader(f, delimiter="\t"):
        g = nz(r.get("Genus"))
        if not g:
            continue
        fams = {x.lower() for x in (nz(r.get("NCBI_Family")), nz(r.get("VMR_Family"))) if x}
        ref[g.lower()] = {"fams": fams, "domain": nz(r.get("Domain")), "name": g}

# ── 工具统计 ──
stat = defaultdict(lambda: {"fam_known": 0, "fam_hit": 0, "realm_known": 0, "realm_hit": 0,
                            "genus_known": 0, "row": 0, "genus_calls": 0, "fam_calls": 0,
                            "species_calls": 0, "sp_genus_ok": 0, "sp_genus_n": 0})

with open(combined, newline="", encoding="utf-8", errors="replace") as f:
    for r in csv.DictReader(f, delimiter="\t"):
        tool = (r.get("tool") or "").strip()
        if not tool:
            continue
        s = stat[tool]
        s["row"] += 1
        g = nz(r.get("Genus"))
        fam = nz(r.get("Family"))
        realm = nz(r.get("Realm"))
        sp = nz(r.get("Species"))
        if fam:
            s["fam_calls"] += 1
        if g:
            s["genus_calls"] += 1
        if sp:
            s["species_calls"] += 1
        rg = ref.get(g.lower()) if g else None
        if rg:
            s["genus_known"] += 1
            if fam:
                s["fam_known"] += 1
                if rg["fams"] and fam.lower() in rg["fams"]:
                    s["fam_hit"] += 1
            if realm and rg["domain"]:
                s["realm_known"] += 1
                if realm.lower() == rg["domain"].lower():
                    s["realm_hit"] += 1
        # 种名首词是否等于所报属（内部一致性）
        if sp and g:
            s["sp_genus_n"] += 1
            if sp.split()[0].lower() == g.lower():
                s["sp_genus_ok"] += 1


def rate(a, b):
    return (100.0 * a / b) if b else float("nan")


rows = []
for tool, s in stat.items():
    fam_r = rate(s["fam_hit"], s["fam_known"])
    realm_r = rate(s["realm_hit"], s["realm_known"])
    sp_r = rate(s["sp_genus_ok"], s["sp_genus_n"])
    rows.append((tool, s, fam_r, realm_r, sp_r))

rows.sort(key=lambda x: (-(x[2] if x[2] == x[2] else -1), -x[1]["fam_known"]))

print("参照表属数 %d" % len(ref))
print()
print("| 工具 | 行数 | 报科数 | 属可查 | 科一致 | 科一致率 | 域可查 | 域一致 | 域一致率 | 种首词=属 | 一致率 |")
print("|---|---|---|---|---|---|---|---|---|---|---|")
for tool, s, fam_r, realm_r, sp_r in rows:
    print("| %s | %d | %d | %d | %d | %.1f%% | %d | %d | %.1f%% | %d/%d | %.1f%% |" % (
        tool, s["row"], s["fam_calls"], s["fam_known"], s["fam_hit"], fam_r,
        s["realm_known"], s["realm_hit"], realm_r, s["sp_genus_ok"], s["sp_genus_n"], sp_r))

print()
print("按「科一致率」降序：")
print(",".join(r[0] for r in rows))
