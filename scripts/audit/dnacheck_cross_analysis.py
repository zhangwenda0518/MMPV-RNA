#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""eve_screen 扫描结果 × 四象限位点 × 病毒组 DNA contig 三向交叉。
坐标归一：ningxia.clean.fa 的 nSX%05d → chr%02d（+1）。"""
import csv, glob, re
from collections import defaultdict

SCAN = "/home/zhangwenda/eve_screen_goji3_20260929/04_Summary"
QUAD = "/home/zhangwenda/biosoft/virus/hi-fever/eve_x_te/caulifinder_overlap"
VERD = "/home/zhangwenda/goji_dnacheck_20260929/final_dna_verdict.tsv"

SPECIES = ["zhonghua", "ningxia", "heiguo"]


def norm_chr(c, sp):
    c = c.strip()
    if sp == "ningxia":
        m = re.match(r"nSX(\d+)", c)
        if m:
            return "chr%02d" % (int(m.group(1)) + 1)
    return c


def parse_locus(loc, sp):
    m = re.match(r"^(.+):(\d+)-(\d+)$", loc)
    if not m:
        return None
    return (norm_chr(m.group(1), sp), int(m.group(2)), int(m.group(3)))


def overlaps(a, b):
    return a[0] == b[0] and a[1] <= b[2] and b[1] <= a[2]


report = []
for sp in SPECIES:
    # 1) 扫描 viral_supported 位点
    scan = []
    with open(f"{SCAN}/{sp}_eve_summary.tsv") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r["verdict"] == "viral_supported":
                p = parse_locus(r["locus"], sp)
                if p:
                    scan.append((p, r))
    # 2) 四象限 v2 位点
    quad = []
    with open(f"{QUAD}/{sp}_v2.bed") as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 4:
                quad.append((norm_chr(p[0], sp), int(p[1]), int(p[2]), p[3]))
    # 3) 病毒组 EVE_strong contig 的基因组命中
    vhits = []
    with open(VERD) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r["genome_evidence"] == "EVE_strong" and r["Cohort"].lower().startswith(
                {"zhonghua": "lycium_chinense", "ningxia": "lycium_barbarum", "heiguo": "lycium_ruthenicum"}[sp]):
                m = re.search(r"(Chr\w+|chr\w+|nSX\d+):(\d+)-(\d+)", r["evidence_detail"])
                if m:
                    vhits.append((norm_chr(m.group(1), sp), int(m.group(2)), int(m.group(3)), r["Species"]))

    # 交叉：四象限位点被扫描覆盖的比例
    q_hit = sum(1 for q in quad if any(overlaps((q[0], q[1], q[2]), s[0]) for s in scan))
    # 扫描位点中被四象限覆盖的比例
    s_hit = sum(1 for s in scan if any(overlaps(s[0], (q[0], q[1], q[2])) for q in quad))
    # 病毒组命中位点落进扫描 viral_supported 的数量
    v_in_scan = sum(1 for v in vhits if any(overlaps((v[0], v[1], v[2]), s[0]) for s in scan))
    report.append((sp, len(scan), len(quad), q_hit, len(scan), s_hit, len(vhits), v_in_scan))
    print(f"[{sp}] scan_viral_supported={len(scan)}  四象限位点={len(quad)}")
    print(f"  四象限位点与扫描重叠: {q_hit}/{len(quad)}")
    print(f"  扫描位点与四象限重叠: {s_hit}/{len(scan)}")
    print(f"  病毒组 EVE_strong 命中: {len(vhits)} 条, 落入扫描位点: {v_in_scan}")
    print()

# 病毒组 392 候选 vs 扫描位点（蛋白质级敏感扫描是否找到候选的基因组对应区）
print("=== 392 候选 contig 的物种队列 vs 各物种扫描位点（坐标级） ===")
# 注意：候选无基因组命中，无法定位坐标——改为按"队列宿主 vs 同物种扫描位点数"给背景
