#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
probe_ref_ws_and_genus.py —— v6.4 逻辑层复核的两个真实数据探针

探针 1（归一化不对称的规模）：参照表 genus_family_ref.tsv 里 NCBI_Family / VMR_Family
        是否含「内部空白」。闸门对行侧 Family（fam_c = gsub("[[:space:]]+","",...)）和
        VMR 侧（vmr_pad 同样去空白）都去了内部空白，唯独 NCBI 侧只做 tolower(tax_norm())
        保留内部空白 —— 若参照表 NCBI_Family 真含空白，两侧容忍度不一致会造成误判置空。

探针 2（harmonize 顺序冲突的规模）：产物 final_integrated_classification.tsv 中
        Species 为双名法（首词以 virus 结尾）时，其首词构成一条「独立属主张」；
        若最终 Genus 与该首词不一致，说明 harmonize_genus_species（由种提属）先跑、
        harmonize_family_genus（科-属校准）后跑把 Genus 改回去了，而两道闸门都不复核
        这层关系；同时看这类行的 Species_agree 是否仍写着 1/1（表面一致的误导信号）。

用法：python3 probe_ref_ws_and_genus.py
"""
import csv
import os
import re
import sys

REF = os.path.expanduser("~/database/taxonomy/genus_family_ref.tsv")

TABLES = [
    ("live_v62产物   ", os.path.expanduser(
        "~/MMPV-paper/goji-virome/02_novel_virus/RNA-Lycium_barbarum_out/"
        "05_Taxonomy/Votus.integrated/final_integrated_classification.tsv")),
    ("v64 legacy     ", "/tmp/rc64b_legacy/final_integrated_classification.tsv"),
    ("v64 cascade    ", "/tmp/rc64b_cascade/final_integrated_classification.tsv"),
]

BAD_SUFFIX = re.compile(r"viridae$|virinae$", re.I)
SP_PREFIX = re.compile(r"^unclassified\b|^uncultured\b|^environmental\b", re.I)


def read_tsv(path):
    with open(path, "r", encoding="utf-8", newline="") as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        return rdr.fieldnames, list(rdr)


def valid(v):
    if v is None:
        return False
    s = v.strip().strip('"')
    return s != "" and s.upper() not in ("NA", "N/A", "-")


def probe_ref():
    print("=" * 78)
    print("探针 1: 参照表内部空白（归一化不对称的规模）")
    print("=" * 78)
    if not os.path.exists(REF):
        print("  [跳过] 参照表不存在: %s" % REF)
        return
    cols, rows = read_tsv(REF)
    print("  文件: %s" % REF)
    print("  行数: %d ; 列: %s" % (len(rows), ",".join(cols or [])))
    stat = {}
    for col in ("NCBI_Family", "VMR_Family"):
        n_inner = 0
        ex = []
        for r in rows:
            v = (r.get(col) or "").strip().strip('"')
            if re.search(r"\s", v):          # trim 后仍有空白 = 内部空白
                n_inner += 1
                if len(ex) < 5:
                    ex.append("%s=[%s]" % (r.get("Genus"), v))
        stat[col] = n_inner
        print("  %s 含内部空白: %d 行 %s" % (col, n_inner, ("例: " + " ; ".join(ex)) if ex else ""))
    # 一属多科（VMR 侧 ; 分隔）规模，用于说明 vmr_pad 去空白的适用面
    multi = sum(1 for r in rows if ";" in (r.get("VMR_Family") or ""))
    print("  VMR_Family 含 ';'（一属多科）: %d 行" % multi)
    # 行侧真正会被误判的形态：NCBI 含空白且去掉空白后可见于产品 Family 取值
    risk = []
    for r in rows:
        nv = (r.get("NCBI_Family") or "").strip().strip('"')
        if re.search(r"\s", nv):
            risk.append((r.get("Genus"), nv, nv.replace(" ", "")))
    if risk:
        print("  受影响属（前 10）: Genus | NCBI_Family(带空白) | 去空白后形态")
        for g, a, b in risk[:10]:
            print("    %s | %s | %s" % (g, a, b))
    if stat.get("NCBI_Family", 0) == 0:
        print("  结论: NCBI_Family 无内部空白 -> 不对称当前无实际触发，属潜伏风险")
    else:
        print("  结论: NCBI_Family 存在内部空白 -> 不对称会实际触发，需统一归一化口径")


def probe_tables():
    print()
    print("=" * 78)
    print("探针 2: 最终 Genus 与 Species 首词（独立属主张）不一致的规模")
    print("=" * 78)
    for tag, path in TABLES:
        if not os.path.exists(path):
            print("  [跳过] %s 不存在: %s" % (tag, path))
            continue
        cols, rows = read_tsv(path)
        need = ("Genus", "Species")
        if not all(c in (cols or []) for c in need):
            print("  [跳过] %s 缺 Genus/Species 列" % tag)
            continue
        has_sagree = "Species_agree" in (cols or [])
        has_gagree = "Genus_agree" in (cols or [])
        n_binom = n_conflict = 0
        n_conf_sagree11 = n_binom_sagree11 = 0
        ex = []
        for r in rows:
            sp = (r.get("Species") or "").strip().strip('"')
            ge = (r.get("Genus") or "").strip().strip('"')
            if " " not in sp:
                continue
            first = sp.split(" ")[0]
            if not first.lower().endswith("virus") or BAD_SUFFIX.search(first) or SP_PREFIX.search(sp):
                continue
            n_binom += 1
            sa = (r.get("Species_agree") or "").strip().strip('"')
            if sa.startswith("1/1"):
                n_binom_sagree11 += 1
            if valid(ge) and first.lower() != ge.lower():
                n_conflict += 1
                if sa.startswith("1/1"):
                    n_conf_sagree11 += 1
                if len(ex) < 5:
                    ex.append((r.get("contig_id"), r.get("Family"), ge, sp,
                               (r.get("Genus_agree") or "").strip('"'),
                               sa, r.get("completeness") or ""))
        print("  %s %s" % (tag, path))
        print("    总行数 %d ; 双名法 Species 行 %d ; 属主张与最终 Genus 不一致 %d 行 (%.2f%% of 双名法)"
              % (len(rows), n_binom, n_conflict,
                 (100.0 * n_conflict / n_binom) if n_binom else 0.0))
        print("    其中 Species_agree 仍为 1/1 的 %d 行（占比 %.1f%%）；双名法行里 Species_agree=1/1 共 %d 行"
              % (n_conf_sagree11,
                 (100.0 * n_conf_sagree11 / n_conflict) if n_conflict else 0.0,
                 n_binom_sagree11))
        if ex:
            print("    样例: contig_id | Family | Genus | Species | Genus_agree | Species_agree | completeness")
            for e in ex:
                print("      %s | %s | %s | %s | %s | %s | %s" % e)
        print()


if __name__ == "__main__":
    sys.stdout.reconfigure(encoding="utf-8")
    probe_ref()
    probe_tables()
