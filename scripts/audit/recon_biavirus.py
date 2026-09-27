#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""465 行 Biavirus 裁决的第二层：属级分辨力与证据来源。

关键背景（实测）：
  - Schizomimiviridae 在 NCBI dmp 下有 2 个属：Biavirus、Kratosvirus
    => 属级 Biavirus 不是由科级自动推出，需要独立证据
  - CAT/diamond_lca/mmseqs 三者 Biavirus 命中数 = 0（只到科级）
    ACVirus 720 行、metabuli 144 行提到 Biavirus

本脚本回答：
  1. 465 行里 rank1 命中到底是 属级Biavirus / 属级Kratosvirus / 族级(无属) / 其他
  2. 目标科(Schizomimiviridae)命中的属级构成
  3. 同一命中的 e 值/一致性，看 Biavirus 属级命中是不是最强证据
  4. ACVirus 对这 465 行的 Best_Contig + Max_Coverage 分布（它自己的证据强度）
  5. 逐工具在 465 行上的 科/属 投票
"""
import csv
import os
import subprocess
import sys
from collections import Counter, defaultdict

BASE = "/home/zhangwenda/MMPV-paper/goji-virome/02_novel_virus/RNA-Lycium_barbarum_out"
VOT = os.path.join(BASE, "05_Taxonomy", "Votus.integrated")
CAL = os.path.join(VOT, "calibration_20260914")
BAD = os.path.join(CAL, "biavirus_arbitration")
CLS = os.path.join(BASE, "05_Taxonomy", "Votus.classed")
SRC = os.path.join(VOT, "final_integrated_classification.tsv")
DMP = "/home/zhangwenda/database/taxonomy/rankedlineage.dmp"
P2T = "/home/zhangwenda/database/virus-db/RVDB-v31/protid2taxid.map"
OUT = sys.stdout


def p(*a):
    print(*a)
    OUT.flush()


def scan_dmp_all():
    """返回 {taxid: (name, species, genus, family, order, class, phylum, kingdom, realm)}"""
    d = {}
    with open(DMP, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not line or line[0] == "#":
                continue
            f = [x.strip().strip("|").strip() for x in line.rstrip("\n").split("\t|\t")]
            if len(f) < 10:
                continue
            d[f[0]] = tuple(f[1:10])
    return d


def load_a2t(need=None):
    a2t = {}
    with open(P2T, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 2:
                a2t[f[0]] = f[1]
    return a2t


def parse_hits(path):
    out = []
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 6:
                continue
            try:
                ev, bs = float(f[-2]), float(f[-1])
            except ValueError:
                continue
            out.append((f[0], f[1], f[2], f[3], ev, bs))
    return out


def main():
    ids = [l.strip() for l in open(os.path.join(BAD, "selected_ids.txt"), encoding="utf-8") if l.strip()]
    idset = set(ids)
    dmp = scan_dmp_all()
    a2t = load_a2t()

    def lin(tx):
        r = dmp.get(tx)
        if not r:
            return None
        name, sp, gen, fam, order, cls, phy, king, realm = r
        return {"name": name, "species": sp, "genus": gen, "family": fam,
                "order": order, "class": cls, "phylum": phy, "realm": realm}

    # ── 1+2+3 用管线自己 ORF 的命中（cat 口径）为主，prodigal 对照 ──
    for tag, hp in (("CAT-ORF", os.path.join(BAD, "hits_cat.tsv")),
                    ("Prodigal-ORF", os.path.join(BAD, "hits.tsv"))):
        rows = parse_hits(hp)
        per = defaultdict(list)
        for r in rows:
            cid = r[0].rsplit("_", 1)[0]
            if cid in idset:
                per[cid].append(r)
        for cid in per:
            per[cid].sort(key=lambda x: (-x[5], x[4]))

        p("\n" + "=" * 70)
        p("【口径 %s】465 行的 rank1 命中细粒度归属" % tag)
        p("=" * 70)
        c_r1 = Counter()
        c_kind = Counter()
        tgt_genus = Counter()
        for cid in ids:
            hs = per.get(cid)
            if not hs:
                c_r1["(无命中)"] += 1
                c_kind["(无命中)"] += 1
                continue
            tx = a2t.get(hs[0][1])
            L = lin(tx)
            if not L:
                c_r1["(taxid未映射)"] += 1
                c_kind["(taxid未映射)"] += 1
                continue
            fam, gen = L["family"], L["genus"]
            if gen:
                c_r1["属 %s (科 %s)" % (gen, fam)] += 1
            elif fam:
                c_r1["科 %s（该节点的属列为空）" % fam] += 1
            else:
                c_r1["%s" % L["name"]] += 1
            c_kind["门 %s" % (L["phylum"] or "(空)")] += 1
            # 目标科命中里出现过的属
            seen = set()
            for h in hs:
                L2 = lin(a2t.get(h[1]))
                if L2 and L2["family"] == "Schizomimiviridae":
                    seen.add(L2["genus"] or "(族级无属)")
            for g in seen:
                tgt_genus[g] += 1
        p("\n-- rank1 细粒度（前 15）--")
        for k, v in c_r1.most_common(15):
            p("   %-58s %4d" % (k, v))
        p("\n-- 命中属级出现过的 Schizomimiviridae 属（按 contig 计）--")
        for k, v in tgt_genus.most_common():
            p("   %-24s %4d" % (k, v))

    # ── 4 ACVirus 自己的证据 ──
    p("\n" + "=" * 70)
    p("【ACVirus 自身证据】465 行的 Best_Contig 与 Max_Coverage")
    p("=" * 70)
    ac_med = os.path.join(CLS, "ACVirus_results", "Votus.acvirus", "medium_result.csv")
    ac_fin = os.path.join(CLS, "ACVirus_results", "Votus.acvirus", "final_result.tsv")
    fin = {}
    if os.path.exists(ac_fin):
        rd = csv.DictReader(open(ac_fin, newline="", encoding="utf-8", errors="replace"), delimiter="\t")
        for d in rd:
            fin[d["Nucleotide"].strip().strip('"')] = d
    if os.path.exists(ac_med):
        rd = csv.DictReader(open(ac_med, newline="", encoding="utf-8", errors="replace"), delimiter="\t")
        cov = []
        best = Counter()
        gen_vote = Counter()
        n = 0
        for d in rd:
            k = d["Nucleotide"].strip().strip('"')
            if k not in idset:
                continue
            n += 1
            bc = (d.get("Best_Contig") or "").strip().strip('"')
            best[bc] += 1
            gen_vote[d.get("Genus", "").strip().strip('"') or "(空)"] += 1
            mc = (d.get("Max_Coverage(%)") or "").replace("%", "").strip()
            try:
                cov.append(float(mc))
            except ValueError:
                pass
        p("ACVirus medium_result 覆盖 %d / 465 行" % n)
        cov.sort()
        if cov:
            def q(f):
                return cov[min(len(cov) - 1, int(len(cov) * f))]
            p("Max_Coverage%%: min %.2f  p10 %.2f  中位 %.2f  p90 %.2f  max %.2f"
              % (cov[0], q(.10), q(.50), q(.90), cov[-1]))
            for lo, hi in ((0, 20), (20, 40), (40, 60), (60, 80), (80, 101)):
                p("   覆盖 %3d-%3d%%: %d" % (lo, hi, sum(1 for x in cov if lo <= x < hi)))
        p("ACVirus 给这 465 行的属投票: %s" % dict(gen_vote.most_common(6)))
        p("ACVirus Best_Contig 前 10（参考序列）:")
        for k, v in best.most_common(10):
            p("   %-24s %4d" % (k, v))
        # 这 465 行在 final_result 里的属分布
        fg = Counter()
        for k in idset:
            if k in fin:
                fg[fin[k].get("Genus", "").strip().strip('"') or "(空)"] += 1
        p("ACVirus final_result 属分布: %s" % dict(fg.most_common(6)))

    # ── 5 逐工具投票 ──
    p("\n" + "=" * 70)
    p("【逐工具投票】465 行上各工具给的 科 / 属")
    p("=" * 70)
    comb = os.path.join(CLS, "Votus_combined_taxonomy.tsv")
    fam_v = defaultdict(Counter)
    gen_v = defaultdict(Counter)
    for d in csv.DictReader(open(comb, newline="", encoding="utf-8", errors="replace"), delimiter="\t"):
        k = d["seq_name"].strip().strip('"')
        if k not in idset:
            continue
        t = d["tool"].strip().strip('"')
        fam_v[t][(d.get("Family") or "").strip().strip('"') or "(空)"] += 1
        gen_v[t][(d.get("Genus") or "").strip().strip('"') or "(空)"] += 1
    for t in sorted(fam_v):
        p("\n-- %s --" % t)
        p("   科: %s" % dict(fam_v[t].most_common(4)))
        p("   属: %s" % dict(gen_v[t].most_common(4)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
