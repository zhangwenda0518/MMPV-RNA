#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Biavirus 自洽行裁决的第二套口径与合并报告。

第一套口径（arbitrate_biavirus.py）：自己对 465 条 centroid 跑 Prodigal（802 ORF）。
第二套口径（本脚本）：改用管线自己的 ORF 预测（Votus.classed/cat_output/
CAT_output.predicted_proteins.faa，CAT 跑 Prodigal 的产物），保证与被裁决的分类结果
用的是同一套 ORF 单元。

两套都做「同一库（RVDB v31）内排名」这件事，然后：
  - 按 目标科(Schizomimiviridae/Biavirus) 在 top10 命中里的占比分层
  - 报告 rank1 的 科/门（门=Nucleocytoviricota 即仍是巨病毒样）
  - 与行内 primary_tool / Genus_agree（哪个工具投了 Biavirus）交叉
  - 两套口径一致性矩阵
输出 ARBITRATION_REPORT.md + tier_table.tsv + hits_cat.tsv。
"""
import csv
import os
import shutil
import subprocess
import sys
from collections import Counter, defaultdict

BASE = "/home/zhangwenda/MMPV-paper/goji-virome/02_novel_virus/RNA-Lycium_barbarum_out"
VOT = os.path.join(BASE, "05_Taxonomy", "Votus.integrated")
OUTDIR = os.path.join(VOT, "calibration_20260914", "biavirus_arbitration")
SRC = os.path.join(VOT, "final_integrated_classification.tsv")
CATFAA = os.path.join(BASE, "05_Taxonomy", "Votus.classed", "cat_output",
                      "CAT_output.predicted_proteins.faa")
DMP = "/home/zhangwenda/database/taxonomy/rankedlineage.dmp"
P2T = "/home/zhangwenda/database/virus-db/RVDB-v31/protid2taxid.map"
DMND = "/home/zhangwenda/database/virus-db/RVDB-v31/RVDB_viroids.diamond_db/U-RVDBv31.0-prot_unique.dmnd"
DIAMOND = "/home/zhangwenda/biosoft/binary/diamond"
PRO = os.path.join(OUTDIR, "proteins.faa")
HITS_A = os.path.join(OUTDIR, "hits.tsv")
CATPRO = os.path.join(OUTDIR, "cat_proteins.faa")
HITS_B = os.path.join(OUTDIR, "hits_cat.tsv")
REP = os.path.join(OUTDIR, "ARBITRATION_REPORT.md")
NCLDV_PHYLUM = "nucleocytoviricota"
TGT_FAM = "schizomimiviridae"
TGT_GEN = "biavirus"


def norm(v):
    if v is None:
        return None
    v = v.strip().strip('"').strip()
    while v.endswith("*"):
        v = v[:-1].strip()
    if not v or v.upper() in ("NA", "N/A"):
        return None
    return v


def read_fasta(path, keep=None):
    seqs, cid, buf = {}, None, []
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                if cid and (keep is None or cid in keep):
                    seqs[cid] = "".join(buf)
                cid = line[1:].split()[0]
                buf = []
            elif cid:
                buf.append(line.strip())
        if cid and (keep is None or cid in keep):
            seqs[cid] = "".join(buf)
    return seqs


def scan_dmp(needed):
    got = {}
    with open(DMP, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not line or line[0] == "#":
                continue
            p = [x.strip().strip("|").strip() for x in line.rstrip("\n").split("\t|\t")]
            if len(p) < 10:
                continue
            if p[0] in needed:
                got[p[0]] = p
    return got


def load_a2t():
    a2t = {}
    with open(P2T, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 2:
                a2t[p[0]] = p[1]
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
    ids = [l.strip() for l in open(os.path.join(OUTDIR, "selected_ids.txt"), encoding="utf-8") if l.strip()]
    idset = set(ids)
    log = open(os.path.join(OUTDIR, "summarize.log"), "w", encoding="utf-8")

    # ── 第二套 ORF：管线自己的 CAT 预测蛋白 ──
    if not os.path.exists(CATPRO) or os.path.getsize(CATPRO) == 0:
        with open(CATFAA, encoding="utf-8", errors="replace") as fh, \
             open(CATPRO, "w", encoding="utf-8") as out:
            cur = 0
            for line in fh:
                if line.startswith(">"):
                    tok = line[1:].split()[0]
                    cid = tok.rsplit("_", 1)[0]
                    cur = 1 if (cid in idset and tok.rsplit("_", 1)[-1].isdigit()) else 0
                if cur:
                    out.write(line)
    cat = read_fasta(CATPRO)
    log.write("CAT 口径 ORF %d 条（覆盖 %d 条 contig）\n"
              % (len(cat), len(set(k.rsplit("_", 1)[0] for k in cat))))

    if not os.path.exists(HITS_B) or os.path.getsize(HITS_B) == 0:
        cmd = [DIAMOND, "blastp", "-q", CATPRO, "-d", DMND, "-o", HITS_B,
               "-k", "60", "-e", "1e-3", "--sensitive", "-f", "6", "-p", "16"]
        log.write("$ %s\n" % " ".join(cmd))
        p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        log.write(p.stdout.decode("utf-8", "replace")[-2000:] + "\n")
        if p.returncode != 0:
            log.write("FATAL diamond rc=%d\n" % p.returncode)
            return 1
    log.flush()

    # ── taxid -> 谱系 ──
    a2t = load_a2t()
    runs = {}
    need = set()
    for tag, hp in (("prodigal", HITS_A), ("cat", HITS_B)):
        rows = parse_hits(hp)
        runs[tag] = rows
        for r in rows:
            tx = a2t.get(r[1])
            if tx:
                need.add(tx)
        log.write("%s: 命中 %d 行\n" % (tag, len(rows)))
    txl = scan_dmp(need)
    log.write("dmp 命中 taxid %d\n" % len(txl))

    def lin(tx):
        if not tx or tx not in txl:
            return ("", "", "", "")
        p = txl[tx]
        gen = p[3] or (p[1] if not p[2] else "")
        return (p[4], gen, p[7], p[9])  # family, genus, phylum, realm

    def per_contig(rows):
        per = defaultdict(list)
        n_unmapped = Counter()
        for r in rows:
            cid = r[0].rsplit("_", 1)[0]
            tx = a2t.get(r[1])
            fam, gen, phy, _ = lin(tx)
            if not tx:
                n_unmapped["nomap"] += 1
            per[cid].append((r[5], r[4], r[1], r[2], fam, gen, phy))
        for cid in per:
            per[cid].sort(key=lambda x: (-x[0], x[1]))
        return per, n_unmapped

    pa, _ = per_contig(runs["prodigal"])
    pb, _ = per_contig(runs["cat"])

    # ── 行内信息 ──
    rows = list(csv.DictReader(open(SRC, newline="", encoding="utf-8", errors="replace"),
                               delimiter="\t"))
    rowinfo = {}
    for d in rows:
        k = d["contig_id"].strip().strip('"')
        if k in idset:
            rowinfo[k] = {kk: (v or "").strip().strip('"') for kk, v in d.items()}

    def tier_of(contig, per, orfcnt):
        hs = per.get(contig, [])
        n_orf = orfcnt.get(contig, 0)
        if n_orf == 0:
            return ("D_noorf", 0, 0, "", "", "", 0, 0, 0, "")
        if not hs:
            return ("E_nohit", 0, 0, "", "", "", 0, 0, 0, "")
        top10 = hs[:10]
        n_tgt = sum(1 for h in top10 if h[4].lower() == TGT_FAM or h[5].lower() == TGT_GEN)
        n_ncldv = sum(1 for h in top10 if h[6].lower() == NCLDV_PHYLUM)
        best_tgt = next((h for h in hs if h[4].lower() == TGT_FAM or h[5].lower() == TGT_GEN), None)
        r1 = hs[0]
        if n_tgt >= 5:
            t = "A_strong"
        elif n_tgt >= 1:
            t = "B_partial"
        else:
            t = "C_negative"
        return (t, n_orf, len(hs), r1[4], r1[5], r1[6], n_tgt, n_ncldv,
                best_tgt[3] if best_tgt else "", best_tgt[4] if best_tgt else "")

    # ORF 计数：分别按各自 ORF 集合统计
    ofa = Counter(k.rsplit("_", 1)[0] for k in read_fasta(PRO).keys())
    ofb = Counter(k.rsplit("_", 1)[0] for k in cat.keys())
    resA = {c: tier_of(c, pa, ofa) for c in ids}
    resB = {c: tier_of(c, pb, ofb) for c in ids}

    def dist(res):
        return Counter(v[0] for v in res.values())

    log.write("\nprodigal 口径分层: %s\n" % dict(dist(resA)))
    log.write("cat 口径分层: %s\n" % dict(dist(resB)))

    # ── 交叉：与工具投票 ──
    tool_by_tier = defaultdict(Counter)
    agree_by_tier = defaultdict(Counter)
    for c in ids:
        ri = rowinfo.get(c, {})
        t = resB[c][0]
        tool_by_tier[t][ri.get("primary_tool", "-")] += 1
        ag = ri.get("Genus_agree", "-")
        ag2 = ag.split(":")[-1] if ":" in ag else ag
        agree_by_tier[t][ag2] += 1

    # ── 一致性（两口径的分层方向）──
    agree_mat = Counter()
    for c in ids:
        a = resA[c][0]
        b = resB[c][0]
        agree_mat[(a, b)] += 1

    # ── 输出 ──
    with open(os.path.join(OUTDIR, "tier_table.tsv"), "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["contig_id", "row_Family", "row_Genus", "primary_tool", "Genus_agree",
                    "Family_agree",
                    "cat_tier", "cat_n_orf", "cat_n_hits", "cat_rank1_family", "cat_rank1_genus",
                    "cat_rank1_phylum", "cat_top10_target", "cat_top10_ncldv", "cat_best_target_rank",
                    "cat_best_target_family",
                    "prodigal_tier", "prodigal_n_orf", "prodigal_n_hits", "prodigal_rank1_family",
                    "prodigal_rank1_genus", "prodigal_rank1_phylum", "prodigal_top10_target",
                    "prodigal_top10_ncldv"])
        for c in ids:
            ri = rowinfo.get(c, {})
            a, b = resA[c], resB[c]
            # 目标科命中排名
            def tgt_rank(per):
                hs = per.get(c, [])
                for i, h in enumerate(hs, 1):
                    if h[4].lower() == TGT_FAM or h[5].lower() == TGT_GEN:
                        return i
                return ""
            w.writerow([c, ri.get("Family", ""), ri.get("Genus", ""), ri.get("primary_tool", ""),
                        ri.get("Genus_agree", ""), ri.get("Family_agree", ""),
                        b[0], b[1], b[2], b[3], b[4], b[5], b[6], b[7], tgt_rank(pb), b[9],
                        a[0], a[1], a[2], a[3], a[4], a[5], a[6], a[7]])

    L = []
    A = L.append
    A("# Biavirus 自洽行 同库排名裁决报告")
    A("")
    A("对象：final_integrated_classification.tsv 中含 Biavirus 词、且在逐级相容性判据下")
    A("自洽（无相邻等级矛盾）的 **465 行**，全部为 Family=Schizomimiviridae + Genus=Biavirus。")
    A("")
    A("判据：这些 contig 的 ORF 在同一参考库 RVDB v31（U-RVDBv31.0-prot_unique）里的排名。")
    A("两套 ORF 口径并行，防止 ORF 预测差异干扰结论：")
    A("  口径1 = 对同一批 centroid 现场跑 Prodigal -p meta；")
    A("  口径2 = 管线自己 CAT 预测的蛋白（Votus.classed/cat_output/CAT_output.predicted_proteins.faa）。")
    A("")
    A("分层定义（目标科 = Schizomimiviridae / 属 = Biavirus）：")
    A("  A_strong  = 目标科占该 contig top10 命中的 ≥5 条")
    A("  B_partial = 目标科在 top10 中出现 1-4 条")
    A("  C_negative= 有命中但 top10 中没有目标科")
    A("  E_nohit   = 该 contig 的 ORF 在库里 1e-3 内无任何命中")
    A("  D_noorf   = 该 contig 没预测出 ORF（无法用此法裁决）")
    A("")
    A("## 一、分层结果")
    A("")
    A("| 分层 | 口径2 CAT-ORF | 口径1 Prodigal-ORF |")
    A("|---|---|---|")
    for t in ("A_strong", "B_partial", "C_negative", "E_nohit", "D_noorf"):
        A("| %s | %d (%.1f%%) | %d (%.1f%%) |"
          % (t, dist(resB)[t], 100.0 * dist(resB)[t] / len(ids),
             dist(resA)[t], 100.0 * dist(resA)[t] / len(ids)))
    A("")
    A("两套口径分层交叉：")
    A("")
    A("| 口径2\\口径1 | A | B | C | E | D |")
    A("|---|---|---|---|---|---|")
    for b in ("A_strong", "B_partial", "C_negative", "E_nohit", "D_noorf"):
        cells = [str(agree_mat[(a, b)]) for a in ("A_strong", "B_partial", "C_negative", "E_nohit", "D_noorf")]
        A("| %s | %s |" % (b, " | ".join(cells)))
    A("")
    A("## 二、rank1 命中的科构成（口径2）")
    A("")
    c1 = Counter(v[3] or "(未定型)" for v in resB.values() if v[0] not in ("E_nohit", "D_noorf"))
    for k, v in c1.most_common():
        A("- %s: %d" % (k, v))
    A("")
    ph = Counter(v[5] or "(未定型)" for v in resB.values() if v[0] not in ("E_nohit", "D_noorf"))
    ncldv = sum(v for k, v in ph.items() if k.lower() == NCLDV_PHYLUM)
    tot = sum(ph.values())
    A("rank1 命中落在门 Nucleocytoviricota（即仍是巨病毒样序列）的：%d / %d (%.1f%%)"
      % (ncldv, tot, 100.0 * ncldv / max(1, tot)))
    A("")
    A("## 三、Biavirus 这个属是谁投的票（按口径2分层交叉）")
    A("")
    A("| 分层 | primary_tool 分布 | Genus_agree 的投票工具 |")
    A("|---|---|---|")
    for t in ("A_strong", "B_partial", "C_negative", "E_nohit", "D_noorf"):
        A("| %s | %s | %s |" % (t,
                                ", ".join("%s:%d" % kv for kv in tool_by_tier[t].most_common(4)),
                                ", ".join("%s:%d" % kv for kv in agree_by_tier[t].most_common(4))))
    A("")
    A("## 四、读法")
    A("")
    A("1. A_strong 层是 Biavirus 归属真正站得住的证据：这些 contig 的 ORF 在同一库里")
    A("   最好的那批命中直接就是 Schizomimiviridae/Biavirus。")
    A("2. C_negative 层值得单独看 rank1 落在哪个门：若仍是 Nucleocytoviricota 的其它科")
    A("   （Mimiviridae/Phycodnaviridae/Pithoviridae 等），说明序列确属巨病毒，只是属级")
    A("   分辨率不够或与 Biavirus 参考距离更远，不宜直接判为错；若落到完全无关的科，")
    A("   则属级/种级赋值可疑。")
    A("3. E_nohit 层不能当作反证：短 contig 上的 ORF 在 1e-3 阈值下无命中，也可能只是")
    A("   分歧度高或 ORF 边界不准。要判错需要额外证据。")
    A("4. 两套 ORF 口径分层方向一致时，结论不依赖 ORF caller。")

    open(REP, "w", encoding="utf-8").write("\n".join(L) + "\n")
    log.write("报告写出 %s\n" % REP)
    log.close()
    return 0


if __name__ == "__main__":
    sys.exit(main())
