#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Biavirus 自洽行的同库排名裁决。

背景：final_integrated_classification.tsv 里带 Biavirus 的行共 838 行，其中一批在
「逐级相容性」判据下自洽（Family=Schizomimiviridae, Genus=Biavirus, Species=Biavirus
raunefjordenense）。这批自洽行不受口径 A 校准影响，但它们出现在植物病毒组里本身可疑
（Biavirus 参考基因组 HG999358 是 1.42 Mb、宿主为原生生物的巨病毒）。

裁决方式（独立于分类表所用的工具链，且只用同一个库 RVDB v31）：
  1. 取这些行的 contig（即 04_CLUSTER/4_centroids/final_centroids.fasta 的 centroid）
  2. Prodigal 预测 ORF
  3. diamond blastp 对 U-RVDBv31.0-prot_unique 取每 query 前 60 命中
  4. sseqid -> protid2taxid.map -> taxid -> rankedlineage.dmp -> 科/属
  5. 看每条的 rank1 命中落在哪个科，以及「第一个 Biavirus/Schizomimiviridae 命中」的排名与 e 值差

只读管线产物，结果写到 calibration_20260914/biavirus_arbitration/。
"""
import csv
import os
import re
import shutil
import subprocess
import sys
from collections import Counter, defaultdict

BASE = "/home/zhangwenda/MMPV-paper/goji-virome/02_novel_virus/RNA-Lycium_barbarum_out"
VOT = os.path.join(BASE, "05_Taxonomy", "Votus.integrated")
SRC = os.path.join(VOT, "final_integrated_classification.tsv")
CENT = os.path.join(BASE, "04_CLUSTER", "4_centroids", "final_centroids.fasta")
RESCUE = os.path.join(BASE, "04_CLUSTER", "rescue_centroids.fasta")
DMP = "/home/zhangwenda/database/taxonomy/rankedlineage.dmp"
P2T = "/home/zhangwenda/database/virus-db/RVDB-v31/protid2taxid.map"
DMND = "/home/zhangwenda/database/virus-db/RVDB-v31/RVDB_viroids.diamond_db/U-RVDBv31.0-prot_unique.dmnd"
OUTDIR = os.path.join(VOT, "calibration_20260914", "biavirus_arbitration")
DIAMOND = "/home/zhangwenda/biosoft/binary/diamond"

RANKS = ["Species", "Genus", "Family", "Order", "Class", "Phylum", "Kingdom", "Realm"]
REF_IDX = {"Species": 2, "Genus": 3, "Family": 4, "Order": 5,
           "Class": 6, "Phylum": 7, "Kingdom": 8, "Realm": 9}
PAIRS = [
    ("Realm", "Kingdom", "Realm"), ("Kingdom", "Phylum", "Kingdom"),
    ("Phylum", "Class", "Phylum"), ("Class", "Order", "Class"),
    ("Order", "Family", "Order"), ("Family", "Genus", "Family"),
    ("Genus", "Species", "Genus"),
]
TARGET = ("schizomimiviridae", "biavirus")


def has_word(v, w):
    """按词边界匹配，避免 Eurybiavirus 内含 biavirus 子串被误选"""
    if not v:
        return False
    return w in [t for t in re.split(r"[^A-Za-z0-9]+", v.lower()) if t]


def norm(v):
    if v is None:
        return None
    v = v.strip().strip('"').strip()
    while v.endswith("*"):
        v = v[:-1].strip()
    if not v or v.upper() in ("NA", "N/A"):
        return None
    return v


def low(v):
    return v.lower() if v else None


def gv(d, k):
    return norm(d.get(k))


def read_table(path):
    with open(path, newline="", encoding="utf-8", errors="replace") as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        return [{k.strip().strip('"'): v for k, v in r.items()} for r in rd]


def scan_dmp_names(needed):
    """按文件顺序取 name 首次出现行；返回 name -> 字段列表（10 列）"""
    got = {}
    with open(DMP, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not line or line[0] == "#":
                continue
            parts = [p.strip().strip("|").strip() for p in line.rstrip("\n").split("\t|\t")]
            if len(parts) < 10:
                continue
            nm = parts[1].lower()
            if nm in needed and nm not in got:
                got[nm] = parts
                if len(got) == len(needed):
                    break
    return got


def scan_dmp_taxids(needed):
    got = {}
    with open(DMP, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not line or line[0] == "#":
                continue
            parts = [p.strip().strip("|").strip() for p in line.rstrip("\n").split("\t|\t")]
            if len(parts) < 10:
                continue
            if parts[0] in needed:
                got[parts[0]] = parts
    return got


def run(cmd, log):
    log.write("$ %s\n" % " ".join(cmd))
    log.flush()
    p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    out = p.stdout.decode("utf-8", "replace")
    log.write(out[-4000:] + "\n")
    log.flush()
    return p.returncode


def read_fasta(path, want=None):
    """want=None 读全部；否则只读 header 第一段在 want 里的序列。返回 (id -> seq, 文件内出现顺序)"""
    seqs = {}
    order = []
    cid, buf = None, []
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                if cid is not None and (want is None or cid in want):
                    seqs[cid] = "".join(buf)
                    order.append(cid)
                cid = line[1:].split()[0]
                buf = []
            elif cid is not None:
                buf.append(line.strip())
        if cid is not None and (want is None or cid in want):
            seqs[cid] = "".join(buf)
            order.append(cid)
    return seqs, order


def main():
    if not os.path.isdir(OUTDIR):
        os.makedirs(OUTDIR)
    log = open(os.path.join(OUTDIR, "run.log"), "w", encoding="utf-8")

    rows = read_table(SRC)
    log.write("原表 %d 行\n" % len(rows))

    # ── 1. 选出所有含 Biavirus 的行并按逐级相容性分类 ──
    involved = []
    needed = set()
    for d in rows:
        vals = [gv(d, r) for r in RANKS]
        if any(has_word(v, "biavirus") for v in vals):
            involved.append(d)
            for v in vals:
                needed.add(v.lower())
    gen_bia = [d for d in involved if has_word(gv(d, "Genus"), "biavirus")]
    log.write("含 Biavirus 词的行 %d（其中 Genus 列就是 Biavirus 的 %d），待查 name %d\n"
              % (len(involved), len(gen_bia), len(needed)))
    ref = scan_dmp_names(needed)
    log.write("dmp 命中 name %d\n" % len(ref))

    bucket = Counter()
    types = Counter()
    selected = []
    detail = []
    for d in involved:
        bad = []
        for hi, lo, refcol in PAIRS:
            vh, vl = gv(d, hi), gv(d, lo)
            if not vh or not vl:
                continue
            rl = ref.get(vl.lower())
            if not rl:
                continue
            exp = rl[REF_IDX[refcol]] or ""
            if exp and exp.lower() != vh.lower():
                bad.append("%s<->%s" % (hi, lo))
        badstr = ";".join(bad)
        fam = gv(d, "Family") or "-"
        gen = gv(d, "Genus") or "-"
        sp = gv(d, "Species") or "-"
        key = "no_contradiction" if not bad else badstr
        bucket[key] += 1
        types[(fam, gen, sp, key)] += 1
        detail.append((d["contig_id"], fam, gen, sp, d.get("primary_tool") or "",
                       d.get("Family_agree") or "", d.get("Genus_agree") or "",
                       d.get("Species_agree") or "", badstr))
        if not bad:
            selected.append(d)

    log.write("\n=== 含 Biavirus 行的相容性分型 ===\n")
    for k, v in bucket.most_common():
        log.write("  %-40s %5d\n" % (k, v))
    log.write("\n=== 自洽行（待裁决）的 科/属/种 构成 ===\n")
    for (fam, gen, sp, k), v in types.most_common():
        if k == "no_contradiction":
            log.write("  %-24s %-20s %-34s %5d\n" % (fam, gen, sp, v))
    log.write("\n选定自洽行 %d 条\n" % len(selected))

    with open(os.path.join(OUTDIR, "select_diagnostics.tsv"), "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["contig_id", "Family", "Genus", "Species", "primary_tool",
                    "Family_agree", "Genus_agree", "Species_agree", "contradictions", "selected"])
        for rec in detail:
            w.writerow(list(rec) + ["1" if not rec[8] else "0"])

    ids = [d["contig_id"] for d in selected]
    with open(os.path.join(OUTDIR, "selected_ids.txt"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(ids) + "\n")

    # ── 2. 抽序列 ──
    want = set(ids)
    seqs, _ = read_fasta(CENT, want)
    miss = [i for i in ids if i not in seqs]
    if miss:
        s2, _ = read_fasta(RESCUE, set(miss))
        seqs.update(s2)
        miss = [i for i in ids if i not in seqs]
    log.write("抽到序列 %d / %d，缺 %d\n" % (len(seqs), len(ids), len(miss)))
    if miss:
        log.write("  缺失样例: %s\n" % ", ".join(miss[:5]))
    SEL = os.path.join(OUTDIR, "biavirus_sel.fasta")
    with open(SEL, "w", encoding="utf-8") as fh:
        for i in ids:
            if i in seqs:
                fh.write(">%s\n%s\n" % (i, seqs[i]))
    tot = sum(len(v) for v in seqs.values())
    log.write("序列总长 %.2f Mb，平均 %.0f bp\n" % (tot / 1e6, tot / max(1, len(seqs))))

    # ── 3. Prodigal ──
    PRO = os.path.join(OUTDIR, "proteins.faa")
    prodigal = shutil.which("prodigal") or "/home/zhangwenda/mambaforge/bin/prodigal"
    if not os.path.exists(PRO) or os.path.getsize(PRO) == 0:
        rc = run([prodigal, "-i", SEL, "-a", PRO, "-p", "meta",
                  "-o", os.path.join(OUTDIR, "prodigal.gbk"),
                  "-d", os.path.join(OUTDIR, "orfs.fna")], log)
        log.write("prodigal rc=%d\n" % rc)
        if rc != 0:
            log.write("FATAL: prodigal 失败\n")
            return 1
    prot, porder = read_fasta(PRO)
    log.write("ORF %d 条\n" % len(prot))

    # ── 4. diamond ──
    HITS = os.path.join(OUTDIR, "hits.tsv")
    if not os.path.exists(HITS) or os.path.getsize(HITS) == 0:
        rc = run([DIAMOND, "blastp", "-q", PRO, "-d", DMND, "-o", HITS,
                  "-k", "60", "-e", "1e-3", "--sensitive", "-f", "6",
                  "-p", "16"], log)
        log.write("diamond rc=%d\n" % rc)
        if rc != 0 or os.path.getsize(HITS) == 0:
            log.write("FATAL: diamond 失败\n")
            return 1

    # ── 5. sseqid -> taxid -> 科属 ──
    a2t = {}
    with open(P2T, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 2:
                a2t[p[0]] = p[1]
    log.write("protid2taxid 载入 %d 条\n" % len(a2t))

    hit_rows = []
    need_tax = set()
    unmapped = 0
    with open(HITS, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            # diamond -f 6 默认 12 列，末两列固定是 evalue / bitscore
            if len(f) < 6:
                continue
            q, s, pid, ln = f[0], f[1], f[2], f[3]
            try:
                ev, bs = float(f[-2]), float(f[-1])
            except ValueError:
                continue
            tx = a2t.get(s)
            if tx is None:
                unmapped += 1
            else:
                need_tax.add(tx)
            hit_rows.append([q, s, pid, ln, ev, bs, tx])
    log.write("命中行 %d，sseqid 无 taxid 映射 %d，涉及 taxid %d\n"
              % (len(hit_rows), unmapped, len(need_tax)))

    txl = scan_dmp_taxids(need_tax)
    log.write("dmp 命中 taxid %d\n" % len(txl))

    def lineage(tx):
        """dmp 列: 0 tax_id,1 tax_name,2 species,3 genus,4 family,5 order,6 class,7 phylum,8 kingdom,9 superkingdom
        属级 taxid 的 genus 列本身为空，用 tax_name 兜底"""
        if tx is None or tx not in txl:
            return ("", "", "", "")
        p = txl[tx]
        gen = p[3] or (p[1] if not p[2] else "")
        return (p[4], gen, p[2], p[9])  # family, genus, species, realm

    with open(os.path.join(OUTDIR, "hit_families.tsv"), "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["qseqid", "sseqid", "pident", "aln_len", "evalue", "bitscore",
                    "taxid", "family", "genus", "species", "realm"])
        for h in hit_rows:
            fam, gen, sp, realm = lineage(h[6])
            w.writerow(h[:6] + [h[6] or ""] + [fam, gen, sp, realm])

    # ── 6. 逐 contig 汇总 ──
    cid_of = {}
    for q in prot:
        cid_of[q] = q.rsplit("_", 1)[0] if q.rsplit("_", 1)[-1].isdigit() else q
    per = defaultdict(list)
    for h in hit_rows:
        c = cid_of.get(h[0])
        if c is None:
            c = h[0]
        per[c].append(h)

    out_rows = []
    verdict = Counter()
    rank1_fam = Counter()
    rowfam_hit = 0
    has_bia = 0
    bia_rank = Counter()
    gap_list = []
    for d in selected:
        cid = d["contig_id"]
        hs = per.get(cid, [])
        rowfam = gv(d, "Family") or "-"
        if not hs:
            verdict["no_hit"] += 1
            out_rows.append([cid, rowfam, 0, "-", "", "", "", "", "", ""])
            continue
        hs.sort(key=lambda x: (-x[5], x[4]))
        top = hs[0]
        tf, tg, ts, tr = lineage(top[6])
        rank1_fam[tf or "(unassigned)"] += 1
        if (tf or "").lower() == rowfam.lower():
            rowfam_hit += 1
        bi = None
        for i, h in enumerate(hs, 1):
            f2, g2, _, _ = lineage(h[6])
            if (f2 or "").lower() in TARGET or (g2 or "").lower() in TARGET:
                bi = (i, h)
                break
        if bi is not None:
            has_bia += 1
            bia_rank[bi[0]] += 1
            gap = (bi[1][4] / top[4]) if top[4] > 0 and bi[1][4] > 0 else 0
            gap_list.append((gap, cid, top[4], bi[1][4], bi[0], tf, rowfam))
        top10 = Counter()
        for h in hs[:10]:
            f2, _, _, _ = lineage(h[6])
            top10[f2 or "(unassigned)"] += 1
        comp = ",".join("%s:%d" % (k, v) for k, v in top10.most_common(3))
        if (tf or "").lower() == rowfam.lower():
            v = "rank1_matches_row_family"
        elif bi is not None and tf and tf.lower() in TARGET:
            v = "rank1_is_target"
        else:
            v = "rank1_other_family"
        verdict[v] += 1
        out_rows.append([cid, rowfam, len(hs), tf or "-", tg or "-", ts or "-",
                         "%.2e" % top[4], top[1], bi[0] if bi else "",
                         "%.3g" % (bi[1][4] / top[4]) if bi and top[4] > 0 else "", comp])

    with open(os.path.join(OUTDIR, "per_contig_summary.tsv"), "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["contig_id", "row_Family", "n_hits", "rank1_family", "rank1_genus",
                    "rank1_species", "rank1_evalue", "rank1_sseqid", "first_target_rank",
                    "evalue_ratio_target_vs_rank1", "top10_family_comp"])
        w.writerows(out_rows)

    log.write("\n=== 裁决结果（自洽行 %d 条）===\n" % len(selected))
    for k, v in verdict.most_common():
        log.write("  %-28s %5d (%.1f%%)\n" % (k, v, 100.0 * v / max(1, len(selected))))
    log.write("  rank1 科 == 行内科: %d (%.1f%%)\n" % (rowfam_hit, 100.0 * rowfam_hit / max(1, len(selected))))
    log.write("  命中里出现过 Biavirus/Schizomimiviridae 的: %d (%.1f%%)\n"
              % (has_bia, 100.0 * has_bia / max(1, len(selected))))
    log.write("\n  rank1 科分布:\n")
    for k, v in rank1_fam.most_common(12):
        log.write("    %-28s %5d\n" % (k, v))
    log.write("\n  第一个目标科命中的排名分布:\n")
    for k, v in sorted(bia_rank.items()):
        log.write("    rank %-4d %5d\n" % (k, v))
    gap_list.sort()
    if gap_list:
        mid = gap_list[len(gap_list) // 2][0]
        log.write("\n  e 值比（目标科最好命中 / rank1）中位数 %.3g，最小 %.3g，最大 %.3g\n"
                  % (mid, gap_list[0][0], gap_list[-1][0]))
        log.write("  样例（差的）:\n")
        for g in gap_list[-3:]:
            log.write("    %s rank1=%s(%.2e) 目标科 rank%d(%.2e) 行内科=%s 比值%.3g\n"
                      % (g[1], g[5], g[2], g[4], g[3], g[6], g[0]))
    log.close()
    return 0


if __name__ == "__main__":
    sys.exit(main())
