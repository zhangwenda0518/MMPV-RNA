#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
===============================================================================
SDT Genus Matrix — ACVirus 结果驱动的属级全基因组相似度矩阵
===============================================================================
适配无参考 GenBank 场景：
  输入不是参考 GenBank + 序列目录，而是 09b_ACVirus_Analysis 的分类结果。
  用户 contigs 的属级归属直接取自 acvirus_classify/final_result_with_confidence.tsv，
  参考序列按 virus.tree.py 的抽样逻辑从 acvirus_db 取同属/同科代表。

流程：
  1. 从 final_result_with_confidence.tsv 选出目标属的用户 contigs（含置信度）
  2. 按 virus.tree.py 口径从 acvirus_db/taxa.txt 抽样参考 accession：
       simple    : 目标属全部参考 + 同科其他属各最多 N 条（外群）
       phylogeny : 目标属按比例 + 同科其他属各最多 N 条
       lineage   : 目标种全部参考（需 --species）
  3. 从 all_virus.fasta 提取参考序列；用户 contigs 做 k-mer 自动定向校正
  4. 调用 virus_auto_pipeline 的 build_mat_sdt_exact（MAFFT 逐对全局比对，
     SDT Get_Similarity 公式）→ 层级聚类排序 → 热图 + 分布图

依赖：与 virus_auto_pipeline.py 同目录运行（服务器 MMPV-RNA/virome_analysis_pipeline）。
外部工具：mafft（PATH）。Python 侧复用 base 环境已有包。
===============================================================================
"""
import sys
import os
import csv
import argparse
import hashlib
import warnings
from pathlib import Path
from collections import Counter, defaultdict

sys.path.insert(0, str(Path(__file__).resolve().parent))

warnings.simplefilter("ignore")

import numpy as np                      # noqa: E402
import pandas as pd                     # noqa: E402
import matplotlib.pyplot as plt         # noqa: E402
from Bio.Seq import Seq                 # noqa: E402
from Bio.SeqRecord import SeqRecord     # noqa: E402
from Bio import SeqIO                   # noqa: E402

import virus_auto_pipeline as vap       # noqa: E402  (forces Agg backend)

plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42

K_ORIENT = 21


# ---------------------------------------------------------------- utilities --
def log(msg):
    print(f"[sdt-genus] {msg}", flush=True)


def sha_seq(s):
    return hashlib.sha1(s.upper().encode()).hexdigest()


def n_percent(seq):
    seq = seq.upper()
    return 100.0 * seq.count("N") / len(seq) if len(seq) else 100.0


def kmer_set(seq, k=K_ORIENT):
    s = seq.upper()
    ok = set("ACGT")
    return {s[i:i + k] for i in range(len(s) - k + 1) if not (set(s[i:i + k]) - ok)}


# ------------------------------------------------------------- input parsing --
def load_classification(tsv_path):
    df = pd.read_csv(tsv_path, sep="\t", dtype=str).fillna("")
    return df


def list_genera(df):
    c = Counter(g for g in df["Genus"] if g.strip())
    for g, n in c.most_common():
        fam = df.loc[df["Genus"] == g, "Family"].mode()
        fam = fam.iloc[0] if len(fam) else "?"
        print(f"{g}\t{n}\t{fam}")
    log(f"total genera: {len(c)}")


def pick_user_contigs(df, contigs_fasta, genus, species=None,
                      min_len=800, max_n=5.0):
    """返回 {contig_id: record}，属级匹配不区分大小写。"""
    sub = df[df["Genus"].str.lower() == genus.lower()]
    if species:
        sub = sub[sub["Species"].str.lower() == species.lower()]
    ids = dict(zip(sub["Nucleotide"], sub["Genus_Confidence"]))
    got = {}
    for rec in SeqIO.parse(contigs_fasta, "fasta"):
        if rec.id in ids:
            seq = str(rec.seq).strip()
            if len(seq) >= min_len and n_percent(seq) <= max_n:
                got[rec.id] = SeqRecord(Seq(seq), id=rec.id, description="")
    dropped = len(ids) - len(got)
    if dropped:
        log(f"[filter] dropped {dropped} user contigs (length < {min_len} bp or N > {max_n}%).")
    return got, sub


# ------------------------------------------------------- reference sampling --
def sample_reference_accessions(taxa_df, genus, mode, species=None,
                                target_ratio=1.0, target_cap=None,
                                other_count=5):
    """virus.tree.py get_sampled_accessions 的矩阵适配版（不含跨科外群）。"""
    hits = taxa_df[taxa_df["Genus"].str.lower() == genus.lower()]
    if hits.empty:
        raise ValueError(f"Genus '{genus}' not found in taxa.txt")
    t_family = hits.iloc[0]["Family"]
    picked = {}

    if mode == "lineage":
        if not species:
            raise ValueError("--species required in lineage mode")
        grp = hits[hits["Species"].str.lower() == species.lower()]
        picked.update({a: (genus, sp) for a, sp in
                       zip(grp["Virus GENBANK accession"], grp["Species"])})
        log(f"[sample] lineage '{species}': kept ALL {len(grp)} reference strains.")
        return picked, t_family

    grp = hits.sample(frac=1.0, random_state=0)          # 打乱后截断，保证可复现
    if mode == "phylogeny":
        n_keep = max(1, int(len(grp) * target_ratio))
    else:                                                 # simple
        n_keep = len(grp)
    if target_cap:
        n_keep = min(n_keep, target_cap)
    for acc, sp in zip(grp.head(n_keep)["Virus GENBANK accession"],
                       grp.head(n_keep)["Species"]):
        picked[str(acc).split(";")[0].strip()] = (genus, sp)
    log(f"[sample] target genus '{genus}' ({mode}): {len(picked)} accessions.")

    if other_count > 0:
        fam_df = taxa_df[taxa_df["Family"] == t_family]
        for g, group in fam_df.groupby("Genus"):
            if g.lower() == genus.lower():
                continue
            take = group.sample(n=min(len(group), other_count), random_state=0)
            for acc, sp in zip(take["Virus GENBANK accession"], take["Species"]):
                picked[str(acc).split(";")[0].strip()] = (g, sp)
        n_other = len(picked) - sum(1 for v in picked.values() if v[0].lower() == genus.lower())
        log(f"[sample] other genera in family '{t_family}': +{n_other} (<= {other_count}/genus).")

    return picked, t_family


def load_ref_records(accessions, db_fasta):
    """accession（无版本）匹配 FASTA 头（常带版本号）；返回 [(rec, genus, species)]。"""
    want = {}
    for acc, (g, sp) in accessions.items():
        want[acc.split(".")[0]] = (acc, g, sp)
    idx = SeqIO.index(db_fasta, "fasta")
    alias = {}
    for k in idx.keys():
        base = k.split(".")[0]
        if base not in alias:
            alias[base] = k
    out = []
    missing = []
    for base, (acc, g, sp) in want.items():
        key = alias.get(base) or (base if base in idx else None)
        if key is None:
            missing.append(acc)
            continue
        rec = idx[key]
        rec.id = acc
        rec.description = ""
        out.append((rec, g, sp))
    idx.close()
    if missing:
        log(f"[warn] {len(missing)} accessions absent from {Path(db_fasta).name}: "
            f"{missing[:5]}{' ...' if len(missing) > 5 else ''}")
    return out


# ------------------------------------------------------------ auto-orientation --
def auto_orient(user_recs, ref_recs, k=K_ORIENT):
    """k-mer 共享计数定向：与最长 3 条同属参考比较正向/反向互补的共享量。"""
    reps = sorted(ref_recs, key=lambda r: len(r.seq), reverse=True)[:3]
    if not reps:
        return {r.id: False for r in user_recs.values()}
    ref_kmers = set()
    for r in reps:
        ref_kmers |= kmer_set(str(r.seq))
    flips = {}
    for rid, rec in user_recs.items():
        ks = kmer_set(str(rec.seq))
        if len(ks) < 5:
            flips[rid] = False
            continue
        fwd = len(ks & ref_kmers)
        rc = kmer_set(str(rec.seq.reverse_complement()))
        rev = len(rc & ref_kmers)
        flips[rid] = rev > fwd * 1.1
        if flips[rid]:
            rec.seq = rec.seq.reverse_complement()
            log(f"[orient] reversed {rid} (fwd={fwd}, rev={rev})")
    return flips


# ------------------------------------------------------------------- main --
def build_argparser():
    p = argparse.ArgumentParser(
        description="Genus-level whole-genome SDT matrix from ACVirus 09b results "
                    "(no reference GenBank needed)",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--analysis_dir", required=True,
                   help="09b_ACVirus_Analysis directory")
    p.add_argument("--genus", help="target genus (e.g. Cilevirus)")
    p.add_argument("--species", help="target species (lineage mode)")
    p.add_argument("--db", default=os.path.expanduser("~/database/virus-db/acvirus_db"),
                   help="acvirus_db directory")
    p.add_argument("-o", "--outdir", default=None)
    p.add_argument("--mode", choices=["simple", "phylogeny", "lineage"],
                   default="simple")
    p.add_argument("--target_ratio", type=float, default=1.0,
                   help="[phylogeny] fraction of target genus refs kept")
    p.add_argument("--target_cap", type=int, default=None,
                   help="hard cap on target genus refs (after ratio)")
    p.add_argument("--other_count", type=int, default=5,
                   help="max refs per OTHER genus in same family (0 disables)")
    p.add_argument("--min_length", type=int, default=800)
    p.add_argument("--max_n", type=float, default=5.0)
    p.add_argument("--no_orient", action="store_true",
                   help="disable k-mer auto-orientation of user contigs")
    p.add_argument("--short_labels", action="store_true",
                   help="compress contig ids to NODE_xxx style for plotting")
    p.add_argument("--aligner", choices=["mafft", "muscle", "clustalw"],
                   default="mafft")
    p.add_argument("--threads", type=int, default=8)
    p.add_argument("--resume", action="store_true")
    p.add_argument("--dpi", type=int, default=300)
    p.add_argument("--plot_format", choices=["png", "pdf"], default="png")
    p.add_argument("--palette", choices=["sdt", "cividis"], default="sdt",
                   help="heatmap colorscale (cividis = colorblind-safe)")
    p.add_argument("--list_genera", action="store_true",
                   help="list genera found in classification result and exit")
    return p


def short_id(name):
    import re
    m = re.search(r"NODE_\d+", name)
    return m.group(0) if m else name[:24]


def main():
    args = build_argparser().parse_args()
    adir = Path(args.analysis_dir)
    db = Path(args.db)
    if args.list_genera:
        df0 = load_classification(adir / "acvirus_classify" / "final_result_with_confidence.tsv")
        list_genera(df0)
        return
    if not args.outdir:
        sys.exit("[ERROR] -o/--outdir is required.")
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    vap.PLOT_DPI = args.dpi
    if args.palette == "cividis":
        vap.SDT_COLORS = [tuple(c[:3]) for c in plt.cm.cividis(np.linspace(0, 1, 10))]
        log("[style] cividis colorblind-safe palette enabled.")

    # ---------- 1. classification ----------
    cls_tsv = adir / "acvirus_classify" / "final_result_with_confidence.tsv"
    df = load_classification(cls_tsv)

    if not args.genus:
        sys.exit("[ERROR] --genus is required (or use --list_genera to inspect).")

    # ---------- 2. user contigs ----------
    contigs_fna = adir / "acvirus_classify" / "contigs.fna"
    user_recs, sub = pick_user_contigs(
        df, contigs_fna, args.genus, args.species,
        min_len=args.min_length, max_n=args.max_n)
    if not user_recs:
        sys.exit(f"[ERROR] no user contigs assigned to Genus '{args.genus}'. "
                 "Try --list_genera.")
    log(f"[input] {len(user_recs)} user contigs in genus '{args.genus}'.")

    # ---------- 3. reference sampling ----------
    taxa_df = pd.read_csv(db / "taxa.txt", dtype=str).fillna("")
    acc_map, t_family = sample_reference_accessions(
        taxa_df, args.genus, args.mode, args.species,
        target_ratio=args.target_ratio, target_cap=args.target_cap,
        other_count=args.other_count)
    ref_recs = load_ref_records(acc_map, db / "all_virus.fasta")
    log(f"[input] {len(ref_recs)} reference genomes loaded from db.")

    # ---------- 4. orientation ----------
    ref_only = [r for r, _, _ in ref_recs]
    flips = {} if args.no_orient else auto_orient(user_recs, ref_only)

    # ---------- 5. dedup + merge ----------
    seen = set()
    used_labels = {}
    records = []          # (label, record, source, genus, species, conf, flipped)
    flips_by_label = {}

    def uniq(label):
        if label not in used_labels:
            used_labels[label] = 1
            return label
        used_labels[label] += 1
        return f"{label}#{used_labels[label]}"

    def add(label, rec, source, g, sp, conf, flipped=False):
        h = sha_seq(str(rec.seq))
        if h in seen:
            log(f"[dedup] skip duplicate sequence: {label}")
            return
        seen.add(h)
        records.append((uniq(label), rec, source, g, sp, conf, flipped))

    for rid, rec in user_recs.items():
        row = sub[sub["Nucleotide"] == rid].iloc[0]
        label = short_id(rid) if args.short_labels else rid
        flips_by_label[label] = flips.get(rid, False)
        add(label, rec, "user_contig", args.genus,
            row["Species"], row["Genus_Confidence"], flips_by_label[label])

    for rec, g, sp in ref_recs:
        add(rec.id, rec, "reference", g, sp, "")

    man = pd.DataFrame([{
        "Label": lb, "Source": src, "Genus": g, "Species": sp,
        "Genus_Confidence": conf, "Length_bp": len(rec.seq),
        "N_pct": round(n_percent(str(rec.seq)), 2),
        "Reversed_for_alignment": flipped,
    } for lb, rec, src, g, sp, conf, flipped in records])
    man.to_csv(outdir / "inclusion_manifest.tsv", sep="\t", index=False)
    log(f"[manifest] {len(man)} sequences -> inclusion_manifest.tsv")

    if len(records) < 3:
        sys.exit("[ERROR] fewer than 3 unique sequences; matrix meaningless.")

    # ---------- 6. SDT matrix ----------
    names = [r[0] for r in records]
    seq_recs = [r[1] for r in records]
    cache_dir = outdir / "cache"
    cache_dir.mkdir(exist_ok=True)
    log(f"[matrix] running sdt_exact ({args.aligner}), {len(names)} sequences -> "
        f"{len(names)*(len(names)-1)//2} pairwise alignments ...")
    mat = vap.build_mat_sdt_exact(
        seq_recs, is_p=False, ign=False,
        aligner_path=args.aligner, aligner_type=args.aligner,
        threads=args.threads, desc="pairwise-identity",
        out_dir=cache_dir, cache_prefix=f"sdt_{args.genus}",
        resume=args.resume)

    # ---------- 7. CSV ----------
    n = len(names)
    full = np.full((n, n), np.nan)
    for i in range(n):
        full[i, i] = 100.0
        for j in range(n):
            if i != j:
                v = mat[i, j] if i < j else mat[j, i]
                full[i, j] = v
    with open(outdir / f"sdt_matrix_{args.genus}.csv", "w", newline="",
              encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow([""] + names)
        for i, nm in enumerate(names):
            w.writerow([nm] + ["NaN" if np.isnan(full[i, j]) else f"{full[i,j]:.1f}"
                               for j in range(n)])

    # ---------- 8. ordering + plots ----------
    ord_idx = vap.get_safe_leaf_order(full, n) if n >= 3 else list(range(n))
    names_o = [names[i] for i in ord_idx]
    mat_o = full[np.ix_(ord_idx, ord_idx)]

    ns = argparse.Namespace(plot_only_csv=False, no_nj=False,
                            align_method=args.aligner, plot=True,
                            plot_format=args.plot_format)
    title = (f"Genus {args.genus}: {sum(1 for r in records if r[2]=='user_contig')} "
             f"user contigs + {sum(1 for r in records if r[2]=='reference')} db refs "
             f"(09b ACVirus)")
    vap.execute_single_output(mat_o, names_o, outdir,
                              f"SDT_{args.genus}", ns, title, is_nt=True)

    log(f"[done] outputs in {outdir.resolve()}")


if __name__ == "__main__":
    main()
