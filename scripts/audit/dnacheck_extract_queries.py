#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""枸杞 DNA 病毒真伪判定 · 第 1 步
从 plant_virus_DNA_RNA_final.tsv 取 DNA 类记录（植物病毒 + 植物病毒?），
按队列从 03b_MergeSamples 的合并 fasta 抽序列，产出 query.fasta + 清单。
同时 join prevalence_full_table.tsv 的 n_samples/prevalence。
只读原表，产物写 ~/goji_dnacheck_20260929/。
"""
import csv, os

FINAL = "/home/zhangwenda/MMPV-paper/goji-virome/00_scripts/_mining_rounds_20260914/outputs/plant_virus_DNA_RNA_final.tsv"
BASE = "/home/zhangwenda/MMPV-paper/goji-virome/02_novel_virus"
OUT = "/home/zhangwenda/goji_dnacheck_20260929"

COHORT_DIR = {
    "Lycium_barbarum": "RNA-Lycium_barbarum_out",
    "Lycium_chinense": "RNA-Lycium_chinense_out",
    "Lycium_ruthenicum": "RNA-Lycium_ruthenicum_out",
    "Lycium_amarum": "RNA-Lycium_amarum_out",
    "Aphis_gossypii": "RNA-Aphis_gossypii_out",
    "Fusarium_nematophilum": "RNA-Fusarium_nematophilum_out",
    "Alternaria_alternata": "RNA-Alternaria_alternata_out",
}


def read_fasta(path):
    seqs = {}
    name = None
    buf = []
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buf)
                name = line[1:].split()[0]
                buf = []
            else:
                buf.append(line.strip())
        if name is not None:
            seqs[name] = "".join(buf)
    return seqs


def main():
    os.makedirs(OUT, exist_ok=True)
    recs = []
    with open(FINAL) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r["Nucleic_acid"] == "DNA" and r["Category"] in ("植物病毒", "植物病毒?"):
                recs.append(r)
    print("DNA 记录:", len(recs))

    by_cohort = {}
    for r in recs:
        by_cohort.setdefault(r["Cohort"], []).append(r)

    # prevalence 表 join
    prev = {}
    for co, d in COHORT_DIR.items():
        fp = f"{BASE}/{d}/prevalence_full_table.tsv"
        if not os.path.isfile(fp):
            continue
        with open(fp) as f:
            for r in csv.DictReader(f, delimiter="\t"):
                prev[(co, r["contig_id"])] = (r.get("n_samples") or "", r.get("prevalence_%") or "",
                                              r.get("checkv_completeness") or "")

    got, missing = [], []
    for co, items in by_cohort.items():
        d = COHORT_DIR[co]
        fa = f"{BASE}/{d}/03b_MergeSamples/all_sample_virus.fasta"
        alt = f"{BASE}/{d}/03b_MergeSamples/all_sample_virus_combined.fasta"
        src = fa if os.path.isfile(fa) else alt
        if not os.path.isfile(src):
            print(f"[warn] {co}: 无合并 fasta")
            missing.extend(items)
            continue
        seqs = read_fasta(src)
        print(f"{co}: 合并 fasta {len(seqs)} 条")
        for it in items:
            cid = it["contig_id"]
            if cid in seqs:
                got.append((co, it, seqs[cid]))
            else:
                missing.append(it)

    print("抽到序列:", len(got), " 未命中:", len(missing))

    with open(f"{OUT}/query.fasta", "w") as fo:
        for co, it, s in got:
            fo.write(f">{co}__{it['contig_id']}\n{s}\n")

    with open(f"{OUT}/dna_manifest.tsv", "w") as fo:
        fo.write("Cohort\tcontig_id\tGenus\tSpecies\tFamily_true\tCategory\tn_samples\tprevalence_pct\tcheckv\tlength\n")
        for co, it, s in got:
            n, p, ck = prev.get((co, it["contig_id"]), ("", "", ""))
            fo.write(f"{co}\t{it['contig_id']}\t{it['Genus']}\t{it['Species']}\t{it['Family_true']}\t"
                     f"{it['Category']}\t{n}\t{p}\t{ck}\t{len(s)}\n")
        for it in missing:
            fo.write(f"{it['Cohort']}\t{it['contig_id']}\t{it['Genus']}\t{it['Species']}\t{it['Family_true']}\t"
                     f"{it['Category']}\t\t\t\tMISSING\n")

    # 未命中清单单独存
    with open(f"{OUT}/missing_ids.tsv", "w") as fo:
        fo.write("Cohort\tcontig_id\n")
        for it in missing:
            fo.write(f"{it['Cohort']}\t{it['contig_id']}\n")
    print("done ->", OUT)


if __name__ == "__main__":
    main()
