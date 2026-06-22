#!/usr/bin/env python3
"""
build_ref_info.py — 生成 auto_known_virus.py 的 --ref_info TSV
=============================================================
从 suvtk taxonomy + 05_Taxonomy R 共识 + all_plant_viruses.fasta 生成

输入:
  09_Virome_Analysis/suvtk_taxonomy/taxonomy.tsv
  05_Taxonomy/Votus.integrated/final_integrated_classification.tsv
  08_Rescue/all_plant_viruses.fasta

输出:
  09_Virome_Analysis/ref_info.tsv

用法:
  python build_ref_info.py -o out/
"""

import argparse, os, sys
from pathlib import Path

def _read_tsv(path):
    rows = []
    p = Path(path)
    if not p.is_file(): return rows
    with open(p) as f:
        hdr = f.readline().strip().split('\t')
        for line in f:
            if not line.strip(): continue
            rows.append(dict(zip(hdr, line.strip().split('\t'))))
    return rows

def main():
    p = argparse.ArgumentParser(description="生成 ref_info TSV")
    p.add_argument("-o", "--output-dir", required=True, help="流水线输出根目录 (out/)")
    args = p.parse_args()
    root = Path(args.output_dir)

    analysis = root / "09_Virome_Analysis"

    # ── 1. 读取多源数据 ──
    # suvtk taxonomy
    suvtk_data = {}
    sv_tax = analysis / "suvtk_taxonomy" / "taxonomy.tsv"
    if sv_tax.is_file():
        for r in _read_tsv(sv_tax):
            cid = r.get("contig_id", r.get("seq_name",""))
            if cid:
                suvtk_data[cid] = {
                    "suvtk_species": r.get("Species", r.get("species","")),
                    "suvtk_genus": r.get("Genus", r.get("genus","")),
                    "suvtk_family": r.get("Family", r.get("family","")),
                    "suvtk_taxid": r.get("taxid", r.get("TaxID","")),
                }

    # R 共识 (05_Taxonomy)
    r_data = {}
    r_tsv = root / "05_Taxonomy" / "Votus.integrated" / "final_integrated_classification.tsv"
    if r_tsv.is_file():
        for r in _read_tsv(r_tsv):
            cid = r.get("contig_id","").strip('"')
            if cid:
                r_data[cid] = {
                    "r_species": r.get("Species","").strip('"'),
                    "r_genus": r.get("Genus","").strip('"'),
                    "r_family": r.get("Family","").strip('"'),
                    "r_realm": r.get("Realm","").strip('"'),
                    "primary_tool": r.get("primary_tool","").strip('"'),
                }

    # all_plant_viruses.fasta (序列长度)
    seq_lens = {}
    all_fa = root / "08_Rescue" / "all_plant_viruses.fasta"
    if all_fa.is_file():
        seq = ""
        for line in open(all_fa):
            if line.startswith('>'):
                if seq:
                    seq_lens[cid] = len(seq)
                cid = line[1:].split()[0]
                seq = ""
            else:
                seq += line.strip()
        if seq: seq_lens[cid] = len(seq)

    # suvtk features (CDS 统计)
    feat_data = {}
    sv_feat = analysis / "suvtk_features" / "featuretable.tbl"
    if sv_feat.is_file():
        cur = None
        for line in open(sv_feat):
            line = line.strip()
            if line.startswith(">Feature"):
                cur = line.split()[-1] if len(line.split()) > 1 else None
                if cur and cur not in feat_data:
                    feat_data[cur] = {"cds": 0, "trna": 0}
            elif cur:
                if "CDS" in line:
                    feat_data.setdefault(cur, {"cds": 0, "trna": 0})["cds"] += 1
                if "tRNA" in line:
                    feat_data.setdefault(cur, {"cds": 0, "trna": 0})["trna"] += 1

    # ── 2. 合并输出 ──
    analysis.mkdir(parents=True, exist_ok=True)
    out_tsv = analysis / "ref_info.tsv"

    all_ids = set(seq_lens.keys())
    if not all_ids:
        all_ids = set(suvtk_data.keys()) | set(r_data.keys())

    with open(out_tsv, "w") as f:
        cols = ["Accession", "Length", "Species", "Genus", "Family", "Realm",
                "suvtk_Species", "suvtk_Genus", "suvtk_Family",
                "CDS_Count", "tRNA_Count", "Primary_Tool"]
        f.write("\t".join(cols) + "\n")

        for cid in sorted(all_ids):
            r = r_data.get(cid, {})
            sv = suvtk_data.get(cid, {})
            ft = feat_data.get(cid, {})
            seq_len = seq_lens.get(cid, "")

            # 最佳 Species: R 共识 > suvtk
            best_sp = r.get("r_species","") or sv.get("suvtk_species","")
            best_ge = r.get("r_genus","") or sv.get("suvtk_genus","")
            best_fa = r.get("r_family","") or sv.get("suvtk_family","")

            vals = [
                cid, seq_len,
                best_sp, best_ge, best_fa, r.get("r_realm",""),
                sv.get("suvtk_species",""), sv.get("suvtk_genus",""), sv.get("suvtk_family",""),
                ft.get("cds", 0), ft.get("trna", 0),
                r.get("primary_tool",""),
            ]
            f.write("\t".join(str(v) for v in vals) + "\n")

    n = len(all_ids)
    print(f"ref_info.tsv: {n} 条 → {out_tsv}")
    print(f"  含 suvtk 分类: {len(suvtk_data)} 条")
    print(f"  含 R 共识:      {len(r_data)} 条")
    print(f"  含 CDS 统计:    {len(feat_data)} 条")
    print(f"  用法: --ref_info {out_tsv} --reference {all_fa}")


if __name__ == "__main__":
    main()
