#!/usr/bin/env python3
"""
integrated_summary.py — 双源交叉验证整合报告 (R共识 + suvtk)
============================================================
输入:
  05_Taxonomy/Votus.integrated/final_integrated_classification.tsv  (R 共识)
  09_Virome_Analysis/suvtk_taxonomy/taxonomy.tsv                     (suvtk)
  09_Virome_Analysis/suvtk_features/featuretable.tbl                 (CDS)
  08_Rescue/all_plant_viruses.fasta                                  (序列)
输出:
  09_Virome_Analysis/integrated_summary.tsv
"""

import argparse, os, sys
from pathlib import Path
from collections import defaultdict

def _read_tsv(path):
    rows = []
    if not Path(path).is_file(): return rows
    with open(path) as f:
        hdr = f.readline().strip().split('\t')
        for line in f:
            if not line.strip(): continue
            rows.append(dict(zip(hdr, line.strip().split('\t'))))
    return rows

def _count_fasta(path):
    if not Path(path).is_file(): return 0
    return sum(1 for _ in open(path) if _.startswith('>'))

def load_r_consensus(path):
    """05_Taxonomy R 共识: contig_id → {rank: value, tools, agrees}"""
    data = {}
    if not Path(path).is_file(): return data
    rows = _read_tsv(path)
    for r in rows:
        cid = r.get("contig_id","").strip().strip('"')
        if not cid: continue
        data[cid] = {
            "Realm": r.get("Realm","").strip('"'),
            "Kingdom": r.get("Kingdom","").strip('"'),
            "Phylum": r.get("Phylum","").strip('"'),
            "Class": r.get("Class","").strip('"'),
            "Order": r.get("Order","").strip('"'),
            "Family": r.get("Family","").strip('"'),
            "Genus": r.get("Genus","").strip('"'),
            "Species": r.get("Species","").strip('"'),
            "n_tools": int(r.get("completeness",0) or 0),
            "confidence": float(r.get("confidence",1) or 1),
            "primary_tool": r.get("primary_tool","").strip('"'),
            "species_agree": r.get("Species_agree","").strip('"') or "",
            "genus_agree": r.get("Genus_agree","").strip('"') or "",
        }
    return data

def load_suvtk_taxonomy(path):
    """suvtk taxonomy.tsv: 实际格式为 contig\ttaxonomy (最低可定级别 '<Taxon> sp.')"""
    data = {}
    if not Path(path).is_file(): return data
    rows = _read_tsv(path)
    for r in rows:
        cid = r.get("contig", r.get("contig_id", r.get("seq_name", "")))
        if not cid: continue
        taxonomy = r.get("taxonomy", "").strip()
        first = taxonomy.split()[0] if taxonomy else ""
        sv_ge = sv_fa = ""
        if first.endswith("idae"):        # 科
            sv_fa = first
        elif first.endswith("virinae"):   # 亚科
            pass
        elif first.endswith("virus"):     # 属
            sv_ge = first
        data[cid] = {
            "suvtk_family": sv_fa,
            "suvtk_genus": sv_ge,
            "suvtk_species": taxonomy,
        }
    return data


def load_miuvig_taxonomy(path):
    """suvtk miuvig_taxonomy.tsv: contig\pred_genome_type\pred_genome_struc"""
    data = {}
    if not Path(path).is_file(): return data
    rows = _read_tsv(path)
    for r in rows:
        cid = r.get("contig", r.get("contig_id", ""))
        if not cid: continue
        data[cid] = {
            "genome_type": r.get("pred_genome_type", ""),
            "genome_struc": r.get("pred_genome_struc", ""),
        }
    return data

def load_suvtk_features(path):
    """suvtk featuretable.tbl: 统计每 contig 的 CDS/tRNA 数量"""
    data = defaultdict(lambda: {"cds_count": 0, "trna_count": 0, "gene_products": []})
    if not Path(path).is_file(): return data
    current_contig = None
    for line in open(path):
        line = line.strip()
        if line.startswith(">Feature"):
            current_contig = line.split()[-1] if len(line.split()) > 1 else None
        elif "CDS" in line and current_contig:
            data[current_contig]["cds_count"] += 1
            # 提取 product 注释
            if "product" in line.lower():
                parts = line.split("\t")
                for p in parts:
                    if p.lower().startswith("product"):
                        data[current_contig]["gene_products"].append(p.split("product")[-1].strip())
        elif "tRNA" in line and current_contig:
            data[current_contig]["trna_count"] += 1
    return data

def compare_taxonomy(r_cons, suvtk, cid):
    """双源分类比较 (R共识 + suvtk): 取最佳"""
    result = {
        "tax_source": "R_consensus",
        "best_species": "", "best_genus": "", "best_family": "",
        "r_species": "", "suvtk_species": "",
        "r_genus": "", "suvtk_genus": "",
        "tax_consensus": 0,
    }
    r = r_cons.get(cid, {})
    sv = suvtk.get(cid, {})

    r_sp = r.get("Species",""); r_ge = r.get("Genus",""); r_fa = r.get("Family","")
    result["r_species"] = r_sp; result["r_genus"] = r_ge; result["r_family"] = r_fa

    sv_sp = sv.get("suvtk_species",""); sv_ge = sv.get("suvtk_genus",""); sv_fa = sv.get("suvtk_family","")
    result["suvtk_species"] = sv_sp; result["suvtk_genus"] = sv_ge

    species_set = {x for x in [r_sp, sv_sp] if x and x != "NA"}
    genus_set = {x for x in [r_ge, sv_ge] if x and x != "NA"}
    result["tax_consensus"] = 2 if len(species_set) == 1 else (1 if len(genus_set) == 1 else 0)

    result["best_family"] = r_fa or sv_fa or ""
    result["best_genus"] = r_ge or sv_ge or ""
    result["best_species"] = r_sp or sv_sp or ""

    for key in ["suvtk_family","suvtk_genus","suvtk_species"]:
        if key not in result: result[key] = ""
    return result


def main():
    p = argparse.ArgumentParser(description="三源交叉验证整合报告")
    p.add_argument("--output-dir", "-o", required=True, help="输出根目录 (通常为 out/)")
    p.add_argument("--analysis-dir", default=None, help="09_Virome_Analysis 目录 (默认: {output_dir}/09_Virome_Analysis)")
    args = p.parse_args()

    root = Path(args.output_dir)
    analysis_dir = Path(args.analysis_dir) if args.analysis_dir else root / "09_Virome_Analysis"

    # 输入 — 自动探测 {sample}.integrated 目录
    r_tsv = root / "05_Taxonomy" / "Votus.integrated" / "final_integrated_classification.tsv"
    if not r_tsv.is_file():
        tax_dir = root / "05_Taxonomy"
        for d in tax_dir.glob("*.integrated"):
            candidate = d / "final_integrated_classification.tsv"
            if candidate.is_file():
                r_tsv = candidate
                break
        else:
            d2 = tax_dir / "integrated"
            c2 = d2 / "final_integrated_classification.tsv"
            if c2.is_file():
                r_tsv = c2
    sv_tax = analysis_dir / "suvtk_taxonomy" / "taxonomy.tsv"
    sv_feat = analysis_dir / "suvtk_features" / "featuretable.tbl"
    all_fa = root / "08_Rescue" / "all_plant_viruses.fasta"

    print("=== 双源交叉验证整合 (R共识 + suvtk) ===\n")

    print(f"[1] 05_Taxonomy R 共识: {r_tsv}")
    r_cons = load_r_consensus(r_tsv)
    print(f"    {len(r_cons)} 条分类记录")

    print(f"[2] suvtk taxonomy: {sv_tax}")
    suvtk_tax = load_suvtk_taxonomy(sv_tax)
    print(f"    {len(suvtk_tax)} 条分类记录")

    print(f"[3] suvtk features: {sv_feat}")
    features = load_suvtk_features(sv_feat)
    print(f"    {len(features)} 条 CDS 特征")

    miuvig_tax_path = analysis_dir / "suvtk_taxonomy" / "miuvig_taxonomy.tsv"
    print(f"[3b] suvtk miuvig taxonomy: {miuvig_tax_path}")
    miuvig_data = load_miuvig_taxonomy(miuvig_tax_path)
    print(f"    {len(miuvig_data)} 条 genome_type 记录")

    n_plant = _count_fasta(all_fa)
    print(f"[4] 植物病毒序列: {all_fa}")
    print(f"    {n_plant} 条序列\n")

    # 4b. CheckV completeness
    checkv_data = {}
    for cv_tsv in [
        root / "08_Rescue" / "checkv" / "Plant" / "completeness.tsv",
        root / "08_Rescue" / "checkv" / "no_rescue" / "completeness.tsv",
        root / "07_Checkv" / "Plant" / "completeness.tsv",
    ]:
        if cv_tsv.is_file():
            for r in _read_tsv(cv_tsv):
                cid = r.get("contig_id", "")
                if cid and cid not in checkv_data:
                    checkv_data[cid] = {
                        "aai_completeness": r.get("aai_completeness", ""),
                        "aai_confidence": r.get("aai_confidence", ""),
                        "aai_expected_length": r.get("aai_expected_length", ""),
                        "aai_error": r.get("aai_error", ""),
                        "aai_num_hits": r.get("aai_num_hits", ""),
                        "aai_top_hit": r.get("aai_top_hit", ""),
                        "aai_id": r.get("aai_id", ""),
                        "aai_af": r.get("aai_af", ""),
                    }
    print(f"[4b] CheckV: {len(checkv_data)} 条评估记录")

    # 整合 — 同时读取序列长度
    plant_ids = set()
    seq_lens = {}
    if all_fa.is_file():
        cur_id, cur_seq = None, []
        for line in open(all_fa):
            if line.startswith('>'):
                if cur_id:
                    seq_lens[cur_id] = len("".join(cur_seq))
                cur_id = line[1:].split()[0]
                plant_ids.add(cur_id)
                cur_seq = []
            else:
                cur_seq.append(line.strip())
        if cur_id:
            seq_lens[cur_id] = len("".join(cur_seq))

    cols = [
        "contig_id", "length",
        "cds_count", "trna_count",
        "genome_type", "genome_struc",
        "aai_completeness", "aai_confidence", "aai_expected_length",
        "aai_error", "aai_num_hits", "aai_top_hit", "aai_id", "aai_af",
        "tax_consensus",
        "best_family", "best_genus", "best_species",
        "r_family", "r_genus", "r_species",
        "suvtk_family", "suvtk_genus", "suvtk_species",
    ]

    analysis_dir.mkdir(parents=True, exist_ok=True)
    out_tsv = analysis_dir / "integrated_summary.tsv"
    with open(out_tsv, "w") as of:
        of.write("\t".join(cols) + "\n")
        for cid in sorted(plant_ids):
            r = r_cons.get(cid, {})
            sv = suvtk_tax.get(cid, {})
            ft = features.get(cid, {})

            tax = compare_taxonomy(r_cons, suvtk_tax, cid)
            miu = miuvig_data.get(cid, {})
            cv = checkv_data.get(cid, {})

            vals = [
                cid, seq_lens.get(cid, ""),
                ft.get("cds_count", 0),
                ft.get("trna_count", 0),
                miu.get("genome_type", ""),
                miu.get("genome_struc", ""),
                cv.get("aai_completeness", ""),
                cv.get("aai_confidence", ""),
                cv.get("aai_expected_length", ""),
                cv.get("aai_error", ""),
                cv.get("aai_num_hits", ""),
                cv.get("aai_top_hit", ""),
                cv.get("aai_id", ""),
                cv.get("aai_af", ""),
                tax["tax_consensus"],
                tax["best_family"], tax["best_genus"], tax["best_species"],
                tax["r_family"], tax["r_genus"], tax["r_species"],
                tax["suvtk_family"], tax["suvtk_genus"], tax["suvtk_species"],
            ]
            of.write("\t".join(str(v) for v in vals) + "\n")

    print(f"整合完成: {out_tsv}")
    print(f"  共 {len(plant_ids)} 条植物病毒")

    # 统计
    n_cds_total = sum(v.get("cds_count",0) for v in features.values())
    print(f"  共 {len(plant_ids)} 条植物病毒, {n_cds_total} CDS")


if __name__ == "__main__":
    main()
