#!/usr/bin/env python3
"""
integrated_summary.py — 三源交叉验证整合报告
============================================
输入:
  05_Taxonomy/Votus.integrated/final_integrated_classification.tsv  (R 共识)
  09_Virome_Analysis/cenote_taker3/run_summary.tsv                  (Cenote-T3)
  09_Virome_Analysis/suvtk_taxonomy/taxonomy.tsv                     (suvtk)
  09_Virome_Analysis/suvtk_features/featuretable.tbl                 (CDS)
  08_Rescue/all_plant_viruses.fasta                                  (序列)
输出:
  09_Virome_Analysis/integrated_summary.tsv
  09_Virome_Analysis/integrated_taxonomy_comparison.tsv
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

def load_cenote_summary(path):
    """Cenote-T3 run_summary.tsv"""
    data = {}
    if not Path(path).is_file(): return data
    rows = _read_tsv(path)
    for r in rows:
        cid = r.get("contig_name", r.get("Contig",""))
        if not cid: continue
        # 标准化 contig_id: 取第一个空格前的部分
        cid = cid.split()[0]
        data[cid] = {
            "hallmark_count": int(r.get("hallmark_count", r.get("hallmark gene count",0)) or 0),
            "circular": str(r.get("circular", "")).lower() in ("true","1","yes"),
            "dtr": str(r.get("dtr", r.get("DTR", "")).strip()),
            "ct3_family": r.get("family", r.get("hallmark_family", "")),
            "ct3_genus": r.get("genus", r.get("hallmark_genus", "")),
            "ct3_taxonomy": r.get("taxonomy", r.get("hallmark_taxonomy", "")),
        }
    return data

def load_suvtk_taxonomy(path):
    """suvtk taxonomy.tsv: contig_id → {rank: value}"""
    data = {}
    if not Path(path).is_file(): return data
    rows = _read_tsv(path)
    for r in rows:
        cid = r.get("contig_id", r.get("seq_name", ""))
        if not cid: continue
        data[cid] = {
            "suvtk_realm": r.get("Realm", r.get("realm", "")),
            "suvtk_kingdom": r.get("Kingdom", r.get("kingdom", "")),
            "suvtk_phylum": r.get("Phylum", r.get("phylum", "")),
            "suvtk_class": r.get("Class", r.get("class", "")),
            "suvtk_order": r.get("Order", r.get("order", "")),
            "suvtk_family": r.get("Family", r.get("family", "")),
            "suvtk_genus": r.get("Genus", r.get("genus", "")),
            "suvtk_species": r.get("Species", r.get("species", "")),
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

def compare_taxonomy(r_cons, cenote, suvtk, cid):
    """三源分类比较: 取最佳"""
    result = {
        "tax_source": "R_consensus",  # 默认 R 共识
        "best_species": "", "best_genus": "", "best_family": "",
        "r_species": "", "ct3_species": "", "suvtk_species": "",
        "r_genus": "", "ct3_genus": "", "suvtk_genus": "",
        "tax_consensus": 0,  # 0-3 sources agree
    }
    r = r_cons.get(cid, {})
    ct = cenote.get(cid, {})
    sv = suvtk.get(cid, {})

    # R 共识 (最高优先级)
    r_sp = r.get("Species","")
    r_ge = r.get("Genus","")
    r_fa = r.get("Family","")
    result["r_species"] = r_sp
    result["r_genus"] = r_ge
    result["r_family"] = r_fa

    # Cenote-T3
    ct_sp = ct.get("ct3_taxonomy","").split(";")[-1] if ct.get("ct3_taxonomy") else ""
    ct_ge = ct.get("ct3_genus","")
    ct_fa = ct.get("ct3_family","")
    result["ct3_species"] = ct_sp
    result["ct3_genus"] = ct_ge

    # suvtk
    sv_sp = sv.get("suvtk_species","")
    sv_ge = sv.get("suvtk_genus","")
    sv_fa = sv.get("suvtk_family","")
    result["suvtk_species"] = sv_sp
    result["suvtk_genus"] = sv_ge

    # 一致性计数
    species_set = {x for x in [r_sp, ct_sp, sv_sp] if x and x != "NA"}
    genus_set = {x for x in [r_ge, ct_ge, sv_ge] if x and x != "NA"}
    result["tax_consensus"] = 1 if len(species_set) == 1 else (2 if len(genus_set) == 1 else 0)

    # 选择最佳: R 共识 > suvtk > Cenote
    result["best_family"] = r_fa or sv_fa or ct_fa or ""
    result["best_genus"] = r_ge or sv_ge or ct_ge or ""
    result["best_species"] = r_sp or sv_sp or ct_sp or ""

    return result


def main():
    p = argparse.ArgumentParser(description="三源交叉验证整合报告")
    p.add_argument("--output-dir", "-o", required=True, help="输出根目录 (通常为 out/)")
    p.add_argument("--analysis-dir", default=None, help="09_Virome_Analysis 目录 (默认: {output_dir}/09_Virome_Analysis)")
    args = p.parse_args()

    root = Path(args.output_dir)
    analysis_dir = Path(args.analysis_dir) if args.analysis_dir else root / "09_Virome_Analysis"

    # 输入
    r_tsv = root / "05_Taxonomy" / "Votus.integrated" / "final_integrated_classification.tsv"
    ct3_tsv = analysis_dir / "cenote_taker3" / "run_summary.tsv"
    sv_tax = analysis_dir / "suvtk_taxonomy" / "taxonomy.tsv"
    sv_feat = analysis_dir / "suvtk_features" / "featuretable.tbl"
    all_fa = root / "08_Rescue" / "all_plant_viruses.fasta"

    print("=== 三源交叉验证整合 ===\n")

    # 加载数据
    print(f"[1] 05_Taxonomy R 共识: {r_tsv}")
    r_cons = load_r_consensus(r_tsv)
    print(f"    {len(r_cons)} 条分类记录")

    print(f"[2] Cenote-Taker3: {ct3_tsv}")
    cenote = load_cenote_summary(ct3_tsv)
    print(f"    {len(cenote)} 条注释记录")

    print(f"[3] suvtk taxonomy: {sv_tax}")
    suvtk_tax = load_suvtk_taxonomy(sv_tax)
    print(f"    {len(suvtk_tax)} 条分类记录")

    print(f"[4] suvtk features: {sv_feat}")
    features = load_suvtk_features(sv_feat)
    print(f"    {len(features)} 条 CDS 特征")

    n_plant = _count_fasta(all_fa)
    print(f"[5] 植物病毒序列: {all_fa}")
    print(f"    {n_plant} 条序列\n")

    # 整合
    plant_ids = set()
    if all_fa.is_file():
        for line in open(all_fa):
            if line.startswith('>'):
                plant_ids.add(line[1:].split()[0])

    cols = [
        "contig_id", "length",
        "hallmark_count", "circular", "dtr",
        "cds_count", "trna_count",
        "tax_source", "tax_consensus",
        "best_family", "best_genus", "best_species",
        "r_family", "r_genus", "r_species",
        "ct3_family", "ct3_genus", "ct3_species",
        "suvtk_family", "suvtk_genus", "suvtk_species",
    ]

    out_tsv = analysis_dir / "integrated_summary.tsv"
    with open(out_tsv, "w") as of:
        of.write("\t".join(cols) + "\n")
        for cid in sorted(plant_ids):
            r = r_cons.get(cid, {})
            ct = cenote.get(cid, {})
            sv = suvtk_tax.get(cid, {})
            ft = features.get(cid, {})

            tax = compare_taxonomy(r_cons, cenote, suvtk_tax, cid)

            # 从 FASTA 获取长度
            seq_len = ""
            # (长度可选, 暂不读)
            vals = [
                cid, seq_len,
                ct.get("hallmark_count", 0),
                "Y" if ct.get("circular") else "N",
                ct.get("dtr", ""),
                ft.get("cds_count", 0),
                ft.get("trna_count", 0),
                tax["tax_consensus"],
                "R_consensus" if tax["best_species"] else ("suvtk" if tax["suvtk_species"] else "Cenote-T3"),
                tax["best_family"], tax["best_genus"], tax["best_species"],
                tax["r_family"], tax["r_genus"], tax["r_species"],
                tax["ct3_family"], tax["ct3_genus"], tax["ct3_species"],
                tax["suvtk_family"], tax["suvtk_genus"], tax["suvtk_species"],
            ]
            of.write("\t".join(str(v) for v in vals) + "\n")

    print(f"整合完成: {out_tsv}")
    print(f"  共 {len(plant_ids)} 条植物病毒")

    # 统计
    n_with_hallmark = sum(1 for v in cenote.values() if v.get("hallmark_count",0) > 0)
    n_circular = sum(1 for v in cenote.values() if v.get("circular"))
    n_cds_total = sum(v.get("cds_count",0) for v in features.values())


if __name__ == "__main__":
    main()
