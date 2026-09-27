#!/usr/bin/env python3
"""
build_ref_info.py — 生成 auto_known_virus.py 的 --ref_info TSV
=============================================================
从 suvtk taxonomy + 05_Taxonomy R 共识 + HQ_plant_viruses.fasta 生成

输入:
  09_Virome_Analysis/suvtk_taxonomy/taxonomy.tsv
  05_Taxonomy/Votus.integrated/final_integrated_classification.tsv
  08_Rescue/HQ_plant_viruses.fasta

输出:
  09_Virome_Analysis/ref_info.tsv

用法:
  python build_ref_info.py -o out/
"""

import argparse, os, sys
from pathlib import Path

# 跨管线统一 I/O 布局 (mmpv_common/, 仓库根): 目录名随 MMPV_IO_LAYOUT 解析
# (编排器已 normalize 环境变量, 子进程导入本模块时快照即正确布局)
_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)
from mmpv_common.io_layout import layout_dirs as _layout_dirs
_D = _layout_dirs(os.environ.get("MMPV_IO_LAYOUT", "legacy"))

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

    analysis = root / _D['d_analysis']
    # 09_ 是 NN[ab]_Name 命名里唯一的例外；改名 09a_Virome_Analysis 后两个名字都认
    if not analysis.is_dir() and (root / _D['d_analysis']).is_dir():
        analysis = root / _D['d_analysis']

    # ── 1. 读取多源数据 ──
    # suvtk taxonomy (多路径查找, 格式: contig\ttaxonomy)
    suvtk_data = {}
    sv_tax = analysis / "suvtk_taxonomy" / "taxonomy.tsv"
    if not sv_tax.is_file():
        sv_tax = root / _D['d_rescue'] / "suvtk.taxonomy_output" / "taxonomy.tsv"
    if sv_tax.is_file():
        for r in _read_tsv(sv_tax):
            cid = r.get("contig_id", r.get("contig", r.get("seq_name","")))
            taxonomy = r.get("taxonomy","").strip()
            # suvtk 输出最低可定级别为 "<Taxon> sp.", 按后缀推断 rank
            first = taxonomy.split()[0] if taxonomy else ""
            sv_ge = sv_fa = ""
            if first.endswith("idae"):        # 科
                sv_fa = first
            elif first.endswith("virinae"):   # 亚科 (无干净科/属)
                pass
            elif first.endswith("virus"):     # 属
                sv_ge = first
            # suvtk_species 保留完整 taxonomy 字符串供双源比对
            if cid:
                suvtk_data[cid] = {
                    "suvtk_species": taxonomy,
                    "suvtk_genus": sv_ge,
                    "suvtk_family": sv_fa,
                }

    # R 共识 (05_Taxonomy)
    r_data = {}
    r_tsv = root / _D['d_taxonomy'] / "Votus.integrated" / "final_integrated_classification.tsv"
    if not r_tsv.is_file():
        tax_dir = root / _D['d_taxonomy']
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

    # HQ_plant_viruses.fasta (序列长度)
    seq_lens = {}
    all_fa = root / _D['d_rescue'] / "HQ_plant_viruses.fasta"
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

    # suvtk features (多路径)
    feat_data = {}
    sv_feat = analysis / "suvtk_features" / "featuretable.tbl"
    if not sv_feat.is_file():
        sv_feat = root / _D['d_rescue'] / "suvtk.features_output" / "featuretable.tbl"
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
