#!/usr/bin/env python3
"""
panvirome_aggregator.py — 泛病毒组突变特征聚合引擎
====================================================
从 S3 (Virus_variants_Results) 中提取所有病毒的 iVar + SnpEff + SNPGenie
结果，合并为统一矩阵，供 panvirome_landscape.py 使用。

输入:
  --variants-dir  病毒变异结果根目录 (含 virus-variants/, virus-SnpEff/, virus-SNPGenie/, virus-annotations/)
输出:
  --outdir        每病毒一个子目录: merged_allele_freq.csv, merged_snpeff.maf, merged_site_pi.csv, gene_coordinates.csv
"""

import argparse
import os
import sys
import glob
import subprocess
import pandas as pd
from Bio import SeqIO


def parse_args():
    p = argparse.ArgumentParser(description="泛病毒组突变特征聚合")
    p.add_argument("--variants-dir", required=True,
                   help="病毒变异结果根目录 (含 virus-variants/, virus-SnpEff/, virus-SNPGenie/, virus-annotations/)")
    p.add_argument("--outdir", default="./Merged_Results_All_Viruses",
                   help="输出目录 (默认: ./Merged_Results_All_Viruses)")
    p.add_argument("--snpeff2maf", default=None,
                   help="snpeff2maf.py 路径 (自动检测: MMPV-RNA/virome_analysis_pipeline/snpeff2maf.py)")
    return p.parse_args()


def extract_sample_id(filename, acc_id):
    base = os.path.basename(filename)
    return base.split(acc_id)[0].strip("._")


def main():
    args = parse_args()
    base = args.variants_dir
    ann_dir = os.path.join(base, "virus-annotations")
    var_dir = os.path.join(base, "virus-variants")
    snpeff_dir = os.path.join(base, "virus-SnpEff")
    snpgenie_dir = os.path.join(base, "virus-SNPGenie")

    # 自动检测 snpeff2maf
    snpeff2maf = args.snpeff2maf
    if not snpeff2maf:
        candidates = [
            os.path.join(os.path.dirname(__file__), "snpeff2maf.py"),
            os.path.expanduser("~/bin/snpeff2maf.py"),
        ]
        for c in candidates:
            if os.path.isfile(c):
                snpeff2maf = c
                break

    gb_files = glob.glob(f"{ann_dir}/*.gb")
    print(f"  检测到 {len(gb_files)} 个 GenBank 注释文件")
    if not gb_files:
        print(f"  [ERROR] {ann_dir} 中无 .gb 文件")
        sys.exit(1)

    os.makedirs(args.outdir, exist_ok=True)

    for gb_file in gb_files:
        acc_id = os.path.basename(gb_file).replace(".gb", "")
        matched = [d for d in os.listdir(var_dir) if acc_id in d]
        if not matched:
            print(f"  [SKIP] {acc_id}: 无变异结果")
            continue

        virus_name = matched[0]
        out_dir = os.path.join(args.outdir, virus_name)
        os.makedirs(out_dir, exist_ok=True)
        print(f"\n  ▶ {virus_name}")

        # 1. 解析 GenBank → 基因坐标
        genes_data = []
        is_viroid = True
        genome_length = 0
        try:
            for rec in SeqIO.parse(gb_file, "genbank"):
                genome_length = len(rec.seq)
                for feat in rec.features:
                    if feat.type == "CDS":
                        is_viroid = False
                        gene = feat.qualifiers.get("gene", [""])[0] or feat.qualifiers.get("product", ["Unk"])[0]
                        genes_data.append({"gene": gene, "start": int(feat.location.start), "end": int(feat.location.end)})
            pd.DataFrame(genes_data).to_csv(os.path.join(out_dir, "gene_coordinates.csv"), index=False)
            tag = "类病毒/非编码" if is_viroid else f"{len(genes_data)} CDS"
            print(f"    GenBank: {genome_length}bp, {tag}")
        except Exception as e:
            print(f"    [ERROR] GenBank 解析: {e}")
            continue

        # 2. 合并 iVar allele frequencies
        ivar_files = glob.glob(f"{var_dir}/{virus_name}/**/*allele_frequencies.tsv", recursive=True)
        df_ivar = []
        for f in ivar_files:
            sid = extract_sample_id(f, acc_id)
            try:
                df = pd.read_csv(f, sep="\t")
                df = df[df["ALT_FREQ"] > 0.05]
                if not df.empty:
                    df["Sample_ID"] = sid
                    df_ivar.append(df)
            except Exception:
                pass
        if df_ivar:
            pd.concat(df_ivar).to_csv(os.path.join(out_dir, "merged_allele_freq.csv"), index=False)
            print(f"    iVar: {len(df_ivar)} 样本合并")

        # 3. 合并 SnpEff (MAF)
        if not is_viroid:
            target_dir = os.path.join(snpeff_dir, virus_name)
            maf_files = glob.glob(f"{target_dir}/**/*.maf", recursive=True)
            if not maf_files:
                vcf_files = glob.glob(f"{target_dir}/**/*.vcf*", recursive=True)
                if vcf_files and snpeff2maf and os.path.isfile(snpeff2maf):
                    print(f"    VCF→MAF 转换 ({len(vcf_files)} 文件)...")
                    subprocess.run([sys.executable, snpeff2maf, target_dir],
                                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
                    maf_files = glob.glob(f"{target_dir}/**/*.maf", recursive=True)
            df_maf = []
            for f in maf_files:
                try:
                    df = pd.read_csv(f, sep="\t", comment="#")
                    if not df.empty:
                        df_maf.append(df)
                except Exception:
                    pass
            if df_maf:
                pd.concat(df_maf).to_csv(os.path.join(out_dir, "merged_snpeff.maf"), sep="\t", index=False)
                print(f"    SnpEff MAF: {len(df_maf)} 文件合并")

        # 4. 合并 SNPGenie product_results
        if not is_viroid:
            prod_files = glob.glob(f"{snpgenie_dir}/{virus_name}/**/product_results.txt", recursive=True)
            df_prod = []
            for f in prod_files:
                sid = os.path.basename(os.path.dirname(f)).replace(f"_{acc_id}", "")
                try:
                    df = pd.read_csv(f, sep="\t")
                    if not df.empty:
                        df["Sample_ID"] = sid
                        df_prod.append(df)
                except Exception:
                    pass
            if df_prod:
                pd.concat(df_prod).to_csv(os.path.join(out_dir, "merged_product_results.csv"), index=False)
                print(f"    SNPGenie product: {len(df_prod)} 样本合并")

        # 5. 合并 site_results (π)
        site_files = glob.glob(f"{snpgenie_dir}/{virus_name}/**/site_results.txt", recursive=True)
        df_site = []
        for f in site_files:
            try:
                df = pd.read_csv(f, sep="\t")
                df = df[pd.to_numeric(df["pi"], errors="coerce").notnull()]
                if not df.empty:
                    df["pi"] = df["pi"].astype(float)
                    df_site.append(df[["site", "pi"]])
            except Exception:
                pass
        if df_site:
            merged = pd.concat(df_site)
            site_mean = merged.groupby("site")["pi"].mean().reset_index()
            site_mean.to_csv(os.path.join(out_dir, "merged_site_pi.csv"), index=False)
            print(f"    SNPGenie π: {len(df_site)} 样本 → {len(site_mean)} 位点")

    print(f"\n  完成 → {args.outdir}/")


if __name__ == "__main__":
    main()
