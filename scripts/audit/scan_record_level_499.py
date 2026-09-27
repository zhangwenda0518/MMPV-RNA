"""记录级残留扫描: 只看 fasta 头部 / TSV 首列, 避免 hit 列撞名。
只统计 contig 专属 id (含 _clean_NODE_) 才是真残留。
"""
import glob
import os

BASE = "/home/zhangwenda/data-test/out"
all_ids = set(l.strip() for l in open("/tmp/stale_499.txt") if l.strip())
ids = {i for i in all_ids if "_clean_NODE_" in i}  # contig 专属
print(f"499 清单: {len(all_ids)} 条, 其中 contig 专属 {len(ids)} 条\n")


def scan_fasta_headers(pattern, tag):
    files = glob.glob(pattern, recursive=True)
    n_file, n_hit = 0, 0
    for fp in files:
        try:
            with open(fp, errors="replace") as f:
                hs = {l[1:].split()[0] for l in f if l.startswith(">")}
        except Exception:
            continue
        ov = hs & ids
        if ov:
            n_file += 1
            n_hit += len(ov)
    print(f"{tag}: {len(files)} 文件, 含残留 {n_file} 文件 / {n_hit} 条头部")
    return n_hit


def scan_tsv_col1(path, tag):
    if not os.path.exists(path):
        print(f"{tag}: 缺")
        return
    c = 0
    with open(path, errors="replace") as f:
        for i, l in enumerate(f):
            if i == 0:
                continue
            if l.split("\t")[0].strip().strip('"') in ids:
                c += 1
    print(f"{tag}: 首列残留 {c} 条")


print("=== 06 HostPrediction ===")
scan_tsv_col1(f"{BASE}/06_HostPrediction/ensemble_host_summary.tsv", "  ensemble_host_summary.tsv")
scan_fasta_headers(f"{BASE}/06_HostPrediction/**/*.fa", "  全部 .fa")
scan_fasta_headers(f"{BASE}/06_HostPrediction/**/*.fasta", "  全部 .fasta")

print("\n=== 07 Checkv ===")
scan_tsv_col1(f"{BASE}/07_Checkv/Plant/completeness.tsv", "  Plant/completeness.tsv")
scan_fasta_headers(f"{BASE}/07_Checkv/**/*.fasta", "  全部 .fasta")

print("\n=== 08 Rescue ===")
scan_fasta_headers(f"{BASE}/08_Rescue/Plant/centroids/*.fasta", "  Plant/centroids/*.fasta")
scan_fasta_headers(f"{BASE}/08_Rescue/**/branch_*/chunk_*/chunk.fasta", "  Plant chunks")
scan_fasta_headers(f"{BASE}/08_Rescue/known/centroids/*.fasta", "  known/centroids/*.fasta")

print("\n=== 09 Virome ===")
scan_tsv_col1(f"{BASE}/09_Virome_Analysis/all_plant_analysis/All_plant.viruses_info.tsv",
              "  All_plant.viruses_info.tsv")
scan_fasta_headers(f"{BASE}/09_Virome_Analysis/**/*.fasta", "  全部 .fasta")

print("\n=== 09b Verify ===")
scan_tsv_col1(f"{BASE}/09b_Analysis_Verify/virus_validation/rescue_evidence_scored.tsv", "  scored.tsv")
scan_tsv_col1(f"{BASE}/09b_Analysis_Verify/virus_validation/cdd_evidence_report.tsv", "  cdd_evidence_report.tsv")
scan_fasta_headers(f"{BASE}/09b_Analysis_Verify/**/*.fasta", "  全部 .fasta")

print("\n=== 10 Reports ===")
for fn in ["final_integrated_classification.tsv", "rescue_report.tsv", "rescue_detection_summary.tsv"]:
    scan_tsv_col1(f"{BASE}/10_Reports/{fn}", f"  {fn}")
