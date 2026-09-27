#!/usr/bin/env python3
"""复刻 report_pipeline.generate_sankey 的 plant_final_taxonomy.tsv 生成逻辑。

原逻辑：读 06_HostPrediction/ensemble_host_summary.tsv，取 Final_Host == "Plant" 的 contig_id 集合，
再从 final_integrated_classification.tsv 逐行按首列（去引号）过滤，写 report_dir/plant_final_taxonomy.tsv。
本脚本不依赖 polars，用 csv 模块等价实现。
"""
import csv
import sys

host_summary, final_tax, out_path = sys.argv[1], sys.argv[2], sys.argv[3]

plant_ids = set()
with open(host_summary, newline="") as fh:
    rd = csv.DictReader(fh, delimiter="\t")
    if "Final_Host" not in (rd.fieldnames or []):
        print("ERROR: Final_Host column not found; header =", rd.fieldnames)
        sys.exit(2)
    for row in rd:
        if (row.get("Final_Host") or "").strip() == "Plant":
            cid = (row.get("contig_id") or "").strip()
            if cid:
                plant_ids.add(cid)
print("plant_ids:", len(plant_ids))

kept = 0
with open(final_tax, newline="") as tf, open(out_path, "w", newline="") as pf:
    header = tf.readline()
    pf.write(header)
    for line in tf:
        cid = line.split("\t")[0].strip().strip('"')
        if cid in plant_ids:
            pf.write(line)
            kept += 1
print("kept_lines:", kept)
