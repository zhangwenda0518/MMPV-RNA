#!/bin/bash
# 等价性验证：新版默认（CASCADE_MIN_SHARE=0.5 + ">= 门槛" 语义）应与 0.499 实验产物逐字节相同。
# 期望 md5 = a3ead5b5bc4fe290a0e9922497d2fc35（/tmp/thr0499_barbarum 用 MMPV_CASCADE_MIN_SHARE=0.499 + ">" 得到）。
set -u
SCRIPT=/home/zhangwenda/MMPV-RNA/virome_discovery_pipeline/virus_classifier_analysis.R
REF=/home/zhangwenda/database/taxonomy/genus_family_ref.tsv
D=/home/zhangwenda/MMPV-paper/goji-virome/02_novel_virus/RNA-Lycium_barbarum_out
C="$D/05_Taxonomy/Votus.classed/Votus_combined_taxonomy.tsv"
O=/tmp/thr_eq_barbarum
rm -rf "$O"; mkdir -p "$O"
echo "=== $(date '+%F %T') START default gate=6.7 >=0.5 ==="
( cd "$O" && MMPV_GENUS_FAMILY_REF="$REF" Rscript "$SCRIPT" --combined "$C" --output "$O" ) > "$O/run.log" 2>&1
echo "rc=$? lines=$(wc -l < "$O/final_integrated_classification.tsv" 2>/dev/null) md5=$(md5sum "$O/final_integrated_classification.tsv" 2>/dev/null | cut -c1-32)"
grep -e 逐级淘汰 -e 平票 "$O/run.log" | tail -5
echo "expected md5=a3ead5b5bc4fe290a0e9922497d2fc35"
echo "=== $(date '+%F %T') DONE ==="
