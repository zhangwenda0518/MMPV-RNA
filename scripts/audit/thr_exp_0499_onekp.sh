#!/bin/bash
# 门槛实验（onekp 规模验证）：MMPV_CASCADE_MIN_SHARE=0.499，只写 /tmp。
set -u
SCRIPT=/home/zhangwenda/MMPV-RNA/virome_discovery_pipeline/virus_classifier_analysis.R
REF=/home/zhangwenda/database/taxonomy/genus_family_ref.tsv
D=/home/zhangwenda/MMPV-paper/onekp-virome/onekp-virus
C="$D/05_Taxonomy/Votus.classed/Votus_combined_taxonomy.tsv"
O=/tmp/thr0499_onekp
mkdir -p "$O"
echo "=== $(date '+%F %T') START onekp thr=0.499 ==="
( cd "$O" && MMPV_GENUS_FAMILY_REF="$REF" MMPV_CASCADE_MIN_SHARE=0.499 Rscript "$SCRIPT" --combined "$C" --output "$O" ) > "$O/run.log" 2>&1
echo "rc=$? lines=$(wc -l < "$O/final_integrated_classification.tsv" 2>/dev/null) md5=$(md5sum "$O/final_integrated_classification.tsv" 2>/dev/null | cut -c1-32)"
grep -e 逐级淘汰 "$O/run.log" | tail -3
echo "=== $(date '+%F %T') DONE ==="
