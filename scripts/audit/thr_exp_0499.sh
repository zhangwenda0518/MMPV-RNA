#!/bin/bash
# 门槛实验：把 CASCADE_MIN_SHARE 从 0.5 降到 0.499，让「半数票」(1/2, 2/4, 3/6) 也触发淘汰。
# 只写 /tmp，不碰线上产物。
set -u
SCRIPT=/home/zhangwenda/MMPV-RNA/virome_discovery_pipeline/virus_classifier_analysis.R
REF=/home/zhangwenda/database/taxonomy/genus_family_ref.tsv
D=/home/zhangwenda/MMPV-paper/goji-virome/02_novel_virus/RNA-Lycium_barbarum_out
C="$D/05_Taxonomy/Votus.classed/Votus_combined_taxonomy.tsv"
O=/tmp/thr0499_barbarum
mkdir -p "$O"
echo "=== $(date '+%F %T') START thr=0.499 ==="
( cd "$O" && MMPV_GENUS_FAMILY_REF="$REF" MMPV_CASCADE_MIN_SHARE=0.499 Rscript "$SCRIPT" --combined "$C" --output "$O" ) > "$O/run.log" 2>&1
echo "rc=$? lines=$(wc -l < "$O/final_integrated_classification.tsv" 2>/dev/null) md5=$(md5sum "$O/final_integrated_classification.tsv" 2>/dev/null | cut -c1-32)"
grep -e 逐级淘汰 -e 平票 "$O/run.log" | tail -5
echo "=== $(date '+%F %T') DONE ==="
