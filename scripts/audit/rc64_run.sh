#!/bin/bash
# v6.4 双模式重跑：legacy（对照 v6.3 legacy 看修复 1/4/5 效果） + cascade
set -u
R=/tmp/virus_classifier_analysis.R.cascade_v64
echo "=== parse check ==="
Rscript -e "invisible(parse('$R')); cat('PARSE OK\n')" 2>&1 | tail -5
echo "=== md5 ==="
md5sum $R
INT=/home/zhangwenda/MMPV-paper/goji-virome/02_novel_virus/RNA-Lycium_barbarum_out/05_Taxonomy
COMB=$INT/Votus.classed/Votus_combined_taxonomy.tsv
echo "=== input ==="
ls -la $COMB
rm -rf /tmp/rc64_legacy /tmp/rc64_cascade
mkdir -p /tmp/rc64_legacy /tmp/rc64_cascade
cd /tmp
echo "=== run legacy ==="
start=$(date +%s)
MMPV_CONSENSUS_MODE=legacy Rscript $R --combined $COMB --output /tmp/rc64_legacy > /tmp/rc64_legacy.log 2>&1
echo "legacy exit=$? elapsed=$(( $(date +%s) - start ))s"
echo "=== run cascade ==="
start=$(date +%s)
MMPV_CONSENSUS_MODE=cascade Rscript $R --combined $COMB --output /tmp/rc64_cascade > /tmp/rc64_cascade.log 2>&1
echo "cascade exit=$? elapsed=$(( $(date +%s) - start ))s"
echo "=== outputs ==="
wc -l /tmp/rc64_legacy/final_integrated_classification.tsv /tmp/rc64_cascade/final_integrated_classification.tsv
ls -la /tmp/rc64_legacy/ /tmp/rc64_cascade/
