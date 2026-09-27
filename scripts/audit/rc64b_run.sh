#!/bin/bash
# v6.4b 双模式重跑：与 /tmp/rc64_* 做产物字节比对，证明「注释与日志串改动不触逻辑」
set -u
R=/tmp/virus_classifier_analysis.R.cascade_v64
echo "=== parse check ==="
Rscript -e "invisible(parse('$R')); cat('PARSE OK\n')" 2>&1 | tail -5
echo "=== md5 ==="
md5sum $R
INT=/home/zhangwenda/MMPV-paper/goji-virome/02_novel_virus/RNA-Lycium_barbarum_out/05_Taxonomy
COMB=$INT/Votus.classed/Votus_combined_taxonomy.tsv
ls -la $COMB
rm -rf /tmp/rc64b_legacy /tmp/rc64b_cascade
mkdir -p /tmp/rc64b_legacy /tmp/rc64b_cascade
cd /tmp
echo "=== run legacy ==="
start=$(date +%s)
MMPV_CONSENSUS_MODE=legacy Rscript $R --combined $COMB --output /tmp/rc64b_legacy > /tmp/rc64b_legacy.log 2>&1
echo "legacy exit=$? elapsed=$(( $(date +%s) - start ))s"
echo "=== run cascade ==="
start=$(date +%s)
MMPV_CONSENSUS_MODE=cascade Rscript $R --combined $COMB --output /tmp/rc64b_cascade > /tmp/rc64b_cascade.log 2>&1
echo "cascade exit=$? elapsed=$(( $(date +%s) - start ))s"
wc -l /tmp/rc64b_legacy/final_integrated_classification.tsv /tmp/rc64b_cascade/final_integrated_classification.tsv
