#!/bin/bash
# barbarum 单项目 A/B 对账：同一 combined 输入分别跑 legacy / cascade，产物写 /tmp，不动盘上正式产物
set -u
SD=/home/zhangwenda/MMPV-RNA/virome_discovery_pipeline
RSCRIPT=/tmp/virus_classifier_analysis.R.cascade_v63
COMB=/home/zhangwenda/MMPV-paper/goji-virome/02_novel_virus/RNA-Lycium_barbarum_out/05_Taxonomy/Votus.classed/Votus_combined_taxonomy.tsv
REF=/home/zhangwenda/database/taxonomy/genus_family_ref.tsv
export MMPV_GENUS_FAMILY_REF=$REF

echo "=== start $(date) ==="
rm -rf /tmp/rc_legacy /tmp/rc_cascade
mkdir -p /tmp/rc_legacy /tmp/rc_cascade

cd /tmp/rc_legacy
echo "--- legacy ---"
MMPV_CONSENSUS_MODE=legacy Rscript $RSCRIPT --combined $COMB --output /tmp/rc_legacy > /tmp/rc_legacy.log 2>&1
echo "legacy exit=$? $(date)"

cd /tmp/rc_cascade
echo "--- cascade ---"
MMPV_CONSENSUS_MODE=cascade Rscript $RSCRIPT --combined $COMB --output /tmp/rc_cascade > /tmp/rc_cascade.log 2>&1
echo "cascade exit=$? $(date)"

echo "=== done $(date) ==="
ls -la /tmp/rc_legacy/final_integrated_classification.tsv /tmp/rc_cascade/final_integrated_classification.tsv
echo DONE > /tmp/rc_done.flag
