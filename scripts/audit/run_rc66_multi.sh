#!/bin/bash
# 对全部数据集重跑 v6.6 共识（只写 /tmp，不碰线上产物）
set -u
SCRIPT=/home/zhangwenda/MMPV-RNA/virome_discovery_pipeline/virus_classifier_analysis.R
REF=/home/zhangwenda/database/taxonomy/genus_family_ref.tsv
ROOT=/tmp/rc66_multi
mkdir -p "$ROOT"
for D in /home/zhangwenda/MMPV-paper/goji-virome/02_novel_virus/RNA-*_out /home/zhangwenda/MMPV-paper/onekp-virome/onekp-virus; do
  [ -d "$D/05_Taxonomy" ] || continue
  B=$(basename "$D")
  C="$D/05_Taxonomy/Votus.classed/Votus_combined_taxonomy.tsv"
  O="$ROOT/$B"
  mkdir -p "$O"
  echo "=== $(date '+%F %T') START $B combined=$(stat -c%s "$C" 2>/dev/null) ==="
  ( cd "$O" && MMPV_GENUS_FAMILY_REF="$REF" Rscript "$SCRIPT" --combined "$C" --output "$O" ) > "$O/run.log" 2>&1
  RC=$?
  N=$(wc -l < "$O/final_integrated_classification.tsv" 2>/dev/null || echo 0)
  echo "=== $(date '+%F %T') DONE $B rc=$RC prod_lines=$N md5=$(md5sum "$O/final_integrated_classification.tsv" 2>/dev/null | cut -c1-32) ==="
done
echo "ALL DONE $(date '+%F %T')"
