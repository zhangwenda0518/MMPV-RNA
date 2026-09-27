#!/bin/bash
# rc64_share_run.sh —— CASCADE_MIN_SHARE 敏感性扫描：同一输入只改票首门槛，看淘汰格次与产物差异
set -u
R=/tmp/virus_classifier_analysis.R.cascade_v64
INT=/home/zhangwenda/MMPV-paper/goji-virome/02_novel_virus/RNA-Lycium_barbarum_out/05_Taxonomy
COMB=$INT/Votus.classed/Votus_combined_taxonomy.tsv
cd /tmp
for s in 0.3 0.4 0.6 0.7; do
  d=/tmp/share_$s
  rm -rf $d; mkdir -p $d
  start=$(date +%s)
  MMPV_CONSENSUS_MODE=cascade MMPV_CASCADE_MIN_SHARE=$s Rscript $R --combined $COMB --output $d > /tmp/share_$s.log 2>&1
  echo "share=$s exit=$? elapsed=$(( $(date +%s) - start ))s  $(grep -h '逐级淘汰' /tmp/share_$s.log | tail -1)"
done
echo "=== 基准 0.5 ==="
grep -h '逐级淘汰' /tmp/rc64b_cascade.log | tail -1
echo "=== 产物比对 ==="
python3 /tmp/cmp_share.py
