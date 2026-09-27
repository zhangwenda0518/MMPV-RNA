#!/bin/bash
# onekp 报告层补跑：两张 sankey + plant_final_taxonomy.tsv
# 先用 report_pipeline.py 的完全相同参数，输出到 /tmp/sankey_fix/out/，确认无误后再回填 10_Reports。
# 不改管线脚本；只读数据集产物。
set -u

D=/home/zhangwenda/MMPV-paper/onekp-virome/onekp-virus
REP=$D/10_Reports
FINAL=$D/05_Taxonomy/Votus.integrated/final_integrated_classification.tsv
HOST=$D/06_HostPrediction/ensemble_host_summary.tsv
SK=/home/zhangwenda/MMPV-RNA/stats/taxonomic_sankey.py
BUILD=/tmp/sankey_fix/build_plant_tax.py
OUT=/tmp/sankey_fix/out
mkdir -p $OUT

echo "=== INPUT"
ls -l --time-style=long-iso $FINAL $SK $HOST
md5sum $FINAL

echo "=== EXISTING PRODUCTS (v6.6 前版本)"
ls -l --time-style=long-iso $REP/classification_sankey.html $REP/classification_sankey_plant.html $REP/plant_final_taxonomy.tsv
md5sum $REP/classification_sankey.html $REP/classification_sankey_plant.html $REP/plant_final_taxonomy.tsv

echo "=== STEP1 PLANT FILTER"
S=$(date +%s)
python3 $BUILD $HOST $FINAL $OUT/plant_final_taxonomy.tsv
echo EXIT_BUILD=$? ELAPSED_BUILD=$(( $(date +%s) - S ))
ls -l --time-style=long-iso $OUT/plant_final_taxonomy.tsv

echo "=== STEP2 ALL SANKEY"
S=$(date +%s)
python3 $SK -i $FINAL -o $OUT/classification_sankey.html --format html --min-flow 1 --min-genus-flow 10 --palette set3 --height 600 --width 800 --node-pad 30 --label-truncate 25 --font-size 9 --title-font-size 14 --title ""
echo EXIT_ALL=$? ELAPSED_ALL=$(( $(date +%s) - S ))
ls -l --time-style=long-iso $OUT/classification_sankey.html

echo "=== STEP3 PLANT SANKEY"
S=$(date +%s)
python3 $SK -i $OUT/plant_final_taxonomy.tsv -o $OUT/classification_sankey_plant.html --format html --min-flow 1 --min-genus-flow 5 --palette set3 --height 600 --width 800 --node-pad 30 --label-truncate 25 --font-size 9 --title-font-size 14 --title ""
echo EXIT_PLANT=$? ELAPSED_PLANT=$(( $(date +%s) - S ))
ls -l --time-style=long-iso $OUT/classification_sankey_plant.html

echo "=== OUTPUT md5"
md5sum $OUT/classification_sankey.html $OUT/classification_sankey_plant.html $OUT/plant_final_taxonomy.tsv
echo "=== DONE"
