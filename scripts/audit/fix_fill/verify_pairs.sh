#!/bin/bash
# 用法: verify_pairs.sh <label> <dataset_dir>
# 输出: 成品按键对比 / 相邻层矛盾(旧vs新) / 成品丢失格与新增格
L=$1
D=$2
I=$D/05_Taxonomy/Votus.integrated
C=$D/05_Taxonomy/Votus.classed
PO=$I.bak_fillfix_20260916/final_integrated_classification.tsv
PN=$I/final_integrated_classification.tsv
CO=$C/Votus_combined_taxonomy.tsv.bak_fillfix_20260916
CN=$C/Votus_combined_taxonomy.tsv

echo "########## $L"
echo "--- 成品按键对比"
python3 /tmp/diff_prod_by_key.py "$PO" "$PN" --top 6 --quiet
echo "--- 成品相邻层矛盾"
echo "OLD:"; python3 /tmp/chimera_rows.py "$PO" --max 0 | tail -2
echo "NEW:"; python3 /tmp/chimera_rows.py "$PN" --max 0 | tail -2
echo "--- 成品单元格丢失/新增"
python3 /tmp/probe_lost_cells.py --prod-old "$PO" --prod-new "$PN" --comb-old "$CO" --comb-new "$CN" --max 0
