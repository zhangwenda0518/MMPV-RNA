#!/bin/bash
# rc64_* 与 rc64b_* 产物字节比对（数据产物必须全 SAME；差异只允许出现在版本戳与日志）
cd /tmp
FILES="final_integrated_classification.tsv fill_stats.tsv tool_weights.tsv species_containment_blanked.tsv taxonomy_gate_check.tsv analysis_summary.txt standardized_ACVirus.tsv standardized_CAT.tsv standardized_VITAP.tsv standardized_diamond_lca.tsv standardized_genomad.tsv standardized_metabuli.tsv standardized_mmseqs.tsv agreement_stats.tsv consistency_summary.tsv consistent_summary.tsv intersection_matrix.tsv common_ids.txt upset_data.json comparison_family.tsv comparison_genus.tsv comparison_species.tsv comparison_realm.tsv"
for m in legacy cascade; do
  echo "=== mode=$m ==="
  for f in $FILES; do
    a=$(md5sum "rc64_$m/$f" 2>/dev/null | cut -d' ' -f1)
    b=$(md5sum "rc64b_$m/$f" 2>/dev/null | cut -d' ' -f1)
    if [ -z "$a" ] && [ -z "$b" ]; then echo "MISSING $f"; continue; fi
    if [ "$a" = "$b" ]; then echo "SAME  $f"; else echo "DIFF  $f  $a -> $b"; fi
  done
done
echo "=== 版本戳 diff（预期仅 script_md5 变） ==="
diff rc64_legacy/taxonomy_gate_stamp.tsv rc64b_legacy/taxonomy_gate_stamp.tsv
echo "--- cascade ---"
diff rc64_cascade/taxonomy_gate_stamp.tsv rc64b_cascade/taxonomy_gate_stamp.tsv
echo "=== 日志 diff（预期仅计时标签行） ==="
diff rc64_legacy.log rc64b_legacy.log | head -20
echo "--- cascade ---"
diff rc64_cascade.log rc64b_cascade.log | head -20
