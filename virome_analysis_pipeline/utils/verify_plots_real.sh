#!/usr/bin/env bash
# 服务器真实数据验证: QST 全图 (PSTVd af_matrix) + vcf_merge LD菱形
set -u
ROOT=~/virus/data-2026/known_virus_all_v2
PIPE=~/MMPV-RNA/virome_analysis_pipeline
cd "$PIPE" || exit 1

echo "== 1. QST 真实数据 (PSTVd af_matrix.tsv, Host 分组)"
rm -rf "$ROOT/test_af_merge/qst_real"
python utils/virus_haplotype_qst.py \
  -m "$ROOT/04_post_analysis/Potato_spindle_tuber_viroid_NC_002030.1/vcf_merge/matrices/af_matrix.tsv" \
  -M "$ROOT/../sra_rna.data4/global_metadata/Global_Unified_Metadata_Core13.tsv" \
  -o "$ROOT/test_af_merge/qst_real" --label PSTVd --hamming-threshold 0.15 \
  > "$ROOT/test_af_merge/qst_real.log" 2>&1
echo "  QST EXIT=$?"
ls "$ROOT/test_af_merge/qst_real/" 2>/dev/null | head -15
tail -6 "$ROOT/test_af_merge/qst_real.log"

echo ""
echo "== 2. vcf_merge 重跑 (PSTVd, 验证 LD 菱形图)"
rm -rf "$ROOT/test_af_merge/PSTVd_plots"
python virus_vcf_pipeline.py \
  -d "$ROOT/03_variants/virus-variants/Potato_spindle_tuber_viroid_NC_002030.1" \
  -o "$ROOT/test_af_merge/PSTVd_plots" --prefix NC_002030.1 \
  --visualize --qc --snp-matrix --tree --dist-metrics both \
  --pca-method genotype --ld --ivar \
  > "$ROOT/test_af_merge/PSTVd_plots.log" 2>&1
echo "  MERGE EXIT=$?"
ls "$ROOT/test_af_merge/PSTVd_plots/figs/" 2>/dev/null | grep -E 'ld|diamond|gene' | head
grep -E 'LD菱形|流程完成|Error|Traceback' "$ROOT/test_af_merge/PSTVd_plots.log" | head -5
echo DONE
