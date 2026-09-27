#!/bin/bash
# TIE_BREAK_ORDER 敏感性测试：同一引擎（v6.6 定稿行为）只换平票兜底顺序
#   S1 = ACVirus 挪到末位（把当前首位工具降权到底）
#   S2 = 完全倒序（最极端）
# 基线 = /tmp/rc66_count（默认顺序 ACVirus,CAT,VITAP,diamond_lca,genomad,metabuli,mmseqs）
set -u
SRC=/tmp/virus_classifier_analysis.R.cascade_v66
INT=/home/zhangwenda/MMPV-paper/goji-virome/02_novel_virus/RNA-Lycium_barbarum_out/05_Taxonomy
COMB=$INT/Votus.classed/Votus_combined_taxonomy.tsv
cd /tmp

echo "=== 源补丁 ==="
md5sum $SRC; grep -n '^TIE_BREAK_ORDER' $SRC

echo "=== 生成变体 ==="
python3 /tmp/make_tiebreak_variant.py $SRC /tmp/virus_classifier_analysis.R.sens_s1 \
  "CAT,VITAP,diamond_lca,genomad,metabuli,mmseqs,ACVirus"
python3 /tmp/make_tiebreak_variant.py $SRC /tmp/virus_classifier_analysis.R.sens_s2 \
  "mmseqs,metabuli,genomad,diamond_lca,VITAP,CAT,ACVirus"
for v in s1 s2; do
  grep -n '^TIE_BREAK_ORDER' /tmp/virus_classifier_analysis.R.sens_$v
  Rscript -e "invisible(parse('/tmp/virus_classifier_analysis.R.sens_$v')); cat('PARSE OK $v\n')" 2>&1 | tail -2
done

rm -rf /tmp/sens_s1 /tmp/sens_s2
mkdir -p /tmp/sens_s1 /tmp/sens_s2

for v in s1 s2; do
  echo "=== 跑 $v ==="
  start=$(date +%s)
  MMPV_CONSENSUS_MODE=cascade Rscript /tmp/virus_classifier_analysis.R.sens_$v \
    --combined $COMB --output /tmp/sens_$v > /tmp/sens_$v.log 2>&1
  echo "exit=$? elapsed=$(( $(date +%s) - start ))s"
  grep -aE "逐级淘汰|平票|科-属校准|自检" /tmp/sens_$v.log
done

echo
echo "=== 基线产物 ==="
ls -la /tmp/rc66_count/final_integrated_classification.tsv
md5sum /tmp/rc66_count/final_integrated_classification.tsv

for v in s1 s2; do
  echo
  echo "########## 敏感性 $v vs 基线（A=基线 B=$v） ##########"
  A=/tmp/rc66_count/final_integrated_classification.tsv \
  B=/tmp/sens_$v/final_integrated_classification.tsv \
  python3 /tmp/cmp_cells_lf.py > /tmp/sens_${v}_cmp.out 2>&1
  head -20 /tmp/sens_${v}_cmp.out
done

echo
echo "=== 变化 contig 数（只取两个 cmp 的合计行） ==="
grep -a "涉及 contig 合计" /tmp/sens_s1_cmp.out /tmp/sens_s2_cmp.out
