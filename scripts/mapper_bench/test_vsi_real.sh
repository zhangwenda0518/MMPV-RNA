#!/bin/bash
# 真实调用 VSI, 验证 minibwa 路径在完整流程中正常工作
# 用 rescue 阶段的原始调用形态
set -e

cd /home/zhangwenda/MMPV-RNA/virome_discovery_pipeline

SCAF=/home/zhangwenda/data-test/out/08_Rescue/Plant/branch_b/scaffolds/ERR2040118_clean_NODE_150_length_3726_cov_29.718708.fasta
R1=/home/zhangwenda/data-test/out/08_Rescue/Plant/branch_b/merged_reads/ERR2040262_clean_NODE_403_length_3351_cov_10.602719_merged_R1.fa.gz
R2=/home/zhangwenda/data-test/out/08_Rescue/Plant/branch_b/merged_reads/ERR2040262_clean_NODE_403_length_3351_cov_10.602719_merged_R2.fa.gz

echo "=== 检查测试输入 ==="
ls -la $SCAF
ls -lh $R1 $R2
echo ""

OUT=/tmp/vsi_real_test
rm -rf $OUT && mkdir -p $OUT

echo "=== 真实调用 VSI (走 minibwa 默认路径) ==="
echo "命令: python3 Virseqimprover.py -1 <R1> -2 <R2> -scaffold <scaf> -o $OUT -t 10"
echo "--- 输出 ---"
timeout 900 python3 Virseqimprover.py \
  -1 $R1 -2 $R2 \
  -scaffold $SCAF \
  -o $OUT \
  -t 10 2>&1 | tail -40

echo ""
echo "=== mapper 决策日志 ==="
grep -E "\[mapper\]" $OUT/*.log 2>/dev/null || echo "(检查 run.sh 与产物)"
echo ""
echo "=== 产物 ==="
ls -la $OUT/ 2>/dev/null | head -20
