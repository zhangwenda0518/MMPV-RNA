#!/bin/bash
# 备份黑名单剔除涉及的产物, 建立回滚点
set -e
OUT=/home/zhangwenda/data-test/out
BAK=/home/zhangwenda/data-test/blacklist_bak_20260902
mkdir -p "$BAK"

echo "=== 备份开始 ==="
FILES=(
  "06_HostPrediction/host_classified_fasta/Plant.classified.fasta"
  "07_Checkv/Plant.fasta"
  "08_Rescue/Plant/input_centroids.fasta"
  "08_Rescue/HQ_plant_viruses.fasta"
  "09_Virome_Analysis/HQ_plant_viruses.fasta"
  "09b_Analysis_Verify/class_KEEP.fasta"
  "09b_Analysis_Verify/acvirus_classify/contigs.fna"
  "09b_Analysis_Verify/virus_validation/keep_candidates.fasta"
  "09b_Analysis_Verify/virus_validation/review_candidates.fasta"
)
for rel in "${FILES[@]}"; do
  src="$OUT/$rel"
  if [ -f "$src" ]; then
    dst="$BAK/$(echo "$rel" | tr '/' '_')"
    cp -p "$src" "$dst"
    echo "  ✓ $rel ($(grep -c '^>' "$src" 2>/dev/null || echo '?') 条)"
  else
    echo "  ✗ 缺失: $rel"
  fi
done

# 受影响科的树目录
echo ""
echo "=== 备份 9 个受影响科的树 ==="
TREES=/home/zhangwenda/data-test/out/09b_Analysis_Verify/acvirus_trees
for fam in Tombusviridae Partitiviridae Aspiviridae Caulimoviridae Metaviridae Geminiviridae Potyviridae Secoviridae Pithoviridae; do
  if [ -d "$TREES/$fam" ]; then
    mkdir -p "$BAK/trees/$fam"
    cp -rp "$TREES/$fam/." "$BAK/trees/$fam/" 2>/dev/null || true
    echo "  ✓ $fam"
  else
    echo "  - $fam (无此目录)"
  fi
done

# SDT 矩阵 (若存在)
echo ""
echo "=== 备份 SDT 矩阵 ==="
SDT=/home/zhangwenda/data-test/out/09b_Analysis_Verify/SDT_matrix
if [ -d "$SDT" ]; then
  du -sh "$SDT" 2>/dev/null || true
  ls "$SDT" | head -20
fi

echo ""
echo "=== 备份完成: $BAK ==="
du -sh "$BAK"
