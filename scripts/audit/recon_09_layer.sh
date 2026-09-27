#!/bin/bash
# 只读：09 分析层的表结构（第 3 条诉求的落点）
for p in ~/MMPV-paper/goji-virome/02_novel_virus/RNA-Lycium_barbarum_out \
         ~/MMPV-paper/onekp-virome/onekp-virus; do
  echo "############ $p"
  for sub in 09_Virome_Analysis/all_plant_analysis 09_Virome_Analysis/HQ_analysis 09_Virome_Analysis/rescue_detection 09_Virome_Analysis/viroid_analysis; do
    echo "=== $sub"
    ls -la "$p/$sub" 2>/dev/null | head -22
  done
  echo "=== 09 下所有 tsv（maxdepth 3）"
  find "$p/09_Virome_Analysis" -maxdepth 3 -type f -name "*.tsv" 2>/dev/null | head -30
done

echo
echo "############ 09 层现有表的表头（抽 all_plant_analysis 与 HQ_analysis 的 tsv）"
for f in $(find ~/MMPV-paper/goji-virome/02_novel_virus/RNA-Lycium_barbarum_out/09_Virome_Analysis -maxdepth 2 -type f -name "*.tsv" 2>/dev/null | head -8); do
  echo "--- $f"
  head -1 "$f" | tr '\t' '\n' | nl
done
