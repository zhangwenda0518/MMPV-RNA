#!/bin/bash
# 只读：全部物种项目的 05 分类表清单 + 逐项目科属不相容复算 + 09 目录布局
echo "=== [1] 全部 05 分类表（find） ==="
find ~/MMPV-paper -maxdepth 7 -type f -name "final_integrated_classification.tsv" 2>/dev/null \
  | grep -v "calibration" | sort > /tmp/all05.txt
while read f; do
  stat -c '%y  %10s  %n' "$f"
done < /tmp/all05.txt
echo "共 $(wc -l < /tmp/all05.txt) 张"

echo
echo "=== [2] 已有校准表 ==="
find ~/MMPV-paper -maxdepth 8 -type f -name "*calibrated*.tsv" 2>/dev/null | sort | while read f; do
  stat -c '%y  %10s  %n' "$f"
done

echo
echo "=== [3] 09 目录布局（抽 3 个项目） ==="
for d in ~/MMPV-paper/goji-virome/02_novel_virus/RNA-Lycium_barbarum_out \
         ~/MMPV-paper/goji-virome/02_novel_virus/RNA-Lycium_chinense_out \
         ~/MMPV-paper/onekp-virome/onekp-virus; do
  echo "--- $d/09_Virome_Analysis"
  ls -la "$d/09_Virome_Analysis" 2>/dev/null | head -30
  echo "--- $d/09b_Analysis_Verify"
  ls -la "$d/09b_Analysis_Verify" 2>/dev/null | head -10
done

echo
echo "=== [4] 逐项目科属不相容复算 ==="
cd /tmp
python3 recheck_famgenus_state.py $(cat /tmp/all05.txt) 2>&1 | grep -E "行数|参照表属数|不存在|跳过"
