#!/bin/bash
# 只读：谁定义了 repair_family_from_genus（R 层现状 vs 备份 vs 本地脚本）
echo "=== [1] 管线目录内所有 .R 里的 repair_family_from_genus ==="
grep -rn "repair_family_from_genus" ~/MMPV-RNA/virome_discovery_pipeline --include="*.R" 2>/dev/null | head -20
echo "(空 = 管线内任何 R 文件都没有这个函数)"
echo
echo "=== [2] 管线目录内 R 文件清单 + md5 ==="
find ~/MMPV-RNA/virome_discovery_pipeline -maxdepth 1 -name "*.R" -o -maxdepth 1 -name "*.R.bak*" 2>/dev/null | sort | while read f; do
  md5sum "$f"; stat -c '        %y  %s bytes' "$f"
done
echo
echo "=== [3] /tmp 里所有提到 repair_family_from_genus 的文件 ==="
grep -rln "repair_family_from_genus" /tmp --include="*.py" --include="*.R" --include="*.md" 2>/dev/null | head -20
echo
echo "=== [4] enforce_rank_containment 函数体（:444 起 40 行） ==="
awk 'NR>=444 && NR<=484 {printf "%4d %s\n", NR, $0}' ~/MMPV-RNA/virome_discovery_pipeline/virus_classifier_analysis.R
