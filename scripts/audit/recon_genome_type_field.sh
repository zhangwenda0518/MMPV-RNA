#!/bin/bash
# 只读：integrated_summary.py 的 genome_type 字段来源与产物落地位置（限制深度，避免全树 grep）
P=~/MMPV-RNA/virome_discovery_pipeline
G=~/MMPV-paper/goji-virome/02_novel_virus/RNA-Lycium_barbarum_out
O=~/MMPV-paper/onekp-virome/onekp-virus

echo "=== [1] integrated_summary.py 关键行 70-110 ==="
awk 'NR>=70 && NR<=110 {printf "%4d %s\n", NR, $0}' $P/integrated_summary.py

echo
echo "=== [2] integrated_summary.py 关键行 190-265 ==="
awk 'NR>=190 && NR<=265 {printf "%4d %s\n", NR, $0}' $P/integrated_summary.py

echo
echo "=== [3] 枸杞树 maxdepth 3 的 tsv 里含 genome_type 的 ==="
find $G -maxdepth 3 -type f -name "*.tsv" 2>/dev/null | while read f; do
  if head -1 "$f" 2>/dev/null | grep -q "genome_type"; then echo "HIT $f"; fi
done
echo "(以上为空 = 该深度内没有 genome_type 列)"

echo
echo "=== [4] 枸杞树 maxdepth 2 目录清单 ==="
find $G -maxdepth 2 -type d 2>/dev/null | sort
