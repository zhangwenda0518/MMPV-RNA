#!/bin/bash
cd /home/zhangwenda/data-test
for d in */; do
  n=$(ls "$d" 2>/dev/null | wc -l)
  sample=$(ls "$d" 2>/dev/null | head -3 | tr '\n' '|')
  sz=$(du -sh --max-depth=0 "$d" 2>/dev/null | cut -f1)
  printf '%-30s %8s  %s  %s\n' "$d" "$sz" "$n entries" "$sample"
done
