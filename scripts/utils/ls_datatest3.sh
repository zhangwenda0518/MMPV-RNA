#!/bin/bash
cd /home/zhangwenda/data-test
for d in out out2 out3 out4 out6 out7 out8 out9 out11 test test2 tmp; do
  if [ -d "$d" ]; then
    sz=$(du -sh --max-depth=0 "$d" 2>/dev/null | cut -f1)
    n=$(ls "$d" 2>/dev/null | wc -l)
    sm=$(ls "$d" 2>/dev/null | head -3 | tr '\n' '|')
    printf '%-8s %8s  %4s entries  %s\n' "$d" "$sz" "$n" "$sm"
  fi
done
