#!/bin/bash
# 09 层文件清单（只读）
set -u
OUT=/tmp/d09_files.txt
: > "$OUT"
grep -v '10_Reports' /tmp/d09_dirs.txt | while read -r d; do
  [ -d "$d" ] || continue
  find "$d" -maxdepth 3 -type f \( -name '*.tsv' -o -name '*.csv' \) -printf '%s\t%p\n' 2>/dev/null
done | sort -k2 >> "$OUT"
echo "文件数: $(wc -l < "$OUT")"
awk -F'\t' '{printf "%10.1fKB  %s\n", $1/1024, $2}' "$OUT"
