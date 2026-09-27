#!/bin/bash
# 1) 把 24 个 09 层目标文件从 .bak_dnarna_20260830 复原（清掉 v1 的错误追加，含 2 个 CRLF 文件）
# 2) 用 v2（字节级、兼容 CRLF）重新落地
set -u
cd /home/zhangwenda/MMPV-paper

mapfile -t FILES < <(find . \( -name 'All_plant.viruses_info.tsv' -o -name 'HQ_plant_viruses_info.tsv' -o -name 'all_plant_viruses_genus_summary.tsv' \) -path '*09_Virome_Analysis/*' | sort)

echo "候选文件数: ${#FILES[@]}"
n=0
for f in "${FILES[@]}"; do
  if [ -f "$f.bak_dnarna_20260830" ]; then
    cp -p "$f.bak_dnarna_20260830" "$f" && n=$((n+1))
  fi
done
echo "已复原: $n"

echo "--- 复原后列数抽样 ---"
for f in "${FILES[@]}"; do
  c=$(head -n 1 "$f" | tr '\t' '\n' | wc -l)
  printf "  %3s  %s\n" "$c" "${f##*MMPV-paper/}"
done

echo
echo "=== 运行 v2 --apply ==="
python3 /tmp/dnarna_apply2.py --apply > /tmp/dnarna_apply2_run.log 2>&1
echo "EXIT=$?"
echo "PASS: $(grep -c 'PASS' /tmp/dnarna_apply2_run.log)  FAIL: $(grep -c 'FAIL' /tmp/dnarna_apply2_run.log)"
tail -n 3 /tmp/dnarna_apply2_run.log
