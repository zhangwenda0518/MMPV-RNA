#!/bin/bash
# 收尾核对：8 数据集 gate 戳 / 成品行数 / 报告层 rescue 同步
G=/home/zhangwenda/MMPV-paper/goji-virome/02_novel_virus
N=/home/zhangwenda/MMPV-paper/onekp-virome/onekp-virus
echo "===== DATE"; date
echo "===== 8 数据集"
printf '%-12s %-5s %-11s %-11s %-10s %-10s %-17s %s\n' DS GATE STAMP_ROWS PROD_ROWS R10_MD5 R08_MD5 R10_MTIME FLAG
for D in $G/RNA-Alternaria_alternata_out $G/RNA-Aphis_gossypii_out $G/RNA-Fusarium_nematophilum_out $G/RNA-Lycium_amarum_out $G/RNA-Lycium_barbarum_out $G/RNA-Lycium_chinense_out $G/RNA-Lycium_ruthenicum_out $N; do
  L=$(basename "$D")
  P="$D/05_Taxonomy/Votus.integrated/final_integrated_classification.tsv"
  GATE="$D/05_Taxonomy/Votus.integrated/taxonomy_gate_stamp.tsv"
  R10="$D/10_Reports/rescue_report.tsv"
  R08="$D/08_Rescue/Plant/rescue_report.tsv"
  gv=$(awk -F'\t' 'NR==2{print $1}' "$GATE" 2>/dev/null)
  rows=$(awk -F'\t' 'NR==2{print $6}' "$GATE" 2>/dev/null)
  nprod=$(( $(wc -l < "$P") - 1 ))
  m10=$(md5sum "$R10" 2>/dev/null | cut -c1-8)
  m08=$(md5sum "$R08" 2>/dev/null | cut -c1-8)
  t10=$(stat -c %y "$R10" 2>/dev/null | cut -c1-16)
  flag=OK
  [ -n "$m10" ] && [ "$m10" != "$m08" ] && flag=R10_08_DIFF
  [ "$rows" != "$nprod" ] && flag="$flag ROWS_MISMATCH"
  printf '%-12s %-5s %-11s %-11s %-10s %-10s %-17s %s\n' "$L" "$gv" "$rows" "$nprod" "$m10" "$m08" "$t10" "$flag"
done
echo "===== onekp 10_Reports 最新 8 个文件"
ls -lt --time-style=+%m-%d_%H:%M "$N/10_Reports" | head -9
echo "===== onekp 09b（应仍为上一轮，未随本轮刷新）"
ls -lt --time-style=+%m-%d_%H:%M "$N/09b_Analysis_Verify" 2>/dev/null | head -5
ls -lt --time-style=+%m-%d_%H:%M "$N/09b_Analysis_Verify/virus_validation" 2>/dev/null | head -5
echo "===== 备份目录"
ls -d $G/RNA-Lycium_barbarum_out/05_Taxonomy/Votus.integrated.bak* $N/10_Reports.bak* 2>/dev/null
