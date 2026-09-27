#!/bin/bash
# mx_ab_test.sh — -Xmx8g vs -Xmx20g A/B: 速度 + 实际峰值 RSS
# 复用 hdist 实验留下的 tips (/tmp/vsi_hdist/ab/*/tips.fasta), hd0 配置, 各 2 次重复
export PATH=/home/zhangwenda/biosoft/ncbi-blast-2.13.0+/bin:/home/zhangwenda/.pixi/bin:/home/zhangwenda/mambaforge/bin:$PATH
export LC_ALL=C
BBDUK=/home/zhangwenda/mambaforge/bin/bbduk.sh
H=/tmp/vsi_mx
rm -rf "$H"; mkdir -p "$H"

declare -A R1=(
  [NODE_132]=/tmp/vsi_snap/ERR2040166_clean_NODE_132_length_4683_cov_23.698795_merged_R1.fa.gz
  [NODE_1088]=/tmp/vsi_snap/ERR2040845_clean_NODE_1088_length_2950_cov_10.158239_merged_R1.fa.gz
)
declare -A R2=(
  [NODE_132]=/tmp/vsi_snap/ERR2040166_clean_NODE_132_length_4683_cov_23.698795_merged_R2.fa.gz
  [NODE_1088]=/tmp/vsi_snap/ERR2040845_clean_NODE_1088_length_2950_cov_10.158239_merged_R2.fa.gz
)
declare -A TIPS=(
  [NODE_132]=/tmp/vsi_hdist/ab/ERR2040166_clean_NODE_132_length_4683_cov_23.698795/tips.fasta
  [NODE_1088]=/tmp/vsi_hdist/ab/ERR2040845_clean_NODE_1088_length_2950_cov_10.158239/tips.fasta
)

for K in NODE_132 NODE_1088; do
(
  D=$H/$K; mkdir -p "$D"
  echo -e "rep\tXmx\t墙钟s\t自报s\t峰值RSS_GB\t对数" > "$D/res.tsv"
  for REP in 1 2; do
    for X in 8g 20g; do
      t0=$(date +%s%3N)
      /usr/bin/time -v bash $BBDUK -da -Xmx$X \
          in="${R1[$K]}" in2="${R2[$K]}" \
          outm="$D/o_${REP}_${X}_1.fa" outm2="$D/o_${REP}_${X}_2.fa" \
          ref="${TIPS[$K]}" k=25 hdist=0 rcomp=t threads=12 overwrite=t prealloc=t \
          2> "$D/t_${REP}_${X}.log"
      rc=$?
      wall=$(( $(date +%s%3N) - t0 ))
      self=$(awk -F'\t' '/Time:/{gsub(/ /,"",$2); sub(/s$/,"",$2); print $2; exit}' "$D/t_${REP}_${X}.log")
      rss=$(awk '/Maximum resident/{printf "%.2f", $6/1048576}' "$D/t_${REP}_${X}.log")
      pairs=$(grep -c '^>' "$D/o_${REP}_${X}_1.fa" 2>/dev/null)
      echo -e "$REP\t$X\t$((wall/1000))\t${self}\t${rss}\t${pairs}\t(rc=$rc)" >> "$D/res.tsv"
    done
  done
) &
done
wait

{
  echo "================ -Xmx8g vs 20g ================"
  for K in NODE_132 NODE_1088; do
    echo "--- $K"
    column -t "$H/$K/res.tsv"
  done
  echo ""
  echo "(机器物理内存约 126G; bbduk 默认自动抓 85%≈107G, 实测负载 190 下曾致 OOM)"
} > "$H/SUMMARY.txt" 2>&1
touch "$H/DONE"
echo done
