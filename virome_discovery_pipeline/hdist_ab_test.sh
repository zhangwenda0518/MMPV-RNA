#!/bin/bash
# hdist_ab_test.sh — BBDuk hdist=0 vs 1 A/B 实验 (速度 + 召回 + 末端闭合)
#
# 阶段1: 每个 contig 复刻 VSI 第一轮 tips, 3 种配置各跑一次全量收割:
#         A hd0      (生产现状基线)
#         B hd0+java (量化 -da -Xmx8g 的速度收益)
#         C hd1+java (目标配置: 速度 + 召回差异)
#         输出: 墙钟/bbduk自报耗时, 捞取对数, C 相对 A 的集合增减
# 阶段2: 3 个 contig 用 VSI_BBDUK_HDIST=1 完整重跑 (原 scaffold + 快照 reads),
#         最终 scaffold 与生产老结果 + 阶段1 hd0 产物 blast 对齐 → 末端闭合判定
export PATH=/home/zhangwenda/biosoft/ncbi-blast-2.13.0+/bin:/home/zhangwenda/.pixi/bin:/home/zhangwenda/mambaforge/bin:$PATH
export LC_ALL=C
BBDUK=/home/zhangwenda/mambaforge/bin/bbduk.sh
SNAP=/tmp/vsi_snap
W=/tmp/vsi_hdist
OLD_B=/home/zhangwenda/data-test/out/08_Rescue/Plant/branch_b
HD0_OUT=$SNAP/out
S=/home/zhangwenda/mambaforge/envs/Virseqimprover/bin/salmon
CK=/home/zhangwenda/database/virus-db/checkv-db-v1.7
CONTIGS=(
  ERR2040166_clean_NODE_132_length_4683_cov_23.698795
  ERR2040845_clean_NODE_1088_length_2950_cov_10.158239
  ERR2041000_clean_NODE_115_length_4426_cov_22.243669
)
rm -rf "$W"; mkdir -p "$W/ab" "$W/hd1_out"

avg_len() { pigz -dc "$1" | head -n 40000 | awk 'NR==1{fq=($0~/^@/)} fq{if(NR%4==2){s+=length($0);n++}} !fq{if($0~/^>/)n++;else s+=length($0)} END{printf "%.1f", s/n}'; }
names()   { awk '/^>/{sub(/^>/,""); print $1}' "$1" | sort -u; }
bp()      { grep -v '>' "$1" 2>/dev/null | tr -d '\n' | wc -c; }

# ── 阶段 1: A/B/C 收割对比 (按 contig 并行, contig 内串行) ──────────
for N in "${CONTIGS[@]}"; do
(
  D=$W/ab/$N; mkdir -p "$D"
  R1=$SNAP/${N}_merged_R1.fa.gz; R2=$SNAP/${N}_merged_R2.fa.gz
  cp "$SNAP/$N.fasta" "$D/ref.fasta"; samtools faidx "$D/ref.fasta" >/dev/null 2>&1
  AV=$(avg_len "$R1"); EDGE=$(awk -v a="$AV" 'BEGIN{print int(a*1.5+0.999)}')
  ID=$(cut -f1 "$D/ref.fasta.fai"); LEN=$(cut -f2 "$D/ref.fasta.fai")
  printf '%s\t0\t%s\n%s\t%s\t%s\n' "$ID" "$EDGE" "$ID" "$((LEN-EDGE))" "$LEN" > "$D/tips.bed"
  bedtools getfasta -fi "$D/ref.fasta" -bed "$D/tips.bed" -fo "$D/tips.fasta" 2>/dev/null
  echo "$N avg=$AV edge=$EDGE" >> "$W/ab/tips.log"

  run_cfg() { # tag hdist jvm(空格分隔, 可空)
    local TAG=$1 HD=$2 JVM=$3
    local t0=$(date +%s%3N)
    bash $BBDUK $JVM in="$R1" in2="$R2" \
        outm="$D/${TAG}_1.fa" outm2="$D/${TAG}_2.fa" \
        ref="$D/tips.fasta" k=25 hdist=$HD rcomp=t threads=12 overwrite=t prealloc=t \
        2> "$D/${TAG}.bbduk.log"
    local rc=$?
    local wall=$(( $(date +%s%3N) - t0 ))
    local self=$(awk -F'\t' '/Time:/{gsub(/ /,"",$2); sub(/s$/,"",$2); print $2; exit}' "$D/${TAG}.bbduk.log")
    echo -e "$TAG\t$HD\t$wall\t$self\t$(grep -c '^>' "$D/${TAG}_1.fa")" >> "$D/cfg.tsv"
    echo "rc=$rc" >> "$D/cfg.tsv"
  }
  run_cfg hd0    0 ""
  run_cfg hd0java 0 "-da -Xmx8g"
  run_cfg hd1java 1 "-da -Xmx8g"

  # 集合差异: hd1 相对 hd0
  names "$D/hd0_1.fa"    > "$D/n.hd0"
  names "$D/hd1java_1.fa" > "$D/n.hd1"
  {
    echo "hd0_pairs=$(wc -l < "$D/n.hd0") hd1_pairs=$(wc -l < "$D/n.hd1")"
    echo "hd1_extra=$(comm -13 "$D/n.hd0" "$D/n.hd1" | wc -l) hd0_only=$(comm -23 "$D/n.hd0" "$D/n.hd1" | wc -l)"
  } >> "$D/cfg.tsv"
) &
done
wait
echo "[ab] phase1 done $(date +%H:%M:%S)"

# ── 阶段 2: hdist=1 完整 VSI 重跑 ───────────────────────────────────
export VSI_BBDUK_HDIST=1
for N in "${CONTIGS[@]}"; do
(
  O=$W/hd1_out/$N; mkdir -p "$O"
  t0=$(date +%s)
  timeout 7200 python -u /tmp/vsi_wire/Virseqimprover.py \
      -1 "$SNAP/${N}_merged_R1.fa.gz" -2 "$SNAP/${N}_merged_R2.fa.gz" \
      -scaffold "$SNAP/$N.fasta" -o "$O" -salmon "$S" -t 10 -checkv_db "$CK" \
      > "$W/hd1_out/${N}.run.log" 2>&1
  echo "$? $(( $(date +%s) - t0 ))" > "$W/hd1_out/${N}.status"
) &
done
wait
unset VSI_BBDUK_HDIST
echo "[ab] phase2 done $(date +%H:%M:%S)"

# ── 汇总 ────────────────────────────────────────────────────────────
{
  echo "================ 阶段1: 速度 + 召回 ================"
  for N in "${CONTIGS[@]}"; do
    D=$W/ab/$N
    echo "--- ${N:8:14} ($(cat "$W/ab/tips.log" | grep "${N:0:20}" | head -1 | sed "s/.*avg=/avg=/"))"
    echo -e "配置\t墙钟ms\tbbduk自报s\t捞取对数"
    ( head -3 "$D/cfg.tsv" ) | column -t
    grep -E "hd0_pairs|hd1_extra" "$D/cfg.tsv"
  done
  echo ""
  echo "================ 阶段2: 末端闭合 (hd1 完整重跑 vs 生产老结果 vs hd0) ================"
  printf "%-14s %7s %7s %8s | %s\n" contig hd0_new hd1_new old_生产
  for N in "${CONTIGS[@]}"; do
    NEW0=$HD0_OUT/$N/scaffold.fasta
    NEW1=$W/hd1_out/$N/scaffold.fasta
    OLD=$OLD_B/out_$N/scaffold.fasta
    printf "%-14s %7s %7s %8s\n" "${N:8:14}" "$(bp "$NEW0")" "$(bp "$NEW1")" "$(bp "$OLD")"
    if [ -s "$NEW1" ] && [ -s "$OLD" ]; then
      L=$(/home/zhangwenda/biosoft/ncbi-blast-2.13.0+/bin/blastn -query "$NEW1" -subject "$OLD" -outfmt '6 pident length qlen' -max_target_seqs 1 2>/dev/null | head -1)
      echo "    hd1_vs_old: pident=$(echo "$L" | cut -d' ' -f1) aligned=$(echo "$L" | cut -d' ' -f2) (hd0基线: NODE_115=99.866/4482, NODE_132=100/4912, NODE_1088=100/2950)"
    fi
    L2=$(/home/zhangwenda/biosoft/ncbi-blast-2.13.0+/bin/blastn -query "$NEW1" -subject "$NEW0" -outfmt '6 pident length qlen' -max_target_seqs 1 2>/dev/null | head -1)
    [ -n "$L2" ] && echo "    hd1_vs_hd0: pident=$(echo "$L2" | cut -d' ' -f1) aligned=$(echo "$L2" | cut -d' ' -f2)"
    COMP=$(awk -F"\t" 'NR==1{for(i=1;i<=NF;i++) if($i=="aai_completeness") C=i} NR==2{printf "%.1f%%", $C}' "$W/hd1_out/$N/checkv_tmp/completeness.tsv" 2>/dev/null)
    read RC T < "$W/hd1_out/${N}.status" 2>/dev/null
    echo "    hd1_checkv=$COMP rc=$RC 耗时=$((T/60))min traceback=$(grep -c Traceback "$W/hd1_out/${N}.run.log" 2>/dev/null)"
  done
} > "$W/SUMMARY.txt" 2>&1
touch "$W/DONE"
echo "[ab] summary ready $(date +%H:%M:%S)"
