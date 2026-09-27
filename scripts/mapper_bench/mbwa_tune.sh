#!/bin/bash
# minibwa 调参逼近 bowtie2 --local 召回的实验
set -e

THREADS=10
MBWA=/home/zhangwenda/bin/minibwa
WORKDIR=/tmp/mbwa_tune
mkdir -p $WORKDIR && cd $WORKDIR

CONTIG=/home/zhangwenda/data-test/out/08_Rescue/Plant/branch_b/scaffolds/ERR2040118_clean_NODE_150_length_3726_cov_29.718708.fasta
R1=/home/zhangwenda/data-test/out/00b_HostDepletion/ERR2040118_clean_1.fa.gz
R2=/home/zhangwenda/data-test/out/00b_HostDepletion/ERR2040118_clean_2.fa.gz
cp $CONTIG ref.fasta

[ -f ref.fasta.mbw ] || $MBWA index -t $THREADS ref.fasta > /dev/null 2>&1

# bowtie2 基线 (上一轮已测: 3193 records, depth 282212)
echo "=== baseline: bowtie2 --local --no-unal = 3193 records / depth 282212 ==="
echo ""

run() {
  local LABEL="$1"; shift
  T0=$(date +%s.%N)
  $MBWA map -u -t $THREADS "$@" ref.fasta $R1 $R2 2> ${LABEL}.err | samtools sort -@ $THREADS - -o ${LABEL}.bam
  T1=$(date +%s.%N)
  local T=$(echo "$T1 - $T0" | bc)
  local REC=$(samtools view -c ${LABEL}.bam)
  local DSUM=$(samtools depth -a ${LABEL}.bam | awk '{s+=$3} END{print s}')
  printf "%-30s %8.2fs  rec=%-8s depth=%s\n" "$LABEL" "$T" "$REC" "$DSUM"
}

echo "--- 默认 ---"
run default

echo "--- 降 -s (DP 分数阈值) ---"
run s10 -s 10
run s5  -s 5

echo "--- 降 -s 和 -m (链分数) ---"
run s10m10 -s 10 -m 10
run s5m5   -s 5 -m 5

echo "--- 降 -k 种子长度 ---"
run s5k15  -s 5 -m 5 -k 15

echo "--- 加宽带宽 ---"
run s5w200 -s 5 -m 5 -w 200

echo "--- 组合激进 ---"
run aggr -s 5 -m 5 -k 15 -w 200 -c 500

echo ""
echo "=== 目标: bowtie2 records=3193 depth=282212 ==="
