#!/bin/bash
# 裸 bowtie2 vs --local 速度对比
set -e

THREADS=10
BT2=/home/zhangwenda/biosoft/bowtie2/bin/bowtie2
BT2BUILD=/home/zhangwenda/biosoft/bowtie2/bin/bowtie2-build
WORKDIR=/tmp/bt2_local_test
mkdir -p $WORKDIR && cd $WORKDIR

CONTIG=/home/zhangwenda/data-test/out/08_Rescue/Plant/branch_b/scaffolds/ERR2040118_clean_NODE_150_length_3726_cov_29.718708.fasta
R1=/home/zhangwenda/data-test/out/00b_HostDepletion/ERR2040118_clean_1.fa.gz
R2=/home/zhangwenda/data-test/out/00b_HostDepletion/ERR2040118_clean_2.fa.gz
cp $CONTIG ref.fasta

echo "=== build index (once) ==="
$BT2BUILD --threads $THREADS ref.fasta bt2idx > /dev/null 2>&1

: > results.txt
timed() {
  local NAME="$1"; shift
  T0=$(date +%s.%N)
  eval "$1" 2> ${NAME}.err | samtools sort -@ $THREADS - -o ${NAME}.bam
  T1=$(date +%s.%N)
  local T=$(echo "$T1 - $T0" | bc)
  local REC=$(samtools view -c ${NAME}.bam)
  local DSUM=$(samtools depth -a ${NAME}.bam | awk '{s+=$3} END{print s}')
  echo "$NAME|$T|$REC|$DSUM" >> results.txt
  echo "[$NAME] time=${T}s records=$REC depth=$DSUM"
}

echo "=== 1) 裸 bowtie2 (原版形态: 无 --local) ==="
timed raw "$BT2 -x bt2idx -1 $R1 -2 $R2"

echo "=== 2) --no-unal only ==="
timed nounal "$BT2 --no-unal -x bt2idx -1 $R1 -2 $R2"

echo "=== 3) --local --no-unal (当前部署) ==="
timed local "$BT2 --local --no-unal -x bt2idx -1 $R1 -2 $R2"

echo ""
echo "############ 纯比对器 A/B (90bp reads, 3726bp contig, 10 threads) ############"
printf "%-18s %10s %12s %14s\n" "mode" "time(s)" "records" "depth_sum"
printf "%-18s %10s %12s %14s\n" "----" "-------" "-------" "---------"
while IFS='|' read n t r d; do
  printf "%-18s %10.2f %12s %14s\n" "$n" "$t" "$r" "$d"
done < results.txt
echo ""
echo "=== depth 一致性 ==="
for n in raw nounal local; do echo "$n: $(md5sum $n.bam 2>/dev/null | cut -d' ' -f1) (bam md5)"; done
