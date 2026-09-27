#!/bin/bash
# 3-way speed benchmark: bowtie2 vs minibwa vs strobealign
# Unified protocol: mapped-only output, samtools sort -> BAM -> depth
set -e

THREADS=10
WORKDIR=/tmp/mbwa_ab_test
BT2=/home/zhangwenda/biosoft/bowtie2/bin/bowtie2
BT2BUILD=/home/zhangwenda/biosoft/bowtie2/bin/bowtie2-build
MBWA=/home/zhangwenda/biosoft/minibwa/minibwa
STROBE=/home/zhangwenda/mambaforge/bin/strobealign

mkdir -p $WORKDIR
cd $WORKDIR

CONTIG=/tmp/bt2_ab_test/ERR2041160_clean_NODE_362_length_3154_cov_1116.864674.fasta
R1=/home/zhangwenda/data-test/out/00b_HostDepletion/ERR2041160_clean_1.fa.gz
R2=/home/zhangwenda/data-test/out/00b_HostDepletion/ERR2041160_clean_2.fa.gz

cp $CONTIG ref.fasta
echo "=== ref: $(stat -c%s ref.fasta) bytes; reads: $(zcat $R1 | head -4 | tail -1 | awk '{print length($0)}')bp-ish ==="
echo ""

run_one() {
  local NAME=$1; shift
  local ODIR=${NAME}_out
  rm -rf $ODIR && mkdir -p $ODIR
  cd $ODIR; cp ../ref.fasta .

  # index
  T0=$(date +%s.%N)
  eval "$IDX_CMD" > idx.log 2>&1 || { echo "$NAME IDX FAIL"; tail -3 idx.log; cd ..; return; }
  T1=$(date +%s.%N)
  local IDX_T=$(echo "$T1 - $T0" | bc)

  # align + sort
  T0=$(date +%s.%N)
  eval "$ALN_CMD" 2> aln.log | samtools sort -@ $THREADS - -o out.bam
  T1=$(date +%s.%N)
  local ALN_T=$(echo "$T1 - $T0" | bc)

  samtools index out.bam 2>/dev/null || true
  samtools depth -a out.bam > depth.txt
  local REC=$(samtools view -c out.bam)
  local DSUM=$(awk '{s+=$3} END{print s}' depth.txt)

  echo "$NAME|$IDX_T|$ALN_T|$REC|$DSUM" >> ../bench_results.txt
  echo "[$NAME] idx=${IDX_T}s aln=${ALN_T}s rec=$REC depth=$DSUM"
  cd ..
}

: > bench_results.txt

echo "=== [1/3] bowtie2 ==="
IDX_CMD='$BT2BUILD --threads $THREADS ref.fasta bt2idx'
ALN_CMD='$BT2 --local --no-unal --threads $THREADS -f -x bt2idx -1 $R1 -2 $R2'
run_one bowtie2

echo "=== [2/3] minibwa ==="
IDX_CMD='$MBWA index -t $THREADS ref.fasta'
ALN_CMD='$MBWA map -u -t $THREADS ref.fasta $R1 $R2'
run_one minibwa

echo "=== [3/3] strobealign ==="
IDX_CMD='$STROBE -t $THREADS -i ref.fasta $R1 $R2'
ALN_CMD='$STROBE -U -t $THREADS --use-index ref.fasta $R1 $R2'
run_one strobealign

echo ""
echo "############ RESULT ############"
printf "%-14s %10s %12s %14s %14s\n" "tool" "idx(s)" "aln(s)" "BAM_records" "depth_sum"
printf "%-14s %10s %12s %14s %14s\n" "----" "------" "------" "-----------" "---------"
while IFS='|' read n i a r d; do
  printf "%-14s %10.2f %12.2f %14s %14s\n" "$n" "$i" "$a" "$r" "$d"
done < bench_results.txt
echo ""
echo "=== depth 一致性 md5 ==="
for n in bowtie2 minibwa strobealign; do
  echo "$n: $(md5sum ${n}_out/depth.txt 2>/dev/null | cut -d' ' -f1)"
done
