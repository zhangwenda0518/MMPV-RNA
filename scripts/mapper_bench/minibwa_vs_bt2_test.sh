#!/bin/bash
# A/B test: minibwa vs bowtie2 on VSI-style alignment task
# Data: same contig + reads as the earlier --no-unal A/B test
set -e

THREADS=10
WORKDIR=/tmp/mbwa_ab_test
BT2=/home/zhangwenda/biosoft/bowtie2/bin/bowtie2
BT2BUILD=/home/zhangwenda/biosoft/bowtie2/bin/bowtie2-build
MBWA=/home/zhangwenda/biosoft/minibwa/minibwa

mkdir -p $WORKDIR
cd $WORKDIR

CONTIG=/tmp/bt2_ab_test/ERR2041160_clean_NODE_362_length_3154_cov_1116.864674.fasta
R1=/home/zhangwenda/data-test/out/00b_HostDepletion/ERR2041160_clean_1.fa.gz
R2=/home/zhangwenda/data-test/out/00b_HostDepletion/ERR2041160_clean_2.fa.gz

for f in $CONTIG $R1 $R2; do
  [ -f "$f" ] || { echo "MISSING: $f"; exit 1; }
done

cp $CONTIG ref.fasta
echo "=== 数据 ==="
echo "contig: $(grep -c '^>' ref.fasta) seq, $(stat -c%s ref.fasta) bytes"
echo "reads: $(zcat $R1 | head -400000 | awk 'END{print NR/4}') (前400k行约计数)"
echo ""

########################
# BOWTIE2
########################
echo "=== [A] bowtie2 ==="
rm -rf bt2_out && mkdir bt2_out
cd bt2_out
cp ../ref.fasta .

T0=$(date +%s.%N)
$BT2BUILD --threads $THREADS ref.fasta bt2idx > build.log 2>&1
T1=$(date +%s.%N)
BT2_BUILD=$(echo "$T1 - $T0" | bc)

T0=$(date +%s.%N)
$BT2 --local --no-unal --threads $THREADS -f -x bt2idx -1 $R1 -2 $R2 2> aln.log \
  | samtools sort -@ $THREADS - -o bt2.bam
T1=$(date +%s.%N)
BT2_ALN=$(echo "$T1 - $T0" | bc)

samtools index bt2.bam
samtools depth -a bt2.bam > bt2.depth.txt
BT2_REC=$(samtools view -c bt2.bam)
BT2_DEPTHSUM=$(awk '{s+=$3} END{print s}' bt2.depth.txt)
cd ..

########################
# MINIBWA
########################
echo "=== [B] minibwa ==="
rm -rf mbwa_out && mkdir mbwa_out
cd mbwa_out
cp ../ref.fasta .

T0=$(date +%s.%N)
$MBWA index -t $THREADS ref.fasta > build.log 2>&1
T1=$(date +%s.%N)
MBWA_BUILD=$(echo "$T1 - $T0" | bc)

T0=$(date +%s.%N)
$MBWA map -t $THREADS ref.fasta $R1 $R2 2> aln.log \
  | samtools sort -@ $THREADS - -o mbwa.bam
T1=$(date +%s.%N)
MBWA_ALN=$(echo "$T1 - $T0" | bc)

samtools index mbwa.bam 2>/dev/null || true
samtools depth -a mbwa.bam > mbwa.depth.txt
MBWA_REC=$(samtools view -c mbwa.bam)
MBWA_DEPTHSUM=$(awk '{s+=$3} END{print s}' mbwa.depth.txt)
cd ..

########################
# REPORT
########################
echo ""
echo "############ RESULT ############"
printf "%-28s %12s %12s\n" "metric" "bowtie2" "minibwa"
printf "%-28s %12s %12s\n" "----" "----" "----"
printf "%-28s %12.2f %12.2f\n" "index time (s)" "$BT2_BUILD" "$MBWA_BUILD"
printf "%-28s %12.2f %12.2f\n" "align+sort time (s)" "$BT2_ALN" "$MBWA_ALN"
printf "%-28s %12s %12s\n" "BAM records" "$BT2_REC" "$MBWA_REC"
printf "%-28s %12s %12s\n" "depth sum" "$BT2_DEPTHSUM" "$MBWA_DEPTHSUM"
echo ""
echo "speedup (align): $(echo "scale=2; $BT2_ALN / $MBWA_ALN" | bc)x"
echo ""
echo "=== bowtie2 aln.log tail ==="
tail -3 bt2_out/aln.log
echo "=== minibwa aln.log tail ==="
tail -3 mbwa_out/aln.log
