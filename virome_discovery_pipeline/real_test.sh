#!/bin/bash
# real_test.sh — 真实样本端到端测试: 新代码(bbduk默认) vs 生产老代码(salmon)刚产出的结果
# 对比标准 (archive/VSI优化方案.md 4.1): 长度差 <1% 且 identity >99%
export PATH=/home/zhangwenda/biosoft/ncbi-blast-2.13.0+/bin:/home/zhangwenda/.pixi/bin:/home/zhangwenda/mambaforge/bin:$PATH
B=/home/zhangwenda/data-test/out/08_Rescue/Plant/branch_b
N=ERR2041195_clean_NODE_1263_length_3739_cov_36.666126
S=/home/zhangwenda/mambaforge/envs/Virseqimprover/bin/salmon
CK=/home/zhangwenda/database/virus-db/checkv-db-v1.7
O=/tmp/vsi_real/out_bbduk
rm -rf /tmp/vsi_real
mkdir -p /tmp/vsi_real

echo "[real] start $(date +%H:%M:%S)"
timeout 7200 python -u /tmp/vsi_wire/Virseqimprover.py \
    -1 "$B/merged_reads/${N}_merged_R1.fa.gz" \
    -2 "$B/merged_reads/${N}_merged_R2.fa.gz" \
    -scaffold "$B/scaffolds/$N.fasta" \
    -o "$O" -salmon "$S" -t 10 -checkv_db "$CK" \
    > /tmp/vsi_real/run.log 2>&1
echo "vsi_rc=$?" >> /tmp/vsi_real/run.log
echo "[real] vsi done $(date +%H:%M:%S)"

NEW=$O/scaffold.fasta
OLD=$B/out_$N/scaffold.fasta
ORIG=$B/scaffolds/$N.fasta
{
  echo "== 长度对比 (bp) =="
  echo "orig: $(grep -v '>' "$ORIG" | tr -d '\n' | wc -c)"
  echo "old(生产salmon): $(grep -v '>' "$OLD" | tr -d '\n' | wc -c)"
  echo "new(本次bbduk):  $(grep -v '>' "$NEW" | tr -d '\n' | wc -c)"
  echo ""
  echo "== blastn 新(query) vs 老(subject) =="
  makeblastdb -in "$OLD" -dbtype nucl > /dev/null 2>&1
  blastn -query "$NEW" -db "$OLD" -outfmt '6 pident length qlen slen' -max_target_seqs 1 | head -3
  echo ""
  echo "== bbduk 收割统计 (每轮) =="
  grep "\[bbduk\]" /tmp/vsi_real/run.log | grep -E "Input:|Time:" | head -10
  echo ""
  echo "== 迭代轨迹 =="
  grep -a "Trying to grow\|CheckV\|VSI finished" /tmp/vsi_real/run.log "$O/output-log.txt" 2>/dev/null | tail -6
} > /tmp/vsi_real/compare.txt 2>&1
touch /tmp/vsi_real/DONE
echo "[real] all done $(date +%H:%M:%S)"
