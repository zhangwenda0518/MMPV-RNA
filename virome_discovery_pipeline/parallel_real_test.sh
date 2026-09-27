#!/bin/bash
# parallel_real_test.sh — 3 个真实样本并行测试 (bbduk 默认路径) vs 生产老结果
export PATH=/home/zhangwenda/biosoft/ncbi-blast-2.13.0+/bin:/home/zhangwenda/.pixi/bin:/home/zhangwenda/mambaforge/bin:$PATH
export LC_ALL=C
B=/home/zhangwenda/data-test/out/08_Rescue/Plant/branch_b
S=/home/zhangwenda/mambaforge/envs/Virseqimprover/bin/salmon
CK=/home/zhangwenda/database/virus-db/checkv-db-v1.7
P=/tmp/vsi_par
rm -rf "$P"; mkdir -p "$P"

CONTIGS=(
  ERR2040798_clean_NODE_118_length_4894_cov_32.897753
  ERR2040900_clean_NODE_876_length_3733_cov_19.598374
  ERR2041191_clean_NODE_1256_length_3327_cov_9.079781
)

bp() { grep -v '>' "$1" 2>/dev/null | tr -d '\n' | wc -c; }

echo "[par] launch $(date +%H:%M:%S): ${#CONTIGS[@]} jobs, each -t 10"
for N in "${CONTIGS[@]}"; do
  (
    O=$P/out_$N
    mkdir -p "$O"
    t0=$(date +%s)
    timeout 7200 python -u /tmp/vsi_wire/Virseqimprover.py \
        -1 "$B/merged_reads/${N}_merged_R1.fa.gz" \
        -2 "$B/merged_reads/${N}_merged_R2.fa.gz" \
        -scaffold "$B/scaffolds/$N.fasta" \
        -o "$O" -salmon "$S" -t 10 -checkv_db "$CK" \
        > "$P/${N}.run.log" 2>&1
    echo "$? $(( $(date +%s) - t0 ))" > "$P/${N}.status"
  ) &
done
wait
echo "[par] all VSI done $(date +%H:%M:%S)"

# ── 汇总 ──
{
  printf "%-12s %8s %8s %8s %7s %9s %8s %8s\n" contig orig old new diff% identity new耗时 rc
  for N in "${CONTIGS[@]}"; do
    NEW=$P/out_$N/scaffold.fasta
    OLD=$B/out_$N/scaffold.fasta
    read RC T < "$P/${N}.status" 2>/dev/null
    O_=$(bp "$B/scaffolds/$N.fasta"); D_=$(bp "$OLD"); W_=$(bp "$NEW")
    if [ -s "$NEW" ] && [ -s "$OLD" ]; then
      LINE=$(blastn -query "$NEW" -subject "$OLD" -outfmt '6 pident length qlen' -max_target_seqs 1 2>/dev/null | head -1)
      PID=$(echo "$LINE" | cut -f1); QLEN=$(echo "$LINE" | cut -f3)
    else
      PID=NA; QLEN=NA
    fi
    DIFF=$(awk -v w=$W_ -v d=$D_ 'BEGIN{if(d>0) printf "%.1f", 100*(w-d)/d; else print "NA"}')
    COMP=$(awk -F"\t" 'NR==1{for(i=1;i<=NF;i++) if($i=="aai_completeness") C=i} NR==2{print int($C)"%"}' "$P/out_$N/checkv_tmp/completeness.tsv" 2>/dev/null)
    printf "%-12s %8s %8s %8s %7s %9s %8s %8s\n" \
      "${N:0:12}" "$O_" "$D_" "${W_:-0}" "$DIFF" "${PID:-NA}" "$((T/60))min" "$RC"
    echo "    checkv_aai_comp=$COMP aligned=$QLEN bp"
  done
} > "$P/SUMMARY.txt" 2>&1
touch "$P/DONE"
echo "[par] summary ready $(date +%H:%M:%S)"
