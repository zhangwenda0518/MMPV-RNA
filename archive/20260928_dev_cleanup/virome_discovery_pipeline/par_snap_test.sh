#!/bin/bash
# par_snap_test.sh — 快照生产正在处理的 contig 输入, 并行跑新代码 (bbduk), 免疫 rescue 清理竞态
export PATH=/home/zhangwenda/biosoft/ncbi-blast-2.13.0+/bin:/home/zhangwenda/.pixi/bin:/home/zhangwenda/mambaforge/bin:$PATH
export LC_ALL=C
B=/home/zhangwenda/data-test/out/08_Rescue/Plant/branch_b
S=/home/zhangwenda/mambaforge/envs/Virseqimprover/bin/salmon
CK=/home/zhangwenda/database/virus-db/checkv-db-v1.7
SNAP=/tmp/vsi_snap
P=/tmp/vsi_snap/out
rm -rf "$SNAP"; mkdir -p "$SNAP" "$P"

CONTIGS=(
  ERR2040166_clean_NODE_132_length_4683_cov_23.698795
  ERR2040845_clean_NODE_1088_length_2950_cov_10.158239
  ERR2041000_clean_NODE_115_length_4426_cov_22.243669
)

# ── 快照输入 (rescue 会删 merged_reads) ──
for N in "${CONTIGS[@]}"; do
  if [ -f "$B/merged_reads/${N}_merged_R1.fa.gz" ] && [ -f "$B/scaffolds/$N.fasta" ]; then
    cp "$B/scaffolds/$N.fasta" "$SNAP/$N.fasta" &
    cp "$B/merged_reads/${N}_merged_R1.fa.gz" "$SNAP/${N}_merged_R1.fa.gz" &
    cp "$B/merged_reads/${N}_merged_R2.fa.gz" "$SNAP/${N}_merged_R2.fa.gz" &
  else
    echo "SKIP $N (input vanished before snapshot)"
  fi
done
wait
echo "[snap] snapshot done $(date +%H:%M:%S)"

# ── 并行跑新代码 ──
for N in "${CONTIGS[@]}"; do
  [ -f "$SNAP/${N}_merged_R1.fa.gz" ] || { echo "$N: no snapshot, skip" >> "$P/skip.txt"; continue; }
  (
    O=$P/$N
    mkdir -p "$O"
    t0=$(date +%s)
    timeout 7200 python -u /tmp/vsi_wire/Virseqimprover.py \
        -1 "$SNAP/${N}_merged_R1.fa.gz" \
        -2 "$SNAP/${N}_merged_R2.fa.gz" \
        -scaffold "$SNAP/$N.fasta" \
        -o "$O" -salmon "$S" -t 10 -checkv_db "$CK" \
        > "$P/${N}.run.log" 2>&1
    echo "$? $(( $(date +%s) - t0 ))" > "$P/${N}.status"
  ) &
done
wait
echo "[par] all VSI done $(date +%H:%M:%S)"

# ── 汇总 ──
bp() { grep -v '>' "$1" 2>/dev/null | tr -d '\n' | wc -c; }
{
  printf "%-14s %8s %10s %8s %9s %8s %8s\n" contig orig new 增长 checkv_comp new耗时 rc
  for N in "${CONTIGS[@]}"; do
    [ -f "$SNAP/${N}_merged_R1.fa.gz" ] || continue
    read RC T < "$P/${N}.status" 2>/dev/null
    O_=$(bp "$SNAP/$N.fasta"); W_=$(bp "$P/$N/scaffold.fasta")
    COMP=$(awk -F"\t" 'NR==1{for(i=1;i<=NF;i++) if($i=="aai_completeness") C=i} NR==2{printf "%.1f%%", $C}' "$P/$N/checkv_tmp/completeness.tsv" 2>/dev/null)
    GROW=$(awk -v w=$W_ -v o=$O_ 'BEGIN{print w-o}')
    printf "%-14s %8s %10s %8s %9s %8s %8s\n" \
      "${N:8:14}" "$O_" "${W_:-0}" "$GROW" "${COMP:-NA}" "$((T/60))min" "$RC"
  done
  echo ""
  echo "== traceback 计数 =="
  for N in "${CONTIGS[@]}"; do
    printf "%s: %s\n" "${N:8:14}" "$(grep -c Traceback "$P/${N}.run.log" 2>/dev/null)"
  done
} > "$P/SUMMARY.txt" 2>&1
touch "$P/DONE"
echo "[par] summary ready $(date +%H:%M:%S)"
