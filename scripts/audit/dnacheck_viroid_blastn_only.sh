#!/bin/bash
# 枸杞三基因组 · 类病毒 word11 分片 blastn + 过滤（分片已就绪）
set -uo pipefail
V=/home/zhangwenda/database/virus-db/viroids-db/viroids.fasta
T=/home/zhangwenda/eve_screen_goji3_20260929/viroid_word11
O=/home/zhangwenda/eve_screen_goji3_20260929

for sp in zhonghua ningxia heiguo; do
  ls "$T/$sp"/chr*.fa 2>/dev/null | xargs -P 48 -I{} bash -c '
    s="{}"; b="${s%.fa}.blastn.tsv"
    blastn -task blastn -word_size 11 -dust no -query '"$V"' -subject "$s" -evalue 1e-4 \
      -outfmt "6 qseqid qlen sseqid slen pident length qstart qend sstart send evalue bitscore qcovs" \
      -out "$b" 2>/dev/null || true
  '
  cat "$T/$sp"/chr*.blastn.tsv 2>/dev/null | awk -F'\t' '$6>=80 && $13>=50' > "$O/${sp}_viroid_word11.tsv"
  raw=$(cat "$T/$sp"/chr*.blastn.tsv 2>/dev/null | wc -l)
  filt=$(wc -l < "$O/${sp}_viroid_word11.tsv")
  echo "[$sp] raw=$raw filtered=$filt"
done
echo VIROID_ALL_DONE
