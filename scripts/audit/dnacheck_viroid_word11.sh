#!/bin/bash
# 枸杞三基因组 · 类病毒 word11 分片快筛（Lycium_EVE 协议：word 11 + dust no，
# 过滤线 aln>=80nt 且 query 覆盖 qcovs>=50%；word 7 精查只对快筛阳性 ±10kb 后补）
# 用法: bash viroid_word11.sh <OUTDIR>
set -euo pipefail
OUT=${1:-/home/zhangwenda/eve_screen_goji3_20260929}
VIROIDS=/home/zhangwenda/database/virus-db/viroids-db/viroids.fasta
TMP=$OUT/viroid_word11
mkdir -p "$TMP"

split_genome() {  # $1=fasta $2=outdir  按序列分片
  python3 - "$1" "$2" << 'PYEOF'
import sys, os
fa, out = sys.argv[1], sys.argv[2]
cur, buf, n = None, "", 0
def flush():
    if cur:
        with open(os.path.join(out, "chr%03d.fa" % n), "w") as w:
            w.write(">" + cur + "\n")
            for i in range(0, len(buf), 60):
                w.write(buf[i:i+60] + "\n")
for line in open(fa):
    if line.startswith(">"):
        flush()
        cur = line[1:].strip().split()[0]
        buf = ""
    else:
        buf += line.strip()
flush()
print("%s: %d shards" % (fa, n))
PYEOF
}

for sp in zhonghua ningxia heiguo; do
  case $sp in
    zhonghua) g=/home/zhangwenda/Lycium_EVE/00_edta/genomes/zhonghua.fa ;;
    ningxia)  g=/home/zhangwenda/Lycium_EVE/00_edta/genomes/ningxia.clean.fa ;;
    heiguo)   g=/home/zhangwenda/Lycium_EVE/00_edta/genomes/heiguo.clean.fa ;;
  esac
  mkdir -p "$TMP/$sp"
  split_genome "$g" "$TMP/$sp"
  ls "$TMP/$sp"/chr*.fa | xargs -P 13 -I{} bash -c '
    s="{}"; b="${s%.fa}.blastn.tsv"
    blastn -task blastn -word_size 11 -dust no \
      -query '"$VIROIDS"' -subject "$s" -evalue 1e-4 \
      -outfmt "6 qseqid qlen sseqid slen pident length qstart qend sstart send evalue bitscore qcovs" \
      -out "$b" 2>/dev/null || true
  '
  cat "$TMP/$sp"/chr*.blastn.tsv 2>/dev/null | awk -F'\t' '$6>=80 && $13>=50' > "$OUT/${sp}_viroid_word11.tsv" || true
  raw=$(cat "$TMP/$sp"/chr*.blastn.tsv 2>/dev/null | wc -l)
  filt=$(wc -l < "$OUT/${sp}_viroid_word11.tsv" || echo 0)
  echo "[$sp] raw=$raw filtered=$filt"
done
echo ALL_DONE
