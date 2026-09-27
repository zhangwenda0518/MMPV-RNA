#!/bin/bash
# 只读：integrated_summary.tsv 的 genome_type 字段实际填充率
G=~/MMPV-paper/goji-virome/02_novel_virus/RNA-Lycium_barbarum_out
O=~/MMPV-paper/onekp-virome/onekp-virus

echo "=== [1] 定位 integrated_summary.tsv ==="
find $G -maxdepth 4 -name "integrated_summary.tsv" 2>/dev/null
find $O -maxdepth 4 -name "integrated_summary.tsv" 2>/dev/null
echo "--- 原始 miuvig_taxonomy.tsv ---"
find $G -maxdepth 5 -name "miuvig_taxonomy.tsv" 2>/dev/null
find $O -maxdepth 5 -name "miuvig_taxonomy.tsv" 2>/dev/null

echo
echo "=== [2] 逐文件统计 genome_type 填充 ==="
for f in $(find $G -maxdepth 4 -name "integrated_summary.tsv" 2>/dev/null) \
         $(find $O -maxdepth 4 -name "integrated_summary.tsv" 2>/dev/null) \
         $(find $G -maxdepth 5 -name "miuvig_taxonomy.tsv" 2>/dev/null) \
         $(find $O -maxdepth 5 -name "miuvig_taxonomy.tsv" 2>/dev/null); do
  echo "--- $f"
  stat -c '    mtime %y  size %s' "$f"
  head -1 "$f"
  python3 - "$f" <<'PY'
import csv, sys, collections
p = sys.argv[1]
with open(p, newline='', encoding='utf-8', errors='replace') as f:
    r = csv.reader(f, delimiter='\t', quotechar='"')
    hdr = [h.strip() for h in next(r)]
    idx = {}
    for k in ('genome_type', 'pred_genome_type', 'genome_struc', 'pred_genome_struc'):
        if k in hdr:
            idx[k] = hdr.index(k)
    if not idx:
        print('    无 genome_type 类列; 表头前 8 列 =', hdr[:8])
        sys.exit()
    n = 0
    cnt = {k: collections.Counter() for k in idx}
    for row in r:
        n += 1
        for k, i in idx.items():
            cnt[k][row[i].strip() if i < len(row) else ''] += 1
    for k in idx:
        nonempty = sum(v for kk, v in cnt[k].items() if kk not in ('', 'NA', 'N/A'))
        print('    %s: 总 %d 行, 非空 %d' % (k, n, nonempty))
        for kk, v in cnt[k].most_common(8):
            print('        %-12s %d' % (repr(kk)[:12], v))
PY
done
