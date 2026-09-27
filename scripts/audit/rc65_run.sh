#!/bin/bash
# v6.5 三次实跑：count 模式正式 / count 模式闸门关闭 / weighted 模式回归对账
set -u
R=/tmp/virus_classifier_analysis.R.cascade_v65
echo "=== parse check ==="
Rscript -e "invisible(parse('$R')); cat('PARSE OK\n')" 2>&1 | tail -5
echo "=== patch md5 / lines ==="
md5sum $R; wc -l $R
INT=/home/zhangwenda/MMPV-paper/goji-virome/02_novel_virus/RNA-Lycium_barbarum_out/05_Taxonomy
COMB=$INT/Votus.classed/Votus_combined_taxonomy.tsv
ls -la $COMB
rm -rf /tmp/rc65_count /tmp/rc65_weighted /tmp/pregate65_count
mkdir -p /tmp/rc65_count /tmp/rc65_weighted /tmp/pregate65_count
cd /tmp

echo "=== run 1/3: v6.5 count（默认，正式版） ==="
start=$(date +%s)
MMPV_CONSENSUS_MODE=cascade Rscript $R --combined $COMB --output /tmp/rc65_count > /tmp/rc65_count.log 2>&1
echo "exit=$? elapsed=$(( $(date +%s) - start ))s"
grep -c . /tmp/rc65_count/final_integrated_classification.tsv
head -2 /tmp/rc65_count/tool_weights.tsv

echo "=== run 2/3: v6.5 count（闸门关闭） ==="
python3 /tmp/make_pregate_variant.py $R /tmp/virus_classifier_analysis.R.pregate65
MMPV_CONSENSUS_MODE=cascade Rscript /tmp/virus_classifier_analysis.R.pregate65 --combined $COMB --output /tmp/pregate65_count > /tmp/pregate65_count.log 2>&1
echo "exit=$? elapsed=$(( $(date +%s) - start ))s"

echo "=== run 3/3: v6.5 weighted（回退口径回归对账） ==="
start=$(date +%s)
MMPV_CONSENSUS_MODE=cascade MMPV_VOTE_WEIGHT=weighted Rscript $R --combined $COMB --output /tmp/rc65_weighted > /tmp/rc65_weighted.log 2>&1
echo "exit=$? elapsed=$(( $(date +%s) - start ))s"

echo "=== weighted 回归：v65(mode=weighted) vs v64 正式产物 ==="
md5sum /tmp/rc65_weighted/final_integrated_classification.tsv /tmp/rc64b_cascade/final_integrated_classification.tsv
md5sum /tmp/rc65_weighted/tool_weights.tsv /tmp/rc64b_cascade/tool_weights.tsv

echo "=== count vs weighted：逐格差异 ==="
python3 - <<'PYEOF'
import csv, os
def load(p):
    with open(p, newline='', encoding='utf-8') as f:
        rd = csv.DictReader(f, delimiter='\t')
        return {r['contig_id']: r for r in rd}, rd.fieldnames
A = '/tmp/rc65_count/final_integrated_classification.tsv'
B = '/tmp/rc65_weighted/final_integrated_classification.tsv'
a, fa = load(A); b, fb = load(B)
print('行数 count=%d weighted=%d 列数 %d/%d 列名一致=%s' % (len(a), len(b), len(fa), len(fb), fa == fb))
for lv in ['Realm','Kingdom','Phylum','Class','Order','Family','Genus','Species']:
    n = sum(1 for c in a if a[c].get(lv,'') != b.get(c,{}).get(lv,''))
    print('  %-9s 差异 %5d 行' % (lv, n))
cell = sum(1 for c in a for lv in ['Realm','Kingdom','Phylum','Class','Order','Family','Genus','Species'] if a[c].get(lv,'') != b.get(c,{}).get(lv,''))
print('  合计差异格次 %d' % cell)
PYEOF
