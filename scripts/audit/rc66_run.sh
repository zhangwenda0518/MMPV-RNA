#!/bin/bash
# v6.6 三次实跑：count 显式两级平票裁决 / weighted 回归 / pregate（闸门关闭）
set -u
R=/tmp/virus_classifier_analysis.R.cascade_v66
echo "=== parse check ==="
Rscript -e "invisible(parse('$R')); cat('PARSE OK\n')" 2>&1 | tail -5
echo "=== patch md5 / lines ==="
md5sum $R; wc -l $R
INT=/home/zhangwenda/MMPV-paper/goji-virome/02_novel_virus/RNA-Lycium_barbarum_out/05_Taxonomy
COMB=$INT/Votus.classed/Votus_combined_taxonomy.tsv
ls -la $COMB
rm -rf /tmp/rc66_count /tmp/rc66_weighted /tmp/pregate66_count
mkdir -p /tmp/rc66_count /tmp/rc66_weighted /tmp/pregate66_count
cd /tmp

echo "=== run 1/3: v6.6 count（默认，正式版） ==="
start=$(date +%s)
MMPV_CONSENSUS_MODE=cascade Rscript $R --combined $COMB --output /tmp/rc66_count > /tmp/rc66_count.log 2>&1
echo "exit=$? elapsed=$(( $(date +%s) - start ))s"

echo "=== run 2/3: v6.6 weighted（回退口径回归对账） ==="
start=$(date +%s)
MMPV_CONSENSUS_MODE=cascade MMPV_VOTE_WEIGHT=weighted Rscript $R --combined $COMB --output /tmp/rc66_weighted > /tmp/rc66_weighted.log 2>&1
echo "exit=$? elapsed=$(( $(date +%s) - start ))s"

echo "=== run 3/3: v6.6 count（闸门关闭，给诊断用） ==="
python3 /tmp/make_pregate_variant.py $R /tmp/virus_classifier_analysis.R.pregate66
MMPV_CONSENSUS_MODE=cascade Rscript /tmp/virus_classifier_analysis.R.pregate66 --combined $COMB --output /tmp/pregate66_count > /tmp/pregate66_count.log 2>&1
echo "exit=$?"

echo "=== 关键日志行（count） ==="
grep -aE "逐级淘汰|校准|置空|相容性约束|自检|计权模式" /tmp/rc66_count.log

echo "=== weighted 回归：v66(mode=weighted) vs v64 正式产物 ==="
md5sum /tmp/rc66_weighted/final_integrated_classification.tsv /tmp/rc64b_cascade/final_integrated_classification.tsv
md5sum /tmp/rc66_weighted/tool_weights.tsv /tmp/rc64b_cascade/tool_weights.tsv

echo "=== v65（隐式行序）vs v66（显式两级裁决）：成品逐格差异 ==="
python3 - <<'PYEOF'
import csv
def load(p):
    with open(p, newline='', encoding='utf-8') as f:
        rd = csv.DictReader(f, delimiter='\t'); return {r['contig_id']: r for r in rd}, rd.fieldnames
LV = ['Realm','Kingdom','Phylum','Class','Order','Family','Genus','Species']
a, fa = load('/tmp/rc65_count/final_integrated_classification.tsv')
b, fb = load('/tmp/rc66_count/final_integrated_classification.tsv')
print('行数 %d/%d 列名一致=%s' % (len(a), len(b), fa == fb))
for lv in LV:
    print('  %-9s 差异 %5d 行' % (lv, sum(1 for c in a if a[c].get(lv,'') != b.get(c,{}).get(lv,''))))
PYEOF
