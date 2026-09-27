#!/bin/bash
# v6.7 三次实跑：count（唯一票首且 >= 门槛）/ weighted 回归 / pregate（闸门关闭）
set -u
R=/tmp/virus_classifier_analysis.R.cascade_v67
echo "=== parse check ==="
Rscript -e "invisible(parse('$R')); cat('PARSE OK\n')" 2>&1 | tail -5
echo "=== patch md5 / lines ==="
md5sum $R; wc -l $R
INT=/home/zhangwenda/MMPV-paper/goji-virome/02_novel_virus/RNA-Lycium_barbarum_out/05_Taxonomy
COMB=$INT/Votus.classed/Votus_combined_taxonomy.tsv
ls -la $COMB
rm -rf /tmp/rc67_count /tmp/rc67_weighted /tmp/pregate67_count
mkdir -p /tmp/rc67_count /tmp/rc67_weighted /tmp/pregate67_count
cd /tmp

echo "=== run 1/3: v6.7 count（默认，正式版） ==="
start=$(date +%s)
MMPV_CONSENSUS_MODE=cascade Rscript $R --combined $COMB --output /tmp/rc67_count > /tmp/rc67_count.log 2>&1
echo "exit=$? elapsed=$(( $(date +%s) - start ))s"

echo "=== run 2/3: v6.7 weighted（回退口径回归对账） ==="
start=$(date +%s)
MMPV_CONSENSUS_MODE=cascade MMPV_VOTE_WEIGHT=weighted Rscript $R --combined $COMB --output /tmp/rc67_weighted > /tmp/rc67_weighted.log 2>&1
echo "exit=$? elapsed=$(( $(date +%s) - start ))s"

echo "=== run 3/3: v6.7 count（闸门关闭，给诊断用） ==="
python3 /tmp/make_pregate_variant.py $R /tmp/virus_classifier_analysis.R.pregate67
MMPV_CONSENSUS_MODE=cascade Rscript /tmp/virus_classifier_analysis.R.pregate67 --combined $COMB --output /tmp/pregate67_count > /tmp/pregate67_count.log 2>&1
echo "exit=$?"

echo "=== 关键日志行（count / weighted） ==="
grep -aE "逐级淘汰|校准|置空|相容性约束|自检|计权模式" /tmp/rc67_count.log
echo "--- weighted ---"
grep -aE "逐级淘汰|校准|置空|自检|计权模式" /tmp/rc67_weighted.log

echo "=== weighted 回归：v67(mode=weighted) vs v64 正式产物 ==="
md5sum /tmp/rc67_weighted/final_integrated_classification.tsv /tmp/rc64b_cascade/final_integrated_classification.tsv
md5sum /tmp/rc67_weighted/tool_weights.tsv /tmp/rc64b_cascade/tool_weights.tsv

echo "=== v66（share>0.5）vs v67（唯一票首且 >=0.5）：成品逐格差异 ==="
python3 - <<'PYEOF'
import csv
def load(p):
    with open(p, newline='', encoding='utf-8') as f:
        rd = csv.DictReader(f, delimiter='\t'); return {r['contig_id']: r for r in rd}, rd.fieldnames
LV = ['Realm','Kingdom','Phylum','Class','Order','Family','Genus','Species']
a, fa = load('/tmp/rc66_count/final_integrated_classification.tsv')
b, fb = load('/tmp/rc67_count/final_integrated_classification.tsv')
print('行数 %d/%d 列名一致=%s' % (len(a), len(b), fa == fb))
tot = set()
for lv in LV:
    d = [c for c in a if a[c].get(lv,'') != b.get(c,{}).get(lv,'')]
    tot.update(d)
    print('  %-9s 差异 %5d 行' % (lv, len(d)))
print('  涉及 contig 合计 %d' % len(tot))
print('  样例（前 12 个发生变化的 contig，Realm..Species 对照）')
n = 0
for c in sorted(tot):
    if n >= 12: break
    n += 1
    print('   ', c)
    print('      v66: ' + ' | '.join(a[c].get(lv,'') or '-' for lv in LV))
    print('      v67: ' + ' | '.join(b[c].get(lv,'') or '-' for lv in LV))
PYEOF
