#!/bin/bash
# v6.4 验收：日志关键行 + 台账 + 四项探针
set -u
cd /tmp
OLD=/home/zhangwenda/MMPV-paper/goji-virome/02_novel_virus/RNA-Lycium_barbarum_out/05_Taxonomy/Votus.integrated/final_integrated_classification.tsv
L63=/tmp/rc_legacy/final_integrated_classification.tsv
C63=/tmp/rc_cascade/final_integrated_classification.tsv
L64=/tmp/rc64_legacy/final_integrated_classification.tsv
C64=/tmp/rc64_cascade/final_integrated_classification.tsv

echo "########## 1. v6.4 legacy 日志关键行"
grep -E "科-种相容性|科-属校准|逐级相容性|占位值清理|后缀|自检|工具权重|完成|ERROR|WARN|阶段|合计" /tmp/rc64_legacy.log | head -40
echo "########## 2. v6.4 cascade 日志关键行"
grep -E "科-种相容性|科-属校准|逐级相容性|逐级淘汰|占位值清理|后缀|自检|完成|ERROR|WARN|阶段|合计" /tmp/rc64_cascade.log | head -40
echo "########## 3. 台账 legacy"
cat /tmp/rc64_legacy/taxonomy_gate_check.tsv
echo "########## 4. 台账 cascade"
cat /tmp/rc64_cascade/taxonomy_gate_check.tsv
echo "########## 5. 工具权重"
cat /tmp/rc64_legacy/tool_weights.tsv
echo "########## 6. 被置空的科种行 样例"
head -6 /tmp/rc64_legacy/species_containment_blanked.tsv
wc -l /tmp/rc64_legacy/species_containment_blanked.tsv /tmp/rc64_cascade/species_containment_blanked.tsv
echo "########## 7. 四方对账 (old -> v63legacy -> v64legacy -> v64cascade)"
python3 -u cmp_products.py $OLD $L63 $L64 $C64
echo "########## 8. 口径陈旧探针"
python3 -u probe_stale.py $L64 $C64
echo "########## 9. 占位/后缀检查"
python3 -u chk_ph.py $L64 $C64
echo "########## 10. 种级闸门影响面"
python3 -u probe_sp_gate.py $L64 $C64
echo "########## DONE"
