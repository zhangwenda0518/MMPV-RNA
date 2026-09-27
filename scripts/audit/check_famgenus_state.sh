#!/bin/bash
# 只读核验：科-属不一致工作流的当前状态（机制层 vs 产物层）
G=~/MMPV-paper/goji-virome/02_novel_virus/RNA-Lycium_barbarum_out
V=$G/05_Taxonomy/Votus.integrated
R=$G/10_Reports

echo "=== [1] R 脚本（科属闸门）==="
cd ~/MMPV-RNA/virome_discovery_pipeline || exit 1
md5sum virus_classifier_analysis.R
stat -c '%y %s B' virus_classifier_analysis.R
echo "备份链条数: $(ls -1 virus_classifier_analysis.R.bak_* 2>/dev/null | wc -l)"

echo
echo "=== [2] 枸杞 05 原表 vs 校准表 ==="
wc -l $V/final_integrated_classification.tsv $V/calibration_20260914/final_integrated_classification.calibrated.tsv
stat -c '%y %n' $V/final_integrated_classification.tsv $V/calibration_20260914/final_integrated_classification.calibrated.tsv
ls -la $V/calibration_20260914/

echo
echo "=== [3] 下游是否引用 calibrated（空=没接）==="
grep -rl calibrated $R 2>/dev/null | head -5
echo "命中数: $(grep -rl calibrated $R 2>/dev/null | wc -l)"

echo
echo "=== [4] 枸杞下游植物表 ==="
wc -l $R/All_plant.viruses_info.tsv $R/All_plant.viruses_info.tsv.bak_kill_20260915
stat -c '%y %n' $R/All_plant.viruses_info.tsv $R/All_plant.viruses_info.tsv.bak_kill_20260915
echo "全部 bak_kill_20260915 文件:"
ls -1 $R/*.bak_kill_20260915 2>/dev/null

echo
echo "=== [5] onekp 下游植物表 ==="
O=~/MMPV-paper/onekp-virome/onekp-virus
ls -la $O/10_Reports/All_plant.viruses_info.tsv 2>/dev/null
wc -l $O/10_Reports/All_plant.viruses_info.tsv 2>/dev/null

echo
echo "=== [6] 参照表与复核脚本 ==="
ls -la ~/database/taxonomy/genus_family_ref.tsv 2>/dev/null
ls -la /tmp/chk_famgenus_state.py 2>/dev/null || echo "chk_famgenus_state.py 已不在 /tmp"
