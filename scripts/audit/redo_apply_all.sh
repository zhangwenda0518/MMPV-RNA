#!/bin/bash
# 复原 Alternaria All_plant 后重跑 apply，补齐第 24 个文件的记录
set -u
A=/home/zhangwenda/MMPV-paper/goji-virome/02_novel_virus/RNA-Alternaria_alternata_out/09_Virome_Analysis/all_plant_analysis/All_plant.viruses_info.tsv
cp -p "$A.bak_dnarna_20260830" "$A" && echo "已复原: $A"
head -n 1 "$A" | tr '\t' '\n' | wc -l | xargs echo "复原后列数(应为20):"
rm -f /tmp/dnarna_apply_20260830.log
python3 /tmp/dnarna_apply.py --apply > /tmp/dnarna_apply_run.log 2>&1
echo "EXIT=$?"
echo "PASS 条数: $(grep -c 'PASS' /tmp/dnarna_apply_run.log)"
tail -n 3 /tmp/dnarna_apply_run.log
