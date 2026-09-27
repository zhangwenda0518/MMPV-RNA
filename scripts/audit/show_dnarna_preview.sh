#!/bin/bash
# 取回预演样例（属级歧义行、无值行），验证追加列格式
D=/tmp/dnarna_preview
A=$D/goji-virome__02_novel_virus__RNA-Lycium_amarum_out__09_Virome_Analysis__all_plant_analysis__All_plant.viruses_info.tsv
B=$D/goji-virome__02_novel_virus__RNA-Lycium_barbarum_out__09_Virome_Analysis__all_plant_analysis__All_plant.viruses_info.tsv
C=$D/onekp-virome__onekp-virus__09_Virome_Analysis__all_plant_analysis__All_plant.viruses_info.tsv

echo "===== 1) All_plant 表头 + 前 3 行（Lycium_amarum） ====="
head -n 4 "$A" | cut -c1-220

echo
echo "===== 2) 行尾列分布（barbarum 全表） ====="
tail -n +2 "$B" | awk -F'\t' '{print $(NF-2)"\t"$(NF-1)"\t"$NF}' | sort | uniq -c | sort -rn | head -20

echo
echo "===== 3) genus_ambiguous 样例（barbarum） ====="
tail -n +2 "$B" | awk -F'\t' '$NF=="genus_ambiguous"' | cut -f5,9,10,11,12,$((21)) | head -6

echo
echo "===== 4) none 样例（onekp，取前 5 行） ====="
tail -n +2 "$C" | awk -F'\t' '$NF=="none"' | cut -f1,10,11,12 | head -5

echo
echo "===== 5) HQ 表头（barbarum） ====="
head -n 3 "$D/goji-virome__02_novel_virus__RNA-Lycium_barbarum_out__09_Virome_Analysis__HQ_analysis__HQ_plant_viruses_info.tsv" | cut -c1-200

echo
echo "===== 6) genus_summary 表头 + 前 3 行（barbarum） ====="
head -n 4 "$D/goji-virome__02_novel_virus__RNA-Lycium_barbarum_out__09_Virome_Analysis__all_plant_analysis__all_plant_viruses_genus_summary.tsv"

echo
echo "===== 7) DNA 行样例（barbarum，判 DNA 的物种/属） ====="
tail -n +2 "$B" | awk -F'\t' '$(NF-2)=="DNA"' | cut -f10,11,12,$((NF-2)),$((NF-1)),$NF | sort | uniq -c | sort -rn | head -10

echo
echo "===== 8) 缺列检查：原表是否曾被改动（对比 09 层原文件行数） ====="
for f in "$B"; do
  n=$(basename "$f")
  orig="/home/zhangwenda/MMPV-paper/goji-virome/02_novel_virus/RNA-Lycium_barbarum_out/09_Virome_Analysis/all_plant_analysis/All_plant.viruses_info.tsv"
  echo "原: $(wc -l < "$orig") 行 / 预演: $(wc -l < "$f") 行"
done
