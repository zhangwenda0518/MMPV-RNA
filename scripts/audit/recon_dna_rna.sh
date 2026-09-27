#!/bin/bash
# 只读侦察：第 3 条诉求（按分类补 DNA/RNA 判断）的当前基础
G=~/MMPV-paper/goji-virome/02_novel_virus/RNA-Lycium_barbarum_out
O=~/MMPV-paper/onekp-virome/onekp-virus
P=~/MMPV-RNA/virome_discovery_pipeline

echo "=== [1] 枸杞 05 原表 表头 ==="
head -1 $G/05_Taxonomy/Votus.integrated/final_integrated_classification.tsv | tr '\t' '\n' | nl
echo
echo "=== [2] onekp 05 原表 表头 ==="
head -1 $O/05_Taxonomy/Votus.integrated/final_integrated_classification.tsv | tr '\t' '\n' | nl
echo
echo "=== [3] 枸杞下游植物表 表头 ==="
head -1 $G/10_Reports/All_plant.viruses_info.tsv | tr '\t' '\n' | nl
echo
echo "=== [4] 管线脚本里提到 Riboviria / DNA / RNA 的文件 ==="
ls -1 $P/*.py $P/*.R 2>/dev/null
echo "--- grep Riboviria ---"
grep -ln "Riboviria" $P/*.py $P/*.R 2>/dev/null
echo "--- grep 基因组类型类字段名 ---"
grep -ln -E "genome_type|Genome_Type|nucleic|Nucleic|dna_rna|DNA_or_RNA|MolType|molecule_type" $P/*.py $P/*.R 2>/dev/null
echo
echo "=== [5] 参照库里有没有 DNA/RNA 字段 ==="
ls -la ~/database/taxonomy/ 2>/dev/null | head -20
echo "--- VMR 表头 ---"
head -1 ~/plant_virus_db/1.virus-host_db/B-ictv/VMR_MSL41.tsv | tr '\t' '\n' | nl | head -40
