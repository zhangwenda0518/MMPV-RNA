#!/bin/bash
# 找输入 fasta 与 SAM，对比 read 名格式
base=~/goji-virome/03_known_virus/known_virus_all_v2/09_dvg/virema_results
d=$base/Cytorhabdovirus_sp._lycii/CRR1126132_OR489165.1
echo "== 对照组目录 =="; ls "$d"
echo; echo "== SAM 头 5 行 =="
head -5 "$d/CRR1126132_OR489165.1_ViReMa.sam" 2>/dev/null | cut -c1-120
echo; echo "== onekp 目录 =="
d2=~/data-test/onekp_analysis/08_dvg/virema_results/Verbena_latent_virus/ERR2040134_PV805050.1
ls -R "$d2" 2>/dev/null | head -20
echo; echo "== onekp sandbox 内容 =="
ls "$d2/virema_sandbox" 2>/dev/null | head -10
