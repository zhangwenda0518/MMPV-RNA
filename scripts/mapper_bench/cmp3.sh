#!/bin/bash
# 看 batch_virema_dvg.py 如何生成输入 fasta，以及 onekp 失败样本的输入文件是否还在
echo "== batch_virema_dvg.py 中 input 生成的代码 =="
grep -n 'input_virema\|fasta\|seqkit\|extract' ~/MMPV-RNA/virome_analysis_pipeline/batch_virema_dvg.py | head -20
echo
echo "== onekp 失败样本 sandbox 是否被清理 =="
ls ~/data-test/onekp_analysis/08_dvg/virema_results/Verbena_latent_virus/ERR2040134_PV805050.1/ -la
echo
echo "== 对照组成功样本有没有 sandbox 残留 =="
find ~/goji-virome/03_known_virus/known_virus_all_v2/09_dvg/virema_results/Cytorhabdovirus_sp._lycii -maxdepth 2 -name '*sandbox*' | head -3
