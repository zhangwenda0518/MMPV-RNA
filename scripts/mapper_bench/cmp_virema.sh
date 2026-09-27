#!/bin/bash
# 对照组成功样本 vs onekp 失败样本：比对率与输入 read 名格式
base=~/goji-virome/03_known_virus/known_virus_all_v2/09_dvg/virema_results
f=$base/Cytorhabdovirus_sp._lycii/CRR1126132_OR489165.1/virema.log
echo "== 对照组成功样本 CRR1126132 =="
grep -E 'reads processed|at least one alignment|failed to align|Reported' "$f" | head -4
echo
echo "== onekp 失败样本 ERR2040134 =="
f2=~/data-test/onekp_analysis/08_dvg/virema_results/Verbena_latent_virus/ERR2040134_PV805050.1/virema.log
grep -E 'reads processed|at least one alignment|failed to align|Reported' "$f2" | head -4
echo
echo "== 对照组输入 fasta 头 2 行 =="
d=$base/Cytorhabdovirus_sp._lycii/CRR1126132_OR489165.1
ls "$d" | head
for fx in "$d"/virema_sandbox/input_virema.fa "$d"/*input*.fa "$d"/*.fa; do
  if [ -f "$fx" ]; then echo "--- $fx"; head -2 "$fx"; break; fi
done
echo
echo "== onekp 输入 fasta 头 2 行 =="
d2=~/data-test/onekp_analysis/08_dvg/virema_results/Verbena_latent_virus/ERR2040134_PV805050.1
ls "$d2" | head
for fx in "$d2"/virema_sandbox/input_virema.fa "$d2"/*input*.fa "$d2"/*.fa; do
  if [ -f "$fx" ]; then echo "--- $fx"; head -2 "$fx"; break; fi
done
