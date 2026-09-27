#!/bin/bash
# 批量 SDT 属级序列一致性矩阵 (接近属平均长度的 KEEP)
cd /home/zhangwenda/MMPV-RNA/virome_analysis_pipeline
OUT=/home/zhangwenda/virus/data-2026/SDT_matrix
mkdir -p $OUT
declare -A DS=(
  [barbarum]=/home/zhangwenda/virus/data-2026/data-test/RNA-Lycium_barbarum_out/09b_ACVirus_Analysis
  [ruthenicum]=/home/zhangwenda/virus/data-2026/data-test/RNA-Lycium_ruthenicum_out/09b_ACVirus_Analysis
  [chinense]=/home/zhangwenda/virus/data-2026/data-test/RNA-Lycium_chinense_out/09b_ACVirus_Analysis
  [amarum]=/home/zhangwenda/virus/data-2026/data-test/RNA-Lycium_amarum_out/09b_ACVirus_Analysis
  [Fusarium]=/home/zhangwenda/virus/data-2026/data-test/RNA-Fusarium_nematophilum_out/09b_ACVirus_Analysis
  [Alternaria]=/home/zhangwenda/virus/data-2026/data-test/RNA-Alternaria_alternata_out/09b_ACVirus_Analysis
  [Aphis]=/home/zhangwenda/virus/data-2026/data-test/RNA-Aphis_gossypii_out/09b_ACVirus_Analysis
  [mix]=/home/zhangwenda/virus/ningxiagouqi/11.merge_assembly/mix/out/09b_ACVirus_Analysis
)
echo "[$(date '+%H:%M:%S')] SDT 批量开始" > /tmp/sdt_batch.log
seen=""
while read genus ds cid len ratio; do
  [ "$genus" = "genus" ] && continue
  key="$ds|$genus"
  case " $seen " in *" $key "*) continue;; esac
  seen="$seen $key "
  d=${DS[$ds]}
  [ -z "$d" ] && continue
  [ -d "$d" ] || continue
  od="$OUT/${ds}__${genus}"
  if [ -f "$od/SDT_${genus}.png" ] || [ -f "$od/SDT_${genus}.pdf" ]; then
    echo "[skip] $ds $genus (已存在)" >> /tmp/sdt_batch.log; continue
  fi
  echo "[$(date '+%H:%M:%S')] === $ds $genus ===" >> /tmp/sdt_batch.log
  python3 sdt_genus_matrix.py --analysis_dir "$d" --genus "$genus" -o "$od" \
    --threads 16 --short_labels --palette cividis --resume --other_count 0 \
    >> /tmp/sdt_batch.log 2>&1
  echo "[$(date '+%H:%M:%S')] $ds $genus rc=$?" >> /tmp/sdt_batch.log
done < /tmp/sdt_candidate_genera.tsv
echo "[$(date '+%H:%M:%S')] SDT 全部完成" >> /tmp/sdt_batch.log
echo ALL_DONE > /tmp/sdt_batch.done