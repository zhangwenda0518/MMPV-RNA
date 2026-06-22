#!/bin/bash
# 批量运行 virome_pipeline.py (跳过 clean+deplete, 从 assembly 开始)
set -e
cd ~/virus/data-2026/data-test

# 物种: 名称 宿主
THREADS=16
JOBS=4
TAX_JOBS=1

run_one() {
    local name="$1"
    local host="$2"
    local out_dir="${name}_out"
    echo "================================================"
    echo "Processing: $name (host=$host)"
    echo "Output: $out_dir"
    echo "================================================"

    python3 ~/MMPV-RNA/virome_discovery_pipeline/virome_pipeline.py \
      --stage assembly identification cobra cluster taxonomy host checkv rescue report \
      --output_dir "$out_dir" \
      --host-filter "$host" \
      --skip_clean --skip_depletion \
      --coassembly \
      --threads "$THREADS" --jobs "$JOBS" --tax_jobs "$TAX_JOBS" \
      --checkv_db /home/zhangwenda/database/virus-db/checkv-db-v1.7 \
      --genomad_db /home/zhangwenda/database/virus-db/genomad_db/ \
      --mmseqs_db /home/zhangwenda/database/virus-db/RVDB-30/RVDB.mmseqs \
      --virus_db /home/zhangwenda/database/virus-db/ \
      --host_db /home/zhangwenda/database/host_db/ \
      --blast_db /home/zhangwenda/database/virus-db/ncbi-virus_ref/ncbi-virus_ref.blast.db \
      --nr_db /home/zhangwenda/database/nr_diamond/nr.dmnd \
      --db-dir /home/zhangwenda/database/virus-db/ \
      2>&1 | tee "${name}.log"

    echo "=== $name DONE ==="
}

# 全部按 Plant 宿主处理
for sp in RNA-Lycium_barbarum RNA-Lycium_chinense RNA-Lycium_ruthenicum \
          RNA-Alternaria_alternata RNA-Fusarium_nematophilum \
          RNA-Aphis_gossypii RNA-Neoceratitis_asiatica METAGENOMIC; do
    run_one "$sp" "Plant"
done

echo "ALL DONE"
