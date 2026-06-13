#!/bin/bash
# ============================================================================
# MMPV-RNA v2.3 — 服务器一键运行脚本
# 用法: bash run_pipeline.sh [stage]
# ============================================================================

set -e

# ── 配置 (修改此处) ──────────────────────────────────────────────
DATA=/data/project
DB=/data/db
OUT=/data/out
THREADS=120
JOBS=20
MEMORY=256

# ── 主控 ──────────────────────────────────────────────────────────
PIPE="python ~/MMPV-RNA/scripts/virome_pipeline.py"

# ── 运行 ──────────────────────────────────────────────────────────
STAGE=${1:-all}

case $STAGE in
    all)
        $PIPE --stage all \
            --input_reads $DATA/raw/ \
            --output_dir $OUT/ \
            --kraken2_db    $DB/kraken2_db/ \
            --host_align_db $DB/host_align_db/ \
            --virus_db      $DB/virus_db/ \
            --checkv_db     $DB/checkv_db/ \
            --phabox-db     $DB/phabox_db_v2_2/ \
            --prob-dir      $DB/cross_analysis/ \
            --ref-genomes   $DB/ICTV_plant_viruses.fasta \
            --assembler all --aligner bowtie2 --seq_type rna-short \
            --rrna --host-filter Plant \
            --min-length 500 --ani 0.95 --qcov 0.85 \
            -t $THREADS -j $JOBS -m $MEMORY
        ;;

    clean)
        $PIPE --stage clean \
            --input_reads $DATA/raw/ --output_dir $OUT/ \
            -t $THREADS -j $JOBS
        ;;

    deplete)
        $PIPE --stage deplete \
            --output_dir $OUT/ \
            --kraken2_db $DB/kraken2_db/ --host_align_db $DB/host_align_db/ \
            --aligner bowtie2 --seq_type rna-short --rrna \
            -t $THREADS -j $JOBS
        ;;

    assembly)
        $PIPE --stage assembly \
            --output_dir $OUT/ --input_reads $OUT/00b_HostDepletion/ \
            --assembler all --virus_db $DB/virus_db/ \
            -t $THREADS -j $JOBS -m $MEMORY
        ;;

    identification)
        $PIPE --stage identification \
            --output_dir $OUT/ --input_reads $OUT/00b_HostDepletion/ \
            --virus_db $DB/virus_db/ -t $THREADS -j $JOBS
        ;;

    cobra)
        $PIPE --stage cobra \
            --output_dir $OUT/ --input_reads $OUT/00b_HostDepletion/ \
            --assembler all -t $THREADS -j 10
        ;;

    cluster)
        $PIPE --stage cluster \
            --output_dir $OUT/ --input_reads $OUT/00b_HostDepletion/ \
            --ref-genomes $DB/ICTV_plant_viruses.fasta \
            --min-length 500 --ani 0.95 --qcov 0.85 -t $THREADS
        ;;

    taxonomy)
        $PIPE --stage taxonomy \
            --output_dir $OUT/ --virus_db $DB/virus_db/ -t $THREADS
        ;;

    host)
        $PIPE --stage host \
            --output_dir $OUT/ --virus_db $DB/virus_db/ \
            --phabox-db $DB/phabox_db_v2_2/ --prob-dir $DB/cross_analysis/ \
            -t $THREADS
        ;;

    checkv)
        $PIPE --stage checkv \
            --output_dir $OUT/ --checkv_db $DB/checkv_db/ -t $THREADS
        ;;

    rescue)
        $PIPE --stage rescue \
            --output_dir $OUT/ --input_reads $OUT/00b_HostDepletion/ \
            --checkv_db $DB/checkv_db/ --virus_db $DB/virus_db/ \
            --host-filter Plant -t $THREADS -j $JOBS
        ;;

    known-detect)
        python ~/MMPV-RNA/scripts/auto_known_virus.py --stage detect \
            --reads_dir $OUT/00b_HostDepletion/ --output_dir $OUT/known_viruses/ \
            --ref_info $DB/ref_info.tsv --reference $DB/ref.fasta \
            --tool salmon --threads 40 --jobs 4
        ;;

    known-variants)
        python ~/MMPV-RNA/scripts/auto_known_virus.py --stage variants \
            --reads_dir $OUT/00b_HostDepletion/ --output_dir $OUT/known_viruses/ \
            --ref_info $DB/ref_info.tsv --reference $DB/ref.fasta \
            --variant_caller freebayes --snpeff --snpgenie --threads 40 --jobs 4
        ;;

    known-full)
        python ~/MMPV-RNA/scripts/auto_known_virus.py --stage full \
            --reads_dir $OUT/00b_HostDepletion/ --output_dir $OUT/known_viruses/ \
            -j 4 -t 40
        ;;

    validate)
        python ~/MMPV-RNA/scripts/validate_novel_viruses.py \
            -i $OUT/08_Rescue/Plant/centroids/final_centroids.fasta \
            --taxonomy $OUT/05_Taxonomy/integrated/final_integrated_classification.tsv \
            --cdhit-known $OUT/04_CLUSTER/centroids/known_association.tsv \
            --clusters-tsv $OUT/04_CLUSTER/3_vclust/vclust_clusters.tsv \
            --host $OUT/06_HostPrediction/ensemble_host_summary.tsv \
            -o $OUT/09_Validation/
        ;;

    *)
        echo "用法: bash run_pipeline.sh [stage]"
        echo ""
        echo "   all              全流程"
        echo "   clean            清洗"
        echo "   deplete          去宿主"
        echo "   assembly         组装"
        echo "   identification   鉴定"
        echo "   cobra            COBRA延伸"
        echo "   cluster          聚类"
        echo "   taxonomy         分类"
        echo "   host             宿主预测"
        echo "   checkv           CheckV预评估"
        echo "   rescue           拯救"
        echo "   known-detect     已知病毒检测"
        echo "   known-variants   已知病毒变异"
        echo "   known-full       已知病毒全长"
        echo "   validate         交叉验证"
        ;;
esac

echo ""
echo "═══════════════════════════════════════"
echo "  MMPV-RNA $STAGE 完成"
echo "  日志: $OUT/orchestrator.log"
echo "═══════════════════════════════════════"
