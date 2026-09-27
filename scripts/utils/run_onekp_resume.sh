#!/bin/bash
cd /home/zhangwenda/data-test/onekp_analysis
nohup python3 /home/zhangwenda/MMPV-RNA/virome_analysis_pipeline/auto_known_virus.py \
  --reads_dir /home/zhangwenda/data-test/out/00b_HostDepletion \
  --output_dir /home/zhangwenda/data-test/onekp_analysis \
  --ref_info /home/zhangwenda/plant_virus_db/3.final-ref-virus.db/final.complete_ref_info.tsv \
  --reference /home/zhangwenda/plant_virus_db/3.final-ref-virus.db/final.complete_ref.fasta \
  --stage all --threads 10 --jobs 10 --align_threads 8 --batch_size 10 >> /home/zhangwenda/data-test/onekp_analysis/resume_all_20260910.log 2>&1
