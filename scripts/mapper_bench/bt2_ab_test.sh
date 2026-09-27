#!/bin/bash
# A/B 测试: 旧管线(view 中间层) vs 新管线(--no-unal) 在同一 reads/contig 上的耗时
set -e
cd /tmp/bt2_ab_test
R1=/home/zhangwenda/data-test/out/00b_HostDepletion/ERR2041160_clean_1.fa.gz
R2=/home/zhangwenda/data-test/out/00b_HostDepletion/ERR2041160_clean_2.fa.gz
FA=ERR2041160_clean_NODE_362_length_3154_cov_1116.864674.fasta

# 索引只需建一次
[ -f idx.1.bt2 ] || bowtie2-build --threads 10 $FA idx > build.log 2>&1

echo "=== A: 旧管线 (view x2) ==="
/usr/bin/time -f "wall=%es cpu=%Ps" bash -c "
bowtie2 --local --threads 10 -f -x idx -1 $R1 -2 $R2 2>/dev/null | \
samtools view -bS - | samtools view -h -F 0x04 -b - | \
samtools sort -@ 10 - -o A.bam 2>/dev/null" || true

echo "=== B: 新管线 (--no-unal) ==="
/usr/bin/time -f "wall=%es cpu=%Ps" bash -c "
bowtie2 --local --no-unal --threads 10 -f -x idx -1 $R1 -2 $R2 2>/dev/null | \
samtools sort -@ 10 - -o B.bam 2>/dev/null" || true

echo "=== 产物对比 ==="
ls -la A.bam B.bam 2>/dev/null | awk '{print $5, $9}'
samtools view -c A.bam 2>/dev/null | xargs echo "A records:"
samtools view -c B.bam 2>/dev/null | xargs echo "B records:"
samtools depth -a A.bam 2>/dev/null | md5sum | cut -c1-12 | xargs echo "A depth md5:"
samtools depth -a B.bam 2>/dev/null | md5sum | cut -c1-12 | xargs echo "B depth md5:"
