#!/bin/bash
# 验证补丁后的 minibwa 路径能正常产出
# 直接跑 VSI 的比对命令形态, 不跑完整 SPAdes 流程
set -e

WORK=/tmp/vsi_patch_verify
rm -rf $WORK && mkdir -p $WORK && cd $WORK

CONTIG=/home/zhangwenda/data-test/out/08_Rescue/Plant/branch_b/scaffolds/ERR2040118_clean_NODE_150_length_3726_cov_29.718708.fasta
R1=/home/zhangwenda/data-test/out/00b_HostDepletion/ERR2040118_clean_1.fa.gz
R2=/home/zhangwenda/data-test/out/00b_HostDepletion/ERR2040118_clean_2.fa.gz
THREADS=10
MBWA=/home/zhangwenda/bin/minibwa

cp $CONTIG scaffold.fasta

echo "=== 模拟补丁后的 minibwa 命令 ==="
cat > run.sh <<EOF
#!/bin/bash
set -o pipefail
cd $WORK
$MBWA index -t $THREADS scaffold.fasta > /dev/null 2>&1
$MBWA map -u -t $THREADS scaffold.fasta $R1 $R2 | samtools sort -@ $THREADS - -o bowtie2-mapped.bam
samtools depth -a bowtie2-mapped.bam > samtools-coverage.txt
EOF
cat run.sh
echo ""
bash run.sh
echo "=== 产物检查 ==="
ls -la bowtie2-mapped.bam samtools-coverage.txt
echo "BAM records: $(samtools view -c bowtie2-mapped.bam)"
echo "coverage lines: $(wc -l < samtools-coverage.txt)"
echo "coverage 前3行:"; head -3 samtools-coverage.txt
echo ""
echo "=== 对比 bowtie2 基线 (3193 records / depth 282212) ==="
samtools depth -a bowtie2-mapped.bam | awk '{s+=$3} END{print "minibwa depth sum: "s}'
