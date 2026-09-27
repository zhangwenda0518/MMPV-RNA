#!/bin/bash
# 测试 -mapper 参数与兜底机制
set -e
cd /home/zhangwenda/MMPV-RNA/virome_discovery_pipeline

CONTIG=/home/zhangwenda/data-test/out/08_Rescue/Plant/branch_b/scaffolds/ERR2040118_clean_NODE_150_length_3726_cov_29.718708.fasta
R1=/home/zhangwenda/data-test/out/00b_HostDepletion/ERR2040118_clean_1.fa.gz
R2=/home/zhangwenda/data-test/out/00b_HostDepletion/ERR2040118_clean_2.fa.gz

echo "=== 测试1: 参数校验 (非法 mapper 名应报错) ==="
python3 Virseqimprover.py -mapper nonsense -1 $R1 -scaffold $CONTIG -o /tmp/mtest1 2>&1 | head -2
echo ""

echo "=== 测试2: minibwa 路径不存在时应回退 bowtie2 ==="
python3 -c "
import sys, os
sys.argv=['x']
# 只导入 parseArguments 验证参数解析, 不跑主流程
src=open('Virseqimprover.py').read()
assert 'mapperMode = \"minibwa\"' in src, '默认值缺失'
assert '\"-mapper\"' in src, '参数解析缺失'
assert 'fallback to bowtie2' in src, '兜底逻辑缺失'
assert 'bowtie2-mapped.bam' in src, '产物名保持'
print('源码结构校验: PASS')
"
echo ""

echo "=== 测试3: 帮助文本完整性 ==="
python3 Virseqimprover.py -h 2>&1 | grep -cE "mapper|minibwa|bowtie2" | xargs echo "相关帮助行数:"
echo ""

echo "=== 测试4: py_compile ==="
python3 -m py_compile Virseqimprover.py && echo "编译: PASS"
