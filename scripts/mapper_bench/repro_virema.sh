#!/bin/bash
# 用复现输入真跑 ViReMa，验证是否复现 KeyError
cd /tmp
rm -rf virema_test && mkdir virema_test && cd virema_test
cp /tmp/repro_input.fa input.fa
head -600 /home/zhangwenda/data-test/out/00b_HostDepletion/ERR2040134_clean_1.fa.gz 2>/dev/null > /dev/null
# 扩展到完整输入（直接用双端合并的完整版）
python3 - << 'PYEOF'
import gzip
with gzip.open('/home/zhangwenda/data-test/out/00b_HostDepletion/ERR2040134_clean_1.fa.gz','rt') as f1, \
     gzip.open('/home/zhangwenda/data-test/out/00b_HostDepletion/ERR2040134_clean_2.fa.gz','rt') as f2, \
     open('input.fa','w') as out:
    for src, tag in [(f1,'/1'),(f2,'/2')]:
        for line in src:
            line=line.rstrip('\n\r')
            if line.startswith('>'):
                if not line.rstrip().endswith(tag):
                    line=line.rstrip()+tag
                out.write(line+'\n')
            else:
                out.write(line+'\n')
print('input built')
PYEOF
ls -la input.fa
/home/zhangwenda/mambaforge/bin/python /home/zhangwenda/MMPV-RNA/biosoft/virema/ViReMa.py \
  /home/zhangwenda/data-test/onekp_analysis/08_dvg/individual_refs/PV805050.1.fasta \
  input.fa ERR_test.sam \
  --Seed 25 --MicroInDel_Length 15 --Defuzz 0 --p 10 \
  -BED -Overwrite --Output_Tag ERR_test --Output_Dir /tmp/virema_test -Fasta 2>&1 | tail -25
