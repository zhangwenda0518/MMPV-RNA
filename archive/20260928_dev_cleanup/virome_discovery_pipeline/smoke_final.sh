#!/bin/bash
# smoke_final.sh — 生产原格式 (.fa.gz) 端到端烟雾测试: 默认 bbduk 路径
# 验证: 格式自适应 getReadLen + tips 生成 + bbduk outm 收割 + SPAdes + 多轮迭代
export PATH=/home/zhangwenda/biosoft/ncbi-blast-2.13.0+/bin:/home/zhangwenda/.pixi/bin:/home/zhangwenda/mambaforge/bin:$PATH
B3=/home/zhangwenda/virus/data-2026/data-test/RNA-Alternaria_alternata_out/08_Rescue/Plant/branch_b
R1=/tmp/bench_bbduk_v3/m_R1.fa.gz
R2=/tmp/bench_bbduk_v3/m_R2.fa.gz
cd /tmp/vsi_wire || exit 1
rm -rf out_bbduk bbduk.log DONE

timeout 1500 python -u /tmp/vsi_wire/Virseqimprover.py -1 "$R1" -2 "$R2" \
    -scaffold "$B3/scaffolds/contig_116.fasta" -o /tmp/vsi_wire/out_bbduk \
    -t 10 -minSuspiciousLen 50 > /tmp/vsi_wire/bbduk.log 2>&1
echo "rc=$?" >> /tmp/vsi_wire/bbduk.log
touch /tmp/vsi_wire/DONE
