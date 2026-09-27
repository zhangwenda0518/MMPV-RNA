#!/bin/bash
# smoke.sh — Virseqimprover.py 双模式分发烟雾测试
# 模式1: 默认 bbduk (应看到 runBbdukHarvest)
# 模式2: -readFrac 0.5 (应看到 runAlignment, 覆盖度过滤回退 salmon)
export PATH=/home/zhangwenda/biosoft/ncbi-blast-2.13.0+/bin:/home/zhangwenda/.pixi/bin:/home/zhangwenda/mambaforge/bin:$PATH
B3=/home/zhangwenda/virus/data-2026/data-test/RNA-Alternaria_alternata_out/08_Rescue/Plant/branch_b
R1=/tmp/bench_bbduk_v3/m_R1.fa.gz
R2=/tmp/bench_bbduk_v3/m_R2.fa.gz
cd /tmp/vsi_wire || exit 1
rm -rf out_bbduk out_salmon bbduk.log salmon.log DONE

echo "[smoke] mode1: default (bbduk)"
timeout 1200 python /tmp/vsi_wire/Virseqimprover.py -1 "$R1" -2 "$R2" \
    -scaffold "$B3/scaffolds/contig_116.fasta" -o /tmp/vsi_wire/out_bbduk \
    -t 10 -minSuspiciousLen 50 > /tmp/vsi_wire/bbduk.log 2>&1
echo "bbduk_rc=$?" >> /tmp/vsi_wire/bbduk.log

echo "[smoke] mode2: -readFrac 0.5 (salmon)"
timeout 1200 python /tmp/vsi_wire/Virseqimprover.py -1 "$R1" -2 "$R2" \
    -scaffold "$B3/scaffolds/contig_116.fasta" -o /tmp/vsi_wire/out_salmon \
    -t 10 -minSuspiciousLen 50 -readFrac 0.5 > /tmp/vsi_wire/salmon.log 2>&1
echo "salmon_rc=$?" >> /tmp/vsi_wire/salmon.log

touch /tmp/vsi_wire/DONE
echo "[smoke] all done"
