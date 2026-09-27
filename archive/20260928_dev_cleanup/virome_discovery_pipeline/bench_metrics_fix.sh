#!/bin/bash
# bench_metrics_fix.sh — 修正 v3 的指标计算 (多行 FASTA 解析 bug), 只重算不重跑工具
# bug: reads FASTA 序列 70 字符折行 (4行/记录), 旧 names_fasta 用 NR%2==1 把序列行当名字
set -uo pipefail
export LC_ALL=C
cd /tmp/bench_bbduk_v3 || exit 1

# 正确解析: header 行以 > 开头
names_fasta() { awk '/^>/{sub(/^>/,""); print $1}' "$1" | sort -u; }
names_sam()   { awk '!/^@/ && int($2/4)%2==0 {print $1}' "$1" | sort -u; }
cnt_sam()     { awk '!/^@/ && int($2/4)%2==0' "$1" | wc -l; }
records()     { grep -c '^>' "$1" 2>/dev/null || echo 0; }

recall() { local miss extra n
    miss=$(comm -13 "$2" "$1" | wc -l); extra=$(comm -23 "$2" "$1" | wc -l); n=$(wc -l < "$1")
    awk -v n=$n -v m=$miss -v e=$extra 'BEGIN{if(n>0) printf "%.3f%% (miss=%d extra=%d truth=%d)", 100*(1-m/n), m, e, n; else print "NA"}'
}

names_sam sm_direct.sam > sal.names
names_sam smtips.sam    > sal_tips.names
names_fasta fb_1.fa     > fb.names
names_fasta bd25_1.fa   > bd25.names
names_fasta bd20_1.fa   > bd20.names
names_fasta bd31_1.fa   > bd31.names
names_fasta bdt_1.fa    > bdt.names
names_fasta bdg_1.fa    > bdg.names
names_fasta fb_2.fa     > fb2.names
names_fasta bd25_2.fa   > bd252.names

echo "================ 修正后指标 ================"
echo "-- 记录数 (对) --"
printf "salmon@全长 mapped records=%s (唯一名 %s)\n" "$(cnt_sam sm_direct.sam)" "$(wc -l < sal.names)"
printf "salmon@TIP  mapped records=%s (唯一名 %s)\n" "$(cnt_sam smtips.sam)" "$(wc -l < sal_tips.names)"
printf "fb=%s bd25=%s bd20=%s bd31=%s bdg=%s bdt=%s\n" \
    "$(records fb_1.fa)" "$(records bd25_1.fa)" "$(records bd20_1.fa)" "$(records bd31_1.fa)" "$(records bdg_1.fa)" "$(records bdt_1.fa)"
echo ""
echo "-- 一致性 --"
echo "bbduk 明文 vs 直读gz 名字集合:  $(cmp -s bd25.names bdg.names && echo 一致 || echo DIFF)"
echo "bbduk 配对完整 (R1 vs R2 名字):  $(cmp -s bd25.names bd252.names && echo OK || echo DESYNC)"
echo "fb 配对完整 (R1 vs R2 名字):     $(cmp -s fb.names fb2.names && echo OK || echo DESYNC)"
echo "salmon 直读gz vs 明文 mapped:    $(cnt_sam sm_direct.sam) vs $(cnt_sam sm_plain.sam)"
echo ""
echo "-- 召回 (真值=filterbyname 从 salmon SAM 提取) --"
FB_N=$(wc -l < fb.names)
echo "recall bbduk k=25 vs fb @全长:   $(recall fb.names bd25.names)"
echo "recall bbduk k=20 vs fb @全长:   $(recall fb.names bd20.names)"
echo "recall bbduk k=31 vs fb @全长:   $(recall fb.names bd31.names)"
echo "recall bbduk k=25 vs salmon @TIP: $(recall sal_tips.names bdt.names)"
echo "fb vs salmon@全长 唯一名:         $(recall sal.names fb.names)"

echo ""
echo "================ 漏捞法证重做 (k=25 全长) ================"
comm -13 bd25.names fb.names > miss25.names
MISSN=$(wc -l < miss25.names)
echo "漏捞总数: $MISSN / $FB_N"
if [ "$MISSN" -gt 0 ]; then
    bash /home/zhangwenda/mambaforge/bin/filterbyname.sh in=p_R1.fa in2=p_R2.fa \
        out=miss_1.fa out2=miss_2.fa names=miss25.names include=t overwrite=t > missfb.log 2>&1 || tail -3 missfb.log
    cat miss_1.fa miss_2.fa > miss_all.fa 2>/dev/null
    echo "取回 records: $(records miss_all.fa)"
    /home/zhangwenda/biosoft/ncbi-blast-2.13.0+/bin/blastn -query miss_all.fa \
        -subject /home/zhangwenda/virus/data-2026/data-test/RNA-Alternaria_alternata_out/08_Rescue/Plant/branch_b/scaffolds/contig_116.fasta \
        -task blastn -dust no -outfmt '6 qseqid pident length' -max_target_seqs 1 -evalue 1e-5 \
        -num_threads 10 > miss.tsv 2>/dev/null || echo "blastn fail"
    awk '{n++; if($3<25) s25++; else if($3<50) s50++; else ge50++; s+=$2}
        END{
          if(n==0){print "漏捞 reads 无 blast 命中"; exit}
          print "有命中的漏捞: " n
          print "平均identity: " s/n "%"
          print "比对长度<25bp (k=25物理捞不到): " s25+0
          print "25-49bp (边缘重叠read): " s50+0
          print ">=50bp: " ge50+0
        }' miss.tsv
fi
touch DONE_FIX
echo "ALL DONE"
