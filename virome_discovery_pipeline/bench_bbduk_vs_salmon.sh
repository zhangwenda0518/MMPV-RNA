#!/bin/bash
# bench_bbduk_vs_salmon.sh v3 — 全面测试 BBDuk 替代 salmon(+filterbyname)
#
# v2 教训 (保留备查):
#   - salmon 0.8.1 不支持逗号多文件 → 本版先 cat 合并成单文件 (等价生产 merged_reads)
#   - reads 是真 FASTA; BBDuk 的 FASTA 读取器要 mmap 真实文件, fifo 不行
#   - salmon 0.8.1 走 pigz fifo 直接报 invoked improperly → fifo 路线废弃,
#     改测「预解压明文」策略 (生产集成候选)
#
# 步骤:
#   [P] 预备(一次性, 单独计时): cat 合并 gz 对; pigz 解压明文对
#   [1] salmon 0.8.1 直读合并gz 全长ref   — 生产现状基线
#   [2] salmon 0.8.1 读明文   全长ref     — 预解压策略速度
#   [3] filterbyname 直读合并gz (生产第二趟)
#   [4] bbduk k=25 明文 全长ref
#   [4b] bbduk k=25 直读gz  全长ref       — 输入模式对比 (定生产集成方案)
#   [5] bbduk k=20 明文; [6] bbduk k=31 明文 — k 敏感性
#   [7] salmon 直读gz TIPref              — 生产真实场景 (两端 225bp)
#   [8] bbduk k=25 明文 TIPref
#   [9] 漏捞法证: blastn 看 miss reads 与 ref 的比对长度分布
set -uo pipefail
export LC_ALL=C

DATA=/home/zhangwenda/virus/data-2026/1.data
B3=/home/zhangwenda/virus/data-2026/data-test/RNA-Alternaria_alternata_out/08_Rescue/Plant/branch_b
R1A=$DATA/SRR23083385_1.fa.gz; R2A=$DATA/SRR23083385_2.fa.gz
R1B=$DATA/SRR26155732_1.fa.gz; R2B=$DATA/SRR26155732_2.fa.gz
REF_FULL=$B3/scaffolds/contig_116.fasta
REF_TIPS=$B3/out_contig_116/scaffold-truncated/scaffold-start-end.fasta
SALMON=${SALMON:-/home/zhangwenda/mambaforge/envs/Virseqimprover/bin/salmon}
BBDUK=${BBDUK:-/home/zhangwenda/mambaforge/bin/bbduk.sh}
FBN=${FBN:-/home/zhangwenda/mambaforge/bin/filterbyname.sh}
BLASTN=${BLASTN:-/home/zhangwenda/biosoft/ncbi-blast-2.13.0+/bin/blastn}
T=${T:-10}
WORK=${WORK:-/tmp/bench_bbduk_v3}

mkdir -p "$WORK" && cd "$WORK"
log()  { echo "[$(date +%H:%M:%S)] $*"; }
load() { echo -n "    (load: "; cut -d' ' -f1 /proc/loadavg; echo ")"; }
now_ms() { date +%s%3N; }
t0=0; step() { log "$*"; load; t0=$(now_ms); }
el()   { echo "    => $(( $(now_ms) - t0 ))ms"; }

names_fasta() { awk 'NR%2==1 {sub(/^>/,""); print $1}' "$1" | sort -u; }
names_sam()   { awk '!/^@/ && int($2/4)%2==0 {print $1}' "$1" | sort -u; }
cnt_sam()     { awk '!/^@/ && int($2/4)%2==0' "$1" | wc -l; }

echo "================ 环境 ================"
echo "salmon: $($SALMON --version 2>&1 | head -1)"
echo "bbduk:  $BBDUK"
echo "threads=$T  cores=$(nproc)"
echo "REF_FULL: $(grep -c '^>' "$REF_FULL") seq;  REF_TIPS: $(grep -c '^>' "$REF_TIPS") seq"
echo ""

# ── [P] 预备 (一次性成本, 不计入每轮) ───────────────────────────
step "[P1] cat 合并 gz 对 (复制 1.6G)"
cat "$R1A" "$R1B" > m_R1.fa.gz; cat "$R2A" "$R2B" > m_R2.fa.gz
el; ls -lh m_R1.fa.gz m_R2.fa.gz | awk '{print "   ", $9, $5}'
step "[P2] pigz 解压明文 (每任务一次的预解压成本)"
pigz -dc -p "$T" m_R1.fa.gz > p_R1.fa; pigz -dc -p "$T" m_R2.fa.gz > p_R2.fa
el; ls -lh p_R1.fa p_R2.fa | awk '{print "   ", $9, $5}'

# ── [1] salmon 直读 gz (生产现状) ──────────────────────────────
rm -rf idx_full res1 sm_direct.sam
step "[1] salmon 0.8.1 直读gz 全长ref (index)"
$SALMON index -t "$REF_FULL" -i idx_full > idx.log 2>&1 && echo "    index ok" || { echo "    INDEX FAIL"; tail -3 idx.log; }
step "[1] salmon 0.8.1 直读gz 全长ref (quant)"
if ! ($SALMON quant -i idx_full -l A -1 m_R1.fa.gz -2 m_R2.fa.gz -o res1 \
        --writeMappings -p "$T" --quasiCoverage 0 2> sal_direct.log \
        | samtools view -h -F 0x04 - > sm_direct.sam); then
    echo "    SALMON DIRECT FAIL"; tail -5 sal_direct.log
fi
el; echo "    mapped records: $(cnt_sam sm_direct.sam)"

# ── [2] salmon 读明文 ──────────────────────────────────────────
step "[2] salmon 0.8.1 明文 全长ref"
if ! ($SALMON quant -i idx_full -l A -1 p_R1.fa -2 p_R2.fa -o res2 \
        --writeMappings -p "$T" --quasiCoverage 0 2> sal_plain.log \
        | samtools view -h -F 0x04 - > sm_plain.sam); then
    echo "    SALMON PLAIN FAIL"; tail -5 sal_plain.log
fi
el; echo "    mapped records: $(cnt_sam sm_plain.sam)"

# ── [3] filterbyname 复刻 ──────────────────────────────────────
step "[3] filterbyname 直读gz (生产第二趟)"
if ! bash "$FBN" in=m_R1.fa.gz in2=m_R2.fa.gz out=fb_1.fa out2=fb_2.fa \
        names=sm_direct.sam include=t overwrite=t > fb.log 2>&1; then
    echo "    FILTERBYNAME FAIL"; tail -5 fb.log
fi
el; echo "    fb_1 records: $(( $(wc -l < fb_1.fa 2>/dev/null || echo 0) / 2 ))"

# ── [4/4b/5/6] bbduk ───────────────────────────────────────────
step "[4] bbduk k=25 明文 全长ref"
bash "$BBDUK" in=p_R1.fa in2=p_R2.fa outm=bd25_1.fa outm2=bd25_2.fa \
     ref="$REF_FULL" k=25 rcomp=t threads=$T overwrite=t > bd25.log 2>&1 \
     || { echo "    BBDUK FAIL"; tail -5 bd25.log; }
el; echo "    bd25_1 records: $(( $(wc -l < bd25_1.fa 2>/dev/null || echo 0) / 2 ))"

step "[4b] bbduk k=25 直读gz 全长ref (输入模式对比)"
bash "$BBDUK" in=m_R1.fa.gz in2=m_R2.fa.gz outm=bdg_1.fa outm2=bdg_2.fa \
     ref="$REF_FULL" k=25 rcomp=t threads=$T overwrite=t > bdg.log 2>&1 \
     || { echo "    BBDUK FAIL"; tail -5 bdg.log; }
el; echo "    bdg_1 records: $(( $(wc -l < bdg_1.fa 2>/dev/null || echo 0) / 2 ))"

for K in 20 31; do
    step "[k$K] bbduk k=$K 明文 全长ref"
    bash "$BBDUK" in=p_R1.fa in2=p_R2.fa outm=bd${K}_1.fa outm2=bd${K}_2.fa \
         ref="$REF_FULL" k=$K rcomp=t threads=$T overwrite=t > bd$K.log 2>&1 \
         || { echo "    BBDUK k=$K FAIL"; tail -5 bd$K.log; }
    el; echo "    bd${K}_1 records: $(( $(wc -l < bd${K}_1.fa 2>/dev/null || echo 0) / 2 ))"
done

# ── [7][8] TIP ref ─────────────────────────────────────────────
rm -rf idx_tips
step "[7] salmon 直读gz TIPref (index+quant)"
$SALMON index -t "$REF_TIPS" -i idx_tips > idx2.log 2>&1 && echo "    index ok"
if ! ($SALMON quant -i idx_tips -l A -1 m_R1.fa.gz -2 m_R2.fa.gz -o res7 \
        --writeMappings -p "$T" --quasiCoverage 0 2> sal_tips.log \
        | samtools view -h -F 0x04 - > smtips.sam); then
    echo "    SALMON TIPS FAIL"; tail -5 sal_tips.log
fi
el; echo "    mapped records: $(cnt_sam smtips.sam)"

step "[8] bbduk k=25 明文 TIPref"
bash "$BBDUK" in=p_R1.fa in2=p_R2.fa outm=bdt_1.fa outm2=bdt_2.fa \
     ref="$REF_TIPS" k=25 rcomp=t threads=$T overwrite=t > bdt.log 2>&1 \
     || { echo "    BBDUK TIPS FAIL"; tail -5 bdt.log; }
el; echo "    bdt_1 records: $(( $(wc -l < bdt_1.fa 2>/dev/null || echo 0) / 2 ))"

# ── 指标 ────────────────────────────────────────────────────────
echo ""
echo "================ 指标 ================"
names_sam sm_direct.sam > sal.names
names_sam smtips.sam    > sal_tips.names
[ -s fb_1.fa ]     && names_fasta fb_1.fa   > fb.names   || : > fb.names
[ -s bd25_1.fa ]   && names_fasta bd25_1.fa > bd25.names || : > bd25.names
[ -s bd20_1.fa ]   && names_fasta bd20_1.fa > bd20.names || : > bd20.names
[ -s bd31_1.fa ]   && names_fasta bd31_1.fa > bd31.names || : > bd31.names
[ -s bdt_1.fa ]    && names_fasta bdt_1.fa  > bdt.names  || : > bdt.names
[ -s bdg_1.fa ]    && names_fasta bdg_1.fa  > bdg.names  || : > bdg.names

recall() { local miss extra n
    miss=$(comm -13 "$2" "$1" | wc -l); extra=$(comm -23 "$2" "$1" | wc -l); n=$(wc -l < "$1")
    awk -v n=$n -v m=$miss -v e=$extra 'BEGIN{if(n>0) printf "%.3f%% (miss=%d extra=%d truth=%d)", 100*(1-m/n), m, e, n; else print "NA"}'
}
PAIR_OK=$([ "$(names_fasta bd25_1.fa 2>/dev/null | md5sum)" = "$(names_fasta bd25_2.fa 2>/dev/null | md5sum)" ] && echo OK || echo DESYNC)
FB_N=$(wc -l < fb.names)
echo "salmon 直读gz vs 明文 mapped 一致:  $(cnt_sam sm_direct.sam) vs $(cnt_sam sm_plain.sam)"
echo "bbduk 明文 vs 直读gz 名字集合一致:   $(cmp -s bd25.names bdg.names && echo OK || echo DIFF)"
echo "bbduk 配对完整 (R1/R2 名字集合):      $PAIR_OK"
echo "fb (filterbyname 真值) 配对数:        $FB_N"
echo "recall bbduk k=25 vs fb @全长:        $(recall fb.names bd25.names)"
echo "recall bbduk k=20 vs fb @全长:        $(recall fb.names bd20.names)"
echo "recall bbduk k=31 vs fb @全长:        $(recall fb.names bd31.names)"
echo "recall bbduk k=25 vs salmon @TIP:     $(recall sal_tips.names bdt.names)"

# ── [9] 漏捞法证 ────────────────────────────────────────────────
echo ""
echo "================ 漏捞法证 (k=25 全长) ================"
comm -13 bd25.names fb.names > miss25.names
MISSN=$(wc -l < miss25.names)
echo "漏捞总数: $MISSN / $FB_N"
if [ "$MISSN" -gt 0 ] && [ -s p_R1.fa ]; then
    bash "$FBN" in=p_R1.fa in2=p_R2.fa out=miss_1.fa out2=miss_2.fa names=miss25.names \
         include=t overwrite=t > missfb.log 2>&1 || tail -3 missfb.log
    cat miss_1.fa miss_2.fa > miss_all.fa 2>/dev/null
    "$BLASTN" -query miss_all.fa -subject "$REF_FULL" -task blastn -dust no \
        -outfmt '6 qseqid pident length' -max_target_seqs 1 -evalue 1e-5 \
        -num_threads "$T" > miss.tsv 2>/dev/null || echo "    blastn fail"
    awk '{n++; L[n]=$3; if($3<25) s25++; else if($3<50) s50++; else ge50++}
        END{
          if(n==0){print "    漏捞 reads 无 blast 命中 (salmon 宽松判定的噪音)"; exit}
          asort(L); med=L[int((n+1)/2)]
          printf "    有命中 %d/%d: 比对长度中位数=%d\n    <25bp (k=25 物理捞不到): %d\n    25-49bp (边缘重叠): %d\n    >=50bp (需检查): %d\n", n,'"$MISSN"',med,s25+0,s50+0,ge50+0
        }' miss.tsv
fi

echo ""
echo "================ 速度汇总 (各步 => 耗时) ================"
uptime
touch "$WORK/DONE"
log "ALL DONE"
