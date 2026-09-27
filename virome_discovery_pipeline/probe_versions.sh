#!/bin/bash
P(){ printf "%-16s " "$1"; shift; "$@" 2>&1 | head -1; }
P fastp /home/zhangwenda/bin/fastp --version
P samtools samtools --version
P blastn /home/zhangwenda/biosoft/ncbi-blast-2.13.0+/bin/blastn -version
P coverm /home/zhangwenda/.pixi/bin/coverm --version
echo "--- virsorter env probe ---"
ls /home/zhangwenda/mambaforge/envs/*/bin/virsorter 2>/dev/null
for f in /home/zhangwenda/mambaforge/envs/*/bin/virsorter; do [ -f "$f" ] && head -1 "$f" && dirname "$f"; done 2>/dev/null
find /home/zhangwenda/mambaforge/envs -maxdepth 5 -name "setup.py" -path "*virsorter*" 2>/dev/null | head -2
grep -r "^__version__" /home/zhangwenda/biosoft/virus/viralVerify/*.py 2>/dev/null
echo "--- salmon ---"
ls /home/zhangwenda/biosoft/RNA-seq/salmon-linux-x86_64/bin 2>/dev/null
/home/zhangwenda/biosoft/RNA-seq/salmon-linux-x86_64/bin/salmon --version 2>&1 | head -1
echo "--- snpeff ---"
ls /home/zhangwenda/biosoft/snpEff*/snpEff.jar 2>/dev/null
java -jar /home/zhangwenda/biosoft/snpEff-*/snpEff.jar -version 2>/dev/null | head -1
java -jar /home/zhangwenda/biosoft/snpEff/snpEff.jar -version 2>/dev/null | head -1
echo "--- clumpify/bbmap ---"
ls /opt/sysoft/bbmap* /home/zhangwenda/biosoft/genome-evaluation/bbmap 2>/dev/null | head -5
grep -oE 'VERSION=[^ ]+' /home/zhangwenda/biosoft/genome-evaluation/bbmap/clumpify.sh 2>/dev/null | head -1
/home/zhangwenda/biosoft/genome-evaluation/bbmap/bbmap.sh --version 2>&1 | tail -1
echo "--- cobra ---"
head -20 /home/zhangwenda/biosoft/virus/cobra/cobra.py 2>/dev/null | grep -iE "version|cobra" | head -3
echo "--- ivar ---"
ls /home/zhangwenda/biosoft/virus/ivar 2>/dev/null | head -5
echo "--- suvtk ---"
which suvtk table2asn tbl2asn 2>/dev/null
find /home/zhangwenda -maxdepth 4 -iname "*suvtk*" -type d 2>/dev/null | head -3
echo "--- RNAVirHost / phabox ---"
ls /home/zhangwenda/database/virus-db 2>/dev/null | grep -iE "phabox|rvh|rnavirhost"
find /home/zhangwenda/biosoft/virus -maxdepth 1 -type d 2>/dev/null | grep -iE "RNAVirHost|phabox|RVH"
echo "--- viroid detect / VirBot ---"
ls /home/zhangwenda/biosoft/VirBot 2>/dev/null | head -5
grep -riE "version" /home/zhangwenda/biosoft/VirBot/VirBot.py 2>/dev/null | head -2
