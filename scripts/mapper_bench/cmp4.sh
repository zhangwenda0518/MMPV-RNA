#!/bin/bash
# 实证检查 onekp 失败样本的输入 fasta：重复名、目标 read、双端文件
d=/home/zhangwenda/data-test/out/00b_HostDepletion
echo "== ERR2040134 相关文件 =="
ls "$d" | grep ERR2040134
echo
echo "== 目标 read 名是否出现及次数 =="
zcat "$d/ERR2040134_clean_1.fa.gz" 2>/dev/null | grep '^>' | grep -c 'ERR2040134.10268640' 
echo
echo "== R1 总条数 vs 唯一名数（查重复）=="
zcat "$d/ERR2040134_clean_1.fa.gz" | grep '^>' | wc -l
zcat "$d/ERR2040134_clean_1.fa.gz" | grep '^>' | sort -u | wc -l
echo
echo "== 头 2 条 read 名 =="
zcat "$d/ERR2040134_clean_1.fa.gz" | grep '^>' | head -2
