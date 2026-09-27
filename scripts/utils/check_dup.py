#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""验证 R1/R2 read 名是否重叠"""
import gzip
r1n, r2n = set(), set()
with gzip.open('/home/zhangwenda/data-test/out/00b_HostDepletion/ERR2040134_clean_1.fa.gz','rt') as f:
    for line in f:
        if line.startswith('>'):
            r1n.add(line.split()[0][1:])
with gzip.open('/home/zhangwenda/data-test/out/00b_HostDepletion/ERR2040134_clean_2.fa.gz','rt') as f:
    for line in f:
        if line.startswith('>'):
            r2n.add(line.split()[0][1:])
print('R1 unique names:', len(r1n))
print('R2 unique names:', len(r2n))
ov = r1n & r2n
print('overlap R1&R2:', len(ov))
print('sample overlap:', list(ov)[:5])
