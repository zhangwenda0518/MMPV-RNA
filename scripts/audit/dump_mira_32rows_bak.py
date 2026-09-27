#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Dump 33 行（Miraophiovirus/Oleurovirus）在 9/4 备份表里的 Decision_Method / Final_Host，
用于把「ad-hoc 按属名删了 32 行」与「模块撤回后只翻转 23 行」的差额 9 行解释清楚。只读。"""
import csv
import os
import collections

SUM = os.path.expanduser(
    '~/MMPV-paper/onekp-virome/onekp-virus/06_HostPrediction/ensemble_host_summary.tsv')
BAK = SUM + '.bak_bl662_20260911'
GEN = {'Miraophiovirus', 'Oleurovirus'}

with open(BAK, encoding='utf-8', errors='replace', newline='') as fh:
    rd = csv.reader(fh, delimiter='\t', quotechar='"')
    header = next(rd)
    low = {h.lower(): i for i, h in enumerate(header)}
    cols = ['contig_id', 'Family', 'Genus', 'Confidence', 'Host_ICTV', 'Final_Host', 'Decision_Method']
    pick = [low.get(c.lower()) for c in cols]
    print('列索引: %s' % dict(zip(cols, pick)))
    rows = []
    for parts in rd:
        if parts[low['genus']] in GEN:
            rows.append({c: (parts[i] if i is not None and i < len(parts) else '') for c, i in zip(cols, pick)})

print('命中 %d 行' % len(rows))
ct = collections.Counter((r['Genus'], r['Final_Host'], r['Decision_Method']) for r in rows)
for k, v in sorted(ct.items()):
    print('  Genus=%-15s Final_Host=%-8s Decision_Method=%-16s -> %d' % (k[0], k[1], k[2], v))

print()
print('逐行（前 40 行）:')
for r in rows[:40]:
    print('  %s | %s | %s | ICTV=%s | %s | %s' % (
        r['contig_id'][:60], r['Family'], r['Genus'], r['Host_ICTV'], r['Final_Host'], r['Decision_Method']))
