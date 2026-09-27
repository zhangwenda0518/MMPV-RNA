#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""对账 Miraophiovirus / Oleurovirus 在三份产物里的行数口径差异。

背景：三个数字互相看着不一致，必须说清各自在数什么。
  A) `ensemble_host_summary.tsv` 的 .bak 里 Miraophiovirus 29 行、现表 1 行 -> 删了 28 行
  B) 本轮 host prediction 逐行对拍：撤回两属后 onekp Plant 净 +23（Miraophiovirus 20 + Oleurovirus 3）
  C) 05 表里 Genus 等于这两个名字的行数（未知）

本脚本只读，不写任何产物。输出每份表的口径定义与计数。
"""
import os
import sys
import collections

ONEKP = os.path.expanduser('~/MMPV-paper/onekp-virome/onekp-virus')
T05 = os.path.join(ONEKP, '05_Taxonomy/Votus.integrated/final_integrated_classification.tsv')
C9 = os.path.join(ONEKP, '06_HostPrediction/C9_ICTV_result/classification_result.tsv')
SUM = os.path.join(ONEKP, '06_HostPrediction/ensemble_host_summary.tsv')
SUMBAK = SUM + '.bak_bl662_20260911'
SUM_BAKS = [
    SUM + '.bak_bl662_20260911',
    SUM + '.bak_ghost_20260912',
    SUM + '.bak_ghost2_20260912',
]
GEN = {'Miraophiovirus', 'Oleurovirus'}


def load(path, wanted):
    """返回 (header, idx, rows)。用 csv 模块解析，字段带引号（05 表全字段带 "）时必须正确去引号。"""
    import csv
    with open(path, encoding='utf-8', errors='replace', newline='') as fh:
        rd = csv.reader(fh, delimiter='\t', quotechar='"')
        header = next(rd)
        low = {h.lower(): i for i, h in enumerate(header)}
        idx = {w: low.get(w.lower()) for w in wanted}
        data = []
        for parts in rd:
            data.append({k: (parts[i] if i is not None and i < len(parts) else '') for k, i in idx.items()})
    return header, idx, data


def report(path, wanted, keycol, hostcol=None):
    if not os.path.exists(path):
        print('  [缺失] %s' % path)
        return None
    header, idx, data = load(path, wanted)
    print('  文件: %s' % os.path.basename(path))
    print('  行数(不含表头): %s' % len(data))
    print('  列映射: %s' % idx)
    hit = [r for r in data if r[keycol] in GEN]
    print('  Genus/名字列命中 %d 行' % len(hit))
    cnt = collections.Counter(r[keycol] for r in hit)
    for k, v in cnt.most_common():
        print('    %s: %d' % (k, v))
    if hostcol and idx.get(hostcol) is not None:
        ct = collections.Counter((r[keycol], r[hostcol]) for r in hit)
        for k, v in sorted(ct.items()):
            print('    %s x %s = %d' % (k[0], repr(k[1]), v))
    return hit


print('=== A) ensemble_host_summary.tsv（现表 + 3 个备份）===')
print('[注意] 05 表与 summary 的字段带双引号，本脚本用 csv 模块解析；'
      '若用 split(\\t) 会得到带引号的值而误判为 0 命中。')
for p in [SUM] + SUM_BAKS:
    report(p, ['contig_id', 'Family', 'Genus', 'Final_Host', 'Decision_Method'], 'Genus', 'Final_Host')

print()
print('=== B) 05 整合表（Family/Genus 口径）===')
report(T05, ['contig_id', 'Family', 'Genus', 'Order'], 'Genus')

print()
print('=== C) C9 分类表（Predicted_Host 口径）===')
if os.path.exists(C9):
    header, idx, data = load(C9, ['contig_id', 'Predicted_Host'])
    print('  列映射: %s' % idx)
    print('  行数: %d' % len(data))
    if idx.get('Predicted_Host') is not None:
        print('  Predicted_Host 取值: %s' % sorted({r['Predicted_Host'] for r in data}))
else:
    print('  [缺失] %s' % C9)

print()
print('=== D) 三表 contig 交集：Genus in GEN 的行，其 C9 Predicted_Host / 05 Family 分布 ===')
if os.path.exists(T05) and os.path.exists(C9):
    _, _, t05 = load(T05, ['contig_id', 'Family', 'Genus'])
    _, _, c9 = load(C9, ['contig_id', 'Predicted_Host'])
    host = {r['contig_id']: r['Predicted_Host'] for r in c9}
    fam = {r['contig_id']: r['Family'] for r in t05}
    hit = [r for r in t05 if r['Genus'] in GEN]
    ct = collections.Counter((r['Genus'], r['Family'], host.get(r['contig_id'], '<无C9行>')) for r in hit)
    for k, v in sorted(ct.items()):
        print('  %s | %s | C9=%s -> %d' % (k[0], k[1], k[2], v))
    print('  合计 %d 行（05 表 Genus 命中）' % len(hit))
    nplant = sum(1 for r in hit if host.get(r['contig_id']) == 'Plant')
    print('  其中 C9 判 Plant 的 %d 行' % nplant)
