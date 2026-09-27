#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""只读：复算科-属相容性（口径同 virus_classifier_analysis.R::enforce_rank_containment）。

判据：某行的 Family 与「该行 Genus 在跨参照表里 NCBI 侧定型科」和「VMR 侧定型科」都不一致
      = 双侧不相容；只有一侧不一致 = 单侧不相容（设计上保留，含分类学版本漂移）。

覆盖 5 张表：枸杞 05 原表 / 枸杞 05 校准表 / OneKP 05 原表 / 枸杞下游植物表 / OneKP 下游植物表。
不写任何文件，只打印统计。args 可传额外表路径。
"""
import os
import sys
import collections

REF = '/home/zhangwenda/database/taxonomy/genus_family_ref.tsv'
G = '/home/zhangwenda/MMPV-paper/goji-virome/02_novel_virus/RNA-Lycium_barbarum_out'
O = '/home/zhangwenda/MMPV-paper/onekp-virome/onekp-virus'

TABLES = [
    ('GOJI 05 原表 08-12', os.path.join(G, '05_Taxonomy/Votus.integrated/final_integrated_classification.tsv')),
    ('GOJI 05 校准表 09-14', os.path.join(
        G, '05_Taxonomy/Votus.integrated/calibration_20260914/final_integrated_classification.calibrated.tsv')),
    ('ONEKP 05 原表 09-04', os.path.join(O, '05_Taxonomy/Votus.integrated/final_integrated_classification.tsv')),
    ('GOJI 下游植物表', os.path.join(G, '10_Reports/All_plant.viruses_info.tsv')),
    ('ONEKP 下游植物表', os.path.join(O, '10_Reports/All_plant.viruses_info.tsv')),
]


def norm(x):
    if x is None:
        return None
    x = x.strip().strip('"').replace('*', '').strip()
    if x == '' or x.upper() in ('NA', 'N/A', '-'):
        return None
    return x.lower()


def load_ref():
    ref = {}
    with open(REF, 'r', encoding='utf-8', errors='replace') as f:
        hdr = f.readline().rstrip('\n').split('\t')
        idx = {h.strip(): i for i, h in enumerate(hdr)}
        gi, ni, vi = idx['Genus'], idx['NCBI_Family'], idx['VMR_Family']
        for line in f:
            p = line.rstrip('\n').split('\t')
            if len(p) <= max(gi, ni, vi):
                continue
            g = norm(p[gi])
            if g is None:
                continue
            nf = norm(p[ni])
            vraw = norm(p[vi])
            e = ref.setdefault(g, [set(), set()])
            if nf:
                e[0].add(nf)
            if vraw:
                for t in vraw.split(';'):
                    t = t.strip()
                    if t:
                        e[1].add(t)
    return ref


def scan(label, path, ref):
    if not os.path.isfile(path):
        print('%-24s 文件不存在: %s' % (label, path))
        return
    st = os.stat(path)
    cnt = collections.Counter()
    bad = []
    n = 0
    with open(path, 'r', encoding='utf-8', errors='replace') as f:
        hdr = f.readline().rstrip('\n').split('\t')
        idx = {h.strip().strip('"'): i for i, h in enumerate(hdr)}
        ci, fi, gi = idx.get('contig_id'), idx.get('Family'), idx.get('Genus')
        if fi is None or gi is None:
            print('%-24s 缺 Family/Genus 列，跳过（表头: %s）' % (label, '|'.join(hdr[:12])))
            return
        for line in f:
            p = line.rstrip('\n').split('\t')
            if len(p) <= max(fi, gi):
                continue
            n += 1
            fam = norm(p[fi])
            gen = norm(p[gi])
            if gen is None:
                cnt['Genus空'] += 1
                continue
            if fam is None:
                cnt['Family空_Genus有值'] += 1
                continue
            if gen not in ref:
                cnt['属无跨参照记录'] += 1
                continue
            nfams, vfams = ref[gen]
            nb = bool(nfams) and fam not in nfams
            vb = bool(vfams) and fam not in vfams
            if nb and vb:
                cnt['双侧不相容'] += 1
                if len(bad) < 5:
                    bad.append((p[ci] if ci is not None else '?', fam, gen,
                                '|'.join(sorted(nfams)) or '-', '|'.join(sorted(vfams)) or '-'))
            elif nb or vb:
                cnt['单侧不相容'] += 1
            else:
                cnt['相容'] += 1
    print('%-24s 行数 %-8d 双侧不相容 %-7d 单侧 %-7d 属无参照 %-7d Genus空 %-6d 相容 %-7d (mtime %s)' % (
        label, n, cnt['双侧不相容'], cnt['单侧不相容'], cnt['属无跨参照记录'],
        cnt['Genus空'], cnt['相容'], __import__('time').strftime('%m-%d %H:%M', __import__('time').localtime(st.st_mtime))))
    for r in bad:
        print('      样例: %s\n        Fam=%s Genus=%s NCBI定型=%s VMR定型=%s' % r)


def main():
    ref = load_ref()
    print('参照表属数 %d  (%s)' % (len(ref), REF))
    print('-' * 150)
    for label, path in TABLES:
        scan(label, path, ref)
    for extra in sys.argv[1:]:
        scan(os.path.basename(extra), extra, ref)


main()
