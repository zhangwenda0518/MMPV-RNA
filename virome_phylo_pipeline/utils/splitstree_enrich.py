#!/usr/bin/env python3
"""
splitstree_enrich.py -- 把 NeighborNet nexus 增强为"GUI-ready"版本, 供 SplitsTree 导入
增强项 (全部可选, 组合使用):
  1. --groups metadata.tsv       分组着色: 写入 BEGIN SETS; TAXSET 块 (SplitsTree 里一键按组着色)
  2. --bootstrap N --alignment x.fasta
                                 bootstrap 支持值: 对比对列重采样, 重跑 NeighborNet N 次,
                                 支持值写入 FORMAT CONFIDENCES=YES 的第二个数值列
  3. --filter planar|weakly|2tree|top:K|minsize:N|minconf:F
                                 过滤: planar=cycle 连续性(平面可画) / weakly=弱兼容贪心 /
                                 2tree=两棵树贪心 / top:K=前K强 / minsize:N / minconf:F

用法示例:
  python3 splitstree_enrich.py -i neighbornet.nexus -o enriched.nexus \
      --groups groups.tsv --alignment aln.fasta --bootstrap 100 \
      --filter minconf:0.5 --min-weight 0
groups.tsv: 两列 TSV (taxon_id<TAB>group), 无表头; 未列出的 taxon 归 "Others"
"""
import argparse
import itertools
import math
import os
import random
import re
import subprocess
import sys
import tempfile

import splitstree_bridge  # 用于 bootstrap 重跑


# ---------- nexus 解析 ----------
def parse_nexus(path):
    text = open(path, encoding='utf-8', errors='ignore').read()
    labels = re.findall(r"\[\d+\]\s+'([^']+)'", text)
    m = re.search(r'CYCLE\s+([\d\s]+?);', text, re.S)
    cycle = [int(x) for x in m.group(1).split()]
    mb = re.search(r'^BEGIN SPLITS;(.*?)^END;', text, re.S | re.M)
    has_conf_hdr = re.search(r'confidences\s*=\s*yes', mb.group(1), re.I) is not None
    splits = []
    for line in mb.group(1).splitlines():
        line = line.strip()
        if has_conf_hdr:
            mm = re.match(r'\[\d+,\s*size=\d+\]\s+([\d.eE+-]+)\s+([\d.eE+-]+)\s+(.*?),?\s*$', line)
            if not mm:
                continue
            w, conf, taxa_txt = float(mm.group(1)), float(mm.group(2)), mm.group(3)
        else:
            mm = re.match(r'\[\d+,\s*size=\d+\]\s+([\d.eE+-]+)\s+(.*?),?\s*$', line)
            if not mm:
                continue
            w, conf, taxa_txt = float(mm.group(1)), None, mm.group(2)
        taxa = frozenset(int(x) for x in taxa_txt.split())
        splits.append({'w': w, 'conf': conf, 'taxa': taxa})
    fit = None
    m = re.search(r'fit=([\d.]+)', text)
    if m:
        fit = float(m.group(1))
    return labels, cycle, splits, fit


def canonical(taxa, ntax):
    """split 的规范表示: 取较小的一侧"""
    if len(taxa) * 2 <= ntax:
        return frozenset(taxa)
    return frozenset(range(1, ntax + 1)) - taxa


# ---------- 过滤器 ----------
def cycle_positions(cycle):
    return {t: k for k, t in enumerate(cycle)}


def is_planar_split(taxa, pos, n):
    """NeighborNet 是 circular split system: split 平面可画 <=> 某一侧在 cycle 上连续"""
    side = taxa if len(taxa) * 2 <= n else set(range(1, n + 1)) - taxa
    if not side or len(side) == n:
        return False
    ps = sorted(pos[t] for t in side)
    gaps = [(ps[(i + 1) % len(ps)] - ps[i]) % n for i in range(len(ps))]
    return sum(1 for g in gaps if g != 1) <= 1 and len(ps) < n


def compatible(a, b):
    """两个 split 兼容: 交集为空 或 一方包含另一方"""
    i = a & b
    return not i or i == a or i == b


def weakly_compatible(new, kept):
    """弱兼容 (三 split 条件): 对任意 kept 中两个 split A,B 与 new C:
       A∩B∩C, A∩B∩C^c, A∩B^c∩C, A^c∩B∩C 至多一个非空"""
    U = None
    for A, B in itertools.combinations(kept, 2):
        i1 = (A & B) & new
        i2 = (A & B) - new
        i3 = (A & new) - B
        i4 = (B & new) - A
        if sum(1 for x in (i1, i2, i3, i4) if x) > 1:
            return False
    return True


def apply_filter(splits, flt, cycle, ntax, support):
    pos = cycle_positions(cycle)
    splits = sorted(splits, key=lambda s: -s['w'])  # 贪心类过滤器需权重降序
    # 规范化取小侧 (兼容/弱兼容条件依赖侧的选择)
    for s in splits:
        s['taxa'] = canonical(s['taxa'], ntax)
    if flt == 'planar':
        return [s for s in splits if is_planar_split(s['taxa'], pos, ntax)]
    if flt == 'weakly':
        kept = []
        for s in splits:
            if weakly_compatible(s['taxa'], [k['taxa'] for k in kept]):
                kept.append(s)
        return kept
    if flt == '2tree':
        t1, t2, rest = [], [], []
        for s in splits:
            if all(compatible(s['taxa'], k['taxa']) for k in t1):
                t1.append(s)
            elif all(compatible(s['taxa'], k['taxa']) for k in t2):
                t2.append(s)
        return t1 + t2
    if flt.startswith('top:'):
        k = int(flt[4:])
        return splits[:k]
    if flt.startswith('minsize:'):
        n = int(flt[8:])
        return [s for s in splits if min(len(s['taxa']), ntax - len(s['taxa'])) >= n]
    if flt.startswith('minconf:'):
        f = float(flt[8:])
        return [s for s in splits if (s['conf'] or 0) >= f]
    raise ValueError('unknown filter: ' + flt)


# ---------- bootstrap ----------
def bootstrap_support(aln_fasta, nrep, tmpdir, max_taxa=200, log=print):
    """列重采样 + 重跑 NeighborNet, 返回 {canonical split -> support 0..1}"""
    # 读比对
    names, seqs = [], []
    cur = None
    for line in open(aln_fasta, encoding='utf-8', errors='replace'):
        line = line.strip()
        if line.startswith('>'):
            names.append(line[1:].split()[0])
            seqs.append([])
        elif seqs:
            seqs[-1].append(line)
    seqs = [''.join(s) for s in seqs]
    L = max(len(s) for s in seqs)
    counts = {}
    rng = random.Random(42)
    for rep in range(nrep):
        cols = [rng.randrange(L) for _ in range(L)]
        rs = [''.join(s[c] for c in cols if c < len(s)) for s in seqs]
        fa = os.path.join(tmpdir, 'bs_%d.fasta' % rep)
        with open(fa, 'w', encoding='utf-8') as f:
            for n, s in zip(names, rs):
                f.write('>%s\n%s\n' % (n, s))
        outd = os.path.join(tmpdir, 'bs_out_%d' % rep)
        argv_bak = sys.argv
        sys.argv = ['splitstree_bridge.py', '-i', fa, '-o', outd,
                    '--network', 'neighbornet', '--max-taxa', str(max_taxa)]
        try:
            import contextlib, io
            with contextlib.redirect_stdout(io.StringIO()):
                rc = splitstree_bridge.main()
        except SystemExit as e:
            rc = e.code
        finally:
            sys.argv = argv_bak
        if rc:
            log('  bootstrap rep %d failed (rc=%s)' % (rep + 1, rc))
            continue
        _, cyc, sp, _ = parse_nexus(os.path.join(outd, 'neighbornet.nexus'))
        nt = len(cyc)
        seen = set()
        for s in sp:
            cs = canonical(s['taxa'], nt)
            if 0 < len(cs) < nt:
                seen.add(cs)
        for cs in seen:
            counts[cs] = counts.get(cs, 0) + 1
        log('  bootstrap %d/%d done (%d splits)' % (rep + 1, nrep, len(seen)))
    n_ok = nrep  # 简化: 按总 rep 数计
    return {k: v / n_ok for k, v in counts.items()}


# ---------- 写增强 nexus ----------
def write_nexus(path, labels, cycle, splits, fit, groups=None):
    ntax = len(cycle)
    with open(path, 'w', encoding='utf-8') as f:
        f.write('#NEXUS\n\n')
        f.write('BEGIN TAXA;\nDIMENSIONS ntax=%d;\nTAXLABELS\n' % ntax)
        for i, lab in enumerate(labels, 1):
            f.write("\t[%d] '%s'\n" % (i, lab))
        f.write(';\nEND; [Taxa]\n\n')
        if groups:
            f.write('BEGIN SETS;\n')
            for g, members in groups.items():
                ids = sorted(i for i, t in enumerate(cycle, 1)
                             if labels[t - 1] in members)
                if ids:
                    f.write('\tTAXSET %s = %s;\n' % (g, ' '.join(map(str, ids))))
            f.write('END; [SETS]\n\n')
        has_conf = any(s['conf'] is not None for s in splits)
        f.write('BEGIN SPLITS;\n')
        f.write('DIMENSIONS ntax=%d nsplits=%d;\n' % (ntax, len(splits)))
        f.write('FORMAT labels=no weights=yes confidences=%s;\n' %
                ('yes' if has_conf else 'no'))
        if fit is not None:
            f.write('PROPERTIES fit=%s;\n' % fit)
        f.write('CYCLE %s;\n' % ' '.join(map(str, cycle)))
        f.write('MATRIX\n')
        for i, s in enumerate(splits, 1):
            side = sorted(canonical(s['taxa'], ntax))
            if has_conf:
                f.write('[%d, size=%d]\t%.8g\t%s\t%s,\n' %
                        (i, len(side), s['w'], '%.4f' % (s['conf'] or 0),
                         ' '.join(map(str, side))))
            else:
                f.write('[%d, size=%d]\t%.8g\t%s,\n' %
                        (i, len(side), s['w'], ' '.join(map(str, side))))
        f.write(';\nEND; [Splits]\n')
        # SplitsTree 打开即用等角布局 (SplitsTree4 方言, SplitsTree6 忽略亦无害)
        f.write('\nBEGIN st_Assumptions;\nuptodate;\n'
                'splitstransform=EqualAngle UseWeights=true RunConvexHull=true '
                'DaylightIterations=0 OptimizeBoxesIterations=0 '
                'SpringEmbedderIterations=0;\n'
                'SplitsPostProcess filter=none;\nEND; [st_Assumptions]\n')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('-i', '--input', required=True)
    ap.add_argument('-o', '--output', required=True)
    ap.add_argument('--groups', help='两列 TSV: taxon_id<TAB>group')
    ap.add_argument('--alignment', help='bootstrap 用的比对 fasta')
    ap.add_argument('--bootstrap', type=int, default=0)
    ap.add_argument('--filter', action='append', default=[],
                    help='可叠加: planar / weakly / 2tree / top:K / minsize:N / minconf:F')
    ap.add_argument('--min-weight', type=float, default=0.0)
    a = ap.parse_args()

    labels, cycle, splits, fit = parse_nexus(a.input)
    ntax = len(cycle)
    print('loaded: %d taxa, %d splits, fit=%s' % (ntax, len(splits), fit))

    # 1. bootstrap
    if a.bootstrap:
        if not a.alignment:
            sys.exit('--bootstrap 需要 --alignment')
        with tempfile.TemporaryDirectory(prefix='st6_bs_') as td:
            support = bootstrap_support(a.alignment, a.bootstrap, td)
        # 应用: 支持值 = 严格匹配 canonical split
        for s in splits:
            cs = canonical(s['taxa'], ntax)
            s['conf'] = support.get(cs, 0.0)
        print('bootstrap %d reps -> support assigned' % a.bootstrap)

    # 2. 过滤 (可叠加, 依次应用)
    for flt in a.filter:
        before = len(splits)
        splits = apply_filter(splits, flt, cycle, ntax, None)
        print('filter %-10s: %d -> %d splits' % (flt, before, len(splits)))
    if a.min_weight > 0:
        splits = [s for s in splits if s['w'] >= a.min_weight]
        print('min-weight %.3g: %d left' % (a.min_weight, len(splits)))

    # 3. 分组
    groups = None
    if a.groups:
        member = {}
        for line in open(a.groups, encoding='utf-8', errors='replace'):
            p = line.rstrip('\n').split('\t')
            if len(p) >= 2 and p[0]:
                member[p[0]] = p[1]
        groups = {}
        for lab in labels:
            groups.setdefault(member.get(lab, 'Others'), set()).add(lab)
        print('groups: ' + ', '.join('%s=%d' % (g, len(m)) for g, m in groups.items()))

    write_nexus(a.output, labels, cycle, splits, fit, groups)
    print('saved:', a.output)


if __name__ == '__main__':
    main()
