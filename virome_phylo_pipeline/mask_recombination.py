#!/usr/bin/env python3
"""mask_recombination.py — 重组区 mask + 稳健性验证（通用）

RDP5 检出重组断点后, mask 掉重组区间, 重算关键统计并对比:
  π / θ / Tajima's D / 变异位点 S / Rm / 按基因 Ka-Ks

用法:
  python -m mask_recombination \
      --fasta aln.fasta --rdp5 RDP5结果.csv [--genes genes.csv] \
      [--min-methods 3] [--mask-char N] [--outdir out]

mask 字符: N (保留位点, 推荐, 可进 BEAST) | - (删列)
"""
import argparse, csv, json, os, re, sys

# ── sys.path 桥接: 本脚本在项目根目录, 依赖模块在 utils/ ──
_here = os.path.dirname(os.path.abspath(__file__))
_utils = os.path.join(_here, "utils")
if _here not in sys.path:
    sys.path.insert(0, _here)
if _utils not in sys.path:
    sys.path.insert(0, _utils)

from pathlib import Path
from collections import Counter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import dnasp
from utils.rdp5_validate import parse_rdp5_csv


def merge_intervals(intervals, aln_len):
    """合并区间 (含跨末端拆分), 返回 [(start, end)] 1-based"""
    segs = []
    for b, e in intervals:
        b, e = max(1, b), min(aln_len, e)
        if b <= e:
            segs.append((b, e))
        else:  # 跨末端: b→len + 1→e
            segs.append((b, aln_len))
            if e >= 1:
                segs.append((1, e))
    segs.sort()
    merged = []
    for b, e in segs:
        if merged and b <= merged[-1][1] + 1:
            merged[-1] = (merged[-1][0], max(merged[-1][1], e))
        else:
            merged.append((b, e))
    return merged


def mask_alignment(fasta_path, intervals, mask_char='N', out_path=None):
    """mask 比对: 区间位点 → mask_char"""
    aln = dnasp.load_alignment(Path(fasta_path))
    seqs = [list(s) for s in aln.seqs]
    n_masked = 0
    for b, e in intervals:
        for i in range(b - 1, min(e, len(seqs[0]))):
            for s in seqs:
                if s[i] not in ('-', 'N'):
                    s[i] = mask_char
            n_masked += 1
    masked = [''.join(s) for s in seqs]
    if out_path:
        # 2026-09-15: 补 encoding
        with open(out_path, 'w', encoding='utf-8') as f:
            for nm, s in zip(aln.names, masked):
                f.write(f'>{nm}\n{s}\n')
    pct = n_masked / len(seqs[0]) * 100
    return masked, n_masked, pct


def stats(aln):
    """dnasp 统计: 多态 + 重组"""
    try:
        res = dnasp.run_analysis(aln, analyses={'polymorphism', 'recombination'})
        poly = res['global']
        recomb = res['recombination']
        return {
            'n': poly.n, 'L': poly.L_net, 'S': poly.S,
            'pi': round(poly.Pi, 5), 'theta': round(poly.ThetaW, 5),
            'tajima_d': round(poly.TajimaD, 4),
            'rm': recomb.Rm, 'n_incompatible_pairs': recomb.n_incompatible_pairs,
        }
    except Exception as e:
        return {'error': str(e)}


def kaks_by_gene(aln, genes_csv, reference=None):
    """按基因 Ka/Ks (复用 dnasp_bridge 逻辑)

    reference: 坐标基准序列名; 缺省取第一条并告警 (2026-09-15 新增)
    """
    import csv as _csv
    if not genes_csv or not os.path.exists(genes_csv):
        return {}
    genes = []
    with open(genes_csv) as f:
        for row in _csv.DictReader(f):
            genes.append((row['gene'], int(row['start']), int(row['end'])))
    seen, genes_u = set(), []
    for g in genes:
        if g[0] not in seen:
            seen.add(g[0]); genes_u.append(g)
    seqs = list(zip(aln.names, aln.seqs))
    # 2026-09-15 修复: 旧版无条件 `ref_seq = seqs[0][1]` —— MAFFT 输出顺序不保证
    # 参考在首位, 基准漂移会让全部基因坐标整体错位且无提示。
    from utils.seq_match import resolve_reference
    _ref_name, ref_seq = resolve_reference(
        seqs, reference, who="mask_recombination.kaks_by_gene")
    ref_pos = {}
    seen_n = 0
    for i, c in enumerate(ref_seq):
        if c != '-':
            seen_n += 1
            ref_pos[seen_n] = i
    out = {}
    for g, s, e in genes_u:
        if s not in ref_pos or e not in ref_pos:
            continue
        cs, ce = ref_pos[s], ref_pos[e] + 1
        gs = [x[cs:ce] for _, x in seqs]
        nc = (ce - cs) // 3 * 3
        gs_c = [x[:nc] for x in gs]
        try:
            r = dnasp.compute_ka_ks(gs_c)
            if r.omega is not None:
                out[g] = {'omega': round(r.omega, 4), 'ka': round(r.Ka, 5),
                          'ks': round(r.Ks, 5)}
        except Exception as e:
            print(f'[mask] 警告: {g} Ka/Ks 计算失败: {e}', file=sys.stderr)
    return out


def drop_recombinants(fasta_path, rdp5_csv, min_methods=3, include_parents=False,
                      out_path=None):
    """保守策略: 删除明确的重组序列 (保留其余序列)

    适用: 重组子少 + 目标建树 + 数据量充足
    返回: (剩余序列数, 删除列表)
    """
    import re as _re
    aln = dnasp.load_alignment(Path(fasta_path))
    drop_codes = set()
    parent_codes = set()
    with open(rdp5_csv) as f:
        for line in f:
            m = _re.match(r'\s*(\d+)\s*,\s*(\d+)([~\^\*\$]*)\s*,\s*(\d+)\s*,\s*(\d+)', line)
            if not m:
                continue
            fields = [x.strip() for x in line.split(',')]
            if len(fields) <= 10:
                continue
            methods = sum(1 for i in range(11, 20) if i < len(fields) and fields[i] not in ('', 'NS'))
            if methods < min_methods:
                continue
            rec = fields[8].split('\n')[0].strip().rstrip(',').lstrip('^')
            minor = fields[9].split('\n')[0].strip().rstrip(',')
            major = fields[10].split('\n')[0].strip().rstrip(',')
            if rec:
                drop_codes.add(rec)
            for p in (minor, major):
                if p and not p.startswith('Unknown'):
                    parent_codes.add(p)
    if include_parents:
        drop_codes |= (parent_codes - drop_codes)

    def code_to_idx(code):
        n = 0
        for c in code:
            n = n * 26 + (ord(c) - ord('A') + 1)
        return n - 1

    names = list(aln.names)
    name_set = set(names)
    # 1) 序列名精确匹配 (RDP5 输出真实序列名, 含数字/下划线时)
    drop_names = {c for c in drop_codes if c in name_set}
    # 2) 字母编码 → 索引匹配 (RDP5 内部 A/B/C 编号, 依赖序列顺序)
    if not drop_names:
        for c in drop_codes:
            if c.isalpha():
                idx = code_to_idx(c)
                if 0 <= idx < len(names):
                    drop_names.add(names[idx])
        if drop_codes and not drop_names:
            # 历史坑: 只认字母编码时, 含数字的真实序列名被 c.isalpha() 过滤掉
            # → 静默删 0 条却报完成。现在明确告警。
            print(f'[mask] 警告: 重组子编码无法匹配任何序列 '
                  f'(样本: {sorted(drop_codes)[:5]}...), 未删除任何序列', file=sys.stderr)

    keep, dropped = [], []
    for nm, s in zip(names, aln.seqs):
        if nm in drop_names:
            dropped.append(nm)
        else:
            keep.append((nm, s))

    if out_path:
        # 2026-09-15: 补 encoding
        with open(out_path, 'w', encoding='utf-8') as f:
            for nm, s in keep:
                f.write(f'>{nm}\n{s}\n')
    return len(keep), dropped


def drop_recombinants_by_validation(fasta_path, validation_json, topology='PASSED',
                                    out_path=None):
    """按 validation_summary 的 topology 筛选, 删除重组子

    适用: 拓扑验证已做, 只删除拓扑确认(PASSED)的重组子
    返回: (剩余序列数, 删除列表)
    """
    import re as _re
    aln = dnasp.load_alignment(Path(fasta_path))
    data = json.load(open(validation_json))
    if isinstance(data, dict):
        data = data.get('events', data.get('validation', []))
    drop_codes = []
    for e in data:
        if topology == 'ALL' or e.get('topology') == topology:
            rec = e.get('recombinant', '')
            m = _re.search(r'(CRR\d+|SRR\d+|ERR\d+|DRR\d+)', rec)
            if m:
                drop_codes.append(m.group(1))
    drop_codes = list(dict.fromkeys(drop_codes))
    # 2026-09-15 修复: 旧版用 `acc in nm` 子串判断 → `CRR1` 会命中 `CRR10`,
    # 删错序列。改为「精确 / 独立 token」判定 (与 utils/rdp5_validate.find_seq 同口径),
    # 不再做最宽松的子串兜底。
    from utils.seq_match import is_same_accession
    keep, dropped = [], []
    for nm, s in zip(aln.names, aln.seqs):
        hit = any(is_same_accession(nm, acc) for acc in drop_codes)
        if hit:
            dropped.append(nm)
        else:
            keep.append((nm, s))
    if out_path:
        # 2026-09-15: 补 encoding
        with open(out_path, 'w', encoding='utf-8') as f:
            for nm, s in keep:
                f.write(f'>{nm}\n{s}\n')
    return len(keep), dropped


def main():
    ap = argparse.ArgumentParser(description='重组区 mask + 稳健性验证')
    ap.add_argument('--fasta', required=True, help='比对 FASTA (RDP5 同一输入)')
    ap.add_argument('--rdp5', required=True, help='RDP5 结果 CSV')
    ap.add_argument('--genes', default=None, help='基因注释 CSV')
    ap.add_argument('--min-methods', type=int, default=3,
                    help='仅 mask 方法数≥此值的事件 (默认3)')
    ap.add_argument('--no-cross-end', action='store_true',
                    help='[已默认开启] 排除跨末端事件')
    ap.add_argument('--keep-cross-end', action='store_true',
                    help='保留跨末端事件 (默认排除: 疑似组装 artifact)')
    ap.add_argument('--mask-char', default='N', choices=['N', '-'],
                    help='mask 字符: N 保留位点 (推荐) | - 删列')
    ap.add_argument('--outdir', default='mask_out')
    ap.add_argument('--validation', default=None,
                    help='validation_summary.json (拓扑验证), 按 topology 筛选删除重组子')
    ap.add_argument('--topology', default='PASSED', choices=['PASSED', 'FAILED', 'ALL'],
                    help='删除哪些 topology 的重组子 (默认 PASSED)')
    ap.add_argument('--mask', action='store_true',
                    help='显式选择 mask 区间模式 (默认删除重组子)')
    ap.add_argument('--drop-recombinants', action='store_true',
                    help='[兼容旧参数] 无 --validation 时按方法数删除 (默认已开启)')
    ap.add_argument('--include-parents', action='store_true',
                    help='删除重组子时连同亲本一起删')
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    # 默认: 删除重组子 (优先 validation PASSED, 回退方法数)
    if not args.mask:
        out_drop = f"{args.outdir}/dropped_recombinants.fasta"
        if args.validation:
            n_keep, dropped = drop_recombinants_by_validation(
                args.fasta, args.validation, args.topology, out_path=out_drop)
            print(f'=== 删除重组子 (validation {args.topology}) ===')
        else:
            n_keep, dropped = drop_recombinants(args.fasta, args.rdp5,
                                                min_methods=args.min_methods,
                                                include_parents=args.include_parents,
                                                out_path=out_drop)
            print(f'=== 删除重组子 (方法数 ≥{args.min_methods}) ===')
        aln = dnasp.load_alignment(Path(args.fasta))
        drop_aln = dnasp.load_alignment(Path(out_drop))
        before, after = stats(aln), stats(drop_aln)
        kaks_b, kaks_a = kaks_by_gene(aln, args.genes), kaks_by_gene(drop_aln, args.genes)
        print(f'事件筛选: {args.topology if args.validation else f"方法数 ≥{args.min_methods}"}')
        print(f'删除 {len(dropped)} 条重组序列 → 剩余 {n_keep} 条 '
              f'({n_keep/(n_keep+len(dropped))*100:.0f}%)')
        for d in dropped:
            print(f'  删除: {d}')
        print()
        hdr = f"{'指标':<14}{'删前':>12}{'删后':>12}{'变化':>10}"
        print(hdr); print('-' * len(hdr))
        for k, label in [('n', '序列数'), ('S', '变异位点'), ('pi', 'π'),
                         ('theta', 'θ'), ('tajima_d', "Tajima's D"), ('rm', 'Rm')]:
            if k in before and k in after:
                b, a = before[k], after[k]
                chg = f"{a - b:+.4f}" if isinstance(b, float) else f"{a - b:+d}"
                print(f"{label:<14}{b:>12}{a:>12}{chg:>10}")
        print(f'\n按基因 Ka/Ks (删前/删后):')
        for g in kaks_b:
            bw = kaks_b[g]['omega']
            aw = kaks_a.get(g, {}).get('omega', 'NA')
            print(f"{g:<6}{bw:>10}{aw:>10}")
        # 2026-09-15: 同上, 句柄泄漏 + 缺 encoding
        with open(f"{args.outdir}/drop_report.json", 'w', encoding='utf-8') as _rf:
            json.dump({'mode': 'drop', 'dropped': dropped, 'n_keep': n_keep,
                       'before': before, 'after': after,
                       'kaks_before': kaks_b, 'kaks_after': kaks_a},
                      _rf, indent=2, default=str)
        print(f'\n删除后比对: {out_drop}')
        print(f'报告: {args.outdir}/drop_report.json')
        print('\n提示: 用删除后比对重跑 IQ-TREE/BEAST 建立可靠系统树')
        return 0

    aln = dnasp.load_alignment(Path(args.fasta))
    aln_len = aln.L
    events = parse_rdp5_csv(args.rdp5)

    # 筛选事件 (默认排除跨末端)
    sel = [ev for ev in events if ev['n_methods'] >= args.min_methods]
    if not args.keep_cross_end:
        n_cross = sum(1 for ev in sel if ev['bp_start'] > ev['bp_end'])
        sel = [ev for ev in sel if ev['bp_start'] <= ev['bp_end']]
        print(f'[默认] 跳过 {n_cross} 个跨末端事件 (疑似组装 artifact, --keep-cross-end 保留)')
    intervals = [(ev['bp_start'], ev['bp_end']) for ev in sel]
    merged = merge_intervals(intervals, aln_len)

    os.makedirs(args.outdir, exist_ok=True)
    out_fasta = f"{args.outdir}/masked_{args.mask_char}.fasta"
    _, n_masked, pct = mask_alignment(args.fasta, merged, args.mask_char, out_fasta)

    # 统计对比
    masked_aln = dnasp.load_alignment(Path(out_fasta))
    before, after = stats(aln), stats(masked_aln)
    kaks_b, kaks_a = kaks_by_gene(aln, args.genes), kaks_by_gene(masked_aln, args.genes)

    print(f'=== 重组区 mask 稳健性验证 ===')
    print(f'事件: {len(events)} 总 / {len(sel)} 高置信 (≥{args.min_methods} 方法)')
    print(f'mask 区间 (合并): {merged}')
    print(f'mask 位点数: {n_masked} ({pct:.1f}% 基因组)\n')

    hdr = f"{'指标':<14}{'mask前':>12}{'mask后':>12}{'变化':>10}"
    print(hdr)
    print('-' * len(hdr))
    for k, label in [('n', '序列数'), ('S', '变异位点'), ('pi', 'π'),
                     ('theta', 'θ'), ('tajima_d', "Tajima's D"), ('rm', 'Rm')]:
        if k in before and k in after:
            b, a = before[k], after[k]
            chg = f"{a - b:+.4f}" if isinstance(b, float) else f"{a - b:+d}"
            print(f"{label:<14}{b:>12}{a:>12}{chg:>10}")

    print(f"\n按基因 Ka/Ks (mask 前后):")
    print(f"{'基因':<6}{'mask前ω':>10}{'mask后ω':>10}")
    for g in kaks_b:
        bw = kaks_b[g]['omega']
        aw = kaks_a.get(g, {}).get('omega', 'NA')
        print(f"{g:<6}{bw:>10}{aw:>10}")

    # 保存报告
    report = {'min_methods': args.min_methods, 'events_total': len(events),
              'events_masked': len(sel), 'mask_intervals': merged,
              'masked_nt': n_masked, 'masked_pct': round(pct, 1),
              'before': before, 'after': after,
              'kaks_before': kaks_b, 'kaks_after': kaks_a}
    with open(f"{args.outdir}/mask_report.json", 'w') as f:
        json.dump(report, f, indent=2, default=str)
    print(f"\nmask 后比对: {out_fasta}")
    print(f"报告: {args.outdir}/mask_report.json")
    return 0


if __name__ == '__main__':
    sys.exit(main())
