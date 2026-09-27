#!/usr/bin/env python3
"""align_qc.py — 重组分析前序列比对 + 质量检查 + 不完整序列去除（通用）

MAFFT 比对 → 逐序列质量检查（长度/gap%/N%/与参考 identity）→ 去除不完整序列 → QC 报告

用法:
  python -m utils.align_qc \
      --input seqs.fasta [--reference REF_ID] [--outdir out] \
      [--max-gap 0.10] [--max-n 0.05] [--min-length 0.90] \
      [--skip-mafft] [--max-identity-drop 0.02]

质量检查维度:
  length    序列长度 / 参考长度 (不完整 = 长度过短)
  gap%      比对后 gap 比例 (比对质量差 = gap 过多)
  N%        简并字符比例 (低质量 = N 过多)
  identity  与参考的 pairwise identity (异常 = 偏离群体)

输出:
  out/align_qc_report.tsv   每序列质量指标 + 判定 + 去除理由
  out/clean.fasta           保留序列 (比对后)
  out/removed.fasta         被去除序列 (原始)
  out/alignment.mafft.fasta MAFFT 比对结果
"""
import argparse
import csv
import os
import subprocess
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def run_mafft(fasta_path, out_path, threads=8):
    """MAFFT 比对"""
    cmd = ['mafft', '--auto', '--thread', str(threads), '--quiet', fasta_path]
    # 2026-09-15: 补 encoding —— Windows 默认 cp936, 路径/内容含中文即乱码或报错
    with open(out_path, 'w', encoding='utf-8') as f:
        subprocess.run(cmd, stdout=f, check=True)
    return out_path


def pairwise_identity(s1, s2):
    """SDT 公式: (1 - dist/denom) * 100, 双方非 gap 才计入"""
    dist, gaps = 0, 0
    for a, b in zip(s1, s2):
        if a != '-' and b != '-':
            if a != b:
                dist += 1
        else:
            gaps += 1
    denom = len(s1) - gaps
    return (1 - dist / denom) * 100 if denom > 0 else 0.0


def main():
    ap = argparse.ArgumentParser(description='比对 + 质量检查 + 不完整序列去除')
    ap.add_argument('--input', required=True, help='输入 fasta (比对或未比对)')
    ap.add_argument('--reference', default=None, help='参考序列 ID (长度基准 + identity 检查)')
    ap.add_argument('--outdir', default='align_qc_out')
    ap.add_argument('--max-gap', type=float, default=0.10, help='最大 gap 比例 (默认 0.10)')
    ap.add_argument('--max-n', type=float, default=0.05, help='最大 N 比例 (默认 0.05)')
    ap.add_argument('--min-length', type=float, default=0.90, help='最小长度/参考 (默认 0.90)')
    ap.add_argument('--max-identity-drop', type=float, default=0.02,
                    help='与参考 identity 允许的最大偏离 (默认 0.02)')
    ap.add_argument('--skip-mafft', action='store_true', help='输入已比对, 跳过 MAFFT')
    ap.add_argument('--threads', type=int, default=8)
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # 1. 读序列
    recs = list(SeqIO.parse(args.input, 'fasta'))
    print(f'输入: {len(recs)} 条序列')
    seq_lens = {r.id: len(str(r.seq).replace("-", "")) for r in recs}

    # 2. 参考长度基准
    if args.reference:
        ref_rec = next((r for r in recs if args.reference in r.id), None)
        if ref_rec:
            ref_len = len(str(ref_rec.seq).replace('-', ''))
        else:
            print(f'[警告] 未找到参考 {args.reference}, 用中位数长度')
            ref_len = sorted(seq_lens.values())[len(recs) // 2]
    else:
        ref_len = sorted(seq_lens.values())[len(recs) // 2]
    print(f'参考长度基准: {ref_len} bp')

    # 3. MAFFT 比对
    if args.skip_mafft:
        aln_path = args.input
        print('跳过 MAFFT (输入已比对)')
    else:
        aln_path = str(outdir / 'alignment.mafft.fasta')
        run_mafft(args.input, aln_path, args.threads)
        print(f'MAFFT 完成: {aln_path}')

    aln = list(SeqIO.parse(aln_path, 'fasta'))
    L = len(aln[0])

    # 4. 逐序列质量检查
    ref_seq = None
    if args.reference:
        ref_rec = next((r for r in aln if args.reference in r.id), None)
        if ref_rec:
            ref_seq = str(ref_rec.seq).upper()

    report = []
    keep, removed = [], []
    reasons = Counter()

    for r in aln:
        s = str(r.seq).upper()
        raw_len = len(s.replace('-', ''))
        n_ratio = s.count('N') / L if L else 1.0
        gap_ratio = s.count('-') / L if L else 1.0
        len_ratio = raw_len / ref_len if ref_len else 1.0

        # identity (与参考)
        identity = pairwise_identity(s, ref_seq) * 1.0 if ref_seq else None

        # 判定
        issues = []
        if len_ratio < args.min_length:
            issues.append(f'长度不足 ({len_ratio:.2f} < {args.min_length})')
        if gap_ratio > args.max_gap:
            issues.append(f'gap过多 ({gap_ratio:.2f} > {args.max_gap})')
        if n_ratio > args.max_n:
            issues.append(f'N过多 ({n_ratio:.2f} > {args.max_n})')
        if identity is not None:
            expected = (1 - args.max_identity_drop) * 100
            if identity < expected:
                issues.append(f'与参考identity过低 ({identity:.1f}% < {expected:.1f}%)')

        status = 'REMOVE' if issues else 'KEEP'
        reason = '; '.join(issues) if issues else ''
        report.append({
            'id': r.id, 'raw_length': raw_len, 'aln_length': L,
            'len_ratio': round(len_ratio, 4), 'gap_ratio': round(gap_ratio, 4),
            'n_ratio': round(n_ratio, 4),
            'identity_to_ref': round(identity, 2) if identity is not None else '',
            'status': status, 'reason': reason,
        })
        if issues:
            removed.append(r)
            for i in issues:
                reasons[i.split(' ')[0]] += 1
        else:
            keep.append(r)

    # 5. 输出报告
    report_path = outdir / 'align_qc_report.tsv'
    with open(report_path, 'w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=list(report[0].keys()), delimiter='\t')
        w.writeheader()
        w.writerows(report)

    # 6. 输出 clean / removed fasta
    keep_path = outdir / 'clean.fasta'
    rem_path = outdir / 'removed.fasta'
    SeqIO.write(keep, keep_path, 'fasta')
    SeqIO.write(removed, rem_path, 'fasta')

    # 7. 汇总打印
    print(f'\n=== 比对质量检查结果 ===')
    print(f'总序列: {len(aln)}, 保留 {len(keep)}, 去除 {len(removed)}')
    print(f'去除理由统计: {dict(reasons)}')
    print(f'\n比对长度: {L} bp (参考 {ref_len} bp)')

    if removed:
        print('\n被去除的序列:')
        for r in report:
            if r['status'] == 'REMOVE':
                print(f"  ✗ {r['id']}: 长度{r['raw_length']} gap={r['gap_ratio']:.0%} "
                      f"N={r['n_ratio']:.0%} → {r['reason']}")

    print(f'\n输出:')
    print(f'  报告: {report_path}')
    print(f'  保留: {keep_path} ({len(keep)} 条)')
    print(f'  去除: {rem_path} ({len(removed)} 条)')
    return 0


if __name__ == '__main__':
    sys.exit(main())
