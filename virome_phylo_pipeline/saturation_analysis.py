#!/usr/bin/env python3
"""
saturation_analysis.py — 替换饱和分析 (Substitution Saturation)

两种判定:
  C 值法  (YR-MPE 参考):  C = std(Ti/Tv) / std(p).  C 低 → 序列信息被饱和覆盖 (高饱和)
  Iss 法  (Xia 2018, 标准判据):  Iss vs Iss.c — Iss < Iss.c → 未饱和 (可用)

输入: 比对 FASTA
输出: saturation_report.csv + 散点图 (p vs Ti/Tv) + 直方图

用法:
  python -m saturation_analysis --fasta mafft.aln.fasta --outdir sat_out
"""

import argparse
import csv
import os
import sys
from pathlib import Path

import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq

# 转换/颠换定义
TRANSITIONS = {('A', 'G'), ('G', 'A'), ('T', 'C'), ('C', 'T')}
TRANSVERSIONS = {('A', 'T'), ('A', 'C'), ('G', 'T'), ('G', 'C'),
                 ('T', 'A'), ('T', 'G'), ('C', 'A'), ('C', 'G')}


def transition_transversion_counts(seq1, seq2, gap_treatment='pairwise'):
    """逐对 Ti/Tv/总位点计数"""
    ti = tv = total = 0
    min_len = min(len(seq1), len(seq2))
    for i in range(min_len):
        b1, b2 = seq1[i].upper(), seq2[i].upper()
        if b1 == '-' or b2 == '-':
            continue
        if b1 not in 'ATCG' or b2 not in 'ATCG':
            continue
        total += 1
        if b1 == b2:
            continue
        if (b1, b2) in TRANSITIONS:
            ti += 1
        elif (b1, b2) in TRANSVERSIONS:
            tv += 1
    return ti, tv, total


def complete_deletion(sequences):
    """保留所有序列均无 gap 且无简并碱基的位点。

    口径说明: 比 YR-MPE 更严（YR-MPE 只删含 '-'/'.' 的列，保留 N/简并碱基列），
    但与 MEGA 的 complete deletion 口径一致（含 N 的列也剔除）。
    有意保留未对齐 YR-MPE：简并碱基在 Ti/Tv 统计中语义不明，剔除更稳。
    """
    if not sequences:
        return sequences
    n = len(sequences)
    seqs = [str(s.seq).upper() for s in sequences]
    keep = []
    for col in range(min(len(x) for x in seqs)):
        bases = [x[col] for x in seqs]
        if '-' not in bases and all(b in 'ATCG' for b in bases):
            keep.append(col)
    out = []
    for s, orig in zip(seqs, sequences):
        new = ''.join(s[i] for i in keep)
        # 历史坑 ①: SeqRecord(seq, id) 第一个位置参数是 seq, 传 orig.id 会把
        #   .seq 变成 ID 字符串、.id 变成序列 → complete_deletion 结果全错。
        #   故必须用关键字传 id/name/description。
        # 历史坑 ② (2026-09-16 上游核对修复): 关键字传参解决了 ①, 但第一个位置参数
        #   传的 `new` 是 str, 而 SeqRecord 只接受 Seq/MutableSeq →
        #   `TypeError: seq argument should be a Seq or MutableSeq object`。
        #   该路径只在 --gap complete 下触发, 默认 pairwise 不走, 故长期未暴露。
        #   实测: `python -m saturation_analysis --fasta X.fa --gap complete` 必崩。
        out.append(orig.__class__(Seq(new), id=orig.id, name=orig.name,
                                  description=orig.description))
    return out


def calculate_saturation_metrics(sequences, gap_treatment='pairwise'):
    """计算饱和度指标 (C 值法, YR-MPE 参考)"""
    n = len(sequences)
    empty = {'p_values': [], 'ti_tv_ratios': [], 'c_value': 0.0,
             'std_ti_tv': 0.0, 'std_p': 0.0, 'n_pairs': 0}
    if n < 2:
        return empty
    if gap_treatment == 'complete':
        sequences = complete_deletion(sequences)
        if len(sequences) < 2:
            return empty

    p_values, ti_tv_ratios = [], []
    for i in range(n):
        for j in range(i + 1, n):
            ti, tv, total = transition_transversion_counts(
                str(sequences[i].seq), str(sequences[j].seq), 'pairwise')
            if total > 0:
                p = (ti + tv) / total
                ratio = (ti / tv) if tv > 0 else 10.0
                p_values.append(p)
                ti_tv_ratios.append(ratio)

    if len(p_values) >= 2:
        std_p = np.std(p_values, ddof=1)
        std_ti_tv = np.std(ti_tv_ratios, ddof=1)
        c_value = std_ti_tv / std_p if std_p > 0 else 0.0
    else:
        std_p = std_ti_tv = c_value = 0.0
    return {'p_values': p_values, 'ti_tv_ratios': ti_tv_ratios,
            'c_value': c_value, 'std_ti_tv': std_ti_tv, 'std_p': std_p,
            'n_pairs': len(p_values)}


def xia_iss(sequences, num_replicates=100, seed=42):
    """简化版饱和度指数 (非 Xia 2003 原版, 命名仅示渊源)

    本实现: Iss = p_obs / p_sat
      p_obs = 平均 pairwise p-distance; p_sat = 1 - Σf_i² (碱基频率随机替换期望差异)
      Iss→1 高饱和; Iss 低未饱和。方向与 Xia 一致, 但口径不同:
      Xia 2003 原版 Iss = Ho/He (位点平均熵 / 随机序列期望熵), 且需对称
      重采样估 Iss.c 临界值做检验 (DAMBE 口径)。本函数不做重采样、
      不产 Iss.c/p_value (返回 None), 饱和判定用经验阈值 0.8。
      需发表级 Xia 检验请用 DAMBE 或按原文补重采样。
    """
    seqs = [str(s.seq).upper() for s in sequences]
    n = len(seqs)
    if n < 2:
        return {'iss': None, 'iss_c': None, 'saturated': None, 'p_value': None}

    min_len = min(len(x) for x in seqs)

    # 观察差异比例 (平均 pairwise p)
    diffs = []
    for i in range(n):
        for j in range(i + 1, n):
            comp = diff = 0
            for k in range(min_len):
                b1, b2 = seqs[i][k], seqs[j][k]
                if b1 == '-' or b2 == '-':
                    continue
                comp += 1
                if b1 != b2:
                    diff += 1
            if comp > 0:
                diffs.append(diff / comp)
    p_obs = float(np.mean(diffs)) if diffs else 0.0

    # 碱基频率 → 饱和期望差异 p_sat = 1 - Σ f_i²
    from collections import Counter
    cnt = Counter(b for s in seqs for b in s if b in 'ACGT')
    tot = sum(cnt.values())
    p_sat = (1 - sum((c / tot) ** 2 for c in cnt.values())) if tot else 0.75
    if p_sat <= 0:
        p_sat = 0.75

    iss = p_obs / p_sat if p_sat > 0 else 0.0
    iss = max(0.0, min(1.0, iss))

    # 经验阈值: Iss > 0.8 → 饱和风险 (接近随机替换水平)
    saturated = iss > 0.8
    return {'iss': iss, 'iss_c': None, 'saturated': bool(saturated), 'p_value': None}


def _norm_cdf(x):
    """标准正态 CDF (无 scipy 依赖)"""
    import math
    return 0.5 * (1 + math.erf(x / np.sqrt(2)))


def plot_saturation(p_values, ti_tv_ratios, out_path):
    """p vs Ti/Tv 散点图 (饱和曲线)"""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.scatter(p_values, ti_tv_ratios, s=12, alpha=0.55, color='#0F4D92', edgecolors='none')
    ax.set_xlabel('p-distance (observed differences)')
    ax.set_ylabel('Ti/Tv ratio')
    ax.set_title('Substitution Saturation (p vs Ti/Tv)')
    if len(p_values) > 1:
        # 平滑趋势线
        order = np.argsort(p_values)
        xs = np.array(p_values)[order]
        ys = np.array(ti_tv_ratios)[order]
        try:
            z = np.polyfit(xs, ys, 3)
            xx = np.linspace(xs.min(), xs.max(), 100)
            ax.plot(xx, np.polyval(z, xx), 'r--', linewidth=1.2, label='trend')
            ax.legend()
        except Exception:
            pass
    fig.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches='tight')
    plt.close(fig)
    return out_path


def run_saturation(fasta_path, outdir='saturation_out', gap_treatment='pairwise',
                   num_replicates=1000, log=None):
    """饱和分析全流程"""
    def emit(msg):
        if log:
            log(msg)
    os.makedirs(outdir, exist_ok=True)
    result = {'success': False}

    if not os.path.exists(fasta_path):
        result['error'] = f'比对文件不存在: {fasta_path}'
        return result

    sequences = list(SeqIO.parse(fasta_path, 'fasta'))
    if len(sequences) < 2:
        result['error'] = '序列数 < 2'
        return result
    emit(f'[saturation] {len(sequences)} 条序列')

    # C 值法
    m = calculate_saturation_metrics(sequences, gap_treatment)
    # Xia Iss
    xia = xia_iss(sequences, num_replicates)

    # 综合判定: 两个判据
    # C 值: 低 (C<5) → 高饱和风险; 高 (C>20) → 无饱和
    # Iss:  Iss = p_obs/p_sat, Iss > 0.8 → 接近饱和 (简化版, 不再重采样估 Iss.c)
    mean_p = float(np.mean(m['p_values'])) if m['p_values'] else 0.0
    c_saturated = m['c_value'] < 5.0
    iss_saturated = bool(xia.get('saturated', False))
    # 综合: C 值法为主 (YR-MPE 参考), Iss 为辅 (Iss>0.8 或 C<5 都提示饱和)
    saturated = c_saturated or iss_saturated
    if c_saturated and not iss_saturated:
        verdict = '饱和 (C 值低, Iss 未达阈值; 建议看散点图)'
    elif not c_saturated and iss_saturated:
        verdict = '偏高 (Iss>0.8; C 值未触发)'
    elif saturated:
        verdict = '饱和'
    else:
        verdict = '未饱和'
    emit(f"[saturation] C={m['c_value']:.4f}, Iss={xia['iss']:.4f}, p={mean_p:.4f} → {verdict}")

    # 报告 CSV
    report = os.path.join(outdir, 'saturation_report.csv')
    with open(report, 'w', newline='', encoding='utf-8') as f:
        w = csv.writer(f)
        w.writerow(['metric', 'value', 'criterion'])
        w.writerow(['n_sequences', len(sequences), ''])
        w.writerow(['n_pairs', m['n_pairs'], ''])
        w.writerow(['C_value', f"{m['c_value']:.4f}", 'C < 5 → 高饱和; C > 20 → 未饱和'])
        w.writerow(['Iss', f"{xia['iss']:.4f}", 'Iss = p_obs/p_sat; >0.8 接近饱和'])
        w.writerow(['saturated', saturated, '综合判定 (C 值法为主 + Iss 辅)'])
        w.writerow(['mean_p_distance', f"{mean_p:.4f}", ''])
        w.writerow(['verdict', verdict, ''])
    emit(f'[saturation] 报告 → {report}')

    # 图
    plot_path = None
    if m['p_values']:
        plot_path = plot_saturation(m['p_values'], m['ti_tv_ratios'],
                                    os.path.join(outdir, 'saturation_plot.png'))
        emit(f'[saturation] 图 → {plot_path}')

    result.update({'success': True, 'report': report, 'plot': plot_path,
                   'c_value': m['c_value'], 'iss': xia['iss'], 'iss_c': xia.get('iss_c'),
                   'saturated': saturated, 'n_sequences': len(sequences),
                   'mean_p': mean_p, 'verdict': verdict})
    return result


def main():
    ap = argparse.ArgumentParser(description='替换饱和分析 (C 值法 + Xia Iss)')
    ap.add_argument('--fasta', required=True, help='比对 FASTA')
    ap.add_argument('--outdir', default='saturation_out')
    ap.add_argument('--gap', default='pairwise', choices=['pairwise', 'complete'])
    ap.add_argument('--replicates', type=int, default=1000, help='Iss 重采样次数')
    args = ap.parse_args()
    r = run_saturation(args.fasta, args.outdir, args.gap, args.replicates,
                       log=lambda m: print(m))
    if r.get('success'):
        # Iss.c 简化版恒为 None (xia_iss 不产临界值), None:.4f 会崩 → 占位符
        iss_c = '—' if r['iss_c'] is None else f"{r['iss_c']:.4f}"
        print(f"✓ 饱和分析完成: C={r['c_value']:.4f}, Iss={r['iss']:.4f}, "
              f"Iss.c={iss_c}, 饱和={'是' if r['saturated'] else '否'}")
        return 0
    print(f"✗ 失败: {r.get('error')}")
    return 1


if __name__ == '__main__':
    sys.exit(main())
