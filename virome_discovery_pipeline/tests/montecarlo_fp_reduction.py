#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
montecarlo_fp_reduction.py — graph 模块假阳性削减量的蒙特卡洛估计
====================================================================
问题: 组装图标签传播模块能减少多少假阳性?

真实数字需一次服务器实跑 (审计表直接可数, 见文末说明); 本脚本在其之前给出
参数校准的模拟区间。参数校准依据:
  * 样本规模: 植物宏转录组去宿主后 MEGAHIT/rnaSPAdes ≥1000bp contigs 数百~千级
  * 社区结构: 宿主链 (大, 强内部连接) / 细菌链 / 病毒小基因组链 (2-6 contigs)
    / 高拷贝质粒 / 无连接单 contig; 跨源噪声边多数 < min_link_pairs 被剪,
    少数强噪声 (重复序列) 与病毒-宿主嵌合边 (整合/LTR 邻接) 存活
  * GC/深度特征: 宿主 0.38/10x, 细菌 0.60/12x, 病毒 0.42/低-中覆盖, 质粒 0.5/高覆盖
  * 投票: 9 工具对真病毒敏感度 ~N(0.62,0.18); 对非病毒 per-tool FP 0.5%-2%
    (候选池定义 = ensemble 并集, 即 ≥1 票 —— 与 virus_identification 一致)

关键分解 (诚实边界):
  FP 候选池按票数拆成三层 —
    ≥2 票: 成为 viral_seed (不可变) → 本模块结构性不削减 (由下游 CDD/rescue 负责)
    =1 票: unknown → 图传播裁决 → 落在强非病毒社区的被削减 (主要收割区)
    (filter 层另有 0-evidence 保留者, 即 CDD/UniProt 双零命中放行 — 同样归图裁决)
  输出: FP 削减率 / TP 误伤 / COBRA 池缩减, 均为 mean±sd

用法: python virome_discovery_pipeline/tests/montecarlo_fp_reduction.py [--reps 15]
"""

import argparse
import random
import sys
from pathlib import Path

import networkx as nx
from networkx.algorithms.community import louvain_communities

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))
sys.path.insert(0, str(REPO_ROOT / "virome_discovery_pipeline"))

import graph_propagation as gp  # noqa: E402

CFG = {'min_link_pairs': 3, 'gc_diff': 0.20, 'depth_ratio': 20.0,
       'resolution': 1.0, 'seed_min_votes': 2, 'nonviral_min_len': 1000,
       'majority_ratio': 10.0}
N_TOOLS = 9
N_CONTIGS = 600


def clip(x, lo, hi):
    if hi is not None:
        x = min(hi, x)
    return max(lo, x)


def build_sample(rng, fp_rate, viral_frac):
    """构造一个模拟样本 → (contigs, links)。truth: viral/nonviral。"""
    contigs, links = {}, []
    n_viral_comm = max(2, round(viral_frac * N_CONTIGS / 4))
    v_quota = max(2, round(viral_frac * N_CONTIGS))
    quota = [('host', round(N_CONTIGS * 0.55)), ('bact', round(N_CONTIGS * 0.20)),
             ('viral', v_quota), ('plasmid', round(N_CONTIGS * 0.05))]
    used = sum(q for _, q in quota)
    cid = 0

    def add_comm(kind, members, gc_mu, gc_sd, dep_mu, dep_sd, chain_w_mu):
        nonlocal cid
        ids = []
        for _ in range(members):
            contigs[f"c{cid}"] = {
                'truth': 'viral' if kind == 'viral' else 'nonviral',
                'len': rng.randrange(1000, 6000),
                'gc': clip(rng.gauss(gc_mu, gc_sd), 0.25, 0.72),
                'depth': clip(rng.gauss(dep_mu, dep_sd), 0.5, None)}
            ids.append(f"c{cid}")
            cid += 1
        rng.shuffle(ids)
        for a, b in zip(ids, ids[1:]):     # 链拓扑
            links.append(((a, b), clip(rng.gauss(chain_w_mu, chain_w_mu / 3), 5, None)))
        return ids

    for kind, q in quota:
        if kind == 'host':
            half = q // 2
            add_comm('host', half, .38, .025, 10, 3, 30)
            add_comm('host', q - half, .40, .025, 10, 3, 30)
        elif kind == 'bact':
            add_comm('bact', q, .60, .04, 12, 4, 25)
        elif kind == 'plasmid':
            add_comm('plasmid', q, .50, .03, 200, 50, 28)
        else:
            per = max(2, q // n_viral_comm)
            for _ in range(n_viral_comm):
                take = min(per, q)
                if take <= 0:
                    break
                add_comm('viral', take, .42, .04, 8, 5, 26)
                q -= take
    while cid < N_CONTIGS:                 # 剩余 = 无连接单 contig (宿主/细菌起源)
        contigs[f"c{cid}"] = {'truth': 'nonviral', 'len': rng.randrange(1000, 6000),
                              'gc': clip(rng.uniform(0.30, 0.68), 0.25, 0.72),
                              'depth': clip(rng.gauss(10, 4), 0.5, None)}
        cid += 1

    # 跨源噪声: 弱 (p=0.03, 权重 1-2 → 被剪) 与 强 (p=0.005, 权重 5-9 → 存活)
    ids = list(contigs)
    for a in ids:
        if rng.random() < 0.03:
            b = rng.choice(ids)
            if b != a:
                links.append(((a, b), rng.uniform(1, 2)))
        if rng.random() < 0.005:
            b = rng.choice(ids)
            if b != a:
                links.append(((a, b), rng.uniform(5, 9)))
    # 病毒-宿主嵌合边 (整合/LTR 邻接): 每病毒 contig p=0.01, 权重 4-8
    viral_ids = [c for c, d in contigs.items() if d['truth'] == 'viral']
    for a in viral_ids:
        if rng.random() < 0.01:
            b = rng.choice([c for c in ids if contigs[c]['truth'] == 'nonviral'])
            links.append(((a, b), rng.uniform(4, 8)))
    return contigs, links


def votes_for(rng, contigs, fp_rate):
    votes = {}
    sens = clip(rng.gauss(0.62, 0.18), 0.25, 0.95)
    for c, d in contigs.items():
        if d['truth'] == 'viral':
            votes[c] = sum(1 for _ in range(N_TOOLS) if rng.random() < sens)
        else:
            p1 = clip(rng.gauss(fp_rate, fp_rate / 2), 0.0005, None)
            votes[c] = sum(1 for _ in range(N_TOOLS) if rng.random() < p1)
    return votes


def run_one(rng, fp_rate, viral_frac):
    contigs, links = build_sample(rng, fp_rate, viral_frac)
    votes = votes_for(rng, contigs, fp_rate)
    lengths = {c: d['len'] for c, d in contigs.items()}
    gc = {c: d['gc'] for c, d in contigs.items()}
    depth = {c: d['depth'] for c, d in contigs.items()}
    seeds = gp.assign_seeds(contigs, lengths, votes, set(), **{
        k: CFG[k] for k in ('seed_min_votes', 'nonviral_min_len')})
    lc = gp.count_cross_links(
        [f"r{i}\t99\t{a}\t1\t60\t150M\t{b}\t2\t0\tA\t*\n"
         for (a, b), w in links for i in range(int(w))], require_proper=True)
    edges = gp.prune_edges(lc, gc, depth, min_link_pairs=CFG['min_link_pairs'],
                           gc_diff=CFG['gc_diff'], depth_ratio=CFG['depth_ratio'],
                           has_depth=True)
    verdicts, _, _ = gp.propagate(edges, contigs, seeds,
                                  resolution=CFG['resolution'],
                                  majority_ratio=CFG['majority_ratio'])
    candidates = {c for c, v in votes.items() if v >= 1}
    fp_pool = {c for c in candidates if contigs[c]['truth'] == 'nonviral'}
    demoted = {c for c, v in verdicts.items() if v == gp.V_GRAPH_NONVIRAL}
    fp_cut = fp_pool & demoted
    tp_loss = {c for c in candidates if contigs[c]['truth'] == 'viral'} & demoted
    protected = {c for c in fp_pool if votes[c] >= CFG['seed_min_votes']}
    weak = fp_pool - protected
    return {'cand': len(candidates), 'fp': len(fp_pool), 'fp_cut': len(fp_cut),
            'protected': len(protected), 'weak': len(weak),
            'weak_cut': len(weak & demoted), 'tp_loss': len(tp_loss)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--reps', type=int, default=15)
    args = ap.parse_args()

    print(f"蒙特卡洛: N_CONTIGS={N_CONTIGS}, 9 工具, 候选池=≥1票, 种子≥2票, "
          f"10x 传播规则, {args.reps} 重复/组合\n")
    header = (f"{'tool FP率':>9} {'病毒占比':>8} | {'候选池':>7} {'FP池':>6} "
              f"{'≥2票保护':>8} {'1票可削':>7} | {'FP削减':>12} {'削减率(全FP)':>11} "
              f"{'削减率(1票)':>11} {'TP误伤':>6}")
    print(header)
    print('-' * len(header))
    for fp_rate in (0.005, 0.01, 0.02):
        for viral_frac in (0.02, 0.05):
            rows = [run_one(random.Random(1000 + i), fp_rate, viral_frac)
                    for i in range(args.reps)]
            mean = lambda k: sum(r[k] for r in rows) / len(rows)
            sd = lambda k: (sum((r[k] - mean(k)) ** 2 for r in rows)
                            / (len(rows) - 1)) ** 0.5 if len(rows) > 1 else 0
            fr_all = [r['fp_cut'] / r['fp'] for r in rows if r['fp']]
            fr_weak = [r['weak_cut'] / r['weak'] for r in rows if r['weak']]
            m = lambda xs: sum(xs) / len(xs) if xs else 0.0
            print(f"{fp_rate:>9.3%} {viral_frac:>8.0%} | "
                  f"{mean('cand'):>7.1f} {mean('fp'):>6.1f} "
                  f"{mean('protected'):>8.1f} {mean('weak'):>7.1f} | "
                  f"{mean('fp_cut'):>5.1f}±{sd('fp_cut'):<5.1f} "
                  f"{m(fr_all):>10.1%} {m(fr_weak):>11.1%} "
                  f"{mean('tp_loss'):>6.2f}")

    print("\n读法:")
    print("  * FP削减率(全FP) = 模块对候选池假阳性的总削减比例; 上限 = 1 - ≥2票占比")
    print("    (≥2票假阳性成为 viral_seed 不可变 — 结构上不由本模块削减, 由下游 CDD/rescue 兜底)")
    print("  * FP削减率(1票) = 可削减区内的实际收割率 (1票弱阳性落在强非病毒社区的比例)")
    print("  * TP误伤 = 真病毒候选被降级的数目 (期望≈0; 病毒-宿主嵌合边是主要风险源)")
    print("真实数字 (服务器实跑后, 直接可数):")
    print("  awk -F'\\t' '$12==\"graph_nonviral\"' d_graph/<S>/<S>.graph_verdict.tsv | wc -l")
    print("  grep -c demoted d_filter/<S>/graph/*_graph_calls.tsv  # filter 层实削数")


if __name__ == '__main__':
    main()
