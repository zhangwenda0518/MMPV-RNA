#!/usr/bin/env python3
"""S4 判别 -> 可执行动作 (filter 表) [判定逻辑最后变更: 模块 v5, 自本文件 v1 起未变].

把 s3_verdict.py 的逐候选 verdict 翻译成下游能直接按列过滤的动作。这一步原本写在
virome_pipeline.py 的阶段后处理里; EVE 判别改为独立后运行脚本后, 表也得跟着搬过来,
否则输出的 verdict 没人翻译成 action, 下游过滤还是得手工做。

verdict -> action 的映射是一对多的, 共 14 种 verdict:

  host_contamination_likely    REMOVE_host_contamination  自身寄主物种高同源且上游无病毒信号
  host_homology_cross_species  REVIEW                     跨物种旁证口径不作自身寄主证据
  host_conflict_review         REVIEW                     自身寄主高同源与上游病毒注释矛盾
  EVE_STRONG_provirus          MOVE_EVE                   真位点齐全 + 密码子退化 = 整合前病毒
  EVE_LTR_TE                   MOVE_EVE                   转座子域/转座子科
  ancient_EVE                  MOVE_EVE                   纯 pol 区分布式退化
  virus_candidate              KEEP_virus                 单位点齐全 + 密码子完好 = 感染中的真病毒
  *review* / review            REVIEW                     证据不足, 不得据此下 EVE 结论

输入: s3 产出的 eve_distinguish_verdict.tsv
输出: dna_vs_eve_filter.tsv (在 verdict 表后面追加 action / action_reason 两列,
 并丢掉只对调试有用的中间列, 让下游 review 时一眼能看到判据)。

v6 注: 后两条 host_* 都是 REVIEW 而不是 REMOVE。`host_contamination_likely` 是唯一
能删数据的判词, 旧版把跨物种命中和上游有病毒科注释的候选也塞进它 —— 1027 条真实候选
上 48 次调用全部打错 (全是带病毒科注释的真病毒)。删数据需要最强的证据, 拿到弱证据时
只能升级人工复核。
"""
import argparse
import csv
import sys
from collections import Counter

ACTIONS = {
    'host_contamination_likely':  ('REMOVE_host_contamination',
                                   '自身寄主物种聚合 >=90%/80%, 无 MP/CP/AP, 上游无病毒信号'),
    # v6: 这两个曾经都写成 host_contamination_likely (会删数据)。现在分开, 都只 REVIEW
    'host_homology_cross_species': ('REVIEW',
                                    '跨物种旁证口径 (样本无自身物种映射), 不作寄主判定依据'),
    'host_conflict_review':       ('REVIEW',
                                   '自身寄主高同源与上游病毒注释矛盾, 需人工定案'),
    'EVE_LTR_TE':                 ('MOVE_EVE', 'LTR 逆转录转座子域/转座子科'),
    'EVE_STRONG_provirus':        ('MOVE_EVE', '真位点 4/5 组件齐全 + 密码子分布式退化 = 整合前病毒'),
    'ancient_EVE':                ('MOVE_EVE', '纯 pol 区分布式退化'),
    'virus_candidate':            ('KEEP_virus', '单位点齐全 + 密码子完好 = 感染中的真病毒'),
    'virus_fragment_review':      ('REVIEW', '架构不全 + 编码完好; RNA-seq 里真 DNA 病毒本就拼不全'),
    'assembly_breakpoint_review': ('REVIEW', '终止集中在单个断裂点, rescue 延伸后复判'),
    'compositional_noise_review': ('REVIEW', '病毒读框不比其他读框更脏, 无退化证据'),
    'EVE_suspect':                ('REVIEW', '架构不全 + 分布式退化, contig 级无法定案'),
    'structure_intact_review':    ('REVIEW', '结构完整/编码保留, Cauli 面板不适用'),
    'structure_decay_review':     ('REVIEW', '分布式退化但无 Cauli 面板证据, 不下 EVE 结论'),
    'review':                     ('REVIEW', '证据不足'),
}

# 逐列原样搬到 filter 表的 verdict 列; 顺序即输出列顺序.
# v5: 加 host_scope —— host_contamination_likely 里"自身物种"和"跨物种旁证"是两种
# 强度的证据, 下游要能分开看 (旧版只有结论没有口径).
# v6.3: 加 rg_call/rg_scope —— 参考基因组通道的证据口径; 旧表没有这两列时按 '' 搬,
# 不影响只跑 s1-s4 的旧目录.
CARRY = ['contig_id', 'tax_family', 'category', 'locus_ncomp', 'locus_completeness',
         'decay_class', 'stop_enrichment', 'host_wpid', 'host_cov', 'host_scope',
         'rg_call', 'rg_scope', 'verdict']
OUT_HDR = CARRY + ['action', 'action_reason']


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('verdict_tsv', help='s3_verdict.py 的产出')
    ap.add_argument('-o', '--out', default='dna_vs_eve_filter.tsv', help='输出路径')
    a = ap.parse_args()

    n_act = Counter()
    n_row = 0
    unknown = Counter()
    with open(a.verdict_tsv, newline='', encoding='utf-8') as fin, \
            open(a.out, 'w', newline='', encoding='utf-8') as fout:
        w = csv.writer(fout, delimiter='\t', lineterminator='\n')
        w.writerow(OUT_HDR)
        for row in csv.DictReader(fin, delimiter='\t'):
            v = row.get('verdict', '')
            act, why = ACTIONS.get(v, (None, ''))
            if act is None:
                # s3 的判定表之外的新 verdict: 不猜, 记未映射并让人工看
                unknown[v] += 1
                act, why = 'REVIEW', '未映射判别 (verdict 值不在 ACTIONS 表内)'
            n_act[act] += 1
            n_row += 1
            w.writerow([row.get(c, '') for c in CARRY] + [act, why])

    print(f"[done] s4 -> {a.out}")
    print(f"候选总数 {n_row}")
    for act, n in n_act.most_common():
        print(f"  {act:28s} {n}")
    if unknown:
        print(f"[warn] {sum(unknown.values())} 条 verdict 未映射, 已按 REVIEW 处理: "
              f"{dict(unknown)}", file=sys.stderr)
        return 2
    return 0


if __name__ == '__main__':
    sys.exit(main())
