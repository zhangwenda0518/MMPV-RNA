#!/usr/bin/env python3
"""S4 判别 -> 可执行动作 (filter 表) v1.

把 s3_verdict.py 的逐候选 verdict 翻译成下游能直接按列过滤的动作。这一步原本写在
virome_pipeline.py 的阶段后处理里; EVE 判别改为独立后运行脚本后, 表也得跟着搬过来,
否则输出的 verdict 没人翻译成 action, 下游过滤还是得手工做。

verdict -> action 的映射是一对多的, 共 12 种 verdict:

  host_contamination_likely    REMOVE_host_contamination  候选就是自身寄主的序列
  EVE_STRONG_provirus          MOVE_EVE                  单位点齐全 + 密码子退化 = 整合前病毒
  EVE_LTR_TE                   MOVE_EVE                  转座子域/转座子科
  ancient_EVE                  MOVE_EVE                  纯 pol 区分布式退化
  virus_candidate              KEEP_virus                单位点齐全 + 密码子完好 = 感染中的真病毒
  *review* / review            REVIEW                    证据不足, 不得据此下 EVE 结论

输入: s3 产出的 eve_distinguish_verdict.tsv
输出: dna_vs_eve_filter.tsv (在 verdict 表后面追加 action / action_reason 两列,
 并丢掉只对调试有用的中间列, 让下游 review 时一眼能看到判据)。
"""
import argparse
import csv
import sys
from collections import Counter

ACTIONS = {
    'host_contamination_likely':  ('REMOVE_host_contamination',
                                   '自身寄主物种聚合 >=90%/80% 且无 MP/CP/AP 病毒结构基因'),
    'EVE_LTR_TE':                 ('MOVE_EVE', 'LTR 逆转录转座子域/转座子科'),
    'EVE_STRONG_provirus':        ('MOVE_EVE', '单位点 4/5 组件齐全 + 密码子分布式退化 = 整合前病毒'),
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

# 逐列原样搬到 filter 表的 verdict 列; 顺序即输出列顺序
CARRY = ['contig_id', 'tax_family', 'category', 'locus_ncomp', 'locus_completeness',
         'decay_class', 'stop_enrichment', 'host_wpid', 'host_cov', 'verdict']
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
