#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
graph_propagation.py — 组装图标签传播 (③ discovery 的 graph 阶段)
====================================================================
定位: 插在 identification (02a_Identification) 与 filter (02b_Filter/03_Filter)
之间, 给组装 contig 的 viral/non-viral 归属做一次图拓扑层面的重估。

思想来源: GMW (Chen et al. 2026, Sci China Life Sci 69:1910-1917) — 在组装图上
做 guilt-by-association 去污染。本实现是它的配对读段改编 (不使用 GMW 代码):
  * 组装图 GFA → 配对读段链接图 (assembler 无关, MEGAHIT/rnaSPAdes/Penguin 通吃)
  * Kraken2 单源种子 → 本管线 10 工具 ensemble 阳性清单做种子 (质量更高)
  * Louvain 社区传播 + GMW >=10 倍多数比规则保持一致
  * 断边规则沿用 GMW: GC 差 >20% / 深度比 >20:1
  * 链接定义: 双端都比对上的主比对 — 跨 contig pair 在 bwa 语义下永远 non-proper
    (0x2 仅同 contig 置位), 故默认不要求 proper flag

输入 (逐样本):
  d_asm/{S}/{S}_{tool}.contig.fasta        组装 contig (节点全集, 不只候选)
  d_hostdep/{S}_1.fa.gz + {S}_2.fa.gz      去宿主 clean reads (PE; SE 优雅退化)
  d_ident/{N}/{N}_virus.{tool}.result.id   ensemble 阳性清单 (N = 去 _clean 样本名)

输出:
  d_graph/{S}/{S}.graph_verdict.tsv        逐 contig 判定 + 审计列 (契约见 VERDICT_COLUMNS)
  d_graph/{S}/{S}.graph_summary.tsv        社区级汇总
  d_graph/{S}/{S}.graph.done               样本级断点标记
  d_graph/candidates/{N}/{N}_virus.all.candidate.fasta
                                           派生候选 = ident 候选 - graph_nonviral
                                           (子目录命名逐镜像 d_ident, 供 COBRA raw
                                           模式按现成模式名发现, cobra_pipeline 零改动)

判定语义:
  viral_seed     被 >= seed_min_votes 个工具判为病毒 (种子, 不可变)
  nonviral_seed  len >= nonviral_min_len 且 0 票且非 viroid 阳性 (种子, 不可变)
  graph_viral    unknown 落在 V >= majority_ratio*N 的社区
  graph_nonviral unknown 落在 N >= majority_ratio*V 的社区
  conflicting    社区两种种子都有但都未达多数比 (保持 unknown)
  unknown        其余 (无链接/无种子/短 contig 未被工具检查过)

消费点:
  * filter 阶段末端 pass 读 {S}.graph_verdict.tsv (rescue/demote, 见 filter_virus.py)
  * 编排器 _cobra_candidate_source 在 apply_to_cobra 开启时把 COBRA raw 候选源
    切到 d_graph/candidates

规模兜底: --max-reads 每端取前 N 条 (确定性前缀截取, 可复现); --jobs 跨样本并行。
"""

import argparse
import gzip
import itertools
import os
import shutil
import subprocess
import sys
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import networkx as nx
from networkx.algorithms.community import louvain_communities

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))
from mmpv_common.io_layout import layout_dirs  # noqa: E402

_D = layout_dirs(os.environ.get("MMPV_IO_LAYOUT"))

# ── 判定值域与 TSV 契约 (filter 消费点按列名读取, 改动需同步 filter_virus.py) ──
V_VIRAL_SEED, V_NONVIRAL_SEED = "viral_seed", "nonviral_seed"
V_GRAPH_VIRAL, V_GRAPH_NONVIRAL = "graph_viral", "graph_nonviral"
V_UNKNOWN, V_CONFLICTING = "unknown", "conflicting"

VERDICT_COLUMNS = ("contig", "source_tool", "length", "gc", "depth", "votes",
                   "seed_label", "community_id", "community_size",
                   "community_viral", "community_nonviral", "verdict", "reason")

SUMMARY_COLUMNS = ("source_tool", "community_id", "size", "viral_seed",
                   "nonviral_seed", "graph_viral", "graph_nonviral",
                   "conflicting", "unknown", "rule")

READ_PATTERNS = ("_1.", "_2.", ".R1.", ".R2.", "_R1.", "_R2.")


# ═══════════════════════════ 纯函数 (单测覆盖) ═══════════════════════════

def read_fasta_lengths_gc(path):
    """流式读取 FASTA → {contig_id: (length, gc_fraction)}。id 取 header 首 token。"""
    out, cur, length, gc = {}, None, 0, 0
    with _open_maybe_gz(path) as fh:
        for line in fh:
            if line.startswith('>'):
                if cur is not None:
                    out[cur] = (length, (gc / length) if length else 0.0)
                cur = line[1:].split()[0] if len(line) > 1 else ''
                length = gc = 0
            elif cur is not None:
                seq = line.strip().upper()
                length += len(seq)
                gc += sum(1 for ch in seq if ch in 'GC')
    if cur is not None:
        out[cur] = (length, (gc / length) if length else 0.0)
    return out


def _open_maybe_gz(path):
    return gzip.open(path, 'rt') if str(path).endswith('.gz') else open(path, 'r')


def load_seed_votes(ident_dir, sample_variants):
    """从 d_ident/{N}/ 读 ensemble 阳性清单 → (votes Counter, viroid_ids set, n, 命中目录)。

    计票文件 = {N}_virus.*.result.id (9 工具阳性并集的来源, 不含分数);
    viroid 阳性 ({N}_viroids.result.id) 不计票, 只作排除集 (类病毒不进病毒并集)。
    sample_variants 依次尝试 (处理 asm 样本名带 _clean 而 ident 目录不带 的差异)。
    """
    ident_dir = Path(ident_dir)
    for n in sample_variants:
        sdir = ident_dir / n
        if not sdir.is_dir():
            continue
        votes, viroid = Counter(), set()
        result_files = sorted(sdir.glob(f"{n}_virus.*.result.id"))
        # 排除 ensemble 并集文件 (all.result.id 是 {tool}.result.id 的派生并集,
        # 不是独立工具 — 计入会使每条候选虚增 1 票, 1 票弱阳性被错误保护成种子;
        # 246 mock benchmark 实测抓出)
        result_files = [f for f in result_files
                        if not os.path.basename(f).startswith(f"{n}_virus.all.")]
        for rf in result_files:
            with open(rf, encoding='utf-8', errors='replace') as fh:
                for line in fh:
                    cid = line.strip()
                    if cid and not cid.startswith('#'):
                        votes[cid] += 1
        viroid_files = sorted(sdir.glob(f"{n}_viroids.result.id"))
        for vf in viroid_files:
            with open(vf, encoding='utf-8', errors='replace') as fh:
                for line in fh:
                    cid = line.strip()
                    if cid and not cid.startswith('#'):
                        viroid.add(cid)
        if result_files or viroid_files:
            return votes, viroid, n, sdir
    return Counter(), set(), None, None


def assign_seeds(contig_ids, lengths, votes, viroid_ids,
                 seed_min_votes=2, nonviral_min_len=1000):
    """contig → 种子标签 (viral_seed / nonviral_seed / unknown)。纯函数。

    票数 < seed_min_votes 的弱阳性 (单工具命中, 假阳性主要来源) 归 unknown,
    交给图传播裁决 — 若钉成 nonviral_seed (不可变), filter 降级与 COBRA 拦截
    两个消费点将永远无法触发。
    """
    seeds = {}
    for cid in contig_ids:
        v = votes.get(cid, 0)
        if v >= seed_min_votes:
            seeds[cid] = V_VIRAL_SEED
        elif v >= 1:
            seeds[cid] = V_UNKNOWN          # 弱阳性: 证据不足, 图裁决
        elif cid in viroid_ids:
            seeds[cid] = V_UNKNOWN          # 类病毒阳性: 不进病毒并集, 也不当非病毒种子
        elif lengths.get(cid, 0) >= nonviral_min_len:
            seeds[cid] = V_NONVIRAL_SEED    # 被工具检查过且全体阴性
        else:
            seeds[cid] = V_UNKNOWN          # 短 contig 未被工具检查过, 不可判阴性
    return seeds


def count_cross_links(sam_lines, require_proper=True):
    """解析 samtools view 文本流 → {frozenset-free (a,b) 归一对: 链接数}。

    过滤: 只收 paired 主比对且双端都比对上 (排除 unmapped/secondary/supplementary);
    require_proper 时再要求 proper pair (0x2)。RNEXT '=' 与 RNAME 相同 → 同 contig,
    不建跨边; RNEXT '*' → mate 未比对, 跳过。
    """
    need = 0x1 | 0x40                    # paired + first-in-pair (每对只从 read1 记一次)
    if require_proper:
        need |= 0x2
    drop = 0x4 | 0x8 | 0x100 | 0x800     # self/mate unmapped, secondary, supplementary
    links = Counter()
    for line in sam_lines:
        if not line or line.startswith('@'):
            continue
        f = line.split('\t', 9)
        if len(f) < 9:
            continue
        flag = int(f[1])
        if (flag & need) != need or (flag & drop):
            continue
        rname, rnext = f[2], f[6]
        if rname == '*' or rnext in ('*', '='):
            continue
        key = (rname, rnext) if rname <= rnext else (rnext, rname)
        links[key] += 1
    return links


def prune_edges(links, gc, depth, min_link_pairs=3,
                gc_diff=0.20, depth_ratio=20.0, has_depth=False):
    """GMW 断边: 权重/GC 差/深度比三规则 → {边: 权重}。纯函数。

    depth 表整体缺失时 (has_depth=False) 跳过深度规则 (调用方负责告警);
    单边一端 depth<=0 视为比值无穷 → 断 (0 深度 contig 实际不会有边, 双保险)。
    """
    kept = {}
    for (a, b), w in links.items():
        if w < min_link_pairs:
            continue
        ga, gb = gc.get(a, 0.0), gc.get(b, 0.0)
        if abs(ga - gb) > gc_diff:
            continue
        if has_depth:
            da, db = depth.get(a, 0.0), depth.get(b, 0.0)
            if da <= 0 or db <= 0:
                continue
            if max(da, db) / min(da, db) > depth_ratio:
                continue
        kept[(a, b)] = w
    return kept


def propagate(edges, contig_ids, seeds, resolution=1.0, majority_ratio=10.0):
    """Louvain 聚类 + GMW >=ratio 多数比传播。

    返回 (node_verdict: {cid: verdict}, community_of: {cid: (tool内社区索引)},
          comm_stats: [{size, viral, nonviral, rule}])。
    种子与 conflicting 社区一律不改标 (GMW 同款语义)。
    """
    G = nx.Graph()
    G.add_nodes_from(sorted(contig_ids))
    for (a, b), w in sorted(edges.items()):
        G.add_edge(a, b, weight=w)
    communities = louvain_communities(G, weight='weight',
                                      resolution=resolution, seed=42)

    node_verdict, community_of, comm_stats = {}, {}, []
    for idx, comm in enumerate(sorted(communities, key=lambda c: min(c))):
        comm = set(comm)
        V = sum(1 for c in comm if seeds.get(c) == V_VIRAL_SEED)
        N = sum(1 for c in comm if seeds.get(c) == V_NONVIRAL_SEED)
        if V > 0 and V >= majority_ratio * N:
            rule, unknown_verdict = 'viral-dominant', V_GRAPH_VIRAL
        elif N > 0 and N >= majority_ratio * V:
            rule, unknown_verdict = 'nonviral-dominant', V_GRAPH_NONVIRAL
        elif V > 0 and N > 0:
            rule, unknown_verdict = 'conflicting', V_CONFLICTING
        else:
            rule, unknown_verdict = 'no-seeds', V_UNKNOWN
        comm_stats.append({'size': len(comm), 'viral': V, 'nonviral': N, 'rule': rule})
        for c in comm:
            community_of[c] = idx
            if c in seeds and seeds[c] in (V_VIRAL_SEED, V_NONVIRAL_SEED):
                node_verdict[c] = seeds[c]          # 种子不可变
            else:
                node_verdict[c] = unknown_verdict
    return node_verdict, community_of, comm_stats


def compute_adjustments(kept_ids, candidate_ids, cdd_unclassified, blastx_hit_ids,
                        verdicts, do_demote=True, do_rescue=True):
    """filter 末端调整的纯函数 (放在本模块以便单测; filter_virus.py 调用)。

    demote: 保留集中 cdd_call=='unclassified' 且无任何 blastx hit (零正面证据,
            仅因全局 no-hit 规则放行) 且图判定 graph_nonviral → 移出。
    rescue: 被整链丢弃 (候选-保留) 且图判定 graph_viral → 加回。
    返回 (final_ids, demoted, rescued)。
    """
    kept, cands = set(kept_ids), set(candidate_ids)
    demoted = set()
    rescued = set()
    if do_demote:
        zero_evidence = (kept & set(cdd_unclassified)) - set(blastx_hit_ids)
        demoted = {c for c in zero_evidence
                   if verdicts.get(c) == V_GRAPH_NONVIRAL}
    if do_rescue:
        rescued = {c for c in (cands - kept)
                   if verdicts.get(c) == V_GRAPH_VIRAL}
    return (kept - demoted) | rescued, demoted, rescued


def write_derived_candidates(ident_dir, n, demoted_ids, candidates_root):
    """派生候选 = ident 候选 - demoted。demoted 为空时逐字节复制。

    输出命名固定 {n}_virus.all.candidate.fasta (COBRA raw 模式第 4 探测模式;
    该目录不含前 3 个模式名, 不会被误认成 filter 产物)。
    返回 (src, dst, removed)。
    """
    src = None
    for name in (f"{n}_virus.all.candidate.fasta", f"{n}.virus.candidate.fasta"):
        p = Path(ident_dir) / n / name
        if p.is_file():
            src = p
            break
    if src is None:
        return None, None, 0
    dst_dir = Path(candidates_root) / n
    dst_dir.mkdir(parents=True, exist_ok=True)
    dst = dst_dir / f"{n}_virus.all.candidate.fasta"
    demoted_ids = set(demoted_ids)
    if not demoted_ids:
        shutil.copyfile(src, dst)
        return src, dst, 0
    removed = 0
    with _open_maybe_gz(src) as fin, open(dst, 'w') as fout:
        write_on = False
        for line in fin:
            if line.startswith('>'):
                cid = line[1:].split()[0] if len(line) > 1 else ''
                write_on = cid not in demoted_ids
                if not write_on:
                    removed += 1
            if write_on:
                fout.write(line)
    return src, dst, removed


# ═══════════════════ filter 输出裁决 (新主链: filter→graph→cobra) ═══════════════════

def _copy_fasta_subset(src, dst, keep_ids):
    """按 header 首 token 过滤复制 FASTA (与 filter_virus._copy_fasta_subset 同语义)"""
    with open(src) as fin, open(dst, 'w') as fout:
        write_on = False
        for line in fin:
            if line.startswith('>'):
                cid = line[1:].split()[0] if len(line) > 1 else ''
                write_on = cid in keep_ids
            if write_on:
                fout.write(line)


def _fasta_ids(path):
    """FASTA 头 ID 列表 (与 filter_virus._fasta_ids 同语义, 独立实现避免交叉导入)"""
    ids = []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                ids.append(line[1:].split()[0] if len(line) > 1 else "")
    return ids


def _resolve_filter_chain(filter_dir, variants):
    """在 d_filter 下按变体解析最终链文件 → (匹配名, 路径) 或 (None, None)"""
    for sv in variants:
        for pat in ("{sv}.virus.candidate_filtered.fasta",
                    "{sv}_virus.all.candidate_filtered.fasta",
                    "{sv}.virus.candidate.fasta"):
            p = os.path.join(filter_dir, sv, pat.format(sv=sv))
            if os.path.isfile(p):
                return sv, p
    return None, None


def adjudicate_sample(sample, out_root, filter_dir, ident_dir,
                      demote=True, rescue=True):
    """对 filter 输出做图裁决 (主链重排后 graph 的核心动作):
      final = filter 链 − (cdd unclassified ∩ 无 blastx hit ∩ graph_nonviral)
              ∪ (被 filter 丢弃 ∩ graph_viral)
    写: ① {out_root}/{sample}/{sample}.virus.candidate.final.fasta
        ② {filter_dir}/candidates/{sv}/{sv}_virus.all.candidate.fasta (覆盖直通镜像,
           COBRA 主链输入)
        ③ {out_root}/{sample}/{sample}.graph_adjust_calls.tsv (审计)
    返回 (n_final, n_demoted, n_rescued); filter 链缺失返回 None (graph 兜底镜像
    由调用方用 ident − graph_nonviral 处理)。"""
    vtsv = Path(out_root) / sample / f"{sample}.graph_verdict.tsv"
    if not vtsv.is_file():
        return None
    verdicts = {}
    with open(vtsv, encoding="utf-8") as fh:
        head = fh.readline().rstrip("\n").split("\t")
        ci = {c: i for i, c in enumerate(head)}
        vi, si = ci.get("verdict"), ci.get("contig")
        if vi is None or si is None:
            return None
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) > vi:
                verdicts[f[si]] = f[vi]

    variants = [sample]
    if sample.endswith("_clean"):
        variants.append(sample[:-len("_clean")])
    else:
        variants.append(sample + "_clean")
    sv, chain = _resolve_filter_chain(filter_dir, variants)
    if chain is None:
        return None
    chain_ids = [i for i in _fasta_ids(chain)]

    cdd_call = {}
    cdd_tsv = os.path.join(filter_dir, sv, "cdd", "cdd_calls.tsv")
    if os.path.isfile(cdd_tsv):
        with open(cdd_tsv) as fh:
            for i, line in enumerate(fh):
                if i == 0 and line.lower().startswith("contig"):
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) >= 2:
                    cdd_call[f[0]] = f[1]
    blastx_hits = set()
    blast_tsv = os.path.join(filter_dir, sv, "uniprot", "diamond_uniprot", "blastx.tsv")
    if os.path.isfile(blast_tsv):
        with open(blast_tsv) as fh:
            for line in fh:
                q = line.split("\t", 1)[0]
                if q:
                    blastx_hits.add(q)

    candidates_fa = None
    if ident_dir:
        for n in variants:
            for name in (f"{n}_virus.all.candidate.fasta", f"{n}.virus.candidate.fasta"):
                p = Path(ident_dir) / n / name
                if p.is_file():
                    candidates_fa = p
                    break
            if candidates_fa:
                break
    cand_ids = set(_fasta_ids(candidates_fa)) if candidates_fa else set(chain_ids)

    final, demoted, rescued = compute_adjustments(
        set(chain_ids), cand_ids,
        {c for c, call in cdd_call.items() if call == "unclassified"},
        blastx_hits, verdicts, do_demote=demote, do_rescue=rescue)

    # 序列源: 有救回时需含全量候选 (rescued 不在链里), 否则链即可
    src = chain
    if rescued:
        if candidates_fa is not None:
            src = str(candidates_fa)
        else:
            print(f"[graph] {sample}: WARN 无法定位候选源, 放弃 {len(rescued)} 条救回",
                  flush=True)
            final = final - rescued
            rescued = set()
    out_fa = Path(out_root) / sample / f"{sample}.virus.candidate.final.fasta"
    _copy_fasta_subset(src, str(out_fa), final)

    mirror_dir = Path(filter_dir) / "candidates" / sv
    mirror_dir.mkdir(parents=True, exist_ok=True)
    _copy_fasta_subset(src, str(mirror_dir / f"{sv}_virus.all.candidate.fasta"), final)

    audit = Path(out_root) / sample / f"{sample}.graph_adjust_calls.tsv"
    with open(audit, "w", encoding="utf-8") as fh:
        fh.write("contig\tcdd_call\tblastx_hit\tgraph_verdict\taction\n")
        for cid in sorted(cand_ids | set(chain_ids)):
            action = ("demoted" if cid in demoted else
                      "rescued" if cid in rescued else
                      "kept" if cid in final else
                      "dropped" if cid in cand_ids else "outside")
            fh.write("%s\t%s\t%s\t%s\t%s\n" % (
                cid, cdd_call.get(cid, ""),
                "yes" if cid in blastx_hits else "no",
                verdicts.get(cid, ""), action))
    return len(final), len(demoted), len(rescued)


# ═══════════════════════════ 外部工具封装 ═══════════════════════════

def _run(cmd, desc):
    """shell 命令 (管道需要), 失败抛 RuntimeError。仿 cobra_pipeline.run_command。"""
    proc = subprocess.run(cmd, shell=True, executable='/bin/bash',
                          capture_output=True, text=True)
    if proc.returncode != 0:
        tail = (proc.stderr or proc.stdout or '')[-2000:]
        raise RuntimeError(f"{desc} 失败 (rc={proc.returncode}): {tail}")
    return proc


def find_read_pair(reads_dir, sample, variants=None):
    """reads_dir 里找 {sample} 的 PE 对 (r1, r2) 或 SE (r1, None)。仿 cobra 命名匹配。"""
    for s in (variants or [sample]):
        r1 = _match_read(reads_dir, s, ('_1.', '.R1.', '_R1.'))
        r2 = _match_read(reads_dir, s, ('_2.', '.R2.', '_R2.'))
        if r1:
            return r1, r2
    return None, None


def _match_read(reads_dir, sample, tags):
    if not Path(reads_dir).is_dir():
        return None
    for f in sorted(os.listdir(reads_dir)):
        for t in tags:
            if f.startswith(sample + t):
                return Path(reads_dir) / f
    return None


def truncate_reads(src, dst, max_reads):
    """确定性前缀截取前 max_reads 条记录 (FASTQ 按固定 4 行组, FASTA 按记录)。

    FASTQ 不能按 '@' 前缀判断记录头 (序列/质量行也可能以 @ 开头), 必须按 4 行组。
    """
    with _open_maybe_gz(src) as fin, open(dst, 'w') as fout:
        head = fin.readline()
        if not head:
            return dst
        if head.startswith('@'):
            n, group = 0, [head]
            for line in fin:
                group.append(line)
                if len(group) == 4:
                    if n >= max_reads:
                        break
                    n += 1
                    fout.writelines(group)
                    group = []
            return dst
        # FASTA: '>' 开头为新记录
        n, keep = 0, True
        for line in itertools.chain([head], fin):
            if line.startswith('>'):
                if n >= max_reads:
                    break
                n += 1
                keep = True
            elif n >= max_reads:
                keep = False
            if keep:
                fout.write(line)
    return dst


# ═══════════════════════════ 样本级流程 ═══════════════════════════

def find_contig_files(asm_sample_dir):
    """{S}/ 顶层 {S}_{tool}.contig.fasta (非递归, 不碰 refineC 子目录产物)。"""
    s = Path(asm_sample_dir).name
    d = Path(asm_sample_dir)
    return sorted(p for p in d.iterdir()
                  if p.is_file() and p.name.startswith(s + '_')
                  and p.name.endswith('.contig.fasta'))


def tool_from_contig_file(path, sample):
    name = Path(path).name
    return name[len(sample) + 1:-len('.contig.fasta')]


def write_filtered_reference(contig_fa, lengths_gc, min_len, workdir):
    """映射参考裁剪: 只保留长度 >= min_len 的 contigs (246 需求: 只分析 unknown 时
    省映射时间 — 短 contig (<1000bp) 不可能是种子、判定不可执行(不注入候选池),
    从参考集剔除可将映射时间省 30-60%, 可执行判定的证据不受损)。"""
    eligible = {c for c, (l, _) in lengths_gc.items() if l >= min_len}
    dst = Path(workdir) / "ref_mapminlen.fa"
    dst.parent.mkdir(parents=True, exist_ok=True)
    kept = 0
    with _open_maybe_gz(contig_fa) as fin, open(dst, "w") as fout:
        keep = False
        for line in fin:
            if line.startswith(">"):
                cid = line[1:].split()[0] if len(line) > 1 else ""
                keep = cid in eligible
                if keep:
                    kept += 1
            if keep:
                fout.write(line)
    return dst, kept


def map_reads_and_links(contig_fa, r1, r2, workdir, threads,
                        require_proper=True, max_reads=0):
    """bwa-mem2 index + mem | samtools sort + index → (cross_links Counter, depth dict)。

    contig 复制进 workdir 再建索引, 不污染 d_asm; BAM 用后由调用方清理。
    """
    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    ref = workdir / (Path(contig_fa).stem + '.contigs.fa')
    shutil.copyfile(contig_fa, ref)

    if max_reads and max_reads > 0:
        r1 = truncate_reads(r1, workdir / (Path(r1).name + '.trunc'), max_reads)
        if r2:
            r2 = truncate_reads(r2, workdir / (Path(r2).name + '.trunc'), max_reads)

    bam = workdir / 'mapped.sorted.bam'
    r2_arg = f"'{r2}'" if r2 else ''
    _run(f"bwa-mem2 index '{ref}'", 'bwa-mem2 index')
    _run(f"bwa-mem2 mem -t {threads} '{ref}' '{r1}' {r2_arg} "
         f"| samtools sort -@ {threads} -o '{bam}' -", 'bwa-mem2 mem | samtools sort')
    _run(f"samtools index '{bam}'", 'samtools index')

    # 跨 contig 检测默认不要求 proper: bwa 的 0x2 仅在同 contig 时置位,
    # 跨 contig pair 永远 non-proper (246 实测 proper 跨 contig = 0/196079)
    flags_f = '0x43' if require_proper else '0x41'
    # 流式解析 (246 实测: XU 样本数千万行 SAM, capture_output 全量入内存耗时 40+ 分钟
    # 且峰值内存数 GB; Popen 逐行迭代常数内存)
    proc = subprocess.Popen(f"samtools view -f {flags_f} -F 0x90C '{bam}'",
                            shell=True, executable='/bin/bash',
                            stdout=subprocess.PIPE, stderr=subprocess.DEVNULL,
                            text=True, bufsize=1 << 20)
    try:
        links = count_cross_links(proc.stdout, require_proper=require_proper)
    finally:
        proc.stdout.close()
        if proc.wait() != 0:
            raise RuntimeError(f"samtools view 失败 (rc={proc.returncode})")

    cov = subprocess.run(f"coverm contig -b '{bam}' -m mean "
                         f"--min-covered-fraction 0 -t {threads}",
                         shell=True, executable='/bin/bash',
                         capture_output=True, text=True)
    depth = {}
    if cov.returncode == 0:
        for line in cov.stdout.splitlines():
            if not line or line.startswith('Contig'):
                continue
            f = line.split('\t')
            if len(f) >= 2:
                try:
                    depth[f[0]] = float(f[1])
                except ValueError:
                    continue
    else:
        print(f"[graph] WARN coverm 失败, 跳过深度断边: {(cov.stderr or '')[-500:]}",
              flush=True)
    return links, depth


def process_tool_graph(contig_fa, sample, tool, lengths_gc, seeds, votes,
                       workdir, threads, r1, r2, cfg):
    """单个 (样本, 组装工具) 的图构建 → (rows, demoted_ids, summary_rows)。"""
    contig_ids = list(lengths_gc)
    if not contig_ids:
        return [], set(), []

    # 映射参考裁剪 (--map-min-len): 短 contig 无可执行判定, 从参考集剔除省映射时间
    map_fa = contig_fa
    if cfg.get('map_min_len', 0) > 0:
        n_short = sum(1 for c in contig_ids
                      if lengths_gc[c][0] < cfg['map_min_len'])
        if n_short == len(contig_ids):
            r1 = None   # 全部低于阈值: 无可执行判定, 跳过映射
            print(f"[graph] {sample}/{tool}: 全部 {n_short} contigs 低于 "
                  f"map-min-len={cfg['map_min_len']}, 跳过映射", flush=True)
        elif n_short > 0:
            map_fa, kept = write_filtered_reference(
                contig_fa, lengths_gc, cfg['map_min_len'], workdir)
            print(f"[graph] {sample}/{tool}: 参考裁剪 {len(contig_ids)}→{kept} "
                  f"(剔除 {n_short} 条 <{cfg['map_min_len']}bp)", flush=True)

    links, depth = Counter(), {}
    if r1:
        links, depth = map_reads_and_links(
            map_fa, r1, r2, workdir, threads,
            require_proper=cfg['require_proper'], max_reads=cfg['max_reads'])
    else:
        print(f"[graph] {sample}/{tool}: 无 PE reads, 跳过建图 (全 unknown)", flush=True)

    edges = prune_edges(links,
                        {c: lengths_gc[c][1] for c in contig_ids}, depth,
                        min_link_pairs=cfg['min_link_pairs'],
                        gc_diff=cfg['gc_diff'], depth_ratio=cfg['depth_ratio'],
                        has_depth=bool(depth))
    verdicts, community_of, comm_stats = propagate(
        edges, contig_ids, seeds, resolution=cfg['resolution'],
        majority_ratio=cfg['majority_ratio'])

    demoted = {c for c, v in verdicts.items() if v == V_GRAPH_NONVIRAL}
    rows = []
    for cid in sorted(contig_ids):
        length, gcf = lengths_gc[cid]
        cidx = community_of.get(cid)
        stats = comm_stats[cidx] if cidx is not None else {
            'size': 1, 'viral': 0, 'nonviral': 0, 'rule': 'no-links'}
        v = verdicts.get(cid, V_UNKNOWN)
        if v == V_VIRAL_SEED:
            reason = f"seed: votes={votes.get(cid, 0)}"
        elif v == V_NONVIRAL_SEED:
            reason = (f"seed: 0 votes, len>={cfg['nonviral_min_len']}, examined")
        elif v == V_GRAPH_VIRAL or v == V_GRAPH_NONVIRAL:
            reason = (f"community {tool}:c{cidx} {stats['rule']}: "
                      f"V={stats['viral']},N={stats['nonviral']}")
        elif v == V_CONFLICTING:
            reason = (f"community {tool}:c{cidx} conflicting: "
                      f"V={stats['viral']},N={stats['nonviral']}")
        elif cid not in community_of:
            reason = "no links"
        else:
            reason = f"community {tool}:c{cidx} {stats['rule']}"
        rows.append({
            'contig': cid, 'source_tool': tool, 'length': length,
            'gc': round(gcf, 4), 'depth': round(depth.get(cid, 0.0), 4),
            'votes': votes.get(cid, 0),
            'seed_label': seeds.get(cid, V_UNKNOWN),
            'community_id': f"{tool}:c{cidx}" if cidx is not None else '',
            'community_size': stats['size'],
            'community_viral': stats['viral'],
            'community_nonviral': stats['nonviral'],
            'verdict': v, 'reason': reason,
        })
    summary_rows = [{
        'source_tool': tool, 'community_id': f"{tool}:c{i}",
        'size': st['size'], 'viral_seed': st['viral'],
        'nonviral_seed': st['nonviral'],
        'graph_viral': sum(1 for r in rows if r['community_id'] == f"{tool}:c{i}"
                           and r['verdict'] == V_GRAPH_VIRAL),
        'graph_nonviral': sum(1 for r in rows if r['community_id'] == f"{tool}:c{i}"
                              and r['verdict'] == V_GRAPH_NONVIRAL),
        'conflicting': sum(1 for r in rows if r['community_id'] == f"{tool}:c{i}"
                           and r['verdict'] == V_CONFLICTING),
        'unknown': sum(1 for r in rows if r['community_id'] == f"{tool}:c{i}"
                       and r['verdict'] == V_UNKNOWN),
        'rule': st['rule'],
    } for i, st in enumerate(comm_stats)]
    return rows, demoted, summary_rows


def process_sample(sample, asm_dir, ident_dir, reads_dir, out_root, cfg,
                   threads, emit_candidates=True):
    """单样本全流程。返回 (verdict_tsv, demoted_total, derived_dst)。"""
    sdir = Path(asm_dir) / sample
    out_sdir = Path(out_root) / sample
    out_sdir.mkdir(parents=True, exist_ok=True)
    done_marker = out_sdir / f"{sample}.graph.done"

    contig_files = find_contig_files(sdir)
    if not contig_files:
        print(f"[graph] {sample}: 无 contig 文件, 跳过", flush=True)
        return None, set(), None

    variants = [sample]
    if sample.endswith('_clean'):
        variants.append(sample[:-len('_clean')])
    variants.append(sample + '_clean')
    votes, viroid, n, ident_sdir = load_seed_votes(ident_dir, variants)
    if n is None:
        print(f"[graph] {sample}: d_ident 无种子清单, 全部按 unknown 处理", flush=True)

    r1, r2 = find_read_pair(reads_dir, sample, variants)
    if not r1:
        print(f"[graph] {sample}: reads_dir 未找到 reads, 全部按 unknown 处理", flush=True)

    all_rows, all_summary, demoted_all = [], [], set()
    for cf in contig_files:
        tool = tool_from_contig_file(cf, sample)
        lengths_gc = read_fasta_lengths_gc(cf)
        if n is None:
            # 无 ident 种子清单: 不能断言任何 contig 非病毒, 全部按 unknown
            seeds = {cid: V_UNKNOWN for cid in lengths_gc}
        else:
            seeds = assign_seeds(lengths_gc.keys(),
                                 {c: l for c, (l, _) in lengths_gc.items()},
                                 votes, viroid,
                                 seed_min_votes=cfg['seed_min_votes'],
                                 nonviral_min_len=cfg['nonviral_min_len'])
        workdir = out_sdir / 'tmp' / tool
        rows, demoted, summary = process_tool_graph(
            cf, sample, tool, lengths_gc, seeds, votes, workdir, threads,
            r1, r2, cfg)
        all_rows.extend(rows)
        all_summary.extend(summary)
        demoted_all |= demoted
        if cfg['keep_bam']:
            for keep in workdir.glob('*.bam*'):
                shutil.move(str(keep), str(out_sdir / keep.name))
        shutil.rmtree(workdir, ignore_errors=True)

    verdict_tsv = out_sdir / f"{sample}.graph_verdict.tsv"
    with open(verdict_tsv, 'w', encoding='utf-8') as fh:
        fh.write('\t'.join(VERDICT_COLUMNS) + '\n')
        for r in all_rows:
            fh.write('\t'.join(str(r[c]) for c in VERDICT_COLUMNS) + '\n')
    summary_tsv = out_sdir / f"{sample}.graph_summary.tsv"
    with open(summary_tsv, 'w', encoding='utf-8') as fh:
        fh.write('\t'.join(SUMMARY_COLUMNS) + '\n')
        for r in all_summary:
            fh.write('\t'.join(str(r[c]) for c in SUMMARY_COLUMNS) + '\n')

    derived_dst = None
    if emit_candidates and n is not None:
        cand_root = Path(out_root) / 'candidates'
        _, derived_dst, removed = write_derived_candidates(
            ident_dir, n, demoted_all, cand_root)
        if derived_dst:
            print(f"[graph] {sample}: 派生候选 {derived_dst.name} "
                  f"(移除 {removed} 条 graph_nonviral)", flush=True)

    done_marker.write_text(
        f"contigs={len(all_rows)} demoted={len(demoted_all)}\n", encoding='utf-8')
    n_comm = len(all_summary)
    print(f"[graph] {sample}: {len(contig_files)} 个组装文件, "
          f"{len(all_rows)} contigs, {n_comm} 社区, "
          f"降级 {len(demoted_all)} 条", flush=True)
    return verdict_tsv, demoted_all, derived_dst


# ═══════════════════════════ 入口 ═══════════════════════════

def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description='组装图标签传播 (graph 阶段): 配对读段链接图 + ensemble 种子 + Louvain',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument('--asm-dir', default=str(_D['d_asm']),
                   help='组装目录 (d_asm), 样本子目录含 {S}_{tool}.contig.fasta')
    p.add_argument('--ident-dir', default=str(_D['d_ident']),
                   help='鉴定目录 (d_ident), 种子来源 *_virus.*.result.id')
    p.add_argument('--reads-dir', default=str(_D['d_hostdep']),
                   help='clean reads 目录 (d_hostdep 或 co-assembly 合并 reads)')
    p.add_argument('-o', '--output-dir', default=str(_D['d_graph']),
                   help='输出目录 (d_graph)')
    p.add_argument('-t', '--threads', type=int, default=8, help='单样本映射线程数')
    p.add_argument('-j', '--jobs', type=int, default=1, help='并行样本数')
    p.add_argument('--min-link-pairs', type=int, default=3)
    p.add_argument('--gc-diff', type=float, default=0.20)
    p.add_argument('--depth-ratio', type=float, default=20.0)
    p.add_argument('--require-proper', action='store_true',
                   help='仅收 proper pair (0x2)。注意: bwa 只对同 contig 且方向/距离正常的 '
                        'pair 置 proper 标志, 跨 contig pair 永远 non-proper — 默认关, '
                        '开启将得不到任何跨 contig 链接 (246 实测 proper 跨 contig = 0/196079)')
    p.add_argument('--resolution', type=float, default=1.0)
    p.add_argument('--seed-min-votes', type=int, default=2,
                   help='病毒种子最少工具票数 (默认 2; 弱阳性票数<阈值归 unknown 交图裁决')
    p.add_argument('--nonviral-min-len', type=int, default=1000)
    p.add_argument('--majority-ratio', type=float, default=10.0)
    p.add_argument('--max-reads', type=int, default=0,
                   help='>0: 每端取前 N 条 (确定性前缀截取, 规模兜底)')
    p.add_argument('--map-min-len', type=int, default=0,
                   help='>0: 映射参考只保留 >= 该长度的 contigs (推荐 1000; 短 contig '
                        '无可执行判定, 剔除可省 30-60%% 映射时间, unknown 裁决不受损)')
    p.add_argument('--keep-bam', action='store_true', help='保留映射 BAM 供调试')
    p.add_argument('--filter-dir', default=str(_D['d_filter']),
                   help='filter 输出目录 (d_filter); 裁决其最终候选集并更新 COBRA 镜像')
    p.add_argument('--no-demote', action='store_true',
                   help='禁用降级 (graph_nonviral 零证据保留 contig 不移出)')
    p.add_argument('--no-rescue', action='store_true',
                   help='禁用救回 (graph_viral 被丢弃 contig 不加回)')
    p.add_argument('--no-candidates', action='store_true',
                   help='不写派生候选目录 (COBRA 拦截将回退 d_ident)')
    p.add_argument('--force', action='store_true', help='忽略 {S}.graph.done 重跑')
    p.add_argument('--samples', nargs='*', default=None,
                   help='只处理指定样本 (默认扫描 asm-dir 全部)')
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    cfg = {'min_link_pairs': args.min_link_pairs, 'gc_diff': args.gc_diff,
           'depth_ratio': args.depth_ratio, 'require_proper': args.require_proper,
           'resolution': args.resolution, 'seed_min_votes': args.seed_min_votes,
           'nonviral_min_len': args.nonviral_min_len,
           'majority_ratio': args.majority_ratio, 'max_reads': args.max_reads,
           'map_min_len': args.map_min_len, 'keep_bam': args.keep_bam,
           'filter_dir': args.filter_dir,
           'demote': not args.no_demote, 'rescue': not args.no_rescue}
    out_root = Path(args.output_dir)
    out_root.mkdir(parents=True, exist_ok=True)

    asm_dir = Path(args.asm_dir)
    samples = args.samples
    if not samples:
        samples = []
        for d in sorted(asm_dir.iterdir()) if asm_dir.is_dir() else []:
            if d.is_dir() and find_contig_files(d):
                samples.append(d.name)
    if not samples:
        print("[graph] asm-dir 下未发现含 contig 的样本目录, 退出", flush=True)
        return

    pending = []
    for s in samples:
        if not args.force and (out_root / s / f"{s}.graph.done").is_file():
            print(f"[graph] {s}: 已完成, 跳过 (--force 重跑)", flush=True)
            continue
        pending.append(s)
    print(f"[graph] 样本 {len(pending)}/{len(samples)} 待处理 "
          f"(threads={args.threads}, jobs={args.jobs})", flush=True)

    failed = []
    with ThreadPoolExecutor(max_workers=max(1, args.jobs)) as ex:
        futs = {ex.submit(process_sample, s, args.asm_dir, args.ident_dir,
                          args.reads_dir, out_root, cfg, args.threads,
                          not args.no_candidates): s for s in pending}
        for fut in as_completed(futs):
            s = futs[fut]
            try:
                fut.result()
            except Exception as exc:  # 样本级容错, 与管线主循环口径一致
                failed.append(s)
                print(f"[graph] {s}: FAIL {exc}", flush=True)

    # 完整性兜底: d_ident 中每个有候选的样本都必须有派生文件 (COBRA 发现零失败)
    if not args.no_candidates:
        cand_root = out_root / 'candidates'
        ident_dir = Path(args.ident_dir)
        if ident_dir.is_dir():
            for d in sorted(ident_dir.iterdir()):
                if not d.is_dir():
                    continue
                has_cand = any((d / f"{d.name}{suffix}").is_file() for suffix in
                               ('_virus.all.candidate.fasta', '.virus.candidate.fasta'))
                if not has_cand:
                    continue
                dst = cand_root / d.name / f"{d.name}_virus.all.candidate.fasta"
                if not dst.is_file():
                    write_derived_candidates(ident_dir, d.name, set(), cand_root)
                    print(f"[graph] 兜底复制派生候选: {d.name}", flush=True)

    # ── filter 输出裁决 (新主链 filter→graph): 救回/降级 + 更新 COBRA 镜像 ──
    if cfg['demote'] or cfg['rescue']:
        filt = cfg.get('filter_dir')
        if filt and os.path.isdir(filt):
            n_adj = n_d = n_r = 0
            for vtsv in sorted(out_root.glob('*/*.graph_verdict.tsv')):
                sname = vtsv.parent.name
                r = adjudicate_sample(sname, str(out_root), str(filt),
                                      str(Path(args.ident_dir)),
                                      demote=cfg['demote'], rescue=cfg['rescue'])
                if r:
                    n_adj += 1
                    n_d += r[1]
                    n_r += r[2]
            print(f"[graph] filter 输出裁决: {n_adj} 样本, 降级 {n_d}, 救回 {n_r}",
                  flush=True)

    produced = list(out_root.glob('*/[!t]*.graph_verdict.tsv')) or \
               list(out_root.glob('*/*.graph_verdict.tsv'))
    if pending and not produced:
        print("[graph] FAIL: 无任何样本产出 verdict TSV", flush=True)
        sys.exit(1)
    if failed:
        print(f"[graph] 部分样本失败: {', '.join(failed)}", flush=True)
    print(f"[graph] 完成: verdict TSV {len(produced)} 个", flush=True)


if __name__ == '__main__':
    main()
