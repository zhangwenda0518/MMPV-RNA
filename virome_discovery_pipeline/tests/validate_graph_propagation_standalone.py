#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
validate_graph_propagation_standalone.py — graph 阶段独立运行验证 (无外部映射工具)
====================================================================================
回答两个问题: ① graph_propagation.py 能否脱离编排器独立运行; ② 判定语义是否成立。

场景 A (自建混合感染社区, ground truth 完全已知, 纯 Python):
  2 个病毒基因组链 + 宿主链 + 细菌链 + 高拷贝质粒 + 弱阳性 contig ×2 + 短病毒片段。
  种子语义对齐真实管线: 真病毒 contig ≥1 票 (弱阳性→unknown 交图裁决), 强种子 ≥2 票;
  0 票长 contig = 全体工具阴性 → nonviral_seed。各埋 3 条跨基因组噪声边, 分别被
  权重/GC差/深度比 剪掉。断言逐 contig verdict + 社区不跨基因组。

场景 B (GMW 论文示例数据 = 真实 SPAdes 4.0 组装图, H3N2+H5N6 混合感染, 377 段):
  解析 input.gfa → 本机 blastn 对两亚型参考定真值 → 沿真实 L 连接模拟跨 contig
  配对 (DP 覆盖度作 depth) → process_sample 全流程。指标:
    ① 剪枝后真实边存活率; ② 短片段拓扑救回率 (unknown 段落进有种子社区 → graph_viral);
    ③ graph_viral 段的 blastn 复核 (纯流感图中应≈100% 确认, 无命中者=拓扑发现的价值);
    ④ 社区亚型纯度 (描述性: 拓扑是否自发分离 H3N2/H5N6)。

用法:
  python validate_graph_propagation_standalone.py [--gmw-repo PATH] [--skip-blast]

依赖: python>=3.10 + networkx; 场景 B 另需 blastn + makeblastdb (本机已装)。
说明: 真实 bwa-mem2/samtools/coverm 映射在 linux 服务器 pixi 环境运行; 本脚本以
      "模拟配对 → SAM 文本 → 真实 count_cross_links" 替代, 模块其余代码路径
      (process_sample/种子加载/剪枝/Louvain/TSV/派生候选) 全部真实执行。
"""

import argparse
import random
import shutil
import subprocess
import sys
import tempfile
from collections import Counter
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))
sys.path.insert(0, str(REPO_ROOT / "virome_discovery_pipeline"))

import graph_propagation as gp  # noqa: E402

CFG = {'min_link_pairs': 3, 'gc_diff': 0.20, 'depth_ratio': 20.0,
       'require_proper': True, 'resolution': 1.0, 'seed_min_votes': 2,
       'nonviral_min_len': 1000, 'majority_ratio': 10.0, 'max_reads': 0,
       'keep_bam': False}


# ── 公共工具 ─────────────────────────────────────────────────────────────

def sam_pair(r1_name, r2_name, count):
    """count 条跨 contig proper pair 的 SAM 行 (flag 99, read1, RNEXT=mate)。"""
    return [f"r{i}\t99\t{r1_name}\t100\t60\t150M\t{r2_name}\t300\t0\tA\t*\n"
            for i in range(count)]


def gc_seq(length, gc, rng):
    """按目标 GC 含量生成随机序列。"""
    n_gc = int(round(gc * length))
    n_at = length - n_gc
    bases = ['G'] * (n_gc // 2) + ['C'] * (n_gc - n_gc // 2) + \
            ['A'] * (n_at // 2) + ['T'] * (n_at - n_at // 2)
    rng.shuffle(bases)
    return ''.join(bases)


def patch_mapper(links, depth):
    """把 map_reads_and_links (bwa+samtools+coverm) 替换为预计算结果。"""
    gp.map_reads_and_links = lambda *a, **kw: (links, depth)


def run_process_sample(sample, asm_dir, ident_dir, reads_dir, out_root):
    tsv, demoted, derived = gp.process_sample(
        sample, asm_dir, ident_dir, reads_dir, out_root, CFG, 1,
        emit_candidates=True)
    verdicts, comm = {}, {}
    with open(tsv, encoding='utf-8') as fh:
        head = fh.readline().rstrip('\n').split('\t')
        ci = {c: i for i, c in enumerate(head)}
        for line in fh:
            f = line.rstrip('\n').split('\t')
            verdicts[f[ci['contig']]] = f[ci['verdict']]
            comm[f[ci['contig']]] = f[ci['community_id']]
    return verdicts, comm, demoted, derived


def make_ident(ident_dir, sample, candidates, strong):
    """result.id ×2 (genomad=全部候选, virbot=强命中) + 派生候选源 FASTA。"""
    (ident_dir / sample).mkdir(parents=True, exist_ok=True)
    with open(ident_dir / sample / f"{sample}_virus.genomad.result.id", 'w') as fh:
        fh.write("".join(c + "\n" for c in candidates))
    with open(ident_dir / sample / f"{sample}_virus.virbot.result.id", 'w') as fh:
        fh.write("".join(c + "\n" for c in strong))
    with open(ident_dir / sample / f"{sample}_virus.all.candidate.fasta", 'w') as fh:
        for c in candidates:
            fh.write(f">{c}\n{'ACGT' * 50}\n")


def setup_scaffold(td, sample):
    (td / "asm" / sample).mkdir(parents=True)
    (td / "ident").mkdir(parents=True)
    reads = td / "reads"
    reads.mkdir()
    (reads / f"{sample}_1.fa").write_text("")   # 占位 (映射被 patch, 不真正读)
    (reads / f"{sample}_2.fa").write_text("")
    return reads


# ── 场景 A: 合成混合感染社区 ─────────────────────────────────────────────

def scenario_a():
    print("=" * 70)
    print("场景 A: 合成混合感染社区 (ground truth 已知, 纯 Python)")
    print("=" * 70)
    rng = random.Random(42)
    td = Path(tempfile.mkdtemp(prefix="gp_val_A_"))
    sample = "S1_clean"
    reads = setup_scaffold(td, sample)

    # (id, len, gc, depth, votes) — 真病毒链 ≥1 票 (弱阳性), 0 票长 contig = 阴性
    contigs = [
        ("V1_c0", 2000, .45, 10, 2), ("V1_c1", 2200, .44, 10, 1),
        ("V1_c2", 1800, .46, 10, 1), ("V1_c3", 2000, .45, 10, 1),
        ("V1_s0", 300, .45, 10, 0),                     # 短病毒片段 (未被工具检查)
        ("V2_c0", 2000, .42, 8, 2), ("V2_c1", 2000, .42, 8, 1),
        ("V2_c2", 2000, .41, 8, 1), ("V2_c3", 1500, .43, 8, 1),
        ("Host_h0", 3000, .38, 12, 0), ("Host_h1", 2500, .38, 12, 0),
        ("Host_h2", 2000, .39, 12, 0), ("Host_hs0", 400, .38, 12, 0),
        ("Bact_b0", 2000, .66, 11, 0), ("Bact_b1", 2000, .65, 11, 0),
        ("Plasmid_p0", 2000, .45, 300, 0),              # 高拷贝污染 (GC 与 V1 相近)
    ]
    fa = td / "asm" / sample / f"{sample}_megahit.contig.fasta"
    with open(fa, 'w') as fh:
        for cid, ln, gc, _, _ in contigs:
            fh.write(f">{cid}\n{gc_seq(ln, gc, rng)}\n")

    candidates = [c for c, _, _, _, v in contigs if v >= 1]
    strong = [c for c, _, _, _, v in contigs if v >= 2]
    make_ident(td / "ident", "S1", candidates, strong)

    # 真实边: 基因组内部链 (权重 30); 噪声边: 三条各被一条规则剪掉
    sam = []
    for a, b in [("V1_c0", "V1_c1"), ("V1_c1", "V1_c2"), ("V1_c2", "V1_c3"),
                 ("V1_c2", "V1_s0"),
                 ("V2_c0", "V2_c1"), ("V2_c1", "V2_c2"), ("V2_c2", "V2_c3"),
                 ("Host_h0", "Host_h1"), ("Host_h1", "Host_h2"),
                 ("Host_h1", "Host_hs0"),
                 ("Bact_b0", "Bact_b1")]:
        sam += sam_pair(a, b, 30)
    sam += sam_pair("V1_c3", "Host_h2", 2)      # 噪声1: 权重 <3 → 剪
    sam += sam_pair("V2_c3", "Bact_b1", 10)     # 噪声2: GC差 |.43-.65|=.22 → 剪
    sam += sam_pair("V1_c2", "Plasmid_p0", 20)  # 噪声3: 深度比 300/10=30 → 剪
    links = gp.count_cross_links(sam, require_proper=True)
    depth = {cid: d for cid, _, _, d, _ in contigs}
    patch_mapper(links, depth)

    verdicts, comm, demoted, derived = run_process_sample(
        sample, str(td / "asm"), str(td / "ident"), str(reads), str(td / "graph"))

    expect = {}
    for cid, ln, _, _, v in contigs:
        if v >= 2:
            expect[cid] = gp.V_VIRAL_SEED
        elif v >= 1:
            expect[cid] = gp.V_GRAPH_VIRAL      # 弱阳性落在自己基因组的 viral-dominant 社区
        elif cid.startswith("Host"):
            expect[cid] = (gp.V_NONVIRAL_SEED if ln >= 1000 else gp.V_GRAPH_NONVIRAL)
        else:
            expect[cid] = gp.V_NONVIRAL_SEED
    expect["V1_s0"] = gp.V_GRAPH_VIRAL          # 短病毒片段靠拓扑救回

    ok = True
    for cid, want in expect.items():
        got = verdicts.get(cid, '<missing>')
        mark = "OK " if got == want else "FAIL"
        if got != want:
            ok = False
        print(f"  [{mark}] {cid:12s} 期望 {want:16s} 实得 {got}")
    # 社区断言: 噪声边全被剪 → 病毒社区与宿主/细菌/质粒社区零交集
    viral_comms = {comm[c] for c in ("V1_c0", "V1_c1", "V1_c2", "V1_c3", "V1_s0",
                                     "V2_c0", "V2_c1", "V2_c2", "V2_c3")}
    other_comms = ({comm[c] for c in ("Host_h0", "Host_h1", "Host_h2", "Host_hs0")} |
                   {comm[c] for c in ("Bact_b0", "Bact_b1")} | {comm["Plasmid_p0"]})
    if viral_comms & other_comms:
        ok = False
        print("  [FAIL] 存在跨基因组社区混合 — 剪枝规则未生效")
    else:
        print("  [OK ] 病毒/宿主/细菌/质粒 社区零混合 (3 条噪声边分别被 权重/GC/深度 剪掉)")
    n_comm = len(set(comm.values()))
    print(f"  demoted(候选剔除)={len(demoted)}, 社区数={n_comm}")
    print(f"  ⇒ 场景 A: {'全部通过' if ok else '存在失败'}\n")
    shutil.rmtree(td, ignore_errors=True)
    return ok


# ── 场景 B: GMW 论文示例数据 ─────────────────────────────────────────────

def parse_gfa(path):
    segs, links = {}, []
    with open(path, encoding='utf-8') as fh:
        for line in fh:
            f = line.rstrip('\n').split('\t')
            if f[0] == 'S':
                dp = 1.0
                for tag in f[3:]:
                    if tag.startswith('DP:f:'):
                        dp = float(tag[5:])
                segs[f[1]] = (f[2], dp)
            elif f[0] == 'L':
                links.append((f[1], f[3]))
    return segs, links


def blastn_subtypes(seqs, ref_dir, tmp):
    """blastn 全部 contig vs merge.fas → {id: (亚型, pident, 是否歧义)}。"""
    fa = tmp / "segs.fa"
    with open(fa, 'w') as fh:
        for sid, (seq, _) in seqs.items():
            fh.write(f">{sid}\n{seq}\n")
    db = tmp / "refdb"
    subprocess.run(["makeblastdb", "-in", str(ref_dir / "merge.fas"),
                    "-dbtype", "nucl", "-out", str(db)],
                   capture_output=True, check=True)
    out = subprocess.run(
        ["blastn", "-query", str(fa), "-db", str(db),
         "-outfmt", "6 qseqid sseqid pident bitscore", "-evalue", "1e-10",
         "-max_target_seqs", "10"], capture_output=True, text=True,
        check=True).stdout

    def st_of(s):
        return "H3N2" if "H3N2" in s else ("H5N6" if "H5N6" in s else "None")

    hits = {}
    for line in out.splitlines():
        q, s, pid, bs = line.split('\t')[:4]
        hits.setdefault(q, []).append((float(bs), st_of(s), float(pid)))
    sub = {}
    for sid in seqs:
        hh = sorted(hits.get(sid, []), reverse=True)
        if not hh:
            sub[sid] = ("None", 0.0, False)
            continue
        bs, st, pid = hh[0]
        near = {s2 for b2, s2, _ in hh if b2 >= bs * 0.95}
        sub[sid] = (st, pid, len(near) > 1)
    return sub


def scenario_b(gmw_repo, use_blast=True):
    print("=" * 70)
    print("场景 B: GMW 论文示例数据 (真实 SPAdes 组装图, H3N2+H5N6 混合感染)")
    print("=" * 70)
    gfa = Path(gmw_repo) / "examples" / "input.gfa"
    ref_dir = Path(gmw_repo) / "examples" / "ref"
    segs, links = parse_gfa(gfa)
    long_segs = {s for s, (q, _) in segs.items() if len(q) >= 1000}
    print(f"  GFA: {len(segs)} 段 / {len(links)} 真实连接; ≥1000bp 段 {len(long_segs)} 个")

    td = Path(tempfile.mkdtemp(prefix="gp_val_B_"))
    sample = "GMWex"
    reads = setup_scaffold(td, sample)

    # 1) 亚型 ground truth (本机 blastn)
    if use_blast:
        try:
            sub = blastn_subtypes(segs, ref_dir, td)
        except (FileNotFoundError, subprocess.CalledProcessError) as exc:
            print(f"  [WARN] blastn 不可用 ({exc}), 亚型真值置 None (仅拓扑验证)")
            sub = {s: ("None", 0.0, False) for s in segs}
    else:
        sub = {s: ("None", 0.0, False) for s in segs}
    by = Counter(st for st, _, _ in sub.values())
    n_amb = sum(1 for _, (_, _, a) in sub.items() if a)
    print(f"  blastn 亚型: {dict(by)}; 歧义段 {n_amb} 个 (双亚型近优命中, 排除出统计)")

    # 2) 组装 FASTA (contig id = GFA 段 id)
    fa = td / "asm" / sample / f"{sample}_spades.contig.fasta"
    with open(fa, 'w') as fh:
        for sid, (seq, _) in segs.items():
            fh.write(f">{sid}\n{seq}\n")

    # 3) 种子 = 全部 ≥1000bp 且 blastn 可分类的段 (模拟 10 工具面板对流感段的阳性判定;
    #    票数 2 = 强种子)。短片段 <1000bp 未被工具检查 → unknown, 交图裁决。
    candidates = sorted(long_segs - {s for s in long_segs if sub[s][0] == "None"
                                     or sub[s][2]})
    make_ident(td / "ident", sample, candidates, candidates)
    print(f"  种子: viral_seed={len(candidates)} (≥1000bp 且 blastn 可分类), "
          f"短片段 unknown={len(segs) - len(long_segs)}")

    # 4) 沿真实 L 连接模拟跨 contig 配对 (每连接 30 对), DP 覆盖度作 depth
    sam = []
    for a, b in links:
        sam += sam_pair(a, b, 30)
    links_c = gp.count_cross_links(sam, require_proper=True)
    depth = {s: dp for s, (_, dp) in segs.items()}
    patch_mapper(links_c, depth)

    # 5) 全流程
    verdicts, comm, demoted, derived = run_process_sample(
        sample, str(td / "asm"), str(td / "ident"), str(reads), str(td / "graph"))
    vc = Counter(verdicts.values())
    print(f"  判定: {dict(vc)}")

    # 指标①  真实边存活率 (真实图连接被 权重/GC/深度 剪掉多少)
    print(f"  指标① 边存活: {len(links_c)}/{len(links)} "
          f"({len(links_c) / len(links):.1%}) — GMW 断边阈值对真实流感图的温和性")

    # 指标② 短片段拓扑救回: unknown 段 → graph_viral
    short_unknown = [s for s in segs if len(segs[s][0]) < 1000]
    rescued = [s for s in short_unknown if verdicts.get(s) == gp.V_GRAPH_VIRAL]
    print(f"  指标② 短片段救回: {len(rescued)}/{len(short_unknown)} "
          f"({len(rescued) / max(1, len(short_unknown)):.1%}) 进入病毒目录")

    # 指标③ graph_viral 段 blastn 复核 (纯流感图: 应≈100% 确认)
    gv = [s for s, v in verdicts.items() if v == gp.V_GRAPH_VIRAL]
    conf = Counter(sub[s][0] for s in gv)
    n_none = conf.get("None", 0) + sum(1 for s in gv if sub[s][2])
    print(f"  指标③ graph_viral={len(gv)} 段的 blastn 复核: {dict(conf)} "
          f"(确认/歧义+无命中={n_none}) — 无命中段即拓扑证据超出序列比对的价值")

    # 指标④ 社区亚型纯度 (描述性: 拓扑是否自发分离 H3N2/H5N6)
    comm_sub = {}
    for s, c in comm.items():
        comm_sub.setdefault(c, []).append(sub[s][0])
    pure = mixed = 0
    for c, sts in comm_sub.items():
        sts = [x for x in sts if x != "None"]
        if len(set(sts)) > 1:
            mixed += 1
        else:
            pure += 1
    print(f"  指标④ 社区亚型纯度: 单亚型社区 {pure}, 跨亚型混合社区 {mixed} "
          f"(混合社区即保守段的重组区, GMW 论文用参考展开解决)")

    print(f"  demoted={len(demoted)}, 派生候选={Path(derived).name if derived else None}")
    shutil.rmtree(td, ignore_errors=True)
    return {"edge_survival": len(links_c) / len(links),
            "rescued": len(rescued), "rescued_total": len(short_unknown),
            "gv": len(gv), "confirmed": conf.get("H3N2", 0) + conf.get("H5N6", 0),
            "pure": pure, "mixed": mixed}


def main():
    ap = argparse.ArgumentParser(description="graph 阶段独立运行验证 (无外部映射工具)")
    ap.add_argument('--gmw-repo',
                    default=r"D:\桌面\文章写作\参考文献\工具、方法与算法\gmw")
    ap.add_argument('--skip-blast', action='store_true')
    args = ap.parse_args()

    ok_a = scenario_a()
    m = scenario_b(args.gmw_repo, use_blast=not args.skip_blast)
    print("=" * 70)
    print(f"结论: 场景A(合成ground truth)={'通过' if ok_a else '失败'}; "
          f"场景B(GMW真实组装图): 边存活 {m['edge_survival']:.1%}, "
          f"短片段救回 {m['rescued']}/{m['rescued_total']}, "
          f"graph_viral {m['gv']} 段 (blastn 确认 {m['confirmed']}), "
          f"单亚型社区 {m['pure']}/混合 {m['mixed']}")
    print("独立运行方式 (服务器, 有 bwa-mem2/samtools/coverm):")
    print("  python graph_propagation.py --asm-dir <d_asm> --ident-dir <d_ident> "
          "--reads-dir <d_hostdep> --output-dir <d_graph> -t 16 -j 4")
    return 0 if ok_a else 1


if __name__ == '__main__':
    sys.exit(main())
