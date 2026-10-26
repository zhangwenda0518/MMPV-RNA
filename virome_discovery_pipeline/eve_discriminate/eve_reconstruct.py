#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""eve_reconstruct.py — 位点蛋白重建 (两级: miniprot 定位 → GeneWise 计数)

为什么需要这一层 (v7.2, 审稿意见): s1 的退化判定直接数"提前终止密码子", 但
终止可能来自**组装 indel** 而不是真实退化 —— 一个单碱基 indel 就让其后全部密码子
错读。hi-fever 用 GeneWise 在核苷酸层面重建蛋白、显式建模移码与内含子, 才能把
"真退化"和"组装错误"分开。本模块补上这一层。

**v7.3 改为两级架构** (借鉴 EasyPseudogene, MIT):
  1. **miniprot** 快速全基因组/位点定位 (PAF, CIGAR 的 F/N 事件) —— 便宜, 用于
     筛出"这个位点有没有同源蛋白、有没有移码/内含子";
  2. **GeneWise** 对 miniprot 给出的候选做精确重建, 在**重建后的蛋白比对串**上
     数终止 (X) 与移码 (!) —— 贵, 但计数语义可靠。
为什么必须两级: 在原始核苷酸序列上数 TAA/TAG/TGA 会被组装 indel 污染 (错读框),
而 GeneWise 的比对串里每个 X 都已核对过对应密码子 (见 `verify_gw_stops`)。
miniprot 单独用于"判定"会漏掉 GeneWise 才能看见的逐密码子事件。

后端选择:
  * `miniprot` 必有 (定位层); `genewise` 可选 (计数层)。
  * genewise 缺失 -> 用 miniprot 的 CIGAR 事件计数, 并在 status 列标
    `miniprot_only` (证据等级较低, 不伪造 GeneWise 语义)。
  * 两者都没有 -> `tool_missing`, 不伪造结果。

关于 miniprot 能否替代 GeneWise:
  * 能替代: 定位、移码/内含子事件的有无、比对质量。
  * 不能完全替代: 逐密码子的终止核验 (GeneWise 会把 X 与显示出的密码子对照)。
    因此**计数层优先 GeneWise**, miniprot 计数只作降级路径。

输出: <outdir>/05_Reconstruct/<NAME>.recon.tsv (逐位点)
  locus  backend  n_frameshift  n_stop  n_intron  aln_aa  identity  coverage
  earliest_event_frac  local_flank_identity  status
  * `earliest_event_frac`: 最早事件在比对查询上的相对位置 (用于"晚期突变"过滤)
  * `local_flank_identity`: 单事件候选的事件两侧局部一致性 (用于判别测序错误)
  * status: ok / no_hit / miniprot_only / tool_missing / error

用法:
  python eve_reconstruct.py -g genome.fa -l loci.bed --loci-fa loci.fa \
      -p proteins.fa -o out/ -n NAME [--backend auto|miniprot|genewise] [-t 16]
"""
import argparse
import os
import re
import shutil
import subprocess
import sys
import tempfile
from bisect import bisect_left
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
# eve_scan_core 在 01_scan/；脚本单独运行时也要能找到
for _p in (SCRIPT_DIR, SCRIPT_DIR.parent / "01_scan"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

from eve_scan_core import safe_dir  # noqa: E402

RECON_DIR = "05_Reconstruct"
RECON_HDR = ("locus\tbackend\tn_frameshift\tn_stop\tn_intron\taln_aa\t"
             "identity\tcoverage\tearliest_event_frac\tlocal_flank_identity\tstatus\n")

# miniprot PAF 的 CIGAR 算子: M/I/D 为比对, N 为内含子, F 为移码 (--fs 时出现)
CIGAR_INTRON = "N"
CIGAR_FRAMESHIFT = "F"

# GeneWise pretty 输出的行结构 (与 EasyPseudogene 同一解析口径):
#   <label> <coord> <query aa>
#   ... (2 行)
#   <genomic DNA 三行>
STOPS = {"TAA", "TAG", "TGA"}
BASES = set("ACGTN")
GW_ROW = re.compile(r'^\S+\s+-?\d+\s+(\S.*)$')


# ---------------- 纯逻辑 (可单测) ----------------
def parse_cigar(cigar):
    """CIGAR -> (n_intron, n_frameshift, aligned_query_bp).

    miniprot PAF 的 cg:Z: 标签用 `[0-9]+[MIDNSHF]` 的算子串, 其中
      N = 内含子 (跳过 reference), F = 移码 (蛋白-基因组比对特有),
      M/I/D = 普通比对, S = soft clip, H = hard clip。
    返回 (内含子事件数, 移码事件数, 比对上的 query 碱基数)。
    解析失败返回 (0, 0, 0) —— 调用方据 status 列区分"真 0"与"没解析出来"。
    """
    if not cigar:
        return 0, 0, 0
    ops = re.findall(r"(\d+)([MIDNSHF])", cigar)
    if not ops:
        return 0, 0, 0
    # 事件数按**算子出现次数**计, 不按长度累加: "500N" 是 1 个内含子 (长 500nt),
    # "2F" 是 1 次移码 (2nt), 判定关心的是"发生了几次", 不是"跨了多少碱基"。
    n_intron = sum(1 for _n, op in ops if op == CIGAR_INTRON)
    n_fs = sum(1 for _n, op in ops if op == CIGAR_FRAMESHIFT)
    aln_bp = sum(int(n) for n, op in ops if op in ("M", "I"))
    return n_intron, n_fs, aln_bp


def parse_miniprot_paf(lines):
    """miniprot --gff 的 ##PAF 行 -> {query: 最佳比对 dict}.

    真实格式 (miniprot 0.18 实测):
      `##PAF\t<q>\t<qlen>\t<qstart>\t<qend>\t<strand>\t<t>\t<tlen>\t<tstart>\t<tend>\
       <n_match>\t<blen>\t<mapq>\t<tags...>`
    注意它**以 `##PAF` 开头** (GFF3 的注释标记), 但内容是标准 PAF 列 ——
    早期实现把所有 `#` 开头的行都跳过, 结果真实输出一条都解析不出来
    (线上表现为"全 no_hit", 而 miniprot 其实有命中)。
    这里只跳过 `##gff-version` / 普通注释, 保留 `##PAF`。

    取每个 query 的**最高比对得分**那条 (AS:i: 标签; 没有就按比对长度)。
    """
    best = {}
    for line in lines:
        if not line.strip():
            continue
        s = line.rstrip("\n")
        if s.startswith("##PAF"):
            # 去掉前缀标记, 后面的列与标准 PAF 一致
            s = s[len("##PAF"):].lstrip("\t")
        elif s.startswith("#"):
            continue
        c = s.split("\t")
        if len(c) < 12:
            continue
        q = c[0]
        qlen = int(c[1]) if c[1].isdigit() else 0
        aln_bp = int(c[3]) if c[3].isdigit() else 0
        # PAF 列 (0 基): 0 qname 1 qlen 2 qstart 3 qend 4 strand 5 tname 6 tstart
        # 7 tend 8 n_match 9 aln_block_len 10 mapq ... (不是 9/10, 是 8/9)
        nmatch = int(c[8]) if c[8].isdigit() else 0
        blen = int(c[9]) if c[9].isdigit() else 0
        tags = {}
        for t in c[12:]:
            parts = t.split(":", 2)
            if len(parts) == 3:
                tags[parts[0]] = parts[2]      # AS:i:300 -> tags["AS"]="300"
        try:
            score = int(tags.get("AS", aln_bp) or 0)
        except ValueError:
            score = aln_bp
        cigar = tags.get("cg", "")
        n_intron, n_fs, cg_bp = parse_cigar(cigar)
        ident = round(nmatch / blen, 4) if blen else 0.0
        rec = {"n_fs": n_fs, "n_intron": n_intron,
               "aln_bp": cg_bp or aln_bp, "identity": ident,
               "aln_aa": max(0, (cg_bp or aln_bp) // 3),
               "target": c[5], "qlen": qlen, "_score": score}
        if q not in best or score > best[q]["_score"]:
            best[q] = rec
    return best


def count_stops_in_frame(seq, frame=0):
    """核苷酸序列按读框数终止密码子 (重建后的对照口径, 与 s1 独立)。"""
    stops = {"TAA", "TAG", "TGA"}
    n = 0
    for i in range(frame, len(seq) - 2, 3):
        if seq[i:i + 3].upper() in stops:
            n += 1
    return n


# ---------------- GeneWise 解析 (移植 EasyPseudogene, MIT) ----------------
def gw_segments(block):
    """GeneWise pretty 输出块 -> [(query_aa, target_verified), ...].

    关键点 (与 EasyPseudogene 一致): GeneWise 的 `X` 可能来自**模糊密码子**
    而不是真终止, 必须用同列显示的三行 DNA 核对 —— 只有密码子确实是
    TAA/TAG/TGA 才记 `X`, 否则记 `?` (不当作终止证据)。
    行结构: `<label> <coord> <query>` 之后隔 2 行是三行 genomic DNA (固定列对齐)。
    """
    lines = block.splitlines()
    for i, line in enumerate(lines[:-5]):
        m = GW_ROW.match(line)
        if not m:
            continue
        offset = m.start(1)
        query = m.group(1).rstrip()
        dna_rows = lines[i + 3:i + 6]
        genomic_row = dna_rows[0]
        if not GW_ROW.match(genomic_row):
            prefix = genomic_row[:offset].strip()
            if ":" not in prefix or not any(
                    all(len(row) > offset + j and row[offset + j].upper() in BASES
                        for row in dna_rows)
                    for j in range(len(query))):
                continue
        target = lines[i + 2][offset:offset + len(query)].ljust(len(query))
        bases = [s.ljust(offset + len(query)) for s in dna_rows]
        verified = []
        for j, aa in enumerate(target):
            codon = "".join(s[offset + j] for s in bases).upper()
            if aa in "X*":
                verified.append("X" if codon in STOPS else "?")
            else:
                verified.append(aa)
        yield query, "".join(verified)


def gw_mutation_counts(block):
    """GeneWise 输出块 -> (n_stop, n_frameshift).

    `!` 是 GeneWise 的移码标记; `X` 只有经过密码子核对才保留 (见 gw_segments)。
    没有 pretty 对齐时退回紧凑输出的保守计数: **不把未核对的 X 当终止**
    (`max(count('*'), count(':X['))` 与 EasyPseudogene 同口径)。
    """
    alns = list(gw_segments(block))
    if alns:
        target = "".join(t for _, t in alns)
        return target.count("X"), target.count("!")
    return max(block.count("*"), block.count(":X[")), block.count("!")


# ---------------- 比对度量 (移植 EasyPseudogene 语义) ----------------
def query_positions(query, target):
    """query 串里真实氨基酸的位置 (跳过内含子空格与 gap 的 '-')。"""
    return [i for i, aa in enumerate(query[:len(target)]) if aa.isalpha()]


def aligned_query_length(query, target):
    return len(query_positions(query, target))


def sequence_identity(query, target, positions=None):
    """按 query 氨基酸位置算一致性 (空格/gap 不计入, 避免虚高)。"""
    if positions is None:
        positions = query_positions(query, target)
    if not positions:
        return 0.0
    return sum(query[i].upper() == target[i].upper() for i in positions) / len(positions)


def local_flank_identity(query, target, column, flank=10):
    """事件位置两侧 (共 2*flank 个 query 氨基酸) 的局部一致性。

    单事件候选的判别依据: 两侧一致性低说明该位置本来就对不齐 (更像测序/组装
    错误), 高说明突变落在保守区 (更像真退化)。移植自 EasyPseudogene。
    """
    positions = query_positions(query, target)
    if not positions:
        return 0.0
    off = bisect_left(positions, column)
    after = off + int(off < len(positions) and positions[off] == column)
    local = positions[max(0, off - flank):off] + positions[after:after + flank]
    return sequence_identity(query, target, local)


def mutation_query_fraction(query, target, column):
    """事件在 query 比对上的相对位置 (0=最前, 1=最后) 与 (covered, total)."""
    positions = query_positions(query, target)
    if not positions:
        return 0.0, 0, 0
    pre = bisect_left(positions, column)
    covered = pre + int(pre < len(positions) and positions[pre] == column)
    return covered / len(positions), covered, len(positions)


def gw_metrics(block, flank=None):
    """GeneWise 块 -> 判定所需的全部度量.

    返回 dict: n_stop / n_frameshift / aln_aa / identity / coverage /
    earliest_event_frac / local_flank_identity / n_events。
    `earliest_event_frac` 与 `local_flank_identity` 分别喂给"晚期突变"与
    "单事件低局部一致性"两个过滤器。
    """
    if flank is None:
        flank = SINGLE_EVENT_FLANK
    alns = list(gw_segments(block))
    if not alns:
        return {"n_stop": 0, "n_frameshift": 0, "aln_aa": 0, "identity": 0.0,
                "coverage": 0.0, "earliest_event_frac": "", 
                "local_flank_identity": "", "n_events": 0}
    # 取最长的那段比对做度量 (多段时按最长代表)
    query, target = max(alns, key=lambda qt: len(qt[0]))
    pos = query_positions(query, target)
    n_stop = target.count("X")
    n_fs = target.count("!")
    events = [i for i, ch in enumerate(target) if ch in ("X", "!")]
    if events:
        first = events[0]
        frac, _c, _t = mutation_query_fraction(query, target, first)
        lfi = local_flank_identity(query, target, first, flank)
    else:
        frac, lfi = "", ""
    return {"n_stop": n_stop, "n_frameshift": n_fs,
            "aln_aa": len(pos), "identity": round(sequence_identity(query, target), 4),
            "coverage": round(len(pos) / max(len(query), 1), 4),
            "earliest_event_frac": round(frac, 4) if frac != "" else "",
            "local_flank_identity": round(lfi, 4) if lfi != "" else "",
            "n_events": n_stop + n_fs}


# ---------------- 四个判定过滤器 (移植 EasyPseudogene 阈值) ----------------
LATE_EVENT_FRAC = 0.9        # 最早事件 >= 0.9 的比对位置 -> 太靠末端, 不可信
SINGLE_FLANK_IDENT = 0.4     # 单事件时两侧局部一致性下限
SINGLE_EVENT_FLANK = 10      # 单事件局部一致性窗口 (aa, 两侧各取; --flank 可改)
MULTI_MIN_COV = 0.7          # 多事件质量闸门
MULTI_MIN_IDENT = 0.05
MULTI_MIN_ALN_AA = 200


def filter_late_event(m):
    """晚期突变过滤: 唯一/最早事件落在比对末端 10% -> 排除。"""
    f = m.get("earliest_event_frac")
    if f == "" or f is None:
        return False, ""
    if f >= LATE_EVENT_FRAC:
        return True, f"最早事件位于比对 {f:.1%} (>= {LATE_EVENT_FRAC:.0%}), 太靠末端"
    return False, ""


def filter_single_event_flank(m):
    """单事件低局部一致性过滤: 只 1 个事件且两侧一致性 < 0.4 -> 疑测序/组装错误。"""
    if m.get("n_events", 0) != 1:
        return False, ""
    lfi = m.get("local_flank_identity")
    if lfi == "" or lfi is None:
        return False, ""
    if lfi < SINGLE_FLANK_IDENT:
        return True, f"单事件但两侧局部一致性仅 {lfi:.2f} (< {SINGLE_FLANK_IDENT}), 疑测序/组装错误"
    return False, ""


def filter_multi_quality(m):
    """多事件质量闸门: 事件 >=2 但覆盖/一致性/比对长度不达标 -> 证据不足。"""
    if m.get("n_events", 0) < 2:
        return False, ""
    bad = []
    if m.get("coverage", 0) < MULTI_MIN_COV:
        bad.append(f"覆盖 {m.get('coverage',0):.2f} < {MULTI_MIN_COV}")
    if m.get("identity", 0) < MULTI_MIN_IDENT:
        bad.append(f"一致性 {m.get('identity',0):.2f} < {MULTI_MIN_IDENT}")
    if m.get("aln_aa", 0) < MULTI_MIN_ALN_AA:
        bad.append(f"比对 {m.get('aln_aa',0)}aa < {MULTI_MIN_ALN_AA}")
    if bad:
        return True, "多事件质量不足: " + "; ".join(bad)
    return False, ""


def filter_same_locus_conflict(mp_rec, gw_m):
    """同一位点冲突: miniprot 报告零 stop/fs 而 GeneWise 报有事件 -> 存疑。

    EasyPseudogene 的 Cross 模式规则: GeneWise-only 的突变若被一条"完整且
    零 stop/fs 标签"的 miniprot 命中反驳, 则排除。这里降级为 REVIEW (不排除),
    因为 EVE 场景下两者差异更可能是参数/边界效应。
    """
    if not mp_rec or not gw_m:
        return False, ""
    if (mp_rec.get("n_fs", 0) == 0 and mp_rec.get("n_intron", 0) == 0
            and gw_m.get("n_events", 0) > 0):
        return True, ("miniprot 报告零移码/零内含子, GeneWise 报告 "
                      f"{gw_m.get('n_events')} 个事件 -> 两通道冲突, 转人工")
    return False, ""


def apply_filters(mp_rec, gw_m):
    """依次跑 4 个过滤器 -> (excluded: bool, reasons: list[str]).

    顺序按 EasyPseudogene: 晚期突变 → 单事件局部 → 多事件质量 → 同位点冲突。
    返回 excluded 表示"该位点不应作为退化证据"; 冲突类只标 reason 不排除。
    """
    reasons = []
    excluded = False
    for fn in (filter_late_event, filter_single_event_flank, filter_multi_quality):
        hit, why = fn(gw_m)
        if hit:
            excluded = True
            reasons.append(why)
    hit, why = filter_same_locus_conflict(mp_rec, gw_m)
    if hit:
        reasons.append(why)
    return excluded, reasons


# ---------------- 后端探测 ----------------
def pick_backend(requested="auto"):
    """探测两级后端 -> (locator, counter).

    locator: 定位层 (miniprot) 可执行文件或 ""
    counter: 计数层 (genewise) 可执行文件或 ""
    requested 只影响是否强制要求某一层:
      auto      -> 有谁用谁 (miniprot 定位 + genewise 计数, 缺 genewise 降级)
      miniprot  -> 只要求 miniprot
      genewise  -> 只要求 genewise (此时也用它做定位)
    """
    mp = shutil.which("miniprot") or ""
    gw = shutil.which("genewise") or ""
    if requested == "miniprot":
        return (mp, "") if mp else ("", "")
    if requested == "genewise":
        return (gw, gw) if gw else ("", "")
    return mp, gw


def miniprot_cmd(exe, genome, proteins, out_paf, threads, fs=True):
    """[参考实现, 当前主流程未调用] 全基因组 miniprot 命令构造。

    保留原因: 这里记录了两个**真实数据实测踩过的坑**, 将来若要做全基因组级
    定位 (而不是逐位点) 会直接用到:
      * 绝对路径时 miniprot 剥扩展名找 `.mpi`, 而 `-d` 建的是 `<genome>.mpi`
        —— 索引明明在也报 "failed to open/build the index"; 用相对名 + chdir。
      * **`-I` 是布尔开关不吃数值**; 写 `-I 200k` 会让 `200k` 变成位置参数,
        报的还是同一个索引错误 (极难定位)。要指定上限用 `-G 200k`。
    当前主流程改为**逐位点** (位点作 ref, 蛋白作 query), 因为全基因组模式下
    PAF 的 query 是蛋白名, 无法映射回位点。
    """
    g = Path(genome).name
    cmd = [exe, "--gff", "--trans", "--aln", "-t", str(max(1, int(threads))),
           "-G", "200k", g, str(Path(proteins).resolve())]
    return cmd, out_paf


def build_index(exe, genome, idx=None, timeout=0):
    """[参考实现, 当前主流程未调用] 预建 miniprot 索引 (`<genome>.mpi`).

    真实数据实测: 1.5G 基因组建索引约 35s / 12.4GB RSS; 目录不可写时抛错
    (不静默降级成"没比对上")。
    """
    g = Path(genome)
    target = Path(idx) if idx else Path(str(g) + ".mpi")
    if not os.access(str(target.parent), os.W_OK):
        raise RuntimeError(f"基因组目录不可写, 无法建 miniprot 索引: {target.parent}"
                           f" (可把基因组软链到可写目录后重试)")
    p = subprocess.run([exe, "-d", str(target), str(g)],
                       capture_output=True, text=True, timeout=timeout or None)
    if p.returncode != 0:
        raise RuntimeError(f"miniprot 建索引失败 (exit {p.returncode}): "
                           f"{(p.stderr or '')[-300:]}")
    return target


def genewise_cmd(exe, protein_fa, locus_fa, block_out, wisecfg_dir=""):
    """GeneWise 单对命令 (与 hi-fever / EasyPseudogene 同参)。

    `-pretty` 输出可解析的逐列比对 (计数层要的就是它); `-both` 双向;
    `-kbyte 2000000` 给大位点留内存; `-silent` 抑制进度。
    wisecfg 目录通过环境变量 WISECONFIGDIR 传入 (GeneWise 自己读它)。
    """
    return [exe, protein_fa, locus_fa, "-both", "-kbyte", "2000000",
            "-matrix", "BLOSUM62.bla", "-pretty", "-silent"]


def run_genewise_on_loci(gw_exe, loci, best_prot, workdir, timeout=0):
    """对每个位点跑 GeneWise -> {locus: 度量 dict}.

    `best_prot`: {locus: 蛋白序列} —— 由 miniprot 定位阶段选出该位点最佳的
    参考蛋白。GeneWise 一次只能比一对, 所以逐位点调用。
    返回的度量来自 `gw_metrics` (含 n_stop/n_frameshift/earliest_event_frac/
    local_flank_identity)。
    """
    out = {}
    for q, seq in loci.items():
        prot = best_prot.get(q)
        if not prot:
            continue
        pf = Path(workdir) / "gw_prot.fa"
        qf = Path(workdir) / "gw_locus.fa"
        pf.write_text(">%s\n%s\n" % (q, prot), encoding="utf-8")
        qf.write_text(">%s\n%s\n" % (q, seq), encoding="utf-8")
        try:
            p = subprocess.run(genewise_cmd(gw_exe, str(pf), str(qf), ""),
                               capture_output=True, text=True,
                               timeout=timeout or None)
        except subprocess.TimeoutExpired:
            continue
        block = (p.stdout or "") + "\n" + (p.stderr or "")
        if "//" not in block and not gw_segments(block):
            continue
        # 多个 // 块: 取第一个 (单对单)
        out[q] = gw_metrics(block.split("//")[0])
    return out


# ---------------- 主流程 (两级) ----------------
def reconstruct(genome, loci_fa, proteins, outdir, name, threads=16,
                backend="auto", timeout=0):
    """两级重建: miniprot 定位 -> (可选) GeneWise 计数.

    返回 (输出路径, 逐位点行数, backend_label)。
    backend_label: `miniprot+genewise` / `miniprot_only` / `genewise` / `unavailable`

    **关键设计** (真实数据实测修正): 不能拿"整个基因组"当 ref、蛋白当 query ——
    那样 miniprot 的 PAF query 是**蛋白名**, 无法映射回位点。正确做法是
    **逐位点**: 位点序列当 ref, 参考蛋白当 query, miniprot 报告的 query 就是
    蛋白名, 而 ref 只有一条 (= 当前位点), 天然一一对应。这也让 GeneWise 的
    第二级只需处理"该位点有蛋白命中"的情形。
    """
    out_d = Path(outdir) / RECON_DIR
    out_d.mkdir(parents=True, exist_ok=True)
    out_tsv = out_d / f"{name}.recon.tsv"
    mp_exe, gw_exe = pick_backend(backend)
    loci = read_fasta(loci_fa)
    prot_seqs = (proteins if isinstance(proteins, dict) else read_fasta(proteins))

    if not mp_exe and not gw_exe:
        rows = [(q, "unavailable", "", "", "", "", "", "", "", "", "tool_missing")
                for q in loci]
        write_rows(out_tsv, rows)
        return out_tsv, len(rows), "unavailable"

    label = ("miniprot+genewise" if mp_exe and gw_exe
             else "miniprot_only" if mp_exe else "genewise")
    rows = []
    with tempfile.TemporaryDirectory(prefix="eve_recon_") as td:
        td = Path(td)
        for q, seq in loci.items():
            # 该位点候选: 位点序列作 ref (单条), 全部参考蛋白作 query
            ref = td / "locus.fa"
            ref.write_text(">%s\n%s\n" % (q, seq), encoding="utf-8")
            prot_fa = td / "proteins.fa"
            if not prot_fa.exists():
                with open(prot_fa, "w", encoding="utf-8") as fo:
                    for pn, ps in prot_seqs.items():
                        fo.write(">%s\n%s\n" % (pn, ps))
            mp_rec = {}
            if mp_exe:
                out_paf = td / "mp.gff"
                idx = Path(str(ref) + ".mpi")
                # 位点很小, 索引瞬时完成; 索引名与 ref 全名一致 (miniprot 的查找规则)
                p = subprocess.run([mp_exe, "-d", str(idx), str(ref)],
                                   capture_output=True, text=True,
                                   timeout=timeout or None)
                if p.returncode != 0:
                    rows.append((q, label, "", "", "", "", "", "", "", "", "error"))
                    continue
                cmd = [mp_exe, "--gff", "--trans", "--aln", "-t",
                       str(max(1, int(threads))), "-G", "200k",
                       ref.name, str(prot_fa.resolve())]
                with open(out_paf, "w", encoding="utf-8") as fo:
                    p = subprocess.run([str(c) for c in cmd], stdout=fo,
                                       stderr=subprocess.PIPE, text=True,
                                       cwd=str(td), timeout=timeout or None)
                if p.returncode != 0:
                    rows.append((q, label, "", "", "", "", "", "", "", "", "error"))
                    continue
                with open(out_paf, encoding="utf-8", errors="replace") as fh:
                    best = parse_miniprot_paf(fh)
                if not best:
                    rows.append((q, label, "", "", "", "", "", "", "", "", "no_hit"))
                    continue
                # 取得分最高的那个蛋白
                pname = max(best, key=lambda k: best[k]["_score"])
                mp_rec = best[pname]
                mp_rec["target"] = pname

            # 第二级: GeneWise 精确计数 (该位点 + 该蛋白)
            gw_m = {}
            if gw_exe:
                pname = mp_rec.get("target") if mp_rec else None
                if pname and pname in prot_seqs:
                    gw_m = run_genewise_on_loci(
                        gw_exe, {q: seq}, {q: prot_seqs[pname]}, td, timeout)

            g = gw_m.get(q)
            if g:
                excluded, reasons = apply_filters(mp_rec, g)
                status = "filtered" if excluded else "ok"
                rows.append((q, label, g["n_frameshift"], g["n_stop"],
                             mp_rec.get("n_intron", ""), g["aln_aa"],
                             g["identity"], g["coverage"],
                             g["earliest_event_frac"],
                             g["local_flank_identity"], status))
            elif mp_rec:
                cov = (round(mp_rec["aln_bp"] / mp_rec["qlen"], 4)
                       if mp_rec.get("qlen") else "")
                rows.append((q, label, mp_rec["n_fs"], "", mp_rec["n_intron"],
                             mp_rec["aln_aa"], mp_rec["identity"], cov,
                             "", "", "miniprot_only"))
            else:
                rows.append((q, label, "", "", "", "", "", "", "", "", "no_hit"))
    write_rows(out_tsv, rows)
    return out_tsv, len(rows), label


def read_fasta(path):
    seqs, name, buf = {}, None, []
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buf)
                name = line[1:].strip().split()[0] if line[1:].strip() else ""
                buf = []
            elif name is not None:
                buf.append(line.strip())
    if name:
        seqs[name] = "".join(buf)
    return seqs


def write_rows(path, rows):
    with open(path, "w", encoding="utf-8", newline="\n") as fo:
        fo.write(RECON_HDR)
        for r in rows:
            fo.write("\t".join(str(x) for x in r) + "\n")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-g", "--genome", required=True, help="宿主基因组 FASTA")
    ap.add_argument("-l", "--loci-bed", required=True, help="Stage1 loci.bed")
    ap.add_argument("--loci-fa", required=True, help="Stage1 loci.fa (位点序列)")
    ap.add_argument("-p", "--proteins", required=True,
                    help="参考蛋白 FASTA (如 Stage1 最佳命中的 ref 序列)")
    ap.add_argument("-o", "--outdir", required=True, help="eve_screen 输出根")
    ap.add_argument("-n", "--name", required=True, help="基因组名")
    ap.add_argument("-t", "--threads", type=int, default=16)
    ap.add_argument("--backend", default="auto",
                    choices=["auto", "miniprot", "genewise"])
    ap.add_argument("--cmd-timeout", type=int, default=0)
    ap.add_argument("--flank", type=int, default=10,
                    help="单事件局部一致性窗口 (aa, 两侧各取这么多; 默认 10)")
    a = ap.parse_args(argv)
    global SINGLE_EVENT_FLANK
    SINGLE_EVENT_FLANK = a.flank
    out, n, b = reconstruct(a.genome, a.loci_fa, a.proteins,
                            safe_dir(a.outdir), a.name, a.threads, a.backend,
                            a.cmd_timeout)
    print(f"[done] reconstruct({b}): {n} 位点 -> {out}")
    if b == "unavailable":
        print("[warn] 未找到 miniprot/genewise: 本层记为 tool_missing, "
              "退化判定退回 s1 的终止计数 (证据强度较低)", file=sys.stderr)
        return 3
    return 0


if __name__ == "__main__":
    sys.exit(main())
