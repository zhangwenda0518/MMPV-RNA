#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
eve_element_db.py — 元件级 (element-level) EVE 后处理层
========================================================
eve_screen.py 收敛到的是**位点**层 (viral_supported / host_like / undetermined);
本层把 viral_supported 位点进一步收敛成**元件**层 (element): 跨基因组同源聚类 →
元件区归并 → consensus 文库 → 拷贝数 / 结构域 / 系统发育输入。

输入契约 (一次跑完的 eve_screen 输出根目录 <outdir>):
  04_Summary/<NAME>_eve_summary.tsv   位点判定汇总 (12 列, 不含 genome 列)
  01_Loci/<NAME>/<NAME>.loci.fa       位点序列 (Stage1 抽取, 头即 locus 串)
  -B/--batch 的 name<TAB>基因组路径    rm 阶段回比用 (与 eve_screen.py -B 同一张表)
输出契约: <outdir>/05_Elements/ —— 文件名沿用源脚本 post/ 下的名字, 下游
  build_eve_db.py 认的 candidates_meta.tsv / cluster_members.tsv / all_candidates.fa
  数据契约原样保留 (列名与列序见各 HDR 常量)。

为什么用 04_Summary 而不是 02_Verdict / 03_RVDB: summary 是唯一把三张表 join 好的单表 ——
Stage1 的 ref_family/ref_sseqid/ref_bitscore、Stage2 的 verdict、Stage3 的
rvdb_stitle/rvdb_bitscore 都在。源脚本的 evidence_bs 取 max(ref_bitscore,
rvdb_bitscore) 两条证据通道的强者, ref_family 又要喂给聚类的家族投票与 phylo 分组;
只读 02_Verdict 会丢 ref_family/rvdb_bitscore, 只读 03_RVDB 会丢 verdict, 且
03 只收 viral_supported + undetermined 两类候选 (host_like 位点根本不在里面)。

locus 与 loci.fa 头两侧都过 locus_key() 再 join: 头是抽序列工具按区域串写的
(chr1:START-END), 历史工具会把它改写成 chr1_START-END —— 哪一侧变了 join 都静默落空
(位点序列全丢、元数据表照样有行、退出码仍是 0), 与 Stage2/3 是同一个坑。

阶段 (断点语义: **该阶段声明的产物文件全部存在即视为完成**, --force 重跑; 不用
  .done 标记 —— 与 eve_screen.py / eve_genome_scan.py 的产物判定保持一致):
  collect    viral_supported 位点序列收集 (加 genome 前缀, 防跨基因组 ID 撞名)
  cluster    MMseqs2 easy-linclust 跨基因组聚类 → 元件主表 (含嵌合体初筛)
  regions    ±5kb 同家族位点归并 → 元件区表 (EPRV 片段化计数校正, Caulifinder 思路)
  consensus  元件 consensus 文库 (默认簇内最长成员; --mafft-consensus 走比对共识)
  rm         RepeatMasker 用 consensus 文库回比各基因组 → 每家族拷贝数/覆盖度
  cdd        mmseqs 翻译搜索 vs CDD profile 库 → 结构域注释 (需 --cdd-db)
  phylo      系统发育输入准备 (代表序列 + 参考集; 命令脚本 phylo_commands.sh)

用法:
  # 只跑不需要额外参考库的四个阶段 (默认)
  python eve_element_db.py -o eve_results --stages collect,cluster,regions,consensus
  # 全量
  python eve_element_db.py -o eve_results --stages all -t 40 -J 8 -B batch.tsv \
      --cdd-db ~/database/cdd/cdd-db/cdd_db --rt-ref RT_refs.fa --run-phylo
  # 断点续跑: 重发同命令 (已完成阶段按产物跳过), 要重跑就加 --force

需要额外依赖的阶段: cdd 要 --cdd-db (mmseqs profile 库前缀), rm 要 -B 与
RepeatMasker (可用 RM_BIN 环境变量或 --rm-bin 指定), phylo 的 --run-phylo 还要
--rt-ref 与 mafft/trimal/iqtree。缺什么开跑前一次性报出来, 不做静默降级。

来源: 移植自 hi-fever/eve_kingdom_post.py (2026-09-19 定版), 算法参数逐项保留:
  easy-linclust --min-seq-id 0.8 -c 0.8 --cov-mode 1; consensus 默认取簇内最长成员,
  --mafft-consensus 走 mafft --auto --anysymbol --thread 1 + 0.4 频率列投票;
  RepeatMasker -no_is -norna -pa <threads>; CDD 走 mmseqs search --search-type 2
  (原 rpstblastn 读不了 msa2profile 建的库)。逐项偏差见函数注释与 doc.md 的移植说明。
"""

import argparse
import os
import shutil
import sys
import time
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))

from eve_genome_scan import setup_logger                    # noqa: E402
from eve_scan_core import (SAFE_CHARS, STAGE_DIRS, clean_name, logmsg,  # noqa: E402
                           locus_key, run, safe_dir, tool_path)

# ---------------- 路径与默认参数 ----------------
ELEM_DIR = "05_Elements"          # 元件层产物目录 (接在 01_Loci/02_Verdict/03_RVDB/04_Summary 之后)
SUMMARY_SUFFIX = "_eve_summary.tsv"
LOCI_SUFFIX = ".loci.fa"

DEFAULT_STAGES = ("collect", "cluster", "regions", "consensus")
ALL_STAGES = ("collect", "cluster", "regions", "consensus", "rm", "cdd", "phylo")
# 数字别名按阶段顺序 (与 eve_genome_scan.parse_stages 的 1/2/3 同风格, 这里 1..7)
STAGE_ALIAS = {"1": "collect", "2": "cluster", "3": "regions", "4": "consensus",
               "5": "rm", "6": "cdd", "7": "phylo"}

# 阶段 -> 完成判定用的产物 (全部存在才算该阶段完成). 与 STAGE_OUTPUTS 逐条对照的用例
# 在 tests/test_eve_element_db.py, 阶段一改名这里就必须跟着改.
STAGE_OUTPUTS = {
    "collect":   ("all_candidates.fa", "candidates_meta.tsv"),
    "cluster":   ("element_table.tsv",),
    "regions":   ("element_regions.tsv",),
    "consensus": ("eve_consensus_library.fa",),
    "rm":        ("family_quant.tsv",),
    "cdd":       ("consensus_domains.tsv",),
    "phylo":     ("phylo_commands.sh",),
}

# 表头常量: 与写盘处共用同一份字面量, 测试拿它跟真实产物第一行逐列比.
# 表头是跨脚本契约 (下游 build_eve_db.py 按列名读), 在别处再抄一份就会悄悄错位.
CAND_META_HDR = "genome\tlocus\tfamily\tverdict\tref_sseqid\tevidence_bs"
ELEMENT_TABLE_HDR = ("cluster_id\trep\tn_members\tn_genomes\tfamily\t"
                     "chimera_flag\tbest_member\tbest_bs")
ELEMENT_REGIONS_HDR = "genome\tregion_id\tcontig\tstart\tend\tn_families\tfamilies"
FAMILY_QUANT_HDR = "family\tgenome\tcopies\tcovered_bp"
# mmseqs convertalis 不写表头, 列序就写在这里 (下游按位置解析这 9 列)
CDD_CONVERT_FIELDS = "query,target,evalue,bits,qstart,qend,tstart,tend,tlen"
PHYLO_SCRIPT = "phylo_commands.sh"

# 算法参数 (与源脚本逐项一致; 不做重新调参)
MIN_BS = 50.0             # collect 的 evidence_bs 下限 (原 --min-bs)
CLUSTER_MIN_SEQ_ID = "0.8"   # easy-linclust --min-seq-id
CLUSTER_COV = "0.8"          # easy-linclust -c
CLUSTER_COV_MODE = "1"       # easy-linclust --cov-mode (覆盖 query)
ER_DIST = 5000            # regions 的元件区归并距离 (原 --er-dist, ±5kb)
MAFFT_CONS_MAX_N = 50     # --mafft-consensus 时参与比对的成员数上限 (源脚本硬编码)
CDD_SENS = 4.0            # mmseqs search -s
CDD_MAX_SEQS = 300        # mmseqs search --max-seqs
MIN_CLUSTERS = 10         # phylo: 家族至少要有这么多个簇才出代表集

REPEATMASKER_DEFAULT = os.environ.get("RM_BIN", "RepeatMasker")


@dataclass
class ElementConfig:
    """元件层全部外部依赖与运行参数.

    与 EveConfig 分开放: 两者面向的工具链完全不同 (diamond/samtools vs
    mmseqs/mafft/RepeatMasker), 合到一起会逼着每一层去填自己根本用不到的字段,
    也说不清哪个阶段的哪个参数才是准的.
    """
    outdir: str = "eve_results"
    elem_dir: str = ""                     # 空 = <outdir>/05_Elements
    threads: int = 40
    jobs: int = 4                          # rm 阶段并行基因组数
    verdicts: tuple = ("viral_supported",)
    min_bs: float = MIN_BS
    min_seq_id: str = CLUSTER_MIN_SEQ_ID
    cov: str = CLUSTER_COV
    er_dist: int = ER_DIST
    mafft_consensus: bool = False
    batch: str = ""                        # name<TAB>genome 批量表 (rm 必需)
    rm_bin: str = REPEATMASKER_DEFAULT     # 也认 RM_BIN 环境变量
    cdd_db: str = ""
    cdd_sens: float = CDD_SENS
    cdd_max_seqs: int = CDD_MAX_SEQS
    rt_ref: str = ""
    run_phylo: bool = False
    min_clusters: int = MIN_CLUSTERS
    mmseqs: str = ""                       # 空 = PATH 里的 mmseqs
    mafft: str = ""
    cmd_timeout: int = 0                   # 单条外部命令超时 (秒); 0=不限
    force: bool = False

    def elem(self):
        """元件层输出目录 (源脚本的 post/ 对应物)."""
        return Path(self.elem_dir) if self.elem_dir else Path(self.outdir) / ELEM_DIR


# ---------------- 通用 IO ----------------
def read_fa(path):
    """FASTA -> {头第一段: 序列}; 与源脚本 read_fa 等价, 但固定 utf-8.

    源脚本用默认编码 (locale): 本机 Windows 默认 GBK, 基因组名/家族名里有非 ASCII
    字符就会 UnicodeDecodeError, 服务器 C locale 下又是另一套行为.
    """
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


def fa_headers(path):
    """只取 FASTA 头 (不读序列) —— 家族名映射只看头, 读整库没必要."""
    heads = []
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith(">"):
                h = line[1:].strip()
                if h:
                    heads.append(h)
    return heads


def write_fa(path, pairs, wrap=80):
    """写 FASTA (显式 newline='\\n').

    Windows 上默认换行翻译会把每个 \\n 写成 \\r\\n, 下游按行 split 的工具会拿到
    带 \\r 的记录 —— 同一份产物在不同平台字节不同, 也没法做逐字节回归.
    """
    with open(path, "w", encoding="utf-8", newline="\n") as fo:
        for name, seq in pairs:
            fo.write(f">{name}\n")
            for i in range(0, len(seq), wrap):
                fo.write(seq[i:i + wrap] + "\n")


def read_table(path):
    """TSV -> (表头 list, 数据行 list[list[str]]); 空文件直接报错.

    返回的行保持原样 (不是 dict): 调用方按**列名**查下标, 历史 13 列 (首列 genome)
    或将来加列都不会错位.
    """
    with open(path, encoding="utf-8") as fh:
        hdr_line = fh.readline().rstrip("\n")
        if not hdr_line.strip():
            raise RuntimeError(f"表是空的 (没有表头): {path}")
        rows = [line.rstrip("\n").split("\t") for line in fh if line.strip()]
    return hdr_line.split("\t"), rows


def write_table(path, hdr, rows):
    """写 TSV: 表头 + 行. 显式 newline='\\n' (理由同 write_fa)."""
    with open(path, "w", encoding="utf-8", newline="\n") as fo:
        fo.write(hdr + "\n")
        for r in rows:
            fo.write("\t".join(str(x) for x in r) + "\n")


def read_members(elem, cid):
    """读 cluster 阶段落的簇成员表 <cid> 的成员 ID 列表 (已排序, 与写盘一致)."""
    p = Path(elem) / f"_members_{cid}.txt"
    if not p.is_file():
        raise RuntimeError(f"缺簇成员表 {p}; 先跑 cluster 阶段 (或 _members_*.txt 被删了)")
    with open(p, encoding="utf-8") as fh:
        return [m for m in fh.read().split() if m]


def resolve_tool(name, configured):
    """解析工具可执行文件: 显式路径按路径处理 (不存在即报错), 裸名字走 PATH.

    与 eve_scan_core.tool_path 的差别只有一点: 那边把任何非空的 configured 都当路径,
    于是 `--rm-bin RepeatMasker` (就是本层默认值) 会报 "配置的工具不存在"。
    默认值表达的是"PATH 里的 RepeatMasker", 必须按 PATH 查, 不能反过来把默认值当成
    用户写错的路径。
    """
    cfg = str(configured or "")
    if cfg and (os.sep in cfg or "/" in cfg):
        return tool_path(name, cfg)          # 含路径分隔符: 交给 tool_path 校验与报错
    found = shutil.which(cfg or name)
    if not found:
        raise SystemExit(f"找不到工具 '{cfg or name}' (PATH 也没有); "
                         f"请安装或用参数显式指定路径")
    return found


# ---------------- 纯逻辑 (可单测, 不碰外部工具) ----------------
def safe_tag(name):
    """家族名/标签清洗: 只保留安全字符, 但**不加** clean_name 的 asm_ 前缀.

    两者用途不同: clean_name 产出的是目录名/文件名主体 (空名要以 asm_ 兜底, 否则
    拼出 <EMPTY>.out 这种路径), safe_tag 产出的是文件名后缀与 FASTA 头里的家族标签,
    空名没有意义也不构成路径逃逸.
    """
    return "".join(ch if ch in SAFE_CHARS else "_" for ch in str(name))


def locus_coords(locus):
    """'chr1:99199-100700' -> (contig, start, end); 解析不出来返回 None.

    取 rsplit(":", 1): 组装的 contig 名里可能本身带 ':' (某些组装把 window 坐标写进
    名字), split(":") 会把 contig 截断成半截名字 —— 与 Stage1 位点命名的口径保持一致.
    """
    parts = str(locus).rsplit(":", 1)
    if len(parts) != 2 or not parts[0]:
        return None
    se = parts[1].split("-")
    if len(se) != 2:
        return None
    try:
        return parts[0], int(se[0]), int(se[1])
    except ValueError:
        return None


def merge_regions(points, er_dist=ER_DIST):
    """元件区归并 (源脚本 p_regions 的核心): points = [(start, end, family), ...].

    同一 contig 上间隔 <= er_dist 的位点并成一个元件区, 家族取并集. 目的与 Caulifinder
    思路一致: EPRV 高度片段化, 一个元件会被切成好几个位点, 归并后再计数才接近真实拷贝数.
    返回 [(start, end, [family, ...]), ...], 家族列表已排序 (输出必须确定).

    与源脚本一致的一处保留: ref_family 为空串时它也算一个"家族" (进入 n_families
    计数与 families 列), 这里不替它过滤 —— 改了就是改口径, 不是修 bug.
    """
    out = []
    if not points:
        return out
    ivs = sorted(points)
    cs, ce, fams = ivs[0][0], ivs[0][1], {ivs[0][2]}
    for s, e, f in ivs[1:]:
        if s <= ce + er_dist:
            ce = max(ce, e)
            fams.add(f)
        else:
            out.append((cs, ce, sorted(fams)))
            cs, ce, fams = s, e, {f}
    out.append((cs, ce, sorted(fams)))
    return out


def candidates_from_rows(name, hdr, rows, verdicts, min_bs=None):
    """位点筛选 (源脚本 p_collect / p_regions 的取行口径), 纯逻辑.

    min_bs=None 时不做证据阈值过滤 (regions 阶段只看 verdict —— 与源脚本一致);
    给数值时 evidence_bs = max(ref_bitscore, rvdb_bitscore), 两条证据通道取强者.
    返回 [(name, locus, family, verdict, ref_sseqid, evidence_bs), ...].
    """
    idx = {h: i for i, h in enumerate(hdr)}
    keep = []

    def get(c, key):
        i = idx.get(key)
        return c[i] if i is not None and i < len(c) else ""

    for c in rows:
        if get(c, "verdict") not in verdicts:
            continue
        bs = 0.0
        for k in ("ref_bitscore", "rvdb_bitscore"):
            try:
                bs = max(bs, float(get(c, k) or 0))
            except ValueError:
                pass
        if min_bs is not None and bs < min_bs:
            continue
        locus = get(c, "locus")
        if not locus:
            continue
        keep.append((name, locus, get(c, "ref_family"), get(c, "verdict"),
                     get(c, "ref_sseqid"), bs))
    return keep


def aggregate_clusters(meta, pairs):
    """聚类聚合 (源脚本 p_cluster 的核心): 家族投票 + 嵌合体初筛 + 最佳成员.

    meta:  member_id -> [genome, locus, family, verdict, ref_sseqid, evidence_bs]
    pairs: cluster_members.tsv 的行 (rep, member)
    返回 (rows, members_by_cid); rows 列序 == ELEMENT_TABLE_HDR.

    两处相对源脚本的**确定性**硬化 (成员集合不变, 只固定顺序):
      1) 成员按名字排序 —— 否则并列长度/并列 bitscore 的胜者跟 mmseqs 输出行序绑定;
      2) 簇按 (-成员数, rep) 排序 —— 源脚本只按成员数降序, 打平时靠输入行序,
         输入一换元件 ID (C000001…) 整体搬家, 产物无法逐字节复现.
    """
    clusters = defaultdict(list)
    for rep, mem in pairs:
        clusters[rep].append(mem)
    for rep in clusters:
        clusters[rep] = sorted(clusters[rep])

    rows, members = [], {}
    for i, rep in enumerate(sorted(clusters, key=lambda r: (-len(clusters[r]), r))):
        ms = clusters[rep]
        cid = f"C{i:06d}"
        fams = Counter(meta.get(m, ["", "", "Other_viral"])[2] for m in ms)
        fam = fams.most_common(1)[0][0] if fams else "Other_viral"
        genomes = {meta[m][0] for m in ms if m in meta}
        minor = [f for f, k in fams.items()
                 if f != fam and k >= max(2, len(ms) // 5)]
        best_m, best_bs = "", 0.0
        for m in ms:
            try:
                bs = float(meta[m][5])
            except (KeyError, ValueError, IndexError):
                continue
            if bs > best_bs:
                best_m, best_bs = m, bs
        rows.append((cid, rep, len(ms), len(genomes), fam,
                     "possible_chimera" if minor else "", best_m, best_bs))
        members[cid] = ms
    return rows, members


def pick_representative(members, seqs):
    """簇内最长成员 (源脚本默认的 consensus 代表口径).

    长度并列时取名字序最小者 —— max() 返回首个最大元素, 成员已排序, 与源脚本
    "按成员顺序取首个最长者"在排序后等价, 且不再随 mmseqs 行序漂移。
    空序列 = 没有序列 (截断/空记录): 源脚本同样不会把它写进文库 (`if seq`), 这里
    提前当"没序列"处理, 好让调用方走同一条告警路径而不是静默丢。
    """
    avail = sorted(m for m in members if seqs.get(m))
    if not avail:
        return None
    return max(avail, key=lambda m: len(seqs[m]))


def consensus_input(members, seqs, mafft_mode):
    """consensus 输入装配, 纯逻辑: 返回 (mode, payload).

    mode="mafft" 时 payload 是参与比对的 {member: seq}; mode="rep" 时是
    [(member, seq)] 单条最长成员.

    源脚本按**成员数** (2 <= n <= 50) 决定是否走比对, 却按"成员序列还在不在
    all_candidates.fa 里"取序列: 上游重跑过 collect 之后成员 ID 会对不上, 源脚本会抓
    出空集合去跑 mafft 直接失败 (mafft 读空文件报错), 报错信息跟真实原因八竿子打不着.
    这里改成按**可用序列数**判断 (空序列不算可用), 不足 2 条就退回最长成员; 只在源脚本
    本来会崩的场景下行为不同, 正常路径逐条一致.
    """
    avail = {m: seqs[m] for m in members if seqs.get(m)}
    if mafft_mode and 2 <= len(members) <= MAFFT_CONS_MAX_N and len(avail) >= 2:
        return "mafft", avail
    rep = pick_representative(members, seqs)
    return "rep", ([(rep, seqs[rep])] if rep else [])


def msa_consensus(records):
    """比对列投票 -> consensus 串 (源脚本 mafft_consensus 的投票部分), 纯逻辑.

    records: [(名字, 比对序列), ...]. 逐列: 非空位碱基少于 2 条的列丢弃; 该列众数
    (统一大写) 得票 >= 非空位数 * 0.4 且不是空位时取该碱基, 否则写 'n' (不确定位).
    """
    cols = defaultdict(Counter)
    for _name, seq in records:
        for i, ch in enumerate(seq):
            cols[i][ch.upper()] += 1
    cons = []
    for i in sorted(cols):
        cnt = cols[i]
        ng = sum(v for k, v in cnt.items() if k != "-")
        if ng < 2:
            continue
        base, freq = cnt.most_common(1)[0]
        cons.append(base if base != "-" and freq >= 0.4 * ng else "n")
    return "".join(cons)


def parse_stages(spec):
    """'collect,cluster' / '1,2' / 'all' -> 阶段集合 (与 eve_genome_scan 同风格)."""
    out = set()
    for tok in str(spec).split(","):
        tok = tok.strip()
        if tok in ("all", ""):
            return set(ALL_STAGES)
        if tok not in STAGE_ALIAS and tok not in ALL_STAGES:
            raise SystemExit(f"未知阶段: {tok} (可选: 1..7 / "
                             f"{', '.join(ALL_STAGES)} / all)")
        out.add(STAGE_ALIAS.get(tok, tok))
    return out


def parse_batch(path):
    """批量表解析: 每行 NAME<TAB>基因组路径; 空行/注释/列数!=2 跳过.

    与 eve_screen.py -B 同一格式同一张表. 名字过 clean_name: 它会落成 rm/<NAME>/ 目录名,
    不清洗 (名字里带 '/' 或 '..') 就是路径逃逸.
    """
    todo = []
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            c = line.rstrip("\n").split("\t")
            if len(c) == 2 and not c[0].startswith("#"):
                todo.append((clean_name(c[0].strip()),
                             os.path.expanduser(c[1].strip())))
    return todo


def pending_rm_jobs(todo, rm_root):
    """待跑的 RM 基因组: 已有 <NAME>.out 的跳过 (断点粒度到单个基因组)."""
    rm_root = Path(rm_root)
    return [(n, g) for n, g in todo
            if not (rm_root / n / f"{n}.out").exists()]


def build_lib_fam_map(headers):
    """consensus 文库头 -> {cluster_id: 家族}, 纯逻辑.

    头形如 'C000123_Caulimoviridae'. 家族名被 safe_tag 洗过, 里面可能本来就含 '_'
    ('Caulimoviridae(endo)' -> 'Caulimoviridae_endo_'), 所以只按第一个 '_' 切.
    """
    m = {}
    for head in headers:
        cid = head.split("_")[0]
        m[cid] = head[len(cid) + 1:] if "_" in head else head
    return m


def parse_rm_out_lines(lines):
    """RepeatMasker .out 数据行 -> (query_begin, query_end, repeat_name).

    .out 的官方列序 (0 基): 0 SW score, 1-3 perc div/del/ins, 4 query 名,
    5 query begin, 6 query end, 7 (left), 8 链 ('+' / 'C'), 9 matching repeat,
    10 repeat class/family, 11-13 repeat begin/end/(left), 14 ID。
    实测多份真 .out (RepeatMasker 4.x) 的 (left) 列都写作 `(1234)`、不带内空格,
    所以按空白切正好是这 15 列; 两行表头与 "There were no repetitive..." 说明行靠
    第 5/6 列不是数字挡掉, 另有"第 8 列必须是链"的守卫 —— 万一哪个版本的列序不同,
    宁可跳过也不要把命计算到一个叫 '+' 的家族头上。

    **源脚本这里取错了一列** (eve_kingdom_post.py:387-391): 它用 c[4]/c[5] 当
    "query begin/end", 还拿它们做 isdigit() 守卫。c[4] 实际是 query 序列名, 对真实
    .out 永远不是数字 —— 于是每一行真实命中都被当成表头跳过, family_quant.tsv
    只剩表头 (退出码 0, 日志还写着 "RM QUANT DONE (0 families x 0 genomes)");
    反过来若 contig 名恰好是纯数字 ('1','2'…), 守卫放行, 又按
    int(c[5]) - int(c[4]) + 1 = 起点 - contig号 + 1 算出垃圾 covered_bp。
    这里按官方列序取 5/6。
    """
    for line in lines:
        c = line.split()
        if (len(c) < 11 or not c[5].isdigit() or not c[6].isdigit()
                or c[8] not in ("+", "C")):
            continue
        yield int(c[5]), int(c[6]), c[9].split("#")[0]


def accumulate_quant(hits, lib_fam, genome):
    """把 (begin, end, name) 命中累加进 {(family, genome): [copies, covered_bp]}.

    一个家族在一个基因组里的"拷贝数"按 RepeatMasker 命中条数算 (碎片化校正交给
    regions 阶段的元件区归并), covered_bp 是该家族全部命中的碱基和.
    """
    quant = defaultdict(lambda: [0, 0])
    for b, e, name in hits:
        cid = name.split("_")[0]
        fam = lib_fam.get(cid, name)
        quant[(fam, genome)][0] += 1
        quant[(fam, genome)][1] += e - b + 1
    return quant


# ---------------- 阶段实现 ----------------
def discover_summaries(outdir):
    """遍历一次筛查输出的基因组: 返回 [(NAME, summary 路径), ...] (按名字排序).

    以 04_Summary/*_eve_summary.tsv 为遍历键, 而不是 01_Loci/<NAME>/ 目录: summary 才是
    本层的输入契约 (verdict 与三条证据都在里面), 名字直接来自文件名, 不会把只跑了
    半截 (没有 summary) 的基因组带进来, 也不会受 01_Loci 目录里遗留垃圾目录影响。
    排序保证表行序与文件系统遍历顺序无关。
    """
    d = Path(outdir) / STAGE_DIRS["summary"]
    if not d.is_dir():
        raise RuntimeError(f"找不到 {STAGE_DIRS['summary']}/: {outdir} 不是一次 EVE "
                           f"筛查的输出根目录 (先跑 eve_screen.py)")
    out = [(p.name[: -len(SUMMARY_SUFFIX)], p)
           for p in sorted(d.glob(f"*{SUMMARY_SUFFIX}"))
           if len(p.name) > len(SUMMARY_SUFFIX)]
    return out


def stage_done(cfg, stage):
    """断点判定: 该阶段声明的产物全都在."""
    elem = cfg.elem()
    return all((elem / f).exists() for f in STAGE_OUTPUTS[stage])


def stage_collect(cfg, log):
    """collect: 04_Summary + 01_Loci -> all_candidates.fa + candidates_meta.tsv.

    ID 规则 <NAME>__<locus>: locus 串 (chr1:99199-100700) 在不同基因组里会重名,
    不加前缀跨基因组聚类就会把不同基因组的同位点当同一条序列。
    """
    elem = cfg.elem()
    verdicts = set(cfg.verdicts)
    pair_list, out_rows, n_genome = [], [], 0
    for name, summ in discover_summaries(cfg.outdir):
        loci_fa = (Path(cfg.outdir) / STAGE_DIRS["loci"] / name
                   / f"{name}{LOCI_SUFFIX}")
        if not loci_fa.is_file():
            logmsg(f"WARN collect: {name} 缺 {LOCI_SUFFIX} 序列文件, 跳过 ({loci_fa})")
            continue
        # 两侧都过 locus_key 再 join (理由见模块 docstring)
        seqs = {locus_key(k): v for k, v in read_fa(loci_fa).items()}
        hdr, rows = read_table(summ)
        n = 0
        for gname, locus, fam, verdict, sseqid, bs in candidates_from_rows(
                name, hdr, rows, verdicts, cfg.min_bs):
            seq = seqs.get(locus_key(locus))
            if not seq:
                continue
            pair_list.append((f"{gname}__{locus}", seq))
            out_rows.append((gname, locus, fam, verdict, sseqid, bs))
            n += 1
        n_genome += 1
        logmsg(f"collect {name}: {n} candidates")
    if not pair_list:
        raise RuntimeError(
            f"没有任何 viral_supported 位点被收集 (verdicts={sorted(verdicts)}, "
            f"min_bs={cfg.min_bs}) —— 检查 04_Summary 是否有该 verdict 的行、"
            f"以及 loci.fa 头能否与 locus 列对上")
    write_fa(elem / "all_candidates.fa", pair_list)
    write_table(elem / "candidates_meta.tsv", CAND_META_HDR, out_rows)
    return len(pair_list), n_genome


def stage_cluster(cfg, log):
    """cluster: all_candidates.fa --easy-linclust--> element_table.tsv (+ 簇成员表)."""
    elem = cfg.elem()
    clu = elem / "cluster_members.tsv"
    if cfg.force or not clu.exists():
        prefix, tmp = elem / "_clu", elem / "_clu_tmp"
        run([resolve_tool("mmseqs", cfg.mmseqs), "easy-linclust",
             str(elem / "all_candidates.fa"), str(prefix), str(tmp),
             "--min-seq-id", cfg.min_seq_id, "-c", cfg.cov,
             "--cov-mode", CLUSTER_COV_MODE, "--threads", str(cfg.threads)],
            log, err_path=elem / "mmseqs_linclust.err", timeout=cfg.cmd_timeout)
        src = elem / "_clu_cluster.tsv"
        if not src.is_file():
            # 源脚本这里直接 src.rename(), 文件不存在就抛 FileNotFoundError 指向中间产物;
            # 退出码 0 却没产出是 mmseqs 侧的问题 (空输入/被包装脚本吞掉输出), 要指名说清
            raise RuntimeError(f"mmseqs easy-linclust 未产出 {src}; "
                               f"看 {elem / 'mmseqs_linclust.err'}")
        src.replace(clu)          # replace 而不是 rename: Windows 上目标存在时 rename 直接失败
    meta = {}
    with open(elem / "candidates_meta.tsv", encoding="utf-8") as fh:
        fh.readline()
        for line in fh:
            c = line.rstrip("\n").split("\t")
            if len(c) >= 6:
                meta[f"{c[0]}__{c[1]}"] = c
    pairs = []
    with open(clu, encoding="utf-8") as fh:
        for line in fh:
            c = line.rstrip("\n").split("\t")
            if len(c) == 2:
                pairs.append((c[0], c[1]))
    rows, members = aggregate_clusters(meta, pairs)
    write_table(elem / "element_table.tsv", ELEMENT_TABLE_HDR, rows)
    for cid, ms in sorted(members.items()):
        with open(elem / f"_members_{cid}.txt", "w", encoding="utf-8",
                  newline="\n") as mf:
            mf.write("\n".join(ms) + "\n")
    return len(rows)


def stage_regions(cfg, log):
    """regions: 每基因组的 ±er_dist 同家族位点归并 -> element_regions.tsv.

    与 collect 的不同: 这里不设证据阈值 (只看 verdict), 因为 regions 是计数校正用的
    元数据表, 不是候选集 —— 阈值一变, 元件区就跟着变, 拷贝数校正也就不可比。
    """
    elem = cfg.elem()
    verdicts = set(cfg.verdicts)
    rows, n_genome = [], 0
    for name, summ in discover_summaries(cfg.outdir):
        hdr, data = read_table(summ)
        pts = defaultdict(list)
        for (_g, locus, fam, _v, _s, _bs) in candidates_from_rows(
                name, hdr, data, verdicts, None):
            co = locus_coords(locus)
            if not co:
                continue
            pts[co[0]].append((co[1], co[2], fam))
        k = 0
        for contig in sorted(pts):
            for cs, ce, fams in merge_regions(pts[contig], cfg.er_dist):
                k += 1
                rows.append((name, f"{name}__ER{k:05d}", contig, cs, ce,
                             len(fams), ";".join(fams)))
        n_genome += 1
    write_table(elem / "element_regions.tsv", ELEMENT_REGIONS_HDR, rows)
    return len(rows), n_genome


def mafft_consensus(cfg, log, seqs):
    """比对共识 (源脚本 --mafft-consensus 路径): mafft --auto + 列投票.

    --thread 1 是源脚本的选择 (比对条数 <= 50, 多线程的收益抵不上 mafft 内部
    线程间的非确定性), 保留不动.
    """
    elem = cfg.elem()
    items = sorted(seqs.items(), key=lambda kv: (-len(kv[1]), kv[0]))
    tmp, aln = elem / "_cons_tmp.fa", elem / "_cons_tmp.aln"
    write_fa(tmp, items)
    try:
        run([resolve_tool("mafft", cfg.mafft), "--auto", "--quiet", "--anysymbol",
             str(tmp), "-o", str(aln), "--thread", "1"], log,
            err_path=elem / "mafft.err", timeout=cfg.cmd_timeout)
        recs = []
        with open(aln, encoding="utf-8") as fh:
            hdr, buf = None, []
            for line in fh:
                if line.startswith(">"):
                    if hdr:
                        recs.append((hdr, "".join(buf)))
                    hdr, buf = line[1:].strip(), []
                elif hdr:
                    buf.append(line.strip())
            if hdr:
                recs.append((hdr, "".join(buf)))
    finally:
        tmp.unlink(missing_ok=True)
        aln.unlink(missing_ok=True)
    return msa_consensus(recs)


def stage_consensus(cfg, log):
    """consensus: element_table + 簇成员 --(最长成员 | mafft 共识)--> 文库.

    文库名 <cid>_<家族>: rm/cdd 阶段要从名字反查家族, 家族名先过 safe_tag。
    """
    elem = cfg.elem()
    seqs_all = read_fa(elem / "all_candidates.fa")
    names = {}
    with open(elem / "element_table.tsv", encoding="utf-8") as fh:
        fh.readline()
        for line in fh:
            c = line.rstrip("\n").split("\t")
            if len(c) >= 5:
                names[c[0]] = c[4]
    rows, n_mafft, n_lost = [], 0, 0
    for cid in sorted(names):
        members = read_members(elem, cid)
        mode, payload = consensus_input(members, seqs_all, cfg.mafft_consensus)
        if mode == "mafft":
            seq = mafft_consensus(cfg, log, payload)
            n_mafft += 1
        else:
            seq = payload[0][1] if payload else ""
        if not seq:
            # 源脚本在这里静默丢簇 —— 文库比 element_table 少几行没人知道, 下游按
            # 元件 ID 回查就找不到序列. 至少让它可追。
            n_lost += 1
            logmsg(f"WARN consensus: 簇 {cid} 没有任何可用成员序列, 跳过")
            continue
        rows.append((f"{cid}_{safe_tag(names[cid])}", seq))
    if n_lost:
        logmsg(f"WARN consensus: {n_lost}/{len(names)} 个簇没出序列")
    write_fa(elem / "eve_consensus_library.fa", rows)
    return len(rows), ("mafft" if n_mafft else "longest-member")


def rm_worker(job):
    """单基因组 RepeatMasker (进程池 worker).

    进池前把要用的东西摊平进 job, 不传 cfg / logger: Windows spawn 下子进程会重新
    import 本模块, 只有可 pickle 的标量才稳; logger 由 worker 自己建 (与
    eve_genome_scan.run_genome 同一套做法)。
    """
    name, genome, elem, threads, log_path, rm_bin, timeout = job
    log = setup_logger(log_path)
    elem = Path(elem)
    rdir = elem / "rm" / name
    rdir.mkdir(parents=True, exist_ok=True)
    marker = rdir / f"{name}.out"
    if not marker.exists():
        run([rm_bin, "-lib", str(elem / "eve_consensus_library.fa"),
             "-pa", str(threads), "-dir", str(rdir), "-no_is", "-norna",
             str(genome)], log, err_path=rdir / "repeatmasker.err",
            timeout=timeout)
        base = Path(genome).name
        for suf in (".out", ".tbl", ".masked", ".cat.gz"):
            src = rdir / (base + suf)
            if src.exists():
                src.replace(rdir / f"{name}{suf}")
    return name, str(marker)


def stage_rm(cfg, log):
    """rm: consensus 文库回比各基因组 -> family_quant.tsv (每家族拷贝数/覆盖度).

    单基因组失败不连带整批: 逐个 future 兜住, 把失败的基因组名攒起来一起报, 成功
    的那些 .out 留在盘上 (下次重发同命令只补失败的)。任何一个失败都不写 family_quant,
    否则会拿半份数据当完整结果。
    """
    elem = cfg.elem()
    lib = elem / "eve_consensus_library.fa"
    todo = parse_batch(cfg.batch)
    if not todo:
        raise RuntimeError(f"批量表里没有可用行 (每行 NAME<TAB>基因组路径): {cfg.batch}")
    if not lib.is_file() or lib.stat().st_size == 0:
        # 空文库时 RepeatMasker 会以无命中/报错收场; 按"零结果"落表更符合断点语义
        logmsg(f"WARN rm: consensus 文库为空 ({lib}), 写空 family_quant.tsv")
        write_table(elem / "family_quant.tsv", FAMILY_QUANT_HDR, [])
        return 0, 0
    rm_root = elem / "rm"
    jobs = todo if cfg.force else pending_rm_jobs(todo, rm_root)
    log_path = str(Path(cfg.outdir) / "logs" / "eve_element_rm.log")
    # 在父进程里把可执行文件解析成绝对路径再进池: Windows 上 subprocess 只按
    # CreateProcess 的规则找 .exe, 裸名字 ("RepeatMasker") 到不了 RepeatMasker.cmd /
    # 无扩展名的包装脚本; 而且 spawn 出来的子进程不该再依赖 PATH 长什么样.
    rm_exe = resolve_tool("RepeatMasker", cfg.rm_bin)
    if jobs:
        logmsg(f"RM: {len(jobs)} 个基因组待屏蔽 (jobs={cfg.jobs}, threads={cfg.threads})")
        joblist = [(n, g, str(elem), cfg.threads, log_path, rm_exe,
                    cfg.cmd_timeout) for n, g in jobs]
        failed = []
        if cfg.jobs > 1:
            with ProcessPoolExecutor(max_workers=cfg.jobs) as ex:
                futs = {ex.submit(rm_worker, j): j[0] for j in joblist}
                for fu in as_completed(futs):
                    try:
                        n, _ = fu.result()
                        logmsg(f"RM done {n}")
                    except Exception as exc:                    # noqa: BLE001
                        failed.append(f"{futs[fu]}({type(exc).__name__})")
        else:
            for j in joblist:
                try:
                    n, _ = rm_worker(j)
                    logmsg(f"RM done {n}")
                except Exception as exc:                        # noqa: BLE001
                    failed.append(f"{j[0]}({type(exc).__name__})")
        if failed:
            raise RuntimeError(f"{len(failed)}/{len(joblist)} 个基因组 RepeatMasker "
                               f"失败: {', '.join(sorted(failed)[:5])}"
                               + (" ..." if len(failed) > 5 else ""))
    lib_fam = build_lib_fam_map(fa_headers(lib))
    quant = defaultdict(lambda: [0, 0])
    for name, _genome in todo:
        outf = rm_root / name / f"{name}.out"
        if not outf.is_file():
            continue
        with open(outf, encoding="utf-8", errors="replace") as fh:
            part = accumulate_quant(parse_rm_out_lines(fh), lib_fam, name)
        for key, v in part.items():
            quant[key][0] += v[0]
            quant[key][1] += v[1]
    fams = sorted({f for f, _ in quant})
    genomes = sorted({g for _, g in quant})
    rows = [(f, g, quant[(f, g)][0], quant[(f, g)][1])
            for f in fams for g in genomes if quant[(f, g)][0]]
    write_table(elem / "family_quant.tsv", FAMILY_QUANT_HDR, rows)
    return len(fams), len(genomes)


def stage_cdd(cfg, log):
    """cdd: consensus 文库 --mmseqs 翻译搜索(六框)--> consensus_domains.tsv.

    --search-type 2 (翻译搜索, 核酸查询 vs 蛋白库) 与 rpstblastn 语义等价 —— 原脚本
    用 rpstblastn, 但本机 CDD 是 msa2profile 建的 mmseqs profile 库, rpstblastn 读不了.
    中间库目录 (_cdd_mm) 每次重跑前清掉: mmseqs createdb 对已存在的库目录报错,
    留着上一次的库会让"改过 consensus 文库之后重跑"拿到旧库的结果.
    """
    elem = cfg.elem()
    lib = elem / "eve_consensus_library.fa"
    out = elem / "consensus_domains.tsv"
    if not lib.is_file() or lib.stat().st_size == 0:
        # 空文库按"零命中"落一个空文件 (与 eve_scan_core.stage3_annotate 对零候选的处理
        # 一致): 这个表本来是 mmseqs convertalis 直出、没有表头的, 补一行表头反而会让
        # 下游按位置解析的第一行变成垃圾
        logmsg(f"WARN cdd: consensus 文库为空 ({lib}), 写空 consensus_domains.tsv")
        out.write_bytes(b"")
        return 0
    mm = elem / "_cdd_mm"
    shutil.rmtree(mm, ignore_errors=True)
    mm.mkdir(parents=True, exist_ok=True)
    qdb, rdb, tmp = mm / "qDB", mm / "rDB", mm / "tmp"
    ms, target = resolve_tool("mmseqs", cfg.mmseqs), os.path.expanduser(cfg.cdd_db)
    err = elem / "mmseqs_cdd.err"
    run([ms, "createdb", str(lib), str(qdb), "-v", "1"], log, err_path=err,
        timeout=cfg.cmd_timeout)
    run([ms, "search", str(qdb), target, str(rdb), str(tmp), "--search-type", "2",
         "-s", str(cfg.cdd_sens), "-v", "1", "--threads", str(cfg.threads),
         "--max-seqs", str(cfg.cdd_max_seqs)], log, err_path=err,
        timeout=cfg.cmd_timeout)
    run([ms, "convertalis", str(qdb), target, str(rdb), str(out),
         "--format-output", CDD_CONVERT_FIELDS, "-v", "1"], log, err_path=err,
        timeout=cfg.cmd_timeout)
    with open(out, encoding="utf-8") as fh:
        n = sum(1 for line in fh if line.strip())
    return n


def stage_phylo(cfg, log):
    """phylo: 每家族代表序列 + RT 参考集 -> phylo/ 与 phylo_commands.sh.

    命令脚本是给人/集群跑的成品 (mafft -> trimal -> iqtree); --run-phylo 会在本机
    直接执行它.
    """
    elem = cfg.elem()
    seqs_all = read_fa(elem / "all_candidates.fa")
    fam_cids = defaultdict(list)
    with open(elem / "element_table.tsv", encoding="utf-8") as fh:
        fh.readline()
        for line in fh:
            c = line.rstrip("\n").split("\t")
            if len(c) >= 5:
                fam_cids[c[4]].append(c[0])
    pdir = elem / "phylo"
    pdir.mkdir(parents=True, exist_ok=True)
    cmd_lines = ["#!/bin/bash", "set -euo pipefail"]
    n_fam = 0
    # 排序键带上家族名: 源脚本只按簇数降序, 并列时脚本里命令的顺序随 dict 插入序漂
    for fam, cids in sorted(fam_cids.items(), key=lambda kv: (-len(kv[1]), kv[0])):
        if len(cids) < cfg.min_clusters:
            continue
        tag = safe_tag(fam)
        reps = []
        for cid in sorted(cids):
            best = pick_representative(read_members(elem, cid), seqs_all)
            if best:
                reps.append((f"{cid}__{best}", seqs_all[best]))
        write_fa(pdir / f"{tag}_reps.fa", reps)
        n_fam += 1
        if cfg.rt_ref:
            ref = os.path.expanduser(cfg.rt_ref)
            cmd_lines += [
                f"# {fam}: {len(cids)} clusters",
                f"cat '{pdir}/{tag}_reps.fa' '{ref}' > '{pdir}/{tag}_query.fa'",
                f"mafft --auto --anysymbol --thread {cfg.threads} "
                f"'{pdir}/{tag}_query.fa' > '{pdir}/{tag}_aln.fa'",
                f"trimal -automated1 -in '{pdir}/{tag}_aln.fa' "
                f"-out '{pdir}/{tag}_trim.fa'",
                f"iqtree -s '{pdir}/{tag}_trim.fa' -m MFP -B 1000 "
                f"-T {cfg.threads} --prefix '{pdir}/{tag}_tree'",
                ""]
        else:
            cmd_lines += [f"# {fam}: {len(cids)} clusters — 准备 RT 参考集 fasta "
                          f"(外源属 + copia/gypsy 外群) 后加 --rt-ref 重跑", ""]
    script = elem / PHYLO_SCRIPT
    with open(script, "w", encoding="utf-8", newline="\n") as fo:
        fo.write("\n".join(cmd_lines) + "\n")
    if cfg.run_phylo:
        # 源脚本是逐行 line.split() 当 argv 执行 (eve_kingdom_post.py:473-477), 那条路
        # 走不通: 'set -euo pipefail' 会被当成可执行文件 (bash 内建, 没有 /usr/bin/set),
        # mafft 那行里的 '> out.fa' 重定向与 cat 的引号也都会被当成普通参数 —— 一跑就崩.
        # 脚本本来就是写给 bash 的, 直接交给 bash 执行 (set -e 保证中途失败即停)。
        bash = shutil.which("bash")
        if not bash:
            raise RuntimeError("--run-phylo 需要 bash 执行生成的脚本")
        run([bash, str(script)], log, err_path=elem / "phylo.err",
            timeout=cfg.cmd_timeout)
    return n_fam


# ---------------- 工具体检 ----------------
def _exe_ok(exe):
    """可执行文件是否可用: 显式路径优先, 否则查 PATH (与 eve_scan_core 同判据).

    调用方传 `cfg.xxx or "<默认名>"` —— 空字符串的意思是"用 PATH 里的同名工具",
    不是"没有工具"; 直接拿空串来判会把默认路径当成缺失, 把正常的 PATH 安装挡在门外.
    """
    if not exe:
        return False
    p = Path(os.path.expanduser(str(exe)))
    return p.is_file() or bool(shutil.which(str(p)))


def check_tools(cfg, stages):
    """开跑前体检: 返回缺失的工具/参数清单 (空 = 全都就绪).

    只查**请求的阶段**真正要用的东西 —— 与 eve_screen.check_tools 一次查全部数据库
    不同: 本层的 cdd/rm/phylo 是可选加装层, 全查会让只想跑四个默认阶段的人被
    --cdd-db 挡住。
    源脚本对应的三处是**静默跳过** (eve_kingdom_post.py:526/529/473: rm 缺 -g、
    cdd 缺 --cdd-db、--run-phylo 缺 --rt-ref 时那一阶段什么都不做, 退出码 0),
    用户以为跑完了去取表, 只有表头。这里改成开跑前报名字。
    """
    missing = []
    if {"cluster", "cdd"} & set(stages):
        if not _exe_ok(cfg.mmseqs or "mmseqs"):
            missing.append(f"mmseqs={cfg.mmseqs or '(PATH 里没有 mmseqs)'}")
    if "consensus" in stages and cfg.mafft_consensus:
        if not _exe_ok(cfg.mafft or "mafft"):
            missing.append(f"mafft={cfg.mafft or '(PATH 里没有 mafft)'} "
                           f"(--mafft-consensus 需要)")
    if "rm" in stages:
        if not _exe_ok(cfg.rm_bin):
            missing.append(f"RepeatMasker={cfg.rm_bin or '(未配置)'} "
                           f"(可用 --rm-bin 或 RM_BIN 指定)")
        if not cfg.batch:
            missing.append("-B/--batch=(未配置; rm 阶段需要 name<TAB>基因组路径 批量表)")
    if "cdd" in stages:
        if not cfg.cdd_db:
            missing.append("--cdd-db=(未配置; cdd 阶段需要 mmseqs profile 格式的 CDD 库前缀)")
        elif not Path(os.path.expanduser(cfg.cdd_db)).exists():
            missing.append(f"--cdd-db={cfg.cdd_db} (路径不存在)")
    if "phylo" in stages and cfg.run_phylo:
        if not cfg.rt_ref:
            missing.append("--rt-ref=(未配置; --run-phylo 需要 RT 参考集 fasta)")
        elif not Path(os.path.expanduser(cfg.rt_ref)).exists():
            missing.append(f"--rt-ref={cfg.rt_ref} (文件不存在)")
        for exe, label in ((cfg.mafft or "mafft", "mafft"), ("trimal", "trimal"),
                           ("iqtree", "iqtree")):
            if not _exe_ok(exe):
                missing.append(f"{label}={exe or '(未配置)'} (--run-phylo 需要)")
    return missing


# ---------------- CLI ----------------
def build_parser():
    p = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="元件级 (element-level) EVE 后处理层: 跨基因组聚类 → 元件区 → "
                    "consensus → 拷贝数/结构域/系统发育")
    g = p.add_argument_group("输入/输出")
    g.add_argument("-o", "--outdir", default="eve_results",
                   help="eve_screen.py 的输出根目录 (也是本层输出根目录)")
    g.add_argument("-B", "--batch", default=None,
                   help="批量 TSV: 每行 NAME<TAB>基因组路径 (rm 阶段必需; 与 "
                        "eve_screen.py -B 同一张表)")
    g.add_argument("--stages", default=",".join(DEFAULT_STAGES),
                   help="阶段: collect/cluster/regions/consensus/rm/cdd/phylo 逗号组合 "
                        "(别名 1..7) 或 all; all 会连带 rm/cdd/phylo, 需要相应参数")
    g.add_argument("--force", action="store_true",
                   help="忽略断点, 重跑请求的阶段 (rm 阶段连带忽略每基因组的 .out 标记)")
    g = p.add_argument_group("运行")
    g.add_argument("-t", "--threads", type=int, default=40)
    g.add_argument("-J", "--jobs", type=int, default=4,
                   help="rm 阶段并行基因组数 (总核数 = jobs x threads)")
    g.add_argument("--cmd-timeout", type=int, default=0,
                   help="单条外部命令超时 (秒); 0=不限 (防工具挂死占核)")
    g = p.add_argument_group("筛选/聚类口径 (与源脚本一致)")
    g.add_argument("--verdicts", default="viral_supported",
                   help="collect/regions 收哪些 verdict (逗号分隔)")
    g.add_argument("--min-bs", type=float, default=MIN_BS,
                   help="evidence_bs 下限 (max(ref_bitscore, rvdb_bitscore))")
    g.add_argument("--id", default=CLUSTER_MIN_SEQ_ID,
                   help="easy-linclust --min-seq-id")
    g.add_argument("--cov", default=CLUSTER_COV, help="easy-linclust -c")
    g.add_argument("--er-dist", type=int, default=ER_DIST,
                   help="元件区归并距离 (nt, ±同家族位点)")
    g.add_argument("--mafft-consensus", action="store_true",
                   help="consensus 走 mafft 比对共识 (默认取簇内最长成员)")
    g = p.add_argument_group("需要额外依赖的阶段")
    g.add_argument("--rm", action="store_true",
                   help="把 rm 阶段加进本次运行 (需 -B 与 RepeatMasker)")
    g.add_argument("--rm-bin", default=REPEATMASKER_DEFAULT,
                   help="RepeatMasker 可执行文件 (默认取 RM_BIN 环境变量, 再默认 PATH)")
    g.add_argument("--cdd-db", default="",
                   help="mmseqs profile 格式的 CDD 库前缀 (给定即把 cdd 阶段加进本次运行)")
    g.add_argument("--cdd-sens", type=float, default=CDD_SENS,
                   help="mmseqs search -s")
    g.add_argument("--cdd-max-seqs", type=int, default=CDD_MAX_SEQS,
                   help="mmseqs search --max-seqs")
    g.add_argument("--rt-ref", default="", help="RT 参考集 fasta (--run-phylo 必需)")
    g.add_argument("--run-phylo", action="store_true",
                   help="本机执行 phylo_commands.sh (需 --rt-ref 与 mafft/trimal/iqtree)")
    g.add_argument("--min-clusters", type=int, default=MIN_CLUSTERS,
                   help="phylo: 家族至少要有这么多个簇才出代表集")
    g = p.add_argument_group("工具")
    g.add_argument("--mmseqs", default=None, help="mmseqs 可执行文件 (默认 PATH)")
    g.add_argument("--mafft", default=None, help="mafft 可执行文件 (默认 PATH)")
    return p


def cfg_from_args(args):
    return ElementConfig(
        outdir=str(safe_dir(args.outdir)), elem_dir="",
        threads=args.threads, jobs=args.jobs,
        verdicts=tuple(v for v in args.verdicts.split(",") if v),
        min_bs=args.min_bs, min_seq_id=args.id, cov=args.cov,
        er_dist=args.er_dist, mafft_consensus=args.mafft_consensus,
        batch=os.path.expanduser(args.batch) if args.batch else "",
        rm_bin=args.rm_bin or "RepeatMasker",
        cdd_db=args.cdd_db, cdd_sens=args.cdd_sens,
        cdd_max_seqs=args.cdd_max_seqs, rt_ref=args.rt_ref,
        run_phylo=args.run_phylo, min_clusters=args.min_clusters,
        mmseqs=args.mmseqs or "", mafft=args.mafft or "",
        cmd_timeout=args.cmd_timeout, force=args.force)


# [done] 行的措辞 (detail 是各阶段返回的小元组; 措辞集中在这里, 免得 [done] 行
# 里的数字单位跟阶段实现各说各话)
STAGE_DONE_MSG = {
    "collect":   lambda d: f"{d[0]} 个位点 (来自 {d[1]} 个基因组) -> "
                           f"all_candidates.fa + candidates_meta.tsv",
    "cluster":   lambda d: f"{d} 个元件 (easy-linclust id={CLUSTER_MIN_SEQ_ID} "
                           f"cov={CLUSTER_COV} cov-mode={CLUSTER_COV_MODE}) -> "
                           f"element_table.tsv",
    "regions":   lambda d: f"{d[0]} 个元件区 (±{ER_DIST}bp 同家族归并, "
                           f"{d[1]} 个基因组) -> element_regions.tsv",
    "consensus": lambda d: f"{d[0]} 条 consensus ({d[1]}) -> "
                           f"eve_consensus_library.fa",
    "rm":        lambda d: f"{d[0]} 个家族 x {d[1]} 个基因组 -> family_quant.tsv",
    "cdd":       lambda d: f"{d} 条结构域命中 -> consensus_domains.tsv",
    "phylo":     lambda d: f"{d} 个家族代表序列 -> phylo/ + {PHYLO_SCRIPT}",
}


def main(argv=None):
    args = build_parser().parse_args(argv)
    cfg = cfg_from_args(args)
    log = setup_logger(str(Path(cfg.outdir) / "logs" / "eve_element_db.log"))
    requested = parse_stages(args.stages)
    if args.rm:
        requested.add("rm")
    if args.cdd_db:
        requested.add("cdd")
    if args.run_phylo:
        requested.add("phylo")

    elem = cfg.elem()
    elem.mkdir(parents=True, exist_ok=True)
    missing = check_tools(cfg, requested)
    if missing:
        logmsg("缺少工具/参数: " + "; ".join(missing))
        return 1

    t0 = time.time()
    logmsg(f"元件层: stages={','.join(s for s in ALL_STAGES if s in requested)}, "
           f"threads={cfg.threads}, jobs={cfg.jobs}, out={elem}")
    # 固定按 ALL_STAGES 的顺序跑: 阶段之间有真实依赖 (cluster 要 collect 的产物,
    # consensus 要 cluster 的, rm/cdd 要 consensus 的), 按用户给的顺序跑会得到
    # 取决于参数写法的结果
    for stage in ALL_STAGES:
        if stage not in requested:
            continue
        if not cfg.force and stage_done(cfg, stage):
            logmsg(f"[done] {stage}: 产物已存在, 跳过 "
                   f"({', '.join(STAGE_OUTPUTS[stage])}); --force 重跑")
            continue
        fn = {"collect": stage_collect, "cluster": stage_cluster,
              "regions": stage_regions, "consensus": stage_consensus,
              "rm": stage_rm, "cdd": stage_cdd, "phylo": stage_phylo}[stage]
        try:
            detail = fn(cfg, log)
        except SystemExit:
            raise
        except Exception as exc:                                # noqa: BLE001
            logmsg(f"[fail] {stage}: {type(exc).__name__}: {exc}")
            return 1
        logmsg(f"[done] {stage}: {STAGE_DONE_MSG[stage](detail)}")
    logmsg(f"[done] all: 请求的阶段全部完成 "
           f"({','.join(s for s in ALL_STAGES if s in requested)}, "
           f"{time.time() - t0:.0f}s) -> {elem}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
