#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
eve_scan_core.py — 内源性病毒元件 (EVE) 筛查核心逻辑
====================================================
对宿主参考基因组做三阶段 EVE 筛查 (内源性病毒化石位点):

  Stage 1  discover  正向发现: 50kb 滑窗 → diamond blastx vs 病毒参考蛋白库
                      → 坐标还原 → 位点合并 (300nt) → 抽序列 → 最佳参考命中 + 科级映射
  Stage 2  verdict   双侧判定: loci vs 植物/病毒拆分库 (id2div 映射)
                      → viral_supported / host_like / undetermined
                      (宿主蛋白误匹配排除: 两侧 bitscore 对比, 防植物基因假阳性)
  Stage 3  annotate  RVDB 深度归属: viral_supported + undetermined 位点 vs RVDB

输出目录结构 (以 --output_dir 为根):
  01_Loci/<NAME>/         Stage1 产物 (fna/chunks/s1_raw/s1_hits/loci.bed/loci.fa/s1_best)
  02_Verdict/<NAME>/      Stage2 产物 (s2_raw/s2_verdict)
  03_RVDB/<NAME>/         Stage3 产物 (cand3.bed/cand3.fa/s3_rvdb)
  04_Summary/             <NAME>_eve_summary.tsv (+ 可选 <NAME>_viroid.tsv)
  logs/                   每基因组日志 (编排器写入)
  --merge 汇总产物写在根目录: kingdom_summary.tsv / kingdom_loci_all.tsv.gz /
                            family_by_genome.tsv

断点续传: 每个阶段产物文件存在即视为完成 (空文件 = 已完成且零结果),
--force 删除旧产物重跑.

来源: 移植自 hi-fever/eve_kingdom.py (2026-09-19 定版, 246 服务器千种植物验证),
数据库路径/工具路径全部参数化 (EveConfig), 适配 MMPV-RNA 配置体系.
数据文件读写统一 utf-8 (服务器 C locale / 本机 GBK 默认编码下不炸);
外部命令支持超时 (cmd_timeout>0) 防工具挂死白占并行核.
"""

import gzip
import os
import re
import shutil
import subprocess
import sys
import time
from contextlib import nullcontext
from dataclasses import dataclass
from pathlib import Path

# ---------------- 默认参数 (CLI / pipeline_config.yaml 可覆盖) ----------------
WINDOW = 50000        # 滑窗宽度 (nt)
EVALUE = "1e-5"       # blastx e-value 阈值
MERGE_D = 300         # 位点合并距离 (nt)
HOST_BS = 50          # Stage2 双侧判定阈值 (bitscore)
VIRAL_BS = 50

OUTFMT1 = ("6 qseqid qlen sseqid slen stitle evalue bitscore "
           "qcovhsp pident length qstart qend")
OUTFMT2 = ("6 qseqid sseqid pident length qstart qend sstart send "
           "evalue bitscore stitle")

# stitle 关键词 -> 科级归属 (轻量 ICTV 映射, 仅用于 ref_family 粗分类列)
# 已知局限 (有意保留): 只匹配属/科级词根, 种名级同义词不展开 ——
#   Tobacco mosaic virus (Tobamovirus/Virgaviridae)、
#   Bean yellow dwarf virus (Babuvirus/Nanoviridae) 等落 Other_viral;
#   不补 "yellow dwarf": Barley yellow dwarf virus 属 Luteoviridae, 会误伤。
# 精确分类由 Stage3 RVDB 注释 + 下游分类工具负责, 本列仅供初筛聚类参考。
FAMILY_KEYWORDS = [
    # 注意: "cauliflower" 必须单独列 — "caulimo" 匹配不到 Cauliflower mosaic virus
    # (cauli-f-lower 无 'm'), 原 eve_kingdom.py 关键词表漏此项导致 CaMV 落 Other_viral
    ("cauliflower", "Caulimoviridae"), ("caulimo", "Caulimoviridae"),
    ("badna", "Caulimoviridae"),
    ("soymo", "Caulimoviridae"), ("tungro", "Caulimoviridae"),
    ("cassava vein", "Caulimoviridae"), ("vein clearing", "Caulimoviridae"),
    ("pararetro", "Caulimoviridae"), ("florendo", "Caulimoviridae(endo)"),
    ("yendo", "Caulimoviridae(endo)"), ("zendo", "Caulimoviridae(endo)"),
    ("gemin", "Geminiviridae"), ("curto", "Geminiviridae"),
    ("becurto", "Geminiviridae"), ("nano", "Nanoviridae"),
    ("mimivir", "NCLDV"), ("megavir", "NCLDV"), ("pandora", "NCLDV"),
    ("phycodna", "NCLDV"), ("chlorovir", "NCLDV"),
    ("irido", "Iridoviridae"), ("afri", "NCLDV"),
    ("partiti", "Partitiviridae"), ("chryso", "Chrysoviridae"),
    ("endorna", "Endornaviridae"), ("ourmia", "Ourmiavirus"),
    ("virga", "Virgaviridae"), ("tobamo", "Virgaviridae"),
    ("tombus", "Tombusviridae"), ("nodo", "Nodaviridae"),
    ("picorna", "Picornavirales"), ("sobemo", "Solemoviridae"),
    ("soil-borne", "Solemoviridae"), ("bromo", "Bromoviridae"),
    ("cucumo", "Bromoviridae"), ("ilar", "Bromoviridae"),
    ("clostero", "Closteroviridae"), ("crinivir", "Closteroviridae"),
    ("poty", "Potyviridae"), ("maclura", "Potyviridae"),
    ("potex", "Alphaflexiviridae"), ("carla", "Betaflexiviridae"),
    ("tricho", "Trichoviridae"), ("beny", "Benyviridae"),
    ("furo", "Furoviridae"), ("oregenvirus", "Nodaviridae"),
    ("quadrivi", "Quadriviridae"), ("amalgavir", "Amalgaviridae"),
    ("reovir", "Sedoreoviridae"), ("cytorhabdo", "Rhabdoviridae"),
    ("nuclrhabdo", "Rhabdoviridae"), ("varicella", "Naldaviricetes"),
    ("nudivir", "Naldaviricetes"), ("baculo", "Naldaviricetes"),
    ("phage", "Phage"), ("mycovir", "Mycovirus"),
]

# 阶段 -> 子目录名 (MMPV 编号目录惯例; 目录名统一由 mmpv_common.io_layout
# 注册表给出 —— legacy/standard 两布局下 EVE 内部编号同名, 输出根约定:
# legacy=用户 -o 任意, standard=<项目>/06_EVE。导入失败时回退内置字面量,
# 保证本模块在脱离仓库环境时仍可独立运行。)
try:
    import sys as _sys
    _EVE_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    if _EVE_REPO_ROOT not in _sys.path:
        _sys.path.insert(0, _EVE_REPO_ROOT)
    from mmpv_common.io_layout import dir_name as _eve_layout_dir
    STAGE_DIRS = {k: _eve_layout_dir("e_" + k)
                  for k in ("loci", "verdict", "rvdb", "summary")}
except Exception:  # pragma: no cover - 脱离仓库环境的兜底
    STAGE_DIRS = {
        "loci":    "01_Loci",
        "verdict": "02_Verdict",
        "rvdb":    "03_RVDB",
        "summary": "04_Summary",
    }

SAFE_CHARS = set("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789._-")


@dataclass
class EveConfig:
    """EVE 筛查全部外部依赖与运行参数 (数据库/工具/阈值)."""
    # 数据库
    ref_dmnd: str = ""        # Stage1: NCBI 病毒参考蛋白 diamond 库
    pv_dmnd: str = ""         # Stage2: 植物+病毒拆分库 (SwissProt 混合)
    id2div: str = ""          # Stage2: sseqid -> viral/plant 归属表
    rvdb_dmnd: str = ""       # Stage3: RVDB diamond 库
    viroids_fa: str = ""      # 可选: 类病毒 blastn 库
    # 工具 (空 = 用 PATH)
    diamond: str = ""
    samtools: str = ""       # 抽序列用 samtools faidx (Stage1/3 都要)
    # 运行参数
    window: int = WINDOW
    merge_d: int = MERGE_D
    evalue: str = EVALUE
    host_bs: float = HOST_BS
    viral_bs: float = VIRAL_BS
    cmd_timeout: int = 0     # 单条外部命令超时 (秒); 0 = 不限 (防工具挂死占核)


def logmsg(msg):
    sys.stderr.write(f"[{time.strftime('%H:%M:%S')}] {msg}\n")
    sys.stderr.flush()


def clean_name(name):
    """组装名清洗: 仅保留安全字符, 防止派生路径逃逸"""
    cleaned = "".join(ch if ch in SAFE_CHARS else "_" for ch in str(name))
    if not cleaned or cleaned.startswith("."):
        cleaned = "asm_" + cleaned
    return cleaned


def safe_dir(p):
    """规范化目录路径 (abspath 后 '..' 已被消解, 无需再查)."""
    return Path(os.path.abspath(str(p)))


def tool_path(name, configured):
    """解析工具可执行文件: 配置路径优先, 否则 PATH 查找"""
    if configured:
        p = Path(os.path.expanduser(configured))
        if p.is_file():
            return str(p)
        raise SystemExit(f"配置的工具不存在: {name} -> {configured}")
    found = shutil.which(name)
    if not found:
        raise SystemExit(f"找不到工具 '{name}' (PATH 也没有); "
                         f"请安装或用参数显式指定路径")
    return found


SAMTOOLS_MIN = (1, 11)      # faidx 的 --fai-idx 从 1.11 起提供


def samtools_too_old(exe):
    """samtools 版本低于 SAMTOOLS_MIN 时返回版本串, 否则 None.

    探测本身失败 (跑不起来/版本行不认识) 一律返回 None — 不因为探测不了就拦下
    任务, 交给 samtools 自己的报错兜底.
    """
    try:
        r = subprocess.run([str(exe), "--version"], capture_output=True,
                           text=True, timeout=60)
        m = re.search(r"(?<![0-9.])(\d+)\.(\d+)", r.stdout or "")
        if not m:
            return None
        ver = (int(m.group(1)), int(m.group(2)))
        return None if ver >= SAMTOOLS_MIN else f"{ver[0]}.{ver[1]}"
    except Exception:
        return None


def check_tools(cfg, need_viroid=False):
    """运行前体检: 返回缺失/不可用的工具与数据库列表 (空 = 全部就绪).

    抽序列只用 samtools faidx (Stage1/3 都要), 所以 samtools 必查且版本要够
    (faidx 的 --fai-idx 从 1.11 起); 版本不够的报错要拦在开跑前 — 否则批处理
    跑到每个基因组的 Stage1 才炸, 几百个基因组白等几小时.
    need_viroid 时 blastn 必查, 否则类病毒层要跑到每基因组末尾才报错,
    批量长跑会白跑几小时.
    """
    missing = []
    for label, path in [("ref_dmnd", cfg.ref_dmnd), ("pv_dmnd", cfg.pv_dmnd),
                        ("id2div", cfg.id2div), ("rvdb_dmnd", cfg.rvdb_dmnd)]:
        if not path or not Path(os.path.expanduser(path)).exists():
            missing.append(f"{label}={path or '(未配置)'}")
    for label, exe in [("diamond", cfg.diamond or "diamond"),
                       ("seqkit", "seqkit"),
                       ("samtools", cfg.samtools or "samtools")]:
        p = Path(os.path.expanduser(exe))
        if not (p.is_file() or shutil.which(str(p))):
            missing.append(f"{label}={exe}")
            continue
        if label == "samtools":
            exe_path = str(p) if p.is_file() else shutil.which(str(p))
            old = samtools_too_old(exe_path)
            if old:
                missing.append(f"samtools={exe} (版本 {old}, 需 >= "
                               f"{SAMTOOLS_MIN[0]}.{SAMTOOLS_MIN[1]}: "
                               f"faidx 缺 --fai-idx)")
    if need_viroid:
        if not shutil.which("blastn"):
            missing.append("blastn=blastn")
        if not cfg.viroids_fa or not Path(os.path.expanduser(cfg.viroids_fa)).exists():
            missing.append(f"viroids_fa={cfg.viroids_fa or '(未配置)'}")
    return missing


def run(cmd, log, err_path=None, timeout=0):
    """列表参数执行, 不经 shell; 失败即终止. stderr/stdout 追加到 err_path 便于排查.

    timeout>0 时超时杀进程并报错 (防 diamond/seqkit 挂死白占并行核);
    timeout=0 不限.
    """
    log.info("RUN %s", " ".join(str(c) for c in cmd))
    with open(err_path, "a", encoding="utf-8") if err_path else nullcontext() as fh:
        try:
            r = subprocess.run([str(c) for c in cmd], stdout=fh,
                               stderr=subprocess.STDOUT, shell=False,
                               timeout=timeout or None)
        except subprocess.TimeoutExpired:
            if err_path and Path(err_path).is_file():
                with open(err_path, encoding="utf-8") as fh2:
                    tail = fh2.readlines()[-15:]
                sys.stderr.write("".join(tail) + "\n")
            raise RuntimeError(f"TIMEOUT ({timeout}s): {cmd[0]}")
    if r.returncode != 0:
        if err_path and Path(err_path).is_file():
            with open(err_path, encoding="utf-8") as fh:
                tail = fh.readlines()[-15:]
            sys.stderr.write("".join(tail) + "\n")
        raise RuntimeError(f"FAILED: {cmd[0]} (exit {r.returncode})")


def diamond_blastx(cfg, db, query, out, threads, sens, outfmt, log,
                   extra=None, err_path=None):
    cmd = [tool_path("diamond", cfg.diamond), "blastx", sens,
           "-F15", "--range-culling", "--max-hsps", "100",
           "-c1", "-b6", "-e", cfg.evalue, "-p", str(threads),
           "-d", db, "-q", query, "-o", out, "--outfmt"] + outfmt.split()
    if extra:
        cmd += [str(x) for x in extra]
    run(cmd, log, err_path=err_path, timeout=cfg.cmd_timeout)


# ---------------- 纯逻辑 (可单测): 坐标还原 / 位点合并 / 双侧判定 / 汇总 ----------------
def restore_coordinates(raw_tsv, hits_bed, window_note=""):
    """diamond blastx 命中 (query=滑窗块) 坐标还原回基因组 + 写 BED.

    raw_tsv 列 (OUTFMT1): qseqid qlen sseqid slen stitle evalue bitscore
                          qcovhsp pident length qstart qend
    qseqid 形如 <contig>_sliding:<块起点1based闭区间>-<块终点> —— seqkit sliding
    的产物, 两端都是 1-based 闭区间 (即 csvtk 的 :START-END 语义), 不是 0-based.
    qstart/qend 也是 1-based 闭区间. 两者都转成 BED 的 0-based 半开:
        gs = 块起点 + qstart - 1        (基因组上的 1-based 位置)
        BED start = gs - 1             BED end = ge  (= 1-based 闭区间终点)
    这里没有 strand 处理: OUTFMT1 不请求 qframe/sframe, 位点是无方向区间.
    输出 BED: contig \t start0 \t end \t <原始行>
    """
    n_hit = 0
    with open(raw_tsv, encoding="utf-8") as fh, open(hits_bed, "w", encoding="utf-8") as fo:
        for line in fh:
            c = line.rstrip("\n").split("\t")
            if len(c) < 12:
                continue
            qseqid = c[0]
            i = qseqid.rfind("_sliding:")
            if i < 0:
                continue
            contig = qseqid[:i]
            try:
                cs = int(qseqid[i + 9:].split("-")[0])
                s, e = int(c[10]), int(c[11])
            except ValueError:
                continue
            if s > e:
                s, e = e, s
            gs, ge = cs + s - 1, cs + e - 1
            fo.write(f"{contig}\t{gs - 1}\t{ge}\t{line.rstrip(chr(10))}\n")
            n_hit += 1
    return n_hit


def merge_loci(hits_bed, loci_bed, merge_d=MERGE_D):
    """区间合并: 同 contig 上间隔 <= merge_d 的命中并为一个位点."""
    ints = {}
    with open(hits_bed, encoding="utf-8") as fh:
        for line in fh:
            c = line.split("\t", 3)
            if len(c) < 3:
                continue
            try:
                s, e = int(c[1]), int(c[2])
            except ValueError:
                continue
            if e > s:
                ints.setdefault(c[0], []).append((s, e))
    n_loci = 0
    with open(loci_bed, "w", encoding="utf-8") as fo:
        for contig in sorted(ints):
            ivs = sorted(ints[contig])
            cur_s, cur_e = ivs[0]
            for s, e in ivs[1:]:
                if s <= cur_e + merge_d:
                    if e > cur_e:
                        cur_e = e
                else:
                    fo.write(f"{contig}\t{cur_s}\t{cur_e}\t"
                             f"{contig}:{cur_s + 1}-{cur_e}\n")
                    n_loci += 1
                    cur_s, cur_e = s, e
            fo.write(f"{contig}\t{cur_s}\t{cur_e}\t"
                     f"{contig}:{cur_s + 1}-{cur_e}\n")
            n_loci += 1
    return n_loci


def assign_best_hits(hits_bed, loci_bed, best_tsv):
    """每位点取最佳 (bitscore 最大) ref 命中. 双指针: 命中已按坐标排序."""
    loci = []
    with open(loci_bed, encoding="utf-8") as fi:
        for line in fi:
            c = line.rstrip("\n").split("\t")
            if len(c) >= 4:
                loci.append((c[0], int(c[1]), int(c[2]), c[3]))
    loci.sort(key=lambda x: (x[0], x[1]))
    hits2 = []
    with open(hits_bed, encoding="utf-8") as fh:
        for line in fh:
            c = line.rstrip("\n").split("\t")
            if len(c) < 15:
                continue
            try:
                hits2.append((c[0], int(c[1]), int(c[2]), c))
            except ValueError:
                continue
    hits2.sort(key=lambda x: (x[0], x[1], x[2]))
    rows = {}
    li, n_loci = 0, len(loci)
    for hc, hs, he, c in hits2:
        while li < n_loci and (loci[li][0] < hc or
                               (loci[li][0] == hc and loci[li][2] <= hs)):
            li += 1
        if li >= n_loci:
            break
        lc, ls, le, lid = loci[li]
        if lc == hc and hs < le and he > ls:
            try:
                bs = float(c[9])
            except ValueError:
                continue
            if lid not in rows or bs > rows[lid][0]:
                rows[lid] = (bs, c[5], family_of(c[7]),
                             float(c[8]), float(c[9]), float(c[10]))
    with open(best_tsv, "w", encoding="utf-8") as fo:
        fo.write("locus\tref_sseqid\tfamily\tref_evalue\tref_bitscore\tref_qcov\n")
        for lid in sorted(rows):
            v = rows[lid]
            fo.write(f"{lid}\t{v[1]}\t{v[2]}\t{v[3]}\t{fmt_num(v[4])}\t"
                     f"{fmt_num(v[5])}\n")
    return len(rows)


def family_of(stitle):
    low = stitle.lower()
    for kw, fam in FAMILY_KEYWORDS:
        if kw in low:
            return fam
    return "Other_viral"


def fmt_num(x):
    """bitscore/覆盖度等数值输出: 整数去 .0 (表格更干净)."""
    f = float(x)
    return str(int(f)) if f.is_integer() else f"{f:g}"


def locus_key(name):
    """位点名归一化 (仅用于跨表 join, 不改变输出): loci.bed 第4列写作
    chr1:50-699, 抽序列工具写 FASTA 头、diamond 回显 qseqid 时可能把 ':'
    换成 '_' (samtools faidx 原样写区域名, 但历史上不同版本/工具行为不保证).
    两侧都过这个函数后 join, 否则 Stage2 判定会被静默丢弃 (verdict/summary
    全变 undetermined, Stage3 候选为 0)."""
    return name.replace(":", "_")


def load_id2div(path):
    """id2div 表: sseqid -> 'viral' / 'plant'"""
    m = {}
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            c = line.rstrip("\n").split("\t")
            if len(c) == 2:
                m[c[0].split()[0]] = c[1]
    return m


def verdict_loci(raw_tsv, id2div_map, verdict_tsv, host_bs=HOST_BS, viral_bs=VIRAL_BS):
    """Stage2 双侧判定: 拆分库内 viral 侧最高 bitscore vs plant 侧最高 bitscore.

    返回 (n_viral_supported, n_host_like, n_undetermined).
    verdict: viral_supported = 病毒侧过阈值且不弱于植物侧 (排除宿主基因假阳性)
             host_like        = 植物侧更强 (宿主基因/重复序列, 非 EVE)
             undetermined     = 两侧都不过阈值 (需 Stage3 深度归属)
    """
    hits = {}
    with open(raw_tsv, encoding="utf-8") as fh:
        for line in fh:
            c = line.rstrip("\n").split("\t")
            if len(c) < 11:
                continue
            try:
                hits.setdefault(c[0], []).append((float(c[9]), c[1]))
            except ValueError:
                continue
    n_v = n_h = n_u = 0
    unmapped = set()
    with open(verdict_tsv, "w", encoding="utf-8") as fo:
        fo.write("locus\tverdict\tbest_viral_id\tbest_viral_bs\t"
                 "best_plant_id\tbest_plant_bs\n")
        for locus in sorted(hits):
            bv_id, bv = "", 0.0
            bp_id, bp = "", 0.0
            for bs, sid in hits[locus]:
                dv = id2div_map.get(sid, "?")
                if dv == "?":
                    unmapped.add(sid)
                if dv == "viral" and bs > bv:
                    bv_id, bv = sid, bs
                elif dv == "plant" and bs > bp:
                    bp_id, bp = sid, bs
            if bv >= viral_bs and bv >= bp:
                v = "viral_supported"
                n_v += 1
            elif bp >= host_bs and bp > bv:
                v = "host_like"
                n_h += 1
            else:
                v = "undetermined"
                n_u += 1
            fo.write(f"{locus}\t{v}\t{bv_id}\t{fmt_num(bv)}\t"
                     f"{bp_id}\t{fmt_num(bp)}\n")
    # 未映射的 sseqid 一律不参与两侧对比: 若 id2div 表分隔符不对或与 pv_dmnd 版本
    # 不匹配, 全部位点会静默落成 undetermined / Stage3 候选为 0 而退出码仍是 0.
    # 这里把它喊出来 —— 阈值内的正常情况是 0 条.
    if unmapped:
        logmsg(f"WARN stage2: {len(unmapped)} 个 sseqid 不在 id2div 表里 "
               f"(最多列出 3 个: {sorted(unmapped)[:3]}); 这些命中不参与 "
               f"viral/plant 两侧对比, 检查 id2div 与 pv_dmnd 是否配套")
    return n_v, n_h, n_u


def summarize_genome(name, outdir):
    """单基因组: s1_best + s2_verdict + s3_rvdb 合并为 _eve_summary.tsv. 返回路径."""
    loci_d = Path(outdir) / STAGE_DIRS["loci"] / name
    verd_d = Path(outdir) / STAGE_DIRS["verdict"] / name
    rvdb_d = Path(outdir) / STAGE_DIRS["rvdb"] / name
    summ_d = Path(outdir) / STAGE_DIRS["summary"]
    summ_d.mkdir(parents=True, exist_ok=True)
    summ = summ_d / f"{name}_eve_summary.tsv"

    best = {}
    bp = loci_d / f"{name}.s1_best.tsv"
    if bp.exists():
        with open(bp, encoding="utf-8") as fh:
            fh.readline()
            for line in fh:
                c = line.rstrip("\n").split("\t")
                if c:
                    best[locus_key(c[0])] = c
    verd = {}
    vp = verd_d / f"{name}.s2_verdict.tsv"
    if vp.exists():
        with open(vp, encoding="utf-8") as fh:
            fh.readline()
            for line in fh:
                c = line.rstrip("\n").split("\t")
                if c:
                    verd[locus_key(c[0])] = c
    rvdb = {}
    s3 = rvdb_d / f"{name}.s3_rvdb.tsv"
    if s3.exists():
        with open(s3, encoding="utf-8") as fh:
            for line in fh:
                c = line.rstrip("\n").split("\t")
                if len(c) < 11:
                    continue
                try:
                    bs = float(c[9])
                except ValueError:
                    continue
                k = locus_key(c[0])
                if k not in rvdb or bs > rvdb[k][0]:
                    rvdb[k] = (bs, c[10] if len(c) > 10 else "")
    with open(loci_d / f"{name}.loci.bed", encoding="utf-8") as fi, open(summ, "w", encoding="utf-8") as fo:
        # 12 列, 不含 genome —— 基因组名在文件名与目录名里, 汇总时由 merge_all 补首列.
        # 这里多写一列会让 merge_all 拼出 14 字段而表头只有 13 个名字, 按列名读的
        # 消费者 (pandas / csv.DictReader) 每一列都右移一位.
        fo.write("locus\tverdict\tref_family\tref_sseqid\tref_bitscore\t"
                 "ref_qcov\tbest_viral_sp\tbest_viral_bs\tbest_plant_sp\t"
                 "best_plant_bs\trvdb_stitle\trvdb_bitscore\n")
        for line in fi:
            c = line.rstrip("\n").split("\t")
            if len(c) < 4:
                continue
            locus = c[3]
            lk = locus_key(locus)
            b = best.get(lk, [""] * 7)
            v = verd.get(lk, ["", "undetermined", "", "", "", ""])
            r = rvdb.get(lk, ("", ""))
            fo.write(f"{locus}\t{v[1]}\t"
                     f"{b[2] if len(b) > 2 else ''}\t{b[1] if len(b) > 1 else ''}\t"
                     f"{b[4] if len(b) > 4 else ''}\t{b[5] if len(b) > 5 else ''}\t"
                     f"{v[2]}\t{v[3]}\t{v[4]}\t{v[5]}\t{r[1]}\t"
                     f"{fmt_num(r[0]) if r[0] != '' else ''}\n")
    return summ


def merge_all(outdir):
    """跨基因组汇总: kingdom_summary / kingdom_loci_all / family_by_genome."""
    outdir = Path(outdir)
    t0 = time.time()
    all_rows = []
    counts = {}
    fam_g = {}
    loci_root = outdir / STAGE_DIRS["loci"]
    if not loci_root.is_dir():
        raise SystemExit(f"无 {STAGE_DIRS['loci']}/ 目录, 先跑筛查阶段: {outdir}")
    for d in sorted(loci_root.iterdir()):
        if not d.is_dir():
            continue
        name = d.name
        s = outdir / STAGE_DIRS["summary"] / f"{name}_eve_summary.tsv"
        if not s.exists():
            continue
        n_v = n_h = n_u = 0
        with open(s, encoding="utf-8") as fh:
            fh.readline()
            for line in fh:
                c = line.rstrip("\n").split("\t")
                if len(c) < 3:
                    continue
                if len(c) == 12 + 1:
                    # 旧布局 (13 列, 首列是重复的 genome) 的 summary: 丢掉首列,
                    # 否则按 12 列索引会整行错位. 新写的 summary 恒为 12 列.
                    # 只认这一个宽度 —— 别的宽度是别的东西 (截断行, 或 eve_kingdom.py
                    # 那套 15 列的旧格式), 一律交给下面的补齐/告警, 不要静默丢首列.
                    c = c[1:]
                    line = "\t".join(c) + "\n"
                all_rows.append((name, line.rstrip("\n")))
                if c[1] == "viral_supported":
                    n_v += 1
                elif c[1] == "host_like":
                    n_h += 1
                else:
                    n_u += 1
                fam = c[2] or "Other_viral"
                if c[1] == "viral_supported":
                    fam_g.setdefault(fam, {})
                    fam_g[fam][name] = fam_g[fam].get(name, 0) + 1
        counts[name] = (n_v, n_h, n_u)
    with open(outdir / "kingdom_summary.tsv", "w", encoding="utf-8") as fo:
        fo.write("genome\tviral_supported\thost_like\tundetermined\n")
        for g in sorted(counts):
            fo.write(f"{g}\t{counts[g][0]}\t{counts[g][1]}\t{counts[g][2]}\n")
    loci_hdr = ("genome\tlocus\tverdict\tref_family\tref_sseqid\tref_bitscore\t"
                "ref_qcov\tbest_viral_sp\tbest_viral_bs\tbest_plant_sp\t"
                "best_plant_bs\trvdb_stitle\trvdb_bitscore")
    width = loci_hdr.count("\t") + 1
    n_ragged = 0
    with gzip.open(outdir / "kingdom_loci_all.tsv.gz", "wt", encoding="utf-8") as z:
        z.write(loci_hdr + "\n")
        for g_, row in all_rows:
            c = row.split("\t")
            if len(c) != width - 1:
                # 单基因组 summary 被截断/多列: 补空或截尾到 12 列, 保证汇总表恒为
                # 矩形 (列名对齐才读得对), 同时计数告警而不是静默错位.
                n_ragged += 1
                c = (c + [""] * 12)[:12]
            z.write(g_ + "\t" + "\t".join(c) + "\n")
    if n_ragged:
        logmsg(f"WARN merge: {n_ragged} 行 summary 列数不是 12, 已补齐/截尾; "
               f"检查对应 <NAME>_eve_summary.tsv 是否被截断")
    with open(outdir / "family_by_genome.tsv", "w", encoding="utf-8") as fo:
        genomes = sorted(counts)
        fo.write("family\t" + "\t".join(genomes) + "\n")
        for fam in sorted(fam_g):
            fo.write(fam + "\t" + "\t".join(
                str(fam_g[fam].get(g, 0)) for g in genomes) + "\n")
    return len(counts), len(all_rows)


# ---------------- 阶段产物路径 ----------------
def genome_dir(outdir, name, stage):
    """该基因组在某阶段的产物目录 (自动创建)."""
    d = Path(outdir) / STAGE_DIRS[stage] / name
    d.mkdir(parents=True, exist_ok=True)
    return d


def loci_paths(outdir, name):
    d = Path(outdir) / STAGE_DIRS["loci"] / name
    return {
        "dir": d,
        "fna": d / f"{name}.fna",
        "chunks": d / f"{name}.chunks.fna",
        "raw": d / f"{name}.s1_raw.tsv",
        "hits": d / f"{name}.s1_hits.bed",
        "loci_bed": d / f"{name}.loci.bed",
        "loci_fa": d / f"{name}.loci.fa",
        "fai": d / f"{name}.genome.fai",
        "best": d / f"{name}.s1_best.tsv",
    }


# ---------------- Stage 1 ----------------
def sliding_cmd(cfg, genome_fa, chunks, step):
    """seqkit sliding 命令 (切块为非重叠滑窗; overlap 时 step 减半).

    -g 是布尔开关 --greedy (导出每个 contig 末尾不足 window 的碎块): 必须写
    成长选项 --greedy, 不能把基因组路径塞给它当值 — 那样路径会被 pflag 漏成
    位置参数, 命令碰巧能跑但意图全丢, 且路径以 '-' 开头即失败.
    基因组为位置参数 (seqkit 从 stdin 或位置参数读输入).
    """
    return ["seqkit", "sliding", "-s", str(step), "-W", str(cfg.window),
            "--greedy", "-o", str(chunks), str(genome_fa)]


def stage1_discover(cfg, genome, name, outdir, threads, fast, overlap,
                    viroid, log, force=False):
    """正向发现. 返回 (基因组 fasta 路径, 阶段目录)."""
    d = genome_dir(outdir, name, "loci")
    p = loci_paths(outdir, name)
    t0 = time.time()
    if force:
        for k in ("chunks", "raw", "hits", "loci_bed", "loci_fa", "fai", "best"):
            p[k].unlink(missing_ok=True)

    # 0) 解压/解包 (支持 .tar.gz 包 / .gz / 裸 fasta)
    g = Path(genome)
    if str(genome).endswith(".tar.gz"):
        if force or not p["fna"].exists():
            log.info("extract genome from tar.gz")
            tdir = d / "_tar"
            tdir.mkdir(exist_ok=True)
            run(["tar", "-xzf", str(genome), "-C", str(tdir)], log,
                timeout=cfg.cmd_timeout)
            cands = [pp for pp in tdir.rglob("*") if pp.is_file()
                     and pp.suffix.lower() in (".fa", ".fasta", ".fna",
                                               ".fa.gz", ".fasta.gz", ".fna.gz")]
            if not cands:
                raise RuntimeError(f"tar.gz 内无 fasta: {genome}")

            def member_pri(pp):
                low = pp.name.lower()
                if "genome" in low:
                    return 0
                if "assembly" in low:
                    return 1
                return 2

            picked = None
            for pp in sorted(cands, key=lambda pp: (member_pri(pp),
                                                    -pp.stat().st_size)):
                if str(pp).endswith(".gz"):
                    probe = subprocess.run(
                        ["gzip", "-cd", str(pp)], capture_output=True,
                        timeout=cfg.cmd_timeout or None).stdout[:5000]
                else:
                    with open(pp, "rb") as fh:
                        probe = fh.read(5000)
                s = probe.upper()
                good = sum(1 for ch in s if ch in b"ACGTN")
                if len(s) > 100 and good / len(s) > 0.9:
                    picked = pp
                    break
            if picked is None:
                raise RuntimeError(
                    f"tar.gz 内无 DNA 基因组 fasta (成员可能是蛋白/CDS 错标): {genome}")
            if str(picked).endswith(".gz"):
                _gzip_to(picked, p["fna"], timeout=cfg.cmd_timeout)
            else:
                shutil.move(str(picked), str(p["fna"]))
            shutil.rmtree(tdir)
        g = p["fna"]
    elif str(genome).endswith(".gz"):
        if force or not p["fna"].exists():
            log.info("decompress genome")
            _gzip_to(genome, p["fna"], timeout=cfg.cmd_timeout)
        g = p["fna"]
    if not Path(g).is_file():
        raise FileNotFoundError(f"基因组文件不存在: {genome}")

    # 1) 切块
    step = cfg.window // 2 if overlap else cfg.window
    if force or not p["chunks"].exists():
        log.info("seqkit sliding W=%d step=%d", cfg.window, step)
        run(sliding_cmd(cfg, g, p["chunks"], step), log,
            err_path=d / "seqkit.err", timeout=cfg.cmd_timeout)
    else:
        log.info("chunks exist, skip")

    # 2) 正向 blastx vs ref
    if force or not p["raw"].exists():
        sens = "--fast" if fast else "--sensitive"
        log.info("stage1 diamond blastx vs ref (%s)", sens)
        diamond_blastx(cfg, cfg.ref_dmnd, str(p["chunks"]), str(p["raw"]),
                       threads, sens, OUTFMT1, log,
                       err_path=d / "diamond_s1.err")
    else:
        log.info("stage1 raw hits exist, skip")

    # 3) 坐标还原
    if force or not p["hits"].exists():
        n_hit = restore_coordinates(p["raw"], p["hits"])
        log.info("hits: %d", n_hit)

    # 4) 位点合并
    if force or not p["loci_bed"].exists():
        n_loci = merge_loci(p["hits"], p["loci_bed"], cfg.merge_d)
        log.info("loci: %d", n_loci)

    # 5) 抽序列 (samtools faidx)
    if force or not p["loci_fa"].exists():
        _extract(cfg, g, p["loci_bed"], p["loci_fa"], p["fai"], log)

    # 6) 每位点最佳 ref 命中
    if force or not p["best"].exists():
        n_best = assign_best_hits(p["hits"], p["loci_bed"], p["best"])
        log.info("best-hit table: %d", n_best)

    log.info("STAGE1 COMPLETE in %.0fs", time.time() - t0)
    return g, d


def _gzip_to(src, dst, timeout=0):
    """gzip -cd src > dst (子进程, 不依赖 Python gzip 模块解大文件)."""
    with open(src, "rb") as fi, open(dst, "wb") as fo:
        subprocess.run(["gzip", "-cd"], stdin=fi, stdout=fo, check=True,
                       timeout=timeout or None)


def bed_to_regions(bed, regions):
    """BED (0 基半开) → samtools 区域文件 (1 基闭区间, 一行一个 chr:FROM-TO).

    位点名 (第 4 列) 不使用: samtools faidx 以区域字符串本身作输出 FASTA 头,
    恰好就是本管线的位点命名 (chr1:50-699); 但下游仍统一过 locus_key()
    归一化, 不依赖这个巧合.
    返回区域行数; 注释/列数不足/起止倒挂的行跳过 (与工具对空 bed 的容忍一致).
    """
    n = 0
    with open(bed, encoding="utf-8") as fi, \
            open(regions, "w", encoding="utf-8") as fo:
        for line in fi:
            s = line.strip()
            if not s or s.startswith(("#", "track", "browser")):
                continue
            c = s.split("\t")
            if len(c) < 3:
                continue
            try:
                start, end = int(c[1]), int(c[2])
            except ValueError:
                continue
            if end <= start:
                continue
            fo.write(f"{c[0]}:{start + 1}-{end}\n")
            n += 1
    return n


def extract_cmd(cfg, genome_fa, fai, regions, out_fa):
    """samtools faidx 按区域文件抽序列.

    --fai-idx 指向输出目录里的索引: 基因组可能是共享只读路径 (如 /db 下的参考),
    samtools 默认把 .fai 写在参考旁边会失败/污染数据库目录.
    选项名必须是 --fai-idx 而不是 --fai: 后者在部分版本上能被无歧义前缀
    匹配到, 但那是未文档化的解析行为, 不能依赖.
    --region-file 而不是把区域逐个当参数: 大基因组位点数以万计,
    命令行长度会撞 ARG_MAX.
    需要 samtools >= 1.11 (--fai-idx 在 1.11 才有; --region-file/--output 1.10 起).
    """
    return [tool_path("samtools", cfg.samtools), "faidx", str(genome_fa),
            "--fai-idx", str(fai), "--region-file", str(regions),
            "--output", str(out_fa)]


def check_extracted(regions, out_fa):
    """核对抽序列产物: 区域文件里每个区域, 输出里都得有一条等长记录.

    samtools faidx 遇到越界区域不报错: 退出码仍是 0, 只在 stderr 写一行
    [faidx] Zero length / Truncated sequence, 而输出里留下空序列或少一截的
    序列 (已实测 1.21). 位点坐标是 restore_coordinates 从滑窗坐标还原来的,
    还原一错, 症状正好落在这里 — 不核对, 脏序列会一路静默喂给 Stage2/3
    和汇总表. run() 只看退出码兜不住, 所以单独查一遍.

    记录的归名统一过 locus_key(): 下游跨表 join 本就不依赖头是 ':' 还是 '_',
    核对用同一把尺子, 免得换个工具写法就误报.
    """
    want = [l.strip() for l in open(regions, encoding="utf-8") if l.strip()]
    got, hdr = {}, None
    for line in open(out_fa, encoding="utf-8"):
        if line.startswith(">"):
            hdr = line[1:].strip()
            got[locus_key(hdr)] = 0
        elif hdr is not None:
            got[locus_key(hdr)] += len(line.strip())
    bad = []
    for r in want:
        s, e = r.rsplit(":", 1)[1].split("-")
        want_len = int(e) - int(s) + 1
        if locus_key(r) not in got:
            bad.append(f"{r} 无记录")
        elif got[locus_key(r)] != want_len:
            bad.append(f"{r} {got[locus_key(r)]}bp != {want_len}bp")
    if bad:
        raise RuntimeError(f"抽序列结果与区域不符 ({len(bad)}/{len(want)} 条): "
                           + "; ".join(bad[:5])
                           + (" ..." if len(bad) > 5 else ""))
    return len(want)


def _extract(cfg, genome_fa, bed, out_fa, fai, log):
    """按 BED 抽序列 (samtools faidx).

    空 bed 时直接写空文件, 不去调 samtools — 空区域文件各版本行为不一,
    且没必要为一个空结果去建整个基因组的索引.
    """
    regions = Path(str(out_fa) + ".regions.txt")
    n = bed_to_regions(bed, regions)
    if n == 0:
        log.info("无位点, 写空序列文件")
        Path(out_fa).write_text("", encoding="utf-8")
        regions.unlink(missing_ok=True)
        return 0
    run(extract_cmd(cfg, genome_fa, fai, regions, out_fa), log,
        timeout=cfg.cmd_timeout)
    check_extracted(regions, out_fa)
    return n


# ---------------- Stage 2 ----------------
def stage2_verdict(cfg, name, outdir, threads, log, force=False):
    d = genome_dir(outdir, name, "verdict")
    p = loci_paths(outdir, name)
    t0 = time.time()
    raw = d / f"{name}.s2_raw.tsv"
    verdict = d / f"{name}.s2_verdict.tsv"
    if force:
        for f in (raw, verdict):
            f.unlink(missing_ok=True)

    if force or not verdict.exists():
        loci_fa = p["loci_fa"]
        if not loci_fa.is_file() or loci_fa.stat().st_size == 0:
            log.info("无 loci 序列, Stage2 输出空表")
            verdict.write_text("locus\tverdict\tbest_viral_id\tbest_viral_bs\t"
                               "best_plant_id\tbest_plant_bs\n")
            return d
        if force or not raw.exists():
            log.info("stage2 diamond blastx vs plant_virus (loci.fa)")
            diamond_blastx(cfg, cfg.pv_dmnd, str(loci_fa), str(raw),
                           threads, "--fast", OUTFMT2, log,
                           extra=["-k", "5"], err_path=d / "diamond_s2.err")
        else:
            log.info("stage2 raw hits exist, skip")
        log.info("load id2div: %s", cfg.id2div)
        div = load_id2div(cfg.id2div)
        n_v, n_h, n_u = verdict_loci(raw, div, verdict,
                                     cfg.host_bs, cfg.viral_bs)
        log.info("STAGE2 COMPLETE in %.0fs (viral_supported=%d host_like=%d "
                 "undetermined=%d)", time.time() - t0, n_v, n_h, n_u)
    else:
        log.info("stage2 verdict exists, skip")
    return d


# ---------------- Stage 3 ----------------
def stage3_annotate(cfg, name, outdir, threads, genome_fa, log, force=False):
    d = genome_dir(outdir, name, "rvdb")
    p = loci_paths(outdir, name)
    t0 = time.time()
    cand_fa = d / f"{name}.cand3.fa"
    cand_bed = d / f"{name}.cand3.bed"
    out = d / f"{name}.s3_rvdb.tsv"
    if force:
        for f in (cand_fa, cand_bed, out):
            f.unlink(missing_ok=True)
    if not force and out.exists():
        log.info("stage3 exists, skip")
        return d

    verdict_f = Path(outdir) / STAGE_DIRS["verdict"] / name / f"{name}.s2_verdict.tsv"
    keep = set()
    if verdict_f.exists():
        with open(verdict_f, encoding="utf-8") as fh:
            fh.readline()
            for line in fh:
                c = line.rstrip("\n").split("\t")
                if len(c) > 1 and c[1] in ("viral_supported", "undetermined"):
                    keep.add(locus_key(c[0]))
    n = 0
    if keep:
        with open(p["loci_bed"], encoding="utf-8") as fi, open(cand_bed, "w", encoding="utf-8") as fo:
            for line in fi:
                c = line.rstrip("\n").split("\t")
                if len(c) > 3 and locus_key(c[3]) in keep:
                    fo.write(line)
                    n += 1
    if n == 0:
        log.info("无 Stage3 候选 (viral_supported + undetermined = 0), 输出空表")
        out.write_text("")
        return d
    _extract(cfg, genome_fa, cand_bed, cand_fa, p["fai"], log)
    log.info("stage3 diamond blastx vs RVDB (%d candidates)", n)
    diamond_blastx(cfg, cfg.rvdb_dmnd, str(cand_fa), str(out),
                   threads, "--fast", OUTFMT2, log,
                   extra=["-k", "5"], err_path=d / "diamond_s3.err")
    log.info("STAGE3 COMPLETE in %.0fs", time.time() - t0)
    return d


# ---------------- 可选类病毒层 ----------------
def stage_viroid(cfg, name, outdir, genome_fa, threads, log, force=False):
    """类病毒 blastn 层 (默认关): viroids.fa vs 基因组."""
    out = Path(outdir) / STAGE_DIRS["summary"] / f"{name}_viroid.tsv"
    out.parent.mkdir(parents=True, exist_ok=True)
    if not force and out.exists():
        log.info("viroid layer exists, skip")
        return out
    log.info("blastn viroids vs genome")
    run(["blastn", "-task", "blastn", "-word_size", "7", "-dust", "no",
         "-query", cfg.viroids_fa, "-subject", str(genome_fa),
         "-num_threads", str(min(int(threads), 4)),
         "-evalue", "1e-4", "-outfmt",
         "6 qseqid qlen sseqid slen pident length qstart qend sstart send evalue bitscore qcovs",
         "-out", str(out)], log, timeout=cfg.cmd_timeout)
    return out


# ---------------- 测序污染筛查 ----------------
# 为什么必须有这一步: 候选里混进 phiX174 (Illumina 建库的 spike-in) 或其它常见细菌/载体
# 序列时, 它们会打中病毒参考库里的噬菌体条目, 然后被当成"病毒信号"一路走到底 ——
# 在千基因组级别的汇总里, 这类位点能占到两成以上 (上游研究实测 Phage 层 11,855/41,887
# 元件, 占 28%, 其质控记录直接写明"Phage 类含 phiX 污染")。我们的管线此前**没有任何
# 污染筛查**, 所以这一层是补上的缺口。
#
# 判定方式与上游研究一致, 用**参考 accession 前缀 + 标题关键词**: 不要求再跑一次比对,
# 只对已有命中做判定, 因此可以事后对任何一次跑完的结果追加筛查。
CONTAM_ACC_PREFIX = (
    "NP_0406", "NP_0407",     # phiX174 蛋白 (NCBI 编号段固定)
    "YP_51237",               # Escherichia phage ID2, phiX 近缘
)
CONTAM_TITLE_KEYWORDS = (
    "phix", "phi x", "escherichia phage",
)
# 输出表头: 测试直接 import CONTAM_HDR 引用, 避免两头各抄一份
CONTAM_HDR = ("query_id\tn_hits\tn_contam_hits\tcontam_frac\ttop_contam_ref\t"
              "top_contam_title\n")


def is_contaminant_hit(sseqid, stitle=""):
    """这条 hit 的参考序列是不是已知测序污染源 (phiX/大肠杆菌噬菌体)."""
    sid = (sseqid or "").strip()
    if any(sid.startswith(p) for p in CONTAM_ACC_PREFIX):
        return True
    low = (stitle or "").lower()
    return any(k in low for k in CONTAM_TITLE_KEYWORDS)


def contamination_scan(raw_tsv, out_tsv):
    """逐 query 统计污染命中占比 -> CONTAM_HDR 表; 返回 (n_query, n_contaminated).

    `contam_frac` 是**该 query 的全部命中里**污染命中占的比例 —— 单看命中数会把
    "打中一次 phiX 的长 contig" 和 "整条都是 phiX" 混为一谈, 占比才分得开。
    下游按 frac 阈值 (配合 n_contam_hits) 过滤; 本函数不下结论, 只报事实。
    """
    tot = {}
    contam = {}
    top = {}
    with open(raw_tsv, encoding="utf-8") as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < 2 or p[0] == "qseqid":
                continue
            q = p[0]
            tot[q] = tot.get(q, 0) + 1
            stitle = p[11] if len(p) > 11 else ""
            if is_contaminant_hit(p[1], stitle):
                contam[q] = contam.get(q, 0) + 1
                # 记录该 query 的污染命中最强的那个参考 (按 bitscore)
                try:
                    bs = float(p[9])
                except (ValueError, IndexError):
                    bs = 0.0
                if q not in top or bs > top[q][0]:
                    top[q] = (bs, p[1], stitle)
    n_contam = 0
    with open(out_tsv, "w", encoding="utf-8") as fo:
        fo.write(CONTAM_HDR)
        for q in sorted(tot):
            nc = contam.get(q, 0)
            if not nc:
                continue
            n_contam += 1
            bs, sid, st = top[q]
            fo.write("%s\t%d\t%d\t%.3f\t%s\t%s\n"
                     % (q, tot[q], nc, nc / tot[q], sid, st[:80]))
    return len(tot), n_contam


# ---------------- 中间文件清理 ----------------
# (victim_stem, guard_stem): guard 存在才允许删 victim (终表永不被删).
# 两侧都必须查它真实的扩展名与阶段目录 —— 早期版本把护栏扩展名硬编码成 .tsv,
# 于是 s1_hits/loci_bed 两条 .bed 护栏永远不存在, 那两条规则从不生效.
_CLEANUP_STEMS = (
    ("chunks", "s1_raw"),
    ("s1_raw", "s1_hits"),
    ("s1_hits", "loci_bed"),
    ("s2_raw", "s2_verdict"),
    ("cand3_fa", "s3_rvdb"),
    ("cand3_bed", "s3_rvdb"),
)
# stem -> (磁盘文件名主体, 扩展名, 所在阶段目录).
# 主体未必等于这里的键: loci.bed 的主体是 "loci" (不是 "loci_bed"), cand3.fa 的主体是
# "cand3" (不是 "cand3_fa"). 早期版本按键直接拼文件名, 且护栏扩展名写死成 .tsv,
# 于是六条规则里有四条指向了根本不存在的路径 —— 清理看着"成功"却什么都没删.
_CLEANUP_META = {
    "chunks":     ("chunks", "fna", "loci"),
    "s1_raw":     ("s1_raw", "tsv", "loci"),
    "s1_hits":    ("s1_hits", "bed", "loci"),
    "loci_bed":   ("loci", "bed", "loci"),
    "s2_raw":     ("s2_raw", "tsv", "verdict"),
    "s2_verdict": ("s2_verdict", "tsv", "verdict"),
    "cand3_fa":   ("cand3", "fa", "rvdb"),
    "cand3_bed":  ("cand3", "bed", "rvdb"),
    "s3_rvdb":    ("s3_rvdb", "tsv", "rvdb"),
}


def cleanup_stage_files(outdir, name):
    """删除已被下游产物取代的中间文件; 终表与 loci.bed/loci.fa 永不删除.

    删除后细粒度断点退化为从 diamond 重跑 (文件不存在即重跑该步).
    返回释放字节数.
    """
    freed = 0
    for victim_stem, guard_stem in _CLEANUP_STEMS:
        v_stem, v_ext, v_stage = _CLEANUP_META[victim_stem]
        g_stem, g_ext, g_stage = _CLEANUP_META[guard_stem]
        victim = (Path(outdir) / STAGE_DIRS[v_stage] / name
                  / f"{name}.{v_stem}.{v_ext}")
        guard = (Path(outdir) / STAGE_DIRS[g_stage] / name
                 / f"{name}.{g_stem}.{g_ext}")
        if victim.exists() and guard.exists():
            try:
                freed += victim.stat().st_size
                victim.unlink()
            except OSError:
                pass
    return freed

