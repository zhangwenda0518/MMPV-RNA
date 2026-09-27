#!/usr/bin/env python3
"""gene_partition_dating.py — 分基因分区 BEAST 定年一键入口 (MMPV-RNA)

把 GCVA/PSTVd 型的「按基因切分 → 逐基因 BEAST skyline → 后台并行 → 自动
解析 MRCA/HPD/拓扑一致性」编排成一条命令，吸收进 virome_phylo_pipeline。

用法:
  python gene_partition_dating.py --virus GCVA \
      --fasta gcva_135.mafft.fasta --metadata phylo_dates_135.csv \
      --gff gcva.gff3 --out gene_beast \
      --chain 15000000 --big-genes L --chain-big 10000000 \
      --alpha-lower 0.5 --threads 24

设计要点（踩坑总结固化）:
  * alpha lower 参数化下发（默认 0.5），规避 BEAST GammaDistribution.pointChi2
    在 df=2*shape 过小时的 `t<0` 崩溃（dr.math.distributions.GammaDistribution）
  * 每条 XML 生成后强制 xml.dom.minidom 语法校验 + 校验 alpha 下界已生效
  * 大基因（如 L）支持独立/更短的 chain 长度；其余用默认 chain
  * 后台启动用 setsid+nohup 可靠 detach；进程检测用对 xml 目录名的模式，
    避免 pkill/grep 匹配到自身命令行（历史教训）
  * watchdog 用「每基因各自完成阈值」判定，逐基因完成即增量解析
"""
import argparse, os, re, shlex, shutil, subprocess, sys, time
from pathlib import Path

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'utils'))

from utils.beast1_bridge import generate_beast1_phylogeo_xml
from utils.virphy_bridge import LogCollector

DEFAULT_ALPHA_LOWER = "0.5"
DEFAULT_CHAIN = 15_000_000


# ─────────────────────────── 坐标/基因分析 ───────────────────────────
def read_gene_coords_gff(gff_path):
    """从 gff3 提取 gene feature 坐标，返回 {gene: (start,end)} 1-based 含端。
    优先 gene_biotype=protein_coding；无 gene 则 fallback CDS。"""
    coords = {}
    path = Path(gff_path)
    if not path.exists():
        return coords
    # 2026-09-15: 补 encoding (Linux LANG=C 下 locale 为 ASCII, gff 含非 ASCII 即崩)
    for line in path.read_text(encoding='utf-8', errors='replace').splitlines():
        if line.startswith('#') or '\t' not in line:
            continue
        f = line.rstrip().split('\t')
        if len(f) < 9:
            continue
        seqid, src, ftype, s, e, score, strand, phase, attrs = f
        if ftype not in ('gene', 'CDS'):
            continue
        m = re.search(r'gene[=;](\S+?)(;|$)', attrs) or re.search(r'Name=(\S+?)(;|$)', attrs)
        if not m:
            continue
        name = m.group(1)
        s, e = int(s), int(e)
        # 取该基因最小-最大 extent（含多 CDS）
        if name in coords:
            coords[name] = (min(coords[name][0], s), max(coords[name][1], e))
        else:
            coords[name] = (s, e)
    return coords


def read_gff_sequence_region(gff_path):
    """从 gff3 的 `##sequence-region <seqid> <start> <end>` 指令读参考序列名与全长。
    返回 (seqid, length)；指令缺失或不可解析时返回 (None, None)。
    2026-09-15 新增: 这是判断"哪条比对序列是坐标基准"最可靠的信号 ——
    gene 坐标的上界只是最靠后的注释终点, 不等于基因组全长, 单靠它无法识别参考。"""
    try:
        with open(gff_path, encoding='utf-8', errors='replace') as fh:
            for line in fh:
                if line.startswith('##sequence-region'):
                    f = line.split()
                    if len(f) >= 4:
                        try:
                            return f[1], int(f[3])
                        except ValueError:
                            continue
                elif line.strip() and not line.startswith('#'):
                    break          # 指令区结束
    except OSError:
        pass
    return None, None


_GAP_CHARS = '-?'
# 视为"无信息"的字符 (全 gap/N/X 的样本在基因区间内丢弃, BEAST 无法读入空序列)
_EMPTY_CHARS = set('-?NnXx')


def _read_fasta_pairs(path):
    """读 fasta -> [(name, seq)]，保持文件顺序。
    2026-09-15: 不再 dict 化 —— 原实现同名序列会互相覆盖 (静默丢样本);
    现显式检测重名并报错 (重名会让下游 BEAST tip 冲突, 早报错好过晚出垃圾)。"""
    pairs, name, buf = [], None, []
    with open(path, encoding='utf-8', errors='replace') as fh:
        for line in fh:
            line = line.rstrip('\r\n')
            if not line.strip():
                continue
            if line.startswith('>'):
                if name is not None:
                    pairs.append((name, ''.join(buf)))
                name = line[1:].split()[0]
                buf = []
            else:
                buf.append(line.strip())
    if name is not None:
        pairs.append((name, ''.join(buf)))

    seen = {}
    for nm, _sq in pairs:
        seen[nm] = seen.get(nm, 0) + 1
    dups = sorted(k for k, v in seen.items() if v > 1)
    if dups:
        raise ValueError(f"比对中序列名重复 ({len(dups)} 个): {dups[:8]}; "
                         f"重名会让下游 BEAST tip 冲突, 请先去重")
    return pairs


def _match_seq(pairs, key):
    """按 精确名 > 前缀 匹配序列；返回 (name, seq) 或 None。"""
    for nm, sq in pairs:
        if nm == key:
            return nm, sq
    for nm, sq in pairs:
        if nm.startswith(key):
            return nm, sq
    return None


def _pick_reference(pairs, genome_len, ref_id=None, ref_hint=None):
    """挑选坐标基准参考序列 -> (name, seq)。

    优先级:
      1. ref_id   —— 用户显式指定 (硬要求, 找不到即报错)
      2. ref_hint —— gff ##sequence-region 的 seqid (软提示, 找不到则忽略)
      3. 去 gap 后长度恰好 == genome_len 且唯一
      4. 长度 >= genome_len 的候选唯一
    无法唯一确定时抛 ValueError —— 宁可停下, 也不静默错切。
    """
    if ref_id:
        hit = _match_seq(pairs, ref_id)
        if hit is None:
            raise ValueError(
                f"--ref-id '{ref_id}' 未在比对中找到; 现有序列: "
                f"{[nm for nm, _ in pairs][:10]}{'...' if len(pairs) > 10 else ''}")
        return hit

    if not pairs:
        raise ValueError("比对为空, 无法确定坐标基准")

    if ref_hint:
        hit = _match_seq(pairs, ref_hint)
        if hit is not None:
            return hit

    ulen = [(nm, sq, sum(1 for c in sq if c not in _GAP_CHARS)) for nm, sq in pairs]
    exact = [(nm, sq) for nm, sq, L in ulen if L == genome_len]
    if len(exact) == 1:
        return exact[0]
    cover = [(nm, sq) for nm, sq, L in ulen if L >= genome_len]
    if len(cover) == 1:
        return cover[0]
    cand = exact or cover
    if not cand:
        raise ValueError(
            f"比对中没有任何序列去 gap 后长度 >= 参考基因组长度 {genome_len} "
            f"(最长 {max(L for _, _, L in ulen)}); 坐标基准不成立, 拒绝切分。"
            f"请确认 --gff 与 --fasta 来自同一参考基因组。")
    names = [c[0] for c in cand]
    raise ValueError(
        f"无法唯一确定参考序列 (候选 {len(names)} 条: {names[:8]}"
        f"{'...' if len(names) > 8 else ''}); 请用 --ref-id 显式指定坐标基准。")


def split_genes(fasta_in, coords, out_dir, ref_len=None, ref_id=None, ref_hint=None):
    """按参考基因组坐标切分 gapped 比对 -> out_dir/<gene>.fasta

    2026-09-15 重写 (P0-6): 旧实现把 GFF 的**参考基因组坐标**直接索引到
    「每条样本各自去 gap 后的序列」上 ——
      · 任一样本在基因上游有插入/缺失, 窗口即整体错位 (仅 `len(ug) < e`
        一道长度检查拦不住, 会静默产出错位的基因比对);
      · 逐样本去 gap 后各序列长度不等, 产出的根本不是比对 (BEAST 读不了);
      · 同名序列在 dict 中互相覆盖, 静默丢样本。
    现按 utils/virphy_bridge.slice_genes_from_alignment 的做法: 先用参考序列建立
    「比对列 -> 基因组坐标」映射, 再反向定位每个基因覆盖的比对列区间, 最后**按列**
    切分 (保留 gap) 得到合法子比对。

    ref_id   : 显式指定坐标基准序列名 (硬要求; 推荐)
    ref_hint : 坐标基准提示 (如 gff ##sequence-region 的 seqid; 软匹配, 找不到则忽略)
    ref_len  : 参考基因组全长 (如 gff ##sequence-region 的 end)。旧签名里是
               "健康标记"用途但从未被使用; 现复用为参考长度, 缺省时退化为
               gene 坐标上界 (偏保守, 可能导致自动推断失败并要求 --ref-id)
    返回 {gene: n_seq}; 参考无法确定 / 坐标越界 / 输入非比对时抛 ValueError。
    """
    pairs = _read_fasta_pairs(fasta_in)
    if not coords:
        raise ValueError("gene 坐标为空, 拒绝切分")

    # 输入必须是等长比对 —— 否则"按列切分"没有意义
    lens = {len(sq) for _nm, sq in pairs}
    if len(lens) > 1:
        raise ValueError(
            f"输入不是比对 (序列长度不一致: {sorted(lens)[:6]}"
            f"{'...' if len(lens) > 6 else ''}); 坐标切分模式需要等长 gapped 比对")

    max_coord = max(int(e) for _s, e in coords.values())
    try:
        genome_len = int(ref_len) if ref_len else max_coord
    except (TypeError, ValueError):
        genome_len = max_coord
    if genome_len < max_coord:
        genome_len = max_coord
    ref_name, ref_seq = _pick_reference(pairs, genome_len, ref_id, ref_hint)

    # 比对列 -> 基因组坐标 (0-based); 参考为 gap 的列不占基因组坐标
    aln_to_genome = {}
    gpos = 0
    for apos, base in enumerate(ref_seq):
        if base not in _GAP_CHARS:
            aln_to_genome[apos] = gpos
            gpos += 1
    if gpos < genome_len:
        raise ValueError(
            f"参考序列 {ref_name} 去 gap 后仅 {gpos} bp < 参考基因组长度 {genome_len}; "
            f"--gff 与 --fasta 可能不是同一参考, 拒绝切分。")

    out_dir = Path(out_dir); out_dir.mkdir(parents=True, exist_ok=True)
    counts = {}
    for gene, (s, e) in sorted(coords.items()):
        s, e = int(s), int(e)
        if s < 1 or e < s:
            counts[gene] = 0
            continue
        # 基因区间 [s, e] (1-based 含端) 覆盖的比对列区间
        a_start = a_end = None
        for apos, gv in aln_to_genome.items():
            if s - 1 <= gv <= e - 1:
                if a_start is None:
                    a_start = apos
                a_end = apos
        if a_start is None:
            counts[gene] = 0
            continue
        n = 0
        with open(out_dir / f"{gene}.fasta", 'w', encoding='utf-8') as fh:
            for nm, sq in pairs:
                seg = sq[a_start:a_end + 1]
                if not seg or all(c in _EMPTY_CHARS for c in seg):
                    continue          # 该区间全 gap/N 的样本丢弃
                fh.write(f">{nm}\n")
                for i in range(0, len(seg), 100):
                    fh.write(seg[i:i + 100] + "\n")
                n += 1
        counts[gene] = n
    return counts


# ─────────────────────────── XML 生成 + 校验 ───────────────────────────
def alpha_stanza(lower):
    """返回 siteModel 的 alpha 行（统一的 exact 串替换/核对目标）。"""
    return f'<parameter id="alpha" value="0.5" lower="{lower}"/>'


def generate_and_validate(gene, fa, metadata, out_dir, chain, log_every,
                          alpha_lower, log, clock_model="ucln", tree_prior="skyline"):
    """生成单基因 xml，先做 alpha 下界替换，再做 xml.dom 语法校验。
    返回 (xml_path, ok, msg)。

    clock_model: "strict" | "ucln"  (收敛困难时改 strict 砍掉 ucld.stdev)
    tree_prior:  "constant" | "exponential" | "skyline" | "bdsky"
                 (skyline 每多一组就多一个难识别的种群参数)"""
    res = generate_beast1_phylogeo_xml(
        fasta_file=str(fa), metadata_csv=metadata, output_dir=str(out_dir),
        chain_length=chain, log_every=log_every,
        clock_model=clock_model, tree_prior=tree_prior,
        discretize_locations=False,   # 分基因定年: 纯分子钟 (无 location CTMC/BSSVS)
        substitution_model="auto", gamma_categories=4,
        xml_name=f"{gene}.xml", log=log)
    # 2026-09-15 (审查 P2-21) 修复: 旧版 `Path(res.get("xml_path", ""))` ——
    # Path("") == Path("."), 而 "." 恒存在, 于是 xml_path 缺失时守卫形同虚设,
    # 后续会把仓库目录当成 XML 去 read_text()/write_text()。
    _xp = (res.get("xml_path") or "").strip() if res.get("success") else ""
    xmlp = Path(_xp) if _xp else None
    if not (xmlp and xmlp.exists()):
        return None, False, res.get("error") or "generate_beast1_phylogeo_xml 未返回有效 xml_path"

    txt = xmlp.read_text()
    # 1) alpha 下界：精确整串替换（源头 bridge 若已改 lower 则 count=0 跳过）
    old = alpha_stanza("0.0")
    new = alpha_stanza(alpha_lower)
    if txt.count(old) >= 1:
        txt = txt.replace(old, new)
        xmlp.write_text(txt)
    elif alpha_stanza(alpha_lower) not in txt:
        return xmlp, False, "alpha lower 未生效"

    # 2) XML 语法校验
    try:
        import xml.dom.minidom
        xml.dom.minidom.parse(str(xmlp))
    except Exception as e:
        return xmlp, False, f"XML 语法错误: {e}"

    # 3) 无污染检查（历史坑：反斜杠、双自闭合）
    bad = ('"/>/>' in txt) or ('\\' in txt)
    if bad:
        return xmlp, False, "XML 存在污染（双闭合/反斜杠）"
    return xmlp, True, f"ok (alpha_lower={alpha_lower})"


# ─────────────────────── 后台提交 / watchdog / 解析 ───────────────────────
WATCHDOG_TMPL = r"""#!/bin/bash
# watchdog: 逐基因按各自 chain 阈值判完成，完成即增量解析
set -u
ROOT="{out_dir}"
targ(){{ case "$1" in {big_genes}) echo {chain_big};; *) echo {chain};; esac; }}
log(){{ echo "$(date '+%F %T') $*" >> "$ROOT/WATCH.log"; }}
parse(){{
  {py} "$ROOT/parse_gene_beast.py" >> "$ROOT/WATCH.log" 2>&1
  log "parser invoked rc=$?"
}}
log "watch started for {genes}"
prev_done=0
for i in $(seq 1 720); do
    done=0
    nprocs=0
    for g in {genes}; do
        last=$(awk -F'\t' 'NR>1{{v=$1}}END{{print v}}' "$ROOT/$g/$g.log" 2>/dev/null)
        t=$(targ "$g")
        [ -n "$last" ] && [ "$last" -ge "$t" ] && done=$((done+1))
        pgrep -f "$ROOT/$g/$g.xml" >/dev/null 2>&1 && nprocs=$((nprocs+1))
    done
    log "check #$i: done=$done/{n_genes} procs=$nprocs"
    if [ "$done" -gt "$prev_done" ]; then parse; prev_done=$done; fi
    if [ "$done" -ge "{n_genes}" ] && [ "$nprocs" -eq 0 ]; then parse; log "ALL_DONE"; break; fi
    sleep 600
done
"""

PARSE_PY = r'''#!/usr/bin/env python3
"""parse_gene_beast.py — 离线解析分基因 BEAST 的 MRCA/HPD/ESS/拓扑 → markdown。"""
import re
from pathlib import Path
import numpy as np

ROOT = Path(__OUTDIR__)
GENES = __GENES__
FULLGENOME_TMRCA = __FULLTMRCA__   # 可比对照（全长定年），None 则跳过
ESS_THRESHOLD = 200.0              # 收敛阈值 (Tracer 惯例)

def hpd(sample, mass=0.95):
    s = np.sort(sample); n = len(s)
    n_in = max(1, int(np.ceil(mass*n)))  # ceil 口径, 对齐 utils.ess.calculate_95hpd/Tracer
    best = (s[-1]-s[0], 0)
    for i in range(n - n_in + 1):
        w = s[i+n_in-1]-s[i]
        if w < best[0]: best = (w, i)
    return s[best[1]], s[best[1]+n_in-1]

def calc_ess(x):
    """ESS (Geyer 1992 initial POSITIVE sequence, Tracer 同算法)
    与 utils/ess.py v3 逐行一致 (2026-08-27 trace.jar 官方库逐字节码对账通过):
    ddof=0 有偏方差, maxLag=min(n-1,2000), 成对和<=0 即停, 无单调性检查。
    独立脚本内联移植 (部署后不能 import utils.ess)。"""
    n = len(x)
    if n < 4:
        return float(n)
    x = np.asarray(x, dtype=float)
    centered = x - np.mean(x)
    var = np.var(centered)  # 有偏 (ddof=0), 同 Tracer gammaStat[0]
    if var == 0:
        return float(n)
    max_lag = min(n - 1, 2000)  # Tracer MAX_LAG = 2000
    acf = np.zeros(max_lag + 1)
    acf[0] = 1.0
    for lag in range(1, max_lag + 1):
        if lag >= n:
            break
        numerator = np.sum(centered[:n - lag] * centered[lag:])
        denominator = var * (n - lag)
        acf[lag] = numerator / denominator if denominator != 0 else 0.0
    # Geyer initial POSITIVE: varStat = gamma0 + 2*sum(成对和), <=0 即停, 无单调检查
    sum_pairs = 0.0
    k = 1
    while 2 * k <= max_lag:
        pair = acf[2 * k - 1] + acf[2 * k]
        if pair <= 0:
            break
        sum_pairs += pair
        k += 1
    tau_hat = 1.0 + 2.0 * sum_pairs
    if tau_hat <= 0:
        return float(n)
    ess = n / tau_hat
    return max(1.0, min(ess, float(n)))

def rootheights(logf, burn=0.2):
    """返回 (after_burn_vals, ess, n_total); 无效返回 None"""
    if not Path(logf).exists(): return None
    lines = Path(logf).read_text().splitlines()
    if len(lines) < 2: return None
    # skip leading '#' comment lines; the column header is the first non-# line
    idx = 0
    while idx < len(lines) and lines[idx].startswith('#'):
        idx += 1
    if idx >= len(lines): return None
    header = [c.strip() for c in lines[idx].split('\t')]
    col = None
    for cand in ("treeModel.rootHeight", "RootHeight"):
        if cand in header: col = header.index(cand); break
    if col is None:
        col = next((i for i,h in enumerate(header) if "rootheight" in h.lower()), None)
    if col is None: return None
    vals = []
    for ln in lines[idx+1:]:
        if not ln.strip() or ln.startswith('#'): continue
        p = ln.split('\t')
        try: vals.append(float(p[col]))
        except (ValueError, IndexError): pass
    if len(vals) < 10: return None  # 最少 10 采样; 采样数低时 ESS 照算并标未收敛
    n_total = len(vals)
    after = vals[int(len(vals)*burn):]
    ess = calc_ess(np.array(after))
    return after, ess, n_total

def main():
    rows, notes = [], []
    for g in GENES:
        d = ROOT/g; rh = rootheights(d/f"{g}.log")
        if not rh:
            rows.append((g, None, None, None, None, None))
            notes.append(f"{g}: 未完成/无 rootHeight"); continue
        vals, ess, n_total = rh
        a = np.array(vals); med=float(np.median(a)); lo,hi=hpd(a)
        conv = "✓" if ess >= ESS_THRESHOLD else "⚠"
        rows.append((g, med, hi, lo, ess, conv))
        if ess < ESS_THRESHOLD:
            notes.append(f"{g}: ESS={ess:.0f} < {ESS_THRESHOLD:.0f} — 未收敛, 建议加长 chain")
    md = ["# 分基因分区 BEAST 定年结果", ""]
    if FULLGENOME_TMRCA is not None:
        md += [f"全长基线 TMRCA: **{FULLGENOME_TMRCA:.1f} yr**", ""]
    md += ["| 基因 | median TMRCA (yr) | 95% HPD (yr) | vs 全长 | ESS | 收敛 |",
           "|---|---|---|---|---|---|"]
    for g,m,hi,lo,ess,conv in rows:
        if m is None:
            md += [f"| {g} | — | — | — | — | — |"]; continue
        vs = "—" if FULLGENOME_TMRCA is None else f"{m-FULLGENOME_TMRCA:+.1f}"
        es = f"{ess:.0f}" if ess is not None else "—"
        md += [f"| {g} | {m:.1f} | {lo:.1f}–{hi:.1f} | {vs} | {es} | {conv} |"]
    md += ["", "## Note", ""] + [f"- {n}" for n in notes]
    md += ["", f"收敛判定: ESS ≥ {ESS_THRESHOLD:.0f} (Tracer 惯例); ⚠ = 未收敛, TMRCA 可能不可靠", ""]
    (ROOT/"GENE_PARTITION_RESULTS.md").write_text("\n".join(md))
    print("\n".join(md))

if __name__ == "__main__": main()
'''

# ────────────────────────────────── 主流程 ──────────────────────────────────
def main():
    # 2026-09-15 修复: 本脚本从不配置 logging, 而 LogCollector.emit() 走
    # logging.getLogger(__name__).info() —— 无 handler 时 INFO 被静默丢弃,
    # 独立运行时全程看不到任何进度 (只有 sys.exit 的 stderr 能看见)。
    import logging as _logging
    _logging.basicConfig(level=_logging.INFO,
                         format='%(asctime)s [%(levelname)s] %(message)s',
                         datefmt='%H:%M:%S')
    ap = argparse.ArgumentParser(description='分基因分区 BEAST 定年一键入口')
    ap.add_argument('--virus', default='GCVA')
    ap.add_argument('--fasta', default=None, help='全长比对 (gapped MAFFT)，坐标切分模式必需')
    ap.add_argument('--metadata', required=True, help='name,date,location CSV')
    ap.add_argument('--gff', default=None, help='参考基因组 gff3 (gene 坐标)，坐标切分模式必需')
    ap.add_argument('--ref-id', default=None,
                    help='坐标切分模式的参考序列名 (gff 坐标基准)。缺省时自动推断, '
                         '推断不唯一会报错并要求显式指定')
    ap.add_argument('--out', default='gene_partition', help='工作目录')
    ap.add_argument('--chain', type=int, default=DEFAULT_CHAIN, help='默认每基因 chain')
    ap.add_argument('--chain-big', type=int, default=10_000_000, help='大基因 chain')
    ap.add_argument('--big-genes', default='L', help='大基因逗号分隔 (默认 L)')
    ap.add_argument('--alpha-lower', default=DEFAULT_ALPHA_LOWER, help='gamma alpha 下界 (默认 0.5)')
    ap.add_argument('--threads', type=int, default=24, help='每基因线程')
    ap.add_argument('--log-every', type=int, default=2000)
    ap.add_argument('--full-tmrca', type=float, default=None, help='全长定年 TMRCA 供对照')
    ap.add_argument('--notify', default=None, help='完成时命令 (可选, 不执行)')
    ap.add_argument('--genes-dir', default=None,
                    help='现成基因比对目录（capheine/cawlign 产物，按 <gene>.fasta 命名）；给出则跳过坐标切分（推荐）')
    ap.add_argument('--submit-only', action='store_true', help='跳过生成只提交')
    ap.add_argument('--clock-model', default='ucln', choices=['strict', 'ucln'],
                    help='分子钟模型 (默认 ucln)。ESS 卡在 ucld.stdev 时可换 strict')
    ap.add_argument('--tree-prior', default='skyline',
                    choices=['constant', 'exponential', 'skyline', 'bdsky'],
                    help='溯祖先验 (默认 skyline)。ESS 卡在 skyline.popSize 时可换 constant')
    args = ap.parse_args()

    out = Path(args.out); out.mkdir(parents=True, exist_ok=True)
    log = LogCollector()
    genes_dir = out / "gene_split"

    big_genes = {g.strip() for g in args.big_genes.split(',') if g.strip()}
    if args.genes_dir:
        # 已有参考锚定基因比对（capheine/cawlign 产物）：跳过坐标切分
        gd_src = Path(args.genes_dir)
        files = sorted(gd_src.glob('*.fasta')) or sorted(gd_src.glob('*.fas'))
        if not files:
            sys.exit(f"--genes-dir 无 *.fasta/*.fas (fd: {args.genes_dir})")
        # capheine/cawlign 命名适配: ref_cds-noStopCodons.part_G-aligned.fasta -> G.fasta
        # (cawlign 产物文件名带前缀与后缀, 逐基因 BEAST 需要干净的 <gene>.fasta)
        for f in files:
            m = re.search(r'\.part_([A-Za-z0-9]+)(?:-aligned|-nodups)?\.(?:fas|fasta)$', f.name)
            if m:
                gname = m.group(1)
                if f.name != f'{gname}.fasta':
                    lnk = gd_src / f'{gname}.fasta'
                    if not lnk.exists():
                        try:
                            lnk.symlink_to(f.name)
                        except OSError:
                            shutil.copy2(f, lnk)
        genes_dir = gd_src
        # ordered 只取干净命名 (排除 capheine 原始命名: .part_/-aligned/-nodups 前缀文件)
        ordered = sorted(
            p.stem for p in gd_src.glob('*.fasta')
            if not re.search(r'(?:\.part_|-aligned|-nodups)', p.stem)
        )
        if not ordered:
            sys.exit(f"--genes-dir 无 *.fasta/*.fas (fd: {args.genes_dir})")
        # 完整性审计 (2026-08-27, GCVA L/G 教训): capheine/cawlign 可能静默丢基因,
        # 对照 seqkit part_* 全集, 缺失则 WARNING (不阻塞, 便于用户决定补跑)
        part_full = {m.group(1) for p in gd_src.glob('*.fasta')
                     if (m := re.search(r'\.part_([A-Za-z0-9]+)', p.stem))}
        missing = part_full - set(ordered)
        if missing:
            log.emit(f"WARNING: 上游比对缺失基因 {sorted(missing)} "
                     f"(seqkit 切分有但 cawlign 未产出比对); "
                     f"分基因定年只覆盖 {sorted(ordered)}, 如需全部基因请检查 capheine 过滤或用 gene_align_generator fallback 补跑")
        log.emit(f"genes-dir 模式: {len(ordered)} 个基因比对 -> {ordered}")
    else:
        if not (args.fasta and args.gff):
            sys.exit("坐标切分模式需 --fasta 与 --gff；或改用 --genes-dir 直接给基因比对")
        coords = read_gene_coords_gff(args.gff)
        if not coords:
            sys.exit("未从 gff3 解析到 gene 坐标；或用 --genes-dir 直接给基因比对")
        ordered = sorted(coords.keys())
        if not args.submit_only:
            # 参考基因组长度/序列名: --ref-id 优先, 其次 gff ##sequence-region
            # 2026-09-15: gene 坐标上界 != 基因组全长, 单靠它无法可靠识别参考序列
            _seqid, _region_len = read_gff_sequence_region(args.gff)
            try:
                counts = split_genes(args.fasta, coords, genes_dir,
                                     ref_len=_region_len,
                                     ref_id=args.ref_id,
                                     ref_hint=_seqid)
            except ValueError as e:
                # 2026-09-15: 坐标基准无法确定时拒绝静默错切, 明确报错并给出解法
                sys.exit(f"基因坐标切分失败: {e}")
            log.emit(f"切分基因坐标: {counts}")

    # 逐基因生成 + 校验（两种模式统一；submit-only 跳过生成只提交）
    if not args.submit_only:
        for gene in ordered:
            fa = genes_dir / f"{gene}.fasta"
            if not fa.exists():
                log.emit(f"{gene}: 无比对文件，跳过"); continue
            gd = out / gene; gd.mkdir(parents=True, exist_ok=True)
            chain = args.chain_big if gene in big_genes else args.chain
            xmlp, ok, msg = generate_and_validate(
                gene, str(fa), args.metadata, str(gd), chain, args.log_every,
                args.alpha_lower, log,
                clock_model=args.clock_model, tree_prior=args.tree_prior)
            log.emit(f"{gene}: {msg} (chain={chain})")

    # 落 watchdog + 解析脚本
    genes_s = ' '.join(ordered) if ordered else 'N P P4 M G L'
    j = ', '.join(f'"{g}"' for g in ordered) if ordered else '[]'
    full = 'None' if args.full_tmrca is None else str(args.full_tmrca)
    # big_genes 逗号 → case pattern 的 | (历史坑: 'L,M' 直接注入 case 永不匹配)
    big_pattern = str(args.big_genes).replace(',', '|')
    # watchdog 文件名中性化 (2026-08-27): 原硬编码 gcva_ 前缀对其他病毒误导
    watch_name = f"{out.name}_watch_parse.sh"
    (out / watch_name).write_text(WATCHDOG_TMPL.format(
        out_dir=out, big_genes=big_pattern, chain=args.chain,
        chain_big=args.chain_big, genes=genes_s, n_genes=len(ordered),
        py=sys.executable))
    # __OUTDIR__ 用 repr 转义 (历史坑: 反斜杠路径里的 \v/\t/\n 被当转义)
    (out / "parse_gene_beast.py").write_text(PARSE_PY
        .replace("__OUTDIR__", repr(str(out)))
        .replace("__GENES__", j if j else '[]')
        .replace("__FULLTMRCA__", full))

    # 后台启动各基因 BEAST（setsid detach；det 到各自目录）
    # BEAST 二进制: 优先 which, 其次 datasets.yaml env, 再次 PHYLO_BEAST_BIN 环境变量。
    # 2026-09-15 (审查 P2-22): 删除硬编码 '/home/zhangwenda/mambaforge/bin/beast'
    # 兜底 —— 换机器后它必然不存在, 却会把"找不到 BEAST"伪装成一个具体的路径问题。
    # 现在找不到就明确报"未找到", 并提示两种指定方式。
    from utils.dataset_config import env as _env
    beast = (shutil.which('beast') or _env('beast_bin')
             or os.environ.get('PHYLO_BEAST_BIN') or '')
    if not beast or not os.path.exists(beast):
        log.emit(f"BEAST 二进制不可用 (解析结果: {beast or '<未找到>'}), 分基因定年无法提交")
        log.emit("  提示: 用 PHYLO_BEAST_BIN=/path/to/beast 或 datasets.yaml env.beast_bin 指定")
        return
    # 确定性种子基准 (2026-09-15 修复): 旧版 os.system 起 BEAST 不传 -seed
    # → 每次 chain 结果不同, 分基因 tMRCA 不可复现。
    try:
        _seed_base = int(os.environ.get("PHYLO_SEED_BASE", "") or 20260915)
    except ValueError:
        _seed_base = 20260915
    for _gi, gene in enumerate(ordered):
        gd = out / gene; xmlp = gd / f'{gene}.xml'
        if not xmlp.exists():
            log.emit(f"{gene}: 无 xml，跳过提交"); continue
        thr = max(8, args.threads)
        seed = _seed_base + (_gi + 1) * 1024
        # 用绝对路径 XML: watchdog 的 pgrep -f "$ROOT/$g/$g.xml" 才能匹配到进程
        # 路径加引号 (2026-09-15 修复): 含空格即被拆成多个参数
        cmd = (f'cd {shlex.quote(str(gd))} && setsid nohup {shlex.quote(beast)} '
               f'-overwrite -threads {int(thr)} -seed {int(seed)} '
               f'{shlex.quote(str(gd / f"{gene}.xml"))} '
               f'>beast_stdout.txt 2>beast_stderr.txt </dev/null &')
        # 2026-09-15: os.system 返回的是 wait status (rc<<8) 而非真实退出码,
        # 且没有超时保护。改 subprocess.run —— 命令自带 `&` 后台化且 stdout/
        # stderr 已重定向到文件, 不会阻塞; 超时仅作为兜底。
        try:
            _rc = subprocess.run(cmd, shell=True, timeout=120).returncode
        except subprocess.TimeoutExpired:
            _rc = -1
        if _rc != 0:
            log.emit(f"{gene}: ⚠️ 提交命令返回码 {_rc}, 请检查 beast_stderr.txt")
        else:
            log.emit(f"{gene}: 已后台提交 (threads={thr}, seed={seed})")
        time.sleep(2)

    log.emit("部署完成。生成文件:")
    for f in (watch_name, 'parse_gene_beast.py'):
        if (out / f).exists():
            log.emit(f"  {out / f}")
    print("DONE")


if __name__ == '__main__':
    main()
