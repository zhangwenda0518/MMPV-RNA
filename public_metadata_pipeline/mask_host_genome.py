#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
mask_host_genome.py — 宿主参考基因组的「病毒样区域」掩码工具 (EVE / EPRV)
============================================================================
把宿主参考基因组中指定的区间掩码成 `N`（或小写），用于**建库前**处理：

    public_metadata_pipeline:  all.genome.uniq.fasta  →  all.genome.masked.fasta
                               ↓ build_hostbase.py (kraken2/bowtie2/hisat2/minimap2)

动机
----
现有去宿主是「宿主基因组全序列」当靶子：`host_depletion.py` 的比对阶段用
`SAM.view_filter(mfilter=True)`（`-f 12 -F 256`，PE）只按「是否比对上」判宿主，
比对到宿主 EVE（内源病毒元件）的外源病毒 reads 会被整段删除，且这些 reads 在
Kraken2 阶段通常并不会被分类为病毒，因此 `--keep-taxids 10239` 的救援**覆盖不到**
（救援发生在比对之前，见 host_depletion.py:651-991 与 SAM.view_filter:491-502）。

把 EVE 区间掩码掉之后，这些位置不再产生 k-mer / 比对种子 → 该区间的 reads
不再被吸收，可进入下游组装与鉴定。

输入
----
1) `--bed`：BED 文件（0 基半开，至少 3 列 `contig start end`），可多个取并集。
2) `--verdict-tsv`：EVE 筛查的判定表，可多个取并集：
   - 单基因组 `06_EVE/02_Verdict/<NAME>/<NAME>.s2_verdict.tsv`
     （列 `locus verdict best_viral_id ...`，`locus` = `contig:start-end`，1 基闭区间）
   - 跨基因组 `kingdom_loci_all.tsv.gz`
     （列 `genome locus verdict ...`）
   由 `--verdict-values`（默认 `viral_supported`）过滤；`undetermined` 默认不进掩码。

输出
----
- 掩码后的 FASTA：**长度不变**、坐标不变（N 化而不是删除，否则下游 GFF/坐标全错）
- 逐 contig 报告 TSV：`contig length n_loci masked_bp masked_frac`
- stderr 摘要（掩码 bp、占比、未命中的 BED contig 数）

设计取舍
--------
- 只用标准库，不 import 仓库内 `utils/mmpv_common`：该脚本要在服务器任意目录直接跑。
- 逐 contig 流式处理（单条 contig 在内存里 mask），不整基因组载入；枸杞基因组
  约 2.2 Gb，单条染色体级别 contig 载入无压力。

示例
----
    # 1) 单基因组：EVE 判定表直接当掩码源（viral_supported 才掩）
    python mask_host_genome.py \\
        -g host_reference/genome/all.genome.uniq.fasta \\
        -v 06_EVE/02_Verdict/ningxia/ningxia.s2_verdict.tsv \\
        --out host_reference/genome/all.genome.masked.fasta \\
        --report host_reference/genome/mask.report.tsv

    # 2) 多个基因组批量 + 外扩 50 bp 吞掉边界错配
    python mask_host_genome.py \\
        -g ningxia.genome.fasta heiguo.genome.fasta zhonghua.genome.fasta \\
        -v 06_EVE/kingdom_loci_all.tsv.gz --pad 50 --out-dir mask-genome/

    # 3) 与 reference_masker.py 的伪读长 BED 叠加（核酸层 + 蛋白层取并集）
    python mask_host_genome.py -g all.genome.uniq.fasta \\
        -b mask-genome/ningxia.masked.bed -v 06_EVE/kingdom_loci_all.tsv.gz \\
        --out all.genome.masked.fasta
"""

from __future__ import annotations

import argparse
import gzip
import os
import re
import sys
from typing import Dict, Iterable, Iterator, List, Sequence, Tuple

if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
        sys.stderr.reconfigure(encoding="utf-8", errors="replace")
    except Exception:  # pragma: no cover - 老环境无 reconfigure
        pass

# `contig:start-end`（1 基闭区间；contig 名里允许出现 ':'，从左贪婪匹配到最后一个
# ':' 前的数字段，EVE 管线自己产出的位点名就是这种形态）
_LOCUS_RE = re.compile(r"^(?P<contig>.+):(?P<start>\d+)-(?P<end>\d+)$")


# ==========================================
# 通用 IO
# ==========================================
def open_maybe_gz(path: str):
    """按魔数判断 gzip，返回文本句柄（FASTA/BED/TSV 都是文本）。"""
    with open(path, "rb") as fh:
        magic = fh.read(2)
    if magic == b"\x1f\x8b":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "r", encoding="utf-8", errors="replace")


def iter_fasta(path: str) -> Iterator[Tuple[str, bytes]]:
    """流式读 FASTA：逐条 yield (完整头行, 序列 bytes)。不假设行宽。"""
    header = None
    chunks: List[bytes] = []
    with open_maybe_gz(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if header is not None:
                    yield header, b"".join(chunks)
                header = line[1:].strip()
                chunks = []
            elif header is not None:
                s = line.strip()
                if s:
                    chunks.append(s.encode("ascii", "replace"))
        if header is not None:
            yield header, b"".join(chunks)


def write_fasta(fh, header: str, seq: bytes, line_width: int = 60) -> None:
    fh.write(">" + header + "\n")
    if line_width <= 0:
        fh.write(seq.decode("ascii") + "\n")
        return
    view = memoryview(seq)
    for i in range(0, len(seq), line_width):
        fh.write(view[i:i + line_width].tobytes().decode("ascii") + "\n")


# ==========================================
# 区间来源
# ==========================================
def parse_bed(paths: Sequence[str]) -> Tuple[Dict[str, List[Tuple[int, int]]], int, int]:
    """解析 BED（0 基半开）。返回 (contig -> 区间列表, 有效行数, 跳过行数)。"""
    out: Dict[str, List[Tuple[int, int]]] = {}
    n_ok = 0
    n_skip = 0
    for p in paths:
        with open_maybe_gz(p) as fh:
            for line in fh:
                s = line.strip()
                if not s or s.startswith(("#", "track", "browser")):
                    continue
                c = s.split("\t")
                if len(c) < 3:
                    c = s.split()
                if len(c) < 3:
                    n_skip += 1
                    continue
                try:
                    start, end = int(c[1]), int(c[2])
                except ValueError:
                    n_skip += 1
                    continue
                if start < 0 or end <= start:
                    n_skip += 1
                    continue
                out.setdefault(c[0].strip(), []).append((start, end))
                n_ok += 1
    return out, n_ok, n_skip


def parse_verdict_tsv(paths: Sequence[str], values: Iterable[str]
                      ) -> Tuple[Dict[str, List[Tuple[int, int]]], int, int]:
    """解析 EVE 判定表，只保留 verdict ∈ values 的位点。

    支持两种形态:
      · 单基因组 `<NAME>.s2_verdict.tsv`：`locus verdict ...`
      · 跨基因组 `kingdom_loci_all.tsv.gz`：`genome locus verdict ...`
    位点名 `contig:start-end` 为 1 基闭区间，转 BED 时 start-1。
    """
    want = {v.strip() for v in values if v and v.strip()}
    out: Dict[str, List[Tuple[int, int]]] = {}
    n_ok = 0
    n_skip = 0
    for p in paths:
        with open_maybe_gz(p) as fh:
            first = True
            i_locus, i_verdict = 0, 1
            for line in fh:
                s = line.rstrip("\n")
                if not s.strip():
                    continue
                c = s.split("\t")
                if first:
                    first = False
                    low = [x.strip().lower() for x in c]
                    if "locus" in low and "verdict" in low:
                        i_locus = low.index("locus")
                        i_verdict = low.index("verdict")
                        continue
                    # 无表头: 退回默认列位
                if len(c) <= max(i_locus, i_verdict):
                    n_skip += 1
                    continue
                if c[i_verdict].strip() not in want:
                    continue
                m = _LOCUS_RE.match(c[i_locus].strip())
                if not m:
                    n_skip += 1
                    continue
                start = int(m.group("start")) - 1
                end = int(m.group("end"))
                if start < 0 or end <= start:
                    n_skip += 1
                    continue
                out.setdefault(m.group("contig"), []).append((start, end))
                n_ok += 1
    return out, n_ok, n_skip


def merge_intervals(ivs: Iterable[Tuple[int, int]], gap: int = 0
                    ) -> List[Tuple[int, int]]:
    """排序合并；`gap` 为允许合并的最大间隙（0 = 仅重叠/相接才合并）。"""
    merged: List[Tuple[int, int]] = []
    for s, e in sorted(ivs):
        if merged and s <= merged[-1][1] + max(0, gap):
            if e > merged[-1][1]:
                merged[-1] = (merged[-1][0], e)
        else:
            merged.append((s, e))
    return merged


def pad_intervals(ivs: Iterable[Tuple[int, int]], pad: int, length: int
                  ) -> List[Tuple[int, int]]:
    """两侧外扩并裁到 [0, length]，随后再合并一次（外扩后可能重叠）。"""
    if pad <= 0:
        return merge_intervals(ivs, 0)
    padded = [(max(0, s - pad), min(length, e + pad)) for s, e in ivs]
    return merge_intervals(padded, 0)


# ==========================================
# 掩码
# ==========================================
def apply_mask(seq: bytes, ivs: Sequence[Tuple[int, int]], mode: str = "N") -> bytes:
    """把 [start, end) 区间替换成 N（或小写）。返回新 bytes，长度不变。"""
    if not ivs:
        return seq
    ba = bytearray(seq)
    n = len(ba)
    if mode == "lower":
        for s, e in ivs:
            cs, ce = max(0, s), min(n, e)
            if ce > cs:
                ba[cs:ce] = bytes(ba[cs:ce]).lower()
        return bytes(ba)
    fill = b"\x4e"  # ASCII 'N'
    for s, e in ivs:
        cs, ce = max(0, s), min(n, e)
        if ce > cs:
            ba[cs:ce] = fill * (ce - cs)
    return bytes(ba)


def mask_genome(genome: str, ivs_by_contig: Dict[str, List[Tuple[int, int]]],
                out_path: str, report_rows: List[List[str]], genome_name: str,
                pad: int, merge_gap: int, mode: str, line_width: int,
                quiet: bool = False) -> Dict[str, float]:
    """单基因组：写出掩码 FASTA + 逐 contig 报告行。返回汇总指标。"""
    used = set()
    tot_len = tot_masked = tot_loci = 0
    os.makedirs(os.path.dirname(os.path.abspath(out_path)) or ".", exist_ok=True)
    with open(out_path, "w", encoding="utf-8") as fo:
        for header, seq in iter_fasta(genome):
            contig = header.split()[0] if header else ""
            seq_len = len(seq)
            raw = ivs_by_contig.get(contig, [])
            if raw:
                used.add(contig)
                ivs = merge_intervals(raw, merge_gap)
                ivs = pad_intervals(ivs, pad, seq_len)
            else:
                ivs = []
            masked = sum(min(e, seq_len) - max(0, s) for s, e in ivs if e > s)
            if ivs:
                seq = apply_mask(seq, ivs, mode)
            write_fasta(fo, header, seq, line_width)
            frac = (masked / seq_len) if seq_len else 0.0
            report_rows.append([genome_name, contig, str(seq_len),
                                str(len(ivs)), str(masked), f"{frac:.6f}"])
            tot_len += seq_len
            tot_masked += masked
            tot_loci += len(ivs)

    missing = [c for c in ivs_by_contig if c not in used]
    if missing and not quiet:
        sys.stderr.write(
            f"[WARN] {genome_name}: {len(missing)} 个区间 contig 在 FASTA 中不存在"
            f"（示例: {', '.join(missing[:3])}）; EVE 表的 contig 名与 FASTA 头不一致时"
            f" 常见，掩码会被静默跳过\n")
    return {
        "genome": genome_name,
        "n_bp": float(tot_len),
        "n_masked_bp": float(tot_masked),
        "n_loci": float(tot_loci),
        "masked_frac": (tot_masked / tot_len) if tot_len else 0.0,
        "n_contigs": float(len(report_rows)),
    }


# ==========================================
# 主流程
# ==========================================
def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="mask_host_genome.py",
        description="把宿主参考基因组的「病毒样区域」(EVE/EPRV) 掩码成 N，供建库使用",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    g = p.add_argument_group("输入")
    g.add_argument("-g", "--genome", nargs="+", required=True,
                   help="宿主参考基因组 FASTA (.fa/.fasta/.fna[.gz])，可多个")
    g.add_argument("-b", "--bed", nargs="*", default=[],
                   help="掩码区间 BED (0 基半开, contig start end)，可多个取并集")
    g.add_argument("-v", "--verdict-tsv", nargs="*", default=[],
                   help="EVE 判定表 (.s2_verdict.tsv / kingdom_loci_all.tsv.gz)，"
                        "可多个取并集")
    g.add_argument("--verdict-values", default="viral_supported",
                   help="判定表中要掩码的 verdict（逗号分隔；默认 viral_supported；"
                        "设为 viral_supported,undetermined 可做「激进」掩码臂）")

    g = p.add_argument_group("掩码参数")
    g.add_argument("--pad", type=int, default=0,
                   help="两侧外扩 bp（默认 0；建议 25-50 吞掉边界错配与短的相邻同源段）")
    g.add_argument("--merge-gap", type=int, default=0,
                   help="区间合并的最大间隙 bp（默认 0）")
    g.add_argument("--mask-mode", choices=["N", "lower"], default="N",
                   help="N: 替换成 N（默认，建库首选）; lower: 转小写（软掩码）")
    g.add_argument("--line-width", type=int, default=60, help="输出 FASTA 行宽（0=不换行）")

    g = p.add_argument_group("输出")
    g.add_argument("-o", "--out", help="单基因组输出 FASTA 路径")
    g.add_argument("--out-dir", help="多基因组输出目录（<stem>.masked.fasta）")
    g.add_argument("--report", help="单基因组逐 contig 报告 TSV（默认 <out>.report.tsv）")
    g.add_argument("--quiet", action="store_true", help="静默（只输出报告文件）")
    return p


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)

    if not args.bed and not args.verdict_tsv:
        sys.stderr.write("[ERROR] 至少需要一个 --bed 或 --verdict-tsv\n")
        return 2
    if len(args.genome) > 1 and not args.out_dir:
        sys.stderr.write("[ERROR] 多基因组必须给 --out-dir\n")
        return 2
    if len(args.genome) == 1 and not (args.out or args.out_dir):
        sys.stderr.write("[ERROR] 需要 --out 或 --out-dir\n")
        return 2

    bed_ivs, n_bed, n_bed_skip = parse_bed(args.bed) if args.bed else ({}, 0, 0)
    vals = [x for x in args.verdict_values.split(",") if x.strip()]
    vd_ivs, n_vd, n_vd_skip = (parse_verdict_tsv(args.verdict_tsv, vals)
                               if args.verdict_tsv else ({}, 0, 0))

    # 两源合并（同一 contig 的区间随后在 mask_genome 里统一 merge/merge_gap/pad）
    all_ivs: Dict[str, List[Tuple[int, int]]] = {}
    for src in (bed_ivs, vd_ivs):
        for contig, ivs in src.items():
            all_ivs.setdefault(contig, []).extend(ivs)

    if not args.quiet:
        sys.stderr.write(
            f"[INFO] 区间来源: BED {n_bed} 行 (跳过 {n_bed_skip}) + "
            f"verdict({','.join(vals)}) {n_vd} 位点 (跳过 {n_vd_skip}); "
            f"合并后 {sum(len(v) for v in all_ivs.values())} 段 / "
            f"{len(all_ivs)} 个 contig\n")

    report_rows: List[List[str]] = []
    summaries: List[Dict[str, float]] = []
    for genome in args.genome:
        stem = os.path.basename(genome)
        for suffix in (".gz", ".fasta", ".fa", ".fna"):
            if stem.lower().endswith(suffix):
                stem = stem[: -len(suffix)]
        if args.out_dir:
            out_path = os.path.join(args.out_dir, f"{stem}.masked.fasta")
            report_path = os.path.join(args.out_dir, f"{stem}.masked.report.tsv")
        else:
            out_path = args.out
            report_path = args.report or (out_path + ".report.tsv")
        summaries.append(mask_genome(
            genome, all_ivs, out_path, report_rows, stem,
            pad=args.pad, merge_gap=args.merge_gap, mode=args.mask_mode,
            line_width=args.line_width, quiet=args.quiet,
        ))
        with open(report_path, "w", encoding="utf-8") as fo:
            fo.write("genome\tcontig\tlength\tn_loci\tmasked_bp\tmasked_frac\n")
            for row in report_rows:
                if row[0] == stem:
                    fo.write("\t".join(row) + "\n")
        if not args.quiet:
            s = summaries[-1]
            sys.stderr.write(
                f"[DONE] {stem}: {int(s['n_bp']):,} bp / {int(s['n_contigs'])} contig; "
                f"掩码 {int(s['n_masked_bp']):,} bp ({s['masked_frac'] * 100:.3f}%), "
                f"{int(s['n_loci'])} 段 -> {out_path}\n")

    if len(summaries) > 1:
        tot = sum(s["n_bp"] for s in summaries)
        m = sum(s["n_masked_bp"] for s in summaries)
        sys.stderr.write(f"[TOTAL] {len(summaries)} 基因组, 掩码 {int(m):,} bp "
                         f"({(m / tot * 100) if tot else 0:.3f}%)\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
