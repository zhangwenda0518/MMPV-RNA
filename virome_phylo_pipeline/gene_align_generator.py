"""gene_align_generator — capheine 产物的 fallback 生成器

当 gene_dating (分基因 BEAST 定年) 找不到 capheine/cawlign 基因比对产物时,
本模块自动生成基因分区比对 (密码子水平), 供 gene_partition_dating.py 消费。

数据流 (对齐 capheine_pipeline 的 cawlign 核心步骤, 只到基因比对层, 不跑 HyPhy):
    reference CDS (ref_cds.fasta / gbk CDS / 用户提供)
        → 去末端终止密码子 (remove_terminal_stop_codon)
        → 按基因拆分 (split_fasta_pure)
        → 每基因 cawlign -t codon -r <gene> -f refmap <unaligned 全序列>
        → time/gene_split/<gene>.fasta   (干净命名, gene_partition_dating 的 --genes-dir 约定)

幂等: 目标目录已有 <gene>.fasta 则复用, 不重复计算。
"""
from __future__ import annotations

import argparse
import logging
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path
from typing import List, Optional

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

logger = logging.getLogger("gene_align_generator")

# ---------------------------------------------------------------- helpers ---

def _find_tool(cmd: str) -> Optional[str]:
    """在 PATH 找可执行工具 (cawlign 等), 找不到返回 None"""
    return shutil.which(cmd)


def remove_terminal_stop_codon(fasta_path: Path, output_path: Path,
                               genetic_code: str = "1") -> Path:
    """去除每个 CDS 末端的终止密码子 (与 capheine_pipeline 同逻辑)"""
    from Bio.Data import CodonTable
    table = CodonTable.unambiguous_dna_by_id[int(genetic_code)]
    stop_codons = set(table.stop_codons)
    records_out = []
    for record in SeqIO.parse(fasta_path, "fasta"):
        seq_str = str(record.seq).upper().replace("U", "T")
        idx = len(seq_str)
        trailing_stops = 0
        while idx >= 3:
            if seq_str[idx - 3: idx] in stop_codons:
                trailing_stops += 1
                idx -= 3
            else:
                break
        for pos in range(0, (idx // 3) * 3, 3):
            if seq_str[pos: pos + 3] in stop_codons:
                raise RuntimeError(f"Internal stop codon in {record.id} at pos {pos}.")
        if trailing_stops > 0:
            record.seq = Seq(seq_str[:idx])
        records_out.append(record)
    SeqIO.write(records_out, output_path, "fasta")
    return output_path


def split_fasta_pure(ref_fasta: Path, outdir: Path) -> List[Path]:
    """按 CDS 记录拆分成单基因 fasta (与 capheine_pipeline 同逻辑)"""
    gene_fastas = []
    prefix = ref_fasta.stem
    for record in SeqIO.parse(ref_fasta, "fasta"):
        safe_id = re.sub(r"[^a-zA-Z0-9_]", "_", record.id)
        out_name = f"{prefix}.part_{safe_id}.fasta"
        out_path = outdir / out_name
        SeqIO.write([record], out_path, "fasta")
        gene_fastas.append(out_path)
    return gene_fastas


def extract_cds_from_gbk(gbk_path: Path, out_path: Path) -> bool:
    """从 GenBank 提取 CDS 序列 → ref_cds.fasta (多 CDS 合并, 每条独立记录)"""
    records_out = []
    for rec in SeqIO.parse(gbk_path, "genbank"):
        for i, feat in enumerate(rec.features):
            if feat.type == "CDS" and feat.location is not None:
                try:
                    cds_seq = feat.extract(rec.seq)
                except Exception:
                    continue
                if len(cds_seq) % 3 != 0:
                    continue  # 非全长 CDS, 跳过
                gene_id = feat.qualifiers.get("gene", [f"CDS{i}"])[0]
                pid = feat.qualifiers.get("protein_id", [f"{rec.id}_{gene_id}"])[0]
                records_out.append(SeqRecord(cds_seq, id=f"{gene_id}_{pid}", description=""))
    if not records_out:
        return False
    SeqIO.write(records_out, out_path, "fasta")
    return True


# ------------------------------------------------------------- core flow ---

def find_reference_cds(work_dir: Path, explicit: Optional[str] = None,
                       accession: Optional[str] = None) -> Optional[Path]:
    """定位参考 CDS fasta, 优先级:
    1) 显式 --ref-cds / datasets.annotation_gb (gb/gbk/fasta 均可, gb 自动提取 CDS)
    2) select/capheine/ref_cds.fasta (capheine 产物)
    3) select/capheine/*noStopCodons.fasta
    4) data/ncbi_ref + work_dir 的 *.gbk/*.gbff/*.gb → 提取 CDS
    5) accession 在线下载 (参考 analysis 的 efetch 逻辑) → 提取 CDS
    """
    if explicit and Path(explicit).exists():
        ex = Path(explicit)
        if ex.suffix.lower() in (".gb", ".gbk", ".gbff"):
            out = work_dir / "time" / "gene_split" / "_ref_cds" / "ref_cds_from_gbk.fasta"
            out.parent.mkdir(parents=True, exist_ok=True)
            if extract_cds_from_gbk(ex, out):
                return out
            return None
        return ex
    for cand in [
        work_dir / "select" / "capheine" / "ref_cds.fasta",
        work_dir / "select" / "capheine" / "ref_cds-noStopCodons.fasta",
        work_dir / "select" / "capheine" / "removeterminalstopcodon",
    ]:
        if cand.is_file():
            return cand
        if cand.is_dir():
            for f in sorted(cand.glob("*.fasta")):
                return f
    # gbk 提取 (含 .gb: virus-annotations 产物是 .gb)
    gbk_dirs = [
        work_dir / "data" / "ncbi_ref",
        work_dir,
    ]
    for d in gbk_dirs:
        for gbk in sorted(d.glob("*.gbk")) + sorted(d.glob("*.gbff")) + sorted(d.glob("*.gb")):
            # 输出到子目录 _ref_cds, 避免被 gene_split 当成基因比对 (glob *.fasta 不递归)
            out = work_dir / "time" / "gene_split" / "_ref_cds" / "ref_cds_from_gbk.fasta"
            out.parent.mkdir(parents=True, exist_ok=True)
            if extract_cds_from_gbk(gbk, out):
                return out
    # 兜底: accession 在线下载 (无任何本地注释时)
    if accession:
        try:
            from utils.gb_fetch import get_or_fetch_cds_fasta
            out = work_dir / "time" / "gene_split" / "_ref_cds" / "ref_cds_from_gbk.fasta"
            cache = work_dir / "data" / "ncbi_ref"
            got = get_or_fetch_cds_fasta(accession, cache, out)
            if got:
                return Path(got)
        except Exception:
            pass
    return None


def find_unaligned_fasta(work_dir: Path, alignment: Optional[str] = None) -> Optional[Path]:
    """定位未比对全序列 fasta:
    1) 显式 --unaligned
    2) work_dir 下的非比对 fasta (排除 mafft.aln / -aligned / 产物目录)
    """
    if alignment and Path(alignment).exists():
        return Path(alignment)
    for f in sorted(work_dir.iterdir()):
        if f.suffix.lower() in (".fasta", ".fa", ".fas"):
            name = f.name.lower()
            if any(k in name for k in ("aln", "aligned", "tree", "ref_cds", "virna", "filtered", "noStop", "-nodups")):
                continue
            # 只取第一条序列就够判断: 含 gap 的是比对, 不含的是未比对
            try:
                rec = next(SeqIO.parse(f, "fasta"))
                if "-" not in str(rec.seq) and "." not in str(rec.seq):
                    return f
            except StopIteration:
                continue
    return None


def generate_gene_alignments(ref_cds: Path, unaligned: Path, out_dir: Path,
                             genetic_code: str = "1") -> List[Path]:
    """核心: 去终止密码子 → 拆分基因 → 每基因 cawlign → 干净命名 <gene>.fasta

    返回生成的基因比对文件列表 (已存在则复用)。
    """
    cawlign_bin = _find_tool("cawlign")
    if not cawlign_bin:
        raise RuntimeError("cawlign 不在 PATH (需要 CAPHEINE 环境的 cawlign 二进制)")

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    tmp = out_dir / "_tmp"
    tmp.mkdir(parents=True, exist_ok=True)

    # 1. 去终止密码子
    ref_nostop = tmp / f"{ref_cds.stem}-noStopCodons.fasta"
    if not ref_nostop.exists():
        remove_terminal_stop_codon(ref_cds, ref_nostop, genetic_code)

    # 2. 拆分基因
    gene_fastas = split_fasta_pure(ref_nostop, tmp)

    # 3. 每基因 cawlign → 干净命名 <gene>.fasta
    from Bio.Data import CodonTable
    code_name = {"1": "Universal", "11": "Bacterial", "4": "Mold"}.get(genetic_code)
    code_arg = f"-c '{code_name}'" if code_name else ""

    generated: List[Path] = []
    for gene_fa in gene_fastas:
        # 从 part_XXX 提取基因名
        m = re.search(r"\.part_([A-Za-z0-9_]+)\.fasta$", gene_fa.name)
        gene_name = m.group(1) if m else gene_fa.stem
        clean_out = out_dir / f"{gene_name}.fasta"
        if clean_out.exists() and clean_out.stat().st_size > 0:
            logger.info(f"  复用已有基因比对: {clean_out.name}")
            generated.append(clean_out)
            continue
        aligned_tmp = tmp / f"{gene_name}-aligned.fasta"
        # 路径加引号 (历史坑: shell=True 下路径含空格被拆词)
        import shlex as _sh
        cmd = (f"cawlign -t codon -r {_sh.quote(str(gene_fa))} -f refmap -s BLOSUM62 {code_arg} "
               f"{_sh.quote(str(unaligned))} > {_sh.quote(str(aligned_tmp))}")
        logger.info(f"  cawlign {gene_name} ...")
        r = subprocess.run(cmd, shell=True, capture_output=True, text=True, encoding='utf-8', errors='replace')
        if r.returncode != 0 or not aligned_tmp.exists():
            logger.warning(f"  cawlign {gene_name} 失败: {(r.stderr or r.stdout)[-200:]}")
            continue
        # 过滤全 gap/N 序列 (cawlign 输出可能含空序列)
        records = []
        for rec in SeqIO.parse(aligned_tmp, "fasta"):
            seq_str = str(rec.seq).upper()
            if set(seq_str) <= {'-', 'N', '?', 'X'}:
                logger.warning(f"  cawlign {gene_name}: 丢弃空序列 {rec.id}")
                continue
            records.append(rec)
        SeqIO.write(records, clean_out, "fasta")
        generated.append(clean_out)
        logger.info(f"  → {clean_out.name}")

    # 清理临时目录
    shutil.rmtree(tmp, ignore_errors=True)
    return generated


# ------------------------------------------------------------------ CLI ---

def main():
    ap = argparse.ArgumentParser(description="capheine 产物 fallback: 自动生成基因分区比对 (供 gene_dating)")
    ap.add_argument("--work-dir", required=True, help="流程 work_dir")
    ap.add_argument("--unaligned", default=None, help="未比对全序列 fasta (默认自动在 work_dir 找)")
    ap.add_argument("--ref-cds", default=None, help="参考 CDS fasta (默认: capheine 产物 → gbk 提取)")
    ap.add_argument("--out", default=None, help="输出目录 (默认 <work_dir>/time/gene_split)")
    ap.add_argument("--code", default="1", help="遗传密码表 (默认 1 Universal)")
    ap.add_argument("--min-genes", type=int, default=1, help="最少成功基因数 (不足则退出码 2)")
    args = ap.parse_args()

    logging.basicConfig(level=logging.INFO, format="[%(name)s] %(message)s")
    work_dir = Path(args.work_dir)
    out_dir = Path(args.out) if args.out else work_dir / "time" / "gene_split"

    ref = find_reference_cds(work_dir, args.ref_cds)
    if not ref:
        logger.error("找不到参考 CDS: 无 ref_cds.fasta / gbk / --ref-cds。无法生成基因比对。")
        sys.exit(2)

    unal = find_unaligned_fasta(work_dir, args.unaligned)
    if not unal:
        logger.error("找不到未比对全序列 fasta (--unaligned 或 work_dir 根目录)。无法生成基因比对。")
        sys.exit(2)

    logger.info(f"参考 CDS: {ref} ({sum(1 for _ in SeqIO.parse(ref, 'fasta'))} 基因)")
    logger.info(f"未比对序列: {unal}")
    logger.info(f"输出: {out_dir}")
    generated = generate_gene_alignments(ref, unal, out_dir, args.code)
    logger.info(f"生成 {len(generated)} 个基因比对 → {out_dir}")
    for g in generated:
        logger.info(f"  {g.name}")
    if len(generated) < args.min_genes:
        logger.error(f"成功基因数 {len(generated)} < min-genes {args.min_genes}")
        sys.exit(2)
    sys.exit(0)


if __name__ == "__main__":
    main()
