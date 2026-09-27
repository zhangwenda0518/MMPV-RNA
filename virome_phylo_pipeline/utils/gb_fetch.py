"""gb_fetch.py — 参考 gb 下载（NCBI），供 capheine/gene_dating 的 CDS 参考兜底

参考 analysis 管线的两种下载实现:
  - fetch_gb_via_efetch (batch_plot_virus_depth.py): efetch CLI (服务器 .pixi 有 v16.2)
  - download_genbank (virus_auto_pipeline.py / consensus_extract.py): Bio.Entrez.efetch

策略: efetch CLI 优先 (快, 无 API key 限制), Entrez 兜底 (CLI 缺失/失败时)。
已存在且非空的 gb 文件直接复用 (幂等)。
"""
from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path

DEFAULT_EMAIL = "mmvp@example.com"


def download_gb(accession: str, out_dir, timeout: int = 120) -> str | None:
    """下载 {accession}.gb 到 out_dir; 已存在非空则复用; 失败返回 None。"""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    out_file = out_dir / f"{accession}.gb"
    if out_file.exists() and out_file.stat().st_size > 500:
        return str(out_file)

    # ① efetch CLI (NCBI E-utilities, 服务器 .pixi 已装 v16.2)
    efetch_bin = shutil.which("efetch")
    if efetch_bin:
        cmd = [efetch_bin, "-db", "nuccore", "-id", accession, "-format", "gb"]
        try:
            # 2026-09-15: subprocess 把字节直接写到 fileno(), 文本模式包装器上
            # 的 encoding 其实不生效 (且与 stdout=f 的语义不一致) -> 改二进制打开。
            with open(out_file, "wb") as f:
                subprocess.run(cmd, stdout=f, stderr=subprocess.PIPE,
                               check=True, timeout=timeout)
            if out_file.exists() and out_file.stat().st_size > 500:
                return str(out_file)
            if out_file.exists():
                out_file.unlink()  # 下载到空/太小文件视为失败
        except (subprocess.CalledProcessError, subprocess.TimeoutExpired, OSError):
            if out_file.exists():
                out_file.unlink()

    # ② Bio.Entrez 兜底 (CLI 缺失/失败)
    try:
        from Bio import Entrez
        Entrez.email = os.environ.get("NCBI_EMAIL", DEFAULT_EMAIL)
        Entrez.tool = "MMPV-PhyloPipeline"
        with Entrez.efetch(db="nucleotide", id=accession,
                           rettype="gb", retmode="text") as handle:
            content = handle.read()
        if "Item not found" in content or content.strip().startswith("Error"):
            return None
        out_file.write_text(content, encoding="utf-8")
        if out_file.stat().st_size > 500:
            return str(out_file)
        out_file.unlink()
    except Exception:
        if out_file.exists():
            out_file.unlink()
    return None


def extract_cds_to_fasta(gb_path: str, out_fasta: str) -> bool:
    """从 gb/gbk 提取全长 CDS → fasta; 成功返回 True。"""
    from Bio import SeqIO
    from Bio.SeqRecord import SeqRecord
    try:
        records_out = []
        for rec in SeqIO.parse(gb_path, "genbank"):
            for i, feat in enumerate(rec.features):
                if feat.type == "CDS" and feat.location is not None:
                    try:
                        cds_seq = feat.extract(rec.seq)
                    except Exception:
                        continue
                    if len(cds_seq) % 3 != 0:
                        continue  # 非全长 CDS, 跳过
                    gene_id = feat.qualifiers.get("gene", [f"CDS{i}"])[0]
                    pid = feat.qualifiers.get("protein_id",
                                              [f"{rec.id}_{gene_id}"])[0]
                    records_out.append(
                        SeqRecord(cds_seq, id=f"{gene_id}_{pid}", description=""))
        if not records_out:
            return False
        Path(out_fasta).parent.mkdir(parents=True, exist_ok=True)
        SeqIO.write(records_out, out_fasta, "fasta")
        return True
    except Exception:
        return False


def get_or_fetch_cds_fasta(accession: str, cache_dir, out_fasta: str) -> str | None:
    """一站式: 下载 gb (缓存) → 提取 CDS → 返回 fasta; 任一环节失败返回 None。"""
    gb = download_gb(accession, cache_dir)
    if not gb:
        return None
    if not extract_cds_to_fasta(gb, out_fasta):
        return None
    return out_fasta


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(description="下载参考 gb 并提取 CDS")
    ap.add_argument("accession", help="NCBI accession, 如 OR489165.1")
    ap.add_argument("--out", default=".", help="输出目录")
    ap.add_argument("--fasta", default=None, help="CDS fasta 输出 (缺省 {out}/{accession}.cds.fasta)")
    args = ap.parse_args()
    fa = args.fasta or os.path.join(args.out, f"{args.accession}.cds.fasta")
    r = get_or_fetch_cds_fasta(args.accession, args.out, fa)
    print(f"OK → {r}" if r else "FAILED")
