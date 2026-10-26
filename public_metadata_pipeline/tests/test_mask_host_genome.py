"""mask_host_genome 的单元测试。

覆盖: BED 解析/合并/外扩、N 化与小写软掩码、verdict 表过滤、gzip 输入输出、
多基因组 --out-dir、报告列、contig 名不匹配告警、长度不变（坐标守恒）。

Run:  python -m pytest public_metadata_pipeline/tests/test_mask_host_genome.py -v
"""
import gzip
import os
import shutil
import sys
import tempfile
import uuid
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import mask_host_genome as mhg  # noqa: E402

# 2 条 contig: 40 bp 全 'A' / 30 bp 全 'C'
FASTA = ">ctg1 desc\n" + "A" * 40 + "\n>ctg2\n" + "C" * 30 + "\n"


@pytest.fixture
def workdir():
    # 默认用系统临时目录; 受限环境 (如 DSH 沙箱: mkdtemp 建出的目录带受限 DACL,
    # 子进程再写会 WinError 5) 可用 MMPV_TEST_TMPDIR 指向一个可写目录。
    base = os.environ.get("MMPV_TEST_TMPDIR")
    if base:
        d = Path(base) / f"mask_host_{uuid.uuid4().hex[:8]}"
        d.mkdir(parents=True, exist_ok=True)
    else:
        d = Path(tempfile.mkdtemp(prefix="mask_host_genome_"))
    yield d
    shutil.rmtree(d, ignore_errors=True)


def _read_fasta(path):
    return [(h, s.decode("ascii")) for h, s in mhg.iter_fasta(str(path))]


def _seqs(path):
    """contig(首个词) -> 序列；完整头行另由 _read_fasta 断言。"""
    return {h.split()[0]: s for h, s in _read_fasta(path)}


def test_mask_bed_basic(workdir):
    fa = workdir / "g.fasta"
    fa.write_text(FASTA, encoding="utf-8")
    bed = workdir / "m.bed"
    bed.write_text("ctg1\t10\t20\n", encoding="utf-8")
    out = workdir / "g.masked.fasta"

    rc = mhg.main(["-g", str(fa), "-b", str(bed), "-o", str(out), "--quiet"])
    assert rc == 0

    recs = _seqs(out)
    assert len(recs["ctg1"]) == 40, "掩码不能改变序列长度（坐标守恒）"
    assert recs["ctg1"][:10] == "A" * 10
    assert recs["ctg1"][10:20] == "N" * 10
    assert recs["ctg1"][20:] == "A" * 20
    assert recs["ctg2"] == "C" * 30, "未被区间覆盖的 contig 原样保留"

    rep = (workdir / "g.masked.fasta.report.tsv").read_text(encoding="utf-8").splitlines()
    assert rep[0].split("\t") == ["genome", "contig", "length", "n_loci",
                                  "masked_bp", "masked_frac"]
    rows = {r.split("\t")[1]: r.split("\t") for r in rep[1:]}
    assert rows["ctg1"][3] == "1" and rows["ctg1"][4] == "10"
    assert rows["ctg2"][3] == "0" and rows["ctg2"][4] == "0"


def test_mask_header_preserved(workdir):
    fa = workdir / "g.fasta"
    fa.write_text(FASTA, encoding="utf-8")
    recs = _read_fasta(fa)
    assert recs[0][0] == "ctg1 desc", "完整头行（含描述）必须原样写回"


def test_pad_and_merge_gap(workdir):
    fa = workdir / "g.fasta"
    fa.write_text(FASTA, encoding="utf-8")
    bed = workdir / "m.bed"
    # 两段相距 3 bp；--merge-gap 5 应合成一段；--pad 2 再外扩
    bed.write_text("ctg1\t10\t15\nctg1\t18\t20\n", encoding="utf-8")
    out = workdir / "o.fasta"
    assert mhg.main(["-g", str(fa), "-b", str(bed), "-o", str(out),
                     "--pad", "2", "--merge-gap", "5", "--quiet"]) == 0
    seq = _seqs(out)["ctg1"]
    assert seq[8:22] == "N" * 14, f"pad=2 + 合并后应为 8..22，实际 {seq}"


def test_lower_soft_mask(workdir):
    fa = workdir / "g.fasta"
    fa.write_text(FASTA, encoding="utf-8")
    bed = workdir / "m.bed"
    bed.write_text("ctg2\t0\t5\n", encoding="utf-8")
    out = workdir / "o.fasta"
    assert mhg.main(["-g", str(fa), "-b", str(bed), "-o", str(out),
                     "--mask-mode", "lower", "--quiet"]) == 0
    seq = _seqs(out)["ctg2"]
    assert seq[:5] == "ccccc" and seq[5:] == "C" * 25


def test_gzip_input_and_output(workdir):
    fa = workdir / "g.fasta.gz"
    with gzip.open(fa, "wt", encoding="utf-8") as fh:
        fh.write(FASTA)
    bed = workdir / "m.bed"
    bed.write_text("ctg1\t0\t4\n", encoding="utf-8")
    out = workdir / "o.fasta"
    assert mhg.main(["-g", str(fa), "-b", str(bed), "-o", str(out), "--quiet"]) == 0
    recs = _seqs(out)
    assert recs["ctg1"].startswith("NNNN") and len(recs["ctg1"]) == 40


def test_verdict_tsv_filters_values(workdir):
    fa = workdir / "g.fasta"
    fa.write_text(FASTA, encoding="utf-8")
    vd = workdir / "s2_verdict.tsv"
    # 1 基闭区间: ctg1:11-20 -> BED 10..20 (10 bp)
    vd.write_text(
        "locus\tverdict\tbest_viral_id\tbest_viral_bs\tbest_plant_id\tbest_plant_bs\n"
        "ctg1:11-20\tviral_supported\tv1\t100\tp1\t60\n"
        "ctg1:21-25\thost_like\tv2\t90\tp2\t200\n"
        "ctg1:31-35\tundetermined\t\t\t\t\n",
        encoding="utf-8")
    out = workdir / "o.fasta"
    assert mhg.main(["-g", str(fa), "-v", str(vd), "-o", str(out), "--quiet"]) == 0
    seq = _seqs(out)["ctg1"]
    assert seq[10:20] == "N" * 10
    assert seq[20:25] == "A" * 5, "host_like 不应进掩码"
    assert seq[30:35] == "A" * 5, "undetermined 默认不掩"


def test_verdict_tsv_aggressive_values(workdir):
    fa = workdir / "g.fasta"
    fa.write_text(FASTA, encoding="utf-8")
    vd = workdir / "s2_verdict.tsv"
    vd.write_text("locus\tverdict\n"
                  "ctg1:11-20\tviral_supported\n"
                  "ctg1:31-35\tundetermined\n", encoding="utf-8")
    out = workdir / "o.fasta"
    assert mhg.main(["-g", str(fa), "-v", str(vd), "-o", str(out),
                     "--verdict-values", "viral_supported,undetermined",
                     "--quiet"]) == 0
    seq = _seqs(out)["ctg1"]
    assert seq[10:20] == "N" * 10 and seq[30:35] == "N" * 5


def test_kingdom_loci_all_layout(workdir):
    """kingdom_loci_all.tsv.gz: 多一列 genome，列位靠表头定位。"""
    fa = workdir / "g.fasta"
    fa.write_text(FASTA, encoding="utf-8")
    vd = workdir / "kingdom_loci_all.tsv.gz"
    with gzip.open(vd, "wt", encoding="utf-8") as fh:
        fh.write("genome\tlocus\tverdict\tref_family\n"
                 "g\tctg2:1-10\tviral_supported\tCaulimoviridae\n"
                 "g\tctg2:11-20\thost_like\tViridiplantae\n")
    out = workdir / "o.fasta"
    assert mhg.main(["-g", str(fa), "-v", str(vd), "-o", str(out), "--quiet"]) == 0
    seq = _seqs(out)["ctg2"]
    assert seq[:10] == "N" * 10 and seq[10:] == "C" * 20


def test_multiple_genomes_out_dir(workdir):
    g1 = workdir / "ningxia.fasta"
    g2 = workdir / "heiguo.fa"
    g1.write_text(FASTA, encoding="utf-8")
    g2.write_text(">ctg1\n" + "G" * 12 + "\n", encoding="utf-8")
    bed = workdir / "m.bed"
    bed.write_text("ctg1\t2\t6\n", encoding="utf-8")
    outdir = workdir / "mask-genome"
    assert mhg.main(["-g", str(g1), str(g2), "-b", str(bed),
                     "--out-dir", str(outdir), "--quiet"]) == 0
    assert (outdir / "ningxia.masked.fasta").exists()
    assert (outdir / "heiguo.masked.fasta").exists()
    assert _seqs(outdir / "ningxia.masked.fasta")["ctg1"][2:6] == "NNNN"
    assert _seqs(outdir / "heiguo.masked.fasta")["ctg1"][2:6] == "NNNN"
    # 逐基因组报告不应串行
    rep = (outdir / "heiguo.masked.report.tsv").read_text(encoding="utf-8")
    assert "ningxia" not in rep


def test_unknown_contig_warns(workdir, capsys):
    fa = workdir / "g.fasta"
    fa.write_text(FASTA, encoding="utf-8")
    bed = workdir / "m.bed"
    bed.write_text("chrX\t0\t10\n", encoding="utf-8")
    out = workdir / "o.fasta"
    assert mhg.main(["-g", str(fa), "-b", str(bed), "-o", str(out)]) == 0
    err = capsys.readouterr().err
    assert "WARN" in err and "chrX" in err
    # 未匹配到时序列原样输出，不崩
    assert _seqs(out)["ctg1"] == "A" * 40


def test_requires_mask_source(workdir):
    fa = workdir / "g.fasta"
    fa.write_text(FASTA, encoding="utf-8")
    assert mhg.main(["-g", str(fa), "-o", str(workdir / "o.fasta")]) == 2


def test_out_required(workdir):
    fa = workdir / "g.fasta"
    fa.write_text(FASTA, encoding="utf-8")
    bed = workdir / "m.bed"
    bed.write_text("ctg1\t0\t5\n", encoding="utf-8")
    assert mhg.main(["-g", str(fa), "-b", str(bed)]) == 2


def test_interval_beyond_contig_is_clipped(workdir):
    fa = workdir / "g.fasta"
    fa.write_text(FASTA, encoding="utf-8")
    bed = workdir / "m.bed"
    bed.write_text("ctg1\t30\t999\n", encoding="utf-8")
    out = workdir / "o.fasta"
    assert mhg.main(["-g", str(fa), "-b", str(bed), "-o", str(out), "--quiet"]) == 0
    seq = _seqs(out)["ctg1"]
    assert seq[:30] == "A" * 30 and seq[30:] == "N" * 10 and len(seq) == 40


def test_merge_intervals_pure():
    assert mhg.merge_intervals([(5, 10), (1, 3)], 0) == [(1, 3), (5, 10)]
    assert mhg.merge_intervals([(1, 5), (5, 9)], 0) == [(1, 9)]
    assert mhg.merge_intervals([(1, 5), (7, 9)], 2) == [(1, 9)]
    assert mhg.merge_intervals([(1, 5), (6, 9)], 0) == [(1, 5), (6, 9)]
    assert mhg.merge_intervals([], 0) == []
