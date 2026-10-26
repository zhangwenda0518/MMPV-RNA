"""build_host_pipeline 的 eve-mask 阶段测试（默认关闭 + 上下游衔接）。

覆盖:
  · 阶段注册顺序 genome-down -> eve-mask -> hostdb -> eve-report
  · 默认关闭：不给 `--eve-mask`/掩码源时 `--stage all` 的计划里没有它（与 v3.1 行为一致）
  · 任一掩码源（BED / verdict / masked-genome / eve-scan）自动启用
  · 显式 `--stage eve-mask` 却无掩码源 → 退出并提示
  · run_eve_mask: N 化、长度守恒、完整头行保留、report/meta 产物、verdict 过滤
  · 区间与 contig 不匹配时**拒绝**生成「假掩码」库
  · resolve_inputs: 原始库 = <work>/hostdb；掩码 = <work>/hostdb_evemask
  · hostdb 阶段确实改用掩码基因组建库，并在库目录写 db_provenance.json

Run:  MMPV_TEST_TMPDIR=<可写目录> python -m pytest \\
          public_metadata_pipeline/tests/test_eve_mask_stage.py -q -p no:cacheprovider
"""
import hashlib
import json
import os
import shutil
import sys
import tempfile
import uuid
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import build_host_pipeline as bhp  # noqa: E402

# chr1 160 bp / chr2 80 bp
GENOME = (
    ">chr1 Lycium barbarum chromosome 1\n"
    + "ACGT" * 40 + "\n"
    + ">chr2\n"
    + "TTTTGGGGCCCCAAAA" * 5 + "\n"
)

VERDICT = (
    "locus\tverdict\tbest_viral_id\tbest_viral_bs\tbest_plant_id\tbest_plant_bs\n"
    "chr1:11-30\tviral_supported\tCaulimoviridae_sp\t120\tViridiplantae\t80\n"
    "chr1:51-60\thost_like\t\tNaN\tViridiplantae\t300\n"
    "chr1:101-130\tviral_supported\tBadnavirus_sp\t150\tViridiplantae\t60\n"
    "chr2:6-15\tundetermined\t\tNaN\tViridiplantae\t90\n"
    "chr2:40-60\tviral_supported\tGeminiviridae_sp\t110\tViridiplantae\t70\n"
)


@pytest.fixture
def workdir():
    # 受限环境 (DSH 沙箱: mkdtemp 建出的目录带受限 DACL) 用 MMPV_TEST_TMPDIR 指到可写目录
    base = os.environ.get("MMPV_TEST_TMPDIR")
    if base:
        d = Path(base) / f"eve_stage_{uuid.uuid4().hex[:8]}"
        d.mkdir(parents=True, exist_ok=True)
    else:
        d = Path(tempfile.mkdtemp(prefix="eve_stage_"))
    yield d
    shutil.rmtree(d, ignore_errors=True)


def _write_inputs(workdir):
    g = workdir / "genome.fa"
    g.write_text(GENOME, encoding="utf-8")
    v = workdir / "verdict.tsv"
    v.write_text(VERDICT, encoding="utf-8")
    return g, v


def _read_fasta(path):
    headers, seqs, cur = [], {}, None
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            cur = line[1:].split()[0]
            headers.append(line[1:])
            seqs[cur] = ""
        elif line.strip():
            seqs[cur] += line.strip()
    return headers, seqs


def make_args(workdir, extra=(), stage="hostdb"):
    parser = bhp.build_parser()
    argv = ["--species", "Lycium barbarum", "--taxid", "112863",
            "--stage", stage, "--work-dir", str(workdir / "work")] + list(extra)
    return parser.parse_args(argv)


def make_pipeline(workdir, extra=(), stage="hostdb"):
    return bhp.BuildHostPipeline(make_args(workdir, extra, stage))


# ── 阶段注册与默认关闭 ─────────────────────────────
def test_stage_registered_between_genome_and_hostdb():
    assert bhp.BuildHostPipeline.STAGES == \
        ["genome-down", "eve-mask", "hostdb", "eve-report"]
    assert "eve-mask" in bhp.BuildHostPipeline.STAGE_DESC
    assert "eve-report" in bhp.BuildHostPipeline.STAGE_DESC


def test_defaults_are_off(workdir):
    args = make_args(workdir, stage="all")
    assert args.eve_mask is False
    assert args.eve_bed == [] and args.eve_verdict == []
    assert args.eve_scan is False and args.masked_genome is None
    assert bhp.eve_mask_enabled(args) is False


def test_stage_all_skips_mask_when_disabled(workdir):
    args = make_args(workdir, stage="all")
    assert bhp.plan_stages(args) == ["genome-down", "hostdb", "eve-report"], \
        "默认关闭时 --stage all 不应包含 eve-mask（eve-report 恒定保留）"


@pytest.mark.parametrize("extra", [
    ["--eve-mask"],
    ["--eve-bed", "x.bed"],
    ["--eve-verdict", "x.tsv"],
    ["--masked-genome", "x.fasta"],
    ["--eve-scan"],
])
def test_any_mask_source_enables_stage(workdir, extra):
    assert bhp.eve_mask_enabled(make_args(workdir, extra, stage="all")) is True


def test_stage_all_includes_mask_when_enabled(workdir):
    args = make_args(workdir, ["--eve-verdict", "x.tsv"], stage="all")
    assert bhp.plan_stages(args) == \
        ["genome-down", "eve-mask", "hostdb", "eve-report"]


def test_explicit_stage_without_source_errors(workdir):
    args = make_args(workdir, stage="eve-mask")
    with pytest.raises(SystemExit) as ei:
        bhp.plan_stages(args)
    assert "需要掩码源" in str(ei.value)


# ── run_eve_mask ──────────────────────────────────
def test_run_eve_mask_writes_masked_fasta_report_and_meta(workdir):
    g, v = _write_inputs(workdir)
    p = make_pipeline(workdir, ["--genome-fasta", str(g), "--eve-verdict", str(v),
                                "--eve-pad", "2"], stage="eve-mask")
    assert p.run_eve_mask() is True

    out = workdir / "work" / "genome" / "all.genome.masked.fasta"
    assert out.is_file()
    headers, seqs = _read_fasta(out)
    assert headers[0].startswith("chr1 Lycium barbarum"), "完整头行要保留"
    assert len(seqs["chr1"]) == 160 and len(seqs["chr2"]) == 80, "掩码不能改变长度"
    # chr1: 11-30 + 101-130，pad=2 → 24 + 34 = 58；host_like 与 undetermined 不掩
    assert seqs["chr1"].count("N") == 58
    # chr2: 40-60 加 pad=2 → 25
    assert seqs["chr2"].count("N") == 25

    rep = (out.parent / (out.name + ".report.tsv")).read_text(encoding="utf-8").splitlines()
    assert rep[0].split("\t") == ["genome", "contig", "length", "n_loci",
                                  "masked_bp", "masked_frac"]
    assert {r.split("\t")[1]: int(r.split("\t")[4]) for r in rep[1:]} == \
        {"chr1": 58, "chr2": 25}

    meta = json.loads((out.parent / (out.name + ".meta.json")).read_text(encoding="utf-8"))
    assert meta["schema"] == "mmpv-eve-mask/1"
    assert meta["totals"]["bp"] == 240
    assert meta["totals"]["masked_bp"] == 83
    assert meta["totals"]["masked_frac"] == pytest.approx(83 / 240, rel=1e-6)
    assert meta["params"]["pad"] == 2
    assert meta["mask_sources"]["verdict_values"] == ["viral_supported"]
    assert meta["advice"], "必须给出 masked_frac 判读建议"
    assert p.ckpt.is_done("eve-mask")


def test_verdict_values_can_include_undetermined(workdir):
    g, v = _write_inputs(workdir)
    p = make_pipeline(workdir, ["--genome-fasta", str(g), "--eve-verdict", str(v),
                                "--eve-verdict-values", "viral_supported,undetermined"],
                      stage="eve-mask")
    assert p.run_eve_mask() is True
    _, seqs = _read_fasta(workdir / "work" / "genome" / "all.genome.masked.fasta")
    # 默认 20+30+21 = 71（chr2:40-60 是 1 基闭区间 = 21 bp）；
    # 加上 undetermined chr2:6-15 (10 bp) 后 81
    assert sum(s.count("N") for s in seqs.values()) == 81


def test_run_eve_mask_accepts_bed(workdir):
    g, _ = _write_inputs(workdir)
    bed = workdir / "m.bed"
    bed.write_text("chr1\t0\t10\nchr2\t0\t5\n", encoding="utf-8")
    p = make_pipeline(workdir, ["--genome-fasta", str(g), "--eve-bed", str(bed)],
                      stage="eve-mask")
    assert p.run_eve_mask() is True
    _, seqs = _read_fasta(workdir / "work" / "genome" / "all.genome.masked.fasta")
    assert seqs["chr1"][:10] == "N" * 10
    assert seqs["chr1"][10:] == ("ACGT" * 40)[10:]
    assert seqs["chr2"][:5] == "NNNNN"
    assert seqs["chr2"][5] != "N"


def test_contig_mismatch_refuses_to_fake_masked_db(workdir):
    g, _ = _write_inputs(workdir)
    v = workdir / "wrong.tsv"
    v.write_text("locus\tverdict\nchrX:1-10\tviral_supported\n", encoding="utf-8")
    p = make_pipeline(workdir, ["--genome-fasta", str(g), "--eve-verdict", str(v)],
                      stage="eve-mask")
    assert p.run_eve_mask() is False, "区间没落在任何 contig 上必须报错"
    out = workdir / "work" / "genome" / "all.genome.masked.fasta"
    assert not out.exists(), "不能留下一个其实没掩码的 masked FASTA 冒充掩码库"
    assert not p.ckpt.is_done("eve-mask")


def test_run_eve_mask_without_sources_fails(workdir):
    g, _ = _write_inputs(workdir)
    p = make_pipeline(workdir, ["--eve-mask", "--genome-fasta", str(g)], stage="eve-mask")
    assert p.run_eve_mask() is False


def test_run_eve_mask_accepts_pregenerated_masked_genome(workdir):
    g, _ = _write_inputs(workdir)
    masked = workdir / "pre.masked.fasta"
    masked.write_text(GENOME, encoding="utf-8")
    p = make_pipeline(workdir, ["--masked-genome", str(masked)], stage="eve-mask")
    assert p.run_eve_mask() is True
    fasta, _, is_masked = p.resolve_inputs()
    assert Path(fasta) == masked.resolve() and is_masked is True
    meta = json.loads((workdir / "pre.masked.fasta.meta.json").read_text(encoding="utf-8"))
    assert meta["mode"] == "provided"


# ── 上下游衔接: 用哪个 FASTA / 哪个库目录 ────────────
def test_resolve_inputs_default_is_raw_hostdb(workdir):
    g, _ = _write_inputs(workdir)
    p = make_pipeline(workdir, ["--genome-fasta", str(g)], stage="hostdb")
    fasta, dbdir, masked = p.resolve_inputs()
    assert Path(fasta) == g
    assert Path(dbdir) == workdir / "work" / "hostdb"
    assert masked is False


def test_resolve_inputs_masked_but_missing_fasta(workdir):
    g, v = _write_inputs(workdir)
    p = make_pipeline(workdir, ["--genome-fasta", str(g), "--eve-verdict", str(v)],
                      stage="hostdb")
    with pytest.raises(FileNotFoundError):
        p.resolve_inputs()


def test_resolve_inputs_uses_hostdb_evemask_when_masked(workdir):
    g, v = _write_inputs(workdir)
    p = make_pipeline(workdir, ["--genome-fasta", str(g), "--eve-verdict", str(v)],
                      stage="hostdb")
    out = workdir / "work" / "genome" / "all.genome.masked.fasta"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(GENOME, encoding="utf-8")
    fasta, dbdir, masked = p.resolve_inputs()
    assert Path(fasta) == out, "掩码文件存在时 hostdb 必须用掩码基因组"
    assert Path(dbdir) == workdir / "work" / "hostdb_evemask"
    assert masked is True


def test_hostdb_stage_builds_from_masked_genome_and_writes_provenance(workdir, monkeypatch):
    g, v = _write_inputs(workdir)
    p = make_pipeline(workdir, ["--genome-fasta", str(g), "--eve-verdict", str(v),
                                "--eve-pad", "2"], stage="eve-mask")
    assert p.run_eve_mask() is True

    calls = []

    def fake_run_cmd(cmd, *a, **k):
        calls.append(cmd)
        return 0

    monkeypatch.setattr(bhp, "run_cmd", fake_run_cmd)

    p2 = make_pipeline(workdir, ["--genome-fasta", str(g), "--eve-verdict", str(v)],
                       stage="hostdb")
    assert p2.run_hostdb() is True

    assert "all.genome.masked.fasta" in calls[0]
    assert "hostdb_evemask" in calls[0]
    assert getattr(p2, "_hostdb_hint", None) == str(workdir / "work" / "hostdb_evemask")

    prov = json.loads((workdir / "work" / "hostdb_evemask" /
                       "db_provenance.json").read_text(encoding="utf-8"))
    assert prov["masked"] is True
    assert prov["species"] == "Lycium barbarum" and prov["taxid"] == 112863
    assert prov["genome"]["path"].endswith("all.genome.masked.fasta")
    assert prov["eve_mask_meta"]["totals"]["masked_bp"] == 83
    assert prov["hostdb_tools"] == "kraken2,bowtie2,hisat2,minimap2"


def test_provenance_records_unmasked_library(workdir):
    g, _ = _write_inputs(workdir)
    p = make_pipeline(workdir, ["--genome-fasta", str(g)], stage="hostdb")
    prov_path = p.write_db_provenance(workdir / "rawdb", str(g), False)
    data = json.loads(Path(prov_path).read_text(encoding="utf-8"))
    assert data["masked"] is False and data["eve_mask_meta"] is None
    assert data["genome"]["md5"] == hashlib.md5(g.read_bytes()).hexdigest()
    assert data["schema"] == "mmpv-hostdb-provenance/1"
