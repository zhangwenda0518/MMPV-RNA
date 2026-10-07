#!/usr/bin/env python3
"""
test_cobra_source_fix.py — COBRA 目录错配与 .ok 旧标记修复的回归测试
====================================================================
背景 (241/246 双机实测):
  1. 编排器把 02b_Filter 传给 COBRA raw 模式 (需要 02a_Identification) → 8/8 样本失败
  2. 失败只 warning → .ok 照写; 空壳样本目录又骗过产物校验 → 下轮旧标记跳过
  3. taxonomy 产物校验查陈旧路径 Votus.integrated/ (现行 integrated/)

运行: python -m pytest virome_discovery_pipeline/tests/test_cobra_source_fix.py -q
"""
from pathlib import Path
import shutil
from types import SimpleNamespace

import pytest

from virome_pipeline import ViromePipeline


def make_pipe(tmp_path: Path, *, ident=True, filt=True, coassembly=False):
    """构造最小 ViromePipeline (绕过 __init__), d 指向合成目录树"""
    root = tmp_path / "proj"
    d = {
        "root": root,
        "ident": root / "02a_Identification",
        "filter": root / "02b_Filter",
        "cobra": root / "03a_COBRA",
        "taxonomy": root / "05_Taxonomy",
        "_cobra_dirs": [root / "03a_COBRA"],
    }
    for k, v in d.items():
        if k not in ("root", "_cobra_dirs"):
            v.mkdir(parents=True, exist_ok=True)
    # 目录存在性必须与旗标一致 (ident=False = 02a 目录整个不存在)
    if not ident:
        shutil.rmtree(d["ident"], ignore_errors=True)
    if not filt:
        shutil.rmtree(d["filter"], ignore_errors=True)
    if ident:
        (d["ident"] / "5.virus_identification").mkdir(parents=True, exist_ok=True)
        (d["ident"] / "5.virus_identification" / "S1_virus.all.candidate.fasta").write_text(">a\nACGT\n")
    if filt:
        (d["filter"] / "S1").mkdir(parents=True, exist_ok=True)
        (d["filter"] / "S1" / "S1.virus.candidate_cdd_filtered.fasta").write_text(">a\nACGT\n")
    pipe = object.__new__(ViromePipeline)
    pipe.d = d
    pipe.args = SimpleNamespace(coassembly=coassembly)
    pipe.layout = "legacy"
    return pipe


# ── 修复 1: 候选来源 ──

def test_cobra_source_prefers_02a(tmp_path):
    """02a 存在时必须选 02a — 即使 02b 也在、即使 virus_mode=filter"""
    pipe = make_pipe(tmp_path, ident=True, filt=True)
    pipe.args.virus_mode = "filter"
    src, warn = pipe._cobra_candidate_source()
    assert src == pipe.d["ident"]
    assert warn is None


def test_cobra_source_falls_back_to_02b_with_warning(tmp_path):
    """02a 缺失时回退 02b 但必须告警 (不能静默错配)"""
    pipe = make_pipe(tmp_path, ident=False, filt=True)
    src, warn = pipe._cobra_candidate_source()
    assert src == pipe.d["filter"]
    assert warn and "02a" in warn


def test_cobra_source_both_missing_warns(tmp_path):
    pipe = make_pipe(tmp_path, ident=False, filt=False)
    src, warn = pipe._cobra_candidate_source()
    assert src == pipe.d["ident"]
    assert warn


# ── 2026-10-07 契约恢复: COBRA 主链吃 filter 结果 (镜像) ──

def make_pipe_v2(tmp_path, *, ident=True, filt=True, mirror=False, graph=False):
    """make_pipe 扩展: 支持 filter 主链镜像与 graph 派生候选目录"""
    import logging
    pipe = make_pipe(tmp_path, ident=ident, filt=filt)
    pipe.log = logging.getLogger("test_cobra_source")
    if mirror:
        m = pipe.d["filter"] / "candidates" / "S1"
        m.mkdir(parents=True, exist_ok=True)
        (m / "S1_virus.all.candidate.fasta").write_text(">a")
    if graph:
        g = pipe.d.setdefault("graph", pipe.d["root"] / "02b_Filter" / "graph")
        gc = g / "candidates" / "S1"
        gc.mkdir(parents=True, exist_ok=True)
        (gc / "S1_virus.all.candidate.fasta").write_text(">a")
    return pipe


def test_cobra_source_prefers_filter_mirror(tmp_path):
    """主链契约恢复 (用户裁定 2026-10-07): filter 镜像存在时 COBRA 必须吃 filter 结果"""
    pipe = make_pipe_v2(tmp_path, mirror=True, graph=True)
    src, warn = pipe._cobra_candidate_source()
    assert src == pipe.d["filter"] / "candidates"
    assert warn is None


def test_cobra_source_graph_fallback_when_no_filter(tmp_path):
    """filter 未跑 (无镜像) 时兜底 graph 派生候选"""
    pipe = make_pipe_v2(tmp_path, filt=False, graph=True)
    src, warn = pipe._cobra_candidate_source()
    assert src == pipe.d["graph"] / "candidates"
    assert warn is None


def test_cobra_source_ident_when_no_mirror_no_graph(tmp_path):
    pipe = make_pipe_v2(tmp_path)
    src, _ = pipe._cobra_candidate_source()
    assert src == pipe.d["ident"]


# ── 修复 2/3: .ok 旧标记不再骗过产物校验 ──

def test_validated_cobra_rejects_empty_sample_dirs(tmp_path):
    """失败残留的空壳样本目录 (无任何文件) 不能让 cobra 校验通过"""
    pipe = make_pipe(tmp_path)
    sample = pipe.d["cobra"] / "S1"
    sample.mkdir()  # 空壳目录 — 失败运行的典型残留
    assert pipe._stage_validated("cobra") is False


def test_validated_cobra_accepts_real_products(tmp_path):
    pipe = make_pipe(tmp_path)
    sample = pipe.d["cobra"] / "S1"
    sample.mkdir()
    (sample / "S1.cobra.extended.fa").write_text(">a\nACGT\n")
    assert pipe._stage_validated("cobra") is True


def test_validated_cobra_true_when_coassembly(tmp_path):
    """co-assembly 模式 COBRA 有意跳过, 信任标记"""
    pipe = make_pipe(tmp_path, coassembly=True)
    assert pipe._stage_validated("cobra") is True


# ── 修复 3b: taxonomy 校验路径 ──

def test_validated_taxonomy_accepts_integrated_path(tmp_path):
    """现行产物路径 integrated/final_integrated_classification.tsv"""
    pipe = make_pipe(tmp_path)
    f = pipe.d["taxonomy"] / "integrated" / "final_integrated_classification.tsv"
    f.parent.mkdir(parents=True, exist_ok=True)
    f.write_text("x" * 200)
    assert pipe._stage_validated("taxonomy") is True


def test_validated_taxonomy_legacy_votus_fallback(tmp_path):
    """历史旧名 Votus.integrated/ 仍兜底可用"""
    pipe = make_pipe(tmp_path)
    f = pipe.d["taxonomy"] / "Votus.integrated" / "final_integrated_classification.tsv"
    f.parent.mkdir(parents=True, exist_ok=True)
    f.write_text("x" * 200)
    assert pipe._stage_validated("taxonomy") is True


def test_validated_taxonomy_rejects_tiny_file(tmp_path):
    pipe = make_pipe(tmp_path)
    f = pipe.d["taxonomy"] / "integrated" / "final_integrated_classification.tsv"
    f.parent.mkdir(parents=True, exist_ok=True)
    f.write_text("x")  # <100 字节, 视为空跑
    assert pipe._stage_validated("taxonomy") is False


def test_build_cobra_mirror_passthrough_of_final_chain(tmp_path):
    """镜像 = filter 最终链原样搬运 (graph 对无证据结果的过滤已由
    apply_graph_adjustments 的 demote 在最终链上完成; 镜像不做额外增删)"""
    import os
    import filter_virus as fv
    outdir = tmp_path / "d_filter" / "S1"
    (outdir / "uniprot").mkdir(parents=True)
    base = "S1_virus.all.candidate"
    (outdir / "uniprot" / f"{base}_uniprot.fasta").write_text(">a\nACGT\n>b\nACGT\n>c\nACGT\n")
    os.symlink(outdir / "uniprot" / f"{base}_uniprot.fasta",
               outdir / f"{base}_filtered.fasta")
    fv.build_cobra_mirror(str(outdir), base, "filter", split_by_prefix=False)
    m = outdir / "candidates" / "S1" / "S1_virus.all.candidate.fasta"
    text = m.read_text()
    assert ">a" in text and ">b" in text and ">c" in text
    assert m.name == "S1_virus.all.candidate.fasta"
