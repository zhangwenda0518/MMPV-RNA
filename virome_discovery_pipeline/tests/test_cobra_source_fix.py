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
