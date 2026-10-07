"""scan_and_pair 配对回归测试。

241 真实数据实测发现: 样本名以 [._-]2 结尾时 (如 GS-2_S2_L001_R1_fastq.gz),
旧的正则直接在整名上搜 R2 标签, 会先命中样本名里的 `-2`, 把 R1 文件同时
配成 R2 (fastp 收到 -i R1 -I R1, PE 数据被静默污染)。
修复后只在样本名之后的部分查找读段标签。
"""
import importlib.util
import os
from pathlib import Path

import pytest

_SPEC = importlib.util.spec_from_file_location(
    "clean_data", Path(__file__).resolve().parents[1] / "clean-data.py")
clean_data = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(clean_data)

FastqScanner = clean_data.FastqScanner


def _pair(tmp_path, names):
    for n in names:
        (tmp_path / n).write_text("@r1\nACGT\n+\nIIII\n")
    paired = {}
    for sample, r1, r2 in FastqScanner.scan_and_pair(str(tmp_path)):
        paired[sample] = (os.path.basename(r1) if r1 else None,
                          os.path.basename(r2) if r2 else None)
    return paired


def test_gs2_name_regression(tmp_path):
    """样本名以 -2 结尾: R1 文件不得被误配成 R2 (241 GS-2 实测 bug)。"""
    paired = _pair(tmp_path, [
        "GS-2_S2_L001_R1_fastq.gz",
        "GS-2_S2_L001_R2_fastq.gz",
    ])
    assert paired == {"GS-2": ("GS-2_S2_L001_R1_fastq.gz",
                               "GS-2_S2_L001_R2_fastq.gz")}


def test_name_with_dash1(tmp_path):
    """样本名以 -1 结尾: r1 搜索同样不得命中样本名里的 -1。"""
    paired = _pair(tmp_path, [
        "QH-1_S2_L001_R1_fastq.gz",
        "QH-1_S2_L001_R2_fastq.gz",
    ])
    assert paired == {"QH-1": ("QH-1_S2_L001_R1_fastq.gz",
                               "QH-1_S2_L001_R2_fastq.gz")}


def test_illumina_full_form(tmp_path):
    """标准 Illumina 全名 (含 S/L 标签): 正确配对, 不误判 `_001` 或 `_S2` 为读段标签。"""
    paired = _pair(tmp_path, [
        "NX-4_S2_L001_R1_001.fastq.gz",
        "NX-4_S2_L001_R2_001.fastq.gz",
    ])
    assert len(paired) == 1
    (r1, r2), = paired.values()
    assert r1 == "NX-4_S2_L001_R1_001.fastq.gz"
    assert r2 == "NX-4_S2_L001_R2_001.fastq.gz"


def test_classic_suffixes(tmp_path):
    paired = _pair(tmp_path, ["S1_R1.fq.gz", "S1_R2.fq.gz"])
    assert paired == {"S1": ("S1_R1.fq.gz", "S1_R2.fq.gz")}


def test_single_end_fallback(tmp_path):
    paired = _pair(tmp_path, ["SE-sample.fastq.gz"])
    assert paired == {"SE-sample": ("SE-sample.fastq.gz", None)}
