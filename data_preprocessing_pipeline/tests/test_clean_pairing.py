"""scan_and_pair 配对回归测试。

241 真实数据实测发现: 样本名以 [._-]2 结尾时 (如 GS-2_S2_L001_R1_fastq.gz),
旧的正则直接在整名上搜 R2 标签, 会先命中样本名里的 `-2`, 把 R1 文件同时
配成 R2 (fastp 收到 -i R1 -I R1, PE 数据被静默污染)。
修复后只在样本名之后的部分查找读段标签。

多 lane 修复 (2026-10-13): 同一样本跨多个 lane 时旧实现只取第一对 R1/R2,
其余 lane 整对丢弃; 只有 R1 多文件时又同样本名拆多任务并发抢写同一输出。
现按 lane 分桶配对, 跨 lane 时任务键加 lane 后缀; 仍有歧义 (同 lane 多个
R1 / 配对后多余文件) → raise ValueError (fail-loud)。
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


# ── 多 lane / 歧义 (2026-10-13) ──────────────────────────────────

def test_multi_lane_pairs_all_kept(tmp_path):
    """多 lane 样本: X_L001/X_L002 各自成对 (旧实现只取第一对, L002 整对丢弃)。"""
    paired = _pair(tmp_path, [
        "X_L001_R1.fastq.gz", "X_L001_R2.fastq.gz",
        "X_L002_R1.fastq.gz", "X_L002_R2.fastq.gz",
    ])
    assert paired == {
        "X_L001": ("X_L001_R1.fastq.gz", "X_L001_R2.fastq.gz"),
        "X_L002": ("X_L002_R1.fastq.gz", "X_L002_R2.fastq.gz"),
    }


def test_multi_lane_r1_only_unique_keys(tmp_path):
    """只有 R1 的多 lane 文件: 加 lane 后缀, 不再同样本名拆多任务抢写同一输出。"""
    paired = _pair(tmp_path, ["X_L001_R1.fastq.gz", "X_L002_R1.fastq.gz"])
    assert paired == {
        "X_L001": ("X_L001_R1.fastq.gz", None),
        "X_L002": ("X_L002_R1.fastq.gz", None),
    }


def test_multi_lane_with_unlabelled_pair(tmp_path):
    """多 lane 且其中一个 lane 无 lane 标记: 两组都保留, 键分别不加/加后缀。"""
    paired = _pair(tmp_path, [
        "X_R1.fastq.gz", "X_R2.fastq.gz",
        "X_L001_R1.fastq.gz", "X_L001_R2.fastq.gz",
    ])
    assert paired == {
        "X": ("X_R1.fastq.gz", "X_R2.fastq.gz"),
        "X_L001": ("X_L001_R1.fastq.gz", "X_L001_R2.fastq.gz"),
    }


def test_same_lane_multiple_r1_raises(tmp_path):
    """同 lane 两个 R1 (R1/R1_001 混合命名) → ValueError, 消息列出冲突文件。"""
    paired_dir = tmp_path
    for n in ("X_L001_R1.fastq.gz", "X_L001_1.fastq.gz", "X_L001_R2.fastq.gz"):
        (paired_dir / n).write_text("@r1\nACGT\n+\nIIII\n")
    with pytest.raises(ValueError) as ei:
        FastqScanner.scan_and_pair(str(paired_dir))
    msg = str(ei.value)
    assert "X_L001_R1.fastq.gz" in msg and "X_L001_1.fastq.gz" in msg


def test_pairing_with_extra_file_raises(tmp_path):
    """同 lane 配对成功但另有无法归属的多余文件 → ValueError (不再静默丢弃)。"""
    for n in ("X_1.fq.gz", "X_2.fq.gz", "X.fq.gz"):
        (tmp_path / n).write_text("@r1\nACGT\n+\nIIII\n")
    with pytest.raises(ValueError) as ei:
        FastqScanner.scan_and_pair(str(tmp_path))
    assert "X.fq.gz" in str(ei.value)
