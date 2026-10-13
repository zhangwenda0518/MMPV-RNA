# -*- coding: utf-8 -*-
"""④ 变异指标统计的数值对拍测试 (2026-10-13, 修 10-10 V03)。

pysam 在部分开发机不可用 → 通过 sys.modules 注入桩后再导入模块,
只测纯计算函数 _site_pi_mean (Pi 公式本体), 不触碰 BAM 路径。
"""
import sys
import types
from collections import Counter
from pathlib import Path

import pytest

_HERE = Path(__file__).resolve().parent
_ANALYSIS = _HERE.parent
_REPO = _ANALYSIS.parent
for p in (str(_ANALYSIS), str(_REPO)):
    if p not in sys.path:
        sys.path.insert(0, p)

if "pysam" not in sys.modules:
    try:
        import pysam  # noqa: F401
    except ImportError:
        sys.modules["pysam"] = types.SimpleNamespace(AlignmentFile=None)

from batch_virus_variants import _site_pi_mean  # noqa: E402


def test_site_pi_matches_hand_calculation():
    """V03: Pi = n/(n-1)·(1−Σp²) 的位点均值 (scikit-allel), 拒绝旧的 1−p_max。

    手算对照:
      位点1 A=6,C=4 (n=10): 10/9 × (1 − 0.36 − 0.16) = 0.5333…
      位点2 A=5,C=5 (n=10): 10/9 × 0.5               = 0.5555…
      均值 = 0.5444…
    旧实现 (1−p_max 均值) = (0.4 + 0.5)/2 = 0.45 —— 必须不等于它。
    """
    counters = [Counter({"A": 6, "C": 4}), Counter({"A": 5, "C": 5})]
    got = _site_pi_mean(counters)
    assert got == pytest.approx((10 / 9 * 0.48 + 10 / 9 * 0.5) / 2, abs=1e-9)
    assert got != pytest.approx(0.45, abs=1e-6)


def test_site_pi_ignores_n_and_low_depth():
    """N 读段按缺失处理; 有效深度 <2 的位点记 NaN 并被跳过。"""
    counters = [Counter({"A": 6, "C": 4}),
                Counter({"N": 5}),            # 只有 N → 无有效等位基因, 剔除
                Counter({"A": 1, "T": 0}),    # 有效深度 1 → NaN, 剔除
                Counter({"A": 8, "G": 2})]    # 10/9 × (1 − 0.64 − 0.04) = 0.3555…
    got = _site_pi_mean(counters)
    assert got == pytest.approx((10 / 9 * 0.48 + 10 / 9 * 0.32) / 2, abs=1e-9)


def test_site_pi_all_uncomputable_is_nan():
    assert _site_pi_mean([Counter({"N": 3})]) != _site_pi_mean([Counter({"N": 3})])  # NaN
    assert _site_pi_mean([]) != _site_pi_mean([])  # NaN
