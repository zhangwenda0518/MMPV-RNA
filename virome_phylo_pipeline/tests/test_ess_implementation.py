# -*- coding: utf-8 -*-
"""⑥ ESS 实现验证 + 收敛诊断出图的回归测试 (2026-10-14)。

原则: 本文件**只断言「金标准性质」与「与已验证数值逐档吻合」**, 不重新实现
任何 ESS / R-hat 公式 —— 公式一律调用 `utils.ess_convergence_fig` 里那几份
(它们本身是「独立第二实现」)。期望值全部来自实测(见 `FIGURES_20261014.md`
与 `ESS_VERIFICATION_20261014.md`), 不是估计值。

覆盖四条历史缺陷, 每条都对应一个真实的错:
  1. **理论 ESS 用错长度** —— AR(1) 序列在 N=27003 上生成、丢弃 WARMUP=500
     暂态, 实际参与估计的长度是 N_EFF=26503。按 N 算 ⇒ 理论值系统性高
     1.8866% ⇒ 5 档偏差**全部**偏负 1.89 个百分点(ρ=0.95 由 +5.1% 变 +3.2%)。
     指纹: 「新旧之差」**恒为同一个数**, 与 rho 无关 —— 故可精确断言。
  2. **只看 i.i.d. 均值会骗人** —— n=100 的 5 次重复里有一次只给 34.2(−66%),
     且**每档都有 2/5 触到 `min(ess, n)` 上限被封顶**。
  3. **合成序列必须逐字复刻 RNG 消耗顺序** —— 第一版忘了复刻 A 层消耗 ⇒
     ρ=0.9 本实现 ESS 算成 1085–1389, 而报告是 1410.9 ⇒ 图与报告互相矛盾。
  4. **`ucld.mean` 与 `ucld.mean.1` 数值逐位相同** ⇒ 不去重会让两张图一个
     16 行、一个 15 行。去重必须**按数值判定**, 不能按名字猜。

⚠ 依赖: numpy / matplotlib(必须)。R + coda(**可选**) —— 缺失时相关用例
`pytest.skip`, 绝不 fail: 出图不该因为环境没装 R 就阻断管线。

双端可跑: 本地 `C:\\Python312\\python.exe -m pytest tests/test_ess_implementation.py`
     服务器 `python -m pytest tests/test_ess_implementation.py`
"""
import ast
import io
import shutil
import sys
from pathlib import Path

import pytest

_HERE = Path(__file__).resolve().parent
_PHYLO = _HERE.parent
_REPO = _PHYLO.parent
for p in (str(_PHYLO), str(_REPO)):
    if p not in sys.path:
        sys.path.insert(0, p)

# ⑥ 的 utils/ 与其他管线 (③④⑤) 的同名包冲突: 全量 pytest 时若其他测试先导入了
# 别的 utils, 这里必须清缓存并把 ⑥ 目录插到最前。
for _m in [m for m in list(sys.modules) if m == "utils" or m.startswith("utils.")]:
    del sys.modules[_m]
sys.path.insert(0, str(_PHYLO))

import numpy as np  # noqa: E402

from utils.ess_convergence_fig import (  # noqa: E402
    AR1_N, AR1_RHOS, AR1_WARMUP, IID_NS, N_REP, SEED,
    calculate_ess, ess_batch_means, ess_ips_independent, generate_gold_standard,
    key_columns, read_log, rhat_python, rhat_via_coda)

_MODULE_PATH = _PHYLO / "utils" / "ess_convergence_fig.py"

# ── 已双端复核的期望值（改动合成规模就会失效 ⇒ 必须同步更新这两个常量） ──
N_EFF = 26503                                   # = AR1_N - AR1_WARMUP
AR1_DEV_PCT = {0.0: -2.93, 0.5: -6.54, 0.9: 1.15, 0.95: 5.12, 0.99: -0.42}
IID_DEV_PCT = {100: -20.56, 1000: -6.78, 10000: -2.91, 27003: -0.50}
# 误用 N 代替 N_EFF 时理论值抬高的比例 —— 5 档**共用同一个数**
LENGTH_FINGERPRINT_PCT = 100.0 * (AR1_N / N_EFF - 1.0)      # = 1.8866…


# ══════════════════════════════════════════════════════════════════
# 夹具
# ══════════════════════════════════════════════════════════════════
@pytest.fixture(scope="module")
def gold():
    """金标准三层数据只生成一次（AR(1) 是 Python 循环，别每例重跑）。"""
    iid, ar1, n_eff = generate_gold_standard()
    return {"iid": iid, "ar1": ar1, "n_eff": n_eff}


def _coda_available(rscript="Rscript"):
    if shutil.which(rscript) is None:
        return False
    import subprocess
    try:
        r = subprocess.run([rscript, "-e", 'library(coda); cat("ok")'],
                           capture_output=True, text=True, timeout=120)
    except Exception:
        return False
    return r.returncode == 0 and "ok" in (r.stdout or "")


CODA = _coda_available()
needs_coda = pytest.mark.skipif(not CODA, reason="Rscript + coda 不可用")


# ══════════════════════════════════════════════════════════════════
# ① 理论 ESS 必须用「实际参与估计的长度」
# ══════════════════════════════════════════════════════════════════
def test_gold_standard_length_is_post_warmup(gold):
    """P0: 理论值分母必须是 N_EFF=26503, 且每条序列实际长度就是这个数。

    这是 2026-10-14 查出的**硬错误**: 原版 `ess_t = N / act_t` 用 N=27003,
    而序列已 `x = x[500:]` ⇒ 理论值系统性偏高。
    """
    assert AR1_N - AR1_WARMUP == N_EFF
    assert gold["n_eff"] == N_EFF

    for rho in AR1_RHOS:
        d = gold["ar1"][rho]
        for i, x in enumerate(d["reps"]):
            assert len(x) == N_EFF, "rho=%.2f rep%d 长度 %d" % (rho, i, len(x))
        # 理论值 = 实际长度 / 理论 ACT（不是 N / ACT）
        assert d["ess_theory"] == pytest.approx(N_EFF / d["act_theory"], rel=1e-12)
        assert d["ess_theory"] != pytest.approx(AR1_N / d["act_theory"], rel=1e-12)


def test_theory_length_error_inflates_theory_by_constant_factor(gold):
    """指纹: 用错长度 ⇒ 理论值**恒高 1.8866%**（比值恒定, 与 rho 无关）。

    ⚠⚠ 这里必须说清**哪个量恒定**：恒定的是**理论值的比值** `N/N_EFF`,
    不是「偏差平移量」。平移量 = `(1+dev)·(1−N_EFF/N)`, 随 dev 微变 ⇒
    实测区间 **1.73–1.96 pp**（均值 ≈1.85）, **不是**严格恒定的 1.85。
    这条守卫就是防「把近似当恒等式」写进文档。
    """
    ratios, shifts = [], []
    for rho in AR1_RHOS:
        d = gold["ar1"][rho]
        ours = float(np.mean(d["ours"]))
        dev_right = 100.0 * (ours / (N_EFF / d["act_theory"]) - 1.0)
        dev_wrong = 100.0 * (ours / (AR1_N / d["act_theory"]) - 1.0)
        ratios.append((AR1_N / d["act_theory"]) / (N_EFF / d["act_theory"]))
        shifts.append(dev_right - dev_wrong)

    # 理论值的比值：严格恒定
    assert all(r == pytest.approx(AR1_N / N_EFF, rel=1e-12) for r in ratios)
    assert LENGTH_FINGERPRINT_PCT == pytest.approx(1.8866, abs=1e-4)

    # 偏差平移量：≈1.85 pp, 但不是常数（随 ρ 微变）
    assert float(np.mean(shifts)) == pytest.approx(1.85, abs=0.02)
    assert all(1.70 <= s <= 2.00 for s in shifts), "shifts=%s" % shifts
    assert max(shifts) - min(shifts) > 0.1, "平移量本应随 ρ 微变, 实测几乎恒定"

    # 方向守卫: 用错长度会让偏差**偏负**（原版就是这样把 +5.1 压成 +3.2 的）
    assert all(s > 0 for s in shifts)


def test_ar1_deviations_match_verified_values(gold):
    """5 档偏差必须与已验证数值逐档吻合（±0.05 个百分点）。

    ⚠ 期望值是**实测**的, 不是从报告抄的 —— 报告早期版本写过「<3.5%」,
    那是另一条(非同源)旧序列的结果, 与这里的同源数据矛盾。
    ⚠⚠ 三个「最大」在**不同** ρ 上, 混用就是静默错误:
        · 最大 |均值偏差|      = ρ=0.5（−6.54）
        · 最大 **带符号** 均值偏差 = ρ=0.95（+5.12）
        · 最大 |单次偏差|      = ρ=0.95（+16.44, rep1）
    """
    for rho, want in AR1_DEV_PCT.items():
        d = gold["ar1"][rho]
        got = 100.0 * (float(np.mean(d["ours"])) / d["ess_theory"] - 1.0)
        assert got == pytest.approx(want, abs=0.05), "rho=%.2f got %+.2f" % (rho, got)

    means = {r: v for r, v in AR1_DEV_PCT.items()}
    assert max(means, key=lambda r: abs(means[r])) == 0.5
    assert max(means, key=means.get) == 0.95

    per = {}
    for rho in AR1_RHOS:
        d = gold["ar1"][rho]
        per[rho] = max(abs(100.0 * (v / d["ess_theory"] - 1.0)) for v in d["ours"])
    assert max(per, key=per.get) == 0.95
    assert per[0.95] == pytest.approx(16.44, abs=0.05)

    # 结论强度: 5 档 |均值偏差| 全部 < 10%, 但**不是**"处处 <3.5%"
    assert max(abs(v) for v in AR1_DEV_PCT.values()) > 3.5
    assert max(abs(v) for v in AR1_DEV_PCT.values()) < 10.0


def test_rng_consumption_order_is_replicated(gold):
    """P0: 合成序列必须逐字复刻审查脚本的 RNG 消耗顺序。

    判据 = 已复核的锚点数值: ρ=0.9 本实现均值必须是 1410.9(报告同值)。
    第一版忘了复刻 A 层消耗 ⇒ 这里是 1085–1389 ⇒ 图与报告矛盾。
    """
    assert float(np.mean(gold["ar1"][0.9]["ours"])) == pytest.approx(1410.9, abs=0.5)
    assert float(np.mean(gold["ar1"][0.0]["ours"])) == pytest.approx(25726.0, abs=1.0)
    # 同 seed 必须完全可复现
    _, ar1b, _ = generate_gold_standard()
    assert np.array_equal(np.asarray(gold["ar1"][0.95]["ours"]),
                          np.asarray(ar1b[0.95]["ours"]))
    assert SEED == 20261014


# ══════════════════════════════════════════════════════════════════
# ② i.i.d. 层：均值收敛 + 触顶封顶 + 单次极差
# ══════════════════════════════════════════════════════════════════
def test_iid_mean_deviation_matches_verified_values(gold):
    for n, want in IID_DEV_PCT.items():
        got = 100.0 * (float(np.mean(gold["iid"][n])) / n - 1.0)
        assert got == pytest.approx(want, abs=0.05), "n=%d got %+.2f" % (n, got)


def test_iid_deviation_shrinks_with_n(gold):
    """白噪声下实现无偏 ⇒ 偏差随 n 收缩（负偏差幅度单调减小）。"""
    devs = [IID_DEV_PCT[n] for n in IID_NS]
    assert all(d < 0 for d in devs), "白噪声应略偏负(IPS 截尾), 实测 %s" % devs
    assert devs == sorted(devs), "|dev| 必须随 n 单调收缩, 实测 %s" % devs
    assert abs(devs[-1]) < 1.0 and abs(devs[0]) > 15.0


def test_iid_every_bucket_hits_the_min_cap(gold):
    """每档都有 2/5 次触到 `min(ess, n)` 上限被封顶。

    这是「只看均值会骗人」的直接证据: 封顶值让均值被**抬高**,
    同时掩盖了另外几次的大幅低估。
    """
    for n in IID_NS:
        v = np.asarray(gold["iid"][n], dtype=float)
        capped = int(np.sum(v >= n - 1e-9))
        assert capped == 2, "n=%d 触顶 %d 次(期望 2)" % (n, capped)
        assert float(v.max()) == pytest.approx(float(n), abs=1e-9)   # 封顶值 = n


def test_iid_single_rep_can_be_badly_wrong(gold):
    """n=100 有一次重复只给 34.2（−65.8%）—— 极差必须被记录, 不能被均值抹平。"""
    v = np.asarray(gold["iid"][100], dtype=float)
    worst = float(v.min())
    assert worst == pytest.approx(34.18, abs=0.05)
    assert 100.0 * (worst / 100.0 - 1.0) < -60.0
    # 同档最好的一次是封顶值 ⇒ 极差跨 3 倍
    assert float(v.max()) / worst > 2.9


# ══════════════════════════════════════════════════════════════════
# ③ 三种算法 / coda 对拍
# ══════════════════════════════════════════════════════════════════
def test_independent_ips_agrees_with_main_implementation(gold):
    """独立第二实现（逐 lag 显式循环）必须与本实现同量级（±2%）。"""
    for rho in AR1_RHOS:
        d = gold["ar1"][rho]
        a, b = float(np.mean(d["ours"])), float(np.mean(d["indep"]))
        assert a == pytest.approx(b, rel=0.02), "rho=%.2f ours=%.1f indep=%.1f" % (rho, a, b)


def test_batch_means_is_only_order_of_magnitude(gold):
    """batch-means 是第三种算法 ⇒ 只断言量级, **不断言位级一致**。"""
    for rho in AR1_RHOS:
        d = gold["ar1"][rho]
        bm = float(np.mean(d["bm"]))
        assert bm == pytest.approx(d["ess_theory"], rel=0.60), \
            "rho=%.2f bm=%.1f theory=%.1f" % (rho, bm, d["ess_theory"])


@needs_coda
def test_coda_crosscheck_matches_reported_panel_c(gold):
    """面板 C 的三个注记数字必须与实测吻合 —— 报告引用的就是这三个。

    报告/图注: median |Δ| = 6.0%、max = 17.9%、max 所在 ρ = 0.95。
    原版硬写 "(ρ=0.99)" ⇒ 改成 argmax 现算。本用例把「现算」钉住。
    ⚠ 容差放宽到 ±1.0 pp: coda 走 **AR 谱密度**、本实现走 **Geyer IPS**,
    属**两种算法**, 只能作量级旁证, 不是位级判据（位级判据见 panel D）。
    """
    from utils.ess_convergence_fig import coda_ess_on_series
    rhos = sorted(gold["ar1"])
    series = [x for r in rhos for x in gold["ar1"][r]["reps"]]
    coda = coda_ess_on_series(series)
    if coda is None:
        pytest.skip("coda 调用失败（环境问题, 不是实现问题）")
    ours = np.array([calculate_ess(x) for x in series])
    c = np.array(coda, dtype=float)
    assert len(c) == len(series)
    rel = np.abs(ours - c) / c * 100
    i_max = int(np.argmax(rel))
    assert float(np.median(rel)) == pytest.approx(6.0, abs=1.0)
    assert float(rel.max()) == pytest.approx(17.9, abs=1.0)
    assert rhos[i_max // N_REP] == 0.95, "max |Δ| 应在 ρ=0.95, 实得 ρ=%.2f" % rhos[i_max // N_REP]
    # 量级旁证: 25 条序列里没有一条差到 2 倍以上
    assert float(rel.max()) < 100.0


def test_ess_implementations_return_within_bounds():
    """所有 ESS 实现必须落在 [1, n] 内（`min(ess, n)` 是硬上界）。"""
    rng = np.random.default_rng(7)
    x = rng.normal(0, 1, 500)
    for f in (ess_ips_independent, ess_batch_means):
        v = f(x)
        assert 1.0 <= v <= len(x)
    # 常量序列: g0=0 ⇒ 退化守卫返回 n（不得 ZeroDivisionError）
    assert ess_ips_independent(np.zeros(300)) == pytest.approx(300.0)
    assert np.isnan(ess_batch_means(np.zeros(300)))


def test_systematic_note_judges_direction_by_sign():
    """方向注记必须**判符号** —— 混合符号不得被写成「一律偏低」。

    旧写法 `abs(rel).min()/max()` 在混合符号时会把「有高有低」静默说成
    「systematic lower」（与判据方向写反同类缺陷）。
    """
    from utils.ess_convergence_fig import systematic_note
    low = systematic_note(np.array([-3.5, -2.6, -1.0]))
    assert "lower" in low and "1.0–3.5" in low   # 幅度区间: |hi|–|lo| 升序
    high = systematic_note(np.array([1.0, 2.6, 3.5]))
    assert "higher" in high and "1.0–3.5" in high
    mixed = systematic_note(np.array([-3.5, 0.4, 2.6]))
    assert "mixed sign" in mixed
    assert "lower" not in mixed and "higher" not in mixed
    # 完全一致（全 0）不得被说成 lower/higher
    same = systematic_note(np.array([0.0, 0.0]))
    assert "no difference" in same
    assert "lower" not in same and "higher" not in same


# ══════════════════════════════════════════════════════════════════
# ④ 读 log / 关键参数去重
# ══════════════════════════════════════════════════════════════════
def _fake_log(tmp_path, n=200, dup=True, name="a.log", seed=3):
    cols = ["posterior", "treeModel.rootHeight", "ucld.mean",
            "ucld.mean" if dup else "skyline.popSize1"]
    p = tmp_path / name
    rng = np.random.default_rng(seed)
    with io.open(p, "w", encoding="utf-8") as f:
        f.write("# BEAST log\n")
        f.write("\t".join(cols) + "\n")
        for i in range(n):
            v = rng.normal(0, 1, len(cols))
            if dup:
                v[3] = v[2]                     # 重复列**数值逐位相同**
            f.write("\t".join("%.6f" % z for z in v) + "\n")
    return str(p)


def test_read_log_suffixes_duplicate_columns_and_applies_burnin(tmp_path):
    names, data = read_log(_fake_log(tmp_path, n=200), burnin_frac=0.10)
    assert names.count("ucld.mean") == 1
    assert "ucld.mean.1" in names                # 与 beast_parser 同口径
    assert len(data["posterior"]) == 180         # 200 × (1 − 0.10)
    assert data["ucld.mean"].shape == data["ucld.mean.1"].shape


def test_key_columns_drops_numerically_identical_columns(tmp_path):
    """P0: 16 → 15 的去重必须**按数值**判定（`ucld.mean.1` 数值与 `ucld.mean` 同）。

    不去重 ⇒ 两张图一个 16 行一个 15 行（历史缺陷）。
    """
    _, d = read_log(_fake_log(tmp_path, n=200, dup=True, name="a.log"))
    _, e = read_log(_fake_log(tmp_path, n=200, dup=True, name="b.log"))
    got = key_columns([d, e])
    assert "ucld.mean" in got and "ucld.mean.1" not in got
    assert "treeModel.rootHeight" in got and "posterior" in got

    # 非重复列不得被误删
    _, f = read_log(_fake_log(tmp_path, n=200, dup=False, name="c.log"))
    _, g = read_log(_fake_log(tmp_path, n=200, dup=False, name="d.log"))
    assert "skyline.popSize1" in key_columns([f, g])

    # 非关键列（不在 KEY_PARAM_PATTERNS 里）不得混入
    h = dict(d)
    h["traceExtra"] = np.zeros(180)
    assert "traceExtra" not in key_columns([h, h])


def test_key_columns_requires_column_in_all_chains(tmp_path):
    _, d = read_log(_fake_log(tmp_path, n=200, name="a.log"))
    e = {k: v for k, v in d.items() if k != "posterior"}
    assert "posterior" not in key_columns([d, e])   # 不是「所有链都有」


# ══════════════════════════════════════════════════════════════════
# ⑤ R-hat
# ══════════════════════════════════════════════════════════════════
def _ar1_chains(m=3, n=800, rho=0.9, seed=11, mean_shift=0.0):
    rng = np.random.default_rng(seed)
    out = []
    for k in range(m):
        e = rng.normal(0, 1, n)
        x = np.empty(n)
        x[0] = e[0]
        for t in range(1, n):
            x[t] = rho * x[t - 1] + e[t]
        out.append(x + mean_shift * k)
    return out


def test_rhat_python_returns_none_on_degenerate_input():
    """退化输入必须返回 None（不得抛异常、不得返回 nan）。"""
    assert rhat_python([np.arange(10.0)]) is None                  # 只有 1 条链
    assert rhat_python([np.ones(5), np.ones(5)]) is None           # W=0（常量链）
    assert rhat_python([np.arange(1.0), np.arange(1.0)]) is None   # n<2
    # ⚠ 非退化: 两条**相同**的非恒定链 W>0 ⇒ R-hat = sqrt((n-1)/n) < 1，不是 None
    r = rhat_python([np.arange(5.0), np.arange(5.0)])
    assert r is not None and r == pytest.approx(np.sqrt(4.0 / 5.0), rel=1e-12)


def test_rhat_python_detects_separated_chains():
    """链均值被拉开 ⇒ R-hat 必须远超 1.1（判据方向: 强信号应当通过）。"""
    good = rhat_python(_ar1_chains(mean_shift=0.0))
    assert good == pytest.approx(1.0, abs=0.15)
    for shift, lo in ((3.0, 1.5), (10.0, 3.0)):
        bad = rhat_python(_ar1_chains(mean_shift=shift))
        assert bad > lo, "shift=%.1f 得 R-hat=%.3f" % (shift, bad)
    # 单调性: 拉开越多, R-hat 越大
    vals = [rhat_python(_ar1_chains(mean_shift=s)) for s in (0.0, 1.0, 3.0, 10.0)]
    assert vals == sorted(vals)


@needs_coda
def test_rhat_python_matches_coda_when_converged():
    """近收敛时两种估计量必须几乎一致（±0.02）—— 证明公式同族。

    ⚠⚠ 这里**只对近收敛链**断言小容差。`rhat_python` 是经典式
    `sqrt(((n−1)/n·W + B/n)/W)`，而 coda 的单变量分支是 **B/W 分解 +
    自由度校正**（`sqrt(df.adj·R2.random)`）⇒ **不是同一个估计量**，
    coda 恒偏大，且混合越差差得越多（实测 R̂≈3–15 时差 17–33%）。
    曾经把这条写成"与 coda 同口径"是错的，见 `test_rhat_python_under_reports_vs_coda`。
    """
    chains = _ar1_chains(m=3, n=800, rho=0.9, seed=5)
    keys = ["posterior"]
    res = rhat_via_coda([{"posterior": c} for c in chains], keys)
    if res is None:
        pytest.skip("coda 调用失败（环境问题, 不是实现问题）")
    got = rhat_python(chains)
    assert got == pytest.approx(res["posterior"]["rhat"], abs=0.02), \
        "python=%.4f coda=%.4f" % (got, res["posterior"]["rhat"])
    # 近收敛 ⇒ 相对差 <0.5%
    assert abs(got - res["posterior"]["rhat"]) / res["posterior"]["rhat"] < 0.005


@needs_coda
def test_rhat_python_under_reports_vs_coda():
    """方向守卫: coda 的校正估计量**恒 ≥** 经典式, 且差距随混合变差放大。

    这条是 2026-10-14 实测查出的实质问题: 回退到纯 Python 会**系统性低估**
    R̂ ⇒ 判据线 1.1 可能被静默放过（实测最差相对差 33%）。
    """
    rows = []
    for rho, shift in ((0.0, 0.0), (0.0, 1.0), (0.0, 3.0), (0.0, 10.0),
                       (0.9, 10.0)):
        chains = _ar1_chains(m=3, n=800, rho=rho, seed=17, mean_shift=shift)
        res = rhat_via_coda([{"posterior": c} for c in chains], ["posterior"])
        if res is None:
            pytest.skip("coda 调用失败（环境问题, 不是实现问题）")
        rows.append((rhat_python(chains), res["posterior"]["rhat"]))

    assert all(p <= c + 1e-9 for p, c in rows), "coda 应恒 ≥ 经典式: %s" % rows
    # 近收敛时几乎相同; 混合差时差距显著（>10% 相对）
    assert (rows[0][1] - rows[0][0]) / rows[0][1] < 0.005
    assert max((c - p) / c for p, c in rows) > 0.10


def test_rhat_via_coda_returns_none_without_rscript():
    """R 不可用时必须返回 None 而不是抛异常 —— 出图不该阻断管线。"""
    chains = _ar1_chains(m=3, n=200, seed=9)
    res = rhat_via_coda([{"posterior": c} for c in chains], ["posterior"],
                        rscript="definitely_not_a_real_rscript_xyz")
    assert res is None


# ══════════════════════════════════════════════════════════════════
# ⑥ 硬规则⑦: 标题/注记里禁硬写「数据结论」数字（AST 判定, 不用字符串匹配）
# ══════════════════════════════════════════════════════════════════
_PLOT_CALLS = {"set_title", "suptitle", "annotate", "text", "figtext",
               "set_xlabel", "set_ylabel", "set_label"}
# 「数字-破折号-数字-单位」= 硬写区间结论的形状（历史 bug: "8–55× loss"）
_HARDCODE_RE = r"\d+\s*[–—]\s*\d+\s*[×xX%]"


def _iter_literal_strings(tree):
    """产出 (行号, 字符串) —— 跳过 docstring，只取真正会被画出来的字面量。"""
    docstrings = set()
    for node in ast.walk(tree):
        if isinstance(node, (ast.Module, ast.FunctionDef, ast.AsyncFunctionDef,
                             ast.ClassDef)):
            ds = ast.get_docstring(node, clean=False)
            if ds:
                docstrings.add(ds)
    for node in ast.walk(tree):
        if not isinstance(node, ast.Constant) or not isinstance(node.value, str):
            continue
        if node.value in docstrings:
            continue
        yield node.lineno, node.value


def _plot_call_strings(tree):
    """AST 判定: 取绘图调用里作为实参出现的字符串字面量（含 f-string 片段）。"""
    out = []
    for node in ast.walk(tree):
        if not isinstance(node, ast.Call):
            continue
        fn = node.func
        name = fn.attr if isinstance(fn, ast.Attribute) else getattr(fn, "id", None)
        if name not in _PLOT_CALLS:
            continue
        for a in list(node.args) + [k.value for k in node.keywords]:
            for n2, s in _iter_literal_strings(a):
                out.append((n2, s))
    return out


def test_no_hardcoded_range_conclusions_in_plot_text():
    """P0: 图里被画出来的文字不得硬写「区间×倍率」这类结论数字。

    历史 bug: `_fig_convergence_diag.py` 曾硬写 `"8–55× loss"` 而实测是
    `8–101×`（现算 `1/max(ra)`、`1/min(ra)` 才对）。**阈值常量**（`1.1`、
    `200`、`1.01`）是判据不是结论, 允许保留; 本用例只拦「区间+单位」形状。
    """
    import re
    tree = ast.parse(_MODULE_PATH.read_text(encoding="utf-8"))
    bad = [(ln, s) for ln, s in _plot_call_strings(tree)
           if re.search(_HARDCODE_RE, s)]
    assert not bad, "绘图文字里硬写了区间结论: %s" % bad


def test_no_hardcoded_decimal_percent_in_plot_text():
    """同上, 拦「硬写小数百分比」（如 `"6.0%"`）—— 结论数字必须 `%` 占位现算。"""
    import re
    tree = ast.parse(_MODULE_PATH.read_text(encoding="utf-8"))
    bad = [(ln, s) for ln, s in _plot_call_strings(tree)
           if re.search(r"\d+\.\d+\s*%", s)]
    assert not bad, "绘图文字里硬写了百分比结论: %s" % bad


def test_module_docstring_numbers_match_measured_values(gold):
    """文档串里写的「已知数值」必须与实测吻合 —— 文档也是产物, 也会过期。

    判据 = 从源码文档串里抓出 ρ=0.95 的**改后**偏差, 与实测比对。
    查不到就 skip（文档措辞变了不等于实现错）。
    """
    import re
    src = _MODULE_PATH.read_text(encoding="utf-8")
    doc = ast.get_docstring(ast.parse(src)) or ""
    m = re.search(r"ρ=0\.95[^\n]{0,20}?→\s*\+?(-?\d+\.\d+)%", doc)
    if not m:
        pytest.skip("文档串未包含 ρ=0.95 的偏差改写记录")
    assert float(m.group(1)) == pytest.approx(AR1_DEV_PCT[0.95], abs=0.1)


# ══════════════════════════════════════════════════════════════════
# ⑦ 出图集成: 两张图都要能出, 且过版式体检
# ══════════════════════════════════════════════════════════════════
def _run_module_capture(argv):
    from utils.check_fig_layout import _run_and_capture
    return _run_and_capture(argv, module="utils.ess_convergence_fig")


def test_verify_only_emits_pdf_and_png(tmp_path):
    """`--verify-only` 必须落盘 pdf + png（HTML 内联靠 pdf，回归判据靠 png）。"""
    figs = _run_module_capture(["prog", "--verify-only", "--outdir", str(tmp_path)])
    assert len(figs) == 1, "应恰好产出 1 张图, 实得 %d" % len(figs)
    for ext in (".pdf", ".png"):
        p = tmp_path / ("ess_implementation_verification" + ext)
        assert p.exists() and p.stat().st_size > 2000, "缺 %s" % p


def test_verification_figure_passes_layout_check(tmp_path):
    """版式体检: 0 对文字重叠 + 无 TEXT-OFF-CANVAS。

    ⚠ 体检会**重渲**目标模块 —— 跑完以体检后产物为准。
    """
    from utils.check_fig_layout import check_figure
    figs = _run_module_capture(["prog", "--verify-only", "--outdir", str(tmp_path)])
    off, ov = check_figure(figs[0], tag="verify", verbose=False)
    assert (off, ov) == (0, 0), "OFF=%d 重叠=%d" % (off, ov)


def test_diagnostics_runs_and_passes_layout_check(tmp_path):
    """两条真实形状的链 log ⇒ 出图②, 且过版式体检。"""
    from utils.check_fig_layout import check_figure
    logs = [_fake_log(tmp_path, n=200, name="a.log", seed=1),
            _fake_log(tmp_path, n=200, name="b.log", seed=2)]
    figs = _run_module_capture(["prog", "--chains", ",".join(logs),
                                "--outdir", str(tmp_path)])
    assert len(figs) == 1
    for ext in (".pdf", ".png"):
        p = tmp_path / ("ess_convergence_diagnostics" + ext)
        assert p.exists() and p.stat().st_size > 2000, "缺 %s" % p
    off, ov = check_figure(figs[0], tag="diag", verbose=False)
    assert (off, ov) == (0, 0), "OFF=%d 重叠=%d" % (off, ov)


def test_cli_skips_diagnostics_with_one_chain(tmp_path, capsys):
    """只有 1 条链 ⇒ 必须**跳过并告警**（rc=2）, 不得抛异常、不得出空图。

    ⚠ R-hat 需要**链间方差**, 用拼接后的 merged.log 算不出 ⇒ 单链场景
    只能跳过, 这是 `merge_results.py` 里 `len(_single) >= 2` 守卫的对应面。
    """
    p = _fake_log(tmp_path, n=200, name="only.log")
    figs = _run_module_capture(["prog", "--chains", p, "--outdir", str(tmp_path)])
    assert figs == [], "单链不应产出任何图"
    assert not (tmp_path / "ess_convergence_diagnostics.pdf").exists()
    out = capsys.readouterr().out
    assert "跳过" in out and "≥2" in out


def test_fallback_path_without_coda_still_passes_layout(tmp_path):
    """无 R/coda 的**回退路径**也必须过版式体检。

    ⚠ 这条分支本地（有 coda）跑不到：回退时副标题变成
    「python fallback (under-reports R-hat)」且 C 面板换成说明框 ——
    文本变长/换形都可能撞版式, 必须单独体检。
    """
    from utils.check_fig_layout import check_figure
    logs = [_fake_log(tmp_path, n=200, name="a.log", seed=1),
            _fake_log(tmp_path, n=200, name="b.log", seed=2)]
    for tag, argv in (
            ("verify", ["prog", "--verify-only", "--rscript", "nope_xyz",
                        "--outdir", str(tmp_path)]),
            ("diag", ["prog", "--chains", ",".join(logs), "--rscript",
                      "nope_xyz", "--outdir", str(tmp_path)])):
        figs = _run_module_capture(argv)
        assert len(figs) == 1
        off, ov = check_figure(figs[0], tag=tag, verbose=False)
        assert (off, ov) == (0, 0), "%s: OFF=%d 重叠=%d" % (tag, off, ov)
