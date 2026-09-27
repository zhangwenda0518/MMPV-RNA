#!/usr/bin/env python3
"""
ess.py — Effective Sample Size (ESS) 计算（Tracer 同算法）

实现 Geyer (1992) initial POSITIVE sequence criterion, 与 Tracer (beast-mcmc
TraceCorrelation.analyseCorrelationNumeric) 严格一致。
金标准对账 2026-08-27: trace.jar 官方库 vs 本实现, 同 log 同 burnin 逐参数核对通过。

算法 (对应 TraceCorrelation.java 逐行):
    gammaStat[lag] = 有偏自协方差 (除以 samples-lag)
    varStat = gammaStat[0] + 2 * sum(成对和 gammaStat[2k-1]+gammaStat[2k], 仅当 >0)
    成对和 <= 0 即停止 (initial positive sequence, Tracer 无单调性检查)
    ACT = stepSize * varStat / gammaStat[0];  ESS = stepSize * samples / ACT
    maxLag = min(samples-1, 2000)   # Tracer MAX_LAG=2000

历史注记:
    v1 逐 lag 遇非正 acf 即截断 → 严重低估 (GCVA 50M run posterior 49.1 被报 6.9)
    v2 (2026-08-27 早) 加了 initial MONOTONE 检查 → 对 acf 零点抖动序列偏高
       (PSTVd posterior: 85 vs Tracer 14.6)
    v3 (当前) initial POSITIVE, 与 Tracer 完全一致; monotone=True 可复现 v2

⚠️ 与 YR-MPE mcmc_utils.calculate_ESS 的有意差异 (2026-09-16 上游核对)
------------------------------------------------------------------------
上游实现 (YR-MPE/YR_MPE/plugins/components/methods/mcmc_utils.py:59-100)
把成对和循环写成 `k = 0 ... P_k = rho_{2k} + rho_{2k+1}`, 于是第一项是
`P_0 = rho_0 + rho_1`。但 `rho_0 ≡ 1` 已经由 `tau_hat = -1 + 2*sum(P_k)` 的
那个 `-1` 代表, 因此 rho_0 被计入两次。

正确配对是 Tracer 的 `for (i = 1; i < maxLag; i += 2)`, 即
`(rho_1, rho_2), (rho_3, rho_4), ...` → `tau = 1 + 2*sum(P_k)`, k 从 1 起。
本模块 (v3) 即此口径。

实测对 Tracer 的中位相对偏差 (20 轮 x 6 档随机序列):
    样本量  我方     YR-MPE
      50    0.090%    1.758%   (上游最坏 26.75%)
     100    0.030%    0.497%   (上游最坏 18.96%)
     250    0.017%    1.503%   (上游最坏  9.86%)
     500    0.000%    1.532%   (上游最坏 18.47%)
    1000    0.000%    0.221%   (上游最坏  7.12%)
    3000    0.001%    1.458%   (上游最坏 35.75%)
6/6 档本实现都更贴近 Tracer。

次级差异: 上游 max_lag = min(n//2, 1000) 且 var 用 ddof=1;
本实现 max_lag = min(n-1, 2000) (Tracer MAX_LAG=2000) 且 ddof=0 (有偏 gamma_0)。

**不要为了"对齐上游"把 k 改成从 0 起** —— 那会引入上游的双计错误。
若历史结论用上游口径复核过 ESS 阈值, 需重跑。
证据脚本: MMPV-RNA/_consistency_check_20260916/ess_stats_robust.py
"""

import numpy as np


def calculate_ess(samples, monotone=False):
    """Effective Sample Size (Geyer initial positive sequence, Tracer 同算法).

    Parameters
    ----------
    samples : array-like, MCMC 样本
    monotone : bool, 额外启用 initial MONOTONE 检查 (成对和递增也截断);
        与 Tracer 不一致, 仅作敏感性分析用。

    Returns
    -------
    float : ESS
    """
    if len(samples) < 4:
        return len(samples)

    samples = np.asarray(samples, dtype=float)
    n = len(samples)

    mean_val = np.mean(samples)
    centered = samples - mean_val
    var = np.var(centered)  # 有偏 (ddof=0), 同 Tracer gammaStat[0]

    if var == 0:
        return n

    max_lag = min(n - 1, 2000)  # Tracer MAX_LAG = 2000
    acf = np.zeros(max_lag + 1)
    acf[0] = 1.0
    for lag in range(1, max_lag + 1):
        if lag >= n:
            break
        numerator = np.sum(centered[:n - lag] * centered[lag:])
        denominator = var * (n - lag)
        acf[lag] = numerator / denominator

    # Geyer initial POSITIVE sequence (与 Tracer TraceCorrelation 逐行一致):
    # varStat = gamma0 + 2*sum_{k>=1}(gamma_{2k-1}+gamma_{2k}), 成对和<=0 即停,
    # 无单调性检查 (Tracer 无此步)。tau = varStat/gamma0 = 1 + 2*sum(成对 rho 和)
    sum_pairs = 0.0
    prev_pair = np.inf
    k = 1
    while 2 * k <= max_lag:
        pair = acf[2 * k - 1] + acf[2 * k]
        if pair <= 0:
            break
        if monotone and pair > prev_pair + 1e-10:
            break
        sum_pairs += pair
        prev_pair = pair
        k += 1

    tau_hat = 1.0 + 2.0 * sum_pairs
    if tau_hat <= 0:
        return float(n)
    ess = n / tau_hat
    return max(1.0, min(ess, float(n)))


def calculate_95hpd(samples):
    """95% Highest Posterior Density interval (最短 95% HPD 区间).

    与 YR-MPE mcmc_utils.calculate_95HPD / Tracer 同算法:
    排序后找含 ceil(0.95*n) 个样本的最短连续区间。
    """
    samples = np.asarray(samples, dtype=float)
    n = len(samples)
    if n == 0:
        return (float("nan"), float("nan"))
    if n == 1:
        return (float(samples[0]), float(samples[0]))
    sorted_s = np.sort(samples)
    hpd_size = int(np.ceil(0.95 * n))
    if hpd_size >= n:
        return (float(sorted_s[0]), float(sorted_s[-1]))
    widths = sorted_s[hpd_size - 1:] - sorted_s[:n - hpd_size + 1]
    idx = int(np.argmin(widths))
    return (float(sorted_s[idx]), float(sorted_s[idx + hpd_size - 1]))
