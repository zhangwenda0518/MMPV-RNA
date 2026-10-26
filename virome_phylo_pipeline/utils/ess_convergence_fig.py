#!/usr/bin/env python3
"""ess_convergence_fig.py — ESS 实现验证 + MCMC 收敛诊断出图 (2026-10-14)

由一次性审查脚本（`reports/_archive/review_dirs/_gcva_time_review_20261013/`）
正式并入管线。回答两个不同的问题：

  ┌ 问题                         ┌ 产物                                 ┌ 依赖
  │ ① ESS 算错了没有？          │ ess_implementation_verification.pdf  │ 无（自造合成序列）
  │ ② 这个 run 的结果能不能用？  │ ess_convergence_diagnostics.pdf      │ 本 run 的链 log
  └                             └                                      └

落盘到 `<merged>/` ⇒ `utils/report_builder.py` 的 `time` 模块 plots glob
含 `**/*.pdf` ⇒ **自动内联进 HTML 报告**，无需改 report_builder。

① 金标准验证（4 面板，通栏 9.32 in）
====================================
  A  i.i.d. 白噪声：理论 ESS = n。给「均值 + 5 次极差 + 触顶标记」。
     ⚠⚠ **只标均值会骗人**：实测 n=100 有一次重复只给 34.2（−66%），且每档
     5 次里有 2 次触到 `min(ess, n)` 上限被封顶。只看均值 20.6% 会误判成
     「实现有偏」，只看 n=2.7e4 的 0.5% 又会误判成「处处都准」。
  B  AR(1) 已知 rho：三方法（本实现 / 独立 IPS / batch-means）vs 理论。
     ⚠⚠ **理论值必须用「实际参与估计的长度」** —— 序列在 N=27003 上生成后
     丢弃 WARMUP=500 暂态，实际长度 N_EFF=26503。原版按 N 算 ⇒ 理论值恒高
     **1.8866%**（= N/N_EFF − 1，**比值恒定，与 ρ 无关**）⇒ 换算到偏差约
     **1.73–1.96 pp**（均值 ≈1.85 pp，非严格恒定：平移量 ∝ (1+dev)）。
     ρ=0.95 由 +3.2% → +5.1%。
     判据 = 第三方算法反推：coda 对白噪声**精确返回 n**（实测 26503）。
  C  与 coda 对拍：**另一种算法**（AR 谱密度 vs Geyer IPS）⇒ 只能算量级旁证，
     不是位级判据。⚠ R/coda 不可用时本面板画成「跳过」说明框（不静默消失）。
  D  真实数据：本实现 vs 独立重写 IPS，逐参数。

② 本 run 收敛诊断（3 面板）
===========================
  A  R-hat 森林图（Gelman-Rubin，**不经过任何 ESS 算法**）：判据线 1.01 / 1.1。
  B  ESS 条形图：逐链求和（coda 口径）vs 拼接 pooled（本项目口径）。
  C  pooled / Σ(per-chain) 比值：= 1 才叫混合良好。

⚠⚠ B 与 C 是**同一组 ESS 数据的两种画法**（C 由 B 导出），**不是两条独立证据**；
本图只有 **2 条**独立证据（R-hat；ESS 口径差）。这条已写进图的副标题，防止
下游把一张图当成两条佐证。

R-hat 来源
==========
优先用 `utils/ess_rhat.R`（coda::gelman.diag，与报告里的数值同源）；R/coda
不可用时退回纯 Python 经典 Gelman-Rubin 点估计（**不画上界**，因 coda 的
97.5% 上界公式是库内实现，不复制一份以免出现两个「上界」口径）。

⚠⚠⚠ **两者不是同一个估计量**（2026-10-14 实测定性，勿再写「同口径」）：
  · 本模块 = 经典式 `sqrt(((n−1)/n·W + B/n)/W)`；
  · coda 的单变量分支 = **B/W 分解 + 自由度校正**，即
    `sqrt(df.adj · R2.random)`（`gelman.diag` 源码 `psrf <- cbind(sqrt(df.adj*R2.estimate), …)`），
    属 Brooks–Gelman 校正版。
  ⇒ 实测（3 链 × 800 步）coda **恒 ≥** 本模块，差值随混合变差迅速放大：
      近收敛（R̂≈1.00）差 ≤0.06% 相对；R̂≈1.1 差约 2%；R̂≈3–15 差 **17–33%**。
  ⇒ **回退到纯 Python 会系统性低估 R̂** ⇒ 判据线 1.1 可能被静默放过。
    故 R 不可用时副标题会显式标注「python fallback (under-reports R-hat)」，
    且下游**不得**把回退值当作与 coda 等价的判据。要判据就装 R + coda。

判据纪律（硬规则 ⑦）
====================
标题/注记里的**一切数字都现算**（`%.2f` / 计数 / `.min()` / `.max()` / argmax），
禁硬写。历史教训：曾硬写「8–55×」而实测为「8–101×」，硬写「(ρ=0.99)」而
实测最大偏差在 ρ=0.95（⚠ 指 **C 面板的 coda 对拍** `|Δ|/coda`，不是理论偏差
——理论偏差的最大值另在 ρ=0.5，两个口径不可混）。**常量阈值**（`MIN_ESS=200`、
`R-hat=1.1`）可以留，因为它们是判据不是结论。

用法
====
  # ① 只用本 run 的链 log 出收敛诊断
  python -m utils.ess_convergence_fig --chains a.log,b.log,c.log --outdir merged

  # ② 金标准验证（可另附真实数据做 panel D）
  python -m utils.ess_convergence_fig --verify --outdir merged \\
      --chains a.log,b.log,c.log

  # ③ 只出验证图（不接真实数据）
  python -m utils.ess_convergence_fig --verify-only --outdir /tmp/x
"""
from __future__ import annotations

# ── 图件四格式导出 (SVG / PNG / PDF / PNG@1200dpi) — 默认开启 ──────────
# 2026-10-20: 所有出图脚本统一口径; 用 MMPV_FIG_EXPORT=0 可整体关闭,
# MMPV_FIG_FORMATS=svg,png,pdf,hires 可裁剪, MMPV_FIG_HIRES_DPI 可改高清 dpi。
try:
    from mmpv_common.figure_export import install_savefig_hook as _install_fig_hook
except ImportError:  # 仓库根不在 sys.path (脚本被单独执行) → 向上找 mmpv_common/
    import os as _os
    import sys as _sys
    _d = _os.path.dirname(_os.path.abspath(__file__))
    while _d and not _os.path.isdir(_os.path.join(_d, "mmpv_common")):
        _parent = _os.path.dirname(_d)
        if _parent == _d:
            break
        _d = _parent
    if _d and _d not in _sys.path:
        _sys.path.insert(0, _d)
    from mmpv_common.figure_export import install_savefig_hook as _install_fig_hook
_install_fig_hook()
# ──────────────────────────────────────────────────────────────────────

import argparse
import io
import os
import shutil
import subprocess
import sys
import tempfile

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ⚠ 双路径 import：既能 `import utils.ess_convergence_fig`（包内），
# 也能 `python utils/ess_convergence_fig.py`（此时 sys.path[0]=utils/）。
try:
    from utils.ess import calculate_ess
    from utils.sci_palette import (
        apply_sci_style, OKABE_ITO, COL_W_1, despine, sci_legend, SEQ_CMAP)
except ImportError:                                    # pragma: no cover
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    from utils.ess import calculate_ess
    from utils.sci_palette import (
        apply_sci_style, OKABE_ITO, COL_W_1, despine, sci_legend, SEQ_CMAP)

# 关键参数集的**唯一来源** —— 与收敛门禁同源，避免两处各写一份而慢慢漂移。
# ⚠ 门禁那边有实测教训：漏 `coefficientOfVariation` 会让最差 ESS 被低估。
try:
    from utils.convergence_gate import _is_key_param
except ImportError:                                    # pragma: no cover
    from convergence_gate import _is_key_param

_CM = plt.get_cmap(SEQ_CMAP)          # 数值型（有序）变量统一走 SEQ_CMAP

C_BAD = OKABE_ITO[5]                  # vermillion — 未达标 / 本实现
C_MID = OKABE_ITO[0]                  # orange     — 临界 / coda
C_OK = OKABE_ITO[2]                   # green      — 达标 / 独立实现
C_BM = OKABE_ITO[4]                   # blue       — batch-means
C_REF = "#444444"

MIN_ESS = 200.0                       # 与 utils/sample_floor.MIN_ESS 同值（判据常量）
RHAT_PRACTICAL = 1.1                  # 操作判据（严格判据 1.01 只写进注记）
RHAT_STRICT = 1.01

# 合成数据的规模（与审查脚本逐字一致 —— 改了就无法与已验证的数值对拍）
SEED = 20261014
IID_NS = (100, 1000, 10000, 27003)
AR1_RHOS = (0.0, 0.5, 0.9, 0.95, 0.99)
AR1_N = 27003
AR1_WARMUP = 500
N_REP = 5


# ══════════════════════════════════════════════════════════════════
# 读 log / 参数筛选
# ══════════════════════════════════════════════════════════════════
def read_log(path, burnin_frac=0.10):
    """读 BEAST tab log，返回 (列名列表, {列名: np.ndarray})，已去 burn-in。

    重名列按出现顺序加 `.1`/`.2` 后缀 —— 与 `utils/beast_parser.py`（pandas）
    同口径。⚠ 真实 BEAST log 里 `ucld.mean` 会出现两次且数值 100% 相同。
    """
    cols, data = None, None
    with io.open(path, encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if cols is None:
                cols, data = p, [[] for _ in p]
                continue
            if len(p) != len(cols):
                continue
            for i, v in enumerate(p):
                data[i].append(v)

    seen, names = {}, []
    for c in cols:
        if c in seen:
            seen[c] += 1
            names.append("%s.%d" % (c, seen[c]))
        else:
            seen[c] = 0
            names.append(c)

    out = {}
    for i, k in enumerate(names):
        v = data[i]
        v = v[int(len(v) * burnin_frac):]
        try:
            out[k] = np.asarray([float(x) for x in v], dtype=float)
        except ValueError:
            continue                      # 非数值列（如 tree 行）跳过
    return names, out


def key_columns(chains):
    """取「所有链都有的关键参数」，并剔除**数值完全重复**的列。

    ⚠ 去重不是可选优化：BEAST 的 `ucld.mean` 与 `ucld.mean.1` 数值逐位相同，
    留着会让人以为有两个参数，且与其它面板的参数数对不上（历史缺陷：
    两张图一个 16 行、一个 15 行）。这里**按数值判定**而不是按名字猜，
    所以任何重复列都会被自动清掉。
    """
    common = set(chains[0])
    for c in chains[1:]:
        common &= set(c)
    keys = [k for k in chains[0] if k in common and _is_key_param(k)]

    kept, seen_vals = [], []
    for k in keys:
        v = chains[0][k]
        if any(v.shape == s.shape and np.array_equal(v, s) for s in seen_vals):
            continue                      # 与已保留列逐位相同 ⇒ 重复列
        kept.append(k)
        seen_vals.append(v)
    return kept


def truncate_common(chains, keys):
    """各链按最短长度等长截断（合并口径要求等长）。"""
    n = min(len(c[keys[0]]) for c in chains)
    return {k: [c[k][:n] for c in chains] for k in keys}, n


# ══════════════════════════════════════════════════════════════════
# ESS 的独立第二实现 / 第三实现（交叉验证用，与 utils.ess 代码路径不同）
# ══════════════════════════════════════════════════════════════════
def ess_ips_independent(x):
    """独立重写的 Geyer IPS（逐 lag 显式循环，不共用 utils.ess 的任何代码）。"""
    x = np.asarray(x, dtype=float)
    n = len(x)
    xc = x - x.mean()
    g0 = np.dot(xc, xc) / n
    if g0 == 0:
        return float(n)
    L = min(n - 1, 2000)
    gam = np.empty(L + 1)
    gam[0] = g0
    for lag in range(1, L + 1):
        gam[lag] = np.dot(xc[:n - lag], xc[lag:]) / n
    s = 0.0
    for k in range(1, L // 2 + 1):
        pair = gam[2 * k - 1] + gam[2 * k]
        if pair <= 0:
            break
        s += pair
    tau = (g0 + 2 * s) / g0
    return max(1.0, min(n / tau, float(n)))


def ess_batch_means(x, nb=50):
    """Batch-means ESS（第三种算法，用于给出「算法间离散度」的量级）。"""
    x = np.asarray(x, dtype=float)
    n = len(x)
    b = n // nb
    if b < 2:
        return float("nan")
    batches = x[:nb * b].reshape(nb, b).mean(axis=1)
    var_b = batches.var(ddof=1)
    var_x = x.var(ddof=1)
    if var_b <= 0 or var_x <= 0:
        return float("nan")
    act = var_b / (var_x / b)
    return max(1.0, min(n / act, float(n)))


# ══════════════════════════════════════════════════════════════════
# R-hat
# ══════════════════════════════════════════════════════════════════
def rhat_python(chains):
    """经典 Gelman-Rubin R-hat（点估计）。**不是 coda 的估计量**（见模块头）。

    R-hat = sqrt( ((n-1)/n * W + B/n) / W )
      B = n * Var(链均值)   —— 链间方差（乘 n）
      W = mean(链内方差)

    ⚠⚠⚠ coda 的单变量分支用 `sqrt(df.adj·R2.random)`（B/W 分解 + 自由度
    校正）⇒ 本函数**恒 ≤** coda 值，且混合越差差得越多（R̂≈1 差 <0.1%，
    R̂≈3–15 差 17–33%）。所以它只是 **coda 不可用时的降级路径**，
    不能当作与 coda 等价的判据。
    ⚠ 只返回点估计：coda 的 97.5% 上界是库内公式，这里不复制一份，
    避免出现两个互相打架的「上界」。上界只在上界来源为 coda 时绘制。
    """
    m = len(chains)
    n = min(len(c) for c in chains)
    if m < 2 or n < 2:
        return None
    x = np.vstack([np.asarray(c[:n], dtype=float) for c in chains])
    xbar = x.mean(axis=1)
    B = n * xbar.var(ddof=1)
    W = x.var(axis=1, ddof=1).mean()
    if not np.isfinite(W) or W <= 0:
        return None
    vh = (n - 1) / n * W + B / n
    return float(np.sqrt(vh / W))


def rhat_via_coda(chain_arrays, keys, rscript="Rscript", timeout=600):
    """调 `utils/ess_rhat.R`（coda::gelman.diag）拿 R-hat + 上界 + coda ESS。

    ⚠ 拿到的是 coda 的 **B/W 分解 + 自由度校正** 估计量（Brooks–Gelman），
    与 `rhat_python` 的经典式**不等价**（coda 恒偏大，混合越差差越多）。

    返回 {param: {"rhat":…, "upper":…, "coda_ess":…}}；任一步失败返回 None
    （调用方退回 `rhat_python`）。⚠ 失败**不抛异常** —— 出图不该阻断管线。
    """
    rfile = os.path.join(os.path.dirname(os.path.abspath(__file__)), "ess_rhat.R")
    if not os.path.exists(rfile) or shutil.which(rscript) is None:
        return None
    tmp = tempfile.mkdtemp(prefix="essrhat_")
    try:
        for idx, arrs in enumerate(chain_arrays, 1):
            with io.open(os.path.join(tmp, "chain%d.csv" % idx), "w",
                         encoding="utf-8") as f:
                f.write(",".join(keys) + "\n")
                n = len(arrs[keys[0]])
                for i in range(n):
                    f.write(",".join("%.10g" % float(arrs[k][i]) for k in keys) + "\n")
        out_tsv = os.path.join(tmp, "rhat.tsv")
        r = subprocess.run([rscript, rfile, tmp, out_tsv],
                           capture_output=True, text=True,
                           encoding="utf-8", errors="replace", timeout=timeout)
        if r.returncode != 0 or not os.path.exists(out_tsv):
            return None
        res = {}
        with io.open(out_tsv, encoding="utf-8") as f:
            hdr = f.readline().rstrip("\n").split("\t")
            for line in f:
                if not line.strip():
                    continue
                d = dict(zip(hdr, line.rstrip("\n").split("\t")))
                res[d["param"]] = {
                    "rhat": float(d["rhat"]),
                    "upper": float(d["upper"]),
                    "coda_ess": float(d["coda_ess"]) if d.get("coda_ess") else None,
                }
        return res or None
    except Exception:
        return None
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


# ══════════════════════════════════════════════════════════════════
# 通用小件
# ══════════════════════════════════════════════════════════════════
def _note_box(ax, txt, xy=(0.97, 0.03), ha="right", va="bottom", fontsize=5.4):
    ax.annotate(txt, xy, xycoords="axes fraction", fontsize=fontsize,
                color="#444444", ha=ha, va=va,
                bbox=dict(boxstyle="round,pad=0.22", fc="white",
                          ec="#cccccc", lw=0.4, alpha=0.92))


def _save(fig, out_pdf, also_png=True):
    os.makedirs(os.path.dirname(os.path.abspath(out_pdf)) or ".", exist_ok=True)
    paths = [out_pdf]
    fig.savefig(out_pdf, dpi=300, bbox_inches="tight")
    if also_png and out_pdf.lower().endswith(".pdf"):
        png = out_pdf[:-4] + ".png"
        fig.savefig(png, dpi=300, bbox_inches="tight")
        paths.append(png)
    plt.close(fig)
    for p in paths:
        print("wrote", p)
    return paths


# ══════════════════════════════════════════════════════════════════
# 图①  ESS 实现验证（金标准）
# ══════════════════════════════════════════════════════════════════
def generate_gold_standard():
    """生成金标准三层数据。**逐字复刻**审查脚本的 RNG 消耗顺序。

    ⚠⚠ 顺序敏感：先 A 层 4 档 × 5 次 `rng.normal`（20 次），再 B 层 5 档 × 5 次。
    任何一处增删都会让序列与已验证的数值对不上（历史缺陷：第一版忘了复刻
    A 层消耗 ⇒ ρ=0.9 本实现 ESS 算成 1085–1389 而报告是 1410.9 ⇒ 图与报告矛盾）。
    """
    rng = np.random.default_rng(SEED)

    iid = {n: [calculate_ess(rng.normal(0, 1, n)) for _ in range(N_REP)]
           for n in IID_NS}

    n_eff = AR1_N - AR1_WARMUP
    ar1 = {}
    for rho in AR1_RHOS:
        reps = []
        for _ in range(N_REP):
            e = rng.normal(0, 1, AR1_N + AR1_WARMUP)
            x = np.empty(AR1_N)
            x[0] = e[0]
            for t in range(1, AR1_N):
                x[t] = rho * x[t - 1] + e[t]
            x = x[AR1_WARMUP:]
            assert len(x) == n_eff, len(x)       # 防回归：长度必须 = N_EFF
            reps.append(x)
        act_t = (1 + rho) / (1 - rho)
        ar1[rho] = {"reps": reps, "act_theory": act_t, "ess_theory": n_eff / act_t,
                    "ours": [calculate_ess(x) for x in reps],
                    "indep": [ess_ips_independent(x) for x in reps],
                    "bm": [ess_batch_means(x) for x in reps]}
    return iid, ar1, n_eff


def coda_ess_on_series(series, rscript="Rscript", timeout=900):
    """对一批序列跑 coda::effectiveSize，返回同长度列表；不可用则 None。"""
    rfile = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                         "ess_coda_synth.R")
    if not os.path.exists(rfile) or shutil.which(rscript) is None:
        return None
    tmp = tempfile.mkdtemp(prefix="esscoda_")
    try:
        for i, x in enumerate(series):
            with io.open(os.path.join(tmp, "s%03d.csv" % i), "w",
                         encoding="utf-8") as f:
                f.write("x\n")
                for v in x:
                    f.write("%.10g\n" % v)
        out_tsv = os.path.join(tmp, "coda.tsv")
        r = subprocess.run([rscript, rfile, tmp, out_tsv],
                           capture_output=True, text=True,
                           encoding="utf-8", errors="replace", timeout=timeout)
        if r.returncode != 0 or not os.path.exists(out_tsv):
            return None
        vals = []
        with io.open(out_tsv, encoding="utf-8") as f:
            hdr = f.readline().rstrip("\n").split("\t")
            for line in f:
                if line.strip():
                    vals.append(float(dict(zip(hdr, line.rstrip("\n").split("\t")))["ess_coda"]))
        return vals if len(vals) == len(series) else None
    except Exception:
        return None
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def _p_iid(ax, iid):
    ns = sorted(iid)
    ax.plot(ns, ns, "-", color=C_REF, lw=1.0, zorder=1, label="theory  ESS = n")
    worst, devs = [], []
    for i, n in enumerate(ns):
        v = np.array(iid[n], dtype=float)
        jit = (np.arange(len(v)) - (len(v) - 1) / 2) * 0.055
        capped = np.isclose(v, float(n), rtol=1e-9)
        ax.vlines(n, v.min(), v.max(), color=C_BAD, lw=0.7, alpha=0.5, zorder=2)
        ax.scatter(np.full(len(v), n) * (1 + jit), v, s=16, facecolor="none",
                   edgecolor=C_BAD, lw=0.85, zorder=3,
                   label="impl. (5 reps)" if i == 0 else None)
        ax.scatter(np.full(int(capped.sum()), n) * (1 + jit[capped]),
                   v[capped], s=15, color=C_BAD, lw=0, zorder=3.5,
                   label="capped at n" if i == 0 else None)
        ax.scatter([n], [v.mean()], s=30, marker="D", color=C_BAD,
                   edgecolor="white", lw=0.5, zorder=4,
                   label="mean" if i == 0 else None)
        worst.append((v.min() - n) / n * 100)
        devs.append((v.mean() - n) / n * 100)
    ax.set_xscale("log")
    ax.set_yscale("log")              # ⚠ 线性 y 会把 n=100 的点压成 0
    ax.set_ylim(bottom=22)
    ax.set_xlabel("sample size  n")
    ax.set_ylabel("ESS")
    n_lo, n_hi = ns[0], ns[-1]
    e_lo, e_hi = int(np.floor(np.log10(n_lo))), int(np.floor(np.log10(n_hi)))
    m_lo, m_hi = n_lo / 10 ** e_lo, n_hi / 10 ** e_hi
    ax.set_title("A  i.i.d. — known truth (ESS = n)", fontsize=7.4, loc="left")
    # ⚠⚠ 逐点百分比**不能**标在点旁（实测压住重复点与理论线）→ 集中到右下空白区
    # ⚠ 刻度/范围文字一律**现算**（禁硬写 "10^2"/"10^4"）：数据规模一改就失配
    ax.annotate("mean dev   " + " / ".join("%.1f%%" % d for d in devs)
                + "   (n = %.0f×10$^{%d}$ … %.1f×10$^{%d}$)"
                % (m_lo, e_lo, m_hi, e_hi)
                + "\nworst rep  %+.0f%%  (n=%.0f×10$^{%d}$)   →   %+.1f%%  (n=%.1f×10$^{%d}$)"
                % (worst[0], m_lo, e_lo, worst[-1], m_hi, e_hi),
                (0.97, 0.03), xycoords="axes fraction", fontsize=5.4,
                color="#444444", ha="right", va="bottom",
                bbox=dict(boxstyle="round,pad=0.22", fc="white",
                          ec="#cccccc", lw=0.4, alpha=0.92))
    sci_legend(ax, loc="upper left")
    despine(ax)


def _p_ar1(ax, ar1):
    rhos = sorted(ar1)
    tv = np.array([ar1[r]["ess_theory"] for r in rhos])
    ax.plot(tv, tv, "-", color=C_REF, lw=1.0, zorder=1, label="theory")
    for key, col, mk, lab in (("ours", C_BAD, "o", "this impl."),
                              ("indep", C_OK, "s", "indep. IPS"),
                              ("bm", C_BM, "^", "batch-means")):
        y = [float(np.mean(ar1[r][key])) for r in rhos]
        e = [float(np.std(ar1[r][key], ddof=1)) if len(ar1[r][key]) > 1 else 0.0
             for r in rhos]
        ax.errorbar(tv, y, yerr=e, fmt=mk, ms=4.2, color=col, lw=0.9,
                    capsize=1.8, elinewidth=0.7, zorder=3, label=lab)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("theoretical ESS")
    ax.set_ylabel("estimated ESS")
    ax.set_title("B  AR(1) — known $\\rho$", fontsize=7.4, loc="left")
    sci_legend(ax, loc="upper left")
    despine(ax)
    txt = []
    for i, r in enumerate(rhos):
        v = float(np.mean(ar1[r]["ours"]))
        # ⚠ 不能直接 "%.0f%%"：−0.42% 会被格式化成 "-0%"，读起来像出错
        iv = int(round((v - tv[i]) / tv[i] * 100))
        txt.append("0%" if iv == 0 else "%+d%%" % iv)
    _note_box(ax, "dev. (this impl.): " + " / ".join(txt))


def _p_coda(ax, ar1, coda_vals):
    if coda_vals is None:
        ax.axis("off")
        # ⚠ 图内文字一律 **ASCII** —— matplotlib 默认字体无 CJK 字形，
        # 中文会渲染成方框（实测有大量 missing-glyph 警告）。
        ax.text(0.5, 0.5,
                "C  coda cross-check\n\nskipped — R / coda not available\n"
                "(see tests/test_ess_implementation.py)",
                ha="center", va="center", fontsize=6.2, color="#888888",
                transform=ax.transAxes,
                bbox=dict(boxstyle="round,pad=0.5", fc="#f7f7f7",
                          ec="#cccccc", lw=0.5))
        return
    rhos = sorted(ar1)
    o, c, rl = [], [], []
    for r in rhos:
        for i, x in enumerate(ar1[r]["reps"]):
            o.append(calculate_ess(x))
            rl.append(r)
    c = np.array(coda_vals, dtype=float)[:len(o)]
    o = np.array(o, dtype=float)[:len(c)]
    rl = np.array(rl, dtype=float)[:len(c)]
    lo, hi = min(o.min(), c.min()) * 0.7, max(o.max(), c.max()) * 1.4
    ax.plot([lo, hi], [lo, hi], "-", color=C_REF, lw=1.0, zorder=1, label="y = x")
    # ⚠⚠ ρ 是**数值型有序**变量 ⇒ 走 SEQ_CMAP，不能用 OKABE_ITO 语义色：
    # 否则 C_BAD 在 panel A/B/D 表示「本实现」、在 panel C 表示「ρ=0.99」
    # —— 同图同色两义，违反「同图颜色要么全语义、要么全数值」。
    for i, (r, mk) in enumerate(zip(rhos, ("o", "s", "^", "D", "v"))):
        m = np.abs(rl - r) < 1e-9
        ax.scatter(o[m], c[m], s=20, color=_CM(0.02 + 0.17 * i), marker=mk,
                   zorder=3, edgecolor="white", lw=0.4, label=r"$\rho$=%.2f" % r)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.set_xlabel("ESS (this impl.)")
    ax.set_ylabel("ESS (coda, AR spectrum)")
    ax.set_title("C  cross-check vs coda (AR spectrum vs IPS)", fontsize=7.4,
                 loc="left")
    sci_legend(ax, loc="upper left", ncol=2)
    despine(ax)
    rel = np.abs(o - c) / c * 100
    i_max = int(np.argmax(rel))
    # ⚠⚠ 最大值所在 ρ **必须现算**（原版硬写 "(ρ=0.99)"，而实测最大值在 ρ=0.95）
    _note_box(ax, "median |Δ| = %.1f%%\nmax = %.1f%%  (ρ=%.2f)"
              % (np.median(rel), rel.max(), rl[i_max]))


def systematic_note(rel):
    """把「本实现 vs 独立实现」的相对差数组变成一句**方向正确**的注记。

    ⚠⚠ 必须判符号：直接用 `abs(rel).min()/max()` 在混合符号时会把
    「有高有低」静默写成「一律偏低」（与判据方向写反同类）。
    """
    lo, hi = float(np.min(rel)), float(np.max(rel))
    if abs(lo) < 1e-9 and abs(hi) < 1e-9:
        return "no difference  (ours == indep. IPS)"
    if hi <= 0:
        return ("systematic:  ours %.1f–%.1f%% lower   "
                "(biased vs unbiased ACF)" % (abs(hi), abs(lo)))
    if lo >= 0:
        return ("systematic:  ours %.1f–%.1f%% higher   "
                "(biased vs unbiased ACF)" % (lo, hi))
    return ("mixed sign:  %+.1f%% … %+.1f%%   "
            "(no systematic direction)" % (lo, hi))


def _p_real(ax, chains, keys, n=None):
    """真实数据：本实现 vs 独立重写 IPS（逐参数，拼接 pooled 口径）。"""
    lab = list(keys)
    pe, pi = [], []
    for k in lab:
        pooled = np.concatenate([np.asarray(c[k], dtype=float) for c in chains])
        pe.append(calculate_ess(pooled))
        pi.append(ess_ips_independent(pooled))
    pe = np.array(pe)
    pi = np.array(pi)
    y = np.arange(len(lab))
    ax.scatter(pe, y, s=26, color=C_BAD, zorder=3, label="this impl.")
    ax.scatter(pi, y, s=13, marker="x", color=C_OK, lw=0.9, zorder=4,
               label="indep. IPS")
    ax.set_yticks(y)
    ax.set_yticklabels(lab, fontsize=5.2)
    ax.invert_yaxis()
    ax.set_xscale("log")
    ax.set_xlabel("pooled ESS")
    ax.set_title("D  real data — impl. cross-check", fontsize=7.4, loc="left")
    # ⚠⚠ 注记**不能**放左上：前两行是 posterior / prior，其 pooled ESS 就在轴
    # 左端，注记框会直接盖住数据点（实测第一版就是这样）。正解 = 顶部额外留
    # 一行空白带（y 已反转，故上界写成负值）再放注记。
    ax.set_ylim(len(lab) - 0.5, -1.5)
    sci_legend(ax, loc="lower right")
    despine(ax)
    rel = (pe - pi) / pi * 100
    # 本实现**系统性低于**独立实现 —— 已知的归一化差异（Tracer 的 (n−lag)
    # 归一 ACF vs (n) 归一），不是随机噪声，故给**区间**而非单个 max。
    ax.annotate(systematic_note(rel), (0.03, 0.985), xycoords="axes fraction",
                fontsize=5.2, color="#444444", ha="left", va="top")


def plot_ess_verification(out_pdf, chains=None, keys=None, rscript="Rscript"):
    """图①：ESS 实现金标准验证。返回统计 dict（供测试/日志断言用）。"""
    iid, ar1, n_eff = generate_gold_standard()

    series = [x for r in sorted(ar1) for x in ar1[r]["reps"]]
    coda_vals = coda_ess_on_series(series, rscript=rscript)

    apply_sci_style(base_size=6.8)
    fig = plt.figure(figsize=(COL_W_1, COL_W_1 * 0.62))
    # ⚠ 主标题用 fig.text 放画布顶缘留白，且 top 下移 —— 否则与面板 A 的
    # axes 标题行重叠（实测：suptitle y=0.985 + top=0.945 会撞）。
    gs = fig.add_gridspec(2, 2, hspace=0.62, wspace=0.30,
                          left=0.075, right=0.985, top=0.905, bottom=0.095)
    _p_iid(fig.add_subplot(gs[0, 0]), iid)
    _p_ar1(fig.add_subplot(gs[0, 1]), ar1)
    _p_coda(fig.add_subplot(gs[1, 0]), ar1, coda_vals)
    if chains is not None and keys:
        _p_real(fig.add_subplot(gs[1, 1]), chains, keys)
    else:
        axd = fig.add_subplot(gs[1, 1])
        axd.axis("off")
        axd.text(0.5, 0.5,
                 "D  real data\n\nnot available\n(no per-chain logs given)",
                 ha="center", va="center", fontsize=6.2, color="#888888",
                 transform=axd.transAxes,
                 bbox=dict(boxstyle="round,pad=0.5", fc="#f7f7f7",
                           ec="#cccccc", lw=0.5))

    fig.text(0.075, 0.972,
             "ESS implementation verification — synthetic truth, external "
             "library, real data",
             fontsize=8.4, ha="left", va="top")

    devs = {r: (float(np.mean(ar1[r]["ours"])) - ar1[r]["ess_theory"])
            / ar1[r]["ess_theory"] * 100 for r in sorted(ar1)}
    stats = {
        "n_eff": n_eff,
        "ar1_dev_pct": devs,
        "iid_mean_dev_pct": {n: float((np.mean(iid[n]) - n) / n * 100) for n in sorted(iid)},
        "coda_available": coda_vals is not None,
    }
    if coda_vals is not None:
        stats["coda_median_abs_pct"] = float(np.median(
            np.abs(np.array(series and [calculate_ess(x) for x in series])[:len(coda_vals)]
                   - np.array(coda_vals)) / np.array(coda_vals) * 100))
    _save(fig, out_pdf)
    return stats


# ══════════════════════════════════════════════════════════════════
# 图②  本 run 收敛诊断
# ══════════════════════════════════════════════════════════════════
def _p_rhat(ax, lab, x, up):
    ax.axvspan(RHAT_PRACTICAL, max(max(up) * 1.15, 2.0), color="#fbe3e0", zorder=0)
    ax.axvline(RHAT_PRACTICAL, color=C_REF, lw=1.0, ls="--", zorder=1)
    ax.annotate("R-hat = %.1f" % RHAT_PRACTICAL, (RHAT_PRACTICAL, 0.995),
                xycoords=("data", "axes fraction"), fontsize=5.2, color=C_REF,
                ha="left", va="top", xytext=(3, 0), textcoords="offset points")
    # ⚠ 1.01 与 1.1 在 xlim=[1, 7.6] 下仅隔 ~1.5% 轴宽，两条线会视觉重合
    # → 只画**操作性判据 1.1**（1.01 在注记里说明），避免假精度。
    ax.annotate("not converged", (0.62, 0.995), xycoords="axes fraction",
                fontsize=5.4, color="#a05050", ha="center", va="top")
    y = np.arange(len(lab))
    if up is not None:
        ax.hlines(y, x, up, color="#bbbbbb", lw=0.8, zorder=2)
    ax.scatter(x, y, s=26, color=C_BAD, zorder=4, edgecolor="white", lw=0.4)
    ax.set_yticks(y)
    ax.set_yticklabels(lab, fontsize=5.2)
    ax.set_xlabel("R-hat  (Gelman–Rubin, independent of ESS)")
    n_ab = int((x > RHAT_PRACTICAL).sum())
    n_st = int((x > RHAT_STRICT).sum())
    ax.set_title("A  R-hat — %d/%d above %.1f" % (n_ab, len(lab), RHAT_PRACTICAL),
                 fontsize=7.4, loc="left")
    # ⚠⚠ 必须用**对数轴**：线性轴下多数点挤在 1.03–1.49（占轴宽 <7%）根本分不开。
    # ⚠ 对数轴默认只给 10^0 一个主刻度 → 必须**显式**给刻度，否则轴上无数值。
    # ⚠ 刻度取值受**标签宽度**约束：实测 [1.0, 1.25, 1.5, …] 会让 "1.25" 与
    # "1.5" 标签重叠 5.3 px²（版式体检抓到，肉眼看不见）→ 用 [1.0, 1.2, 1.5, …]。
    ax.set_xscale("log")
    ax.minorticks_off()
    ax.set_xticks([1.0, 1.2, 1.5, 2.0, 3.0, 5.0, 7.5])
    ax.set_xticklabels(["1.0", "1.2", "1.5", "2", "3", "5", "7.5"])
    ax.set_xlim(1.0, max(max(up) * 1.06, 2.2))
    ax.set_ylim(-1.35, len(lab) - 0.3)     # 底部留白给注记框，避免压住最后一行
    # ⚠ 标题/注记必须与实测一致：最小值可能 **未超 1.1**，不能写 "all params fail"
    _note_box(ax, "worst %.2f  ·  min %.2f\n%d/%d > %.1f  (practical)\n"
                  "%d/%d > %.2f  (strict)"
              % (float(np.max(x)), float(np.min(x)), n_ab, len(lab),
                 RHAT_PRACTICAL, n_st, len(lab), RHAT_STRICT))
    despine(ax)


def _p_ess_bars(ax, lab, sm, po):
    y = np.arange(len(lab))
    h = 0.36
    ax.barh(y + h / 2, sm, height=h, color="#c9d9e8", edgecolor="#8fa8c0",
            lw=0.4, zorder=2, label="sum of per-chain  (coda convention)")
    ax.barh(y - h / 2, po, height=h, color=C_BAD, edgecolor="white",
            lw=0.4, zorder=3, label="pooled / concatenated  (ours)")
    ax.axvline(MIN_ESS, color=C_REF, lw=1.0, ls="--", zorder=4)
    # ⚠ 注记放**轴顶**（轴底有图例；条形顶端在中部）
    ax.annotate("threshold ESS = %g" % MIN_ESS, (MIN_ESS, 0.995),
                xycoords=("data", "axes fraction"), fontsize=5.2, color=C_REF,
                ha="left", va="top", xytext=(3, 0), textcoords="offset points")
    ax.set_yticks(y)
    ax.set_yticklabels(lab, fontsize=5.2)
    ax.set_xscale("log")
    # ⚠ 对数轴上 barh 的矩形是「从 0 到 v」，0 在 log 轴上是 −∞ → matplotlib 把
    # 条裁到轴左界，看起来像从轴左缘起画。若不显式给上界，最长的条会被右轴
    # **静默裁掉**（实测第一版 skyline.popSize2 的蓝条正好顶到右脊线）。
    _lo = min(po.min(), sm.min()) * 0.6
    _hi = max(po.max(), sm.max(), MIN_ESS) * 1.7
    ax.set_xlim(_lo, _hi)
    ax.set_xlabel("ESS  (log scale)")
    ax.set_title("B  ESS — convention gap is the point", fontsize=7.4, loc="left")
    ax.set_ylim(-1.35, len(lab) - 0.3)
    sci_legend(ax, loc="lower right")
    despine(ax)


def _p_pool_ratio(ax, lab, ra):
    y = np.arange(len(lab))
    ax.axvline(1.0, color=C_REF, lw=1.0, ls="--", zorder=1)
    ax.annotate("= 1  good mixing", (0.985, len(lab) - 0.4), fontsize=5.2,
                color=C_REF, ha="right", va="bottom")
    ax.barh(y, ra, height=0.55, color=C_BAD, edgecolor="white", lw=0.4, zorder=3)
    for i, v in enumerate(ra):
        ax.annotate("%.3f" % v, (v, i), xytext=(3, 0), textcoords="offset points",
                    fontsize=5.0, va="center", color="#444444")
    ax.set_yticks(y)
    ax.set_yticklabels(lab, fontsize=5.2)
    ax.set_xlim(0, max(1.22, float(np.max(ra)) * 1.15))
    ax.set_xlabel("pooled ESS / Σ(per-chain ESS)")
    # ⚠ 禁硬写：标题与注记数字一律由数据现算
    ax.set_title("C  mixing collapse — %.2f–%.2f" % (ra.min(), ra.max()),
                 fontsize=7.4, loc="left")
    ax.set_ylim(-1.35, len(lab) - 0.3)
    _note_box(ax, "%.0f–%.0f× loss\nwhen concatenated"
              % (1.0 / ra.max(), 1.0 / ra.min()))
    despine(ax)


def plot_convergence_diagnostics(chain_logs, out_pdf, burnin_frac=0.10,
                                 rscript="Rscript", subtitle=None):
    """图②：本 run 的收敛诊断（R-hat / ESS 口径 / 池化损失）。返回统计 dict。

    Parameters
    ----------
    chain_logs : list[str]  各**单链** BEAST log（不是 merged.log —— R-hat 需要
                            链间方差，用拼接后的 log 算不出链间差异）
    out_pdf : str           输出 PDF（同目录同时写 PNG）
    """
    if len(chain_logs) < 2:
        raise ValueError("R-hat 至少需要 2 条链（当前 %d 条）" % len(chain_logs))

    chains = [read_log(p, burnin_frac)[1] for p in chain_logs]
    keys = key_columns(chains)
    if not keys:
        raise ValueError("各链没有共同的关键参数列，无法出图")
    by_key, n = truncate_common(chains, keys)

    coda = rhat_via_coda(chains, keys, rscript=rscript)
    rows = []
    for k in keys:
        per = [calculate_ess(by_key[k][i]) for i in range(len(chains))]
        pooled = np.concatenate(by_key[k])
        pe = calculate_ess(pooled)
        pi = ess_ips_independent(pooled)
        if coda and k in coda:
            rh, up = coda[k]["rhat"], coda[k]["upper"]
        else:
            rh = rhat_python(by_key[k])
            up = rh                        # ⚠ 无上界来源时不画上界（点估计即上界）
        if rh is None:
            continue
        rows.append({"param": k, "rhat": rh, "upper": up,
                     "sum_perchain": float(sum(per)), "pooled": pe,
                     "pooled_indep": pi, "pool_over_sum": pe / sum(per)})
    if not rows:
        raise ValueError("没有可算 R-hat 的参数列")

    rows.sort(key=lambda r: r["rhat"])
    lab_r = [r["param"] for r in rows]
    x = np.array([r["rhat"] for r in rows])
    up = np.array([r["upper"] for r in rows])

    rows_b = sorted(rows, key=lambda r: r["pooled"])
    lab_b = [r["param"] for r in rows_b]
    sm = np.array([r["sum_perchain"] for r in rows_b])
    po = np.array([r["pooled"] for r in rows_b])

    rows_c = sorted(rows, key=lambda r: r["pool_over_sum"])
    lab_c = [r["param"] for r in rows_c]
    ra = np.array([r["pool_over_sum"] for r in rows_c])

    apply_sci_style(base_size=6.8)
    fig = plt.figure(figsize=(COL_W_1, COL_W_1 * 0.72))
    gs = fig.add_gridspec(1, 3, wspace=0.62, left=0.105, right=0.985,
                          top=0.885, bottom=0.105)
    _p_rhat(fig.add_subplot(gs[0, 0]), lab_r, x, up)
    _p_ess_bars(fig.add_subplot(gs[0, 1]), lab_b, sm, po)
    _p_pool_ratio(fig.add_subplot(gs[0, 2]), lab_c, ra)

    if subtitle is None:
        subtitle = ("%d chains × %d post-burnin samples   ·   "
                    "R-hat source: %s   ·   B and C are the same ESS data "
                    "in two views, not two evidences"
                    % (len(chains), n,
                       "coda::gelman.diag" if coda
                       else "python fallback (under-reports R-hat)"))
    fig.text(0.105, 0.975,
             "MCMC convergence — R-hat and the ESS convention gap",
             fontsize=8.4, ha="left", va="top")
    fig.text(0.105, 0.925, subtitle, fontsize=6.2, color="#555555",
             ha="left", va="top")

    stats = {
        "n_chains": len(chains), "n_samples": n, "n_params": len(rows),
        "rhat_source": "coda" if coda else "python",
        "rhat_min": float(x.min()), "rhat_min_param": lab_r[0],
        "rhat_max": float(x.max()), "rhat_max_param": lab_r[-1],
        "n_gt_practical": int((x > RHAT_PRACTICAL).sum()),
        "n_gt_strict": int((x > RHAT_STRICT).sum()),
        "worst_pooled_ess": float(po.min()),
        "min_pool_over_sum": float(ra.min()), "max_pool_over_sum": float(ra.max()),
        "params": rows,
    }
    _save(fig, out_pdf)
    return stats


# ══════════════════════════════════════════════════════════════════
# CLI
# ══════════════════════════════════════════════════════════════════
def main(argv=None):
    ap = argparse.ArgumentParser(
        description="ESS 实现验证 + MCMC 收敛诊断出图（管线模块）")
    ap.add_argument("--chains", default="",
                    help="逗号分隔的单链 BEAST log 路径（≥2 条才算 R-hat）")
    ap.add_argument("--outdir", default=".", help="输出目录")
    ap.add_argument("--burnin", type=float, default=0.10, help="burn-in 比例")
    ap.add_argument("--verify", action="store_true",
                    help="额外产出金标准验证图（可同时接 --chains 做 panel D）")
    ap.add_argument("--verify-only", action="store_true", help="只出验证图")
    ap.add_argument("--rscript", default="Rscript", help="Rscript 可执行文件")
    ap.add_argument("--prefix", default="ess",
                    help="输出文件名前缀（默认 ess）")
    args = ap.parse_args(argv)

    os.makedirs(args.outdir, exist_ok=True)
    logs = [p for p in args.chains.split(",") if p.strip()]
    rc = 0

    if not args.verify_only:
        if len(logs) < 2:
            print("⚠ 收敛诊断跳过：需要 ≥2 条单链 log（当前 %d 条）" % len(logs))
            rc = 2
        else:
            st = plot_convergence_diagnostics(
                logs, os.path.join(args.outdir,
                                   "%s_convergence_diagnostics.pdf" % args.prefix),
                burnin_frac=args.burnin, rscript=args.rscript)
            print("  R-hat[%s] min=%.4f(%s) max=%.4f(%s)  %d/%d > %.1f  "
                  "worst pooled ESS=%.1f"
                  % (st["rhat_source"], st["rhat_min"], st["rhat_min_param"],
                     st["rhat_max"], st["rhat_max_param"],
                     st["n_gt_practical"], st["n_params"], RHAT_PRACTICAL,
                     st["worst_pooled_ess"]))

    if args.verify or args.verify_only:
        chains = keys = None
        if len(logs) >= 1:
            ch = [read_log(p, args.burnin)[1] for p in logs]
            kk = key_columns(ch)
            if kk:
                chains, keys = ch, kk
        st2 = plot_ess_verification(
            os.path.join(args.outdir, "%s_implementation_verification.pdf" % args.prefix),
            chains=chains, keys=keys, rscript=args.rscript)
        print("  验证: N_EFF=%d  coda=%s  AR(1) dev=%s"
              % (st2["n_eff"], st2["coda_available"],
                 {k: round(v, 2) for k, v in st2["ar1_dev_pct"].items()}))
    return rc


if __name__ == "__main__":
    sys.exit(main())
