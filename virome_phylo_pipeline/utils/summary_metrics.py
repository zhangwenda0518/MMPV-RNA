"""summary_metrics.py — 统一「头条指标」口径（2026-09-15 修复 P0-A / P0-B）

## 为什么需要这个模块

修复前，同一批结果会产出两份产物，而它们对**同名列**取了**不同来源**：

| 指标 | phylo_summary.csv | phylo_report.html |
|---|---|---|
| 速率 | `td_mean_rate` → 退回 `pic_rate` | `metrics['beta']`（TreeTime RTT 斜率） |
| 拟合 | `pic_r` | `metrics['r_squared']`（RTT） |
| 显著性 | `td_relaxed_clock_p` | `metrics['p_value']`（RTT t 检验） |
| TMRCA | `td_tmrca` | `beast_tmrca` → 退回 `td_tmrca` |

列名一模一样、数值不同、且没有任何来源标注 —— 同一病毒两份产物互相打脸。
另外 `Beta` 列写的是 `td_mean_rate or pic_rate`，把 treedater 松弛钟**均值速率**
与 PIC 严格钟**回归斜率**（两个不同估计量、不同假设）混进同一列，
treedater 失败时静默退回 PIC，跨病毒比较分子钟速率时混入异质方法。

## 现在的约定

`headline_metrics(metrics)` 是**唯一**取数入口，CSV 与 HTML 都调它；
返回的每个指标都带 `*_source` 标签，写作产物时必须一起落盘/显示。

优先级（**treedater 优先**，因为它是 BEAST 同族的正式钟模型估计；
PIC 是替代方案；TreeTime RTT 只作交叉核对）：
  速率    td_mean_rate → pic_rate
  拟合    pic_r        → r_squared
  显著性  td_relaxed_clock_p → p_value
  TMRCA   beast_tmrca  → td_tmrca

注意：`headline` 值缺失时返回 `""`（空串），**不返回 0**——
「没算出来」与「算出来是 0」必须可区分。
"""

from __future__ import annotations

from typing import Any, Dict, Tuple

# (metrics 键, 来源标签)，按优先级排列
_RATE_SOURCES: Tuple[Tuple[str, str], ...] = (
    ("td_mean_rate", "treedater"),
    ("pic_rate", "PIC"),
)
_R2_SOURCES: Tuple[Tuple[str, str], ...] = (
    ("pic_r", "PIC"),
    ("r_squared", "RTT"),
)
_P_SOURCES: Tuple[Tuple[str, str], ...] = (
    ("td_relaxed_clock_p", "treedater"),
    ("p_value", "RTT"),
)
_TMRCA_SOURCES: Tuple[Tuple[str, str], ...] = (
    ("beast_tmrca", "BEAST"),
    ("td_tmrca", "treedater"),
)

# 备用列：只落盘供交叉核对，不参与「头条」判定
_ALT_KEYS = (
    ("RTT_beta", "beta", "RTT 斜率"),
    ("RTT_r_squared", "r_squared", "RTT"),
    ("RTT_p_value", "p_value", "RTT"),
    ("PIC_rate", "pic_rate", "PIC"),
    ("treedater_rate", "td_mean_rate", "treedater"),
    ("treedater_p", "td_relaxed_clock_p", "treedater"),
    ("BEAST_tmrca", "beast_tmrca", "BEAST"),
    ("treedater_tmrca", "td_tmrca", "treedater"),
)


def _pick(metrics: Dict[str, Any], sources) -> Tuple[Any, str]:
    """按优先级取第一个"真正有值"的候选。空串/None 都视为缺。"""
    for key, label in sources:
        val = metrics.get(key)
        if val is None or val == "":
            continue
        return val, label
    return "", ""


def headline_metrics(metrics: Dict[str, Any] | None) -> Dict[str, Any]:
    """CSV 汇总与 HTML 报告共用的头条指标提取（含来源标签）。

    返回键：
      rate / rate_source     分子钟速率（-subs/site/year 量级）
      r2   / r2_source       拟合优度
      p    / p_source        钟检验 / 回归显著性
      tmrca / tmrca_source   最近共同祖先时间
      alt                    备用口径 dict {列名: (值, 来源)}，供交叉核对
      any_clock              是否至少拿到一个速率类指标
    """
    m = metrics or {}
    rate, rate_src = _pick(m, _RATE_SOURCES)
    r2, r2_src = _pick(m, _R2_SOURCES)
    p, p_src = _pick(m, _P_SOURCES)
    tmrca, tmrca_src = _pick(m, _TMRCA_SOURCES)
    alt = {}
    for col, key, label in _ALT_KEYS:
        val = m.get(key)
        if val is None or val == "":
            continue
        alt[col] = (val, label)
    return {
        "rate": rate,
        "rate_source": rate_src,
        "r2": r2,
        "r2_source": r2_src,
        "p": p,
        "p_source": p_src,
        "tmrca": tmrca,
        "tmrca_source": tmrca_src,
        "alt": alt,
        "any_clock": bool(rate != "" or r2 != "" or p != ""),
    }


def fmt_metric(value: Any, source: str, dash: str = "—") -> str:
    """渲染成 `值 (来源)`；无值时返回 dash。HTML 与文本日志共用，避免格式漂移。"""
    if value is None or value == "":
        return dash
    txt = f"{value}"
    # 浮点统一保留 6 位有效数字，避免 0.00012300000000000001 之类的噪声
    try:
        f = float(value)
        if f == 0 or abs(f) >= 1e-3:
            txt = f"{f:g}"
        else:
            txt = f"{f:.3g}"
    except (TypeError, ValueError):
        pass
    return f"{txt} ({source})" if source else txt


def headline_columns(metrics: Dict[str, Any] | None) -> Dict[str, Any]:
    """给 phylo_summary.csv 用的扁平列（保持旧列名 + 新增来源列 + 备用口径列）。

    旧列名 Beta / R_squared / P_value / TMRCA **保持不变**（下游脚本与历史
    对照依赖它们），但取值一律来自 headline_metrics，因此与 HTML 报告一致；
    来源写在 *_source 列，备用口径列用 alt 前缀展开。
    """
    hm = headline_metrics(metrics)
    row = {
        "Beta": hm["rate"],
        "Beta_source": hm["rate_source"],
        "R_squared": hm["r2"],
        "R_squared_source": hm["r2_source"],
        "P_value": hm["p"],
        "P_value_source": hm["p_source"],
        "TMRCA": hm["tmrca"],
        "TMRCA_source": hm["tmrca_source"],
    }
    for col, (val, _label) in hm["alt"].items():
        row[f"alt_{col}"] = val
    return row
