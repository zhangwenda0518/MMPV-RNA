#!/usr/bin/env python3
"""convergence_gate.py — 合并产物的 MCMC 收敛门禁 (2026-10-13)

背景 (为什么需要这个模块)
========================
`utils/beast_parser.py` 早已逐参数算好 `converged = ess >= 200`, 并写进
`parameter_summary.csv`。但**下游没有任何一处据此设门禁**:
`merge_results.merge()` 只检查 `parameter_summary.json` **是否存在**就打印 ✅,
于是 ESS=4–9 的未收敛链被静默采信, 一路流进报告与图。

GCVA 实测即为此例: genome 链 `parameter_summary.csv` 里 posterior ESS=9.1、
ucld.mean ESS=4.5, 全部 `converged=False`, 但合并照样"成功"。

本模块把「合并产物能不能信」收敛成**单一判据 + 单一落盘**, 供 merge_results
与下游 (phylo_pipeline / 报告) 共用, 避免各写各的口径。

判据
====
对 `parameter_summary.csv` 的关键参数取最小 ESS, 与 `MIN_ESS`(=200) 比:
  · 全部 key 参数 ESS >= 200  → verdict = "converged"
  · 有任何 key 参数 ESS <  200 → verdict = "not_converged"
  · 无 key 参数行 / 文件缺失   → verdict = "unknown"

⚠ 语义边界 (务必写进下游 caveat):
  verdict 只回答 "MCMC 混合是否充分", **不回答 "数据能否定年"**。
  两者独立 —— 一个零时间信号的数据集即便 ESS 达标, 也不可做分子钟定年
  (反面同理: ESS 不达标的链, 无论数据有无信号都不可采信其数值)。
  故本门禁**只降级/标注**, 不"判死"整个管线 (与 sample_floor 的工程哲学一致:
  要么有结果, 要么有明确告警, 不允许静默产出)。

落盘
====
`<merged>/convergence_gate.json`, 结构:
    {
      "verdict": "converged" | "not_converged" | "unknown",
      "min_ess": 200.0,
      "worst_ess": 4.5,               # 关键参数里的最小 ESS
      "worst_param": "ucld.mean",
      "n_key_params": 3,
      "key_params": {"ucld.mean": 4.5, ...},
      "reason": "…人读说明…",
      "generated_at": "YYYY-MM-DD HH:MM:SS"
    }
"""
import csv
import json
import os
import time

# 收敛判据的唯一来源 (与 beast_parser 的 ess_threshold 默认值一致)
from utils.sample_floor import MIN_ESS

# 关键参数 (定年结论实际依赖的那几个)。用小写子串/精确匹配双通道,
# 与 beast_parser 的 key_names 保持同源语义。
KEY_PARAM_PATTERNS = (
    "posterior",           # 后验 (全局混合)
    "prior",               # 先验
    "likelihood",          # 似然
    "ucld.mean",           # 分子钟速率 (定年核心)
    "clock.rate",
    "clockrate",
    "ucld.stdev",          # 松弛钟偏离度
    "trees.treelength",    # 树长 (与时间尺度耦合)
    "treemodel.rootheight",
    "tmrca",
    "coalescent",
)


def _is_key_param(name):
    """参数名是否属于关键集 (大小写不敏感, 子串匹配)。"""
    n = (name or "").lower()
    return any(pat in n for pat in KEY_PARAM_PATTERNS)


def evaluate(merged_dir, min_ess=None, log=None):
    """评估一个 merged/ 目录的收敛状态, 写 convergence_gate.json 并返回 dict。

    Parameters
    ----------
    merged_dir : str  含 parameter_summary.csv 的目录
    min_ess : float   阈值 (默认 MIN_ESS=200)
    log : optional    给 LogCollector 就写它的 emit/warning

    Returns
    -------
    dict  见模块 docstring 的落盘结构 (始终返回, 不抛异常)
    """
    thr = float(MIN_ESS if min_ess is None else min_ess)
    out = {
        "verdict": "unknown",
        "min_ess": thr,
        "worst_ess": None,
        "worst_param": None,
        "n_key_params": 0,
        "key_params": {},
        "reason": "",
        "generated_at": time.strftime("%Y-%m-%d %H:%M:%S"),
    }

    csv_path = os.path.join(merged_dir, "parameter_summary.csv")
    if not os.path.exists(csv_path):
        out["reason"] = ("parameter_summary.csv 不存在 — 合并未完成或空汇总, "
                         "无法判定收敛")
        _emit(log, f"  ⚠️ [convergence_gate] {out['reason']}")
        _write(merged_dir, out)
        return out

    try:
        with open(csv_path, encoding="utf-8", errors="replace") as f:
            rows = list(csv.DictReader(f))
    except OSError as e:
        out["reason"] = f"读取 parameter_summary.csv 失败: {e}"
        _emit(log, f"  ⚠️ [convergence_gate] {out['reason']}")
        _write(merged_dir, out)
        return out

    key_ess = {}
    for r in rows:
        name = (r.get("name") or "").strip()
        if not name or not _is_key_param(name):
            continue
        try:
            e = float(r.get("ess"))
        except (TypeError, ValueError):
            continue
        # 同名取最小 (防重复行)
        key_ess[name] = min(key_ess.get(name, e), e)

    if not key_ess:
        out["reason"] = ("parameter_summary.csv 无关键参数行 "
                         "(posterior/ucld.mean/… 全缺) — 无法判定收敛")
        _emit(log, f"  ⚠️ [convergence_gate] {out['reason']}")
        _write(merged_dir, out)
        return out

    worst_param = min(key_ess, key=key_ess.get)
    worst_ess = key_ess[worst_param]
    out["key_params"] = {k: round(v, 1) for k, v in sorted(key_ess.items())}
    out["n_key_params"] = len(key_ess)
    out["worst_ess"] = round(worst_ess, 1)
    out["worst_param"] = worst_param

    if worst_ess >= thr:
        out["verdict"] = "converged"
        out["reason"] = (f"关键参数最小 ESS={worst_ess:.1f} ≥ {thr:g} "
                         f"(最差: {worst_param}) — MCMC 混合充分")
        _emit(log, f"  ✅ [convergence_gate] {out['reason']}")
    else:
        out["verdict"] = "not_converged"
        out["reason"] = (f"关键参数最小 ESS={worst_ess:.1f} < {thr:g} "
                         f"(最差: {worst_param}) — MCMC 未收敛, "
                         f"该批参数估计不可采信 (数值仅供探索)")
        _emit(log, f"  ❌ [convergence_gate] {out['reason']}")

    _write(merged_dir, out)
    return out


def _write(merged_dir, payload):
    try:
        with open(os.path.join(merged_dir, "convergence_gate.json"),
                  "w", encoding="utf-8") as f:
            json.dump(payload, f, indent=2, ensure_ascii=False)
    except OSError:
        pass


def _emit(log, msg):
    if log is None:
        print(msg, flush=True)
        return
    for attr in ("emit", "info", "warning"):
        fn = getattr(log, attr, None)
        if callable(fn):
            try:
                fn(msg)
            except Exception:
                print(msg, flush=True)
            return
    print(msg, flush=True)


if __name__ == "__main__":
    import sys
    if len(sys.argv) < 2:
        print("用法: python -m utils.convergence_gate <merged_dir>")
        sys.exit(2)
    r = evaluate(sys.argv[1])
    sys.exit(0 if r["verdict"] == "converged" else 1)
