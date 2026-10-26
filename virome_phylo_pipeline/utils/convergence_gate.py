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
对 `parameter_summary.csv` 的**关键参数**取最小 ESS, 与 `MIN_ESS`(=200) 比:
  · 全部 key 参数 ESS >= 200  → verdict = "converged"
  · 有任何 key 参数 ESS <  200 → verdict = "not_converged"
  · 无 key 参数行 / 文件缺失   → verdict = "unknown"

关键参数集 = `KEY_PARAM_PATTERNS` 子串匹配 (速率 / 松弛度 / 树高 / 树先验 /
全局混合)。⚠ 该集合须覆盖**实际产物里出现的全部定年相关列** —— 漏一条就会
低估最差 ESS (2026-10-14 实测: 曾漏 `coefficientOfVariation`, 上报 4.5 而真值 4.4)。

⚠⚠ `not_converged` 的语义 = **"ESS 不达标"**, 不是 "链未收敛"。两者不可混:
ESS 低既可能源于链未收敛, 也可能仅源于 logEvery 过大导致样本数太少
(ESS 上限被样本数压死)。本门禁**不区分**这两种成因, 故 caveat 措辞
不得断言 "未收敛" (与 §硬规则④ "过度归因" 同一纪律)。

⚠ 语义边界 (下游须理解, 但**不写进 caveat 文案**):
  verdict 只回答 "MCMC 混合是否充分", **不回答 "数据能否定年"**。
  两者独立 —— 一个零时间信号的数据集即便 ESS 达标, 也不做分子钟定年;
  反之 ESS 不达标的链, 其数值也不宜作定量结论。
  ⚠⚠ 但这两句是**给用者的判读指引**, 不是产物该自我宣告的话 (用户 2026-10-14 裁定:
  caveat 只记事实, 不写"不可采信/仅供探索"这类价值判断 —— 见 §判据下方措辞纪律)。
  故本门禁**只记录数值 + 判定**, 不"判死"整个管线 (与 sample_floor 的工程哲学一致:
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

# 关键参数 (定年结论实际依赖的那几个)。用小写子串匹配, 与 beast_parser 的
# key_names 保持同源语义。
#
# ⚠ 2026-10-14 补全 (实测漏检): 原 11 条模式在真实 BEAST2 UCLN+skyline 产物上
# **漏掉**了 meanRate / coefficientOfVariation / skyline / skyline.popSize* ——
# GCVA 实测 `merged/parameter_summary.csv` 里 coefficientOfVariation 的 ESS=4.4
# **低于**当时上报的 worst_ess=4.5, 即门禁**低估了最差 ESS**。已按真实产物补全。
KEY_PARAM_PATTERNS = (
    "posterior",                # 后验 (全局混合)
    "prior",                    # 先验
    "likelihood",               # 似然 (含 treeLikelihood)
    "ucld.mean",                # 分子钟速率 (UCLN, 定年核心)
    "clock.rate",               # 严格钟速率
    "clockrate",
    "meanrate",                 # 平均速率 (BEAST1 输出名)
    "ucld.stdev",               # 松弛钟偏离度 (UCLN)
    "coefficientofvariation",   # 松弛钟 CV (UCLN 的另一写法)
    "trees.treelength",         # 树长 (与时间尺度耦合)
    "treelength",
    "treemodel.rootheight",     # 根高 / TMRCA 相关
    "tmrca",
    "coalescent",               # 树先验 (BEAST2 通用名)
    "skyline",                  # 树先验 (Bayesian Skyline, 含 skyline.popSizeN)
    "demomodel",                # 树先验 (BEAST1 coalescent)
    "popsize",                  # 树先验参数 (skyline.popSizeN / popSize.t)
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
        # ⚠ 措辞纪律 (2026-10-14, 2026-10-14 二次修订): **只陈述实测事实**。
        # 两条约束:
        #   ① **不要**写死 "MCMC 未收敛" —— ESS 不达标有两个成因, 本门禁
        #      **无法区分**, 断言其一即为过度归因:
        #        · 链本身未收敛 (真实问题);
        #        · 采样间隔 logEvery 过大 → 样本数太少 → ESS 上限被样本数压死
        #          (GCVA 实测: time/beast 那批 XML 是 logEvery=50000, 5M 链每链仅
        #           101 个样本, 合并 273 个 → ESS 必然低, 与"是否收敛"无关)。
        #   ② **不要**加"不可采信 / 仅供探索"这类价值判断 (用户 2026-10-14 裁定):
        #      本字段是**事实记录** (数值 + 判定), 不是对结果的自我否决。
        #      用者据数值自行判断; 判据本身 (verdict / worst_ess) 已足够。
        out["reason"] = (f"关键参数最小 ESS={worst_ess:.1f} < {thr:g} "
                         f"(最差: {worst_param}) — ESS 不达标: MCMC 混合不足"
                         f"或采样间隔过大 (logEvery) 使样本数过少; "
                         f"请先复核 logEvery 与每链样本数再判是否真未收敛")
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
