#!/usr/bin/env python3
"""
beast_parser.py — BEAST MCMC 结果解析
======================================
解析 BEAST 1.x/2.x 输出文件，提取参数估计、MCC 树、迁移矩阵。

功能:
  1. parse_beast_log()      — 解析 .log → 参数均值/HPD/ESS 表
  2. parse_beast_trees()    — 解析 .trees → MCC 树 + 节点统计
  3. extract_migration_bf() — 解析 BSSVS indicators → BF 迁移矩阵
  4. build_summary_report() — 生成完整分析报告
  5. plot_beast_diagnostics() — MCMC 收敛诊断图 (trace + density)
"""

import csv, json, os, re, sys, math
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from collections import defaultdict
from datetime import datetime

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from utils.virphy_bridge import LogCollector
# 2026-09-15 (审查 P2-9): BF 阈值与"发散截断值"全管线唯一定义, 见 virphy_bridge
from utils.virphy_bridge import (DEFAULT_BF_THRESHOLD, BF_INFINITE_CAP,
                                 bayes_factor, find_indicator_columns)
# 历史坑: calculate_95hpd 曾只在 _calc_ess 内局部 import, 模块级三处裸用 →
# NameError → parse_beast_log 返回 error, key_params 空, TMRCA 静默丢失。
from utils.ess import calculate_ess, calculate_95hpd


# ═══════════════════════════════════════════════════════════════════
# 1. BEAST Log Parser
# ═══════════════════════════════════════════════════════════════════

def parse_beast_log(
    log_file: str,
    burnin_pct: float = 10.0,
    ess_threshold: float = 200.0,
    log: Optional[LogCollector] = None,
) -> Dict:
    """
    解析 BEAST .log 文件，提取所有参数的统计量。

    Returns
    -------
    dict: {
        success, n_samples, n_burnin, burnin_pct,
        parameters: [{name, mean, median, stdev, hpd_lower, hpd_upper, ess, converged}],
        summary_csv, summary_json
    }
    """
    log = log or LogCollector()
    result = {"success": False, "n_samples": 0, "n_burnin": 0,
              "parameters": [], "converged_count": 0, "error": None}

    if not os.path.exists(log_file):
        result["error"] = f"Log not found: {log_file}"
        return result

    try:
        df = pd.read_csv(log_file, sep="\t", comment="#")
        n_total = len(df)
        if n_total == 0:
            result["error"] = "Empty log"
            return result

        burnin_idx = int(n_total * burnin_pct / 100.0)
        df_burnin = df.iloc[burnin_idx:].reset_index(drop=True)
        n_used = len(df_burnin)

        result["n_samples"] = n_total
        result["n_burnin"] = n_used
        result["burnin_pct"] = burnin_pct

        skip_cols = {"Sample", "sample", "state", "State"}
        params = []

        for col in df_burnin.select_dtypes(include=[np.number]).columns:
            if col in skip_cols:
                continue

            vals = df_burnin[col].dropna().values
            if len(vals) < 10:
                continue

            mean_val = float(np.mean(vals))
            median_val = float(np.median(vals))
            stdev_val = float(np.std(vals))
            hpd_lower, hpd_upper = calculate_95hpd(vals)
            hpd_lower = float(hpd_lower)
            hpd_upper = float(hpd_upper)
            ess_val = _calc_ess(vals)
            converged = ess_val >= ess_threshold

            params.append({
                "name": col,
                "mean": round(mean_val, 6),
                "median": round(median_val, 6),
                "stdev": round(stdev_val, 6),
                "hpd_95_lower": round(hpd_lower, 6),
                "hpd_95_upper": round(hpd_upper, 6),
                "ess": round(ess_val, 1),
                "converged": converged,
            })

        result["parameters"] = sorted(params, key=lambda x: x["ess"])
        result["converged_count"] = sum(1 for p in params if p["converged"])

        # Key parameters
        # 历史坑: 子串匹配 ('prior' in 'posterior', 'likelihood' in 'treeLikelihood')
        # 导致 key_params 项被 ESS 更高的同缀列覆盖, 结果依赖迭代顺序。改精确匹配 + 别名。
        key_names = ["clock.rate", "clockRate", "ucld.mean", "ucld.stdev",
                     "kappa", "constant.popSize", "treeModel.rootHeight",
                     "posterior", "prior", "likelihood",
                     "coefficientOfVariation", "meanRate"]
        _alias = {"likelihood": ["likelihood", "treeLikelihood",
                                  "generalizedSkyLineLikelihood",
                                  "ancestralTreeLikelihood"]}
        result["key_params"] = {}
        for p in params:
            lname = p["name"].lower()
            for kn in key_names:
                kls = kn.lower()
                cands = _alias.get(kn, [kls])
                if lname == kls or (kn in _alias and lname in cands):
                    result["key_params"].setdefault(kn, p)

        # Save
        csv_path = os.path.join(os.path.dirname(log_file),
                                "parameter_summary.csv")
        pd.DataFrame(params).to_csv(csv_path, index=False)
        result["summary_csv"] = csv_path

        log.emit(f"Log parsed: {n_total} samples, {burnin_pct}% burn-in → {n_used}")
        log.emit(f"  {result['converged_count']}/{len(params)} parameters converged (ESS≥{ess_threshold})")

        # Show key parameters
        for kn, kp in result["key_params"].items():
            ess_flag = "✓" if kp["converged"] else "⚠"
            log.emit(f"  {ess_flag} {kn}: {kp['mean']:.4e} [{kp['hpd_95_lower']:.4e}–{kp['hpd_95_upper']:.4e}] ESS={kp['ess']:.0f}")

        result["success"] = True
        return result

    except Exception as e:
        result["error"] = str(e)
        log.error(f"Log parsing failed: {e}")
        return result


def _calc_ess(x: np.ndarray) -> float:
    """Effective Sample Size (Tracer 同算法, Geyer initial positive sequence).

    与 **Tracer** (beast-mcmc TraceCorrelation.analyseCorrelationNumeric) 一致;
    修正旧版逐 lag 截断低估问题。

    ⚠️ 与 YR-MPE mcmc_utils.calculate_ESS **有意不一致** (2026-09-16 上游核对):
       上游 P_k 循环 k 从 0 起, 使 P_0 = rho_0 + rho_1 把 rho_0 ≡ 1 计入,
       而 tau_hat = -1 + 2*sum(P_k) 的 -1 已代表 rho_0 → rho_0 被双计。
       实测 (20 轮 x 6 档, 对 Tracer 的中位偏差): 我方 0.000–0.090%,
       上游 0.22–1.76% (最坏 35.8%)。正确配对是 (rho_1,rho_2),(rho_3,rho_4),...
       即 tau = 1 + 2*sum(P_k), k 从 1 起 (Tracer 的 i 从 1 步进 2)。
       如曾按上游口径复核过 ESS 阈值, 那些结论需重跑。
       证据: MMPV-RNA/_consistency_check_20260916/ess_stats_robust.py
    """
    from utils.ess import calculate_ess, calculate_95hpd
    return calculate_ess(x)


# ═══════════════════════════════════════════════════════════════════
# 2. BEAST Trees Parser
# ═══════════════════════════════════════════════════════════════════

def parse_beast_trees(
    trees_file: str,
    output_dir: str,
    burnin_pct: float = 10.0,
    log: Optional[LogCollector] = None,
) -> Dict:
    """
    解析 BEAST .trees 文件，提取树统计量。

    Returns
    -------
    dict: {success, n_trees, n_burnin, tree_heights: [float],
           mcc_tree_stats: {mean_height, hpd_lower, hpd_upper}}
    """
    log = log or LogCollector()
    result = {"success": False, "n_trees": 0, "n_burnin": 0,
              "tree_heights": [], "mcc_tree_stats": {}, "error": None}

    if not os.path.exists(trees_file):
        result["error"] = f"Trees file not found: {trees_file}"
        return result

    try:
        heights = []
        n_trees = 0

        with open(trees_file, 'r') as f:
            for line in f:
                if not line.startswith('tree'):
                    continue
                n_trees += 1
                # 取树高: 最后一个 height= (根节点; tip 排在最前 height≈0)。
                # 历史坑: re.search 取第一个 → 恒为 tip 高度 → TMRCA 静默归零。
                heights_all = re.findall(r'height=([\d.]+)', line)
                if heights_all:
                    heights.append(float(heights_all[-1]))

        if n_trees == 0:
            result["error"] = "No trees found"
            return result

        burnin_idx = int(n_trees * burnin_pct / 100.0)
        heights_burnin = heights[burnin_idx:] if heights else []
        n_burnin = len(heights_burnin)

        result["n_trees"] = n_trees
        result["n_burnin"] = n_burnin
        result["tree_heights"] = heights

        if heights_burnin:
            _hl, _hh = calculate_95hpd(heights_burnin)
            result["mcc_tree_stats"] = {
                "mean_height": round(np.mean(heights_burnin), 4),
                "median_height": round(np.median(heights_burnin), 4),
                "hpd_lower": round(_hl, 4),
                "hpd_upper": round(_hh, 4),
            }

        log.emit(f"Trees parsed: {n_trees} trees, {burnin_pct}% burn-in → {n_burnin}")
        if result["mcc_tree_stats"]:
            s = result["mcc_tree_stats"]
            log.emit(f"  Root height: {s['mean_height']:.4f} [{s['hpd_lower']:.4f}–{s['hpd_upper']:.4f}]")

        result["success"] = True
        return result

    except Exception as e:
        result["error"] = str(e)
        log.error(f"Trees parsing failed: {e}")
        return result


# ═══════════════════════════════════════════════════════════════════
# 3. Migration BF Matrix
# ═══════════════════════════════════════════════════════════════════

def extract_migration_matrix(
    log_file: str,
    locations: List[str],
    output_dir: str,
    burnin_pct: float = 10.0,
    bf_threshold: float = DEFAULT_BF_THRESHOLD,
    log: Optional[LogCollector] = None,
) -> Dict:
    """
    从 BEAST log 的 BSSVS indicators 提取迁移 BF 矩阵。

    Returns
    -------
    dict: {success, matrix: [[bf]], significant_routes: [{from, to, bf}],
           migration_csv, migration_json}
    """
    log = log or LogCollector()
    result = {"success": False, "matrix": [], "significant_routes": [],
              "migration_csv": None, "error": None}

    if not os.path.exists(log_file):
        result["error"] = f"Log not found: {log_file}"
        return result

    try:
        df = pd.read_csv(log_file, sep="\t", comment="#")
        burnin_idx = int(len(df) * burnin_pct / 100.0)
        df = df.iloc[burnin_idx:]

        n_loc = len(locations)
        if n_loc < 2:
            result["error"] = "Need ≥2 locations"
            return result

        # Find indicator columns (BEAST 1.x: Location.indicators, BEAST 2.x: location.indicators)
        # 列检测口径全管线统一 (virphy_bridge.find_indicator_columns, P2-27 补)
        ind_cols = find_indicator_columns(df.columns)
        rate_cols = [c for c in df.columns
                     if 'rates' in c.lower() and 'indicator' not in c.lower()
                     and 'location' in c.lower()]

        if not ind_cols and not rate_cols:
            log.warning("No BSSVS indicator columns found — no phylogeography in this log")
            result["error"] = "No phylogeography indicators"
            return result

        # Build matrix
        matrix = [[{"bf": 0.0, "posterior": 0.0, "rate": 0.0}
                   for _ in range(n_loc)] for _ in range(n_loc)]

        routes = []

        # ── 列 → (from_i, to_i) 映射 ────────────────────────────────────
        # 2026-09-15 (审查 P2-27): 旧版纯位置索引 (idx // (n_loc-1) ...), 依赖列在
        # DataFrame 里的**出现顺序**。logger 用 sort="smart" (字母序) 时列序是
        # location.indicators1, ...10, ...11, ...2 —— 位置映射整体错位 (n_loc=5 时
        # 实测 19/20 列错), 迁移路线张冠李戴且完全静默。
        #
        # 2026-09-15 (P2-27 补): 仅"命名优先"还不够 —— 本管线 beast1_bridge 写的是
        # **单个向量参数** `location.indicators` (dimension=n_rates), 因此 BEAST1 的
        # log 列名同样是 `location.indicators<k>` (1-based 数字), 命名正则根本命不中,
        # 仍会落到顺序回退上。故补齐三级解析, 逐级降级并显式记录:
        #   ① location.indicators.<from>.<to> —— 按地点名 (与列序无关)
        #   ② location.indicators<k>          —— 按 1-based 数字后缀, 与 XML 声明序
        #                                        一致 (与 phylogeo_bridge.idx_map 同口径)
        #   ③ enumerate 出现顺序              —— 最后兜底, 且**显式告警**
        idx_by_name = {loc: i for i, loc in enumerate(locations)}

        def _resolve_pair(col_name):
            # re.IGNORECASE: BEAST1 的列名前缀是 `Location.indicators` (大写 L),
            # BEAST2 是 `location.indicators` —— 两种都要按名解析 (地点名查表仍精确匹配)。
            m = re.match(r'location\.indicators\.(\w+)\.(\w+)', col_name,
                         re.IGNORECASE)
            if m and m.group(1) in idx_by_name and m.group(2) in idx_by_name:
                return idx_by_name[m.group(1)], idx_by_name[m.group(2)]
            return None

        def _pos_pair(idx):
            """位置映射 (BEAST 参数声明序, row-major 跳对角)"""
            from_i = idx // (n_loc - 1)
            to_i = idx % (n_loc - 1)
            if to_i >= from_i:
                to_i += 1
            return from_i, to_i

        def _numeric_pair(col_name):
            """location.indicators<k> (1-based) → 声明序第 k 对 (row-major 跳对角)。

            与 phylogeo_bridge.extract_migration_bf 的 idx_map 完全同口径:
            idx_map[k-1] = (locations[i], locations[j]), i≠j 按行优先展开。
            """
            m = re.search(r'(\d+)\s*$', col_name)
            if not m:
                return None
            k = int(m.group(1)) - 1
            if not (0 <= k < n_loc * (n_loc - 1)):
                return None
            return _pos_pair(k)

        _by_name, _by_num, _by_pos = [], [], []
        for _i, _c in enumerate(ind_cols):
            _r = _resolve_pair(_c)
            if _r is not None:
                _by_name.append((_c, _r))
                continue
            _r = _numeric_pair(_c)
            if _r is not None:
                _by_num.append((_c, _r))
                continue
            _by_pos.append((_c, _pos_pair(_i)))

        if not _by_pos:
            pairs = _by_name + _by_num
            result["index_source"] = ("column_name" if not _by_num
                                      else "numeric_suffix")
            log.emit(f"  BSSVS 列按{'列名' if not _by_num else '数字后缀'}解析 "
                     f"({len(pairs)} 对 from/to, 与列出现顺序无关)")
        else:
            pairs = _by_name + _by_num + _by_pos
            _unres = [c for c, _ in _by_pos]
            result["index_source"] = "positional"
            result["unnamed_indicator_columns"] = _unres
            log.warning(
                f"  BSSVS 列名无法解析 ({len(_unres)}/{len(ind_cols)} 个): "
                f"{_unres[:4]} → 退回**出现顺序**索引; "
                f"若 BEAST logger 用了 sort=\"smart\" (字母序) 路线会错位")

        # rate 列 → (from_i, to_i): 同样先按名, 再按数字后缀 (与列序无关)
        rate_by_pair = {}
        _rate_by_pos = []
        for _i, c in enumerate(rate_cols):
            m = re.match(r'location\.rates\.(\w+)\.(\w+)', c, re.IGNORECASE)
            if m and m.group(1) in idx_by_name and m.group(2) in idx_by_name:
                rate_by_pair[(idx_by_name[m.group(1)],
                              idx_by_name[m.group(2)])] = c
                continue
            m2 = re.search(r'(\d+)\s*$', c)
            if m2:
                _k2 = int(m2.group(1)) - 1
                if 0 <= _k2 < n_loc * (n_loc - 1):
                    rate_by_pair[_pos_pair(_k2)] = c
                    continue
            _rate_by_pos.append((_i, c))
        for _i, c in _rate_by_pos:
            rate_by_pair.setdefault(_pos_pair(_i), c)

        # 未取整 BF, 供"显著"判据使用 (取整后再比阈值会在边界翻转:
        # 4.999 -> 5.0 会被误判为显著)
        bf_raw = {}

        for _k, (col, (from_i, to_i)) in enumerate(pairs):
            if from_i >= n_loc or to_i >= n_loc or from_i == to_i:
                log.warning(f"  跳过无法映射的 BSSVS 列: {col}")
                continue

            posterior = float(df[col].mean())

            _rcol = rate_by_pair.get((from_i, to_i))
            if _rcol is None and _k < len(rate_cols) and len(rate_cols) == len(ind_cols):
                _rcol = rate_cols[_k]
            mean_rate = float(df[_rcol].mean()) if _rcol is not None else 0.0

            # Bayes Factor (全管线唯一定义, 见 virphy_bridge.bayes_factor)
            bf = bayes_factor(posterior)
            bf_raw[(from_i, to_i)] = bf

            matrix[from_i][to_i] = {
                "bf": round(bf, 2),
                "posterior": round(posterior, 4),
                "rate": round(mean_rate, 6),
            }

            if bf >= bf_threshold:
                routes.append({
                    "from": locations[from_i], "to": locations[to_i],
                    "from_idx": from_i, "to_idx": to_i,
                    "bf": round(bf, 2),
                    "posterior": round(posterior, 4),
                    "rate": round(mean_rate, 6),
                })

        result["matrix"] = matrix
        result["significant_routes"] = sorted(routes, key=lambda x: -x["bf"])

        # Save
        # 2026-09-15 (P2-27 补): 输出目录不存在时旧版 open() 抛 FileNotFoundError,
        # 被外层 except 吞成 success=False + 一个与真实原因无关的报错。自建目录。
        os.makedirs(output_dir, exist_ok=True)
        csv_path = os.path.join(output_dir, "migration_bf.csv")
        with open(csv_path, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(["From", "To", "BF", "Posterior", "Rate", "Significant"])
            for from_i in range(n_loc):
                for to_i in range(n_loc):
                    if from_i != to_i:
                        m = matrix[from_i][to_i]
                        writer.writerow([
                            locations[from_i], locations[to_i],
                            m["bf"], m["posterior"], m["rate"],
                            "Yes" if bf_raw.get((from_i, to_i), 0.0) >= bf_threshold else ""
                        ])
        result["migration_csv"] = csv_path

        log.emit(f"Migration BF matrix: {n_loc}×{n_loc}")
        log.emit(f"  {len(routes)} significant routes (BF≥{bf_threshold})")
        for r in routes[:10]:
            log.emit(f"  {r['from']} → {r['to']}: BF={r['bf']:.1f}")

        result["success"] = True
        return result

    except Exception as e:
        result["error"] = str(e)
        log.error(f"Migration matrix failed: {e}")
        return result


# ═══════════════════════════════════════════════════════════════════
# 4. Summary Report Builder
# ═══════════════════════════════════════════════════════════════════

def build_beast_summary_report(
    log_file: str,
    trees_file: str,
    locations: List[str],
    output_dir: str,
    metadata: Optional[Dict] = None,
    burnin_pct: float = 10.0,
    bf_threshold: float = DEFAULT_BF_THRESHOLD,
    log: Optional[LogCollector] = None,
) -> Dict:
    """
    生成完整的 BEAST 分析总结报告。

    Returns
    -------
    dict: {success, report_path, summary: {clock_rate, tmrca, pop_size, ...}}
    """
    log = log or LogCollector()
    os.makedirs(output_dir, exist_ok=True)
    result = {"success": False, "report_path": None, "summary": {}, "error": None}

    try:
        # Parse log
        log_result = parse_beast_log(log_file, burnin_pct=burnin_pct, log=log)
        # Parse trees
        trees_result = parse_beast_trees(trees_file, output_dir, burnin_pct=burnin_pct, log=log)
        # Extract migrations
        mig_result = extract_migration_matrix(log_file, locations, output_dir,
                                              burnin_pct=burnin_pct,
                                              bf_threshold=bf_threshold, log=log)

        # Extract key summaries
        key = log_result.get("key_params", {})
        summary = {}

        # Clock rate
        for k in ["clock.rate", "clockRate"]:
            if k in key:
                summary["clock_rate"] = key[k]["mean"]
                summary["clock_rate_hpd"] = (key[k]["hpd_95_lower"], key[k]["hpd_95_upper"])
                break

        # Population size
        for k in ["constant.popSize", "popSize"]:
            if k in key:
                summary["pop_size"] = key[k]["mean"]
                break

        # TMRCA (root height)
        for k in ["treeModel.rootHeight", "TreeHeight"]:
            if k in key:
                tmrca = key[k]["mean"]
                summary["tmrca"] = tmrca
                summary["tmrca_hpd"] = (key[k]["hpd_95_lower"], key[k]["hpd_95_upper"])
                break

        # UCLD statistics
        if "ucld.mean" in key:
            summary["ucld_mean"] = key["ucld.mean"]["mean"]
        if "ucld.stdev" in key:
            summary["ucld_stdev"] = key["ucld.stdev"]["mean"]
        if "coefficientOfVariation" in key:
            summary["cov"] = key["coefficientOfVariation"]["mean"]

        # Convergence
        summary["n_params"] = len(log_result.get("parameters", []))
        summary["n_converged"] = log_result.get("converged_count", 0)
        summary["convergence_pct"] = round(
            100 * summary["n_converged"] / max(summary["n_params"], 1), 1)

        # Migration
        summary["n_locations"] = len(locations)
        summary["n_significant_routes"] = len(mig_result.get("significant_routes", []))
        sig_routes = mig_result.get("significant_routes", [])
        summary["top_migration"] = (
            f"{sig_routes[0]['from']} → {sig_routes[0]['to']} (BF={sig_routes[0]['bf']:.1f})"
            if sig_routes else "None"
        )

        # Metadata
        if metadata:
            summary.update(metadata)

        result["summary"] = summary

        # Save summary JSON
        json_path = os.path.join(output_dir, "beast_summary.json")
        with open(json_path, 'w') as f:
            json.dump(summary, f, indent=2, ensure_ascii=False,
                     default=lambda x: str(x) if isinstance(x, (np.floating, np.integer)) else x)

        # Print summary
        log.emit("\n" + "=" * 55)
        log.emit("  BEAST Analysis Summary")
        log.emit("=" * 55)
        if summary.get("clock_rate"):
            log.emit(f"  Clock rate:     {summary['clock_rate']:.2e} subs/site/yr")
        if summary.get("tmrca"):
            log.emit(f"  TMRCA:          {summary['tmrca']:.2f}")
        if summary.get("pop_size"):
            log.emit(f"  Pop size (θ):   {summary['pop_size']:.2f}")
        if summary.get("cov"):
            log.emit(f"  Coeff of var:   {summary['cov']:.3f}")
        log.emit(f"  Convergence:    {summary['n_converged']}/{summary['n_params']} "
                f"({summary['convergence_pct']}%)")
        log.emit(f"  Sig. routes:    {summary['n_significant_routes']}")
        if sig_routes:
            log.emit(f"  Top route:      {summary['top_migration']}")

        # Build HTML report
        report_path = _build_html_report(summary, log_result, mig_result,
                                         locations, output_dir)
        result["report_path"] = report_path
        result["summary_json"] = json_path
        result["success"] = True
        return result

    except Exception as e:
        result["error"] = str(e)
        log.error(f"Report failed: {e}")
        return result


def _build_html_report(
    summary: Dict,
    log_result: Dict,
    mig_result: Dict,
    locations: List[str],
    output_dir: str,
) -> str:
    """Build HTML report for BEAST analysis."""

    params = log_result.get("parameters", [])
    sig_routes = mig_result.get("significant_routes", [])

    param_rows = ""
    for p in params[:30]:
        flag = "✅" if p["converged"] else "⚠️"
        param_rows += (
            f"<tr>"
            f"<td>{flag}</td>"
            f"<td>{p['name']}</td>"
            f"<td>{p['mean']:.4e}</td>"
            f"<td>[{p['hpd_95_lower']:.4e}–{p['hpd_95_upper']:.4e}]</td>"
            f"<td>{p['ess']:.0f}</td>"
            f"</tr>\n"
        )

    route_rows = ""
    for r in sig_routes[:15]:
        route_rows += (
            f"<tr><td>{r['from']}</td><td>{r['to']}</td>"
            f"<td><b>{r['bf']:.1f}</b></td><td>{r['posterior']:.3f}</td></tr>\n"
        )

    report = f'''<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>BEAST Analysis Report</title>
<style>
  * {{ margin:0; padding:0; box-sizing:border-box; }}
  body {{ font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',sans-serif; background:#f5f6fa; color:#2c3e50; }}
  .container {{ max-width:1100px; margin:0 auto; padding:30px; }}
  h1 {{ font-size:24px; margin-bottom:5px; }}
  h2 {{ font-size:18px; margin:25px 0 12px; color:#34495e; border-bottom:2px solid #3498db; padding-bottom:5px; }}
  .meta {{ color:#888; font-size:14px; margin-bottom:20px; }}
  .summary {{ display:grid; grid-template-columns:repeat(auto-fit,minmax(220px,1fr)); gap:15px; margin-bottom:25px; }}
  .card {{ background:white; padding:18px; border-radius:8px; box-shadow:0 1px 4px rgba(0,0,0,0.08); }}
  .card .value {{ font-size:24px; font-weight:700; color:#3498db; }}
  .card .label {{ font-size:12px; color:#888; margin-top:5px; }}
  table {{ width:100%; border-collapse:collapse; background:white; border-radius:8px; overflow:hidden; box-shadow:0 1px 4px rgba(0,0,0,0.08); margin-bottom:20px; }}
  th {{ background:#34495e; color:white; padding:10px 12px; text-align:left; font-size:13px; }}
  td {{ padding:8px 12px; border-bottom:1px solid #ecf0f1; font-size:13px; }}
  tr:hover {{ background:#f8f9fa; }}
  .warn {{ color:#e67e22; }}
  .footer {{ text-align:center; color:#bbb; font-size:12px; margin-top:40px; }}
</style>
</head>
<body>
<div class="container">
<h1>📊 BEAST Bayesian Analysis Report</h1>
<p class="meta">Generated {datetime.now().strftime("%Y-%m-%d %H:%M")}</p>

<div class="summary">
  <div class="card"><div class="value">{summary.get('clock_rate', '?') if summary.get('clock_rate') else 'N/A'}</div><div class="label">Clock Rate (subs/site/yr)</div></div>
  <div class="card"><div class="value">{summary.get('tmrca', '?')}</div><div class="label">TMRCA (years)</div></div>
  <div class="card"><div class="value">{summary.get('pop_size', '?')}</div><div class="label">Pop Size (θ)</div></div>
  <div class="card"><div class="value">{summary.get('n_converged',0)}/{summary.get('n_params',0)}</div><div class="label">Converged Parameters</div></div>
  <div class="card"><div class="value">{summary.get('n_locations',0)}</div><div class="label">Locations</div></div>
  <div class="card"><div class="value">{summary.get('n_significant_routes',0)}</div><div class="label">Significant Routes</div></div>
</div>

<h2>📈 Parameter Estimates</h2>
<table>
<tr><th></th><th>Parameter</th><th>Mean</th><th>95% HPD</th><th>ESS</th></tr>
{param_rows}
</table>

<h2>✈ Migration Routes (BF≥5)</h2>
<table>
<tr><th>From</th><th>To</th><th>Bayes Factor</th><th>Posterior</th></tr>
{route_rows}
</table>

<div class="footer">
  MMPV-RNA virome_phylo_pipeline · BEAST v1.10
</div>
</div>
</body>
</html>'''

    report_path = os.path.join(output_dir, "beast_analysis_report.html")
    with open(report_path, 'w', encoding='utf-8') as f:
        f.write(report)
    return report_path


# ═══════════════════════════════════════════════════════════════════
# 5. MCMC Diagnostics Plot
# ═══════════════════════════════════════════════════════════════════

def plot_beast_diagnostics(
    log_file: str,
    output_dir: str,
    burnin_pct: float = 10.0,
    params_to_plot: Optional[List[str]] = None,
    log: Optional[LogCollector] = None,
) -> str:
    """
    生成 BEAST MCMC 收敛诊断图 (trace + density)。

    Returns path to PDF.
    """
    log = log or LogCollector()

    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        df = pd.read_csv(log_file, sep="\t", comment="#")
        burnin_idx = int(len(df) * burnin_pct / 100.0)

        # Select parameters to plot
        skip = {"Sample", "sample", "state", "State"}
        numeric_cols = [c for c in df.select_dtypes(include=[np.number]).columns
                       if c not in skip]

        if params_to_plot:
            plot_cols = [c for c in params_to_plot if c in numeric_cols]
        else:
            # Auto-select key parameters
            key_patterns = ["posterior", "prior", "likelihood", "clock",
                          "kappa", "popSize", "treeModel.rootHeight",
                          "ucld.mean", "ucld.stdev"]
            plot_cols = []
            for kp in key_patterns:
                for c in numeric_cols:
                    if kp.lower() in c.lower() and c not in plot_cols:
                        plot_cols.append(c)
            if len(plot_cols) > 8:
                plot_cols = plot_cols[:8]

        if not plot_cols:
            plot_cols = numeric_cols[:6]

        n_rows = len(plot_cols)
        fig, axes = plt.subplots(n_rows, 2, figsize=(14, 2.5 * n_rows))
        if n_rows == 1:
            axes = [axes]

        for i, col in enumerate(plot_cols):
            vals = df[col].values

            # Trace plot
            ax_trace = axes[i][0]
            ax_trace.plot(vals, color='#3498db', alpha=0.6, linewidth=0.5)
            ax_trace.axvline(burnin_idx, color='#e74c3c', linestyle='--',
                            alpha=0.7, label=f'Burn-in ({burnin_pct}%)')
            ax_trace.set_ylabel(col, fontsize=10)
            ax_trace.legend(fontsize=8, loc='upper right')

            # Density plot
            ax_dens = axes[i][1]
            post_burnin = vals[burnin_idx:]
            ax_dens.hist(post_burnin, bins=50, color='#2ecc71', alpha=0.7,
                        density=True, edgecolor='white')
            ax_dens.axvline(np.mean(post_burnin), color='#e74c3c',
                           linestyle='--', linewidth=2, label='Mean')
            _vl, _vh = calculate_95hpd(post_burnin)
            ax_dens.axvline(_vl, color='#95a5a6',
                           linestyle=':', linewidth=1)
            ax_dens.axvline(_vh, color='#95a5a6', linestyle=':', linewidth=1, label='95% HPD')
            ax_dens.set_ylabel('Density', fontsize=10)
            ax_dens.legend(fontsize=8)

            if i == 0:
                ax_trace.set_title('MCMC Trace', fontsize=12)
                ax_dens.set_title('Posterior Density', fontsize=12)

        plt.tight_layout()
        diag_path = os.path.join(output_dir, "mcmc_diagnostics.pdf")
        fig.savefig(diag_path, bbox_inches='tight', dpi=150, format='pdf')
        plt.close(fig)

        log.emit(f"MCMC diagnostics: {diag_path} ({n_rows} parameters)")
        return diag_path

    except Exception as e:
        log.warning(f"Diagnostics plot failed: {e}")
        return ""


# ═══════════════════════════════════════════════════════════════════
# 6. Live MCMC Monitor
# ═══════════════════════════════════════════════════════════════════

# monitor_beast_run 在 max_checks<=0 (调用方未设上限) 时的安全轮数上限,
# 防止 check_interval>0 时前台无限挂住 (2026-09-15)
MAX_UNBOUNDED_CHECKS = 240


def monitor_beast_run(
    log_file: str,
    ess_threshold: float = 200.0,
    check_interval: float = 60.0,
    max_checks: int = 0,
    log: Optional[LogCollector] = None,
) -> Dict:
    """
    监控正在运行的 BEAST MCMC，估算当前 ESS 并判断是否可提前停止。
    """
    log = log or LogCollector()
    result = {"success": False, "n_samples": 0, "can_stop": False,
              "parameters": [], "converged_count": 0, "error": None}

    if not os.path.exists(log_file):
        result["error"] = f"Log file not found: {log_file}"
        return result

    import time as _time
    checks_done = 0
    prev_n_samples = 0

    # 2026-09-15 修复: max_checks=0 配 check_interval>0 会构成**无限轮询** ——
    # 调用方不显式传 max_checks 时前台进程永不返回。
    # 语义澄清: max_checks<=0 仍表示"不设人为轮数上限", 但必须给出终止条件:
    #   · check_interval<=0 -> 退化为单次快照 (上限 1 轮)
    #   · check_interval>0  -> 采用安全上限 MAX_UNBOUNDED_CHECKS 轮
    if max_checks <= 0:
        if check_interval <= 0:
            max_checks = 1
        else:
            max_checks = MAX_UNBOUNDED_CHECKS
            log.warning(
                f"monitor_beast_run: max_checks<=0 且 check_interval={check_interval}s "
                f"→ 采用安全上限 {max_checks} 轮 (约 {max_checks * check_interval / 3600:.1f} h); "
                f"需要更久请显式传 max_checks")

    while checks_done < max_checks:
        try:
            df = pd.read_csv(log_file, sep="\t", comment="#")
            n_samples = len(df)
            if n_samples == 0 or n_samples == prev_n_samples:
                if check_interval > 0:
                    _time.sleep(check_interval); checks_done += 1; continue
                break
            prev_n_samples = n_samples

            half_idx = n_samples // 2
            df_recent = df.iloc[half_idx:]
            skip_cols = {"Sample", "sample", "state", "State"}
            params = []
            converged = 0

            for col in df_recent.select_dtypes(include=[np.number]).columns[:20]:
                if col in skip_cols:
                    continue
                vals = df_recent[col].dropna().values
                if len(vals) < 10:
                    continue
                ess_val = _calc_ess(vals)
                is_ok = ess_val >= ess_threshold
                if is_ok:
                    converged += 1
                params.append({"name": col, "ess": round(ess_val, 1), "converged": is_ok})

            params.sort(key=lambda x: x["ess"])
            result["n_samples"] = n_samples
            result["parameters"] = params
            result["converged_count"] = converged

            worst = params[:5]
            log.emit(f"[{checks_done+1}] {n_samples} samples, {converged}/{len(params)} ESS≥{ess_threshold}")
            for p in worst:
                flag = "✓" if p["converged"] else "⚠"
                log.emit(f"  {flag} {p['name'][:40]}: ESS={p['ess']:.0f}")

            # 2026-09-15 修复: 旧版把 success=True 放在收敛 break **之后** →
            # 越收敛 success 越是 False (初值 False), 调用方按 success 判断会得到
            # 完全相反的结论 (最收敛的一次被判为"监控失败")。
            result["success"] = True

            if len(params) > 0 and converged >= len(params) * 0.9:
                result["can_stop"] = True
                result["recommendation"] = f"✓ {converged}/{len(params)} converged. Safe to stop MCMC."
                log.emit(f"\n  {result['recommendation']}")
                break

            if check_interval == 0:
                break
            _time.sleep(check_interval)
            checks_done += 1
        except Exception as e:
            result["error"] = str(e)
            return result

    return result
