#!/usr/bin/env python3
"""
temporal_signal.py — Date-Randomization Test (DRT) 时间信号检验
==============================================================
实现 Duchêne et al. (2015) 的 Date-Randomization Test:
  1. 对真实日期数据跑 TreeTime-RTT，得到真实 R² 和 substitution rate
  2. 将日期随机打乱 N 次，每次重新跑 TreeTime-RTT
  3. 检查真实 R² 是否显著高于随机分布

参考文献:
  Duchêne S, et al. (2015) The performance of the date-randomization test
  in phylogenetic analyses of time-structured virus data.
  Mol Biol Evol 32:1895-1906.

  PMMoV (2026): DRT via BEAST, 20 randomizations, 95% CI non-overlap
  PVS (Frontiers 2025): DRT with randomized dates × 20

两种模式:
  - fast: 用 TreeTime-RTT 作为代理 (秒级)
  - beast: 生成 N 个 BEAST XML 文件用于完整 DRT (需要手动运行 BEAST)
"""

import os, sys, csv, logging, random
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from collections import defaultdict

import numpy as np
import pandas as pd
from scipy import stats

# ── Path setup ──
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from utils.virphy_bridge import LogCollector, run_treetime_rtt


def date_randomization_test(
    fasta_file: str,
    dates_csv: str,
    output_dir: str,
    n_randomizations: int = 20,
    mode: str = "fast",
    threads: int = 4,
    tree_file: Optional[str] = None,
    rng_seed: Optional[int] = 42,
    log: Optional[LogCollector] = None,
) -> Dict:
    """
    Date-Randomization Test.

    Parameters
    ----------
    fasta_file : str — MAFFT 比对 FASTA
    dates_csv : str — CSV (name, date [, location])
    output_dir : str
    n_randomizations : int — 随机打乱次数 (建议 20)
    mode : str — "fast" (TreeTime proxy) 或 "beast" (生成 XML 文件)
    threads : int
    tree_file : str, optional — 已有系统发育树 (推荐传入 IQ-TREE 树, 避免重新建 NJ 树)
    rng_seed : int or None, default 42
        TreeTime 内部 RNG 固定种子，真实数据与每一次随机化都用同一个值。
        2026-08-27 修复: 不给种子时 TreeTime 的根序列抽样每次不同，
        同一输入实测 R² 在 0.0001~0.0667 之间跳，会把抽样噪声当成日期信号。
        固定同一种子后，随机化分布只反映日期打乱效应，DRT 才可比。
    log : LogCollector

    Returns
    -------
    dict: {
        success, passed (bool),
        real_r2, real_beta, real_clock_rate,
        randomized_r2: [float], randomized_beta: [float],
        randomized_clock_rate: [float],
        r2_percentile: float (真实 R² 在随机分布中的百分位),
        rate_overlap: bool (真实 clock rate 是否落在随机化 rate 的中央 95% 区间内;
            True = 与零分布不可区分; 有效随机化 rate < 5 时为 None),
        randomized_rate_ci95: [lo, hi] (随机化 rate 的 2.5%/97.5% 分位),
        conclusion: str,
        summary_csv: str,
        beast_xml_dir: str (仅 beast 模式)
    }
    """
    log = log or LogCollector()
    os.makedirs(output_dir, exist_ok=True)
    result = {
        "success": False, "passed": False,
        "real_r2": None, "real_beta": None, "real_clock_rate": None,
        "randomized_r2": [], "randomized_beta": [],
        "r2_percentile": None, "rate_overlap": None,
        "conclusion": "", "summary_csv": None, "error": None
    }

    try:
        # ── Load dates ──
        # 2026-08-30 审计修复: 旧版硬编码"第 0 列 name、第 1 列 date"，列序不符时静默把错列
        # 当日期。改为列名优先，位置 fallback；同时过滤无日期行 (参考序列等)，避免空值
        # 混进随机化池改变 dated/undated tip 集合 (DRT 应只在有日期 tip 间置换)。
        dates_df = pd.read_csv(dates_csv)
        name_col = next((c for c in dates_df.columns
                         if str(c).strip().lower() in ('name', 'id', 'accession', 'seq', 'taxon')),
                        dates_df.columns[0])
        date_col = next((c for c in dates_df.columns
                         if str(c).strip().lower() in ('date', 'collectiondate', 'collection_date', 'sampling_date')),
                        dates_df.columns[1] if len(dates_df.columns) > 1 else dates_df.columns[-1])
        real_dates = {str(n).strip(): d for n, d in zip(dates_df[name_col], dates_df[date_col])
                      if pd.notna(d) and str(d).strip() != ''}
        if not real_dates:
            result["error"] = f"日期列 {date_col} 全空或未解析到任何日期"
            return result
        log.emit(f"  Dates: {len(real_dates)} dated tips (col={name_col}/{date_col})")

        # ── Step 1: Real data RTT ──
        log.emit("─" * 50)
        log.emit("DRT Step 1: Real data TreeTime-RTT")
        rtt_dir = os.path.join(output_dir, "real")

        # Need tree — use provided tree, or build quick NJ from alignment
        tree_for_drt = None
        iqtree_path = os.path.join(output_dir, "_drt_ref.treefile")
        if tree_file and os.path.exists(tree_file):
            tree_for_drt = tree_file
            log.emit(f"  Using provided tree: {tree_file}")
        elif not os.path.exists(iqtree_path):
            log.emit("  Building NJ reference tree for DRT...")
            try:
                from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
                from Bio import AlignIO
                aln = AlignIO.read(fasta_file, "fasta")
                calculator = DistanceCalculator('identity')
                dm = calculator.get_distance(aln)
                constructor = DistanceTreeConstructor()
                tree = constructor.nj(dm)
                from Bio import Phylo
                Phylo.write(tree, iqtree_path, "newick")
                tree_for_drt = iqtree_path
                log.emit(f"  NJ tree: {iqtree_path}")
            except Exception as e:
                log.warning(f"  NJ tree failed: {e}, will use TreeTime internal")
        else:
            tree_for_drt = iqtree_path
            log.emit(f"  Using cached NJ tree: {iqtree_path}")

        real_rtt = run_treetime_rtt(
            fasta_file=fasta_file,
            tree_file=tree_for_drt,
            dates_file=dates_csv,
            output_dir=rtt_dir,
            rng_seed=rng_seed,
            log=log,
        )

        result["real_r2"] = real_rtt.get("r_squared")
        result["real_beta"] = real_rtt.get("beta")
        # 2026-09-15 (审查 P2-8): real_clock_rate 原先从不赋值 (恒 None)。
        # run_treetime_rtt 现已把 clock_rate 作为正式字段返回。
        result["real_clock_rate"] = real_rtt.get("clock_rate")
        log.emit(f"  Real data: R²={result['real_r2']}, β={result['real_beta']}, "
                 f"rate={result['real_clock_rate']}")

        # ── Step 2: Randomized data ──
        r2_values = []
        beta_values = []
        rate_values = []

        tip_names = list(real_dates.keys())
        date_values = list(real_dates.values())

        # 日期置换专用的独立 RNG
        # 2026-08-27 修复: 旧代码用全局 random.shuffle, 未固定种子 -> 同一输入
        # 两轮 DRT 的随机化日期分配不同 (实测 round1 rand=[0.0516, 7.7e-05, 0.0188,
        # 0.00097, 0.0173] vs round2 rand=[0.0190, 0.00035, 8.3e-05, 0.0143, 0.00038])。
        # 注释里写的「保证随机序列可复现」在实际执行上并未成立。
        # 现改用 default_rng(rng_seed) 顺序抽取, 置换序列与 TreeTime 种子同源可追溯。
        perm_rng = np.random.default_rng(rng_seed)

        def _permuted_dates():
            perm = perm_rng.permutation(len(date_values))
            return {tip_names[k]: date_values[perm[k]] for k in range(len(tip_names))}

        if mode == "fast":
            log.emit(f"\nDRT Step 2: {n_randomizations} randomizations (TreeTime proxy, "
                     f"{min(threads, n_randomizations)} workers)")
            # 预生成全部随机化 CSV (串行, 置换序列由 perm_rng 决定, 可复现)
            rand_jobs = []
            for i in range(n_randomizations):
                rand_dates = _permuted_dates()
                rand_csv = os.path.join(output_dir, f"randomized_{i+1:02d}.csv")
                pd.DataFrame({
                    name_col: list(rand_dates.keys()),
                    date_col: list(rand_dates.values())
                }).to_csv(rand_csv, index=False)
                rand_jobs.append((i + 1, rand_csv))

            def _run_one(job):
                idx, rand_csv = job
                rand_dir = os.path.join(output_dir, f"randomized_{idx:02d}")
                try:
                    rr = run_treetime_rtt(
                        fasta_file=fasta_file,
                        tree_file=tree_for_drt,
                        dates_file=rand_csv,
                        output_dir=rand_dir,
                        rng_seed=rng_seed,
                        log=None,  # silent
                    )
                except KeyboardInterrupt:
                    raise
                except BaseException as e:
                    # 2026-08-27 修复: 原为 except Exception。
                    # 这是跑在 ThreadPoolExecutor worker 里的, 而 TreeTime.run()
                    # 内部失败时抛 SystemExit(2) (BaseException)。只捕 Exception
                    # 会让 SystemExit 被 future.set_exception 存下, 主线程
                    # result() 时重新抛出, 直接终止整个 DRT / 管线进程。
                    # 实测症状: "静默死亡"、round 日志整体丢失。
                    log.warning(f"  [{idx}] TreeTime 失败: {type(e).__name__}: {e}")
                    rr = {"r_squared": None, "beta": None, "clock_rate": None}
                return idx, rr.get("r_squared"), rr.get("beta"), rr.get("clock_rate")

            from concurrent.futures import ThreadPoolExecutor
            n_workers = max(1, min(threads, n_randomizations))
            if n_workers == 1:
                results_seq = [_run_one(j) for j in rand_jobs]
            else:
                with ThreadPoolExecutor(max_workers=n_workers) as ex:
                    results_seq = list(ex.map(_run_one, rand_jobs))
            for idx, r2_val, beta_val, rate_val in results_seq:
                r2_values.append(r2_val)
                beta_values.append(beta_val)
                # 2026-09-15 (审查 P2-8): 一并收集随机化的 clock rate, 供 rate_overlap
                rate_values.append(rate_val)
                log.emit(f"  [{idx:2d}/{n_randomizations}] R²={r2_val}, β={beta_val}")

        elif mode == "beast":
            log.emit(f"\nDRT Step 2: Generating {n_randomizations} BEAST XML files")
            beast_dir = os.path.join(output_dir, "beast_xmls")
            os.makedirs(beast_dir, exist_ok=True)
            from utils.beast_bridge import generate_beast2_xml

            for i in range(n_randomizations):
                rand_dates = _permuted_dates()
                rand_csv = os.path.join(beast_dir, f"randomized_{i+1:02d}.csv")
                pd.DataFrame({
                    name_col: list(rand_dates.keys()),
                    date_col: list(rand_dates.values())
                }).to_csv(rand_csv, index=False)

                rand_dir = os.path.join(beast_dir, f"randomized_{i+1:02d}")
                try:
                    gen = generate_beast2_xml(
                        fasta_file=fasta_file,
                        dates_csv=rand_csv,
                        output_dir=rand_dir,
                        log=None,
                    )
                    log.emit(f"  [{i+1:2d}/{n_randomizations}] BEAST XML: {gen.get('xml_path')}")
                except Exception as e:
                    log.warning(f"  [{i+1:2d}/{n_randomizations}] Failed: {e}")

            result["beast_xml_dir"] = beast_dir
            log.emit(f"\n  BEAST XML files in: {beast_dir}")
            log.emit(f"  Run each with: beast -threads 4 randomized_XX/beast2.xml")

        result["randomized_r2"] = r2_values
        result["randomized_beta"] = beta_values
        result["randomized_clock_rate"] = rate_values
        result["rng_seed"] = rng_seed
        result["deterministic"] = rng_seed is not None

        # ── Step 3: Statistical comparison ──
        valid_r2 = [v for v in r2_values if v is not None] if r2_values else []
        # 2026-08-27 新增: 显式记录零分布的实际样本数。
        # 背景: TreeTime 失败的那次随机化会被记为 None 并被 valid_r2 静默过滤,
        # 零分布实际样本数小于名义 n_randomizations。实测历史产物里 P4 丢 3/10、
        # G 丢 2/10、M 丢 1/10, 而结论里从未体现这件事。
        result["n_randomizations"] = n_randomizations
        result["n_valid_randomizations"] = len(valid_r2)
        if len(valid_r2) < n_randomizations:
            log.warning(f"  {n_randomizations - len(valid_r2)}/{n_randomizations} "
                        f"次随机化失败 (TreeTime 报错), 零分布实际只有 {len(valid_r2)} 个点")
        real_r2 = result.get("real_r2")

        # ── Step 3b: clock rate 是否可与零分布区分 (审查 P2-8) ────────────
        # ⚠⚠ 2026-10-13 重排: rate_overlap 必须在判定 passed 之前算出。
        # 原实现把 passed 定义在 rate_overlap 之前, 判定只用了 r2_percentile,
        # 导致 GCVA 出现「success:false + passed:true + rate_overlap:false」
        # 三者并存的自相矛盾产物 (真实 rate 1.589e-04 落在随机化 95% CI
        # [4.11e-06, 1.318e-04] 之外, 却仍判 "Sufficient temporal signal")。
        # 定义 (与 R² 百分位同构, 不依赖 TreeTime 未暴露的 rate CI):
        #   把 n 次随机化得到的 clock rate 当零分布, 取其中央 95% 区间
        #   [2.5%, 97.5%]; 真实 rate 落在这个区间内 → rate_overlap=True,
        #   即"真实速率与随机化速率不可区分"。
        # 需要 >=5 个有效随机化 rate, 否则写 None (样本不足, 不给结论)。
        _valid_rates = [float(v) for v in rate_values if v is not None]
        result["n_valid_randomized_rates"] = len(_valid_rates)
        rate_overlap = None
        if result.get("real_clock_rate") is not None and len(_valid_rates) >= 5:
            _lo = float(np.percentile(_valid_rates, 2.5))
            _hi = float(np.percentile(_valid_rates, 97.5))
            _real_rate = float(result["real_clock_rate"])
            result["randomized_rate_ci95"] = [_lo, _hi]
            rate_overlap = bool(_lo <= _real_rate <= _hi)
            result["rate_overlap"] = rate_overlap
            log.emit(f"  Clock rate: real={_real_rate:.3e}, "
                     f"randomized 95% CI=[{_lo:.3e}, {_hi:.3e}] "
                     f"→ overlap={rate_overlap}")
        else:
            log.warning(f"  rate_overlap 未计算 (real_clock_rate="
                        f"{result.get('real_clock_rate')}, "
                        f"有效随机化 rate 数={len(_valid_rates)}, 需 >=5)")

        # ── Step 3c: 综合判定 (2026-10-13 修复) ──────────────────────────
        # DRT 判「有时间信号」须同时满足两个条件:
        #   ① r2_percentile >= 95   —— 真实日期的 R² 显著优于随机日期
        #   ② rate_overlap is True  —— 真实速率与随机化速率不可区分
        # 原实现只用 ①。问题: ① 是**分布内相对位置**, 只看"真实是否比随机好",
        # **不看 R² 绝对值**。一个零信号数据集 (绝对 R²≈0.05) 只要随机化后掉得更低,
        # 同样能拿 100 百分位 → 误判 PASSED (GCVA 实测即为此例)。
        # ② 补上了"速率可识别性"这一必要条件。两者缺一不可 (conjunction)。
        # 若 ② 无法计算 (有效随机化 rate < 5) → 保守取 False, 不给"通过"结论。
        if valid_r2 and real_r2 is not None:
            percentile = stats.percentileofscore(valid_r2, real_r2)
            result["r2_percentile"] = float(percentile)

            # 2026-08-27 修复: 必须转成 Python bool。
            # percentileofscore 返回 numpy.float64, 故 `>=` 得到 np.bool_;
            # 而 np.bool_ 不是 Python bool 的子类, 会被写 drt_results.json 时的
            # isinstance(v, (str, int, float, bool, ...)) 过滤器静默剔除。
            # 实测 8/8 个历史 drt_results.json 里都没有 passed 字段。
            r2_ok = bool(percentile >= 95.0)
            rate_ok = bool(rate_overlap) if rate_overlap is not None else False
            passed = bool(r2_ok and rate_ok)
            result["passed"] = passed
            result["r2_criterion_met"] = r2_ok
            result["rate_criterion_met"] = rate_ok

            # 交叉核查: 旧实现 (仅凭 R²) 的判定结果。两者不一致时显式告警,
            # 避免"静默改变历史结论"。
            result["passed_r2_only"] = r2_ok
            if r2_ok and not rate_ok:
                log.warning(
                    f"  ⚠ DRT 判定被推翻: R² 百分位 {percentile:.0f}% ≥95 (旧判据 PASSED), "
                    f"但 rate_overlap={rate_overlap} → 新判据 FAILED。"
                    f"真实速率落在随机化速率分布之外, 时钟不可识别。")

            if passed:
                result["conclusion"] = (
                    f"✓ PASSED: Real R² ({real_r2:.4f}) exceeds "
                    f"{percentile:.0f}% of {len(valid_r2)} randomized values "
                    f"(threshold: 95%), AND clock rate is indistinguishable from "
                    f"the randomized-date null (rate_overlap=True). "
                    f"Sufficient temporal signal for molecular dating."
                )
            elif r2_ok and not rate_ok:
                result["conclusion"] = (
                    f"✗ FAILED: Real R² ({real_r2:.4f}) reaches the "
                    f"{percentile:.0f}th percentile (≥95%), BUT the clock rate "
                    f"does not overlap the randomized-date null "
                    f"(rate_overlap={rate_overlap}). Rate is unidentified → "
                    f"temporal signal insufficient for molecular dating."
                )
            else:
                result["conclusion"] = (
                    f"✗ FAILED: Real R² ({real_r2:.4f}) is only at "
                    f"{percentile:.0f}th percentile of {len(valid_r2)} randomized "
                    f"values (r2_criterion_met={r2_ok}, "
                    f"rate_criterion_met={rate_ok}). "
                    f"Insufficient temporal signal — molecular clock may "
                    f"be unreliable."
                )
            log.emit(f"\n  {result['conclusion']}")
        else:
            result["conclusion"] = "DRT inconclusive: insufficient data"
            real_r2 = None
            result.setdefault("passed", False)

        # ── Step 4: Save summary ──
        summary_csv = os.path.join(output_dir, "drt_summary.csv")
        rows = [{
            "dataset": "real", "R2": result["real_r2"], "beta": result["real_beta"]
        }]
        for i, (r2, beta) in enumerate(zip(r2_values, beta_values)):
            rows.append({
                "dataset": f"randomized_{i+1:02d}", "R2": r2, "beta": beta
            })
        pd.DataFrame(rows).to_csv(summary_csv, index=False)
        result["summary_csv"] = summary_csv

        # Save structured JSON for cross-virus comparison
        drt_json = os.path.join(output_dir, "drt_results.json")
        import json as _json
        _json_safe = {k: v for k, v in result.items()
                      if isinstance(v, (str, int, float, bool, list, type(None)))}
        with open(drt_json, 'w') as _jf:
            _json.dump(_json_safe, _jf, indent=2, default=str)
        result["drt_json"] = drt_json

        # ── Step 5: Visualization ──
        if valid_r2:
            try:
                import matplotlib
                matplotlib.use('Agg')
                import matplotlib.pyplot as plt
                # 2026-10-12: 配色统一 OKABE（原 #E9522E / #00A0E9 / #b0b0b0
                # 均不在项目调色板内）；dpi 150 → 300；补 y 轴网格。
                try:
                    from utils.sci_palette import (
                        apply_sci_style, OKABE_ITO, despine)
                except ImportError:
                    from sci_palette import (
                        apply_sci_style, OKABE_ITO, despine)
                apply_sci_style(base_size=11)

                fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

                # R² distribution
                ax1.hist(valid_r2, bins=max(10, n_randomizations // 2),
                         color=OKABE_ITO[8], edgecolor='#666666', linewidth=0.4,
                         alpha=0.75, label='Randomized dates')
                ax1.axvline(real_r2, color=OKABE_ITO[5], linewidth=2.4, linestyle='--',
                            label=f'Observed ($R^2$ = {real_r2:.4f})')
                ax1.axvline(np.percentile(valid_r2, 95), color=OKABE_ITO[4], linewidth=1.6,
                            linestyle=':', label='95th percentile')
                ax1.set_xlabel('$R^2$', fontsize=13)
                ax1.set_ylabel('Frequency', fontsize=13)
                # 2026-10-12（审稿人视角）：标题写明随机化次数 —— DRT 的检验强度
                # 完全取决于置换次数，读者必须能一眼看到 n。
                ax1.set_title(f'Date-randomization test: $R^2$'
                              f'  (n = {n_randomizations} randomizations)',
                              fontsize=13)
                ax1.grid(alpha=0.3, axis='y')
                ax1.legend(fontsize=10)
                despine(ax1)

                # Beta distribution
                valid_beta = [v for v in beta_values if v is not None]
                if valid_beta and result["real_beta"] is not None:
                    ax2.hist(valid_beta, bins=max(10, n_randomizations // 2),
                             color=OKABE_ITO[8], edgecolor='#666666', linewidth=0.4,
                             alpha=0.75, label='Randomized dates')
                    ax2.axvline(result["real_beta"], color=OKABE_ITO[5], linewidth=2.4,
                                linestyle='--',
                                label=f'Observed ($\\beta$ = {result["real_beta"]:.2e})')
                    ax2.set_xlabel('$\\beta$  (subs/site/yr)', fontsize=13)
                    ax2.set_ylabel('Frequency', fontsize=13)
                    ax2.set_title(f'Date-randomization test: rate'
                                  f'  (n = {n_randomizations} randomizations)',
                                  fontsize=13)
                    # 2026-10-12（审稿人视角）：β ~ 1e-4 量级，旧版 7 个刻度全部
                    # 写成 0.000000/0.000250/… 挤成一串无法辨认。
                    # 改用科学计数法刻度（轴外统一 ×10⁻⁴ 偏移量），刻度只留 4 个。
                    from matplotlib.ticker import MaxNLocator, ScalarFormatter
                    ax2.xaxis.set_major_locator(MaxNLocator(nbins=4))
                    _sf = ScalarFormatter(useMathText=True)
                    _sf.set_powerlimits((0, 0))
                    ax2.xaxis.set_major_formatter(_sf)
                    ax2.ticklabel_format(axis='x', style='sci', scilimits=(0, 0))
                    ax2.grid(alpha=0.3, axis='y')
                    ax2.legend(fontsize=10)
                    despine(ax2)

                drt_plot = os.path.join(output_dir, "drt_results.pdf")
                fig.savefig(drt_plot, bbox_inches='tight', dpi=300, format='pdf')
                plt.close(fig)
                result["drt_plot"] = drt_plot
                log.emit(f"  DRT plot: {drt_plot}")
            except Exception as e:
                log.warning(f"  DRT plot failed: {e}")

        result["success"] = True
        return result

    except Exception as e:
        result["error"] = str(e)
        log.error(f"DRT failed: {e}")
        return result


def check_temporal_signal_simple(
    fasta_file: str,
    dates_csv: str,
    log: Optional[LogCollector] = None,
) -> Dict:
    """
    快速时间信号检查：仅跑一次 TreeTime-RTT，根据 R² 给出建议。

    阈值 (基于 37_Virology_Bioinformatics 课程 + BETS 论文建议):
      R² > 0.85  → 严格分子钟可用
      R² > 0.70  → 放松分子钟 (UCLN)
      R² > 0.10  → 边际信号，必须跑 DRT 确认
      R² > 0.05  → 弱信号，BEAST 可能不收敛，推荐 BETS 检验
      R² > 0.01  → 几乎无信号，不推荐定年
      R² < 0.01  → 无时间信号

    此外，负 β 值 (R² 可能虚高) 提示 TreeTime 自动 rooting 问题，
    应检查 TempEst 图确认根的方向正确。
    """
    log = log or LogCollector()
    import tempfile, shutil
    tmpdir = tempfile.mkdtemp(prefix="temporal_check_")

    try:
        # Ensure we have a tree — build quick NJ if none provided
        tree_for_check = None
        nj_tree_path = os.path.join(tmpdir, "_quick_nj.nwk")
        try:
            from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
            from Bio import AlignIO
            aln = AlignIO.read(fasta_file, "fasta")
            calculator = DistanceCalculator('identity')
            dm = calculator.get_distance(aln)
            constructor = DistanceTreeConstructor()
            tree = constructor.nj(dm)
            from Bio import Phylo
            Phylo.write(tree, nj_tree_path, "newick")
            tree_for_check = nj_tree_path
        except Exception:
            pass

        rtt = run_treetime_rtt(
            fasta_file=fasta_file,
            tree_file=tree_for_check,
            dates_file=dates_csv,
            output_dir=tmpdir,
            log=None,
        )
        r2 = rtt.get("r_squared")
        beta = rtt.get("beta")

        if r2 is None:
            return {"has_signal": False, "level": "unknown",
                    "R2": None, "beta": beta,
                    "recommendation": "Cannot compute R². Check data quality or provide a tree file."}

        if r2 > 0.85:
            level, rec = "strong_strict", "Strict clock viable. BEAST with strict clock recommended."
            clock_rec = "strict"
        elif r2 > 0.70:
            level, rec = "strong", "UCLN relaxed clock recommended. Proceed with BEAST."
            clock_rec = "ucln"
        elif r2 > 0.10:
            level, rec = "moderate", "Marginal signal. Run full DRT before BEAST. Use UCLN clock."
            clock_rec = "ucln"
        elif r2 > 0.05:
            level, rec = "weak", "Weak signal. Run BETS (Bayesian Evaluation of Temporal Signal) to confirm. " \
                      "BEAST may not converge under UCLN."
            clock_rec = "ucln"
        elif r2 > 0.01:
            level, rec = "marginal", "Almost no temporal signal. Bayesian dating not recommended. " \
                      "If you must, use strict clock + exponential popSize prior (BETS recommendation)."
            clock_rec = "strict"
        else:
            level, rec = "none", "No temporal signal. Do NOT attempt molecular dating."
            clock_rec = None

        # Check for negative beta (suspicious rooting)
        beta = rtt.get("beta")
        if beta is not None and beta < 0:
            rec += " ⚠ Negative β detected — possible rooting artifact. Check TempEst root direction."

        return {
            "has_signal": r2 > 0.01,
            "level": level,
            "R2": r2,
            "beta": beta,
            "recommendation": rec,
            "recommended_clock": clock_rec,
        }
    finally:
        import shutil
        shutil.rmtree(tmpdir, ignore_errors=True)
