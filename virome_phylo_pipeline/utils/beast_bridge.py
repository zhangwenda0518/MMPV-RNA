#!/usr/bin/env python3
"""
beast_bridge.py — BEAST2 XML 生成、执行、后处理
================================================
整合 gen_beast_xml.py / gen_beast2.py / auto_beast1.py / prep_beauti.py /
prep_ages.py / run_beast.py / post_beast.py / post_mcc.py 的功能。
"""

import csv
import os
import re
import subprocess
import sys
import time
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
import pandas as pd
from Bio import SeqIO

# ── 将 virphy_bridge.py 的 LogCollector 拿来复用 ──
from utils.virphy_bridge import LogCollector


# ═══════════════════════════════════════════════════════════════════
# BEAST2 XML 生成
# ═══════════════════════════════════════════════════════════════════

def _to_decimal_year(date_str: str, sample: str = "") -> float:
    """将日期字符串转为小数年。支持 YYYY, YYYY-MM, YYYY-MM-DD, YYYY/MM/DD。

    2026-09-15 (审查 P2-26): 实现搬到 utils/decimal_year.py, 全管线统一。
    另修: 旧版解析失败**静默返回 2020.0** —— 一个脏日期会被当成"2020 年采样",
    悄悄污染时间轴与定年结果。现在改为 fail-loud (调用方已有 try/except 剔除逻辑)。
    """
    from utils.decimal_year import to_decimal_year
    return to_decimal_year(date_str, sample=sample, strict=True)


def generate_beast2_xml(
    fasta_file: str,
    dates_csv: str,
    output_dir: str,
    chain_length: int = 5_000_000,
    log_every: int = 5000,
    clock_rate_init: float = 0.001,
    kappa_init: float = 2.0,
    log: Optional[LogCollector] = None,
) -> Dict:
    """
    从 FASTA 多序列比对 + 日期 CSV 生成 BEAST2 XML。

    Parameters
    ----------
    fasta_file : str — MAFFT 比对后的 FASTA
    dates_csv : str — CSV，列: name, date
    output_dir : str
    chain_length : int — MCMC 链长 (默认 5M)
    log_every : int — 每 N 代记录一次
    clock_rate_init : float — 时钟速率初始值
    kappa_init : float — HKY kappa 初始值

    Returns
    -------
    dict: {success, xml_path, n_taxa, seq_len, min_year, max_year}
    """
    log = log or LogCollector()
    os.makedirs(output_dir, exist_ok=True)
    result = {"success": False, "xml_path": None, "n_taxa": 0,
              "seq_len": 0, "min_year": None, "max_year": None, "error": None}

    try:
        # ── 1. Load alignment ──
        aln = list(SeqIO.parse(fasta_file, "fasta"))
        if len(aln) < 3:
            result["error"] = f"Need >=3 sequences, got {len(aln)}"
            return result
        seq_len = len(str(aln[0].seq))
        log.emit(f"Loaded {len(aln)} sequences, {seq_len} bp")

        # ── 2. Load dates ──
        dates = {}
        with open(dates_csv, encoding='utf-8-sig') as f:
            reader = csv.DictReader(f)
            for row in reader:
                name = row.get('name', '').strip()
                date = row.get('date', '').strip()
                if name and date:
                    dates[name] = date

        # Filter to dated sequences
        dated = [(r, dates[r.id]) for r in aln if r.id in dates]
        # 脏日期防御 (2026-09-02): 解析失败剔除并警告, 不再崩溃 (与 beast1_bridge 同款)
        clean = []
        for r, d in dated:
            try:
                _to_decimal_year(d)
                clean.append((r, d))
            except ValueError as e:
                log.warning(f"  ⚠️ 剔除样本 {r.id}: 脏日期 ({e})")
        dropped = len(dated) - len(clean)
        if dropped:
            log.warning(f"  ⚠️ 共剔除 {dropped} 个脏日期样本, 继续 {len(clean)} 条")
        dated = clean
        if len(dated) < 3:
            result["error"] = f"Only {len(dated)} sequences have dates (need >=3)"
            return result

        year_vals = [_to_decimal_year(d) for _, d in dated]
        max_year = max(year_vals)
        min_year = min(year_vals)
        result["n_taxa"] = len(dated)
        result["seq_len"] = seq_len
        result["min_year"] = min_year
        result["max_year"] = max_year
        log.emit(f"{len(dated)} dated sequences, year range: {min_year:.1f}–{max_year:.1f}")

        # ── 3. Generate BEAST2 XML ──
        xml_path = os.path.join(output_dir, "beast2.xml")
        output_prefix = os.path.join(output_dir, "beast2")

        with open(xml_path, 'w', encoding='utf-8') as f:
            f.write('<?xml version="1.0" encoding="UTF-8" standalone="no"?>\n')
            f.write('<beast version="2.7"\n')
            f.write('  namespace="beast.core:beast.evolution.alignment:beast.evolution.tree.coalescent:'
                    'beast.core.util:beast.evolution.operators:beast.evolution.sitemodel:'
                    'beast.evolution.substitutionmodel:beast.evolution.likelihood'
                    ':beast.base.evolution.alignment:beast.base.evolution.tree.coalescent'
                    ':beast.base.inference:beast.base.inference.operator'
                    ':beast.base.evolution.operator:beast.base.evolution.sitemodel'
                    ':beast.base.evolution.substitutionmodel:beast.base.evolution.likelihood'
                    ':beast.base.inference.distribution">\n\n')

            # Alignment
            f.write('  <data id="aln" spec="Alignment" name="alignment">\n')
            for rec, _ in dated:
                f.write(f'    <sequence taxon="{rec.id}" value="{str(rec.seq)}" />\n')
            f.write('  </data>\n\n')

            # Tip dates trait — convert to "years before most recent"
            f.write('  <trait id="dateTrait" spec="beast.base.evolution.tree.TraitSet" '
                    'traitname="date-backward" value="\n')
            date_parts = []
            for rec, d in dated:
                age = max_year - _to_decimal_year(d)
                date_parts.append(f'    {rec.id}={age:.4f}')
            f.write(',\n'.join(date_parts))
            f.write('\n  ">\n')
            f.write('    <taxa id="taxa" spec="TaxonSet" alignment="@aln"/>\n')
            f.write('  </trait>\n\n')

            # Starting tree (simple ladder)
            tree_str = '(' + ','.join([r.id for r, _ in dated]) + ')'

            # MCMC run
            f.write(f'  <run id="mcmc" spec="MCMC" chainLength="{chain_length}" storeEvery="{log_every}">\n')
            f.write('    <state id="state" storeEvery="5000">\n')
            f.write('      <tree id="Tree.t:aln" spec="beast.base.evolution.tree.TreeParser" '
                    'taxa="@aln" IsLabelledNewick="true"\n')
            f.write(f'            newick="{tree_str}" />\n')
            f.write(f'      <parameter id="clockRate.c:aln" spec="RealParameter" '
                    f'lower="1.0E-9" upper="0.1" name="stateNode">{clock_rate_init}</parameter>\n')
            f.write(f'      <parameter id="popSize.t:aln" spec="RealParameter" '
                    f'lower="0.0" upper="1.0E100" name="stateNode">0.3</parameter>\n')
            f.write(f'      <parameter id="kappa.s:aln" spec="RealParameter" '
                    f'lower="0.0" name="stateNode">{kappa_init}</parameter>\n')
            # Base frequencies — HKY needs 4 values
            f.write('      <parameter id="freqParameter.s:aln" spec="RealParameter" '
                    'dimension="4" lower="0.0" upper="1.0" name="stateNode">'
                    '0.25 0.25 0.25 0.25</parameter>\n')
            f.write('    </state>\n\n')

            # Random tree initializer
            f.write('    <init id="RandomTree.t:aln" spec="beast.base.evolution.tree.RandomTree" '
                    'estimate="false" initial="@Tree.t:aln" taxa="@aln">\n')
            f.write('      <populationModel id="ConstantPopulationInit.t:aln" '
                    'spec="beast.evolution.tree.coalescent.ConstantPopulation">\n')
            f.write('        <parameter id="initPopSize.t:aln" spec="RealParameter" '
                    'name="popSize">0.3</parameter>\n')
            f.write('      </populationModel>\n')
            f.write('    </init>\n\n')

            # Posterior distribution
            f.write('    <distribution id="posterior" spec="CompoundDistribution">\n')
            # Prior
            f.write('      <distribution id="prior" spec="CompoundDistribution">\n')
            f.write('        <distribution id="Coalescent" spec="Coalescent">\n')
            f.write('          <populationModel id="ConstantPopulation" spec="ConstantPopulation" '
                    'populationSize="@popSize.t:aln"/>\n')
            f.write('          <treeIntervals id="TreeIntervals" spec="TreeIntervals" '
                    'tree="@Tree.t:aln"/>\n')
            f.write('        </distribution>\n')
            # PopSize prior
            f.write('        <distribution id="popSizePrior" spec="OneOnXPrior">\n')
            f.write('          <parameter idref="popSize.t:aln"/>\n')
            f.write('        </distribution>\n')
            # Kappa prior
            f.write('        <prior id="KappaPrior.s:aln" name="distribution">\n')
            f.write('          <parameter idref="kappa.s:aln"/>\n')
            f.write('          <LogNormal name="distr" M="1.0" S="1.25" offset="0.0"/>\n')
            f.write('        </prior>\n')
            # Clock rate prior
            f.write('        <distribution id="clockPrior" spec="LogNormalDistributionModel" '
                    f'M="{clock_rate_init}" S="1.25" offset="0.0">\n')
            f.write(f'          <parameter spec="RealParameter" estimate="false" name="M">'
                    f'{clock_rate_init}</parameter>\n')
            f.write('          <parameter spec="RealParameter" estimate="false" name="S">'
                    '1.25</parameter>\n')
            f.write('          <parameter idref="clockRate.c:aln"/>\n')
            f.write('        </distribution>\n')
            f.write('      </distribution>\n')
            # Likelihood
            f.write('      <distribution id="likelihood" spec="CompoundDistribution">\n')
            f.write('        <distribution id="treeLikelihood" spec="TreeLikelihood" '
                    'data="@aln" tree="@Tree.t:aln">\n')
            f.write('          <siteModel id="SiteModel" spec="SiteModel">\n')
            f.write('            <parameter spec="RealParameter" estimate="false" '
                    'name="mutationRate">1.0</parameter>\n')
            f.write('            <parameter id="gammaShape" spec="RealParameter" '
                    'estimate="false" name="shape">1.0</parameter>\n')
            f.write('            <parameter id="propInvariant" spec="RealParameter" '
                    'estimate="false" lower="0.0" name="proportionInvariant" '
                    'upper="1.0">0.0</parameter>\n')
            f.write('            <substModel id="hky" spec="HKY" kappa="@kappa.s:aln">\n')
            f.write('              <frequencies id="estimatedFreqs" spec="Frequencies" '
                    'frequencies="@freqParameter.s:aln"/>\n')
            f.write('            </substModel>\n')
            f.write('          </siteModel>\n')
            f.write('          <branchRateModel id="StrictClock.c:aln" spec="StrictClockModel" '
                    'clock.rate="@clockRate.c:aln"/>\n')
            f.write('        </distribution>\n')
            f.write('      </distribution>\n')
            f.write('    </distribution>\n\n')

            # Operators
            ops = [
                ('StrictClockRateScaler.c:aln', 'ScaleOperator',
                 'parameter="@clockRate.c:aln" scaleFactor="0.75" weight="3"'),
                ('strictClockUpDownOperator.c:aln', 'UpDownOperator',
                 'scaleFactor="0.75" weight="3"'),
                ('KappaScaler.s:aln', 'ScaleOperator',
                 'parameter="@kappa.s:aln" scaleFactor="0.75" weight="1"'),
                ('FrequenciesExchanger.s:aln', 'DeltaExchangeOperator',
                 'parameter="@freqParameter.s:aln" delta="0.01" weight="1"'),
                ('CoalescentConstantTreeScaler.t:aln', 'ScaleOperator',
                 'tree="@Tree.t:aln" scaleFactor="0.95" weight="3"'),
                ('CoalescentConstantTreeRootScaler.t:aln', 'ScaleOperator',
                 'tree="@Tree.t:aln" rootOnly="true" scaleFactor="0.95" weight="3"'),
                ('CoalescentConstantUniformOperator.t:aln', 'Uniform',
                 'tree="@Tree.t:aln" weight="30"'),
                ('CoalescentConstantSubtreeSlide.t:aln', 'SubtreeSlide',
                 'tree="@Tree.t:aln" weight="15"'),
                ('CoalescentConstantNarrow.t:aln', 'Exchange',
                 'tree="@Tree.t:aln" isNarrow="true" weight="15"'),
                ('CoalescentConstantWide.t:aln', 'Exchange',
                 'tree="@Tree.t:aln" isNarrow="false" weight="3"'),
                ('CoalescentConstantWilsonBalding.t:aln', 'WilsonBalding',
                 'tree="@Tree.t:aln" weight="3"'),
                ('PopSizeScaler.t:aln', 'ScaleOperator',
                 'parameter="@popSize.t:aln" scaleFactor="0.75" weight="3"'),
            ]

            for op_id, op_spec, op_attrs in ops:
                if op_spec == 'UpDownOperator':
                    f.write(f'    <operator id="{op_id}" spec="{op_spec}" {op_attrs}>\n')
                    f.write('      <up>\n')
                    f.write('        <parameter idref="clockRate.c:aln"/>\n')
                    f.write('      </up>\n')
                    f.write('      <down>\n')
                    f.write('        <tree idref="Tree.t:aln"/>\n')
                    f.write('      </down>\n')
                    f.write('    </operator>\n')
                else:
                    f.write(f'    <operator id="{op_id}" spec="{op_spec}" {op_attrs}/>\n')

            # Loggers
            f.write('\n    <!-- Screen logger -->\n')
            f.write(f'    <logger id="screenlog" spec="Logger" logEvery="{log_every}">\n')
            for ref in ['posterior', 'likelihood', 'prior', 'clockRate.c:aln']:
                f.write(f'      <log idref="{ref}"/>\n')
            f.write('    </logger>\n\n')

            f.write(f'    <logger id="tracelog" spec="Logger" fileName="{output_prefix}.log" '
                    f'logEvery="{log_every}" sort="smart">\n')
            for ref in ['posterior', 'likelihood', 'prior', 'clockRate.c:aln',
                        'popSize.t:aln', 'kappa.s:aln']:
                f.write(f'      <log idref="{ref}"/>\n')
            f.write('      <log id="TreeHeight.t:aln" spec="TreeHeightLogger" '
                    'tree="@Tree.t:aln"/>\n')
            f.write('    </logger>\n\n')

            f.write(f'    <logger id="treelog" spec="Logger" fileName="{output_prefix}.trees" '
                    f'logEvery="{log_every}" mode="tree">\n')
            f.write('      <log id="TreeWithMetaDataLogger" spec="TreeWithMetaDataLogger" '
                    'tree="@Tree.t:aln"/>\n')
            f.write('    </logger>\n')

            f.write('  </run>\n</beast>\n')

        result["xml_path"] = xml_path
        result["success"] = True
        log.emit(f"BEAST2 XML: {xml_path} ({os.path.getsize(xml_path)} bytes)")
        return result

    except Exception as e:
        result["error"] = str(e)
        log.error(f"BEAST2 XML generation failed: {e}")
        return result


# ═══════════════════════════════════════════════════════════════════
# BEAUti 输入准备 (手动模式)
# ═══════════════════════════════════════════════════════════════════

def prepare_beauti_inputs(
    fasta_file: str,
    dates_csv: str,
    output_dir: str,
    max_taxa: int = 50,
    log: Optional[LogCollector] = None,
) -> Dict:
    """
    生成 BEAUti 兼容的简化 FASTA + trait table，用于手动在 Windows BEAUti 中导入。

    如果序列数超过 max_taxa，按时间均匀子抽样。

    Returns
    -------
    dict: {success, fasta_path, traits_path, n_taxa, locations, error}
    """
    import random
    log = log or LogCollector()
    os.makedirs(output_dir, exist_ok=True)
    result = {"success": False, "fasta_path": None, "traits_path": None,
              "n_taxa": 0, "locations": [], "error": None}

    try:
        aln = list(SeqIO.parse(fasta_file, "fasta"))
        dates = {}
        locs = {}
        with open(dates_csv, encoding='utf-8-sig') as f:
            reader = csv.DictReader(f)
            for row in reader:
                name = row.get('name', '').strip()
                d = row.get('date', '').strip()
                loc = row.get('location', '').strip()
                if name and d:
                    dates[name] = d
                if name and loc:
                    locs[name] = loc

        candidates = [(r, dates.get(r.id, '')) for r in aln if r.id in dates]
        candidates.sort(key=lambda x: x[1])

        # Subsampling
        if len(candidates) > max_taxa:
            random.seed(42)
            step = max(1, len(candidates) // max_taxa)
            subset = candidates[::step][:max_taxa]
        else:
            subset = candidates

        n = len(subset)
        log.emit(f"BEAUti prep: {n} taxa (from {len(candidates)} candidates)")

        # Simplify IDs
        def short_id(full: str) -> str:
            return full.split('_')[0] if '_' in full else full

        # Write simplified FASTA
        fasta_path = os.path.join(output_dir, "beauti.fasta")
        with open(fasta_path, 'w') as f:
            for rec, _ in subset:
                f.write(f">{short_id(rec.id)}\n{str(rec.seq)}\n")

        # Write trait table
        traits_path = os.path.join(output_dir, "beauti_traits.txt")
        location_set = set()
        with open(traits_path, 'w') as f:
            headers = ['taxon', 'date']
            has_loc = any(r.id in locs for r, _ in subset)
            if has_loc:
                headers.append('location')
            f.write('\t'.join(headers) + '\n')
            for rec, d in subset:
                sid = short_id(rec.id)
                yr = _to_decimal_year(d)
                parts = [sid, f"{yr:.4f}"]
                if has_loc:
                    raw_loc = locs.get(rec.id, 'Unknown')
                    # Simplify: "China, Ningxia, Yinchuan" → "Ningxia"
                    loc = raw_loc.replace('"', '').strip()
                    if ':' in loc:
                        loc = loc.split(':')[-1].strip()
                    if ',' in loc:
                        loc = [p.strip() for p in loc.split(',')][-1]
                    parts.append(loc)
                    location_set.add(loc)
                f.write('\t'.join(parts) + '\n')

        result["fasta_path"] = fasta_path
        result["traits_path"] = traits_path
        result["n_taxa"] = n
        result["locations"] = sorted(location_set)
        result["success"] = True
        log.emit(f"BEAUti files: fasta={fasta_path}, traits={traits_path}")
        if location_set:
            log.emit(f"  Locations: {location_set}")
        return result

    except Exception as e:
        result["error"] = str(e)
        log.error(f"BEAUti prep failed: {e}")
        return result


# ═══════════════════════════════════════════════════════════════════
# BEAST 执行管理
# ═══════════════════════════════════════════════════════════════════

def run_beast(
    xml_path: str,
    output_dir: str,
    beast_bin: str = "beast",
    threads: int = 8,
    timeout_hours: int = 48,
    log: Optional[LogCollector] = None,
    seed: Optional[int] = None,
) -> Dict:
    """
    运行 BEAST2，实时监控进度。

    Parameters
    ----------
    xml_path : str — BEAST XML 文件路径
    output_dir : str — 工作目录
    beast_bin : str — beast 可执行文件名 (beast 或 beast2)
    threads : int — 线程数
    timeout_hours : int — 超时 (小时)
    log : LogCollector, optional

    Returns
    -------
    dict: {success, log_path, trees_path, exit_code, elapsed_seconds, error}
    """
    log = log or LogCollector()
    os.makedirs(output_dir, exist_ok=True)
    result = {"success": False, "log_path": None, "trees_path": None,
              "exit_code": None, "elapsed_seconds": 0, "error": None}

    if not os.path.exists(xml_path):
        result["error"] = f"XML not found: {xml_path}"
        return result

    # 历史坑: XML 里 overwrite=false, 多链路径带 -overwrite 而单链没带 → --force 重跑被拒
    # 可复现性 (2026-09-15 修复): 旧版不传 -seed, BEAST 默认随机种子 → 同输入两次结果不同。
    # 现默认传确定性种子 (PHYLO_SEED_BASE 环境变量可覆盖), 并写进返回值便于追溯。
    if seed is None:
        try:
            seed = int(os.environ.get("PHYLO_SEED_BASE", "") or 20260915)
        except ValueError:
            seed = 20260915
    cmd = [beast_bin, "-threads", str(threads), "-overwrite",
           "-seed", str(int(seed)), xml_path]
    result["seed"] = int(seed)
    log.emit(f"BEAST: {' '.join(cmd)}")
    log.emit(f"Working dir: {output_dir}")

    t0 = time.time()
    try:
        proc = subprocess.run(
            cmd,
            capture_output=True,
            text=True, encoding='utf-8', errors='replace',
            timeout=timeout_hours * 3600,
            cwd=output_dir,
        )
        elapsed = time.time() - t0
        result["elapsed_seconds"] = elapsed
        result["exit_code"] = proc.returncode

        # Check for SEVERE errors in stderr
        severe_errors = [line for line in proc.stderr.split('\n') if 'SEVERE' in line]
        if severe_errors:
            for err in severe_errors[:5]:
                log.warning(f"BEAST: {err.strip()[:200]}")

        # Look for output files
        for fname in os.listdir(output_dir):
            if fname.endswith('.log') and 'beast' in fname.lower():
                result["log_path"] = os.path.join(output_dir, fname)
            elif fname.endswith('.trees') and 'beast' in fname.lower():
                result["trees_path"] = os.path.join(output_dir, fname)

        if proc.returncode == 0:
            result["success"] = True
            log.emit(f"BEAST completed in {elapsed:.0f}s")
            if result["log_path"]:
                log.emit(f"  Log: {result['log_path']}")
            if result["trees_path"]:
                log.emit(f"  Trees: {result['trees_path']}")
        else:
            result["error"] = f"BEAST exit code {proc.returncode}"
            # Show last stderr lines
            stderr_tail = '\n'.join(proc.stderr.split('\n')[-10:])
            log.error(f"BEAST failed. stderr tail:\n{stderr_tail}")

        return result

    except subprocess.TimeoutExpired:
        elapsed = time.time() - t0
        result["elapsed_seconds"] = elapsed
        result["error"] = f"BEAST timed out after {timeout_hours}h"
        log.error(result["error"])
        return result
    except FileNotFoundError:
        result["error"] = f"BEAST binary not found: {beast_bin}"
        log.error(result["error"])
        return result
    except Exception as e:
        result["error"] = str(e)
        log.error(f"BEAST run failed: {e}")
        return result


# ═══════════════════════════════════════════════════════════════════
# MCMC 收敛诊断 (ESS)
# ═══════════════════════════════════════════════════════════════════

def _effective_sample_size(x: np.ndarray, max_lag: Optional[int] = None) -> float:
    """ESS (Tracer 同算法, Geyer initial positive sequence criterion).

    与 **Tracer** (beast-mcmc TraceCorrelation) 一致。
    旧版 "遇第一个非正自相关即截断" 会低估 ESS, 已修正。

    ⚠️ 与 YR-MPE mcmc_utils.calculate_ESS **有意不一致** (2026-09-16 上游核对):
       上游 P_k 循环 k 从 0 起, 把 rho_0 ≡ 1 双计 (tau_hat = -1 + 2*sum(P_k)
       的 -1 已代表 rho_0)。实测对 Tracer 中位偏差: 我方 0.000–0.090%,
       上游 0.22–1.76% (最坏 35.8%)。详见 utils/ess.py 模块 docstring 与
       MMPV-RNA/_consistency_check_20260916/ess_stats_robust.py
    """
    from utils.ess import calculate_ess
    return calculate_ess(x)


def analyze_beast_ess(
    log_file: str,
    output_dir: str,
    burnin_pct: float = 10.0,
    ess_threshold: float = 200.0,
    log: Optional[LogCollector] = None,
) -> Dict:
    """
    读取 BEAST .log 文件，计算所有参数的 ESS，标记未达标者。

    Parameters
    ----------
    log_file : str — BEAST .log 文件路径
    output_dir : str — ESS 结果 CSV 输出目录
    burnin_pct : float — burn-in 百分比 (默认 10%)
    ess_threshold : float — ESS 达标阈值 (默认 200)

    Returns
    -------
    dict: {success, ess_table, n_params, n_passed, n_warned, worst_params, error}
    """
    log = log or LogCollector()
    result = {"success": False, "ess_table": None, "n_params": 0,
              "n_passed": 0, "n_warned": 0, "worst_params": [], "error": None}

    if not log_file or not os.path.exists(log_file):
        result["error"] = f"Log file not found: {log_file}"
        return result

    try:
        # ── 1. Read BEAST log (tab-separated with # comments) ──
        df = pd.read_csv(log_file, sep="\t", comment="#")
        n_total = len(df)
        if n_total == 0:
            result["error"] = "Empty log file"
            return result

        # ── 2. Apply burn-in ──
        burnin_idx = int(n_total * burnin_pct / 100.0)
        df_burnin = df.iloc[burnin_idx:].reset_index(drop=True)
        n_used = len(df_burnin)
        log.emit(f"ESS analysis: {n_total} samples, {burnin_pct}% burn-in → {n_used} retained")

        # ── 3. Calculate ESS for all numeric columns ──
        skip_cols = {"sample", "Sample", "state", "State"}
        ess_rows = []
        for col in df_burnin.select_dtypes(include=[np.number]).columns:
            if col in skip_cols:
                continue
            ess_val = _effective_sample_size(df_burnin[col].values)
            ess_rows.append({"Parameter": col, "ESS": round(ess_val, 1)})

        if not ess_rows:
            result["error"] = "No numeric columns found in log"
            return result

        ess_df = pd.DataFrame(ess_rows).sort_values("ESS")
        result["n_params"] = len(ess_df)

        # ── 4. Save ESS table ──
        os.makedirs(output_dir, exist_ok=True)
        csv_path = os.path.join(output_dir, "ess_results.csv")
        ess_df.to_csv(csv_path, index=False)
        result["ess_table"] = csv_path

        # ── 5. Check thresholds ──
        passed = ess_df[ess_df["ESS"] >= ess_threshold]
        warned = ess_df[ess_df["ESS"] < ess_threshold]
        result["n_passed"] = len(passed)
        result["n_warned"] = len(warned)

        if len(warned) == 0:
            log.emit(f"ESS: all {len(ess_df)} parameters >= {ess_threshold} ✓")
        else:
            worst = warned.head(10)
            result["worst_params"] = [
                {"param": r["Parameter"], "ess": r["ESS"]}
                for _, r in worst.iterrows()
            ]
            log.warning(f"ESS: {len(warned)}/{len(ess_df)} parameters < {ess_threshold}:")
            for wp in result["worst_params"]:
                log.warning(f"  {wp['param']}: ESS={wp['ess']:.1f}")

        # ── 6. Check key parameters explicitly ──
        key_params = {"posterior", "prior", "likelihood", "clockRate"}
        for kp in key_params:
            matches = ess_df[ess_df["Parameter"].str.contains(kp, case=False, na=False)]
            if not matches.empty:
                min_ess = matches["ESS"].min()
                if min_ess < ess_threshold:
                    log.warning(f"  ⚠ {kp}: min ESS={min_ess:.1f} < {ess_threshold}")

        result["success"] = True
        return result

    except Exception as e:
        result["error"] = str(e)
        log.error(f"ESS analysis failed: {e}")
        return result


# ═══════════════════════════════════════════════════════════════════
# BEAST 后处理 (BSP-Viz + RSPP-Viz + ESS)
# ═══════════════════════════════════════════════════════════════════

def postprocess_beast(
    log_file: str,
    trees_file: str,
    output_dir: str,
    burnin_pct: float = 10.0,
    rscript_path: str = "Rscript",
    log: Optional[LogCollector] = None,
) -> Dict:
    """
    对 BEAST 输出进行后处理：BSP + RSPP 可视化。

    Parameters
    ----------
    log_file : str — .log 文件
    trees_file : str — .trees 文件 (或 MCC 树)
    output_dir : str
    burnin_pct : float — burn-in 百分比
    rscript_path : str — Rscript 路径
    log : LogCollector, optional

    Returns
    -------
    dict: {success, bsp_plot, rspp_plot, rspp_states, error}
    """
    from utils.virphy_bridge import run_bsp_viz, run_rspp_viz

    log = log or LogCollector()
    os.makedirs(output_dir, exist_ok=True)
    result = {"success": False, "bsp_plot": None, "rspp_plot": None,
              "rspp_states": {}, "ess": None, "error": None}

    # ── ESS convergence check ──
    if log_file and os.path.exists(log_file):
        log.emit(f"Post-processing: ESS analysis ({log_file})")
        ess = analyze_beast_ess(log_file, output_dir, burnin_pct=burnin_pct, log=log)
        result["ess"] = ess

    # ── BSP-Viz ──
    if log_file and os.path.exists(log_file):
        try:
            log.emit(f"Post-processing: BSP-Viz ({log_file})")
            bsp = run_bsp_viz(log_file, output_dir, burnin_pct=burnin_pct,
                             rscript_path=rscript_path, log=log)
            if bsp["success"]:
                result["bsp_plot"] = bsp["plot"]
        except Exception as e:
            log.warning(f"BSP-Viz skipped ({e})")

    # ── RSPP-Viz ──
    if trees_file and os.path.exists(trees_file):
        try:
            log.emit(f"Post-processing: RSPP-Viz ({trees_file})")
            rspp = run_rspp_viz(trees_file, output_dir, log=log)
            if rspp["success"]:
                result["rspp_plot"] = rspp["plot"]
                result["rspp_states"] = rspp.get("root_states", {})
        except Exception as e:
            log.warning(f"RSPP-Viz skipped ({e})")

    result["success"] = bool(result["bsp_plot"] or result["rspp_plot"])
    return result


# ═══════════════════════════════════════════════════════════════════
# 一枪头: 生成 → 运行 → 后处理
# ═══════════════════════════════════════════════════════════════════

def run_beast_pipeline(
    fasta_file: str,
    dates_csv: str,
    output_dir: str,
    beast_bin: str = "beast",
    threads: int = 8,
    chain_length: int = 5_000_000,
    burnin_pct: float = 10.0,
    rscript_path: str = "Rscript",
    skip_run: bool = False,
    log: Optional[LogCollector] = None,
) -> Dict:
    """
    完整 BEAST 管线: 生成 XML → 运行 BEAST → 后处理。

    这是整合后可独立调用的便利函数。
    """
    log = log or LogCollector()
    result = {"success": False, "xml": None, "run": None, "postprocess": None}

    # Step 1: Generate XML
    log.emit("=" * 50)
    log.emit("Step 1/3: Generate BEAST2 XML")
    xml_result = generate_beast2_xml(fasta_file, dates_csv, output_dir,
                                     chain_length=chain_length, log=log)
    result["xml"] = xml_result
    if not xml_result["success"]:
        result["error"] = f"XML generation failed: {xml_result.get('error')}"
        return result

    # Step 2: Run BEAST
    if skip_run:
        log.emit("Step 2/3: Skip BEAST run (--skip_run)")
    else:
        log.emit("=" * 50)
        log.emit("Step 2/3: Run BEAST")
        beast_result = run_beast(xml_result["xml_path"], output_dir,
                                beast_bin=beast_bin, threads=threads, log=log)
        result["run"] = beast_result

    # Step 3: Post-process
    log.emit("=" * 50)
    log.emit("Step 3/3: Post-process BEAST output")
    # 历史坑: 硬编码 beast2.log/beast2.trees 丢弃 run_beast 返回的真实路径;
    # fallback 只在 log 缺失时触发 → trees 名字不同时静默指向不存在的文件
    log_path = beast_result.get("log_path") if beast_result else None
    trees_path = beast_result.get("trees_path") if beast_result else None
    log_path = log_path if log_path and os.path.exists(log_path) else os.path.join(output_dir, "beast2.log")
    trees_path = trees_path if trees_path and os.path.exists(trees_path) else os.path.join(output_dir, "beast2.trees")
    if not os.path.exists(log_path) or not os.path.exists(trees_path):
        # 独立探测任意 .log/.trees
        for f in os.listdir(output_dir):
            fp = os.path.join(output_dir, f)
            if f.endswith('.log') and not os.path.exists(log_path):
                log_path = fp
            elif f.endswith('.trees') and not os.path.exists(trees_path):
                trees_path = fp

    post = postprocess_beast(log_path, trees_path, os.path.join(output_dir, "postprocess"),
                            burnin_pct=burnin_pct, rscript_path=rscript_path, log=log)
    result["postprocess"] = post
    # success 综合三步, 不再只看 XML 生成
    result["success"] = bool(xml_result.get("success")
                             and (not beast_result or beast_result.get("success"))
                             and post.get("success", True))
    return result
