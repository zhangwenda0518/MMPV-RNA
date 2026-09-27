#!/usr/bin/env python3
"""
phylogeo_bridge.py — 贝叶斯系统地理学 BEAST2 XML 生成 + 后处理
==============================================================
实现论文标准流程中的核心环节:
  1. 系统地理学 XML 生成 (discrete trait CTMC + BSSVS + UCLN + Skyline)
  2. TreeAnnotator MCC 树提取
  3. Bayes Factor 迁移路径筛选
  4. SpreaD3 JSON 生成

参考文献:
  - Lemey et al. (2009) Bayesian phylogeography finds its roots. PLoS Comput Biol.
  - Bielejec et al. (2016) SpreaD3: interactive visualization of spatiotemporal history.
  - TuMV Silk Road (PNAS 2021), PVS Europe (Virology 2018), ToMV Eurasia (Virology 2021)
"""

import csv, os, re, sys, subprocess, json, time
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from collections import defaultdict

import numpy as np
import pandas as pd
from Bio import SeqIO

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from utils.virphy_bridge import LogCollector
# 2026-09-15 (审查 P2-9): BF 阈值/发散截断/公式 全管线唯一定义
from utils.virphy_bridge import (DEFAULT_BF_THRESHOLD, BF_INFINITE_CAP,
                                 bayes_factor, find_indicator_columns)


# ═══════════════════════════════════════════════════════════════════
# 1. 系统地理学 BEAST2 XML 生成
# ═══════════════════════════════════════════════════════════════════

def generate_phylogeo_xml(
    fasta_file: str,
    metadata_csv: str,
    output_dir: str,
    chain_length: int = 10_000_000,
    log_every: int = 5000,
    clock_model: str = "ucln",          # "strict" or "ucln"
    tree_prior: str = "skyline",         # "constant", "skyline", or "skygrid"
    substitution_model: str = "HKY",     # "HKY", "GTR", or "bModelTest"
    discretize_locations: bool = True,
    location_column: str = "location",
    n_location_groups: int = 10,         # max discrete states
    log: Optional[LogCollector] = None,
) -> Dict:
    """
    生成带有系统地理学 (discrete phylogeography) 的 BEAST2 XML。

    BEAST2 需要的包 (需提前安装):
      BEAST v2.7.x
      BEASTLabs (for BSSVS, Skyline)
      ORC (for Skygrid)
      bModelTest (optional)

    Parameters
    ----------
    fasta_file : str — MAFFT 比对 FASTA
    metadata_csv : str — CSV: name, date, location [, host]
    output_dir : str
    chain_length : int — MCMC 链长 (phylogeography 建议 10M+)
    log_every : int
    clock_model : str — "strict" or "ucln" (uncorrelated lognormal)
    tree_prior : str — "constant", "skyline"
    substitution_model : str — "HKY" or "GTR"
    discretize_locations : bool — 将连续地点 discretize 为区域
    location_column : str — metadata 中的地点列名
    n_location_groups : int — 最大离散状态数
    log : LogCollector

    Returns
    -------
    dict: {success, xml_path, n_taxa, n_locations, locations, seq_len, year_range}
    """
    log = log or LogCollector()
    os.makedirs(output_dir, exist_ok=True)
    result = {"success": False, "xml_path": None, "n_taxa": 0,
              "n_locations": 0, "locations": [], "seq_len": 0,
              "min_year": None, "max_year": None, "error": None}

    # ══════════════════════════════════════════════════════════════════
    # 2026-09-15 修复 (P1-11): BEAST2 生成器整体停用 (fail-loud)
    # ══════════════════════════════════════════════════════════════════
    # 复核结论: 本函数并不是"有几处 bug 的 BEAST2 生成器", 而是**用 BEAST1 语法
    # 套 `spec=` 属性**拼出来的 XML —— 与 BEAST2 的 XML 模型不兼容。逐条对照
    # BEAST2 官方示例 examples/parameterised/RSV2.xml 与 BEAST2 源码:
    #
    #  (1) 先验写法: 本地 `<prior id="KappaPrior" name="distribution">
    #       <parameter idref="kappa.s:aln"/><LogNormal name="distr" .../></prior>`
    #      BEAST2 的 beast.base.inference.distribution.Prior 里 `x` 是
    #      Validate.REQUIRED (Prior.java: m_x), 正确写法是
    #      `<prior ... name="distribution" x="@kappa.s:aln"><LogNormal name="distr">`.
    #      缺 x 会让每个先验直接解析失败。同问题见 clockPrior / ucldMeanPrior /
    #      ucldStdevPrior / popSizePrior (共 5 处)。
    #  (2) `spec="OneOnXPrior"` —— BEAST2 无此类 (应为
    #      beast.base.inference.distribution.OneOnX)。
    #  (3) `spec="UCLN.ClockModel"` —— BEAST2 无此类 (应为
    #      beast.base.evolution.branchratemodel.UCRelaxedClockModel)。
    #  (4) `spec="SumLogger"` —— BEAST2 无此类 (计数统计量应为
    #      beast.base.evolution.Sum, 输入名 arg)。
    #  (5) `spec="TreeHeightLogger"` 未限定包名, 而 namespace 里没有
    #      beast.base.evolution.tree —— 解析不到。
    #  (6) `spec="StrictClockModel"` 同样未限定且缺包
    #      (beast.base.evolution.branchratemodel)。
    #  (7) `<up><parameter idref="..."/></up>` 是 BEAST1 写法; BEAST2 为
    #      `<up idref="..."/>`。
    #  (8) `rateIndicator=` / `symmetric=` 不是 BEAST2 GeneralSubstitutionModel
    #      的合法输入 (源码只声明 rates / eigenSystem; 父类只有 frequencies)。
    #  (9) `beast.math.distributions.Poisson` 类不存在。
    # (10) `namespace="".join(ns)` 无分隔符 (BEAST2 要求冒号分隔)。
    #
    # 修好 (10)(8)(9) 这类"单点"问题不足以让产出可用 —— (1)(2)(3)(4) 属于
    # 整份 XML 的结构性错误, 要正确重写必须**在本机装 BEAST2 逐项验证**。
    # 本机未安装 BEAST2, 仓库内也没有可参照的 BEAST2 DTA 示例。
    # 按"半成品不往外发 / 结论必须站得住"的原则: 明确失败, 不再产出
    # 会崩 (更糟: 静默失效) 的 XML。BEAST1 路径 (管线默认) 不受影响。
    result["error"] = (
        "generate_phylogeo_xml (BEAST2) 已停用: 该生成器使用 BEAST1 语法, "
        "产出的 XML 与 BEAST2 不兼容 —— 至少 10 处独立违规, 其中 5 处是 "
        "Prior 缺 REQUIRED 的 x 属性 (整份 XML 必然解析失败)。"
        "请使用 --beast_version 1 (BEAST1 路径, 管线默认, 其 BSSVS 实现已对照 "
        "biosoft/VirPhyKit/Example/MJRM Generator/PVS/PVS_no_matrix.xml 校验)。"
        "若需要 BEAST2 支持, 请先在本机安装 BEAST2 并按 RSV2.xml 的 XML 习惯重写本函数。")
    log.error(result["error"])
    return result

    try:
        # ── 1. Load data ──
        aln = list(SeqIO.parse(fasta_file, "fasta"))
        if len(aln) < 4:
            result["error"] = f"Need >=4 sequences, got {len(aln)}"
            return result
        seq_len = len(str(aln[0].seq))

        # ── 2. Load metadata ──
        meta = {}
        with open(metadata_csv, encoding='utf-8-sig') as f:
            reader = csv.DictReader(f)
            for row in reader:
                name = row.get('name', '').strip()
                date = row.get('date', '').strip()
                loc = row.get(location_column, '').strip()
                if name and date:
                    meta[name] = {"date": date, "location": loc if loc else "Unknown"}

        # Filter to dated + located sequences
        valid = [(r, meta[r.id]) for r in aln if r.id in meta]
        if len(valid) < 4:
            result["error"] = f"Only {len(valid)} sequences with dates+locations"
            return result

        # ── 3. Discretize locations ──
        all_locs = [m["location"] for _, m in valid]
        loc_counts = defaultdict(int)
        for loc in all_locs:
            simplified = _simplify_location(loc)
            loc_counts[simplified] += 1

        # Sort by frequency, keep top N
        sorted_locs = sorted(loc_counts.items(), key=lambda x: -x[1])
        if len(sorted_locs) > n_location_groups:
            # Merge small groups into "Other"
            top_locs = dict(sorted_locs[:n_location_groups - 1])
            top_locs["Other"] = sum(c for _, c in sorted_locs[n_location_groups - 1:])
            loc_map = {}
            for loc in all_locs:
                simplified = _simplify_location(loc)
                loc_map[loc] = simplified if simplified in top_locs else "Other"
        else:
            top_locs = dict(sorted_locs)
            loc_map = {loc: _simplify_location(loc) for loc in all_locs}

        unique_locs = sorted(top_locs.keys())
        result["locations"] = unique_locs
        result["n_locations"] = len(unique_locs)
        n_loc = len(unique_locs)
        if n_loc < 2:
            log.warning(f"  ⚠ Only {n_loc} location found — skipping phylogeography (need ≥2).")
        else:
            log.emit(f"Discretized locations: {n_loc} states: {unique_locs}")

        # 注: 本函数已在入口处整体 fail-loud 停用 (见函数开头的 (1)-(10) 清单),
        # 以下代码当前**不可达**, 仅为日后按 BEAST2 习惯重写时保留参考。
        # 其中 BSSVS 接线 (GeneralSubstitutionModel 的 rateIndicator / symmetric)
        # 已确认不是 BEAST2 2.7 的合法输入, 相关属性已删除, 详见下方 substModel 处注释。

        # ── 4. Convert dates to decimal years ──
        # 2026-09-15 (审查 P2-26): 旧版 _to_yr 丢弃"日"且解析失败静默返回 2020.0。
        # 统一到 utils/decimal_year (全管线唯一定义)。
        from utils.decimal_year import to_decimal_year as _to_yr

        years = [_to_yr(m["date"]) for _, m in valid]
        max_year = max(years)
        min_year = min(years)
        result["n_taxa"] = len(valid)
        result["seq_len"] = seq_len
        result["min_year"] = min_year
        result["max_year"] = max_year
        log.emit(f"{len(valid)} taxa, {seq_len} bp, years {min_year:.1f}–{max_year:.1f}")

        # ── 5. Generate XML ──
        xml_path = os.path.join(output_dir, "phylogeo.xml")
        output_prefix = os.path.join(output_dir, "phylogeo")

        # Namespace definitions
        # 2026-09-15 修复 (P1-11): 原写法 `"".join(ns)` 把包名**直接首尾相接**,
        # 产出 `beast.corebeast.evolution.alignmentbeast...` —— BEAST2 无法解析出任何包,
        # 于是所有未限定全名的 spec= 都会找不到类。
        # BEAST2 的 namespace 属性是**冒号分隔**的包名列表, 判据来自本仓库自带的
        # 真实 BEAST2 XML 模板 archive/_nonviral_misc/coal_skyline_tpl.xml 以及
        # BEAST2 官方示例 examples/parameterised/RSV2.xml:
        #   namespace="beast.pkgmgmt:beast.base.core:beast.base.inference:..."
        # 同时补齐本文件实际用到但原先缺失的包 (beast.base.evolution.tree:
        # RandomTree/TreeHeightLogger/TreeIntervals 都在这)。
        ns = [
            'beast.pkgmgmt',
            'beast.base.core', 'beast.base.inference',
            'beast.base.inference.distribution', 'beast.base.inference.operator',
            'beast.base.inference.parameter', 'beast.base.inference.util',
            'beast.base.evolution', 'beast.base.evolution.alignment',
            'beast.base.evolution.tree', 'beast.base.evolution.tree.coalescent',
            'beast.base.evolution.branchratemodel',
            'beast.base.evolution.operator', 'beast.base.evolution.sitemodel',
            'beast.base.evolution.substitutionmodel', 'beast.base.evolution.likelihood',
            # 旧包名 (BEAST2 2.x 兼容层, 部分类仍以 beast.evolution.* 解析)
            'beast.evolution.alignment', 'beast.evolution.tree.coalescent',
            'beast.evolution.operators', 'beast.evolution.sitemodel',
            'beast.evolution.substitutionmodel', 'beast.evolution.likelihood',
            # BEASTLabs
            'beastlabs.evolution.likelihood', 'beastlabs.evolution.tree',
        ]
        # 去重但保持顺序 (重复包名虽无害, 但会让生成的 XML 难以比对)
        ns = list(dict.fromkeys(ns))

        with open(xml_path, 'w', encoding='utf-8') as f:
            f.write('<?xml version="1.0" encoding="UTF-8" standalone="no"?>\n')
            f.write(f'<beast version="2.7" namespace="{":".join(ns)}">\n\n')

            # ── Alignment ──
            f.write(f'  <!-- {len(valid)} taxa, {seq_len} bp, {min_year:.1f}–{max_year:.1f} -->\n')
            f.write('  <data id="aln" spec="Alignment" name="alignment">\n')
            for rec, _ in valid:
                f.write(f'    <sequence taxon="{_xml_escape(rec.id)}" '
                        f'value="{str(rec.seq)}" />\n')
            f.write('  </data>\n\n')

            # ── Tip dates trait (date-forward for tip dates) ──
            f.write('  <trait id="dateTrait" spec="beast.base.evolution.tree.TraitSet"\n')
            f.write('         traitname="date-forward"\n')
            f.write('         value="\n')
            date_strs = [f'    {_xml_escape(r.id)}={_to_yr(m["date"]):.4f}'
                        for r, m in valid]
            f.write(',\n'.join(date_strs))
            f.write('\n  ">\n')
            f.write('    <taxa id="taxa" spec="TaxonSet" alignment="@aln"/>\n')
            f.write('  </trait>\n\n')

            # ── Location trait (discrete phylogeography) ──
            if discretize_locations and len(unique_locs) > 1:
                f.write('  <!-- Discrete phylogeography trait -->\n')
                f.write('  <trait id="locationTrait" spec="beast.base.evolution.tree.TraitSet"\n')
                f.write('         traitname="location"\n')
                f.write('         value="\n')
                loc_strs = [f'    {_xml_escape(r.id)}={loc_map[m["location"]]}'
                           for r, m in valid]
                f.write(',\n'.join(loc_strs))
                f.write('\n  ">\n')
                f.write('    <taxa idref="taxa"/>\n')
                f.write('  </trait>\n\n')

            # ── Initial tree ──
            tree_str = '(' + ','.join([_xml_escape(r.id) for r, _ in valid]) + ')'

            # ── MCMC Run ──
            f.write(f'  <run id="mcmc" spec="MCMC" chainLength="{chain_length}" '
                    f'storeEvery="{log_every}">\n\n')
            f.write('    <state id="state" storeEvery="5000">\n')
            f.write(f'      <tree id="Tree.t:aln" spec="beast.base.evolution.tree.TreeParser" '
                    f'taxa="@aln" IsLabelledNewick="true"\n')
            f.write(f'            newick="{tree_str}" />\n')

            # Clock rate
            f.write('      <parameter id="clockRate.c:aln" spec="RealParameter"\n')
            f.write('                 lower="1.0E-9" upper="0.1" name="stateNode">'
                    '0.001</parameter>\n')

            # UCLN parameters
            if clock_model == "ucln":
                f.write('      <parameter id="ucldMean.c:aln" spec="RealParameter"\n')
                f.write('                 lower="0.0" name="stateNode">1.0</parameter>\n')
                f.write('      <parameter id="ucldStdev.c:aln" spec="RealParameter"\n')
                f.write('                 lower="0.0" upper="10.0" name="stateNode">'
                        '0.3333333333333333</parameter>\n')

            # Population size
            f.write('      <parameter id="popSize.t:aln" spec="RealParameter"\n')
            f.write('                 lower="0.0" upper="1.0E100" name="stateNode">'
                    '0.3</parameter>\n')

            # Substitution model params
            f.write('      <parameter id="kappa.s:aln" spec="RealParameter"\n')
            f.write('                 lower="0.0" name="stateNode">2.0</parameter>\n')
            f.write('      <parameter id="freqParameter.s:aln" spec="RealParameter"\n')
            f.write('                 dimension="4" lower="0.0" upper="1.0" name="stateNode">'
                    '0.25 0.25 0.25 0.25</parameter>\n')

            # Location CTMC rate indicators (BSSVS)
            if discretize_locations and len(unique_locs) > 1:
                n_rates = len(unique_locs) * (len(unique_locs) - 1)
                # 2026-09-15 修复 (P1-11): 原实现把 n_rates 算出来却**从未使用**,
                # 两个参数的 dimension 都写死 "2"、初值也只有两个 "1.0"。
                # BEAST2 的 GeneralSubstitutionModel.initAndValidate 会直接抛:
                #   "Dimension of input 'rates' is 2 but a rate matrix of dimension
                #    Nx(N-1)=M was expected"
                # 即 N 个地点必须有 N*(N-1) 个速率 —— 只有 N=2 时写死值才碰巧正确。
                # 判据: BEAST2 源码 GeneralSubstitutionModel.java (dimension 断言) +
                # 本仓库 BEAST1 兄弟实现 utils/beast1_bridge.py:538 用的就是 {n_rates} +
                # 外部参考 biosoft/VirPhyKit/.../PVS_no_matrix.xml: 6 个地区 -> dimension="30"。
                f.write('      <parameter id="location.indicators" spec="parameter.BooleanParameter"\n')
                f.write(f'                 dimension="{n_rates}" name="stateNode">'
                        f'{" ".join(["true"] * n_rates)}</parameter>\n')
                f.write('      <parameter id="location.rates" spec="parameter.RealParameter"\n')
                f.write(f'                 dimension="{n_rates}" lower="0.0" name="stateNode">'
                        f'{" ".join(["1.0"] * n_rates)}</parameter>\n')

            f.write('    </state>\n\n')

            # ── INIT (Random tree) ──
            f.write('    <init id="RandomTree.t:aln" spec="beast.base.evolution.tree.RandomTree"\n')
            f.write('          estimate="false" initial="@Tree.t:aln" taxa="@aln">\n')
            f.write('      <populationModel id="ConstantPopulationInit.t:aln"\n')
            f.write('                        spec="beast.evolution.tree.coalescent.ConstantPopulation">\n')
            f.write('        <parameter id="initPopSize.t:aln" spec="RealParameter" '
                    'name="popSize">0.3</parameter>\n')
            f.write('      </populationModel>\n')
            f.write('    </init>\n\n')

            # ── POSTERIOR ──
            f.write('    <distribution id="posterior" spec="CompoundDistribution">\n')
            f.write('      <distribution id="prior" spec="CompoundDistribution">\n')

            # Coalescent tree prior
            f.write('        <distribution id="Coalescent" spec="Coalescent">\n')
            if tree_prior == "skyline":
                # Bayesian Skyline with 5 groups
                f.write('          <populationModel id="SkylinePopulation"\n')
                f.write('                            spec="beast.evolution.tree.coalescent.BayesianSkyline">\n')
                f.write('            <groupSizes id="groupSizes" spec="RealParameter" '
                        'dimension="5" lower="1">2 2 2 2 24</groupSizes>\n')
                f.write('            <parameter idref="popSize.t:aln"/>\n')
                f.write('          </populationModel>\n')
            else:
                f.write('          <populationModel id="ConstantPopulation" spec="ConstantPopulation"\n')
                f.write('                            populationSize="@popSize.t:aln"/>\n')
            f.write('          <treeIntervals id="TreeIntervals" spec="TreeIntervals" '
                    'tree="@Tree.t:aln"/>\n')
            f.write('        </distribution>\n')

            # 2026-09-15 修复 (P1-11): 以下 5 个先验原先都是 BEAST1 写法 ——
            #   <prior id="X" name="distribution"><parameter idref="p"/><LogNormal name="distr"/></prior>
            # BEAST2 的 beast.base.inference.distribution.Prior 里 `x` 是
            # Validate.REQUIRED, 且目标通过 **属性** `x="@p"` 指定 (不是子元素
            # <parameter idref>)。正确写法见 BEAST2 官方示例 RSV2.xml:
            #   <prior id="KappaPrior.s:RSV2" name="distribution" x="@kappa.s:RSV2">
            #     <LogNormal id="..." name="distr">
            #       <parameter id="..." spec="parameter.RealParameter"
            #                  estimate="false" name="M">1.0</parameter>
            #       ...
            #     </LogNormal>
            #   </prior>
            # 另外 spec="OneOnXPrior" 在 BEAST2 中不存在 (应为 OneOnX)。

            # PopSize prior (OneOnX)
            f.write('        <prior id="popSizePrior" name="distribution" '
                    'x="@popSize.t:aln">\n')
            f.write('          <OneOnX id="popSizeDistr" name="distr"/>\n')
            f.write('        </prior>\n')

            # Kappa prior
            f.write('        <prior id="KappaPrior" name="distribution" '
                    'x="@kappa.s:aln">\n')
            f.write('          <LogNormal id="kappaDistr" name="distr">\n')
            f.write('            <parameter id="kappaDistrM" spec="parameter.RealParameter" '
                    'estimate="false" name="M">1.0</parameter>\n')
            f.write('            <parameter id="kappaDistrS" spec="parameter.RealParameter" '
                    'estimate="false" name="S">1.25</parameter>\n')
            f.write('          </LogNormal>\n')
            f.write('        </prior>\n')

            # Clock rate prior
            f.write('        <prior id="clockPrior" name="distribution" '
                    'x="@clockRate.c:aln">\n')
            f.write('          <LogNormal id="clockDistr" name="distr">\n')
            f.write('            <parameter id="clockDistrM" spec="parameter.RealParameter" '
                    'estimate="false" name="M">0.001</parameter>\n')
            f.write('            <parameter id="clockDistrS" spec="parameter.RealParameter" '
                    'estimate="false" name="S">1.25</parameter>\n')
            f.write('          </LogNormal>\n')
            f.write('        </prior>\n')

            # UCLN priors
            if clock_model == "ucln":
                f.write('        <prior id="ucldMeanPrior" name="distribution" '
                        'x="@ucldMean.c:aln">\n')
                f.write('          <Exponential id="ucldMeanDistr" name="distr">\n')
                f.write('            <parameter id="ucldMeanDistrMean" '
                        'spec="parameter.RealParameter" estimate="false" '
                        'name="mean">0.3333333333333333</parameter>\n')
                f.write('          </Exponential>\n')
                f.write('        </prior>\n')
                f.write('        <prior id="ucldStdevPrior" name="distribution" '
                        'x="@ucldStdev.c:aln">\n')
                f.write('          <Exponential id="ucldStdevDistr" name="distr">\n')
                f.write('            <parameter id="ucldStdevDistrMean" '
                        'spec="parameter.RealParameter" estimate="false" '
                        'name="mean">0.3333333333333333</parameter>\n')
                f.write('          </Exponential>\n')
                f.write('        </prior>\n')

            # Location CTMC prior (BSSVS)
            if discretize_locations and len(unique_locs) > 1:
                f.write('\n        <!-- Phylogeography: CTMC + BSSVS -->\n')
                # 2026-09-15 修复 (P1-11): 原实现三处错误 ——
                #  (1) `beast.math.distributions.Poisson` 在 BEAST2 中不存在
                #      (正确为 beast.base.inference.distribution.Poisson);
                #  (2) Poisson 被当作独立 <distribution> 使用, 但它继承自
                #      ParametricDistribution, **没有 x 输入** —— 必须先包一层
                #      beast.base.inference.distribution.Prior (x= 与 distr= 皆
                #      Validate.REQUIRED, 见 Prior.java), 否则该分布不作用于任何
                #      变量, BSSVS 先验完全空转;
                #  (3) location.lambda 写 lower="0.0" upper="0.0" 而取值 0.6931 ——
                #      区间宽度为 0 且取值越界; RealParameter 会直接报错, 且
                #      Prior.calculateLogP() 对越界的 x 一律返回 -Infinity。
                # 非零计数用 beast.base.evolution.Sum: 其 functionInput 名为 "arg",
                # 对 Boolean/Integer 参数走 integer_mode, 求和即非零个数。
                f.write('        <sumStatistic id="location.nonZeroRates" '
                        'spec="beast.base.evolution.Sum" '
                        'arg="@location.indicators"/>\n')
                f.write('        <prior id="location.indicatorsPrior" name="distribution" '
                        'x="@location.nonZeroRates">\n')
                f.write('          <Poisson id="location.PoissonDistr" name="distr">\n')
                f.write('            <parameter id="location.lambda" '
                        'spec="parameter.RealParameter" estimate="false" '
                        'name="lambda">0.6931471805599453</parameter>\n')
                f.write('          </Poisson>\n')
                f.write('        </prior>\n')
                # location.rates 必须有 proper prior, 否则 MCMC 会把速率推到无穷大
                # -> rate matrix 病态 -> 特征分解不收敛而崩溃。
                # (BEAST1 兄弟实现 utils/beast1_bridge.py 已就此留下教训注释:
                #  "PSTVd/GCVA 均踩过"; BEAST2 版本原先完全没有这个先验。)
                f.write('        <prior id="location.ratesPrior" name="distribution" '
                        'x="@location.rates">\n')
                f.write('          <Exponential id="location.ratesDistr" name="distr">\n')
                f.write('            <parameter id="location.ratesMean" '
                        'spec="parameter.RealParameter" estimate="false" '
                        'name="mean">0.5</parameter>\n')
                f.write('          </Exponential>\n')
                f.write('        </prior>\n')

            f.write('      </distribution>\n')

            # ── LIKELIHOOD ──
            f.write('\n      <distribution id="likelihood" spec="CompoundDistribution">\n')
            f.write('        <distribution id="treeLikelihood.aln" spec="TreeLikelihood"\n')
            f.write('                       data="@aln" tree="@Tree.t:aln">\n')
            f.write('          <siteModel id="SiteModel.aln" spec="SiteModel">\n')
            f.write('            <parameter spec="RealParameter" estimate="false" '
                    'name="mutationRate">1.0</parameter>\n')
            f.write('            <parameter id="gammaShape.aln" spec="RealParameter"\n')
            f.write('                       estimate="false" name="shape">1.0</parameter>\n')
            f.write('            <parameter id="propInvariant.aln" spec="RealParameter"\n')
            f.write('                       lower="0.0" name="proportionInvariant" '
                    'upper="1.0">0.0</parameter>\n')
            f.write('            <substModel id="hky.aln" spec="HKY" '
                    'kappa="@kappa.s:aln">\n')
            f.write('              <frequencies id="estimatedFreqs.aln" spec="Frequencies"\n')
            f.write('                            frequencies="@freqParameter.s:aln"/>\n')
            f.write('            </substModel>\n')
            f.write('          </siteModel>\n')

            # Branch rate model
            # 2026-09-15 修复 (P1-11): 原 spec 值都不是 BEAST2 的类 ——
            #   "StrictClockModel" 未限定包名, 且 namespace 里没有
            #      beast.base.evolution.branchratemodel;
            #   "UCLN.ClockModel" 在 BEAST2 中根本不存在 (BEAST1 专有)。
            # BEAST2 的正确类名 (判据: RSV2.xml 与 BEAST2 源码):
            #   beast.base.evolution.branchratemodel.StrictClockModel
            #   beast.base.evolution.branchratemodel.UCRelaxedClockModel
            # 注意: UCRelaxedClockModel 的参数名是 clock.rate / rates / stdev,
            # 这里保留原有属性写法, 重新启用前需在本机 BEAST2 上验证。
            if clock_model == "strict":
                f.write('          <branchRateModel id="StrictClock.c:aln"\n')
                f.write('                            spec="beast.base.evolution.branchratemodel.StrictClockModel"\n')
                f.write('                            clock.rate="@clockRate.c:aln"/>\n')
            else:
                f.write('          <branchRateModel id="UCLN.c:aln"\n')
                f.write('                            spec="beast.base.evolution.branchratemodel.UCRelaxedClockModel"\n')
                f.write('                            clock.rate="@clockRate.c:aln"\n')
                f.write('                            rate="@ucldMean.c:aln"\n')
                f.write('                            stdev="@ucldStdev.c:aln"/>\n')

            f.write('        </distribution>\n')

            # Location CTMC likelihood
            if discretize_locations and len(unique_locs) > 1:
                f.write('\n        <!-- Phylogeography CTMC likelihood -->\n')
                f.write('        <distribution id="phylogeo.likelihood"\n')
                f.write('                       spec="TreeLikelihood"\n')
                f.write('                       data="@locationTrait"\n')
                f.write('                       tree="@Tree.t:aln">\n')
                f.write('          <siteModel id="location.sitemodel" spec="SiteModel">\n')
                # 2026-09-15 修复 (P1-11): 删掉 rateIndicator= 与 symmetric= 两个属性。
                # BEAST2 2.7 的 GeneralSubstitutionModel 只声明 rates / eigenSystem
                # (父类 SubstitutionModel.Base 只声明 frequencies), 这两个属性不是
                # 合法输入, 会被 XML 解析器拒绝。原本它们承担 BSSVS 的接线作用,
                # 删掉后此处仅剩一个普通(非 BSSVS)离散 CTMC —— 因此上面已在函数入口
                # 对 "BEAST2 + 离散地点" 整体显式拒绝, 本段当前不可达。
                # 保留修正后的代码, 供日后在本机 BEAST2 上验证接线后重新启用。
                f.write('            <substModel id="location.substmodel"\n')
                f.write('                         spec="GeneralSubstitutionModel"\n')
                f.write('                         rates="@location.rates">\n')
                f.write('              <frequencies id="location.freqs" spec="Frequencies"\n')
                f.write('                            data="@locationTrait"/>\n')
                f.write('            </substModel>\n')
                f.write('          </siteModel>\n')
                f.write('          <branchRateModel id="location.branchrate"\n')
                f.write('                            spec="beast.base.evolution.branchratemodel.StrictClockModel"\n')
                f.write('                            clock.rate="@clockRate.c:aln"/>\n')
                f.write('        </distribution>\n')

            f.write('      </distribution>\n')
            f.write('    </distribution>\n\n')

            # ── OPERATORS ──
            _write_phylogeo_operators(f, clock_model, discretize_locations, len(unique_locs))

            # ── LOGGERS ──
            _write_phylogeo_loggers(f, output_prefix, log_every,
                                   discretize_locations, len(unique_locs))

            f.write('  </run>\n</beast>\n')

        result["xml_path"] = xml_path
        result["success"] = True
        log.emit(f"Phylogeography BEAST2 XML: {xml_path} "
                f"({os.path.getsize(xml_path)} bytes)")
        log.emit(f"  Model: {clock_model} clock + {tree_prior} prior + HKY")
        if discretize_locations and len(unique_locs) > 1:
            log.emit(f"  Phylogeography: {len(unique_locs)} locations with BSSVS")
        return result

    except Exception as e:
        result["error"] = str(e)
        log.error(f"Phylogeography XML failed: {e}")
        import traceback
        traceback.print_exc()
        return result


def _write_phylogeo_operators(f, clock_model, has_geo, n_locations):
    """Write MCMC operators for phylogeographic BEAST2 XML"""
    ops = [
        ('StrictClockRateScaler', 'ScaleOperator',
         'parameter="@clockRate.c:aln" scaleFactor="0.75" weight="3"'),
        ('strictClockUpDownOperator', 'UpDownOperator',
         'scaleFactor="0.75" weight="3"'),
        ('KappaScaler', 'ScaleOperator',
         'parameter="@kappa.s:aln" scaleFactor="0.75" weight="1"'),
        ('FrequenciesExchanger', 'DeltaExchangeOperator',
         'parameter="@freqParameter.s:aln" delta="0.01" weight="1"'),
        ('CoalescentTreeScaler', 'ScaleOperator',
         'tree="@Tree.t:aln" scaleFactor="0.95" weight="3"'),
        ('CoalescentTreeRootScaler', 'ScaleOperator',
         'tree="@Tree.t:aln" rootOnly="true" scaleFactor="0.95" weight="3"'),
        ('CoalescentUniformOperator', 'Uniform',
         'tree="@Tree.t:aln" weight="30"'),
        ('CoalescentSubtreeSlide', 'SubtreeSlide',
         'tree="@Tree.t:aln" weight="15"'),
        ('CoalescentNarrow', 'Exchange',
         'tree="@Tree.t:aln" isNarrow="true" weight="15"'),
        ('CoalescentWide', 'Exchange',
         'tree="@Tree.t:aln" isNarrow="false" weight="3"'),
        ('CoalescentWilsonBalding', 'WilsonBalding',
         'tree="@Tree.t:aln" weight="3"'),
        ('PopSizeScaler', 'ScaleOperator',
         'parameter="@popSize.t:aln" scaleFactor="0.75" weight="3"'),
    ]

    if clock_model == "ucln":
        ops += [
            ('UCLDMeanScaler', 'ScaleOperator',
             'parameter="@ucldMean.c:aln" scaleFactor="0.75" weight="3"'),
            ('UCLDStdevScaler', 'ScaleOperator',
             'parameter="@ucldStdev.c:aln" scaleFactor="0.75" weight="3"'),
        ]

    if has_geo and n_locations > 1:
        ops += [
            ('location.ratesScaler', 'ScaleOperator',
             'parameter="@location.rates" scaleFactor="0.75" weight="3"'),
            ('location.indicatorsBitFlip', 'BitFlipOperator',
             'parameter="@location.indicators" weight="3"'),
        ]

    for op_id, op_spec, op_attrs in ops:
        if op_spec == 'UpDownOperator':
            # 2026-09-15 修复 (P1-11): BEAST2 的 UpDownOperator 用 `<up idref="..."/>`
            # 直接引用, 而不是 BEAST1 的 `<up><parameter idref="..."/></up>` 嵌套写法。
            # 判据: BEAST2 官方示例 examples/parameterised/RSV2.xml
            #   <operator id="strictClockUpDownOperator.c:RSV2" spec="UpDownOperator" ...>
            #     <up idref="clockRate.c:RSV2"/>
            #     <down idref="Tree.t:RSV2"/>
            #   </operator>
            f.write(f'    <operator id="{op_id}" spec="{op_spec}" {op_attrs}>\n')
            f.write('      <up idref="clockRate.c:aln"/>\n')
            f.write('      <down idref="Tree.t:aln"/>\n')
            f.write('    </operator>\n')
        else:
            f.write(f'    <operator id="{op_id}" spec="{op_spec}" {op_attrs}/>\n')


def _write_phylogeo_loggers(f, prefix, log_every, has_geo, n_locations):
    """Write loggers for phylogeographic BEAST2 XML"""
    # Screen logger
    f.write(f'\n    <logger id="screenlog" spec="Logger" logEvery="{log_every}">\n')
    for ref in ['posterior', 'likelihood', 'prior', 'clockRate.c:aln']:
        f.write(f'      <log idref="{ref}"/>\n')
    f.write('    </logger>\n\n')

    # Trace log
    f.write(f'    <logger id="tracelog" spec="Logger" fileName="{prefix}.log"\n')
    f.write(f'            logEvery="{log_every}" sort="smart">\n')
    for ref in ['posterior', 'likelihood', 'prior', 'clockRate.c:aln',
                'popSize.t:aln', 'kappa.s:aln']:
        f.write(f'      <log idref="{ref}"/>\n')
    # 2026-09-15 修复 (P1-11): TreeHeightLogger 必须写全包名 —— 原先只写
    # spec="TreeHeightLogger", 而 namespace 里没有 beast.base.evolution.tree,
    # BEAST2 解析不到该类。判据: RSV2.xml 用的是
    #   spec="beast.base.evolution.tree.TreeHeightLogger"
    f.write('      <log id="TreeHeight" '
            'spec="beast.base.evolution.tree.TreeHeightLogger" tree="@Tree.t:aln"/>\n')

    if has_geo and n_locations > 1:
        f.write('      <log idref="location.rates"/>\n')
        f.write('      <log idref="location.indicators"/>\n')
        # 2026-09-15 修复 (P1-11): BEAST2 里没有 SumLogger 这个类。
        # 计数统计量是 beast.base.evolution.Sum, 输入名 arg (不是 <parameter idref>)。
        f.write('      <log id="location.nonZero" '
                'spec="beast.base.evolution.Sum" arg="@location.indicators"/>\n')

    f.write('    </logger>\n\n')

    # Tree log
    f.write(f'    <logger id="treelog" spec="Logger" fileName="{prefix}.trees"\n')
    f.write(f'            logEvery="{log_every}" mode="tree">\n')
    f.write('      <log id="TreeWithMetaDataLogger" spec="TreeWithMetaDataLogger"\n')
    if has_geo and n_locations > 1:
        f.write('           tree="@Tree.t:aln">\n')
        f.write('        <trait idref="locationTrait"/>\n')
    else:
        f.write('           tree="@Tree.t:aln"/>\n')
    f.write('      </log>\n')
    f.write('    </logger>\n')


# ═══════════════════════════════════════════════════════════════════
# 2. TreeAnnotator — MCC 树提取
# ═══════════════════════════════════════════════════════════════════

def run_treeannotator(
    trees_file: str,
    output_dir: str,
    burnin_pct: float = 10.0,
    target_tree: str = "mcc",       # "mcc" or "target" (max clade credibility)
    treeannotator_bin: str = "treeannotator",
    log: Optional[LogCollector] = None,
) -> Dict:
    """
    运行 BEAST TreeAnnotator 从 .trees 文件提取 MCC 树。

    Returns
    -------
    dict: {success, mcc_tree, error}
    """
    log = log or LogCollector()
    result = {"success": False, "mcc_tree": None, "error": None}

    if not os.path.exists(trees_file):
        result["error"] = f"Trees file not found: {trees_file}"
        return result

    os.makedirs(output_dir, exist_ok=True)
    mcc_path = os.path.join(output_dir, "mcc.tree")

    # Calculate burnin as number of trees
    try:
        with open(trees_file) as f:
            n_trees = sum(1 for line in f if line.startswith('tree'))
    except Exception:
        n_trees = 10000  # fallback

    burnin_trees = max(1, int(n_trees * burnin_pct / 100.0))
    log.emit(f"TreeAnnotator: {n_trees} trees, burn-in={burnin_trees} ({burnin_pct}%)")

    cmd = [
        treeannotator_bin,
        "-burnin", str(burnin_trees),
        "-heights", "median",
        trees_file,
        mcc_path,
    ]
    log.emit(f"  {' '.join(cmd)}")

    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, encoding='utf-8', errors='replace', timeout=300, cwd=output_dir)
        if proc.returncode == 0 and os.path.exists(mcc_path):
            result["mcc_tree"] = mcc_path
            result["success"] = True
            log.emit(f"  MCC tree: {mcc_path}")
        else:
            result["error"] = f"TreeAnnotator failed (exit={proc.returncode})\n{proc.stderr[:500]}"
            log.warning(result["error"])
    except subprocess.TimeoutExpired:
        result["error"] = "TreeAnnotator timed out"
    except FileNotFoundError:
        result["error"] = f"TreeAnnotator binary not found: {treeannotator_bin}"
        log.warning(result["error"])
    except Exception as e:
        result["error"] = str(e)

    return result


# ═══════════════════════════════════════════════════════════════════
# 3. BF 迁移路径筛选 + SpreaD3 JSON
# ═══════════════════════════════════════════════════════════════════

def extract_migration_bf(
    log_file: str,
    output_dir: str,
    locations: List[str],
    burnin_pct: float = 10.0,
    bf_threshold: float = DEFAULT_BF_THRESHOLD,
    log: Optional[LogCollector] = None,
) -> Dict:
    """
    从 BEAST trace log 中提取 BSSVS indicators 的 Bayes Factor。

    BF = (posterior prob / prior prob) / ((1 - posterior prob) / (1 - prior prob))
    其中 prior prob = 期望值 (Poisson prior 下约 0.5)

    筛选 BF >= threshold 的显著迁移路径。

    Returns
    -------
    dict: {success, migrations: [{from, to, bf, posterior}], n_significant, bf_table_csv}
    """
    log = log or LogCollector()
    result = {"success": False, "migrations": [], "n_significant": 0,
              "bf_table_csv": None, "error": None}

    if not os.path.exists(log_file):
        result["error"] = f"Log file not found: {log_file}"
        return result

    # 2026-09-15 (审查 P2-27 补): 输出目录不存在时下游 to_csv 会抛错, 被外层
    # except 吞成 success=False + 与真实原因无关的报错 (migrations 其实已算对)。
    os.makedirs(output_dir, exist_ok=True)

    try:
        df = pd.read_csv(log_file, sep="\t", comment="#")
        n_total = len(df)
        burnin_idx = int(n_total * burnin_pct / 100.0)
        df_burnin = df.iloc[burnin_idx:]

        # Find indicator columns (口径全管线统一, P2-27 补)
        ind_cols = find_indicator_columns(df_burnin.columns)
        if not ind_cols:
            result["error"] = "No BSSVS indicator columns found in log"
            return result

        log.emit(f"Found {len(ind_cols)} indicator columns in log")

        # For each indicator, compute posterior mean and BF
        n_loc = len(locations)
        if n_loc < 2:
            result["error"] = "Need >=2 locations for migration analysis"
            return result

        # Build index → (from, to) mapping
        idx_map = {}
        idx = 0
        for i in range(n_loc):
            for j in range(n_loc):
                if i != j:
                    idx_map[idx] = (locations[i], locations[j])
                    idx += 1

        migrations = []
        unmapped = []
        for col in ind_cols:
            posterior = float(df_burnin[col].mean())
            # BF 公式/发散截断: 全管线唯一定义 (virphy_bridge.bayes_factor)。
            # 历史坑: 旧版 posterior>=1.0 时写 float('inf'), 下游 json.dump 产出
            # 非标准 "Infinity" 字面量 (严格 JSON 解析器直接拒绝); 且与
            # beast_parser (截断 999.0) / spread3_viz (截断 999.0) 口径不一致。
            bf = bayes_factor(posterior)

            # ── 列名 → (from, to) ────────────────────────────────────────
            # BEAST1 列名是 location.indicators.<from>.<to> (地点名, 非数字索引)。
            # BEAST2 是 location.indicators<k> (1-based 数字)。
            # 2026-09-15 (审查 P2-3): 旧版两者都解析不出来时 **默认 col_idx=0**,
            # 于是所有无法解析的列都被安到 idx_map[0] = (locations[0], locations[1])
            # —— 地图迁移连线与 heatmap 全部塌到同一对地点, 且完全静默。
            # 现改为: 解析不出就跳过该列并显式记录, 不再伪造一个 loc0→loc1。
            m_col = re.match(r'location\.indicators\.(\w+)\.(\w+)', col,
                             re.IGNORECASE)
            if m_col:
                src, dst = m_col.group(1), m_col.group(2)
            else:
                col_idx_match = re.search(r'(\d+)\s*$', col)
                col_idx = int(col_idx_match.group(1)) - 1 if col_idx_match else None
                if col_idx is not None and col_idx in idx_map:
                    src, dst = idx_map[col_idx]
                else:
                    unmapped.append(col)
                    log.warning(f"  跳过无法映射到地点对的 indicator 列: {col}")
                    continue

            migrations.append({
                "from": src,
                "to": dst,
                "bf": round(bf, 2),
                "posterior": round(posterior, 4),
                "significant": bf >= bf_threshold,
            })

        if unmapped:
            result["unmapped_indicator_columns"] = unmapped
            log.warning(f"  {len(unmapped)}/{len(ind_cols)} 个 indicator 列无法映射到"
                        f"地点对, 已排除 (不计入 n_significant)")
        if not migrations:
            result["error"] = (f"全部 {len(ind_cols)} 个 indicator 列都无法映射到地点对; "
                               f"请检查 locations 列表与 log 列名是否一致")
            log.error(result["error"])
            return result

        # Sort by BF descending
        migrations.sort(key=lambda x: -x["bf"])
        sig_migs = [m for m in migrations if m["significant"]]
        result["migrations"] = migrations
        result["n_significant"] = len(sig_migs)

        log.emit(f"Migration routes: {len(migrations)} total, "
                f"{len(sig_migs)} significant (BF >= {bf_threshold})")
        for m in sig_migs[:15]:
            log.emit(f"  {m['from']} → {m['to']}: BF={m['bf']:.1f}")

        # Save table
        bf_csv = os.path.join(output_dir, "migration_bf.csv")
        pd.DataFrame(migrations).to_csv(bf_csv, index=False)
        result["bf_table_csv"] = bf_csv

        # Generate SpreaD3 JSON
        spread3_json = _generate_spread3_json(migrations, locations, log_file, output_dir, log)
        result["spread3_json"] = spread3_json

        result["success"] = True
        return result

    except Exception as e:
        result["error"] = str(e)
        log.error(f"BF extraction failed: {e}")
        return result


def _generate_spread3_json(
    migrations: List[Dict],
    locations: List[str],
    log_file: str,
    output_dir: str,
    log: LogCollector,
) -> Optional[str]:
    """Generate SpreaD3-compatible JSON for interactive visualization.

    2026-08-30 审计修复 (BUG-1): 旧版此函数写 spread3.json (无坐标 schema B) 会覆写
    spread3_viz.generate_spread3_json 刚写下的带经纬度版本 (schema A)，而地图 HTML
    消费的是 A → 落盘文件与地图不一致且无法复现。现改写 spread3_bf.json 独立命名。
    """
    try:
        # Build nodes
        nodes = []
        for i, loc in enumerate(locations):
            nodes.append({
                "id": loc,
                "label": loc,
                "group": i,
            })

        # Build links (migration routes)
        links = []
        for m in migrations:
            if m["bf"] > 1.0:  # show all with some support
                links.append({
                    "source": locations.index(m["from"]) if m["from"] in locations else 0,
                    "target": locations.index(m["to"]) if m["to"] in locations else 0,
                    "value": m["bf"],
                    "bf": m["bf"],
                    "posterior": m["posterior"],
                })

        spread3 = {
            "nodes": nodes,
            "links": links,
            "metadata": {
                "source": os.path.basename(log_file),
                "description": "Bayesian phylogeographic migration routes (BSSVS)",
                "n_locations": len(locations),
                "n_routes": len(links),
            }
        }

        json_path = os.path.join(output_dir, "spread3_bf.json")  # 不再覆写 spread3_viz 的带坐标版
        with open(json_path, 'w') as f:
            json.dump(spread3, f, indent=2, ensure_ascii=False)
        log.emit(f"SpreaD3 JSON: {json_path}")
        return json_path

    except Exception as e:
        log.warning(f"SpreaD3 JSON failed: {e}")
        return None


# ═══════════════════════════════════════════════════════════════════
# Helpers
# ═══════════════════════════════════════════════════════════════════

def _simplify_location(loc: str) -> str:
    """Simplify location string to region level."""
    loc = loc.replace('"', '').strip()
    parts = [p.strip() for p in loc.split(',')]
    # Take the last meaningful part
    if len(parts) >= 3:
        return parts[-2]  # province/state level
    elif len(parts) == 2:
        return parts[-1]
    return parts[0]


def _xml_escape(s: str) -> str:
    """Escape special XML characters."""
    return s.replace('&', '&amp;').replace('<', '&lt;').replace('>', '&gt;') \
            .replace('"', '&quot;').replace("'", '&apos;')
