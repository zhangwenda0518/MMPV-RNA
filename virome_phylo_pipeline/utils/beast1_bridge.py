#!/usr/bin/env python3
"""
beast1_bridge.py �?BEAST 1.x 离散系统地理�?XML 生成
=====================================================
参照 ebolavirus-gl-phylogeography/xml/DISCRETE.xml 模板�?BEAST 1.x �?16 篇植物病毒贝叶斯系统地理学论文的标准格式�?
BEAST 1.x vs 2.x 关键差异:
  - 日期作为 <taxon><date value="2021.5"/></taxon> 属�?  - 离散性状�?<attr name="Location">Ningxia</attr>
  - CTMC phylogeography: generalSubstitutionModel + BSSVS
  - 树先�? coalescentLikelihood (constant/skyline/exponential)
  - 分子�? strictClockBranchRates / uncorrelatedLogNormal

BEAST 1.10.4+ 需�?
  - BEAST v1.10.4+ (https://github.com/beast-dev/beast-mcmc)
  - BEAGLE library for GPU acceleration
"""

import csv, os, re, sys, textwrap
from pathlib import Path
from typing import Dict, List, Optional
from collections import defaultdict

import numpy as np
import pandas as pd
from Bio import SeqIO

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from utils.virphy_bridge import LogCollector

# 跨管线统一 I/O 布局: ⑤ per-virus 模块目录名随 MMPV_IO_LAYOUT 解析
try:
    from mmpv_common.io_layout import ph_dir
except ImportError:  # 脱离包环境直接运行时自举
    import os as _os, sys as _sys
    _REPO_ROOT = _os.path.dirname(_os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
    if _REPO_ROOT not in _sys.path:
        _sys.path.insert(0, _REPO_ROOT)
    from mmpv_common.io_layout import ph_dir


def _to_decimal_year(date_str: str, sample: str = "") -> float:
    """Convert YYYY / YYYY-MM / YYYY-MM-DD to decimal year.

    Fail loudly: 无效日期直接抛 ValueError, 不再静默默认 2020.0
    (静默默认会污染 TMRCA/定年结果且难以追踪)。

    2026-09-15 (审查 P2-26): 实现搬到 utils/decimal_year.py —— 原先本模块用
    `(dy-1)/365` 而 beast_bridge 用"当月实际天数", phylogeo_bridge 干脆丢掉"日",
    同一条序列在三条路径上得到不同采样时间。现在全管线共用一份实现。
    """
    from utils.decimal_year import to_decimal_year
    return to_decimal_year(date_str, sample=sample, strict=True)


def _simplify_location(loc: str) -> str:
    """Simplify location to province/state level."""
    loc = loc.replace('"', '').strip()
    parts = [p.strip() for p in loc.split(',')]
    if len(parts) >= 3:
        return parts[-2]
    elif len(parts) == 2:
        return parts[-1]
    return parts[0]


def _xml_escape(s: str) -> str:
    return s.replace('&', '&amp;').replace('<', '&lt;').replace('>', '&gt;') \
            .replace('"', '&quot;')


# ══════════════════════════════════════════════════════════════════�?# BEAST 1.x Discrete Phylogeography XML Generator
# ══════════════════════════════════════════════════════════════════�?
def generate_beast1_phylogeo_xml(
    fasta_file: str,
    metadata_csv: str,
    output_dir: str,
    chain_length: int = 50_000_000,
    log_every: int = 5000,
    clock_model: str = "ucln",          # "strict" or "ucln"
    tree_prior: str = "skyline",        # "constant", "skyline" — skyline is the literature standard
    discretize_locations: bool = True,
    location_column: str = "location",
    n_location_groups: int = 12,
    substitution_model: str = "auto",     # "HKY", "GTR", or "auto" (detect from IQ-TREE log)
    gamma_categories: int = 4,
    xml_name: str = "phylogeo_beast1.xml",  # 输出 XML 文件名 (gene_dating 传 <gene>.xml)
    log: Optional[LogCollector] = None,
) -> Dict:
    """
    Generate BEAST 1.x XML for discrete phylogeography with CTMC + BSSVS.

    This replicates the standard workflow used in:
      - TuMV Silk Road (PNAS 2021): BEAST 1.10, UCLN, Skygrid, BSSVS
      - PVY global (Virus Evol 2020): BEAST 1.10, UCLN, Skyline, BSSVS
      - PVS Europe (Virology 2018): BEAST 1.8, UCLN, Skyline, BSSVS
      - RSV Japan (Viruses 2022): BEAST 1.10, UCLN, Skyline, BSSVS
      - TelMV (Phytopath Res): BEAST 1.10, UCLN, Skyline, BSSVS

    BEAST 1.x XML structure (order matters for BEAST 1.x parser):
      taxa �?alignment �?patterns �?startingTree �?treeModel
      �?clockModel �?siteModel �?coalescentLikelihood
      �?locationModel �?operators �?mcmc �?loggers

    Returns
    -------
    dict: {success, xml_path, n_taxa, n_locations, locations, year_range}
    """
    log = log or LogCollector()
    os.makedirs(output_dir, exist_ok=True)
    result = {"success": False, "xml_path": None, "n_taxa": 0,
              "n_locations": 0, "locations": [], "year_range": None, "error": None}

    try:
        # ── 1. Load alignment ──
        aln = list(SeqIO.parse(fasta_file, "fasta"))
        if len(aln) < 4:
            result["error"] = f"Need >=4 sequences, got {len(aln)}"
            return result
        seq_len = len(str(aln[0].seq))

        # ── 2. Load metadata ──
        # 纯定年模式 (discretize_locations=False): location 列可选
        meta = {}
        with open(metadata_csv, encoding='utf-8-sig') as f:
            reader = csv.DictReader(f)
            for row in reader:
                name = row.get('name', '').strip()
                date = row.get('date', '').strip()
                loc = row.get(location_column, '').strip()
                if name and date and (loc or not discretize_locations):
                    meta[name] = {"date": date, "location": loc}

        valid = [(r, meta[r.id]) for r in aln if r.id in meta]
        # 脏日期防御 (2026-09-02): 解析失败的样本剔除并警告, 不再让整个 stage 崩溃
        # (PSTVd 教训: SRR33536645 date='not applicable' 直接 ValueError 干掉 beast/phylogeo)
        clean = []
        for r, m in valid:
            try:
                _to_decimal_year(m["date"], sample=r.id)
                clean.append((r, m))
            except ValueError as e:
                log.warning(f"  ⚠️ 剔除样本 {r.id}: 脏日期 ({e})")
        dropped = len(valid) - len(clean)
        if dropped:
            log.warning(f"  ⚠️ 共剔除 {dropped} 个脏日期样本, 继续 {len(clean)} 条")
        valid = clean
        if len(valid) < 4:
            result["error"] = f"Only {len(valid)} sequences with dates+locations"
            return result

        # ── 3. Discretize locations ──
        # 注意: 地点按实际保留, 不做合并/归并/Unknown 化;
        # 稳定性由 location.rates exponential prior + 频率下界保证
        all_locs_raw = [m["location"] for _, m in valid]
        loc_counts = defaultdict(int)
        for loc in all_locs_raw:
            simplified = _simplify_location(loc)
            loc_counts[simplified] += 1

        sorted_locs = sorted(loc_counts.items(), key=lambda x: -x[1])
        if len(sorted_locs) > n_location_groups:
            top_locs = dict(sorted_locs[:n_location_groups - 1])
            top_locs["Other"] = sum(c for _, c in sorted_locs[n_location_groups - 1:])
            loc_map = {}
            for loc in all_locs_raw:
                s = _simplify_location(loc)
                loc_map[loc] = s if s in top_locs else "Other"
        else:
            top_locs = dict(sorted_locs)
            loc_map = {loc: _simplify_location(loc) for loc in all_locs_raw}

        unique_locs = sorted(top_locs.keys())
        result["locations"] = unique_locs
        result["n_locations"] = len(unique_locs)
        n_loc = len(unique_locs)

        if n_loc < 2:
            log.warning(f"  �?Only {n_loc} location(s) found �?phylogeography requires �? locations.")
            log.warning(f"     Generating clock-only BEAST XML (no BSSVS phylogeography).")
            log.warning(f"     Check your metadata CSV location column.")

        # ── 4. Years ──
        years = [_to_decimal_year(m["date"], sample=r.id) for r, m in valid]
        max_year = max(years)
        min_year = min(years)
        result["n_taxa"] = len(valid)
        result["year_range"] = (min_year, max_year)

        log.emit(f"BEAST 1.x phylogeography: {len(valid)} taxa, {seq_len} bp")
        log.emit(f"  Years: {min_year:.1f}–{max_year:.1f}")
        log.emit(f"  Locations: {n_loc} states: {unique_locs}")

        # Auto-detect substitution model from IQ-TREE log
        # NOTE: GTR 在 BEAST MCMC 中易发 Eigendecomposition 不收敛崩溃,
        # auto 模式统一用 HKY (稳健, 16篇论文常见选择)
        if substitution_model == "auto":
            # 检测 IQ-TREE log。历史坑: 只在 fasta 父目录 glob, work_dir 模式下
            # iqtree 产物在 <out>/phylogeny/ 子目录 → 检测不到 → substitution_model
            # 保持 "auto" → 走 GTR 分支 (无 kappa, 偏离文献标准 HKY)。
            iqtree_candidates = []
            base_dir = Path(fasta_file).parent
            for _depth in range(3):  # 本目录 + 上 2 级
                d = base_dir
                for _ in range(_depth):
                    d = d.parent
                iqtree_candidates += list(d.glob("*.iqtree")) + list(d.glob("*.log"))
                for _sub in (ph_dir('phylogeny'), 'tree'):
                    sd = d / _sub
                    if sd.is_dir():
                        iqtree_candidates += list(sd.glob("*.iqtree")) + list(sd.glob("*.log"))
            _seen = set()
            for cand in iqtree_candidates:
                if str(cand) in _seen:
                    continue
                _seen.add(str(cand))
                with open(cand, encoding='utf-8', errors='replace') as _f:
                    _c = _f.read()
                _m = re.search(r'Best-fit model[:\s]+(\S+)', _c)
                if _m:
                    model_str = _m.group(1)
                    # 统一用 HKY (GTR 易崩, 特征分解不收敛)
                    substitution_model = "HKY"
                    _gm = re.search(r'\+G(\d+)', model_str.upper())
                    if _gm:
                        gamma_categories = int(_gm.group(1))
                    log.emit(f"  Auto-detected model: {model_str} → HKY+G{gamma_categories} (GTR→HKY for robustness)")
                    break

        log.emit(f"  Model: {clock_model} clock + {tree_prior} coalescent + {substitution_model}")

        # ── 5. Write XML ──
        clock_param = "ucld.mean" if clock_model == "ucln" else "clock.rate"
        xml_path = os.path.join(output_dir, xml_name)
        with open(xml_path, 'w', encoding='utf-8') as f:
            _w = f.write

            _w('<?xml version="1.0" encoding="UTF-8" standalone="no"?>\n')
            _w(f'<!-- BEAST 1.x discrete phylogeography XML -->\n')
            _w(f'<!-- {len(valid)} taxa, {seq_len} bp, {n_loc} locations -->\n')
            _w(f'<!-- clock={clock_model}, prior={tree_prior}, model={substitution_model} -->\n')
            _w(f'<!-- Generated by MMPV-RNA virome_phylo_pipeline / beast1_bridge.py -->\n')
            _w('<beast version="1.10.4">\n\n')

            # ── TAXA (BEAST 1.x style: date + attr on taxon) ──
            _w(f'  <taxa id="taxa">\n')
            for rec, m in valid:
                yr = _to_decimal_year(m["date"], sample=rec.id)
                loc = loc_map.get(m["location"], "")
                _w(f'    <taxon id="{_xml_escape(rec.id)}">\n')
                _w(f'      <date value="{yr:.6f}" direction="forwards" units="years"/>\n')
                if discretize_locations and n_loc > 1:
                    _w(f'      <attr name="Location">\n        {_xml_escape(loc)}\n      </attr>\n')
                _w(f'    </taxon>\n')
            _w(f'  </taxa>\n\n')

            # ── ALIGNMENT ──
            _w(f'  <alignment id="alignment" dataType="nucleotide">\n')
            for rec, _ in valid:
                _w(f'    <sequence>\n')
                _w(f'      <taxon idref="{_xml_escape(rec.id)}"/>\n')
                _w(f'      {str(rec.seq)}\n')
                _w(f'    </sequence>\n')
            _w(f'  </alignment>\n\n')

            # ── PATTERNS ──
            _w(f'  <patterns id="patterns" from="1" every="1" strip="false">\n')
            _w(f'    <alignment idref="alignment"/>\n')
            _w(f'  </patterns>\n\n')

            # ── DEMOGRAPHIC MODEL (defined before startingTree for idref) ──
            if tree_prior == "constant":
                _w(f'  <constantSize id="demoModel" units="years">\n')
                _w(f'    <populationSize>\n')
                _w(f'      <parameter id="demoModel.popSize" value="1.0" lower="0.0"/>\n')
                _w(f'    </populationSize>\n')
                _w(f'  </constantSize>\n\n')
            elif tree_prior == "exponential":
                _w(f'  <exponentialGrowth id="demoModel" units="years">\n')
                _w(f'    <populationSize>\n')
                _w(f'      <parameter id="demoModel.popSize" value="1.0" lower="0.0"/>\n')
                _w(f'    </populationSize>\n')
                _w(f'    <growthRate>\n')
                _w(f'      <parameter id="demoModel.growthRate" value="0.0"/>\n')
                _w(f'    </growthRate>\n')
                _w(f'  </exponentialGrowth>\n\n')
            elif tree_prior == "skyline":
                # Skyline 需要起始树也有一个初始 demographic model (用 constantSize)
                _w(f'  <constantSize id="demoModel" units="years">\n')
                _w(f'    <populationSize>\n')
                _w(f'      <parameter id="demoModel.popSize" value="1.0" lower="0.0"/>\n')
                _w(f'    </populationSize>\n')
                _w(f'  </constantSize>\n\n')
            elif tree_prior == "bdsky":
                # bdsky 起始树也需要 constantSize demoModel
                _w(f'  <constantSize id="demoModel" units="years">\n')
                _w(f'    <populationSize>\n')
                _w(f'      <parameter id="demoModel.popSize" value="1.0" lower="0.0"/>\n')
                _w(f'    </populationSize>\n')
                _w(f'  </constantSize>\n\n')

            # ── STARTING TREE ──
            _w(f'  <coalescentSimulator id="startingTree">\n')
            _w(f'    <taxa idref="taxa"/>\n')
            _w(f'    <constantSize idref="demoModel"/>\n')
            _w(f'  </coalescentSimulator>\n\n')

            # ── TREE MODEL ──
            _w(f'  <treeModel id="treeModel">\n')
            _w(f'    <coalescentTree idref="startingTree"/>\n')
            _w(f'    <rootHeight>\n')
            _w(f'      <parameter id="treeModel.rootHeight"/>\n')
            _w(f'    </rootHeight>\n')
            _w(f'    <nodeHeights internalNodes="true">\n')
            _w(f'      <parameter id="treeModel.internalNodeHeights"/>\n')
            _w(f'    </nodeHeights>\n')
            _w(f'    <nodeHeights internalNodes="true" rootNode="true">\n')
            _w(f'      <parameter id="treeModel.allInternalNodeHeights"/>\n')
            _w(f'    </nodeHeights>\n')
            _w(f'  </treeModel>\n\n')

            # ── CLOCK MODEL ──
            if clock_model == "strict":
                _w(f'  <strictClockBranchRates id="default.branchRates">\n')
                _w(f'    <rate>\n')
                _w(f'      <parameter id="clock.rate" value="0.001" lower="0.0"/>\n')
                _w(f'    </rate>\n')
                _w(f'  </strictClockBranchRates>\n\n')
            else:
                _w(f'  <discretizedBranchRates id="default.branchRates">\n')
                _w(f'    <treeModel idref="treeModel"/>\n')
                _w(f'    <distribution>\n')
                _w(f'      <logNormalDistributionModel meanInRealSpace="true">\n')
                _w(f'        <mean>\n')
                _w(f'          <parameter id="ucld.mean" value="0.001" lower="0.0"/>\n')
                _w(f'        </mean>\n')
                _w(f'        <stdev>\n')
                _w(f'          <parameter id="ucld.stdev" value="0.3333" lower="0.0"/>\n')
                _w(f'        </stdev>\n')
                _w(f'      </logNormalDistributionModel>\n')
                _w(f'    </distribution>\n')
                _w(f'    <rateCategories>\n')
                _w(f'      <parameter id="default.branchRates.categories"/>\n')
                _w(f'    </rateCategories>\n')
                _w(f'  </discretizedBranchRates>\n')
                _w(f'  <rateStatistic id="default.meanRate" name="meanRate" mode="mean"\n')
                _w(f'                  internal="true" external="true">\n')
                _w(f'    <treeModel idref="treeModel"/>\n')
                _w(f'    <discretizedBranchRates idref="default.branchRates"/>\n')
                _w(f'  </rateStatistic>\n')
                _w(f'  <rateStatistic id="default.coefficientOfVariation" name="coefficientOfVariation"\n')
                _w(f'                  mode="coefficientOfVariation" internal="true" external="true">\n')
                _w(f'    <treeModel idref="treeModel"/>\n')
                _w(f'    <discretizedBranchRates idref="default.branchRates"/>\n')
                _w(f'  </rateStatistic>\n')

            # ── SITE MODEL (HKY or GTR) ──
            if substitution_model == "HKY":
                _w(f'  <hkyModel id="hky">\n')
                _w(f'    <frequencies>\n')
                _w(f'      <frequencyModel id="estimatedFreqs" dataType="nucleotide">\n')
                _w(f'        <frequencies>\n')
                _w(f'          <parameter id="frequencies" value="0.25 0.25 0.25 0.25" lower="0.01"/>\n')
                _w(f'        </frequencies>\n')
                _w(f'      </frequencyModel>\n')
                _w(f'    </frequencies>\n')
                _w(f'    <kappa>\n')
                _w(f'      <parameter id="kappa" value="2.0" lower="0.0"/>\n')
                _w(f'    </kappa>\n')
                _w(f'  </hkyModel>\n\n')
            else:
                _w(f'  <gtrModel id="gtr">\n')
                _w(f'    <frequencies>\n')
                _w(f'      <frequencyModel dataType="nucleotide">\n')
                _w(f'        <frequencies>\n')
                _w(f'          <parameter id="frequencies" value="0.25 0.25 0.25 0.25"/>\n')
                _w(f'        </frequencies>\n')
                _w(f'      </frequencyModel>\n')
                _w(f'    </frequencies>\n')
                _w(f'    <rateAC><parameter id="ac" value="1.0"/></rateAC>\n')
                _w(f'    <rateAG><parameter id="ag" value="1.0"/></rateAG>\n')
                _w(f'    <rateAT><parameter id="at" value="1.0"/></rateAT>\n')
                _w(f'    <rateCG><parameter id="cg" value="1.0"/></rateCG>\n')
                _w(f'    <rateGT><parameter id="gt" value="1.0"/></rateGT>\n')
                _w(f'  </gtrModel>\n\n')

            # Gamma site heterogeneity
            _w(f'  <siteModel id="siteModel">\n')
            _w(f'    <substitutionModel>\n')
            _w(f'      <hkyModel idref="hky"/>\n' if substitution_model == "HKY"
                else f'      <gtrModel idref="gtr"/>\n')
            _w(f'    </substitutionModel>\n')
            _w(f'    <gammaShape gammaCategories="{gamma_categories}">\n')
            _w(f'      <parameter id="alpha" value="0.5" lower="0.0"/>\n')
            _w(f'    </gammaShape>\n')
            _w(f'  </siteModel>\n\n')

            # ── TREE LIKELIHOOD ──
            _w(f'  <treeLikelihood id="treeLikelihood">\n')
            _w(f'    <patterns idref="patterns"/>\n')
            _w(f'    <treeModel idref="treeModel"/>\n')
            _w(f'    <siteModel idref="siteModel"/>\n')
            _w(f'    <strictClockBranchRates idref="default.branchRates"/>\n' if clock_model == "strict"
                else f'    <discretizedBranchRates idref="default.branchRates"/>\n')
            _w(f'  </treeLikelihood>\n\n')

            # ── COALESCENT TREE PRIOR ──
            if tree_prior == "constant":
                _w(f'  <coalescentLikelihood id="coalescent">\n')
                _w(f'    <model>\n')
                _w(f'      <constantPopulation idref="demoModel"/>\n')
                _w(f'    </model>\n')
                _w(f'    <populationTree>\n')
                _w(f'      <treeModel idref="treeModel"/>\n')
                _w(f'    </populationTree>\n')
                _w(f'  </coalescentLikelihood>\n\n')
            elif tree_prior == "exponential":
                _w(f'  <coalescentLikelihood id="coalescent">\n')
                _w(f'    <model>\n')
                _w(f'      <exponentialGrowth id="exponentialPop" units="years">\n')
                _w(f'        <populationSize>\n')
                _w(f'          <parameter id="exponential.popSize" value="1.0" lower="0.0"/>\n')
                _w(f'        </populationSize>\n')
                _w(f'        <growthRate>\n')
                _w(f'          <parameter id="exponential.growthRate" value="0.0"/>\n')
                _w(f'        </growthRate>\n')
                _w(f'      </exponentialGrowth>\n')
                _w(f'    </model>\n')
                _w(f'    <populationTree>\n')
                _w(f'      <treeModel idref="treeModel"/>\n')
                _w(f'    </populationTree>\n')
                _w(f'  </coalescentLikelihood>\n\n')
            elif tree_prior == "skyline":
                # BEAST 1.x Bayesian Skyline: generalizedSkyLineLikelihood 是自包含似然 (含 populationTree)
                # 不包进 coalescentLikelihood, 直接作为独立 likelihood
                # 历史坑: groupSize 不给 value → 默认全 0, skyline 时间轴失真。
                # 按 BEAUti 惯例初始化: sum(groupSize) = n_taxa - 1, 最近组 (group 0) 取余。
                _n_int = max(len(valid) - 1, 1)
                _g = _n_int // 5
                _r = _n_int - _g * 5
                _gsizes = [(_g + _r)] + [_g] * 4
                _gs_str = " ".join(str(x) for x in _gsizes)
                _w(f'  <generalizedSkyLineLikelihood id="skyline" linear="false">\n')
                _w(f'    <populationSizes>\n')
                _w(f'      <parameter id="skyline.popSize" dimension="5" value="1.0" lower="0.0"/>\n')
                _w(f'    </populationSizes>\n')
                _w(f'    <groupSizes>\n')
                _w(f'      <parameter id="skyline.groupSize" dimension="5" value="{_gs_str}"/>\n')
                _w(f'    </groupSizes>\n')
                _w(f'    <populationTree>\n')
                _w(f'      <treeModel idref="treeModel"/>\n')
                _w(f'    </populationTree>\n')
                _w(f'  </generalizedSkyLineLikelihood>\n\n')
                # EML prior on skyline pop sizes (smoothing)
                _w(f'  <exponentialMarkovLikelihood id="skyline.eml" jeffreys="true">\n')
                _w(f'    <chainParameter>\n')
                _w(f'      <parameter idref="skyline.popSize"/>\n')
                _w(f'    </chainParameter>\n')
                _w(f'  </exponentialMarkovLikelihood>\n\n')
                # 起始树用 coalescentSimulator + demoModel (已在前面定义)
                # skyline 不生成 coalescentLikelihood, 直接用 generalizedSkyLineLikelihood 作为似然
            elif tree_prior == "bdsky":
                # BEAST 1.x Birth-Death Skyline (BDSKY): Stadler et al. 2013
                # 关键结构 (已用小数据集验证跑通):
                #  1. birthDeathSkyline 元素: timesStartFromOrigin + units 属性必需
                #  2. 子元素: times/birthRate/deathRate/psi (dimension=n_epochs) + sampleProbability/origin (标量)
                #  3. 模型源码 bug: p.addBounds(origin.getSize()) → sampleProbability 和 origin 必须 dimension=1
                #  4. 需 speciationLikelihood 包装 (model + speciesTree), prior 引用它
                n_epochs = 5
                _w(f'  <birthDeathSkyline id="bdsky" timesStartFromOrigin="false" units="years">\n')
                _w(f'    <times>\n')
                _w(f'      <parameter id="bdsky.times" dimension="{n_epochs}" '
                    f'value="0.0 20.0 40.0 60.0 100.0"/>\n')
                _w(f'    </times>\n')
                _w(f'    <birthRate>\n')
                _w(f'      <parameter id="bdsky.lambda" dimension="{n_epochs}" '
                    f'value="{" ".join(["0.1"] * n_epochs)}"/>\n')
                _w(f'    </birthRate>\n')
                _w(f'    <deathRate>\n')
                _w(f'      <parameter id="bdsky.mu" dimension="{n_epochs}" '
                    f'value="{" ".join(["0.05"] * n_epochs)}"/>\n')
                _w(f'    </deathRate>\n')
                _w(f'    <psi>\n')
                _w(f'      <parameter id="bdsky.psi" dimension="{n_epochs}" '
                    f'value="{" ".join(["0.01"] * n_epochs)}"/>\n')
                _w(f'    </psi>\n')
                _w(f'    <sampleProbability>\n')
                _w(f'      <parameter id="bdsky.sampleProbability" value="1.0" lower="0.0" upper="1.0"/>\n')
                _w(f'    </sampleProbability>\n')
                _w(f'    <origin>\n')
                _w(f'      <parameter id="bdsky.origin" value="50.0" lower="0.0"/>\n')
                _w(f'    </origin>\n')
                _w(f'  </birthDeathSkyline>\n\n')
                # speciationLikelihood 包装 (树先验)
                _w(f'  <speciationLikelihood id="speciation">\n')
                _w(f'    <model>\n')
                _w(f'      <birthDeathSkyline idref="bdsky"/>\n')
                _w(f'    </model>\n')
                _w(f'    <speciesTree>\n')
                _w(f'      <treeModel idref="treeModel"/>\n')
                _w(f'    </speciesTree>\n')
                _w(f'  </speciationLikelihood>\n\n')

            # ── LOCATION CTMC + BSSVS (discrete phylogeography) ──
            if discretize_locations and n_loc > 1:
                _w(f'\n  <!-- ═══ Discrete Phylogeography: CTMC + BSSVS ═══ -->\n\n')
                n_rates = n_loc * (n_loc - 1)

                # General data type for location trait (with state definitions)
                _w(f'  <generalDataType id="location.dataType">\n')
                for loc in unique_locs:
                    _w(f'    <state code="{_xml_escape(loc)}"/>\n')
                _w(f'  </generalDataType>\n\n')

                # Attribute patterns for discrete trait
                _w(f'  <attributePatterns id="location.pattern" attribute="Location">\n')
                _w(f'    <taxa idref="taxa"/>\n')
                _w(f'    <generalDataType idref="location.dataType"/>\n')
                _w(f'  </attributePatterns>\n\n')

                # Frequency model (normalized from data)
                _w(f'  <frequencyModel id="location.frequencyModel" normalize="true">\n')
                _w(f'    <generalDataType idref="location.dataType"/>\n')
                _w(f'    <frequencies>\n')
                _w(f'      <parameter id="location.frequencies" dimension="{n_loc}" lower="0.01"/>\n')
                _w(f'    </frequencies>\n')
                _w(f'  </frequencyModel>\n\n')

                # CTMC substitution model for locations (BSSVS)
                _w(f'  <generalSubstitutionModel id="location.model">\n')
                _w(f'    <frequencies>\n')
                _w(f'      <frequencyModel idref="location.frequencyModel"/>\n')
                _w(f'    </frequencies>\n')
                _w(f'    <rates>\n')
                _w(f'      <parameter id="location.rates" dimension="{n_rates}" '
                    f'value="{" ".join(["1.0"] * n_rates)}" lower="0.0"/>\n')
                _w(f'    </rates>\n')
                _w(f'    <rateIndicator>\n')
                _w(f'      <parameter id="location.indicators" dimension="{n_rates}" '
                    f'value="{" ".join(["1"] * n_rates)}"/>\n')
                _w(f'    </rateIndicator>\n')
                _w(f'  </generalSubstitutionModel>\n\n')

                # Site model for location trait
                _w(f'  <siteModel id="location.siteModel">\n')
                _w(f'    <substitutionModel>\n')
                _w(f'      <generalSubstitutionModel idref="location.model"/>\n')
                _w(f'    </substitutionModel>\n')
                _w(f'  </siteModel>\n\n')

                # Branch rate model for location trait (fixed at 1.0)
                _w(f'  <strictClockBranchRates id="location.branchRates">\n')
                _w(f'    <rate>\n')
                _w(f'      <parameter id="location.clockRate" value="1.0" lower="0.0"/>\n')
                _w(f'    </rate>\n')
                _w(f'  </strictClockBranchRates>\n\n')

                # Ancestral state reconstruction (CTMC on tree, BEAUti format)
                _w(f'  <ancestralTreeLikelihood id="location.treeLikelihood"\n')
                _w(f'    stateTagName="location.states" useUniformization="true"\n')
                _w(f'    saveCompleteHistory="false" logCompleteHistory="false">\n')
                _w(f'    <attributePatterns idref="location.pattern"/>\n')
                _w(f'    <treeModel idref="treeModel"/>\n')
                _w(f'    <siteModel idref="location.siteModel"/>\n')
                _w(f'    <generalSubstitutionModel idref="location.model"/>\n')
                _w(f'    <strictClockBranchRates idref="location.branchRates"/>\n')

                # MJRM: Markov jump indicators + rewards (counts migration events)
                _w(f'    <!-- Markov jump configuration (MJRM) -->\n')
                # Indicator matrices for each migration pair
                for si, src in enumerate(unique_locs):
                    for di, dst in enumerate(unique_locs):
                        if src == dst:
                            continue
                        _w(f'    <parameter id="{_xml_escape(src)}-to-{_xml_escape(dst)}" value="\n')
                        for row in range(n_loc):
                            row_vals = ['0'] * n_loc
                            if row == si:
                                row_vals[di] = '1'
                            _w(f'      {" ".join(row_vals)}\n')
                        _w(f'    "/>\n')
                # Reward vectors for time-in-location
                _w(f'    <rewards>\n')
                for i, loc in enumerate(unique_locs):
                    vals = ['0.0'] * n_loc
                    vals[i] = '1.0'
                    _w(f'      <parameter id="{_xml_escape(loc)}_reward" value="{" ".join(vals)}" />\n')
                _w(f'    </rewards>\n')
                _w(f'  </ancestralTreeLikelihood>\n\n')

            # ── PRIORS ──
            _w(f'  <!-- Priors -->\n')

            # BSSVS statistic (defined before priors, referenced by Poisson prior)
            if discretize_locations and n_loc > 1:
                _w(f'  <sumStatistic id="location.nonZeroRates" elementwise="true">\n')
                _w(f'    <parameter idref="location.indicators"/>\n')
                _w(f'  </sumStatistic>\n\n')

            _w(f'  <prior id="prior">\n')
            # BDSKY tree prior (Speciation model - 放 prior 而非 likelihood)
            if tree_prior == "bdsky":
                _w(f'    <speciationLikelihood idref="speciation"/>\n')
            # Clock rate prior (BETS recommendation: broad CTMC-rate reference prior)
            # NOTE: bdsky 模式下不加 ctmcScalePrior — 小树高+短序列时其密度溢出为 Infinity
            # (PSTVd 359bp/树高12年实测踩坑; ucld.mean 已有 exponential prior 兜底)
            if tree_prior != "bdsky":
                _w(f'    <ctmcScalePrior>\n')
                _w(f'      <ctmcScale>\n')
                _w(f'        <parameter idref="{clock_param}"/>\n')
                _w(f'      </ctmcScale>\n')
                _w(f'      <treeModel idref="treeModel"/>\n')
                _w(f'    </ctmcScalePrior>\n')
            # Population size prior (BETS: exponential is safest)
            if tree_prior == "constant":
                _w(f'    <exponentialPrior mean="1.0" offset="0.0">\n')
                _w(f'      <parameter idref="demoModel.popSize"/>\n')
                _w(f'    </exponentialPrior>\n')
            elif tree_prior == "exponential":
                _w(f'    <exponentialPrior mean="1.0" offset="0.0">\n')
                _w(f'      <parameter idref="exponential.popSize"/>\n')
                _w(f'    </exponentialPrior>\n')
            elif tree_prior == "skyline":
                _w(f'    <exponentialPrior mean="1.0" offset="0.0">\n')
                _w(f'      <parameter idref="skyline.popSize"/>\n')
                _w(f'    </exponentialPrior>\n')
            elif tree_prior == "bdsky":
                # BDSKY priors: lambda/mu exponential, psi uniform, origin lognormal
                _w(f'    <exponentialPrior mean="0.1" offset="0.0">\n')
                _w(f'      <parameter idref="bdsky.lambda"/>\n')
                _w(f'    </exponentialPrior>\n')
                _w(f'    <exponentialPrior mean="0.05" offset="0.0">\n')
                _w(f'      <parameter idref="bdsky.mu"/>\n')
                _w(f'    </exponentialPrior>\n')
                _w(f'    <uniformPrior lower="0.0" upper="1.0">\n')
                _w(f'      <parameter idref="bdsky.psi"/>\n')
                _w(f'    </uniformPrior>\n')
                _w(f'    <logNormalPrior mean="100.0" stdev="1.0" offset="0.0">\n')
                _w(f'      <parameter idref="bdsky.origin"/>\n')
                _w(f'    </logNormalPrior>\n')
            # Kappa prior (only for HKY)
            if substitution_model == "HKY":
                _w(f'    <logNormalPrior mean="1.0" stdev="1.25" offset="0.0">\n')
                _w(f'      <parameter idref="kappa"/>\n')
                _w(f'    </logNormalPrior>\n')
            # UCLD priors
            if clock_model == "ucln":
                _w(f'    <exponentialPrior mean="0.3333" offset="0.0">\n')
                _w(f'      <parameter idref="ucld.mean"/>\n')
                _w(f'    </exponentialPrior>\n')
                _w(f'    <exponentialPrior mean="0.3333" offset="0.0">\n')
                _w(f'      <parameter idref="ucld.stdev"/>\n')
                _w(f'    </exponentialPrior>\n')
            # BSSVS priors: location rates exponential + Poisson on non-zero rates
            if discretize_locations and n_loc > 1:
                # CRITICAL: location.rates 必须有 proper prior, 否则 MCMC 将 rate 推到无穷大
                # → rate matrix 病态 → Eigendecomposition 不收敛崩溃 (PSTVd/GCVA 均踩过)
                _w(f'    <exponentialPrior mean="0.5" offset="0.0">\n')
                _w(f'      <parameter idref="location.rates"/>\n')
                _w(f'    </exponentialPrior>\n')
                _w(f'    <poissonPrior mean="0.693" offset="0.0">\n')
                _w(f'      <statistic idref="location.nonZeroRates"/>\n')
                _w(f'    </poissonPrior>\n')
            _w(f'  </prior>\n\n')

            # ── OPERATORS ──
            _w(f'  <operators id="operators" optimizationSchedule="log">\n')
            _w(f'    <scaleOperator scaleFactor="0.75" weight="3">\n')
            _w(f'      <parameter idref="{clock_param}"/>\n')
            _w(f'    </scaleOperator>\n')
            # Substitution model operators
            if substitution_model == "HKY":
                _w(f'    <scaleOperator scaleFactor="0.75" weight="1">\n')
                _w(f'      <parameter idref="kappa"/>\n')
                _w(f'    </scaleOperator>\n')
            else:
                for _r in ['ac', 'ag', 'at', 'cg', 'gt']:
                    _w(f'    <scaleOperator scaleFactor="0.75" weight="1">\n')
                    _w(f'      <parameter idref="{_r}"/>\n')
                    _w(f'    </scaleOperator>\n')
            _w(f'    <deltaExchange delta="0.01" weight="1">\n')
            _w(f'      <parameter idref="frequencies"/>\n')
            _w(f'    </deltaExchange>\n')
            _w(f'    <scaleOperator scaleFactor="0.75" weight="1">\n')
            _w(f'      <parameter idref="alpha"/>\n')
            _w(f'    </scaleOperator>\n')
            # Tree operators
            _w(f'    <subtreeSlide size="1.0" gaussian="true" weight="15">\n')
            _w(f'      <treeModel idref="treeModel"/>\n')
            _w(f'    </subtreeSlide>\n')
            _w(f'    <narrowExchange weight="15">\n')
            _w(f'      <treeModel idref="treeModel"/>\n')
            _w(f'    </narrowExchange>\n')
            _w(f'    <wideExchange weight="3">\n')
            _w(f'      <treeModel idref="treeModel"/>\n')
            _w(f'    </wideExchange>\n')
            _w(f'    <wilsonBalding weight="3">\n')
            _w(f'      <treeModel idref="treeModel"/>\n')
            _w(f'    </wilsonBalding>\n')
            # Node height operators
            _w(f'    <scaleOperator scaleFactor="0.75" weight="3">\n')
            _w(f'      <parameter idref="treeModel.rootHeight"/>\n')
            _w(f'    </scaleOperator>\n')
            _w(f'    <uniformOperator weight="30">\n')
            _w(f'      <parameter idref="treeModel.internalNodeHeights"/>\n')
            _w(f'    </uniformOperator>\n')
            # Population size operators
            if tree_prior == "constant":
                _w(f'    <scaleOperator scaleFactor="0.75" weight="3">\n')
                _w(f'      <parameter idref="demoModel.popSize"/>\n')
                _w(f'    </scaleOperator>\n')
            elif tree_prior == "skyline":
                _w(f'    <scaleOperator scaleFactor="0.75" weight="3">\n')
                _w(f'      <parameter idref="skyline.popSize"/>\n')
                _w(f'    </scaleOperator>\n')
            # UCLD operators
            if clock_model == "ucln":
                _w(f'    <scaleOperator scaleFactor="0.75" weight="3">\n')
                _w(f'      <parameter idref="ucld.mean"/>\n')
                _w(f'    </scaleOperator>\n')
                _w(f'    <scaleOperator scaleFactor="0.75" weight="3">\n')
                _w(f'      <parameter idref="ucld.stdev"/>\n')
                _w(f'    </scaleOperator>\n')
            # Up-down operator (clock rate vs tree heights)
            _w(f'    <upDownOperator scaleFactor="0.75" weight="3">\n')
            _w(f'      <up>\n')
            _w(f'        <parameter idref="{clock_param}"/>\n')
            _w(f'      </up>\n')
            _w(f'      <down>\n')
            _w(f'        <parameter idref="treeModel.allInternalNodeHeights"/>\n')
            _w(f'      </down>\n')
            _w(f'    </upDownOperator>\n')
            # BSSVS operators
            if discretize_locations and n_loc > 1:
                _w(f'    <scaleOperator scaleFactor="0.75" weight="3">\n')
                _w(f'      <parameter idref="location.rates"/>\n')
                _w(f'    </scaleOperator>\n')
                _w(f'    <bitFlipOperator weight="3">\n')
                _w(f'      <parameter idref="location.indicators"/>\n')
                _w(f'    </bitFlipOperator>\n')
            # BDSKY operators (2维参数用 randomWalk, 有界标量 origin 用 scale)
            if tree_prior == "bdsky":
                _w(f'    <randomWalkOperator windowSize="0.02" weight="3">\n')
                _w(f'      <parameter idref="bdsky.lambda"/>\n')
                _w(f'    </randomWalkOperator>\n')
                _w(f'    <randomWalkOperator windowSize="0.02" weight="3">\n')
                _w(f'      <parameter idref="bdsky.mu"/>\n')
                _w(f'    </randomWalkOperator>\n')
                _w(f'    <randomWalkOperator windowSize="0.01" weight="3">\n')
                _w(f'      <parameter idref="bdsky.psi"/>\n')
                _w(f'    </randomWalkOperator>\n')
                _w(f'    <scaleOperator scaleFactor="0.75" weight="3">\n')
                _w(f'      <parameter idref="bdsky.origin"/>\n')
                _w(f'    </scaleOperator>\n')
                _w(f'    <randomWalkOperator windowSize="2.0" weight="2">\n')
                _w(f'      <parameter idref="bdsky.times"/>\n')
                _w(f'    </randomWalkOperator>\n')
            _w(f'  </operators>\n\n')

            # ── MCMC ──
            _w(f'  <mcmc id="mcmc" chainLength="{chain_length}" autoOptimize="true"\n')
            _w(f'        operatorAnalysis="{os.path.join(output_dir, xml_name.replace(".xml", ".ops.txt"))}">\n')
            _w(f'    <posterior id="posterior">\n')
            _w(f'      <prior idref="prior"/>\n')
            _w(f'      <likelihood id="likelihood">\n')
            _w(f'        <treeLikelihood idref="treeLikelihood"/>\n')
            if tree_prior == "skyline":
                _w(f'        <generalizedSkyLineLikelihood idref="skyline"/>\n')
                _w(f'        <exponentialMarkovLikelihood idref="skyline.eml"/>\n')
            elif tree_prior != "bdsky":
                _w(f'        <coalescentLikelihood idref="coalescent"/>\n')
            if discretize_locations and n_loc > 1:
                _w(f'        <ancestralTreeLikelihood idref="location.treeLikelihood"/>\n')
            _w(f'      </likelihood>\n')
            _w(f'    </posterior>\n')
            _w(f'    <operators idref="operators"/>\n\n')

            # ── LOGGERS ──
            # Screen log
            _w(f'    <log id="screenLog" logEvery="{log_every}">\n')
            _w(f'      <column label="Posterior" dp="4" width="12">\n')
            _w(f'        <posterior idref="posterior"/>\n')
            _w(f'      </column>\n')
            _w(f'      <column label="Prior" dp="4" width="12">\n')
            _w(f'        <prior idref="prior"/>\n')
            _w(f'      </column>\n')
            _w(f'      <column label="Likelihood" dp="4" width="12">\n')
            _w(f'        <likelihood idref="likelihood"/>\n')
            _w(f'      </column>\n')
            _w(f'      <column label="ClockRate" dp="6" width="12">\n')
            _w(f'        <parameter idref="{clock_param}"/>\n')
            _w(f'      </column>\n')
            _w(f'    </log>\n\n')

            # File log
            log_path = os.path.join(output_dir, xml_name.replace('.xml', '.log'))
            _w(f'    <log id="fileLog" logEvery="{log_every}" fileName="{log_path}" overwrite="false">\n')
            _w(f'      <posterior idref="posterior"/>\n')
            _w(f'      <prior idref="prior"/>\n')
            _w(f'      <likelihood idref="likelihood"/>\n')
            _w(f'      <parameter idref="{clock_param}"/>\n')
            if substitution_model == "HKY":
                _w(f'      <parameter idref="kappa"/>\n')
            _w(f'      <parameter idref="frequencies"/>\n')
            _w(f'      <parameter idref="alpha"/>\n')
            _w(f'      <treeLikelihood idref="treeLikelihood"/>\n')
            if tree_prior == "skyline":
                _w(f'      <generalizedSkyLineLikelihood idref="skyline"/>\n')
            elif tree_prior == "bdsky":
                _w(f'      <speciationLikelihood idref="speciation"/>\n')
            else:
                _w(f'      <coalescentLikelihood idref="coalescent"/>\n')
            _w(f'      <parameter idref="treeModel.rootHeight"/>\n')
            if tree_prior == "constant":
                _w(f'      <parameter idref="demoModel.popSize"/>\n')
            elif tree_prior == "skyline":
                _w(f'      <parameter idref="skyline.popSize"/>\n')
            if clock_model == "ucln":
                _w(f'      <parameter idref="ucld.mean"/>\n')
                _w(f'      <parameter idref="ucld.stdev"/>\n')
                _w(f'      <rateStatistic idref="default.meanRate"/>\n')
                _w(f'      <rateStatistic idref="default.coefficientOfVariation"/>\n')
            if discretize_locations and n_loc > 1:
                _w(f'      <parameter idref="location.rates"/>\n')
                _w(f'      <parameter idref="location.indicators"/>\n')
                _w(f'      <ancestralTreeLikelihood idref="location.treeLikelihood"/>\n')
            _w(f'    </log>\n\n')

            # Location rate matrix log
            if discretize_locations and n_loc > 1:
                rate_log = os.path.join(output_dir, xml_name.replace('.xml', '.location.rates.log'))
                _w(f'    <log id="location.rateMatrixLog" logEvery="{log_every}"\n')
                _w(f'         fileName="{rate_log}">\n')
                _w(f'      <parameter idref="location.rates"/>\n')
                _w(f'      <parameter idref="location.indicators"/>\n')
                _w(f'    </log>\n\n')

            # Tree log
            tree_log = os.path.join(output_dir, xml_name.replace('.xml', '.trees'))
            _w(f'    <logTree id="treeLog" logEvery="{log_every}" fileName="{tree_log}"\n')
            _w(f'             nexusFormat="false" sortTranslationTable="true">\n')
            _w(f'      <treeModel idref="treeModel"/>\n')
            _w(f'      <strictClockBranchRates idref="default.branchRates"/>\n' if clock_model == "strict"
                else f'      <discretizedBranchRates idref="default.branchRates"/>\n')
            if discretize_locations and n_loc > 1:
                _w(f'      <trait name="location.states" tag="Location">\n')
                _w(f'        <ancestralTreeLikelihood idref="location.treeLikelihood"/>\n')
                _w(f'      </trait>\n')
            _w(f'    </logTree>\n')

            _w(f'  </mcmc>\n')
            _w(f'</beast>\n')

        result["xml_path"] = xml_path
        result["success"] = True
        log.emit(f"BEAST 1.x XML: {xml_path} ({os.path.getsize(xml_path)} bytes)")
        return result

    except Exception as e:
        result["error"] = str(e)
        log.error(f"BEAST 1.x XML generation failed: {e}")
        import traceback
        traceback.print_exc()
        return result
