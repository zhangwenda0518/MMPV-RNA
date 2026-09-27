#!/usr/bin/env python3
"""
path_sampling.py — Path Sampling / Stepping Stone 边际似然估计
=============================================================
实现 BEAST model-selection 包的标准 path sampling 工作流:
  - 生成 N 个 power posterior XML (beta 从 0→1)
  - 支持 BEAST 1.x 和 2.x
  - 自动收集 likelihood，计算 log marginal likelihood
  - 支持 BETS 时间信号检验 (heterochronous vs isochronous 比较)

参考文献:
  - Baele et al. (2012) Improving the accuracy of demographic and molecular clock
    model comparison while accommodating phylogenetic uncertainty. Mol Biol Evol.
  - Duchene et al. (2020) Bayesian Evaluation of Temporal Signal (BETS).
  - BETS_prior_sensitivity_ms (PLOS manuscript): 指数先验 + stepping stone 推荐

用法:
  # 单次 path sampling
  python path_sampling.py --xml phylogeo_beast1.xml --steps 20 --alpha 0.3

  # BETS: 比较有时间 vs 无时间模型
  python path_sampling.py --bets --fasta mafft.aln.fasta --meta phylo_dates.csv
"""

import os, sys, csv, re, math, glob, subprocess, time, shutil
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from collections import defaultdict

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from utils.virphy_bridge import LogCollector


# ═══════════════════════════════════════════════════════════════════
# Path Sampling XML Generator
# ═══════════════════════════════════════════════════════════════════

def generate_path_sampling_xmls(
    source_xml: str,
    output_dir: str,
    num_steps: int = 20,
    alpha: float = 0.3,
    chain_length: int = 1_000_000,
    pre_burnin: int = 100_000,
    log_every: int = 1000,
    posterior_to_prior: bool = True,
    num_replicates: int = 1,
    beast_version: str = "1",
    root_height_upper: Optional[float] = None,
    log: Optional[LogCollector] = None,
) -> Dict:
    """
    从 BEAST XML 生成 path sampling 所需的多个 power posterior XML。

    Beta 值是 power posterior 的指数:
      q(θ|β) = p(D|θ)^β × p(θ)
      当 β=0 时为 prior，β=1 时为 posterior

    Parameters
    ----------
    source_xml : str — BEAST XML 文件 (1.x 或 2.x)
    output_dir : str
    num_steps : int — power posterior 步数 (推荐 20-100)
    alpha : float — beta 分布的形状参数 (0 < α < 1，推荐 0.3)
    chain_length : int — 每步 MCMC 链长
    pre_burnin : int — 每步 burn-in
    log_every : int — 日志频率
    posterior_to_prior : bool — True: posterior→prior (β: 1→0), False: prior→posterior
    num_replicates : int — 重复次数 (推荐 1-2)
    beast_version : str — "1" or "2"

    Returns
    -------
    dict: {success, xml_dir, steps: [{beta, xml_path}], betas_file}
    """
    log = log or LogCollector()
    os.makedirs(output_dir, exist_ok=True)
    result = {"success": False, "xml_dir": output_dir, "steps": [], "error": None}

    if not os.path.exists(source_xml):
        result["error"] = f"Source XML not found: {source_xml}"
        return result

    try:
        # Read source XML
        with open(source_xml, 'r', encoding='utf-8') as f:
            xml_content = f.read()

        # Compute beta values
        betas = []
        for i in range(num_steps):
            if posterior_to_prior:
                p = (num_steps - 1 - i) / (num_steps - 1)
            else:
                p = i / (num_steps - 1)

            if p <= 0.0:
                beta = 0.0
            elif p >= 1.0:
                beta = 1.0
            else:
                beta = p ** (1.0 / alpha)
            betas.append(beta)

        log.emit(f"Path sampling: {num_steps} steps, α={alpha}")
        log.emit(f"  β range: {betas[0]:.4f} → {betas[-1]:.4f}")
        log.emit(f"  Chain: {chain_length:,} + {pre_burnin:,} burn-in")

        # Generate XMLs for each replicate × step
        # 历史坑: BEAST1 与 BEAST2 架构不同 — BEAST1 原生引擎在单一 XML 内部遍历
        # beta 网格并直接输出 log ML, 不需要每个 beta 一个 XML;
        # 逐 beta 架构只适配 BEAST2。
        for rep in range(1, num_replicates + 1):
            rep_dir = os.path.join(output_dir, f"run{rep}")
            os.makedirs(rep_dir, exist_ok=True)

            betas_file = os.path.join(rep_dir, "betas.txt")
            with open(betas_file, 'w') as bf:
                bf.write("step\tbeta\n")
                for i, beta in enumerate(betas):
                    bf.write(f"{i}\t{beta:.16f}\n")

            if str(beast_version) == "1":
                # BEAST1 原生 MLE: 单一 XML, 引擎内部遍历 beta 网格
                step_xml_path = os.path.join(rep_dir, "marginal_likelihood.xml")
                modified_xml = _inject_beast1_mle(
                    xml_content, chain_length, num_steps, alpha, pre_burnin,
                    root_height_upper=root_height_upper, log=log)
                with open(step_xml_path, 'w', encoding='utf-8') as f:
                    f.write(modified_xml)
                result["steps"].append({
                    "replicate": rep, "step": 0, "beta": 1.0,
                    "xml_path": step_xml_path, "dir": rep_dir,
                    "mode": "beast1_native_mle",
                })
                log.emit(f"  BEAST1 native MLE XML: {step_xml_path} "
                         f"({num_steps} path steps in one run)")
                continue

            for i, beta in enumerate(betas):
                step_dir = os.path.join(rep_dir, f"step{i}")
                os.makedirs(step_dir, exist_ok=True)

                step_xml_path = os.path.join(step_dir, "pathsampling.xml")

                # Modify the XML: replace <run>/<mcmc> with PathSamplingStep
                modified_xml = _inject_path_sampling(xml_content, beta,
                                                     chain_length, pre_burnin,
                                                     log_every, beast_version)

                with open(step_xml_path, 'w', encoding='utf-8') as f:
                    f.write(modified_xml)

                result["steps"].append({
                    "replicate": rep,
                    "step": i,
                    "beta": beta,
                    "xml_path": step_xml_path,
                    "dir": step_dir,
                })

        result["success"] = True
        log.emit(f"Generated {len(result['steps'])} path sampling XMLs")
        return result

    except Exception as e:
        result["error"] = str(e)
        log.error(f"Path sampling XML generation failed: {e}")
        return result


def _inject_path_sampling(
    xml_content: str,
    beta: float,
    chain_length: int,
    pre_burnin: int,
    log_every: int,
    beast_version: str,
    path_steps: int = 0,
    alpha: float = 0.3,
    root_height_upper: Optional[float] = None,
    log: Optional[LogCollector] = None,
) -> str:
    """Inject PathSamplingStep into BEAST XML, replacing the MCMC run element."""

    # BEAST 1.x: 原生 <marginalLikelihoodEstimator> + <pathLikelihood> 语法
    # (引擎 dr.inference.mcmc.MarginalLikelihoodEstimator 内置于 beast.jar,
    #  XML 语法规则已从 Parser.class 字节码完整还原, 2026-08-27):
    #   <marginalLikelihoodEstimator chainLength="N" pathSteps="N" alpha="0.3">
    #     <pathLikelihood>
    #       <source>  ...posterior compound model... </source>
    #       <destination> ...prior compound model... </destination>
    #     </pathLikelihood>
    #     <mcmc idref="mcmc"/>            ← parser 收集其 operators/loggers
    #   </marginalLikelihoodEstimator>
    # 与 BEAST2 不同: BEAST1 引擎在一个进程内遍历 beta 网格并直接输出 log ML,
    # 每个 beta 一个 XML 的架构只适配 BEAST2。
    if str(beast_version) == "1":
        return _inject_beast1_mle(xml_content, chain_length, path_steps,
                                  alpha, pre_burnin,
                                  root_height_upper=root_height_upper, log=log)

    # BEAST 2.x: replace <run id="mcmc" ...>
    else:
        likelihood_logger = (
            f'\n        <logger id="likelihoodLog" spec="Logger" '
            f'fileName="likelihood.log" logEvery="{log_every}">\n'
            f'            <log idref="likelihood"/>\n'
            f'        </logger>\n'
        )

        if 'likelihoodLog' not in xml_content:
            xml_content = xml_content.replace(
                '  </run>',
                likelihood_logger + '  </run>'
            )

        xml_content = re.sub(
            r'<run\s+id="mcmc"\s+spec="MCMC"\s+chainLength="\d+"[^>]*>',
            f'<run id="PathSamplingStep" spec="modelselection.inference.PathSamplingStep"\n'
            f'     chainLength="{chain_length}" preBurnin="{pre_burnin}" beta="{beta:.16f}">',
            xml_content
        )

    return xml_content


# ═══════════════════════════════════════════════════════════════════
# Path Sampling Result Collection
# ═══════════════════════════════════════════════════════════════════

def _inject_beast1_mle(
    xml_content: str,
    chain_length: int,
    path_steps: int,
    alpha: float,
    pre_burnin: int,
    root_height_upper: Optional[float] = None,
    log: Optional[LogCollector] = None,
) -> str:
    """BEAST1 原生 MLE 注入: 在原 XML 上添加 pathLikelihood + marginalLikelihoodEstimator。

    语法依据 (dr.inference.mcmc.MarginalLikelihoodEstimator$2 Parser 字节码还原, 2026-08-27):
      元素名: marginalLikelihoodEstimator
      必需属性: chainLength (int); 可选: pathSteps/burnin/prerun (int),
                alpha/beta (double, 默认 0.5), pathScheme (string),
                linear/lacing/spawn/printOperatorAnalysis (bool), fixedValues (double[])
      子元素: <pathLikelihood> (必需) + <mcmc> (0+) + <samplers> + MCLogger (0+)
    parser 行为: 从内层 <mcmc> 收集 OperatorSchedule 与 loggers。

    source = posterior (prior + likelihood 复合), destination = prior,
    两者都是对原 XML 已定义组件的 idref 重组合。
    注意: 原生引擎自动遍历 beta 网格并直接输出 log ML, 每个 beta 一个 XML
    的架构只适配 BEAST2, 这里产出的是单一 XML。
    """
    # 提取 posterior 块内容 → source 复合似然
    m_post = re.search(
        r'<posterior\s+id="posterior"[^>]*>(.*?)</posterior>',
        xml_content, re.DOTALL)
    if not m_post:
        raise ValueError(
            'BEAST1 MLE 注入失败: 源 XML 中找不到 <posterior id="posterior"> 块。'
            'path_sampling 只支持 beast1_bridge.py 生成的 phylogeo XML。')
    posterior_body = m_post.group(1).strip()

    # 历史坑: posterior body 内的 <likelihood id="likelihood"> 是定义块,
    # 直接复制进 source 会报 "Object with Id=likelihood already exists"
    # (2026-08-29 服务器实跑复现)。改为纯 idref 重组。
    likelihood_ids = re.findall(
        r'<(\w+Likelihood|likelihood)\s+id="([^"]+)"', posterior_body)
    if not likelihood_ids:
        raise ValueError(
            'BEAST1 MLE 注入失败: posterior 块内未找到 likelihood 定义。')
    source_refs = '\n'.join(
        f'      <{tag} idref="{lid}"/>' for tag, lid in likelihood_ids)

    if not re.search(r'<mcmc\s+id="mcmc"', xml_content):
        raise ValueError('BEAST1 MLE 注入失败: 找不到 <mcmc id="mcmc"> 块。')

    # 历史坑: gammaShape alpha 的 lower="0.0" 允许 alpha==0, MLE 在 β→0 纯先验端
    # 采样时无先验约束的参数会漂到边界 → GammaDistribution.quantile 抛
    # "Arguments out of range: t < 0" (2026-08-29 服务器实跑复现)。
    # clock.rate 等其它 lower="0.0" 参数同理, 统一顶开到 1e-3 避免退化。
    # 仅动 likelihood 链条上的参数 (clock.rate / alpha), 不动 popSize 等树先验参数。
    for _pid in ('alpha', 'clock.rate'):
        xml_content = re.sub(
            r'(<parameter\s+id="' + re.escape(_pid) +
            r'"[^>]*lower=")0\.0(")',
            r'\g<1>1e-3\g<2>', xml_content)

    # 硬边界补丁 (v6 实验废弃 2026-08-30): scaleOperator 拒绝有限上界参数
    # ("Scale operator can only be used on parameters with an infinite upper
    #  or lower bound"), BEAUti 默认 XML 无硬上界 + 全维度 working prior 锚定
    # (Duchene 原作者路线) 才是正解。仅保留 alpha/clock.rate 的 lower 顶开。

    # === 官方 GSS 三段式 (beast.community/model_selection_2, 2026-08-30 逐字对照) ===
    # 全面复查发现的 4 处偏差修正:
    #   A. prior 块必须含 <coalescentLikelihood idref="coalescent"/> (官方示例)
    #   B. source 直接 <posterior idref="posterior"/> (官方原文, 不手工 compoundLikelihood)
    #   C. estimator 的 MLE log 必须带 fileName 写文件 (官方 MLE.con.log)
    #   D. XML 末尾必须追加 steppingStoneSamplingAnalysis + pathSamplingAnalysis
    #      runnable —— 官方靠它们读 MLE log 算出并打印 log marginal likelihood;
    #      缺了它们 estimator 跑完也不出数 (此前所有"无输出"的直接原因)

    # A. prior 块补 coalescentLikelihood idref (若缺)
    # A2. 历史坑 (v4 实锢 2026-08-30): 官方 model_selection 教程反复强调
    #     "所有被估计参数必须配 proper priors, 否则 PS/SS 出数值问题"。
    #     base XML 的 prior 块只覆盖 clock.rate/popSize, GTR 五速率与 alpha 裸奔
    #     → β→0 纯先验端无约束漂移 → 似然发散 (-4500→-11700, SS 得 -1e121)。
    #     补官方 kappa 同款 logNormal(1, 1.25) 先验。
    m_prior_blk = re.search(r'<prior\s+id="prior"[^>]*>(.*?)</prior>', xml_content, re.DOTALL)
    if not m_prior_blk:
        raise ValueError('BEAST1 MLE 注入失败: 找不到 <prior id="prior"> 块。')
    _missing_priors = ''
    for _pid in ('ac', 'ag', 'at', 'cg', 'gt', 'alpha'):
        if f'idref="{_pid}"' not in m_prior_blk.group(1):
            _missing_priors += (
                f'    <logNormalPrior mean="1.0" stdev="1.25" offset="0.0" meanInRealSpace="false">\n'
                f'      <parameter idref="{_pid}"/>\n'
                '    </logNormalPrior>\n')
    if 'coalescentLikelihood idref="coalescent"' not in m_prior_blk.group(1):
        _missing_priors += '    <coalescentLikelihood idref="coalescent"/>\n'
    if _missing_priors:
        _pb = m_prior_blk.group(0)
        _pb_new = _pb.replace('</prior>', _missing_priors + '  </prior>')
        xml_content = xml_content.replace(_pb, _pb_new)

    # fileLog 文件名 + 确保所有需要 working prior 的参数都在 fileLog 里
    m_filelog = re.search(
        r'<log\s+id="fileLog"[^>]*fileName="([^"]+)"', xml_content)
    if not m_filelog:
        raise ValueError('BEAST1 MLE 注入失败: 找不到 fileLog 的 fileName。')
    log_file = m_filelog.group(1)
    # 官方模板 fileName 全是纯文件名 (5M.expolog.log), 无目录。
    # BEAST 以 XML 所在目录为 cwd 解析相对路径; 给绝对路径会被二次拼接成
    # 错误路径 "<dir><abs_path>" (2026-08-30 v3 实跑: can not be opened for
    # logTransformedNormalReferencePrior)。统一降为 basename。
    log_file = os.path.basename(log_file)
    # GTR 速率列 (ac/ag/at/cg/gt) 默认不在 fileLog, 注入补上
    filelog_body_end = xml_content.find('</log>', m_filelog.start())
    filelog_inner = xml_content[m_filelog.start():filelog_body_end]
    extra_cols = ''
    for _pid in ('ac', 'ag', 'at', 'cg', 'gt'):
        if f'<parameter idref="{_pid}"/>' not in filelog_inner:
            extra_cols += f'      <parameter idref="{_pid}"/>' + '\n'
    if extra_cols:
        xml_content = (xml_content[:filelog_body_end] + extra_cols +
                       xml_content[filelog_body_end:])

    # working prior: 正值参数用 logTransformedNormal, 频率向量用 normal
    _burnin = max(chain_length // 20, 100)  # 读 log 时弃前 95%, 官方建议后 50%
    wp_parts = []
    for _pid in ('clock.rate', 'alpha', 'ac', 'ag', 'at', 'cg', 'gt'):
        if f'<parameter id="{_pid}"' in xml_content or f'idref="{_pid}"' in xml_content:
            wp_parts.append(
                f'          <logTransformedNormalReferencePrior fileName="{log_file}"\n'
                f'                parameterColumn="{_pid}" burnin="{_burnin}">\n'
                f'            <parameter idref="{_pid}"/>\n'
                '          </logTransformedNormalReferencePrior>')
    # 历史坑 (v3c 实锢 2026-08-30): 多维参数 frequencies 不配 normalWorkingPrior ——
    # fileLog 中列名为 frequencies1-4 (逐维展开), 无裸 frequencies 列,
    # NormalKDEDistribution 报 "Column 'frequencies' can not be found"。
    # 官方模板 working prior 仅覆盖标量参数; frequencies 有 Dirichlet 边界约束,
    # 不配参考分布。
    # 【v7 重大修正 2026-08-30, Duchene BETS_code_data 原作者 XML 对照】
    # working prior 必须覆盖全部维度, 否则 β→0 端无锚:
    #   ① frequencies: logitTransformedNormalReferencePrior dimension="4"
    #      (原作者同款; v3c 的裸列名报错是因为用错了 normal 变体)
    #   ② 树维度: strictClockBranchRates idref="branchRates"
    #      (树高/拓扑的参考, 没有它 18-taxa 树在纯先验端漂移 → GTR 似然溢出
    #       double 上限, v4/v5 source 列 -1e308 实锢)
    if ('idref="frequencies"' in xml_content
            or '<parameter id="frequencies"' in xml_content):
        wp_parts.append(
            f'          <logitTransformedNormalReferencePrior fileName="{log_file}"\n'
            f'                parameterColumn="frequencies" dimension="4" burnin="{_burnin}">\n'
            '            <parameter idref="frequencies"/>\n'
            '          </logitTransformedNormalReferencePrior>')
    _m_brates = re.search(r'<strictClockBranchRates\s+id="([^"]+)"', xml_content)
    if _m_brates:
        wp_parts.append(
            f'          <strictClockBranchRates idref="{_m_brates.group(1)}"/>')
    # matching coalescent: constantSize 参数从 log 读 (官方 exponentialLogistic
    # 模板的同构做法; 不同树先验需对应替换)
    _coalescent_ref = ''
    if '<constantSize id="demoModel"' in xml_content:
        _coalescent_ref = (
            f'    <constantSize id="demoModelReference" units="years">\n'
            '      <populationSize>\n'
            f'        <parameter id="demoModelReference.popSize" parameterColumn="demoModel.popSize"\n'
            f'                      fileName="{log_file}" burnin="{_burnin}"/>\n'
            '      </populationSize>\n'
            '    </constantSize>\n'
            f'    <coalescentLikelihood id="coalescentReference">\n'
            '      <model>\n'
            '        <constantSize idref="demoModelReference"/>\n'
            '      </model>\n'
            '      <populationTree>\n'
            '        <treeModel idref="treeModel"/>\n'
            '      </populationTree>\n'
            '    </coalescentLikelihood>\n')
        wp_parts.append('          <coalescentLikelihood idref="coalescentReference"/>')

    dest_inner = (
        '        <workingPrior id="workingPrior">\n'
        + '\n'.join(wp_parts) + '\n'
        '        </workingPrior>')

    source_inner = ('      <prior idref="prior"/>\n' + source_refs)


    # 【v8 修正 2026-08-30, Duchene 逐字重读】两处硬差异:
    #   (1) treeModel.rootHeight 必须在 prior 内有 uniformPrior 硬界
    #      (原作者 cholera years 尺度用 0~300; 无界时 beta->0 端树高自由膨胀,
    #       coalescent+GTR 似然一起爆 -> v7b source/destination 两列同爆 -1e308)
    #   (2) prior 块补 strictClockBranchRates idref (速率维度 prior 贡献)
    #
    # 2026-09-15 (审查 P2-5) 修复: 旧版把 _rh_upper 硬编码成 300.0
    # (注释自承"源自 ice_viruses 取样 1947-2014")。这个上界对时间深度 >300 年
    # 的数据集是**硬截断 TMRCA** —— rootHeight 被卡在 300 以内, 定年结果系统性偏年轻,
    # 而且不报错。
    # 现改为从 XML 里的 taxon 采样日期反推数据跨度 span, 取 5×span 作上界
    # (对 ice_viruses 的 67 年跨度 → 335 年, 与旧值同量级; 对深时间数据自动放宽)。
    # 需要固定值时显式传 root_height_upper。
    if root_height_upper is None:
        _dates = [float(x) for x in
                  re.findall(r'<date\s+value="([-+0-9.eE]+)"', xml_content)]
        _span = (max(_dates) - min(_dates)) if len(_dates) >= 2 else 0.0
        # 下界 10.0 防止 span≈0 (同日采样) 时上界退化到 0
        _rh_upper = max(10.0, 5.0 * _span)
        _rh_src = f"数据跨度 {_span:g} 年 × 5"
    else:
        _rh_upper = float(root_height_upper)
        _rh_src = "调用方指定"
    if log is not None:
        log.emit(f"  [path_sampling] rootHeight uniform 上界 = {_rh_upper:g} ({_rh_src})")
    if 'treeModel.rootHeight' in xml_content:
        _m_pr = re.search(r'<prior\s+id="prior"[^>]*>(.*?)</prior>', xml_content, re.DOTALL)
        if _m_pr and 'treeModel.rootHeight' not in _m_pr.group(1):
            _add_rh = (
                f'    <uniformPrior lower="0.0" upper="{_rh_upper}">' + chr(10) +
                '      <parameter idref="treeModel.rootHeight"/>' + chr(10) +
                '    </uniformPrior>' + chr(10))
            _pb0 = _m_pr.group(0)
            xml_content = xml_content.replace(
                _pb0, _pb0.replace('</prior>', _add_rh + '  </prior>'))
    _m_br2 = re.search(r'<strictClockBranchRates\s+id="([^"]+)"', xml_content)
    if _m_br2:
        _m_pr2 = re.search(r'<prior\s+id="prior"[^>]*>(.*?)</prior>', xml_content, re.DOTALL)
        if _m_pr2 and f'strictClockBranchRates idref="{_m_br2.group(1)}"' not in _m_pr2.group(1):
            _pb1 = _m_pr2.group(0)
            xml_content = xml_content.replace(
                _pb1, _pb1.replace(
                    '</prior>',
                    f'    <strictClockBranchRates idref="{_m_br2.group(1)}"/>' + chr(10) + '  </prior>'))

    # 历史坑: pathLikelihood 必须是 marginalLikelihoodEstimator 的子元素
    # (兄弟位置报 "Exactly one ELEMENT of type PathLikelihood REQUIRED",
    # 2026-08-29 实跑; 之前报 samplers 错误是校验顺序掩盖了这一点)。
    # 历史坑: GSS 必须配 pseudo-prior (BEAUti 的 "Tree working prior: Matching
    # coalescent model"), 否则 β→0 端参数无约束漂移 → GammaDistribution.quantile
    # 抛 "Arguments out of range: t < 0" 带伤运行 (2026-08-29 实跑 792 次异常,
    # 参考 Raindy BETS 教程确认官方路径)。PS (非 GSS) 不需要。
    # sourcePseudoPrior = prior 原样复制 (树先验+参数先验都在里面, 起 β 端锚定作用)。
    _m_prior = re.search(
        r'(<prior id="prior".*?</prior>)', xml_content, re.DOTALL)
    pseudo_block = ''
    if _m_prior:
        pseudo_block = (
            '      <sourcePseudoPrior>\n'
            '        <compoundLikelihood id="mle.sourcePseudoPrior">\n'
            '          <prior idref="prior"/>\n'
            '        </compoundLikelihood>\n'
            '      </sourcePseudoPrior>\n'
        )

    # B. source = 官方原文写法: 直接 posterior idref
    # C. MLE log 带 fileName; D. 末尾追加 steppingStone/pathSampling analysis runnable
    mle_log_file = log_file + '.mle'
    mle_base = os.path.basename(mle_log_file)

    mle_block = (
        '\n  <!-- BEAST1 native marginal likelihood estimation (GSS, official 3-stage) -->\n'
        f'{_coalescent_ref}'
        f'  <marginalLikelihoodEstimator id="mle" chainLength="{chain_length}"\n'
        f'        pathSteps="{path_steps}" alpha="{alpha}" burnin="{pre_burnin}"\n'
        f'        pathScheme="betaQuantile" printOperatorAnalysis="false">\n'
        '    <pathLikelihood id="pathLikelihood">\n'
        '      <source>\n'
        '        <posterior idref="posterior"/>\n'
        '      </source>\n'
        '      <destination>\n'
        f'{dest_inner}\n'
        '      </destination>\n'
        '    </pathLikelihood>\n'
        '    <samplers>\n'
        '      <mcmc idref="mcmc"/>\n'
        '    </samplers>\n'
        f'    <log id="MLE.con" logEvery="{max(chain_length // 10, 100)}" fileName="{mle_base}">\n'
        '      <pathLikelihood idref="pathLikelihood"/>\n'
        '    </log>\n'
        '  </marginalLikelihoodEstimator>\n\n'
        f'  <steppingStoneSamplingAnalysis fileName="{mle_base}">\n'
        '    <likelihoodColumn name="pathLikelihood.delta"/>\n'
        '    <thetaColumn name="pathLikelihood.theta"/>\n'
        '  </steppingStoneSamplingAnalysis>\n'
        f'  <pathSamplingAnalysis fileName="{mle_base}">\n'
        '    <likelihoodColumn name="pathLikelihood.delta"/>\n'
        '    <thetaColumn name="pathLikelihood.theta"/>\n'
        '  </pathSamplingAnalysis>\n'
    )
    return xml_content.replace('</beast>', mle_block + '</beast>')


def _collect_beast1_mle_results(ps_dir: str, log: Optional[LogCollector] = None) -> Dict:
    """收集 BEAST1 原生 MLE 结果。

    引擎输出位置: 每次 run 的 stdout 及 run 目录下的 .mle.out (由 runner 重定向);
    MarginalLikelihoodEstimator 输出格式 (字节码还原):
      "Path Sampling Marginal Likelihood Estimator:\n\tEquilibrating chain ...\n"
      以及 integrator 最终的 "log marginal likelihood = <value>" 行。
    数值行由 runner 端 (run_bets_test / 用户脚本) 重定向捕获; 这里读所有 *.out / stdout.log。
    """
    log = log or LogCollector()
    result = {"success": False, "log_marginal_likelihood": None,
              "log_marginal_likelihood_ss": None, "steps": [], "error": None,
              "mode": "beast1_native_mle"}

    ml_values = []
    for root, _dirs, files in os.walk(ps_dir):
        for fn in files:
            if not (fn.endswith('.out') or fn.endswith('.stdout.log') or fn == 'mle_output.txt'):
                continue
            path = os.path.join(root, fn)
            try:
                with open(path, encoding='utf-8', errors='replace') as f:
                    content = f.read()
            except OSError:
                continue
            # 典型行: "log marginal likelihood = -1234.5678" 或
            #          "Marginal likelihood: -1234.5678"
            for m in re.finditer(
                    r'log marginal likelihood\s*[=:]\s*(-?\d+\.?\d*(?:[eE][-+]?\d+)?)',
                    content, re.IGNORECASE):
                ml_values.append((path, float(m.group(1))))
            for m in re.finditer(
                    r'Marginal likelihood:\s*(-?\d+\.?\d*(?:[eE][-+]?\d+)?)',
                    content):
                ml_values.append((path, float(m.group(1))))

    if not ml_values:
        result["error"] = (
            "BEAST1 native MLE: 未找到 log marginal likelihood 输出。"
            "确认 run 已完成且 stdout 已重定向到 run 目录下的 .out 文件。")
        log.error(result["error"])
        return result

    # 多个 run 取均值, 记录来源
    vals = [v for _, v in ml_values]
    result["log_marginal_likelihood"] = round(float(np.mean(vals)), 4)
    result["steps"] = [{"source": p, "log_ml": v} for p, v in ml_values]
    result["success"] = True
    log.emit(f"BEAST1 native MLE: log ML = {result['log_marginal_likelihood']} "
             f"(from {len(vals)} run(s))")
    return result


def collect_path_sampling_results(
    ps_dir: str,
    burnin_pct: float = 50.0,
    alpha: float = 0.3,
    log: Optional[LogCollector] = None,
) -> Dict:
    """
    收集 path sampling 各步的 likelihood，计算 log marginal likelihood。

    两个估计量 (都要求 β **升序** 积分):

      Thermodynamic Integration (TI, 梯形法):
        log ML = ∫₀¹ E_β[log p(D|θ)] dβ
               ≈ Σ_k (β_{k+1} - β_k) · (f_k + f_{k+1}) / 2,  f_k = E_{β_k}[log p(D|θ)]

      Stepping Stone (SS, Xie et al. 2011):
        log ML ≈ Σ_k (β_{k+1} - β_k) · log( (1/n) Σ_i p(D|θ_{k,i})^{β_{k+1}} )

    `likelihood.log` 里的 `likelihood` 列是**未加幂**的 log p(D|θ)（BEAST2
    model-selection 的 PathSampleAnalyser 同样这样解读），β 幂次由本函数施加。

    Returns
    -------
    dict: {
        success, log_marginal_likelihood, log_marginal_likelihood_ss,
        steps: [{step, beta, mean_likelihood, n_replicates,
                 ti_contribution, ss_contribution}],
        summary_csv
    }
    """
    log = log or LogCollector()
    result = {"success": False, "log_marginal_likelihood": None,
              "log_marginal_likelihood_ss": None, "steps": [], "error": None}

    try:
        # BEAST1 原生 MLE: 引擎直接在 stdout/日志输出 log ML,
        # 不产出逐 beta likelihood.log。先探测 marginal_likelihood.xml 模式。
        _mle_xmls = []
        for _root, _dirs, _files in os.walk(ps_dir):
            if 'marginal_likelihood.xml' in _files:
                _mle_xmls.append(_root)
        if _mle_xmls:
            return _collect_beast1_mle_results(ps_dir, log=log)

        # Find all likelihood.log files
        likelihood_files = []
        for root, dirs, files in os.walk(ps_dir):
            if 'likelihood.log' in files:
                step_match = re.search(r'step(\d+)', root)
                run_match = re.search(r'run(\d+)', root)
                if step_match:
                    likelihood_files.append({
                        "path": os.path.join(root, 'likelihood.log'),
                        "step": int(step_match.group(1)),
                        "run": int(run_match.group(1)) if run_match else 1,
                        "dir": root,
                    })

        if not likelihood_files:
            result["error"] = "No likelihood.log files found"
            return result

        # Group by step, average across replicates
        step_data = defaultdict(list)

        for lf in likelihood_files:
            try:
                df = pd.read_csv(lf["path"], sep="\t", comment="#")
                if len(df) == 0:
                    continue
                burnin_idx = int(len(df) * burnin_pct / 100.0)
                df_burnin = df.iloc[burnin_idx:]

                if 'likelihood' in df_burnin.columns:
                    mean_ll = df_burnin['likelihood'].mean()
                elif 'Likelihood' in df_burnin.columns:
                    mean_ll = df_burnin['Likelihood'].mean()
                else:
                    # Try first numeric column
                    numeric_cols = df_burnin.select_dtypes(include=[np.number]).columns
                    mean_ll = df_burnin[numeric_cols[0]].mean() if len(numeric_cols) > 0 else 0.0

                step_data[lf["step"]].append(mean_ll)
            except Exception:
                pass

        if not step_data:
            result["error"] = "Could not parse any likelihood files"
            return result

        # ── β 网格 ──
        # 历史坑: collect 阶段曾用默认 alpha 重新算 beta, 若 generate 用了不同 alpha 或
        # posterior_to_prior=False 则 beta 映射错位; 某步 likelihood 缺失时 max_step 压缩,
        # 整条 beta 轴再错位。优先读 generate 落盘的 betas.txt。
        # 2026-09-15 修复 (P1-10):
        #   · 原 `os.path.join(output_dir, "run1", "betas.txt")` 里 output_dir 是
        #     **未定义变量** (函数形参叫 ps_dir) -> 该分支必然 NameError;
        #   · 只在 run1 找 betas.txt, num_replicates>1 时其余 replicate 的 β 轴无从核对;
        #   · 找不到 betas.txt 时的兜底网格假设 "step 序号递增 => β 递增", 而
        #     generate_path_sampling_xmls 默认 posterior_to_prior=True (β: 1→0),
        #     方向正好相反 —— 见下面"按 β 升序积分"的说明。
        betas = {}
        _betas_file = None
        for _rep in sorted(glob.glob(os.path.join(ps_dir, "run*"))):
            _cand = os.path.join(_rep, "betas.txt")
            if os.path.exists(_cand):
                _betas_file = _cand
                break
        if _betas_file:
            try:
                with open(_betas_file, encoding='utf-8') as _bf:
                    next(_bf)                       # 跳过 "step\tbeta" 表头
                    for _line in _bf:
                        _parts = _line.strip().split('\t')
                        if len(_parts) == 2:
                            betas[int(_parts[0])] = float(_parts[1])
            except (OSError, ValueError) as _e:
                log.warning(f"  读取 betas.txt 失败 ({_e}), 回退内部网格")
                betas = {}
        if betas:
            log.emit(f"  β 网格来自 {os.path.basename(os.path.dirname(_betas_file))}"
                     f"/betas.txt ({len(betas)} 值, "
                     f"β={min(betas.values()):.4g}→{max(betas.values()):.4g})")
        else:
            _mx = max(step_data.keys())
            for s in sorted(step_data.keys()):
                p = s / _mx if _mx > 0 else 0.0
                betas[s] = p ** (1.0 / alpha) if p > 0 else 0.0
            log.warning("  未找到 betas.txt, 用内部 alpha 重建 β 网格 —— "
                        "若 generate 用了不同 alpha / posterior_to_prior, 结果会偏")

        # ── 积分网格: 必须按 β 升序 ──
        # 2026-09-15 修复 (P1-10, 核心): TI 的定义是 log ML = ∫₀¹ E_β[log p(D|θ)] dβ,
        # 梯形展开为 Σ (β_{k+1} - β_k)·(f_k + f_{k+1})/2, **要求 β 升序**。
        # 旧代码按 step 序号遍历, 而 betas.txt 默认 posterior_to_prior=True (β: 1→0),
        # 于是每个 (β_{k+1} - β_k) 都是负数 -> 积分整体**符号反转**, 得到 -log ML。
        # 现显式按 β 升序排序后再积分, 与 β 轴的书写方向解耦。
        grid = []
        for s in sorted(step_data.keys()):
            grid.append((betas.get(s, 0.0), s,
                         float(np.mean(step_data[s])), len(step_data[s])))
        grid.sort(key=lambda t: t[0])

        log_ml = 0.0
        log_ml_ss = 0.0
        steps_out = []

        for i, (beta_val, s, mean_ll, n_rep) in enumerate(grid):
            ti_contrib = 0.0
            ss_contrib = 0.0
            if i > 0:
                prev_beta, _prev_s, prev_ll, _prev_n = grid[i - 1]
                d_beta = beta_val - prev_beta
                # Thermodynamic integration (trapezoidal rule)
                ti_contrib = d_beta * (mean_ll + prev_ll) / 2.0
                log_ml += ti_contrib
                # Stepping stone (Xie et al. 2011):
                #   Σ_k (β_{k+1} - β_k) · log( (1/n) Σ_i p(D|θ_{k,i})^{β_{k+1}} )
                # 2026-09-15 修复 (P1-10):
                #   · 旧代码 `np.exp(step_data[s])` 对**对数似然**直接取指数 ——
                #     真实数据 log L ≈ -1e3 ~ -1e5, exp 必然下溢为 0.0,
                #     紧接着 math.log(0.0) 抛 "math domain error" 让整个函数失败;
                #   · 且漏掉 β 幂次 (应为 p^{β_{k+1}}, 不是 p^1);
                #   · step_data[s] 是各 replicate 的均值列表, 对它再取 exp 语义也不对。
                #   现用 logsumexp 稳定形式, 并把 β 幂次放回去。
                _lls = np.asarray(step_data[s], dtype=float) * beta_val
                if _lls.size:
                    _m = float(_lls.max())
                    log_ratio = _m + math.log(float(np.mean(np.exp(_lls - _m))))
                else:
                    log_ratio = 0.0
                ss_contrib = d_beta * log_ratio
                log_ml_ss += ss_contrib

            steps_out.append({
                "step": s, "beta": round(beta_val, 6),
                "mean_likelihood": round(mean_ll, 4),
                "n_replicates": n_rep,
                "ti_contribution": round(float(ti_contrib), 6),
                "ss_contribution": round(float(ss_contrib), 6),
            })

        # 网格按 β 升序输出, 便于人工核对 (step 列保留原始步号)
        result["log_marginal_likelihood"] = round(log_ml, 4)
        result["log_marginal_likelihood_ss"] = round(log_ml_ss, 4)
        result["steps"] = steps_out

        # Save summary
        summary_csv = os.path.join(ps_dir, "path_sampling_summary.csv")
        pd.DataFrame(steps_out).to_csv(summary_csv, index=False)
        result["summary_csv"] = summary_csv

        log.emit(f"Log marginal likelihood (TI): {log_ml:.4f}")
        log.emit(f"Log marginal likelihood (SS): {log_ml_ss:.4f}")
        result["success"] = True
        return result

    except Exception as e:
        result["error"] = str(e)
        log.error(f"Path sampling collection failed: {e}")
        return result


# ═══════════════════════════════════════════════════════════════════
# BETS: Bayesian Evaluation of Temporal Signal
# ═══════════════════════════════════════════════════════════════════

def run_beast1_mle(
    xml_path: str,
    beast_bin: str = "beast",
    threads: int = 4,
    timeout: int = 3600,
    log: Optional[LogCollector] = None,
    seed: Optional[int] = None,
) -> Dict:
    """执行单个 BEAST1 原生 MLE XML, stdout 重定向到 run 目录 mle_output.txt。

    历史坑: MarginalLikelihoodEstimator 的 log ML 结果只输出到 stdout
    (不写 .log 文件), 必须重定向才能被 _collect_beast1_mle_results 读到。
    可复现性 (2026-09-15 修复): 补 -seed (旧版不传 → BEAST 随机种子, log ML 不可复现)。
    """
    log = log or LogCollector()
    result = {"success": False, "log_ml": None, "stdout_file": None, "error": None}
    run_dir = os.path.dirname(xml_path)
    out_file = os.path.join(run_dir, "mle_output.txt")
    if seed is None:
        try:
            seed = int(os.environ.get("PHYLO_SEED_BASE", "") or 20260915)
        except ValueError:
            seed = 20260915
    cmd = [beast_bin, "-threads", str(threads), "-seed", str(int(seed)), xml_path]
    result["seed"] = int(seed)
    log.emit(f"  BEAST1 MLE: {' '.join(cmd)}")
    try:
        with open(out_file, 'w', encoding='utf-8', errors='replace') as out_f:
            proc = subprocess.run(cmd, stdout=out_f, stderr=subprocess.STDOUT,
                                  timeout=timeout, cwd=run_dir)
        if proc.returncode != 0:
            result["error"] = f"beast exit code {proc.returncode}, see {out_file}"
            return result
    except subprocess.TimeoutExpired:
        result["error"] = f"beast timeout after {timeout}s"
        return result
    except FileNotFoundError:
        result["error"] = f"beast binary not found: {beast_bin}"
        return result

    coll = _collect_beast1_mle_results(run_dir, log=log)
    result["stdout_file"] = out_file
    if coll["success"]:
        result["success"] = True
        result["log_ml"] = coll["log_marginal_likelihood"]
    else:
        result["error"] = coll["error"]
    return result


def run_bets_test(
    fasta_file: str,
    metadata_csv: str,
    output_dir: str,
    beast_version: str = "1",
    clock_model: str = "ucln",
    tree_prior: str = "constant",
    path_steps: int = 20,
    path_alpha: float = 0.3,
    path_chain: int = 1_000_000,
    execute: bool = False,
    beast_bin: str = "beast",
    threads: int = 4,
    mle_timeout: int = 7200,
    log: Optional[LogCollector] = None,
) -> Dict:
    """
    BETS: Bayesian Evaluation of Temporal Signal.

    比较两个模型的 log marginal likelihood:
      1. Heterochronous (有时间信息) — 采样时间不同
      2. Isochronous (无时间信息) — 所有样本设相同日期

    log BF > 5  → very strong evidence for temporal signal
    log BF > 3  → strong evidence
    log BF < -3 → strong evidence against
    中间值       → inconclusive

    References:
      Duchene et al. (2020) BETS. Mol Biol Evol.
      BETS_prior_sensitivity_ms: exponential popSize prior recommended.

    Returns
    -------
    dict: {
        success, log_bf, has_temporal_signal, conclusion,
        heterochronous_ml, isochronous_ml,
        hetero_xml, iso_xml, hetero_ps_dir, iso_ps_dir
    }
    """
    log = log or LogCollector()
    result = {"success": False, "log_bf": None, "has_temporal_signal": None,
              "conclusion": "", "error": None}

    os.makedirs(output_dir, exist_ok=True)

    try:
        from utils.beast1_bridge import generate_beast1_phylogeo_xml
        from utils.phylogeo_bridge import generate_phylogeo_xml

        # ── Step 1: Generate heterochronous XML ──
        log.emit("─" * 50)
        log.emit("BETS Step 1: Heterochronous model (with sampling times)")

        hetero_dir = os.path.join(output_dir, "heterochronous")
        if str(beast_version) == "1":
            hetero_xml_r = generate_beast1_phylogeo_xml(
                fasta_file=fasta_file, metadata_csv=metadata_csv,
                output_dir=hetero_dir, clock_model=clock_model,
                tree_prior=tree_prior, discretize_locations=False,
                log=log,
            )
        else:
            hetero_xml_r = generate_phylogeo_xml(
                fasta_file=fasta_file, metadata_csv=metadata_csv,
                output_dir=hetero_dir, clock_model=clock_model,
                tree_prior=tree_prior, discretize_locations=False,
                log=log,
            )

        if not hetero_xml_r["success"]:
            result["error"] = f"Heterochronous XML failed: {hetero_xml_r.get('error')}"
            return result

        hetero_xml = hetero_xml_r["xml_path"]
        result["hetero_xml"] = hetero_xml

        # ── Step 2: Generate isochronous XML (all dates = max_date) ──
        log.emit("─" * 50)
        log.emit("BETS Step 2: Isochronous model (without sampling times)")

        # Read metadata, set all dates to same value
        iso_dir = os.path.join(output_dir, "isochronous")
        os.makedirs(iso_dir, exist_ok=True)

        df = pd.read_csv(metadata_csv)
        date_col = df.columns[1]
        max_date = df[date_col].max() if len(df) > 0 else "2020"
        df_iso = df.copy()
        df_iso[date_col] = max_date
        iso_csv = os.path.join(iso_dir, "isochronous_dates.csv")
        df_iso.to_csv(iso_csv, index=False)

        if str(beast_version) == "1":
            iso_xml_r = generate_beast1_phylogeo_xml(
                fasta_file=fasta_file, metadata_csv=iso_csv,
                output_dir=iso_dir, clock_model=clock_model,
                tree_prior=tree_prior, discretize_locations=False,
                log=log,
            )
        else:
            iso_xml_r = generate_phylogeo_xml(
                fasta_file=fasta_file, metadata_csv=iso_csv,
                output_dir=iso_dir, clock_model=clock_model,
                tree_prior=tree_prior, discretize_locations=False,
                log=log,
            )

        if not iso_xml_r["success"]:
            result["error"] = f"Isochronous XML failed: {iso_xml_r.get('error')}"
            return result

        iso_xml = iso_xml_r["xml_path"]
        result["iso_xml"] = iso_xml

        # ── Step 3: Generate path sampling XMLs for both ──
        log.emit("─" * 50)
        log.emit("BETS Step 3: Path sampling XMLs")

        hetero_ps_dir = os.path.join(hetero_dir, "path_sampling")
        iso_ps_dir = os.path.join(iso_dir, "path_sampling")

        ps_hetero = generate_path_sampling_xmls(
            source_xml=hetero_xml, output_dir=hetero_ps_dir,
            num_steps=path_steps, alpha=path_alpha,
            chain_length=path_chain, beast_version=beast_version,
            log=log,
        )
        ps_iso = generate_path_sampling_xmls(
            source_xml=iso_xml, output_dir=iso_ps_dir,
            num_steps=path_steps, alpha=path_alpha,
            chain_length=path_chain, beast_version=beast_version,
            log=log,
        )

        result["hetero_ps_dir"] = hetero_ps_dir
        result["iso_ps_dir"] = iso_ps_dir

        # ── Step 4: Instructions for running ──
        log.emit("─" * 50)
        log.emit("BETS Step 4: Run path sampling MCMC")
        log.emit(f"  Heterochronous: {hetero_ps_dir}/run1/step*/pathsampling.xml")
        log.emit(f"  Isochronous:   {iso_ps_dir}/run1/step*/pathsampling.xml")
        log.emit(f"\n  After MCMC completes, run:")
        log.emit(f"  python path_sampling.py --collect {hetero_ps_dir}")
        log.emit(f"  python path_sampling.py --collect {iso_ps_dir}")

        # ── Step 5: Run + collect ──
        # execute=True 时直接跑 BEAST (BEAST1 原生 MLE: 每模型单一 XML);
        # 否则只收集已有结果 (外部自行跑 MCMC 后调用)。
        hetero_ml = collect_path_sampling_results(hetero_ps_dir, log=log)
        iso_ml = collect_path_sampling_results(iso_ps_dir, log=log)

        if execute:
            # BEAST1 原生 MLE 模式: 找 marginal_likelihood.xml 并执行
            for label, ps_dir in [("heterochronous", hetero_ps_dir),
                                  ("isochronous", iso_ps_dir)]:
                run_dir = os.path.join(ps_dir, "run1")
                xml_f = os.path.join(run_dir, "marginal_likelihood.xml")
                if os.path.exists(xml_f):
                    log.emit(f"  Running BEAST1 MLE ({label})...")
                    run_beast1_mle(xml_f, beast_bin=beast_bin, threads=threads,
                                   timeout=mle_timeout, log=log)
            # 重收集 (现在有 mle_output.txt 了)
            hetero_ml = collect_path_sampling_results(hetero_ps_dir, log=log)
            iso_ml = collect_path_sampling_results(iso_ps_dir, log=log)

        if hetero_ml["success"] and iso_ml["success"]:
            log_bf = (hetero_ml["log_marginal_likelihood"] or 0) - \
                     (iso_ml["log_marginal_likelihood"] or 0)
            result["log_bf"] = round(log_bf, 4)
            result["heterochronous_ml"] = hetero_ml["log_marginal_likelihood"]
            result["isochronous_ml"] = iso_ml["log_marginal_likelihood"]

            if log_bf >= 5:
                result["has_temporal_signal"] = True
                result["conclusion"] = (
                    f"✓ Very strong temporal signal (log BF={log_bf:.1f} ≥ 5). "
                    f"Proceed with molecular dating."
                )
            elif log_bf >= 3:
                result["has_temporal_signal"] = True
                result["conclusion"] = (
                    f"✓ Strong temporal signal (log BF={log_bf:.1f} ≥ 3). "
                    f"Molecular dating is supported."
                )
            elif log_bf <= -3:
                result["has_temporal_signal"] = False
                result["conclusion"] = (
                    f"✗ Strong evidence against temporal signal (log BF={log_bf:.1f} ≤ -3). "
                    f"Do NOT perform molecular dating."
                )
            else:
                result["has_temporal_signal"] = None
                result["conclusion"] = (
                    f"⚠ Inconclusive (log BF={log_bf:.1f}). "
                    f"Increase path sampling steps or chain length."
                )
            log.emit(f"\n  {result['conclusion']}")
        else:
            # 2026-09-15 (审查 P2-10) 修复: 旧版无论两个 ML 是否算出, 一律
            # result["success"]=True —— 调用方拿到 success=True + log_bf=None,
            # 分不清"log BF 不显著"和"根本没算出来", 会把一个空转的 BETS 当成功
            # (并可能被 checkpoint 记为完成)。
            _why = []
            if not hetero_ml.get("success"):
                _why.append(f"heterochronous ML 失败 ({hetero_ml.get('error') or '未解析到结果'})")
            if not iso_ml.get("success"):
                _why.append(f"isochronous ML 失败 ({iso_ml.get('error') or '未解析到结果'})")
            result["error"] = "BETS 未完成: " + "; ".join(_why)
            result["heterochronous_ml"] = hetero_ml.get("log_marginal_likelihood")
            result["isochronous_ml"] = iso_ml.get("log_marginal_likelihood")
            result["conclusion"] = ("⚠ BETS 未完成 —— 未得到 log BF, 不能据此判断时间信号 "
                                    "(这不是'无时间信号')")
            log.error(f"  {result['error']}")
            log.error(f"  {result['conclusion']}")
            return result

        result["success"] = True
        return result

    except Exception as e:
        result["error"] = str(e)
        log.error(f"BETS failed: {e}")
        import traceback
        traceback.print_exc()
        return result


# ═══════════════════════════════════════════════════════════════════
# BEAST Job Manager (SLURM/PBS/本地)
# ═══════════════════════════════════════════════════════════════════

def generate_slurm_submit_script(
    xml_path: str,
    output_dir: str,
    job_name: str = "beast_job",
    beast_bin: str = "beast",
    threads: int = 8,
    time_limit: str = "48:00:00",
    memory: str = "16G",
    partition: str = "normal",
    modules: List[str] = None,
    log: Optional[LogCollector] = None,
) -> str:
    """
    Generate SLURM submission script for BEAST job.

    Returns path to the .sh script.
    """
    log = log or LogCollector()
    os.makedirs(output_dir, exist_ok=True)

    script_path = os.path.join(output_dir, f"submit_{job_name}.sh")
    xml_basename = os.path.basename(xml_path)

    # 2026-09-15 (审查 P2-1) 修复两处 Windows 致命问题:
    #   (1) 文本模式 open(...,'w') 在 Windows 下把 \n 写成 \r\n →
    #       shebang 变成 "#!/bin/bash\r", 内核找不到解释器, sbatch 必失败;
    #   (2) 直接把 Windows 路径 (反斜杠) 写进 #SBATCH --output 与 cd ——
    #       SLURM 在 Linux 上解析不了 "D:\a\b"。
    # 统一: 换行强制 LF; 路径统一转 POSIX 分隔符并加引号 (防空格)。
    def _posix(p: str) -> str:
        return str(p).replace('\\', '/')

    _out_posix = _posix(output_dir)

    module_lines = ""
    if modules:
        module_lines = "\n".join(f"module load {m}" for m in modules) + "\n"

    script = f"""#!/bin/bash
#SBATCH --job-name={job_name}
#SBATCH --output={_out_posix}/{job_name}_%j.out
#SBATCH --error={_out_posix}/{job_name}_%j.err
#SBATCH --time={time_limit}
#SBATCH --ntasks=1
#SBATCH --cpus-per-task={threads}
#SBATCH --mem={memory}
#SBATCH --partition={partition}

{module_lines}
echo "Starting BEAST job: {job_name}"
echo "XML: {xml_basename}"
echo "Date: $(date)"

cd "{_out_posix}"

{beast_bin} -threads {threads} "{xml_basename}"

EXIT_CODE=$?
echo "BEAST finished with exit code: $EXIT_CODE"
echo "Date: $(date)"
exit $EXIT_CODE
"""

    with open(script_path, 'w', encoding='utf-8', newline='\n') as f:
        f.write(script)
    try:
        os.chmod(script_path, 0o755)
    except OSError:
        pass  # Windows 上 chmod 语义有限, 非致命

    log.emit(f"SLURM script: {script_path}")
    log.emit(f"  Submit with: sbatch {script_path}")
    return script_path


def batch_submit_path_sampling(
    ps_dir: str,
    beast_bin: str = "beast",
    threads: int = 4,
    time_limit: str = "04:00:00",
    memory: str = "8G",
    partition: str = "normal",
    log: Optional[LogCollector] = None,
) -> List[str]:
    """
    Generate SLURM submission scripts for all path sampling steps.

    Returns list of script paths.
    """
    log = log or LogCollector()
    scripts = []

    for root, dirs, files in os.walk(ps_dir):
        if 'pathsampling.xml' in files:
            xml_path = os.path.join(root, 'pathsampling.xml')
            step_name = os.path.basename(os.path.dirname(root))
            run_name = os.path.basename(os.path.dirname(os.path.dirname(root)))
            job_name = f"ps_{run_name}_{step_name}"

            script = generate_slurm_submit_script(
                xml_path=xml_path,
                output_dir=root,
                job_name=job_name,
                beast_bin=beast_bin,
                threads=threads,
                time_limit=time_limit,
                memory=memory,
                partition=partition,
                log=None,
            )
            scripts.append(script)

    log.emit(f"Generated {len(scripts)} SLURM scripts")
    return scripts


# ═══════════════════════════════════════════════════════════════════
# Auto Tree Prior Selection (runs quick path sampling MCMC)
# ═══════════════════════════════════════════════════════════════════

def auto_select_tree_prior(
    fasta_file: str,
    metadata_csv: str,
    output_dir: str,
    clock_model: str = "ucln",
    quick_steps: int = 10,
    quick_chain: int = 200_000,
    beast_bin: str = "beast",
    threads: int = 4,
    beast_version: str = "1",
    mle_timeout: int = 3600,
    log: Optional[LogCollector] = None,
) -> Dict:
    """
    自动选择树先验：跑快速 path sampling，比较 constant vs skyline。

    流程:
      1. 生成 constant 和 skyline 两份 XML
      2. 各跑 quick_steps 步 × quick_chain 链长的 path sampling
      3. 收集 marginal likelihood
      4. 计算 log BF = log ML(skyline) - log ML(constant)
      5. 决策:
         log BF > 5  → 选 skyline (11/16论文标准)
         log BF < -5 → 选 constant (群体恒定)
         |log BF| ≤ 5 → 无法区分 → 选 skyline (文献标准, 保守)

    总耗时: ~20-40 分钟 (取决于 taxa 数量)

    Returns
    -------
    dict: {selected, skyline_ml, constant_ml, log_bf, reason, error,
           comparison_failed, run_errors}

    `comparison_failed=True` 表示**先验比较没有真正完成** (marginal likelihood
    算不出来), 此时 `selected` 只是回退默认值, `log_bf` 无效 —— 调用方必须
    把它当失败上报, 不能当比较结论使用。
    """
    log = log or LogCollector()
    os.makedirs(output_dir, exist_ok=True)
    result = {"selected": "skyline", "skyline_ml": None, "constant_ml": None,
              "log_bf": None, "reason": "", "error": None,
              "comparison_failed": False, "run_errors": []}

    try:
        # 2026-09-15 修复 (P1-9): 原实现无条件调用 BEAST1 的 XML 生成器,
        # 与下面 beast_version="1" 的硬编码一致; 现改为按 beast_version 分派,
        # 否则 BEAST2 配置下会生成 BEAST1 XML 再按 BEAST2 布局去找 step*。
        if str(beast_version) == "1":
            from utils.beast1_bridge import generate_beast1_phylogeo_xml as _gen_xml
        else:
            from utils.phylogeo_bridge import generate_phylogeo_xml as _gen_xml

        # BEAST 二进制可用性预检 (P1-9): 先确认 beast 能跑, 避免生成完 XML、
        # 跑到最后才发现二进制不在 —— 那种失败最难定位。
        if not os.path.isabs(beast_bin) and shutil.which(beast_bin) is None:
            _env_bin = os.environ.get("PHYLO_BEAST_BIN") or ""
            if _env_bin and (os.path.exists(_env_bin) or shutil.which(_env_bin)):
                log.warning(f"Auto tree prior: '{beast_bin}' 不在 PATH, "
                            f"改用 PHYLO_BEAST_BIN={_env_bin}")
                beast_bin = _env_bin
            else:
                result["comparison_failed"] = True
                result["error"] = (f"beast 二进制不可用 (PATH 中无 '{beast_bin}', "
                                   f"PHYLO_BEAST_BIN 未设置); 未做先验比较")
                result["reason"] = "beast 不可用, 回退默认 skyline (非比较结论)"
                log.error(result["error"])
                return result

        # ── 1. Generate XMLs ──
        log.emit("Auto tree prior: generating models...")

        const_dir = os.path.join(output_dir, "constant")
        const_xml_r = _gen_xml(
            fasta_file, metadata_csv, const_dir,
            clock_model=clock_model, tree_prior="constant",
            discretize_locations=False, log=log)
        if not const_xml_r["success"]:
            result["comparison_failed"] = True
            result["error"] = f"Constant XML failed: {const_xml_r.get('error')}"
            result["reason"] = "constant XML 生成失败, 回退默认 skyline (非比较结论)"
            log.error(result["error"])
            return result

        sky_dir = os.path.join(output_dir, "skyline")
        sky_xml_r = _gen_xml(
            fasta_file, metadata_csv, sky_dir,
            clock_model=clock_model, tree_prior="skyline",
            discretize_locations=False, log=log)
        if not sky_xml_r["success"]:
            result["comparison_failed"] = True
            result["error"] = f"Skyline XML failed: {sky_xml_r.get('error')}"
            result["reason"] = "skyline XML 生成失败, 回退默认 skyline (非比较结论)"
            log.error(result["error"])
            return result

        # ── 2. Generate + run path sampling ──
        const_ps = generate_path_sampling_xmls(
            const_xml_r["xml_path"], os.path.join(const_dir, "ps"),
            num_steps=quick_steps, chain_length=quick_chain,
            beast_version=beast_version, log=log)
        sky_ps = generate_path_sampling_xmls(
            sky_xml_r["xml_path"], os.path.join(sky_dir, "ps"),
            num_steps=quick_steps, chain_length=quick_chain,
            beast_version=beast_version, log=log)

        # 2026-09-15 (P1-9): XML 生成失败必须显式上报, 否则下面会在空目录上
        # 空转一轮再报一个含糊的 "ML 未解析到"。
        for _label, _ps in (("constant", const_ps), ("skyline", sky_ps)):
            if not _ps.get("success") or not _ps.get("steps"):
                result["comparison_failed"] = True
                result["error"] = (f"{_label} path sampling XML 生成失败: "
                                   f"{_ps.get('error') or '无 steps 产出'}")
                result["reason"] = "path sampling XML 生成失败, 回退默认 skyline (非比较结论)"
                log.error(result["error"])
                return result

        log.emit(f"\nRunning path sampling "
                 f"({'BEAST1 native MLE' if str(beast_version) == '1' else 'BEAST2 per-β'}, "
                 f"{quick_steps} β steps × {quick_chain:,} chain)...")

        # 2026-09-15 修复 (P1-9): 原实现无条件枚举 run1/step*/pathsampling.xml。
        # 但 BEAST1 原生 MLE 只产 run1/marginal_likelihood.xml, **不建 step* 目录**
        # (见 generate_path_sampling_xmls 的 BEAST1 分支 `continue`) —— 于是这个循环
        # 一次都不执行, 一步 MCMC 都没跑, 之后 collect 必然拿不到 ML, 最终静默
        # "defaulting to skyline"。用户以为做了先验比较, 实际什么都没算。
        # 现按版本分派: BEAST1 直接跑 run1/marginal_likelihood.xml (经 run_beast1_mle,
        # 带确定性种子与 stdout 重定向), BEAST2 才走 step* 逐 β。
        run_errors = []
        for label, ps_dir in [("constant", const_dir), ("skyline", sky_dir)]:
            ps_root = os.path.join(ps_dir, "ps")

            if str(beast_version) == "1":
                run_dir = os.path.join(ps_root, "run1")
                xml_file = os.path.join(run_dir, "marginal_likelihood.xml")
                if not os.path.exists(xml_file):
                    _e = f"{label}: 缺少 BEAST1 MLE XML ({xml_file})"
                    run_errors.append(_e)
                    log.error(f"  [{label}] {_e}")
                    continue
                log.emit(f"  [{label}] BEAST1 native MLE ...")
                _r = run_beast1_mle(xml_file, beast_bin=beast_bin, threads=threads,
                                    timeout=mle_timeout, log=log)
                if not _r.get("success"):
                    _e = f"{label}: {_r.get('error')}"
                    run_errors.append(_e)
                    log.error(f"  [{label}] BEAST1 MLE 失败: {_r.get('error')}")
                continue

            # BEAST2: 逐 β XML
            run_dir = os.path.join(ps_root, "run1")
            if not os.path.isdir(run_dir):
                _e = f"{label}: 缺少 BEAST2 path sampling 目录 ({run_dir})"
                run_errors.append(_e)
                log.error(f"  [{label}] {_e}")
                continue
            step_dirs = sorted(
                d for d in os.listdir(run_dir)
                if d.startswith('step') and os.path.isdir(os.path.join(run_dir, d)))
            if not step_dirs:
                _e = f"{label}: {run_dir} 下没有任何 step* 目录"
                run_errors.append(_e)
                log.error(f"  [{label}] {_e}")
                continue
            for step in step_dirs:
                step_dir = os.path.join(run_dir, step)
                xml_file = os.path.join(step_dir, "pathsampling.xml")
                if not os.path.exists(xml_file):
                    continue
                # 2026-09-15 修复 (P1-9 连带): 原写法 `log.emit(..., end=" ")` 会直接
                # 抛 TypeError —— LogCollector.emit(self, msg) 不接受 end 参数
                # (见 utils/virphy_bridge.py:56)。该异常会冒泡到函数最外层 except,
                # 让整个 auto prior 直接失败。改为先跑后记, 单条消息带结果。
                _outcome = "done"
                try:
                    _p = subprocess.run(
                        [beast_bin, "-threads", str(threads), xml_file],
                        capture_output=True, text=True, encoding='utf-8', errors='replace',
                        timeout=600, cwd=step_dir)
                    if _p.returncode != 0:
                        run_errors.append(f"{label}/{step}: exit {_p.returncode}")
                        _outcome = f"exit {_p.returncode}"
                except subprocess.TimeoutExpired:
                    run_errors.append(f"{label}/{step}: timeout")
                    _outcome = "timeout"
                except Exception as e:
                    run_errors.append(f"{label}/{step}: {e}")
                    _outcome = f"failed: {e}"
                log.emit(f"  [{label}] {step} ... {_outcome}")

        # ── 3. Collect marginal likelihoods ──
        log.emit("\nCollecting marginal likelihoods...")
        const_ml_r = collect_path_sampling_results(
            os.path.join(const_dir, "ps"), log=log)
        sky_ml_r = collect_path_sampling_results(
            os.path.join(sky_dir, "ps"), log=log)

        const_ml = const_ml_r.get("log_marginal_likelihood")
        sky_ml = sky_ml_r.get("log_marginal_likelihood")

        # 2026-09-15 修复 (P1-9): ML 缺失时不再伪装成"正常选择了 skyline"。
        # 原实现只写一条 warning 就 return, selected 保持 "skyline", 调用方
        # (run_phase_phylogeo) 会把它当正常比较结果写进 metrics.tree_prior_selected /
        # tree_prior_log_bf —— 用户看到"比较过了, 选了 skyline", 实际一步都没跑。
        # 现在显式标记 comparison_failed, 让上层能区分
        # "算出来选了 skyline" 与 "没算出来只能退回 skyline"。
        if const_ml is None or sky_ml is None:
            _missing = []
            if const_ml is None:
                _missing.append(f"constant ({const_ml_r.get('error') or '未解析到 ML'})")
            if sky_ml is None:
                _missing.append(f"skyline ({sky_ml_r.get('error') or '未解析到 ML'})")
            if run_errors:
                _missing.append("运行错误: " + "; ".join(run_errors[:4]))
            result["comparison_failed"] = True
            result["run_errors"] = run_errors
            result["error"] = "marginal likelihood 未算出, 先验比较未完成: " + " | ".join(_missing)
            result["reason"] = ("先验比较未完成 (marginal likelihood 计算失败), "
                                "回退默认 skyline —— 这不是比较结论, log BF 无效")
            result["selected"] = "skyline"
            log.error(f"Auto tree prior: {result['error']}")
            return result

        result["run_errors"] = run_errors

        log_bf = sky_ml - const_ml
        result["skyline_ml"] = sky_ml
        result["constant_ml"] = const_ml
        result["log_bf"] = round(log_bf, 2)

        # ── 4. Decision ──
        if log_bf > 5:
            result["selected"] = "skyline"
            result["reason"] = f"skyline significantly better (log BF={log_bf:.1f} > 5)"
        elif log_bf < -5:
            result["selected"] = "constant"
            result["reason"] = f"constant significantly better (log BF={log_bf:.1f} < -5)"
        else:
            result["selected"] = "skyline"
            result["reason"] = (f"models indistinguishable (|log BF|={abs(log_bf):.1f} ≤ 5); "
                               f"selecting skyline per literature standard (11/16 papers)")

        log.emit(f"\n  Decision: {result['selected']}")
        log.emit(f"  {result['reason']}")
        return result

    except Exception as e:
        # 2026-09-15 (P1-9): 异常路径同样是"比较未完成", 必须打标记,
        # 否则上层会把回退的 skyline 当成比较结论。
        result["comparison_failed"] = True
        result["error"] = str(e)
        log.error(f"Auto tree prior failed: {e}")
        result["reason"] = f"先验比较异常 ({e}), 回退默认 skyline (非比较结论)"
        return result


# ═══════════════════════════════════════════════════════════════════
# Tree Prior Comparison: Constant vs Skyline
# ═══════════════════════════════════════════════════════════════════

def compare_tree_priors(
    fasta_file: str,
    metadata_csv: str,
    output_dir: str,
    beast_version: str = "1",
    clock_model: str = "ucln",
    path_steps: int = 20,
    path_chain: int = 1_000_000,
    log: Optional[LogCollector] = None,
) -> Dict:
    """
    自动比较 constant vs skyline 树先验的边际似然。

    论文中少数做了此比较 (PVY China-Japan 2017, BETS paper)。
    多数直接使用 skyline (11/16 papers)。

    如果 log BF > 5 支持 skyline → 用 skyline (群体大小确实变化)
    如果 |log BF| ≤ 5 → 模型无差异 → 用 skyline (更灵活, 论文标准)
    如果 log BF < -5 → constant 更好 → 用 constant (群体大小恒定)

    Returns
    -------
    dict: {success, constant_dir, skyline_dir, constant_xml, skyline_xml}
    """
    log = log or LogCollector()
    os.makedirs(output_dir, exist_ok=True)
    result = {"success": False, "constant_dir": None, "skyline_dir": None,
              "constant_xml": None, "skyline_xml": None, "error": None}

    try:
        from utils.beast1_bridge import generate_beast1_phylogeo_xml

        # ── 1. Constant ──
        log.emit("=" * 55)
        log.emit("Model 1: Constant population size")
        const_dir = os.path.join(output_dir, "constant_prior")
        const_xml_r = generate_beast1_phylogeo_xml(
            fasta_file=fasta_file, metadata_csv=metadata_csv,
            output_dir=const_dir, clock_model=clock_model,
            tree_prior="constant", discretize_locations=False, log=log)
        if not const_xml_r["success"]:
            result["error"] = "Constant XML failed"
            return result

        const_ps_dir = os.path.join(const_dir, "path_sampling")
        generate_path_sampling_xmls(const_xml_r["xml_path"], const_ps_dir,
                                    num_steps=path_steps, chain_length=path_chain,
                                    beast_version=beast_version, log=log)
        result["constant_dir"] = const_ps_dir
        result["constant_xml"] = const_xml_r["xml_path"]

        # ── 2. Skyline ──
        log.emit("\n" + "=" * 55)
        log.emit("Model 2: Bayesian Skyline (5 groups)")
        sky_dir = os.path.join(output_dir, "skyline_prior")
        sky_xml_r = generate_beast1_phylogeo_xml(
            fasta_file=fasta_file, metadata_csv=metadata_csv,
            output_dir=sky_dir, clock_model=clock_model,
            tree_prior="skyline", discretize_locations=False, log=log)
        if not sky_xml_r["success"]:
            result["error"] = "Skyline XML failed"
            return result

        sky_ps_dir = os.path.join(sky_dir, "path_sampling")
        generate_path_sampling_xmls(sky_xml_r["xml_path"], sky_ps_dir,
                                    num_steps=path_steps, chain_length=path_chain,
                                    beast_version=beast_version, log=log)
        result["skyline_dir"] = sky_ps_dir
        result["skyline_xml"] = sky_xml_r["xml_path"]

        # ── 3. Instructions ──
        log.emit("\n" + "=" * 55)
        log.emit("Tree Prior Comparison Instructions")
        log.emit("=" * 55)
        log.emit(f"  1. Run path sampling for constant:")
        log.emit(f"     cd {const_ps_dir}/run1")
        log.emit(f"     for f in step*/pathsampling.xml; do beast -threads 4 $f; done")
        log.emit(f"  2. Run path sampling for skyline:")
        log.emit(f"     cd {sky_ps_dir}/run1")
        log.emit(f"     for f in step*/pathsampling.xml; do beast -threads 4 $f; done")
        log.emit(f"  3. Collect and compare:")
        log.emit(f"     python path_sampling.py --collect {const_ps_dir}")
        log.emit(f"     python path_sampling.py --collect {sky_ps_dir}")
        log.emit(f"")
        log.emit(f"  Decision rule:")
        log.emit(f"    log BF = log ML(skyline) - log ML(constant)")
        log.emit(f"    BF > 5  → skyline better → use skyline")
        log.emit(f"    BF < -5 → constant better → use constant")
        log.emit(f"    |BF| ≤ 5 → indistinguishable → use skyline (literature standard)")

        result["success"] = True
        return result

    except Exception as e:
        result["error"] = str(e)
        log.error(f"Tree prior comparison failed: {e}")
        return result


# ═══════════════════════════════════════════════════════════════════
# Clock Model Comparison: Strict vs UCLN
# ═══════════════════════════════════════════════════════════════════

def compare_clock_models(
    fasta_file: str,
    metadata_csv: str,
    output_dir: str,
    beast_version: str = "1",
    tree_prior: str = "constant",
    path_steps: int = 20,
    path_chain: int = 1_000_000,
    log: Optional[LogCollector] = None,
) -> Dict:
    """
    自动比较 strict vs UCLN 分子钟模型的边际似然。

    流程:
      1. 生成 strict clock BEAST XML + path sampling XMLs
      2. 生成 UCLN clock BEAST XML + path sampling XMLs
      3. 打印运行指令和比较方法

    论文参考: 16篇植物病毒文献中12篇做了此比较。
    BETS 论文建议: 优先用 UCLN，但若 strict 的 log ML 更高且 log BF > 5，用 strict。

    Returns
    -------
    dict: {success, strict_dir, ucln_dir, strict_xml, ucln_xml, comparison_guide}
    """
    log = log or LogCollector()
    os.makedirs(output_dir, exist_ok=True)
    result = {"success": False, "strict_dir": None, "ucln_dir": None,
              "strict_xml": None, "ucln_xml": None, "error": None}

    try:
        from utils.beast1_bridge import generate_beast1_phylogeo_xml

        # ── 1. Strict clock ──
        log.emit("=" * 55)
        log.emit("Model 1: Strict molecular clock")
        log.emit("=" * 55)
        strict_dir = os.path.join(output_dir, "strict_clock")
        strict_xml_r = generate_beast1_phylogeo_xml(
            fasta_file=fasta_file, metadata_csv=metadata_csv,
            output_dir=strict_dir, clock_model="strict",
            tree_prior=tree_prior, discretize_locations=False,
            log=log)
        if not strict_xml_r["success"]:
            result["error"] = "Strict clock XML failed"
            return result

        strict_ps_dir = os.path.join(strict_dir, "path_sampling")
        ps_strict = generate_path_sampling_xmls(
            source_xml=strict_xml_r["xml_path"], output_dir=strict_ps_dir,
            num_steps=path_steps, chain_length=path_chain,
            beast_version=beast_version, log=log)
        result["strict_dir"] = strict_ps_dir
        result["strict_xml"] = strict_xml_r["xml_path"]

        # ── 2. UCLN clock ──
        log.emit("\n" + "=" * 55)
        log.emit("Model 2: UCLN relaxed clock")
        log.emit("=" * 55)
        ucln_dir = os.path.join(output_dir, "ucln_clock")
        ucln_xml_r = generate_beast1_phylogeo_xml(
            fasta_file=fasta_file, metadata_csv=metadata_csv,
            output_dir=ucln_dir, clock_model="ucln",
            tree_prior=tree_prior, discretize_locations=False,
            log=log)
        if not ucln_xml_r["success"]:
            result["error"] = "UCLN clock XML failed"
            return result

        ucln_ps_dir = os.path.join(ucln_dir, "path_sampling")
        ps_ucln = generate_path_sampling_xmls(
            source_xml=ucln_xml_r["xml_path"], output_dir=ucln_ps_dir,
            num_steps=path_steps, chain_length=path_chain,
            beast_version=beast_version, log=log)
        result["ucln_dir"] = ucln_ps_dir
        result["ucln_xml"] = ucln_xml_r["xml_path"]

        # ── 3. Instructions ──
        log.emit("\n" + "=" * 55)
        log.emit("Model Comparison Instructions")
        log.emit("=" * 55)
        log.emit(f"  1. Run path sampling for strict clock:")
        log.emit(f"     cd {strict_ps_dir}/run1")
        log.emit(f"     for f in step*/pathsampling.xml; do beast -threads 4 $f; done")
        log.emit(f"  2. Run path sampling for UCLN clock:")
        log.emit(f"     cd {ucln_ps_dir}/run1")
        log.emit(f"     for f in step*/pathsampling.xml; do beast -threads 4 $f; done")
        log.emit(f"  3. Collect and compare:")
        log.emit(f"     python path_sampling.py --collect {strict_ps_dir}")
        log.emit(f"     python path_sampling.py --collect {ucln_ps_dir}")
        log.emit(f"")
        log.emit(f"  Decision rule (BETS paper):")
        log.emit(f"    log BF = log ML(UCLN) - log ML(strict)")
        log.emit(f"    BF > 5  → UCLN significantly better → use UCLN")
        log.emit(f"    BF < -5 → Strict significantly better → use strict")
        log.emit(f"    |BF| ≤ 5 → models indistinguishable → use UCLN (more flexible)")

        result["success"] = True
        return result

    except Exception as e:
        result["error"] = str(e)
        log.error(f"Clock comparison failed: {e}")
        return result
