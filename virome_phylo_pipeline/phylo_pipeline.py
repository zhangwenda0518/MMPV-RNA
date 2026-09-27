#!/usr/bin/env python3
"""
phylo_pipeline.py — MMPV-RNA 系统发育分析管线 (v2, 流程化 stage 编排)
======================================================================
主运行脚本（位于 pipeline 根目录），依赖模块在 utils/。

stage 体系（--stage，可逗号组合；自动补齐依赖；顺序: 数据→比对→诊断→树→群体/选择→时间→地理→报告）:
  prep        数据收集 (序列+采样日期+元数据)
  online      NCBI 参考补充 (SeqHarvester) → data/ncbi_ref
  align       MAFFT 全长比对 + 替换饱和分析 Iss → phylogeny/mafft.aln.fasta
  splitstree  分裂网络 NeighborNet (SplitsTree6)
  rdp5        RDP5 重组检测 → 默认删整条重组子+重新 MAFFT (rdp5_dropped.realn.fasta);
              --rdp5_mask 改用旧区间置N策略 (masked_N.fasta, 不重比对)
  tree        IQ-TREE 建树 → phylogeny/iqtree.treefile (仅作下游供树, 已去超长 bootstrap)
  popgen      群体遗传学 (pypopart: π/θ/Tajima/Fst/MJN) — 只吃比对
  host        宿主分化 (比对必需, 树可选增强)
  capheine    正选择分析 (CAWLign→IQ-TREE→HyPhy→DRHIP) — 吃 CDS 参考
  clock       分子钟联合分析 (2026-08-27 由 rtt+temporal+genes 三 stage 合并):
              整体+分基因 RTT/LTT + DRT 时间信号, 末尾自动汇总 time/clock_summary.tsv
              (DRT 默认 10 次随机化; 可独立 CLI: python -m utils.clock_analysis)
  beast       全长 BEAST 定年
  gene_dating 分基因 BEAST 定年 (吃 capheine 产物, full_tmrca 取自 beast)
  phylogeo    贝叶斯系统地理 CTMC+BSSVS + SpreaD3 (重计算, 多链后台)
  geo         地理统计 (Mantel/树地理/VirSpaceTime) — phylogeo 的快速预检
  report      汇总 CSV + HTML

已移除 (2026-08-27): saturation(并入 align), rtt/temporal/genes(并入 clock),
  concat(无下游消费者; 工具函数保留 multigene_tools.concatenate_alignments)
密码子偏好性 (ENC/RSCU/Fop): 独立工具目录 codon_usage/ (不在 stage 体系内)

默认参数: config.yaml (CLI 显式参数 > config.yaml > 代码默认)

用法:
  python phylo_pipeline.py --extract_dir 06_extraction/ --metadata metadata.tsv --stage all -t 20
  python phylo_pipeline.py --work_dir some_virus/ --alignment mafft.aln.fasta --tree tree.nwk \\
      --meta dates.csv --stage tree,clock,beast
  python phylo_pipeline.py --config config.yaml --stage all
"""

import argparse
import csv
import hashlib
import json
import logging
import os
import re
import subprocess
import sys
import time
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Optional

# ── Add project root ──
SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent
sys.path.insert(0, str(SCRIPT_DIR))
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

# 跨管线统一 I/O 布局 (mmpv_common/): ⑤ per-virus 模块目录名布局感知,
# legacy=data/phylogeny/... standard=01_data/02_phylogeny/...(按输出顺序编号)
from mmpv_common.io_layout import ph_dir, normalize_layout_env, dir_name  # noqa: E402

try:
    import yaml
except ImportError:
    yaml = None

from utils.data_collector import (
    load_sample_metadata,
    scan_extraction_dir,
    prepare_virus_inputs,
    find_ref_fasta,
)
from utils.virphy_bridge import (
    LogCollector,
    run_treetime_rtt,
    run_treedater_ltt,
    run_mantel_test,
    run_tree_geography,
    run_virspacetime,
    run_tree_host,
    run_host_differentiation,
    slice_genes_from_alignment,
    run_gene_rtt_batch,
    check_dependencies,
)
from utils.seqharvester_bridge import run_seqharvester as seq_harvester
from utils.seqharvester_bridge import get_taxid
from utils.seqgrouper_bridge import run_seqgrouper as seq_grouper
from utils.report_builder import build_html_report
# 2026-09-15 (P0-A/P0-B): 汇总 CSV 与 HTML 报告共用的头条指标口径
from utils.summary_metrics import headline_columns
from utils.beast_bridge import (
    prepare_beauti_inputs,
    run_beast,
    postprocess_beast,
)
from utils.phylogeo_bridge import (
    generate_phylogeo_xml,
    extract_migration_bf,
)
from utils.beast1_bridge import (
    generate_beast1_phylogeo_xml,
)
from run_phylogeo import submit_phylogeo_chains, submit_mcmc_chains
from utils.beast_parser import parse_beast_log
from utils.temporal_signal import (
    date_randomization_test,
)
# 2026-09-15 (审查 P3) 清理: 删除 9 个从未被引用的导入 ——
#   beast_bridge.generate_beast2_xml / phylogeo_bridge.run_treeannotator /
#   temporal_signal.check_temporal_signal_simple /
#   path_sampling.{generate_path_sampling_xmls, collect_path_sampling_results,
#                  run_bets_test, generate_slurm_submit_script, compare_tree_priors} /
#   spread3_viz.build_phylogeo_report
# (AST 逐名核对 used_count==0; 这些函数本身仍在各模块内保留, 需要时可显式导入)
from utils.path_sampling import (
    auto_select_tree_prior,
)
from utils.spread3_viz import (
    generate_spread3_json,
    generate_interactive_map_html,
)
from utils.dataset_config import dataset as get_dataset, known_dataset_names

# ═══════════════════════════════════════════════════════════════════
# Stage registry
# ═══════════════════════════════════════════════════════════════════

STAGE_ORDER = [
    'prep', 'metadata', 'online', 'clean', 'align', 'splitstree', 'rdp5', 'tree',
    'popgen', 'host', 'capheine',
    'clock', 'beast', 'gene_dating',
    'phylogeo', 'geo', 'geo_paths',
    'report',
]

STAGE_DEPS = {
    'prep':         [],
    'metadata':     ['prep'],
    'online':       ['prep', 'metadata'],
    'clean':        ['prep'],
    'align':        ['prep', 'clean'],
    'splitstree':   ['align'],
    'rdp5':         ['align'],
    'tree':         ['align'],
    'clock':        ['tree'],
    'popgen':       ['align'],
    'capheine':     ['prep'],
    'beast':        ['tree'],
    'gene_dating':  ['clock', 'capheine'],
    'phylogeo':     ['tree'],
    'geo':          ['tree'],
    'geo_paths':    ['phylogeo'],   # 2026-09-16 新增: 吃 phylogeo 的 merged/mcc.tree
    'host':         ['tree'],
    'report':       [],
}

# rtt/temporal/genes 已合并为 clock (子步骤名仅用于 checkpoint 兼容)
_LEGACY_CLOCK_SUBSTAGES = ('rtt', 'temporal', 'genes')

# stage → 模块目录 (产物归组, 对应 8 大模块)
# 2026-09-15 (审查 P3): 键集合与 STAGE_ORDER 对齐。原先含 'genes' / 'rtt' / 'temporal'
# 三个**不存在于 STAGE_ORDER 的 stage** —— 它们永远不会被 resolve_stages 选中,
# 留着只会让人以为"有这几个 stage"。
# 2026-09-16 修: 注释声称"与 STAGE_ORDER 对齐"但实际**漏了 'metadata'** ——
# metadata stage 的产物落在 `data/metadata/`, 模块就是 data。已补上, 现真对齐。
# 注: 本常量目前无代码引用 (纯信息), 保留供文档/外部脚本查阅;
# 由 `_consistency_check_20260916/verify_docs_vs_code.py` A4 断言守门。
STAGE_MODULE_DIR = {
    'prep': 'data', 'metadata': 'data', 'online': 'data', 'clean': 'data',
    'align': 'phylogeny', 'splitstree': 'phylogeny', 'tree': 'phylogeny',
    'rdp5': 'recomb',
    'capheine': 'select',
    'beast': 'time', 'gene_dating': 'time', 'clock': 'time',
    'phylogeo': 'geography', 'geo': 'geography', 'geo_paths': 'geography',
    'popgen': 'popgen', 'host': 'popgen',
    'report': 'report',
}

# 8 大模块分组 (--stage 顶层名) → 子步骤 (细粒度仍可单独跑)
STAGE_GROUPS = {
    'data':      ['prep', 'online', 'clean'],
    'popgen':    ['popgen', 'host'],
    'recomb':    ['rdp5'],
    'select':    ['capheine'],
    'phylogeny': ['align', 'splitstree', 'tree'],
    'time':      ['clock', 'beast', 'gene_dating'],
    'geography': ['phylogeo', 'geo', 'geo_paths'],
    'report':    ['report'],
}

STAGE_DESC = {
    'prep':         '数据收集 (序列+日期+元数据)',
    'metadata':     'metadata 治理 (时间+地理 检查与矫正) + 三通道拆分',
    'online':       'NCBI 参考补充 (SeqHarvester)',
    'clean':        '分析前序列清洗 (长度/N/大小写/非法字符/元数据缺失)',
    'align':        'MAFFT 全长比对 + 替换饱和分析 (Iss, 数据质量)',
    'splitstree':   '分裂网络 (SplitsTree6 NeighborNet)',
    'rdp5':         'RDP5 重组检测 + mask',
    'tree':         'IQ-TREE 建树',
    'clock':        '分子钟联合分析 (整体+分基因 RTT/LTT + DRT 时间信号)',
    'popgen':       '群体遗传学 (pypopart: π/θ/Tajima/Fst/MJN)',
    'capheine':     '正选择分析 (CAWLign→IQ-TREE→HyPhy→DRHIP)',
    'beast':        '全长 BEAST 定年',
    'gene_dating':  '分基因 BEAST 定年 (吃 capheine 产物)',
    'phylogeo':     '系统地理 (CTMC+BSSVS)',
    'geo':          '地理分析 (Mantel/树地理)',
    'geo_paths':    '事件级传播路径 (pathways/episode/置信分级 + 扩散统计 + 传播图)',
    'host':         '宿主分化',
    'report':       '汇总 CSV + HTML 报告',
}

# ── Checkpoint / 断点续传 ──
CHECKPOINT_DIR = ".checkpoints"

# 每 stage 的完成标志产物（相对 virus work_dir；.done 标记优先，产物存在兼容旧数据）
STAGE_CHECKPOINT_FILES = {
    # 2026-09-18 补: 原先漏了 metadata / clean（09-16 新增的两个 stage）。
    # 缺条目只是"没有标志产物回退" —— 只有 .done 标记一条判定路径, 老数据（无 .done）
    # 会被判未完成而每次重跑。产物路径见 RUN_GUIDE §2/§4。
    'metadata':    'data/metadata/metadata_std.csv',
    'clean':       'data/clean/clean.fasta',
    'online':      'data/ncbi_ref',
    'align':       'phylogeny/mafft.aln.fasta',
    'splitstree':  'phylogeny/splitstree/neighbornet.nexus',
    'rdp5':        'recomb/rdp5/rdp5_dropped.realn.fasta',  # 2026-08-27 改默认 drop+realign; --rdp5_mask 时为 masked/masked_N.fasta
    'tree':        'phylogeny/iqtree.treefile',
    'clock':       'time/treedater_ltt/Phylogeny_dated.pdf',
    'popgen':      'popgen/popgen_report.txt',
    'capheine':    'select/capheine/drhip',
    'beast':       'time/beast/beast1.xml',
    'gene_dating': 'time/gene_dating/gcva_watch_parse.sh',
    'phylogeo':    'geography/phylogeography/merged/mcc.tree',
    'geo':         'geography/geo_analysis/mantel_test.pdf',
    'geo_paths':   'geography/geo_analysis/pathways/pathway_qc.json',  # 2026-09-16 新增
    'host':        'popgen/host_analysis/tree_host.pdf',
    'report':      'report/phylo_summary.csv',
    'prep':        None,   # 数据收集不参与 checkpoint
}


def _cp_file(virus_dir, stage):
    return os.path.join(virus_dir, CHECKPOINT_DIR, f"{stage}.done")


def _inputs_fingerprint(paths, params: Optional[str] = None) -> str:
    """输入文件指纹 (路径 + 大小 + mtime 的 sha1 前 16 位)。

    2026-09-15 新增 (P1-1): checkpoint 原先只判断"标志文件是否存在且非空",
    **没有任何输入指纹校验** —— 换了 MAFFT 参数 / 加了新序列 / 换了 RDP5 策略
    (drop 的序列集不同 → 重比对得到新坐标) 之后, 只要 iqtree.treefile 还在,
    tree 就会拿旧比对建出的旧树静默复用。现把输入指纹记进 .done 文件。
    缺失文件记 "MISSING" (稳定值), 不会因文件不存在而反复触发重跑。

    2026-09-16 新增 params: **工具参数**也参与指纹。动机: 给 IQ-TREE 补自举后,
    比对文件没变但建树参数变了, 旧 checkpoint 会让 tree 整段跳过 → 修复不生效。
    """
    h = hashlib.sha1()
    for p in paths:
        if not p:
            continue
        try:
            st = os.stat(p)
            h.update(f"{os.path.abspath(p)}|{st.st_size}|{int(st.st_mtime)}".encode())
        except OSError:
            h.update(f"{p}|MISSING".encode())
    if params:
        h.update(f"|params:{params}".encode())
    return h.hexdigest()[:16]


def stage_done(virus_dir, stage, force=False, inputs=None) -> bool:
    """stage 是否已完成。

    2026-09-15 (P1-1) 新增 inputs 指纹校验:
      · force=True -> 恒 False (强制重跑)
      · .done 存在且**记录了** inputs 指纹, 与当前 inputs 不符 -> 视为未完成 (重跑)
      · .done 存在但没有 inputs 记录 (升级前的旧 checkpoint) -> 沿用旧行为 (True),
        并在本次 mark 时补记指纹 —— 保证首次升级不会触发全量重跑
      · 否则回落到"完成标志产物存在"
    """
    if force:
        return False
    cp = _cp_file(virus_dir, stage)
    if os.path.exists(cp):
        if inputs is not None:
            try:
                with open(cp, encoding='utf-8') as f:
                    data = json.load(f)
            except Exception:
                data = {}
            prev = data.get('inputs')
            if prev is not None and prev != inputs:
                return False          # 输入变了 -> 旧产物已过期
        return True
    marker = STAGE_CHECKPOINT_FILES.get(stage)
    if marker:
        p = os.path.join(virus_dir, marker)
        if os.path.exists(p) and os.path.getsize(p) > 0:
            return True
    return False


def mark_stage_done(virus_dir, stage, meta=None) -> None:
    cp = _cp_file(virus_dir, stage)
    os.makedirs(os.path.dirname(cp), exist_ok=True)
    data = {"stage": stage, "time": datetime.now().isoformat(timespec="seconds"), "status": "ok"}
    if meta:
        data.update(meta)
    with open(cp, "w", encoding="utf-8") as f:
        json.dump(data, f, ensure_ascii=False, indent=1)


def clear_stage_checkpoint(virus_dir, stage=None) -> None:
    """清除 checkpoint (stage=None 清全部)"""
    cp_dir = os.path.join(virus_dir, CHECKPOINT_DIR)
    if stage:
        p = _cp_file(virus_dir, stage)
        if os.path.exists(p):
            os.remove(p)
        return
    if os.path.isdir(cp_dir):
        import shutil
        shutil.rmtree(cp_dir, ignore_errors=True)

# config.yaml section → argparse dest 映射
CONFIG_SECTION_ARG = {
    'runtime': {},  # key 与 dest 同名
    'rdp5': {'enabled': 'rdp5', 'script': 'rdp5_script', 'genes': 'rdp5_genes',
             'mask_min_methods': 'mask_min_methods', 'stage': 'rdp5_stage'},
    'beast': {'enabled': 'beast', 'bin': 'beast_bin', 'chain': 'beast_chain',
              'burnin': 'beast_burnin', 'timeout_hours': 'beast_timeout',
              'skip_run': 'skip_beast_run', 'version': 'beast_version',
              # 2026-09-16 修 (与 phylogeo.chains 同一类): `chains` 原先**不在映射里**
              # → config.yaml 里写的 `beast.chains: 5` 被静默忽略, 实际一直用
              # argparse 默认 3。文档收敛时才发现 (见 _docconv_20260916/)。
              'chains': 'beast_chains'},
    'gene_dating': {'enabled': 'gene_dating', 'genes_dir': 'genes_dir',
                    'chain': 'gene_chain', 'chain_big': 'gene_chain_big',
                    'big_genes': 'big_genes'},
    # 2026-09-22 新增: clean 段此前**整段未登记** ——
    # datasets.yaml 的 fasta 是**已比对文件**(含 '-' gap), 而 clean 默认把 gap 当
    # 非法字符整条剔除 → 全部 REMOVE → clean.fasta 为空 → 下游无输入 (PSTVd 实测
    # 111/111 全被剔除, 比对文件 0 字节)。登记本段后可在 config.yaml 里持久化
    # `clean.allow_gap_chars: true`, 不必每次命令行加 --clean_allow_gap_chars。
    'clean': {'allow_gap_chars': 'clean_allow_gap_chars',
              'strip_illegal': 'clean_strip_illegal',
              'min_length_ratio': 'clean_min_length_ratio',
              'max_length_ratio': 'clean_max_length_ratio',
              'drop_no_date': 'clean_drop_no_date',
              'drop_no_location': 'clean_drop_no_location'},
    'popgen': {'group': 'popgen_group', 'min_n': 'popgen_min_n', 'perm': 'popgen_perm',
               'exclude_reference': 'popgen_exclude_reference',
               'exclude_ids': 'popgen_exclude_ids'},
    'capheine': {'ref': 'capheine_ref', 'unaligned': 'capheine_unaligned',
                 'code': 'capheine_code', 'workers': 'capheine_workers',
                 'cpus_iqtree': 'capheine_cpus_iqtree', 'cpus_hyphy': 'capheine_cpus_hyphy',
                 'use_mpi': 'capheine_mpi'},
    'phylogeo': {'enabled': 'phylogeo', 'clock': 'phylogeo_clock',
                 'prior': 'phylogeo_prior', 'chain': 'phylogeo_chain',
                 'bf': 'phylogeo_bf', 'treeannotator_bin': 'treeannotator_bin',
                 # 2026-09-16 修: `chains` 原先**不在映射里** → config.yaml 里写的
                 # `phylogeo.chains: 5` 被静默忽略, 实际一直用 argparse 默认 3。
                 # 这类"键没登记就悄悄丢"的坑, 加键时必须同步登记。
                 'chains': 'phylogeo_chains',
                 # RRT 树注释法(VirPhyKit 口径)的随机化 MCC 树目录; 不填=不跑(默认关闭)
                 'rrt_randomized_dir': 'rrt_randomized_dir'},
    # 2026-09-16 新增: geo_paths 子阶段 (事件级传播路径)。键必须与 g_geo 的
    # --no_pathways/--direct_snp/--indirect_snp/--censor_years/--n_perm/--geo_map/--geo_gif
    # 同步 (chains 键静默丢失的教训)。
    'geo_paths': {'no_pathways': 'no_pathways', 'direct_snp': 'direct_snp',
                  'indirect_snp': 'indirect_snp', 'censor_years': 'censor_years',
                  'n_perm': 'n_perm', 'map': 'geo_map', 'gif': 'geo_gif',
                  'gif_frames': 'geo_gif_frames', 'tree_map': 'geo_tree_map',
                  'calibrate': 'calibrate'},
    'temporal': {"check": "check_temporal", "drt_randomizations": "drt_randomizations",
                 'bets': 'bets', 'bets_steps': 'bets_steps', 'bets_chain': 'bets_chain'},
    'saturation': {'replicates': 'sat_replicates'},
    'report': {'enabled': 'report_enabled'},
}


def validate_stage_graph() -> None:
    """启动自检 (2026-09-18): 依赖必须排在**使用者之前**，否则拒绝启动。

    ⚠️ 为什么需要这道闸门
    --------------------
    `process_virus` 里**没有 `for stage in stages` 循环** —— 每个 stage 是一段
    `if '<stage>' in stages:` 硬编码块，所以**真实执行顺序 = 这些块在文件里的物理顺序**，
    而 `resolve_stages` 的依赖闭包只保证"被选中"，**不保证先后**。

    历史事故: `geo_paths` 声明依赖 `phylogeo`（要吃它的 `merged/mcc.tree`），却因为块被写在
    `phylogeo` 之前 → 首次全新跑时 mcc.tree 还不存在 → 整段静默跳过；而 `_skip` 不计入
    `stages_failed` → **跑完仍然报成功**，要跑到第二/三遍才出结果。
    2026-09-18 已把物理块序对齐到 STAGE_ORDER（STAGE_ORDER 现即真实执行序）；本函数再钉一道:
    有人日后改 STAGE_ORDER 时, 只要把使用者排到其依赖前面, 立刻在此报错而不是静默跳过。
    """
    pos = {s: i for i, s in enumerate(STAGE_ORDER)}
    bad = []
    for s in STAGE_ORDER:
        for d in STAGE_DEPS.get(s, ()):
            if d in pos and pos[d] > pos[s]:
                bad.append(f"{s}(第{pos[s]+1}位) 依赖 {d}(第{pos[d]+1}位)")
    if bad:
        raise SystemExit(
            "[stage] STAGE_ORDER 违反依赖序（依赖必须排在使用者之前）: " + "; ".join(bad)
            + "\n        （执行顺序 = process_virus 里 if 块的物理顺序，改 STAGE_ORDER 时"
              "必须同步移动对应代码块）")


def resolve_stages(stage_arg: str) -> List[str]:
    """解析 --stage: 'all' 展开全部; 逗号列表自动补齐依赖; 返回按 STAGE_ORDER 的执行清单"""
    validate_stage_graph()
    if stage_arg in (None, '', 'all'):
        return list(STAGE_ORDER)
    requested = [s.strip() for s in str(stage_arg).split(',') if s.strip()]
    # 展开大模块名 → 子步骤
    expanded = []
    for s in requested:
        if s in STAGE_GROUPS:
            expanded.extend(STAGE_GROUPS[s])
        else:
            expanded.append(s)
    for s in expanded:
        if s not in STAGE_DEPS:
            raise SystemExit(f"[stage] 未知 stage: '{s}'. 大模块: {', '.join(STAGE_GROUPS)} | 子步骤: {', '.join(STAGE_ORDER)}")
    full = set()
    for s in expanded:
        stack = [s]
        while stack:
            cur = stack.pop()
            if cur in full:
                continue
            full.add(cur)
            for d in STAGE_DEPS[cur]:
                stack.append(d)
    return [s for s in STAGE_ORDER if s in full]


def load_config(path: Optional[str] = None) -> Dict:
    """读取 config.yaml → dict (缺文件/无 yaml 时返回 {})"""
    if yaml is None:
        return {}
    p = path or os.path.join(str(SCRIPT_DIR), 'config.yaml')
    if not os.path.exists(p):
        return {}
    try:
        with open(p, encoding='utf-8') as f:
            return yaml.safe_load(f) or {}
    except Exception as e:
        print(f"[config] 读取失败: {e}")
        return {}


def config_defaults(cfg: Dict) -> Dict:
    """config dict → {argparse dest: value}"""
    flat = {}
    for section, mapping in CONFIG_SECTION_ARG.items():
        sec = cfg.get(section, {})
        for key, value in sec.items():
            dest = mapping.get(key, key if section == 'runtime' else None)
            if dest is None:
                continue
            if key == 'stage' and section == 'runtime':
                dest = 'stage'
            flat[dest] = value
    return flat


def config_orphan_keys(cfg: Dict) -> List[str]:
    """列出 config.yaml 里**登记不到 dest、会被静默丢弃**的键。

    2026-09-16 新增 (文档收敛时发现 `beast.chains` 就是这个坑):
    `config_defaults()` 对非 runtime 分区取 `mapping.get(key)` —— 没登记的键直接
    `continue`，**不报错、不警告**，于是 `config.yaml` 里写了也等于没写。
    历史上 `phylogeo.chains` 与 `beast.chains` 都中过招 (用户以为跑了 5 条链,
    实际是 argparse 默认 3)。本函数把这类"静默失效"变成**可检测**的，
    由 `verify_docs_vs_code.py` 与启动期自检共同消费。

    Returns
    -------
    List[str]: `section.key` 形式的孤儿键列表 (空 = 全部已登记)
    """
    orphans: List[str] = []
    for section, sec in (cfg or {}).items():
        if not isinstance(sec, dict):
            continue
        if section not in CONFIG_SECTION_ARG:
            orphans.append(f"{section}.* (整段未登记)")
            continue
        mapping = CONFIG_SECTION_ARG[section]
        for key in sec:
            if section == 'runtime' or key in mapping:
                continue
            orphans.append(f"{section}.{key}")
    return orphans


def _explicit_cli_dests(parser) -> set:
    """返回命令行上**真正显式出现**过的参数 dest 集合。

    2026-09-15 新增 (修 P1-5): 用于让 `--config` 不覆盖 CLI 显式传参。
    做法是扫描 sys.argv 里出现的选项字符串并映射回 dest —— 比 deepcopy parser
    再把 default 换成 SUPPRESS 更稳 (argparse 对象并不保证可深拷贝)。
    支持 `--opt val` / `--opt=val` / 合并短选项 `-tv` 三种写法。
    """
    opt2dest = {}
    for a in getattr(parser, "_actions", []):
        for o in getattr(a, "option_strings", []) or []:
            opt2dest[o] = a.dest
    explicit = set()
    for tok in sys.argv[1:]:
        if tok == "--":
            break
        if not tok.startswith("-") or tok == "-":
            continue
        key = tok.split("=", 1)[0]
        if key in opt2dest:
            explicit.add(opt2dest[key])
            continue
        # 合并短选项: -tv == -t -v
        if len(key) > 2 and key[1] != "-":
            for ch in key[1:]:
                if "-" + ch in opt2dest:
                    explicit.add(opt2dest["-" + ch])
    return explicit


# ═══════════════════════════════════════════════════════════════════
# Logger / tools
# ═══════════════════════════════════════════════════════════════════

def setup_logger(output_dir: str, level: str = "INFO") -> logging.Logger:
    logger = logging.getLogger("PhyloPipeline")
    logger.setLevel(logging.DEBUG)
    logger.handlers.clear()

    os.makedirs(output_dir, exist_ok=True)

    ch = logging.StreamHandler()
    ch.setLevel(getattr(logging, level.upper(), logging.INFO))
    ch.setFormatter(logging.Formatter('[%(asctime)s] %(message)s', datefmt='%H:%M:%S'))
    logger.addHandler(ch)

    fh = logging.FileHandler(os.path.join(output_dir, 'phylo_pipeline.log'), encoding='utf-8')
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(logging.Formatter('[%(asctime)s] %(levelname)s %(message)s'))
    logger.addHandler(fh)

    return logger


def find_tool(name: str) -> Optional[str]:
    result = subprocess.run(['where', name] if os.name == 'nt' else ['which', name],
                           capture_output=True, text=True, encoding='utf-8', errors='replace')
    if result.returncode == 0:
        return result.stdout.strip().split('\n')[0]
    return None


def check_tools(logger: logging.Logger) -> Dict[str, Optional[str]]:
    """启动时探测外部工具。

    2026-09-18 修 (自检盲区): 原先只查 mafft/iqtree/iqtree2/Rscript，而代码实际还会调用
    beast (24 处) / hyphy (13) / lsd2 / raxml-ng —— 启动看到「工具 ✓」并不代表这些阶段
    能跑，要到那一步才失败或静默跳过。现把「默认链要用的」拆成两组分别报告:
      · 必需: 本地链 (prep→clean→align→tree) 缺任一项就跑不动
      · 可选: 对应阶段开关打开时才需要 (beast / capheine / clock 的 LSD2 对照 / 交叉验证)
    注: rdp5 不是可执行文件 (走 `python -m recombination_analysis`)，treetime 是 pip 包
        (已在 check_dependencies 里查)，两者都不该出现在这里。
    """
    REQUIRED = ['mafft', 'iqtree', 'iqtree2', 'Rscript']
    OPTIONAL = ['beast', 'hyphy', 'lsd2', 'raxml-ng']
    tools: Dict[str, Optional[str]] = {}
    for name in REQUIRED + OPTIONAL:
        tools[name] = find_tool(name)
    for name, found in tools.items():
        kind = '必需' if name in REQUIRED else '可选'
        status = '✓' if found else '✗ NOT FOUND'
        logger.info(f"  {name} [{kind}]: {status}")
    miss_opt = [n for n in OPTIONAL if not tools.get(n)]
    if miss_opt:
        logger.info(f"  (可选工具缺失: {', '.join(miss_opt)} —— 只要不开对应阶段就无影响)")
    return tools


def run_mafft(input_fasta, output_aln, threads=8, logger=None) -> bool:
    log = logger or logging.getLogger(__name__)
    cmd = f'mafft --auto --thread {threads} --reorder "{input_fasta}" > "{output_aln}"'
    log.info(f"  MAFFT: {Path(input_fasta).name} → {Path(output_aln).name}")
    try:
        subprocess.run(cmd, shell=True, check=True, timeout=1800)
        return True
    except subprocess.CalledProcessError as e:
        log.error(f"  MAFFT failed (exit={e.returncode})")
        return False
    except subprocess.TimeoutExpired:
        log.error(f"  MAFFT timed out")
        return False


# ═══════════════════════════════════════════════════════════════════
# IQ-TREE 支持值 (2026-09-16 补: 主树原不带自举 → 论文里没有分支支持值)
# ═══════════════════════════════════════════════════════════════════
# 历史: v1 (2026-08-27) 曾写 `-m -m utils.MFP -bb 1000`（`-m` 重复 + `utils.MFP`
# 并非合法模型名 → 实际跑不起来），后为控耗时改成 `-m MFP --quiet` 并把 bootstrap
# 整个去掉（文件头注释原话"已去超长 bootstrap"）→ 主树没有任何支持值。
#
# 房屋口径（与平台侧一致，勿各自发明）: 平台 `Virus_Platform_Core/phylo.py::_run_iqtree`
# 用的是 `-m MFP -B 1000 -alrt 1000`（UFBoot + SH-aLRT 双支持值）。
# 本处按同一口径补齐；参考项目 zika_Vietnam 的建树步骤同样带支持值
# （RAxML-NG `--bs-trees 1000 --bs-metric fbp`；IQ-TREE `-b 100`）。
IQTREE_BOOT_DEFAULT = 1000      # UFBoot 重复数 (IQ-TREE 硬性要求 >=1000)
IQTREE_ALRT_DEFAULT = 1000      # SH-aLRT 重复数 (同样要求 >=1000)
IQTREE_TIMEOUT_DEFAULT = 3600   # **单次尝试**超时(秒); 带 UFBoot 后需按数据规模上调
IQTREE_SUPPORT_META = "iqtree.support.json"   # 支持值元数据 sidecar (供下游/报告取证)


def _iqtree_support_flags(boot, alrt, bnni=False) -> List[str]:
    """把支持值参数翻成 IQ-TREE 旗标 (自动过滤不合法值)。

    ⚠️ `-B` / `-alrt` 在 IQ-TREE 里都**要求 >=1000** —— 想"降级跑少一点重复"
    是不行的 (会被直接拒绝)。所以降级路径只能换成"另一种支持值"或"不要支持值"，
    见 `run_iqtree` 的三段阶梯。
    """
    flags: List[str] = []
    try:
        boot = int(boot or 0)
    except (TypeError, ValueError):
        boot = 0
    try:
        alrt = int(alrt or 0)
    except (TypeError, ValueError):
        alrt = 0
    if boot >= 1000:
        flags += ["-B", str(boot)]
        if bnni:
            # UFBoot2 的 NNI 优化: 对"近缘序列多"的数据集能压低 UFBoot 的高估
            flags.append("--bnni")
    if alrt >= 1000:
        flags += ["-alrt", str(alrt)]
    return flags


def compose_iqtree_cmd(iqtree_bin, alignment, prefix, threads=8, model="MFP",
                       flags=None, boot=IQTREE_BOOT_DEFAULT,
                       alrt=IQTREE_ALRT_DEFAULT, bnni=False) -> List[str]:
    """构造 IQ-TREE 命令行 (**列表形式**, 不走 shell —— 中文/空格路径不会被拆)。

    `flags` 显式给出时直接用它(供降级阶梯复用); 否则由 boot/alrt/bnni 推出。
    """
    if flags is None:
        flags = _iqtree_support_flags(boot, alrt, bnni)
    return [str(iqtree_bin), "-s", str(alignment), "-pre", str(prefix),
            "-nt", str(int(threads)), "-m", str(model)] + list(flags) + ["--quiet"]


# 兼容旧名 (文档/测试里引用过)
build_iqtree_cmd = compose_iqtree_cmd


def support_type_from_flags(flags: List[str]) -> str:
    """由**实际用过的旗标**判定支持值类型 (权威来源, 优于从树标签反猜)。

    `-alrt 1000` + `-B 1000` → `SH-aLRT/UFBoot`（IQ-TREE 写标签为 `aLRT/UFBoot`）；
    只有 UFBoot → `UFBoot`；只有 `-b` → `bootstrap`；都没有 → `none`。
    """
    has_uf = "-B" in flags or "-bb" in flags or "--ufboot" in flags
    has_alrt = "-alrt" in flags or "--alrt" in flags
    has_boot = "-b" in flags or "--boot" in flags
    if has_alrt and has_uf:
        return "SH-aLRT/UFBoot"
    if has_uf:
        return "UFBoot"
    if has_alrt:
        return "SH-aLRT"
    if has_boot:
        return "bootstrap"
    return "none"


def tree_support_stats(treefile: str) -> Dict:
    """统计树的分支支持值 (内部节点数 / 带标签数 / 标签形态)。

    只做**只读**解析, 用于: ① 判断已有树要不要因"缺支持值"而重建;
    ② 把 `n_labeled` 落进 metrics 供报告引用。
    IQ-TREE 的标签形如 `0/56`( aLRT/UFBoot )、`100/`(仅 aLRT )、`98`(单一支持值)。

    ⚠️ **必须同时看 `clade.name` 和 `clade.confidence`** (2026-09-16 踩):
    Bio.Phylo 的 Newick 解析器把**纯数字**标签放进 `confidence`、`name` 留 None,
    而带斜杠的 `95/100` 才进 `name`。只看 `name` 会把"只有单一支持值"的树
    (如 `-b 100` 的传统自举、或只跑 UFBoot) 误判成"无支持值" → `run_stage_tree`
    每次都判定要重建, 白跑一遍 IQ-TREE。
    """
    out = {"n_internal": 0, "n_labeled": 0, "has_support": False, "examples": []}
    try:
        from Bio import Phylo
        with open(treefile, encoding="utf-8", errors="replace") as f:
            tree = next(iter(Phylo.parse(f, "newick")), None)
        if tree is None:
            return out
        for clade in tree.get_nonterminals():
            if clade is tree.root:
                continue                      # 根节点没有可检验的分支支持值
            out["n_internal"] += 1
            nm = (clade.name or "").strip()
            lab = ""
            if nm and re.fullmatch(r"[\d.]+(?:/[\d.]*)?|[\d.]*/[\d.]+", nm):
                lab = nm
            elif clade.confidence is not None:
                try:
                    lab = "%g" % float(clade.confidence)
                except (TypeError, ValueError):
                    lab = str(clade.confidence)
            if lab:
                out["n_labeled"] += 1
                if len(out["examples"]) < 3:
                    out["examples"].append(lab)
        out["has_support"] = out["n_labeled"] > 0
    except Exception:
        pass
    return out


def write_support_meta(output_dir: str, payload: Dict) -> str:
    p = os.path.join(output_dir, IQTREE_SUPPORT_META)
    try:
        with open(p, "w", encoding="utf-8") as f:
            json.dump(payload, f, ensure_ascii=False, indent=2)
    except OSError:
        pass
    return p


def read_support_meta(output_dir: str) -> Optional[Dict]:
    p = os.path.join(output_dir, IQTREE_SUPPORT_META)
    try:
        with open(p, encoding="utf-8") as f:
            return json.load(f)
    except Exception:
        return None


def run_iqtree(alignment, output_dir, threads=8, iqtree_bin="iqtree2", logger=None,
               boot=IQTREE_BOOT_DEFAULT, alrt=IQTREE_ALRT_DEFAULT, bnni=False,
               model="MFP", timeout=IQTREE_TIMEOUT_DEFAULT) -> Optional[str]:
    """IQ-TREE 建树, **默认带 UFBoot + SH-aLRT 双支持值**。

    2026-09-16 修复: 原实现是 `-m MFP --quiet`（无任何自举）→ 主树没有分支支持值,
    论文里无法报告节点置信度。现按平台侧口径补 `-B 1000 -alrt 1000`。

    三段阶梯（**任何降级都必须留痕，不允许静默丢支持值**）：
      ① `-B {boot} -alrt {alrt} [--bnni]` —— 期望路径
      ② ①失败 → 退 `-alrt {alrt}`（SH-aLRT 便宜且稳; UFBoot 在小比对/近缘序列上
         本就会失败 —— v1 的注释"ultrafast bootstrap may fail for small alignments"
         说的就是这件事），**仍带支持值**
      ③ ②也失败 → 退纯 `-m MFP`，`log.warning` 明确说"支持值缺失"，并且
         sidecar 里记 `degraded: true` + 失败原因, 供报告端如实披露

    无论走哪条, 都把**实际命令行 / 支持值类型 / 降级状态**写进
    `<output_dir>/iqtree.support.json`, 并在日志里报 `n_labeled/n_internal` ——
    这样"有没有支持值"可被机械核对, 不靠人回忆。

    Parameters
    ----------
    boot, alrt : int
        UFBoot / SH-aLRT 重复数; 传 0 表示**不要**该支持值 (两者都 0 = 回到修复前行为)。
    bnni : bool
        给 UFBoot 加 `--bnni` (对近缘序列压低 UFBoot 高估), 以少数时间换更稳的支持值。
    timeout : int
        **单次尝试**的超时秒数。带 UFBoot 后耗时明显上升, 大比对需上调
        (命令行 `--iqtree_timeout`)。

    Returns
    -------
    str | None: treefile 路径; 失败返回 None
    """
    log = logger or logging.getLogger(__name__)
    os.makedirs(output_dir, exist_ok=True)
    prefix = os.path.join(output_dir, "iqtree")

    attempts = [
        ("full", _iqtree_support_flags(boot, alrt, bnni)),
    ]
    if alrt and int(alrt or 0) >= 1000:
        attempts.append(("alrt_only", _iqtree_support_flags(0, alrt, False)))
    attempts.append(("no_support", []))

    tried, last_err = [], None
    for kind, flags in attempts:
        cmd = compose_iqtree_cmd(iqtree_bin, alignment, prefix, threads, model,
                                 flags=flags)
        if kind != "full":
            log.warning(f"  IQ-TREE 重试({kind}): {' '.join(cmd)}")
        else:
            log.info(f"  IQ-TREE: {Path(alignment).name} → {prefix}.treefile")
        try:
            r = subprocess.run(cmd, timeout=timeout, cwd=output_dir,
                               capture_output=True, text=True,
                               encoding="utf-8", errors="replace")
            rc = r.returncode
        except subprocess.TimeoutExpired:
            rc, last_err = None, f"timeout({timeout}s)"
        except OSError as e:
            rc, last_err = None, str(e)
        tried.append({"attempt": kind, "flags": flags, "rc": rc, "err": last_err})

        treefile = prefix + ".treefile"
        if rc == 0 and os.path.exists(treefile) and os.path.getsize(treefile) > 0:
            st = tree_support_stats(treefile)
            sup_type = support_type_from_flags(flags)
            degraded = (kind == "no_support") and bool(
                (_iqtree_support_flags(boot, alrt, bnni)))
            # `degraded` 只表示"最终没有支持值"; "退到单支持值"另用 fallback 标出 ——
            # 否则 aLRT-only 会被误当成"拿到了请求的双支持值"。
            fallback = kind != "full"
            payload = {
                "tree": treefile,
                "tool": str(iqtree_bin),
                "requested": {"boot": int(boot or 0), "alrt": int(alrt or 0), "bnni": bool(bnni)},
                "actual": {"flags": flags, "support_type": sup_type,
                           "boot": int(boot or 0) if "-B" in flags else 0,
                           "alrt": int(alrt or 0) if ("-alrt" in flags) else 0},
                "support_type": sup_type,
                "degraded": bool(degraded),
                "success_attempt": kind,
                "fallback": bool(fallback),
                "n_internal": st["n_internal"],
                "n_labeled": st["n_labeled"],
                "attempts": tried,
                "cmd": " ".join(cmd),
                "time": datetime.now().isoformat(timespec="seconds"),
            }
            if degraded:
                payload["degrade_reason"] = (
                    "支持值计算两次尝试均失败 (UFBoot 易在小比对/近缘序列上失败); "
                    "已退到无支持值建树 —— 报告里不得声称有分支支持值")
                log.warning("  ⚠ IQ-TREE 未能算出支持值, 已退到 ``-m MFP``; "
                            "详见 iqtree.support.json")
            elif fallback:
                first = tried[0] if tried else {}
                rsn = ("①(-B/-alrt 双支持值)失败: %s; 已退到 %s"
                       % (first.get("err") or ("exit=%s" % first.get("rc")), kind))
                payload["fallback_reason"] = rsn
                log.warning("  ⚠ IQ-TREE 请求的双支持值未能算出 (%s); "
                            "本次支持值为 %s (已记入 iqtree.support.json)" % (rsn, sup_type))
            write_support_meta(output_dir, payload)
            if sup_type == "none":
                log.warning(f"  IQ-TREE 建树完成但**无支持值**: {treefile}")
            else:
                log.info(f"  IQ-TREE 完成, 支持值 {sup_type} "
                         f"({st['n_labeled']}/{st['n_internal']} 内部节点有标签)")
            return treefile

        if rc is None:
            log.warning(f"  IQ-TREE 尝试失败({kind}): {last_err}")
        else:
            tail = (r.stderr or r.stdout or "").strip().splitlines()[-1:] 
            log.warning(f"  IQ-TREE 尝试失败({kind}, exit={rc}): "
                        f"{tail[0] if tail else ''}")

    log.error("  IQ-TREE 三次尝试均失败 (含无支持值兜底), 无 treefile")
    write_support_meta(output_dir, {
        "tree": None, "tool": str(iqtree_bin),
        "requested": {"boot": int(boot or 0), "alrt": int(alrt or 0), "bnni": bool(bnni)},
        "support_type": "none", "degraded": True,
        "degrade_reason": "IQ-TREE 全部尝试失败, 没有产物",
        "attempts": tried, "time": datetime.now().isoformat(timespec="seconds"),
    })
    return None



# ═══════════════════════════════════════════════════════════════════
# Stage runners
# ═══════════════════════════════════════════════════════════════════

def run_stage_metadata(prep, args, logger) -> Dict:
    """metadata (入口治理): 时间 + 地理 的检查与矫正 —— 产出标准化 metadata + 体检报告

    2026-09-16 新增 (老师要求「时间和地理的检查和矫正，补充到我们的流程中对 metadata
    信息进行格式化，以免下游错误」; 老师拍板: 另存新文件 / Unknown 层级删除 /
    virphy_bridge 用 mode='start')。

    定位: 全流程最上游的 metadata 把关。**只做归一 + 标注, 不改任何判定阈值**。
    产出 (输出到 data/metadata/):
      · metadata_std.csv   —— 标准化后 metadata (下游唯一输入)
      · metadata_report.tsv —— 逐行体检报告 (日期精度/问题, 地理层级/删除, 列内小数年)
      · metadata_summary.txt —— 汇总统计 (人类可读)

    原件**不动**。下游若需要, 由编排层把 metadata 路径改写为 metadata_std.csv。
    """
    from utils import metadata_governance as _mg

    vname = prep["virus_name"]
    out_dir = prep["output_dir"]
    log = logger

    # 2026-09-16: 来源回退 —— 此前只认 `--metadata` (batch 入口),
    # 于是 work_dir 模式 (`--meta`) 与 --virus 模式 (datasets.yaml) 的 metadata
    # 虽然给了, 治理却整段跳过 (实测日志 "[metadata] 无 metadata 输入 (), 跳过")。
    #
    # 三种入口统一优先级:
    #   --metadata (batch) > prep['metadata'] > prep['meta_csv'] / ['dates_csv']
    #   > --meta / --host_meta (兜底)
    #
    # ⚠ 不能简单用 `or` 串联: work_dir 模式下 `--meta meta.csv` 是**未经解析的
    # 相对字符串**, 而 prep['meta_csv'] 才是 `_resolve()` 过的绝对路径。若把
    # args.meta 排在前面, `os.path.exists("meta.csv")` 会以 **CWD** 为基准判定
    # 失败 → 表现成"明明给了 metadata 却整段跳过"(2026-09-16 实测:
    # `[metadata] _e2e_ch: 无 metadata 输入 (meta.csv), 跳过`)。
    # 故改为「按优先级取第一个**真实存在**的路径」, 相对/绝对一律以存在性裁决。
    _cand = (getattr(args, "metadata", None), prep.get("metadata", ""),
             prep.get("meta_csv", ""), prep.get("dates_csv", ""),
             getattr(args, "meta", None), getattr(args, "host_meta", None))
    src = next((c for c in _cand if c and os.path.exists(c)), "")
    if not src:
        _shown = next((c for c in _cand if c), "")
        log.warning(f"  [metadata] {vname}: 无 metadata 输入 ({_shown}), 跳过")
        return {"success": False, "error": "no metadata input"}

    # 幂等: 若输入已是本 stage 的产物 (例如编排层回写后重复调用), 直接复用,
    # 避免"自己治理自己的产物"造成二次处理 / 永久改名。
    if os.path.basename(src) == "metadata_std.csv" and os.path.basename(
            os.path.dirname(src)) == "metadata":
        log.info(f"  [metadata] {vname}: 输入已是治理产物, 跳过重复治理")
        return {"success": True, "std_csv": src, "report_tsv": "",
                "summary": "", "metrics": {}, "result": {}, "outputs": {}}

    md_dir = os.path.join(out_dir, ph_dir('data'), "metadata")
    os.makedirs(md_dir, exist_ok=True)
    std_csv = os.path.join(md_dir, "metadata_std.csv")
    rep_tsv = os.path.join(md_dir, "metadata_report.tsv")
    summary = os.path.join(md_dir, "metadata_summary.txt")

    mode = getattr(args, "metadata_mode", "mid")
    try:
        res = _mg.govern_table(
            src,
            date_col=getattr(args, "metadata_date_col", None) or None,
            loc_col=getattr(args, "metadata_loc_col", None) or None,
            name_col=getattr(args, "metadata_name_col", None) or None,
            mode=mode,
            drop_unknown_geo=not getattr(args, "metadata_keep_unknown_geo", False),
            derive_coords=not getattr(args, "metadata_no_derive_coords", False),
            blank_placeholders=not getattr(args, "metadata_no_blank_placeholders", False),
            drop_placeholder_rows=getattr(args, "metadata_drop_placeholder_rows", False),
            drop_no_date=getattr(args, "metadata_drop_no_date", False),
            drop_no_location=getattr(args, "metadata_drop_no_location", False),
        )
        paths = _mg.write_outputs(res, std_csv, rep_tsv)
    except Exception as e:  # 治理失败不应让整条流程崩 —— 降级用原 metadata
        log.warning(f"  [metadata] {vname}: 治理失败 ({e}), 降级使用原 metadata")
        return {"success": False, "error": str(e)}

    s = res["stats"]
    lines = [
        f"metadata 治理报告 — {vname}",
        f"源文件: {src}",
        f"行数: {s['n_total']}  保留: {s['kept']}",
        f"识别列: date={s['date_col']!r}  loc={s['loc_col']!r}  name={s['name_col']!r}",
        f"坐标列: {s['coord_cols']}",
        f"日期口径 mode: {mode}  ({'月中/年中' if mode == 'mid' else '月初/年初 (VirPhyKit)'})",
        "",
        "日期精度分布:",
    ]
    for k, v in s["date"].items():
        lines.append(f"  {k:<8} {v}")
    if s["date_issues"]:
        lines.append(f"日期问题: {s['date_issues']}")
    if s["geo_issues"]:
        lines.append(f"地理问题: {s['geo_issues']}")
    lines.append(f"已删除 Unknown 层级: {s['geo_drop_unknown']} 行 / {s['geo_dropped_levels']} 级")
    if s["coords_derived"]:
        lines.append(f"坐标补算成功: {s['coords_derived']}  来源={s['coords_derived_by_source']}")
    if s["coords_unrecoverable"]:
        lines.append(f"坐标补算不出 (已清空): {s['coords_unrecoverable']}")
    if s["coord_placeholder_fields"]:
        lines.append(f"坐标列占位符 (补算后仍为空): {s['coord_placeholder_fields']}")
    if s["blanked_fields"]:
        lines.append(f"字段清空统计 (占位符 → 空): {s['blanked_fields']}")
        lines.append(f"  含占位符字段的行: {s['n_rows_with_placeholder']}")
    lines.append(f"列内已是小数年(未重复转换): {s['already_decimal_in_col']} 条")
    lines.append(f"剔除: 占位符行 {s['dropped_placeholder']};  "
                 f"无日期 {s['dropped_no_date']};  无地理 {s['dropped_no_location']}")
    with open(summary, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")

    n_drop = s["dropped_no_date"] + s["dropped_no_location"] + s["dropped_placeholder"]
    log.info(f"  [metadata] {vname}: {s['n_total']} 行 → 保留 {s['kept']}"
             + (f" (剔除 {n_drop})" if n_drop else "")
             + f"; 日期 day/month/year = {s['date'].get('day', 0)}/{s['date'].get('month', 0)}/{s['date'].get('year', 0)}"
             + f"; 删 Unknown 地理层级 {s['geo_drop_unknown']} 行")
    if s["blanked_fields"]:
        log.warning(f"  [metadata] {vname}: 占位符字段已清空 {s['blanked_fields']} "
                    f"(共 {s['n_rows_with_placeholder']} 行受影响; 详见 {rep_tsv})")
    if s["date_issues"]:
        log.warning(f"  [metadata] {vname}: 日期问题 {s['date_issues']} (详见 {rep_tsv})")
    if s["geo_issues"]:
        log.warning(f"  [metadata] {vname}: 地理问题 {s['geo_issues']} (详见 {rep_tsv})")
    if s["kept"] == 0:
        log.error(f"  [metadata] {vname}: 治理后无可用行, 降级使用原 metadata")
        return {"success": False, "error": "no rows after governance",
                "outputs": {"std_csv": std_csv, "report_tsv": rep_tsv, "summary": summary},
                "metrics": s}

    return {"success": True, "std_csv": std_csv, "report_tsv": rep_tsv,
            "summary": summary, "metrics": s, "result": res,
            "outputs": {"std_csv": std_csv, "report_tsv": rep_tsv, "summary": summary}}


def run_stage_online(prep, args, logger) -> Dict:
    """online: NCBI 参考补充 → ncbi_ref/ (efetch 单参考) + ncbi_harvest/ (SeqHarvester 按物种名)"""
    vname = prep["virus_name"]
    out_dir = prep["output_dir"]
    log = logger
    log.info(f"\n[Stage online] {vname} — NCBI 参考补充")
    results = {"virus": vname, "success": False, "ncbi_sequences": 0, "ncbi_fasta": None,
               "metrics": {}, "outputs": {}}

    # ── ① efetch 单病毒参考 (优先 datasets.yaml 的 accession, 回退到目录名解析) ──
    acc_match = prep.get("accession", "")
    if not acc_match:
        parts = vname.split('_')
        for p in parts:
            if re.search(r'^[A-Z]{1,2}_\d+\.\d+$', p) or re.search(r'^[A-Z]{2}\d+\.\d+$', p):
                acc_match = p
                break
    if acc_match:
        ncbi_dir = os.path.join(out_dir, ph_dir('data'), "ncbi_ref")
        os.makedirs(ncbi_dir, exist_ok=True)
        try:
            from Bio import Entrez, SeqIO
            Entrez.email = "mmvp@example.com"
            Entrez.tool = "MMPV-PhyloPipeline"
            log.info(f"  Fetching NCBI reference: {acc_match}")
            handle = Entrez.efetch(db="nucleotide", id=acc_match, rettype="gb", retmode="text")
            gb_file = os.path.join(ncbi_dir, f"{acc_match}.gbk")
            with open(gb_file, 'w') as f:
                f.write(handle.read())
            handle.close()
            for rec in SeqIO.parse(gb_file, "genbank"):
                fasta_file = os.path.join(ncbi_dir, f"{acc_match}.ref.fasta")
                SeqIO.write(rec, fasta_file, "fasta")
                results["ncbi_fasta"] = fasta_file
                log.info(f"  Reference: {rec.description[:80]}...")
                break
            results["ncbi_sequences"] = 1
            results["success"] = True
        except Exception as e:
            log.warning(f"  NCBI fetch failed: {e}")
    else:
        log.info("  无法从目录名解析 NCBI accession, 跳过 efetch 单参考")

    # ── ② SeqHarvester 按 taxid/物种名补充参考序列 (taxid 优先, 全长过滤) ──
    species_name = prep.get("species", "")
    if not species_name:
        species_name = re.sub(r'_?(?:NC_|OR|MW|MK|MT|MN|LC|AB|DQ|AY|AF)\d+[\.\d]*$', '', vname)
        species_name = species_name.replace('_', ' ').strip()
    if species_name or acc_match:
        harv_dir = os.path.join(out_dir, ph_dir('data'), "ncbi_harvest")
        taxid = get_taxid(accession=acc_match or prep.get("accession", ""), species=species_name)
        if taxid:
            log.info(f"  [online] SeqHarvester: taxid={taxid} (全长过滤)")
        else:
            log.info(f"  [online] SeqHarvester: '{species_name}' (全长过滤)")
        try:
            ncbi = seq_harvester(species_name, harv_dir, max_total=0,
                                 taxid=taxid, full_length=True, log=LogCollector(log))
            if ncbi and ncbi.get("success"):
                n_seq = ncbi.get("n_sequences", results.get("ncbi_sequences", 0))
                results["ncbi_sequences"] = n_seq
                results["metrics"]["ncbi_sequences"] = n_seq
                if ncbi.get("fasta_file"):
                    results["outputs"]["ncbi_harvest_fasta"] = ncbi["fasta_file"]
                if ncbi.get("metadata_csv") and os.path.exists(ncbi["metadata_csv"]):
                    results["outputs"]["ncbi_metadata"] = ncbi["metadata_csv"]
                    sg = seq_grouper(ncbi["metadata_csv"], harv_dir, "location", log=LogCollector(log))
                    if sg and sg.get("success"):
                        results["outputs"]["ncbi_group_loc"] = sg["pie_chart"]
                    sg2 = seq_grouper(ncbi["metadata_csv"], harv_dir, "host", log=LogCollector(log))
                    if sg2 and sg2.get("success"):
                        results["outputs"]["ncbi_group_host"] = sg2["pie_chart"]
                results["success"] = True
                # 2026-09-18 明确化 (原先是"分析死胡同"却没人说明):
                # 这是**收集/报告**通路, 不是分析输入 —— 产物只在 report 层被引用。
                log.info("  [online] 注: ncbi_harvest/ 的产物仅供报告（分组图 + 序列清单），"
                         "**不进** combined.fasta / 比对 / 建树 / 定年。")
                log.info("  [online]     要让它真扩样本: 把 data/ncbi_harvest/ncbi_sequences.fasta "
                         "并入 06_extraction/<病毒目录>/ 后重跑 prep（并确认对应 metadata 也在表里）。")
        except Exception as e:
            log.warning(f"  SeqHarvester failed: {e}")
    else:
        log.info("  无法从目录名解析物种名, 跳过 SeqHarvester")
    return results


def run_stage_clean(prep, args, logger) -> Dict:
    """clean (第一道, 比对前): 序列本体清洗 —— 长度/N/大小写/非法字符/元数据缺失

    2026-09-16 新增 (老师要求「分析前，再次长度检查，去除长度异常，N 碱基异常等，
    修正大小写，去除时间和地理都没记录的序列等」)。

    与后续 align_qc (第二道, 比对后) 的分工:
      · 本道: 只看序列本体 —— 必须在 MAFFT 之前, 因为空格/星号/小写
        会让 MAFFT 本身报错或产出错位比对, 事后 QC 只能"亡羊补牢"
      · align_qc: 看需要比对坐标系的指标 (gap% / identity / 比对坐标系 N%)

    产出 data/clean/clean.fasta (+ removed.fasta + clean_report.tsv)。
    调用方需把 prep['combined_fasta'] 指向 clean.fasta。
    """
    import importlib
    from utils import seq_clean as _sc

    vname = prep["virus_name"]
    out_dir = prep["output_dir"]
    log = logger
    src = prep.get("combined_fasta", "")
    if not src or not os.path.exists(src):
        log.warning(f"  [clean] {vname}: 无输入序列 ({src}), 跳过")
        return {"success": False, "error": "no input fasta"}

    clean_dir = os.path.join(out_dir, ph_dir('data'), "clean")
    os.makedirs(clean_dir, exist_ok=True)
    clean_fa = os.path.join(clean_dir, "clean.fasta")
    report = os.path.join(clean_dir, "clean_report.tsv")

    # 2026-09-16: metadata 缺省回退。
    # 此前默认 None → 用户开了 --clean_drop_no_date/-location 却忘了 --clean_metadata,
    # 会走到 seq_clean 的 SystemExit, 被本函数捕获 → **整个 clean stage 降级**,
    # 等于"去除时间和地理都没记录的序列"这条要求静默失效 (老师原话的那条)。
    # 回退顺序刻意只取**本病毒的两张通道表**: 它们与 fasta 的样品集一致;
    # 若退到全局 metadata, 凡"表里没这一行"的序列都会被判无日期/无地理而删掉。
    meta_path = (getattr(args, "clean_metadata", None)
                 or prep.get("meta_csv")
                 or None)
    if meta_path and not os.path.exists(meta_path):
        log.warning(f"  [clean] metadata 不存在: {meta_path} (时间/地理门槛不生效)")
        meta_path = None
    try:
        r = _sc.clean_sequences(
            src,
            reference=getattr(args, "clean_reference", None) or None,
            min_length_ratio=getattr(args, "clean_min_length_ratio", 0.9),
            max_length_ratio=getattr(args, "clean_max_length_ratio", 1.5),
            max_n=getattr(args, "clean_max_n", 0.05),
            allow_rna=getattr(args, "clean_allow_rna", False),
            strip_illegal=getattr(args, "clean_strip_illegal", False),
            metadata_path=meta_path,
            drop_no_date=getattr(args, "clean_drop_no_date", False),
            drop_no_location=getattr(args, "clean_drop_no_location", False),
            strict_gap_chars=not getattr(args, "clean_allow_gap_chars", False),
            log=LogCollector(log),
        )
        produced = _sc.write_outputs(r, clean_fa, report)
    except SystemExit as e:
        log.warning(f"  [clean] {vname}: {e}")
        return {"success": False, "error": str(e)}

    st = r["stats"]
    log.info(f"  [clean] {vname}: 保留 {st['n_keep']} / 剔除 {st['n_removed']} "
             f"(共 {st['n_input']}); 基准 {st['basis']}bp ({st['basis_src']})")

    # 全部被剔光 → 明确失败, 不要让下游拿到空文件
    if st["n_keep"] == 0:
        log.error(f"  [clean] {vname}: 所有序列均被剔除, 无法继续")
        return {"success": False, "error": "all sequences removed by clean",
                "outputs": produced, "metrics": st}

    return {"success": True, "clean_fasta": clean_fa,
            "outputs": produced, "metrics": st, "result": r}


def run_stage_align(prep, tools, args, logger, force=False) -> Dict:
    """align: MAFFT 全长比对 → mafft.aln.fasta"""
    vname = prep["virus_name"]
    out_dir = prep["output_dir"]
    log = logger
    aln_file = os.path.join(out_dir, ph_dir('phylogeny'), "mafft.aln.fasta")

    # 产物已存在: 直接复用 (work-dir 提供的比对或上次成功产物)
    # 2026-09-15 (P1-1): 加 force 判断 —— 旧版即使 --force 也复用已存在产物,
    # 用户以为强制重跑, 实际什么都没变。
    if not force and os.path.exists(aln_file) and os.path.getsize(aln_file) > 0:
        log.info(f"  [align] {vname}: MAFFT alignment exists, use")
        return {"success": True, "alignment": aln_file}
    if force and os.path.exists(aln_file):
        log.info(f"  [align] {vname}: --force, 重新计算比对 (忽略已有 {Path(aln_file).name})")

    provided = prep.get("combined_fasta", "")
    if provided and os.path.exists(provided) and provided != aln_file:
        log.info(f"  [align] {vname}: using provided alignment: {Path(provided).name}")
        return {"success": True, "alignment": provided}

    if args.skip_msa:
        log.warning(f"  [align] {vname}: --skip_msa but no alignment available, skip")
        return {"success": False, "error": "no alignment available (skip_msa)"}

    if not provided or not os.path.exists(provided):
        log.warning(f"  [align] {vname}: input sequences missing: {provided}")
        return {"success": False, "error": "input sequences missing"}

    # work-dir 模式: combined_fasta 即比对文件本身, 空/损坏时无原始输入可重算, 拒绝 MAFFT
    if os.path.abspath(provided) == os.path.abspath(aln_file):
        return {"success": False,
                "error": "alignment file empty/corrupt (work-dir mode, no raw input to recompute)"}

    os.makedirs(os.path.join(out_dir, ph_dir('phylogeny')), exist_ok=True)
    log.info(f"  [align] {vname}: running MAFFT...")
    ok = run_mafft(provided, aln_file, args.threads, log)
    if not ok:
        return {"success": False, "error": "MAFFT failed"}
    return {"success": True, "alignment": aln_file}


def run_stage_rdp5(aln_file, prep, tools, args, logger) -> Optional[str]:
    """rdp5: RDP5 重组检测 → 删除重组子 (默认) 或 mask 区间 (--rdp5_mask)

    默认策略 (2026-08-27 改): drop_recombinants 删整条重组序列
      → 对剩余序列重新 MAFFT → 后续 tree/clock/beast 等全部吃重比对
    保留策略: --rdp5_mask 时走旧 mask 区间置 N (不重比对, 保样本量)
    返回: 后续 stage 应使用的比对路径"""
    vname = prep["virus_name"]
    out_dir = prep["output_dir"]
    log = logger
    rdp5_dir = os.path.join(out_dir, ph_dir('recomb'), "rdp5")
    os.makedirs(rdp5_dir, exist_ok=True)

    # 短序列防御 (2026-09-02): 类病毒 (如 PSTVd ~359nt) 太短, RDP5 窗口方法无足够位点,
    # wine 下静默跑不出 CSV → 记 skip 而非 fail
    try:
        with open(aln_file, errors="replace") as _af:
            _first_seq = []
            for _line in _af:
                if _line.startswith(">"):
                    if _first_seq:
                        break
                else:
                    _first_seq.append(_line.strip())
        _aln_len = max((len(s) for s in _first_seq), default=0)
    except OSError:
        _aln_len = 0
    if 0 < _aln_len < 500:
        log.info(f"  [rdp5] 比对仅 {_aln_len} bp (<500, 类病毒/超短序列), RDP5 窗口方法不适用, 记 skip")
        return aln_file

    log.info(f"  [rdp5] {vname}: RDP5 recombination detection (stage={args.rdp5_stage})...")

    rdp5_script = os.path.expanduser(args.rdp5_script)
    rec_cmd = [sys.executable, "-m", "recombination_analysis",
               "--fasta", aln_file, "--prefix", vname, "--outdir", rdp5_dir,
               "--rdp5-script", rdp5_script, "--stage", args.rdp5_stage,
               "--threads", str(args.threads)]
    if args.rdp5_genes:
        rec_cmd += ["--genes", args.rdp5_genes]
    rr = subprocess.run(rec_cmd, capture_output=True, text=True, encoding='utf-8', errors='replace', cwd=str(SCRIPT_DIR))
    if rr.returncode != 0:
        log.warning(f"  [rdp5] RDP5 failed (exit={rr.returncode}): {(rr.stderr or rr.stdout)[-400:]}")
        return None

    rdp5_csv = os.path.join(rdp5_dir, f"{vname}.csv")
    if not os.path.exists(rdp5_csv):
        log.warning(f"  [rdp5] RDP5 CSV not found: {rdp5_csv}")
        return None
    log.info(f"  [rdp5] RDP5 events CSV: {rdp5_csv}")

    # ── 策略 A (默认): 删重组子 → 重新 MAFFT ──
    if not getattr(args, "rdp5_mask", False):
        drop_dir = os.path.join(rdp5_dir, "dropped")
        drop_cmd = [sys.executable, "-m", "mask_recombination",
                    "--fasta", aln_file, "--rdp5", rdp5_csv,
                    "--min-methods", str(args.mask_min_methods),
                    "--outdir", drop_dir]
        if args.rdp5_genes:
            drop_cmd += ["--genes", args.rdp5_genes]
        dr = subprocess.run(drop_cmd, capture_output=True, text=True, encoding='utf-8', errors='replace', cwd=str(SCRIPT_DIR))
        if dr.returncode != 0:
            log.warning(f"  [rdp5] drop_recombinants failed: {(dr.stderr or dr.stdout)[-300:]}")
            return None
        dropped_fa = os.path.join(drop_dir, "dropped_recombinants.fasta")
        if not os.path.exists(dropped_fa):
            log.warning(f"  [rdp5] drop output not found: {dropped_fa}")
            return None
        # 删除后为纯序列 (含旧 gap 列), 重新 MAFFT 建立新比对坐标
        realn = os.path.join(rdp5_dir, "rdp5_dropped.realn.fasta")
        log.info(f"  [rdp5] realigning {Path(dropped_fa).name} (MAFFT) → {Path(realn).name}...")
        if not run_mafft(dropped_fa, realn, threads=args.threads, logger=log):
            log.warning("  [rdp5] realign failed, fallback to masked-strategy alignment")
            return None
        log.info(f"  [rdp5] recombinant-dropped + realigned: {realn}")
        return realn

    # ── 策略 B (--rdp5_mask): 旧 mask 区间置 N (不重比对) ──
    mask_cmd = [sys.executable, "-m", "mask_recombination",
                "--fasta", aln_file, "--rdp5", rdp5_csv,
                "--min-methods", str(args.mask_min_methods),
                "--outdir", os.path.join(rdp5_dir, "masked"),
                "--mask"]
    if args.rdp5_genes:
        mask_cmd += ["--genes", args.rdp5_genes]
    mr = subprocess.run(mask_cmd, capture_output=True, text=True, encoding='utf-8', errors='replace', cwd=str(SCRIPT_DIR))
    if mr.returncode != 0:
        log.warning(f"  [rdp5] mask_recombination failed: {(mr.stderr or mr.stdout)[-300:]}")
        return None
    masked_fa = os.path.join(rdp5_dir, "masked", "masked_N.fasta")
    if os.path.exists(masked_fa):
        log.info(f"  [rdp5] masked alignment: {masked_fa}")
        return masked_fa
    log.warning(f"  [rdp5] mask output not found: {masked_fa}")
    return None


def _iqtree_params_from_args(args) -> Dict:
    """从命令行取出 IQ-TREE 支持值/超时参数 (带默认值, 兼容老 args 对象)。"""
    def _i(name, default):
        try:
            return int(getattr(args, name, default) or 0)
        except (TypeError, ValueError):
            return default
    return {
        "boot": _i("iqtree_boot", IQTREE_BOOT_DEFAULT),
        "alrt": _i("iqtree_alrt", IQTREE_ALRT_DEFAULT),
        "bnni": bool(getattr(args, "iqtree_bnni", False)),
        "timeout": _i("iqtree_timeout", IQTREE_TIMEOUT_DEFAULT) or IQTREE_TIMEOUT_DEFAULT,
    }


def iqtree_param_token(args) -> str:
    """IQ-TREE **参数指纹** —— 拼进 tree stage 的 checkpoint inputs。

    为什么需要: `stage_done()` 原先只比"输入文件指纹"。补上支持值后, 比对文件没变
    但**建树参数变了**, 旧 checkpoint 会让 tree stage 整个跳过 → 旧的无支持值树被
    静默沿用, 修复等于没生效。把它并进 inputs 指纹, 参数一变就自动重跑一次
    (之后稳定, 不会反复触发)。
    """
    p = _iqtree_params_from_args(args)
    return "iqtree|bin=%s|model=MFP|boot=%s|alrt=%s|bnni=%s" % (
        getattr(args, "iqtree_bin", "iqtree2"), p["boot"], p["alrt"], p["bnni"])


def run_stage_tree(prep, aln_file, tools, args, logger, force=False) -> Dict:
    """tree: IQ-TREE 建树 → iqtree.treefile (rdp5 开启时 aln_file 已为 masked)

    2026-09-16: 主树默认带 UFBoot + SH-aLRT 支持值 (原为纯 `-m MFP`, 无支持值)。
    复用旧树前会**核对支持值**是否齐全 —— 否则"修复"在已有产物上不会生效。
    """
    vname = prep["virus_name"]
    out_dir = prep["output_dir"]
    log = logger
    phylo_dir = os.path.join(out_dir, ph_dir('phylogeny'))
    tree_file = os.path.join(phylo_dir, "iqtree.treefile")
    iq = _iqtree_params_from_args(args)
    want_support = iq["boot"] >= 1000 or iq["alrt"] >= 1000

    # 产物已存在: 复用 (但要求支持值满足本次请求)
    # 2026-09-15 (P1-1): 加 force 判断 (旧版 --force 也复用旧树)
    if not force and os.path.exists(tree_file) and os.path.getsize(tree_file) > 0:
        st = tree_support_stats(tree_file)
        if (not want_support) or st["has_support"]:
            log.info(f"  [tree] {vname}: IQ-TREE tree exists, use"
                     f" ({st['n_labeled']}/{st['n_internal']} 节点有支持值)" if st["has_support"]
                     else f"  [tree] {vname}: IQ-TREE tree exists, use (未请求支持值)")
            return {"success": True, "tree": tree_file,
                    "support": {"type": support_type_from_flags(
                        _iqtree_support_flags(iq['boot'], iq['alrt'], iq['bnni']))
                        if want_support else "none",
                        "n_labeled": st["n_labeled"],
                        "n_internal": st["n_internal"]}}
        log.warning(f"  [tree] {vname}: 已有主树**无分支支持值** "
                    f"({st['n_labeled']}/{st['n_internal']}), 但本次要求 "
                    f"UFBoot/alrt 支持值 → 重建")
    if force and os.path.exists(tree_file):
        log.info(f"  [tree] {vname}: --force, 重新建树 (忽略已有 {Path(tree_file).name})")

    provided = prep.get("tree_file", "")
    if provided and os.path.exists(provided) and provided != tree_file:
        log.info(f"  [tree] {vname}: using provided tree: {Path(provided).name}")
        st = tree_support_stats(provided)
        if want_support and not st["has_support"]:
            # 外部提供的树不可重建 → 必须说清"这篇论文里这棵树没有支持值"
            log.warning(f"  [tree] {vname}: 外部提供的树无分支支持值, 且无法重建; "
                        f"支持值将缺失 (如需支持值请去掉 --tree 让本管线自建)")
        return {"success": True, "tree": provided,
                "support": {"type": "provided", "n_labeled": st["n_labeled"],
                            "n_internal": st["n_internal"]}}

    if args.skip_tree:
        log.warning(f"  [tree] {vname}: --skip_tree but no tree available, skip")
        return {"success": False, "error": "no tree available (skip_tree)"}

    if not aln_file or not os.path.exists(aln_file):
        log.warning(f"  [tree] {vname}: no alignment for IQ-TREE, skip")
        return {"success": False, "error": "no alignment for IQ-TREE"}

    os.makedirs(phylo_dir, exist_ok=True)
    log.info(f"  [tree] {vname}: running IQ-TREE on {Path(aln_file).name} "
             f"(MFP + boot={iq['boot']} alrt={iq['alrt']}"
             f"{' +bnni' if iq['bnni'] else ''}, timeout={iq['timeout']}s)...")
    t = run_iqtree(aln_file, phylo_dir, args.threads, args.iqtree_bin, log,
                   boot=iq["boot"], alrt=iq["alrt"], bnni=iq["bnni"],
                   timeout=iq["timeout"])
    if not t:
        return {"success": False, "error": "IQ-TREE failed"}
    meta = read_support_meta(phylo_dir) or {}
    if meta.get("degraded"):
        log.warning(f"  [tree] {vname}: ⚠ 本次主树**没有**支持值 —— "
                    f"报告/正文不得声称有分支支持值 ({meta.get('degrade_reason', '')})")
    elif meta.get("fallback"):
        log.warning(f"  [tree] {vname}: ⚠ 本次主树只拿到 {meta.get('support_type')} "
                    f"(请求的双支持值降级; 详见 iqtree.support.json)")
    return {"success": True, "tree": t,
            "support": {"type": meta.get("support_type", "unknown"),
                        "n_labeled": meta.get("n_labeled", 0),
                        "n_internal": meta.get("n_internal", 0),
                        "degraded": bool(meta.get("degraded")),
                        "fallback": bool(meta.get("fallback")),
                        "meta_file": os.path.join(phylo_dir, IQTREE_SUPPORT_META)}}



def run_stage_rtt(prep, aln_file, tree_file, tools, args, logger) -> Dict:
    """rtt: TreeTime-RTT + TreeDater-LTT

    ⚠️ 2026-09-15 (审查 P3) 状态标注: **本函数当前无任何调用点** ——
    对应能力已并入其它 stage (rtt/temporal/genes → clock; beast → run_phase_beast)。
    保留仅为历史参考与可能的 --stage 兼容扩展; 修改它**不会**影响在跑的管线。
    真正在用的入口见 process_virus() 的 stage 分发。
    """

    vname = prep["virus_name"]
    out_dir = prep["output_dir"]
    log = logger
    results = {"virus": vname, "success": False, "outputs": {}, "metrics": {}}

    # Step: TreeTime-RTT
    rtt_dir = os.path.join(out_dir, ph_dir('time'), "treetime_rtt")
    rtt_checkpoint = os.path.join(rtt_dir, "RootToTip_regression.pdf")
    if not args.force and os.path.exists(rtt_checkpoint):
        log.info(f"  [rtt] {vname}: TreeTime-RTT checkpoint exists, skip")
        results["outputs"]["rtt_plot"] = rtt_checkpoint
    else:
        log.info(f"  [rtt] {vname}: running TreeTime-RTT...")
        rtt = run_treetime_rtt(
            fasta_file=aln_file, tree_file=tree_file,
            dates_file=prep["dates_csv"], output_dir=rtt_dir,
            log=LogCollector(log))
        if rtt["success"]:
            results["metrics"]["beta"] = rtt["beta"]
            results["metrics"]["r_squared"] = rtt["r_squared"]
            results["metrics"]["p_value"] = rtt["p_value"]
            results["outputs"]["rtt_plot"] = rtt["plots"].get("rtt")
            results["outputs"]["timetree"] = rtt["output_tree"]
            results["outputs"]["timetree_plot"] = rtt["plots"].get("tree")
        else:
            log.warning(f"  [rtt] TreeTime-RTT: {rtt.get('error', 'unknown error')}")

    # Step: TreeDater-LTT
    dater_dir = os.path.join(out_dir, ph_dir('time'), "treedater_ltt")
    dater_checkpoint = os.path.join(dater_dir, "Phylogeny_dated.pdf")
    if not args.force and os.path.exists(dater_checkpoint):
        log.info(f"  [rtt] {vname}: TreeDater-LTT checkpoint exists, skip")
        results["outputs"]["dated_phylo"] = dater_checkpoint
    else:
        rscript = tools.get("Rscript") or "Rscript"
        seq_len = prep.get("avg_seq_len", 10000)
        dater_tree = results["outputs"].get("timetree") or tree_file
        log.info(f"  [rtt] {vname}: running TreeDater-LTT (seq_len={seq_len})...")
        dater = run_treedater_ltt(
            tree_file=dater_tree, metadata_file=prep["dates_csv"],
            seq_len=seq_len, output_dir=dater_dir, plot_ltt=True,
            rscript_path=rscript, threads=args.threads)
        if dater["success"]:
            results["outputs"]["dated_phylo"] = dater["plots"].get("phylogeny")
            results["outputs"]["ltt_plot"] = dater["plots"].get("ltt")
            results["metrics"]["treedater_output"] = dater["output"]
            # ── 2026-08-31 treedater 四功能 metrics ──
            _tdm = dater.get("metrics") or {}
            for _k, _lbl in [("cov_rate", "td_cov_rate"),
                             ("relaxed_clock_p", "td_relaxed_clock_p"),
                             ("tmrca", "td_tmrca"),
                             ("mean_rate", "td_mean_rate"),
                             ("n_outlier_tips", "td_n_outlier_tips")]:
                if _k in _tdm:
                    results["metrics"][_lbl] = _tdm[_k]
            if "tmrca_ci" in _tdm:
                results["metrics"]["td_tmrca_ci95"] = (
                    f"{_tdm['tmrca_ci'][0]:.3f}-{_tdm['tmrca_ci'][1]:.3f}")
            if "rate_ci" in _tdm:
                results["metrics"]["td_rate_ci95"] = (
                    f"{_tdm['rate_ci'][0]:.3e}-{_tdm['rate_ci'][1]:.3e}")
            if dater.get("outlier_tips"):
                results["outputs"]["treedater_outlier_tips"] = dater["outlier_tips"]
        else:
            log.warning(f"  [rtt] TreeDater-LTT: {dater.get('error', 'unknown error')}")

    results["success"] = tree_file is not None
    return results


def run_stage_saturation(prep, aln_file, args, logger) -> Dict:
    """saturation: 替换饱和分析 (C 值法 + Xia Iss) → saturation_out/"""
    vname = prep["virus_name"]
    out_dir = prep["output_dir"]
    results = {"virus": vname, "success": False, "outputs": {}, "metrics": {}}
    if not aln_file or not os.path.exists(aln_file):
        logger.warning(f"  [saturation] 无比对, 跳过")
        return results
    from saturation_analysis import run_saturation
    logger.info(f"  [saturation] {vname}: 替换饱和分析...")
    r = run_saturation(aln_file, outdir=os.path.join(out_dir, ph_dir('phylogeny'), "saturation_out"),
                       gap_treatment="pairwise", num_replicates=args.sat_replicates,
                       log=lambda m: logger.info(f"    {m}"))
    if r.get("success"):
        results["success"] = True
        results["outputs"]["saturation_report"] = r["report"]
        results["outputs"]["saturation_plot"] = r.get("plot")
        results["metrics"]["saturation_c"] = r["c_value"]
        results["metrics"]["saturation_iss"] = r["iss"]
        results["metrics"]["saturation_saturated"] = r["saturated"]
        logger.info(f"  [saturation] C={r['c_value']:.4f}, Iss={r['iss']:.4f}, 饱和={r['saturated']}")
    else:
        logger.warning(f"  [saturation] 失败: {r.get('error')}")
    return results


def run_stage_splitstree(prep, aln_file, args, logger) -> Dict:
    """splitstree: SplitsTree6 NeighborNet 分裂网络 (检验树状信号/网状进化)"""
    vname = prep["virus_name"]
    out_dir = prep["output_dir"]
    results = {"virus": vname, "success": False, "outputs": {}, "metrics": {}}
    if not aln_file or not os.path.exists(aln_file):
        logger.warning(f"  [splitstree] 无比对, 跳过")
        return results
    outdir = os.path.join(out_dir, ph_dir('phylogeny'), "splitstree")
    try:
        from splitstree_bridge import main as _bridge_main
    except ImportError:
        logger.warning("  [splitstree] splitstree_bridge 不可用, 跳过")
        return results
    logger.info(f"  [splitstree] {vname}: NeighborNet 分裂网络...")
    import contextlib, io as _io
    argv_backup = sys.argv
    sys.argv = ['splitstree_bridge.py', '-i', aln_file, '-o', outdir,
                '--network', 'neighbornet', '--max-taxa', str(getattr(args, 'splitstree_max_taxa', 200))]
    try:
        with contextlib.redirect_stdout(_io.StringIO()):
            rc = _bridge_main()
    except SystemExit as e:
        rc = e.code
    except Exception as e:
        rc = 1
        logger.warning(f"  [splitstree] 异常: {e}")
    finally:
        sys.argv = argv_backup
    nexus = os.path.join(outdir, 'neighbornet.nexus')
    if rc == 0 and os.path.exists(nexus) and os.path.getsize(nexus) > 100:
        results["success"] = True
        results["outputs"]["splitstree_nexus"] = nexus
        # 从 nexus 提取 splits 数作指标
        try:
            with open(nexus, errors='ignore') as f:
                for line in f:
                    if 'nsplits=' in line:
                        ns = int(line.split('nsplits=')[1].split()[0].rstrip(';'))
                        results["metrics"]["splitstree_nsplits"] = ns
                        logger.info(f"  [splitstree] nsplits={ns}")
                        break
        except Exception:
            pass
    else:
        logger.warning("  [splitstree] 失败或未产出")
    return results


# run_stage_concat 已移除 (2026-08-27, stage 删除; 工具函数在 multigene_tools.concatenate_alignments)

def run_stage_popgen(prep, aln_file, args, logger) -> Dict:
    """popgen: pypopart 群体遗传学 → popgen_out/popgen_report.txt"""
    vname = prep["virus_name"]
    out_dir = prep["output_dir"]
    popgen_dir = os.path.join(out_dir, ph_dir('popgen'))
    os.makedirs(popgen_dir, exist_ok=True)
    results = {"virus": vname, "success": False, "outputs": {}}

    if not aln_file or not os.path.exists(aln_file):
        logger.warning(f"  [popgen] 无比对, 跳过")
        return results

    cmd = [sys.executable, "-m", "popgen_analysis",
           "--fasta", aln_file,
           "--group", args.popgen_group,
           "--min-n", str(args.popgen_min_n),
           "--perm", str(args.popgen_perm)]
    if getattr(args, "popgen_exclude_reference", False):
        cmd.append("--exclude-reference")
    if getattr(args, "popgen_exclude_ids", None):
        cmd += ["--exclude-ids", str(args.popgen_exclude_ids)]
    meta = prep.get("meta_csv", "")
    if meta and os.path.exists(meta):
        cmd += ["--metadata", meta]

    logger.info(f"  [popgen] {vname}: pypopart 群体遗传 (group={args.popgen_group})...")
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, encoding='utf-8', errors='replace', cwd=str(SCRIPT_DIR), timeout=3600)
    except subprocess.TimeoutExpired:
        logger.warning("  [popgen] 超时 (>3600s)")
        return results

    # 成功才写报告 (失败时写 error 文件, 避免 checkpoint 误判完成)
    if r.returncode == 0:
        out_txt = os.path.join(popgen_dir, "popgen_report.txt")
        with open(out_txt, "w", encoding="utf-8") as f:
            f.write(r.stdout)
        results["outputs"]["popgen_report"] = out_txt
        results["success"] = True
        logger.info(f"  [popgen] 完成 → {out_txt}")
        # 解析关键指标进 metrics
        m = re.search(r"Tajima's D = ([0-9.\-]+)", r.stdout)
        if m:
            results["metrics"] = {"popgen_tajima_d": float(m.group(1))}
    else:
        logger.warning(f"  [popgen] 失败: {(r.stderr or r.stdout)[-300:]}")
    return results


def run_stage_capheine(prep, args, logger) -> Dict:
    """capheine: 正选择分析 (CAWLign→IQ-TREE→HyPhy→DRHIP→MultiQC) → capheine/"""
    vname = prep["virus_name"]
    out_dir = prep["output_dir"]
    capheine_dir = os.path.join(out_dir, ph_dir('select'), "capheine")
    os.makedirs(capheine_dir, exist_ok=True)
    results = {"virus": vname, "success": False, "outputs": {}, "metrics": {}}

    # ── 1. reference CDS: 用户提供 → datasets annotation_gb → work_dir gb/gbk → NCBI 下载 ──
    ref = args.capheine_ref
    ref_src = "--capheine_ref"
    if not ref and prep.get("annotation_gb") and os.path.exists(prep["annotation_gb"]):
        ref = prep["annotation_gb"]
        ref_src = "datasets.annotation_gb"
    if not ref:
        gbks = sorted(Path(out_dir).glob("*.gbk")) + sorted(Path(out_dir).glob("*.gb"))
        if not gbks:
            gbks = sorted(Path(out_dir, "data", "ncbi_ref").glob("*.gbk")) + \
                   sorted(Path(out_dir, "data", "ncbi_ref").glob("*.gb"))  # online 产物
        if gbks:
            ref = os.path.join(capheine_dir, "ref_cds.fasta")
            gx = subprocess.run(
                [sys.executable, "-m", "gbk_extractor", "-i", str(gbks[0]), "-n", ref],
                capture_output=True, text=True, encoding='utf-8', errors='replace', cwd=str(SCRIPT_DIR))
            if not (os.path.exists(ref) and os.path.getsize(ref) > 0):
                # viroid 无 CDS 判定 (2026-09-02): 与后置提取路径同款逻辑
                _has_cds = False
                try:
                    with open(gbks[0], errors="replace") as _gf:
                        _gb_head = _gf.read(4000)
                    _has_cds = " CDS " in _gb_head or "CDS\t" in _gb_head or "codon" in _gb_head.lower()
                except OSError:
                    pass
                if not _has_cds:
                    results["skip_reason"] = "viroid/no-CDS reference: 正选择 (密码子方法) 不适用"
                    logger.info(f"  [capheine] 参考 GBK 无 CDS (类病毒?), 记 skip: {gbks[0]}")
                else:
                    logger.warning(f"  [capheine] gbk_extractor 未产出 CDS 参考: {ref}")
                return results
            logger.info(f"  [capheine] 从 {gbks[0].name} 提取 CDS 参考 → ref_cds.fasta")
            ref_src = os.path.basename(gbks[0])
        elif prep.get("accession"):
            # 兜底: 无本地注释 → 按 accession 在线下载 (参考 analysis 的 efetch 逻辑)
            from utils.gb_fetch import get_or_fetch_cds_fasta
            ref = os.path.join(capheine_dir, "ref_cds.fasta")
            got = get_or_fetch_cds_fasta(prep["accession"],
                                         os.path.join(out_dir, ph_dir('data'), "ncbi_ref"), ref)
            if not got:
                logger.warning(f"  [capheine] 在线下载参考 gb 失败: {prep['accession']}")
                return results
            logger.info(f"  [capheine] 在线下载 NCBI 注释 → 提取 CDS 参考 → ref_cds.fasta")
            ref_src = f"NCBI:{prep['accession']}"
        else:
            logger.warning("  [capheine] 无本地注释/accession, 请传 --capheine_ref 或配置 datasets.annotation_gb")
            return results
    if not os.path.exists(ref):
        logger.warning(f"  [capheine] reference 不存在: {ref}")
        return results

    # annotation_gb / 下载的 gb 是 gb/gbk 文件时需先提取 CDS → ref_cds.fasta
    if not ref.endswith(".fasta") and not ref.endswith(".fa"):
        ref_fa = os.path.join(capheine_dir, "ref_cds.fasta")
        gx = subprocess.run(
            [sys.executable, "-m", "gbk_extractor", "-i", ref, "-n", ref_fa],
            capture_output=True, text=True, encoding='utf-8', errors='replace', cwd=str(SCRIPT_DIR))
        if os.path.exists(ref_fa) and os.path.getsize(ref_fa) > 0:
            ref = ref_fa
            logger.info(f"  [capheine] 从 {Path(ref).name} 提取 CDS 参考 → ref_cds.fasta")
        else:
            # 类病毒/无 CDS 参考 (2026-09-02): GBK 提不出 CDS 大多是 viroid (无编码区),
            # 密码子层面的正选择分析不适用 → 记 skip 而非 fail
            with open(ref, errors="replace") as _gf:
                _gb_head = _gf.read(4000)
            _has_cds = " CDS " in _gb_head or "CDS\u0009" in _gb_head or "codon" in _gb_head.lower()
            if not _has_cds:
                results["skip_reason"] = "viroid/no-CDS reference: 正选择 (密码子方法) 不适用"
                logger.info(f"  [capheine] 参考 GBK 无 CDS (类病毒?), 记 skip: {ref}")
            else:
                logger.warning(f"  [capheine] 提取 CDS 失败: {ref}")
            return results

    # ── 2. unaligned: 用户提供或 prep 序列 ──
    unaligned = args.capheine_unaligned or prep.get("combined_fasta", "")
    if not unaligned or not os.path.exists(unaligned):
        logger.warning(f"  [capheine] unaligned 序列不存在: {unaligned}")
        return results

    # ── 3. 调 capheine_pipeline.py ──
    cmd = [sys.executable, str(SCRIPT_DIR / "capheine_pipeline.py"),
           "--reference", ref, "--unaligned", unaligned,
           "--outdir", capheine_dir,
           "--workers", str(args.capheine_workers),
           "--cpus_iqtree", str(args.capheine_cpus_iqtree),
           "--cpus_hyphy", str(args.capheine_cpus_hyphy),
           "--code", str(args.capheine_code)]
    if args.capheine_mpi:
        cmd.append("--use_mpi")

    logger.info(f"  [capheine] {vname}: 正选择分析 (workers={args.capheine_workers})...")
    logger.info(f"    reference={Path(ref).name}, unaligned={Path(unaligned).name}")
    try:
        # 不设时限: 134+ 条全长序列的 HyPhy 五件套可能远超 2h, 由用户/集群层面控制 (用户拍板 2026-09-02)
        r = subprocess.run(cmd, capture_output=True, text=True, encoding='utf-8', errors='replace', cwd=str(SCRIPT_DIR))
    except subprocess.TimeoutExpired:
        logger.warning("  [capheine] 超时 (>7200s)")
        return results

    drhip_dir = os.path.join(capheine_dir, "drhip")
    if r.returncode == 0 and os.path.isdir(drhip_dir):
        results["success"] = True
        results["outputs"]["capheine_dir"] = capheine_dir
        results["outputs"]["capheine_drhip"] = drhip_dir
        results["metrics"]["capheine_genes"] = len(list(Path(drhip_dir).glob("*.csv")))
        logger.info(f"  [capheine] 完成 → {capheine_dir}/drhip")
    else:
        logger.warning(f"  [capheine] 失败: {(r.stderr or r.stdout)[-400:]}")
    return results


def run_stage_clock(prep, aln_file, tree_file, tools, args, logger,
                    results=None) -> Dict:
    """clock: 分子钟联合分析 (整体+分基因 RTT/LTT/DRT), 实现在 utils/clock_analysis.py,
    支持独立命令行运行; 此处仅做包装 (找 gbk/dates/工具路径 + 汇总 metrics)。"""
    out = {"outputs": {}, "metrics": {}}
    vname = prep["virus_name"]
    out_dir = prep["output_dir"]
    gbk = prep.get("annotation_gb")
    if gbk and not os.path.exists(gbk):
        gbk = None
    from clock_analysis import run_clock_analysis
    cr = run_clock_analysis(
        alignment=aln_file, tree=tree_file,
        dates_csv=prep.get("dates_csv", ""), out_dir=out_dir,
        meta_csv=prep.get("meta_csv", ""), gbk=gbk,
        seq_len=prep.get("avg_seq_len", 10000),
        threads=args.threads,
        iqtree_bin=tools.get("iqtree2") or tools.get("iqtree") or "iqtree2",
        rscript=tools.get("Rscript") or "Rscript",
        drt_randomizations=getattr(args, "drt_randomizations", 10),
        gene_drt_randomizations=10, gene_ltt=False,
        log=LogCollector(logger))
    out["outputs"].update(cr.get("outputs", {}))
    out["metrics"].update({k: v for k, v in cr.get("metrics", {}).items()
                            if v is not None})
    if results is not None and cr.get("summary_rows"):
        results["gene_results"] = [r for r in cr["summary_rows"] if r.get("level") == "gene"]
    out["success"] = bool(cr.get("summary_rows"))
    return out


def run_stage_temporal(prep, aln_file, tree_file, args, logger) -> Dict:
    """temporal: DRT (--check_temporal) + BETS (--bets)

    ⚠️ 2026-09-15 (审查 P3) 状态标注: **本函数当前无任何调用点** ——
    对应能力已并入其它 stage (rtt/temporal/genes → clock; beast → run_phase_beast)。
    保留仅为历史参考与可能的 --stage 兼容扩展; 修改它**不会**影响在跑的管线。
    真正在用的入口见 process_virus() 的 stage 分发。
    """

    vname = prep["virus_name"]
    out_dir = prep["output_dir"]
    log = logger
    results = {"virus": vname, "success": True, "metrics": {}}
    dates_csv = prep.get("dates_csv", "")
    meta_csv = prep.get("meta_csv", "")

    if args.check_temporal and aln_file and os.path.exists(aln_file) and dates_csv and os.path.exists(dates_csv):
        temporal_dir = os.path.join(out_dir, ph_dir('time'), "temporal_signal")
        log.info(f"  [temporal] {vname}: DRT temporal signal check...")
        drt = date_randomization_test(
            fasta_file=aln_file, dates_csv=dates_csv,
            output_dir=temporal_dir, n_randomizations=args.drt_randomizations,
            tree_file=tree_file, log=LogCollector(log))
        results["metrics"]["temporal_signal"] = drt.get("conclusion", "")
        # 2026-09-15 (审查 P2-7): 只有 DRT 真正产出结论才写 bool; 崩溃时写 None
        # ("没算出来" ≠ "没通过")
        _drt_ok = bool(drt.get("success"))
        results["metrics"]["drt_passed"] = drt.get("passed") if _drt_ok else None
        results["metrics"]["drt_status"] = "ok" if _drt_ok else "error"
        if not _drt_ok:
            results["metrics"]["drt_error"] = drt.get("error") or "(无错误信息)"

    if args.bets and aln_file and os.path.exists(aln_file) and meta_csv and os.path.exists(meta_csv):
        # 2026-08-30 审计定性: BETS/GSS 路线已归档 (见 phylo_results/BETS_GSS_ARCHIVE.md)。
        # 历史坑: 旧版 execute=False 只生成 XML、log_bf=None 却仍写 metrics → temporal 被
        # _mark 假完成 (checkpoint 落盘, 无任何 BF 结果)。现改为显式禁用哨兵 + 不写 metrics。
        log.warning(f"  [temporal] {vname}: BETS 已归档止损 (GSS 注入器发散, 见 "
                    f"phylo_results/BETS_GSS_ARCHIVE.md), 本次不执行。如需重启用 --bets 时手动取消哨兵")
    return results


def run_stage_beast(prep, phase_results, args, tools, logger) -> Dict:
    """beast: 全长 BEAST 定年 (Phase 5)

    ⚠️ 2026-09-15 (审查 P3) 状态标注: **本函数当前无任何调用点** ——
    对应能力已并入其它 stage (rtt/temporal/genes → clock; beast → run_phase_beast)。
    保留仅为历史参考与可能的 --stage 兼容扩展; 修改它**不会**影响在跑的管线。
    真正在用的入口见 process_virus() 的 stage 分发。
    """

    return run_phase_beast(prep, phase_results, args, tools, logger)


def run_stage_gene_dating(prep, results, args, logger) -> None:
    """gene_dating: 分基因 BEAST 定年 (gene_partition_dating.py, 吃 capheine 产物)"""
    out_dir = prep["output_dir"]
    logger.info(f"  [gene_dating] {prep['virus_name']}: per-gene BEAST dating...")
    genes_dir = args.genes_dir
    if not genes_dir or not os.path.exists(genes_dir):
        for cand in [os.path.join(out_dir, ph_dir('select'), "capheine", "cawlign"),
                     os.path.join(out_dir, ph_dir('select'), "capheine"),
                     os.path.join(out_dir, ph_dir('time'), "gene_split")]:
            if os.path.isdir(cand):
                genes_dir = cand
                break
    if not genes_dir or not os.path.exists(genes_dir):
        # fallback: capheine 产物缺失时自动生成基因分区比对 (cawlign)
        logger.info("  [gene_dating] 未找到 capheine 产物, 自动生成基因比对 (cawlign fallback)...")
        try:
            from gene_align_generator import generate_gene_alignments, \
                find_reference_cds, find_unaligned_fasta
            gene_split = os.path.join(out_dir, ph_dir('time'), "gene_split")
            ref = find_reference_cds(Path(out_dir), args.capheine_ref or prep.get("annotation_gb"),
                                     accession=prep.get("accession"))
            # 优先用 prep 的原始序列 (capheine 同源), 其次自动找
            unal = None
            if prep.get("combined_fasta") and os.path.exists(prep["combined_fasta"]):
                unal = Path(prep["combined_fasta"])
            else:
                unal = find_unaligned_fasta(Path(out_dir), args.alignment)
            if ref and unal:
                generated = generate_gene_alignments(ref, unal, Path(gene_split), args.capheine_code)
                if generated:
                    genes_dir = gene_split
                    logger.info(f"  [gene_dating] 自动生成 {len(generated)} 个基因比对 → {gene_split}")
                else:
                    logger.warning("  [gene_dating] fallback 生成 0 个基因比对, 跳过分基因定年")
            else:
                logger.warning(f"  [gene_dating] fallback 缺参考 CDS (ref={ref}) 或未比对序列 (unal={unal}), 跳过")
        except Exception as e:
            logger.warning(f"  [gene_dating] fallback 生成失败: {e}")
    if not genes_dir or not os.path.exists(genes_dir):
        logger.warning("  [gene_dating] 未找到基因比对目录 (传 --genes_dir 或先跑 capheine)。跳过分基因定年。")
        return

    # 空目录防御 (2026-09-02, PSTVd 教训): capheine 失败后目录存在但无基因比对,
    # 旧逻辑直接部署空任务并写 checkpoint。必须验证 alignments 真的有文件。
    def _count_gene_aligns(gd):
        n = 0
        for sub in ["cawlign/alignments", "alignments", "gene_split", ""]:
            d = os.path.join(gd, sub) if sub else gd
            if os.path.isdir(d):
                n = sum(1 for f in os.listdir(d) if f.endswith((".fasta", ".fa", ".aln")))
                if n:
                    break
        return n
    n_aligns = _count_gene_aligns(genes_dir)
    if n_aligns == 0:
        # viroid 无 CDS → 分基因定年不适用, 记 skip 语义 (返回不写 checkpoint)
        logger.info("  [gene_dating] 基因比对目录为空 (viroid 无 CDS 或 capheine 未产出), 分基因定年不适用, 跳过")
        return

    gene_dating_dir = os.path.join(out_dir, ph_dir('time'), "gene_dating")
    meta_for_gene = prep["meta_csv"] if (prep.get("meta_csv") and os.path.exists(prep["meta_csv"])) else prep["dates_csv"]
    full_tmrca = results.get("metrics", {}).get("beast_tmrca")
    gd_cmd = [sys.executable, str(SCRIPT_DIR / "gene_partition_dating.py"),
              "--metadata", meta_for_gene, "--out", gene_dating_dir,
              "--genes-dir", genes_dir,
              "--chain", str(args.gene_chain),
              "--chain-big", str(args.gene_chain_big),
              "--big-genes", args.big_genes,
              "--threads", str(args.threads)]
    if full_tmrca:
        gd_cmd += ["--full-tmrca", str(full_tmrca)]
    logger.info(f"  [gene_dating] 提交分基因 BEAST (genes-dir={genes_dir})...")
    try:
        # 无时限 (2026-09-02 拆除 900s 部署上限): 分基因 BEAST 部署本身可能慢, 由集群层面控制
        gdr = subprocess.run(gd_cmd, capture_output=True, text=True, encoding='utf-8', errors='replace', cwd=str(SCRIPT_DIR))
    except subprocess.TimeoutExpired:
        logger.warning("  [gene_dating] 部署超时 (>900s)")
        return
    if gdr.returncode == 0:
        results["outputs"]["gene_dating_dir"] = gene_dating_dir
        results["metrics"]["gene_dating_genes_dir"] = genes_dir
        logger.info(f"  [gene_dating] 部署完成 → {gene_dating_dir} (后台 BEAST 逐基因运行)")
    else:
        logger.warning(f"  [gene_dating] 失败: {(gdr.stderr or gdr.stdout)[-400:]}")


def _find_phylogeo_mcc(out_dir: str) -> Optional[str]:
    """定位 phylogeo 产出的 MCC 树。

    2026-09-16 抽成函数: 原先 phylogeo 阶段内联写了一遍 glob, geo 阶段要用同一份
    候选顺序 (优先 `merged/mcc.tree`, 再退化为 `*mcc*.tree`); 两处各写一份必漂移。
    """
    phylogeo_dir = os.path.join(out_dir, ph_dir('geography'), "phylogeography")
    cands = list(Path(phylogeo_dir).glob("merged/mcc.tree")) \
        + list(Path(phylogeo_dir).glob("*mcc*.tree"))
    return str(cands[0]) if cands else None


def run_stage_geo(prep, results, aln_file, tree_file, args, logger) -> None:
    """geo: 地理分析 (Mantel / 树地理 / VirSpaceTime)"""
    out_dir = prep["output_dir"]
    meta_csv = prep["meta_csv"]
    if not prep.get("has_location"):
        return
    geo_dir = os.path.join(out_dir, ph_dir('geography'), "geo_analysis")
    os.makedirs(geo_dir, exist_ok=True)

    # 前处理: 地理均衡抽样 (--geo_subsample N, 输出子集到 geo 目录, 不覆盖全量)
    geo_aln = aln_file
    if getattr(args, "geo_subsample", 0) and aln_file and os.path.exists(aln_file):
        try:
            from utils.geo_subsampler import subsample_fasta
            sub_dir = os.path.join(geo_dir, "subsampled")
            os.makedirs(sub_dir, exist_ok=True)
            r = subsample_fasta(aln_file, num_seqs=args.geo_subsample,
                                output_dir=sub_dir)
            # 历史坑: 模块默认输出 extract.fas, 调用处却检查 subsampled.fasta → 永远失效
            out_fa = (r or {}).get("extract") or os.path.join(sub_dir, "subsampled.fasta")
            if os.path.exists(out_fa):
                geo_aln = out_fa
                # 2026-09-15: 记录实际种子, 便于复现同一子集
                results["metrics"]["geo_subsample_seed"] = (r or {}).get("seed", "")
                logger.info(f"  [geo] 地理均衡抽样: {args.geo_subsample} 条 → "
                            f"{Path(out_fa).name} (seed={(r or {}).get('seed')})")
        except Exception as e:
            logger.warning(f"  [geo] 子抽样跳过: {e}")

    logger.info(f"  [geo] {prep['virus_name']}: geographic analysis...")
    tree = results.get("outputs", {}).get("tree") or tree_file
    if geo_aln and os.path.exists(geo_aln):
        mantel = run_mantel_test(geo_aln, meta_csv, geo_dir, log=LogCollector(logger))
        if mantel["success"]:
            results["outputs"]["mantel_plot"] = mantel["plot"]
            results["metrics"]["mantel_r"] = mantel["mantel_r"]
    if tree and os.path.exists(tree) and geo_aln and os.path.exists(geo_aln):
        geo_tree = run_tree_geography(tree, meta_csv, geo_aln, geo_dir, log=LogCollector(logger))
        if geo_tree["success"]:
            results["outputs"]["geo_tree"] = geo_tree["plot"]
    vst = run_virspacetime(meta_csv, geo_dir, log=LogCollector(logger))
    if vst["success"]:
        results["outputs"]["spacetime_plot"] = vst["plot"]

    # 后处理: geo_analysis 扩展 (智能子抽样 + RRT 距离矩阵法 [+ RRT 树注释法])
    if aln_file and os.path.exists(aln_file) and meta_csv and os.path.exists(meta_csv):
        extras = ["--fasta", aln_file, "--metadata", meta_csv,
                  "--outdir", os.path.join(geo_dir, "extended")]
        # ── RRT 树注释法 (VirPhyKit 同口径, 2026-09-16 接线; 默认关闭) ──
        # 只在用户显式给了 --rrt_randomized_dir 时才跑; 注意 geo 阶段在 STAGE_ORDER 里
        # 排在 phylogeo **之后**(两者都只依赖 tree), 所以此时 MCC 树通常已合并好。
        rrt_rand = getattr(args, "rrt_randomized_dir", None)
        if rrt_rand:
            mcc = _find_phylogeo_mcc(out_dir)
            if not mcc:
                logger.warning(
                    "  [geo] 已给 --rrt_randomized_dir, 但未找到 phylogeo MCC 树 "
                    "(需 geography/phylogeography/**/mcc.tree) → 跳过 RRT 树注释法。"
                    "请先跑 phylogeo 阶段, 或改用 --stage phylogeo,geo 让两者同批完成")
            else:
                # 注意: 不传 --no-rrt —— 距离矩阵法(RRT 方法一)是既有结论, 不能被
                # "新加了个参数"顺手关掉。两种方法口径不同, 各自独立报告。
                extras += ["--rrt-tree", mcc, "--rrt-randomized", str(rrt_rand)]
                logger.info(f"  [geo] RRT 树注释法 (VirPhyKit 口径): MCC={os.path.basename(mcc)}, "
                            f"随机化树源={rrt_rand}")
        run_postprocess("geo_analysis", logger, *extras)
        # 把 RRT 树注释法的结论回读进 results (供 report 阶段汇总)
        _rrt_sum = os.path.join(geo_dir, "extended", "rrt_tree_summary.json")
        if os.path.exists(_rrt_sum):
            try:
                import json as _json
                with open(_rrt_sum, encoding="utf-8") as _f:
                    _d = _json.load(_f)
                if _d.get("success"):
                    results["metrics"]["rrt_tree_passed"] = _d.get("passed")
                    results["metrics"]["rrt_tree_target"] = _d.get("target_state")
                    results["metrics"]["rrt_tree_margin"] = _d.get("margin")
                    results["outputs"]["rrt_tree_summary"] = _rrt_sum
                    if _d.get("plot"):
                        results["outputs"]["rrt_tree_plot"] = _d["plot"]
            except Exception as e:
                logger.warning(f"  [geo] 读取 rrt_tree_summary.json 失败: {e}")


def run_stage_geo_paths(prep, results, aln_file, tree_file, args, logger) -> None:
    """geo_paths: 事件级传播路径 (pathways/episode/三级置信 + LTL + 扩散统计 [+ 传播图])

    2026-09-16 新增子阶段 (吸收 phymapr/seraphim/ggphylogeo/flu_d_project 方法, 出处见
    utils/transmission_paths.py 与 utils/transmission_map.py docstring)。

    输入:  phylogeo 阶段的 merged/mcc.tree (STAGE_CHECKPOINT_FILES['phylogeo'])
    输出:  geography/geo_analysis/pathways/{branch_table,pathways,episodes,ltl_summary}.csv
           geography/geo_analysis/pathways/pathway_qc.json   ← stage checkpoint
           geography/geo_analysis/dispersion/{velocity_table,wavefront_series,summary}
           [可选 --geo_map] geography/visualization/transmission_map.png(+gif)

    定位:  MCC **描述层** (与 RRT/TEMPMIG 并列), 不改动 BSSVS BF 推断结论;
           阈值 direct/indirect 期望替换数默认沿用 phymapr (2/5), 进论文前须按
           本数据 E[S] 分布标定 (见 INTEGRATION_20260916.md §2)。
    """
    out_dir = prep["output_dir"]
    meta_csv = prep["meta_csv"]
    if not prep.get("has_location"):
        return
    if getattr(args, "no_pathways", False):
        logger.info("  [geo_paths] --no-pathways 显式关闭, 跳过")
        return
    mcc = _find_phylogeo_mcc(out_dir)
    if not mcc:
        logger.warning(
            "  [geo_paths] 未找到 phylogeo MCC 树 (geography/phylogeography/**/mcc.tree) "
            "→ 跳过。请先跑 phylogeo 阶段 (其 checkpoint 即 merged/mcc.tree)")
        return

    geo_dir = os.path.join(out_dir, ph_dir('geography'), "geo_analysis")
    os.makedirs(geo_dir, exist_ok=True)
    from utils.transmission_paths import run_pathways

    logger.info(f"  [geo_paths] {prep['virus_name']}: event-level pathways "
                f"(MCC={os.path.basename(mcc)})")
    pr = run_pathways(
        mcc, meta_csv, os.path.join(geo_dir, "pathways"),
        fasta=aln_file,
        direct_snp=getattr(args, "direct_snp", 2.0),
        indirect_snp=getattr(args, "indirect_snp", 5.0),
        censor_years=getattr(args, "censor_years", 0.5),
        calibrate=getattr(args, "calibrate", None),
        log=LogCollector(logger),
    )
    if not pr.get("success"):
        logger.warning(f"  [geo_paths] pathways 失败: {pr.get('error')}")
        return

    results["outputs"]["pathway_qc"] = os.path.join(geo_dir, "pathways", "pathway_qc.json")
    results["outputs"]["pathways_csv"] = pr["pathways"]
    results["metrics"]["n_pathways"] = pr["n_pathways"]
    results["metrics"]["n_episodes"] = pr["n_episodes"]
    qc = pr.get("qc_data") or {}
    conf = qc.get("confidence_counts") or {}
    for k in ("direct", "indirect", "distant_import", "unclassified"):
        if k in conf:
            results["metrics"][f"pathway_conf_{k}"] = conf[k]

    # 扩散统计 (默认开; 坐标走 geo_resolver 离线词典, 解析失败则跳过不致命)
    try:
        from utils.diffusion_stats import run_dispersion
        dr = run_dispersion(
            pr["branch_table"],
            output_dir=os.path.join(geo_dir, "dispersion"),
            n_perm=getattr(args, "n_perm", 200),
            log=LogCollector(logger),
        )
        if dr.get("success"):
            results["outputs"]["dispersion_summary"] = dr["summary"]
            # 关键数字回填 metrics → report 汇总列 (2026-09-16)
            sd = dr.get("summary_data") or {}
            results["metrics"]["geo_mean_velocity_km_yr"] = sd.get("mean_velocity_km_yr")
            results["metrics"]["geo_max_wavefront_km"] = sd.get("max_wavefront_km")
            _null = sd.get("permutation_null") or {}
            results["metrics"]["geo_perm_p_displacement"] = _null.get("p_displacement")
    except SystemExit as e:   # resolve_coordinates 缺坐标时 raise SystemExit
        logger.warning(f"  [geo_paths] 扩散统计跳过: {e}")
    except Exception as e:    # noqa: BLE001
        logger.warning(f"  [geo_paths] 扩散统计失败: {e}")

    # 静态传播图 (+可选 GIF) — 写到 visualization/, 与 spread3 产物并列
    if getattr(args, "geo_map", False):
        try:
            from utils.transmission_map import draw_animation, draw_static_map
            import pandas as _pd
            bt = _pd.read_csv(pr["branch_table"])
            vis_dir = os.path.join(out_dir, ph_dir('geography'), "visualization")
            os.makedirs(vis_dir, exist_ok=True)
            locs = sorted(set(bt["start_location"].dropna()) | set(bt["end_location"].dropna()))
            from utils.transmission_map import resolve_coordinates
            coords = resolve_coordinates(locs, log=LogCollector(logger))
            missing = [l for l in locs if coords.get(l) is None]
            if missing:
                logger.warning(f"  [geo_paths] 传播图跳过: 地点缺坐标 {missing} "
                               f"(可用 geo_analysis CLI --coords 提供)")
            else:
                png = draw_static_map(bt, coords, os.path.join(vis_dir, "transmission_map.png"),
                                      log=LogCollector(logger))
                results["outputs"]["transmission_map"] = png
                if getattr(args, "geo_gif", False):
                    gif = draw_animation(bt, coords,
                                         os.path.join(vis_dir, "transmission_map.gif"),
                                         n_frames=getattr(args, "geo_gif_frames", 36),
                                         log=LogCollector(logger))
                    results["outputs"]["transmission_gif"] = gif
        except Exception as e:  # noqa: BLE001
            logger.warning(f"  [geo_paths] 传播图失败 (不影响 stage): {e}")

    # 后验 Markov jumps (SpreaD3 同名统计; MCC 单树计数的后验校正; 2026-09-17 新增)
    try:
        from utils.diffusion_stats import find_posterior_trees, markov_jumps_posterior
        _pt = find_posterior_trees(mcc)
        if _pt:
            mj = markov_jumps_posterior(_pt, log=LogCollector(logger))
            if mj.get("success"):
                mj_path = os.path.join(geo_dir, "markov_jumps.json")
                with open(mj_path, "w", encoding="utf-8") as _f:
                    json.dump(mj, _f, ensure_ascii=False, indent=2)
                results["outputs"]["markov_jumps"] = mj_path
                results["metrics"]["geo_markov_jumps_mean"] = mj["mean"]
                results["metrics"]["geo_markov_jumps_hpd95"] = str(mj["hpd95"])
        else:
            logger.warning("  [geo_paths] 未找到后验树文件 → 跳过 Markov jumps")
    except Exception as e:  # noqa: BLE001
        logger.warning(f"  [geo_paths] Markov jumps 失败 (不影响 stage): {e}")

    # 双流一致性闸门 (时间流 RTT/DRT vs 地理流 MCC 根高/速率; 2026-09-17 新增)
    try:
        from utils.transmission_paths import timescale_consistency
        clk = os.path.join(out_dir, ph_dir('time'), "clock_summary.tsv")
        tc_path = os.path.join(geo_dir, "timescale_consistency.json")
        if os.path.exists(clk):
            tc = timescale_consistency(mcc, clk, tc_path, log=LogCollector(logger))
            results["outputs"]["timescale_consistency"] = tc_path
            results["metrics"]["geo_root_height_yr"] = (tc.get("geo_stream") or {}).get("root_height_median_yr")
            results["metrics"]["timescale_rate_ratio"] = tc.get("rate_ratio_geo_vs_time")
            results["metrics"]["timescale_flag"] = tc.get("flag")
        else:
            logger.warning("  [geo_paths] 双流一致性: 无 time/clock_summary.tsv (时间流未跑) → 跳过")
    except Exception as e:  # noqa: BLE001
        logger.warning(f"  [geo_paths] 双流一致性检查失败 (不影响 stage): {e}")

    # 树-图联动双面板 (phymapr Tree|Map 图版的 matplotlib 替代; 2026-09-16 补)
    if getattr(args, "geo_tree_map", False):
        try:
            from utils.transmission_map import draw_tree_map
            draw_tree_map(mcc, prep["meta_csv"],
                          coords=None, out_path=os.path.join(
                              out_dir, ph_dir('geography'), "visualization",
                              "transmission_tree_map.png"),
                          pathways_csv=os.path.join(geo_dir, "pathways", "pathways.csv"),
                          log=LogCollector(logger))
            results["outputs"]["tree_map"] = os.path.join(
                out_dir, ph_dir('geography'), "visualization", "transmission_tree_map.png")
        except SystemExit as e:
            logger.warning(f"  [geo_paths] 树图联动跳过: {e}")
        except Exception as e:  # noqa: BLE001
            logger.warning(f"  [geo_paths] 树图联动失败 (不影响 stage): {e}")


def _host_channel_rows(path: str) -> "list":
    """读宿主通道 → [(name, host), ...]。文件不存在/无该列 → []。

    抽成函数是为了让「有没有宿主」的判据只有一处实现 (与写入端
    `has_host()` 同源), 避免"写表时算一次、分析时又算一遍"的漂移。
    """
    import csv as _csv
    from utils.metadata_governance import has_host, normalize_host
    if not (path and os.path.exists(path)):
        return []
    try:
        with open(path, encoding="utf-8-sig", errors="replace", newline="") as f:
            head = f.readline()
            f.seek(0)
            delim = "\t" if head.count("\t") > head.count(",") else ","
            rdr = _csv.DictReader(f, delimiter=delim)
            ncol = None
            for c in (rdr.fieldnames or []):
                if (c or "").strip().lower() in ("host", "host_species",
                                                 "host_scientific_name"):
                    ncol = c
                    break
            if not ncol:
                return []
            out = []
            for r in rdr:
                nm = str(r.get("name", "") or "").strip()
                hv = normalize_host(r.get(ncol))
                if nm and has_host(hv):
                    out.append((nm, hv))
            return out
    except OSError:
        return []


def run_stage_host(prep, results, aln_file, tree_file, args, logger) -> None:
    """host: 宿主分化分析

    2026-09-16 第三轮 (老师: 「host 也独立出来, 仅对有 host 的进行分析」):
    输入由**地理通道**改为**宿主通道** `data/host.csv`。该表只收有宿主的样品,
    「仅对有 host 的进行分析」因此是**结构性保证**, 而不是靠每个分析脚本各自过滤。
    """
    out_dir = prep["output_dir"]
    host_csv = prep.get("host_csv") or ""
    if not (host_csv and os.path.exists(host_csv)):
        logger.info("  [host] 无宿主通道 (host.csv), 跳过 (无非空宿主信息)")
        return
    if not (aln_file and os.path.exists(aln_file)):
        return
    rows = _host_channel_rows(host_csv)
    n_hosts_total = len({h for _, h in rows})
    if not rows:
        logger.info(f"  [host] 宿主通道为空 ({os.path.basename(host_csv)}), "
                    f"跳过 (没有带宿主的样品)")
        return
    logger.info(f"  [host] 宿主通道: {len(rows)} 个带宿主的样品 / "
                f"{n_hosts_total} 个宿主")
    host_dir = os.path.join(out_dir, ph_dir('popgen'), "host_analysis")
    os.makedirs(host_dir, exist_ok=True)
    logger.info(f"  [host] {prep['virus_name']}: host differentiation (仅宿主子集)...")
    results["metrics"]["n_host_samples"] = len(rows)
    results["metrics"]["n_hosts_expected"] = n_hosts_total
    if tree_file and os.path.exists(tree_file):
        h_tree = run_tree_host(tree_file, host_csv, aln_file, host_dir,
                               log=LogCollector(logger))
        if h_tree["success"]:
            results["outputs"]["host_tree"] = h_tree["plot"]
            results["metrics"]["n_hosts"] = h_tree["n_hosts"]
    h_diff = run_host_differentiation(aln_file, host_csv, host_dir,
                                      log=LogCollector(logger))
    if h_diff["success"]:
        results["outputs"]["host_diff_plot"] = h_diff["plot"]
        results["metrics"]["host_between_mean"] = h_diff["between_mean"]
        results["metrics"]["host_within_mean"] = h_diff["within_mean"]
        results["metrics"]["host_p_value"] = h_diff["p_value"]


def run_stage_genes(prep, results, aln_file, tools, args, logger) -> None:
    """genes: 分基因 TreeTime-RTT (gbk 坐标切分)

    ⚠️ 2026-09-15 (审查 P3) 状态标注: **本函数当前无任何调用点** ——
    对应能力已并入其它 stage (rtt/temporal/genes → clock; beast → run_phase_beast)。
    保留仅为历史参考与可能的 --stage 兼容扩展; 修改它**不会**影响在跑的管线。
    真正在用的入口见 process_virus() 的 stage 分发。
    """

    out_dir = prep["output_dir"]
    if not aln_file or not os.path.exists(aln_file):
        return
    # 兼容 .gbk 和 .gb 后缀 (datasets.annotation_gb 常用 .gb)
    gbk_files = sorted(set(Path(out_dir).glob("*.gbk")) | set(Path(out_dir).glob("*.gb")))
    gbk_files = [f for f in gbk_files if f.is_file()]
    if not gbk_files:
        return
    gene_dir = os.path.join(out_dir, ph_dir('phylogeny'), "gene_analysis")
    logger.info(f"  [genes] {prep['virus_name']}: per-gene TreeTime-RTT...")
    gene_alns = slice_genes_from_alignment(
        str(gbk_files[0]), aln_file, os.path.join(gene_dir, "alignments"),
        log=LogCollector(logger))
    if gene_alns:
        gene_results = run_gene_rtt_batch(
            gene_alns, prep["meta_csv"], os.path.join(gene_dir, "trees"),
            iqtree_bin=tools.get("iqtree2") or tools.get("iqtree") or "iqtree2",
            threads=args.threads, log=LogCollector(logger))
        results["gene_results"] = gene_results
        for gr in gene_results:
            results["metrics"][f"R2_{gr['gene']}"] = gr.get('r_squared')

        # 2026-08-30 新增: 逐基因 DRT 时间信号检验 (自动化)。
        # 动机: 全基因组 DRT FAILED 会掩盏慢进化基因 (如 CP) 的真实信号;
        # 分基因 RTT 只有 R² 无随机化对照 → 无法判定单基因时间信号真伪。
        # 每基因随机化次数降为 10 (控制耗时; 全基因组级仍用 --drt_randomizations)。
        # 检验不阻断 genes stage; 通过与否均诚实落盘 metrics。
        dates_csv_g = prep.get("dates_csv", "")
        if dates_csv_g and os.path.exists(dates_csv_g) and getattr(args, 'check_temporal', False):
            from utils.temporal_signal import date_randomization_test
            drt_dir = os.path.join(gene_dir, "drt")
            for gr in gene_results:
                gname = gr.get('gene', '')
                g_aln = os.path.join(gene_dir, "alignments", f"{gname}.fasta")
                if not (gname and os.path.exists(g_aln)):
                    continue
                n_seq = sum(1 for l in open(g_aln) if l.startswith('>'))
                if n_seq < 10 or n_seq > 500:
                    logger.info(f"  [genes] DRT skip {gname}: 序列数 {n_seq} 超出 [10,500] 合理区间")
                    continue
                g_tree = os.path.join(gene_dir, "trees", f"{gname}.treefile")
                try:
                    gd = date_randomization_test(
                        fasta_file=g_aln, dates_csv=dates_csv_g,
                        output_dir=os.path.join(drt_dir, gname),
                        n_randomizations=10,
                        tree_file=g_tree if os.path.exists(g_tree) else None,
                        log=LogCollector(logger))
                    if gd.get("conclusion"):
                        _gd_ok = bool(gd.get("success"))
                        results["metrics"][f"drt_{gname}_passed"] = (
                            gd.get("passed") if _gd_ok else None)
                        results["metrics"][f"drt_{gname}_status"] = (
                            "ok" if _gd_ok else "error")
                        results["metrics"][f"drt_{gname}_pct"] = gd.get("r2_percentile")
                        logger.info(f"  [genes] DRT {gname}: passed={gd.get('passed')} "
                                    f"pct={gd.get('r2_percentile')}")
                except Exception as e:
                    logger.warning(f"  [genes] DRT {gname} 失败 (不阻断): {e}")
        elif not getattr(args, 'check_temporal', False):
            logger.info("  [genes] 逐基因 DRT 未开启 (--check_temporal 未加), 跳过")


def _wait_chains_and_merge(work_dir, logger, label="beast", merge_fn=None,
                           poll_interval=60, stale_minutes=30):
    """多链前台等待闭环: 轮询链进度直到全部完成 → 合并 (logcombiner+treeannotator)。

    - 每 poll_interval 秒读 run_status.json + 各链 log 末行 state, 打印进度
    - 链全部完成 (all_done) 时调 merge_fn(work_dir, st) 合并并返回 True
    - 无 run_status.json / 链死亡 (log 超过 stale_minutes 无增长且无 beast 进程) → 返回 False 不阻断
    - 用户拍板 (2026-09-02): 多链提交后前台等待闭环, 不再后台挂起等二次重跑
    """
    from utils.merge_results import load_task, all_done, read_states, merge as _default_merge
    st = load_task(work_dir)
    if not st:
        logger.warning(f"  [{label}] 无 run_status.json, 无法等待 (视为后台模式)")
        return False
    merge_fn = merge_fn or (lambda wd, s: _default_merge(wd, s, background=False))
    prefixes = st.get("prefixes", [])
    chain_len = st.get("chain_length", 0) or 1
    logger.info(f"  [{label}] 前台等待 {len(prefixes)} 链完成 (chain_length={chain_len/1e6:.0f}M, "
                f"每 {poll_interval}s 轮询; Ctrl-C 不影响后台链)")
    import time as _time, subprocess as _sp
    last_sizes = {}
    last_progress_time = _time.time()
    while True:
        _time.sleep(poll_interval)
        states = read_states(work_dir, prefixes)
        done_n = sum(1 for s in states if s >= chain_len * 0.999)
        prog = ", ".join(f"{p}: {s/1e6:.2f}M" for p, s in zip(prefixes, states))
        r = _sp.run('ps aux | grep "[b]east" | grep java | wc -l',
                    shell=True, capture_output=True, text=True)
        n_proc = int(r.stdout.strip() or 0)
        logger.info(f"  [{label}] 进度 [{prog}] 完成 {done_n}/{len(prefixes)}, 进程 {n_proc}")
        if all_done(work_dir, st):
            logger.info(f"  [{label}] 全部链完成, 合并 (logcombiner + treeannotator)...")
            merge_fn(work_dir, st)
            return True
        # 停滞检测: 所有 log 大小不变 且 无 beast 进程 → 链已死
        cur_sizes = {p: os.path.getsize(os.path.join(work_dir, f"{p}.log"))
                     if os.path.exists(os.path.join(work_dir, f"{p}.log")) else 0
                     for p in prefixes}
        grew = any(cur_sizes.get(p, 0) > last_sizes.get(p, 0) for p in prefixes)
        if grew:
            last_progress_time = _time.time()
        last_sizes = cur_sizes
        if not grew and n_proc == 0:
            # 再给一次机会: 可能刚好在 all_done 边界
            if all_done(work_dir, st):
                merge_fn(work_dir, st)
                return True
            logger.warning(f"  [{label}] 链停滞且无进程 (疑似死亡), 放弃等待")
            return False
        if _time.time() - last_progress_time > stale_minutes * 60:
            logger.warning(f"  [{label}] 链超过 {stale_minutes} 分钟无进展, 放弃等待")
            return False


def _meta_locations(meta_csv):
    """从 metadata CSV 提取去重后的 location 列表 (XML 复用分支给 SpreaD3 用)

    注意: 返回的是 **原始** location 名 (CSV 顺序), 而 XML 里实际写入的是
    `_simplify_location()` 后的省名 (sorted 顺序)。两者不等价 —— 只做兜底用,
    优先走 `_locations_from_xml()`。
    """
    locs = []
    if meta_csv and os.path.exists(meta_csv):
        import csv as _csv
        try:
            with open(meta_csv, encoding='utf-8-sig') as f:
                for row in _csv.DictReader(f):
                    loc = (row.get('location') or '').strip()
                    if loc and loc not in locs:
                        locs.append(loc)
        except Exception:
            pass
    return locs


def _locations_from_xml(xml_path):
    """从已生成的 BEAST1 phylogeo XML 里读出 location 状态列表 (顺序即声明序)。

    2026-09-15 (审查 P2-4) 修复: 复用已有 XML 的分支原先调 `_meta_locations()`,
    拿到的是 **原始名 + CSV 顺序**; 而新建 XML 的分支用生成器返回的
    `locations`(**简化省名 + sorted 顺序**)。两条路径给 SpreaD3 / migration_bf
    的 locations 口径不同 → BF 索引映射与地图节点张冠李戴。
    正确做法是直接读 XML 里真正声明的那份 (generator 写的就是它):

        <generalDataType id="location.dataType">
          <state code="Ningxia"/>
          ...
    """
    try:
        with open(xml_path, encoding='utf-8', errors='replace') as f:
            txt = f.read()
    except OSError:
        return []
    m = re.search(r'<generalDataType\s+id="location\.dataType"[^>]*>(.*?)'
                  r'</generalDataType>', txt, re.DOTALL)
    if not m:
        return []
    return [c for c in re.findall(r'<state\s+code="([^"]*)"', m.group(1)) if c]


def _resolve_locations(xml_path, meta_csv, logger):
    """位置列表解析: 优先读 XML 声明 (与 log/tree 注释同源), 兜底 metadata。"""
    locs = _locations_from_xml(xml_path)
    if locs:
        return locs
    locs = _meta_locations(meta_csv)
    if locs:
        logger.warning(f"  [phylogeo] 无法从 XML 读出 location 状态, 回退 metadata 原始名 "
                       f"({len(locs)} 个) —— SpreaD3 节点名可能与 MCC 树注释不一致")
    return locs


def run_stage_phylogeo(prep, results, aln_file, args, logger) -> bool:
    """phylogeo: 贝叶斯系统地理 (CTMC+BSSVS) + SpreaD3 可视化
    返回 True = 多链后台运行中 (pending, 不写 checkpoint, 下次重跑检查合并) """
    out_dir = prep["output_dir"]
    meta_csv = prep["meta_csv"]
    if not (aln_file and os.path.exists(aln_file) and meta_csv and os.path.exists(meta_csv)):
        logger.warning("  [phylogeo] 缺少 alignment 或 metadata, 跳过")
        return False

    tree_prior = args.phylogeo_prior
    if tree_prior == "auto":
        prior_cmp_dir = os.path.join(out_dir, ph_dir('geography'), "prior_comparison")
        logger.info(f"  [phylogeo] Auto prior: constant vs skyline...")
        try:
            auto_r = auto_select_tree_prior(
                fasta_file=aln_file, metadata_csv=meta_csv,
                output_dir=prior_cmp_dir, clock_model=args.phylogeo_clock,
                quick_steps=10, quick_chain=200_000,
                beast_bin=getattr(args, "beast_bin", None) or "beast",
                threads=min(args.threads, 8),
                beast_version=args.beast_version,
                log=LogCollector(logger))
        except NotImplementedError as e:
            # 历史坑: BEAST1 path sampling 未实现 (fail-fast), 降级用默认 skyline
            logger.warning(f"  [phylogeo] auto prior skipped — {e}; fallback skyline")
            auto_r = {"comparison_failed": True, "error": str(e),
                      "reason": f"auto prior 未实现: {e}; 回退默认 skyline (非比较结论)"}
        tree_prior = auto_r.get("selected", "skyline")
        results["metrics"]["tree_prior_selected"] = tree_prior
        results["metrics"]["tree_prior_log_bf"] = auto_r.get("log_bf")
        results["metrics"]["tree_prior_reason"] = auto_r.get("reason", "")
        # 2026-09-15 修复 (P1-9): 先验比较若未真正完成, 必须显式标记, 不能把回退的
        # skyline 当作"比较结论"混进汇总。旧行为下 BEAST1 的 auto 模式一步 MCMC
        # 都没跑, 却仍写 tree_prior_selected="skyline" + log_bf=None, 用户无法区分。
        if auto_r.get("comparison_failed"):
            results["metrics"]["tree_prior_comparison_failed"] = True
            results["metrics"]["tree_prior_reason"] = auto_r.get("reason", "")
            logger.warning(f"  [phylogeo] 先验自动比较未完成, 使用默认 skyline: "
                           f"{auto_r.get('error') or auto_r.get('reason')}")
        else:
            results["metrics"]["tree_prior_comparison_failed"] = False

    phylogeo_dir = os.path.join(out_dir, ph_dir('geography'), "phylogeography")
    xml_1 = os.path.join(phylogeo_dir, "phylogeo_beast1.xml")
    xml_2 = os.path.join(phylogeo_dir, "phylogeo.xml")
    geo_xml = None
    if not args.force and os.path.exists(xml_1):
        logger.info("  [phylogeo] XML exists, skip (BEAST1)")
        locs = _resolve_locations(xml_1, meta_csv, logger)
        geo_xml = {"xml_path": xml_1, "n_locations": len(locs), "locations": locs}
    elif not args.force and os.path.exists(xml_2):
        logger.info("  [phylogeo] XML exists, skip (BEAST2)")
        locs = _resolve_locations(xml_2, meta_csv, logger)
        geo_xml = {"xml_path": xml_2, "n_locations": len(locs), "locations": locs}
    else:
        logger.info(f"  [phylogeo] BEAST {args.beast_version}.x XML (prior={tree_prior})...")
        if str(args.beast_version) == "1":
            geo_xml = generate_beast1_phylogeo_xml(
                fasta_file=aln_file, metadata_csv=meta_csv,
                output_dir=phylogeo_dir,
                chain_length=args.phylogeo_chain,
                clock_model=args.phylogeo_clock,
                tree_prior=tree_prior,
                log=LogCollector(logger))
        else:
            geo_xml = generate_phylogeo_xml(
                fasta_file=aln_file, metadata_csv=meta_csv,
                output_dir=phylogeo_dir, chain_length=args.phylogeo_chain,
                clock_model=args.phylogeo_clock, tree_prior=tree_prior,
                log=LogCollector(logger))
            # 2026-09-15 修复 (P1-11): generate_phylogeo_xml (BEAST2) 现已 fail-loud
            # (该生成器用的是 BEAST1 语法, 产出与 BEAST2 不兼容)。原代码不检查
            # success, 会把 xml_path=None 一路带下去, 最后 submit 一个不存在的
            # XML —— 报错点离真正原因很远。这里显式失败并给出可执行的替代方案。
            if not geo_xml.get("success") or not geo_xml.get("xml_path"):
                results["error"] = geo_xml.get("error") or "BEAST2 phylogeo XML 生成失败"
                results["metrics"]["phylogeo_xml_failed"] = True
                logger.error(f"  [phylogeo] {results['error']}")
                logger.error("  [phylogeo] 建议改用 --beast_version 1 (BEAST1 路径)")
                return False
    results["outputs"]["phylogeo_xml"] = geo_xml.get("xml_path")
    results["metrics"]["phylogeo_locations"] = geo_xml.get("n_locations", 0)

    # ── 多链 BEAST 正式运行 (默认提交; --skip_phylogeo_run 跳过; 链 1 已存在则防重复) ──
    pending = False
    chain1_log = os.path.join(phylogeo_dir, "phylogeo_sky1.log")
    if not getattr(args, "skip_phylogeo_run", False):
        if os.path.exists(chain1_log) and not args.force:
            logger.info("  [phylogeo] 多链输出存在, 跳过提交")
        else:
            submit_phylogeo_chains(
                xml_path=geo_xml.get("xml_path") or os.path.join(phylogeo_dir, "phylogeo_beast1.xml"),
                work_dir=phylogeo_dir,
                chains=getattr(args, "phylogeo_chains", 3) or 3,
                threads=max(4, args.threads // max(getattr(args, "phylogeo_chains", 3) or 3, 1)),
                prefix="phylogeo_sky",
                chain_length=args.phylogeo_chain,
                logger=logger)

    # ── 多链闭环: 链全部完成后 logcombiner + treeannotator → MCC → TempMig ──
    try:
        from utils.merge_results import load_task, all_done, merge as merge_chains
        st = load_task(phylogeo_dir)
        if st and all_done(phylogeo_dir, st):
            merged_log = os.path.join(phylogeo_dir, "merged", "merged.log")
            mcc_tree = os.path.join(phylogeo_dir, "merged", "mcc.tree")
            if os.path.exists(merged_log) and os.path.exists(mcc_tree) and not args.force:
                logger.info("  [phylogeo] 合并产物已存在, 复用")
            else:
                logger.info("  [phylogeo] 全部链完成, 合并 (logcombiner + treeannotator)...")
                merge_chains(phylogeo_dir, st, background=False)
        elif st:
            # 前台等待闭环 (2026-09-02): 轮询到链完成再合并, 不再 pending 等二次重跑
            ok = _wait_chains_and_merge(
                phylogeo_dir, logger, label="phylogeo",
                merge_fn=lambda wd, s: merge_chains(wd, s, background=False))
            if not ok:
                pending = True
    except Exception as e:
        logger.warning(f"  [phylogeo] 多链合并步骤异常: {e}")

    # ── TempMig: 逐年迁移矩阵 (VirPhyKit TempMig 等效; 链合并后有 MCC 时) ──
    # 2026-09-16: 改用 _find_phylogeo_mcc (与 geo 阶段的 RRT 树注释法共用同一候选顺序)
    mcc_tree_path = _find_phylogeo_mcc(out_dir)
    if mcc_tree_path and meta_csv:
        tempmig_dir = os.path.join(out_dir, ph_dir('geography'), "tempmig")
        logger.info(f"  [phylogeo] TempMig migration matrix "
                    f"(MCC: {os.path.basename(mcc_tree_path)})...")
        run_postprocess("tempmig_full", logger,
                        "--mcc", mcc_tree_path,
                        "--dates", meta_csv,
                        "--outdir", tempmig_dir)

    # ── SpreaD3 可视化 + 发表级多面板图 (链完成 + merge 后) ──
    # 链产物名是 phylogeo_sky*.log/trees (prefix=phylogeo_sky), 不是 phylogeo_beast1*。
    # 2026-08-30 审计修复 (viz 门控): 旧版只 glob 链原始 log，链后台跑完 → 重跑本 stage
    # 时若 checkpoint 已 done 直接 skip，viz 块从未执行。加入 merged/merged.log 候选，
    # 且优先用合并链 (burnin 已处理、样本更完整)。
    log_files = list(Path(phylogeo_dir).glob("merged/merged.log")) \
        + list(Path(phylogeo_dir).glob("phylogeo_sky*.log")) \
        + list(Path(phylogeo_dir).glob("phylogeo_beast1*.log"))
    locations = geo_xml.get("locations", []) if geo_xml else []
    if log_files and locations:
        viz_dir = os.path.join(out_dir, ph_dir('geography'), "visualization")
        os.makedirs(viz_dir, exist_ok=True)
        logger.info(f"  [phylogeo] SpreaD3 map ({len(locations)} locations)...")
        spread3 = generate_spread3_json(
            # 2026-09-16 审查 P3 修复: 旧版恒传 ""，而上面第 1872 行明明已找到 MCC
            # 只喂给了 TempMig → 时空散点 (mcc_events) 整体缺失且只落一条 warning。
            log_file=str(log_files[0]), mcc_tree_file=(mcc_tree_path or ""),
            locations=locations, output_dir=viz_dir,
            bf_threshold=args.phylogeo_bf, log=LogCollector(logger))
        if spread3.get("success"):
            results["outputs"]["spread3_json"] = spread3["json_path"]
            results["outputs"]["spread3_map"] = generate_interactive_map_html(
                spread3["json_path"], f"{viz_dir}/phylogeo_map.html",
                title=f"{Path(out_dir).name} — Phylogeography", log=None)

        # 生成 migration_bf.csv (历史坑: generate_spread3_json 不写 bf_csv,
        # 只有 extract_migration_bf 写; 缺文件 → 多面板图 Panel B 永远 not available)
        try:
            bf_r = extract_migration_bf(
                log_file=str(log_files[0]), output_dir=viz_dir,
                locations=locations, bf_threshold=args.phylogeo_bf,
                log=LogCollector(logger))
            if bf_r.get("bf_table_csv"):
                logger.info(f"  [phylogeo] migration_bf.csv "
                            f"({bf_r.get('n_significant', 0)} 显著路线)")
        except Exception as e:
            logger.warning(f"  [phylogeo] migration_bf 提取失败: {e}")

        # 发表级系统地理多面板图 (trees 优先用合并后的 merged.trees)
        merged_trees = os.path.join(phylogeo_dir, "merged", "merged.trees")
        trees_files = ([merged_trees] if os.path.exists(merged_trees) else [])\
            + list(Path(phylogeo_dir).glob("phylogeo_sky*.trees"))
        if trees_files:
            loc_str = ",".join(locations)
            run_postprocess("phylogeo_figures", logger,
                            "--log", str(log_files[0]), "--trees", str(trees_files[0]),
                            "--locations", loc_str,
                            "--migration", os.path.join(viz_dir, "migration_bf.csv"),
                            "--output", os.path.join(viz_dir, "phylogeo_multipanel.pdf"))

        # T5 (2026-09-16): skyline 95% HPD 带 → 单文件自包含 HTML (日历年)。
        # 报告里只放相对路径链接 (见 report_builder PLOT_LINK_EXTS), 不内联数 MB 的 plotly。
        try:
            from phylogeo_figures import build_skyline_html, max_tip_year_from_metadata
            _mty = None
            _yr = (geo_xml or {}).get("year_range")
            if _yr:
                _mty = float(_yr[1])
            if _mty is None:
                _mty = max_tip_year_from_metadata(meta_csv, log=LogCollector(logger))
            sky_html = build_skyline_html(
                str(log_files[0]), os.path.join(viz_dir, "phylogeo_skyline.html"),
                max_tip_year=_mty, log=LogCollector(logger))
            if sky_html:
                results["outputs"]["skyline_html"] = sky_html
                logger.info(f"  [phylogeo] skyline HTML: {os.path.basename(sky_html)} "
                            f"(日历年{'可用' if _mty else '未知, 用距今年代'})")
        except Exception as e:
            logger.warning(f"  [phylogeo] skyline HTML 生成失败: {type(e).__name__}: {e}")
    return pending


# ═══════════════════════════════════════════════════════════════════
# Phase 5: BEAST (full-length)
# ═══════════════════════════════════════════════════════════════════

def run_phase_beast(prep, phase_results, args, tools, logger) -> Dict:
    vname = prep["virus_name"]
    out_dir = prep["output_dir"]
    log = logger

    beast_dir = os.path.join(out_dir, ph_dir('time'), "beast")
    beast_xml = os.path.join(beast_dir, "beast1.xml")
    beast_log = os.path.join(beast_dir, "beast1.log")

    results = {"virus": vname, "phase": "beast", "success": False,
               "outputs": {}, "metrics": {}}

    aln_file = phase_results.get("outputs", {}).get("alignment")
    dates_csv = prep.get("dates_csv", "")
    if not aln_file or not os.path.exists(aln_file):
        results["error"] = "No alignment file, skip BEAST"
        log.warning(f"  [beast] {results['error']}")
        return results
    if not dates_csv or not os.path.exists(dates_csv):
        results["error"] = "No dates file, skip BEAST"
        log.warning(f"  [beast] {results['error']}")
        return results

    log.info(f"\n[Stage beast] {vname} — 全长 BEAST 分子定年")

    if not args.force and os.path.exists(beast_xml):
        log.info(f"  BEAST1 XML exists, skip")
        results["outputs"]["beast_xml"] = beast_xml
    else:
        log.info(f"  Generating BEAST1 XML (UCLN + skyline, clock-only dating)...")
        xml_r = generate_beast1_phylogeo_xml(
            fasta_file=aln_file, metadata_csv=dates_csv, output_dir=beast_dir,
            chain_length=args.beast_chain, log_every=max(1000, args.beast_chain // 100),
            clock_model="ucln", tree_prior="skyline",
            discretize_locations=False,   # 全长分子定年: 纯分子钟 (无 location CTMC/BSSVS)
            substitution_model="auto", gamma_categories=4,
            xml_name="beast1.xml", log=LogCollector(log))
        if xml_r["success"]:
            results["outputs"]["beast_xml"] = xml_r["xml_path"]
            results["metrics"]["beast_n_taxa"] = xml_r["n_taxa"]
            if xml_r.get("year_range"):
                y0, y1 = xml_r["year_range"]
                results["metrics"]["year_range"] = f"{y0:.1f}–{y1:.1f}"
        else:
            log.error(f"  BEAST XML failed: {xml_r.get('error')}")
            return results

    beauti_dir = os.path.join(beast_dir, "beauti")
    beauti_r = prepare_beauti_inputs(aln_file, dates_csv, beauti_dir, log=LogCollector(log))
    if beauti_r["success"]:
        results["outputs"]["beauti_fasta"] = beauti_r["fasta_path"]
        results["outputs"]["beauti_traits"] = beauti_r["traits_path"]
        results["metrics"]["beauti_locations"] = beauti_r.get("locations", [])

    if args.skip_beast_run:
        log.info(f"  BEAST run skipped (--skip_beast_run)")
    elif not os.path.exists(beast_xml):
        log.warning(f"  No XML to run")
        return results
    else:
        beast_chains = getattr(args, "beast_chains", 3) or 1
        # 多链模式下链 1 的产物名 (token 替换: beast1.log → beast1{1}.log = beast11.log)
        chain1_log = os.path.join(beast_dir, "beast11.log")
        if beast_chains <= 1:
            # 单链: 前台等待 (原逻辑)
            if os.path.exists(beast_log) and not args.force:
                log.info(f"  BEAST output exists, skip run")
                results["outputs"]["beast_log"] = beast_log
                results["outputs"]["beast_trees"] = os.path.join(beast_dir, "beast1.trees")
            else:
                beast_bin = args.beast_bin or "beast"
                if not beast_bin:
                    log.warning(f"  No BEAST binary configured, skip run")
                else:
                    log.info(f"  Running BEAST1 (chain={args.beast_chain:,}, threads={args.threads})...")
                    beast_r = run_beast(
                        xml_path=beast_xml, output_dir=beast_dir, beast_bin=beast_bin,
                        threads=args.threads, timeout_hours=args.beast_timeout,
                        log=LogCollector(log))
                    if beast_r["success"]:
                        results["outputs"]["beast_log"] = beast_r.get("log_path")
                        results["outputs"]["beast_trees"] = beast_r.get("trees_path")
                        results["metrics"]["beast_elapsed"] = f"{beast_r['elapsed_seconds']:.0f}s"
                    else:
                        log.warning(f"  BEAST run failed: {beast_r.get('error')}")
        else:
            # 多链: 默认提交 N 链 (后台, 复用 submit_mcmc_chains)
            if os.path.exists(chain1_log) and not args.force:
                log.info(f"  BEAST 多链输出存在, 跳过提交")
            else:
                log.info(f"  [beast] 提交 {beast_chains} 链 × {max(4, args.threads // beast_chains)} 线程 (多链, 默认)...")
                submit_mcmc_chains(
                    xml_path=beast_xml, work_dir=beast_dir,
                    chains=beast_chains,
                    threads=max(4, args.threads // beast_chains),
                    prefix="beast1", chain_length=args.beast_chain,
                    logger=log)
                results["metrics"]["beast_chains_submitted"] = beast_chains

            # 多链闭环: 全部链完成 → merge_results 合并 (logcombiner + treeannotator)
            merged_log = os.path.join(beast_dir, "merged", "merged.log")
            try:
                from utils.merge_results import load_task, all_done, merge as merge_chains
                st = load_task(beast_dir)
                if st and all_done(beast_dir, st):
                    if os.path.exists(merged_log) and not args.force:
                        log.info(f"  [beast] 合并产物已存在, 复用 merged.log")
                    else:
                        log.info(f"  [beast] 全部链完成, 合并 (logcombiner + treeannotator)...")
                        merge_chains(beast_dir, st, background=True)
                    if os.path.exists(merged_log):
                        results["outputs"]["beast_log"] = merged_log
                        results["outputs"]["beast_trees"] = os.path.join(beast_dir, "merged", "mcc.tree")
                    else:
                        log.warning(f"  [beast] merged.log 未生成 (logcombiner 失败?)")
                else:
                    # 前台等待闭环 (用户拍板 2026-09-02 不放后台): 轮询到链完成再合并
                    if st:
                        ok = _wait_chains_and_merge(
                            beast_dir, log, label="beast",
                            merge_fn=lambda wd, s: merge_chains(wd, s, background=True))
                        if ok and os.path.exists(merged_log):
                            results["outputs"]["beast_log"] = merged_log
                            results["outputs"]["beast_trees"] = os.path.join(beast_dir, "merged", "mcc.tree")
                        elif not ok:
                            results["beast_pending"] = True   # 等待失败退回 pending 语义
                            if os.path.exists(chain1_log):
                                results["outputs"]["beast_log"] = chain1_log
                                results["outputs"]["beast_trees"] = os.path.join(beast_dir, "beast11.trees")
            except Exception as e:
                log.warning(f"  [beast] 多链合并步骤异常: {e}")
                if os.path.exists(chain1_log):
                    results["outputs"]["beast_log"] = chain1_log
                    results["outputs"]["beast_trees"] = os.path.join(beast_dir, "beast11.trees")

    pp_dir = os.path.join(beast_dir, "postprocess")
    log_file = results["outputs"].get("beast_log")
    trees_file = results["outputs"].get("beast_trees")
    if log_file or trees_file:
        log.info(f"  Post-processing BEAST output...")
        rscript = tools.get("Rscript") or "Rscript"
        post = postprocess_beast(
            log_file=log_file or "", trees_file=trees_file or "",
            output_dir=pp_dir, burnin_pct=args.beast_burnin,
            rscript_path=rscript, log=LogCollector(log))
        if post.get("bsp_plot"):
            results["outputs"]["bsp_plot"] = post["bsp_plot"]
        if post.get("rspp_plot"):
            results["outputs"]["rspp_plot"] = post["rspp_plot"]
        if post.get("rspp_states"):
            results["metrics"]["rspp_root_states"] = post["rspp_states"]
        if post.get("ess"):
            ess = post["ess"]
            results["metrics"]["ess_n_params"] = ess.get("n_params", 0)
            results["metrics"]["ess_passed"] = ess.get("n_passed", 0)
            results["metrics"]["ess_warned"] = ess.get("n_warned", 0)
            results["outputs"]["ess_table"] = ess.get("ess_table")
            if ess.get("n_warned", 0) > 0:
                log.warning(f"  ⚠ BEAST ESS: {ess['n_warned']} parameters below threshold!")

    # ── 全长 TMRCA: parse_beast_log (BEAST 1.x treeModel.rootHeight) → 喂 gene_dating ──
    if log_file and os.path.exists(log_file):
        try:
            parsed = parse_beast_log(log_file, burnin_pct=args.beast_burnin, log=LogCollector(log))
            kp = parsed.get("key_params", {})
            rh = kp.get("treeModel.rootHeight")
            if rh:
                results["metrics"]["beast_tmrca"] = rh["median"]
                results["metrics"]["beast_tmrca_hpd"] = f"{rh['hpd_95_lower']:.1f}–{rh['hpd_95_upper']:.1f}"
                results["metrics"]["beast_tmrca_ess"] = rh["ess"]
                log.info(f"  [beast] 全长 TMRCA: {rh['median']:.1f} yr "
                         f"[{rh['hpd_95_lower']:.1f}–{rh['hpd_95_upper']:.1f}] ESS={rh['ess']:.0f}")
            for k, v in [("ucld.mean", "beast_ucld_mean"), ("ucld.stdev", "beast_ucld_stdev")]:
                if k in kp:
                    results["metrics"][v] = kp[k]["mean"]
        except Exception as e:
            log.warning(f"  [beast] TMRCA 解析失败: {e}")

    results["success"] = not results.get("beast_pending", False)
    return results


# ═══════════════════════════════════════════════════════════════════
# Per-virus driver
# ═══════════════════════════════════════════════════════════════════

def run_postprocess(script_name: str, logger, *extra_args) -> None:
    """调用后处理脚本 (-m, 可带额外参数); 失败仅告警, 不阻断主流程。"""
    if not script_name:
        return
    try:
        logger.info(f"  [postprocess] {script_name} {' '.join(extra_args[:4])}...")
        # utils/ 下的后处理模块用 -m 调用时需加 utils. 前缀
        # (否则 ModuleNotFoundError, 如 sample_plot / tempmig_full)
        mod = script_name
        if os.path.exists(os.path.join(SCRIPT_DIR, "utils", f"{script_name}.py")):
            mod = f"utils.{script_name}"
        r = subprocess.run(
            [sys.executable, "-m", mod] + list(extra_args),
            capture_output=True, text=True, encoding='utf-8', errors='replace', cwd=str(SCRIPT_DIR), timeout=3600)
        if r.returncode != 0:
            tail = (r.stderr or r.stdout).strip().splitlines()
            logger.warning(f"  [postprocess] {script_name} 失败: {tail[-1][:200] if tail else '?'}")
        else:
            logger.info(f"  [postprocess] {script_name} 完成")
    except Exception as e:
        logger.warning(f"  [postprocess] {script_name} 异常: {e}")


# ═══════════════════════════════════════════════════════════════════
# metadata 三通道拆分 (2026-09-16)
# ═══════════════════════════════════════════════════════════════════
# 背景: 要求「在元数据清洗检查后, 定年和地理分开去研究」。
#   · 计算层早已分开 —— `time` 组 (clock/beast/gene_dating) 与
#     `geography` 组 (phylogeo/geo) 在 STAGE_DEPS 里都只依赖 tree, 互不依赖;
#     beast 侧 `discretize_locations=False` (纯分子钟), geo 侧 Mantel/RRT
#     根本不需要 BEAST。
#   · 但**数据层**此前是断的: 通道表 (dates.csv / sample_metadata.csv) 由
#     `prepare_virus_inputs()` 在 prep 阶段写出, 而 metadata 治理发生在 prep
#     **之后** → 治理结果只回写到 `args.metadata`, 而它的唯一消费点
#     (`load_sample_metadata`) 早就跑完了 → 治理形同只影响日志。
# 本处把治理产物**按样品名回填**到通道表, 使各条分析线各吃各的表。
#
# 2026-09-16 第三轮 (老师: 「host 也独立出来, 仅对有 host 的进行分析」):
# 通道由两张扩为**三张** —— `host` 从地理表的附属列独立成 `data/host.csv`,
# 且**只收有宿主的样品** (判据 `has_host()`), 于是宿主分析天然只跑在有宿主的子集上。
#   dates.csv            name, date[, location]  → 时间通道
#   sample_metadata.csv  name, date, location    → 地理通道
#   host.csv             name, host, tissue      → 宿主通道 (仅含宿主样品)
def _channel_paths(out_dir: str) -> Dict[str, str]:
    """三张通道表 + 留痕报告的约定位置 (与 prep 的产物同目录)。"""
    d = os.path.join(out_dir, ph_dir('data'))
    return {"dates": os.path.join(d, "dates.csv"),
            "meta": os.path.join(d, "sample_metadata.csv"),
            "host": os.path.join(d, "host.csv"),
            "report": os.path.join(d, "metadata", "channels.json")}


def _refresh_metadata_channels(prep, args, logger) -> Dict:
    """把 metadata (治理产物优先) 回填到 时间 / 地理 / 宿主 三张通道表。

    三种入口统一口径:
      · batch     —— prep 已按样品写出三张表 → **就地回填** (覆盖前自动备份)
      · --virus   —— 骨架 = datasets.yaml 指的那份 metadata → 派生到 <work>/data/
      · work_dir  —— 骨架 = --meta / --host_meta → 派生到 <work>/data/

    来源优先级: 治理产物 (prep['metadata_std'] 或 <out>/data/metadata/metadata_std.csv)
                > 当前通道表自身 (无治理时只做"拆分", 不改值)。

    幂等: 来源指纹未变且产物在 → 跳过 (checkpoint 语义)。
    失败**不阻断**流程: 告警并沿用原表 (与 metadata/clean 同款 fail-soft)。
    """
    from utils import metadata_channels as _mc

    out_dir = prep.get("output_dir")
    if not out_dir:
        return {}
    tgt = _channel_paths(out_dir)
    cur_dates = prep.get("dates_csv") or ""
    cur_meta = prep.get("meta_csv") or ""
    cur_host = prep.get("host_csv") or ""

    gov = prep.get("metadata_std") or os.path.join(
        out_dir, ph_dir('data'), "metadata", "metadata_std.csv")
    gov = gov if os.path.exists(gov) else None

    if gov:
        src_date = src_geo = gov
    else:
        src_date = cur_dates if os.path.exists(cur_dates) else ""
        src_geo = cur_meta if os.path.exists(cur_meta) else src_date
        if not src_date:
            return {}
        if os.path.abspath(cur_dates) == os.path.abspath(tgt["dates"]):
            # batch 且无治理产物: 三张表本就在约定位置, 无需拆分
            logger.info("  [channels] 无治理产物且通道表已就位, 跳过拆分")
            return {}

    if not _mc.needs_rebuild(tgt["report"], src_date, src_geo):
        prep["dates_csv"] = tgt["dates"]
        if os.path.exists(tgt["meta"]):
            prep["meta_csv"] = tgt["meta"]
        if os.path.exists(tgt["host"]):
            prep["host_csv"] = tgt["host"]
        logger.info("  [channels] checkpoint OK, skip (来源未变)")
        return {"skipped": True, "dates_csv": prep["dates_csv"],
                "meta_csv": prep.get("meta_csv"),
                "host_csv": prep.get("host_csv")}

    r = _mc.rebuild_channels(
        cur_dates, cur_meta or None, host_csv=cur_host or None,
        meta_source=src_date, geo_source=src_geo,
        out_dates=tgt["dates"], out_meta=tgt["meta"], out_host=tgt["host"],
        logger=logger)
    if not r.get("success"):
        logger.warning(f"  [channels] 拆分失败 ({r.get('error')}), 沿用原表")
        return r
    # 关键: 把下游 (beast/phylogeo/geo/host/popgen/clean) 指到新表
    prep["dates_csv"] = r["dates_csv"]
    if r.get("meta_csv") and os.path.exists(r["meta_csv"]):
        prep["meta_csv"] = r["meta_csv"]
    # 宿主通道: 落地了才改指向 (未落地 → 保留原值, 不清空)
    if r.get("host_csv") and os.path.exists(r["host_csv"]):
        prep["host_csv"] = r["host_csv"]
    prep["_channels_source"] = src_date
    return r


def _governed_meta(prep, fallback_dates: str = "",
                   fallback_meta: str = "") -> Dict[str, str]:
    """收尾后处理 (publication_figures / sample_plot) 该吃的 metadata 路径。

    2026-09-16 与 stage 侧同源 —— **治理产物优先**。

    此前后处理 (单病毒 --virus / work_dir 两个分支) 直接用 `main()` 里在治理**之前**
    捕获的 `metadata` 局部变量, 或 work_dir 的 `_resolve(args.meta)` → 于是
    "治理只影响日志": TempEst 的日期与 sample_plot 的地理仍是脏值。
    (batch 分支早已吃 `vq['meta_csv']`, 因 `vq` 就是 prep 同一对象,
     `process_virus` 就地改写后自然指向治理回填表 —— 属"巧合正确", 一并在此显式化。)

    三张通道表在 `_refresh_metadata_channels()` 之后已指向治理回填后的文件,
    直接复用; 不存在时逐级回退到传入的原表, 保证不会因缺文件而整体漏跑。
    """
    d = prep.get("dates_csv") or ""
    g = prep.get("meta_csv") or ""
    if not (d and os.path.exists(d)) and fallback_dates and os.path.exists(fallback_dates):
        d = fallback_dates
    if not (g and os.path.exists(g)):
        if fallback_meta and os.path.exists(fallback_meta):
            g = fallback_meta
        else:
            g = g or d
    # 宿主通道: **不回退**到地理表 —— 曾靠"地理表里凑巧有 host 列"运转,
    # 现在 host 已独立; 没有 host.csv 就是"没有宿主信息", 由消费方按缺省处理。
    h = prep.get("host_csv") or ""
    if not (h and os.path.exists(h)):
        h = ""
    return {"dates": d, "meta": g, "host": h}


def process_virus(prep, stages: List[str], args, tools, logger) -> Dict:
    """单个病毒按 stage 清单顺序执行, 返回汇总 results (含 stage 状态追踪)"""
    vname = prep["virus_name"]
    out_dir = prep["output_dir"]
    log = logger
    st = {"done": [], "skipped": [], "failed": []}
    results = {"virus": vname, "success": False, "outputs": {}, "metrics": {},
               "stages_done": [], "stages_skipped": [], "stages_failed": [],
               "virus_dir": out_dir}  # virus_dir: report_builder 按模块扫产物用 (2026-09-02)

    log.info(f"\n{'='*60}")
    log.info(f"Processing {vname} | stages: {','.join(stages)}")
    log.info(f"{'='*60}")

    force = getattr(args, 'force', False)
    aln_file = None
    tree_file = None

    def _inputs_of(stage):
        """该 stage 的输入指纹 (仅对"输入明确且稳定"的 stage 启用; 2026-09-15 P1-1)。

        align 的输入 = 上游合并 fasta; tree 的输入 = 比对文件 + **建树参数**。
        其余 stage 输入来源复杂 (BEAST XML / 外部工具配置等), 暂不启用,
        避免误判导致不必要的全量重跑。
        """
        try:
            if stage == 'align':
                return _inputs_fingerprint([prep.get('combined_fasta')])
            if stage == 'tree':
                # 2026-09-16: 并入 iqtree 参数指纹 —— 加自举后参数一变就重建,
                # 否则旧的无支持值主树会被 checkpoint 静默沿用。
                return _inputs_fingerprint([aln_file or prep.get('combined_fasta')],
                                           params=iqtree_param_token(args))
        except Exception:
            return None
        return None

    def _done(stage):
        return stage_done(out_dir, stage, force, inputs=_inputs_of(stage))

    def _mark(stage, meta=None):
        fp = _inputs_of(stage)
        if fp is not None:
            meta = dict(meta or {})
            meta.setdefault('inputs', fp)
        mark_stage_done(out_dir, stage, meta)
        st["done"].append(stage)

    def _skip(stage):
        st["skipped"].append(stage)

    def _fail(stage, msg=""):
        st["failed"].append(stage)
        if msg:
            log.warning(f"  [{stage}] {msg}")

    # ── metadata: 入口治理 (时间 + 地理 的检查与矫正) ──
    # 2026-09-16 新增 (老师拍板)。必须在 online 之前: 先把已有 metadata 归一化,
    # 下游 (clean 的 drop_no_date / beast 取小数年 / phylogeo 取坐标) 才不会被
    # 占位符日期、Unknown 地理层级、越界日期静默污染。
    if 'metadata' in stages:
        if _done('metadata'):
            log.info("  [metadata] checkpoint OK, skip"); _skip('metadata')
            sc = os.path.join(out_dir, ph_dir('data'), "metadata", "metadata_std.csv")
            if os.path.exists(sc):
                prep["metadata_std"] = sc
                prep["_metadata_original"] = getattr(args, "metadata", None)
                args.metadata = sc          # 下游一律消费标准化产物
                results["outputs"]["metadata_std"] = sc
        else:
            mr = run_stage_metadata(prep, args, log)
            results["outputs"].update(mr.get("outputs", {}))
            results["metrics"]["metadata"] = mr.get("metrics", {})
            if mr.get("success") and mr.get("std_csv"):
                # 关键: 把下游 metadata 输入换成标准化产物
                prep["metadata_std"] = mr["std_csv"]
                prep["_metadata_original"] = getattr(args, "metadata", None)
                args.metadata = mr["std_csv"]
                _mark('metadata')
            else:
                # 治理失败不阻断流程, 显式告警 (降级为直接用原 metadata)
                log.warning(f"  [metadata] 失败, 降级使用原 metadata: {mr.get('error')}")
                _skip('metadata')

    # ── metadata 三通道拆分 (2026-09-16): 治理产物 → 时间表 / 地理表 / 宿主表 ──
    # 位置很关键: 必须在 clean 之前 (clean 的 --clean_drop_no_date/-location 要吃它),
    # 也必须在 clock/beast (吃 dates_csv) 与 phylogeo/geo (吃 meta_csv)、
    # host (吃 host_csv) 之前。
    # 无治理产物时本调用是廉价的 no-op 或"纯拆分"; 失败不阻断。
    _ch = _refresh_metadata_channels(prep, args, log)
    if _ch.get("dates_csv"):
        results["outputs"]["dates_csv"] = _ch["dates_csv"]
    if _ch.get("meta_csv"):
        results["outputs"]["sample_metadata"] = _ch["meta_csv"]
    if _ch.get("host_csv"):
        results["outputs"]["host_metadata"] = _ch["host_csv"]
    if _ch.get("success") and not _ch.get("skipped"):
        results["metrics"]["metadata_channels"] = _ch.get("stats", {})

    # ── online: NCBI 参考补充 ──
    if 'online' in stages:
        if _done('online'):
            log.info("  [online] checkpoint OK, skip"); _skip('online')
        else:
            r = run_stage_online(prep, args, log)
            if r.get("ncbi_fasta"):
                results["outputs"]["ncbi_fasta"] = r["ncbi_fasta"]
                _mark('online')
            else:
                # 在线参考补充是可选项 (work_dir 模式无 accession 时必然无产物)
                _skip('online')

    # ── clean: 第一道过滤 (比对前, 序列本体清洗) ──
    # 2026-09-16 新增。必须在 align 之前: 空格/星号/小写会让 MAFFT 出错或错位比对,
    # 事后 QC (align_qc) 只能发现问题, 不能修复比对本身。
    if 'clean' in stages:
        if _done('clean'):
            log.info("  [clean] checkpoint OK, skip"); _skip('clean')
            cf = os.path.join(out_dir, ph_dir('data'), "clean", "clean.fasta")
            if os.path.exists(cf):
                prep["combined_fasta"] = cf
                prep["_raw_combined_fasta"] = prep.get("combined_fasta")
                results["outputs"]["clean_fasta"] = cf
        else:
            cr = run_stage_clean(prep, args, log)
            results["outputs"].update(cr.get("outputs", {}))
            results["metrics"]["clean"] = cr.get("metrics", {})
            if cr.get("success") and cr.get("clean_fasta"):
                # 关键: 把比对输入换成清洗后的序列
                prep["_raw_combined_fasta"] = prep.get("combined_fasta")
                prep["combined_fasta"] = cr["clean_fasta"]
                _mark('clean')
            else:
                # 清洗失败不阻断流程, 但要显式告警 (降级为直接用原序列)
                log.warning(f"  [clean] 失败, 降级使用原始序列: {cr.get('error')}")
                _skip('clean')

    # ── align: MAFFT 比对 + 内联 saturation (纯Python秒级, 合并运行) ──
    if 'align' in stages:
        if _done('align'):
            log.info("  [align] checkpoint OK, skip"); _skip('align')
            aln_file = os.path.join(out_dir, ph_dir('phylogeny'), "mafft.aln.fasta")
            if not os.path.exists(aln_file):
                aln_file = prep.get("combined_fasta")
            results["outputs"]["alignment"] = aln_file
            # saturation 随 align 产出; 老数据补跑
            if not _done('saturation') and aln_file and os.path.exists(aln_file):
                sr = run_stage_saturation(prep, aln_file, args, log)
                results["outputs"].update(sr.get("outputs", {}))
                results["metrics"].update(sr.get("metrics", {}))
                if sr.get("success"):
                    _mark('saturation')
        else:
            r = run_stage_align(prep, tools, args, log, force=force)
            if r.get("success") and r.get("alignment"):
                aln_file = r["alignment"]
                results["outputs"]["alignment"] = aln_file
                _mark('align')
                # 内联 saturation: 比对一旦产出立即验饱和 (闸门逻辑)
                sr = run_stage_saturation(prep, aln_file, args, log)
                results["outputs"].update(sr.get("outputs", {}))
                results["metrics"].update(sr.get("metrics", {}))
                if sr.get("success"):
                    _mark('saturation')
                # 后处理: 比对 QC 报告 (长度/gap%/N%/identity, 不覆盖原比对)
                if aln_file and os.path.exists(aln_file):
                    run_postprocess("align_qc", log,
                                    "--input", aln_file,
                                    "--outdir", os.path.join(out_dir, ph_dir('phylogeny'), "qc"),
                                    "--skip-mafft")
                    # 比对 QC 图 (长度分布 + gap 比例, SCI 风格)
                    try:
                        from utils.stage_plots import plot_align_qc
                        _p = plot_align_qc(aln_file, os.path.join(out_dir, ph_dir('phylogeny'), "qc"))
                        if _p:
                            results["outputs"]["align_qc_plot"] = _p
                    except Exception as _e:
                        log.warning(f"  [align_qc_plot] skip: {_e}")
            else:
                _fail('align', r.get("error", "failed"))

    # ── splitstree: SplitsTree6 分裂网络 (树状/网状信号) ──
    if 'splitstree' in stages and aln_file:
        if _done('splitstree'):
            log.info("  [splitstree] checkpoint OK, skip"); _skip('splitstree')
        else:
            ssr = run_stage_splitstree(prep, aln_file, args, log)
            results["outputs"].update(ssr.get("outputs", {}))
            results["metrics"].update(ssr.get("metrics", {}))
            if ssr.get("success"):
                _mark('splitstree')
            else:
                _fail('splitstree', '分裂网络未产出')

    # ── rdp5: 重组检测 → 删重组子+重比对 (默认) / mask (可选) ──
    if 'rdp5' in stages and aln_file:
        if _done('rdp5'):
            log.info("  [rdp5] checkpoint OK, skip"); _skip('rdp5')
            # 恢复比对: 优先新策略产物, 回退旧 masked (老数据兼容)
            r_aln = os.path.join(out_dir, ph_dir('recomb'), "rdp5", "rdp5_dropped.realn.fasta")
            if not os.path.exists(r_aln):
                r_aln = os.path.join(out_dir, ph_dir('recomb'), "rdp5", "masked", "masked_N.fasta")
            if os.path.exists(r_aln):
                aln_file = r_aln
                results["outputs"]["masked_alignment"] = r_aln
        else:
            masked = run_stage_rdp5(aln_file, prep, tools, args, log)
            if masked:
                aln_file = masked
                results["outputs"]["masked_alignment"] = masked
                _mark('rdp5')
                # 后处理: 重组事件验证层 (validate_recombinants, 生成验证任务 + 汇总表)
                rdp5_csv = os.path.join(out_dir, ph_dir('recomb'), "rdp5", f"{prep['virus_name']}.csv")
                if os.path.exists(rdp5_csv):
                    run_postprocess("validate_recombinants", log,
                                    "--events", rdp5_csv,
                                    "--alignment", aln_file,
                                    "--outdir", os.path.join(out_dir, ph_dir('recomb'), "rdp5", "validation"))
            else:
                _fail('rdp5', 'RDP5/mask 未产出干净比对')

    # ── tree: IQ-TREE 建树 ──
    if 'tree' in stages and aln_file:
        if _done('tree'):
            log.info("  [tree] checkpoint OK, skip"); _skip('tree')
            tree_file = os.path.join(out_dir, ph_dir('phylogeny'), "iqtree.treefile")
            if not os.path.exists(tree_file):
                tree_file = prep.get("tree_file")
            results["outputs"]["tree"] = tree_file
        else:
            r = run_stage_tree(prep, aln_file, tools, args, log, force=force)
            if r.get("success") and r.get("tree"):
                tree_file = r["tree"]
                results["outputs"]["tree"] = tree_file
                _mark('tree')
                # 后处理: ML 树环形可视化 (纯 matplotlib, 无 ETE 依赖)
                try:
                    from utils.stage_plots import plot_tree_circular
                    _tp = plot_tree_circular(tree_file, os.path.join(out_dir, ph_dir('phylogeny')))
                    if _tp:
                        results["outputs"]["tree_plot"] = _tp
                except Exception as _e:
                    log.warning(f"  [tree_plot] skip: {_e}")
                # 后处理: 交叉验证 (RAxML-NG 与 IQ-TREE 对照 + 拓扑比较, 工具缺失优雅跳过)
                if aln_file and os.path.exists(aln_file) and tree_file and os.path.exists(tree_file):
                    run_postprocess("phylo_cross", log,
                                    "raxml", "--aln", aln_file,
                                    "--outdir", os.path.join(out_dir, ph_dir('phylogeny'), "cross_validation"),
                                    "--ref-tree", tree_file)
            else:
                _fail('tree', r.get("error", "failed"))

    # concat stage 已移除 (2026-08-27: 无下游消费者, 研究不涉拓扑叙事;
    # 工具函数 concatenate_alignments 保留在 multigene_tools)

    # ── popgen: 群体遗传学 (pypopart) ──
    if 'popgen' in stages:
        if _done('popgen'):
            log.info("  [popgen] checkpoint OK, skip"); _skip('popgen')
        else:
            pr = run_stage_popgen(prep, aln_file, args, log)
            results["outputs"].update(pr.get("outputs", {}))
            results["metrics"].update(pr.get("metrics", {}))
            if pr.get("success"):
                _mark('popgen')
                # 后处理: popgen 指标柱状图 (π/θ/Tajima's D/Fst, SCI 风格)
                try:
                    from utils.stage_plots import plot_popgen
                    _pp = os.path.join(out_dir, ph_dir('popgen'), "popgen_report.txt")
                    if os.path.exists(_pp):
                        _pg = plot_popgen(_pp, os.path.join(out_dir, ph_dir('popgen')))
                        if _pg:
                            results["outputs"]["popgen_plot"] = _pg
                except Exception as _e:
                    log.warning(f"  [popgen_plot] skip: {_e}")
                # 后处理: 遗传距离矩阵 + NJ 树 + 热图 (popgen 多样性可视化)
                if aln_file and os.path.exists(aln_file):
                    run_postprocess("distance_matrix", log,
                                    "--fasta", aln_file,
                                    "--outdir", os.path.join(out_dir, ph_dir('popgen'), "distance"))
                # 后处理: 单倍型网络图 + 变异热点图 (datasets 配置驱动, 无配置跳过)
                vkey = prep.get("dataset_key", "")
                if vkey:
                    # 历史坑: prep 字典从不含 haplo/popgen_out/genes_csv (单病毒手工构造、
                    # 批量 prepare_virus_inputs 均无) → 钩子永不触发, 配置静默失效。
                    # 改为从 datasets.yaml 直接读。
                    try:
                        from utils.dataset_config import dataset as _ds_get
                        _dcfg = _ds_get(vkey) or {}
                    except Exception:
                        _dcfg = {}
                    if prep.get("haplo") or _dcfg.get("haplo"):
                        run_postprocess("build_haplo_outputs", log, "--virus", vkey)
                    if prep.get("popgen_out") or prep.get("genes_csv") \
                            or _dcfg.get("popgen_out") or _dcfg.get("genes_csv"):
                        run_postprocess("popgen_figures", log, "--virus", vkey)
                # 后处理: DnaSP 补充检验 (Fu&Li D*/F* / Rm / Ka-Ks, 需 pop 文件)
                meta_csv = prep.get("meta_csv") or prep.get("dates_csv", "")
                if meta_csv and os.path.exists(meta_csv):
                    pop_file = os.path.join(out_dir, ph_dir('popgen'), "location_pops.tsv")
                    try:
                        import csv as _csv
                        rows = []
                        with open(meta_csv, encoding='utf-8-sig') as _f:
                            for _row in _csv.DictReader(_f):
                                _n = _row.get('name', '').strip()
                                _l = _row.get('location', '').strip()
                                if _n and _l:
                                    rows.append(f"{_n}\t{_l}")
                        if rows:
                            with open(pop_file, 'w', encoding='utf-8') as _pf:
                                _pf.write("#seq\tpop\n" + "\n".join(rows))
                            # 2026-08-30 审计修复 B3: 全基因组比对当单 CDS 算 Ka/Ks 是方法学错误
                            # (基因间区混入 → ω 无意义; 长度非 3 倍数时直接跳过)。改为
                            # 默认分析集去掉 kaks，Ka/Ks 走 capheine 分基因路径 (codeml/cawlign)。
                            run_postprocess("dnasp_bridge", log,
                                            "--fasta", aln_file, "--pop-file", pop_file,
                                            "--analyses", "polymorphism,fufs,recombination,ld",
                                            "--out", os.path.join(out_dir, ph_dir('popgen'), "dnasp"),
                                            "--cli-report")
                    except Exception as _e:
                        log.warning(f"  [popgen] dnasp 补充检验跳过: {_e}")
            else:
                _fail('popgen', 'pypopart 未产出报告 (可能缺 pypopart/元数据)')

    # ── host: 宿主分化（数据不适用时视为 skipped）──
    if 'host' in stages:
        if _done('host'):
            log.info("  [host] checkpoint OK, skip"); _skip('host')
        else:
            run_stage_host(prep, results, aln_file, tree_file, args, log)
            if results.get("outputs", {}).get("host_tree") or results.get("outputs", {}).get("host_diff_plot"):
                _mark('host')
            else:
                _skip('host')

    # ── capheine: 正选择分析 (Phase S) ──
    if 'capheine' in stages:
        if _done('capheine'):
            log.info("  [capheine] checkpoint OK, skip"); _skip('capheine')
        else:
            cr = run_stage_capheine(prep, args, log)
            results["outputs"].update(cr.get("outputs", {}))
            results["metrics"].update(cr.get("metrics", {}))
            if cr.get("success"):
                _mark('capheine')
                # 后处理: 密码子位点可视化 (visual_codon_miner, drhip 产物驱动)
                drhip_csv = os.path.join(out_dir, ph_dir('select'), "capheine", "drhip", "combined_sites.csv")
                cln_dir = os.path.join(out_dir, ph_dir('select'), "capheine", "hyphy", "CLN")
                if os.path.exists(drhip_csv) and os.path.isdir(cln_dir):
                    run_postprocess("visual_codon_miner", log,
                                    "--drhip", drhip_csv, "--clndir", cln_dir,
                                    "--outdir", os.path.join(out_dir, ph_dir('select'), "capheine", "codon_miner"))
                # 后处理: 正选择位点 3D 结构可视化 (batch_draw_pymol) 已改为按需手动运行, 不进管线
                # (基因名与 AlphaFold 资产库命名不匹配时全部跳过, 每次白跑; 需要时:
                #  python3 -m batch_draw_pymol --input <drhip_csv> --output <dir> --pdb_dir <dir>)
                # 后处理: 选择分析套件 (selection_suite FUBAR + selection_analysis)
                capheine_dir = os.path.join(out_dir, ph_dir('select'), "capheine")
                meta_for_sel = prep.get("meta_csv") or prep.get("dates_csv", "")
                run_postprocess("selection_suite", log,
                                "fubar", "--capheine-dir", capheine_dir,
                                "--metadata", meta_for_sel,
                                "--outdir", os.path.join(capheine_dir, "selection_out"))
                cawlign_dir = os.path.join(capheine_dir, "cawlign")
                if os.path.isdir(cawlign_dir):
                    # selection_analysis 用子命令式 CLI (parse/run/...); capheine 已有结果 → parse。
                    # 历史坑: 缺 parse 位置参数时 argparse 把 --capheine-dir 的值当 cmd → 永远失败
                    run_postprocess("selection_analysis", log,
                                    "parse",
                                    "--capheine-dir", capheine_dir,
                                    "--outdir", os.path.join(capheine_dir, "selection_out"))
                # 2026-08-30 决策: 正选择证据链以 CAPHEINE 为主 (用户拍板),
                # codeml 不进自动链, 且唯一实现 codeml_bridge.py 已随之归档到
                # _archive_codeml_20260830/ (不在 utils/ 下), selection_suite 的
                # codeml 子命令也已移除。如后续需复活, 先取回归档实现再重新接线。
            else:
                if cr.get("skip_reason"):
                    # 内层已定性不适用 (viroid 无 CDS 等), 记 skip 不记 fail (2026-09-02)
                    log.info(f"  [capheine] skip: {cr['skip_reason']}")
                    _skip('capheine')
                else:
                    _fail('capheine', 'capheine 未产出 (缺 GBK/参考 CDS 或计算失败)')

    # ── clock: 分子钟联合 (整体+分基因 RTT/LTT + DRT) ──
    if 'clock' in stages and aln_file and tree_file:
        if _done('clock'):
            log.info("  [clock] checkpoint OK, skip"); _skip('clock')
        else:
            prep["_cp"] = lambda name: _done(name)
            cr = run_stage_clock(prep, aln_file, tree_file, tools, args, log,
                                 results=results)
            results["outputs"].update(cr.get("outputs", {}))
            results["metrics"].update(cr.get("metrics", {}))
            if cr.get("outputs") or results.get("gene_results"):
                _mark('clock')
            else:
                _fail('clock', 'RTT/LTT/DRT 无输出')

    # rtt/temporal/genes 已并入 clock (此处不再单独跑; 保留函数供 --stage 兼容)

    # ── beast: 全长 BEAST 定年 ──
    if 'beast' in stages:
        if _done('beast'):
            log.info("  [beast] checkpoint OK, skip"); _skip('beast')
            beast_xml = os.path.join(out_dir, ph_dir('time'), "beast", "beast1.xml")
            if os.path.exists(beast_xml):
                results["outputs"]["beast_xml"] = beast_xml
            # skip 时也恢复 TMRCA → metrics (供 gene_dating --full-tmrca)
            beast_log = os.path.join(out_dir, ph_dir('time'), "beast", "beast1.log")
            if os.path.exists(beast_log):
                try:
                    parsed = parse_beast_log(beast_log, burnin_pct=getattr(args, 'beast_burnin', 10),
                                              log=LogCollector(log))
                    rh = parsed.get("key_params", {}).get("treeModel.rootHeight")
                    if rh:
                        results["metrics"]["beast_tmrca"] = rh["median"]
                        results["metrics"]["beast_tmrca_hpd"] = f"{rh['hpd_95_lower']:.1f}–{rh['hpd_95_upper']:.1f}"
                        results["metrics"]["beast_tmrca_ess"] = rh["ess"]
                        log.info(f"  [beast] 全长 TMRCA (复用): {rh['median']:.1f} yr "
                                 f"[{rh['hpd_95_lower']:.1f}–{rh['hpd_95_upper']:.1f}] ESS={rh['ess']:.0f}")
                except Exception as e:
                    log.warning(f"  [beast] skip 分支 TMRCA 解析失败: {e}")
        else:
            beast_results = run_phase_beast(prep, results, args, tools, log)
            # 2026-09-15 修复: run_phase_beast 返回的 dict 含同名 "outputs"/"metrics" 键,
            # 原 `results.update(beast_results)` 是整体替换 → beast 之前累积的
            # clock/popgen/capheine/geo/host 全部指标被清空 (汇总 CSV/HTML 只剩 beast 之后的数据)。
            # 改为深合并: 只并入新增键, outputs/metrics 逐键更新。
            for _k in ("outputs", "metrics"):
                _bv = beast_results.get(_k)
                if isinstance(_bv, dict):
                    results.setdefault(_k, {}).update(_bv)
            for _k in ("error", "beast_pending"):
                if _k in beast_results:
                    results[_k] = beast_results[_k]
            if beast_results.get("success"):
                _mark('beast')
            elif beast_results.get("beast_pending"):
                # 多链后台运行中: 不写 checkpoint, 下次跑检查合并
                _skip('beast')
            else:
                _fail('beast', beast_results.get("error", "failed"))

    # ── gene_dating: 分基因 BEAST 定年 ──
    if 'gene_dating' in stages:
        if _done('gene_dating'):
            log.info("  [gene_dating] checkpoint OK, skip"); _skip('gene_dating')
        else:
            run_stage_gene_dating(prep, results, args, log)
            if results.get("outputs", {}).get("gene_dating_dir"):
                _mark('gene_dating')
            else:
                # 无 capheine 前置产物 = 数据依赖缺失, 非计算失败
                _skip('gene_dating')

    # ── phylogeo: 系统地理 ──
    if 'phylogeo' in stages:
        if _done('phylogeo'):
            log.info("  [phylogeo] checkpoint OK, skip"); _skip('phylogeo')
        else:
            geo_pending = run_stage_phylogeo(prep, results, aln_file, args, log)
            if results.get("outputs", {}).get("phylogeo_xml"):
                if geo_pending:
                    # 链后台运行中: 不写 checkpoint, 下次跑检查合并 (同 beast pending 逻辑)
                    _skip('phylogeo')
                else:
                    _mark('phylogeo')
            else:
                _fail('phylogeo', '未生成 XML (可能缺数据)')
    # ── geo: 地理分析 (Mantel / 树地理)。位置关键：必须排在 phylogeo 之后 ——
    #    其 --rrt_randomized_dir 分支要吃 phylogeo 的 merged/mcc.tree ──
    if 'geo' in stages:
        if _done('geo'):
            log.info("  [geo] checkpoint OK, skip"); _skip('geo')
        else:
            run_stage_geo(prep, results, aln_file, tree_file, args, log)
            if results.get("outputs", {}).get("mantel_plot") or results.get("outputs", {}).get("geo_tree"):
                _mark('geo')
            else:
                _skip('geo')

    # ── geo_paths: 事件级传播路径 (2026-09-16 新增; 吃 phylogeo 的 mcc.tree) ──
    #    位置关键：**必须排在 phylogeo 之后**，否则首跑 mcc.tree 尚未生成 → 整段跳过。
    #    2026-09-18: 原先误排在 phylogeo 之前（STAGE_DEPS 声明了依赖却没人执行它）。
    if 'geo_paths' in stages:
        if _done('geo_paths'):
            log.info("  [geo_paths] checkpoint OK, skip"); _skip('geo_paths')
        else:
            run_stage_geo_paths(prep, results, aln_file, tree_file, args, log)
            if results.get("outputs", {}).get("pathway_qc"):
                _mark('geo_paths')
            else:
                _skip('geo_paths')


    results["stages_done"] = st["done"]
    results["stages_skipped"] = st["skipped"]
    results["stages_failed"] = st["failed"]
    # 2026-09-15 修复 (P1-B): 旧版 `success = bool(tree_file) or bool(outputs)`,
    # 只要任一 stage 往 results["outputs"] 写过东西就算成功 —— 于是 align/tree
    # 双双失败、仅 online 阶段写出 ncbi_fasta 时也会被判 success=True, 汇总表
    # 把"只跑完可选的序列补充"报成通过。现要求「无任何 stage 失败」且「有实质产物」。
    results["success"] = (not st["failed"]) and bool(
        tree_file or aln_file or results.get("outputs"))
    log.info(f"  [{vname}] stages: done={st['done']} skipped={st['skipped']} failed={st['failed']}")
    return results


def generate_summary(all_results, output_dir, logger) -> None:
    """汇总所有病毒的 stage 结果 -> phylo_summary.csv

    列语义 (2026-09-15 二次修订, 修 P0-A / P0-B):
      Beta / R_squared / P_value / TMRCA —— 取值一律经 utils.summary_metrics.
      headline_metrics(), 与 HTML 报告**同源同值**; 每个值配一列 *_source 标注
      实际来源(treedater / PIC / RTT / BEAST), 并额外落 alt_* 备用口径列。

      优先级: 速率 td_mean_rate(treedater) → pic_rate(PIC)
              R²  pic_r(PIC) → r_squared(RTT)
              P   td_relaxed_clock_p(treedater) → p_value(RTT)
              TMRCA beast_tmrca(BEAST) → td_tmrca(treedater)

      修订原因: ① 旧版 CSV 与 HTML 对同名列取不同来源(TreeTime RTT vs
      treedater/PIC), 同一病毒两份产物互相打脸且无来源标注; ② 旧版 Beta 列
      写 `td_mean_rate or pic_rate`, 把 treedater 松弛钟均值速率与 PIC 严格钟
      回归斜率两个不同估计量混进一列, 跨病毒比较时混入异质方法。现由
      headline_metrics 单点决定, 来源显式落列。
    """
    summary_path = os.path.join(output_dir, "phylo_summary.csv")
    rows = []
    for r in all_results:
        metrics = r.get("metrics", {})
        row = {
            "Virus": r.get("virus", ""),
            "Phase": r.get("phase", ""),
            "Success": r.get("success", False),
            "Stages_done": ";".join(r.get("stages_done", [])),
            "Stages_skipped": ";".join(r.get("stages_skipped", [])),
            "Stages_failed": ";".join(r.get("stages_failed", [])),
        }
        row.update(headline_columns(metrics))
        row.update({
            "Geo_n_pathways": metrics.get("n_pathways", ""),
            "Geo_n_episodes": metrics.get("n_episodes", ""),
            "Geo_conf_direct": metrics.get("pathway_conf_direct", ""),
            "Geo_conf_indirect": metrics.get("pathway_conf_indirect", ""),
            "Geo_conf_distant": metrics.get("pathway_conf_distant_import", ""),
            "Geo_velocity_km_yr": metrics.get("geo_mean_velocity_km_yr", ""),
            "Geo_wavefront_km": metrics.get("geo_max_wavefront_km", ""),
            "Geo_perm_p": metrics.get("geo_perm_p_displacement", ""),
            "Geo_root_height_yr": metrics.get("geo_root_height_yr", ""),
            "Geo_markov_jumps_mean": metrics.get("geo_markov_jumps_mean", ""),
            "Geo_markov_jumps_hpd95": metrics.get("geo_markov_jumps_hpd95", ""),
            "Timescale_ratio": metrics.get("timescale_rate_ratio", ""),
            "Timescale_flag": metrics.get("timescale_flag", ""),
            "Temporal_signal": metrics.get("temporal_signal", ""),
            # 2026-09-15 (P2-7): None (DRT 未产出结论) 渲染成 NA, 不再显示成 False
            "DRT_passed": ("NA" if metrics.get("drt_passed") is None
                           else metrics.get("drt_passed", "")),
            "DRT_status": metrics.get("drt_status", ""),
        })
        rows.append(row)
    if rows:
        with open(summary_path, 'w', newline='', encoding='utf-8') as f:
            writer = csv.DictWriter(f, fieldnames=rows[0].keys())
            writer.writeheader()
            writer.writerows(rows)
        logger.info(f"\nSummary: {summary_path}")


# ═══════════════════════════════════════════════════════════════════
# Main
# ═══════════════════════════════════════════════════════════════════

def main():
    cfg = load_config()
    defaults = config_defaults(cfg)
    # 启动期自检: config.yaml 里"登记不到 dest"的键会被静默丢弃 (历史踩坑
    # phylogeo.chains / beast.chains) —— 这里显式喊出来, 不让它继续悄悄失效。
    _orphans = config_orphan_keys(cfg)
    if _orphans:
        print("[config] ⚠ 以下键未登记映射, 会被静默忽略: "
              + ", ".join(_orphans)
              + "  (请在 CONFIG_SECTION_ARG 里登记)", file=sys.stderr)

    parser = argparse.ArgumentParser(
        description="MMPV-RNA PhyloPipeline v2 — 流程化系统发育分析 (VirPhyKit)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # 全流程批量分析 (从 06_extraction)
  python phylo_pipeline.py --extract_dir 06_extraction/ --metadata metadata.tsv --stage all -t 20

  # 指定 stage (自动补依赖: beast -> prep,align,tree,beast)
  python phylo_pipeline.py --extract_dir 06_extraction/ --metadata metadata.tsv --stage tree,beast

  # work-dir 模式 (已有 alignment+tree, 直接补时钟/定年/地理)
  python phylo_pipeline.py --work_dir gcva_results/ --alignment mafft.aln.fasta \\
      --tree iqtree.treefile --meta dates.csv --stage clock,beast

  # config 文件 (默认参数)
  python phylo_pipeline.py --config config.yaml --stage all

Stages: """ + ", ".join(f"{k}({v})" for k, v in STAGE_DESC.items()),
    )

    # ── Config ──
    g_cfg = parser.add_argument_group("Config")
    g_cfg.add_argument("--config", default=None, help="config.yaml 路径 (默认: pipeline 根目录 config.yaml)")

    # ── Input sources ──
    g_input = parser.add_argument_group("Input")
    g_input.add_argument("--extract_dir", default=None, help="06_extraction 目录 (批量模式; standard 布局为 04_Analysis/06_Extraction)")
    g_input.add_argument("--io-layout", choices=["legacy", "standard"], default=None,
                         help="I/O 目录布局: legacy=v3.0 现行名 (默认); standard=管线独立根"
                              "(05_Phylo 内模块编号 01_data..08_report, 见 doc/IO_LAYOUT_DESIGN.md); "
                              "亦可 MMPV_IO_LAYOUT 环境变量")
    g_input.add_argument("--metadata", default=None, help="Global_Unified_Metadata_Core14.tsv 路径")
    g_input.add_argument("--info_dir", default=None, help="public_metadata_pipeline 的 info/ 输出目录")
    g_input.add_argument("--variants_dir", default=None, help="03_variants 目录 (查找参考序列)")
    g_input.add_argument("--virus", default=None,
                        help="单个病毒名 (datasets.yaml 注册名, 如 GCVA/PSTVD; 自动定位比对/元数据)")
    g_input.add_argument("--work_dir", default=None, help="单个病毒工作目录 (已有 alignment+tree 的模式)")

    # ── Work dir mode ──
    g_work = parser.add_argument_group("Work-dir mode (已有 alignment+tree)")
    g_work.add_argument("--alignment", default=None, help="已有 MAFFT 比对文件")
    g_work.add_argument("--tree", default=None, help="已有 Newick 树文件")
    g_work.add_argument("--meta", default=None, help="已有 dates.csv (用于时间分析)")
    g_work.add_argument("--host_meta", default=None,
                        help="宿主元数据 CSV (name,host[,tissue]) —— 会独立成宿主通道 "
                             "data/host.csv, 只收有 host 的样品; 省略则由 --meta 拆出")

    # ── Output ──
    g_out = parser.add_argument_group("Output")
    g_out.add_argument("--output_dir", "-o", default=defaults.get("output_dir", None),
                       help="输出根目录 (默认: phylo_results/)")

    # ── Runtime ──
    g_run = parser.add_argument_group("Runtime")
    g_run.add_argument("--stage", default=defaults.get("stage", "all"),
                       help="运行 stage: 大模块 " + ",".join(STAGE_GROUPS) + " (或子步骤 " + ",".join(STAGE_ORDER) + ")")
    g_run.add_argument("--disable", default=defaults.get("disable", ""),
                       help="真正关闭 stage (逗号分隔, 从最终清单里剔除)。"
                            "用于覆盖 --stage all / 大模块展开: 如 --disable beast,rdp5")
    g_run.add_argument("--threads", "-t", type=int, default=defaults.get("threads", 40))
    g_run.add_argument("--min_samples", type=int, default=defaults.get("min_samples", 5))
    g_run.add_argument("--dedup", action=argparse.BooleanOptionalAction,
                       default=defaults.get("dedup", True),
                       help="data 阶段 (date, location, seq) 三元组克隆去重 (--no-dedup 关闭)")
    g_run.add_argument("--force", action="store_true", default=defaults.get("force", False))
    g_run.add_argument("--skip_msa", action="store_true", default=defaults.get("skip_msa", False))
    g_run.add_argument("--skip_tree", action="store_true", default=defaults.get("skip_tree", False))
    g_run.add_argument("--max_viruses", type=int, default=defaults.get("max_viruses", 0))
    g_run.add_argument("--iqtree_bin", default=defaults.get("iqtree_bin", "iqtree2"))
    # ── IQ-TREE 分支支持值 (2026-09-16 补) ──
    # 默认与平台侧 `Virus_Platform_Core/phylo.py::_run_iqtree` 同口径。
    # ⚠️ `-B`/`-alrt` 在 IQ-TREE 里硬性要求 >=1000, "少跑几次" 会被拒绝;
    # 传 0 才表示不要该支持值 (两者都 0 = 回到修复前行为)。
    g_run.add_argument("--iqtree_boot", type=int,
                       default=defaults.get("iqtree_boot", IQTREE_BOOT_DEFAULT),
                       help="UFBoot 重复数 (-B); 默认 1000, 传 0 关闭 (IQ-TREE 要求 >=1000)")
    g_run.add_argument("--iqtree_alrt", type=int,
                       default=defaults.get("iqtree_alrt", IQTREE_ALRT_DEFAULT),
                       help="SH-aLRT 重复数 (-alrt); 默认 1000, 传 0 关闭 (要求 >=1000)")
    g_run.add_argument("--iqtree_bnni", action="store_true",
                       default=defaults.get("iqtree_bnni", False),
                       help="给 UFBoot 加 --bnni (近缘序列多时压低 UFBoot 高估, 稍慢)")
    g_run.add_argument("--iqtree_timeout", type=int,
                       default=defaults.get("iqtree_timeout", IQTREE_TIMEOUT_DEFAULT),
                       help="IQ-TREE **单次尝试**超时秒数; 默认 3600, 大比对需上调")
    g_run.add_argument("--dry", action="store_true", help="只打印执行计划, 不运行")

    # ── 分析前序列清洗 (clean stage, 比对前第一道) ──
    # 2026-09-16 新增。设计原则: 默认值保守且**不改行为** —— 所有阈值都给成
    # "刚好拦住明显异常" 的水平, 且元数据筛除默认关闭, 避免静默改变既有结果。
    g_clean = parser.add_argument_group(
        "分析前序列清洗 (clean stage, 比对前；与比对后的 align_qc 互补)")
    g_clean.add_argument("--clean_min_length_ratio", type=float,
                         default=defaults.get("clean_min_length_ratio", 0.9),
                         help="长度 / 基准 的下限 (默认 0.9)")
    g_clean.add_argument("--clean_max_length_ratio", type=float,
                         default=defaults.get("clean_max_length_ratio", 1.5),
                         help="长度 / 基准 的上限, 抓异常长嵌合序列 (默认 1.5)")
    g_clean.add_argument("--clean_max_n", type=float,
                         default=defaults.get("clean_max_n", 0.05),
                         help="N + 简并码占比上限 (默认 0.05)")
    g_clean.add_argument("--clean_reference", default=defaults.get("clean_reference", None),
                         help="参考序列 ID; 给了则用它作长度基准, 否则用中位数")
    g_clean.add_argument("--clean_allow_rna", action="store_true",
                         default=defaults.get("clean_allow_rna", False),
                         help="允许 U (RNA 病毒)")
    g_clean.add_argument("--clean_strip_illegal", action="store_true",
                         default=defaults.get("clean_strip_illegal", False),
                         help="删除非法字符而非整条剔除 (默认整条剔除, 更保守)")
    g_clean.add_argument("--clean_allow_gap_chars", action="store_true",
                         default=defaults.get("clean_allow_gap_chars", False),
                         help="允许未比对序列中出现 '-'/'.' (默认视为异常剔除)")
    g_clean.add_argument("--clean_metadata", default=defaults.get("clean_metadata", None),
                         help="metadata CSV/TSV (用于时间/地理筛除)")
    g_clean.add_argument("--clean_drop_no_date", action="store_true",
                         default=defaults.get("clean_drop_no_date", False),
                         help="剔除日期为空的序列")
    g_clean.add_argument("--clean_drop_no_location", action="store_true",
                         default=defaults.get("clean_drop_no_location", False),
                         help="剔除地理为空的序列")

    # ── metadata 入口治理 (metadata stage, 全流程最上游) ──
    # 2026-09-16 新增 (老师拍板)。产出另存新文件, 原件不动; 失败降级用原文件。
    g_md = parser.add_argument_group(
        "metadata 入口治理 (metadata stage；时间+地理 的检查与矫正)")
    g_md.add_argument("--metadata_mode", choices=("mid", "start"),
                      default=defaults.get("metadata_mode", "mid"),
                      help="日期口径: mid=月中/年中 (decimal_year); "
                           "start=月初/年初 (virphy_bridge/VirPhyKit 原口径)。"
                           "注意: 本 stage 的产物供下游消费, 与 virphy_bridge 内部"
                           "的 mode='start' 是两个独立选择, 默认 mid 保持不变")
    g_md.add_argument("--metadata_date_col", default=defaults.get("metadata_date_col", None),
                      help="强制指定日期列名 (默认自动识别)")
    g_md.add_argument("--metadata_loc_col", default=defaults.get("metadata_loc_col", None),
                      help="强制指定地理列名 (默认自动识别)")
    g_md.add_argument("--metadata_name_col", default=defaults.get("metadata_name_col", None),
                      help="强制指定名称列名 (默认自动识别)")
    g_md.add_argument("--metadata_keep_unknown_geo", action="store_true",
                      default=defaults.get("metadata_keep_unknown_geo", False),
                      help="保留地理中的 Unknown 层级 (默认删除, 老师拍板)")
    g_md.add_argument("--metadata_no_derive_coords", action="store_true",
                      default=defaults.get("metadata_no_derive_coords", False),
                      help="**不补算**坐标 (默认: 坐标占位符先用地名推算, 推不出才清空)")
    g_md.add_argument("--metadata_no_blank_placeholders", action="store_true",
                      default=defaults.get("metadata_no_blank_placeholders", False),
                      help="**不清空**占位符字段 (默认清空: 老师拍板「仅清空该字段」)")
    g_md.add_argument("--metadata_drop_placeholder_rows", action="store_true",
                      default=defaults.get("metadata_drop_placeholder_rows", False),
                      help="含占位符字段的行整行剔除 (默认不剔, 只清空字段)")
    g_md.add_argument("--metadata_drop_no_date", action="store_true",
                      default=defaults.get("metadata_drop_no_date", False),
                      help="治理时剔除日期不可解析的行")
    g_md.add_argument("--metadata_drop_no_location", action="store_true",
                      default=defaults.get("metadata_drop_no_location", False),
                      help="治理时剔除地理为空的行")

    # ── Recombination (RDP5, Phase R) ──
    g_rdp5 = parser.add_argument_group("Recombination (RDP5, Phase R)")
    g_rdp5.add_argument("--rdp5", action="store_true", default=defaults.get("rdp5", False),
                        help="MAFFT 后、IQ-TREE 前运行 RDP5 重组检测 + mask 重组区")
    g_rdp5.add_argument("--rdp5_script", default=defaults.get("rdp5_script", "~/MMPV-RNA/biosoft/rdp5/run_rdp5.sh"))
    g_rdp5.add_argument("--rdp5_genes", default=defaults.get("rdp5_genes", None))
    g_rdp5.add_argument("--mask_min_methods", type=int, default=defaults.get("mask_min_methods", 3))
    g_rdp5.add_argument("--rdp5_mask", action="store_true", default=defaults.get("rdp5_mask", False),
                        help="重组处理改用 mask 区间置N (不删序列不重比对); 默认删重组子+重新MAFFT")
    g_rdp5.add_argument("--rdp5_stage", default=defaults.get("rdp5_stage", "run"), choices=["run", "all"])

    # ── BEAST (Phase 5) ──
    g_beast = parser.add_argument_group("BEAST (Phase 5)")
    g_beast.add_argument("--beast", action="store_true", default=defaults.get("beast", False),
                         help="运行全长 BEAST 分子定年 (= stage 含 beast)")
    g_beast.add_argument("--beast_bin", default=defaults.get("beast_bin", None))
    g_beast.add_argument("--beast_chain", type=int, default=defaults.get("beast_chain", 5_000_000))
    g_beast.add_argument("--beast_burnin", type=float, default=defaults.get("beast_burnin", 10.0))
    g_beast.add_argument("--beast_timeout", type=int, default=defaults.get("beast_timeout", 48))
    g_beast.add_argument("--beast_chains", type=int, default=defaults.get("beast_chains", 3),
                         help="全长 BEAST 并行链数 (默认 3, 发表级多链; 1 = 单链前台等待)")
    g_beast.add_argument("--skip_beast_run", action="store_true", default=defaults.get("skip_beast_run", False))

    # ── Per-gene BEAST (Phase G) ──
    g_gene = parser.add_argument_group("Per-gene BEAST dating (Phase G)")
    g_gene.add_argument("--gene_dating", action="store_true", default=defaults.get("gene_dating", False),
                        help="运行分基因 BEAST 定年 (吃 capheine 产物)")
    g_gene.add_argument("--genes_dir", default=defaults.get("genes_dir", None),
                        help="基因比对目录 (capheine/cawlign 产物; 默认自动找)")
    g_gene.add_argument("--gene_chain", type=int, default=defaults.get("gene_chain", 15_000_000))
    g_gene.add_argument("--gene_chain_big", type=int, default=defaults.get("gene_chain_big", 10_000_000))
    g_gene.add_argument("--big_genes", default=defaults.get("big_genes", "L"))

    # ── Saturation (数据质量) ──
    g_sat = parser.add_argument_group("Saturation (替换饱和)")
    # 历史坑: 全文件唯一不走 config 默认值的参数, 改 saturation.replicates 静默失效
    g_sat.add_argument("--sat_replicates", type=int,
                       default=defaults.get("sat_replicates", 1000),
                       help="Xia Iss 重采样次数 (默认: config saturation.replicates, 再默认 1000)")

    # ── Popgen (Phase P) ──
    g_pop = parser.add_argument_group("Popgen (Phase P, pypopart)")
    g_pop.add_argument("--popgen_group", default=defaults.get("popgen_group", "location"),
                       help="分组列 (默认: location)")
    g_pop.add_argument("--popgen_min_n", type=int, default=defaults.get("popgen_min_n", 5))
    g_pop.add_argument("--popgen_perm", type=int, default=defaults.get("popgen_perm", 999))
    # 2026-09-15 (审查 P2-15): 参考/外群默认参与 π/S/θ/Hd 会把多样性系统性抬高。
    # 默认行为保持不变 (不排除), 但提供开关; popgen 侧对"检测到但未排除"会告警。
    g_pop.add_argument("--popgen_exclude_reference", action="store_true",
                       default=defaults.get("popgen_exclude_reference", False),
                       help="群体遗传学: 自动排除参考基因组/外群序列")
    g_pop.add_argument("--popgen_exclude_ids", default=defaults.get("popgen_exclude_ids", None),
                       help="群体遗传学: 额外排除的序列 id (逗号分隔)")

    # ── Capheine (Phase S) ──
    g_cap = parser.add_argument_group("Capheine (Phase S, positive selection)")
    g_cap.add_argument("--capheine_ref", default=defaults.get("capheine_ref", None),
                       help="参考 CDS FASTA (缺省自动从 work_dir GBK 提取)")
    g_cap.add_argument("--capheine_unaligned", default=defaults.get("capheine_unaligned", None),
                       help="未比对序列 FASTA (缺省用 prep 序列)")
    g_cap.add_argument("--capheine_code", default=defaults.get("capheine_code", "1"),
                       help="遗传密码 (默认: 1=Universal)")
    g_cap.add_argument("--capheine_workers", type=int, default=defaults.get("capheine_workers", 4),
                       help="并行基因数 (默认: 4, 2026-09-03 提升; 246 服务器 256 核充裕)")
    g_cap.add_argument("--capheine_cpus_iqtree", type=int, default=defaults.get("capheine_cpus_iqtree", 6))
    g_cap.add_argument("--capheine_cpus_hyphy", type=int, default=defaults.get("capheine_cpus_hyphy", 32),
                       help="HyPhy MPI rank 数 (默认: 32, 2026-09-03 提升; PRIME 按位点切分近线性加速, 单 rank 内存小)")
    g_cap.add_argument("--capheine_mpi", action=argparse.BooleanOptionalAction, default=defaults.get("capheine_mpi", True),
                       help="HyPhy MPI")

    # ── Phylogeography (Phase 6) ──
    g_geo = parser.add_argument_group("Phylogeography (Phase 6)")
    g_geo.add_argument("--phylogeo", action="store_true", default=defaults.get("phylogeo", False))
    g_geo.add_argument("--phylogeo_clock", default=defaults.get("phylogeo_clock", "ucln"), choices=["strict", "ucln"])
    g_geo.add_argument("--phylogeo_prior", default=defaults.get("phylogeo_prior", "skyline"),
                       choices=["constant", "skyline", "bdsky", "auto"])
    g_geo.add_argument("--phylogeo_chain", type=int, default=defaults.get("phylogeo_chain", 10_000_000))
    g_geo.add_argument("--phylogeo_bf", type=float, default=defaults.get("phylogeo_bf", 5.0))
    g_geo.add_argument("--phylogeo_chains", type=int, default=defaults.get("phylogeo_chains", 3),
                       help="phylogeo BEAST 并行链数 (默认 3, 发表级多链)")
    g_geo.add_argument("--skip_phylogeo_run", action="store_true",
                       default=defaults.get("skip_phylogeo_run", False),
                       help="跳过 phylogeo 多链 BEAST 提交 (默认生成 XML 后即提交)")
    g_geo.add_argument("--geo_subsample", type=int, default=defaults.get("geo_subsample", 0),
                       help="geo 前按地点均衡抽样数 (0=不抽样, 用全量)")
    g_geo.add_argument("--treeannotator_bin", default=defaults.get("treeannotator_bin", "treeannotator"))
    g_geo.add_argument("--beast_version", default=defaults.get("beast_version", "1"), choices=["1", "2"])
    # ── RRT 树注释法 (VirPhyKit 同口径, 2026-09-16 接线) ──
    # 需外部准备好 N 棵"区域随机化"MCC 树 (上游口径 N=20, 需 ~21 次完整 BEAST 重跑),
    # 因此**默认关闭**: 不给该参数就完全不跑, 也不会改变既有地理结论。
    g_geo.add_argument("--rrt_randomized_dir", default=defaults.get("rrt_randomized_dir", None),
                       help="随机化 MCC 树目录/glob; 给了才跑 VirPhyKit 口径 RRT(树注释法). "
                            "树需自备: 打乱地点标签重跑 phylogeo BEAST + treeannotator 取 MCC, 重复 N 次")
    # ── geo_paths 子阶段 (2026-09-16 新增): 事件级传播路径 + 扩散统计 + 传播图 ──
    # ⚠ config.yaml 键必须同步登记 CONFIG_SECTION_ARG['geo_paths'] (chains 键静默丢失的教训)
    g_geo.add_argument("--no_pathways", action="store_true",
                       default=defaults.get("no_pathways", False),
                       help="关闭 geo_paths 子阶段 (默认在 phylogeo MCC 存在时自动运行)")
    g_geo.add_argument("--direct_snp", type=float, default=defaults.get("direct_snp", 2.0),
                       help="三级置信 direct 阈值(期望替换数; 默认沿用 phymapr, 论文前须标定)")
    g_geo.add_argument("--indirect_snp", type=float, default=defaults.get("indirect_snp", 5.0),
                       help="三级置信 indirect 阈值(期望替换数)")
    g_geo.add_argument("--censor_years", type=float, default=defaults.get("censor_years", 0.5),
                       help="LTL 持久性删失窗(年; 默认 0.5, polio-wpv1 同款)")
    g_geo.add_argument("--n_perm", type=int, default=defaults.get("n_perm", 200),
                       help="扩散统计置换零模型次数 (坐标重标注方案, 种子固定 42)")
    g_geo.add_argument("--calibrate", default=defaults.get("calibrate", "manual"),
                       choices=["auto", "manual"],
                       help="置信阈值标定: auto=本数据 E[S] 分布 q25/q75 (推荐, 阈值随数据尺度自适应); "
                            "manual=用 --direct_snp/--indirect_snp 显式值")
    g_geo.add_argument("--geo_map", action="store_true",
                       default=defaults.get("geo_map", False),
                       help="geo_paths 附带发表级静态传播图 (PNG+PDF+SVG, 默认纯矢量合规底图)")
    g_geo.add_argument("--geo_gif", action="store_true",
                       default=defaults.get("geo_gif", False),
                       help="geo_paths 附带事件级时间揭示动画 GIF (需 --geo_map)")
    g_geo.add_argument("--geo_gif_frames", type=int,
                       default=defaults.get("geo_gif_frames", 36),
                       help="动画 GIF 帧数 (默认 36)")
    g_geo.add_argument("--geo_tree_map", action="store_true",
                       default=defaults.get("geo_tree_map", False),
                       help="geo_paths 附带树-图联动双面板 (phymapr Tree|Map 图版的 matplotlib 替代)")

    # ── Temporal signal (Phase 0) ──
    g_drt = parser.add_argument_group("Temporal signal (Phase 0)")
    g_drt.add_argument("--check_temporal", action="store_true", default=defaults.get("check_temporal", False),
                       help="运行时间信号检验 (DRT)")
    g_drt.add_argument("--drt_randomizations", type=int, default=defaults.get("drt_randomizations", 20))
    g_drt.add_argument("--bets", action="store_true", default=defaults.get("bets", False),
                       help="运行 BETS 时间信号检验")
    g_drt.add_argument("--bets_steps", type=int, default=defaults.get("bets_steps", 20))
    g_drt.add_argument("--bets_chain", type=int, default=defaults.get("bets_chain", 1_000_000))

    # ── Report (汇总 CSV + HTML) ──
    # 2026-09-15 (P1-7): config 的 report.enabled 此前根本没注册成 argparse 参数
    # (只有 CONFIG_SECTION_ARG 里的映射), 因此这个配置项是死的。现补上开关。
    g_rep = parser.add_argument_group("Report")
    g_rep.add_argument("--report_enabled", action="store_true",
                       default=defaults.get("report_enabled", True),
                       help="生成 phylo_summary.csv + phylo_report.html (config: report.enabled)")
    g_rep.add_argument("--no-report", dest="report_enabled", action="store_false",
                       help="不生成汇总 CSV / HTML 报告")

    args = parser.parse_args()

    # 布局解析 (CLI > MMPV_IO_LAYOUT > legacy) 并写回环境供子进程继承
    normalize_layout_env(getattr(args, "io_layout", None))
    if not getattr(args, "output_dir", None):
        args.output_dir = dir_name("p_root")   # legacy=phylo_results / standard=05_Phylo

    # ── 重载 config (若 --config 显式给出) ──
    # 2026-09-15 修复 (P1-5): 旧版无条件 `setattr(args, dest, val)`, 注释却写着
    # "CLI 显式传参优先" —— 判断根本没实现。只要显式传 --config, 命令行上的
    # -t / --stage / --force / --output_dir 等全部被 yaml 静默盖掉。
    # 现先算出"命令行上真正显式出现过"的 dest 集合, 这些 dest 不被 config 覆盖。
    explicit_cli = _explicit_cli_dests(parser)
    if args.config:
        cfg2 = load_config(args.config)
        if cfg2:
            _overridden = []
            for dest, val in config_defaults(cfg2).items():
                if not hasattr(args, dest):
                    continue
                if dest in explicit_cli:
                    _overridden.append(dest)
                    continue          # CLI 显式传参优先
                setattr(args, dest, val)
            if _overridden:
                print(f"[config] {args.config}: 以下参数在命令行上显式给出, "
                      f"不被 config 覆盖: {', '.join(sorted(_overridden))}")

    # ── 解析 stage 清单 (开关参数合并) ──
    # 2026-09-15 修复 (P1-6): config 的 <module>.enabled 现在是**真正的关闭开关**。
    # 旧版它只被用来"把 stage 追加进清单", 而 `--stage all` 展开的 STAGE_ORDER
    # 本身已含 rdp5/beast/gene_dating/phylogeo —— 于是 enabled:false 完全不起作用,
    # 用户读 config 以为关掉了、实际照跑, 只在末尾打一行告警。
    # 现在: enabled=false 且用户既没在 --stage 里点名、也没在命令行显式打开时,
    # 该 stage 会被移出执行清单 (并打印被移除的项, 不静默)。
    # ⚠️ 配套改动: config.yaml 里这四项的 enabled 已同步改为 true, 以保持
    #    `--stage all` 的实际行为与本次改动前完全一致。
    _enable_dests = ('rdp5', 'beast', 'gene_dating', 'phylogeo')

    # 用户通过 --stage 点名的 stage (展开大模块名), 用于豁免 enabled 门控
    _user_stages = set()
    _raw_stage = '' if args.stage is None else str(args.stage)
    if _raw_stage.strip() not in ('', 'all'):
        for _s in _raw_stage.split(','):
            _s = _s.strip()
            if not _s:
                continue
            if _s in STAGE_GROUPS:
                _user_stages.update(STAGE_GROUPS[_s])
            else:
                _user_stages.add(_s)

    # 2026-09-15 (P1-6): 只有**命令行上显式给出**的开关才追加 stage。
    # 旧版直接看 `args.X`, 而 `args.X` 的默认值来自 config 的 <module>.enabled ——
    # 于是 config 写 enabled:true 会把该 stage 隐式塞进**任何**部分运行的清单
    # (`--stage clock` 也会顺带跑 beast), 把"默认全集成员资格"和"本次显式请求"
    # 混成了一个开关, 也正是 P1-6 那个"看得见但语义矛盾"的根源。拆分后:
    #   config <module>.enabled      -> 该 stage 是否属于 `--stage all` 的默认全集
    #   --rdp5/--beast/... 命令行开关 -> 本次显式请求, 追加进执行清单
    stage_extra = []
    for _mod in _enable_dests:
        if _mod in explicit_cli and getattr(args, _mod, False):
            stage_extra.append(_mod)
    if args.check_temporal or args.bets:
        stage_extra.append('clock')  # 2026-08-27: rtt/temporal/genes 已合并为 clock

    stage_arg = args.stage
    if stage_extra:
        if stage_arg in (None, '', 'all'):
            stage_arg = 'all'
        else:
            stage_arg = stage_arg + ',' + ','.join(stage_extra)
    stages = resolve_stages(stage_arg)

    # ── --disable: 真正的关闭开关 ──
    _disabled = [s.strip() for s in (getattr(args, 'disable', '') or '').split(',') if s.strip()]
    if _disabled:
        _unknown = [s for s in _disabled if s not in STAGE_ORDER]
        if _unknown:
            sys.exit(f"--disable 含未知 stage: {_unknown}; "
                     f"可用: {', '.join(STAGE_ORDER)}")
        _before = list(stages)
        stages = [s for s in stages if s not in _disabled]
        _removed = [s for s in _before if s not in stages]
        print(f"[stage] --disable 生效, 已移除: {', '.join(_removed) or '(无)'}")

    # ── config <module>.enabled:false -> 真门控 (P1-6) ──
    _cfg_off = []
    for _mod in _enable_dests:
        if _mod not in stages:
            continue
        if getattr(args, _mod, False):
            continue                       # config / 命令行显式打开
        if _mod in explicit_cli:
            continue                       # 命令行上显式给了该开关
        if _mod in _user_stages:
            continue                       # --stage 里点名要求
        _cfg_off.append(_mod)
    if _cfg_off:
        stages = [s for s in stages if s not in _cfg_off]
        print(f"[stage] config 里 {', '.join(_cfg_off)}.enabled=false -> 已移出执行清单。"
              f"如需运行: 把 config 改为 enabled: true, 或在 --stage 中显式点名, "
              f"或加对应的命令行开关。")

    # ── Setup ──
    output_dir = Path(args.output_dir).resolve()
    logger = setup_logger(str(output_dir))
    logger.info("=" * 60)
    logger.info("MMPV-RNA PhyloPipeline v2 — 流程化 stage 编排")
    logger.info(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    logger.info(f"Stages: {', '.join(stages)}")
    logger.info("=" * 60)

    if args.dry:
        print(f"[dry] stages 执行计划: {', '.join(stages)}")
        for s in stages:
            print(f"  {s:<12} {STAGE_DESC[s]}")
        return

    # ── Check dependencies ──
    logger.info("\n[Check] VirPhyKit dependencies:")
    for name, ok in check_dependencies().items():
        logger.info(f"  {name}: {'✓' if ok else '✗'}")
    logger.info("\n[Check] External tools:")
    tools = check_tools(logger)

    # ── Single-virus mode (datasets.yaml 自动定位) ──
    if args.virus:
        key = args.virus.upper()
        ds = get_dataset(key)
        if not ds or not ds.get("dir"):
            avail = ', '.join(known_dataset_names()) or '(无注册数据集)'
            logger.error(f"未知病毒 '{args.virus}'. 已知: {avail}")
            sys.exit(1)
        data_dir = str(ds["dir"])
        fasta = os.path.join(data_dir, ds.get("fasta", ""))
        metadata = os.path.join(data_dir, ds.get("metadata", ""))
        # fasta 空/缺失时 fallback 到同目录最大的有效比对
        if not os.path.exists(fasta) or os.path.getsize(fasta) == 0:
            cands = sorted(
                [f for f in os.listdir(data_dir)
                 if f.endswith(('.aln.fasta', '.fasta', '.fas'))],
                key=lambda f: os.path.getsize(os.path.join(data_dir, f)), reverse=True)
            for c in cands:
                p = os.path.join(data_dir, c)
                if os.path.getsize(p) > 1000:
                    logger.warning(f"datasets.yaml fasta 空/缺失 ({Path(fasta).name}), 用 fallback: {c}")
                    fasta = p
                    break
        if not os.path.exists(fasta) or os.path.getsize(fasta) == 0:
            logger.error(f"比对文件不存在或为空: {fasta}")
            sys.exit(1)
        if not os.path.exists(metadata):
            logger.error(f"元数据不存在: {metadata}")
            sys.exit(1)
        virus_name = os.path.basename(os.path.dirname(data_dir.rstrip('/\\')))
        tree_file = None
        for cand in [os.path.join(data_dir, ph_dir('phylogeny'), "iqtree.treefile"),
                     os.path.join(data_dir, "iqtree.treefile")]:
            if os.path.exists(cand):
                tree_file = cand
                break

        logger.info(f"\n[Mode] Single virus (datasets.yaml): {key} → {data_dir}")
        logger.info(f"  alignment: {Path(fasta).name}, metadata: {Path(metadata).name}, "
                    f"tree: {Path(tree_file).name if tree_file else '(待建)'}")

        prep = {
            "virus_name": virus_name,
            "dataset_key": key,
            "combined_fasta": fasta,
            "dates_csv": metadata,
            "meta_csv": metadata,
            # 宿主通道: datasets.yaml 只指一份 metadata → 先留空, 由
            # `_refresh_metadata_channels` 从它拆出 `data/host.csv` (含 host 列才落地)
            "host_csv": "",
            "tree_file": tree_file,
            "output_dir": data_dir,
            "accession": ds.get("accession", ""),
            "species": ds.get("species", ""),
            "annotation_gb": ds.get("annotation_gb", ""),
            "pdb_dir": ds.get("pdb_dir", ""),
            "avg_seq_len": 10000,
            "n_samples": 0,
            "n_with_dates": 0,
            "has_location": True,
        }
        # ── 去重闸 (2026-09-03): 本地+线上合并后按 (date, location, seq) 去除克隆重复 ──
        # 参考 geo_analysis smart subsample 风格; 仅在去重产物不存在时执行 (checkpoint 语义)。
        # 注意: 已有 mafft 比对时不替换 combined_fasta (避免 checkpoint 失配), 报告仍生成。
        if getattr(args, 'dedup', True) and not os.path.exists(os.path.join(data_dir, ph_dir('data'), 'dedup_report.csv')):
            _aln_exists = os.path.exists(os.path.join(data_dir, ph_dir('phylogeny'), 'mafft.aln.fasta'))
            try:
                from utils.data_collector import dedup_triplet
                _dd = dedup_triplet(fasta, metadata, os.path.join(data_dir, ph_dir('data')), logger=logger)
                if _dd.get('n_removed', 0) > 0:
                    if not _aln_exists:
                        prep['combined_fasta'] = _dd['fasta']
                        logger.info(f"  [dedup] 下游改用去重后序列: {Path(_dd['fasta']).name}")
                    else:
                        logger.warning("  [dedup] 检出克隆重复但比对已存在, 不替换 (重跑需删比对产物)")
            except Exception as _e:
                logger.warning(f"  [dedup] 去重失败 ({_e}), 使用原序列")

        results = process_virus(prep, stages, args, tools, logger)

        report_dir = os.path.join(data_dir, ph_dir('report'))
        os.makedirs(report_dir, exist_ok=True)
        # 2026-09-15 (P1-7): report 由 main() 收尾统一产出, 不是 process_virus 里的
        # stage 分支 —— 故它是"总是产出"而非门控。这里接上 config 的 report.enabled
        # (--report_enabled / --no-report), 让该配置项真正生效。
        if getattr(args, "report_enabled", True):
            generate_summary([results], report_dir, logger)
            if results:
                report_path = build_html_report([results], report_dir)
                logger.info(f"  HTML report: {report_path}")
        else:
            logger.info("  [report] report.enabled=false, 跳过 phylo_summary.csv / HTML")
        # 后处理: 发表级三件套 (TempEst 图 + MCC 树 + 参数附表)
        # 历史坑: beast_log 硬编码 beast1.log, 多链模式 (beast_chains>1) 下该文件不存在 →
        # Table S1 静默缺失; 且从未传 --mcc_tree → Fig 2 MCC 树静默缺失
        _beast_dir = os.path.join(data_dir, ph_dir('time'), "beast")
        _beast_log = None
        for _cand in (os.path.join(_beast_dir, "merged", "merged.log"),
                      os.path.join(_beast_dir, "beast1.log"),
                      os.path.join(_beast_dir, "beast11.log")):
            if os.path.exists(_cand):
                _beast_log = _cand; break
        _mcc_tree = os.path.join(_beast_dir, "merged", "mcc.tree")
        # 后处理 metadata 一律走「治理产物优先」(2026-09-16): 此前直接用治理前
        # 捕获的 `metadata`, 等于 TempEst 日期 / 采样地图仍吃脏值。
        _pp = _governed_meta(prep, metadata)
        _pf_cmd = ["--rtt_fasta", fasta,
                        "--rtt_dates", _pp["dates"],
                        "--rtt_tree", tree_file or ""]
        if _beast_log:
            _pf_cmd += ["--beast_log", _beast_log]
        if os.path.exists(_mcc_tree):
            _pf_cmd += ["--mcc_tree", _mcc_tree]
        _pf_cmd += ["--output_dir", os.path.join(report_dir, "publication")]
        run_postprocess("publication_figures", logger, *_pf_cmd)
        # 后处理: 采样分布图 (时间线 + 地理地图)
        run_postprocess("sample_plot", logger,
                        "--metadata", _pp["meta"],
                        "--outdir", os.path.join(report_dir, "sample_plot"))
        logger.info(f"\nDone. Results: {data_dir}")
        return

    # ── Work-dir mode (单病毒) ──
    if args.work_dir and args.alignment and args.tree and args.meta:
        work_dir = Path(args.work_dir).resolve()
        logger.info(f"\n[Mode] Single virus work-dir: {work_dir}")

        def _resolve(p):
            pp = Path(p)
            return str(pp) if pp.is_absolute() else str(work_dir / pp)

        prep = {
            "virus_name": work_dir.name,
            "combined_fasta": _resolve(args.alignment),
            "dates_csv": _resolve(args.meta),
            "meta_csv": _resolve(args.host_meta) if args.host_meta else _resolve(args.meta),
            # 宿主通道: `--host_meta` 若已给出, 它天然含 host 列 → 正好当宿主骨架;
            # 否则让 `_refresh_metadata_channels` 从 --meta 拆出 host.csv。
            "host_csv": _resolve(args.host_meta) if args.host_meta else "",
            "tree_file": _resolve(args.tree),
            "output_dir": str(work_dir),
            "avg_seq_len": 10000,
            "n_samples": 0,
            "n_with_dates": 0,
            "has_location": True,
        }
        results = process_virus(prep, stages, args, tools, logger)

        report_dir = os.path.join(str(work_dir), "report")
        os.makedirs(report_dir, exist_ok=True)
        if getattr(args, "report_enabled", True):
            generate_summary([results], report_dir, logger)
            if results:
                report_path = build_html_report([results], report_dir)
                logger.info(f"  HTML report: {report_path}")
        else:
            logger.info("  [report] report.enabled=false, 跳过 phylo_summary.csv / HTML")
        # 后处理: 发表级三件套 (同单病毒模式修复: beast_log 多链探测 + 补 mcc_tree)
        _beast_dir_w = os.path.join(str(work_dir), "time", "beast")
        _beast_log_w = None
        for _cand in (os.path.join(_beast_dir_w, "merged", "merged.log"),
                      os.path.join(_beast_dir_w, "beast1.log"),
                      os.path.join(_beast_dir_w, "beast11.log")):
            if os.path.exists(_cand):
                _beast_log_w = _cand; break
        _mcc_tree_w = os.path.join(_beast_dir_w, "merged", "mcc.tree")
        # 同 --virus 分支: 后处理走治理产物优先, 而非未解析/未治理的 _resolve(args.meta)
        _pp_w = _governed_meta(prep, _resolve(args.meta))
        _pf_cmd_w = ["--rtt_fasta", _resolve(args.alignment),
                        "--rtt_dates", _pp_w["dates"],
                        "--rtt_tree", _resolve(args.tree)]
        if _beast_log_w:
            _pf_cmd_w += ["--beast_log", _beast_log_w]
        if os.path.exists(_mcc_tree_w):
            _pf_cmd_w += ["--mcc_tree", _mcc_tree_w]
        _pf_cmd_w += ["--output_dir", os.path.join(report_dir, "publication")]
        run_postprocess("publication_figures", logger, *_pf_cmd_w)
        # 后处理: 采样分布图
        run_postprocess("sample_plot", logger,
                        "--metadata", _pp_w["meta"],
                        "--outdir", os.path.join(report_dir, "sample_plot"))
        logger.info(f"\nDone. Results: {work_dir}")
        return

    # ── Batch mode ──
    if not args.extract_dir:
        logger.error("Need --extract_dir (batch mode) or --work_dir (single mode)")
        sys.exit(1)

    extract_dir = Path(args.extract_dir).resolve()
    if not extract_dir.exists():
        logger.error(f"Extract dir not found: {extract_dir}")
        sys.exit(1)

    # ── Stage prep: Phase 1 数据收集 ──
    logger.info(f"\n{'='*60}")
    logger.info("Stage prep: Data Collection")
    logger.info(f"{'='*60}")
    logger.info(f"  Extract dir: {extract_dir}")

    sample_meta = load_sample_metadata(
        info_dir=args.info_dir or str(extract_dir.parent),
        metadata_csv=args.metadata, logger=logger)
    logger.info(f"  Sample metadata: {len(sample_meta)} entries")

    viruses = scan_extraction_dir(str(extract_dir), logger)
    logger.info(f"  Viruses found: {len(viruses)}")

    virus_queue = []
    for vkey, vinfo in viruses.items():
        vout = os.path.join(str(output_dir), vkey)
        ref_fa = find_ref_fasta(vkey, args.variants_dir)
        prep = prepare_virus_inputs(vinfo, vout, sample_meta, ref_fa, args.min_samples, logger)
        if prep:
            virus_queue.append(prep)
    virus_queue.sort(key=lambda x: x["n_with_dates"], reverse=True)
    if args.max_viruses > 0:
        virus_queue = virus_queue[:args.max_viruses]

    logger.info(f"\n  Viruses queued: {len(virus_queue)}")
    for i, vq in enumerate(virus_queue):
        logger.info(f"  [{i+1}] {vq['virus_name']}: {vq['n_with_dates']}/{vq['n_samples']} samples with dates")

    if not virus_queue:
        logger.warning("No viruses meet the minimum sample threshold. Exiting.")
        sys.exit(0)

    # ── 每病毒按 stage 执行 ──
    all_results = []
    for i, vq in enumerate(virus_queue):
        logger.info(f"\n--- [{i+1}/{len(virus_queue)}] ---")
        t_start = time.time()
        results = process_virus(vq, stages, args, tools, logger)
        t_elapsed = time.time() - t_start
        status = "✓" if results.get("success") else "✗"
        logger.info(f"  [{i+1}/{len(virus_queue)}] {status} {vq['virus_name']} ({t_elapsed:.0f}s)")
        all_results.append(results)

    # ── 汇总 CSV + HTML ──
    # 说明 (P1-7): 'report' 在 STAGE_ORDER 里, 但汇总实际由这里(process_virus 之外)
    # 无条件产出 —— 也就是 report 从来不是一个"阶段分支", 只做收尾。
    # 之所以不能改成 `if 'report' in stages`, 是因为 --stage clock 这类局部重跑
    # 也必须刷新汇总表。开关用 config 的 report.enabled, 而不是 stage 清单。
    report_dir = os.path.join(str(output_dir), "report")
    os.makedirs(report_dir, exist_ok=True)
    if getattr(args, "report_enabled", True):
        generate_summary(all_results, report_dir, logger)
        if all_results:
            report_path = build_html_report(all_results, report_dir)
            logger.info(f"  HTML report: {report_path}")
    else:
        logger.info("  [report] report.enabled=false, 跳过 phylo_summary.csv / HTML")

    # ── 每病毒后处理: 采样分布图 + 发表级三件套 ──
    # 历史坑: sample_plot / publication_figures 曾只挂在单病毒与 work-dir 模式,
    # batch 模式 (主路径) 从未调用 → 采样时间线/地图与 MCC 树图静默缺失
    for vq in virus_queue:
        vname = vq["virus_name"]
        vdir = os.path.join(str(output_dir), vname)
        v_meta = vq.get("meta_csv") or vq.get("dates_csv", "")
        v_report = os.path.join(vdir, ph_dir('report'))
        # 采样分布图 (时间线 + 地理地图)
        if v_meta and os.path.exists(v_meta):
            run_postprocess("sample_plot", logger,
                            "--metadata", v_meta,
                            "--outdir", os.path.join(v_report, "sample_plot"))
        # 分基因 tMRCA 汇总 (森林图 + TSV; gene_dating 后台 BEAST 完成后才有数据)
        _gd_dir = os.path.join(vdir, ph_dir('time'), "gene_dating")
        if os.path.isdir(_gd_dir):
            try:
                from utils.stage_plots import plot_gene_dating_summary
                _tsv, _fig = plot_gene_dating_summary(_gd_dir, _gd_dir)
                if _fig:
                    logger.info(f"  [gene_tmrca] {vname}: {_fig}")
            except Exception as _e:
                logger.warning(f"  [gene_tmrca] {vname} skip: {_e}")
        # 发表级三件套 (beast_log 多链探测 + mcc_tree)
        _bd = os.path.join(vdir, ph_dir('time'), "beast")
        _bl = None
        for _cand in (os.path.join(_bd, "merged", "merged.log"),
                      os.path.join(_bd, "beast1.log"),
                      os.path.join(_bd, "beast11.log")):
            if os.path.exists(_cand):
                _bl = _cand; break
        _mcc = os.path.join(_bd, "merged", "mcc.tree")
        _cmd = ["--rtt_fasta", vq.get("combined_fasta", ""),
                "--rtt_dates", v_meta,
                "--rtt_tree", vq.get("tree_file", "")]
        if _bl:
            _cmd += ["--beast_log", _bl]
        if os.path.exists(_mcc):
            _cmd += ["--mcc_tree", _mcc]
        _cmd += ["--output_dir", os.path.join(v_report, "publication")]
        run_postprocess("publication_figures", logger, *_cmd)

    logger.info(f"\n{'='*60}")
    logger.info(f"PhyloPipeline complete!")
    logger.info(f"Output: {output_dir}")
    logger.info(f"{'='*60}")


if __name__ == "__main__":
    main()
