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
    'treetime',
    'phylogeo', 'geo', 'geo_paths',
    'phymapr',
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
    # 2026-10-12 新增 (STAGE_REFERENCE §19/§20): TreeTime 路线的两环。
    # treetime 只吃 tree (+比对/日期); phymapr 吃 treetime 的 timetree.nexus,
    # **不依赖 BEAST** —— 与 BEAST 路线 (phylogeo→geo_paths) 完全并列。
    'treetime':     ['tree'],
    'phymapr':      ['treetime'],
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
    'treetime': 'time',
    'phylogeo': 'geography', 'geo': 'geography', 'geo_paths': 'geography',
    'phymapr': 'geography',
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
    'time':      ['clock', 'beast', 'gene_dating', 'treetime'],
    'geography': ['phylogeo', 'geo', 'geo_paths', 'phymapr'],
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
    'treetime':     'TreeTime CLI 完整时间树 (timetree.nexus; phymapr 路线前置)',
    'phylogeo':     '系统地理 (CTMC+BSSVS)',
    'geo':          '地理分析 (Mantel/树地理)',
    'geo_paths':    '事件级传播路径 (pathways/episode/置信分级 + 扩散统计 + 传播图)',
    'phymapr':      'phymapr R BM-ASR 系统地理 + 传播地图 (TreeTime 路线)',
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
    # 2026-10-12 修 (P1): 旧值 'time/beast/beast1.xml' 是**输入 XML** —— 它由
    # generate_beast1_phylogeo_xml() 在 BEAST 开跑**之前**写出。于是只要 XML 一落盘,
    # stage_done('beast') 就为真, 哪怕 BEAST 一条链都没跑、merge 一次都没做:
    #   · run_phase_beast 的 `beast_pending` 分支注释写着"不写 checkpoint, 下次跑
    #     检查合并" —— 但下次 _done('beast') 依然为真 → 整段跳过 → merge_chains()
    #     永不执行 → time/beast/merged/ 永久缺失。(实测: GCVA 结果树 merged/ 缺失,
    #     只能手工 `python -m utils.merge_results --work-dir ...` 补。)
    #   · 与同类 'phylogeo' 的标志产物 (merged/mcc.tree) 口径必须一致: 标志产物
    #     是**产物**, 不是输入。
    # 取 merged/mcc.tree: 多链模式的终态产物 (logcombiner + treeannotator 之后)。
    # 单链模式没有 merged/, 首次重入会走 run_phase_beast 的"输出已存在, skip run"
    # 分支并 _mark('beast') 写下 .done —— 自限, 不会反复重跑。
    'beast':       'time/beast/merged/mcc.tree',
    'gene_dating': 'time/gene_dating/gcva_watch_parse.sh',
    # 2026-10-12 新增 (STAGE_REFERENCE §19/§20; TreeTime 路线两 stage)
    'treetime':    'time/treetime/timetree.nexus',
    'phymapr':     'geography/phymapr/phymapr_summary.json',
    'phylogeo':    'geography/phylogeography/merged/mcc.tree',
    'geo':         'geography/geo_analysis/mantel_test.pdf',
    'geo_paths':   'geography/geo_analysis/pathways/pathway_qc.json',  # 2026-09-16 新增
    'host':        'popgen/host_analysis/tree_host.pdf',
    'report':      'report/phylo_summary.csv',
    'prep':        None,   # 数据收集不参与 checkpoint
}


def _cp_file(virus_dir, stage):
    return os.path.join(virus_dir, CHECKPOINT_DIR, f"{stage}.done")


# STAGE_CHECKPOINT_FILES 里的相对路径第一段是**模块名**(legacy 口径)。
# standard (默认) 布局下模块目录带编号 (time→06_time, report→08_report …),
# 直接 join 会永远不命中 —— 探测时把第一段换算成 ph_dir(), 并与 legacy 字面量
# 一并尝试 (存量结果树兼容)。
_MARKER_MODULE_KEYS = ("data", "phylogeny", "popgen", "recomb", "select",
                       "time", "geography", "report")


def _marker_candidates(virus_dir, marker: str) -> List[str]:
    """checkpoint 标志产物 (legacy 相对字面量) → 候选绝对路径 (standard 名优先)。"""
    parts = marker.split("/")
    cands = []
    if parts[0] in _MARKER_MODULE_KEYS and len(parts) > 1:
        cands.append(os.path.join(virus_dir, ph_dir(parts[0]), *parts[1:]))
    legacy = os.path.join(virus_dir, marker)
    if legacy not in cands:
        cands.append(legacy)
    return cands


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
        # 2026-10-12 修 (P1-6): 表内路径是 legacy 字面量 (report/phylo_summary.csv),
        # 而 standard (默认) 布局的模块目录叫 08_report/… —— 旧实现直接 join
        # 永远不命中, "标志产物回退"这条判定路径在 standard 下形同虚设。
        # 现按当前布局换算第一段模块名, 并与 legacy 字面量一并探测 (存量树兼容)。
        for p in _marker_candidates(virus_dir, marker):
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


def _usable_alignment(path) -> bool:
    """比对文件是否"可用" —— 存在 + 非空 + 至少有一条 FASTA 记录。

    2026-10-11 新增 (修 P0)。此前全管线对输入只用 `os.path.exists()` 判定:
    PSTVd 工作副本 `pstvd_mafft.aln.fasta` 被截成 0 字节后, align 的 checkpoint
    分支照样把它当成合法比对 → rdp5/tree/clock/beast 全部吃空文件, 各自以
    "没有产出"的形式失败, 而**真正的原因(输入为空)从未被喊出来**。
    「0 字节不是合法输入」必须在唯一的入口处断言。
    """
    if not path:
        return False
    try:
        if not os.path.isfile(path) or os.path.getsize(path) == 0:
            return False
        with open(path, "rb") as f:
            head = f.read(4096)
    except OSError:
        return False
    # 必须是 FASTA: 首部出现 '>' 记录符
    return b">" in head


def _usable_file(path) -> bool:
    """存在且非空 (用于 tree/dates 等非 FASTA 输入)。"""
    if not path:
        return False
    try:
        return os.path.isfile(path) and os.path.getsize(path) > 0
    except OSError:
        return False


def _derive_virus_name(data_dir, ds) -> str:
    """从 dataset 条目推出**病毒名** (用于汇总表 Virus 列 / 产物文件名前缀)。

    2026-10-11 修 (P1)。原写法是
        os.path.basename(os.path.dirname(data_dir.rstrip('/\\\\')))
    即"病毒名 = data_dir 的**父目录名**" —— 这要求 datasets.yaml 的 dir 恰好是
    `<phase>/<virus>` 两层。实测 PSTVD_FULL 的 dir 是
        {base}/PSTVd/phylo_pstvd_full_20260831/Potato_spindle_tuber_viroid_NC_002030.1
    父目录名 = `phylo_pstvd_full_20260831` → 汇总表 Virus 列写成了**阶段目录名**,
    而 GCVA (dir 多一层 `gcva_analysis`) 恰好推对。同一张表里两个病毒的 Virus 列
    口径不同。

    新规则(按证据强度):
      ① 路径中**含 accession 的那一层** —— 这是病毒目录的命名契约
         (`<物种>_<accession>`), 与 dir 嵌套层数无关;
      ② 退而用 species+accession 拼 (与 ① 同构);
      ③ 再退: 从深到浅找**第一个非通用目录名** (跳过 `phylo_*` 阶段目录与
         `gcva_analysis`/`analysis_v2`/`data` 之类的工作目录名) —— 这一步修正了
         旧写法"取父目录名"的错误 (父目录恰好是阶段目录)。
    """
    acc = str(ds.get("accession") or "").strip()
    parts = [p for p in Path(data_dir).parts if p not in ("/", "\\")]
    if acc:
        for p in reversed(parts):
            if acc in p:
                return p
    sp = str(ds.get("species") or "").strip()
    if sp and acc:
        return f"{sp.replace(' ', '_')}_{acc}"
    # ③ 跳过通用/阶段目录名, 取最深的"像病毒名"的那一层
    for p in reversed(parts):
        if p in _GENERIC_DIR_NAMES or p.startswith("phylo_"):
            continue
        return p
    return parts[-1] if parts else "unknown_virus"


# 非病毒名的通用目录 (阶段目录 / 工作目录)。_derive_virus_name 的兜底规则用它。
_GENERIC_DIR_NAMES = {
    "gcva_analysis", "analysis", "analysis_v2", "analysis_v3", "phylo_results",
    "results", "data", "report", "time", "phylogeny", "popgen", "geography",
    "recomb", "select", "assembly", "extraction", "phylogeny_out",
}

# 各布局下"模块目录"的代表名, 用于从既有工作目录**反推**它用的是哪套布局。
_LAYOUT_MARKERS = {
    "standard": ("01_data", "08_report"),
    "legacy": ("data", "report"),
}


def detect_layout_from_tree(work_dir: str) -> Optional[str]:
    """从既有工作目录反推 I/O 布局; 判不出 (空目录 / 两套都在) 返回 None。

    2026-10-11 新增。布局解析优先级是 `--io-layout` > `MMPV_IO_LAYOUT` >
    **standard**(默认), 而存量结果树可能是 legacy。此时按默认跑会**静默新建一套
    并行目录树**, 老目录纹丝不动 —— 用户去看老报告发现"什么都没变"。
    实测 2026-10-11: 不带 --io-layout 跑 PSTVd 的 report, 直接在病毒目录下新建了
    `01_data/` 与 `08_report/`, 而全部存量结果都在 `data/` 与 `report/` 里。
    """
    if not work_dir or not os.path.isdir(work_dir):
        return None
    hits = [lay for lay, marks in _LAYOUT_MARKERS.items()
            if any(os.path.isdir(os.path.join(work_dir, m)) for m in marks)]
    return hits[0] if len(hits) == 1 else None


def clear_stage_checkpoint(virus_dir, stage=None) -> None:
    """清除 checkpoint (stage=None 清全部)"""
    cp_dir = os.path.join(virus_dir, CHECKPOINT_DIR)
    if stage:
        p = _cp_file(virus_dir, stage)
        if os.path.exists(p):
            os.remove(p)
        # 2026-10-11: 连同状态快照一起清 —— 否则清了 .done 却留着 <stage>.state.json,
        # 重跑若中途失败, 下游 _restore 会把**上一轮**的 outputs/metrics 当成
        # 本轮结果回填 (幽灵数据)。状态快照与 .done 必须同生共死。
        sp = _state_file(virus_dir, stage)
        if os.path.exists(sp):
            os.remove(sp)
        return
    if os.path.isdir(cp_dir):
        import shutil
        shutil.rmtree(cp_dir, ignore_errors=True)


# ═══════════════════════════════════════════════════════════════════
# Stage 状态快照 (2026-10-11 新增, 修 P0: "重跑即空汇总")
# ═══════════════════════════════════════════════════════════════════
# 【问题】每个 stage 的 `if _done(x): log(...); _skip(x)` 分支只打一行日志,
#   **不恢复该 stage 的 outputs/metrics**。于是任何「checkpoint 已存在, 只重跑
#   report / beast / gene_dating」的增量重跑里, results["metrics"] 只含**本次
#   请求的 stage** 写过的东西 → headline_metrics() 取不到任何速率/R²/p/TMRCA。
#
#   实测 (PSTVd, 2026-09-22 活副本 report/phylo_summary.csv):
#     Virus,Phase,Success,Stages_done,...,Beta,R_squared,P_value,TMRCA
#     phylo_pstvd_full_20260831,,True,clean;align;saturation,,,,,,
#   即整行只有 3 个 stage 名, 全部头条指标为空 —— 而磁盘上
#   time/clock_summary.tsv / time/temporal_signal/drt_results.json /
#   time/beast/merged/ 的结果一直都在。
#
# 【修法】stage 成功时把状态落盘:
#     <virus_dir>/.checkpoints/<stage>.state.json
#       {"stage":…, "time":…, "outputs":{…}, "metrics":{…}, "extra":{…}}
#   skip 分支统一调 _restore(stage) 把它读回来 → 汇总/报告与"当次真跑"完全等价。
#   extra 承载 stage 间的局部变量 (align→aln_file, tree→tree_file,
#   clean→combined_fasta, metadata→metadata_std, rdp5→aln_file)。
#
# 【为什么不直接从产物反推】产物反推 (如从 clock_summary.tsv 读 β/R²) 需要每个
#   stage 各写一套解析器, 且口径易漂移; 落盘快照是**同一份内存对象**的副本,
#   与当次运行的 results 逐键一致, 不会出现"报告与产物对不上"。
STAGE_STATE_SUFFIX = ".state.json"

# 整个 results 的快照 (供 main() 在"只跑部分 stage"时重建汇总)
STAGE_RESULTS_FILE = "stage_results.json"


def _state_file(virus_dir, stage):
    return os.path.join(virus_dir, CHECKPOINT_DIR, f"{stage}{STAGE_STATE_SUFFIX}")


def _jsonable(obj):
    """把任意对象压成 JSON 可序列化形式 (不可序列化的退化成 str)。

    2026-10-11: metrics 里混有 numpy 标量 / Path / set, 直接 json.dump 会抛
    TypeError 且**冒泡打断整个 stage** —— 状态落盘是"锦上添花"的能力,
    绝不允许它把主流程搞挂, 故先做一次显式降级。
    """
    if obj is None or isinstance(obj, (bool, int, float, str)):
        return obj
    if isinstance(obj, dict):
        return {str(k): _jsonable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_jsonable(v) for v in obj]
    if isinstance(obj, set):
        return sorted(str(v) for v in obj)
    try:
        import numpy as _np
        if isinstance(obj, _np.generic):
            return obj.item()
        if isinstance(obj, _np.ndarray):
            return _jsonable(obj.tolist())
    except Exception:
        pass
    return str(obj)


def save_stage_state(virus_dir, stage, outputs=None, metrics=None, extra=None) -> None:
    """落盘单 stage 状态快照 (失败仅告警, 绝不阻断主流程)。"""
    try:
        cp_dir = os.path.join(virus_dir, CHECKPOINT_DIR)
        os.makedirs(cp_dir, exist_ok=True)
        payload = {
            "stage": stage,
            "time": datetime.now().isoformat(timespec="seconds"),
            "outputs": _jsonable(outputs or {}),
            "metrics": _jsonable(metrics or {}),
            "extra": _jsonable(extra or {}),
        }
        with open(_state_file(virus_dir, stage), "w", encoding="utf-8") as f:
            json.dump(payload, f, ensure_ascii=False, indent=1)
    except Exception as e:                      # noqa: BLE001 — 兜底不阻断
        try:
            print(f"[checkpoint] ⚠ {stage} 状态快照落盘失败 (非阻断): {e}",
                  file=sys.stderr)
        except Exception:
            pass


def load_stage_state(virus_dir, stage) -> Dict:
    """读回单 stage 状态快照; 不存在/损坏返回 {}。"""
    p = _state_file(virus_dir, stage)
    if not os.path.exists(p):
        return {}
    try:
        with open(p, encoding="utf-8") as f:
            data = json.load(f)
        return data if isinstance(data, dict) else {}
    except Exception:
        return {}


def save_stage_results(virus_dir, results) -> None:
    """落盘整个 results (供 main() 在部分 stage 重跑时重建完整汇总)。"""
    try:
        os.makedirs(virus_dir, exist_ok=True)
        payload = {
            "time": datetime.now().isoformat(timespec="seconds"),
            "results": _jsonable(results),
        }
        with open(os.path.join(virus_dir, STAGE_RESULTS_FILE), "w",
                  encoding="utf-8") as f:
            json.dump(payload, f, ensure_ascii=False, indent=1)
    except Exception as e:                      # noqa: BLE001
        try:
            print(f"[checkpoint] ⚠ results 快照落盘失败 (非阻断): {e}",
                  file=sys.stderr)
        except Exception:
            pass


def load_stage_results(virus_dir) -> Dict:
    """读回整个 results 快照; 不存在/损坏返回 {}。"""
    p = os.path.join(virus_dir, STAGE_RESULTS_FILE)
    if not os.path.exists(p):
        return {}
    try:
        with open(p, encoding="utf-8") as f:
            data = json.load(f)
    except Exception:
        return {}
    res = data.get("results") if isinstance(data, dict) else None
    return res if isinstance(res, dict) else {}


def _stage_state_files(virus_dir) -> Dict[str, Dict]:
    """扫描 .checkpoints/*.state.json → {stage: payload} (按 stage 名排序)。"""
    out: Dict[str, Dict] = {}
    cp_dir = os.path.join(virus_dir, CHECKPOINT_DIR)
    if not os.path.isdir(cp_dir):
        return out
    for fn in sorted(os.listdir(cp_dir)):
        if not fn.endswith(STAGE_STATE_SUFFIX):
            continue
        stage = fn[: -len(STAGE_STATE_SUFFIX)]
        data = load_stage_state(virus_dir, stage)
        if data:
            out[stage] = data
    return out


def merge_persisted_results(results: Dict, virus_dir: str) -> Dict:
    """把磁盘持久化的 stage 状态并入当次 results (当次值优先, 磁盘只补缺)。

    用途: 「只跑 report」「只跑 beast」这类增量重跑里, process_virus 只记录了
    **本次请求的 stage** —— 汇总表/HTML 报告因此退化成只有本次 stage 的空壳。
    这里把 <virus_dir>/stage_results.json 与 .checkpoints/*.state.json 里的
    outputs/metrics/stage 清单回填, 使"重跑一次 report"能产出与全量跑等价的汇总。

    返回**新的** results dict (不改入参), 便于调用方显式接管。
    """
    if not virus_dir or not os.path.isdir(virus_dir):
        return results
    merged = dict(results or {})
    # ⚠ 必须**深一层拷贝** outputs/metrics —— `dict(results)` 是浅拷贝, 内层
    # 仍是同一个 dict 对象; 后面往 merged["metrics"][k] 写值会**改到调用方的
    # results**, 破坏"纯函数"契约 (本文件单测 [6] 的「merge 不改入参」断言
    # 就是专门守这一点的)。
    merged["outputs"] = dict((results or {}).get("outputs") or {})
    merged["metrics"] = dict((results or {}).get("metrics") or {})

    # ① 整份 results 快照 (上一轮全量/增量跑的最终态)
    prev = load_stage_results(virus_dir)

    # ② 逐 stage 快照 (比整份快照更细, 且能覆盖"整份快照还没写"的旧数据)
    per_stage = _stage_state_files(virus_dir)

    def _fill_from(src_outputs, src_metrics):
        """当次值优先; 只补磁盘上有、当次没有的键。"""
        for k, v in (src_outputs or {}).items():
            if k not in merged["outputs"] and v not in (None, "", {}):
                merged["outputs"][k] = v
        for k, v in (src_metrics or {}).items():
            # 指标缺值判定: 键不存在 **或** 值为 None/"" —— 与
            # summary_metrics._pick 的"空串/None 都视为缺"保持一致
            if merged["metrics"].get(k) in (None, "") and v not in (None, ""):
                merged["metrics"][k] = v

    # 先并入整份快照, 再用逐 stage 快照补齐(逐 stage 通常更完整)
    sources = []
    if prev:
        sources.append((prev.get("outputs"), prev.get("metrics")))
    for stage, payload in per_stage.items():
        sources.append((payload.get("outputs"), payload.get("metrics")))

    for so, sm in sources:
        _fill_from(so, sm)

    # 2026-10-11: 对"有 <stage>.done 但没有 state 快照"的 stage 补一次**产物重建**。
    # 【为什么必须在 merge 里做, 而不是只在 process_virus 的 skip 分支里做】
    #   ① 老结果树的 *.done **全部**没有 state 快照 (实测活副本 23/23), 而跳过分支
    #      只覆盖"本次请求的 stage" —— 跑 `--stage report` 时 align/clock/beast 的
    #      代码块根本不会进入, 汇总依旧是空的;
    #   ② merge 是三种模式 (--virus / --work_dir / batch) 汇总前的唯一公共入口。
    # 【优先级】本次运行写下的值 > *.state.json 快照 > 磁盘产物重建
    #   (_fill_from 只在当前值为空时才写, 故顺序天然满足)。
    # 重建成功后落盘快照 → 下次不必重复解析 (自愈)。
    try:
        _cp_entries = sorted(os.listdir(os.path.join(virus_dir, CHECKPOINT_DIR)))
    except OSError:
        _cp_entries = []
    for _fn in _cp_entries:
        if not _fn.endswith(".done"):
            continue
        _st = _fn[:-len(".done")]
        if _st in per_stage:          # 已有快照, 上面已并入
            continue
        _rd = rediscover_stage_state(virus_dir, _st)
        if not _rd:
            continue
        _fill_from(_rd.get("outputs"), _rd.get("metrics"))
        save_stage_state(virus_dir, _st, outputs=_rd.get("outputs"),
                         metrics=_rd.get("metrics"), extra=_rd.get("extra"))

    # stage 清单取并集 (当次在前, 历史在后, 保序去重)
    for field in ("stages_done", "stages_skipped", "stages_failed"):
        cur = [str(x) for x in (merged.get(field) or [])]
        for src in [prev] + list(per_stage.values()):
            if field not in src:
                continue
            for x in (src.get(field) or []):
                x = str(x)
                if x not in cur:
                    cur.append(x)
        merged[field] = cur
    # 显式 `<stage>.skipped` 标记 (内层已定性"不适用") —— 对 done 有**否决权**。
    _skipped_marks = {_fn[:-len(".skipped")] for _fn in _cp_entries
                      if _fn.endswith(".skipped")}

    # (a) 用 <stage>.done 的存在补齐 stages_done。
    # 2026-10-11 修: 原实现只遍历 `per_stage` (有 *.state.json 的那些) —— 老结果树
    # 一个快照都没有, 于是 Stages_done 只剩"本次请求的 stage"(PSTVd 实测只写了
    # `clean;align;saturation`), 报告里所有模块都被判成"未运行"。改为认**全部** .done。
    # 注意: 必须保留这条 —— 有些 stage 不在 STAGE_ORDER 里 (saturation / splitstree),
    # 只能靠 *.done 发现。
    for _fn in _cp_entries:
        if not _fn.endswith(".done"):
            continue
        _st = _fn[:-len(".done")]
        if _st in _skipped_marks:      # 不适用标记优先 (遗留 .done 可能是旧的虚报)
            continue
        if _st not in merged["stages_done"]:
            merged["stages_done"].append(_st)

    # (b) 2026-10-11 二次修: 只认 *.done 会**漏报** —— `stage_done()` 还有一条
    # 「STAGE_CHECKPOINT_FILES 声明的标志产物存在」的回退路径。实测 GCVA:
    # `time/beast/beast1.xml` 存在 → `stage_done('beast')` 为真 (于是 beast 整段被
    # checkpoint 跳过), 但 `beast.done` 从未写过 → 汇总表 Stages_done 里没有 beast,
    # 与管线自己的判断**自相矛盾**。这里补一次 stage_done() 判定, 与管线同源。
    for _st in STAGE_ORDER:
        if _st in ("prep", "report") or _st in _skipped_marks:
            continue
        if _st in merged["stages_done"]:
            continue
        try:
            if stage_done(virus_dir, _st):
                merged["stages_done"].append(_st)
        except Exception:
            pass

    # (c) 2026-10-11 三次修: `<stage>.skipped` → stages_skipped, 并从 done 里剔除。
    # 此前"不适用"只在**本次运行请求了该 stage** 时才进 st['skipped'] —— 而
    # `--stage report` 不进任何 stage 代码块, 于是 rdp5/capheine 的"不适用"在
    # 汇总表里**彻底消失** (既不在 done 也不在 skipped)。标记文件是权威来源。
    if _skipped_marks:
        merged["stages_done"] = [x for x in merged["stages_done"]
                                 if x not in _skipped_marks]
        for _st in sorted(_skipped_marks):
            if _st not in merged["stages_skipped"]:
                merged["stages_skipped"].append(_st)

    if prev.get("virus") and not merged.get("virus"):
        merged["virus"] = prev["virus"]
    if prev.get("virus_dir") and not merged.get("virus_dir"):
        merged["virus_dir"] = prev["virus_dir"]
    # success: 只要**没有失败 stage** 且**有实质产物**就算成功 (与 process_virus 同判据)
    merged["success"] = (not merged["stages_failed"]) and bool(
        merged["outputs"] or merged["metrics"] or merged["stages_done"])
    return merged


# ═══════════════════════════════════════════════════════════════════
# Stage 状态重建 (rediscover) —— 状态快照缺失时从磁盘产物反推
# ═══════════════════════════════════════════════════════════════════
# 【为什么还需要这一层】
#   <stage>.state.json 只在**新代码跑过的 stage** 里才会产生。以下情况它必然缺失,
#   而 skip 分支若只认快照, 汇总表照旧为空:
#     ① 结果树由**旧代码**产出 (PSTVd 2026-09-22 / GCVA 2026-09-02 都是);
#     ② 结果树被 `cp -a` 到别处 —— 只剩 *.done, 没有 *.state.json;
#     ③ 手工删过 state 文件。
#   实测: 活副本 results_PSTVd_GCVA_20260905 下 23 个 *.done, 带 state 快照的 **0 个**。
#
# 【设计原则】
#   · 只认**唯一权威产物** → 指标, 口径与"当次真跑"逐字一致 (同一文件 / 同一列 /
#     同一解析函数)。例: beast 的 TMRCA 一律取 parse_beast_log 的
#     treeModel.rootHeight.**median**; merged/parameter_summary.json 里的 tmrca.mean
#     是另一个估计量 —— 混用会让"跳过"与"真跑"对同一份数据给出不同数字。
#   · 不二次推断 / 不插值 / 不编造默认值: 产物缺失或解析失败就**留空**
#     (与 headline_metrics "空 = 没算出来" 的语义一致, 绝不用 0 冒充)。
#   · 重建成功即**落盘快照** → 下次直接读快照, 不再重复解析 (自愈)。
#
# 【覆盖范围为什么只到"指标 + 路径"】
#   HTML 报告的图件由 report_builder 的 MODULES 表按**目录 glob 扫盘**收集,
#   不读 results["outputs"]; 汇总 CSV 只读 results["metrics"]。故重建只需保证
#   ① metrics 齐 (clock/beast 是 headline_metrics 的唯一来源);
#   ② extra 齐 (align/clean/tree/rdp5/metadata 的产物路径 —— 下游 stage 靠它定位输入)。
def _rd_read_tsv(path) -> List[Dict[str, str]]:
    """读 TSV → [{列名: 值}]; 读不了返回 []。"""
    if not _usable_file(path):
        return []
    try:
        with open(path, encoding="utf-8-sig", newline="") as f:
            return [dict(r) for r in csv.DictReader(f, delimiter="\t")]
    except (OSError, UnicodeDecodeError, csv.Error):
        return []


def _rd_read_json(path) -> Dict:
    """读 JSON → dict; 读不了或顶层不是 dict 返回 {}。"""
    if not _usable_file(path):
        return {}
    try:
        with open(path, encoding="utf-8") as f:
            d = json.load(f)
    except (OSError, ValueError):
        return {}
    return d if isinstance(d, dict) else {}


def _rd_num(v) -> Optional[float]:
    """→ float; 空 / NA / 布尔 / 非有限 → None (区分"没算出来"与"算出 0")。"""
    if v is None or v == "" or isinstance(v, bool):
        return None
    try:
        f = float(v)
    except (TypeError, ValueError):
        return None
    return f if f == f and f not in (float("inf"), float("-inf")) else None


def _rd_bool(v) -> Optional[bool]:
    """'True'/'False'/'1'/'0'/… → bool; 认不出 → None。"""
    if isinstance(v, bool):
        return v
    s = str(v or "").strip().lower()
    if s in ("true", "1", "yes", "y", "t"):
        return True
    if s in ("false", "0", "no", "n", "f"):
        return False
    return None


def _rd_beast_burnin_pct(virus_dir) -> float:
    """从 run_status.json 反推本次 BEAST 的 burnin 百分比 (与真跑同源)。

    parse_beast_log 的 burnin_pct 直接决定 TMRCA 的取样区间 —— 重建时若硬编码
    10% 而实跑是别的值, 两边数字就对不上。run_status.json 记着 merge_burnin 与
    chain_length, 相除即真值; 文件缺失时退回 10.0 (config/argparse 默认)。
    """
    st = _rd_read_json(os.path.join(virus_dir, ph_dir('time'), "beast", "run_status.json"))
    b, n = _rd_num(st.get("merge_burnin")), _rd_num(st.get("chain_length"))
    if b is not None and n:
        return 100.0 * b / n
    return 10.0


def rediscover_stage_state(virus_dir: str, stage: str) -> Dict:
    """从磁盘产物重建单 stage 的 {outputs, metrics, extra}; 无法重建返回 {}。

    仅在 <stage>.state.json 缺失时被 _restore() 调用 (见 process_virus)。
    返回空 dict 表示"该 stage 的产物不足以反推" —— 调用方照常告警。
    """
    if not virus_dir or not os.path.isdir(virus_dir):
        return {}
    outs: Dict = {}
    mets: Dict = {}
    extra: Dict = {}

    def _reg(key, *parts, fasta=False):
        """产物存在且可用 → 登记 outputs 并返回路径; 否则 None。"""
        fp = os.path.join(virus_dir, *parts)
        if (_usable_alignment(fp) if fasta else _usable_file(fp)):
            outs[key] = fp
            return fp
        return None

    if stage == "metadata":
        fp = _reg("metadata_std", ph_dir('data'), "metadata", "metadata_std.csv")
        if fp:
            extra["metadata_std"] = fp

    elif stage == "clean":
        fp = (_reg("clean_fasta", ph_dir('data'), "clean", "clean.fasta", fasta=True)
              or _reg("clean_fasta", ph_dir('phylogeny'), "qc", "clean.fasta", fasta=True))
        if fp:
            extra["clean_fasta"] = fp

    elif stage == "align":
        # 真跑写 phylogeny/mafft.aln.fasta。但 work-dir 模式 / clean 已产出时,
        # run_stage_align 会**直接沿用 provided 而不落盘**该文件 (实测 PSTVd 活副本
        # 根目录只有 pstvd_mafft.aln.fasta, phylogeny/ 下无 mafft.aln.fasta)
        # → 依次回退到 clean 产物, 这正是下游 stage 实际吃的那份。
        fp = (_reg("alignment", ph_dir('phylogeny'), "mafft.aln.fasta", fasta=True)
              or _reg("alignment", ph_dir('data'), "clean", "clean.fasta", fasta=True)
              or _reg("alignment", ph_dir('phylogeny'), "qc", "clean.fasta", fasta=True))
        if fp:
            extra["aln_file"] = fp

    elif stage == "rdp5":
        fp = (_reg("masked_alignment", ph_dir('recomb'), "rdp5",
                   "rdp5_dropped.realn.fasta", fasta=True)
              or _reg("masked_alignment", ph_dir('recomb'), "rdp5", "masked",
                      "masked_N.fasta", fasta=True))
        if fp:
            extra["aln_file"] = fp
        _reg("rdp5_report", ph_dir('recomb'), "rdp5", "recombination_report.md")

    elif stage == "tree":
        fp = (_reg("tree", ph_dir('phylogeny'), "iqtree.treefile")
              or _reg("tree", ph_dir('phylogeny'), "tree.nwk"))
        if fp:
            extra["tree_file"] = fp

    elif stage == "popgen":
        _reg("popgen_report", ph_dir('popgen'), "popgen_report.txt")
        _reg("popgen_plot", ph_dir('popgen'), "popgen_metrics.png")

    elif stage == "capheine":
        _reg("capheine_sites", ph_dir('select'), "capheine", "drhip", "combined_sites.csv")

    elif stage == "host":
        _reg("host_diff_plot", ph_dir('popgen'), "host_analysis", "host_diff.png")

    elif stage == "clock":
        t = ph_dir('time')
        # ① clock_summary.tsv —— β/R²/P 与 drt_passed 的权威落盘 (clock_analysis 写)
        cs = _reg("clock_summary", t, "clock_summary.tsv")
        if cs:
            extra["clock_summary"] = cs
            g = next((r for r in _rd_read_tsv(cs)
                      if (r.get("level") or "").strip() == "genome"), None)
            if g:
                for k in ("beta", "r_squared", "p_value"):
                    v = _rd_num(g.get(k))
                    if v is not None:
                        mets[k] = v
                b = _rd_bool(g.get("drt_passed"))
                if b is not None:
                    mets["drt_passed"] = b
        # ② PIC 严格钟 (clock_analysis 写 robust_report.json 的 pic_strict_clock)
        rr = _reg("robust_report", t, "treetime_rtt", "robust", "robust_report.json")
        if rr:
            pic = _rd_read_json(rr).get("pic_strict_clock")
            if isinstance(pic, dict):
                for src, dst in (("rate", "pic_rate"), ("r_pic", "pic_r")):
                    v = _rd_num(pic.get(src))
                    if v is not None:
                        mets[dst] = v
        # ③ DRT 结论串 (与 metrics['temporal_signal'] 同一字段)
        drt = _reg("drt_results", t, "temporal_signal", "drt_results.json")
        if drt:
            d = _rd_read_json(drt)
            if d.get("conclusion"):
                mets["temporal_signal"] = d["conclusion"]
            v = _rd_num(d.get("r2_percentile"))
            if v is not None:
                mets["drt_r2_percentile"] = v
        # ④ treedater 四功能 (2026-10-11 起 virphy_bridge 落盘 metrics.json;
        #    旧数据没有该文件 → 留空, 不追认)
        td = _reg("treedater_metrics", t, "treedater_ltt", "metrics.json")
        if td:
            d = _rd_read_json(td)
            for src, dst in (("cov_rate", "td_cov_rate"),
                             ("relaxed_clock_p", "td_relaxed_clock_p"),
                             ("tmrca", "td_tmrca"),
                             ("mean_rate", "td_mean_rate"),
                             ("n_outlier_tips", "td_n_outlier_tips")):
                v = _rd_num(d.get(src))
                if v is not None:
                    mets[dst] = v
        # ⑤ 图件 (报告靠扫盘收集, 这里只登记便于下游引用)
        _reg("dated_phylo", t, "treedater_ltt", "Phylogeny_dated.pdf")
        _reg("ltt_plot", t, "treedater_ltt", "LTT.pdf")
        _reg("treedater_outlier_tips", t, "treedater_ltt", "outlier_tips.tsv")
        _reg("rtt_plot", t, "treetime_rtt", "RootToTip_regression.pdf")
        _reg("timetree", t, "treetime_rtt", "timetree_inferred.nwk")
        _reg("robust_rtt_plot", t, "treetime_rtt", "robust", "robust_rtt_compare.png")

    elif stage == "beast":
        bd = os.path.join(virus_dir, ph_dir('time'), "beast")
        _reg("beast_xml", ph_dir('time'), "beast", "beast1.xml")
        _reg("beast_mcc_tree", ph_dir('time'), "beast", "merged", "mcc.tree")
        _reg("beast_parameter_summary", ph_dir('time'), "beast", "merged",
             "parameter_summary.json")
        _reg("beast_ess", ph_dir('time'), "beast", "postprocess", "ess_results.csv")
        # TMRCA: 复用 parse_beast_log (与真跑、与 skip 分支同一入口同一口径)
        log_path = None
        for cand in (os.path.join(bd, "merged", "merged.log"),
                     os.path.join(bd, "beast1.log"),
                     os.path.join(bd, "beast11.log")):
            if _usable_file(cand):
                log_path = cand
                break
        if log_path:
            try:
                parsed = parse_beast_log(log_path,
                                         burnin_pct=_rd_beast_burnin_pct(virus_dir),
                                         log=LogCollector())
                kp = parsed.get("key_params") if isinstance(parsed, dict) else None
                rh = (kp or {}).get("treeModel.rootHeight") if isinstance(kp, dict) else None
                if isinstance(rh, dict) and _rd_num(rh.get("median")) is not None:
                    mets["beast_tmrca"] = rh["median"]
                    mets["beast_tmrca_hpd"] = (f"{rh['hpd_95_lower']:.1f}"
                                               f"–{rh['hpd_95_upper']:.1f}")
                    mets["beast_tmrca_ess"] = rh["ess"]
                if isinstance(kp, dict):
                    for k, dst in (("ucld.mean", "beast_ucld_mean"),
                                   ("ucld.stdev", "beast_ucld_stdev")):
                        v = kp.get(k)
                        if isinstance(v, dict) and _rd_num(v.get("mean")) is not None:
                            mets[dst] = v["mean"]
            except Exception:      # noqa: BLE001 —— 重建失败不得影响主流程
                pass

    elif stage == "gene_dating":
        _reg("gene_dating_dir", ph_dir('time'), "gene_dating")

    elif stage == "treetime":
        # 2026-10-12 (TreeTime 路线, STAGE_REFERENCE §19): checkpoint 标志物 = timetree.nexus
        tt = _reg("treetime_nexus", ph_dir('time'), "treetime", "timetree.nexus")
        if tt:
            extra["treetime_nexus"] = tt

    elif stage == "phymapr":
        # 2026-10-12 (STAGE_REFERENCE §20): 标志物 = phymapr_summary.json (桥落盘)
        pj = _reg("phymapr_summary", ph_dir('geography'), "phymapr", "phymapr_summary.json")
        if pj:
            d = _rd_read_json(pj)
            _m = d.get("metrics") or {}
            for k in ("n_pathways", "n_episodes", "n_locations",
                      "n_tips_used", "n_tips_total", "n_tips_dropped"):
                v = _rd_num(_m.get(k))
                if v is not None:
                    mets[f"phymapr_{k}"] = v
            for k, v in (_m.get("transmission_counts") or {}).items():
                mets[f"phymapr_conf_{k}"] = v
        for _k, _rel in (("tree_plot", "tree_plot.png"),
                         ("gif", "transmission_map.gif"),
                         ("static_map", "transmission_map_static.png"),
                         ("episodes_csv", "episodes.csv"),
                         ("pathways_csv", "pathways.csv")):
            _reg(_k, ph_dir('geography'), "phymapr", _rel)

    elif stage == "phylogeo":
        _reg("phylogeo_xml", ph_dir('geography'), "phylogeography", "merged", "mcc.tree")

    elif stage == "geo":
        _reg("mantel_plot", ph_dir('geography'), "geo_analysis", "mantel_test.pdf")
        _reg("geo_tree", ph_dir('geography'), "geo_analysis", "tree_geography.pdf")

    elif stage == "geo_paths":
        g = ph_dir('geography')
        _reg("pathway_qc", g, "geo_analysis", "pathways", "pathway_qc.json")
        _reg("pathways_csv", g, "geo_analysis", "pathways", "pathways.csv")
        _reg("dispersion_summary", g, "geo_analysis", "dispersion", "dispersion_summary.json")
        _reg("timescale_consistency", g, "geo_analysis", "timescale_consistency.json")
        _reg("markov_jumps", g, "geo_analysis", "markov_jumps.json")

    elif stage == "splitstree":
        _reg("splitstree_nexus", ph_dir('phylogeny'), "splitstree", "splits.nex")
        _reg("splitstree_pdf", ph_dir('phylogeny'), "splitstree", "neighbornet.pdf")
        _reg("splitstree_png", ph_dir('phylogeny'), "splitstree", "neighbornet.png")

    elif stage == "online":
        _reg("ncbi_fasta", ph_dir('data'), "ncbi_ref", "ncbi_ref.fasta")

    if not (outs or mets or extra):
        return {}
    return {"outputs": outs, "metrics": mets, "extra": extra}


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
    # 2026-10-12 新增 (TreeTime 路线; STAGE_REFERENCE §19/§20)。键必须与
    # argparse 的 --treetime_*/--phymapr_* 同步登记 —— (chains 键静默丢失的教训)。
    'treetime': {'enabled': 'treetime', 'clock_filter': 'treetime_clock_filter',
                 'coalescent': 'treetime_coalescent', 'clock_rate': 'treetime_clock_rate'},
    'phymapr': {'enabled': 'phymapr', 'map_style': 'phymapr_map_style',
                'n_prune_early': 'phymapr_n_prune_early',
                'direct_snp': 'phymapr_direct_snp',
                'indirect_snp': 'phymapr_indirect_snp',
                'distant_style': 'phymapr_distant_style',
                'pace': 'phymapr_pace'},
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
    # 2026-10-12: splitstree 段登记 (group_col 控制外圈分组弧线, RUN_GUIDE 已承诺)
    'splitstree': {'max_taxa': 'splitstree_max_taxa',
                   'group_col': 'splitstree_group_col'},
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

# RDP5 窗口方法的最小比对长度 (bp)。低于此值 RDP5 无足够位点, 静默跑不出 CSV
# —— 类病毒 (PSTVd 379 bp) 就属于这种情形, 必须记"不适用"而非"失败"。
# 2026-10-11: 原先这个 500 是 run_stage_rdp5 里的魔法数字, 现提为常量并写进
# skip_reason, 便于报告与日志引用同一口径。
RDP5_MIN_ALN_LEN = 500
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
        f"metadata governance report — {vname}",
        f"source file: {src}",
        f"rows: {s['n_total']}  kept: {s['kept']}",
        f"detected columns: date={s['date_col']!r}  loc={s['loc_col']!r}  name={s['name_col']!r}",
        f"coordinate columns: {s['coord_cols']}",
        f"date convention mode: {mode}  ({'mid-month/mid-year' if mode == 'mid' else 'start of month/year (VirPhyKit)'})",
        "",
        "date precision distribution:",
    ]
    for k, v in s["date"].items():
        lines.append(f"  {k:<8} {v}")
    if s["date_issues"]:
        lines.append(f"date issues: {s['date_issues']}")
    if s["geo_issues"]:
        lines.append(f"geo issues: {s['geo_issues']}")
    lines.append(f"dropped Unknown levels: {s['geo_drop_unknown']} rows / {s['geo_dropped_levels']} levels")
    if s["coords_derived"]:
        lines.append(f"coordinates derived: {s['coords_derived']}  source={s['coords_derived_by_source']}")
    if s["coords_unrecoverable"]:
        lines.append(f"coordinates not recoverable (blanked): {s['coords_unrecoverable']}")
    if s["coord_placeholder_fields"]:
        lines.append(f"coordinate-column placeholders (still empty after deriving): {s['coord_placeholder_fields']}")
    if s["blanked_fields"]:
        lines.append(f"blanked fields (placeholder → empty): {s['blanked_fields']}")
        lines.append(f"  rows with placeholder fields: {s['n_rows_with_placeholder']}")
    lines.append(f"already decimal years in column (not re-converted): {s['already_decimal_in_col']}")
    lines.append(f"dropped: placeholder rows {s['dropped_placeholder']};  "
                 f"no date {s['dropped_no_date']};  no location {s['dropped_no_location']}")
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
    # 2026-10-11 (P0): 复用条件由 `size > 0` 收紧为 `_usable_alignment`
    # (存在 + 非空 + 含 FASTA 记录符)。0 字节/半截文件不得被当成"已有比对"。
    if not force and _usable_alignment(aln_file):
        log.info(f"  [align] {vname}: MAFFT alignment exists, use")
        return {"success": True, "alignment": aln_file}
    if force and os.path.exists(aln_file):
        log.info(f"  [align] {vname}: --force, 重新计算比对 (忽略已有 {Path(aln_file).name})")

    provided = prep.get("combined_fasta", "")
    # 2026-10-11 修 (P0): 原判定只有 `os.path.exists(provided)` —— **0 字节的
    # "provided alignment" 会被当作成功返回**, 于是下游 rdp5/tree/clock 全部吃
    # 空文件。PSTVd 工作副本的 pstvd_mafft.aln.fasta 被截成 0 字节后, 整条链
    # 就是这样一路"成功"到各个 stage 各自失败, 真正的原因从未被喊出来。
    # 2026-10-13 修 (P0, AUDIT-1009 P0-1): 「provided 早退」只允许 prep 显式声明
    # provided_alignment=True 的模式 (work-dir 的 --alignment)。batch/datasets
    # 模式的 combined_fasta 已被 clean 就地改写成**未比对** clean.fasta —— 它非空
    # 且是合法 FASTA, 旧判定恒命中早退 → MAFFT 从不执行, rdp5/tree/clock/beast/
    # gene_dating/splitstree 全链吃未比对输入。
    if (provided and prep.get("provided_alignment")
            and provided != aln_file and _usable_alignment(provided)):
        log.info(f"  [align] {vname}: using provided alignment: {Path(provided).name}")
        return {"success": True, "alignment": provided}
    if (provided and prep.get("provided_alignment")
            and provided != aln_file and os.path.exists(provided)):
        # 存在但不可用 → 明确报错, 不要静默降级
        log.error(f"  [align] {vname}: provided alignment 不可用 (0 字节/非 FASTA): {provided}")
        return {"success": False,
                "error": f"provided alignment unusable (empty or not FASTA): {provided}"}

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


def run_stage_rdp5(aln_file, prep, tools, args, logger) -> Dict:
    """rdp5: RDP5 重组检测 → 删除重组子 (默认) 或 mask 区间 (--rdp5_mask)

    默认策略 (2026-08-27 改): drop_recombinants 删整条重组序列
      → 对剩余序列重新 MAFFT → 后续 tree/clock/beast 等全部吃重比对
    保留策略: --rdp5_mask 时走旧 mask 区间置 N (不重比对, 保样本量)

    ── 返回值 (2026-10-11 改为结构化) ──
    旧签名 `-> Optional[str]` 只有"真值/None"两态, 调用方 `if masked:` 无法区分
    「产出了干净比对」与「方法不适用, 原样返回输入」—— 于是 PSTVd 的短序列 skip
    被记成了 rdp5 **done**(写下 rdp5.done + 把原比对登记为 masked_alignment)。
    现返回三态:
        {"status": "ok",      "alignment": <干净比对路径>, "metrics": {...}}
        {"status": "skipped", "skip_reason": <原因>}
        {"status": "failed",  "error": <原因>}
    调用方据此分别 _mark / _mark_skipped / _fail。
    """
    vname = prep["virus_name"]
    out_dir = prep["output_dir"]
    log = logger
    rdp5_dir = os.path.join(out_dir, ph_dir('recomb'), "rdp5")
    os.makedirs(rdp5_dir, exist_ok=True)

    # ── 输入合法性: 空/损坏比对必须显式失败 ──
    # 2026-10-11 修 (P0): 原实现里 _aln_len==0 时 `0 < _aln_len < 500` 为假 →
    # 会**带着空文件去跑 RDP5**, 结果只是"CSV not found"这种误导性报错。
    if not _usable_alignment(aln_file):
        return {"status": "failed",
                "error": f"输入比对不可用 (缺失/0 字节/非 FASTA): {aln_file}"}

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
        # 2026-10-11 修 (P1): 原为 `max(len(s) for s in _first_seq)` —— 那是**最大行宽**,
        # 只在"一条序列一行"(MAFFT/iqtree 输出) 时等于序列长。若比对被换行(如 60 列),
        # 500 bp 的序列会被算成 60 → 误判为"类病毒超短序列"而 skip 整个 RDP5。
        # 改为**求和**, 对未换行文件结果不变, 对换行文件才正确。
        _aln_len = sum(len(s) for s in _first_seq)
    except OSError as e:
        return {"status": "failed", "error": f"读取比对失败: {e}"}
    if _aln_len == 0:
        return {"status": "failed", "error": "比对序列长度为 0 (空记录)"}
    if _aln_len < RDP5_MIN_ALN_LEN:
        reason = (f"alignment only {_aln_len} bp (<{RDP5_MIN_ALN_LEN}, viroid/very short sequence), "
                  f"RDP5 window methods lack enough sites, not applicable")
        log.info(f"  [rdp5] {reason}")
        return {"status": "skipped", "skip_reason": reason,
                "metrics": {"rdp5_aln_len": _aln_len, "rdp5_skipped": True}}

    log.info(f"  [rdp5] {vname}: RDP5 recombination detection (stage={args.rdp5_stage})...")

    rdp5_script = os.path.expanduser(args.rdp5_script)
    # 2026-10-12 (test E): 脚本缺失守卫 (直调也生效) —— 环境无 RDP5 时**明确报错并
    # 返回 None**, 绝不启动子进程把误导性的子进程报错当成本 stage 的结论。
    # 调用方 (process_virus) 预检后仍会归一化为 _skip('rdp5', 'na')。
    if not os.path.exists(rdp5_script):
        log.warning(f"  [rdp5] RDP5 脚本不存在: {rdp5_script} —— 环境无 RDP5, 跳过")
        return None
    rec_cmd = [sys.executable, "-m", "recombination_analysis",
               "--fasta", aln_file, "--prefix", vname, "--outdir", rdp5_dir,
               "--rdp5-script", rdp5_script, "--stage", args.rdp5_stage,
               "--threads", str(args.threads)]
    if args.rdp5_genes:
        rec_cmd += ["--genes", args.rdp5_genes]
    rr = subprocess.run(rec_cmd, capture_output=True, text=True, encoding='utf-8', errors='replace', cwd=str(SCRIPT_DIR))
    if rr.returncode != 0:
        msg = f"RDP5 failed (exit={rr.returncode}): {(rr.stderr or rr.stdout)[-400:]}"
        log.warning(f"  [rdp5] {msg}")
        return {"status": "failed", "error": msg}

    rdp5_csv = os.path.join(rdp5_dir, f"{vname}.csv")
    if not os.path.exists(rdp5_csv):
        msg = f"RDP5 CSV not found: {rdp5_csv}"
        log.warning(f"  [rdp5] {msg}")
        return {"status": "failed", "error": msg}
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
            msg = f"drop_recombinants failed: {(dr.stderr or dr.stdout)[-300:]}"
            log.warning(f"  [rdp5] {msg}")
            return {"status": "failed", "error": msg}
        dropped_fa = os.path.join(drop_dir, "dropped_recombinants.fasta")
        if not _usable_alignment(dropped_fa):
            msg = f"drop 产物缺失或为空: {dropped_fa}"
            log.warning(f"  [rdp5] {msg}")
            return {"status": "failed", "error": msg}
        # 删除后为纯序列 (含旧 gap 列), 重新 MAFFT 建立新比对坐标
        realn = os.path.join(rdp5_dir, "rdp5_dropped.realn.fasta")
        log.info(f"  [rdp5] realigning {Path(dropped_fa).name} (MAFFT) → {Path(realn).name}...")
        if not run_mafft(dropped_fa, realn, threads=args.threads, logger=log) \
                or not _usable_alignment(realn):
            msg = "realign 失败或产出为空 (fallback 到 masked 策略亦不可用)"
            log.warning(f"  [rdp5] {msg}")
            return {"status": "failed", "error": msg}
        log.info(f"  [rdp5] recombinant-dropped + realigned: {realn}")
        n_events = _count_rdp5_events(rdp5_csv)
        return {"status": "ok", "alignment": realn,
                "metrics": {"rdp5_events": n_events,
                            "rdp5_aln_len": _aln_len,
                            "rdp5_strategy": "drop_realign"}}

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
        msg = f"mask_recombination failed: {(mr.stderr or mr.stdout)[-300:]}"
        log.warning(f"  [rdp5] {msg}")
        return {"status": "failed", "error": msg}
    masked_fa = os.path.join(rdp5_dir, "masked", "masked_N.fasta")
    if _usable_alignment(masked_fa):
        log.info(f"  [rdp5] masked alignment: {masked_fa}")
        return {"status": "ok", "alignment": masked_fa,
                "metrics": {"rdp5_events": _count_rdp5_events(rdp5_csv),
                            "rdp5_aln_len": _aln_len,
                            "rdp5_strategy": "mask"}}
    msg = f"mask 产物缺失或为空: {masked_fa}"
    log.warning(f"  [rdp5] {msg}")
    return {"status": "failed", "error": msg}


def _count_rdp5_events(csv_path) -> int:
    """统计 RDP5 事件 CSV 的数据行数 (去掉表头); 解析失败返回 0。"""
    try:
        with open(csv_path, encoding="utf-8-sig", errors="replace") as f:
            lines = [ln for ln in f.read().splitlines() if ln.strip()]
        return max(0, len(lines) - 1)
    except OSError:
        return 0


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


def _trimal_trim(aln_file, out_fasta, logger, mode="gt", value=0.1):
    """用 trimAl 修剪比对 (剔除高 gap 列 / 低质量区), 供 NeighborNet 建网。

    ⚠ 为什么要做 (2026-10-12 第 11 轮) —— 喂给 Splitstree 的是**未修剪**的 mafft 输出,
      145x16310 里有 **9.5%** 的列 gap 率 >10% (421 个碎片段, 含末端残缺/引物区)。
      这些列的 P-distance 噪声会把 NeighborNet 的冲突信号堆在中心 → 出一团解不开
      的乱网。修剪后 (保留 90.9% 列) 网络明显收敛, 且 **aspect 从随机二选一
      稳定为 1.63** (不再出现 449/514 两种实现) —— 因为高变列正是 `EqualAngle` 的
      `isCircular` 判定翻转的来源。

    参数:
      mode='gt'   → `-gt <value>`  (gap 阈值; value=0.1 保留 gap<=10% 的列)
      mode='auto' → `-automated1`  (启发式; 本数据实测**过度修剪**, 退化成三叉星芒, 勿用)
      mode='nogaps' → `-nogaps`    (剪得最狠, 也偏退化)

    实测对照 (145x16310, P-distance + NeighborNet):
      | 参数        | 保留列 | 节点 | 边   | aspect |
      | 原始        | 16310  |449/514|625/829|1.00/1.62|
      | -gt 0.1     | 14829  | 516  | 835  | 1.63  |  ← 采用
      | -gt 0.2     | 14813  | 508  | 820  | 1.63  |
      | -automated1 | 14812  | 437  | 603  | 0.74  |  ← 退化
      | -nogaps     | 13596  | 460  | 650  | 1.03  |  ← 退化

    返回修剪后的 fasta 路径; 失败 (无 trimAl / 报错 / 结果为空) 返回 None (调用方回退原比对)。
    """
    try:
        trimal = find_tool("trimal")
    except Exception:                        # noqa: BLE001
        trimal = None
    if not trimal:
        import shutil as _sh
        trimal = _sh.which("trimal")
    if not trimal:
        logger.warning("  [splitstree] trimAl 不可用, 用未修剪比对建网")
        return None
    if mode == "auto":
        opt = ["-automated1"]
    elif mode == "nogaps":
        opt = ["-nogaps"]
    else:
        opt = ["-gt", str(value)]
    os.makedirs(os.path.dirname(os.path.abspath(out_fasta)), exist_ok=True)
    cmd = [trimal, "-in", aln_file] + opt + ["-fasta", "-out", out_fasta]
    try:
        env = dict(os.environ)
        # trimAl 动态链接 libgomp 等; 服务器上需 mambaforge/lib (见 MEMORY.md 环境节)
        _ml = os.path.expanduser("~/mambaforge/lib")
        if os.path.isdir(_ml):
            env["LD_LIBRARY_PATH"] = _ml + ":" + env.get("LD_LIBRARY_PATH", "")
        r = subprocess.run(cmd, capture_output=True, timeout=1800, env=env)
        if r.returncode != 0 or not os.path.exists(out_fasta):
            logger.warning("  [splitstree] trimAl 失败 (rc=%s): %s"
                           % (r.returncode, r.stderr.decode(errors="replace")[:200]))
            return None
        if os.path.getsize(out_fasta) < 100:
            logger.warning("  [splitstree] trimAl 输出为空, 回退原比对")
            return None
        # 记一行修剪摘要 (列数变化)
        try:
            def _ncol(p):
                _L = [l for l in open(p, errors="ignore") if not l.startswith(">")]
                return len(_L[0].strip()) if _L else 0
            logger.info("  [splitstree] trimAl %s: %d -> %d 列"
                        % (" ".join(opt), _ncol(aln_file), _ncol(out_fasta)))
        except Exception:                     # noqa: BLE001
            pass
        return out_fasta
    except Exception as e:                    # noqa: BLE001
        logger.warning("  [splitstree] trimAl 异常: %s" % e)
        return None


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
    # ── trimAl 修剪 (第 11 轮新增): 未修剪比对的高 gap 列会让 NeighborNet 中心
    #    堆成一团乱网, 且 aspect 在 449/514 两种实现间随机跳。修剪后收敛。
    #    `--splitstree-trim {gt,auto,nogaps,none}` 控制, 默认 gt=0.1。
    _tm = str(getattr(args, "splitstree_trim", "gt") or "gt").lower()
    _aln_use = aln_file
    if _tm != "none":
        _tf = os.path.join(outdir, "neighbornet_trimmed.aln.fasta")
        _r = _trimal_trim(aln_file, _tf, logger, mode=_tm,
                          value=float(getattr(args, "splitstree_trim_gap", 0.1)))
        if _r:
            _aln_use = _r
            results["outputs"]["splitstree_trimmed_aln"] = _r
    logger.info(f"  [splitstree] {vname}: NeighborNet 分裂网络...")
    import contextlib, io as _io
    argv_backup = sys.argv
    sys.argv = ['splitstree_bridge.py', '-i', _aln_use, '-o', outdir,
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
        # ── 出版级渲染 (2026-10-12 第 8 轮重做, 第 9 轮加分组弧线): 用 **SplitsTree6
        #    自己的 SplitNetworkLayout** 出图。旧版是 Python 自绘 equal-angle, 用户实测
        #    "根本不是 splitstree 的结果图" (自绘出星芒扇, SplitsTree 出宽扁网状)。
        #    坐标走 webserver /draw_distances (见 splitstree_bridge.native_coords;
        #    同一请求会随机返回两种实现, 第 9 轮起按 `select='radial'` 取**辐射版**)。
        #    外圈分组弧线 + 组名由 splitstree_native_viz.draw_group_arcs() 自绘 ——
        #    SplitsTree6 本体**没有**这个功能 (第 9 轮源码审查, 见 FIX_REPORT §⑬)。
        #    取不到坐标时回退自绘版。渲染失败不影响 nexus 产物。
        try:
            from utils.splitstree_viz import load_group_map
            _gc = str(getattr(args, "splitstree_group_col", "auto") or "auto")
            _groups, _gnote = (None, None)
            if _gc != "none":
                # ⚠ 2026-10-13: 分组列可由 metadata_std（只有 location）之外的
                #   **第二张表**提供（host.csv 只有 host）。`load_group_map` 是
                #   单表接口 → 这里按候选表逐个试，取第一个命中的。
                _metas = [prep.get("metadata_std") or prep.get("meta_csv") or ""]
                _droot = os.path.dirname(os.path.abspath(outdir))
                for _rel in ("data/host.csv", "data/sample_metadata.csv"):
                    _p = os.path.join(_droot, _rel)
                    if os.path.exists(_p) and _p not in _metas:
                        _metas.append(_p)
                for _m in _metas:
                    if _m and os.path.exists(_m):
                        _groups, _gnote = load_group_map(_m, _gc)
                        if _groups:
                            break
            _pdf = os.path.join(outdir, "neighbornet.pdf")
            _png = os.path.join(outdir, "neighbornet.png")
            _coords = os.path.join(outdir, "_st6_native", "nnet_radial.coords")
            _done = False
            if str(getattr(args, "splitstree_layout", "native") or "native") != "equal_angle":
                try:
                    from splitstree_bridge import native_coords
                    os.makedirs(os.path.dirname(_coords), exist_ok=True)
                    native_coords(os.path.join(outdir, "neighbornet.stree6"),
                                  _coords, width=2886, height=1496)
                    # ── 等长 (unweighted) 几何: 把每条边权重当 1 重放 DFS ──────
                    # 用户 2026-10-12: "取消枝长, 等长" —— 参考图 (SplitsTree GUI
                    # 截图) 是等长的规则网格, 而加权布局边长跨 344 倍 (放射长枝)。
                    # 转换器 `_unweight.py` 自检 = "原权重重放 == 原坐标" (残差 < 0.5)。
                    _unw = str(getattr(args, "splitstree_unweight", "on") or "on")
                    if _unw in ("on", "auto"):
                        try:
                            # 用 `splitstree_native_viz` 自己的 `unweight_coords()`
                            # (2026-10-13 起该能力已内联进脚本, 这里不再重复实现)。
                            from utils.splitstree_native_viz import (
                                unweight_coords as _uw_fn)
                            _coords_unw = _coords.replace(
                                ".coords", "_unweighted.coords")
                            _ok_uw, _new_uw = _uw_fn(_coords, _coords_unw)
                            if _ok_uw:
                                _coords = _new_uw
                                logger.info("  [splitstree] 等长布局 (unweighted) OK")
                            else:
                                logger.warning("  [splitstree] unweight 自检未通过, "
                                               "沿用加权布局")
                        except Exception as _ue:  # noqa: BLE001
                            logger.warning(f"  [splitstree] unweight 失败, 沿用加权: {_ue}")
                    from utils.splitstree_native_viz import render as _native_render
                    # show_legend=False: 组名已在图内 (分组弧线), 不再叠图例
                    # label_mode='ring': 145 个 tip 标签在中心密集区**必然**互压 ——
                    #   旧 'tip' 模式实测 **478 对重叠** (不可交付)。'ring' 把标签推到
                    #   网络外圈多层环上 + 细引线连回, 并做**环间避让**, 实测 **0 对**。
                    #   (第 11 轮定型: 每字符宽系数 0.8656 + 逐对精确间距 + 环间消解。)
                    _native_render(_coords, _pdf, _png, nexus=nexus,
                                   title=vname, groups=_groups or None,
                                   group_note=_gnote, label_mode="ring",
                                   label_rings=int(getattr(args, "splitstree_label_rings", 9) or 9),
                                   show_legend=False, lw_preset="auto")
                    # 面板版: 半栏 4.58 in, 无标题/图例/题注/标签 (Fig3 面板 D)
                    #   —— 符合「面板宽必须匹配源画布宽」, 否则字号打 5 折不可读。
                    #   保留分组弧线 + 组名 (这是第 9 轮新增的图内信息), 字号压到 6.5。
                    #   ⚠ lw_preset="auto": 等长坐标 → 自动切 flat (等粗)。用户口径
                    #     「useWeights = false」= 等长 + 等粗, 两者必须同时生效。
                    #   ⚠ canvas_aspect=1.0: 面板 D 与同排面板 C
                    #     (`phylogeny/tree_circular.pdf`, 4.580 x 4.580, aspect 1.0)
                    #     并排 → 画布比必须一致, 否则 D 比 C 矮 → 字母 D 错位
                    #     (实测等长版 D 的自然比 0.9105, 字母低 0.41 in)。
                    _panel = os.path.join(outdir, "neighbornet_panel.pdf")
                    _native_render(_coords, _panel, None, groups=_groups or None,
                                   label_mode="none", show_legend=False,
                                   fit_text=False, width_in=4.58,
                                   group_arcs=True, arc_label_fontsize=6.5,
                                   lw_preset="auto", canvas_aspect=1.0)
                    _done = True
                    results["outputs"]["splitstree_coords"] = _coords
                    logger.info("  [splitstree] 原生布局渲染 (SplitsTree6 SplitNetworkLayout)")
                except Exception as _e:  # noqa: BLE001
                    logger.warning(f"  [splitstree] 原生布局渲染失败, 回退自绘: {_e}")
            if not _done:
                from utils.splitstree_viz import render as _ss_render
                _ss_render(nexus, _pdf, _png,
                           groups=_groups or None, group_note=_gnote)
                # 面板版兜底: Fig3 面板 D 指向 neighbornet_panel.pdf,
                # 原生路径失败时退化为整页图副本 (不阻断 Fig3)
                _panel = os.path.join(outdir, "neighbornet_panel.pdf")
                if not os.path.exists(_panel) and os.path.exists(_pdf):
                    import shutil as _sh
                    _sh.copyfile(_pdf, _panel)
            if os.path.exists(_pdf):
                results["outputs"]["splitstree_pdf"] = _pdf
            if os.path.exists(_png):
                results["outputs"]["splitstree_png"] = _png
            if _groups:
                logger.info(f"  [splitstree] 分组着色 ({_gnote}): "
                            f"{len(set(_groups.values()))} 组")
        except Exception as _e:  # noqa: BLE001 — 渲染失败不拖垮网络产物
            logger.warning(f"  [splitstree] 出版级渲染失败 (nexus 保留): {_e}")
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
           "--perm", str(args.popgen_perm),
           # 2026-10-12: 结构化统计 JSON (全局+按群体+Fst 矩阵), 供多面板统计图渲染
           "--stats-json", os.path.join(popgen_dir, "popgen_stats.json")]
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
        # ── 统计可视化 (2026-10-12): 结构化 JSON → 多面板统计图 ──
        # π/θ/Hd/Tajima's D/Fu's Fs/Fst 此前只有纯文本报告没有图; pypopart 的
        # visualization 模块只画网络 (static_plot) 不画统计, 故自研渲染。
        # 渲染失败只降级 (报告文本与 checkpoint 不受影响)。
        try:
            stats_json = os.path.join(popgen_dir, "popgen_stats.json")
            if not os.path.exists(stats_json):
                logger.warning("  [popgen] popgen_stats.json 未产出, 跳过统计绘图")
            else:
                from utils.popgen_plots import render_popgen_stats
                pdf_p = os.path.join(popgen_dir, "popgen_stats.pdf")
                png_p = os.path.join(popgen_dir, "popgen_stats.png")
                with open(stats_json, encoding="utf-8") as _sj:
                    _stats = json.load(_sj)
                render_popgen_stats(_stats, pdf_p, png_p, title=vname)
                if os.path.exists(png_p):
                    results["outputs"]["popgen_stats_pdf"] = pdf_p
                    results["outputs"]["popgen_stats_png"] = png_p
                    n_pops = len(_stats.get("populations") or [])
                    _fp = (_stats.get("fst_pairwise") or {}).get("pops") or []
                    results["metrics"]["popgen_n_pops_plotted"] = n_pops
                    results["metrics"]["popgen_n_fst_pairs"] = len(_fp) * (len(_fp) - 1) // 2
                    logger.info(f"  [popgen] 统计图 → {os.path.basename(png_p)} "
                                f"({n_pops} 个群体)")
        except Exception as e:  # noqa: BLE001
            logger.warning(f"  [popgen] 统计可视化失败 (报告不受影响): {e}")
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
            # 2026-10-13 修 (P1): 硬编码 "data"/"ncbi_ref" 在 standard 布局下失配
            # (online 写的是 ph_dir('data')/ncbi_ref) → capheine 找不到 GBK 参考。
            gbks = sorted(Path(out_dir, ph_dir('data'), "ncbi_ref").glob("*.gbk")) + \
                   sorted(Path(out_dir, ph_dir('data'), "ncbi_ref").glob("*.gb"))  # online 产物
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
                    results["skip_reason"] = "viroid/no-CDS reference: positive selection (codon methods) not applicable"
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
                results["skip_reason"] = "viroid/no-CDS reference: positive selection (codon methods) not applicable"
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
    if gdr.returncode == 3:
        # 2026-10-13 修 (P1): 子进程约定退出码 3 = BEAST 工具不可用 (旧版 exit 0
        # 假成功 → checkpoint 被永久跳过)。
        logger.warning("  [gene_dating] BEAST 不可用, 部署未提交 (装好 BEAST 后重跑本 stage)")
        return
    if gdr.returncode == 0:
        # 成功判据加产物验证 (2026-10-13): 退出码 0 但目录空 = 部署脚本异常, 不写
        # checkpoint (否则下轮 _done('gene_dating') 恒真、永久跳过)。
        _gd_products = ([f for f in os.listdir(gene_dating_dir)
                         if f.endswith((".sh", ".xml", ".tsv"))]
                        if os.path.isdir(gene_dating_dir) else [])
        if not _gd_products:
            logger.warning("  [gene_dating] 部署退出码 0 但目录无产物, 不写 checkpoint")
            return
        results["outputs"]["gene_dating_dir"] = gene_dating_dir
        results["metrics"]["gene_dating_genes_dir"] = genes_dir
        logger.info(f"  [gene_dating] 部署完成 → {gene_dating_dir} (后台 BEAST 逐基因运行)")
    else:
        logger.warning(f"  [gene_dating] 失败: {(gdr.stderr or gdr.stdout)[-400:]}")


def run_stage_treetime(prep, aln_file, tree_file, args, logger, results) -> None:
    """treetime: TreeTime CLI 完整时间树 (timetree.nexus) — phymapr 路线前置。

    2026-10-12 新增 (STAGE_REFERENCE §19; utils/treetime_bridge.py)。
    与 clock 阶段的关系: clock 的 run_treetime_rtt 走 **Python API** 只做 RTT 回归口径
    (无注解 newick); 本 stage 走 **`python -m treetime` CLI**: clock-filter(IQD) +
    coalescent(skyline) + --confidence + 固定种子 → 带 `[&date=...]` 注解的
    timetree.nexus (phymapr `treeio::read.beast` 的硬性输入)。

    外部工具缺失 (phylo-treetime 未装/日期不足/CLI 失败) → **只告警不抛异常**:
    由调用方按 outputs 判定为"优雅 skip"(default 关闭的 stage, 缺失不阻断主线)。
    """
    out_dir = prep["output_dir"]
    out = os.path.join(out_dir, ph_dir('time'), "treetime")
    from utils.treetime_bridge import run_treetime_full
    logger.info(f"  [treetime] {prep['virus_name']}: TreeTime CLI 完整时间树 "
                f"(clock-filter={getattr(args, 'treetime_clock_filter', 4.0)}, "
                f"coalescent={getattr(args, 'treetime_coalescent', 'skyline')})...")
    tt = run_treetime_full(
        tree_file=tree_file, aln_file=aln_file,
        dates_csv=prep.get("dates_csv", ""), output_dir=out,
        clock_filter=getattr(args, "treetime_clock_filter", 4.0),
        coalescent=getattr(args, "treetime_coalescent", "skyline"),
        clock_rate=getattr(args, "treetime_clock_rate", None),
        log=LogCollector(logger))
    if not tt.get("success"):
        logger.warning(f"  [treetime] 未产出 (treetime 未装/日期不足/CLI 失败), "
                       f"skipped: {tt.get('error')}")
        return
    results["outputs"]["treetime_nexus"] = tt["timetree_nexus"]
    _m = tt.get("metrics") or {}
    for _k, _dst in (("clock_rate", "treetime_clock_rate"),
                     ("r_squared", "treetime_r_squared"),
                     ("n_tips_total", "treetime_n_tips_total"),
                     ("n_tips_used", "treetime_n_tips_used"),
                     ("n_tips_pruned_undated", "treetime_n_tips_pruned_undated"),
                     ("clock_filter_iqd", "treetime_clock_filter_iqd"),
                     ("coalescent", "treetime_coalescent"),
                     ("rng_seed", "treetime_rng_seed")):
        if _m.get(_k) is not None:
            results["metrics"][_dst] = _m[_k]
    # DRT 未通过 → 时间树仅供探索, 不进论文定年结论 (与 clock 阶段同一纪律:
    # "PSTVd 不报 TMRCA" 同源)
    if results["metrics"].get("drt_passed") is False:
        results["metrics"]["treetime_temporal_caveat"] = (
            "clock 阶段 DRT 未通过: 时间树仅供探索, 不做定年结论")
    logger.info(f"  [treetime] timetree.nexus → {tt['timetree_nexus']}")


def run_stage_phymapr(prep, results, args, tools, logger) -> None:
    """phymapr: 系统地理 TreeTime 路线 (BM-ASR 祖先重建 + 传播地图)。

    2026-10-12 新增 (STAGE_REFERENCE §20; utils/phymapr_bridge.py)。
    吃 treetime stage 的 timetree.nexus + 管线三通道元数据 (name/date/location),
    **不依赖 BEAST** —— 与本管线 BEAST 路线 (phylogeo→geo_paths) 方法不同源。
    外部依赖: Rscript + phymapr R 包 (缺 → 只告警, 由调用方判 skip, 不抛异常)。
    """
    out_dir = prep["output_dir"]
    timetree = (results.get("outputs", {}).get("treetime_nexus")
                or os.path.join(out_dir, ph_dir('time'), "treetime", "timetree.nexus"))
    meta = prep.get("metadata_std") or prep.get("meta_csv") or ""
    if not _usable_file(timetree):
        logger.warning(f"  [phymapr] 缺 timetree.nexus (先跑 treetime 阶段): {timetree} → skipped")
        return
    if not _usable_file(meta):
        logger.warning(f"  [phymapr] 缺元数据 (name/date/location): {meta} → skipped")
        return

    phy_dir = os.path.join(out_dir, ph_dir('geography'), "phymapr")
    rscript = (tools or {}).get("Rscript") or "Rscript"
    from utils.phymapr_bridge import build_phymapr_metadata, run_phymapr
    logger.info(f"  [phymapr] {prep['virus_name']}: phymapr BM-ASR 系统地理 "
                f"(map_style={getattr(args, 'phymapr_map_style', 'polygon')})...")
    bridge = build_phymapr_metadata(
        meta, os.path.join(phy_dir, "phymapr_metadata.tsv"), log=LogCollector(logger))
    if not bridge.get("success"):
        logger.warning(f"  [phymapr] 元数据桥失败 → skipped: {bridge.get('error')}")
        return
    pr = run_phymapr(
        timetree_nexus=timetree, meta_tsv=bridge["out_tsv"], output_dir=phy_dir,
        rscript_path=rscript,
        map_style=getattr(args, "phymapr_map_style", "polygon"),
        n_prune_early=getattr(args, "phymapr_n_prune_early", 4),
        direct_snp=getattr(args, "phymapr_direct_snp", 2.0),
        indirect_snp=getattr(args, "phymapr_indirect_snp", 5.0),
        distant_style=getattr(args, "phymapr_distant_style", "import_arrow"),
        pace=getattr(args, "phymapr_pace", 0.5),
        log=LogCollector(logger))
    if not pr.get("success"):
        logger.warning(f"  [phymapr] 未产出 (R 依赖缺失/计算失败) → skipped: "
                       f"{pr.get('error')}")
        return
    results["outputs"].update(pr.get("outputs") or {})
    if pr.get("summary_json"):
        results["outputs"]["phymapr_summary"] = pr["summary_json"]
    _m = pr.get("metrics") or {}
    for _k in ("n_pathways", "n_episodes", "n_locations",
               "n_tips_used", "n_tips_total", "n_tips_dropped"):
        if _m.get(_k) is not None:
            results["metrics"][f"phymapr_{_k}"] = _m[_k]
    for _k, _v in (_m.get("transmission_counts") or {}).items():
        results["metrics"][f"phymapr_conf_{_k}"] = _v
    logger.info(f"  [phymapr] 完成 → {phy_dir} "
                f"(pathways={_m.get('n_pathways')}, episodes={_m.get('n_episodes')})")


def _find_phylogeo_mcc(out_dir: str) -> Optional[str]:
    """定位 phylogeo 产出的 MCC 树。

    2026-09-16 抽成函数: 原先 phylogeo 阶段内联写了一遍 glob, geo 阶段要用同一份
    候选顺序 (优先 `merged/mcc.tree`, 再退化为 `*mcc*.tree`); 两处各写一份必漂移。
    """
    phylogeo_dir = os.path.join(out_dir, ph_dir('geography'), "phylogeography")
    cands = list(Path(phylogeo_dir).glob("merged/mcc.tree")) \
        + sorted(Path(phylogeo_dir).glob("*mcc*.tree"))
    return str(cands[0]) if cands else None


def _phylogeo_log_candidates(phylogeo_dir: str) -> List[str]:
    """phylogeo 的 BEAST log 候选, **按可信度排序** (最佳在前)。

    2026-10-12 修 (P1): 旧版直接
        list(Path(d).glob("merged/merged.log"))
        + list(Path(d).glob("phylogeo_sky*.log"))
        + list(Path(d).glob("phylogeo_beast1*.log"))
    然后取 `[0]`。问题: `phylogeo_sky*.log` 会同时命中 BEAST **控制台输出**
    `phylogeo_sky1_run.log` —— 它只有 1 列、**0 个** `location.indicators` 列。
    文件系统枚举顺序下它常排在最前, 于是**没有** merged/merged.log 时 (实测 GCVA
    快照正是这种形态):
      · generate_spread3_json → "log 中未解析到任何 BSSVS 迁移路线" → spread3.json 不生成
      · extract_migration_bf  → migration_bf.csv 不生成 → 多面板图 Panel B 缺失
    实测同目录三份 log 的 indicator 列数: phylogeo_sky1.log=30,
    phylogeo_sky1.location.rates.log=30, phylogeo_sky1_run.log=**0**。

    现口径: ① 优先 merged/merged.log (burnin 已处理, 样本最全);
            ② 主链 log `<prefix><N>.log` (正则精确匹配, 排除任何二级后缀);
            ③ 其余候选 (仍排除 *_run.log);
    且**全部排序** —— 文件系统枚举顺序不可靠, 会导致同一目录两次跑出不同结果。
    """
    out: List[str] = []
    merged = os.path.join(phylogeo_dir, "merged", "merged.log")
    if os.path.exists(merged) and os.path.getsize(merged) > 0:
        out.append(merged)
    cands = []
    for pat in ("phylogeo_sky", "phylogeo_beast1"):
        for p in Path(phylogeo_dir).glob(f"{pat}*.log"):
            if p.name.endswith("_run.log"):
                continue            # BEAST 控制台输出, 不是主 log
            cands.append(p)
    main = sorted(p for p in cands
                  if re.fullmatch(r"phylogeo_(?:sky|beast1)\d+\.log", p.name))
    main_set = set(main)
    rest = sorted(p for p in cands if p not in main_set)
    out += [str(p) for p in main] + [str(p) for p in rest]
    return out


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
    merge_fn = merge_fn or (lambda wd, s: _default_merge(wd, s))
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
                # 2026-10-13: 原先不传 log_every → 吃 bridge 默认 5000。
                # 10M 链只落 2000 样本/链 (GCVA 实测 merged n=5403), ucld.mean
                # 合并 ESS 仅 48。改为按目标样本数反算 (10M → 1000)。
                log_every=beast_log_every(args.phylogeo_chain),
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
                # 2026-10-12: merge() 的 background 形参已删 (恒同步), 不再传。
                merge_chains(phylogeo_dir, st)
        elif st:
            # 前台等待闭环 (2026-09-02): 轮询到链完成再合并, 不再 pending 等二次重跑
            ok = _wait_chains_and_merge(
                phylogeo_dir, logger, label="phylogeo",
                merge_fn=lambda wd, s: merge_chains(wd, s))
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
    # 2026-10-12 修 (P1): 候选改由 _phylogeo_log_candidates 给出 —— 旧版裸 glob
    # 会把控制台输出 phylogeo_sky1_run.log 排在首位 (0 个 indicator 列) → 迁移矩阵
    # 解析成 0 条 → spread3.json / migration_bf.csv 双双不生成。
    log_files = _phylogeo_log_candidates(phylogeo_dir)
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
            + sorted(Path(phylogeo_dir).glob("phylogeo_sky*.trees"))
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

def beast_log_every(chain_length: int, target_samples: int = 10_000) -> int:
    """BEAST logEvery: 按「每链目标样本数」反算采样间隔。

    ⚠ 历史缺陷 (2026-10-13 修): 原调用处写的是
        log_every = max(1000, args.beast_chain // 100)
    即 **无论链多长都只取 ~100 个样本**。GCVA 实测 chain_length=5,000,000
    → logEvery=50,000 → 每链仅 101 个样本 (烧 10% burnin 后 91), 三链合并 273。
    后果:
      · **无法诊断收敛**: 合并序列 ACF 在 lag 30 仍 >0.7, ACT 被估到 7e6 世代
        (比链长还大), ess_results.csv 全线 4–10 —— 粗采样把「采样间隔」当成了
        「自相关时间」;
      · **MCC 树只有 273 棵**, treeannotator 的 height_95%_HPD 宽度退化到 ~1e-14;
      · **skyline 只有 273 点**, BSP 阶梯图剧烈锯齿。
    现按 target_samples (默认 10,000) 反算, 使链长与样本数解耦:
      5M → 500;  20M → 2,000;  50M → 5,000;  200M → 20,000。
    下限 100 防 chain_length 过小产生 0 间隔。

    注: 严格说 ESS ≈ chain_length / ACT_generations, **与 logEvery 无关** ——
    本修复不直接抬 ESS, 但它是「能否测量 ESS」「MCC/HPD 是否可用」的前提。
    """
    cl = int(chain_length or 0)
    ts = max(1, int(target_samples))
    return max(100, cl // ts)


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
            chain_length=args.beast_chain, log_every=beast_log_every(args.beast_chain),
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
                        # 2026-10-12: 原先传 background=True —— 而 merge() 从未使用该参数
                        # (恒同步), 名字却暗示"已转后台可退出", 是误导。参数已删。
                        merge_chains(beast_dir, st)
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
                            merge_fn=lambda wd, s: merge_chains(wd, s))
                        if ok and os.path.exists(merged_log):
                            results["outputs"]["beast_log"] = merged_log
                            results["outputs"]["beast_trees"] = os.path.join(beast_dir, "merged", "mcc.tree")
                        elif not ok:
                            results["beast_pending"] = True   # 等待失败退回 pending 语义
                            if os.path.exists(chain1_log):
                                results["outputs"]["beast_log"] = chain1_log
                                results["outputs"]["beast_trees"] = os.path.join(beast_dir, "beast11.trees")
                # 收敛门禁: 把 merged/convergence_gate.json 的结论带入 metrics,
                # 供报告层据实标注 (未收敛 → 数值仅供探索)。⚠ 门禁只判"MCMC 混合",
                # 与"数据有无时间信号"独立 (后者由 DRT/R² 管)。
                _gate = os.path.join(beast_dir, "merged", "convergence_gate.json")
                if os.path.exists(_gate):
                    try:
                        with open(_gate, encoding="utf-8") as _gf:
                            _g = json.load(_gf)
                        results["metrics"]["beast_convergence_verdict"] = _g.get("verdict")
                        results["metrics"]["beast_worst_ess"] = _g.get("worst_ess")
                        if _g.get("verdict") != "converged":
                            results["metrics"]["beast_convergence_caveat"] = (
                                "MCMC 未收敛 (关键参数最小 ESS "
                                f"{_g.get('worst_ess')} < 200, 最差: {_g.get('worst_param')}): "
                                "参数估计仅供探索, 不做定量结论")
                            log.warning(f"  [beast] ⚠ 收敛门禁: {results['metrics']['beast_convergence_caveat']}")
                    except Exception as _ge:
                        log.warning(f"  [beast] 读取 convergence_gate.json 失败: {_ge}")
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


def count_beast_log_rows(path: str) -> int:
    """BEAST .log 的**数据行数** (跳过 '#' 注释与 'state' 表头行)。

    2026-10-03 引入: 退化日志闸门共用 (链未跑完 / logcombiner burnin 超界 →
    近空 merged.log; 该日志会让 ESS 走"样本数"兜底、BSP 画空白图)。
    2026-10-13 恢复: 随记账机制一并丢失。
    返回 0 = 空/仅表头/不可读。
    """
    try:
        with open(path, encoding="utf-8", errors="replace") as f:
            return sum(1 for L in f
                       if L.strip() and not L.startswith("#")
                       and not L.startswith("state"))
    except OSError:
        return 0


def pipeline_exit_code(all_results) -> tuple:
    """退出码约定 (2026-10-03, P1-1): failed→1; blocked/pending→2; 干净→0。

    旧版 main 无论如何退出 0 → 自动化/`&&` 链/per-virus driver 一律"成功"。
    2026-10-13 恢复: 该函数在 10-12 工作区重写中丢失 (RUN_GUIDE:298-301 仍按本
    约定验收)。返回 (code, 说明)。工作目录模式与批量模式共用 (两处出口都要调)。
    """
    n_fail = sum(len(r.get("stages_failed", [])) for r in all_results)
    n_blocked = sum(len(r.get("stages_blocked", [])) for r in all_results)
    n_pending = sum(1 for r in all_results
                    if r.get("incomplete") and not r.get("stages_blocked"))
    if n_fail:
        return 1, f"{n_fail} 个 stage 失败 — 见 Stages_failed"
    if n_blocked or n_pending:
        extra = f", 另有 {n_pending} 个病毒有 pending 链" if n_pending else ""
        return 2, f"结果不完整 (受阻 stage {n_blocked} 个{extra}) — 修好前置后重跑"
    return 0, "全部 stage 干净"


def process_virus(prep, stages: List[str], args, tools, logger) -> Dict:
    """单个病毒按 stage 清单顺序执行, 返回汇总 results (含 stage 状态追踪)"""
    vname = prep["virus_name"]
    out_dir = prep["output_dir"]
    log = logger
    # blocked = 真受阻 (前置缺失/本 stage 未产出/降级继续); skip_reasons 记录每个跳过原因
    # (2026-10-13 恢复: 五态记账在 10-12 重写中丢失 → stages_blocked/Incomplete/退出码全失效)
    st = {"done": [], "skipped": [], "failed": [], "blocked": [], "skip_reasons": {}}
    results = {"virus": vname, "success": False, "outputs": {}, "metrics": {},
               "stages_done": [], "stages_skipped": [], "stages_failed": [],
               "virus_dir": out_dir}  # virus_dir: report_builder 按模块扫产物用 (2026-09-02)

    log.info(f"\n{'='*60}")
    log.info(f"Processing {vname} | stages: {','.join(stages)}")
    log.info(f"{'='*60}")

    # 2026-10-11: 布局与既有结果树不一致的**静默**陷阱。默认布局是 standard, 而存量
    # 结果树可能是 legacy —— 按默认跑不会报错, 只是把产物写到另一套目录 (01_data/
    # 08_report), 老目录一个字节都不变, 用户以为"跑了但什么都没变"。
    # 实测同一天就踩到: 不带 --io-layout 跑 PSTVd report → 新建 01_data/ 与 08_report/。
    _lay_hit = detect_layout_from_tree(out_dir)
    _lay_active = os.environ.get("MMPV_IO_LAYOUT") or "standard"
    if _lay_hit and _lay_hit != _lay_active:
        log.error(f"  [layout] ⚠ 工作目录的既有产物是 **{_lay_hit}** 布局, 但本次生效的是 "
                  f"**{_lay_active}** —— 产物会写到另一套目录里, 既有结果不会被更新。"
                  f"请改用 `--io-layout {_lay_hit}` 或 `export MMPV_IO_LAYOUT={_lay_hit}` 重跑。"
                  f" (目录: {out_dir})")

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

    def _mark(stage, meta=None, state=None):
        fp = _inputs_of(stage)
        if fp is not None:
            meta = dict(meta or {})
            meta.setdefault('inputs', fp)
        mark_stage_done(out_dir, stage, meta)
        # 2026-10-11: 同步落盘该 stage 的 outputs/metrics 快照, 供 skip 分支回填。
        # state 形如 {"outputs": {...}, "metrics": {...}, "extra": {...}}
        if state:
            save_stage_state(out_dir, stage,
                             outputs=state.get("outputs"),
                             metrics=state.get("metrics"),
                             extra=state.get("extra"))
        st["done"].append(stage)

    def _restore(stage):
        """把 <stage>.state.json 里的 outputs/metrics 并回 results, 返回 extra。

        2026-10-11 新增 (修 P0 "重跑即空汇总")。所有 checkpoint-skip 分支都调它,
        使"跳过"与"真跑"对下游**完全等价** —— 汇总表 Beta/R²/P/TMRCA 不再为空。
        返回 extra dict, 供调用方恢复 stage 间局部变量 (aln_file/tree_file/...)。
        """
        data = load_stage_state(out_dir, stage)
        if not data:
            # 2026-10-11: 快照缺失 (旧代码跑出的老数据 / 结果树被 cp -a / 手工删过)
            # → 从磁盘产物重建, 并**顺手落盘快照**实现自愈 (下次直接读快照)。
            # 不这么做的话, 存量结果永远补不出指标 —— 而重建用的都是当次真跑写下的
            # 同一份权威产物 (见 rediscover_stage_state 的 docstring)。
            data = rediscover_stage_state(out_dir, stage)
            if data:
                log.info(f"  [{stage}] 无状态快照 → 已从磁盘产物重建 "
                         f"({len(data.get('metrics') or {})} 项指标, "
                         f"{len(data.get('outputs') or {})} 项产物), 并落盘快照")
                save_stage_state(out_dir, stage, outputs=data.get("outputs"),
                                 metrics=data.get("metrics"), extra=data.get("extra"))
        if not data:
            log.warning(f"  [{stage}] ⚠ checkpoint 已跳过但无状态快照 "
                        f"(.checkpoints/{stage}{STAGE_STATE_SUFFIX} 缺失), "
                        f"且磁盘产物不足以重建 —— "
                        f"本次运行的汇总表将缺少该 stage 的指标。"
                        f"如需补全: 删除 .checkpoints/{stage}.done 后重跑该 stage。")
            return {}
        for k, v in (data.get("outputs") or {}).items():
            results["outputs"].setdefault(k, v)
        for k, v in (data.get("metrics") or {}).items():
            # 与 summary_metrics._pick 一致: None/"" 视为缺
            if results["metrics"].get(k) in (None, ""):
                results["metrics"][k] = v
        return data.get("extra") or {}

    def _skip(stage, reason="checkpoint", detail=""):
        """stage 未执行/未产出。reason 取值 (2026-10-03 引入, 区分「正常跳过」与「真受阻」;
        2026-10-13 恢复: 该五态语义在 10-12 工作区重写中丢失 → blocked 不再记账):

        · checkpoint — 已完成, 本轮跳过 (正常)
        · na         — 数据/配置不适用 (如 viroid 无 CDS、无宿主通道; 正常)
        · pending    — 后台链条运行中, 下次跑继续 (正常但结果未定)
        · degraded   — 本 stage 失败后降级继续 (有问题, 进 blocked)
        · blocked    — 前置产物缺失 / 本 stage 未产出 (有问题, 进 blocked)

        历史事故: 旧版 `_skip` 不计入 failed → `success` 仍为 True,
        geo/geo_paths/gene_dating 等整块缺席也报成功 (要跑到第二/三遍才出结果)。
        """
        st["skipped"].append(stage)
        st["skip_reasons"][stage] = reason
        if detail:
            log.info(f"  [{stage}] skip({reason}): {detail}")
        elif reason and reason != "checkpoint":
            log.info(f"  [{stage}] 跳过 ({reason})")
        if reason in ("degraded", "blocked"):
            st["blocked"].append(f"{stage}({reason})")

    def _fail(stage, msg=""):
        st["failed"].append(stage)
        if msg:
            log.warning(f"  [{stage}] {msg}")

    def _snap():
        """记录当前 outputs/metrics 快照, 供 _delta 求增量。"""
        return (dict(results.get("outputs") or {}),
                dict(results.get("metrics") or {}))

    def _delta(before):
        """相对 _snap() 的增量 {outputs, metrics}。

        2026-10-11: host / geo / geo_paths / phylogeo / gene_dating 这些 stage 是
        **直接往 results 里写** (签名吃 results 参数), 不返回 dict —— 无法像
        align/tree 那样拿到"本次产出"。用增量快照把它们也纳入状态持久化,
        否则这些 stage 一被 checkpoint 跳过, 其指标就从汇总里消失。
        """
        bo, bm = before
        no = {k: v for k, v in (results.get("outputs") or {}).items()
              if k not in bo}
        nm = {k: v for k, v in (results.get("metrics") or {}).items()
              if k not in bm or bm.get(k) in (None, "")}
        return {"outputs": no, "metrics": nm}

    def _mark_skipped(stage, reason="", state=None):
        """记录「不适用而跳过」—— 写 <stage>.skipped 让跳过也具备幂等性。

        2026-10-11 新增。此前 capheine 的"类病毒无 CDS"分支只调 _skip(),
        **不写任何标记** → 每次重跑都要重新跑一遍 gbk_extractor 才能再次得出
        "不适用"的结论 (慢且日志噪音大); 而 rdp5 的短序列分支更糟 —— 它直接
        return 原比对被当成成功, 写下了 rdp5.done (见 run_stage_rdp5 的修复)。
        """
        try:
            cp_dir = os.path.join(out_dir, CHECKPOINT_DIR)
            os.makedirs(cp_dir, exist_ok=True)
            with open(os.path.join(cp_dir, f"{stage}.skipped"), "w",
                      encoding="utf-8") as f:
                json.dump({"stage": stage, "status": "skipped",
                           "reason": reason,
                           "time": datetime.now().isoformat(timespec="seconds")},
                          f, ensure_ascii=False, indent=1)
        except OSError as e:
            log.warning(f"  [{stage}] 跳过标记写入失败 (非阻断): {e}")
        if state:
            save_stage_state(out_dir, stage,
                             outputs=state.get("outputs"),
                             metrics=state.get("metrics"),
                             extra=state.get("extra"))
        _skip(stage, 'na', reason)

    def _skipped_marker(stage) -> str:
        """读 <stage>.skipped 的 reason; 不存在返回 ""。"""
        p = os.path.join(out_dir, CHECKPOINT_DIR, f"{stage}.skipped")
        if not os.path.exists(p):
            return ""
        try:
            with open(p, encoding="utf-8") as f:
                return json.load(f).get("reason", "") or ""
        except Exception:
            return ""

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
            # 2026-10-11: 回填 outputs/metrics 快照 (修 P0 "重跑即空汇总")
            _ex = _restore('metadata')
            if not prep.get("metadata_std") and _ex.get("metadata_std"):
                prep["metadata_std"] = _ex["metadata_std"]
                args.metadata = _ex["metadata_std"]
        else:
            mr = run_stage_metadata(prep, args, log)
            results["outputs"].update(mr.get("outputs", {}))
            results["metrics"]["metadata"] = mr.get("metrics", {})
            if mr.get("success") and mr.get("std_csv"):
                # 关键: 把下游 metadata 输入换成标准化产物
                prep["metadata_std"] = mr["std_csv"]
                prep["_metadata_original"] = getattr(args, "metadata", None)
                args.metadata = mr["std_csv"]
                _mark('metadata', state={"outputs": mr.get("outputs"),
                                         "metrics": {"metadata": mr.get("metrics", {})},
                                         "extra": {"metadata_std": mr["std_csv"]}})
            else:
                # 治理失败不阻断流程, 显式告警 (降级为直接用原 metadata)
                _skip('metadata', 'degraded',
                      f"治理失败, 降级使用原 metadata: {mr.get('error')}")

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
            _restore('online')            # 2026-10-11: 回填 outputs/metrics 快照
        else:
            r = run_stage_online(prep, args, log)
            if r.get("ncbi_fasta"):
                results["outputs"]["ncbi_fasta"] = r["ncbi_fasta"]
                _mark('online', state={"outputs": {"ncbi_fasta": r["ncbi_fasta"]}})
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
            if os.path.exists(cf) and os.path.getsize(cf) > 0:
                prep["_raw_combined_fasta"] = prep.get("combined_fasta")
                prep["combined_fasta"] = cf
                results["outputs"]["clean_fasta"] = cf
            else:
                # 2026-10-11: 产物缺失/空 → 不能静默沿用上游原序列, 必须喊出来。
                log.error(f"  [clean] ⚠ checkpoint 标记完成但产物缺失或为空: {cf} —— "
                          f"下游将退回未清洗序列 (可能含空格/小写/非法字符)。"
                          f"请删除 .checkpoints/clean.done 后重跑 clean。")
            _ex = _restore('clean')       # 2026-10-11: 回填 outputs/metrics 快照
            if _ex.get("clean_fasta") and os.path.exists(_ex["clean_fasta"]):
                prep["combined_fasta"] = _ex["clean_fasta"]
        else:
            cr = run_stage_clean(prep, args, log)
            results["outputs"].update(cr.get("outputs", {}))
            results["metrics"]["clean"] = cr.get("metrics", {})
            if cr.get("success") and cr.get("clean_fasta"):
                # 关键: 把比对输入换成清洗后的序列
                prep["_raw_combined_fasta"] = prep.get("combined_fasta")
                prep["combined_fasta"] = cr["clean_fasta"]
                _mark('clean', state={"outputs": cr.get("outputs"),
                                      "metrics": {"clean": cr.get("metrics", {})},
                                      "extra": {"clean_fasta": cr["clean_fasta"]}})
            else:
                # 清洗失败不阻断流程, 但要显式告警 (降级为直接用原序列)
                _skip('clean', 'degraded',
                      f"清洗失败, 降级使用原始序列: {cr.get('error')}")

    # ── align: MAFFT 比对 + 内联 saturation (纯Python秒级, 合并运行) ──
    if 'align' in stages:
        if _done('align'):
            log.info("  [align] checkpoint OK, skip"); _skip('align')
            aln_file = os.path.join(out_dir, ph_dir('phylogeny'), "mafft.aln.fasta")
            # 2026-10-11 修 (P0): 原判定只有 `os.path.exists` —— 0 字节文件被当成
            # 合法比对沿用。PSTVd 实测: 工作副本 pstvd_mafft.aln.fasta 被截成
            # 0 字节, 该分支照样 `results["outputs"]["alignment"]=…` 并让
            # rdp5/tree/clock 全部吃空文件 → 下游一串莫名其妙的失败。
            # 现要求"存在且非空", 且回退候选也必须非空。
            if not _usable_alignment(aln_file):
                _ex0 = _restore('align')
                cand = _ex0.get("aln_file")
                aln_file = cand if _usable_alignment(cand) else None
                if aln_file:
                    log.warning(f"  [align] 主比对 {Path(aln_file).name} 缺失/为空, "
                                f"按状态快照回退到: {aln_file}")
            if aln_file is None:
                _pf = prep.get("combined_fasta")
                if _usable_alignment(_pf):
                    aln_file = _pf
                else:
                    log.error("  [align] ⚠ checkpoint 标记完成但比对缺失或为空 "
                              f"(主路径与 combined_fasta={_pf} 均不可用) —— "
                              "比对是全部下游 stage 的唯一输入, 拒绝以空输入继续。"
                              "请恢复比对文件后删除 .checkpoints/align.done 重跑。")
            results["outputs"]["alignment"] = aln_file
            _restore('align')             # 2026-10-11: 回填 outputs/metrics 快照
            # saturation 随 align 产出; 老数据补跑
            if not _done('saturation') and aln_file and os.path.exists(aln_file):
                sr = run_stage_saturation(prep, aln_file, args, log)
                results["outputs"].update(sr.get("outputs", {}))
                results["metrics"].update(sr.get("metrics", {}))
                if sr.get("success"):
                    _mark('saturation', state={"outputs": sr.get("outputs"),
                                               "metrics": sr.get("metrics")})
        else:
            r = run_stage_align(prep, tools, args, log, force=force)
            if r.get("success") and r.get("alignment"):
                aln_file = r["alignment"]
                results["outputs"]["alignment"] = aln_file
                _mark('align', state={"outputs": {"alignment": aln_file},
                                      "extra": {"aln_file": aln_file}})
                # 内联 saturation: 比对一旦产出立即验饱和 (闸门逻辑)
                sr = run_stage_saturation(prep, aln_file, args, log)
                results["outputs"].update(sr.get("outputs", {}))
                results["metrics"].update(sr.get("metrics", {}))
                if sr.get("success"):
                    _mark('saturation', state={"outputs": sr.get("outputs"),
                                               "metrics": sr.get("metrics")})
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
            _restore('splitstree')        # 2026-10-11: 回填 outputs/metrics 快照
        else:
            ssr = run_stage_splitstree(prep, aln_file, args, log)
            results["outputs"].update(ssr.get("outputs", {}))
            results["metrics"].update(ssr.get("metrics", {}))
            if ssr.get("success"):
                _mark('splitstree', state={"outputs": ssr.get("outputs"),
                                           "metrics": ssr.get("metrics")})
            else:
                _fail('splitstree', '分裂网络未产出')

    # ── rdp5: 重组检测 → 删重组子+重比对 (默认) / mask (可选) ──
    # 2026-10-11 修 (P0): run_stage_rdp5 改为返回**结构化状态**, 以区分三态:
    #   ok(产出干净比对) / skipped(方法不适用, 如类病毒 379bp) / failed(真失败)。
    # 旧版把"短序列不适用"也 return 原比对 → 调用方 `if masked:` 为真 →
    # 写下了 rdp5.done 并把**原始比对**登记成 masked_alignment —— 于是
    #   · 汇总表 Stages_done 里出现 rdp5 (实际没做)
    #   · 下游以为重组已被处理, 拿到的是未处理比对
    #   · 下次重跑走 skip 分支去找 rdp5_dropped.realn.fasta (不存在) → 前后矛盾
    # PSTVd 实测正是如此 (rdp5.done 落于 2026-09-02 20:52, 而 recomb/rdp5/ 里
    # 只有一次 09-01 的失败残留, 没有任何 CSV)。
    if 'rdp5' in stages and aln_file:
        _rdp_skip_reason = _skipped_marker('rdp5')
        if _rdp_skip_reason:
            log.info(f"  [rdp5] 已记录不适用, skip: {_rdp_skip_reason}")
            _skip('rdp5')
            _restore('rdp5')
        elif _done('rdp5'):
            log.info("  [rdp5] checkpoint OK, skip"); _skip('rdp5')
            # 恢复比对: 优先新策略产物, 回退旧 masked (老数据兼容)
            r_aln = os.path.join(out_dir, ph_dir('recomb'), "rdp5", "rdp5_dropped.realn.fasta")
            if not _usable_alignment(r_aln):
                r_aln = os.path.join(out_dir, ph_dir('recomb'), "rdp5", "masked", "masked_N.fasta")
            if _usable_alignment(r_aln):
                aln_file = r_aln
                results["outputs"]["masked_alignment"] = r_aln
            else:
                # checkpoint 说做过了, 但产物不在 —— 不能静默沿用未处理比对
                log.error("  [rdp5] ⚠ checkpoint 标记完成但既无 rdp5_dropped.realn.fasta "
                          "也无 masked/masked_N.fasta —— 本次将沿用**未去重组**的比对。"
                          "如需重建请删除 .checkpoints/rdp5.done 后重跑 rdp5。")
            _restore('rdp5')
        else:
            # 2026-10-12 (test E): 环境无 RDP5 (服务器未装 / wine 缺失) → **预检脚本**,
            # 缺失即优雅跳过: _skip('rdp5', 'na') —— 不 fail (非数据问题)、不写
            # checkpoint (装好 RDP5 后重跑该 stage 即生效)。旧实现要跑到
            # recombination_analysis 里才失败, 报出的是 "RDP5 failed (exit=...)"
            # 这种误导信息。
            rdp5_script = os.path.expanduser(getattr(args, "rdp5_script", "") or "")
            if not os.path.exists(rdp5_script):
                log.warning(f"  [rdp5] RDP5 脚本不存在 ({rdp5_script}) —— 本环境无 RDP5, "
                            f"跳过该 stage (不 fail)")
                _skip('rdp5', 'na')
            else:
                rr = run_stage_rdp5(aln_file, prep, tools, args, log)
                if rr is None:
                    # 直调守卫同款: 脚本缺失 (预检与调用之间被删) → 同判 skip
                    _skip('rdp5', 'na')
                    rr = {}
                _st = rr.get("status")
                if _st == "ok":
                    masked = rr.get("alignment")
                    aln_file = masked
                    results["outputs"]["masked_alignment"] = masked
                    _mark('rdp5', state={"outputs": {"masked_alignment": masked},
                                         "metrics": rr.get("metrics"),
                                         "extra": {"aln_file": masked}})
                    # 后处理: 重组事件验证层 (validate_recombinants, 生成验证任务 + 汇总表)
                    rdp5_csv = os.path.join(out_dir, ph_dir('recomb'), "rdp5", f"{prep['virus_name']}.csv")
                    if os.path.exists(rdp5_csv):
                        run_postprocess("validate_recombinants", log,
                                        "--events", rdp5_csv,
                                        "--alignment", aln_file,
                                        "--outdir", os.path.join(out_dir, ph_dir('recomb'), "rdp5", "validation"))
                elif _st == "skipped":
                    # 方法不适用 (短序列/类病毒): 记 skip 且**写 .skipped 标记**保证幂等,
                    # 绝不写 .done、绝不把原比对冒充成去重组产物。
                    log.info(f"  [rdp5] 不适用, 记 skip: {rr.get('skip_reason', '')}")
                    _mark_skipped('rdp5', rr.get('skip_reason', ''),
                                  state={"metrics": rr.get("metrics")})
                else:
                    _fail('rdp5', rr.get("error") or 'RDP5/mask 未产出干净比对')

    # ── tree: IQ-TREE 建树 ──
    if 'tree' in stages and aln_file:
        if _done('tree'):
            log.info("  [tree] checkpoint OK, skip"); _skip('tree')
            tree_file = os.path.join(out_dir, ph_dir('phylogeny'), "iqtree.treefile")
            if not _usable_file(tree_file):
                tree_file = prep.get("tree_file")
            results["outputs"]["tree"] = tree_file
            if not _usable_file(tree_file):
                log.error("  [tree] ⚠ checkpoint 标记完成但主树缺失或为空 "
                          f"({tree_file}) —— 分子钟/定年/地理全部依赖它。")
            _restore('tree')              # 2026-10-11: 回填 outputs/metrics 快照
        else:
            r = run_stage_tree(prep, aln_file, tools, args, log, force=force)
            if r.get("success") and r.get("tree"):
                tree_file = r["tree"]
                results["outputs"]["tree"] = tree_file
                _mark('tree', state={"outputs": {"tree": tree_file},
                                     "extra": {"tree_file": tree_file}})
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
            _restore('popgen')            # 2026-10-11: 回填 outputs/metrics 快照
        else:
            pr = run_stage_popgen(prep, aln_file, args, log)
            results["outputs"].update(pr.get("outputs", {}))
            results["metrics"].update(pr.get("metrics", {}))
            if pr.get("success"):
                _mark('popgen', state={"outputs": pr.get("outputs"),
                                       "metrics": pr.get("metrics")})
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
            _restore('host')              # 2026-10-11: 回填 outputs/metrics 快照
        else:
            _b = _snap()
            run_stage_host(prep, results, aln_file, tree_file, args, log)
            if results.get("outputs", {}).get("host_tree") or results.get("outputs", {}).get("host_diff_plot"):
                _mark('host', state=_delta(_b))
            else:
                _hc = prep.get("host_csv") or ""
                if not (_hc and os.path.exists(_hc)):
                    _skip('host', 'na', '无宿主通道 (host.csv), 数据不适用')
                else:
                    _skip('host', 'blocked', '有宿主通道但未产出 (见上方 host 日志)')

    # ── capheine: 正选择分析 (Phase S) ──
    if 'capheine' in stages:
        _cap_skip = _skipped_marker('capheine')
        if _cap_skip:
            # 2026-10-11: "不适用"也具备幂等性 —— 类病毒无 CDS 时不必每次重跑
            # gbk_extractor 才能再次得出同一结论。
            log.info(f"  [capheine] 已记录不适用, skip: {_cap_skip}")
            _skip('capheine')
            _restore('capheine')
        elif _done('capheine'):
            log.info("  [capheine] checkpoint OK, skip"); _skip('capheine')
            _restore('capheine')          # 2026-10-11: 回填 outputs/metrics 快照
        else:
            cr = run_stage_capheine(prep, args, log)
            results["outputs"].update(cr.get("outputs", {}))
            results["metrics"].update(cr.get("metrics", {}))
            if cr.get("success"):
                _mark('capheine', state={"outputs": cr.get("outputs"),
                                         "metrics": cr.get("metrics")})
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
                    # 2026-10-11: 同时写 <stage>.skipped 标记 → 跳过具备幂等性,
                    # 且原因落盘 (此前只在日志里出现过一次, 报告里看不到"为什么没有")。
                    log.info(f"  [capheine] skip: {cr['skip_reason']}")
                    _mark_skipped('capheine', cr['skip_reason'],
                                  state={"metrics": cr.get("metrics")})
                else:
                    _fail('capheine', 'capheine 未产出 (缺 GBK/参考 CDS 或计算失败)')

    # ── clock: 分子钟联合 (整体+分基因 RTT/LTT + DRT) ──
    if 'clock' in stages and aln_file and tree_file:
        if _done('clock'):
            log.info("  [clock] checkpoint OK, skip"); _skip('clock')
            # 2026-10-11 修 (P0): 原分支**只打日志**。clock 是 headline_metrics
            # 的首选来源 (td_mean_rate / td_relaxed_clock_p / td_tmrca 来自
            # treedater; pic_rate / pic_r 来自 PIC; beta/r_squared/p_value 来自
            # TreeTime RTT), 一旦被 checkpoint 跳过, 汇总表 Beta/R²/P/TMRCA
            # 立刻全空 —— 这正是"重跑一次 report 就得到空汇总"的直接原因。
            _ex = _restore('clock')
            if _ex.get("clock_summary") and os.path.exists(_ex["clock_summary"]):
                results["outputs"].setdefault("clock_summary", _ex["clock_summary"])
        else:
            prep["_cp"] = lambda name: _done(name)
            cr = run_stage_clock(prep, aln_file, tree_file, tools, args, log,
                                 results=results)
            results["outputs"].update(cr.get("outputs", {}))
            results["metrics"].update(cr.get("metrics", {}))
            if cr.get("outputs") or results.get("gene_results"):
                _mark('clock', state={"outputs": cr.get("outputs"),
                                      "metrics": cr.get("metrics"),
                                      "extra": {"clock_summary":
                                                (cr.get("outputs") or {}).get("clock_summary")}})
            else:
                _fail('clock', 'RTT/LTT/DRT 无输出')

    # rtt/temporal/genes 已并入 clock (此处不再单独跑; 保留函数供 --stage 兼容)

    # ── beast: 全长 BEAST 定年 ──
    def _beast_log_candidates():
        """后验日志候选 (合并优先; 与 main 末尾 publication_figures 的探测顺序一致)。"""
        bd = os.path.join(out_dir, ph_dir('time'), "beast")
        return [os.path.join(bd, "merged", "merged.log"),
                os.path.join(bd, "beast1.log"),
                os.path.join(bd, "beast11.log")]

    def _beast_ran() -> bool:
        """是否真跑过后验 (链日志/合并产物非空)。--skip_beast_run 只出 XML → False。

        2026-10-13 恢复 (P1): 该闸门在 10-12 工作区重写中丢失 → --skip_beast_run
        写下 .done 后正式跑被 checkpoint 静默跳过, 形成"已跑过"假象。
        """
        _cands = _beast_log_candidates() + [os.path.join(
            out_dir, ph_dir('time'), "beast", "merged", "mcc.tree")]
        return any(os.path.exists(c) and os.path.getsize(c) > 0 for c in _cands)

    if 'beast' in stages:
        _beast_done = _done('beast')
        if _beast_done and not _beast_ran():
            # 2026-10-03 (P1-2): 旧 checkpoint 标志物只看 .done, 而 --skip_beast_run
            # (只出 XML 的用法) 也会写下它 → 正式跑被整段跳过、形成"已跑过"假象。
            log.warning("  [beast] checkpoint 命中但无任何后验日志/合并产物 "
                        "(只出过 XML 或产物被清理?) → 视为未跑, 重新执行")
            _beast_done = False
        if _beast_done:
            log.info("  [beast] checkpoint OK, skip"); _skip('beast')
            beast_xml = os.path.join(out_dir, ph_dir('time'), "beast", "beast1.xml")
            if os.path.exists(beast_xml):
                results["outputs"]["beast_xml"] = beast_xml
            # skip 时也恢复 TMRCA → metrics (供 gene_dating --full-tmrca)
            # 2026-10-11: 日志路径改为**多链感知** —— 原来硬编码 beast1.log,
            # 而多链模式 (beast_chains>1) 下只有 beast11/12/13.log, 该文件不存在
            # → TMRCA 在 skip 分支静默恢复失败 (PSTVd/GCVA 实测双双如此)。
            _bd = os.path.join(out_dir, ph_dir('time'), "beast")
            beast_log = None
            for _c in _beast_log_candidates():
                if not _usable_file(_c):
                    continue
                # 2026-10-03 (P1-4, 2026-10-13 恢复): 近空日志 (<10 行) 会被解析出
                # 不可信的 TMRCA 进汇总 —— 拒绝复用。
                _n_rows = count_beast_log_rows(_c)
                if _n_rows < 10:
                    log.warning(f"  [beast] 复用分支: {os.path.basename(_c)} 仅 "
                                f"{_n_rows} 行数据 → 跳过 TMRCA 复用 (不可用)")
                    continue
                beast_log = _c
                break
            if beast_log:
                try:
                    parsed = parse_beast_log(beast_log, burnin_pct=getattr(args, 'beast_burnin', 10),
                                              log=LogCollector(log))
                    rh = parsed.get("key_params", {}).get("treeModel.rootHeight")
                    if rh:
                        results["metrics"]["beast_tmrca"] = rh["median"]
                        results["metrics"]["beast_tmrca_hpd"] = f"{rh['hpd_95_lower']:.1f}–{rh['hpd_95_upper']:.1f}"
                        results["metrics"]["beast_tmrca_ess"] = rh["ess"]
                        log.info(f"  [beast] 全长 TMRCA (复用 {Path(beast_log).name}): "
                                 f"{rh['median']:.1f} yr "
                                 f"[{rh['hpd_95_lower']:.1f}–{rh['hpd_95_upper']:.1f}] ESS={rh['ess']:.0f}")
                except Exception as e:
                    log.warning(f"  [beast] skip 分支 TMRCA 解析失败: {e}")
            else:
                log.warning(f"  [beast] skip 分支未找到可用 BEAST log (查过 merged/merged.log, "
                            f"beast1.log, beast11.log) → TMRCA 无法复用")
            _restore('beast')             # 2026-10-11: 回填 outputs/metrics 快照
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
                _mark('beast', state={"outputs": beast_results.get("outputs"),
                                      "metrics": beast_results.get("metrics")})
            elif beast_results.get("beast_pending"):
                # 多链后台运行中: 不写 checkpoint, 下次跑检查合并
                _skip('beast', 'pending', '多链后台运行中, 下次跑检查合并')
            else:
                _fail('beast', beast_results.get("error", "failed"))

    # ── gene_dating: 分基因 BEAST 定年 ──
    if 'gene_dating' in stages:
        if _done('gene_dating'):
            log.info("  [gene_dating] checkpoint OK, skip"); _skip('gene_dating')
            _restore('gene_dating')       # 2026-10-11: 回填 outputs/metrics 快照
        else:
            _b = _snap()
            run_stage_gene_dating(prep, results, args, log)
            if results.get("outputs", {}).get("gene_dating_dir"):
                _mark('gene_dating', state=_delta(_b))
            else:
                _skip('gene_dating', 'blocked',
                      '无 capheine 前置产物 (先跑 select/capheine 或传 --genes_dir)')

    # ── treetime: TreeTime CLI 完整时间树 (phymapr 路线前置; 2026-10-12 新增) ──
    #    位置关键: 依赖 tree (与 clock 后的 gene_dating 同层), 必须在 phymapr 之前。
    if 'treetime' in stages:
        if not (aln_file and tree_file):
            _fail('treetime', '缺比对/树输入 (需 align/tree 产物)')
        elif _done('treetime'):
            log.info("  [treetime] checkpoint OK, skip"); _skip('treetime')
            _restore('treetime')          # 回填 outputs/metrics 快照
        else:
            _b = _snap()
            run_stage_treetime(prep, aln_file, tree_file, args, log, results)
            if results.get("outputs", {}).get("treetime_nexus"):
                _mark('treetime', state=_delta(_b))
            else:
                # 工具缺失 (phylo-treetime 未装) / 日期不足 / CLI 失败 → 优雅 skip:
                # 本 stage 默认关闭, 缺失不阻断主线 (装好后重跑即生效)。
                _skip('treetime', 'na')

    # ── phylogeo: 系统地理 ──
    if 'phylogeo' in stages:
        if _done('phylogeo'):
            log.info("  [phylogeo] checkpoint OK, skip"); _skip('phylogeo')
            _restore('phylogeo')          # 2026-10-11: 回填 outputs/metrics 快照
        else:
            _b = _snap()
            geo_pending = run_stage_phylogeo(prep, results, aln_file, args, log)
            if results.get("outputs", {}).get("phylogeo_xml"):
                if geo_pending:
                    # 链后台运行中: 不写 checkpoint, 下次跑检查合并 (同 beast pending 逻辑)
                    _skip('phylogeo', 'pending', '链后台运行中, 下次跑检查合并')
                else:
                    _mark('phylogeo', state=_delta(_b))
            else:
                _fail('phylogeo', '未生成 XML (可能缺数据)')
    # ── geo: 地理分析 (Mantel / 树地理)。位置关键：必须排在 phylogeo 之后 ——
    #    其 --rrt_randomized_dir 分支要吃 phylogeo 的 merged/mcc.tree ──
    if 'geo' in stages:
        if _done('geo'):
            log.info("  [geo] checkpoint OK, skip"); _skip('geo')
            _restore('geo')               # 2026-10-11: 回填 outputs/metrics 快照
        else:
            _b = _snap()
            run_stage_geo(prep, results, aln_file, tree_file, args, log)
            if results.get("outputs", {}).get("mantel_plot") or results.get("outputs", {}).get("geo_tree"):
                _mark('geo', state=_delta(_b))
            else:
                _skip('geo', 'blocked', 'Mantel/树地理均未产出')

    # ── geo_paths: 事件级传播路径 (2026-09-16 新增; 吃 phylogeo 的 mcc.tree) ──
    #    位置关键：**必须排在 phylogeo 之后**，否则首跑 mcc.tree 尚未生成 → 整段跳过。
    #    2026-09-18: 原先误排在 phylogeo 之前（STAGE_DEPS 声明了依赖却没人执行它）。
    if 'geo_paths' in stages:
        if _done('geo_paths'):
            log.info("  [geo_paths] checkpoint OK, skip"); _skip('geo_paths')
            _restore('geo_paths')         # 2026-10-11: 回填 outputs/metrics 快照
        else:
            _b = _snap()
            run_stage_geo_paths(prep, results, aln_file, tree_file, args, log)
            if results.get("outputs", {}).get("pathway_qc"):
                _mark('geo_paths', state=_delta(_b))
            else:
                _skip('geo_paths', 'blocked',
                      'pathway_qc 未产出 (常见: 缺 phylogeo merged/mcc.tree)')

    # ── phymapr: 系统地理 TreeTime 路线 (BM-ASR + 传播地图; 2026-10-12 新增) ──
    #    位置关键: 吃 treetime 的 timetree.nexus, **不依赖 BEAST/phylogeo** ——
    #    与 geo_paths (BEAST 路线) 并列; 排在 geo_paths 之后保持与 STAGE_ORDER 一致。
    if 'phymapr' in stages:
        if _done('phymapr'):
            log.info("  [phymapr] checkpoint OK, skip"); _skip('phymapr')
            _restore('phymapr')           # 回填 outputs/metrics 快照
        else:
            _b = _snap()
            run_stage_phymapr(prep, results, args, tools, log)
            if results.get("outputs", {}).get("phymapr_summary"):
                _mark('phymapr', state=_delta(_b))
            else:
                # R/依赖缺失 或 timetree 未产出 → 优雅 skip (不 fail)
                _skip('phymapr', 'na')


    # 2026-10-13 恢复 (P0-1): 未被任何块记账的 stage 兜底记 blocked。
    # 触发场景: `if 'splitstree' in stages and aln_file:` —— align 未产出时整块不执行,
    # 旧版既不 done 也不 skipped/failed, 在三个清单里全部消失 (静默)。
    # 注: 'report' 由 main 末尾统一执行 (不按病毒), 不参与单病毒记账。
    _accounted = set(st["done"]) | set(st["skipped"]) | set(st["failed"])
    for _s in stages:
        if _s in ('prep', 'report') or _s in _accounted:
            continue
        _skip(_s, 'blocked', '未执行: 前置产物缺失')

    results["stages_done"] = st["done"]
    results["stages_skipped"] = st["skipped"]
    results["stages_failed"] = st["failed"]
    results["stages_blocked"] = st["blocked"]
    results["stages_skip_reasons"] = st["skip_reasons"]
    # incomplete = 有受阻 stage 或仍有链在后台 (结果不完整/未定)
    results["incomplete"] = bool(st["blocked"]) or any(
        r == 'pending' for r in st["skip_reasons"].values())
    # 2026-09-15 修复 (P1-B): 旧版 `success = bool(tree_file) or bool(outputs)`,
    # 只要任一 stage 往 results["outputs"] 写过东西就算成功 —— 于是 align/tree
    # 双双失败、仅 online 阶段写出 ncbi_fasta 时也会被判 success=True, 汇总表
    # 把"只跑完可选的序列补充"报成通过。现要求「无任何 stage 失败」且「有实质产物」。
    results["success"] = (not st["failed"]) and bool(
        tree_file or aln_file or results.get("outputs"))
    log.info(f"  [{vname}] stages: done={st['done']} skipped={st['skipped']} failed={st['failed']}")
    if st["blocked"]:
        log.warning(f"  ⚠ [{vname}] 受阻 stage (结果不完整, 勿直接引用): {st['blocked']} — "
                    f"原因: { {k: v for k, v in st['skip_reasons'].items() if k in st['blocked'] or v in ('blocked', 'degraded')} }; "
                    f"修好前置后重跑")
    if any(r == 'pending' for r in st["skip_reasons"].values()):
        log.warning(f"  ⏳ [{vname}] 有后台链条运行中 (pending) — 结果未定, 下次重跑自动合并")
    # 2026-10-11: 落盘整份 results —— 供 main() 在"只跑部分 stage"时重建完整汇总
    # (修 P0 "重跑一次 report 就得到空汇总")。
    save_stage_results(out_dir, results)
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
            # 2026-10-03 (P1-1): 结果完整性 —— Incomplete=True 表示有受阻 stage 或链仍在
            # 后台 (Success 只代表"无硬失败", 不再等价于"结果完整")。
            # 2026-10-13 恢复: 两列在 10-12 工作区重写中随记账机制一并丢失。
            "Incomplete": r.get("incomplete", False),
            "Stages_done": ";".join(r.get("stages_done", [])),
            "Stages_skipped": ";".join(r.get("stages_skipped", [])),
            "Stages_failed": ";".join(r.get("stages_failed", [])),
            "Stages_blocked": ";".join(r.get("stages_blocked", [])),
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

    # ── 分裂网络 (SplitsTree6; 2026-10-12 补登记 —— RUN_GUIDE/STAGE_REFERENCE 已承诺) ──
    g_ss = parser.add_argument_group("Split network (SplitsTree6)")
    g_ss.add_argument("--splitstree_max_taxa", type=int,
                      default=defaults.get("splitstree_max_taxa", 200),
                      help="NeighborNet 最多取样 tip 数 (默认 200)")
    g_ss.add_argument("--splitstree_group_col",
                      # ⚠ 2026-10-13 默认由 "auto" 改为 "host"：用户指定
                      #   「改成宿主着色」（与 ML 树口径统一，见 stage_plots
                      #   ::plot_tree_circular 的 color_by 默认值）。
                      #   ⚠ "auto" 是依次试 location → host，本数据两张元数据表
                      #   各缺一列，auto 会命中 location → 仍是地区着色。
                      default=defaults.get("splitstree_group_col", "host"),
                      choices=["auto", "location", "host", "none"],
                      help="NeighborNet 分组列 (默认 host=宿主着色; auto=location→host 回退; "
                           "none=只画网络); 用于 tip 着色 + 外圈分组弧线 + 组名")
    g_ss.add_argument("--splitstree_layout",
                      default=defaults.get("splitstree_layout", "native"),
                      choices=["native", "equal_angle"],
                      help="neighbornet.pdf/png 的布局来源: native=SplitsTree6 自己的 "
                           "SplitNetworkLayout (默认, 需 splitstree6-tools); "
                           "equal_angle=旧版 Python 自绘 equal-angle (回退用)")
    g_ss.add_argument("--splitstree_unweight",
                      default=defaults.get("splitstree_unweight", "on"),
                      choices=["on", "off"],
                      help="native 布局是否用**等长**几何 (on=默认). on 时把每条边的"
                           "权重当 1 重放 DFS (源码 EqualAngle), 得到所有边长一致的"
                           "规则网格 —— 与用户参考图 (SplitsTree GUI 截图) 一致; "
                           "线宽同时自动切等粗 (flat), 避免用粗细重复表达枝长. "
                           "off=旧加权布局 (边长跨 344 倍, 放射状长枝).")
    g_ss.add_argument("--splitstree_trim",
                      default=defaults.get("splitstree_trim", "gt"),
                      choices=["gt", "auto", "nogaps", "none"],
                      help="建网前用 trimAl 修剪比对 (默认 gt=按 gap 阈值剪). "
                           "gt=保留 gap<=阈值 的列; auto=-automated1; nogaps=-nogaps; "
                           "none=用未修剪比对. ⚠ 未修剪时高 gap 列会让 NeighborNet "
                           "中心堆成一团乱网 (实测), 且布局在 449/514 两种实现间随机跳。")
    g_ss.add_argument("--splitstree_trim_gap", type=float,
                      default=defaults.get("splitstree_trim_gap", 0.1),
                      help="--splitstree_trim gt 时的 gap 阈值 (0-1, 默认 0.1). "
                           "实测 0.1 保留 90.9%% 列、aspect 稳定 1.63; "
                           "0.2 接近; auto/nogaps 会过度修剪退化成三叉星芒。")
    g_ss.add_argument("--splitstree_label_rings", type=int,
                      default=defaults.get("splitstree_label_rings", 9),
                      help="tip 标签的外围环数 (默认 9). 145 个标签在中心密集区会互压, "
                           "推到外围环上 + 引线; 实测 9 环 + 环间避让 = 0 对重叠。")

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

    # ── TreeTime 路线 (2026-10-12 新增; STAGE_REFERENCE §19/§20) ──
    # 与 BEAST 路线并列的第二条时间-地理方法: treetime (CLI 完整时间树) →
    # phymapr (R BM-ASR 祖先重建 + 传播地图)。两条默认关闭 (config enabled: false),
    # 开启方式: config 改 true / 命令行开关 / --stage 点名 (点名豁免门控)。
    g_tt = parser.add_argument_group("TreeTime route (Phase 6T)")
    g_tt.add_argument("--treetime", action="store_true", default=defaults.get("treetime", False),
                      help="运行 TreeTime CLI 完整时间树 → time/treetime/timetree.nexus")
    g_tt.add_argument("--treetime_clock_filter", type=float,
                      default=defaults.get("treetime_clock_filter", 4.0),
                      help="IQD 离群 tip 剔除阈值 (0=关; 默认 4.0)")
    g_tt.add_argument("--treetime_coalescent", default=defaults.get("treetime_coalescent", "skyline"),
                      choices=["const", "skyline", "optimal", "0.0"],
                      help="coalescent 先验 (默认 skyline)")
    g_tt.add_argument("--treetime_clock_rate", type=float,
                      default=defaults.get("treetime_clock_rate", None),
                      help="手动钟速率 (慢进化病原用; 默认 None=TreeTime 回归自估)")
    g_tt.add_argument("--phymapr", action="store_true", default=defaults.get("phymapr", False),
                      help="运行 phymapr R 系统地理 (吃 treetime 的 timetree.nexus, 不依赖 BEAST)")
    g_tt.add_argument("--phymapr_map_style", default=defaults.get("phymapr_map_style", "polygon"),
                      choices=["polygon", "tile"],
                      help="地图样式 (polygon=离线矢量, 默认; tile=暗色瓦片需联网)")
    g_tt.add_argument("--phymapr_n_prune_early", type=int,
                      default=defaults.get("phymapr_n_prune_early", 4),
                      help="动画丢弃最早 N 条 pathway (默认 4)")
    g_tt.add_argument("--phymapr_direct_snp", type=float,
                      default=defaults.get("phymapr_direct_snp", 2.0),
                      help="Direct 置信阈值 (phymapr 内部 snp_distance 语义; 默认 2.0)")
    g_tt.add_argument("--phymapr_indirect_snp", type=float,
                      default=defaults.get("phymapr_indirect_snp", 5.0),
                      help="Indirect 置信阈值 (默认 5.0)")
    g_tt.add_argument("--phymapr_distant_style",
                      default=defaults.get("phymapr_distant_style", "import_arrow"),
                      choices=["import_arrow", "pulse", "curve", "hide"],
                      help="distant 路径渲染样式 (默认 import_arrow)")
    g_tt.add_argument("--phymapr_pace", type=float,
                      default=defaults.get("phymapr_pace", 0.5),
                      help="动画节奏: 1=纯日历时间 0=事件均速 (默认 0.5)")

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
    # 2026-10-12: treetime/phymapr (TreeTime 路线) 并入门控 —— config 里两者
    # enabled: false, 同 rdp5/beast/gene_dating/phylogeo 一套门控语义。
    _enable_dests = ('rdp5', 'beast', 'gene_dating', 'phylogeo',
                     'treetime', 'phymapr')

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
        # ── 依赖保护 (2026-10-12) ──
        # 门控不得把**被保留 stage 的依赖**删掉: 例 `--stage geography` 点名
        # phymapr (豁免门控), 而 treetime 的 config enabled=false 会被移出清单 ——
        # 于是 phymapr 跑起来却拿不到 timetree.nexus, 整条 TreeTime 路线形同虚设。
        # 被保留 stage 的依赖闭包一律保留 (只是"随依赖被带上", 不算用户点名)。
        _kept_deps = set()
        _stack = [s for s in stages if s not in _cfg_off]
        while _stack:
            _cur = _stack.pop()
            for _d in STAGE_DEPS.get(_cur, ()):
                if _d in _cfg_off and _d not in _kept_deps:
                    _kept_deps.add(_d)
                    _stack.append(_d)
        _removed = [s for s in _cfg_off if s not in _kept_deps]
        if _kept_deps:
            print(f"[stage] 依赖保护: {', '.join(sorted(_kept_deps))} 虽 enabled=false, "
                  f"但被保留 stage 依赖 → 随依赖保留")
        stages = [s for s in stages if s not in _removed]
        if _removed:
            print(f"[stage] config 里 {', '.join(_removed)}.enabled=false -> 已移出执行清单。"
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
        virus_name = _derive_virus_name(data_dir, ds)
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
        # 2026-10-11 修 (P0): 「只跑 report」时 process_virus 的 results 只含本次
        # 请求的 stage → 汇总表/HTML 退化成空壳 (实测 PSTVd 2026-09-03 00:44 与
        # 2026-09-22 23:13 两次 report-only 重跑: Beta/R²/P/TMRCA 全空)。
        # 这里把磁盘上的 stage_results.json + .checkpoints/*.state.json 回填进来,
        # 使"重跑一次 report"与"全量跑"产出等价。
        results = merge_persisted_results(results, data_dir)
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
        # 2026-10-13 修 (P0-3 同类): datasets.yaml 的 fasta 是未比对原料, AlignIO.read
        # 对长度不齐输入直接抛 ValueError → 优先用已产出的 MAFFT 比对。
        _aln_r = os.path.join(data_dir, ph_dir('phylogeny'), "mafft.aln.fasta")
        _rtt_fa = (_aln_r if (os.path.exists(_aln_r) and os.path.getsize(_aln_r) > 0)
                   else fasta)
        _pf_cmd = ["--rtt_fasta", _rtt_fa,
                        "--rtt_dates", _pp["dates"],
                        "--rtt_tree", tree_file or ""]
        if _beast_log:
            _pf_cmd += ["--beast_log", _beast_log]
        if os.path.exists(_mcc_tree):
            _pf_cmd += ["--mcc_tree", _mcc_tree]
        # 2026-10-11: 把病毒名传给 publication_figures —— 其 Fig1/Fig2 标题原先
        # **硬编码 "GCVA"**, 跑 PSTVd 时两张发表级图的标题都写成 GCVA。
        # 优先用 datasets.yaml 的 species (可读名, 如 "Potato spindle tuber viroid");
        # 没有则退化为规范名 (目录名, 如 Potato_spindle_tuber_viroid_NC_002030.1)。
        _pf_cmd += ["--virus_name", (prep.get("species") or prep.get("virus_name", "")),
                    "--output_dir", os.path.join(report_dir, "publication")]
        run_postprocess("publication_figures", logger, *_pf_cmd)
        # 后处理: 采样分布图 (时间线 + 地理地图)
        # 2026-10-13: 追加宿主表 → 时间分布图的气泡内切宿主饼 (用户「图 B 把物种
        #   也加上」)。`_governed_meta` 已带 host 通道 (无 host.csv 时为空串,
        #   sample_plot 会自动退回纯地点着色圆)。
        _sp_cmd = ["--metadata", _pp["meta"]]
        if _pp.get("host"):
            _sp_cmd += ["--host-metadata", _pp["host"]]
        run_postprocess("sample_plot", logger, *_sp_cmd,
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
            # --alignment 是用户给的**比对文件**: 允许 align 早退复用 (2026-10-13 P0)
            "provided_alignment": True,
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

        report_dir = os.path.join(str(work_dir), ph_dir('report'))
        os.makedirs(report_dir, exist_ok=True)
        # 2026-10-11 修 (P0): 同 --virus 分支 —— 回填磁盘状态, 避免 report-only 空壳
        results = merge_persisted_results(results, str(work_dir))
        if getattr(args, "report_enabled", True):
            generate_summary([results], report_dir, logger)
            if results:
                report_path = build_html_report([results], report_dir)
                logger.info(f"  HTML report: {report_path}")
        else:
            logger.info("  [report] report.enabled=false, 跳过 phylo_summary.csv / HTML")
        # 后处理: 发表级三件套 (同单病毒模式修复: beast_log 多链探测 + 补 mcc_tree)
        # 2026-10-13 修 (P1): 硬编码 "time" 与 batch 分支 ph_dir('time') 口径不一致,
        # standard 布局下探测不到 beast 日志/MCC → Table S1/Fig2/3 面板静默缺失。
        _beast_dir_w = os.path.join(str(work_dir), ph_dir('time'), "beast")
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
        # 2026-10-11: 同 --virus 分支 —— 传病毒名, 避免图题硬编码 GCVA。
        _pf_cmd_w += ["--virus_name", (prep.get("species") or prep.get("virus_name", "")),
                      "--output_dir", os.path.join(report_dir, "publication")]
        run_postprocess("publication_figures", logger, *_pf_cmd_w)
        # 后处理: 采样分布图
        _sp_cmd_w = ["--metadata", _pp_w["meta"]]
        if _pp_w.get("host"):
            _sp_cmd_w += ["--host-metadata", _pp_w["host"]]
        run_postprocess("sample_plot", logger, *_sp_cmd_w,
                        "--outdir", os.path.join(report_dir, "sample_plot"))
        logger.info(f"\nDone. Results: {work_dir}")
        # 退出码 (2026-10-03, 2026-10-13 恢复): work-dir 分支有独立出口,
        # 同样要按记账结果定退出码 (failed→1 / blocked·pending→2 / 干净→0)
        _wcode, _wwhy = pipeline_exit_code([results])
        if _wcode:
            logger.warning(f"EXIT {_wcode}: {_wwhy}")
        raise SystemExit(_wcode)
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
    # 2026-10-12 修 (P2-13): 目录名原先写死 "report", 与单病毒/work-dir 分支的
    # ph_dir('report') 口径不一致 —— standard 布局下汇总会落到 05_Phylo/report,
    # 而 per-virus 产物在 08_report/, 报告与产物分家。统一走 ph_dir()。
    report_dir = os.path.join(str(output_dir), ph_dir('report'))
    os.makedirs(report_dir, exist_ok=True)
    # 2026-10-11 修 (P0): 局部重跑时用磁盘状态回填, 避免汇总退化成"只有本次 stage"
    all_results = [merge_persisted_results(r, r.get("virus_dir", ""))
                   for r in all_results]
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
    # 2026-10-13 修 (P0): publication_figures 必须吃 **比对后序列与树** ——
    # prep.combined_fasta 已被 clean 就地改写成未比对 clean.fasta (AlignIO.read
    # 对长度不齐输入直接抛 ValueError, Fig2/Fig3/TableS1 全灭), 且 prep 没有
    # tree_file 键 (恒 "")。现按病毒名取记账结果里 align/tree 的真实产物。
    _res_by_virus = {r.get("virus"): r for r in all_results}
    for vq in virus_queue:
        vname = vq["virus_name"]
        vdir = os.path.join(str(output_dir), vname)
        v_meta = vq.get("meta_csv") or vq.get("dates_csv", "")
        v_report = os.path.join(vdir, ph_dir('report'))
        _vout = (_res_by_virus.get(vname) or {}).get("outputs", {})
        # 采样分布图 (时间线 + 地理地图)
        if v_meta and os.path.exists(v_meta):
            # 2026-10-13: 同单病毒分支 —— 有宿主表就传给 sample_plot 画宿主饼。
            _sp_cmd_b = ["--metadata", v_meta]
            _vhost = vq.get("host_csv") or ""
            if _vhost and os.path.exists(_vhost):
                _sp_cmd_b += ["--host-metadata", _vhost]
            run_postprocess("sample_plot", logger, *_sp_cmd_b,
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
        _cmd = ["--rtt_fasta", _vout.get("alignment") or vq.get("combined_fasta", ""),
                "--rtt_dates", v_meta,
                "--rtt_tree", _vout.get("tree") or vq.get("tree_file", "")]
        if _bl:
            _cmd += ["--beast_log", _bl]
        if os.path.exists(_mcc):
            _cmd += ["--mcc_tree", _mcc]
        # 2026-10-11: batch 模式同样传病毒名, 避免图题硬编码 GCVA
        # (优先 species, 与单病毒分支口径一致)。
        _cmd += ["--virus_name", (vq.get("species") or vname),
                 "--output_dir", os.path.join(v_report, "publication")]
        run_postprocess("publication_figures", logger, *_cmd)

    logger.info(f"\n{'='*60}")
    logger.info(f"PhyloPipeline complete!")
    logger.info(f"Output: {output_dir}")
    logger.info(f"{'='*60}")

    # 退出码 (2026-10-03, P1-1; 2026-10-13 恢复): 旧版无论有无失败都退出 0 →
    # 自动化/&&链/per-virus driver 一律"成功"。现: 有 failed → 1; 无失败但有
    # blocked/pending (结果不完整/未定) → 2; 全部干净 → 0。
    _n_fail = sum(len(r.get("stages_failed", [])) for r in all_results)
    _n_blocked = sum(len(r.get("stages_blocked", [])) for r in all_results)
    _n_pending = sum(1 for r in all_results if r.get("incomplete") and
                     not r.get("stages_blocked"))
    if _n_fail:
        logger.warning(f"EXIT 1: {_n_fail} 个 stage 失败 — 见上表 Stages_failed")
        raise SystemExit(1)
    if _n_blocked or _n_pending:
        logger.warning(f"EXIT 2: 结果不完整 (受阻 stage {_n_blocked} 个"
                       f"{f', 另有 {_n_pending} 个病毒有 pending 链' if _n_pending else ''})"
                       f" — 修好前置后重跑")
        raise SystemExit(2)
    raise SystemExit(0)


if __name__ == "__main__":
    main()
