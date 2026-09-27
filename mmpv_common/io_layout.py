#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
io_layout.py — MMPV 五管线统一 I/O 目录布局 (legacy / standard 双布局)
=======================================================================

MMPV 的 ①公共数据 ②数据清洗 ③病毒发现 ④已知病毒分析 ⑤系统发育 五条管线
各自维护一套输出目录名, 边界交接靠人记路径。本模块把「逻辑目录键 → 相对
目录名」收敛到一处, 提供两套布局:

  legacy   现行目录名 (00a_CleanData / 02a_Identification / 09a_Virome_Analysis
           / 01_detection ...), 与 v3.0 及之前逐字节一致 —— 默认值,
           已有 checkpoint / 正在跑的任务零影响。
  standard 管线独立根布局 (设计见 doc/IO_LAYOUT_DESIGN.md):
           每条管线一个独立输出根 (01_PublicData / 02_Preprocessing /
           03_Discovery / 04_Analysis / 05_Phylo), 根内部各自从 01 独立编号;
           编排器的 --output_dir/--work-dir 即指向本管线自己的根。
           跨管线交接物统一登记在项目根的 90_Handoff/。

布局解析优先级:  CLI --io-layout  >  环境变量 MMPV_IO_LAYOUT  >  legacy。
环境变量同时是向子进程 (② 被 ③ 调用等) 传播布局的机制 —— 编排器入口
统一调用 normalize_layout_env() 把解析结果写回环境。

边界产物定位 (跨布局): locate_boundary() 按「布局各目录名 → 产物相对路径」
扫描项目根, 供桥脚本 (discovery2analysis --from-discovery) 与文档引用。

纯标准库, 无第三方依赖; 所有编排器可安全 import。
"""

import os
from pathlib import Path
from typing import Dict, List, Optional

LAYOUT_ENV = "MMPV_IO_LAYOUT"
LAYOUTS = ("legacy", "standard")
DEFAULT_LAYOUT = "legacy"

# ═══════════════════════════════════════════════════════════════════
# 目录名注册表: 逻辑键 → 相对路径 (相对各管线的输出根目录)
#
# 键名约定 (跨管线不重号):
#   <pipeline prefix>  d_  = ②③ 共用发现链;  a_  = ④ 分析;  m_  = ① 元数据;
#                      h_  = ① 宿主参考;   p_  = ⑤ 进化根;   x_  = 交接层
# ═══════════════════════════════════════════════════════════════════

_DIR_TABLE: Dict[str, Dict[str, str]] = {
    # ── legacy: 与 v3.0 逐字节一致 (含 ③ 的字母后缀体系) ──
    "legacy": {
        # ② 数据清洗 (③ 的 clean/deplete 阶段复用同一对目录)
        "d_clean":        "00a_CleanData",
        "d_hostdep":      "00b_HostDepletion",
        "d_bbnorm":       "00c_BBnorm",
        # ③ 病毒发现
        "d_asm":          "01_Assembly",
        "d_ident":        "02a_Identification",
        "d_filter":       "02b_Filter",
        "d_cobra":        "03a_COBRA",
        "d_merge":        "03b_MergeSamples",
        "d_cluster":      "04_CLUSTER",
        "d_centroids":    "04_CLUSTER/4_centroids",
        "d_taxonomy":     "05_Taxonomy",
        "d_host_pred":    "06_HostPrediction",
        "d_checkv":       "07_Checkv",
        "d_rescue":       "08_Rescue",
        "d_analysis":     "09a_Virome_Analysis",
        "d_verify":       "09b_Analysis_Verify",
        "d_reports":      "10_Reports",
        # ④ 已知病毒分析
        "a_detect":       "01_detection",
        "a_filter":       "02_filtering",
        "a_variants":     "03_variants",
        "a_post":         "04_post_analysis",
        "a_assembly":     "05_assembly",
        "a_extract":      "06_extraction",
        "a_similarity":   "07_similarity",
        "a_dvg":          "08_dvg",
        "a_report":       "09_report",
        # ① 公共数据
        "m_search":       "search",
        "m_info":         "info",
        "m_plot":         "plot",
        "m_down":         "down",
        "h_genome":       "host_reference/genome",
        "h_hostdb":       "host_reference/hostdb",
        # ⑤ 进化 (默认输出根, 仅作登记; ⑤ 内部 per-virus 结构两布局一致)
        "p_root":         "phylo_results",
        "ph_data":        "data",
        "ph_phylogeny":   "phylogeny",
        "ph_recomb":      "recomb",
        "ph_popgen":      "popgen",
        "ph_select":      "select",
        "ph_time":        "time",
        "ph_geography":   "geography",
        "ph_report":      "report",
        # ⑥ EVE (内部编号两布局一致)
        "e_loci":         "01_Loci",
        "e_verdict":      "02_Verdict",
        "e_rvdb":         "03_RVDB",
        "e_summary":      "04_Summary",
        # 交接层
        "x_handoff":      "90_Handoff",
    },

    # ── standard: 管线独立根 + 根内独立编号 (doc/IO_LAYOUT_DESIGN.md) ──
    # 相对路径均相对「本管线自己的输出根」(编排器 --output_dir/--work-dir 指向它):
    #   ① --work-dir    = <项目>/01_PublicData
    #   ② --output_dir  = <项目>/02_Preprocessing
    #   ③ --output_dir  = <项目>/03_Discovery
    #   ④ --output_dir  = <项目>/04_Analysis
    #   ⑤ --output_dir  = <项目>/05_Phylo
    "standard": {
        # ① 公共数据 (根: 01_PublicData, 按 stage 顺序编号)
        "m_search":       "01_Search",
        "m_info":         "02_Unified",
        "m_down":         "03_RawData",
        "m_plot":         "04_Plots",
        "h_genome":       "05_HostRef/genome",
        "h_hostdb":       "05_HostRef/hostdb",
        # ② 数据清洗 (根: 02_Preprocessing)
        "d_clean":        "01_CleanData",
        "d_hostdep":      "02_HostDepleted",
        "d_bbnorm":       "03_BBnorm",
        # ③ 病毒发现 (根: 03_Discovery, 15 阶段独立编号)
        "d_asm":          "01_Assembly",
        "d_ident":        "02_Identification",
        "d_filter":        "03_Filter",
        "d_cobra":        "04_COBRA",
        "d_merge":        "05_CoAssembly",
        "d_cluster":      "06_CLUSTER",
        "d_centroids":    "06_CLUSTER/centroids",
        "d_taxonomy":     "07_Taxonomy",
        "d_host_pred":    "08_HostPrediction",
        "d_checkv":       "09_CheckV",
        "d_rescue":       "10_Rescue",
        "d_analysis":     "11_ViromeAnalysis",
        "d_verify":       "12_AnalysisVerify",
        "d_reports":      "13_Report",
        # ④ 已知病毒分析 (根: 04_Analysis, 9 阶段独立编号)
        "a_detect":       "01_Detection",
        "a_filter":       "02_Filtering",
        "a_variants":     "03_Variants",
        "a_post":         "04_PostAnalysis",
        "a_assembly":     "05_Assembly",
        "a_extract":      "06_Extraction",
        "a_similarity":   "07_Similarity",
        "a_dvg":          "08_DVG",
        "a_report":       "09_Report",
        # ⑤ 进化 (--output_dir 指向本键; per-virus 内部模块按输出顺序编号)
        "p_root":         "05_Phylo",
        "ph_data":        "01_data",
        "ph_phylogeny":   "02_phylogeny",
        "ph_popgen":      "03_popgen",
        "ph_recomb":      "04_recomb",
        "ph_select":      "05_select",
        "ph_time":        "06_time",
        "ph_geography":   "07_geography",
        "ph_report":      "08_report",
        # ⑥ EVE (输出根约定 06_EVE; 内部编号两布局一致, 本就符合编号惯例)
        "e_loci":         "01_Loci",
        "e_verdict":      "02_Verdict",
        "e_rvdb":         "03_RVDB",
        "e_summary":      "04_Summary",
        # 交接层 (相对项目根, 跨管线共有)
        "x_handoff":      "90_Handoff",
    },
}

# ③→④ 边界产物: 相对各布局键目录的固定文件名
CENTROIDS_BASENAME = "final_centroids.fasta"
TAXONOMY_RELPATH = "integrated/final_integrated_classification.tsv"

# ③ centroids 目录的历史名 (legacy 内部三代命名, 定位时全部尝试)
_LEGACY_CENTROID_DIRS = ("04_CLUSTER/4_centroids",
                         "04_CLUSTER/04_centroids",
                         "04_CLUSTER/4.centroids")


def resolve_layout(explicit: Optional[str] = None) -> str:
    """布局名解析: 显式参数 > 环境变量 > legacy。非法值报错而非静默回退。"""
    name = explicit or os.environ.get(LAYOUT_ENV) or DEFAULT_LAYOUT
    if name not in LAYOUTS:
        raise ValueError(
            f"未知 I/O 布局: {name!r} (可选: {', '.join(LAYOUTS)}; "
            f"环境变量 {LAYOUT_ENV})")
    return name


def normalize_layout_env(explicit: Optional[str] = None) -> str:
    """解析布局并把结果写回环境变量, 返回布局名。

    编排器入口必须调用本函数: 一处解析, 子进程 (② 被 ③ 调用等) 经环境
    继承同一布局, 不需要每个下游脚本都加 CLI 参数。
    """
    name = resolve_layout(explicit)
    os.environ[LAYOUT_ENV] = name
    return name


def layout_dirs(layout: Optional[str] = None) -> Dict[str, str]:
    """返回布局的 {逻辑键: 相对目录名} 只读副本。"""
    return dict(_DIR_TABLE[resolve_layout(layout)])


def dir_name(key: str, layout: Optional[str] = None) -> str:
    """取单个逻辑键的目录名。"""
    try:
        return _DIR_TABLE[resolve_layout(layout)][key]
    except KeyError:
        raise KeyError(f"io_layout: 未知目录键 {key!r} (布局 {resolve_layout(layout)})")


# ═══════════════════════════════════════════════════════════════════
# 各编排器的目录字典构造器
# ═══════════════════════════════════════════════════════════════════

def build_discovery_dirs(out: Path, raw: Path, layout: Optional[str] = None) -> Dict:
    """③ virome_pipeline.py 的 self.d 目录字典。

    legacy 分支逐键复刻 v3.0 字面量 (含 _ 下划线回退键), 保证既有
    checkpoint / 结果目录解析行为不变。
    """
    name = resolve_layout(layout)
    if name == "legacy":
        return {
            'raw':       raw,
            'root':      out,
            'clean':     out / '00a_CleanData',
            'hostdep':   out / '00b_HostDepletion',
            'asm':       out / '01_Assembly',
            'ident':     out / '02a_Identification',
            'filter':    out / '02b_Filter',
            'cobra':         out / '03_COBRA',
            '_cobra_dirs':   [out / '03a_COBRA', out / '03_COBRA'],
            'merge_samples': out / '03b_MergeSamples',
            'cluster':       out / '04_CLUSTER',
            'centroids':     out / '04_CLUSTER' / '4_centroids',   # 新集群用 4_centroids
            '_centroids_v1': out / '04_CLUSTER' / '04_centroids',  # 旧回退
            '_centroids_v2': out / '04_CLUSTER' / '4.centroids',   # 旧回退
            'taxonomy':  out / '05_Taxonomy',
            'host_pred':  out / '06_HostPrediction',
            'checkv_dir': out / '07_Checkv',
            'rescue_dir': out / '08_Rescue',
            'analysis':  out / '09_Virome_Analysis',
            # 阶段目录命名统一为 NN[ab]_Name（00a/00b、03a/03b、09b），09_ 是唯一
            # 例外。改名为 09a_Virome_Analysis 与 09b_ 配对时，下游按这张表解析
            # 新旧两个名字（同 _cobra_dirs / _centroids_v1 的写法）。
            '_analysis_dirs': [out / '09a_Virome_Analysis',
                               out / '09_Virome_Analysis'],
            'analysis_verify': out / '09b_Analysis_Verify',
            'reports':    out / '10_Reports',
        }
    d = _DIR_TABLE["standard"]
    return {
        'raw':       raw,
        'root':      out,
        'clean':     out / d['d_clean'],
        'hostdep':   out / d['d_hostdep'],
        'asm':       out / d['d_asm'],
        'ident':     out / d['d_ident'],
        'filter':    out / d['d_filter'],
        'cobra':         out / d['d_cobra'],
        '_cobra_dirs':   [out / d['d_cobra']],
        'merge_samples': out / d['d_merge'],
        'cluster':       out / d['d_cluster'],
        'centroids':     out / d['d_centroids'],
        '_centroids_v1': out / d['d_centroids'],   # standard 无历史命名, 回退键同指
        '_centroids_v2': out / d['d_centroids'],
        'taxonomy':  out / d['d_taxonomy'],
        'host_pred':  out / d['d_host_pred'],
        'checkv_dir': out / d['d_checkv'],
        'rescue_dir': out / d['d_rescue'],
        'analysis':  out / d['d_analysis'],
        '_analysis_dirs': [out / d['d_analysis']],
        'analysis_verify': out / d['d_verify'],
        'reports':    out / d['d_reports'],
    }


def build_analysis_dirs(out: Path, layout: Optional[str] = None) -> Dict[str, Path]:
    """④ auto_known_virus.py 的 9 个阶段输出目录 {键: Path}。"""
    d = _DIR_TABLE[resolve_layout(layout)]
    return {
        'detect':     out / d['a_detect'],
        'filter':     out / d['a_filter'],
        'variants':   out / d['a_variants'],
        'post':       out / d['a_post'],
        'assembly':   out / d['a_assembly'],
        'extract':    out / d['a_extract'],
        'similarity': out / d['a_similarity'],
        'dvg':        out / d['a_dvg'],
        'report':     out / d['a_report'],
    }


def build_meta_dirs(work_dir: Path, layout: Optional[str] = None) -> Dict[str, str]:
    """① public_data_pipeline.py 的 search/info/down/plot 目录 (str 路径)。"""
    d = _DIR_TABLE[resolve_layout(layout)]
    return {
        'search': str(work_dir / d['m_search']),
        'info':   str(work_dir / d['m_info']),
        'plot':   str(work_dir / d['m_plot']),
        'down':   str(work_dir / d['m_down']),
    }


# ⑤ per-virus 模块 (编号顺序 = 主线分析顺序, 与 STAGE_REFERENCE 8 大模块表一致:
#   data[prep→clean] → phylogeny[align,splitstree,tree] → popgen[popgen,host] →
#   recomb[rdp5] → select[capheine] → time[clock,beast,gene_dating] →
#   geography[phylogeo,geo,geo_paths] → report)
PHYLO_MODULES = ("data", "phylogeny", "popgen", "recomb",
                 "select", "time", "geography", "report")


def ph_dir(module: str, layout: Optional[str] = None) -> str:
    """⑤ 单个 per-virus 模块目录名 (布局感知)。

    用法: os.path.join(work, ph_dir('phylogeny'), 'mafft.aln.fasta')
    """
    if module not in PHYLO_MODULES:
        raise KeyError(f"io_layout: 未知 ⑤ 模块 {module!r} (可选: {', '.join(PHYLO_MODULES)})")
    return dir_name('ph_' + module, layout)


def build_phylo_dirs(work_dir: Path, layout: Optional[str] = None) -> Dict[str, Path]:
    """⑤ per-virus 工作目录的 8 个模块子目录 {模块名: Path}。"""
    d = _DIR_TABLE[resolve_layout(layout)]
    return {m: Path(work_dir) / d['ph_' + m] for m in PHYLO_MODULES}


def build_eve_dirs(out_dir: Path, layout: Optional[str] = None) -> Dict[str, Path]:
    """⑥ eve_screen.py 的四阶段输出目录 {阶段: Path} (两布局同名, 已编号)。"""
    d = _DIR_TABLE[resolve_layout(layout)]
    return {k: Path(out_dir) / d['e_' + k]
            for k in ("loci", "verdict", "rvdb", "summary")}


# ═══════════════════════════════════════════════════════════════════
# 边界产物定位与登记 (桥脚本 / 文档共用)
# ═══════════════════════════════════════════════════════════════════

def locate_centroids(root: Path) -> Optional[Path]:
    """在发现管线输出根下定位 final_centroids.fasta (两布局 + 三代旧名)。"""
    root = Path(root)
    cands: List[Path] = []
    for layout in LAYOUTS:
        cands.append(root / _DIR_TABLE[layout]['d_centroids'] / CENTROIDS_BASENAME)
    for sub in _LEGACY_CENTROID_DIRS:                     # 旧名兜底
        cands.append(root / sub / CENTROIDS_BASENAME)
    for p in cands:
        if p.is_file():
            return p
    if root.is_dir():                                     # 最后兜底: 限深递归
        hits = sorted(root.glob(f"*/*/{CENTROIDS_BASENAME}"))
        if hits:
            return hits[0]
    return None


def locate_taxonomy(root: Path) -> Optional[Path]:
    """定位 8 工具加权投票的整合分类表 (两布局)。"""
    root = Path(root)
    for layout in LAYOUTS:
        p = root / _DIR_TABLE[layout]['d_taxonomy'] / TAXONOMY_RELPATH
        if p.is_file():
            return p
    if root.is_dir():
        hits = sorted(root.glob(f"*/*/integrated/{os.path.basename(TAXONOMY_RELPATH)}"))
        if hits:
            return hits[0]
    return None


def write_boundary(root: Path, boundary: str, entries: Dict[str, str],
                   produced_by: str = "") -> Path:
    """把一次边界交接登记到 <项目根>/90_Handoff/handoff_manifest.tsv。

    90_Handoff 属于项目根 (跨管线共有): 若传入的 root 本身是某个管线独立根
    (01_PublicData/02_Preprocessing/03_Discovery/04_Analysis/05_Phylo),
    自动上提到其父目录再落 90_Handoff/。共享根布局 (legacy) 下 root 即项目根,
    行为不变。幂等: 同 (boundary, key) 重写其行, 其余行保留。
    """
    root = Path(root)
    _PIPELINE_ROOTS = {"01_PublicData", "02_Preprocessing", "03_Discovery",
                       "04_Analysis", "05_Phylo"}
    if root.name in _PIPELINE_ROOTS:
        root = root.parent
    handoff_dir = root / _DIR_TABLE["standard"]["x_handoff"]
    handoff_dir.mkdir(parents=True, exist_ok=True)
    manifest = handoff_dir / "handoff_manifest.tsv"

    rows: "Dict[tuple, Dict[str, str]]" = {}
    if manifest.is_file():
        with open(manifest, encoding="utf-8") as f:
            header = f.readline().rstrip("\n").split("\t")
            for line in f:
                parts = line.rstrip("\n").split("\t")
                rec = dict(zip(header, parts))
                rows[(rec.get("boundary", ""), rec.get("key", ""))] = rec

    import datetime
    now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    for key, path in entries.items():
        rows[(boundary, key)] = {"boundary": boundary, "key": key,
                                 "path": str(path), "produced_by": produced_by,
                                 "produced_at": now}

    header = ["boundary", "key", "path", "produced_by", "produced_at"]
    with open(manifest, "w", encoding="utf-8", newline="") as f:
        f.write("\t".join(header) + "\n")
        for (bnd, _key), rec in rows.items():
            f.write("\t".join(rec.get(h, "") for h in header) + "\n")
    return manifest


if __name__ == "__main__":
    for lay in LAYOUTS:
        print(f"── {lay} ──")
        for k, v in _DIR_TABLE[lay].items():
            print(f"  {k:14s} {v}")
