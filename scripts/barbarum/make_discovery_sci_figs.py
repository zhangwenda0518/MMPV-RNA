#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""make_discovery_sci_figs.py — 宁夏枸杞 (Lycium barbarum) 新病毒发现模块 SCI 主图

背景
----
`virome_discovery_pipeline` 的产物（`02_novel_virus/RNA-Lycium_barbarum_out/`）此前只有
一份交互式 HTML 报告（`report_pipeline.py` → `10_Reports/pipeline_report.html`），
**没有任何静态出版级图件**。已知病毒模块有 `S9_SCI_Composite/Fig_KnownVirome_Composite.*`
（`scripts/barbarum/make_sci_composite_figs.py`），发现模块缺一张对应主图 —— 本脚本补这张。

数据口径（全部来自 `10_Reports/` 的汇总表，脚本不重算任何统计量）
--------------------------------------------------------------
    分母 run 数        : 01_Assembly 下的样本目录数（权威）
    组装 contigs       : 10_Reports/assembly_summary.tsv 求和（或 stage_summary 的 Key_Metric）
    鉴定序列           : 10_Reports/ident_summary.tsv 的 TOTAL 行（多工具并集）
    合并 / 聚类        : 10_Reports/cluster_pipeline_reduction.tsv
    分类 / 宿主 / 质检 : taxonomy_composition.tsv / host_distribution.tsv / checkv_summary.tsv
    植物病毒           : all_plant_viruses_genus_summary.tsv / final_judgement_table.tsv
    类病毒             : Viroid.species_info.tsv

⚠ 口径纪律：**不要**把「全部 20,889 条 vOTU」当作植物病毒 —— 枸杞 RNA-seq 的病毒组
主体是原生生物/藻类/细菌/真菌相关的大 DNA 病毒（Nucleocytoviricota 10,280 条）。
植物相关子集只有 `all_plant_viruses_*` / `final_judgement_table.tsv` 那一支。图上必须标清。

配色
----
统一走 `virome_discovery_pipeline/utils/sci_palette.py`（与 phylo 侧同源）。
图内**不出现中文**（`report_pipeline.py` 的 `_en()` 词表只用于 HTML 报告，静态图全英文）。

用法
----
    python scripts/barbarum/make_discovery_sci_figs.py \
      --root  <RNA-Lycium_barbarum_out> \
      --outdir <...>/10_Reports/S10_SCI_Discovery \
      --virus-name "Lycium barbarum" \
      [--panels-only] [--composite-only] [--dpi 300]

产出
----
    <outdir>/Fig_GojiVirome_Discovery.{png,pdf}   组合主图 (A–I, v13 九面板)
    <outdir>/panels/panel_A..I_*.{png,pdf}        单面板（供排版/审稿人索取）
    <outdir>/_build/fig_manifest.json             面板清单 + 每个面板的数值来源
    <outdir>/_assets/sankey_flow.png              桑基预渲染图（供无 Chrome 的服务器复用）

v13 版式（9 面板 / 4 种图型，2026-10-12 重构）
---------------------------------------------
    A 漏斗(条形) B 鉴定工具(条形) C 簇大小(条形) **D 分类桑基(6 层)**
    E 新病毒分层 **F 宿主 vs HQ(成对条)** G 植物属(条形) H 终判(条形)
    **I 基因组型 + 类病毒(并排)**
    原 E(分类填充率)下线；原 G+H 合并为 F；原 K+L 合并为 I。

v14（2026-10-12，借 github.com/gabrielvpina/dataViz 参考图补 2 面板）
-------------------------------------------------------------------
    **J** nt identity × coverage 散点（+ 边缘箱线，着 = verdict，参考图 C）
    **K** 按科的 nt identity 箱线（参考图 F）

v15（2026-10-12，对齐参考论文 Viruses 2025, 17, 624 的 Fig2 版式）
------------------------------------------------------------------
    **L** 保守域完整性饼图（Complete / Fragment / absent，对齐论文 **Fig2-E**）
    至此 12 面板 A–L。
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import os
import re
import sys
import traceback
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Patch
import matplotlib.patheffects as pe
import numpy as np

# ── 统一配色（唯一真源） ───────────────────────────────────────────
# ⚠ **不要**写成 `from utils.sci_palette import ...`：服务器
# `~/mambaforge/lib/python3.10/site-packages/utils/` 里有一个**第三方同名包**
# （有 `__init__.py` 的常规包）。常规包**优先于命名空间包**，与 sys.path 顺序无关，
# 于是 `virome_discovery_pipeline/utils/`（无 `__init__.py`）会被静默盖掉 →
# `ModuleNotFoundError: No module named 'utils.sci_palette'`（2026-10-12 实测）。
# 这里按**显式文件路径**加载，既不依赖包名，也不受 site-packages 干扰。
import importlib.util as _ilu

_HERE = Path(__file__).resolve()
_REPO = _HERE.parents[2]                       # MMPV-RNA/


def _load_palette():
    for rel in ("virome_discovery_pipeline/utils/sci_palette.py",
                "virome_phylo_pipeline/utils/sci_palette.py"):
        cand = _REPO / rel
        if cand.is_file():
            spec = _ilu.spec_from_file_location("_sci_palette", cand)
            mod = _ilu.module_from_spec(spec)
            sys.modules["_sci_palette"] = mod
            spec.loader.exec_module(mod)
            return mod
    raise ImportError(f"sci_palette.py not found under {_REPO}")


_P = _load_palette()
apply_sci_style = _P.apply_sci_style
colors = _P.colors
OKABE_ITO = _P.OKABE_ITO
SEMANTIC = _P.SEMANTIC
SEQ_CMAP = _P.SEQ_CMAP
COL_W_1 = _P.COL_W_1
sci_legend = _P.sci_legend
despine = _P.despine

DEFAULT_DPI = 300


# ══════════════════════════════════════════════════════════════════
# I/O 辅助
# ══════════════════════════════════════════════════════════════════
def read_tsv(path: Path):
    """读 TSV → list[dict]（缺文件返回 None，便于上层做退化守卫）。"""
    if not path.is_file():
        return None
    with path.open(encoding="utf-8", newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    return rows


def _num(s, default=None):
    """宽松取数：'1,234' / '12.5%' / 'NA' / '""' → float。"""
    if s is None:
        return default
    t = str(s).strip().strip('"')
    if t in ("", "NA", "N/A", "-", "nan", "None"):
        return default
    t = t.replace(",", "").rstrip("%")
    try:
        return float(t)
    except ValueError:
        return default


def _int(s, default=None):
    v = _num(s)
    return default if v is None else int(round(v))


def _float(s, default=None):
    """宽松取数 → float（`_num` 的语义别名，便于表达「这里要小数」）。"""
    return _num(s, default)


def _find_reports(root: Path) -> Path:
    """在 root 下定位 10_Reports（兼容 09_Reports / reports）。"""
    for name in ("10_Reports", "09_Reports", "09_report", "reports"):
        p = root / name
        if p.is_dir():
            return p
    return root


def parse_stage_summary(rep: Path):
    """把 stage_summary.tsv 的 Key_Metric 解析成 {stage: {k: v}}。

    ⚠ Key_Metric 是**给人看的字符串**（"446 样本, 14,387,168 contigs, 5357.5 Mb"），
    不是稳定 schema。这里只提取本图要用的字段，取不到就留 None（面板自行退化）。
    """
    rows = read_tsv(rep / "stage_summary.tsv")
    out = {}
    if not rows:
        return out
    for r in rows:
        stage = (r.get("Stage") or "").strip()
        if not stage or stage.startswith("└"):      # 样本级明细行跳过
            continue
        km = r.get("Key_Metric") or ""
        d = {"status": (r.get("Status") or "").strip(), "raw": km}
        m = re.search(r"([\d,]+)\s*contigs", km)
        if m:
            d["contigs"] = _int(m.group(1))
        m = re.search(r"([\d,]+)\s*Mb", km)
        if m:
            d["mb"] = _num(m.group(1))
        m = re.search(r"([\d,]+)\s*(?:samples|样本)", km)
        if m:
            d["samples"] = _int(m.group(1))
        m = re.search(r"([\d,]+)\s*viral sequences", km)
        if m:
            d["viral_seqs"] = _int(m.group(1))
        m = re.search(r"combined=([\d,]+)", km)
        if m:
            d["combined"] = _int(m.group(1))
        # ⚠ Key_Metric 是**本地化字符串**：`report_pipeline.py`（有 `_en()` 中译英）
        #   重跑会输出英文（"20935 seqs"），旧版产物是中文（"20935 条"）。
        #   只认一种会让面板 A 静默丢环节（实测丢 2 条：分类 / 宿主预测）。
        m = re.search(r"([\d,]+)\s*(?:条|seqs)", km)
        if m:
            d["classified"] = _int(m.group(1))
        m = re.search(r"(\d+)\s*HQ\s*/\s*([\d,]+)\s*total", km)
        if m:
            d["hq"], d["checkv_total"] = _int(m.group(1)), _int(m.group(2))
        m = re.search(r"(\d+)\s*viruses\s*\(known=(\d+)\s*\+\s*rescued=(\d+)\)", km)
        if m:
            d["rescue_total"] = _int(m.group(1))
            d["rescue_known"] = _int(m.group(2))
            d["rescue_rescued"] = _int(m.group(3))
        # 06_HostPrediction 的 Key_Metric 形如
        #   "20892 seqs, Protist=7595 | Animal=6204 | ... | Plant=1063 | ..."
        # → 植物宿主子集数（**权威口径**；注意它 ≠ Plant.classified.fasta 的 1,061，
        #   后者是 2026-09-15 旧快照，未随宿主判定重建）。
        m = re.search(r"Plant=([\d,]+)", km)
        if m:
            d["plant_host"] = _int(m.group(1))
        m = re.search(r"HQ viruses:\s*([\d,]+)\s*seqs", km)
        if m:
            d["hq_plant_seqs"] = _int(m.group(1))
        out[stage] = d
    return out


# ══════════════════════════════════════════════════════════════════
# 数据装配
# ══════════════════════════════════════════════════════════════════
class Data:
    """所有面板的数据源。每个属性要么是可用结构，要么是 None（退化守卫）。"""

    def __init__(self, root: Path, virus_name: str):
        self.root = root
        self.virus_name = virus_name
        self.rep = _find_reports(root)
        self.manifest = {}          # panel -> {"source":..., "note":...}
        self._outdir = None         # 由 main() 注入，供 _find_asset 找预渲染资产
        self._sankey_plant_warn = None   # 植物过滤不可用时记一条（不静默）
        self._load()

    # ── 加载 ───────────────────────────────────────────────────
    def _load(self):
        rep = self.rep
        self.stage = parse_stage_summary(rep)
        self.assembly = read_tsv(rep / "assembly_summary.tsv")
        self.ident = read_tsv(rep / "ident_summary.tsv")
        self.cluster_red = read_tsv(rep / "cluster_pipeline_reduction.tsv")
        self.cluster_size = read_tsv(rep / "cluster_size_distribution.tsv")
        # 面板 A 的「cobra」环节：逐样本摘要（`Total_Queries` 汇总 = 35,348）
        self.cobra = read_tsv(rep / "cobra_summary.tsv")
        self.tax_comp = read_tsv(rep / "taxonomy_composition.tsv")
        self.tax_fill = read_tsv(rep / "taxonomy_fill_stats.tsv")
        self.host_dist = read_tsv(rep / "host_distribution.tsv")
        self.checkv = read_tsv(rep / "checkv_summary.tsv")
        self.hq_genome = read_tsv(rep / "hq_genome_type_summary.tsv")
        self.hq_family = read_tsv(rep / "hq_genome_type_family.tsv")
        self.plant_genus = read_tsv(rep / "all_plant_viruses_genus_summary.tsv")
        self.viroid_sp = read_tsv(rep / "Viroid.species_info.tsv")
        self.rescue = read_tsv(rep / "rescue_report.tsv")
        # Fig2-E「结构域完整性」的唯一正确来源：`cdd_evidence_report.tsv` 的 `evidence`
        #   列（逐条 CDD 保守域证据：PASS_CORE/PASS_CORE_FAM/PASS_VIRAL/AMBIGUOUS/REVIEW）。
        self.cdd = read_tsv(rep / "cdd_evidence_report.tsv")
        # Fig2-B「比对靶标」：HMM 命中的最佳病毒域描述（`ct3_best_viral_desc`）。
        self.ct3 = read_tsv(rep / "hmm_ct3_hits.tsv")
        # ⚠ 参考图（dataViz `data_report.R`）的散点图 C 与箱线图 F 需要**逐 contig 的
        #   identity / coverage / 科归属**，只有 `rescue_evidence_scored.tsv` 同时有
        #   `nt_pident`（成对相似度 %）+ `nt_qcovs`（覆盖度 %）+ `tax_family` + `verdict`。
        #   不要改用 `All_plant.viruses_info.tsv`：它只有 `confidence`（1/0.88/0.75… 是
        #   证据分而非 % identity），画成"identity 散点"会张冠李戴。
        self.evidence = read_tsv(self._locate("rescue_evidence_scored.tsv"))
        self.judge = read_tsv(self._locate("final_judgement_table.tsv"))
        self.keep = read_tsv(rep / "keep_summary.tsv")
        # ⚠ **面板 A 的「过滤」环节**：`02b_Filter` 的高置信过滤产物。
        #   管线数据流（`virome_pipeline.py` 顶部箭头图）：
        #     `02a_Identification → 02b_Filter → 03a_COBRA → 03b_MergeSamples → 04_CLUSTER`
        #   ⚠ `filter_summary.tsv` 在本项目里**是空的（396 行全 0）**：`ident_stats.py`
        #     按 `uniprot_filter_output*` 目录名找产物，与实际目录名（`uniprot/`）
        #     不匹配 → 静默全 0。**不能**用它。
        #   故 `filter_per_sample.tsv` 由 `02b_Filter/<样本>/uniprot/
        #   all_candidates_uniprot.fasta` 的 `>` 计数汇总而来（权威、可复核）。
        self.filter_ps = read_tsv(rep / "filter_per_sample.tsv")
        # ⚠ **桑基图的唯一正确数据源**：逐 contig 的完整分类谱系（Realm→…→Genus，
        #   20,935 行，各阶元填充率与 taxonomy_fill_stats.tsv 逐位吻合）。
        #   **不能**用 taxonomy_composition.tsv —— 它只有 (`Rank`,`Taxon`,`Count`) 三列、
        #   各阶元是**独立汇总**、没有 lineage → 桑基的层级边无从谈起（会凭空编造归属）。
        self.integrated = read_tsv(rep / "final_integrated_classification.tsv")
        # ⚠ **面板 D 换植物桑基后需要的宿主映射**：逐 contig 的权威宿主判定
        #   （`Final_Host`）。锚点与 `plant_denom()` 同源 → 过滤出的植物子集
        #   必然 = 1,063，不会与面板 A/G/I 的口径打架。
        #   查找：报告目录 → `06_HostPrediction/`（**不能**用会漂移的
        #   `plant_final_taxonomy.tsv` 中间产物，它是旧快照、且列结构会变）。
        self.host_map = (read_tsv(rep / "ensemble_host_summary.tsv")
                         or read_tsv(self.root / "06_HostPrediction" / "ensemble_host_summary.tsv")
                         or read_tsv(self._locate("ensemble_host_summary.tsv")))
        # 权威 run 数：01_Assembly 下的样本目录数
        asm_root = self.root / "01_Assembly"
        self.n_runs = None
        if asm_root.is_dir():
            self.n_runs = len([d for d in asm_root.iterdir()
                               if d.is_dir() and d.name.lower() not in ("logs", "log")])
        if self.n_runs is None and self.stage.get("01_Assembly", {}).get("samples"):
            self.n_runs = self.stage["01_Assembly"]["samples"]

    def _locate(self, name: str) -> Path:
        for cand in (self.rep / name, self.root / name,
                     self.rep.parent / name, self.root.parent / name):
            if cand.is_file():
                return cand
        return self.rep / name

    # ── 派生量 ─────────────────────────────────────────────────
    def funnel(self):
        """发现管线主干计数（v18：按用户指令补齐为 13 环节全链）。

        用户给的链条（2026-10-13）：
          组装 → 病毒 → 过滤 → cobra → merge → cluster → 分类 → 宿主 → 植物
          → HQ → rescue → high → final HQ

        ⚠ **口径分两块，不可读成单一漏斗**（用户要求「全部纳入 + 逐项标注来源」）：
          【A】全宿主块（前 8 项）= 枸杞 RNA-seq 病毒组的主体是原生生物/藻/菌/
               细菌相关大 DNA 病毒，故"宿主"= 全 8 类宿主之和。
          【B】植物块（后 5 项）= 从【A】的 `Host-assigned` 里切出的植物子集，
               它是【A】的**切片**而非下一级 —— 图上用同色 + 括号注明。

        ⚠ 各环节来源（逐项可复核）：
          · `Assembled contigs`  = `01_Assembly` contigs（14,387,168）
          · `Viral sequences`    = `02_Identification` viral_seqs（556,773）
          · `Filtered candidates`= `02b_Filter/<样本>/uniprot/all_candidates_uniprot.fasta`
             的序列数汇总（**94,826**）。**不能**读 `filter_summary.tsv` —— 本项目里
             它是 **396 行全 0 的空表**（`ident_stats.py` 按 `uniprot_filter_output*`
             目录名找产物，与实际 `uniprot/` 不匹配 → 静默全 0）。
          · `COBRA queries`      = `03a_COBRA` 逐样本 `cobra_summary.tsv` 的
             `Total_Queries` 汇总 = **35,348**（延伸后的 query 数，可 > 输入）。
          · `Merged / re-assembled` = `03b_Merge` combined（80,637）
          · `vOTU clusters`      = `cluster_pipeline_reduction.tsv`（37,803）
          · `Taxonomically classified` = `05_Taxonomy` classified（20,935）
          · `Host-assigned`      = `06_HostPrediction` classified（20,892）
          【B】植物块 ──
          · `Plant-host (subset)`= 同上 `Plant=` 行（**1,063**）
          · `Plant HQ vOTUs`     = `checkv_summary.tsv` **Plant 行**的 `HQ` 列
             （**63** = Complete 33 + High-quality 30）。用户 2026-10-13 明确选定
             此口径（不取 `07_CheckV` 全宿主 200）。
             ⚠ 该表 Plant 分母是 **1,090**（≠ 本图其它处的 1,063）——属已知 R6
               分歧（两套植物口径），图上按用户要求**逐项标注来源**。
          · `Rescued (VSI)`      = `08_Rescue` 的 rescued（**276**）
          · `CheckV high-confidence` = `keep_summary.tsv` 中
             `checkv_confidence == "high"` 条数（**45**）
          · `Final HQ plant viruses` = `10_Virome_Analysis` 的 `HQ viruses`（**291**）
        """
        s = self.stage
        items = []
        # ──【A】全宿主块 ──
        a = s.get("01_Assembly", {})
        if a.get("contigs"):
            items.append(("Assembled contigs", a["contigs"]))
        i = s.get("02_Identification", {})
        if i.get("viral_seqs"):
            items.append(("Viral sequences (10 tools)", i["viral_seqs"]))
        nfilt = self.filter_total()
        if nfilt:
            items.append(("Filtered candidates", nfilt))
        ncobra = self.cobra_queries()
        if ncobra:
            items.append(("COBRA queries (extended)", ncobra))
        m = s.get("03b_Merge", {})
        if m.get("combined"):
            items.append(("Merged / re-assembled", m["combined"]))
        if self.cluster_red:
            d = {r["Stage"]: _int(r["Count"]) for r in self.cluster_red}
            if d.get("vOTU clusters"):
                items.append(("vOTU clusters", d["vOTU clusters"]))
        c = s.get("05_Taxonomy", {})
        if c.get("classified"):
            items.append(("Taxonomically classified", c["classified"]))
        h = s.get("06_HostPrediction", {})
        if h.get("classified"):
            items.append(("Host-assigned", h["classified"]))
        # ──【B】植物子集块（与【A】不同轨，图上标清）──
        if h.get("plant_host"):
            items.append(("Plant-host (subset)", h["plant_host"]))
        nplant_hq = self.plant_hq()
        if nplant_hq:
            items.append(("Plant HQ vOTUs", nplant_hq))
        r = s.get("08_Rescue", {})
        if r.get("rescue_rescued"):
            items.append(("Rescued (VSI, plant subset)", r["rescue_rescued"]))
        nhigh = self.keep_high()
        if nhigh:
            items.append(("CheckV high-confidence", nhigh))
        nhq = self.final_hq()
        if nhq:
            items.append(("Final HQ plant viruses", nhq))
        return items

    def filter_total(self):
        """02b_Filter 的高置信过滤总量（逐样本 uniprot fasta 计数汇总）。"""
        if not self.filter_ps:
            return None
        n = 0
        for r in self.filter_ps:
            n += _int(r.get("Uniprot_Filtered") or 0)
        return n or None

    def cobra_queries(self):
        """03a_COBRA 的延伸 query 总量（`cobra_summary.tsv` 的 `Total_Queries` 汇总）。"""
        if not self.cobra:
            return None
        n = 0
        for r in self.cobra:
            n += _int(r.get("Total_Queries") or 0)
        return n or None

    def keep_high(self):
        """`keep_summary.tsv` 中 `checkv_confidence == "high"` 的条数。"""
        if not self.keep:
            return None
        return sum(1 for r in self.keep
                   if (r.get("checkv_confidence") or "").strip().lower() == "high") or None

    def final_hq(self):
        """`10_Virome_Analysis` 的 `HQ viruses`（最终结果，植物 HQ 病毒）。

        ⚠ 该值来自 `stage_summary.tsv` 的 `10_Virome_Analysis` 行，形如
          `HQ viruses: 291 seqs | All plant: 1061 seqs | …`。
        """
        raw = self.stage.get("10_Virome_Analysis", {}).get("raw") or ""
        m = re.search(r"HQ\s+viruses:\s*([\d,]+)", raw)
        return _int(m.group(1)) if m else None

    def plant_hq(self):
        """植物子集的 CheckV HQ（Complete + High-quality）——**用户 2026-10-13 选定口径**。

        ⚠ 来源 = `checkv_summary.tsv` 的 **Plant 行** `HQ` 列（= 63）。
          **不取** `07_CheckV` 的全宿主 HQ（= 200）—— 用户明确要植物子集的口径。
        ⚠ 该表 Plant 行的分母是 **1,090**，与面板 A/D/I 其它处的 **1,063**
          不同源（属已知 R6「两套植物口径」分歧）。图上必须逐项标注来源，
          否则读者会拿 63/1,063 与 63/1,090 对不上。
        """
        if not self.checkv:
            return None
        for r in self.checkv:
            if (r.get("Host") or "").strip().lower() == "plant":
                return _int(r.get("HQ")) or None
        return None

    def ident_tools(self):
        """ident_summary.tsv 的 TOTAL 行 → [(tool, n)]，降序。

        ⚠ 第 2 列 `All_Candidate` 是**多工具并集**（= 556,773），不是工具；
        且 virsorter2 / virhunter / viralm 三列为 0（未产出）。三者都要排除，
        否则并集会与工具并列，读成"某个工具贡献了 55 万条"。
        """
        if not self.ident:
            return None
        cols = list(self.ident[0].keys())
        tot = None
        for r in self.ident:
            if (r.get(cols[0]) or "").strip() == "TOTAL":
                tot = r
                break
        if tot is None:
            return None
        skip = {cols[0], "All_Candidate", "summary_plots"}
        vals = []
        for c in cols[1:]:
            if c in skip:
                continue
            n = _int(tot.get(c))
            if n and n > 0:
                vals.append((c, n))
        vals.sort(key=lambda x: -x[1])
        return vals or None

    def ident_union(self):
        """多工具并集（All_Candidate）—— 画成参考线，不当作工具。"""
        if not self.ident:
            return None
        cols = list(self.ident[0].keys())
        for r in self.ident:
            if (r.get(cols[0]) or "").strip() == "TOTAL":
                return _int(r.get("All_Candidate"))
        return None

    def novelty(self):
        """05_Taxonomy 的 Novel Rank → [(tier, n)]。

        ⚠ 该字段在 `stage_summary.tsv` 里是 **`└ Novel Rank` 子行**（前缀 `└`），
        不是 05_Taxonomy 主行 —— 主行只有 "20935 条, ★17282 已知, ★★1125 新种"。
        所以必须扫全表（含子行），不能只看非 `└` 行。
        """
        rows = read_tsv(self.rep / "stage_summary.tsv") or []
        for r in rows:
            km = r.get("Key_Metric") or ""
            m = re.search(r"Known=(\d+)\s*NewSp=(\d+)\s*NewGe=(\d+)\s*NewFa=(\d+)", km)
            if m:
                return [("Known", _int(m.group(1))), ("Novel species", _int(m.group(2))),
                        ("Novel genus", _int(m.group(3))), ("Novel family", _int(m.group(4)))]
        return None

    def checkv_matrix(self):
        """checkv_summary.tsv → (hosts, qualities, matrix)。去掉 TOTAL 行。"""
        if not self.checkv:
            return None
        quals = ["Complete", "High-quality", "Medium-quality",
                 "Low-quality", "Not-determined"]
        hosts, mat = [], []
        for r in self.checkv:
            h = (r.get("Host") or "").strip()
            if not h or h.upper() == "TOTAL":
                continue
            hosts.append(h)
            mat.append([_int(r.get(q), 0) or 0 for q in quals])
        if not hosts:
            return None
        order = np.argsort([-sum(x) for x in mat])
        return ([hosts[i] for i in order], quals,
                np.array([mat[i] for i in order], dtype=float))

    def plant_genera(self, top=18):
        if not self.plant_genus:
            return None
        rows = []
        for r in self.plant_genus:
            g = (r.get("Genus") or "").strip()
            n = _int(r.get("n_contigs"), 0) or 0
            if not g or g in ("NA", "-") or n <= 0:
                continue
            rows.append((g, n, (r.get("Genome_Type") or "none").strip()))
        rows.sort(key=lambda x: -x[1])
        return rows[:top] or None

    def plant_denom(self):
        """全图**唯一的植物宿主分母**（权威口径）。

        ⚠ 项目里有 5 套植物口径（详见 MEMORY_full「植物口径 5 套」）：
            · **1,063** = `host_distribution.tsv` 的 Plant 行 = `ensemble_host_summary.tsv`
              的 `Final_Host=="Plant"` 计数 = `stage_summary.tsv` 的 `Plant=`（**三处一致**）
              → **权威**，本函数返回它；
            · 1,090 = `checkv_summary.tsv` 的 Plant 行 `Total`（更宽/更旧，混入 27 条非植物）；
            · 1,061 = `Plant.classified.fasta`（2026-09-15 旧快照）；
            · 291 = HQ plant viruses；310 = plant-virus vOTUs。
        混用会让同一张组合图里出现两个植物分母（实测 G 用 1,063、H/I 用 1,090）→ 本函数
        统一取 1,063，H 面板的 checkv 分母另由 `checkv_matrix()` 内部按宿主组自洽处理，
        并在面板内标注来源，避免读者拿两个数对不上。
        """
        if self.host_dist:
            for r in self.host_dist:
                if (r.get("Host") or "").strip().lower() == "plant":
                    v = _int(r.get("Count"))
                    if v:
                        return v
        s = self.stage.get("06_HostPrediction", {})
        return s.get("plant_host")

    def plant_contig_ids(self):
        """`ensemble_host_summary.tsv` 里 `Final_Host == "Plant"` 的 contig_id 集合。

        与 `plant_denom()` **同源**（1,063）→ 用它过滤出的子集必然与面板 A/G/I
        标注的植物分母一致，不会出现"图上两个植物分母"。

        ⚠ 取不到 `host_map` 时返回 **None**（不是空集！）—— 空集会被下游当成
          "植物一条都没有"而画出**空桑基**，属静默错误；None 让下游**退回全量**并
          在 manifest 里记一条警告。
        """
        if not self.host_map:
            return None
        ids = {(r.get("contig_id") or "").strip() for r in self.host_map
               if (r.get("Final_Host") or "").strip() == "Plant"}
        ids.discard("")
        return ids or None

    def sankey_flows(self, topn=None, plant_only=False):
        """Realm → Kingdom → Phylum → Class → Family → Genus 六列的流量（供桑基图）。

        返回 `(nodes, links)`：nodes = [(col, name)]，links = [(src_idx, dst_idx, value)]。
        每级只保留 top-N，其余归 "Other"；空值统一记 "Unclassified"（**不丢行**，
        否则各列总量对不上，桑基会画成"越流越少"的假象）。

        plant_only=True → 只统计 `Final_Host == "Plant"` 的 contig（1,063 条），
        对齐管线报告里 `Plant Virus Taxonomy Sankey`（`classification_sankey_plant.html`）
        的口径。⚠ 过滤键是 `contig_id`；`final_integrated_classification.tsv` 与
        `ensemble_host_summary.tsv` 的 contig_id 空间一致（已验证 plant ⊂ integrated）。
        """
        if not self.integrated:
            return None
        rows = self.integrated
        if plant_only:
            pids = self.plant_contig_ids()
            if pids is None:
                # 不静默：退回全量，把原因记进 manifest（由 panel_D_sankey 落盘）
                self._sankey_plant_warn = (
                    "ensemble_host_summary.tsv not found -> plant-only filter skipped, "
                    "falling back to ALL hosts")
            else:
                rows = [r for r in self.integrated
                        if (r.get("contig_id") or "").strip() in pids]
        cols = ("Realm", "Kingdom", "Phylum", "Class", "Family", "Genus")
        # 每层保留的节点数：上两层少（本来就没几个），下三层多（种类多、要收敛）
        if topn is None:
            topn = {"Realm": 8, "Kingdom": 8, "Phylum": 10,
                    "Class": 10, "Family": 10, "Genus": 10}
        NA = "Unclassified"

        def _norm(x):
            t = (x or "").strip()
            return t if t and t not in ("NA", "-", "None", "nan") else NA

        # 每列统计 + 取 top-N（Unclassified 永远保留，否则最大的一类会被归到 Other）
        tops = {}
        for c in cols:
            n = topn.get(c, 10)
            cnt = {}
            for r in rows:
                k = _norm(r.get(c))
                cnt[k] = cnt.get(k, 0) + 1
            keep = [k for k, _ in sorted(cnt.items(), key=lambda x: -x[1])[:n]]
            tops[c] = keep

        def _collapse(c, k):
            return k if (k in tops[c] or k == NA) else "Other"

        # 逐级聚合（只连相邻两列）
        nodes, idx = [], {}

        def nid(col, name):
            key = (col, name)
            if key not in idx:
                idx[key] = len(nodes)
                nodes.append((col, name))
            return idx[key]

        links = {}
        for a, b in zip(cols, cols[1:]):
            for r in rows:
                ka, kb = _collapse(a, _norm(r.get(a))), _collapse(b, _norm(r.get(b)))
                s, d = nid(a, ka), nid(b, kb)
                links[(s, d)] = links.get((s, d), 0) + 1
        if not links:
            return None
        return nodes, [(s, d, v) for (s, d), v in links.items()]

    def host_vs_hq(self):
        """宿主组：预测总数 vs 其中 HQ（面板 F 用，把原 G+H 合成一张）。

        ⚠ 两个数**不同源**且**故意如此**：
          · 总数 = `host_distribution.tsv`（权威宿主判定，Plant=1,063）；
          · HQ   = `checkv_summary.tsv` 的 Complete+High-quality（Plant=63）。
        两者不要求相等（HQ 是 CheckV 对基因组的质检，不是宿主判定），但**类目必须对齐**
        → 按 `host_distribution.tsv` 的类目为基准，CheckV 里缺的补 0。
        """
        if not self.host_dist:
            return None
        hq = {}
        if self.checkv:
            for r in self.checkv:
                h = (r.get("Host") or "").strip()
                if not h or h.upper() == "TOTAL":
                    continue
                hq[h] = (_int(r.get("Complete"), 0) or 0) + (_int(r.get("High-quality"), 0) or 0)
        rows = []
        for r in self.host_dist:
            h = (r.get("Host") or "").strip()
            n = _int(r.get("Count"), 0) or 0
            if not h or n <= 0:
                continue
            rows.append((h, n, hq.get(h, 0)))
        rows.sort(key=lambda x: -x[1])
        return rows or None

    def genome_viroid(self):
        """面板 I：基因组型（HQ） + 类病毒种，合成一对横向条。

        ⚠ **术语撞车**：`hq_genome_type_summary.tsv` 的 "HQ" 指**高质量基因组型构成**
        （合计 895 = DNA 44 + RNA 237 + Viroid 604 + Unassigned 10），**不是** A 面板那个
        "HQ vOTUs = 200"（CheckV Complete+High-quality）。面板标题必须写明
        "high-quality genomes"，否则读者会把 895 和 200 当成同一个数。
        """
        out = {}
        if self.hq_genome:
            rows = [(r["Category"], _int(r.get("n_contigs"), 0) or 0)
                    for r in self.hq_genome
                    if (r.get("Category") or "").strip().upper() != "TOTAL"]
            rows = [r for r in rows if r[1] > 0]
            if rows:
                out["genome_type"] = rows
        if self.viroid_sp:
            vr = []
            for r in self.viroid_sp:
                sp = (r.get("species") or "").strip()
                n_c = _int(r.get("n_contigs"), 0) or 0
                n_s = _int(r.get("n_samples"), 0) or 0
                if not sp or (n_c == 0 and n_s == 0):
                    continue
                vr.append((sp, n_c, n_s))
            vr.sort(key=lambda x: -x[1])
            if vr:
                out["viroid"] = vr
        return out or None

    def verdicts(self):
        if not self.judge:
            return None
        from collections import Counter
        cnt = Counter()
        for r in self.judge:
            cat = (list(r.values())[-1] or "").strip()
            if cat:
                cnt[cat] += 1
        if not cnt:
            return None
        return cnt.most_common()

    # ── 参考图（dataViz）借用的两个面板的数据 ────────────────────
    def ident_vs_cov(self, min_pts=8):
        """面板：nt identity (%) × coverage (%) 散点（参考图 C）。

        返回 `(pts, xlab, ylab)`，`pts = [(ident, cov, verdict), ...]`。
        ⚠ 只取**两个坐标都有效**（>0 且非空）的行：`nt_pident=0` 表示无成对命中，
          画进散点会在 x=0 堆一条假竖线，把"低相似度"当成真实分布。
        """
        if not self.evidence:
            return None
        pts = []
        for r in self.evidence:
            x = _float(r.get("nt_pident"))
            y = _float(r.get("nt_qcovs"))
            if x is None or y is None or x <= 0 or y <= 0:
                continue
            v = (r.get("verdict") or "").strip() or "undetermined"
            pts.append((x, y, v))
        if len(pts) < min_pts:
            return None
        return pts, "Nucleotide identity (%)", "Coverage (%)"

    def family_identity(self, min_n=3, top=14):
        """面板：按科的 nucleotide identity 箱线（参考图 F）。

        返回 `[(family, [ident, ...]), ...]`，按中位数降序，只保留 `n >= min_n` 的科
        （`n < 3` 的箱线只有一条线、无信息，且会把 y 轴标签挤满 → 退化守卫剔除）。
        """
        if not self.evidence:
            return None
        from collections import defaultdict
        g = defaultdict(list)
        for r in self.evidence:
            f = (r.get("tax_family") or "").strip()
            x = _float(r.get("nt_pident"))
            if not f or f in ("NA", "-") or x is None or x <= 0:
                continue
            g[f].append(x)
        rows = [(f, v) for f, v in g.items() if len(v) >= min_n]
        if len(rows) < 3:
            return None
        rows.sort(key=lambda r: -float(np.median(r[1])))
        return rows[:top] or None

    # ── 参考论文（Viruses 2025, 17, 624 / PMC12115555）Fig2-E 的数据 ──
    def domain_completeness(self):
        """Fig2-E「complete / fragmented / absent」结构域完整性三分。

        参考论文的口径（`dataViz/data_report.R`）：
            `domains <- ifelse(filo$Contig.Domains == filo$BestHit.Domains,
                               "Complete", "Fragment")` —— 即**逐条**比对
            「我们的 contig 命中域集合」与「最佳命中（参考）的域集合」，完全相同 = Complete，
            否则 = Fragment。

        我们这一侧最接近、且**可核验**的口径是 `cdd_evidence_report.tsv` 的
        `evidence` 列（逐条 CDD/保守域证据判定）：
            · PASS_CORE / PASS_CORE_FAM  = 命中该科的**核心**域全集 → **Complete**
            · PASS_VIRAL                 = 命中泛病毒域但缺科核心 → **Fragment**
            · AMBIGUOUS / REVIEW         = 未命中任何可用域证据     → **No conserved domain**
        ⚠ 不要改用 `checkv_summary.tsv`：CheckV 的 Complete/High-quality 是**基因组完整性**
          （基因组长度的完整度），不是「保守域齐不齐」，画成同一个饼会张冠李戴。
        """
        if not getattr(self, "cdd", None):
            return None
        groups = {"Complete": 0, "Fragment": 0, "No conserved domain": 0}
        seen = 0
        for r in self.cdd:
            ev = (r.get("evidence") or "").strip().upper()
            if not ev:
                continue
            seen += 1
            if ev in ("PASS_CORE", "PASS_CORE_FAM"):
                groups["Complete"] += 1
            elif ev == "PASS_VIRAL":
                groups["Fragment"] += 1
            else:
                groups["No conserved domain"] += 1
        if seen < 8:
            return None
        # 固定顺序（Complete → Fragment → absent），颜色走 SEMANTIC["quality"]
        order = ["Complete", "Fragment", "No conserved domain"]
        return [(k, groups[k]) for k in order if groups[k] > 0], seen

    # ── 参考论文 Fig2-B「alignment targets of viral contigs」──────
    _TARGET_RULES = (
        # 顺序敏感：先匹配更专一的（RdRp 先于泛 polyprotein）
        ("RdRp (RNA-dependent RNA polymerase)", ("rna-dependent-rna-polymerase", "rdrp")),
        ("RT / RNase H", ("reverse-transcriptase", "rnase-h")),
        ("Capsid / coat protein", ("capsid", "-coat-protein", "coat-protein")),
        ("Nucleoprotein", ("nucleoprotein",)),
        ("Helicase", ("helicase",)),
        ("Protease / helper component", ("proteinase", "helper-component", "protease")),
        ("Movement protein", ("movement-protein", "cell-to-cell")),
        ("Replicase", ("replicase",)),
        ("Integrase", ("integrase",)),
    )

    def protein_targets(self, min_n=2):
        """Fig2-B：病毒 contig 的**比对靶标**（RdRp / Capsid / 多聚蛋白 / …）分布。

        参考论文 Fig2-B: `Alignment targets of viral contigs`（蛋白功能类目直方图）。
        ⚠ 论文的类目是「polyprotein / Capsid / RdRp / …」；我们这一侧**没有**现成的
          蛋白功能类目表，只能从 `hmm_ct3_hits.tsv` 的 `ct3_best_viral_desc`
          （HMM 命中的**最佳病毒域**描述串，如 `PFAM:PF00978.21-RNA-dependent-
          RNA-polymerase`）里**按关键词归类**。这是**派生量、非原始统计量**，
          图上与 manifest 都写明归类规则，避免被当成权威注释。
        """
        if not getattr(self, "ct3", None):
            return None
        from collections import Counter
        cnt = Counter()
        unmapped = 0
        seen = 0
        for r in self.ct3:
            d = (r.get("ct3_best_viral_desc") or "").strip().lower()
            if not d:
                continue
            seen += 1
            for label, keys in self._TARGET_RULES:
                if any(k in d for k in keys):
                    cnt[label] += 1
                    break
            else:
                unmapped += 1
        if seen < 8:
            return None
        rows = [(k, v) for k, v in cnt.most_common() if v >= min_n]
        if len(rows) < 2:
            return None
        return rows, seen, unmapped


# ══════════════════════════════════════════════════════════════════
# 绘图辅助
# ══════════════════════════════════════════════════════════════════
def _empty(ax, msg="(no data)"):
    ax.axis("off")
    ax.text(0.5, 0.5, msg, ha="center", va="center", color="#999999", fontsize=7,
            transform=ax.transAxes)


def _box(ax, data, horizontal=False, **kw):
    """版本无关的箱线调用。

    ⚠ matplotlib 3.11 起 `vert=` 弃用（3.13 移除）→ 改用 `orientation=`；但服务器
      mpl 版本更旧、**不认** `orientation` → 直接写任一个都会在另一端报错/告警。
      故按版本探测后再传（本图要在本地 3.11 与服务器两处都能跑）。
    """
    kw = dict(kw)
    try:
        major, minor = (int(x) for x in matplotlib.__version__.split(".")[:2])
    except Exception:
        major, minor = 3, 10
    if (major, minor) >= (3, 11):
        kw["orientation"] = "horizontal" if horizontal else "vertical"
    else:
        kw["vert"] = not horizontal
    return ax.boxplot(data, **kw)


def _log_ticks(ax, lo, hi):
    """给 log 轴按数量级放刻度，并用千分位纯数字标签（不要 10^n）。"""
    if lo <= 0 or hi <= lo:
        return
    dec_lo, dec_hi = math.floor(math.log10(lo)), math.ceil(math.log10(hi))
    ticks = [10 ** k for k in range(dec_lo, dec_hi + 1)]
    ax.set_xticks(ticks)
    ax.set_xticklabels([f"{t:,}" if t < 1e6 else f"{t/1e6:g}M" for t in ticks],
                       fontsize=5.4)
    ax.set_xlim(lo * 0.45, hi * 2.2)


def _hbar_labels(ax, labels, vals, bar_colors, fmt="{:,}", log=False,
                 title=None, xlabel=None, label_fs=6.2, val_fs=5.6,
                 xlabel_fs=6.4):
    y = np.arange(len(labels))[::-1]
    ax.barh(y, vals, color=bar_colors, height=0.72,
            edgecolor="#333333", linewidth=0.4)
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=label_fs)
    if log:
        ax.set_xscale("log")
        _log_ticks(ax, min(vals), max(vals))
        for yy, v in zip(y, vals):
            ax.text(v * 1.25, yy, fmt.format(v), va="center", ha="left",
                    fontsize=val_fs, color="#333333")
    else:
        ax.margins(x=0.22)
        xmax = max(vals) if vals else 1
        for yy, v in zip(y, vals):
            ax.text(v + xmax * 0.02, yy, fmt.format(v), va="center", ha="left",
                    fontsize=val_fs, color="#333333")
    if title:
        ax.set_title(title, fontsize=7.2, fontweight="bold", loc="left")
    if xlabel:
        ax.set_xlabel(xlabel, fontsize=xlabel_fs)
    despine(ax)
    ax.tick_params(length=2)


def panel_letter(ax, letter, dx=-0.075, dy=1.045):
    ax.text(dx, dy, letter, transform=ax.transAxes, fontsize=9.5,
            fontweight="bold", va="bottom", ha="left", color="#111111")


# ══════════════════════════════════════════════════════════════════
# 各面板
# ══════════════════════════════════════════════════════════════════
def panel_A(ax, D: Data):
    """A 发现管线主干 13 环节（log 计数）：组装 → … → 宿主 │ 植物 → … → final HQ。

    v18（2026-10-13，用户指令）：补 `过滤 / cobra / 植物 HQ / CheckV high / 最终 HQ`
    五项，使链条完整覆盖：
      组装 → 病毒 → 过滤 → cobra → merge → cluster → 分类 → 宿主 → 植物 →
      HQ → rescue → high → final HQ
    ⚠ **两组口径**：前 8 项全宿主（all hosts），后 5 项植物子集（plant subset）。
      植物块是「宿主」那一档的**切片**、不是它的下一级 —— 图上用同色 + 右侧
      scope 括号 + 括注说清，避免被读成"20,892 → 1,063 掉了 95 %"的假漏斗。
    """
    items = D.funnel()
    if not items:
        return _empty(ax)
    labels = [k for k, _ in items]
    vals = [v for _, v in items]
    sc = SEMANTIC["stage"]
    keymap = {"Assembled contigs": "assembly",
              "Viral sequences (10 tools)": "identification",
              "Filtered candidates": "filter",
              "COBRA queries (extended)": "cobra",
              "Merged / re-assembled": "cobra",
              "vOTU clusters": "cluster",
              "Taxonomically classified": "taxonomy",
              "Host-assigned": "host",
              # ⚠ `Plant-host (subset)` 故意复用 `host` 色 = 与 `Host-assigned` 同色，
              #   让读者一眼看出"这是它的切片"，而不是又一条独立口径。
              "Plant-host (subset)": "host",
              "Plant HQ vOTUs": "checkv",
              "Rescued (VSI, plant subset)": "rescue",
              "CheckV high-confidence": "checkv",
              "Final HQ plant viruses": "final"}
    bc = [sc.get(keymap.get(k, ""), "#999999") for k in labels]
    _hbar_labels(ax, labels, vals, bc, log=True, xlabel="Count (log scale)",
                 title="Discovery pipeline reduction", label_fs=6.0)
    # ⚠ `_hbar_labels` 把 xlim 收到 `(min*0.45, max*2.2)`（值标签画在 `v*1.25`，
    #   故右侧余量只够放标签）。本面板还要在右侧放两块口径的 scope 括号 + 括注。
    #   做法：**括号与注记一律用 `ax.transAxes`（轴分数坐标）**，与数据范围无关
    #   —— 比"按 max(vals) 的倍数估数据坐标"稳健得多（后者在 log 轴上要反复
    #   反解 xlim，实测会把注记顶出画布 → 体检报 `TEXT-OFF-CANVAS`）。
    #   唯一的数据坐标动作：把右界从 2.2× 撑到 4.5×，给条形右侧留出横向净空。
    if vals:
        ax.set_xlim(min(vals) * 0.45, max(vals) * 4.5)
    # ⚠⚠ 面板变高后（v21 上排 3.7 in）默认 y 余量会撑出上下空白 → 显式收紧到条外 0.4 单位
    #   （13 条横条在 3.7 in 里每条约 0.26 in，留 0.4 单位 ≈ 0.11 in 边距足够）。
    ax.set_ylim(-0.6, len(labels) - 0.4)
    # ── 两块口径的 scope 括号 + 括注（轴分数坐标，全部落在轴内）──
    y = np.arange(len(labels))[::-1]
    PLANT_ITEMS = ("Plant-host (subset)", "Plant HQ vOTUs",
                   "Rescued (VSI, plant subset)", "CheckV high-confidence",
                   "Final HQ plant viruses")
    ia = [i for i, k in enumerate(labels) if k not in PLANT_ITEMS]   # all-host
    ip = [i for i, k in enumerate(labels) if k in PLANT_ITEMS]       # plant
    ymin, ymax = y.min() - 0.5, y.max() + 0.5
    def _af(yy):                       # data-y → axes fraction
        return (yy - ymin) / (ymax - ymin)
    x_brace = 0.655                    # 括号线的轴分数位置
    # ⚠⚠ **括号 + note 的水平位置在半栏版式下必须重算**（v21 实测）：
    #   v21 上排是 A|J 并排 → A 轴宽仅 **4.12 in**（v20 通栏 7.61 in）。
    #   实测 A 轴内文本占用（`get_window_extent` 换算轴分数）：
    #     值标签最右 = `14,387,168` → **0.914–1.026**（本身已伸出轴外）
    #     `80,637` 0.568–0.636、`35,348` 0.513–0.581
    #   → 括号放 0.735 会把组名压在值标签上；note 放 0.755 会伸出轴右缘撞 J 的刻度。
    #   正解 = ① 括号移到 **0.655**（值标签右侧的净空带）；
    #           ② note **改用极短文案**（`plant ⊂ Host-assigned`，原 `plant = slice of
    #              Host-assigned` 在 4.12 in 轴里占 0.39 轴分数 → 必然越界）并左起 0.675。
    #   口径说明不丢：完整表述写在 `D.manifest["A"]["note"]` 与主图题注里。
    for idx, txt in ((ia, "all hosts"), (ip, "plant subset")):
        if not idx:
            continue
        f0, f1 = _af(y[idx[0]] + 0.40), _af(y[idx[-1]] - 0.40)
        ax.plot([x_brace, x_brace], [f0, f1], color="#888888", lw=0.8,
                solid_capstyle="butt", zorder=6, transform=ax.transAxes)
        for ff in (f0, f1):
            ax.plot([x_brace, x_brace - 0.020], [ff, ff], color="#888888", lw=0.8,
                    solid_capstyle="butt", zorder=6, transform=ax.transAxes)
        # 组名竖排贴在括号线**右侧**（v21 改）：左侧 0.627 会撞 `80,637` 的值标签
        #   （实测 7.7 × 2.3 px OVERLAP）→ 右侧 0.663 是净空带。
        ax.text(x_brace + 0.008, (f0 + f1) / 2, txt, fontsize=5.0, color="#666666",
                ha="left", va="center", rotation=90, zorder=6,
                transform=ax.transAxes)
    # 关键数写成可核对的等式（不用形容词），放在 scope 括号右侧、轴内
    #   ⚠ **文案必须够短**（半栏 A 轴只有 4.12 in）：见上方 `x_brace` 注释。
    ph = D.stage.get("06_HostPrediction", {}).get("plant_host")
    nphq, nres = D.plant_hq(), D.stage.get("08_Rescue", {}).get("rescue_rescued")
    nhigh, nhq = D.keep_high(), D.final_hq()
    note = "plant \u2282 Host-assigned"
    if ph:
        note += f"\n{ph:,} plant-host"
    if nphq:
        note += f"\n\u2192 {nphq} plant HQ"
    if nres:
        note += f"\n\u2192 {nres} rescued"
    if nhigh:
        note += f"\n\u2192 {nhigh} CheckV high"
    if nhq:
        note += f"\n\u2192 {nhq} final HQ"
    if ip:
        f_anchor = (_af(y[ip[0]] + 0.40) + _af(y[ip[-1]] - 0.40)) / 2
        ax.text(0.700, f_anchor, note, fontsize=4.8, color="#666666",
                ha="left", va="center", zorder=6, transform=ax.transAxes)
    D.manifest["A"] = {
        "source": "10_Reports/stage_summary.tsv + cluster_pipeline_reduction.tsv "
                  "+ cobra_summary.tsv + keep_summary.tsv + filter_per_sample.tsv "
                  "+ checkv_summary.tsv (Plant row)",
        "note": "13-step chain (v18). ALL-HOST block: Assembled → Viral sequences → "
                "Filtered candidates → COBRA queries → Merged → vOTU clusters → "
                "Taxonomically classified → Host-assigned. "
                "PLANT-subset block: Plant-host → Plant HQ (CheckV) → Rescued (VSI) "
                "→ CheckV high-confidence → Final HQ plant viruses. "
                "NOT a single-scale funnel: the plant block is a SLICE of Host-assigned "
                "(same colour on purpose), annotated by a right-side scope brace. "
                "⚠ 'Filtered candidates' from 02b_Filter/<sample>/uniprot/*.fasta counts "
                "(filter_summary.tsv is an all-zero dud in this run). "
                "⚠ 'Plant HQ' uses checkv_summary.tsv Plant row (denominator 1,090, "
                "≠ the 1,063 used elsewhere → known R6 split).",
        "n_items": len(items),
        "values": {k: v for k, v in items}}


def panel_B(ax, D: Data):
    """B 鉴定环节的多工具贡献（并集画成参考线）。

    ⚠ **已下线（v17, 2026-10-13，用户指令「BC 也删了呢」）**：B 不再是主图面板。
    下线理由 = 工具贡献属**方法学过程量**，与 A 漏斗的 `Viral sequences (10 tools)`
    同源，读者已由 A 得知规模；保留函数体供单跑/复用。当前 `PANELS` 无调用点。
    """
    tools = D.ident_tools()
    if not tools:
        return _empty(ax)
    labels = [t for t, _ in tools]
    vals = [n for _, n in tools]
    bc = colors(len(labels))
    _hbar_labels(ax, labels, vals, bc, log=True,
                 title="Identification tools", xlabel="Viral sequences (log)")
    union = D.ident_union()
    if union:
        ax.axvline(union, color="#444444", lw=0.9, ls="--", zorder=5)
        ax.annotate(f"union of all tools\n{union:,}",
                    xy=(union, 0.15), xytext=(union * 0.42, 0.55),
                    fontsize=5.0, color="#444444", ha="right", va="center",
                    arrowprops=dict(arrowstyle="->", lw=0.5, color="#666666"))
    D.manifest["B"] = {"source": "10_Reports/ident_summary.tsv (TOTAL row)",
                       "note": "All_Candidate column and zero-count tools excluded"}


def panel_C(ax, D: Data):
    """C vOTU 簇大小分布（附组装→聚类的缩减注记）。

    ⚠ **已下线（v17, 2026-10-13，用户指令「BC 也删了呢」）**：C 不再是主图面板。
    下线理由 = 簇统计属**中间产物 QC**，不承载结论（singleton 占比高是宏病毒组常规）。
    保留函数体供单跑/复用。当前 `PANELS` 无调用点。
    """
    if not D.cluster_size:
        return _empty(ax)
    rows = [(r["Size"], _int(r["Count"], 0) or 0) for r in D.cluster_size]
    rows = [r for r in rows if r[1] > 0]
    if not rows:
        return _empty(ax)
    labels = [k for k, _ in rows]
    vals = [v for _, v in rows]
    red = ""
    if D.cluster_red:
        d = {r["Stage"]: _int(r["Count"]) for r in D.cluster_red}
        red = (f" · {d.get('Input contigs', 0):,} contigs → "
               f"{d.get('vOTU clusters', 0):,} clusters")
    _hbar_labels(ax, labels, vals, [OKABE_ITO[4]] * len(labels), log=True,
                 title="vOTU cluster size distribution",
                 xlabel="Clusters (log)" + red, label_fs=6.0, xlabel_fs=5.2)
    D.manifest["C"] = {"source": "10_Reports/cluster_size_distribution.tsv "
                                 "+ cluster_pipeline_reduction.tsv"}


def panel_F(ax, D: Data):
    """⚠ **已下线（v16, 2026-10-13，用户指令）**：E「Novelty tier」不再进主图。

    不在 `PANELS`/`LAYOUT` 里 → **不会被调用**。函数体保留仅为便于将来恢复或单独出图。
    下线理由：新颖性分层与面板 A 漏斗的 novelty 分段语义重复，信息量低。
    """
    """F 新病毒等级（Novel Rank）。"""
    nov = D.novelty()
    if not nov:
        return _empty(ax)
    labels = [k for k, _ in nov]
    vals = [v for _, v in nov]
    sc = SEMANTIC["novelty"]
    keymap = {"Known": "Known", "Novel species": "NewSp",
              "Novel genus": "NewGe", "Novel family": "NewFa"}
    bc = [sc.get(keymap.get(l, ""), "#999999") for l in labels]
    x = np.arange(len(labels))
    ax.bar(x, vals, color=bc, edgecolor="#333333", linewidth=0.4, width=0.68)
    ax.set_xticks(x)
    ax.set_xticklabels([l.replace(" ", "\n") for l in labels], fontsize=6.0)
    ax.set_yscale("log")
    lo, hi = min(vals), max(vals)
    ticks = [t for t in (10, 100, 1000, 10000) if lo * 0.5 <= t <= hi * 2.5]
    if ticks:
        ax.set_yticks(ticks)
        ax.set_yticklabels([f"{t:,}" for t in ticks], fontsize=5.6)
    ax.set_ylim(lo * 0.45, hi * 3.0)
    ax.set_ylabel("vOTUs (log)", fontsize=6.4)
    for xx, v in zip(x, vals):
        ax.text(xx, v * 1.18, f"{v:,}", ha="center", fontsize=5.4, color="#333333")
    ax.set_title("Novelty tier", fontsize=7.2, fontweight="bold", loc="left")
    despine(ax)
    D.manifest["E"] = {"source": "stage_summary.tsv `└ Novel Rank` row"}


def panel_F_pair(ax, D: Data):
    """⚠ **已下线（v16, 2026-10-13，用户指令）**：F「Host group: predicted vs HQ」不再进主图。

    不在 `PANELS`/`LAYOUT` 里 → **不会被调用**；函数体保留便于恢复。
    """
    """F 宿主组：预测总数 vs 其中 HQ（合并原 G + H）。

    ⚠ 原 G 与 H 是**同一组 9 个类目、同一单位**的两张横条图，且 HQ 是总数的切片 →
    属硬重复，读者要看两遍、还得自己心算比例。这里合成**成对横条**：浅色 = 该宿主组
    预测总数，深色 = 其中 Complete+High-quality，右侧标 HQ 占本组的百分比。
    """
    rows = D.host_vs_hq()
    if not rows:
        return _empty(ax)
    labels = [h for h, _, _ in rows]
    tot = [n for _, n, _ in rows]
    hq = [q for _, _, q in rows]
    y = np.arange(len(labels))[::-1]
    hgt = 0.36
    ax.barh(y + hgt / 2, tot, height=hgt, color="#CCCCCC",
            edgecolor="#333333", linewidth=0.4, label="Predicted vOTUs")
    ax.barh(y - hgt / 2, hq, height=hgt,
            color=SEMANTIC["quality"]["High-quality"],
            edgecolor="#333333", linewidth=0.4, label="HQ (Complete + High-quality)")
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=6.0)
    xmax = max(tot) if tot else 1
    # ⚠ 三坑叠加（实测均报 OVERLAP）：
    #   ① 小百分比标签若画在**条形末端**（`q+eps`）→ 全挤在起点叠字；
    #   ② 改成"统一右列"后，若右列取 1.02·xmax，又会撞上最长条的值标签（`7,595`）
    #      —— 因为值标签从条末 `n+eps` 向右延伸，正好压到 1.02·xmax 那一列；
    #   ③ 图例放图内右下/中右都会碰到末排百分比标签。
    #   修法：xlim 显式拉宽到 1.75·xmax（给值标签 + 百比分列都留出净空），
    #        百分比列固定在 1.46·xmax；图例放到**图下方外侧**（不与任何文字同区）。
    ax.set_xlim(0, xmax * 1.75)
    x_pct = xmax * 1.46
    for yy, (h, n, q) in zip(y, rows):
        ax.text(n + xmax * 0.02, yy + hgt / 2, f"{n:,}", va="center", ha="left",
                fontsize=5.0, color="#333333")
        lab = f"{q}  ({q / n * 100:.1f}%)" if q > 0 else "0"
        ax.text(x_pct, yy - hgt / 2, lab, va="center", ha="left",
                fontsize=5.0, color="#333333")
    ax.set_xlabel("vOTUs", fontsize=6.4)
    # ⚠ 顶部留白：最上一行的值标签（如 Protist 的 7,595）会顶出轴框、压到标题。
    ax.set_ylim(-0.65, len(labels) - 0.35)
    ax.set_title("Host group: predicted vOTUs vs high-quality", fontsize=7.2,
                 fontweight="bold", loc="left")
    # 图例放**图外下方**：与任何标签都不同区，彻底避让（frameon=False 不占视觉）。
    ax.legend(loc="upper left", bbox_to_anchor=(0.0, -0.22), ncol=2,
              fontsize=4.8, frameon=False, borderpad=0.2, handlelength=1.4,
              columnspacing=1.2)
    despine(ax)
    pd_ = D.plant_denom()
    D.manifest["F"] = {
        "source": "host_distribution.tsv (total) + checkv_summary.tsv (HQ)",
        "note": "total from host_distribution.tsv; HQ = Complete+High-quality from "
                "checkv_summary.tsv; the two files use different plant denominators "
                f"(1,063 vs 1,090) — panel labels the total with the authoritative "
                f"host_distribution value (Plant={pd_})"}


def _sankey_data_dirs():
    """桑基数据层 `novelty_sankey_data/` 的候选目录（供读 manifest / contigs）。

    ⚠ `_HERE = Path(__file__).resolve()` 是**文件路径**不是目录 → 必须 `.parent`。
      本项目已因这个坑静默漏过资产一次（见 `_find_asset` docstring）。
    """
    return [
        _HERE.parent / "novelty_sankey_data",          # scripts/barbarum/novelty_sankey_data
    ]


def _find_asset(name: str, D: "Data", extra=None):
    """在 outdir / _build / _assets / 数据目录里找一个预渲染资产（png）。

    ⚠ 末位候选是 `<outdir>/data/06_HostPrediction/` —— 本地构建时把宿主映射与
      预渲染桑基 PNG 一起放在那儿，这样 **服务器上也能直接复用**（服务器无 Chrome）。
    ⚠⚠ 曾经这里写 `_HERE / "_assets" / name` —— 但 `_HERE = Path(__file__).resolve()`
      是**文件路径**不是目录 → 拼出 `<script>.py/_assets/<name>`（永不存在），
      末位候选形同虚设，新资产会被**静默漏掉**并退回旧路径（2026-10-13 实测）。
      → 改为 `_HERE.parent / "_assets"`（脚本同级 `_assets/`，与 docstring 一致）。
    """
    cands = []
    od = getattr(D, "_outdir", None)
    if od:
        cands += [od / name, od / "_build" / name, od / "_assets" / name,
                  od / "data" / "06_HostPrediction" / name]
    cands.append(_HERE.parent / "_assets" / name)
    for d in (extra or []):
        try:
            cands.append(Path(d) / name)
        except Exception:
            pass
    for c in cands:
        if c.is_file():
            return c
    return None


def panel_D_sankey(ax, D: Data):
    """D 分类谱系桑基。

    ⭐⭐ **现行口径（2026-10-13，用户「主图改通栏」）= 09b 五层桑基（121 单位）**：
        `Known/Novel → DNA/RNA/Viroid → family → genus → species`
        · 唯一权威数据源 = `scripts/barbarum/make_novelty_sankey_data.py`
          （`novelty_sankey_data/contigs.tsv`，121 条 = Known 40 + Novel 81）；
        · 唯一绘图实现 = `scripts/barbarum/fig_novelty_sankey.py`
          （`--embed-asset-plain` 产出的 `sankey_novelty_embed.png`）。
      ⚠⚠ **本面板不自行重画桑基**：历史上本项目吃过「同一面板两份标题副本/两份实现
        静默漂移」的亏（见 §9 硬规则）。这里只负责**贴图**，与旧 `sankey_flow_plant.png`
        旁路同构 —— 桑基几何、L5 着色、加粗规则全部由 fig_novelty_sankey.py 单点决定。
      ⚠ 为什么必须用 **plain（无发生率后缀）** 资产：发生率后缀 `(xx.x%)` 会把最长标签
        从 0.753 in 撑到 0.993 in，在通栏面板（可用 1.058 in）里余量只剩 6 %，
        服务器字形宽 ~12–14 %（硬规则）→ 必然溢出。发生率信息保留在**独立图**的注记里。
      ⚠ **与面板 A / I 是不同口径，不得混读**：A = 445 calls 漏斗、I = 1,063 植物宿主池、
        本面板 = 121 高置信单位（含 5 个类病毒物种，且类病毒按物种去重）。

    ── 历史保留：旧 6 层植物宿主桑基（Realm→…→Genus，1,063）──
    退化 / 回退路径仍保留在原处（`sankey_flow_plant.png`）：
      · 口径 = 植物宿主子集（1,063），对齐 `classification_sankey_plant.html`；
      · matplotlib 无原生桑基 → plotly 生成后 kaleido 导出静态 PNG 再贴图。
      · plotly/kaleido 缺失时**退化**为 Phylum 横条（`_panel_D_bars`），图面不空。
    """
    # ── ① 现行路径：09b 五层桑基（121）──
    #   ⚠ 资产查找顺序（`_find_asset`）：outdir → _build → _assets → 数据目录；
    #     生产端（fig_novelty_sankey.py --embed-asset-plain）把资产写到 outdir 下。
    pre09 = _find_asset("sankey_novelty_embed.png", D)
    if pre09 is not None:
        try:
            from PIL import Image
            img = np.asarray(Image.open(pre09))
            ax.imshow(img, aspect="auto")
            ax.axis("off")
            # ⚠⚠ 2026-10-13：标题数字**不得硬编码**（旧版写死 121 → 口径改成 61 后
            #    面板 D 会静默说谎）。改从桑基数据层 manifest 读实际单位数。
            _n_units = None
            _manp = _find_asset("manifest.json", D, extra=_sankey_data_dirs())
            if _manp is not None:
                try:
                    _n_units = int(json.loads(_manp.read_text(encoding="utf-8"))["total"])
                except Exception:
                    _n_units = None
            if _n_units is None:
                _n_units = 61   # 兜底（与数据层默认口径一致）
            ax.set_title(f"Novelty taxonomy flow \u00b7 {_n_units} high-confidence units",
                         fontsize=7.2, fontweight="bold", loc="left")
            ax.text(0.0, -0.055,
                    "Known/Novel \u2192 DNA/RNA/Viroid \u2192 family \u2192 genus \u2192 species"
                    "    (separate unit pool from panel A = 445 calls and"
                    " plant-host pool = 1,063)",
                    transform=ax.transAxes, fontsize=5.0, color="#555555",
                    ha="left", va="top")
            D.manifest["D"] = {
                "source": "scripts/barbarum/novelty_sankey_data/contigs.tsv "
                          "(via fig_novelty_sankey.py --embed-asset-plain)",
                "render": f"pre-rendered PNG: {pre09.name} (matplotlib, single "
                          f"implementation)",
                "scope": f"{_n_units} high-confidence units (Known + Novel); "
                         "\u26a0 2026-10-13 \u53e3\u5f84\uff1aNovel \u5df2\u540c\u5c5e\u53bb\u5197\u4f59"
                         "\uff08\u6bcf\u5c5e\u53ea\u5c55\u793a 1 \u6761\u4ee3\u8868\uff09"
                         "\uff0c\u5b8c\u6574\u660e\u7ec6\u89c1\u9644\u56fe novel_dedup_full.tsv",
                "note": "5-level flow (tier \u2192 genome type \u2192 family \u2192 genus \u2192 "
                        "species). Same unit pool as Fig_GojiVirome_NoveltySankey. "
                        "\u26a0 separate from panel A (445 calls) and the 1,063 plant-host "
                        "pool used by G/I; \u26a0 prevalence suffix intentionally omitted "
                        "here (kept in the standalone figure's caption)"}
            return
        except Exception:
            pass

    flow = D.sankey_flows(plant_only=True)
    if not flow:
        return _empty(ax)
    nodes, links = flow
    n_flow = sum(v for s, d, v in links if s < len(nodes) and nodes[s][0] == "Realm") \
        or sum(v for s, d, v in links if nodes[s][0] == "Realm")
    pd_ = D.plant_denom()
    scope = f"plant-host subset (n = {pd_:,})" if pd_ else "plant-host subset"
    warn = getattr(D, "_sankey_plant_warn", None)
    if warn:
        # 取不到宿主映射 → 退回全量，但**必须**在图面与 manifest 同时暴露口径变化
        scope = "ALL hosts (plant filter unavailable)"

    # ── 旁路：优先读预渲染的桑基 PNG ──────────────────────────────
    # ⚠ 服务器（202.119.189.246）的 kaleido **没有 Chrome** → plotly `to_image` 抛
    #   `RuntimeError: Kaleido requires Google Chrome`，且无外网装不上。
    #   故允许把**本地渲好的** `sankey_flow_plant.png` 放到 outdir 里，优先贴它。
    #   查找顺序同 `_locate`：outdir → <outdir>/_build → <outdir>/_assets → 脚本同级 _assets。
    pre = _find_asset("sankey_flow_plant.png", D)
    if pre is not None:
        try:
            from PIL import Image
            img = np.asarray(Image.open(pre))
            ax.imshow(img)
            ax.axis("off")
            # ⚠ v17：旁路分支的标题**必须与直渲分支逐字一致**（历史上这里是第二份
            #    副本，改标题时只改了直渲那处 → 图上仍是旧标题）。层级链改写成
            #    下方 text 注记，避免跨 2 列时侵入右侧 L 面板标题。
            ax.set_title(f"Taxonomy flow · {scope}", fontsize=7.2, fontweight="bold",
                         loc="left")
            ax.text(0.0, -0.045,
                    "Realm → Kingdom → Phylum → Class → Family → Genus",
                    transform=ax.transAxes, fontsize=5.4, color="#555555",
                    ha="left", va="top")
            man = {
                "source": "final_integrated_classification.tsv , filtered to Plant by "
                          "ensemble_host_summary.tsv (`Final_Host`)",
                "render": f"pre-rendered PNG: {pre.name} (plotly+kaleido on build host)",
                "note": "Realm→…→Genus; top taxa kept, rest 'Other'; empty -> "
                        "'Unclassified'. plant-only filter mirrors pipeline's "
                        "classification_sankey_plant.html. ⚠ NOT taxonomy_composition.tsv "
                        "(no lineage)"}
            if warn:
                man["warning"] = warn
            D.manifest["D"] = man
            return
        except Exception:
            pass

    try:
        import io

        import plotly.graph_objects as go
        from PIL import Image
    except Exception:
        return _panel_D_bars(ax, D)

    col_order = ["Realm", "Kingdom", "Phylum", "Class", "Family", "Genus"]
    xmap = {c: i for i, c in enumerate(col_order)}

    def _flow_of(i):
        return sum(v for s, d, v in links if s == i or d == i)

    pos = {}
    for c in col_order:
        ns = [i for i, (col, _) in enumerate(nodes) if col == c]
        ns.sort(key=lambda i: -_flow_of(i))
        for k, i in enumerate(ns):
            pos[i] = (xmap[c], -k)

    # ── 节点配色：每层一个主色系（OKABE_ITO），色内按名次微调明度 ──
    _BASE = {"Realm": "#0072B2", "Kingdom": "#009E73", "Phylum": "#D55E00",
             "Class": "#CC79A7", "Family": "#E69F00", "Genus": "#56B4E9"}
    GREY = "#B0B0B0"

    def _shade(hexc, k, n):
        """把基色按层内名次 k/n 做明度渐变（越靠前越深）。"""
        h = hexc.lstrip("#")
        r, g, b = int(h[0:2], 16), int(h[2:4], 16), int(h[4:6], 16)
        t = 0.0 if n <= 1 else 0.55 * (k / max(n - 1, 1))
        r = int(r + (255 - r) * t)
        g = int(g + (255 - g) * t)
        b = int(b + (255 - b) * t)
        return f"#{r:02X}{g:02X}{b:02X}"

    rank_in_col, col_counts = {}, {}
    for c in col_order:
        ns = [i for i, (col, _) in enumerate(nodes) if col == c]
        ns.sort(key=lambda i: -_flow_of(i))
        col_counts[c] = len(ns)
        for k, i in enumerate(ns):
            rank_in_col[i] = k

    node_col = []
    for i, (col, name) in enumerate(nodes):
        if name in ("Unclassified", "Other"):
            node_col.append(GREY)
        else:
            node_col.append(_shade(_BASE[col], rank_in_col[i], col_counts[col]))

    node_x = [xmap[col] / (len(col_order) - 1) for col, _ in nodes]
    node_y = [-pos[i][1] for i in range(len(nodes))]
    node_lab = [name for _, name in nodes]

    # ── 链接：渐变带（源色→目标色），同色系时用半透明 ──
    def _rgba(hexc, a):
        h = hexc.lstrip("#")
        return f"rgba({int(h[0:2],16)},{int(h[2:4],16)},{int(h[4:6],16)},{a})"

    src = [s for s, d, v in links]
    dst = [d for s, d, v in links]
    val = [v for s, d, v in links]
    link_col = []
    for s, d, _v in links:
        cs, cd = node_col[s], node_col[d]
        if cs == GREY or cd == GREY:
            link_col.append(_rgba(GREY, 0.30))
        else:
            link_col.append(_rgba(cs, 0.28))

    fig = go.Figure(go.Sankey(
        arrangement="snap",
        # ⚠ **不要手写 node.x / node.y**：plotly 的 `arrangement="snap"` 会在给了显式
        #   坐标后**忽略自动去重叠**，实测节点挤成一列、标签互相压住且首末列被裁。
        #   交给 plotly 自动布局（只按秩给顺序），颜色仍按层保留。
        node=dict(label=node_lab, pad=16, thickness=12,
                  line=dict(color="rgba(60,60,60,0.35)", width=0.4),
                  color=node_col),
        link=dict(source=src, target=dst, value=val, color=link_col),
        textfont=dict(size=11, family="Arial", color="#222222"),
    ))
    # ⚠ 标签必须留边距：plotly 把节点名画在节点**外侧**，margin 太小会把第一列/最后一列
    #   的名字整个裁掉（实测只看到色块、看不到名字）。
    #   v17：末列（Genus）名字较短且已在面板内 → 右边距收到 60，减少贴图后的右侧留白。
    fig.update_layout(margin=dict(l=95, r=60, t=10, b=10),
                      font=dict(size=11), paper_bgcolor="white",
                      plot_bgcolor="white")
    # ⚠ 画布宽高比要贴近**面板**宽高比，否则 `ax.imshow` 会留大片空白。
    #   v17：D 跨 2 列，面板实测 4.70 × 1.63 in → 比 2.89；故源图直接取 2.89
    #   （1760×610 ≈ 2.885），让 `imshow(aspect='equal')` 两头都不留白。
    png = fig.to_image(format="png", width=1760, height=610, scale=2)
    img = np.asarray(Image.open(io.BytesIO(png)))
    # 落盘：供**无 Chrome 的服务器**复用（见 panel_D_sankey 顶部旁路）。
    # ⚠ 文件名带 `_plant`：这是**植物子集**版，不能覆盖旧的全宿主版。
    try:
        if getattr(D, "_outdir", None):
            (D._outdir / "_assets").mkdir(parents=True, exist_ok=True)
            (D._outdir / "_assets" / "sankey_flow_plant.png").write_bytes(png)
    except Exception:
        pass
    ax.imshow(img)
    ax.axis("off")
    # ⚠ v17：原标题含 `(Realm → Kingdom → … → Genus)` 全程，跨 2 列时**右端会侵入
    #    右侧 L 面板的标题**（实测重叠）。层级信息桑基本身已用 6 列表达 → 标题精简为
    #    `Taxonomy flow · <scope>`，层级链改写在 text 注记里（不占标题行宽）。
    ax.set_title(f"Taxonomy flow · {scope}", fontsize=7.2, fontweight="bold",
                 loc="left")
    ax.text(0.0, -0.045, "Realm → Kingdom → Phylum → Class → Family → Genus",
            transform=ax.transAxes, fontsize=5.4, color="#555555",
            ha="left", va="top")
    man = {
        "source": "final_integrated_classification.tsv , filtered to Plant by "
                  "ensemble_host_summary.tsv (`Final_Host`)",
        "note": "Realm→Kingdom→Phylum→Class→Family→Genus; top taxa kept, rest 'Other'; "
                "empty cells -> 'Unclassified'. plant-only filter mirrors pipeline's "
                "classification_sankey_plant.html. ⚠ NOT taxonomy_composition.tsv "
                "(that file has no lineage, only per-rank totals)"}
    if warn:
        man["warning"] = warn
    D.manifest["D"] = man


def _panel_D_bars(ax, D: Data):
    """D 的退化版：Phylum 级 top 12 横条（plotly 不可用时）。

    ⚠ 与主路径同口径：**植物宿主子集**（1,063）——退化不代表换口径。
    """
    if not D.integrated:
        return _empty(ax)
    from collections import Counter
    pids = D.plant_contig_ids()
    rows = D.integrated if pids is None else \
        [r for r in D.integrated if (r.get("contig_id") or "").strip() in pids]
    scope = "all hosts" if pids is None else "plant-host subset"
    ph = Counter()
    for r in rows:
        t = (r.get("Phylum") or "").strip()
        ph[t if t and t not in ("NA", "-") else "Unclassified"] += 1
    top = ph.most_common(12)
    labels = [k for k, _ in top]
    vals = [v for _, v in top]
    _hbar_labels(ax, labels, vals, colors(len(labels)), log=True,
                 title="Taxonomy (Phylum, top 12)",
                 xlabel=f"vOTUs (log, {scope})",
                 label_fs=5.6)


def panel_G_plants(ax, D: Data):
    """⚠ **已下线（v16, 2026-10-13，用户指令）**：G「Plant-associated virus genera」不再进主图。

    不在 `PANELS`/`LAYOUT` 里 → **不会被调用**；函数体保留便于恢复。
    """
    """G 植物病毒属组成（按基因组型着色）。原 I 面板，版面重排后为 G。

    ⚠ **分母双写**：条值来自 `all_plant_viruses_genus_summary.tsv`（合计 1,090），
    而权威植物宿主池是 `plant_denom()` = 1,063（见遗留 R6）。两个数就是不一样，
    此处如实写进 xlabel，不掩盖、也不擅自改名掩盖。
    """
    rows = D.plant_genera(top=15)
    if not rows:
        return _empty(ax)
    labels = [g for g, _, _ in rows]
    vals = [n for _, n, _ in rows]
    gc = SEMANTIC["genome"]
    bc = [gc.get(gt, "#999999") for _, _, gt in rows]
    n_pool = sum(_int(r.get("n_contigs"), 0) or 0 for r in (D.plant_genus or []))
    pd_ = D.plant_denom()
    xlab = f"Contigs · genus table (n = {n_pool:,})"
    if pd_ and pd_ != n_pool:
        xlab += f"; plant-host pool = {pd_:,}"
    # ⚠ 值标签字号是**跨平台版式硬约束**：G 有 15 行、行距只有 ~7 px，`val_fs` 取默认
    #   5.6 时本机（matplotlib 3.11 / 不同 hinting）不碰，但服务器上字形略高 →
    #   相邻条的数值标签竖向叠 **0.8 px**（实测 14 对 OVERLAP，全在 ax6/G）。
    #   → G 的值标签单独降到 4.8，给跨平台字形留余量。
    _hbar_labels(ax, labels, vals, bc, label_fs=5.4, val_fs=4.8,
                 title="Plant-associated virus genera (top 15)",
                 xlabel=xlab)
    present = []
    for _, _, gt in rows:
        g = gt if gt else "none"
        if g not in present:
            present.append(g)
    handles = [Patch(facecolor=gc.get(g, "#999999"), edgecolor="#333333",
                     linewidth=0.4, label=(g if g != "none" else "unassigned"))
               for g in present]
    if handles:
        # ⚠ 图例放 ``lower right`` 时右缘会伸出轴框、撞到**相邻面板 H 的 y 刻度**
        #   （实测 `unassigned` ↔ `Fungal virus/contam.` 18.4×1.8px 跨面板重叠）。
        #   改放**图内右下偏内**（bbox_to_anchor x<1），保证整块不越轴框。
        ax.legend(handles=handles, title="Genome", loc="lower right",
                  bbox_to_anchor=(0.995, 0.02),
                  fontsize=5.0, title_fontsize=5.4, frameon=True,
                  borderpad=0.3, labelspacing=0.25)
    D.manifest["G"] = {"source": "10_Reports/all_plant_viruses_genus_summary.tsv",
                       "scope": "plant-host subset (contrast with panel D = all hosts)",
                       "n_shown": len(rows), "n_pool": n_pool,
                       "n_plant_host_authoritative": pd_}


def panel_H(ax, D: Data):
    """⚠ **已下线（v16, 2026-10-13，用户指令）**：H「Final verdict」不再进主图。

    不在 `PANELS`/`LAYOUT` 里 → **不会被调用**；函数体保留便于恢复。
    """
    """H 终判类别。原 J 面板，版面重排后为 H。"""
    v = D.verdicts()
    if not v:
        return _empty(ax)
    # 英文标签（图内不出现中文）
    EN = {"已知植物病毒": "Known plant virus", "新种候选": "Novel species cand.",
          "远缘病毒候选": "Distant virus cand.", "远缘(未定)": "Distant (undet.)",
          "未定": "Undetermined", "宿主基因污染": "Host-gene contam.",
          "昆虫病毒/污染": "Insect virus/contam.", "真菌病毒/污染": "Fungal virus/contam.",
          "噬菌体污染": "Phage contam."}
    labels = [EN.get(k, k) for k, _ in v]
    vals = [n for _, n in v]
    sc = SEMANTIC["verdict"]
    bc = [sc.get(k, "#999999") for k, _ in v]
    _hbar_labels(ax, labels, vals, bc, title="Final verdict (plant-virus candidates)",
                 xlabel="Contigs")
    D.manifest["H"] = {"source": "final_judgement_table.tsv (category column)"}


def panel_I_pair(ax, D: Data):
    """I 高质量基因组型 + 类病毒种（合并原 K + L）。

    ⚠ **"HQ" 一词在本图有两义**，极易误读：
      · 面板 A 的 `HQ vOTUs (all hosts)` = **200**（CheckV Complete + High-quality）；
      · 本面板的 "high-quality genomes" 合计 = **895**（DNA 44 + RNA 237 + Viroid 604
        + Unassigned 10）—— 它统计的是**基因组型构成**，不是同一类东西。
    故标题写全 "high-quality genomes"，与 A 的 "HQ vOTUs" 明确区分。

    左半：高质量基因组型（横向条）；右半：类病毒种（contigs vs samples 双条）。
    原 K（5 个类病毒种）单独成面板太单薄 → 与 L 并排，共享一个面板、各占一半宽度。
    """
    d = D.genome_viroid()
    if not d:
        return _empty(ax)
    # 用 inset 把一张轴分成左右两块
    from matplotlib.transforms import Bbox
    pos = ax.get_position()
    ax.set_position(pos)                      # 占位（inset 另起）
    ax.axis("off")

    # ⚠ v19：I 在通栏版式下只占 1 列（1.79 in 宽）→ 左右两块各约 0.80 in。
    #   ⚠⚠ 两块之间**必须留出 y 轴标签的净空**：右块的 y 标签是类病毒种名
    #     （`Posplviroid fusituberis` 等，@4.0pt 约 0.45 in 宽），若右块起点太靠左
    #     会直接压在左块的柱子上（2026-10-13 实测重叠）。
    #   → 左块给 0.00–0.42（不含右侧标签），右块从 0.52 起（留 0.10 轴宽 ≈ 0.18 in）。
    left = ax.inset_axes([0.00, 0.14, 0.42, 0.70])
    right = ax.inset_axes([0.62, 0.14, 0.38, 0.70])

    # ── 左：基因组型 ──
    gt = d.get("genome_type")
    if gt:
        labels = [k for k, _ in gt]
        vals = [v for _, v in gt]
        gc = SEMANTIC["genome"]
        bc = [gc.get(l, "#999999") for l in labels]
        x = np.arange(len(labels))
        left.bar(x, vals, color=bc, edgecolor="#333333", linewidth=0.4, width=0.66)
        left.set_xticks(x)
        left.set_xticklabels(labels, fontsize=4.2, rotation=0)
        left.set_yscale("log")
        lo, hi = min(vals), max(vals)
        tk = [t for t in (1, 10, 100, 1000) if lo * 0.5 <= t <= hi * 2.5]
        if tk:
            left.set_yticks(tk)
            left.set_yticklabels([f"{t:,}" for t in tk], fontsize=4.0)
        left.set_ylim(lo * 0.5, hi * 3.0)
        left.set_ylabel("genomes (log)", fontsize=4.6)
        for xx, v in zip(x, vals):
            left.text(xx, v * 1.2, f"{v:,}", ha="center", fontsize=4.0, color="#333333")
        left.set_title("Genome type", fontsize=5.2, fontweight="bold", loc="left")
        despine(left)
    else:
        left.axis("off")

    # ── 右：类病毒种 ──
    vr = d.get("viroid")
    if vr:
        # ⚠ 面板宽仅 1.79 in → 类病毒种名必须**大幅截断**（原 24 字符 @5.4pt 就超出）。
        #   截到 14 字符 + 省略号；完整名留给独立面板图。
        short = [(s if len(s) <= 14 else s[:13] + "\u2026") for s, _, _ in vr]
        y = np.arange(len(vr))[::-1]
        h = 0.34
        right.barh(y + h / 2, [r[1] for r in vr], height=h, color=OKABE_ITO[2],
                   edgecolor="#333333", linewidth=0.4, label="Contigs")
        right.barh(y - h / 2, [r[2] for r in vr], height=h, color=OKABE_ITO[1],
                   edgecolor="#333333", linewidth=0.4, label="Samples")
        right.set_yticks(y)
        right.set_yticklabels(short, fontsize=3.8)
        right.margins(x=0.20)
        right.set_xlabel("Count", fontsize=4.6)
        right.set_title("Viroid species", fontsize=5.2, fontweight="bold", loc="left")
        sci_legend(right, loc="lower right", fontsize=4.0)
        despine(right)
    else:
        right.axis("off")

    D.manifest["I"] = {"source": "hq_genome_type_summary.tsv + Viroid.species_info.tsv",
                       "note": "merged original K+L. ⚠ 'high-quality genomes' here "
                               "totals 895 and is NOT the same as panel A's "
                               "'HQ vOTUs = 200' (CheckV Complete+High-quality)"}


# ══════════════════════════════════════════════════════════════════
# 参考图（dataViz `data_report.R`）借用的两个面板
#   · J = 参考图 C：nt identity × coverage 散点（+ 边缘分布）
#   · K = 参考图 F：按科的 nt identity 箱线（横向）
# 判据：这两个面板补的是**我们原本 9 个面板全缺的"证据质量分布"**——
#   原图只有计数类条形图，读者看不到"这些判定到底有多可靠"。
# ══════════════════════════════════════════════════════════════════
def panel_J_scatter(ax, D: Data):
    """J identity vs coverage 散点（参考图 C）。

    参考图用 `geom_point(alpha=0.5, size=4)` + `ggMarginal(type="boxplot")`。
    matplotlib 无 ggMarginal → 用 `inset_axes` 手工加顶部/右侧边缘箱线，
    否则"两维分布"读不出单维集中度。
    """
    got = D.ident_vs_cov()
    if not got:
        return _empty(ax)
    pts, xlab, ylab = got
    x = np.array([p[0] for p in pts], dtype=float)
    y = np.array([p[1] for p in pts], dtype=float)

    # ── 按 verdict 上色（同 SEMANTIC.verdict，与面板 H 同源）──
    # ⚠ `rescue_evidence_scored.tsv` 的 `verdict` 列取值是**英文**（KEEP/DROP/REVIEW），
    #   而 `final_judgement_table.category`（面板 H）是中文。两者都要能落到同一色表，
    #   否则同一张图里"保留"会有两种颜色。
    vc = SEMANTIC["verdict"]
    vmap = {"KEEP": "KEEP", "keep": "KEEP", "已知植物病毒": "KEEP", "新种候选": "KEEP",
            "DROP": "DROP", "drop": "DROP", "宿主基因污染": "DROP",
            "REVIEW": "REVIEW", "review": "REVIEW", "远缘病毒候选": "REVIEW",
            "远缘(未定)": "REVIEW", "未定": "REVIEW"}
    groups = {}
    for xx, yy, v in pts:
        k = vmap.get(v, "REVIEW")
        groups.setdefault(k, ([], []))
        groups[k][0].append(xx)
        groups[k][1].append(yy)

    # 边缘箱线的轴（在主轴上开 inset，避免改变主轴刻度范围）
    # ⚠⚠ **inset 的尺寸必须按「英寸」算、不能写死轴分数**（2026-10-13 v21 实测）：
    #   轴高从 v19 的 1.6 in → v21 的 3.74 in 后，写死的 `h=0.12` 轴分数 = **0.45 in**，
    #   边缘箱线变成一条横贯整图顶部的粗带，**直接压住 suptitle 与 J 标题**。
    #   → 顶部条高固定 **0.22 in**、右侧条宽固定 **0.30 in**，再换算成轴分数。
    # ⚠ 三档宽度策略：
    #   · **通栏（≥ 6.0 in）** → 外置（参考图 `ggMarginal` 观感，v20 通栏 J 用过）；
    #   · **半栏（3.0–6.0 in）** → 内嵌（v21 的 J = 3.44 in）；
    #   · **窄（< 3.0 in）** → 内嵌更紧凑（v19 的 J 仅 1.8 in）。
    try:
        _pos = ax.get_position()
        _fw, _fh = ax.figure.get_size_inches()
        _w_in = _pos.width * _fw
        _h_in = _pos.height * _fh
    except Exception:
        _w_in, _h_in = 5.0, 3.0
    _outer = _w_in >= 6.0
    _bar_h_ax = 0.22 / max(_h_in, 0.1) if _h_in else 0.10     # 顶部条高（轴分数）
    _bar_w_ax = 0.30 / max(_w_in, 0.1) if _w_in else 0.10     # 右侧条宽（轴分数）
    # ⚠⚠ **内嵌档改为「轴外贴顶」（非「轴内压顶」）**（2026-10-13 v21 实测）：
    #   首版内嵌把箱线放进轴内 0.885–0.985 → **压住 y=100 的散点**（coverage 满值那排）。
    #   现在统一贴轴外，高度只有 0.22 in（v20 通栏时写死的 0.14 轴分数 = 0.45 in 才会盖 suptitle），
    #   配合 `top=0.915` 的净空 → 既不影响点云、也不撞标题。
    ax_top = ax.inset_axes([0.0, 1.015, 0.86, _bar_h_ax])
    ax_rgt = ax.inset_axes([1.015, 0.0, _bar_w_ax, 1.0])

    for k in sorted(groups):
        gx, gy = groups[k]
        ax.scatter(gx, gy, s=13, alpha=0.72, edgecolors="#333333", linewidths=0.25,
                   color=vc.get(k, "#999999"), label=f"{k} (n={len(gx)})", zorder=3)
    # 参考图的"可靠区"参考线（identity 80 %）——**只画竖线**：参考图里那两条灰线
    # 是 x/y 各一条，但本数据 y（coverage）在 0–100 全域铺开，画横线会横穿整个点云、
    # 与数据点抢视觉；identity 的 80 % 才是判"是不是同种"的判据 → 只留竖线。
    #   （竖线本体在 set_xlim 之后画，见下方；这里只记判据。）

    ax.set_xlabel(xlab, fontsize=6.2)
    ax.set_ylabel(ylab, fontsize=6.2)
    lo_x, hi_x = max(0, x.min() - 6), min(101, x.max() + 6)
    lo_y, hi_y = max(0, y.min() - 6), min(101, y.max() + 6)
    ax.set_xlim(lo_x, hi_x)
    ax.set_ylim(lo_y, hi_y)
    ax.tick_params(labelsize=5.4)
    # 参考线（必须在 set_xlim/set_ylim 之后画，见上方注释）
    ax.axvline(80, color="#999999", lw=0.5, ls=":", zorder=1)
    # ⚠ "80% identity" 竖排注记：**必须放在 set_ylim 之后**，否则 ylim 还是默认
    #   (0,1) → 注记落在画布外（实测被裁）。且注记 y 取 lo_y 而非 0，与点云底对齐。
    ax.text(80.7, lo_y, "80% identity", fontsize=4.2, color="#999999",
            rotation=90, va="bottom", ha="left")
    ax.set_title(f"Identity vs coverage (n = {len(pts):,})", fontsize=7.2,
                 fontweight="bold", loc="left")
    ax.legend(loc="lower left", fontsize=4.6, frameon=True, borderpad=0.25,
              labelspacing=0.22, handletextpad=0.3, framealpha=0.9)
    despine(ax)

    # ── 边缘箱线（展示单维分布；与主图同坐标轴范围）──
    for a, vals, orient, lo, hi in ((ax_top, x, "h", lo_x, hi_x),
                                    (ax_rgt, y, "v", lo_y, hi_y)):
        _box(a, [vals], horizontal=(orient == "h"), widths=0.55,
             patch_artist=True, showfliers=False,
             medianprops=dict(color="#333333", lw=0.7),
             boxprops=dict(facecolor="#B0B0B0", edgecolor="#666666", lw=0.4),
             whiskerprops=dict(color="#666666", lw=0.4),
             capprops=dict(color="#666666", lw=0.4))
        if orient == "h":
            a.set_xlim(lo, hi)
            a.set_xticks([])
            a.set_yticks([])
        else:
            a.set_ylim(lo, hi)
            a.set_xticks([])
            a.set_yticks([])
        for sp in a.spines.values():
            sp.set_visible(False)
    # 说明文字：贴**箱线右端外侧**（与箱线同 y 带），**绝不用 `set_title`**
    #   （inset 起于轴顶上方 → 其 title 会与主轴标题撞，v20/v21 两次实测）。
    ax_top.text(1.02, 0.5, "Distribution (box = IQR)", fontsize=4.4,
                color="#666666", ha="left", va="center",
                transform=ax_top.transAxes)

    D.manifest["J"] = {
        "source": "10_Reports/rescue_evidence_scored.tsv",
        "note": "nt_pident × nt_qcovs, only rows with BOTH > 0; colour = verdict. "
                "⚠ NOT All_plant.viruses_info.tsv (its `confidence` is an evidence "
                "score, not a % identity)"}


def panel_K_family(ax, D: Data):
    """K 按科的 nucleotide identity 箱线（参考图 F）。

    参考图 `geom_boxplot()` + `coord_flip()`，按中位数重排。
    ⚠ 退化守卫：`n < 3` 的科不画（箱线退化成一条线、还会把 y 轴标签挤满）。
    """
    rows = D.family_identity(min_n=3, top=14)
    if not rows:
        return _empty(ax)
    labels = [f for f, _ in rows]
    data = [v for _, v in rows]
    y = np.arange(len(labels))[::-1]
    bp = _box(ax, data, horizontal=True, positions=y, widths=0.62,
              patch_artist=True, showfliers=False,
              medianprops=dict(color="#333333", lw=0.9),
              boxprops=dict(facecolor=OKABE_ITO[6], edgecolor="#333333", lw=0.4),
              whiskerprops=dict(color="#333333", lw=0.5),
              capprops=dict(color="#333333", lw=0.5))
    # 每个箱左侧标 n（读者要判断箱线的可信度）
    for yy, (f, v) in zip(y, rows):
        ax.text(max(0.5, min(v) - 3.0), yy, f"n={len(v)}", va="center", ha="right",
                fontsize=4.0, color="#666666")
    ax.set_yticks(y)
    # ⚠ v19：K 在通栏版式下只占 1 列（1.79 in 宽）→ 科名字号收到 4.6（原 5.6 会挤）
    ax.set_yticklabels(labels, fontsize=4.6)
    ax.set_xlabel("Nucleotide identity (%)", fontsize=5.2)
    ax.tick_params(axis="x", labelsize=4.6)
    # ⚠ **不要写死 `set_xlim(0, 102)`**：本数据所有科的中位数都在 78–100 %，
    #   固定 0–102 会让左半张图全是空白（实测左 2/3 无墨迹），箱线被压成细条。
    #   改为按实际数据范围留边（但仍**必须在轴上保留 0 % 附近的锚**：只画 75–100
    #   容易读成"全部高相似"，而 80 % 参考线才是判据）。
    allv = [v for _, vs in rows for v in vs]
    lo_d, hi_d = min(allv), max(allv)
    span = max(hi_d - lo_d, 8.0)
    x0 = max(0.0, lo_d - 0.22 * span)
    x1 = min(102.0, hi_d + 0.10 * span)
    if x1 - x0 < 30:                       # 数据极窄时退回半宽视图（不然看不出差异）
        x0, x1 = max(0.0, (lo_d + hi_d) / 2 - 15), min(102.0, (lo_d + hi_d) / 2 + 15)
    ax.set_xlim(x0, x1)
    # 参考图的对照：80 % 常被当作"同种"阈值
    if x0 < 80 < x1:
        ax.axvline(80, color="#999999", lw=0.5, ls=":", zorder=1)
    ax.set_ylim(-0.6, len(labels) - 0.4)
    ax.set_title("Identity by virus family (n ≥ 3)", fontsize=7.2, fontweight="bold",
                 loc="left")
    despine(ax)
    D.manifest["K"] = {
        "source": "10_Reports/rescue_evidence_scored.tsv (tax_family × nt_pident)",
        "note": "families with n<3 dropped (degenerate boxplot); sorted by median; "
                "dotted line = 80% identity reference"}


def panel_L_domains(ax, D: Data):
    """L 保守域完整性饼图（对齐参考论文 Fig2-E）。

    参考论文 Fig2-E: `Distribution of complete, fragmented, and absent viral
    conserved domains in the characterized sequences`（`coord_polar` 饼图）。
    本面板用 `cdd_evidence_report.tsv` 的 `evidence` 列三分：
        Complete（PASS_CORE/PASS_CORE_FAM）/ Fragment（PASS_VIRAL）/ absent（其余）。
    ⚠ 退化守卫：只剩一类时饼图无信息 → 不画（`_empty`）。
    """
    got = D.domain_completeness()
    if not got:
        return _empty(ax)
    rows, n_tot = got
    if len(rows) < 2:
        return _empty(ax)
    labels = [k for k, _ in rows]
    vals = [v for _, v in rows]
    sc = SEMANTIC["quality"]
    cols = [sc.get(k, "#999999") for k in labels]
    wedges, _texts, autotexts = ax.pie(
        vals, colors=cols, startangle=90, counterclock=False,
        autopct=lambda p: f"{p:.0f}%", pctdistance=0.68,
        wedgeprops=dict(edgecolor="white", linewidth=0.8))
    for t in autotexts:
        t.set_fontsize(5.2)
        t.set_color("white")
        t.set_fontweight("bold")
    ax.axis("equal")
    # ⚠⚠ 图例**不能挂饼图右侧**：面板宽仅 1.79 in（通栏版式下 L 只占 1 列），
    #   右侧图例会直接压到相邻 K 面板的 y 轴标签上（2026-10-13 实测重叠）。
    #   → 改挂**饼图下方**、两列排（`ncol=2`），轴已在 `build_panels` 里上移让出底部。
    leg = [f"{k}  ({v:,})" for k, v in zip(labels, vals)]
    ax.legend(wedges, leg, loc="upper center", bbox_to_anchor=(0.5, -0.02),
              fontsize=5.0, frameon=False, handlelength=1.0, handleheight=1.0,
              labelspacing=0.25, borderpad=0.0, columnspacing=0.8, ncol=2)
    ax.set_title(f"Conserved-domain completeness (n = {n_tot:,})",
                 fontsize=6.4, fontweight="bold", loc="left", x=0.0, pad=3)
    D.manifest["L"] = {
        "source": "10_Reports/cdd_evidence_report.tsv (`evidence` column)",
        "note": "Complete = PASS_CORE/PASS_CORE_FAM (family core domains all hit); "
                "Fragment = PASS_VIRAL (viral domain hit, family core missing); "
                "No conserved domain = AMBIGUOUS/REVIEW. "
                "⚠ NOT checkv_summary.tsv (that is GENOME completeness, not domain "
                "completeness)"}


def panel_M_targets(ax, D: Data):
    """M 病毒 contig 的比对靶标（对齐参考论文 Fig2-B）。

    参考论文 Fig2-B: `Alignment targets of viral contigs`。
    ⚠ 我们的类目是从 `hmm_ct3_hits.tsv` 的 `ct3_best_viral_desc` **按关键词归类**得到
      （派生量），归不上的记为 unmapped 并在 xlabel 里如实报数。
    """
    got = D.protein_targets()
    if not got:
        return _empty(ax)
    rows, seen, unmapped = got
    labels = [k for k, _ in rows]
    vals = [v for _, v in rows]
    _hbar_labels(ax, labels, vals, colors(len(labels)),
                 title="Alignment targets of viral contigs",
                 xlabel=f"HMM best-hit domain (n = {seen:,}; "
                        f"{unmapped:,} unclassified)",
                 label_fs=5.6, val_fs=5.2)
    D.manifest["M"] = {
        "source": "10_Reports/hmm_ct3_hits.tsv (`ct3_best_viral_desc`)",
        "note": "derived by keyword rule (RdRp/RT/Capsid/Nucleoprotein/Helicase/"
                "Protease/Movement/Replicase/Integrase); NOT a curated annotation. "
                "targets with n<2 dropped; unclassified count reported in xlabel"}


PANELS = [
    ("A", panel_A),
    ("D", panel_D_sankey),
    ("I", panel_I_pair),
    ("J", panel_J_scatter), ("K", panel_K_family), ("L", panel_L_domains),
    ("M", panel_M_targets),
]

# 组合图布局（**v21, 3 面板 / 2 行 —— A 与 J 并排一行、D 独占一行**）。
#   用户 2026-10-13 晚（二次指令）：「AJ一行 呢 D 一行」
#   → 在 v20（A/J/D 三行全通栏）基础上改为 **上排 A|J 并排、下排 D 通栏**。
#   ⚠⚠ **列宽比必须偏向 A（漏斗）**：A 的 13 条横条 + 右侧 scope 括号 + note 极吃横向空间，
#     J 是散点、3 in 宽即舒展。实测 20 列网格下：
#       11:9 → A 4.34 × 3.74 in、J 3.44 × 3.74 in（J 略过 3.2 in → 走外置箱线）
#       1:1  → A 3.40 in（**太方**，13 条横条会被挤扁、右侧 note 无净空）✗
#     故取 **A:J = 11:9**（≈55 %:45 %）。网格用 **20 列**（11+9）保证比例可精确表达。
#   ⚠ D 行权单独反解（历史判据 = L5 标签 ≥ 4.2 pt）；本轮 61 单位版本 L5 37 个节点。
#   ⚠ **字母标识保持稳定（跳号 A/D/I/J/K/L/M）**：文档 §6.2c 对照表、manifest key、
#     `WIDE` 集合、`panel_letter` 分支全部锚定这些字母 → 删减/重排**不重编字母**。
_GRID_COLS = 20                     # 上排 11(A) + 9(J)；下排 D 跨满 20
LAYOUT = [
    ("A", 0, 0, 1, 11),
    ("J", 0, 11, 1, 9),
    ("D", 1, 0, 1, _GRID_COLS),
]
# ⚠⚠ **`LAYOUT` 的行号必须与 `height_ratios` 的下标严格对应**（2026-10-13 实测踩坑）：
#   首版把 5 行权写成 `[1.05, 1.30, 1.15, 0.90, 2.10]`（以为末位 2.10 给通栏的 D），
#   但 `LAYOUT` 里 D 在 **row=1**（权 1.30）、row=4 **根本没人用** →
#   D 面板只拿到 **1.51 in**（不是 2.97）、末行空白被 `bbox_inches="tight"` 裁掉，
#   主图高度从 12.6 回落到 9.31 in，而**没有任何报错**（静默！）。
#   v21：上排（A+J）需要足够行高（A 13 条横条、J 散点+边缘箱线）；D 行按 L5 可读下限。
#   ⭐⭐ v21.1（2026-10-13，用户「再宽高些」→ 澄清「D 面板 L1–L2 块再高些，不用宽了」）：
#     行权 `[1.28, 1.00]` → **`[0.95, 1.00]`**、画布高 `10.0` → **`11.5` in**。
#     ⚠⚠ **块高只能靠「加高 D 行」来放大，不能靠字号/块宽** —— 桑基的块高恒 =
#       `n × (画布高 / Y)`，而 `Y` 由 gap 地板（21.94 px，物理量）决定、**与画布高无关**
#       ⇒ **画布越高，单位像素线性增长**（详见 `fig_novelty_sankey.py::solve_geometry`
#       的实测结论）。故「让 L1–L2 更高」= 抬高 D 行本身。
#     ⚠ 行权刻意**压到 0.95**（<1.00）：这样 D 行增幅远大于上排 —— 实测
#       上排 4.40 → **4.32 in**（A/J 基本不变）、D 3.43 → **4.55 in（+33 %）**。
#       若沿用 1.28，D 要 4.55 in 就得把画布撑到 17.4 in（上排跟着涨到 6.9 in，浪费）。
#   ⭐⭐ v21.2（2026-10-13 21:32，用户「D图都看不清字体大些」）：
#     行权 `[0.95, 1.00]` → **`[0.75, 1.00]`**、画布高 `11.5` → **`13.0` in**。
#     ⚠⚠⚠ **本轮查清了「D 面板字号物理上限」（先量后判，结论写死在 README）**：
#       · D 面板的**桑基带高**实测仅 **2.48 in**（= D 行 4.55 in × 0.545，余下是标题/注记）；
#       · L5 有 **37 行** → 行距上限 `2.48×300/37 = 20.1 px = 4.83 pt`
#         → 字号封顶 `4.83/1.30 = `**`3.71 pt`**（实测渲出仅 **~2.0 pt**，还叠加了各向异性缩放）。
#       ⚠ **抬画布 / 抬字号都不解决**：抬 D 行只线性抬带高；抬 `--all-fs` 反会抬 `Y`
#         （gap 地板）→ 单位像素更小（实测资产 6.5 pt 贴进面板只等效 2.0 pt）。
#       ⇒ **要 6 pt 需桑基带 ≥ 4.01 in**，而 D 行即便 6.18 in（figH=14）也只 3.37 in ⇒ **不可达**。
#     ⇒ 用户选「**只抬 D 行（有限改善）**」：行权压到 0.75、画布 13.0 in
#       → D 行 4.55 → **5.25 in（+15 %）**、桑基带 2.48 → **~2.86 in**；上排 4.32 → 3.94 in。
LAYOUT_ROW_RATIOS = [0.6, 1.0]
# ⚠ 画布高 13.0 in 是**反解**出来的（v21.2）：可用高度 `top-bottom = 0.893-0.038 = 0.855` →
#   上排权 0.75 + 下排 1.00 + hspace 0.42 → A/J 高 ≈ 3.94 in、D 高 ≈ 5.25 in。
#   （v21.1 为 11.5 in → 上排 4.32 / D 4.55；v21 为 10.0 in → 4.40 / 3.43。）
#   A 从 v19 的 ~1.6 in 提上来（用户要「A 展开高些」），J 从 1.80 in 宽提到 3.44 in 宽，
#   两者在 v21 里**同处上排、顶部对齐**，视觉上是一组。
COMPOSITE_FIGSIZE = (8.55, 14.0)


# ══════════════════════════════════════════════════════════════════
# 组装
# ══════════════════════════════════════════════════════════════════
def save(fig, outdir: Path, stem: str, dpi: int):
    outdir.mkdir(parents=True, exist_ok=True)
    png, pdf = outdir / f"{stem}.png", outdir / f"{stem}.pdf"
    fig.savefig(png, dpi=dpi, bbox_inches="tight", pad_inches=0.05,
                facecolor="white")
    fig.savefig(pdf, dpi=dpi, bbox_inches="tight", pad_inches=0.05,
                facecolor="white")
    return png, pdf


def build_composite(D: Data, outdir: Path, dpi: int):
    # v21：上排 **A|J 并排**（11:9，用户「AJ一行」）、下排 **D 通栏**（用户「D 一行」）。
    #   ⚠ 网格用 20 列（11+9）→ 比例可精确表达，D 跨满 20 列。
    #   ⚠ 桑基面板 D 高须满足 L5 标签可读下限（历史判据 ≥ 4.2 pt；本轮 37 个节点）。
    #   ⚠ **不要用「面板比 ≈ 源图比」来定行高比**：现行 D 的源资产是**扁平化**桑基
    #     （--embed-size 给足宽度），本体比与面板比**故意不等**，贴图端 `imshow(aspect="auto")`。
    fig = plt.figure(figsize=COMPOSITE_FIGSIZE)
    # ⚠⚠ **`top` 必须同时给「2 行 suptitle」和「行内面板标题」留净空**（2026-10-13 v21 实测）：
    #   v21 上排是 A|J 并排，J 的标题 `loc="left"` 从 **J 轴左缘**（画布中部）起排 →
    #   它正好落在 suptitle **第二行**的尾部 → 两者重叠（实测 "Identity vs coverage (n = 82)"
    #   压住 "host prediction and plant-virus recovery"）。
    #   修法：`top` 从 0.945 降到 **0.915**，把面板标题整体下移让出 suptitle 两行的净空；
    #   suptitle `y` 提到 0.986（两行从 0.986 往下排，末行底约在 0.945）。
    #   ⚠⚠ **0.915 仍不够**（`get_window_extent` 实测）：J 标题落在 **y frac 0.953**，
    #     与 suptitle 末行（0.945–0.985）只差 0.007 图高（≈20 px）→ 仍重叠。
    #   → 最终 `top = 0.893`（J 标题 → ≈0.931），留 0.027 图高净空。
    # ⚠⚠ v21.3：`hspace` 0.42 → **0.28**。原 0.42 是为「上排 x 轴标签 + D 标题」留空，
    #   实测 0.42 × 平均行高 = **2.08 in** 的空带（图上肉眼可见大片留白）→ 浪费。
    #   压到 0.28 后同页高下 D 行 6.18 → **6.56 in**、上排 3.71 → **3.94 in**（两者都变好）。
    #   ⚠ 判据：A/J 的 x 轴标签与 D 的面板标题不得相撞（已跑 `_check_discovery_overlap` 0 对）。
    gs = GridSpec(2, _GRID_COLS, figure=fig, hspace=0.28, wspace=0.60,
                  height_ratios=LAYOUT_ROW_RATIOS,
                  left=0.085, right=0.975, top=0.893, bottom=0.038)
    made = []
    for letter, r, c, rs, cs in LAYOUT:
        ax = fig.add_subplot(gs[r:r + rs, c:c + cs])
        fn = dict(PANELS)[letter]
        try:
            fn(ax, D)
        except Exception as e:                       # 单面板失败不拖垮整图
            traceback.print_exc()
            _empty(ax, f"panel {letter} failed: {type(e).__name__}")
        if ax.has_data() or ax.texts or letter == "I":
            # 面板 I 走 axis("off")+inset → has_data() 恒 False，须显式补字母
            panel_letter(ax, letter)
        made.append(letter)
    fig.suptitle(
        f"Virus discovery in the {D.virus_name} RNA virome: assembly, multi-tool "
        f"identification, taxonomy,\nhost prediction and plant-virus recovery",
        fontsize=8.5, fontweight="bold", y=0.986)
    png, pdf = save(fig, outdir, "Fig_GojiVirome_Discovery", dpi)
    plt.close(fig)
    return png, pdf, made


def build_panels(D: Data, outdir: Path, dpi: int):
    pdir = outdir / "panels"
    pdir.mkdir(parents=True, exist_ok=True)
    made = []
    # 横向条形图 / 桑基 / 散点 / 箱线 / 比对靶标给宽画布（≈ 通栏），纵向柱状图/饼图给半栏
    # v16：删掉 E/F/G/H → 同步移除已不存在的 "F"/"G"。
    # v17：删掉 B/C → 再移除 "B"/"C"（保留字母稳定，仅清理集合）。
    # v19：D 升通栏（3 列 7.61 in 宽）、M/I 同置底行（M 跨 2 列 4.66 in、I 跨 1 列 1.70 in）
    #      → 单面板画布尺寸必须**与组合图里的面板比例同源**，否则组合图与单面板观感不一致
    #        （且 D 的单面板若仍用旧 9.32x8.40 会与贴图的扁平资产不匹配）。
    WIDE = {"A", "D", "I", "J", "K", "M"}
    for letter, fn in PANELS:
        wide = letter in WIDE
        if letter == "D":
            # 通栏桑基：与组合图 D 面板同比例。
            # ⚠⚠ 2026-10-13 v21.3：组合图 D 行 5.25 → 6.18 in，**本尺寸必须同步**
            #   （组合图 D 的 **axes 高 = GridSpec 行高**，故单面板 axes 也要 6.18 in
            #    → canvas = 6.18/0.90 = 6.87 in）。
            #   ⚠ 上一版（3.35 in）是按错误的「桑基带」口径设的，比组合图矮一半。
            fig = plt.figure(figsize=(7.61, 6.87))
            ax = fig.add_axes([0.005, 0.035, 0.990, 0.900])
        elif letter == "M":
            # 比对靶标：组合图里跨 2 列（4.66 in 宽）；7 条横条不需通栏
            fig = plt.figure(figsize=(4.66, 1.90))
            ax = fig.add_axes([0.30, 0.22, 0.66, 0.68])
        elif letter == "K":
            # 科级箱线：组合图里占 1 列（1.79 in 宽）→ 画布同步收窄。
            # ⚠ y 轴要放 14 个科名 → 左侧留白给足（轴只占右侧 0.60）。
            fig = plt.figure(figsize=(2.70, 1.77))
            ax = fig.add_axes([0.40, 0.26, 0.58, 0.56])
        elif letter == "L":
            # 饼图：组合图里占 1 列（1.79 in 宽）。
            # ⚠⚠ 旧版 5.0 in 画布 + 图例挂 `bbox_to_anchor=(1.0, 0.5)`（饼**右侧**）
            #   → 缩到 1.79 in 后面板必被挤爆（2026-10-13 实测：图例压到 K 的 y 标签上）。
            #   → 图例改挂**饼下方**（两列排），轴上移收窄让出底部空间。
            fig = plt.figure(figsize=(2.05, 1.77))
            ax = fig.add_axes([0.02, 0.36, 0.96, 0.56])
        elif letter == "I":
            # 两轴并排：组合图里占 1 列（1.79 in 宽）→ 画布同步收窄。
            #   左右各约 0.90 in：左 = 4 条基因组型柱、右 = 5 个类病毒种双条。
            fig = plt.figure(figsize=(1.85, 1.77))
            ax = fig.add_axes([0.02, 0.14, 0.96, 0.78])
        elif letter == "J":
            # 散点带边缘箱线：主轴要留出右侧/顶部 inset 的空间
            fig = plt.figure(figsize=(5.6, 3.6))
            ax = fig.add_axes([0.13, 0.15, 0.72, 0.68])
        else:
            fig = plt.figure(figsize=((7.0, 3.2) if wide else (4.4, 3.2)))
            ax = fig.add_axes([0.30 if wide else 0.20, 0.20, 0.66, 0.70])
        try:
            fn(ax, D)
        except Exception as e:
            traceback.print_exc()
            _empty(ax, f"panel {letter} failed: {type(e).__name__}")
        if ax.has_data() or ax.texts or letter == "I":
            # ⚠ 面板 I 走 `axis("off")` + inset_axes，`has_data()` 恒为 False
            #   → 必须显式补字母，否则组合图与单面板都缺 "I" 标号。
            #   ⚠ v19：I 在组合图里只占 **1 列**（宽 1.70 in），轴左缘就在画布 2 % 处
            #     → 字母偏移必须很小，否则 bbox x0 为负（判 TEXT-OFF-CANVAS）。
            #   ⚠ v19：D 升为通栏（轴左缘在画布 8.5 % 处，与 A 同列起点）
            #     → 用默认偏移即可；M 跨 2 列（左缘同 A）同理。
            if letter == "I":
                panel_letter(ax, letter, dx=0.0, dy=1.05)
            elif letter in ("J", "K"):
                panel_letter(ax, letter, dx=-0.13 if letter == "K" else -0.145, dy=1.03)
            elif letter == "L":
                # ⚠⚠ 饼图轴几乎占满画布左缘（轴 x0 = 0.02 画布宽 ≈ 0.041 in，宽 1.968 in）
                #   → 字母偏移必须**接近 0**；且标题是 `loc="left"` 从轴左缘起排，
                #     "Conserved-domain completeness (n = …" 在 1.79 in 面板里几乎占满整行
                #     → 字母若与标题同 y 必然相撞（2026-10-13 实测 4.1 x 2.3 px OVERLAP）。
                #   修法：dy 提到标题**上方**（标题在 dy≈1.06 一带），字母给 1.17。
                panel_letter(ax, letter, dx=-0.005, dy=1.17)
            elif letter == "D":
                # ⚠⚠ v19：D 轴占满通栏画布（左缘 0.005 画布宽 ≈ 0.038 in，宽 7.53 in）
                #   → 偏移超过 `-0.005` 就会让 bbox x0 为负（判 TEXT-OFF-CANVAS）；
                #   ⚠ 但 dx=0 时字母又落回标题行（标题 `loc="left"` 从轴左缘起排）
                #     → 实测 11.0 x 2.9 px OVERLAP。
                #   修法：dx 贴轴左缘 + dy 抬到标题**上方**。
                panel_letter(ax, letter, dx=0.0, dy=1.12)
            elif letter == "M":
                panel_letter(ax, letter, dx=-0.34, dy=1.05)
            else:
                panel_letter(ax, letter, dx=-0.34 if wide else -0.20, dy=1.03)
        stem = f"panel_{letter}"
        save(fig, pdir, stem, dpi)
        plt.close(fig)
        made.append(letter)
    return made


def prune_stale_panels(outdir: Path, keep: list[str]) -> list[Path]:
    """删除 panels/ 下不属于 `keep` 的陈旧单面板图。

    ⚠ 必要性：`build_panels()` **只写不删** → 面板一旦下线（如 v16 删 E/F/G/H），
    旧 png/pdf 会原地留存，后续极易被误读为「仍在用」。与项目「退化守卫：跳过作图
    须删旧图」同一条纪律。
    仅删 `panel_<X>.{png,pdf}` 形态的文件，不动 `phylo_*` 等其它产物。
    """
    pdir = outdir / "panels"
    if not pdir.is_dir():
        return []
    valid = {f"panel_{x}.{ext}" for x in keep for ext in ("png", "pdf")}
    removed = []
    for f in sorted(pdir.glob("panel_*.*")):
        if f.name not in valid:
            f.unlink()
            removed.append(f)
    return removed


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", required=True,
                    help="discovery 输出根目录（含 10_Reports/）")
    ap.add_argument("--outdir", default=None,
                    help="输出目录（默认 <root>/10_Reports/S10_SCI_Discovery）")
    ap.add_argument("--virus-name", default="Lycium barbarum")
    ap.add_argument("--dpi", type=int, default=DEFAULT_DPI)
    ap.add_argument("--panels-only", action="store_true")
    ap.add_argument("--composite-only", action="store_true")
    a = ap.parse_args()

    root = Path(a.root).expanduser().resolve()
    if not root.is_dir():
        sys.exit(f"[fatal] --root not a directory: {root}")
    outdir = Path(a.outdir).expanduser().resolve() if a.outdir \
        else (root / "10_Reports" / "S10_SCI_Discovery")

    apply_sci_style(base_size=7.5, grid=False)

    D = Data(root, a.virus_name)
    D._outdir = outdir                       # 供 _find_asset 找预渲染资产
    print(f"[in]  root        = {root}")
    print(f"[in]  reports dir = {D.rep}")
    print(f"[in]  runs (denom)= {D.n_runs}")

    made = []
    if not a.panels_only:
        png, pdf, m = build_composite(D, outdir, a.dpi)
        print(f"[out] composite   = {png.name} / {pdf.name}  panels={''.join(m)}")
        made += m
    if not a.composite_only:
        m = build_panels(D, outdir, a.dpi)
        print(f"[out] panels      = {len(m)} files in {outdir / 'panels'}")
        # ⚠ 删除已下线面板的陈旧图（见 prune_stale_panels docstring）
        stale = prune_stale_panels(outdir, [k for k, _ in PANELS])
        if stale:
            print(f"[prune] removed {len(stale)} stale panel file(s): "
                  f"{', '.join(f.name for f in stale)}")

        # ⭐ 附图：新颖性分类**完整明细**（用户 2026-10-13：「详细的结果放附图」）
        #   主图桑基每属只展示 1 条代表 → 85 条 Novel 里 64 条被折叠；
        #   本附图把全部 85 条可视化（数据源 = 数据层产出的 novel_dedup_full.tsv）。
        #   ⚠ 失败不拖垮主流程（主图/面板已经落盘），只打印告警。
        try:
            _supp = _HERE.parent / "fig_novelty_sankey_supp.py"
            _spec2 = _ilu.spec_from_file_location("fig_ns_supp", _supp)
            _mod2 = _ilu.module_from_spec(_spec2)
            _spec2.loader.exec_module(_mod2)
            _dsrc = _HERE.parent / "novelty_sankey_data" / "novel_dedup_full.tsv"
            _mod2.build(_dsrc, outdir)
            print("[out] supplementary = Fig_GojiVirome_NoveltySankey_Supp.png/.pdf")
        except SystemExit as e:
            print(f"[warn] supplementary figure skipped: {e}")
        except Exception as e:                       # noqa: BLE001
            print(f"[warn] supplementary figure failed: {type(e).__name__}: {e}")

        # ⭐ 附图 2：L/K/I/M 四面板（用户 2026-10-13 晚：「LKIM 放附图吧」）
        #   v20 主图只留 A/J/D；L/K/I/M 移入 `Fig_GojiVirome_Discovery_Supp`。
        #   ⚠ 该附图**复用本模块的面板函数** → 必须传 root/outdir，由它自己
        #     `apply_sci_style` 后重跑（口径零漂移）。
        #   ⚠ 失败不拖垮主流程（主图/面板已经落盘），只打印告警。
        try:
            _supp2 = _HERE.parent / "fig_discovery_supp.py"
            _spec3 = _ilu.spec_from_file_location("fig_disc_supp", _supp2)
            _mod3 = _ilu.module_from_spec(_spec3)
            _spec3.loader.exec_module(_mod3)
            _mod3.build(root, outdir, a.dpi, a.virus_name)
            print("[out] supplementary = Fig_GojiVirome_Discovery_Supp.png/.pdf")
        except SystemExit as e:
            print(f"[warn] discovery supplementary skipped: {e}")
        except Exception as e:                       # noqa: BLE001
            print(f"[warn] discovery supplementary failed: {type(e).__name__}: {e}")

    man = outdir / "_build" / "fig_manifest.json"
    man.parent.mkdir(parents=True, exist_ok=True)
    man.write_text(json.dumps({"root": str(root), "virus_name": a.virus_name,
                               "n_runs": D.n_runs, "panels": D.manifest},
                              ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"[out] manifest    = {man}")


if __name__ == "__main__":
    main()
