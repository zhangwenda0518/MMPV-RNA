#!/usr/bin/env python3
"""
report_builder.py — 模块化系统发育分析汇总 HTML 报告

按 8 大模块 (data → phylogeny → popgen → recomb → select → time → geography → report)
分节, 顺序与 mmpv_common.io_layout.PHYLO_MODULES / phylo_pipeline.STAGE_MODULES 对齐。
每节展示: stage 状态徽章 + 该模块 summary 表 (csv/tsv/json → HTML table) + plots (png base64 内嵌)。
单文件零依赖输出 phylo_report.html (常驻左侧边栏导航: 模块状态 dot + 病毒锚点,
任何窗口宽度都不隐藏, 仅打印时隐藏)。
"""

import os
import json
import csv
import base64
import html as _html
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from datetime import datetime

# 2026-09-15 修复 (P0-A): HTML 报告与 phylo_summary.csv 共用同一取数入口。
# 修复前 HTML 读 metrics['beta'/'r_squared'/'p_value'] (TreeTime RTT),
# 而 CSV 读 td_mean_rate/pic_r/td_relaxed_clock_p (treedater/PIC) —— 同名列
# 不同来源, 同一病毒两份产物数字互相打脸。现两处都走 headline_metrics()。
from utils.summary_metrics import headline_metrics, fmt_metric

# 跨管线统一 I/O 布局: ⑤ per-virus 模块目录名随 MMPV_IO_LAYOUT 解析
try:
    from mmpv_common.io_layout import ph_dir
except ImportError:  # 脱离包环境直接运行时自举
    import os as _os, sys as _sys
    _REPO_ROOT = _os.path.dirname(_os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
    if _REPO_ROOT not in _sys.path:
        _sys.path.insert(0, _REPO_ROOT)
    from mmpv_common.io_layout import ph_dir

try:
    from mmpv_common.table_widgets import TW_BLOCK as _TW
except Exception:
    _TW = ""   # 共享表格控件不可达时静默降级为原生表格

# ─────────────────────────────────────────────────────────────
# 模块定义: 每模块的目录、状态 stage 集、summary 文件、plots glob
# 相对每个病毒 work_dir 解析
# ─────────────────────────────────────────────────────────────
MODULES = [
    {
        "key": "data", "name": "Data (prep + metadata + online + clean)",
        "dir": ph_dir("data"),
        # 2026-10-16 修 (#5): 徽章原先只有 prep/online —— 09-16 新增的 metadata
        # (data/metadata/metadata_std.csv) 与 clean (data/clean/clean.fasta) 两个
        # stage 产物一直在 data 模块里, 但报告上完全不显示其状态 (审计确认)。
        "stages": ["prep", "metadata", "online", "clean"],
        # data 模块本身无汇总表, 但 seqgrouper 的分组饼图/统计表落在这里
        "summaries": ["ncbi_harvest/seqgroup_*.csv"],
        "plots": ["ncbi_harvest/*.png", "ncbi_harvest/*.pdf", "ncbi_harvest/*.svg",
                  "**/*.pdf"],
        "notes": ["combined.fasta samples + reference sequences", "dates.csv / sample_metadata.csv", "ncbi_harvest public sequence collection"],
    },
    {
        "key": "phylogeny", "name": "Phylogeny (align + saturation + tree + splitstree)",
        "dir": ph_dir("phylogeny"),
        "stages": ["align", "saturation", "tree", "splitstree"],
        "summaries": ["saturation_out/saturation_report.csv", "qc/align_qc_report.tsv"],
        "plots": ["saturation_out/*.png", "qc/*.png", "tree_circular.png", "splitstree/*.png",
                  "gene_analysis/*.png", "cross_validation/*.png",
                  "qc/*.pdf", "tree_circular.pdf", "**/*.pdf"],
        # 2026-10-16 修 (#7): 多图模块 (比对/饱和/建树/分基因/交叉验证...) 给 60 张上限,
        # 原先统一 20 张会静默截断; 被截断的候选仍在报告里带链接列出。
        "max_imgs": 60,
        "notes": ["MAFFT alignment + IQ-TREE ML tree", "Xia substitution saturation test", "alignment QC"],
    },
    {
        "key": "popgen", "name": "Population genetics (popgen + host)",
        "dir": ph_dir("popgen"),
        "stages": ["popgen", "host"],
        "summaries": ["dnasp/results.tsv", "dnasp/ld_pairs.tsv", "location_pops.tsv",
                      # popgen_stats (2026-10-12): 多面板统计图的数据源 + EasyHap 借鉴三表
                      "popgen_stats.json",
                      "popgen_stats.private_haplotypes.tsv",
                      "popgen_stats.haplotype_enrichment.tsv",
                      "popgen_stats.similarity_clusters.tsv"],
        "plots": ["dnasp/figures/*.png", "host_analysis/*.png", "distance/*.png",
                  "popgen_metrics.png", "popgen_stats.png",
                  "*hap_state_heatmap.png", "**/*.pdf"],
        # extra_dirs (2026-10-12): 模块目录外的三棵产物树 (EasyHap 实跑 + MSN 网络)
        "extra_dirs": [
            {"dir": "haplo_outputs",
             "plots": ["*pub.png", "*map.png", "*.html"], "max_imgs": 8,
             "tables": ["network_stats.json"]},
            {"dir": "EasyHap_results/analyze",
             "plots": ["*HaplotypeNetwork.png", "*GroupStackedBar.png", "*GroupPieChart.png"],
             "max_imgs": 8,
             "tables": ["*HapSummary.tsv", "*PrivateHaplotypes.tsv"]},
            {"dir": "EasyHap_results/analyze_genes",
             "plots": ["*HaplotypeHeatmap.png", "*LD_r2_Heatmap.png",
                       "*HaplotypeNetwork.png", "*GroupStackedBar.png"],
             "max_imgs": 30},
        ],
        "notes": ["haplotype diversity + Tajima's D / Fu's Fs", "host differentiation analysis"],
    },
    {
        "key": "recomb", "name": "Recombination (rdp5)",
        "dir": ph_dir("recomb"),
        "stages": ["rdp5"],
        "summaries": ["rdp5/recombination_report.md", "rdp5/*.csv"],
        "plots": ["rdp5/*.png", "rdp5/figures/*.png", "**/*.pdf"],
        "notes": ["RDP5 seven methods + recombinant masking"],
    },
    {
        "key": "select", "name": "Selection (capheine)",
        "dir": ph_dir("select"),
        "stages": ["capheine"],
        # 2026-10-16 修 (#2): summaries 原先为空数组 —— capheine 跑完一张表都不渲染。
        # 现声明两个当前真实落盘口 (相对 select 模块目录):
        #   · capheine/drhip/combined_summary.csv —— DRHIP 基因级汇总
        #     (visual_codon_miner.py:249-256 读它; 管线登记 phylo_pipeline.py:917)
        #   · capheine/selection_out/selection_results.json —— selection_analysis
        #     parse 产物 (selection_analysis.py:552; 由 phylo_pipeline.py:4958-4971
        #     的 capheine 后处理写出)
        # 不声明: drhip/combined_sites.csv (逐位点, 行数大且表格渲染价值低)、
        # multiqc/capheine_summary.html (整页 HTML, 不属于表格/图件渲染路径)。
        # 告警门控见 missing_declared_summaries: 模块没跑 (无 stages_done 且无
        # .checkpoints/capheine.done) 时不告警。
        "summaries": ["capheine/drhip/combined_summary.csv",
                      "capheine/selection_out/selection_results.json"],
        "plots": ["capheine/**/*.png", "capheine/**/*.pdf"],
        # 2026-10-16 修 (#7): capheine 每基因一套图, 20 张必然截断 → 60。
        "max_imgs": 60,
        "notes": ["HyPhy FEL/MEME/BUSTED etc. + DRHIP site integration"],
    },
    {
        "key": "time", "name": "Time (clock + beast + gene_dating + treetime)",
        "dir": ph_dir("time"),
        # 2026-10-16 修 (#4): 徽章原先写 legacy 子步骤名 rtt/temporal —— 真实 stage
        # (STAGE_ORDER / STAGE_GROUPS) 是 clock: rtt/temporal/genes 已合并为 clock
        # (phylo_pipeline._LEGACY_CLOCK_SUBSTAGES), 老名字永远匹配不上 stages_done,
        # 报告里 clock 完成了也显示灰色。.done 回退判定对旧标记的兼容见
        # stage_badge 的 _LEGACY_STAGE_ALIASES。
        "stages": ["clock", "beast", "gene_dating", "treetime"],
        "summaries": ["clock_summary.tsv", "temporal_signal/drt_summary.csv",
                      "beast/merged/parameter_summary.csv",
                      "treetime_rtt/robust/robust_report.json"],
                      # 注: treetime stage 的 molecular_clock.txt 不进 summaries ——
                      # 该 stage 默认关, 没跑过的病毒会在报告里常驻 "declared but
                      # not found" 噪音; 产物本身仍被 plots/**.pdf 与文件树覆盖。
        "plots": ["treetime_rtt/*.png", "treetime_rtt/robust/*.png", "treedater_ltt/*.png",
                  "temporal_signal/real/*.png", "beast/postprocess/*.png",
                  "beast/merged/*.png", "gene_dating/*.png", "**/*.pdf"],
        "notes": ["Root-to-tip + DRT temporal signal", "BEAST Bayesian dating", "clock_summary — rate estimates across methods",
                  "treetime: TreeTime CLI timetree.nexus (prerequisite for the phymapr route, 2026-10-12)"],
    },
    {
        "key": "geography", "name": "Geography (phylogeo + geo + phymapr)",
        "dir": ph_dir("geography"),
        "stages": ["phylogeo", "geo", "geo_paths", "phymapr"],
        "summaries": ["phylogeography/merged/parameter_summary.json",
                      # 2026-10-12 实测对齐: tempmig 是 visualization/migration_bf.csv;
                      # phymapr 默认关 (未跑则不声明, 免 not-found) —— 跑了 phymapr/
                      # treetime 的病毒自行可见。
                      # 2026-10-16 修 (#3): geo_paths 产物已从历史运行目录
                      # absorb_test_20260916_full/ 归位到 geo_analysis/ —— 落盘口见
                      # phylo_pipeline.py:1052-1058 的产物登记 (_reg "pathway_qc" /
                      # "pathways_csv" / "dispersion_summary" / "timescale_consistency"),
                      # pathways/episodes 由 utils/transmission_paths.run_pathways 写,
                      # dispersion/velocity_table.csv 由 utils/diffusion_stats 写。
                      # 旧路径声明导致"表永远 not-found + 新表永远不出现"。
                      "visualization/migration_bf.csv",
                      "geo_analysis/spatiotemporal_counts.csv",
                      "geo_analysis/pathways/pathway_qc.json",
                      "geo_analysis/pathways/pathways.csv",
                      "geo_analysis/pathways/episodes.csv",
                      "geo_analysis/dispersion/dispersion_summary.json",
                      "geo_analysis/dispersion/velocity_table.csv",
                      "geo_analysis/timescale_consistency.json"],
        # visualization/ 里既有 SpreaD3 地图 (html) 也有发表级多面板图 (pdf)
        # 2026-10-16 修 (#6): 补 visualization/*.svg —— transmission_map.py:315 与
        # 多面板静态图在 png/pdf 之外还落 .svg, 原先无 glob 命中 → 静默漏收。
        "plots": ["visualization/*.png", "visualization/*.svg", "geo_analysis/*.png",
                  "tempmig/*.png",
                  "visualization/*.gif",
                  "phymapr/*.png", "phymapr/*.gif",
                  "**/*.pdf", "visualization/*.html", "geo_analysis/*.html"],
        "notes": ["CTMC+BSSVS ancestral geographic reconstruction", "Mantel genetic-geographic correlation", "temporal migration",
                  "geo_paths: event-level transmission pathways/episodes/three-tier confidence + dispersal stats (2026-09-16)",
                  "phymapr: TreeTime route BM-ASR + transmission map GIF (2026-10-12)"],
    },
    {
        # 2026-10-16 新增 (#1): 管道真实输出还有第 8 个模块 report/ —— 汇总表 +
        # 发表级三件套 + 采样分布图, 原先不在 MODULES 里, 报告永远看不到它们。
        # 落盘口 (相对 ⑤ per-virus work_dir, 目录名经 ph_dir('report') 解析,
        # standard 布局 = 08_report/):
        #   · phylo_summary.csv          ← generate_summary() (phylo_pipeline.py:6131-6142,
        #                                    root 级 = 本模块目录根, 不是 per-virus work_dir 根)
        #   · publication/*.pdf|png      ← publication_figures 后处理 (Fig1 TempEst /
        #                                    Fig2 MCC / Fig3 / Fig4 / Table S1; phylo_pipeline.py:6200)
        #   · sample_plot/*.pdf|png      ← sample_plot 后处理 (采样时间线 + 采样地图; :6170)
        # stages 留空: prep/report 不写 checkpoint (STAGE_CHECKPOINT_FILES['prep']=None,
        # report 由 main() 收尾无条件产出, 不是真实阶段分支) → 不给徽章;
        # 同时使 missing_declared_summaries 对 report 永不告警 (非 stage 产物)。
        "key": "report", "name": "Report",
        "dir": ph_dir("report"),
        "stages": [],
        "summaries": ["phylo_summary.csv"],
        "plots": ["publication/*.pdf", "publication/*.png",
                  "sample_plot/*.pdf", "sample_plot/*.png"],
        # 2026-10-16 修 (#7): 发表级三件套 + 采样图全是图, 20 张偏紧 → 60。
        "max_imgs": 60,
        "notes": ["phylo_summary.csv — headline metrics across viruses (same source as this page)",
                  "publication/: Fig1 TempEst / Fig2 MCC tree / Fig3 / Fig4 / Table S1",
                  "sample_plot/: sampling timeline + sampling map"],
    },
]

# ─────────────────────────────────────────────────────────────
# 可渲染产物白名单 (E1, 2026-09-16)
# ─────────────────────────────────────────────────────────────
# 历史坑: collect_plots 只认 `.suffix == ".png"`, 而管线地理/系统地理类产物
# 几乎全是 PDF (phylogeo_multipanel.pdf / mantel_test.pdf / tempmig.pdf /
# rrt_results.pdf / VirSpaceTime.pdf / sample_map.pdf / mcmc_diagnostics.pdf …)
# → 报告 Geography 段长期「有产物但空图」且不报任何错。
# 同一条 `.png` 过滤还把 MODULES['data'] 里写的 `ncbi_harvest/*.pdf` 变成死 glob。
#
# 现在按扩展名分派: 图片 → <img>, PDF → <embed>+下载链接 (浏览器原生 PDF 阅读器)。
PLOT_IMG_EXTS = (".png", ".jpg", ".jpeg", ".gif", ".svg", ".webp")
PLOT_PDF_EXTS = (".pdf",)
# 自包含交互图 (phylogeo_map.html / skyline.html …): 体积可达数 MB (内联 plotly.min.js),
# 不做内联, 报告里给**相对路径链接**——报告与产物同在 work_dir 下一棵树内, 链接有效。
PLOT_LINK_EXTS = (".html",)
PLOT_EXTS = PLOT_IMG_EXTS + PLOT_PDF_EXTS

# 同一 stem 多格式共存时选谁内联 (png 最通用; svg 次之; pdf 兜底, 体积最大; html 只做链接)
_EXT_PRIORITY = {".png": 0, ".svg": 1, ".jpg": 2, ".jpeg": 3,
                 ".gif": 4, ".webp": 5, ".pdf": 6, ".html": 9}
_MIME = {".png": "image/png", ".jpg": "image/jpeg", ".jpeg": "image/jpeg",
         ".gif": "image/gif", ".svg": "image/svg+xml", ".webp": "image/webp",
         ".pdf": "application/pdf"}

# 内联体积上限 (base64 膨胀 ~33%, 报告是单文件, 不能无限膨胀)
# 2026-10-17 改 (B): 图片上限 500 → 1200 KB —— 300dpi 的发表级 PNG 普遍
# 600-1100 KB, 卡在 500 时全被踢出内联分支 (报告只剩 PDF embed/下载链接)。
# PDF 上限 6000 KB 不变; 超限产物照旧进「未内联/被排除」清单, 不静默丢弃。
MAX_IMG_KB = 1200
MAX_PDF_KB = 6000

# 2026-10-16 修 (#7): 每模块内联图件上限改为可按模块覆盖 —— 原先全模块硬编码 20 张,
# 多图模块 (phylogeny / select / report) 稳定被截断。MODULES 条目可用 max_imgs
# 覆盖; 超出上限的候选不静默丢弃, 在报告里带下载链接列出 (上限 MAX_SKIP_LIST 条)。
MAX_IMGS_PER_MODULE = 20
MAX_SKIP_LIST = 30


# ─────────────────────────────────────────────────────────────
# 工具函数
# ─────────────────────────────────────────────────────────────

def embed_asset(path: str,
                max_img_kb: int = MAX_IMG_KB,
                max_pdf_kb: int = MAX_PDF_KB) -> Tuple[str, str, str]:
    """产物 → (kind, data_uri, reason)。

    kind ∈ {'img', 'pdf', ''}；不内联时 data_uri 为空且 **reason 非空** —— 
    调用方据此在报告里显式写出「未内联原因」，不静默丢产物。
    """
    if not path:
        return "", "", "empty path"
    p = Path(path)
    if not p.is_file():
        return "", "", "file not found"
    ext = p.suffix.lower()
    if ext not in PLOT_EXTS:
        return "", "", f"extension {ext or '(none)'} not in renderable whitelist {PLOT_EXTS}"
    kind = "pdf" if ext in PLOT_PDF_EXTS else "img"
    size_kb = p.stat().st_size / 1024
    cap = max_pdf_kb if kind == "pdf" else max_img_kb
    if size_kb > cap:
        return kind, "", f"{size_kb:.0f} KB > inline limit {cap} KB"
    try:
        with open(path, "rb") as f:
            b64 = base64.b64encode(f.read()).decode()
    except OSError as e:
        return kind, "", f"read failed: {e}"
    return kind, f"data:{_MIME[ext]};base64,{b64}", ""


def _rel_href(path: str, output_dir: Optional[str]) -> str:
    if output_dir:
        try:
            return os.path.relpath(path, output_dir).replace("\\", "/")
        except ValueError:
            pass
    return Path(path).name


def asset_card(path: str, output_dir: Optional[str] = None) -> Tuple[str, str]:
    """产物 → (HTML 卡片, 跳过原因)。不内联时仍给相对路径下载链接。"""
    name = Path(path).name
    label = _html.escape(Path(path).stem.replace("_", " ").title())
    if Path(path).suffix.lower() in PLOT_LINK_EXTS:
        href = _html.escape(_rel_href(path, output_dir), quote=True)
        return (f'<div class="img-card"><div class="img-label">{label}'
                f'<span class="tag">Interactive HTML</span></div>'
                f'<p class="asset-note">↗ <a href="{href}">Open in new tab {_html.escape(name, quote=True)}</a>'
                f'(single self-contained file, no external links)</p></div>'), ""
    kind, uri, reason = embed_asset(path)
    if uri:
        if kind == "pdf":
            return (f'<div class="img-card"><div class="img-label">{label}'
                    f'<span class="tag">PDF</span></div>'
                    f'<embed class="pdf-embed" type="application/pdf" src="{uri}" /></div>'), ""
        return (f'<div class="img-card"><div class="img-label">{label}</div>'
                f'<img src="{uri}" /></div>'), ""
    href = _html.escape(_rel_href(path, output_dir), quote=True)
    return (f'<div class="img-card off"><div class="img-label">{label}'
            f'<span class="tag warn">not inlined</span></div>'
            f'<p class="asset-note">⚠ {_html.escape(reason, quote=True)}<br>'
            f'<a href="{href}" download>Download {_html.escape(name, quote=True)}</a></p></div>'), reason


def collect_plots(vdir: str, patterns: List[str], max_imgs: int = MAX_IMGS_PER_MODULE,
                  skipped: Optional[List[Tuple[str, str]]] = None) -> List[str]:
    """按 glob 收集模块可渲染产物 (png/svg/jpg/gif/webp/pdf/html), 同 stem 取优先级最高者。

    max_imgs: 本模块内联上限 (MODULES 条目的 max_imgs 字段覆盖默认值, 见 #7)。
    skipped: 传入 list 时把每个「被排除的候选 + 原因」追加进去 (上限之外、
            同 stem 低优先级、扩展名不在白名单), 由调用方写进报告 —— 排除不等于没发生。
    """
    by_stem: Dict[Tuple[str, str], Path] = {}
    for pat in patterns:
        for p in sorted(Path(vdir).glob(pat)):
            if not p.is_file():
                continue
            ext = p.suffix.lower()
            if ext not in PLOT_EXTS and ext not in PLOT_LINK_EXTS:
                if skipped is not None:
                    skipped.append((str(p), f"extension {ext or '(none)'} not in renderable whitelist"))
                continue
            key = (p.parent.as_posix(), p.stem.lower())
            cur = by_stem.get(key)
            if cur is None:
                by_stem[key] = p
            elif _EXT_PRIORITY[ext] < _EXT_PRIORITY[cur.suffix.lower()]:
                if skipped is not None:
                    skipped.append((str(cur), f"same-stem format {ext} has higher priority, not embedded twice"))
                by_stem[key] = p
            elif skipped is not None:
                skipped.append((str(p), f"same-stem format {cur.suffix.lower()} has higher priority, not embedded twice"))
    out = [str(p) for p in sorted(by_stem.values(),
                                  key=lambda x: (x.parent.as_posix(), x.stem.lower()))]
    if len(out) > max_imgs:
        if skipped is not None:
            for _p in out[max_imgs:]:
                skipped.append((_p, f"beyond per-module limit max_imgs={max_imgs}"))
        out = out[:max_imgs]
    return out


def _stage_done_on_disk(vdir: str, stage: str) -> bool:
    """该 stage 是否在磁盘上留下 .done 标记 (与 stage_badge 的回退判定同一口径)。"""
    return bool(vdir) and os.path.exists(
        os.path.join(vdir, ".checkpoints", f"{stage}.done"))


def missing_declared_summaries(vdir: str, m: Dict, stages_done) -> List[str]:
    """模块已跑、但声明产物在磁盘上找不到的条目 (E2/G3 复发护栏)。

    历史坑: report_builder 声明 `geo_analysis/spatiotemporal_counts.csv`, 产出端
    若写在别的目录, 报告只会「静默少一张表」。此处把「声明了但没找到」显式列出。

    2026-10-16 修 (#2/#3): 门控由「本轮 stages_done 命中」扩为「或 .checkpoints/<stage>.done
    在盘上」。原先 stages_done 为空 (report-only 重跑等) 直接整段静默, 声明漂移
    查不出来; 而模块整个没跑 (既无 stages_done 也无 .done) 时仍不告警 —— 不制造噪音。
    """
    if not vdir:
        return []
    stages_done = stages_done or []
    if not (any(s in stages_done for s in m["stages"])
            or any(_stage_done_on_disk(vdir, s) for s in m["stages"])):
        return []
    missing = []
    for sp in m["summaries"]:
        if any(ch in sp for ch in "*?["):
            hit = list(Path(vdir).glob(f"{m['dir']}/{sp}"))
        else:
            hit = [p for p in [os.path.join(vdir, m["dir"], sp)] if os.path.exists(p)]
        if not hit:
            missing.append(f"{m['dir']}/{sp}")
    return missing


def _fmt_cell(v) -> str:
    s = str(v)
    # 科学计数/长小数紧凑化
    if len(s) > 12:
        try:
            f = float(s)
            s = f"{f:.4g}"
        except ValueError:
            pass
    # 2026-10-12 修 (P2-7): 表格单元格一律 HTML 转义 —— 产物里的序列名/病毒名/
    # 备注可能含 `<`/`&`(如 `<script>`), 旧实现直接内插进 <td>, 报告页可被注入
    # 裸标签。所有单元格只经本函数这一个出口 (csv/json/kv 同源), 故只在这里
    # 转义一次 (render_md 自己转义, 不走本函数, 不会二次转义)。
    return _html.escape(s)


def render_table(path: str, max_rows: int = 50) -> str:
    """csv/tsv → HTML 表 (最多 max_rows 行)"""
    if not path or not os.path.exists(path):
        return ""
    ext = Path(path).suffix.lower()
    try:
        if ext in (".tsv", ".txt"):
            rows = list(csv.reader(open(path, encoding="utf-8", errors="replace"), delimiter="\t"))
        else:
            rows = list(csv.reader(open(path, encoding="utf-8", errors="replace")))
    except Exception:
        return ""
    if not rows:
        return ""
    # 去全空行
    rows = [r for r in rows if any(c.strip() for c in r)]
    header, body = rows[0], rows[1:1 + max_rows]
    thead = "".join(f"<th>{_fmt_cell(c)}</th>" for c in header)
    trs = []
    for r in body:
        trs.append("<tr>" + "".join(f"<td>{_fmt_cell(c)}</td>" for c in r) + "</tr>")
    extra = f'<p class="tbl-more">… {len(rows)-1-max_rows} more rows (see source file)</p>' if len(rows) - 1 > max_rows else ""
    return (f'<table class="mod-table"><thead><tr>{thead}</tr></thead>'
            f'<tbody>{"".join(trs)}</tbody></table>{extra}')


def render_json(path: str, max_rows: int = 40) -> str:
    """json → key/value 两列表 (嵌套 dict 展开)"""
    if not path or not os.path.exists(path):
        return ""
    try:
        data = json.load(open(path, encoding="utf-8", errors="replace"))
    except Exception:
        return ""
    if isinstance(data, dict):
        items = list(data.items())[:max_rows]
        trs = "".join(f'<tr><td class="k">{_fmt_cell(k)}</td><td>{_fmt_cell(v)}</td></tr>'
                      for k, v in items)
        return f'<table class="mod-table kv"><thead><tr><th>Parameter</th><th>Value</th></tr></thead><tbody>{trs}</tbody></table>'
    if isinstance(data, list):
        return render_kv_list(data, max_rows)
    return ""


def render_kv_list(items: List, max_rows: int = 30) -> str:
    rows = []
    for it in items[:max_rows]:
        if isinstance(it, dict):
            rows.append("<tr>" + "".join(f"<td>{_fmt_cell(v)}</td>" for v in it.values()) + "</tr>")
    if not rows:
        return ""
    cols = list(items[0].keys())
    thead = "".join(f"<th>{_fmt_cell(c)}</th>" for c in cols)
    return f'<table class="mod-table"><thead><tr>{thead}</tr></thead><tbody>{"".join(rows)}</tbody></table>'


def render_md(path: str, max_chars: int = 3000) -> str:
    """md/txt → <pre> 简单展示"""
    if not path or not os.path.exists(path):
        return ""
    try:
        text = open(path, encoding="utf-8", errors="replace").read(max_chars)
    except Exception:
        return ""
    text = text.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
    return f'<pre class="md-pre">{text}</pre>'


# 2026-10-16 修 (#4): clock 合并了旧 rtt/temporal/genes 三个子步骤
# (phylo_pipeline._LEGACY_CLOCK_SUBSTAGES)。徽章一律按新 stage 名 (clock) 显示,
# 但 .done 回退判定要兼容老病毒只写的 rtt.done / temporal.done / genes.done ——
# 只认 clock.done 会把已完成的老病毒显示成灰色未跑。
_LEGACY_STAGE_ALIASES = {"clock": ("rtt", "temporal", "genes")}

# 2026-10-17 改 (A): 侧栏状态体系 —— dot 颜色与 .badge ok/fail/skip/off 同色,
# 左边条 class 走 ③ 的 sb-pass/sb-fail/sb-skip (外加 sb-off 表示"未跑")。
_SB_COLORS = {"pass": "#27ae60", "fail": "#e74c3c", "skip": "#f39c12", "off": "#95a5a6"}
_SB_PRIORITY = {"off": 0, "skip": 1, "pass": 2, "fail": 3}


def module_state(m: Dict, r: Dict, has_dir: bool, n_rendered: int) -> str:
    """模块在该病毒下的聚合状态 (侧栏 dot 用): fail > pass > skip > off。

    口径直接复用 stage_badge 的判定 (同一份徽章 HTML 里读 class), 不另写一套
    完成度逻辑 —— 徽章显示 ok 侧栏就必须是绿点, 两处不许打架。
    无 stage 的模块 (report: 不写 checkpoint) 有产物目录/渲染物即算 pass。
    """
    state = "off"
    for s in m["stages"]:
        b = stage_badge(s, r)
        if 'class="badge fail"' in b:
            state = _merge_state(state, "fail")
        elif 'class="badge ok"' in b:
            state = _merge_state(state, "pass")
        elif 'class="badge skip"' in b:
            state = _merge_state(state, "skip")
    if not m["stages"] and (has_dir or n_rendered):
        state = _merge_state(state, "pass")
    return state


def _merge_state(a: str, b: str) -> str:
    """两状态取优先级更高者 (多病毒聚合 / 模块内多 stage 聚合共用)。"""
    return a if _SB_PRIORITY.get(a, 0) >= _SB_PRIORITY.get(b, 0) else b


def stage_badge(stage: str, r: Dict) -> str:
    vdir = r.get("virus_dir", "")
    done = stage in r.get("stages_done", [])
    skipped = stage in r.get("stages_skipped", [])
    failed = stage in r.get("stages_failed", [])
    # 2026-10-03: 受阻 stage (前置缺失/未产出/降级) 用 fail 色但文案区分;
    # 跳过原因 (checkpoint/na/pending) 显示在 badge 里, 不再一律 "skip"。
    blocked = stage in r.get("stages_blocked", []) or any(
        b.startswith(f"{stage}(") for b in r.get("stages_blocked", []))
    reason = (r.get("stages_skip_reasons") or {}).get(stage, "")
    # 本轮未涉及该 stage 时回退看 checkpoint (report 单独重跑场景);
    # 老 stage 名 (clock ← rtt/temporal/genes) 的 .done 标记同样认。
    if not (done or skipped or failed) and vdir:
        for _mark in (stage,) + _LEGACY_STAGE_ALIASES.get(stage, ()):
            if os.path.exists(os.path.join(vdir, ".checkpoints", f"{_mark}.done")):
                done = True
                break
    if done:
        return f'<span class="badge ok">{stage}</span>'
    if blocked:
        return f'<span class="badge fail">{stage} · blocked</span>'
    if skipped:
        return f'<span class="badge skip">{stage} · {reason or "skip"}</span>'
    if failed:
        return f'<span class="badge fail">{stage} · fail</span>'
    return f'<span class="badge off">{stage}</span>'


# ─────────────────────────────────────────────────────────────
# 主函数
# ─────────────────────────────────────────────────────────────

def build_html_report(
    results: List[Dict],
    output_dir: str,
    title: str = "PhyloPipeline Report",
) -> str:
    """
    从分析结果列表生成模块化 HTML 汇总报告。

    results 每个元素包含 virus, phase, success, outputs, metrics,
    stages_done, stages_skipped, stages_failed, 以及病毒 work_dir (virus_dir)。
    """
    summary_rows = []
    virus_sections = []
    # 2026-10-16 (#8) / 2026-10-17 改 (A): 侧栏数据 —— virus_nav: (index, 已转义病毒名,
    # 总状态); mod_state: 模块 key → 跨病毒聚合状态 (多病毒取优先级最高者)。
    virus_nav = []
    mod_state: Dict[str, str] = {}

    for i, r in enumerate(results):
        # 2026-10-12 修 (P2-7): 病毒名进 HTML (总表 <td> 与卡片 <h3>) 前必须转义 ——
        # 病毒名来自目录/数据集元数据, 可能含 `<`/`&`。
        vname = _html.escape(str(r.get("virus", f"Virus_{i}")))
        metrics = r.get("metrics", {})
        success = r.get("success", False)
        vdir = r.get("virus_dir") or r.get("output_dir") or ""

        # ── 总表行: 加 Mantel / TMRCA / stage 计数 ──
        # 2026-09-15 (P0-A): 头条指标统一走 headline_metrics, 与 phylo_summary.csv
        # 完全同源; 每个值后面带来源标签, 读者不必猜是 treedater / PIC / RTT。
        hm = headline_metrics(metrics)
        mantel = metrics.get("mantel_r", "")
        tmrca_cell = fmt_metric(hm["tmrca"], hm["tmrca_source"])
        n_done = len(r.get("stages_done", []))
        n_skip = len(r.get("stages_skipped", []))
        n_fail = len(r.get("stages_failed", []))
        n_blocked = len(r.get("stages_blocked", []))
        # report 单独重跑: 本轮 stages 空时以 checkpoint 为准
        if not (n_done or n_skip or n_fail) and vdir:
            cp_dir = os.path.join(vdir, ".checkpoints")
            if os.path.isdir(cp_dir):
                n_done = len([f for f in os.listdir(cp_dir) if f.endswith(".done")])
        # 2026-10-03 (P1-1): 受阻/未定 (blocked|pending) 一律降为 PARTIAL, 不再报 PASS
        if n_fail:
            overall = "FAIL"
        elif n_blocked or r.get("incomplete"):
            overall = "PARTIAL"
        elif n_done > 0 or success:
            overall = "PASS"
        else:
            overall = "PARTIAL"

        # ── 每模块一节 ──
        mod_html = []
        n_rendered = 0        # 实际渲染出的表/图件数 (P2-8 闸门用)
        n_decl_missing = 0    # "声明但未找到"的产物数
        for m in MODULES:
            mdir = os.path.join(vdir, m["dir"]) if vdir else ""
            badges = "".join(stage_badge(s, r) for s in m["stages"])
            has_dir = mdir and os.path.isdir(mdir)

            # summary 内容
            tables = []
            for sp in m["summaries"]:
                matched = []
                if any(ch in sp for ch in "*?["):
                    if vdir:
                        matched = [str(p) for p in sorted(Path(vdir).glob(f"{m['dir']}/{sp}"))][:5]
                else:
                    sp_full = os.path.join(vdir, m["dir"], sp) if vdir else ""
                    if sp_full and os.path.exists(sp_full):
                        matched = [sp_full]
                for sp_full in matched:
                    ext = Path(sp_full).suffix.lower()
                    label = os.path.relpath(sp_full, os.path.join(vdir, m["dir"])) if vdir else Path(sp).name
                    if ext == ".json":
                        t = render_json(sp_full)
                    elif ext == ".md":
                        t = render_md(sp_full)
                    else:
                        t = render_table(sp_full)
                    if t:
                        tables.append(f'<div class="tbl-wrap"><div class="tbl-title">{label}</div>{t}</div>')

            # plots (E1: png/svg/jpg/gif/webp/pdf; 同 stem 取高优先级; 被排除的显式列出)
            plot_skips = []
            # 2026-10-16 修 (#7): 逐模块内联上限 —— MODULES 条目的 max_imgs 覆盖默认值
            # (phylogeny/select/report = 60), 超出部分照旧进 plot_skips 并带链接列出。
            imgs = collect_plots(mdir, m["plots"],
                                 max_imgs=m.get("max_imgs", MAX_IMGS_PER_MODULE),
                                 skipped=plot_skips) if has_dir else []
            img_tags, not_inlined = [], []
            for p in imgs:
                card, reason = asset_card(p, output_dir)
                img_tags.append(card)
                if reason:
                    not_inlined.append((p, reason))

            # extra_dirs (2026-10-12): 模块目录之外的产物树 (EasyHap_results/、
            # haplo_outputs/… 按 vdir 相对路径收集; 图走 asset_card, 表走
            # render_table/render_json; 缺目录静默跳过)
            for xd in m.get("extra_dirs", []):
                xdir = os.path.join(vdir, xd["dir"]) if vdir else xd["dir"]
                if not os.path.isdir(xdir):
                    continue
                import glob as _glob
                for sp in xd.get("tables", []):
                    # 支持 glob (文件名常带区域前缀, 如 GCVA_1_16310.HapSummary.tsv)
                    _hits = sorted(_glob.glob(os.path.join(xdir, sp))) if any(
                        ch in sp for ch in "*?[") else [
                        p for p in [os.path.join(xdir, sp)] if os.path.isfile(p)]
                    for _p in _hits:
                        _ext = Path(_p).suffix.lower()
                        _t = render_json(_p) if _ext == ".json" else render_table(_p)
                        if _t:
                            _label = os.path.relpath(_p, vdir) if vdir else _p
                            tables.append(f'<div class="tbl-wrap"><div class="tbl-title">'
                                          f'{_label}</div>{_t}</div>')
                _xsk = []
                _xi = collect_plots(xdir, xd["plots"],
                                    max_imgs=xd.get("max_imgs", 20), skipped=_xsk)
                for p in _xi:
                    card, reason = asset_card(p, output_dir)
                    img_tags.append(card)
                    if reason:
                        not_inlined.append((p, reason))

            # E2: stage 已完成但声明的 summary 产物不在磁盘 → 显式列出 (G3 复发护栏)
            missing = missing_declared_summaries(vdir, m, r.get("stages_done", []))
            n_decl_missing += len(missing)
            n_rendered += len(tables) + len(img_tags)
            warn_html = ""
            if missing:
                warn_html += ('<p class="asset-note">⚠ declared but not found (output path may '
                              'not match the report glob): ' + _html.escape(", ".join(missing)) + '</p>')
            if not_inlined:
                warn_html += ('<p class="asset-note">⚠ assets not inlined: ' + _html.escape(
                    "; ".join(f"{Path(p).name} ({why})" for p, why in not_inlined)) + '</p>')
            if plot_skips:
                # 2026-10-16 修 (#7): 被排除候选原先只列前 8 个文件名且无链接 ——
                # 读者看到"被排除"却点不开产物。现改为带相对路径下载链接 (与
                # asset_card 的 Download 写法一致), 上限 MAX_SKIP_LIST 条, 余量只报数量。
                _sk_items = "; ".join(
                    f'<a href="{_html.escape(_rel_href(p, output_dir), quote=True)}" download>'
                    f'{_html.escape(Path(p).name, quote=True)}</a> ({_html.escape(why, quote=True)})'
                    for p, why in plot_skips[:MAX_SKIP_LIST])
                _sk_more = (f'; … {len(plot_skips) - MAX_SKIP_LIST} more not listed'
                            if len(plot_skips) > MAX_SKIP_LIST else "")
                warn_html += ('<p class="asset-note">· candidates not embedded in this module ('
                              f'{len(plot_skips)}): {_sk_items}{_sk_more}</p>')

            # 未运行且无目录 → 整节灰显
            state_cls = "mod" if (has_dir or tables or img_tags) else "mod off"
            # 2026-10-17 改 (A): 侧栏 dot 的模块状态 (跨病毒聚合)
            mod_state[m["key"]] = _merge_state(
                mod_state.get(m["key"], "off"),
                module_state(m, r, has_dir, len(tables) + len(img_tags)))
            # 2026-10-16 新增 (#8): 模块锚点供顶栏导航跳转。HTML 要求 id 唯一,
            # 首个病毒用 mod-{key} (顶栏链接目标), 其余病毒加 -{i} 后缀。
            _mod_id = f"mod-{m['key']}" if i == 0 else f"mod-{m['key']}-{i}"
            notes = "".join(f'<li>{n}</li>' for n in m.get("notes", []))
            body = ""
            if tables:
                body += "".join(tables)
            if img_tags:
                body += f'<div class="img-grid">{"".join(img_tags)}</div>'
            if warn_html:
                body += warn_html
            if not body:
                body = f'<p class="mod-empty">{"not run / not applicable (skip)" if not has_dir else "no table or figure assets in this module"}</p>'
            mod_html.append(f"""
            <div class="{state_cls}" id="{_mod_id}">
                <div class="mod-head">
                    <h4>{m['name']}</h4>
                    <div class="badges">{badges}</div>
                </div>
                {f'<ul class="mod-notes">{notes}</ul>' if notes and not tables else ''}
                {body}
            </div>""")

        # 2026-10-12 修 (P2-8): 「零产物仍 PASS」闸门 —— overall 原先只看
        # stages_done/success: stage 全部"完成"但**一张表/一张图都没渲染出来**
        # (声明的 summary 产物全部未找到) 时, 报告照打 PASS。现改为:
        # 渲染产物计数为 0 且存在"声明但未找到"条目 → 降为 PARTIAL 并给出原因。
        overall_note = ""
        if n_rendered == 0 and n_decl_missing and overall == "PASS":
            overall = "PARTIAL"
            overall_note = (f'<span class="asset-note">PARTIAL reason: zero rendered assets — '
                            f'{n_decl_missing} declared asset(s) not found</span>')

        summary_rows.append(f"""
        <tr>
            <td>{vname}</td>
            <td><span class="badge {'ok' if overall == 'PASS' else ('fail' if overall == 'FAIL' else 'skip')}">{overall}</span>{overall_note}</td>
            <td>{n_done} / {n_skip} / {n_fail} / {n_blocked}</td>
            <td>{fmt_metric(hm['rate'], hm['rate_source'])}</td>
            <td>{fmt_metric(hm['r2'], hm['r2_source'])}</td>
            <td>{tmrca_cell}</td>
            <td>{mantel if mantel != '' else '—'}</td>
            <td>{metrics.get("n_pathways", '') if metrics.get("n_pathways", '') != '' else '—'}</td>
            <td>{metrics.get("geo_mean_velocity_km_yr", '') if metrics.get("geo_mean_velocity_km_yr", '') != '' else '—'}</td>
        </tr>""")

        virus_sections.append(f"""
        <div class="virus-card" id="virus-{i}">
            <h3>{vname} <span class="badge {'ok' if overall == 'PASS' else ('fail' if overall == 'FAIL' else 'skip')}">{overall}</span>{overall_note}</h3>
            <div class="metrics">
                <span>β={fmt_metric(hm['rate'], hm['rate_source'])}</span>
                <span>R²={fmt_metric(hm['r2'], hm['r2_source'])}</span>
                <span>P={fmt_metric(hm['p'], hm['p_source'])}</span>
                <span>TMRCA={fmt_metric(hm['tmrca'], hm['tmrca_source'])}</span>
            </div>
            {''.join(mod_html)}
        </div>""")
        virus_nav.append((i, vname, overall))   # 侧栏病毒跳转清单 (含总状态)

    # ── 常驻左侧边栏 (2026-10-17 改 A: 取代顶部 sticky chips) ──
    # 视觉与 ③ virome_discovery_pipeline/report_pipeline.py 的 sb 体系同款;
    # 结构: "Pipeline Modules" (8 模块, 带各自状态 dot) + "Viruses" (各病毒锚点)。
    # 侧栏常驻 (见 CSS 注释), scroll-spy 高亮当前所在 section。
    _sb_mods = "".join(
        f'<a href="#mod-{m["key"]}" class="sb-item sb-{mod_state.get(m["key"], "off")}"'
        f' title="{_html.escape(m["name"], quote=True)}">'
        f'<span class="sb-dot" style="background:{_SB_COLORS[mod_state.get(m["key"], "off")]}"></span>'
        f'<span class="sb-label">{_html.escape(m["key"].title())}</span></a>'
        for m in MODULES)
    _sb_viruses = ""
    for _i, _vn, _ov in virus_nav:
        _cls = {"PASS": "pass", "FAIL": "fail"}.get(_ov, "skip")   # PARTIAL → skip 色
        _sb_viruses += (
            f'<a href="#virus-{_i}" class="sb-item sb-{_cls}" title="{_vn}">'
            f'<span class="sb-dot" style="background:{_SB_COLORS[_cls]}"></span>'
            f'<span class="sb-label">{_vn}</span></a>')
    sidebar = f"""
<div class="sidebar">
  <div class="sb-title">Pipeline Modules</div>
  {_sb_mods}
  <div class="sb-title" style="margin-top:10px;">Viruses</div>
  {_sb_viruses}
</div>"""

    # Full HTML
    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<!-- impeccable-disable numbered-section-markers,em-dash-overuse: module order & CJK scientific copy -->
<title>{title}</title>
<style>
    * {{ margin: 0; padding: 0; box-sizing: border-box; }}
    body {{ font-family: -apple-system, 'Segoe UI', Arial, sans-serif; background: #f5f6fa; color: #2c3e50; padding: 20px 20px 20px 0; }}
    /* 2026-10-17 impeccable: 标题/hero 用衬线, 正文保持系统 sans —— 双字体层级,
       零外部资源 (EB Garamond 缺失时回退 Georgia / Times)。 */
    h1, h2, h3, h4 {{ font-family: 'EB Garamond', Georgia, 'Times New Roman', serif; font-weight: 600; letter-spacing: .01em; }}
    h1 {{ text-align: center; margin: 20px 0; font-size: 28px; color: #1a1a2e; }}
    .summary-table {{ width: 100%; max-width: 1100px; margin: 20px auto; border-collapse: collapse; background: white; border-radius: 8px; overflow: hidden; box-shadow: 0 2px 8px rgba(0,0,0,0.08); }}
    .summary-table th {{ background: #1a1a2e; color: white; padding: 12px 16px; text-align: left; font-size: 13px; }}
    .summary-table td {{ padding: 10px 16px; border-bottom: 1px solid #eee; font-size: 13px; }}
    .summary-table tr:hover {{ background: #f8f9ff; }}
    .virus-card {{ background: white; border-radius: 10px; padding: 24px; margin: 20px auto; max-width: 1200px; box-shadow: 0 2px 12px rgba(0,0,0,0.06); }}
    .virus-card h3 {{ font-size: 18px; margin-bottom: 12px; display: flex; align-items: center; gap: 10px; }}
    .badge {{ font-size: 11px; padding: 3px 10px; border-radius: 12px; font-weight: 600; white-space: nowrap; }}
    .badge.ok {{ background: #27ae6018; color: #27ae60; }}
    .badge.fail {{ background: #e74c3c18; color: #e74c3c; }}
    .badge.skip {{ background: #f39c1218; color: #b9770e; }}
    .badge.off {{ background: #95a5a622; color: #7f8c8d; }}
    .metrics {{ display: flex; gap: 16px; margin-bottom: 16px; font-size: 13px; color: #666; }}
    .metrics span {{ background: #f0f0f5; padding: 4px 10px; border-radius: 6px; font-family: monospace; }}
    .mod {{ background: #fafafc; border: 1px solid #e8e8ef; border-radius: 8px; padding: 16px; margin: 14px 0; }}
    .mod.off {{ opacity: 0.55; }}
    .mod-head {{ display: flex; justify-content: space-between; align-items: center; flex-wrap: wrap; gap: 8px; margin-bottom: 8px; }}
    .mod-head h4 {{ font-size: 15px; color: #1a1a2e; }}
    .badges {{ display: flex; gap: 6px; flex-wrap: wrap; }}
    .mod-notes {{ margin: 4px 0 8px 18px; font-size: 12px; color: #888; }}
    .mod-empty {{ font-size: 13px; color: #aaa; font-style: italic; padding: 6px 0; }}
    .tbl-wrap {{ margin: 10px 0; }}
    .tbl-title {{ font-size: 12px; font-weight: 600; color: #555; margin-bottom: 4px; font-family: monospace; }}
    .mod-table {{ border-collapse: collapse; width: 100%; max-height: 400px; display: block; overflow: auto; background: white; }}
    .mod-table th {{ background: #2c3e50; color: white; padding: 6px 10px; font-size: 12px; text-align: left; position: sticky; top: 0; }}
    .mod-table td {{ padding: 5px 10px; border-bottom: 1px solid #eee; font-size: 12px; font-family: monospace; }}
    .mod-table.kv td.k {{ font-weight: 600; width: 40%; }}
    .tbl-more {{ font-size: 11px; color: #999; padding: 4px; }}
    .md-pre {{ background: white; border: 1px solid #eee; border-radius: 6px; padding: 12px; font-size: 12px; max-height: 300px; overflow: auto; white-space: pre-wrap; }}
    .img-grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(320px, 1fr)); gap: 16px; margin-top: 12px; }}
    .img-card {{ background: white; border-radius: 8px; overflow: hidden; border: 1px solid #eee; }}
    .img-card.off {{ border-style: dashed; }}
    .img-label {{ font-size: 12px; padding: 8px 12px; background: #1a1a2e; color: white; font-weight: 500; display: flex; justify-content: space-between; align-items: center; gap: 8px; }}
    .img-card img {{ width: 100%; display: block; }}
    .pdf-embed {{ width: 100%; height: 460px; display: block; background: #f0f0f5; }}
    .tag {{ font-size: 10px; padding: 2px 7px; border-radius: 9px; background: #ffffff2e; font-weight: 600; }}
    .tag.warn {{ background: #e74c3c; }}
    .asset-note {{ font-size: 11px; color: #8a5a00; background: #fdf6e3; border: 1px solid #f0d9a8; border-radius: 6px; padding: 6px 10px; margin-top: 8px; line-height: 1.6; }}
    .asset-note a {{ color: #2c3e50; }}
    .footer {{ text-align: center; margin-top: 40px; padding: 20px; color: #999; font-size: 12px; }}
    /* 2026-10-17 改 (A): 顶部 sticky chips 换成常驻左侧边栏 —— 视觉拷贝 ③
       (virome_discovery_pipeline/report_pipeline.py) 的 sb 体系:
       .sidebar/.sb-title/.sb-item/.sb-dot/.sb-label/.sb-active + 左边条状态色。
       ⚠ 用户裁定 (两次投诉侧栏消失): **侧栏常驻, 任何窗口宽度都不隐藏** ——
       不许在 max-width 断点里给侧栏写 display:none;
       只在 @media print 隐藏 (存档/打印不需要导航)。 */
    .sidebar {{ position: fixed; top: 0; left: 0; width: 224px; height: 100vh;
                background: #ffffff; border-right: 1px solid #e8e8ef;
                box-shadow: 2px 0 8px rgba(0,0,0,0.04); z-index: 200;
                display: flex; flex-direction: column; padding-top: 12px; overflow-y: auto; }}
    .sb-title {{ font-size: 11px; font-weight: 700; color: #7f8c8d; text-transform: uppercase;
                 letter-spacing: .8px; padding: 8px 16px 6px; margin: 0; }}
    .sb-item {{ display: flex; align-items: center; gap: 10px; padding: 9px 16px;
                text-decoration: none; color: #2c3e50; font-size: 12.5px; font-weight: 500;
                border: 0; transition: background .15s, color .15s; }}
    .sb-item:hover {{ background: rgba(83,125,150,.08); }}
    .sb-dot {{ width: 10px; height: 10px; border-radius: 50%; flex-shrink: 0; }}
    .sb-label {{ overflow: hidden; text-overflow: ellipsis; white-space: nowrap; }}
    /* 状态只由 sb-dot 圆点承载 (颜色内联自 _SB_COLORS); 不再用左侧彩条。 */
    .sb-active {{ background: rgba(83,125,150,.1); font-weight: 600; }}
    /* 主内容让出侧栏宽度 (body 左 padding 归零, 偏移全部由这里承担 → 内容起点 240px) */
    .container {{ padding-left: 240px; }}
    /* 非 sticky 布局: 锚点只需小位移, 不再为吸顶导航留 64px */
    .virus-card, .mod {{ scroll-margin-top: 16px; }}
    @media print {{ .sidebar {{ display: none; }} .container {{ padding-left: 20px; }} }}
</style>
</head>
<body>
{sidebar}
<div class="container">
<h1>{title}</h1>
<p style="text-align:center;color:#888;margin-bottom:30px">
    Generated: {datetime.now().strftime('%Y-%m-%d %H:%M')} · MMPV-RNA PhyloPipeline · VirPhyKit
</p>

<h2 style="max-width:1100px;margin:0 auto 10px;font-size:16px;">Summary</h2>
<table class="summary-table">
    <thead>
        <tr><th>Virus</th><th>Status</th><th>done/skip/fail/blocked</th><th>β (rate)</th><th>R²</th><th>TMRCA</th><th>Mantel r</th><th>Geo paths</th><th>Geo velocity</th></tr>
    </thead>
    <tbody>
        {''.join(summary_rows)}
    </tbody>
</table>
<p style="max-width:1100px;margin:8px auto 0;color:#888;font-size:12px;">
    Metric sources are given in parentheses. Priority: rate td_mean_rate(treedater) → pic_rate(PIC);
    R²  pic_r(PIC) → r_squared(RTT); P  td_relaxed_clock_p(treedater) → p_value(RTT);
    TMRCA  beast_tmrca(BEAST) → td_tmrca(treedater).
    Same source and same values as <code>phylo_summary.csv</code> (that table also stores <code>*_source</code> and <code>alt_*</code> alternate-estimator columns).
</p>

<h2 style="max-width:1200px;margin:30px auto 10px;font-size:16px;">Per-virus Results by Module</h2>
{''.join(virus_sections) if virus_sections else '<p style="text-align:center;color:#999;">No results to display</p>'}

<div class="footer">MMPV-RNA PhyloPipeline · Powered by VirPhyKit (Yin &amp; Gao, 2025)</div>
</div><!-- /.container -->
<script>
// 侧栏 scroll-spy (2026-10-17 改 A; 从 ③ report_pipeline.py 的同款脚本拷贝,
// 目标选择器换成本报告的锚点: .virus-card[id] / .mod[id])。纯内联, 零外部资源。
(function(){{
  var items=document.querySelectorAll('.sb-item');
  var targets=document.querySelectorAll('.virus-card[id], .mod[id]');
  if(!items.length||!targets.length)return;
  function onScroll(){{
    var current='';
    targets.forEach(function(t){{ if(t.getBoundingClientRect().top<=160) current=t.id; }});
    items.forEach(function(a){{
      a.classList.toggle('sb-active', a.getAttribute('href')==='#'+current);
    }});
  }}
  window.addEventListener('scroll',onScroll,{{passive:true}});
  onScroll();
}})();
</script>
{_TW}</body>
</html>"""

    report_path = os.path.join(output_dir, "phylo_report.html")
    with open(report_path, 'w', encoding='utf-8') as f:
        f.write(html)

    return report_path
