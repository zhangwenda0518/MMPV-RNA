#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""host_expression.py — ⑧ 宿主表达管线 (可选补充模块) 编排器
=========================================================================
从公共转录组的宿主侧做"表达定量 → 差异表达"分析。抽取自 goji-virome
差异表达分析项目 (枸杞病毒宿主响应) 的核心链路, 以**可选补充模块**形式
并入 MMPV: 不在 ①→⑤ 主链上、不改变任何既有管线行为; 输入全部取自主链
现有产物, 输出落独立根 08_HostExpression/。

阶段 (内部编号 = 运行顺序):
  01_HostRef  宿主 CDS 提取 (genome + GFF3 → gffread -x; 或 --reference 直接给)
  02_Quant    宿主转录组 salmon 定量 (salmon_auto_pipeline.py, FASTA/FASTQ 均可)
  03_QC       定量 QC (比对率/reads 数双阈值, 剔除 16S 垃圾样本)
  04_Labels   病毒标签推导 (④ 检出表 → 阳性/阴性 + 连续载量; 合并 Tissue 元数据)
  05_Coexp    TEA-GCN 差异共表达 (主分析): 阳/阴双网络 → ΔzScore(MR) 显著对
              (teagcn/ 引擎; 闸门: 每组>=30 样本且失衡<=3×; mp 入口仅 Linux)
  06_DE       DESeq2 差异表达 (交叉参考; 可选 --de-strata Tissue 组织分层)
  07_Enrich   功能富集 ORA (紧跟 06_DE: 消费 DE/hub/模块基因集 +
              用户注释 eggNOG-mapper/通用 gene-term 映射; 缺注释时跳过)
  08_LoadScan 病毒载量×全基因组 Spearman 扫描 (可选剂量验证翼, 默认跑,
              --no-with-loadscan 可关; Tissue 存在时自动组内扫描)
  09_Report   HTML 汇总报告

输入边界 (可用 --project-root 自动定位, 也可显式给路径):
  ② 02_Preprocessing/01_CleanData/{3.clumpify,2.fasta}  → --reads_dir
  ① 01_PublicData/05_HostRef/genome/                    → --genome / --gff
  ④ 04_Analysis/01_Detection/summary/all_viruses.best.summary.tsv → --detect
  ① 01_PublicData/02_Unified/Global_Unified_Metadata_Core14.{csv,tsv} → --metadata

用法:
  # 全流程 (输入自动定位; 参考由 genome+GFF3 现场提取)
  python host_expression_pipeline/host_expression.py --stage all \\
      --project-root <项目根> -t 60 -j 8

  # 已有宿主 CDS 时跳过提取
  python host_expression_pipeline/host_expression.py --stage all \\
      --output_dir <项目>/08_HostExpression --reads_dir clean/ \\
      --reference host_cds.fasta --detect all_viruses.best.summary.tsv -t 60

  # 单阶段重跑 / 预览
  python host_expression_pipeline/host_expression.py --stage de --output_dir ... --force
  python host_expression_pipeline/host_expression.py --stage all --dry-run ...

依赖: salmon (pixi 已含), gffread (仅 hostref 阶段), pydeseq2/组学库 + einsumt
      (coexp/de/loadscan 阶段; pixi [pypi-dependencies] 已带, 见模块 README)
"""

import argparse
import logging
import os
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Optional, Sequence

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

from mmpv_common.io_layout import (CLEAN_READS_SUBDIRS, PIPELINE_ROOT_NAMES,
                                   build_hostexpr_dirs, dir_name,
                                   locate_detection_summary,
                                   normalize_layout_env, write_boundary)

SCRIPT_DIR = Path(__file__).resolve().parent
TEAGCN_DIR = SCRIPT_DIR / "teagcn"
# 阶段顺序 = 编号顺序; 05_Coexp (TEA-GCN Δz) 为主分析, 06/07 为交叉参考/剂量验证,
# 07_Enrich 功能富集 (消费 05/06 基因集 + 用户注释), 09_Report 收尾
STAGES = ("hostref", "quant", "qc", "labels", "coexp", "de", "enrich",
          "loadscan", "report")

# 各阶段的产物 (断点续传判定: .ok 标记 + 产物存在 双条件)
STAGE_PRODUCTS = {
    "hostref":  (),          # 产物即参考序列本身, 单独判定
    "quant":    ("matrix", "NumReads_matrix.tsv"),
    "qc":       ("qc_summary.tsv",),
    "labels":   ("sample_labels.tsv",),
    "coexp":    ("diff_coexp_significant.tsv", "diff_coexp_summary.tsv"),
    "de":       ("deseq2_all_results.tsv",),
    "loadscan": ("load_scan_summary.tsv",),
    "enrich":   ("enrich_summary.tsv",),
    "report":   ("host_expression_report.html",),
}
# coexp 内部子步 (昂贵: 双网络构建数小时; 逐子步断点)
#   qc = 网络 z 分布闸门 (KingdomStress Fig S8); modules = DEG 限定模块检测 (Fig 5)
COEXP_SUBSTEPS = ("matrices", "partitions", "networks", "ranks", "qc", "diff",
                  "modules")


def setup_logger(output_dir: Path, log_file: bool = True) -> logging.Logger:
    logger = logging.getLogger("HostExpression")
    logger.setLevel(logging.DEBUG)
    logger.handlers.clear()
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    ch.setFormatter(logging.Formatter("[%(asctime)s] %(levelname)s %(message)s",
                                      datefmt="%H:%M:%S"))
    logger.addHandler(ch)
    if log_file:                       # dry-run 时不落盘、不建目录
        output_dir.mkdir(parents=True, exist_ok=True)
        fh = logging.FileHandler(output_dir / "host_expression.log",
                                 encoding="utf-8")
        fh.setLevel(logging.DEBUG)
        fh.setFormatter(logging.Formatter("[%(asctime)s] %(levelname)s %(message)s"))
        logger.addHandler(fh)
    return logger


def run_cmd(cmd: Sequence[str], log: logging.Logger, step: str,
            cwd: Optional[Path] = None) -> None:
    log.info("[%s] %s", step, " ".join(str(c) for c in cmd)
             + (f"  (cwd={cwd})" if cwd else ""))
    subprocess.run([str(c) for c in cmd], check=True, cwd=str(cwd) if cwd else None)


def _dry_or_raise(args, log: logging.Logger, label: str, func, placeholder: str):
    """dry-run 下输入定位失败不中断 (用占位符展示命令); 实跑时原样抛错。"""
    try:
        return func()
    except FileNotFoundError as exc:
        if args.dry_run:
            log.warning("[dry-run] %s 未定位 (%s)", label, exc)
            return Path(placeholder)
        raise


def _resolve_located(explicit: Optional[str], locator, project_root) -> Optional[Path]:
    """显式路径优先; 否则经 locator 从项目根定位 (统一 resolve 为绝对路径)。"""
    if explicit:
        return Path(explicit).resolve()
    if not project_root:
        return None
    hit = locator(Path(project_root))
    return hit.resolve() if hit else None


# ═══════════════════════════════════════════════════════════════════
# 输入自动定位 (standard 优先, legacy 兜底; 全部经 io_layout 键名)
# ═══════════════════════════════════════════════════════════════════

def locate_clean_reads(project_root: Path) -> Optional[Path]:
    """② clean reads 最终子目录 (3.clumpify 优先, 回退 2.fasta)。"""
    roots = [project_root / PIPELINE_ROOT_NAMES["preprocessing"],
             project_root / PIPELINE_ROOT_NAMES["discovery"],
             project_root]
    for base_root in roots:
        for layout in ("standard", "legacy"):
            base = base_root / dir_name("d_clean", layout)
            for sub in CLEAN_READS_SUBDIRS:
                cand = base / sub
                if cand.is_dir() and any(cand.iterdir()):
                    return cand.resolve()
    return None


def locate_host_genome(project_root: Path) -> Optional[Path]:
    """① 宿主基因组 FASTA (all.genome.uniq.fasta 优先)。"""
    bases = [project_root / PIPELINE_ROOT_NAMES["publicdata"]
             / dir_name("h_genome", "standard"),
             project_root / dir_name("h_genome", "legacy")]
    for base in bases:
        if not base.is_dir():
            continue
        first = base / "all.genome.uniq.fasta"
        if first.is_file():
            return first
        hits = sorted(p for p in base.glob("*")
                      if p.suffix.lower() in (".fa", ".fasta", ".fna", ".gz"))
        if hits:
            return hits[0]
    return None


def locate_host_gff(project_root: Path) -> Optional[Path]:
    """① 宿主注释 GFF3。"""
    bases = [project_root / PIPELINE_ROOT_NAMES["publicdata"]
             / dir_name("h_genome", "standard"),
             project_root / dir_name("h_genome", "legacy")]
    for base in bases:
        if not base.is_dir():
            continue
        hits = sorted(p for p in base.iterdir()
                      if p.name.lower().endswith((".gff3", ".gff", ".gff3.gz")))
        if hits:
            return hits[0]
    return None


def locate_metadata_core14(project_root: Path) -> Optional[Path]:
    """① Core14 元数据 (Global_Unified_Metadata_Core14.{csv,tsv})。"""
    bases = [project_root / PIPELINE_ROOT_NAMES["publicdata"]
             / dir_name("m_info", "standard"),
             project_root / dir_name("m_info", "legacy")]
    for base in bases:
        if not base.is_dir():
            continue
        for name in ("Global_Unified_Metadata_Core14.csv",
                     "Global_Unified_Metadata_Core14.tsv"):
            f = base / name
            if f.is_file():
                return f
        hits = sorted(base.glob("Global_Unified_Metadata_Core14*"))
        if hits:
            return hits[0]
    return None


def locate_detect_table(project_root: Path) -> Optional[Path]:
    """④ 检出表 (04_Analysis/01_Detection/summary/all_viruses.best.summary.tsv)。"""
    roots = [project_root / PIPELINE_ROOT_NAMES["analysis"], project_root]
    for root in roots:
        hit = locate_detection_summary(root)
        if hit is not None:
            return hit.resolve()
    return None


# ═══════════════════════════════════════════════════════════════════
# 断点续传 (标记 + 产物核验)
# ═══════════════════════════════════════════════════════════════════

def stage_products(dirs: Dict[str, Path], stage: str,
                   reference: Optional[Path]) -> List[Path]:
    if stage == "hostref":
        return [reference] if reference else []
    base = dirs[stage]
    rel = STAGE_PRODUCTS[stage]
    if stage == "quant":
        return [base / "matrix" / "NumReads_matrix.tsv"]
    return [base / p for p in rel]


def stage_done(output_dir: Path, stage: str, products: Sequence[Path]) -> bool:
    if not products or not all(Path(p).is_file() for p in products):
        return False
    return (output_dir / ".checkpoints" / f"{stage}.ok").is_file()


def mark_done(output_dir: Path, stage: str) -> None:
    ck = output_dir / ".checkpoints"
    ck.mkdir(parents=True, exist_ok=True)
    (ck / f"{stage}.ok").write_text("ok\n", encoding="utf-8")


# ═══════════════════════════════════════════════════════════════════
# 阶段实现
# ═══════════════════════════════════════════════════════════════════

def stage_hostref(args, dirs, log) -> Optional[Path]:
    """宿主 CDS 提取: 用户参考优先, 否则 genome+GFF3 → gffread -x。"""
    if args.reference:
        ref = Path(args.reference).resolve()
        if not ref.is_file():
            raise FileNotFoundError(f"--reference 不存在: {ref}")
        log.info("使用用户提供的参考: %s", ref)
        return ref

    ref_out = dirs["hostref"] / "host_cds.fasta"
    if ref_out.is_file() and not args.force:
        log.info("参考已存在, 跳过提取: %s", ref_out)
        return ref_out

    genome = _resolve_located(args.genome, locate_host_genome, args.project_root)
    gff = _resolve_located(args.gff, locate_host_gff, args.project_root)
    if genome is None or gff is None:
        msg = ("hostref 需要 --reference (现成 CDS/转录本 FASTA), "
               "或 --genome + --gff (可用 --project-root 从 ① 宿主参考产物 "
               "05_HostRef/genome 自动定位)")
        if args.dry_run:
            log.warning("[dry-run] %s", msg)
            return ref_out
        raise FileNotFoundError(msg)

    dirs["hostref"].mkdir(parents=True, exist_ok=True)
    cmd = ["gffread", str(gff), "-g", str(genome), "-x", str(ref_out)]
    if args.dry_run:
        if shutil.which("gffread") is None:
            log.warning("[dry-run] gffread 不在 PATH (实际运行时需 pixi 环境)")
        log.info("[dry-run] %s", " ".join(cmd))
        return ref_out
    if shutil.which("gffread") is None:
        raise EnvironmentError("未找到 gffread (pixi 环境应已含); "
                               "或改用 --reference 直接提供 CDS FASTA")
    run_cmd(cmd, log, "hostref")
    if not ref_out.is_file() or ref_out.stat().st_size == 0:
        raise RuntimeError(f"gffread 未产出 CDS: {ref_out}")
    with ref_out.open(encoding="utf-8", errors="replace") as fh:
        n = sum(1 for l in fh if l.startswith(">"))
    log.info("宿主 CDS 提取完成: %s (%d 条序列)", ref_out, n)
    return ref_out


def resolve_reference(args, dirs, log) -> Path:
    if args.reference:
        return Path(args.reference).resolve()
    cand = dirs["hostref"] / "host_cds.fasta"
    if cand.is_file():
        return cand
    raise FileNotFoundError(
        "无可用参考序列: 请先跑 hostref 阶段 (或 --reference 指定 CDS FASTA)")


def resolve_reads(args) -> Path:
    if args.reads_dir:
        reads = Path(args.reads_dir).resolve()
        if not reads.is_dir():
            raise FileNotFoundError(f"--reads_dir 不存在: {reads}")
        return reads
    if args.project_root:
        reads = locate_clean_reads(Path(args.project_root))
        if reads is not None:
            return reads.resolve()
    raise FileNotFoundError(
        "未找到 clean reads: 请 --reads_dir 指定 (② 产物 "
        "02_Preprocessing/01_CleanData/{3.clumpify,2.fasta}), "
        "或用 --project-root 自动定位")


def resolve_detect(args, log) -> Path:
    if args.detect:
        p = Path(args.detect).resolve()
        if not p.exists():
            raise FileNotFoundError(f"--detect 不存在: {p}")
        return p
    if args.project_root:
        hit = locate_detect_table(Path(args.project_root))
        if hit is not None:
            return hit.resolve()
    raise FileNotFoundError(
        "未找到病毒检出表: 请 --detect 指定 (④ 产物 "
        "04_Analysis/01_Detection/summary/all_viruses.best.summary.tsv), "
        "或用 --project-root 自动定位")


def resolve_samples_list(dirs) -> Optional[Path]:
    """labels 阶段的样本名单: 优先 QC 汇总 (含全样本), 回退 quant 的 sample_info。"""
    qc_summary = dirs["qc"] / "qc_summary.tsv"
    if qc_summary.is_file():
        return qc_summary
    info = dirs["quant"] / "metadata" / "sample_info.tsv"
    if info.is_file():
        return info
    return None


def stage_quant(args, dirs, log) -> None:
    reads = _dry_or_raise(args, log, "clean reads",
                          lambda: resolve_reads(args), "<reads_dir>")
    ref = _dry_or_raise(args, log, "宿主参考",
                        lambda: resolve_reference(args, dirs, log), "<reference>")
    threads_per_sample = max(1, args.threads // max(1, args.jobs))
    cmd = [sys.executable, str(SCRIPT_DIR / "salmon_auto_pipeline.py"),
           "-i", str(reads), "-r", str(ref), "-o", str(dirs["quant"]),
           "--salmon", args.salmon,
           "-j", str(args.jobs), "--threads-per-sample", str(threads_per_sample)]
    if args.force:
        cmd.append("--force-quant")
    if args.dry_run:
        log.info("[dry-run] %s", " ".join(cmd))
        return
    run_cmd(cmd, log, "quant")


def stage_qc(args, dirs, log) -> None:
    cmd = [sys.executable, str(SCRIPT_DIR / "quant_qc.py"),
           "--quant_dir", str(dirs["quant"]), "--output_dir", str(dirs["qc"]),
           "--min-mapped", f"{args.min_mapped:g}",
           "--min-processed", f"{args.min_processed:g}"]
    if args.dry_run:
        log.info("[dry-run] %s", " ".join(cmd))
        return
    run_cmd(cmd, log, "qc")


def stage_labels(args, dirs, log) -> None:
    out = dirs["labels"] / "sample_labels.tsv"
    # 通用模式: 用户自备分组表 (任意二分组列; 无需 ④ 病毒检出)
    if args.labels_file:
        src = Path(args.labels_file)
        if args.dry_run:
            log.info("[dry-run] cp %s %s   (外部分组表)", src, out)
            return
        if not src.is_file():
            raise FileNotFoundError(f"--labels-file 不存在: {src.resolve()}")
        dirs["labels"].mkdir(parents=True, exist_ok=True)
        shutil.copyfile(src, out)
        log.info("外部分组表: %s → %s (通用模式; loadscan 无病毒载量将自动跳过)",
                 src.resolve(), out)
        return
    # 病毒组学模式: 从 ④ 检出表推导 (默认)
    detect = _dry_or_raise(args, log, "病毒检出表",
                           lambda: resolve_detect(args, log), "<detect_table>")
    cmd = [sys.executable, str(SCRIPT_DIR / "derive_virus_labels.py"),
           "--detect", str(detect), "--output_dir", str(dirs["labels"]),
           "--cov-threshold", str(args.cov_threshold)]
    samples = resolve_samples_list(dirs)
    if samples is not None:
        cmd += ["--samples", str(samples)]
    if args.metadata:
        cmd += ["--metadata", str(args.metadata)]
    elif args.project_root:
        meta = locate_metadata_core14(Path(args.project_root))
        if meta is not None:
            cmd += ["--metadata", str(meta.resolve())]
    if args.dry_run:
        log.info("[dry-run] %s", " ".join(cmd))
        return
    run_cmd(cmd, log, "labels")


# ═══════════════════════════════════════════════════════════════════
# coexp 子步 (TEA-GCN Δz 主分析; 05_Coexp)
# ═══════════════════════════════════════════════════════════════════

def _coexp_ck(output_dir: Path, name: str) -> Path:
    return output_dir / ".checkpoints" / f"coexp_{name}.ok"


def _coexp_substep(name: str, products: Sequence[Path], fn, args,
                   output_dir: Path, log) -> bool:
    """coexp 子步断点: .ok + 产物双条件; 返回是否实际执行。

    dry-run 下结构性不执行 (预览由 stage_coexp 的预览块负责)。"""
    if args.dry_run:
        return False
    if not args.force and _coexp_ck(output_dir, name).is_file() \
            and all(Path(p).exists() for p in products):
        log.info("[coexp/%s] 已完成, 跳过", name)
        return False
    fn()
    _coexp_ck(output_dir, name).parent.mkdir(parents=True, exist_ok=True)
    _coexp_ck(output_dir, name).write_text("ok\n", encoding="utf-8")
    return True


def _find_completed_gcn(work: Path) -> Path:
    """定位 TEA-GCN 工作目录下唯一的 Completed_GCN/<网络名>/。"""
    base = work / "Completed_GCN"
    hits = [p for p in base.iterdir() if p.is_dir()] if base.is_dir() else []
    if len(hits) != 1:
        raise RuntimeError(f"{base} 下应有且仅有一个网络目录, 实际 {len(hits)}")
    return hits[0]


def _coexp_group_spec(args, labels_path: Path):
    """coexp 分组口径: 默认 --group-col +/−; --virus 给定时切 <病毒>_positive 1/0。

    病毒列缺失时列出 sample_labels 中可用的 *_positive 列并报错。"""
    if not args.virus:
        return args.group_col, args.pos_level, args.neg_level
    import derive_virus_labels as dvl
    vcol = dvl.sanitize_virus(args.virus) + "_positive"
    with open(labels_path, encoding="utf-8") as fh:
        header = fh.readline().rstrip("\n").split("\t")
    if vcol not in header:
        avail = [c for c in header if c.endswith("_positive")]
        raise FileNotFoundError(
            f"--virus {args.virus!r}: sample_labels 无 {vcol} 列; 可用: "
            f"{avail or '(无 — 04_Labels 需来自 ④ 检出表推导, 外部分组表请自带该列)'}")
    return vcol, "1", "0"


def _run_modules_substep(args, dirs, log) -> None:
    """coexp/modules 子步: DEG 限定共表达模块检测 (KingdomStress Fig5 流程)。

    需要 06_DE 的 DEG 表 — coexp 先于 de 运行时此处优雅跳过, stage_de
    完成后会自动补跑; 自带独立断点 (coexp_modules.ok)。"""
    cdir = dirs["coexp"]
    ck = cdir / ".checkpoints" / "coexp_modules.ok"
    ck.parent.mkdir(parents=True, exist_ok=True)
    prod = cdir / "modules" / "module_summary.tsv"
    deg = Path(args.module_deg) if args.module_deg \
        else dirs["de"] / "deseq2_degs.tsv"
    if not deg.is_file():
        log.warning("[coexp/modules] 无 DEG 表 (%s); 模块检测跳过 "
                    "(06_DE 完成后会自动补跑, 或 --module-deg 指定)", deg)
        return
    if not args.force and ck.is_file() and prod.exists():
        log.info("[coexp/modules] 已完成, 跳过")
        return
    cmd = [sys.executable, str(SCRIPT_DIR / "coexp_modules.py"),
           "--deg", str(deg),
           "--z-threshold", f"{args.module_z:g}",
           "--resolution", f"{args.module_resolution:g}",
           "--output_dir", str(cdir / "modules")]
    for grp, work in (("positive", cdir / "teagcn_positive"),
                      ("negative", cdir / "teagcn_negative")):
        cmd += ["--gcn-dir", f"{grp}={_find_completed_gcn(work)}"]
    qc_path = cdir / "network_qc.tsv"
    if qc_path.is_file():
        cmd += ["--network-qc", str(qc_path)]
    vlabels = dirs["labels"] / "virus_labels.tsv"
    if vlabels.is_file():
        cmd += ["--tpm", str(dirs["quant"] / "matrix" / "TPM_matrix.tsv"),
                "--virus-labels", str(vlabels)]
        passed = dirs["qc"] / "qc_passed_samples.txt"
        if args.use_qc and passed.is_file():
            cmd += ["--samples", str(passed)]
    run_cmd(cmd, log, "coexp/modules")
    ck.write_text("ok\n", encoding="utf-8")


def stage_coexp(args, dirs, log) -> None:
    """TEA-GCN 差异共表达 (主分析): 切矩阵 → 双网络 ×2 → 秩变换 ×2 → z-QC →
    Δz → DEG 限定模块检测 (DEG 缺时 modules 由 06_DE 后补跑)。"""
    if args.teagcn_entry == "mp" and not hasattr(os, "fork"):
        msg = ("teagcn mp 入口依赖 fork (仅 Linux); 本平台请用 "
               "--teagcn-entry serial 或 fast")
        if args.dry_run:
            log.warning("[dry-run] %s (预览继续)", msg)
        else:
            raise EnvironmentError(msg)
    py = args.teagcn_python or sys.executable
    cdir = dirs["coexp"]
    output_dir = cdir.parent
    gate_fail = "skip" if "all" in args.stage else "error"

    # 分组口径: 默认 Virus_any +/-; --virus 切换为该病毒 1 vs 0
    gcol, plev, nlev = args.group_col, args.pos_level, args.neg_level
    if args.virus:
        if args.dry_run:
            try:
                gcol, plev, nlev = _coexp_group_spec(
                    args, dirs["labels"] / "sample_labels.tsv")
            except FileNotFoundError as exc:
                log.warning("[dry-run] %s", exc)
                gcol, plev, nlev = f"{args.virus}_positive", "1", "0"
        else:
            gcol, plev, nlev = _coexp_group_spec(
                args, dirs["labels"] / "sample_labels.tsv")

    matrices = cdir / "matrices"
    work_pos, work_neg = cdir / "teagcn_positive", cdir / "teagcn_negative"
    entry = {"mp": "main/Run_TEA-GCN_mp.py", "serial": "main/Run_TEA-GCN.py",
             "fast": "main/Run_TEA-GCN_fast.py"}[args.teagcn_entry]

    def build_matrices():
        cmd = [py, str(SCRIPT_DIR / "build_coexp_matrices.py"),
               "--tpm", str(dirs["quant"] / "matrix" / "TPM_matrix.tsv"),
               "--labels", str(dirs["labels"] / "sample_labels.tsv"),
               "--output_dir", str(matrices),
               "--group-col", gcol,
               "--pos-level", plev, "--neg-level", nlev,
               "--min-group", str(args.coexp_min_group),
               "--max-imbalance", str(args.coexp_max_imbalance),
               "--on-gate-fail", gate_fail]
        if args.use_qc and (dirs["qc"] / "qc_passed_samples.txt").is_file():
            cmd += ["--samples", str(dirs["qc"] / "qc_passed_samples.txt")]
        if args.force:
            cmd.append("--force")
        run_cmd(cmd, log, "coexp/matrices")

    if args.dry_run:
        # 预览: 依次打印 7 个子步计划 (k/pc 的实际收缩在运行时按组规模定)
        log.info("[dry-run] 1/7 matrices: build_coexp_matrices.py "
                 "(group-col=%s %s/%s, 闸门 >=%d 且 <=%.1f×, %s 模式)",
                 gcol, plev, nlev,
                 args.coexp_min_group, args.coexp_max_imbalance, gate_fail)
        log.info("[dry-run] 2/7 partitions: Generate_partitions.py ×2 "
                 "(k 扫描 %d–%d, 运行时按组规模收缩)",
                 args.coexp_ks, args.coexp_ke)
        log.info("[dry-run] 3/7 networks: %s ×2 (--resume, workers=%d)",
                 entry, args.coexp_workers)
        log.info("[dry-run] 4/7 ranks: rank_transform_mp.py ×2 "
                 "(fork 并行 ×%d; 输出与串行逐字节一致)", args.rank_workers)
        log.info("[dry-run] 5/7 qc: coexp_network_qc.py ×2 网络 "
                 "(z 分布: FAIL<1.5 / LOW<1.7; Δz 闸门=%s)", args.coexp_qc_gate)
        log.info("[dry-run] 6/7 diff: diff_coexp.py (top-n=%d)", args.coexp_top_n)
        log.info("[dry-run] 7/7 modules: coexp_modules.py "
                 "(z>=%.1f ∩ DEG → Louvain res=%.1f; DEG 缺则跳过, "
                 "06_DE 完成后自动补跑)", args.module_z, args.module_resolution)
        return

    _coexp_substep(
        "matrices", [matrices / "positive.tsv", matrices / "group_summary.tsv"],
        build_matrices, args, output_dir, log)

    import csv as _csv
    summary = {}
    sum_path = matrices / "group_summary.tsv"
    if sum_path.is_file():
        with sum_path.open(encoding="utf-8") as fh:
            rdr = list(_csv.DictReader(fh, delimiter="\t"))
            summary = rdr[0] if rdr else {}
    if summary.get("status") == "SKIPPED":
        log.warning("[coexp] 闸门未过, 阶段跳过 (%s); 06/07 阶段继续",
                    summary.get("reason", ""))
        return
    if not (matrices / "positive.tsv").is_file():
        raise FileNotFoundError(f"缺 {matrices / 'positive.tsv'}; 请先跑 matrices 子步")

    with (matrices / "positive.tsv").open(encoding="utf-8") as fh:
        n_genes = sum(1 for _ in fh) - 1
    with sum_path.open(encoding="utf-8") as fh:
        n_min = min(int(summary.get("n_positive") or 0),
                    int(summary.get("n_negative") or 0)) or 2
    ke = max(2, min(args.coexp_ke, n_min - 1))
    ks = max(2, min(args.coexp_ks, ke))
    pc = max(2, min(100, n_min - 1, n_genes))
    blas_t = min(args.threads, 64)
    est_gb = 2 * 70 * n_genes ** 2 / 1e9
    if est_gb > 50:
        log.warning("[coexp] 预计磁盘 ~%.0f GB (双网络 RAW+Completed)", est_gb)
    log.info("[coexp] 参数: 分组 %s %s/%s, 基因 %d, 组规模 pos/neg (min %d), "
             "k %d–%d, pc %d, entry=%s, workers=%d", gcol, plev, nlev,
             n_genes, n_min, ks, ke, pc, args.teagcn_entry, args.coexp_workers)

    def run_partitions():
        for tag, work, mat in (("pos", work_pos, matrices / "positive.tsv"),
                               ("neg", work_neg, matrices / "negative.tsv")):
            run_cmd([py, "main/Generate_partitions.py", "-o", work, "-im", mat,
                     "-ks", str(ks), "-ke", str(ke), "-st", "1",
                     "-t", str(blas_t), "-pc", str(pc)],
                    log, f"coexp/partitions[{tag}]", cwd=TEAGCN_DIR)

    _coexp_substep("partitions",
                   [work_pos / "Generate_partitions" / "selected_k.pkl",
                    work_neg / "Generate_partitions" / "selected_k.pkl"],
                   run_partitions, args, output_dir, log)

    def run_networks():
        for tag, work, mat in (("pos", work_pos, matrices / "positive.tsv"),
                               ("neg", work_neg, matrices / "negative.tsv")):
            run_cmd([py, entry, "-o", work, "-im", mat, "-k", "0",
                     "-t", str(blas_t), "-w", str(args.coexp_workers), "--resume"],
                    log, f"coexp/networks[{tag}]", cwd=TEAGCN_DIR)

    _coexp_substep("networks",
                   [work_pos / "RAW_GCN", work_neg / "RAW_GCN"],
                   run_networks, args, output_dir, log)

    def run_ranks():
        # 并行版 (fork 分片, 输出与串行逐字节一致, 见 rank_transform_mp.py);
        # 无 fork 平台脚本内部自动退化为单进程
        for tag, work in (("pos", work_pos), ("neg", work_neg)):
            run_cmd([py, str(SCRIPT_DIR / "rank_transform_mp.py"),
                     "-o", work, "-w", str(args.rank_workers)],
                    log, f"coexp/ranks[{tag}]", cwd=TEAGCN_DIR)

    _coexp_substep("ranks",
                   [work_pos / "Completed_GCN", work_neg / "Completed_GCN"],
                   run_ranks, args, output_dir, log)

    pos_gcn = _find_completed_gcn(work_pos)
    neg_gcn = _find_completed_gcn(work_neg)

    # 网络 z 分布 QC (KingdomStress Fig S8 物种剔除规则的工程化)
    qc_path = cdir / "network_qc.tsv"

    def run_qc():
        run_cmd([py, str(SCRIPT_DIR / "coexp_network_qc.py"),
                 "--gcn-dir", f"positive={pos_gcn}",
                 "--gcn-dir", f"negative={neg_gcn}",
                 "--group-summary", str(matrices / "group_summary.tsv"),
                 "--output_dir", str(cdir)], log, "coexp/qc")

    _coexp_substep("qc", [qc_path], run_qc, args, output_dir, log)

    def run_diff():
        run_cmd([py, str(SCRIPT_DIR / "diff_coexp.py"),
                 "--pos-gcn", pos_gcn, "--neg-gcn", neg_gcn,
                 "--output_dir", str(cdir), "--top-n", str(args.coexp_top_n)],
                log, "coexp/diff")

    # 网络 z-QC 闸门 (KingdomStress Fig S8): FAIL 网络的 Δz 跨网络比较不可信
    skip_diff = False
    if qc_path.is_file() and args.coexp_qc_gate != "off":
        with qc_path.open(encoding="utf-8") as fh:
            qc_rows = list(_csv.DictReader(fh, delimiter="\t"))
        fails = [r.get("Group", "?") for r in qc_rows if r.get("Status") == "FAIL"]
        if fails:
            msg = (f"网络 z 分布截断 FAIL: {fails} (max_z<1.5, Fig S8 — "
                   f"样本量不足压低 MR z 分布)")
            if args.coexp_qc_gate == "error":
                raise RuntimeError(msg + "; 修复组规模/ QC 后重跑, "
                                        "或 --coexp-qc-gate warn/skip")
            log.warning("[coexp] %s — 闸门=%s", msg, args.coexp_qc_gate)
            if args.coexp_qc_gate == "skip":
                log.warning("[coexp] Δz 按闸门跳过 (下次运行自动重试)")
                skip_diff = True
    if not skip_diff:
        _coexp_substep("diff",
                       [cdir / "diff_coexp_significant.tsv",
                        cdir / "diff_coexp_summary.tsv"],
                       run_diff, args, output_dir, log)

    _run_modules_substep(args, dirs, log)


def stage_de(args, dirs, log) -> None:
    counts = dirs["quant"] / "matrix" / "NumReads_matrix.tsv"
    labels = dirs["labels"] / "sample_labels.tsv"
    cmd = [sys.executable, str(SCRIPT_DIR / "run_deseq2.py"),
           "--counts", str(counts), "--labels", str(labels),
           "--output_dir", str(dirs["de"]),
           "--group-col", args.group_col,
           "--ref-level", args.ref_level, "--alt-level", args.alt_level,
           "--padj", str(args.padj), "--lfc", str(args.lfc),
           "-t", str(args.threads)]
    if args.de_strata:
        cmd += ["--strata", args.de_strata]
    if args.use_qc:
        passed = dirs["qc"] / "qc_passed_samples.txt"
        if passed.is_file():
            cmd += ["--samples", str(passed)]
        else:
            log.warning("--use-qc 已开但 %s 不存在; 用全部有标签样本", passed)
    if args.dry_run:
        log.info("[dry-run] %s", " ".join(cmd))
        return
    run_cmd(cmd, log, "de")

    # coexp 先于 de 运行时 modules 因缺 DEG 被跳过 — DE 就绪后自动补跑
    if (dirs["coexp"] / "teagcn_positive" / "Completed_GCN").is_dir():
        log.info("[de] DE 就绪 — 补跑 coexp/modules (DEG 限定模块检测)")
        try:
            _run_modules_substep(args, dirs, log)
        except (subprocess.CalledProcessError, RuntimeError,
                FileNotFoundError) as exc:
            log.warning("[coexp/modules] 补跑失败 (不影响 DE 产物): %s", exc)


def stage_loadscan(args, dirs, log) -> None:
    tpm = dirs["quant"] / "matrix" / "TPM_matrix.tsv"
    labels = dirs["labels"] / "virus_labels.tsv"
    if not labels.is_file() and not args.dry_run:
        log.warning("[loadscan] 无 virus_labels.tsv (通用分组模式或 ④ 未跑); "
                    "该阶段为病毒载量专属, 跳过 — 05_Coexp/06_DE 不受影响")
        return
    cmd = [sys.executable, str(SCRIPT_DIR / "viral_load_scan.py"),
           "--tpm", str(tpm), "--labels", str(labels),
           "--output_dir", str(dirs["loadscan"]),
           "--strata-col", args.strata_col]
    slabels = dirs["labels"] / "sample_labels.tsv"
    if slabels.is_file():
        cmd += ["--meta", str(slabels)]
    if args.use_qc:
        passed = dirs["qc"] / "qc_passed_samples.txt"
        if passed.is_file():
            cmd += ["--samples", str(passed)]
    if args.tissue:
        cmd += ["--tissue", args.tissue]
    if args.dry_run:
        log.info("[dry-run] %s", " ".join(cmd))
        return
    run_cmd(cmd, log, "loadscan")


def _collect_enrich_sets(dirs, out_dir: Path, hub_n: int) -> dict:
    """从 ⑧ 产物自动收集基因集, 归一化为 07_Enrich/sets/<name>.txt。

    来源: 06_DE 显著 DEG (全量) / 05_Coexp Δz hub (Top-N, aggregate 已按
    显著对数排序) / modules 每模块成员。返回 name→文件路径 (可能为空)。"""
    sets_dir = Path(out_dir) / "sets"
    sets_dir.mkdir(parents=True, exist_ok=True)
    found = {}

    degs = dirs["de"] / "deseq2_degs.tsv"
    if degs.is_file():
        genes, up, down = [], [], []
        with degs.open(encoding="utf-8") as fh:
            hdr = fh.readline().rstrip("\n").split("\t")
            try:
                ilfc = hdr.index("log2FoldChange")
            except ValueError:
                ilfc = -1
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if not parts or not parts[0]:
                    continue
                genes.append(parts[0])
                if ilfc >= 0 and len(parts) > ilfc:
                    try:
                        v = float(parts[ilfc])
                    except ValueError:
                        v = 0.0
                    (up if v > 0 else down).append(parts[0])
        for name, gs in (("DE", genes), ("DE_up", up), ("DE_down", down)):
            if gs:
                p = sets_dir / f"{name}.txt"
                p.write_text("\n".join(gs) + "\n", encoding="utf-8")
                found[name] = p

    hubs = dirs["coexp"] / "diff_coexp_gene_aggregate.tsv"
    if hubs.is_file():
        genes = [l.split("\t", 1)[0] for l in
                 hubs.read_text(encoding="utf-8").splitlines()[1:]
                 if l.strip()][:hub_n]
        if genes:
            p = sets_dir / "dz_hub.txt"
            p.write_text("\n".join(genes) + "\n", encoding="utf-8")
            found["dz_hub"] = p

    mod_dir = dirs["coexp"] / "modules"
    if mod_dir.is_dir():
        for f in sorted(mod_dir.glob("module_table_*.tsv")):
            group = f.stem.replace("module_table_", "")
            members = {}
            with f.open(encoding="utf-8") as fh:
                next(fh, None)
                for line in fh:
                    parts = line.rstrip("\n").split("\t")
                    if len(parts) >= 2 and parts[0]:
                        members.setdefault(parts[1], []).append(parts[0])
            for mid, genes in sorted(members.items()):
                name = f"module_{group}_{mid}"
                p = sets_dir / f"{name}.txt"
                p.write_text("\n".join(genes) + "\n", encoding="utf-8")
                found[name] = p
    return found


def stage_enrich(args, dirs, log) -> None:
    """07_Enrich 功能富集: 自动收集 05/06 基因集 + TPM 背景 → 超几何 ORA。

    注释源 (--annotation, eggNOG-mapper 或通用映射) 是物种相关的外部输入,
    未提供时优雅跳过 (与 loadscan 无 virus_labels 同策略)。"""
    if not args.annotation:
        log.warning("[enrich] 未提供 --annotation (eggNOG-mapper 注释或通用 "
                    "gene/term 映射表); 功能富集跳过 — 见模块 README §9")
        return
    tpm = dirs["quant"] / "matrix" / "TPM_matrix.tsv"
    if not tpm.is_file() and not args.dry_run:
        raise FileNotFoundError(f"缺 TPM 矩阵: {tpm} (请先跑 quant)")

    edir = dirs["enrich"]
    background = edir / "background_genes.txt"
    gene_sets = _collect_enrich_sets(dirs, edir, args.enrich_hub_n)
    if not gene_sets:
        log.warning("[enrich] 无可用基因集 (需 06_DE 的 DEG 或 05_Coexp 的 "
                    "hub/模块产物); 功能富集跳过")
        return

    cmd = [sys.executable, str(SCRIPT_DIR / "enrich_analysis.py"),
           "--annotation", str(Path(args.annotation).resolve()),
           "--background", str(background),
           "--output_dir", str(edir)]
    for name, path in gene_sets.items():
        cmd += ["--gene-set", f"{name}={path}"]
    if args.dry_run:
        log.info("[dry-run] %s", " ".join(cmd))
        return
    edir.mkdir(parents=True, exist_ok=True)
    with tpm.open(encoding="utf-8") as fin, \
            background.open("w", encoding="utf-8") as fout:
        next(fin, None)
        for line in fin:
            gene = line.split("\t", 1)[0].strip()
            if gene:
                fout.write(gene + "\n")
    run_cmd(cmd, log, "enrich")


def stage_report(args, dirs, log) -> None:
    # FIX(§5 P0-3): 原编排器从不调用 host_expression_plots.py，导致 `--stage all`
    #   不会生成 09_Report/figures/fig01–fig05，报告内嵌图（fig01–fig05）全部缺失。
    #   在生成 HTML 报告前先出图。出图失败仅告警、不阻断报告（HTML 仍可用）。
    plots_cmd = [sys.executable, str(SCRIPT_DIR / "host_expression_plots.py"),
                 "--root", str(args.output_dir)]
    if args.dry_run:
        log.info("[dry-run] %s", " ".join(plots_cmd))
    else:
        try:
            run_cmd(plots_cmd, log, "plots")
        except Exception as exc:  # noqa: BLE001 - 出图失败不应使整条流水线失败
            log.warning("plots 阶段失败（报告仍将生成，图可能缺失）: %s", exc)

    cmd = [sys.executable, str(SCRIPT_DIR / "host_expression_report.py"),
           "--root", str(args.output_dir)]
    if args.dry_run:
        log.info("[dry-run] %s", " ".join(cmd))
        return
    run_cmd(cmd, log, "report")


STAGE_FUNCS = {
    "hostref": stage_hostref, "quant": stage_quant, "qc": stage_qc,
    "labels": stage_labels, "coexp": stage_coexp, "de": stage_de,
    "loadscan": stage_loadscan, "enrich": stage_enrich, "report": stage_report,
}


# ═══════════════════════════════════════════════════════════════════

def _expand_stages(requested: Sequence[str], with_loadscan: bool) -> List[str]:
    """展开 --stage 列表 (all → STAGES, 去重保序)。

    08_LoadScan 剂量验证层默认含于 'all'; --no-with-loadscan 可关
    (显式点名 loadscan 永远生效)。"""
    stages: List[str] = []
    for s in requested:
        stages.extend(STAGES if s == "all" else [s])
    if "all" in requested and not with_loadscan and "loadscan" not in requested:
        stages = [s for s in stages if s != "loadscan"]
    seen, ordered = set(), []
    for s in stages:
        if s not in seen:
            seen.add(s)
            ordered.append(s)
    return ordered


def parse_args():
    p = argparse.ArgumentParser(
        description="⑧ 宿主表达管线 (可选补充模块) 编排器",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--stage", nargs="+", default=["all"],
                   choices=list(STAGES) + ["all"],
                   help="运行阶段 (all = 全流程, 可多选)")
    p.add_argument("--output_dir", "-o",
                   help="⑧ 输出根 (08_HostExpression); 缺省时用 "
                        "<--project-root>/08_HostExpression")
    p.add_argument("--project-root",
                   help="项目根: 用于自动定位 ① 宿主参考/元数据、② clean reads、"
                        "④ 检出表 (standard 优先, legacy 兜底)")
    p.add_argument("--io-layout", choices=["legacy", "standard"], default=None,
                   help="I/O 布局 (缺省: 环境变量 MMPV_IO_LAYOUT 或 standard)")

    g = p.add_argument_group("输入 (缺省时经 --project-root 自动定位)")
    g.add_argument("--reads_dir", help="clean reads 目录 (② 产物)")
    g.add_argument("--reference", help="宿主 CDS/转录本 FASTA (跳过提取)")
    g.add_argument("--genome", help="宿主基因组 FASTA (hostref 阶段)")
    g.add_argument("--gff", help="宿主注释 GFF3 (hostref 阶段)")
    g.add_argument("--detect", help="病毒检出表或目录 (④ 产物)")
    g.add_argument("--metadata", help="Core14 元数据 (合并 Tissue)")

    g = p.add_argument_group("工具与资源")
    g.add_argument("--salmon", default=os.environ.get("MMPV_SALMON", "salmon"),
                   help="salmon 可执行文件")
    g.add_argument("--threads", "-t", type=int, default=60, help="总线程数")
    g.add_argument("--jobs", "-j", type=int, default=8, help="并行样本数")

    g = p.add_argument_group("分析参数")
    g.add_argument("--cov-threshold", type=float, default=90.0,
                   help="阳性判定: 基因组覆盖度阈值 (%%)")
    g.add_argument("--min-mapped", type=float, default=50.0,
                   help="QC: percent_mapped 阈值 (%%)")
    g.add_argument("--min-processed", type=float, default=1e6,
                   help="QC: 处理 reads 数阈值")
    g.add_argument("--use-qc", action=argparse.BooleanOptionalAction,
                   default=True,
                   help="下游分析只用 QC 通过样本 (默认开)")
    g.add_argument("--group-col", default="Virus_any",
                   help="分组列 (coexp/de 共用; 默认病毒口径 Virus_any, "
                        "也可用 <病毒>_positive 或任意二分列如 Condition)")
    g.add_argument("--labels-file", default=None,
                   help="通用模式: 直接提供分组表 (Sample+分组列 的 TSV, "
                        "可含 Tissue 等元数据列), 跳过 ④ 检出表推导 — "
                        "非病毒项目/任意二分组比较用")
    g.add_argument("--ref-level", default="-",
                   help="DE 对照组水平 (默认 '-' = 病毒阴性)")
    g.add_argument("--alt-level", default="+",
                   help="DE 处理组水平 (默认 '+' = 病毒阳性)")
    g.add_argument("--pos-level", default="+",
                   help="coexp 阳性组水平 (默认 '+')")
    g.add_argument("--neg-level", default="-",
                   help="coexp 阴性组水平 (默认 '-')")
    g.add_argument("--strata-col", default="Tissue",
                   help="分层列名 (DE --de-strata 与 LoadScan 组织内扫描共用; "
                        "默认 Core14 的 Tissue, 可换任意元数据列)")
    g.add_argument("--de-strata", default=None,
                   help="DE 分层开关: 给定则按 --strata-col 列分层内各做一次对比, "
                        "公共队列建议开启 (组织混杂防护)")
    g.add_argument("--padj", type=float, default=0.05, help="DE padj 阈值")
    g.add_argument("--lfc", type=float, default=1.0, help="DE |log2FC| 阈值")
    g.add_argument("--tissue", help="LoadScan 范围过滤 (按 --strata-col 列子串匹配, 逗号分隔; "
                                    "不给时若该列存在则自动加做组内扫描)")

    g = p.add_argument_group("coexp 参数 (TEA-GCN Δz 主分析)")
    g.add_argument("--coexp-min-group", type=int, default=30,
                   help="闸门: 每组样本数下限 (默认 30)")
    g.add_argument("--virus", default=None,
                   help="coexp 指定病毒: 分组切为 <病毒>_positive 1 vs 0 "
                        "(默认 Virus_any 阳性 vs 阴性; 06_DE 不受影响, "
                        "仍走 --group-col; 切换口径会覆盖 05_Coexp 现有产物)")
    g.add_argument("--coexp-max-imbalance", type=float, default=3.0,
                   help="闸门: 两组样本量比上限 (默认 3.0)")
    g.add_argument("--coexp-ks", type=int, default=5,
                   help="TEA-GCN k 扫描起点 (自动收缩到组规模允许的上限)")
    g.add_argument("--coexp-ke", type=int, default=30,
                   help="TEA-GCN k 扫描终点")
    g.add_argument("--coexp-top-n", type=int, default=20,
                   help="Δz: 每基因取 |Δz| Top-N 邻居")
    g.add_argument("--coexp-workers", type=int, default=0,
                   help="TEA-GCN mp 进程数 (0 = threads//8)")
    g.add_argument("--rank-workers", type=int, default=0,
                   help="Rank_transform 并行进程数 (0 = max(8, threads//4); "
                        "输出与串行版逐字节一致)")
    g.add_argument("--teagcn-entry", choices=["mp", "serial", "fast"], default="mp",
                   help="TEA-GCN 网络构建入口 (mp=fork 多进程, 仅 Linux)")
    g.add_argument("--teagcn-python", default=None,
                   help="运行 teagcn 子脚本的 Python (默认当前解释器; "
                        "需含 einsumt — pixi 环境已带)")
    g.add_argument("--coexp-qc-gate", choices=["warn", "skip", "error", "off"],
                   default="warn",
                   help="Δz 前网络 z-QC 闸门 (max_z<1.5 判 FAIL, Fig S8): "
                        "warn=警告后继续 / skip=跳过 Δz / error=中止 / off=不检查")
    g.add_argument("--module-z", type=float, default=1.7,
                   help="模块检测边阈值 (KingdomStress Fig5: z>=1.7)")
    g.add_argument("--module-resolution", type=float, default=1.0,
                   help="Louvain resolution (论文 Fig5: 1.0)")
    g.add_argument("--module-deg", default=None,
                   help="模块检测 DEG 表 (缺省 06_DE/deseq2_degs.tsv)")
    g.add_argument("--annotation", default=None,
                   help="07_Enrich 功能富集注释: eggNOG-mapper "
                        ".emapper.annotations (GOs/KEGG_ko) 或通用 gene/term "
                        "映射 TSV (MapMan/Mercator 导出亦可); 缺省跳过富集")
    g.add_argument("--enrich-hub-n", type=int, default=50,
                   help="富集基因集: Δz hub 取前 N (默认 50)")

    g = p.add_argument_group("流程控制")
    g.add_argument("--force", action="store_true", help="强制重跑 (忽略断点)")
    g.add_argument("--with-loadscan", action=argparse.BooleanOptionalAction,
                   default=True,
                   help="08_LoadScan 剂量验证层 (默认含于 --stage all; "
                        "--no-with-loadscan 可关)")
    g.add_argument("--dry-run", action="store_true", help="只打印将执行的命令")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    normalize_layout_env(args.io_layout)
    if not args.coexp_workers:
        args.coexp_workers = max(1, args.threads // 8)
    if not args.rank_workers:
        args.rank_workers = max(8, args.threads // 4)

    if not args.output_dir:
        if not args.project_root:
            print("[ERROR] 需要 --output_dir, 或提供 --project-root 以推导 "
                  "<项目根>/08_HostExpression", file=sys.stderr)
            sys.exit(2)
        args.output_dir = str(Path(args.project_root)
                              / PIPELINE_ROOT_NAMES["hostexpr"])
    output_dir = Path(args.output_dir).resolve()

    log = setup_logger(output_dir, log_file=not args.dry_run)
    log.info("⑧ 宿主表达管线 (可选补充模块) | 输出根: %s", output_dir)

    stages = _expand_stages(args.stage, args.with_loadscan)
    if "all" in args.stage and "loadscan" not in stages:
        log.info("[loadscan] 已通过 --no-with-loadscan 关闭")
    log.info("阶段: %s", " → ".join(stages))

    dirs = build_hostexpr_dirs(output_dir)
    resolved_reference: Optional[Path] = None

    failed = False
    for stage in stages:
        products = stage_products(dirs, stage, resolved_reference)
        if not args.dry_run and not args.force \
                and stage != "hostref" \
                and stage_done(output_dir, stage, products):
            log.info("[%s] 已完成 (产物齐全), 跳过; --force 可重跑", stage)
            continue
        log.info("── [%s] ──", stage)
        try:
            result = STAGE_FUNCS[stage](args, dirs, log)
        except (subprocess.CalledProcessError, FileNotFoundError,
                EnvironmentError, RuntimeError) as exc:
            log.error("[%s] 失败: %s", stage, exc)
            failed = True
            break
        if stage == "hostref":
            resolved_reference = result if isinstance(result, Path) else None
        elif resolved_reference is None:
            try:
                resolved_reference = resolve_reference(args, dirs, log)
            except FileNotFoundError:
                pass
        if not args.dry_run and stage != "hostref":
            mark_done(output_dir, stage)

    if failed:
        if args.dry_run:
            log.error("dry-run 预览中断 (存在输入定位失败, 见上方 ERROR)")
            sys.exit(1)
        log.error("管线中断; 修复后重发同命令即可续跑 (已完成阶段自动跳过)")
        sys.exit(1)
    if args.dry_run:
        log.info("dry-run 结束 (未执行任何命令)")
        return

    # 边界留痕: 记录本模块消费的主链产物与关键输出 (幂等; 路径一律绝对)
    reads_entry = Path(args.reads_dir).resolve() if args.reads_dir else None
    if reads_entry is None and args.project_root:
        reads_entry = locate_clean_reads(Path(args.project_root))
    detect_entry = Path(args.detect).resolve() if args.detect else None
    if detect_entry is None and args.project_root:
        detect_entry = locate_detect_table(Path(args.project_root))
    reference_entry = resolved_reference
    if reference_entry is None:
        cand = dirs["hostref"] / "host_cds.fasta"
        reference_entry = cand if cand.is_file() else None

    entries = {}
    for key, path in (("reads_dir", reads_entry),
                      ("reference", reference_entry),
                      ("detect_table", detect_entry),
                      ("tpm_matrix", dirs["quant"] / "matrix" / "TPM_matrix.tsv"),
                      ("sample_labels", dirs["labels"] / "sample_labels.tsv"),
                      ("deseq2_all_results", dirs["de"] / "deseq2_all_results.tsv")):
        if path and Path(path).exists():
            entries[key] = str(path)
    if entries:
        write_boundary(output_dir, "hostexpr", entries, "host_expression.py")

    log.info("全部阶段完成 ✓  报告: %s",
             dirs["report"] / "host_expression_report.html")


if __name__ == "__main__":
    main()
