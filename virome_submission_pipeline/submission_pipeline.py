#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
submission_pipeline.py — MMPV-RNA → NCBI 提交统一入口 v1.0
==========================================================

接收上游 MMPV-RNA 管道 9b (analysis_verify) 的 class_KEEP.fasta, 自动检测已有文件,
按需运行 suvtk taxonomy/features/hypothetical,
最终生成两种格式的提交文件。

上游输入:
  09b_Analysis_Verify/class_KEEP.fasta     (默认输入, KEEP 高置信子集 ≥1000bp)
  08_Rescue/
  ├── HQ_plant_viruses.fasta          (回退全集, 未运行 9b 时)
  ├── suvtk.taxonomy_output/          (可选, 已存在则跳过)
  ├── suvtk.features_output/          (可选, 已存在则跳过)
  ├── analyze_hypothetical/           (可选, 已存在则跳过)
  └── metadata_output/                (可选, 自动填充 source.src)

两种输出模式:
  A) suvtk 模式: source.src + suvtk comments → table2asn → .sqn
  B) Sequin 模式: .fsa + .tbl + .cmt → tbl2asn → .sqn (Cenote-Taker3 风格)

用法:
  # 从 08_Rescue 一键运行
  python submission_pipeline.py \\
      --work-dir $OUT/08_Rescue/ \\
      --run-title ningxiagouqi_plant_virome \\
      --mode both \\
      --suvtk-db ~/database/virus-db/suvtk_db/ \\
      --rvdb-db ~/database/virus-db/RVDB-v31/ \\
      -t 40
"""

import argparse
import os
import re
import sys
import shutil
import subprocess
import logging

# 跨管线统一 I/O 布局 (mmpv_common/, 仓库根): ③ 输出目录名随 MMPV_IO_LAYOUT 解析
_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)
from mmpv_common.io_layout import layout_dirs as _layout_dirs, normalize_layout_env
_D = _layout_dirs(os.environ.get("MMPV_IO_LAYOUT"))
from pathlib import Path
from datetime import datetime
from Bio import SeqIO

SCRIPT_DIR = Path(__file__).resolve().parent

# bash 解析: Linux= /bin/bash; Windows(Git Bash) 经 shutil.which 找到 bash.exe。
# Git 常装在带空格路径, 故统一 [BASH, '-c', cmd] 列表式调用 (shell=True+executable 会断)。
BASH = shutil.which('bash') or '/bin/bash'


def _read_tsv_simple(path):
    """读取 TSV 文件, 返回 [{col: val}]"""
    rows = []
    if not Path(path).is_file():
        return rows
    with open(path) as f:
        hdr = f.readline().strip().split('\t')
        for line in f:
            if not line.strip(): continue
            rows.append(dict(zip(hdr, line.strip().split('\t'))))
    return rows


def _count_lines(path):
    """计算文件行数"""
    if not Path(path).is_file(): return 0
    with open(path) as f:
        return sum(1 for _ in f)


# ══════════════════════════════════════════════════════════════
# 管道参数自动检测 (从 pipeline_config.yaml)
# ══════════════════════════════════════════════════════════════

def load_pipeline_params(work_dir):
    """从 pipeline_config.yaml 自动提取组装/测序/富集参数"""
    import yaml as _yaml
    params = {
        'assembler': 'MEGAHIT;1.2.9;default parameters',
        'sequencer': 'Illumina NovaSeq 6000',
        'enrichment': 'none',
        'source_uvig': 'metatranscriptome (not viral targeted)',  # 默认公共转录组
        'metagenome_source': 'plant virome',
        'pipeline': 'MMPV-RNA v2.3 + suvtk v0.1.1',
    }

    # 从 work_dir 往上找 pipeline_config.yaml
    for parent in [work_dir, work_dir.parent, work_dir.parent.parent,
                   work_dir.parent.parent.parent]:
        config_path = parent / 'pipeline_config.yaml'
        if not config_path.exists():
            continue
        try:
            with open(config_path) as f:
                config = _yaml.safe_load(f)
            profiles = config.get('profiles', {})
            default = profiles.get('default', {})
            # 组装
            assembly = default.get('assembly', {})
            if assembly.get('assembler'):
                params['assembler'] = f"{assembly['assembler']};auto;default parameters"
            # 工具
            tools = default.get('tools', {})
            if tools.get('diamond'):
                pass  # 已全局使用
            # 运行时
            runtime = default.get('runtime', {})
            # 默认用 public data 的 source_uvig
            break
        except Exception:
            continue

    return params


# ══════════════════════════════════════════════════════════════
# 日志
# ══════════════════════════════════════════════════════════════

def setup_logger(out_dir):
    logger = logging.getLogger("SubPipeline")
    logger.setLevel(logging.DEBUG)
    logger.handlers.clear()
    os.makedirs(out_dir, exist_ok=True)

    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    ch.setFormatter(logging.Formatter('[%(asctime)s] %(message)s', datefmt='%H:%M:%S'))
    logger.addHandler(ch)

    fh = logging.FileHandler(os.path.join(out_dir, 'submission_pipeline.log'), encoding='utf-8')
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(logging.Formatter('[%(asctime)s] %(levelname)s %(message)s'))
    logger.addHandler(fh)

    return logger


def run(cmd, log, step_name, check=False):
    """执行命令, 容忍失败

    返回 True 仅当子进程 exit code == 0。check=True 时非零退出抛
    CalledProcessError (行为保持); check=False 时不再无条件返回 True,
    避免 table2asn/suvtk 失败却被当成成功 (stderr 尾部会打印到日志)。
    """
    log.info("[%s] %s", step_name, cmd[:120])
    try:
        p = subprocess.run([BASH, '-c', cmd], check=check,
                           stderr=subprocess.PIPE, text=True)
    except subprocess.CalledProcessError as e:
        log.error("[%s] ✗ (exit=%d)", step_name, e.returncode)
        for line in (e.stderr or '').strip().splitlines()[-10:]:
            log.error("[%s]   %s", step_name, line)
        return False
    except OSError as e:
        # 例如 /bin/bash 不存在 (非 Linux 环境): 按失败处理而不是崩溃
        log.error("[%s] ✗ 无法执行: %s", step_name, e)
        return False
    if p.returncode != 0:
        log.error("[%s] ✗ (exit=%d)", step_name, p.returncode)
        for line in (p.stderr or '').strip().splitlines()[-10:]:
            log.error("[%s]   %s", step_name, line)
        return False
    if p.stderr:
        log.debug("[%s] stderr:\n%s", step_name, p.stderr[-2000:])
    log.info("[%s] ✓", step_name)
    return True


def run_argv(argv, log, step_name):
    """执行 argv 命令列表 (不经 shell): Windows 中文/空格/反斜杠路径安全。

    与 run() 同语义 (返回 True 仅当 exit code == 0), 供 report_html 等
    内部 python 步骤使用 (bash -c 会吃反斜杠路径 → 报告步骤静默失败)。
    """
    log.info("[%s] %s", step_name, " ".join(str(a) for a in argv[:6]))
    try:
        p = subprocess.run([str(a) for a in argv], stderr=subprocess.PIPE,
                           text=True, errors='replace')
    except OSError as e:
        log.error("[%s] ✗ 无法执行: %s", step_name, e)
        return False
    if p.returncode != 0:
        log.error("[%s] ✗ (exit=%d)", step_name, p.returncode)
        for line in (p.stderr or '').strip().splitlines()[-10:]:
            log.error("[%s]   %s", step_name, line)
        return False
    if p.stderr:
        log.debug("[%s] stderr:\n%s", step_name, p.stderr[-2000:])
    log.info("[%s] ✓", step_name)
    return True


def _first_existing(*paths):
    """返回第一个存在的文件路径 (全部不存在/None 时返回 None)。"""
    for p in paths:
        if p is not None and Path(p).is_file():
            return p
    return None


def _ncbi_email(log=None):
    """NCBI 联系邮箱: 优先环境变量 NCBI_EMAIL; 取不到返回空串并告警

    不在代码里硬编码任何邮箱 (含示例/真实邮箱)。
    """
    email = os.environ.get("NCBI_EMAIL", "").strip()
    if not email and log is not None:
        log.warning("  未设置 NCBI_EMAIL 环境变量: template.sbt 联系邮箱留空, 提交前必须补全")
    return email


def _find_suvtk_output(work_dir, rel_path, min_size=100):
    """在 work_dir 及其上级目录的 ③ 分析目录 (布局感知) 中查找已有 suvtk 输出"""
    # 1. work_dir 直接子目录
    candidate = work_dir / rel_path
    if candidate.exists() and os.path.getsize(candidate) > min_size:
        return candidate
    # 2. 上级目录的 09_Virome_Analysis (主流水线输出位置)
    for parent in [work_dir.parent, work_dir.parent.parent]:
        alt = parent / _D["d_analysis"] / rel_path
        if alt.exists() and os.path.getsize(alt) > min_size:
            return alt
    return None


def _detect_keep_fasta(work_dir):
    """定位 9b 的 class_KEEP.fasta (KEEP 高置信子集, ≥1000bp)。

    work_dir 可能是 $OUT 或 $OUT/08_Rescue, 逐层探测:
      <base>/09b_Analysis_Verify/class_KEEP.fasta   (新目录名)
      <base>/09b_ACVirus_Analysis/class_KEEP.fasta  (旧目录名)
    base ∈ {work_dir, work_dir.parent, work_dir.parent.parent}
    全部未命中时回退到 08_Rescue/HQ_plant_viruses.fasta 全集。
    """
    acv_names = [_D["d_verify"], "09b_ACVirus_Analysis"]   # 后者为历史旧名兜底
    bases = [work_dir, work_dir.parent, work_dir.parent.parent]
    for base in bases:
        for acv in acv_names:
            f = base / (_D["d_verify"] if acv == "09b_Analysis_Verify" else acv) / "class_KEEP.fasta"
            if f.is_file() and f.stat().st_size > 100:
                print(f"[*] 自动检测 KEEP FASTA: {f}")
                return str(f)
    # 兜底: 08_Rescue 的 HQ 全集 (未运行 9b 时)
    for base in bases:
        for hq_rel in [Path("HQ_plant_viruses.fasta"),
                       Path(_D["d_rescue"]) / "HQ_plant_viruses.fasta"]:
            hq = base / hq_rel
            if hq.is_file() and hq.stat().st_size > 100:
                print(f"[*] 回退 HQ 全集 FASTA: {hq}")
                return str(hq)
    return None


# ══════════════════════════════════════════════════════════════
# 出 .sqn 前的硬校验 gate (S1-⑦)
# ══════════════════════════════════════════════════════════════

# 位置写法: 常规两列 a..b (可含 < > 部分标记) / 旧版 complement(a..b) / join(a..b,c..d)
_COMPLEMENT_LOC_RE = re.compile(r'^complement\(\s*(.*?)\s*\)$', re.IGNORECASE)
_SEGMENT_LOC_RE = re.compile(r'[<>]?(\d+)\s*\.\.\s*[<>]?(\d+)')

# source.src 必填列 (与 suvtk table2asn 前既有校验保持一致)
_SRC_REQUIRED_COLS = ("Collection_date", "geo_loc_name", "Lat_Lon",
                      "Bioproject", "Biosample", "SRA")

_REVCOMP_TABLE = str.maketrans('ACGTNRYKMSWBDHV', 'TGCANYRMKSWVHDB')


def _revcomp(seq):
    """反向互补 (支持 IUPAC 简并碱基)"""
    return seq.translate(_REVCOMP_TABLE)[::-1]


def _placeholder_pattern():
    """占位符正则: 优先跨文件 import unified_metadata.PLACEHOLDER_RE

    导入失败 (如 SCRIPT_DIR 不在 sys.path) 时回退等价副本。
    """
    try:
        if str(SCRIPT_DIR) not in sys.path:
            sys.path.insert(0, str(SCRIPT_DIR))
        from unified_metadata import PLACEHOLDER_RE
        return PLACEHOLDER_RE
    except Exception:
        return (r'XXXX|YYYY-MM-DD|DD-Mmm-YYYY|Country:Region|PRJNAXXXX|SAMNXXXXXXXX'
                r'|XX\.[Xx0-9]{1,3}\s*N|XXX\.[Xx0-9]{1,3}\s*E|XX\.[Xx]{1,3}')


def _gate_locate(workdir, src_path, names, rel_dirs):
    """定位 gate 输入: 显式路径优先, 否则按常规目录 + 递归查找

    第二形参名与调用方 (validate_submission_gate 的 tbl_path/fna_path/src_path)
    保持一致 —— 关键字调用 `src_path=` 可用, 不再只有位置传参一条路。
    """
    if src_path:
        p = Path(src_path)
        return p if p.is_file() else None
    wd = Path(workdir)
    for rel in rel_dirs:
        for n in names:
            p = wd / rel / n
            if p.is_file():
                return p
    for n in names:
        hits = sorted(wd.rglob(n))
        if hits:
            return hits[0]
    return None


def _parse_cds_location(col1, col2=None):
    """解析 .tbl 位置列 → (segments, strand, is_partial); 无法解析返回 None

    - 常规: col1, col2 两列 (列1 > 列2 = 负链反转区间; < > 为部分标记)
    - 旧版: 单列 complement(a..b)
    - join/order(a..b,c..d): 按片段累加 (仅校验长度)
    """
    text = col1.strip() if col2 is None else f"{col1.strip()}..{col2.strip()}"
    strand = 1
    m = _COMPLEMENT_LOC_RE.match(text)
    if m:
        strand = -1
        text = m.group(1).strip()
    elif col2 is not None:
        # 5 列格式负链 = 反转区间 (列1 > 列2)
        try:
            n1 = int(col1.strip().replace('>', '').replace('<', ''))
            n2 = int(col2.strip().replace('>', '').replace('<', ''))
            if n1 > n2:
                strand = -1
        except ValueError:
            pass
    is_partial = ('<' in text) or ('>' in text)
    low = text.lower()
    if low.startswith('join(') or low.startswith('order('):
        inner = text[text.find('(') + 1: text.rfind(')')]
        seg_texts = [s for s in re.split(r'[,\s]+', inner) if s]
    else:
        seg_texts = [text]
    segments = []
    for st in seg_texts:
        m2 = _SEGMENT_LOC_RE.fullmatch(st)
        if not m2:
            return None
        a, b = int(m2.group(1)), int(m2.group(2))
        segments.append((min(a, b), max(a, b)))
    if not segments:
        return None
    return segments, strand, is_partial


def validate_submission_gate(workdir, logger, tbl_path=None, fna_path=None,
                             src_path=None):
    """出 .sqn 前的硬校验 gate (S1-⑦)

    在 suvtk table2asn 运行前调用; 检查:
      1) featuretable.tbl 每个 CDS: 坐标可解析 (含负链反转区间 / complement(a..b)),
         长度 %3==0、坐标不越界 → ERROR; 起始/终止密码子非规范 → WARN
      2) qualifier 值 (gene/product/locus_tag/note 等) 含内嵌 tab/换行 → ERROR
      3) source.src 必填列 (_SRC_REQUIRED_COLS) 非空且不命中占位符 → ERROR

    任一 ERROR → 打印汇总并返回 False, 调用方应中止该步骤 (不标记成功);
    全部通过返回 True。--skip-gate 可在调用方显式跳过。
    """
    log = logger
    errors, warnings = [], []
    n_cds = 0

    tbl = _gate_locate(workdir, tbl_path,
                       ("featuretable_updated.tbl", "featuretable.tbl"),
                       ("analyze_hypothetical", "suvtk.features_output", "2_features"))
    fna = _gate_locate(workdir, fna_path, ("reoriented_nucleotide_sequences.fna",),
                       ("suvtk.features_output", "2_features"))
    src = _gate_locate(workdir, src_path, ("source.src",),
                       ("suvtk_submission", "submission", "metadata_output", "3_metadata"))

    log.info("─" * 40 + "\n  提交前硬校验 (gate)\n" + "─" * 40)
    if not tbl:
        errors.append("未找到 featuretable(.tbl): 无法校验 CDS 坐标/qualifier")
    if not fna:
        errors.append("未找到 .fna 核酸序列: 无法校验 CDS 长度/密码子")
    if not src:
        errors.append("未找到 source.src: 必填列无法校验")

    # 载入核酸序列 (校验坐标/密码子用)
    seqs = {}
    if fna:
        try:
            for rec in SeqIO.parse(str(fna), "fasta"):
                seqs[rec.id] = str(rec.seq).upper().replace('U', 'T')
        except Exception as e:                       # 文件解析失败按 ERROR 收集
            errors.append(f"{fna.name} 解析失败: {e}")

    # ── 1+2. featuretable.tbl ──
    if tbl:
        try:
            tbl_text = tbl.read_text(encoding='utf-8', errors='replace')
        except OSError as e:
            tbl_text = ''
            errors.append(f"{tbl.name} 读取失败: {e}")
        cur_id = None
        for i, raw in enumerate(tbl_text.split('\n'), 1):
            line = raw.rstrip('\r')
            if not line.strip():
                continue
            if line.startswith('>Feature '):
                ids = line.split('>Feature ', 1)[1].split()
                cur_id = ids[0] if ids else None
                if len(ids) > 1:
                    warnings.append(f"{tbl.name}:{i} >Feature 行含多个 ID {ids[:3]} "
                                    "→ 按第一个 ID 归属 (旧版语法)")
                continue
            stripped = line.strip()
            if line[0] in '\t ':
                # qualifier 行: 值含内嵌 tab/换行 → ERROR (S1-⑦②)
                kv = stripped.split('\t')
                if len(kv) > 2:
                    errors.append(f"{tbl.name}:{i} qualifier {kv[0]!r} 值含内嵌 tab/换行: "
                                  f"{stripped[:60]!r}")
                continue
            # 位置行: 仅 CDS 参与坐标校验
            parts = stripped.split()
            if not parts or parts[-1] != 'CDS':
                continue
            n_cds += 1
            loc_cols = parts[:-1]
            if len(loc_cols) == 1:
                parsed = _parse_cds_location(loc_cols[0])
            elif len(loc_cols) == 2:
                parsed = _parse_cds_location(loc_cols[0], loc_cols[1])
            else:
                parsed = None
            if parsed is None:
                errors.append(f"{tbl.name}:{i} CDS 位置无法解析: {stripped[:60]!r}")
                continue
            segments, strand, is_partial = parsed
            cds_len = sum(e - s + 1 for s, e in segments)
            if cds_len % 3 != 0:
                errors.append(f"{tbl.name}:{i} CDS 长度 {cds_len} 非 3 倍数 "
                              f"(位置 {stripped[:40]!r})")
            seq = seqs.get(cur_id) if cur_id else None
            if seq is None:
                errors.append(f"{tbl.name}:{i} 特征 ID {cur_id!r} 在 "
                              f"{fna.name if fna else '.fna'} 中不存在")
                continue
            lo = min(s for s, _ in segments)
            hi = max(e for _, e in segments)
            if lo < 1 or hi > len(seq):
                errors.append(f"{tbl.name}:{i} 坐标 {lo}..{hi} 超出 {cur_id} "
                              f"长度 {len(seq)}")
                continue
            # 起始/终止密码子 (单片段、非部分标记时检查; 异常仅 WARN)
            if len(segments) == 1 and not is_partial:
                s0, e0 = segments[0]
                sub = seq[s0 - 1:e0]
                if strand == -1:
                    sub = _revcomp(sub)
                if len(sub) >= 6:
                    if sub[:3] != 'ATG':
                        warnings.append(f"{tbl.name}:{i} {cur_id} 起始密码子 {sub[:3]} "
                                        "非 ATG")
                    if sub[-3:] not in ('TAA', 'TAG', 'TGA'):
                        warnings.append(f"{tbl.name}:{i} {cur_id} 终止密码子 {sub[-3:]} "
                                        "非 TAA/TAG/TGA")

    # ── 3. source.src 必填列 ──
    if src:
        ph_re = re.compile(_placeholder_pattern())
        try:
            with open(src, encoding='utf-8-sig') as vf:
                header = vf.readline().rstrip('\n').rstrip('\r').split('\t')
                missing_cols = [c for c in _SRC_REQUIRED_COLS if c not in header]
                if missing_cols:
                    errors.append(f"{src.name} 缺少必填列: {', '.join(missing_cols)}")
                n_rows, n_bad, bad_hits = 0, 0, []
                for ln, vline in enumerate(vf, 2):
                    if not vline.strip():
                        continue
                    n_rows += 1
                    vmap = dict(zip(header, vline.rstrip('\n').rstrip('\r').split('\t')))
                    for c in _SRC_REQUIRED_COLS:
                        if c not in header:
                            continue
                        val = str(vmap.get(c, '')).strip()
                        if not val or ph_re.search(val):
                            n_bad += 1
                            if len(bad_hits) < 5:
                                bad_hits.append(f"L{ln} {c}={val!r}")
                            break
                if n_bad:
                    errors.append(f"{src.name} 必填列空值/占位符: {n_bad}/{n_rows} 行 "
                                  f"({'; '.join(bad_hits)}) — 提交前必须补齐")
                else:
                    log.info("  source.src 必填列校验通过 (%d 行)", n_rows)
        except OSError as e:
            errors.append(f"{src.name} 读取失败: {e}")

    # ── 汇总 ──
    for w in warnings[:10]:
        log.warning("  [WARN] %s", w)
    if len(warnings) > 10:
        log.warning("  ... 另有 %d 条 WARN", len(warnings) - 10)
    if errors:
        log.error("  提交前硬校验未通过: %d ERROR, %d WARN (CDS %d 个)",
                  len(errors), len(warnings), n_cds)
        for e in errors[:20]:
            log.error("  [ERROR] %s", e)
        if len(errors) > 20:
            log.error("  ... 另有 %d 条 ERROR", len(errors) - 20)
        log.error("  修复后重跑; 确需跳过请显式加 --skip-gate")
        return False
    log.info("  ✓ 提交前硬校验通过 (CDS %d 个, WARN %d 条)", n_cds, len(warnings))
    return True


# ══════════════════════════════════════════════════════════════
# 自动检测 + 补运行缺失步骤
# ══════════════════════════════════════════════════════════════

def ensure_taxonomy(work_dir, fasta, suvtk_db, threads, log):
    """确保 suvtk taxonomy 已完成"""
    tax_out = work_dir / "suvtk.taxonomy_output"
    tax_tsv = tax_out / "taxonomy.tsv"
    # 检查是否已有 suvtk 输出 (可能在 09_Virome_Analysis)
    existing = _find_suvtk_output(work_dir, "suvtk.taxonomy_output/taxonomy.tsv")
    if existing and existing != tax_tsv:
        log.info("[1/5] taxonomy — 在 %s 找到, 复用", existing.parent)
        os.makedirs(tax_out, exist_ok=True)
        for fname in ["taxonomy.tsv", "miuvig_taxonomy.tsv", "miuvig_taxonomy_enriched.tsv",
                      "segmented_viruses_info.tsv", "taxresults_lca.tsv"]:
            src = existing.parent / fname
            dst = tax_out / fname
            if src.exists():
                if dst.exists() or dst.is_symlink():
                    dst.unlink()
                os.symlink(src, dst)
        return tax_out
    if tax_tsv.exists() and os.path.getsize(tax_tsv) > 100:
        log.info("[1/5] taxonomy — 已存在, 跳过")
        return tax_out

    log.info("[1/5] 运行 suvtk taxonomy ...")
    run(
        f"suvtk taxonomy -i {fasta} -o {tax_out} -d {suvtk_db} -s 0.7 -t {threads}",
        log, "suvtk taxonomy"
    )
    return tax_out if tax_tsv.exists() else None


def ensure_features(work_dir, fasta, tax_dir, suvtk_db, threads, log):
    """确保 suvtk features 已完成"""
    feat_out = work_dir / "suvtk.features_output"
    tbl = feat_out / "featuretable.tbl"
    existing = _find_suvtk_output(work_dir, "suvtk.features_output/featuretable.tbl")
    if existing and existing != tbl:
        log.info("[2/5] features — 在 %s 找到, 复用", existing.parent)
        os.makedirs(feat_out, exist_ok=True)
        # 链核心文件
        for fname in ["featuretable.tbl", "proteins.faa", "proteins_updated.faa",
                      "reoriented_nucleotide_sequences.fna", "miuvig_features.tsv",
                      "alignment.m8", "tophit_info.tsv"]:
            src = existing.parent / fname
            dst = feat_out / fname
            if src.exists():
                if dst.exists() or dst.is_symlink():
                    dst.unlink()
                os.symlink(src, dst)
        return feat_out
    if tbl.exists() and os.path.getsize(tbl) > 100:
        log.info("[2/5] features — 已存在, 跳过")
        return feat_out

    tax_tsv = tax_dir / "taxonomy.tsv" if tax_dir else None
    cmd = f"suvtk features -i {fasta} -o {feat_out} -d {suvtk_db} --coding-complete -t {threads}"
    if tax_tsv and tax_tsv.exists():
        cmd += f" --taxonomy {tax_tsv}"

    log.info("[2/5] 运行 suvtk features ...")
    run(cmd, log, "suvtk features")
    return feat_out if tbl.exists() else None


def ensure_hypothetical(work_dir, feat_dir, rvdb_dir, threads, log):
    """确保 analyze_hypothetical 已完成"""
    hypo_out = work_dir / "analyze_hypothetical"
    updated_tbl = hypo_out / "featuretable_updated.tbl"
    existing = _find_suvtk_output(work_dir, "analyze_hypothetical/featuretable_updated.tbl")
    if existing and existing != updated_tbl:
        log.info("[2.5/5] hypothetical — 在 %s 找到, 复用", existing.parent)
        os.makedirs(hypo_out, exist_ok=True)
        os.symlink(existing, updated_tbl)
        return hypo_out
    if updated_tbl.exists() and os.path.getsize(updated_tbl) > 100:
        log.info("[2.5/5] hypothetical — 已存在, 跳过")
        return hypo_out

    tbl = feat_dir / "featuretable.tbl"
    faa = feat_dir / "proteins.faa"
    if not tbl.exists() or not faa.exists():
        log.warning("[2.5/5] hypothetical — 缺少 input, 跳过")
        return None

    # 自动检测 RVDB 数据库
    rvdb_fasta = Path(rvdb_dir) / "U-RVDBv31.0-prot.EX.acc.fasta" if rvdb_dir else None
    rvdb_hmm = Path(rvdb_dir) / "U-RVDBv31.0-prot.hmm" if rvdb_dir else None
    rvdb_annot = Path(rvdb_dir) / "U-RVDBv31.0-prot.info.tab" if rvdb_dir else None

    script = SCRIPT_DIR / "analyze_hypothetical.py"

    cmd = (
        f"python {script} "
        f"-t {tbl} -f {faa} "
        f"-o {hypo_out} "
        f"--blast {hypo_out}/merged_blast.txt "
        f"--diamond -d {rvdb_fasta} "
        f"--hmmer --hmmer-db {rvdb_hmm} "
        f"--annot {rvdb_annot} "
        f"--threads {threads}"
    )

    log.info("[2.5/5] 运行 analyze_hypothetical (diamond + HMM) ...")
    run(cmd, log, "analyze_hypothetical")
    return hypo_out if updated_tbl.exists() else None


# ══════════════════════════════════════════════════════════════
# 已知病毒管道 (INFO → suvtk)
# ══════════════════════════════════════════════════════════════

def gather_metadata_from_summary(all_summary_path, output_dir, metadata_dir=None, log=None):
    """INFO 步骤: 从 all_summary.tsv 提取所有 SRA/CRR, 运行 gsa_sra.info.py 生成元数据。

    流程:
      1. 从 all_summary.tsv 提取唯一 SRA/CRR 编号
      2. 运行 gsa_sra.info.py -i <id_list> -o <output_dir>/0_info/
         → 生成 Global_Unified_Metadata_Full.tsv / Core14.tsv
      3. 合并额外的 Corrected_Metadata.tsv 和 runinfo.csv (如存在)
      4. 返回合并后的元数据 lookup

    Returns: dict {SRA_run: {collection_date, geo_loc_name, bioproject,
                              biosample, platform, host, tissue, source, center}}
    """
    import pandas as pd

    # 1. 提取唯一 SRA/CRR
    summary = pd.read_csv(all_summary_path, sep='\t')
    all_runs = set()
    for _, row in summary.iterrows():
        sra = str(row.get('Sample', '')).strip()
        if sra and sra.lower() not in ('nan', ''):
            all_runs.add(sra)

    if log:
        log.info("[INFO] 从 all_summary 提取 %d 个唯一 SRA/CRR", len(all_runs))

    # 2. 运行 gsa_sra.info.py 生成元数据
    info_out = Path(output_dir) / '0_info'
    info_out.mkdir(parents=True, exist_ok=True)
    full_tsv = info_out / 'Global_Unified_Metadata_Full.tsv'

    if not (full_tsv.exists() and os.path.getsize(full_tsv) > 1000):
        id_file = info_out / 'run_ids.txt'
        with open(id_file, 'w') as f:
            for rid in sorted(all_runs):
                f.write(rid + '\n')

        gsa_script = SCRIPT_DIR.parent / 'public_metadata_pipeline' / 'gsa_sra.info.py'
        if not gsa_script.exists():
            gsa_script = SCRIPT_DIR.parent.parent / 'public_metadata_pipeline' / 'gsa_sra.info.py'

        if gsa_script.exists():
            if log:
                log.info("  运行 gsa_sra.info.py — 获取 %d 个 Run 的元数据 ...", len(all_runs))
            cmd = f'python {gsa_script} -i {id_file} -o {info_out} -m local'
            run(cmd, log, 'gsa_sra.info')
        else:
            log.warning("  未找到 gsa_sra.info.py, 跳过在线元数据获取")

    # 3. 读取生成的元数据
    meta_lookup = {}

    def _ingest(df, source_name, log):
        added = 0
        for _, row in df.iterrows():
            run = str(row.get('Run', '')).strip()
            if not run or run not in all_runs:
                continue
            if run not in meta_lookup:
                meta_lookup[run] = {}

            for out_key, col_name in [
                ('collection_date', 'CollectionDate'),
                ('geo_loc_name', 'Location'),
                ('bioproject', 'BioProject'),
                ('biosample', 'BioSample'),
                ('platform', 'Platform'),
                ('host', 'ScientificName'),
                ('tissue', 'Tissue'),
                ('source', 'Source'),
                ('center', 'CenterName'),
            ]:
                if col_name in df.columns:
                    v = str(row.get(col_name, '')).strip()
                    if v and v.lower() not in ('nan', 'not_provided', '', ' '):
                        meta_lookup[run][out_key] = v
                        added += 1
        if log:
            log.info("  %s: 贡献 %d 字段", source_name, added)

    # Load Full.tsv (from gsa_sra.info.py output)
    if full_tsv.exists():
        df_full = pd.read_csv(full_tsv, sep=None, engine='python')
        df_full.columns = [c.lstrip('\ufeff') for c in df_full.columns]
        _ingest(df_full, 'Full', log)

    # Core14
    core14 = info_out / 'Global_Unified_Metadata_Core14.tsv'
    if core14.exists():
        df_core = pd.read_csv(core14, sep=None, engine='python')
        df_core.columns = [c.lstrip('\ufeff') for c in df_core.columns]
        _ingest(df_core, 'Core14', log)

    # 4. 额外合并 Corrected_Metadata.tsv 和 runinfo.csv (优先覆盖)
    if metadata_dir and os.path.isdir(metadata_dir):
        corrected = os.path.join(metadata_dir, 'Corrected_Metadata.tsv')
        if os.path.exists(corrected):
            df_cm = pd.read_csv(corrected, sep='\t')
            _ingest(df_cm, 'Corrected_Metadata', log)

        runinfo = os.path.join(metadata_dir, 'runinfo.csv')
        if os.path.exists(runinfo):
            df_ri = pd.read_csv(runinfo, sep=',')
            _ingest(df_ri, 'runinfo', log)

    if log:
        n_with_date = sum(1 for v in meta_lookup.values() if v.get('collection_date'))
        n_with_bs = sum(1 for v in meta_lookup.values() if v.get('biosample'))
        n_with_plat = sum(1 for v in meta_lookup.values() if v.get('platform'))
        log.info("[INFO] 元数据合并: %d/%d 个 Run, 含日期=%d, BioSample=%d, Platform=%d",
                 len(meta_lookup), len(all_runs), n_with_date, n_with_bs, n_with_plat)

    return meta_lookup


def build_taxonomy_from_summary(all_summary_path, fasta_dir=None, log=None):
    """从 all_summary.tsv 构建 contig-level taxonomy, 按 FASTA 过滤。

    Returns:
        tax_tsv_path, n_seqs, n_viruses
    """
    import pandas as pd

    df = pd.read_csv(all_summary_path, sep='\t')
    contig_tax = {}
    for _, row in df.iterrows():
        sample = str(row.get('Sample', '')).strip()
        acc = str(row.get('Accession', '')).strip()
        taxonomy = str(row.get('Taxonomy', '')).strip()
        if not sample or not acc:
            continue
        contig = f'{sample}_{acc}'
        if taxonomy and contig not in contig_tax:
            contig_tax[contig] = taxonomy

    # Filter by FASTA
    if fasta_dir and os.path.isdir(fasta_dir):
        fa_contigs = set()
        from pathlib import Path
        for fpath in Path(fasta_dir).rglob('*.fasta'):
            for line in open(fpath):
                if line.startswith('>'):
                    fa_contigs.add(line[1:].strip())
        before = len(contig_tax)
        contig_tax = {k: v for k, v in contig_tax.items() if k in fa_contigs}
        if log:
            log.info("[INFO] FASTA 过滤: %d → %d contigs (匹配 %d FASTA 头)",
                     before, len(contig_tax), len(fa_contigs))

    # Write taxonomy.tsv
    out_dir = Path(all_summary_path).parent.parent  # go up from summary/
    tax_path = out_dir / 'taxonomy.tsv'
    with open(tax_path, 'w') as f:
        f.write('contig\ttaxonomy\n')
        for contig, tax in sorted(contig_tax.items()):
            f.write(f'{contig}\t{tax}\n')

    return str(tax_path), len(contig_tax)


def build_inputs_from_extraction(extraction_dir, out_dir, log=None):
    """从 ④ 的 06_Extraction/ 全长提取产物构建 known 模式输入 (推荐口径)。

    布局 (B4 边界): <extraction_dir>/<病毒物种_登录号>/<样本>.<登录号>.full.fasta
      例: 06_Extraction/Cytorhabdovirus_sp_lycii_OR489165.1/GS-3.OR489165.1.full.fasta
          头: >GS-3_OR489165.1

    taxonomy.tsv 的物种名由病毒目录名解析 (下划线名 → 空格名, `_sp_` → ` sp. `,
    尾部登录号剥离); 后续 suvtk taxonomy 的 LCA 会再校准。

    Returns:
        (tax_tsv_path, combined_fasta_path, n_seqs, n_viruses)
    """
    edir = Path(extraction_dir)
    if not edir.is_dir():
        raise SystemExit(f'[!] extraction 目录不存在: {edir}')
    fasta_files = sorted(edir.glob('*/*.full.fasta')) or sorted(edir.rglob('*.full.fasta'))
    if not fasta_files:
        raise SystemExit(f'[!] {edir} 下未找到 */*.full.fasta 全长产物')

    # 目录名 → 物种名: 剥离尾部登录号 (NC_002030.1 / OR489165.1 等), _sp_ 还原为 sp.
    _ACC_RE = re.compile(r'^(?P<species>.+)_(?P<acc>[A-Za-z]{1,2}_?\d+(?:\.\d+)?)$')

    def dir_to_taxonomy(dirname):
        m = _ACC_RE.match(dirname)
        species = m.group('species') if m else dirname
        species = species.replace('_sp_', ' sp. ')
        return species.replace('_', ' ')

    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    combined_fa = out / 'combined_known_viruses.fasta'
    tax_path = out / 'taxonomy.tsv'

    n_seqs = 0
    virus_seen = set()
    with open(combined_fa, 'w', encoding='utf-8') as cf, \
            open(tax_path, 'w', encoding='utf-8') as tf:
        tf.write('contig\ttaxonomy\n')
        for fpath in fasta_files:
            virus_tax = dir_to_taxonomy(fpath.parent.name)
            virus_seen.add(fpath.parent.name)
            with open(fpath, encoding='utf-8') as inf:
                cf.write(inf.read())
            for line in open(fpath, encoding='utf-8'):
                if line.startswith('>'):
                    contig = line[1:].strip().split()[0]
                    tf.write(f'{contig}\t{virus_tax}\n')
                    n_seqs += 1
    if log:
        log.info("  全长提取产物: %d 条序列 / %d 个病毒参考 → %s",
                 n_seqs, len(virus_seen), combined_fa)
    return str(tax_path), str(combined_fa), n_seqs, len(virus_seen)


def prepare_known_viruses(all_summary_path, fasta_dir, metadata_dir, output_dir,
                          run_title, suvtk_db, threads, sequencer, assembler,
                          authors, pipeline, enrichment, source_uvig, metagenome_source,
                          log, extraction_dir=None):
    """已知病毒全流程: INFO → suvtk taxonomy → suvtk features → unified_metadata → 导出。

    参数:
      all_summary_path : 2_Virus_variants_Results/summary/all_summary.tsv
                         (与 extraction_dir 二选一; extraction 模式优先)
      extraction_dir   : ④ 的 06_Extraction/ 全长提取产物 (<病毒>/<样本>.full.fasta);
                         给出时序列直接取全长产物, 不需要 all_summary/fasta_dir
      fasta_dir        : 4_assemblies_clean/
      metadata_dir     : 6.Virus_variants_Results/summary/ (含 Full.tsv, runinfo.csv 等)
      output_dir       : 提交输出目录
    """
    import pandas as pd
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    log.info("=" * 60)
    log.info("已知病毒 GenBank 提交准备 (Known Virus Pipeline)")
    log.info("  all_summary: %s", all_summary_path)
    log.info("  fasta_dir:   %s", fasta_dir)
    log.info("  output:      %s", out)
    log.info("=" * 60)

    if extraction_dir:
        # ── 06_Extraction 全长模式: 序列取全长产物, 无 all_summary/Run 对照 ──
        log.info("[0/4] INFO — extraction 模式: 无 all_summary, 元数据走占位符/后续人工补")
        meta_lookup = {}
        log.info("[1/4] 从 06_Extraction 构建 taxonomy + 合并全长 FASTA ...")
        tax_path, combined_fa, n_contigs, n_viruses = build_inputs_from_extraction(
            extraction_dir, str(out), log)
        log.info("  合并全长 FASTA: %d 条序列 (%d 个病毒参考)", n_contigs, n_viruses)
    else:
        # ── all_summary 模式 (旧口径) ──
        # ── Step 0: INFO — gather metadata ──
        log.info("[0/4] INFO — 收集元数据 ...")
        meta_lookup = gather_metadata_from_summary(
            all_summary_path, str(out), metadata_dir, log)

        # ── Step 1: Build taxonomy + combine FASTA ──
        log.info("[1/4] 构建 taxonomy + 合并 FASTA ...")
        tax_path, n_contigs = build_taxonomy_from_summary(all_summary_path, fasta_dir, log)

        # Combine FASTA
        combined_fa = out / 'combined_known_viruses.fasta'
        with open(combined_fa, 'w') as cf:
            for fpath in Path(fasta_dir).rglob('*.fasta'):
                with open(fpath) as inf:
                    cf.write(inf.read())
        n_fa = sum(1 for _ in open(combined_fa) if _.startswith('>'))
        log.info("  合并 FASTA: %d 条序列 → %s", n_fa, combined_fa)

    # ── Step 2: suvtk taxonomy ──
    tax_dir = feat_dir = None
    if suvtk_db:
        log.info("[2/4] suvtk taxonomy ...")
        tax_dir = out / '1_taxonomy'
        tax_dir.mkdir(exist_ok=True)
        tax_out_tsv = tax_dir / 'taxonomy.tsv'
        if not tax_out_tsv.exists() or os.path.getsize(tax_out_tsv) < 100:
            cmd = (f'suvtk taxonomy -i {combined_fa} -o {tax_dir} '
                   f'-d {suvtk_db} -s 0.7 -t {threads}')
            run(cmd, log, 'suvtk taxonomy')
        else:
            log.info("  taxonomy 已存在, 跳过")
    else:
        log.info("[2/4] suvtk taxonomy — 跳过 (未指定 --suvtk-db)")

    # ── Step 3: suvtk features ──
    if suvtk_db:
        log.info("[3/4] suvtk features ...")
        feat_dir = out / '2_features'
        feat_dir.mkdir(exist_ok=True)
        tbl = feat_dir / 'featuretable.tbl'
        if not tbl.exists() or os.path.getsize(tbl) < 100:
            tax_tsv = out / '1_taxonomy' / 'taxonomy.tsv'
            cmd = (f'suvtk features -i {combined_fa} -o {feat_dir} '
                   f'-d {suvtk_db} --coding-complete -t {threads}')
            if tax_tsv.exists():
                cmd += f' --taxonomy {tax_tsv}'
            run(cmd, log, 'suvtk features')
        else:
            log.info("  features 已存在, 跳过")
    else:
        log.info("[3/4] suvtk features — 跳过 (未指定 --suvtk-db)")

    # ── Step 4: unified_metadata + 导出 ──
    log.info("[4/4] 生成统一元数据 + 导出 ...")
    meta_dir = out / 'submission'
    meta_dir.mkdir(exist_ok=True)

    meta_script = SCRIPT_DIR / 'unified_metadata.py'

    # Use the contig-taxonomy from summary (not suvtk's taxonomy)
    # suvtk taxonomy.tsv might have different format
    summary_tax = tax_path  # from build_taxonomy_from_summary

    cmd = (
        f'python {meta_script} '
        f'--taxonomy {summary_tax} '
        f'--run-title {run_title} '
        f'--assembler "{assembler}" '
        f'--sequencer "{sequencer}" '
        f'--enrichment "{enrichment}" '
        f'--source-uvig "{source_uvig}" '
        f'--metagenome-source "{metagenome_source}" '
        f'--pipeline "{pipeline}" '
        f'--authors "{authors}" '
        f'-o {meta_dir}/'
    )
    if metadata_dir:
        # Point to metadata dir so load_metadata can auto-detect Full.tsv
        pass
    run(cmd, log, 'unified_metadata')

    # ── Enrich CSV with metadata + abundance ──
    csv_path = meta_dir / 'unified_metadata.csv'
    if csv_path.exists():
        log.info("  充实 CSV 元数据 ...")
        df = pd.read_csv(csv_path, dtype=str, keep_default_na=False)
        enriched = 0

        # Build abundance lookup from all_summary.tsv (contig → {MeanDepth, Covered%, ...})
        abu_lookup = {}
        if all_summary_path and os.path.exists(all_summary_path):
            from_summary = pd.read_csv(all_summary_path, sep='\t')
            for _, r2 in from_summary.iterrows():
                s = str(r2.get('Sample', '')).strip()
                a = str(r2.get('Accession', '')).strip()
                if s and a:
                    contig_key = f'{s}_{a}'
                    abu_lookup[contig_key] = {
                        'depth': str(r2.get('Rep_MeanDepth', '')).strip(),
                        'covered': str(r2.get('Covered%', '')).strip(),
                        'reads': str(r2.get('Reads', '')).strip(),
                        'rel_abund': str(r2.get('Asm_Rel_Abund(%)', '')).strip(),
                    }
        elif extraction_dir:
            # extraction 模式: 丰度信息从 06_extraction 的 assembly_stats.csv 取 (若在)
            _stats = Path(extraction_dir) / 'assembly_stats.csv'
            if _stats.exists():
                try:
                    from_summary = pd.read_csv(_stats, sep=None, engine='python')
                    log.info("  丰度参考: %s", _stats.name)
                except Exception as _e:
                    log.info("  assembly_stats.csv 读取失败 (%s), 跳过丰度回填", _e)
                    from_summary = None
                if from_summary is not None:
                    for _, r2 in from_summary.iterrows():
                        s = str(r2.get('Sample', '')).strip()
                        a = str(r2.get('Accession', '')).strip()
                        if s and a:
                            abu_lookup[f'{s}_{a}'] = {
                                'depth': str(r2.get('Rep_MeanDepth', r2.get('MeanDepth', ''))).strip(),
                                'covered': str(r2.get('Covered%', '')).strip(),
                                'reads': str(r2.get('Reads', '')).strip(),
                                'rel_abund': str(r2.get('Asm_Rel_Abund(%)', '')).strip(),
                            }

        for idx, row in df.iterrows():
            contig = str(row.get('sequence_name', '')).strip()
            sra = re.match(r'([SC]RR\d+)', contig)
            sra = sra.group(1) if sra else contig.split('_')[0] if '_' in contig else contig

            # ── Fill from gsa_sra.info.py metadata ──
            if sra in meta_lookup:
                m = meta_lookup[sra]
                for csv_col, meta_key in [
                    ('collection_date', 'collection_date'),
                    ('src-geo_loc_name', 'geo_loc_name'),
                    ('bioproject', 'bioproject'),
                    ('biosample', 'biosample'),
                    ('src-Host', 'host'),
                    ('src-Tissue_type', 'tissue'),
                    ('src-Cultivar', 'source'),
                    ('src-Collected_by', 'center'),
                ]:
                    if csv_col in df.columns and m.get(meta_key):
                        current = str(row.get(csv_col, '')).strip()
                        is_ph = (not current or current.lower() in ('nan', 'not_provided', '') or
                                 'XXXX' in current or 'YYYY' in current or
                                 'Country:Region' in current or 'PRJNAXXXX' in current or
                                 'SAMNXXXXXXXX' in current or 'XX.' in current)
                        if is_ph:
                            df.at[idx, csv_col] = m[meta_key]
                            enriched += 1
                # Platform → Sequencing Technology
                if m.get('platform') and 'cmt-Sequencing_Technology' in df.columns:
                    cur = str(row.get('cmt-Sequencing_Technology', '')).strip()
                    if cur == sequencer or not cur or 'Illumina NovaSeq' in cur:
                        df.at[idx, 'cmt-Sequencing_Technology'] = m['platform']
                        enriched += 1

            # ── Fill abundance from all_summary ──
            if contig in abu_lookup and 'cmt-Genome_Coverage' in df.columns:
                a = abu_lookup[contig]
                cur_cov = str(row.get('cmt-Genome_Coverage', '')).strip()
                if not cur_cov:
                    depth = a['depth']
                    covered = a['covered']
                    parts = []
                    if depth and depth not in ('nan', ''):
                        parts.append(f'{depth}x')
                    if covered and covered not in ('nan', ''):
                        parts.append(f'{covered}%')
                    if parts:
                        df.at[idx, 'cmt-Genome_Coverage'] = ' '.join(parts)
                        enriched += 1

        df.to_csv(csv_path, index=False, encoding='utf-8-sig')
        log.info("  充实 %d 个单元格", enriched)

    # ── Topology ──
    # feat_dir 可缺 (未给 --suvtk-db): 显式 None 判定, 不再 NameError/TypeError
    fna = (feat_dir / 'reoriented_nucleotide_sequences.fna') if feat_dir else None
    if fna and fna.exists():
        topo_script = SCRIPT_DIR / 'viral_topology.py'
        topo_out = meta_dir / 'topology.tsv'
        run(f'python {topo_script} --taxonomy {summary_tax} --fasta {fna} -o {topo_out}',
            log, 'viral_topology')

    # ── HTML report (argv 直连, 不经 shell; 存在的产物全部透传) ──
    if csv_path.exists():
        report_script = SCRIPT_DIR / 'report_html.py'
        report_out = meta_dir / 'report.html'
        argv = [sys.executable, str(report_script), '--csv', str(csv_path),
                '--run-title', run_title, '-o', str(report_out)]
        suvtk_tax = _first_existing(
            tax_dir / 'miuvig_taxonomy.tsv', tax_dir / 'taxonomy.tsv') if tax_dir else None
        suvtk_feat = _first_existing(
            feat_dir / 'featuretable.tbl', feat_dir / 'miuvig_features.tsv') if feat_dir else None
        for flag, p in (('--miuvig', _first_existing(meta_dir / 'miuvig.tsv')),
                        ('--assembly', _first_existing(meta_dir / 'assembly.tsv')),
                        ('--topology', _first_existing(meta_dir / 'topology.tsv')),
                        ('--suvtk-tax', suvtk_tax),
                        ('--suvtk-feat', suvtk_feat)):
            if p is not None:
                argv += [flag, str(p)]
        run_argv(argv, log, 'report_html')

    log.info("=" * 60)
    log.info("完成! 输出: %s", meta_dir)
    log.info("  序列: %d, 元数据覆盖: %d/%d", n_contigs, len(meta_lookup), n_contigs)
    log.info("=" * 60)

    return str(meta_dir)


# ══════════════════════════════════════════════════════════════
# 主流程
# ══════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(
        description="submission_pipeline.py — MMPV-RNA → NCBI 提交统一入口",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument('--io-layout', choices=['standard', 'legacy'], default=None,
                        help='I/O 布局 (缺省: MMPV_IO_LAYOUT 环境变量, 再缺省 standard; '
                             '2026-10-15 新增 — 原先只吃环境变量, 独立运行 legacy 项目会被静默按 standard 解析)')
    sub = parser.add_subparsers(dest='command', help='管道模式')

    # ── known 模式: 已知病毒 (INFO → suvtk) ──
    p_known = sub.add_parser('known', help='已知病毒提交 (INFO → suvtk)',
        epilog="""
示例 (06_Extraction 全长模式, 推荐 — 序列直接取 ④ 提取后的全长产物):
  python submission_pipeline.py known \\
      --extraction-dir $OUT/06_Extraction/ \\
      --suvtk-db ~/database/virus-db/suvtk_db/ \\
      -o genbank_submission/test_known/ \\
      -t 40

示例 (all_summary 模式, 旧口径):
  python submission_pipeline.py known \\
      --all-summary known_virus_pipeline/2_Virus_variants_Results/summary/all_summary.tsv \\
      --fasta-dir known_virus_pipeline/4_assemblies_clean/ \\
      --metadata-dir 6.Virus_variants_Results/summary/ \\
      --suvtk-db ~/database/virus-db/suvtk_db/ \\
      -o genbank_submission/test_known/ \\
      -t 40
""")
    p_known.add_argument('--all-summary',
                         help='all_summary.tsv 路径 (2_Virus_variants_Results/summary/); '
                              '与 --extraction-dir 二选一')
    p_known.add_argument('--extraction-dir',
                         help='④ 的 06_Extraction/ 全长提取产物目录 (<病毒>/<样本>.full.fasta); '
                              '给出时序列直接取全长产物, 不需要 --all-summary/--fasta-dir')
    p_known.add_argument('--fasta-dir',
                         help='4_assemblies_clean/ 目录 (all_summary 模式用; extraction 模式忽略)')
    p_known.add_argument('--metadata-dir',
                         help='元数据目录 (默认自动从 all_summary 推断: 同 data-* 父级下的 6.Virus_variants_Results/summary/)')
    p_known.add_argument('--suvtk-db', help='suvtk 数据库路径')
    p_known.add_argument('-o', '--output', default='./genbank_submission/known',
                         help='输出目录')
    p_known.add_argument('--run-title', default='known_virus_submission',
                         help='运行标题')
    p_known.add_argument('-t', '--threads', type=int, default=40)
    p_known.add_argument('--assembler', default='MEGAHIT;1.2.9;default parameters')
    p_known.add_argument('--sequencer', default='Illumina NovaSeq 6000')
    p_known.add_argument('--enrichment', default='none')
    p_known.add_argument('--pipeline', default='MMPV-RNA v2.3 + suvtk v0.1.1')
    p_known.add_argument('--source-uvig', default='metatranscriptome (not viral targeted)')
    p_known.add_argument('--metagenome-source', default='plant virome')
    p_known.add_argument('--authors', default='Zhang, Wenda')

    # ── 默认模式: discovery pipeline (08_Rescue) ──
    p_disc = sub.add_parser('discovery', help='新病毒提交 (discovery pipeline, 08_Rescue)',
        epilog="""
示例:
  python submission_pipeline.py discovery --work-dir $OUT/08_Rescue/ --suvtk-db ~/.../suvtk_db/ -t 40
""")
    p_disc.add_argument('--work-dir', required=True,
                        help='上游输出目录 ($OUT 或 $OUT/08_Rescue, 自动检测 9b 的 class_KEEP.fasta)')
    p_disc.add_argument('--fasta',
                        help='输入 FASTA (默认自动检测 9b 的 class_KEEP.fasta)')
    p_disc.add_argument('--suvtk-db', help='suvtk 数据库路径')
    p_disc.add_argument('--rvdb-db', help='RVDB 数据库目录')
    p_disc.add_argument('--run-title', default='viral_submission', help='运行标题')
    p_disc.add_argument('--mode', choices=['suvtk', 'sequin', 'both'], default='both')
    p_disc.add_argument('-t', '--threads', type=int, default=40)
    p_disc.add_argument('--fetch-metadata', action='store_true')
    p_disc.add_argument('--metadata-dir', help='Core14 元数据目录')
    p_disc.add_argument('--config', help='YAML 配置文件')
    p_disc.add_argument('--interactive', action='store_true')
    p_disc.add_argument('--pipeline-type', choices=['auto', 'discovery', 'analysis'],
                        default='auto')
    p_disc.add_argument('--ref-info', help='已知病毒参考信息')
    p_disc.add_argument('--force-taxonomy', action='store_true')
    p_disc.add_argument('--force-features', action='store_true')
    p_disc.add_argument('--force-hypo', action='store_true')
    p_disc.add_argument('--export-metadata', action='store_true')
    p_disc.add_argument('--skip-gate', action='store_true',
                        help='显式跳过出 .sqn 前的硬校验 gate (CDS 长度/密码子/qualifier/source.src)')

    args = parser.parse_args()

    # 布局解析 (2026-10-15): CLI > 环境变量 > standard。解析结果写回环境供子进程
    # 继承, 并重绑模块级 _D —— _D 原先只在 import 时按环境变量算一次。
    global _D
    _D = _layout_dirs(normalize_layout_env(getattr(args, 'io_layout', None)))

    # ── known 模式 ──
    if args.command == 'known':
        out = Path(args.output)
        log = setup_logger(str(out))

        if not args.all_summary and not args.extraction_dir:
            parser.error("known 模式需要 --extraction-dir (④ 06_Extraction 全长产物) 或 "
                         "--all-summary (2_Virus_variants_Results/summary/) 二选一")

        # Auto-detect metadata_dir from all_summary path (extraction 模式无 Run 对照, 跳过)
        metadata_dir = args.metadata_dir
        if metadata_dir is None and args.all_summary:
            summary_dir = Path(args.all_summary).resolve().parent  # .../2_Virus_variants_Results/summary
            # Go up 3 levels to data-{year}/sra_rna.dataN/
            for up_to in [summary_dir.parent.parent.parent,  # data-*/sra_rna.*
                          summary_dir.parent.parent.parent.parent]:
                candidate = up_to / '6.Virus_variants_Results' / 'summary'
                if candidate.is_dir():
                    metadata_dir = str(candidate)
                    log.info("  自动检测 metadata-dir: %s", metadata_dir)
                    break
            if not metadata_dir:
                log.warning("  未找到 6.Virus_variants_Results/summary/, 元数据将为空")

        prepare_known_viruses(
            all_summary_path=args.all_summary,
            fasta_dir=args.fasta_dir,
            extraction_dir=args.extraction_dir,
            metadata_dir=metadata_dir or '',
            output_dir=str(out),
            run_title=args.run_title,
            suvtk_db=args.suvtk_db or '',
            threads=args.threads,
            sequencer=args.sequencer,
            assembler=args.assembler,
            authors=args.authors,
            pipeline=args.pipeline,
            enrichment=args.enrichment,
            source_uvig=args.source_uvig,
            metagenome_source=args.metagenome_source,
            log=log,
        )
        return

    # ── discovery 模式 ──
    if args.command != 'discovery':
        parser.print_help()
        sys.exit(0)

    work_dir = Path(args.work_dir)
    fasta = args.fasta or _detect_keep_fasta(work_dir)

    if not fasta or not os.path.exists(fasta):
        print(f"[!] 未找到 KEEP FASTA 文件: {fasta or '(空)'}")
        print("    请用 --fasta 显式指定, 或先运行 discovery 的 9b (analysis_verify) 阶段")
        sys.exit(1)

    log = setup_logger(str(work_dir))

    # 提交未完成项登记 (2026-10-15): gate 未过 / .sqn 未生成 / 缺输入不再只打日志,
    # main 末尾聚合后非零退出 —— 旧行为可"gate 未过、无 .sqn"仍退出 0 (假成功)。
    submission_failures = []

    # ── 自动检测管道类型 ──
    pipeline_type = args.pipeline_type
    if pipeline_type == 'auto':
        # 检测 discovery pipeline 特征 (Plant/, centroids/)
        is_discovery = (work_dir / "Plant").exists() or (work_dir / "centroids").exists()
        # 检测 analysis pipeline 特征 (4_assemblies_clean/, known_viruses/)
        is_analysis = (work_dir.parent / "4_assemblies_clean").exists() or \
                      (work_dir / "4_assemblies_clean").exists() or \
                      (len(list(work_dir.glob("*_assemblies_clean"))) > 0)

        if is_discovery and not is_analysis:
            pipeline_type = 'discovery'
        elif is_analysis and not is_discovery:
            pipeline_type = 'analysis'
        elif is_discovery and is_analysis:
            pipeline_type = 'discovery'
        else:
            pipeline_type = 'discovery'

    # ── 收集病毒列表 ──
    virus_list = []  # [{name, accession, fastas: [path, ...], n_samples}]
    if pipeline_type == 'analysis':
        # 检测 analysis pipeline 输出: 4_assemblies_clean/
        analysis_dir = work_dir.parent / "4_assemblies_clean"
        if not analysis_dir.exists():
            analysis_dir = work_dir / "4_assemblies_clean"
        if not analysis_dir.exists():
            # 尝试从 work_dir 的兄弟目录找
            for parent in [work_dir.parent, work_dir.parent.parent]:
                candidate = parent / "4_assemblies_clean"
                if candidate.exists():
                    analysis_dir = candidate
                    break

        if analysis_dir.exists():
            for virus_dir in sorted(analysis_dir.iterdir()):
                if virus_dir.is_dir() and not virus_dir.name.startswith('.'):
                    fastas = sorted(virus_dir.glob("*.full.fasta")) + sorted(virus_dir.glob("*.fasta"))
                    # 从目录名解析病毒名和 accession
                    parts = virus_dir.name.rsplit('_', 1)
                    virus_name = parts[0] if len(parts) > 1 else virus_dir.name
                    accession = parts[1] if len(parts) > 1 else ''
                    virus_list.append({
                        'name': virus_name,
                        'accession': accession,
                        'dir': str(virus_dir),
                        'fastas': [str(f) for f in fastas],
                        'n_samples': len(fastas),
                    })
            log.info("  检测到 analysis pipeline: 1 种病毒 × %d 个样本 (全长序列)", len(virus_list) if virus_list else 0)
            for v in virus_list:
                log.info("    %s (%s): %d samples", v['name'], v['accession'], v['n_samples'])
        else:
            log.warning("  未找到 4_assemblies_clean/ 目录")

    elif pipeline_type == 'discovery':
        log.info("  检测到 discovery pipeline: %s", fasta)

    log.info("=" * 60)
    log.info("MMPV-RNA → NCBI 提交流水线 v1.0")
    log.info("  工作目录: %s", work_dir)
    log.info("  管道类型: %s — %s",
             pipeline_type,
             f"1种病毒 × {virus_list[0]['n_samples'] if virus_list else '?'} 个样本" if pipeline_type == 'analysis'
             else "多病毒, 每病毒1条代表序列")
    log.info("  输出模式: %s", args.mode)
    log.info("=" * 60)

    # ═══ 步骤 1-2.5: 自动补缺失 ═══

    # 1. taxonomy
    tax_dir = None
    if args.suvtk_db:
        tax_dir = work_dir / "suvtk.taxonomy_output"
        if args.force_taxonomy and tax_dir.exists():
            import shutil
            shutil.rmtree(tax_dir)
        tax_dir = ensure_taxonomy(work_dir, fasta, args.suvtk_db, args.threads, log)

    # 2. features
    feat_dir = None
    if args.suvtk_db:
        feat_dir = work_dir / "suvtk.features_output"
        if args.force_features and feat_dir.exists():
            import shutil
            shutil.rmtree(feat_dir)
        feat_dir = ensure_features(work_dir, fasta, tax_dir, args.suvtk_db, args.threads, log)

    # 2.5. hypothetical (仅 discovery pipeline 需要)
    hypo_dir = None
    if pipeline_type == 'discovery' and feat_dir and args.rvdb_db:
        hypo_dir = work_dir / "analyze_hypothetical"
        if args.force_hypo and hypo_dir.exists():
            import shutil
            shutil.rmtree(hypo_dir)
        hypo_dir = ensure_hypothetical(work_dir, feat_dir, args.rvdb_db, args.threads, log)
    elif pipeline_type == 'analysis':
        log.info("[2.5/5] hypothetical — 已知病毒跳过 (analysis pipeline 使用参考基因组注释)")

    # 2.6. 公共元数据获取 (public_metadata_pipeline)
    meta_dir = args.metadata_dir or str(work_dir / "metadata_output")
    core14 = Path(meta_dir) / "Global_Unified_Metadata_Full.tsv"
    if not core14.exists():
        core14 = Path(meta_dir) / "Global_Unified_Metadata_Core14.tsv"

    # taxonomy.tsv 由步骤 1 产出; 在本函数后段 (步骤 3) 才正式赋值, 此处提前解析
    tax_tsv_early = (tax_dir / "taxonomy.tsv") if tax_dir else None
    if args.fetch_metadata and tax_tsv_early and tax_tsv_early.exists():
        tax_tsv = tax_tsv_early
        log.info("[2.6/5] 获取公共元数据 (NCBI SRA → 自动填 source.src)")

        # 提取 SRA/CRR 编号
        srr_list = work_dir / "srr_list.txt"
        if not srr_list.exists() or args.force_hypo:
            log.info("  提取 SRA/CRR 编号 ...")
            import subprocess as _sp
            try:
                result = _sp.run([BASH, '-c',
                                  f"grep -oP '[SC]RR\\d+' {tax_tsv} | sort -u > {srr_list}"])
            except OSError as _e:
                # 非 POSIX 环境 (无 bash): 按失败处理, srr_list 缺失时下游已兜底
                log.warning("  SRA/CRR 编号提取失败 (shell 不可用): %s", _e)
            n_sra = sum(1 for _ in open(srr_list)) if srr_list.exists() else 0
            log.info("  找到 %d 个唯一 SRA/CRR", n_sra)

        if not core14.exists():
            log.info("  运行 gsa_sra.info.py 获取元数据 ...")
            gsa_script = work_dir.parent.parent.parent / "MMPV-RNA" / "public_metadata_pipeline" / "gsa_sra.info.py"
            if not gsa_script.exists():
                gsa_script = Path.home() / "MMPV-RNA" / "public_metadata_pipeline" / "gsa_sra.info.py"
            cmd = (
                f"python {gsa_script} "
                f"-i {srr_list} "
                f"-o {meta_dir} "
                f"-m local "
                f"-t {args.threads}"
            )
            run(cmd, log, "gsa_sra.info", check=False)

    if core14.exists():
        log.info("  元数据已就绪: %s", core14)
    else:
        log.info("  无公共元数据文件")
        if pipeline_type == 'analysis':
            log.info("  ═══════════════════════════════════════════════")
            log.info("  📋 新测序数据提交指南:")
            log.info("  1. unified_metadata.csv 已生成 (占位符待填)")
            log.info("  2. 去 NCBI 注册 BioProject → 拿到 PRJNA-xxx")
            log.info("     https://submit.ncbi.nlm.nih.gov/")
            log.info("  3. 用 biosample_template.tsv 批量注册 BioSample")
            log.info("  4. 编辑 unified_metadata.csv 填入:")
            log.info("     collection_date, geo_loc_name, lat_lon,")
            log.info("     bioproject (你的 PRJNA), biosample (你的 SAMN)")
            log.info("  5. 重新运行 --export-metadata 不带 --fetch-metadata")
            log.info("  ═══════════════════════════════════════════════")
        else:
            log.info("  公共数据: source.src 从 Core14 自动填, 样本信息已有")
        meta_dir = None

    # ═══ 步骤 3-5: 生成提交文件 ═══

    # 确定使用的 featuretable
    updated_tbl = None
    updated_faa = None
    if hypo_dir and (hypo_dir / "featuretable_updated.tbl").exists():
        updated_tbl = str(hypo_dir / "featuretable_updated.tbl")
        updated_faa = str(hypo_dir / "proteins_updated.faa")
    elif feat_dir:
        updated_tbl = str(feat_dir / "featuretable.tbl")
        updated_faa = str(feat_dir / "proteins.faa")

    tax_tsv = tax_dir / "taxonomy.tsv" if tax_dir else None
    miuvig_tax = tax_dir / "miuvig_taxonomy.tsv" if tax_dir else None
    fna = feat_dir / "reoriented_nucleotide_sequences.fna" if feat_dir else None
    miuvig_feat = feat_dir / "miuvig_features.tsv" if feat_dir else None

    # 统一输出目录 (work_dir 平级，不污染分析目录)
    sub_out = work_dir.parent / f"submission_{args.run_title}"
    sub_out.mkdir(parents=True, exist_ok=True)
    report_path = sub_out / "report.html"   # ⑦ 报告实际落点 (尾段判存在打印)

    if args.mode in ('suvtk', 'both'):
        log.info("─" * 40 + "\n  suvtk 模式: 生成 source.src + .cmt → table2asn\n" + "─" * 40)

        sv_dir = sub_out / "suvtk_submission"
        sv_dir.mkdir(parents=True, exist_ok=True)

        miuvig_tax = tax_dir / "miuvig_taxonomy.tsv" if tax_dir else None
        miuvig_feat = feat_dir / "miuvig_features.tsv" if feat_dir else None
        reoriented_fna = feat_dir / "reoriented_nucleotide_sequences.fna" if feat_dir else None

        # 1. 生成 miuvig.tsv (MIUVIG global parameters, key-value)
        #    ref: https://github.com/NCBI-Codeathons/suvtk
        miuvig_tsv = sv_dir / "miuvig.tsv"
        if not miuvig_tsv.exists():
            pipe_params = load_pipeline_params(work_dir)
            with open(miuvig_tsv, "w") as mf:
                mf.write("MIUVIG_parameter\tvalue\n")
                mf.write(f"source_uvig\t{pipe_params.get('source_uvig','metatranscriptome (not viral targeted)')}\n")
                mf.write(f"assembly_software\t{pipe_params.get('assembler','MEGAHIT v1.2.9')}\n")
                mf.write("vir_ident_software\tGenomad;1.7.0;default parameters\n")
                mf.write(f"nucl_acid_ext\t{pipe_params.get('enrichment','RNA extracted from leaf tissue')}\n")
            log.info("  生成 miuvig.tsv")

        # 2. 生成 assembly.tsv (Assembly-Data structured comment, key-value)
        asm_tsv = sv_dir / "assembly.tsv"
        if not asm_tsv.exists():
            pipe_params = load_pipeline_params(work_dir)
            with open(asm_tsv, "w") as af:
                af.write("Assembly_parameter\tvalue\n")
                af.write("StructuredCommentPrefix\tAssembly-Data\n")
                af.write(f"Assembly Method\t{pipe_params.get('assembler','MEGAHIT v1.2.9')}\n")
                af.write(f"Sequencing Technology\t{pipe_params.get('sequencer','Illumina NovaSeq 6000')}\n")
            log.info("  生成 assembly.tsv")

        # 2.8 tbl2gb: 生成 GenBank .gb 文件 (供 SnpEff/capheine/similarity 下游)
        gb_dir = sv_dir / "virus-annotations"
        gb_dir.mkdir(exist_ok=True)
        tbl2gb_script = SCRIPT_DIR.parent / "virome_analysis_pipeline" / "utils" / "tbl2gb.py"
        tax_tsv_path = tax_dir / "taxonomy.tsv" if tax_dir else None
        if tbl2gb_script.exists() and reoriented_fna and updated_tbl and tax_tsv_path and tax_tsv_path.exists():
            if not list(gb_dir.glob("*.gb")):
                tbl2gb_cmd = (
                    f"python {tbl2gb_script} "
                    f"--tbl {updated_tbl} "
                    f"--fasta {reoriented_fna} "
                    f"--taxonomy {tax_tsv_path} "
                    f"-o {gb_dir}"
                )
                run(tbl2gb_cmd, log, "tbl2gb → .gb")
            else:
                log.info("[2.8/5] tbl2gb — .gb 文件已存在, 跳过")
        else:
            missing = []
            if not tbl2gb_script.exists(): missing.append("tbl2gb.py")
            if not reoriented_fna or not reoriented_fna.exists(): missing.append("reoriented.fna")
            if not updated_tbl or not os.path.exists(updated_tbl): missing.append("featuretable.tbl")
            if missing:
                log.info("[2.8/5] tbl2gb — 跳过 (%s)", ", ".join(missing))

        # 3. suvtk comments: 生成 MIUVIG 结构化注释 .cmt (suvtk 会自动加 .cmt 后缀)
        comments_cmt = sv_dir / "comments"
        if miuvig_tax and miuvig_tax.exists() and miuvig_feat and miuvig_feat.exists():
            comments_cmd = (
                f"suvtk comments "
                f"-t {miuvig_tax} "
                f"-f {miuvig_feat} "
                f"-m {miuvig_tsv} "
                f"-a {asm_tsv} "
                f"-o {comments_cmt}"
            )
            run(comments_cmd, log, "suvtk comments")

        # 4. suvtk table2asn → .sqn
        src_file = sv_dir / "source.src"
        sbt_file = sv_dir / "template.sbt"

        # 生成 source.src (TSV 格式, suvtk table2asn 要求)
        # ref: suvtk docs — Sequence_ID, Organism, Isolate, Collection_date, geo_loc_name, Lat_Lon, ...
        # NOTE: Collection_date/geo_loc_name/Lat_Lon 不写死默认值 (旧版本硬编码 2024 /
        #       China:Ningxia, 把推测写成观测); 必须由用户补齐后重新生成。
        if src_file.exists():
            log.info("  复用已存在的 source.src: %s (如需重新生成请先删除)", src_file)
        else:
            # 从 taxonomy 和 ref_info 读取实际数据
            # ref_info 定位 (2026-10-15 修): ref_info.tsv 在 ③ 输出根的 d_analysis
            # (standard=11_ViromeAnalysis / legacy=09a_Virome_Analysis) 下 —— 按文档
            # 姿势 --work-dir <③根> 时应查 work_dir 本身; 旧的 work_dir.parent 只在
            # 把 work-dir 指到 rescue 子目录时才命中, 两姿势互斥且失败时静默落
            # 'unclassified plant virus' 占位。现两布局 × 两级全探测, 未命中告警。
            _da_names = []
            for _da in (_D["d_analysis"],
                        _layout_dirs("standard")["d_analysis"],
                        _layout_dirs("legacy")["d_analysis"]):
                if _da not in _da_names:
                    _da_names.append(_da)
            _ref_info_cands = [base / da / "ref_info.tsv"
                               for da in _da_names
                               for base in (work_dir, work_dir.parent)]
            ref_info_path = next((c for c in _ref_info_cands if c.is_file()),
                                 _ref_info_cands[0])
            has_ref_info = ref_info_path.is_file()
            if not has_ref_info:
                log.warning(
                    "  未找到 ref_info.tsv (已探测: %s) — source.src 的 Organism 将落 "
                    "'unclassified plant virus' 占位; 如需真实物种名, 请把 --work-dir "
                    "指向 ③ 输出根或手工放置 ref_info.tsv",
                    "; ".join(str(c) for c in _ref_info_cands))
            # 输入可能是 9b 的 class_KEEP.fasta (HQ 子集), 而 ref_info 覆盖 HQ 全集,
            # 只保留输入 fasta 中实际存在的序列 ID, 避免 source.src 出现多余行导致 table2asn 对不上
            fasta_ids = {rec.id for rec in SeqIO.parse(str(fasta), "fasta")}
            with open(src_file, "w") as sf:
                cols = ["Sequence_ID","Organism","Isolate","Collection_date","geo_loc_name",
                        "Lat_Lon","Bioproject","Biosample","SRA","Segment","Metagenomic","Metagenome_source","Host"]
                sf.write("\t".join(cols) + "\n")
                # Organism 用 suvtk taxonomy 结果; isolate 以 (物种, 样本) 为组:
                # 同物种同样本的多条 contig (分段病毒各片段) 共享同一 isolate —
                # NCBI/suvtk 按 Isolate 归组分段病毒 (旧版每序列一个 isolate,
                # 分段无法归组); 组哈希与 FASTA 行序无关, 跨数据集合并也不重名
                import hashlib

                def _isolate_for(organism, cid):
                    m = re.search(r'([SC]RR\d+)', cid)
                    sample = m.group(1) if m else cid.split('_')[0]
                    key = f'{organism}|{sample}'
                    if key not in _isolate_for._seen:
                        _isolate_for._seen[key] = 'NXGQ' + hashlib.md5(
                            key.encode('utf-8')).hexdigest()[:5].upper()
                    return _isolate_for._seen[key]

                _isolate_for._seen = {}
                n_rows = 0
                for row in _read_tsv_simple(ref_info_path) if has_ref_info else []:
                    cid = row.get("Accession", row.get("contig_id", ""))
                    if not cid or cid not in fasta_ids: continue
                    organism = row.get("suvtk_Species", row.get("Species", "unclassified plant virus"))
                    if organism in ("nan", "", "NA", None):
                        organism = "unclassified plant virus"
                    isolate = _isolate_for(organism, cid)
                    sf.write("\t".join([cid, organism, isolate,
                                        "", "", "", "", "", "", "", "TRUE",
                                        "plant virome", "Lycium chinense"]) + "\n")
                    n_rows += 1
                if not has_ref_info:
                    # fallback: generate from FASTA
                    for rec in SeqIO.parse(str(fasta), "fasta"):
                        isolate = _isolate_for("unclassified plant virus", rec.id)
                        sf.write("\t".join([rec.id, "unclassified plant virus", isolate,
                                            "", "", "", "", "", "", "", "TRUE",
                                            "plant virome", "Lycium chinense"]) + "\n")
                        n_rows += 1
            log.info("  生成 source.src (%d 条)", n_rows)

        # source.src 必填列校验: Collection_date/Lat_Lon/Bioproject/Biosample/SRA 为空时显式告警
        _required_src_cols = ("Collection_date", "geo_loc_name", "Lat_Lon",
                              "Bioproject", "Biosample", "SRA")
        try:
            with open(src_file) as vf:
                vh = vf.readline().rstrip("\n").split("\t")
                n_rows, n_bad = 0, 0
                empty_cols = set()
                for vline in vf:
                    if not vline.strip():
                        continue
                    n_rows += 1
                    vals = vline.rstrip("\n").split("\t")
                    vmap = dict(zip(vh, vals))
                    for c in _required_src_cols:
                        if not vmap.get(c, "").strip():
                            empty_cols.add(c)
                    if any(not vmap.get(c, "").strip() for c in _required_src_cols):
                        n_bad += 1
                if empty_cols:
                    log.warning("  ⚠️ source.src 必填列未填: %s (%d/%d 行不完整) — "
                                "提交前必须人工补齐, 否则 table2asn 会失败",
                                ", ".join(sorted(empty_cols)), n_bad, n_rows)
                else:
                    log.info("  source.src 必填列校验通过 (%d 行)", n_rows)
        except OSError as e:
            log.warning("  source.src 校验读取失败: %s", e)

        # 生成 template.sbt (NCBI ASN.1 格式, table2asn 兼容)
        # NOTE: ASN.1 解析对格式极其敏感, 每个字段必须单独一行
        # 联系邮箱来自 NCBI_EMAIL 环境变量, 取不到留空并告警 (不硬编码邮箱)
        email = _ncbi_email(log)
        if not sbt_file.exists():
            with open(sbt_file, "w") as sf:
                sf.write(
                    'Submit-block ::= {\n'
                    '  contact {\n'
                    '    contact {\n'
                    '      name name {\n'
                    '        last "Zhang",\n'
                    '        first "Wenda",\n'
                    '        middle "",\n'
                    '        initials "",\n'
                    '        suffix "",\n'
                    '        title ""\n'
                    '      },\n'
                    '      affil std {\n'
                    '        affil "Ningxia University",\n'
                    '        div "College of Life Sciences",\n'
                    '        city "Yinchuan",\n'
                    '        sub "Ningxia",\n'
                    '        country "China",\n'
                    '        street "Ningxia University",\n'
                    '        email "' + email + '",\n'
                    '        postal-code "750021"\n'
                    '      }\n'
                    '    }\n'
                    '  },\n'
                    '  cit {\n'
                    '    authors {\n'
                    '      names std {\n'
                    '        {\n'
                    '          name name {\n'
                    '            last "Zhang",\n'
                    '            first "Wenda",\n'
                    '            middle "",\n'
                    '            initials "",\n'
                    '            suffix "",\n'
                    '            title ""\n'
                    '          }\n'
                    '        }\n'
                    '      },\n'
                    '      affil std {\n'
                    '        affil "Ningxia University",\n'
                    '        div "College of Life Sciences",\n'
                    '        city "Yinchuan",\n'
                    '        sub "Ningxia",\n'
                    '        country "China",\n'
                    '        street "Ningxia University",\n'
                    '        postal-code "750021"\n'
                    '      }\n'
                    '    }\n'
                    '  },\n'
                    '  subtype new\n'
                    '}\n'
                    'Seqdesc ::= pub {\n'
                    '  pub {\n'
                    '    gen {\n'
                    '      cit "unpublished",\n'
                    '      authors {\n'
                    '        names std {\n'
                    '          {\n'
                    '            name name {\n'
                    '              last "Zhang",\n'
                    '              first "Wenda",\n'
                    '              middle "",\n'
                    '              initials "",\n'
                    '              suffix "",\n'
                    '              title ""\n'
                    '            }\n'
                    '          }\n'
                    '        }\n'
                    '      },\n'
                    '      title "Plant virome of Lycium chinense in Ningxia, China"\n'
                    '    }\n'
                    '  }\n'
                    '}\n'
                    + (('Seqdesc ::= user {\n'
                        '  type str "Submission",\n'
                        '  data {\n'
                        '    {\n'
                        '      label str "AdditionalComment",\n'
                        '      data str "ALT EMAIL:' + email + '"\n'
                        '    }\n'
                        '  }\n'
                        '}\n') if email else '')
                    + 'Seqdesc ::= user {\n'
                    '  type str "Submission",\n'
                    '  data {\n'
                    '    {\n'
                    '      label str "AdditionalComment",\n'
                    '      data str "Submission Title:None"\n'
                    '    }\n'
                    '  }\n'
                    '}\n'
                )
            log.info("  生成 template.sbt")

        comments_file = Path(str(comments_cmt) + ".cmt")  # suvtk 自动加 .cmt

        # 4.5 出 .sqn 前硬校验 gate (S1-⑦): 任一 ERROR → 中止 table2asn (不标记成功)
        gate_ok = True
        if args.skip_gate:
            log.warning("  --skip-gate: 跳过提交前硬校验 gate")
        else:
            gate_ok = validate_submission_gate(
                work_dir, log,
                tbl_path=updated_tbl,
                fna_path=str(reoriented_fna) if reoriented_fna else None,
                src_path=str(src_file))

        # 5. suvtk table2asn → .sqn
        if not gate_ok:
            log.error("  提交前硬校验未通过 → 跳过 suvtk table2asn (确需跳过请加 --skip-gate)")
            submission_failures.append("提交前硬校验 gate 未通过 (已跳过 table2asn)")
        elif reoriented_fna and reoriented_fna.exists() and updated_tbl and os.path.exists(updated_tbl):
            if not comments_file.exists():
                log.warning("  comments.cmt 不存在, 跳过 suvtk table2asn")
                submission_failures.append("comments.cmt 缺失 (跳过 table2asn)")
            elif not sbt_file.exists():
                log.warning("  template.sbt 不存在, 跳过 suvtk table2asn")
                submission_failures.append("template.sbt 缺失 (跳过 table2asn)")
            else:
                # 直连 NCBI table2asn (argv 列表, 不经 shell — Windows 中文路径安全);
                # suvtk 0.1.x/0.2.x 不可用或直连失败时回退 suvtk 包装器
                sqn_prefix = sv_dir / "submission"
                produced = False
                try:
                    from virome_submission import _run_table2asn_chain as _t2chain
                    produced = _t2chain(reoriented_fna, updated_tbl, str(src_file),
                                        str(sbt_file), str(comments_file),
                                        str(sqn_prefix), log)
                except Exception as e:
                    log.warning("  直连 table2asn 链异常 (%s), 回退 suvtk 包装器", e)
                if not produced:
                    sqn_cmd = (
                        f"suvtk table2asn "
                        f"-i {reoriented_fna} "
                        f"-o {sv_dir}/submission "
                        f"-s {src_file} "
                        f"-f {updated_tbl} "
                        f"-t {sbt_file} "
                        f"-c {comments_file}"
                    )
                    if run(sqn_cmd, log, "suvtk table2asn"):
                        produced = (sv_dir / "submission.sqn").is_file() and \
                                   (sv_dir / "submission.sqn").stat().st_size > 0
                if produced:
                    log.info("  ✓ submission.sqn 已生成 (%d KB)",
                             (sv_dir / "submission.sqn").stat().st_size // 1024)
                else:
                    log.warning("  submission.sqn 未生成, 请检查")
                    submission_failures.append("submission.sqn 未生成 (table2asn 两级均未产出)")
        else:
            missing = []
            if not reoriented_fna or not reoriented_fna.exists(): missing.append("reoriented.fna")
            if not updated_tbl or not os.path.exists(updated_tbl): missing.append("featuretable.tbl")
            log.warning("  缺少输入文件, 跳过 table2asn: %s", ", ".join(missing))
            submission_failures.append("缺输入文件: " + ", ".join(missing))

    # ── 统一生成 topology.tsv (circular/linear, 供 sequin .fsa 头部与报告复用) ──
    topo_out = sub_out / "topology.tsv"
    topo_script = SCRIPT_DIR / "viral_topology.py"
    if tax_tsv and fna and topo_script.exists() and not topo_out.exists():
        run(f"python {topo_script} --taxonomy {tax_tsv} --fasta {fna} -o {topo_out}",
            log, "viral_topology")

    if args.mode in ('sequin', 'both'):
        log.info("─" * 40 + "\n  Sequin 模式: 生成 .fsa + .tbl + .cmt (Cenote-Taker3 风格)\n" + "─" * 40)

        sequin_script = SCRIPT_DIR / "sequin_builder.py"

        if tax_tsv and miuvig_tax and fna and updated_tbl and miuvig_feat:
            cmd = (
                f"python {sequin_script} "
                f"--taxonomy {tax_tsv} "
                f"--miuvig-tax {miuvig_tax} "
                f"--features {updated_tbl} "
                f"--proteins {updated_faa} "
                f"--fasta {fna} "
                f"--miuvig-feat {miuvig_feat} "
                f"--run-title {args.run_title} "
            )
            if topo_out.exists():
                cmd += f" --topology {topo_out} "
            cmd += f"-o {sub_out}/"
            run(cmd, log, "sequin_builder")

    # ── 统一元数据 CSV 导出 (SeqSender 兼容) ──
    if args.export_metadata and tax_tsv:
        log.info("─" * 40 + "\n  导出统一元数据 CSV (SeqSender 兼容)\n" + "─" * 40)

        pipe_params = load_pipeline_params(work_dir)

        # 根据管道类型自动决定 source_uvig 和 enrichment
        if pipeline_type == 'discovery':
            params = pipe_params.copy()
        else:
            # analysis pipeline: 自己测序, 可能做了病毒富集
            params = pipe_params.copy()

        # sub_out 已在上方定义

        meta_script = SCRIPT_DIR / "unified_metadata.py"
        cmd = (
            f"python {meta_script} "
            f"--taxonomy {tax_tsv} "
            f"--run-title {args.run_title} "
            f"--assembler \"{params['assembler']}\" "
            f"--sequencer \"{params['sequencer']}\" "
            f"--enrichment \"{params['enrichment']}\" "
            f"--source-uvig \"{params['source_uvig']}\" "
            f"--metagenome-source \"{params['metagenome_source']}\" "
            f"--pipeline \"{params['pipeline']}\" "
            f"-o {sub_out}/"
        )

        if core14.exists():
            cmd += f" --reference-metadata {core14}"
        if args.config:
            cmd += f" --config {args.config}"
        if args.interactive:
            cmd += " --interactive"

        run(cmd, log, "unified_metadata")

        # ── 病毒拓扑判断 (复用上方统一生成的 topology.tsv) ──

        # ── 交互式 HTML 报告 (argv 直连, 不经 shell; 存在的产物全部透传) ──
        csv_path = sub_out / "unified_metadata.csv"
        miuvig_path = sub_out / "miuvig.tsv"
        asm_path = sub_out / "assembly.tsv"
        report_script = SCRIPT_DIR / "report_html.py"
        if csv_path.exists():
            argv = [sys.executable, str(report_script), "--csv", str(csv_path),
                    "--run-title", args.run_title, "-o", str(report_path)]
            suvtk_tax = _first_existing(miuvig_tax, tax_tsv)
            suvtk_feat = _first_existing(
                (feat_dir / "featuretable.tbl") if feat_dir else None, miuvig_feat)
            hypo_tbl = _first_existing(
                (hypo_dir / "featuretable_updated.tbl") if hypo_dir else None)
            for flag, p in (("--miuvig", _first_existing(miuvig_path)),
                            ("--assembly", _first_existing(asm_path)),
                            ("--topology", _first_existing(topo_out)),
                            ("--suvtk-tax", suvtk_tax),
                            ("--suvtk-feat", suvtk_feat),
                            ("--hypo", hypo_tbl),
                            ("--gb-dir", sub_out / "suvtk_submission" / "virus-annotations")):
                if p is None:
                    continue
                if flag == "--gb-dir":
                    if os.path.isdir(p):
                        argv += [flag, str(p)]
                elif Path(p).is_file():
                    argv += [flag, str(p)]
            run_argv(argv, log, "report_html")

    # ═══ 最终报告 ═══
    log.info("=" * 60)
    log.info("完成!")
    log.info("  分析目录: %s", work_dir)
    log.info("  提交输出: %s", sub_out)
    if report_path.exists():
        log.info("  📊 交互报告:  %s", report_path)
    if tax_dir:
        log.info("  [1] taxonomy:    %s", tax_dir)
    if feat_dir:
        log.info("  [2] features:    %s", feat_dir)
    if hypo_dir:
        log.info("  [2.5] hypo:      %s", hypo_dir)
    if args.mode in ('suvtk', 'both'):
        log.info("  [suvtk]  report:  %s/report/", work_dir)
    if args.mode in ('sequin', 'both'):
        gb_dir = work_dir / "genbank_submission" / args.run_title
        log.info("  [sequin] output:  %s/", gb_dir)
    log.info("=" * 60)

    # 聚合退出码 (2026-10-15): 登记的未完成项非空 → 非零退出。
    # 旧行为 gate 失败 / .sqn 未生成只打日志, main 正常返回 0 —— 批处理把
    # 失败提交当成功 (与 ⑤ 的 0/1/2 语义、③ 的 failed_stages 口径对齐)。
    if submission_failures:
        log.error("提交流程存在未完成项 (退出码 1): %s", "; ".join(submission_failures))
        sys.exit(1)


if __name__ == '__main__':
    main()
