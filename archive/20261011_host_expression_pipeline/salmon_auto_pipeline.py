#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Salmon RNA-seq one-click quantification workflow — 宿主转录组定量引擎 (⑧ 宿主表达模块)

来源与改动 (provenance)
------------------------------------------------------------------
本脚本自 goji-virome 项目的同名脚本原样抽取 (枸杞病毒宿主表达分析,
项目管理目录 D:\\桌面\\差异表达分析\\salmon_auto_pipeline.py), 逻辑与
原版逐行一致, 仅两处适配:

  1. 后缀识别增加 FASTA 扩展名 (.fa/.fna/.fasta, 含 .gz) —— MMPV ② 管线
     clean-data 的输出为 FASTA (seqkit fq2fa, 质量值已剥除), 原版只认
     FASTQ 后缀会识别不到样本; salmon 不使用质量值, 定量不受影响。
  2. 常量 READ_SUFFIXES 更名 READ_SUFFIXES (语义已不限于 FASTQ)。
  3. 新增 salmon 定量后 QC 闸门 (quant_qc.tsv / qc_failed_samples.tsv;
     默认 min_library_size=1e6, 未通过者从主矩阵剔除, 另存 *_matrix_all_samples.tsv)。
     动因: 枸杞队列 19 个近空文库样本 (最小 1,065 reads) 曾混入下游,
     污染共表达网络与分区结论。

功能：
1. 扫描输入目录中的 reads 文件 (FASTQ / FASTA, 均支持 .gz)
2. 自动识别单端（SE）/双端（PE）数据
3. 自动检测已有 Salmon index；若无索引，则根据参考 transcript/CDS FASTA 自动建立
4. 自动并行执行 salmon quant
5. 自动整理 quant.sf
6. 自动汇总 TPM / NumReads(count) / FPKM 表达矩阵
7. 输出样本信息、命令和日志

支持常见文件名：
    sample_R1.fastq.gz / sample_R2.fastq.gz
    sample_R1_001.fastq.gz / sample_R2_001.fastq.gz
    sample_L001_R1_001.fastq.gz / sample_L001_R2_001.fastq.gz
    sample_1.fastq.gz / sample_2.fastq.gz          (② clean data: sample_1.fa.gz / sample_2.fa.gz)
    sample.R1.fq.gz / sample.R2.fq.gz
    单端文件：sample.fastq.gz

示例：
    python salmon_auto_pipeline.py \
        -i /data/clean_reads \
        -r /data/host_cds.fa \
        -o /data/host_expression/01_Quant

若已有索引：
    python salmon_auto_pipeline.py \
        -i /data/clean_reads \
        --index /data/salmon_index \
        -o /data/host_expression/01_Quant

依赖：
    Salmon 必须已安装并能在命令行直接运行：
        salmon --version

本脚本仅使用 Python 标准库，不依赖 pandas。
"""

from __future__ import annotations

import argparse
import csv
import json
import logging
import os
import re
import shutil
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple

# 可选依赖：tqdm 进度条（未安装时自动降级，不影响功能）
try:
    from tqdm import tqdm as _tqdm

    _HAS_TQDM = True
except ImportError:
    _HAS_TQDM = False


def make_pbar(total: int, desc: str):
    """返回 tqdm 进度条；无 tqdm 或非交互终端（如 nohup 重定向）时返回 None。"""
    if not _HAS_TQDM:
        return None
    return _tqdm(
        total=total,
        desc=desc,
        unit="sample",
        disable=not sys.stderr.isatty(),
    )


READ_SUFFIXES = (
    ".fastq.gz",
    ".fq.gz",
    ".fastq",
    ".fq",
    ".fasta.gz",
    ".fa.gz",
    ".fna.gz",
    ".fasta",
    ".fa",
    ".fna",
)


# ---------------------------------------------------------------------------
# Salmon 定量后 QC 默认门槛
#
# 依据（宁夏枸杞 RNA-Lycium_barbarum，445 样本实测分布）：
#   文库量中位 ~10.6M；detected_genes 中位 ~26,000；percent_mapped 中位 ~61%
#   关键事实：lib < 1e5 的 54 个样本与 lib >= 2.4e6 的正常样本之间存在
#   **约 100 倍的空档**（最差正常样本 2,438,517，最好问题样本 15,860）。
#   因此门槛取 1e6 落在空档中央，在 [1e5, 2.4e6] 区间内取任何值结果完全一致。
#
# 历史教训：原先只用「lib == 0」做过滤，导致 19 个近空文库样本
# （最小 1,065 reads，全部为单端、全部病毒阴性）进入下游分析，
# 污染了表达矩阵、共表达网络、模块与分区结论。
# ---------------------------------------------------------------------------
DEFAULT_MIN_LIBRARY_SIZE = 1_000_000     # 硬门槛：ΣNumReads
DEFAULT_WARN_LIBRARY_SIZE = 5_000_000    # 软门槛（仅告警，不剔除）
DEFAULT_MIN_DETECTED_GENES = 5_000       # 硬门槛：NumReads >= 1 的基因数
DEFAULT_WARN_DETECTED_GENES = 20_000
DEFAULT_MIN_MAPPING_RATE = 10.0          # 硬门槛：percent_mapped
DEFAULT_WARN_MAPPING_RATE = 30.0


@dataclass
class Sample:
    name: str
    layout: str  # PE or SE
    read1: List[Path] = field(default_factory=list)
    read2: List[Path] = field(default_factory=list)
    single: List[Path] = field(default_factory=list)
    missing_read: Optional[str] = None  # 缺端的 read 标记：'R1' 或 'R2'；PE/正常 SE 为 None


def setup_logging(log_file: Path, verbose: bool = False) -> None:
    log_file.parent.mkdir(parents=True, exist_ok=True)
    level = logging.DEBUG if verbose else logging.INFO

    formatter = logging.Formatter(
        "%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    root = logging.getLogger()
    root.setLevel(level)
    root.handlers.clear()

    console = logging.StreamHandler(sys.stdout)
    console.setFormatter(formatter)
    root.addHandler(console)

    file_handler = logging.FileHandler(log_file, encoding="utf-8")
    file_handler.setFormatter(formatter)
    root.addHandler(file_handler)


def strip_fastq_suffix(filename: str) -> Optional[str]:
    lower = filename.lower()
    for suffix in READ_SUFFIXES:
        if lower.endswith(suffix):
            return filename[: -len(suffix)]
    return None


def natural_key(text: str):
    return [
        int(x) if x.isdigit() else x.lower()
        for x in re.split(r"(\d+)", text)
    ]


def fastq_files(input_dir: Path, recursive: bool) -> List[Path]:
    iterator = input_dir.rglob("*") if recursive else input_dir.glob("*")
    files = []

    for path in iterator:
        if not path.is_file():
            continue

        low = path.name.lower()
        if any(low.endswith(s) for s in READ_SUFFIXES):
            # 注意：不要 path.resolve()。
            # 输入目录中的 FASTQ 可能是符号链接（如 CRR1126132_R1.fq.gz -> x.fa.gz），
            # resolve 会换成真实文件名，导致后缀匹配失败、样本无法识别。
            # input_dir 在 main() 中已 resolve 为绝对路径，此处 path 即为绝对路径。
            files.append(path)

    return sorted(files, key=lambda x: natural_key(str(x)))


def parse_read_name(stem: str) -> Optional[Tuple[str, str]]:
    """
    返回 (sample_name, read_number) 或 None。

    优先匹配明确的 R1/R2：
        sample_R1
        sample_R1_001
        sample_L001_R1_001
        sample.R1
        sample-R1

    再匹配 _1/_2、.1/.2、-1/-2。
    """

    patterns = [
        # sample_L001_R1_001
        re.compile(
            r"^(?P<sample>.+?)(?:[_\.-]L\d{3})?[_\.-]R(?P<read>[12])"
            r"(?:[_\.-]\d{3})?$",
            re.IGNORECASE,
        ),
        # sample_R1 / sample.R1 / sample-R1
        re.compile(
            r"^(?P<sample>.+?)[_\.-]R(?P<read>[12])$",
            re.IGNORECASE,
        ),
        # sample_1_001 / sample_1
        re.compile(
            r"^(?P<sample>.+?)(?:[_\.-]L\d{3})?[_\.-](?P<read>[12])"
            r"(?:[_\.-]\d{3})?$",
            re.IGNORECASE,
        ),
    ]

    for pattern in patterns:
        m = pattern.match(stem)
        if m:
            sample = m.group("sample").rstrip("._-")
            read = m.group("read")
            if sample:
                return sample, read

    return None


def detect_samples(input_dir: Path, recursive: bool = False) -> List[Sample]:
    """
    自动检测样本。

    规则：
    - 能找到匹配 R1/R2 的样本 -> PE
    - 只有一侧 read 的文件 -> 当作 SE，并发出警告
    - 文件名中完全没有 R1/R2/1/2 标记 -> SE
    - 同一样本多个 lane 会自动合并为一个 Salmon 样本
    """

    files = fastq_files(input_dir, recursive)

    if not files:
        raise RuntimeError(
            f"输入目录中没有发现 FASTQ 文件：{input_dir}"
        )

    paired_groups: Dict[str, Dict[str, List[Path]]] = {}
    unclassified: List[Path] = []

    for path in files:
        stem = strip_fastq_suffix(path.name)
        if stem is None:
            continue

        parsed = parse_read_name(stem)

        if parsed is None:
            unclassified.append(path)
            continue

        sample_name, read_num = parsed

        if sample_name not in paired_groups:
            paired_groups[sample_name] = {"1": [], "2": []}

        paired_groups[sample_name][read_num].append(path)

    samples: List[Sample] = []

    # 处理有 R1/R2 标记的数据
    for sample_name in sorted(paired_groups, key=natural_key):
        group = paired_groups[sample_name]

        r1 = sorted(group["1"], key=lambda x: natural_key(str(x)))
        r2 = sorted(group["2"], key=lambda x: natural_key(str(x)))

        if r1 and r2:
            if len(r1) != len(r2):
                raise RuntimeError(
                    f"样本 {sample_name} 的 R1/R2 文件数量不一致："
                    f"R1={len(r1)}, R2={len(r2)}"
                )

            samples.append(
                Sample(
                    name=sample_name,
                    layout="PE",
                    read1=r1,
                    read2=r2,
                )
            )

        elif r1:
            logging.warning(
                "样本 %s 仅检测到 R1，没有检测到 R2；"
                "将按单端 SE 数据处理（缺失端已标记）。",
                sample_name,
            )
            samples.append(
                Sample(
                    name=sample_name,
                    layout="SE",
                    single=r1,
                    missing_read="R2",
                )
            )

        elif r2:
            logging.warning(
                "样本 %s 仅检测到 R2，没有检测到 R1；"
                "将按单端 SE 数据处理（缺失端已标记）。",
                sample_name,
            )
            samples.append(
                Sample(
                    name=sample_name,
                    layout="SE",
                    single=r2,
                    missing_read="R1",
                )
            )

    # 无 R1/R2 标记的文件，每个文件作为一个单端样本
    for path in unclassified:
        stem = strip_fastq_suffix(path.name)
        sample_name = stem if stem else path.stem

        # 防止与已检测样本重名
        existing = {s.name for s in samples}
        original_name = sample_name
        n = 2

        while sample_name in existing:
            sample_name = f"{original_name}_{n}"
            n += 1

        samples.append(
            Sample(
                name=sample_name,
                layout="SE",
                single=[path],
            )
        )

    samples.sort(key=lambda s: natural_key(s.name))

    if not samples:
        raise RuntimeError("未能从输入目录识别任何样本。")

    return samples


def samples_from_quant_dir(quant_root: Path) -> List[Sample]:
    """
    --qc-only 模式：直接从已有 quant 目录重建样本清单。

    用于给历史结果补做质控（原始 FASTQ 可能已删除/迁移）。
    Layout 由各样本的 cmd_info.json 判断：
      含 "mates"       -> PE
      含 "unmatedReads" -> SE
    """

    if not quant_root.is_dir():
        raise RuntimeError(f"quant 目录不存在：{quant_root}")

    samples: List[Sample] = []

    for child in sorted(
        (p for p in quant_root.iterdir() if p.is_dir()),
        key=lambda p: natural_key(p.name),
    ):
        if not (child / "quant.sf").exists():
            continue

        layout = "SE"
        cmd_info = child / "cmd_info.json"

        if cmd_info.exists():
            try:
                info = json.loads(cmd_info.read_text(encoding="utf-8"))
                # salmon 的 cmd_info.json：PE 记 mates1/mates2，SE 记 unmatedReads
                if info.get("mates1") or info.get("mates2"):
                    layout = "PE"
            except (ValueError, OSError):
                pass

        samples.append(Sample(name=child.name, layout=layout))

    if not samples:
        raise RuntimeError(f"quant 目录中没有可用样本：{quant_root}")

    return samples


def salmon_exists(salmon_bin: str) -> str:
    path = shutil.which(salmon_bin)

    if path is None:
        raise RuntimeError(
            f"找不到 Salmon 可执行程序：{salmon_bin}\n"
            "请先安装 Salmon，并确认 `salmon --version` 可以正常运行。"
        )

    return path


def salmon_version(salmon_bin: str) -> str:
    result = subprocess.run(
        [salmon_bin, "--version"],
        capture_output=True,
        text=True,
        check=False,
    )

    text = (result.stdout or result.stderr or "").strip()
    return text if text else "unknown"


def valid_salmon_index(index_dir: Path) -> bool:
    if not index_dir.is_dir():
        return False

    # Salmon/pufferfish index 通常包含 versionInfo.json。
    # 同时要求目录不为空，避免误把普通目录当作 index。
    return (
        (index_dir / "versionInfo.json").exists()
        and any(index_dir.iterdir())
    )


def find_existing_index(
    explicit_index: Optional[Path],
    output_dir: Path,
    input_dir: Path,
    reference: Optional[Path],
) -> Optional[Path]:

    candidates: List[Path] = []

    if explicit_index is not None:
        candidates.append(explicit_index)

    candidates.append(output_dir / "salmon_index")
    candidates.append(input_dir / "salmon_index")

    if reference is not None:
        candidates.append(reference.parent / "salmon_index")

    # 去重
    seen = set()
    unique_candidates = []

    for c in candidates:
        c = c.expanduser().resolve()
        if str(c) not in seen:
            unique_candidates.append(c)
            seen.add(str(c))

    for candidate in unique_candidates:
        if valid_salmon_index(candidate):
            return candidate

    # 再扫描上述几个位置的一级子目录
    search_parents = [output_dir, input_dir]

    if reference is not None:
        search_parents.append(reference.parent)

    for parent in search_parents:
        if not parent.exists() or not parent.is_dir():
            continue

        try:
            for child in parent.iterdir():
                if child.is_dir() and valid_salmon_index(child):
                    logging.info("自动发现 Salmon index：%s", child)
                    return child.resolve()
        except PermissionError:
            continue

    return None


def build_index(
    salmon_bin: str,
    reference: Path,
    index_dir: Path,
    threads: int,
    force: bool = False,
) -> Path:

    if valid_salmon_index(index_dir) and not force:
        logging.info("Salmon index 已存在：%s", index_dir)
        return index_dir

    if not reference.exists():
        raise RuntimeError(f"参考 FASTA 不存在：{reference}")

    if force and index_dir.exists():
        logging.warning("强制重建 index，删除旧目录：%s", index_dir)
        shutil.rmtree(index_dir)

    index_dir.parent.mkdir(parents=True, exist_ok=True)

    cmd = [
        salmon_bin,
        "index",
        "-t",
        str(reference),
        "-i",
        str(index_dir),
        "-p",
        str(threads),
        # 默认保留重复转录本（salmon index 默认会移除序列重复的转录本，
        # 对含多拷贝基因的转录组会低估；与 quant 阶段的 --keepDuplicates 配套）。
        "--keepDuplicates",
    ]

    logging.info("开始建立 Salmon index。")
    logging.info("命令：%s", " ".join(cmd))

    subprocess.run(cmd, check=True)

    if not valid_salmon_index(index_dir):
        raise RuntimeError(
            "Salmon index 命令执行结束，但没有检测到有效 index。"
        )

    return index_dir.resolve()


def sanitize_sample_name(name: str) -> str:
    safe = re.sub(r"[^A-Za-z0-9._-]+", "_", name)
    safe = safe.strip("._-")
    return safe or "sample"


def run_quant(
    salmon_bin: str,
    index_dir: Path,
    sample: Sample,
    quant_root: Path,
    log_root: Path,
    threads: int,
    fld_mean: float,
    fld_sd: float,
    force: bool,
    validate_mappings: bool,
) -> Tuple[str, str]:

    sample_name = sanitize_sample_name(sample.name)
    sample_out = quant_root / sample_name
    quant_sf = sample_out / "quant.sf"
    log_file = log_root / f"{sample_name}.salmon.log"

    if quant_sf.exists() and quant_sf.stat().st_size > 0 and not force:
        logging.info("[%s] quant.sf 已存在，跳过。", sample_name)
        return sample_name, "SKIPPED"

    if force and sample_out.exists():
        shutil.rmtree(sample_out)

    sample_out.mkdir(parents=True, exist_ok=True)
    log_root.mkdir(parents=True, exist_ok=True)

    cmd = [
        salmon_bin,
        "quant",
        "-i",
        str(index_dir),
        "-l",
        "A",
        "-p",
        str(threads),
        # 注意：--keepDuplicates 仅存在于 salmon index 阶段（保留重复转录本）；
        # salmon quant 没有该选项，加在 quant 命令会导致 unrecognised option 直接退出。
    ]

    if validate_mappings:
        cmd.append("--validateMappings")

    if sample.layout == "PE":
        left = ",".join(str(x) for x in sample.read1)
        right = ",".join(str(x) for x in sample.read2)

        cmd.extend(
            [
                "-1",
                left,
                "-2",
                right,
            ]
        )

    else:
        reads = ",".join(str(x) for x in sample.single)

        cmd.extend(
            [
                "-r",
                reads,
                "--fldMean",
                str(fld_mean),
                "--fldSD",
                str(fld_sd),
            ]
        )

    cmd.extend(
        [
            "-o",
            str(sample_out),
        ]
    )

    logging.info(
        "[%s] 开始 Salmon quant（%s）。",
        sample_name,
        sample.layout,
    )

    with log_file.open("w", encoding="utf-8") as log:
        log.write("COMMAND:\n")
        log.write(" ".join(cmd) + "\n\n")
        log.flush()

        process = subprocess.run(
            cmd,
            stdout=log,
            stderr=subprocess.STDOUT,
            text=True,
            check=False,
        )

    if process.returncode != 0:
        raise RuntimeError(
            f"样本 {sample_name} Salmon quant 失败，"
            f"退出码={process.returncode}。日志：{log_file}"
        )

    if not quant_sf.exists() or quant_sf.stat().st_size == 0:
        raise RuntimeError(
            f"样本 {sample_name} Salmon quant 未产生有效 quant.sf。"
        )

    logging.info("[%s] 完成。", sample_name)
    return sample_name, "DONE"


def sample_status(sample: Sample) -> str:
    """样本识别状态：PE / SE / SE_missing_R1 / SE_missing_R2"""
    if sample.layout == "PE":
        return "PE"
    if sample.missing_read is not None:
        return f"SE_missing_{sample.missing_read}"
    return "SE"


def write_sample_info(samples: List[Sample], output_file: Path) -> None:
    output_file.parent.mkdir(parents=True, exist_ok=True)

    with output_file.open("w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f, delimiter="\t")

        writer.writerow(
            [
                "Sample",
                "Layout",
                "Status",
                "Read1_or_SE",
                "Read2",
                "NumFiles",
            ]
        )

        for sample in samples:
            if sample.layout == "PE":
                left = ",".join(str(p) for p in sample.read1)
                right = ",".join(str(p) for p in sample.read2)
                n = len(sample.read1) + len(sample.read2)
            else:
                left = ",".join(str(p) for p in sample.single)
                right = ""
                n = len(sample.single)

            writer.writerow(
                [
                    sample.name,
                    sample.layout,
                    sample_status(sample),
                    left,
                    right,
                    n,
                ]
            )


def read_quant_sf(path: Path):
    records = {}

    with path.open("r", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")

        required = {
            "Name",
            "Length",
            "EffectiveLength",
            "TPM",
            "NumReads",
        }

        missing = required - set(reader.fieldnames or [])

        if missing:
            raise RuntimeError(
                f"{path} 缺少 quant.sf 必需字段："
                + ", ".join(sorted(missing))
            )

        for row in reader:
            name = row["Name"]

            records[name] = {
                "Length": float(row["Length"]),
                "EffectiveLength": float(row["EffectiveLength"]),
                "TPM": float(row["TPM"]),
                "NumReads": float(row["NumReads"]),
            }

    return records


def write_matrix(
    output_file: Path,
    feature_order: List[str],
    sample_names: List[str],
    values: Dict[str, Dict[str, float]],
) -> None:

    output_file.parent.mkdir(parents=True, exist_ok=True)

    with output_file.open("w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["Name"] + sample_names)

        for feature in feature_order:
            row = [feature]

            for sample in sample_names:
                value = values.get(feature, {}).get(sample, 0.0)

                # 保留足够精度，避免表达量被过度四舍五入
                row.append(f"{value:.10g}")

            writer.writerow(row)


def collect_quant_metrics(
    quant_root: Path,
    samples: List[Sample],
) -> List[Dict[str, object]]:
    """
    读取每个样本的 salmon quant 产出，汇总定量 QC 指标。

    指标来源：
      - library_size / detected_genes：quant.sf 的 NumReads 列
      - percent_mapped：aux_info/meta_info.json（salmon 自带比对率）
    """

    metrics: List[Dict[str, object]] = []

    for sample in samples:
        sample_name = sanitize_sample_name(sample.name)
        sample_out = quant_root / sample_name
        quant_sf = sample_out / "quant.sf"

        if not quant_sf.exists() or quant_sf.stat().st_size == 0:
            logging.warning(
                "[%s] 缺少有效 quant.sf，QC 时跳过。",
                sample_name,
            )
            continue

        library_size = 0
        detected_genes = 0

        with quant_sf.open("r", encoding="utf-8", newline="") as f:
            reader = csv.DictReader(f, delimiter="\t")

            for row in reader:
                num_reads = float(row.get("NumReads", 0.0) or 0.0)
                library_size += num_reads

                if num_reads >= 1.0:
                    detected_genes += 1

        percent_mapped = float("nan")
        meta_info = sample_out / "aux_info" / "meta_info.json"

        if meta_info.exists():
            try:
                meta = json.loads(meta_info.read_text(encoding="utf-8"))
                percent_mapped = float(meta.get("percent_mapped", "nan"))
            except (ValueError, OSError, TypeError):
                percent_mapped = float("nan")

        metrics.append(
            {
                "Sample": sample_name,
                "Layout": sample.layout,
                "Status": sample_status(sample),
                "library_size": int(library_size),
                "detected_genes": detected_genes,
                "percent_mapped": percent_mapped,
            }
        )

    return metrics


def evaluate_quant_qc(
    metrics: List[Dict[str, object]],
    min_library_size: int = DEFAULT_MIN_LIBRARY_SIZE,
    warn_library_size: int = DEFAULT_WARN_LIBRARY_SIZE,
    min_detected_genes: int = DEFAULT_MIN_DETECTED_GENES,
    warn_detected_genes: int = DEFAULT_WARN_DETECTED_GENES,
    min_mapping_rate: float = DEFAULT_MIN_MAPPING_RATE,
    warn_mapping_rate: float = DEFAULT_WARN_MAPPING_RATE,
) -> List[Dict[str, object]]:
    """
    按门槛给每个样本打 PASS / WARN / FAIL 并记录原因。

    FAIL 的样本会从主表达矩阵中剔除；WARN 仅记录，不剔除。
    """

    evaluated: List[Dict[str, object]] = []

    for m in metrics:
        reasons: List[str] = []
        hard_fail = False

        lib = int(m["library_size"])

        if lib < min_library_size:
            reasons.append(f"library_size={lib}<{min_library_size}")
            hard_fail = True
        elif lib < warn_library_size:
            reasons.append(f"library_size={lib}<{warn_library_size}")

        det = int(m["detected_genes"])

        if det < min_detected_genes:
            reasons.append(f"detected_genes={det}<{min_detected_genes}")
            hard_fail = True
        elif det < warn_detected_genes:
            reasons.append(f"detected_genes={det}<{warn_detected_genes}")

        pct = m["percent_mapped"]

        if isinstance(pct, float) and pct == pct:  # 非 NaN
            if pct < min_mapping_rate:
                reasons.append(f"percent_mapped={pct:.1f}<{min_mapping_rate}")
                hard_fail = True
            elif pct < warn_mapping_rate:
                reasons.append(f"percent_mapped={pct:.1f}<{warn_mapping_rate}")

        # 缺端样本（R2/R1 缺失、被迫按单端跑）本身是数据缺陷，记入原因但不单独判死
        if str(m.get("Status", "")).startswith("SE_missing_"):
            reasons.append(str(m["Status"]))

        if hard_fail:
            status = "FAIL"
        elif reasons:
            status = "WARN"
        else:
            status = "PASS"

        evaluated.append(
            {**m, "QC_status": status, "QC_reason": ";".join(reasons)}
        )

    return evaluated


def write_quant_qc_report(
    evaluated: List[Dict[str, object]],
    metadata_root: Path,
) -> Tuple[Path, List[str], List[str]]:
    """
    写出定量 QC 报告，返回 (报告路径, 通过样本名单, 未通过样本名单)。

    产物：
      - metadata/quant_qc.tsv           全部样本的指标 + PASS/WARN/FAIL + 原因
      - metadata/qc_pass_samples.txt    可进入下游的样本（PASS + WARN）
      - metadata/qc_failed_samples.tsv  被剔除的样本及原因
    """

    metadata_root.mkdir(parents=True, exist_ok=True)

    qc_path = metadata_root / "quant_qc.tsv"

    columns = [
        "Sample",
        "Layout",
        "Status",
        "library_size",
        "detected_genes",
        "percent_mapped",
        "QC_status",
        "QC_reason",
    ]

    with qc_path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(columns)

        for m in evaluated:
            pct = m["percent_mapped"]
            pct_text = f"{pct:.2f}" if isinstance(pct, float) and pct == pct else "NA"

            writer.writerow(
                [
                    m["Sample"],
                    m["Layout"],
                    m["Status"],
                    m["library_size"],
                    m["detected_genes"],
                    pct_text,
                    m["QC_status"],
                    m["QC_reason"],
                ]
            )

    passed = [str(m["Sample"]) for m in evaluated if m["QC_status"] != "FAIL"]
    failed = [m for m in evaluated if m["QC_status"] == "FAIL"]

    with (metadata_root / "qc_pass_samples.txt").open("w", encoding="utf-8") as f:
        for name in passed:
            f.write(f"{name}\n")

    with (metadata_root / "qc_failed_samples.tsv").open(
        "w", encoding="utf-8", newline=""
    ) as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(
            [
                "Sample",
                "library_size",
                "detected_genes",
                "percent_mapped",
                "QC_reason",
            ]
        )

        for m in failed:
            pct = m["percent_mapped"]
            pct_text = f"{pct:.2f}" if isinstance(pct, float) and pct == pct else "NA"

            writer.writerow(
                [
                    m["Sample"],
                    m["library_size"],
                    m["detected_genes"],
                    pct_text,
                    m["QC_reason"],
                ]
            )

    n_pass = sum(1 for m in evaluated if m["QC_status"] == "PASS")
    n_warn = sum(1 for m in evaluated if m["QC_status"] == "WARN")
    n_fail = len(failed)

    logging.info(
        "定量 QC 完成：PASS=%d，WARN=%d，FAIL=%d（共 %d 个样本）。",
        n_pass,
        n_warn,
        n_fail,
        len(evaluated),
    )

    if n_fail:
        logging.warning(
            "%d 个样本未通过定量 QC，已从主表达矩阵剔除；"
            "明细见：%s",
            n_fail,
            metadata_root / "qc_failed_samples.tsv",
        )

    return qc_path, passed, [str(m["Sample"]) for m in failed]


def merge_quant_results(
    quant_root: Path,
    matrix_root: Path,
    expected_samples: List[str],
    sample_filter: Optional[set] = None,
    write_all_samples: bool = True,
) -> List[str]:
    """
    汇总 quant.sf 为表达矩阵。

    sample_filter 为 None 时，所有样本写入主矩阵；
    否则主矩阵只含 sample_filter 内的样本，并另写
    ``*_matrix_all_samples.tsv`` 保存全部样本（便于回溯与对比）。

    返回实际写入主矩阵的样本名单。
    """

    logging.info("开始汇总 quant.sf 为表达矩阵。")

    sample_records = {}
    sample_names = []
    feature_order = []
    feature_seen = set()

    for raw_name in expected_samples:
        sample_name = sanitize_sample_name(raw_name)
        quant_sf = quant_root / sample_name / "quant.sf"

        if not quant_sf.exists():
            logging.warning(
                "样本 %s 缺少 quant.sf，矩阵汇总时跳过。",
                sample_name,
            )
            continue

        records = read_quant_sf(quant_sf)
        sample_records[sample_name] = records
        sample_names.append(sample_name)

        for feature in records:
            if feature not in feature_seen:
                feature_seen.add(feature)
                feature_order.append(feature)

    if not sample_records:
        raise RuntimeError("没有可用于汇总的 quant.sf。")

    tpm_values: Dict[str, Dict[str, float]] = {}
    count_values: Dict[str, Dict[str, float]] = {}
    fpkm_values: Dict[str, Dict[str, float]] = {}
    effective_length_values: Dict[str, Dict[str, float]] = {}

    for sample_name, records in sample_records.items():
        total_numreads = sum(
            max(0.0, x["NumReads"])
            for x in records.values()
        )

        for feature, rec in records.items():
            length = rec["Length"]
            tpm = rec["TPM"]
            count = rec["NumReads"]
            effective_length = rec["EffectiveLength"]

            tpm_values.setdefault(feature, {})[sample_name] = tpm
            count_values.setdefault(feature, {})[sample_name] = count
            effective_length_values.setdefault(feature, {})[
                sample_name
            ] = effective_length

            # transcript-level FPKM：
            # FPKM = NumReads * 1e9 / (transcript_length_bp * total_NumReads)
            if total_numreads > 0 and length > 0:
                fpkm = count * 1e9 / (length * total_numreads)
            else:
                fpkm = 0.0

            fpkm_values.setdefault(feature, {})[sample_name] = fpkm

    if sample_filter is None:
        kept_samples = sample_names
    else:
        kept_samples = [s for s in sample_names if s in sample_filter]

        dropped = len(sample_names) - len(kept_samples)

        if dropped:
            logging.warning(
                "QC 过滤：%d 个样本未通过门槛，已从主矩阵剔除"
                "（仍保留在 *_matrix_all_samples.tsv）。",
                dropped,
            )

    if not kept_samples:
        raise RuntimeError(
            "QC 过滤后没有剩余样本可用于汇总表达矩阵。"
        )

    matrix_root.mkdir(parents=True, exist_ok=True)

    def write_pair(stem: str, values: Dict[str, Dict[str, float]]) -> None:
        write_matrix(
            matrix_root / f"{stem}_matrix.tsv",
            feature_order,
            kept_samples,
            values,
        )

        if write_all_samples and sample_filter is not None:
            write_matrix(
                matrix_root / f"{stem}_matrix_all_samples.tsv",
                feature_order,
                sample_names,
                values,
            )

    write_pair("TPM", tpm_values)
    write_pair("NumReads", count_values)
    write_pair("FPKM", fpkm_values)
    write_pair("EffectiveLength", effective_length_values)

    # 同时输出 sample.list，兼容后续 kallisto2matrix 等工具
    sample_list = matrix_root / "sample.list"

    with sample_list.open("w", encoding="utf-8") as f:
        for sample_name in kept_samples:
            quant_sf = quant_root / sample_name / "quant.sf"
            f.write(f"{quant_sf},{sample_name}\n")

    logging.info(
        "表达矩阵汇总完成：%s（%d 个样本）",
        matrix_root,
        len(kept_samples),
    )

    return kept_samples


def write_run_summary(
    summary_file: Path,
    salmon_ver: str,
    input_dir: Path,
    output_dir: Path,
    index_dir: Path,
    reference: Optional[Path],
    samples: List[Sample],
    jobs: int,
    threads_per_sample: int,
    extra: Optional[Dict[str, object]] = None,
) -> None:

    pe = sum(1 for s in samples if s.layout == "PE")
    se = sum(1 for s in samples if s.layout == "SE")

    data = {
        "salmon_version": salmon_ver,
        "input_dir": str(input_dir) if input_dir else None,
        "output_dir": str(output_dir),
        "reference": str(reference) if reference else None,
        "index_dir": str(index_dir),
        "samples_total": len(samples),
        "paired_end_samples": pe,
        "single_end_samples": se,
        "parallel_jobs": jobs,
        "threads_per_sample": threads_per_sample,
        "estimated_max_threads": jobs * threads_per_sample,
    }

    if extra:
        data.update(extra)

    summary_file.parent.mkdir(parents=True, exist_ok=True)

    with summary_file.open("w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, ensure_ascii=False)


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "自动识别单端/双端、自动检测或建立 Salmon index、"
            "批量 quant 并汇总 TPM/Count/FPKM。"
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "-i",
        "--input-dir",
        type=Path,
        default=None,
        help="FASTQ 输入目录（--qc-only 时可省略，届时从 quant 目录自举样本清单）",
    )

    parser.add_argument(
        "-o",
        "--output-dir",
        required=True,
        type=Path,
        help="结果输出目录",
    )

    parser.add_argument(
        "-r",
        "--reference",
        type=Path,
        default=None,
        help=(
            "transcript/CDS FASTA。"
            "若找不到已有 Salmon index，则必须提供。"
        ),
    )

    parser.add_argument(
        "--index",
        type=Path,
        default=None,
        help="已有 Salmon index 目录；不提供时自动搜索。",
    )

    parser.add_argument(
        "--salmon",
        default="salmon",
        help="Salmon 可执行程序名称或完整路径",
    )

    parser.add_argument(
        "--threads-per-sample",
        type=int,
        default=8,
        help="每个 salmon quant 使用的线程数",
    )

    parser.add_argument(
        "-j",
        "--jobs",
        type=int,
        default=0,
        help=(
            "同时运行的样本数。0=根据 CPU 自动计算，"
            "自动保证 jobs × threads_per_sample 不明显超过 CPU 核数。"
        ),
    )

    parser.add_argument(
        "--index-threads",
        type=int,
        default=0,
        help="建立 index 使用的线程数；0=使用全部可用 CPU",
    )

    parser.add_argument(
        "--fld-mean",
        type=float,
        default=200.0,
        help="单端数据 fragment length mean",
    )

    parser.add_argument(
        "--fld-sd",
        type=float,
        default=80.0,
        help="单端数据 fragment length standard deviation",
    )

    parser.add_argument(
        "--recursive",
        action="store_true",
        help="递归扫描输入目录中的子目录",
    )

    parser.add_argument(
        "--force-index",
        action="store_true",
        help="强制重新建立 Salmon index",
    )

    parser.add_argument(
        "--force-quant",
        action="store_true",
        help="强制重新运行已有样本的 quant",
    )

    parser.add_argument(
        "--no-validate-mappings",
        action="store_true",
        help="关闭 Salmon --validateMappings",
    )

    parser.add_argument(
        "--require-paired",
        action="store_true",
        help="严格模式：只定量 R1/R2 完整的双端样本，缺一端的样本跳过"
        "（不参与定量，名单写入 metadata/missing_read_samples.tsv）",
    )

    # ---- 定量后 QC ----
    qc_group = parser.add_argument_group(
        "定量后 QC（默认开启，未通过门槛的样本从主矩阵剔除）"
    )

    qc_group.add_argument(
        "--min-library-size",
        type=int,
        default=DEFAULT_MIN_LIBRARY_SIZE,
        help="硬门槛：单样本 ΣNumReads 低于此值判 FAIL 并从主矩阵剔除",
    )

    qc_group.add_argument(
        "--warn-library-size",
        type=int,
        default=DEFAULT_WARN_LIBRARY_SIZE,
        help="软门槛：单样本 ΣNumReads 低于此值只告警（WARN），不剔除",
    )

    qc_group.add_argument(
        "--min-detected-genes",
        type=int,
        default=DEFAULT_MIN_DETECTED_GENES,
        help="硬门槛：NumReads>=1 的基因数低于此值判 FAIL",
    )

    qc_group.add_argument(
        "--warn-detected-genes",
        type=int,
        default=DEFAULT_WARN_DETECTED_GENES,
        help="软门槛：检测到的基因数低于此值只告警",
    )

    qc_group.add_argument(
        "--min-mapping-rate",
        type=float,
        default=DEFAULT_MIN_MAPPING_RATE,
        help="硬门槛：percent_mapped（salmon meta_info）低于此值判 FAIL",
    )

    qc_group.add_argument(
        "--warn-mapping-rate",
        type=float,
        default=DEFAULT_WARN_MAPPING_RATE,
        help="软门槛：percent_mapped 低于此值只告警",
    )

    qc_group.add_argument(
        "--no-quant-qc",
        action="store_true",
        help="关闭定量后 QC（不生成报告、不做任何剔除）",
    )

    qc_group.add_argument(
        "--qc-report-only",
        action="store_true",
        help="只出 QC 报告，不按 QC 结果剔除样本（主矩阵仍含全部样本）",
    )

    qc_group.add_argument(
        "--qc-only",
        action="store_true",
        help="跳过 index/quant，直接对已有 quant 目录做 QC + 矩阵汇总"
        "（用于给历史结果补做质控）",
    )

    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="输出更多调试信息",
    )

    return parser.parse_args()


def main():
    args = parse_args()

    input_dir = (
        args.input_dir.expanduser().resolve()
        if args.input_dir is not None
        else None
    )
    output_dir = args.output_dir.expanduser().resolve()
    reference = (
        args.reference.expanduser().resolve()
        if args.reference is not None
        else None
    )
    explicit_index = (
        args.index.expanduser().resolve()
        if args.index is not None
        else None
    )

    if args.qc_only:
        # --qc-only：不读 FASTQ，样本清单从 quant 目录自举
        if not (output_dir / "quant").is_dir():
            print(
                f"[ERROR] --qc-only 需要已存在的 quant 目录：{output_dir / 'quant'}",
                file=sys.stderr,
            )
            sys.exit(1)

    else:
        if input_dir is None:
            print(
                "[ERROR] 缺少 -i/--input-dir（仅 --qc-only 时可省略）。",
                file=sys.stderr,
            )
            sys.exit(1)

        if not input_dir.is_dir():
            print(
                f"[ERROR] 输入目录不存在：{input_dir}",
                file=sys.stderr,
            )
            sys.exit(1)

    output_dir.mkdir(parents=True, exist_ok=True)

    log_root = output_dir / "logs"
    quant_root = output_dir / "quant"
    matrix_root = output_dir / "matrix"
    metadata_root = output_dir / "metadata"

    setup_logging(
        output_dir / "workflow.log",
        verbose=args.verbose,
    )

    try:
        if args.qc_only:
            # --qc-only 不需要 salmon 可执行程序
            salmon_bin = args.salmon
            salmon_ver = "skipped (--qc-only)"
        else:
            salmon_bin = salmon_exists(args.salmon)
            salmon_ver = salmon_version(salmon_bin)

            logging.info("Salmon：%s", salmon_bin)
            logging.info("版本：%s", salmon_ver)

        # 1. 检测样本
        if args.qc_only:
            samples = samples_from_quant_dir(quant_root)
            logging.info(
                "--qc-only：从 quant 目录重建样本清单，共 %d 个：%s",
                len(samples),
                quant_root,
            )
        else:
            logging.info("扫描 FASTQ：%s", input_dir)
            samples = detect_samples(
                input_dir,
                recursive=args.recursive,
            )

        pe_count = sum(1 for x in samples if x.layout == "PE")
        se_count = sum(1 for x in samples if x.layout == "SE")
        missing_count = sum(1 for x in samples if x.missing_read is not None)

        logging.info(
            "共识别 %d 个样本：PE=%d，SE=%d（其中 %d 个缺一端已标记）。",
            len(samples),
            pe_count,
            se_count,
            missing_count,
        )

        # 缺端样本报告（无论是否 --require-paired 都写）
        missing_samples = [s for s in samples if s.missing_read is not None]
        if missing_samples:
            missing_path = metadata_root / "missing_read_samples.tsv"
            missing_path.parent.mkdir(parents=True, exist_ok=True)
            with missing_path.open("w", encoding="utf-8", newline="") as f:
                writer = csv.writer(f, delimiter="\t")
                writer.writerow(["Sample", "Missing", "PresentReads"])
                for s in missing_samples:
                    reads = ",".join(str(p) for p in s.single)
                    writer.writerow([s.name, s.missing_read, reads])
            logging.warning(
                "%d 个样本缺一端，名单已写入：%s",
                len(missing_samples),
                missing_path,
            )

        # --require-paired：剔除缺端样本
        if args.require_paired:
            skipped = [s for s in samples if s.missing_read is not None]
            if skipped:
                logging.warning(
                    "--require-paired 生效：跳过 %d 个缺端样本。",
                    len(skipped),
                )
            samples = [s for s in samples if s.missing_read is None]

        if not samples:
            raise RuntimeError(
                "没有可定量的样本（--require-paired 过滤后为空，"
                "或输入目录无可识别样本）。"
            )

        if not args.qc_only:
            for sample in samples:
                if sample.layout == "PE":
                    logging.info(
                        "样本 %-30s PE  R1=%d  R2=%d",
                        sample.name,
                        len(sample.read1),
                        len(sample.read2),
                    )
                else:
                    logging.info(
                        "样本 %-30s SE  FASTQ=%d",
                        sample.name,
                        len(sample.single),
                    )

        # --qc-only 时 sample_info.tsv 来自 quant 目录重建（无 FASTQ 路径），
        # 不覆盖历史产物，改写到 sample_info_from_quant.tsv
        write_sample_info(
            samples,
            metadata_root
            / ("sample_info_from_quant.tsv" if args.qc_only else "sample_info.tsv"),
        )

        # 2-4. 检测 index / 执行 salmon quant（--qc-only 时整段跳过）
        if args.qc_only:
            logging.info(
                "--qc-only：跳过 index 建立与 salmon quant，"
                "直接对已有 quant 目录做定量 QC 与矩阵汇总。"
            )
            index_dir = output_dir / "salmon_index"
            failures = []
            jobs = args.jobs if args.jobs > 0 else 1
            threads_per_sample = args.threads_per_sample
            quant_root.mkdir(parents=True, exist_ok=True)
        else:
            # 2. 自动检测 / 建立 index
            found_index = None

            if not args.force_index:
                found_index = find_existing_index(
                    explicit_index=explicit_index,
                    output_dir=output_dir,
                    input_dir=input_dir,
                    reference=reference,
                )

            if found_index is not None:
                index_dir = found_index
                logging.info(
                    "使用已有 Salmon index：%s",
                    index_dir,
                )

            else:
                if reference is None:
                    raise RuntimeError(
                        "没有检测到有效 Salmon index，且未提供参考 FASTA。\n"
                        "请增加：-r /path/to/cds.fa\n"
                        "或者：--index /path/to/salmon_index"
                    )

                index_dir = (
                    explicit_index
                    if explicit_index is not None
                    else output_dir / "salmon_index"
                )

                total_cpu = os.cpu_count() or 1

                index_threads = (
                    args.index_threads
                    if args.index_threads > 0
                    else total_cpu
                )

                index_dir = build_index(
                    salmon_bin=salmon_bin,
                    reference=reference,
                    index_dir=index_dir,
                    threads=index_threads,
                    force=args.force_index,
                )

            # 3. 自动计算并行参数
            total_cpu = os.cpu_count() or 1

            if args.threads_per_sample < 1:
                raise RuntimeError("--threads-per-sample 必须 >= 1")

            threads_per_sample = args.threads_per_sample

            if args.jobs > 0:
                jobs = args.jobs
            else:
                jobs = max(
                    1,
                    total_cpu // threads_per_sample,
                )
                jobs = min(jobs, len(samples))

            logging.info(
                "CPU=%d；并行样本数=%d；每样本线程=%d；"
                "最大理论线程=%d。",
                total_cpu,
                jobs,
                threads_per_sample,
                jobs * threads_per_sample,
            )

            quant_root.mkdir(parents=True, exist_ok=True)
            log_root.mkdir(parents=True, exist_ok=True)

            # 4. Salmon quant
            failures = []
            done_count = 0
            pbar = make_pbar(len(samples), "salmon quant")

            with ThreadPoolExecutor(max_workers=jobs) as executor:
                future_map = {
                    executor.submit(
                        run_quant,
                        salmon_bin,
                        index_dir,
                        sample,
                        quant_root,
                        log_root,
                        threads_per_sample,
                        args.fld_mean,
                        args.fld_sd,
                        args.force_quant,
                        not args.no_validate_mappings,
                    ): sample
                    for sample in samples
                }

                for future in as_completed(future_map):
                    sample = future_map[future]

                    try:
                        name, status = future.result()
                        logging.info(
                            "[%s] 状态：%s",
                            name,
                            status,
                        )
                    except Exception as exc:
                        failures.append(
                            (sample.name, str(exc))
                        )
                        logging.error(
                            "[%s] 失败：%s",
                            sample.name,
                            exc,
                        )

                    if pbar is not None:
                        pbar.update(1)
                        done_count += 1
                        pbar.set_postfix(
                            done=done_count,
                            fail=len(failures),
                        )

            if pbar is not None:
                pbar.close()

            if failures:
                logging.error(
                    "共有 %d 个样本定量失败。",
                    len(failures),
                )

                failure_file = metadata_root / "failed_samples.tsv"
                failure_file.parent.mkdir(
                    parents=True,
                    exist_ok=True,
                )

                with failure_file.open(
                    "w",
                    encoding="utf-8",
                    newline="",
                ) as f:
                    writer = csv.writer(f, delimiter="\t")
                    writer.writerow(["Sample", "Error"])

                    for name, err in failures:
                        writer.writerow([name, err])

                logging.warning(
                    "失败样本已写入：%s",
                    failure_file,
                )

        # 5. 定量后 QC（补充质量步骤）+ 合并表达矩阵
        sample_filter = None

        if not args.no_quant_qc:
            metrics = collect_quant_metrics(
                quant_root=quant_root,
                samples=samples,
            )

            evaluated = evaluate_quant_qc(
                metrics,
                min_library_size=args.min_library_size,
                warn_library_size=args.warn_library_size,
                min_detected_genes=args.min_detected_genes,
                warn_detected_genes=args.warn_detected_genes,
                min_mapping_rate=args.min_mapping_rate,
                warn_mapping_rate=args.warn_mapping_rate,
            )

            qc_path, qc_passed, qc_failed = write_quant_qc_report(
                evaluated,
                metadata_root=metadata_root,
            )

            logging.info("定量 QC 报告：%s", qc_path)

            if args.qc_report_only:
                logging.warning(
                    "--qc-report-only 生效：不按 QC 结果剔除样本。"
                )
            else:
                sample_filter = set(qc_passed)
        else:
            logging.info("--no-quant-qc：跳过定量后 QC。")

        kept_samples = merge_quant_results(
            quant_root=quant_root,
            matrix_root=matrix_root,
            expected_samples=[s.name for s in samples],
            sample_filter=sample_filter,
            write_all_samples=True,
        )

        # 6. 汇总信息
        qc_extra: Dict[str, object] = {
            "quant_qc_enabled": not args.no_quant_qc,
            "quant_qc_report_only": bool(args.qc_report_only),
            "matrix_samples": len(kept_samples),
            "quant_qc_thresholds": {
                "min_library_size": args.min_library_size,
                "warn_library_size": args.warn_library_size,
                "min_detected_genes": args.min_detected_genes,
                "warn_detected_genes": args.warn_detected_genes,
                "min_mapping_rate": args.min_mapping_rate,
                "warn_mapping_rate": args.warn_mapping_rate,
            },
        }

        if not args.no_quant_qc:
            qc_extra["quant_qc_pass"] = sum(
                1 for m in evaluated if m["QC_status"] != "FAIL"
            )
            qc_extra["quant_qc_fail"] = sum(
                1 for m in evaluated if m["QC_status"] == "FAIL"
            )
            qc_extra["quant_qc_warn"] = sum(
                1 for m in evaluated if m["QC_status"] == "WARN"
            )

        write_run_summary(
            metadata_root / "run_summary.json",
            salmon_ver=salmon_ver,
            input_dir=input_dir,
            output_dir=output_dir,
            index_dir=index_dir,
            reference=reference,
            samples=samples,
            jobs=jobs,
            threads_per_sample=threads_per_sample,
            extra=qc_extra,
        )

        logging.info("=" * 70)
        logging.info("Salmon RNA-seq workflow 完成。")
        logging.info("输出目录：%s", output_dir)
        logging.info("单样本 quant：%s", quant_root)
        logging.info(
            "TPM 矩阵：%s（%d 个样本）",
            matrix_root / "TPM_matrix.tsv",
            len(kept_samples),
        )
        logging.info(
            "Count 矩阵：%s",
            matrix_root / "NumReads_matrix.tsv",
        )
        logging.info(
            "FPKM 矩阵：%s",
            matrix_root / "FPKM_matrix.tsv",
        )
        if not args.no_quant_qc:
            logging.info(
                "定量 QC 报告：%s",
                metadata_root / "quant_qc.tsv",
            )
        logging.info(
            "样本信息：%s",
            metadata_root / "sample_info.tsv",
        )
        logging.info(
            "总日志：%s",
            output_dir / "workflow.log",
        )
        logging.info("=" * 70)

        if failures:
            logging.warning(
                "流程已完成，但存在失败样本。"
                "请检查 logs/ 和 metadata/failed_samples.tsv。"
            )
            sys.exit(2)

    except KeyboardInterrupt:
        logging.error("用户终止流程。")
        sys.exit(130)

    except Exception as exc:
        logging.exception("流程失败：%s", exc)
        sys.exit(1)


if __name__ == "__main__":
    main()
