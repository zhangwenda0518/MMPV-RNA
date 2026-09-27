#!/usr/bin/env python3
"""
data_collector.py — 从 MMPV-RNA 分析管线输出收集系统发育分析所需数据
====================================================================
四阶段的数据准备层：
  1. 从 06_extraction 收集各病毒的全长序列
  2. 关联 Global_Unified_Metadata_Core14 获取日期 + 地理 + 宿主
  3. 按病毒组织为 VirPhyKit 兼容的输入格式
  4. 输出 combined.fasta + dates.csv + sample_metadata.csv

去重 (2026-09-03): 本地 + 线上合并后按 (date, location, sequence) 三元组去重,
参考 geo_analysis.py smart subsample 的报告风格, 产物 data/dedup_report.csv。
"""


def dedup_triplet(fasta_path, meta_path, out_dir, logger=None):
    """
    按 (date, location, sequence-identical) 三元组去重。

    仅当两条样本 捕样日期相同 且 采样地点相同 且 序列完全一致 时,
    才视为克隆重复, 组内保留第一条 (或参考序列), 其余剔除。
    不同时间/地点的相同序列不剔除 (跨时空同序列 = 传播事件, 有流行病学信息)。

    Returns
    -------
    dict: {n_original, n_kept, n_removed, fasta, removed_csv}
          无重复时 fasta 指回原文件, 不写新文件。
    """
    import csv as _csv
    import os as _os
    from collections import OrderedDict
    from Bio import SeqIO

    log = logger
    recs = list(SeqIO.parse(fasta_path, "fasta"))
    if not recs:
        return {"n_original": 0, "n_kept": 0, "n_removed": 0,
                "fasta": fasta_path, "removed_csv": None}

    # 读元数据 (兼容 csv/tsv)
    delim = "\t" if str(meta_path).endswith((".tsv", ".txt")) else ","
    meta = {}
    try:
        with open(meta_path, encoding="utf-8", errors="replace") as f:
            for row in _csv.DictReader(f, delimiter=delim):
                name = (row.get("name") or row.get("sample") or row.get("id") or "").strip()
                if name:
                    meta[name] = {
                        "date": (row.get("date") or "").strip(),
                        "location": (row.get("location") or "").strip(),
                    }
    except OSError as e:
        if log:
            log.warning(f"  [dedup] 元数据读取失败 ({e}), 仅按序列一致性去重")

    # 三元组聚类
    groups = OrderedDict()
    for r in recs:
        sid = base_sample_id(r.id)   # 2026-09-15: 统一口径 (见 utils/seq_ids.py)
        m = meta.get(r.id) or meta.get(sid) or {}
        key = (m.get("date", ""), m.get("location", ""), str(r.seq).upper())
        groups.setdefault(key, []).append(r)

    kept, removed_rows = [], []
    for (dt, loc, _seq), grp in groups.items():
        # 参考序列 (description 含 reference) 优先保留
        grp_sorted = sorted(grp, key=lambda r: 0 if "reference" in (r.description or "") else 1)
        keep = grp_sorted[0]
        kept.append(keep)
        for r in grp_sorted[1:]:
            removed_rows.append({
                "removed_id": r.id, "kept_id": keep.id,
                "date": dt, "location": loc,
                "reason": "identical sequence, same date & location",
            })

    n_removed = len(removed_rows)
    if n_removed == 0:
        if log:
            log.info(f"  [dedup] 无克隆重复 (date+location+seq), 保留全部 {len(recs)} 条")
        return {"n_original": len(recs), "n_kept": len(recs), "n_removed": 0,
                "fasta": fasta_path, "removed_csv": None}

    # 写去重后 fasta + 报告 (原文件不动, 产 data/combined_dedup.fasta)
    _os.makedirs(out_dir, exist_ok=True)
    out_fasta = _os.path.join(out_dir, "combined_dedup.fasta")
    SeqIO.write(kept, out_fasta, "fasta")
    removed_csv = _os.path.join(out_dir, "dedup_report.csv")
    with open(removed_csv, "w", newline="", encoding="utf-8") as f:
        w = _csv.DictWriter(f, fieldnames=["removed_id", "kept_id", "date", "location", "reason"])
        w.writeheader(); w.writerows(removed_rows)
    if log:
        loc_stat = {}
        for row in removed_rows:
            loc_stat[row["location"] or "Unknown"] = loc_stat.get(row["location"] or "Unknown", 0) + 1
        top = ", ".join(f"{k}:{v}" for k, v in sorted(loc_stat.items(), key=lambda x: -x[1])[:5])
        log.info(f"  [dedup] 克隆重复去除: {len(recs)} → {len(kept)} (去除 {n_removed}; {top})")
        log.info(f"  [dedup] 明细: {removed_csv}")
    return {"n_original": len(recs), "n_kept": len(kept), "n_removed": n_removed,
            "fasta": out_fasta, "removed_csv": removed_csv}

import csv
import os
import re
import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from collections import defaultdict

from utils.seq_ids import base_sample_id
from utils.metadata_governance import has_host, normalize_host

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# 跨管线统一 I/O 布局: ⑤ per-virus 模块目录名随 MMPV_IO_LAYOUT 解析
try:
    from mmpv_common.io_layout import ph_dir
except ImportError:  # 脱离包环境直接运行时自举
    import os as _os, sys as _sys
    _REPO_ROOT = _os.path.dirname(_os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
    if _REPO_ROOT not in _sys.path:
        _sys.path.insert(0, _REPO_ROOT)
    from mmpv_common.io_layout import ph_dir


def load_sample_metadata(
    info_dir: str,
    metadata_csv: Optional[str] = None,
    logger: Optional[logging.Logger] = None,
) -> Dict[str, Dict[str, str]]:
    """
    从多种来源加载 Run → {date, location, host, tissue} 映射。

    优先级:
      1. metadata_csv (Global_Unified_Metadata_Core14.tsv) — 最完整
      2. SRA_GSA_Merged_Final.csv — fallback (ReleaseDate + Location)

    Returns
    -------
    dict: {run_id: {"date": "YYYY-MM-DD", "location": "...", "host": "...", "tissue": "..."}}
    """
    log = logger or logging.getLogger(__name__)
    meta: Dict[str, Dict[str, str]] = {}

    # ── Source 1: Global_Unified_Metadata_Core14 (full fields) ──
    if metadata_csv and os.path.exists(metadata_csv):
        delimiter = '\t' if str(metadata_csv).endswith('.tsv') else ','
        try:
            with open(metadata_csv, encoding='utf-8-sig') as f:
                reader = csv.DictReader(f, delimiter=delimiter)
                for row in reader:
                    run = row.get('Run', '').strip()
                    if not run:
                        continue
                    entry = {}
                    # Date: CollectionDate > ReleaseDate
                    cd = row.get('CollectionDate', '').strip()
                    rd = row.get('ReleaseDate', '').strip()
                    if cd and cd not in ('Not_Provided', '', 'None', 'NA'):
                        entry['date'] = cd
                    elif rd and rd not in ('Not_Provided', '', 'None', 'NA'):
                        entry['date'] = rd[:10] if len(rd) >= 10 else rd
                    # Location
                    loc = row.get('Location', '').strip()
                    if loc and loc not in ('Not_Provided', '', 'None', 'NA', 'missing'):
                        entry['location'] = loc
                    # Host
                    host = row.get('ScientificName', '').strip()
                    if host and host not in ('Not_Provided', '', 'None', 'NA'):
                        entry['host'] = host
                    # Tissue
                    tissue = row.get('Tissue', '').strip()
                    if tissue and tissue not in ('Not_Provided', '', 'None', 'NA'):
                        entry['tissue'] = tissue

                    if entry:
                        meta[run] = entry
            log.info(f"Loaded {len(meta)} sample metadata entries from Core14")
        except Exception as e:
            log.warning(f"Failed to load Core14 metadata: {e}")

    # ── Source 2: SRA_GSA_Merged_Final.csv (fallback) ──
    merged_csv = os.path.join(info_dir, "..", "search", "SRA_GSA_Merged_Final.csv")
    merged_csv = os.path.normpath(merged_csv)
    if os.path.exists(merged_csv):
        try:
            with open(merged_csv, encoding='utf-8-sig') as f:
                reader = csv.DictReader(f)
                for row in reader:
                    run = row.get('Run', '').strip()
                    if not run or run in meta:
                        continue
                    entry = {}
                    rel = row.get('ReleaseDate', '').strip()
                    if rel and rel not in ('missing', 'Not_Provided', ''):
                        entry['date'] = rel[:10] if len(rel) >= 10 else rel
                    loc = row.get('Location', '').strip()
                    if loc and loc not in ('missing', 'Not_Provided', ''):
                        entry['location'] = loc
                    if entry:
                        meta[run] = entry
            already = len(meta)
            log.info(f"Sample metadata: {already} total entries")
        except Exception as e:
            log.warning(f"Failed to load Merged_Final: {e}")

    return meta


def scan_extraction_dir(
    extract_dir: str,
    logger: Optional[logging.Logger] = None,
) -> Dict[str, Dict]:
    """
    扫描 06_extraction 目录，按病毒整理序列和样本。

    Returns
    -------
    { "Virus_Accession": {"dir": ..., "name": ..., "samples": [...], "n_samples": N}, ... }
    """
    log = logger or logging.getLogger(__name__)
    extract_path = Path(extract_dir)
    viruses: Dict[str, Dict] = {}

    if not extract_path.exists():
        log.error(f"Extraction directory not found: {extract_dir}")
        return viruses

    for virus_dir in sorted(extract_path.iterdir()):
        if not virus_dir.is_dir():
            continue

        acc_match = re.search(r'([A-Z]{1,2}_\d+\.\d+|[A-Z]{2}\d+\.\d+)$', virus_dir.name)
        virus_key = acc_match.group(1) if acc_match else virus_dir.name

        fasta_files = sorted(virus_dir.glob("*.full.fasta"))
        if not fasta_files:
            fasta_files = sorted(virus_dir.glob("*.fasta"))

        samples = []
        for fa in fasta_files:
            sample_id = fa.stem
            seq_len = 0
            n_ratio = 0.0
            try:
                for rec in SeqIO.parse(str(fa), "fasta"):
                    s = str(rec.seq).upper()
                    seq_len = len(s)
                    n_ratio = (s.count('N') / seq_len * 100) if seq_len > 0 else 0
                    break
            except Exception:
                pass
            samples.append({
                "path": str(fa), "sample_id": sample_id,
                "seq_len": seq_len, "n_ratio": n_ratio,
            })

        if samples:
            viruses[virus_key] = {
                "dir": str(virus_dir), "name": virus_dir.name,
                "samples": samples, "n_samples": len(samples),
            }

    log.info(f"Scanned extraction dir: {len(viruses)} viruses, "
             f"total {sum(v['n_samples'] for v in viruses.values())} samples")
    return viruses


def prepare_virus_inputs(
    virus_info: Dict,
    output_dir: str,
    sample_meta: Dict[str, Dict[str, str]],
    ref_fasta: Optional[str] = None,
    min_samples: int = 5,
    logger: Optional[logging.Logger] = None,
) -> Optional[Dict]:
    """
    为单个病毒准备 VirPhyKit 兼容的输入文件。

    产出:
      combined.fasta       — 所有样本全长序列 + 参考序列
      dates.csv            — name, date [, location]        时间通道
      sample_metadata.csv  — name, date, location           地理通道
      host.csv             — name, host, tissue             宿主通道 (**只含有宿主的样品**)

    2026-09-16 第三轮: 三条通道各管一个研究维度 —— host/tissue 从地理表移出,
    独立成 `host.csv`, 并只收 `has_host()` 为真的样品 (老师:
    「host 也独立出来, 仅对有 host 的进行分析」)。

    Returns
    -------
    dict or None
    """
    log = logger or logging.getLogger(__name__)
    os.makedirs(output_dir, exist_ok=True)
    data_dir = os.path.join(output_dir, ph_dir("data"))
    os.makedirs(data_dir, exist_ok=True)

    samples = virus_info["samples"]
    if len(samples) < min_samples:
        log.info(f"  [{virus_info['name']}] {len(samples)} samples < {min_samples}, skipping")
        return None

    # ── Match samples with metadata ──
    matched = []
    for s in samples:
        sid = s["sample_id"]
        # 规范化：去 accession 后缀 (CRR1126135.OR489165.1 → CRR1126135)
        base_sid = base_sample_id(sid)   # 2026-09-15: 统一口径
        info = sample_meta.get(sid) or sample_meta.get(base_sid)
        if info:
            matched.append({**s, **info})

    n_with_dates = sum(1 for s in matched if s.get('date'))
    n_with_loc = sum(1 for s in matched if s.get('location'))
    # ⚠ 必须用 `has_host()` 而不是 `bool(s.get('host'))` —— 否则日志报的
    #   "hosts=N" 会把 `Unknown` / `uncultured bacterium` 也算进去, 与宿主通道
    #   实际行数对不上 (实测: 报 5、通道里 3)。判据只能有一处实现。
    n_with_host = sum(1 for s in matched if has_host(s.get('host')))
    has_location = n_with_loc >= min_samples

    if n_with_dates < min_samples:
        log.info(f"  [{virus_info['name']}] {n_with_dates} samples with dates < {min_samples}, skipping")
        return None

    log.info(f"  [{virus_info['name']}] dates={n_with_dates} locations={n_with_loc} hosts={n_with_host} / {len(samples)}")

    # ── Write combined.fasta ──
    combined_fasta = os.path.join(data_dir, "combined.fasta")
    records = []
    for s in matched:
        try:
            for rec in SeqIO.parse(s["path"], "fasta"):
                rec.id = s["sample_id"]
                rec.description = ""
                records.append(rec)
                break
        except Exception as e:
            log.warning(f"    Failed to read {s['path']}: {e}")

    if ref_fasta and os.path.exists(ref_fasta):
        try:
            for rec in SeqIO.parse(ref_fasta, "fasta"):
                rec.id = rec.id.split('.')[0] if '.' in rec.id else rec.id
                rec.description = "reference"
                records.append(rec)
        except Exception as e:
            log.warning(f"    Failed to read reference {ref_fasta}: {e}")

    SeqIO.write(records, combined_fasta, "fasta")
    log.info(f"    combined.fasta: {len(records)} sequences")

    # ── Write dates.csv (with location if available) ──
    dates_csv = os.path.join(data_dir, "dates.csv")
    with open(dates_csv, 'w', newline='') as f:
        writer = csv.writer(f)
        if has_location:
            writer.writerow(["name", "date", "location"])
            for s in matched:
                writer.writerow([s["sample_id"], s.get('date', ''), s.get('location', '')])
        else:
            writer.writerow(["name", "date"])
            for s in matched:
                writer.writerow([s["sample_id"], s.get('date', '')])

    # ── Write sample_metadata.csv (for geo analysis, VirSpaceTime) ──
    #
    # 2026-09-16 第三轮 (老师: 「host 也独立出来, 仅对有 host 的进行分析」):
    # 地理通道**只留地理** (`name, date, location`), host/tissue 挪到宿主通道。
    # 这样三条通道各管一个研究维度, 不再出现"同一列两处各写一份"的漂移风险:
    #
    #   dates.csv            name, date[, location]  → 时间通道 (定年)
    #   sample_metadata.csv  name, date, location    → 地理通道 (系统地理)
    #   host.csv             name, host, tissue      → 宿主通道 (宿主分化)
    meta_csv = os.path.join(data_dir, "sample_metadata.csv")
    with open(meta_csv, 'w', newline='') as f:
        fields = ["name", "date", "location"]
        writer = csv.DictWriter(f, fieldnames=fields, extrasaction='ignore')
        writer.writeheader()
        for s in matched:
            writer.writerow({
                "name": s["sample_id"],
                "date": s.get('date', ''),
                "location": s.get('location', ''),
            })

    # ── Write host.csv (宿主通道) —— **只收有宿主的样品** ──
    # 「有没有宿主」用 `has_host()` 统一判据 (未填 / `uncultured bacterium` /
    # `metagenome` 这类无信息宿主一律不算) → 下游宿主分析天然只跑在这个子集上,
    # 不需要各分析脚本各自再过滤一遍 (重复判据 = 必然漂移)。
    host_csv = os.path.join(data_dir, "host.csv")
    host_rows = [s for s in matched if has_host(s.get('host'))]
    with open(host_csv, 'w', newline='') as f:
        hfields = ["name", "host", "tissue"]
        writer = csv.DictWriter(f, fieldnames=hfields, extrasaction='ignore')
        writer.writeheader()
        for s in host_rows:
            writer.writerow({
                "name": s["sample_id"],
                "host": normalize_host(s.get('host')),
                "tissue": str(s.get('tissue', '') or '').strip(),
            })
    n_host_rows = len(host_rows)

    avg_seq_len = int(sum(s["seq_len"] for s in matched) / len(matched)) if matched else 0

    return {
        "virus_name": virus_info["name"],
        "virus_key": virus_info.get("dir", ""),
        "combined_fasta": combined_fasta,
        "dates_csv": dates_csv,
        "meta_csv": meta_csv,
        "host_csv": host_csv,
        "n_host_rows": n_host_rows,
        "n_samples": len(samples),
        "n_with_dates": n_with_dates,
        "n_with_location": n_with_loc,
        "n_with_host": n_with_host,
        "has_location": has_location,
        "has_host": n_host_rows > 0,
        "avg_seq_len": avg_seq_len,
        "output_dir": output_dir,
    }


def find_ref_fasta(virus_key: str, variants_dir: Optional[str] = None) -> Optional[str]:
    """查找病毒的参考序列 FASTA"""
    if not variants_dir:
        return None
    candidates = [
        os.path.join(variants_dir, "virus-fasta", f"ref_{virus_key}", f"ref_{virus_key}.ref.fasta"),
        os.path.join(variants_dir, f"ref_{virus_key}.fasta"),
    ]
    base_key = virus_key.split('.')[0] if '.' in virus_key else virus_key
    for root, dirs, files in os.walk(variants_dir):
        for f in files:
            if base_key in f and f.endswith('.fasta') and 'ref' in f.lower():
                candidates.append(os.path.join(root, f))
    for c in candidates:
        if os.path.exists(c):
            return c
    return None
