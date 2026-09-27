#!/usr/bin/env python3
"""
seqharvester_bridge.py — NCBI 序列收获器（自研，等效 VirPhyKit SeqHarvester）
=============================================================================
Entrez 搜索 → 批量下载 → 元数据解析 → 保存。不依赖 VirPhyKit。
"""

import os, sys, re, io, logging, time, csv
from typing import Dict, List, Optional


from Bio import Entrez, SeqIO
Entrez.email = "mmvp@example.com"
Entrez.tool = "MMPV-PhyloPipeline"


class LogCollector:
    def __init__(self, logger=None):
        self.logger = logger or logging.getLogger(__name__)
        self.messages = []
    def emit(self, msg):
        self.messages.append(msg); self.logger.info(msg)
    def warning(self, msg):
        self.messages.append(f"WARN: {msg}"); self.logger.warning(msg)  # 对齐 virphy_bridge 版
    def error(self, msg):
        self.messages.append(f"ERROR: {msg}"); self.logger.error(msg)


# ═══════════════════════════════════════════════════════════════════
# 以下函数直接提取自 VirPhyKit AnalysisWorker 类方法
# (src/SeqHarvester/function_SeqHarvester.py, lines 189-400+)
# ====================================================================
# 原类: class AnalysisWorker(QThread)
# 提取方法: search_sequences, fetch_records, parse_records
# 改动: self.virus_name → virus_name 参数
#        self.status_update.emit() → print/log
#        self.progress_update.emit() → print
# ═══════════════════════════════════════════════════════════════════

def get_taxid(accession=None, species=None):
    """病毒名/accession → NCBI taxonomy ID (taxid)

    方法1: 从 accession efetch GBK 的 source feature 解析 db_xref=taxon:xxx
    方法2: 从学名 esearch taxonomy 数据库
    """
    if accession:
        try:
            h = Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="text")
            for rec in SeqIO.parse(h, "gb"):
                for feat in rec.features:
                    if feat.type == "source":
                        for xref in feat.qualifiers.get("db_xref", []):
                            if xref.startswith("taxon:"):
                                return xref.split(":")[1]
                break
        except Exception:
            pass
    if species:
        try:
            h = Entrez.esearch(db="taxonomy", term=species)
            r = Entrez.read(h)
            if r.get("IdList"):
                return r["IdList"][0]
        except Exception:
            pass
    return None


def _search_sequences(virus_name: str, max_total: int = 0, taxid: str = None,
                      full_length: bool = False) -> List[str]:
    """检索: taxid 优先 (最准), 学名 fallback; 可选全长过滤。

    taxid:    txid{taxid}[Organism] — 精确到 taxonomy 节点, 避免学名歧义
    full_length: AND (complete genome/cds/sequence) — 过滤片段/部分序列
    """
    ids = []
    retmax = 10000
    retstart = 0
    if taxid:
        term = f'txid{taxid}[Organism]'
    else:
        term = f'"{virus_name}"[Organism]'
    term += ' NOT patent[Title]'
    if full_length:
        term += ' AND (complete genome[Title] OR complete cds[Title] OR complete sequence[Title])'
    while True:
        handle = Entrez.esearch(
            db="nucleotide",
            term=term,
            retmax=retmax, retstart=retstart)
        record = Entrez.read(handle)
        ids.extend(record["IdList"])
        retstart += retmax
        if retstart >= int(record["Count"]):
            break
    if max_total > 0 and len(ids) > max_total:
        ids = ids[:max_total]
    return ids


def _fetch_records(ids: List[str], log=None) -> List[str]:
    """NCBI 批量下载 (自研: 小批 + 重试 + 限流抗断流)。"""
    records = []
    batch_size = 100          # 降批: 500 批在大响应时易 IncompleteRead
    total_b = (len(ids) + batch_size - 1) // batch_size
    for i in range(0, len(ids), batch_size):
        batch_ids = ids[i:i + batch_size]
        batch_records = None
        for attempt in range(3):          # 断流重试
            try:
                handle = Entrez.efetch(
                    db="nucleotide", id=",".join(batch_ids),
                    rettype="gb", retmode="text")
                data = handle.read()
                handle.close()
                if data.strip():
                    batch_records = data.split("\n//\n")
                    break
            except Exception as e:
                if log:
                    log.warning(f"  Batch {i//batch_size+1} 尝试{attempt+1}失败: {str(e)[:120]}")
                time.sleep(1.5 * (attempt + 1))
        if batch_records:
            records.extend([
                r + "\n//\n" if not r.endswith("\n//\n") else r
                for r in batch_records if r.strip()
            ])
            if log:
                log.emit(f"  Batch {i//batch_size+1}/{total_b} OK")
        time.sleep(0.35)      # NCBI 限流 (~3 req/s)
    return records


def _parse_single_record(record_str: str, seq_id: str) -> dict:
    """
    VirPhyKit 风格的记录解析，结合:
    - SeqGrouper 的 parse_single_record (function_group.py line 77-111)
    - SeqHarvester 的 parse_records 分类逻辑 (function_SeqHarvester.py line 249+)
    """
    from datetime import datetime

    record = SeqIO.read(io.StringIO(record_str), "genbank")

    entry = {
        "Isolate": "N/A",
        "ID": record.id,
        "Organism": record.annotations.get("organism", "N/A"),
        "Length": len(record.seq),
        "Host": "N/A",
        "Geo Location": "N/A",
        "Collection Date": "N/A",
        "accession": seq_id,
    }

    cds_features = [f for f in record.features if f.type == "CDS"]
    desc_l = record.description.lower()
    # 完整性：标题关键词 + 排除 partial (与 VirPhyKit 原版一致)
    is_complete = (("complete genome" in desc_l or
                    "complete sequence" in desc_l or
                    "whole genome" in desc_l) and
                   "partial" not in desc_l)

    for feature in record.features:
        if feature.type == "source":
            qualifiers = feature.qualifiers
            # 位置 (VirPhyKit 原版逻辑: geo_loc_name > country)
            if "geo_loc_name" in qualifiers:
                geo = qualifiers["geo_loc_name"][0]
                entry["Geo Location"] = geo
            elif "country" in qualifiers:
                geo = qualifiers["country"][0]
                entry["Geo Location"] = geo
            # 宿主
            if "host" in qualifiers:
                entry["Host"] = qualifiers["host"][0]
            # 日期 (VirPhyKit 原版多格式解析)
            if "collection_date" in qualifiers:
                date_str = qualifiers["collection_date"][0]
                for fmt in ("%d-%b-%Y", "%d-%B-%Y", "%b-%Y", "%B-%Y",
                            "%Y-%m-%d", "%Y-%m", "%Y"):
                    try:
                        dt = datetime.strptime(date_str, fmt)
                        entry["Collection Date"] = dt.strftime("%Y-%m-%d")
                        break
                    except ValueError:
                        entry["Collection Date"] = date_str
            if "isolate" in qualifiers:
                entry["Isolate"] = qualifiers["isolate"][0]
            break

    entry["IsComplete"] = is_complete
    entry["CDS_Count"] = len(cds_features)

    return entry


def run_seqharvester(
    virus_name: str,
    output_dir: str,
    max_total: int = 0,
    taxid: str = None,
    full_length: bool = True,
    log: Optional[LogCollector] = None,
) -> Dict:
    """
    桥接 VirPhyKit SeqHarvester 全流程: 搜索 → 下载 → 解析 → 保存。

    taxid: NCBI taxonomy ID (优先用 txid[Organism] 检索, 最准)
    full_length: 只收完整基因组/序列 (过滤片段)
    """
    log = log or LogCollector()
    os.makedirs(output_dir, exist_ok=True)

    try:
        # Step 1: 搜索 (taxid 优先, 学名 fallback, 可选全长过滤)
        if taxid:
            log.emit(f"SeqHarvester: searching txid{taxid} (full_length={full_length})...")
        else:
            log.emit(f"SeqHarvester: searching '{virus_name}' (full_length={full_length})...")
        ids = _search_sequences(virus_name, max_total, taxid=taxid, full_length=full_length)
        log.emit(f"  Found {len(ids)} records")
        if not ids:
            return {"success": False, "error": "No sequences found", "n_sequences": 0}

        # Step 2: 下载 (VirPhyKit fetch_records)
        log.emit(f"  Downloading {len(ids)} records...")
        records = _fetch_records(ids, log)
        log.emit(f"  Fetched {len(records)} records")
        if not records:
            return {"success": False, "error": "No records fetched", "n_sequences": 0}

        # Step 3: 解析 (VirPhyKit parse_records 逻辑)
        log.emit(f"  Parsing records...")
        metadata_rows = []
        seq_records = []
        parsed = 0
        for idx, rs in enumerate(records):
            try:
                rec = SeqIO.read(io.StringIO(rs), "genbank")
                # 历史坑: batch 下载失败会丢弃整批 records → ids 与 records 错位,
                # accession 张冠李戴。record.id 才是真实来源。
                entry = _parse_single_record(rs, rec.id)
                seq_records.append(rec)
                metadata_rows.append(entry)
                parsed += 1
            except Exception:
                continue

        log.emit(f"  Parsed {parsed} records")

        # Step 4: 保存
        fasta_path = os.path.join(output_dir, "ncbi_sequences.fasta")
        meta_path = os.path.join(output_dir, "ncbi_metadata.csv")

        SeqIO.write(seq_records, fasta_path, 'fasta')

        with open(meta_path, 'w', newline='') as f:
            w = csv.DictWriter(f, fieldnames=[
                'ID', 'accession', 'Collection Date', 'Geo Location', 'Host',
                'Isolate', 'Organism', 'Length', 'IsComplete', 'CDS_Count'])
            w.writeheader()
            w.writerows(metadata_rows)

        n_dates = sum(1 for m in metadata_rows
                      if m['Collection Date'] and m['Collection Date'] != 'N/A')
        n_loc = sum(1 for m in metadata_rows
                    if m['Geo Location'] and m['Geo Location'] != 'N/A')
        n_host = sum(1 for m in metadata_rows
                     if m['Host'] and m['Host'] != 'N/A')

        log.emit(f"  Done: {parsed} seqs ({n_dates} dates, {n_loc} locations, {n_host} hosts)")

        return {
            "success": True,
            "n_sequences": parsed,
            "fasta_file": fasta_path,
            "metadata_csv": meta_path,
            "n_with_dates": n_dates,
            "n_with_location": n_loc,
            "n_with_host": n_host,
        }

    except Exception as e:
        log.error(f"SeqHarvester failed: {e}")
        return {"success": False, "error": str(e), "n_sequences": 0}
