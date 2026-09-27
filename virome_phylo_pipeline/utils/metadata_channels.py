#!/usr/bin/env python3
"""metadata_channels.py — 把 metadata 拆成「时间 / 地理 / 宿主」三张通道表

背景 (2026-09-16 要求「在元数据清洗检查后, 定年和地理分开去研究」)
=================================================================
`utils/metadata_governance.py` 产出 `metadata_std.csv` (日期结构化 + 地理分层归一),
但下游真正消费的是三张分工不同的表:

    data/dates.csv            name, date [, location]      → 时间通道
                              clock (TreeTime-RTT/TreeDater) / beast / gene_dating
    data/sample_metadata.csv  name, date, location         → 地理通道
                              phylogeo / geo / popgen
    data/host.csv             name, host, tissue           → 宿主通道
                              host (宿主分化) / 单倍型着色

这三张表由 `utils/data_collector.prepare_virus_inputs()` 在 **prep 阶段**写出,
而 metadata 治理发生在 prep **之后** → 治理结果此前无法回流 (断链)。
本模块负责把治理后的值**按样品名回填**到这三张表, 使三条分析线各吃各的表。

宿主通道独立 (2026-09-16 第三轮, 老师: 「host 也独立出来, 仅对有 host 的进行分析」)
====================================================================================
`host` 原本只是地理表的附属列 → 与 `location` 混在一张表里, 无法表达
「只对有宿主的样品做宿主分析」。现在:

* **独立成 `data/host.csv`**, 地理表只留地理 (`name, date, location`)。
* **只收 `has_host()` 为真的样品** —— 未填 (`Unknown`/`NA`) 与**无信息宿主**
  (`uncultured bacterium` / `metagenome` / `synthetic construct`) 都不进表。
  于是下游宿主分析天然跑在这个子集上, 不必各脚本再各写一遍过滤 (重复判据必漂移)。
* 判据唯一来源: `utils.metadata_governance.has_host()`。

回填规则 (保守: 只改治理确实提供了内容的字段)
================================================
* **骨架** = 原表 `name` 列及其行顺序。**不增不删行** —— 行的增减会改变比对输入与
  checkpoint 语义, 属另一件事 (那由治理阶段的 `--drop-no-date/-location` 决定)。
  唯一例外: 宿主通道的**行集本身就是过滤器** (只收有宿主的样品), 这是设计意图。
* `date`     —— 治理表有日期列就覆盖 (含「占位符被清空 → 该行变空」, 这正是要修的)。
* `location` —— 治理表有地名列就覆盖。
* `host` / `tissue` —— 取值优先级 (宿主值常常**不在**治理范围内, 来自
  public_metadata_pipeline 的 info_dir / Core14):
    ① 治理来源表有 host 列  → 用它 (权威)
    ② 现有 host.csv 有 host  → 保留原值 (迁移期自保持, 保证幂等)
    ③ 旧版地理表仍带 host 列 → 用它 (迁移前写出的 work_dir 也能平滑升级)
  三级都取不到 → 宿主通道为空表, 并在日志里说明 (留痕, 不静默)。
* 治理表里**找不到**的样品名 → 保留原值, 计入 `unmatched` (留痕, 不静默)。
* 覆盖前把原表另存 `<表名>.bak_pre_governance` (**只在不存在时存一次**,
  避免第二次运行把「已治理值」当成原件备份掉)。
* 输出 `data/metadata/channels.json` 记录来源/指纹/计数 —— 便于机械核对。

样品名匹配口径
================
与 `prepare_virus_inputs()` 保持一致: 先原样命中, 再退到 `base_sample_id()`
(去 accession 版本号等)。若两个不同原名归一到同一键, 记为 `collisions` 并让
**首个命中**生效 —— 冲突可见, 不静默。
"""
import argparse
import csv
import hashlib
import json
import logging
import os
import sys
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple

# 脚本方式运行 (python utils/metadata_channels.py) 时 sys.path[0] 是 utils/,
# 而 `from utils.X import ...` 需要管线根目录在路径上 → 与 import_export.py 同款 bootstrap。
_ROOT = str(Path(__file__).resolve().parent.parent)
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

from utils.metadata_governance import (  # noqa: E402
    DATE_HINTS, LOC_HINTS, NAME_HINTS, has_host, normalize_host,
)
from utils.seq_ids import base_sample_id  # noqa: E402

# host / tissue 的列名提示 (治理模块不管这两列, 本处自行识别)
#
# ⚠ `scientificname` / `scientific_name` 必须在内, 且放**末位**:
#   本项目 `Global_Unified_Metadata_Core14.tsv` 的宿主列就叫 `ScientificName`
#   (见 `data_collector.load_sample_metadata`: `host = row.get('ScientificName')`)。
#   此前它不在 hints 里 → 只认 prep 直接写出的 host.csv; 一旦改成"以 Core14 为
#   宿主来源"(批处理/--virus 的常见形态), 宿主列识别不到 → 宿主通道空 → fail-safe
#   跳过 → 宿主分析静默不跑。放末位保证真有 `host` 列时仍优先精确命中 `host`。
HOST_HINTS = ("host", "host_species", "host_scientific_name", "isolation_host",
              "scientificname", "scientific_name")
TISSUE_HINTS = ("tissue", "isolation_source", "source", "organ")

DATES_FIELDS = ("name", "date", "location")
META_FIELDS = ("name", "date", "location")       # 地理通道: 只留地理
HOST_FIELDS = ("name", "host", "tissue")         # 宿主通道: 只收有宿主的样品

# 通道布局版本 —— **必须**参与 checkpoint 判据。
#
# 为什么需要它: `needs_rebuild()` 只看「来源指纹 + 产物是否在」。当通道**布局**
# 变了 (2026-09-16: 2 张 → 3 张, `host` 从地理表迁出), 旧 work_dir 里的
# `channels.json` 既没记 `out_host`、来源文件又没变 → 判定"无需重建" →
# **宿主通道永远不会被拆出来**, 下游主机分析静默不跑。这不属于"来源变了",
# 属于"口径变了", 必须单独一个版本号来戳破。
CHANNELS_SCHEMA = 3   # 2 = 时间+地理; 3 = 时间+地理+宿主 (2026-09-16 第三轮)


# ══════════════════════════════════════════════════════════════════
# 读表 / 列识别
# ══════════════════════════════════════════════════════════════════
def read_table(path: str) -> Tuple[List[str], List[Dict[str, str]]]:
    """读 CSV/TSV → (列名, 行)。带 BOM 容忍、坏字符替换 (与治理模块同口径)。"""
    with open(path, encoding="utf-8-sig", errors="replace", newline="") as f:
        head = f.readline()
        f.seek(0)
        delim = "\t" if head.count("\t") > head.count(",") else ","
        reader = csv.DictReader(f, delimiter=delim)
        cols = [c for c in (reader.fieldnames or []) if c]
        rows = [r for r in reader]
    return cols, rows


def _pick(cols: List[str], hints, taken: set) -> Optional[str]:
    """精确 > 前缀 > 包含 (与 metadata_governance._pick_col 同规则)。"""
    low = [(c, c.lower().strip()) for c in cols]
    for h in hints:
        for c, lc in low:
            if lc == h and c not in taken:
                return c
    for h in hints:
        for c, lc in low:
            if lc.startswith(h) and c not in taken:
                return c
    for h in hints:
        for c, lc in low:
            if h in lc and c not in taken:
                return c
    return None


def pick_columns(cols: List[str]) -> Dict[str, Optional[str]]:
    """识别 name / date / location / host / tissue 列 (与治理模块同一套 hints)。"""
    taken: set = set()
    out: Dict[str, Optional[str]] = {}
    for key, hints in (("date", DATE_HINTS), ("location", LOC_HINTS),
                       ("name", NAME_HINTS), ("host", HOST_HINTS),
                       ("tissue", TISSUE_HINTS)):
        c = _pick(cols, hints, taken)
        out[key] = c
        if c:
            taken.add(c)
    return out


# ══════════════════════════════════════════════════════════════════
# 按名索引
# ══════════════════════════════════════════════════════════════════
def index_by_name(rows: List[Dict[str, str]], name_col: Optional[str]):
    """建 name → 行 的索引。原样名与 base_sample_id 归一化名都进索引。

    Returns (index, collisions)
      collisions: 归一化键被不同原名重复占用的次数 (首个命中生效)。
    """
    idx: Dict[str, Dict[str, str]] = {}
    collisions = 0
    for r in rows:
        raw = (r.get(name_col, "") if name_col else "") or ""
        raw = str(raw).strip()
        if not raw:
            continue
        for k in {raw, base_sample_id(raw)}:
            if not k:
                continue
            if k in idx:
                prev = str(idx[k].get(name_col, "")).strip()
                if prev != raw:
                    collisions += 1
                continue
            idx[k] = r
    return idx, collisions


def lookup(idx: Dict[str, Dict[str, str]], name: str) -> Optional[Dict[str, str]]:
    """先原样、再归一化命中。"""
    n = str(name).strip()
    return idx.get(n) or idx.get(base_sample_id(n))


# ══════════════════════════════════════════════════════════════════
# 回填
# ══════════════════════════════════════════════════════════════════
def _fingerprint(paths: List[str]) -> str:
    h = hashlib.sha1()
    for p in paths:
        try:
            st = os.stat(p)
            h.update(f"{os.path.basename(p)}|{st.st_size}|{int(st.st_mtime)}".encode())
        except OSError:
            h.update(f"{os.path.basename(p)}|missing".encode())
    return h.hexdigest()[:16]


def source_fingerprint(meta_source: str, geo_source: Optional[str] = None) -> str:
    """**只**对"值来源"取指纹。

    ⚠ 不能把骨架表 (dates.csv / sample_metadata.csv) 纳入指纹 —— 就地回填会改写
    骨架表自身的 mtime, 纳入就会每次判定"需要重建", 永不稳定。骨架的产物存在性
    另由 `needs_rebuild()` 检查。
    """
    ps = [meta_source]
    if geo_source and geo_source != meta_source:
        ps.append(geo_source)
    return _fingerprint(ps)


def _backup_once(path: str, log) -> Optional[str]:
    """覆盖前备份, 只在备份不存在时写 (避免把已治理值当原件备份)。"""
    bak = path + ".bak_pre_governance"
    if os.path.exists(path) and not os.path.exists(bak):
        with open(path, "rb") as src, open(bak, "wb") as dst:
            dst.write(src.read())
        log.info(f"  [channels] 原表已备份: {os.path.basename(bak)}")
        return bak
    return bak if os.path.exists(bak) else None


def _write_table(path: str, fields: Tuple[str, ...], rows: List[Dict[str, str]]):
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(fields), extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in fields})


def rebuild_channels(dates_csv: str,
                     meta_csv: Optional[str] = None,
                     *,
                     host_csv: Optional[str] = None,
                     meta_source: Optional[str] = None,
                     geo_source: Optional[str] = None,
                     out_dates: Optional[str] = None,
                     out_meta: Optional[str] = None,
                     out_host: Optional[str] = None,
                     backup: bool = True,
                     logger: Optional[logging.Logger] = None) -> Dict:
    """把治理后的 metadata 回填到「时间 / 地理 / 宿主」三张通道表。

    Parameters
    ----------
    dates_csv, meta_csv, host_csv : str
        现有三张通道表 (提供 name 骨架与行顺序)。`meta_csv` / `host_csv` 可省略。
    meta_source : str
        **日期值**的来源表 (通常是治理产物 `metadata_std.csv`)。给 None 则原样返回。
    geo_source : str
        **地理/宿主值**的来源表; 省略时同 `meta_source`。
    out_dates, out_meta, out_host : str
        输出路径; 省略时**就地覆盖** `dates_csv` / `meta_csv` / `host_csv`。

    Returns dict(success, dates_csv, meta_csv, host_csv, report_path, stats...)
    """
    log = logger or logging.getLogger(__name__)
    if not meta_source or not os.path.exists(meta_source):
        return {"success": False, "error": f"meta_source 不存在: {meta_source}"}
    if not dates_csv or not os.path.exists(dates_csv):
        return {"success": False, "error": f"骨架表不存在: {dates_csv}"}

    geo_source = geo_source or meta_source
    out_dates = out_dates or dates_csv
    out_meta = out_meta or meta_csv
    # 宿主通道输出: 显式给 > 已有 host.csv > **由地理表同目录派生** `host.csv`。
    # 派生是刻意的 —— 本模块的职责就是"一张表拆三条通道", 不给输出路径也应落地
    # (否则 host 会静默留在原地, 等于没独立出来)。
    if not out_host:
        out_host = host_csv or (os.path.join(os.path.dirname(out_meta), "host.csv")
                                if out_meta else None)

    # ── 骨架 ──
    d_cols, d_rows = read_table(dates_csv)
    d_pick = pick_columns(d_cols)
    if not d_pick["name"]:
        return {"success": False, "error": f"骨架表无 name 列: {dates_csv}"}
    m_rows: List[Dict[str, str]] = []
    m_cols: List[str] = []
    m_pick: Dict[str, Optional[str]] = {}
    if meta_csv and os.path.exists(meta_csv):
        m_cols, m_rows = read_table(meta_csv)
        m_pick = pick_columns(m_cols)
        if not m_pick["name"]:
            m_rows, m_cols, m_pick = [], [], {}

    # 迁移期: 旧版地理表 (name,date,location,host,tissue) 仍带 host 列
    legacy_host_rows = m_rows if m_pick.get("host") else []
    legacy_host_pick = m_pick if legacy_host_rows else {}

    h_cols: List[str] = []
    h_rows: List[Dict[str, str]] = []
    h_pick: Dict[str, Optional[str]] = {}
    h_idx: Dict[str, Dict[str, str]] = {}
    if host_csv and os.path.exists(host_csv):
        h_cols, h_rows = read_table(host_csv)
        h_pick = pick_columns(h_cols)
        if not h_pick["name"]:
            h_cols, h_rows, h_pick = [], [], {}
        else:
            h_idx, _ = index_by_name(h_rows, h_pick["name"])

    # ── 值来源 ──
    s_cols, s_rows = read_table(meta_source)
    s_pick = pick_columns(s_cols)
    idx, collisions = index_by_name(s_rows, s_pick["name"])
    if geo_source != meta_source and os.path.exists(geo_source):
        g_cols, g_rows = read_table(geo_source)
        g_pick = pick_columns(g_cols)
        gidx, gcol = index_by_name(g_rows, g_pick["name"])
        collisions += gcol
    else:
        g_pick, gidx = s_pick, idx

    if not s_pick["name"]:
        return {"success": False, "error": f"来源表无 name 列: {meta_source}"}

    # ── 宿主值来源: **有序取值链** (不是"选一个") ──
    #   ① 治理来源表有 host 列 —— 权威 (老师要求治理后的值必须生效)
    #   ② 现有 host.csv        —— 自保持。**必须留在链里**: 治理来源可能是
    #      部分表 (只含部分样品), 若只认 ①, 没被治理覆盖到的样品会**丢掉宿主**,
    #      宿主通道凭空变小 (selftest 实测: 3 样品 → 1 样品)。
    #   ③ 旧版地理表 host 列    —— 迁移期, 让治理前写出的 work_dir 平滑升级
    #   三级都取不到 → 该样品确无宿主信息 → 不进通道; 全链为空 → 不落地 (fail-safe)。
    host_chain: List[Tuple[Dict[str, Dict[str, str]], Dict[str, Optional[str]], str]] = []
    if g_pick.get("host"):
        host_chain.append((gidx, g_pick, "governance"))
    if s_pick.get("host") and s_pick is not g_pick:
        host_chain.append((idx, s_pick, "governance"))
    if h_pick.get("host"):
        host_chain.append((h_idx, h_pick, "host_channel"))
    if legacy_host_pick.get("host"):
        _lh_idx, _ = index_by_name(legacy_host_rows, legacy_host_pick["name"])
        host_chain.append((_lh_idx, legacy_host_pick, "legacy_geo_channel"))
    host_src_kind = host_chain[0][2] if host_chain else "none"

    # tissue 取值链: 与 host 同款思路 (逐级取第一个非空), 只是不要"可用性"判据
    tissue_chain: List[Tuple[Dict[str, Dict[str, str]], Dict[str, Optional[str]], str]] = []
    for _idx, _pick, _kind in ((gidx, g_pick, "governance"), (idx, s_pick, "governance"),
                               (h_idx, h_pick, "host_channel")):
        if _pick.get("tissue") and all(_idx is not c[0] or _pick is not c[1]
                                       for c in tissue_chain):
            tissue_chain.append((_idx, _pick, _kind))

    stats = {"n_skeleton": len(d_rows), "n_source": len(s_rows),
             "n_matched": 0, "n_unmatched": 0, "n_date_changed": 0,
             "n_loc_changed": 0, "n_host_changed": 0, "collisions": collisions,
             "n_host_rows": 0, "n_host_dropped": 0, "n_host_skeleton": 0,
             "host_src_kind": host_src_kind, "host_src_used": {},
             "src_date_col": s_pick["date"], "src_loc_col": s_pick["location"],
             "src_host_col": (host_chain[0][1].get("host") if host_chain else None),
             "src_tissue_col": s_pick["tissue"]}

    def _apply(row: Dict[str, str], src: Optional[Dict[str, str]],
               src_pick: Dict[str, Optional[str]], do_date: bool) -> Dict[str, str]:
        """按规则回填一行 (时间/地理通道)。src=None → 原样。"""
        out = dict(row)
        if src is None:
            return out
        if do_date and s_pick["date"]:
            new = (src.get(s_pick["date"], "") or "")
            new = "" if new is None else str(new).strip()
            if new != str(out.get("date", "") or "").strip():
                stats["n_date_changed"] += 1
            out["date"] = new
        if src_pick["location"]:
            new = str(src.get(src_pick["location"], "") or "").strip()
            if new != str(out.get("location", "") or "").strip():
                stats["n_loc_changed"] += 1
            out["location"] = new
        return out

    # ── 时间通道 (dates.csv): 只回填 date (+ 若原表本就有 location 列则一并) ──
    #
    # 2026-09-16 第三轮: 骨架若是一张**整表** (work_dir 的 `--meta` 常是
    # name,date,location,host,tissue 全字段表), 时间通道会把 host/tissue 一并继承
    # → 会出现**第三份** host (时间表 + 宿主表), 违背"一处一义"。
    # 故时间通道**剔掉属别的通道的列**: 时间表只留时间 (+ 地理, 兼容原形状)。
    # 其它未识别列保留 (避免把用户自己带的辅助列莫名删掉)。
    _OTHER_CHANNEL = {"host", "host_species", "host_scientific_name", "isolation_host",
                      "tissue", "isolation_source", "organ"}
    if d_cols:
        d_fields = tuple(c for c in d_cols if c.strip().lower() not in _OTHER_CHANNEL)
    else:
        d_fields = ("name", "date", "location")
    out_d: List[Dict[str, str]] = []
    for row in d_rows:
        nm = row.get(d_pick["name"], "")
        src = lookup(idx, nm)
        if src is None:
            stats["n_unmatched"] += 1
        else:
            stats["n_matched"] += 1
        out_d.append(_apply(row, src, s_pick, do_date=True))

    # ── 地理通道 (sample_metadata.csv): name / date / location ──
    out_m: List[Dict[str, str]] = []
    for row in (m_rows or []):
        nm = row.get(m_pick["name"], "")
        # 地理值优先地理专用来源, 缺失再退日期来源
        src = lookup(gidx, nm) or lookup(idx, nm)
        out_m.append(_apply(row, src, g_pick, do_date=True))

    # ── 宿主通道 (host.csv): **只收有宿主的样品** ──
    # 行集 = 地理骨架 (样品全集) 的顺序, 附加 host.csv 里独有的名字 (防御性, 不丢行)。
    host_skel_rows = m_rows or d_rows
    host_skel_pick = m_pick if m_rows else d_pick
    out_h: List[Dict[str, str]] = []
    seen: set = set()
    for row in host_skel_rows:
        nm = str(row.get(host_skel_pick["name"], "") or "").strip()
        if not nm or nm in seen:
            continue
        seen.add(nm)
        stats["n_host_skeleton"] += 1
        # host: 沿取值链逐级取第一个「真实可用」的值 (空/无信息宿主不算命中)
        hv, used = "", ""
        for _idx, _pick, _kind in host_chain:
            _r = lookup(_idx, nm)
            if not _r:
                continue
            _v = normalize_host(_r.get(_pick["host"]))
            if has_host(_v):
                hv, used = _v, _kind
                break
        if not hv:
            # 链上全无 → 退回骨架自身 (兼容调用方直接把旧表当骨架传进来的用法)
            _v = normalize_host(row.get("host", ""))
            if has_host(_v):
                hv, used = _v, "skeleton"
        # tissue 同样走取值链, 但**必须在本行内重新查表**。
        # ⚠ 别复用上面地理循环遗留的 `src`: 它指向**最后一行**的地理来源,
        #   会让所有样品都拿到末行的 tissue (实测 bug: S1/S2 全变成 S3 的 `stem`)。
        tv = ""
        for _idx, _pick, _kind in tissue_chain:
            _r = lookup(_idx, nm)
            if not _r:
                continue
            _v = str(_r.get(_pick["tissue"], "") or "").strip()
            if _v:
                tv = _v
                break
        if not tv:
            tv = str(row.get("tissue", "") or "").strip()
        if not has_host(hv):
            stats["n_host_dropped"] += 1
            continue
        stats["host_src_used"][used] = stats["host_src_used"].get(used, 0) + 1
        # 变更计数: 与**上一轮宿主通道**比 (而非与地理骨架比 —— 后者已无 host 列,
        # 每轮都会误计成"变了", 使计数器失去意义)
        _prev_row = h_idx.get(nm) or h_idx.get(base_sample_id(nm)) or {}
        prev = normalize_host(_prev_row.get("host", "") or row.get("host", ""))
        if prev != hv:
            stats["n_host_changed"] += 1
        out_h.append({"name": nm, "host": hv, "tissue": tv})
        stats["n_host_rows"] += 1
    # host.csv 里独有、地理骨架里没有的名字 → 追加 (不静默丢)
    for row in h_rows:
        nm = str(row.get(h_pick["name"], "") or "").strip()
        if not nm or nm in seen:
            continue
        hv = normalize_host(row.get("host", ""))
        if not has_host(hv):
            continue
        seen.add(nm)
        out_h.append({"name": nm, "host": hv,
                      "tissue": str(row.get("tissue", "") or "").strip()})
        stats["n_host_rows"] += 1

    if backup:
        _backup_once(out_dates, log)
        if out_meta and os.path.exists(out_meta):
            _backup_once(out_meta, log)
        if out_host and os.path.exists(out_host):
            _backup_once(out_host, log)

    _write_table(out_dates, d_fields, out_d)
    if out_meta and (m_rows or not os.path.exists(out_meta)):
        _write_table(out_meta, META_FIELDS, out_m)
    # ⚠ fail-safe: 只有当**确实拿到了 host 信息**时才写宿主通道。
    #   `host_src_kind == "none"` 意味着三处都找不到 host 列 (治理来源没有、
    #   调用方没把现有 host.csv 传进来、地理表也已是新版无 host) —— 此时若照写,
    #   就会用一个**空表覆盖掉原有 host.csv**, 静默丢掉全部宿主信息 (selftest 实测)。
    #   宁可"不落地 + 告警", 也不做不可逆的销毁。
    host_written = bool(out_host) and host_src_kind != "none"
    if host_written:
        _write_table(out_host, HOST_FIELDS, out_h)
    elif out_host:
        log.warning(f"  [channels] ⚠ 无任何 host 来源 (治理表无 host 列, 且未传入现有 "
                    f"{os.path.basename(out_host)}) → **跳过**宿主通道写出, "
                    f"保留既有内容不动")

    rep = {
        "time": datetime.now().isoformat(timespec="seconds"),
        # 布局版本: 见 CHANNELS_SCHEMA 注释 —— 旧版 channels.json 会被
        # needs_rebuild 判为"需重建", 从而把新增的宿主通道补齐。
        "channels_schema": CHANNELS_SCHEMA,
        "meta_source": meta_source, "geo_source": geo_source,
        "skeleton_dates": dates_csv, "skeleton_meta": meta_csv,
        "skeleton_host": host_csv,
        "out_dates": out_dates, "out_meta": out_meta,
        # 未落地时记 None → needs_rebuild 不会永远要求一个刻意不存在的产物
        "out_host": out_host if host_written else None,
        "fingerprint": source_fingerprint(meta_source, geo_source),
        "stats": stats,
    }
    rep_path = os.path.join(os.path.dirname(out_dates), "metadata", "channels.json")
    try:
        os.makedirs(os.path.dirname(rep_path), exist_ok=True)
        with open(rep_path, "w", encoding="utf-8") as f:
            json.dump(rep, f, ensure_ascii=False, indent=2)
    except OSError as e:
        log.warning(f"  [channels] channels.json 写出失败: {e}")
        rep_path = None

    hk = {"governance": "治理来源表", "host_channel": "现有 host.csv (自保持)",
          "legacy_geo_channel": "旧版地理表 (迁移)", "none": "**无** (跳过宿主通道)"}
    log.info(f"  [channels] 时间通道 {len(out_d)} 行 / 地理通道 {len(out_m)} 行 / "
             f"宿主通道 {len(out_h)} 行 (屏弃无宿主 {stats['n_host_dropped']}); "
             f"host 取值来源: {hk.get(host_src_kind, host_src_kind)}; "
             f"未匹配保留原值 {stats['n_unmatched']}")
    return {"success": True, "dates_csv": out_dates, "meta_csv": out_meta,
            "host_csv": out_host if host_written else None,
            "report_path": rep_path, "stats": stats,
            "host_written": host_written}


def needs_rebuild(report_path: str, meta_source: str,
                  geo_source: Optional[str] = None) -> bool:
    """来源指纹一致 + 布局版本一致 + 产物在 → 不必重建 (checkpoint 语义)。

    参数与 `rebuild_channels()` 同源, 避免两处各写一份清单而对不上。
    """
    if not os.path.exists(report_path):
        return True
    try:
        with open(report_path, encoding="utf-8") as f:
            rep = json.load(f)
    except (OSError, ValueError):
        return True
    # 布局版本不一致 → 必须重建 (否则新增的宿主通道永远不会出现, 见 CHANNELS_SCHEMA)
    if rep.get("channels_schema") != CHANNELS_SCHEMA:
        return True
    if rep.get("fingerprint") != source_fingerprint(meta_source, geo_source):
        return True
    for p in (rep.get("out_dates"), rep.get("out_meta"), rep.get("out_host")):
        if p and not os.path.exists(p):
            return True
    return False


# ══════════════════════════════════════════════════════════════════
# CLI
# ══════════════════════════════════════════════════════════════════
def _main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description="把治理后的 metadata 回填到 时间/地理/宿主 三张通道表")
    ap.add_argument("--dates", default=None, help="骨架: 现有 dates.csv")
    ap.add_argument("--meta", default=None, help="骨架: 现有 sample_metadata.csv")
    ap.add_argument("--host", default=None, help="骨架: 现有 host.csv (可选)")
    ap.add_argument("--source", default=None, help="值来源 (通常 metadata_std.csv)")
    ap.add_argument("--geo-source", default=None, help="地理值专用来源 (省略=同 --source)")
    ap.add_argument("--out-dates", default=None, help="输出 (省略=就地覆盖)")
    ap.add_argument("--out-meta", default=None, help="输出 (省略=就地覆盖)")
    ap.add_argument("--out-host", default=None, help="输出 (省略=就地覆盖)")
    ap.add_argument("--no-backup", action="store_true", help="不写 .bak_pre_governance")
    ap.add_argument("--selftest", action="store_true", help="跑内置自测")
    args = ap.parse_args(argv)

    logging.basicConfig(level=logging.INFO, format="%(message)s")
    if args.selftest:
        return _selftest()
    if not args.dates or not args.source:
        ap.error("非 --selftest 时必须给 --dates 与 --source")
    r = rebuild_channels(args.dates, args.meta, host_csv=args.host,
                         meta_source=args.source,
                         geo_source=args.geo_source, out_dates=args.out_dates,
                         out_meta=args.out_meta, out_host=args.out_host,
                         backup=not args.no_backup)
    print(json.dumps({k: v for k, v in r.items() if k != "stats"},
                     ensure_ascii=False, indent=2))
    return 0 if r.get("success") else 1


def _selftest() -> int:
    """内置自测: 临时目录造「脏骨架 + 治理后来源」, 断言回填规则逐条成立。"""
    import tempfile
    tmp = tempfile.mkdtemp(prefix="channels_selftest_")
    ok, fail = 0, 0

    def check(cond, msg):
        nonlocal ok, fail
        if cond:
            ok += 1
            print(f"  [ok  {ok}] {msg}")
        else:
            fail += 1
            print(f"  [FAIL] {msg}")

    d = os.path.join(tmp, "data")
    os.makedirs(d, exist_ok=True)
    # 骨架: 3 个样品, 日期含占位符与欧式格式, 无 location 列
    with open(os.path.join(d, "dates.csv"), "w", encoding="utf-8", newline="") as f:
        f.write("name,date\nA1,2020-01-01\nA2,15/07/2016\nA3,YYYY-MM-DD\n")
    # 地理骨架 (迁移期旧版式样: 仍带 host/tissue 列) → 地理通道只留地理,
    # host/tissue 由三级优先里的第 ③ 级 (legacy_geo_channel) 接走
    with open(os.path.join(d, "sample_metadata.csv"), "w", encoding="utf-8", newline="") as f:
        f.write("name,date,location,host,tissue\n"
                "A1,2020-01-01,Ningxia,Lycium barbarum,leaf\n"
                "A2,15/07/2016,China:Unknown:Yinchuan,Lycium ruthenicum,fruit\n"
                "A3,YYYY-MM-DD,Country:Region,Unknown,stem\n"
                "A4,2021-03-01,China:Beijing,uncultured bacterium,root\n")
    # 治理产物: A2 归一化, A3 占位符清空, location 删 Unknown 层; 只有 host 列缺失
    src = os.path.join(tmp, "metadata_std.csv")
    with open(src, "w", encoding="utf-8", newline="") as f:
        f.write("name,collection_date,geo_loc_name\n"
                "A1,2020-01-01,Ningxia\n"
                "A2,2016-07-15,China:Yinchuan\n"
                "A3,,\n")          # 占位符 → 清空

    r = rebuild_channels(os.path.join(d, "dates.csv"),
                         os.path.join(d, "sample_metadata.csv"), meta_source=src)
    check(r.get("success"), "回填成功返回 success=True")
    st = r["stats"]
    check(st["src_date_col"] == "collection_date", "识别治理表日期列 collection_date")
    check(st["src_loc_col"] == "geo_loc_name", "识别治理表地名列 geo_loc_name")
    check(st["n_matched"] == 3, "3 行全部按名匹配")

    with open(os.path.join(d, "dates.csv"), encoding="utf-8") as f:
        dt = {r_["name"]: r_ for r_ in csv.DictReader(f)}
    check(dt["A2"]["date"] == "2016-07-15", "欧式日期 15/07/2016 → 2016-07-15")
    check(dt["A3"]["date"] == "", "占位符日期被清空 (正是要修的)")
    check(list(dt) == ["A1", "A2", "A3"], "行顺序不变 (不增不删行)")

    with open(os.path.join(d, "sample_metadata.csv"), encoding="utf-8") as f:
        _mh = csv.DictReader(f)
        _mcols = _mh.fieldnames
        mt = {r_["name"]: r_ for r_ in _mh}
    check(mt["A2"]["location"] == "China:Yinchuan", "Unknown 地理层级已删 (回填生效)")
    check(mt["A3"]["location"] == "", "占位符地理被清空")
    # ── 地理通道**只留地理**: host/tissue 列已迁出 (老师「host 也独立出来」) ──
    check(_mcols == ["name", "date", "location"],
          f"地理通道列 = name/date/location (实际 {_mcols})")

    # ── 宿主通道: 独立文件 + 只收有宿主的样品 ──
    _hp = os.path.join(d, "host.csv")
    check(os.path.exists(_hp), "独立写出 host.csv (宿主通道)")
    with open(_hp, encoding="utf-8") as f:
        _hh = csv.DictReader(f)
        _hcols = _hh.fieldnames
        ht = {r_["name"]: r_ for r_ in _hh}
    check(_hcols == ["name", "host", "tissue"], f"宿主通道列 (实际 {_hcols})")
    check(sorted(ht) == ["A1", "A2"],
          f"只收有宿主的样品 (A3=Unknown / A4=uncultured bacterium 均被屏弃; 实际 {sorted(ht)})")
    check(ht["A1"]["host"] == "Lycium barbarum", "宿主值保留 (来自旧版地理表)")
    check(ht["A1"]["tissue"] == "leaf", "tissue 随宿主一并迁入")
    check(st["host_src_kind"] == "legacy_geo_channel",
          f"host 取值来源判定 = legacy_geo_channel (实际 {st['host_src_kind']})")
    check(st["n_host_dropped"] == 2, f"屏弃无宿主 2 条 (实际 {st['n_host_dropped']})")

    check(os.path.exists(os.path.join(d, "dates.csv.bak_pre_governance")),
          "原件已备份 .bak_pre_governance")
    _rep = os.path.join(d, "metadata", "channels.json")
    check(os.path.exists(_rep), "写出 channels.json 留痕")
    check(not needs_rebuild(_rep, src),
          "指纹一致时 needs_rebuild=False (checkpoint 语义)")
    src2 = os.path.join(tmp, "metadata_std_v2.csv")
    with open(src2, "w", encoding="utf-8", newline="") as f:
        f.write("name,collection_date,geo_loc_name\nA1,2021-01-01,Ningxia\n")
    check(needs_rebuild(_rep, src2), "换了来源文件 → needs_rebuild=True")

    # ── 布局版本 (CHANNELS_SCHEMA): 旧版 channels.json 必须被判"需重建" ──
    # 否则 2→3 张表的升级永远不会发生 (旧 report 无 out_host, 来源又没变 → 误判 skip)
    _snap = json.load(open(_rep, encoding="utf-8"))
    _no_schema = {k: v for k, v in _snap.items() if k != "channels_schema"}
    json.dump(_no_schema, open(_rep, "w", encoding="utf-8"), ensure_ascii=False)
    check(needs_rebuild(_rep, src), "旧版 report 无 channels_schema → needs_rebuild=True (强制升级)")
    json.dump(_snap, open(_rep, "w", encoding="utf-8"), ensure_ascii=False)
    check(not needs_rebuild(_rep, src), "版本补齐后恢复 skip (不空转)")

    # ── 宿主通道幂等: host.csv 自身作为下一轮骨架, host_src_kind 应降为自保持 ──
    r_h = rebuild_channels(os.path.join(d, "dates.csv"),
                           os.path.join(d, "sample_metadata.csv"),
                           host_csv=_hp, meta_source=src, backup=False)
    with open(_hp, encoding="utf-8") as f:
        ht2 = {r_["name"]: r_ for r_ in csv.DictReader(f)}
    check(sorted(ht2) == ["A1", "A2"] and ht2["A1"]["host"] == "Lycium barbarum",
          "宿主通道二次回填幂等 (行集与值不变)")

    # ── 治理来源自带 host 列 → 优先级 ① 生效 (权威覆盖) ──
    src3 = os.path.join(tmp, "metadata_std_v3.csv")
    with open(src3, "w", encoding="utf-8", newline="") as f:
        f.write("name,collection_date,geo_loc_name,host\n"
                "A1,2020-01-01,Ningxia,Lycium chinense\n")
    r3 = rebuild_channels(os.path.join(d, "dates.csv"),
                          os.path.join(d, "sample_metadata.csv"),
                          host_csv=_hp, meta_source=src3, backup=False)
    check(r3["stats"]["host_src_kind"] == "governance",
          f"治理来源带 host 列 → 优先级 ① 生效 (实际 {r3['stats']['host_src_kind']})")
    with open(_hp, encoding="utf-8") as f:
        ht3 = {r_["name"]: r_ for r_ in csv.DictReader(f)}
    check(ht3["A1"]["host"] == "Lycium chinense", "治理来源的 host 覆盖了骨架原值")
    # 取值链兜底: src3 只含 A1 → A2 在治理来源里"查无此人", 不能因此丢掉宿主
    check(sorted(ht3) == ["A1", "A2"],
          f"治理来源只含部分样品时, 未覆盖样品仍保留宿主 (实际 {sorted(ht3)})")
    check(r3["stats"]["host_src_used"].get("host_channel") == 1,
          f"兜底确实命中 host.csv 一级 (实际 {r3['stats']['host_src_used']})")

    # 未匹配: 治理表里没有的名字保留原值
    with open(src, "a", encoding="utf-8", newline="") as f:
        f.write("")
    with open(src, encoding="utf-8", newline="") as f:
        body = f.read()
    with open(src, "w", encoding="utf-8", newline="") as f:
        f.write(body.replace("A3,,\n", ""))          # 去掉 A3
    r2 = rebuild_channels(os.path.join(d, "dates.csv"),
                          os.path.join(d, "sample_metadata.csv"), meta_source=src,
                          backup=False)
    check(r2["stats"]["n_unmatched"] == 1, "来源缺 A3 → unmatched=1")
    with open(os.path.join(d, "dates.csv"), encoding="utf-8") as f:
        dt2 = {r_["name"]: r_ for r_ in csv.DictReader(f)}
    check(dt2["A3"]["date"] == "", "未匹配行保留上一轮 (空) 值, 不报错")

    # ── fail-safe: 三处都拿不到 host 时, **不得用空表覆盖**既有 host.csv ──
    # (上面这次 r2 调用刻意没传 host_csv, 而地理表已是新版无 host 列 →
    #  host_src_kind == "none"。若照写就会静默清空宿主通道。)
    with open(_hp, encoding="utf-8") as f:
        _after_r2 = {r_["name"]: r_ for r_ in csv.DictReader(f)}
    check(r2["stats"]["host_src_kind"] == "none",
          f"无 host 来源时判定为 none (实际 {r2['stats']['host_src_kind']})")
    check(not r2["host_written"] and r2.get("host_csv") is None,
          "无 host 来源 → 不落地宿主通道 (host_written=False)")
    check(sorted(_after_r2) == ["A1", "A2"],
          f"既有 host.csv 未被空表覆盖 (实际 {sorted(_after_r2)})")
    with open(_rep, encoding="utf-8") as f:
        check(json.load(f).get("out_host") is None,
              "未落地时 channels.json 记 out_host=None (needs_rebuild 不空转)")

    print()
    if fail:
        print(f"SELFTEST FAILED — {ok} passed, {fail} failed")
        return 1
    print(f"SELFTEST PASSED — {ok} 项断言全部通过")
    return 0


if __name__ == "__main__":
    sys.exit(_main())
