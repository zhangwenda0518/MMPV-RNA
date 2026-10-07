#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
stat_prescreen.py — 基于 NCBI STAT 的 SRA 病毒预筛选器

用途
----
给定一批 SRA Run accession（SRR / ERR / DRR / CRR），向 NCBI Trace 与 CNCB GSA
取回该 run 的**预计算物种分类组成**，同时回答两个问题：

  (1) 这条 run 的物种记录对吗？  —— 用实测主要物种对拍元信息 ScientificName
  (2) 这条 run 里有病毒吗？      —— 病毒占比 / 归一化占比 / 主要病毒

输出（按用途分层）
------------------
主表
  sample_decision.tsv       每 run 一行，规范 schema：状态 / 物种对拍 / 病毒预鉴定 / 决策
派生清单（可直接喂下游，每行一个 accession）
  prescreen_run.list        建议下载的 run（decision = DOWNLOAD / REVIEW）
  prescreen_local_check.list 无预计算分类，只能本地比对的 run（decision = LOCAL_CHECK）
  prescreen_skip.list       明确无病毒信号的 run（decision = SKIP）
明细与兼容
  stat_run_summary.tsv      兼容旧版的宽表（含全部列，列名未变）
  stat_viral_taxa.tsv       病毒类 taxon 长表（供下游富集 / 系统发育取样）
  prescreen_pass.tsv        通过阈值的清单（按病毒占比降序）
  unresolved.tsv / .list    无分类数据的 run 及原因
元信息
  prescreen_meta.json       参数快照 / 输入指纹 / 判定依据（可复现性）
  STAT_PRESCREEN_REPORT.md  人读摘要

决策语义（decision 列）
-----------------------
  DOWNLOAD     检出病毒信号且物种对拍一致  -> 值得下载
  REVIEW       检出病毒但主要物种与元信息不符 -> 仍下载，但人工复核（疑似元信息记错）
  SKIP         有分类数据但无病毒信号      -> 可不下
  LOCAL_CHECK  拿不到预计算分类           -> 只能下载 reads 本地比对

核心价值：**不下载任何测序数据**即可把上万条候选 run 先筛一遍，
只把真正含病毒的样本送进 virome_discovery_pipeline / virome_analysis_pipeline；
同时用实测物种反向校验元数据，避免"元信息记错物种"一路带到下游。

数据来源
--------
    GET https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/run_taxonomy
          ?acc={RUN}&cluster_name=public

    -> [{"tax_totals": {...}, "tax_table": [{"tax_id","org","self_count","total_count"}, ...]}]

    ⚠ 必须带 cluster_name=public，否则返回 HTTP 200 + 空 body（这个坑很隐蔽）。

口径说明（来自 NCBI 官方答复 ngs-tools#43）
-------------------------------------------
    analysed     : 分母。被赋予 taxid 的 spot 总数。
    total_count  : 该节点 **含全部后代** 的累计计数（所以 "cellular organisms" 极大）。
    identified   : 至少命中 1 个 k-mer 的 spot 数。
    占比 pct     = total_count / analysed * 100

覆盖范围（本机实测，2026-10-01）
--------------------------------
    SRR (NCBI)      ✅ 有数据
    DRR (DDBJ)      ✅ 有数据（DRR839754 -> 1280 个 taxon）
    ERR (ENA)       ✅ 有数据
    CRR (CNCB/GSA)  ❌ 不在 SRA 内，接口返回空 -> 状态 NO_STAT
    新增/未分析完的 run 同样返回空 -> NO_STAT

    => CRR 样本无法用本方法预筛，需走本地比对路线（见 --help 末尾说明）。

用法
----
    # 1) 从文本列表
    python stat_prescreen.py --runs sra.list -o ./stat_out

    # 2) 直接从 pipeline 的元数据 CSV 取 Run 列
    python stat_prescreen.py --csv Global_Unified_Metadata_Core14.csv \
                             --run-col Run -o ./stat_out

    # 3) 只看病毒占比 > 0.5% 的样本，并把通过的清单单独导出
    python stat_prescreen.py --runs sra.list -o ./stat_out --virus-pct 0.5

    # 4) 精确定位病毒（按 NCBI Taxonomy 子树判定，而非名字猜）
    python stat_prescreen.py --runs sra.list -o ./stat_out --taxdump ./taxdump

    # 5) 断点续跑（默认开启，已下载的 JSON 直接复用）
    python stat_prescreen.py --runs sra.list -o ./stat_out
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import random
import re
import sys
import threading
import time
import urllib.error
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime
from typing import Any, Dict, Iterable, List, Optional, Set, Tuple

if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass

# ---- 复用 pipeline 的 UI，保持输出风格一致；脱离 pipeline 时降级为内置实现 ----
try:
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "utils"))
    from pipeline_utils import UI  # type: ignore
except Exception:  # pragma: no cover
    class UI:  # type: ignore
        C = {"cyan": "\033[96m", "green": "\033[92m", "yellow": "\033[93m",
             "red": "\033[91m", "gray": "\033[90m", "purple": "\033[95m",
             "bold": "\033[1m", "reset": "\033[0m"}

        @staticmethod
        def banner(t): print(f"\n{UI.C['purple']}{UI.C['bold']}  {t}{UI.C['reset']}\n")
        @classmethod
        def stage(cls, n, s="start"): print(f"{cls.C['cyan']}[>>>] {n}{cls.C['reset']}")
        @classmethod
        def ok(cls, m): print(f"  {cls.C['green']}[OK]{cls.C['reset']} {m}")
        @classmethod
        def warn(cls, m): print(f"  {cls.C['yellow']}[WARN]{cls.C['reset']} {m}")
        @classmethod
        def err(cls, m): print(f"  {cls.C['red']}[ERR]{cls.C['reset']} {m}")
        @classmethod
        def info(cls, m): print(f"  {cls.C['gray']}> {m}{cls.C['reset']}")


# ============================================================================
# 常量
# ============================================================================
# 版本号：输出 schema 或判定规则变更时 +1。写入 prescreen_meta.json，
# 下游可据此判断自己读到的表是哪一代格式。
SCRIPT_VERSION = "2.0"
SCHEMA_VERSION = "2"

STAT_URL = "https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/run_taxonomy?acc={acc}&cluster_name=public"

TAXID_VIRUSES = 10239          # NCBI: Viruses
TAXID_ROOT = 1
TAXID_CELLULAR = 131567        # cellular organisms

# 名字兜底启发式：当没有 --taxdump 时，用它从 tax_table 里挑出病毒行。
# 这是「够用但不严谨」的判据，仅供预览；严格判定请用 --taxdump。
VIRAL_NAME_RE = re.compile(
    r"(viridae|virales|virinae|virus|viruses|viroid|virophage|"
    r"\bphage\b|bacteriophage|retroelement|Caulimoviridae)",
    re.IGNORECASE,
)

USER_AGENT = "stat-prescreen/1.0 (research; contact: local pipeline)"

# ---- GSA (CNCB / NGDC) ----
# GSA 对 CRR 也做了一套独立的分类分析，但格式与 NCBI STAT 完全不同：
#   * 结构：嵌套树（name / hits / percent / children），不是扁平表
#   * 没有 tax_id，只有名字  -> 只能按名字判定病毒
#   * percent 的**分母是总 reads**（identified + unidentified），
#     而 STAT 的分母是 analysed（已鉴定 reads）。
#     两者不可直接比较，因此本脚本额外给出归一化列 virus_pct_ident。
GSA_BASE = "https://ngdc.cncb.ac.cn"
GSA_SEARCH = GSA_BASE + "/gsa/search?searchTerm={acc}"
GSA_TREE = GSA_BASE + "/gsa/taxonomyAnalysis/treePath/{cra}/{acc}/{acc}.json"
# 搜索页模板里恒定出现的噪声 CRA，必须排除，否则会把它当成候选
GSA_CRA_NOISE = {"CRA000112"}
# GSA 缓存格式版本。改动缓存结构时 +1，旧缓存会被自动作废并重抓。
GSA_CACHE_V = 2

# GSA 树里的病毒节点名（GSA 有顶层 "Viruses" 节点，名字判定比 STAT 更可靠）
GSA_VIRAL_RE = re.compile(
    r"(^Viruses$|viridae|virales|virinae|virus|viruses|viroid|virophage|phage)",
    re.IGNORECASE,
)


# ============================================================================
# Taxonomy（可选）：精确判定 taxid 是否落在 Viruses(10239) 子树内
# ============================================================================
class Taxonomy:
    """从 NCBI taxdump 的 nodes.dmp 构建 parent 映射，支持祖先查询。"""

    def __init__(self, taxdump_dir: str) -> None:
        nodes = os.path.join(taxdump_dir, "nodes.dmp")
        if not os.path.exists(nodes):
            raise FileNotFoundError(f"未找到 {nodes}")
        self.parent: Dict[int, int] = {}
        with open(nodes, "r", encoding="utf-8", errors="replace") as f:
            for line in f:
                # 格式: tax_id \t|\t parent_tax_id \t|\t rank ...
                parts = line.split("\t|\t")
                if len(parts) < 2:
                    continue
                try:
                    self.parent[int(parts[0])] = int(parts[1])
                except ValueError:
                    continue
        self._cache: Dict[int, bool] = {}

    def is_under(self, taxid: int, ancestor: int, _guard: int = 64) -> bool:
        """判断 taxid 是否位于 ancestor 的子树内（含自身）。"""
        if taxid == ancestor:
            return True
        if taxid in self._cache:
            return self._cache[taxid]
        cur, hops, seen = taxid, 0, set()
        while hops < _guard:
            p = self.parent.get(cur)
            if p is None or p == cur or p in seen:
                break
            if p == ancestor:
                self._cache[taxid] = True
                return True
            seen.add(cur)
            cur, hops = p, hops + 1
            if cur == TAXID_ROOT:
                break
        self._cache[taxid] = False
        return False


# ============================================================================
# 抓取
# ============================================================================
class Fetcher:
    """带重试/退避/限速的 STAT 抓取器，并落盘缓存以支持断点续跑。"""

    def __init__(self, raw_dir: str, retries: int = 4, timeout: int = 90,
                 min_interval: float = 0.0) -> None:
        self.raw_dir = raw_dir
        self.retries = retries
        self.timeout = timeout
        self.min_interval = min_interval
        self._lock = threading.Lock()
        self._last = 0.0
        # 读缓存失败 (损坏/形态不认识) 的计数 —— 旧版静默 pass, 无法发现缓存全失效
        self.n_cache_corrupt = 0
        os.makedirs(raw_dir, exist_ok=True)

    def _throttle(self) -> None:
        if self.min_interval <= 0:
            return
        with self._lock:
            gap = time.time() - self._last
            if gap < self.min_interval:
                time.sleep(self.min_interval - gap)
            self._last = time.time()

    def cache_path(self, acc: str) -> str:
        return os.path.join(self.raw_dir, f"{acc}.json")

    def fetch(self, acc: str, use_cache: bool = True) -> Tuple[str, Optional[Dict[str, Any]], str, float]:
        """返回 (acc, payload|None, status, seconds)。

        status ∈ {OK, NO_STAT, FAILED}
        """
        cp = self.cache_path(acc)
        if use_cache and os.path.exists(cp):
            try:
                with open(cp, "r", encoding="utf-8") as f:
                    txt = f.read()
                if txt.strip() == "":          # 缓存了「无数据」这一事实
                    return acc, None, "NO_STAT", 0.0
                # 兼容两种历史形态: 早期写端 json.dump(payload) 落 dict, 后统一为
                # json.dump([payload]) 落 list。旧读端固定取 [0] —— 对 dict 恒
                # KeyError, 被 except 吞掉, 于是缓存永不命中 (磁盘上几百个 dict
                # 缓存全废, 每次重抓)。两种形态都要认。
                obj = json.loads(txt)
                return acc, (obj[0] if isinstance(obj, list) else obj), "OK", 0.0
            except Exception:
                self.n_cache_corrupt += 1
                UI.warn(f"缓存不可用 (已重抓): {os.path.basename(cp)}")

        url = STAT_URL.format(acc=urllib.parse.quote(acc))
        last_err = ""
        t0 = time.time()
        for attempt in range(self.retries):
            self._throttle()
            try:
                req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
                with urllib.request.urlopen(req, timeout=self.timeout) as resp:
                    raw = resp.read()
                if not raw:
                    # 空 body = 该 accession 无 STAT 结果（CRR / 新 run / 不存在）
                    with open(cp, "w", encoding="utf-8") as f:
                        f.write("")
                    return acc, None, "NO_STAT", time.time() - t0
                payload = json.loads(raw)[0]
                with open(cp, "w", encoding="utf-8") as f:
                    json.dump([payload], f, ensure_ascii=False)   # list 形态, 与读端/上游一致
                return acc, payload, "OK", time.time() - t0
            except urllib.error.HTTPError as e:
                last_err = f"HTTP {e.code}"
                if e.code in (429, 500, 502, 503, 504):
                    time.sleep(min(2 ** attempt + random.random(), 30))
                    continue
                break
            except Exception as e:                      # 超时 / 连接重置
                last_err = f"{type(e).__name__}: {e}"
                time.sleep(min(2 ** attempt + random.random(), 30))
        return acc, None, f"FAILED ({last_err})", time.time() - t0


# ============================================================================
# GSA 抓取（CRR）
# ============================================================================
class GsaFetcher:
    """抓取 NGDC/GSA 的 CRR 分类树。

    与 NCBI 不同，GSA 的 URL 需要 (CRA, CRR) 两个 accession，
    而输入通常只有 CRR，因此要先解析 CRA：
      1. --cra-map 提供的映射（最快，推荐用于大批量）
      2. 抓 gsa/search?searchTerm={CRR}，正则出候选 CRA（排除模板噪声），
         再逐个用 treePath 验证 —— 验证过才算数，不靠猜。
    """

    def __init__(self, raw_dir: str, retries: int = 4, timeout: int = 60,
                 min_interval: float = 0.4, cra_map: Optional[Dict[str, str]] = None) -> None:
        self.raw_dir = raw_dir
        self.retries = retries
        self.timeout = timeout
        self.min_interval = min_interval
        self.cra_map = cra_map or {}
        self._lock = threading.Lock()
        self._last = 0.0
        self._cra_cache: Dict[str, List[str]] = {}
        os.makedirs(raw_dir, exist_ok=True)

    def _throttle(self) -> None:
        with self._lock:
            gap = time.time() - self._last
            if gap < self.min_interval:
                time.sleep(self.min_interval - gap)
            self._last = time.time()

    def _get(self, url: str) -> str:
        last_err: Optional[Exception] = None
        for attempt in range(self.retries):
            self._throttle()
            try:
                req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
                with urllib.request.urlopen(req, timeout=self.timeout) as r:
                    return r.read().decode("utf-8", errors="replace")
            except urllib.error.HTTPError as e:
                last_err = e
                if e.code in (429, 500, 502, 503, 504):
                    # NGDC 对高频访问会 429，退避要够长
                    time.sleep(min(3 * (attempt + 1) + random.random(), 45))
                    continue
                raise
            except Exception as e:
                last_err = e
                time.sleep(min(2 ** attempt + random.random(), 30))
        raise last_err if last_err else RuntimeError("unknown")

    def cra_candidates(self, acc: str) -> List[str]:
        if acc in self._cra_cache:
            return self._cra_cache[acc]
        if acc in self.cra_map:
            cands = [self.cra_map[acc]]
        else:
            html = self._get(GSA_SEARCH.format(acc=acc))
            cands = [c for c in sorted(set(re.findall(r"CRA\d{6}", html)))
                     if c not in GSA_CRA_NOISE]
        self._cra_cache[acc] = cands
        return cands

    def cache_path(self, acc: str) -> str:
        return os.path.join(self.raw_dir, f"{acc}.json")

    def _write_cache(self, cp: str, status: str, cra: str = "",
                     detail: str = "", tree: Optional[Dict[str, Any]] = None) -> None:
        """统一缓存包装。

        带 `_v` 版本号是有意的：旧格式（空文件 / 裸 tree）一律视为不可信并重抓，
        否则旧版把「临时故障」误记成「无数据」的结果会被一直沿用。
        """
        payload: Dict[str, Any] = {"_v": GSA_CACHE_V, "_status": status,
                                   "_cra": cra, "_detail": detail}
        if tree is not None:
            payload["_tree"] = tree
        with open(cp, "w", encoding="utf-8") as f:
            json.dump(payload, f, ensure_ascii=False)

    def fetch(self, acc: str, use_cache: bool = True,
              recheck_nodata: bool = False) -> Tuple[str, Optional[Dict[str, Any]], str, float, str, str]:
        """返回 (acc, tree|None, status, seconds, cra, detail)。

        status 取值：
          OK           拿到分类树
          NO_ANALYSIS  CRA 解析成功，但该 run 没有分类分析结果（treePath 404）—— 确认态，可缓存
          NOT_FOUND    搜索页找不到任何 CRA 候选，accession 不在 GSA —— 确认态，可缓存
          FAILED       429 / 网络 / 超时等**临时**故障 —— 绝不缓存，下次重试

        只有「确认态」才写缓存。临时故障一旦被写成 NO_DATA，
        就会永久污染结果（这正是本脚本早期版本的一个 bug）。
        """
        cp = self.cache_path(acc)

        if use_cache and os.path.exists(cp):
            try:
                obj = json.loads(open(cp, "r", encoding="utf-8").read() or "null")
            except Exception:
                obj = None
            # 只认当前版本的缓存包装；旧格式（空文件 / 裸 tree）一律重抓
            if isinstance(obj, dict) and obj.get("_v") == GSA_CACHE_V:
                st = obj.get("_status", "")
                if st == "OK":
                    tree = obj.get("_tree")
                    if isinstance(tree, dict):
                        return acc, tree, "OK", 0.0, obj.get("_cra", ""), ""
                elif not recheck_nodata:
                    # --recheck-nodata：忽略「无数据」标记，重新问一次 GSA
                    # （GSA 会陆续补做分析，昨天的无数据今天可能就有了）
                    return (acc, None, st, 0.0,
                            obj.get("_cra", ""), obj.get("_detail", ""))

        t0 = time.time()

        # ---- 1) 解析 CRA ----
        try:
            cands = self.cra_candidates(acc)
        except Exception as e:
            # 解析失败多半是限流/网络 —— 临时故障，不缓存
            return acc, None, f"FAILED (CRA解析 {type(e).__name__})", time.time() - t0, "", ""

        if not cands:
            d = "搜索页无 CRA 候选（accession 不在 GSA）"
            self._write_cache(cp, "NOT_FOUND", detail=d)
            return acc, None, "NOT_FOUND", time.time() - t0, "", d

        # ---- 2) 逐个候选 CRA 取树 ----
        confirmed_missing = 0     # 明确 404 的候选数
        transient: List[str] = []  # 临时故障的候选
        for cra in cands:
            try:
                raw = self._get(GSA_TREE.format(cra=cra, acc=acc))
                if not raw.strip():
                    confirmed_missing += 1
                    continue
                tree = json.loads(raw)
                if not isinstance(tree, dict):
                    confirmed_missing += 1
                    continue
                self._write_cache(cp, "OK", cra=cra, tree=tree)
                return acc, tree, "OK", time.time() - t0, cra, ""
            except urllib.error.HTTPError as e:
                if e.code == 404:
                    confirmed_missing += 1
                    continue
                transient.append(f"{cra}:HTTP{e.code}")
            except json.JSONDecodeError:
                confirmed_missing += 1
            except Exception as e:
                transient.append(f"{cra}:{type(e).__name__}")

        # 只要有一个候选是**临时**故障，就不能断定「没有数据」
        if transient:
            return (acc, None, f"FAILED ({';'.join(transient)})", time.time() - t0, "",
                    "候选 CRA 全部请求失败（临时故障，未缓存）")

        cra_str = "/".join(cands)
        detail = (f"候选 CRA {cra_str} 均返回 404，"
                  f"GSA 未对该 run 做分类分析（可日后用 --recheck-nodata 复查）")
        self._write_cache(cp, "NO_ANALYSIS", cra=cra_str, detail=detail)
        return acc, None, "NO_ANALYSIS", time.time() - t0, cra_str, detail


# ============================================================================
# 解析
# ============================================================================
def parse_payload(payload: Dict[str, Any], tax: Optional[Taxonomy],
                  viral_taxids: Set[int]) -> Dict[str, Any]:
    """把 STAT payload 提炼成一行汇总 + 病毒明细。"""
    tot = payload.get("tax_totals", {}) or {}
    table = payload.get("tax_table", []) or []

    analysed = int(tot.get("analysed") or tot.get("total") or 0)
    identified = int(tot.get("identified") or 0)
    denom = analysed if analysed > 0 else 0

    def pct(c: int) -> float:
        return round(100.0 * c / denom, 6) if denom else 0.0

    row: Dict[str, Any] = {
        "source": "NCBI_STAT",
        "analysed": analysed,
        "identified": identified,
        "identified_pct": pct(identified),
        "n_taxa": len(table),
        "virus_count": 0,
        "virus_pct": 0.0,
        "virus_pct_ident": 0.0,
        "top_organism": "",
        "top_organism_pct": 0.0,
    }

    # --- 病毒判定 ---
    viral_rows: List[Dict[str, Any]] = []
    virus_node = next((r for r in table if r.get("tax_id") == TAXID_VIRUSES), None)
    if virus_node:
        row["virus_count"] = int(virus_node.get("total_count", 0))
        row["virus_pct"] = pct(row["virus_count"])

    for r in table:
        tid = r.get("tax_id")
        if tid is None or tid in (TAXID_ROOT, TAXID_CELLULAR, TAXID_VIRUSES):
            continue
        if tax is not None:
            hit = tax.is_under(int(tid), TAXID_VIRUSES)
        else:
            hit = bool(VIRAL_NAME_RE.search(r.get("org", "")))
        if hit:
            viral_rows.append({
                "tax_id": int(tid),
                "org": r.get("org", ""),
                "path": "",
                "self_count": int(r.get("self_count", 0)),
                "total_count": int(r.get("total_count", 0)),
                "pct": pct(int(r.get("total_count", 0))),
            })

    # 病毒明细按 total_count 降序（父节点在前，便于看层级）
    viral_rows.sort(key=lambda x: -x["total_count"])
    row["_viral_rows"] = viral_rows

    # 若没有显式的 Viruses(10239) 行，用检出的病毒节点里 total_count 最大者兜底
    if not virus_node and viral_rows:
        row["virus_count"] = viral_rows[0]["total_count"]
        row["virus_pct"] = viral_rows[0]["pct"]

    # 归一化口径：病毒命中 / 已鉴定命中。用于和 GSA 的结果跨来源比较
    # （GSA 的 percent 分母是总 reads，STAT 的是 analysed，直接比会低估 GSA）
    if identified > 0:
        row["virus_pct_ident"] = round(100.0 * row["virus_count"] / identified, 6)

    # 「最特异的病毒」：按 **self_count**（直接命中该节点的读数，不含后代继承）选。
    # 用 total_count 会永远选中 "Viruses"/"ssRNA viruses" 这类高层节点，没有信息量。
    if viral_rows:
        best_v = max(viral_rows, key=lambda x: (x["self_count"], x["total_count"]))
        row["top_virus"] = best_v["org"]
        row["top_virus_pct"] = best_v["pct"]
    elif row["virus_count"] > 0:
        # 有病毒信号但没有细分到具名节点 —— 别留空，否则下游会以为没检出
        row["top_virus"] = "Viruses (未细分)"
        row["top_virus_pct"] = row["virus_pct"]
    else:
        row["top_virus"] = ""
        row["top_virus_pct"] = 0.0

    # 非病毒的主要物种：同样排除病毒子树，并用 self_count 选最特异的那个
    non_viral = [
        r for r in table
        if r.get("tax_id") not in (TAXID_ROOT, TAXID_CELLULAR, TAXID_VIRUSES)
        and not (tax.is_under(int(r["tax_id"]), TAXID_VIRUSES) if tax else
                 VIRAL_NAME_RE.search(r.get("org", "")))
    ]
    if non_viral:
        best_o = max(non_viral, key=lambda x: (int(x.get("self_count", 0)),
                                               int(x.get("total_count", 0))))
        row["top_organism"] = best_o.get("org", "")
        row["top_organism_pct"] = pct(int(best_o.get("total_count", 0)))
    return row


def _gsa_flat(node: Dict[str, Any], path: Tuple[str, ...] = ()) -> List[Dict[str, Any]]:
    """把 GSA 的嵌套树拍平成列表，保留完整路径（GSA 没有 tax_id，路径即层级证据）。"""
    p = path + (node.get("name", ""),)
    out = [{"name": node.get("name", ""), "hits": node.get("hits"),
            "percent": node.get("percent"), "path": " > ".join(p)}]
    for c in node.get("children") or []:
        out.extend(_gsa_flat(c, p))
    return out


def _gsa_pct(v: Any) -> float:
    """'0.04%' -> 0.04；无效返回 0。"""
    if not isinstance(v, str):
        return 0.0
    try:
        return float(v.rstrip("%").strip())
    except ValueError:
        return 0.0


def _gsa_hits(v: Any) -> int:
    try:
        return int(str(v).strip())
    except (TypeError, ValueError):
        return 0


def parse_gsa_tree(tree: Dict[str, Any]) -> Dict[str, Any]:
    """把 GSA 分类树提炼成与 STAT 路径同构的一行汇总。

    注意口径差异：
      * GSA 的 percent 分母是**总 reads**（identified + unidentified）
      * STAT 的分母是 analysed（已鉴定 reads）
    所以这里额外算 virus_pct_ident = 病毒命中 / 已鉴定命中，用于跨来源比较。
    """
    rows = _gsa_flat(tree)
    root = tree
    ident_pct = _gsa_pct(root.get("identified"))
    unident_pct = _gsa_pct(root.get("unidentified"))

    # 用命中最多的节点反推总 reads： total = hits / (percent/100)
    total_reads = 0.0
    for r in rows:
        h, p = _gsa_hits(r["hits"]), _gsa_pct(r["percent"])
        if h > 0 and p > 0:
            total_reads = max(total_reads, h / (p / 100.0))
    ident_hits = int(round(total_reads * ident_pct / 100.0)) if total_reads else 0

    row: Dict[str, Any] = {
        "source": "GSA",
        # analysed 统一语义 = virus_pct 的分母。GSA 的分母是总 reads。
        "analysed": int(round(total_reads)) if total_reads else 0,
        "identified": ident_hits,
        "identified_pct": round(ident_pct, 4),
        "n_taxa": len(rows) - 1,
        "virus_count": 0,
        "virus_pct": 0.0,
        "virus_pct_ident": 0.0,
        "top_organism": "",
        "top_organism_pct": 0.0,
    }

    # --- 病毒节点 ---
    viral: List[Dict[str, Any]] = []
    virus_root = next((r for r in rows if (r["name"] or "").strip().lower() == "viruses"), None)
    for r in rows:
        nm = r["name"] or ""
        if nm.strip().lower() == "viruses":
            continue
        if not GSA_VIRAL_RE.search(nm):
            continue
        hits = _gsa_hits(r["hits"])
        viral.append({
            "tax_id": "",
            "org": nm,
            "path": r["path"],
            "self_count": "",
            "total_count": hits,
            "pct": round(_gsa_pct(r["percent"]), 6),
        })

    if virus_root is not None:
        row["virus_count"] = _gsa_hits(virus_root["hits"])
        row["virus_pct"] = round(_gsa_pct(virus_root["percent"]), 6)

    # 没有顶层 Viruses 节点时，用病毒节点里的最大值兜底
    if not row["virus_count"] and viral:
        row["virus_count"] = max(v["total_count"] for v in viral)
        row["virus_pct"] = max(v["pct"] for v in viral)

    if ident_hits > 0:
        row["virus_pct_ident"] = round(100.0 * row["virus_count"] / ident_hits, 6)

    # 最特异的病毒 = 「没有病毒后代」的那些节点里命中最多的
    viral_names = {v["org"] for v in viral}
    leaves = [v for v in viral
              if not any(o["path"].startswith(v["path"] + " > ") and o["org"] in viral_names
                         for o in viral)]
    pool = leaves or viral
    if pool:
        best_v = max(pool, key=lambda x: x["total_count"])
        row["top_virus"] = best_v["org"]
        row["top_virus_pct"] = best_v["pct"]
    elif row["virus_count"] > 0:
        # 有病毒信号但没有细分到具名节点 —— 别留空，否则下游会以为没检出
        row["top_virus"] = "Viruses (未细分)"
        row["top_virus_pct"] = row["virus_pct"]

    viral.sort(key=lambda x: -x["total_count"])
    row["_viral_rows"] = viral

    # 主要非病毒物种：最深层且命中最多的非病毒、非根节点
    non_viral = [r for r in rows
                 if (r["name"] or "").strip().lower() not in ("", "viruses")
                 and not GSA_VIRAL_RE.search(r["name"] or "")
                 and r["name"] not in ("cellular organisms",)]
    if non_viral:
        maxdepth = max(r["path"].count(" > ") for r in non_viral)
        deep = [r for r in non_viral if r["path"].count(" > ") == maxdepth]
        best_o = max(deep, key=lambda x: _gsa_hits(x["hits"]))
        row["top_organism"] = best_o["name"]
        row["top_organism_pct"] = round(_gsa_pct(best_o["percent"]), 6)
    return row


# ============================================================================
# 物种鉴定：元信息记录 vs 预筛实测 对拍
# ============================================================================
# 动机：元数据里的 ScientificName 是提交者填的，可能张冠李戴、只填到属、
# 或干脆填成另一个物种。而 STAT / GSA 的分类是**从 reads 直接测出来的**，
# 两者对拍能自动揪出可疑记录 —— 这是「植物物种鉴定」这一半价值。
#
# 判定规则（保守，宁可不下结论也不要误报）：
#   MATCH    全名相等，或一方包含另一方（"Lycium barbarum" ⊂ "Lycium barbarum x ..."）
#   PARTIAL  同属不同种（近缘种 / 杂交种），且占比达到证据门槛
#   MISMATCH 属名都不同，且占比达到证据门槛 —— 元信息很可能记错了
#   COARSE   实测物种只到「属以上」阶元（如 Solanaceae 科）
#            —— reads 分辨率不足，**不能**据此判断元信息对错
#   WEAK     实测物种占比低于证据门槛 —— 多半是污染或数据库最近邻，证据不足
#   UNKNOWN  任一侧为空（无分类数据 / 元信息缺物种），不下结论
#
# ⚠ COARSE / WEAK 是实测倒逼出来的：最初只有 MATCH/PARTIAL/MISMATCH 三态时，
#   470 条真实数据报出 303 条 MISMATCH —— 其中 229 条主要物种是 "Solanaceae"
#   （枸杞所属的科，即宿主自己的祖先节点），40 条是占比 0.36% 的近缘种噪声。
#   把「分辨率不足」当成「元信息错误」会直接毁掉这个功能的可信度。
_GENUS_RE = re.compile(r"([A-Za-z][A-Za-z\-]{2,})")

# 这些不是学名属名，出现时说明 observed 侧不是正经物种名（如 GSA 的
# "Viruses (未细分)"、STAT 的 "unidentified"），此时应判 UNKNOWN 而非 MISMATCH
_GENUS_STOP = {
    "viruses", "virus", "viroid", "viroids", "unidentified", "unclassified",
    "environmental", "uncultured", "cellular", "organisms", "other", "synthetic",
    "unknown", "unassigned", "mixed", "samples", "metagenome",
}

# 高于「属」的分类阶元 —— observed 落在这些阶元时属于分辨率不足
_COARSE_RANKS = {"family", "subfamily", "order", "suborder", "class",
                 "phylum", "kingdom", "domain", "high"}

# 学名后缀 -> 阶元（植物 / 动物 / 微生物 / 病毒的常见命名规范）
_RANK_SUFFIX = (
    ("viridae", "family"), ("virinae", "subfamily"), ("virales", "order"),
    ("aceae", "family"), ("ales", "order"), ("idae", "family"),
    ("inae", "subfamily"), ("phyceae", "class"), ("opsida", "class"),
    ("phyta", "phylum"), ("mycota", "phylum"), ("bacteria", "class"),
)

# 单名且无标准后缀的高阶类群（Eukaryota / Pentapetalae / Proteobacteria ...）
_HIGH_RANK_NAMES = {
    "eukaryota", "bacteria", "archaea", "viruses", "metazoa", "animalia",
    "plantae", "viridiplantae", "streptophyta", "embryophyta", "tracheophyta",
    "spermatophyta", "magnoliophyta", "magnoliopsida", "eudicots", "rosids",
    "asterids", "lamiids", "campanulids", "pentapetalae", "superasterids",
    "superrosids", "fabids", "malvids", "gunneridae", "mesangiospermae",
    "proteobacteria", "firmicutes", "actinobacteria", "bacteroidetes",
    "cyanobacteria", "ascomycota", "basidiomycota", "mucoromycota",
    "arthropoda", "chordata", "mammalia", "insecta", "nematoda",
}

# 表示「未鉴定到种」的种加词
_INFRASUB_RE = re.compile(r"\b(?:subgen|subg|subsp|ssp|var|cf|aff|sp|spp)\.", re.I)


def genus_of(name: str) -> str:
    """取拉丁学名的属名（第一个词）。

    两条规则用来过滤非学名字符串，避免产生假的 MISMATCH：
      * 方括号是 NCBI 的「分类位置不确定」标记（'[Nicotiana] tabacum'），
        属名仍在括号内，所以只去符号不去内容；
      * 学名属名首字母必大写，且不在停用词表里（'viruses' 之类不算属名）。
    """
    if not name:
        return ""
    s = re.sub(r"\([^)]*\)", " ", name)   # 去掉 "(未细分)" 这类括注
    s = re.sub(r"[\[\]]", " ", s)         # 去方括号符号，保留 [Nicotiana] 的属名
    m = _GENUS_RE.search(s)
    if not m:
        return ""
    g = m.group(1)
    if not g[0].isupper():
        return ""
    gl = g.lower()
    return "" if gl in _GENUS_STOP else gl


def _norm_name(s: str) -> str:
    """规范化学名用于比较：去括注 / 方括号、压缩空白、转小写。

    这样 'Nicotiana tabacum' 与 NCBI 的 '[Nicotiana] tabacum' 能对上。
    """
    s = re.sub(r"\([^)]*\)", " ", s or "")
    s = re.sub(r"[\[\]]", " ", s)
    return re.sub(r"\s+", " ", s).strip().lower()


def rank_guess(name: str) -> str:
    """从学名字形猜测分类阶元（轻量启发式，不依赖 taxdump）。

    返回 species / genus / family / order / class / phylum / high /
    below_genus / unknown。仅用于判断 observed 是否**细到足以对拍物种**，
    不做严谨分类学断言。
    """
    n = _norm_name(name)
    if not n:
        return "unknown"
    if _INFRASUB_RE.search(n):
        return "below_genus"          # "Solanum subgen. Lycopersicon"
    if len(n.split()) >= 2:
        return "species"              # "Lycium barbarum"
    for suf, rk in _RANK_SUFFIX:
        if n.endswith(suf):
            return rk                 # "Solanaceae" -> family
    if n in _HIGH_RANK_NAMES:
        return "high"                 # "Eukaryota"
    return "genus"                    # "Lycium"（单名无后缀，多半是属）


def species_match(expected: str, observed: str, observed_pct: float = 0.0,
                  min_pct: float = 5.0) -> str:
    """判定元信息物种 (expected) 与实测主要物种 (observed) 是否一致。

    observed_pct 是实测主要物种在样本中的占比；低于 min_pct 视为证据不足，
    不据此断言元信息有误（低占比多半是污染或数据库最近邻）。
    """
    e, o = (expected or "").strip(), (observed or "").strip()
    if not e or not o:
        return "UNKNOWN"
    el, ol = _norm_name(e), _norm_name(o)
    if not el or not ol:
        return "UNKNOWN"
    # 1) 直接一致（含一方包含另一方：'Lycium barbarum' ⊂ 'Lycium barbarum x ...'）
    if el == ol or el in ol or ol in el:
        return "MATCH"
    # 2) 实测只到属以上阶元 -> 分辨率不足，不能据此判断元信息对错
    if rank_guess(o) in _COARSE_RANKS:
        return "COARSE"
    # 3) 实测占比过低 -> 证据不足
    if observed_pct < min_pct:
        return "WEAK"
    # 4) 属级比较
    ge, go = genus_of(e), genus_of(o)
    if ge and go:
        return "PARTIAL" if ge == go else "MISMATCH"
    return "UNKNOWN"


def organism_confidence(match: str, pct: float) -> str:
    """对拍结论的可信度 —— 主要物种占比越高，结论越硬。

    注意这是「对拍结论」的置信度，不是「物种鉴定」的置信度：
    若 MISMATCH 且实测物种占比很高，说明元信息记错的可能性也很大。
    """
    if match in ("UNKNOWN", "COARSE", "WEAK"):
        return "NA"
    if pct >= 20.0:
        return "HIGH"
    if pct >= 5.0:
        return "MEDIUM"
    return "LOW"


def virus_signal_level(pct_ident: float) -> str:
    """病毒信号强度。用归一化口径 virus_pct_ident 分级，跨来源可比。

    阈值取 0.01%（万分之一）：低于此值多半是 k-mer 误命中的背景噪声，
    但**不丢弃**，只降级为 WEAK 交给下游判断。
    """
    if pct_ident <= 0:
        return "NONE"
    if pct_ident < 0.01:
        return "WEAK"
    if pct_ident < 1.0:
        return "MODERATE"
    return "STRONG"


# ============================================================================
# 输入读取
# ============================================================================
def fmt_cell(v: Any) -> Any:
    """TSV 里把浮点写成定点小数，避免 4e-06 这种科学计数法影响可读性。"""
    if isinstance(v, float):
        s = f"{v:.8f}".rstrip("0").rstrip(".")
        return s if s not in ("", "-0") else "0"
    return v


def read_accessions(args: argparse.Namespace) -> List[str]:
    accs: List[str] = []
    if args.runs:
        with open(args.runs, "r", encoding="utf-8-sig", errors="replace") as f:
            for line in f:
                s = line.strip().split()[0] if line.strip() else ""
                if s and not s.startswith("#"):
                    accs.append(s)
    if args.csv:
        with open(args.csv, "r", encoding="utf-8-sig", errors="replace") as f:
            for rec in csv.DictReader(f):
                v = (rec.get(args.run_col) or "").strip()
                if v:
                    accs.append(v)
    if args.acc:
        accs.extend(args.acc)

    # 去重保序
    seen: Set[str] = set()
    out: List[str] = []
    for a in accs:
        a = a.strip()
        if a and a not in seen:
            seen.add(a)
            out.append(a)
    return out


def read_meta_species(path: str, run_col: str = "Run",
                      species_col: str = "ScientificName",
                      taxid_col: str = "TaxID"
                      ) -> Tuple[Dict[str, str], Dict[str, str]]:
    """从元数据表读 Run -> ScientificName / TaxID 映射，供物种对拍。

    兼容 pipeline 的 Core14（Global_Unified_Metadata_Core14.csv，列名固定），
    也接受任意含 Run 列的表。**列名找不到时返回空映射而不报错** ——
    对拍是增强项，不该让主流程挂掉。
    """
    sp: Dict[str, str] = {}
    tx: Dict[str, str] = {}
    if not path or not os.path.exists(path):
        return sp, tx
    with open(path, "r", encoding="utf-8-sig", errors="replace") as f:
        rd = csv.DictReader(f)
        if not rd.fieldnames:
            return sp, tx
        names = {n.strip().lstrip("\ufeff"): n for n in rd.fieldnames if n}
        rc = names.get(run_col) or names.get("Run") or rd.fieldnames[0]
        sc = names.get(species_col)
        tc = names.get(taxid_col)
        for rec in rd:
            r = (rec.get(rc) or "").strip()
            if not r:
                continue
            if sc and (rec.get(sc) or "").strip():
                sp[r] = (rec.get(sc) or "").strip()
            if tc and (rec.get(tc) or "").strip():
                tx[r] = (rec.get(tc) or "").strip()
    return sp, tx


# ============================================================================
# 主流程
# ============================================================================
def main() -> int:
    p = argparse.ArgumentParser(
        description="SRA / GSA 病毒预筛选器（基于 NCBI STAT 与 CNCB GSA 预计算分类结果）",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
双通路自动路由
--------------
按 accession 前缀自动分流，无需手工指定：

  SRR / ERR / DRR  ->  NCBI STAT
      GET trace.ncbi.nlm.nih.gov/Traces/sra-db-be/run_taxonomy?acc=X&cluster_name=public
      扁平表，带 tax_id；分母 analysed（已鉴定 reads）

  CRR              ->  CNCB GSA
      GET ngdc.cncb.ac.cn/gsa/taxonomyAnalysis/treePath/{CRA}/{CRR}/{CRR}.json
      嵌套树，**无 tax_id**（只能按名字判病毒）；分母是**总 reads**

  => 两者分母不同，直接比较会低估 GSA。故输出多一列 virus_pct_ident
     （病毒命中 / 已鉴定命中），跨来源比较请用这一列。

覆盖率（实测 2026-10）
---------------------
* NCBI：SRR/ERR/DRR 基本都有
* GSA ：**只有一部分 CRR 做过分类分析**。实测 148 条 CRR 中 60 条有数据、88 条没有。
        缺失**主要但不完全**按项目聚集：10 个项目里 8 个整批没有，
        1 个（CRA007417）是项目内部分缺失（21 条中 12 条无），2 个全覆盖。

无数据的 run 会被单列出来（unresolved.tsv / unresolved.list），
状态区分 NO_ANALYSIS（确认没做分析）、NOT_FOUND（accession 不在 GSA）、
FAILED（临时故障，**不缓存**，重跑即可补）。

* 想要「反向查询」——即"全库里有谁含某病毒"，本接口做不到。
  那需要 BigQuery / Athena 上的 sra_tax_analysis_tool 表做 SQL 过滤
  （Google Cloud 在国内需代理；AWS Athena 可达但需账号）。

* 两者都是 **k-mer 快速分类**，不是严谨比对。低丰度病毒可能漏检，
  也可能把宿主内源病毒序列（EVM）算进去。它的定位是 **便宜的预筛**，
  结论请用下游比对/组装复核。
""",
    )
    p.add_argument("--runs", help="accession 列表文件（每行一个）")
    p.add_argument("--csv", help="含 Run 列的 CSV/TSV")
    p.add_argument("--run-col", default="Run", help="CSV 中的 Run 列名（默认 Run）")
    p.add_argument("--acc", nargs="*", help="直接给 accession")
    p.add_argument("-o", "--outdir", default="./stat_prescreen_out", help="输出目录")
    p.add_argument("-j", "--jobs", type=int, default=6, help="并发数（默认 6，勿过高）")
    p.add_argument("--gsa-jobs", type=int, default=3,
                   help="GSA 并发数（NGDC 有反爬，默认 3，勿调高）")
    p.add_argument("--min-interval", type=float, default=0.0,
                   help="NCBI 最小请求间隔秒（0=不限）")
    p.add_argument("--gsa-min-interval", type=float, default=0.4,
                   help="GSA 最小请求间隔秒（默认 0.4，太低会 429）")
    p.add_argument("--retries", type=int, default=4, help="失败重试次数")
    p.add_argument("--timeout", type=int, default=90, help="单请求超时秒")
    p.add_argument("--virus-pct", type=float, default=None,
                   help="病毒占比阈值(%%); 达到则标记 pass 并导出 prescreen_pass.tsv")
    p.add_argument("--virus-pct-ident", type=float, default=None,
                   help="归一化口径阈值(%%): 病毒命中/已鉴定命中。**推荐用这个**，"
                        "跨 NCBI/GSA 来源可比；同时驱动 decision 列")
    p.add_argument("--core-csv", default=None,
                   help="元数据表（如 Global_Unified_Metadata_Core14.csv），"
                        "读其中 Run/ScientificName/TaxID 作为物种对拍基准")
    p.add_argument("--expect-species", default=None,
                   help="期望物种拉丁学名（无 --core-csv 时的对拍基准，如 'Lycium barbarum'）")
    p.add_argument("--taxids", help="自定义目标 taxid 列表文件（每行一个），替代 Viruses(10239)")
    p.add_argument("--taxdump", help="NCBI taxdump 目录（含 nodes.dmp），用于精确子树判定")
    p.add_argument("--cra-map", help="CRR->CRA 映射文件（TSV 两列），可跳过 GSA 网页解析")
    p.add_argument("--no-gsa", action="store_true", help="禁用 GSA 通路，只跑 NCBI")
    p.add_argument("--recheck-nodata", action="store_true",
                   help="忽略「无数据」标记，重新问一次 GSA（GSA 会陆续补做分析，可定期复查）")
    p.add_argument("--no-cache", action="store_true", help="忽略已有 JSON 缓存，强制重抓")
    args = p.parse_args()

    accs = read_accessions(args)
    if not accs:
        UI.err("没有读到任何 accession。用 --runs / --csv / --acc 指定输入。")
        return 2

    os.makedirs(args.outdir, exist_ok=True)
    raw_dir = os.path.join(args.outdir, "stat_raw")
    gsa_raw_dir = os.path.join(args.outdir, "gsa_raw")

    tax: Optional[Taxonomy] = None
    if args.taxdump:
        try:
            tax = Taxonomy(args.taxdump)
            UI.ok(f"已加载 taxonomy: {len(tax.parent):,} 个节点（精确病毒子树判定）")
        except Exception as e:
            UI.warn(f"taxdump 加载失败({e})，回退到名称启发式")

    viral_taxids: Set[int] = {TAXID_VIRUSES}
    if args.taxids:
        with open(args.taxids, "r", encoding="utf-8-sig") as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith("#"):
                    try:
                        viral_taxids.add(int(line.split()[0]))
                    except ValueError:
                        pass
        UI.info(f"目标 taxid 集合: {sorted(viral_taxids)[:10]}{' ...' if len(viral_taxids) > 10 else ''}")

    # ---- 按前缀路由 ----
    gsa_accs = [a for a in accs if a.upper().startswith("CRR")]
    stat_accs = [a for a in accs if a not in set(gsa_accs)]
    if args.no_gsa:
        stat_accs, gsa_accs = accs, []

    UI.banner("SRA / GSA 病毒预筛选")
    UI.info(f"待处理 run 数 : {len(accs):,}  (NCBI STAT {len(stat_accs)} / GSA {len(gsa_accs)})")
    UI.info(f"并发          : NCBI {args.jobs} 线程 (间隔 {args.min_interval}s) | "
            f"GSA {args.gsa_jobs} 线程 (间隔 {args.gsa_min_interval}s)")
    UI.info(f"输出目录      : {os.path.abspath(args.outdir)}")
    if not tax and stat_accs:
        UI.warn("未提供 --taxdump：NCBI 病毒明细用名称启发式（够用但不严谨）")
    if gsa_accs:
        UI.info("GSA 通路：无 tax_id，病毒按名字判定；分母是总 reads（见 virus_pct_ident）")

    cra_map: Dict[str, str] = {}
    if args.cra_map and os.path.exists(args.cra_map):
        with open(args.cra_map, "r", encoding="utf-8-sig", errors="replace") as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 2:
                    cra_map[parts[0]] = parts[1]
        UI.ok(f"已加载 CRR->CRA 映射 {len(cra_map)} 条")

    # ---- 物种对拍基准：逐 run 的元信息优先，其次全局 --expect-species ----
    meta_species: Dict[str, str] = {}
    meta_taxid: Dict[str, str] = {}
    if args.core_csv:
        meta_species, meta_taxid = read_meta_species(args.core_csv)
        if meta_species:
            UI.ok(f"已加载元信息物种 {len(meta_species):,} 条 "
                  f"（{os.path.basename(args.core_csv)}）-> 启用物种对拍")
        else:
            UI.warn(f"{os.path.basename(args.core_csv)} 里没读到 ScientificName 列，"
                    f"物种对拍将退化为 UNKNOWN")
    if not meta_species and args.expect_species:
        UI.info(f"物种对拍基准（全局）: {args.expect_species}")
    if not meta_species and not args.expect_species:
        UI.warn("未提供 --core-csv / --expect-species：跳过物种对拍"
                "（organism_match 全为 UNKNOWN）")

    fetcher = Fetcher(raw_dir, retries=args.retries, timeout=args.timeout,
                      min_interval=args.min_interval)
    gsa = GsaFetcher(gsa_raw_dir, retries=args.retries, timeout=args.timeout,
                     min_interval=args.gsa_min_interval, cra_map=cra_map)
    gsa_sem = threading.Semaphore(max(1, args.gsa_jobs))

    summary_rows: List[Dict[str, Any]] = []
    viral_rows_all: List[Dict[str, Any]] = []
    status_count: Dict[str, int] = {}
    t_start = time.time()
    done = 0

    gsa_set = set(gsa_accs)

    def work(acc: str) -> Tuple[str, str, Optional[Dict[str, Any]], str, float, str, str]:
        if acc in gsa_set:
            with gsa_sem:
                a, payload, status, dt, cra, detail = gsa.fetch(
                    acc, use_cache=not args.no_cache, recheck_nodata=args.recheck_nodata)
            return a, "GSA", payload, status, dt, detail, cra
        a, payload, status, dt = fetcher.fetch(acc, use_cache=not args.no_cache)
        return a, "NCBI_STAT", payload, status, dt, "", ""

    with ThreadPoolExecutor(max_workers=max(1, args.jobs)) as ex:
        futs = {ex.submit(work, a): a for a in accs}
        for fut in as_completed(futs):
            acc, source, payload, status, dt, detail, cra = fut.result()
            done += 1
            key = f"{source}:{status.split(' ')[0]}"
            status_count[key] = status_count.get(key, 0) + 1

            if payload is None:
                summary_rows.append({
                    "acc": acc, "source": source, "status": status, "detail": detail,
                    "cra": cra,
                    "analysed": 0, "identified": 0, "identified_pct": 0.0,
                    "n_taxa": 0, "virus_count": 0, "virus_pct": 0.0, "virus_pct_ident": 0.0,
                    "top_virus": "", "top_virus_pct": 0.0,
                    "top_organism": "", "top_organism_pct": 0.0,
                    "fetch_s": round(dt, 2), "pass": "",
                })
            else:
                r = (parse_gsa_tree(payload) if source == "GSA"
                     else parse_payload(payload, tax, viral_taxids))
                vr = r.pop("_viral_rows", [])
                for x in vr:
                    viral_rows_all.append({"acc": acc, "source": source, **x})
                r["acc"] = acc
                r["source"] = source
                r["status"] = status
                r["detail"] = detail
                r["cra"] = cra
                r["fetch_s"] = round(dt, 2)
                if args.virus_pct is not None:
                    r["pass"] = "PASS" if r["virus_pct"] >= args.virus_pct else ""
                else:
                    r["pass"] = ""
                summary_rows.append(r)

            if done % 25 == 0 or done == len(accs):
                el = time.time() - t_start
                rate = done / el if el > 0 else 0
                eta = (len(accs) - done) / rate if rate > 0 else 0
                print(f"  {UI.C['gray']}[{done}/{len(accs)}] {rate:.2f} run/s "
                      f"ETA {eta/60:.1f} min{UI.C['reset']}")

    summary_rows.sort(key=lambda x: x["acc"])

    # ==================== 规范化派生字段 ====================
    # 索引病毒明细，供主表回填「检出几个病毒类群 / 主要病毒的 taxid」
    viral_n: Dict[str, int] = {}
    viral_taxid_of: Dict[Tuple[str, str], str] = {}
    for x in viral_rows_all:
        viral_n[x["acc"]] = viral_n.get(x["acc"], 0) + 1
        viral_taxid_of.setdefault((x["acc"], x["org"]), x.get("tax_id", ""))

    thr_ident, thr_pct = args.virus_pct_ident, args.virus_pct

    for r in summary_rows:
        # --- 状态：拆成「机器可读枚举」+「人读详情」，信息不丢 ---
        st_raw = r.get("status") or ""
        st = st_raw.split(" ")[0] or "UNKNOWN"
        if st_raw != st:
            extra = st_raw[len(st):].strip().strip("()")
            if extra:
                r["detail"] = f"{r['detail']}; {extra}" if r.get("detail") else extra
        r["status"] = st
        r["retryable"] = "yes" if st == "FAILED" else "no"

        # --- 植物物种鉴定：元信息 vs 实测 ---
        msp = meta_species.get(r["acc"], "") or (args.expect_species or "")
        r["meta_species"] = msp
        r["meta_taxid"] = meta_taxid.get(r["acc"], "")
        obs = r.get("top_organism", "") or ""
        obs_pct = float(r.get("top_organism_pct") or 0)
        mt = species_match(msp, obs, obs_pct)
        r["organism_match"] = mt
        r["organism_conf"] = organism_confidence(mt, obs_pct)

        # --- 病毒预鉴定 ---
        vi = float(r.get("virus_pct_ident") or 0)
        vp = float(r.get("virus_pct") or 0)
        r["virus_signal"] = virus_signal_level(vi) if st == "OK" else "NA"
        r["n_viral_taxa"] = viral_n.get(r["acc"], 0)
        r["top_virus_taxid"] = viral_taxid_of.get(
            (r["acc"], r.get("top_virus", "")), "")
        r["pass_ident"] = "PASS" if (thr_ident is not None and vi >= thr_ident) else ""

        # --- 决策 ---
        # 阈值优先级：归一化口径（跨来源可比）> 旧口径 > 仅看是否 >0
        if thr_ident is not None:
            passed, thr_note = vi >= thr_ident, f"≥{thr_ident}%(已鉴定口径)"
        elif thr_pct is not None:
            passed, thr_note = vp >= thr_pct, f"≥{thr_pct}%"
        else:
            passed, thr_note = vi > 0, ">0%"

        if st != "OK":
            r["decision"] = "LOCAL_CHECK"
            r["decision_reason"] = "无预计算分类结果，需下载 reads 本地比对"
        elif not passed:
            r["decision"] = "SKIP"
            r["decision_reason"] = f"病毒信号未达阈值({thr_note})"
        elif mt == "MISMATCH":
            r["decision"] = "REVIEW"
            r["decision_reason"] = (f"检出病毒，但实测主要物种「{obs}」"
                                    f"与元信息「{msp}」不符，疑似元信息有误")
        else:
            r["decision"] = "DOWNLOAD"
            r["decision_reason"] = "检出病毒信号且物种对拍一致"

    # ---------------- 输出 ----------------
    # virus_pct      : 各来源自身口径（STAT=/analysed；GSA=/总 reads）
    # virus_pct_ident: 归一化口径（病毒命中/已鉴定命中），**跨来源比较用这一列**

    el = time.time() - t_start
    ok_rows = [r for r in summary_rows if r["status"] == "OK"]
    with_virus = [r for r in ok_rows if float(r["virus_pct_ident"] or 0) > 0]

    dec_count: Dict[str, int] = {}
    for r in summary_rows:
        dec_count[r["decision"]] = dec_count.get(r["decision"], 0) + 1
    mt_count: Dict[str, int] = {}
    for r in ok_rows:
        mt_count[r["organism_match"]] = mt_count.get(r["organism_match"], 0) + 1

    DEC_DESC = {
        "DOWNLOAD":    "检出病毒信号且物种对拍一致 —— 值得下载",
        "REVIEW":      "检出病毒但主要物种与元信息不符 —— 仍下载，人工复核",
        "SKIP":        "有分类数据但无病毒信号 —— 可不下",
        "LOCAL_CHECK": "拿不到预计算分类 —— 只能下载 reads 本地比对",
    }

    # ---- 1) 主表：规范 schema，一表回答「物种对不对 + 有没有病毒 + 该怎么办」----
    DEC_COLS = [
        "acc", "source", "cra", "status", "retryable", "detail",
        "meta_species", "meta_taxid",
        "top_organism", "top_organism_pct", "organism_match", "organism_conf",
        "virus_signal", "virus_pct", "virus_pct_ident",
        "top_virus", "top_virus_taxid", "top_virus_pct", "n_viral_taxa",
        "analysed", "identified", "identified_pct", "n_taxa",
        "decision", "decision_reason", "pass", "pass_ident", "fetch_s",
    ]
    dec_path = os.path.join(args.outdir, "sample_decision.tsv")
    with open(dec_path, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=DEC_COLS, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in summary_rows:
            w.writerow({k: fmt_cell(v) for k, v in r.items()})

    # ---- 2) 兼容表：列与旧版完全一致，老下游脚本无需改动 ----
    cols = ["acc", "source", "cra", "status", "detail", "virus_pct", "virus_pct_ident",
            "virus_count", "top_virus", "top_virus_pct",
            "analysed", "identified", "identified_pct", "n_taxa",
            "top_organism", "top_organism_pct", "pass", "fetch_s"]
    sum_path = os.path.join(args.outdir, "stat_run_summary.tsv")
    with open(sum_path, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in summary_rows:
            w.writerow({k: fmt_cell(v) for k, v in r.items()})

    # ---- 3) 病毒明细长表 ----
    vir_path = os.path.join(args.outdir, "stat_viral_taxa.tsv")
    with open(vir_path, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["acc", "source", "tax_id", "org", "path",
                                          "self_count", "total_count", "pct"],
                           delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in viral_rows_all:
            w.writerow({k: fmt_cell(v) for k, v in r.items()})

    # ---- 4) 派生清单：每行一个 accession，可直接喂下游 ----
    def _write_list(fname: str, rows_: List[Dict[str, Any]]) -> str:
        p_ = os.path.join(args.outdir, fname)
        with open(p_, "w", encoding="utf-8") as fh:
            for r_ in sorted(rows_, key=lambda x: x["acc"]):
                fh.write(r_["acc"] + "\n")
        return p_

    dl_rows = [r for r in summary_rows if r["decision"] in ("DOWNLOAD", "REVIEW")]
    lc_rows = [r for r in summary_rows if r["decision"] == "LOCAL_CHECK"]
    sk_rows = [r for r in summary_rows if r["decision"] == "SKIP"]
    dl_list = _write_list("prescreen_run.list", dl_rows)
    lc_list = _write_list("prescreen_local_check.list", lc_rows)
    sk_list = _write_list("prescreen_skip.list", sk_rows)

    # ---- 5) 通过阈值清单（按归一化口径降序，跨来源可比）----
    pass_path = None
    if thr_ident is not None or thr_pct is not None:
        passed = [r for r in summary_rows
                  if r.get("pass_ident") == "PASS" or r.get("pass") == "PASS"]
        pass_path = os.path.join(args.outdir, "prescreen_pass.tsv")
        with open(pass_path, "w", encoding="utf-8", newline="") as f:
            w = csv.DictWriter(f, fieldnames=DEC_COLS, delimiter="\t", extrasaction="ignore")
            w.writeheader()
            for r in sorted(passed, key=lambda x: -float(x["virus_pct_ident"] or 0)):
                w.writerow({k: fmt_cell(v) for k, v in r.items()})

    # ---- 6) 未决清单：没有分类数据的 run，需要走本地比对路线 ----
    unresolved = [r for r in summary_rows if r["status"] != "OK"]
    unrel_path = os.path.join(args.outdir, "unresolved.tsv")
    with open(unrel_path, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["acc", "source", "cra", "status", "detail"],
                           delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in sorted(unresolved, key=lambda x: (x["source"], x["acc"])):
            w.writerow({k: fmt_cell(v) for k, v in r.items()})
    unrel_list = os.path.join(args.outdir, "unresolved.list")
    with open(unrel_list, "w", encoding="utf-8") as f:
        for r in sorted(unresolved, key=lambda x: x["acc"]):
            f.write(r["acc"] + "\n")

    # ---- 7) 运行元信息：参数快照 + 输入指纹，保证可复现与可追溯 ----
    acc_fp = hashlib.sha1("\n".join(sorted(accs)).encode("utf-8")).hexdigest()[:12]
    if meta_species:
        baseline = f"core_csv:{os.path.basename(args.core_csv)}({len(meta_species)})"
    elif args.expect_species:
        baseline = f"expect_species:{args.expect_species}"
    else:
        baseline = "none"
    meta = {
        "script": "stat_prescreen.py",
        "script_version": SCRIPT_VERSION,
        "schema_version": SCHEMA_VERSION,
        "generated_at": datetime.now().isoformat(timespec="seconds"),
        "elapsed_min": round(el / 60, 2),
        "inputs": {
            "runs_file": args.runs, "csv_file": args.csv, "core_csv": args.core_csv,
            "n_acc": len(accs), "acc_fingerprint": acc_fp,
        },
        "params": {
            "jobs": args.jobs, "gsa_jobs": args.gsa_jobs,
            "min_interval": args.min_interval, "gsa_min_interval": args.gsa_min_interval,
            "retries": args.retries, "timeout": args.timeout,
            "virus_pct": args.virus_pct, "virus_pct_ident": args.virus_pct_ident,
            "taxdump": args.taxdump, "no_gsa": args.no_gsa,
        },
        "judgement": {
            "viral_taxid_root": TAXID_VIRUSES,
            "ncbi_method": "taxdump_subtree" if tax else "name_heuristic",
            "gsa_method": "name_heuristic",
            "species_baseline": baseline,
        },
        "counts": {
            "total": len(accs), "ok": len(ok_rows), "with_virus": len(with_virus),
            "decision": dec_count, "status": status_count, "organism_match": mt_count,
        },
    }
    meta_path = os.path.join(args.outdir, "prescreen_meta.json")
    with open(meta_path, "w", encoding="utf-8") as f:
        json.dump(meta, f, ensure_ascii=False, indent=2)

    # ---- 8) 人读摘要报告 ----
    rep = os.path.join(args.outdir, "STAT_PRESCREEN_REPORT.md")
    with open(rep, "w", encoding="utf-8") as f:
        f.write("# SRA / GSA 病毒预筛选 & 物种鉴定报告\n\n")
        f.write(f"- 脚本版本: stat_prescreen.py v{SCRIPT_VERSION} (schema {SCHEMA_VERSION})\n")
        f.write(f"- 生成时间: {datetime.now().isoformat(timespec='seconds')}\n")
        f.write(f"- 输入 run 数: {len(accs)}  （指纹 {acc_fp}）\n")
        f.write(f"- 有分类数据: {len(ok_rows)}\n")
        f.write(f"- 检出病毒(>0%): {len(with_virus)}\n")
        f.write(f"- 状态分布: {status_count}\n")
        f.write(f"- 无分类数据: {len(unresolved)}（见 unresolved.tsv）\n")
        f.write(f"- NCBI 病毒判定: {'taxdump 精确子树' if tax else '名称启发式(近似)'}\n")
        f.write("- GSA 病毒判定: 名称判定（GSA 分类树不含 tax_id）\n")
        f.write(f"- 物种对拍基准: {baseline}\n\n")

        f.write("## 决策分布\n\n")
        f.write("| decision | 条数 | 含义 |\n|---|---|---|\n")
        for k in ("DOWNLOAD", "REVIEW", "SKIP", "LOCAL_CHECK"):
            f.write(f"| {k} | {dec_count.get(k, 0)} | {DEC_DESC[k]} |\n")
        f.write("\n")

        if baseline != "none":
            f.write("## 植物物种鉴定（元信息 vs 实测主要物种）\n\n")
            f.write("| 判定 | 条数 | 含义 |\n|---|---|---|\n")
            for k, d in (("MATCH", "一致"),
                         ("PARTIAL", "同属不同种（近缘 / 杂交）"),
                         ("MISMATCH", "属名都不同且占比足够 —— 疑似元信息有误"),
                         ("COARSE", "实测只到属以上阶元（分辨率不足，不下结论）"),
                         ("WEAK", "实测占比过低（证据不足，不下结论）"),
                         ("UNKNOWN", "一侧缺失，不下结论")):
                f.write(f"| {k} | {mt_count.get(k, 0)} | {d} |\n")
            f.write("\n")
            mism = [r for r in ok_rows if r["organism_match"] == "MISMATCH"]
            if mism:
                f.write("### 疑似元信息记录有误（前 40）\n\n")
                f.write("| run | 元信息记录 | 实测主要物种 | 实测占比% | 置信度 | 主要病毒 |\n")
                f.write("|---|---|---|---|---|---|\n")
                for r in sorted(mism, key=lambda x: -float(x["top_organism_pct"] or 0))[:40]:
                    f.write(f"| {r['acc']} | {r['meta_species']} | {r['top_organism']} | "
                            f"{float(r['top_organism_pct'] or 0):.2f} | {r['organism_conf']} | "
                            f"{r['top_virus']} |\n")
                if len(mism) > 40:
                    f.write(f"\n_（共 {len(mism)} 条，此处仅列前 40）_\n")
                f.write("\n")

        if unresolved:
            f.write("## 无分类数据的 run（需本地比对）\n\n")
            f.write("这些 run 拿不到预计算分类结果，**不是分析失败**，"
                    "而是上游没有这份数据。要判断物种只能下载 reads 本地比对。\n\n")
            f.write(f"完整清单见 `unresolved.tsv`（含原因）与 `unresolved.list`"
                    f"（纯 accession，可直接喂下游）。共 {len(unresolved)} 条。\n\n")

            # 按项目聚合 —— 缺失主要按项目聚集，聚合后一眼看出该等哪个项目
            by_proj: Dict[str, List[str]] = {}
            for r in unresolved:
                by_proj.setdefault(r.get("cra") or "(未解析)", []).append(r["acc"])
            f.write("### 按项目（CRA）聚合\n\n")
            f.write("| CRA | 无数据 run 数 | 状态 |\n|---|---|---|\n")
            for k in sorted(by_proj, key=lambda x: -len(by_proj[x])):
                st = next((r["status"] for r in unresolved
                           if (r.get("cra") or "(未解析)") == k), "")
                f.write(f"| {k} | {len(by_proj[k])} | {st} |\n")
            f.write("\n")

            f.write("### 明细（前 40 条）\n\n")
            f.write("| run | 来源 | CRA | 状态 | 原因 |\n|---|---|---|---|---|\n")
            for r in sorted(unresolved, key=lambda x: (x["source"], x["acc"]))[:40]:
                f.write(f"| {r['acc']} | {r['source']} | {r.get('cra', '')} | "
                        f"{r['status']} | {r['detail']} |\n")
            if len(unresolved) > 40:
                f.write(f"\n_（共 {len(unresolved)} 条，此处仅列前 40）_\n")
            f.write("\n")

        f.write("## 口径警告\n\n")
        f.write("- `virus_pct` 是各来源**自身口径**：NCBI STAT 分母为 analysed（已鉴定 reads），\n")
        f.write("  GSA 分母为**总 reads**。直接比较会低估 GSA。\n")
        f.write("- 跨来源比较请用 `virus_pct_ident`（病毒命中 / 已鉴定命中）。\n")
        f.write("- k-mer 分类是**便宜的预筛**，会漏低丰度病毒、可能算进宿主内源病毒序列(EVM)，\n")
        f.write("  结论请用下游比对/组装复核。\n\n")
        f.write("## 病毒占比 Top 30\n\n")
        f.write("| run | 来源 | 病毒% | 病毒%(已鉴定口径) | 主要病毒 | 测序量 | 主要非病毒物种 |\n")
        f.write("|---|---|---|---|---|---|---|\n")
        for r in sorted(with_virus, key=lambda x: -float(x["virus_pct_ident"] or 0))[:30]:
            f.write(f"| {r['acc']} | {r['source']} | {float(r['virus_pct'] or 0):.4f} | "
                    f"{float(r['virus_pct_ident'] or 0):.4f} | {r['top_virus']} | "
                    f"{int(r['analysed'] or 0):,} | {r['top_organism']} |\n")

    # ---------------- 控制台小结 ----------------
    UI.stage("完成", "done")
    UI.ok(f"总耗时 {el/60:.1f} min，平均 {len(accs)/el:.2f} run/s")
    UI.ok(f"有分类数据: {len(ok_rows)}/{len(accs)}")
    for src in ("NCBI_STAT", "GSA"):
        sub = [r for r in summary_rows if r["source"] == src]
        if sub:
            okn = sum(1 for r in sub if r["status"] == "OK")
            UI.info(f"  {src}: {okn}/{len(sub)} 有数据")
    UI.ok(f"检出病毒(>0%): {len(with_virus)}")
    if fetcher.n_cache_corrupt:
        UI.warn(f"STAT 缓存有 {fetcher.n_cache_corrupt} 个不可用 (已重抓); "
                "若数量异常多请检查 raw_dir 是否被写坏")
    UI.info("决策分布: " + " | ".join(
        f"{k}={dec_count.get(k, 0)}" for k in ("DOWNLOAD", "REVIEW", "SKIP", "LOCAL_CHECK")))
    if baseline != "none":
        UI.info("物种对拍: " + " | ".join(
            f"{k}={mt_count.get(k, 0)}"
            for k in ("MATCH", "PARTIAL", "MISMATCH", "COARSE", "WEAK", "UNKNOWN")))
        if mt_count.get("MISMATCH"):
            UI.warn(f"{mt_count['MISMATCH']} 条实测物种与元信息不符 -> 见 sample_decision.tsv"
                    f" 的 organism_match=MISMATCH")
    UI.info(f"主表     -> {dec_path}")
    UI.info(f"下载清单 -> {dl_list}  ({len(dl_rows)} 条)")
    UI.info(f"本地比对 -> {lc_list}  ({len(lc_rows)} 条)")
    UI.info(f"可跳过   -> {sk_list}  ({len(sk_rows)} 条)")
    if pass_path:
        UI.info(f"通过清单 -> {pass_path}")
    UI.info(f"兼容表   -> {sum_path}")
    UI.info(f"病毒明细 -> {vir_path}")
    UI.info(f"元信息   -> {meta_path}")
    UI.info(f"报告     -> {rep}")

    if unresolved:
        UI.warn(f"{len(unresolved)} 个 run 没有分类数据 -> {unrel_path}")
        UI.info(f"  （纯 accession 列表: {unrel_list}，可直接喂给下游本地比对）")
        by_status: Dict[str, int] = {}
        for r in unresolved:
            k = f"{r['source']}:{r['status']}"
            by_status[k] = by_status.get(k, 0) + 1
        for k in sorted(by_status):
            UI.info(f"    {k}: {by_status[k]}")
        if any(r["status"] == "FAILED" for r in unresolved):
            UI.warn("其中有 FAILED —— 属于临时故障，**未写入缓存**，直接重跑即可补上")

    return 0


if __name__ == "__main__":
    sys.exit(main())
