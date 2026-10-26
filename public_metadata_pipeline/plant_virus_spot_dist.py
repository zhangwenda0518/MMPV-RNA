#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
plant_virus_spot_dist.py — 每个 run 的「植物病毒科」spot 分布
============================================================

数据来源与 `stat_prescreen.py` 同口径
------------------------------------
  SRR / ERR / DRR : NCBI STAT
      GET https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/run_taxonomy
          ?acc={acc}&cluster_name=public
      ⚠ 必须带 cluster_name=public，否则 HTTP 200 + 空 body。
  CRR             : CNCB GSA（NCBI STAT 不含 CRR）
      GET https://ngdc.cncb.ac.cn/gsa/taxonomyAnalysis/treePath/{CRA}/{CRR}/{CRR}.json

spot 口径（两条通路不可直接比较）
--------------------------------
  STAT  self_count  该 taxon **直接命中**的 spot 数（不含后代）。
                    对 tax_table 每一行沿 parent 链找「最近的植物病毒科祖先」，
                    按该科归组求和 -> 每个 run 每个科的 spot。
  GSA   hits        只有**含后代**的累计 hits；科级节点即该科总命中（不再累加子节点）。
                    ⚠ GSA 的 hits 是 reads 量级，STAT 的 spot 是 k-mer 量级。

植物病毒科名单
--------------
  默认读 `00_data/02_references/plant_virus_families_classified.tsv` 第一列
  （51 个科，含类病毒 Pospiviroidae / Avsunviroidae）。
  类病毒单独成列，因为它们的 spot 量级会主导分布。

输出
----
  per_sample.tsv          每 run 一行：植物病毒 spot 合计 / 类病毒 spot / 主要科
  per_sample_family.tsv   每 run × 科 长表
  DISTRIBUTION_REPORT.md  分位数、分档、按项目分组
  stat_raw/ gsa_raw/      原始 JSON 缓存（断点续跑）

用法
----
    python plant_virus_spot_dist.py --runs samples.txt -o ./pv_spot
    python plant_virus_spot_dist.py --runs a.txt --runs b.txt -o ./pv_spot --label onekp --label goji
    python plant_virus_spot_dist.py --runs samples.txt -o ./pv_spot --no-gsa
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import random
import re
import statistics
import sys
import threading
import time
import urllib.error
import urllib.parse
import urllib.request
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed

STAT_URL = ("https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/run_taxonomy"
            "?acc={acc}&cluster_name=public")
GSA_BASE = "https://ngdc.cncb.ac.cn"
GSA_SEARCH = GSA_BASE + "/gsa/search?searchTerm={acc}"
GSA_TREE = GSA_BASE + "/gsa/taxonomyAnalysis/treePath/{cra}/{acc}/{acc}.json"
GSA_CRA_NOISE = {"CRA000112"}
GSA_CACHE_V = 1
UA = "pv-spot-dist/1.0 (research)"

TAXID_ROOT = 1
TAXID_CELLULAR = 131567
TAXID_VIRUSES = 10239

# 类病毒两个科（NCBI family 级）
VIROID_FAMILIES = {"Pospiviroidae", "Avsunviroidae"}


# ============================================================================
# Taxonomy：nodes.dmp + names.dmp + merged.dmp
# ============================================================================
class Tax:
    def __init__(self, d: str) -> None:
        self.parent: dict[int, int] = {}
        self.rank: dict[int, str] = {}
        self.merged: dict[int, int] = {}
        nodes = os.path.join(d, "nodes.dmp")
        if not os.path.exists(nodes):
            raise FileNotFoundError(f"未找到 {nodes}")
        with open(nodes, "r", encoding="utf-8", errors="replace") as f:
            for line in f:
                p = line.split("\t|\t")
                if len(p) < 3:
                    continue
                try:
                    tid, par = int(p[0]), int(p[1])
                except ValueError:
                    continue
                self.parent[tid] = par
                self.rank[tid] = p[2].strip()
        mp = os.path.join(d, "merged.dmp")
        if os.path.exists(mp):
            with open(mp, "r", encoding="utf-8", errors="replace") as f:
                for line in f:
                    p = line.split("\t|\t")
                    if len(p) >= 2:
                        try:
                            self.merged[int(p[0])] = int(p[1])
                        except ValueError:
                            pass
        self._cache: dict[int, int | None] = {}

    def norm(self, tid: int) -> int:
        seen = set()
        while tid in self.merged and tid not in seen:
            seen.add(tid)
            tid = self.merged[tid]
        return tid

    def family_of(self, tid: int, fam_set: set[int]) -> int | None:
        """沿 parent 链找最近的、在 fam_set 里的祖先（含自身）。"""
        t = self.norm(tid)
        if t in self._cache:
            return self._cache[t]
        cur, hops, out = t, 0, None
        while hops < 80:
            if cur in fam_set:
                out = cur
                break
            p = self.parent.get(cur)
            if p is None or p == cur or cur == TAXID_ROOT:
                break
            cur = p
            hops += 1
        self._cache[t] = out
        return out



def load_family_taxids(taxdump: str, fam_names: set[str]) -> tuple[dict[str, int], set[int]]:
    """科名 -> taxid；返回 (name2tid, tid_set)。只认 scientific name。"""
    name2tid: dict[str, int] = {}
    names = os.path.join(taxdump, "names.dmp")
    with open(names, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            p = line.split("\t|\t")
            if len(p) < 4:
                continue
            nm = p[1].strip()
            if nm not in fam_names:
                continue
            if p[3].rstrip("\t|\n").strip() != "scientific name":
                continue
            try:
                name2tid.setdefault(nm, int(p[0]))
            except ValueError:
                pass
    return name2tid, set(name2tid.values())


# ============================================================================
# 抓取
# ============================================================================
class Fetcher:
    def __init__(self, raw_dir: str, retries: int = 4, timeout: int = 90) -> None:
        self.raw_dir = raw_dir
        self.retries = retries
        self.timeout = timeout
        os.makedirs(raw_dir, exist_ok=True)

    def path(self, acc: str) -> str:
        return os.path.join(self.raw_dir, f"{acc}.json")

    def fetch(self, acc: str) -> tuple[str, dict | None, str]:
        cp = self.path(acc)
        if os.path.exists(cp):
            try:
                txt = open(cp, "r", encoding="utf-8").read()
                if not txt.strip():
                    return acc, None, "NO_STAT"
                obj = json.loads(txt)
                return acc, (obj[0] if isinstance(obj, list) else obj), "OK"
            except Exception:
                pass
        url = STAT_URL.format(acc=urllib.parse.quote(acc))
        last = ""
        for attempt in range(self.retries):
            try:
                req = urllib.request.Request(url, headers={"User-Agent": UA})
                with urllib.request.urlopen(req, timeout=self.timeout) as r:
                    raw = r.read()
                if not raw:
                    open(cp, "w", encoding="utf-8").write("")
                    return acc, None, "NO_STAT"
                payload = json.loads(raw)[0]
                json.dump([payload], open(cp, "w", encoding="utf-8"), ensure_ascii=False)
                return acc, payload, "OK"
            except urllib.error.HTTPError as e:
                last = f"HTTP {e.code}"
                if e.code in (429, 500, 502, 503, 504):
                    time.sleep(min(2 ** attempt + random.random(), 30))
                    continue
                break
            except Exception as e:
                last = f"{type(e).__name__}: {e}"
                time.sleep(min(2 ** attempt + random.random(), 30))
        return acc, None, f"FAILED ({last})"


class GsaFetcher:
    def __init__(self, raw_dir: str, retries: int = 4, timeout: int = 60,
                 min_interval: float = 0.4) -> None:
        self.raw_dir = raw_dir
        self.retries = retries
        self.timeout = timeout
        self.min_interval = min_interval
        self._lock = threading.Lock()
        self._last = 0.0
        os.makedirs(raw_dir, exist_ok=True)

    def _throttle(self) -> None:
        with self._lock:
            gap = time.time() - self._last
            if gap < self.min_interval:
                time.sleep(self.min_interval - gap)
            self._last = time.time()

    def _get(self, url: str) -> str:
        last: Exception | None = None
        for attempt in range(self.retries):
            self._throttle()
            try:
                req = urllib.request.Request(url, headers={"User-Agent": UA})
                with urllib.request.urlopen(req, timeout=self.timeout) as r:
                    return r.read().decode("utf-8", errors="replace")
            except urllib.error.HTTPError as e:
                last = e
                if e.code in (429, 500, 502, 503, 504):
                    time.sleep(min(3 * (attempt + 1) + random.random(), 45))
                    continue
                raise
            except Exception as e:
                last = e
                time.sleep(min(2 ** attempt + random.random(), 30))
        raise last if last else RuntimeError("unknown")

    def path(self, acc: str) -> str:
        return os.path.join(self.raw_dir, f"{acc}.json")

    def fetch(self, acc: str) -> tuple[str, dict | None, str, str]:
        cp = self.path(acc)
        if os.path.exists(cp):
            try:
                obj = json.loads(open(cp, "r", encoding="utf-8").read() or "null")
                if isinstance(obj, dict) and obj.get("_v") == GSA_CACHE_V:
                    if obj.get("_status") == "OK" and isinstance(obj.get("_tree"), dict):
                        return acc, obj["_tree"], "OK", obj.get("_cra", "")
                    return acc, None, obj.get("_status", "NO_DATA"), obj.get("_cra", "")
            except Exception:
                pass
        try:
            html = self._get(GSA_SEARCH.format(acc=acc))
            cands = [c for c in sorted(set(re.findall(r"CRA\d{6}", html)))
                     if c not in GSA_CRA_NOISE]
        except Exception as e:
            return acc, None, f"FAILED (CRA {type(e).__name__})", ""
        if not cands:
            json.dump({"_v": GSA_CACHE_V, "_status": "NOT_FOUND", "_cra": ""},
                      open(cp, "w", encoding="utf-8"))
            return acc, None, "NOT_FOUND", ""
        transient = []
        for cra in cands:
            try:
                raw = self._get(GSA_TREE.format(cra=cra, acc=acc))
                if not raw.strip():
                    continue
                tree = json.loads(raw)
                if not isinstance(tree, dict):
                    continue
                json.dump({"_v": GSA_CACHE_V, "_status": "OK", "_cra": cra, "_tree": tree},
                          open(cp, "w", encoding="utf-8"), ensure_ascii=False)
                return acc, tree, "OK", cra
            except urllib.error.HTTPError as e:
                if e.code == 404:
                    continue
                transient.append(f"{cra}:HTTP{e.code}")
            except json.JSONDecodeError:
                continue
            except Exception as e:
                transient.append(f"{cra}:{type(e).__name__}")
        if transient:
            return acc, None, f"FAILED ({';'.join(transient)})", ""
        json.dump({"_v": GSA_CACHE_V, "_status": "NO_ANALYSIS",
                   "_cra": "/".join(cands)}, open(cp, "w", encoding="utf-8"))
        return acc, None, "NO_ANALYSIS", "/".join(cands)


# ============================================================================
# 解析
# ============================================================================
def stat_family_spots(payload: dict, tax: Tax, fam_set: set[int],
                      tid2name: dict[int, str]) -> tuple[dict[str, float], dict]:
    """STAT -> {科名: spot}，按最近植物病毒科祖先归组 self_count。"""
    tot = payload.get("tax_totals", {}) or {}
    table = payload.get("tax_table", []) or []
    analysed = int(tot.get("analysed") or tot.get("total") or 0)
    identified = int(tot.get("identified") or 0)

    fam_spots: dict[str, float] = defaultdict(float)
    for r in table:
        try:
            tid = int(r.get("tax_id"))
        except (TypeError, ValueError):
            continue
        if tid in (TAXID_ROOT, TAXID_CELLULAR, TAXID_VIRUSES):
            continue
        fam = tax.family_of(tid, fam_set)
        if fam is None:
            continue
        nm = tid2name.get(fam)
        if nm is None:
            continue
        try:
            sc = float(r.get("self_count") or 0)
        except (TypeError, ValueError):
            sc = 0.0
        if sc:
            fam_spots[nm] += sc

    meta = {"analysed": analysed, "identified": identified, "n_taxa": len(table)}
    return dict(fam_spots), meta


def _gsa_flat(node: dict, path: tuple[str, ...] = ()):
    p = path + (node.get("name", ""),)
    yield node, p
    for c in node.get("children") or []:
        yield from _gsa_flat(c, p)


GSA_ALL_VIRUSES = "__ALL_VIRUSES__"


def gsa_family_spots(tree: dict, fam_names: set[str]) -> tuple[dict[str, float], dict]:
    """GSA -> {科名: hits} + {__ALL_VIRUSES__: hits}。

    ⚠ GSA 的分类树**普遍只细到属**（实测 67 条 OK 里仅 5 条有 Viruses 节点，
    且最细是 `Viruses > Pospiviroid > Tomato chlorotic dwarf viroid`，
    没有科级节点），所以科级匹配大多是空的。真正能用的是 Viruses 根节点的
    累计 hits —— 它含后代，不能再累加子节点。
    """
    fam_spots: dict[str, float] = {}
    n_nodes = 0
    for node, _p in _gsa_flat(tree):
        n_nodes += 1
        nm = (node.get("name") or "").strip()
        try:
            hits = float(str(node.get("hits", 0)).strip() or 0)
        except (TypeError, ValueError):
            hits = 0.0
        if nm.lower() == "viruses" and hits > 0:
            fam_spots[GSA_ALL_VIRUSES] = max(fam_spots.get(GSA_ALL_VIRUSES, 0.0), hits)
        if nm not in fam_names:
            continue
        if hits > 0:
            fam_spots[nm] = fam_spots.get(nm, 0.0) + hits

    # 反推总 reads / 已鉴定 reads（GSA percent 分母是总 reads）
    ident_pct = 0.0
    try:
        ident_pct = float(str(tree.get("identified", "0")).rstrip("%"))
    except ValueError:
        pass
    total_reads = 0.0
    for node, _p in _gsa_flat(tree):
        try:
            h = float(str(node.get("hits", 0)).strip() or 0)
            p_ = float(str(node.get("percent", "0")).rstrip("%"))
        except ValueError:
            continue
        if h > 0 and p_ > 0:
            total_reads = max(total_reads, h / (p_ / 100.0))
    meta = {"analysed": int(round(total_reads)), "identified": int(round(total_reads * ident_pct / 100.0)),
            "n_taxa": n_nodes}
    return fam_spots, meta


# ============================================================================
# 分布统计
# ============================================================================
def quantiles(vals: list[float]) -> dict[str, float]:
    if not vals:
        return {}
    v = sorted(vals)
    n = len(v)

    def q(p: float) -> float:
        if n == 1:
            return v[0]
        i = (n - 1) * p
        lo, hi = int(i), min(int(i) + 1, n - 1)
        return v[lo] + (v[hi] - v[lo]) * (i - lo)

    return {
        "n": n, "min": v[0], "p05": q(0.05), "p25": q(0.25), "median": q(0.5),
        "p75": q(0.75), "p90": q(0.90), "p95": q(0.95), "p99": q(0.99),
        "max": v[-1], "mean": statistics.fmean(v),
    }


BINS = [(0, 0), (1, 9), (10, 49), (50, 99), (100, 499), (500, 999),
        (1000, 4999), (5000, 9999), (10000, None)]


def bin_label(lo: int, hi: int | None) -> str:
    return f"{lo}" if hi == lo else (f"{lo}-{hi}" if hi is not None else f">={lo}")


def bin_counts(vals: list[float]) -> list[tuple[str, int]]:
    out = []
    for lo, hi in BINS:
        if hi is None:
            c = sum(1 for x in vals if x >= lo)
        elif hi == lo:
            c = sum(1 for x in vals if x == lo)
        else:
            c = sum(1 for x in vals if lo <= x <= hi)
        out.append((bin_label(lo, hi), c))
    return out


def fmt_q(q: dict[str, float]) -> str:
    if not q:
        return "（无数据）"
    def f(x):
        return f"{x:,.0f}" if abs(x) >= 100 else f"{x:,.2f}"
    return (f"n={q['n']} | min={f(q['min'])} | p5={f(q['p05'])} | p25={f(q['p25'])} | "
            f"**median={f(q['median'])}** | p75={f(q['p75'])} | p90={f(q['p90'])} | "
            f"p95={f(q['p95'])} | p99={f(q['p99'])} | max={f(q['max'])}")


# ============================================================================
# 主流程
# ============================================================================
def main() -> int:
    ap = argparse.ArgumentParser(
        description="每个 run 的植物病毒科 spot 分布（NCBI STAT / CNCB GSA）")
    ap.add_argument("--runs", action="append", required=True,
                    help="样本列表文件，每行一个 acc；可多次给，配合 --label")
    ap.add_argument("--label", action="append", default=None,
                    help="与 --runs 一一对应的分组名（默认 file1/file2/...）")
    ap.add_argument("-o", "--outdir", required=True)
    ap.add_argument("--taxdump", default=os.path.expanduser(
        "~/plant_virus_mining/00_data/02_references/taxdump_20261021"))
    ap.add_argument("--families", default=os.path.expanduser(
        "~/plant_virus_mining/00_data/02_references/plant_virus_families_classified.tsv"))
    ap.add_argument("-j", "--jobs", type=int, default=6)
    ap.add_argument("--gsa-jobs", type=int, default=3)
    ap.add_argument("--gsa-min-interval", type=float, default=0.4)
    ap.add_argument("--no-gsa", action="store_true")
    ap.add_argument("--only-gsa", action="store_true")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    # ---- 分组 ----
    groups: list[tuple[str, list[str]]] = []
    for i, path in enumerate(args.runs):
        lbl = (args.label[i] if args.label and i < len(args.label)
               else os.path.splitext(os.path.basename(path))[0])
        accs = []
        with open(path, "r", encoding="utf-8-sig", errors="replace") as f:
            for line in f:
                s = line.strip().split()[0] if line.strip() else ""
                if s and not s.startswith("#"):
                    accs.append(s)
        groups.append((lbl, accs))

    # ---- 科名单 ----
    fam_names: set[str] = set()
    with open(args.families, "r", encoding="utf-8-sig", errors="replace") as f:
        for i, rec in enumerate(csv.reader(f, delimiter="\t")):
            if i == 0:
                continue
            if rec and rec[0].strip():
                fam_names.add(rec[0].strip())
    fam_names |= VIROID_FAMILIES
    print(f"[fam] 植物病毒科名单 {len(fam_names)} 个（含类病毒 2 个）")

    # ---- taxonomy ----
    print(f"[tax] 加载 {args.taxdump} ...", flush=True)
    t0 = time.time()
    tax = Tax(args.taxdump)
    print(f"[tax] {len(tax.parent):,} 节点，{time.time()-t0:.1f}s", flush=True)
    name2tid, fam_set = load_family_taxids(args.taxdump, fam_names)
    missing = sorted(fam_names - set(name2tid))
    print(f"[tax] 解析到科 taxid {len(name2tid)}/{len(fam_names)}")
    if missing:
        print(f"[tax] ⚠ taxdump 里找不到（按名字保留，仅 GSA 通路可用）：{', '.join(missing)}")
    tid2name = {v: k for k, v in name2tid.items()}

    # ---- 抓取 ----
    stat_accs, gsa_accs = [], []
    for lbl, accs in groups:
        for a in accs:
            (gsa_accs if a.upper().startswith("CRR") else stat_accs).append(a)
    if args.no_gsa:
        stat_accs += gsa_accs
        gsa_accs = []
    if args.only_gsa:
        gsa_accs += stat_accs
        stat_accs = []

    f_stat = Fetcher(os.path.join(args.outdir, "stat_raw"))
    f_gsa = GsaFetcher(os.path.join(args.outdir, "gsa_raw"),
                       min_interval=args.gsa_min_interval)

    print(f"[get] NCBI STAT {len(stat_accs)} | GSA {len(gsa_accs)}", flush=True)
    results: dict[str, dict] = {}
    lock = threading.Lock()
    done = 0
    total = len(stat_accs) + len(gsa_accs)

    def work_stat(a: str):
        acc, payload, status = f_stat.fetch(a)
        fam_spots, meta = ({}, {}) if payload is None else stat_family_spots(
            payload, tax, fam_set, tid2name)
        return acc, "NCBI_STAT", status, fam_spots, meta, ""

    def work_gsa(a: str):
        acc, tree, status, cra = f_gsa.fetch(a)
        fam_spots, meta = ({}, {}) if tree is None else gsa_family_spots(tree, fam_names)
        return acc, "GSA", status, fam_spots, meta, cra

    with ThreadPoolExecutor(max_workers=max(1, args.jobs)) as ex:
        futs = [ex.submit(work_stat, a) for a in stat_accs]
        futs += [ex.submit(work_gsa, a) for a in gsa_accs]
        for fut in as_completed(futs):
            acc, src, status, fam_spots, meta, cra = fut.result()
            with lock:
                results[acc] = {"source": src, "status": status, "fam": fam_spots,
                                "meta": meta, "cra": cra}
                done += 1
                if done % 25 == 0 or done == total:
                    print(f"  [{done}/{total}]", flush=True)

    # ---- 汇总 ----
    rows = []
    long_rows = []
    for lbl, accs in groups:
        for a in accs:
            r = results.get(a)
            if r is None:
                rows.append({"group": lbl, "acc": a, "source": "", "status": "MISSING",
                             "analysed": 0, "identified": 0, "n_pv_families": 0,
                             "pv_spots": 0.0, "viroid_spots": 0.0, "pv_spots_no_viroid": 0.0,
                             "pv_family_resolved": "no",
                             "top_pv_family": "", "top_pv_family_spots": 0.0,
                             "top_viroid_family": "", "top_viroid_spots": 0.0})
                continue
            fam = r["fam"]
            # GSA 只有 Viruses 根节点可用（树不细到科），单列，不混进科合计
            all_v = fam.pop(GSA_ALL_VIRUSES, 0.0) if r["source"] == "GSA" else 0.0
            viroid = {k: v for k, v in fam.items() if k in VIROID_FAMILIES}
            noviroid = {k: v for k, v in fam.items() if k not in VIROID_FAMILIES}
            pv_total = sum(fam.values()) if fam else all_v
            v_total = sum(viroid.values())
            # GSA 未落科时，Viruses 根节点值含类病毒，无法拆分 -> 扣除类病毒口径记 NA
            unresolved = (r["source"] == "GSA" and not fam)
            novir = "" if unresolved else pv_total - v_total
            top = max(fam.items(), key=lambda x: x[1]) if fam else (
                ("Viruses(仅根节点)" if all_v > 0 else ""), all_v)
            top_v = max(viroid.items(), key=lambda x: x[1]) if viroid else ("", 0.0)
            rows.append({
                "group": lbl, "acc": a, "source": r["source"], "status": r["status"],
                "analysed": r["meta"].get("analysed", 0),
                "identified": r["meta"].get("identified", 0),
                "n_pv_families": len(fam),
                "pv_spots": pv_total,
                "viroid_spots": v_total,
                "pv_spots_no_viroid": novir,
                "pv_family_resolved": "no" if unresolved else "yes",
                "top_pv_family": top[0], "top_pv_family_spots": top[1],
                "top_viroid_family": top_v[0], "top_viroid_spots": top_v[1],
            })
            for k, v in sorted(fam.items(), key=lambda x: -x[1]):
                long_rows.append({"group": lbl, "acc": a, "source": r["source"],
                                  "family": k, "spots": v,
                                  "is_viroid": "yes" if k in VIROID_FAMILIES else "no"})

    # ---- 写表 ----
    cols = ["group", "acc", "source", "status", "analysed", "identified",
            "n_pv_families", "pv_spots", "pv_spots_no_viroid", "viroid_spots",
            "pv_family_resolved", "top_pv_family", "top_pv_family_spots",
            "top_viroid_family", "top_viroid_spots"]
    per_path = os.path.join(args.outdir, "per_sample.tsv")
    with open(per_path, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in sorted(rows, key=lambda x: (x["group"], -float(x["pv_spots"] or 0))):
            w.writerow(r)

    long_path = os.path.join(args.outdir, "per_sample_family.tsv")
    with open(long_path, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["group", "acc", "source", "family",
                                          "spots", "is_viroid"],
                           delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in sorted(long_rows, key=lambda x: (x["group"], -x["spots"])):
            w.writerow(r)

    # ---- 报告 ----
    rep = os.path.join(args.outdir, "DISTRIBUTION_REPORT.md")
    with open(rep, "w", encoding="utf-8") as f:
        f.write("# 植物病毒科 spot 分布（每 run）\n\n")
        f.write(f"- 生成时间：{time.strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"- 植物病毒科名单：{len(fam_names)} 个（含类病毒 Pospiviroidae / Avsunviroidae）\n")
        f.write(f"- spot 口径：STAT `self_count`（直接命中，不含后代），按最近植物病毒科祖先归组\n")
        f.write("- ⚠ GSA（CRR）的 `hits` 是含后代的累计 reads 量级，与 STAT 的 k-mer spot **不可直接比较**\n\n")

        for lbl, accs in groups:
            sub = [r for r in rows if r["group"] == lbl]
            ok = [r for r in sub if r["status"] == "OK"]
            f.write(f"## {lbl}\n\n")
            f.write(f"- 样本数 {len(sub)}；有分类数据 {len(ok)}；"
                    f"状态 {dict(sorted(_count(r['status'] for r in sub).items()))}\n")
            if not ok:
                f.write("\n（无可用数据）\n\n")
                continue

            pv = [float(r["pv_spots"]) for r in ok]
            nv = [float(r["pv_spots_no_viroid"]) for r in ok
                  if str(r["pv_spots_no_viroid"]) != ""]
            vv = [float(r["viroid_spots"]) for r in ok]
            f.write(f"\n**全部植物病毒科合计**：{fmt_q(quantiles(pv))}\n\n")
            f.write(f"**扣除类病毒后**（n={len(nv)}，GSA 未落科的 "
                    f"{len(ok)-len(nv)} 条不计入）：{fmt_q(quantiles(nv))}\n\n")
            f.write(f"**类病毒（两科）**：{fmt_q(quantiles(vv))}\n\n")

            # 按数据来源分层 —— STAT 的 spot 与 GSA 的 hits 不是同一量级
            by_src: dict[str, list[dict]] = defaultdict(list)
            for r in ok:
                by_src[r["source"]].append(r)
            if len(by_src) > 1 or "GSA" in by_src:
                f.write("### 按数据来源分层（⚠ 两种口径不可合并）\n\n")
                f.write("| 来源 | n | 中位数 | p75 | p90 | p95 | p99 | max |\n")
                f.write("|---|---:|---:|---:|---:|---:|---:|---:|\n")
                for src in ("NCBI_STAT", "GSA"):
                    sub_s = by_src.get(src) or []
                    if not sub_s:
                        continue
                    q = quantiles([float(r["pv_spots"]) for r in sub_s])
                    f.write(f"| {src} | {q['n']} | {q['median']:,.0f} | {q['p75']:,.0f} | "
                            f"{q['p90']:,.0f} | {q['p95']:,.0f} | {q['p99']:,.0f} | "
                            f"{q['max']:,.0f} |\n")
                f.write("\n")
                gsa_rows = by_src.get("GSA") or []
                if gsa_rows:
                    unres = sum(1 for r in gsa_rows if r["pv_family_resolved"] == "no")
                    f.write(f"- GSA 通路 {len(gsa_rows)} 条中 **{unres}** 条无法落到科级"
                            f"（GSA 分类树普遍只细到属，无科节点），其 spot 值为 Viruses "
                            f"根节点累计 hits。\n")
                    f.write("- STAT 的 spot 是 k-mer 计数；GSA 的 hits 是 reads 计数，"
                            "**量级差约 1-2 个数量级，不可直接比较**。\n\n")

            f.write("### 分档（全部植物病毒科合计）\n\n")
            f.write("| spot 区间 | run 数 | 占比 |\n|---|---:|---:|\n")
            for lab, c in bin_counts(pv):
                f.write(f"| {lab} | {c} | {100.0*c/len(ok):.1f}% |\n")
            f.write("\n### 分档（扣除类病毒）\n\n")
            f.write("| spot 区间 | run 数 | 占比 |\n|---|---:|---:|\n")
            for lab, c in bin_counts(nv):
                f.write(f"| {lab} | {c} | {100.0*c/len(ok):.1f}% |\n")
            f.write("\n")

            # 科层面：出现频次 + spot 合计
            fam_freq: dict[str, int] = defaultdict(int)
            fam_sum: dict[str, float] = defaultdict(float)
            for r in long_rows:
                if r["group"] != lbl:
                    continue
                fam_freq[r["family"]] += 1
                fam_sum[r["family"]] += float(r["spots"])
            f.write("### 检出科 Top 20（按出现 run 数）\n\n")
            f.write("| 科 | run 数 | run 占比 | spot 合计 |\n|---|---:|---:|---:|\n")
            for k in sorted(fam_freq, key=lambda x: -fam_freq[x])[:20]:
                f.write(f"| {k} | {fam_freq[k]} | {100.0*fam_freq[k]/len(ok):.1f}% | "
                        f"{fam_sum[k]:,.0f} |\n")
            f.write("\n")

            # 类病毒 top 样本
            vt = sorted(ok, key=lambda r: -float(r["viroid_spots"]))[:10]
            if any(float(r["viroid_spots"]) > 0 for r in vt):
                f.write("### 类病毒 spot Top 10\n\n")
                f.write("| run | 类病毒 spot | 主要类病毒科 | 全部植物病毒 spot |\n")
                f.write("|---|---:|---|---:|\n")
                for r in vt:
                    if float(r["viroid_spots"]) <= 0:
                        continue
                    f.write(f"| {r['acc']} | {float(r['viroid_spots']):,.0f} | "
                            f"{r['top_viroid_family']} | {float(r['pv_spots']):,.0f} |\n")
                f.write("\n")

            # 最高 spot 样本
            f.write("### 植物病毒 spot Top 10\n\n")
            f.write("| run | 来源 | 植物病毒 spot | 扣除类病毒 | 主要科 | analysed |\n")
            f.write("|---|---|---:|---:|---|---:|\n")
            for r in sorted(ok, key=lambda x: -float(x["pv_spots"]))[:10]:
                nvs = (f"{float(r['pv_spots_no_viroid']):,.0f}"
                       if str(r["pv_spots_no_viroid"]) != "" else "NA")
                f.write(f"| {r['acc']} | {r['source']} | {float(r['pv_spots']):,.0f} | "
                        f"{nvs} | {r['top_pv_family']} | "
                        f"{int(r['analysed']):,} |\n")
            f.write("\n")

    print(f"\n[done] {per_path}\n[done] {long_path}\n[done] {rep}")
    return 0


def _count(it):
    d = defaultdict(int)
    for x in it:
        d[x] += 1
    return d


if __name__ == "__main__":
    sys.exit(main())
