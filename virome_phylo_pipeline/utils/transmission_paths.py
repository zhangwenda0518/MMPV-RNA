#!/usr/bin/env python3
"""
utils/transmission_paths.py — 事件级传播路径表 + episode 折叠 + 三级置信分类 + LTL
================================================================================
把 BEAST1 离散系统地理学产物（MCC 树 + tip 元数据）折叠成"事件级"解释层：

  1. branch_table  — 逐分支时空表（seraphim treeExtractions 同构字段）
  2. pathways      — 只保留"地点发生变化"的分支 = 一次传播事件（起止时间/地点/置信度）
  3. episodes      — 同地点连续节点折叠成"驻留时段"（回答：单次长驻留还是多次反复引入）
  4. classify      — 三级置信分类 Direct / Indirect / Distant-Import（框架借自 phymapr，
                     判据重写：SNP 维度用"期望替换数"= rate × Δheight × L，不用 mutations 注解）
  5. ltl           — 本地传播单元(LTL)持久性摘要（借自 polio-wpv1-phylodynamics，MIT）

方法学出处（逐条可追溯）
----------------------
* 分支级时空表 schema：
  git-repo/wnv_north_america/Scripts_&_data/Landscape_phylogeographic_analyses/
  R_script_analyses.r §1 (seraphim::treeExtractions; 另见 mccTreeExtraction.r)。
  字段对齐 startLon/startLat/endLon/endLat/startYear/endYear 语义。
* episode 折叠：phymapr R/epidemiologic_inference.R:103-146 的**思想**，
  但其实现按日期升序遍历（同日期样本会先于父节点处理 → episode_id 永久 NA 的幽灵组），
  本模块改为**自根向下 DFS（拓扑序）**，父节点必然先于子节点处理，同日期 bug 结构性不存在。
* 三级置信分类框架：phymapr R/epidemiologic_inference.R:237-250 的多维定级框架；
  判据**不照抄**——原版 SNP 维度读 `mutations` 树注解，经实测（解包 BEAST 1.10.4 jar）
  BEAST1/BEAST2/TreeTime 均不产出该注解，判据在一切合法输入下退化为时间单维。
  本模块 SNP 维度用期望替换数 E[S] = rate × Δheight × L（MCC 自带 rate/height 注解）。
  阈值默认沿用 phymapr 的 2/5 SNP，**必须在 METHODS.md 标定后才能进论文**。
* LTL/集群持久性：git-repo/polio-wpv1-phylodynamics/R/ltl.R:18-71 (MIT License) 的
  时间删失窗思想（默认 6 个月），简化为 episode 级持久性汇总。

自包含性：只依赖 numpy/pandas/标准库 + Bio.Phylo（管线既有依赖）；
注释解析复用 utils.tempmig_full.parse_beast_annotation（项目"单一判据只实现一次"）。
"""

from __future__ import annotations

import json
import os
import re
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

try:
    from utils.virphy_bridge import LogCollector
except ImportError:  # 允许 `python utils/transmission_paths.py` 直接运行
    from .virphy_bridge import LogCollector  # type: ignore


# ─────────────────────────────────────────────────────────────────────
# 树解析
# ─────────────────────────────────────────────────────────────────────

def _read_tree(tree_file: str):
    """Bio.Phylo 读树，NEXUS/Newick 自动探测（与 run_tempmig 同款探测逻辑）。"""
    from Bio import Phylo
    with open(tree_file, encoding="utf-8", errors="replace") as f:
        head = f.read(1024)
    fmt = "nexus" if head.lstrip().upper().startswith("#NEXUS") else "newick"
    return Phylo.read(tree_file, fmt)


_RATE_RE = re.compile(r"(?:^|[,&\s])rate=([0-9.eE+-]+)")
_POST_RE = re.compile(r"(?:^|[,&\s])posterior=([0-9.eE+-]+)")
_LOC_POINT_RE = re.compile(r'(?:^|[,&\s])Location="([^"]+)"')
_LOC_POINT_PROB_RE = re.compile(r"(?:^|[,&\s])Location\.prob=([0-9.eE+-]+)")


def _node_annotation(clade) -> Dict:
    """单节点注解 → location(dict state→prob) + height + rate + posterior。

    复用 parse_beast_annotation（覆盖 Location={X=p} / X.set+X.set.prob / Location="X"
    三种格式）；rate/posterior 用带边界正则补齐（避免误吃 rateCat/rate_median/rate_range）。
    """
    from utils.tempmig_full import parse_beast_annotation
    info = parse_beast_annotation(clade)
    comment = getattr(clade, "comment", None) or ""
    m = _RATE_RE.search(comment)
    if m:
        try:
            info["rate"] = float(m.group(1))
        except ValueError:
            pass
    m = _POST_RE.search(comment)
    if m:
        try:
            info["posterior"] = float(m.group(1))
        except ValueError:
            pass
    return info


def _load_metadata_locations(metadata_csv: str, name_col: str = "name",
                             location_col: str = "location") -> Dict[str, str]:
    df = pd.read_csv(metadata_csv, encoding="utf-8-sig")
    out = {}
    for _, row in df.iterrows():
        name = str(row.get(name_col, "") or "").strip()
        loc = str(row.get(location_col, "") or "").strip()
        if name and loc:
            out[name] = loc
    return out


def _match_known_location(meta_loc: str, known: List[str]) -> Optional[str]:
    """元数据长地址（如 'China, Ningxia, Yinchuan'）→ 树内地点词表（省名）。
    规则：已知地点名作为子串出现即命中；多命中取最长名（更具体）。"""
    hits = [k for k in known if k and k in meta_loc]
    if not hits:
        return None
    return max(hits, key=len)


def extract_tree_events(
    tree_file: str,
    metadata_csv: Optional[str] = None,
    name_col: str = "name",
    location_col: str = "location",
    log: Optional[LogCollector] = None,
) -> Dict:
    """解析 MCC 树 → 节点表(nodes) + 边表(edges)，含地点/高度/速率/后验。

    年龄约定：BEAST/TreeAnnotator 的 `height` 注解 = 距今年代（root 最大、tip=0），
    单位与定年一致（年）。缺 height 注解时退化为自根累加 branch length，
    并在 QC 里把 age_source 标为 'branch_length'（单位存疑，仅拓扑用途）。

    Returns
    -------
    dict: {nodes: List[dict], edges: List[dict], known_locations: List[str], qc: dict}
    """
    log = log or LogCollector()
    tree = _read_tree(tree_file)

    meta_loc = _load_metadata_locations(metadata_csv, name_col, location_col) if metadata_csv else {}

    # ── 第一遍：收集全部注解状态词表（与 spread3_viz 同约定：内部 argmax 必须落在已知地点内）──
    known: List[str] = []
    seen = set()
    for clade in tree.find_clades():
        info = _node_annotation(clade)
        for st in (info.get("location") or {}):
            if st not in seen:
                seen.add(st)
                known.append(st)

    # ── 第二遍：DFS 自根向下（拓扑序），父先于子 ──
    nodes: List[Dict] = []
    edges: List[Dict] = []
    n_height_annot = 0
    n_rate_annot = 0

    def _nodename(clade, idx: int) -> str:
        return clade.name if clade.name else f"node{idx}"

    counter = {"i": 0}

    def walk(clade, parent: Optional[Dict]):
        nonlocal n_height_annot, n_rate_annot
        counter["i"] += 1
        idx = counter["i"]
        info = _node_annotation(clade)
        locdist = info.get("location") or {}
        if locdist:
            loc = max(locdist, key=locdist.get)
            loc_prob = float(locdist[loc])
            loc_set_size = len(locdist)
        else:
            loc, loc_prob, loc_set_size = None, None, 0
        height = info.get("height")
        if height is not None:
            n_height_annot += 1
        rate = info.get("rate")
        if rate is not None:
            n_rate_annot += 1
        posterior = info.get("posterior")
        if posterior is None and loc_prob is not None:
            posterior = loc_prob  # 与 spread3_viz 同款回填

        node = {
            "node_id": f"n{idx}",
            "name": clade.name or "",
            "is_tip": bool(clade.is_terminal()),
            "location": loc,
            "location_prob": (round(loc_prob, 4) if loc_prob is not None else None),
            "location_set_size": loc_set_size,
            "age": (float(height) if height is not None else None),
            "rate": rate,
            "posterior": (round(posterior, 4) if posterior is not None else None),
            "parent_id": (parent["node_id"] if parent else ""),
            "_clade": clade,
            "_branch_len": (clade.branch_length or 0.0),
        }
        # tip 地点优先元数据（地面真值），注解兜底；内部节点反之
        if node["is_tip"] and clade.name and clade.name in meta_loc:
            mloc = _match_known_location(meta_loc[clade.name], known)
            if mloc:
                node["location"] = mloc
                node["location_source"] = "metadata"
            elif node["location"]:
                node["location_source"] = "annotation"
        elif node["location"]:
            node["location_source"] = "annotation"
        else:
            node["location_source"] = ""
        nodes.append(node)

        if parent is not None:
            edges.append({
                "branch_id": f"b{len(edges) + 1}",
                "parent_id": parent["node_id"],
                "child_id": node["node_id"],
            })
        for child in clade.clades:
            walk(child, node)

    root = tree.root
    walk(root, None)

    by_id = {n["node_id"]: n for n in nodes}

    # ── 年龄兜底：无 height 注解 → 自根累加 branch length（单位=subs，仅拓扑语义）──
    age_source = "height_annotation"
    n_missing_height = sum(1 for n in nodes if n["age"] is None)
    if n_missing_height == len(nodes):
        age_source = "branch_length_subs"
        acc: Dict[str, float] = {}

        def accum(node: Dict, base: float):
            base += node["_branch_len"]
            acc[node["node_id"]] = base
            for e in edges:
                if e["parent_id"] == node["node_id"]:
                    accum(by_id[e["child_id"]], base)

        accum(nodes[0], 0.0)
        for n in nodes:
            n["age"] = acc[n["node_id"]]
    elif n_missing_height:
        # 个别缺 height（TreeAnnotator 正常不会缺）→ parent_age - branch_len 近似
        for n in nodes:
            if n["age"] is None and n["parent_id"]:
                n["age"] = (by_id[n["parent_id"]]["age"] or 0.0) - n["_branch_len"]

    # ── 边表补全：起止年龄/地点/期望替换数 ──
    rates_present = [n["rate"] for n in nodes if n["rate"] is not None]
    median_rate = float(pd.Series(rates_present).median()) if rates_present else None

    for e in edges:
        p, c = by_id[e["parent_id"]], by_id[e["child_id"]]
        e.update({
            "parent_name": p["name"],
            "child_name": p["name"] if False else c["name"],
            "child_is_tip": c["is_tip"],
            "start_age": p["age"],
            "end_age": c["age"],
            "duration": (p["age"] - c["age"]) if (p["age"] is not None and c["age"] is not None) else None,
            "start_location": p["location"],
            "end_location": c["location"],
            "start_location_prob": p["location_prob"],
            "end_location_prob": c["location_prob"],
            "rate": (0.5 * (p["rate"] + c["rate"])
                     if (p["rate"] is not None and c["rate"] is not None)
                     else (p["rate"] if p["rate"] is not None else c["rate"])),
            "posterior": c["posterior"],
        })
        e["location_change"] = bool(
            e["start_location"] and e["end_location"]
            and e["start_location"] != e["end_location"]
        )
        e["rate_imputed"] = e["rate"] is None and median_rate is not None
        if e["rate"] is None and median_rate is not None:
            e["rate"] = median_rate

    qc = {
        "n_nodes": len(nodes),
        "n_tips": sum(1 for n in nodes if n["is_tip"]),
        "n_edges": len(edges),
        "known_locations": known,
        "age_source": age_source,
        "n_height_annot": n_height_annot,
        "n_rate_annot": n_rate_annot,
        "median_rate": median_rate,
        "tip_loc_coverage": round(sum(1 for n in nodes if n["is_tip"] and n["location"])
                                  / max(1, sum(1 for n in nodes if n["is_tip"])), 4),
        "internal_loc_coverage": round(sum(1 for n in nodes if not n["is_tip"] and n["location"])
                                       / max(1, sum(1 for n in nodes if not n["is_tip"])), 4),
        "n_rate_imputed": sum(1 for e in edges if e["rate_imputed"]),
    }
    # 解除对 Bio.Phylo 对象的引用，保证 nodes/edges 可安全 JSON/CSV 化
    for n in nodes:
        n.pop("_clade", None)
        n.pop("_branch_len", None)
    return {"nodes": nodes, "edges": edges, "known_locations": known, "qc": qc}


# ─────────────────────────────────────────────────────────────────────
# episode 折叠（拓扑序，父先于子）
# ─────────────────────────────────────────────────────────────────────

def fold_episodes(nodes: List[Dict], edges: List[Dict]) -> List[Dict]:
    """同地点连续节点折叠为驻留时段(episode)。

    遍历方式：nodes 已由 extract_tree_events 按自根 DFS 生成（拓扑序），
    逐节点继承父节点 episode 或开新 episode —— **不按日期排序**，
    phymapr 同日期 NA 幽灵组 bug（epidemiologic_inference.R:107-126）结构性规避。
    """
    by_id = {n["node_id"]: n for n in nodes}
    child_edges = {e["child_id"]: e for e in edges}
    episode_of: Dict[str, int] = {}
    episodes: List[Dict] = []

    for n in nodes:  # 拓扑序：父节点必在子节点之前
        pid = n["parent_id"]
        loc = n["location"] or "Unknown"
        p_ep = episode_of.get(pid) if pid else None
        if p_ep is not None:
            parent = by_id.get(pid, {})
            if (parent.get("location") or "Unknown") == loc:
                episode_of[n["node_id"]] = p_ep
                ep = episodes[p_ep]
                ep["n_nodes"] += 1
                # age：大=旧，小=新；episode 跨度取 [max(start_age), min(end_age)]
                if n["age"] is not None:
                    ep["start_age"] = max(ep["start_age"], n["age"])
                    ep["end_age"] = min(ep["end_age"], n["age"])
                continue
        episodes.append({
            "episode_id": len(episodes) + 1,
            "location": loc,
            "start_age": (n["age"] if n["age"] is not None else float("nan")),
            "end_age": (n["age"] if n["age"] is not None else float("nan")),
            "n_nodes": 1,
        })
        episode_of[n["node_id"]] = len(episodes) - 1

    for ep in episodes:
        ep["duration"] = ep["start_age"] - ep["end_age"]
    return episodes


# ─────────────────────────────────────────────────────────────────────
# 三级置信分类
# ─────────────────────────────────────────────────────────────────────

def classify_pathways(
    pathways: List[Dict],
    seq_len: Optional[float] = None,
    direct_snp: float = 2.0,
    indirect_snp: float = 5.0,
    direct_time: Optional[float] = None,
    indirect_time: Optional[float] = None,
    calibrate: Optional[str] = None,
) -> Tuple[List[Dict], Dict]:
    """三级分类 Direct / Indirect / Distant-Import（框架 phymapr，判据重写）。

    SNP 维度：期望替换数 E[S] = rate × duration × seq_len（期望值非实测 SNP，
    进论文须在 METHODS.md 声明）。seq_len 缺失 → 该维度不参与，tier 记 'unclassified'。
    时间维度：direct_time/indirect_time 显式给出才参与（默认不启用——
    phymapr 的 0.12/0.35 年魔法数无出处，borrow 报告 §3.3/§4.3 已论证）。

    calibrate:
      None / "manual" — 使用显式 direct_snp/indirect_snp（默认 2/5，phymapr 原值）。
      "auto"          — 数据驱动标定：取本数据 E[S] 分布的 q25/q75 作分界。
                        ⚠ 2/5 这类固定数没有跨尺度通用性（E[S] ∝ L×rate×时长；
                        实测同一阈值在 GCVA 16kb 上"全部 Distant"、放到 359bp
                        类病毒尺度上会"几乎全 Direct"）→ auto 让分界线跟随数据尺度。
                        事件数 <4 或 q75<=q25（分布过窄）时自动回退手动阈值，
                        回退原因写入 qc。
    """
    # 第一遍: 统一算 E[S]
    for p in pathways:
        p["expected_subs"] = None
        if seq_len and p.get("rate") is not None and p.get("duration") is not None:
            p["expected_subs"] = round(p["rate"] * p["duration"] * seq_len, 2)

    eff_direct, eff_indirect = float(direct_snp), float(indirect_snp)
    mode = "manual"
    calib_note = ""
    if calibrate == "auto":
        vals = [p["expected_subs"] for p in pathways if p.get("expected_subs") is not None]
        if len(vals) < 4:
            mode = "manual(fallback)"
            calib_note = f"auto 需 >=4 条含 E[S] 的事件, 实际 {len(vals)} -> 回退手动阈值"
        else:
            q25 = float(np.quantile(vals, 0.25))
            q75 = float(np.quantile(vals, 0.75))
            if q75 <= q25:
                mode = "manual(fallback)"
                calib_note = "E[S] 分布过窄 (q75<=q25) -> 回退手动阈值"
            else:
                eff_direct, eff_indirect = q25, q75
                mode = "auto(q25/q75)"
                calib_note = (f"E[S] 分位数标定: direct<=q25={q25:.4g}, "
                              f"indirect<=q75={q75:.4g}")

    used = {"snp_dim": seq_len is not None, "time_dim": direct_time is not None or indirect_time is not None}
    out = []
    for p in pathways:
        es = p.get("expected_subs")
        tier = "unclassified"
        if es is not None:
            snp_ok_direct = es <= eff_direct
            snp_ok_indirect = es <= eff_indirect
            t_ok_direct = (duration_ok(p["duration"], direct_time) if direct_time is not None else True)
            t_ok_indirect = (duration_ok(p["duration"], indirect_time) if indirect_time is not None else True)
            if snp_ok_direct and t_ok_direct:
                tier = "direct"
            elif snp_ok_indirect and t_ok_indirect:
                tier = "indirect"
            else:
                tier = "distant_import"
        p["confidence"] = tier
        out.append(p)
    qc = {
        "thresholds": {"direct": eff_direct, "indirect": eff_indirect,
                       "direct_time": direct_time, "indirect_time": indirect_time,
                       "requested_direct_snp": direct_snp,
                       "requested_indirect_snp": indirect_snp},
        "calibration": {"mode": mode, "note": calib_note},
        "dims_used": used,
        "note": "SNP 维度=期望替换数(rate×duration×seq_len)，非实测；阈值无论手动或 auto，"
                "进论文都必须在 METHODS.md 写明取值与依据",
    }
    return out, qc


def duration_ok(duration: Optional[float], limit: float) -> bool:
    return duration is not None and duration <= limit


# ─────────────────────────────────────────────────────────────────────
# 双流一致性闸门（时间流 vs 地理流）
# ─────────────────────────────────────────────────────────────────────

_RATE_ALL_RE = re.compile(r"(?:^|[,&\s])rate=([0-9.eE+-]+)")


def timescale_consistency(
    mcc_tree: str,
    clock_summary_tsv: Optional[str],
    out_json: str,
    log: Optional[LogCollector] = None,
) -> Dict:
    """时间流 vs 地理流的时间尺度并排比对（双流闸门，2026-09-17 新增）。

    时间流: clock_summary.tsv 的 genome 级 RTT 速率 + DRT 判据（clock stage）。
    地理流: phylogeo MCC 的根高（中位/95%HPD）与节点速率中位（BEAST 联合后验）。

    背景（为什么需要）：两条流是**独立推断**，各自内部自洽，但估计量不同
    （ML 树上的 RTT 回归 vs BEAST 联合后验），速率/根高可系统偏离——
    实测 GCVA 全量：时间流 9.3e-05 vs 地理流 3.3e-04（3.5 倍），且时间流
    DRT 阴性。闸门把两个口径**摆在同一张纸上**并打旗标，防止论文混引。

    flag:
      OK   — 速率比在 [0.5, 2] 且时间流 DRT 通过
      WARN — 速率比超出 2 倍带 **或** 时间流 DRT 阴性（地理流绝对年代须谨慎引用）
    """
    log = log or LogCollector()
    # ── 地理流 ──
    text = open(mcc_tree, encoding="utf-8", errors="replace").read()
    # TreeAnnotator 产出 height_median= 与 height= 同值并存; 手工/合成树可能只有 height=
    h_meds = [float(x) for x in re.findall(r"height(?:_median)?=([\d.eE+-]+)", text)]
    root_h = max(h_meds) if h_meds else None
    root_hpd = None
    for a, b in re.findall(r"height_95%_HPD=\{([\d.eE+-]+),([\d.eE+-]+)\}", text):
        a, b = float(a), float(b)
        if root_hpd is None or b > root_hpd[1]:   # 上界最高的 HPD = 根（最老节点）
            root_hpd = [a, b]
    rates = [float(x) for x in _RATE_ALL_RE.findall(text)]
    geo_rate = (sorted(rates)[len(rates) // 2] if rates else None)

    # ── 时间流 ──
    time_rate, drt_passed = None, None
    if clock_summary_tsv and os.path.exists(clock_summary_tsv):
        import csv as _csv
        with open(clock_summary_tsv, encoding="utf-8-sig") as f:
            for row in _csv.DictReader(f, delimiter="\t"):
                if str(row.get("level", "")).strip().lower() == "genome":
                    try:
                        time_rate = float(row.get("beta"))
                    except (TypeError, ValueError):
                        time_rate = None
                    dp = str(row.get("drt_passed", "")).strip().lower()
                    drt_passed = {"true": True, "false": False}.get(dp)
                    break

    ratio = (geo_rate / time_rate) if (geo_rate and time_rate) else None
    flags = []
    if ratio is not None and (ratio > 2 or ratio < 0.5):
        flags.append(f"速率比 {ratio:.2f}x 超出 2 倍带 — 两流口径必须分开引用")
    if drt_passed is False:
        flags.append("时间流 DRT 阴性 — 绝对 TMRCA 不报；地理流绝对年代引用须附 HPD 并谨慎")
    if not flags:
        flags.append("速率比在 2 倍带内且 DRT 通过")
    flag = "WARN" if any(f.startswith(("速率比", "时间流")) for f in flags) else "OK"

    out = {
        "time_stream": {"rate_rtt": time_rate, "drt_passed": drt_passed,
                        "source": os.path.basename(clock_summary_tsv) if clock_summary_tsv else None},
        "geo_stream": {"root_height_median_yr": root_h,
                       "root_height_95hpd_yr": root_hpd,
                       "median_rate": geo_rate,
                       "source": os.path.basename(mcc_tree)},
        "rate_ratio_geo_vs_time": (round(ratio, 3) if ratio else None),
        "flag": flag,
        "flags": flags,
    }
    os.makedirs(os.path.dirname(os.path.abspath(out_json)), exist_ok=True)
    with open(out_json, "w", encoding="utf-8") as f:
        json.dump(out, f, ensure_ascii=False, indent=2)
    for fl in flags:
        log.emit(f"  [timescale] {flag}: {fl}")
    return out


# ─────────────────────────────────────────────────────────────────────
# LTL 本地传播单元持久性（polio-wpv1 ltl.R 思想，episode 级简化）
# ─────────────────────────────────────────────────────────────────────

def build_ltl(episodes: List[Dict], censor_years: float = 0.5) -> pd.DataFrame:
    """按地点汇总 episode → LTL 持久性表。

    persistent 判据：该地点存在时长 ≥ censor_years（默认 0.5 年，polio 默认删失窗）
    的驻留时段。Unknown 地点排除。
    """
    df = pd.DataFrame([e for e in episodes if e["location"] != "Unknown"])
    if df.empty:
        return pd.DataFrame(columns=["location", "n_episodes", "total_persistent_years",
                                     "max_episode_years", "first_age", "last_age", "persistence"])
    g = df.groupby("location")
    rows = []
    for loc, sub in g:
        durs = sub["duration"].dropna()
        rows.append({
            "location": loc,
            "n_episodes": int(len(sub)),
            "total_persistent_years": round(float(durs.sum()), 4) if len(durs) else 0.0,
            "max_episode_years": round(float(durs.max()), 4) if len(durs) else 0.0,
            "first_age": round(float(sub["start_age"].max()), 4),   # 最大 age = 最早出现
            "last_age": round(float(sub["end_age"].min()), 4),      # 最小 age = 最晚出现
            "persistence": ("persistent" if (len(durs) and durs.max() >= censor_years)
                            else "transient"),
        })
    return pd.DataFrame(rows).sort_values(["persistence", "n_episodes"],
                                          ascending=[False, False])


# ─────────────────────────────────────────────────────────────────────
# 编排入口
# ─────────────────────────────────────────────────────────────────────

def run_pathways(
    mcc_tree: str,
    metadata_csv: Optional[str] = None,
    output_dir: str = "pathways",
    fasta: Optional[str] = None,
    seq_len: Optional[float] = None,
    name_col: str = "name",
    location_col: str = "location",
    direct_snp: float = 2.0,
    indirect_snp: float = 5.0,
    direct_time: Optional[float] = None,
    indirect_time: Optional[float] = None,
    censor_years: float = 0.5,
    calibrate: Optional[str] = None,
    log: Optional[LogCollector] = None,
) -> Dict:
    """geo_stage 旁挂子分析：MCC 树 → 事件级传播路径全套产物。

    calibrate="auto": 用本数据 E[S] 分布 q25/q75 自动标定 Direct/Indirect 分界
    （事件数<4 或分布过窄自动回退手动值；选中的阈值与原因写入 pathway_qc.json）。

    输出（output_dir/）：
      branch_table.csv   逐分支时空表（seraphim 同构，供 diffusion_stats/地图复用）
      pathways.csv       事件级传播表（仅地点变化的分支，含三级置信）
      episodes.csv       驻留时段表
      ltl_summary.csv    地点持久性摘要
      pathway_qc.json    覆盖率/来源/阈值 QC
    """
    log = log or LogCollector()
    os.makedirs(output_dir, exist_ok=True)
    res: Dict = {"success": False, "error": None}

    # seq_len：显式给出 > 从 fasta 推平均长度
    if seq_len is None and fasta:
        from Bio import SeqIO
        lens = [len(r.seq) for r in SeqIO.parse(fasta, "fasta")]
        if lens:
            seq_len = float(sum(lens) / len(lens))
            log.emit(f"seq_len 从 fasta 推得: {seq_len:.0f} bp (平均 {len(lens)} 条)")

    ev = extract_tree_events(mcc_tree, metadata_csv, name_col, location_col, log=log)
    nodes, edges, qc = ev["nodes"], ev["edges"], ev["qc"]

    # "Unknown" 是 BEAST 对无地点 tip 的真实祖先状态, 不是可映射地点 → 不进路径表
    # (2026-09-17 修: PSTVd 全量 MCC 实测含 Unknown 节点, 原实现让其进入 pathways
    #  导致下游地图/扩散统计因无法解析坐标而失败)
    n_unknown_excluded = sum(1 for e in edges if e["location_change"]
                             and "Unknown" in (e["start_location"], e["end_location"]))
    pathways = [dict(e) for e in edges if e["location_change"]
                and "Unknown" not in (e["start_location"], e["end_location"])]
    episodes = fold_episodes(nodes, edges)
    pathways, cls_qc = classify_pathways(
        pathways, seq_len=seq_len, direct_snp=direct_snp, indirect_snp=indirect_snp,
        direct_time=direct_time, indirect_time=indirect_time, calibrate=calibrate)
    ltl = build_ltl(episodes, censor_years)

    # 清理内部标记列后落盘
    def _clean(rows: List[Dict]) -> pd.DataFrame:
        df = pd.DataFrame(rows)
        return df.drop(columns=["rate_imputed"], errors="ignore")

    branch_csv = os.path.join(output_dir, "branch_table.csv")
    _clean(edges).to_csv(branch_csv, index=False)
    path_csv = os.path.join(output_dir, "pathways.csv")
    _clean(pathways).to_csv(path_csv, index=False)
    epi_csv = os.path.join(output_dir, "episodes.csv")
    pd.DataFrame(episodes).to_csv(epi_csv, index=False)
    ltl_csv = os.path.join(output_dir, "ltl_summary.csv")
    ltl.to_csv(ltl_csv, index=False)

    qc.update({
        "n_pathways": len(pathways),
        "n_unknown_excluded": n_unknown_excluded,
        "n_episodes": len(episodes),
        "seq_len": seq_len,
        "classification": cls_qc,
        "confidence_counts": pd.Series([p["confidence"] for p in pathways]).value_counts().to_dict()
        if pathways else {},
        "inputs": {"mcc_tree": os.path.abspath(mcc_tree),
                   "metadata_csv": os.path.abspath(metadata_csv) if metadata_csv else None,
                   "fasta": os.path.abspath(fasta) if fasta else None},
    })
    qc_path = os.path.join(output_dir, "pathway_qc.json")
    with open(qc_path, "w", encoding="utf-8") as f:
        json.dump(qc, f, ensure_ascii=False, indent=2)

    log.emit(f"pathways: {len(pathways)} 传播事件 / {len(episodes)} 驻留时段 / "
             f"tip 地点覆盖 {qc['tip_loc_coverage']:.0%} / 内部节点覆盖 {qc['internal_loc_coverage']:.0%}")
    for t, c in sorted(qc["confidence_counts"].items()):
        log.emit(f"  confidence={t}: {c}")

    res.update({"success": True, "branch_table": branch_csv, "pathways": path_csv,
                "episodes": epi_csv, "ltl_summary": ltl_csv, "qc": qc_path,
                # 2026-09-16 补: "qc" 键是**路径字符串**(易误用为 dict, stage 接线时踩过),
                # 这里同时给出解析好的 dict, 调用方一律用 qc_data
                "qc_data": qc,
                "n_pathways": len(pathways), "n_episodes": len(episodes),
                "nodes": nodes, "edges": edges})
    return res


# ─────────────────────────────────────────────────────────────────────
# CLI（独立运行）
# ─────────────────────────────────────────────────────────────────────

def main():
    import argparse
    ap = argparse.ArgumentParser(
        description="事件级传播路径表：MCC 树 → branch_table/pathways/episodes/LTL")
    ap.add_argument("--mcc-tree", required=True, help="BEAST MCC 树 (NEXUS/Newick)")
    ap.add_argument("--metadata", default=None, help="name,location CSV (tip 地点兜底)")
    ap.add_argument("--fasta", default=None, help="比对 fasta (推 seq_len 供期望替换数)")
    ap.add_argument("--seq-len", type=float, default=None, help="显式序列长度 (覆盖 fasta 推断)")
    ap.add_argument("-o", "--outdir", default="pathways")
    ap.add_argument("--direct-snp", type=float, default=2.0)
    ap.add_argument("--indirect-snp", type=float, default=5.0)
    ap.add_argument("--calibrate", default="manual", choices=["auto", "manual"],
                    help="auto=按本数据 E[S] 分布 q25/q75 自动标定置信阈值 (manual 时用上面的显式值)")
    ap.add_argument("--direct-time", type=float, default=None,
                    help="可选：direct 时间上限(年)；默认不启用时间维")
    ap.add_argument("--indirect-time", type=float, default=None)
    ap.add_argument("--censor-years", type=float, default=0.5, help="LTL 持久性删失窗(年)")
    args = ap.parse_args()
    log = LogCollector()
    r = run_pathways(args.mcc_tree, args.metadata, args.outdir, fasta=args.fasta,
                     seq_len=args.seq_len, direct_snp=args.direct_snp,
                     indirect_snp=args.indirect_snp, direct_time=args.direct_time,
                     indirect_time=args.indirect_time, censor_years=args.censor_years,
                     calibrate=args.calibrate,
                     log=log)
    if not r["success"]:
        raise SystemExit(f"pathways failed: {r['error']}")
    print("PATHWAYS_DONE")


if __name__ == "__main__":
    main()
