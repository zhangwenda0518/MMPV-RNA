#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
run_digest.py — MMPV 运行状态数字化 (四层递进架构 · L1 状态层)
================================================================

把一个项目目录里各管线的「跑没跑、跑到哪、关键计数是多少」抽成一份
机器可读的 run digest (JSON), 落 <项目>/90_Handoff/state/run_digest.json。

这是 AI 驱动四层架构的地基: L2 解读 (ai.summarize_run)、L3 顾问
(advisor)、L4 copilot/MCP 读的都是它 —— AI 只见摘要, 永不见原始
fastq/大 TSV, 成本可控且不出平台。

设计:
  * 纯标准库 + mmpv_common.io_layout (目录名以 io_layout 注册表为准,
    standard / legacy 双布局通吃; 目录缺失一律记 "absent" 而不报错 ——
    项目本来就是渐进跑出来的);
  * 解析失败的摘要文件记 "unparsable", 不让单个坏文件拖垮整个 digest;
  * schema_version 从 1 起, 字段只加不改。

CLI:
  python -m mmpv_common.run_digest <项目根> [--write] [--compact]
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Optional

try:
    from mmpv_common import io_layout
except ImportError:  # 直接以脚本运行时的兜底
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
    from mmpv_common import io_layout  # type: ignore

# ═══════════════════════════════════════════════════════════════════
# 常量
# ═══════════════════════════════════════════════════════════════════

SCHEMA_VERSION = 1

# standard 布局: 各管线输出根 (doc/IO_LAYOUT_DESIGN.md); legacy 布局目录
# 全部平铺在项目根, 由 io_layout.layout_dirs() 的映射表给出目录名。
PIPELINE_ROOTS: Dict[str, str] = {
    "public":          "01_PublicData",
    "preprocessing":   "02_Preprocessing",
    "discovery":       "03_Discovery",
    "analysis":        "04_Analysis",
    "phylo":           "05_Phylo",
    # ⑥ 已归档 2026-10-11 → archive/20261011_endogenous_virus_pipeline; 键保留以兼容历史输出
    "eve":             "06_EVE",
    # ⑧ 已归档 2026-10-11、2026-10-27 移出主仓至 `差异表达分析` 项目归档; 键保留以兼容历史输出
    "host_expression": "08_HostExpression",
}

HANDOFF_DIR = "90_Handoff"
STATE_SUBDIR = "state"
DIGEST_NAME = "run_digest.json"

# io_layout 逻辑键 → 所属管线 (digest 只统计这些键, 其余忽略)
KEYS_BY_PIPELINE: Dict[str, List[str]] = {
    "public":        ["m_search", "m_info", "m_prescreen", "m_down", "m_plot",
                      "h_genome", "h_hostdb"],
    "preprocessing": ["d_clean", "d_hostdep", "d_bbnorm"],
    "discovery":     ["d_asm", "d_ident", "d_filter", "d_graph", "d_cobra",
                      "d_merge", "d_cluster", "d_taxonomy", "d_host_pred",
                      "d_checkv", "d_rescue", "d_analysis", "d_verify",
                      "d_reports"],
    # [阶段重排 2026-11] a_align = 03_Extract (比对与制备); a_variants 与 a_post
    # 现共用 04_Variants —— 由 _digest_analysis 按真实路径去重, 体积只计一次
    "analysis":      ["a_detect", "a_filter", "a_align", "a_variants", "a_post",
                      "a_assembly", "a_extract", "a_similarity", "a_dvg",
                      "a_report"],
    "phylo":         ["ph_data", "ph_phylogeny", "ph_popgen", "ph_recomb",
                      "ph_select", "ph_time", "ph_geography", "ph_report"],
    "eve":           ["e_loci", "e_verdict", "e_rvdb", "e_summary"],
    "host_expression": ["he_hostref", "he_quant", "he_qc", "he_labels",
                        "he_coexp", "he_de", "he_enrich", "he_loadscan",
                        "he_report"],
}

READ_EXTS = {".fastq", ".gz", ".fq"}
SUMMARY_CAP = 8   # 文件名清单最多列几个, 防 digest 膨胀


# ═══════════════════════════════════════════════════════════════════
# 底层统计
# ═══════════════════════════════════════════════════════════════════

def _dir_stats(path: Path) -> Dict[str, Any]:
    """目录统计: 文件数/字节/子目录名。不存在 → {exists: False}。"""
    if not path.is_dir():
        return {"exists": False}
    files = 0
    total = 0
    subdirs: List[str] = []
    try:
        for entry in path.iterdir():
            if entry.is_dir():
                subdirs.append(entry.name)
            elif entry.is_file():
                files += 1
                try:
                    total += entry.stat().st_size
                except OSError:
                    pass
    except OSError as e:
        return {"exists": True, "error": str(e)}
    return {"exists": True, "files": files, "bytes": total,
            "subdirs": sorted(subdirs)[:SUMMARY_CAP],
            "n_subdirs": len(subdirs)}


def _count_rows(path: Path) -> Optional[int]:
    """数据行数 (减表头); 读不动返回 None。"""
    try:
        with open(path, "r", encoding="utf-8", errors="replace") as f:
            n = sum(1 for line in f if line.strip())
        return max(n - 1, 0)
    except OSError:
        return None


def _read_lines(path: Path, cap: int = 100000) -> Optional[int]:
    try:
        with open(path, "r", encoding="utf-8", errors="replace") as f:
            return sum(1 for _ in zip(range(cap), f))
    except OSError:
        return None


def _glob_first(root: Path, pattern: str) -> Optional[Path]:
    if not root.is_dir():
        return None
    try:
        for p in sorted(root.rglob(pattern)):
            return p
    except OSError:
        pass
    return None


def _list_html(root: Path) -> List[str]:
    if not root.is_dir():
        return []
    try:
        return sorted(p.name for p in root.rglob("*.html"))[:SUMMARY_CAP]
    except OSError:
        return []


# ═══════════════════════════════════════════════════════════════════
# 逐管线摘要 (每个都独立容错)
# ═══════════════════════════════════════════════════════════════════

def _count_fastq(path: Path) -> int:
    if not path.is_dir():
        return 0
    n = 0
    for p in path.rglob("*"):
        if p.is_file() and (p.suffix in READ_EXTS or p.name.endswith((".fastq.gz", ".fq.gz"))):
            n += 1
    return n


def _fastq_sizes(path: Path, cap: int = 200) -> List[int]:
    """reads 文件大小样本 (字节), 供深度差异规则用。"""
    sizes: List[int] = []
    if not path.is_dir():
        return sizes
    for p in sorted(path.rglob("*")):
        if p.is_file() and p.name.endswith((".fastq.gz", ".fq.gz")):
            try:
                sizes.append(p.stat().st_size)
            except OSError:
                pass
            if len(sizes) >= cap:
                break
    return sizes


def _digest_public(base: Path, dirs: Dict[str, str]) -> Dict[str, Any]:
    out: Dict[str, Any] = {"stages": {}}
    for key in KEYS_BY_PIPELINE["public"]:
        out["stages"][key] = _dir_stats(base / dirs[key])
    core = _glob_first(base / dirs["m_info"], "*Unified_Metadata_Core14*.csv") \
        if out["stages"]["m_info"].get("exists") else None
    out["core14_rows"] = _count_rows(core) if core else None
    out["core14_file"] = core.name if core else None
    raw = base / dirs["m_down"]
    out["raw_reads"] = _count_fastq(raw) if raw.is_dir() else 0
    out["raw_sra"] = sum(1 for p in raw.rglob("*.sra")) if raw.is_dir() else 0
    return out


def _digest_preprocessing(base: Path, dirs: Dict[str, str]) -> Dict[str, Any]:
    out: Dict[str, Any] = {"stages": {}}
    for key in KEYS_BY_PIPELINE["preprocessing"]:
        out["stages"][key] = _dir_stats(base / dirs[key])
    clean = base / dirs["d_clean"]
    dep = base / dirs["d_hostdep"]
    out["clean_reads"] = _count_fastq(clean) if clean.is_dir() else 0
    out["depleted_reads"] = _count_fastq(dep) if dep.is_dir() else 0
    out["clean_read_sizes"] = _fastq_sizes(clean)
    summary = _glob_first(base, "preprocessing_summary.tsv")
    out["summary_rows"] = _count_rows(summary) if summary else None
    ready = _glob_first(base, "assembly_ready.list")
    out["assembly_ready"] = _read_lines(ready) if ready else None
    return out


def _digest_discovery(base: Path, dirs: Dict[str, str]) -> Dict[str, Any]:
    out: Dict[str, Any] = {"stages": {}}
    for key in KEYS_BY_PIPELINE["discovery"]:
        out["stages"][key] = _dir_stats(base / dirs[key])
    asm = base / dirs["d_asm"]
    # 组装样本数 = 组装目录下的样本子目录数
    out["assembly_samples"] = (
        len([d for d in asm.iterdir() if d.is_dir()]) if asm.is_dir() else 0)
    ident = out["stages"]["d_ident"]
    out["identification_tools"] = ident.get("n_subdirs", 0) if ident.get("exists") else 0
    rescue = base / dirs["d_rescue"]
    out["rescue_report_html"] = (
        _list_html(base / dirs["d_reports"]) if rescue.is_dir() else [])
    out["report_html"] = _list_html(base / dirs["d_reports"])
    # HQ vOTU 粗计: cluster 目录下的 fasta 计数交由 deeper parse, 这里只看目录在不在
    out["cluster_present"] = bool(out["stages"]["d_cluster"].get("exists"))
    return out


def _digest_analysis(base: Path, dirs: Dict[str, str]) -> Dict[str, Any]:
    out: Dict[str, Any] = {"stages": {}}
    seen: Dict[str, str] = {}          # 真实路径 → 首个统计它的逻辑键
    for key in KEYS_BY_PIPELINE["analysis"]:
        path = base / dirs[key]
        try:
            sig = str(path.resolve())
        except OSError:
            sig = str(path)
        if sig in seen:
            # a_variants 与 a_post 共用 04_Variants: 体积/文件数不得重复计入
            out["stages"][key] = {"exists": out["stages"][seen[sig]].get("exists", False),
                                  "aliased_to": seen[sig]}
            continue
        seen[sig] = key
        out["stages"][key] = _dir_stats(path)
    det = _glob_first(base / dirs["a_detect"], "*summary*.tsv")
    out["detection_summary_rows"] = _count_rows(det) if det else None
    out["report_html"] = _list_html(base / dirs["a_report"])
    return out


def _digest_phylo(base: Path, dirs: Dict[str, str]) -> Dict[str, Any]:
    out: Dict[str, Any] = {"stages": {}}
    for key in KEYS_BY_PIPELINE["phylo"]:
        out["stages"][key] = _dir_stats(base / dirs[key])
    out["report_html"] = _list_html(base / dirs["ph_report"])
    return out


def _digest_eve(base: Path, dirs: Dict[str, str]) -> Dict[str, Any]:
    out: Dict[str, Any] = {"stages": {}}
    for key in KEYS_BY_PIPELINE["eve"]:
        out["stages"][key] = _dir_stats(base / dirs[key])
    loci = base / dirs["e_loci"]
    out["genomes"] = len([d for d in loci.iterdir() if d.is_dir()]) if loci.is_dir() else 0
    ks = _glob_first(base / dirs["e_summary"], "kingdom_summary.tsv")
    out["kingdom_summary_rows"] = _count_rows(ks) if ks else None
    return out


def _digest_hostexpr(base: Path, dirs: Dict[str, str]) -> Dict[str, Any]:
    out: Dict[str, Any] = {"stages": {}}
    for key in KEYS_BY_PIPELINE["host_expression"]:
        out["stages"][key] = _dir_stats(base / dirs[key])
    out["report_html"] = _list_html(base / dirs["he_report"])
    return out


_DIGESTERS = {
    "public": _digest_public,
    "preprocessing": _digest_preprocessing,
    "discovery": _digest_discovery,
    "analysis": _digest_analysis,
    "phylo": _digest_phylo,
    "eve": _digest_eve,
    "host_expression": _digest_hostexpr,
}


# ═══════════════════════════════════════════════════════════════════
# 主入口
# ═══════════════════════════════════════════════════════════════════

def state_dir(project_dir: Path) -> Path:
    return Path(project_dir) / HANDOFF_DIR / STATE_SUBDIR


def guess_state_dir(output_root: Path) -> Optional[Path]:
    """从管线输出根反推项目 90_Handoff/state。

    仅当父目录已有 90_Handoff/ 才返回 (保守: 不为单跑一个报告的项目凭空
    造目录); 否则 None (调用方跳过溯源日志)。"""
    proj = Path(output_root).resolve().parent
    if (proj / HANDOFF_DIR).is_dir():
        return proj / HANDOFF_DIR / STATE_SUBDIR
    return None


def find_state_root(start: Path, max_up: int = 3) -> Optional[Path]:
    """从任意起点向上找含 90_Handoff/ 的项目根 (最多 max_up 层)。

    供过程面板等"报告内嵌全项目信息"的场景: 报告文件可能在项目任意深度,
    找不到返回 None, 调用方优雅放弃 —— 绝不臆造项目根。"""
    p = Path(start).resolve()
    for _ in range(max_up + 1):
        if (p / HANDOFF_DIR).is_dir():
            return p
        if p.parent == p:
            break
        p = p.parent
    return None


def digest_path(project_dir: Path) -> Path:
    return state_dir(project_dir) / DIGEST_NAME


def build_digest(project_dir: Path, layout: Optional[str] = None) -> Dict[str, Any]:
    """扫描项目根, 产出 run digest。任何单点失败都不中断 (记 absent/error)。"""
    project_dir = Path(project_dir).resolve()
    lay = io_layout.resolve_layout(layout)
    dirs = io_layout.layout_dirs(lay)
    digest: Dict[str, Any] = {
        "schema_version": SCHEMA_VERSION,
        "generated_at": time.strftime("%Y-%m-%d %H:%M:%S"),
        "project_dir": str(project_dir),
        "layout": lay,
        "pipelines": {},
    }
    for pk in _DIGESTERS:
        if lay == "standard":
            base = project_dir / PIPELINE_ROOTS[pk]
        else:
            # legacy: 目录平铺在项目根 (键名映射不变)
            base = project_dir
        if lay == "standard" and not base.is_dir():
            digest["pipelines"][pk] = {"present": False}
            continue
        try:
            section = _DIGESTERS[pk](base, dirs)
            section["present"] = any(
                s.get("exists") for s in section.get("stages", {}).values())
            digest["pipelines"][pk] = section
        except Exception as e:  # 单管线解析失败不拖垮整体
            digest["pipelines"][pk] = {"present": False, "error": f"{type(e).__name__}: {e}"}
    present = [pk for pk, v in digest["pipelines"].items() if v.get("present")]
    digest["totals"] = {
        "pipelines_present": len(present),
        "pipelines_present_list": present,
        "project_exists": project_dir.is_dir(),
    }
    return digest


HISTORY_NAME = "digest_history.jsonl"


def write_digest(digest: Dict[str, Any], project_dir: Path) -> Path:
    out = digest_path(project_dir)
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", encoding="utf-8") as f:
        json.dump(digest, f, ensure_ascii=False, indent=2)
    # 过程时间线: 每次落盘追加一行快照 (过程报告的数据源)
    try:
        snap = {
            "ts": digest.get("generated_at") or time.strftime("%Y-%m-%d %H:%M:%S"),
            "pipelines_present": digest.get("totals", {}).get("pipelines_present", 0),
            "pipelines": {k: bool(v.get("present"))
                          for k, v in digest.get("pipelines", {}).items()},
        }
        with open(out.parent / HISTORY_NAME, "a", encoding="utf-8") as hf:
            hf.write(json.dumps(snap, ensure_ascii=False) + "\n")
    except OSError:
        pass
    return out


def load_history(project_dir: Path) -> List[Dict[str, Any]]:
    """digest 落盘历史 (时间线), 最早→最近。无历史返回空表。"""
    path = state_dir(project_dir) / HISTORY_NAME
    out: List[Dict[str, Any]] = []
    try:
        with open(path, "r", encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if line:
                    try:
                        out.append(json.loads(line))
                    except json.JSONDecodeError:
                        continue
    except OSError:
        pass
    return out


def load_digest(path: Path) -> Dict[str, Any]:
    with open(path, "r", encoding="utf-8") as f:
        d = json.load(f)
    if not isinstance(d, dict) or "schema_version" not in d:
        raise ValueError(f"不是合法的 run digest: {path}")
    return d


def compact_for_ai(digest: Dict[str, Any]) -> Dict[str, Any]:
    """给 AI 提示词用的瘦身版: 去子目录清单, 阶段只留 exists/files/bytes。"""
    slim = {k: v for k, v in digest.items() if k != "pipelines"}
    slim["pipelines"] = {}
    for pk, sec in digest.get("pipelines", {}).items():
        s: Dict[str, Any] = {"present": sec.get("present", False)}
        for k, v in sec.items():
            if k == "stages":
                continue
            if isinstance(v, list) and len(v) > 5:
                v = v[:5] + ["..."]
            s[k] = v
        stages = {}
        for sk, sv in sec.get("stages", {}).items():
            if not sv.get("exists"):
                stages[sk] = "absent"
            else:
                stages[sk] = {k2: v2 for k2, v2 in sv.items() if k2 in ("files", "bytes", "n_subdirs")}
        s["stages"] = stages
        slim["pipelines"][pk] = s
    return slim


# ═══════════════════════════════════════════════════════════════════

def main(argv: Optional[List[str]] = None) -> int:
    ap = argparse.ArgumentParser(
        description="MMPV run digest — 项目运行状态数字化 (L1)")
    ap.add_argument("project", help="项目根目录")
    ap.add_argument("--io-layout", choices=["standard", "legacy"], default=None)
    ap.add_argument("--write", action="store_true",
                    help=f"写入 {HANDOFF_DIR}/{STATE_SUBDIR}/{DIGEST_NAME}")
    ap.add_argument("--compact", action="store_true", help="输出 AI 瘦身版")
    args = ap.parse_args(argv)

    d = build_digest(Path(args.project), layout=args.io_layout)
    if args.compact:
        d = compact_for_ai(d)
    payload = json.dumps(d, ensure_ascii=False, indent=2)
    if args.write:
        out = write_digest(build_digest(Path(args.project), layout=args.io_layout),
                           Path(args.project))
        print(f"已写入: {out}", file=sys.stderr)
        print(payload)
    else:
        print(payload)
    return 0


if __name__ == "__main__":
    sys.exit(main())
