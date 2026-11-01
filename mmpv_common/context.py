#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
context.py — MMPV AI 上下文组装 (Phase 1 · 摘要层地基)
======================================================

把「一个阶段此刻长什么样」收敛为 :class:`StageContext`: run digest、阶段卡片、
``90_Handoff/ai/stage_state`` 记录、日志尾部、管线输出目录统计 —— 一份
JSON-serializable 的上下文, 供摘要 (summaries.py)、doctor、watch 复用。
设计纪律沿用 run_digest:

  * **只读**: 本模块不创建、不修改项目里的任何文件;
  * **不抛异常**: 目录 / JSON / 日志缺失或 OSError 一律降级为
    ``None`` / 空串 / ``{"exists": False}`` —— 摘要层永远有东西可写;
  * **布局感知**: standard / legacy 双布局经 io_layout 解析; stage 名到
    digest 目录键的映射见 :data:`STAGE_IO_KEYS` (近似映射已注明);
  * **JSON 安全**: :meth:`StageContext.to_dict` 保证 ``json.dumps`` 可序列化,
    日志超长时只保留最后 ``max_chars`` 字符并标记 ``truncated``。

术语约定 (与各编排器一致):
  ``output_dir``   管线输出根 (③ 即 ``<项目>/03_Discovery``; legacy 为项目根);
  ``stage_dir``    单个阶段的产物目录 (digest 目录键对应的目录);
  ``stage_state``  ``<项目>/90_Handoff/ai/stage_state/<pipeline>.<stage>.json``。

纯标准库, 可被任何模块安全 import。
"""

from __future__ import annotations

import json
import os
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional

try:
    from mmpv_common import io_layout, run_digest, stage_cards
except ImportError:  # 直接以脚本路径运行时的兜底
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
    from mmpv_common import io_layout, run_digest, stage_cards  # type: ignore


# ═══════════════════════════════════════════════════════════════════
# 常量: AI 交接目录约定
# ═══════════════════════════════════════════════════════════════════

AI_HANDOFF_DIR = "ai"
STAGE_STATE_SUBDIR = "stage_state"
AI_LOG_SUBDIR = "logs"
SUMMARIES_SUBDIR = "summaries"

# 阶段日志约定位置 (相对管线输出根): standard ③ 报告目录 / legacy 旧名
REPORT_LOG_DIR = "13_Report"
LEGACY_REPORT_LOG_DIR = "10_Reports"

DEFAULT_MAX_CHARS = 12000
LOG_TAIL_LINES = 100


# ═══════════════════════════════════════════════════════════════════
# stage 名 → io_layout 逻辑目录键
#
# 为什么需要这张表: stage_cards 的 stage 名是编排器 CLI 的 choices
# (assembly / identification / ...), digest 的键是 io_layout 逻辑键
# (d_asm / d_ident / ...) —— 两套命名服务不同读者, 这里做一次收敛。
# 一个 stage 复用多个目录键时取主产物目录 (public.convert 与 down 共用
# m_down); 没有独立目录的 stage (各类 report / submission) 留空, 摘要层
# 会明确写"无 digest 目录键"而不是臆造路径。
# ═══════════════════════════════════════════════════════════════════

STAGE_IO_KEYS: Dict[str, Dict[str, str]] = {
    "public": {
        "search": "m_search",
        "info": "m_info",
        "prescreen": "m_prescreen",
        "down": "m_down",
        "convert": "m_down",          # 近似: 转换产物仍在 03_RawData
        "plot": "m_plot",
        "hostref": "h_genome",
        "hostdb": "h_hostdb",
    },
    "preprocessing": {
        "clean": "d_clean",
        "deplete": "d_hostdep",
        "bbnorm": "d_bbnorm",
    },
    "discovery": {
        "clean": "d_clean",
        "deplete": "d_hostdep",
        "assembly": "d_asm",
        "identification": "d_ident",
        "filter": "d_filter",
        "cobra": "d_cobra",
        "merge": "d_merge",
        "cluster": "d_cluster",
        "taxonomy": "d_taxonomy",
        "host": "d_host_pred",
        "checkv": "d_checkv",
        "rescue": "d_rescue",
        "analysis": "d_analysis",
        "analysis_verify": "d_verify",
        "report": "d_reports",
    },
    "analysis": {
        "detect": "a_detect",
        "filter": "a_filter",
        "variants": "a_variants",
        "post": "a_post",
        "full": "a_assembly",
        # [阶段重排] extract = 比对与制备, 产物在 03_Extract
        "extract": "a_align",
        "similarity": "a_similarity",
        "dvg": "a_dvg",
        "report": "a_report",
    },
    "phylo": {
        "prep": "ph_data",
        "metadata": "ph_data",
        "online": "ph_data",
        "clean": "ph_data",
        "align": "ph_phylogeny",
        "splitstree": "ph_recomb",
        "rdp5": "ph_recomb",
        "tree": "ph_phylogeny",
        "popgen": "ph_popgen",
        "host": "ph_popgen",
        "capheine": "ph_select",
        "clock": "ph_time",
        "beast": "ph_time",
        "gene_dating": "ph_time",
        "phylogeo": "ph_geography",
        "geo": "ph_geography",
        "geo_paths": "ph_geography",
        "report": "ph_report",
    },
    "eve": {
        "discover": "e_loci",
        "verdict": "e_verdict",
        "annotate": "e_rvdb",
        "merge": "e_summary",
    },
    # ⑦ submission 消费 ④/③ 产物, 无独立输出根 → 无 digest 目录键
    "submission": {},
    "host_expression": {
        "quant": "he_quant",
        "labels": "he_labels",
        "coexp": "he_coexp",
        "de": "he_de",
        "loadscan": "he_loadscan",
        "enrich": "he_enrich",
        "report": "he_report",
    },
}


# ═══════════════════════════════════════════════════════════════════
# 路径解析 (全部容错)
# ═══════════════════════════════════════════════════════════════════

def resolve_layout(layout: Optional[str] = None) -> str:
    """布局名解析; 非法值回退默认布局 (上下文层永不因布局名崩)。"""
    try:
        return io_layout.resolve_layout(layout)
    except (ValueError, KeyError):
        return io_layout.DEFAULT_LAYOUT


def stage_io_key(pipeline: str, stage: str) -> Optional[str]:
    """stage 名 → io_layout 逻辑目录键; 未登记返回 None。"""
    return STAGE_IO_KEYS.get(pipeline, {}).get(stage)


def pipeline_root(project: Path, pipeline: str,
                  layout: Optional[str] = None) -> Optional[Path]:
    """管线输出根。

    standard: ``<项目>/<PIPELINE_ROOTS[pipeline]>``; legacy: 项目根 (目录平铺);
    无输出根的管线 (⑦ submission) 返回 None。
    """
    lay = resolve_layout(layout)
    project = Path(project)
    if pipeline not in run_digest.PIPELINE_ROOTS:
        return None
    if lay == "standard":
        sub = run_digest.PIPELINE_ROOTS[pipeline]
        return project / sub
    return project


def pipeline_output_dir(project: Path, pipeline: str,
                        layout: Optional[str] = None,
                        stage_state: Optional[Dict[str, Any]] = None) -> Optional[Path]:
    """编排器 ``--output_dir`` 语义的管线输出根。

    优先 stage_state 里登记的 ``output_dir`` (Phase 2 会写入实际路径), 其次
    按布局推导。这是 StageContext.output_stats / 13_Report 日志的锚点。
    """
    state = stage_state if isinstance(stage_state, dict) else {}
    raw = state.get("output_dir")
    if isinstance(raw, (str, Path)) and str(raw).strip():
        p = Path(raw)
        return p if p.is_absolute() else Path(project) / p
    return pipeline_root(project, pipeline, layout)


def stage_dir(project: Path, pipeline: str, stage: str,
              layout: Optional[str] = None) -> Optional[Path]:
    """阶段产物目录 (纯 digest 键映射, 不被 stage_state.output_dir 覆盖)。

    watch 用它对单个 stage 做 mtime 检测; 摘要语境下的 output_stats 走
    :func:`pipeline_output_dir`。
    """
    key = stage_io_key(pipeline, stage)
    if key is None:
        return None
    root = pipeline_root(project, pipeline, layout)
    if root is None:
        return None
    try:
        return root / io_layout.dir_name(key, layout)
    except (KeyError, ValueError):
        return None


def handoff_ai_dir(project: Path) -> Path:
    """``<项目>/90_Handoff/ai``。"""
    return Path(project) / run_digest.HANDOFF_DIR / AI_HANDOFF_DIR


def stage_state_dir(project: Path) -> Path:
    return handoff_ai_dir(project) / STAGE_STATE_SUBDIR


def ai_logs_dir(project: Path) -> Path:
    return handoff_ai_dir(project) / AI_LOG_SUBDIR


def summaries_dir(project: Path) -> Path:
    return handoff_ai_dir(project) / SUMMARIES_SUBDIR


def stage_state_path(project: Path, pipeline: str, stage: str) -> Path:
    return stage_state_dir(project) / f"{pipeline}.{stage}.json"


# ═══════════════════════════════════════════════════════════════════
# 安全读取
# ═══════════════════════════════════════════════════════════════════

def _read_json_dict(path: Path) -> Optional[Dict[str, Any]]:
    """读 JSON 对象; 不存在 / 非对象 / 解析失败一律 None。"""
    try:
        with open(path, "r", encoding="utf-8") as f:
            data = json.load(f)
    except (OSError, json.JSONDecodeError, UnicodeDecodeError, ValueError):
        return None
    return data if isinstance(data, dict) else None


def read_stage_state(project: Path, pipeline: str,
                     stage: str) -> Optional[Dict[str, Any]]:
    """读 ``90_Handoff/ai/stage_state/<pipeline>.<stage>.json``; 缺失返回 None。"""
    return _read_json_dict(stage_state_path(project, pipeline, stage))


def read_log_tail(path: Optional[Path], lines: int = LOG_TAIL_LINES) -> str:
    """安全读取文件最后 ``lines`` 行 (按字节从尾部回溯, 不整文件入内存)。

    文件不存在 / 不可读 / 为空 → ``""``。编码用 utf-8 + replace, 坏字节
    不会抛异常。返回文本不含末尾换行。
    """
    if lines <= 0 or path is None:
        return ""
    p = Path(path)
    try:
        if not p.is_file():
            return ""
        with open(p, "rb") as f:
            f.seek(0, os.SEEK_END)
            pos = f.tell()
            chunk = 16 * 1024
            data = b""
            while pos > 0 and data.count(b"\n") <= lines:
                step = min(chunk, pos)
                pos -= step
                f.seek(pos)
                data = f.read(step) + data
                chunk = min(chunk * 2, 1024 * 1024)
    except OSError:
        return ""
    text = data.decode("utf-8", errors="replace")
    if pos > 0:
        # 从文件中间开始读: 首行可能是半行, 丢弃
        nl = text.find("\n")
        text = text[nl + 1:] if nl != -1 else ""
    all_lines = text.splitlines()
    return "\n".join(all_lines[-lines:])


def dir_stats(path: Optional[Path], max_subdirs: int = 32) -> Dict[str, Any]:
    """目录统计: 存在性 / 文件数 / 字节 / 子目录名 (直接子项, 不递归)。"""
    if path is None:
        return {"exists": False, "root": ""}
    p = Path(path)
    out: Dict[str, Any] = {"exists": False, "root": str(p)}
    try:
        if not p.is_dir():
            return out
        files = 0
        total = 0
        subdirs: List[str] = []
        for entry in p.iterdir():
            if entry.is_dir():
                subdirs.append(entry.name)
            elif entry.is_file():
                files += 1
                try:
                    total += entry.stat().st_size
                except OSError:
                    pass
    except OSError as e:
        out.update({"exists": True, "error": str(e)})
        return out
    out.update({"exists": True, "files": files, "bytes": total,
                "subdirs": sorted(subdirs)[:max_subdirs],
                "n_subdirs": len(subdirs)})
    return out


def _safe_glob(directory: Path, pattern: str) -> List[Path]:
    try:
        return sorted(directory.glob(pattern))
    except OSError:
        return []


# ═══════════════════════════════════════════════════════════════════
# StageContext
# ═══════════════════════════════════════════════════════════════════

@dataclass
class StageContext:
    """一个阶段的可摘要上下文 (纯数据, JSON 可序列化)。

    output_stats 统计的是**管线输出根** (编排器 --output_dir 语义), 不是
    单个 stage 目录 —— 阶段级文件数在 ``digest["pipelines"][pipeline]
    ["stages"][io_key]`` 里, 模板摘要以 digest 为准, 这里只补充根目录地貌。
    """

    project: Path
    pipeline: str
    stage: str
    digest: Dict[str, Any]
    card: Optional[Dict[str, Any]]
    stage_state: Optional[Dict[str, Any]]
    log_tail: str
    output_stats: Dict[str, Any]

    def to_dict(self, max_chars: int = DEFAULT_MAX_CHARS) -> Dict[str, Any]:
        """JSON-serializable 视图; 日志超长截断并标注 ``truncated``。"""
        tail = self.log_tail or ""
        if max_chars <= 0:
            truncated = bool(tail)
            tail = ""
        else:
            truncated = len(tail) > max_chars
            if truncated:
                tail = tail[-max_chars:]   # 尾部 (最新日志) 优先
        return {
            "project": str(self.project),
            "pipeline": self.pipeline,
            "stage": self.stage,
            "digest": self.digest,
            "card": self.card,
            "stage_state": self.stage_state,
            "log_tail": tail,
            "output_stats": self.output_stats,
            "truncated": truncated,
        }


# ═══════════════════════════════════════════════════════════════════
# ContextBuilder
# ═══════════════════════════════════════════════════════════════════

class ContextBuilder:
    """按需组装 digest / stage / project 上下文。所有方法不抛异常。"""

    def __init__(self, project: Path, layout: Optional[str] = None):
        self.project_dir = Path(project)
        self._layout_arg = layout

    # ── 布局 ────────────────────────────────────────────────────

    def _layout(self) -> str:
        return resolve_layout(self._layout_arg)

    # ── digest ──────────────────────────────────────────────────

    def digest(self, compact: bool = True) -> Dict[str, Any]:
        """项目 run digest; compact=True 时用 ``compact_for_ai`` 瘦身。

        build_digest 自身已逐管线容错; 这里再兜一层, 保证任何情况下都
        返回带 ``pipelines`` 的 dict (失败时附 ``error``)。
        """
        lay = self._layout()
        try:
            d = run_digest.build_digest(self.project_dir, layout=lay)
        except Exception as e:  # noqa: BLE001 — 上下文层永不抛
            d = {
                "schema_version": run_digest.SCHEMA_VERSION,
                "generated_at": time.strftime("%Y-%m-%d %H:%M:%S"),
                "project_dir": str(self.project_dir),
                "layout": lay,
                "pipelines": {},
                "totals": {"pipelines_present": 0, "pipelines_present_list": [],
                           "project_exists": False},
                "error": f"{type(e).__name__}: {e}",
            }
        if compact:
            try:
                return run_digest.compact_for_ai(d)
            except Exception:  # noqa: BLE001 — 瘦身失败就用全量
                return d
        return d

    # ── 单个 stage ──────────────────────────────────────────────

    def _log_path_candidates(self, pipeline: str, stage: str,
                             output_dir: Optional[Path],
                             state: Optional[Dict[str, Any]]) -> List[Path]:
        """日志候选, 顺序即优先级 (契约):
        stage_state["log_path"] → <output_dir>/13_Report/logs/<stage>.log →
        <项目>/90_Handoff/ai/logs/<pipeline>.<stage>.log; 末尾追加 legacy
        报告目录 10_Reports 兜底 (不改变上面三者的优先级)。
        """
        cands: List[Path] = []
        if isinstance(state, dict):
            raw = state.get("log_path")
            if isinstance(raw, (str, Path)) and str(raw).strip():
                p = Path(raw)
                cands.append(p if p.is_absolute() else self.project_dir / p)
        if output_dir is not None:
            cands.append(Path(output_dir) / REPORT_LOG_DIR / "logs" / f"{stage}.log")
        cands.append(ai_logs_dir(self.project_dir) / f"{pipeline}.{stage}.log")
        if output_dir is not None:
            cands.append(Path(output_dir) / LEGACY_REPORT_LOG_DIR / "logs" / f"{stage}.log")
        return cands

    def stage(self, pipeline: str, stage: str,
              log_lines: int = LOG_TAIL_LINES) -> StageContext:
        """组装一个阶段的上下文 (digest + 卡片 + stage_state + 日志 + 目录统计)。"""
        digest = self.digest(compact=False)
        layout = self._layout()
        try:
            card = stage_cards.get_card(pipeline, stage)
        except Exception:  # noqa: BLE001 — 卡片查找失败按未知处理
            card = None
        state = read_stage_state(self.project_dir, pipeline, stage)
        out_dir = pipeline_output_dir(self.project_dir, pipeline, layout, state)
        stats = dir_stats(out_dir)
        sdir = stage_dir(self.project_dir, pipeline, stage, layout)
        stats["stage_dir"] = str(sdir) if sdir is not None else ""

        log_tail = ""
        for cand in self._log_path_candidates(pipeline, stage, out_dir, state):
            text = read_log_tail(cand, log_lines)
            if text:
                log_tail = text
                break
        return StageContext(
            project=self.project_dir, pipeline=pipeline, stage=stage,
            digest=digest, card=card, stage_state=state,
            log_tail=log_tail, output_stats=stats,
        )

    # ── 项目级 ──────────────────────────────────────────────────

    def project(self) -> Dict[str, Any]:
        """项目级上下文: digest + 各管线进度 + AI 交接文件清单。"""
        d = self.digest(compact=False)
        pipelines: Dict[str, Any] = {}
        for pk, sec in (d.get("pipelines") or {}).items():
            if not isinstance(sec, dict):
                continue
            stages = sec.get("stages") if isinstance(sec.get("stages"), dict) else {}
            present = sorted(k for k, v in stages.items()
                             if isinstance(v, dict) and v.get("exists"))
            counts: Dict[str, Any] = {}
            for k, v in sec.items():
                if k in ("stages", "present"):
                    continue
                if isinstance(v, bool) or not isinstance(v, (int, float)):
                    continue
                counts[k] = v
            pipelines[pk] = {
                "present": bool(sec.get("present")),
                "n_stages_present": len(present),
                "stages_present": present,
                "counts": counts,
            }
        return {
            "project": str(self.project_dir),
            "layout": d.get("layout") or self._layout(),
            "generated_at": d.get("generated_at") or "",
            "digest": d,
            "totals": d.get("totals") or {},
            "pipelines": pipelines,
            "stage_states": [p.name for p in _safe_glob(stage_state_dir(self.project_dir), "*.json")],
            "summaries": [p.name for p in _safe_glob(summaries_dir(self.project_dir), "*.json")],
        }


# 兼容别名 (集成分支可能按模块级函数取上下文)
def build_digest(project: Path, compact: bool = True,
                 layout: Optional[str] = None) -> Dict[str, Any]:
    """便捷函数: ``ContextBuilder(project, layout).digest(compact)``。"""
    return ContextBuilder(project, layout=layout).digest(compact=compact)


if __name__ == "__main__":  # 手工排查: python -m mmpv_common.context <项目> <管线> <stage>
    import argparse

    ap = argparse.ArgumentParser(description="MMPV 阶段上下文自检")
    ap.add_argument("project")
    ap.add_argument("pipeline", nargs="?", default=None)
    ap.add_argument("stage", nargs="?", default=None)
    ap.add_argument("--io-layout", choices=list(io_layout.LAYOUTS), default=None)
    ns = ap.parse_args()
    b = ContextBuilder(Path(ns.project), layout=ns.io_layout)
    payload = (b.stage(ns.pipeline, ns.stage).to_dict()
               if ns.pipeline and ns.stage else b.project())
    print(json.dumps(payload, ensure_ascii=False, indent=2))
