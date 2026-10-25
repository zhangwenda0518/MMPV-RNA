#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""元数据管线的 AI 溯源日志 (与 mmpv_common.ai.log_ai_event 同一 schema)。

为什么需要
----------
平台 AI 驱动层 (mmpv_common/ai.py) 把所有 AI 调用记进
``<项目>/90_Handoff/state/ai_log.jsonl``, ``mmpv_ai.py provenance`` 与
过程报告据此统计「AI 参与」。但 public_metadata_pipeline 的 AI 调用是内嵌
的 (search 洗髓 / 检索回填 / summary, info 推断 / 仲裁 / 正交化 / Location
规范化), 从未写这份日志 —— 元数据管线的 AI 参与在溯源里完全隐形。

本模块以同一份事件 schema 打点, 写入位置按项目布局自动反推:

* ``<outdir>`` 或向上 3 级父目录存在 ``90_Handoff/`` -> ``<项目>/90_Handoff/state``
  (这样 provenance / 过程报告能直接聚合到);
* 否则退化写到 ``<outdir>/ai_log.jsonl`` 本身 (仍可离线审计, 只是不计入
  项目级溯源, 命令行会打一行提示)。

纪律
----
* 打点函数 **绝不抛异常** (不能因为溯源失败拖垮主流程);
* 事件只记「元信息」(op / 模型 / prompt 版本与哈希 / 计数 / token),
  不记原始数据内容 —— 与平台层「AI 只见摘要」的红线一致。
"""

from __future__ import annotations

import hashlib
import json
import os
import sys
import time
from pathlib import Path
from typing import Any, Dict, Optional

#: 内嵌提示词版本 (改措辞必须升版本; 事件里同时记 sha 便于逐字核对)
PROMPT_VERSIONS: Dict[str, str] = {
    "metadata_sanitizer": "1.1.0",
    "metadata_search_backfill": "1.0.0",
    "metadata_summary": "1.0.0",
    "metadata_sra_infer": "1.0.0",
    "metadata_sra_arbitrate": "1.0.0",
    "metadata_gsa_enhance": "1.0.0",
    "metadata_location_normalize": "1.0.0",
    "metadata_orthogonalize": "1.1.0",
}

_DEFAULT_STATE: Optional[Path] = None


# ═══════════════════════════════════════════════════════════════════
# 目录解析
# ═══════════════════════════════════════════════════════════════════

def resolve_state_dir(outdir: Any) -> Optional[Path]:
    """从产物目录反推项目 ``90_Handoff/state``; 反推不到则返回 outdir 本身。"""
    if outdir is None:
        return None
    try:
        p = Path(str(outdir)).resolve()
    except (OSError, ValueError):
        return None
    for cand in (p, *list(p.parents)[:3]):
        try:
            if (cand / "90_Handoff").is_dir():
                return cand / "90_Handoff" / "state"
        except OSError:
            continue
    return p


def set_default_state_dir(path: Optional[Path]) -> None:
    global _DEFAULT_STATE
    _DEFAULT_STATE = Path(path) if path is not None else None


def get_state_dir() -> Optional[Path]:
    return _DEFAULT_STATE


def auto_state(outdir: Any) -> Optional[Path]:
    """按产物目录 (或 MMPV_AI_STATE_DIR 环境变量) 设定默认溯源目录."""
    env = os.environ.get("MMPV_AI_STATE_DIR")
    if env:
        try:
            set_default_state_dir(Path(env))
            return Path(env)
        except (OSError, ValueError):
            pass
    d = resolve_state_dir(outdir)
    set_default_state_dir(d)
    return d


# ═══════════════════════════════════════════════════════════════════
# 写出
# ═══════════════════════════════════════════════════════════════════

def prompt_sha(text: Optional[str]) -> str:
    if not text:
        return ""
    return hashlib.sha256(text.encode("utf-8")).hexdigest()[:16]


def _platform_logger():
    """优先复用平台层 log_ai_event (保证 schema 完全一致)。"""
    try:
        repo = Path(__file__).resolve().parents[1]
        if str(repo) not in sys.path:
            sys.path.insert(0, str(repo))
        from mmpv_common.ai import log_ai_event  # type: ignore
        return log_ai_event
    except Exception:
        return None


def _local_log(path: Path, body: Dict[str, Any]) -> Optional[Path]:
    try:
        path.parent.mkdir(parents=True, exist_ok=True)
        with open(path, "a", encoding="utf-8") as f:
            f.write(json.dumps(body, ensure_ascii=False) + "\n")
        return path
    except OSError:
        return None


def trace(op: str, *, model: str = "", ok: bool = True, latency_ms: int = 0,
          usage: Optional[Dict[str, Any]] = None, cached: bool = False,
          prompt_text: Optional[str] = None, version: Optional[str] = None,
          state_dir: Optional[Path] = None, **extra: Any) -> Optional[Path]:
    """追加一条 AI 调用事件; 任何失败都静默返回 None。"""
    try:
        target = Path(state_dir) if state_dir is not None else _DEFAULT_STATE
        if target is None:
            return None
        event: Dict[str, Any] = {
            "kind": "call",
            "op": op,
            "model": model or "",
            "ok": bool(ok),
            "latency_ms": int(latency_ms or 0),
            "cached": bool(cached),
            "source": "public_metadata_pipeline",
            "prompt": {
                "name": op,
                "version": version or PROMPT_VERSIONS.get(op, "0.0.0"),
            },
            "usage": dict(usage or {}),
        }
        sha = prompt_sha(prompt_text)
        if sha:
            event["prompt"]["sha"] = sha
        for k, v in extra.items():
            if v is not None:
                event[k] = v
        logger = _platform_logger()
        if logger is not None:
            return logger(target, event)
        return _local_log(Path(target) / "ai_log.jsonl",
                          {"ts": time.strftime("%Y-%m-%d %H:%M:%S"), **event})
    except Exception:
        return None


def usage_tokens(res: Any) -> Dict[str, Any]:
    """从 openai 响应对象里取 token 用量 (取不到返回空 dict)。"""
    try:
        u = getattr(res, "usage", None)
        if u is None:
            return {}
        total = getattr(u, "total_tokens", None)
        if total is None:
            return {}
        return {"total_tokens": int(total)}
    except Exception:
        return {}
