#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ai_hooks.py — 阶段级 AI 摘要钩子 (Phase 2)
==========================================

每条编排管线跑完一个 stage 后调用一次 :func:`after_stage`, 把「这一阶段
到底产出了什么、成功还是失败、花了多久、日志在哪」固化成:

  <项目>/90_Handoff/ai/stage_state/<pipeline>.<stage>.json   机器可读状态
  <项目>/90_Handoff/ai/stage_summaries/<pipeline>.<stage>.md 人读摘要
  <项目>/90_Handoff/ai/index.jsonl                           追加式时间线

摘要由 ``mmpv_common.summaries`` 生成 (lazy import); 该模块缺失、接口
不匹配或抛异常时, 本模块自己写一份**离线确定性模板** —— 钩子绝不调用
LLM 之外的任何东西, 也绝不抛异常打断主管线。

三种模式 (环境变量 ``MMPV_AI_STAGE_SUMMARY``, 大小写不敏感):

  off       完全不写任何文件 (想彻底关闭时用)
  template  默认。只写 stage_state + 模板摘要 + index, 零 AI 调用
  ai        在 template 基础上再让 summaries 生成 AI 摘要 (需已配置 API key;
            生成失败时 summaries 或本模块回退到模板, 永不阻断)

设计纪律:
  * best-effort: 任何异常都吞掉并反映在返回值里, 不改管线退出码/stdout;
  * 所有 import 尽量 lazy —— 管线在跑子进程前 import 本模块, 不能因为
    某个可选依赖缺失而启动失败;
  * 路径规范化: 传入各管线自己的输出根 (如 ``03_Discovery``) 时自动上提
    到项目根, 保证 90_Handoff 落在项目根而不是管线根里。

CLI (调试/回填):
  python -m mmpv_common.ai_hooks --backfill <项目根> [--pipeline discovery]
"""

from __future__ import annotations

import datetime as _dt
import json
import os
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

# ═══════════════════════════════════════════════════════════════════
# 常量
# ═══════════════════════════════════════════════════════════════════

SCHEMA_VERSION = 1
ENV_MODE = "MMPV_AI_STAGE_SUMMARY"
MODE_OFF = "off"
MODE_TEMPLATE = "template"
MODE_AI = "ai"
MODES = (MODE_OFF, MODE_TEMPLATE, MODE_AI)
DEFAULT_MODE = MODE_TEMPLATE

HANDOFF_DIR = "90_Handoff"
AI_SUBDIR = "ai"
STATE_SUBDIR = "stage_state"
SUMMARY_SUBDIR = "stage_summaries"
INDEX_NAME = "index.jsonl"

# 单阶段目录扫描上限 (病态目录不拖垮管线)
MAX_SCAN_FILES = 200_000

# 管线根名 → 项目根上提用 (与 run_digest.PIPELINE_ROOTS 同源; 这里复制一份
# 仅作 import 失败时的兜底, 正常路径用 run_digest 的权威表)
_PIPELINE_ROOT_FALLBACK = {
    "public": "01_PublicData",
    "preprocessing": "02_Preprocessing",
    "discovery": "03_Discovery",
    "analysis": "04_Analysis",
    "phylo": "05_Phylo",
    "eve": "06_EVE",
    "host_expression": "08_HostExpression",
}

# 编排器 stage 名 → io_layout 逻辑键 (digest_keys 用; 与各管线 CLI 的真实
# stage 名对齐, 键不存在于 KEYS_BY_PIPELINE 的管线会被忽略)
STAGE_KEYS: Dict[str, Dict[str, str]] = {
    "public": {
        "search": "m_search", "info": "m_info", "prescreen": "m_prescreen",
        "down": "m_down", "download": "m_down", "plot": "m_plot",
        "hostref": "h_genome", "hostref_genome": "h_genome",
        "hostdb": "h_hostdb",
    },
    "preprocessing": {
        "clean": "d_clean", "deplete": "d_hostdep", "bbnorm": "d_bbnorm",
    },
    "discovery": {
        "clean": "d_clean", "deplete": "d_hostdep", "assembly": "d_asm",
        "identification": "d_ident", "filter": "d_filter", "graph": "d_graph",
        "cobra": "d_cobra", "merge": "d_merge", "cluster": "d_cluster",
        "taxonomy": "d_taxonomy", "host": "d_host_pred", "checkv": "d_checkv",
        "rescue": "d_rescue", "analysis": "d_analysis",
        "analysis_verify": "d_verify", "report": "d_reports",
    },
    "analysis": {
        "detect": "a_detect", "filter": "a_filter", "variants": "a_variants",
        "post": "a_post", "assembly": "a_assembly", "full": "a_assembly",
        # [阶段重排] extract = 靶向 reads 提取, 产物在 03_Variants/virus_reads
        "extract": "a_variants", "similarity": "a_similarity", "dvg": "a_dvg",
        "report": "a_report",
    },
    "phylo": {
        "data": "ph_data", "phylogeny": "ph_phylogeny", "popgen": "ph_popgen",
        "recomb": "ph_recomb", "select": "ph_select", "time": "ph_time",
        "geography": "ph_geography", "report": "ph_report",
    },
    "eve": {
        "loci": "e_loci", "verdict": "e_verdict", "rvdb": "e_rvdb",
        "summary": "e_summary",
    },
    "host_expression": {
        "hostref": "he_hostref", "quant": "he_quant", "qc": "he_qc",
        "labels": "he_labels", "coexp": "he_coexp", "de": "he_de",
        "enrich": "he_enrich", "loadscan": "he_loadscan", "report": "he_report",
    },
}

# 反向表: 单键 → 规范 stage 名 (首次出现者胜; 用于 backfill 枚举)
_KEY_TO_STAGE: Dict[str, Dict[str, str]] = {}
for _pk, _m in STAGE_KEYS.items():
    _rev: Dict[str, str] = {}
    for _st, _key in _m.items():
        _rev.setdefault(_key, _st)
    _KEY_TO_STAGE[_pk] = _rev

_PREFIXES = ("d_", "a_", "m_", "h_", "ph_", "e_", "he_")


# ═══════════════════════════════════════════════════════════════════
# 基础工具
# ═══════════════════════════════════════════════════════════════════

def resolve_mode(env: Optional[Dict[str, str]] = None) -> str:
    """解析 ``MMPV_AI_STAGE_SUMMARY``; 大小写不敏感, 非法/空值回退 template。"""
    e = env if env is not None else os.environ
    raw = str(e.get(ENV_MODE) or "").strip().lower()
    return raw if raw in MODES else DEFAULT_MODE


def project_root(project_dir: Any) -> Path:
    """把「管线输出根」规范化为项目根。

    standard 布局下编排器 --output_dir 指向 ``03_Discovery`` 这类独立根;
    90_Handoff 属于项目根, 因此名称命中管线根名时上提一层。
    legacy 布局下输出根就是项目根, 行为不变。
    """
    p = Path(project_dir).expanduser().resolve()
    try:
        from mmpv_common.run_digest import PIPELINE_ROOTS as _roots
        names = set(_roots.values())
    except Exception:  # pragma: no cover - run_digest 始终存在, 纯兜底
        names = set(_PIPELINE_ROOT_FALLBACK.values())
    if p.name in names and p.parent != p:
        return p.parent
    return p


def ai_dir(project: Path) -> Path:
    """``<项目>/90_Handoff/ai``。"""
    return Path(project) / HANDOFF_DIR / AI_SUBDIR


def state_path(project: Path, pipeline: str, stage: str) -> Path:
    """阶段状态 JSON 的规范路径。"""
    return ai_dir(project) / STATE_SUBDIR / f"{_safe(pipeline)}.{_safe(stage)}.json"


def summary_path(project: Path, pipeline: str, stage: str) -> Path:
    """本模块回退模板摘要的规范路径。"""
    return ai_dir(project) / SUMMARY_SUBDIR / f"{_safe(pipeline)}.{_safe(stage)}.md"


def io_keys_for_stage(pipeline: str, stage: str) -> List[str]:
    """stage 名 → 该管线 digest 里的 io_layout 逻辑键 (可能为空)。"""
    if not pipeline or not stage:
        return []
    table = STAGE_KEYS.get(str(pipeline), {})
    key = table.get(str(stage).strip().lower())
    if key:
        return [key]
    # 允许直接传逻辑键 (如 d_clean)
    try:
        from mmpv_common.run_digest import KEYS_BY_PIPELINE
        if str(stage) in KEYS_BY_PIPELINE.get(str(pipeline), []):
            return [str(stage)]
    except Exception:
        pass
    return []


def stage_for_io_key(pipeline: str, key: str) -> str:
    """io_layout 逻辑键 → 规范 stage 名 (未登记时剥前缀兜底)。"""
    rev = _KEY_TO_STAGE.get(str(pipeline), {})
    if key in rev:
        return rev[key]
    k = str(key)
    for pref in _PREFIXES:
        if k.startswith(pref):
            return k[len(pref):]
    return k


def scan_output_dir(path: Any, max_files: int = MAX_SCAN_FILES) -> Dict[str, Any]:
    """递归统计目录: 文件数/字节/最后修改时间。不存在 → exists=False。

    单文件路径按 1 个文件处理; 任何 OSError 只跳过该条目, 绝不抛。
    """
    out: Dict[str, Any] = {"exists": False, "n_files": 0, "bytes": 0,
                           "last_modified": None}
    if path is None:
        return out
    p = Path(path)
    try:
        if not p.exists():
            return out
        out["exists"] = True
        if p.is_file():
            st = p.stat()
            out.update(n_files=1, bytes=st.st_size,
                       last_modified=_fmt_ts(st.st_mtime))
            return out
        latest = 0.0
        stack: List[Path] = [p]
        scanned = 0
        while stack:
            d = stack.pop()
            try:
                with os.scandir(d) as it:
                    for entry in it:
                        try:
                            if entry.is_dir(follow_symlinks=False):
                                stack.append(Path(entry.path))
                            elif entry.is_file(follow_symlinks=False):
                                st = entry.stat(follow_symlinks=False)
                                out["n_files"] += 1
                                out["bytes"] += int(st.st_size)
                                latest = max(latest, float(st.st_mtime))
                                scanned += 1
                                if scanned >= max_files:
                                    stack.clear()
                                    break
                        except OSError:
                            continue
            except OSError:
                continue
        if latest:
            out["last_modified"] = _fmt_ts(latest)
    except OSError:
        pass
    return out


def call_supported(fn: Any, candidates: Dict[str, Any]) -> Any:
    """按目标函数签名挑它能接受的参数调用 (供 summaries/diagnosis 等并行模块)。

    只传签名里真实存在且提供了候选值的参数; 必需参数既不在候选值里也
    没有默认值时抛 TypeError, 由调用方决定是否降级。
    """
    import inspect

    sig = inspect.signature(fn)
    args: List[Any] = []
    kwargs: Dict[str, Any] = {}
    for name, param in sig.parameters.items():
        if param.kind in (param.VAR_POSITIONAL, param.VAR_KEYWORD):
            continue
        if name in candidates:
            val = candidates[name]
            if param.kind is param.POSITIONAL_ONLY:
                args.append(val)
            else:
                kwargs[name] = val
        elif param.default is inspect.Parameter.empty:
            raise TypeError(f"{getattr(fn, '__name__', 'fn')} 缺少必需参数: {name}")
    return fn(*args, **kwargs)


def _safe(name: str) -> str:
    """文件名安全化 (stage/pipeline 不应带路径分隔符)。"""
    s = str(name or "unknown").strip() or "unknown"
    for ch in ("/", "\\", ":", "*", "?", '"', "<", ">", "|"):
        s = s.replace(ch, "_")
    return s


def _fmt_ts(ts: Optional[float]) -> Optional[str]:
    if not ts:
        return None
    try:
        return _dt.datetime.fromtimestamp(float(ts)).strftime("%Y-%m-%d %H:%M:%S")
    except (OSError, OverflowError, ValueError):
        return None


def _now() -> str:
    return _dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def _as_float(value: Any) -> Optional[float]:
    if value is None:
        return None
    try:
        return round(float(value), 3)
    except (TypeError, ValueError):
        return None


def _error_fields(error: Any) -> Tuple[Optional[str], Optional[str]]:
    """异常/字符串 → (error_type, error_message)。"""
    if error is None:
        return None, None
    if isinstance(error, BaseException):
        msg = str(error) or repr(error)
        return type(error).__name__, msg
    return "Error", str(error)


def _write_json(path: Path, payload: Dict[str, Any]) -> None:
    """原子写 JSON (先写临时文件再替换)。"""
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(json.dumps(payload, ensure_ascii=False, indent=2),
                   encoding="utf-8")
    os.replace(tmp, path)


def _write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(text, encoding="utf-8")
    os.replace(tmp, path)


def _append_index(project: Path, record: Dict[str, Any]) -> Path:
    p = ai_dir(project) / INDEX_NAME
    p.parent.mkdir(parents=True, exist_ok=True)
    with open(p, "a", encoding="utf-8") as f:
        f.write(json.dumps(record, ensure_ascii=False) + "\n")
    return p


# ═══════════════════════════════════════════════════════════════════
# 摘要生成 (summaries 优先, 模板兜底)
# ═══════════════════════════════════════════════════════════════════

_SUMMARY_TEXT_KEYS = ("markdown", "md", "text", "summary", "content", "body")
_SUMMARY_PATH_KEYS = ("summary_path", "path", "output_path", "md_path", "file")


def _extract_text(payload: Any) -> Optional[str]:
    """从 summaries 返回值里提取 Markdown 文本 (str 或 dict 常见键)。"""
    if isinstance(payload, str):
        return payload.strip() or None
    if isinstance(payload, dict):
        for key in _SUMMARY_TEXT_KEYS:
            value = payload.get(key)
            if isinstance(value, str) and value.strip():
                return value.strip()
    return None


def _extract_path(payload: Any) -> Optional[Path]:
    if isinstance(payload, dict):
        for key in _SUMMARY_PATH_KEYS:
            value = payload.get(key)
            if isinstance(value, (str, os.PathLike)) and str(value).strip():
                p = Path(str(value))
                # summaries 返回的是 JSON 路径时, 优先指向同名 .md (人读摘要)
                if p.suffix.lower() == ".json":
                    md = p.with_suffix(".md")
                    if md.is_file():
                        return md
                return p
    return None


def _try_summaries(project: Path, pipeline: str, stage: str, *, mode: str,
                   state: Dict[str, Any], output_dir: Optional[Path],
                   config: Any, client: Any) -> Tuple[Optional[str], Optional[Path], str]:
    """尝试调用 ``mmpv_common.summaries.summarize_stage``。

    返回 (markdown, 已存在的摘要路径, 来源标记)。任何失败都返回 (None, None, why)。
    """
    try:
        from mmpv_common import summaries  # lazy import: 并行开发中, 可能不存在
    except Exception as e:  # ImportError / 语法错误 / 循环导入
        return None, None, f"fallback:summaries_unavailable({type(e).__name__})"
    fn = getattr(summaries, "summarize_stage", None)
    if not callable(fn):
        return None, None, "fallback:summarize_stage_missing"

    candidates: Dict[str, Any] = {
        "project_dir": project,
        "project": project,
        "project_path": project,
        "root": project,
        "pipeline": pipeline,
        "stage": stage,
        "level": "stage",
        "mode": mode,
        "ai": mode == MODE_AI,
        "use_ai": mode == MODE_AI,
        "output_dir": output_dir,
        "state": state,
        "stage_state": state,
        "state_obj": state,
        "config": config,
        "client": client,
    }
    try:
        payload = call_supported(fn, candidates)
    except TypeError:
        # 宽松重试: 按最保守的位置参数形态调用
        try:
            payload = fn(project, pipeline, stage)
        except Exception as e:  # noqa: BLE001 - 摘要失败必须降级
            return None, None, f"fallback:{type(e).__name__}"
    except Exception as e:  # noqa: BLE001
        return None, None, f"fallback:{type(e).__name__}"

    if hasattr(payload, "__await__"):  # 异步接口不支持, 降级
        return None, None, "fallback:async_summary_unsupported"

    text = _extract_text(payload)
    path = _extract_path(payload)
    source = "ai" if mode == MODE_AI else "template:summaries"
    if path is not None and text is not None and not path.exists():
        # summaries 返回了文本 + 目标路径但尚未落盘: 由本模块写入
        try:
            _write_text(path if path.is_absolute() else project / path, text)
        except OSError:
            path = None
    if text is None and path is not None and path.exists():
        try:
            text = path.read_text(encoding="utf-8", errors="replace")
        except OSError:
            text = None
    if text is None and path is None:
        return None, None, "fallback:empty_summary"
    return text, path, source


def _fallback_markdown(pipeline: str, stage: str, state: Dict[str, Any],
                       n_files: int, n_bytes: int, log_path: Optional[Path],
                       error_type: Optional[str], error_message: Optional[str]) -> str:
    """离线确定性模板摘要 (零 AI 调用)。"""
    lines: List[str] = [
        f"# 阶段摘要: {pipeline} / {stage}",
        "",
        f"- **状态**: `{state.get('status')}`",
        f"- **完成时间**: {state.get('finished_at') or '未知'}",
        f"- **耗时**: {state.get('elapsed_s') if state.get('elapsed_s') is not None else '未知'} s",
        f"- **输出目录**: `{state.get('output_dir') or '未指定'}`",
        f"- **产物规模**: {n_files} 个文件 / {_fmt_bytes(n_bytes)}",
        f"- **日志**: `{log_path}`" if log_path else "- **日志**: 未定位到阶段日志",
    ]
    if error_type or error_message:
        lines += ["", f"- **错误类型**: `{error_type or 'Error'}`",
                  f"- **错误信息**: {error_message or ''}"]
    keys = state.get("digest_keys") or []
    if keys:
        lines.append(f"- **digest 键**: {', '.join(keys)}")
    lines += [
        "",
        "> 本摘要由 `mmpv_common.ai_hooks` 离线模板生成 (未调用 LLM)。"
        "配置 `MMPV_AI_STAGE_SUMMARY=ai` 且 `mmpv_common.summaries` 可用时，"
        "同一阶段会追加 AI 解读。",
        "",
    ]
    return "\n".join(lines)


def _fmt_bytes(n: int) -> str:
    value = float(max(int(n or 0), 0))
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if value < 1024.0 or unit == "TB":
            return f"{value:.1f} {unit}" if unit != "B" else f"{int(value)} B"
        value /= 1024.0
    return f"{int(n)} B"


# ═══════════════════════════════════════════════════════════════════
# 对外接口 1: after_stage
# ═══════════════════════════════════════════════════════════════════

def after_stage(project_dir: Any, pipeline: str, stage: str, *, status: str,
                output_dir: Any = None, log_path: Any = None,
                elapsed_s: Any = None, error: Any = None,
                config: Any = None, client: Any = None) -> Dict[str, Any]:
    """记录一个 stage 的结束状态并生成摘要 (best-effort, 永不抛)。

    参数:
        project_dir: 项目根, 或某管线自己的输出根 (如 03_Discovery; 自动上提)
        pipeline:    管线键 (public/preprocessing/discovery/analysis/phylo/...)
        stage:       阶段名 (clean/assembly/rescue/...)
        status:      调用方给的结局标签 (ok/failed/error/skipped/unknown ...)
        output_dir:  该阶段的输出目录 (统计文件数/字节)
        log_path:    该阶段的日志文件 (写进状态, 供 MCP 读日志尾)
        elapsed_s:   阶段耗时 (秒)
        error:       失败时的异常对象或字符串
        config/client: 传给 summaries 的可选 AI 连接 (template 模式不使用)

    返回 {"status": "ok|skipped|error", "mode": ..., "state_path": ...,
          "summary_path": ..., "summary_source": ...}
    """
    mode = resolve_mode()
    result: Dict[str, Any] = {"status": "skipped", "mode": mode,
                              "state_path": None, "summary_path": None,
                              "summary_source": None}
    if mode == MODE_OFF:
        result["reason"] = f"{ENV_MODE}=off (钩子关闭)"
        return result
    try:
        if project_dir is None or str(project_dir).strip() == "":
            result.update(status="error", reason="project_dir 为空")
            return result
        project = project_root(project_dir)
        if not project.is_dir():
            result.update(status="error",
                          reason=f"项目目录不存在: {project}")
            return result

        pipeline = str(pipeline or "unknown").strip() or "unknown"
        stage = str(stage or "unknown").strip() or "unknown"
        status = str(status or "unknown").strip() or "unknown"

        finished = _now()
        elapsed = _as_float(elapsed_s)
        started = None
        if elapsed is not None:
            try:
                started = (_dt.datetime.now()
                           - _dt.timedelta(seconds=elapsed)).strftime("%Y-%m-%d %H:%M:%S")
            except (OverflowError, ValueError):
                started = None

        out_path = Path(output_dir) if output_dir else None
        stats = scan_output_dir(out_path)
        logp = Path(log_path) if log_path else _guess_log_path(out_path, stage)
        err_type, err_msg = _error_fields(error)
        keys = io_keys_for_stage(pipeline, stage)

        state: Dict[str, Any] = {
            "schema_version": SCHEMA_VERSION,
            "pipeline": pipeline,
            "stage": stage,
            "status": status,
            "started_at": started,
            "finished_at": finished,
            "elapsed_s": elapsed,
            "output_dir": str(out_path) if out_path else None,
            "n_files": int(stats["n_files"]),
            "bytes": int(stats["bytes"]),
            "log_path": str(logp) if logp else None,
            "error_type": err_type,
            "error_message": err_msg,
            "digest_keys": keys,
        }

        spath = state_path(project, pipeline, stage)
        _write_json(spath, state)

        # ── 摘要: summaries 优先, 模板兜底 ──
        # skipped/unknown 类状态不触发 AI 调用 (续跑/回填时不产生额外成本)
        summary_mode = MODE_TEMPLATE if (mode == MODE_AI
                                         and status in ("skipped", "unknown")) else mode
        markdown, ext_path, source = _try_summaries(
            project, pipeline, stage, mode=summary_mode, state=state,
            output_dir=out_path, config=config, client=client)
        final_summary: Optional[Path] = ext_path
        if markdown is not None:
            default = summary_path(project, pipeline, stage)
            target = ext_path if ext_path is not None else default
            if not target.is_absolute():
                target = project / target
            if not target.exists() or ext_path is None:
                _write_text(target, markdown)
            final_summary = target
        elif ext_path is None:
            markdown = _fallback_markdown(pipeline, stage, state,
                                          stats["n_files"], stats["bytes"],
                                          logp, err_type, err_msg)
            final_summary = summary_path(project, pipeline, stage)
            _write_text(final_summary, markdown)
            source = "fallback:template"

        _append_index(project, {
            "ts": finished,
            "pipeline": pipeline,
            "stage": stage,
            "status": status,
            "summary_path": str(final_summary) if final_summary else None,
            "state_path": str(spath),
            "error": f"{err_type}: {err_msg}" if err_type else None,
        })

        result.update(status="ok", state_path=str(spath),
                      summary_path=str(final_summary) if final_summary else None,
                      summary_source=source)
        return result
    except SystemExit as e:  # 摘要层内部若触发退出, 也只记录不扩散
        result.update(status="error", reason=f"SystemExit: {e}")
        return result
    except Exception as e:  # noqa: BLE001 - 钩子绝不让主管线崩溃
        result.update(status="error", reason=f"{type(e).__name__}: {e}")
        return result


def _guess_log_path(output_dir: Optional[Path], stage: str) -> Optional[Path]:
    """在阶段输出目录里猜日志位置 (调用方给了 log_path 就不会走到这里)。"""
    if output_dir is None:
        return None
    base = Path(output_dir)
    cands = [base / f"{stage}.log", base / "logs" / f"{stage}.log",
             base / f"run_{stage}.log"]
    for c in cands:
        try:
            if c.is_file():
                return c
        except OSError:
            continue
    try:
        for c in sorted(base.glob(f"run_*{stage}*.log"))[:1]:
            return c
    except OSError:
        pass
    return None


# ═══════════════════════════════════════════════════════════════════
# 对外接口 2: backfill
# ═══════════════════════════════════════════════════════════════════

def backfill(project_dir: Any, *, pipeline: Optional[str] = None) -> Dict[str, Any]:
    """为存量项目补写阶段状态/摘要 (按 digest 目录枚举, 不跑任何管线)。

    * 枚举 ``run_digest.PIPELINE_ROOTS`` × ``KEYS_BY_PIPELINE`` 中真实存在
      的 stage 目录; 有文件 → status=ok, 空目录 → status=unknown;
    * 已存在 state 文件默认跳过 (可重复调用, 幂等);
    * 返回 counts (scanned/written/skipped_existing/unknown/errors)。
    """
    mode = resolve_mode()
    out: Dict[str, Any] = {"status": "ok", "mode": mode, "pipeline": pipeline,
                           "scanned": 0, "written": 0,
                           "skipped_existing": 0, "unknown": 0,
                           "errors": 0, "details": []}
    if mode == MODE_OFF:
        out.update(status="skipped", reason=f"{ENV_MODE}=off (钩子关闭)")
        return out
    try:
        project = project_root(project_dir)
        if not project.is_dir():
            out.update(status="error", reason=f"项目目录不存在: {project}")
            return out
        from mmpv_common import io_layout, run_digest  # lazy: 避免顶层耦合

        digest = run_digest.build_digest(project)
        layout = digest.get("layout") or io_layout.resolve_layout(None)
        dirs = io_layout.layout_dirs(layout)
        roots = run_digest.PIPELINE_ROOTS
        keys_by = run_digest.KEYS_BY_PIPELINE

        pipelines = [pipeline] if pipeline else list(keys_by)
        for pk in pipelines:
            if pk not in keys_by:
                out["errors"] += 1
                out["details"].append({"pipeline": pk, "status": "error",
                                       "reason": f"未知 pipeline: {pk}"})
                continue
            if layout == "legacy":
                base = project
            else:
                base = project / roots.get(pk, "")
            if not base.is_dir():
                continue
            for key in keys_by.get(pk, []):
                dirname = dirs.get(key)
                if not dirname:
                    continue
                stage_dir = base / dirname
                if not stage_dir.is_dir():
                    continue
                stage = stage_for_io_key(pk, key)
                out["scanned"] += 1
                spath = state_path(project, pk, stage)
                if spath.is_file():
                    out["skipped_existing"] += 1
                    continue
                stats = scan_output_dir(stage_dir)
                status = "ok" if stats["n_files"] > 0 else "unknown"
                if status == "unknown":
                    out["unknown"] += 1
                res = after_stage(project, pk, stage, status=status,
                                  output_dir=stage_dir)
                if res.get("status") == "ok":
                    out["written"] += 1
                    out["details"].append({"pipeline": pk, "stage": stage,
                                           "status": status,
                                           "state_path": res.get("state_path")})
                else:
                    out["errors"] += 1
                    out["details"].append({"pipeline": pk, "stage": stage,
                                           "status": "error",
                                           "reason": res.get("reason")})
        return out
    except SystemExit as e:  # pragma: no cover - 理论上不会发生
        out.update(status="error", reason=f"SystemExit: {e}")
        return out
    except Exception as e:  # noqa: BLE001
        out.update(status="error", reason=f"{type(e).__name__}: {e}")
        return out


# ═══════════════════════════════════════════════════════════════════
# CLI (调试/回填用; 管线不依赖)
# ═══════════════════════════════════════════════════════════════════

def main(argv: Optional[List[str]] = None) -> int:
    import argparse

    ap = argparse.ArgumentParser(description="MMPV 阶段级 AI 摘要钩子工具")
    ap.add_argument("project", help="项目根目录")
    ap.add_argument("--backfill", action="store_true",
                    help="按现有 stage 目录补写 state/摘要")
    ap.add_argument("--pipeline", default=None, help="仅处理指定管线")
    ap.add_argument("--mode", choices=list(MODES), default=None,
                    help=f"覆盖 {ENV_MODE} (仅本次进程)")
    args = ap.parse_args(argv)
    if args.mode:
        os.environ[ENV_MODE] = args.mode
    if args.backfill:
        report = backfill(Path(args.project), pipeline=args.pipeline)
    else:
        report = {"status": "ok",
                  "hint": "加 --backfill 才写文件; 本工具默认只报告模式",
                  "mode": resolve_mode()}
    print(json.dumps(report, ensure_ascii=False, indent=2))
    return 0 if report.get("status") in ("ok", "skipped") else 1


if __name__ == "__main__":
    raise SystemExit(main())
