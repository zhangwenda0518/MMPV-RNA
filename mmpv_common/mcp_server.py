#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
mcp_server.py — MMPV MCP 服务器 (四层递进架构 · L4 编排层)
===========================================================

把 MMPV 的 AI 能力 (L1 digest / L3 顾问 / L1 卡片 / L4 命令校验 / L2 问答)
以 MCP (Model Context Protocol) 工具的形式开放给任意 AI 客户端
(ZCode / Claude / Cursor / 任何支持 MCP 的 IDE)。

零第三方依赖: 手写 stdio 上的 JSON-RPC 2.0 (MCP 规范的 newline-delimited
JSON 传输), 只 import 本包模块 —— 与 io_layout 同级的"可安全 import"纪律。

【红线】刻意【不提供】任何执行类工具: 没有 run_pipeline / run_command。
AI 客户端最多拿到校验过的命令草稿 (mmpv_validate_command), 执行权在
人的终端 (doc/AI_DRIVEN_DESIGN.md §红线 4)。

接入 (以 ZCode 为例, ~/.zcode/workspace/<ws>/.zcode/config.json):
  "mcp": {"servers": {"mmpv": {"type": "stdio", "command": "python",
          "args": ["<仓库根>/mmpv_common/mcp_server.py"]}}}

工具清单:
  mmpv_project_digest     项目 run digest (compact 建议 true)
  mmpv_advise             规则顾问 (+ 可选 AI 叙事, 需环境变量已配置)
  mmpv_stage_cards        阶段卡片查询
  mmpv_validate_command   命令白名单校验 + dry-run 预览 (不执行)
  mmpv_ask                带 digest 上下文的 AI 问答 (离线时明确告知)
"""

from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from mmpv_common import stage_cards                        # noqa: E402
from mmpv_common import run_digest                         # noqa: E402
from mmpv_common import stage_tools                        # noqa: E402
from mmpv_common.advisor import advise, advise_with_ai     # noqa: E402
from mmpv_common.copilot import validate                   # noqa: E402
from mmpv_common.ai import ai_status, ask                  # noqa: E402

SERVER_NAME = "mmpv"
SERVER_VERSION = "0.1.0"
PROTOCOL_VERSION = "2024-11-05"


# ═══════════════════════════════════════════════════════════════════
# 工具实现 (输入 dict → 文本)
# ═══════════════════════════════════════════════════════════════════

def _tool_digest(args: Dict[str, Any]) -> str:
    project = Path(args["project_dir"])
    digest = run_digest.build_digest(project)      # 只构建一次: 落盘与返回共用同一份
    # write 默认关, 且仅显式 true 才落盘 (字符串 "true"/1 等一律不触发写文件)
    out = run_digest.write_digest(digest, project) if args.get("write") is True else None
    if args.get("compact", True):
        digest = run_digest.compact_for_ai(digest)
    if out is not None:
        digest["written_to"] = str(out)
    return json.dumps(digest, ensure_ascii=False, indent=2)


def _tool_advise(args: Dict[str, Any]) -> str:
    digest = run_digest.build_digest(Path(args["project_dir"]))
    if args.get("ai_narrative", False):
        payload = advise_with_ai(digest)
        return json.dumps(payload, ensure_ascii=False, indent=2)
    recs = advise(digest)
    return json.dumps({"recommendations": [r.to_dict() for r in recs]},
                      ensure_ascii=False, indent=2)


def _tool_cards(args: Dict[str, Any]) -> str:
    pipeline = args.get("pipeline")
    stage = args.get("stage")
    if pipeline:
        card = stage_cards.get_card(pipeline, stage)
        if card is None:
            return json.dumps({"error": f"未知 pipeline/stage: {pipeline}/{stage}",
                               "pipelines": stage_cards.pipelines()},
                              ensure_ascii=False)
        return json.dumps(card, ensure_ascii=False, indent=2)
    return stage_cards.summary_text()


def _tool_validate(args: Dict[str, Any]) -> str:
    result = validate(str(args.get("command", "")))
    body = result.to_dict()
    body["note"] = "本工具只校验不执行; 请把 preview 中的完整命令复制到终端人工执行"
    return json.dumps(body, ensure_ascii=False, indent=2)


def _tool_stage_tools(args: Dict[str, Any]) -> str:
    """stages 工具化清单: 参数面/前置/后置/退出码语义 (只描述, 不执行)。"""
    tools = stage_tools.tools_json()
    key = args.get("tool")
    if key:
        tools = [t for t in tools if t["name"] == key or t["name"] == f"mmpv_run_{key}"]
        if not tools:
            return json.dumps({"error": f"未知工具: {key}"}, ensure_ascii=False)
    return json.dumps({"tools": tools,
                       "exit_semantics": stage_tools.EXIT_SEMANTICS,
                       "note": "执行须经 agent.py 受控循环 (人在终端, 人工确认门)"},
                      ensure_ascii=False, indent=2)


def _tool_plan_chain(args: Dict[str, Any]) -> str:
    """基于 digest 的主链执行计划 (只规划, 不执行)。"""
    digest = run_digest.build_digest(Path(args["project_dir"]))
    inputs = args.get("inputs") or {}
    pipelines = args.get("pipelines")
    steps = stage_tools.plan_chain(digest, inputs=inputs, pipelines=pipelines)
    return json.dumps({"project_dir": str(args["project_dir"]),
                       "n_steps": len(steps),
                       "steps": [s.to_dict() for s in steps],
                       "note": "missing_inputs 非空的步骤需要人工补参 (--set), 路径不臆造"},
                      ensure_ascii=False, indent=2)


def _tool_ask(args: Dict[str, Any]) -> str:
    question = str(args.get("question", "")).strip()
    if not question:
        return json.dumps({"error": "question 不能为空"}, ensure_ascii=False)
    context = ""
    project = args.get("project_dir")
    if project:
        slim = run_digest.compact_for_ai(run_digest.build_digest(Path(project)))
        context = "项目 run digest:\n" + json.dumps(slim, ensure_ascii=False)
    ok, why = ai_status()
    if not ok and not args.get("_client"):
        return json.dumps({"status": "unavailable", "reason": why,
                           "hint": "export MMPV_AI_API_KEY / DEEPSEEK_API_KEY 后重启 MCP"},                           ensure_ascii=False)
    res = ask(question, context=context, client=args.get("_client"))
    return json.dumps(res, ensure_ascii=False, indent=2)


TOOLS: List[Dict[str, Any]] = [
    {"name": "mmpv_project_digest",
     "description": "提取 MMPV 项目 run digest: 各管线跑到哪、关键计数。AI 只读此摘要, 不读原始数据。",
     "inputSchema": {"type": "object", "properties": {
         "project_dir": {"type": "string", "description": "MMPV 项目根目录"},
         "compact": {"type": "boolean", "default": True},
         "write": {"type": "boolean", "default": False}},
         "required": ["project_dir"]}},
    {"name": "mmpv_advise",
     "description": "基于 digest 的下一步建议 (规则引擎; ai_narrative=true 时附 AI 总体叙事)。",
     "inputSchema": {"type": "object", "properties": {
         "project_dir": {"type": "string"},
         "ai_narrative": {"type": "boolean", "default": False}},
         "required": ["project_dir"]}},
    {"name": "mmpv_stage_cards",
     "description": "MMPV 阶段卡片: 各管线各 stage 的输入/输出/风险/前置。不传参数返回全平台总览。",
     "inputSchema": {"type": "object", "properties": {
         "pipeline": {"type": "string", "enum": stage_cards.pipelines()},
         "stage": {"type": "string"}}}},
    {"name": "mmpv_validate_command",
     "description": "把一条 MMPV 命令草稿过白名单校验并给 dry-run 预览。只校验, 绝不执行。",
     "inputSchema": {"type": "object", "properties": {
         "command": {"type": "string"}}, "required": ["command"]}},
    {"name": "mmpv_stage_tools",
     "description": "stages 工具化清单: 8 个编排器工具的参数面/前置条件/后置条件/退出码语义。只描述, 不执行。",
     "inputSchema": {"type": "object", "properties": {
         "tool": {"type": "string", "description": "可选, 过滤单个工具 (如 discovery)"}}}},
    {"name": "mmpv_plan_chain",
     "description": "基于项目 digest 生成按主链 ①→⑧ 排序的执行计划 (advisor 规则驱动, 不执行)。",
     "inputSchema": {"type": "object", "properties": {
         "project_dir": {"type": "string"},
         "pipelines": {"type": "array", "items": {"type": "string"}},
         "inputs": {"type": "object", "description": "占位符→值的替换表, 如 {\"RAW_FASTQ_DIR\": \"D:/raw\"}"}},
         "required": ["project_dir"]}},
    {"name": "mmpv_ask",
     "description": "带项目 digest 上下文的 AI 问答 (需已配置 MMPV_AI_API_KEY/DEEPSEEK_API_KEY)。",
     "inputSchema": {"type": "object", "properties": {
         "question": {"type": "string"},
         "project_dir": {"type": "string"}},
         "required": ["question"]}},
]

_DISPATCH = {
    "mmpv_project_digest": _tool_digest,
    "mmpv_advise": _tool_advise,
    "mmpv_stage_cards": _tool_cards,
    "mmpv_validate_command": _tool_validate,
    "mmpv_stage_tools": _tool_stage_tools,
    "mmpv_plan_chain": _tool_plan_chain,
    "mmpv_ask": _tool_ask,
}


# ═══════════════════════════════════════════════════════════════════
# JSON-RPC / MCP 骨架
# ═══════════════════════════════════════════════════════════════════

def _resp(rid: Any, result: Dict[str, Any]) -> Dict[str, Any]:
    return {"jsonrpc": "2.0", "id": rid, "result": result}


def _err(rid: Any, code: int, message: str) -> Dict[str, Any]:
    return {"jsonrpc": "2.0", "id": rid, "error": {"code": code, "message": message}}


def handle(msg: Dict[str, Any]) -> Optional[Dict[str, Any]]:
    """处理一条 JSON-RPC 消息。通知 (无 id) 返回 None。"""
    method = msg.get("method", "")
    rid = msg.get("id")
    params = msg.get("params") or {}

    if method == "initialize":
        return _resp(rid, {
            "protocolVersion": params.get("protocolVersion") or PROTOCOL_VERSION,
            "capabilities": {"tools": {}},
            "serverInfo": {"name": SERVER_NAME, "version": SERVER_VERSION},
        })
    if method == "ping":
        return _resp(rid, {})
    if method == "tools/list":
        return _resp(rid, {"tools": TOOLS})
    if method == "tools/call":
        name = params.get("name", "")
        fn = _DISPATCH.get(name)
        if fn is None:
            return _err(rid, -32602, f"未知工具: {name}")
        try:
            text = fn(params.get("arguments") or {})
            return _resp(rid, {"content": [{"type": "text", "text": text}]})
        except Exception as e:
            return _resp(rid, {"content": [{"type": "text",
                                            "text": f"{type(e).__name__}: {e}"}],
                               "isError": True})
    if method.startswith("notifications/"):
        return None
    if rid is not None:
        return _err(rid, -32601, f"method not found: {method}")
    return None


def serve(stdin=None, stdout=None) -> int:
    """stdio 主循环。Ctrl+C / EOF 退出。日志走 stderr。"""
    stdin = stdin or sys.stdin
    stdout = stdout or sys.stdout
    print(f"[mmpv-mcp] {SERVER_NAME} v{SERVER_VERSION} ready "
          f"({len(TOOLS)} tools)", file=sys.stderr)
    for line in stdin:
        line = line.strip()
        if not line:
            continue
        try:
            msg = json.loads(line)
        except json.JSONDecodeError as e:
            out = _err(None, -32700, f"parse error: {e}")
        else:
            out = handle(msg)
        if out is not None:
            stdout.write(json.dumps(out, ensure_ascii=False) + "\n")
            stdout.flush()
    return 0


if __name__ == "__main__":
    sys.exit(serve())
