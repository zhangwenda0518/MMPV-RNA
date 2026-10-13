#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
stage_tools.py — stages 工具化 (四层递进架构 · L4 编排层)
==========================================================

把 MMPV 的 8 个编排器包装成带 JSON Schema 的「工具」: 每个工具声明

  * inputSchema   合法参数面 (来自 copilot.COMMAND_SPECS, 即真实 argparse);
  * preconditions 前置条件 (run digest 谓词 —— 交接契约的机器可查形态);
  * postconditions 后置条件 (跑完后 digest 里应出现的产物键);
  * exit semantics 退出码语义 (沿 phylo_pipeline 2026-10-03 约定:
    0=干净完成, 1=failed(下游不得继续), 2=blocked(前置缺失, 可跳过))。

plan_chain() 把 advisor 的规则建议升格为按主链 ①→⑧ 排序的执行计划:
每步带命令草稿与显式 missing_inputs (<> 占位符, 大写 ASCII 键) —— **路径不臆造**,
缺的输入列出来等人补 (agent.py 的 --set / 交互提示), 绝不自动填。

本模块是纯规划层: 不执行任何命令 (执行见 agent.py 的受控循环,
确认门与红线见 doc/AI_DRIVEN_DESIGN.md §4)。
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional

try:
    from mmpv_common import stage_cards
    from mmpv_common.advisor import advise
    from mmpv_common.copilot import COMMAND_SPECS
except ImportError:
    from mmpv_common import stage_cards  # type: ignore
    from mmpv_common.advisor import advise  # type: ignore
    from mmpv_common.copilot import COMMAND_SPECS  # type: ignore

# 主链顺序 (README 七管线图; hostref 归 ①, host_expression 为可选补充)
# 【⑥ EVE、⑧ 宿主表达均已归档 2026-10-11 → archive/20261011_{endogenous_virus,host_expression}_pipeline;
#   键保留以兼容历史状态】
PIPELINE_ORDER = ["public", "hostref", "preprocessing", "discovery",
                  "analysis", "phylo", "eve", "submission", "host_expression"]

EXIT_SEMANTICS = {
    0: "干净完成 (继续下一步)",
    1: "failed — 下游不得继续, 停链排查",
    2: "blocked — 前置缺失/依赖未满足, 可跳过该步继续",
}

_PLACEHOLDER = re.compile(r"<([^<>]+)>")

# 前置/后置条件: 工具键 → (digest 谓词描述, 检查函数 / 产物键列表)
# 谓词收 digest (build_digest 产物); 保守判定: 缺证据 = 不满足。
_TOOL_TITLES = {
    "public": "① 公共数据获取",
    "hostref": "① 宿主参考建库",
    "preprocessing": "② 数据清洗 (质控+去宿主)",
    "discovery": "③ 病毒发现",
    "analysis": "④ 已知病毒深度分析",
    "phylo": "⑤ 系统发育/进化",
    "eve": "⑥ EVE 内源病毒筛查",
    "submission": "⑦ 数据提交",
}


def _d(digest: Dict[str, Any], pipeline: str, key: str) -> Any:
    return digest.get("pipelines", {}).get(pipeline, {}).get(key)


def _stage_present(digest: Dict[str, Any], pipeline: str, io_key: str) -> bool:
    return bool(digest.get("pipelines", {}).get(pipeline, {})
                .get("stages", {}).get(io_key, {}).get("exists"))


# ── 前置条件 (机器谓词) ─────────────────────────────────────────────

def _pre_public(d: Dict[str, Any]) -> List[str]:
    return []  # 数据入口, 无前置


def _pre_hostref(d: Dict[str, Any]) -> List[str]:
    return []  # species/taxid 或 --genome-fasta 由参数承担


def _pre_preprocessing(d: Dict[str, Any]) -> List[str]:
    raw = _d(d, "public", "raw_reads") or 0
    sra = _d(d, "public", "raw_sra") or 0
    if not (raw or sra or _stage_present(d, "preprocessing", "d_clean")):
        return ["无原始数据 (① m_down 为空) — 先跑 public 检索下载, 或自备 FASTQ"]
    return []


def _pre_discovery(d: Dict[str, Any]) -> List[str]:
    if not (_stage_present(d, "preprocessing", "d_hostdep")
            or _stage_present(d, "discovery", "d_asm")):
        return ["无去宿主 reads (② 02_HostDepleted 缺失) — 先跑 preprocessing"]
    return []


def _pre_analysis(d: Dict[str, Any]) -> List[str]:
    miss = []
    if not (_stage_present(d, "preprocessing", "d_hostdep")
            or _stage_present(d, "analysis", "a_detect")):
        miss.append("无去宿主 reads — 先跑 preprocessing")
    if not _d(d, "discovery", "present"):
        miss.append("④ 依赖病毒参考库 (ref.fasta + ref_info.tsv), 请确认已备妥"
                    " — 不属于本项目 digest, 工具执行前自查")
    return miss


def _pre_phylo(d: Dict[str, Any]) -> List[str]:
    if not (_stage_present(d, "analysis", "a_extract")
            or _stage_present(d, "phylo", "ph_data")):
        return ["无 ④ 06_Extraction 全长序列 — 先跑 analysis full/extract"]
    return []


def _pre_eve(d: Dict[str, Any]) -> List[str]:
    if not _stage_present(d, "public", "h_genome"):
        return ["无宿主基因组 (① 05_HostRef/genome 缺失) — 先跑 hostref"]
    return []


def _pre_submission(d: Dict[str, Any]) -> List[str]:
    if not _stage_present(d, "discovery", "d_rescue"):
        return ["无 ③ 10_Rescue 产物 — 先跑 discovery 到 rescue"]
    return []


_PRECONDITIONS = {
    "public": _pre_public,
    "hostref": _pre_hostref,
    "preprocessing": _pre_preprocessing,
    "discovery": _pre_discovery,
    "analysis": _pre_analysis,
    "phylo": _pre_phylo,
    "eve": _pre_eve,
    "submission": _pre_submission,
}

# ── 后置条件 (跑完后 digest 应出现的键) ────────────────────────────

_POSTCONDITIONS = {
    "public":       [("public", "raw_reads"), ("public", "core14_rows")],
    "hostref":      [("public", "stages.h_hostdb.exists")],
    "preprocessing": [("preprocessing", "depleted_reads")],
    "discovery":    [("discovery", "stages.d_rescue.exists")],
    "analysis":     [("analysis", "detection_summary_rows")],
    "phylo":        [("phylo", "present")],
    "eve":          [("eve", "kingdom_summary_rows")],
    "submission":   [("discovery", "stages.d_rescue.exists")],  # .sqn 落 rescue 目录旁
}


def check_preconditions(tool_key: str, digest: Dict[str, Any]) -> List[str]:
    """返回未满足的前置条件列表 (空 = 可跑)。"""
    fn = _PRECONDITIONS.get(tool_key)
    if fn is None:
        return [f"未知工具: {tool_key}"]
    try:
        return list(fn(digest))
    except Exception as e:  # 谓词自身异常按未满足处理 (保守)
        return [f"前置检查异常: {e}"]


def check_postconditions(tool_key: str, digest: Dict[str, Any]) -> List[str]:
    """返回未达成的后置条件 (点路径逐段下钻)。仅作运行后核对提示。"""
    out = []
    for pipe, keypath in _POSTCONDITIONS.get(tool_key, []):
        node: Any = digest.get("pipelines", {}).get(pipe, {})
        for seg in keypath.split("."):
            if isinstance(node, dict):
                node = node.get(seg)
            else:
                node = None
                break
        if not node:
            out.append(f"{pipe}.{keypath} 未出现")
    return out


# ═══════════════════════════════════════════════════════════════════
# 工具清单 (MCP tools/list 兼容的 JSON Schema)
# ═══════════════════════════════════════════════════════════════════

def tools_json() -> List[Dict[str, Any]]:
    """8 个管线工具的 {name, description, inputSchema, mmpv: 元数据} 列表。

    注意: 这里只描述参数面, **不提供执行入口** — 工具的调用者是
    agent.py 的受控循环 (人在终端启动), 不是 MCP 远端。
    ⑦ submission 这类带子命令的编排器: stage enum 为空, 另给 subcommand enum。
    """
    out = []
    for key in PIPELINE_ORDER:
        spec = COMMAND_SPECS.get(key)
        if spec is None:
            continue
        card = stage_cards.get_card(spec.key)
        title = _TOOL_TITLES.get(key, spec.key)
        allowed = list(spec.optional)
        for sub in spec.subcommands:
            allowed += [f"{sub.name}:{f}" for f in sub.optional]
        props: Dict[str, Any] = {
            "project_dir": {"type": "string", "description": "项目根目录"},
            "stage": {"type": "string", "enum": list(spec.stages),
                      "description": "要跑的阶段 (可多选语义见编排器)"},
            "args": {"type": "object",
                     "description": (f"额外 flag→值; 允许: {', '.join(allowed)}"
                                     f"; 必填: {', '.join(spec.required)}")},
        }
        if spec.subcommands:
            props["subcommand"] = {
                "type": "string", "enum": [s.name for s in spec.subcommands],
                "description": "顶层子命令 (必选其一; 该编排器无 --stage)"}
        out.append({
            "name": f"mmpv_run_{key}",
            "title": title,
            "description": (f"{title} — {spec.desc}。前置: {'; '.join(check_preconditions(key, {})) or '无 (入口工具)'}"
                            f"。退出码: 0=完成 1=failed停链 2=blocked可跳。"
                            f"本工具仅描述参数面, 执行须经 agent.py 受控循环 (人工确认门)。"),
            "inputSchema": {"type": "object", "properties": props},
            "mmpv": {
                "script": spec.script,
                "stages": list(spec.stages),
                "subcommands": [s.name for s in spec.subcommands],
                "required_args": list(spec.required),
                "dangerous_args": list(spec.dangerous),
                "precondition_keys": list(_POSTCONDITIONS.get(key, [])),
                "example": spec.hint,
            },
        })
    return out


# ═══════════════════════════════════════════════════════════════════
# 计划器: advisor 建议 → 主链排序的可执行计划
# ═══════════════════════════════════════════════════════════════════

@dataclass
class ChainStep:
    order: int
    tool: str                      # COMMAND_SPECS 键
    pipeline: str                  # stage_cards 管线键
    title: str
    rationale: str                 # 来自 advisor 规则 (title + detail)
    command: str                   # 占位符替换后的命令草稿
    missing_inputs: List[str] = field(default_factory=list)   # 未填 <> 占位符
    unmet_preconditions: List[str] = field(default_factory=list)

    def to_dict(self) -> Dict[str, Any]:
        return {"order": self.order, "tool": self.tool, "pipeline": self.pipeline,
                "title": self.title, "rationale": self.rationale,
                "command": self.command, "missing_inputs": self.missing_inputs,
                "unmet_preconditions": self.unmet_preconditions}


def _map_rec_pipeline_to_tool(pipeline: str) -> str:
    """advisor 用 stage_cards 管线键; 工具键基本同名, public 分检索/建库两个。"""
    return pipeline  # public 建库建议不出现在 advisor (无对应规则), 直通即可


def substitute_placeholders(command: str, project_dir: str,
                            inputs: Optional[Dict[str, str]] = None) -> tuple:
    """<PROJECT_ROOT> 自动替换为 project_dir; 其余占位符用 inputs 补, 补不上的返回 missing。"""
    inputs = inputs or {}
    missing: List[str] = []

    def _sub(m: "re.Match[str]") -> str:
        name = m.group(1)
        if name == "PROJECT_ROOT":
            return project_dir
        if name in inputs and inputs[name]:
            return inputs[name]
        missing.append(name)
        return m.group(0)

    return _PLACEHOLDER.sub(_sub, command), missing


def plan_chain(digest: Dict[str, Any], *,
               inputs: Optional[Dict[str, str]] = None,
               pipelines: Optional[List[str]] = None) -> List[ChainStep]:
    """advisor 规则建议 → 按主链排序的计划。纯函数, 不执行。"""
    project_dir = digest.get("project_dir", "")
    recs = advise(digest)
    steps: List[ChainStep] = []
    for rec in recs:
        tool = _map_rec_pipeline_to_tool(rec.pipeline)
        if tool not in COMMAND_SPECS or tool not in PIPELINE_ORDER:
            continue
        if pipelines and tool not in pipelines:
            continue
        command, missing = substitute_placeholders(rec.command, project_dir, inputs)
        unmet = check_preconditions(tool, digest)
        steps.append(ChainStep(
            order=PIPELINE_ORDER.index(tool), tool=tool, pipeline=rec.pipeline,
            title=rec.title, rationale=f"{rec.title}: {rec.detail}",
            command=command, missing_inputs=missing, unmet_preconditions=unmet))
    steps.sort(key=lambda s: s.order)
    return steps
