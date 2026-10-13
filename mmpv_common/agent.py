#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
agent.py — Agent 全链路受控编排 (四层递进架构 · L4 编排层)
===========================================================

落地顺序第四阶段的最后一环: plan → (确认门) → execute → read_digest →
iterate 的**受控**循环。与"全自动 Agent"的本质区别 (红线, 见
doc/AI_DRIVEN_DESIGN.md §4):

  * 计划只来自规则 (advisor + stage_tools.plan_chain), AI 不决定跑什么;
  * 每步命令先过 copilot 白名单三道防线, 以 token 列表直接 subprocess
    (不经 shell, 无注入面);
  * 危险命令 (--force 类) 即使 --yes 也必须交互确认; 非交互一律跳过;
  * 路径不臆造: 计划里的 <> 占位符缺输入就跳过该步并明说;
  * 预算封顶: --max-steps; 退出码语义: 1=failed 停链, 2=blocked 跳过;
  * 全程留痕: 90_Handoff/state/agent_runs/<ts>.json + ai_log.jsonl。

用法 (经 mmpv_ai.py):
  python mmpv_ai.py agent <项目根> --dry-run                  # 只看计划
  python mmpv_ai.py agent <项目根> --set "RAW_FASTQ_DIR=D:/raw" --yes
  python mmpv_ai.py agent <项目根> --only preprocessing,discovery
"""

from __future__ import annotations

import json
import subprocess
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional

try:
    from mmpv_common import run_digest
    from mmpv_common.ai import log_ai_event
    from mmpv_common.copilot import validate
    from mmpv_common.stage_tools import (PIPELINE_ORDER, ChainStep,
                                         check_postconditions, plan_chain)
except ImportError:
    from mmpv_common import run_digest  # type: ignore
    from mmpv_common.ai import log_ai_event  # type: ignore
    from mmpv_common.copilot import validate  # type: ignore
    from mmpv_common.stage_tools import (PIPELINE_ORDER, ChainStep,  # type: ignore
                                         check_postconditions, plan_chain)

# 步骤状态
ST_PLANNED = "planned"            # dry-run: 只计划未执行
ST_RAN = "ran"                    # 退出码 0
ST_FAILED = "failed"              # 退出码 1 → 停链
ST_BLOCKED = "blocked"            # 退出码 2 → 跳过继续
ST_SKIPPED = "skipped"            # 缺输入/白名单不过/前置不满足
ST_UNCONFIRMED = "unconfirmed"    # 危险命令未获确认

STATUS_RANK = {ST_RAN: 0, ST_PLANNED: 1, ST_BLOCKED: 2, ST_SKIPPED: 3,
               ST_UNCONFIRMED: 4, ST_FAILED: 5}


@dataclass
class AgentStep:
    index: int
    tool: str
    title: str
    command: str
    status: str
    detail: str = ""
    exit_code: Optional[int] = None
    postcondition_missing: List[str] = field(default_factory=list)

    def to_dict(self) -> Dict[str, Any]:
        return {"index": self.index, "tool": self.tool, "title": self.title,
                "command": self.command, "status": self.status,
                "detail": self.detail, "exit_code": self.exit_code,
                "postcondition_missing": self.postcondition_missing}


@dataclass
class AgentReport:
    project_dir: str
    steps: List[AgentStep] = field(default_factory=list)
    started_at: str = ""
    finished_at: str = ""
    journal_path: str = ""
    stop_reason: str = ""          # "" | "failed" | "max_steps" | "计划为空"

    def summary_line(self) -> str:
        ran = sum(1 for s in self.steps if s.status == ST_RAN)
        failed = sum(1 for s in self.steps if s.status == ST_FAILED)
        skipped = sum(1 for s in self.steps
                      if s.status in (ST_SKIPPED, ST_UNCONFIRMED, ST_BLOCKED))
        return f"ran {ran} · failed {failed} · skipped {skipped} · {len(self.steps)} step(s)" \
               + (f" · stop: {self.stop_reason}" if self.stop_reason else "")

    def to_dict(self) -> Dict[str, Any]:
        return {"project_dir": self.project_dir, "started_at": self.started_at,
                "finished_at": self.finished_at, "journal_path": self.journal_path,
                "stop_reason": self.stop_reason, "summary": self.summary_line(),
                "steps": [s.to_dict() for s in self.steps]}


# ═══════════════════════════════════════════════════════════════════

def _default_confirm(prompt: str) -> bool:
    try:
        return input(prompt).strip().lower() in ("y", "yes")
    except (EOFError, KeyboardInterrupt):
        return False


def run_agent(project_dir: Path, *,
              inputs: Optional[Dict[str, str]] = None,
              pipelines: Optional[List[str]] = None,
              max_steps: int = 8,
              dry_run: bool = False,
              yes: bool = False,
              confirm_fn: Optional[Callable[[str], bool]] = None,
              executor: Optional[Callable[..., Any]] = None) -> AgentReport:
    """受控 agent 循环。executor/confirm_fn 可注入 (测试桩)。

    语义约定:
      * dry_run=True  → 全部步骤只产出计划 (status=planned), 不执行;
      * yes=True      → 非危险步骤免确认; 危险步骤仍须 confirm_fn 同意,
                        无 confirm_fn (非交互) 则跳过 (红线: 危险不自动);
      * yes=False     → 有 confirm_fn 则逐步确认, 没有则只读计划
                        (等价 dry_run 的执行面)。
    """
    project_dir = Path(project_dir)
    inputs = dict(inputs or {})
    executor = executor or (lambda cmd, **kw: subprocess.run(
        cmd, capture_output=True, text=True, **kw))
    confirm = confirm_fn or _default_confirm

    started = time.strftime("%Y-%m-%d %H:%M:%S")
    digest = run_digest.build_digest(project_dir)
    chain = plan_chain(digest, inputs=inputs, pipelines=pipelines)
    report = AgentReport(project_dir=str(project_dir), started_at=started)

    if not chain:
        report.stop_reason = ("empty plan (advisor has no recommendation — project may be "
                              "fully advanced or missing raw data)")
        _finish(report, project_dir, started)
        return report

    ran_steps = 0
    for i, step in enumerate(chain):
        if ran_steps >= max_steps:
            report.stop_reason = f"max_steps={max_steps} budget exhausted"
            break
        idx = len(report.steps)
        # 缺输入 → 明说跳过 (路径不臆造)
        if step.missing_inputs:
            report.steps.append(AgentStep(
                idx, step.tool, step.title, step.command, ST_SKIPPED,
                detail="missing inputs " + ", ".join(f"<{m}>" for m in step.missing_inputs)
                       + " — fill with --set \"PLACEHOLDER=value\" and rerun"))
            continue
        # 白名单三道防线
        v = validate(step.command)
        if not v.ok:
            report.steps.append(AgentStep(
                idx, step.tool, step.title, step.command, ST_SKIPPED,
                detail="白名单校验未过: " + "; ".join(v.errors)))
            continue
        # 前置条件 (以最新 digest 复查)
        fresh = run_digest.build_digest(project_dir)
        unmet = [c for c in step.unmet_preconditions] or []
        if not unmet:
            from mmpv_common.stage_tools import check_preconditions
            unmet = check_preconditions(step.tool, fresh)
        if unmet:
            report.steps.append(AgentStep(
                idx, step.tool, step.title, v.preview.splitlines()[0] if v.preview else step.command,
                ST_BLOCKED, detail="前置未满足: " + "; ".join(unmet)))
            continue
        # 确认门 (危险命令红线)
        if not dry_run:
            need_confirm = v.requires_confirmation or not yes
            if need_confirm:
                if confirm_fn is None and v.requires_confirmation:
                    report.steps.append(AgentStep(
                        idx, step.tool, step.title, step.command, ST_UNCONFIRMED,
                        detail="危险命令且非交互 (--yes 不豁免), 已跳过"))
                    continue
                if not confirm(f"[agent] 第 {idx + 1} 步 {step.title}\n  {step.command}\n  确认执行? [y/N]: "):
                    report.steps.append(AgentStep(
                        idx, step.tool, step.title, step.command, ST_UNCONFIRMED,
                        detail="人工未确认"))
                    continue
        else:
            report.steps.append(AgentStep(
                idx, step.tool, step.title, step.command, ST_PLANNED,
                detail="[dry-run] " + (v.preview.splitlines()[0] if v.preview else "")))
            ran_steps += 1
            continue

        # ── 执行 (token 列表直 subprocess, 不经 shell) ──
        try:
            tokens = _tokens_from_validation(step.command)
            proc = executor(tokens)
            code = getattr(proc, "returncode", 1)
        except Exception as e:
            report.steps.append(AgentStep(
                idx, step.tool, step.title, step.command, ST_FAILED,
                detail=f"执行异常: {type(e).__name__}: {e}", exit_code=-1))
            report.stop_reason = "failed"
            break
        if code == 0:
            post_miss = check_postconditions(step.tool,
                                             run_digest.build_digest(project_dir))
            report.steps.append(AgentStep(
                idx, step.tool, step.title, step.command, ST_RAN,
                detail="完成" + ("" if not post_miss else
                                 f" (注意后置未达: {'; '.join(post_miss)})"),
                exit_code=0, postcondition_missing=post_miss))
            ran_steps += 1
        elif code == 2:
            report.steps.append(AgentStep(
                idx, step.tool, step.title, step.command, ST_BLOCKED,
                detail="编排器报 blocked (前置缺失), 跳过继续", exit_code=2))
            ran_steps += 1
        else:
            report.steps.append(AgentStep(
                idx, step.tool, step.title, step.command, ST_FAILED,
                detail="退出码 1 — 按停链语义终止, 排查后再续跑", exit_code=int(code)))
            report.stop_reason = "failed"
            break

    _finish(report, project_dir, started)
    return report


def _tokens_from_validation(command: str) -> List[str]:
    """从白名单校验后的原始命令取 token 列表 (与 validate 同一套切词)。"""
    from mmpv_common.copilot import parse_tokens
    return parse_tokens(command)


def _finish(report: AgentReport, project_dir: Path, started: str) -> None:
    report.finished_at = time.strftime("%Y-%m-%d %H:%M:%S")
    state = run_digest.state_dir(project_dir)
    try:
        journal_dir = state / "agent_runs"
        journal_dir.mkdir(parents=True, exist_ok=True)
        path = journal_dir / (time.strftime("%Y%m%d_%H%M%S") + ".json")
        path.write_text(json.dumps(report.to_dict(), ensure_ascii=False, indent=2),
                        encoding="utf-8")
        report.journal_path = str(path)
        log_ai_event(state, {"kind": "agent", "op": "run_chain",
                             "summary": report.summary_line(),
                             "stop_reason": report.stop_reason,
                             "n_steps": len(report.steps),
                             "journal": str(path)})
    except OSError:
        pass  # 留痕失败不吞执行结果
