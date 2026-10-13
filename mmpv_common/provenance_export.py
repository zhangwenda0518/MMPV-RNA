#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
provenance_export.py — AI 参与溯源一键导出 (P3)
================================================

论文方法章的"AI 参与声明"一键生成器。聚合项目里全部 AI 决策留痕
(ai_log.jsonl / agent_runs/ / digest_history.jsonl / *.triage.tsv) +
平台元数据 (mmpv_common 版本, 布局, 生成时间), 产出两份:

  90_Handoff/state/provenance_export.json   机器可读全量记录
  90_Handoff/AI_PROVENANCE.md               人读方法学文档 (含中英文
                                            可直接粘贴的方法章段落)

红线不变: 这里只**汇总**已发生的留痕, 不产生任何新的 AI 调用 ——
导出本身就是确定性的。

CLI: python mmpv_ai.py provenance <项目根>
"""

from __future__ import annotations

import json
import time
from pathlib import Path
from typing import Any, Dict, List, Optional

try:
    import mmpv_common
    from mmpv_common import process_report, run_digest
except ImportError:
    import sys
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
    import mmpv_common  # type: ignore
    from mmpv_common import process_report, run_digest  # type: ignore

EXPORT_JSON = "provenance_export.json"
EXPORT_MD = "AI_PROVENANCE.md"
SCHEMA_VERSION = 1


def build_provenance(project_dir: Path, layout: Optional[str] = None) -> Dict[str, Any]:
    """聚合 → 结构化溯源记录 (纯确定性, 零 AI 调用)。"""
    project = Path(project_dir)
    record = process_report.collect_process(project, layout=layout)
    state = run_digest.state_dir(project)
    ai_log_events: List[Dict[str, Any]] = []
    log_path = state / "ai_log.jsonl"
    if log_path.is_file():
        for line in log_path.read_text(encoding="utf-8").splitlines():
            line = line.strip()
            if not line:
                continue
            try:
                ai_log_events.append(json.loads(line))
            except json.JSONDecodeError:
                continue
    return {
        "schema_version": SCHEMA_VERSION,
        "exported_at": time.strftime("%Y-%m-%d %H:%M:%S"),
        "project_dir": str(project),
        "layout": record.get("layout"),
        "platform": {"name": "MMPV (Massive Meta-mining of Plant Viruses)",
                     "ai_layer_version": mmpv_common.__version__,
                     "design_doc": "doc/AI_DRIVEN_DESIGN.md"},
        "digest_totals": record["digest"].get("totals", {}),
        "ai_summary": record["ai"],
        "ai_log_events": ai_log_events,
        "agent_runs": record["agent_runs"],
        "timeline": record["timeline"],
        "triage": record["triage"],
        "governance": (
            "AI 仅用于结果解读、参数/流程建议、候选二值标注 (KEEP/REVIEW) 与命令白名单校验; "
            "不执行管线命令, 不改写任何统计算法结论 (分类=多工具加权投票, 完整性=CheckV, "
            "差异=DESeq2/TEA-GCN); 全部 AI 调用含模型、提示词版本与输入哈希, 留痕于 "
            "90_Handoff/state/ai_log.jsonl, 可离线复现规则部分。"),
    }


def _md_table(rows: List[List[Any]], header: List[str]) -> str:
    out = ["| " + " | ".join(header) + " |",
           "|" + "|".join("---" for _ in header) + "|"]
    for r in rows:
        out.append("| " + " | ".join(str(c) for c in r) + " |")
    return "\n".join(out)


def build_markdown(prov: Dict[str, Any]) -> str:
    ai = prov.get("ai_summary", {})
    lines: List[str] = [
        "# MMPV AI 参与溯源 (AI Provenance)",
        "",
        f"> 导出: {prov['exported_at']} · 项目: `{prov['project_dir']}` · "
        f"布局: {prov.get('layout')} · AI 层版本: mmpv_common "
        f"{prov['platform']['ai_layer_version']}",
        "",
        "本文件由 `python mmpv_ai.py provenance <项目根>` 自动生成 (确定性聚合, "
        "导出过程零 AI 调用)。机器可读全量记录: "
        "`90_Handoff/state/provenance_export.json`。",
        "",
        "## 1. 项目推进",
        "",
        f"管线在线: **{prov['digest_totals'].get('pipelines_present', 0)}** 条; "
        f"digest 时间线快照 {len(prov.get('timeline', []))} 条。",
        "",
        "## 2. AI 参与总览",
        "",
    ]
    if ai.get("n_calls"):
        lines.append(_md_table(
            [[op, n, ai.get("by_op_tokens", {}).get(op, 0)]
             for op, n in ai["by_op"].items()],
            ["操作 (op)", "调用次数", "tokens"]))
        lines += ["",
                  f"- 模型: {', '.join(ai['models'])}",
                  f"- 提示词版本: {', '.join(ai['prompt_versions'])}",
                  f"- token 合计: ~{ai['total_tokens']:,}"
                  + (f" (缓存命中 {ai['n_cached']} 次, 未重复计费)" if ai.get("n_cached") else ""),
                  f"- 首次调用: {ai['first_ts']} · 最近调用: {ai['last_ts']}"]
    else:
        lines.append("**无 AI 调用记录 — 本项目为纯规则运行。**")
    agents = prov.get("agent_runs", [])
    lines += ["", "## 3. 受控 agent 运行", ""]
    if agents:
        lines.append(_md_table(
            [[a.get("finished_at", ""), a.get("summary", ""), a.get("stop_reason", "") or "正常收尾"]
             for a in agents], ["结束时间", "结局", "停因"]))
    else:
        lines.append("无 agent 运行记录。")
    triage = prov.get("triage", [])
    lines += ["", "## 4. 候选 triage 标注", ""]
    if triage:
        lines.append(_md_table(
            [[t["file"], t["rows"], t["keep"], t["review"]] for t in triage],
            ["标注文件", "行数", "KEEP", "REVIEW"]))
        lines.append("")
        lines.append("标注只追加 `ai_action/ai_reason/ai_source` 三列; "
                     "原始数据列逐字节未改; REVIEW 为待人工裁决项。")
    else:
        lines.append("无 triage 标注文件。")
    lines += ["", "## 5. 治理声明 (红线)", "", prov["governance"], ""]

    # 方法学段落 (可直接粘贴论文)
    if ai.get("n_calls"):
        zh = (f"数据分析由 MMPV 平台完成, 其中大语言模型 (LLM) 仅作为辅助层参与: "
              f"共调用 {ai['n_calls']} 次 (模型: {', '.join(ai['models'])}), 用途限于结果解读、"
              f"流程建议与候选序列的二值标注 (KEEP/REVIEW, 不改变任何计算结果); 所有病毒分类"
              f"判定均由多工具加权投票等确定性算法给出, 基因组完整性由 CheckV 评估, "
              f"差异表达由 DESeq2 与 TEA-GCN 完成。全部 LLM 调用的模型、提示词版本与"
              f"输入摘要哈希均已留痕 ({'90_Handoff/state/ai_log.jsonl'}), 涉及 AI 标注的"
              f"候选均经人工复核 REVIEW 项后采用。")
        en = (f"Data analysis was performed with the MMPV platform, in which large "
              f"language models served only as an assistive layer: {ai['n_calls']} calls "
              f"(models: {', '.join(ai['models'])}) were made for result interpretation, "
              f"workflow advice, and binary candidate labelling (KEEP/REVIEW) that never "
              f"altered any computational result. All taxonomic assignments were produced "
              f"by deterministic weighted voting across identification tools; genome "
              f"completeness was assessed with CheckV; differential expression was "
              f"computed with DESeq2 and TEA-GCN. Every LLM call (model, prompt version, "
              f"input digest) is logged in 90_Handoff/state/ai_log.jsonl, and all "
              f"AI-labelled candidates underwent manual review of REVIEW items before use.")
    else:
        zh = ("数据分析由 MMPV 平台完成, 全程为确定性规则计算, 未使用大语言模型; "
              "病毒分类由多工具加权投票给出, 完整性由 CheckV 评估, 差异表达由 "
              "DESeq2 与 TEA-GCN 完成。")
        en = ("Data analysis was performed with the MMPV platform entirely by "
              "deterministic rule-based computation, with no large language model "
              "involved; taxonomic assignments were produced by weighted voting across "
              "identification tools, completeness was assessed with CheckV, and "
              "differential expression was computed with DESeq2 and TEA-GCN.")
    lines += ["## 6. 方法学段落 (可直接引用)", "",
              "**中文:**", "", zh, "", "**English:**", "", en, ""]
    return "\n".join(lines)


def export_provenance(project_dir: Path, layout: Optional[str] = None,
                      md_out: Optional[Path] = None) -> Dict[str, Any]:
    """生成 JSON + MD 两份溯源文件, 返回路径与摘要。"""
    project = Path(project_dir)
    prov = build_provenance(project, layout)
    state = run_digest.state_dir(project)
    state.mkdir(parents=True, exist_ok=True)
    jpath = state / EXPORT_JSON
    jpath.write_text(json.dumps(prov, ensure_ascii=False, indent=2), encoding="utf-8")
    mpath = Path(md_out) if md_out else project / run_digest.HANDOFF_DIR / EXPORT_MD
    mpath.parent.mkdir(parents=True, exist_ok=True)
    mpath.write_text(build_markdown(prov), encoding="utf-8")
    return {"status": "ok", "json": str(jpath), "md": str(mpath),
            "n_ai_calls": prov["ai_summary"].get("n_calls", 0),
            "n_agent_runs": len(prov["agent_runs"]),
            "n_triage_files": len(prov["triage"])}
