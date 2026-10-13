#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
advisor.py — MMPV 参数/流程顾问 (四层递进架构 · L3 决策层)
===========================================================

输入 L1 的 run digest, 输出结构化建议列表。两级:

  * advise()          纯规则: 确定性、零成本、离线可用 —— 一切建议的骨架;
  * advise_with_ai()  在规则骨架之上让 LLM 写总体叙事 (只准串规则建议,
                      不准新增建议; 数字只准引用 digest) —— 补"为什么"。

建议统一是「建议」: 落 execution("command") 供人/copilot 复核,
本模块绝不执行任何命令 (红线, 见 doc/AI_DRIVEN_DESIGN.md §红线)。
"""

from __future__ import annotations

import json
import statistics
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional

try:
    from mmpv_common import stage_cards
    from mmpv_common.ai import render_prompt, chat, log_ai_event
except ImportError:
    from mmpv_common import stage_cards  # type: ignore
    from mmpv_common.ai import render_prompt, chat, log_ai_event  # type: ignore

# ═══════════════════════════════════════════════════════════════════
# 建议数据结构
# ═══════════════════════════════════════════════════════════════════

LEVELS = ("info", "suggest", "warn")


@dataclass
class Recommendation:
    id: str
    level: str                    # info | suggest | warn
    title: str
    detail: str
    source: str                   # rule:<id> | ai-narrative
    command: str = ""             # 建议执行链的命令草稿 (需人工/copilot 复核)
    pipeline: str = ""            # 关联管线 key (stage_cards)
    meta: Dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> Dict[str, Any]:
        return {"id": self.id, "level": self.level, "title": self.title,
                "detail": self.detail, "command": self.command,
                "pipeline": self.pipeline, "source": self.source, "meta": self.meta}


def _d(digest: Dict[str, Any], pipeline: str, key: str, default: Any = None) -> Any:
    """从 digest 里安全取数: digest["pipelines"][pipeline][key]。"""
    return digest.get("pipelines", {}).get(pipeline, {}).get(key, default)


def _stage(digest: Dict[str, Any], pipeline: str, io_key: str) -> Dict[str, Any]:
    return digest.get("pipelines", {}).get(pipeline, {}).get("stages", {}).get(io_key, {})


def _present(digest: Dict[str, Any], pipeline: str, io_key: str) -> bool:
    return bool(_stage(digest, pipeline, io_key).get("exists"))


# ═══════════════════════════════════════════════════════════════════
# 规则库 (每条规则一个函数, 返回 None 或 Recommendation)
# ═══════════════════════════════════════════════════════════════════

def _rule_start_clean(d: Dict[str, Any]) -> Optional[Recommendation]:
    clean = _present(d, "preprocessing", "d_clean")
    raw = _d(d, "public", "raw_reads", 0) or _d(d, "preprocessing", "clean_reads", 0)
    if not clean and (raw or _d(d, "public", "raw_sra", 0)):
        return Recommendation(
            id="start_preprocessing", level="suggest",
            title="Raw data present but preprocessing not run",
            detail=f"Raw reads/SRA detected (raw_reads={_d(d, 'public', 'raw_reads', 0)}, "
                   f"raw_sra={_d(d, 'public', 'raw_sra', 0)}) while 02_Preprocessing has "
                   "not been run. Run QC + host depletion before entering the discovery pipeline.",
            pipeline="preprocessing",
            command="python data_preprocessing_pipeline/data_preprocessing.py "
                    "--stage clean deplete --input_reads <RAW_FASTQ_DIR> "
                    "--output_dir <PROJECT_ROOT>/02_Preprocessing --host_db <HOST_DB>",
            source="rule:start_preprocessing")
    return None


def _rule_deplete_missing(d: Dict[str, Any]) -> Optional[Recommendation]:
    if _present(d, "preprocessing", "d_clean") and not _present(d, "preprocessing", "d_hostdep"):
        n_clean = _d(d, "preprocessing", "clean_reads", 0)
        return Recommendation(
            id="run_deplete", level="suggest",
            title="QC done but host depletion not run",
            detail=f"{n_clean} clean-read file(s) ready; the discovery pipeline's "
                   "--input_reads expects host-depleted reads — run deplete first.",
            pipeline="preprocessing",
            command="python data_preprocessing_pipeline/data_preprocessing.py "
                    "--stage deplete --output_dir <PROJECT_ROOT>/02_Preprocessing "
                    "--host_db <HOST_DB>",
            source="rule:run_deplete")
    return None


def _rule_depth_variance(d: Dict[str, Any]) -> Optional[Recommendation]:
    sizes = _d(d, "preprocessing", "clean_read_sizes") or []
    if len(sizes) < 4:
        return None
    mn, mx = min(sizes), max(sizes)
    if mn > 0 and mx / mn >= 10:
        return Recommendation(
            id="suggest_bbnorm", level="info",
            title="Large depth variance — consider normalization before co-assembly",
            detail=f"Clean-read file sizes span {mx / mn:.0f}x "
                   f"(min={mn / 1e6:.0f}MB, max={mx / 1e6:.0f}MB). Co-assembly would be "
                   "dominated by high-depth samples; if co-assembling, run --stage bbnorm "
                   "first. Note: bbnorm alters abundance profiles — keep using "
                   "non-normalized reads for the ④ quantification pipeline.",
            pipeline="preprocessing",
            command="python data_preprocessing_pipeline/data_preprocessing.py "
                    "--stage bbnorm --output_dir <PROJECT_ROOT>/02_Preprocessing",
            source="rule:suggest_bbnorm")
    return None


def _rule_single_sample(d: Dict[str, Any]) -> Optional[Recommendation]:
    n_asm = _d(d, "discovery", "assembly_samples", 0) or 0
    n_clean = _d(d, "preprocessing", "clean_reads", 0) or 0
    if n_asm == 0 and 0 < n_clean <= 2:
        return Recommendation(
            id="single_sample_coassembly", level="info",
            title="Very few samples — merge (co-assembly) stage can be skipped",
            detail=f"Only {n_clean} clean sample(s) (≤ 2); co-assembly (Flye merge) "
                   "adds little value. The orchestrator skips it automatically — "
                   "this is expected behavior.",
            pipeline="discovery", source="rule:single_sample_coassembly")
    return None


def _rule_discovery_chain(d: Dict[str, Any]) -> Optional[Recommendation]:
    stages = _stage(d, "discovery", "d_rescue")
    ident = _present(d, "discovery", "d_ident")
    if ident and not stages.get("exists"):
        missing = [k for k in ("d_filter", "d_cobra", "d_cluster", "d_taxonomy",
                               "d_checkv") if not _present(d, "discovery", k)]
        hint = ", ".join(missing) if missing else "intermediate stages"
        return Recommendation(
            id="continue_discovery", level="suggest",
            title="Identification done but downstream discovery stages not run",
            detail=f"Identification outputs exist while {hint} not run yet. The discovery "
                   "pipeline auto-fills missing stages by dependency — just resume "
                   "(checkpoint-safe).",
            pipeline="discovery",
            command="python virome_discovery_pipeline/virome_pipeline.py "
                    "--stage all --output_dir <PROJECT_ROOT>/03_Discovery "
                    "--input_reads <DEPLETED_READS_DIR>",
            source="rule:continue_discovery")
    return None


def _rule_analysis_report(d: Dict[str, Any]) -> Optional[Recommendation]:
    if _present(d, "analysis", "a_detect") and not _present(d, "analysis", "a_report"):
        return Recommendation(
            id="finish_analysis", level="suggest",
            title="④ Quantification done but report not generated",
            detail="Detection outputs exist while 09_Report is missing; the report stage "
                   "also produces the structured stats used for AI interpretation — "
                   "worth running.",
            pipeline="analysis",
            # ④ 手工校验 reads_dir/reference/ref_info/output_dir 四件套
            # (auto_known_virus.py), 只跑 report 也必须带上 —— 占位符缺值由 plan/agent 提示补
            command="python virome_analysis_pipeline/auto_known_virus.py "
                    "--stage report --reads_dir <DEPLETED_READS_DIR> "
                    "--reference <REF_FASTA> --ref_info <REF_INFO_TSV> "
                    "--output_dir <PROJECT_ROOT>/04_Analysis",
            source="rule:finish_analysis")
    return None


def _rule_eve_optional(d: Dict[str, Any]) -> Optional[Recommendation]:
    # ⑥ EVE 已归档 2026-10-11 → archive/20261011_endogenous_virus_pipeline
    # (方法不成熟): 不再向用户推荐已归档管线; 恢复模块时删除下面一行即可。
    return None
    genome = _present(d, "public", "h_genome")
    eve = _d(d, "eve", "present", False)
    if genome and not eve:
        return Recommendation(
            id="consider_eve", level="info",
            title="Host genome ready — optional EVE endogenous-virus screen",
            detail="EVE screening takes a host-genome FASTA and is independent of the "
                   "reads pipelines; run it when interested in endogenous viral "
                   "fossils / DNA-virus candidates.",
            pipeline="eve",
            command="python endogenous_virus_pipeline/eve_screen.py "
                    "-g <HOST_GENOME_FASTA> -n <SPECIES_NAME> -o <PROJECT_ROOT>/06_EVE -t 40",
            source="rule:consider_eve")
    return None


def _rule_submission_ready(d: Dict[str, Any]) -> Optional[Recommendation]:
    rescue = _present(d, "discovery", "d_rescue")
    if rescue and not _d(d, "submission", "present", False):
        return Recommendation(
            id="submission_ready", level="info",
            title="Rescue outputs ready — NCBI submission files can be generated",
            detail="The discovery pipeline has produced 10_Rescue; pipeline ⑦ can go "
                   "from rescue straight to .sqn in one step.",
            pipeline="submission",
            command="python virome_submission_pipeline/submission_pipeline.py discovery "
                    "--work-dir <PROJECT_ROOT>/03_Discovery/10_Rescue "
                    "--run-title <RUN_TITLE> --mode both",
            source="rule:submission_ready")
    return None


def _rule_nothing_present(d: Dict[str, Any]) -> Optional[Recommendation]:
    if not d.get("totals", {}).get("pipelines_present"):
        return Recommendation(
            id="empty_project", level="warn",
            title="Project directory empty or layout not recognized",
            detail="No pipeline outputs scanned. If the directory is truly empty, start "
                   "from ① public-data retrieval; for legacy projects, re-scan with "
                   "--io-layout legacy.",
            pipeline="", source="rule:empty_project")
    return None


RULES = [
    _rule_nothing_present,
    _rule_start_clean,
    _rule_deplete_missing,
    _rule_depth_variance,
    _rule_single_sample,
    _rule_discovery_chain,
    _rule_analysis_report,
    _rule_eve_optional,
    _rule_submission_ready,
]


# ═══════════════════════════════════════════════════════════════════
# 对外入口
# ═══════════════════════════════════════════════════════════════════

def advise(digest: Dict[str, Any]) -> List[Recommendation]:
    """纯规则顾问: 确定性、离线。返回按 warn > suggest > info 排序的建议。"""
    out: List[Recommendation] = []
    for rule in RULES:
        try:
            rec = rule(digest)
        except Exception:  # 单规则异常不影响其余
            continue
        if rec is not None:
            out.append(rec)
    rank = {"warn": 0, "suggest": 1, "info": 2}
    out.sort(key=lambda r: rank.get(r.level, 3))
    return out


def advise_with_ai(digest: Dict[str, Any], *,
                   state_dir: Optional[Path] = None,
                   config: Any = None,
                   client: Any = None,
                   compact_digest: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
    """规则骨架 + AI 叙事。返回 {"recommendations": [...], "narrative": {...}}。

    AI 只负责把规则建议串成执行路径叙事; 提示词明令不得新增建议、
    数字只准引用给定 JSON。离线时 narrative.status = "unavailable"。
    """
    recs = advise(digest)
    payload: Dict[str, Any] = {
        "recommendations": [r.to_dict() for r in recs],
        "narrative": {"status": "skipped", "reason": "AI not requested"},
    }
    if client is None:
        try:
            from mmpv_common.ai import ai_status
            ok, _ = ai_status(config)
        except Exception:
            ok = False
        if not ok:
            payload["narrative"] = {"status": "unavailable",
                                    "reason": "AI not configured; rule-based "
                                              "recommendations remain fully available"}
            return payload
    slim = compact_digest
    if slim is None:
        try:
            from mmpv_common.run_digest import compact_for_ai
            slim = compact_for_ai(digest)
        except Exception:
            slim = {"totals": digest.get("totals", {})}
    prompt_text, meta = render_prompt(
        "advisor_narrative",
        rules_json=json.dumps([r.to_dict() for r in recs], ensure_ascii=False)[:15000],
        digest_json=json.dumps(slim, ensure_ascii=False)[:15000])
    res = chat([{"role": "user", "content": prompt_text}], config=config, client=client)
    if res.status != "ok":
        payload["narrative"] = {"status": res.status, "reason": res.reason}
        return payload
    log_ai_event(state_dir, {"kind": "call", "op": "advisor_narrative",
                             "model": res.model, "prompt": meta.to_dict(),
                             "usage": res.usage,
                             "latency_s": round(res.latency_s, 2),
                             "n_rules": len(recs)})
    payload["narrative"] = {"status": "ok", "text": res.text.strip(),
                            "meta": {**res.to_dict(), "prompt": meta.to_dict()}}
    return payload
