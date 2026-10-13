#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
triage.py — MMPV 候选三级 triage 写入器 (四层递进架构 · L3 决策层)
===================================================================

把 eve_distinguish/s4_filter.py 的 verdict→action (REMOVE/MOVE_EVE/
KEEP_virus/REVIEW) 模式推广为通用标注器: 对任意证据 TSV 的每一行给出

    KEEP    证据一致支持保留
    REVIEW  证据冲突 / 信息不足 / AI 输出不可解析 (保守回退)

【宪法】沿袭 gsa_sra.info.py 的"只回填空白"纪律并进一步收紧:
  * 只「追加」ai_action / ai_reason / ai_source 三列;
  * 原有列逐字节原样回写 (tests/test_triage.py 有不变量断言);
  * 本模块永远不执行删除/移动 —— REMOVE 类动作必须由人或在
    eve_distinguish 的规则链里做, AI 只标 REVIEW。
  * 列名冲突时拒绝写入 (输出到新文件), 绝不静默覆盖。

用法:
  from mmpv_common.triage import triage
  report = triage("candidates.tsv", out_path="annotated.tsv",
                  evidence_cols=["genomad","virsorter2","viralverify"],
                  id_col="contig_id", mode="rules")        # 离线
  report = triage(..., mode="ai", state_dir=...)           # AI 批注 (需配置)

CLI: python -m mmpv_common.triage <in.tsv> [-o out.tsv] [--mode rules|ai]
     [--evidence colA,colB] [--id-col contig_id] [--state-dir DIR]
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional

try:
    from mmpv_common.ai import AIConfig, ai_status, chat, log_ai_event, parse_json, render_prompt
except ImportError:
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
    from mmpv_common.ai import AIConfig, ai_status, chat, log_ai_event, parse_json, render_prompt  # type: ignore

TRIAGE_COLS = ("ai_action", "ai_reason", "ai_source")
VALID_ACTIONS = ("KEEP", "REVIEW")
DEFAULT_BATCH = 20


# ═══════════════════════════════════════════════════════════════════
# 规则模式 (确定性)
# ═══════════════════════════════════════════════════════════════════

def _rule_annotate(row: Dict[str, str], evidence_cols: List[str]) -> Dict[str, str]:
    vals = [str(row.get(c, "")).strip().lower() for c in evidence_cols]
    vals = [v for v in vals if v not in ("", "na", "n/a", "none", "nan", "-")]
    if not vals:
        return {"ai_action": "REVIEW", "ai_reason": "evidence columns all empty — cannot assess",
                "ai_source": "rule:unanimous"}
    unanimous = all(v == vals[0] for v in vals)
    if unanimous:
        if vals[0] in ("0", "1"):
            return {"ai_action": "KEEP",
                    "ai_reason": f"{vals.count('1')}/{len(vals)} evidence sources unanimously support",
                    "ai_source": "rule:unanimous"}
        return {"ai_action": "KEEP",
                "ai_reason": f"evidence unanimous: {vals[0]}", "ai_source": "rule:unanimous"}
    if all(v in ("0", "1") for v in vals):
        return {"ai_action": "REVIEW",
                "ai_reason": f"evidence split: {vals.count('1')}/{len(vals)} support",
                "ai_source": "rule:conflict"}
    return {"ai_action": "REVIEW",
            "ai_reason": f"evidence conflict: {'/'.join(vals)}",
            "ai_source": "rule:conflict"}


# ═══════════════════════════════════════════════════════════════════
# 实战预设 (P2): 对 MMPV 真实产物表定位/构建证据, 不臆造列名
# ═══════════════════════════════════════════════════════════════════

# ② 鉴定阶段产出 *.result.id 的 9 个工具 (report_pipeline.py _collect_identification 同源;
#    第 10 工具 Viroid BLASTN 走 blastn 管线, 不产 result.id)
IDENT_TOOLS = ['genomad', 'blast', 'metabuli', 'virsorter2', 'viralverify',
               'virhunter', 'virbot', 'viralm', 'rdrpcatch']

PRESETS = {
    "ident": "③ 鉴定证据矩阵: 现场合并 02_Identification 各工具 *.result.id → contig×工具 1/0 矩阵",
    # ⑥ 已归档 2026-10-11 → archive/20261011_endogenous_virus_pipeline; 预设保留,
    # 仅供消费历史项目的 dna_vs_eve_filter.tsv (新项目不再产出该文件)
    "eve_distinguish": "⑥ eve_distinguish【已归档】: 定位历史 dna_vs_eve_filter.tsv, 证据列自动侦测",
}


def locate_ident_dir(project: Path, layout: Optional[str] = None) -> Optional[Path]:
    """定位鉴定目录 (standard <proj>/03_Discovery/02_Identification; legacy 平铺)。"""
    from mmpv_common import io_layout
    from mmpv_common.run_digest import PIPELINE_ROOTS
    dirs = io_layout.layout_dirs(io_layout.resolve_layout(layout))
    cands = [Path(project) / PIPELINE_ROOTS["discovery"] / dirs["d_ident"],
             Path(project) / dirs["d_ident"]]
    for c in cands:
        if c.is_dir():
            return c
    return None


def build_ident_evidence(ident_dir: Path, out_path: Optional[Path] = None,
                         tools: Optional[List[str]] = None) -> Dict[str, Any]:
    """从各工具 *.result.id + all.candidate.fasta 现场构建 contig×工具 1/0 证据矩阵。

    宪法: 值只有 1 (该工具判病毒) / 0 (未判); "全票"分母 = 实际产出了
    result.id 的工具 (tools_present), 不用 9 工具全集 —— 缺库未跑的工具
    不应稀释一致性; 零工具命中的 contig 不进矩阵 (对 triage 无意义)。
    """
    tools = tools or IDENT_TOOLS
    ident_dir = Path(ident_dir)
    sample_dirs = [d for d in sorted(ident_dir.iterdir()) if d.is_dir()]
    # 第一遍: 哪些工具真实在跑 (有 result.id 文件)
    tools_present = [t for t in tools
                     if any((d / f"{d.name}_virus.{t}.result.id").is_file()
                            for d in sample_dirs)]
    if not tools_present:
        tools_present = list(tools)   # 极端: 只有 fasta 无工具结果 → 全 0 列兜底
    hits: Dict[str, Dict[str, int]] = {}
    sample_of: Dict[str, str] = {}
    for sample_dir in sample_dirs:
        for tool in tools_present:
            idf = sample_dir / f"{sample_dir.name}_virus.{tool}.result.id"
            if not idf.is_file():
                continue
            try:
                for line in idf.read_text(encoding="utf-8", errors="replace").splitlines():
                    cid = line.strip().lstrip(">").split()[0] if line.strip() else ""
                    if cid:
                        hits.setdefault(cid, {})[tool] = 1
                        sample_of.setdefault(cid, sample_dir.name)
            except OSError:
                continue
    if out_path is None:
        out_path = ident_dir / "ident_evidence.tsv"
    out_path = Path(out_path)
    with open(out_path, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(["contig_id", "sample"] + tools_present)
        for cid in sorted(hits):
            row = hits[cid]
            w.writerow([cid, sample_of.get(cid, "")] +
                       [row.get(t, 0) for t in tools_present])
    n_unanimous = sum(1 for h in hits.values() if len(h) == len(tools_present))
    return {"status": "ok", "evidence": str(out_path),
            "n_samples": len(sample_dirs), "n_contigs": len(hits),
            "tools_present": tools_present,
            "n_unanimous": n_unanimous,
            "n_partial": len(hits) - n_unanimous}


def _locate_eve_distinguish(project: Path) -> Optional[Path]:
    """定位 dna_vs_eve_filter.tsv (06_EVE 优先, 兼 legacy/其他深度)。"""
    proj = Path(project)
    for base in (proj / "06_EVE", proj / "03_Discovery", proj):
        if not base.is_dir():
            continue
        try:
            hits = sorted(base.rglob("dna_vs_eve_filter.tsv"))
        except OSError:
            hits = []
        if hits:
            return hits[0]
    return None


def run_preset(preset: str, project: Path, *, mode: str = "rules",
               out: Optional[Path] = None, batch_size: int = DEFAULT_BATCH,
               state_dir: Optional[Path] = None,
               config: Any = None, client: Any = None) -> Dict[str, Any]:
    """实战预设入口: 定位/构建证据表 → 走同一套 triage() (宪法不变)。

    返回 dict (含 triage 报告或 error 说明); 找不到目标文件时 status=error,
    不抛异常 —— CLI/agent 可优雅提示。
    """
    project = Path(project)
    if preset == "ident":
        ident = locate_ident_dir(project)
        if ident is None:
            return {"status": "error",
                    "reason": "未找到鉴定目录 (03_Discovery/02_Identification); 先跑 discovery identification"}
        built = build_ident_evidence(ident)
        ev = Path(built["evidence"])
        report = triage(ev, out, evidence_cols=built["tools_present"],
                        id_col="contig_id", mode=mode, state_dir=state_dir,
                        config=config, client=client, batch_size=batch_size)
        report["preset"] = "ident"
        report["evidence_built"] = built
        return report
    if preset == "eve_distinguish":
        tsv = _locate_eve_distinguish(project)
        if tsv is None:
            return {"status": "error",
                    "reason": "未找到 dna_vs_eve_filter.tsv; ⑥ 已归档 (2026-10-11 → "
                              "archive/20261011_endogenous_virus_pipeline), 历史项目请确认 "
                              "06_EVE 输出存在, 新项目无法再产出该文件"}
        report = triage(tsv, out, mode=mode, state_dir=state_dir,
                        config=config, client=client, batch_size=batch_size)
        report["preset"] = "eve_distinguish"
        report["source"] = str(tsv)
        return report
    return {"status": "error", "reason": f"未知预设: {preset} (可用: {', '.join(PRESETS)})"}


# ═══════════════════════════════════════════════════════════════════
# AI 模式 (批处理 + schema 强制 + 保守回退)
# ═══════════════════════════════════════════════════════════════════

def _ai_batch(rows: List[Dict[str, str]], id_col: str, evidence_cols: List[str],
              config: Any, client: Any, model_used: List[str]) -> Dict[str, Dict[str, str]]:
    """一批行 → {行id: 标注}。schema 不合的行回退 REVIEW。"""
    compact = []
    for r in rows:
        compact.append({"id": r.get(id_col, ""),
                        **{c: r.get(c, "") for c in evidence_cols}})
    prompt_text, meta = render_prompt("row_triage",
                                      rows_json=json.dumps(compact, ensure_ascii=False))
    res = chat([{"role": "user", "content": prompt_text}],
               config=config, client=client, json_mode=True, temperature=0.0)
    fallback = lambda why: {"ai_action": "REVIEW", "ai_reason": why,
                            "ai_source": "ai-fallback"}  # noqa: E731
    if res.status != "ok":
        return {r.get(id_col, ""): fallback(f"AI unavailable: {res.reason[:60]}") for r in rows}
    model_used.append(res.model)
    try:
        parsed = parse_json(res.text)
    except ValueError as e:
        return {r.get(id_col, ""): fallback(str(e)[:60]) for r in rows}
    if not isinstance(parsed, list):
        parsed = parsed.get("rows", []) if isinstance(parsed, dict) else []
    out: Dict[str, Dict[str, str]] = {}
    seen = set()
    for item in parsed:
        if not isinstance(item, dict):
            continue
        rid = str(item.get("id", ""))
        action = str(item.get("action", "")).upper()
        reason = str(item.get("reason", ""))[:120]
        if rid and action in VALID_ACTIONS and rid not in seen:
            out[rid] = {"ai_action": action, "ai_reason": reason or "AI gave no reason",
                        "ai_source": f"model:{res.model}"}
            seen.add(rid)
    for r in rows:  # 漏标/错schema → 保守 REVIEW
        rid = r.get(id_col, "")
        if rid not in out:
            out[rid] = fallback("row not covered by AI output; conservative REVIEW fallback")
    return out


# ═══════════════════════════════════════════════════════════════════
# 主入口
# ═══════════════════════════════════════════════════════════════════

def triage(input_tsv: Path, out_path: Optional[Path] = None, *,
           evidence_cols: Optional[List[str]] = None,
           id_col: Optional[str] = None,
           mode: str = "rules",
           state_dir: Optional[Path] = None,
           config: Any = None,
           client: Any = None,
           batch_size: int = DEFAULT_BATCH) -> Dict[str, Any]:
    """对 TSV 逐行标注, 输出新文件。返回报告 dict (counts / output / status)。

    不变量 (tests/test_triage.py 断言): 输出文件的原有列与值同输入逐字节
    一致, 仅在尾部追加 TRIAGE_COLS 三列。
    """
    input_tsv = Path(input_tsv)
    if mode not in ("rules", "ai"):
        raise ValueError(f"mode 必须是 rules|ai, 得到 {mode!r}")
    if out_path is None:
        out_path = input_tsv.with_name(input_tsv.stem + ".triage.tsv")
    out_path = Path(out_path)
    if out_path.resolve() == input_tsv.resolve():
        out_path = input_tsv.with_name(input_tsv.stem + ".triage.tsv")

    with open(input_tsv, "r", encoding="utf-8", errors="replace", newline="") as f:
        reader = csv.reader(f, delimiter="\t")
        try:
            header = next(reader)
        except StopIteration:
            raise ValueError(f"空 TSV: {input_tsv}")
        rows = [dict(zip(header, r + [""] * (len(header) - len(r)))) for r in reader]

    # 宪法检查: 三列已在表中 → 拒绝, 要求换输出路径
    clash = [c for c in TRIAGE_COLS if c in header]
    if clash:
        raise ValueError(f"输入已含标注列 {clash}, 拒绝覆盖; 请用 -o 指定新输出文件")

    # 列名冲突防御: 数据里恰有 ai_action 之外的匹配已由 clash 覆盖。
    if id_col is None:
        for cand in ("contig_id", "id", "name", "seq_id", "record_id"):
            if cand in header:
                id_col = cand
                break
        else:
            id_col = header[0]
    if evidence_cols is None:
        evidence_cols = [c for c in header if c not in (id_col,) and not c.startswith("ai_")][:6]

    counts = {"KEEP": 0, "REVIEW": 0}
    models: List[str] = []
    annotated: List[Dict[str, str]] = []
    if mode == "rules":
        for row in rows:
            mark = _rule_annotate(row, evidence_cols)
            counts[mark["ai_action"]] += 1
            annotated.append({**row, **mark})
    else:
        ok, why = ai_status(config)
        if not ok and client is None:
            raise SystemExit(f"AI 模式不可用: {why} (可先 --mode rules 离线跑)")
        for i in range(0, len(rows), batch_size):
            batch = rows[i:i + batch_size]
            marks = _ai_batch(batch, id_col, evidence_cols, config, client, models)
            for row in batch:
                mark = marks.get(row.get(id_col, ""),
                                 {"ai_action": "REVIEW", "ai_reason": "unknown row",
                                  "ai_source": "ai-fallback"})
                counts[mark["ai_action"]] += 1
                annotated.append({**row, **mark})
        if state_dir is not None and models:
            log_ai_event(state_dir, {
                "kind": "decision", "op": "triage",
                "model": models[0], "mode": "ai",
                "input": str(input_tsv), "n_rows": len(rows),
                "counts": counts,
                "prompt": "row_triage@%s" % _prompt_version(),
                "evidence_cols": evidence_cols})

    with open(out_path, "w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f, delimiter="\t", lineterminator="\n")
        writer.writerow(header + list(TRIAGE_COLS))
        for row in annotated:
            writer.writerow([row.get(c, "") for c in header] +
                            [row.get(c, "") for c in TRIAGE_COLS])

    return {"status": "ok", "mode": mode, "input": str(input_tsv),
            "output": str(out_path), "n_rows": len(rows),
            "id_col": id_col, "evidence_cols": evidence_cols,
            "counts": counts, "models": sorted(set(models))}


def _prompt_version() -> str:
    from mmpv_common.ai import PROMPTS
    return sorted(PROMPTS.get("row_triage", {}))[-1] if PROMPTS.get("row_triage") else "?"


# ═══════════════════════════════════════════════════════════════════

def main(argv: Optional[List[str]] = None) -> int:
    ap = argparse.ArgumentParser(description="MMPV 通用 triage 标注器 (只加 ai_* 列)")
    ap.add_argument("input", nargs="?", default=None, help="输入 TSV (或改用 --preset)")
    ap.add_argument("-o", "--output", default=None, help="输出 TSV (默认 <stem>.triage.tsv)")
    ap.add_argument("--evidence", default=None,
                    help="证据列 (逗号分隔; 默认自动取前 6 个非 id 列)")
    ap.add_argument("--id-col", default=None, help="行 id 列 (默认自动识别)")
    ap.add_argument("--mode", choices=["rules", "ai"], default="rules")
    ap.add_argument("--state-dir", default=None, help="AI 溯源日志目录 (90_Handoff/state)")
    ap.add_argument("--batch", type=int, default=DEFAULT_BATCH)
    ap.add_argument("--preset", choices=sorted(PRESETS), default=None,
                    help="实战预设: 定位/构建 MMPV 真实证据表 (ident / eve_distinguish)")
    ap.add_argument("--project", default=None, help="预设用的项目根目录")
    args = ap.parse_args(argv)

    if args.preset:
        if not args.project:
            ap.error("--preset 需要 --project <项目根>")
        report = run_preset(args.preset, Path(args.project),
                            mode=args.mode,
                            out=Path(args.output) if args.output else None,
                            batch_size=args.batch,
                            state_dir=Path(args.state_dir) if args.state_dir else None)
        print(json.dumps(report, ensure_ascii=False, indent=2))
        return 0 if report.get("status") == "ok" else 2
    if not args.input:
        ap.error("需要输入 TSV, 或 --preset <name> --project <项目根>")

    evidence = args.evidence.split(",") if args.evidence else None
    report = triage(Path(args.input), Path(args.output) if args.output else None,
                    evidence_cols=evidence, id_col=args.id_col, mode=args.mode,
                    state_dir=Path(args.state_dir) if args.state_dir else None,
                    batch_size=args.batch)
    print(json.dumps(report, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
