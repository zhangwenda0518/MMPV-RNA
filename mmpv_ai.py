#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
mmpv_ai.py — MMPV AI 驱动层统一入口 (四层递进架构)
====================================================

  L1 状态   digest   项目 run digest 提取/落盘 (mmpv_common.run_digest)
  L1 知识   cards    阶段卡片查阅 (mmpv_common.stage_cards)
  L2 解读   ask      带 digest 上下文的一次性 AI 问答
  L3 决策   advise   规则顾问 (+ --ai 叙事增强)
  L3 决策   triage   候选三级标注 (只加 ai_* 列, 不动数据)
  L4 编排   chat     自然语言副驾驶 (NL → 白名单校验命令草稿, 不执行)
  L4 编排   validate 命令白名单校验 + dry-run 预览
  L4 编排   mcp      启动 MCP server (接入 ZCode/Claude/Cursor 等)

红线 (doc/AI_DRIVEN_DESIGN.md): AI 只解读/建议/标注/校验, 本入口
不执行任何管线命令; 执行权永远在人的终端。

示例:
  python mmpv_ai.py digest out/ --write
  python mmpv_ai.py advise out/ --ai
  python mmpv_ai.py triage out/03_Discovery/02_Identification/candidates.tsv \\
      --evidence genomad,virsorter2,viralverify --mode rules
  python mmpv_ai.py chat
  python mmpv_ai.py mcp
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from mmpv_common import run_digest, stage_cards                    # noqa: E402
from mmpv_common.ai import ai_status                               # noqa: E402
from mmpv_common.advisor import advise, advise_with_ai             # noqa: E402
from mmpv_common.copilot import plan_with_ai, suggest_intents, validate  # noqa: E402
from mmpv_common.triage import triage                              # noqa: E402


# ═══════════════════════════════════════════════════════════════════
# 子命令实现
# ═══════════════════════════════════════════════════════════════════

def _print_ai_banner() -> None:
    ok, why = ai_status()
    print(f"[AI] {'可用: ' + why if ok else '不可用 — ' + why}")


def cmd_digest(args) -> int:
    d = run_digest.build_digest(Path(args.project), layout=args.io_layout)
    if args.write:
        out = run_digest.write_digest(d, Path(args.project))
        print(f"[digest] 已写入 {out}", file=sys.stderr)
    payload = run_digest.compact_for_ai(d) if args.compact else d
    print(json.dumps(payload, ensure_ascii=False, indent=2))
    return 0


def cmd_cards(args) -> int:
    if args.pipeline:
        card = stage_cards.get_card(args.pipeline, args.stage)
        if card is None:
            print(f"未知 pipeline/stage: {args.pipeline}/{args.stage}\n"
                  f"可用管线: {', '.join(stage_cards.pipelines())}", file=sys.stderr)
            return 2
        print(json.dumps(card, ensure_ascii=False, indent=2))
    else:
        print(stage_cards.summary_text())
    return 0


def cmd_advise(args) -> int:
    d = run_digest.build_digest(Path(args.project), layout=args.io_layout)
    if args.ai:
        payload = advise_with_ai(d, state_dir=run_digest.state_dir(Path(args.project)))
        for r in payload["recommendations"]:
            print(f"[{r['level'].upper():7s}] {r['title']}\n    {r['detail']}")
            if r.get("command"):
                print(f"    建议: {r['command']}")
        n = payload["narrative"]
        print("\n── AI 总体叙事 " +
              (f"({n['status']})" if n["status"] != "ok" else "") + " ──")
        print(n.get("text") or n.get("reason", ""))
    else:
        for r in advise(d):
            print(f"[{r.level.upper():7s}] {r.title}")
            print(f"    {r.detail}")
            if r.command:
                print(f"    建议: {r.command}")
    return 0


def cmd_triage(args) -> int:
    if args.preset:
        from mmpv_common.triage import run_preset
        report = run_preset(args.preset, Path(args.project),
                            mode=args.mode,
                            out=Path(args.output) if args.output else None,
                            batch_size=args.batch,
                            state_dir=Path(args.state_dir) if args.state_dir else None)
        print(json.dumps(report, ensure_ascii=False, indent=2))
        return 0 if report.get("status") == "ok" else 2
    report = triage(Path(args.input),
                    Path(args.output) if args.output else None,
                    evidence_cols=args.evidence.split(",") if args.evidence else None,
                    id_col=args.id_col, mode=args.mode,
                    state_dir=Path(args.state_dir) if args.state_dir else None,
                    batch_size=args.batch)
    print(json.dumps(report, ensure_ascii=False, indent=2))
    return 0


def cmd_provenance(args) -> int:
    """AI 参与溯源一键导出 (确定性聚合, 零 AI 调用)。"""
    from mmpv_common.provenance_export import export_provenance
    res = export_provenance(Path(args.project), layout=args.io_layout)
    print(f"[provenance] JSON → {res['json']}")
    print(f"[provenance] MD   → {res['md']}")
    print(f"[provenance] AI 调用 {res['n_ai_calls']} 次 · agent {res['n_agent_runs']} 轮 · "
          f"triage {res['n_triage_files']} 份")
    return 0


def cmd_ask(args) -> int:
    context = ""
    if args.project:
        slim = run_digest.compact_for_ai(
            run_digest.build_digest(Path(args.project), layout=args.io_layout))
        context = "项目 run digest:\n" + json.dumps(slim, ensure_ascii=False)
    from mmpv_common.ai import ask as ai_ask
    res = ai_ask(args.question, context=context,
                 state_dir=run_digest.state_dir(Path(args.project)) if args.project else None)
    if res["status"] != "ok":
        print(f"[AI 不可用] {res.get('reason', '')}", file=sys.stderr)
        return 3
    print(res["answer"])
    return 0


def cmd_chat(args) -> int:
    """自然语言副驾驶 REPL: 意图 → 命令草稿 → 白名单校验 → dry-run 预览。"""
    _print_ai_banner()
    print("""MMPV 副驾驶 (不执行命令, 只产出校验过的命令草稿)
  试: "帮我对 out/ 项目跑发现管线全部阶段" / "质控一下 raw_fastqs"
  子命令: /digest <项目>  /advise <项目>  /quit""")
    while True:
        try:
            line = input("\nmmpv> ").strip()
        except (EOFError, KeyboardInterrupt):
            print()
            return 0
        if not line:
            continue
        if line in ("/quit", "/exit", "q"):
            return 0
        if line.startswith("/digest "):
            d = run_digest.build_digest(Path(line.split(None, 1)[1]))
            print(json.dumps(run_digest.compact_for_ai(d), ensure_ascii=False, indent=2)[:3000])
            continue
        if line.startswith("/advise "):
            d = run_digest.build_digest(Path(line.split(None, 1)[1]))
            for r in advise(d):
                print(f"[{r.level.upper():7s}] {r.title}\n    {r.detail}")
            continue
        # 意图规划: AI 优先, 规则兜底
        planned = False
        ok, _ = ai_status()
        if ok:
            res = plan_with_ai(line, state_dir=None)
            if res.get("status") == "ok":
                v = res["validation"]
                print(f"意图: {res.get('explain', '')}")
                if res.get("missing"):
                    print("缺参数: " + "; ".join(res["missing"]))
                print(v["preview"])
                for w in v["warnings"]:
                    print(f"  ⚠ {w}")
                if not v["ok"]:
                    print("  ✗ 校验未过: " + "; ".join(v["errors"]))
                planned = True
        if not planned:
            sugs = suggest_intents(line)
            if not sugs:
                print("没听懂。可以试试: /advise <项目>, 或描述管线名"
                      "(质控/发现/定量/进化/提交/下载; EVE/宿主表达已归档)")
                continue
            print("猜测你想做 (AI 未配置, 关键词匹配; 详细度有限):")
            for s in sugs:
                print(f"  ◆ {s['pipeline']} — {s['desc']}")
                print(f"    示例: {s['example']}")
            print("补全参数后可用 `validate` 子命令校验完整命令。")


def cmd_validate(args) -> int:
    res = validate(args.command)
    print(res.preview if res.ok else "校验未通过:")
    for e in res.errors:
        print(f"  ✗ {e}")
    for w in res.warnings:
        print(f"  ⚠ {w}")
    return 0 if res.ok else 1


def cmd_agent(args) -> int:
    """受控 agent 全链路编排: plan → 确认门 → execute → read_digest → iterate。"""
    from mmpv_common.agent import run_agent
    inputs = {}
    for pair in (args.set or []):
        if "=" not in pair:
            print(f"--set 格式应为 \"占位符=值\": {pair}", file=sys.stderr)
            return 2
        k, v = pair.split("=", 1)
        inputs[k.strip()] = v.strip()
    pipelines = args.only.split(",") if args.only else None
    report = run_agent(
        Path(args.project), inputs=inputs, pipelines=pipelines,
        max_steps=args.max_steps, dry_run=args.dry_run, yes=args.yes,
        confirm_fn=(lambda p: input(p).strip().lower() in ("y", "yes"))
        if not args.dry_run else None,
    )
    print(f"[agent] 项目: {report.project_dir}")
    for s in report.steps:
        mark = {"ran": "✓", "planned": "◦", "blocked": "⊘", "skipped": "–",
                "unconfirmed": "?", "failed": "✗"}.get(s.status, "·")
        print(f"  {mark} [{s.status:12s}] {s.title}")
        print(f"      {s.command}")
        if s.detail:
            print(f"      ↳ {s.detail}")
    print(f"[agent] {report.summary_line()}")
    if report.journal_path:
        print(f"[agent] 留痕: {report.journal_path}")
    # 退出码契约 (与 agent.py 一致): 1=failed 停链; 2=有 blocked 步骤 (前置缺失, 可继续)
    if report.stop_reason == "failed":
        return 1
    if any(s.status == "blocked" for s in report.steps):
        return 2
    return 0


def cmd_report(args) -> int:
    """项目级主报告: 结果 × 过程 × AI 参与 (90_Handoff/MMPV_master_report.html)。"""
    from mmpv_common.process_report import generate_master_report
    res = generate_master_report(Path(args.project), layout=args.io_layout,
                                 out_path=Path(args.out) if args.out else None,
                                 want_ai=not args.no_ai)
    print(f"[report] 主报告 → {res['path']}")
    print(f"[report] AI 过程叙事: {res['ai_narrative']}"
          + (" (未配置 key, 过程面板/时间线/建议不受影响)" if res["ai_narrative"] != "ok" else ""))
    print(f"[report] 规则建议 {res['n_recommendations']} 条已收录")
    return 0


def cmd_enrich(args) -> int:
    """事后补 AI: 对已生成的报告 HTML 补插 AI 解读 + 过程面板 (幂等, 无需重跑管线)。"""
    from mmpv_common.process_report import enrich_project_reports
    res = enrich_project_reports(Path(args.project), layout=args.io_layout,
                                 want_ai=not args.no_ai)
    for r in res["reports"]:
        print(f"  {r['report']}: AI={r['ai']} · 面板={r['panel']}")
    m = res["master"]
    print(f"[enrich] 主报告 → {m['path']} (AI 叙事: {m['ai_narrative']})")
    if not res["reports"] and m.get("ai_narrative") not in ("ok",):
        print("[enrich] 未发现可补的报告 (先跑管线报告), 或 AI 未配置 (--no-ai 可只补面板)")
    return 0


def cmd_mcp(args) -> int:
    from mmpv_common import mcp_server
    return mcp_server.serve()


# ═══════════════════════════════════════════════════════════════════

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="mmpv_ai",
                                description="MMPV AI 驱动层统一入口 (四层架构)")
    sub = p.add_subparsers(dest="cmd", required=True)

    sp = sub.add_parser("digest", help="项目 run digest 提取 (L1)")
    sp.add_argument("project")
    sp.add_argument("--io-layout", choices=["standard", "legacy"], default=None)
    sp.add_argument("--write", action="store_true", help="落盘 90_Handoff/state/")
    sp.add_argument("--compact", action="store_true", help="AI 瘦身版")
    sp.set_defaults(fn=cmd_digest)

    sp = sub.add_parser("cards", help="阶段卡片查阅 (L1)")
    sp.add_argument("pipeline", nargs="?", default=None)
    sp.add_argument("stage", nargs="?", default=None)
    sp.set_defaults(fn=cmd_cards)

    sp = sub.add_parser("advise", help="下一步建议 (L3)")
    sp.add_argument("project")
    sp.add_argument("--ai", action="store_true", help="附加 AI 总体叙事 (需配置)")
    sp.add_argument("--io-layout", choices=["standard", "legacy"], default=None)
    sp.set_defaults(fn=cmd_advise)

    sp = sub.add_parser("triage", help="候选三级标注 (L3, 只加 ai_* 列; 或 --preset 实战预设)")
    sp.add_argument("input", nargs="?", default=None, help="输入 TSV (或改用 --preset)")
    sp.add_argument("-o", "--output", default=None)
    sp.add_argument("--evidence", default=None, help="证据列, 逗号分隔")
    sp.add_argument("--id-col", default=None)
    sp.add_argument("--mode", choices=["rules", "ai"], default="rules")
    sp.add_argument("--state-dir", default=None)
    sp.add_argument("--batch", type=int, default=20)
    sp.add_argument("--preset", choices=["ident", "eve_distinguish"], default=None,
                    help="实战预设: 定位/构建 MMPV 真实证据表")
    sp.add_argument("--project", default=None, help="预设用的项目根目录")
    sp.set_defaults(fn=cmd_triage)

    sp = sub.add_parser("provenance", help="AI 参与溯源一键导出 (JSON+MD, 方法章一键引用)")
    sp.add_argument("project")
    sp.add_argument("--io-layout", choices=["standard", "legacy"], default=None)
    sp.set_defaults(fn=cmd_provenance)

    sp = sub.add_parser("ask", help="一次性 AI 问答 (L2, 可带项目上下文)")
    sp.add_argument("question")
    sp.add_argument("--project", default=None)
    sp.add_argument("--io-layout", choices=["standard", "legacy"], default=None)
    sp.set_defaults(fn=cmd_ask)

    sp = sub.add_parser("chat", help="自然语言副驾驶 REPL (L4)")
    sp.set_defaults(fn=cmd_chat)

    sp = sub.add_parser("validate", help="命令白名单校验 + dry-run (L4)")
    sp.add_argument("command")
    sp.set_defaults(fn=cmd_validate)

    sp = sub.add_parser("agent", help="受控 agent 全链路编排 (L4): plan→确认→execute→digest→iterate")
    sp.add_argument("project")
    sp.add_argument("--dry-run", action="store_true", help="只产出计划, 不执行")
    sp.add_argument("--yes", action="store_true",
                    help="非危险步骤免确认 (危险步骤 --force 类仍须交互确认, 红线)")
    sp.add_argument("--max-steps", type=int, default=8, help="单次预算上限 (默认 8)")
    sp.add_argument("--only", default=None, help="限定管线, 逗号分隔 (如 preprocessing,discovery)")
    sp.add_argument("--set", action="append", default=[],
                    help="占位符赋值, 可多次: --set \"RAW_FASTQ_DIR=D:/raw\"")
    sp.set_defaults(fn=cmd_agent)

    sp = sub.add_parser("report", help="项目级主报告: 结果×过程×AI 参与 (L2+过程层)")
    sp.add_argument("project")
    sp.add_argument("--io-layout", choices=["standard", "legacy"], default=None)
    sp.add_argument("--out", default=None, help="输出 HTML (默认 90_Handoff/MMPV_master_report.html)")
    sp.add_argument("--no-ai", action="store_true", help="跳过 AI 过程叙事 (面板/时间线/建议仍生成)")
    sp.set_defaults(fn=cmd_report)

    sp = sub.add_parser("enrich", help="事后补 AI: 已生成报告补插 AI 解读+过程面板 (幂等)")
    sp.add_argument("project")
    sp.add_argument("--io-layout", choices=["standard", "legacy"], default=None)
    sp.add_argument("--no-ai", action="store_true", help="只补过程面板, 不调 AI")
    sp.set_defaults(fn=cmd_enrich)

    sp = sub.add_parser("mcp", help="启动 MCP server (L4)")
    sp.set_defaults(fn=cmd_mcp)
    return p


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    return args.fn(args)


if __name__ == "__main__":
    sys.exit(main())
