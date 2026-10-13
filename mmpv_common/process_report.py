#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
process_report.py — 过程报告层: 让最终报告承载整个过程 (不只读数据)
====================================================================

回应的需求: 报告此前只解读"本管线的结果数字"; 本模块把**过程维度**注入报告:

  * digest 时间线   项目怎么一步步长出来的 (digest_history.jsonl 快照序列);
  * AI 参与清单     哪些环节调过 AI、做了什么、用了哪个模型/提示词版本、
                    花了多少 token (ai_log.jsonl 聚合);
  * agent 留痕      受控编排跑过几轮、每轮结局 (90_Handoff/state/agent_runs/);
  * 治理声明        AI 只建议/解读/标注, 科学判定均为确定性算法 (审稿友好)。

三件交付:
  1. build_process_panel_html()   确定性 HTML 面板 —— **不需要 AI、不需要网络**,
                                  任何报告嵌入后立刻能展示全过程骨架;
  2. append_process_panel_to_html()  写盘点: 在已生成 HTML 的 </body> 前嵌入面板
                                  (幂等: 已嵌过不再重复; 无 90_Handoff 则原样不动);
  3. generate_master_report()     项目级主报告 (90_Handoff/MMPV_master_report.html):
                                  管线全景 + 过程时间线 + AI 参与治理 + 当前建议
                                  + (可选) AI 过程叙事 (process_narrative 提示词)。

纯标准库 + mmpv_common; HTML 渲染走 mmpv_common/report_template.py (Jinja2 模板,
templates/process_panel.html.j2 / master_report.html.j2)。AI 叙事离线时明确标注
"未配置", 主报告照常生成。
"""

from __future__ import annotations

import html as _html
import json
import time
from pathlib import Path
from typing import Any, Dict, List, Optional

try:
    from mmpv_common import run_digest, stage_cards
    from mmpv_common.advisor import advise
    from mmpv_common.ai import (_ai_cache_file, chat, log_ai_event,
                                render_prompt)
    from mmpv_common.report_style import PALETTE
    from mmpv_common.report_template import render
except ImportError:
    import sys
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
    from mmpv_common import run_digest, stage_cards  # type: ignore
    from mmpv_common.advisor import advise  # type: ignore
    from mmpv_common.ai import _ai_cache_file, chat, log_ai_event, render_prompt  # type: ignore
    from mmpv_common.report_style import PALETTE  # type: ignore
    from mmpv_common.report_template import render  # type: ignore

PANEL_ID = "mmpv-process-panel"
MASTER_REPORT_NAME = "MMPV_master_report.html"

# 面板/主报告的管线展示顺序 (主链 ①→⑧; host_expression 为可选补充殿后)
# 【⑥ EVE、⑧ 宿主表达均已归档 2026-10-11 → archive/20261011_{endogenous_virus,host_expression}_pipeline;
#   入口保留, 历史项目有对应输出根时仍能挂面板, 新项目自然空缺】
_DISPLAY_ORDER = ["public", "preprocessing", "discovery", "analysis",
                  "phylo", "eve", "submission", "host_expression"]

_LEVEL_COLOR = {"warn": PALETTE["danger"], "suggest": PALETTE["warning"],
                "info": PALETTE["info"]}


def _esc(v: Any) -> str:
    return _html.escape(str(v) if v is not None else "")


# ═══════════════════════════════════════════════════════════════════
# 过程记录聚合
# ═══════════════════════════════════════════════════════════════════

def _summarize_ai_log(state: Path) -> Dict[str, Any]:
    """ai_log.jsonl → 参与摘要 (调用数/操作分布/模型/提示词版本/token/按月用量)。"""
    out: Dict[str, Any] = {"n_calls": 0, "by_op": {}, "by_op_tokens": {}, "by_month": {},
                           "models": [], "prompt_versions": [], "total_tokens": 0,
                           "n_cached": 0, "first_ts": None, "last_ts": None,
                           "agent_runs_logged": 0}
    p = state / "ai_log.jsonl"
    if not p.is_file():
        return out
    ops: Dict[str, int] = {}
    op_tokens: Dict[str, int] = {}
    months: Dict[str, Dict[str, int]] = {}
    models: List[str] = []
    versions: List[str] = []
    for line in p.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line:
            continue
        try:
            rec = json.loads(line)
        except json.JSONDecodeError:
            continue
        if rec.get("kind") == "agent":
            out["agent_runs_logged"] += 1
            continue
        out["n_calls"] += 1
        op = rec.get("op", "?")
        ops[op] = ops.get(op, 0) + 1
        usage = rec.get("usage") or {}
        tok = usage.get("total_tokens")
        n_tok = int(tok) if isinstance(tok, (int, float)) else 0
        if n_tok:
            op_tokens[op] = op_tokens.get(op, 0) + n_tok
            out["total_tokens"] += n_tok
        ts = rec.get("ts", "")
        ym = ts[:7] if len(ts) >= 7 else "?"
        slot = months.setdefault(ym, {"calls": 0, "tokens": 0})
        slot["calls"] += 1
        slot["tokens"] += n_tok
        if rec.get("model") and rec["model"] not in models:
            models.append(rec["model"])
        pm = rec.get("prompt") or {}
        pv = f"{pm.get('name', '?')}@{pm.get('version', '?')}"
        if pm and pv not in versions:
            versions.append(pv)
        if rec.get("cached"):
            out["n_cached"] += 1
        if ts:
            if out["first_ts"] is None:
                out["first_ts"] = ts
            out["last_ts"] = ts
    out["by_op"] = dict(sorted(ops.items(), key=lambda kv: -kv[1]))
    out["by_op_tokens"] = op_tokens
    out["by_month"] = dict(sorted(months.items()))
    out["models"] = models
    out["prompt_versions"] = versions
    return out


def _summarize_agent_runs(state: Path) -> List[Dict[str, Any]]:
    """agent_runs/*.json → [{finished_at, summary, stop_reason}] (最多 10 条)。"""
    d = state / "agent_runs"
    if not d.is_dir():
        return []
    runs = []
    for p in sorted(d.glob("*.json")):
        try:
            rec = json.loads(p.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            continue
        runs.append({"finished_at": rec.get("finished_at", p.stem),
                     "summary": rec.get("summary", ""),
                     "stop_reason": rec.get("stop_reason", "")})
    return runs[-10:]


def collect_triage(project: Path, layout: Optional[str] = None) -> List[Dict[str, Any]]:
    """扫描项目各管线根下的 *.triage.tsv → KEEP/REVIEW 计数 (面板/主报告数据源)。"""
    from mmpv_common import io_layout
    from mmpv_common.run_digest import PIPELINE_ROOTS
    proj = Path(project)
    lay = io_layout.resolve_layout(layout)
    out: List[Dict[str, Any]] = []
    seen = set()
    for pk, root_name in PIPELINE_ROOTS.items():
        base = proj / root_name if lay == "standard" else proj
        d = base / "."
        if not d.is_dir():
            continue
        try:
            files = sorted(base.rglob("*.triage.tsv"))
        except OSError:
            continue
        for f in files:
            key = str(f.resolve())
            if key in seen:
                continue
            seen.add(key)
            keep = review = 0
            rows = 0
            try:
                with open(f, "r", encoding="utf-8", errors="replace") as fh:
                    header = fh.readline().rstrip("\n").split("\t")
                    if "ai_action" in header:
                        idx = header.index("ai_action")
                        for line in fh:
                            parts = line.rstrip("\n").split("\t")
                            if len(parts) > idx:
                                rows += 1
                                if parts[idx] == "KEEP":
                                    keep += 1
                                elif parts[idx] == "REVIEW":
                                    review += 1
            except OSError:
                continue
            try:
                rel = f.resolve().relative_to(proj.resolve()).as_posix()
            except ValueError:
                rel = str(f)
            out.append({"file": rel, "pipeline": pk, "rows": rows,
                        "keep": keep, "review": review})
    return out


def collect_pipeline_reports(project: Path, layout: Optional[str] = None) -> List[Dict[str, str]]:
    """已知落点的管线 HTML 报告清单 (主报告导航/面板互链用)。

    只登记确定路径 (各报告 main 的默认输出), 不做模糊搜索 —— 链接指向的
    文件存在与否由 exists 字段如实反映, 绝不给不存在的报告造链接。
    """
    from mmpv_common import io_layout
    from mmpv_common.run_digest import PIPELINE_ROOTS
    proj = Path(project)
    lay = io_layout.resolve_layout(layout)
    dirs = io_layout.layout_dirs(lay)
    std = lay == "standard"

    def _root(pk: str) -> Path:
        return proj / PIPELINE_ROOTS[pk] if std else proj

    def _first_html(d: Path) -> Optional[Path]:
        try:
            hits = sorted(d.glob("*.html"))
        except OSError:
            return None
        return hits[0] if hits else None

    def _glob_first(base: Path, pattern: str) -> Optional[Path]:
        """受 pattern 约束的首个匹配 (目录不存在 → None; 不做无约束模糊搜索)。"""
        try:
            hits = sorted(base.glob(pattern)) if base.is_dir() else []
        except OSError:
            return None
        return hits[0] if hits else None

    def _rglob_first(base: Path, name: str) -> Optional[Path]:
        """固定文件名在管线根下的首个匹配 (⑤ 报告落 per-virus 子目录)。"""
        try:
            hits = sorted(base.rglob(name)) if base.is_dir() else []
        except OSError:
            return None
        return hits[0] if hits else None

    cand = [
        # ① 报告 = generate_report.py 默认输出 (work_dir = 01_PublicData 根)
        ("public", "① Public Data Report",
         _root("public") / "Pipeline_Summary_Report.html"),
        ("discovery", "③ Discovery Report",
         _root("discovery") / dirs["d_reports"] / "pipeline_report.html"),
        ("analysis", "④ Analysis Report", _first_html(_root("analysis") / dirs["a_report"])),
        # ⑤ 报告名固定 phylo_report.html, 落 per-virus 子目录 (report_builder)
        ("phylo", "⑤ Phylogeny Report", _rglob_first(_root("phylo"), "phylo_report.html")),
        ("eve", "⑥ EVE Report", _root("eve") / dirs["e_summary"] / "eve_report.html"),
        # ⑦ 提交产物落项目根 submission_<run_title>/ (⑦ 无独立输出根)
        ("submission", "⑦ Submission Report", _glob_first(proj, "submission_*/report.html")),
        ("host_expression", "⑧ Host Expression Report",
         _first_html(_root("host_expression") / dirs["he_report"])),
        ("preprocessing", "② Preprocessing Report",
         _root("preprocessing") / "preprocessing_report.html"),
    ]
    out = []
    for pk, title, path in cand:
        if path is None:
            continue
        try:
            rel = path.resolve().relative_to(proj.resolve()).as_posix()
        except ValueError:
            rel = str(path)
        out.append({"pipeline": pk, "title": title, "rel": rel,
                    "exists": "yes" if path.is_file() else "no"})
    return out


def _timeline_svg(timeline: List[Dict[str, Any]]) -> str:
    """digest 演化图: 在线管线数随快照序列的阶梯折线 (纯内联 SVG, 零依赖)。

    少于 2 个快照时返回空串 (主报告自动省略)。"""
    pts_data = [s for s in timeline if isinstance(s.get("pipelines_present"), int)]
    if len(pts_data) < 2:
        return ""
    W, H = 640, 130
    pad_l, pad_r, pad_t, pad_b = 34, 100, 14, 24
    ys = [s["pipelines_present"] for s in pts_data]
    ymax = max(ys) or 1
    n = len(pts_data)

    def X(i: float) -> float:
        return pad_l + i * (W - pad_l - pad_r) / (n - 1)

    def Y(v: float) -> float:
        return H - pad_b - (v / ymax) * (H - pad_t - pad_b)

    grid = []
    step = max(1, ymax // 4)
    for gy in range(0, ymax + 1, step):
        yy = Y(gy)
        grid.append(f'<line x1="{pad_l}" y1="{yy:.1f}" x2="{W - pad_r}" y2="{yy:.1f}" '
                    f'stroke="{PALETTE["border"]}" stroke-width="1"/>'
                    f'<text x="{pad_l - 6}" y="{yy + 4:.1f}" text-anchor="end" '
                    f'font-size="10" fill="{PALETTE["ink-secondary"]}">{gy}</text>')
    poly = " ".join(f"{X(i):.1f},{Y(v):.1f}" for i, v in enumerate(ys))
    dots = "".join(
        f'<circle cx="{X(i):.1f}" cy="{Y(v):.1f}" r="3" fill="{PALETTE["primary"]}"/>'
        for i, v in enumerate(ys))
    t0 = _esc(pts_data[0].get("ts", ""))[5:16]
    t1 = _esc(pts_data[-1].get("ts", ""))[5:16]
    return (f'<svg viewBox="0 0 {W} {H}" style="width:100%;max-width:640px" '
            f'role="img" aria-label="digest timeline">{" ".join(grid)}'
            f'<polyline points="{poly}" fill="none" stroke="{PALETTE["primary"]}" '
            f'stroke-width="2"/>{dots}'
            f'<text x="{pad_l}" y="{H - 6}" font-size="10" fill="{PALETTE["ink-secondary"]}">{t0}</text>'
            f'<text x="{W - pad_r}" y="{H - 6}" text-anchor="end" font-size="10" '
            f'fill="{PALETTE["ink-secondary"]}">{t1}</text>'
            f'<text x="{W - pad_r + 8}" y="{Y(ys[-1]) + 4:.1f}" font-size="11" '
            f'font-weight="bold" fill="{PALETTE["primary"]}">{ys[-1]}</text></svg>')


def collect_process(project_dir: Path, layout: Optional[str] = None) -> Dict[str, Any]:
    """聚合全过程记录: 当前 digest + 时间线 + AI 参与 + agent 运行 + triage。"""
    project = Path(project_dir)
    state = run_digest.state_dir(project)
    digest = run_digest.build_digest(project, layout=layout)
    return {
        "project_dir": str(project),
        "generated_at": time.strftime("%Y-%m-%d %H:%M:%S"),
        "layout": digest.get("layout"),
        "digest": digest,
        "timeline": run_digest.load_history(project),
        "ai": _summarize_ai_log(state),
        "agent_runs": _summarize_agent_runs(state),
        "triage": collect_triage(project, layout),
        "pipeline_reports": collect_pipeline_reports(project, layout),
    }


def compact_for_ai(record: Dict[str, Any]) -> Dict[str, Any]:
    """过程记录瘦身版 (喂 process_narrative 提示词; 上限可控)。

    刻意剔除 generated_at 等易变字段: 同一项目两次采集的 prompt 相同 →
    叙事命中 ai_cache, 事后补解读/enrich 重跑零计费且文案稳定。
    """
    digest = record.get("digest", {})
    slim = {k: v for k, v in digest.items()
            if k not in ("pipelines", "generated_at")}
    slim["pipelines"] = {}
    for pk, sec in digest.get("pipelines", {}).items():
        stages = sec.get("stages", {})
        present_keys = [k for k, v in stages.items() if isinstance(v, dict) and v.get("exists")]
        slim["pipelines"][pk] = {"present": sec.get("present", False),
                                 "stages_done": present_keys}
        for k, v in sec.items():
            if k not in ("stages", "present") and not isinstance(v, (dict, list)):
                slim["pipelines"][pk][k] = v
    tl = record.get("timeline", [])
    slim["timeline_tail"] = tl[-5:]
    slim["ai"] = {k: v for k, v in record.get("ai", {}).items()
                  if k in ("n_calls", "by_op", "models", "total_tokens")}
    slim["agent_runs"] = [{"summary": r.get("summary"),
                           "stop_reason": r.get("stop_reason")}
                          for r in record.get("agent_runs", [])]
    return slim


# ═══════════════════════════════════════════════════════════════════
# 确定性过程面板 (零 AI 依赖)
# ═══════════════════════════════════════════════════════════════════

def _timeline_ctx(timeline: List[Dict[str, Any]], tail: int = 6) -> List[Dict[str, Any]]:
    """digest 时间线 → 模板行上下文 [{ts, n, on}] (取末尾 tail 条;
    on = 在线管线名逗号串。循环/条件渲染在 process_panel/master_report 模板)。"""
    out: List[Dict[str, Any]] = []
    for snap in (timeline or [])[-tail:]:
        on = [k for k, v in (snap.get("pipelines") or {}).items() if v]
        out.append({"ts": snap.get("ts", "?"),
                    "n": snap.get("pipelines_present", 0),
                    "on": ", ".join(on)})
    return out


def _rel_href(base_dir: Path, target: Path) -> str:
    """面板宿主报告 → 目标的相对链接 (posix 风格; 转义由模板 autoescape 负责)。"""
    import os
    return Path(os.path.relpath(target.resolve(),
                                Path(base_dir).resolve())).as_posix()


def build_process_panel_html(record: Dict[str, Any], *,
                             links_base: Optional[Path] = None,
                             project_root: Optional[Path] = None) -> str:
    """确定性过程面板 HTML (无 AI、无网络、幂等标记 PANEL_ID)。

    links_base/project_root 给定时追加"相关报告"互链行 (主报告 ↔ 各管线报告);
    链接只指向真实存在的文件 —— 不存在的报告绝不造链接。
    渲染: templates/process_panel.html.j2 (autoescape + StrictUndefined)。
    """
    digest = record.get("digest", {})
    ai = record.get("ai", {})
    triage_list = record.get("triage") or []

    chips: List[Dict[str, Any]] = []
    for pk in _DISPLAY_ORDER:
        card = stage_cards.get_card(pk)
        if card is None:
            continue
        sec = digest.get("pipelines", {}).get(pk, {})
        chips.append({"num": card["num"], "name": card["name"],
                      "on": bool(sec.get("present"))})

    links: List[Dict[str, str]] = []
    if links_base is not None and project_root is not None:
        master = Path(project_root) / run_digest.HANDOFF_DIR / MASTER_REPORT_NAME
        # 面板可能嵌在主报告里 — 自己链自己是噪声, 排除
        cur = Path(links_base).resolve()
        if master.is_file() and master.resolve().parent != cur:
            links.append({"href": _rel_href(links_base, master),
                          "label": "Master Report"})
        for rep in record.get("pipeline_reports") or []:
            if rep.get("exists") != "yes":
                continue
            p = Path(project_root) / rep["rel"]
            if p.resolve().parent == cur:
                continue   # 面板宿主报告自身不链
            links.append({"href": _rel_href(links_base, p),
                          "label": rep["title"]})

    return render(
        "process_panel.html.j2",
        panel_id=PANEL_ID,
        chips=chips,
        timeline_tail=_timeline_ctx(record.get("timeline", []), tail=6),
        ai=ai,
        ai_ops_str=" · ".join(f"{k}×{v}"
                              for k, v in ai.get("by_op", {}).items()),
        agents=record.get("agent_runs", []),
        triage_n=len(triage_list),
        triage_keep=sum(t["keep"] for t in triage_list),
        triage_review=sum(t["review"] for t in triage_list),
        links=links)


def append_process_panel_to_html(html_path: Path, start_dir: Path) -> Dict[str, Any]:
    """把过程面板嵌入已生成的报告 HTML (幂等; 无项目状态时原样不动)。"""
    html_path = Path(html_path)
    root = run_digest.find_state_root(Path(start_dir))
    if root is None:
        return {"status": "no_state", "written": False,
                "reason": "未找到 90_Handoff/ (独立运行?), 面板跳过"}
    try:
        content = html_path.read_text(encoding="utf-8")
    except OSError as e:
        return {"status": "error", "written": False, "reason": str(e)}
    if PANEL_ID in content:
        return {"status": "already", "written": False}
    if "</body>" not in content:
        return {"status": "error", "written": False, "reason": "无 </body>"}
    record = collect_process(root)
    panel = build_process_panel_html(record, links_base=html_path.parent,
                                     project_root=root)
    content = content.replace("</body>", panel + "\n</body>", 1)
    try:
        html_path.write_text(content, encoding="utf-8")
    except OSError as e:
        return {"status": "error", "written": False, "reason": str(e)}
    return {"status": "ok", "written": True, "project_root": str(root)}


# ═══════════════════════════════════════════════════════════════════
# AI 过程叙事 (离线降级)
# ═══════════════════════════════════════════════════════════════════

def process_ai_narrative(record: Dict[str, Any], *,
                         state_dir_path: Optional[Path] = None,
                         config: Any = None, client: Any = None) -> Dict[str, Any]:
    """过程记录 → AI 三段叙事。纪律: 数字只引用给定 JSON; 离线返回 unavailable。"""
    prompt_text, meta = render_prompt(
        "process_narrative",
        process_json=json.dumps(compact_for_ai(record), ensure_ascii=False)[:20000])
    res = chat([{"role": "user", "content": prompt_text}], config=config,
               client=client,
               cache_path=_ai_cache_file(state_dir_path))
    if res.status != "ok":
        return {"status": res.status, "text": "", "reason": res.reason}
    log_ai_event(state_dir_path, {
        "kind": "call", "op": "process_narrative", "model": res.model,
        "cached": res.cached, "prompt": meta.to_dict(), "usage": res.usage,
        "latency_s": round(res.latency_s, 2)})
    return {"status": "ok", "text": res.text.strip(),
            "meta": {**res.to_dict(), "prompt": meta.to_dict()}}


# ═══════════════════════════════════════════════════════════════════
# 项目级主报告 (渲染: templates/master_report.html.j2)
# ═══════════════════════════════════════════════════════════════════

def _pipeline_details_ctx(digest: Dict[str, Any]) -> List[Dict[str, Any]]:
    """管线全景卡上下文: 在线管线 → [{num, name, n_done, n_stages, extra}]。"""
    out: List[Dict[str, Any]] = []
    for pk in _DISPLAY_ORDER:
        card = stage_cards.get_card(pk)
        if card is None:
            continue
        sec = digest.get("pipelines", {}).get(pk, {})
        if not sec.get("present"):
            continue
        stages = sec.get("stages", {})
        done = [k for k, v in stages.items() if isinstance(v, dict) and v.get("exists")]
        extra = " · ".join(
            f"{k}={v}" for k, v in sec.items()
            if k not in ("stages", "present") and not isinstance(v, (dict, list))
            and v not in (None, 0, ""))
        out.append({"num": card["num"], "name": card["name"],
                    "n_done": len(done), "n_stages": len(stages), "extra": extra})
    return out


def generate_master_report(project_dir: Path, *,
                           layout: Optional[str] = None,
                           out_path: Optional[Path] = None,
                           want_ai: bool = True,
                           client: Any = None) -> Dict[str, Any]:
    """项目级主报告: 全景 + 过程 + AI 参与治理 + 建议 + (可选)AI 过程叙事。"""
    project = Path(project_dir)
    record = collect_process(project, layout=layout)
    digest = record["digest"]
    recs = advise(digest)
    state = run_digest.state_dir(project)

    narrative = {"status": "skipped"}
    if want_ai:
        narrative = process_ai_narrative(record, state_dir_path=state,
                                         client=client)
    out = Path(out_path) if out_path else project / run_digest.HANDOFF_DIR / MASTER_REPORT_NAME

    ai = record["ai"]
    n_pipelines = digest.get("totals", {}).get("pipelines_present", 0)
    recs_ctx = [{"level": r.level, "level_upper": r.level.upper(),
                 "title": r.title, "detail": r.detail, "command": r.command or "",
                 "color": _LEVEL_COLOR.get(r.level, PALETTE["info"])}
                for r in recs]

    html_doc = render(
        "master_report.html.j2",
        project_dir=record["project_dir"],
        layout=record["layout"],
        generated_at=record["generated_at"],
        nav_reps=[r for r in record.get("pipeline_reports", [])
                  if r.get("exists") == "yes"],
        kpi_pipelines=f"{n_pipelines}/{len(_DISPLAY_ORDER)}",
        kpi_ai_calls=ai["n_calls"],
        kpi_tokens=ai["total_tokens"],
        kpi_agent_runs=len(record["agent_runs"]),
        kpi_triage_review=sum(t["review"] for t in record.get("triage", [])),
        pipeline_details=_pipeline_details_ctx(digest),
        svg=_timeline_svg(record.get("timeline", [])),
        timeline_rows=_timeline_ctx(record.get("timeline", []), tail=12),
        ai=ai,
        triage=record.get("triage") or [],
        recs=recs_ctx,
        narrative_status=narrative.get("status"),
        narrative_text=narrative.get("text", ""),
        narrative_reason=narrative.get("reason", "AI not configured"),
        panel_html=build_process_panel_html(record, links_base=out.parent,
                                            project_root=project))

    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(html_doc, encoding="utf-8")
    return {"status": "ok", "path": str(out),
            "ai_narrative": narrative.get("status", "skipped"),
            "n_recommendations": len(recs)}


def enrich_project_reports(project_dir: Path, *, layout: Optional[str] = None,
                           want_ai: bool = True, include_master: bool = True,
                           client: Any = None) -> Dict[str, Any]:
    """事后补 AI: 对项目里【已生成】的报告 HTML 补插 AI 解读段 + 过程面板。

    典型场景: 无 key 时管线照常跑、报告照常出 (离线零影响); 之后配好 key
    再跑本函数, 五份管线报告与主报告一次性补齐 AI 内容 —— 无需重跑任何
    管线阶段, 无需重新生成报告。

    幂等: AI 段 (id=mmpv-ai-section) 与过程面板 (id=PANEL_ID) 已在的文件
    跳过, 可反复执行; 无 key 时 AI 段保持缺席 (文件不动), 面板照补。
    补入解读的数据来源 = 该项目 run digest (机器状态, 事后唯一可信源)。
    """
    from mmpv_common.ai import AIConfig, append_ai_section_to_html
    project = Path(project_dir)
    state = run_digest.state_dir(project)
    slim = run_digest.compact_for_ai(run_digest.build_digest(project, layout=layout))
    results: List[Dict[str, str]] = []
    for rep in collect_pipeline_reports(project, layout):
        if rep.get("exists") != "yes":
            continue
        html = project / rep["rel"]
        ai_status = "skipped"
        if want_ai:
            stats = {"pipeline": rep["pipeline"], "report": rep["title"],
                     "note": "post-hoc interpretation: based on the project run digest (machine state)",
                     "pipeline_digest": slim["pipelines"].get(rep["pipeline"], {}),
                     "project_totals": slim.get("totals", {})}
            r = append_ai_section_to_html(
                html, stats, state_dir=state, client=client,
                note="post-hoc interpretation · source: project run digest (90_Handoff/state/run_digest.json)")
            ai_status = r.get("status", "?")
        p = append_process_panel_to_html(html, project)
        results.append({"report": rep["rel"], "ai": ai_status,
                        "panel": p.get("status", "?")})
    master: Dict[str, Any] = {"status": "skipped"}
    if include_master:
        master = generate_master_report(project, layout=layout,
                                        want_ai=want_ai, client=client)
    return {"status": "ok", "project_dir": str(project),
            "reports": results,
            "master": {"path": master.get("path", ""),
                       "ai_narrative": master.get("ai_narrative", "skipped")}}
