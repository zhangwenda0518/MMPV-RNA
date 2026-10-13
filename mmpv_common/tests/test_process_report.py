# -*- coding: utf-8 -*-
"""过程报告层测试 — digest 时间线 / AI 参与聚合 / 过程面板 / 主报告 / AI 缓存。

全部离线: AI 与 client 均桩注入或明确未配置。
"""

import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from mmpv_common import process_report, run_digest
from mmpv_common.ai import AIConfig, chat
from mmpv_common.process_report import (PANEL_ID, append_process_panel_to_html,
                                        build_process_panel_html,
                                        collect_process, generate_master_report)

from test_ai_driven import make_project


@pytest.fixture()
def project(tmp_path):
    """假项目 + 过程数据 (digest 历史 / ai_log / agent_runs)。"""
    root = make_project(tmp_path / "proj")
    # 时间线: 两次 digest 落盘
    run_digest.write_digest(run_digest.build_digest(root), root)
    run_digest.write_digest(run_digest.build_digest(root), root)
    state = run_digest.state_dir(root)
    # AI 参与记录
    state.joinpath("ai_log.jsonl").write_text(
        "\n".join([
            json.dumps({"kind": "call", "op": "ai_interpretation", "model": "stub-m",
                        "cached": False, "ts": "2026-10-08 01:00:00",
                        "prompt": {"name": "ai_interpretation", "version": "1.0.0"},
                        "usage": {"total_tokens": 120}}),
            json.dumps({"kind": "call", "op": "process_narrative", "model": "stub-m",
                        "cached": True, "ts": "2026-10-08 02:00:00",
                        "prompt": {"name": "process_narrative", "version": "1.0.0"},
                        "usage": {"total_tokens": 80}}),
            json.dumps({"kind": "agent", "op": "run_chain", "ts": "2026-10-08 02:30:00",
                        "summary": "执行 1 · 失败 0 · 跳过 0 · 共 1 步", "n_steps": 1}),
        ]) + "\n", encoding="utf-8")
    # agent 留痕
    ar = state / "agent_runs"
    ar.mkdir(parents=True, exist_ok=True)
    (ar / "20260101_000000.json").write_text(json.dumps({
        "finished_at": "2026-10-08 02:30:00",
        "summary": "执行 1 · 失败 0 · 跳过 0 · 共 1 步", "stop_reason": ""}),
        encoding="utf-8")
    return root


# ── digest 时间线 ──────────────────────────────────────────────────

def test_digest_history_appended(project):
    hist = run_digest.load_history(project)
    assert len(hist) >= 2
    assert hist[-1]["pipelines_present"] >= 5
    assert hist[-1]["pipelines"]["preprocessing"] is True


def test_find_state_root(project):
    assert run_digest.find_state_root(project / "03_Discovery") == project.resolve()
    assert run_digest.find_state_root(project / "03_Discovery" / "01_Assembly") == project.resolve()
    assert run_digest.find_state_root(project.parent) is None  # 项目外 → 不臆造


# ── 过程记录聚合 ────────────────────────────────────────────────────

def test_collect_process_aggregates(project):
    rec = collect_process(project)
    assert rec["ai"]["n_calls"] == 2
    assert rec["ai"]["by_op"]["ai_interpretation"] == 1
    assert rec["ai"]["n_cached"] == 1
    assert rec["ai"]["total_tokens"] == 200
    assert rec["ai"]["models"] == ["stub-m"]
    assert "process_narrative@1.0.0" in rec["ai"]["prompt_versions"]
    assert len(rec["agent_runs"]) == 1
    assert len(rec["timeline"]) >= 2
    assert rec["digest"]["totals"]["pipelines_present"] >= 5


def test_collect_process_empty_state(tmp_path):
    root = make_project(tmp_path / "bare")     # 无 90_Handoff/state 内容
    rec = collect_process(root)
    assert rec["ai"]["n_calls"] == 0
    assert rec["agent_runs"] == []


# ── 确定性过程面板 ──────────────────────────────────────────────────

def test_panel_html_deterministic_content(project):
    rec = collect_process(project)
    panel = build_process_panel_html(rec)
    assert PANEL_ID in panel
    assert "Project Process Overview" in panel and "AI Involvement Trace" in panel
    assert "AI calls: <b>2</b>" in panel
    assert "No agent runs recorded" not in panel        # 有 agent 记录
    assert "Governance statement" in panel              # 审稿友好声明恒在
    assert "①" in panel and "⑦" in panel              # 管线 chips


def test_panel_no_ai_records(tmp_path):
    root = make_project(tmp_path / "bare")
    panel = build_process_panel_html(collect_process(root))
    assert "No AI calls recorded" in panel
    assert "No timeline yet" in panel


def test_append_panel_idempotent(project):
    html = project / "03_Discovery" / "13_Report" / "r.html"
    html.parent.mkdir(parents=True, exist_ok=True)
    html.write_text("<html><body>报告</body></html>", encoding="utf-8")
    r1 = append_process_panel_to_html(html, project / "03_Discovery")
    assert r1["status"] == "ok" and r1["written"]
    body = html.read_text(encoding="utf-8")
    assert PANEL_ID in body and body.index("报告") < body.index(PANEL_ID)
    r2 = append_process_panel_to_html(html, project / "03_Discovery")
    assert r2["status"] == "already" and not r2["written"]
    assert html.read_text(encoding="utf-8").count(PANEL_ID) == 1


def test_append_panel_no_state_leaves_file_alone(tmp_path):
    root = tmp_path / "lonely"
    root.mkdir()
    html = root / "r.html"
    html.write_text("<html><body>x</body></html>", encoding="utf-8")
    res = append_process_panel_to_html(html, root)
    assert res["status"] == "no_state"
    assert html.read_text(encoding="utf-8") == "<html><body>x</body></html>"


# ── 项目级主报告 ────────────────────────────────────────────────────

def test_master_report_offline(project):
    res = generate_master_report(project, want_ai=True)   # 无 key → 叙事降级
    out = Path(res["path"])
    assert out.exists() and out.name == "MMPV_master_report.html"
    assert out.parent.name == "90_Handoff"
    body = out.read_text(encoding="utf-8")
    assert "Results × Process × AI Involvement" in body
    assert "Pipeline Overview" in body and "Process Timeline" in body
    assert "Current Recommendations" in body   # advisor 规则建议
    assert PANEL_ID in body                   # 过程面板
    assert res["ai_narrative"] in ("unavailable", "error")
    assert "unaffected" in body               # 离线降级说明文案


def test_master_report_ai_narrative_stub(project):
    resp = SimpleNamespace(
        choices=[SimpleNamespace(message=SimpleNamespace(
            content="项目已推进到 ③ 发现与 ④ 分析。过程中 AI 参与解读 1 次。建议下一步跑 ⑦ 提交。"))],
        usage=SimpleNamespace(prompt_tokens=10, completion_tokens=5, total_tokens=15))

    class Stub:
        class chat:  # noqa: N801
            class completions:  # noqa: N801
                @staticmethod
                def create(**kw):
                    return resp

    res = generate_master_report(project, want_ai=True, client=Stub())
    assert res["ai_narrative"] == "ok"
    body = Path(res["path"]).read_text(encoding="utf-8")
    assert "AI Process Narrative" in body and "发现与 ④ 分析" in body
    # 叙事调用也入溯源日志
    log = run_digest.state_dir(project).joinpath("ai_log.jsonl").read_text(encoding="utf-8")
    assert "process_narrative" in log


# ── AI 缓存 ─────────────────────────────────────────────────────────

def _counting_client(text, counter):
    resp = SimpleNamespace(
        choices=[SimpleNamespace(message=SimpleNamespace(content=text))],
        usage=SimpleNamespace(prompt_tokens=1, completion_tokens=1, total_tokens=2))

    class Stub:
        class chat:  # noqa: N801
            class completions:  # noqa: N801
                @staticmethod
                def create(**kw):
                    counter.append(1)
                    return resp
    return Stub()


def test_chat_cache_hit_and_miss(tmp_path, monkeypatch):
    monkeypatch.delenv("MMPV_AI_CACHE", raising=False)   # 默认开启
    cache = tmp_path / "ai_cache.json"
    cfg = AIConfig(api_key="k")
    c1 = _counting_client("答案A", [])
    msgs = [{"role": "user", "content": "同一问题"}]
    r1 = chat(msgs, config=cfg, client=c1, cache_path=cache)
    assert r1.status == "ok" and not r1.cached
    # 第二次同输入: 换一个新桩, 若缓存生效则桩不会被调用
    c2 = _counting_client("答案B", [])
    r2 = chat(msgs, config=cfg, client=c2, cache_path=cache)
    assert r2.status == "ok" and r2.cached and r2.text == "答案A"
    # 不同输入 → miss, 走桩
    c3 = _counting_client("答案C", [])
    r3 = chat([{"role": "user", "content": "另一个问题"}], config=cfg,
              client=c3, cache_path=cache)
    assert r3.status == "ok" and not r3.cached and r3.text == "答案C"


def test_chat_cache_disabled_by_env(tmp_path, monkeypatch):
    monkeypatch.setenv("MMPV_AI_CACHE", "0")
    from mmpv_common.ai import _ai_cache_file
    assert _ai_cache_file(tmp_path) is None
    monkeypatch.delenv("MMPV_AI_CACHE")
    assert _ai_cache_file(tmp_path) == tmp_path / "ai_cache.json"
    assert _ai_cache_file(None) is None


def test_cache_key_salted_by_api_key(tmp_path):
    """P2-29e: 缓存键掺 api_key 盐 (只入 8 位哈希) —— 换 key 不命中旧缓存。"""
    from mmpv_common.ai import _cache_key
    msgs = [{"role": "user", "content": "同一问题"}]
    k1 = _cache_key(AIConfig(api_key="key-one"), msgs, 0.2, False)
    assert k1 == _cache_key(AIConfig(api_key="key-one"), msgs, 0.2, False)   # 同 key 稳定
    assert k1 != _cache_key(AIConfig(api_key="key-two"), msgs, 0.2, False)   # 换 key 变键
    assert k1 != _cache_key(AIConfig(api_key=""), msgs, 0.2, False)          # 离线键独立
    # 端到端: 同一 cache 文件, 换 key 后不命中 (重新走桩)
    cache = tmp_path / "ai_cache.json"
    r1 = chat(msgs, config=AIConfig(api_key="k1"),
              client=_counting_client("答案A", []), cache_path=cache)
    r2 = chat(msgs, config=AIConfig(api_key="k2"),
              client=_counting_client("答案B", []), cache_path=cache)
    assert r1.status == "ok" and not r1.cached
    assert r2.status == "ok" and not r2.cached and r2.text == "答案B"
