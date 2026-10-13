# -*- coding: utf-8 -*-
"""ai.py 单测 — L2 核心: 提示词版本化 / 离线降级 / 溯源日志 / JSON 容错解析。

全部离线: 不发网络请求, LLM 用桩对象 (SimpleNamespace) 注入。
"""

import json
from types import SimpleNamespace

import pytest

from mmpv_common import ai
from mmpv_common.ai import (AIConfig, AIResult, ai_status, ask, chat,
                            log_ai_event, parse_json, register_prompt,
                            render_prompt, summarize_run)


# ── 桩: 模拟 openai client ──────────────────────────────────────────

def _stub_client(text='{"ok": true}', tokens=(10, 5)):
    resp = SimpleNamespace(
        choices=[SimpleNamespace(message=SimpleNamespace(content=text))],
        usage=SimpleNamespace(prompt_tokens=tokens[0], completion_tokens=tokens[1],
                              total_tokens=tokens[0] + tokens[1]))
    calls = []

    def create(**kwargs):
        calls.append(kwargs)
        return resp

    return SimpleNamespace(chat=SimpleNamespace(
        completions=SimpleNamespace(create=create))), calls


# ── 配置与离线降级 ──────────────────────────────────────────────────

def test_config_from_env_priority(tmp_path, monkeypatch):
    monkeypatch.setenv("MMPV_AI_API_KEY", "k-mmpv")
    monkeypatch.setenv("DEEPSEEK_API_KEY", "k-deepseek")
    cfg = AIConfig.from_env()
    assert cfg.api_key == "k-mmpv"          # MMPV_ 优先
    monkeypatch.delenv("MMPV_AI_API_KEY")
    assert AIConfig.from_env().api_key == "k-deepseek"   # 回退 DEEPSEEK (gsa_sra 同口径)


def test_ai_status_offline_reason(monkeypatch):
    monkeypatch.delenv("MMPV_AI_API_KEY", raising=False)
    monkeypatch.delenv("DEEPSEEK_API_KEY", raising=False)
    ok, why = ai_status(AIConfig(api_key=""))
    assert not ok
    assert "MMPV_AI_API_KEY" in why


def test_chat_offline_returns_unavailable(monkeypatch):
    monkeypatch.delenv("MMPV_AI_API_KEY", raising=False)
    monkeypatch.delenv("DEEPSEEK_API_KEY", raising=False)
    res = chat([{"role": "user", "content": "hi"}], config=AIConfig(api_key=""))
    assert res.status == "unavailable"
    assert res.reason


# ── 提示词注册表 ────────────────────────────────────────────────────

def test_render_prompt_versioning_and_hash():
    text, meta = render_prompt("run_summary", digest_json='{"a": 1}')
    assert "digest" in text and '{"a": 1}' in text
    assert meta.name == "run_summary" and meta.version
    assert len(meta.sha256) == 64
    text2, meta2 = render_prompt("run_summary", digest_json='{"a": 1}')
    assert meta2.sha256 == meta.sha256        # 同输入同哈希 → 可复现
    _, meta3 = render_prompt("run_summary", "1.0.0", digest_json="{}")
    assert meta3.version == "1.0.0"


def test_register_prompt_conflict_rejected():
    register_prompt("_t_conflict", "1.0.0", "A {x}")
    with pytest.raises(ValueError):
        register_prompt("_t_conflict", "1.0.0", "B {x}")
    register_prompt("_t_conflict", "1.0.0", "A {x}")  # 同内容重复登记允许
    register_prompt("_t_conflict", "2.0.0", "B {x}")  # 升版本允许


def test_render_unknown_prompt_raises():
    with pytest.raises(KeyError):
        render_prompt("_no_such_prompt_xyz")


def test_render_missing_field_strict_raises():
    with pytest.raises(KeyError):
        render_prompt("run_summary")   # digest_json 缺失 → strict 抛 KeyError


# ── 调用层 (桩注入) ─────────────────────────────────────────────────

def test_chat_with_stub_extracts_text_and_usage():
    client, calls = _stub_client("hello")
    cfg = AIConfig(api_key="k")
    res = chat([{"role": "user", "content": "hi"}], config=cfg, client=client)
    assert res.status == "ok" and res.text == "hello"
    assert res.usage["total_tokens"] == 15
    assert calls[0]["model"] == cfg.model


def test_chat_stub_error_becomes_error_result():
    class Boom:
        class chat:  # noqa: N801
            class completions:  # noqa: N801
                @staticmethod
                def create(**kw):
                    raise RuntimeError("network down")
    res = chat([{"role": "user", "content": "hi"}], config=AIConfig(api_key="k"),
               client=Boom())
    assert res.status == "error" and "network down" in res.reason


# ── JSON 容错解析 ───────────────────────────────────────────────────

def test_parse_json_variants():
    assert parse_json('{"a": 1}') == {"a": 1}
    assert parse_json('```json\n{"a": 1}\n```') == {"a": 1}
    assert parse_json('结果如下: {"a": [1,2]} 请查收') == {"a": [1, 2]}
    assert parse_json('前言 [1, 2] 后记') == [1, 2]
    with pytest.raises(ValueError):
        parse_json("完全不是 JSON")


# ── 溯源日志 ────────────────────────────────────────────────────────

def test_log_ai_event_jsonl(tmp_path):
    p1 = log_ai_event(tmp_path, {"kind": "call", "op": "t"})
    p2 = log_ai_event(tmp_path, {"kind": "call", "op": "t2"})
    assert p1 == p2 and p1.exists()
    lines = p1.read_text(encoding="utf-8").strip().splitlines()
    assert len(lines) == 2
    rec = json.loads(lines[0])
    assert rec["kind"] == "call" and "ts" in rec


def test_log_ai_event_none_is_noop():
    assert log_ai_event(None, {"kind": "call"}) is None


# ── 高层封装 ────────────────────────────────────────────────────────

def test_summarize_run_ok_and_logged(tmp_path):
    client, _ = _stub_client("## 运行状态\n一切正常")
    res = summarize_run({"totals": {"pipelines_present": 1}},
                        state_dir=tmp_path, client=client)
    assert res["status"] == "ok" and "运行状态" in res["summary"]
    assert res["meta"]["prompt"]["name"] == "run_summary"
    log = (tmp_path / "ai_log.jsonl").read_text(encoding="utf-8")
    assert "summarize_run" in log


def test_summarize_run_offline_graceful():
    res = summarize_run({}, client=None,
                        config=AIConfig(api_key=""))
    assert res["status"] == "unavailable" and res["summary"] == ""


def test_ai_interpretation_section_offline_empty_html():
    from mmpv_common.ai import ai_interpretation_section
    res = ai_interpretation_section({"stage": 1}, client=None,
                                    config=AIConfig(api_key=""))
    assert res["status"] == "unavailable" and res["html"] == ""


def test_ai_interpretation_section_ok(tmp_path):
    from mmpv_common.ai import ai_interpretation_section
    client, _ = _stub_client("解读文本")
    res = ai_interpretation_section({"stage": 1}, state_dir=tmp_path, client=client)
    assert res["status"] == "ok"
    assert "AI Interpretation" in res["html"] and "解读文本" in res["html"]


def test_ask_offline_and_ok():
    res = ask("问题", client=None, config=AIConfig(api_key=""))
    assert res["status"] == "unavailable"
    client, _ = _stub_client("答案")
    res2 = ask("问题", context="背景", client=client)
    assert res2["status"] == "ok" and res2["answer"] == "答案"


def test_ai_result_defaults():
    r = AIResult(status="ok")
    assert r.to_dict()["status"] == "ok"
