# -*- coding: utf-8 -*-
"""事后补 AI (enrich) 测试 — 幂等 AI 段 / 离线只补面板 / 桩注入补齐 / 缓存稳定。"""

from pathlib import Path
from types import SimpleNamespace

import pytest

from mmpv_common import run_digest
from mmpv_common.ai import AIConfig, append_ai_section_to_html
from mmpv_common.process_report import (compact_for_ai, collect_process,
                                        enrich_project_reports)

from test_ai_driven import make_project


@pytest.fixture()
def project(tmp_path):
    """假项目 + 三份"离线时代生成"的管线报告 (无 AI 段、无面板)。"""
    root = make_project(tmp_path / "proj")
    d = root / "03_Discovery" / "13_Report"
    d.mkdir(parents=True, exist_ok=True)
    (d / "pipeline_report.html").write_text(
        "<html><body><h1>发现报告</h1></body></html>", encoding="utf-8")
    (root / "02_Preprocessing" / "preprocessing_report.html").write_text(
        "<html><body><h1>清洗报告</h1></body></html>", encoding="utf-8")
    # make_project 的 ④ 报告是残缺 <html/> (无 </body>), 换成规整报告
    (root / "04_Analysis" / "09_Report" / "report.html").write_text(
        "<html><body><h1>分析报告</h1></body></html>", encoding="utf-8")
    return root


def _stub_client(text="事后补的 AI 解读内容"):
    resp = SimpleNamespace(
        choices=[SimpleNamespace(message=SimpleNamespace(content=text))],
        usage=SimpleNamespace(prompt_tokens=1, completion_tokens=1, total_tokens=2))

    class Stub:
        class chat:  # noqa: N801
            class completions:  # noqa: N801
                @staticmethod
                def create(**kw):
                    return resp
    return Stub()


# ── AI 段幂等 ───────────────────────────────────────────────────────

def test_append_ai_section_idempotent(tmp_path):
    html = tmp_path / "r.html"
    html.write_text("<html><body>x</body></html>", encoding="utf-8")
    r1 = append_ai_section_to_html(html, {"pipeline": "t"}, client=_stub_client())
    assert r1["status"] == "ok" and r1["written"]
    r2 = append_ai_section_to_html(html, {"pipeline": "t"}, client=_stub_client())
    assert r2["status"] == "already" and not r2["written"]
    body = html.read_text(encoding="utf-8")
    assert body.count('id="mmpv-ai-section"') == 1


def test_append_ai_section_checks_file_before_calling_ai(tmp_path):
    """已插过的文件不应再触发一次 AI 调用 (省 token)。"""
    calls = []

    class Counting:
        class chat:  # noqa: N801
            class completions:  # noqa: N801
                @staticmethod
                def create(**kw):
                    calls.append(1)
                    return SimpleNamespace(
                        choices=[SimpleNamespace(message=SimpleNamespace(content="x"))],
                        usage=None)
    html = tmp_path / "r.html"
    html.write_text("<html><body>x</body></html>", encoding="utf-8")
    append_ai_section_to_html(html, {"pipeline": "t"}, client=Counting())
    append_ai_section_to_html(html, {"pipeline": "t"}, client=Counting())
    assert len(calls) == 1


# ── 缓存稳定 (事后补叙事不重复计费的前提) ──────────────────────────

def test_compact_for_ai_drops_volatile_fields(project):
    r1 = collect_process(project)
    r2 = collect_process(project)          # generated_at 必不同
    c1, c2 = compact_for_ai(r1), compact_for_ai(r2)
    assert "generated_at" not in c1
    assert c1 == c2                          # 同项目两次采集 → 同 prompt 输入


# ── enrich 三种形态 ────────────────────────────────────────────────

def test_enrich_offline_panels_only(project):
    res = enrich_project_reports(project, want_ai=True)   # 无 key
    by = {r["report"]: r for r in res["reports"]}
    assert by["03_Discovery/13_Report/pipeline_report.html"]["panel"] == "ok"
    assert by["03_Discovery/13_Report/pipeline_report.html"]["ai"] == "unavailable"
    body = (project / "03_Discovery" / "13_Report" / "pipeline_report.html").read_text(encoding="utf-8")
    assert "发现报告" in body                      # 报告本体还在
    assert "mmpv-process-panel" in body           # 面板补上了
    assert "mmpv-ai-section" not in body          # AI 段如实缺席
    assert res["master"]["ai_narrative"] == "unavailable"


def test_enrich_with_stub_injects_everything(project):
    res = enrich_project_reports(project, want_ai=True, client=_stub_client())
    for r in res["reports"]:
        assert r["ai"] == "ok" and r["panel"] == "ok"
    body = (project / "02_Preprocessing" / "preprocessing_report.html").read_text(encoding="utf-8")
    assert "事后补的 AI 解读内容" in body
    assert "run digest" in body                   # 注明解读来源
    assert res["master"]["ai_narrative"] == "ok"
    master = (project / "90_Handoff" / "MMPV_master_report.html").read_text(encoding="utf-8")
    assert "AI Process Narrative" in master


def test_enrich_idempotent_second_run(project):
    enrich_project_reports(project, want_ai=True, client=_stub_client())
    res2 = enrich_project_reports(project, want_ai=True, client=_stub_client())
    for r in res2["reports"]:
        assert r["ai"] == "already" and r["panel"] == "already"
    body = (project / "03_Discovery" / "13_Report" / "pipeline_report.html").read_text(encoding="utf-8")
    assert body.count('id="mmpv-ai-section"') == 1
    assert body.count('id="mmpv-process-panel"') == 1


def test_enrich_no_reports(tmp_path):
    root = make_project(tmp_path / "bare")        # 无报告文件
    (root / "04_Analysis" / "09_Report" / "report.html").unlink()
    res = enrich_project_reports(root, want_ai=False, include_master=True)
    assert res["reports"] == []
    assert Path(res["master"]["path"]).exists()   # 主报告仍会生成
