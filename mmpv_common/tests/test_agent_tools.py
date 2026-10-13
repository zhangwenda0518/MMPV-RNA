# -*- coding: utf-8 -*-
"""落地顺序四阶段收尾测试 — stage_tools / agent 受控循环 / 报告 AI 接线。

全部离线: LLM 桩注入, executor/confirm 注入, 不真跑编排器。
"""

import json
import shutil
from pathlib import Path
from types import SimpleNamespace

import pytest

from mmpv_common import run_digest, stage_tools
from mmpv_common.agent import (ST_BLOCKED, ST_FAILED, ST_PLANNED, ST_RAN,
                               ST_SKIPPED, ST_UNCONFIRMED, run_agent)
from mmpv_common.ai import AIConfig, append_ai_section_to_html

from test_ai_driven import make_project   # 复用假项目构造器

# 假项目的 advisor 建议链只含"提交"一步, 其命令带 <RUN_TITLE> 占位符
FILL = {"RUN_TITLE": "demo"}


@pytest.fixture()
def project(tmp_path):
    return make_project(tmp_path / "proj")


def _raw_only(parent: Path) -> Path:
    """只有原始数据的项目 (触发 preprocessing 建议)。"""
    root = parent / "raw_only"
    shutil.rmtree(root, ignore_errors=True)
    raw = root / "01_PublicData" / "03_RawData"
    raw.mkdir(parents=True)
    (raw / "s.fastq.gz").write_bytes(b"x")
    return root


def _exec(code_book: dict, calls: list):
    """按命令片段返回退出码的桩 executor。"""
    def _run(tokens, **kw):
        joined = " ".join(tokens)
        calls.append(joined)
        for frag, code in code_book.items():
            if frag in joined:
                return SimpleNamespace(returncode=code)
        return SimpleNamespace(returncode=0)
    return _run


# ═══════════════════════════════════════════════════════════════════
# stages 工具化
# ═══════════════════════════════════════════════════════════════════

def test_tools_json_covers_pipelines():
    tools = stage_tools.tools_json()
    names = {t["name"] for t in tools}
    assert {"mmpv_run_public", "mmpv_run_hostref", "mmpv_run_preprocessing",
            "mmpv_run_discovery", "mmpv_run_analysis", "mmpv_run_phylo",
            "mmpv_run_eve", "mmpv_run_submission"} <= names
    for t in tools:
        assert t["inputSchema"]["type"] == "object"
        assert "stage" in t["inputSchema"]["properties"]
        assert "退出码" in t["description"]


def test_preconditions_and_postconditions(project):
    d = run_digest.build_digest(project)
    assert stage_tools.check_preconditions("public", d) == []
    assert stage_tools.check_preconditions("discovery", d) == []   # 假项目有去宿主
    assert stage_tools.check_preconditions("submission", d) == []  # 有 10_Rescue
    empty = run_digest.build_digest(project.parent / "nope")
    assert any("去宿主" in m or "原始数据" in m
               for m in stage_tools.check_preconditions("discovery", empty))
    assert stage_tools.check_postconditions("discovery", d) == []


def test_plan_chain_order_and_substitution(project):
    d = run_digest.build_digest(project)
    steps = stage_tools.plan_chain(d)
    assert steps, "假项目 (有 rescue 无提交) 应产出提交建议步"
    orders = [s.order for s in steps]
    assert orders == sorted(orders)                     # 主链序
    for s in steps:
        assert "<PROJECT_ROOT>" not in s.command            # 项目根已自动替换
        assert str(project) in s.command or s.missing_inputs


def test_plan_chain_missing_inputs_not_invented(project):
    """占位符必须显式进 missing_inputs, 不臆造路径。"""
    raw_root = _raw_only(project.parent)
    d = run_digest.build_digest(raw_root)
    steps = stage_tools.plan_chain(d, pipelines=["preprocessing"])
    assert steps and steps[0].tool == "preprocessing"
    assert any("RAW_FASTQ_DIR" in m for m in steps[0].missing_inputs)
    assert "<RAW_FASTQ_DIR>" in steps[0].command


def test_plan_chain_pipeline_filter_and_inputs(project):
    d = run_digest.build_digest(project)
    steps = stage_tools.plan_chain(d, pipelines=["submission"],
                                   inputs={"RUN_TITLE": "tomato"})
    assert all(s.tool == "submission" for s in steps)
    assert all("tomato" in s.command or s.missing_inputs for s in steps)


# ═══════════════════════════════════════════════════════════════════
# agent 受控循环 (executor/confirm 注入)
# ═══════════════════════════════════════════════════════════════════

def test_agent_dry_run_no_execution(project, tmp_path):
    calls = []
    rep = run_agent(project, inputs=FILL, dry_run=True, yes=True,
                    executor=_exec({}, calls), confirm_fn=lambda p: False)
    assert calls == []                                  # dry-run 不执行
    assert rep.steps and all(s.status == ST_PLANNED for s in rep.steps)
    assert rep.journal_path and Path(rep.journal_path).exists()


def test_agent_ran_and_journal(project, tmp_path):
    calls = []
    rep = run_agent(project, inputs=FILL, yes=True,
                    executor=_exec({}, calls), confirm_fn=lambda p: True)
    assert rep.steps
    ran = [s for s in rep.steps if s.status == ST_RAN]
    assert ran, [s.status for s in rep.steps]
    assert rep.journal_path
    body = json.loads(Path(rep.journal_path).read_text(encoding="utf-8"))
    assert body["steps"]


def test_agent_exit1_stops_chain(project, tmp_path):
    calls = []
    book = {"submission_pipeline.py": 1}   # 提交步必败
    rep = run_agent(project, inputs=FILL, yes=True,
                    executor=_exec(book, calls), confirm_fn=lambda p: True)
    assert rep.stop_reason == "failed"
    failed = [s for s in rep.steps if s.status == ST_FAILED]
    assert failed and failed[0].exit_code == 1
    idx_fail = rep.steps.index(failed[0])
    assert all(s.status != ST_RAN for s in rep.steps[idx_fail + 1:])


def test_agent_exit2_blocked_continues(project, tmp_path):
    calls = []
    book = {"submission_pipeline.py": 2}
    rep = run_agent(project, inputs=FILL, yes=True,
                    executor=_exec(book, calls), confirm_fn=lambda p: True)
    blocked = [s for s in rep.steps if s.status == ST_BLOCKED]
    assert blocked and blocked[0].exit_code == 2
    assert rep.stop_reason != "failed"


def test_agent_skipped_when_inputs_missing(project, tmp_path):
    """不填 <RUN_TITLE> → 该步显式 skipped, 绝不臆造路径执行。"""
    calls = []
    rep = run_agent(project, yes=True, executor=_exec({}, calls),
                    confirm_fn=lambda p: True)
    skipped = [s for s in rep.steps if s.status == ST_SKIPPED]
    assert skipped and any("RUN_TITLE" in s.detail for s in skipped)
    assert calls == []


def test_agent_dangerous_never_auto_confirmed(project, tmp_path):
    """红线: 危险命令即使 --yes, 非交互 (confirm_fn=None) 也必须跳过。"""
    from mmpv_common.copilot import validate
    v = validate("python endogenous_virus_pipeline/eve_screen.py -g g.fa -n sp "
                 "-o out/ --force")
    assert v.requires_confirmation
    calls = []
    rep = run_agent(project, inputs=FILL, yes=True,
                    executor=_exec({}, calls), confirm_fn=None)
    for s in rep.steps:
        if s.status == ST_RAN:
            assert "--force" not in s.command


def test_agent_unconfirmed_when_user_says_no(project, tmp_path):
    calls = []
    rep = run_agent(project, inputs=FILL, yes=False,
                    executor=_exec({}, calls), confirm_fn=lambda p: False)
    assert all(s.status != ST_RAN for s in rep.steps)   # 全被拒
    assert any(s.status == ST_UNCONFIRMED for s in rep.steps)


def test_agent_max_steps_budget(project, tmp_path):
    calls = []
    rep = run_agent(project, inputs=FILL, yes=True, max_steps=1,
                    executor=_exec({}, calls), confirm_fn=lambda p: True)
    ran = [s for s in rep.steps if s.status == ST_RAN]
    assert len(ran) <= 1
    assert rep.stop_reason == "max_steps=1 budget exhausted" or ran


def test_agent_empty_plan(tmp_path):
    rep = run_agent(tmp_path / "nothing", dry_run=True)
    assert rep.steps == [] and "empty plan" in rep.stop_reason


# ═══════════════════════════════════════════════════════════════════
# 报告 AI 接线 (阶段一尾)
# ═══════════════════════════════════════════════════════════════════

def _stub_client(text="解读内容"):
    resp = SimpleNamespace(
        choices=[SimpleNamespace(message=SimpleNamespace(content=text))],
        usage=None)

    class Stub:
        class chat:  # noqa: N801
            class completions:  # noqa: N801
                @staticmethod
                def create(**kw):
                    return resp
    return Stub()


def test_append_ai_section_offline_no_change(tmp_path):
    html = tmp_path / "r.html"
    html.write_text("<html><body>报告本体</body></html>", encoding="utf-8")
    before = html.read_text(encoding="utf-8")
    res = append_ai_section_to_html(html, {"pipeline": "t"},
                                    config=AIConfig(api_key=""))
    assert res["status"] == "unavailable" and not res["written"]
    assert html.read_text(encoding="utf-8") == before   # 离线文件原样


def test_append_ai_section_inserts_before_body(tmp_path):
    html = tmp_path / "r.html"
    html.write_text("<html><body>报告本体</body></html>", encoding="utf-8")
    res = append_ai_section_to_html(html, {"pipeline": "t"},
                                    client=_stub_client("AI 说: 一切正常"))
    assert res["status"] == "ok" and res["written"]
    out = html.read_text(encoding="utf-8")
    assert "报告本体" in out and "AI 说: 一切正常" in out
    assert out.index("报告本体") < out.index("AI Interpretation") < out.rindex("</body>")


def test_append_ai_section_no_body_tag(tmp_path):
    html = tmp_path / "r.html"
    html.write_text("<html>残缺", encoding="utf-8")
    before = html.read_text(encoding="utf-8")
    res = append_ai_section_to_html(html, {"pipeline": "t"},
                                    client=_stub_client())
    assert res["status"] == "error" and not res["written"]
    assert html.read_text(encoding="utf-8") == before
