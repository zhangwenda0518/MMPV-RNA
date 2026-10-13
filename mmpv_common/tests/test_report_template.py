#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""test_report_template.py — 报告模板层 (Jinja2) 守门测试
========================================================

覆盖 (对应报告层模板化的安全基线):
  1. 共享模板文件真实存在、全部可编译 (审计脚本不扫 .j2, 由本测试兜底);
  2. StrictUndefined: 缺上下文变量直接抛错, 不产出半成品报告;
  3. autoescape: 元数据值里的 HTML 不再破坏标记;
  4. 离线自包含: 渲染产物无外链资源引用 (src=/href= 指向 http);
  5. report_style 组件包装与宏输出一致 (单一实现来源);
  6. 面板/主报告/提交报告模板的回归契约子串。

跑法: python -m pytest mmpv_common/tests/test_report_template.py -q
"""

import re

import pytest
from jinja2 import UndefinedError

from mmpv_common import report_template as rt
from mmpv_common.report_style import badge, data_table, kpi_card, section_header

EXPECTED_TEMPLATES = {
    "base.html.j2", "macros.html.j2", "process_panel.html.j2",
    "master_report.html.j2", "submission_report.html.j2",
}

# 离线自包含: 资源型外链 (纯文本提及 URL 不算)
_EXT_REF_RE = re.compile(r"""(?:src|href)\s*=\s*["'](?:https?:)?//""", re.IGNORECASE)


# ── 1. 模板资产 ─────────────────────────────────────────────────────────

def test_template_files_present():
    names = set(rt.list_templates())
    missing = EXPECTED_TEMPLATES - names
    assert not missing, f"共享模板缺失: {missing} (mmpv_common/templates/)"


def test_all_templates_compile():
    env = rt.get_env()
    for name in rt.list_templates():
        env.get_template(name)  # 语法错误直接抛


# ── 2/3. 安全基线 ───────────────────────────────────────────────────────

def test_strict_undefined_fails_loud():
    with pytest.raises(UndefinedError):
        rt.render("process_panel.html.j2")  # 缺全部上下文 → 必须抛错


def test_autoescape_escapes_values():
    html = data_table(["col"], [["<script>alert(1)</script>"]])
    assert "<script>" not in html
    assert "&lt;script&gt;" in html


def test_base_css_not_escaped():
    """base_css 是可信静态样式 (内含 CSS 单引号), 不得被 autoescape 打成 &#39;。"""
    html = rt.render("base.html.j2")
    style = html.split("<style>", 1)[1].split("</style>", 1)[0]
    assert "&#39;" not in style, "base_css 被转义, 字体名将失效"
    assert "--primary:" in style


def test_base_render_offline():
    assert not _EXT_REF_RE.search(rt.render("base.html.j2")), "骨架模板出现外链资源"


# ── 4. 组件: 包装与宏单一来源 ──────────────────────────────────────────

def test_style_wrappers_match_macros():
    mod = rt.get_env().get_template("macros.html.j2").module
    assert kpi_card("检测数", 47, accent="success") == str(
        mod.kpi_card(label="检测数", value=47, accent="success"))
    assert data_table(["a"], [["1"]], table_id="tt") == str(
        mod.data_table(headers=["a"], rows=[["1"]], table_id="tt"))


def test_component_class_contract():
    assert 'class="kpi-card success"' in kpi_card("标签", 47, accent="success")
    assert 'class="kpi-value"' in kpi_card("标签", 47)
    assert 'class="data-table" id="tt"' in data_table(["a"], [["1"]], table_id="tt")
    assert 'class="badge badge-success"' in badge("PASS", "success")
    assert "section-header" in section_header("概览", icon="X")


# ── 5. 面板/主报告/提交报告回归契约 ────────────────────────────────────

def _minimal_record():
    return {
        "digest": {"pipelines": {"public": {"present": True},
                                 "submission": {"present": False}}},
        "ai": {"n_calls": 2, "by_op": {"interpret": 2},
               "by_op_tokens": {"interpret": 1200}, "total_tokens": 1200,
               "n_cached": 0, "models": ["test-model"], "prompt_versions": [],
               "by_month": {"2026-09": {"calls": 2, "tokens": 1200}}},
        "agent_runs": [],
        "triage": [],
        "timeline": [{"ts": "2026-09-01T10:00:00", "pipelines_present": 1,
                      "pipelines": {"public": True}}],
        "pipeline_reports": [],
    }


def test_process_panel_contract():
    from mmpv_common.process_report import PANEL_ID, build_process_panel_html
    html = build_process_panel_html(_minimal_record())
    assert f'id="{PANEL_ID}"' in html
    assert "Project Process Overview · AI Involvement Trace" in html
    assert "AI calls: <b>2</b>" in html            # 字面 <b>, autoescape 不得转义
    assert "tokens ~1,200" in html                  # fmt_int 千分位
    assert "①" in html and "⑦" in html              # 管线 chips
    assert "Governance statement" in html
    assert not _EXT_REF_RE.search(html)


def test_master_report_template_contract():
    ctx = {
        "project_dir": "D:/proj", "layout": "standard",
        "generated_at": "2026-10-08 12:00:00", "nav_reps": [],
        "kpi_pipelines": "1/8", "kpi_ai_calls": 2, "kpi_tokens": 1200,
        "kpi_agent_runs": 0, "kpi_triage_review": 0,
        "pipeline_details": [{"num": "①", "name": "公共数据", "n_done": 3,
                              "n_stages": 5, "extra": "reads=1200"}],
        "svg": "<svg viewBox='0 0 10 10'></svg>",
        "timeline_rows": [], "ai": _minimal_record()["ai"], "triage": [],
        "recs": [], "narrative_status": "unavailable",
        "narrative_text": "", "narrative_reason": "未配置 AI",
        "panel_html": "<section id='x'>panel</section>",
    }
    html = rt.render("master_report.html.j2", **ctx)
    assert "MMPV Master Report — Results × Process × AI Involvement" in html
    assert "Pipeline Overview" in html and "AI Usage" in html
    assert "unaffected" in html                        # 离线降级文案
    assert "Operation (by op)" in html and "interpret" in html
    assert "1,200" in html                             # tokens 千分位
    assert "<section id='x'>panel</section>" in html   # panel 片段不转义注入
    assert not _EXT_REF_RE.search(html)


def test_submission_report_template_contract():
    ctx = {
        "run_title": "demo_run", "timestamp": "2026-10-08 12:00:00",
        "pipeline": "MMPV-RNA v2.3 + suvtk v0.1.1", "n_seqs": 3, "n_viruses": 2,
        "suvtk_parts": ["topology: 2 linear / 1 circular"],
        "src_headers": ["SeqID", "organism"],
        "src_rows": [["v1", "Betaflexiviridae sp."], ["v2", "<uncultured>"]],
        "cmt_headers": ["Parameter", "Value"],
        "cmt_rows": [["Assembly_Method", "MEGAHIT;1.2.9"]],
        "bs_headers": ["sample_name", "bioproject"],
        "bs_rows": [["v1", "PRJNA123456"]],
        "miuvig_rows": [["miuvig_completeness", "high"]],
        "asm_rows": [["Assembly Method", "MEGAHIT v1.2.9"]],
        "topo_rows": [{"contig": "v1", "taxonomy": "Viruses",
                       "rule_topology": "linear", "final_topology": "linear",
                       "cls": "filled", "seq_length": "8123", "dtr": "✓",
                       "evidence": "DTR"}],
        "topo_summary": "1 seqs: 1 linear / 0 circular",
        "checklist": [("done", "suvtk taxonomy — ICTV classification complete"),
                      ("todo", "Register BioProject")],
        "js": "function downloadTSV(){}",
    }
    html = rt.render("submission_report.html.j2", **ctx)
    assert "GenBank Submission Report — demo_run" in html
    for tid in ("srcTable", "cmtTable", "bsTable", "miuvigTable",
                "topoTable", "asmTable", "completeBar"):
        assert tid in html, f"交互 JS 依赖的 id={tid} 缺失"
    assert "<script>function downloadTSV(){}</script>" in html  # JS 不转义注入
    assert "&lt;uncultured&gt;" in html                          # 元数据值被转义
    assert "✓" in html
    # 修复的断链 class: 专属 CSS 层必须真正定义 .toolbar/.section-title/.toast
    style = html.split("<style>", 1)[1].split("</style>", 1)[0]
    for cls in (".toolbar", ".btn", ".section-title", ".section-body",
                ".toast", ".legend", ".checklist", ".stat", ".progress"):
        assert cls in style, f"submission 专属 CSS 缺 {cls} 定义"
    assert not _EXT_REF_RE.search(html)
