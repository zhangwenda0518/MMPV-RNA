# -*- coding: utf-8 -*-
"""优化提升收尾测试 — triage 实战预设 / 报告互链与演化图 / provenance 导出。

全部离线。
"""

import json
from pathlib import Path

import pytest

from mmpv_common import run_digest
from mmpv_common.process_report import (collect_pipeline_reports,
                                        collect_process, collect_triage,
                                        generate_master_report)
from mmpv_common.provenance_export import export_provenance
from mmpv_common.triage import (build_ident_evidence, locate_ident_dir,
                                run_preset)

from test_ai_driven import make_project
from test_process_report import *  # noqa: F401,F403  (复用 project fixture)


# ═══════════════════════════════════════════════════════════════════
# triage 实战预设
# ═══════════════════════════════════════════════════════════════════

def _make_ident(project: Path) -> Path:
    """仿 02_Identification/<sample>/ 结构: result.id (每行一个 contig id)
    + all.candidate.fasta (候选全集, 含零命中 contig c4)。"""
    ident = project / "03_Discovery" / "02_Identification" / "S1"
    ident.mkdir(parents=True)
    (ident / "S1_virus.all.candidate.fasta").write_text(
        ">c1 desc\nACGT\n>c2 desc\nACGT\n>c3 desc\nACGT\n>c4 desc\nACGT\n",
        encoding="utf-8")
    (ident / "S1_virus.genomad.result.id").write_text("c1\nc2\nc3\n", encoding="utf-8")
    (ident / "S1_virus.virsorter2.result.id").write_text("c1\nc2\n", encoding="utf-8")
    (ident / "S1_virus.blast.result.id").write_text("c1\n", encoding="utf-8")
    return ident.parent


def test_locate_ident_dir_standard_and_missing(tmp_path):
    proj = tmp_path / "p"
    proj.mkdir()
    assert locate_ident_dir(proj) is None
    _make_ident(proj)
    assert locate_ident_dir(proj) == proj / "03_Discovery" / "02_Identification"


def test_build_ident_evidence_matrix(tmp_path):
    proj = tmp_path / "p"
    proj.mkdir()
    ident = _make_ident(proj)
    built = build_ident_evidence(ident)
    assert built["status"] == "ok" and built["n_contigs"] == 3  # c4 零命中不进
    assert built["n_unanimous"] == 1 and built["n_partial"] == 2
    lines = Path(built["evidence"]).read_text(encoding="utf-8").strip().splitlines()
    header = lines[0].split("\t")
    assert header[:2] == ["contig_id", "sample"]
    row = dict(zip(header, lines[1].split("\t")))   # c1 (排序首行)
    assert row["contig_id"] == "c1"
    assert row["genomad"] == "1" and row["virsorter2"] == "1" and row["blast"] == "1"


def test_run_preset_ident_rules(tmp_path):
    proj = tmp_path / "p"
    proj.mkdir()
    _make_ident(proj)
    report = run_preset("ident", proj, mode="rules")
    assert report["status"] == "ok" and report["preset"] == "ident"
    assert report["counts"] == {"KEEP": 1, "REVIEW": 2}   # c1 全票 / c2,c3 分歧
    out = Path(report["output"])
    assert out.name == "ident_evidence.triage.tsv"
    body = [l.split("\t") for l in out.read_text(encoding="utf-8").strip().splitlines()[1:]]
    by_id = {r[0]: r for r in body}
    assert by_id["c1"][-3] == "KEEP" and "3/3" in by_id["c1"][-2]
    assert by_id["c2"][-3] == "REVIEW" and "2/3" in by_id["c2"][-2]


def test_run_preset_ident_missing(tmp_path):
    res = run_preset("ident", tmp_path / "empty", mode="rules")
    assert res["status"] == "error" and "identification" in res["reason"]


def test_run_preset_eve_distinguish(tmp_path):
    proj = tmp_path / "p"
    proj.mkdir()
    tsv = proj / "06_EVE" / "eve_distinguish" / "dna_vs_eve_filter.tsv"
    tsv.parent.mkdir(parents=True)
    tsv.write_text("contig_id\tdecay\tdomain\taction\n"
                   "e1\thigh\tfull\tKEEP_virus\n"
                   "e2\tlow\tpartial\tREVIEW\n", encoding="utf-8")
    report = run_preset("eve_distinguish", proj, mode="rules")
    assert report["status"] == "ok"
    assert report["counts"]["REVIEW"] == 2   # decay/domain 两列值不一致 → REVIEW
    assert Path(report["output"]).exists()


def test_run_preset_unknown():
    assert run_preset("nope", Path(".")).get("status") == "error"


# ═══════════════════════════════════════════════════════════════════
# 用量聚合 / triage 扫描 / 报告互链 / 演化图
# ═══════════════════════════════════════════════════════════════════

def test_ai_by_month_and_op_tokens(project):
    rec = collect_process(project)
    ai = rec["ai"]
    assert ai["by_month"]["2026-10"]["calls"] == 2
    assert ai["by_month"]["2026-10"]["tokens"] == 200
    assert ai["by_op_tokens"]["ai_interpretation"] == 120


def test_collect_triage_and_panel_line(project):
    # 在发现管线根下放一份 triage 输出
    t = project / "03_Discovery" / "02_Identification" / "demo.triage.tsv"
    t.parent.mkdir(parents=True, exist_ok=True)
    t.write_text("contig_id\tgenomad\tai_action\tai_reason\tai_source\n"
                 "c1\tvirus\tKEEP\tok\trule\n"
                 "c2\tx\tREVIEW\tbad\trule\n"
                 "c3\ty\tREVIEW\tbad2\trule\n", encoding="utf-8")
    items = collect_triage(project)
    assert len(items) == 1
    assert items[0]["keep"] == 1 and items[0]["review"] == 2
    assert items[0]["file"].startswith("03_Discovery/")
    from mmpv_common.process_report import build_process_panel_html
    rec = collect_process(project)
    panel = build_process_panel_html(rec)
    assert "Triage annotations: <b>1</b> file(s) (KEEP 1 / REVIEW 2)" in panel


def test_collect_pipeline_reports_links(project):
    # 造两份已知落点的报告
    d = project / "03_Discovery" / "13_Report"
    d.mkdir(parents=True, exist_ok=True)
    (d / "pipeline_report.html").write_text("<html/>", encoding="utf-8")
    (project / "02_Preprocessing" / "preprocessing_report.html").write_text(
        "<html/>", encoding="utf-8")
    reps = collect_pipeline_reports(project)
    by_pk = {r["pipeline"]: r for r in reps}
    assert by_pk["discovery"]["exists"] == "yes"
    assert by_pk["discovery"]["rel"] == "03_Discovery/13_Report/pipeline_report.html"
    assert by_pk["preprocessing"]["exists"] == "yes"
    assert by_pk["eve"]["exists"] == "no"    # 不存在的如实标 no, 不造链接


def test_master_report_nav_svg_usage_triage(project):
    d = project / "03_Discovery" / "13_Report"
    d.mkdir(parents=True, exist_ok=True)
    (d / "pipeline_report.html").write_text("<html/>", encoding="utf-8")
    # triage 文件 + digest 历史(≥2) fixture 已具备
    t = project / "03_Discovery" / "02_Identification" / "demo.triage.tsv"
    t.parent.mkdir(parents=True, exist_ok=True)
    t.write_text("contig_id\tgenomad\tai_action\tai_reason\tai_source\n"
                 "c1\tvirus\tKEEP\tok\trule\n", encoding="utf-8")
    res = generate_master_report(project, want_ai=False)
    body = Path(res["path"]).read_text(encoding="utf-8")
    assert 'href="03_Discovery/13_Report/pipeline_report.html"' in body   # 导航 (相对 90_Handoff)
    assert "<svg" in body and "polyline" in body                          # 演化图
    assert "AI Usage" in body and "2026-10" in body                       # 按月用量
    assert "Triage Candidate Annotation Summary" in body and "demo.triage.tsv" in body    # triage 分区
    assert "KEEP</th>" in body


def test_master_report_panel_links_back(project):
    """主报告生成后, 再嵌面板的管线报告应反链主报告 (相对路径正确)。"""
    from mmpv_common.process_report import append_process_panel_to_html
    generate_master_report(project, want_ai=False)
    d = project / "03_Discovery" / "13_Report"
    d.mkdir(parents=True, exist_ok=True)
    html = d / "pipeline_report.html"
    html.write_text("<html><body>x</body></html>", encoding="utf-8")
    r = append_process_panel_to_html(html, project / "03_Discovery")
    assert r["status"] == "ok"
    body = html.read_text(encoding="utf-8")
    assert 'href="../../90_Handoff/MMPV_master_report.html"' in body


# ═══════════════════════════════════════════════════════════════════
# provenance 导出
# ═══════════════════════════════════════════════════════════════════

def test_export_provenance(project):
    res = export_provenance(project)
    assert res["status"] == "ok"
    assert Path(res["json"]).exists() and Path(res["md"]).exists()
    assert Path(res["md"]).parent.name == "90_Handoff"
    prov = json.loads(Path(res["json"]).read_text(encoding="utf-8"))
    assert prov["schema_version"] == 1
    assert prov["ai_summary"]["n_calls"] == 2
    assert len(prov["ai_log_events"]) == 3      # 2 call + 1 agent 事件全量在内
    assert "governance" in prov
    md = Path(res["md"]).read_text(encoding="utf-8")
    assert "方法学段落" in md and "ai_log.jsonl" in md
    assert "大语言模型" in md and "deterministic" in md   # 中英文段落齐备
    assert "| 操作 (op)" in md                             # 参与总览表


def test_export_provenance_no_ai(tmp_path):
    root = make_project(tmp_path / "bare")
    res = export_provenance(root)
    md = Path(res["md"]).read_text(encoding="utf-8")
    assert "纯规则运行" in md
    assert "未使用大语言模型" in md
