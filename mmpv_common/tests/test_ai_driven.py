# -*- coding: utf-8 -*-
"""四层架构集成单测 — L1 digest/卡片 · L3 顾问/triage · L4 copilot/MCP。

全部离线: run_digest 只读目录结构; advisor 只跑规则; triage 用 rules 模式
或桩 client; copilot/MCP 只校验不执行。
"""

import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from mmpv_common import mcp_server, run_digest, stage_cards
from mmpv_common.advisor import advise, advise_with_ai
from mmpv_common.copilot import suggest_intents, validate
from mmpv_common.triage import TRIAGE_COLS, triage


# ═══════════════════════════════════════════════════════════════════
# 造一个迷你 standard 布局项目
# ═══════════════════════════════════════════════════════════════════

def make_project(root, *, deep_diff=False, with_rescue=True):
    pub = root / "01_PublicData"
    (pub / "02_Unified").mkdir(parents=True)
    (pub / "02_Unified" / "Global_Unified_Metadata_Core14.csv").write_text(
        "Run,Location\nSRR1,China, Beijing, Beijing_AI\nSRR2,X\n", encoding="utf-8")
    raw = pub / "03_RawData"
    raw.mkdir()
    (raw / "s1.fastq.gz").write_bytes(b"x" * 100)

    pre = root / "02_Preprocessing"
    clean = pre / "01_CleanData"
    clean.mkdir(parents=True)
    if deep_diff:
        (clean / "a.fastq.gz").write_bytes(b"x" * 10_000_000)
        (clean / "b.fastq.gz").write_bytes(b"x" * 100_000)
    else:
        (clean / "a.fastq.gz").write_bytes(b"x" * 1000)
        (clean / "b.fastq.gz").write_bytes(b"x" * 900)
    (pre / "02_HostDepleted").mkdir()
    (pre / "preprocessing_summary.tsv").write_text("sample\treads\na\t1\nb\t2\n",
                                                   encoding="utf-8")

    disc = root / "03_Discovery"
    (disc / "01_Assembly" / "S1").mkdir(parents=True)
    (disc / "02_Identification" / "genomad").mkdir(parents=True)
    if with_rescue:
        (disc / "10_Rescue").mkdir(parents=True)

    ana = root / "04_Analysis"
    (ana / "01_Detection").mkdir(parents=True)
    (ana / "01_Detection" / "summary.tsv").write_text("virus\tcnt\nv1\t5\n", encoding="utf-8")
    (ana / "09_Report").mkdir()
    (ana / "09_Report" / "report.html").write_text("<html/>", encoding="utf-8")

    (root / "06_EVE" / "01_Loci" / "GenomeA").mkdir(parents=True)
    (root / "90_Handoff").mkdir()
    return root


@pytest.fixture()
def project(tmp_path):
    return make_project(tmp_path / "proj")


# ═══════════════════════════════════════════════════════════════════
# L1 stage_cards
# ═══════════════════════════════════════════════════════════════════

def test_stage_cards_self_check_clean():
    assert stage_cards.self_check() == []


def test_stage_cards_structure():
    assert "discovery" in stage_cards.pipelines()
    assert "rescue" in stage_cards.stages_of("discovery")
    card = stage_cards.get_card("discovery", "rescue")
    assert card["name"] and card["inputs"] and card["risks"] is not None
    assert stage_cards.get_card("nope") is None


# ═══════════════════════════════════════════════════════════════════
# L1 run_digest
# ═══════════════════════════════════════════════════════════════════

def test_digest_empty_project(tmp_path):
    d = run_digest.build_digest(tmp_path / "empty_proj")
    assert d["totals"]["pipelines_present"] == 0
    assert d["layout"] == "standard"


def test_digest_extracts_facts(project):
    d = run_digest.build_digest(project)
    pre = d["pipelines"]["preprocessing"]
    assert pre["present"] and pre["clean_reads"] == 2
    assert pre["summary_rows"] == 2
    disc = d["pipelines"]["discovery"]
    assert disc["assembly_samples"] == 1
    assert disc["stages"]["d_rescue"]["exists"]
    assert d["pipelines"]["eve"]["genomes"] == 1
    assert d["pipelines"]["analysis"]["detection_summary_rows"] == 1
    assert "report.html" in d["pipelines"]["analysis"]["report_html"][0] or \
           d["pipelines"]["analysis"]["report_html"]
    assert d["totals"]["pipelines_present"] >= 5


def test_digest_write_load_roundtrip(project):
    d = run_digest.build_digest(project)
    out = run_digest.write_digest(d, project)
    assert out.exists() and out.parent.name == "state" and out.parent.parent.name == "90_Handoff"
    d2 = run_digest.load_digest(out)
    assert d2["generated_at"] == d["generated_at"]
    bad = out.parent / "bad.json"
    bad.write_text("{}", encoding="utf-8")
    with pytest.raises(ValueError):
        run_digest.load_digest(bad)


def test_digest_compact_for_ai(project):
    d = run_digest.build_digest(project)
    slim = run_digest.compact_for_ai(d)
    pre = slim["pipelines"]["preprocessing"]
    assert pre["clean_reads"] == 2
    disc_json = json.dumps(slim["pipelines"]["discovery"]["stages"])
    assert "S1" not in disc_json and "genomad" not in disc_json  # 子目录名清单被剥掉


# ═══════════════════════════════════════════════════════════════════
# L3 advisor
# ═══════════════════════════════════════════════════════════════════

def _ids(recs):
    return [r.id for r in recs]


def test_advise_empty_project_warns(tmp_path):
    recs = advise(run_digest.build_digest(tmp_path / "nope"))
    assert "empty_project" in _ids(recs)
    assert recs[0].level == "warn"


def test_advise_raw_without_clean(project):
    # 造一个只有原始数据的项目
    import shutil
    src = project.parent / "raw_only"
    shutil.rmtree(src, ignore_errors=True)
    src.mkdir()
    raw = src / "01_PublicData" / "03_RawData"
    raw.mkdir(parents=True)
    (raw / "s.fastq.gz").write_bytes(b"x")
    recs = advise(run_digest.build_digest(src))
    assert "start_preprocessing" in _ids(recs)


def test_advise_clean_without_deplete(project):
    import shutil
    src = project.parent / "clean_only"
    shutil.rmtree(src, ignore_errors=True)
    shutil.copytree(project, src)
    shutil.rmtree(src / "02_Preprocessing" / "02_HostDepleted")
    recs = advise(run_digest.build_digest(src))
    assert "run_deplete" in _ids(recs)


def test_advise_depth_variance(project):
    import shutil
    src = project.parent / "deepdiff"
    shutil.rmtree(src, ignore_errors=True)
    shutil.copytree(project, src)
    cd = src / "02_Preprocessing" / "01_CleanData"
    # 4 个文件 (≥2 样本口径), 极差 100 倍 → 触发 bbnorm 建议
    (cd / "a_1.fastq.gz").write_bytes(b"x" * 10_000_000)
    (cd / "a_2.fastq.gz").write_bytes(b"x" * 10_000_000)
    (cd / "b_1.fastq.gz").write_bytes(b"x" * 100_000)
    (cd / "b_2.fastq.gz").write_bytes(b"x" * 100_000)
    recs = advise(run_digest.build_digest(src))
    assert "suggest_bbnorm" in _ids(recs)
    r = [x for x in recs if x.id == "suggest_bbnorm"][0]
    assert "abundance" in r.detail   # 提醒定量管线别用归一化 reads (detail 已英文化)


def test_advise_submission_ready(project):
    recs = advise(run_digest.build_digest(project))
    assert "submission_ready" in _ids(recs)


def test_advise_with_ai_offline_narrative_unavailable(project):
    payload = advise_with_ai(run_digest.build_digest(project),
                             config=SimpleNamespace(api_key=""))
    assert payload["recommendations"]
    assert payload["narrative"]["status"] in ("unavailable", "skipped")


def test_advise_with_ai_stub_narrative(project, tmp_path):
    resp = SimpleNamespace(
        choices=[SimpleNamespace(message=SimpleNamespace(content="先跑去宿主, 再进发现管线。"))],
        usage=None)

    class Stub:
        class chat:  # noqa: N801
            class completions:  # noqa: N801
                @staticmethod
                def create(**kw):
                    return resp

    payload = advise_with_ai(run_digest.build_digest(project),
                             state_dir=tmp_path, client=Stub())
    assert payload["narrative"]["status"] == "ok"
    assert "去宿主" in payload["narrative"]["text"]
    assert (tmp_path / "ai_log.jsonl").exists()   # 叙事也要溯源


# ═══════════════════════════════════════════════════════════════════
# L3 triage
# ═══════════════════════════════════════════════════════════════════

def _write_tsv(path, rows, header):
    lines = ["\t".join(header)] + ["\t".join(r) for r in rows]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def test_triage_rules_unanimous_and_conflict(tmp_path):
    tsv = tmp_path / "c.tsv"
    _write_tsv(tsv,
               [("c1", "virus", "virus"), ("c2", "virus", "plasmid"),
                ("c3", "", "")],
               ["contig_id", "genomad", "virsorter2"])
    out = tmp_path / "c.triage.tsv"
    report = triage(tsv, out, evidence_cols=["genomad", "virsorter2"],
                    id_col="contig_id", mode="rules")
    assert report["counts"]["KEEP"] == 1 and report["counts"]["REVIEW"] == 2
    lines = out.read_text(encoding="utf-8").strip().splitlines()
    assert lines[0].split("\t") == ["contig_id", "genomad", "virsorter2", *TRIAGE_COLS]
    body = [l.split("\t") for l in lines[1:]]
    assert body[0][3] == "KEEP" and "unanimous" in body[0][4]
    assert body[1][3] == "REVIEW" and "conflict" in body[1][4]
    assert body[2][3] == "REVIEW"   # 证据全空 → 保守


def test_triage_never_touches_original_columns(tmp_path):
    """宪法不变量: 原有列逐字节保留, 仅追加 ai_* 三列。"""
    tsv = tmp_path / "c.tsv"
    orig_rows = [("c1", "virus", "virus", "extra,值"), ("c2", "x", "y", "z")]
    _write_tsv(tsv, orig_rows, ["id", "e1", "e2", "note"])
    out = tmp_path / "out.tsv"
    triage(tsv, out, evidence_cols=["e1", "e2"], id_col="id", mode="rules")
    orig_lines = tsv.read_text(encoding="utf-8").strip().splitlines()
    out_lines = out.read_text(encoding="utf-8").strip().splitlines()
    for o, n in zip(orig_lines, out_lines):
        assert n.startswith(o + "\t")     # 原内容是输出的前缀 → 未被改动


def test_triage_rejects_existing_ai_columns(tmp_path):
    tsv = tmp_path / "c.tsv"
    _write_tsv(tsv, [("c1", "KEEP")], ["id", "ai_action"])
    with pytest.raises(ValueError):
        triage(tsv, tmp_path / "o.tsv", mode="rules")


def test_triage_ai_mode_schema_fallback(tmp_path):
    tsv = tmp_path / "c.tsv"
    _write_tsv(tsv, [("c1", "virus"), ("c2", "phage")], ["id", "tool"])
    # 桩: 输出只有一行且 action 非法 → 两行都应保守回退 REVIEW
    resp = SimpleNamespace(
        choices=[SimpleNamespace(message=SimpleNamespace(
            content='[{"id": "c1", "action": "REMOVE", "reason": "x"}]'))],
        usage=None)

    class Stub:
        class chat:  # noqa: N801
            class completions:  # noqa: N801
                @staticmethod
                def create(**kw):
                    return resp

    report = triage(tsv, tmp_path / "o.tsv", evidence_cols=["tool"],
                    id_col="id", mode="ai", client=Stub())
    assert report["counts"]["REVIEW"] == 2   # REMOVE 不在合法动作集 → 拒绝并回退


def test_triage_ai_mode_ok(tmp_path):
    tsv = tmp_path / "c.tsv"
    _write_tsv(tsv, [("c1", "virus"), ("c2", "phage")], ["id", "tool"])
    content = json.dumps([{"id": "c1", "action": "KEEP", "reason": "工具判病毒"},
                          {"id": "c2", "action": "REVIEW", "reason": "噬菌体混入植物样本"}],
                         ensure_ascii=False)
    resp = SimpleNamespace(
        choices=[SimpleNamespace(message=SimpleNamespace(content=content))],
        usage=SimpleNamespace(prompt_tokens=1, completion_tokens=1, total_tokens=2))

    class Stub:
        class chat:  # noqa: N801
            class completions:  # noqa: N801
                @staticmethod
                def create(**kw):
                    return resp

    report = triage(tsv, tmp_path / "o.tsv", evidence_cols=["tool"],
                    id_col="id", mode="ai", state_dir=tmp_path, client=Stub())
    assert report["counts"] == {"KEEP": 1, "REVIEW": 1}
    assert report["models"] and report["models"][0]
    assert (tmp_path / "ai_log.jsonl").exists()


# ═══════════════════════════════════════════════════════════════════
# L4 copilot
# ═══════════════════════════════════════════════════════════════════

GOOD_DISCOVERY = ("python virome_discovery_pipeline/virome_pipeline.py --stage rescue "
                  "--output_dir D:/proj/03_Discovery -t 40")


def test_validate_good_discovery_command():
    res = validate(GOOD_DISCOVERY)
    assert res.ok, res.errors
    assert "dry-run" in res.preview and "virome_pipeline.py" in res.preview
    assert not res.requires_confirmation


def test_validate_rejects_unknown_script():
    res = validate("python evil.py --stage all")
    assert not res.ok and "白名单" in res.errors[0]


def test_validate_rejects_shell_metachars():
    res = validate(GOOD_DISCOVERY + " && echo pwned")
    assert not res.ok and "元字符" in res.errors[0]


def test_validate_rejects_unknown_flag_and_stage():
    res = validate("python virome_discovery_pipeline/virome_pipeline.py "
                   "--stage not_a_stage --output_dir x --zzz 1")
    assert not res.ok
    assert any("not_a_stage" in e for e in res.errors)
    assert any("--zzz" in e for e in res.errors)


def test_validate_missing_required():
    res = validate("python virome_discovery_pipeline/virome_pipeline.py --stage rescue")
    assert not res.ok and any("--output_dir" in e for e in res.errors)


def test_validate_dangerous_force_needs_confirmation():
    res = validate("python endogenous_virus_pipeline/eve_screen.py -g g.fasta "
                   "-n sp -o out/ --force")
    assert res.ok and res.requires_confirmation
    assert any("--force" in w for w in res.warnings)


def test_validate_conditional_required():
    res = validate("python data_preprocessing_pipeline/data_preprocessing.py "
                   "--stage clean --output_dir out/")
    assert not res.ok and any("--input_reads" in e for e in res.errors)


def test_validate_analysis_needs_reference():
    res = validate("python virome_analysis_pipeline/auto_known_virus.py "
                   "--stage all --reads_dir rd/")
    assert not res.ok
    assert any("--reference" in e for e in res.errors)


def test_validate_public_inline_key_warns():
    # --stage 是 public_data_pipeline 的 argparse 必填 (spec 已按真实 CLI 收紧),
    # 本用例只锁"内联 key 告警"这一行为
    res = validate("python public_metadata_pipeline/public_data_pipeline.py "
                   "--species Tomato --taxid 4081 --stage search "
                   "--deepseek-api sk-secret123")
    assert res.ok and any("环境变量" in w for w in res.warnings)


def test_validate_submission_subcommands():
    """⑦ 是子命令编排器 (known/discovery), 无 --stage; 必填按子命令分开。"""
    base = "python virome_submission_pipeline/submission_pipeline.py "
    ok = validate(base + "discovery --work-dir D:/proj/03_Discovery/10_Rescue "
                         "--run-title demo --mode both")
    assert ok.ok, ok.errors
    assert "子命令: discovery" in ok.preview
    assert not validate(base + "--work-dir x --run-title t").ok        # 缺子命令
    assert not validate(base + "bogus --work-dir x").ok               # 未知子命令
    assert not validate(base + "discovery --fasta f.fa").ok           # 缺 --work-dir
    assert not validate(base + "discovery --work-dir x --stage all").ok   # 无 --stage
    known_ok = validate(base + "known --extraction-dir E --suvtk-db D -o out/")
    assert known_ok.ok, known_ok.errors                              # -o 是 known 的 --output
    assert not validate(base + "known --suvtk-db D").ok              # known 缺序列来源


def test_validate_hostref_build_alias():
    """P2-28: build_host_pipeline 由 "hostref"/"hostref_build" 两个键共同匹配。"""
    from mmpv_common.copilot import COMMAND_SPECS
    assert COMMAND_SPECS["hostref_build"] is COMMAND_SPECS["hostref"]
    cmd = ("python public_metadata_pipeline/build_host_pipeline.py "
           "--species S --taxid 1 --stage all --threads 8")
    assert validate(cmd).ok
    assert not validate("python public_metadata_pipeline/build_host_pipeline.py "
                        "--species S --taxid 1").ok      # --stage 必填


def test_interpreter_whitelist_strict():
    """P1-21: pythonXxx/pythonevil 不是解释器; python3.11 / python.exe 是。"""
    base = "virome_discovery_pipeline/virome_pipeline.py --stage rescue --output_dir o"
    assert validate(f"python3.11 {base}").ok
    assert validate(f"python.exe {base}").ok
    assert validate(f"C:/Python312/python.exe {base}").ok
    assert not validate(f"pythonXxx {base}").ok
    assert not validate(f"pythonevil {base}").ok


def test_validate_errors_are_rejections_not_confirmations():
    """P2-29b: 校验失败就是拒绝 —— 不能退化成"确认后可跑"。"""
    res = validate("python virome_discovery_pipeline/virome_pipeline.py --stage rescue")
    assert not res.ok and not res.requires_confirmation


def test_plan_with_ai_logs_structure_not_free_text(tmp_path):
    """P2-29d: ai_log 只记管线名 + 子命令 + flag 键名, 不记用户原话/AI 原文。"""
    from mmpv_common.copilot import plan_with_ai
    secret = "我的样本在 D:/保密/raw 目录, 请帮我跑 ⑦"
    answer = json.dumps({"pipeline": "submission", "subcommand": "discovery",
                         "stages": [], "args": {"--work-dir": "D:/proj/10_Rescue",
                                                "--run-title": "demo", "--bogus": "x"},
                         "missing": [], "explain": "submit"},
                        ensure_ascii=False)
    resp = SimpleNamespace(choices=[SimpleNamespace(
        message=SimpleNamespace(content=answer))], usage=None)

    class Stub:
        class chat:  # noqa: N801
            class completions:  # noqa: N801
                @staticmethod
                def create(**kw):
                    return resp

    out = plan_with_ai(secret, state_dir=tmp_path, client=Stub())
    assert out["status"] == "ok" and out["validation"]["ok"]
    assert "discovery" in out["command"] and "--work-dir" in out["command"]
    assert "--bogus" not in out["command"]          # 未登记 flag 不拼进草稿
    log = (tmp_path / "ai_log.jsonl").read_text(encoding="utf-8")
    assert "plan_with_ai" in log and '"pipeline": "submission"' in log
    assert "保密" not in log and "请帮我跑" not in log       # 用户自由文本不落盘
    assert '"--work-dir"' in log and "D:/proj/10_Rescue" not in log   # 只记键名不记值


# ═══════════════════════════════════════════════════════════════════
# spec × 真实 argparse 一致性 (spec 漂移回归: 只查"幻影 flag")
# ═══════════════════════════════════════════════════════════════════

_REPO = Path(__file__).resolve().parents[2]


def _real_flags(path: Path) -> set:
    """AST 提取 add_argument 的 flag 面 (BooleanOptionalAction 补隐式 --no-XXX)。"""
    import ast
    tree = ast.parse(path.read_text(encoding="utf-8"))
    out = set()
    for n in ast.walk(tree):
        if isinstance(n, ast.Call) and isinstance(n.func, ast.Attribute) \
                and n.func.attr == "add_argument":
            names = [a.value for a in n.args
                     if isinstance(a, ast.Constant) and isinstance(a.value, str)]
            out.update(names)
            kw = {k.arg: k.value for k in n.keywords}
            if kw.get("action") is not None and "BooleanOptionalAction" in ast.dump(kw["action"]):
                base = next((x for x in names if x.startswith("--")), None)
                if base:
                    out.add("--no-" + base[2:])
    return {f for f in out if f not in ("-h", "--help")}


@pytest.mark.parametrize("key", ["preprocessing", "discovery", "analysis", "public",
                                 "hostref", "phylo", "submission"])
def test_spec_flags_exist_in_real_argparse(key):
    """spec 不许登记脚本里不存在的 flag (P1-19/P2-28 的漂移回归)。

    反向 (脚本新加的 flag 未登记) 不断言 —— 那只让草稿被拒, 不会产生幻影放行。
    """
    from mmpv_common.copilot import COMMAND_SPECS
    spec = COMMAND_SPECS[key]
    path = _REPO / spec.script
    if not path.exists():
        pytest.skip(f"编排器不在仓库: {spec.script} (归档/移动)")
    real = _real_flags(path)
    declared = set(spec.required) | set(spec.optional) | set(spec.dangerous)
    for sub in spec.subcommands:
        declared |= set(sub.required) | set(sub.optional)
    phantom = sorted(f for f in declared if f.startswith("--") and f not in real)
    assert not phantom, f"{key}: spec 登记了脚本没有的 flag: {phantom}"


def test_suggest_intents_offline():
    hits = suggest_intents("帮我做病毒发现和组装")
    assert hits and hits[0]["pipeline"] == "discovery"
    assert suggest_intents("天气预报") == []


# ═══════════════════════════════════════════════════════════════════
# L4 MCP server (进程内, 不起子进程)
# ═══════════════════════════════════════════════════════════════════

def _rpc(msg):
    return mcp_server.handle(msg)


def test_mcp_initialize():
    out = _rpc({"jsonrpc": "2.0", "id": 1, "method": "initialize",
                "params": {"protocolVersion": "2024-11-05"}})
    assert out["result"]["serverInfo"]["name"] == "mmpv"
    assert "tools" in out["result"]["capabilities"]


def test_mcp_notification_no_response():
    assert _rpc({"jsonrpc": "2.0", "method": "notifications/initialized"}) is None


def test_mcp_tools_list():
    out = _rpc({"jsonrpc": "2.0", "id": 2, "method": "tools/list"})
    names = {t["name"] for t in out["result"]["tools"]}
    assert names == {"mmpv_project_digest", "mmpv_advise", "mmpv_stage_cards",
                     "mmpv_validate_command", "mmpv_stage_tools",
                     "mmpv_plan_chain", "mmpv_ask"}


def test_mcp_plan_chain_readonly(project):
    out = _rpc({"jsonrpc": "2.0", "id": 7, "method": "tools/call",
                "params": {"name": "mmpv_plan_chain",
                           "arguments": {"project_dir": str(project)}}})
    body = json.loads(out["result"]["content"][0]["text"])
    assert body["n_steps"] >= 1
    assert all("command" in s for s in body["steps"])
    # 不存在任何执行类工具 — 红线
    names = {t["name"] for t in _rpc({"jsonrpc": "2.0", "id": 8,
             "method": "tools/list"})["result"]["tools"]}
    assert not any("exec" in n or "run_pipeline" in n or n.endswith("_execute")
                   for n in names)


def test_mcp_call_validate_tool():
    out = _rpc({"jsonrpc": "2.0", "id": 3, "method": "tools/call",
                "params": {"name": "mmpv_validate_command",
                           "arguments": {"command": GOOD_DISCOVERY}}})
    body = json.loads(out["result"]["content"][0]["text"])
    assert body["ok"] is True and "dry-run" in body["preview"]


def test_mcp_call_unknown_tool_is_error():
    out = _rpc({"jsonrpc": "2.0", "id": 4, "method": "tools/call",
                "params": {"name": "run_pipeline", "arguments": {}}})
    assert out["error"]["code"] == -32602


def test_mcp_digest_and_cards(project):
    out = _rpc({"jsonrpc": "2.0", "id": 5, "method": "tools/call",
                "params": {"name": "mmpv_project_digest",
                           "arguments": {"project_dir": str(project)}}})
    body = json.loads(out["result"]["content"][0]["text"])
    assert body["pipelines"]["preprocessing"]["present"]
    out2 = _rpc({"jsonrpc": "2.0", "id": 6, "method": "tools/call",
                 "params": {"name": "mmpv_stage_cards",
                            "arguments": {"pipeline": "eve", "stage": "verdict"}}})
    card = json.loads(out2["result"]["content"][0]["text"])
    assert card["name"]


def _digest_rpc(project, **arguments):
    out = _rpc({"jsonrpc": "2.0", "id": 9, "method": "tools/call",
                "params": {"name": "mmpv_project_digest",
                           "arguments": {"project_dir": str(project), **arguments}}})
    return json.loads(out["result"]["content"][0]["text"])


def test_mcp_digest_write_is_opt_in_and_single_build(project):
    """P2-29a: write 默认关且仅显式 true 生效; 落盘内容 = 返回的那一份 digest。"""
    from mmpv_common import run_digest
    digest_file = run_digest.digest_path(project)
    assert "written_to" not in _digest_rpc(project)                      # 默认不落盘
    assert not digest_file.exists()
    assert "written_to" not in _digest_rpc(project, write="true")        # 字符串不算显式 true
    assert not digest_file.exists()
    body = _digest_rpc(project, write=True, compact=False)
    snap = json.loads(digest_file.read_text(encoding="utf-8"))
    assert body["written_to"] == str(digest_file)
    assert snap == {k: v for k, v in body.items() if k != "written_to"}  # 同一份, 非二次构建
