"""
Tests for shared pipeline utilities: Checkpoint, UI, run_cmd, sanitizers.

Run:  cd public_metadata_pipeline && python -m pytest tests/ -v
"""

import os
import sys
import json
import tempfile
import shutil
from pathlib import Path

import pytest

# Ensure the package root is on sys.path so imports work
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from utils.pipeline_utils import UI, Checkpoint, run_cmd


# ── Checkpoint ──────────────────────────────────────────────────────────────

class TestCheckpoint:
    """JSON checkpoint state-machine tests."""

    @pytest.fixture
    def tmp_work_dir(self):
        d = tempfile.mkdtemp(prefix="ckpt_test_")
        yield d
        shutil.rmtree(d, ignore_errors=True)

    def test_initial_state(self, tmp_work_dir):
        ckpt = Checkpoint(tmp_work_dir)
        assert ckpt.is_done("search") is False
        assert ckpt.state == {"stages": {}}

    def test_mark_start_and_is_done(self, tmp_work_dir):
        ckpt = Checkpoint(tmp_work_dir)
        ckpt.mark_start("search")
        assert ckpt.is_done("search") is False
        assert ckpt.state["stages"]["search"]["status"] == "running"

    def test_mark_done(self, tmp_work_dir):
        ckpt = Checkpoint(tmp_work_dir)
        ckpt.mark_start("search")
        ckpt.mark_done("search")
        assert ckpt.is_done("search") is True
        assert "completed" in ckpt.state["stages"]["search"]

    def test_mark_fail(self, tmp_work_dir):
        ckpt = Checkpoint(tmp_work_dir)
        ckpt.mark_start("search")
        ckpt.mark_fail("search", "network timeout")
        assert ckpt.is_done("search") is False
        assert ckpt.state["stages"]["search"]["status"] == "failed"
        assert "network timeout" in ckpt.state["stages"]["search"]["error"]

    def test_reset(self, tmp_work_dir):
        ckpt = Checkpoint(tmp_work_dir)
        ckpt.mark_done("search")
        ckpt.mark_done("info")
        ckpt.reset()
        assert ckpt.state == {"stages": {}}
        assert ckpt.is_done("search") is False

    def test_persistence(self, tmp_work_dir):
        """State survives a second Checkpoint instantiation (real disk write)."""
        ckpt1 = Checkpoint(tmp_work_dir)
        ckpt1.mark_done("search")
        ckpt2 = Checkpoint(tmp_work_dir)
        assert ckpt2.is_done("search") is True

    def test_summary(self, tmp_work_dir):
        ckpt = Checkpoint(tmp_work_dir)
        ckpt.mark_done("search")
        ckpt.mark_start("info")
        summary = ckpt.summary(["search", "info", "down"])
        assert "[OK]" in summary
        assert "[..]" in summary
        assert "[  ]" in summary  # down — untouched = pending

    def test_fail_error_truncation(self, tmp_work_dir):
        ckpt = Checkpoint(tmp_work_dir)
        long_err = "x" * 300
        ckpt.mark_start("down")
        ckpt.mark_fail("down", long_err)
        stored = ckpt.state["stages"]["down"]["error"]
        assert len(stored) <= 200

    def test_unknown_stage_is_not_done(self, tmp_work_dir):
        ckpt = Checkpoint(tmp_work_dir)
        assert ckpt.is_done("nonexistent") is False


# ── UI ──────────────────────────────────────────────────────────────────────

class TestUI:
    def test_color_map_keys(self):
        expected = {"cyan", "green", "yellow", "red", "purple", "gray", "bold", "reset"}
        assert set(UI.C.keys()) == expected

    def test_color_values_are_ansi(self):
        for c in UI.C.values():
            assert c.startswith("\033[") or c == ""

    def test_stage_status_sym(self, capsys):
        UI.stage("test stage", "skip")
        captured = capsys.readouterr().out
        assert "[---]" in captured
        assert "test stage" in captured


# ── run_cmd ─────────────────────────────────────────────────────────────────

class TestRunCmd:

    def test_successful_command(self, tmp_path):
        log_dir = str(tmp_path / "logs")
        rc = run_cmd("echo hello", "echo_test", log_dir, timeout=10)
        assert rc == 0
        log_file = os.path.join(log_dir, "echo_test.log")
        assert os.path.isfile(log_file)
        with open(log_file, "r", encoding="utf-8") as f:
            content = f.read()
        assert "hello" in content

    def test_failed_command(self, tmp_path):
        log_dir = str(tmp_path / "logs")
        # On Windows 'false' may not exist; use a guaranteed failing command
        cmd = 'python -c "import sys; sys.exit(42)"'
        rc = run_cmd(cmd, "fail_test", log_dir, timeout=10)
        assert rc == 42

    def test_secret_masking(self, tmp_path):
        log_dir = str(tmp_path / "logs")
        secret = "sk-secret-key-123"
        cmd = f"echo {secret}"
        rc = run_cmd(cmd, "mask_test", log_dir, timeout=10, secrets=[secret])
        assert rc == 0
        log_file = os.path.join(log_dir, "mask_test.log")
        with open(log_file, "r", encoding="utf-8") as f:
            content = f.read()
        assert secret not in content
        assert "***" in content

    def test_log_dir_auto_created(self, tmp_path):
        log_dir = str(tmp_path / "nonexistent" / "logs")
        run_cmd("echo ok", "create_test", log_dir, timeout=10)
        assert os.path.isdir(log_dir)


# ── sanitize_dirname (from gsa_sra.down.py) ─────────────────────────────────

# Import the function directly
import importlib.util
_down_spec = importlib.util.spec_from_file_location(
    "gsa_sra_down",
    Path(__file__).resolve().parent.parent / "gsa_sra.down.py",
)
_down = importlib.util.module_from_spec(_down_spec)
_down_spec.loader.exec_module(_down)


class TestSanitizeDirname:
    def test_empty(self):
        assert _down.sanitize_dirname("") == "Uncategorized"

    def test_normal(self):
        assert _down.sanitize_dirname("leaf_sample") == "leaf_sample"

    def test_special_chars(self):
        result = _down.sanitize_dirname("leaf: sample<1>?")
        assert ":" not in result
        assert "<" not in result
        assert ">" not in result
        assert "?" not in result

    def test_whitespace_replaced(self):
        result = _down.sanitize_dirname("leaf sample 1")
        assert " " not in result
        assert "_" in result


# ── clean_series (from gsa_sra.plot.py) ────────────────────────────────────

import pandas as pd

_plot_spec = importlib.util.spec_from_file_location(
    "gsa_sra_plot",
    Path(__file__).resolve().parent.parent / "gsa_sra.plot.py",
)
_plot = importlib.util.module_from_spec(_plot_spec)
_plot_spec.loader.exec_module(_plot)


class TestCleanSeries:
    def test_removes_nan(self):
        s = pd.Series(["leaf", None, "root", pd.NA])
        result = _plot.clean_series(s)
        assert list(result) == ["leaf", "root"]

    def test_removes_not_provided(self):
        s = pd.Series(["leaf", "not_provided", "Not_Provided", "root"])
        result = _plot.clean_series(s)
        assert list(result) == ["leaf", "root"]

    def test_removes_missing(self):
        s = pd.Series(["leaf", "missing", "not collected"])
        result = _plot.clean_series(s)
        assert list(result) == ["leaf"]

    def test_removes_empty_string(self):
        s = pd.Series(["leaf", "", "root"])
        result = _plot.clean_series(s)
        assert list(result) == ["leaf", "root"]

    def test_all_invalid_returns_empty(self):
        s = pd.Series(["not_provided", "nan", "missing"])
        result = _plot.clean_series(s)
        assert len(result) == 0

    def test_case_insensitive(self):
        s = pd.Series(["leaf", "NOT_PROVIDED", "UNKNOWN"])
        result = _plot.clean_series(s)
        assert list(result) == ["leaf"]


# ── sanitize_for_json (from gsa_sra.info.py) ────────────────────────────────

_info_spec = importlib.util.spec_from_file_location(
    "gsa_sra_info",
    Path(__file__).resolve().parent.parent / "gsa_sra.info.py",
)
_info = importlib.util.module_from_spec(_info_spec)
_info_spec.loader.exec_module(_info)


class TestSanitizeForJSON:
    def test_plain_dict(self):
        # Create a minimal GSAPipeline to access sanitize_for_json
        # GSAPipeline.__init__ requires (mode, api_client, ai_model, out_dir, fill_date)
        gsa = _info.GSAPipeline(
            mode="local",
            api_client=None,
            ai_model=None,
            out_dir=tempfile.mkdtemp(prefix="gsa_test_"),
            fill_date=False,
        )
        result = gsa.sanitize_for_json({"a": 1, "b": "hello"})
        assert result == {"a": 1, "b": "hello"}

    def test_handles_pd_na(self):
        gsa = _info.GSAPipeline(
            mode="local", api_client=None, ai_model=None,
            out_dir=tempfile.mkdtemp(prefix="gsa_test_"), fill_date=False,
        )
        result = gsa.sanitize_for_json({"a": 1, "b": pd.NA})
        assert result == {"a": 1, "b": None}

    def test_handles_none(self):
        gsa = _info.GSAPipeline(
            mode="local", api_client=None, ai_model=None,
            out_dir=tempfile.mkdtemp(prefix="gsa_test_"), fill_date=False,
        )
        result = gsa.sanitize_for_json({"x": None})
        assert result == {"x": None}

    def test_handles_nested_dict(self):
        gsa = _info.GSAPipeline(
            mode="local", api_client=None, ai_model=None,
            out_dir=tempfile.mkdtemp(prefix="gsa_test_"), fill_date=False,
        )
        result = gsa.sanitize_for_json({"outer": {"inner": pd.NA}})
        assert result == {"outer": {"inner": None}}

    def test_json_serializable(self):
        import json
        gsa = _info.GSAPipeline(
            mode="local", api_client=None, ai_model=None,
            out_dir=tempfile.mkdtemp(prefix="gsa_test_"), fill_date=False,
        )
        data = {"name": "test", "tissue": pd.NA, "count": 5}
        result = gsa.sanitize_for_json(data)
        json.dumps(result)  # must not raise


# ── metadata backfill (from gsa_sra.search.py, DeepSeek Search API) ─────────

_search_spec = importlib.util.spec_from_file_location(
    "gsa_sra_search",
    Path(__file__).resolve().parent.parent / "gsa_sra.search.py",
)
_search = importlib.util.module_from_spec(_search_spec)
_search_spec.loader.exec_module(_search)


class TestBuildSearchQuery:
    def test_contains_accession_species_and_keyword(self):
        q = _search.build_search_query("Lycium barbarum", "SAMC1234567", "PRJCA002345",
                                       ["CollectionDate"])
        assert "SAMC1234567" in q
        assert "Lycium barbarum" in q
        assert "collection date" in q

    def test_truncated_to_128_chars(self):
        q = _search.build_search_query("x" * 200, "SAMC1", "PRJCA1",
                                       ["CollectionDate", "Location", "Tissue"])
        assert len(q) <= 128

    def test_without_accessions(self):
        q = _search.build_search_query("Arabidopsis thaliana", "", "", ["Tissue"])
        assert "Arabidopsis thaliana" in q
        assert "tissue" in q


class TestValidBackfillValue:
    def test_accepts_location(self):
        v = "China, Qinghai, Qaidam Basin"
        assert _search._valid_backfill_value("Location", v) == v

    def test_accepts_date_with_year(self):
        assert _search._valid_backfill_value("CollectionDate", "collected May 2021") == "collected May 2021"

    def test_rejects_date_without_year(self):
        assert _search._valid_backfill_value("CollectionDate", "collected recently") is None

    def test_rejects_junk(self):
        for junk in ("unknown", "N/A", "not collected", "missing", "", "x"):
            assert _search._valid_backfill_value("Location", junk) is None

    def test_rejects_overlong(self):
        assert _search._valid_backfill_value("Location", "y" * 300) is None

    def test_strips_quotes(self):
        assert _search._valid_backfill_value("Tissue", '"leaf"') == "leaf"


class TestIsMissingMeta:
    def test_na_and_junk(self):
        assert _search._is_missing_meta(pd.NA) is True
        assert _search._is_missing_meta(None) is True
        assert _search._is_missing_meta("Not_Provided") is True
        assert _search._is_missing_meta("not collected") is True

    def test_real_values(self):
        assert _search._is_missing_meta("leaf") is False
        assert _search._is_missing_meta("China, Ningxia, Yinchuan_AI") is False


class _FakeMessage:
    content = ('{"collection_date": "2021-05", "location": "China, Qinghai, Qaidam Basin", '
               '"tissue": "root", "growth_stage": "5-year-old plant", '
               '"evidence_url": "https://example.org/paper"}')


class _FakeChoice:
    message = _FakeMessage()


class _FakeResponse:
    choices = [_FakeChoice()]


class _FakeCompletions:
    def create(self, **kwargs):
        return _FakeResponse()


class _FakeChat:
    completions = _FakeCompletions()


class _FakeClient:
    chat = _FakeChat()


def _fake_search_results(calls):
    def _fake(query, api_key, api_base, count=8, ori_question=None, need_content=False):
        calls.append(query)
        return [{
            "title": "Study on wolfberry root samples",
            "url": "https://example.org/paper",
            "snippet": "roots collected in May 2021 in the Qaidam Basin, Qinghai, China",
            "summary": "Roots of 5-year-old Lycium barbarum plants collected in May 2021 "
                       "from the Qaidam Basin, Qinghai, China.",
            "rerank_score": 2.0,
            "date": None,
        }]
    return _fake


class TestRunAiMetadataBackfill:
    def _make_df(self, tissues):
        rows = []
        for i, tissue in enumerate(tissues):
            rows.append({"Run": f"CRR{i+1}", "Database": "GSA", "BioProject": "PRJCA1",
                         "BioSample": "SAMC1", "ScientificName": "Lycium barbarum",
                         "Location": pd.NA, "Tissue": tissue, "Age_GrowthStage": pd.NA,
                         "CollectionDate": pd.NA})
        return pd.DataFrame(rows)

    def test_backfills_missing_fields_and_logs(self, tmp_path, monkeypatch):
        calls = []
        monkeypatch.setattr(_search, "_ds_search", _fake_search_results(calls))
        df = self._make_df([pd.NA, "leaf"])

        out = _search.run_ai_metadata_backfill(df, "sk-test", "https://api.deepseek.com",
                                               "deepseek-chat", str(tmp_path), client=_FakeClient())
        # 同一样本的两个 Run 共享一次检索
        assert len(calls) == 1
        assert "SAMC1" in calls[0]

        row1 = out[out["Run"] == "CRR1"].iloc[0]
        row2 = out[out["Run"] == "CRR2"].iloc[0]
        assert row1["CollectionDate"] == "2021-05"
        # Location 统一成「国家, 省/州, 市/县_AI」一种写法 (与洗髓/建库一致)
        assert row1["Location"] == "China, Qinghai, Qaidam Basin_AI"
        assert row1["Tissue"] == "root"
        assert row1["Age_GrowthStage"] == "5-year-old plant"
        assert "CollectionDate" in str(row1["AI_Backfill"])
        # 已有值不被覆盖
        assert row2["Tissue"] == "leaf"
        assert "Tissue" not in str(row2["AI_Backfill"])

        log_file = tmp_path / "metadata_backfill_log.json"
        assert log_file.is_file()
        log = json.loads(log_file.read_text(encoding="utf-8"))
        assert log[0]["filled"]
        assert log[0]["urls"]

    def test_no_missing_skips_search(self, tmp_path, monkeypatch):
        calls = []
        monkeypatch.setattr(_search, "_ds_search", _fake_search_results(calls))
        df = self._make_df(["leaf"])
        df.loc[0, ["Location", "Age_GrowthStage", "CollectionDate"]] = \
            ["China, Qinghai, X", "5-year-old plant", "2021-05"]

        out = _search.run_ai_metadata_backfill(df, "sk-test", "https://api.deepseek.com",
                                               "deepseek-chat", str(tmp_path), client=_FakeClient())
        assert len(calls) == 0
        assert out["Location"].iloc[0] == "China, Qinghai, X"

    def test_extract_failure_leaves_df_untouched(self, tmp_path, monkeypatch):
        monkeypatch.setattr(_search, "_ds_search",
                            lambda *a, **k: [{"title": "t", "url": "u", "summary": "s",
                                              "rerank_score": 2.0, "date": None}])

        class _BadCompletions:
            def create(self, **kwargs):
                raise RuntimeError("boom")

        class _BadChat:
            completions = _BadCompletions()

        class _BadClient:
            chat = _BadChat()

        df = self._make_df([pd.NA])
        out = _search.run_ai_metadata_backfill(df, "sk-test", "https://api.deepseek.com",
                                               "deepseek-chat", str(tmp_path), client=_BadClient())
        assert out["Location"].isna().all()
        assert out["CollectionDate"].isna().all()


# ── search → info 检索回填贯通 (merge_search_backfill, 修复 Core14 丢回填) ────

class TestMergeSearchBackfill:
    """search 阶段 --deepseek-search 的回填值必须能进 Core14, 且只补空值。"""

    COLS = ["Run", "CollectionDate", "Location", "Tissue", "Age_GrowthStage", "AI_Backfill"]

    def _write_backfill_csv(self, tmp_path, rows):
        import csv as _csv
        p = tmp_path / "SRA_GSA_Merged_Final.csv"
        with open(p, "w", encoding="utf-8-sig", newline="") as f:
            w = _csv.DictWriter(f, fieldnames=self.COLS)
            w.writeheader()
            for r in rows:
                w.writerow(r)
        return str(p)

    def _make_core(self):
        return pd.DataFrame([
            {"Run": "CRR1", "ReleaseDate": "2023", "CollectionDate": pd.NA,
             "Location": pd.NA, "Tissue": pd.NA, "Age_GrowthStage": pd.NA,
             "ScientificName": "Lycium barbarum"},
            {"Run": "CRR2", "ReleaseDate": "2023", "CollectionDate": "2021-05",
             "Location": "China, Ningxia, Yinchuan_AI", "Tissue": "leaf",
             "Age_GrowthStage": pd.NA, "ScientificName": "Lycium barbarum"},
        ])

    def test_fills_only_blanks_and_records_provenance(self, tmp_path):
        csv_p = self._write_backfill_csv(tmp_path, [
            {"Run": "CRR1", "CollectionDate": "2021-05",
             "Location": "China, Qinghai, Qaidam Basin", "Tissue": "Roots",
             "Age_GrowthStage": "5-year-old plant",
             "AI_Backfill": "CollectionDate|Location|Tissue|Age_GrowthStage"},
            {"Run": "CRR2", "CollectionDate": "2018-01",
             "Location": "China, Gansu, Lanzhou", "Tissue": "fruit",
             "Age_GrowthStage": "3-year-old plant",
             "AI_Backfill": "CollectionDate|Location|Tissue|Age_GrowthStage"},
        ])
        out, stats = _info.merge_search_backfill(self._make_core(), csv_p, str(tmp_path))

        r1 = out[out["Run"] == "CRR1"].iloc[0]
        assert r1["CollectionDate"] == "2021-05"
        assert r1["Location"] == "China, Qinghai, Qaidam Basin_AI"  # 统一 _AI 来源标记
        assert r1["Tissue"] == "root"  # 部位归一化 (Roots -> root)
        assert r1["Age_GrowthStage"] == "5-year-old plant"
        assert set(r1["AI_Backfill_Search"].split("|")) == {
            "CollectionDate", "Location", "Tissue", "Age_GrowthStage"}

        # 已有真值不被检索回填覆盖
        r2 = out[out["Run"] == "CRR2"].iloc[0]
        assert r2["CollectionDate"] == "2021-05"
        assert r2["Tissue"] == "leaf"
        assert r2["Age_GrowthStage"] == "3-year-old plant"
        assert r2["AI_Backfill_Search"] == "Age_GrowthStage"

        assert stats == {"filled_runs": 2, "filled_cells": 5}
        assert (tmp_path / "search_backfill_merge.json").is_file()

    def test_ignores_invalid_date_and_untagged_fields(self, tmp_path):
        csv_p = self._write_backfill_csv(tmp_path, [
            {"Run": "CRR1", "CollectionDate": "recently", "Location": "China, X",
             "Tissue": pd.NA, "Age_GrowthStage": pd.NA,
             "AI_Backfill": "CollectionDate"},  # Location 未打标 → 不采纳
        ])
        out, stats = _info.merge_search_backfill(self._make_core(), csv_p)
        r1 = out[out["Run"] == "CRR1"].iloc[0]
        assert pd.isna(r1["CollectionDate"])  # 无 4 位年份 → 丢弃
        assert pd.isna(r1["Location"])         # 未声明证据 → 不采纳
        assert stats["filled_runs"] == 0

    def test_missing_source_file_is_noop(self, tmp_path):
        df = self._make_core()
        out, stats = _info.merge_search_backfill(df, str(tmp_path / "nope.csv"))
        assert stats == {"filled_runs": 0, "filled_cells": 0}
        assert out.equals(df)

# ── AI 值裁决 (ai_value_policy): 非破坏性清洗 + Location 统一写法 ──────────

import ai_value_policy as _policy
import ai_trace as _trace


class TestAiValuePolicy:
    def test_blank_caliber(self):
        for v in (None, pd.NA, "", "  ", "nan", "Not_Provided", "not collected", "未知"):
            assert _policy.is_blank(v) is True
        for v in ("leaf", "China, Ningxia, Yinchuan_AI", "2021-05"):
            assert _policy.is_blank(v) is False

    def test_format_location_single_suffix_idempotent(self):
        assert _policy.format_location("China:Yinchuan") == "China, Yinchuan_AI"
        assert _policy.format_location("China:Ningxia:Yinchuan_AI") == "China, Ningxia, Yinchuan_AI"
        assert _policy.format_location("China, Ningxia, Yinchuan_AI") == "China, Ningxia, Yinchuan_AI"

    def test_ai_value_never_deletes_existing_truth(self):
        df = pd.DataFrame([{"Location": pd.NA, "Tissue": "leaf", "Age_GrowthStage": pd.NA,
                            "CollectionDate": "2021-05"}])
        adopted, cf = _policy.apply_updates(df, 0, {
            "Location": "China, Ningxia, Yinchuan",
            "Tissue": "not_provided",          # AI 想清空 -> 必须忽略
            "Age_GrowthStage": "5 years, young fruit stage",
            "CollectionDate": "2022-03",
        })
        assert df.at[0, "Location"] == "China, Ningxia, Yinchuan_AI"
        assert df.at[0, "Tissue"] == "leaf"                     # 保留原值 (AI 只给占位词 = 无建议)
        assert df.at[0, "Age_GrowthStage"] == "5 years, young fruit stage"
        assert df.at[0, "CollectionDate"] == "2021-05"          # 冲突不覆盖
        assert "Tissue" not in str(df.at[0, "AI_Conflict"])     # 占位词不是冲突
        assert "CollectionDate" in str(df.at[0, "AI_Conflict"])
        assert set(adopted) == {"Location", "Age_GrowthStage"}
        assert {c["field"] for c in cf} == {"CollectionDate"}

    def test_enrichment_adopted_then_kept(self):
        # 旧值只是新值的子集 (AI 补全省份) -> 采纳更完整的写法
        df = pd.DataFrame([{"Location": "China, Yinchuan_AI"}])
        _policy.apply_updates(df, 0, {"Location": "China, Ningxia, Yinchuan"})
        assert df.at[0, "Location"] == "China, Ningxia, Yinchuan_AI"
        # 再次应用同一结论 -> 幂等, 且不新增冲突
        adopted, cf = _policy.apply_updates(df, 0, {"Location": "China, Ningxia, Yinchuan"})
        assert adopted == []
        assert cf == []


class _SanitizeMessage:
    def __init__(self, payload):
        self.content = json.dumps(payload, ensure_ascii=False)


class _SanitizeChoice:
    def __init__(self, payload):
        self.message = _SanitizeMessage(payload)


class _SanitizeCompletions:
    def __init__(self, payload):
        self.payload = payload

    def create(self, **kwargs):
        return type("R", (), {"choices": [_SanitizeChoice(self.payload)],
                              "usage": type("U", (), {"total_tokens": 42})()})()


class _SanitizeChat:
    def __init__(self, payload):
        self.completions = _SanitizeCompletions(payload)


class _SanitizeClient:
    def __init__(self, payload):
        self.chat = _SanitizeChat(payload)


class TestRunAiSanitizer:
    def _df(self):
        return pd.DataFrame([
            {"Run": "CRR1", "Location": pd.NA, "Tissue": "leaf",
             "Age_GrowthStage": "3 | Young Fruit", "CollectionDate": "2021-05"},
            {"Run": "CRR2", "Location": "China:Yinchuan", "Tissue": pd.NA,
             "Age_GrowthStage": pd.NA, "CollectionDate": pd.NA},
        ])

    def test_conflicts_keep_original_and_are_logged(self, tmp_path):
        payload = {"Location": "China, Ningxia, Yinchuan", "Tissue": "root",
                   "Age_GrowthStage": "3 years, young fruit stage"}
        out = _search.run_ai_sanitizer(self._df(), "sk-test", "https://api.deepseek.com",
                                       "deepseek-chat", out_dir=str(tmp_path),
                                       client=_SanitizeClient(payload))
        r1 = out[out["Run"] == "CRR1"].iloc[0]
        r2 = out[out["Run"] == "CRR2"].iloc[0]
        assert r1["Location"] == "China, Ningxia, Yinchuan_AI"
        assert r1["Tissue"] == "leaf"                       # AI 想改 root -> 保留 leaf 并记冲突
        assert "Tissue" in str(r1["AI_Conflict"])
        assert r1["CollectionDate"] == "2021-05"            # 洗髓不碰采集日期
        assert r2["Location"] == "China, Ningxia, Yinchuan_AI"  # 省市补全 + 单后缀
        assert r2["Tissue"] == "root"
        lines = (tmp_path / "ai_conflicts.jsonl").read_text(encoding="utf-8").splitlines()
        assert any("CRR1" in ln or "Tissue" in ln for ln in lines)


# ── AI 溯源日志 (ai_trace): 与平台 ai_log.jsonl 同 schema ──────────────────

class TestAiTrace:
    def teardown_method(self):
        _trace.set_default_state_dir(None)

    def test_trace_writes_platform_compatible_event(self, tmp_path):
        _trace.set_default_state_dir(tmp_path)
        p = _trace.trace("metadata_sanitizer", model="deepseek-chat",
                         prompt_text="prompt-x", rows=3, changed_rows=2)
        assert p is not None
        rec = json.loads(Path(p).read_text(encoding="utf-8").splitlines()[-1])
        assert rec["op"] == "metadata_sanitizer"
        assert rec["model"] == "deepseek-chat"
        assert rec["prompt"]["version"] == _trace.PROMPT_VERSIONS["metadata_sanitizer"]
        assert rec["source"] == "public_metadata_pipeline"
        assert rec["rows"] == 3
        assert rec["ok"] is True

    def test_auto_state_resolves_project_90handoff(self, tmp_path):
        (tmp_path / "90_Handoff").mkdir()
        out = tmp_path / "data" / "work" / "search"
        out.mkdir(parents=True)
        assert _trace.auto_state(out) == tmp_path / "90_Handoff" / "state"

    def test_auto_state_falls_back_to_outdir(self, tmp_path):
        out = tmp_path / "standalone"
        out.mkdir()
        assert _trace.auto_state(out) == out.resolve()

    def test_no_state_dir_is_noop(self):
        assert _trace.trace("metadata_summary") is None

    def test_prompt_sha_is_stable(self):
        assert _trace.prompt_sha("abc") == _trace.prompt_sha("abc")
        assert _trace.prompt_sha("abc") != _trace.prompt_sha("abd")
        assert _trace.prompt_sha(None) == ""


# ── fill_ai_blanks: 只填空 + Location 统一 _AI ─────────────────────────────

class TestFillAiBlanks:
    def test_fills_blanks_only_and_marks_location(self):
        df = pd.DataFrame([{"Run": "CRR1", "Location": "China, Ningxia, Yinchuan_AI",
                            "Tissue": "leaf", "Age_GrowthStage": pd.NA}])
        filled = _info.fill_ai_blanks(df, 0, {
            "Location": "China, Sichuan, Chengdu",
            "Tissue": "root",
            "Age_GrowthStage": "seedling",
        }, ["Location", "Tissue", "Age_GrowthStage"])
        assert filled == ["Age_GrowthStage"]
        assert df.at[0, "Location"] == "China, Ningxia, Yinchuan_AI"   # 已有真值不动
        assert df.at[0, "Tissue"] == "leaf"
        assert df.at[0, "Age_GrowthStage"] == "seedling"

    def test_new_location_gets_single_ai_suffix(self):
        df = pd.DataFrame([{"Run": "CRR1", "Location": pd.NA}])
        _info.fill_ai_blanks(df, 0, {"Location": "China:Yunnan:kunming"}, ["Location"])
        assert df.at[0, "Location"] == "China, Yunnan, kunming_AI"


# ── 冲突落盘 + 正交化防火墙 (非破坏写入) ──────────────────────────────────

class TestConflictLogAndFirewall:
    def test_write_conflicts_writes_csv_and_jsonl(self, tmp_path):
        n = _policy.write_conflicts(tmp_path, [
            {"run": "SRR1", "field": "CollectionDate", "old": "2021-05",
             "new": "2022-03", "action": "conflict", "stage": "info_firewall"}])
        assert n == 1
        assert "CollectionDate" in (tmp_path / "ai_conflicts.csv").read_text(encoding="utf-8-sig")
        rec = json.loads((tmp_path / "ai_conflicts.jsonl").read_text(encoding="utf-8").splitlines()[-1])
        assert rec["field"] == "CollectionDate" and rec["run"] == "SRR1"

    def test_write_conflicts_is_noop_without_records(self, tmp_path):
        assert _policy.write_conflicts(tmp_path, []) == 0
        assert not (tmp_path / "ai_conflicts.csv").exists()

    def test_firewall_keeps_truth_and_logs_conflict(self, tmp_path):
        payload = {"CollectionDate": "2022-03", "Location": "China, Yunnan, Kunming",
                   "ScientificName": "Malus pumila var. domestica"}
        df = pd.DataFrame([{"Run": "SRR9", "CollectionDate": "2021-05",
                            "Location": "China, Yunnan_AI",
                            "ScientificName": "Malus pumila"}])
        stats = _info.orthogonalize_ai_columns(df, _SanitizeClient(payload),
                                               "deepseek-chat", out_dir=str(tmp_path))
        assert df.at[0, "CollectionDate"] == "2021-05"              # 采集时间不被覆盖
        assert df.at[0, "Location"] == "China, Yunnan, Kunming_AI"  # 补全 + 单后缀
        assert df.at[0, "ScientificName"] == "Malus pumila"         # 学名冲突保留旧值
        assert "CollectionDate" in str(df.at[0, "AI_Conflict"])
        assert stats["conflict_count"] == 2
        assert stats["errors"] == 0
        assert (tmp_path / "ai_conflicts.csv").exists()

    def test_firewall_tolerates_client_errors(self, tmp_path):
        class _Boom:
            class chat:
                class completions:
                    @staticmethod
                    def create(**kwargs):
                        raise RuntimeError("api down")
        df = pd.DataFrame([{"Run": "SRR1", "CollectionDate": "2021-05"}])
        stats = _info.orthogonalize_ai_columns(df, _Boom(), "deepseek-chat",
                                               out_dir=str(tmp_path))
        assert stats["errors"] == 1 and stats["changed_rows"] == 0
        assert df.at[0, "CollectionDate"] == "2021-05"
