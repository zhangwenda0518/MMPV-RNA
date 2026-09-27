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
