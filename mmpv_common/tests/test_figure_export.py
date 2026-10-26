# -*- coding: utf-8 -*-
"""图件四格式导出层测试 (2026-10-20)。

覆盖 `mmpv_common.figure_export`:
  · target_path / figure_paths —— 四格式路径规则 (hires 进独立子目录)
  · export_figure —— SVG/PNG/PDF/hires 全部落盘, 尺寸与 dpi 正确
  · 全局 savefig 钩子 —— 老写法 (只写 png / 只写 pdf / 双写) 都补齐四格式,
    且同一 stem 不重复渲染; PdfPages 与文件对象原样透传
  · 环境变量开关 —— MMPV_FIG_EXPORT=0 关闭, MMPV_FIG_FORMATS 裁剪
  · 失败隔离 —— 某格式写失败不抛错, 其他格式照常落盘
"""
import os
import sys

import pytest

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(_HERE)
if _REPO not in sys.path:
    sys.path.insert(0, _REPO)

matplotlib = pytest.importorskip("matplotlib")
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from mmpv_common import figure_export as fe  # noqa: E402

PNG_MAGIC = b"\x89PNG\r\n\x1a\n"
PDF_MAGIC = b"%PDF-"


@pytest.fixture(autouse=True)
def _clean_env(monkeypatch):
    """每个用例从默认配置出发, 并保证钩子状态干净。"""
    for k in ("MMPV_FIG_EXPORT", "MMPV_FIG_FORMATS", "MMPV_FIG_PNG_DPI",
              "MMPV_FIG_HIRES_DPI", "MMPV_FIG_HIRES_DIR"):
        monkeypatch.delenv(k, raising=False)
    yield
    plt.close("all")


def _fig(w=2.0, h=2.0):
    fig, ax = plt.subplots(figsize=(w, h))
    ax.plot([1, 2, 3], [1, 4, 9])
    return fig


def _head(path, n=8):
    with open(path, "rb") as f:
        return f.read(n)


# ═══════════════════════════════════════════════════════════════════
# 路径规则
# ═══════════════════════════════════════════════════════════════════

def test_target_path_four_formats(tmp_path):
    base = str(tmp_path / "sub" / "fig1.pdf")
    assert fe.target_path(base, "svg") == str(tmp_path / "sub" / "fig1.svg")
    assert fe.target_path(base, "png") == str(tmp_path / "sub" / "fig1.png")
    assert fe.target_path(base, "pdf") == str(tmp_path / "sub" / "fig1.pdf")
    assert fe.target_path(base, "hires") == str(
        tmp_path / "sub" / fe.HIRES_SUBDIR / "fig1.png")


def test_target_path_is_extension_agnostic(tmp_path):
    """传 .png / .svg / 无扩展名, 都应得到同一组主名。"""
    for name in ("fig9.png", "fig9.svg", "fig9"):
        got = fe.figure_paths(str(tmp_path / name))
        assert got["pdf"] == str(tmp_path / "fig9.pdf")
        assert got["hires"] == str(tmp_path / fe.HIRES_SUBDIR / "fig9.png")


def test_target_path_rejects_unknown_format(tmp_path):
    with pytest.raises(ValueError):
        fe.target_path(str(tmp_path / "a.pdf"), "tiff")


# ═══════════════════════════════════════════════════════════════════
# export_figure
# ═══════════════════════════════════════════════════════════════════

def test_export_figure_writes_all_four(tmp_path):
    out = fe.export_figure(_fig(), str(tmp_path / "e.pdf"))
    for fmt in fe.ALL_FORMATS:
        p = out[fmt]
        assert p and os.path.exists(p), f"缺 {fmt}"
        assert os.path.getsize(p) > 1000, f"{p} 过小"
    assert _head(out["png"]) == PNG_MAGIC
    assert _head(out["hires"]) == PNG_MAGIC
    assert _head(out["pdf"], 5) == PDF_MAGIC
    assert out["svg"].endswith(".svg")
    assert os.path.basename(os.path.dirname(out["hires"])) == fe.HIRES_SUBDIR


def test_export_figure_dpi_defaults(tmp_path):
    """PNG=300 dpi, hires=1200 dpi (2 in 画布 → 600 / 2400 px)。"""
    from matplotlib.image import imread
    out = fe.export_figure(_fig(2, 2), str(tmp_path / "d.pdf"))
    assert imread(out["png"]).shape[:2] == (600, 600)
    assert imread(out["hires"]).shape[:2] == (2400, 2400)


def test_export_figure_dpi_override(tmp_path):
    from matplotlib.image import imread
    out = fe.export_figure(_fig(1, 1), str(tmp_path / "o.pdf"),
                           png_dpi=150, hires_dpi=600)
    assert imread(out["png"]).shape[:2] == (150, 150)
    assert imread(out["hires"]).shape[:2] == (600, 600)


def test_export_figure_pdf_keeps_caller_dpi(tmp_path):
    """PDF 的 dpi 用调用方参数 (历史口径不变), 不被 png_dpi 污染。"""
    out = fe.export_figure(_fig(), str(tmp_path / "p.pdf"), dpi=150)
    assert out["pdf"] and os.path.exists(out["pdf"])


def test_export_figure_formats_subset(tmp_path):
    out = fe.export_figure(_fig(), str(tmp_path / "s.pdf"),
                           formats=("png", "pdf"))
    assert out["png"] and out["pdf"]
    assert out["svg"] is None and out["hires"] is None
    assert not os.path.exists(str(tmp_path / fe.HIRES_SUBDIR / "s.png"))


def test_export_figure_path_overrides(tmp_path):
    custom = str(tmp_path / "自定义.png")
    out = fe.export_figure(_fig(), str(tmp_path / "x.pdf"),
                           path_overrides={"png": custom})
    assert out["png"] == custom and os.path.exists(custom)


def test_export_figure_overwrite_false_skips(tmp_path):
    p = str(tmp_path / "k.pdf")
    fe.export_figure(_fig(), p)
    mtime = os.path.getmtime(str(tmp_path / "k.png"))
    fe.export_figure(_fig(), p, overwrite=False)
    assert os.path.getmtime(str(tmp_path / "k.png")) == mtime


def test_export_figure_hires_degrades_on_huge_figure(tmp_path, monkeypatch):
    """超大画布 → 高清 dpi 自动降级, 不 OOM。"""
    monkeypatch.setattr(fe, "MAX_HIRES_PIXELS", 1_000_000)
    assert fe._effective_hires_dpi(_fig(10, 10), 1200) <= fe.HIRES_FALLBACK_DPI


def test_export_figure_failure_is_isolated(tmp_path, monkeypatch):
    """单格式失败只告警, 其余格式照常落盘。"""
    real = fe._save_one

    def flaky(fig, path, fmt, *a, **kw):
        if fmt == "svg":
            raise RuntimeError("模拟 SVG 失败")
        return real(fig, path, fmt, *a, **kw)

    monkeypatch.setattr(fe, "_save_one", flaky)
    out = fe.export_figure(_fig(), str(tmp_path / "f.pdf"))
    assert out["svg"] is None
    assert out["png"] and out["pdf"] and out["hires"]


# ═══════════════════════════════════════════════════════════════════
# 环境变量开关
# ═══════════════════════════════════════════════════════════════════

def test_env_disable_only_affects_hook(tmp_path, monkeypatch):
    """`MMPV_FIG_EXPORT=0` 关闭的是**自动补齐**; 显式 export_figure 仍然照写。"""
    monkeypatch.setenv("MMPV_FIG_EXPORT", "0")
    out = fe.export_figure(_fig(), str(tmp_path / "off.pdf"))
    assert out["png"] and os.path.exists(out["png"])

    # 显式 formats=() 才是「什么都不写」
    out2 = fe.export_figure(_fig(), str(tmp_path / "none.pdf"), formats=())
    assert all(v is None for v in out2.values())
    assert not os.path.exists(str(tmp_path / "none.png"))


def test_env_formats(tmp_path, monkeypatch):
    monkeypatch.setenv("MMPV_FIG_FORMATS", "png,pdf")
    out = fe.export_figure(_fig(), str(tmp_path / "env.pdf"))
    assert out["png"] and out["pdf"]
    assert out["svg"] is None and out["hires"] is None


def test_env_hires_dir(tmp_path, monkeypatch):
    monkeypatch.setenv("MMPV_FIG_HIRES_DIR", "_print")
    out = fe.export_figure(_fig(), str(tmp_path / "hd.pdf"))
    assert os.path.basename(os.path.dirname(out["hires"])) == "_print"


# ═══════════════════════════════════════════════════════════════════
# 全局 savefig 钩子
# ═══════════════════════════════════════════════════════════════════

@pytest.fixture
def hooked():
    """安装钩子, 用例结束恢复原 savefig。"""
    fe.install_savefig_hook()
    assert fe.hook_installed()
    yield
    fe.uninstall_savefig_hook()


def test_hook_png_only_fills_others(tmp_path, hooked):
    fig = _fig()
    fig.savefig(str(tmp_path / "h.png"), dpi=300, bbox_inches="tight")
    for fmt in fe.ALL_FORMATS:
        assert os.path.exists(fe.target_path(str(tmp_path / "h.png"), fmt)), fmt


def test_hook_pdf_only_fills_others(tmp_path, hooked):
    fig = _fig()
    fig.savefig(str(tmp_path / "p.pdf"), dpi=300)
    for fmt in fe.ALL_FORMATS:
        assert os.path.exists(fe.target_path(str(tmp_path / "p.pdf"), fmt)), fmt


def test_hook_png_pdf_pair_does_not_rerender(tmp_path, hooked):
    """老写法 png + pdf 双写: 第二次调用不应重写已存在的兄弟文件。"""
    fig = _fig()
    fig.savefig(str(tmp_path / "pair.png"), dpi=300)
    svg_mtime = os.path.getmtime(str(tmp_path / "pair.svg"))
    hires_mtime = os.path.getmtime(
        str(tmp_path / fe.HIRES_SUBDIR / "pair.png"))
    fig.savefig(str(tmp_path / "pair.pdf"), dpi=300)
    assert os.path.getmtime(str(tmp_path / "pair.svg")) == svg_mtime
    assert os.path.getmtime(
        str(tmp_path / fe.HIRES_SUBDIR / "pair.png")) == hires_mtime


def test_hook_passthrough_for_pdfpages(tmp_path, hooked):
    """PdfPages 是文件对象 → 原样透传, 不得多出任何格式。"""
    from matplotlib.backends.backend_pdf import PdfPages
    import fitz

    fig = _fig()
    with PdfPages(str(tmp_path / "multi.pdf")) as pp:
        pp.savefig(fig)
        pp.savefig(fig)

    assert sorted(os.listdir(tmp_path)) == ["multi.pdf"]
    assert fitz.open(str(tmp_path / "multi.pdf")).page_count == 2


def test_hook_passthrough_for_file_object(tmp_path, hooked):
    fig = _fig()
    buf = tmp_path / "buf.png"
    with open(buf, "wb") as fh:
        fig.savefig(fh, format="png", dpi=100)
    assert buf.stat().st_size > 1000
    assert sorted(os.listdir(tmp_path)) == ["buf.png"]


def test_hook_disabled_by_env(tmp_path, hooked, monkeypatch):
    monkeypatch.setenv("MMPV_FIG_EXPORT", "0")
    fig = _fig()
    fig.savefig(str(tmp_path / "d.png"), dpi=100)
    assert sorted(os.listdir(tmp_path)) == ["d.png"]


def test_hook_install_is_idempotent():
    assert fe.install_savefig_hook() is True
    assert fe.install_savefig_hook() is True
    fe.uninstall_savefig_hook()


# ═══════════════════════════════════════════════════════════════════
# 与报告层的兼容性
# ═══════════════════════════════════════════════════════════════════

def test_hires_subdir_not_collected_by_report_layer(tmp_path):
    """报告层按 (目录, 主名) 去重; 高清 PNG 在独立子目录 → 不得重复内嵌。

    这条约束是 `HIRES_SUBDIR` 存在的理由: 若把 1200 dpi 副本写成
    `<stem>_1200dpi.png` 放在同目录, 报告会把它当第二张图。
    """
    out = fe.export_figure(_fig(), str(tmp_path / "fig1.pdf"))
    assert os.path.exists(out["hires"])

    # 复刻 report_builder.collect_plots 的去重口径: key = (parent, stem)
    by_key = {}
    for name in sorted(os.listdir(tmp_path)):
        p = tmp_path / name
        if not p.is_file() or p.suffix.lower() not in (".png", ".pdf", ".svg"):
            continue
        by_key.setdefault((str(p.parent), p.stem.lower()), []).append(name)

    assert len(by_key) == 1, f"同目录同主名被拆成多组: {by_key}"
    assert list(by_key.values())[0] == ["fig1.pdf", "fig1.png", "fig1.svg"]
    assert os.path.dirname(out["hires"]) != str(tmp_path)
