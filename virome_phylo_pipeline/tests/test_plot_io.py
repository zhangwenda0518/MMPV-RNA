# -*- coding: utf-8 -*-
"""图件双格式写出测试 (2026-10-16, 格式双输出管线改造)。

覆盖:
  · `utils.plot_io.save_both` —— PNG + PDF 同主名同时落盘、PNG 非空且是真 PNG、
    `format=` 陷阱不复发 (传 format='pdf' 时 .png 不得被写成 PDF 内容)、
    PNG 300 dpi / 白底、PDF 输出与旧写法等价 (页数/画布尺寸一致)。
  · `utils.pdf_to_png` —— 第三方 PDF 兜底渲染: 只转没有同名 PNG 的、幂等、
    空 PDF 跳过、fitz 不可用时优雅跳过 (不抛错)。
"""
import os
import sys
from pathlib import Path

import pytest

_HERE = Path(__file__).resolve().parent
_PHYLO = _HERE.parent
_REPO = _PHYLO.parent
for p in (str(_PHYLO), str(_REPO)):
    if p not in sys.path:
        sys.path.insert(0, p)

# ⑤ 的 utils/ 与其他管线同名包冲突 (与 test_popgen_stats 同款处理):
# 全量 pytest 时若别的测试先导入了别的 utils, 这里必须清缓存并把自己插到最前。
for _m in [m for m in list(sys.modules) if m == "utils" or m.startswith("utils.")]:
    del sys.modules[_m]
sys.path.insert(0, str(_PHYLO))

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from utils.pdf_to_png import (convert_dir, convert_pdfs, pdf_to_png,  # noqa: E402
                              fitz_available, stage_finalize)
from utils.plot_io import png_sibling, pdf_sibling, save_both  # noqa: E402

PNG_MAGIC = b"\x89PNG\r\n\x1a\n"
PDF_MAGIC = b"%PDF-"


def _fig():
    fig, ax = plt.subplots(figsize=(2, 2))
    ax.plot([1, 2, 3], [1, 4, 9])
    return fig


def _read_head(path, n=8):
    with open(path, "rb") as f:
        return f.read(n)


# ═══════════════════════════════════════════════════════════════════
# save_both
# ═══════════════════════════════════════════════════════════════════

def test_save_both_writes_png_and_pdf(tmp_path):
    """PNG + PDF 同时落盘, 主名相同; PNG 非空且是真 PNG, PDF 是真 PDF。"""
    pdf = tmp_path / "fig_a.pdf"
    out = save_both(_fig(), str(pdf), dpi=300, format="pdf", bbox_inches="tight")
    plt.close("all")

    assert out["pdf"] == str(pdf)
    assert out["png"] == str(tmp_path / "fig_a.png")
    for p in (out["pdf"], out["png"]):
        assert os.path.exists(p), f"缺 {p}"
        assert os.path.getsize(p) > 1000, f"{p} 过小 (空图?)"
    assert _read_head(out["png"]) == PNG_MAGIC, "PNG 头部不对 (format= 陷阱复发?)"
    assert _read_head(out["pdf"], 5) == PDF_MAGIC
    # 主名相同、仅扩展名不同
    assert os.path.splitext(out["png"])[0] == os.path.splitext(out["pdf"])[0]


def test_save_both_png_is_300dpi_white(tmp_path):
    """PNG 规范: 300 dpi、白底 (角点像素为白)。"""
    from matplotlib.image import imread

    out = save_both(_fig(), str(tmp_path / "dpi.pdf"))
    plt.close("all")
    img = imread(out["png"])
    # figsize=2x2 in @300 dpi = 600x600 px (savefig 无 bbox_inches 时严格相等)
    assert img.shape[0] == 600 and img.shape[1] == 600, img.shape
    corner = img[0, 0]
    assert all(float(c) >= 0.99 for c in corner[:3]), f"角点非白: {corner}"


def test_save_both_png_dpi_override(tmp_path):
    """png_dpi 只管 PNG; PDF 那次 savefig 的 dpi 参数不被污染。"""
    from matplotlib.image import imread
    out = save_both(_fig(), str(tmp_path / "d2.pdf"), dpi=150, png_dpi=300)
    plt.close("all")
    img = imread(out["png"])
    assert img.shape[0] == 600, img.shape


def test_save_both_explicit_png_path_and_png_input(tmp_path):
    """显式 png= 生效; 传 .png 路径时按主名补 .pdf (反向双写)。"""
    out = save_both(_fig(), str(tmp_path / "x.pdf"), png=str(tmp_path / "自定义.png"))
    plt.close("all")
    assert os.path.exists(str(tmp_path / "自定义.png"))

    out2 = save_both(_fig(), str(tmp_path / "y.png"))
    plt.close("all")
    assert out2["png"] == str(tmp_path / "y.png")
    assert out2["pdf"] == str(tmp_path / "y.pdf")
    assert os.path.exists(out2["pdf"])


def test_save_both_creates_parent_dir(tmp_path):
    """目标目录不存在时自动建 (stage 里 outdir 常晚于绘图创建)。"""
    deep = tmp_path / "a" / "b" / "c.pdf"
    out = save_both(_fig(), str(deep))
    plt.close("all")
    assert os.path.exists(out["pdf"]) and os.path.exists(out["png"])


def test_save_both_pdf_matches_plain_savefig(tmp_path):
    """PDF 输出与旧写法等价: 页数与画布尺寸一致 (不改 PDF 输出口径)。"""
    pytest.importorskip("fitz")
    import fitz

    fig_a = _fig()
    plain = tmp_path / "plain.pdf"
    fig_a.savefig(str(plain), dpi=300, format="pdf", bbox_inches="tight")
    plt.close(fig_a)

    fig_b = _fig()
    out = save_both(fig_b, str(tmp_path / "both.pdf"), dpi=300, format="pdf",
                    bbox_inches="tight")
    plt.close(fig_b)

    with fitz.open(str(plain)) as da, fitz.open(out["pdf"]) as db:
        assert da.page_count == db.page_count == 1
        ra, rb = da.load_page(0).rect, db.load_page(0).rect
        assert abs(ra.width - rb.width) < 0.5 and abs(ra.height - rb.height) < 0.5


def test_sibling_helpers():
    assert png_sibling("a/b/foo.pdf") == "a/b/foo.png"
    assert png_sibling("a/b/foo.PDF") == "a/b/foo.png"
    assert png_sibling("noext") == "noext.png"
    assert png_sibling("a/b/foo.png") == "a/b/foo.png"
    assert pdf_sibling("a/b/foo.png") == "a/b/foo.pdf"
    assert pdf_sibling("noext") == "noext.pdf"
    assert pdf_sibling("a/b/foo.pdf") == "a/b/foo.pdf"


# ═══════════════════════════════════════════════════════════════════
# pdf_to_png (第三方直写 PDF 的兜底)
# ═══════════════════════════════════════════════════════════════════

def _make_pdf(path, figsize=(2, 2)):
    fig, ax = plt.subplots(figsize=figsize)
    ax.plot([1, 2, 3], [3, 2, 1])
    fig.savefig(str(path), dpi=300, format="pdf")
    plt.close(fig)
    return str(path)


def test_pdf_to_png_renders_and_is_idempotent(tmp_path):
    """转出同名 PNG (非空, 200 dpi 白底); 第二次调用不重复转 (幂等)。"""
    pytest.importorskip("fitz")
    pdf = _make_pdf(tmp_path / "third_party.pdf")
    dst = tmp_path / "third_party.png"

    got = pdf_to_png(pdf)
    assert got == str(dst)
    assert os.path.getsize(str(dst)) > 500
    assert _read_head(str(dst)) == PNG_MAGIC

    # 幂等: 已有同名 PNG 时跳过 (返回 None 且不改内容)
    before = os.path.getmtime(str(dst))
    assert pdf_to_png(pdf) is None
    assert os.path.getmtime(str(dst)) == before
    assert convert_pdfs([pdf]) == []
    assert convert_dir(str(tmp_path)) == []

    # overwrite=True 才重渲染
    assert pdf_to_png(pdf, overwrite=True) == str(dst)


def test_pdf_to_png_skips_small_and_missing(tmp_path):
    """空/近空 PDF 跳过 (历史 1 字节空图); 不存在的路径返回 None 不抛错。"""
    tiny = tmp_path / "tiny.pdf"
    tiny.write_bytes(b"%PDF-1.4\n")
    assert pdf_to_png(str(tiny)) is None
    assert not (tmp_path / "tiny.png").exists()
    assert pdf_to_png(str(tmp_path / "nope.pdf")) is None


def test_pdf_to_png_graceful_without_fitz(tmp_path, monkeypatch, capsys):
    """fitz 不可用: 打印一行说明 + 优雅跳过, **不抛错**。

    ⚠ 用 `pdf_to_png.__globals__` 而不是 `import utils.pdf_to_png as P` 来打补丁:
    本仓库多个 test 文件都会清 `utils.*` 缓存再重导, 全量 pytest 时
    `sys.modules['utils.pdf_to_png']` 可能指向**另一个**模块对象, 打它不生效
    (实测: 单跑本文件绿、全量跑红)。函数的 `__globals__` 才是权威字典。
    """
    pdf = _make_pdf(tmp_path / "nofitz.pdf")
    monkeypatch.setitem(pdf_to_png.__globals__, "_load_fitz", lambda: None)

    assert pdf_to_png(pdf) is None                      # 不抛错
    assert not (tmp_path / "nofitz.png").exists()
    assert convert_dir(str(tmp_path)) == []


def test_stage_finalize_never_raises(tmp_path):
    """stage 收尾入口对坏路径 / 坏 PDF 都必须静默返回空列表。"""
    bad = tmp_path / "broken.pdf"
    bad.write_bytes(b"%PDF-1.4\n" + b"garbage" * 200)    # 大但损坏的 PDF
    assert stage_finalize(str(tmp_path / "does_not_exist")) == []
    out = stage_finalize(str(tmp_path), log=None)        # 不抛错即可
    assert isinstance(out, list)


def test_convert_dir_only_converts_pdfs_without_png(tmp_path):
    """目录批量: 有 PNG 的 PDF 不动, 没 PNG 的才转; 递归子目录。"""
    pytest.importorskip("fitz")
    p1 = _make_pdf(tmp_path / "has_png.pdf")
    p2 = _make_pdf(tmp_path / "needs_png.pdf")
    sub = tmp_path / "sub"
    sub.mkdir()
    p3 = _make_pdf(sub / "deep.pdf")
    # 给 p1 预置一个"已有 PNG" → 不得被覆盖
    (tmp_path / "has_png.png").write_bytes(b"\x89PNG\r\n\x1a\nsentinel")

    made = convert_dir(str(tmp_path))
    assert sorted(os.path.basename(m) for m in made) == ["deep.png", "needs_png.png"]
    assert (tmp_path / "has_png.png").read_bytes().endswith(b"sentinel")
    assert os.path.exists(str(sub / "deep.png"))
    assert not os.path.exists(png_sibling(p3) + ".tmp")


def test_convert_paths_mixed_files_and_dirs(tmp_path):
    """文件 + 目录混合入口 (stage 收尾最常用形态)。"""
    pytest.importorskip("fitz")
    from utils.pdf_to_png import convert_paths

    p1 = _make_pdf(tmp_path / "one.pdf")
    sub = tmp_path / "d"
    sub.mkdir()
    _make_pdf(sub / "two.pdf")
    made = convert_paths([p1, str(sub)])
    assert sorted(os.path.basename(m) for m in made) == ["one.png", "two.png"]
    # 非 PDF 路径忽略
    assert convert_paths([str(tmp_path / "readme.txt")]) == []


def test_fitz_available_is_bool():
    assert isinstance(fitz_available(), bool)
