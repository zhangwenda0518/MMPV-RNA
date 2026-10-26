#!/usr/bin/env python3
"""plot_io.py — 图件双格式写出统一收口 (PNG + PDF 同主名, 仅扩展名不同)。

口径 (2026-10-16 裁定):
  · 管道所有图件一律**同时**产出 PNG + PDF, 主名相同、仅扩展名不同;
    报告优先内嵌 PNG (report_builder 既有 png>pdf 优先级), 矢量 PDF 供出版。
  · PNG 规范: 300 dpi、白底 (matplotlib 直出 PNG, 不嵌 JPEG)。
  · **不改图的视觉设计**, 也**不改既有 PDF 输出** —— save_both 对 PDF 完全按
    调用方原参数写出 (逐字节口径与原 `fig.savefig(out_pdf, ...)` 一致),
    只是把同一 figure 再写一份 PNG。多页 PdfPages 保持 PDF-only, 不走这里。

用法:
    from utils.plot_io import save_both          # 脚本直跑时 `from plot_io import ...`

    # 旧: fig.savefig(out_pdf, dpi=300, format='pdf', bbox_inches='tight')
    # 新:
    save_both(fig, out_pdf, dpi=300, format='pdf', bbox_inches='tight')

⚠ **format= 陷阱**: `fig.savefig("x.png", format="pdf")` 会写出一个**内容是
PDF、扩展名却是 .png** 的文件 (报告会把它当图片内嵌 → 空白/解析失败)。
save_both 写 PNG 时会剥掉调用方的 `format=`, 强制 `format="png"`。

⚠ **dpi 口径**: 调用方 kw 里的 `dpi` 只作用于 PDF (保持历史输出不变);
PNG 单独用 `png_dpi` (默认 300)。历史上有 dpi=150 的 PDF, 不代表 PNG 也要 150。
"""
from __future__ import annotations

import os
import sys

__all__ = ["png_sibling", "pdf_sibling", "save_both", "ensure_png_pair"]

DEFAULT_PNG_DPI = 300
DEFAULT_PNG_FACECOLOR = "white"

# 仓库根 (本文件在 <root>/virome_phylo_pipeline/utils/) —— 供共享导出层导入
_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

_UNSET = object()
_EXPORT_FN = _UNSET


def png_sibling(path: str) -> str:
    """`a/b/foo.pdf` → `a/b/foo.png`; 无扩展名 → 直接加 `.png`。"""
    root, ext = os.path.splitext(str(path))
    if ext.lower() == ".png":
        return str(path)
    if not ext:
        return str(path) + ".png"
    return root + ".png"


def pdf_sibling(path: str) -> str:
    """`a/b/foo.png` → `a/b/foo.pdf`; 无扩展名 → 直接加 `.pdf`。"""
    root, ext = os.path.splitext(str(path))
    if ext.lower() == ".pdf":
        return str(path)
    if not ext:
        return str(path) + ".pdf"
    return root + ".pdf"


def save_both(fig, path, png=None, *, png_dpi=DEFAULT_PNG_DPI,
              png_facecolor=DEFAULT_PNG_FACECOLOR, **savefig_kw):
    """把同一 figure 写成 PDF + PNG (同主名, 仅扩展名不同)。

    Parameters
    ----------
    fig : matplotlib.figure.Figure
    path : str — PDF 路径 (或任意路径: `.png` 会给主名补 `.pdf`)
    png : str, optional — 显式 PNG 路径; 默认 = path 主名 + `.png`
    png_dpi : int — PNG 的 dpi (默认 300, 与 scientific-visualization 规范一致)
    png_facecolor : str — PNG 底色 (默认 white)
    **savefig_kw — 原样传给 PDF 那次 savefig (dpi / bbox_inches / format /
        facecolor / transparent …) —— 保证既有 PDF 输出不变。

    Returns
    -------
    dict — {"pdf": <pdf 路径>, "png": <png 路径>}; PNG 写失败时 "png" 为 None
        (只告警, 不抛错: PDF 已落盘, 不能因为 PNG 拖垮整个 stage)。
    """
    pdf_path = str(path)
    if os.path.splitext(pdf_path)[1].lower() == ".png":
        pdf_path = pdf_sibling(pdf_path)
    elif not os.path.splitext(pdf_path)[1]:
        pdf_path = pdf_path + ".pdf"

    png_path = str(png) if png else png_sibling(pdf_path)

    for p in (pdf_path, png_path):
        d = os.path.dirname(os.path.abspath(p))
        if d:
            os.makedirs(d, exist_ok=True)

    # 2026-10-20: 统一走 mmpv_common.figure_export —— 除 PDF/PNG 外再补
    # SVG 与 _hires_1200dpi/PNG。PDF/PNG 的路径与参数口径**逐字不变**
    # (pdf 用调用方原参, png 用 png_dpi/facecolor), 只是把两次 savefig
    # 收进导出层, 由它按 formats 顺序写出四种格式。
    export = _load_exporter()
    if export is not None:
        out = export(
            fig, pdf_path,
            path_overrides={"png": png_path},
            png_dpi=png_dpi,
            png_facecolor=png_facecolor,
            **savefig_kw,
        )
        return {"pdf": out.get("pdf"), "png": out.get("png"),
                "svg": out.get("svg"), "hires": out.get("hires")}

    # 兜底: 仓库根不可导入 mmpv_common 时保持旧行为 (PNG + PDF 双写)
    # ① PDF: 完全按调用方参数 (含历史 format='pdf' / dpi=150 等) —— 口径不变
    fig.savefig(pdf_path, **savefig_kw)

    # ② PNG: 剥掉 format= 强制 png (见模块 docstring 的 format 陷阱),
    #    dpi/facecolor 用 PNG 专属值, 其余 (bbox_inches 等) 与 PDF 同参。
    png_kw = {k: v for k, v in savefig_kw.items() if k != "format"}
    png_kw["format"] = "png"
    png_kw["dpi"] = png_dpi
    png_kw.setdefault("facecolor", png_facecolor)
    try:
        fig.savefig(png_path, **png_kw)
    except Exception as e:  # noqa: BLE001 — PNG 失败不得中断管线
        print(f"[plot_io] 警告: PNG 写出失败 ({png_path}): {e}", file=sys.stderr)
        return {"pdf": pdf_path, "png": None}
    return {"pdf": pdf_path, "png": png_path}


def _load_exporter():
    """惰性取 mmpv_common.figure_export.export_figure; 不可用则返回 None。"""
    global _EXPORT_FN
    if _EXPORT_FN is not _UNSET:
        return _EXPORT_FN
    try:
        if _REPO_ROOT not in sys.path:
            sys.path.insert(0, _REPO_ROOT)
        from mmpv_common.figure_export import export_figure
        _EXPORT_FN = export_figure
    except Exception:  # noqa: BLE001 — 缺共享层不致命, 走旧双写兜底
        _EXPORT_FN = None
    return _EXPORT_FN


def ensure_png_pair(pdf_path: str, png_path: str = None) -> dict:
    """仅补算缺的 PNG 路径 (不渲染); 供调用方自行判断/记录。

    第三方工具直写 PDF 的兜底渲染请用 `utils.pdf_to_png.pdf_to_png`。
    """
    pdf_path = pdf_sibling(pdf_path) if not os.path.splitext(str(pdf_path))[1] \
        else str(pdf_path)
    png_path = str(png_path) if png_path else png_sibling(pdf_path)
    return {"pdf": pdf_path, "png": png_path,
            "has_pdf": os.path.exists(pdf_path), "has_png": os.path.exists(png_path)}
