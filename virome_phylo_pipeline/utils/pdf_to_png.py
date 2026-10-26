#!/usr/bin/env python3
"""pdf_to_png.py — 第三方工具直写 PDF 的 PNG 兜底渲染 (PyMuPDF / fitz)。

背景: treetime / SpreaD3 / BEAST / 外部 R 脚本等自己写 PDF, 我们改不了它内部;
而报告要内嵌 PNG (`report_builder` 的 png>pdf 优先级)。本工具把这些 PDF
**渲染成同主名 PNG** (默认 200 dpi、白底)。

口径:
  · **只转「没有同名 PNG」的 PDF** —— 幂等; 已有 PNG 的默认**不覆盖**
    (已有的 PNG 多为 matplotlib 直出 300 dpi, 比 PDF 渲染更清晰)。
  · PyMuPDF (`import fitz`) 不可用时打印一行说明并**优雅跳过**, 绝不抛错、
    绝不中断管线。
  · 只取第 1 页 (多页 PDF 的完整原件仍保留, 报告层无需逐页 PNG)。
  · 空/近空 PDF (< 256 B, 历史上有 1 字节空图事故) 直接跳过。
  · 本模块对**第三方兜底**; 我们自己的 matplotlib 图请用 `utils.plot_io.save_both`
    双写 (矢量直出, 不经过栅格化)。

CLI:
    python -m utils.pdf_to_png IMG.pdf [DIR ...] [--dpi 200] [--overwrite] [-q]
"""
from __future__ import annotations

import os
import sys

__all__ = ["fitz_available", "pdf_to_png", "convert_pdfs", "convert_dir",
           "convert_paths", "main"]

DEFAULT_DPI = 200
MIN_PDF_BYTES = 256          # 空/近空 PDF 跳过阈值 (历史事故: 1 字节 PDF)

_FITZ_STATE = {"tried": False, "mod": None}


def _load_fitz():
    """惰性导入 fitz; 不可用时打印一行说明并返回 None (只打一次)。"""
    if not _FITZ_STATE["tried"]:
        _FITZ_STATE["tried"] = True
        try:
            import fitz  # noqa: F401  (PyMuPDF)
            _FITZ_STATE["mod"] = fitz
        except Exception:  # noqa: BLE001 — 缺依赖不是错误, 只是跳过
            _FITZ_STATE["mod"] = None
    return _FITZ_STATE["mod"]


def fitz_available() -> bool:
    """PyMuPDF 是否可用 (供调用方在收尾处判断, 不产生副作用)。"""
    return _load_fitz() is not None


def _png_sibling(pdf_path: str) -> str:
    return os.path.splitext(str(pdf_path))[0] + ".png"


def _emitter(log):
    """把 log 归一成一个 callable —— 兼容 None / 函数 / LogCollector(.emit/.info)。"""
    if log is None:
        return None
    if callable(log):
        return log
    for attr in ("emit", "info", "warning"):
        fn = getattr(log, attr, None)
        if callable(fn):
            return fn
    return None


def _resolve_emitter(log, quiet: bool):
    """显式给了 log 就听它的 (管线日志里要留痕); 没给且 quiet=True 则静默。"""
    e = _emitter(log)
    if e is not None:
        return e
    return (lambda _m: None) if quiet else (lambda m: print(m))


def pdf_to_png(pdf_path: str, png_path: str = None, dpi: int = DEFAULT_DPI,
               overwrite: bool = False, log=None, quiet: bool = False) -> str:
    """把一个 PDF 第 1 页渲染成 PNG。

    返回 PNG 路径; 跳过或失败返回 None (**不抛错**)。
    """
    emit = _resolve_emitter(log, quiet)
    pdf_path = str(pdf_path)
    png_path = str(png_path) if png_path else _png_sibling(pdf_path)

    if not os.path.exists(pdf_path):
        return None
    if os.path.exists(png_path) and not overwrite:
        return None                                  # 幂等: 已有 PNG, 不覆盖
    try:
        if os.path.getsize(pdf_path) < MIN_PDF_BYTES:
            emit(f"[pdf_to_png] 跳过 (空/近空 PDF, <{MIN_PDF_BYTES} B): {pdf_path}")
            return None
    except OSError:
        return None

    fitz = _load_fitz()
    if fitz is None:
        emit("[pdf_to_png] PyMuPDF (fitz) 不可用 → 跳过 PDF→PNG 兜底渲染 "
             "(安装: pip install PyMuPDF)。其余管线不受影响。")
        return None

    try:
        d = os.path.dirname(os.path.abspath(png_path))
        if d:
            os.makedirs(d, exist_ok=True)
        zoom = max(float(dpi), 1.0) / 72.0
        mat = fitz.Matrix(zoom, zoom)                # 白底: PDF 页自身背景
        with fitz.open(pdf_path) as doc:
            if doc.page_count < 1:                   # pragma: no cover
                emit(f"[pdf_to_png] 跳过 (0 页): {pdf_path}")
                return None
            page = doc.load_page(0)
            pix = page.get_pixmap(matrix=mat, alpha=False)
            pix.save(png_path)
        if os.path.getsize(png_path) <= 0:           # pragma: no cover
            os.remove(png_path)
            return None
        return png_path
    except Exception as e:  # noqa: BLE001 — 兜底渲染失败只告警, 不中断管线
        emit(f"[pdf_to_png] 渲染失败 (跳过): {pdf_path} — {e}")
        return None


def has_png(pdf_path: str) -> bool:
    """同名 PNG 是否已存在。"""
    return os.path.exists(_png_sibling(pdf_path))


def convert_pdfs(pdfs, dpi: int = DEFAULT_DPI, overwrite: bool = False,
                 log=None, quiet: bool = False) -> list:
    """批量转换一组 PDF 路径; **只转没有同名 PNG 的**。返回新建 PNG 列表。"""
    emit = _resolve_emitter(log, quiet)
    made = []
    for p in pdfs:
        if has_png(p) and not overwrite:
            continue
        got = pdf_to_png(p, dpi=dpi, overwrite=overwrite, log=log, quiet=quiet)
        if got:
            made.append(got)
            emit(f"[pdf_to_png] {os.path.basename(p)} → {os.path.basename(got)} "
                 f"({dpi} dpi)")
    return made


def convert_dir(directory: str, dpi: int = DEFAULT_DPI, recursive: bool = True,
                overwrite: bool = False, log=None, quiet: bool = True,
                patterns=None) -> list:
    """把目录下「没有同名 PNG」的 PDF 全部兜底渲染。

    幂等: 重复调用第二次返回空列表。任何单文件失败都只告警, 不抛出。
    """
    if not directory or not os.path.isdir(directory):
        return []
    pats = tuple(patterns or (".pdf",))
    pdfs = []
    if recursive:
        for dp, _dn, fn in os.walk(directory):
            for f in fn:
                if f.lower().endswith(pats):
                    pdfs.append(os.path.join(dp, f))
    else:
        for f in os.listdir(directory):
            if f.lower().endswith(pats):
                pdfs.append(os.path.join(directory, f))
    return convert_pdfs(sorted(pdfs), dpi=dpi, overwrite=overwrite,
                        log=log, quiet=quiet)


def convert_paths(paths, dpi: int = DEFAULT_DPI, recursive: bool = True,
                  overwrite: bool = False, log=None, quiet: bool = True) -> list:
    """路径 (文件或目录) 混合入口 —— stage 收尾最省事的调用方式。"""
    made = []
    pdfs = []
    for p in paths:
        p = str(p)
        if os.path.isdir(p):
            made += convert_dir(p, dpi=dpi, recursive=recursive,
                                overwrite=overwrite, log=log, quiet=quiet)
        elif p.lower().endswith(".pdf"):
            if has_png(p) and not overwrite:
                continue
            pdfs.append(p)
    if pdfs:
        made += convert_pdfs(pdfs, dpi=dpi, overwrite=overwrite,
                             log=log, quiet=quiet)
    return made


def stage_finalize(output_dir: str, log=None, dpi: int = DEFAULT_DPI,
                   quiet: bool = True) -> list:
    """stage 收尾约定入口: 兜底渲染 output_dir 下缺 PNG 的 PDF。

    永不抛错 —— 收尾环节挂掉不该让 stage 失败。返回新建 PNG 列表。
    log 可为 None / callable / LogCollector (带 .emit)。
    """
    emit = _emitter(log)
    try:
        return convert_dir(output_dir, dpi=dpi, log=emit, quiet=quiet)
    except Exception as e:  # noqa: BLE001
        if emit:
            emit(f"[pdf_to_png] 收尾兜底异常 (忽略): {e}")
        return []


def main(argv=None) -> int:
    import argparse
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("paths", nargs="+", help="PDF 文件或目录 (目录则递归兜底)")
    ap.add_argument("--dpi", type=int, default=DEFAULT_DPI)
    ap.add_argument("--overwrite", action="store_true",
                    help="已存在同名 PNG 时也重渲染 (默认不覆盖)")
    ap.add_argument("-q", "--quiet", action="store_true")
    a = ap.parse_args(argv)
    if not fitz_available():
        print("PyMuPDF (fitz) 不可用 → 跳过 PDF→PNG 兜底渲染。"
              "安装: pip install PyMuPDF", file=sys.stderr)
        return 0                                     # 不是错误, 不返回非零
    made = convert_paths(a.paths, dpi=a.dpi, overwrite=a.overwrite, quiet=a.quiet)
    if not a.quiet:
        print(f"[pdf_to_png] 新建 {len(made)} 个 PNG")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
