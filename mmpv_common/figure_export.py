#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""figure_export.py — 图件四格式统一导出 (SVG / PNG / PDF / PNG@1200 dpi)。

口径 (2026-10-20 裁定)
----------------------
仓库内**所有**出图脚本默认同时产出四种文件, 主名相同:

    <stem>.svg          矢量, 文本可编辑 (`svg.fonttype='none'`)
    <stem>.png          300 dpi, 白底 (报告内嵌用)
    <stem>.pdf          矢量, 出版投稿用 (保持调用方原有 dpi/参数)
    _hires_1200dpi/<stem>.png   1200 dpi, 白底 (印刷/审稿人索图用)

`_hires_1200dpi/` 是**独立子目录**而非同目录改名: 报告层 `collect_plots`
按 `(parent, stem)` 去重 (png>pdf>svg), 同目录多一张 PNG 会被当成第二张图
重复内嵌; 放进独立目录后报告 glob 不会命中, 版面不变。

两种用法
--------
1. **显式调用** (推荐新代码)::

       from mmpv_common.figure_export import export_figure
       export_figure(fig, "out/fig1.pdf", bbox_inches="tight")

2. **全局钩子** (零改动覆盖既有 `fig.savefig(...)` / `plt.savefig(...)`)::

       from mmpv_common.figure_export import install_savefig_hook
       install_savefig_hook()          # matplotlib 尚未导入也可调用 (延迟安装)

   安装后任何 `Figure.savefig("x.png")` 都会自动补齐其余三种格式。
   同一 stem 的兄弟文件**已存在则跳过**, 所以脚本里既写 png 又写 pdf 的
   老写法不会重复渲染 (总渲染次数 = 每种格式一次)。

环境变量开关
------------
    MMPV_FIG_EXPORT=0            完全关闭 (回退到原生 savefig)
    MMPV_FIG_FORMATS=png,pdf,... 只产出列出的格式; 默认 svg,png,pdf,hires
    MMPV_FIG_PNG_DPI=300         普通 PNG dpi
    MMPV_FIG_HIRES_DPI=1200      高清 PNG dpi
    MMPV_FIG_HIRES_DIR=_hires_1200dpi

安全边界
--------
* `PdfPages` / 文件对象 / BytesIO 目标**原样透传** (多页 PDF 不是单图, 不补格式)。
* 单个兄弟格式写失败只告警, 绝不中断管线 (PDF 已落盘, 不能因 SVG 失败拖垮 stage)。
* 预估像素数超过 `MAX_HIRES_PIXELS` 时自动降级到 600 dpi 并告警, 避免 OOM。
"""
from __future__ import annotations

import os
import sys
import threading

__all__ = [
    "ALL_FORMATS", "DEFAULT_PNG_DPI", "DEFAULT_HIRES_DPI", "HIRES_SUBDIR",
    "export_figure", "install_savefig_hook", "uninstall_savefig_hook",
    "hook_installed", "target_path", "figure_paths",
]

DEFAULT_PNG_DPI = 300
DEFAULT_HIRES_DPI = 1200
HIRES_SUBDIR = "_hires_1200dpi"
ALL_FORMATS = ("svg", "png", "pdf", "hires")

# 1200 dpi 下超过这个像素数就降级 (20000x20000 ≈ 4 亿像素, 约 1.6 GB RGBA)
MAX_HIRES_PIXELS = 400_000_000
HIRES_FALLBACK_DPI = 600

_LOCK = threading.RLock()
_STATE = {"installed": False, "finder": None, "orig": None, "suppress": 0}


# ═══════════════════════════════════════════════════════════════════
# 配置
# ═══════════════════════════════════════════════════════════════════

def _env_flag(name: str, default: str) -> str:
    v = os.environ.get(name)
    return default if v is None or v == "" else v


def _enabled() -> bool:
    return _env_flag("MMPV_FIG_EXPORT", "1").strip().lower() not in (
        "0", "false", "no", "off")


def _formats() -> tuple:
    raw = _env_flag("MMPV_FIG_FORMATS", ",".join(ALL_FORMATS))
    want = tuple(x.strip().lower() for x in raw.split(",") if x.strip())
    return tuple(f for f in want if f in ALL_FORMATS)


def _png_dpi() -> int:
    try:
        return max(1, int(_env_flag("MMPV_FIG_PNG_DPI", str(DEFAULT_PNG_DPI))))
    except ValueError:
        return DEFAULT_PNG_DPI


def _hires_dpi() -> int:
    try:
        return max(1, int(_env_flag("MMPV_FIG_HIRES_DPI", str(DEFAULT_HIRES_DPI))))
    except ValueError:
        return DEFAULT_HIRES_DPI


def _hires_dirname() -> str:
    return _env_flag("MMPV_FIG_HIRES_DIR", HIRES_SUBDIR)


# ═══════════════════════════════════════════════════════════════════
# 尺寸安全 (matplotlib 渲染器硬限: 每边 < 2^16 = 65536 px)
# 246 实测: 1234 样本/1150 病毒规模下 occurrence 柱 130514px、共丰度热图
# 75480px、filter 柱 143550px —— 超限即 ValueError 且此前无人兜底。
# ═══════════════════════════════════════════════════════════════════

RENDER_PX_CAP = 60000  # 留 5k 余量


def clamp_dpi(fig, dpi, cap_px: int = RENDER_PX_CAP) -> int:
    """按当前 figsize 返回不超过 cap_px 的安全 dpi。"""
    try:
        w, h = fig.get_size_inches()
        m = max(float(w), float(h))
    except Exception:
        return int(dpi)
    if m <= 0:
        return int(dpi)
    return max(10, min(int(dpi), int(cap_px // m)))


def safe_savefig(fig, path, dpi=300, logger=None, **kw) -> bool:
    """clamp dpi 的 savefig: 仍超限或渲染失败时降级为告警而不是崩掉整个绘图阶段。"""
    try:
        fig.savefig(path, dpi=clamp_dpi(fig, dpi), **kw)
        return True
    except ValueError as e:
        _msg = f"[figure_export] 跳过 {os.path.basename(str(path))}: {e}"
        if logger is not None:
            logger.warning(_msg)
        else:
            print(_msg, file=sys.stderr, flush=True)
        return False


# ═══════════════════════════════════════════════════════════════════
# 路径
# ═══════════════════════════════════════════════════════════════════

def _split(path) -> tuple:
    """返回 (目录, 主名); 目录为空时用 '.'。"""
    p = os.fspath(path)
    d, base = os.path.split(p)
    stem, _ext = os.path.splitext(base)
    if not stem:                      # 传入的是目录/空名 → 不猜
        raise ValueError(f"figure_export: 无法从 {p!r} 推断主名")
    return (d or "."), stem


def target_path(path, fmt: str) -> str:
    """给定任意一个目标路径, 返回 `fmt` 格式的落盘路径。

    fmt ∈ {"svg", "png", "pdf", "hires"}。
    """
    fmt = fmt.lower()
    if fmt not in ALL_FORMATS:
        raise ValueError(f"未知格式 {fmt!r}; 可选 {ALL_FORMATS}")
    d, stem = _split(path)
    if fmt == "hires":
        return os.path.join(d, _hires_dirname(), stem + ".png")
    return os.path.join(d, f"{stem}.{fmt}")


def figure_paths(path) -> dict:
    """返回四种格式的完整路径 (不渲染)。"""
    return {f: target_path(path, f) for f in ALL_FORMATS}


# ═══════════════════════════════════════════════════════════════════
# 渲染
# ═══════════════════════════════════════════════════════════════════

def _base_kwargs(savefig_kw: dict) -> dict:
    """剥离 format/dpi, 保留 bbox_inches / transparent / facecolor / pad_inches 等。"""
    kw = {k: v for k, v in savefig_kw.items() if k not in ("format", "dpi")}
    return kw


def _save_one(fig, path: str, fmt: str, savefig_kw: dict,
              png_dpi: int, hires_dpi: int, png_facecolor: str = "white") -> bool:
    """用**原生** savefig 写单个格式; 返回是否成功。"""
    from matplotlib.figure import Figure

    kw = _base_kwargs(savefig_kw)
    d = os.path.dirname(os.path.abspath(path))
    if d:
        os.makedirs(d, exist_ok=True)

    if fmt == "svg":
        kw["format"] = "svg"
        # SVG 是矢量; dpi 只影响栅格化元素 (imshow), 沿用调用方或 PNG 口径
        kw.setdefault("dpi", savefig_kw.get("dpi") or png_dpi)
    elif fmt == "pdf":
        kw["format"] = "pdf"
        if savefig_kw.get("dpi") is not None:
            kw["dpi"] = savefig_kw["dpi"]
    elif fmt == "png":
        kw["format"] = "png"
        kw["dpi"] = png_dpi
        kw.setdefault("facecolor", png_facecolor)
    elif fmt == "hires":
        kw["format"] = "png"
        kw["dpi"] = _effective_hires_dpi(fig, hires_dpi)
        kw.setdefault("facecolor", "white")
    else:
        raise ValueError(fmt)

    orig = _STATE["orig"] or Figure.savefig
    _STATE["suppress"] += 1
    try:
        orig(fig, path, **kw)
    finally:
        _STATE["suppress"] -= 1
    return True


def _effective_hires_dpi(fig, hires_dpi: int) -> int:
    """按 figure 物理尺寸估算像素数, 过大则降级并告警 (防 OOM)。"""
    try:
        w_in, h_in = fig.get_size_inches()
        px = float(w_in) * hires_dpi * float(h_in) * hires_dpi
    except Exception:                                  # noqa: BLE001
        return hires_dpi
    if px <= MAX_HIRES_PIXELS:
        return hires_dpi
    scale = (MAX_HIRES_PIXELS / px) ** 0.5
    fallback = max(1, int(hires_dpi * scale))
    # [尺寸安全] 面积达标还不够: 单边仍受 2^16 硬限 (246 实测 600in 宽在
    # 面积达标 dpi 下单边仍 99600px 崩溃)
    side = max(float(w_in), float(h_in))
    if side > 0:
        side_cap = max(1, int(65535 / side))
        if side_cap < fallback:
            fallback = side_cap
    fallback = min(fallback, HIRES_FALLBACK_DPI)
    print(f"[figure_export] 警告: {w_in:.1f}x{h_in:.1f} in @ {hires_dpi} dpi "
          f"≈ {px/1e6:.0f} Mpx 超限, 高清 PNG 降级到 {fallback} dpi",
          file=sys.stderr)
    return fallback


def export_figure(fig, path, *, formats=None, png_dpi=None, hires_dpi=None,
                  overwrite=True, path_overrides=None, png_facecolor=None,
                  **savefig_kw) -> dict:
    """把同一 figure 写成四种格式。

    Parameters
    ----------
    fig : matplotlib.figure.Figure
    path : str | os.PathLike
        任意一个目标路径; 只取**目录 + 主名**, 扩展名由 formats 决定。
    formats : iterable, optional
        默认 `ALL_FORMATS` (svg/png/pdf/hires), 可用 `MMPV_FIG_FORMATS` 覆盖。
    png_dpi, hires_dpi : int, optional
        默认 300 / 1200, 可用环境变量覆盖。
    overwrite : bool
        False 时已存在的文件跳过。
    path_overrides : dict, optional
        指定某格式的自定义落盘路径, 如 ``{"png": "out/自定义.png"}``。
    png_facecolor : str, optional
        普通 PNG 的底色 (默认 white); 高清 PNG 恒用 white。
    **savefig_kw
        透传给 savefig (bbox_inches / transparent / pad_inches / facecolor …)。
        `dpi` 只作用于 PDF (保持历史输出), PNG 用 png_dpi, 高清用 hires_dpi。

    Returns
    -------
    dict — {"svg": path|None, "png": ..., "pdf": ..., "hires": ...};
        写失败的格式值为 None (只告警, 不抛错)。

    Notes
    -----
    `MMPV_FIG_EXPORT=0` 只关闭**自动补齐** (全局 savefig 钩子), 不关闭显式
    调用 —— 显式点名要写四格式的代码, 调用方意图明确, 不该被环境开关静默吞掉。
    需要完全关闭时改用 `formats=()`。
    """
    formats = tuple(formats) if formats is not None else _formats()
    png_dpi = _png_dpi() if png_dpi is None else int(png_dpi)
    hires_dpi = _hires_dpi() if hires_dpi is None else int(hires_dpi)
    overrides = {str(k).lower(): v for k, v in (path_overrides or {}).items()}
    png_fc = png_facecolor or "white"

    out = {f: None for f in ALL_FORMATS}
    for fmt in ALL_FORMATS:
        if fmt not in formats:
            continue
        dst = str(overrides[fmt]) if fmt in overrides else target_path(path, fmt)
        if not overwrite and os.path.exists(dst):
            out[fmt] = dst
            continue
        try:
            _save_one(fig, dst, fmt, savefig_kw, png_dpi, hires_dpi,
                      png_facecolor=png_fc)
            out[fmt] = dst
        except Exception as e:                          # noqa: BLE001
            print(f"[figure_export] 警告: {fmt} 写出失败 ({dst}): {e}",
                  file=sys.stderr)
    return out


# ═══════════════════════════════════════════════════════════════════
# 全局 savefig 钩子
# ═══════════════════════════════════════════════════════════════════

def _sibling_already_written(path) -> bool:
    """本次 savefig 的目标是否是这组兄弟文件之一 (避免重复渲染)。"""
    try:
        return os.path.exists(os.fspath(path))
    except Exception:                                   # noqa: BLE001
        return False


def _fill_siblings(fig, fname, savefig_kw: dict) -> None:
    """给一次原生 savefig 补齐其余格式 (已存在则跳过)。"""
    if not _enabled():
        return
    fmt_of_target = os.path.splitext(os.fspath(fname))[1].lower().lstrip(".")
    formats = _formats()
    png_dpi = _png_dpi()
    hires_dpi = _hires_dpi()

    for fmt in formats:
        dst = target_path(fname, fmt)
        # 本次调用刚写过的那个格式跳过 (它已由原生 savefig 落盘)
        if fmt == fmt_of_target and fmt != "hires":
            continue
        if _sibling_already_written(dst):
            continue
        try:
            _save_one(fig, dst, fmt, savefig_kw, png_dpi, hires_dpi)
        except Exception as e:                          # noqa: BLE001
            print(f"[figure_export] 警告: 补 {fmt} 失败 ({dst}): {e}",
                  file=sys.stderr)


def _make_patched(orig):
    def _patched_savefig(self, fname, *args, **kwargs):
        # ① 非路径目标 (PdfPages / 文件对象 / BytesIO) → 原样透传, 不补格式
        if _STATE["suppress"] or not isinstance(fname, (str, os.PathLike)):
            return orig(self, fname, *args, **kwargs)

        # [尺寸安全] 栅格目标在渲染前按最长边钳制 dpi (matplotlib 硬限 2^16/边);
        # 矢量 (pdf/svg) 无此限制, 不动 dpi
        _sfx = str(fname).lower().rsplit(".", 1)[-1]
        if _sfx in ("png", "jpg", "jpeg", "tif", "tiff", "webp"):
            kwargs["dpi"] = clamp_dpi(self, kwargs.get("dpi") or self.dpi)

        result = orig(self, fname, *args, **kwargs)

        # ② 补齐兄弟格式; 任何异常都不得影响主调用
        try:
            _fill_siblings(self, fname, kwargs)
        except Exception as e:                          # noqa: BLE001
            print(f"[figure_export] 警告: 补齐兄弟格式失败: {e}", file=sys.stderr)
        return result

    _patched_savefig.__name__ = "savefig"
    _patched_savefig.__doc__ = orig.__doc__
    _patched_savefig.__wrapped__ = orig
    return _patched_savefig


def _patch_now() -> bool:
    """matplotlib.figure 已加载时直接替换 Figure.savefig。"""
    from matplotlib.figure import Figure
    with _LOCK:
        if _STATE["installed"]:
            return True
        orig = getattr(Figure.savefig, "__wrapped__", None)
        _STATE["orig"] = orig or Figure.savefig
        Figure.savefig = _make_patched(_STATE["orig"])
        _STATE["installed"] = True
    return True


class _DeferredFinder:
    """matplotlib 尚未导入时, 在 `matplotlib.figure` 执行完后安装钩子。"""

    def __init__(self):
        self._done = False
        self._resolving = False

    def find_spec(self, fullname, path=None, target=None):
        if fullname != "matplotlib.figure" or self._done or self._resolving:
            return None
        import importlib.abc
        import importlib.util

        # importlib.util.find_spec 会重走 meta_path → 用 _resolving 挡住自递归。
        self._resolving = True
        try:
            spec = importlib.util.find_spec(fullname)
        finally:
            self._resolving = False
        if spec is None or spec.loader is None:
            self._done = True
            return None
        inner = spec.loader

        class _Loader(importlib.abc.Loader):
            def create_module(self, s):
                return inner.create_module(s)

            def exec_module(self, module):
                inner.exec_module(module)
                try:
                    _patch_now()
                except Exception as e:                  # noqa: BLE001
                    print(f"[figure_export] 警告: 钩子安装失败: {e}",
                          file=sys.stderr)

        spec.loader = _Loader()
        self._done = True
        return spec

    def invalidate_caches(self):
        pass


def install_savefig_hook() -> bool:
    """安装全局 savefig 钩子 (幂等; matplotlib 未导入时延迟安装)。"""
    if not _enabled():
        return False
    with _LOCK:
        if _STATE["installed"]:
            return True
    if "matplotlib.figure" in sys.modules:
        try:
            return _patch_now()
        except Exception as e:                          # noqa: BLE001
            print(f"[figure_export] 警告: 钩子安装失败: {e}", file=sys.stderr)
            return False
    with _LOCK:
        if _STATE["finder"] is None:
            finder = _DeferredFinder()
            sys.meta_path.insert(0, finder)
            _STATE["finder"] = finder
    return True


def uninstall_savefig_hook() -> bool:
    """撤销钩子 (测试/调试用)。"""
    with _LOCK:
        finder = _STATE.get("finder")
        if finder is not None and finder in sys.meta_path:
            sys.meta_path.remove(finder)
        _STATE["finder"] = None
        if not _STATE["installed"]:
            return False
        try:
            from matplotlib.figure import Figure
            Figure.savefig = _STATE["orig"] or Figure.savefig
        except Exception:                               # noqa: BLE001
            pass
        _STATE["installed"] = False
        _STATE["orig"] = None
    return True


def hook_installed() -> bool:
    return bool(_STATE["installed"])
