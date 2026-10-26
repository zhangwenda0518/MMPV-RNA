#!/usr/bin/env python3
"""check_fig_layout.py — 通用版式体检：文字重叠 + 文字出画布 (2026-10-14)

由一次性审查脚本（`reports/_archive/review_dirs/_gcva_time_review_20261013/`）
正式并入管线，并把入参泛化（原版写死 `--datadir/--outdir`，只能体检那两个脚本）。

判据（与项目既有 `_check_*_overlap.py` 同源）
============================================
  1. **0 对文字重叠**（TOL = 0.5 px）
  2. **无 TEXT-OFF-CANVAS**（文字 bbox 越出画布边界）

口径要点（都是踩过的坑）
========================
  · 必须用 `get_window_extent(renderer)` 量**实际渲染** bbox —— `tight_layout()`
    只调 axes、**不搬 `fig.text`**，用它的推算值会漏掉主标题溢出。
  · 只对**可见且有内容**的 Text 计入（跳过空串与 `set_visible(False)`）。
  · ⚠⚠ **对数轴有假阳性**：locator 会生成**视域外**的刻度 Text（如
    xlim=[22, 40000] 仍产出 10^1 / 10^5），这些 label `get_visible()` 为 True
    但绘制时被裁 → 直接量 bbox 会误报 OFF-CANVAS（实测 5 个假阳性）。
    权威判据 = 刻度**位置**是否落在 view interval 内（`_tick_in_range`）。
  · ⚠⚠ 本脚本会**重渲**目标脚本 —— 它不是只读工具，跑完以**体检后产物**为准。
  · ⚠ 不适用：图片嵌入型面板 / ComplexUpset / 桑基（`axis("off")`）——
    这类图的"文字"不是 matplotlib Text，须目视。本脚本对它们会**假通过**。

用法
====
  # 直接体检某个绘图脚本（透传它自己的参数）
  python -m utils.check_fig_layout --script path/to/fig.py --args "--outdir /tmp/x"

  # 体检管线模块
  python -m utils.check_fig_layout --module utils.ess_convergence_fig \
      --args "--chains a.log,b.log --outdir /tmp/x"

  # 也可当库用
  from utils.check_fig_layout import check_figure
  off, ov = check_figure(fig)
"""
from __future__ import annotations

import argparse
import runpy
import shlex
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.figure as mfig
import matplotlib.pyplot as plt

TOL = 0.5          # px


def _tick_in_range(ax, axis_name):
    """返回该轴上**落在视域内**的刻度位置集合（避开对数轴假阳性）。"""
    axis = ax.xaxis if axis_name == "x" else ax.yaxis
    lo, hi = sorted(axis.get_view_interval())
    return [loc for loc, _ in zip(axis.get_majorticklocs(),
                                  axis.get_majorticklabels())
            if lo - 1e-9 <= loc <= hi + 1e-9]


def _iter_texts(fig):
    for ax in fig.get_axes():
        if not ax.get_visible():
            continue
        for t in ax.texts + [ax.title, ax.xaxis.label, ax.yaxis.label]:
            yield t
        xin = set(_tick_in_range(ax, "x"))
        yin = set(_tick_in_range(ax, "y"))
        for loc, t in zip(ax.xaxis.get_majorticklocs(), ax.get_xticklabels()):
            if loc in xin:
                yield t
        for loc, t in zip(ax.yaxis.get_majorticklocs(), ax.get_yticklabels()):
            if loc in yin:
                yield t
        lg = ax.get_legend()
        if lg is not None and lg.get_visible():
            for t in lg.get_texts():
                yield t
    for t in fig.texts:
        yield t


def _bbox(t, renderer):
    try:
        bb = t.get_window_extent(renderer=renderer)
    except Exception:
        return None
    if bb.width <= 0 or bb.height <= 0:
        return None
    return bb


def _overlap(a, b):
    dx = min(a.x1, b.x1) - max(a.x0, b.x0)
    dy = min(a.y1, b.y1) - max(a.y0, b.y0)
    if dx <= TOL or dy <= TOL:
        return None
    return dx * dy


def check_figure(fig, tag="fig", verbose=True):
    """体检单个 Figure，返回 (n_off_canvas, n_overlap)。"""
    # ⚠ runpy 场景下 canvas 可能是 FigureCanvasBase（无 get_renderer）
    # → 显式挂 Agg canvas 再 draw。
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    if not hasattr(fig.canvas, "get_renderer"):
        fig.set_canvas(FigureCanvasAgg(fig))
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    W, H = fig.canvas.get_width_height()

    items = []
    for t in _iter_texts(fig):
        if not t.get_visible():
            continue
        s = t.get_text()
        if not s or not s.strip():
            continue
        bb = _bbox(t, r)
        if bb is None:
            continue
        items.append((bb, s))

    off = [(s, (bb.x0, bb.y0, bb.x1, bb.y1)) for bb, s in items
           if bb.x0 < -TOL or bb.y0 < -TOL or bb.x1 > W + TOL or bb.y1 > H + TOL]

    ov = []
    for i in range(len(items)):
        for j in range(i + 1, len(items)):
            area = _overlap(items[i][0], items[j][0])
            if area and area > TOL:
                ov.append((items[i][1], items[j][1], area))

    if verbose:
        print("[%s] 画布 %dx%d px | 计入文字 %d 个" % (tag, W, H, len(items)))
        print("   TEXT-OFF-CANVAS : %d  %s" % (len(off), "OK" if not off else "✗"))
        for s, bb in off[:12]:
            print("     ✗ %-42s bbox=(%.1f, %.1f, %.1f, %.1f)"
                  % (s[:42].replace("\n", "|"), bb[0], bb[1], bb[2], bb[3]))
        print("   重叠对数         : %d  %s" % (len(ov), "OK" if not ov else "✗"))
        for a, b, area in ov[:15]:
            print("     ✗ %-34s × %-34s  area=%.1f px²"
                  % (a[:34].replace("\n", "|"), b[:34].replace("\n", "|"), area))
    return len(off), len(ov)


def _run_and_capture(argv, module=None, script=None):
    """跑目标脚本并截获它 savefig 过的所有 Figure。

    ⚠ 为什么用 savefig 钩子而不是 `plt.gcf()`：绘图脚本普遍在保存后
    `plt.close(fig)`，事后拿不到对象；钩子在 savefig 之前把引用留住。
    """
    captured = []
    orig = mfig.Figure.savefig

    def _hook(self, *a, **k):
        # ⚠ 去重：同一 Figure 存 PDF 再存 PNG 会触发两次钩子，
        # 不去重就会把同一张图体检两遍（输出里出现重复的 fig1/fig2）。
        if not any(c is self for c in captured):
            captured.append(self)
        return orig(self, *a, **k)

    mfig.Figure.savefig = _hook
    old_argv = sys.argv
    sys.argv = argv
    try:
        if module:
            runpy.run_module(module, run_name="__main__", alter_sys=True)
        else:
            runpy.run_path(script, run_name="__main__")
    except SystemExit:
        # ⚠ 目标脚本普遍以 `sys.exit(main())` 结尾 —— SystemExit 会一路穿过
        # runpy 冒到本模块，把体检结果一起带走。这里必须吞掉。
        pass
    finally:
        mfig.Figure.savefig = orig
        sys.argv = old_argv
    return captured


def main(argv=None):
    ap = argparse.ArgumentParser(description="通用版式体检（0 重叠 + 无 OFF-CANVAS）")
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--script", help="绘图脚本路径")
    src.add_argument("--module", help="以 -m 方式运行的模块名，如 utils.ess_convergence_fig")
    ap.add_argument("--args", default="",
                    help="透传给目标脚本的参数串（会被 shlex 切分）")
    args = ap.parse_args(argv)

    extra = shlex.split(args.args)
    if args.module:
        argv_run = [args.module] + extra
        captured = _run_and_capture(argv_run, module=args.module)
    else:
        argv_run = [args.script] + extra
        captured = _run_and_capture(argv_run, script=args.script)

    if not captured:
        print("⚠ 没有截获到任何 Figure —— 目标脚本可能没有 savefig，"
              "或参数不对（检查 --args）")
        return 2

    tot_off = tot_ov = 0
    for i, f in enumerate(captured):
        o, v = check_figure(f, "fig%d" % (i + 1))
        tot_off += o
        tot_ov += v
    print()
    print("=" * 62)
    ok = (tot_off == 0 and tot_ov == 0)
    print("总判定: %s   (OFF=%d, 重叠=%d)"
          % ("✓ 通过" if ok else "✗ 不通过", tot_off, tot_ov))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
