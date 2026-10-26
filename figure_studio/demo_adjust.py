# -*- coding: utf-8 -*-
"""Figure Studio · 调色算子对照图（同布局，改的是 adjust）。

版式与 ``demo_types.py`` **完全一致**（同样 8 块、同样的行/列/尺寸），
只有每块的 ``adjust`` 不同 —— 方便与 ``demo_types.png`` 直接对照。

每块的算子都按**该面板的真实问题**选，不是随机套：
* A 环形树里 ``#0172b2`` 比项目 ``#0072B2`` 差 1 个色阶（抗锯齿偏色）→ 容差替换拉回
* C 桑基节点是**α 混合**出来的浅色（``#add1e6``），精确替换抓不住 → 只能整体着色
* E 是 R 出的调色板 PNG、带灰底 → 近白替换成纯白 + 自动对比度
* G/H 各有一个不在项目色板里的主色 → 换成项目语义色

运行::

    C:/Python312/python.exe demo_adjust.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import composer  # noqa: E402
import demo_types as D  # noqa: E402  （复用面板清单与行结构）

# letter -> (adjust, 题注)
ADJUSTS: dict[str, tuple[dict, str]] = {
    "A": ({"remap": [{"from": "#0172b2", "to": "#0072B2", "tol": 8},
                     {"from": "#e6a617", "to": "#E69F00", "tol": 10}]},
          "colour remap: near-miss tints snapped back to the project palette"),
    "B": ({"grayscale": True},
          "grayscale (monochrome print)"),
    "C": ({"tint": {"color": "#0072B2", "amount": 0.30}},
          "single-colour tint (alpha-blended nodes cannot be remapped)"),
    "D": ({"hue": 150.0},
          "hue rotation +150 deg"),
    "E": ({"flatten": "#ffffff", "flatten_tol": 18, "autocontrast": True},
          "near-white to pure white + autocontrast"),
    "F": ({"brightness": 0.85, "contrast": 1.45},
          "brightness 0.85 / contrast 1.45 (deepen the heat map)"),
    "G": ({"remap": [{"from": "#cc0000", "to": "#D55E00", "tol": 12}]},
          "colour remap: off-palette red to project vermillion"),
    "H": ({"remap": [{"from": "#59aeb6", "to": "#56B4E9", "tol": 14}],
           "saturation": 1.15},
          "colour remap: off-palette teal to project sky blue"),
}


def main() -> int:
    layout = D.build()                       # 版式完全复用 demo_types
    for it in layout["items"]:
        L = it["letter"]
        adj, cap = ADJUSTS[L]
        it["adjust"] = adj
        it["caption"] = f"{L} · {cap}"
    layout["title"]["text"] = "Figure Studio — colour adjustments (same layout as demo_types)"
    layout["title"]["subtitle"] = ("image-level ops on real panels: colour remap / grayscale / "
                                   "tint / hue / flatten / brightness-contrast")

    res = composer.render_layout(layout, str(HERE / "outputs"), "demo_adjust",
                                formats=("png", "pdf"))
    (HERE / "layouts").mkdir(exist_ok=True)
    (HERE / "layouts" / "demo_adjust.json").write_text(
        json.dumps(layout, ensure_ascii=False, indent=2), encoding="utf-8")

    for it in layout["items"]:
        print(f"  {it['letter']}  {it['w_in']:6.3f} x {it['h_in']:6.3f} in   "
              f"adjust={json.dumps(it['adjust'], ensure_ascii=False)}")
    print("告警:", res["warnings"] or "无")
    print("PNG:", res["png"], res["size_px"])
    print("PDF:", res["pdf"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
