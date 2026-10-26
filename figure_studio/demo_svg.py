# -*- coding: utf-8 -*-
"""SVG 面板实测：矢量级调色 / 图例隐藏 / 字号缩放 + 与 PNG 混排。

与 demo_types.py / demo_adjust.py 的区别：这里用 **SVG 源图**，
编辑发生在 XML 层（svg_edit.apply_svg_edit）而非像素层。
产物 demo_svg_mixed.png 同时含 3 个 SVG 面板 + 1 个 PNG 面板。
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import composer as C  # noqa: E402

OUT = HERE / "outputs"
LAY = HERE / "layouts"

# ── 面板来源（真实成品图）────────────────────────────────────
PANELS = {
    "A": HERE / "outputs/_svgdemo/panel_net.svg",
    "B": HERE / "outputs/_svgdemo/panel_tracks.svg",
    "C": HERE / "outputs/_svgdemo/panel_pca.svg",
    "D": Path("D:/桌面/MMPV主图档案库_20261014/13_HostPrediction_宿主预测/"
              "成品图/Fig_GojiVirome_HostPrediction_APanel.png"),
}

CAPS = {
    "A": "SVG: recoloured to the project palette, text enlarged 1.35x",
    "B": "SVG: both legends hidden, text enlarged 1.2x",
    "C": "SVG (R output): recoloured, text enlarged 1.25x",
    "D": "PNG: grayscale image-level adjust after a vector step",
}

EDITS = {
    "A": {"svg": {"colors": [{"from": "#0f4d92", "to": "#0072B2"},
                             {"from": "#d9544d", "to": "#D55E00"}],
                  "font_scale": 1.35}},
    "B": {"svg": {"legend": {"action": "hide"}, "font_scale": 1.2}},
    "C": {"svg": {"colors": [{"from": "#1b9e77", "to": "#0072B2"}],
                  "font_scale": 1.25}},
    "D": {"svg": {"font_scale": 1.15}, "grayscale": True},
}

ROWS = [["A", "B"], ["C", "D"]]
MARGIN, GAP, TOP, CAP_ALLOW, TITLE_ALLOW = 0.10, 0.10, 0.62, 0.30, 0.10
W = 10.0


def build() -> dict:
    items, y = [], TOP + TITLE_ALLOW
    for row in ROWS:
        aspects = {}
        for L in row:
            pr = C.probe(str(PANELS[L]), dpi=300, adjust=EDITS[L])
            aspects[L] = pr["aspect"]
        h = (W - 2 * MARGIN - (len(row) - 1) * GAP) / sum(
            1.0 / aspects[L] for L in row)
        x = MARGIN
        for L in row:
            w = h / aspects[L]
            items.append({"id": "i_" + L, "src": str(PANELS[L]),
                          "x_in": round(x, 4), "y_in": round(y, 4),
                          "w_in": round(w, 4), "h_in": round(h, 4),
                          "letter": L, "caption": CAPS[L], "z": 0,
                          "adjust": EDITS[L]})
            x += w + GAP
        y += h + GAP + CAP_ALLOW
    # ⚠ 不能写成 y - GAP - CAP_ALLOW + MARGIN：那样会把末行题注空间抵消掉，
    #   底部只剩 MARGIN - gap = 0.04 in，题注必被裁（且零报错）。
    #   demo_types.py 用的就是这个正确式。
    H = y - GAP + MARGIN
    return {
        "canvas": {"width_in": W, "height_in": round(H, 4), "dpi": 300},
        "style": {"bg": "white", "frame": True, "frame_color": "#d9d9d9",
                  "frame_lw": 0.6, "letter_fs": 11, "letter_dx": 0.0,
                  "letter_dy": 1.02, "show_letters": True, "caption_fs": 7,
                  "caption_gap_in": 0.06, "caption_linespacing": 1.2,
                  "margin_in": MARGIN, "pack_margin_in": MARGIN,
                  "gap_in": GAP},
        "title": {"text": "SVG panels: vector-level editing",
                  "subtitle": "colour / legend / font-size edited in the XML, "
                              "then composited with a raster PNG",
                  "fs": 13, "sub_fs": 8.5, "x_in": MARGIN, "y_in": 0.30},
        "items": items,
    }


def main() -> int:
    OUT.mkdir(exist_ok=True)
    LAY.mkdir(exist_ok=True)
    layout = build()
    r = C.render_layout(layout, str(OUT), "demo_svg_mixed")
    worst = 0.0
    for it in layout["items"]:
        pr = C.probe(it["src"], dpi=300, adjust=it.get("adjust"))
        worst = max(worst, abs(it["h_in"] / it["w_in"] - pr["aspect"]) / pr["aspect"])
    print("[demo_svg] canvas %.3f x %.3f in, png %s"
          % (layout["canvas"]["width_in"], layout["canvas"]["height_in"],
             r["size_px"]))
    print("[demo_svg] max aspect deviation %.2e" % worst)
    for w in r["warnings"]:
        print("  -", w)
    (LAY / "demo_svg_mixed.json").write_text(
        json.dumps(layout, ensure_ascii=False, indent=1), encoding="utf-8")
    print("[demo_svg] layout -> layouts/demo_svg_mixed.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
