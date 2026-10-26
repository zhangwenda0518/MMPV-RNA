# -*- coding: utf-8 -*-
"""Figure Studio · 异构面板实测（不是条形图的那些）。

用本项目档案库里的**真实成品图**（环形 ML 树 / NeighborNet 网络 / 桑基 /
UpSet / 宿主预测调色板 PNG / SDT 热图 / SDT 分布 / 覆盖度折线）拼一张通栏图，
验证合成器对图型完全无感 —— 它只认像素和宽高比。

同时把算好的 layout 落到 ``layouts/demo_types.json``，
可以直接在网页端「加载布局」里打开继续拖。

运行::

    C:/Python312/python.exe demo_types.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import composer  # noqa: E402

ARCHIVE = Path(r"D:/桌面/MMPV主图档案库_20261014")
DISCO = Path(r"D:/桌面/延伸基因组/MMPV-RNA/reports/_archive/review_dirs"
             r"/_discovery_sci_20261012/_lt/10_Reports/S10_SCI_Discovery")

# letter -> (源图, 题注)
PANELS = {
    "A": (ARCHIVE / "03_MLtree_circular/成品图/7_MLtree_circular_host.png",
          "Maximum-likelihood phylogeny of GCVA (circular, tips coloured by host)"),
    "B": (ARCHIVE / "04_NeighborNet/成品图/8_NeighborNet_host.png",
          "NeighborNet split network of GCVA genomes"),
    "C": (ARCHIVE / "07_Fig_GojiVirome_NoveltySankey/成品图/4_Fig_GojiVirome_NoveltySankey.png",
          "Flow of known / novel virus detections (Sankey)"),
    "D": (ARCHIVE / "09_UpSet_工具分类宿主/成品图/Fig_GojiVirome_UpsetClassification.png",
          "Tool x host classification intersections (UpSet)"),
    "E": (ARCHIVE / "13_HostPrediction_宿主预测/成品图/Fig_GojiVirome_HostPrediction_APanel.png",
          "Host prediction, panel A (palette-mode PNG from R)"),
    "F": (DISCO / "Fig_GojiVirome_Potexvirus_SDT.png",
          "Sequence demarcation tool heat map (Potexvirus)"),
    "G": (DISCO / "Fig_GojiVirome_Potexvirus_SDTdist.png",
          "Pairwise SDT distance distribution"),
    "H": (DISCO / "Fig_GojiVirome_Potexvirus_Coverage.png",
          "Read coverage profile along the reference"),
}

ROWS = [["A", "B", "F"], ["C", "D", "E"], ["G", "H"]]

CANVAS_W = 9.32          # 通栏（项目惯例）
DPI = 300
MARGIN = 0.10
GAP = 0.12
CAPTION_ALLOW = 0.40     # 行间给题注留的高度（题注会按面板宽折行，按 3 行算）
TITLE_ALLOW = 0.86       # 顶部标题块 + 面板字母的高度（字母画在面板上缘之外）


def build() -> dict:
    info = {}
    for L, (src, _cap) in PANELS.items():
        if not src.is_file():
            raise SystemExit(f"面板 {L} 源图缺失: {src}")
        info[L] = composer.probe(str(src))

    avail = CANVAS_W - 2 * MARGIN
    items, y = [], MARGIN + TITLE_ALLOW
    for row in ROWS:
        # 行高 = (可用宽 - 行内 gap 总和) / Σ(1/aspect)  → 行内等高、整行刚好铺满
        # ⚠ 必须先扣 (n-1) 个 gap，否则整行会超宽 (n-1)*gap
        row_avail = avail - (len(row) - 1) * GAP
        h = row_avail / sum(1.0 / info[L]["aspect"] for L in row)
        x = MARGIN
        for L in row:
            w = h / info[L]["aspect"]
            items.append({
                "id": f"d_{L}", "src": str(PANELS[L][0]),
                "x_in": round(x, 4), "y_in": round(y, 4),
                "w_in": round(w, 4), "h_in": round(h, 4),
                "letter": L, "caption": PANELS[L][1], "z": 0,
            })
            x += w + GAP
        y += h + GAP + CAPTION_ALLOW
    # ⚠ 末行也要留 CAPTION_ALLOW：H 里扣掉它 → 最后一行题注被下缘裁掉（零报错）
    H = y - GAP + MARGIN

    return {
        "canvas": {"width_in": CANVAS_W, "height_in": round(H, 4), "dpi": DPI},
        "style": {"bg": "white", "frame": False, "frame_color": "#bbbbbb",
                  "frame_lw": 0.6, "letter_fs": 12, "letter_dx": 0.0,
                  "letter_dy": 1.01, "caption_fs": 6.6, "caption_gap_in": 0.07,
                  "caption_linespacing": 1.15, "show_letters": True},
        "title": {"text": "Figure Studio — heterogeneous panels",
                  "subtitle": "circular tree / network / Sankey / UpSet / heat map / line, "
                              "each placed at its native aspect ratio",
                  "fs": 13, "sub_fs": 8.5, "x_in": MARGIN, "y_in": 0.26},
        "items": items,
    }


def main() -> int:
    layout = build()
    out_dir = HERE / "outputs"
    res = composer.render_layout(layout, str(out_dir), "demo_types",
                                formats=("png", "pdf"))

    (HERE / "layouts").mkdir(exist_ok=True)
    (HERE / "layouts" / "demo_types.json").write_text(
        json.dumps(layout, ensure_ascii=False, indent=2), encoding="utf-8")

    print(f"画布  {layout['canvas']['width_in']} x {layout['canvas']['height_in']} in")
    for it in layout["items"]:
        print(f"  {it['letter']}  {it['w_in']:6.3f} x {it['h_in']:6.3f} in"
              f"  @ ({it['x_in']:5.2f}, {it['y_in']:5.2f})")
    print("形变告警:", res["warnings"] or "无")
    print("PNG:", res["png"], res["size_px"])
    print("PDF:", res["pdf"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
