#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
check_ef_library.py
===================
体检 `make_sci_composite_figs.py::EF_LIBRARY` 里每一张候选图：
裁白边 + 切掉源图标题后，**在组合图里占多高、印刷后字号还剩多少**。

为什么需要它（两条实测约束，详见 reports/virome_analysis_pipeline/figures/README_SCI_composite.md §七）：

* 源图是整页画布（3000 px = 300 dpi 下 10 in；4800 px = 16 in），缩到组合图的
  通栏（9.775 in）会按比例缩字号。印刷后的字号是

      font_print = 源字号 × (通栏宽 / 源画布宽) × (期刊印刷宽 / 组合图设计宽)
                 = 12 pt × (D/10) × (7.09/D) = 12 × 0.709 ≈ 8.5 pt   （源画布 10 in）

  **设计宽 D 被约掉了** —— 这是"比性质"，把组合图整体放大并不会让嵌图字号变大。
* 于是：通栏 ≈ 7.4-7.9 pt（可进正文）、半栏 ≈ 3.5 pt（不行）。**E/F 只能各占一个通栏**，
  总高因此涨到 22 in 量级 —— 这个代价必须提前算出来。

用法（在服务器上）：
    export LD_LIBRARY_PATH=/home/zhangwenda/mambaforge/lib:$LD_LIBRARY_PATH
    /home/zhangwenda/mambaforge/bin/python reports/virome_analysis_pipeline/check_ef_library.py
"""
from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import make_sci_composite_figs as M          # noqa: E402

# 组合图几何（与 make_sci_composite_figs.py 的 main() 保持一致）
DESIGN_W = 11.5                              # 设计宽 (in)
PANEL_L, PANEL_R = 0.095, 0.945              # B~F 区左右边界 (figure 分数)
EF_W = (PANEL_R - PANEL_L) * DESIGN_W        # 通栏可用宽度 = 9.775 in
PRINT_W = 7.09                               # 期刊双栏 180 mm
SRC_PT = 12.0                                # 源图正文基准字号 (pt)
SRC_DPI = 300.0                              # 源图出图分辨率

DEFAULT_POST = ("/home/zhangwenda/MMPV-paper/goji-virome/03_known_virus/"
                "known_virus_barbarum/04_post_analysis/"
                "Cytorhabdovirus_sp._lycii_OR489165.1")


def main(argv=None):
    post = (argv[0] if argv else None) or os.environ.get("POST_DIR") or DEFAULT_POST
    if not os.path.isdir(post):
        print(f"[err] 目录不存在: {post}", file=sys.stderr)
        return 2
    print(f"[post] {post}")
    print(f"[geom] 通栏 = {EF_W:.3f} in  设计宽 {DESIGN_W} in  印刷宽 {PRINT_W} in\n")
    print(f"{'key':16s} {'裁后 px':>12s} {'宽高比':>7s} {'通栏高':>8s} "
          f"{'印刷字号':>9s}  判定")
    print("-" * 62)
    bad = 0
    for k in M.EF_IMG_KEYS:
        try:
            p, _cap = M.resolve_ef("E", k, post)
            im = M._crop_white(p)
        except Exception as exc:                       # 缺文件 / 全白
            print(f"{k:16s}   !! {exc}")
            bad += 1
            continue
        ar = im.width / im.height
        h = EF_W / ar
        # 源图按 SRC_DPI 出 → 裁后像素宽对应的"英寸"；再乘印前缩放
        pt = SRC_PT * (EF_W / (im.width / SRC_DPI)) * (PRINT_W / DESIGN_W)
        flag = "OK" if pt >= 6.0 else "** 偏小 **"
        if pt < 6.0:
            bad += 1
        print(f"{k:16s} {im.width:5d}x{im.height:<6d} {ar:7.2f} {h:7.2f}in "
              f"{pt:8.1f}pt  {flag}")
    print("-" * 62)
    print(f"[sum] {len(M.EF_IMG_KEYS)} 张候选，其中 {bad} 张印刷字号 < 6 pt（不建议进正文）")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
