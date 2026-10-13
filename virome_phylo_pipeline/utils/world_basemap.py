#!/usr/bin/env python3
"""
world_basemap.py — 离线世界海岸线/国界底图 (数据借自 EasyHap)
==============================================================

2026-10-12 新增。数据源: biosoft/EasyHap/easyhap/world_map.py 的
`get_world_lines()` (base85+zlib 内嵌, 零网络/零 geo 依赖):
  {"coastlines": [147 条折线], "countries": [450 条折线]}, 每条折线为
  [[lon, lat], ...] 的列表的列表 (多段)。

借鉴出处 (按项目规范注明): Guangqi-He/EasyHap (GPL-3.0),
easyhap/world_map.py; 本模块只做 matplotlib 绘制适配, 不复制其绘制代码。

用法:
    from utils.world_basemap import draw_world
    draw_world(ax)                      # 只画海岸线
    draw_world(ax, borders=True)        # 海岸线 + 国界
EasyHap 不可导入时 draw_world 静默返回 False (调用方自行降级),
绝不联网 —— 与 transmission_map 的「离线合规」红线一致。
"""

import os
import sys
from typing import Optional

_EASYHAP_PATH = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "..", "biosoft", "EasyHap")


def _world_lines():
    """加载 EasyHap 内嵌的世界折线; 失败返回 None。"""
    try:
        p = os.path.normpath(_EASYHAP_PATH)
        if p not in sys.path and os.path.isdir(p):
            sys.path.insert(0, p)
        from easyhap.world_map import get_world_lines
        return get_world_lines()
    except Exception:  # noqa: BLE001 — 缺 EasyHap / 数据损坏一律降级
        return None


def draw_world(ax, borders: bool = False, color: str = "#9a9a9a",
               border_color: str = "#c0c0c0", lw: float = 0.5,
               zorder: int = 1) -> bool:
    """在 ax 上画世界海岸线 (可选国界)。成功 True / EasyHap 不可用 False。

    调用方先设好经纬度 xlim/ylim (本函数不 set_xlim, 只加线, 不改视图)。
    """
    wl = _world_lines()
    if not wl:
        return False

    def _add(key, c, w):
        for segs in wl.get(key, []):
            # 每条折线可含多段 (MultiLineString 风格): [[ [lon,lat],... ], ...]
            first = segs[0] if segs and isinstance(segs[0], (list, tuple)) else None
            if first and isinstance(first[0], (list, tuple)):
                for sub in segs:
                    if len(sub) >= 2:
                        xs, ys = zip(*sub)
                        ax.plot(xs, ys, color=c, lw=w, zorder=zorder)
            elif segs and len(segs) >= 2:
                xs, ys = zip(*segs)
                ax.plot(xs, ys, color=c, lw=w, zorder=zorder)

    _add("coastlines", color, lw)
    if borders:
        _add("countries", border_color, max(0.3, lw - 0.2))
    return True


def available() -> bool:
    """EasyHap 世界折线是否可加载 (供调用方预先决定是否标注海岸线来源)。"""
    return _world_lines() is not None


if __name__ == "__main__":
    # 烟测: 画全球海岸线存 PNG
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(10, 5))
    ok = draw_world(ax, borders=True)
    ax.set_xlim(-180, 180)
    ax.set_ylim(-60, 80)
    ax.set_aspect("equal")
    fig.savefig("world_basemap_smoke.png", dpi=110)
    print("draw:", ok, "->", os.path.abspath("world_basemap_smoke.png"))
