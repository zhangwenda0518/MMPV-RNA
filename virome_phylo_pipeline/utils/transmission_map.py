#!/usr/bin/env python3
"""
utils/transmission_map.py — 发表级静态传播图 + 事件级时间动画
==============================================================
吃 utils/transmission_paths.py 产出的 branch_table.csv（或现场解析 MCC 树），
在真实经纬度上画"贝塞尔弧线 + 末段箭头 + 时间色带"的传播图，并可输出
时间轴推进 GIF（残影拖尾 + blended pacing）。

方法学出处（逐条可追溯）
----------------------
* 二次贝塞尔弧线（垂直偏移中点控制点，沿曲线插值 age）：
  git-repo/ggphylogeo/R/geom_phylo_branches.R:33-72。
* 自适应曲率（近距离放大弯曲，避免短弧视觉重叠）：
  git-repo/flu_d_project/analysis/migration_map.ipynb cell1
  `adjust_bezier_control`: k/sqrt(sqrt(|dx|)) 思想。
* 末段箭头（只给曲线末段画箭头，避免整条曲线多段箭头）：
  phymapr R/mapping.R:262-303 技巧 + phylomovie R/plot_mcc_tree_with_hpds.R:155
  的 graphics::arrows 语义（本项目用 matplotlib annotate 重写）。
* 时间色带（HSV 随节点时间渐变）：phylomovie plot_mcc_tree_with_hpds.R:106-117，
  本项目用 matplotlib colormap 等价实现。
* 残影拖尾（当前帧实心、前 lagg 个切片空心淡点叠加）：
  phylomovie vignettes/make_phylomovie.Rmd:158-195。
* blended pacing（时间/秩混合步进，解决采样时间不均时大部分帧无事件）：
  phymapr R/epidemiologic_inference.R:274-313 apply_blended_pacing。
* 帧时长 50ms 下限钳制（卡顿时慢放而非跳帧）：
  git-repo/spreadgl2.github.io/src/features/timeline/playback.ts:18-29
  MAX_FRAME_DELTA_MS 语义。

底图合规（红线）
--------------
默认**无任何在线瓦片/底图**（纯矢量点线，离线合规）。需要国界/省界时，
用 `--boundaries` 传入**自有审图号来源**的 GeoJSON（自然资源部/天地图下载），
本模块只做矢量绘制，不内置、不下载任何第三方边界或瓦片（不使用
Natural Earth / CartoDB / OSM —— 无审图号，不合中国发表要求）。

自包含性：numpy/pandas/matplotlib/Pillow + 标准库；坐标解析优先用户 CSV，
否则走 utils.geo_resolver 离线词典（11960 条，默认离线）。
"""

from __future__ import annotations

import io
import json
import math
import os
import re
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

try:
    from utils.virphy_bridge import LogCollector
except ImportError:  # 允许直接运行
    from .virphy_bridge import LogCollector  # type: ignore


# ─────────────────────────────────────────────────────────────────────
# 坐标
# ─────────────────────────────────────────────────────────────────────

def haversine_km(lat1: float, lon1: float, lat2: float, lon2: float) -> float:
    """球面大圆距离 km（IUGG 平均半径 6371.0088；与 geopy geodesic 1% 内一致）。"""
    r = 6371.0088
    p1, p2 = math.radians(lat1), math.radians(lat2)
    dp = p2 - p1
    dl = math.radians(lon2 - lon1)
    a = math.sin(dp / 2) ** 2 + math.cos(p1) * math.cos(p2) * math.sin(dl / 2) ** 2
    return 2 * r * math.asin(math.sqrt(a))


def resolve_coordinates(
    locations: List[str],
    coord_csv: Optional[str] = None,
    log: Optional[LogCollector] = None,
) -> Dict[str, Optional[Tuple[float, float]]]:
    """地点 → (lat, lon)。优先 coord_csv（列: location,lat,lon），否则离线 geo_resolver。"""
    log = log or LogCollector()
    out: Dict[str, Optional[Tuple[float, float]]] = {}
    csv_map: Dict[str, Tuple[float, float]] = {}
    if coord_csv:
        df = pd.read_csv(coord_csv, encoding="utf-8-sig")
        for _, row in df.iterrows():
            k = str(row.get("location", "") or "").strip()
            try:
                csv_map[k] = (float(row["lat"]), float(row["lon"]))
            except (TypeError, ValueError, KeyError):
                continue
    use_resolver = False
    for loc in locations:
        if loc in csv_map:
            out[loc] = csv_map[loc]
            continue
        if not use_resolver:
            try:
                from utils.geo_resolver import resolve_location  # 离线词典，默认不走网络
                use_resolver = True
            except Exception as e:  # noqa: BLE001
                log.warning(f"geo_resolver 不可用 ({e})，无 CSV 的地点将无坐标")
        try:
            from utils.geo_resolver import resolve_location
            coord = resolve_location(loc)
        except Exception:  # noqa: BLE001
            coord = None
        out[loc] = coord
        if coord is None:
            log.warning(f"地点无坐标: {loc} (请用 --coords 提供 location,lat,lon)")
    return out


# ─────────────────────────────────────────────────────────────────────
# 弧线数学
# ─────────────────────────────────────────────────────────────────────

def bezier_arc_points(
    x0: float, y0: float, x1: float, y1: float,
    curvature: float = 0.2,
    adaptive: bool = True,
    n: int = 60,
) -> Tuple[np.ndarray, np.ndarray]:
    """二次贝塞尔弧线采样点（ggphylogeo geom_phylo_branches.R:33-72）。

    控制点 = 中点 + curvature × 垂直方向（cx=midx-c·dy, cy=midy+c·dx）。
    adaptive=True 时按 flu_d_project 思想对近距离弧放大曲率：c ∝ k/√√|dx|。
    """
    dx, dy = x1 - x0, y1 - y0
    c = curvature
    if adaptive:
        scale = math.sqrt(math.sqrt(max(abs(dx), 1e-9)))
        c = curvature / scale * math.sqrt(math.sqrt(4.0))  # 归一基准: |dx|=4° 时曲率=curvature
    cx = (x0 + x1) / 2 - c * dy
    cy = (y0 + y1) / 2 + c * dx
    t = np.linspace(0.0, 1.0, n)
    bx = (1 - t) ** 2 * x0 + 2 * (1 - t) * t * cx + t ** 2 * x1
    by = (1 - t) ** 2 * y0 + 2 * (1 - t) * t * cy + t ** 2 * y1
    return bx, by


def _age_colormap(ages: np.ndarray):
    """时间色带（早=深，晚=亮），phylomovie HSV 时间渐变的 matplotlib 等价。"""
    import matplotlib.pyplot as plt
    return plt.get_cmap("viridis")


# ─────────────────────────────────────────────────────────────────────
# 边界 GeoJSON（用户自备审图号来源；本模块不内置不下载数据）
# ─────────────────────────────────────────────────────────────────────

def _draw_boundaries(ax, geojson_path: str):
    with open(geojson_path, encoding="utf-8") as f:
        gj = json.load(f)
    from matplotlib.patches import Polygon as MplPolygon
    from matplotlib.collections import PatchCollection

    def add_ring(ring):
        pts = [(p[0], p[1]) for p in ring]
        if len(pts) >= 3:
            patches.append(MplPolygon(pts, closed=True))

    patches: list = []
    feats = gj.get("features", []) if isinstance(gj, dict) else gj
    for feat in feats:
        geom = feat.get("geometry", feat) if isinstance(feat, dict) else {}
        gtype = geom.get("type")
        coords = geom.get("coordinates")
        if gtype == "Polygon":
            for ring in coords or []:
                add_ring(ring)
        elif gtype == "MultiPolygon":
            for poly in coords or []:
                for ring in poly:
                    add_ring(ring)
    if patches:
        coll = PatchCollection(patches, facecolor="none", edgecolor="#8a8a8a",
                               linewidths=0.5, zorder=1)
        ax.add_collection(coll)


# ─────────────────────────────────────────────────────────────────────
# 静态图
# ─────────────────────────────────────────────────────────────────────

def draw_static_map(
    branch_table: pd.DataFrame,
    coords: Dict[str, Optional[Tuple[float, float]]],
    out_path: str,
    boundaries_geojson: Optional[str] = None,
    title: str = "Phylogeographic transmission",
    curvature: float = 0.2,
    stream: bool = True,
    dpi: int = 300,
    log: Optional[LogCollector] = None,
) -> str:
    """发表级静态传播图（PNG + 按扩展名另存 PDF/SVG 向量版）。

    弧线 = 地点变化的分支；颜色 = 分支中点年龄时间色带；线宽 ∝ 期望替换数（若有）；
    stream=True（默认）沿方向渐显（ggphylogeo "stream" 式）；箭头只画末段。
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    log = log or LogCollector()

    pw = branch_table[branch_table.get("location_change", pd.Series(dtype=bool)) == True]  # noqa: E712
    pw = pw[pw["start_location"].map(lambda x: coords.get(x) is not None)
            & pw["end_location"].map(lambda x: coords.get(x) is not None)].copy()
    if pw.empty:
        raise SystemExit("无可绘制的传播事件（地点变化分支为 0 或坐标缺失）")

    # 上色年龄 = 分支**中点**（2026-09-16 修：原用分支起点 47.4 年前 → 视觉上读成
    # "最早迁往内蒙古"，用户指正。地点变迁时刻只被约束在 [child_age, parent_age]
    # 区间内，中点是其最朴素的点估计；区间宽度本身记在 pathways.csv 的 duration）
    pw = pw.assign(_color_age=(pw["start_age"].astype(float)
                               + pw["end_age"].astype(float)) / 2.0)

    fig, ax = plt.subplots(figsize=(11, 9))
    if boundaries_geojson:
        _draw_boundaries(ax, boundaries_geojson)

    ages = pw["_color_age"].values
    amin, amax = float(np.nanmin(ages)), float(np.nanmax(ages))
    cmap = _age_colormap(ages)
    norm = plt.Normalize(vmin=amin, vmax=max(amax, amin + 1e-9))

    subs = pw["expected_subs"] if "expected_subs" in pw else None
    if subs is not None and subs.notna().any():
        w = 0.8 + 2.2 * (subs.fillna(subs.median()) / max(float(subs.median()), 1e-9)).clip(upper=3)
    else:
        w = pd.Series([1.6] * len(pw), index=pw.index)

    # 地点散点
    seen: set = set()
    for _, row in pw.iterrows():
        for loc in (row["start_location"], row["end_location"]):
            if loc not in seen:
                seen.add(loc)
    for loc in sorted(seen):
        lat, lon = coords[loc]
        ax.scatter(lon, lat, s=90, c="#d64541", edgecolors="#2c3e50",
                   linewidths=1.2, zorder=5)
        ax.annotate(loc, (lon, lat), textcoords="offset points",
                    xytext=(6, 6), fontsize=10, zorder=6)

    # 弧线 + 末段箭头
    for (_, row), lw in zip(pw.iterrows(), w):
        la0, lo0 = coords[row["start_location"]]
        la1, lo1 = coords[row["end_location"]]
        bx, by = bezier_arc_points(lo0, la0, lo1, la1, curvature=curvature)
        color = cmap(norm(float(row["_color_age"])))
        if stream:
            # ggphylogeo "stream" 式渐显: 沿方向 alpha/线宽渐增 (方向线索), 末段箭头收尾
            n_seg = max(8, len(bx) // 4)
            for i in range(1, n_seg + 1):
                t0, t1 = (i - 1) / n_seg, i / n_seg
                j0, j1 = int(t0 * (len(bx) - 1)), int(t1 * (len(bx) - 1))
                ax.plot(bx[j0:j1 + 1], by[j0:j1 + 1], color=color,
                        linewidth=float(lw) * (0.5 + 0.5 * t1),
                        alpha=0.25 + 0.65 * t1, zorder=3, solid_capstyle="round")
        else:
            ax.plot(bx, by, color=color, linewidth=float(lw), alpha=0.85, zorder=3)
        # 末段箭头：只取曲线最后 2 个采样点（phymapr mapping.R 技巧）
        ax.annotate("", xy=(bx[-1], by[-1]), xytext=(bx[-3], by[-3]),
                    arrowprops=dict(arrowstyle="-|>", color=color, lw=0, mutation_scale=16),
                    zorder=4)

    # 地点置信环 (phylowood 饼图思路的静态化): 祖先后验均值 < 0.9 的地点画警示环,
    # 后验越低环越实 —— 把 location.states.set.prob 的不确定度画到图上
    probs: Dict[str, list] = {}
    for _, row in pw.iterrows():
        for lc, pb in ((row["start_location"], row.get("start_location_prob")),
                       (row["end_location"], row.get("end_location_prob"))):
            if pb is not None and pd.notna(pb):
                probs.setdefault(lc, []).append(float(pb))
    for loc, pv in probs.items():
        p = float(np.mean(pv))
        if p < 0.9:
            lat, lon = coords[loc]
            ax.scatter(lon, lat, s=340, facecolors="none", edgecolors="#2c3e50",
                       linewidths=1.6, linestyles="dashed", alpha=min(1.0, 1.1 - p),
                       zorder=4)
            ax.annotate(f"p={p:.2f}", (lon, lat), textcoords="offset points",
                        xytext=(8, -12), fontsize=8, color="#555", zorder=6)

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, shrink=0.75, pad=0.02)
    cbar.set_label("Age of transmission branch midpoint (years before present)\n"
                   "(transition time constrained within branch span)")
    ax.set_title(title, fontsize=13, fontweight="bold")
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.set_aspect("equal", adjustable="datalim")
    ax.autoscale()
    ax.margins(0.12)
    fig.tight_layout()

    os.makedirs(os.path.dirname(os.path.abspath(out_path)), exist_ok=True)
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
    root, ext = os.path.splitext(out_path)
    for ve in (".pdf", ".svg"):
        if ext.lower() != ve:
            fig.savefig(root + ve, bbox_inches="tight")
    plt.close(fig)
    log.emit(f"静态传播图: {out_path} (+pdf/svg), {len(pw)} 条弧线")
    return out_path


# ─────────────────────────────────────────────────────────────────────
# 动画
# ─────────────────────────────────────────────────────────────────────

def blended_pacing(times: np.ndarray, lam: float = 0.5) -> np.ndarray:
    """时间/秩混合步进（phymapr apply_blended_pacing, L274-313）。

    u=线性时间位置, r=秩位置, t' = t_min + (λu + (1-λ)r)(t_max - t_min)。
    λ=1 纯时间；λ=0 纯秩（每帧必出事件）。
    """
    order = np.argsort(times)                      # 早→晚
    t_sorted = times[order]
    t_min, t_max = float(t_sorted.min()), float(t_sorted.max())
    span = max(t_max - t_min, 1e-9)
    u = (t_sorted - t_min) / span
    r = np.linspace(0.0, 1.0, len(t_sorted)) if len(t_sorted) > 1 else np.array([0.5])
    mixed = lam * u + (1 - lam) * r
    out = np.empty_like(mixed)
    out[order] = t_min + mixed * span
    return out


def draw_animation(
    branch_table: pd.DataFrame,
    coords: Dict[str, Optional[Tuple[float, float]]],
    out_gif: str,
    n_frames: int = 48,
    trail: Optional[int] = None,
    pacing_lambda: float = 0.5,
    fps: int = 8,
    boundaries_geojson: Optional[str] = None,
    title: str = "Phylogeographic transmission (time reveal)",
    curvature: float = 0.2,
    log: Optional[LogCollector] = None,
) -> str:
    """事件级时间揭示 GIF（phymapr blended pacing + 可选 phylomovie 残影拖尾）。

    trail 语义（2026-09-16 修订，默认改为**累积模式**）:
      None / -1  — 累积：路线按时间依次出现后**永久保留**（发表默认。原默认
                   trail=5 会把最早——往往最重要——的扩散事件在末帧淡到
                   alpha=0.15 近乎不可见，真实 gcva32 末帧实测丢失内蒙古弧线，
                   用户指出后改定）。
      int > 0    — 残影：只高亮最近 trail 条，更早的按衰减叠加（观感用）。
    帧时长 50ms 下限钳制（spreadgl2 playback.ts 语义）。
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from PIL import Image
    log = log or LogCollector()

    cumulative = trail is None or trail < 0

    pw = branch_table[branch_table.get("location_change", pd.Series(dtype=bool)) == True]  # noqa: E712
    pw = pw[pw["start_location"].map(lambda x: coords.get(x) is not None)
            & pw["end_location"].map(lambda x: coords.get(x) is not None)].copy()
    if pw.empty:
        raise SystemExit("无可动画的传播事件")
    pw = pw.sort_values("start_age", ascending=False).reset_index(drop=True)  # 早→晚

    reveal = blended_pacing(pw["start_age"].astype(float).values, lam=pacing_lambda)
    t_min, t_max = float(reveal.min()), float(reveal.max())

    cmap = _age_colormap(pw["start_age"].astype(float).values)
    amin, amax = float(pw["start_age"].min()), float(pw["start_age"].max())
    norm = plt.Normalize(vmin=amin, vmax=max(amax, amin + 1e-9))

    # 静态坐标框架
    all_lats = [coords[l][0] for l in set(pw["start_location"]) | set(pw["end_location"])]
    all_lons = [coords[l][1] for l in set(pw["start_location"]) | set(pw["end_location"])]
    lat_pad = max(1.0, (max(all_lats) - min(all_lats)) * 0.15)
    lon_pad = max(1.0, (max(all_lons) - min(all_lons)) * 0.15)

    frames: List[Image.Image] = []
    for k in range(n_frames):
        tau = t_min + (t_max - t_min) * k / max(1, n_frames - 1)
        fig, ax = plt.subplots(figsize=(10, 8))
        if boundaries_geojson:
            _draw_boundaries(ax, boundaries_geojson)
        for loc in sorted(set(pw["start_location"]) | set(pw["end_location"])):
            lat, lon = coords[loc]
            ax.scatter(lon, lat, s=70, c="#d64541", edgecolors="#2c3e50",
                       linewidths=1.0, zorder=5)
            ax.annotate(loc, (lon, lat), textcoords="offset points",
                        xytext=(5, 5), fontsize=9, zorder=6)

        revealed = np.where(reveal <= tau + 1e-9)[0]
        newest = revealed[-1] if len(revealed) else None
        for j in revealed:
            row = pw.iloc[j]
            la0, lo0 = coords[row["start_location"]]
            la1, lo1 = coords[row["end_location"]]
            bx, by = bezier_arc_points(lo0, la0, lo1, la1, curvature=curvature)
            color = cmap(norm(float(row["start_age"])))
            slot = len(revealed) - 1 - int(np.where(revealed == j)[0][0])  # 0=最新
            if cumulative:
                alpha, lw = (1.0 if slot == 0 else 0.85), (2.2 if slot == 0 else 1.6)
            else:
                alpha = 1.0 if slot == 0 else max(0.3, 0.8 - 0.12 * slot)  # 残影衰减(下限抬高)
                lw = 2.2 if slot == 0 else max(1.0, 1.8 - 0.15 * slot)
            ax.plot(bx, by, color=color, linewidth=lw,
                    alpha=alpha, zorder=3)
            if slot == 0:
                ax.annotate("", xy=(bx[-1], by[-1]), xytext=(bx[-3], by[-3]),
                            arrowprops=dict(arrowstyle="-|>", color=color, lw=0,
                                            mutation_scale=15), zorder=4)
                ax.scatter([bx[-1]], [by[-1]], s=(60 + 40 * math.sin(k / 2.0)),
                           facecolors="none", edgecolors=color, linewidths=1.6, zorder=6)

        ax.set_xlim(min(all_lons) - lon_pad, max(all_lons) + lon_pad)
        ax.set_ylim(min(all_lats) - lat_pad, max(all_lats) + lat_pad)
        ax.set_title(f"{title}   t = {tau:.2f} yr", fontsize=12, fontweight="bold")
        ax.set_xlabel("Longitude")
        ax.set_ylabel("Latitude")
        ax.set_aspect("equal", adjustable="box")
        buf = io.BytesIO()
        fig.savefig(buf, format="png", dpi=110, bbox_inches="tight")
        plt.close(fig)
        buf.seek(0)
        frames.append(Image.open(buf).convert("P", palette=Image.ADAPTIVE))

    # 帧时长: 总时长≈n_frames/fps, 单帧 50ms 下限钳制（spreadgl2 playback.ts 语义）
    frame_ms = max(50, int(1000 / max(1, fps)))
    os.makedirs(os.path.dirname(os.path.abspath(out_gif)), exist_ok=True)
    frames[0].save(out_gif, save_all=True, append_images=frames[1:],
                   duration=frame_ms, loop=0, optimize=True)
    mode = "cumulative" if cumulative else f"ghost-trail({trail})"
    log.emit(f"动画 GIF: {out_gif} ({n_frames} 帧, {frame_ms}ms/帧, mode={mode}, "
             f"pacing λ={pacing_lambda})")
    return out_gif


# ─────────────────────────────────────────────────────────────────────
# 扩散进展多边形 + 时间快照面板 (seraphim R_script_analyses.r §2 / phylomovie 静态化)
# ─────────────────────────────────────────────────────────────────────

def _convex_hull_poly(pts: List[Tuple[float, float]]) -> Optional[List[Tuple[float, float]]]:
    """凸包顶点（Andrew 单调链, 纯实现避免 scipy 依赖差异）；点数<3 返回 None。"""
    pts = sorted(set(pts))
    if len(pts) < 3:
        return None

    def cross(o, a, b):
        return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])

    lower: List[Tuple[float, float]] = []
    for p in pts:
        while len(lower) >= 2 and cross(lower[-2], lower[-1], p) <= 0:
            lower.pop()
        lower.append(p)
    upper: List[Tuple[float, float]] = []
    for p in reversed(pts):
        while len(upper) >= 2 and cross(upper[-2], upper[-1], p) <= 0:
            upper.pop()
        upper.append(p)
    return lower[:-1] + upper[:-1]


def _active_positions(bt: pd.DataFrame, coords: Dict, age_hi: float) -> List[Tuple[float, float]]:
    """截至年龄 age_hi（含）已出现的传播事件的端点坐标集合（累积口径）。

    返回 (lon, lat) 顺序 —— 与 plt.Polygon 的 (x, y) 对齐。
    （2026-09-16 修：原返回 (lat, lon) 直接喂 Polygon，凸包画到经纬倒置的
    错误位置——真实数据上表现为多边形出现在欧洲上空，用户可见即修。）"""
    pts: List[Tuple[float, float]] = []
    for _, r in bt.iterrows():
        mid = (float(r["start_age"]) + float(r["end_age"])) / 2.0
        if mid >= age_hi:  # age 大 = 早; mid >= 阈值 = 已发生
            for lc in (r["start_location"], r["end_location"]):
                c = coords.get(lc)
                if c:
                    pts.append((c[1], c[0]))
    return pts


def draw_progress_polys(
    branch_table: pd.DataFrame,
    coords: Dict[str, Optional[Tuple[float, float]]],
    out_path: str,
    n_steps: int = 4,
    boundaries_geojson: Optional[str] = None,
    dpi: int = 300,
    log: Optional[LogCollector] = None,
) -> str:
    """扩散进展多边形叠加图（seraphim wnv R_script_analyses.r §2 凸包进展的静态版）。

    按"传播事件中点年龄"取 n_steps 个累积时点；每个时点对**已到达位置集合**画
    凸包半透明多边形（时间色带嵌套）—— 一图看清波及面如何随时间扩张。
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    log = log or LogCollector()

    pw = branch_table[branch_table.get("location_change", pd.Series(dtype=bool)) == True]  # noqa: E712
    pw = pw[pw["start_location"].map(lambda x: coords.get(x) is not None)
            & pw["end_location"].map(lambda x: coords.get(x) is not None)].copy()
    if pw.empty:
        raise SystemExit("无可绘制事件")
    pw = pw.assign(_mid=(pw["start_age"].astype(float) + pw["end_age"].astype(float)) / 2.0)
    mids = sorted(pw["_mid"], reverse=True)  # 早→晚
    cuts = [mids[int(len(mids) * i / n_steps)] for i in range(n_steps)]

    cmap = _age_colormap(np.array(cuts + [cuts[-1]]))
    amin, amax = float(min(cuts)), float(max(cuts))
    norm = plt.Normalize(vmin=amin, vmax=max(amax, amin + 1e-9))

    fig, ax = plt.subplots(figsize=(11, 9))
    if boundaries_geojson:
        _draw_boundaries(ax, boundaries_geojson)
    all_pts = _active_positions(pw, coords, amin - 1e9)
    for loc in sorted({l for l in pw["start_location"]} | {l for l in pw["end_location"]}):
        lat, lon = coords[loc]
        ax.scatter(lon, lat, s=90, c="#d64541", edgecolors="#2c3e50",
                   linewidths=1.2, zorder=5)
        ax.annotate(loc, (lon, lat), textcoords="offset points", xytext=(6, 6),
                    fontsize=10, zorder=6)
    for i, cut in enumerate(cuts):
        pts = _active_positions(pw, coords, cut)
        hull = _convex_hull_poly(pts)
        if not hull:
            continue
        color = cmap(norm(cut))
        poly = plt.Polygon(hull, closed=True, facecolor=color, alpha=0.14 + 0.06 * i,
                           edgecolor=color, linewidth=1.4, linestyle="--", zorder=2 + i)
        ax.add_patch(poly)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, shrink=0.75, pad=0.02)
    cbar.set_label("Spread-stage age (branch-midpoint, years before present)")
    ax.set_title("Progressive spread area (convex hull per stage)", fontsize=13,
                 fontweight="bold")
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.set_aspect("equal", adjustable="datalim")
    ax.autoscale()
    ax.margins(0.12)
    fig.tight_layout()
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
    root, _ = os.path.splitext(out_path)
    fig.savefig(root + ".pdf", bbox_inches="tight")
    plt.close(fig)
    log.emit(f"进展多边形图: {out_path} (+pdf), {len(cuts)} 个扩散阶段")
    return out_path


def draw_snapshot_panels(
    branch_table: pd.DataFrame,
    coords: Dict[str, Optional[Tuple[float, float]]],
    out_path: str,
    n_panels: int = 4,
    boundaries_geojson: Optional[str] = None,
    dpi: int = 300,
    log: Optional[LogCollector] = None,
) -> str:
    """时间快照多面板（phylomovie 时间切片的静态发表版）：每格 = 截至该时刻
    已发生的全部路线（累积口径），标题标注时段。GIF 的可嵌入替代。"""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    log = log or LogCollector()

    pw = branch_table[branch_table.get("location_change", pd.Series(dtype=bool)) == True]  # noqa: E712
    pw = pw[pw["start_location"].map(lambda x: coords.get(x) is not None)
            & pw["end_location"].map(lambda x: coords.get(x) is not None)].copy()
    if pw.empty:
        raise SystemExit("无可绘制事件")
    pw = pw.assign(_mid=(pw["start_age"].astype(float) + pw["end_age"].astype(float)) / 2.0)
    pw = pw.sort_values("_mid", ascending=False).reset_index(drop=True)  # 早→晚
    mids = pw["_mid"].tolist()
    cuts = [mids[min(int(len(mids) * (i + 1) / n_panels), len(mids) - 1)]
            for i in range(n_panels)]
    cmap = _age_colormap(pw["_mid"].values)
    amin, amax = float(pw["_mid"].min()), float(pw["_mid"].max())
    norm = plt.Normalize(vmin=amin, vmax=max(amax, amin + 1e-9))

    all_lats = [coords[l][0] for l in set(pw["start_location"]) | set(pw["end_location"])]
    all_lons = [coords[l][1] for l in set(pw["start_location"]) | set(pw["end_location"])]
    lat_pad = max(1.0, (max(all_lats) - min(all_lats)) * 0.18)
    lon_pad = max(1.0, (max(all_lons) - min(all_lons)) * 0.18)

    ncol = min(n_panels, 2)
    nrow = (n_panels + ncol - 1) // ncol
    fig, axes = plt.subplots(nrow, ncol, figsize=(6.5 * ncol, 5.6 * nrow), squeeze=False)
    for k in range(n_panels):
        ax = axes[k // ncol][k % ncol]
        cut = cuts[k]
        if boundaries_geojson:
            _draw_boundaries(ax, boundaries_geojson)
        sub = pw[pw["_mid"] >= cut]
        for loc in sorted(set(pw["start_location"]) | set(pw["end_location"])):
            lat, lon = coords[loc]
            shown = (loc in set(sub["start_location"]) | set(sub["end_location"]))
            ax.scatter(lon, lat, s=60 if shown else 14,
                       c="#d64541" if shown else "#c9c9c9",
                       edgecolors="#2c3e50" if shown else "#bdbdbd",
                       linewidths=0.8, zorder=5)
            if shown:
                ax.annotate(loc, (lon, lat), textcoords="offset points", xytext=(5, 5),
                            fontsize=8, zorder=6)
        for _, row in sub.iterrows():
            la0, lo0 = coords[row["start_location"]]
            la1, lo1 = coords[row["end_location"]]
            bx, by = bezier_arc_points(lo0, la0, lo1, la1, curvature=0.2)
            color = cmap(norm(float(row["_mid"])))
            ax.plot(bx, by, color=color, linewidth=1.8, alpha=0.9, zorder=3)
            ax.annotate("", xy=(bx[-1], by[-1]), xytext=(bx[-3], by[-3]),
                        arrowprops=dict(arrowstyle="-|>", color=color, lw=0,
                                        mutation_scale=13), zorder=4)
        hull = _convex_hull_poly(_active_positions(pw, coords, cut))
        if hull:
            ax.add_patch(plt.Polygon(hull, closed=True, facecolor="#34495e",
                                     alpha=0.06, edgecolor="#7f8c8d",
                                     linewidth=1.0, linestyle=":", zorder=1))
        ax.set_xlim(min(all_lons) - lon_pad, max(all_lons) + lon_pad)
        ax.set_ylim(min(all_lats) - lat_pad, max(all_lats) + lat_pad)
        prev_cut = cuts[k - 1] if k else float("inf")
        n_new = int((pw["_mid"] <= prev_cut).sum() - (pw["_mid"] <= cut).sum())
        ax.set_title(f"by ~{cut:.1f} yr  ({n_new} new event"
                     f"{'s' if n_new != 1 else ''})", fontsize=10)
        ax.set_aspect("equal", adjustable="box")
    for k in range(n_panels, nrow * ncol):
        axes[k // ncol][k % ncol].axis("off")
    fig.suptitle("Transmission snapshots (cumulative by branch-midpoint age)",
                 fontsize=13, fontweight="bold")
    fig.tight_layout()
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
    root, _ = os.path.splitext(out_path)
    fig.savefig(root + ".pdf", bbox_inches="tight")
    plt.close(fig)
    log.emit(f"快照多面板: {out_path} (+pdf), {n_panels} 格")
    return out_path


# ─────────────────────────────────────────────────────────────────────
# 树-图联动双面板 (phymapr generate_phylo_transmission 图版替代, 纯 matplotlib)
# ─────────────────────────────────────────────────────────────────────

def _to_decimal_year(val) -> Optional[float]:
    """YYYY / YYYY-MM / YYYY-MM-DD → 十进制年（YYYY-MM 取月中）。"""
    import calendar as _cal
    import datetime as _dt
    val = str(val or "").strip()
    if not val or val.lower() in ("nan", "none", "null"):
        return None
    val = val.replace("/", "-")
    parts = val.split("-")
    try:
        if len(parts) >= 3:
            y, m, d = int(parts[0]), int(parts[1]), int(parts[2])
            doy = (_dt.date(y, m, d) - _dt.date(y, 1, 1)).days
            return y + doy / (366.0 if _cal.isleap(y) else 365.0)
        if len(parts) == 2 and len(parts[0]) == 4:
            y, m = int(parts[0]), int(parts[1])
            return y + (m - 0.5) / 12.0
        if len(parts[0]) == 4:
            return float(parts[0])
    except (ValueError, TypeError):
        return None
    return None


def _categorical_palette(locations: List[str]) -> Dict[str, tuple]:
    """地点 → 分类色（tab20；两面板共用的 Region 配色）。"""
    import matplotlib.pyplot as plt
    cmap = plt.get_cmap("tab20")
    return {loc: cmap(i % 20) for i, loc in enumerate(sorted(locations))}


def _ladderize_children(nodes: List[Dict], edges: List[Dict]) -> Dict[str, List[str]]:
    """children 表 + 按子树 tip 数排序（ladderize，图面整齐）。"""
    from collections import defaultdict
    children = defaultdict(list)
    for e in edges:
        children[e["parent_id"]].append(e["child_id"])

    size: Dict[str, int] = {}

    def subtree_size(nid: str) -> int:
        if nid in size:
            return size[nid]
        s = 1 + sum(subtree_size(c) for c in children.get(nid, []))
        size[nid] = s
        return s

    root = nodes[0]["node_id"]
    subtree_size(root)
    for k in children:
        children[k].sort(key=lambda c: size.get(c, 0))
    return children


def _load_rename_map(x) -> Dict[str, str]:
    """重命名映射：dict 直接用；文件按行解析（制表/逗号/⇥ 分隔，# 注释跳过）。

    格式: 旧ID <sep> 新ID（与平台 🏷 序列 ID 重命名工具的「旧ID⇥新ID」一致）。
    """
    if not x:
        return {}
    if isinstance(x, dict):
        return dict(x)
    m: Dict[str, str] = {}
    for line in Path(x).read_text(encoding="utf-8-sig").splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        for sep in ("\t", "⇥", ",", " "):
            if sep in line:
                old, _, new = line.partition(sep)
                old, new = old.strip(), new.strip()
                if old and new:
                    m[old] = new
                break
    return m


def draw_tree_map(
    mcc_tree: str,
    metadata_csv: str,
    coords: Optional[Dict[str, Optional[Tuple[float, float]]]],
    out_path: str,
    pathways_csv: Optional[str] = None,
    boundaries_geojson: Optional[str] = None,
    rename_map: Optional[object] = None,
    curvature: float = 0.2,
    dpi: int = 300,
    log: Optional[LogCollector] = None,
) -> str:
    """树-图联动双面板（phymapr "Tree|Map" 图版的 matplotlib 替代）。

    rename_map: dict 或映射文件路径（旧ID→新ID）——tip 标签换成短名
    （发表图美化；未映射的 tip 保持原名）。

    coords=None 时自动走离线 geo_resolver 解析树内地点。
    左panel: 时间标尺 MCC 树（x=采样日期；tip/内部节点按地点配色，内部节点
             透明度=地点后验）；右panel: 传播地图（弧线按终点地点配色、
             线型=三级置信: 实线 Direct / 虚线 Indirect / 点线 Distant-Import）。
    两panel共用 Region 配色图例。出处: phymapr R/mapping.R+plot_frontpage 的
    图版设计（ggtree/ggplot 实现），本函数为 matplotlib 重写。
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.gridspec import GridSpec
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    from collections import defaultdict
    log = log or LogCollector()

    from utils.transmission_paths import extract_tree_events
    ev = extract_tree_events(mcc_tree, metadata_csv)
    nodes, edges = ev["nodes"], ev["edges"]
    by_id = {n["node_id"]: n for n in nodes}
    children = _ladderize_children(nodes, edges)
    if coords is None:
        locs = sorted({n["location"] for n in nodes if n.get("location")})
        coords = resolve_coordinates(locs, log=log)
    rmap = _load_rename_map(rename_map)

    # ── present = 元数据最晚采样日期（十进制年）──
    meta = pd.read_csv(metadata_csv, encoding="utf-8-sig")
    dcol = next((c for c in meta.columns if c.strip().lower() == "date"), None)
    dates: Dict[str, float] = {}
    if dcol is not None:
        for _, r in meta.iterrows():
            dy = _to_decimal_year(r.get(dcol))
            nm = str(r.get("name", "") or "").strip()
            if dy and nm:
                dates[nm] = dy
    tip_dates = [dates.get(n["name"]) for n in nodes if n["is_tip"] and n["name"]]
    tip_dates = [d for d in tip_dates if d is not None]
    present = max(tip_dates) if tip_dates else None
    if present is None:
        raise SystemExit("元数据无可用日期 → 无法标定时间轴（需要 date 列）")

    # ── y 布局: tip 按 DFS 序编号, 内部节点=子节点均值 ──
    tip_y: Dict[str, float] = {}
    y0 = 0
    for n in nodes:
        if n["is_tip"]:
            tip_y[n["node_id"]] = float(y0)
            y0 += 1
    node_y: Dict[str, float] = {}

    def place(nid: str) -> float:
        if nid in node_y:
            return node_y[nid]
        if nid in tip_y:
            node_y[nid] = tip_y[nid]
            return tip_y[nid]
        ys = [place(c) for c in children.get(nid, [])]
        node_y[nid] = sum(ys) / len(ys) if ys else 0.0
        return node_y[nid]

    place(nodes[0]["node_id"])
    node_x = {n["node_id"]: (present - n["age"] if n["age"] is not None else None)
              for n in nodes}

    # ── 置信分级映射（pathways.csv 的 branch_id → confidence）──
    conf_by_branch: Dict[str, str] = {}
    if pathways_csv and os.path.exists(pathways_csv):
        pwc = pd.read_csv(pathways_csv)
        if {"branch_id", "confidence"}.issubset(pwc.columns):
            conf_by_branch = dict(zip(pwc["branch_id"], pwc["confidence"]))

    locs: set = set()
    for n in nodes:
        if n.get("location"):
            locs.add(n["location"])
    for ed in edges:
        for nid in (ed["parent_id"], ed["child_id"]):
            loc = by_id[nid].get("location")
            if loc:
                locs.add(loc)
    palette = _categorical_palette(sorted(locs))

    conf_style = {"direct": ("-", 2.4), "indirect": ("--", 1.9),
                  "distant_import": (":", 1.7), "unclassified": ("-", 1.4)}

    fig = plt.figure(figsize=(21, 9))
    gs = GridSpec(1, 2, width_ratios=[1.5, 1.0], wspace=0.22,
                  left=0.05, right=0.96, top=0.90, bottom=0.08)
    axT = fig.add_subplot(gs[0, 0])
    axM = fig.add_subplot(gs[0, 1])

    # ── 左: 树 ──
    for e in edges:
        p, c = by_id[e["parent_id"]], by_id[e["child_id"]]
        x0, x1 = node_x[p["node_id"]], node_x[c["node_id"]]
        if x0 is None or x1 is None:
            continue
        yy0, yy1 = node_y[p["node_id"]], node_y[c["node_id"]]
        axT.plot([x0, x1, x1], [yy0, yy0, yy1], color="#4a4a4a", lw=1.0,
                 solid_capstyle="round", zorder=2)
    for n in nodes:
        if not n.get("location") or n["age"] is None:
            continue
        x, yy = node_x[n["node_id"]], node_y[n["node_id"]]
        if n["is_tip"]:
            axT.scatter(x, yy, s=52, c=palette[n["location"]],
                        edgecolors="#2c3e50", linewidths=0.7, zorder=4)
            if n["name"]:
                lbl = rmap.get(n["name"], n["name"])
                lbl = lbl if len(lbl) <= 24 else lbl[:23] + "…"
                axT.annotate(lbl, (x, yy), xytext=(6, -2),
                             textcoords="offset points", fontsize=6.5,
                             color="#444", zorder=5)
        else:
            pr = n.get("location_prob") or 0.0
            axT.scatter(x, yy, s=24, c=palette[n["location"]],
                        alpha=0.35 + 0.6 * pr, edgecolors="none", zorder=3)
    axT.set_yticks([])
    axT.set_xlabel("Date", fontsize=11)
    axT.set_title("Tree", fontsize=15, fontweight="bold", pad=10)
    axT.spines["top"].set_visible(False)
    axT.spines["right"].set_visible(False)
    axT.spines["left"].set_visible(False)
    axT.tick_params(axis="x", labelsize=9)
    import matplotlib.ticker as mticker
    axT.xaxis.set_major_locator(mticker.MaxNLocator(integer=True, nbins=8))
    axT.xaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: f"{int(round(v))}"))

    # ── 右: 地图（弧线按终点地点配色 + 置信线型）──
    pw = pd.read_csv(pathways_csv) if (pathways_csv and os.path.exists(pathways_csv)) \
        else pd.DataFrame(columns=["branch_id", "start_location", "end_location",
                                   "start_age", "end_age", "confidence"])
    if boundaries_geojson:
        _draw_boundaries(axM, boundaries_geojson)
    drawn_conf = set()
    for _, r in pw.iterrows():
        if coords.get(r["start_location"]) is None or coords.get(r["end_location"]) is None:
            continue
        la0, lo0 = coords[r["start_location"]]
        la1, lo1 = coords[r["end_location"]]
        bx, by = bezier_arc_points(lo0, la0, lo1, la1, curvature=curvature)
        conf = r.get("confidence", "unclassified")
        ls, lw = conf_style.get(conf, conf_style["unclassified"])
        drawn_conf.add(conf)
        axM.plot(bx, by, color=palette.get(r["end_location"], "#555"),
                 linestyle=ls, linewidth=lw, alpha=0.9, zorder=3)
        axM.annotate("", xy=(bx[-1], by[-1]), xytext=(bx[-3], by[-3]),
                     arrowprops=dict(arrowstyle="-|>",
                                     color=palette.get(r["end_location"], "#555"),
                                     lw=0, mutation_scale=14), zorder=4)
    for loc in sorted(palette):
        if coords.get(loc) is None:
            continue
        lat, lon = coords[loc]
        axM.scatter(lon, lat, s=80, c=palette[loc], edgecolors="#2c3e50",
                    linewidths=1.0, zorder=5)
        axM.annotate(loc, (lon, lat), textcoords="offset points", xytext=(5, 5),
                     fontsize=9, zorder=6)
    axM.set_title("Map", fontsize=15, fontweight="bold", pad=10)
    axM.set_xlabel("Longitude")
    axM.set_ylabel("Latitude")
    axM.set_aspect("equal", adjustable="datalim")
    axM.autoscale()
    axM.margins(0.12)

    # ── 共用图例: Region（两panel配色）居中; Pathway Type 在地图内 ──
    region_handles = [Patch(facecolor=palette[l], edgecolor="#2c3e50", label=l)
                      for l in sorted(palette)]
    fig.legend(handles=region_handles, loc="center", bbox_to_anchor=(0.615, 0.28),
               title="Region", fontsize=8, title_fontsize=9, frameon=True)
    type_handles = [Line2D([], [], color="#555", ls=ls, lw=lw_, label=lbl)
                    for lbl, (ls, lw_) in
                    [("Direct (likely)", conf_style["direct"]),
                     ("Indirect (missing intermediates)", conf_style["indirect"]),
                     ("Distant / Import", conf_style["distant_import"])]]
    if drawn_conf:
        axM.legend(handles=type_handles, loc="lower left", fontsize=8,
                   title="Pathway Type", title_fontsize=9, frameon=True)

    os.makedirs(os.path.dirname(os.path.abspath(out_path)), exist_ok=True)
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
    root_, _ = os.path.splitext(out_path)
    fig.savefig(root_ + ".pdf", bbox_inches="tight")
    plt.close(fig)
    log.emit(f"树-图联动面板: {out_path} (+pdf), {len(nodes)} 节点 / {len(pw)} 事件")
    return out_path


# ─────────────────────────────────────────────────────────────────────
# CLI
# ─────────────────────────────────────────────────────────────────────

def main():
    import argparse
    ap = argparse.ArgumentParser(
        description="发表级传播图 + 事件动画（底图合规：默认纯矢量，无在线瓦片）")
    ap.add_argument("--branch-table", default=None, help="transmission_paths 的 branch_table.csv")
    ap.add_argument("--mcc-tree", default=None, help="或直接给 MCC 树（现场解析）")
    ap.add_argument("--metadata", default=None, help="name,location CSV（--mcc-tree 时用）")
    ap.add_argument("--coords", default=None, help="location,lat,lon CSV（缺省走离线 geo_resolver）")
    ap.add_argument("--boundaries", default=None,
                    help="国界/省界 GeoJSON（须自备审图号来源；缺省不画边界）")
    ap.add_argument("-o", "--out", default="transmission_map.png", help="静态图输出(PNG)")
    ap.add_argument("--gif", default=None, help="动画 GIF 输出路径（给了才生成）")
    ap.add_argument("--gif-frames", type=int, default=48)
    ap.add_argument("--trail", type=int, default=-1,
                    help="-1/负数=累积模式(默认,已揭示路线永久保留); 正整数=残影模式(只高亮最近 N 条)")
    ap.add_argument("--pacing-lambda", type=float, default=0.5,
                    help="blended pacing: 1=纯时间 0=纯秩")
    ap.add_argument("--curvature", type=float, default=0.2)
    ap.add_argument("--no-stream", action="store_true", help="关闭弧线渐显(用统一透明度)")
    ap.add_argument("--progress", default=None,
                    help="扩散进展多边形叠加图输出路径 (seraphim 式凸包嵌套)")
    ap.add_argument("--snapshots", default=None,
                    help="时间快照多面板输出路径 (GIF 的静态发表替代)")
    ap.add_argument("--snap-panels", type=int, default=4)
    ap.add_argument("--tree-map", default=None,
                    help="树-图联动双面板输出路径 (phymapr Tree|Map 图版替代; 需 --mcc-tree+--metadata)")
    ap.add_argument("--rename", default=None,
                    help="tip 短名映射文件 (旧ID<TAB/逗号/⇥>新ID, # 注释) — 树图标签美化")
    args = ap.parse_args()

    if args.branch_table:
        bt = pd.read_csv(args.branch_table)
    elif args.mcc_tree:
        from utils.transmission_paths import extract_tree_events
        ev = extract_tree_events(args.mcc_tree, args.metadata)
        bt = pd.DataFrame(ev["edges"])
    else:
        raise SystemExit("需要 --branch-table 或 --mcc-tree")

    log = LogCollector()
    locs = sorted(set(bt["start_location"].dropna()) | set(bt["end_location"].dropna()))
    coords = resolve_coordinates(locs, args.coords, log=log)
    missing = [l for l in locs if coords.get(l) is None]
    if missing:
        # 2026-09-17 修: 缺坐标地点相关事件跳行+警告 (原实现硬退出;
        # 真实 PSTVd MCC 含 Unknown 地点节点时整张图出不来)
        log.warning(f"地点缺坐标, 相关事件跳过: {missing}")
        bt = bt[(~bt["start_location"].isin(missing)) & (~bt["end_location"].isin(missing))]
        if bt.empty:
            raise SystemExit(f"全部事件都涉及缺坐标地点: {missing}")

    draw_static_map(bt, coords, args.out, boundaries_geojson=args.boundaries,
                    curvature=args.curvature, stream=not args.no_stream, log=log)
    if args.progress:
        draw_progress_polys(bt, coords, args.progress, boundaries_geojson=args.boundaries,
                            log=log)
    if args.snapshots:
        draw_snapshot_panels(bt, coords, args.snapshots, n_panels=args.snap_panels,
                             boundaries_geojson=args.boundaries, log=log)
    if args.tree_map:
        if not (args.mcc_tree and args.metadata):
            raise SystemExit("--tree-map 需要 --mcc-tree + --metadata（树面板要吃树和日期）")
        # 置信分级在同目录 pathways.csv 里（transmission_paths 的产物布局）
        pw_csv = (args.branch_table.replace("branch_table.csv", "pathways.csv")
                  if args.branch_table else None)
        draw_tree_map(args.mcc_tree, args.metadata, coords, args.tree_map,
                      pathways_csv=pw_csv, rename_map=args.rename,
                      boundaries_geojson=args.boundaries, log=log)
    if args.gif:
        draw_animation(bt, coords, args.gif, n_frames=args.gif_frames,
                       trail=args.trail,
                       pacing_lambda=args.pacing_lambda, boundaries_geojson=args.boundaries,
                       curvature=args.curvature, log=log)
    print("TRANSMISSION_MAP_DONE")


if __name__ == "__main__":
    main()
