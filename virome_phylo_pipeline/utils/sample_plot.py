#!/usr/bin/env python3
"""sample_plot.py — VirPhyKit SamplePlot 复刻 (采样分布可视化)

复刻 VirPhyKit SamplePlot 模块的两个功能:
  1. 采样时间分布图 (generate_plot.R): Year/Total/地点 → 时间线折线 + 地点分布点线
  2. 采样地理地图 (generate_map.R): region/Longitude/Latitude → 地图 + 采样点

纯 Python (matplotlib + cartopy), 无 R 依赖。

口径说明 (2026-09-16 上游核对修复)
----------------------------------
上游 generate_map.R 要求输入表**已经存在** `Longitude/Latitude/region` 三列 ——
经纬度与大区是在上游流程的前置步骤里预先算好的, 不是从 location 现场猜。

初版本管线改成从 `location` 字符串硬匹配一张写死的 12 个省名表, 后果实测:
用上游自带示例 `Example/TreeTime-RTT/H3N2`(476 样本, 27 个唯一 location) 运行,
**476/476 全部落进 'Other' 桶**, 时间分布图地点面板只剩一行、地图只剩一个点,
且全程无任何告警。对全球数据等于该图报废。

现改为:
  1. 优先用输入表已有的 `Longitude/Latitude`(+可选 `region`) 列 —— 与上游接口对齐;
  2. 无坐标列时, 调 `utils.geo_resolver` 分层解析地名 (中国全行政区划 + 常见国家 +
     逗号段拆分 + 持久缓存 + 可选在线兜底), 取代写死的省名表;
  3. 解析失败/落 Other 占比超过阈值时**显式告警**(不再静默),
     并打印未识别地名清单, 便于人工补 geo 缓存。

用法:
  python -m utils.sample_plot --metadata X.csv \
      [--date-col date] [--location-col location] [--outdir out]
"""
import argparse
import csv
import os
import re
import sys
from collections import Counter, defaultdict

# 2026-10-16: 图件双格式 (PNG+PDF 同主名); 统一由 plot_io.save_both 收口。
# 双路径 import: `-m utils.sample_plot` 与 `python utils/sample_plot.py` 均可用。
try:
    from utils.plot_io import save_both
except ImportError:                              # pragma: no cover
    from plot_io import save_both

# 未识别地名占比超过该阈值即告警 (0.2 = 20%)
UNRESOLVED_WARN_FRACTION = 0.2


def _get_resolver():
    """惰性获取 geo_resolver (允许 import 失败时退化, 但必须告警)。"""
    try:
        from utils.geo_resolver import get_resolver
        return get_resolver()
    except Exception as e:                       # pragma: no cover
        print(f"[sample_plot] 警告: geo_resolver 不可用 ({e}), "
              f"地图将无法解析地名", file=sys.stderr)
        return None


def load_coord_columns(rows):
    """若输入表自带 Longitude/Latitude (+region) 列, 直接采用 (与上游接口对齐)。

    返回 {row_index: (lat, lon, region)} 或 None。
    """
    if not rows:
        return None
    cols = {c.lower(): c for c in rows[0].keys()}
    lat_c = cols.get('latitude')
    lon_c = cols.get('longitude')
    if not (lat_c and lon_c):
        return None
    reg_c = cols.get('region')
    out = {}
    for i, r in enumerate(rows):
        try:
            lat = float(str(r.get(lat_c, '')).strip())
            lon = float(str(r.get(lon_c, '')).strip())
        except (TypeError, ValueError):
            continue
        if not (-90 <= lat <= 90 and -180 <= lon <= 180):
            continue                      # 坐标越界 → 丢弃, 不伪造
        out[i] = (lat, lon, (r.get(reg_c) or '').strip() if reg_c else '')
    return out or None


def _use_cjk_font_if_available():
    """地点标签常是中文 (metadata 的 location/region 列) —— 探测系统 CJK 字体并
    插到 sans-serif 首位; 没有则保持默认。

    2026-10-12 实测: 不探测时中文地名全出方框 (DejaVu Sans 无 CJK 字形) ——
    A/B 两面板的图例与行标签全部不可读。必须在 apply_sci_style() **之后**调用
    (它会重置 rcParams)。与 utils/splitstree_viz.py::_setup_font 同一探测清单。
    """
    try:
        import matplotlib
        from matplotlib import font_manager
        names = {f.name for f in font_manager.fontManager.ttflist}
        for cand in ("Microsoft YaHei", "SimHei", "Noto Sans CJK SC",
                     "Source Han Sans SC", "WenQuanYi Zen Hei", "PingFang SC"):
            if cand in names:
                cur = list(matplotlib.rcParams["font.sans-serif"])
                if cand not in cur:
                    matplotlib.rcParams["font.sans-serif"] = [cand] + cur
                break
        matplotlib.rcParams["axes.unicode_minus"] = False
    except Exception:  # noqa: BLE001 — 字体探测失败不致命 (退回默认字体)
        pass


def extract_year(date_str):
    m = re.search(r'(\d{4})', str(date_str or ''))
    return int(m.group(1)) if m else None


def _canonical_label(loc_str):
    """把原始地名压成显示用标签: 取逗号段的最后一段 (最具体), 保留原拼写。"""
    parts = [p.strip() for p in re.split(r'[,;]', str(loc_str)) if p.strip()]
    if not parts:
        return str(loc_str).strip()
    # 末段若是 'China'/'中国' 这类宽泛词, 往前取一段
    for p in reversed(parts):
        if p.lower() not in ('china', '中国', 'unknown', 'na'):
            return p
    return parts[0]


# ── 国家级/无信息地名（不能落到具体采样点）──────────────────────────
# 2026-10-12（审稿人视角修复）：metadata 里有 1 条 location 只写 'China' 且无日期。
# 旧版会把它解析成"中国质心"并在地图上打一个气泡 —— 那不是一个采样点，
# 而且它只出现在地图面板、不出现于时间面板（无日期被跳过），
# 造成 A/B 两地点的集合与配色双双不一致。现统一剔除并记数。
_COUNTRY_LEVEL = {'china', '中国', 'unknown', 'na', 'n/a', ''}


def _is_specific_label(label: str) -> bool:
    return str(label).strip().lower() not in _COUNTRY_LEVEL


def _row_region(row) -> str:
    """取该行的 region 列值 (大小写不敏感); 无该列返回 ""。"""
    for k, v in row.items():
        if str(k).strip().lower() == 'region':
            return str(v or '').strip()
    return ''


def resolve_label(row, location_col='location', region='', resolver=None):
    """**统一**地点标签解析 —— 时间面板与地图面板共用同一入口。

    2026-10-12 修 (P1-5): 旧实现两面板各走一路 —— 时间面板只用 geo_resolver
    解析 location (完全无视表内 region 列), 地图面板优先 region 列; 同一份数据
    会得出两套标签 (实测: 时间面板 1 个 'Other', 地图 8 个真实地点), 位置面板
    的"同名同色"承诺在一开始就不成立。
    优先级 (与上游 generate_map.R 的 region 口径对齐):
      1) 输入表 region 列非空 → 直接采用 (上游预置口径, 原样信任);
      2) geo_resolver 能解析 location 字符串 → 取最具体的一段 (_canonical_label);
      3) 均无 → None (调用方统一落 'Other')。
    第 2 步解析出的**只到国家级**标签(如 'China')同样返回 None —— 它不是可定位
    的采样点 (地图侧另有 dropped_imprecise 计数, 不静默丢数据)。
    """
    name = str(region or '').strip()
    if name:
        return name
    loc = (row.get(location_col, '') or '').strip()
    if not loc or resolver is None:
        return None
    coord, _src = resolver.resolve(loc)
    if coord is None:
        return None
    name = _canonical_label(loc)
    return name if _is_specific_label(name) else None


def loc_color_map(names):
    """地点名集合 → {name: hex}。

    **跨面板稳定的配色映射**：按地点名排序后取 OKABE_ITO。
    旧版三个面板各自按"枚举顺序"取色（A 用 dict 插入序、B 用 sorted、
    V 用 R 侧独立顺序），同一地点在 A 是绿色、在 B 却是橙色 —— 多面板图里
    这是硬伤。统一走本函数后 A/B/V 颜色一一对应。
    """
    from utils.sci_palette import colors
    ordered = sorted({str(n) for n in names if str(n).strip()})
    pal = colors(max(len(ordered), 1))
    return {n: pal[i % len(pal)] for i, n in enumerate(ordered)}


def _size_legend(ax, counts, scale, *, title='isolates', loc='lower left',
                 bbox=None, fontsize=5.6, ncol=1):
    """在 ax 上追加"气泡大小 = 样本数"的图例（审稿人必问项）。

    旧版两个气泡图都只解释了颜色、没解释大小 —— 读者无从判断圆的大小对应多少条序列。
    标记尺寸换算：matplotlib 的 `s` 是 pt²，故 markersize = sqrt(s)。

    两处必须处理（2026-10-12 实测）：
      · 参考级别取 **3 档**（小/中/大）而非 4 档 —— 档位越多图例越高，容易顶到数据；
      · `handleheight/handlelength` 要按**最大色块**反推，否则大圆会**溢出图例框**，
        看起来像凭空多了一个数据点（vmax=69 时色块直径 ≈ 37 pt，实测溢出）。
    """
    from matplotlib.lines import Line2D
    vmax = int(max(counts)) if counts else 1
    vmin = int(min(counts)) if counts else 1
    lv = sorted({max(1, vmin), max(1, int(round(vmax * 0.5))), max(1, vmax)})
    dmax = max((max(v * scale, 12) ** 0.5) for v in lv)
    handles = [Line2D([], [], marker='o', linestyle='none',
                      markersize=max(2.2, max(v * scale, 12) ** 0.5),
                      markerfacecolor='#BBBBBB', markeredgecolor='#555555',
                      markeredgewidth=0.4, label=str(v)) for v in lv]
    hh = max(0.7, dmax / float(fontsize))
    # ⚠ 必须挂到 **figure** 级（fig.legend + bbox_transform=ax.transAxes），
    #   不能走 ax.legend()：后者会把同一个 ax 上**先加的地点颜色图例删掉**
    #   （实测 A/B 的"颜色=地点"图例就是这样消失的）。figure 级图例不碰
    #   ax.legend_，两条图例可以共存，且仍按 ax 的坐标分数定位。
    leg = ax.figure.legend(
        handles=handles, loc=loc, fontsize=fontsize, frameon=True,
        title=title, labelspacing=0.62, borderpad=0.4,
        handletextpad=0.7, bbox_to_anchor=bbox,
        bbox_transform=ax.transAxes,
        handleheight=hh, handlelength=hh, ncol=ncol)
    if title:
        leg.get_title().set_fontsize(fontsize)
    return leg


def _marker_bboxes(ax, fig, extra=None):
    """收集 ax 上所有散点标记的**显示坐标** bbox（含 inset 等额外避让框）。

    用于图例自动避让：`Collection.get_window_extent()` 对 scatter 不可靠，
    故从 offsets（数据坐标）经 transData 换到显示坐标，半径由标记面积 s 反推
    （matplotlib 里 s 是 pt²，圆标记直径 ≈ sqrt(s) pt）。
    """
    from matplotlib.transforms import Bbox
    boxes = list(extra or [])
    try:
        rend = fig.canvas.get_renderer()
    except Exception:
        fig.canvas.draw()
        rend = getattr(fig.canvas, "renderer", None)
    if rend is None:
        return boxes
    for c in ax.collections:
        # ⚠ cartopy 的 FeatureArtist（COASTLINE/BORDERS）也在 ax.collections 里，
        #   但它没有 get_sizes()/get_offsets() —— 必须逐个 try，否则整段崩。
        try:
            off = c.get_offsets()
            sizes = list(c.get_sizes())
        except Exception:
            continue
        if off is None or len(off) == 0:
            continue
        r_pt = 0.5 * ((max(sizes) if sizes else 12.0) ** 0.5)
        r_px = r_pt * fig.dpi / 72.0
        for x, y in ax.transData.transform(off):
            boxes.append(Bbox([[x - r_px, y - r_px], [x + r_px, y + r_px]]))
    return boxes


def _place_size_legend(ax, counts, scale, fig, candidates, *, title='isolates',
                       fontsize=5.6, extra_boxes=None, ncol=1, outside=None):
    """把"气泡大小 = 样本数"图例放到**第一个既不出框、又不压气泡**的候选位置。

    尺寸图例的最大色块是按真实比例画的（vmax=69 → 直径 ≈ 37 pt ≈ 0.51 in），
    一旦落在有数据的角落就会把气泡整个盖住 —— 2026-10-12 实测：右上角的
    '69' 色块大到**溢出图例框**，压住地图上最大的采样点。故不能写死角落，
    必须按当前数据逐个候选位置试算。

    ⚠ 判据必须**同时**满足两条（只查气泡会出事）：
      (a) 不与任何气泡相交；
      (b) 图例整体落在坐标区内 —— 否则会越出坐标区压住刻度标签、甚至被裁掉
          （实测：B 面板放 (0, 0.16) 时 '36' 色块掉到 x 轴标签上并被裁一半）。

    矮而密的面板（如 B）可能**四个角全被气泡占满** → 此时退到 `outside`
    给出的"坐标区外"候选（只要求不压气泡，允许出框；`bbox_inches='tight'`
    会把画布撑大，不会被裁）。
    """
    boxes = _marker_bboxes(ax, fig, extra=extra_boxes)
    ax_bb = ax.get_window_extent()
    _dbg = bool(__import__("os").environ.get("SP_DEBUG"))
    if _dbg:
        print("[sizelegend] n_marker_boxes=%d  ax_bb=%s  inside=%s  outside=%s"
              % (len(boxes), tuple(round(v, 1) for v in ax_bb.extents),
                 [c[0] for c in candidates], [c[0] for c in (outside or [])]),
              file=__import__("sys").stderr)
    chosen = None
    plan = [(c, True) for c in candidates] + [(c, False) for c in (outside or [])]
    for (loc, bbox), need_inside in plan:
        leg = _size_legend(ax, counts, scale, title=title, loc=loc, bbox=bbox,
                           fontsize=fontsize, ncol=ncol)
        try:
            fig.canvas.draw()
            lb = leg.get_window_extent(fig.canvas.get_renderer())
            inside = (lb.x0 >= ax_bb.x0 - 1 and lb.x1 <= ax_bb.x1 + 1 and
                      lb.y0 >= ax_bb.y0 - 1 and lb.y1 <= ax_bb.y1 + 1)
            nov = sum(1 for b in boxes if lb.overlaps(b))
            if _dbg:
                print("   try %-14s inside=%-5s overlap=%d  lb=%s"
                      % (loc, inside, nov,
                         tuple(round(v, 1) for v in lb.extents)),
                      file=__import__("sys").stderr)
            if nov == 0 and (inside or not need_inside):
                chosen = (loc, bbox)
                break
        except Exception as e:
            if _dbg:
                print("   try %-14s EXC %s" % (loc, e), file=__import__("sys").stderr)
            chosen = (loc, bbox)
            break
        leg.remove()
    if chosen is None:
        # 候选位置全被占：宁可重叠也不能没有图例，退回最后一个候选
        chosen = (outside or candidates)[-1]
        _size_legend(ax, counts, scale, title=title, loc=chosen[0],
                     bbox=chosen[1], fontsize=fontsize, ncol=ncol)
    return chosen


def build_temporal_data(rows, date_col, location_col, log=None, host_col=None):
    """从 metadata 生成 Year/Total/地点 宽表。

    地点标签统一走 resolve_label() (region 列 → geo_resolver → 'Other'),
    与地图面板同源 —— 两面板的标签集合/配色严格一致 (2026-10-12 修 P1-5)。

    2026-10-13: 加 `host_col` → 额外统计 (year, location, host) 三元计数,
    供气泡内切饼用 (用户「图 B 把物种也加上」)。

    返回 (years, locs, year_loc, unresolved) ；给了 host_col 时返回 5 元组,
    第 5 项 = year_loc_host[year][loc][host]。
    """
    resolver = _get_resolver()
    coord_cols = load_coord_columns(rows) or {}
    year_loc = defaultdict(Counter)
    year_loc_host = defaultdict(lambda: defaultdict(Counter))
    unresolved = Counter()          # 未解析成功的原始地名 → 出现次数

    for i, r in enumerate(rows):
        year = extract_year(r.get(date_col, ''))
        loc = (r.get(location_col, '') or '').strip()
        reg = coord_cols.get(i, (None, None, ''))[2] if i in coord_cols \
            else _row_region(r)
        if not (year and (loc or reg)):
            continue
        label = resolve_label(r, location_col, region=reg, resolver=resolver)
        if label is None:
            unresolved[loc or reg] += 1
            label = 'Other'
        year_loc[year][label] += 1
        if host_col:
            # ⚠ 宿主是**双名法**，可能写成 "Lycium barbarum, Nicotiana ..."
            #   → 取**首段**（主宿主），与 stage_plots 的 _load_host_map 同口径。
            _h = (r.get(host_col) or '').split(',')[0].strip()
            if _h:
                year_loc_host[year][label][_h] += 1

    years = sorted(year_loc)
    locs = sorted({l for c in year_loc.values() for l in c})
    if log is not None and unresolved:
        n_unres = sum(unresolved.values())
        n_tot = sum(sum(c.values()) for c in year_loc.values())
        frac = n_unres / n_tot if n_tot else 0.0
        msg = (f"[sample_plot] 地名解析: {n_unres}/{n_tot} "
               f"({frac*100:.1f}%) 未识别, 已归入 'Other'")
        if frac >= UNRESOLVED_WARN_FRACTION:
            log.warning(msg + " —— 占比过高, 地点面板/地图可能失真。"
                              "未识别清单: " + ", ".join(
                                  f"{k}({v})" for k, v in unresolved.most_common(10)))
        else:
            log.emit(msg)
    if host_col:
        return years, locs, year_loc, unresolved, year_loc_host
    return years, locs, year_loc, unresolved


def plot_sample_temporal(rows, date_col, location_col, out_path, log=None,
                         host_col=None):
    """时间线折线 + 地点分布气泡 (两面板)。

    2026-10-13（用户「图 B 把物种也加上，饼图？图例在空白处」）：
      · `host_col` 给了 → 每个气泡内**切饼**表示该(年,地点)的宿主构成；
      · 两组图例统一**移到面板内左上空白区**（原在坐标区下方，占地方）。
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from utils.sci_palette import (apply_sci_style, colors, despine,
                                   sci_legend, COL_W_2, OKABE_ITO)
    apply_sci_style(base_size=8)
    _use_cjk_font_if_available()   # 中文地名防方框 (P1-5 验证需要)

    _bt = build_temporal_data(rows, date_col, location_col, log=log,
                              host_col=host_col)
    if host_col:
        years, locs, year_loc, _, year_loc_host = _bt
    else:
        years, locs, year_loc, _ = _bt
        year_loc_host = {}
    totals = {y: sum(year_loc[y].values()) for y in years}

    # ⚠⚠ 2026-10-13：布局改用 **3 行 GridSpec**（折线 / 图例带 / 气泡）。
    #   图例从"坐标区下方"挪到"两面板之间"后，若仍用 2 行 + tight_layout + hspace，
    #   会触发 `UserWarning: This figure includes Axes that are not compatible with
    #   tight_layout`（因为 `fig.legend(bbox_transform=fig.transFigure)` 让
    #   tight_layout 无法正确估算），且几何不可控。3 行 GridSpec 把"图例带"变成
    #   一行真实的、高度可配的区域，几何全由 `height_ratios` 决定，无需 tight_layout。
    fig = plt.figure(figsize=(COL_W_2, COL_W_2 * 1.05))
    # ⚠ 图例带（中间行）高度要够：放不下会与上排的年份刻度、下排的标题**文字相撞**
    #   （实测 0.62 时 2018/2020 等年份标签直接压在图例标题上 —— 零报错）。
    _gs = fig.add_gridspec(3, 1, height_ratios=[3, 1.05, 2], hspace=0.06,
                           left=0.145, right=0.985, top=0.975, bottom=0.12)
    ax1 = fig.add_subplot(_gs[0])
    _gap_ax = fig.add_subplot(_gs[1])           # 只作图例容器：全隐藏
    ax2 = fig.add_subplot(_gs[2], sharex=ax1)
    _gap_ax.set_axis_off()
    _gap_ax.set_xlim(0, 1)
    _gap_ax.set_ylim(0, 1)
    # 空/无有效记录守卫 (历史坑: location 列全空或列名不符时 min/max 直接 ValueError)
    if not years or not totals:
        print("  [!] 无有效 年份×地点 记录, 跳过采样分布图")
        plt.close(fig)
        return
    # 上: 总样本时间线 (复刻 p1) —— 单序列用 OKABE 蓝, 可读性最好
    ax1.plot(years, [totals[y] for y in years],
             color=OKABE_ITO[4], lw=1.2, marker='o', ms=2.6, zorder=3)
    ax1.set_ylabel('No. of isolates')
    # 2026-10-12: 两端各留 0.6 年 —— 旧版 xlim 恰好卡在 min/max 年上, 首末两个
    # 数据点被坐标区边缘切掉一半 (实测 2024 那个点只画出一半)。
    ax1.set_xlim(min(years) - 0.6, max(years) + 0.6)
    ax1.set_ylim(0, max(totals.values()) * 1.12 + 1)
    ax1.grid(axis='y', color='#DDDDDD', lw=0.5)
    ax1.grid(axis='x', visible=False)
    despine(ax1)
    ax1.tick_params(axis='x', labelrotation=30)

    # 下: 各地点分布点线 —— 统一 OKABE_ITO 分类色
    # 2026-10-12: 颜色改走 loc_color_map()（按地点名排序），与地图面板 A 一一对应。
    cmap = loc_color_map(locs)
    # 气泡尺度：按最大计数自适应，保证最大圆不超过一行行高（旧版固定 v*26，
    # 最大地点 n=69 → 直径 42 pt ≈ 两倍行高，圆会盖住相邻行的标签）。
    _vmax = max([year_loc[y].get(l, 0) for y in years for l in locs] or [1])
    _scale = min(26.0, 380.0 / max(_vmax, 1))
    _all_counts = []

    # 宿主色板（给了 host_col 才用）—— 与 ML 树 / NeighborNet 的宿主着色同源，
    # 保证整组主图「橙 = L. barbarum / 蓝 = L. chinense / 绿 = L. ruthenicum」一致。
    _hosts = sorted({h for y in year_loc_host for l in year_loc_host[y]
                     for h in year_loc_host[y][l]}) if year_loc_host else []
    _host_cmap = {h: OKABE_ITO[i % len(OKABE_ITO)] for i, h in enumerate(_hosts)}
    _pie_specs = []      # 多宿主气泡的饼规格，等 tight_layout 之后再落多边形

    if not year_loc_host:
        # ── 旧路径：纯地区着色圆 ──
        for i, loc in enumerate(locs):
            xs = [y for y in years if year_loc[y].get(loc, 0) > 0]
            ys = [year_loc[y].get(loc, 0) for y in xs]
            _all_counts.extend(ys)
            ax2.scatter(xs, [i] * len(xs), s=[max(v * _scale, 12) for v in ys],
                        color=cmap[loc], alpha=0.85,
                        edgecolor='white', linewidth=0.4, zorder=3)
            ax2.plot(xs, [i] * len(xs), color=cmap[loc], lw=0.5, alpha=0.7, zorder=2)
    else:
        # ── 新路径：气泡内切饼表示宿主构成（用户 2026-10-13 要求）──
        # ⚠⚠ 2026-10-13 实测**两个致命坑**（都零报错）：
        #   ① `Wedge(center, r, ...)` / `Polygon` 的半径是**数据坐标**，而本图 x（年）
        #      与 y（行）两轴尺度完全不同 → 直接传"点半径"会把扇形画成铺满面板的
        #      巨大色块；且若在 `tight_layout()` **之前**算换算比例，拿到的还是旧
        #      axis box → 量纲和几何全错（实测饼缩成小点、还被拉成椭圆）。
        #   ② `transData.transform` 返回**像素**（含 dpi），而 `s` / `r_pt` 是**点** →
        #      必须 ×(72/dpi) 才是同量纲。
        #   正解 = 分两趟：本趟只**登记**饼的规格（中心 + 半径点值 + 扇区占比），
        #      把真正的多边形放到 `tight_layout()` **之后**（几何已定）再生成。
        for i, loc in enumerate(locs):
            xs = [y for y in years if year_loc[y].get(loc, 0) > 0]
            for x in xs:
                v = year_loc[x].get(loc, 0)
                _all_counts.append(v)
                s_pt = max(v * _scale, 12)           # scatter 的 s（点²）
                r_pt = (s_pt ** 0.5) / 2.0           # 点半径
                comp = year_loc_host.get(x, {}).get(loc, {})
                if not comp:
                    ax2.scatter([x], [i], s=[s_pt], color=cmap[loc], alpha=0.85,
                                edgecolor='white', linewidth=0.4, zorder=3)
                    continue
                tot = sum(comp.values()) or 1
                if len(comp) == 1:                   # 单宿主 → 整圆（免扇形缝）
                    ax2.scatter([x], [i], s=[s_pt],
                                color=_host_cmap[next(iter(comp))], alpha=0.85,
                                edgecolor='white', linewidth=0.4, zorder=3)
                    continue
                _pie_specs.append((x, i, r_pt,
                                   [(h, n / tot) for h, n
                                    in sorted(comp.items(), key=lambda kv: -kv[1])]))
            # 2026-10-15（用户裁定）：Y 轴刻度已**逐行写出地区名** → 行连线不再用
            #   地点色，改**中性灰**，彻底消除「没有图例解释的颜色」。地点色带随之
            #   删掉（见下方图例带；无宿主表的那条路径仍然保留地点图例）。
            ax2.plot(xs, [i] * len(xs), color='#AAAAAA', lw=0.5, alpha=0.55, zorder=2)
    ax2.set_yticks(range(len(locs)))
    ax2.set_yticklabels(locs)
    ax2.set_xlabel('Year')
    # 气泡上下各留 0.6 行、左右各留 0.6 年, 否则首末行的圆会被坐标区切掉半个
    ax2.set_ylim(-0.6, len(locs) - 0.4)
    ax2.set_xlim(min(years) - 0.6, max(years) + 0.6)
    ax2.grid(axis='x', color='#DDDDDD', lw=0.5)
    ax2.grid(axis='y', visible=False)
    despine(ax2)
    ax2.tick_params(axis='x', labelrotation=30)

    # ── 图例（2026-10-13 用户：「图例在空白处，不要在下面占地方」）──
    # 旧版把地区图例挂在 `bbox_to_anchor=(0.5, -0.46)`（坐标区**下方**），既占版面
    # 又把画布撑高。改为压缩成**一条横向图例带**，塞进折线图与气泡图之间的空隙。
    #
    # ⚠⚠ 2026-10-13 实测：**不能塞进气泡图内部**。真实数据 9 地点 + 3 宿主 + 尺寸
    #   图例共 13 项，竖排在面板里会从顶贯到底、把银川/北京那几行的气泡**全部盖住**，
    #   而且框会溢到上面的折线图里（首版实测，且**零报错**）。用户看过后选"放两面板
    #   之间的空隙" → 只有这一种排法能在不遮数据、不占下方版面的前提下放下全部项。
    #
    # 做法：GridSpec 里给"图例带"留一行（见函数开头的 3 行 GridSpec），图例直接挂在
    #   那一行的 `_gap_ax` 上（用它自己的 transAxes，几何完全可控、无需 tight_layout）。
    from matplotlib.lines import Line2D
    # 2026-10-12（审稿人视角）：与 plot_sample_map 同理 —— 气泡 alpha=0.85，
    # 图例色块必须同透明度，否则同一地点在图例与气泡上呈两种色调。
    handles = [Line2D([], [], marker='o', linestyle='none', markersize=3.4,
                      markerfacecolor=cmap[loc], alpha=0.85,
                      markeredgecolor='white', markeredgewidth=0.4, label=loc)
               for loc in locs]
    fig.canvas.draw()

    # ── 饼图落盘（⚠ 必须在 tight_layout + canvas.draw 之后）──
    # 此时 axis box 已定，`transData` 反映真实几何 → 换算比例才对。
    # `transData` 给**像素**，`r_pt` 是**点** → ×(72/dpi) 才是同量纲。
    if _pie_specs:
        from matplotlib.patches import Polygon
        import math as _m
        _px2pt = 72.0 / max(fig.dpi, 1e-9)
        _xlo, _xhi = ax2.get_xlim()
        _ylo, _yhi = ax2.get_ylim()
        (x0p, y0p), (x1p, y1p) = ax2.transData.transform(
            [(_xlo, _ylo), (_xhi, _yhi)])
        _sx = abs(x1p - x0p) / max(abs(_xhi - _xlo), 1e-9) * _px2pt   # 点 / 年
        _sy = abs(y1p - y0p) / max(abs(_yhi - _ylo), 1e-9) * _px2pt   # 点 / 行
        for _x, _i, _r_pt, _fracs in _pie_specs:
            _rx = _r_pt / max(_sx, 1e-9)      # 水平数据半径（年）
            _ry = _r_pt / max(_sy, 1e-9)      # 垂直数据半径（行）
            _a0 = 90.0
            for _h, _f in _fracs:
                _a1 = _a0 - 360.0 * _f
                _k = max(3, int(abs(_a1 - _a0) / 3.0) + 2)
                _ang = [_m.radians(_a0 + (_a1 - _a0) * t / (_k - 1))
                        for t in range(_k)]
                _pts = [(_x, _i)] + [(_x + _rx * _m.cos(a), _i + _ry * _m.sin(a))
                                     for a in _ang]
                ax2.add_patch(Polygon(_pts, closed=True,
                                      facecolor=_host_cmap[_h], alpha=0.85,
                                      edgecolor='white', linewidth=0.4,
                                      zorder=3))
                _a0 = _a1

    # ── 图例带：挂在 `_gap_ax`（3 行 GridSpec 的中间那行）上 ──
    # 位置全用 `_gap_ax` 的轴分数，几何由 height_ratios 决定 → 无需 tight_layout，
    # 也不会"图例溢出到面板里"。
    # 排法（2026-10-15 起）：**只有**宿主图例 + 气泡尺度说明。
    #   ⚠⚠ 地点图例**只在「没有宿主表」时才画** —— 那条路径下气泡是**按地点着色**
    #     的（见 `if not year_loc_host:` 分支），删掉就没法读图；有宿主表时气泡改按
    #     宿主着色、行连线改中性灰 → 地点色带**纯冗余**（Y 轴刻度已逐行写出地区名，
    #     面板 A 的环形树外圈也直接标注了地区名）→ 用户 2026-10-15 裁定删除。
    #   ⚠ 宿主图例**不能**横排 ncol=3 —— 三个双名法名横排会顶到最右、把
    #     `L. ruthenicum` 挤出画布（实测被裁掉且零报错）。
    # ⚠⚠ 2026-10-15 实测（PyMuPDF 逐 span bbox 对拍，非目视）：图例带**竖直居中
    #   (y=0.5)** 时，4 行的宿主图例块顶边（y≈166.2 pt）会**压进上排年份刻度**的
    #   45° 旋转标签 —— '2016' 的 bbox 一直到 y≈171.2 pt，二者实测重叠 **63 pt²**
    #   （**零报错**，只能量像素）。根因：图例带那一行只有 ~49 pt 高，而 ax1 的旋转
    #   年份标签比行间距 (hspace=0.06) 还高，必然越界到图例带顶部。
    #   正解 = 把整条图例带**下移**（锚点 y 0.5 → 0.28），**不动任何坐标区尺寸**：
    #   下移后距 '2016' 底 ~6 pt、距下排坐标区顶 ~11 pt，上下两清。
    #   ⚠ 改这个值必须**同时**改下面气泡尺度说明 `_gap_ax.text` 的 y，否则两者错位。
    _BAND_Y = 0.28
    _fs_leg = 4.9
    _gbb = None
    _prev_r = 0.0
    if not year_loc_host:
        _leg_loc = _gap_ax.legend(handles=handles, loc='center left',
                                  bbox_to_anchor=(0.0, _BAND_Y),
                                  bbox_transform=_gap_ax.transAxes,
                                  ncol=4, fontsize=_fs_leg,
                                  title='location', frameon=False, handletextpad=0.3,
                                  columnspacing=0.85, labelspacing=0.32,
                                  borderpad=0.0, borderaxespad=0.0)
        _leg_loc.get_title().set_fontsize(_fs_leg)
        _leg_loc.get_title().set_fontweight('bold')
        # ⚠ 不要再 `_gap_ax.add_artist(_leg_loc)`！`ax.legend()` 已把图例挂在
        #   `ax.legend_` 上、由 Axes 正常绘制；再 `add_artist` 会让**同一条图例被
        #   绘制两次**（2026-10-15 实测：PDF 里每个图例标签都出现两遍，bbox 完全
        #   重合 → 文字重叠对拍时表现为一堆「自重叠」，把真正的碰撞淹没掉）。
        _prev_r = 0.60
        try:
            fig.canvas.draw()
            _gbb = _gap_ax.get_window_extent()
            _lb = _leg_loc.get_window_extent(fig.canvas.get_renderer())
            _loc_w = (_lb.x1 - _gbb.x0) / max(_gbb.width, 1e-9)
            _prev_r = min(_loc_w + 0.015, 0.72)
        except Exception:
            pass

    if year_loc_host:
        _h_handles = [Line2D([], [], marker='o', linestyle='none', markersize=3.4,
                             markerfacecolor=_host_cmap[h], alpha=0.85,
                             markeredgecolor='white', markeredgewidth=0.4, label=h)
                      for h in _hosts]
        _leg_h = _gap_ax.legend(handles=_h_handles, loc='center left',
                                bbox_to_anchor=(_prev_r, _BAND_Y),
                                bbox_transform=_gap_ax.transAxes,
                                ncol=1, fontsize=_fs_leg,
                                title='host plant', frameon=False,
                                handletextpad=0.3, columnspacing=0.85,
                                labelspacing=0.28, borderpad=0.0,
                                borderaxespad=0.0)
        _leg_h.get_title().set_fontsize(_fs_leg)
        _leg_h.get_title().set_fontweight('bold')
        # ⚠ 同上：不要 `add_artist(_leg_h)`（会双绘）。
        try:
            fig.canvas.draw()
            if _gbb is None:                 # 无宿主表时地点图例已算过；有宿主表时补算
                _gbb = _gap_ax.get_window_extent()
            _lb2 = _leg_h.get_window_extent(fig.canvas.get_renderer())
            _prev_r = (_lb2.x1 - _gbb.x0) / max(_gbb.width, 1e-9)
        except Exception:
            _prev_r += 0.22

    # 气泡尺寸：不放色块图例（最大色块 ≈ 37 pt，在这条窄带里必然溢出）
    # → 用**一句文字**说明尺度，既不遮数据也不撑高画布。
    # ⚠ 文案要短：这条带子右侧余量有限，长句会折行（实测 "(max = 36)" 被折到下一行）。
    _vmax_txt = int(max(_all_counts or [1]))
    _note = f'area \u221d isolates; max = {_vmax_txt}'
    if _prev_r < 0.78:
        _gap_ax.text(min(_prev_r + 0.02, 0.99), _BAND_Y, _note,
                     ha='left', va='center', fontsize=_fs_leg, color='#444444',
                     transform=_gap_ax.transAxes)
    else:
        _gap_ax.text(0.0, _BAND_Y, _note, ha='left', va='center',
                     fontsize=_fs_leg, color='#444444',
                     transform=_gap_ax.transAxes)

    save_both(fig, out_path, dpi=300, bbox_inches='tight')
    plt.close(fig)
    return out_path


def _separate_bubbles(ax, fig, coords, counts, scale, *, pad_pt=1.6, iters=80):
    """把相互重叠的采样气泡在**显示坐标**里推开。

    为什么需要: 同一省内的两个采样点 (如银川 69 条 与 中宁 21 条, 相距仅约 110 km)
    在"华北—全国"尺度的地图上, 气泡圆必然相交 —— 图上看就是一团色块, 读者分不清
    那里有几个采样点。旧版直接照真实坐标画, 审稿人会读成"采样点重叠 / 坐标错误"。

    做法: 在显示坐标系里把相交的圆沿连心线等量对推, 直到两两不交 (含 pad);
    被移动的点由调用方画一条细引线连回真实坐标, 保证位置可追溯。

    返回 (adjusted: {loc: (lon, lat)}, moved: set(loc))。坐标转换用 ax.transData,
    因此**必须在最终几何 (tight_layout + canvas.draw) 之后再调用**。
    """
    import numpy as np
    locs = list(coords)
    if not locs:
        return {}, set()
    true_px = np.asarray([ax.transData.transform(coords[l]) for l in locs], float)
    rad_px = np.asarray([(max(counts[l] * scale, 12) ** 0.5) / 2.0
                         * fig.dpi / 72.0 for l in locs], float)
    P = true_px.copy()
    pad = pad_pt * fig.dpi / 72.0
    try:
        bb = ax.get_window_extent()
    except Exception:  # noqa: BLE001
        bb = None
    for _ in range(iters):
        moved = False
        for i in range(len(P)):
            for j in range(i + 1, len(P)):
                d = P[j] - P[i]
                dist = float(np.hypot(d[0], d[1]))
                need = rad_px[i] + rad_px[j] + pad
                if dist >= need:
                    continue
                if dist < 1e-6:
                    d = np.array([1.0, 0.0])
                    dist = 1.0
                u = d / dist
                push = (need - dist) / 2.0
                P[i] -= u * push
                P[j] += u * push
                moved = True
        if not moved:
            break
    if bb is not None:
        for k in range(len(P)):
            r = rad_px[k]
            P[k, 0] = min(max(P[k, 0], bb.x0 + r), bb.x1 - r)
            P[k, 1] = min(max(P[k, 1], bb.y0 + r), bb.y1 - r)
    inv = ax.transData.inverted()
    adjusted, moved_set = {}, set()
    for k, l in enumerate(locs):
        lon, lat = inv.transform(P[k])
        adjusted[l] = (float(lon), float(lat))
        # ⚠ 只有"真实点已经落到画出来的圆之外"才需要引线。位移小于半径时
        #   真实点仍在圆内，画引线只会在圆里留一道像划痕的斜线（实测如此）。
        if float(np.hypot(*(P[k] - true_px[k]))) > rad_px[k] + 1.0:
            moved_set.add(l)
    return adjusted, moved_set


def plot_sample_map(rows, location_col, out_path, log=None):
    """复刻 generate_map.R: 地图 + 采样点。

    坐标来源优先级:
      1) 输入表自带 Longitude/Latitude 列 (与上游 generate_map.R 接口一致)
      2) utils.geo_resolver 分层解析 location (取代原写死省名表)
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from utils.sci_palette import (apply_sci_style, sci_legend, COL_W_2)
    apply_sci_style(base_size=8)
    _use_cjk_font_if_available()   # 中文地名防方框 (P1-5 验证需要)

    loc_counts = Counter()
    loc_coords = {}
    unresolved = Counter()
    dropped_imprecise = Counter()

    # 地点标签统一走 resolve_label() (region 列 → geo_resolver → 'Other'),
    # 与时间面板同源 (2026-10-12 修 P1-5: 旧版地图优先 region、时间面板无视
    # region, 同一数据两面板标签/配色不一致)。
    resolver = _get_resolver()

    # ── 路径 1: 自带坐标列 (上游口径) ──
    coord_cols = load_coord_columns(rows)
    if coord_cols is not None:
        if log is not None:
            log.emit(f"[sample_plot] 使用输入表自带 Longitude/Latitude "
                     f"({len(coord_cols)}/{len(rows)} 行有效)")
        n_other = 0
        for i, (lat, lon, region) in coord_cols.items():
            name = resolve_label(rows[i], location_col, region=region,
                                 resolver=resolver)
            if name is None:
                n_other += 1
                name = 'Other'
            loc_counts[name] += 1
            loc_coords.setdefault(name, (lon, lat))
        if log is not None and n_other:
            log.emit(f"[sample_plot] 地图: {n_other} 行 region 为空且地名无法解析, "
                     f"已归入 'Other' (与时间面板同一口径)")
    else:
        # ── 路径 2: geo_resolver (上游要求预置列; 这里做增强) ──
        for r in rows:
            loc = (r.get(location_col, '') or '').strip()
            reg = _row_region(r)
            if not (loc or reg):
                continue
            name = resolve_label(r, location_col, region=reg, resolver=resolver)
            if name is None:
                # 无标签: 只到国家级 (如 'China') 与真·未识别分开记数
                # (国家级无法定位到采样点, 记录在日志里, 不静默丢数据)
                if loc and not _is_specific_label(_canonical_label(loc)):
                    dropped_imprecise[loc] += 1
                else:
                    unresolved[loc or reg] += 1
                continue
            coord = None
            if resolver is not None and loc:
                coord, _src = resolver.resolve(loc)
            if coord is None:
                unresolved[loc or reg] += 1
                continue
            lat, lon = coord
            loc_counts[name] += 1
            loc_coords.setdefault(name, (lon, lat))

        if log is not None and dropped_imprecise:
            log.emit(f"[sample_plot] 地图: 剔除 {sum(dropped_imprecise.values())} 条"
                     f"仅到国家级的记录（无法定位到采样点）: "
                     + ", ".join(f"{k}({v})" for k, v in dropped_imprecise.items()))
        if log is not None and unresolved:
            n_unres = sum(unresolved.values())
            frac = n_unres / max(len(rows), 1)
            msg = (f"[sample_plot] 地图坐标解析: {n_unres}/{len(rows)} "
                   f"({frac*100:.1f}%) 未识别, 已从地图剔除")
            if frac >= UNRESOLVED_WARN_FRACTION:
                log.warning(msg + " —— 占比过高, 地图可能严重缺采样点。"
                                  "未识别清单: " + ", ".join(
                                      f"{k}({v})" for k, v in unresolved.most_common(10)))
            else:
                log.emit(msg)

    fig = plt.figure(figsize=(COL_W_2, COL_W_2 * 1.02))

    # 有真实坐标时按数据范围自适应; 否则退回全球视图 (上游 bbox = 全球)
    if loc_coords:
        lons = [c[0] for c in loc_coords.values()]
        lats = [c[1] for c in loc_coords.values()]
        pad_x = max(5.0, (max(lons) - min(lons)) * 0.15)
        pad_y = max(5.0, (max(lats) - min(lats)) * 0.15)
        extent = [min(lons) - pad_x, max(lons) + pad_x,
                  min(lats) - pad_y, max(lats) + pad_y]
    else:
        extent = [-180, 180, -90, 90]

    try:
        import cartopy.crs as ccrs
        import cartopy.feature as cfeature
        ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
        ax.add_feature(cfeature.COASTLINE, lw=0.5)
        ax.add_feature(cfeature.BORDERS, lw=0.4)
        ax.set_extent(extent, crs=ccrs.PlateCarree())
        ax.gridlines(draw_labels=True, lw=0.3, alpha=0.5)
        transform = ccrs.PlateCarree()
    except Exception:
        ax = fig.add_subplot(1, 1, 1)
        ax.set_xlim(extent[0], extent[1])
        ax.set_ylim(extent[2], extent[3])
        transform = None

    # 2026-10-12（审稿人视角）：颜色改走 loc_color_map()，与时间面板 B 一一对应。
    cmap = loc_color_map(loc_counts.keys())
    # 气泡面积 ∝ 分离株数 (s = count × BUBBLE_S)。
    # 2026-10-12 二次下调 20 → 8 → 5：vmax=69 时直径 2·√(69×5/π) ≈ 21 pt ≈ 0.29 in，
    # 占 4.58 in 图宽的 6.3%（20 时为 12.7%，8 时仍与中宁的圆相交）。缩小后
    # 银川(69)/中宁(21) 两圆不再相切，同时仍是清晰的视觉层级。
    BUBBLE_S = 5                       # pt² / isolate
    ax.set_title('Sampling locations', fontsize=9, fontweight='bold')

    # ── 定位插图（合规 + 审稿人可读性）─────────────────────────────
    # 主图只覆盖研究区（西北—华北），必须让读者看到研究区在全国的位置；
    # 插图把国土全境（含中国台湾与南海诸岛所在海域）一并框入，避免出现
    # "只画一半"的残缺地图。仅画海岸线（物理地理），不画国界。
    if transform is not None:
        try:
            import cartopy.crs as ccrs
            import matplotlib.patches as mpatches
            axins = ax.inset_axes([0.015, 0.015, 0.30, 0.27],
                                  projection=ccrs.PlateCarree())
            axins.set_extent([73, 136, 2, 55], crs=ccrs.PlateCarree())
            axins.add_feature(cfeature.COASTLINE, lw=0.25)
            axins.add_patch(mpatches.Rectangle(
                (extent[0], extent[2]), extent[1] - extent[0],
                extent[3] - extent[2], fill=False, edgecolor='#D55E00',
                linewidth=0.8, transform=ccrs.PlateCarree(), zorder=6))
            axins.set_xticks([])
            axins.set_yticks([])
            axins.set_facecolor('white')
            for sp in axins.spines.values():
                sp.set_edgecolor('#888888')
                sp.set_linewidth(0.5)
            axins.set_title('Study region', fontsize=5.2, pad=1.0)
        except Exception:
            pass

    if loc_counts:
        # 图例放到坐标区外下方，横向排布 —— 避免压住采样点（旧版压在地图左下角）
        # 2026-10-12: 图例标记必须**统一小尺寸**。直接让 legend 复用 scatter 句柄会
        # 继承 `s=max(count*20,12)` —— 最大地点的圆大到把相邻标签全遮住（实测已发生）。
        from matplotlib.lines import Line2D
        # 2026-10-12（审稿人视角）：图例色块必须与地图气泡**同透明度**。
        # 气泡用 alpha=0.9 叠在白底上（实测 #D55E00 → RGB(217,110,25)），
        # 图例若用全不透明色，同一地点在图例与地图上会呈现两种色调 —— 审稿人
        # 会当成配色不一致。这里显式给 alpha，让"单层气泡"的颜色与图例逐像素相同。
        handles = [Line2D([], [], marker='o', linestyle='none', markersize=4.5,
                          markerfacecolor=cmap[loc], alpha=0.9,
                          markeredgecolor='white', markeredgewidth=0.5, label=loc)
                   for loc in loc_counts]
        sci_legend(ax, handles=handles, loc='upper center',
                   bbox_to_anchor=(0.5, -0.06),
                   ncol=min(4, max(1, len(loc_counts))), fontsize=6.5,
                   title=None, frameon=True,
                   handletextpad=0.4, columnspacing=1.0)
    # ⚠ 先定版、再放尺寸图例：tight_layout 会改变坐标区几何，若在图例放好之后
    #   才 tight_layout，避让计算就基于旧几何、可能失效。
    fig.tight_layout(pad=0.4)

    # ── 采样气泡：先分离重叠、再落点 ─────────────────────────────
    # 必须在 tight_layout + canvas.draw 之后做：分离算法用的是 ax.transData
    # 显示坐标，几何未定版时算出的推挤量是错的。未被移动的点 adj 等于真实坐标，
    # 行为与旧版完全一致；只有真正相交的点才会被推开并补一条引线。
    adj, moved = {}, set()
    if loc_counts and transform is not None:
        try:
            fig.canvas.draw()
            adj, moved = _separate_bubbles(ax, fig, loc_coords, loc_counts,
                                           BUBBLE_S)
        except Exception:      # noqa: BLE001 — 分离失败即退回真实坐标
            adj, moved = {}, set()
    for loc, count in loc_counts.items():
        lon, lat = adj.get(loc, loc_coords[loc])
        kw = dict(s=max(count * BUBBLE_S, 12), color=cmap[loc], label=loc,
                  edgecolor='white', linewidth=0.6, alpha=0.9, zorder=3)
        if transform is not None:
            kw['transform'] = transform
        ax.scatter(lon, lat, **kw)
        if loc in moved:
            tlon, tlat = loc_coords[loc]
            lk = dict(color='#AAAAAA', lw=0.5, zorder=2)
            if transform is not None:
                lk['transform'] = transform
            ax.plot([tlon, lon], [tlat, lat], **lk)

    if loc_counts:
        # 气泡大小图例（审稿人必问：圆的大小代表多少条序列）
        # 自动避让：尺寸图例的最大色块是按真实比例画的（vmax=69 → 直径≈37 pt），
        # 写死角落会压住采样点。这里按当前数据逐个候选位置试算，并把定位插图
        # 也作为避让框（插图在左下角，恰好是候选之一）。
        from matplotlib.transforms import Bbox
        try:
            (_ax0, _ay0), (_ax1, _ay1) = ax.transAxes.transform(
                [[0.015, 0.015], [0.315, 0.285]])
            avoid = [Bbox([[_ax0, _ay0], [_ax1, _ay1]])]
        except Exception:
            avoid = []
        _place_size_legend(ax, list(loc_counts.values()), BUBBLE_S, fig,
                           [('upper left', (0.012, 0.99)),
                            ('lower right', (0.99, 0.012)),
                            ('upper right', (0.99, 0.99))],
                           title='isolates', fontsize=5.6, extra_boxes=avoid)
    save_both(fig, out_path, dpi=300, bbox_inches='tight')
    plt.close(fig)
    return out_path


def main():
    ap = argparse.ArgumentParser(description='VirPhyKit SamplePlot 复刻')
    ap.add_argument('--metadata', required=True, help='元数据 CSV (含日期+地点)')
    ap.add_argument('--date-col', default='date')
    ap.add_argument('--location-col', default='location')
    ap.add_argument('--host-col', default='host')
    ap.add_argument('--host-metadata', default=None,
                    help='第二张元数据表 (含 name + 宿主列)。给了就在时间分布图的'
                         '气泡里按宿主切饼 (2026-10-13 用户要求「物种也加上」)。')
    ap.add_argument('--outdir', default='sample_plot_out')
    args = ap.parse_args()

    class _L:
        def emit(self, m): print(m)
        def warning(self, m): print(m, file=sys.stderr)
    log = _L()

    # 自动检测分隔符 (tab 或逗号)
    with open(args.metadata, encoding='utf-8-sig') as f:
        sample = f.read(2048)
        delim = '\t' if sample.count('\t') > sample.count(',') else ','
    with open(args.metadata, encoding='utf-8-sig') as f:
        rows = list(csv.DictReader(f, delimiter=delim))
    print(f'样本数: {len(rows)} (分隔符: {"tab" if delim == chr(9) else "逗号"})')

    # 宿主表按 name 并入 (列名冲突时以主表为准, 宿主列另起名)
    host_col = None
    if args.host_metadata and os.path.exists(args.host_metadata):
        with open(args.host_metadata, encoding='utf-8-sig') as f:
            hs = f.read(2048)
            hd = '\t' if hs.count('\t') > hs.count(',') else ','
        with open(args.host_metadata, encoding='utf-8-sig') as f:
            hrows = {r.get('name', '').strip(): r for r in csv.DictReader(f, delimiter=hd)}
        _n = 0
        for r in rows:
            h = hrows.get((r.get('name') or '').strip())
            if h:
                v = (h.get(args.host_col) or h.get('host') or '').strip()
                if v:
                    r['__host__'] = v
                    _n += 1
        host_col = '__host__' if _n else None
        print(f'宿主表并入: {_n}/{len(rows)} 条命中 ({args.host_metadata})')
    else:
        print('宿主表未提供 → 时间分布图不画宿主饼')

    os.makedirs(args.outdir, exist_ok=True)

    # 1. 时间分布图
    t_path = os.path.join(args.outdir, 'sample_temporal.pdf')
    plot_sample_temporal(rows, args.date_col, args.location_col, t_path,
                         log=log, host_col=host_col)
    print(f'采样时间分布图: {t_path}')

    # 2. 地理地图
    m_path = os.path.join(args.outdir, 'sample_map.pdf')
    plot_sample_map(rows, args.location_col, m_path, log=log)
    print(f'采样地理地图: {m_path}')

    # 落盘 geo 缓存 (在线解析结果供下次秒查)
    try:
        from utils.geo_resolver import get_resolver
        get_resolver().flush()
    except Exception:
        pass

    return 0


if __name__ == '__main__':
    sys.exit(main())
