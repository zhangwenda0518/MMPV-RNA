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


def build_temporal_data(rows, date_col, location_col, log=None):
    """从 metadata 生成 Year/Total/地点 宽表。

    地点归并优先用 geo_resolver (分层解析), 解析不到才落 'Other'。
    返回 (years, locs, year_loc, unresolved_places)
    """
    resolver = _get_resolver()
    year_loc = defaultdict(Counter)
    unresolved = Counter()          # 未解析成功的原始地名 → 出现次数

    for r in rows:
        year = extract_year(r.get(date_col, ''))
        loc = (r.get(location_col, '') or '').strip()
        if not (year and loc):
            continue
        label = None
        if resolver is not None:
            coord, src = resolver.resolve(loc)
            if coord is not None:
                label = _canonical_label(loc)
        if label is None:
            unresolved[loc] += 1
            label = 'Other'
        year_loc[year][label] += 1

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
    return years, locs, year_loc, unresolved


def plot_sample_temporal(rows, date_col, location_col, out_path, log=None):
    """复刻 generate_plot.R: 时间线折线 + 地点分布点线 (两面板)"""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    years, locs, year_loc, _ = build_temporal_data(
        rows, date_col, location_col, log=log)
    totals = {y: sum(year_loc[y].values()) for y in years}

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 6), sharex=True,
                                   gridspec_kw={'height_ratios': [3, 2]})
    # 空/无有效记录守卫 (历史坑: location 列全空或列名不符时 min/max 直接 ValueError)
    if not years or not totals:
        print("  [!] 无有效 年份×地点 记录, 跳过采样分布图")
        plt.close(fig)
        return
    # 上: 总样本时间线 (复刻 p1)
    ax1.plot(years, [totals[y] for y in years], color='red', lw=1.5)
    ax1.set_ylabel('Isolates', fontsize=12)
    ax1.set_xlim(min(years), max(years))
    ax1.set_ylim(0, max(totals.values()) + 2)
    ax1.grid(False)
    ax1.tick_params(axis='x', labelrotation=25)

    # 下: 各地点分布点线 (复刻 p2)
    colors = plt.cm.Set1(range(max(len(locs), 1)))
    for i, loc in enumerate(locs):
        xs = [y for y in years if year_loc[y].get(loc, 0) > 0]
        ys = [year_loc[y].get(loc, 0) for y in xs]
        ax2.scatter(xs, [i] * len(xs), s=[v * 30 for v in ys],
                    color=colors[i % len(colors)], label=loc, alpha=0.7)
        ax2.plot(xs, [i] * len(xs), color=colors[i % len(colors)], lw=0.5)
    ax2.set_yticks(range(len(locs)))
    ax2.set_yticklabels(locs, fontsize=8)
    ax2.set_xlabel('Year', fontsize=12)
    ax2.tick_params(axis='x', labelrotation=25)

    fig.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches='tight')
    plt.close(fig)
    return out_path


def plot_sample_map(rows, location_col, out_path, log=None):
    """复刻 generate_map.R: 地图 + 采样点。

    坐标来源优先级:
      1) 输入表自带 Longitude/Latitude 列 (与上游 generate_map.R 接口一致)
      2) utils.geo_resolver 分层解析 location (取代原写死省名表)
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    loc_counts = Counter()
    loc_coords = {}
    unresolved = Counter()

    # ── 路径 1: 自带坐标列 (上游口径) ──
    coord_cols = load_coord_columns(rows)
    if coord_cols is not None:
        if log is not None:
            log.emit(f"[sample_plot] 使用输入表自带 Longitude/Latitude "
                     f"({len(coord_cols)}/{len(rows)} 行有效)")
        for i, (lat, lon, region) in coord_cols.items():
            name = region or _canonical_label(rows[i].get(location_col, '')) \
                   or f"row{i}"
            loc_counts[name] += 1
            loc_coords.setdefault(name, (lon, lat))
    else:
        # ── 路径 2: geo_resolver (上游要求预置列; 这里做增强) ──
        resolver = _get_resolver()
        for r in rows:
            loc = (r.get(location_col, '') or '').strip()
            if not loc:
                continue
            name = _canonical_label(loc)
            if resolver is None:
                unresolved[loc] += 1
                continue
            coord, src = resolver.resolve(loc)
            if coord is None:
                unresolved[loc] += 1
                continue
            lat, lon = coord
            loc_counts[name] += 1
            loc_coords.setdefault(name, (lon, lat))

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

    fig = plt.figure(figsize=(10, 8))

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

    colors = plt.cm.Set1(range(max(len(loc_counts), 1)))
    for i, (loc, count) in enumerate(loc_counts.items()):
        lon, lat = loc_coords[loc]
        kw = dict(s=count * 40, color=colors[i % len(colors)], label=loc,
                  edgecolor='white', linewidth=0.8)
        if transform is not None:
            kw['transform'] = transform
        ax.scatter(lon, lat, **kw)
    ax.set_title('Sampling Locations', fontsize=14)
    if loc_counts:
        ax.legend(loc='lower left', fontsize=9, frameon=False)
    fig.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches='tight')
    plt.close(fig)
    return out_path


def main():
    ap = argparse.ArgumentParser(description='VirPhyKit SamplePlot 复刻')
    ap.add_argument('--metadata', required=True, help='元数据 CSV (含日期+地点)')
    ap.add_argument('--date-col', default='date')
    ap.add_argument('--location-col', default='location')
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

    os.makedirs(args.outdir, exist_ok=True)

    # 1. 时间分布图
    t_path = os.path.join(args.outdir, 'sample_temporal.pdf')
    plot_sample_temporal(rows, args.date_col, args.location_col, t_path, log=log)
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
