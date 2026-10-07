#!/usr/bin/env python3
"""coord_lookup.py — 地理名称 → 坐标 的**唯一来源** (2026-10-03, S3-⑤)

背景
----
同一批地名在 virome_phylo_pipeline 里存在四/五份坐标表, 且互不一致
(审计实测: 'Ningxia' 五表五样; 'Qinghai' 一表给省质心 (36.28, 95.24)、
另一表给省会西宁 (36.62, 101.78) —— 相差约 586 km)。Mantel / SpreaD3 /
metadata clean 各读各的表 → 同一份数据换个模块就得到不同的地理距离。

本模块把全部表合并为**一份权威查找**:

  基底 (权威):  utils/geo_resolver.PROVINCE_LOOKUP / CITY_LOOKUP
  合并 (补缺):  virphy_bridge.COORD_LOOKUP       (含裸省份名 + 复合地名)
               spread3_viz._DEFAULT_LAT_LON      (含 Europe/Africa/Other 等占位)
               metadata_governance._CITY_COORDS  (30 个中国城市手校表)
               transmission_map                  (无自有表, 本就走 geo_resolver)

冲突消解规则 (S3-⑤ 约定):
  · 省份键 → geo_resolver **省质心**;  城市键 → 城市表 (CITY_LOOKUP)
  · 同一键多表不一致 → 取 geo_resolver 值, 差异记入差异清单 (diff_lines)
  · 消费方旧表与合并表取值不同的键 → 消费方打印 WARN (键/旧值/新值/来源)

⚠ 不改 submission 侧 unified_metadata.py 的 30 城表 (跨模块保留, 不在本管线的
  坐标源收敛范围内)。

API
---
  lookup(name) → (lat, lon, source, granularity) | None
      归一化匹配 (大小写/空白/逗号/别名); 静态合并表未命中时回落到
      geo_resolver 的离线分层解析 (段拆分/模糊/缓存; 在线仅显式开启时)。
  diff_lines(label) → [str, ...]
      label ∈ {'virphy_bridge','spread3_viz','metadata_governance'};
      返回该消费方旧表与合并表取值不同的差异行, 供消费方打印 WARN。
  get_resolver() / flush()   透传 geo_resolver (缓存落盘由消费方 flush)。
"""
from __future__ import annotations

import os
import re
from typing import Dict, List, Optional, Tuple

__all__ = ["lookup", "diff_lines", "flush", "get_resolver", "MERGED", "DIFF_REPORT"]


# ═══════════════════════════════════════════════════════════════════
# geo_resolver 加载 (包内导入 → 顶层导入 → 按文件路径裸加载, 三级兜底)
# ═══════════════════════════════════════════════════════════════════
def _load_geo_resolver():
    """加载 geo_resolver 模块。

    virphy_bridge 有 spec_from_file_location 裸加载路径 (无包上下文), 本模块
    也必须能在同样条件下工作 → 不依赖 utils 包已在 sys.path。
    """
    try:
        from . import geo_resolver as _g          # 包内正常导入
        return _g
    except ImportError:
        pass
    try:
        import geo_resolver as _g                 # utils/ 已在 sys.path
        return _g
    except ImportError:
        pass
    try:
        from utils import geo_resolver as _g      # 仓库根在 sys.path
        return _g
    except ImportError:
        pass
    import importlib.util
    _path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'geo_resolver.py')
    try:
        _spec = importlib.util.spec_from_file_location('_coord_lookup_geo_resolver', _path)
        _mod = importlib.util.module_from_spec(_spec)
        _spec.loader.exec_module(_mod)
        return _mod
    except Exception as _e:  # noqa: BLE001 — 兜底: 不阻断消费方, 退化为旧表-only
        print(f"[coord_lookup] WARN: geo_resolver 加载失败 ({_e}); "
              f"权威基底为空表, 仅旧表条目可用")
        return None


_GEO = _load_geo_resolver()
if _GEO is None:                     # 极端部署异常: 空基底, 不崩
    class _GeoStub:
        PROVINCE_LOOKUP = {}
        CITY_LOOKUP = {}
    _GEO = _GeoStub()

_PROVINCE_BASE: Dict[str, Tuple[float, float]] = dict(getattr(_GEO, 'PROVINCE_LOOKUP', {}))
_CITY_BASE: Dict[str, Tuple[float, float]] = dict(getattr(_GEO, 'CITY_LOOKUP', {}))


# ── 归一化 / 占位符判据 (与 geo_resolver 保持一致; 加载失败时用本地等价实现) ──
def _norm_fallback(s: str) -> str:
    s = str(s or '').strip().lower()
    s = re.sub(r'[\s,]+', ' ', s)
    s = re.sub(r'\b(province|sheng|region|autonomous)\b', '', s).strip()
    s = s.replace('province', '').strip()
    return s


def _placeholder_fallback(key: str) -> bool:
    _PH = {'', 'unknown', 'na', 'n/a', 'not applicable', 'none', 'null',
           'missing', '其他', '未知'}
    k = str(key or '').strip().lower()
    if not k or k in _PH:
        return True
    core = k.strip(':').strip()
    if not core or core in _PH:
        return True
    return bool(re.fullmatch(r"[x\s.]+", core) or re.fullmatch(r"[<\-\s]+", k))


_norm = getattr(_GEO, '_norm', None) or _norm_fallback
_is_placeholder_key = getattr(_GEO, '_is_placeholder_key', None) or _placeholder_fallback


# ═══════════════════════════════════════════════════════════════════
# 旧消费方表 (合并进唯一源; 也是 diff_lines 的对比基准)
# ═══════════════════════════════════════════════════════════════════

# ── virphy_bridge.COORD_LOOKUP (2026-09-02 Mantel 用, 现已移除本地字面量) ──
_LEGACY_VIRPHY: Dict[str, Tuple[float, float]] = {
    'China, Ningxia, Yinchuan': (38.47, 106.27),
    'China, Ningxia, Zhongning': (37.49, 105.68),
    'China, Neimenggu': (40.82, 111.75),
    'China, Beijing, Beijing': (39.90, 116.40),
    'China, Beijing': (39.90, 116.40),
    'China, Beijing, Unknown': (39.90, 116.40),
    'China, Guangdong, Guangzhou': (23.13, 113.26),
    'China, Qinghai, Xining': (36.62, 101.78),
    'China, Gansu, Wuwei': (37.93, 102.64),
    'China, Liaoning, Shenyang': (41.80, 123.43),
    'China, ningxia, Yinchuan': (38.47, 106.27),
    # 裸省份名 (PSTVd metadata 格式)
    'Ningxia': (38.47, 106.27),
    'Gansu': (36.06, 103.83),
    'Qinghai': (36.62, 101.78),
    'Neimenggu': (40.82, 111.75),
    'Inner Mongolia': (40.82, 111.75),
    'Beijing': (39.90, 116.40),
    'Guangdong': (23.13, 113.26),
    'Jiangsu': (32.06, 118.78),
    'Shandong': (36.67, 117.02),
    'Hebei': (38.04, 114.51),
    'Xinjiang': (43.79, 87.63),
    'Yunnan': (25.04, 102.71),
    'Sichuan': (30.66, 104.07),
    'Hunan': (28.23, 112.94),
    'Hubei': (30.59, 114.31),
    'Henan': (34.75, 113.63),
    'Shaanxi': (34.27, 108.95),
    'Shanxi': (37.87, 112.55),
    'Liaoning': (41.80, 123.43),
    'Fujian': (26.08, 119.30),
    'Zhejiang': (30.27, 120.16),
    'Anhui': (31.86, 117.28),
    'Heilongjiang': (45.80, 126.53),
    'Jilin': (43.88, 125.32),
    'Guangxi': (22.82, 108.32),
    'Guizhou': (26.65, 106.63),
    'Hainan': (20.04, 110.32),
    'Chongqing': (29.56, 106.55),
    'Tianjin': (39.08, 117.20),
    'Shanghai': (31.23, 121.47),
    'Tibet': (29.65, 91.13),
}

# ── spread3_viz._DEFAULT_LAT_LON (SpreaD3 JSON 内置表) ──
_LEGACY_SPREAD3: Dict[str, Tuple[float, float]] = {
    "Beijing": (39.9, 116.4),
    "Ningxia": (37.3, 106.0),
    "Neimenggu": (44.0, 113.0),
    "China": (35.0, 105.0),
    "Japan": (36.0, 138.0),
    "Korea": (37.0, 127.5),
    "India": (20.0, 78.0),
    "Europe": (50.0, 10.0),
    "USA": (40.0, -100.0),
    "South America": (-15.0, -60.0),
    "Africa": (0.0, 25.0),
    "Other": (30.0, 80.0),
}

# ── metadata_governance._CITY_COORDS (30 个中国城市手校表; NCBI 串格式已转 float) ──
_LEGACY_METADATA_GOV: Dict[str, Tuple[float, float]] = {
    'yinchuan': (38.47, 106.27),
    'ningxia': (37.48, 105.68),
    'zhongning': (37.48, 105.68),
    'beijing': (39.90, 116.40),
    'shanghai': (31.23, 121.47),
    'guangzhou': (23.13, 113.26),
    'wuhan': (30.59, 114.31),
    'nanjing': (32.06, 118.79),
    'hangzhou': (30.27, 120.15),
    'chengdu': (30.57, 104.07),
    'xian': (34.26, 108.94),
    'kunming': (25.04, 102.68),
    'harbin': (45.80, 126.53),
    'zhengzhou': (34.75, 113.62),
    'jinan': (36.65, 116.98),
    'taiyuan': (37.87, 112.55),
    'changsha': (28.23, 112.94),
    'fuzhou': (26.07, 119.30),
    'guiyang': (26.65, 106.63),
    'lanzhou': (36.06, 103.79),
    'xining': (36.62, 101.77),
    'urumqi': (43.79, 87.58),
    'lhasa': (29.65, 91.10),
    'shenyang': (41.80, 123.43),
    'dalian': (38.91, 121.61),
    'qingdao': (36.07, 120.38),
    'suzhou': (31.30, 120.62),
    'shenzhen': (22.54, 114.06),
    'chongqing': (29.56, 106.55),
    'tianjin': (39.13, 117.18),
}

_LEGACY_TABLES = {
    'virphy_bridge': _LEGACY_VIRPHY,
    'spread3_viz': _LEGACY_SPREAD3,
    'metadata_governance': _LEGACY_METADATA_GOV,
}

# ── geo_resolver PROVINCE_LOOKUP 里的"国家"块 (其余为省级行政区; 港澳台按
#    中国一级行政区计, granularity='province') ──
_COUNTRY_KEYS = frozenset({
    'china', '中国', 'japan', 'korea', 'india', 'vietnam', 'thailand', 'usa',
    'united states', 'brazil', 'australia', 'russia', 'germany', 'france',
    'uk', 'united kingdom', 'unitedkingdom', 'newzealand', 'new zealand',
    'peru', 'singapore', 'denmark', 'malaysia', 'netherlands', 'nicaragua',
    'canada', 'mexico', 'iran', 'guam', 'philippines', 'indonesia', 'italy',
    'spain', 'sweden', 'norway', 'finland', 'poland', 'belgium', 'switzerland',
    'austria', 'portugal', 'ireland', 'greece', 'turkey', 'egypt',
    'south africa', 'argentina', 'chile', 'colombia', 'cambodia', 'laos',
    'myanmar', 'bangladesh', 'pakistan', 'sri lanka', 'nepal', 'mongolia',
    'kazakhstan', 'saudi arabia', 'israel', 'kenya', 'qatar', 'serbia',
    'croatia', 'czech', 'hungary', 'romania', 'bulgaria', 'kuwait', 'oman',
    'jordan', 'lebanon', 'morocco', 'tunisia', 'nigeria', 'ethiopia', 'ghana',
    'tanzania', 'madagascar', 'uruguay', 'paraguay', 'bolivia', 'ecuador',
    'venezuela', 'cuba', 'jamaica', 'panama', 'costa rica', 'guatemala',
    'honduras', 'brunei', 'fiji', 'papua new guinea',
})
_CONTINENT_KEYS = frozenset({'europe', 'south america', 'africa', 'asia',
                             'north america', 'oceania'})

_PROVINCE_BASE_NORM = {_norm(k) for k in _PROVINCE_BASE}
_CITY_BASE_NORM = {_norm(k) for k in _CITY_BASE}


def _granularity_for_province(nk: str) -> str:
    if nk in _CITY_BASE_NORM:
        return 'city'          # 基表同键既是省又是市时按市 (更具体)
    return 'province' if nk not in _COUNTRY_KEYS else 'country'


def _granularity_for_legacy(nk: str, label: str) -> str:
    """旧表键的粒度推断: 命中省/市基表 → 对应粒度; 复合地名看最后一段。"""
    if label == 'metadata_governance':
        # 该表按定义全部是城市 (30 个中国城市手校表); 部分城市 (wuhan 等)
        # 不在 geo_resolver.CITY_LOOKUP, 不能因此降级成 'other'。
        return 'city'
    if nk in _CITY_BASE_NORM:
        return 'city'
    if nk in _PROVINCE_BASE_NORM:
        return _granularity_for_province(nk)
    if nk in _CONTINENT_KEYS:
        return 'continent'
    if nk == 'other':
        return 'other'
    last = nk.split(' ')[-1]
    if last in _CITY_BASE_NORM:
        return 'city'
    if last in _PROVINCE_BASE_NORM:
        return _granularity_for_province(last)
    return 'other'


# ═══════════════════════════════════════════════════════════════════
# 合并表构建 (geo_resolver 权威; 旧表补缺; 冲突记录)
# ═══════════════════════════════════════════════════════════════════
def _build_merged():
    merged: Dict[str, Tuple[float, float, str, str]] = {}
    conflicts: List[str] = []

    # 1) geo_resolver 基底: 省/市表 (权威)
    for k, v in _PROVINCE_BASE.items():
        nk = _norm(k)
        merged[nk] = (float(v[0]), float(v[1]), 'geo_resolver:province',
                      _granularity_for_province(nk))
    for k, v in _CITY_BASE.items():
        nk = _norm(k)
        # 城市表晚于省表 → 同键时城市表覆盖 (更具体; 如 'beijing city' 不与 'beijing' 冲突)
        merged[nk] = (float(v[0]), float(v[1]), 'geo_resolver:city', 'city')

    # 2) 旧表: 补缺 + 冲突记录 (geo_resolver 值胜出)
    for label, table in _LEGACY_TABLES.items():
        for k, v in table.items():
            nk = _norm(k)
            if not nk:
                continue
            v = (float(v[0]), float(v[1]))
            if nk in merged:
                cur = merged[nk]
                if (round(cur[0], 4), round(cur[1], 4)) != (round(v[0], 4), round(v[1], 4)):
                    conflicts.append(
                        f"{label}: {k!r} 旧值 ({v[0]}, {v[1]}) → 新值 "
                        f"({cur[0]}, {cur[1]}) ← {cur[2]}")
            else:
                merged[nk] = (v[0], v[1], f'legacy:{label}',
                              _granularity_for_legacy(nk, label))
    return merged, conflicts


MERGED, DIFF_REPORT = _build_merged()


# ═══════════════════════════════════════════════════════════════════
# 查询
# ═══════════════════════════════════════════════════════════════════
def lookup(name) -> Optional[Tuple[float, float, str, str]]:
    """地名 → (lat, lon, source, granularity); 解析失败返回 None。

    source: 'geo_resolver:province' / 'geo_resolver:city' /
            'legacy:<消费方>' / 'resolver:<exact|city|province|fuzzy|cache|online>'
    granularity: province / city / country / continent / other / resolved
    """
    s = str(name or '')
    nk = _norm(s)
    if not nk or _is_placeholder_key(nk):
        return None
    hit = MERGED.get(nk)
    if hit:
        return (hit[0], hit[1], hit[2], hit[3])
    # 静态表未命中 → geo_resolver 离线分层解析 (段拆分/模糊/缓存; 在线需显式开启)
    try:
        coord, src = get_resolver().resolve(s)
    except Exception:  # noqa: BLE001 — 解析器异常不致命, 视为未命中
        coord, src = None, 'miss'
    if coord:
        return (float(coord[0]), float(coord[1]), f'resolver:{src}', 'resolved')
    return None


def diff_lines(label: str) -> List[str]:
    """消费方旧表 vs 合并表的差异行 (键/旧值/新值/来源), 供打印 WARN。"""
    table = _LEGACY_TABLES.get(label)
    if not table:
        return []
    lines = []
    for k, v in table.items():
        nk = _norm(k)
        cur = MERGED.get(nk)
        if cur is None:
            lines.append(f"  坐标差异({label}): 键 {k!r} 不在合并表 (内部错误)")
            continue
        if (round(cur[0], 4), round(cur[1], 4)) != (round(float(v[0]), 4), round(float(v[1]), 4)):
            lines.append(f"  坐标差异({label}): {k!r} 旧值 ({v[0]}, {v[1]}) → "
                         f"新值 ({cur[0]}, {cur[1]}) ← {cur[2]}")
    return lines


_RESOLVER = None


def get_resolver():
    """geo_resolver 单例 (惰性; 与本模块同一模块实例, flush 共享同一缓存)。

    geo_resolver 不可用时返回 None (静态合并表仍可用; 不崩消费方)。
    """
    global _RESOLVER
    if _RESOLVER is None:
        try:
            _RESOLVER = _GEO.get_resolver()
        except AttributeError:
            return None
    return _RESOLVER


def flush() -> None:
    """把 geo_resolver 缓存落盘 (消费方循环结束后调用)。"""
    _r = get_resolver()
    if _r is not None:
        _r.flush()
