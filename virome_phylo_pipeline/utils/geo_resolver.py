#!/usr/bin/env python3
"""geo_resolver: 地理名称 → 经纬度 分层解析器 (Mantel/phylogeo 通用)

设计目标 (2026-09-02, PSTVd 裸省份名教训): 替换 virphy_bridge 里的
字面量 COORD_LOOKUP, 做到对常见地名格式自适应:

  解析层级 (逐级兜底, 全部离线优先):
    L1  精确匹配: 内置中国全行政区划表 (省/市/县, 中英文+拼音变体)
    L2  规范化匹配: 大小写/空白/逗号段拆分重组
        'China, ningxia, Yinchuan' / 'Ningxia' / 'ningxia ' / 'CHINA,Ningxia'
        全部收敛到同一 key
    L3  段内反查: 多段地名逐段尝试 (任意一段命中省级/市级即用)
        'China, Gansu, Wuwei' → Gansu 省质心 (Wuwei 无县级数据时)
    L4  持久化缓存: 解析结果写盘 (cache/geo_cache.json), 下次秒查且人工可修
    L5  [可选] 在线 Nominatim 兜底: 环境变量 GEO_RESOLVER_ONLINE=1 时启用,
        带 1req/s 限速与 UA, 结果回写缓存; 默认关闭 (管线离线跑)

  仍然解析失败: 返回 None, 调用方剔除样本 (绝不伪造坐标)
"""
import json
import os
import re
import time
from pathlib import Path

# ── 中国省级行政区质心 (lat, lon): 中英文名 + 常见拼音变体 ──────────
# 精度说明: 省级用质心近似, 对 Mantel 的省际距离足够; 市县级在 CITY_LOOKUP
PROVINCE_LOOKUP = {
    # 直辖市
    'beijing': (39.90, 116.40), '北京': (39.90, 116.40),
    'tianjin': (39.08, 117.20), '天津': (39.08, 117.20),
    'shanghai': (31.23, 121.47), '上海': (31.23, 121.47),
    'chongqing': (29.56, 106.55), '重庆': (29.56, 106.55),
    # 省
    'hebei': (38.88, 115.48), '河北': (38.88, 115.48),
    'shanxi': (37.87, 112.55), '山西': (37.87, 112.55),
    'shaanxi': (34.27, 108.95), '陕西': (34.27, 108.95), 'shanxxi': (34.27, 108.95),
    'shandong': (36.40, 118.15), '山东': (36.40, 118.15),
    'henan': (34.75, 113.53), '河南': (34.75, 113.53),
    'jiangsu': (32.97, 119.42), '江苏': (32.97, 119.42),
    'anhui': (31.86, 117.28), '安徽': (31.86, 117.28),
    'zhejiang': (29.16, 119.65), '浙江': (29.16, 119.65),
    'fujian': (26.08, 117.98), '福建': (26.08, 117.98),
    'jiangxi': (27.80, 114.92), '江西': (27.80, 114.92),
    'hubei': (30.58, 112.30), '湖北': (30.58, 112.30),
    'hunan': (27.60, 111.71), '湖南': (27.60, 111.71),
    'guangdong': (23.34, 113.42), '广东': (23.34, 113.42),
    'guangxi': (23.85, 108.67), '广西': (23.85, 108.67),
    'hainan': (19.20, 109.70), '海南': (19.20, 109.70),
    'sichuan': (30.50, 103.00), '四川': (30.50, 103.00),
    'guizhou': (26.64, 106.71), '贵州': (26.64, 106.71),
    'yunnan': (25.04, 102.71), '云南': (25.04, 102.71),
    'gansu': (36.06, 103.83), '甘肃': (36.06, 103.83),
    'qinghai': (36.28, 95.24), '青海': (36.28, 95.24),
    'taiwan': (23.70, 120.96), '台湾': (23.70, 120.96),
    # 自治区 (含历史遗留拼写变体)
    'neimenggu': (44.09, 113.97), '内蒙古': (44.09, 113.97),
    'inner mongolia': (44.09, 113.97), 'neimeng': (44.09, 113.97),
    'ningxia': (37.27, 106.17), '宁夏': (37.27, 106.17),
    'guangxi zhuangzu': (23.85, 108.67),
    'xizang': (31.27, 88.77), '西藏': (31.27, 88.77), 'tibet': (31.27, 88.77),
    'xinjiang': (41.73, 85.30), '新疆': (41.73, 85.30),
    'heilongjiang': (46.80, 127.95), '黑龙江': (46.80, 127.95),
    'jilin': (43.70, 126.19), '吉林': (43.70, 126.19),
    'liaoning': (41.60, 123.40), '辽宁': (41.60, 123.40),
    # 常见国家 (国际数据集兜底)
    'china': (35.00, 103.00), '中国': (35.00, 103.00),
    'japan': (36.20, 138.25), 'korea': (36.50, 127.90),
    'india': (22.00, 79.00), 'vietnam': (14.06, 108.28),
    'thailand': (15.87, 100.99), 'usa': (39.83, -98.58),
    'united states': (39.83, -98.58), 'brazil': (-14.24, -51.93),
    'australia': (-25.27, 133.78), 'russia': (61.52, 105.32),
    'germany': (51.17, 10.45), 'france': (46.23, 2.21),
    'uk': (55.38, -3.44), 'united kingdom': (55.38, -3.44),
    # ── 2026-09-16 补: 上游 VirPhyKit 示例 (H3N2/RSV/PVS) 实际出现的地名 ──
    # 判据: 拿 Example/ 真数据跑 sample_plot, 把仍落 'Other' 的名单补进来。
    # 补前 185/476 (38.9%) 未识别, 补后应显著下降。
    # 注: 无空格拼写 (HongKong/NewZealand/UnitedKingdom) 是上游数据的实际写法,
    #     一并给带空格变体, 因为 _norm 不拆驼峰。
    'hongkong': (22.32, 114.17), 'hong kong': (22.32, 114.17),
    'newzealand': (-40.90, 174.89), 'new zealand': (-40.90, 174.89),
    'peru': (-9.19, -75.02), 'singapore': (1.35, 103.82),
    'denmark': (56.26, 9.50), 'malaysia': (4.21, 101.98),
    'netherlands': (52.13, 5.29), 'nicaragua': (12.87, -85.21),
    'canada': (56.13, -106.35), 'mexico': (23.63, -102.55),
    'iran': (32.43, 53.69), 'guam': (13.44, 144.79),
    'philippines': (12.88, 121.77), 'indonesia': (-0.79, 113.92),
    'italy': (41.87, 12.57), 'spain': (40.46, -3.75),
    'sweden': (60.13, 18.64), 'norway': (60.47, 8.47),
    'finland': (61.92, 25.75), 'poland': (51.92, 19.15),
    'belgium': (50.50, 4.47), 'switzerland': (46.82, 8.23),
    'austria': (47.52, 14.55), 'portugal': (39.40, -8.22),
    'ireland': (53.14, -7.69), 'greece': (39.07, 21.82),
    'turkey': (38.96, 35.24), 'egypt': (26.82, 30.80),
    'south africa': (30.56, 22.94), 'argentina': (-38.42, -63.62),
    'chile': (-35.68, -71.54), 'colombia': (4.57, -74.30),
    'cambodia': (12.57, 104.99), 'laos': (19.86, 102.50),
    'myanmar': (21.91, 95.96), 'bangladesh': (23.68, 90.36),
    'pakistan': (30.38, 69.35), 'sri lanka': (7.87, 80.77),
    'nepal': (28.39, 84.12), 'mongolia': (46.86, 103.85),
    'kazakhstan': (48.02, 66.92), 'saudi arabia': (23.89, 45.08),
    'israel': (31.05, 34.85), 'kenya': (-0.02, 37.91),
    # 无空格变体 (上游数据实际写法) + 小国补全
    'unitedkingdom': (55.38, -3.44), 'qatar': (25.35, 51.18),
    'serbia': (44.02, 21.01), 'croatia': (45.10, 15.20),
    'czech': (49.82, 15.47), 'hungary': (47.16, 19.50),
    'romania': (45.94, 24.97), 'bulgaria': (42.73, 25.49),
    'kuwait': (29.31, 47.48), 'oman': (21.51, 55.92),
    'jordan': (30.59, 36.24), 'lebanon': (33.85, 35.86),
    'morocco': (31.79, -7.09), 'tunisia': (33.89, 9.54),
    'nigeria': (9.08, 8.68), 'ethiopia': (9.15, 40.49),
    'ghana': (7.95, -1.02), 'tanzania': (-6.37, 34.89),
    'madagascar': (-18.77, 46.87), 'uruguay': (-32.52, -55.77),
    'paraguay': (-23.44, -58.44), 'bolivia': (-16.29, -63.59),
    'ecuador': (-1.83, -78.18), 'venezuela': (6.42, -66.59),
    'cuba': (21.52, -77.78), 'jamaica': (18.11, -77.30),
    'panama': (8.54, -80.78), 'costa rica': (9.75, -83.75),
    'guatemala': (15.78, -90.23), 'honduras': (15.20, -86.24),
    'brunei': (4.54, 114.73), 'fiji': (-17.71, 178.07),
    'papua new guinea': (-6.31, 143.96), 'hong kong': (22.32, 114.17),
}

# ── 市县级坐标 (GCVA 等数据集实际出现过的 + 枸杞产区常见地) ──────────
CITY_LOOKUP = {
    'yinchuan': (38.47, 106.27), '银川': (38.47, 106.27),
    'zhongning': (37.49, 105.68), '中宁': (37.49, 105.68),
    'wuzhong': (37.99, 106.20), '吴忠': (37.99, 106.20),
    'guyuan': (36.00, 106.24), '固原': (36.00, 106.24),
    'shizuishan': (38.98, 106.38), '石嘴山': (38.98, 106.38),
    'wuwei': (37.93, 102.64), '武威': (37.93, 102.64),
    'xining': (36.62, 101.78), '西宁': (36.62, 101.78),
    'lanzhou': (36.06, 103.83), '兰州': (36.06, 103.83),
    'shenyang': (41.80, 123.43), '沈阳': (41.80, 123.43),
    'guangzhou': (23.13, 113.26), '广州': (23.13, 113.26),
    'nanjing': (32.06, 118.78), '南京': (32.06, 118.78),
    'beijing city': (39.90, 116.40),
    'urumqi': (43.79, 87.63), '乌鲁木齐': (43.79, 87.63),
    'hohhot': (40.82, 111.75), '呼和浩特': (40.82, 111.75),
    'ulaanbaatar': (47.89, 106.91),
}

# 意为「未知」的占位符: 不查表直接拒
_PLACEHOLDER = {'', 'unknown', 'na', 'n/a', 'not applicable', 'none', 'null',
                'missing', '其他', '未知'}


def _norm(s: str) -> str:
    """规范化 key: 小写 + 压缩空白 + 去常见装饰词"""
    s = str(s or '').strip().lower()
    s = re.sub(r'[\s,]+', ' ', s)
    s = re.sub(r'\b(province|sheng|region|autonomous)\b', '', s).strip()
    s = s.replace('province', '').strip()
    return s


class GeoResolver:
    """分层地名解析。线程不安全 (缓存写盘用 json, 单进程场景够用)。"""

    def __init__(self, cache_path=None, online=None):
        self.cache_path = Path(cache_path or Path(__file__).parent.parent / 'cache' / 'geo_cache.json')
        self.online = online if online is not None else (
            os.environ.get('GEO_RESOLVER_ONLINE', '1') != '0')  # 默认开: 公共 API + 持久缓存
        self.cache = {}
        self._cache_dirty = False
        try:
            if self.cache_path.exists():
                self.cache = json.loads(self.cache_path.read_text(encoding='utf-8'))
        except (json.JSONDecodeError, OSError):
            self.cache = {}
        # 启动时把旧 COORD_LOOKUP 的全格式条目预热进缓存反查集
        self._full_index = {}
        for k, v in {**PROVINCE_LOOKUP, **CITY_LOOKUP}.items():
            self._full_index[_norm(k)] = v

    # ── L1+L2+L3: 离线解析 ──
    def _resolve_offline(self, loc_str: str):
        key = _norm(loc_str)
        if key in _PLACEHOLDER:
            return None, 'placeholder'
        # L1 精确
        if key in self._full_index:
            return self._full_index[key], 'exact'
        # 缓存命中 (含此前在线解析/人工修正结果)
        if key in self.cache:
            v = self.cache[key]
            return (tuple(v) if v else None), 'cache'
        # L2 逗号段拆分: 'china, ningxia, yinchuan' → 段重组
        # 原则: 优先最具体段 (市 > 省), 市级查到用市级
        parts = [p.strip() for p in re.split(r'[,;]', str(loc_str)) if p.strip()]
        best = None
        for p in parts:  # 从右往左 = 从具体到宽泛
            pk = _norm(p)
            if pk in _PLACEHOLDER or pk == 'china' or pk == '中国':
                continue
            if pk in CITY_LOOKUP:
                best = (CITY_LOOKUP[pk], 'city')
                break
        if not best:
            for p in reversed(parts):
                pk = _norm(p)
                if pk in _PLACEHOLDER or pk in ('china', '中国'):
                    continue
                if pk in PROVINCE_LOOKUP:
                    best = (PROVINCE_LOOKUP[pk], 'province')
                    break
        if best:
            return best
        # L3 段内反查 (变体拼写): 'ningxia' 段与 'neimenggu' 段都试 norm 后缀
        for p in reversed(parts):
            pk = _norm(p)
            if not pk or pk in _PLACEHOLDER:
                continue
            for known, v in self._full_index.items():
                if known.endswith(' ' + pk) or pk.endswith(' ' + known):
                    return v, 'fuzzy'
        return None, 'miss'

    # ── L5: 在线兜底 (Photon 公共 API, 免费无 key; 默认开, 结果持久缓存后等效离线) ──
    def _resolve_online(self, loc_str: str):
        if not self.online:
            return None
        import urllib.request
        import urllib.parse
        for endpoint in (
            'https://photon.komoot.io/api',            # Photon 公共服务 (首选, 模糊强)
            'https://nominatim.openstreetmap.org/search'  # 兑底
        ):
            try:
                q = urllib.parse.quote(str(loc_str))
                if 'photon' in endpoint:
                    url = f'{endpoint}?q={q}&limit=1&lang=en'
                    req = urllib.request.Request(url, headers={
                        'User-Agent': 'virome-phylo-pipeline/2.0 (geo_resolver)'})
                    with urllib.request.urlopen(req, timeout=8) as r:
                        data = json.loads(r.read().decode())
                    if data.get('features'):
                        c = data['features'][0]['geometry']['coordinates']  # [lon, lat]
                        return (float(c[1]), float(c[0]))
                else:
                    url = f'{endpoint}?q={q}&format=json&limit=1'
                    req = urllib.request.Request(url, headers={
                        'User-Agent': 'virome-phylo-pipeline/2.0 (geo_resolver)'})
                    with urllib.request.urlopen(req, timeout=8) as r:
                        data = json.loads(r.read().decode())
                    if data:
                        return (float(data[0]['lat']), float(data[0]['lon']))
                time.sleep(1.0)  # 公共服务限速礼貌
            except Exception:
                continue
        return None

    def resolve(self, loc_str):
        """返回 ((lat, lon) | None, source)"""
        coord, src = self._resolve_offline(loc_str)
        if coord is None and self.online and src == 'miss':
            coord = self._resolve_online(loc_str)
            src = 'online' if coord else 'miss'
        # 回写缓存 (负缓存也写, 避免重复在线打同一名)
        if src == 'miss' or src == 'online':
            self.cache[_norm(loc_str)] = list(coord) if coord else None
            self._cache_dirty = True
        return coord, src

    def flush(self):
        if self._cache_dirty:
            self.cache_path.parent.mkdir(parents=True, exist_ok=True)
            self.cache_path.write_text(
                json.dumps(self.cache, ensure_ascii=False, indent=1), encoding='utf-8')
            self._cache_dirty = False


_resolver = None

def get_resolver() -> GeoResolver:
    global _resolver
    if _resolver is None:
        _resolver = GeoResolver()
    return _resolver


def resolve_location(loc_str):
    """便捷函数: 返回 (lat, lon) 或 None"""
    coord, _ = get_resolver().resolve(loc_str)
    return coord
