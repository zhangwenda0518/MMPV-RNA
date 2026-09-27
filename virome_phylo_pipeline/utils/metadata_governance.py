#!/usr/bin/env python3
"""metadata_governance.py — metadata 入口治理 (时间 + 地理 的检查与矫正)

2026-09-16 建立
---------------
**背景**. 管线入口的 metadata 来自多个上游 (GenBank 抓取 / 公共库合并 / 用户自备),
字段格式不统一。历史上已发生过的「静默出错」共四类 (均为全量扫描实测, 见
`_consistency_check_20260916/diagnose_metadata_20260916.py`):

  A. 坐标占位符     `XX.XX N XXX.XX E` 出现 450/450 —— 坐标列未填,
                    **但 408 条本可由地名推算** (geo_to_latlon / gazetteer);
                    2026-09-16 老师拍板改为「先补算, 再清空」
  B. 日期占位符     `collection_date` 字面量 `YYYY-MM-DD` 42 条 —— 形似合法日期,
                    任何"格式校验"都会放行
  C. 地理层级 Unknown  `China, Unknown, Yinchuan_AI` —— geo_resolver 会 skip
                    `Unknown` 段, 静默降级成 China 国家质心 (35,103), 误差上千公里
  D. 只给年 / GenBank 日期  39 条 `YYYY` + 36 条 `DD-Mon-YYYY`; 后者在
                    virphy_bridge 侧被 `dropna` 静默剔除

**本模块的定位**. 只做「归一 + 标注 + 体检」, **不改变**任何现有判定阈值;
在流程入口产出标准化 metadata, 作为下游唯一输入。原有模块 (geo_resolver /
decimal_year / virphy_bridge) 一律只读不改。

**三条原则**:
  1. 不猜   —— 解析不出就标 none + 记 issue, 绝不套默认值
  2. 不静默 —— 任何降精度 / 任何丢弃都进报告, 可统计条数
  3. 不改语义 —— 只归一 + 标注; 口径由调用方显式指定

**老师拍板 (2026-09-16)**:
  · 模块范围: 体检 + 格式化, **另存新文件** (原件不动)
  · 地理 Unknown 层级: **删除**该层级
  · virphy_bridge 月初口径: **接入统一入口**, 用 `mode='start'`
"""
from __future__ import annotations

import calendar
import csv
import io
import json
import pathlib
import re
import unicodedata
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

__all__ = [
    "DateInfo", "parse_date", "decimal_year", "DATE_FORMAT_TABLE",
    "split_location", "is_placeholder", "GeoInfo", "parse_location",
    "coord_from_location", "govern_table", "write_outputs",
    "has_host", "HOST_PLACEHOLDER_TOKENS", "normalize_host",
]

# ══════════════════════════════════════════════════════════════════
# 宿主 (host) —— 「有没有宿主」的唯一判据
# ══════════════════════════════════════════════════════════════════
# 2026-09-16 第三轮 (老师: 「host 也独立出来, 仅对有 host 的进行分析」)。
#
# host 从「地理通道的附属列」独立成**第三张通道表** (`data/host.csv`),
# 且该表**只收有宿主的样品** —— 于是宿主分析天然只跑在有宿主的子集上。
# 「有没有宿主」这个判定必须**全局唯一**, 否则「写表时算一次、分析时又算一次」
# 必然漂移 (本次修复的 C7 幂等 bug 就是同款问题)。
#
# 哪些值算「没有宿主」: 上游未填的占位字面量 + NCBI 常见的**无信息宿主**取值。
# 后者尤其重要 —— `uncultured bacterium` / `metagenome` 这类字符串看着像宿主,
# 放进宿主分化分析会凭空造出一个"物种组", 使组间/组内距离失去意义。
HOST_PLACEHOLDER_TOKENS = frozenset({
    "", "na", "n/a", "nan", "none", "null", "missing", "unknown", "not applicable",
    "not_provided", "not provided", "tbd", "未", "未知", "待定",
    # NCBI 常见「无信息宿主」—— 不是宿主物种, 是"宿主未知"的另一种写法
    "unidentified", "unidentified organism", "uncultured", "uncultured bacterium",
    "unclassified", "other", "synthetic construct", "metagenome",
    "mixed culture", "environmental sample", "laboratory culture",
})


def normalize_host(value) -> str:
    """宿主名归一: 去首尾空白 + 压平内部空白 + 去尾部分隔符。"""
    return re.sub(r"\s+", " ", str(value or "").strip()).strip(" ,;|/")


def has_host(value) -> bool:
    """该样品的 host 是否是**可用于宿主分化分析**的真实宿主。

    返回 False 的取值分两类, 都不应进入 `data/host.csv`:
      · 未填: `''` / `Unknown` / `Not_Provided` / `NA` …
      · 无信息: `uncultured bacterium` / `metagenome` / `synthetic construct` …
        —— 它们会在地理/宿主维度里造出虚假的"物种组"。
    """
    h = normalize_host(value)
    if not h:
        return False
    if h.lower() in HOST_PLACEHOLDER_TOKENS:
        return False
    if is_placeholder(h):          # 复用日期/地理那套字面量与形似格式拦截
        return False
    return True


# ══════════════════════════════════════════════════════════════════
# 占位符 (上游未填值的字面量) —— 一律拦截, 不进入下游
# ══════════════════════════════════════════════════════════════════
PLACEHOLDER_LITERALS = {
    "", "na", "n/a", "nan", "none", "null", "missing", "unknown", "not applicable",
    "not_provided", "tbd", "未", "未知", "待定",
}
# 形似格式的占位符 (最危险: 通过一切"格式校验")
# ⚠ 教训 (2026-09-16): 曾用泛化正则 `^[a-z]+\s*[:：]\s*[a-z]+$` 想兜住任意
#   "Word:Word" 占位符, 结果**误伤真实地理** `China:Ningxia` / `China:Beijing`
#   (72 条全被打成 placeholder)。改为只精确匹配已知占位字面量。
PLACEHOLDER_PATTERNS = [
    re.compile(r"^y{4}[-/.]m{2}[-/.]d{2}$", re.I),   # YYYY-MM-DD
    re.compile(r"^y{4}[-/.]m{2}$", re.I),            # YYYY-MM
    re.compile(r"^y{4}$", re.I),                     # YYYY
    re.compile(r"^x+(\.[x]+)?\s*[ns]\s*x+(\.[x]+)?\s*[ew]$", re.I),  # XX.XX N XXX.XX E
    re.compile(r"^country\s*[:：]\s*region$", re.I),  # Country:Region (仅这一组)
]

_MONTHS = {
    "jan": 1, "feb": 2, "mar": 3, "apr": 4, "may": 5, "jun": 6,
    "jul": 7, "aug": 8, "sep": 9, "sept": 9, "oct": 10, "nov": 11, "dec": 12,
}


def is_placeholder(v) -> bool:
    """判断是否为「未填值」占位符 (含形似格式的假日期/假地理)。"""
    s = str(v or "").strip()
    if s.lower() in PLACEHOLDER_LITERALS:
        return True
    for p in PLACEHOLDER_PATTERNS:
        if p.match(s):
            return True
    return False


# ══════════════════════════════════════════════════════════════════
# 日期
# ══════════════════════════════════════════════════════════════════

# 显式白名单登记表: (名称, 说明, 归一化后 precision)
DATE_FORMAT_TABLE = [
    ("iso_day",      "YYYY-MM-DD",   "day"),
    ("iso_month",    "YYYY-MM",      "month"),
    ("iso_year",     "YYYY",         "year"),
    ("genbank",      "DD-Mon-YYYY",  "day"),
    ("euro_day",     "DD/MM/YYYY 或 DD.MM.YYYY", "day"),
    ("mon_year",     "Mon-YYYY",     "month"),
    ("decimal_year", "已转好的小数年 (如 2007.569473)", "day*"),
]


@dataclass
class DateInfo:
    """一次日期解析的完整结果 (结构化, 不是裸 float)。"""
    raw: str = ""
    iso: Optional[str] = None          # 归一化后 ISO: YYYY-MM-DD / YYYY-MM / YYYY
    decimal: Optional[float] = None    # 按指定 mode 算出的小数年
    precision: str = "none"            # day | month | year | none
    confidence: float = 0.0            # 1.0=精确命中白名单; <1 表示推断
    source_format: str = "unknown"     # 命中的 DATE_FORMAT_TABLE 名称
    issue: Optional[str] = None        # placeholder | out_of_range | unparsable | ambiguous
    consumed_as: str = "fresh"         # fresh=从字符串解析; already_decimal=输入本就是小数年

    def to_row(self) -> Dict[str, str]:
        d = asdict(self)
        d["decimal"] = "" if self.decimal is None else f"{self.decimal:.10f}"
        d["raw"] = self.raw
        return {k: ("" if v is None else str(v)) for k, v in d.items()}


def _split_seps(s: str) -> List[str]:
    """吸收 BioAider 的唯一优点: `.` `/` `_` `-` 四种分隔符统一处理。

    与 BioAider 的 `re.split("[./_-]", s)` 等价, 但只在**归一化阶段**用,
    解析与判定仍走我们自己的严格逻辑 (不继承它的静默回显 / 越界强解)。
    """
    return [p for p in re.split(r"[./_\-]", s) if p != ""]


def _classify_iso(parts: List[str], raw: str) -> Optional[Tuple[str, str, str]]:
    """判断是否 YYYY 在前。返回 (iso, precision, fmt)，否则 None。"""
    try:
        y = int(parts[0])
    except ValueError:
        return None
    if not (1900 <= y <= 2100):
        return None
    if len(parts) == 1 and len(parts[0]) == 4:
        return f"{y:04d}", "year", "iso_year"
    if len(parts) == 2:
        m = int(parts[1])
        if 1 <= m <= 12 and len(parts[0]) == 4:
            return f"{y:04d}-{m:02d}", "month", "iso_month"
        return None
    if len(parts) == 3 and len(parts[0]) == 4:
        m, d = int(parts[1]), int(parts[2])
        if 1 <= m <= 12 and 1 <= d <= 31:
            dim = calendar.monthrange(y, m)[1]
            if d > dim:
                return None
            return f"{y:04d}-{m:02d}-{d:02d}", "day", "iso_day"
        return None
    return None


def _classify_euro(parts: List[str]) -> Optional[Tuple[str, str, str, str]]:
    """判断是否 DD 在前 (DD/MM/YYYY 或 DD.MM.YYYY)。返回 (iso, precision, fmt, note)。"""
    if len(parts) != 3:
        return None
    try:
        a, b, y = int(parts[0]), int(parts[1]), int(parts[2])
    except ValueError:
        return None
    if len(parts[2]) != 4 or not (1900 <= y <= 2100):
        return None
    if not (1 <= a <= 31 and 1 <= b <= 12):
        return None
    dim = calendar.monthrange(y, b)[1]
    if a > dim:
        return None
    note = "ambiguous" if a <= 12 else "clear"   # a<=12 时 DD 与 MM 都可解释
    return f"{y:04d}-{b:02d}-{a:02d}", "day", "euro_day", note


def _classify_genbank(s: str) -> Optional[Tuple[str, str, str]]:
    """DD-Mon-YYYY / Mon-YYYY / Mon-YYYY (GenBank 标准)。"""
    m = re.match(r"^(\d{1,2})[-\s]?([A-Za-z]{3,4})[-\s]?(\d{4})$", s)
    if m:
        mo = _MONTHS.get(m.group(2).lower())
        y, d = int(m.group(3)), int(m.group(1))
        if mo and 1900 <= y <= 2100 and 1 <= d <= calendar.monthrange(y, mo)[1]:
            return f"{y:04d}-{mo:02d}-{d:02d}", "day", "genbank"
    m = re.match(r"^([A-Za-z]{3,4})[-\s]?(\d{4})$", s)
    if m:
        mo = _MONTHS.get(m.group(1).lower())
        y = int(m.group(2))
        if mo and 1900 <= y <= 2100:
            return f"{y:04d}-{mo:02d}", "month", "mon_year"
    return None


def parse_date(raw, *, mode: str = "mid") -> DateInfo:
    """解析日期字符串 → DateInfo。

    mode:
      'mid'   → 只给月/只给年时取该月 15 日 / 6 月 15 日 (utils/decimal_year 现口径)
      'start' → 只给月/只给年时取月初 / 年初       (virphy_bridge 现口径)

    ⚠ 'start' 与 'mid' 在**完整日期 (YYYY-MM-DD)** 上结果相同, 只在
      「只给年 / 只给年月」两级有差异 —— 这正是上游 VirPhyKit 与本管线
      的**有意差异**, 由调用方显式选择, 不藏在实现里。
    """
    assert mode in ("mid", "start"), f"mode must be 'mid'|'start', got {mode!r}"
    info = DateInfo(raw=str(raw or "").strip())

    if is_placeholder(info.raw):
        info.issue = "placeholder"
        return info

    s = info.raw

    # ① 已是小数年 (VirPhyKit Example 的 date 列) —— 必须识别, 否则会重复转换。
    #    消歧判据 (实测): 真实小数年的**小数位数 ≥ 4** (VirPhyKit h3n2_na_500
    #    实测 431 条 6 位 / 42 条 5 位 / 3 条 4 位); 而 `年.月` 写法的小数位
    #    **恒为 2** (`2014.05` / `2022.11`, BioAider GUI 示例)。故 2 位让给
    #    「年月」分支, ≥4 位才判小数年, 两者不冲突。
    m_frac = re.fullmatch(r"(\d{4})\.(\d{4,})", s)
    if m_frac:
        try:
            v = float(s)
        except ValueError:
            info.issue = "unparsable"
            return info
        y = int(v)
        if 1900 <= y <= 2100:
            import datetime as _dt
            doy = max(1, min(366, round((v - y) * 365) + 1))
            base = _dt.date(y, 1, 1) + _dt.timedelta(days=doy - 1)
            info.iso = base.isoformat()
            info.decimal = v
            info.precision = "day"
            info.confidence = 0.9
            info.source_format = "decimal_year"
            info.consumed_as = "already_decimal"
            return info
        info.issue = "out_of_range"
        return info

    # ② GenBank 文字月 (必须先于分隔符归一, 否则 'Jun' 会被拆散)
    gb = _classify_genbank(s)
    if gb:
        info.iso, info.precision, info.source_format = gb[0], gb[1], gb[2]
        info.confidence = 1.0
        info.decimal = _iso_to_decimal(info.iso, mode)
        return info

    parts = _split_seps(s)
    if not parts:
        info.issue = "unparsable"
        return info

    # ③ YYYY 在前
    r = _classify_iso(parts, s)
    if r:
        info.iso, info.precision, info.source_format = r
        info.confidence = 1.0
        info.decimal = _iso_to_decimal(info.iso, mode)
        return info

    # ④ DD 在前 (仅 3 段), 记录歧义
    r = _classify_euro(parts)
    if r:
        info.iso, info.precision, info.source_format, note = r
        info.confidence = 0.7 if note == "ambiguous" else 1.0
        if note == "ambiguous":
            info.issue = "ambiguous"
        info.decimal = _iso_to_decimal(info.iso, mode)
        return info

    # ⑤ 4 段及以上 / 数字越界
    nums = [p for p in parts if p.isdigit()]
    if nums and any(len(p) == 4 and p.startswith(("19", "20")) for p in nums):
        info.issue = "out_of_range"
    else:
        info.issue = "unparsable"
    return info


def _iso_to_decimal(iso: str, mode: str) -> float:
    """ISO → 小数年。mode='mid' 与 utils/decimal_year 逐位一致;
    mode='start' 与 virphy_bridge._to_decimal 逐位一致 (只给年→年初, 只给月→月初)。"""
    p = iso.split("-")
    y = int(p[0])
    if len(p) == 1:
        # 只给年
        return float(y) if mode == "start" else y + (6.0 - 1.0 + (15.0 - 1.0) / 30.0) / 12.0
    mo = int(p[1])
    if len(p) == 2:
        # 只给年月
        if mode == "start":
            return float(y) + (float(mo) - 1) / 12.0
        dim = calendar.monthrange(y, mo)[1]
        return y + (mo - 1.0 + (15.0 - 1.0) / dim) / 12.0
    dy = int(p[2])
    if mode == "start":
        # ⚠ virphy_bridge 的 3 段走的是 utils/decimal_year (统一口径), 非"月初"
        return _mid_decimal(y, mo, dy)
    return _mid_decimal(y, mo, dy)


def _mid_decimal(y: int, mo: int, dy: int) -> float:
    """utils/decimal_year.to_decimal_year 的月内插值公式 (逐位复刻)。"""
    dim = calendar.monthrange(y, mo)[1]
    d = min(dy, dim)
    return y + (mo - 1.0 + (d - 1.0) / dim) / 12.0


def decimal_year(raw, *, mode: str = "mid") -> Optional[float]:
    """便捷入口: 直接拿小数年 (拿不到返回 None, 绝不编造)。"""
    return parse_date(raw, mode=mode).decimal


# ══════════════════════════════════════════════════════════════════
# 地理
# ══════════════════════════════════════════════════════════════════

_UNKNOWN_TOKENS = {"unknown", "unk", "未知", "na", "n/a", "none", "null", ""}


def _is_unknown_level(p: str) -> bool:
    """该层级是否为「未知」。

    判定: 完全等于占位词, 或**以占位词开头** (处理 `Unknown_AI` / `Unknown-Region`
    这类「占位词 + 来源标记」的写法)。注意排除 `unknown` 出现在中间的真实地名。
    """
    t = p.strip().lower()
    if t in _UNKNOWN_TOKENS:
        return True
    for k in _UNKNOWN_TOKENS:
        if k and (t.startswith(k + "_") or t.startswith(k + "-") or t.startswith(k + " ")):
            return True
    return False


def _is_country_only_term(t: str) -> bool:
    """该层级是否为「纯国家名」—— 即**比国家更粗、没有地理分辨率**的层级。

    用于判断「删掉 Unknown 层级后是否还剩可用的地理信息」:
      · `China` / `巴西` / `United States` → True (只有国家级)
      · `Taiwan, China` / `Macao` / `香港`   → False (在中国主权下是**一级行政区**,
        属「比国家更细的层级」, 有独立坐标)

    ⚠ 必须排除 `_CN_SUBREGION_KEYS`(台港澳): 它们为消歧被登记进了
      `_COUNTRY_NAME_TO_ISO2`, 但**不是纯国家名**。若一并判为 True,
      `Macao, Unknown` 会被当成"无地理信息"整条作废, 而它本可拿到
      `22.13 N 113.556 E`。
    """
    k = str(t or "").strip().lower()
    return k in _COUNTRY_NAME_TO_ISO2 and k not in _CN_SUBREGION_KEYS


@dataclass
class GeoInfo:
    raw: str = ""
    normalized: Optional[str] = None   # 归一化后, 如 "China, Ningxia, Yinchuan"
    components: List[str] = field(default_factory=list)  # 保留的层级
    dropped_levels: List[str] = field(default_factory=list)  # **被删除**的 Unknown 层级
    unresolved: bool = False           # 是否含未识别层级
    issue: Optional[str] = None        # placeholder | unknown_level | country_only | unparsable

    @property
    def usable(self) -> bool:
        """本记录是否**还有可用于下游分析的地理信息**。

        与 `components` 的区别 (`country_only` 一类):
          · `China, Unknown, Unknown_AI` —— `China` 确实是「保留的层级」(不是占位符),
            故 `components == ['China']` 是**如实**的; 但删掉占位层级后只剩国家名,
            采样省/市本来就是 Unknown, 该记录**没有任何地理分辨率**可用 ——
            沿用旧实现赋"中国中心点"属**编造数据** (误差可达上千公里),
            还会撞幂等 (产物写成 `China` 后再治理 → 命中中心点, 两次坐标不同)。

        故下游判"有没有地理"必须用 `usable` (→ `govern_table` 的 `no_loc`),
        不能再用 `not components`, 否则这类记录会带着空坐标存活到系统地理分析里。

        为何不直接清空 `components`: 那会让报告里的 `geo_components` 丢失
        「国家其实是知道的」这一事实, 排障时看不出被判废的原因。
        """
        return bool(self.components) and self.issue != "country_only"

    def to_row(self) -> Dict[str, str]:
        return {
            "geo_raw": self.raw,
            "geo_normalized": self.normalized or "",
            "geo_components": "|".join(self.components),
            "geo_dropped": "|".join(self.dropped_levels),
            "geo_unresolved": str(self.unresolved),
            "geo_issue": self.issue or "",
        }


def split_location(raw) -> List[str]:
    """拆地理层级。**冒号与逗号都拆** (修 geo_resolver._norm 只拆 [,;] 的缺口):
        'China:Ningxia'            → ['China', 'Ningxia']
        'China, Ningxia, Yinchuan_AI' → ['China', 'Ningxia', 'Yinchuan_AI']
    """
    return [p.strip() for p in re.split(r"[,;:：]", str(raw or "")) if p.strip()]


def parse_location(raw, *, drop_unknown: bool = True) -> GeoInfo:
    """解析地理字符串。

    drop_unknown=True (老师拍板): **删除** Unknown 层级, 并在 dropped_levels
    中留痕。这样 `China, Unknown, Yinchuan_AI` → `China, Yinchuan`, 不会像
    geo_resolver 那样静默降级成 China 国家质心却仍标记为"已解析"。

    2026-09-16 第三轮 (老师: 「Unknown:Unknown 这些在清理过程应该去掉」):
    删掉占位层级后**只剩纯国家名**时 (`China, Unknown, Unknown_AI`),
    该记录已无任何地理分辨率 —— issue 标 `country_only`, `normalized=None`,
    由 `govern_table` 清空地名格并按 `no_loc` 走 `--drop-no-location` 剔除。
    判据见 `GeoInfo.usable`。
    """
    info = GeoInfo(raw=str(raw or "").strip())
    if is_placeholder(info.raw):
        info.issue = "placeholder"
        info.unresolved = True
        return info

    parts = split_location(info.raw)
    kept, dropped = [], []
    for p in parts:
        if _is_unknown_level(p):
            if drop_unknown:
                dropped.append(p)
            else:
                kept.append(p)
        else:
            kept.append(p)

    info.components = kept
    info.dropped_levels = dropped
    if dropped:
        info.issue = "unknown_level"
        info.unresolved = True

    # 只剩纯国家名 + 占位层级被删 → 无可用地理信息 (不得赋国家中心点)
    if dropped and kept and all(_is_country_only_term(p) for p in kept):
        info.issue = "country_only"
        info.normalized = None
        info.unresolved = True
        return info

    info.normalized = ", ".join(kept) if kept else (None if drop_unknown else ", ".join(parts))
    if info.normalized is None and not dropped:
        info.issue = info.issue or "unparsable"
        info.unresolved = True
    return info


# ══════════════════════════════════════════════════════════════════
# 表级治理
# ══════════════════════════════════════════════════════════════════

# `run` 是 Core14 (`Global_Unified_Metadata_Core14.tsv`) 的样品主键列名 ——
# 缺了它, 批处理入口的治理产物就没有 name 列, 下游"按名回填两张通道表"直接失败
# (`来源表无 name 列`)。放**末位**: _pick_col 是 精确>前缀>包含 且按 hints 顺序,
# 放最后可确保既有 name/tip/accession 等更明确的列名优先命中。
NAME_HINTS = ("name", "tip", "sequence_name", "seq", "id", "accession", "isolate",
              "strain", "run")
DATE_HINTS = ("collection_date", "collectiondate", "date", "sampling_date", "time", "year")
LOC_HINTS = ("location", "geo_loc_name", "geolocname", "region", "locality",
             "country", "state", "province", "place", "geo_loc")
COORD_HINTS = ("lat_lon", "latlon", "lat-lon", "coordinates", "coords", "latitude", "longitude")


# ══════════════════════════════════════════════════════════════════
# 地名 → 坐标 推算 (2026-09-16 老师拍板: 坐标占位符「先补算, 再清空」)
# ══════════════════════════════════════════════════════════════════
# 2026-09-16 第二轮 (老师: 「不能只是中国的省市呀, 全世界的地理怎么做」):
# 由「30 个中国城市硬编码字典」升级为**多级全球查找**。
#
# 查找优先级 (高 → 低):
#   1. `_CITY_COORDS`  —— 30 个中国城市手校字典。**保留且最高优先**,
#      保证与 submission pipeline 的 `geo_to_latlon` 对拍结果一字不变。
#   2. admin1 省/州级 —— 全球 4410 条 (Natural Earth 1:10m)
#   3. admin0 国家级 —— 全球 255 条 (Natural Earth 1:10m)
#   4. 推不出 → 返回 '' , 由调用方决定清空
#
# 数据来源: `_geo_build_20260916/dist/gazetteer.json`
#   由 Natural Earth 1:10m (public domain) 用 shapely `representative_point()`
#   构建, 脚本 `_geo_build_20260916/scripts/build_gazetteer.py` 可复跑。
#   参考实现: 植物病毒分析平台 git-repo/spreadgl2.github.io (MIT) 的
#   `src/lib/format/gazetteer.ts` 及其公开构建仓库 spreadgl2/spreadgl2-gazetteer。
#
# ⚠ 中国主权相关: 上游 Natural Earth 把台湾标为 `Sovereign country`、
#   中文名写作「中华民国」, 与中华人民共和国官方立场不一致。构建脚本已在
#   生成时规范化: 台湾/香港/澳门 一律归入 `countryIso2='CN'`,
#   规范名 `Taiwan, China` / `Hong Kong, China` / `Macao, China`,
#   中文名 `中国台湾` / `中国香港` / `中国澳门`。
_CITY_COORDS = {
    'yinchuan':     '38.47 N 106.27 E',
    'ningxia':      '37.48 N 105.68 E',
    'zhongning':    '37.48 N 105.68 E',
    'beijing':      '39.90 N 116.40 E',
    'shanghai':     '31.23 N 121.47 E',
    'guangzhou':    '23.13 N 113.26 E',
    'wuhan':        '30.59 N 114.31 E',
    'nanjing':      '32.06 N 118.79 E',
    'hangzhou':     '30.27 N 120.15 E',
    'chengdu':      '30.57 N 104.07 E',
    'xian':         '34.26 N 108.94 E',
    'kunming':      '25.04 N 102.68 E',
    'harbin':       '45.80 N 126.53 E',
    'zhengzhou':    '34.75 N 113.62 E',
    'jinan':        '36.65 N 116.98 E',
    'taiyuan':      '37.87 N 112.55 E',
    'changsha':     '28.23 N 112.94 E',
    'fuzhou':       '26.07 N 119.30 E',
    'guiyang':      '26.65 N 106.63 E',
    'lanzhou':      '36.06 N 103.79 E',
    'xining':       '36.62 N 101.77 E',
    'urumqi':       '43.79 N 87.58 E',
    'lhasa':        '29.65 N 91.10 E',
    'shenyang':     '41.80 N 123.43 E',
    'dalian':       '38.91 N 121.61 E',
    'qingdao':      '36.07 N 120.38 E',
    'suzhou':       '31.30 N 120.62 E',
    'shenzhen':     '22.54 N 114.06 E',
    'chongqing':    '29.56 N 106.55 E',
    'tianjin':      '39.13 N 117.18 E',
}
_COUNTRY_CENTROID_CN = '35.86 N 104.19 E'

# ── 全球 gazetteer 数据路径 ──
# 本文件位于 virome_phylo_pipeline/utils/ 下:
#   parents[0] = utils/                  parents[1] = virome_phylo_pipeline/
#   parents[2] = MMPV-RNA/  (项目根)     → 根/_geo_build_20260916/dist/gazetteer.json
_GAZETTEER_REL = ("_geo_build_20260916", "dist", "gazetteer.json")
_GAZETTEER_PATH = pathlib.Path(__file__).resolve().parents[2].joinpath(*_GAZETTEER_REL)

# 查找键归一化用的转写表 (与构建脚本保持一致)
_TRANSLITERATION = {
    "Æ": "AE", "æ": "ae", "Œ": "OE", "œ": "oe", "Ø": "O", "ø": "o",
    "Ð": "D", "ð": "d", "Đ": "D", "đ": "d", "Þ": "Th", "þ": "th",
    "Ł": "L", "ł": "l", "Ŋ": "N", "ŋ": "n", "Ħ": "H", "ħ": "h",
    "Ə": "E", "ə": "e", "Ŧ": "T", "ŧ": "t", "ß": "ss", "ẞ": "SS",
    "İ": "I", "ı": "i", "’": "'", "‘": "'", "ʼ": "'", "ʻ": "'",
    "ʿ": "'", "“": '"', "”": '"', "–": "-", "—": "-", "−": "-",
}
_NORM_RE = re.compile(r"[^\w]+", re.UNICODE)

# 国家级别名 → ISO2, 用于把 `United States:Maryland` 这类
# `国家:一级行政区` 写法里的国家前缀解析成消歧用的 ISO2
_COUNTRY_NAME_TO_ISO2 = {
    "united states": "US", "usa": "US", "us": "US", "america": "US",
    "united states of america": "US",
    "china": "CN", "prc": "CN", "peoples republic of china": "CN",
    "中华人民共和国": "CN", "中国": "CN",
    "united kingdom": "GB", "uk": "GB", "great britain": "GB",
    "russia": "RU", "russian federation": "RU",
    "south korea": "KR", "korea": "KR", "republic of korea": "KR",
    "north korea": "KP", "japan": "JP", "india": "IN", "brazil": "BR",
    "germany": "DE", "france": "FR", "italy": "IT", "spain": "ES",
    "canada": "CA", "australia": "AU", "mexico": "MX", "argentina": "AR",
    "south africa": "ZA", "nigeria": "NG", "kenya": "KE", "egypt": "EG",
    "vietnam": "VN", "thailand": "TH", "indonesia": "ID",
    "philippines": "PH", "malaysia": "MY", "pakistan": "PK",
    "bangladesh": "BD", "turkey": "TR", "iran": "IR", "iraq": "IQ",
    "saudi arabia": "SA", "israel": "IL", "poland": "PL",
    "netherlands": "NL", "belgium": "BE", "switzerland": "CH",
    "sweden": "SE", "norway": "NO", "denmark": "DK", "finland": "FI",
    "portugal": "PT", "greece": "GR", "ireland": "IE", "austria": "AT",
    "new zealand": "NZ", "chile": "CL", "colombia": "CO", "peru": "PE",
    "venezuela": "VE", "ecuador": "EC", "bolivia": "BO", "uruguay": "UY",
    "paraguay": "PY", "cuba": "CU", "ethiopia": "ET", "ghana": "GH",
    "tanzania": "TZ", "uganda": "UG", "morocco": "MA", "algeria": "DZ",
    "tunisia": "TN", "sudan": "SD", "senegal": "SN", "zimbabwe": "ZW",
    "zambia": "ZM", "mozambique": "MZ", "angola": "AO", "cameroun": "CM",
    "cameroon": "CM", "ukraine": "UA", "romania": "RO", "hungary": "HU",
    "czech republic": "CZ", "czechia": "CZ", "slovakia": "SK",
    "bulgaria": "BG", "croatia": "HR", "serbia": "RS", "slovenia": "SI",
    "lithuania": "LT", "latvia": "LV", "estonia": "EE", "belarus": "BY",
    "kazakhstan": "KZ", "uzbekistan": "UZ", "mongolia": "MN",
    "nepal": "NP", "sri lanka": "LK", "myanmar": "MM", "cambodia": "KH",
    "laos": "LA", "singapore": "SG", "taiwan": "CN", "taiwan, china": "CN",
    "hong kong": "CN", "hong kong, china": "CN",
    "macao": "CN", "macau": "CN", "macao, china": "CN", "中国台湾": "CN",
    "中国香港": "CN", "中国澳门": "CN",
}

_GAZ = None  # (exact_index, normalized_index), 惰性加载
_GAZ_LOAD_FAILED = False

# 中国台湾/香港/澳门: 在中国主权下是**一级行政区**, 不是国家。
# 这些词条在 `_COUNTRY_NAME_TO_ISO2` 里映射到 CN 只是为了消歧提示,
# 但 `_has_finer_level()` 判定时必须把它们视为「比国家更细的层级」,
# 以便交给 gazetteer 拿到真实坐标, 而不是降级成中国中心点。
_CN_SUBREGION_KEYS = frozenset({
    "taiwan", "taiwan, china", "台湾", "中国台湾",
    "hong kong", "hong kong, china", "香港", "中国香港",
    "macao", "macau", "macao, china", "macau, china", "澳门", "中国澳门",
})


def _norm_key(value) -> str:
    """查找键归一化: 转写 → NFKD 去音符 → 小写 → 仅保留字母数字。

    与构建脚本 `build_gazetteer.norm_key` / 参考实现
    `gazetteer.ts: normalizedLookupKey` 等价。
    """
    s = "".join(_TRANSLITERATION.get(ch, ch) for ch in str(value))
    s = unicodedata.normalize("NFKD", s)
    s = "".join(ch for ch in s if not unicodedata.combining(ch))
    return _NORM_RE.sub("", s.lower())


def _load_gazetteer():
    """惰性加载全球 gazetteer。加载失败返回 None (降级为仅中国城市字典)。"""
    global _GAZ, _GAZ_LOAD_FAILED
    if _GAZ is not None or _GAZ_LOAD_FAILED:
        return _GAZ
    try:
        with open(_GAZETTEER_PATH, encoding="utf-8") as fh:
            entries = json.load(fh)
    except (OSError, ValueError):
        _GAZ_LOAD_FAILED = True
        return None

    exact, normalized = {}, {}
    for e in entries:
        cand = {
            "name": e["name"], "iso2": e.get("iso2", ""),
            "lat": e["lat"], "lon": e["lon"], "kind": e.get("kind", ""),
            "countryIso2": (e.get("countryIso2") or "").upper(),
            "countryCode3": (e.get("countryCode3") or "").upper(),
            "nameZh": e.get("nameZh"),
            # 城市层: 同名消歧的次级判据 (Alexandria → 取埃及那只而非美国的)
            "popMax": int(e.get("popMax") or 0),
        }
        keys = [e.get("name"), e.get("nameZh") or "", e.get("iso2", "")]
        keys += list(e.get("aliases") or [])
        for k in keys:
            if not k:
                continue
            for bucket, key in ((exact, str(k).strip().lower()), (normalized, _norm_key(k))):
                if key:
                    bucket.setdefault(key, []).append(cand)
    _GAZ = (exact, normalized)
    return _GAZ


def _fmt_coord(lat: float, lon: float) -> str:
    """输出 NCBI 格式 `XX.XX N XXX.XX E` (纬度 2 位小数, 经度 3 位)。"""
    ns = "N" if lat >= 0 else "S"
    ew = "E" if lon >= 0 else "W"
    return f"{abs(lat):.2f} {ns} {abs(lon):.3f} {ew}"


def _lookup_gazetteer(term: str, prefer_iso2: str = ""):
    """在 gazetteer 里查一个地名。

    返回候选 dict 或 None。`prefer_iso2` 用于消歧: 同名地名优先取该国的
    (如 `Maryland` 同时属美国 US-MD 与利比里亚 LR-MY)。
    优先级: admin1 (省/州) > country (国家) > city (城市) > region;
    同级内先按 `prefer_iso2` 命中, 再按**人口**从大到小 (同名小城镇很多,
    取人口最大者最可能是被引用那个)。

    ⚠ 为什么 city 排在 admin1/country **之后** (而不是按"越细越优先"):
      · `Brazil`(国名) 与 `Brazil, Indiana`(同名市) → 必须取**国家**,
        否则 C11a 判据 (country 层级) 会失守;
      · `Brazil:Sao Paulo` 与圣保罗市同名 → 必须取 **admin1 BR-SP**,
        否则 C11a 判据 (admin1 层级) 会失守。
      把 city 放在末尾, 城市层就只填补 admin1/country 都命中不了的空白,
      对既有解析结果是**纯增量**(实测: 全库 15 个不同地名解析结果一条未变)。
    """
    gaz = _load_gazetteer()
    if not gaz:
        return None
    exact, normalized = gaz
    key_exact = str(term).strip().lower()
    cands = exact.get(key_exact) or normalized.get(_norm_key(term)) or []
    if not cands:
        return None
    # 层级精确度: admin1 > country > city > region
    rank = {"admin1": 0, "country": 1, "city": 2, "region": 3}
    ordered = sorted(cands, key=lambda c: (rank.get(c["kind"], 9), -c["popMax"]))
    if prefer_iso2:
        p = prefer_iso2.upper()
        for c in ordered:
            if c["countryIso2"] == p or c["iso2"].upper() == p:
                return c
    return ordered[0]


# 地名里要剔除的噪声后缀（数据里常见 `_AI` 之类的来源标记）。
#
# ⚠ 必须要求**分隔符存在**（`_AI` / `-AI` / `.AI`）或整体就是个噪声标记，
#   绝不能只按「结尾是 ai」判断 —— 否则会误伤以 ai 结尾的真实地名：
#   `Shanghai`→`Shangh`、`Qinghai`→`Qingh`、`Chennai`→`Chenn`、`Zhanjiang` 等。
#   这是 2026-09-16 第二轮实测发现并修掉的 bug。
_LOC_NOISE_RE = re.compile(
    r"[_\-](?:ai|v[1-9]|new|final|rev\d*|copy)$", re.IGNORECASE
)


def _split_loc(loc: str):
    """把地名拆成候选词条与消歧提示。

    支持的分隔与写法:
      · `China, Ningxia, Yinchuan_AI`  逗号分段
      · `China:Ningxia`                冒号分段 (NCBI 国:区)
      · `United States:Maryland`       国:省
      · `China > Ningxia > Yinchuan`   层级符

    返回 `(parts, country_hint)`:
      parts        —— 各段词条, **保持原顺序**（不再反转）
      country_hint —— 国家前缀的 ISO2, 用于同名地名消歧

    ⚠ 不在此处决定「由粗到细」还是「由细到粗」—— 两种写法都真实存在:
      · `China:Ningxia`(国:省)    → 细的在**后**
      · `Macao, China`(地名,国)   → 细的在**前**
    故本函数只负责切分, 顺序策略交给 `_coord_from_gazetteer` 两级都试。
    """
    s = str(loc).strip()
    for sep in ("|", ">", ";", ":"):
        s = s.replace(sep, ",")
    parts = [p.strip() for p in s.split(",") if p.strip()]
    parts = [_LOC_NOISE_RE.sub("", p).strip() for p in parts]
    parts = [p for p in parts if p]
    if not parts:
        return [], ""

    # 消歧提示: 任一段若是国家名/ISO 码即采用
    country_hint = ""
    for p in parts:
        key = p.strip().lower()
        if key in _COUNTRY_NAME_TO_ISO2:
            country_hint = _COUNTRY_NAME_TO_ISO2[key]
            break
        if len(p.strip()) == 2 and p.strip().isalpha():
            country_hint = p.strip().upper()
            break
    return parts, country_hint


def _coord_from_gazetteer(loc: str):
    """用全球 gazetteer 解析地名, 返回 (coord, source) 或 None。

    source 形如 `admin1:CN-NX` / `country:BR` / `city:US` /
    `admin1:US-MD (from Maryland)`。
    两轮尝试（细粒度优先, 兼顾两种写法）:
      1. **由细到粗**（反转顺序）—— 处理 `China:Ningxia`、`United States:Maryland`
      2. **由粗到细**（原顺序）  —— 处理 `Macao, China`、`Taiwan, China`
    每轮内跳过占位符与国家名段；命中即返回。这样 `Macao, China` 不会被
    首段 `China` 抢先解析成国家质心。
    """
    parts, country_hint = _split_loc(loc)
    if not parts:
        return None

    def _try(term: str):
        if _is_placeholder_text(term):
            return None
        hint = country_hint
        if not hint:
            hint = _COUNTRY_NAME_TO_ISO2.get(term.strip().lower(), "")
        cand = _lookup_gazetteer(term, hint)
        if not cand:
            return None
        # 城市没有 ISO 3166-2 (`iso2` 留空), 退化用 countryIso2 / countryCode3,
        # 使来源标记仍是可读的 `city:US` 而不是 `city:`。
        tag = f"{cand['kind']}:{cand['iso2'] or cand['countryIso2'] or cand.get('countryCode3', '')}"
        if hint and cand["countryIso2"] == hint and cand["iso2"].upper() != hint:
            tag += f" (from {term})"
        return _fmt_coord(cand["lat"], cand["lon"]), tag

    def _is_pure_country_term(t: str) -> bool:
        """是"纯国家名"吗 —— 即**可以整段忽略**的层级。

        ⚠ 必须排除 `_CN_SUBREGION_KEYS` (台湾/香港/澳门): 它们为消歧被登记进了
          `_COUNTRY_NAME_TO_ISO2`, 但在中国主权下是**一级行政区**, 属于
          「比国家更细的层级」。若一并跳过, `Macao, China` 就没有任何段落可解析
          → 回落成中国中心点 (实测 bug: 22.13 N 113.56 E 被换成 35.86 N 104.19 E)。

        实现统一走模块级 `_is_country_only_term()` (同一个判据同时服务
        `parse_location` 的 `country_only` 判定, 避免两处口径漂移)。
        """
        return _is_country_only_term(t)

    # 轮 1: 由细到粗。跳过纯国家名段（除非整串只有国家名）
    for term in reversed(parts):
        if _is_pure_country_term(term) and len(parts) > 1:
            continue
        got = _try(term)
        if got:
            return got
    # 轮 2: 由粗到细（处理 `Macao, China` 这类「地名 + 国家」写法）
    #
    # ⚠ 2026-09-16 第三轮: 轮 2 **也必须**跳过纯国家名段 ——
    #   否则 `China, Unknown, Unknown_AI` 会在这里被 "China" 命中, 拿到
    #   Natural Earth 的国家标注点 (32.50 N 106.337 E), 等于编造地理信息
    #   (实测泄漏: 第 2 级已拦住手校中心点, 却被这里绕过)。
    #   纯国家名的回落**统一**由 `coord_from_location` 的第 2/4 级负责
    #   (`China` → 手校中心点 35.86 N 104.19 E, 与原实现逐字一致)。
    for term in parts:
        if _is_pure_country_term(term) and len(parts) > 1:
            continue
        got = _try(term)
        if got:
            return got
    return None


def _is_placeholder_text(s: str) -> bool:
    """判断是否占位符文本 (不含坐标/日期语义, 仅占位标记)。"""
    t = str(s).strip().lower()
    if not t:
        return True
    if t in ("country", "region", "unknown", "na", "nan", "not_provided",
             "none", "unknown:unknown", "country:region"):
        return True
    # 全部是 X/Y 占位或尖括号
    if re.fullmatch(r"[x\s.]+", t) or re.fullmatch(r"[<\-\s]+", t):
        return True
    return False


def coord_from_location(loc) -> Tuple[str, str]:
    """从地名推算经纬度字符串 (NCBI 格式 `XX.XX N XXX.XX E`)。

    返回 `(coord, source)`:
      coord  —— 坐标字符串; 推不出时为 `''`
      source —— 来源标记, 取值:
                  `'city'`              命中 30 个中国城市手校字典 (最精确)
                  `'country_centroid'`  仅命中 "China" 时的中国中心点
                  `'admin1:<ISO>'`      全球一级行政区 (省/州)
                  `'country:<ISO>'`     全球国家
                  `''`                  推不出

    查找优先级:
      1. `_CITY_COORDS` (30 个中国城市, 手校) —— **最高优先**,
         保证与 submission pipeline `geo_to_latlon` 的对拍结果不变
      2. 全球 admin1 (省/州) —— 由细到粗先试, 兼顾 `China:Ningxia` 与
         `Macao, China` 两种写法
      3. 全球 admin0 (国家)
      4. `''`

    行为要点:
      · `'China'` 单独出现时仍返回原中国中心点 `country_centroid`,
        与原 `geo_to_latlon` **完全一致** (避免对拍出现差异)
      · 全球查找是对中国字典的**纯增量**: 命中中国城市时行为一字不变
      · 中国台湾/香港/澳门 由构建层规范化为中国的一部分, 不会输出
        与国家主权立场不符的表述 (`Taiwan` → `admin1:CN-TW`)
    """
    if loc is None:
        return "", ""
    s = str(loc).strip()
    if not s or _is_placeholder_text(s):
        return "", ""

    low = s.lower()
    # ── 第 1 级: 中国城市手校字典 (最高优先, 与 geo_to_latlon 逐字一致) ──
    for city, latlon in _CITY_COORDS.items():
        if city in low:
            return latlon, "city"

    # ── 第 2 级: **纯国家名** → 原中国中心点 ──
    #   `China` 单独出现时, 原 `geo_to_latlon` 就是返回手校的 `35.86 N 104.19 E`。
    #   **保持逐字一致**, 不引入 Natural Earth 的 China 标注点
    #   (32.50 N 106.337 E) —— 后者虽也合理, 但会让本可零变更的场景产生漂移。
    #
    #   ⚠ 2026-09-16 第三轮修 (老师: 「Unknown:Unknown 这些在清理过程应该去掉」):
    #   若地名里**显式带了占位层级** (`China, Unknown, Unknown_AI` / `China:Unknown`),
    #   说明提交者根本不知道采样省份 —— 此时赋"中国中心点"是**编造地理信息**
    #   (误差可达上千公里), 会污染系统地理分析。改为返回 `''` (不可解析),
    #   由 `clean` 的 `--drop-no-location` 连同其序列一并剔除。
    #   判据区分: `China`(只知道国家, 合法回落) vs `China, Unknown, Unknown`(无地理信息)。
    if ("china" in low and not _has_finer_level(s)
            and not _has_placeholder_level(s)):
        return _COUNTRY_CENTROID_CN, "country_centroid"

    # ── 第 3 级: 全球 gazetteer (admin1 省州 → admin0 国家) ──
    # ⚠ 放在「中国中心点兜底」**之后**但**之前于**最终兜底:
    #   · `China:Ningxia` / `China:Qinghai` → 交给 gazetteer 拿真实省坐标
    #   · `Macao, China` / `Taiwan, China`  → 首段是地名, 交给 gazetteer
    hit = _coord_from_gazetteer(s)
    if hit:
        return hit

    # ── 第 4 级: 兜底 —— 含 "China" 但更细层级推不出 → 中国中心点 ──
    #   ⚠ 与第 2 级同款限制: 显式占位层级存在时不编造坐标 (见第 2 级注释)。
    if "china" in low and not _has_placeholder_level(s):
        return _COUNTRY_CENTROID_CN, "country_centroid"
    return "", ""


def _has_placeholder_level(s: str) -> bool:
    """地名里是否**显式**出现占位层级 (`Unknown` / `Country` / `Region` / `NA` …)。

    用于把「只知道国家」(→ `China`) 与「知道国家但更细层级是占位符」
    (→ `China, Unknown, Unknown_AI`) 区分开:

      · 前者是**合法**的国家级记录 —— 允许回落到国家中心点 (与原实现逐字一致);
      · 后者**没有任何地理信息** —— 赋国家中心点等于编造数据 (误差上千公里),
        必须返回不可解析, 交由 `clean` 剔除。

    这是 2026-09-16 第三轮针对「Unknown:Unknown 这类记录应在清理过程去掉」的判据。
    """
    parts, _ = _split_loc(s)
    return any(_is_placeholder_text(p) for p in parts)


def _has_finer_level(s: str) -> bool:
    """判断 `China, X` 里是否有比「国家」更细的层级可供 gazetteer 解析。

    只有「国家名 + 占位段」才算无更细层级 (→ 用原中国中心点)。

    ⚠ 中国台湾/香港/澳门**必须**算作「更细层级」—— 它们在中国主权下是
      一级行政区, 而非"国家名"。否则 `Macao, China` / `Taiwan, China`
      这类写法会被判为"无更细层级", 降级成中国中心点 (实测 bug)。
    """
    for p in _split_loc(s)[0]:
        if _is_placeholder_text(p):
            continue
        key = p.strip().lower()
        if key in _CN_SUBREGION_KEYS:
            return True          # 港澳台: 是更细层级, 交给 gazetteer
        if key in _COUNTRY_NAME_TO_ISO2:
            continue             # 纯国家名: 不算更细层级
        return True              # 其他地名 (省/市): 是更细层级
    return False


def _pick_col(cols: List[str], hints, taken: set) -> Optional[str]:
    """列名识别: 精确 > 前缀 > 包含; 排除已占用列 (修历史"三 pick 撞同一列" bug)。"""
    low = [(c, c.lower().strip()) for c in cols]
    for h in hints:
        for c, lc in low:
            if lc == h and c not in taken:
                return c
    for h in hints:
        for c, lc in low:
            if lc.startswith(h) and c not in taken:
                return c
    for h in hints:
        for c, lc in low:
            if h in lc and c not in taken:
                return c
    return None


def _pick_cols_all(cols: List[str], hints) -> List[str]:
    """识别**所有**匹配的列 (坐标列可能有多份: src-Lat_Lon / bs-geo_loc_name 等)。"""
    out = []
    low = [(c, c.lower().strip()) for c in cols]
    for c, lc in low:
        if any(h == lc or lc.startswith(h) or h in lc for h in hints):
            out.append(c)
    return out


def govern_table(path, *, date_col: Optional[str] = None,
                 loc_col: Optional[str] = None, name_col: Optional[str] = None,
                 coord_cols: Optional[List[str]] = None,
                 mode: str = "mid", drop_unknown_geo: bool = True,
                 blank_placeholders: bool = True,
                 derive_coords: bool = True,
                 drop_no_date: bool = False, drop_no_location: bool = False,
                 drop_placeholder_rows: bool = False,
                 encoding: str = "utf-8-sig") -> Dict:
    """读一张 metadata 表 → 治理结果。

    坐标处置 (2026-09-16 老师拍板「**先补算, 再清空**」):
      · derive_coords=True (默认) —— 坐标列若为占位符/空, 先用**地名列**的值
        调 `coord_from_location()` 推算坐标 (与 submission pipeline 的
        `geo_to_latlon` 同逻辑); 推算成功则**填入**并在报告标 `coord_source`
        (`city` / `country_centroid`); 推算不出才清空。
      · derive_coords=False —— 回到旧行为: 占位符直接清空, 不推算。

      ⚠ 补算值是**推算约值**, 不是原始观测。`coord_source` 必须一并消费,
        不可与真实坐标混同。`country_centroid` 误差可达上千公里。

    数据质量处置 (2026-09-16 老师拍板: 质量问题的处置为「仅清空该字段」):
      · blank_placeholders=True (默认) —— 把占位符字段**清空**, **不丢行**;
        每个被清空的字段在报告里留痕 (issue='placeholder')。
      · drop_placeholder_rows=True —— 若某行含占位符字段, 则**整行剔除**。
      · drop_no_date / drop_no_location —— 按「该行是否缺日期/缺地理」剔除。

    ⚠ 「坐标列为空 → 判该行无地理」: 坐标列全为占位符时等价于没有坐标,
      地理信息以地名列 (loc_col) 为准。故 coord_cols 的占位符**不单独触发剔除**,
      只清空并计数, 便于报告。

    返回 dict:
      rows        : 标准化后的行 (list[dict])
      report      : 逐行体检记录 (list[dict])
      stats       : 统计
      date_col / loc_col / name_col / coord_cols : 实际使用的列名
      columns     : 原表列名
    """
    p = Path(path)
    with open(p, encoding=encoding, errors="replace", newline="") as f:
        # 分隔符嗅探 (2026-09-16 修): 批处理入口的规范源是
        # `Global_Unified_Metadata_Core14.tsv` —— **TSV**; 而此前 csv.DictReader
        # 只按逗号切 → 整行被当成**一个**列名
        # (`识别列: date='Run\tCollectionDate\tLocation\tScientificName\tTissue'`),
        # 于是列识别全灭、日期 100% unparsable, 治理形同空转 (batch 入口实测)。
        head = f.readline()
        f.seek(0)
        delim = "\t" if head.count("\t") > head.count(",") else ","
        reader = csv.DictReader(f, delimiter=delim)
        cols = [c for c in (reader.fieldnames or []) if c]
        rows_in = list(reader)

    taken: set = set()
    d_c = date_col or _pick_col(cols, DATE_HINTS, taken)
    if d_c:
        taken.add(d_c)
    l_c = loc_col or _pick_col(cols, LOC_HINTS, taken)
    if l_c:
        taken.add(l_c)
    n_c = name_col or _pick_col(cols, NAME_HINTS, taken)

    # 坐标列: 显式给则用; 否则识别所有匹配列 (可能多份), 且排除已占用的日期/地名列
    if coord_cols:
        co_cs = [c for c in coord_cols if c in cols]
    else:
        co_cs = [c for c in _pick_cols_all(cols, COORD_HINTS)
                 if c not in (d_c, l_c, n_c)]

    out_rows, report = [], []
    stats = {
        "n_total": len(rows_in),
        "date": {"day": 0, "month": 0, "year": 0, "none": 0},
        "date_issues": {},
        "geo_issues": {},
        "geo_dropped_levels": 0,
        "geo_drop_unknown": 0,
        # 地理"判废"计数 (2026-09-16 第三轮):
        #   geo_country_only —— 删 Unknown 后只剩纯国家名 → 无分辨率
        #   geo_unusable     —— 含 country_only 在内的全部"无可用地理"行数
        "geo_country_only": 0,
        "geo_unusable": 0,
        "already_decimal_in_col": 0,
        "kept": 0,
        "dropped_no_date": 0,
        "dropped_no_location": 0,
        "dropped_placeholder": 0,
        # 字段级清空计数: {列名: 清空条数}
        "blanked_fields": {},
        "coord_placeholder_fields": {},
        # 坐标补算 (2026-09-16 老师拍板「先补算再清空」):
        #   coords_derived            {列名: 成功补算条数}
        #   coords_derived_by_source  {"city": n, "country_centroid": m}
        #   coords_unrecoverable      {列名: 补不出而清空的条数}
        "coords_derived": {},
        "coords_derived_by_source": {},
        "coords_unrecoverable": {},
        "n_rows_with_placeholder": 0,
        "date_col": d_c, "loc_col": l_c, "name_col": n_c,
        "coord_cols": list(co_cs),
    }

    def _bump(d, k):
        d[k] = d.get(k, 0) + 1

    for i, r in enumerate(rows_in, start=2):  # 行号从 2 起 (含表头)
        rec = {"row": i, "name": (r.get(n_c) or "") if n_c else ""}
        no_date = no_loc = False
        blanked = []          # 本行被清空的字段
        coord_derived = {}    # {坐标列: 本行补算出的坐标}
        coord_src = {}        # {坐标列: 'city' | 'country_centroid'}

        if d_c:
            di = parse_date(r.get(d_c), mode=mode)
            rec.update({
                "date_raw": di.raw, "date_iso": di.iso or "",
                "date_decimal": "" if di.decimal is None else f"{di.decimal:.10f}",
                "date_precision": di.precision, "date_confidence": di.confidence,
                "date_format": di.source_format, "date_issue": di.issue or "",
                "date_consumed_as": di.consumed_as,
            })
            stats["date"][di.precision] = stats["date"].get(di.precision, 0) + 1
            if di.issue:
                _bump(stats["date_issues"], di.issue)
            if di.consumed_as == "already_decimal":
                stats["already_decimal_in_col"] += 1
            no_date = di.decimal is None
            if di.issue == "placeholder" and blank_placeholders:
                blanked.append(d_c)

        if l_c:
            gi = parse_location(r.get(l_c), drop_unknown=drop_unknown_geo)
            rec.update(gi.to_row())
            if gi.issue:
                _bump(stats["geo_issues"], gi.issue)
            if gi.dropped_levels:
                stats["geo_dropped_levels"] += len(gi.dropped_levels)
                stats["geo_drop_unknown"] += 1
            # `country_only`: 删掉占位层级后只剩纯国家名 → 无可用地理分辨率。
            # 单列计数, 便于报告里一眼看到"被判废了多少条"(老师:
            # 「Unknown:Unknown 这些在清理过程应该去掉」)。
            if gi.issue == "country_only":
                stats["geo_country_only"] += 1
            if not gi.usable:
                stats["geo_unusable"] += 1
            # ⚠ 用 `usable` 而非 `not components`: 后者会把 `China, Unknown,
            #   Unknown_AI` 当"有地理"放行, 于是它带着空坐标活到系统地理分析里。
            no_loc = not gi.usable
            # `unknown_level`(仍有更细层级) 不在此列 —— 它有可用的省/市, 只是部分层级缺失。
            if gi.issue in ("placeholder", "country_only") and blank_placeholders:
                blanked.append(l_c)

        # 坐标列: 「**先补算, 再清空**」(2026-09-16 老师拍板)。
        # 占位符/空 → 先用**地名列** (loc_col) 的值推算坐标; 推算成功则填入,
        # 推不出才清空。补算值是约值, 来源记入 rec['coord_source'] 供下游分辨。
        # ⚠ 老师拍板「坐标列为空 → 判该行无地理」—— 坐标占位符**不得**触发
        #   `drop_placeholder_rows`, 否则 src-Lat_Lon 450/450 占位符会把整表剔光。
        #   坐标占位符只清空 + 计数; 是否"无地理"由地名列 (loc_col) 决定。
        for cc in co_cs:
            v = r.get(cc)
            if is_placeholder(v) or not str(v or "").strip():
                dval, dsrc = ("", "")
                if derive_coords and l_c:
                    dval, dsrc = coord_from_location(r.get(l_c))
                if dval:
                    coord_derived[cc] = dval
                    coord_src[cc] = dsrc
                    _bump(stats["coords_derived"], cc)
                    _bump(stats["coords_derived_by_source"], dsrc)
                else:
                    if blank_placeholders:
                        blanked.append(cc)
                    _bump(stats["coord_placeholder_fields"], cc)
                    _bump(stats["coords_unrecoverable"], cc)

        blanked = sorted(set(blanked))
        # 触发整行剔除的字段: 排除坐标列 (坐标占位符不等于"该行不可用")
        blanked_for_drop = [b for b in blanked if b not in co_cs]
        if blanked:
            stats["n_rows_with_placeholder"] += 1
            rec["blanked_fields"] = "|".join(blanked)
            for b in blanked:
                _bump(stats["blanked_fields"], b)
        else:
            rec["blanked_fields"] = ""
        # 坐标补算来源: `列名:city` / `列名:country_centroid`; 空 = 本行无需补算
        rec["coord_source"] = "|".join(f"{k}:{v}" for k, v in coord_src.items())

        # ── 剔除判定 ──
        drop_why = None
        if drop_placeholder_rows and blanked_for_drop:
            drop_why = "placeholder"
        elif drop_no_date and no_date:
            drop_why = "no_date"
        elif drop_no_location and no_loc:
            drop_why = "no_location"

        if drop_why == "placeholder":
            stats["dropped_placeholder"] += 1
            rec["action"] = "dropped:placeholder"
            report.append(rec)
            continue
        if drop_why == "no_date":
            stats["dropped_no_date"] += 1
            rec["action"] = "dropped:no_date"
            report.append(rec)
            continue
        if drop_why == "no_location":
            stats["dropped_no_location"] += 1
            rec["action"] = "dropped:no_location"
            report.append(rec)
            continue

        rec["action"] = "kept"
        stats["kept"] += 1
        report.append(rec)

        new = dict(r)
        # 日期: 能解析则写归一 ISO; 占位符则清空
        if d_c:
            if rec.get("date_iso"):
                new[d_c] = rec["date_iso"]
            elif rec.get("date_issue") == "placeholder" and blank_placeholders:
                new[d_c] = ""
        if l_c:
            if rec.get("geo_normalized"):
                new[l_c] = rec["geo_normalized"]
            elif rec.get("geo_issue") == "placeholder" and blank_placeholders:
                new[l_c] = ""
        # 坐标: 补算值优先写入; 其余占位符 (推不出) 才清空
        for cc, dv in coord_derived.items():
            if cc in new:
                new[cc] = dv
        for cc in blanked:
            if cc in new and cc not in coord_derived:
                new[cc] = ""
        new["_date_decimal"] = rec.get("date_decimal", "")
        new["_date_precision"] = rec.get("date_precision", "")
        new["_date_issue"] = rec.get("date_issue", "")
        new["_geo_unresolved"] = rec.get("geo_unresolved", "")
        out_rows.append(new)

    stats.setdefault("coord_placeholder_fields", {})
    return {"rows": out_rows, "report": report, "stats": stats,
            "date_col": d_c, "loc_col": l_c, "name_col": n_c,
            "coord_cols": co_cs, "columns": cols}


def write_outputs(result: Dict, out_csv, report_tsv, *, quiet: bool = False) -> Dict:
    """写标准化 metadata (out_csv) 与体检报告 (report_tsv)。"""
    out_csv, report_tsv = Path(out_csv), Path(report_tsv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)

    fieldnames = list(result["columns"])
    for extra in ("_date_decimal", "_date_precision", "_date_issue", "_geo_unresolved"):
        if result["rows"] and extra in result["rows"][0]:
            fieldnames.append(extra)
    with open(out_csv, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        w.writeheader()
        w.writerows(result["rows"])

    rep_cols = ["row", "name", "action", "blanked_fields", "coord_source",
                "date_raw", "date_iso", "date_decimal",
                "date_precision", "date_format", "date_confidence", "date_issue",
                "date_consumed_as", "geo_raw", "geo_normalized", "geo_components",
                "geo_dropped", "geo_unresolved", "geo_issue"]
    with open(report_tsv, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=rep_cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for rec in result["report"]:
            w.writerow({k: rec.get(k, "") for k in rep_cols})

    return {"std_csv": str(out_csv), "report_tsv": str(report_tsv)}


# ══════════════════════════════════════════════════════════════════
# CLI
# ══════════════════════════════════════════════════════════════════

def _main(argv=None) -> int:
    import argparse
    ap = argparse.ArgumentParser(
        description="metadata 入口治理: 时间/地理 归一 + 标注 + 体检")
    ap.add_argument("--metadata", required=True, help="输入 metadata CSV")
    ap.add_argument("--out", required=True, help="输出标准化 CSV (另存新文件)")
    ap.add_argument("--report", required=True, help="输出体检报告 TSV")
    ap.add_argument("--date-col", default=None)
    ap.add_argument("--loc-col", default=None)
    ap.add_argument("--name-col", default=None)
    ap.add_argument("--coord-cols", default=None,
                    help="坐标列名, 逗号分隔 (默认自动识别 lat_lon 类列)")
    ap.add_argument("--mode", choices=("mid", "start"), default="mid",
                    help="mid=月中/年中 (decimal_year); start=月初/年初 (virphy_bridge)")
    ap.add_argument("--keep-unknown-geo", action="store_true",
                    help="保留地理中的 Unknown 层级 (默认删除, 老师拍板)")
    ap.add_argument("--no-derive-coords", action="store_true",
                    help="**不补算**坐标 (默认: 坐标占位符先用地名推算, 推不出才清空)")
    ap.add_argument("--no-blank-placeholders", action="store_true",
                    help="**不清空**占位符字段 (默认清空: 老师拍板「仅清空该字段」)")
    ap.add_argument("--drop-placeholder-rows", action="store_true",
                    help="含占位符字段的行**整行剔除** (默认不剔, 只清空)")
    ap.add_argument("--drop-no-date", action="store_true",
                    help="剔除日期不可解析的行")
    ap.add_argument("--drop-no-location", action="store_true",
                    help="剔除地理为空的行")
    args = ap.parse_args(argv)

    co_cs = [s.strip() for s in args.coord_cols.split(",")] if args.coord_cols else None
    res = govern_table(args.metadata, date_col=args.date_col, loc_col=args.loc_col,
                       name_col=args.name_col, coord_cols=co_cs, mode=args.mode,
                       drop_unknown_geo=not args.keep_unknown_geo,
                       derive_coords=not args.no_derive_coords,
                       blank_placeholders=not args.no_blank_placeholders,
                       drop_placeholder_rows=args.drop_placeholder_rows,
                       drop_no_date=args.drop_no_date,
                       drop_no_location=args.drop_no_location)
    paths = write_outputs(res, args.out, args.report)
    s = res["stats"]
    print(f"输入: {args.metadata}  ({s['n_total']} 行)")
    print(f"  日期列={s['date_col']!r}  地理列={s['loc_col']!r}  名称列={s['name_col']!r}  mode={args.mode}")
    print(f"  坐标列={s['coord_cols']}")
    print(f"  日期精度: " + "  ".join(f"{k}={v}" for k, v in s["date"].items()))
    if s["date_issues"]:
        print(f"  日期问题: {s['date_issues']}")
    if s["geo_issues"]:
        print(f"  地理问题: {s['geo_issues']}")
    print(f"  已删除 Unknown 层级: {s['geo_drop_unknown']} 行 / {s['geo_dropped_levels']} 级")
    if s["coords_derived"]:
        print(f"  坐标补算成功: {s['coords_derived']}  来源={s['coords_derived_by_source']}")
    if s["coords_unrecoverable"]:
        print(f"  坐标补算不出 (已清空): {s['coords_unrecoverable']}")
    if s["coord_placeholder_fields"]:
        print(f"  坐标占位符: {s['coord_placeholder_fields']}")
    if s["blanked_fields"]:
        print(f"  字段清空统计: {s['blanked_fields']}  (含占位符的行 {s['n_rows_with_placeholder']} 条)")
    print(f"  列内已是小数年: {s['already_decimal_in_col']} 条")
    print(f"  保留 {s['kept']} 行;  剔占位符 {s['dropped_placeholder']};  "
          f"无日期 {s['dropped_no_date']};  无地理 {s['dropped_no_location']}")
    print(f"  → {paths['std_csv']}")
    print(f"  → {paths['report_tsv']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(_main())
