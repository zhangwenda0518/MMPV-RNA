#!/usr/bin/env python3
"""decimal_year.py — 日期字符串 → 小数年 (全管线唯一定义)

2026-09-15 (审查 P2-26) 背景
----------------------------
同一份 metadata 在管线里被三个模块各自换算成小数年, 而三个实现的约定都不同:

    beast1_bridge.py   yr + (mo-1)/12 + (dy-1)/365          ← 用 365 天年
    beast_bridge.py    yr + (mo-1 + (min(dy,dim)-1)/dim)/12 ← 用当月实际天数
    phylogeo_bridge.py yr + (mo-1)/12                       ← 直接丢弃"日"

后果: BEAST1 与 BEAST2 两条路径给同一条序列算出不同的采样时间 (2020-12-31 相差
约 0.0015 年 ≈ 0.55 天), 而 `YYYY-MM-DD` 在 phylogeo 路径上干脆退化成月初。
对以"年"为尺度的时间信号检验/定年而言, 这是同一份数据两个答案。

本模块给出唯一实现, 三个调用方一律复用。
"""
import calendar

__all__ = ["to_decimal_year", "to_decimal_year_safe"]


def to_decimal_year(date_str: str, sample: str = "", strict: bool = True) -> float:
    """把 YYYY / YYYY-MM / YYYY-MM-DD (也接受 '/' 分隔) 转成小数年。

    约定 (全管线统一, 与 beast_bridge 的"当月实际天数"口径一致):
        yr + (mo - 1 + (dy - 1) / days_in_month(mo)) / 12

    · 只给年月 → 取该月 15 日; 只给年 → 取 6 月 15 日
    · 日期落在月初时结果 == yr + (mo-1)/12 (与旧 phylogeo 口径兼容)
    · 闰年自动由 calendar.monthrange 处理, 无 365/366 漂移

    Parameters
    ----------
    date_str : str
    sample : str
        序列名, 仅用于错误信息定位。
    strict : bool, default True
        True  → 无法解析/越界时抛 ValueError (fail loudly);
        False → 返回 None (调用方自行决定如何处理, 绝不静默造一个年份)。
    """
    if date_str is None or not str(date_str).strip():
        if strict:
            raise ValueError(f"Empty date for sample '{sample}'")
        return None

    s = str(date_str).replace('/', '-').strip()
    parts = s.split('-')
    try:
        yr = float(parts[0])
        mo = float(parts[1]) if len(parts) > 1 else 6.0
        dy = float(parts[2]) if len(parts) > 2 else 15.0
        if not (1 <= mo <= 12):
            raise ValueError(f"month={mo} out of range [1,12]")
        if not (1 <= dy <= 31):
            raise ValueError(f"day={dy} out of range [1,31]")
        if not (1900 <= yr <= 2100):
            raise ValueError(f"year={yr} out of plausible range [1900,2100]")
        dim = calendar.monthrange(int(yr), int(mo))[1]
        dy = min(dy, float(dim))
        return yr + (mo - 1.0 + (dy - 1.0) / dim) / 12.0
    except (ValueError, IndexError, TypeError) as e:
        if strict:
            raise ValueError(
                f"Cannot parse date '{date_str}' for sample '{sample}': {e}. "
                f"Fix the metadata CSV date column before running BEAST."
            ) from e
        return None


def to_decimal_year_safe(date_str: str, sample: str = ""):
    """宽松版: 解析失败返回 None (不抛异常)。"""
    return to_decimal_year(date_str, sample=sample, strict=False)
