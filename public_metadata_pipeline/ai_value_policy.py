#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""AI 值的非破坏性采纳策略 (元数据管线共用)。

背景
----
search 阶段的 AI 洗髓、info 阶段的「9 列正交化自检防火墙」、以及 Location
批量规范化都会拿到「AI 建议值」。旧实现一律无条件覆盖表格 (df.at[idx, k] = val),
与提示词里「只做减法 / 不许改已有事实」的约束自相矛盾 —— 模型一次幻觉就能
把数据库里的真实采集时间、真实地点改写掉, 且没有任何痕迹。

本模块把「采纳什么、丢弃什么」变成一条显式、可测、可审计的规则:

  fill       旧值为空/占位词          -> 采纳
  normalize  两者语义等价, 仅规范写法 -> 采纳 (如 roots -> root, 补 `_AI` 后缀)
  enrich     新值包含旧值信息         -> 采纳 (如 2021 -> 2021-05, 补全省/市层级)
  keep       新值更弱或不可用         -> 保留旧值, 不写
  conflict   两者都非空且语义相斥     -> **保留旧值**, 冲突记入 AI_Conflict 待人工复核

原则: 永远不静默改写既有事实; AI 只允许「补空、规范、补全」, 不允许「替换」。
"""

from __future__ import annotations

import json
import os
import re
from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Sequence, Tuple

# ═══════════════════════════════════════════════════════════════════
# 常量
# ═══════════════════════════════════════════════════════════════════

AI_SUFFIX = "_AI"

#: 视为「空」的占位词 (含中英文, 大小写不敏感)
BLANK_TOKENS = {
    "", "nan", "none", "null", "na", "n/a", "nat", "<na>", "pd.na",
    "not_provided", "not provided", "not applicable", "not collected",
    "not available", "missing", "unknown", "unknow", "no evidence",
    "未知", "未提供", "未采集", "缺", "-",
}

#: 组织部位单数归一 (与 gsa_sra.info 历史映射保持一致)
TISSUE_SINGULAR = {
    "leaves": "leaf", "fruits": "fruit", "roots": "root", "stems": "stem",
    "flowers": "flower", "seeds": "seed", "myceliums": "mycelium",
}

_DATE_RE = re.compile(r"((?:19|20)\d{2})(?:[-/.年](\d{1,2}))?(?:[-/.月](\d{1,2}))?")
_LOC_SEP_RE = re.compile(r"[,，;；:：/|、]+")
_WS_RE = re.compile(r"\s+")


# ═══════════════════════════════════════════════════════════════════
# 基础判定 / 规整
# ═══════════════════════════════════════════════════════════════════

def is_blank(value: Any) -> bool:
    """空值判定: None / NaN / pd.NA / 占位词。"""
    if value is None:
        return True
    if not isinstance(value, str):
        try:
            import pandas as pd
            if pd.isna(value):
                return True
        except (TypeError, ValueError):
            pass
    s = str(value).strip()
    return s == "" or s.lower() in BLANK_TOKENS


def strip_ai(value: Any) -> str:
    """去掉末尾 `_AI` 来源标记 (大小写不敏感, 幂等)。"""
    return re.sub(re.escape(AI_SUFFIX) + r"$", "", str(value).strip(), flags=re.IGNORECASE)


def with_ai(value: Any) -> str:
    """规范为「恰好一个 `_AI` 后缀」的字符串 (幂等, 不会出现 `_AI_AI`)。"""
    core = strip_ai(value)
    return f"{core}{AI_SUFFIX}" if core else ""


def _space(value: Any) -> str:
    return _WS_RE.sub(" ", str(value).strip())


def format_location(value: Any) -> str:
    """Location 规范写法: `国家, 省/州, 市/县` + 单一 `_AI` 后缀。"""
    core = strip_ai(value)
    parts = [p.strip() for p in _LOC_SEP_RE.split(core) if p.strip()]
    return with_ai(", ".join(parts))


def format_tissue(value: Any) -> str:
    """Tissue 规范写法: 小写单数核心组织词 (剥离修饰语由提示词负责)。"""
    low = _space(value).lower().strip(".'\" ")
    return TISSUE_SINGULAR.get(low, low)


def format_date(value: Any) -> str:
    return _space(value).replace("  ", " ")


def format_age(value: Any) -> str:
    return _space(value).strip("|;- ")


def _canon_location(value: Any) -> Optional[Tuple[str, ...]]:
    core = strip_ai(value)
    parts = [p.strip().lower() for p in _LOC_SEP_RE.split(core) if p.strip()]
    return tuple(parts) if parts else None


def _canon_date(value: Any) -> Optional[Tuple[int, ...]]:
    m = _DATE_RE.search(str(value))
    if not m:
        return None
    y = int(m.group(1))
    mo = int(m.group(2)) if m.group(2) else None
    d = int(m.group(3)) if m.group(3) else None
    if mo is not None and not 1 <= mo <= 12:
        mo = None
    if d is not None and not 1 <= d <= 31:
        d = None
    if mo is None:
        return (y,)
    if d is None:
        return (y, mo)
    return (y, mo, d)


def _canon_tissue(value: Any) -> Optional[Tuple[str, ...]]:
    low = format_tissue(value)
    toks = tuple(t for t in re.split(r"[\s,;/|]+", low) if t)
    return toks or None


def _canon_age(value: Any) -> Optional[Tuple[str, ...]]:
    toks = tuple(t for t in re.split(r"[\s,;/|]+", _space(value).lower()) if t)
    return toks or None


def _subset(a: Tuple[str, ...], b: Tuple[str, ...]) -> bool:
    return bool(set(a)) and set(a) <= set(b)


# ═══════════════════════════════════════════════════════════════════
# 裁决
# ═══════════════════════════════════════════════════════════════════

@dataclass(frozen=True)
class Decision:
    """单个字段的裁决结果。value 为最终应写入的值 (未采纳时即旧值)。"""
    field: str
    old: Any
    new: Any
    value: Any
    action: str  # fill | normalize | enrich | keep | conflict

    @property
    def adopted(self) -> bool:
        return self.action in ("fill", "normalize", "enrich")


def resolve_value(field: str, old: Any, new: Any) -> Decision:
    """按字段类型裁决 AI 建议值是否可写入。

    字段支持: CollectionDate / Location / Tissue / Age_GrowthStage (其它字段
    按「空则填、等则留、异则冲突」的通用规则处理)。
    """
    if new is None or is_blank(new):
        # AI 没给出可用值 (或明确说不确定) -> 绝不删除旧值
        return Decision(field, old, new, old, "keep")

    if field == "Location":
        formatted = format_location(new)
        if not strip_ai(formatted):
            return Decision(field, old, new, old, "keep")
        if is_blank(old):
            return Decision(field, old, new, formatted, "fill")
        oc, nc = _canon_location(old), _canon_location(formatted)
        if not nc:
            return Decision(field, old, new, old, "keep")
        if not oc:
            return Decision(field, old, new, formatted, "normalize")
        if oc == nc:
            return Decision(field, old, new, formatted, "normalize")
        if _subset(oc, nc):
            # 补全层级: China -> China, Ningxia, Yinchuan
            return Decision(field, old, new, formatted, "enrich")
        if _subset(nc, oc):
            # 去掉占位/噪声层级: China, Unknown, Yinchuan -> China, Yinchuan
            return Decision(field, old, new, formatted, "normalize")
        return Decision(field, old, new, old, "conflict")

    if field == "CollectionDate":
        if is_blank(old):
            return Decision(field, old, new, format_date(new), "fill")
        oc, nc = _canon_date(old), _canon_date(new)
        if not nc:
            return Decision(field, old, new, old, "keep")
        if not oc:
            return Decision(field, old, new, format_date(new), "normalize")
        if oc == nc:
            return Decision(field, old, new, format_date(new), "normalize")
        if oc == nc[: len(oc)]:
            # 2021 -> 2021-05 (更细, 且与旧值一致)
            return Decision(field, old, new, format_date(new), "enrich")
        if nc == oc[: len(nc)]:
            # 2021-05 -> 2021 (更粗) 不采纳, 保留更精确的旧值
            return Decision(field, old, new, old, "keep")
        return Decision(field, old, new, old, "conflict")

    if field == "Tissue":
        formatted = format_tissue(new)
        if not formatted:
            return Decision(field, old, new, old, "keep")
        if is_blank(old):
            return Decision(field, old, new, formatted, "fill")
        oc, nc = _canon_tissue(old), _canon_tissue(formatted)
        if not nc:
            return Decision(field, old, new, old, "keep")
        if not oc:
            return Decision(field, old, new, formatted, "normalize")
        if oc == nc:
            return Decision(field, old, new, formatted, "normalize")
        if _subset(nc, oc):
            # 剥离修饰语: young leaves -> leaf
            return Decision(field, old, new, formatted, "normalize")
        if _subset(oc, nc):
            return Decision(field, old, new, formatted, "enrich")
        return Decision(field, old, new, old, "conflict")

    if field == "Age_GrowthStage":
        formatted = format_age(new)
        if not formatted:
            return Decision(field, old, new, old, "keep")
        if is_blank(old):
            return Decision(field, old, new, formatted, "fill")
        oc, nc = _canon_age(old), _canon_age(formatted)
        if not nc:
            return Decision(field, old, new, old, "keep")
        if not oc:
            return Decision(field, old, new, formatted, "normalize")
        if oc == nc or _subset(nc, oc) or _subset(oc, nc):
            # 合并/凝练 (3 | Young Fruit -> 3 years, young fruit stage) 只改写法
            return Decision(field, old, new, formatted, "normalize")
        return Decision(field, old, new, old, "conflict")

    # 通用字段: 空则填; 语义相同不写; 其余冲突
    if is_blank(old):
        return Decision(field, old, new, str(new).strip(), "fill")
    if _space(old).lower() == _space(new).lower():
        return Decision(field, old, new, old, "normalize")
    return Decision(field, old, new, old, "conflict")


# ═══════════════════════════════════════════════════════════════════
# 批量写入
# ═══════════════════════════════════════════════════════════════════

def _conflict_note(field: str, dec: Decision, limit: int = 60) -> str:
    old = _space(dec.old)[:limit]
    new = _space(dec.new)[:limit]
    return f"{field}: {old} ⇒ {new}"


def apply_updates(df, idx: Any, updates: Dict[str, Any], *,
                  fields: Optional[Sequence[str]] = None,
                  conflict_col: str = "AI_Conflict",
                  stage: str = "") -> Tuple[List[str], List[Dict[str, Any]]]:
    """把 AI 建议值按策略写入 df 的 idx 行。

    返回 (已采纳字段列表, 冲突记录列表)。冲突既不写值也不丢证据:
    记入 conflict_col (一行内 ` | ` 拼接) 并返回结构化记录供落盘审计。
    ``fields`` 省略时按 ``updates`` 的键逐个裁决。
    """
    adopted: List[str] = []
    conflicts: List[Dict[str, Any]] = []
    for field in (fields if fields is not None else list(updates)):
        if field not in updates:
            continue
        old = df.at[idx, field] if field in df.columns else None
        dec = resolve_value(field, old, updates[field])
        if dec.adopted:
            # 幂等: 规范化写法与旧值一致时不重复记"采纳"
            try:
                if not is_blank(old) and _space(old) == _space(dec.value):
                    continue
            except Exception:
                pass
            if field not in df.columns:
                df[field] = None
            try:
                if df[field].dtype != object:
                    df[field] = df[field].astype(object)
            except Exception:
                pass
            df.at[idx, field] = dec.value
            adopted.append(field)
        elif dec.action == "conflict":
            note = _conflict_note(field, dec)
            if conflict_col not in df.columns:
                df[conflict_col] = None
            try:
                if df[conflict_col].dtype != object:
                    df[conflict_col] = df[conflict_col].astype(object)
            except Exception:
                pass
            prev = df.at[idx, conflict_col]
            merged = note if is_blank(prev) else f"{prev} | {note}"
            df.at[idx, conflict_col] = merged
            conflicts.append({
                "row": str(df.at[idx, "Run"]) if "Run" in df.columns else "",
                "field": field,
                "old": _space(dec.old),
                "new": _space(dec.new),
                "action": dec.action,
                "stage": stage,
            })
    return adopted, conflicts


def write_conflicts(out_dir: Any, records: Sequence[Dict[str, Any]], *,
                    csv_name: str = "ai_conflicts.csv",
                    jsonl_name: str = "ai_conflicts.jsonl") -> int:
    """把冲突记录落盘: CSV 汇总 (覆盖写) + JSONL 追加; 返回写入条数。

    落盘失败不影响主流程 (审计尽力而为), 因此任何 IO 异常都吞掉。
    """
    records = [r for r in (records or []) if r]
    if not out_dir or not records:
        return 0
    try:
        os.makedirs(out_dir, exist_ok=True)
        import csv as _csv
        cols = ["run", "row", "field", "old", "new", "action", "stage"]
        with open(os.path.join(out_dir, csv_name), "w", newline="",
                  encoding="utf-8-sig") as f:
            w = _csv.DictWriter(f, fieldnames=cols, extrasaction="ignore")
            w.writeheader()
            for r in records:
                w.writerow({k: r.get(k, "") for k in cols})
        with open(os.path.join(out_dir, jsonl_name), "a", encoding="utf-8") as f:
            for r in records:
                f.write(json.dumps(r, ensure_ascii=False) + "\n")
        return len(records)
    except OSError:
        return 0
