#!/usr/bin/env python3
"""seq_match.py — 序列名匹配 / 参考序列解析 (全管线统一口径)

背景 (2026-09-15 复审):
  · 多处用 `acc in name` 做子串判断 → `CRR1` 会命中 `CRR10`, 删错序列 / 选错参考;
  · 多处直接把 `seqs[0]` 当参考序列且无任何校验 → MAFFT 输出顺序不保证参考在首位,
    坐标基准会静默漂移, 下游 CDS 坐标整体错位且没有任何提示。

本模块提供统一的「精确 → token → 子串」三级匹配与参考解析, 供
popgen_analysis / mask_recombination / dnasp_bridge / virphy_bridge 复用,
避免每个模块各写一份、口径不一致。
"""
from __future__ import annotations

import re
import sys
from typing import List, Optional, Sequence, Tuple

_TOKEN_SPLIT = re.compile(r'[_\-.|:/\s]+')


def _tokens(s: str) -> set:
    return {t for t in _TOKEN_SPLIT.split(str(s)) if t}


def match_by_name(pairs: Sequence[Tuple[str, str]], key: Optional[str],
                  allow_substring: bool = True):
    """在 [(name, seq), ...] 中按 精确 > token > 子串 三级匹配 key。

    返回 (name, seq) 或 None。allow_substring=False 时禁用最宽松的子串兜底
    (用于"必须精确对应"的场景, 如重组序列删除)。
    """
    if key is None:
        return None
    key = str(key)
    for nm, sq in pairs:
        if nm == key:
            return nm, sq
    kt = _tokens(key)
    if kt:
        for nm, sq in pairs:
            nt = _tokens(nm)
            # 双向 token 包含: 兼容 "key 比 name 长" (X1_extra → X1)
            # 与 "name 比 key 长" (CRR1_x → CRR1) 两种常见写法
            if nt and (kt <= nt or nt <= kt):
                return nm, sq
    if allow_substring:
        for nm, sq in pairs:
            if key in nm:
                return nm, sq
    return None


def is_same_accession(name: str, acc: str) -> bool:
    """判断序列名是否"就是"某个 accession (精确或作为独立 token 出现)。

    与 `acc in name` 的区别: `CRR1` 不会命中 `CRR10` (后者不是独立 token)。
    """
    if name == acc:
        return True
    return acc in _tokens(name)


def resolve_reference(pairs: Sequence[Tuple[str, str]],
                      reference: Optional[str] = None,
                      log=None, who: str = "reference"):
    """解析参考序列 -> (name, seq)。

    reference 显式给出: 三级匹配; 找不到抛 ValueError (拒绝静默换基准)。
    未给出: 退回第一条并**明确告警** (MAFFT 输出顺序不保证参考在首位)。
    """
    pairs = list(pairs)
    if not pairs:
        raise ValueError(f"{who}: 序列列表为空, 无法确定参考")
    if reference:
        hit = match_by_name(pairs, reference)
        if hit is None:
            names = [nm for nm, _ in pairs]
            raise ValueError(
                f"{who}: 参考序列 '{reference}' 不在输入中 "
                f"(共 {len(names)} 条: {names[:8]}{'...' if len(names) > 8 else ''})")
        return hit

    msg = (f"{who}: 未指定参考序列 → 默认取第一条 '{pairs[0][0]}' 作为坐标基准; "
           f"MAFFT 输出顺序不保证参考在首位, 坐标可能整体错位, 建议显式指定")
    if log is not None:
        for m in ("warning", "emit", "error"):
            f = getattr(log, m, None)
            if callable(f):
                f(msg)
                break
        else:
            print(f"[警告] {msg}", file=sys.stderr)
    else:
        print(f"[警告] {msg}", file=sys.stderr)
    return pairs[0]
