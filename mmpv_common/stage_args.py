"""跨管线共享的 --stage 多阶段解析 (MMPV 全部管线入口统一接入)。

支持空格与逗号两种分隔, 可混用:
  --stage clean deplete report          # 空格分隔 (argparse nargs='+')
  --stage clean,deplete,report          # 逗号分隔 (单 token)
  --stage clean,deplete report          # 混合
行为约定:
  - 去重保序 (执行顺序仍由各管线自身的 STAGE_ORDER/依赖逻辑决定);
  - 未知阶段名 SystemExit 并列出合法值;
  - 'all' 不在本模块展开 —— 展开语义管线各异 (如 data_preprocessing 的
    all 不含 bbnorm, phylo 的 all 受 config 开关影响), 由各管线原逻辑处理。
"""
from __future__ import annotations

import sys

__all__ = ["split_stage_tokens", "validate_stage_names", "parse_multi_stage"]


def split_stage_tokens(raw):
    """把 --stage 的原始值 (str / list / tuple, nargs='+' 或单值) 拆成阶段名列表。

    逗号与空白都视为分隔符:
      "clean,deplete" / ["clean", "deplete"] / "clean, deplete" / ["clean deplete"]
    """
    tokens = []
    for item in raw if isinstance(raw, (list, tuple)) else [raw]:
        for chunk in str(item).split(','):
            tokens.extend(chunk.split())
    return tokens


def validate_stage_names(tokens, valid, what="--stage"):
    """校验阶段名并去重保序; 未知值 SystemExit 并列出合法值。"""
    seen, ordered = set(), []
    for t in tokens:
        if t not in seen:
            seen.add(t)
            ordered.append(t)
    unknown = [t for t in ordered if t not in valid]
    if unknown:
        sys.exit("ERROR: unknown %s value(s): %s. Valid: %s"
                 % (what, ", ".join(unknown), "/".join(valid)))
    return ordered


def parse_multi_stage(raw, valid, what="--stage"):
    """split_stage_tokens + validate_stage_names 一步到位, 返回干净的阶段列表。"""
    return validate_stage_names(split_stage_tokens(raw), valid, what=what)
