#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
report_template.py — 报告层 Jinja2 模板统一出口
================================================

把各管线 HTML 报告的 f-string/str.format 拼接迁移到共享模板
(mmpv_common/templates/*.html.j2), 实现组件与数据模板分离:

  * base.html.j2     整页骨架 (head/CSS/页脚/时间戳)
  * macros.html.j2   规范组件宏 (kpi_card/data_table/badge/section_header)
  * 其余模板         各报告专属页 (master_report / process_panel / submission_report)

纪律 (与 report_style.py 一致): 纯本地模板、零外部资源、可离线渲染、打印友好。
安全基线: autoescape=True (元数据值里的 <>& 不再破坏标记) + StrictUndefined
(缺变量立刻抛错, 不产出"看着像成功"的半成品报告 —— fail-loud)。

用法:

    from mmpv_common.report_template import render
    html = render("master_report.html.j2", **context)

render() 自动注入两个上下文默认值 (可被同名入参覆盖):
  base_css      共享基础样式 (Markup, 不转义)
  generated_at  生成时间戳 (%Y-%m-%d %H:%M:%S)

管线专属模板 (后续迁移) 可把自己的模板目录传给 make_env(extra_dirs=[...]),
ChoiceLoader 优先查管线目录、共享目录殿后。
"""

from __future__ import annotations

import time
from functools import lru_cache
from pathlib import Path
from typing import Any, Optional, List

try:
    from jinja2 import ChoiceLoader, Environment, FileSystemLoader, StrictUndefined
    from markupsafe import Markup
except ImportError as _e:  # pragma: no cover - 环境缺 jinja2
    raise ImportError(
        "报告模板层需要 jinja2 (pixi [dependencies] 已声明; 手动环境: pip install jinja2)"
    ) from _e

try:
    from mmpv_common.report_style import PALETTE, get_base_css
except ImportError:
    import sys
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
    from mmpv_common.report_style import PALETTE, get_base_css  # type: ignore

TEMPLATES_DIR = Path(__file__).resolve().parent / "templates"


def _fmt_int(value: Any) -> str:
    """千分位整数 (None/异常值原样返回字符串, 模板里数字列统一走它)。"""
    try:
        return f"{int(value):,}"
    except (TypeError, ValueError):
        return str(value if value is not None else "")


def make_env(extra_dirs: Optional[List[Path]] = None) -> Environment:
    """构建报告模板环境。extra_dirs: 管线专属模板目录 (优先于共享目录)。"""
    loaders = []
    for d in extra_dirs or []:
        p = Path(d)
        if not p.is_dir():
            raise FileNotFoundError(f"extra_dirs 模板目录不存在: {p}")
        loaders.append(FileSystemLoader(str(p)))
    loaders.append(FileSystemLoader(str(TEMPLATES_DIR)))
    env = Environment(
        loader=ChoiceLoader(loaders),
        autoescape=True,
        undefined=StrictUndefined,
        trim_blocks=True,
        lstrip_blocks=True,
    )
    # PALETTE 键带连字符 (如 success-light), Jinja 点语法取不到 —— 注入下划线
    # 别名视图, 模板里 palette.success_light 与 palette['success-light'] 均可用。
    env.globals["palette"] = {**PALETTE,
                              **{k.replace("-", "_"): v for k, v in PALETTE.items()}}
    # base_css 是可信静态样式 (内含 CSS 单引号, 不能走 autoescape), 以 Markup 注入
    env.globals["base_css"] = Markup(get_base_css())
    env.filters["fmt_int"] = _fmt_int
    return env


@lru_cache(maxsize=1)
def get_env() -> Environment:
    """共享目录模板环境 (进程内复用; 管线 extra_dirs 场景用 make_env 自建)。"""
    return make_env()


def render(template_name: str, **context: Any) -> str:
    """按名渲染共享目录下的模板。缺变量/模板不存在均直接抛错 (fail-loud)。"""
    context.setdefault("generated_at", time.strftime("%Y-%m-%d %H:%M:%S"))
    return get_env().get_template(template_name).render(**context)


def list_templates() -> "dict[str, Path]":
    """共享模板目录清单 (测试/审计用: 确认模板文件真实存在)。"""
    if not TEMPLATES_DIR.is_dir():
        return {}
    return {p.name: p for p in sorted(TEMPLATES_DIR.glob("*.j2"))}
