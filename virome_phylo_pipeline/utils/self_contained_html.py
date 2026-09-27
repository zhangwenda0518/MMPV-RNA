#!/usr/bin/env python3
"""self_contained_html.py — 单文件自包含 HTML 图 (零外链) 统一出口

E3 (2026-09-16)。背景: 管线此前的交互式 HTML 靠外部 CDN (unpkg leaflet、
plotly CDN、CartoDB 瓦片), 离线机器/内网直接白屏 —— 而且白屏**不报任何错**,
属典型静默降级。本模块提供统一出口, 并把「有没有残留外链」变成可检查的事实。

要点 (踩过的坑):
  · plotly 图: 必须 `include_plotlyjs=True` —— 从**已安装的 plotly 包**里把
    plotly.min.js 正文内联。写成 `include_plotlyjs="<某路径>"` 时 plotly 会生成
    `src="D:\\...\\plotly.min.js"` 绝对路径引用, 换台机器即白屏, 不算自包含。
  · 写完必须自检: `external_refs()` 扫 src/href/url() 里的 http(s):// 与协议相对
    //host 引用; 有残留逐条 log.warning (或 strict=True 直接抛错)。

无 plotly 时 fail-loud (ImportError + 修复提示), 不退化成 CDN 引用。
"""

import os
import re
from typing import List, Optional

# 协议相对 (//cdn.example.com/x.js) 也要抓 —— 离线时同样是白屏原因
_EXT_REF_RE = re.compile(
    r"""(?:(?:src|href)\s*=\s*["'](?P<attr>(?:https?:)?//[^"']+)["'])"""
    r"""|(?:url\(\s*["']?(?P<css>(?:https?:)?//[^"')\s]+))""",
    re.IGNORECASE)

# 本机绝对路径 src/href (最常见: plotly 写成 src="D:\...\plotly.min.js") —— 换机即白屏
_LOCAL_REF_RE = re.compile(
    r"""(?:src|href)\s*=\s*["'](?P<local>(?:[A-Za-z]:[\\/]|file://|\\\\)[^"']*)["']""",
    re.IGNORECASE)

# 内联 plotly 正文里含 100+ 处「惰性」网址 (attribution / maplibre 图标 / mapbox 样式),
# 它们不会在页面加载时发请求, 扫外链时必须先把正文剔掉 —— 否则自检永远不通过,
# 反而会让人把「外链检查」当噪声关掉。见 external_refs(strip=...)。
def external_refs(html_text: str, strip: Optional[str] = None) -> List[str]:
    """返回 html 里所有**加载期**外部引用 (去重, 保持出现顺序)。空 = 真自包含。

    strip: 先把这段文本 (通常 = 内联的 plotly.min.js 正文) 从 html 里删掉再扫,
           用于排除库内部惰性网址的干扰。
    同时抓本机绝对路径引用 (src="D:\\...\\plotly.min.js" 这种「假自包含」)。
    """
    body = html_text or ""
    if strip:
        body = body.replace(strip, "", 1)
    out, seen = [], set()
    for m in _EXT_REF_RE.finditer(body):
        ref = m.group("attr") or m.group("css") or ""
        if ref and ref not in seen:
            seen.add(ref)
            out.append(ref)
    for m in _LOCAL_REF_RE.finditer(body):
        ref = "(本机绝对路径) " + m.group("local")
        if ref not in seen:
            seen.add(ref)
            out.append(ref)
    return out


def write_html(html_text: str, output_path: str, log=None,
               strict: bool = False, what: str = "HTML",
               strip: Optional[str] = None) -> str:
    """落盘 + 外链自检。strict=True 时发现残留外链直接抛错, 否则逐条告警。"""
    os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
    with open(output_path, "w", encoding="utf-8") as f:
        f.write(html_text)
    refs = external_refs(html_text, strip=strip)
    if refs:
        msg = (f"{what} 仍含 {len(refs)} 处外部引用, 离线会白屏: "
               + "; ".join(refs[:6]) + ("…" if len(refs) > 6 else ""))
        if strict:
            raise RuntimeError(msg)
        if log is not None:
            log.warning(f"  ⚠ {msg}")
    elif log is not None:
        log.emit(f"  ✓ {what} 自包含 (零外链): {os.path.basename(output_path)}")
    return output_path


def plotlyjs_source() -> str:
    """内联用的 plotly.min.js 正文。无 plotly 包时 fail-loud。"""
    try:
        from plotly.offline.offline import get_plotlyjs
    except ImportError as e:  # pragma: no cover - 环境缺 plotly
        raise ImportError(
            "自包含 HTML 需要 plotly 包 (用于内联 plotly.min.js)。"
            f"当前环境导入失败: {e}。修复: pip install plotly"
        ) from e
    js = get_plotlyjs()
    if "</script" in js:  # pragma: no cover - plotly 包异常才会发生
        raise RuntimeError("plotly.min.js 含 </script 片段, 无法安全内联")
    return js


def figure_to_html(fig, title: str = "", note: str = "", height: int = None,
                   config: Optional[dict] = None) -> str:
    """plotly Figure → 单文件自包含 HTML 文本 (plotly.min.js 正文内联)。

    note: 追加到图下方的一段说明 (数据自检/口径提示, 见各调用点)。
    """
    cfg = {"displayModeBar": False, "responsive": True}
    cfg.update(config or {})
    html = fig.to_html(include_plotlyjs=True, full_html=True, config=cfg,
                       default_height=height or None)
    if note:
        banner = (f'<div style="max-width:1180px;margin:8px auto 24px;padding:8px 14px;'
                  f'background:#fdf6e3;border-left:3px solid #f39c12;color:#7f5a04;'
                  f'font:12px/1.7 -apple-system,Segoe UI,Arial,sans-serif">{note}</div>')
        html = html.replace("</body>", banner + "</body>", 1)
    return html


def write_figure(fig, output_path: str, title: str = "", note: str = "",
                 height: int = None, config: Optional[dict] = None,
                 log=None, strict: bool = False) -> str:
    """plotly Figure → 自包含 HTML 落盘 (+外链自检)。返回输出路径。"""
    html = figure_to_html(fig, title=title, note=note, height=height, config=config)
    js = plotlyjs_source()
    if js not in html:  # include_plotlyjs=True 失效 → 立刻报, 不产出「看着像自包含」的文件
        raise RuntimeError("plotly.js 未被内联 (include_plotlyjs=True 未生效)")
    return write_html(html, output_path, log=log, strict=strict,
                      what=f"图 {os.path.basename(output_path)}", strip=js)
