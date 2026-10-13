#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
report_style.py — MMPV 七（八）管线 HTML 报告共享样式框架
==========================================================

所有管线的 HTML 报告通过 import 本模块获取统一的 CSS 与 HTML 组件，
确保跨管线视觉一致。用法：

    from mmpv_common.report_style import get_base_css, kpi_card, section_header

    css = get_base_css()
    html = f"<style>{css}</style>"
    html += kpi_card("检测病毒种数", "47", accent="#2563eb")
    html += section_header("检测 overview", icon="🔍")

组件的规范实现在 mmpv_common/templates/macros.html.j2 (Jinja2 宏,
报告模板内用 {% import "macros.html.j2" as m %} 直接调用);
本模块的组件函数是它的懒加载薄包装 —— Python 侧与模板侧单一实现来源,
输出一致。模板渲染统一走 mmpv_common/report_template.py。

色板：MMPV Blue (#2563eb) 为主色，中性灰底，语义色用于状态标记。
设计原则：纯 CSS（零外部依赖），可离线渲染，打印友好。
"""

# ── 色板 ────────────────────────────────────────────────────────────────
PALETTE = {
    "primary":       "#2563eb",
    "primary-hover": "#1d4ed8",
    "primary-light": "#eff4ff",
    "success":       "#059669",
    "success-light": "#ecfdf5",
    "warning":       "#d97706",
    "warning-light": "#fffbeb",
    "danger":        "#dc2626",
    "danger-light":  "#fef2f2",
    "info":          "#0891b2",
    "info-light":    "#ecfeff",
    "ink":           "#1a2332",
    "ink-secondary": "#5a6a7e",
    "surface":       "#ffffff",
    "surface-alt":   "#f6f8fb",
    "border":        "#e1e5eb",
    "border-light":  "#eef0f4",
    "sidebar-bg":    "#f3f5f8",
}

def _css_vars():
    pairs = "\n".join(f"  --{k}: {v};" for k, v in PALETTE.items())
    return f":root {{\n{pairs}\n}}"

# ── 共享 CSS ───────────────────────────────────────────────────────────
BASE_CSS = _css_vars() + """

/* ── Reset & base ── */
*, *::before, *::after { box-sizing: border-box; margin: 0; padding: 0; }
body {
  font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, 'Helvetica Neue', Arial, sans-serif;
  color: var(--ink); background: var(--surface-alt);
  line-height: 1.6; font-size: 14px;
}
a { color: var(--primary); text-decoration: none; }
a:hover { text-decoration: underline; }

/* ── Layout ── */
.container { max-width: 1200px; margin: 0 auto; padding: 24px 16px; }
.two-col { display: grid; grid-template-columns: 1fr 1fr; gap: 16px; }
.three-col { display: grid; grid-template-columns: 1fr 1fr 1fr; gap: 16px; }
@media (max-width: 768px) {
  .two-col, .three-col { grid-template-columns: 1fr; }
}

/* ── Cards ── */
.card {
  background: var(--surface); border: 1px solid var(--border);
  border-radius: 8px; padding: 20px; margin-bottom: 16px;
  box-shadow: 0 1px 3px rgba(0,0,0,.06);
}
.card-header {
  font-size: 16px; font-weight: 600; color: var(--ink);
  border-bottom: 1px solid var(--border-light); padding-bottom: 10px; margin-bottom: 14px;
}
.card-body { font-size: 14px; }

/* ── KPI stat cards ── */
.kpi-grid {
  display: grid; grid-template-columns: repeat(auto-fit, minmax(180px, 1fr));
  gap: 16px; margin-bottom: 24px;
}
.kpi-card {
  background: var(--surface); border: 1px solid var(--border);
  border-radius: 8px; padding: 20px 16px; text-align: center;
  box-shadow: 0 1px 3px rgba(0,0,0,.06);
  border-top: 3px solid var(--primary);
}
.kpi-card .kpi-value {
  font-size: 36px; font-weight: 700; color: var(--primary); line-height: 1.1;
}
.kpi-card .kpi-label {
  font-size: 12px; color: var(--ink-secondary); margin-top: 6px;
  text-transform: uppercase; letter-spacing: .5px;
}
.kpi-card.success { border-top-color: var(--success); }
.kpi-card.success .kpi-value { color: var(--success); }
.kpi-card.warning { border-top-color: var(--warning); }
.kpi-card.warning .kpi-value { color: var(--warning); }
.kpi-card.danger  { border-top-color: var(--danger); }
.kpi-card.danger  .kpi-value { color: var(--danger); }

/* ── Sidebar ── */
.layout { display: flex; min-height: 100vh; }
.sidebar {
  width: 220px; background: var(--sidebar-bg);
  border-right: 1px solid var(--border); padding: 16px 0;
  flex-shrink: 0;
}
.sidebar .brand {
  font-size: 16px; font-weight: 700; color: var(--primary);
  padding: 0 16px 14px; border-bottom: 1px solid var(--border); margin-bottom: 8px;
}
.sidebar a {
  display: block; padding: 8px 16px; color: var(--ink-secondary);
  font-size: 13px; border-left: 3px solid transparent;
}
.sidebar a:hover, .sidebar a.active {
  color: var(--primary); background: var(--primary-light);
  border-left-color: var(--primary);
}
.main-content { flex: 1; padding: 24px; overflow-x: auto; }
@media (max-width: 768px) {
  .layout { flex-direction: column; }
  .sidebar { width: 100%; border-right: none; border-bottom: 1px solid var(--border); }
}

/* ── Tables ── */
.table-wrap { overflow-x: auto; margin-bottom: 16px; }
table.data-table {
  width: 100%; border-collapse: collapse; font-size: 13px;
}
table.data-table th {
  background: var(--surface-alt); padding: 10px 12px; text-align: left;
  font-weight: 600; color: var(--ink-secondary);
  border-bottom: 2px solid var(--border); white-space: nowrap;
}
table.data-table td {
  padding: 8px 12px; border-bottom: 1px solid var(--border-light);
}
table.data-table tr:hover td { background: var(--primary-light); }

/* ── Badges / labels ── */
.badge {
  display: inline-block; padding: 3px 8px; border-radius: 4px;
  font-size: 11px; font-weight: 600; line-height: 1.4;
}
.badge-primary { background: var(--primary-light); color: var(--primary); }
.badge-success { background: var(--success-light); color: var(--success); }
.badge-warning { background: var(--warning-light); color: var(--warning); }
.badge-danger  { background: var(--danger-light);  color: var(--danger); }
.badge-info    { background: #ecfeff; color: var(--info); }

/* ── Chart containers ── */
.chart-container {
  position: relative; margin-bottom: 20px;
  background: var(--surface); border: 1px solid var(--border);
  border-radius: 8px; padding: 16px;
}
.chart-container .chart-title {
  font-size: 14px; font-weight: 600; margin-bottom: 12px; color: var(--ink);
}

/* ── Section headers ── */
.section-header {
  display: flex; align-items: center; gap: 8px;
  margin: 28px 0 16px; padding-bottom: 8px;
  border-bottom: 2px solid var(--primary);
}
.section-header .icon { font-size: 20px; }
.section-header h2 {
  font-size: 18px; font-weight: 700; color: var(--ink); margin: 0;
}

/* ── Progress bars ── */
.progress-bar {
  height: 8px; background: var(--border-light);
  border-radius: 4px; overflow: hidden; margin-top: 6px;
}
.progress-bar .fill {
  height: 100%; background: var(--primary); border-radius: 4px;
  transition: width .3s ease;
}

/* ── Print friendly ── */
@media print {
  .sidebar { display: none; }
  .card, .kpi-card { box-shadow: none; border: 1px solid #ccc; }
  body { font-size: 12px; }
}
"""

def get_base_css():
    """返回完整的基础 CSS 字符串（含 CSS 变量定义）。"""
    return BASE_CSS


# ── HTML 组件快捷函数 (规范实现 = templates/macros.html.j2, 此处为懒加载薄包装) ──
# 懒加载: report_template 顶层 import 本模块取 PALETTE, 函数内反向 import 避免循环。

def _macros():
    from mmpv_common.report_template import get_env
    return get_env().get_template("macros.html.j2").module


def kpi_card(label: str, value, accent: str = "", extra: str = "") -> str:
    """生成一个 KPI 大数字卡片。
    accent: "", "success", "warning", "danger"
    extra:  附加可信 HTML（如趋势箭头、副标签; 原样注入不转义）
    """
    return str(_macros().kpi_card(label=label, value=value, accent=accent,
                                  extra_html=extra))


def kpi_grid(cards_html: str) -> str:
    """将多张 KPI 卡包裹在网格容器中。cards_html 为已构建的可信 HTML。"""
    return str(_macros().kpi_grid(cards_html))


def section_header(title: str, icon: str = "") -> str:
    return str(_macros().section_header(title=title, icon=icon))


def data_table(headers: list, rows: list, caption: str = "",
               table_id: str = "") -> str:
    """生成标准数据表。headers = 列名列表，rows = 值列表的列表。
    值经 autoescape 转义; table_id 供交互 JS (getElementById) 挂钩。"""
    return str(_macros().data_table(headers=headers, rows=rows, caption=caption,
                                    table_id=table_id))


def badge(text: str, kind: str = "primary") -> str:
    return str(_macros().badge(text=text, kind=kind))
