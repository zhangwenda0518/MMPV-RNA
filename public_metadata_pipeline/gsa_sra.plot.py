#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
📊 SCI 级宏观组学数据可视化引擎 (2026-10-14 与平台版 landscape_plot 对拍同步)

本次从植物病毒分析平台 Virus_Platform_Core/public_meta/landscape_plot.py 回移
（该文件 2026-10 由本脚本移植出去后又硬化演化，实测优于本版，故整体对齐）：

1. 面板函数化：合并总览图按 PANELS 顺序排布 3x2（只出主图，不出独立小图
   ——2026-10-14 按需求裁掉，平台版保留该能力供其 GUI 分开展示）。
2. Panel A 自适应：观测年份 ≤12 画柱状图（离散年份、标注零碰撞——折线版
   实测会被线穿过标注）；>12 保留折线+填充体现趋势，刻度交给 MaxNLocator。
3. Panel B 补回 total==0 保护（0 行 CSV 时 pie 抛 "All wedge sizes are zero"，
   空表对拍实抓平台版丢过这道保护）；单一来源改印环心——印在 360° 圆环上
   会被坐标区裁掉 (n=xxx) 尾巴（三七 836 条实测）。
4. CSV 读取改 utf-8-sig：平台导出的 CSV 带 BOM，不指定编码时首列列名带
   '\ufeff' 前缀，列匹配静默失败。
5. 配色升级 Okabe-Ito（色盲友好、期刊惯例），中文字体仍置前（机构/部位
   等类目含中文）。

与平台版的差异（有意保留）：不依赖 Virus_Platform_Core.sci_palette，色板
常量内联；rcParams 在本文件内置，独立脚本可直跑。
"""

# ── 图件四格式导出 (SVG / PNG / PDF / PNG@1200dpi) — 默认开启 ──────────
# 2026-10-20: 所有出图脚本统一口径; 用 MMPV_FIG_EXPORT=0 可整体关闭,
# MMPV_FIG_FORMATS=svg,png,pdf,hires 可裁剪, MMPV_FIG_HIRES_DPI 可改高清 dpi。
try:
    from mmpv_common.figure_export import install_savefig_hook as _install_fig_hook
except ImportError:  # 仓库根不在 sys.path (脚本被单独执行) → 向上找 mmpv_common/
    import os as _os
    import sys as _sys
    _d = _os.path.dirname(_os.path.abspath(__file__))
    while _d and not _os.path.isdir(_os.path.join(_d, "mmpv_common")):
        _parent = _os.path.dirname(_d)
        if _parent == _d:
            break
        _d = _parent
    if _d and _d not in _sys.path:
        _sys.path.insert(0, _d)
    from mmpv_common.figure_export import install_savefig_hook as _install_fig_hook
_install_fig_hook()
# ──────────────────────────────────────────────────────────────────────

import os
import sys
if sys.platform == 'win32' and sys.stdout is not None:
    try:
        sys.stdout.reconfigure(encoding='utf-8', errors='replace')
    except Exception:
        pass

import re
import math
import argparse
import textwrap
import pandas as pd
import matplotlib
matplotlib.use('Agg')  # non-interactive backend, avoids Tk memory issues
import matplotlib.pyplot as plt
import matplotlib as mpl

# ==========================================
# 0. SCI 期刊全局格式设置（自包含，无外部样式依赖）
# ==========================================
mpl.rcParams['pdf.fonttype'] = 42   # TrueType 嵌入，PDF 文字可编辑（期刊要求）
mpl.rcParams['ps.fonttype'] = 42
# CJK 字体置前 (机构/部位等类目含中文), 原西文字体保留作回退
mpl.rcParams['font.sans-serif'] = ["Microsoft YaHei", "SimHei", "Arial", "Helvetica", "DejaVu Sans", "sans-serif"]
mpl.rcParams['font.family'] = "sans-serif"
mpl.rcParams['axes.unicode_minus'] = False
mpl.rcParams['axes.spines.right'] = False
mpl.rcParams['axes.spines.top'] = False

NULL_WORDS = ['not_provided', 'nan', 'none', 'unknown', 'missing', '', 'not applicable', 'not collected']

# Okabe-Ito 色盲友好色板（与平台 sci_palette 同源）
OKABE_BLUE = '#0072B2'
OKABE_VERMILLION = '#D55E00'
OKABE_GREEN = '#009E73'
_GREY = '#999999'
# 数据库来源色：SRA 蓝 / GSA 朱红
COLORS_DB = {'SRA': OKABE_BLUE, 'GSA': OKABE_VERMILLION}


def clean_series(series):
    return series.dropna()[~series.dropna().astype(str).str.strip().str.lower().isin(NULL_WORDS)]

def wrap_labels_capped(labels, width=26, max_lines=2, fontsize=10):
    """行数封顶的标签包装：前 max_lines-1 行取整行，末行按词截断加 …。

    修复点（check_figures 体检实抓）：Top10 机构名不换行直接放进 y 轴会
    把整列宽度让给单行标签、六个子图被压窄；且不能用 shorten(width×
    max_lines) 预截——按词折行的边界损耗仍可能折出第 3 行，必须按折行
    结果强制截断末行才可靠。
    """
    out = []
    for s in labels:
        s = str(s).strip() or '(unknown)'
        lines = textwrap.wrap(s, width=width)
        if len(lines) > max_lines:
            kept = lines[:max_lines - 1]
            rest = s[len(' '.join(kept)):].strip()
            if len(rest) > width - 1:
                rest = rest[:width - 2].rstrip() + '…'
            lines = kept + [rest]
        out.append('\n'.join(lines))
    need_small = any(l.count('\n') >= max_lines for l in out)
    return out, (fontsize - 1 if need_small else fontsize)

def _note_empty(ax, msg='No valid data'):
    """面板无可用数据时给出居中占位提示，而不是留一张空坐标轴。"""
    ax.set_xticks([])
    ax.set_yticks([])
    ax.text(0.5, 0.5, msg, ha='center', va='center', transform=ax.transAxes,
            fontsize=14, color='#888888', style='italic')

def _prep(df):
    """列名兼容与深度清洗（时间/机构/组织/地理/发育时期），返回副本。"""
    df = df.copy()
    # 兼容不同来源的列名 (info Core14 / search 旧版 / search 新版)
    if 'Database' not in df.columns and 'Run' in df.columns:
        df['Database'] = df['Run'].astype(str).str.extract(r'^([A-Za-z]+)')[0].map(
            {'SRR': 'SRA', 'ERR': 'SRA', 'DRR': 'SRA'}).fillna('GSA')
    if 'CenterName' in df.columns:
        df['Organization_CenterName'] = df['CenterName']
    elif 'Organization_CenterName' not in df.columns:
        df['Organization_CenterName'] = pd.NA

    # 深度数据聚合与纠错清洗
    if 'ReleaseDate' in df.columns:
        df['ReleaseDate'] = pd.to_datetime(df['ReleaseDate'], errors='coerce')
        df['Year'] = df['ReleaseDate'].dt.year.fillna(0).astype(int)

    if 'Organization_CenterName' in df.columns:
        orgs = df['Organization_CenterName'].astype(str).str.strip().str.title()
        orgs = orgs.str.replace('Unversity', 'University', flags=re.IGNORECASE)
        orgs = orgs.str.replace('&Amp;', '&', flags=re.IGNORECASE)
        orgs = orgs.str.replace(r'\s+', ' ', regex=True)  # collapse whitespace
        df['Organization_CenterName'] = orgs

    if 'Tissue' in df.columns:
        df['Tissue'] = df['Tissue'].astype(str).str.strip().str.title()

    if 'Location' in df.columns:
        def clean_loc(val):
            val = str(val).strip(' "')
            if val.lower() in NULL_WORDS or 'missing' in val.lower():
                return pd.NA
            val = val.replace(':', ', ')
            val = re.sub(r'_Ai$', '', val, flags=re.IGNORECASE)
            val = re.sub(r'\s+', ' ', val)
            return val.title().strip()
        df['Location_Clean'] = df['Location'].apply(clean_loc)

    if 'Age_GrowthStage' in df.columns:
        def clean_age(val):
            if pd.isna(val): return pd.NA
            parts = [p.strip().title() for p in str(val).split('|') if p.strip().lower() not in NULL_WORDS]
            return " | ".join(parts) if parts else pd.NA
        df['Age_GrowthStage_Clean'] = df['Age_GrowthStage'].apply(clean_age)
    return df

# ==========================================
# 1. 六个面板的独立绘制函数（合并图与单图共用）
# ==========================================
def _pnl_a_temporal(ax, df):
    """A. 时间分布趋势。年份少（≤12）画柱状图——稀疏折线加大面积填充
    会把整个坐标区糊成色块，且数值标注会被折线穿过；年份多时保留
    折线+填充体现趋势。"""
    ax.set_title('A. Temporal Distribution of Sequencing Data', loc='left',
                 fontsize=18, fontweight='bold')
    ax.set_xlabel('Release Year', fontsize=14)
    ax.set_ylabel('Number of Runs', fontsize=14)
    if 'Year' not in df.columns:
        _note_empty(ax, 'No release dates')
        return
    year_counts = df[df['Year'] > 2000]['Year'].value_counts().sort_index()
    if year_counts.empty:
        _note_empty(ax, 'No valid release dates')
        return
    years = list(year_counts.index)
    vals = list(year_counts.values)
    if len(years) <= 12:
        ax.bar(range(len(years)), vals, width=0.6, color=OKABE_BLUE, zorder=3)
        for i, v in enumerate(vals):
            ax.text(i, v + max(vals) * 0.02, str(v), ha='center', va='bottom',
                    fontsize=13, fontweight='bold', color='#333333', zorder=4)
        ax.set_xticks(range(len(years)))
        ax.set_xticklabels([str(y) for y in years], fontsize=12)
        ax.set_ylim(0, max(vals) * 1.15)
    else:
        ax.plot(years, vals, marker='o', linestyle='-', linewidth=3.5,
                markersize=10, color=OKABE_BLUE, markerfacecolor='white',
                markeredgewidth=2.5, zorder=3)
        ax.fill_between(years, vals, color=OKABE_BLUE, alpha=0.15, zorder=2)
        for x, y in zip(years, vals):
            ax.text(x, y + (max(vals) * 0.03), str(y), ha='center', va='bottom',
                    fontsize=13, fontweight='bold', color='#333333')
        # 年份多时交给自动刻度定位（每年一 ticks 会挤成一团并把绘图区压窄）
        ax.xaxis.set_major_locator(mpl.ticker.MaxNLocator(nbins=8, integer=True))
        ax.tick_params(axis='x', labelsize=12)
    ax.grid(axis='y', linestyle='--', alpha=0.6, color='#D3D3D3', zorder=1)

def _pnl_b_database(ax, df):
    """B. 数据来源比例环形图。"""
    db_counts = df['Database'].value_counts()
    total = db_counts.sum()
    if total == 0:
        # 0 行 CSV / 无有效 Database 值时 pie 抛 "All wedge sizes are zero"
        _note_empty(ax, 'No database metadata')
        ax.set_title('B. Proportion of Data Origin', loc='left',
                     fontsize=18, fontweight='bold', pad=30)
        return
    wedges, _ = ax.pie(
        db_counts,
        labels=None,
        startangle=140,
        colors=[COLORS_DB.get(x, _GREY) for x in db_counts.index],
        wedgeprops=dict(width=0.55, edgecolor='w', linewidth=3),
        radius=1.0,
    )
    if len(db_counts) == 1:
        # 单一来源：比例文字印在环心（印在 360° 圆环上会被坐标区裁掉尾巴）
        name = db_counts.index[0]
        ax.text(0, 0, f'{name}\n100%\n(n={total})', ha='center', va='center',
                fontsize=15, fontweight='bold', color=COLORS_DB.get(name, _GREY),
                linespacing=1.4)
    else:
        ring_center = 1.0 - 0.55 / 2
        label_r = ring_center * 0.95
        for i, wedge in enumerate(wedges):
            ang = math.radians((wedge.theta1 + wedge.theta2) / 2)
            arc_deg = abs(wedge.theta2 - wedge.theta1)
            fs = max(9, min(14, arc_deg * 0.06))
            name = db_counts.index[i]
            val = db_counts.values[i]
            pct = val / total * 100
            label = f"{name}\n{pct:.1f}%  (n={val})"
            ax.text(label_r * math.cos(ang), label_r * math.sin(ang),
                    label, ha='center', va='center',
                    fontsize=fs, fontweight='bold', color='white',
                    linespacing=1.2)
    ax.set_title('B. Proportion of Data Origin', loc='left', fontsize=18, fontweight='bold', pad=30)

def _pnl_c_organizations(ax, df):
    """C. 核心贡献机构 Top10。"""
    org_counts = clean_series(df['Organization_CenterName']).value_counts().head(10).sort_values(ascending=True)
    if org_counts.empty:
        _note_empty(ax, 'No organization metadata')
    else:
        n = len(org_counts)
        y_pos = range(n)
        ax.barh(y_pos, org_counts.values, height=0.65, color=OKABE_BLUE)
        for i, v in enumerate(org_counts.values):
            ax.text(v + max(org_counts.values)*0.02, i, str(v), va='center', fontsize=11)
        ax.set_yticks(list(y_pos))
        _lbls, _fs = wrap_labels_capped(org_counts.index, width=26, fontsize=10)
        ax.set_yticklabels(_lbls, fontsize=_fs, linespacing=1.3)
        ax.invert_yaxis()
    top_n = f'Top {len(org_counts)}' if 0 < len(org_counts) < 10 else 'Top 10'
    ax.set_title(f'C. {top_n} Contributing Organizations', loc='left', fontsize=16, fontweight='bold')
    ax.set_xlabel('Number of Runs', fontsize=14)

def _pnl_d_tissues(ax, df):
    """D. 研究部位偏好 Top10。"""
    if 'Tissue' in df.columns:
        tissue_counts = clean_series(df['Tissue']).value_counts().head(10).sort_values(ascending=True)
    else:
        tissue_counts = pd.Series(dtype=int)
    if tissue_counts.empty:
        _note_empty(ax, 'No tissue metadata')
    else:
        n = len(tissue_counts)
        y_pos = range(n)
        # 分类条目统一主色 (不用顺序色图映射分类, 期刊惯例)
        ax.barh(y_pos, tissue_counts.values, height=0.65, color=OKABE_BLUE, edgecolor='white')
        for i, v in enumerate(tissue_counts.values):
            ax.text(v + max(tissue_counts.values)*0.02, i, str(v), va='center', fontsize=11)
        ax.set_yticks(list(y_pos))
        _lbls, _fs = wrap_labels_capped(tissue_counts.index, width=30, fontsize=10)
        ax.set_yticklabels(_lbls, fontsize=_fs)
        ax.invert_yaxis()
    ax.set_title('D. Top Investigated Biological Tissues', loc='left', fontsize=16, fontweight='bold')
    ax.set_xlabel('Number of Runs', fontsize=14)

def _pnl_e_regions(ax, df):
    """E. 采样地理分布 Top10。"""
    if 'Location_Clean' in df.columns:
        loc_counts = clean_series(df['Location_Clean']).value_counts().head(10).sort_values(ascending=True)
    else:
        loc_counts = pd.Series(dtype=int)
    if loc_counts.empty:
        _note_empty(ax, 'No location metadata')
    else:
        n = len(loc_counts)
        y_pos = range(n)
        ax.barh(y_pos, loc_counts.values, height=0.65, color=OKABE_VERMILLION)
        for i, v in enumerate(loc_counts.values):
            ax.text(v + max(loc_counts.values)*0.02, i, str(v), va='center', fontsize=11)
        ax.set_yticks(list(y_pos))
        _lbls, _fs = wrap_labels_capped(loc_counts.index, width=24, fontsize=10)
        ax.set_yticklabels(_lbls, fontsize=_fs)
        ax.invert_yaxis()
    ax.set_title('E. Top Sample Collection Regions', loc='left', fontsize=16, fontweight='bold')
    ax.set_xlabel('Number of Runs', fontsize=14)

def _pnl_f_stages(ax, df):
    """F. 发育时期偏好 Top10。"""
    if 'Age_GrowthStage_Clean' in df.columns:
        stage_counts = clean_series(df['Age_GrowthStage_Clean']).value_counts().head(10).sort_values(ascending=True)
    else:
        stage_counts = pd.Series(dtype=int)
    if stage_counts.empty:
        _note_empty(ax, 'No developmental stage metadata')
    else:
        n = len(stage_counts)
        y_pos = range(n)
        ax.barh(y_pos, stage_counts.values, height=0.65, color=OKABE_GREEN, edgecolor='white')
        for i, v in enumerate(stage_counts.values):
            ax.text(v + max(stage_counts.values)*0.02, i, str(v), va='center', fontsize=11)
        ax.set_yticks(list(y_pos))
        _lbls, _fs = wrap_labels_capped(stage_counts.index, width=30, fontsize=10)
        ax.set_yticklabels(_lbls, fontsize=_fs, linespacing=1.2)
        ax.invert_yaxis()
    ax.set_title('F. Top Developmental Stages / Ages', loc='left', fontsize=16, fontweight='bold')
    ax.set_xlabel('Number of Runs', fontsize=14)

# 六个面板的绘制顺序（合并总览图按此排布 3x2）
PANELS = [_pnl_a_temporal, _pnl_b_database, _pnl_c_organizations,
          _pnl_d_tissues, _pnl_e_regions, _pnl_f_stages]

def plot_sci_landscape(csv_path, output_dir="SCI_Figures_Output"):
    """主入口：输出 3x2 合并总览图（PNG 600dpi + PDF 矢量）。"""
    os.makedirs(output_dir, exist_ok=True)
    print(f"📥 正在读取数据: {csv_path}")
    # utf-8-sig：平台导出的 CSV 带 BOM（Excel 直开不乱码），不指定编码时首列
    # 列名会带 '\ufeff' 前缀，导致列匹配静默失败
    df = _prep(pd.read_csv(csv_path, encoding='utf-8-sig'))

    # ==========================================
    # 2. 合并总画布 (3x2 完美布局)
    #    显示 dpi 150; 保存时 PNG 600 / PDF 矢量 (出版要求 ≥300dpi)
    # ==========================================
    fig, axes = plt.subplots(3, 2, figsize=(16, 20), dpi=150)
    ax_list = axes.flatten()
    for i, draw in enumerate(PANELS):
        print(f"📊 绘制 {chr(65 + i)}: ...")
        draw(ax_list[i], df)
    plt.tight_layout(pad=2.0, h_pad=3.0, w_pad=3.0)

    # 3. 输出总拼图 (PDF 矢量 + PNG 600dpi, 满足期刊 ≥300dpi 要求)
    out_pdf = os.path.join(output_dir, "Combined_Landscape_Full.pdf")
    out_png = os.path.join(output_dir, "Combined_Landscape_Full.png")
    fig.savefig(out_pdf, format='pdf', bbox_inches='tight')
    fig.savefig(out_png, format='png', dpi=600, bbox_inches='tight')
    plt.close(fig)
    print(f"✅ 合并总览图: {out_png}")

    print(f"\n🎉 完美收工！所有图表已存放至文件夹: [ {output_dir}/ ]")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="SCI 级宏观组学数据可视化引擎 (与平台版同步)")
    parser.add_argument("-i", "--input", required=True, help="输入的 SRA_GSA_Merged_Final.csv 文件路径")
    parser.add_argument("-o", "--outdir", default="SCI_Figures_Output", help="图表输出的文件夹名称")

    args = parser.parse_args()

    if not os.path.exists(args.input):
        print(f"❌ 找不到文件: {args.input}")
    else:
        plot_sci_landscape(args.input, args.outdir)
