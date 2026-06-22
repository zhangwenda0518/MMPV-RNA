#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
📊 SCI 级宏观组学数据可视化引擎 (精修美化终极版)
核心优化：
1. 折线图 (Panel A) 审美重构：科学蓝主色调，白底空心标记点，通透感增强。
2. 环形图 (Panel B) 标签重构：强制将 数据库名、比例、数量 三合一居中印在彩色圆环上。
3. 文本换行与防遮挡设计完美保留。
"""

import os
import sys
if sys.platform == 'win32':
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')

import re
import math
import argparse
import textwrap
import pandas as pd
import matplotlib
matplotlib.use('Agg')  # non-interactive backend, avoids Tk memory issues
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib as mpl

# ==========================================
# 0. SCI 期刊全局格式设置与安全回退
# ==========================================
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.sans-serif'] = ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"]
mpl.rcParams['font.family'] = "sans-serif"
mpl.rcParams['axes.spines.right'] = False
mpl.rcParams['axes.spines.top'] = False

NULL_WORDS = ['not_provided', 'nan', 'none', 'unknown', 'missing', '', 'not applicable', 'not collected']

def clean_series(series):
    return series.dropna()[~series.dropna().astype(str).str.strip().str.lower().isin(NULL_WORDS)]

def wrap_labels(labels, width=35):
    return [textwrap.fill(str(label), width=width) for label in labels]

def add_bar_labels(ax, values, is_horizontal=True):
    max_val = max(values) if len(values) > 0 else 1
    offset = max_val * 0.02 
    for i, v in enumerate(values):
        if is_horizontal:
            ax.text(v + offset, i, str(v), va='center', ha='left', fontsize=12, color='black')
        else:
            ax.text(i, v + offset, str(v), va='bottom', ha='center', fontsize=12, color='black')

def plot_sci_landscape(csv_path, output_dir="SCI_Figures_Output"):
    os.makedirs(output_dir, exist_ok=True)
    print(f"📥 正在读取数据: {csv_path}")
    df = pd.read_csv(csv_path)

    # 兼容 info Core13 输出 (CenterName / 无 Database)
    if 'Database' not in df.columns and 'Run' in df.columns:
        df['Database'] = df['Run'].astype(str).str.extract(r'^([A-Za-z]+)')[0].map(
            {'SRR': 'SRA', 'ERR': 'SRA', 'DRR': 'SRA'}).fillna('GSA')
    if 'CenterName' in df.columns and 'Organization_CenterName' not in df.columns:
        df['Organization_CenterName'] = df['CenterName']

    # 1. 深度数据聚合与纠错清洗
    df['ReleaseDate'] = pd.to_datetime(df['ReleaseDate'], errors='coerce')
    df['Year'] = df['ReleaseDate'].dt.year.fillna(0).astype(int)
    
    if 'Organization_CenterName' in df.columns:
        orgs = df['Organization_CenterName'].astype(str).str.strip().str.title()
        orgs = orgs.str.replace('Unversity', 'University', flags=re.IGNORECASE)
        orgs = orgs.str.replace('&Amp;', '&', flags=re.IGNORECASE)
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

    # ==========================================
    # 2. 动态构建总画布 (3x2 完美布局)
    # ==========================================
    fig, axes = plt.subplots(3, 2, figsize=(14, 16), dpi=150)
    ax_list = axes.flatten()
    colors_db = {'SRA': '#4C72B0', 'GSA': '#C44E52'} 
    plot_idx = 0

    # ------------------------------------------
    # Panel A: Temporal Distribution (高颜值美化版)
    # ------------------------------------------
    print("📊 绘制 A: 时间分布趋势图...")
    ax = ax_list[plot_idx]
    df_valid_years = df[df['Year'] > 2000] 
    year_counts = df_valid_years['Year'].value_counts().sort_index()
    
    # 使用经典科学蓝 (#0072B2)，并将数据点改为“白底蓝边”的空心样式，提升通透感
    ax.plot(year_counts.index, year_counts.values, marker='o', linestyle='-', linewidth=3.5, 
            markersize=10, color='#0072B2', markerfacecolor='white', markeredgewidth=2.5, zorder=3)
    # 淡蓝色面积填充
    ax.fill_between(year_counts.index, year_counts.values, color='#0072B2', alpha=0.15, zorder=2)
    ax.grid(axis='y', linestyle='--', alpha=0.6, color='#D3D3D3', zorder=1)
    
    for x, y in zip(year_counts.index, year_counts.values):
        ax.text(x, y + (max(year_counts.values)*0.03), str(y), ha='center', va='bottom', fontsize=13, fontweight='bold', color='#333333')
        
    ax.set_title('A. Temporal Distribution of Sequencing Data', loc='left', fontsize=20, fontweight='bold')
    ax.set_xlabel('Release Year', fontsize=16)
    ax.set_ylabel('Number of Runs', fontsize=16)
    ax.tick_params(axis='x', rotation=45, labelsize=13)
    ax.set_xticks(year_counts.index)
    plot_idx += 1

    # ------------------------------------------
    # Panel B: Database Proportion (强制圈内标记版)
    # ------------------------------------------
    print("📊 绘制 B: 数据库比例环形图...")
    ax = ax_list[plot_idx]
    db_counts = df['Database'].value_counts()
    
    # 手动拼装多行文字：名称 + 比例 + 数量
    custom_labels = []
    total = db_counts.sum()
    for name, val in db_counts.items():
        pct = val / total * 100
        custom_labels.append(f"{name}\n{pct:.1f}%\n(n={val})")

    # 绘制甜甜圈（不放外部标签，手动定位到环形中心）
    wedges, _ = ax.pie(
        db_counts,
        labels=None,
        startangle=140,
        colors=[colors_db.get(x, '#555555') for x in db_counts.index],
        wedgeprops=dict(width=0.55, edgecolor='w', linewidth=3),
        radius=1.3, center=(0, 0),
    )

    # 将标签精确放置在每个楔形的环形中心
    ring_center = 1.0 - 0.55 / 2  # 环形中心 = 外半径 - 环宽/2
    label_r = ring_center * 1.3    # 实际径向距离

    for i, wedge in enumerate(wedges):
        ang = math.radians((wedge.theta1 + wedge.theta2) / 2)
        x = label_r * math.cos(ang)
        y = label_r * math.sin(ang)
        ax.text(x, y, custom_labels[i], ha='center', va='center',
                fontsize=15, fontweight='bold', color='white')

    ax.set_title('B. Proportion of Data Origin', loc='left', fontsize=20, fontweight='bold', pad=40)
    plot_idx += 1

    # ------------------------------------------
    # Panel C: Top Organizations
    # ------------------------------------------
    print("📊 绘制 C: 核心贡献机构...")
    ax = ax_list[plot_idx]
    org_counts = clean_series(df['Organization_CenterName']).value_counts().head(10).sort_values(ascending=True)
    wrapped_org_index = wrap_labels(org_counts.index, width=35)
    
    ax.hlines(y=wrapped_org_index, xmin=0, xmax=org_counts.values, color='skyblue', linewidth=4)
    ax.plot(org_counts.values, wrapped_org_index, "o", markersize=12, color='royalblue')
    add_bar_labels(ax, org_counts.values, is_horizontal=True)
    
    ax.set_title('C. Top 10 Contributing Organizations', loc='left', fontsize=20, fontweight='bold')
    ax.set_xlabel('Number of Runs', fontsize=16)
    ax.tick_params(axis='y', labelsize=12)
    plot_idx += 1

    # ------------------------------------------
    # Panel D: Tissues
    # ------------------------------------------
    print("📊 绘制 D: 研究部位偏好...")
    ax = ax_list[plot_idx]
    if 'Tissue' in df.columns:
        tissue_counts = clean_series(df['Tissue']).value_counts().head(10).sort_values(ascending=False)
        if not tissue_counts.empty:
            wrapped_tissue_index = wrap_labels(tissue_counts.index, width=25)
            sns.barplot(x=tissue_counts.values, y=wrapped_tissue_index, hue=wrapped_tissue_index, palette="magma", legend=False, ax=ax)
            add_bar_labels(ax, tissue_counts.values, is_horizontal=True)
    
    ax.set_title('D. Top Investigated Biological Tissues', loc='left', fontsize=20, fontweight='bold')
    ax.set_xlabel('Number of Runs', fontsize=16)
    ax.tick_params(axis='y', labelsize=13)
    plot_idx += 1

    # ------------------------------------------
    # Panel E: Geographical Distribution
    # ------------------------------------------
    print("📊 绘制 E: 采样地理分布...")
    ax = ax_list[plot_idx]
    if 'Location_Clean' in df.columns:
        country_counts = clean_series(df['Location_Clean']).value_counts().head(10).sort_values(ascending=True)
        if not country_counts.empty:
            wrapped_country_index = wrap_labels(country_counts.index, width=30)
            ax.hlines(y=wrapped_country_index, xmin=0, xmax=country_counts.values, color='lightcoral', linewidth=4)
            ax.plot(country_counts.values, wrapped_country_index, "o", markersize=12, color='firebrick')
            add_bar_labels(ax, country_counts.values, is_horizontal=True)
    
    ax.set_title('E. Top Sample Collection Regions', loc='left', fontsize=20, fontweight='bold')
    ax.set_xlabel('Number of Runs', fontsize=16)
    ax.tick_params(axis='y', labelsize=13)
    plot_idx += 1
    
    # ------------------------------------------
    # Panel F: Age & Growth Stage
    # ------------------------------------------
    print("📊 绘制 F: 发育时期偏好...")
    ax = ax_list[plot_idx]
    if 'Age_GrowthStage_Clean' in df.columns:
        stage_counts = clean_series(df['Age_GrowthStage_Clean']).value_counts().head(10).sort_values(ascending=False)
        if not stage_counts.empty:
            wrapped_stage_index = wrap_labels(stage_counts.index, width=35)
            sns.barplot(x=stage_counts.values, y=wrapped_stage_index, hue=wrapped_stage_index, palette="crest", legend=False, ax=ax)
            add_bar_labels(ax, stage_counts.values, is_horizontal=True)
    
    ax.set_title('F. Top Developmental Stages / Ages', loc='left', fontsize=20, fontweight='bold')
    ax.set_xlabel('Number of Runs', fontsize=16)
    ax.tick_params(axis='y', labelsize=12)

    # 整体排版调优
    plt.tight_layout(pad=4.0, h_pad=6.0, w_pad=4.0)

    # 4. 输出总拼图
    out_pdf = os.path.join(output_dir, "Combined_Landscape_Full.pdf")
    out_png = os.path.join(output_dir, "Combined_Landscape_Full.png")
    
    fig.savefig(out_pdf, format='pdf', bbox_inches='tight')
    fig.savefig(out_png, format='png', dpi=150, bbox_inches='tight')
    
    print(f"\n🎉 完美收工！所有图表已存放至文件夹: [ {output_dir}/ ]")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="SCI 级宏观组学数据可视化引擎 (精修版)")
    parser.add_argument("-i", "--input", required=True, help="输入的 SRA_GSA_Merged_Final.csv 文件路径")
    parser.add_argument("-o", "--outdir", default="SCI_Figures_Output", help="图表输出的文件夹名称")
    
    args = parser.parse_args()
    
    if not os.path.exists(args.input):
        print(f"❌ 找不到文件: {args.input}")
    else:
        plot_sci_landscape(args.input, args.outdir)
