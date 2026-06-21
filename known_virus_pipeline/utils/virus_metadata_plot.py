import os
import sys
import argparse
import warnings
import re
import subprocess
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

warnings.filterwarnings("ignore", category=UserWarning, module="scipy")
warnings.filterwarnings("ignore", category=FutureWarning, module="seaborn")
try:
    warnings.filterwarnings("ignore", category=pd.errors.SettingWithCopyWarning)
except AttributeError:
    pass

try:
    # 兼容多平台字体，解决负号显示问题
    plt.rcParams['font.sans-serif'] = ['Arial Unicode MS', 'SimHei', 'Microsoft YaHei', 'DejaVu Sans']
    plt.rcParams['axes.unicode_minus'] = False
except: pass

def safe_filename(name):
    return re.sub(r'[^a-zA-Z0-9_\-]', '_', str(name)).strip('_')

def clean_ai_labels(val):
    """【新增】为了图表美观，去除 _AI 后缀，并统一缺失值为 Unknown"""
    s = str(val).strip()
    s = re.sub(r'_AI$', '', s, flags=re.IGNORECASE)
    if s.lower() in ['not_provided', 'unknown', 'nan', 'none', '', '<na>']:
        return 'Unknown'
    return s

class VirusMetadataAnalyzer:
    def __init__(self, virus_file, meta_file, out_dir):
        self.virus_file = virus_file
        self.meta_file = meta_file
        self.out_dir = out_dir
        os.makedirs(self.out_dir, exist_ok=True)
        
        # 全面包含 11 列核心元数据
        self.meta_features = [
            'ScientificName', 'BioProject', 'CenterName', 'Tissue', 'Source', 
            'Location', 'Age_GrowthStage', 'CollectionDate', 'LibrarySource', 'TaxID'
        ]
        
    def _clean_dataframe(self, df):
        # 【核心修复 1】: 强制剔除表头中可能存在的不可见 BOM 字符 (\ufeff) 和首尾空格
        df.columns = [str(c).replace('\ufeff', '').strip() for c in df.columns]
        for col in df.columns:
            if df[col].dtype == 'object':
                df[col] = df[col].astype(str).str.strip()
        return df

    def load_and_merge(self):
        print(f"📥 读取病毒鉴定结果: {self.virus_file}")
        df_virus = pd.read_csv(self.virus_file, sep=None, engine='python', skipinitialspace=True)
        df_virus = self._clean_dataframe(df_virus)
        
        if 'Covered%' in df_virus.columns:
            df_virus['Covered%'] = pd.to_numeric(df_virus['Covered%'], errors='coerce')
            df_virus = df_virus[df_virus['Covered%'] > 0]
            print(f"  ✂️ 已剔除 'Covered%' 为 0 的记录: 剩余 {len(df_virus)} 条。")

        # 【核心修复 2】: 病毒表自动主键对齐 (容错一切可能的名字)
        virus_rename_dict = {c: 'Run' for c in df_virus.columns if c.lower() in ['sample', 'sample_id', 'query', 'query_id', 'id', 'srr', 'run_accession']}
        if virus_rename_dict:
            df_virus = df_virus.rename(columns=virus_rename_dict)
            
        # 如果还是没有 Run，直接将第一列强制霸占为 Run
        if 'Run' not in df_virus.columns:
            first_col = df_virus.columns[0]
            print(f"  ⚠️ 病毒表未找到标准 'Run' 列，强制将第一列 '{first_col}' 设为主键 Run。")
            df_virus = df_virus.rename(columns={first_col: 'Run'})
            
        print(f"📥 读取 11维度大一统元数据: {self.meta_file}")
        df_meta = pd.read_csv(self.meta_file, sep=None, engine='python', skipinitialspace=True)
        df_meta = self._clean_dataframe(df_meta)
        
        # 【核心修复 3】: 元数据表主键容错
        meta_rename_dict = {c: 'Run' for c in df_meta.columns if c.lower() in ['run', 'run_accession', 'query_id', 'sample']}
        if meta_rename_dict:
            df_meta = df_meta.rename(columns=meta_rename_dict)
            
        if 'Run' not in df_meta.columns:
            first_col = df_meta.columns[0]
            print(f"  ⚠️ 元数据表未找到标准 'Run' 列，强制将第一列 '{first_col}' 设为主键 Run。")
            df_meta = df_meta.rename(columns={first_col: 'Run'})
            
        # 清洗 AI 标签和空值
        for col in self.meta_features:
            if col in df_meta.columns:
                df_meta[col] = df_meta[col].apply(clean_ai_labels)
            else: 
                df_meta[col] = 'Unknown'
                
        print("🔗 正在进行强力联表匹配...")
        
        # 如果此时仍然报错，那将是不可能的（因为首列已经被强制接管为 Run）
        merged_by_query = pd.DataFrame()
        if 'query_id' in df_meta.columns:
            merged_by_query = pd.merge(df_virus, df_meta, left_on='Run', right_on='query_id', how='inner')
        
        if not merged_by_query.empty:
            unmatched_virus = df_virus[~df_virus['Run'].isin(merged_by_query['Run'])]
        else:
            unmatched_virus = df_virus.copy()
            
        merged_by_run = pd.merge(unmatched_virus, df_meta, on='Run', how='inner')
        self.df_merged = pd.concat([merged_by_query, merged_by_run], ignore_index=True)
        
        final_unmatched = df_virus[~df_virus['Run'].isin(self.df_merged['Run'])].copy()
        if not final_unmatched.empty:
            for col in df_meta.columns:
                if col not in final_unmatched.columns:
                    final_unmatched[col] = 'Unknown'
            self.df_merged = pd.concat([self.df_merged, final_unmatched], ignore_index=True)

        self.df_merged['ScientificName'] = self.df_merged['ScientificName'].replace('Unknown', 'Unknown_Host')
        
        matched_count = len(self.df_merged[self.df_merged['ScientificName'] != 'Unknown_Host'])
        unmatched_count = len(self.df_merged[self.df_merged['ScientificName'] == 'Unknown_Host'])
        print(f"✅ 合并完成！总关联记录: {len(self.df_merged)} (成功找到元数据: {matched_count}, 依然丢失: {unmatched_count})")

    def _save_plot_and_table(self, ct, out_dir, filename_prefix, fig, ax, kind='bar'):
        if kind == 'bar':
            for c in ax.containers:
                labels = [int(v.get_height()) if v.get_height() > (ct.max().max()*0.05) else '' for v in c]
                ax.bar_label(c, labels=labels, label_type='center', fontsize=9, color='white', weight='bold')
            totals = ct.sum(axis=1)
            for i, total in enumerate(totals):
                if total > 0: ax.text(i, total + (totals.max() * 0.01), f"{int(total)}", ha='center', va='bottom', weight='bold')
        elif kind == 'barh':
            for c in ax.containers:
                labels = [int(v.get_width()) if v.get_width() > (ct.max().max()*0.05) else '' for v in c]
                ax.bar_label(c, labels=labels, label_type='center', fontsize=9, color='white', weight='bold')
            totals = ct.sum(axis=1)
            for i, total in enumerate(totals):
                if total > 0: ax.text(total + (totals.max() * 0.01), i, f" {int(total)}", va='center', weight='bold')

        plt.tight_layout()
        plt.savefig(os.path.join(out_dir, f"{filename_prefix}.png"), dpi=300, bbox_inches='tight')
        plt.close()
        
        ct_export = ct.copy()
        ct_export['Total'] = ct_export.sum(axis=1)
        ct_export.to_csv(os.path.join(out_dir, f"{filename_prefix}.csv"), encoding='utf-8-sig')

    def plot_overall_summary(self):
        print("\n📊 模块 1: 正在绘制全局概览图...")
        summary_dir = os.path.join(self.out_dir, "01_Global_Summary")
        os.makedirs(summary_dir, exist_ok=True)

        sample_vc = self.df_merged.groupby(['Run', 'ScientificName'])['Taxonomy'].nunique().reset_index()
        sample_vc.columns = ['Run', 'ScientificName', 'Virus_Count']
        sample_vc['Count_Label'] = sample_vc['Virus_Count'].astype(str) + " Virus(es)"
        
        ct = pd.crosstab(sample_vc['ScientificName'], sample_vc['Count_Label'])
        fig, ax = plt.subplots(figsize=(12, 7))
        ct.plot(kind='bar', stacked=True, colormap='Set2', ax=ax, edgecolor='black', linewidth=0.5)
        plt.title('Global Co-infection Overview across Hosts', fontsize=16, pad=15)
        plt.xticks(rotation=45, ha='right')
        plt.legend(title='Concurrent Infections', bbox_to_anchor=(1.05, 1), loc='upper left')
        
        self._save_plot_and_table(ct, summary_dir, "01_Coinfection_Counts_by_Host", fig, ax, kind='bar')

        combos = self.df_merged.groupby('Run')['Taxonomy'].apply(lambda x: ' + \n'.join(sorted(x))).reset_index()
        combo_counts = combos['Taxonomy'].value_counts().head(15) 
        
        plt.figure(figsize=(12, 8))
        ax = sns.barplot(x=combo_counts.values, y=combo_counts.index, hue=combo_counts.index, legend=False, palette='Spectral')
        plt.title('Top 15 Viral Infection Combinations (Global)', fontsize=16, pad=15)
        for i, v in enumerate(combo_counts.values):
            ax.text(v + (combo_counts.max() * 0.01), i, str(int(v)), color='black', va='center', weight='bold', fontsize=11)
        plt.tight_layout()
        plt.savefig(os.path.join(summary_dir, "02_Top_Viral_Combinations.png"), dpi=300, bbox_inches='tight')
        plt.close()
        
        combo_counts.reset_index().rename(columns={'Taxonomy': 'Combination', 'count': 'Sample_Count'}).to_csv(
            os.path.join(summary_dir, "02_Top_Viral_Combinations.csv"), index=False, encoding='utf-8-sig'
        )

    def plot_global_features_vs_virus(self):
        print("📊 模块 2: 正在绘制 全局特征 vs 病毒谱...")
        feature_dir = os.path.join(self.out_dir, "02_Global_Features_vs_Viruses")
        os.makedirs(feature_dir, exist_ok=True)
        
        top_viruses = self.df_merged['Taxonomy'].value_counts().head(15).index
        df_filtered = self.df_merged[self.df_merged['Taxonomy'].isin(top_viruses)].copy()
        
        for feature_col in self.meta_features:
            if feature_col not in df_filtered.columns: continue
            
            top_features = df_filtered[feature_col].value_counts().head(10).index
            df_plot = df_filtered.copy()
            df_plot[feature_col] = df_plot[feature_col].apply(lambda x: x if x in top_features else 'Other')
            
            ct = pd.crosstab(df_plot['Taxonomy'], df_plot[feature_col])
            ct['Total'] = ct.sum(axis=1)
            ct = ct.sort_values('Total', ascending=True).drop(columns='Total') 
            
            fig, ax = plt.subplots(figsize=(14, 8))
            ct.plot(kind='barh', stacked=True, ax=ax, colormap='tab20', edgecolor='black', linewidth=0.5)
            
            plt.title(f'Viral Taxonomy Distribution across {feature_col}', fontsize=16, pad=15)
            plt.ylabel('Viral Taxonomy (Top 15)', fontsize=14)
            plt.xlabel('Number of Detected Infections', fontsize=14)
            plt.legend(title=feature_col, bbox_to_anchor=(1.05, 1), loc='upper left')
            
            idx_num = str(self.meta_features.index(feature_col) + 1).zfill(2)
            self._save_plot_and_table(ct, feature_dir, f"{idx_num}_Virus_vs_{feature_col}", fig, ax, kind='barh')

    def plot_infection_complexity_breakdown(self):
        print("📊 模块 3: 正在精细拆解 1种、2种、3种+ 病毒组合情况...")
        breakdown_dir = os.path.join(self.out_dir, "03_Infection_Complexity_Breakdown")
        os.makedirs(breakdown_dir, exist_ok=True)
        
        sample_info = self.df_merged.groupby('Run').agg({
            'Taxonomy': lambda x: sorted(list(set(x))),
            'ScientificName': 'first'
        }).reset_index()
        
        sample_info['Complexity'] = sample_info['Taxonomy'].apply(len)
        sample_info['Combo_Name'] = sample_info['Taxonomy'].apply(lambda x: ' + \n'.join(x))
        sample_info['Level'] = sample_info['Complexity'].apply(lambda x: x if x <= 2 else '3+')
        
        for level in [1, 2, '3+']:
            df_lvl = sample_info[sample_info['Level'] == level]
            if df_lvl.empty: continue
            
            top_combos = df_lvl['Combo_Name'].value_counts().head(15).index
            df_lvl_top = df_lvl[df_lvl['Combo_Name'].isin(top_combos)]
            if df_lvl_top.empty: continue
            
            ct = pd.crosstab(df_lvl_top['Combo_Name'], df_lvl_top['ScientificName'])
            ct['Total'] = ct.sum(axis=1)
            ct = ct.sort_values('Total', ascending=True).drop(columns='Total')
            
            fig, ax = plt.subplots(figsize=(14, max(6, len(ct)*0.8)))
            ct.plot(kind='barh', stacked=True, ax=ax, colormap='tab20', edgecolor='black', linewidth=0.5)
            
            title_prefix = "Single Infections" if level == 1 else f"{level}-Virus Co-infections"
            plt.title(f'{title_prefix} Breakdown by Host', fontsize=16, pad=15)
            plt.ylabel('Viral Combination', fontsize=14)
            plt.xlabel('Number of Samples', fontsize=14)
            plt.legend(title='Host (ScientificName)', bbox_to_anchor=(1.05, 1), loc='upper left')
            
            self._save_plot_and_table(ct, breakdown_dir, f"Level_{level}_Breakdown", fig, ax, kind='barh')

    def plot_virus_specific_profiles(self):
        print("📊 模块 4: 正在为 Top 8 核心病毒生成专属特写画像及数据表...")
        profile_dir = os.path.join(self.out_dir, "04_Virus_Specific_Profiles")
        os.makedirs(profile_dir, exist_ok=True)
        
        top_viruses = self.df_merged['Taxonomy'].value_counts().head(8).index
        host_stack_hues = ['Source', 'Location', 'Tissue', 'CenterName', 'CollectionDate', 'LibrarySource']
        
        for virus in top_viruses:
            print(f"  -> 深度裂变绘制: {virus}")
            v_dir = os.path.join(profile_dir, f"Profile_{safe_filename(virus)}")
            os.makedirs(v_dir, exist_ok=True)
            
            df_v = self.df_merged[self.df_merged['Taxonomy'] == virus].copy()
            
            # 【A】常规透视
            for feature in self.meta_features:
                if feature == 'ScientificName' or feature not in df_v.columns: continue
                
                top_features = df_v[feature].value_counts().head(8).index
                df_v[f'{feature}_clean'] = df_v[feature].apply(lambda x: x if x in top_features else 'Other')
                
                ct = pd.crosstab(df_v[f'{feature}_clean'], df_v['ScientificName'])
                
                if feature == 'CollectionDate':
                    others_row = ct[ct.index == 'Other']
                    times_row = ct[ct.index != 'Other'].sort_index(ascending=True)
                    ct = pd.concat([times_row, others_row])
                else:
                    ct['Total'] = ct.sum(axis=1)
                    ct = ct.sort_values('Total', ascending=False).drop(columns='Total')
                    
                if ct.empty: continue
                
                fig, ax = plt.subplots(figsize=(10, 6))
                ct.plot(kind='bar', stacked=True, ax=ax, colormap='Set1', edgecolor='black', linewidth=0.5)
                
                plt.title(f"'{virus}' Infection across {feature}", fontsize=14, pad=15)
                plt.xlabel(feature, fontsize=12)
                plt.ylabel('Number of Detected Infections', fontsize=12)
                plt.xticks(rotation=45, ha='right')
                plt.legend(title='Host (ScientificName)', bbox_to_anchor=(1.05, 1), loc='upper left')
                
                self._save_plot_and_table(ct, v_dir, f"vs_{feature}", fig, ax, kind='bar')

            # 【B】宿主专属透视
            if 'ScientificName' in df_v.columns:
                top_sci = df_v['ScientificName'].value_counts().head(8).index
                df_v['Sci_clean'] = df_v['ScientificName'].apply(lambda x: x if x in top_sci else 'Other')
                
                for hue_col in host_stack_hues:
                    if hue_col not in df_v.columns: continue
                    
                    top_hues = df_v[hue_col].value_counts().head(8).index
                    df_v[f'{hue_col}_clean'] = df_v[hue_col].apply(lambda x: x if x in top_hues else 'Other')
                    
                    ct = pd.crosstab(df_v['Sci_clean'], df_v[f'{hue_col}_clean'])
                    ct['Total'] = ct.sum(axis=1)
                    ct = ct.sort_values('Total', ascending=False).drop(columns='Total')
                    if ct.empty: continue
                    
                    fig, ax = plt.subplots(figsize=(10, 6))
                    ct.plot(kind='bar', stacked=True, ax=ax, colormap='Set2', edgecolor='black', linewidth=0.5)
                    
                    plt.title(f"'{virus}' Hosts Stacked by {hue_col}", fontsize=14, pad=15)
                    plt.xlabel('Host (ScientificName)', fontsize=12)
                    plt.ylabel('Number of Detected Infections', fontsize=12)
                    plt.xticks(rotation=45, ha='right')
                    plt.legend(title=hue_col, bbox_to_anchor=(1.05, 1), loc='upper left')
                    
                    filename = f"ScientificName_stacked_by_{hue_col}"
                    self._save_plot_and_table(ct, v_dir, filename, fig, ax, kind='bar')

    def run_all(self):
        self.load_and_merge()
        self.plot_overall_summary()
        self.plot_global_features_vs_virus()
        self.plot_infection_complexity_breakdown()
        self.plot_virus_specific_profiles()
        
        out_table = os.path.join(self.out_dir, "Viral_Infection_with_Metadata_Summary.tsv")
        self.df_merged.to_csv(out_table, sep='\t', index=False, encoding='utf-8-sig')
        print(f"\n🎉 大一统完美收工！所有模块可视化图表已保存至: {os.path.abspath(self.out_dir)}")

def main():
    parser = argparse.ArgumentParser(description="🦠 病毒组学 vs 11维宿主元数据 联合分析与绘图引擎")
    parser.add_argument("-v", "--virus", required=True, help="病毒鉴定汇总表 (.csv/.tsv)")
    parser.add_argument("-m", "--meta", default=None, help="元数据表路径 (缺省自动查找/生成)")
    parser.add_argument("-o", "--outdir", default="./virus_multi_analysis", help="图表输出目录")
    parser.add_argument("--sra_list", default=None, help="SRA编号列表文件 (触发自动生成元数据)")
    parser.add_argument("--meta_outdir", default=None, help="元数据输出目录 (缺省: virus_multi_analysis/../sra_results)")
    parser.add_argument("--deepseek_api", default=None, help="DeepSeek API Key")
    parser.add_argument("--ncbi_api", default=None, help="NCBI API Key")
    parser.add_argument("--email", default=None, help="NCBI 联系邮箱")
    parser.add_argument("--meta_threads", type=int, default=4, help="元数据解析线程数")
    args = parser.parse_args()

    # ── 解析最终使用的元数据文件 ──
    meta_file = args.meta

    if meta_file and os.path.exists(meta_file):
        pass
    elif meta_file is None or not os.path.exists(meta_file):
        script_dir = Path(__file__).resolve().parent
        gsa_script = script_dir.parent / "metadata" / "gsa_sra.info.py"
        default_outdir = args.meta_outdir or os.path.join(os.path.dirname(args.outdir), "sra_results")
        default_meta = os.path.join(default_outdir, "Global_Unified_Metadata_Core13.tsv")

        if os.path.exists(default_meta):
            meta_file = default_meta
        elif args.sra_list and os.path.exists(args.sra_list) and gsa_script.exists():
            print(f"⏳ 未找到元数据文件，自动运行 metadata/gsa_sra.info.py 生成...")
            cmd = [
                sys.executable, str(gsa_script),
                "-i", args.sra_list,
                "-o", default_outdir,
                "-m", "both",
                "-t", str(args.meta_threads),
            ]
            if args.deepseek_api:
                cmd += ["--deepseek-api", args.deepseek_api]
            if args.ncbi_api:
                cmd += ["--ncbi-api", args.ncbi_api]
            if args.email:
                cmd += ["--email", args.email]
            result = subprocess.run(cmd)
            if result.returncode != 0:
                print(f"❌ 元数据生成失败 (exit={result.returncode})")
                sys.exit(1)
            meta_file = default_meta
        else:
            print("❌ 未找到元数据文件。请提供 --meta 或 --sra_list 参数。")
            print(f"   期望路径: {default_meta}")
            print(f"   或手动运行: python {gsa_script} -i <sra.list> -o {default_outdir}")
            sys.exit(1)

    if not meta_file or not os.path.exists(meta_file):
        print(f"❌ 元数据文件不存在: {meta_file}")
        sys.exit(1)

    print(f"📋 元数据: {meta_file}")
    analyzer = VirusMetadataAnalyzer(args.virus, meta_file, args.outdir)
    analyzer.run_all()

if __name__ == "__main__":
    main()
