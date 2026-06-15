#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
🧬 GSA & SRA Global Species-Targeted Retrieval Engine (AI 极致提纯终极版)
核心修复：
1. GSA Excel 解析现在会正确抓取并拼接 'Age unit'。
2. 重写了 AI 提示词，强制剔除 Tissue 中的环境/状态修饰语 (如 under stress, young)。
3. 强化了 Location 的三级自动推理补全 (如 China:Yinchuan -> China, Ningxia, Yinchuan_AI)。
4. 清洗了各种形式的 missing/not collected 等无效信息。
"""

import os
import re
import time
import json
import random
import argparse
import requests
import pandas as pd
from io import StringIO
from tqdm import tqdm
from urllib.parse import quote
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
from openai import OpenAI
import warnings

warnings.filterwarnings("ignore")

HEADERS = {
    "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36",
    "Accept": "text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8",
    "Connection": "keep-alive"
}

def get_retry_session(retries=3, backoff_factor=1.5):
    session = requests.Session()
    retry_strategy = Retry(
        total=retries, backoff_factor=backoff_factor, 
        status_forcelist=[403, 429, 500, 502, 503, 504]
    )
    adapter = HTTPAdapter(max_retries=retry_strategy, pool_connections=20, pool_maxsize=20)
    session.mount("https://", adapter)
    session.mount("http://", adapter)
    session.headers.update(HEADERS)
    return session

# ==========================================
# 🟢 SRA 检索引擎
# ==========================================
class SRAEngine:
    def __init__(self, query, source, out_dir, detailed=False, ncbi_api=None):
        self.query = query
        self.source = source
        self.out_dir = out_dir
        self.detailed = detailed
        self.ncbi_api = ncbi_api
        self.session = get_retry_session()

    def fetch_runinfo(self):
        print("\n" + "="*65)
        mode_str = "详细模式 (XML解析)" if self.detailed else "极速模式 (基础信息)"
        print(f"🟢 启动 SRA 检索引擎 [{mode_str}]")
        print("="*65)
        
        term = f'"{self.query}"[Organism]'
        if self.source:
            if self.source.upper() == 'TRANSCRIPTOMIC':
                term += ' AND "biomol rna"[Properties]'
            elif self.source.upper() == 'GENOMIC':
                term += ' AND "biomol dna"[Properties]'
            else:
                term += f' AND "{self.source}"[Properties]'
                
        print(f"🧩 构建 SRA 检索逻辑: {term}")
        
        try:
            esearch_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
            params = {"db": "sra", "term": term, "usehistory": "y", "retmode": "json"}
            if self.ncbi_api: params["api_key"] = self.ncbi_api
            res = self.session.get(esearch_url, params=params, timeout=30).json()
            count = int(res.get('esearchresult', {}).get('count', 0))
            if count == 0:
                print("⚠️ 未在 SRA 找到相关数据。")
                return pd.DataFrame()
                
            print(f"🎯 锁定 {count} 个 Run。正在下载 RunInfo 基础表...")
            webenv = res['esearchresult']['webenv']
            query_key = res['esearchresult']['querykey']
            
            efetch_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
            params_csv = {"db": "sra", "query_key": query_key, "WebEnv": webenv, "rettype": "runinfo", "retmode": "text"}
            if self.ncbi_api: params_csv["api_key"] = self.ncbi_api
            res_csv = self.session.get(efetch_url, params=params_csv, timeout=60)
            df_sra = pd.read_csv(StringIO(res_csv.text))
            
            if self.detailed:
                print("🧬 [详细模式] 正在批量拉取 SRA XML 获取深层生物学特征...")
                params_xml = {"db": "sra", "query_key": query_key, "WebEnv": webenv, "rettype": "xml", "retmode": "text"}
                if self.ncbi_api: params_xml["api_key"] = self.ncbi_api
                res_xml = self.session.get(efetch_url, params=params_xml, timeout=120)
                xml_text = res_xml.text
                
                meta_dict = {}
                for pkg_match in re.finditer(r'<EXPERIMENT_PACKAGE>([\s\S]*?)</EXPERIMENT_PACKAGE>', xml_text):
                    pkg = pkg_match.group(1)
                    runs = re.findall(r'<RUN[^>]*accession="([E|S|D]RR\d+)"', pkg)
                    
                    attrs = re.findall(r'<SAMPLE_ATTRIBUTE>\s*<TAG>([\s\S]*?)</TAG>\s*<VALUE>([\s\S]*?)</VALUE>', pkg, re.I)
                    attr_dict = {t.strip().lower(): v.strip() for t, v in attrs}
                    
                    tissue = attr_dict.get('tissue', attr_dict.get('cell type', attr_dict.get('tissue type', pd.NA)))
                    loc = attr_dict.get('geo_loc_name', attr_dict.get('country', attr_dict.get('geographic location', pd.NA)))
                    
                    age_parts = []
                    for k in ['age', 'dev_stage', 'development stage', 'growth stage']:
                        if k in attr_dict and attr_dict[k]:
                            age_parts.append(attr_dict[k])
                    stage = " | ".join(age_parts) if age_parts else pd.NA
                    
                    for r in runs:
                        meta_dict[r] = {'Tissue': tissue, 'Age_GrowthStage': stage, 'Location': loc}
                
                df_sra['Tissue'] = df_sra['Run'].map(lambda x: meta_dict.get(x, {}).get('Tissue', pd.NA))
                df_sra['Age_GrowthStage'] = df_sra['Run'].map(lambda x: meta_dict.get(x, {}).get('Age_GrowthStage', pd.NA))
                df_sra['Location'] = df_sra['Run'].map(lambda x: meta_dict.get(x, {}).get('Location', pd.NA))

            return df_sra
            
        except Exception as e:
            print(f"❌ SRA 获取失败: {e}")
            return pd.DataFrame()

# ==========================================
# 🔵 GSA 检索引擎
# ==========================================
class GSAEngine:
    def __init__(self, query, source, out_dir, detailed=False):
        self.query = str(query).strip()
        self.source_filter = str(source).strip() if source else None
        self.center_filter = "NGDC" 
        self.detailed = detailed
        self.base_dir = os.path.join(out_dir, "GSA_Results")
        self.d_web = os.path.join(self.base_dir, "0_web_cache")
        self.d_xls = os.path.join(self.base_dir, "1_xls_cache")
        os.makedirs(self.d_web, exist_ok=True)
        if self.detailed:
            os.makedirs(self.d_xls, exist_ok=True)
        self.session = get_retry_session()

    def get_accession_list(self):
        url = "https://ngdc.cncb.ac.cn/gsa/search/getAccessionList"
        term = f'&quot;{self.query}&quot;[organism]'
        if self.source_filter: term = f'({term} AND &quot;{self.source_filter}&quot;[source])'
        term = f'({term} AND &quot;{self.center_filter}&quot;[center])'

        clean_term = term.replace('&quot;', '"')
        print(f"🧩 构建 GSA 检索逻辑: {clean_term}")
        
        payload = f"searchField=&searchTerm={quote(term)}&totalDatas=99999"
        headers = self.session.headers.copy()
        headers["Content-Type"] = "application/x-www-form-urlencoded"

        try:
            res = self.session.post(url, data=payload, headers=headers, timeout=45)
            return sorted(list(set(re.findall(r'(CR[RPAX]\d{6,})', res.text))))
        except Exception as e:
            return []

    def parse_all_from_html(self, html, acc):
        feat = {
            "ReleaseDate": pd.NA, "Organization": pd.NA, "CRA": pd.NA, "CRX": pd.NA, 
            "PRJ": pd.NA, "SAMC": pd.NA, "Platform": pd.NA, "TaxID": pd.NA, "ScientificName": pd.NA,
            "LibraryStrategy": "Not_Provided", "LibrarySource": "Not_Provided", "RunRecords": []
        }
        if not html: return feat

        def extract_text(pat):
            m = re.search(pat, html, re.I)
            return re.sub(r'<[^>]+>', '', m.group(1)).strip() if m else pd.NA

        feat["ReleaseDate"] = extract_text(r'<th[^>]*>(?:发布日期|Release date)</th>\s*<td[^>]*>([\s\S]*?)</td>')
        feat["Organization"] = extract_text(r'<th[^>]*>(?:所属单位|Organization)</th>\s*<td[^>]*>([\s\S]*?)</td>')
        feat["ScientificName"] = extract_text(r'wwwtax\.cgi\?id=\d+"[^>]*>([\s\S]*?)</a>')
        feat["CRX"] = extract_text(r'<th[^>]*>(?:实验编号|Accession)</th>\s*<td[^>]*>.*?(CRX\d+)[\s\S]*?</td>')
        feat["PRJ"] = extract_text(r'<th[^>]*>(?:项目编号|BioProject)</th>\s*<td[^>]*>.*?(PRJ[A-Z]*\d+)[\s\S]*?</td>')
        feat["SAMC"] = extract_text(r'<th[^>]*>(?:样本编号|BioSample)</th>\s*<td[^>]*>.*?(SAM[A-Z]*\d+)[\s\S]*?</td>')
        feat["Platform"] = extract_text(r'<th[^>]*>(?:测序平台|Platform)</th>\s*<td[^>]*>([\s\S]*?)</td>')
        
        cra_m = re.search(r'(CRA\d{6,})', html)
        if cra_m: feat["CRA"] = cra_m.group(1)

        lib_table_m = re.search(r'<th[^>]*>(?:建库信息|Library)</th>\s*<td[^>]*>[\s\S]*?<table[^>]*>([\s\S]*?)</table>', html, re.I)
        if lib_table_m:
            tds = re.findall(r'<td[^>]*>([\s\S]*?)</td>', lib_table_m.group(1), re.I)
            if len(tds) >= 4:
                feat["LibraryStrategy"] = re.sub(r'<[^>]+>', '', tds[2]).strip()
                feat["LibrarySource"] = re.sub(r'<[^>]+>', '', tds[3]).strip()

        file_block_m = re.findall(r'<a href="[^"]*browse/(CRA\d+)/(CRR\d+)"[^>]*>[\s\S]*?</a>\s*</td>\s*<td[^>]*>([\s\S]*?)</td>\s*<td[^>]*>([\s\S]*?)</td>', html, re.I)
        if file_block_m:
            feat["CRA"] = file_block_m[0][0] 
            for match in file_block_m:
                if re.findall(r'\.(?:fq|fastq|bam|sra)', match[2], re.I): feat["RunRecords"].append({"Run": match[1]})
        else:
            file_m = re.search(rf'<a[^>]*>{acc}</a>[\s\S]*?</td>\s*<td[^>]*>([\s\S]*?)</td>\s*<td[^>]*>([\s\S]*?)</td>', html, re.I)
            if file_m and re.findall(rf'\.(?:fq|fastq|bam|sra)', file_m.group(1), re.I): feat["RunRecords"].append({"Run": acc})

        return feat

    def fetch_excel_attributes(self, cra_list):
        meta_dict = {}
        print("\n📊 [详细模式] 正在向 CNCB 提交 Excel 级联解析请求...")
        for cra in tqdm(cra_list, desc="解析 Excel 附件"):
            if pd.isna(cra): continue
            xls_path = os.path.join(self.d_xls, f"{cra}.xlsx")
            
            if not os.path.exists(xls_path):
                try:
                    rx = self.session.post("https://ngdc.cncb.ac.cn/gsa/file/exportExcelFile", data={"type": 3, "dlAcession": cra}, timeout=30)
                    if len(rx.content) > 1000:
                        with open(xls_path, "wb") as f: f.write(rx.content)
                except: pass
            
            if os.path.exists(xls_path):
                try:
                    xls = pd.ExcelFile(xls_path, engine='openpyxl')
                    if 'Sample' in xls.sheet_names:
                        df_samp = pd.read_excel(xls, sheet_name='Sample')
                        cols = df_samp.columns.str.lower().str.strip()
                        df_samp.columns = cols
                        
                        acc_col = 'accession' if 'accession' in cols else 'biosample accession'
                        tissue_col = 'tissue' if 'tissue' in cols else 'organism part'
                        
                        for _, row in df_samp.iterrows():
                            sam = row.get(acc_col)
                            if pd.isna(sam): continue
                            
                            # 【核心修复】智能拼接 age 和 age unit
                            age_val = row.get('age')
                            age_unit = row.get('age unit')
                            stage_val = row.get('dev stage')
                            
                            age_str = ""
                            if pd.notna(age_val):
                                age_str = str(age_val).strip()
                                if pd.notna(age_unit):
                                    age_str += f" {str(age_unit).strip()}"
                                    
                            stage_parts = []
                            if age_str: stage_parts.append(age_str)
                            if pd.notna(stage_val): stage_parts.append(str(stage_val).strip())
                            
                            meta_dict[sam] = {
                                'Tissue': row.get(tissue_col, pd.NA),
                                'Age_GrowthStage': " | ".join(stage_parts) if stage_parts else pd.NA,
                                'Location': row.get('geographic location', pd.NA)
                            }
                except: pass
        return meta_dict

    def fetch_gsa(self):
        print("\n" + "="*65)
        mode_str = "详细模式 (含Excel解析)" if self.detailed else "极速模式 (基础信息)"
        print(f"🔵 启动 GSA 检索引擎 [{mode_str}]")
        print("="*65)
        
        acc_list = self.get_accession_list()
        if not acc_list: return pd.DataFrame()

        all_records = []
        unique_cras = set()
        
        for acc in tqdm(acc_list, desc="🔵 解析 HTML"):
            web_cache_f = os.path.join(self.d_web, f"{acc}.html")
            try:
                if not os.path.exists(web_cache_f):
                    time.sleep(random.uniform(1.0, 2.0))
                    res_web = self.session.get(f"https://ngdc.cncb.ac.cn/gsa/search?searchTerm={acc}", timeout=30)
                    with open(web_cache_f, 'w', encoding='utf-8') as f: f.write(res_web.text)
                with open(web_cache_f, 'r', encoding='utf-8') as f: html = f.read()

                wf = self.parse_all_from_html(html, acc)
                if self.query.lower() not in str(wf.get("ScientificName", "")).lower(): continue
                
                cra_id = wf.get("CRA", "UNKNOWN_CRA")
                unique_cras.add(cra_id)
                
                for r in wf.get("RunRecords", [{"Run": acc}]):
                    all_records.append({
                        "Run": r["Run"], "ReleaseDate": wf.get("ReleaseDate", pd.NA),
                        "LibraryStrategy": wf.get("LibraryStrategy"), "LibrarySource": wf.get("LibrarySource"),
                        "BioProject": wf.get("PRJ", pd.NA), "BioSample": wf.get("SAMC", pd.NA), 
                        "Platform": wf.get("Platform", pd.NA), "Organization": wf.get("Organization", pd.NA)
                    })
            except Exception as e: pass

        df_gsa = pd.DataFrame(all_records)
        
        if self.detailed and not df_gsa.empty and unique_cras:
            excel_meta = self.fetch_excel_attributes(unique_cras)
            df_gsa['Tissue'] = df_gsa['BioSample'].map(lambda x: excel_meta.get(x, {}).get('Tissue', pd.NA))
            df_gsa['Age_GrowthStage'] = df_gsa['BioSample'].map(lambda x: excel_meta.get(x, {}).get('Age_GrowthStage', pd.NA))
            df_gsa['Location'] = df_gsa['BioSample'].map(lambda x: excel_meta.get(x, {}).get('Location', pd.NA))

        return df_gsa

# ==========================================
# 🤖 AI 智能清洗与规范化引擎 (全新提纯规则)
# ==========================================
def run_ai_sanitizer(df, api_key, api_base, model):
    if df.empty or not api_key: return df
    print("\n🤖 启动 AI 洗髓引擎 (深度规范 Tissue, Age_GrowthStage, Location)...")
    
    client = OpenAI(api_key=api_key, base_url=api_base)
    cols = ['Location', 'Tissue', 'Age_GrowthStage']
    for c in cols: 
        if c not in df.columns: df[c] = pd.NA

    prompt = """你是一个极其严苛的生命科学数据清理程序。请严格按以下规则处理输入的JSON，返回清洗后的JSON：
1. **全局去杂**：遇到 'not applicable', 'not collected', 'missing', 'N/A', 'nan', 或只包含无意义编号(如 'tissue9')，一律清空替换为 'Not_Provided'。
2. **Tissue (核心剥离)**：强制剔除所有的状态、环境、时间、物种等修饰语（例如：将 "flower under drought stress in the stage 2", "leaf of Wolfberry", "young leaves" 统统剥离）。你【只能】输出最纯粹的单数英文核心组织词，例如：leaf, root, stem, flower, fruit, seed, anther, stamen, pistil, whole plant 等。绝对不要保留多余信息！
3. **Age_GrowthStage (凝练)**：剥离多余符号和无效重复，合并年龄和时期信息。例如将 "3 | Young Fruit" 规范为 "3 years, young fruit stage"；将 "not applicable | missing" 清除为 "Not_Provided"。
4. **Location (三级正交化)**：必须凭借你的地理知识，强制补全为【国家, 省/州, 市/县_AI】的标准格式（如 "China:Yinchuan" 或 "China:yinchuan" 必须补全为 "China, Ningxia, Yinchuan_AI"；"China:Qaidam Basin" 补全为 "China, Qinghai, Qaidam Basin_AI"）。实在缺失的层级用 Unknown 补全。
仅返回合法 JSON，勿带 ```json 标记。绝对禁止凭空捏造数据！"""

    for idx, row in tqdm(df.iterrows(), total=len(df), desc="AI Sanitization"):
        if pd.isna(row.get('Location')) and pd.isna(row.get('Tissue')) and pd.isna(row.get('Age_GrowthStage')): continue
        target = {k: str(row.get(k, "Not_Provided")) for k in cols}
        
        try:
            res = client.chat.completions.create(
                model=model,
                messages=[{"role": "system", "content": prompt}, {"role": "user", "content": json.dumps(target, ensure_ascii=False)}],
                temperature=0.1
            )
            clean_str = re.sub(r"^```(?:json)?\s*|\s*```$", "", res.choices[0].message.content.strip(), flags=re.IGNORECASE)
            clean_data = json.loads(clean_str)
            
            for k in cols:
                if k in clean_data:
                    val = str(clean_data[k]).strip()
                    if val.lower() in ["", "none", "nan", "not_provided", "unknown"]:
                        df.at[idx, k] = pd.NA
                    else:
                        df.at[idx, k] = val
        except: pass
    return df

# ==========================================
# 🌍 数据大一统合并
# ==========================================
def merge_results(df_sra, df_gsa, out_dir, detailed, api_key, api_base, model):
    if not df_sra.empty:
        df_sra = df_sra.rename(columns={'CenterName': 'Organization_CenterName'})
        df_sra['Database'] = 'SRA'
    else: df_sra = pd.DataFrame()

    if not df_gsa.empty:
        df_gsa = df_gsa.rename(columns={'Organization': 'Organization_CenterName'})
        df_gsa['Database'] = 'GSA'
    else: df_gsa = pd.DataFrame()

    df_merged = pd.concat([df_sra, df_gsa], ignore_index=True)
    if df_merged.empty: return

    if detailed and api_key:
        df_merged = run_ai_sanitizer(df_merged, api_key, api_base, model)

    target_cols = ['Database', 'Run', 'BioProject', 'BioSample']
    if detailed:
        target_cols.extend(['Tissue', 'Age_GrowthStage', 'Location'])
    target_cols.extend(['LibraryStrategy', 'LibrarySource', 'Platform', 'Organization_CenterName', 'ReleaseDate'])
    
    final_cols = [c for c in target_cols if c in df_merged.columns]
    
    df_final = df_merged[final_cols].drop_duplicates(subset=['Run'])
    
    out_file = os.path.join(out_dir, "SRA_GSA_Merged_Final.csv")
    df_final.to_csv(out_file, index=False, encoding='utf-8-sig')
    
    print(f"\n🎉 大一统合并成功！共汇总 {len(df_final)} 条 Run 记录。")
    print(f"📁 最终输出路径: {out_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="SRA + GSA 双引擎物种检索与提纯工具 (双模式版)")
    parser.add_argument("-q", "--query", required=True, help="输入物种拉丁名 (如 'Lycium barbarum')")
    parser.add_argument("-s", "--source", help="限制测序类型 (如 'TRANSCRIPTOMIC')")
    parser.add_argument("-o", "--outdir", default="./Global_Species_Results", help="输出根目录")
    
    parser.add_argument("--detailed", action="store_true", help="开启详细模式，抓取并解析 Tissue/Stage/Location 等深层特征")
    parser.add_argument("--ncbi-api", help="NCBI E-utilities API Key (提升速率)")
    parser.add_argument("--deepseek-api", help="DeepSeek API Key (AI 智能清洗)")
    parser.add_argument("--deepseek-model", default="deepseek-chat", help="DeepSeek 模型名称")

    args = parser.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    sra_engine = SRAEngine(args.query, args.source, args.outdir, detailed=args.detailed, ncbi_api=args.ncbi_api)
    df_sra = sra_engine.fetch_runinfo()

    gsa_engine = GSAEngine(args.query, args.source, args.outdir, detailed=args.detailed)
    df_gsa = gsa_engine.fetch_gsa()

    merge_results(df_sra, df_gsa, args.outdir, args.detailed, args.deepseek_api, "https://api.deepseek.com", args.deepseek_model)
