#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
stage_cards.py — MMPV 阶段卡片 (四层递进架构 · L1 知识层)
==========================================================

把「每条管线的每个 stage 是什么、吃什么、吐什么、会怎么失败」收敛为
结构化卡片 —— 一份真相源同时喂三张嘴:

  * 文档 (--help / wiki 引用本模块核对, 防文档漂移);
  * AI (run digest / advisor / copilot 的提示词上下文);
  * 人 (mmpv_ai.py cards 子命令直接查阅)。

内容来源: 各编排器源码 STAGE_ORDER / docstring 与 doc/*.md (2026-10 整理)。
阶段名与编排器 choices 逐字对齐; 编排器改动 stage 时请同步本表
(结构性校验见 tests/test_stage_cards.py)。

纯标准库, 所有模块可安全 import。
"""

from __future__ import annotations

from typing import Dict, List, Optional

# 卡片字段:
#   name     中文短名
#   desc     一句话职责
#   inputs / outputs  关键输入输出 (语义级, 不锁死绝对路径)
#   risks    失败模式 / 踩坑要点 (供 AI 解读与排障提示)
#   requires 前置 stage (同管线内)
STAGE_CARDS: Dict[str, Dict] = {
    "public": {
        "num": "①", "name": "Public Data Retrieval",
        "script": "public_metadata_pipeline/public_data_pipeline.py",
        "desc": "SRA/GSA 检索 → 元数据清洗 → 批量下载 → 转换 → 可视化",
        "stages": {
            "search":   {"name": "检索", "desc": "NCBI SRA + CNCB GSA 双引擎物种检索",
                         "inputs": ["物种名/taxid"], "outputs": ["SRA_GSA_Merged_Final.csv"],
                         "risks": ["NCBI 限流 (建议配 NCBI_API_KEY)"], "requires": []},
            "info":     {"name": "元数据清洗", "desc": "XML 解析 + GSA 爬虫 + AI 元数据清洗 → Core14",
                         "inputs": ["检索结果"], "outputs": ["Global_Unified_Metadata_Core14.csv"],
                         "risks": ["AI 只回填空白单元格, 已有真值不覆盖",
                                   "Location 强制三级格式 (国家, 省/州, 市/县_AI)"],
                         "requires": ["search"]},
            "prescreen": {"name": "预筛", "desc": "下载前按元数据规则预筛 (--prescreen-filter)",
                          "inputs": ["Core14"], "outputs": ["prescreen 名单"],
                          "risks": [], "requires": ["info"]},
            "down":     {"name": "下载", "desc": "aria2c/wget/prefetch 双协议批量下载",
                         "inputs": ["名单"], "outputs": ["原始 FASTQ/SRA (03_RawData)"],
                         "risks": ["大批量断点续传依赖 --resume"], "requires": ["info"]},
            "convert":  {"name": "转换", "desc": "sracha/fasterq-dump 批量 SRA → FASTQ.GZ (自动后端 + 回退)",
                         "inputs": [".sra"], "outputs": ["FASTQ.GZ"],
                         "risks": ["磁盘配额", "sracha 严格校验失败会降级 --no-strict 再回退 sra-tools"],
                         "requires": ["down"]},
            "plot":     {"name": "可视化", "desc": "时间/组织/地区/机构 6 面板 SCI 级图",
                         "inputs": ["Core14"], "outputs": ["04_Plots"], "risks": [], "requires": ["info"]},
            "hostref":  {"name": "宿主参考", "desc": "宿主基因组多通道下载 (datasets/ngd/FTP/gget; 或 --genome-fasta)",
                         "inputs": ["物种/taxid 或 FASTA"], "outputs": ["05_HostRef/genome"],
                         "risks": ["非模式植物 NCBI 缺种 → gget 通道"], "requires": []},
            "hostdb":   {"name": "宿主建库", "desc": "Kraken2/Bowtie2/HISAT2/Minimap2 四索引构建",
                         "inputs": ["宿主基因组+GFF3"], "outputs": ["05_HostRef/hostdb"],
                         "risks": ["建库耗时与磁盘大"], "requires": ["hostref"]},
            "report":   {"name": "报告", "desc": "交互式 HTML + sample_handoff.csv 交接清单",
                         "inputs": ["全阶段产物"], "outputs": ["HTML", "sample_handoff.csv"],
                         "risks": [], "requires": []},
        },
    },
    "preprocessing": {
        "num": "②", "name": "Data Preprocessing",
        "script": "data_preprocessing_pipeline/data_preprocessing.py",
        "desc": "质控 (Fastp→Seqkit→Clumpify) → 去宿主去 rRNA (Kraken2→比对→Ribodetector)",
        "stages": {
            "clean":   {"name": "质控", "desc": "Fastp QC + Seqkit 转换 + Clumpify 聚类重排",
                        "inputs": ["原始 FASTQ"], "outputs": ["01_CleanData (standard)"],
                        "risks": ["双端配对完整性靠脚本校验"], "requires": []},
            "deplete": {"name": "去宿主", "desc": "分类初筛 (Kraken2 / Kun-peng) → 比对 (bowtie2/hisat2/minimap2) → Ribodetector 去 rRNA",
                        "inputs": ["clean reads", "宿主四索引"], "outputs": ["02_HostDepleted"],
                        "risks": ["分类库必须与比对索引同物种口径",
                                  "keep_taxids 默认救回病毒 reads (10239)"], "requires": ["clean"]},
            "bbnorm":  {"name": "归一化", "desc": "覆盖度归一化 (可选, 默认不跑; 共组装前深度差异大才需要)",
                        "inputs": ["去宿主 reads"], "outputs": ["03_BBnorm"],
                        "risks": ["会改变丰度分布, 定量管线慎用"], "requires": ["deplete"]},
            "report":  {"name": "报告", "desc": "preprocessing_summary.tsv + assembly_ready.list + HTML",
                        "inputs": ["clean/deplete 产物"], "outputs": ["汇总 TSV", "交接清单", "HTML"],
                        "risks": [], "requires": ["clean"]},
        },
    },
    "discovery": {
        "num": "③", "name": "Virus Discovery (de novo)",
        "script": "virome_discovery_pipeline/virome_pipeline.py",
        "desc": "组装 → 10 工具鉴定 → 分层过滤 → 延伸 → 聚类 → 分类 → 宿主 → 评估 → 拯救 → 报告 (15 阶段)",
        "stages": {
            "clean":           {"name": "质控", "desc": "复用 ② 的 clean-data.py",
                                "inputs": ["原始 FASTQ"], "outputs": ["CleanData"],
                                "risks": ["与 ② 输出目录同键"], "requires": []},
            "deplete":         {"name": "去宿主", "desc": "复用 ② 的 host_depletion.py",
                                "inputs": ["clean reads", "宿主库"], "outputs": ["HostDepleted"],
                                "risks": [], "requires": ["clean"]},
            "assembly":        {"name": "组装", "desc": "MEGAHIT / rnaviralSPAdes / Penguin 组装",
                                "inputs": ["去宿主 reads"], "outputs": ["各工具 contig"],
                                "risks": ["RNA 病毒优先 rnaviralSPAdes"], "requires": ["deplete"]},
            "identification":  {"name": "鉴定", "desc": "Genomad+Diamond+RdRpCatch+ViraLM+VirBot+VirSorter2+ViralVerify+VirHunter+Metabuli+Viroid 10 工具并行",
                                "inputs": ["contig"], "outputs": ["逐工具判定表"],
                                "risks": ["10 库全路径必须在 config 配齐"], "requires": ["assembly"]},
            "filter":          {"name": "过滤", "desc": "UniProt-strict 高置信过滤 (CDD 层已移除); 可选 graph 子步骤二次裁决",
                                "inputs": ["逐工具判定表"], "outputs": ["高置信病毒 contig"],
                                "risks": ["graph.enabled 默认关闭"], "requires": ["identification"]},
            "cobra":           {"name": "延伸", "desc": "BWA-MEM2 → CoverM → COBRA-Meta 单样本延伸",
                                "inputs": ["候选 contig", "reads"], "outputs": ["延伸 contig"],
                                "risks": ["输入=filter 主链镜像"], "requires": ["filter"]},
            "merge":           {"name": "共组装", "desc": "多样本 Flye 共组装 (co-assembly 模式)",
                                "inputs": ["延伸 contig"], "outputs": ["共组装 contig"],
                                "risks": ["单样本项目自动跳过"], "requires": ["cobra"]},
            "cluster":         {"name": "聚类", "desc": "CD-HIT 参考引导预聚类 + vclust Leiden → vOTU",
                                "inputs": ["合并 contig"], "outputs": ["vOTU 代表序列"],
                                "risks": [], "requires": ["merge"]},
            "taxonomy":        {"name": "分类", "desc": "8 工具分类加权投票 → 8 级 taxonomy",
                                "inputs": ["vOTU"], "outputs": ["taxonomy.tsv"],
                                "risks": ["AI 不得改写投票结论, 只可解读"], "requires": ["cluster"]},
            "host":            {"name": "宿主预测", "desc": "ICTV > RNAVirHost > PhaBOX2 三级 + 植物专属后过滤",
                                "inputs": ["vOTU"], "outputs": ["宿主预测表"],
                                "risks": [], "requires": ["cluster"]},
            "checkv":          {"name": "完整性", "desc": "CheckV 完整性评估",
                                "inputs": ["vOTU"], "outputs": ["completeness 表"],
                                "risks": ["阈值 checkv_threshold 默认 90"], "requires": ["cluster"]},
            "rescue":          {"name": "拯救", "desc": "四支路级联: A CheckV≥90 直通 → B VSI 延伸 → C BLASTN 参考重建 → D genus-length 兜底",
                                "inputs": ["checkv + reads + 参考库"], "outputs": ["HQ vOTU catalog"],
                                "risks": ["--force 重跑代价大"], "requires": ["checkv"]},
            "analysis":        {"name": "下游分析", "desc": "整合拯救证据 + 病毒组下游分析",
                                "inputs": ["rescue 产物"], "outputs": ["分析汇总"], "risks": [],
                                "requires": ["rescue"]},
            "analysis_verify": {"name": "复核", "desc": "HMM/CT3 证据扫描复核",
                                "inputs": ["analysis 产物"], "outputs": ["复核表"], "risks": [],
                                "requires": ["analysis"]},
            "report":          {"name": "报告", "desc": "TSV 汇总 + Sankey + 交互式 HTML",
                                "inputs": ["全阶段产物"], "outputs": ["13_Report HTML"], "risks": [],
                                "requires": []},
        },
    },
    "analysis": {
        "num": "④", "name": "Known-Virus Deep Analysis",
        "script": "virome_analysis_pipeline/auto_known_virus.py",
        "desc": "对已知参考病毒: 定量 → 变异 → 全长重建 → 相似性 → DVG → 报告",
        "stages": {
            "detect":     {"name": "定量", "desc": "Salmon/Kallisto/Bowtie2 快速定量 + Poisson 过滤",
                           "inputs": ["去宿主 reads", "病毒参考+ref_info"],
                           "outputs": ["定量表"], "risks": ["参考须配套 ref_info.tsv"],
                           "requires": []},
            "filter":     {"name": "过滤", "desc": "高置信度过滤",
                           "inputs": ["定量表"], "outputs": ["过滤后定量"], "risks": [],
                           "requires": ["detect"]},
            "variants":   {"name": "变异", "desc": "FreeBayes/iVar/LoFreq 变异检出 + 共识序列",
                           "inputs": ["reads", "参考"], "outputs": ["VCF + 共识 FASTA"],
                           "risks": ["ivar 需特定 BAM 排序"], "requires": ["detect"]},
            "post":       {"name": "变异后处理", "desc": "VCF 可视化 + SnpEff 注释 + MAF 瀑布 + SnpGenie dN/dS",
                           "inputs": ["VCF"], "outputs": ["注释与图表"], "risks": [],
                           "requires": ["variants"]},
            "full":       {"name": "全长构建", "desc": "OmniVirusAssembler 12 步全长基因组重建 + 尾段集成全长序列提取",
                           "inputs": ["靶向 reads", "参考"], "outputs": ["全长基因组", "06_Extraction FASTA"],
                           "risks": ["计算重; RNA 注意 N<5% 才 N-fill"], "requires": ["extract"]},
            "extract":    {"name": "靶向 reads 提取", "desc": "从 virus-bam 按病毒×样本提取 R1/R2/single (batch_extract_reads.py)",
                           "inputs": ["virus-bam"], "outputs": ["virus_reads/"], "risks": [],
                           "requires": ["variants"]},
            "similarity": {"name": "相似性", "desc": "全长相似性热图 + 层次聚类 (MAFFT)",
                           "inputs": ["全长 FASTA"], "outputs": ["热图/聚类"], "risks": [],
                           "requires": ["full"]},
            "dvg":        {"name": "DVG", "desc": "ViReMa DVG/重组检测 + Circos",
                           "inputs": ["reads", "参考"], "outputs": ["DVG 表 + Circos 图"],
                           "risks": ["ViReMa 对嵌合读长敏感"], "requires": ["detect"]},
            "report":     {"name": "报告", "desc": "交互式 HTML 综合报告 + AI 解读面板",
                           "inputs": ["全阶段产物"], "outputs": ["09_Report HTML"], "risks": [],
                           "requires": []},
        },
    },
    "phylo": {
        "num": "⑤", "name": "Phylogeny / Evolution / Population Genetics",
        "script": "virome_phylo_pipeline/phylo_pipeline.py",
        "desc": "18 stage / 8 大模块: 数据→比对→树→重组/选择→时间→地理→报告",
        "stages": {
            "prep":        {"name": "数据收集", "desc": "序列+日期+元数据汇总",
                            "inputs": ["④ 06_Extraction", "元数据"], "outputs": ["data/"], "risks": [],
                            "requires": []},
            "metadata":    {"name": "元数据治理", "desc": "时间+地理检查与矫正 + 三通道拆分",
                            "inputs": ["元数据"], "outputs": ["治理后元数据"],
                            "risks": ["坐标不臆造原则 (S1 系列)"], "requires": ["prep"]},
            "online":      {"name": "在线补充", "desc": "NCBI 参考序列补充 (SeqHarvester)",
                            "inputs": ["物种清单"], "outputs": ["参考序列"], "risks": ["需网络"],
                            "requires": ["prep"]},
            "clean":       {"name": "序列清洗", "desc": "长度/N/大小写/非法字符/元数据缺失过滤",
                            "inputs": ["序列"], "outputs": ["清洗序列"], "risks": [],
                            "requires": ["prep"]},
            "align":       {"name": "比对", "desc": "MAFFT 全长比对 + 替换饱和分析 (Iss)",
                            "inputs": ["序列"], "outputs": ["比对 + 饱和诊断"],
                            "risks": ["Iss 高 → 警惕随机化"], "requires": ["clean"]},
            "splitstree":  {"name": "网络", "desc": "SplitsTree 重组网络", "inputs": ["比对"],
                            "outputs": ["网络图"], "risks": [], "requires": ["align"]},
            "rdp5":        {"name": "重组", "desc": "RDP5 重组检测", "inputs": ["比对"],
                            "outputs": ["重组事件表"], "risks": ["Windows GUI 依赖注意"], "requires": ["align"]},
            "tree":        {"name": "进化树", "desc": "ML 树构建", "inputs": ["比对"], "outputs": ["newick"],
                            "risks": [], "requires": ["align"]},
            "popgen":      {"name": "群体遗传", "desc": "群体遗传学统计", "inputs": ["比对"],
                            "outputs": ["popgen 表"], "risks": [], "requires": ["align"]},
            "host":        {"name": "宿主关联", "desc": "宿主跳跃/关联分析", "inputs": ["树+宿主"],
                            "outputs": ["宿主关联表"], "risks": [], "requires": ["tree"]},
            "capheine":    {"name": "正选择", "desc": "HyPhy FEL/MEME/BUSTED/PRIME (Capheine)",
                            "inputs": ["比对"], "outputs": ["选择压力表"],
                            "risks": ["v3.0 起归属本管线"], "requires": ["prep"]},
            "clock":       {"name": "分子钟", "desc": "根-支离散分子钟 (rtt/temporal/genes 合并)",
                            "inputs": ["树+日期"], "outputs": ["进化速率"],
                            "risks": ["时间信号不足会失败, 看 date 拉散度"], "requires": ["tree"]},
            "beast":       {"name": "BEAST 定年", "desc": "BEAST 定年树", "inputs": ["比对+日期"],
                            "outputs": ["定年树"], "risks": ["单独部署, 见依赖文档 A12"],
                            "requires": ["tree"]},
            "gene_dating": {"name": "基因定年", "desc": "分段基因定年", "inputs": ["clock+capheine"],
                            "outputs": ["基因定年表"], "risks": [], "requires": ["clock", "capheine"]},
            "phylogeo":    {"name": "系统地理", "desc": "系统地理分析", "inputs": ["树+地理"],
                            "outputs": ["地理传播结构"], "risks": [], "requires": ["tree"]},
            "geo":         {"name": "地理解析", "desc": "geo_resolver 坐标/地名解析", "inputs": ["元数据"],
                            "outputs": ["解析坐标表"], "risks": ["坐标不臆造: 查不到留空"],
                            "requires": ["tree"]},
            "geo_paths":   {"name": "传播路径", "desc": "transmission_map 传播路径图",
                            "inputs": ["phylogeo merged/mcc.tree"], "outputs": ["传播地图"], "risks": [],
                            "requires": ["phylogeo"]},
            "report":      {"name": "报告", "desc": "HTML 汇总报告", "inputs": ["全模块产物"],
                            "outputs": ["08_report"], "risks": [], "requires": []},
        },
    },
    "eve": {
        "num": "⑥", "name": "Endogenous Viruses (EVE) (archived 2026-10-11, method immature)",
        "script": "archive/20261011_endogenous_virus_pipeline/eve_screen.py",
        "desc": "【已归档, 不可运行】宿主基因组 EVE 正向筛查 (输入是基因组 FASTA, 与 reads 管线相互独立); 归档原因与恢复方法见 archive/20261011_endogenous_virus_pipeline/ARCHIVED_20261011.md",
        "stages": {
            "discover": {"name": "位点发现", "desc": "滑窗 blastx vs 病毒参考蛋白 → 位点合并",
                         "inputs": ["宿主基因组 FASTA"], "outputs": ["01_Loci/<NAME>"],
                         "risks": ["--fast 快 6 倍但位点数约 -54%"], "requires": []},
            "verdict":  {"name": "双侧判定", "desc": "植物/病毒双侧 bitscore 判定",
                         "inputs": ["loci"], "outputs": ["02_Verdict"],
                         "risks": ["SwissProt 拆分库 + id2div 映射须配齐"], "requires": ["discover"]},
            "annotate": {"name": "RVDB 归属", "desc": "RVDB 深度归属注释",
                         "inputs": ["verdict"], "outputs": ["03_RVDB"], "risks": [],
                         "requires": ["verdict"]},
            "merge":    {"name": "汇总", "desc": "kingdom_summary.tsv / family_by_genome.tsv (--merge 单独重跑)",
                         "inputs": ["全基因组产物"], "outputs": ["04_Summary"], "risks": [],
                         "requires": ["discover"]},
        },
    },
    "submission": {
        "num": "⑦", "name": "Data Submission",
        "script": "virome_submission_pipeline/submission_pipeline.py",
        "desc": "拓扑 → 元数据模板 → 假定蛋白注释 → Sequin/tbl2asn → .sqn",
        "stages": {
            "topology":    {"name": "拓扑判断", "desc": "BWA-MEM2 末端比对 → circular/linear",
                            "inputs": ["vOTU/全长"], "outputs": ["拓扑表"], "risks": [],
                            "requires": []},
            "metadata":    {"name": "元数据模板", "desc": "统一元数据模板 (source.src + features + organism)",
                            "inputs": ["检出/分类表"], "outputs": ["模板"],
                            "risks": ["与 ⑦ 列定义同步 (unified_metadata)"], "requires": ["topology"]},
            "hypothetical": {"name": "假定蛋白", "desc": "HHsuite/Diamond/DeepLoc/PSORTb/TMHMM 5 工具功能注释",
                             "inputs": ["蛋白"], "outputs": ["注释表"], "risks": ["HHsuite 库大"],
                             "requires": ["metadata"]},
            "sequin":      {"name": "Sequin 构建", "desc": "Sequin .tbl 格式 (Cenote-Taker3 风格)",
                            "inputs": ["模板"], "outputs": [".tbl"], "risks": [], "requires": ["metadata"]},
            "submit":      {"name": "端到端", "desc": "suvtk / Sequin tbl2asn 双模式 → .sqn",
                            "inputs": [".tbl"], "outputs": [".sqn 提交文件"],
                            "risks": ["tbl2asn 需单独部署"], "requires": ["sequin"]},
            "report":      {"name": "报告", "desc": "交互式 HTML 全表编辑报告",
                            "inputs": ["全阶段"], "outputs": ["HTML"], "risks": [], "requires": []},
        },
    },
    "host_expression": {
        "num": "⑧", "name": "Host Expression Analysis (archived 2026-10-11, moved out 2026-10-27)",
        "script": "(已移出主仓 → 差异表达分析 项目归档)",
        "desc": "【已归档 2026-10-11、2026-10-27 移出主仓, 不可运行】病毒阳性/阴性分组宿主侧表达: 定量 → 标签 → TEA-GCN/DESeq2 → 富集 → 报告; 归档原因与恢复方法见 `差异表达分析` 项目归档内的 ARCHIVED_20261011.md",
        "stages": {
            "quant":    {"name": "定量", "desc": "宿主 CDS 提取 (gffread) + Salmon 定量 + QC",
                         "inputs": ["① 宿主参考", "② clean reads"], "outputs": ["02_Quant"],
                         "risks": ["输入全部取自主链现有产物"], "requires": []},
            "labels":   {"name": "病毒标签", "desc": "④ 检出表 → 阳/阴分组 + 连续载量 (合并 Core14 组织)",
                         "inputs": ["④ 检出表"], "outputs": ["04_Labels"], "risks": [],
                         "requires": ["quant"]},
            "coexp":    {"name": "差异共表达", "desc": "TEA-GCN 主分析: ΔzScore + 网络 z-QC 闸门 + DEG 限定 Louvain 模块 (内置样本量/失衡闸门)",
                         "inputs": ["定量+标签"], "outputs": ["05_Coexp"],
                         "risks": ["样本量闸门可能 blocked —— 属预期保护"], "requires": ["labels"]},
            "de":       {"name": "差异表达", "desc": "DESeq2 交叉参考 (可选组织分层)",
                         "inputs": ["定量+标签"], "outputs": ["06_DE"], "risks": [],
                         "requires": ["labels"]},
            "loadscan": {"name": "载量扫描", "desc": "载量 × 全基因组剂量验证",
                         "inputs": ["载量+定量"], "outputs": ["08_LoadScan"], "risks": [],
                         "requires": ["labels"]},
            "enrich":   {"name": "富集", "desc": "功能富集分析", "inputs": ["DE/模块基因"],
                         "outputs": ["07_Enrich"], "risks": [], "requires": ["de"]},
            "report":   {"name": "报告", "desc": "HTML 汇总报告", "inputs": ["全模块"],
                         "outputs": ["09_Report"], "risks": [], "requires": []},
        },
    },
}


# ═══════════════════════════════════════════════════════════════════
# 查询接口
# ═══════════════════════════════════════════════════════════════════

def pipelines() -> List[str]:
    """全部管线 key (顺序即主链顺序)。"""
    return list(STAGE_CARDS)


def get_card(pipeline: str, stage: Optional[str] = None) -> Optional[Dict]:
    """取管线卡片; 给 stage 则取单张阶段卡。找不到返回 None。"""
    p = STAGE_CARDS.get(pipeline)
    if p is None:
        return None
    if stage is None:
        return p
    return p["stages"].get(stage)


def stages_of(pipeline: str) -> List[str]:
    p = STAGE_CARDS.get(pipeline)
    return list(p["stages"]) if p else []


def summary_text(max_stage_desc_chars: int = 60) -> str:
    """紧凑总览 (喂 AI 提示词 / chat 意图识别的上下文)。"""
    lines: List[str] = []
    for key, p in STAGE_CARDS.items():
        lines.append(f"{p['num']} {key} ({p['name']}): {p['desc']}")
        for s, c in p["stages"].items():
            desc = c["desc"][:max_stage_desc_chars]
            lines.append(f"  - {s}: {desc}")
    return "\n".join(lines)


def self_check() -> List[str]:
    """结构自检: 卡片缺字段/前置 stage 不存在等。返回问题列表 (空=通过)。"""
    problems: List[str] = []
    for pk, p in STAGE_CARDS.items():
        for field in ("num", "name", "script", "desc", "stages"):
            if field not in p:
                problems.append(f"{pk}: 缺字段 {field}")
        for sk, c in p.get("stages", {}).items():
            for field in ("name", "desc", "inputs", "outputs", "risks", "requires"):
                if field not in c:
                    problems.append(f"{pk}.{sk}: 缺字段 {field}")
            for req in c.get("requires", []):
                if req not in p.get("stages", {}):
                    problems.append(f"{pk}.{sk}: requires 未知的 stage {req}")
    return problems
