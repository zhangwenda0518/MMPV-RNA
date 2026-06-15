# MMPV-RNA v2.3 四大管道技术路线图

---

## 管道 1: public_data_pipeline.py — 公共数据获取管道

```mermaid
flowchart LR
    A["📂 物种+TaxID"] --> B["search<br/>GSA+SRA 双引擎检索"]
    B --> C["info<br/>元数据深度解析<br/>+ 文献溯源"]
    C --> D["down<br/>aria2c/prefetch<br/>高通量下载"]
    D --> E["plot<br/>SCI 六图可视化"]
    E --> F["✅ FASTQs<br/>+ 元数据表<br/>+ PDF 图表"]

    style A fill:#607d8b,color:#fff
    style F fill:#2e7d32,color:#fff
```

| 阶段 | 脚本 | 功能 |
|------|------|------|
| search | `gsa_sra.search.py` | NGDC GSA + NCBI SRA 双引擎物种检索，DeepSeek AI 辅助生成检索词，合并去重 → `SRA_GSA_Merged_Final.csv` |
| info | `gsa_sra.info.py` | Run 号元数据批量获取 (13 核心字段)，可选 DeepSeek 文献溯源补全 → `Global_Unified_Metadata_Core13.csv` |
| down | `gsa_sra.down.py` | aria2c (NGDC) + prefetch/fasterq-dump (NCBI) 双通道下载，断点续传，并发控制 |
| plot | `gsa_sra.plot.py` | 六张 SCI 级统计图: 平台/地理/时间/组织/来源/类型 → `Combined_Landscape_Full.pdf` |

```bash
python public_data_pipeline.py --species "Lycium barbarum" --taxid 112863 \
    --deepseek-api "sk-xxx" --ncbi-api "xxx" --stage all
```

---

## 管道 2: data_preprocessing.py — 数据预处理管道

```mermaid
flowchart LR
    A["📂 Raw FASTQs"] --> B

    subgraph S1["00a CleanData"]
        direction LR
        B["Fastp QC<br/>接头+质量过滤"] --> C["Seqkit 统计"] --> D["Clumpify 去重"]
    end

    D --> E{"宿主去除?"}

    subgraph S2["00b HostDepletion"]
        direction LR
        E -->|执行| F["Kraken2<br/>快速标记"] --> G["Bowtie2<br/>严格比对"] --> H["RefineC<br/>去污染"] --> I["Ribodetector<br/>rRNA过滤"]
    end

    E -->|跳过| J
    I --> J["✅ Clean FASTQs<br/>去宿主/去rRNA"]

    style A fill:#607d8b,color:#fff
    style J fill:#2e7d32,color:#fff
```

| 步骤 | 工具 | 功能 |
|------|------|------|
| Fastp | fastp | 切除 Illumina 接头，Q<15 过滤，<50bp 丢弃，生成 JSON 质检报告 |
| Seqkit | seqkit stats | 统计每步 reads 数和碱基数 |
| Clumpify | clumpify.sh | 去除光学/PCR 重复 reads |
| Kraken2 | kraken2 | 快速分类标记宿主 reads |
| Bowtie2 | bowtie2 | 严格比对宿主参考基因组 |
| RefineC | refineC | 去除比对边界的宿主污染 |
| Ribodetector | ribodetector | 去除残留核糖体 RNA |

```bash
python data_preprocessing.py -i raw_fastqs/ -o out/ --host_ref host.fa --threads 64
```

---

## 管道 3: virome_pipeline.py — 宏病毒组端到端发现管道

```mermaid
flowchart LR
    A["📂 Clean FASTQs"] --> B

    subgraph G1["数据准备"]
        direction LR
        B["01 Assembly<br/>MEGAHIT/SPAdes"] --> C["02 Identification<br/>9工具+UniProt"] --> D["03 COBRA<br/>末端延伸"] --> E["04 CLUSTER<br/>CD-HIT+vclust"]
    end

    E --> F

    subgraph G2["注释评估"]
        direction LR
        F["05 Taxonomy<br/>9工具+R共识"] --> G["06 Host<br/>Plant靶向"] --> H["07 CheckV<br/>完整性评估"]
    end

    H --> I

    subgraph G3["拯救与报告"]
        direction LR
        I["08 Rescue<br/>A:CheckV<br/>B:VSI reads<br/>C:ragtag参考"] --> J["vclust去重"] --> K["09 Report<br/>HTML+旭日图"]
    end

    K --> L["✅ all_plant_viruses.fasta<br/>+ pipeline_report.html"]

    style A fill:#607d8b,color:#fff
    style L fill:#2e7d32,color:#fff
```

### 三支路 Rescue 细节

```mermaid
flowchart LR
    A["Plant novel<br/>centroids"] --> B["分支 A<br/>CheckV 直接评估<br/>≥90% 免拯救"]
    A --> C["分支 B<br/>VSI reads 延伸<br/>SPAdes+salmon+bowtie2<br/>scaffold-truncated"]
    A --> D["分支 C<br/>ragtag 参考引导<br/>BLASTN→blastdbcmd<br/>→ragtag scaffold"]
    B --> E["CheckV 复评"]
    C --> E
    D --> E
    E --> F["vclust 最终去重"]
    F --> G["all_plant_viruses.fasta<br/>免拯救 + rescued"]

    style A fill:#ef6c00,color:#fff
    style G fill:#2e7d32,color:#fff
```

| 阶段 | 核心功能 |
|------|----------|
| 01 Assembly | MEGAHIT/rnaviralSPAdes/Penguin de novo 组装，输出 N50/N90 统计 |
| 02 Identification | 9 工具并行鉴定 + UniProt 蛋白级过滤 (raw/filter/strict) |
| 03 COBRA | BWA-MEM2 末端延伸，统计延伸率/孤儿率 |
| 04 CLUSTER | CD-HIT 参考预聚类 + vclust (Leiden/ANI) 去冗余 |
| 05 Taxonomy | 9 工具并行分类 (mmseqs/CAT/DIAMOND 等) + R 共识筛选 |
| 06 Host | RNAVirHost + iPHoP + ICTV ensemble 宿主预测 |
| 07 CheckV | 按宿主分组完整性评估 (Complete/High/Medium/Low/NA) |
| 08 Rescue | 三支路拯救 (A: 直接 / B: VSI / C: ragtag), `--checkv_threshold` 可调 |
| 09 Report | 调用独立 `report_pipeline.py` 生成期刊级交互式 HTML |

```bash
python virome_pipeline.py --input_reads clean_fastqs/ --output_dir out/ --host-filter Plant
python virome_pipeline.py --stage rescue --output_dir out/ --checkv_threshold 90
```

---

## 管道 4: auto_known_virus.py — 已知病毒定量与变异管道

```mermaid
flowchart LR
    A["📂 FASTQs<br/>+ vOTU参考"] --> B

    subgraph S1["detect 检测"]
        direction LR
        B["Salmon<br/>索引"] --> C["伪比对<br/>TPM定量"] --> D["Poisson<br/>打假"] --> E["双轨<br/>过滤"] --> F["depth_summary.tsv"]
    end

    F --> G

    subgraph S2["variants 变异"]
        direction LR
        G["Reads<br/>提取"] --> H["共识<br/>序列"] --> I["FreeBayes<br/>SNP/INDEL"] --> J["SnpEff<br/>注释"] --> K["SnpGenie<br/>dN/dS"] --> L["variant_summary.tsv"]
    end

    L --> M

    subgraph S3["full 单倍型"]
        direction LR
        M["virus-full.py<br/>批量调度"] --> N["全基因组<br/>单倍型"]
    end

    N --> O["✅ TPM矩阵<br/>变异报告<br/>可视化图表"]

    style A fill:#607d8b,color:#fff
    style O fill:#2e7d32,color:#fff
```

| 阶段 | 工具 | 功能 |
|------|------|------|
| detect | Salmon | 伪比对快速定量 → Poisson 检验打假 → 双轨过滤 (RNA-seq 仅保留 RNA 病毒) |
| variants | FreeBayes + SnpEff + SnpGenie | SNP/INDEL 检出 → 功能注释 (同义/错义/移码) → 群体遗传 (dN/dS, π, θ) |
| full | virus-full.py | 批量单倍型全长组装 |

> **注意**: 输入 reads 建议使用仅经 Fastp QC、**未宿主去除**的 reads（宿主去除可能误删病毒 reads）。

```bash
python auto_known_virus.py --stage all --ref all_plant_viruses.fasta \
    -i clean_fastqs/ -o out_known/ --threads 64
```

---

## 四管道整体数据流

```mermaid
flowchart LR
    P1["管道 ①<br/>public_data_pipeline<br/>━━━━━━━━━━━━<br/>search→info<br/>→down→plot"] -->|"FASTQ"| P2["管道 ②<br/>data_preprocessing<br/>━━━━━━━━━━━━<br/>clean→deplete"]

    P2 -->|"Clean FASTQ"| P3["管道 ③<br/>virome_pipeline<br/>━━━━━━━━━━━━<br/>Assembly→ID→COBRA<br/>→Cluster→Tax→Host<br/>→CheckV→Rescue<br/>→Report"]

    P3 -->|"参考 FASTA"| P4["管道 ④<br/>auto_known_virus<br/>━━━━━━━━━━━━<br/>detect→variants→full"]

    P3 -->|"|"| OUT

    subgraph OUT["最终产出"]
        direction TB
        O1["all_plant_viruses.fasta<br/>完整植物病毒基因组集"]
        O2["pipeline_report.html<br/>期刊级交互式报告"]
        O3["depth_summary.tsv<br/>TPM/深度/覆盖度矩阵"]
        O4["variant_summary.tsv<br/>群体变异统计"]
    end

    style P1 fill:#eceff1,stroke:#607d8b,color:#333
    style P2 fill:#e3f2fd,stroke:#1565c0,color:#333
    style P3 fill:#e8eaf6,stroke:#1a237e,color:#fff
    style P4 fill:#e8f5e9,stroke:#2e7d32,color:#333
```
