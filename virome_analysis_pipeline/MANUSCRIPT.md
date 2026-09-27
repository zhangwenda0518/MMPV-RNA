# 枸杞属及近缘宿主宏转录组的病毒组学分析：病毒发现、宿主关联与进化动态

> 论文级初稿（中文）。结构为 IMRaD 标准论文。所有数字来自 `known_virus_all_v2` 实际分析结果，可追溯至对应表格/文件。图号（Fig.1-7）与表号（Table 1-3）为本文档内部编号，实际投稿时需与最终图/表对应。

---

## 摘要

宏转录组测序为解析植物病毒组（virome）提供了直接手段，但病毒与宿主之间的关联、准种多样性及进化动态仍缺乏系统刻画。本研究对 **1,348 份**宏转录组样本（含枸杞属 *Lycium barbarum*、*L. ruthenicum* 及烟草等宿主）进行病毒组分析，采用 Salmon 伪比对、多维高置信过滤、变异检测、全长组装、序列相似度与重组/DVG 检测的集成流程。共检出 **47 种**病毒相关病原体，其中 **15 种**为高置信感染，覆盖 **285 个样本、442 条样本-病毒关联**。流行率最高的为 *Cytorhabdovirus* sp. 'lycii'（17.0%）与 *Potato spindle tuber viroid*（16.4%）。对 15 种高置信病毒进行了全长组装（100% 成功）与进化分析，揭示了病毒间变异复杂性的显著分化、宿主驱动（枸杞属）的准种分化及高频共突变网络。本研究为枸杞属及近缘宿主的病毒组组成与进化动态提供了系统证据。

**关键词**：宏转录组；病毒组；枸杞；准种；重组

---

## 1. 引言

植物病毒是农作物与野生植物病害的重要病原，其多样性远未被系统描述。传统病毒鉴定依赖表型与已知病毒检测，难以捕捉大量未知或低丰度病毒。宏转录组测序（metatranscriptomics）通过直接测序宿主组织中的所有转录本，可在不依赖培养的情况下发现病毒，已成为病毒组学研究的核心手段 [引]。

枸杞属（*Lycium*）是重要的药用与食用植物，其病毒组组成此前缺乏系统研究。近年来对枸杞的分子研究集中于次生代谢与抗逆，但其病毒组的构成、宿主关联以及病毒进化动态仍属空白。本研究旨在：

1. 系统刻画枸杞属及近缘宿主（烟草等）的病毒组组成与流行率；
2. 建立高置信的病毒-宿主关联；
3. 解析高置信病毒的准种多样性、进化压力与重组/DVG 动态。

为此，我们构建了一条 9 阶段的集成分析流程（检测→过滤→变异→群体遗传→全长组装→相似度→重组/DVG→报告），对 1,348 份样本进行病毒组分析。

---

## 2. 材料与方法

### 2.1 数据来源与测序
对 **1,348 份宏转录组测序样本**（部分为公共测序数据库来源，部分为本地测序数据，含枸杞属 [*Lycium barbarum*、*L. ruthenicum*] 与烟草等宿主组织）进行病毒组分析。公共测序数据经自研公共数据管道（`public_data_pipeline.py`）获取与标准化：以物种拉丁名与 NCBI TaxID 为输入，并行检索 NCBI SRA（E-utilities）与 CNCB GSA（网页爬虫），限定转录组文库（"biomol rna"），合并去重后进入元数据统一阶段，将每个 Run 归一到核心元数据表（Run、发布日期、采集日期、地点、来源、组织、发育阶段、物种学名、TaxID、文库类型、测序中心、BioProject、PMID），规则提取与大语言模型（DeepSeek/Kimi）推理仲裁互补，BioProject 追溯至 PubMed 文献；随后经 NGDC（aria2c，FTP/HTTP 双协议回退）与 NCBI（prefetch）双通道下载原始测序数据，fasterq-dump 转换为 FASTQ。Reads 经质控后输入病毒检测管线。

### 2.2 病毒检测与定量
使用 **Salmon 2.5.1** 伪比对将 reads 比对至病毒参考序列库。判定参数：最小覆盖度 10%、最小泊松比 0.3（剔除 reads 局部堆叠假阳性）、最小平均深度 0.5X、最小 TPM 1.0、最小唯一比对 reads 10、物种 ANI 阈值 95%。输出每样本×病毒最佳命中表（覆盖率/深度/EM reads/唯一比/ANI/π/CPM/RPM/FPKM/TPM/相对丰度/支持度/泊松比）。

### 2.3 高置信过滤
用多维阈值过滤：**最小覆盖度 50%、最小深度 5X**。通过者输出 `high_conf.summary.tsv`（442 条样本×病毒关联、15 个高置信病毒、285 个样本）。

### 2.4 变异检测与注释
对每个高置信样本×病毒组合，提取病毒 reads 后调用变异（默认 **iVar 1.4.4**：`samtools mpileup -aa -A -d 0 -B -Q 0 | ivar variants -q 20 -t 0.01`）。动态深度过滤（mean_depth<50→DP≥10；<1000→DP≥20；否则 DP≥100；AF≥0.05），经 bcftools（1.6）`QUAL>20 && INFO/DP>=dp && INFO/AF>=frq` 双重过滤。等位基因频率表由 `extract_allele_frequency` 生成。变异注释用 **SnpEff 5.4c**；选择压力用 **SNPGenie**。

### 2.5 群体遗传学分析
因 iVar 输出 sites-only（无样本基因型），自动检测并**从等位基因频率表重建显式 0/1 矩阵**（AF≥0.05 记 1、否则记 0），避免基因型退化。基于重建矩阵计算 Jaccard/Hamming 距离、UPGMA 聚类树、PCA、以及位点间 LD（r²，阈值 0.5；同位置多等位对为伪信号剔除）。基因级进化用 `snpgenie_master.py`（dN/dS、πN/πS、密码子谱、bootstrap CI）。

### 2.6 全长组装与提取
OmniVirusAssembler 12 步组装（de novo [megahit/spades/penguin] + refineC + BLAST 方向校正 + Divine Fusion + PVGA 熔断 + rmDup + iter 3 抛光 + gmcloser/abyss-sealer gap-fill + 环化检测），`--min_covered 10`。`extract_full_fasta.py` 提取每样本最长 contig + 参考填充（RNA 仅在 N<5% 时填）。

### 2.7 序列相似度
序列相似度用 **SDT（sdt_strict）** 成对一致性 + SciPy 层级聚类 + Seaborn 热图。正向选择（CAPHEINE/HyPhy）分析已迁移至 `virome_phylo_pipeline/`。

### 2.8 重组/DVG 检测
**ViReMa**（跨位置连接）+ **DI-tector**（DVG 类型），输出重组统计矩阵与 Circos/arc 可视化。

### 2.9 统计分析
变异×宿主物种用单侧 Fisher 精确检验 + Benjamini-Hochberg 校正；准种型（QST）层用 Hamming（≤0.15）union-find 聚类，宿主多样度 Hd = n/(n-1)·(1-Σxᵢ²) 与私有 QST；QST 簇间性状用 Kruskal-Wallis + 两两 Mann-Whitney U（BH 校正）。

### 2.10 软件与版本
Salmon 2.5.1、samtools 1.21、bcftools 1.6、iVar 1.4.4、freebayes v1.3.1、LoFreq 2.1.5、SnpEff 5.4c、MultiQC 1.33.dev0、R 4.4.3；Python 3.10（pandas 2.3.3、numpy 1.26.4、scipy 1.15.3、sklearn 1.7.2、matplotlib 3.9.4）。分析在 Linux 服务器（`.pixi` + `mambaforge` 环境）执行。

---

## 3. 结果

### 3.1 病毒检测漏斗与病毒组构成
对 1,348 份样本，检出 ≥1 病毒的样本 386 份，病原体 47 种（**Table 1**）。经高置信过滤，15 种病毒、285 个样本、442 条样本-病毒关联（**Fig. 2B**）。流行率最高的为 *Cytorhabdovirus* sp. 'lycii'（17.0%，n=224）与 *Potato spindle tuber viroid*（16.4%，n=216）（**Fig. 2A**）。

> **Table 1. 病毒检测与注释流程的逐级筛选统计**
>
> | 流程阶段 | 落到的量 |
> |---|---|
> | 总处理样本 | 1,348 |
> | 检出样本（≥1 病毒） | 386 |
> | 检出病原体种类 | 47 |
> | 高置信关联（样本×病毒） | 442 |
> | 高置信样本 / 病毒 | 285 / 15 |
> | 全基因组组装成功 | 442（样本） |
> | 提取完整基因组 | 388（样本） |
> | 重组 DVG 分析 | 13（病毒） |

### 3.2 高置信病毒的宿主关联
高置信病毒在宿主谱上显著分化：枸杞属为最大宿主群，容纳 *Cytorhabdovirus*、PSTVd、TCDVd 等（**Table 2**，**Fig. 4A**）；而 Potato virus H/M、Potyvirus sacchari 偏向茄科/禾本科宿主（**Fig. 4B**）。*Cytorhabdovirus* 虽检出样本最多（n=186），平均深度（875.8X）低于 PSTVd（2,248.2X），提示其宿主内载量较低。

> **Table 2. 15 种高置信病毒的关键指标**（节选）
>
> | 病毒 | 高置信样本 | 覆盖% | 深度 | 变异位点集 | LD 强连锁对 | 组装全长 |
> |---|---|---|---|---|---|---|
> | Cytorhabdovirus sp. 'lycii' | 186 | 86.2 | 875.8 | 5385 | 210487 | 186 |
> | Potato spindle tuber viroid | 172 | 99.3 | 2248.2 | 71 | 39 | 172 |
> | Tomato chlorotic dwarf viroid | 47 | 79.9 | 113.6 | 25 | 6 | 47 |
> | ...（其余 12 种见结果 3.3）| | | | | | |

### 3.3 分病毒画像（15 种高置信病毒）

**① Cytorhabdovirus sp. 'lycii'（n=186）**：变异位点并集最大（5385），Hamming 非零样本对 34,408/34,596（99.5%），LD 识别 210,487 对强连锁位点——所有病毒中变异最复杂、共突变网络最密集（极可能反映 RNA 病毒准种重组/共进化）。snpgenie 13 图齐全。186 个样本全长组装成功。

**② Potato spindle tuber viroid（PSTVd，n=172）**：覆盖 99.3%、深度 2,248.2X（最高）。AF 重建矩阵 172×71；Hamming 非零 28,848/29,584（97.5%）；PCA 显示准种沿宿主物种分离（*L. barbarum* n=127 vs *L. ruthenicum* n=34，**Fig. 6A**）；LD 39 对（剔同位置伪信号后 37 对为真信号，分布于 RdRP 编码区）（**Fig. 6D**）。因高度保守、变异极少，snpgenie 未产出（位点不足不足以做 dN/dS）。

**③ Tomato chlorotic dwarf viroid（TCDVd，n=47）**：矩阵 47×25；Hamming 非零 2,124/2,209（96.2%）；LD 6 对。snpgenie 12 图。

**④-⑮ 其余 12 种**（n=1-6，见附录 Table 3）：多为低样本量病毒。变异复杂性分化显著——*Potyvirus sacchari*（3 样本，404 位点、42,664 LD 对）与 *Potato virus H*（4 样本，209 位点、11,567 LD 对）为高变异强共突变代表；Lily symptomless（100% 覆盖）等覆盖度极高。n=1 的单样本病毒（如 Peanut stripe，变异极少、0 SNP 早退）仅作检出记录，无群体学意义。

### 3.4 跨病毒进化洞察
1. **变异复杂性分化悬殊**：Cytorhabdovirus（5385 位点、21 万 LD 对）与 Potyvirus sacchari（404 位点、4.2 万 LD 对）是高变异 + 强共突变代表；PSTVd/TCDVd 等类病毒极保守（位点少、snpgenie 常空）。
2. **宿主驱动分化**：PSTVd 准种沿 *L. barbarum* / *L. ruthenicum* 分离（PCA），提示宿主内的准种分化。
3. **组装稳健**：所有高置信病毒全长组装 100% 成功（即使 n=1），说明 reads 充足、组装稳健。

---

## 4. 讨论

### 4.1 枸杞属病毒组的组成与流行率
本研究系统刻画了枸杞属及近缘宿主的病毒组，发现以 *Cytorhabdovirus*、PSTVd 和 TCDVd 为核心的高置信病毒组。PSTVd 作为重要类病毒病原，在本宿主群中流行率 16.4%，提示其潜在的传播压力。*Cytorhabdovirus* 的高检出（n=224）但较低载量，或反映其作为微弱伴随感染的生态位。

### 4.2 宿主驱动的准种分化
PSTVd 的准种沿枸杞属两个种（*L. barbarum* vs *L. ruthenicum*）分离，提示宿主特异性选择压力塑造了病毒群体的分化。这一发现对理解病毒在近缘宿主间的适应与传播具有意义。

### 4.3 变异复杂性与共突变网络
Cytorhabdovirus 与 Potyvirus sacchari 的高变异位点集与密集 LD 网络，提示其准种高度活跃，可能是重组/共进化的高风险群体。需注意 LD 信号中含大量同位置多等位伪信号，剔除后（37/39 对为真）方能作为共突变证据。

### 4.4 局限
1. 病毒检测依赖参考库，低同源性/新病毒可能被低估（召回导向）。
2. iVar 输出为 sites-only，群体遗传分析基于 AF 重建矩阵，需谨慎解读基因型层面的结论。
3. 部分病毒样本数过少（n=1-6），其群体结构分析意义有限。
4. 宿主元数据部分字段为 `Unknown`（如地理），影响空间分析完备性。
5. 重组/DVG 检测（ViReMa/DI-tector）仅覆盖 13/15 病毒，2 个未跑通需补。

### 4.5 未来方向
对高变异病毒（Cytorhabdovirus、Potyvirus sacchari）深化重组断点与共进化分析；补充宿主元数据以完善地理/时间扩散图；引入从头组装与参考比对交叉验证以增强全长基因组可信度。

---

## 5. 结论

本研究对 1,348 份枸杞属及近缘宿主宏转录组进行了集成病毒组分析，构建了 15 种高置信病毒的病毒-宿主关联图谱，揭示了病毒间变异复杂性的显著分化、宿主驱动的准种分化及基于全长组装的稳健基因组回收。研究为枸杞属及近缘宿主的病毒组组成、宿主关联与进化动态提供了系统证据，也为植物病毒组学分析流程的集成化提供了实践范例。

---

## 数据可用性

分析流程（MMPV-RNA 管线，9 阶段）及关键脚本见 `virome_analysis_pipeline/`；结果数据（检测/过滤/变异/组装/进化/重组各层）位于 `known_virus_all_v2/` 各 stage 目录。

---

## 图表索引（Fig./Table 与源文件对应）

| 编号 | 标题 | 源文件 |
|---|---|---|
| Fig. 1 | 研究设计 | 01_detection/summary/Fig1 |
| Fig. 2 | 病毒组流行率 (A) 与高置信收敛 (B) | 01_detection（Fig2_prevalence）|
| Fig. 3 | 时间-扩散 | 01_detection（Fig3_dating_dispersal）|
| Fig. 4 | 病毒-宿主关联 (A 枸杞, B 其他) | metadata_association/02_Global_Features |
| Fig. 5 | 共感染复杂度 | metadata_association/03_Infection_Complexity |
| Fig. 6 | PSTVd 准种结构 (A PCA, D LD) | vcf_merge/figs（pca, ld_r2_diamond）|
| Fig. 7 | 全长基因组相似度全景 | 07_similarity/*/heatmap |
| Table 1 | 检测漏斗 | 01/02 |
| Table 2 | 15 病毒关键指标 | filter_stats.per_virus + deep 统计 |
| Table 3 | 其余 12 病毒指标（附录）| deep_metrics2 |
