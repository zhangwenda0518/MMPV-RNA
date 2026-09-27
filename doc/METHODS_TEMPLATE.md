# MMPV-RNA 方法章总模板（论文级 · 中文版）

> **用途**：本文档汇总 MMPV-RNA 各管道的论文级方法段，供使用 MMPV-RNA 的研究者在投稿时直接取用。按论文实际使用的管道组合挑选对应小节粘贴进 Methods，将文末占位符替换为自身数据即可。
>
> **如何使用**：按研究设计组合选取小节——①仅做病毒发现 → 管道一＋管道二（如需提交序列再加管道五）；②已知病毒定量/变异/群体遗传 → 管道三；③系统发育/定年/重组 → 管道四；④序列提交 → 管道五。管道一（公共数据获取与宿主库构建）为公共前置，仅在使用了公共数据或需声明宿主去除策略时写入。粘贴后统一替换文末"占位符清单"中的条目。
>
> **来源文件**：本模板为各管道独立方法章的合并精简版，逐管道源文件见 `virome_discovery_pipeline/METHODS_Virome_Discovery.md`（管道一/二）、`virome_analysis_pipeline/METHODS.md`（管道三）、`virome_phylo_pipeline/METHODS.md`（管道四）、`virome_submission_pipeline/METHODS.md`（管道五）。
>
> **数据口径**：所有工具名、参数、阈值均取自各管道脚本的实际调用（argparse 默认值 / 代码逻辑）及各管道独立方法章，非编造。软件版本号标 ⊙ 者待从服务器运行环境确认（本机 `SOFTWARE_VERSIONS.txt` 及管道独立方法章中已核实的直接写明）。文末参考文献按字母排序，覆盖全部管道引用的软件。

---

## 引用声明（Citation / Acknowledgment）

建议在 Methods 或 Acknowledgment 中按以下形式声明：

> Viral genomes were identified, assembled, and analyzed using MMPV-RNA (v2.3) [GitHub/DOI 待填], a modular metatranscriptomic virome analysis platform.

（GitHub 地址与 DOI 按实际发布信息填写；v2.3 与代码内版本一致——`virome_pipeline.py`、`pixi.toml` 均为 v2.3。注意仓库 README 曾以 "MMPV v3.0" 作为平台展示名，投稿引用版本号以实际发布信息为准。）

---

## 管道一：公共测序数据获取与宿主参考库构建

### 公共测序数据获取

公共宏转录组数据集经五阶段元数据管道（`public_data_pipeline.py`，v3.1）获取，以宿主物种拉丁名与 NCBI TaxID 为输入。检索阶段并行检索 NCBI SRA（E-utilities esearch/efetch）与 CNCB GSA（网页爬虫），限定转录组文库（"biomol rna"），合并去重为统一登录号表；详细模式解析每个 Run 的 SRA XML 与 GSA Excel 记录，提取组织、发育阶段与地点等深层特征，并用大语言模型（DeepSeek）辅助清洗自由文本字段。元数据统一阶段将每个 Run 归一到核心表（Run、发布日期、采集日期、地点、来源、组织、发育阶段、物种学名、TaxID、文库类型、测序中心、BioProject、PMID），规则提取与大语言模型推理（DeepSeek/Kimi）仲裁互补，BioProject 追溯至 PubMed 文献。下载阶段经 NGDC（CRR 登录号，aria2c，FTP/HTTP 回退）与 NCBI（SRR/ERR/DRR 登录号，SRA Toolkit prefetch）双通道下载，均支持断点续传。转换阶段用 fasterq-dump 将 SRA 档案转为 FASTQ。可视化阶段生成 3×2 面板出版图，汇总样本队列的时间、数据库、机构、组织、地理与发育阶段分布。

### 宿主参考库构建

为支撑宿主 reads 去除，为每个宿主物种构建竞争性宿主参考库（`build_host_pipeline.py`，v3.1）。用 NCBI datasets CLI 下载参考基因组组装（核基因组、GFF3 注释、序列报告），合并为单一去重 FASTA。随后编译四套互补索引：Kraken2 数据库（含 taxonomy 下载或软链接，标准库含古菌、细菌、质粒、真菌、原生动物与 UniVec，非模式物种宿主 TaxID 注入序列头）、Bowtie2 索引、HISAT2 索引，以及 Minimap2 索引（dna-short、rna-short、nanopore、pacbio 四套 read 预设）。

---

## 管道二：宏转录组病毒组发现与注释

### 流程概述

MMPV-RNA 为端到端宏转录组病毒组发现管线（`virome_pipeline.py`，v2.3），以 Python 编排，依次经历十五个阶段：质量控制（clean）、宿主去除（deplete）、组装（assembly）、多工具病毒鉴定（identification）、层级假阳性过滤（filter）、延伸（cobra）、跨样本合并（merge）、聚类（cluster）、分类学分配（taxonomy）、宿主预测（host）、质量评估（checkv）、拯救（rescue）、下游分析（analysis）、证据整合验证（analysis_verify）与报告（report），支持检查点续跑与全程溯源日志。

### 数据预处理与宿主 reads 去除

原始 reads 用 fastp（Chen et al., 2018）做接头修剪与质量过滤，seqkit（Shen et al., 2016）转为 FASTA，Clumpify（BBTools）做 k-mer 聚类重排，以提升压缩率并加速下游重叠组装（管线不做 reads 去重；fastp 自带去重默认未开启）。宿主 reads 分三步去除：(i) Kraken2（Wood et al., 2019）粗分类过滤（confidence 0.2，保留病毒 taxid 10239）；(ii) Bowtie2（Langmead & Salzberg, 2012）对宿主基因组精细比对（可替换 HISAT2 或 minimap2），丢弃一致比对上的 read 对；(iii) RiboDetector（Deng et al., 2022）去除核糖体 RNA。（可选的 k-mer 丰度归一化 BBNorm 已移至数据预处理管线，默认不执行。）

### 从头组装与候选病毒鉴定

组装用 rnaSPAdes（rnaviralspades 模式；Bushmanova et al., 2019）、MEGAHIT（Li et al., 2015）或 Penguin 逐样本生成，丢弃短于最小长度阈值的 contig。病毒 contig 由十个互补工具组成的共识面板鉴定：geNomad（Camargo et al., 2024）、DIAMOND BLASTX（对病毒蛋白、nr、UniProt 三个参考库；Buchfink et al., 2021）、RdRp-Catch、ViraLM（Peng et al., 2024）、VirBot（Chen et al., 2023）、VirSorter2（Guo et al., 2021；dsDNAphage、NCLDV、RNA、ssDNA、Lavidaviridae 组）、ViralVerify、VirHunter（Sukhorukov et al., 2022）、Metabuli（Kim & Steinegger, 2024），以及针对类病毒样序列的 BLASTN 模块。200–1000 bp 的 contig 归入类病毒候选集，≥1000 bp 进入病毒目录，命中阈值 E-value ≤ 1e−5。候选序列经两层假阳性筛查：(i) DIAMOND BLASTX 对 UniRef90（Suzek et al., 2015；E-value ≤ 1e−3），top-hit 多数投票限定病毒 taxid 与策展病毒关键词；(ii) 翻译后 MMseqs2（Steinegger & Söding, 2017）检索保守结构域数据库（CDD；Marchler-Bauer et al., 2017），按结构域证据将 contig 分为病毒结构域（tier 1）、宿主蛋白（tier 2）或根层级（tier 3），在默认（拯救导向）、strict、raw 三种模式应用。

### 序列延伸、跨样本合并与聚类

候选 contig 经 BWA-MEM2 比对与 CoverM 深度估计后，用 COBRA（Chen & Banfield, 2024；k-mer 范围 21–141，linkage mismatch 2）对原始 reads 延伸。逐样本延伸结果合并，片段化延伸经 metaFlye（`--subassemblies`；Kolmogorov et al., 2020）脚手架化，并记录 read 到样本溯源。去冗余分两层：CD-HIT（Fu et al., 2012）对合并的 ICTV/NCBI 参考基因组集做参考引导预聚类（ANI ≥ 0.95，query coverage ≥ 0.85），随后仅对新序列做基于 Leiden 的图聚类（vclust），产出代表质心（final_centroids.fasta）。

### 分类学分配与宿主预测

病毒分类学由八工具集成分配：geNomad、Metabuli、CAT/BAT（von Meijenfeldt et al., 2019）、DIAMOND LCA、VITAP、MMseqs2、ACVirus 与 vConTACT3，各自对标策展数据库（RVDB v31、VMR MSL40、ACVirus 库）。各工具八阶元谱系经加权投票共识引擎（R）整合，输出带一致性与置信度的最终分配。宿主预测遵循优先级决策树：ICTV 认可的宿主注释 > RNAVirHost > PhaBOX2/CHERRY。

### 基因组完整性评估与高质量目录构建

完整性用 CheckV（Nayfach et al., 2021）评估。四分支级联拯救不完整基因组：(A) CheckV 完整性 ≥ 90% 直接通过；(B) 失败基因组用 Virseqimprover（内部 Salmon 定量加 RagTag 脚手架；Alonge et al., 2022）迭代 read 延伸后复评；(C) dc-megablast BLASTN 对自定义病毒库做参考引导延伸并补洞；(D) 属级长度恢复，推断基因组落在同属平均长度 ±15% 内。高质量（HQ）植物病毒操作分类单元（vOTU）跨分支用 vclust 去重，产出最终 HQ 目录。

### 下游分析、注释验证与证据整合

HQ 目录经 suvtk 做结构注释（分类学推断、特征表提取），假定蛋白经 DIAMOND/HMMER（Eddy, 2011）对 RVDB 检索验证。病毒丰度与流行率用 Salmon（Patro et al., 2017）跨全部输入文库定量。新病毒声明经五层证据链裁决：(i) 拯救序列的 CDD 结构域验证；(ii) 核苷酸/氨基酸一致性聚合；(iii) 加权证据打分（score = 0.30·BLASTN + 0.30·BLASTX + 0.40·domain，domain = max(CDD, CT3-HMM 病毒全长 profile 证据)，互补 profile 证据填补 CDD 盲区且不双重计分）给出 KEEP/REVIEW/DROP 判定；(iv) 保留过滤（≥1000 bp）；(v) ACVirus 再分类辅以科级系统发育定位。GenBank 就绪提交包（.sqn）用 tbl2asn 生成，采用 NCBI 验证器要求的本地序列 ID 格式。类病毒样候选环化经末端自 BLASTN（identity ≥ 95%，query coverage ≥ 90%）推断。所有阶段输出结构化 TSV 汇总，整合进带 Sankey 流可视化的交互式 HTML 报告。

---

## 管道三：已知病毒定量、变异与群体遗传分析

### 病毒检测与定量

用 Salmon v2.5.1（Patro et al., 2017）伪比对将 reads 比对至病毒参考序列库，计算每样本×病毒的丰度与覆盖指标。判定参数：最小覆盖度 10%、最小泊松比 0.3、最小平均深度 0.5X、最小 TPM 1.0、最小唯一比对 reads 10、物种判定 ANI 阈值 95%。输出每样本×病毒最佳命中表，字段含覆盖率、平均深度、EM 分配 reads、唯一 reads 比、reads-参考一致性（ANI）、核苷酸多样性（π）、CPM/RPM/FPKM/TPM、相对丰度、预测支持度与泊松打假得分。

### 高置信过滤

对最佳命中表做多维阈值过滤，保留高置信阳性：最小覆盖度 50%、最小平均深度 5X。通过者输出 high_conf.summary.tsv，逐病毒/逐样本聚合统计。

### 变异检测与功能注释

对每个高置信样本×病毒组合，从 BAM 提取病毒 reads，比对后调用变异（默认 iVar v1.4.4；Grubaugh et al., 2019：`samtools mpileup -aa -A -d 0 -B -Q 0 | ivar variants -q 20 -t 0.01`），TSV 经 ivar_tsv_to_vcf 转为 8 列 sites-only VCF。备选 caller 为 freebayes（Garrison & Marth, 2012；`-p 1 --pooled-continuous`）或 LoFreq（Wilm et al., 2012）。动态深度过滤：按平均深度分档 DP 阈值（<50 → DP ≥ 10；<1000 → DP ≥ 20；否则 DP ≥ 100），等位基因频率 AF ≥ 0.05，经 bcftools（Danecek et al., 2021）`QUAL>20 && INFO/DP>=dp && INFO/AF>=frq` 双重过滤。变异注释用 SnpEff v5.4c（Cingolani et al., 2012），种群内选择压力用 SNPGenie（Nelson et al., 2015）。

### 群体遗传学分析

因 iVar 输出 sites-only（无样本基因型），管线自动从等位基因频率表（allele_frequencies.tsv）重建显式 0/1 矩阵（AF ≥ 0.05 记 1、否则记 0）与连续 AF 矩阵，避免基因型退化。基于重建矩阵计算 Jaccard 与 Hamming 距离、UPGMA/NJ 树、主成分分析、位点间连锁不平衡（r²，阈值 0.5，同位置多等位对视为伪信号剔除）、滑动窗群体遗传（π 与 Tajima's D），以及基因级突变宏（oncoprint、MAF 突变谱、dN/dS、πN/πS、密码子谱、bootstrap CI）。

### 全长组装与完整基因组提取

对每个组装样本做 de novo 全长病毒基因组组装（OmniVirusAssembler，12 步精炼：de novo [megahit/spades/penguin] + refineC + BLAST 方向校正 + Divine Fusion + PVGA 熔断 + rmDup + 迭代抛光 + gmcloser/abyss-sealer 补洞 + 环化检测），`--min_covered 10` 过滤低丰度靶标。提取每样本最长 contig，可选参考填充缺口（RNA 仅在 N < 5% 时填、DNA 总是填）。

### 正向选择与序列相似度

对病毒编码序列做密码子水平正选择检测：基因树用 IQ-TREE v3.1.2（Minh et al., 2020）构建，选择检测用 HyPhy v2.5.97（Kosakovsky Pond et al., 2020）系列（FEL、MEME、BUSTED、RELAX、CONTRASTFEL、PRIME、CLN、LABELTREE），多方法（FEL + BUSTED + MEME）交叉验证。序列相似度用 SDT（Muhire et al., 2014）计算成对一致性，SciPy 层级聚类加 Seaborn 热图呈现序列空间结构。

### 重组与缺陷干扰颗粒检测

用 ViReMa v0.29（Routh & Johnson, 2014）检测跨位置连接（重组/DVG 断点），DI-tector（Beauclair et al., 2018）判定 DVG 类型，结果汇总为重组统计矩阵，R 版 Circos / arc / top-event 图可视化。

### 统计分析

变异×宿主物种关联用单侧 Fisher 精确检验（超几何尾概率），跨位点 Benjamini-Hochberg 校正，效应量以携带率比值表示。准种型（QST）单倍型层：将样本×位点 AF 矩阵 token 化，唯一模式经 Hamming 距离（≤ 0.15）union-find 连通聚类为 QST 簇，计算各宿主物种单倍型多样性（Hd = n/(n−1)·(1−Σxᵢ²)）与私有 QST；QST 簇间性状比较用 Kruskal-Wallis 加两两 Mann-Whitney U（BH 校正）。

### 软件与版本

| 工具 | 版本 | 工具 | 版本 |
|---|---|---|---|
| Salmon | 2.5.1 | samtools/bcftools | 1.21 / 1.6 |
| iVar | 1.4.4 | freebayes | v1.3.1 |
| LoFreq | 2.1.5 | SnpEff | 5.4c |
| IQ-TREE | 3.1.2 | HyPhy | 2.5.97(MP) |
| ViReMa | 0.29 | MultiQC | 1.33.dev0 |

Python 3.10（pandas 2.3.3、numpy 1.26.4、scipy 1.15.3、sklearn 1.7.2、matplotlib 3.9.4、networkx 3.4.2、polars 1.39.3），R 4.4.3。分析在 Linux 计算服务器（246）执行（`.pixi` + `mambaforge` 环境），并行度由 `--jobs`/`--threads` 控制。

---

## 管道四：系统发育与分子进化分析

### 序列比对与建树

对每个高置信病毒，将全长基因组序列用 MAFFT v7.525（Katoh & Standley, 2013）做多序列比对（`--auto --reorder`），IQ-TREE v3.1.2（Minh et al., 2020）构建最大似然树（ModelFinder 自动选模 `-m MFP`；Kalyaanamoorthy et al., 2017，ultrafast bootstrap 1000 次）。比对前对低覆盖序列质控（coverage ≥ 85%）。

### 时间信号验证

定年前先检验时间信号：TreeTime（Sagulenko et al., 2018）root-to-tip 回归计算 R² 与进化速率 β；日期随机化检验（DRT；Duchêne et al., 2015）将采样日期随机打乱 20 次，真实 R² 需超过随机分布 95% 百分位；BETS（Duchêne et al., 2020）以 path sampling 比较异时与同步模型。定年决策：仅当 root-to-tip R² > 0.3 才进行 BEAST 定年并报告 TMRCA。

### 分子钟定年与系统地理学

对通过时间信号检验的病毒，用 BEAST v1.10.4（Suchard et al., 2018）配合 BEAGLE（Ayres et al., 2019）贝叶斯定年：uncorrelated lognormal relaxed clock（UCLN），树先验为 coalescent skyline（默认）/ constant / birth-death skyline，替代模型 HKY 或 GTR（+G4），多链独立运行（5 链 × 8 线程，每链 5×10⁷ states），ESS > 200 评估收敛，burnin 10% 后 LogCombiner 合并、TreeAnnotator 生成 MCC 树。系统地理学以采样地点为离散性状，CTMC 模型加 BSSVS（Lemey et al., 2009）推断迁移率与显著迁移路径（BF ≥ 5），辅以 RRT 与 TempMig，SpreaD3（Bielejec et al., 2016）交互式地图可视化。

### 群体遗传学

用 pypopart（Leigh & Bryant, 2015）做单倍型与多态性分析：coverage ≥ 85% 质控、单倍型鉴定、多样性指标（π、θ、S、Tajima's D、Fu's Fs）、单倍型网络（MJN 无向与 VirNA 最小生成网络有向）、Fst（位点法加 permutation）。

### 选择压力与重组

选择压力用 DnaSP 6（Rozas et al., 2017）计算 Ka/Ks、Fu & Li's D*/F*、Rm、LD 与位点频谱，codeml（PAML 4；Yang, 2007）与 HyPhy v2.5.97（Kosakovsky Pond et al., 2020）做正选择检验（FEL、MEME、BUSTED）。重组用 RDP5（Martin et al., 2021）检测断点，事件解析、IQ-TREE 拓扑切换与 SimPlot 验证，下游建树采用 mask 重组区或删除高置信（≥ 3 方法）重组序列两种策略，DnaSP Rm 交叉验证。

### 软件与版本

MAFFT 7.525、IQ-TREE 3.1.2、BEAST 1.10.4、BEAGLE ⊙、TreeTime ⊙、TempEst ⊙、TreeDater ⊙、pypopart ⊙、DnaSP 6、RDP5 ⊙、codeml ⊙、HyPhy v2.5.97(MP) ⊙（与分析管道实测一致）、VirNA ⊙、SpreaD3 ⊙。分析在 Linux 计算服务器（246）执行，Python 3.10（numpy、pandas、scipy、Biopython）；BEAST 经 `java -Xmx2048m` 调用，并行链数与线程数由 `--chains`/`--threads` 控制。

---

## 管道五：GenBank 序列提交

对 MMPV-RNA 鉴定的新病毒与已知病毒序列，用 suvtk v0.1.1（GitHub: LanderDC/suvtk，无正式发表文献）工具链准备 GenBank 提交文件（.sqn），由 `virome_submission.py` 编排并支持断点续传。新病毒序列来自发现管道 rescue 阶段的输出（final_centroids.fasta），已知病毒序列来自分析管道的全长组装阶段。流程分五步：(i) suvtk taxonomy 以 MMseqs2 LCA 分配 ICTV 分类学并预测基因组类型（`-s 0.7`）；(ii) suvtk features 做序列方向校正、pyrodigal-gv（Larralde, 2022）ORF 预测（coding-complete，CDS > 50% 基因组）、BFVD（Kim et al., 2025）蛋白功能注释，输出五列特征表 .tbl，无命中 ORF 标注 hypothetical protein 并可经 BLASTP 补充注释；(iii) 序列 ID 统一加 `lcl|` 前缀规范化；(iv) 按 GenBank 与 MIUVIG 标准准备 source.src、miuvig.tsv、assembly.tsv 与 template.sbt；(v) suvtk comments 整合结构化注释，suvtk table2asn 生成 Sequin 提交文件 submission.sqn 与验证报告 submission.val。提交前核对 taxonomy 无非 NA、CDS 无内部终止密码子、占位符已替换、分段病毒 isolate 一致、metagenomic 均为 TRUE、验证报告无 ERROR。

### 软件与版本

| 工具 | 版本 | 引用 |
|---|---|---|
| suvtk | 0.1.1 | GitHub（无正式文献） |
| MMseqs2 | bd01c22 | Steinegger & Söding, 2017 |
| pyrodigal-gv | ⊙ | Larralde, 2022 |
| BFVD | ⊙ | Kim et al., 2025 |
| table2asn | ⊙ | NCBI（无文献） |
| MEGAHIT | 1.2.9 | Li et al., 2015 |

分析在 Linux 计算服务器（246）执行，Python 3.10；suvtk 数据库（约 5 GB）经 `suvtk download-database` 部署于 `~/database/virus-db/suvtk_db/`。

---

## 参考文献

- Alonge M, Lebeigle L, Kirsche M, et al. 2022. Automated assembly scaffolding using RagTag elevates a new tomato system for high-throughput genome editing. Genome Biology 23:258.
- Ayres DL, Cummings MP, Baele G, et al. 2019. BEAGLE 3: Improved performance, scaling, and usability for a high-performance computing library for statistical phylogenetics. Systematic Biology 68:1052–1061.
- Beauclair G, Mura C, Combredet C, et al. 2018. DI-tector: defective interfering viral genomes' detector for next-generation sequencing data. RNA 24:1285–1296.
- Bielejec F, Baele G, Vrancken B, Suchard MA, Rambaut A, Lemey P. 2016. SpreaD3: Interactive visualization of spatiotemporal history and trait evolutionary processes. Molecular Biology and Evolution 33:2167–2169.
- Bin Jang H, Bolduc B, Zablocki O, et al. 2019. Taxonomic assignment of uncultivated prokaryotic virus genomes is enabled by gene-sharing networks. Nature Biotechnology 37:632–639.
- Buchfink B, Reuter K, Drost HG. 2021. Sensitive protein alignments at tree-of-life scale using DIAMOND. Nature Methods 18:366–368.
- Bushmanova E, Antipov D, Lapidus A, Prjibelski AD. 2019. rnaSPAdes: a de novo transcriptome assembler and its application to RNA-Seq data. GigaScience 8:giz100.
- Camargo AP, Roux S, Schulz F, et al. 2024. Identification of mobile genetic elements with geNomad. Nature Biotechnology 42:1303–1312.
- Chen G, Tang X, Shi M, Sun Y. 2023. VirBot: an RNA viral contig detector for metagenomic data. Bioinformatics 39:btad093.
- Chen LX, Banfield JF. 2024. COBRA improves the completeness and contiguity of viral genomes assembled from metagenomes. Nature Microbiology 9:737–750.
- Chen S, Zhou Y, Chen Y, Gu J. 2018. fastp: an ultra-fast all-in-one FASTQ preprocessor. Bioinformatics 34:i884–i890.
- Cingolani P, Platts A, Wang LL, et al. 2012. A program for annotating and predicting the effects of single nucleotide polymorphisms, SnpEff. Fly 6:80–92.
- Danecek P, Bonfield JK, Liddle J, et al. 2021. Twelve years of SAMtools and BCFtools. GigaScience 10:giab008.
- Deng ZL, Münch PC, Mreches R, McHardy AC. 2022. Rapid and accurate identification of ribosomal RNA sequences via deep learning. Nucleic Acids Research 50:e60.
- Duchêne S, Duchêne D, Holmes EC, Ho SYW. 2015. The performance of the date-randomization test in phylogenetic analyses of time-structured virus data. Molecular Biology and Evolution 32:1895–1906.
- Duchêne S, Lemey P, Stadler T, et al. 2020. Bayesian evaluation of temporal signal in measurably evolving populations. Molecular Biology and Evolution 37:3363–3379.
- Eddy SR. 2011. Accelerated profile HMM searches. PLoS Computational Biology 7:e1002195.
- Ewels P, Magnusson M, Lundin S, Käller M. 2016. MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics 32:3047–3048.
- Fu L, Niu B, Zhu Z, Wu S, Li W. 2012. CD-HIT: accelerated for clustering the next-generation sequencing data. Bioinformatics 28:3150–3152.
- Garrison E, Marth G. 2012. Haplotype-based variant detection from short-read sequencing. arXiv:1207.3907.
- Grubaugh ND, Gangavarapu K, Quick J, et al. 2019. An amplicon-based sequencing framework for accurately measuring intrahost virus diversity using PrimalSeq and iVar. Genome Biology 20:8.
- Guo J, Bolduc B, Zayed AA, et al. 2021. VirSorter2: a multi-classifier, expert-guided approach to detect diverse DNA and RNA viruses. Microbiome 9:37.
- Kalyaanamoorthy S, Minh BQ, Wong TKF, von Haeseler A, Jermiin LS. 2017. ModelFinder: fast model selection for accurate phylogenetic estimates. Nature Methods 14:587–589.
- Katoh K, Standley DM. 2013. MAFFT multiple sequence alignment software version 7: improvements in performance and usability. Molecular Biology and Evolution 30:772–780.
- Kim D, Paggi JM, Park C, Bennett C, Salzberg SL. 2019. Graph-based genome alignment and genotyping with HISAT2 and HISAT-genotype. Nature Biotechnology 37:907–915.
- Kim J, Steinegger M. 2024. Metabuli: sensitive and specific metagenomic classification via joint analysis of amino acid and DNA. Nature Methods 21:971–973.
- Kim RS, Levy Karin E, Mirdita M, Chikhi R, Steinegger M. 2025. BFVD—a large repository of predicted viral protein structures. Nucleic Acids Research 53:D340–D347.
- Kolmogorov M, Bickhart DM, Behsaz B, et al. 2020. metaFlye: scalable long-read metagenome assembly using repeat graphs. Nature Methods 17:1103–1110.
- Kosakovsky Pond SL, Poon AFY, Velazquez R, et al. 2020. HyPhy 2.5—A customizable platform for evolutionary hypothesis testing using phylogenies. Molecular Biology and Evolution 37:295–299.
- Langmead B, Salzberg SL. 2012. Fast gapped-read alignment with Bowtie 2. Nature Methods 9:357–359.
- Larralde M. 2022. Pyrodigal: Python bindings and interface to Prodigal, an efficient method for gene prediction in prokaryotes. Journal of Open Source Software 7:4296.
- Leigh JW, Bryant D. 2015. PopART: Full-feature software for haplotype network construction. Methods in Ecology and Evolution 6:1110–1116.
- Lemey P, Rambaut A, Drummond AJ, Suchard MA. 2009. Bayesian phylogeography finds its roots. PLoS Computational Biology 5:e1000520.
- Li D, Liu CM, Luo R, Sadakane K, Lam TW. 2015. MEGAHIT: an ultra-fast single-node solution for large and complex metagenomics assembly via succinct de Bruijn graph. Bioinformatics 31:1674–1676.
- Li H. 2018. Minimap2: pairwise alignment for nucleotide sequences. Bioinformatics 34:3094–3100.
- Marchler-Bauer A, Bo Y, Han L, et al. 2017. CDD/SPARCLE: functional classification of proteins via subfamily domain architectures. Nucleic Acids Research 45:D200–D203.
- Martin DP, Varsani A, Roumagnac P, et al. 2021. RDP5: a computer program for analyzing recombination in, and removing signals of recombination from, nucleotide sequence datasets. Virus Evolution 7:veaa087.
- Minh BQ, Schmidt HA, Chernomor O, et al. 2020. IQ-TREE 2: New models and efficient methods for phylogenetic inference in the genomic era. Molecular Biology and Evolution 37:1530–1534.
- Muhire BM, Varsani A, Martin DP. 2014. SDT: a virus classification tool based on pairwise sequence alignment and identity calculation. PLoS One 9:e108277.
- Nayfach S, Camargo AP, Schulz F, et al. 2021. CheckV assesses the quality and completeness of metagenome-assembled viral genomes. Nature Biotechnology 39:578–585.
- Nelson CW, Moncla LH, Hughes AL. 2015. SNPGenie: estimating evolutionary parameters to detect natural selection using pooled next-generation sequencing data. Bioinformatics 31:3709–3711.
- Peng C, Shang J, Guan J, Wang D, Sun Y. 2024. ViraLM: empowering virus discovery through the genome foundation model. Bioinformatics 40:btae704.
- Patro R, Duggal G, Love MI, Irizarry RA, Kingsford C. 2017. Salmon provides fast and bias-aware quantification of transcript expression. Nature Methods 14:417–419.
- Rambaut A, Lam TT, Max Carvalho L, Pybus OG. 2016. Exploring the temporal structure of heterochronous sequences using TempEst (formerly Path-O-Gen). Virus Evolution 2:vew007.
- Routh A, Johnson JE. 2014. Discovery of functional genomic motifs in viruses with ViReMa—a Virus Recombination Mapper—for analysis of next-generation sequencing data. Nucleic Acids Research 42:e11.
- Rozas J, Ferrer-Mata A, Sánchez-DelBarrio JC, et al. 2017. DnaSP 6: DNA sequence polymorphism analysis of large data sets. Molecular Biology and Evolution 34:3299–3302.
- Sagulenko P, Puller V, Neher RA. 2018. TreeTime: Maximum-likelihood phylodynamic analysis. Virus Evolution 4:vex042.
- Shen W, Le S, Li Y, Hu F. 2016. SeqKit: A Cross-Platform and Ultrafast Toolkit for FASTA/Q File Manipulation. PLoS One 11:e0163962.
- Steinegger M, Söding J. 2017. MMseqs2 enables sensitive protein sequence searching for the analysis of massive data sets. Nature Biotechnology 35:1026–1028.
- Suchard MA, Lemey P, Baele G, Ayres DL, Drummond AJ, Rambaut A. 2018. Bayesian phylogenetic and phylodynamic data integration using BEAST 1.10. Virus Evolution 4:vey016.
- Sukhorukov G, Khalili M, Gascuel O, et al. 2022. VirHunter: A Deep Learning-Based Method for Detection of Novel RNA Viruses in Plant Sequencing Data. Frontiers in Bioinformatics 2:867111.
- Suzek BE, Wang Y, Huang H, McGarvey PB, Wu CH. 2015. UniRef clusters: a comprehensive and scalable alternative for improving sequence similarity searches. Bioinformatics 31:926–932.
- Volz EM, Frost SDW. 2017. Scalable relaxed clock phylogenetic dating. Virus Evolution 3:vex025.
- von Meijenfeldt FAB, Arkhipova K, Cambuy DD, Coutinho FH, Dutilh BE. 2019. Robust taxonomic classification of uncharted microbial sequences and bins with CAT and BAT. Genome Biology 20:217.
- Wilm A, Aw PPK, Bertrand D, et al. 2012. LoFreq: a sequence-quality aware, ultra-sensitive variant caller for uncovering cell-population heterogeneity from high-throughput sequencing datasets. Nucleic Acids Research 40:11189–11201.
- Wood DE, Lu J, Langmead B. 2019. Improved metagenomic analysis with Kraken 2. Genome Biology 20:257.
- Yang Z. 2007. PAML 4: phylogenetic analysis by maximum likelihood. Molecular Biology and Evolution 24:1586–1591.

---

## 无正式发表文献的工具（引用方式）

**第三方工具（无同行评议论文）**——投稿时建议在软件表或脚注中标注 GitHub 仓库与访问日期：suvtk（github.com/LanderDC/suvtk）、CAPHEINE（github.com/veg/CAPHEINE）、VirNA、RdRp-Catch、ViralVerify、VITAP、ACVirus、RNAVirHost（github.com/terry-wagner/RNAvirhost）、vConTACT3（其前身 vConTACT2 文献见 Bin Jang et al., 2019）、BBTools/Clumpify（JGI，github.com/BioInfoTools/BBMap）、BWA-MEM2（github.com/bwa-mem2/bwa-mem2）、CoverM（github.com/wwood/CoverM）、Penguin、PhaBOX2（github.com/KennthShang/PhaBOX2）、CHERRY（github.com/KennthShang/CHERRY）、RVDB（fzer.github.io/rvdbtools）。NCBI 官方工具（datasets CLI、SRA Toolkit [prefetch/fasterq-dump]、tbl2asn）直接标注 NCBI 官网 URL 与访问日期即可；pypopart 为 PopART 的 Python 实现，方法文献见 Leigh & Bryant, 2015（仓库 github.com/Adamtaranto/pypopart）。

**自研脚本**（Virseqimprover、OmniVirusAssembler、vclust 等）在正文中直接引用 MMPV-RNA 平台即可，无需单独文献。

---

## 占位符清单（投稿前需替换）

**全局**
- [ ] GitHub 地址 / DOI（引用声明处）
- [ ] 样本数与宿主物种（如 "1,348 份宏转录组样本，含 *Lycium barbarum*、*L. ruthenicum* 与烟草"）
- [ ] 标 ⊙ 的软件版本号

**管道一**
- [ ] 检索物种与 TaxID、最终纳入的 Run 数

**管道二**
- [ ] 实际启用的鉴定工具子集
- [ ] 病毒/类病毒候选数、HQ vOTU 数

**管道三**
- [ ] 高置信病毒数、样本数、关联数（如 "15 个高置信病毒、285 个样本、442 条关联"）
- [ ] 各病毒流行率与组装成功率

**管道四**
- [ ] 各病毒定年数值（R²、进化速率、TMRCA 及 95% HPD）
- [ ] MCMC 链长与 burnin 比例（若与默认不同）
- [ ] 重组事件数量

**管道五**
- [ ] BioProject / BioSample / SRA 登录号
- [ ] GenBank 登录号区间（或 "will be deposited"）
