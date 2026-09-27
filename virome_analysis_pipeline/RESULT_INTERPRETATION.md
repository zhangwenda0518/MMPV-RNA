# MMPV 管线结果解读说明文档

> 用途：供论文写作理解。本文档基于 `virome_analysis_pipeline` 各 stage 的实际脚本逻辑（README/docstring/源码）与服务器 `~/virus/data-2026/known_virus_all_v2` 的真实产出，逐 stage 说明「做了什么、输出什么东西、每个结果字段/图的统计含义是什么、论文里如何解读」。凡标注「代码声明」者为脚本里写死的输出路径，标注「实测」者为已在服务器核对过的文件。

---

## 0. 管线总览：9 stage 数据流

```
FASTQ reads
  │
  ├─ 01_detection   快速病毒定量（Salmon/Kallisto/Bowtie2 伪比对/比对）
  │        └─ all_viruses.best.summary.tsv（每样本×病毒的最佳命中 + 丰度/覆盖/ANI）
  ├─ 02_filtering   高置信双向过滤（覆盖度/深度/Reads/TPM/ANI/泊松比多维阈值）
  │        └─ high_conf.summary.tsv（442 条高置信「样本×病毒」关联，论文核心结果）
  ├─ 03_variants    变异检测 + SnpEff 注释 + SNPGenie 选择压力（逐样本逐病毒）
  │        └─ virus-variants/（filtered.vcf + allele_frequencies.tsv）+ SnpEff + SNPGenie
  ├─ 04_post_analysis  六路并行后处理（vcf_viz / vcf_merge / snpeff_macro / maftools / snpgenie / metadata_association）
  ├─ 05_assembly    每样本 de novo 全基因组组装（OmniVirusAssembler 12 步精炼）
  │        └─ 11.Ultimate_Circular_Result.fasta（组装全长）
  ├─ 06_extraction  提取组装最长 contig + N 缺口参考填充
  │        └─ {sample}.full.fasta（每样本完整基因组序列）
  ├─ 07_similarity  病毒-病毒相似度全景（层级聚类 + 热图）
  ├─ 08_dvg         DVG 缺陷干扰颗粒与重组检测（ViReMa + DI-tector + R 版 Circos）
  └─ 09_report      汇总 HTML 报告
```

> **2026-08-27 迁移**：密码子水平正向选择分析（原 Stage 7 `07_capheine`，CAPHEINE/HyPhy 全家桶）已整体迁移至 `virome_phylo_pipeline/`，本管线重编号为 9 stage（原 08_similarity→07、09_dvg→08、10_report→09）。其历史结果（15 病毒）记录见 RESULT_15VIRUS.md 与 MANUSCRIPT.md。

**数据流主线**：01 检测（谁在谁体内，丰度多高）→ 02 过滤（哪些是可靠阳性）→ 03+04 变异与进化（群体遗传）→ 05+06 全长基因组（序列层面）→ 07 进化压力（dN/dS）→ 08 跨毒株相似度 → 09 重组/DVG（准种动态）。

---

## 1. Stage 1：快速病毒定量（01_detection）

**脚本**：`batch_virus_depth.py`（主引擎）+ `batch_plot_virus_depth.py`（编排器末尾出图）

**做什么**：把每条测序样本的 reads 比对/伪比对到病毒参考库，量化每条病毒在该样本里的相对丰度、覆盖度、平均深度，并用「泊松比（Poisson Ratio）」剔除 reads 局部堆叠造成的假阳性。

**核心逻辑要点**（docstring）：
- 智能索引复用；多引擎支持（传统比对 Bowtie2/BWA 等 + 伪比对 Kallisto/Salmon，本项目用 Salmon）。
- 泊松打假：`Poisson_Ratio` 建模剔除 reads 局部堆叠假阳性。
- 双轨过滤：`genes_cov` 活跃转录区挽救机制。
- Parquet 断点续传 + Spawn 多进程 + 资源监控报表（Max RSS）。

**关键输出（实测）**：
| 文件 | 含义 |
|---|---|
| `summary/all_viruses.best.summary.tsv` | **核心结果表**：每样本×病毒的最佳命中，含丰度/覆盖/深度/ANI/泊松比。实测 754 条确诊白名单 |
| `summary/all_viruses.raw.tsv` | 全量原始命中 |
| `summary/all_viruses.summary.tsv` | 过滤后 |
| `summary/analysis_report.txt` | 人类可读的「宏病毒定量全景报告」（实测含各病毒检出率、平均 CPM/FPKM/覆盖度/泊松得分） |
| `bam/*.sorted.bam` | 每样本比对 bam |
| `stat/*.SiteDepth.gz` | 每样本逐位点深度 |
| `sample_distribution/` | 样本×病毒计数 bar/pie、病毒出现频次 bar |
| `coabundance/` | 病毒共丰度热图（abundance/relabund/cooccurrence/spearman） |
| `depth_plots/{sample}_{tax}_{vname}_depth.pdf` | 每样本每病毒深度覆盖图 |
| `summary/Fig1_study_design.*`、`Fig2_prevalence.*`、`Fig3_dating_dispersal.*` | **论文开篇三图**：研究设计、流行率、时间-扩散 |

**字段统计含义**（best.summary.tsv，实测列）：
- `Rep_Coverage(%)`：参考序列被 reads 覆盖的比例（测序深度覆盖度）。
- `Rep_MeanDepth`：平均测序深度（X）。
- `Asm_EM_Reads`：EM 算法分配给该病毒的唯一 reads 数（丰度代理）。
- `Unique(%)`：唯一比对 reads 占比。
- `Avg_Read_ANI`：reads 与参考的平均核酸一致性（物种判定）。
- `Avg_Pi`：平均核苷酸多样性（π）。
- `Asm_CPM/RPM/FPKM/TPM`：多种标准化丰度（TPM 最常用，占该样本总转录的百万分之几）。
- `Asm_Rel_Abund(%)`：相对丰度（占所有病毒 reads 的百分比）。
- `Predicted_Support`：支持度（0-1，判断可信度）。
- `Poisson_Ratio`：泊松打假得分（越高越不像随机堆叠假阳性）。
- `Segment_Accessions`、`Molecule_type`：分段/分子类型（ssRNA 等）。

**实测检出数据（analysis_report.txt）**：1348 样本，确诊白名单 754。各病毒检出率：Cytorhabdovirus sp. 'lycii' 17.0%（229 例）、Potato spindle tuber viroid 16.4%（221 例）、Tomato chlorotic dwarf viroid 8.5%（115 例）、Citrus exocortis Yucatan viroid 1.9%、Grapevine-associated RNA virus 4 1.9%、Potato virus H 1.0% 等。

**论文解读**：这一节讲「病毒组构成与流行率」——在 1348 份（枸杞属/烟草等宿主）样本里，哪些病毒最普遍、丰度/覆盖多高、是否构成共感染。Fig1 研究设计、Fig2 流行率、Fig3 时间-扩散。

---

## 2. Stage 2：高置信过滤（02_filtering）

**脚本**：`filter_summary.py`（utils/）

**做什么**：对 01 的 best.summary.tsv 做**多维双向过滤**，剔除低覆盖、低深度、低 reads、低 TPM、低 ANI、低泊松分的记录，只保留「高置信」阳性。

**过滤维度**（`--min_cov/--min_depth/--min_reads/--min_tpm/--min_ani/--min_abund/--min_support/--min_poisson/--min_unique`）→ 映射到 KNOWN_COLS：`coverage→Rep_Coverage(%)`、`depth→Rep_MeanDepth`、`reads→Asm_EM_Reads`、`tpm→Asm_TPM`、`ani→Avg_Read_ANI`、`abund→Asm_Rel_Abund(%)`、`support→Predicted_Support`、`poisson→Poisson_Ratio`、`unique_pct→Unique(%)`。

**关键输出（实测）**：
| 文件 | 含义 |
|---|---|
| `high_conf.summary.tsv` | **核心结果表**：通过全部过滤的高置信「样本×病毒」关联。**实测 442 条记录**（443 行含表头） |
| `high_conf.summary.discarded.tsv` | 被剔除的记录 |
| `filter_stats.per_sample.tsv` | 每样本聚合统计（该样本检出多少病毒、平均覆盖/深度等） |
| `filter_stats.per_virus.tsv` | **每病毒聚合统计**：`n_samples`（该病毒阳性的样本数）、`avg_cov_%`、`avg_depth`、`total_EM_reads`、`total_TPM` |
| `filter_summary_plot.pdf/png` | 过滤前后分布图 |
| `metadata_association/` | **病毒×宿主元数据关联**（见下） |

**实测 per_virus 示例**：Cytorhabdovirus n_samples=224、avg_cov 86.2%、avg_depth 875.8；PSTVd n_samples=216、avg_cov 99.3%、avg_depth 2248.2。

**metadata_association/（病毒×宿主关联，置顶记忆核心）**：
- `Viral_Infection_with_Metadata_Summary.tsv`：442 条高置信感染×宿主元数据关联（Run/Taxonomy/Location/Tissue/Host=含 Lycium barbarum 等/BioProject 15 列）。
- `01_Global_Summary`、`02_Global_Features_vs_Viruses`、`03_Infection_Complexity_Breakdown`、`04_Virus_Specific_Profiles`：四组图。

**论文解读**：这一节确立「确定感染的病毒-宿主关系」——把这些病毒在哪些宿主（枸杞属 4 种 + 烟草等）哪些组织位点检出、是否多病毒共感染、宿主宿主间病毒谱差异讲清楚。**这是论文的结果主干**（442 条高置信关联）。

---

## 3. Stage 3：变异检测（03_variants）

**脚本**：`batch_virus_variants.py`

**做什么**：对每个「样本×病毒」阳性组合，提取该病毒的 reads → 比对 → 变异检测 → SnpEff 注释 → SNPGenie 选择压力。**逐样本逐病毒独立分析**。

**调用流程**（worker_variants）：
1. `samtools view` 提取病毒 reads（`--extract_reads`）。
2. **变异检测 caller**（`--variant_caller`，本项目默认 iVar，可选 freebayes/lofreq）：
   - freebayes：`-p 1 --pooled-continuous --min-alternate-fraction 0.01`
   - lofreq：`indelqual --dindel` + `call --call-indels`
   - iVar：`samtools mpileup -aa -A -d 0 -B -Q 0 | ivar variants -q 20 -t 0.01`，再用 `ivar_tsv_to_vcf` 转 VCF
3. **动态 DP 过滤**（`mean_depth<50→DP≥10，<1000→DP≥20，否则 DP≥100；AF≥0.05`）+ bcftools soft/hard 双滤（`QUAL>20 && INFO/DP≥dp && INFO/AF≥frq`）。
4. SnpEff 注释（本地库 + NCBI efetch 兜底）。
5. SNPGenie（`dN/dS`，V13.2 OMEGA）。
6. `worker_calc_metrics` 算 Pi（`1-max/tot`）、Shannon、ANI（NM tag）。

**关键输出（目录结构，实测）**：`03_variants/virus-variants/{virus}/{sample}/`：
| 文件 | 含义 |
|---|---|
| `{sample}.filtered.vcf` | 过滤后变异 VCF（iVar 转换，8 列 sites-only，含 INFO/DP/AF） |
| `{sample}.variants.vcf` | 变异原始 VCF |
| `{sample}.{contig}.allele_frequencies.tsv` | **每样本每病毒等位基因频率表**（CHROM/POS/REF/ALT/ALT_FREQ），AF 兜底方案的输入 |
| `{sample}.variants.tsv` | 变异表 |
| `virus-SnpEff/{sample}.{contig}.ann.vcf` + 摘要 TSV | SnpEff 注释 |
| `virus-SNPGenie/` | SNPGenie 选择压力结果 |

**注意**：03 的输出主体在 `virus-variants/`（交给 04 做 vcf_merge）、`virus-SnpEff/`（交给 04 的 snpeff_macro/maftools）、`virus-SNPGenie/`。

**论文解读**：这一节建立「病毒准种/microvariome」——每个样本里的病毒不是单一序列而是一群突变体（quasispecies）。allele_frequencies.tsv 记录每个位点的等位基因频率，是后续 QST/关联/选择分析的基础。

---

## 4. Stage 4：post-hoc 六路（04_post_analysis）

**脚本**：`auto_known_virus.py` 的 `process_one_virus`，对每个高置信病毒并行跑六路。

**六路产物**（每病毒一个子目录 `04_post_analysis/{virus}/`）：
| 子目录 | 脚本 | 产出 |
|---|---|---|
| `vcf_viz/` | `virus_variants_analyzer.py` | 单样本 VCF 可视化（9+1 图引擎：AF/深度/变异分布等） |
| `vcf_merge/` | `virus_vcf_pipeline.py` | 跨样本 merge → QC → SNP 矩阵 → Jaccard/Hamming 距离 → UPGMA 树 → PCA → LD（AF 兜底版，见下） |
| `snpeff_macro/` | `snpeff_analysis.py` | 六矩阵（含 Matrix_05 连续 AF）+ 5 图 + R ComplexHeatmap oncoprint |
| `maftools/` | `snpeff2maf.py` + `viral_maftools.R` | MAF 格式 + 突变数/oncoplot/somaticInteractions |
| `snpgenie/` | `snpgenie_master.py` | 11 图统计引擎（dN/dS、πN/πS、密码子谱） |
| `metadata_association/` | `utils/virus_metadata_plot.py` | 病毒级×元数据（宿主/组织/地理） |

**vcf_merge（变异模块重点）**——本项目已修复的关键点：
- 输入 iVar 8 列 sites-only → `virus_vcf_pipeline.py` 自动检测（GT 无 0 态）→ 改吃 `allele_frequencies.tsv` 重建**显式 0/1 矩阵**（AF≥0.05 记 1，否则记 0），并输出连续 `af_matrix.tsv`。
- 产出：`snp_matrix.tsv`（0/1）、`af_matrix.tsv`（连续 AF）、`distance_jaccard/hamming.tsv`、`tree.newick`、`qc_summary.tsv`、`epistatic_co_mutations.tsv`（LD 共突变）、`pca.png`、`distance_clustermap.png`、`afs.png`、`dendrogram.png`、`ld_r2_triangle/diamond.png`、`gene_variant_links.png`（有 GTF 时）。
- **实测**：PSTVd 重跑后 Hamming 非零对 28848/29584，LD 39 对。修复前矩阵 1/NA 退化、Hamming 全零。

**论文解读**：这一节是「病毒多样性/群体遗传/选择压力」的解读层——同一病毒在不同样本里变异谱如何、是否形成亚群（树/PCA 聚类）、哪些位点强连锁（共突变，上位效应）、每个基因的 dN/dS 是否暗示正选择。snpeff_macro 的 Matrix_05 连续 AF 矩阵是 QST/关联脚本的输入。

---

## 5. Stage 5：全基因组组装（05_assembly）

**脚本**：`batch_virus_full.py`（桥接调度）→ `virus-full.py`（OmniVirusAssembler 组装引擎）

**做什么**：对每个高置信病毒的 reads 做 de novo 全长组装（含 12 步精炼：de novo + refineC split-merge + Shiver-like BLAST 方向校正 + Divine Fusion 骨架融合 + PVGA 假连接熔断 + rmDup + 迭代抛光 + gmcloser/abyss-sealer 双引擎 Gap filling + 环化检测）。

**关键输出（实测）**：`05_assembly/{病毒}/{样本}/`：
| 文件 | 含义 |
|---|---|
| `11.Ultimate_Circular_Result.fasta` | **组装全长（环化检测）**结果 |
| `{sample}_assembly.log` | 每样本组装日志 |
| 中间 12 步（`2.refinec_merged_raw` … `10.gap_filled_final`） | 逐步精炼产物 |
| `{sample}_plot_1_base_composition.png`、`{sample}_plot_2_total_depth.png`、`{prefix}_Compact/Stacked.pdf/png` | 组装图 |

**论文解读**：这一节提供「全长病毒基因组」——用组装出的全长序列做后续相似度（08）、系统发育、序列特征分析。对比 04 的变异（基于参考）和 05 的从头组装（独立序列），可交叉验证。

---

## 6. Stage 6：提取最长 contig + 参考填充（06_extraction）

**脚本**：`extract_full_fasta.py`

**做什么**：从每个样本的组装结果里取**最长 contig**，可选**用参考序列填补 N 缺口**（RNA 仅在 N<5% 时填，DNA 总是填），产出每样本一条完整基因组 fasta。

**关键输出（实测）**：`06_extraction/{sample}.full.fasta`（如 `CRR1126134.NC.002030.1.full.fasta`），`--plot` 时 `assembly_stats.png/csv`。

**论文解读**：提供每样本的**完整基因组序列**（最长 contig + 缺口修复），用于下游的序列比对、系统发育树、跨样本全长序列比较。与 05 的组装中间产物相比，06 是「去芜存菁的一根最长最完整序列」。
对比记忆点：用户偏好「Extraction」为"extract longest contig with N-fill"，产出完整 fasta 供建树。

---

## 7. Stage 7：相似度全景（07_similarity）

**脚本**：`virus_auto_pipeline.py`

**做什么**：对提取的病毒全长/基因序列做**跨毒株相似度矩阵**（用 sdt_strict 等比对），层级聚类 + 热图 + 去重报告，回答「这些病毒序列彼此多像、如何分组」——即序列空间全景。

**关键输出（代码声明）**：`07_similarity/{病毒}/`：
| 文件 | 含义 |
|---|---|
| `{prefix}_{align_method}.csv` | 相似度矩阵（默认 sdt_strict） |
| `{prefix}_heatmap.pdf/png`、`{prefix}_distribution.pdf/png` | 相似度热图 + 分布图 |
| `04_filtered_sequences/`（`gene_{g}_NT/AA.fasta`、`concatenated_CDS_NT/AA.fasta`、`overall_WholeGenome_NT.fasta`） | 去重/过滤后的序列集 |
| `03_deduplication_reports/{gene}_{g}_NT/AA/Composite_Clustering_Report.tsv` | 去重聚类报告 |

**论文解读**：这一节放在「病毒组多样性/序列空间」——同种病毒在不同样本/宿主的序列有多大分化（94-100% ANI？），是否形成离散的基因型簇，用热图+聚类直观呈现。是「群体遗传」和「系统发育」之间的桥梁。

---

## 8. Stage 8：DVG 与重组检测（08_dvg）

**脚本**：`batch_virema_dvg.py`（整合 virema/src 全部能力：ViReMa + Compiler + Visualize + DI-tector）

**做什么**：检测**缺陷干扰颗粒（DVG）**与**病毒重组**——ViReMa 找跨位点连接（重组/DVG 断点），DI-tector 找 DVG 类型（假缺陷/复合等），R 版 Circos 可视化。

**关键输出（实测）**：`08_dvg/`：
| 文件 | 含义 |
|---|---|
| `virema_results/{病毒}/{样本}/` | ViReMa 原始结果（`*_Results.txt`、`*_Insertions.txt`、`*_Substitutions.txt`、`*.bed`、`*.bedgraph`） |
| `Summary_Analysis_Report/Aggregated_Sequence_Information.csv` | 汇总序列信息 |
| `Summary_Analysis_Report/Matrix_Sample_Wise_Recombination_Statistics.csv` | **样本级重组统计矩阵** |
| `Summary_Analysis_Report/Virus_Specific_Plots/` | 每病毒专属图（Circos/arc/top_events） |
| 图 `arc_{chrom}.png/pdf`、`top_events_{chrom}.png/pdf` | 重组/DVG 弧线图、top 事件图 |

**论文解读**：这一节回答「病毒准种里的重组/缺陷变异」——重组断点在哪、哪些样本有 DVG、DVG 类型分布。呼应记忆点：用户长期搞 RDP5/重组检测算法，08 是管线内建的重组/DVG 模块（VirPhyKit 之外的独立工具）。对 RNA 病毒进化、准种动态、致病性解释有价值。

---

## 9. Stage 9：汇总报告（09_report）

**脚本**：`auto_known_virus.py` 的 report 段 + `generate_pipeline_report.py`

**做什么**：把 01-09 的关键结果收集成一份 HTML 总报告（`Pipeline_Summary_Report.html`），含各 stage 的图表与汇总。

**关键输出（实测）**：`09_report/`：`Pipeline_Summary_Report.html`、`S1_Detection`、`S2_Filter`、`S3_Variants`、`S5_Assembly_Stats`、`S6_Post_hoc`、`S7_Similarity`、`S8_DVG` 各 stage 的产物副本。

**注意**：当前 `Pipeline_Summary_Report.html` 为 8 月 2 日版，早于本轮 vcf_merge AF 兜底修复与 15 个病毒重跑，**未收录新结果**；重跑 `--stage report` 可刷新。

---

## 10. 论文写作对照速查（哪节引用哪些结果）

| 论文章节 | 用到的 stage | 关键文件/图 |
|---|---|---|
| 引言（背景） | — | — |
| 结果一：病毒组构成与流行率 | 01 | Fig1 研究设计、Fig2 流行率、analysis_report 检出率 |
| 结果二：高置信病毒-宿主关联 | 02 | high_conf.summary.tsv（442 条）、metadata_association 四组图、per_virus 统计 |
| 结果三：感染者内病毒多样性与准种 | 03 + 04 | allele_frequencies.tsv、vcf_merge（树/PCA/LD）、snpeff_macro Matrix_05 |
| 结果四：全长基因组与序列空间 | 05 + 06 + 08 | full.fasta、similarity 热图/聚类 |
| 结果五：进化选择压力 | phylo | capheine/hyphy（FEL/MEME/BUSTED/RELAX，已迁移 virome_phylo_pipeline） |
| 结果六：重组与 DVG 动态 | 09 | recombination 矩阵、DVG 类型、Circos/arc 图 |
| 讨论 | 跨 stage | 交叉验证（参考变异 vs 从头组装） |

---

## 11. 重要技术说明（写作时避免踩坑）

1. **样本×病毒 vs 每病毒样本数**：`high_conf.summary.tsv` 是「样本×病毒」对（442 条）；`filter_stats.per_virus.tsv` 是「每病毒阳性样本数」（如 Cytorhabdovirus 224、PSTVd 216）。两者不冲突，写作时明确所指层级。
2. **参考基因组 vs 从头组装**：03/04 基于参考（freebayes/iVar 比对 → 变异）；05/06 是 de novo 组装（无参考）。两者可交叉验证一致性。
3. **变异 caller 输出差异**：iVar 转 VCF 为 8 列 sites-only（仅变异位点+AF，无样本基因型）；这在 vcf_merge 里已被 AF 兜底修复（改吃 allele_frequencies.tsv 重建显式 0/1 矩阵）。写作时不要直接用旧版 GT 矩阵描述样本基因型。
4. **宿主列名**：元数据 `Global_Unified_Metadata_Core14.tsv` 的宿主信息在 `ScientificName` 列（非 `Host`），两个新脚本已做自动检测（Host→ScientificName→Species）。
5. **产出层级**：01/02 是全局汇总（跨所有样本）；03/04/05/06/07/08 是逐病毒（甚至逐样本）子目录；08 是全局 + 每病毒。写作时按需取层级。

---

## 12. 快速引用（每个 stage 的一句话定性）

- **01 检测**：Salmon 伪比对 + 泊松打假，1348 样本 → 754 确诊白名单，给出各病毒流行率与丰度。
- **02 过滤**：多维阈值筛出 442 条高置信「样本×病毒」关联，附宿主/组织/地理元数据。
- **03 变异**：逐样本逐病毒 iVar 变异检测 + SnpEff 注释 + SNPGenie。
- **04 post**：vcf_merge（AF 兜底重建）+ snpeff_macro + maftools + snpgenie + vcf_viz + metadata_association。
- **05 组装**：OmniVirusAssembler 12 步全长组装。
- **06 提取**：最长 contig + N-fill 完整基因组。
- **选择压力（原 Stage 7，已迁移 virome_phylo_pipeline）**：CAPHEINE + HyPhy 密码子正向选择。
- **07 相似度**：SDT 相似度矩阵 + 层级聚类去重，病毒序列空间全景。
- **08 DVG**：ViReMa + DI-tector 重组/DVG 检测 + Circos 可视化。
- **09 报告**：HTML 汇总（8/2 版待刷新）。
