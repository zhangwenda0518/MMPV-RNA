# Methods —— 病毒组学管线（MMPV-RNA）

> 论文级 Methods 章节。所有工具名、参数、阈值、统计方法均来自管线脚本的实际调用（`auto_known_virus.py` 及各 stage 脚本的 argparse/逻辑），非编造。软件具体版本号标 ⊙ 待从运行环境/日志确认（服务器 mambaforge / biosoft 环境）。

---

## 1. 数据来源与测序

对 **1,348 份宏转录组测序样本**（部分为公共测序数据库来源，部分为本地测序数据，含枸杞属 [*Lycium barbarum*、*L. ruthenicum*] 与烟草等宿主组织）进行拼接与病毒组分析。样本来源为公共测序数据库及本地测序，Reads 经质控后输入病毒检测管线。

## 2. 病毒检测与定量（Stage 1）

使用 Salmon（⊙版本）伪比对，将 reads 比对至病毒参考序列库（`ref_*` 各病毒参考，含完整基因组）。每样本×病毒计算丰度与覆盖指标。判定参数（`batch_virus_depth.py`）：

- 最小覆盖度 `--coverage` = **10%**
- 最小泊松比 `--ratio` = **0.3**（剔除 reads 局部堆叠假阳性）
- 最小平均深度 `--meandepth` = **0.5X**
- 最小 TPM `--min_tpm` = **1.0**
- 最小唯一比对 reads `--min_uniq_reads` = **10**
- 物种判定 ANI 阈值 `--sp_thresh` = **95%**

每样本×病毒的"最佳命中"记录写入 `all_viruses.best.summary.tsv`，字段含覆盖率（`Rep_Coverage`）、平均深度（`Rep_MeanDepth`）、EM 分配 reads（`Asm_EM_Reads`）、唯一reads比（`Unique`）、reads-参考一致性（`Avg_Read_ANI`）、核苷酸多样性（`Avg_Pi`）、归一化丰度（CPM/RPM/FPKM/TPM）、相对丰度（`Asm_Rel_Abund`）、预测支持度（`Predicted_Support`）、泊松打假得分（`Poisson_Ratio`）。

## 3. 高置信过滤（Stage 2）

用 `filter_summary.py` 对最佳命中表做多维阈值过滤，仅保留高置信阳性。默认阈值（`--filter_cov` / `--filter_depth`）：

- 最小覆盖度 `filter_cov` = **50%**
- 最小平均深度 `filter_depth` = **5X**

记录通过所有阈值者输出至 `high_conf.summary.tsv`（本数据集 442 条样本×病毒关联、15 个高置信病毒、285 个样本），未通过者至 `high_conf.summary.discarded.tsv`。逐病毒/逐样本聚合统计写入 `filter_stats.per_virus.tsv` / `per_sample.tsv`。

## 4. 变异检测与功能注释（Stage 3）

对每个高置信"样本×病毒"组合，从 BAM 提取该病毒 reads，比对后调用变异（**默认 caller 为 iVar**）：`samtools mpileup -aa -A -d 0 -B -Q 0 | ivar variants -q 20 -t 0.01 -r ref`，TSV 经 `ivar_tsv_to_vcf` 转为 8 列 sites-only VCF（INFO 记录 DP/AF）。备选 caller 为 freebayes（`-p 1 --pooled-continuous`）或 LoFreq（`indelqual --dindel` + `call --call-indels`）。

**动态深度过滤**：依据平均深度分档 DP 阈值（`mean_depth<50→DP≥10；<1000→DP≥20；否则 DP≥100`），等位基因频率 `AF≥0.05`；经 `bcftools filter`（`QUAL>20 && INFO/DP>=dp && INFO/AF>=frq`）软/硬双重过滤得 `filtered.vcf`。

等位基因频率表 `{sample}.{contig}.allele_frequencies.tsv` 由 `extract_allele_frequency` 生成（解析 INFO，优先 `AF=`、兜底 `AO/DP`，多等位位点每 ALT 一行）。

变异功能注释用 **SnpEff**（本地库 + NCBI efetch 兜底）。种群内选择压力用 **SNPGenie**（dN/dS 计算）。

## 5. 群体遗传学分析（Stage 4）

对每个高置信病毒，跨样本构建群体遗传学图像（`virus_vcf_pipeline.py`）。由于 iVar 输出为 sites-only（无样本基因型），管线自动检测并**从 `allele_frequencies.tsv` 重建显式 0/1 矩阵**（`snp_matrix.tsv`，AF≥0.05 记 1、否则记 0）与连续 AF 矩阵（`af_matrix.tsv`），避免 sites-only 输入导致的基因型退化。

基于重建矩阵计算：

- **距离矩阵**：Jaccard 与 Hamming 距离（`distance_{jaccard,hamming}.tsv`）
- **系统树**：UPGMA 聚类树（`tree.newick`）
- **主成分分析**：基于基因型矩阵（`pca.png`）
- **连锁不平衡**：位点间 r²（`epistatic_co_mutations.tsv`），阈值 `min_r2=0.5`，同位置多等位对（`Distance_bp==0`）视为伪信号剔除
- **滑动窗群体遗传**：π 与 Tajima's D（可选用）
- 基因级突变宏（`snpeff_analysis.py`）与 oncoprint、MAF 突变谱（`snpeff2maf.py` + `viral_maftools.R`）、选择压力（`snpgenie_master.py`，dN/dS、πN/πS、密码子谱、bootstrapped CI）。

## 6. 全长组装与完整基因组提取（Stage 5/6）

对每个组装样本进行 de novo 全长病毒基因组组装（`batch_virus_full.py` → `virus-full.py`，OmniVirusAssembler）。12 步精炼：de novo（megahit/spades/penguin，`--assembly_tools all`）→ refineC split-merge → Shiver-like BLAST 方向校正 → Divine Fusion 骨架融合 → PVGA 假连接熔断 → rmDup → 迭代抛光（`--extra_args "--iter 3 --vc-min-depth 1"`）→ gmcloser/abyss-sealer 双引擎 Gap 填充 → 环化检测。`--min_covered 10` 过滤低丰度靶标，输出 `11.Ultimate_Circular_Result.fasta`。

`extract_full_fasta.py` 提取每样本最长 contig，可选参考填充缺口（`--fill`：RNA 仅在 N<5% 时填、DNA 总是填），输出 `{sample}.{contig}.full.fasta` + `assembly_stats`。

## 7. 序列相似度全景（Stage 7）

对提取的病毒全长/基因序列做跨样本/跨宿主相似度矩阵（`virus_auto_pipeline.py`），比对方法默认 **SDT（sdt_strict）**，计算成对一致性（％ identity）。用 SciPy 层级聚类 + Seaborn 热图呈现序列空间结构（`{prefix}_heatmap`）；成对一致性分布用于判断离散基因型簇 vs 连续谱系；去重报告（`03_deduplication_reports/`）用于剔除冗余序列。

## 8. 重组与缺陷干扰颗粒（DVG）检测（Stage 8）

用 **ViReMa**（⊙版本，`batch_virema_dvg.py`）检测跨位置连接（重组/DVG 断点），用 **DI-tector** 判定 DVG 类型。参数：`--seed`、`--mindel`、`--min_cov`（⊙具体值待确认）。ViReMa 结果（`virema_results/`，含 `*_Results.txt`、`*.bed`、`*.bedgraph`）经汇总生成 `Matrix_Sample_Wise_Recombination_Statistics.csv`，并用 R 版 Circos / arc / top-event 图可视化（`Summary_Analysis_Report/`）。

## 9. 统计分析

- **变异×宿主物种关联**：对每位点做单侧 Fisher 精确检验（超几何尾概率），跨位点用 **Benjamini-Hochberg** 校正（`variant_species_association.py`），效应量以携带率比值（Rate Ratio）表示。
- **准种型（QST）单倍型层**：将样本×位点 AF 矩阵 token 化（二值或分级剂量），唯一模式经 Hamming 距离（≤0.15，共享非缺失位点计算）union-find 连通聚类为 QST 簇，按携带样本数降序编号；计算各宿主物种的单倍型多样性（Hd = n/(n-1)·(1-Σxᵢ²)）与私有 QST（`virus_haplotype_qst.py`）。
- 变异位点与宿主关联的显著性（Fisher+BH）；QST 簇间性状比较用 Kruskal-Wallis + 两两 Mann-Whitney U（BH 校正）。

## 10. 软件与版本

流程在 Linux 计算服务器执行，Python 3.10 环境。关键工具版本（从服务器运行环境确认）：

| 工具 | 版本 | 环境/位置 |
|---|---|---|
| Salmon | 2.5.1 | `.pixi` |
| samtools | 1.21 | `mambaforge` |
| bcftools | 1.6 | `mambaforge` |
| iVar | 1.4.4 | `biosoft/binary` |
| freebayes | v1.3.1 | `biosoft` (MaSuRCA) |
| LoFreq | 2.1.5 | `mambaforge` |
| SnpEff | 5.4c (2026-02-23) | `biosoft/snpEff` |
| ViReMa | (biosoft/virema) | `biosoft/virema/ViReMa.py` |
| SNPGenie | (OMEGA) | `snpgenie_master.py` 封装 |
| CAPHEINE | (biosoft/CAPHEINE) | 已迁移至 `virome_phylo_pipeline/`（2026-08-27） |
| MultiQC | 1.33.dev0 | `mambaforge` |
| Rscript | 4.4.3 | `mambaforge` |
| pandas / numpy | 2.3.3 / 1.26.4 | Python 3.10 |
| scipy / sklearn | 1.15.3 / 1.7.2 | Python 3.10 |
| matplotlib / networkx / polars | 3.9.4 / 3.4.2 / 1.39.3 | Python 3.10 |

**运行环境**：分析在 Linux 计算服务器（246）执行，含 `.pixi` 与 `mambaforge` 两个环境；并行度由 `--jobs`/`--threads` 控制。
