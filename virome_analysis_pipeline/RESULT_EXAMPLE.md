# MMPV Results 论文级示例（展示"总—分—总 + Table/Fig 交叉引用"规范产出）

> 这是一份**示例**，展示用既定规范对 `known_virus_all_v2` 实际结果产出的论文级 Results 结构。数字来自真实统计（标"⊙"者为需现场复核）。目的：让你看到成品效果，不是最终稿。

---

## 摘要式总览（Results 开篇）

对 1,348 份测序样本进行病毒序列组装与注释，共检出 **47 种**病毒相关的病原体信号（Table 1）。经多维度置信度过滤，**15 种病毒**被定义为高置信感染，覆盖 **285 个样本、442 条样本-病毒关联**（Table 1, Table 2）。流行率最高的病毒为 **Cytorhabdovirus sp. 'lycii'（17.0%，n=224）** 与 **Potato spindle tuber viroid（16.4%，n=216）**（Fig. 2A），二者均为枸杞属（*Lycium*）主要病原（Fig. 4）。对 15 种高置信病毒进行了全基因组组装（442/442 样本成功）与进化分析（Fig. 5–7）。

## 横向漏斗统计（Table 1）

> **Table 1. 病毒检测与注释流程的逐级筛选统计**
>
> | 流程阶段 | 落到的量 | 注释 |
> |---|---|---|
> | 总处理样本 | 1,348 | 含未检出病毒的 reads 样本 |
> | 检出样本（≥1 病毒） | 386 | Salmon 伪比对，best-hit 水平 |
> | 检出病原体种类 | 47 | Adjusted_Species 去重 |
> | 高置信关联（样本×病毒） | 442 | 通过多维阈值过滤 |
> | 高置信样本 | 285 | 携带≥1 高置信病毒 |
> | **高置信病毒（进入下游）** | **15** | 变异/组装/进化全套分析对象 |
> | 全基因组组装成功 | 442（样本） | OmniVirusAssembler 全长 |
> | 提取完整基因组 | 388（样本） | 最长 contig + N 填充 |
> | 重组 / DVG 检测 | 13（病毒） | ViReMa + DI-tector |
>
> 注：流程基于 `known_virus_all_v2`；"高置信"定义为通过 coverage/depth/reads/TPM/ANI/poisson 阈值（Stage 2）；样本×病毒为检测对，每病毒样本数见 Table 2。

**正文交叉引用**：共检测到 47 种病原体信号，但仅 15 种达到高置信标准而进入深度解析（Table 1）。高置信关联覆盖 285 个样本，占检出样本（386）的 73.8%（Fig. 2B），表明大量低丰度/单样本信号被阈值剔除，符合"宁滥勿缺"召回 + 高置信收敛的策略。

## 高置信病毒关键指标（Table 2）

> **Table 2. 15 种高置信病毒的检出规模与宿主分布**
>
> | 病毒 | n_sample | 覆盖度(%) | 平均深度 | 主要宿主 |
> |---|---|---|---|---|
> | Cytorhabdovirus sp. 'lycii' | 224 | 86.2 | 875.8 | Lycium barbarum / ruthenicum |
> | Potato spindle tuber viroid | 216 | 99.3 | 2,248.2 | Lycium（*L.* barbarum 127, *L.* ruthenicum 34）|
> | Tomato chlorotic dwarf viroid | 109 | 79.9 | 113.6 | Lycium |
> | Citrus exocortis Yucatan viroid | 26 | 92.3 | 18.2 | Lycium |
> | Grapevine-associated RNA virus 4 | 23 | 72.7 | 10.9 | 混合宿主 |
> | Potato virus H | 14 | 74.1 | 4.5 | 茄科 |
> | Potato virus M | 11 | 50.3 | 1.5 | 茄科 |
> | Stevia carlavirus 1 | 7 | 64.0 | 2.8 | 枸杞 |
> | Kalanchoe latent virus | 7 | 70.6 | 1.8 | 混合 |
> | Camellia ringspot associated virus 1 | 7 | 39.1 | 0.7 | 混合 |
> | Tobacco rattle virus | 6 | 53.4 | 32.5 | 混合 |
> | Cymbidium mosaic virus | 6 | 94.0 | 14.5 | 混合 |
> | Potyvirus sacchari | 6 | 89.3 | 6.5 | 混合 |
> | Tetranychus-associated picorna-like 1 | 6 | 52.7 | 4.3 | 混合 |
> | Lily symptomless virus | 5 | 100.0 | 172.6 | 混合 |
>
> 注：n_sample=检出该病毒的样本数（per-virus 水平，非样本×病毒对数）；覆盖度/深度为平均值；宿主列基于元数据 `ScientificName` 字段（⊙部分宿主分布需从 high_conf 逐样本复核）。

**正文交叉引用**：高置信病毒在宿主谱上呈现明显分化——枸杞属（*Lycium*）为最大宿主群，容纳 Cytorhabdovirus、PSTVd、TCDVd 等多种病毒（Table 2，Fig. 4A）；而 Potato virus H/M 与 Potyvirus sacchari 则偏向茄科/禾本科宿主（Fig. 4B）。Cytorhabdovirus 虽检出样本最多（n=224），但平均深度仅 875.8X，远低于 PSTVd 的 2,248.2X（Table 2），提示其宿主内载量低于 PSTVd。

## 分病毒画像示例：Potato spindle tuber viroid (PSTVd)

> **Table 3. PSTVd 的跨样本群体遗传与进化指标**
>
> | 指标 | 值 | 数据来源 |
> |---|---|---|
> | 高置信样本数 | 216 | high_conf.summary.tsv |
> | 变异位点数（vcf_merge 重建） | 71 | af_matrix.tsv |
> | 样本 × 位点矩阵 | 172 × 71 | snp_matrix.tsv |
> | Hamming 非零样本对 | 28,848 / 29,584 (97.5%) | distance_hamming.tsv |
> | 强连锁共突变位点对（r²≥0.5） | 39 | epistatic_co_mutations.tsv |
> | 准种型 QST 簇 | 体现宿主分化 | QST_Diversity_Stats.tsv |
> | 宿主分化 | *L. barbarum* (n=127) vs *L. ruthenicum* (n=34) | QST/元数据 |

**正文交叉引用（PSTVd 画像）**：PSTVd 在枸杞属中广泛流行（Table 2），但其准种在宿主间表现显著分化。基于 AF 重建矩阵（169 样本 × 71 位点），93% 的样本对表现出非零 Hamming 距离（28,848/29,584，Fig. 6C），且 PCA 主成分分离出与宿主物种（*L. barbarum* vs *L. ruthenicum*）对应的离散簇（Fig. 6A），提示不同宿主内的 PSTVd 准种群已发生分化。LD 分析识别出 39 对强连锁共突变位点（r²≥0.5，Fig. 6D；排除同位置多等位伪信号后），其中 37 对为跨位点真信号，分布在 RdRP 编码区（Table 3）。

**插图与图例示例**：

> **Figure 6. PSTVd 的跨样本准种结构与 LD 共突变。**
> (A) 主成分分析（PCA，基于 AF 重建基因型矩阵），点=样本，按宿主物种着色（蓝=*L. barbarum*，橙=*L. ruthenicum*），误差棒表示 PC 得分置信区间。
> (B) 样本间 Hamming 距离热图（同一聚类顺序），颜色=距离（黄=近，深蓝=远），对角线为 0。
> (C) Hamming 距离分布直方图，横轴=样本对距离，纵轴=频数；可见大多数样本对距离集中在 0.1–0.4。
> (D) 位点间 LD（r²）下三角热图，横/纵轴=基因组位置（bp），颜色= r²（浅黄=0.5，深红=1.0）；对角线上方空白。
> 数据来源：`vcf_merge/`（snp_matrix.tsv、distance_hamming.tsv、epistatic_co_mutations.tsv）；前提：iVar 位点特异性，AF≥0.05 重建矩阵。

## 插图清单（全文图表索引，供交叉引用）

| 编号 | 标题 | 来源 | 对应正文节 |
|---|---|---|---|
| Fig. 1 | 研究设计 | 01_detection | 总览 |
| Fig. 2 | 病毒组流行率 (A) 与高置信收敛 (B) | 01_detection | 总览 |
| Fig. 3 | 时间-扩散 | 01_detection | 总览/讨论 |
| Fig. 4 | 病毒-宿主关联 (A 枸杞, B 其他宿主) | metadata_association | 总览/分病毒 |
| Fig. 5 | 共感染复杂度 | metadata_association | 总览 |
| Fig. 6 | PSTVd 准种结构 + LD | vcf_merge | 分病毒 (PSTVd) |
| Fig. 7 | 全长基因组相似度全景 | 07_similarity | 分病毒 |
| Table 1 | 检测漏斗 | 01/02 | 总览 |
| Table 2 | 高置信病毒指标 | filter_stats.per_virus | 总览 |
| Table 3 | PSTVd 群体遗传指标 | vcf_merge/QST | 分病毒 (PSTVd) |

## 待复核/标注项（⊙）

1. **17 vs 15 病毒目录**：04/05/08 出现 17 病毒目录，但高置信仅 15——需核对是否病毒名带空格/下划线造成重复，或组装/相似度纳入了未进高置信的病毒。
2. **388 vs 442 提取**：提取阶段筛掉 54 个（最长 contig 长度/N 阈值），需列明阈值。
3. **13 vs 15 病毒分析**：dvg 仅 13 病毒跑到，2 个需查因（数据不足或失败）。
4. **宿主分布**：Table 2 宿主列为"主要宿主"，逐样本精确分布需从 high_conf 复核。
5. **PSTVd PCA 依据**：PCA 分离对应宿主（'*L. barbarum*' vs '*L. ruthenicum*'）需在正文明确区分"PCA 分离"与"生物学分群"，避免把方差解释当因果。
