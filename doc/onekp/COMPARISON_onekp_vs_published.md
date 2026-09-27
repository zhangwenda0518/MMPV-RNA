# OneKP 病毒组结果 vs. 已发表千种植物转录组研究 — 对比分析

对比对象：`D:\桌面\文章写作\参考文献\千种转录组` 下 5 篇文献 vs. 我们的 MMPV-RNA v2.3 OneKP 全量分析（1,342 转录组，`~/onekp-virome/`）。文献要点直接提取自 PDF 原文（页码可溯）；我方数字出自 `RESULTS_virus_discovery.md` / `RESULTS_known_virus_onekp.md`（均含来源 TSV）。

## 1. 五篇文献定位

| 文献 | 期刊/年 | 一句话定位 |
|---|---|---|
| One Thousand Plant Transcriptomes… (Leebens-Mack et al.) | Nature 2019 | OneKP 数据集本体：**1,342 个转录组 / 1,124 物种**，植物系统发育，不做病毒 |
| Mifsud et al., Transcriptome mining expands knowledge of RNA viruses across the plant kingdom | J. Virol. 2022 | **最直接对标**：挖 1KP 公开转录组找 RNA 病毒，RdRp 系统发育为核心 |
| Vogrinec et al., Wild and globally traded ornamental aquatic plants… | Environ. Microbiome 2025 | 子集深挖：79 种水生植物（复用 1KP 组装 contigs）+ 14 种自测观赏水生植物 |
| Higuita et al., An integrated analysis of the Passifloraceae virome | 预印本/53页 | 单科聚焦：114 个 SRA 数据集（12 种 Passiflora），文献+数据整合 |
| Zhang et al., An atlas of plant viruses… (PVD 数据库) | Plant Disease | 参考数据库：3,353 病毒物种、9,010 病毒–宿主对（可作我们参考库的对照基准） |

## 2. 数据集口径对照

| | 我们 (MMPV) | Mifsud 2022 | Vogrinec 2025 | Higuita (Passifloraceae) |
|---|---|---|---|---|
| 样本范围 | **全量 1,342 转录组**（=OneKP 原始口径，与 Nature 2019 的 1,342 完全一致） | 1,079 个公开转录组 / 960 物种（100 个未公开被弃） | 79 种水生植物的 1KP 组装 contigs（未自行组装）+ 14 种自测 | 114 个 SRA 数据集 / 13 BioProject / 12 物种 |
| 测序量 | raw 5.83 Tbp（中位 ~4.3 Gbp/样本，data_summary.tsv） | 平均 1.99 Gbp/样本 | 复用 1KP contigs，无 reads 层 | 未统一（SRA 各异） |
| 组装 | rnaviralSPAdes 逐样本 + co-assembly，共 **1.80 亿条 contigs**，COBRA+四分支 rescue | de novo 组装，**4,126 万条 contigs**（中位 36,015/库），主版本组装 | 直接下载 1KP contigs，BLASTX 复筛（≥500 nt） | SRA 逐数据集组装 |
| 鉴定 | 10 工具面板（geNomad/DIAMOND/Metabuli/RdRp-Catch/ViralVerify…）+ 两层去假阳 + CDD 分层 | DIAMOND BLASTx 同源检索 + RdRp ML 系统发育（人工策展） | BLASTX + 家族级归类 + CLC 去冗余 | 文献 + GenBank + SRA 组装 |

**要点**：我们是唯一跑全量 1,342 个转录组、且覆盖 DNA+RNA+类病毒全谱的流水线；Mifsud 的 RNA-only 策略与我们的植物 RNA 病毒子集最具可比性。

## 3. 结果规模对照

| 指标 | 我们 | Mifsud 2022 | Vogrinec 2025 | Higuita | Atlas/PVD |
|---|---|---|---|---|---|
| 病毒阳性范围 | 原始检出 887/1,341 样本（66.1%）；高置信 251 样本（18.7%） | 病毒样转录物见于 **603/960 物种**（62.8%） | 水生植物文献此前仅 31 病毒/18 物种；本研究扩至多科 | 12 种 Passiflora | 全库 3,353 物种 |
| 检出分类群 | 植物 vOTU 20,054；高置信 71 taxa / 118 对 | ~30 科病毒样转录物；+ssRNA 占 61%（Betaflexi 30%/Poty 19%/Seco 16%/Alphaflexi 10%），RT-dsDNA 22% | 已知作物病原（TuYV、CMV、Lettuce chlorosis virus）+ 多科新种 contigs | 12 已知 + 6 新关联 + 6 推定新种 | 1,986/3,353 可感染作物 |
| **新病毒** | **184 新种候选**（contig 级，证据链）+ 843 远缘候选 | **104 条新病毒 RdRp contigs**（物种级策展，已存 GenBank；40% +ssRNA、1/3 dsRNA、其余 −ssRNA 含新科 Viridisbunyaviridae） | 多个推定新种（含 Sagittaria 新 potyvirus、新 begomovirus） | 6 推定新种 | — |
| 独有维度 | **变异/群体遗传（22,377 变异）、全长恢复（76 条）、DVG 筛查、共感染矩阵、类病毒（CEVd 398 样本）** | 共进化分析：跨物种传播 65%（46–79%）为主 | 贸易流通风险视角 | 混合感染+跨国传播网络 | 宿主域比较（dsDNA 宿主域最广） |

**口径警示**：我们的"184 新种候选"是 contig 级（未逐条人工策展、未提交 GenBank），Mifsud 的"104"是 RdRp 物种级策展单元——**不能直接比大小**；论文中应写为"contig-level candidates"对照"curated species-level novelties"。

## 4. 分类学交叉验证（同一数据集，两条独立路线）

- **Deltapartitivirus 双重验证**：Mifsud 在低等植物发现 14 条 plant-associated partiti 序列（10 条属 Deltapartitivirus，裸子植物具 dsRNA3）；我们植物子集第一大科恰为 Partitiviridae（2,608 条 contigs）、第一大属 Deltapartitivirus（1,608），且高置信集含 *Deltapartitivirus* sp. 'yunanense'——两条路线在 dsRNA/Durnavirales 上完全收敛。
- **Potyviridae / Flexiviridae 一致性**：他们的 +ssRNA 前四科（Betaflexi/Poty/Seco/Alphaflexi）与我们植物子集 Top（Potyvirus 1,712、Caulimovirus 类 2,117）和高置信集（PVY、Lily mottle、CymMV、PVX、PVS、Schlumbergera virus X、Arabis mosaic virus/Nepovirus）一致。
- **Varicosavirus**：Mifsud 首报蕨类/苔类双段 varicosavirus-like（TfVV/MgVV）；我们高置信集含 *Varicosavirus lactucae*（2 对）——同科同形态。
- **Amalgaviridae**：我们检出 Cannabis sativa amalgavirus 1；Vogrinec 引述海草 amalgavirus——低丰度一致性类群。
- **Begomovirus**：我们 441 条 contigs；Vogrinec 在观赏水生植物发现已知+新 begomovirus。
- **Lettuce chlorosis virus / CMV**：与我们高置信集逐字同名（Vogrinec 亦报）——三库独立互证。
- **CMV satellite CARNA-5**：我们单样本检出且进入全长组装；文献少有同期报道。

## 5. 污染信号互证（重要方法论共识）

Mifsud 明确量化了 1KP 文库污染：influenza A（16 库）、human mastadenovirus C（30 库）、HIV（15 库）、PIV-5（3 库），另有 11 库 18S "worrisome contamination"；并因疑似污染剔除 Phycodnaviridae 与大部分 Mimiviridae。

我们的全局分类（taxonomy_composition.tsv）独立重现了同一景观：**Lentivirus humimdef1（=HIV）18,271 条 vOTU**、Adenoviridae 13,154、Orthoherpesviridae 50,585（Simplexvirus 5,837 / Varicellovirus 5,397）、Mimiviridae 120,124、Protist 宿主 vOTU 185,751 居首。**两条路线在同一数据集上看到同一批 vertebrate/protist 污染信号**——写论文时这是现成的"污染可控性/可重复性"论据，但同时也要求我们对 Plant 之外的分析层保持与 Mifsud 同等保守。

## 6. 差异与各自局限

| 维度 | 他们的局限（=我们的相对优势） | 我们的局限（=他们的优势） |
|---|---|---|
| 覆盖面 | Mifsud 仅 RNA 病毒（polyA 偏倚+RdRp 框架），主动剔除 Phycodna/Mimivi；未做类病毒 | 我们同样继承 polyA 偏倚，但流水线覆盖 DNA 病毒与类病毒（viroid 4 种，CEVd 398 样本） |
| 深度 | 平均 1.99 Gbp/样本，浅；宿主去除可能误删 RT 病毒 reads | 我们做了完整宿主去除+定量+变异层；但 depth≥5× 高置信规则对低载量样本偏严（CEVd 104 检出→16 确认） |
| 可复现 | 人工策展为主，流水线难复现；无标准化阈值链 | 全流程参数化（10 工具、CDD 分层、5 层证据、Salmon/iVar），可复跑 |
| 新病毒确认 | 104 条已交 GenBank、系统发育定级（含新科建议） | 我们的 184 候选未逐条 Rdp 策展/提交——这是发表前必须补的活 |
| 进化分析 | 共进化（COPAP/Cophylogeny：65% 跨种传播）、30K MP 起源等演化叙事 | 我们暂无共进化分析层（可引他们的结论作为背景） |
| 定量/变异 | 无定量、无变异、无全长恢复率 | — |

## 7. 对论文写作的直接建议

1. **引言/讨论引 Mifsud 2022 为唯一同数据集前驱**：差异点写"全量 1,342 vs 公开 1,079、全谱（DNA/RNA/viroid）vs RNA-only、流水线可复现 vs 人工策展、定量+变异层 vs 仅检出现"。
2. **数字并列时固定口径**：603/960 物种（Mifsud）↔ 887/1,341 样本原始检出（我们）；104 curated species ↔ 184 contig-level candidates（勿混写）。
3. **Partitiviridae/Deltapartitivirus、Varicosavirus、CMV+satRNA、Lettuce chlorosis virus 等重叠分类群**是天然的独立验证素材，可在 Results/Discussion 显式互引。
4. **污染段落**引 Mifsud 的定量（flu A/adenovirus/HIV）与我们 humimdef1/Adenoviridae/Herpes 的对应计数，证明该数据集的污染景观已被两条路线刻画。
5. **PVD（3,353 物种）**可作参考数据库对照（我们用 ICTV MSL41）写入 Methods 或 Table S。
6. 发表前需补：新种候选的 RdRp 策展+系统发育定位+GenBank 提交（对齐 Mifsud 的物种级标准），把"184 候选"升级为可被 ICTV 对话的物种级清单。

## 附：五篇文献关键数字备查

- OneKP 本体：1,342 transcriptomes / 1,124 species（Nature 2019）
- Mifsud：1,079 libraries/960 spp；median 25,187,714 PE reads/库；median 82% host reads；41,256,176 contigs；603 spp 病毒阳性；+ssRNA 61%/RT-dsDNA 22%；104 novel RdRp；13 科宿主域扩至低等植物、4 科扩至藻类；共进化 cross-species transmission median 65%
- Vogrinec：79 种水生（1KP contigs）+14 种自测；此前文献仅 31 病毒/18 种水生植物；检出 TuYV、CMV、LCV、新 Sagittaria potyvirus、begomoviruses
- Higuita：114 SRA 数据集/13 BioProject/12 种；12 已知 + 6 新关联 + 6 推定新种
- PVD/Atlas：3,353 病毒种、9,010 病毒–宿主对、1,986 作物病毒；dsDNA 宿主域最广
