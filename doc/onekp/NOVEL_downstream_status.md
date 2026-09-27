# OneKP 新种候选（184）画像与下游分析缺口 — 2026-09-14

## 一、数量口径

| 口径 | n | 来源 |
|---|---|---|
| **新种候选（主口径）** | **184** | final_judgement_table.tsv `category=新种候选`；证据分中位 75.3/100；KEEP 96 / REVIEW 88 |
| 新病毒证据标签（含 REVIEW/DROP） | 3,847 | rescue_evidence_scored.tsv novelty 列（blast 蛋白级 461 + 仅域证据 2,095 + blast 无域 1,291） |
| 远缘候选 | 843 | 同表 category=Distantly-related virus candidate |
| ⚠ 口径警示 | — | Mifsud 2022 的 104 为 RdRp 策展级物种（已交 GenBank）；本表为 contig 级候选，不可直接互比 |

## 二、184 候选画像

- **科**：Partitiviridae 72（39%）≫ Caulimoviridae 27、Rhabdoviridae 25、Betaflexiviridae 15、Geminiviridae 9、Secoviridae/Tymoviridae 各 5、Potyviridae/Amalgaviridae 各 4、NA 7
- **属**（52 个 NA 待补）：Deltapartitivirus 17、Solendovirus 14、Betanucleorhabdovirus 12、Carlavirus 10、Betapartitivirus 8、Alphapartitivirus 7、Badnavirus 6、Soymovirus/Welwivirus 各 5
- **长度**：中位 1,550 bp，最长 13,225 bp（多为片段）
- **流行度**：中位 2 样本；≥2 样本 107 个；单样本 77 个
- **内部不一致**：16 条 novelty 标 "Known virus(blast蛋白)" 却被归为新种候选（与审计 #17 同族的判定矛盾，需复核）

## 三、已完成的分析（发现线内）

1. 四分支 rescue 级联（定级对象即 rescue 后序列）
2. 五层证据打分（0.30 BLASTN + 0.30 BLASTX + 0.40 CDD）
3. ACVirus 独立分类验证（1,969 条池）+ SDT 矩阵（09b）
4. vOTU 聚类（进 20,054 植物集）+ 流行度表
5. class_KEEP.fasta 保留（96 条 KEEP 新候选在内）

## 四、五项下游缺口（全部未做）

| # | 缺口 | 现状 | 优先级 |
|---|---|---|---|
| 1 | **CheckV 完整度** | 184 条 aai_completeness/confidence **全部 NA** | ★★★ |
| 2 | **系统发育定级** | **基建已在**：ACVirus 已对 class_KEEP（2,322 条）按科建 **64 棵 ML 树**（prodigal 预测蛋白 → DIAMOND 对 acvirus_db 定科提参考 → MAFFT 比对 → IQ-TREE，.treefile 含 bootstrap）。缺口：96 KEEP 新候选已在树中但未做"独立分支=新种"的正式解读；88 条 REVIEW 候选与 52 个 NA-属候选未进树 | ★★★ |
| 3 | **变异/群体遗传** | 管道三只覆盖 71 种已知病毒（118 对）；新候选 0 | ★★（≥2 样本的 107 个值得做） |
| 4 | **全长恢复** | 05–07 阶段仅对 73 个已知 accession；新候选 0 进组装/提取/相似性 | ★★ |
| 5 | **GenBank 提交** | 0 条 | ★★★（发表硬要求） |
| — | 宿主物种归属 | 被 final_judgement `host_species=out` 失效问题阻塞 | ★★ |

## 五、RdRp-Catch 是否需要对延伸后结果重跑？——需要，且数字支撑

**已做的**：02a 阶段 RdRp-Catch 跑在**原始组装 contig**上（9 工具面板之一，原始阳性 55,838 条），用途是"识别"。

**缺口**：03a COBRA + 08 四分支 rescue（尤其 C 分支 9,212 条新 contig）产生的延伸序列**从未被 RdRp 层面重新分析**。延伸后 CDD 层的现状：

- 8,915 条植物证据序列中仅 **160 条**带 RdRp 域（pfam00946 系），其中 104 条有属标签；
- **184 个新候选中仅 35 条**有 RdRp 域证据 → 149 条的 RdRp 状态未知。

**重跑目的（不是重新识别，是定级）**：
1. 确认延伸后 RdRp 区 motif I–VII 是否完整 → 决定哪些候选可进系统发育树（缺口 #2 的前置）；
2. RdRp 最近邻分类 → 给 52 个 NA-属新候选 + 远缘候选补属级标签；
3. 产出 Mifsud 式策展与 GenBank 提交所需的 RdRp 比对集。

**建议范围**：无需重跑全目录（610,535 vOTU）；先对 **184 新候选 + 843 远缘候选（≈1,027 条）**重跑 RdRp-Catch，再视结果扩展到 8,915 条植物证据集。计算成本低（千条级 HMM 扫描）。

## 六、建议执行顺序

1. 对 1,027 条（新候选+远缘）重跑 RdRp-Catch → 补属标签 + RdRp 完整性清单
2. 对有完整 RdRp 的候选建树定级（对照 MSL41/RVDB 参考集）
3. 对定级候选跑 CheckV（当前全 NA）
4. ≥2 样本的 107 个做 consensus + reads 重比对的变异分析（复用管道三）
5. 全长恢复（复用 05–07 流程）→ GenBank 提交包
