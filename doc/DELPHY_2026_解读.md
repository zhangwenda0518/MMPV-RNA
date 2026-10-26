# Delphy（Nature 2026）全文解读：对本项目「时间 + 地理分析」的可用性

> **原文**：Varilly P, Schifferli M, Yang K, et al. *Scalable near-real-time Bayesian phylogenetics for outbreaks with Delphy.* **Nature 658, 747–754 (15 October 2026)**. DOI `10.1038/s41586-026-11012-6`｜PMID 42749814｜Open access｜投稿 2025-04-01，接收 2026-08-10，online 2026-09-16。
> **本地文件**：`E:/谷歌下载/s41586-026-11012-6.pdf`（22 页）→ 已抽全文 `_refs/delphy/_nature2026.txt`（113 K 字符）。
> **代码**：`_refs/delphy/`（HEAD `5d5989f`，2026-10-01，MIT）｜**数据**：`_refs/delphy-2026-paper-data/`（tree:0 惰性克隆）。
> **第四来源**（二手，仅用于交叉核对）：微信公众号「基因组学研究前沿」罗林焘，2026-09-22，`mp.weixin.qq.com/s/yIewFJeMXRPjpX2Oo23o_Q`。该文提供了**版本演进对照**与**两处原文不一致**的线索（已核，见 §1.5 / §9）；凡二手文与我核实结果不同处，**以原文为准并注明**。
> 本文只回答一件事：**这条链里的哪些环节能用 Delphy，哪些不能用，判据是什么。**

---

## 0. 结论（读完全文后的修正版）

**Delphy 补的是「时间轴」的速度与独立验证，补不了「地理轴」——而且对本项目最关键的限制不是模型，是「日期精度」。**

1. ✅ **精确**，不是近似：数据似然被写成「匹配则 1，否则 0」（测序错误被忽略），全部信息压在 **genetic prior** 上；把 root 序列与突变史积分掉后**恰好等于 Felsenstein 树似然**（原文 Methods + Supplementary）。
2. ✅ **快 2–3 个数量级**（vs BEAST2）/ **约 2 个数量级**（vs BEAST X）；Ebola 数据集接近 **3 个数量级**；H5N1 实测 **42 分钟 vs BEAST X 跑了 21 天还没到一半**。⚠️ **但这个倍率随样本量增长**——757 / 1,610 条才拿得到；**81 / 174 条级别原文只写 "modestly more efficient"**（见 §1.5）。本项目 GCVA 135 / PSTVd 145 条正好落在后者。
3. ❌ **地理分析在本版里没有模型**。原文把 phylogeography 列为「同一数学基础可以支持」的**未来方向**，初始版本只有单分区 HKY + coalescent。所有地点着色都是 **MCC 上的朴素简约法**（`ties resolved arbitrarily`）。
4. ⚠️⚠️ **年精度日期是硬伤**。原文 Discussion 明写：「**year-long tip-date uncertainties are problematic for convergence…we recommend…excluding tips with date uncertainties above 1 year**」。H5N1 案例里，他们**直接把日期不到天的序列全部剔掉**（4,080 → 3,339），因为「runs with the larger set that includes those with uncertain dates exhibited severe convergence problems」。
5. ⚠️ **只有严格钟**（`Q = μν(ℓ)q_ab` 且 `This form implies a strict molecular clock`）。非严格钟被明确列为 future work。

→ 定位：**`clock` 阶段的分钟级独立定年 + `beast` 的双引擎一致性检验**。`phylogeo` / `geo` / `geo_paths` / `phymapr` 一律不动。

---

## 1. 论文到底做了什么（按原文结构）

### 1.1 核心：EMAT

- 带时间的二叉树上，**每条枝挂显式的、带时间戳的 mutation** `(a, ℓ, b, t)` 与 **missation** `(a, ℓ)`；只显式存一条参考序列，其余全是增量。
- 相对 UShER 的 MAT 三处扩展：① 节点与突变**都有时间**；② 允许同一枝同一站点多次突变；③ **显式编码缺失数据**。
- 关键代价/收益：`posterior ratio` 只需「与受影响突变数成正比」的少量局部运算，而 Felsenstein 侧要重算 attachment point 到 root 之间**每个节点、每个可变位点**的 partial likelihood（**几十万次运算**）。原文自陈代价是**统计效率更低**（相邻样本间需要更多步），但单步便宜得多。

### 1.2 三个关键机制

| 机制 | 作用 | 原文要点 |
|---|---|---|
| **mdSPR**（mutation-directed SPR） | 大范围重排的提议分布 | 用 Jukes–Cantor 近似提议 P′，Metropolis–Hastings 纠偏；提议分布「closely tracks」真实条件后验 → 远距离移动也能高接受率 |
| **N-pruning** | 缺失数据 | 不在 tip 附近插补缺失；对以缺失结尾的位点子树上做 Felsenstein 剪枝（贡献恒为 1）→ genetic prior 的**形式不变**，局部性与可并行性都保住 |
| **增广 coalescent** | 并行化 | 把 `k(t)` 换成「分区内 `k_p(t)` + 高斯耦合场」，积分掉耦合场后严格还原 Kingman coalescent |

> ⚠️ 原文明确：**mdSPR 的稳健性会随 tip 日期不确定性变差**——「we speculate that the entropy across many uncertain branch lengths overwhelms the constraints from tip sequences… **For now, Delphy is unsuitable for such datasets**」。这是 §0 第 4 条的技术根源。

### 1.3 精度证据：与 BEAST 的一致性

- SARS-CoV-2（757 条）：MCC 结构/时间一致；**clade correlation graph 覆盖 8,712 个 clade**，support 与 mean tMRCA 都落在对角线附近；少数离群点归因于「site-rate heterogeneity 处理不同」+「rugged tree landscapes」。
- Ebola（1,610 条，Skygrid）：MCC、clade 图、参数分布均一致；统计效率**接近 BEAST X 的 1,000 倍**。
- 唯一系统性差异：**α（位点异质）**。Delphy 用**连续 Gamma**（K=∞），BEAST 用 **K 类离散近似**；原文用 K=2 / 4 / 8 三档证明 α 随 K 收敛到 Delphy 的结果。**关掉异质后无显著差异**。
- ML 路线（IQ-TREE 2 + TreeTime）：「captures the essential lineage structure, but infers tMRCAs of major lineages less accurately」；对**群体曲线**及其连带的 μ、tree height 更不可靠。

### 1.4 速度实测（可直接引用的数字）

| 数据集 | Delphy | 对手 | 倍率 |
|---|---|---|---|
| SARS-CoV-2（757） | — | BEAST2 | **2–3 个数量级** |
| SARS-CoV-2（757） | — | BEAST X | **约 2 个数量级** |
| Ebola（1,610，Skygrid） | — | BEAST X | **接近 3 个数量级** |
| H5N1（3,339） | **42 min** | BEAST X **跑 21 天被中止**；IQ-TREE 数小时 | — |
| sims（100,000） | **14.80 h，ESS 248**（192 线程 / 96 vCPU） | 此前最大贝叶斯运行约 40,000 条、固定拓扑、约 3 周 | — |
| **小数据集（81 Ebola / 174 Zika）** | — | BEAST X / BEAST2 | **「较为有限」**（正式版措辞） |

**sims 的完整规模曲线（可外推本项目）**：

| 序列数 | 运行时间 | 后验 ESS | 线程 |
|---|---|---|---|
| 100 | 0.06 h | 1,234 | 5 |
| 1,000 | 0.28 h | 1,934 | 50 |
| 10,000 | 1.48 h | 694 | 192 |
| 100,000 | **14.80 h** | **248** | 192 |

> ⚠️⚠️ **两条必须遵守的引用纪律**（微信解读文专门点了这两条，论文本身没这么强调）：
> 1. **不能把「约 1,000 倍」当通用加速率**。正式版把 81 条 Ebola / 174 条 Zika 这类小数据集的优势从预印版的「约 10 倍」**改写为「较为有限」**。本项目 GCVA 135 条 / PSTVd 145 条正好落在这个区间——**指望 2–3 个数量级是错的，指望「几分钟 vs 几十小时」也要打问号**。
> 2. **H5N1 的「42 min vs 21 天」不是两个完整运行的耗时比**。BEAST X 是**跑到一半被中止**的，原文写 "interrupted before reaching halfway after 21 days"。引用时必须带这个限定。

> ⚠️ **Zika 是反例**：采样密度最低（约 1,200 个突变、**平均 4 个突变/枝**），Delphy 仍然精确但**变慢**，拓扑 ESS/min/sequence 比其他基准低约 **7 倍**。原文 Discussion 直接说：对「thinly sampled data spanning evolutionary timescales」，显式突变法**可能收益很小，甚至为负**。

### 1.5 版本演进：预印版 → 正式版（读二手材料前必看）

这篇工作的口径在四个版本里改过：白皮书 v1（2024-07-01）→ 白皮书 v1.1（2024-09-19）→ 预印版 → 正式版（2026-09-16 online）。**网上流传的中文解读多数基于预印版，照抄会误判加速倍率。**

| 项目 | 早期版本（预印版 / 白皮书） | 正式版（Nature 2026） | 本地核实 |
|---|---|---|---|
| 小数据集（81 Ebola / 174 Zika） | 优势「约 10 倍」 | 「**only modestly more efficient** than BEAST 2 or BEAST X」 | ✅ `_nature2026.txt:532` |
| Skygrid | 白皮书 v1.1 **全文 0 次提及** | **已实现**：Ebola 1,610 与 H5N1 3,339 两个基准都跑了 Skygrid | ✅ 白皮书 grep 无命中；`_nature2026.txt:2312` |
| 基准数据集 | SARS-CoV-2 757 / Ebola / Zika 174 | **新增 H5N1 3,339**（4,080 中筛出天级日期） | ✅ `_nature2026.txt:2297` |
| 校准研究 | 无 | **新增 WCSS**：200 重复 × 20 个研究（`wcss/` 已下载） | ✅ `_refs/delphy-2026-paper-data/wcss/` |
| 非严格钟 / 地理模型 | future work | **仍是 future work** | ✅ 结论未变 |

> ⚠️ 上表第 1 行的「约 10 倍」出自微信解读文对预印版的对照，**我在本地文件中无法直接核实**（白皮书正文不含加速倍率表）。可信的是**方向**：正式版对小数据集明确降级为 "modestly"，并在同页给出原因——「Delphy's advantage over Felsenstein pruning **grows with dataset size and sampling density**」（`_nature2026.txt:525`）。
> → **对本项目的直接后果**：GCVA 135 条 / PSTVd 145 条与 174 条 Zika 同量级，**必须按 "modestly" 预期**，即「比 BEAST 快一些」，而不是「快两个数量级」。

> 💡 **白皮书 v1.1 给出了机理**（`_whitepaper.txt:295-299`，正式版没这么直白）：小数据集慢不是算法问题，是**并行度上不去**——
> 「the Delphy runs for the **81-sequence Ebola** dataset can only sustain a load of **just over 1 vCPU** (resulting in ∼0.5 million steps per second), whereas the SARS-CoV-2 runs can keep **3–4 vCPUs** in use (∼2.5 million steps per second). **The implementation of the serial part of the MCMC in Delphy is currently particularly inefficient, and this impacts the Ebola runs most seriously.**」
> → **工程含义**：在 246 上开 96 线程对 135 条序列**没有意义**，吞吐大约卡在 1 个核；`--v0-threads` 调再大也不会有线性收益。真要提速只能加长链或换指数增长模型省维度。

---

## 2. 模型面板（Methods 原文口径，比 CLI help 更权威）

| 维度 | 原文表述 | 备注 |
|---|---|---|
| 数据似然 | 「**1 if the EMAT is consistent, 0 otherwise**」 | **不建模测序错误**；tip 只能是 A/C/G/T 或全 N |
| 替换模型 | 单分区 **HKY**；`Q = μν(ℓ)q_ab` | 部分简并（Y = C 或 T）**必须当缺失处理** |
| 分子钟 | 「**This form implies a strict molecular clock**」 | 非严格钟 = future work |
| 位点异质 | `ν(ℓ) ~ Gamma(α, α)`，`α ~ Expo(1)`；**连续**，非 K 类离散 | 因此 α 与 BEAST 不可比 |
| 根序列 | `π_a ~ Dir(1,1,1,1)` | |
| κ | `log κ ~ N(1, 1.25²)` | |
| tip 时间 | 「a priori uniformly distributed between their minimum and maximum times」 | 支持 `YYYY` / `YYYY-MM` / `YYYY-MM-DD` / **`YYYY-MM-DD/YYYY-MM-DD` 区间** |
| 祖先先验 ① | 指数增长 `N(t) = n0 e^{g(t−t0)}`；n0 默认 improper 1/x；g ~ Laplace(mean 0.001/yr, scale 30.701135/yr) | |
| 祖先先验 ② | **Skygrid**：staircase（标准）或 **log-linear（Delphy 特有）**；默认**固定** τ 而不推断；加「1 天以下 coalescence 时间」惩罚 | 原文：比 Gill 2012 原版更稳 |
| 并行 | ~400 个时间格离散 coalescent | 100k 数据集用了 10,000 格 |

**原文对离散化的判据（重要）**：当采样强度高到 `k(t)/N(t)` 与全基因组突变率**同量级**时，400 格会出现可察觉的数值伪影（Extended Data Fig. 4）。→ 本项目若样本密集，需要调大 `--v0-target-coal-prior-cells`。

---

## 3. 能力边界：与本项目现状逐条对照

| 本项目阶段 | 现用 | Delphy | 判断 |
|---|---|---|---|
| `clock`（RTT / DRT） | `temporal_signal.py` + 20 次随机化 | 可给贝叶斯 μ 区间 | ✅ **第二判据** |
| `beast`（全长定年） | BEAST1 `ucln` + `skyline`，5 链 × 5e6 | strict + skygrid/exponential | ⚠️ **只能交叉验证** |
| `gene_dating` | 逐基因 BEAST1 | CLI 无分区 | ❌ |
| `treetime` | TreeTime 0.11.x | 不同性质 | ❌（但原文证明 ML 的 tMRCA 不可靠，可作为**反向论据**） |
| `phylogeo` | BEAST1 + CTMC + BSSVS | **无地理模型** | ❌ **不可替代** |
| `geo` / `geo_paths` | RRT / TempMig / 传播路径 | 产物无 `location` 注释 | ❌ 吃不了 |
| `phymapr` | TreeTime 的 `timetree.nexus` | 可喂 Delphy MCC | ⚠️ 需实测解析 |
| MCC | TreeAnnotator | `delphy-mcc`（毫秒） | ⚠️ **只作自检** |

### 关于地理：必须准确表述

原文有两处容易误读，必须分清：

1. **「Delphy 的数学基础可以支持 phylogeography」**（"A flexible computational core" 一节）——这是**能力声明**，不是已实现功能。同段紧接着写：初始版本只有单分区 HKY + coalescent，「Future updates will relax these restrictions」。
2. **mpox / H5N1 图上的地理着色**——原文图注写明是 **"naive parsimony on the MCC"**，方法里写 **"ties resolved arbitrarily"**。作者自己给了限度：

   > 对 mpox，「Delphy's MCC parsimony reproduces the ref. 4 central conclusion at regional resolution… **However, in contrast to their phylogeographic reconstruction, Delphy's simple parsimony cannot establish a finer geographical state-level assignment to Rivers with high certainty.**」

3. **源码层面的确证（2026-10-17 复核，比论文更硬）**——直接对 `_refs/delphy/` 全仓检索：

   | 检索 | 结果 |
   |---|---|
   | `grep -rn -i "phylogeo\|location model\|discrete trait\|migration" core/ tools/` | **0 命中** |
   | `core/cmdline.cpp` 里全部 `location` 命中 | 只有 3 处，全是 Laplace 先验的 **location 参数 μ**，与地理位置无关 |
   | 全仓 `fitch` 命中 | `core/utree.cpp`（**核苷酸** Fitch set，用于建树/SPR）、`tools/beast_trees_to_dphy.cpp`（BEAST 树转 EMAT 时重建**祖先碱基**）、tests、plans |
   | 仓内 `.js/.ts/.html` | **0 个文件**——Web UI 未开源 |

   → 结论：**没有地理模型，也没有任何地理相关命令行开关**。论文图上那层地点着色由 **Web UI 的 Customize Panel** 做（Extended Data Fig. 1D：「A naive Fitch parsimony algorithm propagates metadata to interior MCC nodes」），**不在 CLI 产物里**。若走 CLI 做管线，`.dphy` / BEAST 格式输出**不带 location 注释**，Fitch 传播得自己实现（算法本身很简单，但没有后验含义）。

**所以**：不能写「Delphy 也能做系统地理学」。能写的是「Delphy 的 MCC 简约法可在区域层面复现地理结论，但**无法给出 state 级的高置信分配**」。这恰好是 BSSVS 存在的理由。

**三层能力分级（判断时按这个对号）**：

| 层 | 能力 | Delphy | 产物性质 |
|---|---|---|---|
| **L0** | MCC 上 Fitch 简约传播 + 着色 | ✅ 仅 UI 层 | **点估计**，`ties resolved arbitrarily`，无置信度 |
| **L1** | 离散性状 CTMC（速率矩阵 `Q`、平衡频率） | ❌ 无实现 | — |
| **L2** | BSSVS（迁移次数、BF）+ 路径后验 + 祖先地点概率 | ❌ 无实现 | — |

本项目 `phylogeo` / `geo` / `geo_paths` / `phymapr` 要的产物**全部在 L1/L2**——这不是「效果差一点」，是**输入结构上就不存在**（`geo_analysis.py::parse_mcc_states()` 读的是 MCC 上的 location state，Delphy MCC 里没有）。

---

## 4. 对本项目数据的适配度评估（按论文自己的判据）

| 判据（原文） | 本项目情况 | 结论 |
|---|---|---|
| **突变稀疏度**（1–2 突变/枝最佳；Zika 4 突变/枝已掉 7 倍） | GCVA 13 kb / 135 条；需实测「每个 tip 与其最近邻的差异位数」 | ⚠️ **先算再决定**；若 >4，加速比会显著缩水 |
| **日期精度**（>1 年建议剔除；年精度本身有问题） | 枸杞样本多为**年份**精度 | ⚠️⚠️ **最大障碍** |
| **分子钟** | 需要 ucln | ❌ 只能做敏感性分析 |
| **群体模型** | 多年定点采集，非暴发 | ⚠️ 必须显式换 `skygrid`，默认的指数增长会强行拟合 |
| **缺失数据** | 有 N/gap | ✅ N-pruning 正是为此设计 |
| **简并碱基** | 若比对含 IUPAC | ⚠️ 会被当缺失，需先统计占比 |

> **实操建议**：先用一天时间做「可行性三算」——① 每条序列与最近邻的差异位数中位数；② 日期精度分布（天/月/年各多少）；③ 比对中 IUPAC 简并位点占比。三项都过关再谈接入。

---

## 5. 五个可落地用途

### P0 — 定年前 gate（分钟级）

BEAST 5 链 × 5e6 步 × 48 h 超时，跑完才知道 μ / rootHeight 是否离谱。Delphy 几分钟给后验。

**换算**：原文明确「**10 Delphy steps ≈ 1 BEAST step**」（Methods，用于给 BEAST 设步数）。反过来：想对齐 5e6 步 BEAST，Delphy 跑 ~5e7 步即可，不必照搬。

### P0 — 双引擎一致性检验

同一模型下二者统计等价（原文：「samples the same tree distribution… once explicit mutation events are projected out」），所以 `treeHeight` / `clockRate` 的 95% HPD 重叠 = 一条强独立证据。

⚠️ **写法**：比的是「**严格钟假设下两者一致**」。`ucld.mean` vs Delphy 的 μ 不同名不同义，Methods 里要写成**模型一致性检查**，不是模型选择。

### P1 — WCSS 校准模板（`wcss/` 已下载）

200 重复 / 200–500 tip / L=30,000 / 三条判据（95% HPD 覆盖率、rank 均匀、clade 概率良校准）。

study 11（tip 日期不确定）实测：186 个有效重复，`rootHeight` 覆盖 **0.968**、`mu` **0.968**、`alpha` **0.946**，全在二项 95% 区间内；tip 日期本身月精度 94.05%（略低于下界 94.42%）、年精度 94.14%。

> ⚠️ 200 个重复里 **14 个因低 ESS 被剔除**（`excluded_replicates.tsv`）。先剔再算覆盖率会**偏乐观**，引用时必须说明。

### P1 — 树空间精度指标

`tree_ess/`（Rust，37 文件已下载）：`calc-tree-ess`（Frechet correlation ESS，Magee et al. 2024）、`compare-clades`（clade correlation graph）、`calc-clade-coverage`。
原文用它们对 **8,712 个 clade**（SARS-CoV-2）/ **43,845 个 clade**（H5N1）做了比对。**本项目现在只报参数 ESS，回答不了「拓扑收敛没有」。**

### P2 — 快速 MCC 自检

`delphy-mcc <in.trees> <out.nexus>`，毫秒级；原文说它能处理 100,000 条的模拟，「TreeAnnotator 在大基准上撞到运行时与内存上限」。

⚠️ 三个硬事实：
1. **burn-in 写死 30%**（原文 Methods 明说），本项目 `MERGE_BURNIN_FRAC = 0.10` → 喂之前先按 10% 裁好。
2. **等价于 TreeAnnotator 的 `--heights ca`**：节点时间是**全体后验树**上下游 tip 的 mean tMRCA，不是只在单系树里取。
3. **不输出 HPD、不输出 location 注释** → `geo_analysis.parse_mcc_states` 会崩。

### ❌ 明确不做 — 用 Delphy 做地理

绕法（固定 Delphy 树 + 单独离散 ASR）产出的不是 BSSVS，与本项目 `phylogeo` 的 BF 检验**口径不同**，并排即错。原文自己也没有这么做：mpox 的 BEAST 对照实验里，他们**主动移除了 phylogeography** 才做比较（"we modified the original BEAST X XML… and removed phylogeography, spillover detection and pre/post-spillover partitioning"）。

**如果只是要一张「扩散示意」**（不写进统计结论）：Delphy MCC + 自己写 20 行 Fitch 传播（或直接用 `phangorn::parsimony` / baltic）可以出区域级着色图。但它给的是 **单棵 MCC 上的点估计**，没有迁移速率、没有路径概率、没有 BF，**不能喂 `parse_mcc_states()` / `geo_paths`**，也**不能与 BEAST 的 BSSVS 结果并排比较**。定位只能是「配图」或「定性 sanity check」。

**真正值得做的地理相关用法**（交叉验证，非替代）：把 Delphy 的**时间标尺**与 BEAST 的**地理模型**组合成敏感性分析——若 BEAST 的 CTMC 区域级结论与 Delphy 简约法结论一致，说明地理推断不依赖具体树形；不一致的节点恰是模型在起作用的地方。⚠️ 但这属于**讨论段**能写的话，不能当新分析出图。

---

## 6. 命令模板

### 6.1 部署（246）

```bash
# 预编译二进制（Linux x86_64 / ARM64）或 Docker
# ghcr.io/broadinstitute/delphy:latest
```
⚠️ 网页版比原生命令行**慢 2–3 倍**（WebAssembly，原文 Methods），**论文数字都来自原生命令行**。

### 6.2 输入（与现有 `phylo_dates_135.csv` 不同）

ID 行用 `|` 分隔：**第一段 = 短 ID（与 metadata 的 `id`/`accession` 列对应）**，**最后一段 = 日期**，中间字段忽略。

```
>GCVA_001|China|2019-05-14
>GCVA_002|China|2019-05
>GCVA_003|China|2019
>GCVA_004|China|2019-05-10/2019-05-20
```

- 日期四档：天 / 月 / 年 / **区间**。
- **ID 里不能有空格**（整行都当 ID，与 BEAST/IQ-TREE 一致）。
- metadata：CSV/TSV，含 `id` 或 `accession` 列；缺失值写空、`-`、`noknown` 或 `none`。

### 6.3 跑法

```bash
# ① 快速预检（严格钟 + 指数增长，对应 BEAST 的 exponential）
delphy --v0-in-fasta gcva_delphy.fasta \
  --v0-pop-model exponential \
  --v0-steps 50000000 \
  --v0-threads 8 --v0-seed 42 \
  --v0-out-log-file gcva.log --v0-out-trees-file gcva.trees \
  --v0-out-delphy-file gcva.dphy

# ② 对齐本项目 skyline 口径：skygrid，24 个区间，cutoff 1 年
delphy --v0-in-fasta gcva_delphy.fasta \
  --v0-pop-model skygrid \
  --v0-skygrid-type staircase \
  --v0-skygrid-num-parameters 24 \
  --v0-skygrid-cutoff 1.0 \
  --v0-target-coal-prior-cells 400 \
  --v0-out-log-file gcva_skygrid.log --v0-out-trees-file gcva_skygrid.trees

# ③ MCC（先按 10% 裁树，再喂它——它内部还会再切 30%）
delphy-mcc gcva_10pct.trees gcva_mcc.nexus

# ④ 导出等价 BEAST XML（用于确认「比的是什么模型」）
delphy --v0-in-fasta gcva_delphy.fasta --v0-out-beast-xml gcva.xml \
  --v0-out-beast-version 2.7.7 -n
```

### 6.4 与 BEAST 对齐时的三个换算

| 项目 | Delphy | BEAST1（本项目） | 对齐方式 |
|---|---|---|---|
| 步数 | ~100,000/tip 默认 | 5e6 / 15e6 | **10:1**（原文启发式） |
| 群体模型 | `skygrid`（M+1 个 knot） | `skyline`（`skyline.popSize` dimension=5） | 区间数需换算，**不是 1:1** |
| 时间单位 | log 输出用**线性年** | `decimal_year` | 自写解析器**禁用 `days/365`** |

---

## 7. 踩坑清单（13 条，每条都会静默出错）

1. **⚠️⚠️ 年精度日期 → 收敛问题。** 原文建议**剔除不确定性 > 1 年的 tip**；H5N1 里 4,080 条中只保留 3,339 条有**天级**日期的。技术根源是 mdSPR 在大量不确定枝长下失效。
2. **⚠️⚠️ 稀疏采样可能「越用越慢」。** Zika（4 突变/枝）拓扑 ESS/min/seq 低 7 倍；原文说对跨演化尺度的稀疏数据「may offer smaller efficiency gains, or even a loss」。
3. **⚠️⚠️ `delphy-mcc` burn-in 硬编码 30%，且等价 `--heights ca`。** 与本项目 10% 口径冲突；无 HPD、无 location。
4. **⚠️⚠️ 连续 Gamma ≠ 离散 K 类 → α 不可比。** 原文专门用 K=2/4/8 证明 α 随 K 收敛到 Delphy 值；关掉异质才无差异。
5. **⚠️⚠️ 时间单位。** `core/beasty_output.cpp` 顶部注释：早期用 `/365` 输出导致后验日历年系统性偏低 1–2 天，被 WCSS 抓出；正解是 `to_linear_year()`。
6. **⚠️ Skygrid 参数差一位。** `--v0-skygrid-num-parameters` = BEAUTi 的 "Number of parameters" = Gill 的 `M`+1。写 24 = **25 个 knot**。
7. **⚠️ 简并碱基当缺失。** 「partially ambiguous states (such as Y = C or T) must be treated as missing data」。若比对里有 IUPAC，先统计占比，否则等于人为删数据。
8. **⚠️ 不建模测序错误。** 数据似然是 0/1。低质量碱基会被当成真突变，BEAST 至少还有误差项可吸收。
9. **⚠️ coalescent 离散化的密度判据。** `k(t)/N(t)` 与全基因组突变率同量级时 400 格会出伪影（Extended Data Fig. 4）；需调大 `--v0-target-coal-prior-cells`。
10. **⚠️ 只有严格钟。** `ucln` 无法复现；非严格钟是 future work。
11. **⚠️ `--v0-` 前缀 = 不稳定语法**（论文期间改名多轮）。钉版本 + 钉命令行。
12. **⚠️ 版本错配。** 论文基准用 **v1.1.4 (build 2044)**；`sims` / `gisaid-week-by-week` 用 **v1.0（无 Skygrid）**。引用数字时必须对应版本。
13. **⚠️ `--v0-paranoid` 很慢。** 只在怀疑树完整性时开。

---

## 8. 建议的接入顺序

| 步 | 动作 | 验收判据 |
|---|---|---|
| 0 | **可行性三算**：最近邻差异位数中位数 / 日期精度分布 / IUPAC 占比 | 差异位数 ≤2 才指望 2–3 个数量级；年精度占比高则先补日期 |
| 1 | 246 装二进制；写 `utils/delphy_io.py`（FASTA+CSV → delphy FASTA） | 条数、序列 md5 与源一致；ID 与 `seq_ids.py` 映射一致 |
| 2 | GCVA / PSTVd 各跑一次 `exponential`，比 `treeHeight` / `clockRate` | 两个 95% HPD 重叠；不重叠先查日期与掩蔽 |
| 3 | `skygrid`（区间数对齐 skyline 维度），比群体曲线趋势 | 趋势一致即可，**绝对 N(t) 不比** |
| 4 | `tree_ess/` 三个指标接进 `merge_results.py` | 树空间 ESS 与参数 ESS 趋势一致 |
| 5 | 缩小版 WCSS（20–50 组） | `rootHeight` / `μ` 的 95% HPD 覆盖率 90–98% |
| 不做 | 任何「Delphy 做地理」的尝试 | — |

**落盘**：脚本 `virome_phylo_pipeline/utils/delphy_io.py` + 新 stage `delphy`（排在 `clock` 之后、`beast` 之前）。
⚠️ 新 stage 必须遵守 **`--io-layout legacy`**（硬规则 ①），否则会静默写到 `06_time/`。

### 8.1 成本视角：接 Delphy 的边际收益到底有多大

先算清楚本管线的钱花在哪（数据来自 `config.yaml` 与本数据 `datasets.yaml`，2026-10-18 核）：

| 环节 | 配置 | 相对成本 | 说明 |
|---|---|---|---|
| `clock`（RTT / DRT） | 20 次随机化 | **≈0** | 分钟级 |
| **`beast`（时间分析）** | `chain: 5,000,000` × `chains: 5` = 25 M states | **1×（基准）** | GCVA 145 条 × 16.3 kb |
| `phylogeo`（地理推断） | `chain: 10,000,000` × `chains: 5` = 50 M states | **≈2.5×** | 8 个地点状态 → `n_rates = 8×7 = 56`（`beast1_bridge.py:505`） |
| `geo` / `geo_paths` | MCC 解析 + 200 次置换 + 36 帧 GIF | **≈0** | 秒到分钟级 |
| **RRT 树注释法** | 需 **20 棵随机化 MCC** = 20 次完整 `phylogeo` 重跑 | **≈50×** | `config.yaml:83-88` 自己写明「相当于 ~21 次完整 BEAST 重跑」 |

**两个反直觉的点**：

1. **地理模型的「单步开销」其实很小（<5%）**。trait 似然是 `(2n−1) × L²` ≈ 289×64 ≈ 1.8 万次运算；而 GCVA 的树似然是 `(2n−1) × 位点模式数 × 16` ≈ 数百万次。所以 `phylogeo` 慢**不是因为每步更慢**，而是因为 ①链长写成了 2 倍，②BSSVS 的 56 个速率 + 56 个指示子让后验更多峰、ESS 更低，**这 2 倍是收敛需要的，不是浪费**。
2. **真正的时间黑洞是 RRT，不是 phylogeo**。`geo` 那一堆下游分析（pathways / dispersion / 速度场）都是秒级；但 RRT 树法要 20 次「打乱地点标签 → 重跑 phylogeo」，直接乘 20。

→ **对 Delphy 的结论**：它能压缩的只有上表第 2 行的 **1×**，而且按 §1.5 对 145 条这个量级只 "modestly"。相对总成本（1 + 2.5 + 50），**边际收益 ≈ 2%**。若 RRT 是必做项，接 Delphy 省下的时间在总量里几乎看不见。
→ 反过来说：**若只跑「时间分析 + 地理推断」、不做 RRT**，地理部分占 2.5/(1+2.5) ≈ **71%**，此时把 `phylogeo` 的链长压下来（或先只跑 2 链探路）比接 Delphy 更划算。

---

## 9. 引用前必读：数字口径与两处「看起来矛盾」的地方

**凡是把 Delphy 的数字写进论文/报告，先过这一节。**

### 9.1 `3,339` 与 `3,383` 不是矛盾，是两套分母

微信解读文把这里标为「原文不一致、PDF 未解释两者的对应关系」。我核了 Methods 原文，**对应关系其实可推**：

| 数字 | 含义 | 出处 |
|---|---|---|
| **4,080** | 全部候选：11,541 个 unique SRR → 剔撤回 → cattle host + B3.13 基因型 | `_nature2026.txt:2297` |
| **3,383** | 其中**有 GenBank 登录号**的全集（含只有年精度的）——**地理分析的分母** | `:2309` |
| **3,339** | 3,383 中**日期精确到天**的子集——**定年的分母** | `:2304` |
| **3,194** | 3,383 中 `geo_loc_name` 能解析到美国州级的——**地理分析的实际可用数** | `:2309` |

- 差额 **3,383 − 3,339 = 44 条**，正是「有登录号但日期只到年」的那批。
- 原文用 "the larger set that includes those with uncertain dates" 指代 3,383，**但没写出这个减法**——所以严格说不是「矛盾」，是**没显式写等式**。
- ⚠️ **引用纪律**：写「H5N1 3,339 条」必须带定语「**日期精确到天**」；写地理分析用 **3,383 / 3,194**。两个数字混用会被审稿人抓。

### 9.2 图与正文的数值版本错配

微信文指出：图 3 存在「图片内嵌旧数值 / 旧图例」与正式版正文数值不一致的情况。**我无法从抽出的纯文本验证图内嵌字**（图是位图），因此采取保守规则：

> **一切数字以正文 + Methods 为准，图注 / 图内嵌字不作引用来源。**

### 9.3 三条容易被二手解读带偏的口径

1. **「快 1,000 倍」只对 Ebola 1,610 / SARS-CoV-2 757 成立**，且原文限定为「statistical efficiency」（ESS 除以墙钟时间），不是端到端墙钟时间。小数据集是 "modestly"（§1.5）。
2. **H5N1 的「42 min vs 21 天」不是两个完整运行之比**——BEAST X 是 "interrupted before reaching halfway"（§1.4）。
3. **α（位点异质）与 BEAST 的 α 不可比**（连续 Gamma vs K 类离散，§2）。要比就先把 Delphy 的位点异质关掉。

### 9.4 两个方法学细节，值得直接抄进本项目

| 细节 | 原文做法 | 本项目可借用处 |
|---|---|---|
| **自动 burn-in 判据** | 取轨迹**后半段**的均值与标准差 → 找最后一次波动**超过 5 sd** 的位置 → 再找其后**首次落回 2 sd 以内**的位置；对 `log-posterior` / `μ` / `total evolutionary time` 三条轨迹取**最晚**的那个；用户不覆盖则持续更新 | 可替代 `run_phylogeo.py::MERGE_BURNIN_FRAC = 0.10` 的固定 10%。⚠️ **不冲突**：论文基准另行统一用 **30%**（"verified to suffice by visual inspection… higher than necessary"），30% 是「安全过量」，自动判据是「够用即可」 |
| **Minimum ESS 的口径** | 汇总所有观测量 ESS 取最小，但**排除 mutation count**——它只取少数离散值，相关高斯随机游走假设不成立（Extended Data Fig. 1 注） | 本项目报 ESS 时也应排除同类「离散 / 计数型」观测量，否则最小 ESS 会被它拖死 |

> 顺带一条可直接引用的「为什么该自动化」：原文说早期用户「requiring manual cut-off selection led to either **needless friction or no cut-off at all, biasing the results**」（`_nature2026.txt:2069`）。

### 9.5 mdSPR 的 99% / 1% 扫描比例

原文 Methods：提议新拓扑 P′ 时，「**scanning points up to 1 mutation away from P in 99% of cases, the whole tree (partition) in 1% of cases**」（`_nature2026.txt:1996`）。

→ 解读：**默认是极局部搜索**，全局跳跃只占 1%。这正好解释两件事：① 为什么「稀疏采样 / 长枝多」会变慢（Zika 4 突变/枝，§1.4）——局部邻域里找不到好提议，1% 的全局跳跃不够用；② 为什么 `--v0-paranoid` 这类「更彻底」的开关代价很高。

---

## 10. 可直接引用的原文句子

> "Delphy, an **exact reformulation** of Bayesian phylogenetics designed to transform its speed, scalability and accessibility while retaining Bayesian state-of-the-art accuracy."

> "we demonstrate state-of-the-art accuracy with up to **2–3 orders of magnitude improvements in speed**."

> "a simulated dataset of **100,000 sequences can be analysed within a day**."

> "**For now, Delphy is unsuitable for such datasets**"（指大量 tip 日期不确定性大的数据）

> "**year-long tip-date uncertainties are problematic for convergence**… we recommend corroborating any conclusion that depends on exact sequence placement through multiple independent runs and **excluding tips with date uncertainties above 1 year**."

> "Delphy's simple parsimony **cannot establish a finer geographical state-level assignment** to Rivers with high certainty."（对比 ref. 4 的 phylogeographic reconstruction）

> "For datasets that are far from parsimonious, locally or globally, Delphy **degrades predictably: it remains exact but slows down**, rather than silently introducing artifacts."

> "for the 81-sequence Ebola and 174-sequence Zika datasets, Delphy is **only modestly more efficient** than BEAST 2 or BEAST X."（小数据集口径，§1.5）

> "Delphy's advantage over Felsenstein pruning **grows with dataset size and sampling density**."（同一页，解释了上一条）

> "**scanning points up to 1 mutation away from P in 99% of cases, the whole tree (partition) in 1% of cases**."（mdSPR 的局部性，§9.5）

> "experience with early users showed that requiring manual cut-off selection led to either **needless friction or no cut-off at all, biasing the results**."（自动 burn-in 的理由，§9.4）

---

## 11. 文件索引

```
_refs/delphy/
├── _nature2026.txt               # ★ 论文全文（22 页，已抽文本）
├── delphywp.pdf / _whitepaper.txt# 白皮书 v1.1（2024-09-19，28 页）——注意比论文旧
├── doc/dphy_file_format.md       # .dphy FlatBuffers 字段表
├── core/                         # 引擎（41 个 .h）
│   ├── cmdline.cpp               # ★ 全部 CLI 选项与默认值
│   ├── dates.cpp/.h              # 2020-01-01=day0；parse_iso_date/month/year
│   ├── io.cpp                    # FASTA 解析；missation 区间语法
│   ├── beasty_output.cpp         # ★ BEAST 格式输出；顶部有 /365 偏 2 天的注释
│   ├── coalescent.cpp / scalable_coalescent.cpp / very_scalable_coalescent.cpp
│   ├── pop_model.cpp/.h          # exponential / skygrid
│   └── mcc_tree.cpp/.h           # 指纹 MCC
├── tools/                        # delphy / delphy_mcc / delphy_ui / delphy_wasm
├── plans/                        # 30 篇设计文档（先验推导细节）
└── tests/                        # 28 个 gtest

_refs/delphy-2026-paper-data/     # tree:0 惰性克隆
├── README.md                     # 11 个数据集 + 软件版本 + AWS 配置
├── tree_ess/                     # Rust：calc-tree-ess / compare-clades / calc-clade-coverage
└── wcss/                         # 200 重复 × 20 个研究
    ├── plans/*.md
    └── 11_tip_date_uncertainty/  # ★ 与「年份精度日期」最相关（含实测覆盖率 TSV）
```

**数据/代码存档**：数据 Zenodo `10.5281/zenodo.17899935`（论文口径；README 另给 `19819574`）；代码 `10.5281/zenodo.11922813`；Sapling `10.5281/zenodo.20717608`。
