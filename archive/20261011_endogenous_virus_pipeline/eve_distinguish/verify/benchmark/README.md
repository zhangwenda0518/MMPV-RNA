# benchmark/ — 策展基准集（判别器的灵敏度/召回）

`verify/` 那四个工具只能测**特异性**：端到端空模型把信号全破坏了，它证明"随机序列不会
触发"，证明不了"真病毒会被抓住"。这个目录补最后那一环——**标签已知**的序列。

## 四类标签与来源

| 类 | 含义 | n | 来源 |
|---|---|---|---|
| `A_extant_virus` | 完整基因组、编码完好、**未整合**——生物上不是 EVE | 486 | NCBI RefSeq 完整基因组：Caulimoviridae 七属 102 条（Caulimovirus 12 / Badnavirus 77 / Soymovirus 5 / Cavemovirus 3 / Tungrovirus 2 / Solendovirus 2 / Petuvirus 1）+ Geminiviridae 384 条 |
| `B1_eprv_clean` | 文献报道已整合的元件（**完整/大片段**档） | 123 | NCBI 查询 `endogenous[Title] AND Caulimoviridae`、`endogenous pararetrovirus`、`Musa balbisiana BAC`（含 BSGFV 整合等位 EPRV-7/9），加登录号 AP009325 / AP009326 |
| `B1b_eprv_decayed` | **退化型个体基因组拷贝（目标类别）** | 103 | 糖甜菜 (*Beta vulgaris*) EPRV 研究的个体基因组拷贝，Zenodo 3888270 的 S2–S5 五个比对文件，**去 gap 后** 869–8,066 bp，命名带染色体位置与方向 |
| `B2_eprv_activatable` | **可激活型** EPRV（可侵染）——标签有歧义 | 12 | TVCV（NC_003378）、PVCV、BSV 等 complete genome 记录（去重后 12 条独立序列） |
| `C_host_viral_domain` | 带病毒样结构域的**寄主基因** | 1023 | 植物端粒酶 TERT（含 RT 域，630）、植物 RNase H 域蛋白（392）、dUTPase（1），全部 RefSeq |

`B1b` 是**加进去专门为了分开两种难度**的：`B1` 大多来自 RefSeq complete genome 提交与 BAC 克隆，
属于"完好/可激活"档；而本方法的设计目标是**退化型**元件（"不全"不算 EVE 证据，是 v4 那个倒置的
核心）。把两种难度混在一个分母里，灵敏度就没有意义。糖甜菜的个体拷贝是已发表研究在真实基因组里
全基因组筛出的、带侧翼上下文的复制子——这才是"目标类别"。

第四类"转座子"不在这里下载：服务器上已有策展库（`Dfam3.9_plants.fa` 16,304 个植物重复
家族、`REXdb/GyDB2.hmm` 的 LTR 反转座子域模型），而且 pipeline 自己的 baits 就是
Caulifinder banks——直接用它们更贴近实际运行条件。

`B2` 单列、**不混进任何分母**：可激活 EPRV 既是整合元件又编码完好，判成"像病毒"在生物学
上并不算错，拿它当"必须判成 EVE"的分母是错的。

## 怎么跑

```bash
bash verify/benchmark/build_benchmark.sh --out BENCH        # 需要 efetch/esearch + 网络
bash verify/benchmark/run_benchmark.sh --bench BENCH --work WORK \
     --mod <模块目录> --evidence rescue_evidence_scored.tsv --threads 32
```

已验证可复现：2026-09-24 实测去重后 **1,747 条**（合并时丢掉 7 条跨查询重复）。

## 实测结果（2026-09-24，寄主通道关闭）

```
A_extant_virus (486)     review 79.2% / EVE_suspect 12.1% / virus_fragment_review 7.6% / ancient_EVE 1.0%
B1_eprv_clean  (123)     virus_fragment_review 30.1% / review 21.1% / EVE_suspect 17.1%
                         / ancient_EVE 10.6% / compositional_noise_review 8.9% / EVE_LTR_TE 8.9% / 其它 3.2%
B1b_eprv_decayed (103)   EVE_suspect 45.6% / virus_fragment_review 34.0%
                         / compositional_noise_review 14.6% / ancient_EVE 2.9% / 其它 2.9%
B2_eprv_activatable (12) EVE_suspect 83.3% / virus_fragment_review 16.7%
C_host_viral_domain(1023) review 99.2% / structure_intact_review 0.4% / EVE_LTR_TE 0.4%
```

**分两个严重度看**（混在一起会把风险讲大）：

| 类 | 会真的动数据 | 仅判成 EVE 侧、转人工复核（不删不移） |
|---|---|---|
| `A_extant_virus` | **5 / 486 = 1.0%**（`ancient_EVE` → MOVE_EVE） | 59 / 486 = 12.1%（`EVE_suspect` → REVIEW） |
| `C_host_viral_domain` | **4 / 1023 = 0.4%**（`EVE_LTR_TE` → MOVE_EVE） | 0 / 1023 = 0% |
| `B1_eprv_clean` | 灵敏度 **24 / 123 = 19.5%** | 其余 99 条留在病毒集里（漏清，但不删错） |
| `B1b_eprv_decayed` | 灵敏度 **3 / 103 = 2.9%** | 其余 100 条留在病毒集里 |
| `B2_eprv_activatable` | 12 条全部转人工，无一被移出 | — |

### 结论：特异性好，灵敏度差——而这是本基准集最重要的产出

- **特异性好**：现存病毒被误移 **1.0%**、寄主 RT 域基因被误移 **0.4%**。流水线真正会造成
  损失的那类错误比例很低。可信。
- **C 类那 4 条漏网机制可指名**：植物端粒酶/RNase H 基因打中 **TE baits**，走 `EVE_LTR_TE`
  被移出。TE 通道一个已知泄漏点，量小但明确。
- **灵敏度差，且目标类别上更差**：退化型 EPRV（`B1b`）只有 **2.9%** 被判 EVE 移出，
  低于完整型（19.5%）。**我此前的假设（"19.5% 低估了目标类别上的灵敏度"）被这一测直接
  推翻——方向相反。** 91% 的退化型拷贝落在 `EVE_suspect`(45.6%) 与
  `virus_fragment_review`(34.0%)，也就是**留在病毒集里不清**。
- **卡点是架构而不是退化**：`EVE_suspect` 的定义是"架构不全 + 分布式退化"。退化通道确实
  在开火（空模型也显示它 45% 触发），**真正卡住 MOVE_EVE 的是组件齐全度要求**——
  单一位点上要凑够 2–4 个规范组件才够格。这与敏感度扫描独立吻合：`COMP_AA_MIN`(8.2%)
  与 `MIN_PID`(11.5%) 是仅有的两个能动结论的阈值，而 `PROVIRUS_*` 已不参与判定。

**所以对外表述必须改口径**：这条流水线不是"能识别 EVE"，而是**"高特异性、低灵敏度"的保守
过滤器**——它很少误删（1.0% / 0.4%），但也基本不主动判出 EVE（2.9%）。说它"识别出了真 EVE"
是错的；说它"不删错、但会大量漏清"才对。

### 可执行的下一步（由本结果直接推出）

卡点是组件齐全度 ⇒ 把 `COMP_AA_MIN` / `MIN_PID` 放宽后重跑本基准集，看灵敏度是否上升、
以及**代价**（`A`/`C` 两类的假阳性升多少）。这正是一个 ROC 式的取舍曲线，而且现在有标签了，
可以真的画出来。这是本基准集最该被用去的地方。

## 必须一起读的三条限制

1. **`B1b` 的输入分布与方法的预期输入不完全匹配。** 糖甜菜拷贝是**基因组位点**（含寄主
   侧翼，869–8,066 bp），而本方法的设计输入是**病毒富集的组装 contig**。EPRV 片段在整条
   序列里可能只占一部分，这会压低命中率。所以 2.9% 要读成"在基因组窗口上"的灵敏度；
   在真正的病毒富集 contig 上会更高——但**高多少，本基准集测不出来**，需要 `B1b` 那类
   位点先做病毒富集（或改用只含元件的片段）。
2. **寄主通道关闭**：基准序列不在任何 OneKP 样本里，`host_*` 三类结论结构上不出现。
   本矩阵**不测寄主否决**（它只会移除/复核，不会创造 EVE 结论，所以也不影响上面那个灵敏度）。
3. **上游证据通道未测**：基准集 ID 不在证据表里，`tax_family`/CheckV/各条同源通道全为空。
   测的是 **S1×S2 的判定逻辑**。

另外 `A_extant_virus` 里 Geminiviridae 占 384/486，类群分布不均衡；Geminiviridae 的
"complete genome" 记录里有双组分基因组的单组分，本身近似片段。分属统计是下一步的细化方向。
