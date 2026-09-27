# 原生（Protist/Algae）与植物病毒能否按科属机械分开 —— 实测报告

日期：2026-08-30
数据源：
- `database/cross_analysis/{order,family,genus,species}_host_probability.tsv`（C7 表，构建脚本 `scripts/3_host_probability.py`，参照 `VMR_MSL41.v1.20260320.xlsx` + `ICTV_MSL41_*.tsv`）
- `06_HostPrediction/C9_ICTV_result/Plant.classified.tsv`（1,065 行，植物分支实际输入）
- `06_HostPrediction/C9_ICTV_result/classification_result.tsv`（全库 20,892 行）
- `05_Taxonomy/Votus.integrated/final_integrated_classification.tsv`

审计脚本（只读，未改任何产物）：
- `scripts/audit/plant_origin_split.py`
- `scripts/audit/plant_native_leak.py`

---

## 结论一句话

科级分不开，属级能分开一部分，而**恰好分不出最需要分的那一批**（新病毒、跨宿主科、科属种嵌合行）。
根源是宿主标签的生成方式：它不是 ICTV 权威 host range，而是"宿主记录数统计 + 最高概率"派生；记录数又与病毒的研究热度正相关，新病毒天然没有宿主证据。

---

## 一、判据底座的实际规模

| 层级 | 条目数 | Plant | Protist | Algae | 原生+藻 |
|---|---|---|---|---|---|
| Family | 337 | 35 | 13 | 3 | 16（4.7%） |
| Genus | 2,786 | 196 | 35 | 6 | 41（1.5%） |
| Species | 55,568 | — | — | — | — |
| Order | 90 | — | — | — | — |

原生/藻类科名单（16）：Allomimiviridae, Alvernaviridae, Bacilladnaviridae, Circoviridae, Leishbuviridae, Mamonoviridae, Marnaviridae, Marseilleviridae, Maviroviridae, Mimiviridae, Ouroboviridae, Phycodnaviridae, Phypoliviridae, Pseudototiviridae, Sputniviroviridae, Yaraviridae

原生/藻类属名单（41，节选）：Chlorovirus, Mimivirus, Megavirus, Pandoravirus, Phaeovirus, Prasinovirus, Tupanvirus, Fadolivirus, Bacillarnavirus, Marnavirus, Marseillevirus, Medusavirus, Cryspovirus, Leishmaniavirus …

**名单本身有两处值得质疑**：`Circoviridae → Protist`（P=0.872，仅 16 条记录）与 `Pseudototiviridae → Protist`（P=0.578，542 条记录）都不符合常识，属于记录统计噪声。

---

## 二、植物表的分层能力（Plant.classified.tsv，1,065 行）

宿主判定层级构成：

| Determination_Level | 行数 | 占比 |
|---|---|---|
| Species(via Species) | 736 | 69.1% |
| Family(via Family) | 155 | 14.6% |
| Genus(via Genus) | 145 | 13.6% |
| Order(via Order) | 29 | 2.7% |

证据强度构成：High 724 / Medium 177 / Low(Singleton/Rare) 164

可判性：

| 指标 | 植物表 | 全库 20,892 行 |
|---|---|---|
| 属级有条目（可按属判） | 811（76.2%） | 12,794（61.2%） |
| 只能科级兜底 | 254（23.8%） | 8,098（38.8%） |
| 有科级条目 | 1,004（94.3%） | — |

**属级覆盖率就是"能分开"的上限。** 有 23.8% 的植物结果行拿不到属级宿主证据，只能继承科的判定。

---

## 三、科级为什么不足以切开：跨宿主科

植物表涉及的科里，11 个科的记录同时含植物与原生/藻/真菌：

| 植物表行数 | 科 | 科级判定 | P(Plant) | P(原生+藻) | P(Fungi) | 记录数 | Shannon |
|---|---|---|---|---|---|---|---|
| 104 | Partitiviridae | Plant | 0.401 | 0.045 | **0.399** | 4,245 | 2.043 |
| 93 | Rhabdoviridae | **Mammalia** | 0.034 | 0 | 0 | 45,026 | 1.189 |
| 74 | Alphaflexiviridae | Plant | 0.998 | 0 | 0.002 | 5,642 | 0.027 |
| 51 | Tymoviridae | Plant | 0.966 | 0 | 0.006 | 1,796 | 0.278 |
| 17 | Tospoviridae | Plant | 0.968 | 0 | 0 | 3,771 | 0.210 |
| 11 | Endornaviridae | Plant | 0.755 | 0.009 | 0.155 | 1,060 | 1.199 |
| 10 | Phenuiviridae | **Human** | 0.114 | 0 | 0.001 | 19,779 | 2.010 |
| 9 | Botourmiaviridae | **Fungi** | 0.132 | 0 | 0.327 | 897 | 2.318 |
| 8 | Spinareoviridae | **Aves** | 0.125 | 0 | 0.004 | 18,807 | 2.312 |
| 5 | Aspiviridae | Plant | 0.988 | 0 | 0.002 | 899 | 0.117 |
| 1 | Amalgaviridae | Plant | 0.936 | 0 | 0.056 | 481 | 0.387 |

合计 11 科、植物表 383 行（36.0%）。

两个极端：
- **Partitiviridae** 是最坏情形。P(Plant)=0.4009 与 P(Fungi)=0.3986 几乎并列，Shannon 2.04。科级标 Plant 只是 40% 微弱多数。
- **Rhabdoviridae** 科级判 Mammalia（45,026 条记录里植物只占 3.4%），却被植物表继承。

能切与不能切的真正标尺是 **P_Max 与 Shannon**，不是科名。

---

## 四、漏洞 A：级别优先于证据强度

`utils/classify_contigs.py` 的候选排序键是 `(level_rank, Integrated_Confidence)` 降序，级别（Species=4 > Genus=3 > Family=2 > Order=1）先于置信度。
后果：一个记录数极少的物种级条目能压过属级/科级的强证据。

植物表里带原生科/属名字的行，共 10 行，全部靠 `Species(via Species)`（9 行）或 `Genus(via Genus)`（1 行）压过原生命中：

| contig_id | 科 | 属 | 种 | 判定 | 种级证据 |
|---|---|---|---|---|---|
| SRR33536652_clean_NODE_329 | Mimiviridae | Cotonvirus | Trichovirus mali | Species | Plant, 1151, High |
| SRR8191328_clean_NODE_509 | Mimiviridae | Fadolivirus | Potexvirus ecspotati | Species | Plant, 759, High |
| SRR23215176_clean_NODE_605 | Mimiviridae | Ilarvirus | Ilarvirus PNRSV | Species | Plant, 953, High |
| CRR527041_clean_NODE_1995 | Mimiviridae | Metavirus | Olive virus O | Species | Plant, **5**, Medium |
| SRR6206200_clean_NODE_10855 | Mimiviridae | Metavirus | Fabavirus vitis | Species | Plant, 102, High |
| CRR527041_clean_NODE_1593 | Mimiviridae | Potyvirus | Potyvirus rapae | Species | Plant, 2455, High |
| SRR23587951_clean_NODE_458 | Mimiviridae | Tupanvirus | Trichovirus mali | Species | Plant, 1151, High |
| contig_3578 | Pithoviridae | Alphacedratvirus | Potexvirus triasparagi | Species | Plant, **6**, Medium |
| SRR22742702_clean_NODE_2222 | Tymoviridae | Pandoravirus | Maculavirus vitis | Species | Plant, 405, High |
| SRR23215177_clean_NODE_19 | Marseilleviridae | Betanucleorhabdovirus | Coptis betanucleorhabdovirus 1 | Genus | 无种级条目 |

这些行的科、属、种三者互不相干（Mimiviridae + Tupanvirus + Trichovirus mali），是典型嵌合行。
注意方向性：不能因为"科=原生"就断言该行是原生病毒。有可能是真植物病毒被工具误挂巨病毒科，也有可能是污染序列恰好撞上植物病毒种名。判据必须能区分"哪一级可信"。

另有 29 行属级判为非植物却进了植物表（Magoulivirus→Oomycetes 5、Iflavirus→Insecta 4、Orthodiscovirus→Fungi 3、Metavirus→Fungi 2、Haifavirus→Bacteria 2、Pandoravirus/Tupanvirus/Fadolivirus→Protist 各 1 …），同样是更细层级压过属级的结果。

---

## 五、漏洞 B：科标签来自跨树投票，同一属能分裂到 7 个宿主阵营

全库 `Genus=Biavirus` 共 745 行，按 Family × Predicted_Host 分裂：

| 行数 | 科 | 宿主 | 判定层级 |
|---|---|---|---|
| 450 | Schizomimiviridae | Protist | Order(via Order) |
| 70 | Mimiviridae | Protist | Family(via Family) |
| 54 | Partitiviridae | Insecta | Species(via Species) |
| 25 | Schizomimiviridae | Protist | Species(via Species) |
| 23 | Marseilleviridae | Protist | Family(via Family) |
| 14 | Phycodnaviridae | Algae | Family(via Family) |
| **14** | **Partitiviridae** | **Plant** | **Family(via Family)** |
| 10 | Schizomimiviridae | Human | Species(via Species) |
| 9 | Adenoviridae | Human | Family(via Family) |
| 8 | Mimiviridae | Protist | Species(via Species) |
| 7 | Ascoviridae | Insecta | Family(via Family) |
| 5 | Pithoviridae | Animal_other | Order(via Order) |
| 4 | Hepaciviridae | Human | Order(via Order) |
| … | | | |

同一个属名在 10 余个不同科标签下被分到 7 个不同宿主阵营。这说明判据事实上由"科标签"这一个变量决定，而科标签来自跨树投票（C9 的 `Family_agree` 常为部分一致），本身不自洽。

**Biavirus 是这套判据的完美盲区标本**：
- `Biavirus` 在属级表**不存在**（2,786 属里没有）
- 它的真科 `Schizomimiviridae` 在科级表**不存在**
- `Biavirus raunefjordenense`（VMR MSL41 名）与 `polar freshwater` 在种级表**不存在**（表里只有 `Heliosvirus raunefjordenense|Algae|18`）

于是这 14 行的宿主只能由 `Family=Partitiviridae` 兜底，而 Partitiviridae 恰好是 P(Plant)=0.4009 / P(Fungi)=0.3986 的混合科 → 判 Plant → 经 C9 级联落进植物成果集，再被 rescue 按 Genus 分组继承。

按属切不到它，按科切会把它算成植物，按种也切不到。它没有任何一级宿主证据。

---

## 六、全库宿主构成（20,892 行）

| Predicted_Host | 行数 | 占比 | Species | Genus | Family | Order |
|---|---|---|---|---|---|---|
| Protist | 7,555 | 36.2% | 3563 | 1088 | 1801 | 1103 |
| Human | 2,611 | 12.5% | 1863 | 210 | 391 | 147 |
| Mammalia | 1,590 | 7.6% | 1359 | 118 | 100 | 13 |
| Animal_other | 1,570 | 7.5% | 523 | 191 | 82 | 774 |
| Insecta | 1,481 | 7.1% | 1224 | 132 | 111 | 14 |
| Bacteria | 1,254 | 6.0% | 963 | 158 | 118 | 15 |
| Algae | 1,235 | 5.9% | 723 | 262 | 239 | 11 |
| Plant | 1,065 | 5.1% | 736 | 145 | 155 | 29 |
| Unknown | 1,001 | 4.8% | — | — | — | — |
| Fungi | 912 | 4.4% | 289 | 275 | 344 | 4 |
| Aves | 360 | 1.7% | 290 | 70 | 0 | 0 |
| Oomycetes | 162 | 0.8% | 143 | 19 | 0 | 0 |
| Arachnida | 78 | 0.4% | 59 | 11 | 8 | 0 |
| Archaea | 18 | 0.1% | 15 | 1 | 0 | 2 |

判为原生+藻类共 **8,790 行（42.1%）**，植物 5.1%。
原生/藻类的科分布：Mimiviridae 4413、Phycodnaviridae 1281、NA 821、Schizomimiviridae 652、Marseilleviridae 416、Pithoviridae 331、Mesomimiviridae 300、Allomimiviridae 164 …
原生/藻类的属分布：NA 1447、Tupanvirus 876、Fadolivirus 774、Biavirus 606、Pandoravirus 583、Chlorovirus 582、Phaeovirus 409 …

注：原生占比 42% 这一结构本身更可能反映参照库组成与分类系统偏倚，而非真实生态。若宿主判定的科目的是论文表，这个偏斜必须单独交代。

---

## 七、四条可落地的机械判据（待定夺）

- **P1 属级优先、科级仅兜底**（现状已如此），但须叠加记录数门槛：属级条目 `Total_Records < N`（建议 15，与 C7 自身的 95% 置信拐点一致）时不参与判定。
- **P2 证据强度优先于级别**：把 `classify_contig` 排序键从 `(level_rank, IC)` 改为 `(IC, level_rank)`，或对更细层级条目设准入条件（`Confidence_Level != Low` 且 `Total_Records >= 15`）才能压过上一级。可直接修掉第四节 10 行与 29 行两类漏网。
- **P3 科级混合度门控**：科级 `Shannon_Entropy > 1.0` 或 `P_Max < 0.6` 判为 mixed family，其下无属级证据的行不进植物成果表，落 `review`。按此会拦下 Partitiviridae 104 行、Rhabdoviridae 93 行、Botourmiaviridae 9 行、Endornaviridae 11 行、Phenuiviridae 10 行、Spinareoviridae 8 行。
- **P4 原生名单剥离**：属 41 + 科 16 名单硬切，但必须与 P2 联用，否则会误伤第四节那 10 行（真植物病毒被误挂巨病毒科的情形）。

## 八、盲区（无法靠机械规则解决）

1. **新病毒无宿主记录**：Biavirus 属级/科级/种级三层全缺，任何基于"概率表查表"的方法都判不了它。只能靠外挂权威表（ICTV VMR 的 host range 字段、NCBI Virus 宿主注释）单独补一层。
2. **真菌-植物双宿主科**：Partitiviridae（40.1% / 39.9%）、Endornaviridae、Botourmiaviridae 在生物学上确实跨宿主，机械切分必然产生假阳/假阴。
3. **嵌合行**：科属种来自不同工具投票，三者互不相干时，判据无法从单行内部判断该信哪一级。
4. **记录偏斜**：Phycodnaviridae 1,078 条记录 vs Biavirus 0 条，判定能力与病毒研究热度正相关，冷门类群系统性失去宿主证据。
5. **概率表来源链已断**：C7 构建脚本 `3_host_probability.py` 的输入 `classified_clean/` 在服务器上不存在（已 find 无果），上游 C5/C6 脚本注释还写着"排除 HVT 科"（含 Partitiviridae）与 C7 的口径不一致。表是怎么生成的目前无法从服务器上复现。

---

# 附：原生/藻类与植物的科属交叉程度（层面实测）

日期：2026-08-30　脚本：`scripts/audit/plant_native_cross.py`、`scripts/audit/plant_cross_rows.py`

## 层面 1　名称层面：零交叉

科表 337 条、属表 2,786 条，每条只携带一个 `Predicted_Host`。判 Plant 的科（35）与判原生/藻的科（16）交集 **0**；判 Plant 的属（196）与判原生/藻的属（41）交集 **0**。

即不存在"同一个科名或属名同时归属植物与原生"的歧义。按科属名单剥离在语法上完全可行，剩下的问题只是标签本身对不对。

## 层面 2　记录层面：交叉很窄

337 科中同时含 Plant 与原生/藻记录者 **8 个**：

| 记录数 | 科 | 科标签 | P(Plant) | P(原生藻) | P(Fungi) | Shannon |
|---|---|---|---|---|---|---|
| 171,507 | Sedoreoviridae | Human | 0.001 | 0.000 | 0.000 | 1.109 |
| 19,779 | Phenuiviridae | Human | 0.114 | 0.000 | 0.001 | 2.010 |
| 9,429 | Peribunyaviridae | Human | 0.003 | 0.000 | 0.001 | 1.821 |
| 4,245 | Partitiviridae | Plant | 0.401 | 0.045 | 0.399 | 2.043 |
| 3,771 | Tospoviridae | Plant | 0.968 | 0.000 | 0.000 | 0.210 |
| 1,060 | Endornaviridae | Plant | 0.755 | 0.009 | 0.155 | 1.199 |
| 824 | Orthototiviridae | Fungi | 0.107 | 0.024 | 0.333 | 2.712 |
| 101 | Qinviridae | Animal_other | 0.020 | 0.089 | 0.010 | 1.841 |

2,786 属中只有 **1 个**：Orthotospovirus（Plant，P=0.969）。

要点：**植物 × 真菌的交叉远大于植物 × 原生藻**。Partitiviridae P(Fungi)=0.399 vs P(Protist)=0.042；Endornaviridae 0.155 vs 0.009。原生/藻的直接渗透面很窄。

## 层面 3　层级层面：存在交叉

属→科映射 2,929 属（源 `database/virus-db/acvirus_db/taxa.txt`，VMR MSL41）。

- 科标签=Plant 但科下有原生/藻属：**1 例**，Partitiviridae → Cryspovirus(Protist)。
- 科标签=原生/藻但科下有 Plant 属：**0 例**。
- 全表可比对属（属级+科级都有条目）1,951 个，其中属级标签与科级标签不一致 **506 个（25.9%）**。

涉原生/藻的属级/科级标签冲突（全部 10 条）：

| 科 | 科标签 | 属 | 属标签 |
|---|---|---|---|
| Betaormycoviridae | Fungi | Stormycovirus | Protist |
| Circoviridae | Protist | Circovirus | Fungi |
| Naryaviridae | Mammalia | Nimphelosvirus | Protist |
| Naryaviridae | Mammalia | Phialvirus | Protist |
| Nenyaviridae | Human | Mazarbulvirus | Protist |
| Partitiviridae | Plant | Cryspovirus | Protist |
| Phycodnaviridae | Algae | Phaeovirus | Protist |
| Phycodnaviridae | Algae | Raphidovirus | Protist |
| Pseudototiviridae | Protist | Victorivirus | Fungi |
| Sedoreoviridae | Human | Mimoreovirus | Algae |

注意后四条：`Phycodnaviridae`(Algae) 下挂 `Phaeovirus`/`Raphidovirus`(Protist)，`Sedoreoviridae`(Human) 下挂 `Mimoreovirus`(Algae)。**"原生生物"与"藻类"这条界线在表里本身就不统一**。

## 层面 4　数据行层面

全库 20,892 行：科属标签一致 8,621（41.3%）／不一致 2,687（12.9%）／**缺一侧条目 9,584（45.9%）**；不一致且涉原生/藻 943 行（4.5%）。

植物表 1,065 行：可比对 806 行，一致 667，不一致 **139（17.2%）**，其中涉原生/藻仅 **8 行**。

植物表里科级标签与属级标签不一致的主要组合：

| 行数 | 科标签 | 属标签 | 典型属 | 谁对 |
|---|---|---|---|---|
| 90 | Mammalia | Plant | Betacytorhabdovirus, Betanucleorhabdovirus | **属级对**（植物弹状病毒属） |
| 14 | Human | Plant | Rubodvirus, Coguvirus | **属级对** |
| 8 | Aves | Plant | — | **属级对** |
| 8 | Fungi | Oomycetes | Magoulivirus, Betabotoulivirus | 待定 |
| 4 | Fungi | Plant | — | 待定 |
| 3 | Protist | Plant | — | 属级对 |
| 193 | Plant | 缺属条目 | — | 无法比对 |
| 254 | （任意） | 缺属条目合计 | — | 无法比对 |

全库中科标签=Plant(Partitiviridae) 而属标签=Protist 的 13 行，属名是 Theiavirus／Tethysvirus／Tupanvirus／Alphahydrivirus／Phaeovirus 这类原生/藻巨病毒属，仍是同一批嵌合行（`Cryspovirus` 在 20,892 行里 0 次出现）。

## 交叉的形式与本问题的影响

1. **不是"同名跨阵营"，而是"同科内记录跨界 + 科属标签互不背书"。** 名称层零交叉，硬切名单可行；危险全部来自标签与证据不匹配。
2. **两个方向的错误同时存在。** 植物表 104 行（Mammalia/Human→Plant）里，科级判哺乳动物是错的、属级判植物是对的；Biavirus 那批里，科级判 Plant 是错的、属级（真值）应是原生。所以"科级优先"与"属级优先"两个简单规则都不成立。
3. **剥离原生污染的可操作面很窄。** 植物表里直接涉原生/藻的交叉只有 8 行，靠属级原生名单（41 个）+ 科级原生名单（16 个）足以覆盖。真正的大头是 Partitiviridae 那 104 行，属于**真菌-植物双宿主**问题，不是原生问题。
4. **45.9% 的行缺一侧条目**，这是科属双重佐证方法的覆盖上限。

---

# 附二　宿主判定机制与黑白名单（代码级实测，2026-09-14）

本节全部结论来自实读两端源码与在服务器上重放真实表格，脚本留存于 `scripts/audit/`，可复跑。

## 2.1 判定链路的三段

| 段 | 文件 | 位置 | 做什么 |
|---|---|---|---|
| C9 查表 | `utils/classify_contigs.py` | 服务器 `utils/`，两端 md5 `02056c31` 一致，233 行，mtime 2026-08-01 17:37 | 查 C7 四表得分 → 写 `Predicted_Host` + `Determination_Level` |
| 融合决策 | `run_host_prediction.py` | 服务器 478 行 / 23,617 B，md5 `2366d1fd` | 级联 `Rule(噬菌体纲) → ICTV(C9) → RVH → PB2 → Unknown`，另有旧并行树保留 |
| 名单 | 同上，`NON_PLANT_GENERA` / `NON_PLANT_FAMILIES_FALLBACK` | :48 / :51 | 只作用于 ICTV 这一票 |

C7 表规模（实测行数）：order 90、family 337、genus 2,786、species 55,568。C9 内部**没有任何黑白名单、没有任何置信度阈值、没有宿主先验**；候选按 `(级别, Integrated_Confidence)` 排序，Species > Genus > Family > Order，级别压倒证据强度，并含 12 种「错位匹配」（例如 `Genus(via Species*)`）。

## 2.2 名单清单（服务器现版，逐条实测）

- `NON_PLANT_GENERA` = **66 条**；`NON_PLANT_FAMILIES_FALLBACK` = **38 条**；`TRUSTED_LEVELS` = {genus, species, species\*}；`PHAGE_CLASSES` = 4 纲；`RVH_MAP` = 7 键。
- 门控原文（`is_blacklisted`），三条结构性盲区都藏在这里：

```python
def is_blacklisted(family, genus, det_level):
    if is_trusted_level(det_level):          # 盲区① 属/种级永久豁免
        return False
    gen = str(genus).strip() if genus is not None else ''
    fam = str(family).strip() if family is not None else ''
    if gen and gen not in ('NA', 'nan'):     # 盲区② 有属名则科兜底被跳过
        return gen in NON_PLANT_GENERA
    return bool(fam and fam in NON_PLANT_FAMILIES_FALLBACK)
```

- 盲区③：挂钩动作只是 `if h_ictv=='Plant' and is_blacklisted(...): h_ictv='Unknown'`，级联继续走 RVH → PB2，**下游可把否决行拉回 Plant**（下方有实证）。
- **白名单侧：不存在任何宿主白名单。** 实测到的近似物有四类：CDD 域白名单 `~/MMPV-RNA/database/cdd/cdd_virus_final_v4.txt`（`filter_virus.py:22`、`validate_rescue_cdd.py:49`，管病毒蛋白结构域，与宿主无关）；`TRUSTED_LEVELS` 豁免；`PHAGE_CLASSES` 规则；`virus_classifier.py:757` 的 `--blacklist ''`（空值）。
- `~/plant_virus_db/1.virus-host_db/C-host_classify/classified_clean/Plant.tsv` 确实存在：**201,488 行宿主-病毒对，6,247 个唯一病毒**（唯一 taxid 6,247、唯一 Accession 201,455、`Host_Category` 全为 Plant），62.2 MB。它是 C7 的补丁数据源，不参与运行时判定。审查稿写的「6,247 行记录」实为 6,247 个唯一病毒。

## 2.3 版本链：两端差异的真相

| 文件 | 行数 / 字节 | md5 | 定性 |
|---|---|---|---|
| 本地 `run_host_prediction.py` | 431（非空行 375）/ 19,980 | `0e794d41` | = 服务器 `.bak_blacklist_20260902/2/3`，**黑名单前基线（已含级联树）** |
| `.bak_cascade_20260904` | 285 / 13,295 | `07971a63` | 级联改造**前**，`normalize_c9` 为 `pd.isna` 列表版 |
| `.bak_blacklist4_20260902` | 478 / 23,459 | `3ad1b172` | 黑名单版 |
| 服务器现版 | 478 / 23,617 | `2366d1fd` | = blacklist4 + `run_ensemble` 增加 Family/Genus/Determination_Level 合并 |

审查稿「本地是没打补丁的版本、四次补丁只改服务器」成立；行数口径上，本地文件 431 行是总行数，脚本内非空行 375。

## 2.4 定量：黑名单实际动了什么

服务器 `/tmp/bl_replay.py` 用现版函数重放两棵树真实表格（`ensemble_host_summary.tsv` 加 C9 的 `Determination_Level`）：

**枸杞树（20,892 行）**：现版名单开火 **27 行**，其中历史上判 Plant 的 26 行：25 行离开 Plant（17 → Unknown、5 → Animal、2 → Fungi、1 → Bacteria），**1 行被 RVH 拉回 Plant**（盲区③的复活实证）；第 27 行历史上已是 Animal。占 Plant 行（1,092）的 2.4%。

**OneKP 树（534,438 行）**：C9 表上命中 654 行；**现表上黑名单净影响 0 行**（残留的 2 行命中行先被噬菌体纲规则判为 Bacteria）。已删 662 条 **100% 是属名单命中**（Biavirus 485 条 = 73%，Sylvanvirus 39、Rimosavirus 35、Miraophiovirus 28），判定层级 Family 638 / Order 18 / None 6，科兜底一次未用。

**补原生/藻名单的增量 = 0**：

| 树 | 现有名单已覆盖 | 漏登 | 补 10 科 + 32 属后新增拦截 |
|---|---|---|---|
| 枸杞 | 科 6/16、属 9/41 | 10 科 32 属 | **+0 行** |
| OneKP | 同上 | 同上 | **+0 行**（341 行残留全部在豁免级） |

**那 341 行（枸杞 10 行同类）为什么拦不住**：`/tmp/bl_vote.py` 逐行判读宿主票来源，341 行 = Species 293 + Genus 48，宿主票**全部来自 species/genus 表命中，且该表判 Plant**；物种名如 `Fabavirus pruni`、`Ipomovirus manihotis`、`Potexvirus papayae`、`Tritimovirus bromi`，而同行 Family 字段是 Mimiviridae / Phycodnaviridae / Marseilleviridae / Allomimiviridae。**科属打架更像科标签坏，宿主票未必错**，这批行不能据科名字就判为污染。

## 2.5 一致性红旗（本轮新发现）

用现版代码重放，两棵树的 `Final_Host` **不能复现**：枸杞变动 123 行、OneKP 变动 461 行。逐行归因后：

- OneKP 的 461 行 **全部**走 `ICTV_Preferred`，形态是 `Host_ICTV` 与 `Final_Host` 在同一行内不相容（如 `Host_ICTV='Animal'` 而 `Final_Host='Plant'`）；方向集中在 Mammalia / Protist / Algae → Animal 的**粗化**。
- 枸杞 123 行中，98 行同类（`ICTV_Preferred` 内部矛盾），25 行才是黑名单驱动。
- 时间线：OneKP 的 host stage 落在 2026-09-04 10:59（`host.log` 与 `C9_ICTV_result/` 同步），而 `.bak_cascade_20260904` 为该日 16:11；此后 `ensemble_host_summary.tsv` 又被 09-11、09-12 两次改写（`.bak_ghost_20260912`、`.bak_ghost2_20260912`，现档 09-12 15:33）。

结论：**历史宿主标签对应的代码版本无法定版，文件内部存在自相矛盾行**。任何后续「重跑对齐 / 重算宿主」之前，必须先冻结一版代码与一版 C9 输入，再全量重算，否则新旧数字不可比。

## 2.6 对「补原生/藻黑名单进行去除」的回答

1. **可以做，但不解决主要问题。** 补齐漏掉的 10 科 32 属，两棵树新增拦截均为 0 行。
2. **主因是三条结构性盲区**，不在名字表：属/种级永久豁免；有属名时科兜底被跳过；RVH/PB2 可复活被否决行。
3. **真正被现名单干掉的是 662 条（OneKP），形态与原生/藻无关**：Biavirus 占 73%，靠的是属名单 + 科级判定层级。
4. **若要拦住 341 / 10 这类行，正确的闸门是「科属自洽性检查」**：科判原生/藻而属判植物病毒属（或科为巨型病毒科而属为植物病毒属）时标记复核，属**标记**层动作，不动名字表。
5. 注意口径：OneKP 现档已删 662 条，kill 清单与 `blacklist_bak_20260902` 已不存，仅 `.bak_bl662_20260911` 等备份可反推删除范围（备份 20,716 条 vs 现 20,054 条，差 662，反向 0）。

