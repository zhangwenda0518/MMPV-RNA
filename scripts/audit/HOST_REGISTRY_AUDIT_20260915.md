# 宿主判定白/黑名单与嵌合行审计（2026-09-15）

对象：`virome_discovery_pipeline/run_host_prediction.py`（06_HostPrediction 的决策树）
关联：05 闸门（`TAX_GATE_VERSION=6.2`，科-属一致性）、`utils/classify_contigs.py`（C9 ICTV 查表）
数据：goji 树 `MMPV-paper/goji-virome/02_novel_virus/RNA-Lycium_barbarum_out`，05 表 md5 `12e8f1715922ee967db9b10cb70e523a`（20,892 contig）

---

## 一、用户三行的判定

| 询问 | 结论 | 证据 |
|---|---|---|
| Family=Mimiviridae / Genus=Potyvirus | 非植物科，**已在**名单 | Plant.tsv（201,488 条）中 Mimiviridae 零记录；原始 38 科名单含此项 |
| Family=Marseilleviridae / Genus=Betanucleorhabdovirus | 非植物科，**已在**名单 | 同上，Marseilleviridae 零记录 |
| Family=Pithoviridae / Genus=Potexvirus | 非植物科，**已在**名单 | 同上，Pithoviridae 零记录 |

三行仍出现在植物结果，原因不是漏登，而是名单在现网代码里**从未生效**（见第四节），
且属/种级判定在设计上直接绕过科名单（见第五节）。

## 二、名单规模与来源

| 名单 | 原值 | 现值 | 判据 |
|---|---|---|---|
| `NON_PLANT_FAMILIES_FALLBACK` | 38 | **227**（+199 登记，−10 撤回） | 该科在 Plant.tsv 零记录，且在 goji/OneKP 的 05 表或下游植物表中出现 |
| `NON_PLANT_GENERA` | 66 | **64**（−2 撤回） | 属级门控：C9 只给到科级时以属为判据，属在黑名单则否决 Plant |
| `PLANT_FAMILIES_WHITELIST` | — | **59**（新增） | Plant.tsv 位置口径 50 科 + ICTV 含植物但库零覆盖 9 科；只做断言，不参与决策（见第十节） |
| `PLANT_GENERA_WHITELIST` | — | **202**（新增） | Plant.tsv 位置口径全部属；同上 |

登记候选 199 科由 `scripts/audit/gen_family_registry_append.py` 生成，
逐科行数与代表属见 `scripts/audit/non_plant_family_registry.tsv`。
199 科中出现在下游植物表的 8 科：Metaviridae、Draupnirviridae、Iflaviridae、Barnaviridae、
Discoviridae、Hydriviridae、Kanorauviridae、Narnaviridae。

## 三、ICTV 全量宿主核对与撤回（10 科 + 2 属）

判定源：ICTV VMR **MSL41**（VMR_MSL41.v1.20260729）+ 官方 *Virus Properties by Family*
（含 Host 列，全表 **427 科**，`items_per_page=100` 逐页抓取）。逐科明细与 URL：
`scripts/audit/ICTV_HOST_SWEEP_233.tsv`；含植物/未定/非现行三类汇总：`ICTV_PLANT_HOST_FLAG.tsv`。
判据「Plant.tsv 零记录」在**改名/别名**与**宿主列含植物**两种情形下都会失效。

| 科 | ICTV Host 原文 | 撤回理由 |
|---|---|---|
| Ourmiaviridae | invertebrates, plants | 含植物病毒属 Ourmiavirus；原在 38 科名单内，属误登 |
| Mitoviridae | fungi, plants | 植物线粒体 mitovirus 有独立侵染证据；原在 38 科名单内 |
| Artoviridae | plants, invertebrates | 宿主列含 plants |
| Chrysoviridae | fungi, plants, invertebrates | 宿主列含 plants（成员为真菌病毒） |
| Genomoviridae | fungi, plants, invertebrates, vertebrates | 宿主列含 plants（CRESS-DNA） |
| Pestiviridae | vertebrates, invertebrates, plants | 宿主列含 plants |
| Spiciviridae | plants, invertebrates | 宿主列含 plants |
| Pseudoviridae | protists, fungi, plants, invertebrates | 宿主列含 plants，成员为 Ty1/copia 型内源反转座子 |
| Metaviridae | protists, fungi, plants, invertebrates, vertebrates | 宿主列含 plants，成员为 Ty3/gypsy 型内源反转座子 |
| Kanorauviridae | 官方表列 plants；ViralZone 记宿主未知（污水/粪便） | 两处冲突，撤回待核 |

属级复核（植物库属记录 + 物种名前缀匹配）：

| 属 | 植物库证据 | 处理 |
|---|---|---|
| Betapartitivirus | 属级记录 79 条 + 物种名 8 条（Betapartitivirus primulae / trifolii） | **撤回** |
| Geminivirus | 物种名 1 条；该名为非 ICTV 统称，Geminiviridae 在植物库 30,001 条 | **撤回** |

撤回实测影响：撤回宿主含植物的 6 科后 goji Plant 由 1028 → **1029**（翻转 1 行，
Pseudoviridae/Sirevirus/Alfalfa cytorhabdovirus 2，该行本身也是科-属-种三方不一致）；
撤回 2 属对 goji 树 0 行影响。

留存但需标注的三类（均不含植物，不影响植物召回）：
- 宿主未定/预测值：Hydriviridae（`soil (S)`）、Oomyviridae（`predicted protists`）、Smacoviridae（`predicted archaea`）。
- 官方表未收录的科名（原因已用 MSL41 `Taxa Renamed or Abolished` 表 + 提案原文核实）：

  | 科名 | MSL41 处置 |
  |---|---|
  | Reoviridae | 2021（MSL36）拆分为 Sedoreoviridae 与 Spinareoviridae（二者均为现行科，宿主含 plants，但不在本名单内） |
  | Microviridae | 2025 年废止并升格为纲 Microviricetes（提案 2025.043B），原亚科升为目 Bullavirales / Gokushovirales |
  | Autographiviridae | 2025 年由科提升为目 Autographivirales（提案 2024.045B） |
  | Adomaviridae / Cruciviridae / Pandoraviridae | 均为文献 proposed 科名，从未被 ICTV 采纳 |

  这 5 个非现行名称仍在名单内（无植物风险，保留可保护同名结果行）：
  goji 05 表中 Microviridae 2 行、Autographiviridae 5 行；onekp 05 表中 Microviridae 1 行、Autographiviridae 703 行。
- 拼写核对：任务转写文本里出现过 `Fervensiviridae`，复核代码名单与两棵树的 05 表，
  实际均为官方名 `Fervensviridae`（goji 10 行，宿主 archaea，Plants_Included=否，已在名单内），
  错拼仅存在于文档转写链路，未进入代码与数据。
- 官方分页表与 MSL 计数一致（427 = 427）；Microviridae 缺席不是抓漏，而是该科已废止（见上）。

未发现其他冲突：227 科中无任何科的成员属出现在植物库属表（位置口径 **202 属**，2026-09-15 第二轮复验，见第十节）；
64 属黑名单中除已撤回 2 条外，植物库属记录与物种名前缀均为 0。
同一断言已机器化落进代码：导入 `run_host_prediction.py` 时若白名单与黑名单相交则直接抛错。

## 四、缺陷一：C9 字段碰撞使名单恒不生效（已修）

`run_ensemble()` 把 C9 表的 `Family`/`Genus`/`Determination_Level` 并入 05 表，
而 05 表本身已有 `Family`/`Genus`。pandas merge 同名冲突列自动加后缀：

```
merge 后列: ... Family_x, Genus_x, Species, ..., Host_ICTV, Family_y, Genus_y, Determination_Level
'Family' in merged = False  →  row.get('Family') = None
```

`is_blacklisted(None, None, det_level)` 恒返回 False。实测：**命中 0 / 20,892 行**。

- 引入版本：`run_host_prediction.py.bak_blacklist4_20260902`（Sep 10 21:24）新增 `_c9_extra`；
  上一版 `bak_blacklist3_20260902` 仅取 `['contig_id','Predicted_Host']`，无碰撞。
- 修复：`_c9_extra` 只保留 `Determination_Level`。实测 C9 的 Family/Genus 与 05 表
  **20,892/20,892 行完全一致**（classify_contigs.py 以 05 表为输入直通），故 Family/Genus 冗余。
- 修复后实测：否决命中 **3,954** 行（其中 C9 原判 Plant 的 27 行）；Plant 行 1087 → **1061**。

## 五、缺陷二：属/种级判定绕过科名单（已改为科级优先）

`is_blacklisted` 原逻辑：`is_trusted_level(det_level)` 为真（genus/species 级）→ 立即放行；
属非空时也只查属名单，不查科名单。后果是同一 (科,属) 组合在属缺失时被否决、属冲突时反而不否决。
例：Barnaviridae/Sobemovirus（科=真菌巴纳病毒科，属=植物南方豆花叶病毒属）落在 Plant。

`FAMILY_FIRST_VETO = True`（现装开关，一行可回退）：科已判为非植物科即否决 Plant，与判定层无关。

三口径实测（goji，其余条件相同）：

| 口径 | 否决命中 | Plant 行 |
|---|---|---|
| 现网（名单恒不生效，等价无否决） | 0 | 1087 |
| 仅修碰撞（C9 字段） | 3954 | 1061 |
| 修碰撞 + 科级优先（现装） | 14472 | **1028** |
| 参考：改用「科级支持数 ≥ 深层支持数」才否决 | — | 1043 |

被否决的 65 行中，4 行经 RVH（pred|L1=Viridiplantae）重新落回 Plant。
**不建议**把否决扩展到 RVH/PB2 分支：其中 Aspiviridae/Miraophiovirus 一行是植物库
已收录科（Aspiviridae）下的属，RVH 的 Plant 判断与植物库一致，扩展否决会误杀。

## 六、嵌合行证据（三科 10 行）

05 表逐 rank 独立投票、无跨 rank 一致性约束，产生分类学不可能的组合：

| contig | Family / Family_agree | Genus / Genus_agree | Species / Species_agree |
|---|---|---|---|
| CRR527041_..._1593 | Mimiviridae 3/4: CAT,diamond_lca,mmseqs | Potyvirus 1/1: metabuli | Potyvirus rapae 1/3: metabuli |
| CRR527041_..._210 | Pithoviridae 2/3: CAT,diamond_lca | Potexvirus 1/1: mmseqs | NA |
| SRR23215177_..._19 | Marseilleviridae 3/4 | Betanucleorhabdovirus 1/1: metabuli | Coptis betanucleorhabdovirus 1 1/3 |
| contig_3578 | Pithoviridae 1/4: ACVirus | Alphacedratvirus 1/2: ACVirus | Potexvirus triasparagi 1/3: metabuli |

同时出现 `Realm=Riboviria + Kingdom=Bamfordvirae` 这类不可能组合。
逐 contig 追踪：9 / 302 行（goji / OneKP 植物表）全部保留在同一巨病毒科下，
植物属被替换或置空，零 contig 丢失。
C9 侧这 9 行 `Determination_Level=Species/Genus(via …)` → `is_trusted_level` 为真 →
黑名单短路；且 05 闸门只规范 Genus，Species 列仍留植物病毒种名，故闸门后仍是 Plant。

## 七、独立发现：物种级宿主参照表存在错误标注

`database/cross_analysis/species_host_probability.tsv`（55,568 行）：

| 物种 | Plant_Records | Insecta_Records | Predicted_Host | 置信 |
|---|---|---|---|---|
| Bactericera cockerelli picorna-like virus（马铃薯木虱病毒，昆虫） | 2 | 0 | Plant | Low (Singleton/Rare) |
| Potyvirus rapae | 2455 | 0 | Plant | High |

物种宿主来自样本宿主关联，并非病毒生物学宿主；昆虫病毒的物种记录被标成 Plant。
影响：Iflaviridae/Iflavirus 4 行（goji）因此落 Plant。

## 八、备份与回滚

| 备份 | md5 | 说明 |
|---|---|---|
| `run_host_prediction.py.bak_famregistry_20260915` | 2366d1fd1a61ee7b563200834e299c64 | 补丁前 |
| `run_host_prediction.py.bak_c9collision_20260915` | 53a95ee0715660e756b36e597d6a56d0 | 199 科登记后 |
| `run_host_prediction.py.bak_registryretract_20260915` | 8887b9128cc9e80a18e1c1ab622015fb | C9 碰撞修复后 |
| `run_host_prediction.py.bak_familyfirst_20260915` | 6089af9f12bd914972ef0bd4c3d95c68 | 撤回 4 科后 |
| `run_host_prediction.py.bak_genusretract_20260915` | 8b13f3574fa569d3c219ec5bc7d5fbe2 | 撤回 2 属后 |
| `run_host_prediction.py.bak_plantwl_20260915` | 96072985332d58fe68a9cfb988382f46 | 加白名单前（227 科 + 64 属） |
| `run_host_prediction.py.bak_wlgenus_20260830`（= 白名单已落、开关未开的状态） | 4c1e14a2fa7720af03c10ce9a9d842d6 | 开 W3 前（`WHITELIST_OVERRIDES_FAMILY_VETO=False`） |
| `run_host_prediction.py.bak_genusretract_20260830` | 2911b44d28c7d006914b5f82bdd5efed | 属级撤回 + `normalize_c9` 合并前（= W3 已开状态，第三轮） |
| 现装 | **a42ae81330a01958c4ce000e15d9daad** | 227 科 + **62 属**（64 − 2 撤回）+ 白名单 59 科/202 属，`FAMILY_FIRST_VETO=True`，`WHITELIST_OVERRIDES_FAMILY_VETO=True`，单一定义 `normalize_c9`（第四轮，见 10.9） |

本地留档（md5 与服务器逐一对齐）：`scripts/audit/server_snapshots/run_host_prediction.py.4c1e14a2`（开 W3 前）、
`.2911b44d`（开 W3 后、属级撤回前）、`.bak_genusretract_20260830`（内容 = `.2911b44d`）、`.a42ae813`（现装）。
回滚到第四轮前 = `cp run_host_prediction.py.bak_genusretract_20260830 run_host_prediction.py`；
回滚到 W3 开启前 = `cp run_host_prediction.py.bak_wlgenus_20260830 run_host_prediction.py`。

## 九、待办

1. 05 闸门版成品落地两棵树，再重跑 06 / 09 / 10（现树 06 产物为旧口径）。
2. 物种级宿主表可疑条目治理（昆虫/真菌病毒被标 Plant，见第七节）。
3. 内源反转座子类单列处置：Metaviridae（05 表 133 行）与 Pseudoviridae（05 表 239 行）
   在植物病毒组中的存在多为宿主基因组内源元件，不宜用「非植物科」表述，
   宜新建「宿主基因组内源元件」类别，并在论文植物病毒表中标注排除。
4. ~~对 64 属黑名单做 ICTV 属级宿主复核~~ **已闭合**（第四轮，见 10.9）：改用 ICTV VMR MSL41
   `Host source` 列逐属复核 64 属，5 属含 plant（2 属不带 (S) → **已撤回**，3 属仅 (S)/真菌 → 保留），
   7 属 VMR 无记录（保留但标注）。第一轮已完成 12 个 ICTV 含植物科的属级细查（第十节）：黑名单**不新增**，白名单已落代码。
5. 非现行科名（Microviridae / Autographiviridae / Cruciviridae / Adomaviridae / Pandoraviridae）
   在名单内保留但需在下游报告中按 MSL41 口径改名（如 Autographiviridae → 目 Autographivirales）。
6. ~~决策 `WHITELIST_OVERRIDES_FAMILY_VETO`~~ **已决策：开启**（第三轮，见 10.6 / 10.8）。
7. V2 口径（科名单完全让位于属/种级判定）**已否**，但其额外 102 行（goji 19 + onekp 83）
   作为「科/属标签与种名互斥」的待人工复核清单留档（`modes_v2minusv1_*.tsv`，见 10.8）。
8. **待大王定**：争议条目 4 条的处置（第四轮未动，证据与影响量已量出，见 10.9 表）：
   Sylvanvirus（onekp 688 行 / 35 行翻转，名字在 ICTV 属表、管线参考库、VMR 三处皆 0 命中）、
   Crucivirus（60 / 2）、Rimosavirus（onekp 35 / 33、goji 3 / 3；VMR 无确认植物宿主、仅 (S)）、
   科级 Discoviridae（goji 剥离后 +3 行 Plant，VMR 有 1 条 plants 记录）/ Ouroboviridae（+0）。
9. **待确认**：VMR `Host source` 列 `(S)` 后缀的书面定义。按取值形态推断为「样本/环境来源，非确认宿主」
   （(S) 只出现在 soil / air / sewage / freshwater / marine / phytobiome / unknown / plants / protists / invertebrates 上），
   但未取到 ICTV 书面说明；10.9 的撤回判据依赖该语义，若将来定义不同需重跑 `vmr_genus_host_audit.py`。

## 十、12 个 ICTV 含植物科的属级细查与白名单落地（2026-09-15 第二轮）

### 10.1 解析口径
Plant.tsv 按**位置口径**：`Virus_lineage` 恰 9 段（201,488/201,488 行），
index5=Family / index6=Genus，与 `3_host_probability.py` 的 `levels` 完全一致 → **50 科 / 202 属**。
旧「锚定口径」（取科名后一个词）漏掉 Alphasatellitidae / Avsunviroidae / Pospiviroidae /
Tolecusatellitidae，只得到 46 科 / 175 属，已弃用。

### 10.2 属级证据分层（12 科共 116 个 (科,属) 组，`PLANT_FAMILY_GENUS_TRIAGE.tsv`）

| 判定 | 组数 | 定义 |
|---|---|---|
| E1/E2 植物属-白名单 | 21 | Plant.tsv 有该属记录，或属级概率表 Plant_Records>0 |
| E3 非植物属-候选 | 58 | 属级概率表有条目、Plant_Records=0、Predicted_Host≠Plant |
| E4 证据不足-不动 | 37 | Plant.tsv 与属级概率表都没有该属 |

科内构成（新出现的 3 科全部为植物属）：

| 科 | Plant.tsv 行/属 | 属级结果 |
|---|---|---|
| Tomosaviridae | 53 / 1 属 | Virtovirus 53/53 全植物，无杂属 |
| Virgaviridae | 6,077 / 7 属 | Tobamovirus、Tobravirus、Pomovirus、Furovirus、Pecluvirus、Hordeivirus、Goravirus 全 Plant；另现 13 个非植物属，均为嵌合行（Pandoravirus、Negevirus、Betabaculovirus、Lentivirus 等） |
| Solemoviridae | 7,715 / 4 属 | Polerovirus、Sobemovirus、Enamovirus、Polemovirus 全 Plant，另有 Potyvirus（Plant.tsv 归 Potyviridae，30,348 条）落在本科行上；非植物属仅 Pandoravirus、Hubsclerovirus |

9 个既有科（Artoviridae / Chrysoviridae / Genomoviridae / Kanorauviridae / Mitoviridae /
Ourmiaviridae / Pestiviridae / Pseudoviridae / Spiciviridae）Plant.tsv 覆盖均为 0（实测 `plant_tsv_family_coverage.py`）：
库能认得的是少数植物属，它们作为**属列**落在这些科（及其它非植物科）的行上：

| 属 | Plant.tsv 记录数 | Plant.tsv 里的科 | 出现在（科级标签） |
|---|---|---|---|
| Begomovirus | 25,168 | Geminiviridae | Genomoviridae |
| Potyvirus | 30,348 | Potyviridae | Solemoviridae / Orthoherpesviridae 等 |
| Orthotospovirus | 3,519 | Tospoviridae | Mitoviridae |
| Caulimovirus | 1,080 | Caulimoviridae | Pseudoviridae |
| Cilevirus | 775 | Kitaviridae | Virgaviridae |
| Maculavirus | 552 | Tymoviridae | Mitoviridae |
| Tungrovirus | 159 | Caulimoviridae | Pseudoviridae |
| Totivirus | 41 | Orthototiviridae | Mitoviridae |
| Ourmiavirus | 40 | Botourmiaviridae | Ourmiaviridae |

同一行的科级标签与属级来源分属不同植物科，即两个 rank 被投给了不同科（嵌合/错配），
属级证据才是真实来源，因此这些行必须保持 Plant，不能被属黑名单拦掉。

### 10.3 黑名单**不新增**（实测依据）

保守候选 40 个（宿主 ∈ Insecta/Human/Mammalia/Aves/Animal_other/Algae/Protist/Bacteria，
且 Plant.tsv 零记录），其中 32 个不在现装名单。把这 32 属一次性加入后对拍：

| 树 | 行数 | 现状 Plant | 加入后 Plant | 翻转 |
|---|---|---|---|---|
| goji | 20,892 | 1,029 | 1,029 | **0** |
| onekp | 535,103 | 19,560 | 19,560 | **0** |

原因：属级否决只在判定层为科/目级（`is_trusted_level` 为假）时才会被查到，
而候选属所在行的判定层是属/种级（C9 依植物属直接给 Plant）或 C9 本就不判 Plant。
加了不生效，还会把「库缺口型」属（Mitovirus 等）拖回误否决，故**不新增**。

### 10.4 为何排除 Fungi 宿主属

ICTV 把这 9 科的宿主列写成含 plants，而 Plant.tsv 对这些科零覆盖。
对 Alphachrysovirus / Betachrysovirus / Chrysovirus / Circovirus / Duamitovirus /
Gemycircularvirus / Mitovirus / Triamitovirus / Unuamitovirus 而言，
「属级宿主=Fungi」更可能是**库缺植物记录**造成的假非植物，一律不拉黑。

### 10.5 白名单已落代码（断言型，零行为变更）

- `PLANT_FAMILIES_WHITELIST`：59 = Plant.tsv 50 科 + 零覆盖的 9 个 ICTV 含植物科（即第三节撤回的全部科名）。
- `PLANT_GENERA_WHITELIST`：202（Plant.tsv 位置口径全部属）。
- `ICTV_PLANT_INCLUDING_FAMILIES`：12 科 → ICTV Host 列原文。
- 导入时 fail-loud 断言「白名单 ∩ 黑名单 = ∅」，报错文案直接指向本文档。
- 实测：科交集空 / 属交集空 / 12 科均未被拉黑 / 12 科均已在白名单；
  两棵树 Plant 计数不变（goji 1029、onekp 19560）。
- 白名单**不覆盖科级否决**：科级否决是刻意的嵌合拦截（Mimiviridae+Potyvirus 等），
  越权会反向制造假阳性。如需越权，走 10.6 的显式开关。

### 10.6 白名单越权（W3）：**已决策开启**（2026-09-15 第三轮）

背景：现状 `FAMILY_FIRST_VETO=True` 让科级否决先于判定层检查，因此 C9 依植物属给出 Plant 的嵌合行仍被否决：

| 树 | 科级否决行 | 其中属命中植物白名单 | 属级否决行 | C9=Plant 总行 |
|---|---|---|---|---|
| goji | 44 | 15 | 20 | 1,065 |
| onekp | 579 | 463 | 636 | 19,210 |

W3 口径（科命中黑名单时，若判定层为属/种级且属在植物白名单 → 放行）对拍：

| 树 | 现状 Plant | W3 Plant | 净增 |
|---|---|---|---|
| goji | 1,029 | 1,042 | +13 |
| onekp | 19,560 | 19,967 | **+407** |

被救回行示例：Mimiviridae+Potyvirus（det=Species）、Orthoherpesviridae+Macluravirus（43 行）、
Mimiviridae+Potexvirus（13 行）、Barnaviridae+Sobemovirus、Pithoviridae+Potexvirus、Nimaviridae+Allexivirus。
代价：若某行科名正确而属名是假的（真为非植物病毒），会引入 Plant 假阳性。
明细 `w3_rescued_goji.tsv` / `w3_rescued_onekp.tsv`。

**决策**（依据 10.8 的票型/内在一致性证据）：开启 W3。
`WHITELIST_OVERRIDES_FAMILY_VETO = True`（415 行），现装 md5 `2911b44d28c7d006914b5f82bdd5efed`，
备份 `run_host_prediction.py.bak_wlgenus_20260830`（md5 `4c1e14a2fa7720af03c10ce9a9d842d6`）。
开启后实测（`verify_w3_on.py`）：goji Plant 1,029 → **1,042**、onekp 19,560 → **19,967**；
与开启前的 V1 预测逐行差异 **0**，与 V0 的差集恰好 13 / 407 行且每行属均在 `PLANT_GENERA_WHITELIST` 内（PASS）。

同时**否决 V2**（把 `FAMILY_FIRST_VETO` 置否、让科名单完全让位于属/种级判定）：
多出的 goji 19 / onekp 83 行里，属级证据不可信（onekp 83 行中 55 行属为空、其余为非植物属
Chlorovirus / Arenavirus / Simplexvirus / Alphabaculovirus 等），种名又出现节肢动物相关病毒
（Tetranychus urticae-associated picorna-like virus、Bactericera cockerelli picorna-like virus），
属于「科名与种名互相矛盾的弱证据行」，不宜自动放行。该 102 行已落盘为待人工复核清单。

### 10.7 本轮改动的验证（逐行等价 + 独立复算）

验的是「这次改动到底动了什么、能不能动到判定结果」：

1. **文本级**（`diff_check_plantwl.py`）：备份（743 行）vs 现装（857 行）差异 = 新增 3 块共 114 行、
   删除 0 行、替换 0 块。三个插入点：白名单常量前、`FAMILY_FIRST_VETO = True` 后、
   `is_blacklisted` 科级否决分支内。黑名单定义区未被动过。
2. **逐行等价**（`regress_plantwl_equivalence.py`）：同一份输入分别交给备份版与现装版判完所有行。
   goji 20,892 行逐行差异 **0**（Plant 两侧均 1029）；onekp 535,103 行逐行差异 **0**（两侧均 19560）；
   两棵树出现过的 (Family, Genus, Determination_Level) 三元组共 1667 / 5075 个，
   两版 `is_blacklisted` 判定不一致 **0** 个。
3. **白名单独立复算**（`verify_whitelist_vs_planttsv.py`）：用 pandas 的 `split` 与纯文本逐行切分两种实现
   重算 Plant.tsv → 均为 50 科 / 202 属且互相一致；Plant.tsv 里的科/属全部在白名单内，
   属白名单无多余条目，科白名单多出的 9 个正是零覆盖的 ICTV 含植物科。
4. **附带面**（`check_names_and_importers.sh`）：新引入的 4 个名字在文件里只出现在预期位置；
   全管线无其他文件 `import run_host_prediction`，故导入期断言不会连坐其他 stage；
   全管线只有本文件解析 Plant.tsv，不存在需同步的第二份名单。
5. ICTV 12 科宿主原文已用本地 `ICTV_HOST_SWEEP_233.tsv` 逐行复核（Plants_Included 均为是，与代码常量一致）。

结论：本轮改动对两棵树 **零行为变更**（逐行可验），只增加了「防自相矛盾的断言」与「一个默认关闭的开关」。
其中唯一带判断性的两点需显式知情：
- 导入期断言是 **fail-loud**（白黑名单相交时直接抛错，不降级为警告）。好处是不会静默出错，
  代价是若将来要故意把某个植物库在册的科拉黑，会先卡在导入上，需改这一行或改文档。
- `FAMILY_FIRST_VETO=True`（上一轮引入）才是真正影响面最大的一条：onekp 579 行科级否决中
  463 行的属是植物库认证的植物病毒属。要不要改成 10.6 的 W3，是召回与科名优先之间的取舍。
  **第三轮已定：开 W3（`WHITELIST_OVERRIDES_FAMILY_VETO=True`），`FAMILY_FIRST_VETO` 保持 True**，证据见 10.8。

### 10.8 第三轮证据链：为何开 W3、为何否 V2

三口径对拍（`compare_veto_modes2.py`，逐行重新判完两棵树）：

| 口径 | 定义 | goji Plant | onekp Plant |
|---|---|---|---|
| V0 | 现状（科级否决优先，W3 关） | 1,029 | 19,560 |
| V1 | 科级命中黑名单时，属/种级判定 + 属在植物白名单 → 放行（= W3 开） | 1,042 | 19,967 |
| V2 | 再加：科名单完全让位于属/种级判定 | 1,061 | 20,050 |

**V1 放行行（goji 13 / onekp 407）的构成**：属命中植物白名单 **407/407 与 13/13（100%）**。
- 科标签（那一个异类）：Mimiviridae 174、Orthoherpesviridae 44、Phycodnaviridae 34、Pithoviridae 20、Nimaviridae 19、Marseilleviridae 18。
- 属标签：Macluravirus 45、Potyvirus 38、Maculavirus 26、Potexvirus 22、Nepovirus 19、Ampelovirus 19。
- 种名：Rumex macluravirus 1（44）、Maculavirus vitis（26）、Allexivirus deltallii（11）、Potexvirus papayae（10）、Blueberry virus L（9）、Sequivirus stellatum、Arepavirus arecae。
- 判定层：Species(via Species) 342 + Genus(via Genus) 65；Confidence High 294。
- 票型（`rank_agreement_delta.py`）：科级票 **407/407 均为非一致票**（形如 3/4、2/3），属级票 192 行一致。

**V2 额外行（goji 19 / onekp 83）的构成**：属命中白名单 **0**。
- onekp 83 行中 **55 行属为空**，其余为非植物属：Chlorovirus 5、Lymphocystivirus 2、Klosneuvirus 2、Alphabaculovirus 2、Arenavirus 2、Prasinovirus 1、Coccolithovirus 1、Iltovirus 1、Simplexvirus 1、Tupanvirus 1。
- 种名里仍有植物病毒种（Iris potyvirus A 16、Ullucus tymovirus 1 7、Vaccinium-associated virus C 4、Tomato necrotic ring virus 4），
  但 goji 19 行出现节肢动物相关种名（Bactericera cockerelli picorna-like virus 4、
  Grapevine wood holobiome associated mycobunyavirales-like virus 2 3），属标签也有 Iflavirus 4。
- 置信度：Low (Singleton/Rare) 59/83。

**判据（为什么 V1 行可信、V2 行不可信）**：这两类行都是「同一行内科名与属/种名在 ICTV 树系上互斥」（Mimiviridae + Potyvirus 不可能同属一条谱系），
所以至少有一个是错的。区别在：
- V1 行里**属名与种名前缀自洽**（Sequivirus + Sequivirus stellatum、Arepavirus + Arepavirus arecae、Maculavirus + Maculavirus vitis），
  科名是那个孤立异类；且属级证据可由 Plant.tsv 独立认证。
- V2 行里**属名与种名也互斥**（Chlorovirus + Paris capillovirus 1、Varicellovirus + Ullucus tymovirus 1，或属干脆缺失，Genus_agree 形如 0/0），
  即属、科两个标签都不可信，只剩一个低置信度的种级票支撑，故不自动放行。

**残留不确定性（已知情）**：V2 额外行的种名 **82/82** 都能在 Plant.tsv 里查到（C9 的 Plant 不是凭空来的）；
即「种名在植物库里查得到」是否已足以采信，是一个独立判断。若将来定下这条，把 `FAMILY_FIRST_VETO` 置否即可，
脚本与清单已就位（`modes_v2minusv1_goji_19.tsv` / `modes_v2minusv1_onekp_83.tsv`）。
另：`/tmp/labels_<tree>_V0.tsv`、`_V1.tsv`、`_V2.tsv` 为三口径的逐行标签快照，可用于任何后续对拍。

### 10.9 第四轮：属级名单的 ICTV VMR 权威复核与 E1/E2 修复（2026-08-30）

背景：64 属黑名单的原始判据是「植物库（Plant.tsv）零记录」，这是**库缺口敏感**的判据（库没收录 ≠ 宿主不是植物）。
本轮改用 ICTV 自己的宿主记录做交叉复核：`~/plant_virus_db/1.virus-host_db/B-ictv/VMR_MSL41.tsv`
（19,272 行 × 34 列，col26 = `Host source`）。

**`Host source` 取值实测**（`probe_vmr_hostsource_semantics.py`）：24 个不同取值，含 plant 的只有三种写法
`plants`(2,926) / `invertebrates, plants`(85) / `plants (S)`(68)，共 3,079 条。
带 `(S)` 的取值全部是 soil/unknown/invertebrates/freshwater/plants/protists/marine/sewage/phytobiome/air，
即**环境或样本来源**型取值；不带 (S) 的取值是 bacteria/vertebrates/plants/invertebrates/fungi/archaea/protists，
即**宿主界**。故本轮判据：`Host source` 写 plants 且**不带 (S)** = ICTV 确认的植物宿主；只带 `(S)` 不计。
（该语义为按取值形态推断，书面定义待确认，见第九节待办 9。）

**64 属逐属复核结果**（`vmr_genus_host_audit.py`）：

| 类型 | 属数 | 明细 |
|---|---|---|
| VMR 宿主列含 plants 且**不带 (S)** | **2** | Miraophiovirus `plants`(13 条，Aspiviridae)、Oleurovirus `plants`(1 条，Geminiviridae) |
| 只有 (S) 或真菌 | 3 | Alphambiguivirus 真菌 21 + plants(S) 1、Betambiguivirus 真菌 5 + plants(S) 2、Rimosavirus plants(S) 3 + 无脊椎(S) 4 + 真菌 1 + 脊椎 1 |
| VMR 无该属记录 | 7 | Arenavirus、Clandestinovirus、Crucivirus、Hokovirus、Indivirus、Klosneuvirus、Sylvanvirus |

名字溯源（`grep -w`，`ICTV_MSL41_Genus.tsv` 与管线参考库 `final.cluster.ref_info.tsv`）：
Miraophiovirus / Oleurovirus 两处**均在册**（ICTV 属表 1 命中；参考库 17 / 1 命中）→ 确为 ICTV 现行属；
同一属表里 Alphambiguivirus / Betambiguivirus / Rimosavirus 也各 1 命中（Rimosavirus 参考库 3 命中）→ 为现行属、但宿主记录仅 (S)。
而 Arenavirus / Sylvanvirus / Crucivirus / Klosneuvirus / Hokovirus / Indivirus / Clandestinovirus 在
**ICTV 属表与管线参考库均 0 命中**（Arenavirus 为已废止名，现行对应属为 Mammarenavirus 等），
其中后 6 个只在 05 整合表的 `Genus` 列出现（provenance 待查，见待决条目）。

**处置：撤回 2 属**（`patch_genusretract_e1e2.py`，备份 `run_host_prediction.py.bak_genusretract_20260830`，
现装 md5 `2911b44d…` → `a42ae813…`）。撤回后内存对拍实测：

| 树 | 撤回前 Plant | 撤回后 Plant | 增量 | 增量行构成 |
|---|---|---|---|---|
| goji | 1,042 | 1,042 | **+0** | 该树 Miraophiovirus 1 行已被 RVH 复活为 Plant，撤回不再改变 |
| onekp | 19,967 | 19,990 | **+23** | Miraophiovirus 20 行（Aspiviridae + C9 判 Plant，原被属黑名单否决）、Oleurovirus 3 行（Geminiviridae） |

上线验收（`verify_genusretract_on.py`，PASS）：两属已不在 `NON_PLANT_GENERA`、已入 `NON_PLANT_GENERA_UNDER_REVIEW`（4 项）；
白名单∩黑名单 = ∅；`def normalize_c9` 仅 1 处且 `'NA'/'nan'/''/None` → `'Unknown'`；
与 V1 基线逐行对拍 goji **0** 差异、onekp **23** 行且差异行属集合恰为 {Miraophiovirus, Oleurovirus}。

**E2 修复：`normalize_c9` 双定义合并**。文件里原有两处定义（旧 497 行 = 带 strip/NA 归一；655 行 = 生效版，无 strip）。
输入域实测（`probe_c9_domain.py`）：C9 的 `Predicted_Host` 是 **14 类封闭取值**（Protist/Human/Mammalia/`Animal_other`/
Insecta/Bacteria/Algae/Plant/Unknown/Fungi/Aves/Oomycetes/Arachnida/Archaea），两棵树**对 normalize 敏感的原值 0 行**；
两版整树对拍判定差异 **0 / 0**。即两版在真实数据上等价，差别只在「未来出现 NA/空白时是否把无效值当有效宿主类别」。
处置：删掉死副本，保留**带 strip + NA 归一的版本**（更强、语义与 `read_ictv_hits` 一致），并用上述 docstring 说明合并理由。

**E3 结论：原「`'species*'` 不可达」的说法被证伪，不改代码**。C9 的 `Determination_Level` 真实取值含
`Species(via Species)`、`Family(via Family)`、`Genus(via Genus)`、`Order(via Order)`、`Species*(via Genus)` 与 NaN；
`Species*(via Genus)` 出现在 **goji 69 行 / onekp 2,181 行**，`is_trusted_level` 切分 `(` 后得 `species*` → 命中 `TRUSTED_LEVELS`。
删掉 `'species*'` 会把这 2,250 行的判定层从「可信」降级，属**引入新错误**，故保留。

**S1 盲区计数**（否决非终局，RVH/PB2 可复活）：goji 科级否决命中 14,162 行，最终仍判 Plant **2** 行；
onekp 371,378 行 → 最终 Plant **148** 行（首因 Partitiviridae + Biavirus 的 RVH 复活 102 行，其次 Aspiviridae + Miraophiovirus 6 行 RVH + 2 行 PB2）。
即「否决可被复活」这条通道的规模在 0.04%–0.9% 量级，不构成结构性风险，但需在论文方法里如实描述。

**产物侧待对账（非本轮改动引入，但受本轮改动影响）**：onekp 树的 `06_HostPrediction/ensemble_host_summary.tsv`
被一次性后处理删过一轮（备份 `ensemble_host_summary.tsv.bak_bl662_20260911`，另有 `.bak_ghost_20260912` /
`.bak_ghost2_20260912`）。行数（不含表头）：9/4 备份 535,103 → 现 534,438（少 665）。

本轮用 csv 引号口径逐列对账（`reconcile_mira_counts.py` / `dump_mira_32rows_bak.py`，两份脚本均为只读）：

| 口径 | Miraophiovirus | Oleurovirus | 说明 |
|---|---|---|---|
| onekp 05 表 `Genus` 列（=备份表全量） | 29 | 4 | 28 行 Aspiviridae + 1 行 Mimiviridae（嵌合）；4 行 Geminiviridae |
| 9/4 备份表里 `Final_Host=Plant` | 28 | 4 | 按 Decision_Method 拆：ICTV_Preferred 20+3、ICTV_RVH_Agree 6+1、ICTV_PB2_Agree_Tiebreaker 2 |
| **现装模块撤回后 label 翻转** | **20** | **3** | = 上表里 ICTV_Preferred 那批（它们是旧名单下**真被否决**的行） |
| 现表/ghost/ghost2 保留 | 1 | 0 | 保留的 1 行是 Mimiviridae + C9=Protist 那行 |

差额 9 行（ICTV_RVH_Agree 6+1 与 ICTV_PB2_Agree_Tiebreaker 2）的解释：这些行的 Plant 判定走的是
RVH/PB2 一致或平局裁决路径，**模块从来没用属黑名单否决过它们**，但 9/11 那次后处理是按属名无条件删行，
于是把它们一并删了。即：那次清理比模块自己的规则更狠，多删了 9 行真 Plant。
这属独立缺陷（后处理脚本已不在库内，只剩备份可供对账），不手工改表，归入第九节待办 1，
重跑 06/09/10 后以管线产物为准；备份仅用于对账，不能当回滚点（回滚会重新引入旧名单口径）。

方法学提醒：05 表与 summary 的字段**全带双引号**，任何按列取值的核对脚本必须用 csv 引号解析，
直接用 `line.split('\t')` 会拿到 `"Miraophiovirus"` 而误报 0 命中（本轮实测踩过一次）。
`grep -c` 类全字符串统计不受影响。

**本轮未动的待决条目**（证据与影响量已量出，等大王定）：

| 条目 | 类型 | 树内行数 / 剥离后翻转 | 现有证据 |
|---|---|---|---|
| Sylvanvirus | 属黑名单 | onekp 688 / **35 行 → Plant**（其中 30 行判定层为科级、科名为 Caulimoviridae/Potyviridae/Secoviridae 等植物专有科）；goji 0 | 名字在 ICTV 属表、管线参考库、VMR 三处 **0 命中**，provenance 不明 |
| Crucivirus | 属黑名单 | 60 / 2 行 → Plant（Geminiviridae 科级） | 同上三处 0 命中；该名是 CRESS-DNA 病毒的非 ICTV 惯用名 |
| Rimosavirus | 属黑名单 | onekp 35 / 33 行 → Plant；goji 3 / 3 | VMR 只有 plants(S) 3 条，无确认植物宿主 |
| Discoviridae | 科兜底黑名单 | goji 17 行 / 剥离后 **+3 行 Plant**；onekp 2 行 / +0 | VMR 有 1 条非 (S) plants 记录：rice dwarf-associated bunya-like virus（Orthodiscovirus） |
| Ouroboviridae | 科兜底黑名单 | goji 0 行；onekp 4 行 / +0 | VMR 有 2 条非 (S) plants 记录：Trichosanthes kirilowii CRESS virus pt11（Demetevirus）、Forsythia suspensa CRESS virus pt110-nan（Persevirus） |
| Klosneuvirus / Hokovirus / Indivirus | 属黑名单 | 各 3 行 → Plant（均为 Partitiviridae 嵌合行） | 三处 0 命中；属名来自巨型病毒惯用名，保留否决更安全 |
| Clandestinovirus / Arenavirus | 属黑名单 | 0 / 0 与 22 / 0 | 无翻转，保留 |

判定建议（待大王裁）：Sylvanvirus 的 35 行里 30 行落在植物专有科，若这个属名本身是上游参考库的伪名，
则属黑名单条目在「科级判定」分支上**反向惩罚了召回**（属名缺失时反而会走科兜底放行）；
建议先查 05 表 `Genus` 列的 provenance（哪个中间文件、哪条参考序列），再决定撤回或保留。
Crucivirus / Rimosavirus / Klosneuvirus 三组影响量 ≤ 3 行，改动收益小，建议保持现状并记录。
Discoviridae 的 +3 行与 Ouroboviridae 的 +0 属「科兜底条目 与 该科例示病毒的宿主记录」冲突：
两科的 **ICTV 科级 Host 列不含 plants**（含则第三节就已按同一规则撤回），但 VMR 的例示病毒 Host source
出现**非 (S) 的 plants 记录**，即 ICTV 自己的两处数据自相矛盾。按哪一处为准，需大王定口径；
若定「以 VMR 例示病毒宿主为准」，撤回这两科的相应改动 = `NON_PLANT_FAMILIES_FALLBACK` 删 2 项，
需同时改 `PLANT_FAMILIES_WHITELIST`（两科目前不在 59 科白名单内，撤回后应同步入白名单，否则断言会失败）。

## 十一、产物与脚本索引

ICTV 宿主核对产物（可复现）：
- `scripts/audit/fetch_ictv_virus_properties.py`（抓取脚本）
- `scripts/audit/ICTV_HOST_SWEEP_233.tsv`（237 科逐科宿主，md5 `18B4746CD6D9F5502B124791B7B1EC76`）
- `scripts/audit/ICTV_PLANT_HOST_FLAG.tsv`（21 条：含植物 12 / 宿主未定 3 / 非现行科名 6）

`three_family_presence.py`、`three_family_trace.py`、`tool_vote_trace.py`、`list_giantfam_rows.py`、
`dump_tax_rows.py`、`dump_c9_rows.py`、`dump_agree_cols.py`、`diag_merge_columns.py`、
`diag_host_blacklist_effect.py`、`verify_c9fix_and_gate.py`、`verify_c9fix_gate_focus.py`、
`quant_parked_families.py`、`dump_eight_family_rows.py`、`quant_family_first_rule.py`、
`quant_bl_rule_variants.py`、`final_effect_report.py`、`probe_residual_rvh_plant.py`、
`audit_registry_lists.py`、`sweep_registry_vs_plantdb.py`、`sweep_genera_vs_plantdb.py`、
`sweep_genera_alias.py`、`sweep_genera_family_health.py`、`dump_genera_list.py`、
`verify_registry_state.py`、`quant_genus_retract.py`

12 科属级细查与白名单（2026-09-15 第二轮，可复现）：
- `scripts/audit/triage_v2.py` → `PLANT_FAMILY_GENUS_TRIAGE.tsv`（116 行属级分层，见第十节）
- `scripts/audit/gen_plant_whitelist_block.py` → 生成白名单代码块（`plant_whitelist_block_generated.py`）
- `scripts/audit/patch_insert_whitelist.py`、`scripts/audit/patch_add_wl_switch.py`（幂等插入补丁）
- `scripts/audit/verify_registry.py`（导入即断言：白黑不相交 / 12 科未被拉黑）
- `scripts/audit/effect_genus_additions.py`（加入 32 候选属的影响实测：两棵树翻转 0 行）
- `scripts/audit/diag_veto_vs_plantgenus.py`（科级否决命中与植物属代价面）
- `scripts/audit/diag_whitelist_override.py` → `w3_rescued_goji.tsv` / `w3_rescued_onekp.tsv`
- `scripts/audit/diff_check_plantwl.py`（文本级：差异必须全部为新增行）
- `scripts/audit/regress_plantwl_equivalence.py`（逐行等价回归：备份版 vs 现装版，两棵树）
- `scripts/audit/verify_whitelist_vs_planttsv.py`（白名单 vs Plant.tsv 两种独立实现复算）
- `scripts/audit/check_names_and_importers.sh`（命名冲突 / 导入方 / 同类名单重复维护 检查）
- `scripts/audit/plant_tsv_family_coverage.py`、`scripts/audit/plant_tsv_genus_family_dist.py`
- `scripts/audit/plant_tsv_family_coverage.py`（12 科在 Plant.tsv 的科级行数/属数，复现 10.2 表）
- `scripts/audit/plant_tsv_genus_family_dist.py`（植物属在 Plant.tsv 的科归属，复现 10.2 下表的嵌合溯源）
- 服务器快照：`scripts/audit/server_snapshots/run_host_prediction.py.4c1e14a2`

第三轮（W3 决策与开启，2026-09-15）产物与脚本：
- `scripts/audit/compare_veto_modes2.py`（V0/V1/V2 三口径逐行对拍；同时落盘 `/tmp/labels_<tree>_<mode>.tsv`）
- `scripts/audit/rank_agreement_delta.py`、`scripts/audit/rank_agreement_delta2.py`（放行行的科/属/种票型与 Plant.tsv 查证）
- `scripts/audit/trace_v2_c9_rows.py`（差集行回接 C9 原始判定；`trace_v2_rows_upstream.py` 为早期版本，会超时，已弃用）
- `scripts/audit/enable_w3.sh`（备份 + 把开关置 True + py_compile + md5）
- `scripts/audit/verify_w3_on.py`（开完后的验收：与 V1 预测逐行 0 差异；与 V0 差集全在植物属白名单内）
- 明细清单：`modes_v1minusv0_goji_13.tsv`、`modes_v1minusv0_onekp_407.tsv`（已放行）
  与 `modes_v2minusv1_goji_19.tsv`、`modes_v2minusv1_onekp_83.tsv`（未放行，待人工复核）
- 服务器快照：`scripts/audit/server_snapshots/run_host_prediction.py.2911b44d`（开 W3 后）

第四轮（VMR 属级复核 + E1/E2 修复，2026-08-30）产物与脚本：
- `scripts/audit/probe_vmr_hostsource_semantics.py`（VMR `Host source` 取值分布与 (S) 形态，复现 10.9 的取值实测）
- `scripts/audit/vmr_genus_host_audit.py`（64 属 + 227 科逐条对 VMR 宿主列复核）
- `scripts/audit/vmr_plant_entry_detail.py`（VMR 含植物记录逐条展开到种名/属名：Miraophiovirus / Oleurovirus /
  Ouroboviridae / Discoviridae / Rimosavirus）
- `scripts/audit/vmr_contested_final.py`、`scripts/audit/trace_novmr_genera.py`（无 VMR 记录属的埋名溯源与逐属撤回影响）
- `scripts/audit/measure_genus_retract_e1.py`（内存撤回 Miraophiovirus / Oleurovirus 的两棵树增量实证）
- `scripts/audit/probe_c9_domain.py`（C9 `Predicted_Host` / `Determination_Level` 取值域 + 两版 `normalize_c9` 对拍）
- `scripts/audit/patch_genusretract_e1e2.py`（撤回 2 属 + 合并 `normalize_c9`，带备份/断言/py_compile/md5）
- `scripts/audit/verify_genusretract_on.py`（上线验收：名单状态 / 单一定义语义 / 两树计数 / 与 V1 逐行差集）
- `scripts/audit/measure_family_retract_contended.py`（争议科 Ouroboviridae / Discoviridae / Ambiguiviridae /
  Steitzviridae / Fiersviridae 的树内行数与剥离后翻转量）
- `scripts/audit/reconcile_mira_counts.py`（05 表 / C9 / summary 三口径对账，csv 引号解析）
- `scripts/audit/dump_mira_32rows_bak.py`（33 行逐行 Decision_Method 拆解，证明 32 vs 23 的差额来源）
- 服务器快照：`scripts/audit/server_snapshots/run_host_prediction.py.a42ae813`（现装）、
  `…py.2911b44d`（W3 开启后）、`…py.bak_genusretract_20260830`（= `2911b44d`）、`…py.4c1e14a2`（W3 开启前）
- 第三口径明细清单：`modes_v1minusv0_goji_13.tsv`、`modes_v1minusv0_onekp_407.tsv`、
  `modes_v2minusv1_goji_19.tsv`、`modes_v2minusv1_onekp_83.tsv`
