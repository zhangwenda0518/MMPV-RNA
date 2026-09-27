# phylo_pipeline Stage 完整参考（8 大模块 + 细粒度子步骤）

> **📌 文档定位：本文件是「每个 stage 干什么」的权威** —— 各 stage 的输入 / 处理 / 输出 /
> 参数 / 依赖 / checkpoint 逐项细节以本文件为准。
> 「怎么跑」（CLI 参数、config 分区、输出目录、checkpoint 规则）见 `RUN_GUIDE.md`；
> 模块归属与数据流见 `PIPELINE_FLOW.md`；论文方法学表述见 `METHODS.md`；文档总索引见 `README.md`。
> 机器校验：`_consistency_check_20260916/verify_docs_vs_code.py`。
>
> **顶层 8 大模块**（`--stage` 分组，自动展开）+
> **细粒度子步骤**（每个可单独运行，`--stage <子步骤>` 自动补依赖）。

> ⚠ **已废弃 stage（2026-08-27）**：`saturation` / `rtt` / `temporal` / `genes` / `concat`
> 这 5 个 stage 已从 `STAGE_ORDER` 移除，不再参与自动编排；单独调用仍可运行（保留兼容），
> 但**不会**被任何大模块自动展开。下文对应章节保留作历史参考，标题已加废弃标注。
> 现行 `STAGE_ORDER` 共 18 个（脚本清点，见附录 B）：
> `prep → metadata → online → clean → align → splitstree → rdp5 → tree → popgen → host
> → capheine → clock → beast → gene_dating → phylogeo → geo → geo_paths → report`
> （其中 `splitstree` 为新增的可选网络分析，详见 §6）

> **同步状态**：本文件已于 **2026-09-16** 与代码逐项复核并修正 9 处过期内容
> （rdp5 默认策略/路径、tree 命令与路径与支持值、popgen 路径、phylogeo checkpoint 与
> TempMig/RRT、geo 的 RRT 两种口径、report 路径与收尾语义；tree 一节于同日**再次**更新 ——
> 补支持值后「无自举」表述已作废）。
> 备份：`STAGE_REFERENCE.md.bak_stagesync_20260916`。
>
> ⚠ **两处同名死代码**（`2026-09-15` 审查 P3 标注，勿改错地方）：
> `run_stage_rtt()` 与 `run_stage_beast()` **当前无任何调用点**，能力已分别并入
> `clock`（实现在 `utils/clock_analysis.py`）与 `run_phase_beast()`。
> 改它们**不会**影响在跑的管线；真正在用的入口是 `process_virus()` 的 stage 分发。

## 8 大模块（--stage 顶层名）

> 下表与 `phylo_pipeline.py:174` 的 `STAGE_GROUPS` **逐字对应**。
> 注意：`metadata` **不在**任何大模块展开清单里（`data` 组只含 `prep/online/clean`），
> 因为它由 `--metadata` 入口参数驱动、并在 `prep` 之后自动调用；需要单独重跑时直接写 `--stage metadata`。
> 同理 `rdp5` 只在 `recomb` 组，不在 `phylogeny` 组。

| 大模块 | 展开的子步骤 | 对应你的骨架 |
|---|---|---|
| `data` | prep + online + clean | 数据准备（Phase 1/2/3） |
| `popgen` | popgen + host | 群体遗传学（Phase P，变异/宿主分化） |
| `recomb` | rdp5 | 重组检测（Phase R） |
| `select` | capheine | 选择压力（Phase S） |
| `phylogeny` | align + splitstree + tree | 序列比对+建树（Phase 3） |
| `time` | clock + beast + gene_dating | 时间信号+定年（Phase 0a/0b/5/G） |
| `geography` | phylogeo + geo + geo_paths | 系统地理（Phase 6，含遗传-地理相关性） |
| `report` | report | 可视化+报告（Phase 7） |

```bash
# 大模块方式 (推荐)
python phylo_pipeline.py --extract_dir 06_extraction/ --metadata meta.tsv --stage data,phylogeny,time,geography

# 子步骤方式 (精细控制, 每个都能单独跑)
python phylo_pipeline.py --work_dir gcva/ --alignment mafft.aln.fasta --tree iqtree.treefile --meta dates.csv --stage beast
```

## 目录结构（按 8 大模块）

```
<病毒>/                              ← 病毒 work_dir
├── data/                           ← data 模块（prep + online + clean）
│   ├── combined.fasta               ← prep 样本+参考序列
│   ├── dates.csv / sample_metadata.csv
│   ├── ncbi_ref/                    ← online efetch 单参考
│   ├── ncbi_harvest/                ← online SeqHarvester 公共序列
│   ├── metadata/                    ← metadata（入口治理：metadata_std.csv + 报告）
│   └── clean/                       ← clean（clean.fasta + removed.fasta + clean_report.tsv）
├── popgen/                         ← popgen 模块（popgen + host）
│   ├── popgen_report.txt
│   └── host_analysis/
├── recomb/                         ← recomb 模块（rdp5）
│   └── rdp5/（事件 CSV + masked/masked_N.fasta）
├── select/                         ← select 模块（capheine）
│   └── capheine/（cawlign/iqtree/hyphy/drhip/multiqc）
├── phylogeny/                      ← phylogeny 模块（align + splitstree + tree）
│   ├── mafft.aln.fasta              ← align
│   ├── saturation_out/              ← align 内置饱和分析（Iss）
│   ├── splitstree/                  ← splitstree（NeighborNet 网络）
│   └── iqtree.treefile (+ .log/.ckp.gz 等)  ← tree
├── time/                           ← time 模块（clock + beast + gene_dating）
│   ├── clock/                       ← clock（整合 rtt + temporal + genes 三个旧子步骤）
│   │   ├── treetime_rtt/ treedater_ltt/ ← 旧 rtt
│   │   ├── temporal_signal/ bets/       ← 旧 temporal
│   │   └── gene_analysis/               ← 旧 genes（分基因 RTT）
│   ├── beast/（beast1.xml/log/trees/beauti/postprocess）← beast（BEAST 1.x 全长定年）
│   └── gene_dating/                 ← gene_dating（分基因 BEAST）
├── geography/                      ← geography 模块（phylogeo+geo）
│   ├── phylogeography/ + visualization/  ← phylogeo
│   └── geo_analysis/ + prior_comparison/ ← geo
├── report/                         ← report 模块
│   ├── phylo_summary.csv
│   └── phylo_report.html
└── .checkpoints/                   ← 断点续传（隐藏，不归模块）
```

> 输入文件（work_dir 模式的 --alignment/--tree/--meta、GBK 注释）在病毒目录根，不归模块。

---

## 细粒度子步骤参考（每个可单独运行）

## 1. `prep` — 数据收集

| 项 | 内容 |
|---|---|
| **输入** | `--extract_dir`（06_extraction 目录）、`--metadata`（Core14.tsv）、`--info_dir`、`--variants_dir`（参考序列） |
| **处理** | load_sample_metadata（Core14 + Merged 双源）→ scan_extraction_dir（按病毒整理）→ prepare_virus_inputs（样本×元数据匹配） |
| **输出** | `<virus>/combined.fasta`、`dates.csv`、`sample_metadata.csv` |
| **参数** | `--min_samples`(5)、`--max_viruses`(0) |
| **依赖** | 无（batch 模式起始） |
| **Checkpoint** | 不参与（每次重新扫描构建队列） |
| **单独运行** | `python phylo_pipeline.py --extract_dir 06_extraction/ --metadata meta.tsv --stage prep` |

### 去重闸（datasets.yaml 单病毒模式，2026-09-03）

单病毒模式在 prep 构建后、process_virus 前执行 `dedup_triplet()`（utils/data_collector.py）：
本地+线上合并后的 fasta 按 **(date, location, 序列完全一致)** 三元组去除克隆重复；跨时空相同序列保留（传播事件有信息量）。参考序列（description 含 reference）优先保留。

- 产物：`data/combined_dedup.fasta` + `data/dedup_report.csv`（removed_id/kept_id/date/location/reason），无重复时不写文件、fasta 指回原文件
- Checkpoint：`dedup_report.csv` 存在则跳过
- 保护：`phylogeny/mafft.aln.fasta` 已存在时只警告不替换（避免下游 checkpoint 失配）
- 参数：`--dedup`（默认开）/ `--no-dedup` 关闭
- 实测：PSTVd 134 条零克隆重复；合成数据（同 date+loc+seq）精确去除并出报告

---

## 2. `metadata` — 入口 metadata 治理（时间 + 地理 的检查与矫正）

> 2026-09-16 新增（老师要求「时间和地理的检查和矫正，补充到我们的流程中对 metadata 信息
> 进行格式化，以免下游错误」）。位于**全流程最上游**，`online` 之前。

| 项 | 内容 |
|---|---|
| **输入** | `--metadata`（原始 metadata CSV/TSV，多来源混杂） |
| **处理** | `utils/metadata_governance.py`：日期结构化解析（格式白名单）→ 精度分档 → 小数年；地理层级拆分 → `Unknown` 层删除 → 编码归一（冒号/逗号）；占位符识别与处置 |
| **输出** | `<virus>/data/metadata/`：`metadata_std.csv`（下游唯一输入）、`metadata_report.tsv`（逐行体检 18 列）、`metadata_summary.txt`（汇总） |
| **参数** | `--metadata_mode`(mid)、`--metadata_date_col`/`--metadata_loc_col`/`--metadata_name_col`、`--metadata_keep_unknown_geo`、`--metadata_no_derive_coords`、`--metadata_no_blank_placeholders`、`--metadata_drop_placeholder_rows`、`--metadata_drop_no_date`、`--metadata_drop_no_location` |
| **依赖** | prep |
| **Checkpoint** | `data/metadata/metadata_std.csv` 存在 |
| **单独运行** | `--stage metadata` |

### 2.1 日期：格式白名单 + 精度分档

| 格式 | 例 | 归一化 | precision |
|---|---|---|---|
| ISO 全 | `2024-04-24` | 原样 | day |
| ISO 到月 | `2019-06` | 补 `-15` | month |
| 仅年 | `2024` | 保留（小数年取年中/年初，由 mode 定） | year |
| GenBank | `20-Jun-2019` | → `2019-06-20` | day |
| 欧式 | `15/07/2016`、`15.07.2016` | → `2016-07-15` | day（`01/07` 类日月皆 ≤12 时标 `ambiguous`） |
| 已是小数年 | `2007.569473` | 反解为日期 | day（`consumed_as=already_decimal`） |

- **分隔符归一化**：吸收 BioAider 的 `re.split("[./_-]")`，`. / _ -` 四种混用都能解析
  （如 `2014_05`、`2022.11-11`）。**这是从 BioAider 唯一吸收的一点**——它其余行为
  （脏值静默回显、整百年闰年判错、月份越界静默强解）都比我方弱，未采用。
- **小数年消歧判据**：小数位 **≥4 位** → 小数年（VirPhyKit 实测 431×6位 / 42×5位 /
  3×4位）；**恰好 2 位** → 年月（`2014.05`）。避免重复转换。
- **两条口径**：`mode='mid'`（月中/年中，同 `utils/decimal_year`，默认）；
  `mode='start'`（月初/年初，同 `virphy_bridge` / VirPhyKit）。**逐位对拍验证**。

### 2.2 地理：层级校验 + `Unknown` 删除

- **冒号也拆**：`China:Ningxia` → `['China','Ningxia']`（补 `geo_resolver._norm` 只拆
  `[,;]` 的缺口，原先靠 L3 段内反查侥幸命中）
- **`Unknown` 层级删除**（老师拍板）：`China, Unknown, Yinchuan_AI` →
  `China, Yinchuan_AI`，删除项记入 `dropped_levels` 并在报告留痕。
  覆盖 `Unknown_AI` 这类「占位词 + 来源标记」写法。
  > 与 `geo_resolver` 的区别：后者 `skip` 掉 `Unknown` 段后**仍返回国家质心并标为
  > `'province'`**，下游无从分辨；本模块删除该层并显式标 `unknown_level`。

### 2.3 占位符处置（2026-09-16 老师拍板：**坐标先补算，再清空**）

| 情形 | 默认行为 |
|---|---|
| 日期列占位符（字面量 `YYYY-MM-DD`、`NA` 等） | **清空该字段**，行保留 |
| 地名列占位符（`Country:Region`） | **清空该字段**，行保留 |
| 坐标列占位符（`XX.XX N XXX.XX E`） | **先用地名补算** → 成功则填入；推不出才清空 |

- `--metadata_no_blank_placeholders`：不清空，保留占位符原值
- `--metadata_no_derive_coords`：**不补算**坐标，回到"直接清空"旧行为
  （补算本身默认开启，覆盖**全球**地名，见 §2.3.1）
- `--metadata_drop_placeholder_rows`：含占位符字段的行**整行剔除**
  （坐标占位符**不触发**剔除——否则 450/450 全占位符的坐标列会把整表剔光）
- 全部剔除时 `kept=0` → 编排层判定失败并**降级使用原 metadata**，不阻断流程
- **幂等**：输入已是 `metadata_std.csv` 时直接复用，不二次处理
- 原件**完全不动**（另存新文件）

#### 2.3.1 坐标补算（`coord_from_location`）—— 全球多级查找

**2026-09-16 第二轮**（老师：「不能只是中国的省市呀，全世界的地理怎么做」）：
由「30 个中国城市硬编码字典」升级为**多级全球查找**。

查找优先级（高 → 低）：

| 级 | 来源 | 覆盖 | `coord_source` | 误差量级 |
|---|---|---|---|---|
| 1 | `_CITY_COORDS` 手校字典 | 30 个中国城市 | `city` | 约几十公里 |
| 2 | 中国中心点兜底（仅「China 无更细层级」） | 1 | `country_centroid` | 上千公里 |
| 3 | 全球 admin1 省/州 | **4410** 条 | `admin1:<ISO>` | 几十~几百公里 |
| 4 | 全球 admin0 国家 | **255** 条（236 国） | `country:<ISO>` | 数百~上千公里 |
| 5 | 推不出 | — | `''` | 清空 |

**数据来源**（可复跑、可审计）：

- `_geo_build_20260916/dist/gazetteer.json`，由 **Natural Earth 1:10m**
  （**public domain**，`naturalearthdata.com/about/terms-of-use`）构建
- 构建脚本 `_geo_build_20260916/scripts/build_gazetteer.py`，
  用 shapely `representative_point()` 取几何**内部**点
  （不用 centroid：凹多边形/多部件会落到陆地外）
- 参考实现：`D:\桌面\植物病毒分析平台\git-repo\spreadgl2.github.io`（MIT）的
  `src/lib/format/gazetteer.ts` + 其公开构建仓库 `spreadgl2/spreadgl2-gazetteer`
- 输出确定性 + 源文件 `sha256` 记于 `dist/sources.json`

**支持的地名写法**：`China:Ningxia`、`China, Ningxia, Yinchuan_AI`、
`United States:Maryland`、`Macao, China`、`China > Ningxia`、中文名（`宁夏`/`中国台湾`）。

**同名消歧**：`Maryland` 同属美国(`US-MD`)与利比里亚(`LR-MY`)，
按国家前缀选对（`United States:Maryland` → `US-MD`）。
Natural Earth 里共 **95 组**重名 admin1，故消歧是必需而非可选。

> ⚠ **中国主权相关（重要）**：上游 Natural Earth 把台湾标为
> `Sovereign country`、中文名写作「中华民国」，与中华人民共和国官方立场不一致。
> 构建脚本在生成时**已规范化**：台湾/香港/澳门 一律归入 `countryIso2='CN'`，
> 规范名 `Taiwan, China` / `Hong Kong, China` / `Macao, China`，
> 中文名 `中国台湾` / `中国香港` / `中国澳门`，**kind 由 country 降为 admin1**。
> 解析层另设 `_CN_SUBREGION_KEYS`，确保港澳台判为「比国家更细的层级」，
> 不会输出国家层级表述。

> ⚠ 补算值是**推算约值**，不是原始观测。`coord_source` 必须一并消费。

实测：

| 文件 | 坐标占位符 | 补算成功 | 推不出 | 覆盖率 |
|---|---|---|---|---|
| `unified_metadata.csv` | 450 | **408**（city 406 / centroid 2） | 42 | **0% → 90.7%** |
| `unified_metadata_barbarum_real.csv` | 9 | **8**（city 5 / US-MD×2 / US-MA×1） | 1 | — |

推不出的 42 行：地名恰好也是 `Country:Region`（与 42 条日期占位符**同一批行**）。
barbarum 推不出的 1 行：`Unknown:Unknown`（确无地理信息）。

**与旧行为的关系**（2026-09-16 老师逐条确认）：

| 地名 | 旧 `geo_to_latlon` | 新实现 | 说明 |
|---|---|---|---|
| 30 个中国城市 | ✓ | **逐字一致** | 回归硬约束 |
| `China`（无更细层级） | 中国中心点 | **逐字一致** | 零漂移 |
| `China:Gansu` / `:Qinghai` | 中国中心点（误差上千公里） | **真实省坐标** | 老师选「精度提升」 |
| `United States:Maryland` / `:Massachusetts` | **推不出** | 能算 | 老师选「采用」 |

全库对拍：原实现可解析 **479** → 新实现 **482**，**一条未丢**，
差异仅限上述白名单 4 个地名。

### 2.4 实测（真实数据）

| 文件 | 行数 | 结果 |
|---|---|---|
| `unified_metadata.csv` | 450 | 日期 day=408 / placeholder=42；**坐标补算 408**（推不出 42）；删 Unknown 层级 105 行 |
| `unified_metadata_barbarum_real.csv` | 75 | 日期 day=36 / year=39；**坐标补算 8**（余 1 条 `Unknown:Unknown` 推不出）；删层留痕 |

---

## 3. `online` — NCBI 参考补充

| 项 | 内容 |
|---|---|
| **输入** | 病毒名（accession + 物种名） |
| **处理** | efetch 单参考 GBK + SeqHarvester 按物种名批量收集公共序列 + seq_grouper 分组图 |
| **输出** | `<virus>/ncbi_ref/`（单参考）、`ncbi_harvest/`（公共序列 + metadata + 分组图） |
| **参数** | 无专属（网络依赖） |
| **依赖** | prep、metadata |
| **Checkpoint** | `ncbi_ref` 或 `ncbi_harvest` 目录存在 |
| **单独运行** | `--stage online` |

---

## 4. `clean` — 比对前序列清洗（第一道过滤）

> 2026-09-16 新增（老师要求「分析前，再次长度检查，去除长度异常，N 碱基异常等，
> 修正大小写，去除时间和地理都没记录的序列」）。在 `align` **之前**，与比对后的
> `align_qc` 互补。

| 项 | 内容 |
|---|---|
| **输入** | prep 的 `combined.fasta`（或 work-dir 模式 `--alignment`） |
| **处理** | `utils/seq_clean.py`：大小写归一 / 非法字符 / gap 字符 / 长度上下限 / N+简并码占比 / RNA(U) / 元数据缺失 |
| **输出** | `<virus>/data/clean/`：`clean.fasta`（下游比对输入）、`removed.fasta`（被剔除序列，可追溯）、`clean_report.tsv` |
| **参数** | `--clean_min_length_ratio`(0.9)、`--clean_max_length_ratio`(1.5)、`--clean_max_n`(0.05)、`--clean_reference`、`--clean_allow_rna`、`--clean_strip_illegal`、`--clean_allow_gap_chars`、`--clean_metadata`、`--clean_drop_no_date`、`--clean_drop_no_location` |
| **依赖** | prep |
| **Checkpoint** | `data/clean/clean.fasta` 存在 |
| **单独运行** | `--stage clean` |

### 4.1 与 `align_qc` 的分工

| | `clean`（本 stage） | `align_qc` |
|---|---|---|
| **时机** | 比对**前** | 比对**后** |
| **看什么** | 序列本体（大小写、非法字符、长度、N） | 需比对坐标系的指标（gap%、identity、比对坐标系 N%） |
| **为何必须在 MAFFT 前** | 空格/星号/小写会让 MAFFT 报错或产出**错位比对**，事后 QC 只能发现问题、不能修复比对本身 | — |

- 编排放行：`prep['combined_fasta']` 被改写为 `clean.fasta`，原值存 `_raw_combined_fasta`
- **全部剔除时** `success=False` 并**降级用原序列**，不阻断流程
- 实测：12 条脏数据对照，`align_qc` 放行 10 条，`clean` 只放行 5 条（多拦住 5 类本体问题）

---

## 5. `align` — MAFFT 全长比对

| 项 | 内容 |
|---|---|
| **输入** | `clean` 的 `clean.fasta`（work-dir 模式用 `--alignment`） |
| **处理** | `mafft --auto --thread N --reorder`，并内联 saturation（纯 Python，秒级，已并入本 stage） |
| **输出** | `<virus>/phylogeny/mafft.aln.fasta` |
| **参数** | `--skip_msa`、`--threads` |
| **依赖** | prep、clean |
| **Checkpoint** | `mafft.aln.fasta` 存在且非空 |
| **单独运行** | `--stage align` |

---

## 6. `splitstree` — 分裂网络（NeighborNet，树状/网状信号判别）

| 项 | 内容 |
|---|---|
| **输入** | `align` 的 `mafft.aln.fasta`（work-dir 模式用 `--alignment`） |
| **处理** | `utils/splitstree_bridge.py` → SplitsTree6 命令行，`--network neighbornet` 构建分裂网络（split decomposition） |
| **输出** | `<virus>/phylogeny/splitstree/neighbornet.nexus`（+ 内部 nsplits 指标） |
| **参数** | `--splitstree_max_taxa`(200) |
| **依赖** | align |
| **Checkpoint** | `phylogeny/splitstree/neighbornet.nexus` 存在 |
| **单独运行** | `--stage splitstree` |
| **定位** | 与 `tree` 并列：`tree` 强制给出**一棵树**（隐含"数据是树状"假设），`splitstree` 不假设树状，用分裂图暴露**冲突信号**（重组、不完全谱系分选）。两者互不依赖，可同时开 |
| **降级** | 无比对 / bridge 不可用 / 未产出 → 记 warning 跳过，不阻断流程 |

> 判读：NeighborNet 图上出现**箱形（box）结构**而非干净树状分支时，提示存在网状进化——此时应回看 `rdp5` 的重组检测结果。

---

## 7. `rdp5` — RDP5 重组检测 + mask

| 项 | 内容 |
|---|---|
| **输入** | mafft.aln.fasta（+ 可选 --rdp5_genes 基因注释） |
| **处理** | `recombination_analysis`（RDP5 7 方法）→ **默认 `drop_recombinants`：删整条重组序列后重新 MAFFT**，后续 tree/clock/beast 全部吃重比对；`--rdp5_mask` 时走旧策略（重组区置 N，**不重比对**，保样本量） |
| **输出** | `recomb/rdp5/<virus>.csv` + 默认 `recomb/rdp5/rdp5_dropped.realn.fasta`（干净比对喂 tree）；`--rdp5_mask` 时为 `recomb/rdp5/masked/masked_N.fasta` |
| **参数** | `--rdp5`、`--rdp5_script`、`--rdp5_genes`、`--mask_min_methods`(3)、`--rdp5_stage`(run\|all)、`--rdp5_mask` |
| **依赖** | align |
| **Checkpoint** | `recomb/rdp5/rdp5_dropped.realn.fasta`（`--rdp5_mask` 时为 `recomb/rdp5/masked/masked_N.fasta`） |
| **单独运行** | `--stage rdp5` |
| **短序列防御** (2026-09-02) | 比对 <500 bp（如类病毒 ~359nt）RDP5 窗口方法无足够位点，记 skip 返回原比对，不再空跑 11 分钟后 fail |

---

## 8. `tree` — IQ-TREE 建树

| 项 | 内容 |
|---|---|
| **输入** | `phylogeny/mafft.aln.fasta`，或 rdp5 开启时其 drop/mask 后的比对 |
| **处理** | `iqtree2 -s <aln> -pre <phylogeny>/iqtree -nt N -m MFP -B 1000 -alrt 1000 --quiet`（ModelFinder 自动选模型 + **UFBoot + SH-aLRT 双支持值**） |
| **输出** | `<virus>/phylogeny/iqtree.treefile`（+ 同名 `.iqtree` / `.log` / `.ckp.gz`）+ `phylogeny/iqtree.support.json`（支持值元数据 sidecar） |
| **参数** | `--skip_tree`、`--iqtree_bin`、`--iqtree_boot`、`--iqtree_alrt`、`--iqtree_bnni`、`--iqtree_timeout`、`--threads` |
| **依赖** | align |
| **Checkpoint** | `phylogeny/iqtree.treefile`；**输入指纹含 IQ-TREE 参数 token**（`iqtree_param_token()`），建树参数一变即重跑 |
| **单独运行** | `--stage tree` |

> ✅ **分支支持值已补齐（2026-09-16）**：默认拼 `-B 1000 -alrt 1000`，与平台侧
> `Virus_Platform_Core/phylo.py::_run_iqtree` 同口径。此前为控耗时写成 `-m MFP --quiet`，
> 主树无任何支持值 —— **该「不跑自举」表述自本日作废**（本文件头已同步声明）。

**三段降级阶梯**（`run_iqtree()`，`phylo_pipeline.py:634-759`）：

| 次序 | 实际旗标 | 说明 |
|---|---|---|
| ① `full` | `-B <boot> -alrt <alrt>` [ `--bnni` ] | 期望路径（双支持值） |
| ② `alrt_only` | `-alrt <alrt>` | 仅当 `alrt >= 1000` 才挂；UFBoot 在小比对/近缘序列上更易失败 |
| ③ `no_support` | （空）→ 等价 `-m MFP --quiet` | 最终兜底；此时树**无支持值**，报告不得声称有 |

> ⚠️ **`-B` / `-alrt` 在 IQ-TREE 里硬性要求 `>= 1000`** —— 想"少跑几轮"来降级是**不行的**
> （会被直接拒绝，不会自动打折）。要关闭支持值须把对应项写 **0**（两项都 0 = 回到修复前行为）。
> `--iqtree_boot`/`--iqtree_alrt` 传小于 1000 的值会被 `_iqtree_support_flags()` 静默过滤。

**结果取证**：每次建树都写 `phylogeny/iqtree.support.json`，字段含 `requested` / `actual`
（`flags`、`support_type`、实际生效的 `boot`/`alrt`）、`support_type`
（`SH-aLRT/UFBoot` / `UFBoot` / `SH-aLRT` / `bootstrap` / `none`）、`degraded` / `fallback`
（含原因）、`attempts`（每次尝试的 rc/err）、`n_internal` / `n_labeled`。
`support_type` 由**实际用过的旗标**判定（`support_type_from_flags()`），优于从树标签反猜。

> ⚠️ **解析已有树的支持值必须同时看 `clade.name` 与 `clade.confidence`**
> （`tree_support_stats()`）：Bio.Phylo 把**纯数字**标签放进 `confidence`、`name` 留 None，
> 带斜杠的 `95/100` 才进 `name`。只看 `name` 会把"只有单一支持值"的树误判成"无支持值"
> → 每次运行都白跑一遍 IQ-TREE。

> 仓库内唯一仍用旧旗标 `-bb 1000` 的地方是 `validate_recombinants.py:178`
> （重组验证的子区段建树），与主树无关；若要口径统一，改这里。

---

## 9. `popgen` — 群体遗传（pypopart）

| 项 | 内容 |
|---|---|
| **输入** | mafft.aln.fasta + sample_metadata.csv |
| **处理** | `popgen_analysis`：π/θ/Tajima D/Fu's Fs/Fst/MJN 网络 |
| **输出** | `<virus>/popgen/popgen_report.txt`（stdout 捕获）+ metrics；下游另有 `popgen/distance/`（距离矩阵+NJ树）、`popgen/dnasp/`、可选 `popgen/haplo`/`popgen/figures` |
| **参数** | `--popgen_group`(location)、`--popgen_min_n`(5)、`--popgen_perm`(999) |
| **依赖** | align |
| **Checkpoint** | `popgen/popgen_report.txt` |
| **单独运行** | `--stage popgen` |

---

## 10. `capheine` — 正选择分析（Phase S）

| 项 | 内容 |
|---|---|
| **输入** | GBK 注释（work_dir *.gbk，自动 gbk_extractor 提取 CDS 参考）+ 序列 |
| **处理** | capheine_pipeline：CAWLign 密码子比对 → IQ-TREE 基因树 → HyPhy FEL/MEME/PRIME/BUSTED/CONTRASTFEL/RELAX → DRHIP → MultiQC |
| **输出** | `capheine/`（cawlign/iqtree/hyphy/drhip/multiqc） |
| **参数** | `--capheine_ref`、`--capheine_unaligned`、`--capheine_code`(1)、`--capheine_workers`(4)、`--capheine_cpus_iqtree`(6)、`--capheine_cpus_hyphy`(32)、`--capheine_mpi` |
| **依赖** | prep |
| **Checkpoint** | `capheine/drhip` |
| **单独运行** | `--stage capheine` |
| **viroid skip** (2026-09-02) | GBK 提不出 CDS 时查文件头 CDS 特征：无 CDS（类病毒）记 skip（`参考 GBK 无 CDS (类病毒?), 记 skip`），有 CDS 但提取失败才 fail。注意 run_stage_capheine 内有**两条** GBK 提取路径（前置 datasets annotation_gb 与后置 gbk 搜索），两处均已接入；外层归账读 `cr['skip_reason']` 决定 skip/fail |
| **MPI 默认开** (2026-08-31) | `--capheine_mpi` 默认 True（BooleanOptionalAction），FEL/MEME/PRIME/CONTRASTFEL 走 mpirun；BUSTED/RELAX 不走 |
| **LPT 调度** (2026-09-03) | 基因按参考长度降序提交（长基因优先占槽，如 L 先跑），短基因填充空档；同基因内 FEL/MEME/PRIME/BUSTED（±CONTRASTFEL/RELAX）无相互依赖，全部并行 |

---

## 11. `clock` — 分子钟联合分析（整体 + 分基因 RTT/LTT + DRT）

> 本 stage 由 2026-08-27 合并旧 `rtt` / `temporal` / `genes` 而来（见附录 A），
> 是全流程**唯一**的时间信号/分子钟入口。

| 项 | 内容 |
|---|---|
| **输入** | `mafft.aln.fasta` + `iqtree.treefile` + `dates.csv`（+ 可选 `--meta` 供 BETS、`--gbk` 供分基因拆分） |
| **处理** | 实现在 `utils/clock_analysis.py`（`run_clock_analysis`，支持独立命令行）；本 stage 只做包装（找 gbk/dates/工具路径 + 汇总 metrics）。内含：① 整体 TreeTime-RTT（β/R²/p）② TreeDater-LTT 定年 ③ 分基因 RTT ④ 可选 DRT 时间信号 |
| **输出** | `time/clock/treetime_rtt/`、`time/clock/treedater_ltt/`、`time/clock/gene_analysis/`、可选 `time/clock/temporal_signal/` |
| **参数** | `--drt_randomizations`(10)、`--threads`、`--iqtree_bin`、`--rscript` |
| **依赖** | tree |
| **Checkpoint** | `time/treedater_ltt/Phylogeny_dated.pdf`（STAGE_CHECKPOINT_FILES 用的旧路径） |
| **单独运行** | `--stage clock` |
| **失败语义** | `success = bool(cr.get("summary_rows"))` —— 无任何汇总行即视为失败 |
| **metrics 过滤** | 只回填非 `None` 的 metrics（`{k:v for ... if v is not None}`），避免空值污染报告 |
| **子步骤记账** | `summary_rows` 里 `level == "gene"` 的行回填到 `results["gene_results"]`，供 `gene_dating` 取用 |

---

## 12. `beast` — 全长 BEAST 定年

| 项 | 内容 |
|---|---|
| **输入** | mafft.aln.fasta + dates.csv |
| **处理** | generate_beast1_phylogeo_xml（BEAST 1.x，UCLN + skyline，纯分子钟无 location）→ beauti 输入 → run_beast（--skip_beast_run 只出 XML）→ postprocess（ESS 检查，BSP/RSPP 降级）→ parse_beast_log 出全长 TMRCA（median/HPD/ESS，喂 gene_dating --full-tmrca） |
| **输出** | `beast/beast1.xml` + `beast1.log` + `beast1.trees` + `postprocess/` |
| **参数** | `--beast`、`--beast_bin`、`--beast_chain`(5M)、`--beast_burnin`(10)、`--beast_timeout`(48h)、`--skip_beast_run` |
| **依赖** | tree |
| **Checkpoint** | `beast/beast1.xml` |
| **单独运行** | `--stage beast --skip_beast_run`（先出 XML） |
| **注意** | 服务器只有 BEAST 1.10.4；XML 必须是 BEAST 1.x 格式（BEAST 2.7 XML 会被 1.10.4 拒收） |
| **脏日期防御** (2026-09-02) | 构建样本列表时预解析日期，解析失败（如 `not applicable`）的样本剔除并警告（`⚠️ 剔除样本 SRRxxx: 脏日期`），不再让整个 stage 崩溃。beast_bridge（BEAST2 路径）与 beast1_bridge 同款逻辑 |
| **多链闭环** (2026-08-31) | 多链模式前台等待：60s 轮询链 state，all_done 自动 merge；链停滞+无进程或 30 分钟无进展退回 pending。config.yaml `chains: 5` 为以后默认 |
| **ESS 算法对齐 Tracer** (2026-08-27) | utils/ess.py 改为 Geyer initial POSITIVE sequence（beast-mcmc TraceCorrelation 逐行一致：有偏自协方差 + 成对和>0 截断，maxLag 2000，无单调性检查；`monotone=True` 可复现旧敏感性口径）。trace.jar 官方库逐参数对账通过（残余偏差<15% 来自 burnin 取整）。历史：v1 逐 lag 截断严重低估（GCVA 50M run posterior 实际 49.1 被报 6.9）；v2 单调版对 acf 零点抖动序列虚高（PSTVd posterior 85 vs Tracer 14.6）。GCVA/PSTVd 两份 ess_results.csv 已按 v3 重算（旧版备份 .bak_v2_20260827） |

---

## 13. `gene_dating` — 分基因 BEAST 定年

| 项 | 内容 |
|---|---|
| **输入** | capheine 产物（--genes_dir 或自动找 `<virus>/capheine/cawlign`）+ meta |
| **处理** | gene_partition_dating：基因命名适配 → 逐基因 BEAST XML → setsid 后台提交 → watchdog 解析 |
| **输出** | `gene_dating/`（gene_split/ + 每基因 XML + watchdog） |
| **参数** | `--gene_dating`、`--genes_dir`、`--gene_chain`(15M)、`--gene_chain_big`(10M)、`--big_genes`(L) |
| **依赖** | capheine 产物（外置） |
| **Checkpoint** | `gene_dating/gcva_watch_parse.sh` |
| **单独运行** | `--stage gene_dating --genes_dir ...` |
| **空目录防御** (2026-09-02) | `_count_gene_aligns()` 验证 genes-dir 内比对文件数；空目录（viroid 无 CDS 或 capheine 未产出）记 skip 不写 checkpoint，不再空部署假 done |

---

## 14. `phylogeo` — 系统地理（CTMC+BSSVS）

| 项 | 内容 |
|---|---|
| **输入** | mafft.aln.fasta + meta（含 location） |
| **处理** | ① 生成 BEAST 离散性状 CTMC+BSSVS XML（`--phylogeo_prior auto` 时先做 constant vs skyline 先验比较）② 多链 BEAST 提交（`phylogeo_sky{i}.xml`，确定性种子）③ 链完成后 logcombiner + Python burnin 过滤 + treeannotator → **MCC 树** ④ TempMig 逐年迁移矩阵 ⑤ SpreaD3 交互地图 + 迁移 Bayes factor |
| **输出** | `<virus>/geography/phylogeography/`：`phylogeo_beast1.xml`、`phylogeo_sky{i}.log/.trees`、`merged/{merged.log,merged.trees,mcc.tree}`、`visualization/`（SpreaD3 JSON + HTML + `migration_bf.csv`）；`<virus>/geography/tempmig/` |
| **参数** | `--phylogeo`、`--phylogeo_clock`(ucln)、`--phylogeo_prior`(skyline)、`--phylogeo_chain`(10M)、`--phylogeo_chains`(3)、`--phylogeo_bf`(5)、`--beast_version`(1)、`--skip_phylogeo_run`、`--treeannotator_bin` |
| **依赖** | tree |
| **Checkpoint** | `geography/phylogeography/merged/mcc.tree`（链还在后台跑时**不写** checkpoint，下次重跑再检查合并 → pending 语义） |
| **单独运行** | `--stage phylogeo` |
| **config 对齐坑** (2026-09-16) | `config.yaml` 里的 `phylogeo.chains` 原先**不在** `CONFIG_SECTION_ARG['phylogeo']` 映射里 → 被静默忽略、实际一直用 argparse 默认 3。现已登记为 `phylogeo_chains`。**加 config 键必须同步登记映射，否则静默失效** |

---

## 15. `geo` — 地理分析

| 项 | 内容 |
|---|---|
| **输入** | mafft.aln.fasta + iqtree.treefile + meta（含 location） |
| **处理** | ① 可选地理均衡抽样（`--geo_subsample N`，只写子集不覆盖全量）② Mantel 检验（遗传×地理）③ 树地理着色 ④ VirSpaceTime 时空图 ⑤ 后处理 `geo_analysis`：**RRT 距离矩阵法**（1000 次置换，固定种子）⑥ 可选 **RRT 树注释法**（VirPhyKit 同口径，默认关闭） |
| **输出** | `<virus>/geography/geo_analysis/`：`mantel_test.pdf`、`tree_geography.pdf`、`VirSpaceTime.pdf`、`extended/rrt_results.pdf`；可选 `extended/rrt_tree/{rrt_tree_table.csv,rrt_tree.pdf}` + `extended/rrt_tree_summary.json` |
| **参数** | `--geo_subsample`(0=不抽样)、`--rrt_randomized_dir`（**默认无 → 不跑树注释法**） |
| **依赖** | tree |
| **Checkpoint** | `geography/geo_analysis/mantel_test.pdf` |
| **RRT 两种口径**（2026-09-16 接线） | ⑤ 与 ⑥ 回答的是**两个不同问题**，不可互相引用：<br>· **距离矩阵法** `run_rrt`：组间遗传距离随机化 → "地理聚集是否显著"（不需 BEAST，一直在跑）<br>· **树注释法** `run_rrt_tree`：MCC 树根状态后验 vs 区域随机化 → "地理状态是否优于随机"（**= VirPhyKit RRT 口径**，已接线但**默认关闭**）<br>⚠️ 树注释法需要**自备 N 棵区域随机化 MCC 树**（上游口径 N=20 ⇒ ~21 次完整 BEAST 重跑）。<br>✅ **输入格式不是障碍**（2026-09-16 更正）：`列名.set` + `列名.set.prob` 由 **TreeAnnotator 汇总时自动生成**（对字符串型节点属性做频次统计），所以本管线 phylogeo 的 `merged/mcc.tree` 本来就带 `location.states.set`/`.set.prob`，**直接喂即可**；只有误传逐样本 `.trees` 才会看到"只有点状态"。详见 `UPSTREAM_CONSISTENCY_20260916.md` §5.4 |
| **单独运行** | `--stage geo` |
| **地名解析** (2026-09-02) | 坐标解析接入 `utils/geo_resolver.py` 分层解析器：L0 遗留 COORD_LOOKUP → L1 省/市中英文+拼音表 → L2 规范化段拆分 → L3 模糊 → L4 持久缓存 `cache/geo_cache.json` → L5 在线兑底（Photon 公共 API → Nominatim，默认开，`GEO_RESOLVER_ONLINE=0` 关）。解析失败仍诚实剔除（绝不伪造坐标），日志打印解析来源统计（如 `{'legacy': 132, 'province': 5}`） |
| **维度修复** (2026-09-02) | 坐标过滤后 geo_dist 与 gen_dist 均用新 n 重建（旧代码 geo_dist 残留旧 n 会 ValueError 维度不匹配） |

---

## 16. `geo_paths` — 事件级传播路径（2026-09-16 新增）

| 项 | 内容 |
|---|---|
| **输入** | phylogeo 的 `merged/mcc.tree`（吃 `location.states.set/.set.prob` 后验）+ meta（name,location）+ mafft.aln.fasta（推 seq_len） |
| **处理** | ① `utils/transmission_paths.py::run_pathways`：逐分支时空表（seraphim `treeExtractions` 同构）→ 只留地点变化分支 = 事件级 pathway → 同地点连续节点折叠为 episode（**自根 DFS 拓扑序**，免疫 phymapr 按日期排序的同日期 NA bug）→ 三级置信分类 Direct/Indirect/Distant-Import（SNP 维度=期望替换数 `rate×Δheight×L`，**非实测**；阈值 2/5 沿用 phymapr，**论文前须按本数据 E[S] 分布标定**）→ LTL 持久性摘要（删失窗 0.5 年，polio-wpv1 同款）② `utils/diffusion_stats.py`：分支扩散速度/扩散系数/wavefront 时间序列 + **坐标重标注置换零模型**（种子 42；初版端点洗牌在单源网络下无检验力已废弃）③ 可选 `--geo_map`：`utils/transmission_map.py` 贝塞尔弧线+末段箭头+时间色带静态图（**默认纯矢量无第三方瓦片，合规**）+ 可选 `--geo_gif` 动画 |
| **输出** | `geography/geo_analysis/pathways/{branch_table,pathways,episodes,ltl_summary}.csv + pathway_qc.json`；`geography/geo_analysis/dispersion/{velocity_table,wavefront_series}.csv + dispersion_summary.json`、`geo_analysis/timescale_consistency.json`（双流一致性闸门：时间流 RTT/DRT vs 地理流根高/速率并排 + WARN 旗标，2026-09-17）；可选 `geography/visualization/transmission_map.png(+.pdf/.svg/.gif)` |
| **参数** | `--no_pathways`(关)、`--direct_snp`(2.0)、`--indirect_snp`(5.0)、`--censor_years`(0.5)、`--n_perm`(200)、`--geo_map`、`--geo_gif`、`--geo_gif_frames`(36)、`--geo_tree_map`(树-图联动双面板)、`--calibrate`(auto=E[S] 分位数 q25/q75 自标定, 默认经 config) |
| **依赖** | phylogeo（checkpoint 即其 MCC） |
| **Checkpoint** | `geography/geo_analysis/pathways/pathway_qc.json` |
| **单独运行** | `--stage geo_paths`（geography 组自动包含） |
| **出图默认** (2026-09-17) | config `geo_paths.map/tree_map: true` → 每次 run 自动产出 `visualization/transmission_map.*` 与 `transmission_tree_map.*`；GIF 仍按需 (`--geo_gif`) |
| **方法定位** | MCC **描述层**（与 RRT/TEMPMIG 并列），不改动 BSSVS BF 推断结论；置换零模型与 seraphim RRW 参数化零模型**不同**，不可互换引用。方法出处逐条见各模块 docstring 与 `INTEGRATION_20260916.md` |

---

## 17. `host` — 宿主分化

| 项 | 内容 |
|---|---|
| **输入** | mafft.aln.fasta + iqtree.treefile + meta（含 host/ScientificName 列） |
| **处理** | 宿主树着色 + 组间/组内遗传距离检验 |
| **输出** | `host_analysis/`（tree_host.pdf、host_diff.pdf） |
| **参数** | 无（自动） |
| **依赖** | tree |
| **Checkpoint** | `host_analysis` 目录 |
| **单独运行** | `--stage host` |

---

---

## 18. `report` — 汇总报告

| 项 | 内容 |
|---|---|
| **输入** | 各 stage 产物 + results 汇总 |
| **处理** | generate_summary（CSV，含 Stages_done/Skipped/Failed）+ build_html_report |
| **输出** | `report/phylo_summary.csv` + `report/phylo_report.html`（batch 模式落 `<output_dir>/report/`；单病毒模式落 `<data_dir>/report/`） |
| **参数** | 无 |
| **依赖** | 已完成 stage |
| **Checkpoint** | `phylo_summary.csv` |
| **单独运行** | `--stage report`（batch 模式需已有产物） |
| **收尾语义**（2026-09-16 澄清） | 汇总**不是** `process_virus` 里的 stage 分支，而是 `main()` 收尾时**无条件**产出 —— 这样 `--stage clock` 这类局部重跑也会刷新汇总表。开关用 config `report.enabled`（`--report_enabled` / `--no-report`），**不是** stage 清单 |

---

## 附录 A：已废弃 stage（2026-08-27，保留作历史参考）

> 以下 5 个 stage 已从 `STAGE_ORDER` 移除，**不参与自动编排、不被任何大模块展开**。
> 代码里 `_LEGACY_CLOCK_SUBSTAGES = ('rtt','temporal','genes')` 保留旧名字仅用于
> **checkpoint 兼容**（老数据目录仍能被识别）。新流程请用合并后的替代 stage。

| 废弃 stage | 合并去向 | 说明 |
|---|---|---|
| `saturation` | → `align` 内置 | 替换饱和分析（Iss）在 MAFFT 比对后**内联**执行，秒级 |
| `rtt` | → `clock` | TreeTime-RTT + TreeDater-LTT |
| `temporal` | → `clock` | DRT / BETS 时间信号 |
| `genes` | → `clock` | 分基因 TreeTime-RTT |
| `concat` | **移除** | 多基因拼接（当前不做 supermatrix，分基因各自建树/定年） |

### A.1 `saturation` — 替换饱和分析（数据质量）

| 项 | 内容 |
|---|---|
| **输入** | mafft.aln.fasta |
| **处理** | 逐对 Ti/Tv/p → C 值（YR-MPE 法）+ Xia Iss/Iss.c（全局频率随机重采样）；p<0.05 时以 C 值法为准 |
| **输出** | `saturation_out/saturation_report.csv` + `saturation_plot.png` |
| **参数** | `--sat_replicates`(1000) |
| **现状** | 已并入 `align`，产物仍落在 `phylogeny/saturation_out/` |

### A.2 `rtt` — TreeTime-RTT + TreeDater-LTT

| 项 | 内容 |
|---|---|
| **输入** | mafft.aln.fasta + iqtree.treefile + dates.csv |
| **处理** | TreeTime 根到梢回归（β/R²/p）→ TreeDater LTT 定年 |
| **输出** | `time/clock/treetime_rtt/`、`time/clock/treedater_ltt/` |
| **现状** | 已并入 `clock`；子步骤名 `rtt` 仅用于 checkpoint 兼容 |

### A.3 `temporal` — 时间信号（DRT/BETS）

| 项 | 内容 |
|---|---|
| **输入** | mafft.aln.fasta + iqtree.treefile + dates.csv（DRT）/ meta（BETS） |
| **处理** | DRT：真实 RTT vs N 次日期打乱；BETS：path sampling 异时/同时比较 |
| **输出** | `time/clock/temporal_signal/drt_results.pdf`、`time/clock/bets/` |
| **参数** | `--check_temporal`、`--drt_randomizations`(20)、`--bets`、`--bets_steps`(20)、`--bets_chain`(1M) |
| **现状** | 已并入 `clock`；BETS 需在同一次运行里带上 `--bets` |

### A.4 `genes` — 分基因 TreeTime-RTT

| 项 | 内容 |
|---|---|
| **输入** | 各基因比对 + 树 + dates |
| **处理** | 逐基因 RTT（与全长 RTT 同法） |
| **输出** | `time/clock/gene_analysis/` |
| **现状** | 已并入 `clock` |

### A.5 `concat` — 多基因拼接（超矩阵）

| 项 | 内容 |
|---|---|
| **输入** | 各基因比对 + 分区定义 |
| **处理** | 拼接 supermatrix + 分区文件 |
| **输出** | `phylogeny/concat_out/` |
| **现状** | **已彻底移除**。当前流程按基因分别建树、分别定年，不做超矩阵串联 |

---

## 附录 B：stage 速查表（与 `STAGE_ORDER` 逐字对应）

`STAGE_ORDER` 共 **18 项**（`phylo_pipeline.py:126-132`，逐字抄录并于 2026-09-16 用脚本清点）：

```text
STAGE_ORDER = [
    'prep', 'metadata', 'online', 'clean', 'align', 'splitstree', 'rdp5', 'tree',
    'popgen', 'host', 'capheine',
    'clock', 'beast', 'gene_dating',
    'phylogeo', 'geo', 'geo_paths',
    'report',
]   # len == 18
```

| `STAGE_ORDER` 序号 | stage | 本文档章节 |
|---|---|---|
| 1 | prep | §1 |
| 2 | metadata | §2 |
| 3 | online | §3 |
| 4 | clean | §4 |
| 5 | align | §5（饱和分析内联，见附录 A.1） |
| 6 | splitstree | §6 |
| 7 | rdp5 | §7 |
| 8 | tree | §8 |
| 9 | popgen | §9 |
| 10 | host | §17 |
| 11 | capheine | §10 |
| 12 | clock | §11（整合旧 rtt/temporal/genes） |
| 13 | beast | §12 |
| 14 | gene_dating | §13 |
| 15 | phylogeo | §14 |
| 16 | geo | §15 |
| 17 | geo_paths | §16 |
| 18 | report | §18 |

> ⚠ **章节号 ≠ 执行序**：本文档把 `capheine`/`clock` 提前（按"选择压力→分子钟"的逻辑
> 归组更易读），`host` 后置；`STAGE_ORDER` 里 `host` 紧跟 `popgen`。
> **执行顺序一律以 `STAGE_ORDER` 为准，本章节号仅供检索。**
>
> 校验方法（可复现）：`python -c "import re,io; s=io.open('phylo_pipeline.py',encoding='utf-8').read(); m=re.search(r'STAGE_ORDER = \[(.*?)\]',s,re.S); print(len(re.findall(r\"'(\w+)'\",m.group(1))))"` → `17`

---

## 分组别名（`--stage` 顶层名，与 `STAGE_GROUPS` 一致）

| 别名 | 展开 |
|---|---|
| `data` | prep, online, clean |
| `popgen` | popgen, host |
| `recomb` | rdp5 |
| `select` | capheine |
| `phylogeny` | align, splitstree, tree |
| `time` | clock, beast, gene_dating |
| `geography` | phylogeo, geo, geo_paths |
| `report` | report |

> `metadata` 不在任何别名里：它由 `--metadata` 入口参数驱动、在 `prep` 之后自动执行；
> 需单独重跑时写 `--stage metadata`。
