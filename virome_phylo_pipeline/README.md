# MMPV-RNA virome_phylo_pipeline — 病毒系统地理学 + 群体遗传学工具链

> 枸杞属病毒组（GCVA / PSTVd）分析沉淀。双轨道：BEAST 1.x 系统地理学 + pypopart 群体遗传学。

---

## 0. 文档权威（**先读这一节**）

六份文档各有**唯一职责**，交叉内容只保留一处，冲突时按下表判定：

| 文档 | 权威范围 | 一句话 |
|---|---|---|
| **`README.md`**（本文件） | **入口索引** + BEAST 踩坑 + 服务器路径 | 「从哪读起、哪里有坑」 |
| **`RUN_GUIDE.md`** | **CLI 参数 · config 分区 · 输出目录 · checkpoint** | 「怎么跑」唯一权威 |
| **`STAGE_REFERENCE.md`** | **18 个 stage 的逐项细节**（输入/处理/输出/参数/依赖/checkpoint） | 「每个 stage 干什么」唯一权威 |
| **`PIPELINE_FLOW.md`** | **模块归属总图 / 数据流分叉**（概念层） | 「模块坐在数据流哪里」唯一权威 |
| **`METHODS.md`** | **论文方法章**（科学表述、阈值、参考文献） | 「论文里怎么写」唯一权威 |
| **`RUNNING_PROCESS.md`** | **端到端叙事**（跑起来发生什么） | 叙事版，**表格一律指向 `RUN_GUIDE.md`** |

**冲突判定规则**：
1. 机器可核事实（stage 列表、依赖、参数名、路径、config 键、checkpoint）→ **代码为最终真相**；
   `RUN_GUIDE.md` / `STAGE_REFERENCE.md` 必须与代码一致，由
   `_consistency_check_20260916/verify_docs_vs_code.py` 逐项对拍（不符即非零退出）。
2. 科学表述（阈值、判据、参考文献）→ **`METHODS.md`**。
3. 概念归属（谁吃谁的产物）→ **`PIPELINE_FLOW.md`**。
4. 叙事文档（`RUNNING_PROCESS.md`）**不得**再复制表格；一旦复制就会二次漂移。

> 报告类文档（`AUDIT_REVIEW_*` / `FIX_REPORT_*` / `TEST_REPORT_*` / `AUDIT_FIX_REPORT_*` /
> `AUDIT_STATUS_*`）是**时点快照**，不随代码更新；只作历史证据，不作运行依据。

---

## 1. 工具链架构

```
virome_phylo_pipeline/
├── phylo_pipeline.py           # 总管线（96 个 CLI 参数；三入口：--extract_dir / --virus / --work_dir）
├── run_phylogeo.py             # BEAST 统一入口（XML 生成 → 提交 → 状态登记）
├── capheine_pipeline.py        # 正选择（stage capheine 的实现）
├── gene_partition_dating.py    # 分基因 BEAST 定年（stage gene_dating 的实现）
├── popgen_analysis.py          # 群体遗传（stage popgen 的实现）
├── recombination_analysis.py   # RDP5 重组检测（stage rdp5 的实现）
├── selection_suite.py          # 正选择下游（独立 CLI，未进总管线）
├── gbk_extractor.py            # GBK → CDS（分基因流上游）
├── config.yaml                 # 默认参数（优先级：CLI 显式 > config.yaml > 代码默认）
├── datasets.yaml               # 数据集注册表（--virus 模式）
├── utils/                      # 依赖模块（51 个）
│   ├── metadata_governance.py  # metadata 治理（日期/地理判定 + has_host 宿主判据）
│   ├── metadata_channels.py    # 三通道回填（dates.csv / sample_metadata.csv / host.csv）
│   ├── beast1_bridge.py        # BEAST 1.x XML 生成（skyline/bdsky/constant）
│   ├── popgen_analysis.py      # pypopart: π/θ/Tajima's D/Fu's Fs/Fst
│   ├── build_haplo_outputs.py  # 单倍型全套输出（发表图/交互 HTML/地图/统计）
│   ├── clock_analysis.py       # 分子钟联合分析（RTT/LTT + DRT）
│   ├── virphy_bridge.py        # Mantel/RRT/VirSpaceTime
│   ├── merge_results.py        # 多链合并（logcombiner + treeannotator）
│   └── ...（共 51 个模块）
└── *.md                        # 文档（权威关系见 §0）
```

> ⚠️ 历史脚本 `merge_results.py` / `clean_workdir.py` / `batch_phylogeo.py` / `pstvd_qc.py`
> **已不在根目录**：`merge_results.py` 迁入 `utils/`（改 `python -m utils.merge_results`），
> `clean_workdir.py` 与 `pstvd_qc.py` 已移除，`batch_phylogeo.py` 的功能由
> `--stage phylogeo` 批量模式取代。旧文档里的根目录调用全部作废。

---

## 2. BEAST 1.x 踩坑记录（务必先读）

### 2.1 parser 未注册 → "Object with idref=bdsky has not been parsed"
**现象**: `<birthDeathSkyline>` 元素被当未知跳过, idref 找不到。
**根因**: BEAST 1.10.4 的 `release_parsers.properties` 没登记
`BirthDeathSerialSkylineModelParser` (仅 BEAUti 可用, 运行时未注册)。
**解法**: classpath 前缀放自定义 `release_parsers.properties` (原 507 类 + SerialSkyline):
```
java -cp /home/zhangwenda/bdsky-override:beast.jar:beast-beagle.jar \
     dr.app.beast.BeastMain -overwrite file.xml
```

### 2.2 模型源码 bug → "Incorrect dimension of bounds"
`BirthDeathSerialSkylineModel` 构造器里
`p.addBounds(DefaultBounds(∞,0, origin.getSize()))` 误把 origin 的 bounds 加到
sampleProbability 上 → 维度冲突。**解法**: sampleProbability 和 origin 用标量
(dimension=1), times/birthRate/deathRate/psi 用 n 维。

### 2.3 ctmcScalePrior 溢出 → "Total = Infinity"
小树高 + 短序列 (PSTVd 359bp) 时 ctmcScalePrior 密度溢出。**解法**: bdsky 模式
不加 ctmcScalePrior (ucld.mean 有 exponential prior 兜底)。

### 2.4 operator 限制
- `scaleOperator` 只接受无界参数; 有界 (lower/upper) 用 `randomWalkOperator`
- 2 维参数不能加 XML lower 属性 (与模型 DefaultBounds 维度冲突)
- `-overwrite` 是 BEAST 应用参数, 必须在类名后

### 2.5 logcombiner trees 极慢 (1.2GB 1小时未果)
**解法**: trees 文件是纯文本行, 用 python 按行跳过 burnin 拼接 (2 秒完成,
49510 棵树), 然后 sumtrees 出 MCC。

### 2.6 treeannotator 单线程巨慢 (49510 棵树 2 小时未果)
**解法**: 抽样 (每 10 棵取 1 → 4951 棵) + DendroPy sumtrees (-s mcct), 2 分钟完成。

---

## 3. 群体遗传学流程 (pypopart, 无 R 依赖)

```
① 数据质控: coverage≥85% 过滤低质量序列 (FastViromeExplorer 低深度组装
            会产生 40-85% 覆盖的嵌合序列, 差异虚高)
② 单倍型:   pypopart identify_haplotypes
③ 多样性:   pypopart stats.popgen (π/θ/S/Tajima's D/Fu's Fs)
④ 网络:     pypopart MJN (无向) + VirNA MSN (有向, HC 突变累积约束)
⑤ Fst:      python 位点法 + permutation
⑥ 输出:     create_publication_figure (SVG/PDF) + 交互 HTML + 地图
```

**注意事项**:
- pypopart 单倍型的 population 需构建后手动 `hap.add_sample(sid, prov)`,
  sample_id 格式 `Location/ID/date` 取 `split('/')[1]`
- pypopart 的 Tajima's D 用"过滤 gap 列"口径 (与 pegas pairwise-deletion 有差异,
  有 gap 数据注意标注)
- **全唯一单倍型** (Hd=1.0) 是保守病毒的特征不是缺陷: 星状网络如实呈现,
  不强行聚类; 单倍型网络视觉分区 ≠ Fst 显著 (频率差异才显著)

---

## 4. 定年决策标准 (方法学亮点)

```
定年判据: DRT 经验 p 值 (主判据, 百分位 ≥ 95 即通过, 见 utils/temporal_signal.py)
          root-to-tip R² 只作参考 (经验参考线 R² > 0.3, 非硬门槛, 不作定年开关)
GCVA (16kb): R²=0.023 显著 → TMRCA 33.9 年 ✅
PSTVd (359bp): R²=0.019 不显著 → 不报 TMRCA ✅ (文献标准)
类病毒 359bp + 位点独立突变 + 弱时间信号 = 分子钟不可靠
```

> 完整方法学表述（阈值分档、阴性结论）见 `METHODS.md`。

### 4.1 DRT 方法学诊断脚本（`phylo_results/`，只读，手动运行）

DRT 主判据由 `utils/temporal_signal.py` 产出（经 `phylo_pipeline.py` 编排）；下列诊断脚本在
`phylo_results/` 下**手动运行**（与本地既有 100+ 诊断脚本同一模式），用于支撑结论的统计严谨性：

| 脚本 | 回答的问题 | 输出 |
|---|---|---|
| `_drt_input_diag.py` | 采样日期构成多少簇、各簇地点/宿主构成；共享日期的序列在 ML 树上是否单系（无根 Steiner 树口径） | `DRT_INPUT_STRUCTURE.md` |
| `_drt_cluster_diag.py` | 簇级置换的排列空间上限（Murray 2015：p 下限 = 1/唯一排列数）；R² 判据与 pos_root 统计量的实际区分力 | `DRT_CLUSTER_DIAGNOSIS.md` |
| `_drt_mantel.py` | **时间×地理混淆**：Mantel 与部分 Mantel（控制地点后时间还剩多少解释力，反之亦然），1000 次置换 | `DRT_MANTEL.md` |

> **注 (2026-09-28)**: 本节所列诊断脚本已整体归档至 `archive/20260928_phylo_results_diag/`
> (含 _drt_* 全系、RunGD.java、_drt_inputs/ 输入数据；结论已沉淀于 DRT_*.md 文档)。
> 另有一批口径核查脚本（`_drt_tt3.py` 主引擎、`_drt_tips.py` 输入口径、`_drt_verify.py` 数据错配核查、
> `_drt_compare_all.py` 口径汇总等）同样位于 `phylo_results/`。
> 审稿人关注的「时间与地理共线」问题由 `_drt_mantel.py` 给出带 p 值的结论（文献：Murray 2015 / Mao 2019 / Guan 2018 / Long 2024）。

---

## 5. 常用命令

运行命令的**完整参考见 `RUN_GUIDE.md`**（本处只留高频速查）：

```bash
# 总管线（先 --dry 预览执行清单）
python phylo_pipeline.py --stage all --dry
python phylo_pipeline.py --extract_dir 06_extraction/ --metadata metadata.tsv --stage all -t 20
python phylo_pipeline.py --virus GCVA --stage phylogeny,time

# BEAST 单病毒入口
python run_phylogeo.py --virus GCVA --prior skyline --chains 5 --threads 8
python -m utils.merge_results --virus GCVA          # 多链合并（原根目录 merge_results.py 已迁入 utils/）

# 群体遗传学
python -m utils.popgen_analysis --fasta X.fasta --metadata Y.csv --group location

# DnaSP 全套 (Fu&Li/Rm/Ka-Ks/LD/SFS...)
python -m utils.dnasp_bridge --fasta X.fasta --out dir

# 按基因 Ka/Ks (需基因注释 CSV + 参考序列名)
python -m utils.dnasp_bridge --fasta X.fasta --out dir --genes genes.csv \
       --kaks-by-gene --reference REF_ACCESSION

# 重组断点分析 (Rm + 断点分布 + 基因映射)
python -m utils.dnasp_bridge --fasta X.fasta --out dir --genes genes.csv --recomb-profile

# 重组分析 (RDP5 完整流程)
# 1) RDP5 检测 (服务器, wine 封装)
bash ~/MMPV-RNA/biosoft/rdp5/run_rdp5.sh alignment.fas my_prefix ~/rdp5_results
# 2) 验证 (事件解析 + IQ-TREE 拓扑切换 + SimPlot + 树图)
python -m utils.rdp5_validate -r ~/rdp5_results/my_prefix.csv -a alignment.fas
# 3) 断点分布图 (需 R + RDP5_RBDP_Rgrapher)
python -m utils.rdp5_bdp_plot -b "my_prefix breakpoint distribution.csv" -g orf_coords.csv -p ...BreakpointPositions.csv -n Virus
# 4) 统一入口 (串联 + dnasp Rm 交叉验证)
python -m utils.recombination_analysis --fasta aln.fas --prefix my --rdp5-csv my.csv \
       --genes genes.csv --outdir out

# 重组区处理 (mask_recombination.py, 默认排除跨末端事件)
# 策略 A: mask 重组区间 → N (区域层面稳健性验证)
python -m utils.mask_recombination --fasta aln.fasta --rdp5 my.csv --genes genes.csv --outdir mask_out
# 策略 B: 删除明确重组序列 (保守, 建树推荐, 默认仅删高置信≥3方法重组子)
python -m utils.mask_recombination --fasta aln.fasta --rdp5 my.csv --genes genes.csv --outdir drop_out --drop-recombinants

# 单倍型全套输出 (着色: location/host)
python -m utils.build_haplo_outputs --virus GCVA --color-by host
```

> ⚠️ `python pstvd_qc.py`（旧的 coverage≥85% 质控构建脚本）**已不存在**；
> 该质控现由 `--clean_*` 系列参数（stage `clean`）承担，见 `RUN_GUIDE.md` §四。

---

## 6. 服务器路径

```
代码:     /home/zhangwenda/MMPV-RNA/virome_phylo_pipeline/
GCVA:     .../8.Virus_full.result/Cytorhabdovirus_sp_lycii_OR489165.1/gcva_analysis/
PSTVd:    .../8.Virus_full.result/Potato_spindle_tuber_viroid_NC_002030.1/analysis_v2/
BDSKY:    /home/zhangwenda/bdsky-override/ (classpath override)
工具:     VirNA (/tmp/VirNA), pypopart (/tmp/pypopart, pip install -e)
```

---

## 7. 关键结果存档

| 结果 | 位置 |
|------|------|
| GCVA 主图 (MCC/skyline/迁移) | gcva_analysis/figures_main/ |
| 单倍型网络图 (发表级) | *analysis/haplo_outputs/ |
| 群体遗传统计 | popgen_analysis CLI 输出 |
| 基因变异映射 | gcva_gene_variation.json |
| 共检出宿主分层 | coabundance_host_strata_results.tsv |
