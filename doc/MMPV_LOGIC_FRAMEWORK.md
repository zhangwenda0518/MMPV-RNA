# MMPV 逻辑框架与核心完整性说明

> 2026-09-05 整理。回应两个问题：① MMPV 的模块逻辑框架；② 昨日合并/清理后核心内容是否受损。
> 结论先行：**核心目录零文件缺失（git index 与磁盘逐一比对通过）；管线间 270 处脚本引用全部核查，真实断裂仅 3 处，且均为本会话之前已知的迁移残留，不影响任何活跃调用链。**

---

## 1. 平台定位

MMPV (Massive Meta-mining of Plant Viruses)：从**公共测序档案**（NCBI SRA / CNCB GSA）出发，经**病毒发现 → 已知病毒深度分析 → 系统发育/进化 → GenBank 提交**的完整闭环平台。Python 编排 + 外部 bioconda 工具 + 少量 R，全部阶段带检查点（`.ok` 标记 + 产物核验）与断点续跑。

## 2. 模块与数据流（核心 8 目录 + 5 根配置）

```
                     ┌────────────────────────────────────┐
                     │  ① public_metadata_pipeline/       │  数据入口
                     │  public_data_pipeline.py           │
                     │   search → info → down → convert → plot
                     │  build_host_pipeline.py            │
                     │   genome-down → hostdb (kraken2/bowtie2/hisat2/minimap2)
                     └──────────┬─────────────────────────┘
              Global_Unified_Metadata_Core14.tsv │  去宿主 reads (00b_HostDepletion)
                                │                │
          ┌─────────────────────┼────────────────┼──────────────────────┐
          ▼                     ▼                │                      │
┌──────────────────────┐  ┌──────────────────────────────────────────┐   │
│ ③ metadata_gui/      │  │  ② virome_discovery_pipeline/  病毒发现 │   │
│  (PySide6 桌面)      │  │  virome_pipeline.py  15 阶段编排:         │   │
│  search_bridge 动态  │  │   clean→deplete→assembly→          │   │
│  加载 gsa_sra.*.py   │  │   identification(10工具)→filter(UniRef90+ │   │
│  viz_panel 调 plot   │  │   CDD 两层)→cobra延伸→merge→cluster       │   │
│  Core14 表浏览/编辑  │  │   →taxonomy(8工具R投票)→host→checkv→      │   │
└──────────────────────┘  │   rescue(四分支A/B/C/D)→analysis→         │   │
                          │   analysis_verify(5层证据)→report          │   │
                          │  产出: 08_Rescue/HQ_plant_viruses.fasta    │   │
                          │        final_judgement_table.tsv           │   │
                          └───────┬───────────────────────┬───────────┘   │
                                  │ final_centroids.fasta │               │
                                  ▼ (known参考)           ▼ (新病毒)      │
              ┌──────────────────────────────┐  ┌───────────────────────┐  │
              │ ④ virome_analysis_pipeline/  │  │ ⑥ virome_submission_  │  │
              │  auto_known_virus.py 9阶段:  │  │    pipeline/          │  │
              │   detect(Salmon)→filter→     │  │  submission_pipeline  │  │
              │   variants(ivar/fb/lofreq+   │  │   .py known|discovery │  │
              │   SnpEff+SNPGenie)→post(6套件│  │   taxonomy→features→  │  │
              │   VCF/PCA/LD/MAF)→full(12步  │  │   hypothetical→       │  │
              │   组装)→extract→similarity   │  │   suvtk table2asn→    │  │
              │   (SDT)→dvg(ViReMa)→report   │  │   .sqn + sequin 模式  │  │
              │  产出: 01..09 目录 + HTML报告 │  └──────────┬────────────┘  │
              └──────┬───────────────────────┘             │               │
                     │ 全长基因组                            │ .sqn 提交件   │
                     ▼                                       ▼               │
              ┌──────────────────────────────────────────────────┐         │
              │  ⑤ virome_phylo_pipeline/  系统发育与分子进化      │         │
              │  phylo_pipeline.py 阶段化:                        │         │
              │   MAFFT→IQ-TREE(MFP)→时间信号门禁(RTT/DRT/BETS)→  │         │
              │   BEAST 定年(仅过门禁者)→CTMC系统地理→popgen      │         │
              │   (pypopart)→选择(FEL/MEME/BUSTED/codeml)→       │         │
              │   重组(RDP5)→发表图                               │         │
              └──────────────────────────────────────────────────┘         │
                                                                           │
              ┌──────────────────────────────────────────────────┐         │
              │  ⑦ submission_gui/  (PySide6)                    │◄────────┘
              │  unified_metadata.csv 多标签编辑/占位符高亮/      │
              │  template.sbt 生成/验证报告 → 与⑥经 sync 闭环     │
              └──────────────────────────────────────────────────┘

  支撑层:
    biosoft/    第三方免编译工具 (snpEff.jar, ViReMa, VirBot, VirHunter权重, ...)
    database/   ICTV MSL41 分类表 + VMR + 参考库种子 (大型DB在服务器 ~/database)
    doc/        17 篇阶段文档 + METHODS_TEMPLATE(五管道方法章) + 依赖清单
    pixi.toml / pipeline_config.yaml / DATABASE_SETUP.md / SOFTWARE_VERSIONS.txt
```

## 3. 关键交接文件（管线间契约）

| 交接物 | 生产者 → 消费者 | 说明 |
|---|---|---|
| `Global_Unified_Metadata_{Core14,Full}.tsv` | ① → ③④⑥⑦ | 14 列核心元数据，GUI 编辑与提交元数据的唯一来源 |
| `00b_HostDepletion/*.fa.gz` | ①preprocess/②deplete → ②③ | 去宿主 reads，两管线共用 |
| `04_CLUSTER/4_centroids/final_centroids.fasta` | ② → ②host/④ | HQ vOTU 代表序列 |
| `08_Rescue/HQ_plant_viruses.fasta` + `final_judgement_table.tsv` | ② → ⑥ | 新病毒提交输入 + KEEP/REVIEW/DROP 裁决表 |
| `01_detection/.../all_viruses.best.summary.tsv` | ④S1 → ④S2-S9 | STAGE1_SCHEMA 校验的样本×病毒主表 |
| `02_filtering/high_conf.summary.tsv` | ④S2 → ④S3+ | 高置信过滤后主表（S2 之后 get_summary 自动切换） |
| `05_assembly→06_extraction` 全长序列 | ④ → ⑤⑥ | 系统发育与提交的序列来源 |
| `submission/` (source.src, miuvig.tsv, assembly.tsv, template.sbt) | ⑥ ↔ ⑦ | GUI 编辑后经 sync 回写服务器重出 .sqn |

## 4. 跨模块桥接（代码级）

- `metadata_gui/controllers/search_bridge.py`：importlib 动态加载 ① 的 `gsa_sra.search.py`（同进程 QThread）；`deep_extract` 走 subprocess 调 `gsa_sra.info.py`；`viz_panel` 动态加载 `gsa_sra.plot.py`
- ② `virome_pipeline.py` 的 analysis 阶段直接 subprocess 调 ④ 的 `batch_virus_depth.py`；analysis_verify 调 `~/bin/acvirus_tree_pro.py`（回退）与 ④ 的 `sdt_genus_matrix.py`
- ⑥ `submission_pipeline.py` 步骤 2.6 调 ① 的 `gsa_sra.info.py -m local` 回填元数据；`tbl2gb` 借用 ④/utils/tbl2gb.py
- ⑥⑦ 经 `sync_client.py`（SSH）+ 服务器端 `sync_sqn_from_csv.py` 闭环
- ⑤ `datasets.yaml` 以绝对路径挂载数据集（PSTVD_FULL、GCVA 等），结果统一落 `results_PSTVd_GCVA_20260905/`（原路径符号链接兼容）

## 5. 核心完整性审计（2026-09-05，回应对归档的关切）

**方法**：git index ↔ 磁盘逐文件比对 + 全管线脚本引用扫描（270 处 `*.py/*.R/*.sh` 引用逐一验证目标存在）。

**昨日改动中涉及核心目录的删除，全部有冗余替代**：
| 被删/移动 | 替代 | 状态 |
|---|---|---|
| `virome_discovery_pipeline/utils/analyze_viroid.py`（损坏：引号被剥） | `analyze_viroid_batch.py`（编排器实际调用者） | ✅ 无损 |
| `virome_analysis_pipeline/utils/batch_plot_virus_depth.py`（坏副本：close后存PNG） | 根目录 `batch_plot_virus_depth.py`（编排器调用者） | ✅ 无损 |
| `virome_analysis_pipeline/utils/snpeff_build.py`（旧副本缺关键参数） | `batch_virus_variants.py` 内联构建（含 -noCheckCds） | ✅ 无损 |
| `run_all_species.sh` | 功能已在管线参数化（--species/--taxid） | ✅ 无损 |

**真实引用断裂仅 3 处，均为本会话之前已存在**（git 记录佐证）：
1. `utils/auto_known_virus.py`（旧版 10-stage 编排器，无任何调用方）引用的 `utils/gbk_extractor.py`、`utils/visual_codon_miner.py` —— 这两文件在**更早的 capheine 迁移**中移至 `virome_phylo_pipeline/`（服务器端 R100 记录一致）。旧版编排器保留原地未归档，仅作历史参考。
2. `submission_gui/sync_client.py` 引用 `sync_sqn_from_csv.py` —— 该文件仅存在于服务器端（本来就要在服务器跑）。
3. `build_virus_db.py` 的 3 个辅助脚本（SearchAccessionIdToTaxId.py 等）—— 仅在服务器，该工具为独立建库工具，不在主编排链。

**检查器报告的其余 ~95 处"断裂"经逐一核实均为误报**：跨管线动态定位（SCRIPT_DIR/'utils'/'x.py' 拼写）、运行时自生成脚本（gene_partition_dating.py 生成 parse_gene_beast.py）、外部命令（rnaviralspades/ragtag）、报告 HTML 的文字提及（capheine）、PyInstaller dist 产物等。

**核心目录现状**：`doc database metadata_gui public_metadata_pipeline submission_gui virome_analysis_pipeline virome_discovery_pipeline virome_phylo_pipeline virome_submission_pipeline + pixi.toml pipeline_config.yaml DATABASE_SETUP.md SOFTWARE_VERSIONS.txt README.md` —— git index 内每一个文件磁盘都在，无一缺失。
