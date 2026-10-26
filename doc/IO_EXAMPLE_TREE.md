# MMPV 全管线完整输入/输出目录示例结构（standard 布局）

> 2026-10-07 生成：对 ①→⑧ 全部管线脚本做**代码字面量级重扫**后的完整示例树。
> 与 `doc/IO_LAYOUT_DESIGN.md`（规范文档）的关系：该文档定义布局契约与映射表；
> 本文是**逐脚本核实过的实际产物形态**（含每一层的文件名、每样本占位符、断点标记），
> 两者冲突时以本文记录的代码现状为准，并在 §13 列出契约与实现的偏差清单。
>
> 占位符：`{sample}` = 样本名；`{NAME}` = ⑥ 的基因组名；`{virus}` / `{vkey}` = 病毒名；
> `{S}` = ③ 的 contig 前缀样本名；`{acc}` = SRA/GSA accession。
> 目录键由 `mmpv_common/io_layout.py` 注册表解析；legacy 旧名对照见设计文档 §4。
>
> **【2026-10-11 归档说明】** ⑥ `endogenous_virus_pipeline`（EVE）与 ⑧ `host_expression_pipeline`（宿主表达）
> 已因"方法不成熟"归档至 `archive/20261011_endogenous_virus_pipeline/` 与
> `archive/20261011_host_expression_pipeline/`。下文 **§8（`06_EVE/`）与 §10（`08_HostExpression/`）
> 的目录树为归档前的产物形态存档，命令不可再运行**；历史项目已有的对应输出仍可解读，恢复方法见各归档目录内说明。

---

## 1. 项目根总览（一屏）

```
<project_root>/
├── 01_PublicData/            ← ① 公共数据（public_data_pipeline.py --work-dir）
├── 02_Preprocessing/         ← ② 数据清洗（data_preprocessing.py --output_dir）
├── 03_Discovery/             ← ③ 病毒发现（virome_pipeline.py --output_dir；⑦ 默认 --work-dir 也指这里）
├── 04_Analysis/              ← ④ 已知病毒分析（auto_known_virus.py --output_dir）
├── 05_Phylo/                 ← ⑤ 系统发育（phylo_pipeline.py --output_dir）
├── 06_EVE/                   ← ⑥ 内源病毒（eve_screen.py -o）【⑥ 已归档 2026-10-11，不可运行】
├── 07_*                      ← （空缺）⑦ 提交管线无独立输出根，产物落 <项目>/submission_<run_title>/
├── 08_HostExpression/        ← ⑧ 宿主表达可选模块（host_expression.py -o/--project-root）【⑧ 已归档 2026-10-11，不可运行】
├── 90_Handoff/               ← 跨管线交接层（③→④ 桥 + 边界台账）
└── submission_<run_title>/   ← ⑦ 提交产物（work_dir 的父目录下）
```

---

## 2. 输入侧总清单（开跑前要准备的东西）

| 类别 | 内容 | 形态 |
|---|---|---|
| 必填 CLI | ① 的 `--species`（拉丁名）+ `--taxid`（int）；各编排器 `--work-dir/--output_dir` 指向 §1 各自的根 | — |
| 凭据 | `DEEPSEEK_API_KEY` / `NCBI_API_KEY` 环境变量或 `--deepseek-api` / `--ncbi-api`（缺失时 ① info 降级 local 模式） | env 或 CLI |
| 检索控制（可选） | `--source-type`（默认 TRANSCRIPTOMIC）、`--prescreen-filter {off,safe,strict}`、`--taxdump`（含 nodes.dmp 的目录）、`--skip-list`、`--cra-map` | 文件/目录 |
| 宿主参考（可选） | `--host-fasta`（已有基因组 FASTA 则跳过下载）；否则 ① 自动走 datasets/ngd/FTP/gget 四通道 | FASTA |
| 外部数据库 | 宿主库（kraken2 hash.k2d + bowtie2/hisat2/minimap2 索引，① hostdb 阶段或 ② `--host_db` 自动建）；③ 鉴定库 `--virus_db`、`--uniprot_db`、`--nr_db`、CheckV 库、8 套分类库（geNomad/CAT/MMseqs/VITAP/ACVirus/VContact3…）；⑥ 的 ref/pv/RVDB/id2div/viroids DIAMOND 库；⑦ 的 suvtk 库 + RVDB v31 目录；④ 的 SnpEff 库（可自动建） | 目录 |
| 外部工具 | fasterq-dump/prefetch/pigz、esearch/efetch(entrez-direct)、datavzrd、datasets CLI、kraken2-build、megahit/rnaviralSPAdes/penguin、10 鉴定工具、COBRA、Flye、vclust/cd-hit、IQ-TREE、RDP5、BEAST1、salmon、suvtk、table2asn 等（各管线 README 有版本表） | PATH |
| 跨管线交接 | 全部由上游产物自动衔接，见 §12 边界契约；无需人工搬运 | — |

---

## 3. ① 公共数据 — `01_PublicData/`

**编排**：`public_data_pipeline.py --species <拉丁名> --taxid <int> --stage all --work-dir <项目>/01_PublicData`
**stage 顺序**：`search → info → prescreen → down → convert → plot → hostref → hostdb → report`（`.checkpoints/state.json` 记状态，`logs/{stage}.log` 逐阶段日志）

```
01_PublicData/
├── .checkpoints/state.json                 # utils/pipeline_utils.Checkpoint：stage 状态 JSON
├── logs/{search,info,prescreen,down,convert,plot,hostref,hostdb,report}.log
├── 01_Search/                              # gsa_sra.search.py
│   ├── SRA_GSA_Merged_Final.csv            # ★检索合并总表（按 Run 去重）
│   │                                       #   列: Database,Run,BioProject,BioSample,ScientificName
│   │                                       #   [,Tissue,Age_GrowthStage,Location(detailed)],
│   │                                       #   LibraryStrategy,LibrarySource,Platform,CenterName,ReleaseDate,FileSize_GB
│   ├── download_links.csv                  # 列 Run,Database,url_1,md5_1,url_2,md5_2（down 阶段零步输入）
│   ├── AI_Summary.txt                      # 仅 --summary + DeepSeek key
│   └── GSA_Results/
│       ├── 0_web_cache/{acc}.html          # GSA 页面缓存（acc=CRR/CRX…）
│       └── 1_xls_cache/{cra}.xlsx          # CRA 项目表缓存（仅 detailed）
├── 02_Unified/                             # gsa_sra.info.py；★全平台元数据权威
│   ├── sra.list                            # 编排器从 Merged csv Run 列提取（prescreen/down 输入）
│   ├── Global_Unified_Metadata_Core14.csv/.tsv   # ★Core14（⑤ --metadata 直接消费）
│   │                                       #   14 列: Run,ReleaseDate,CollectionDate,Location,Source,Tissue,
│   │                                       #   Age_GrowthStage,ScientificName,TaxID,LibrarySource,CenterName,
│   │                                       #   BioProject,BioSample,PMID
│   ├── Global_Unified_Metadata_Full.csv/.tsv     # 全列版（⑦ 元数据自动填充用）
│   ├── Config_…yaml + Report_…/index.html  # datavzrd 交互页
│   ├── SRA_Results/
│   │   ├── 1_raw_xml/{srr}.xml → 2_full_json/{srr}_full.json
│   │   ├── 3_1_local_parsed/{srr}_local.json → 3_2_api_inferred/{srr}_api.json → 3_3_ai_arbitrated/{srr}_arb.json
│   │   └── SRA_Ultimate_Merged.csv/.tsv
│   ├── GSA_Results/
│   │   ├── 0_web_cache/{acc}.json → 1_csv/{acc}.csv → 2_xlsx/{cra}.xlsx + {acc}_merged.csv
│   │   ├── 3_api_enhanced/{q_id}_api.json
│   │   └── GSA_Ultimate_Merged.csv/.tsv    # 固定 28 列
│   └── BioProject_Results/{1_cache/{bp}_trace.json, BioProject_Trace_Details.csv}  # 缺 PMID 项目溯源
├── 02b_Prescreen/                          # stat_prescreen.py（不下载数据的预筛）
│   ├── sample_decision.tsv                 # ★逐库决策: acc,source,cra,status,…,decision,decision_reason,pass…
│   │                                       #   decision ∈ DOWNLOAD/REVIEW/SKIP/LOCAL_CHECK
│   ├── stat_run_summary.tsv                # 兼容宽表
│   ├── stat_viral_taxa.tsv                 # 长表: acc,source,tax_id,org,path,self_count,total_count,pct
│   ├── prescreen_run.list / _local_check.list / _skip.list   # 按 decision 分流的 acc 清单
│   ├── prescreen_filtered_{safe|strict}.list  # --prescreen-filter 生效时（down 按此下载）
│   ├── prescreen_meta.json / STAT_PRESCREEN_REPORT.md
│   ├── stat_raw/{acc}.json  gsa_raw/{crr}.json   # 抓取缓存（空文件=NO_STAT 事实缓存）
│   └── virus_profile/                      # stat_virus_profile.py（有 stat_viral_taxa.tsv 时）
│       ├── virus_prevalence.tsv  virus_sample_matrix.tsv  virus_top_by_sample.tsv
│       ├── virus_profile_meta.json  VIRUS_PROFILE_REPORT.md
│       └── figures/01_virus_prevalence_bar.png … Combined_Virus_Profile.png（仅 --plot）
├── 03_RawData/                             # gsa_sra.down.py + sra2fastx.py
│   ├── {acc}.sra                           # 下载产物（平铺；convert 后删除）
│   ├── download_report_{YYYYmmdd_HHMMSS}.csv     # Accession,Sample,Method,Status
│   ├── failed_sra_{YYYYmmdd_HHMMSS}.txt    # 失败清单（可作 --list 重跑）
│   └── {acc}_1.fastq.gz / {acc}_2.fastq.gz / {acc}.fastq.gz   # ★convert 后（PE/SE）——② 的输入
├── 04_Plots/
│   ├── from_search/Combined_Landscape_Full.pdf + .png    # 3x2 六面板（PNG 600dpi）
│   └── from_info/Combined_Landscape_Full.pdf + .png
├── 05_HostRef/
│   ├── genome/
│   │   ├── all.genome.uniq.fasta           # ★宿主基因组去重合并（⑥/⑧ 的输入；⑥⑧ 均已归档 2026-10-11）
│   │   ├── all.genome.uniq.fasta.meta.txt  genome_report.txt
│   │   ├── extracted/  genome_down.zip  gget_ref.json    # 各下载通道中间物
│   │   └── {chloroplast,mitochondrion}/*.fasta           # 仅 --include-organelles
│   └── hostdb/
│       ├── kraken2/{hash,opts,taxo}.k2d  host_with_taxid.fasta
│       ├── bowtie2/host.{1..}.bt2[l]
│       ├── hisat2/host.{1..}.ht2
│       ├── minimap2/host_{dna-short|rna-short|nanopore|pacbio}.mmi
│       └── .checkpoints                    # build_hostbase 分步断点（纯文本 step_id）
├── Pipeline_Summary_Report.html            # generate_report.py（report 阶段）
└── sample_handoff.csv                      # 列 Run,Layout(PE/SE),R1_Path,R2_Path,SE_Path,Size_MB,
                                            #   ScientificName,Tissue,Location,BioProject,Status(ready|missing)
```

---

## 4. ② 数据清洗 — `02_Preprocessing/`

**编排**：`data_preprocessing.py --input_reads <①>/03_RawData --host_db <①>/05_HostRef/hostdb --output_dir <项目>/02_Preprocessing --stage all`
**stage 顺序**：`clean → deplete → [bbnorm 仅显式] → report`；子脚本 `.clean_checkpoints` / `.checkpoints` 逐样本断点

```
02_Preprocessing/
├── preprocessing.log                       # 编排器 DEBUG 日志
├── preprocessing_summary.tsv               # preprocess_report.py；列: Sample Raw_Reads Fastp_Clean_Reads
│                                           #   Fastp_Q30(%) Depletion_Input After_Kraken2 After_HostAlign
│                                           #   After_rRNA Retained(%) Aligner_Rate(%) Final_Files Status
│                                           #   （末行 TOTAL；Status ∈ PASS/LOW_RETAINED/NO_OUTPUT/NO_DATA）
├── assembly_ready.list                     # ★最终 reads 绝对路径清单（人工核对/下游取数）
├── preprocessing_report.html               # Chart.js 内嵌离线报告
├── host_reference/                         # 仅现场构建宿主库时（host_fasta / host_species+host_taxid）
│   ├── genome/all.genome.uniq.fasta
│   └── hostdb/{kraken2/, bowtie2|hisat2|minimap2/}
├── 01_CleanData/                           # clean-data.py（fastp→seqkit→clumpify）
│   ├── .clean_checkpoints                  # 已完成样本名逐行（--force 重建）
│   ├── logs/{sample}_fastp_report.html|.json     # fastp QC（JSON 是 report 阶段统计源）
│   │      {sample}.01_fastp.0_fastp.log  {sample}.02_seqkit.{0,1}_seqkit.log
│   │      {sample}.03_clumpify.0_clumpify.log
│   ├── 1.fastp_tmp/{sample}_1.fq.gz        # fastp 出；seqkit 成功后即删（目录常残留为空）
│   ├── 2.fasta/{sample}_1.fa.gz            # seqkit fq2fa -w 70；clumpify 成功后删
│   └── 3.clumpify/{sample}_1.fa.gz         # ★clumpify reorder（②→⑧ 取 clean reads 的首选子目录）
│       # 注意: --stage all 跑完后 deplete 成功的样本数据文件会被清理，
│       # 只剩 logs/ 与 fastp JSON —— 这是编排器 cleanup_00a_after_deplete 的预期行为
├── 02_HostDepleted/                        # host_depletion.py（kraken2→比对去宿主→去 rRNA）
│   ├── .checkpoints                        # 已完成样本名逐行
│   ├── {sample}_clean_1.fa.gz / _2.fa.gz   # ★★PE 去宿主去 rRNA 最终 reads —— ③ --input_reads 的正式落点
│   ├── {sample}_clean.fa.gz                # SE 形态
│   ├── rrna/                               # 仅 --keep_rrna
│   │   ├── {sample}.rrna.{1,2}.fq.gz       # ribodetector 分支（SE {sample}.rrna.fq.gz）
│   │   └── {sample}.silva_rrna.{1,2}.fq.gz # silva 分支
│   └── logs/
│       ├── host_depletion_seqkit_summary.tsv   # Sample/Stage/…；Stage=1_Raw/2_Kraken2_Filtered/
│       │                                       #   3_Host_Filtered/4_rRNA_Filtered
│       ├── host_depletion_resource_usage.tsv   # Sample/Elapsed_Seconds/Peak_Mem_MB/CPU_Time_Seconds/Avg_CPUs
│       ├── ribodetector.report.txt             # 仅 --rrna（Sample/Total_sequences/non_rRNA/rRNA/…%）
│       ├── host_depletion_plot_line.png / _bar.png
│       ├── kraken2_report/{sample}_kraken2_report.txt   # 唯一保留的 kraken2 产物
│       └── logs/{sample}.step1_kraken2.0_kraken2.log    # 二层 logs/：逐样本逐步骤工具日志
│                 {sample}.step2_align_filter.{0_bowtie2,1_samtools}.log
│                 {sample}.step{3_sort,4_extract,5_ribodetector,5_silva_*,6_*,7_compress_*}.*.log
│                 {sample}.resource.log
└── 03_BBnorm/                              # 仅 --stage bbnorm（run_bbnorm.py）
    ├── {sample}_clean_norm_R1.fq.gz / _R2.fq.gz   # PE（恒 .fq.gz）
    └── {sample}_clean_norm_SE.fq.gz        # SE/无标签
```

---

## 5. ③ 病毒发现 — `03_Discovery/`

**编排**：`virome_pipeline.py --input_reads <②>/02_HostDepleted --output_dir <项目>/03_Discovery --stage all`
**stage_order**：`clean → deplete → assembly → identification → filter → graph → cobra → merge → cluster → taxonomy → host → checkv → rescue → analysis → analysis_verify → report`
**断点**：`03_Discovery/.{stage}.ok`（16 个阶段标记）；阶段日志统一 `13_Report/logs/{stage}.log`

```
03_Discovery/
├── .{stage}.ok                             # 阶段断点（clean/deplete/assembly/…/report）
├── final_judgement_table.tsv               # ★根级判定表（gen_final_judgement.py）
│                                           #   11 列: contig_id,host_species,length,prevalence_%,
│                                           #   nt_pident,nt_qcovs,nt_species,aa_pident,aa_species,cdd_top,category
├── prevalence_full_table.tsv               # ★根级流行率表（panvirome_prevalence.py；⑦/9b 也读）
│
├── 01_CleanData/ 02_HostDepleted/          # 仅 ③ 内跑 clean/deplete 时（推荐 ② 独立跑）
│
├── 01_Assembly/                            # assembly_pipeline.py（多工具）
│   ├── logs/ resource_usage_summary.tsv  assembly.log
│   └── {sample}/
│       ├── {sample}_megahit.contig.fasta   # 每工具一套（头已加 {sample}_ 前缀）
│       ├── {sample}_rnaviralspades.contig.fasta (+ .scaffolds.fasta)
│       ├── {sample}_penguin.contig.fasta
│       ├── {sample}.{tool}.log  {sample}.{tool}.time.mem.log
│       ├── {sample}_all_tools_refineC_merge/*.merged.fasta   # 多工具 refineC 合并（最终移回样本根）
│       └── {sample}.all_tools.split_merged.fasta
│   └── coassembly_merged/ALL_merged_R1.fq.gz, _R2.fq.gz      # 仅 --coassembly
│
├── 02_Identification/                      # virus_identification.py（10 工具并集）
│   ├── ident.log  resource_usage.tsv  resource_summary.tsv
│   ├── summary_plots/All_Samples_Comparison.png
│   └── {sample}/
│       ├── {sample}_virus.all.candidate.fasta        # ★★10 工具并集候选（02a 候选，主链锚点）
│       ├── {sample}_virus.all.result.id
│       ├── {sample}_virus.{blast|genomad|rdrpcatch|viralm|virbot|virsorter2|viralverify|virhunter|metabuli|viroid}.result.id
│       ├── {sample}.virus.candidate.fasta / {sample}.viroids.candidate.fasta   # 长度/类病毒预筛
│       ├── {sample}.processing.done        # 样本级断点
│       ├── {sample}_resource.tsv  seqkit_filter.log
│       ├── blast_output/ genomad_output/ virsorter2_output/ viralverify_output/
│       │   virhunter_output/ metabuli_output/ viroid_output/ rdrpcatch_output/
│       │   viralm_output/ virbot_output/             # 每工具原始输出
│       ├── comparison_plots/{sample}_{venn_diagram,upset_plot,bar_charts}.png
│       ├── uniprot_filter_output[_filter|_strict]/    # UniProt strict 过滤
│       └── nr_filter_output[_filter|_strict]/
│
├── 03_Filter/                              # filter_virus.py（UniProt+CDD 分层）+ graph_propagation.py
│   ├── filter_merged.log  all_candidates.fasta   # 合并中间件（拆完即删）
│   ├── candidates/{S}/{S}_virus.all.candidate.fasta  # ★COBRA 主链镜像（filter 写，graph 裁决后覆写）
│   ├── {S}/
│   │   ├── {S}.virus.candidate_filtered.fasta        # 最终过滤产物（uniprot→cdd→02a 符号链）
│   │   ├── {S}.virus.candidate_uniprot_filtered.fasta / _cdd_filtered.fasta / _virus.all.candidate.fasta
│   │   ├── filter.log
│   │   ├── cdd/{…}_cdd_{filter|strict}.fasta  cdd_calls.tsv  cdd.unclassified.id
│   │   └── uniprot/{…}_uniprot.fasta  uniprot.unclassified.id  diamond_uniprot/blastx.tsv
│   └── graph/                              # graph 阶段（配对读段链接图 + Louvain 裁决）
│       ├── graph.log
│       ├── candidates/{N}/{N}_virus.all.candidate.fasta   # ident 候选 − graph 非病毒
│       └── {S}/
│           ├── {S}.graph_verdict.tsv       # contig,source_tool,length,gc,depth,votes,seed_label,
│           │                               #   community_*,verdict,reason
│           ├── {S}.graph_summary.tsv  {S}.graph.done
│           ├── {S}.virus.candidate.final.fasta       # filter 链 − demote ∪ rescue
│           └── {S}.graph_adjust_calls.tsv  # contig,cdd_call,blastx_hit,graph_verdict,action
│
├── 04_COBRA/                               # cobra_pipeline.py --mode mix（单样本延伸）
│   ├── cobra.log  checkpoint_status.json
│   └── {base_sample}/{cobra_{megahit|rnaviralspades|penguin}_result}/
│       ├── {base}.mix.{asm}.virus.fa       # 归一化候选
│       ├── {base}.mix.{asm}.coverage.txt   # CoverM contig\trpkm
│       ├── {base}.mix.{asm}.cobra.fa       # ★延伸结果（merge 阶段 rglob *.cobra.fa 收集）
│       ├── {base}.mix.{asm}.COBRA/COBRA_category_*.fasta   # cobra-meta 原始输出
│       └── {base}.mix.{asm}.CoverM.txt  {base}.mix.{asm}.log
│
├── 05_CoAssembly/                          # merge 阶段
│   ├── all_sample_virus.fasta              # 全样本 COBRA 合并
│   ├── all_sample_virus_combined.fasta     # ★Flye 延伸版拼接 —— cluster 输入
│   └── flye_coassembly/{flye.log, flye_output/assembly.fasta, flye_mapping.tsv,
│                        flye_sample_map.json, flye_summary.txt, trace.log}
│
├── 06_CLUSTER/                             # cluster_pipeline.py
│   ├── cluster.log  global_summary.txt  ref_merged.fasta  rescue_centroids.fasta
│   ├── 1_seqkit/virus.candidate.fasta
│   ├── 2_cdhit/{cdhit_combined.fasta, known_centroids.fasta, known_linked_centroids.fasta,
│   │           known_association.tsv, novel_contigs.fasta, vclust_cdhit/,
│   │           known_clusters/cluster_*.all.fasta}          # 仅 --ref-genomes
│   ├── 3_vclust/{vclust_prefilter.txt, vclust_ani.tsv, vclust_clusters.tsv, vclust_ani.ids.tsv,
│   │            all.cluster.ref.fasta, logs/, split_fastas/cluster_X.{ref,all}.fasta}
│   └── 4_centroids/                        # ⚠ 实写目录名（契约名 centroids/，见 §13-D1）
│       ├── final_centroids.fasta           # ★★★③→④ 交接物（known 直通 + vclust novel 代表）
│       ├── known_ids.txt                   # CD-HIT known 簇代表 ID（免拯救白名单）
│       └── known_association.tsv           # contig_id\treference_accession
│
├── 07_Taxonomy/                            # virus_classifier.py + virus_classifier_analysis.R
│   ├── taxonomy.log  r_consensus.log
│   ├── {Votus}.classed/                    # 8 工具逐级分类（{Votus}=--tax_sample_name，默认 Votus）
│   │   ├── {Votus}_combined_taxonomy.tsv   # 8 工具合并 8 级表（seq_name,tool,Realm..Species）
│   │   ├── {Votus}_{genomad|metabuli|diamond_lca|CAT|mmseqs|VITAP|ACVirus|vcontact3|PhaGCN3}_taxonomy.tsv
│   │   ├── {Votus}_resource_usage.tsv
│   │   └── genomad_annotate/ VITAP_results/ ACVirus_results/ mmseqs_results/ …
│   └── {Votus}.integrated/                 # ★★R 共识（taxonomy 逐样本契约：{sample}.integrated/）
│       ├── final_integrated_classification.tsv   # ★★★整合分类表
│       │   # 列: contig_id,primary_tool,completeness,confidence,Realm,…,Species,{Realm..Species}_agree
│       ├── fill_stats.tsv  analysis_summary.txt  consensus_summary.pdf|png
│       ├── intersection_matrix.tsv  common_ids.txt  upset_data.json  standardized_*.tsv
│       ├── agreement_stats.tsv  consistency_summary.tsv  consistent_summary.tsv
│       ├── comparison_{level}.tsv  consistent_{level}.tsv|/
│       ├── classification_counts.pdf|png  agreement_rates.pdf|png
│       ├── intersection_venn.pdf|png  intersection_upset.pdf|png  combined_tax_level_venn.pdf|png
│       └── taxonomy_gate_stamp.tsv + taxonomy_gate_check.tsv   # 闸门戳（当前树内 R 脚本未写，见 §13-D4）
│
├── 08_HostPrediction/                      # run_host_prediction.py（三方法 ensemble）
│   ├── ensemble_host_summary.tsv           # ★contig_id,Final_Host,Decision_Method,Host…
│   ├── host_classified_fasta/{host}.classified.fasta     # 如 Plant.classified.fasta
│   ├── C9_ICTV_result/{classification_result.tsv, classification_summary.tsv, confidence_report.tsv}
│   ├── RVH_result/{result.csv, RVH_taxa.csv}
│   ├── phabox2_output/final_prediction/cherry_prediction.tsv
│   └── cascade_sub{1_ictv,2_rvh}_miss.fasta
│
├── 09_CheckV/                              # checkv 阶段
│   ├── checkv.log
│   ├── {safe_host}.fasta                   # 按 Final_Host 拆分
│   ├── {safe_host}/completeness.tsv        # checkv 原始输出
│   └── checkv_pass_ids.txt                 # ≥90% completeness 的 contig ID
│
├── 10_Rescue/                              # run_rescue + rescue_pipeline.py（两分支拯救）
│   ├── rescue.log
│   ├── HQ_plant_viruses.fasta              # ★★★HQ 病毒集（免拯救 known + 各宿主 rescue 拼接）
│   ├── known/centroids/final_centroids.fasta   # CD-HIT known + CheckV≥90% 直通
│   └── {Plant}/                            # = _".join(host_filter)，默认 Plant
│       ├── branch_a/{branchA_{pass,fail}.fasta, branch_a_prepass.txt, checkv_out/completeness.tsv}
│       ├── branch_b/{branchB_{pass,fail}.fasta, logs/, scaffolds/, merged_reads/,
│       │            {scaffold.fasta, pilon_out.fasta, run.ok, .vsi_failed, .vsi_error}}
│       │            # (2026-10-11) branch_c/branch_d 已归档 (archive/20261011_rescue_branch_cd/),
│       │            #  历史输出目录中的 branch_c/ branch_d/ 仍可存在于旧运行里
│       ├── merged/all_HQ.fasta             # A+B pass 拼接
│       ├── centroids/final_centroids.fasta # ★HQ vOTU（最终 vclust 去重）
│       ├── rescue_report.tsv               # contig_id,input_len,final_len,branch,method,genus,
│       │                                   #   genus_avg_len,pct_of_genus,checkv,note
│       ├── rescue_summary.md
│       └── checkv/{no_rescue,unknown,all,skipped_*,<host>}/completeness.tsv
│
├── 11_ViromeAnalysis/                      # analysis 阶段（⑦ 也从这里拿 ref_info）
│   ├── HQ_plant_viruses.fasta              # 从 10_Rescue 拷贝
│   ├── All_plant.viruses.fasta             # symlink → 08_HostPrediction/host_classified_fasta/Plant.classified.fasta
│   ├── Viroid.all.fasta                    # symlink → 02_Identification/*viroids.candidate.fasta
│   ├── logs/{nucleic_acid,plant_virus_info,viroid_analysis,rescue_detection}.log
│   ├── HQ_analysis/HQ_plant_viruses_info.tsv
│   ├── all_plant_analysis/{All_plant.viruses_info.tsv, all_plant_viruses_genus_summary.tsv,
│   │                       plant_virus_cluster_summary.tsv}
│   ├── viroid_analysis/{Viroid.all_info.tsv, Viroid.species_info.tsv, Viroid.per_virus.tsv,
│   │                    Viroid.per_sample.tsv}
│   ├── rescue_detection/{rescue_ref_info.tsv, summary/all_viruses.raw.tsv, frequency_table.tsv}
│   ├── ref_info.tsv                        # build_ref_info.py 12 列（⑦ source.src 的 ID/物种来源）
│   └── (可选 suvtk 链: suvtk_taxonomy/ suvtk_features/ hypothetical/ topology.tsv
│        submission/{sequences.fsa, featuretable.tbl, comments.cmt, submission.sqn})
│
├── 12_AnalysisVerify/                      # analysis_verify（9b）阶段
│   ├── virus_validation/                   # CDD + HMM-CT3 + 加权评分
│   │   ├── cdd_hits.tsv  cdd_evidence_report.tsv
│   │   ├── queryDB/ resDB/ tmp/
│   │   ├── keep|review|drop_candidates.{txt,fasta}   # 仅 --mode filter
│   │   ├── hmm_ct3_hits.tsv  hmm_ct3_hits_raw.tsv
│   │   ├── rescue_evidence_scored.tsv      # ★30/30/40 加权评分 + verdict KEEP/REVIEW/DROP
│   │   └── hmm_ct3_rescue_candidates.tsv
│   ├── class_KEEP.fasta                    # ★⑦ 提交管线默认主输入
│   ├── KEEP_len_filter_table.tsv  keep_summary.tsv（只读，生产者在外部）
│   ├── acvirus_classify/final_result_with_confidence.tsv
│   ├── acvirus_identify/identification_result_score.csv
│   ├── acvirus_trees/<Family>/{keep_<Family>.fasta, *.treefile}
│   └── SDT_matrix/<genus>/SDT_<genus>.{png,pdf}
│
├── 13_Report/                              # report_pipeline.py
│   ├── pipeline_report.html  stage_summary.tsv  directory_tree.txt
│   ├── {data,hostdep,assembly,ident,filter,tool_filter,cobra,checkv}_summary.tsv
│   ├── cluster_size_distribution.tsv  cluster_pipeline_reduction.tsv  cluster_skip_by_host.tsv
│   ├── taxonomy_{fill_stats,agreement_stats,consistency_summary,composition}.tsv
│   ├── final_integrated_classification.tsv  upset_data.json  taxonomy_images.json
│   ├── host_distribution.tsv  host_decision_method.tsv  host_ictv_{classification_summary,confidence}.tsv
│   ├── rescue_report.tsv  plant_virus_summary.tsv
│   ├── classification_sankey{,_plant}.html  taxonomy_sunburst.html
│   ├── {HQ_plant_viruses_info, All_plant.viruses_info, …}.tsv
│   ├── rescue_evidence_scored.tsv  cdd_evidence_report.tsv  hmm_ct3_hits.tsv
│   ├── keep_summary.tsv  frequency_table.tsv  rescue_detection_summary.tsv
│   └── logs/{clean,deplete,assembly,identification,filter,graph,cobra,merge,cluster,
│              taxonomy,host,checkv,rescue,analysis,analysis_verify,report}.log
│
└── （eve_distinguish 判别产物 2026-10-07 起归口 06_EVE/eve_distinguish/，见 §8 与 §13-D3；
     存量运行与 -D 指向非 03_Discovery 目录名时仍落 <DISCOVERY>/08b_EVE_Distinguish/；
     ⑥ 已归档 2026-10-11，此判别不再新跑）
```

**③→④ 桥**（运行 `discovery2analysis.py --from-discovery 03_Discovery --output_prefix 90_Handoff/analysis`）：

```
90_Handoff/
├── analysis.reference.fasta   # centroids 按 min_length 过滤后 fasta（contig_id=头首 token）
├── analysis.ref_info.tsv      # 17 列: Accession Taxid Species_NCBI Species_ICTV Segment Molecule_type
│                              #   Molecule_Type2 Genome_Length Taxonomic_Completeness
│                              #   Realm Kingdom Phylum Class Order Family Genus Species
│                              #   （Taxid=9000000+md5 派生 9,000,000–9,899,999）
└── handoff_manifest.tsv       # 台账: boundary key path produced_by produced_at（幂等主键=(boundary,key)）
```

---

## 6. ④ 已知病毒分析 — `04_Analysis/`

**编排**：`auto_known_virus.py --reference <90_Handoff>/analysis.reference.fasta --ref_info <90_Handoff>/analysis.ref_info.tsv --reads_dir <②>/02_HostDepleted --output_dir <项目>/04_Analysis --stage all`
**stage 顺序**：`detect → filter → variants → post → full(assembly) → extract → similarity → dvg → report`
**断点**：各阶段目录 `run_{1..9}_{stage}.log` + `run.ok`（Stage9 report 不写 run.ok，每次重跑）

```
04_Analysis/
├── known_virus.log                         # 编排器根日志
├── tmp/                                    # 编排器 TMPDIR
│
├── 01_Detection/                           # batch_virus_depth.py（salmon/bowtie2 检出）
│   ├── bam/{sample}.sorted.bam  stat/{sample}.chr.stat.gz  index/  batches/batch_0001.parquet
│   ├── pseudo_quant/{sample}/abundance.tsv          # 仅 salmon/kallisto
│   ├── logs/{pipeline_*.log, sample_{s}.log, sample_resource_usage.tsv}
│   ├── plots/virus_analysis/freq_multi_metrics_log10.{png,pdf}
│   ├── sample_distribution/                # 编排器 --mode sample 出图
│   ├── coabundance/                        # --mode coabundance 出图
│   ├── depth_plots/{sample}/{sample}_{tax}_{virus}_depth.{png,pdf}   # --mode depth 出图
│   └── summary/
│       ├── all_viruses.raw.tsv             # 全量原始命中
│       ├── all_viruses.summary.tsv         # 泊松/覆盖过滤后
│       ├── all_viruses.best.summary.tsv    # ★★★确诊白名单（④→⑧ 交接物；STAGE1_SCHEMA 校验:
│       │                                   #   Sample,Rep_Accession,Rep_Coverage(%),Rep_MeanDepth,
│       │                                   #   Asm_TPM,Poisson_Ratio,taxid,Adjusted_Species + 扩展列）
│       ├── all_viruses.unclassified.tsv    # 新种（ANI<95%）
│       └── analysis_report.txt
│
├── 02_Filtering/                           # utils/filter_summary.py
│   ├── high_conf.summary.tsv               # ★通过过滤（Stage3/4/6/8 的 get_summary() 首选输入）
│   ├── high_conf.discarded.summary.tsv  filter_stats.per_sample.tsv  filter_stats.per_virus.tsv
│   └── filter_summary_plot.{pdf,png}       # --plot
│
├── 03_Variants/                            # batch_virus_variants.py
│   ├── virus-fasta/ref_{safe_virus}/ref_{safe_virus}.ref.fasta    # 逐病毒归一参考
│   ├── virus-index/harmonized_references.fasta + harmonized_bowtie2*
│   ├── virus-bam/{sample}.sorted.bam
│   ├── virus_reads/{sample}.{safe_virus}*.fastq.gz          # --extract_reads 靶向 reads
│   ├── virus-consensus/{L1}/{L2}/{L2}.consensus.fasta + {L2}.fixed.bam
│   │       # L1={safe(tax)}_{safe(virus)}; L2={safe(sample)}_{safe(virus)}
│   ├── virus-variants/{L1}/{L2}/{L2}.{variants,filtered}.vcf + {L2}.fixed.bam
│   ├── virus-SnpEff/{L1}/{L2}/{L2}.ann.vcf + {L2}.annotation_summary.tsv
│   ├── virus-SNPGenie/…  virus-annotations/（GTF/GBK）
│   ├── summary/all_summary.tsv             # ★Stage4/5/6 的 summary 输入
│   ├── summary/Coinfection_Matrix_Reads.tsv
│   ├── summary/individual_virus_reports/{Sample}_{Taxonomy}.summary.tsv
│   └── logs/master_commands.log + {L2}_*.log
│
├── 04_PostAnalysis/                        # post 阶段（每病毒四件套）
│   ├── {Species_safe}_{Acc}/
│   │   ├── vcf_viz/                        # virus_variants_analyzer.py
│   │   │   ├── Table1_Unified_HighQuality_Variants.csv
│   │   │   ├── Table2_Sliding_Window_PopGen_Stats.csv
│   │   │   └── Figure1A…Figure9(_PCA_Lineages)(_3D).{png,pdf} + *_Data.csv
│   │   ├── vcf_merge/                      # virus_vcf_pipeline.py
│   │   │   ├── {prefix}.vcf.gz  stats/  snpgenie_merged/
│   │   │   ├── matrices/{snp_matrix.tsv, af_matrix.tsv, distance_*.tsv}
│   │   │   ├── trees/{newick, dendrogram.png}
│   │   │   └── figs/{pca.png, distance_clustermap.png, afs.png, ld_r2_*.png, …}
│   │   ├── snpeff_macro/                   # snpeff_analysis.py --miner
│   │   │   ├── Matrix_01_Global_Melted_Variants.csv … Matrix_06_Ultimate_OncoPrint_State.csv
│   │   │   └── Figure_1..4.{pdf,png}  Figure_5_Ultimate_Quasispecies_OncoPrint.pdf
│   │   ├── maftools/                       # snpeff2maf.py + viral_maftools.R
│   │   └── snpgenie/                       # snpgenie_master.py（Fig01…Fig11b + Matrix_ML_Clusters_3D.csv）
│   ├── sra.list
│   └── metadata_association/               # ★编排器写在这里（v2: figures/fig01…fig07 + tables/；
│                                           #   v1: 01_Global_Summary/…+ Viral_Infection_with_Metadata_Summary.tsv）
│                                           #   ⚠ 报告读的是 02_Filtering/metadata_association，见 §13-D2
│
├── 05_Assembly/                            # batch_virus_full.py → virus-full.py
│   └── {NormTax}_{safe_virus}/{safe_sample}_{safe_virus}/     # + 同名 _assembly.log
│       ├── 1.DeNovo_Assembly/ … 10.Gap_Filling/              # 12 个编号步骤目录
│       ├── 11.Ultimate_Circular_Result.fasta     # ★extract 目标产物
│       ├── 12.Coverage_Visualization/  13.Evolution_Stats.tsv  Global_Evolution_Stats.tsv
│       ├── Fully-assembled.ok
│       └── logs/
│
├── 06_Extraction/                          # utils/extract_full_fasta.py
│   ├── extraction_status.tsv               # virus/sample/status/reason/seq_len/n_pct
│   ├── assembly_stats.{pdf,png,csv}        # --plot
│   └── {病毒目录名}/                        # = 05_Assembly 一级目录名原样复制
│       └── {sample 中 '_'→'.'}.full.fasta  # ★★★④→⑤ 交接物（如 CRR527042_MN386961.1 →
│                                           #   CRR527042.MN386961.1.full.fasta）
│
├── 07_Similarity/                          # virus_auto_pipeline.py（逐病毒）
│   └── {vname}/                            # vname = 03_Variants/virus-consensus 子目录名
│       ├── 01_extraction_summary.tsv
│       ├── 01_extracted_genes/{prefix}_{feature}_{nucl,prot}.fasta
│       └── Mode_{Strict|Filter|Fill}/…/02_similarity_matrices/{prefix}_{align_method}.csv,
│         {prefix}_heatmap.{pdf,png}, 03_deduplication_reports/…, 04_filtered_sequences/*.fasta
│
├── 08_DVG/                                 # batch_virema_dvg.py + virema_summary_report.R
│   ├── individual_refs/{safe(acc)}.fasta
│   ├── virema_results/{safe(tax)}/{safe(sample)}_{safe(virus)}/   # ViReMa .sam + bed_results/
│   ├── all_gathered_beds/  r_summary_run.log
│   └── Summary_Analysis_Report/
│       ├── Aggregated_Sequence_Information.csv  Matrix_Sample_Wise_Recombination_Statistics.csv
│       └── Virus_Specific_Plots/{Acc}/{prefix}_{SV_Profile_Matrix.csv, Scatterplot.pdf,
│                                   Deletions.pdf, Duplications.pdf, Circos_4Track.pdf}
│
└── 09_Report/                              # generate_pipeline_report.py
    ├── Pipeline_Summary_Report.html        # ★报告文件名
    ├── stage_funnel/{stage_funnel_summary.tsv, virus_stage_matrix.tsv, *.png|pdf}
    │                                       #   由 stage_funnel_summary.py 生成
    └── run_9_report.log
```

**panvirome 聚合层**（不在 04_Analysis 树内，作用于多物种聚合根 `--base`）：

```
<聚合根>/                                 # --base（多物种项目共用目录）
├── prevalence_{sp}_{tag}/                 # 每物种一个 ④ 根（panvirome_prevalence 读它 + 对应 ③ 根）
├── {prefix}_{sp}_out/                     # 每物种一个 ③ 根（读 07_Taxonomy/10_Rescue/09_CheckV）
├── prevalence_full_table.tsv              # panvirome_prevalence.py --out
├── novel_vs_known_final.tsv               # panvirome_novel_known.py --out（中间物 all_rescue_centroids.fasta、
│                                          #   blastn/mmseqs best_hit.tsv、cdd_hits.tsv）
├── novel_vs_known_report.html             # panvirome_novel_known_report.py
└── Merged_Results_All_Viruses/{virus}/    # panvirome_aggregator.py --outdir
    ├── {gene_coordinates.csv, merged_allele_freq.csv, merged_snpeff.maf,
    │   merged_product_results.csv, merged_site_pi.csv}
    ├── Plot1_Mutation_Density.pdf … Plot5_Diversity_Sliding_Window.pdf   # panvirome_suite.R
    └── Plot6_Super_Evolutionary_Landscape.{pdf,png}                      # panvirome_landscape.py
```

---

## 7. ⑤ 系统发育/进化 — `05_Phylo/`

**编排**：`phylo_pipeline.py --extract_dir <④>/06_Extraction --metadata <①>/02_Unified/Global_Unified_Metadata_Core14.csv --output_dir <项目>/05_Phylo --stage all`
**19 个细粒度 stage**：`prep,metadata,online,clean,align,splitstree,rdp5,tree,popgen,host,capheine,clock,beast,gene_dating,phylogeo,geo,geo_paths,report`（`--stage` 大组名 = 8 模块）
**断点**：`<virus>/.checkpoints/<stage>.done`（JSON 含输入指纹；align/tree 有指纹校验）

```
05_Phylo/
├── report/                                 # batch 模式根级汇总（main 末尾无条件写）
│   ├── phylo_summary.csv                   # Virus/Phase/Success/Incomplete/Stages_* + headline 列
│   └── phylo_report.html                   # report_builder（单文件零依赖）
└── <vkey>/                                 # vkey = ④ 病毒目录名的 accession（正则提取，否则全名）
    ├── .checkpoints/<stage>.done
    ├── 01_data/
    │   ├── combined.fasta                  # prep：样本+参考（reference description）
    │   ├── dates.csv                       # name,date[,location]（时间通道）
    │   ├── sample_metadata.csv             # name,date,location（地理通道）
    │   ├── host.csv                        # name,host,tissue（宿主通道，仅 has_host 样本）
    │   ├── combined_dedup.fasta + dedup_report.csv    # --virus 模式去重命中时
    │   ├── metadata/{metadata_std.csv, metadata_report.tsv, metadata_summary.txt, channels.json}
    │   │                                   # metadata_std.csv 是下游唯一元数据输入
    │   ├── clean/{clean.fasta, removed.fasta, clean_report.tsv}
    │   ├── ncbi_ref/{<ACC>.gbk, <ACC>.ref.fasta}
    │   └── ncbi_harvest/{ncbi_sequences.fasta, ncbi_metadata.csv,
    │                      seqgroup_location.{pdf,csv}, seqgroup_host.{pdf,csv}}   # 仅报告用
    ├── 02_phylogeny/
    │   ├── mafft.aln.fasta                 # align（rdp5 开启时 tree 吃 04_recomb 的 realn 产物）
    │   ├── iqtree.treefile (+ .iqtree/.log/.ckp.gz/iqtree.support.json)
    │   ├── tree_circular.{pdf,png}
    │   ├── saturation_out/{saturation_report.csv, saturation_plot.png}
    │   ├── qc/{align_qc_report.tsv, clean.fasta, removed.fasta, align_qc.{pdf,png}}
    │   ├── splitstree/neighbornet.nexus
    │   ├── cross_validation/raxml.raxml.{bestTree,support,log}
    │   └── gene_analysis/{alignments/<gene>.mafft.fasta, trees/<gene>/, drt/<gene>/, ltt/<gene>/}
    ├── 03_popgen/
    │   ├── popgen_report.txt               # pypopart（π/θ/S/Hd/Tajima/Fst）
    │   ├── popgen_metrics.{pdf,png}  location_pops.tsv
    │   ├── distance/{distance_matrix.csv, nj_tree.nwk, distance_heatmap.png, nj_tree.png}
    │   ├── host_analysis/{tree_host.pdf, host_differentiation.pdf}
    │   └── dnasp/{dnasp_results.json, results.tsv, ld_pairs.tsv, figures/*.png}
    ├── 04_recomb/rdp5/
    │   ├── <virus>.csv                     # RDP5 七方法事件表
    │   ├── <virus> breakpoint distribution.csv  rdp5_events.json  recombination_report.md
    │   ├── dropped/{dropped_recombinants.fasta, drop_report.json}
    │   ├── rdp5_dropped.realn.fasta        # 默认策略产物（tree/clock/beast 实际消费）
    │   ├── masked/masked_N.fasta + mask_report.json    # --rdp5_mask 策略
    │   └── validation/{master_summary.tsv, refined_events.csv, events.rdp5ML,
    │                  run_all.sh|.ps1, <event>/{region_*.fasta, iqtree/, virusrecom/}}
    ├── 05_select/capheine/
    │   ├── ref_cds.fasta                   # gbk_extractor 提取
    │   ├── removeterminalstopcodon/ seqkit/ cawlign/<gene>-aligned.fasta
    │   ├── removeambigseqs/<gene>-clean.fasta  iqtree/
    │   ├── hyphy/CLN/<gene>-nodups.fasta  hyphy/LABELTREE/
    │   │   hyphy/{FEL,MEME,PRIME,BUSTED,CONTRASTFEL,RELAX}/<gene>.<TOOL>.json
    │   ├── drhip/{combined_sites.csv, combined_summary.csv, hyphy/}   # checkpoint 标志物
    │   ├── multiqc/{capheine_summary.html, multiqc_report.html}
    │   ├── selection_out/{selection_results.json, fubar/<gene>.FUBAR.json + fubar_summary.json,
    │   │                  relax/, contrastfel/}
    │   └── codon_miner/{<gene>_Protein_Site_<site>_Evolution.png + .csv,
    │                    selection_per_gene.{pdf,png}}
    ├── 06_time/
    │   ├── clock_summary.tsv               # 整体+分基因 β/R²/p/DRT 并排
    │   ├── treetime_rtt/{RootToTip_regression.pdf, TimeTree.pdf, timetree_inferred.nwk,
    │   │                 robust/{robust_report.json, robust_rtt_compare.{pdf,png}}}
    │   ├── treedater_ltt/{Phylogeny_dated.pdf, LTT.pdf, outlier_tips.tsv,
    │   │                  timetree_pruned.nwk, dates_decimal.csv}     # clock checkpoint 标志物
    │   ├── temporal_signal/{drt_summary.csv, drt_results.{json,pdf}, real/, randomized_NN.*,
    │   │                    beast_xmls/, _drt_ref.treefile}
    │   ├── beast/
    │   │   ├── beast1.xml (+ beast11.xml beast12.xml …)   # BEAST1 UCLN+skyline
    │   │   ├── beast1.log/.trees、beast11.log/.trees、*_run.log、run_status.json
    │   │   ├── beauti/{beauti.fasta, beauti_traits.txt}
    │   │   ├── merged/{merged.log, merged.trees, merged_sampled.trees, mcc.tree,
    │   │   │           parameter_summary.json, parameter_summary.csv, mcmc_diagnostics.pdf}
    │   │   └── postprocess/{ess_results.csv, BSP.pdf, RSPP_<chart>.pdf}
    │   ├── gene_split/{<gene>.fasta, _ref_cds/, _tmp/}    # gene_align_generator fallback
    │   └── gene_dating/{<gene>/{<gene>.xml,.log,.trees}, gene_dating_watch_parse.sh,
    │                     parse_gene_beast.py, gene_tmrca_summary.tsv, gene_tmrca_forest.{pdf,png}}
    ├── 07_geography/
    │   ├── prior_comparison/               # --phylogeo_prior auto
    │   ├── phylogeography/
    │   │   ├── phylogeo_beast1.xml + phylogeo_sky{1..N}.{xml,log,trees,location.rates.log}
    │   │   ├── *_run.log  run_status.json
    │   │   └── merged/{merged.log, merged.trees, mcc.tree, parameter_summary.json,
    │   │                 mcmc_diagnostics.pdf}                    # phylogeo checkpoint=mcc.tree
    │   ├── visualization/{spread3.json, phylogeo_map.html, migration_bf.csv,
    │   │                  phylogeo_multipanel.pdf, phylogeo_skyline.html,
    │   │                  transmission_map.{png,pdf,svg}, transmission_map.gif,
    │   │                  transmission_tree_map.{png,pdf,svg}}
    │   ├── tempmig/{migration_over_time.csv, migration_yearly_summary.csv, migration_over_time.pdf}
    │   └── geo_analysis/
    │       ├── mantel_test.pdf             # geo checkpoint 标志物
    │       ├── geo_distance.csv  genetic_distance.csv  tree_geography.pdf
    │       ├── VirSpaceTime.pdf  spatiotemporal_counts.csv
    │       ├── pathways/{branch_table.csv, pathways.csv, episodes.csv, ltl_summary.csv,
    │       │            pathway_qc.json}                          # geo_paths checkpoint
    │       ├── dispersion/{velocity_table.csv, wavefront_series.csv, dispersion_summary.json}
    │       ├── markov_jumps.json  timescale_consistency.json
    │       ├── subsampled/                 # --geo_subsample
    │       └── extended/{rrt_results.{pdf,csv}, rrt_plot.pdf, rrt_tree/{…,rrt_tree_table.csv},
    │                      tempmig.pdf, tempmig_summary.csv, genetic_subsample_sensitivity.csv}
    └── 08_report/
        ├── phylo_summary.csv + phylo_report.html   # 仅 --virus 单病毒模式（batch 模式在根 report/）
        ├── sample_plot/{sample_temporal.pdf, sample_map.pdf}
        └── publication/{Fig1_temprest.pdf, Fig2_mcc_tree.pdf, FigS1_parameters.{csv,tex}}
```

---

## 8. ⑥ 内源性病毒 — `06_EVE/`【已归档 2026-10-11，命令不可运行；下为归档前产物形态】

**编排**：`eve_screen.py -g <宿主基因组>.fasta -n <NAME> -o <项目>/06_EVE -t 40`（批量 `-B batch.tsv`，TSV 两列 `NAME<TAB>路径`）
**stage 顺序**：`discover(1) → verdict(2) → annotate(3) → [viroid 可选] → merge(4)`；断点=按产物存在性判定（空文件=已完成零结果）

```
06_EVE/
├── logs/{eve_screen.log, <NAME>.log, eve_element_db.log}
├── 01_Loci/<NAME>/                         # Stage1 滑窗 blastx 发现
│   ├── <NAME>.fna                          # 归一化/解包后基因组（tar.gz 自动挑 genome）
│   ├── <NAME>.chunks.fna                   # seqkit sliding -W 50000
│   ├── <NAME>.s1_raw.tsv → <NAME>.s1_hits.bed → <NAME>.loci.bed   # diamond→坐标还原→300nt 合并
│   ├── <NAME>.loci.fa + <NAME>.loci.fa.regions.txt + <NAME>.genome.fai
│   ├── <NAME>.s1_best.tsv                  # 每位点最佳命中（locus ref_sseqid family …）
│   ├── seqkit.err  diamond_s1.err
│   └── (chunks/s1_raw/s1_hits 被 cleanup 规则按下游产物逐步清理)
├── 02_Verdict/<NAME>/                      # Stage2 植物/病毒双侧 bitscore 判定
│   ├── <NAME>.s2_raw.tsv
│   └── <NAME>.s2_verdict.tsv               # locus verdict best_viral_{id,bs} best_plant_{id,bs}
├── 03_RVDB/<NAME>/                         # Stage3 RVDB 深度归属
│   ├── <NAME>.cand3.{bed,fa}
│   └── <NAME>.s3_rvdb.tsv                  # 空文件=已完成零候选（断点标记）
├── 04_Summary/
│   ├── <NAME>_eve_summary.tsv              # ★12 列汇总（断点标记）
│   ├── <NAME>_viroid.tsv                   # --blastn-viroid
│   └── eve_report.html                     # eve_report.py（--merge 后）
├── 05_Elements/                            # eve_element_db.py 元件级后处理
│   ├── all_candidates.fa  candidates_meta.tsv  element_table.tsv  element_regions.tsv
│   ├── eve_consensus_library.fa  family_quant.tsv  consensus_domains.tsv
│   ├── cluster_members.tsv  _clu_*(_seq.fasta, _tmp/)  _members_C%06d.txt
│   ├── rm/<NAME>/<NAME>.out  _cdd_mm/  phylo/<家族>_reps.fa + phylo_commands.sh
├── qc_contamination.tsv + QC_CONTAMINATION.md      # eve_qc_contam.py 事后 QC
├── kingdom_summary.tsv                     # merge：genome viral_supported host_like undetermined
├── kingdom_loci_all.tsv.gz                 # merge：13 列全位点
├── family_by_genome.tsv                    # merge：family × genome 计数矩阵
└── eve_distinguish/                        # ③↔⑥ 后运行判别（run_all.sh -D 指向 03_Discovery 时归口这里）
    ├── panel.fasta + panel.dmnd  baits.dmnd
    ├── panel_hits.tsv  baits_hits.tsv  all_hits.tsv
    ├── s1_decay.tsv  s2_domains.tsv  locus_blastn.tsv  locus_errors.log
    ├── locus_architecture.tsv  (refgenome_evidence.tsv + refgenome_blastn.tsv 可选)
    ├── eve_distinguish_verdict.tsv         # 22 列裁决表
    └── dna_vs_eve_filter.tsv               # ★最终 15 列过滤表（CARRY 13 列 + action + action_reason）
```

---

## 9. ⑦ 提交管线（无独立输出根）

**编排**：`submission_pipeline.py --work-dir <项目>/03_Discovery --run-title <title> --mode both`
产物落在 **`work_dir` 内部**（中间物）+ **`work_dir.parent/submission_<run_title>/`**（最终提交包）：

```
03_Discovery/                              # --work-dir（中间物直接写进 ③ 根）
├── submission_pipeline.log
├── suvtk.taxonomy_output/{taxonomy.tsv, miuvig_taxonomy.tsv, miuvig_taxonomy_enriched.tsv,
│                          segmented_viruses_info.tsv, taxresults_lca.tsv}
├── suvtk.features_output/{featuretable.tbl, proteins.faa, proteins_updated.faa,
│                          reoriented_nucleotide_sequences.fna, miuvig_features.tsv,
│                          alignment.m8, tophit_info.tsv}
├── analyze_hypothetical/{featuretable_updated.tbl, proteins_updated.faa, merged_blast.txt}
├── srr_list.txt                           # --fetch-metadata（taxonomy.tsv 提取 [SC]RR\d+）
└── metadata_output/{Global_Unified_Metadata_Full.tsv, Core14.tsv}   # gsa_sra.info.py -m local

<project_root>/submission_<run_title>/     # sub_out = work_dir.parent / f"submission_{run_title}"
├── suvtk_submission/
│   ├── miuvig.tsv  assembly.tsv  comments.cmt  template.sbt
│   ├── source.src                         # 列: Sequence_ID Organism Isolate Collection_date
│   │                                      #   geo_loc_name Lat_Lon Bioproject Biosample SRA
│   │                                      #   Segment Metagenomic Metagenome_source Host
│   ├── virus-annotations/*.gb             # tbl2gb.py
│   └── submission.sqn                     # ★table2asn 终产物（出前过 validate_submission_gate）
├── topology.tsv                           # viral_topology.py（DTR/ITR → circular/linear）
├── <run_title>/…                          # sequin 模式（sequin_builder.py）
├── unified_metadata.csv (+.bak)  source.src  source_individual/source_<virus>.src
├── biosample_template.tsv  miuvig.tsv  assembly.tsv  submission_log.csv
├── report.html                            # report_html.py
└── id_norm/{source_{filled,norm}.src, comments_*, featuretable_*, sequences_*}
                                           # fix_sqn.py / sync_sqn_from_csv.py（服务器同步链路）
```

> ⑦ 输入侧：默认吃 `12_AnalysisVerify/class_KEEP.fasta`（缺则回退 `10_Rescue/HQ_plant_viruses.fasta`）；
> 物种/宿主来自 `11_ViromeAnalysis/ref_info.tsv`（或 90_Handoff 桥表）；元数据来自 ① Full/Core14；
> 外部依赖 suvtk + table2asn + RVDB v31。

---

## 10. ⑧ 宿主表达（可选模块）— `08_HostExpression/`【已归档 2026-10-11，命令不可运行；下为归档前产物形态】

**编排**：`host_expression.py --project-root <项目根>`（四条输入自动定位，standard 优先 legacy 兜底）
**执行序**：`hostref → quant → qc → labels → coexp → de → enrich → loadscan → report`
（⚠ 目录号与执行序不同：**07_Enrich 先于 08_LoadScan 执行**）
**断点**：`.checkpoints/<stage>.ok` 且 `STAGE_PRODUCTS[stage]` 产物齐备

```
08_HostExpression/
├── host_expression.log  .checkpoints/{<stage>.ok, coexp_<子步>.ok}
├── 01_HostRef/host_cds.fasta               # gffread -x（genome+GFF3 → CDS）
├── 02_Quant/                               # salmon_auto_pipeline.py
│   ├── quant/{sample}/quant.sf
│   ├── matrix/{TPM,NumReads,FPKM,EffectiveLength}_matrix.tsv + sample.list
│   ├── metadata/{sample_info.tsv, run_summary.json, missing_read_samples.tsv, failed_samples.tsv}
│   ├── logs/  workflow.log  salmon_index/
├── 03_QC/{qc_summary.tsv, qc_passed_samples.txt}    # sample num_processed num_mapped percent_mapped pass note
├── 04_Labels/{virus_labels.tsv, sample_labels.tsv, virus_columns.tsv}
│       # virus_labels: Sample Virus Covered_pct MeanDepth Load Reads Positive
│       # sample_labels: Sample Virus_any(+/-) N_virus + 每病毒 <sanitize>_positive/_load + Tissue 等
├── 05_Coexp/                               # TEA-GCN Δz 主分析（七子步）
│   ├── matrices/{positive.tsv, negative.tsv, group_summary.tsv}
│   ├── teagcn_{positive,negative}/{Generate_partitions/, RAW_GCN/, Completed_GCN/}
│   ├── network_qc.tsv                      # Group … Max_z … Status(PASS≥1.7/LOW≥1.5/FAIL)
│   ├── diff_coexp_{candidates,significant,gene_aggregate,summary}.tsv   # ★断点产物=significant+summary
│   └── modules/{module_table_<组>.tsv, module_edges_<组>.tsv, module_summary.tsv,
│                module_load_assoc.tsv}
├── 06_DE/{deseq2_all_results.tsv, deseq2_degs.tsv, volcano_plot.png, deseq2_summary.tsv,
│          deseq2_strata/<组织>.tsv}
├── 07_Enrich/{sets/{DE.txt, DE_up.txt, DE_down.txt, dz_hub.txt, module_<组>_<id>.txt},
│              background_genes.txt, enrich_<Set>.tsv, enrich_summary.tsv}
├── 08_LoadScan/{<病毒>_genome_scan.tsv, <病毒>__<组织>_genome_scan.tsv, load_scan_summary.tsv}
└── 09_Report/{host_expression_report.html, figures/fig01_delta_z_distribution.… fig05_circular_pathway.…}
```

---

## 11. 90_Handoff（交接层，两布局统一）

```
<project_root>/90_Handoff/
├── analysis.reference.fasta               # ③→④ 桥（discovery2analysis.py）
├── analysis.ref_info.tsv                  # 17 列（见 §5 末尾）
└── handoff_manifest.tsv                   # 边界台账: boundary key path produced_by produced_at
    # 现有写入方:
    #   discovery2analysis.py  → boundary=discovery->analysis, keys=centroids,taxonomy,taxonomy_files,reference,ref_info
    #   host_expression.py     → boundary=hostexpr, keys=reads_dir,reference,detect_table,tpm_matrix,sample_labels,deseq2_all_results
    #   ① generate_report.py   → 独立实现（不调 io_layout）
```

---

## 12. 边界契约一览（生产者 → 消费者）

| # | 边界 | 产物（standard 实际路径） | 消费入口 |
|---|---|---|---|
| B1 | ①→② | `01_PublicData/03_RawData/{acc}_{1,2}.fastq.gz` + `05_HostRef/hostdb/` | `data_preprocessing.py --input_reads --host_db` |
| B2 | ②→③ | `02_Preprocessing/02_HostDepleted/{sample}_clean_{1,2}.fa.gz`（+ `assembly_ready.list`） | `virome_pipeline.py --input_reads` |
| B3 | ③→④ | `03_Discovery/06_CLUSTER/4_centroids/final_centroids.fasta` + `07_Taxonomy/{Votus}.integrated/final_integrated_classification.tsv` → 桥 → `90_Handoff/analysis.reference.{fasta,tsv}` | `auto_known_virus.py --reference --ref_info --reads_dir <②>/02_HostDepleted` |
| B4 | ④→⑤ | `04_Analysis/06_Extraction/{病毒目录}/{样本}.{acc}.full.fasta` | `phylo_pipeline.py --extract_dir` |
| B5 | ①→⑤ | `01_PublicData/02_Unified/Global_Unified_Metadata_Core14.csv` | `phylo_pipeline.py --metadata`（data_collector 原生直读） |
| B6 | ①→⑥ | `01_PublicData/05_HostRef/genome/all.genome.uniq.fasta` | `eve_screen.py -g`（或 `-B batch.tsv`）（**⑥ 已归档 2026-10-11，边界停用**） |
| B7 | ③↔⑥ | ⑥ 判别器吃 ③ 的 `evescan_query.fasta`/`rescue_evidence_scored.tsv`，产出 `dna_vs_eve_filter.tsv` 归口 `06_EVE/eve_distinguish/`（-D 目录名非 `03_Discovery` 时旧默认 `<③根>/08b_EVE_Distinguish/`） | `eve_distinguish/run_all.sh -D <③根> -A <1kp 组装根>`（**⑥ 已归档 2026-10-11，边界停用；历史产物仍可被 `triage --preset eve_distinguish` 消费**） |
| B8 | ②①④→⑧ | ①`05_HostRef/genome`(+GFF3)、①`02_Unified/Core14`、②`01_CleanData/{3.clumpify,2.fasta}`、④`01_Detection/summary/all_viruses.best.summary.tsv` | `host_expression.py --project-root`（四条自动定位）（**⑧ 已归档 2026-10-11，边界停用**） |
| B9 | ③→⑦ | `12_AnalysisVerify/class_KEEP.fasta`（主）→ `10_Rescue/HQ_plant_viruses.fasta`（回退）+ `11_ViromeAnalysis/ref_info.tsv` | `submission_pipeline.py --work-dir <③根>` |
| B10 | ③④→panvirome | `{base}/prevalence_{sp}_{tag}/`（④ 根）+ `{base}/{prefix}_{sp}_out/`（③ 根） | `panvirome_prevalence.py --base` |

---

## 13. 重扫发现的契约/实现偏差清单（2026-10-07）

重扫确认以下几处 **代码现状 ≠ 设计文档/注册表**，均已核实到行号：

| # | 位置 | 偏差 | 影响 |
|---|---|---|---|
| D1 | `virome_discovery_pipeline/cluster_pipeline.py:617` | cluster 实写 `06_CLUSTER/4_centroids/`，而 io_layout standard 契约名是 `06_CLUSTER/centroids/` | 编排器靠 `_locate_centroids_dir()` 重绑 + `locate_centroids()` 的 `*/*` 限深兜底命中，主链不断；但按契约路径直接取文件会 miss |
| D2 | ~~`auto_known_virus.py:735,750` vs `generate_pipeline_report.py:837,2142`~~ **已解决（2026-10 核实）** | 原：post 阶段把 metadata_association 写到 `04_PostAnalysis/metadata_association`，报告只读 `02_Filtering/metadata_association`。现状：报告三处触点均已改为「`04_PostAnalysis/metadata_association` 优先、`02_Filtering` 兜底回退」（`generate_pipeline_report.py:937-942` Overview 指标、`:976-982` S2 阶段文件、`:2306-2310` 拷贝归档），错配消除 | 报告可正常取到 metadata 汇总表；本条仅存档 |
| D3 | ~~`eve_distinguish/run_all.sh:165`~~ **已解决（2026-10-07）** | 原：判别产物落 `03_Discovery/08b_EVE_Distinguish/`，与文档 B7 的 `14_EVE_Distinguish/` 互不相认。现已按用户裁定归口 **`06_EVE/eve_distinguish/`**（`-D` 目录名为 `03_Discovery` 即 standard 布局时生效；其余目录名与存量运行维持旧默认 08b，`-o` 始终可显式覆盖）。测试：eve_distinguish 21/21（⑥ 已归档 2026-10-11） | 判别产物归 ⑥ 名下，不再混入 ③ 的编号体系 |
| D4 | `virome_pipeline.py`（`TAX_GATE_REQUIRED="6.2"`）vs `virus_classifier_analysis.R` | 编排器强制校验 `taxonomy_gate_stamp.tsv`/`taxonomy_gate_check.tsv`，但当前树内的 R 共识脚本不写这两个文件（写闸门戳的版本在 `archive/20261007_audit_cleanup/snapshot/live_virus_classifier_analysis.R`，原 scripts/audit/snapshot 已归档） | taxonomy 阶段闸门恒判 FAIL → 每次重跑 R 共识并记 error 日志（不会静默放行） |
| D5 | `host_expression.py:69`（⑧ 已归档 2026-10-11） | ⑧ 执行序为 …coexp→de→**enrich→loadscan**…，目录号 07/08 与执行序相反 | 仅需知晓，非缺陷 |
| D6 | ② `data_preprocessing.py` `cleanup_00a_after_deplete()` | `--stage all` 跑完后 `01_CleanData/{1.fastp_tmp,2.fasta,3.clumpify}` 的样本数据文件被清理（留 logs/JSON） | ②→⑧ 的 `3.clumpify` 取 clean reads 仅在**只跑 clean 未跑 deplete 清理**时成立；跑完全链后 ⑧ 需从 `02_HostDepleted` 取（或保留 ② 单独跑 clean 的副本） |
| D7 | `scripts/audit/check_io_handoffs.py` SCOPES | 审计范围缺 `virome_submission_pipeline/`（它 import 了 layout_dirs） | ⑦ 的 legacy 字面量不在审计覆盖内 |

其余全部脚本核实**与契约一致**（②/④/⑤/⑥/⑧ 的目录名均经 `io_layout` 注册表解析，20/24 项 io_layout 冻结测试通过，`check_io_handoffs` 审计 0 FAIL）。
