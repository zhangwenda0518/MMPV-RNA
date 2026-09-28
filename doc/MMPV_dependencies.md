# MMPV-RNA 管线外部依赖清单

> 来源：`pipeline_config.yaml`、`pixi.toml`、`SOFTWARE_VERSIONS.txt`、`DATABASE_SETUP.md`、`envs/viralm.yaml` 及各管道脚本的 subprocess/import 扫描。
> 路径为服务器路径（`/home/zhangwenda/...`），本机映射为 `~/...`。
> 2026-09-06 复核：全部外部命令工具/Python/R 包与脚本实际调用逐项比对；pixi.toml 中 bioconda 不存在的包已移入手动安装区。

---

## A. 外部命令工具

### A1. 数据清洗 & 宿主去除
| 工具 | 用途 | 出现脚本 |
|---|---|---|
| `fastp` | QC + trimming | clean-data.py, preprocess.py, virome_pipeline.py |
| `seqkit` | FASTA/Q 操作（split2/grep/seq/stats） | clean-data.py, host_depletion.py, virus_identification.py, cluster_pipeline.py, build_virus_db.py |
| `clumpify.sh`(bbmap) | k-mer 聚类重排 (压缩/组装加速) | clean-data.py, preprocess.py |
| `bbnorm.sh/bbduk.sh/bbmap.sh`(bbmap) | 归一化/去污染/比对 | run_bbnorm.py, preprocess.py, host_depletion.py |
| `bowtie2` | 宿主比对（RNA-short） | host_depletion.py, batch_virus_depth*.py |
| `hisat2` | 宿主比对（splice-aware, 三选一） | host_depletion.py, pixi.toml |
| `minimap2` | 长读段/宿主比对 | host_depletion.py, pixi.toml |
| `kraken2` | 宿主分类去除 | host_depletion.py, build_virus_db.py |
| `ribodetector` | rRNA 检测去除 | pixi.toml |
| `seqtk` | subseq 提取 | host_depletion.py |
| `pigz/gzip` | 并行压缩/解压 | host_depletion.py, batch_virus_variants.py |
| `vmtouch` | Kraken2 库锁 page cache | host_depletion.py |

### A2. 组装
| 工具 | 用途 | 出现脚本 |
|---|---|---|
| `megahit` | 宏基因组组装 | assembly_pipeline.py, virome_pipeline.py |
| `rnaviralspades.py/spades` | 病毒引导组装 | assembly_pipeline.py, virome_pipeline.py |
| `penguin` | 引导组装（**源码自编译, bioconda 无包**） | assembly_pipeline.py |
| `refineC` | 组装优化（split/merge）（**源码编译, bioconda 无包**） | assembly_pipeline.py, virome_pipeline.py, rescue_pipeline.py |
| `flye` | 长读/共组装（可选, biosoft 源码版 2.9-b1779） | virome_pipeline.py, utils/flye_trace_native.py |
| `bioawk` | rnaviralspades 依赖 | assembly_pipeline.py |

### A3. 病毒鉴定（Stage 2）
| 工具 | 用途 | 出现脚本 |
|---|---|---|
| `diamond` | BLASTX 蛋白搜索 | virus_identification.py, build_virus_db.py, rescue_pipeline.py |
| `blastn/blastx` | 核酸/蛋白比对 | virus_identification.py, pyhmmer_viralverify.py |
| `makeblastdb/blastdbcmd` | BLAST 建库/取序列 | rescue_pipeline.py, build_virus_db.py |
| `genomad` | 病毒鉴定 | virus_identification.py |
| `virsorter` | VirSorter2 分类 | virus_identification.py |
| `viralverify` | 病毒序列验证 | virus_identification.py |
| `metabuli` | 宏基因组分类 | virus_identification.py, build_virus_db.py |
| `virhunter`(biosoft/virhunter) | DL 病毒检测(keras/tf) | virus_identification.py |
| `virbot`(biosoft/VirBot) | DL 病毒检测 | virus_identification.py |
| `viralm`(utils/viralm_cpu.py) | DNABERT-2 病毒鉴定(env viralm) | virus_identification.py |
| `rdrpcatch` | RdRp 检测(env rdrpcatch) | virus_identification.py |
| `checkv` | 完整性/质量评估 | virome_pipeline.py |
| `taxonkit` | NCBI TaxID 谱系 | virus_classifier.py |

### A4. 病毒分类（Stage 5，8 工具）
| 工具 | 用途 | 出现脚本 |
|---|---|---|
| `mmseqs`(easy-taxonomy/easy-search/createdb/createtaxdb/convertalis) | 快速分类搜索 + CDD 验证 | virus_classifier.py, validate_rescue_cdd.py, build_virus_db.py |
| `genomad` classify | 病毒分类 | virus_classifier.py |
| `metabuli` classify | 分类 | virus_classifier.py |
| `CAT`(CAT_pack) | 分类 | virus_classifier.py, classify_contigs.py |
| `VITAP`(assignment) | 病毒分类 | virus_classifier.py |
| `ACVirus`(classify/cli.py) | 科/属级分类 + 新病毒鉴定 | virus_classifier.py, virome_pipeline.py |
| `vcontact3`(run) | 基因共享聚类（**opt-in：默认 `--tools all` 不包含，需显式指定才跑**） | virus_classifier.py |
| `diamond_lca` | 蛋白 LCA | virus_classifier.py |
| `suvtk`(taxonomy) | ICTV 分类学注释 | virome_submission.py, virome_analysis.py |

### A5. 聚类（Stage 4）
| 工具 | 用途 | 出现脚本 |
|---|---|---|
| `vclust`(prefilter/align/cluster/deduplicate) | 病毒聚类去冗余(Leiden+cd-hit) | cluster_pipeline.py |
| `cd-hit` | 参考预聚类 | cluster_pipeline.py, virome_pipeline.py |
| `genome_rmDuplicates.pl`(Perl) | 组装 contig 去重 | cluster_pipeline.py |

### A6. COBRA 延伸（Stage 3，实际 = bwa-mem2 组装）
| 工具 | 用途 | 出现脚本 |
|---|---|---|
| `bwa-mem2` | 短读比对组装 | cobra_pipeline.py, rescue_pipeline.py |
| `samtools`(sort/index/view/mpileup/faidx/idxstats) | BAM 处理 | cobra_pipeline.py, rescue_pipeline.py |
| `coverm` | 覆盖度计算 | cobra_pipeline.py |

### A7. 宿主预测（Stage 6）
| 工具 | 用途 | 出现脚本 |
|---|---|---|
| `rnavirhost`(classify_order/predict) | RNA 病毒宿主（conda 1.0.5 + 用户修改版, mambaforge） | run_host_prediction.py |
| `phabox2`(--task cherry) | 噬菌体宿主（**GitHub 手动安装, bioconda 无包**） | run_host_prediction.py |
| `cherry` | CHERRY 宿主预测（**GitHub 手动安装, bioconda 无包**） | run_host_prediction.py |

### A8. Rescue / 共识
| 工具 | 用途 | 出现脚本 |
|---|---|---|
| `ragtag.py`(scaffold) | 参考引导 scaffolding | rescue_pipeline.py |
| `viral_consensus` | 共识序列 | rescue_pipeline.py, batch_virus_depth*.py |
| `blastdbcmd/makeblastdb/blastn` | BLAST 拯救 | rescue_pipeline.py |
| `bwa-mem2/bwa` | 比对 | rescue_pipeline.py |
| `salmon` | VSI 检测 | rescue_pipeline.py, Virseqimprover.py |
| `bedtools`(getfasta) | VSI scaffold 截取（**Virseqimprover 拼入 shell 命令真实执行, 硬依赖**） | Virseqimprover.py |

### A9. 变异检测 & 注释（Analysis）
| 工具 | 用途 | 出现脚本 |
|---|---|---|
| `ivar`(variants) | 病毒变异检测(默认) | batch_virus_variants.py, batch_virus_depth*.py |
| `freebayes` | 变异检测(备选) | batch_virus_depth*.py |
| `lofreq` | 低频变异检测 | batch_virus_variants.py |
| `bcftools` | VCF 处理 | batch_virus_variants.py |
| `snpEff.jar/SnpSift.jar`(JRE≥21) | 变异注释 | biosoft/snpEff/, snpeff_build.py |
| `snpgenie`(Perl) | 群体遗传 dN/dS | snpgenie_master.py, auto_known_virus.py |
| `pandepth` | 位点深度 | batch_virus_depth*.py |
| `bgzip/tabix` | VCF 压缩索引 | pixi.toml |

### A10. DVG / 重组
| 工具 | 用途 | 出现脚本 |
|---|---|---|
| `ViReMa.py`(biosoft/ViReMa) | DVG/重组检测 | batch_virema_dvg.py |
| `DI-tector_06.py`(biosoft/ViReMa) | DVG 检测(备选) | batch_virema_dvg.py |
| `bowtie`(v1, bowtie-build) | DVG 比对 | batch_virema_dvg.py |

### A11. 选择压力（CAPHEINE）
| 工具 | 用途 | 出现脚本 |
|---|---|---|
| `cawlign` | 密码子感知比对 | capheine_pipeline.py |
| `hyphy`(busted/fel/meme/fubar/relax/absrel/slac/cln) | 正选择分析 | capheine_pipeline.py, selection_suite.py |
| `iqtree` | ML 树 | capheine_pipeline.py |
| `drhip` | DRHIP 分析 | capheine_pipeline.py |

### A12. 系统发育 & 谱系地理（Phylo）
| 工具 | 用途 | 出现脚本 |
|---|---|---|
| `mafft` | 多序列比对 | phylo_pipeline.py, sdt_genus_matrix.py |
| `iqtree/iqtree2` | ML 建树 | phylo_pipeline.py, rdp5_validate.py |
| `beast`(v1/v2) | 贝叶斯定年+谱系地理 | phylo_pipeline.py, batch_phylogeo.py |
| `logcombiner/treeannotator` | BEAST 链合并/MCC | merge_results.py |
| `treetime`(RTT) | 根到尖回归定年 | phylo_pipeline.py |
| `treedater`(LTT) | 分子定年 | phylo_pipeline.py |
| `codeml`(PAML) | M7/M8 LRT 正选择 | codeml_bridge.py |
| `dnasp`(DnaSP 6 重实现) | 群体遗传 | dnasp_bridge.py |
| `rdp5`(逆向 rdp5_ini.py + rdp5.sh) | 重组检测(RDP/GENECONV/BOOTSCAN/MaxChi/Chimaera/SiScan/3Seq) | rdp5_ini.py, rdp5_validate.py |
| `tempmig`(TempMig) | 时序迁移追踪 | geo_analysis.py, virphy_bridge.py |
| `spread3` | 谱系地理可视化 | spread3_viz.py |
| `seqharvester/seqgrouper` | 序列采集/分组 | seqharvester_bridge.py, seqgrouper_bridge.py |
| `pypopart` | 单倍型网络 | build_haplo_outputs.py |
| `SDT`(sdt_genus_matrix.py 内调 mafft) | 属级序列一致性矩阵 | sdt_genus_matrix.py |
| `Rscript` | 运行 R 脚本 | ViPhyKit 等 |

> ✅ **pypopart 已固定到修复版 `4562cba`（2026-09-14，含 indel 修复）** —— 本地副本 `biosoft/pypopart/`
> （含 `PROVENANCE.md`），服务器运行环境 `~/MMPV-RNA/biosoft/pypopart`（git 仓库，2026-09-22 从
> `97d76b1` 升级至 `4562cba`）。
>
> **修复的 bug**：v0.1.0（≤ `97d76b1`, 2025-11-19）的 `core/haplotype.py: identify_haplotypes_from_alignment`
> 用 `remove_gaps()` 后的序列**同时**作分组 key 和单倍型代表序列，而下游 `core/distance.py: hamming_distance`
> 强制 `len(seq1) == len(seq2)` 否则抛 `ValueError`。病毒比对含 indel → 各序列去 gap 后长度不等 →
> **MJN / MSN / MST / TCS / TSW / parsimony 全部网络算法崩溃**（这些算法都 import 该函数）。
> 上游 `555ee04`（2026-09-07，作为 "Phase 2: shared infrastructure" 重构的一部分）改为按**含 gap 的完整
> 比对序列**分组（`key = seq.data`、`sequence_map[key] = seq`），与 PopART 的 `condenseSeqs` 语义一致，
> 保证所有单倍型共享比对长度。**升级验证**：含 indel 比对在旧版抛
> `ValueError: Sequences must have same length: 11 vs 10`，新版通过（10 节点建网成功）。
>
> **本仓库补丁 `virome_phylo_pipeline/patch_pypopart_haplotype.py` 的角色**：2026-08-28 的独立修复
> （早于上游 10 天），现降级为「无法升级时的权宜方案」，并且**分组语义与上游不同**——补丁保留去 gap
> 序列作 key，会把 gap 位置不同但碱基相同的序列归为同一单倍型（实测：该输入旧版+补丁得 1 个单倍型，
> 上游得 2 个）。**论文级结果请用上游版本**，保证与 PopART 可比。

### A13. 公共数据获取
`prefetch/fasterq-dump`(sra-tools) · `datasets`(ncbi-datasets-cli) · `esearch/efetch`(entrez-direct) · `aria2c` · `wget/curl`

### A14. 数据库构建（build_virus_db.py 专用）
`aria2c, rapidgzip, seqkit, SearchAccessionIdToTaxId.py, db_seqid2taxid_add_legth.py, make_ktaxonomy.py, metabuli build, centrifuger-build, kraken2-build, bracken-build, krakenuniq-build, kmcp compute/index, ganon build-custom, sylph sketch, kma index, salmon index, kallisto index, kaiju-mkbwt/mkfmi, diamond makedb, mmseqs createdb/createtaxdb, CAT_pack prepare, taxonkit`

> ⚠️ 3 个自写辅助脚本（SearchAccessionIdToTaxId.py / db_seqid2taxid_add_legth.py / make_ktaxonomy.py）已建收编目录
> `virome_discovery_pipeline/utils/db_build/`：`build_virus_db.py` 自动把该目录插到 PATH 最前，从服务器拷入即用（命令见其 README.md）。

### A15. 其它
`conda`（conda run -n 环境）、`python3/sys.executable`、`taskset`、`gffread`、`multiqc`、POSIX shell 工具

### A16. biosoft/ 中的参考工具（⚠️ 非运行时依赖）
以下目录**仅作为算法参考/移植来源**，管线代码中均为注释性提及（"复刻/借鉴/移植自"）或已纯 Python 重实现，**没有任何 subprocess 调用**；目录已 gitignore，不进 pixi、不进部署清单：

| 目录 | 在代码中的角色 |
|---|---|
| `VirPhyKit` | geo_analysis.py 复刻其 RRT/TempMig/GeoSubsampler 逻辑 |
| `YR-MPE` | multigene_tools.py / saturation_analysis.py 参考其拼接与 C 值法口径 |
| `EasyHap` | virus_vcf_pipeline.py / virus_haplotype_qst.py 的 LD 热图与单倍型命名思路 |
| `shinyTempSignal` | clock_analysis.py 稳健 RTT 层借鉴其设计（纯 Python 重实现） |
| `vfam_trees` | acvirus_tree_pro.py 的定根/LCA/配色代码移植来源 |
| `viralclust` | 无任何引用 |

---

## B. Python 第三方包

| 包 | 用途 | 出现脚本 |
|---|---|---|
| `Bio/biopython` | 序列/比对/Entrez | 几乎所有脚本 |
| `pandas` | 表格处理 | 大量 |
| `polars` | 高性能表格 | virome_pipeline, cluster_pipeline, filter_virus, batch_virus_depth* |
| `numpy` / `scipy` | 数值/统计 | 大量（phylo, analysis, 绘图） |
| `matplotlib` / `seaborn` | 绘图 | 大量 |
| `plotly` + `kaleido` | 交互图/导出 | report_pipeline, taxonomic_sankey |
| `pysam` | BAM 处理 | batch_virus_depth*, batch_virus_variants |
| `scikit-learn` | PCA/KMeans | virus_variants_analyzer, snpgenie_master |
| `tqdm` / `psutil` | 进度/资源 | 大量 |
| `pyhmmer` | HMMER 蛋白域搜索 | virus_identification, biosoft/pyhmmer_* |
| `torch` + `transformers` | DNABERT-2 病毒鉴定 | viralm_cpu.py, envs/viralm.yaml |
| `venn` / `upsetplot` | 韦恩/UpSet 图 | virus_identification.py |
| `yaml`(PyYAML) | 配置解析 | virome_pipeline, submission_pipeline |
| `pymol` | 蛋白结构渲染 | batch_draw_pymol.py |
| `networkx` / `pypopart` | 图/单倍型网络 | phylo/utils/build_haplo_outputs.py |
| `openai` | AI 报告 | pixi 声明 |
| `requests` / `urllib3` / `xmltodict` | HTTP/XML | pixi 声明 |
| `pyside6` | 桌面 GUI (metadata_gui / submission_gui / virome_submission_pipeline) | GUI 全部入口 |
| `scikit-learn` | PCA/KMeans | virus_variants_analyzer, snpgenie_master, virus_vcf_pipeline |
| `openpyxl` | 读 xlsx (元数据/VMR 表) | gsa_sra.search/info, data_store, submission_gui |
| `treetime` | RTT 定年门禁 (Python 包, bioconda) | phylo_pipeline.py |
| `cartopy` + `geopy` | 地理坐标/地图可视化 | phylo utils/sample_plot.py, virphy_bridge.py |
| `marsilea` | 复合图 | phylo utils/pub_plots.py |
| `pyarrow` | Arrow 数据 | envs/viralm.yaml |
| huggingface_hub / safetensors / gdown / accelerate / peft / datasets / tokenizers | ViraLM 模型 | envs/viralm.yaml |

---

## C. R 包

| 包 | 用途 | 出现脚本 |
|---|---|---|
| `optparse` | CLI 解析 | 所有 *.R |
| `data.table` | 快速表格 | virus_classifier_analysis.R, panvirome_suite.R |
| `ggplot2` | 绘图 | 所有 *.R |
| `dplyr/tidyr` | 数据整理 | virema_summary_report.R, panvirome_suite.R |
| `stringr` | 字符串 | virus_classifier_analysis.R |
| `VennDiagram/ggVennDiagram` | 韦恩图 | virus_classifier_analysis.R |
| `grid/gridExtra` | 布局 | virus_classifier_analysis.R, panvirome_oncoprint.R |
| `RColorBrewer` / `scales` / `viridis` | 配色/刻度 | virus_classifier_analysis.R, virus_frequency_plot.R |
| `cowplot/patchwork` | 拼版 | virus_classifier_analysis.R, panvirome_suite.R |
| `ComplexHeatmap` | 热图/oncoprint | panvirome_oncoprint.R |
| `ggrepel/ggpubr/zoo/fs` | 标注/统计 | panvirome_suite.R |
| `maftools`(Bioconductor) | 突变注释可视化 | panvirome_suite.R |
| `circlize/httr` | 环形图/HTTP | virema_summary_report.R |
| `ape/beautier/ggbeast/bdskytools/XML/seqinr` | 系统发育/定年(历史) | phylo/archive/*.R |

---

## D. 数据库 / 参考文件

### D1. 鉴定数据库
| 数据库 | 默认路径 |
|---|---|
| geNomad DB | `~/database/virus-db/genomad_db/` |
| RVDB v31(diamond) | `~/database/virus-db/RVDB-v31/RVDB_viroids.diamond_db/U-RVDBv31.0-prot_unique.dmnd` |
| RVDB v31(hmm/fasta/annot) | `U-RVDBv31.0-prot.hmm` / `.EX.acc.fasta` / `.info.tab` |
| UniRef90 | `~/database/uniport_db/uniref90/uniref90.dmnd` |
| NCBI nr | `~/database/nr_db/nr.dmnd` |
| VirSorter2 DB | `~/database/virus-db/virsorter2_db` |
| ViralVerify HMM | `~/database/virus-db/viralverify_db/nbc_hmms.hmm` |
| Metabuli DB | `~/database/virus-db/RVDB-v31/RVDB_viroids.metabuli_db` |
| Viroids DB | `~/database/virus-db/viroids-db/viroids.fasta` |
| MMseqs2 DB | `~/database/virus-db/RVDB-30/RVDB.mmseqs` |
| ncbi-virus_ref | `~/database/virus-db/ncbi-virus_ref/ncbi-virus_ref.blast.db` |

### D2. 宿主/组装/拯救
| 数据库 | 默认路径 |
|---|---|
| host_db | `~/database/host_db/` |
| Kraken2 宿主库 | `~/database/kraken2/k2_pluspfp_20260626`（标准预建索引, 下载见 DATABASE_SETUP §6） |
| CheckV DB | `~/database/virus-db/checkv-db-v1.7/` |
| genus_lens | `~/database/virus-db/db/genus_lens`（+ `genus_lens_no_prefix.json`；构建脚本 make_genus-length.py 已收编 utils/db_build/） |
| 植物病毒参考 | `~/plant_virus_db/3.final-ref-virus.db/final.complete_ref.fasta` 等 |

### D3. 分类/宿主预测
| 数据库 | 默认路径 |
|---|---|
| NCBI taxonomy | `~/database/taxonomy/` |
| ICTV MSL41 | `~/database/virus-db/VMR_MSL41.v1.20260320.xlsx` + `database/ICTV_MSL41_*.tsv` |
| ACVirus DB | `~/database/virus-db/acvirus_db` |
| vConTACT3 / VITAP DB | vcontact3 `--db-path` / VITAP `-d` |
| PhaBOX2 DB | `~/database/virus-db/phabox_db_v2_2/` |
| ICTV 宿主概率表 | `database/cross_analysis/` |
| suvtk DB | `~/database/virus-db/suvtk_db/`（`suvtk download-database`, 备用 Zenodo 10.5281/zenodo.15374439, 见 DATABASE_SETUP §3.8） |

### D4. 深度学习/结构域
| 数据库 | 默认路径 |
|---|---|
| ViraLM 模型 | `~/database/virus-db/viralm_db/`(DNABERT-2, gdown) |
| VirHunter 权重 | `biosoft/virhunter/weights/generalistic` |
| VirBot 参考 | `biosoft/VirBot/ref/VirBot.dmnd` |
| RdRpCatch DB | `~/database/virus-db/rdrp-db/rdrpcatch_dbs/`（`rdrpcatch databases` 命令下载, 见 DATABASE_SETUP §2.9） |
| CDD 库 | `~/database/cdd/cdd-db/cdd_db` + `cdd_classified_taxid.tsv` + `cdd_virus_final_v4.txt`(MMPV-RNA/database/cdd)；**构建脚本 `~/database/cdd/build_cdd_databases.sh` 已收编 utils/db_build/**；两个产物文件被 `database/*` gitignore，本地工作副本缺失需从服务器拷回 |
| pyhmmer HMM | `{db_dir}/hmm/viral/combined.hmm`、`hmm/pfam/Pfam-A-*.hmm` |

---

## 部署口径（三类来源）

1. **pixi 可装**（`pixi.toml`，conda-forge+bioconda）：约 120 包，含全部 Python/R 库 + 主流生信工具（2026-09-06 逐包 API 核查过渠道存在性）。
2. **手动拷贝/源码编译**（`biosoft/`+`utils/`）：ViReMa.py(biosoft/ViReMa)、VirBot.py、VirHunter/predict_cpu.py、snpEff.jar/SnpSift.jar、viralm_cpu.py、genome_rmDuplicates.pl、penguin(自编译)、refineC(编译)、flye(biosoft 源码版)、PhaBOX2、CHERRY。
3. **独立 conda 环境**：`viralm`（envs/viralm.yaml，Py3.8+Torch2.0.1+CUDA11.8）、`rdrpcatch`、`virhunter`。

> ⚠️ **virome_phylo_pipeline 不在 pixi 覆盖范围**：BEAST/treeannotator/logcombiner、RDP5、PAML(codeml)、DnaSP、TempMig、Spread3、SeqHarvester/SeqGrouper、pypopart、treedater、drhip、shinyTempSignal 等需按 A12 单独部署（Python 侧的 treetime/cartopy/geopy/marsilea/networkx 已入 pixi）。

## 路径配置：环境变量与占位符（2026-09-22 引入）

**设计目标：换环境不改代码、尽量不改配置文件。** 三层解耦：

| 层 | 机制 | 位置 |
|---|---|---|
| 代码层 | 占位符 `{repo}`（仓库根）/ `{base}`（数据根）/ `{results}`（产物根）/ `{dir}`（条目自身目录）；环境变量 `${VAR}` / `${VAR:-默认}`（支持嵌套，迭代展开） | `virome_phylo_pipeline/utils/dataset_config.py`、`virome_discovery_pipeline/virome_pipeline.py: _expand_env_config()` |
| 配置层 | 单一配置点 + 环境变量兜底；profile 用 YAML 锚点继承（只写差异） | `pipeline_config.yaml`、`virome_phylo_pipeline/datasets.yaml` |
| 脚本层 | 兜底路径一律按 `Path(__file__)` / 仓库相对定位，**禁止写死机器绝对路径** | 各脚本 |

### 环境变量清单

**主流程（`virome_pipeline.py` / `pipeline_config.yaml`）**

| 变量 | 作用 | 默认 |
|---|---|---|
| `MMPV_DB_ROOT` | 数据库根 | `/home/zhangwenda/database` |
| `MMPV_VIRUS_DB` | 病毒库根 | `{MMPV_DB_ROOT}/virus-db` |
| `MMPV_HOST_DB` | 宿主库根 | `{MMPV_DB_ROOT}/host_db/` |
| `MMPV_KRAKEN2_DB` | Kraken2 宿主库 | `{MMPV_DB_ROOT}/kraken2/k2_pluspfp_20260626` |
| `MMPV_CHECKV_DB` | CheckV 数据库 | `{MMPV_DB_ROOT}/virus-db/checkv-db-v1.7` |
| `MMPV_PLANT_VIRUS_DB` | 植物病毒参考库根 | `/home/zhangwenda/plant_virus_db` |
| `MMPV_SALMON` / `MMPV_DIAMOND` / `MMPV_RAGTAG` / `MMPV_VIRBOT` | 工具可执行文件 | 见配置文件 |

**系统发育（`virome_phylo_pipeline/datasets.yaml`）**

| 变量 | 作用 | 默认 |
|---|---|---|
| `MMPV_DATA_BASE` | 数据根 | `/home/zhangwenda/MMPV-paper/results_PSTVd_GCVA_20260905` |
| `MMPV_RESULTS` | 运行产物根 | `{repo}/virome_phylo_pipeline/phylo_results` |
| `MMPV_PYPOPART` | pypopart src | `{repo}/biosoft/pypopart/src` |
| `MMPV_BEAST` / `MMPV_BEAST_CP` / `MMPV_BEAST_JLP` | BEAST 可执行与 classpath | 见配置文件 |
| `MMPV_ANNOT_DIR` | 注释 GB 目录 | 空（缺失则相关步骤跳过） |
| `MMPV_VMR` / `MMPV_PLANT_VIRUS_DB` | VMR 表 / 植物病毒库 | 见 `annotate_nucleic_acid.py` 查找链 |

**提交管线（`virome_submission_pipeline/`）**

| 变量 | 作用 | 默认 |
|---|---|---|
| `MMPV_SUBMIT_BASE` | 提交数据根 | `/home/zhangwenda/virus/data-2026/data-test` |
| `MMPV_SUBMIT_MIX` | 混合样品目录 | `/home/zhangwenda/virus/ningxiagouqi/11.merge_assembly/mix/out` |
| `MMPV_SUBMIT_CORE13` | Core13 元数据表 | `.../global_metadata/Global_Unified_Metadata_Core13.tsv` |
| `MMPV_T2ASN` | table2asn 可执行 | `/home/zhangwenda/.pixi/bin/table2asn` |

### 换环境的三种做法（按侵入性递增）

```bash
# ① 只设环境变量（推荐，零改动）
export MMPV_DB_ROOT=/data/db MMPV_DATA_BASE=/data/results
python virome_pipeline.py --profile plant --output_dir out/

# ② 改配置里的 base / results 两行（datasets.yaml 顶部）
# ③ 复制配置文件并 --config my_config.yaml（最灵活）
```

### 本次修复的历史问题

- `datasets.yaml` 原称"改 base 即可"，实则需改 6 处（base + env×4 + PSTVD_FULL 绝对路径），且 `haplo.*` 不复用 `dir` 而写绝对长路径 → 已改为 `${VAR:-默认}` + 占位符复用
- `pipeline_config.yaml` 的 `plant` profile 曾**完整复制** default 且**漏了 10 个键**（`deplete.*`、`tools.ragtag/virbot`、`taxonomy.*`、`cobra.*` 等，运行时会静默用 argparse 默认值）→ 已改为锚点继承（48 行 → 15 行，补齐全部）
- 2 处 `D:\桌面\...` Windows 路径写死（`orphan_audit.py`、`dependency_audit.py`）→ 改为 `Path(__file__)` 自解析
- 6 处兜底绝对路径（`build_haplo_outputs.py` 的 pypopart、`annotate_nucleic_acid.py` 的 plant_virus_db、submission 的 4 处）→ 环境变量 + 仓库相对

> ⚠️ **YAML 浅合并陷阱**：`<<: *anchor` 只合并**顶层键**，嵌套字典会被整体替换。给 `plant` 加锚点继承时若只写 `cluster: {min_length: 300}`，default 的 `cluster.ani/qcov` 会**丢失**（实测踩到）。正确做法：嵌套段各自加锚点（如 `cluster: &cluster_default` + `<<: *cluster_default`）。

## 易踩坑
- **vcontact3 是 opt-in 工具**：`virus_classifier.py` 里 `optin_tools=["vcontact3"]`，`--tools all` 只含 7 工具（genomad/metabuli/CAT/diamond_lca/VITAP/mmseqs/ACVirus）；vConTACT3_db 数据库虽装（14G）但流程默认不调用，需显式 `--tools vcontact3` 才运行。
- **sambamba 全程没用到**（统一 `samtools`）；**COBRA 阶段实际是 `bwa-mem2 + samtools + coverm`**，cobra 只是目录/阶段名。
- **SDT 无独立可执行**，`sdt_genus_matrix.py` 内部调 `mafft` 两两比对（支持 muscle/clustalw）。
- **RDP5 是 Windows GUI**，管线用 `rdp5_ini.py`(逆向 ini 生成器)+`rdp5.sh` 在 Linux 跑。
- **genus_lens 有两种格式**：原始 `genus_lens`(g__Genus) 与 `genus_lens_no_prefix.json`。
- **SnpEff 需 JRE≥21**（pixi 已改为 java-openjdk>=21）。
- **bioconda 无包名单（2026-09-06 API 逐一核查）**：penguin、refineC、flye(conda版非 biosoft 2.9)、PhaBOX2、CHERRY、VirHunter、VirBot、ViraLM —— 这些在 pixi.toml 里只能进注释/手动区，加回 `[dependencies]` 会导致 `pixi install` 整体求解失败。
- **kaleido 的 conda-forge 包名是 `python-kaleido`**（直接写 kaleido 会 404）；`datavzrd` 在 conda-forge 而非 bioconda。
- **biosoft/ 中 VirPhyKit/YR-MPE/EasyHap/shinyTempSignal/vfam_trees/viralclust 是参考工具**（见 A16），不要写进任何部署清单。
- **Virseqimprover 硬依赖 bedtools**：VSI scaffold 截取直接执行 `bedtools getfasta`，rescue 分支 B 的机器必须有 bedtools。