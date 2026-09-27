# GenBank 提交准备（virome_submission_pipeline）—— 论文级方法章

> 本方法章描述病毒序列的 GenBank 提交准备流程，可直接粘贴进投稿 Methods。工具名、参数取自 `virome_submission_pipeline/` 脚本实际调用。软件版本号标 ⊙ 者待确认。文末附参考文献与占位符清单。

---

## 1. 序列提交准备流程

对 MMPV-RNA 鉴定的新病毒与已知病毒序列，用 suvtk v0.1.1（GitHub: LanderDC/suvtk，无正式发表文献）工具链准备 GenBank 提交文件（.sqn）。流程分五步，由 `virome_submission.py` 编排并支持断点续传。新病毒序列来自发现管道的 rescue 阶段（final_centroids.fasta），已知病毒序列来自分析管道的全长组装阶段。

## 2. 分类学分配（suvtk taxonomy）

用 MMseqs2（Steinegger & Söding, 2017）对每条序列做 LCA 检索，比对 ICTV 参考蛋白库，分配病毒分类学（门到种各阶元）并预测基因组类型与结构（`-s 0.7` 灵敏度）。输出 taxonomy.tsv 与 miuvig_taxonomy.tsv。无法分配至已知分类阶元的序列按实际阶元标记为 novel。

## 3. 特征表构建（suvtk features）

先做序列方向校正（负链 RNA 病毒反转至编码方向），再用 pyrodigal-gv（Larralde, 2022）预测开放阅读框，筛选 coding-complete 序列（CDS 占基因组 > 50%）。ORF 产物经 MMseqs2 比对 BFVD 病毒蛋白结构库（Kim et al., 2025）做功能注释，无命中者标注 hypothetical protein，可经 BLASTP 在线检索补充注释。输出五列特征表 .tbl（CDS 位置、产物、推断证据）与 proteins.faa。

## 4. 序列标识规范化

GenBank 要求序列 ID 为 NCBI 本地 ID 格式。对非标准 ID（如组装器产生的含点号长 ID）统一加 `lcl|` 前缀（仅保留字母数字下划线），并同步规范化 FASTA、.tbl、source.src、comments.cmt 四个文件中的 ID，避免 table2asn 报 "Malformatted ID" 解析失败。

## 5. 元数据准备

按 GenBank 与 MIUVIG 标准准备三类元数据：

- **source.src**（逐序列样本信息）：Sequence_ID、Organism、Isolate（分段病毒各片段共用同一 isolate）、Collection_date（DD-Mmm-YYYY）、geo_loc_name（Country:Region）、Lat_Lon、BioProject、BioSample、SRA 登录号、Metagenomic（TRUE）、Metagenome_source。样本字段可从 `gsa_sra.info.py` 产出的 Core14 元数据表自动填充。
- **miuvig.tsv**（全局 MIUVIG 参数）：source_uvig（如 viral fraction RNA metagenome）、assembly_software（MEGAHIT 1.2.9）、assembly_method、sequencing_platform、viral_enrichment。
- **assembly.tsv**（组装信息）与 **template.sbt**（从 NCBI 提交门户下载的作者模板）。

## 6. 结构化注释整合与 .sqn 生成（suvtk comments + table2asn）

`suvtk comments` 将 miuvig_taxonomy、特征表、miuvig.tsv 与 assembly.tsv 整合为结构化注释文件 output.cmt；`suvtk table2asn` 将重定向后的核苷酸序列、特征表、source.src、comments 与 template.sbt 打包生成 Sequin 提交文件 submission.sqn，并输出验证报告 submission.val。提交前核对：taxonomy 无非 NA 记录、CDS 无内部终止密码子、占位符已替换、分段病毒 isolate 一致、metagenomic 均为 TRUE、验证报告无 ERROR。

## 7. 软件与版本

| 工具 | 版本 | 引用 |
|---|---|---|
| suvtk | 0.1.1 | GitHub（无正式文献） |
| MMseqs2 | bd01c22 | Steinegger & Söding, 2017 |
| pyrodigal-gv | ⊙ | Larralde, 2022 |
| BFVD | ⊙ | Kim et al., 2025 |
| table2asn | ⊙ | NCBI（无文献） |
| MEGAHIT | 1.2.9 | Li et al., 2015 |

分析在 Linux 计算服务器（246）执行，Python 3.10；suvtk 数据库（约 5 GB）经 `suvtk download-database` 部署于 `~/database/virus-db/suvtk_db/`。

---

## 参考文献

- Kim RS, Levy Karin E, Mirdita M, Chikhi R, Steinegger M. 2025. BFVD—a large repository of predicted viral protein structures. Nucleic Acids Research 53:D340–D347.
- Larralde M. 2022. Pyrodigal: Python bindings and interface to Prodigal, an efficient method for gene prediction in prokaryotes. Journal of Open Source Software 7:4296.
- Li D, Liu CM, Luo R, Sadakane K, Lam TW. 2015. MEGAHIT: an ultra-fast single-node solution for large and complex metagenomics assembly via succinct de Bruijn graph. Bioinformatics 31:1674–1676.
- Steinegger M, Söding J. 2017. MMseqs2 enables sensitive protein sequence searching for the analysis of massive data sets. Nature Biotechnology 35:1026–1028.

---

## 占位符清单（投稿前需替换）

- [ ] 提交的序列数与病毒种类数
- [ ] BioProject / BioSample / SRA 登录号
- [ ] 标 ⊙ 的软件版本号
- [ ] 若改用其他组装软件或自建 suvtk 数据库，更新 miuvig.tsv / assembly.tsv 相应字段
- [ ] GenBank 登录号区间（或 "will be deposited"）
