# GenBank 提交工作流 / GenBank Submission Workflow

> 基于 suvtk 工具链 | MMPV-RNA → GenBank .sqn 提交文件

---

## 一、概述 / Overview

suvtk 是一个专门用于向 GenBank 提交病毒序列的命令行工具，自动化 ICTV 分类、ORF 预测、功能注释、MIUVIG 元数据整合和 .sqn 生成。

| 步骤 | 命令 | 输入 | 输出 |
|------|------|------|------|
| 1 | `suvtk taxonomy` | FASTA 序列 | `taxonomy.tsv` + `miuvig_taxonomy.tsv` |
| 2 | `suvtk features` | FASTA + taxonomy.tsv | `.tbl` 特征表 + `.fna` + `.faa` |
| 3 | 手动准备 | — | `source.src` + `miuvig.tsv` + `assembly.tsv` |
| 4 | `suvtk comments` | 步骤 1/2/3 输出 | `output.cmt` |
| 5 | `suvtk table2asn` | 步骤 2/3/4 输出 | `submission.sqn` ✅ |

---

## 二、安装与数据库 / Installation & Database

```bash
# 安装 suvtk
pip install suvtk

# 下载数据库 (~5GB)
aria2c -x 16 -s 16 https://zenodo.org/records/15423947/files/suvtk_db_v0.1.1.tar.gz
tar -xzf suvtk_db_v0.1.1.tar.gz -C ~/database/virus-db/

# 或使用 suvtk 内置下载
suvtk download-database -o ~/database/virus-db/suvtk_db/
```

---

## 三、自动化提交脚本 / Automated Submission Script

`virome_submission_pipeline/virome_submission.py` 将上述 5 步集成为一个命令。

### 3.1 新病毒 (Novel Viruses)

```bash
# 一键运行 (步骤1-4)
python virome_submission_pipeline/virome_submission.py novel \
    --fasta $OUT/08_Rescue/Plant/centroids/final_centroids.fasta \
    --taxonomy $OUT/05_Taxonomy/integrated/final_integrated_classification.tsv \
    --host $OUT/06_HostPrediction/ensemble_host_summary.tsv \
    --checkv $OUT/08_Rescue/checkv/ \
    --suvtk-db ~/database/virus-db/suvtk_db/ \
    --output ./genbank_submission/novel/ \
    -t 40

# 编辑 source.src 后, 分步运行 table2asn
python virome_submission_pipeline/virome_submission.py step \
    --step table2asn \
    --work-dir ./genbank_submission/novel/ \
    --suvtk-db ~/database/virus-db/suvtk_db/
```

### 3.2 已知病毒 (Known Viruses)

```bash
# 一键运行
python virome_submission_pipeline/virome_submission.py known \
    --fasta $OUT/known_viruses/3_Virus_assemblies_final/ \
    --summary $OUT/known_viruses/1_FastViromeExplorer/summary/best.summary.tsv \
    --ref-info /db/ref_info.tsv \
    --suvtk-db ~/database/virus-db/suvtk_db/ \
    --output ./genbank_submission/known/ \
    -t 40
```

### 3.3 分步模式 (断点续传)

```bash
python virome_submission_pipeline/virome_submission.py step --step taxonomy    --fasta seqs.fasta --work-dir ./work/ --suvtk-db ~/db/ -t 40
python virome_submission_pipeline/virome_submission.py step --step features    --fasta seqs.fasta --work-dir ./work/ --suvtk-db ~/db/ -t 40
python virome_submission_pipeline/virome_submission.py step --step metadata    --work-dir ./work/
python virome_submission_pipeline/virome_submission.py step --step comments    --work-dir ./work/ --suvtk-db ~/db/
python virome_submission_pipeline/virome_submission.py step --step table2asn   --work-dir ./work/ --suvtk-db ~/db/
```

---

## 四、手动逐步执行 / Manual Step-by-Step

### Step 1: suvtk taxonomy

```bash
suvtk taxonomy \
    -i sequences.fasta \
    -o 1_taxonomy/ \
    -d ~/database/virus-db/suvtk_db/ \
    -s 0.7 \
    -t 40
```

**输出文件:**
- `taxonomy.tsv`: 主分类表 (`contig`, `taxonomy`, `taxid`)
- `miuvig_taxonomy.tsv`: MIUVIG 格式 (`contig`, `pred_genome_type`, `pred_genome_struc`)
- `taxonomy.log`: MMseqs2 LCA 比对日志

**关键参数:**
- `-s 0.7`: MMseqs2 灵敏度 (0-1, 越高越敏感但越慢)
- `-t 40`: 线程数

---

### Step 2: suvtk features

```bash
suvtk features \
    -i sequences.fasta \
    -o 2_features/ \
    -d ~/database/virus-db/suvtk_db/ \
    --coding-complete \
    --taxonomy 1_taxonomy/taxonomy.tsv \
    -t 40
```

**输出文件:**
- `*.tbl`: NCBI 5-column feature table (CDS 位置 + 产物 + 推断证据)
- `reoriented_nucleotide_sequences.fna`: 可能被重新定向的序列
- `proteins.faa`: 预测的所有蛋白序列
- `miuvig_features.tsv`: 软件/数据库版本信息

**处理逻辑:**
1. 序列定向 (如 Negarnaviricota 负链 → 反转)
2. pyrodigal-gv ORF 预测 (筛选 coding-complete: CDS > 50% 基因组)
3. MMseqs2 对比 BFVD → 功能注释 (如 "RNA-directed RNA polymerase")
4. 无命中的 → "hypothetical protein"
5. 格式化为 .tbl

---

### Step 2.5: 假定蛋白分析 (可选)

```bash
# 对 hypothetical proteins 进一步分析
python ~/bin/analyze_hypothetical.py \
    -t 2_features/sequences.tbl \
    -f 2_features/proteins.faa \
    -o hypothetical.faa \
    --blast -r blast_results.txt
```

---

### Step 3: 准备元数据文件 / Prepare Metadata

#### 3.1 source.src (每条序列的样本信息)

```tsv
Sequence_ID	Organism	Isolate	Collection_date	geo_loc_name	Lat_Lon	Bioproject	Biosample	SRA	Metagenomic	Metagenome_source	Segment
contig_001	Chrysoviridae sp.	SAMPLE001_isolate	15-Jun-2024	China:Jiangsu	32.06 N 118.79 E	PRJNA123456	SAMN12345678	SRR12345678	TRUE	soil metagenome
contig_002	Partitiviridae sp.	SAMPLE001_isolate	15-Jun-2024	China:Jiangsu	32.06 N 118.79 E	PRJNA123456	SAMN12345678	SRR12345678	TRUE	soil metagenome
```

**关键字段说明:**

| 字段 | 必须 | 说明 |
|------|------|------|
| `Sequence_ID` | ✅ | 与 taxonomy.tsv contig 列一致 |
| `Organism` | ✅ | 来自 taxonomy.tsv taxonomy 列 |
| `Isolate` | ✅ | 唯一标识符; 分段病毒同一病毒必须相同 |
| `Collection_date` | ✅ | 格式 DD-Mmm-YYYY (如 15-Jun-2024) |
| `geo_loc_name` | ✅ | 格式 Country:Region (如 China:Jiangsu) |
| `Lat_Lon` | ✅ | 格式 XX.XX N/S XXX.XX E/W |
| `Bioproject` | — | PRJNA 登录号 |
| `Biosample` | — | SAMN 登录号 |
| `SRA` | — | SRR 登录号 |
| `Metagenomic` | — | 必须为 TRUE |
| `Metagenome_source` | — | 如 "soil metagenome" |
| `Segment` | — | 分段病毒的片段编号 |

#### 3.2 miuvig.tsv (全局 MIUVIG 参数)

```tsv
sample_id	viral_enrichment	sequencing_platform	sequencing_method	assembly_software	assembly_method	quality_check_software
ALL	rRNA_depletion	Illumina_NovaSeq	RNA-SHORT	MEGAHIT	metatranscriptomic	CheckV
```

#### 3.3 assembly.tsv (组装信息)

```tsv
Sequencing_Technology	Assembly_Method	Assembly_Name	Assembly_Software	Coverage
Illumina_NovaSeq	metatranscriptomic_assembly	MMPV-RNA_v2.3	MEGAHIT	NOT_PROVIDED
```

#### 3.4 template.sbt (NCBI 作者模板)

从 NCBI 网站生成: https://submit.ncbi.nlm.nih.gov/genbank/template/submission/

---

### Step 4: suvtk comments

```bash
suvtk comments \
    --taxonomy 1_taxonomy/miuvig_taxonomy.tsv \
    --features 2_features/miuvig_features.tsv \
    --miuvig 3_metadata/miuvig.tsv \
    --assembly 3_metadata/assembly.tsv \
    --quality checkv/completeness.tsv \
    -o 4_comments/
```

**输出:** `output.cmt` — 结构化 MIUVIG 注释文件

---

### Step 5: suvtk table2asn → .sqn

```bash
suvtk table2asn \
    --fasta 2_features/reoriented_nucleotide_sequences.fna \
    --features 2_features/reoriented_nucleotide_sequences.tbl \
    --source 3_metadata/source.src \
    --comments 4_comments/output.cmt \
    --template template.sbt \
    -o 5_submission/
```

**输出:**
- `submission.sqn`: 可直接上传 GenBank 的 Sequin 文件 ✅
- `submission.val`: 验证报告 (错误/警告/信息)

---

## 五、辅助工具 / Utility Tools

### 5.1 co-occurrence (分段病毒关联)

```bash
# 已知病毒: 通过丰度相关性找出分段病毒的共现模式
suvtk co-occurrence \
    --abundance best.summary.tsv \
    -o cooccurrence_out/
```

### 5.2 gbk2tbl (格式转换)

```bash
# 将 phold 等工具生成的 .gbk 转为 .tbl
suvtk gbk2tbl -i input.gbk -o output_dir/
```

### 5.3 virus-info (分段病毒提示)

```bash
# 检查哪些序列可能属于已知的分段病毒科
suvtk virus-info -t taxonomy.tsv
```

---

## 六、提交前检查清单 / Pre-Submission Checklist

- [ ] `taxonomy.tsv` 中每条序列都有非 NA 的 taxonomy
- [ ] `.tbl` 特征表中 CDS 位置正确 (无内部终止密码子)
- [ ] `source.src` 中所有占位符已替换为真实值
- [ ] 分段病毒的所有片段使用相同的 `Isolate` 值
- [ ] `metagenomic` 列均为 TRUE
- [ ] `collection_date` 格式正确 (DD-Mmm-YYYY)
- [ ] `geo_loc_name` 格式正确 (Country:Region)
- [ ] template.sbt 已从 NCBI 下载
- [ ] `submission.val` 验证报告无 ERROR (WARNING 可接受)
- [ ] BioProject / BioSample / SRA 登录号已在 NCBI 注册

---

## 七、本次实际数据提交方案 / Actual Submission Plan

### 7.1 管道产出概览 / Pipeline Output Summary

| 类别 | 数量 | 来源 |
|------|------|------|
| 总 centroids | 2,062 | `04_CLUSTER/centroids/final_centroids.fasta` |
| ★ known (Species≠NA) | 783 | taxonomy |
| ★★ novel_species (Genus≠NA) | 84 | taxonomy |
| ★★ novel_genus (Family≠NA) | 41 | taxonomy |
| ★★★ truly_novel (all NA) | 43 | taxonomy |
| Plant-related viruses | ~60 | Virgaviridae, Bromoviridae, Partitiviridae, etc. |
| CD-HIT known linked | 3 | `04_CLUSTER/2_cdhit/known_linked_centroids.fasta` |
| CheckV pass (≥90%) | 6 | `07_Checkv/checkv_pass_ids.txt` |
| 已知病毒检测 | — | `known_viruses/1_FastViromeExplorer/summary/best.summary.tsv` |
| 已知病毒变异 | — | `known_viruses/2_Virus_variants_Results/` |

### 7.2 实际运行命令 / Actual Commands

```bash
# 设定变量
OUT=/home/zhangwenda/data-test2/out
SUVTK_DB=~/database/virus-db/suvtk_db

# ═══ 新病毒提交 (novel) — 使用全部 centroids ═══
python virome_submission_pipeline/virome_submission.py novel \
    --fasta $OUT/04_CLUSTER/centroids/final_centroids.fasta \
    --taxonomy $OUT/05_Taxonomy/integrated/final_integrated_classification.tsv \
    --host $OUT/06_HostPrediction/ensemble_host_summary.tsv \
    --suvtk-db $SUVTK_DB \
    --output ./genbank_submission/novel/ \
    -t 40

# 如果只需要植物相关病毒子集 (推荐用于植物病毒组):
# 先从 taxonomy 提取植物病毒科的 contig ID, 然后子集化 FASTA
python -c "
import polars as pl
# 已知植物病毒科
plant_families = ['Partitiviridae','Chrysoviridae','Totiviridae','Megatotiviridae',
    'Orthototiviridae','Endornaviridae','Amalgaviridae','Botourmiaviridae',
    'Bromoviridae','Virgaviridae','Tombusviridae','Sobemoviridae',
    'Potyviridae','Closteroviridae','Luteoviridae','Nanoviridae',
    'Geminiviridae','Tymoviridae','Secoviridae','Alphaflexiviridae',
    'Betaflexiviridae','Gammaflexiviridae','Mitoviridae','Narnaviridae',
    'Ourmiavirus','Fimoviridae','Phenuiviridae','Tospoviridae',
    'Rhabdoviridae','Ophioviridae','Aspiviridae','Benyviridae',
    'Furovirus','Pecluvirus','Pomovirus','Tobravirus','Hordeivirus',
    'Mayoviridae','Kitaviridae','Tepovirus','Carlavirus','Foveavirus',
    'Capillovirus','Trichovirus','Vitivirus','Ampelovirus']
tax = pl.read_csv('$OUT/05_Taxonomy/integrated/final_integrated_classification.tsv', separator='\t')
plant_ids = tax.filter(pl.col('Family').is_in(plant_families))
print(f'Plant virus contigs: {plant_ids.height}')
plant_ids.select('contig_id').write_csv('plant_virus_ids.txt', include_header=False)
"

# 然后用 seqkit 提取子集
# seqkit grep -f plant_virus_ids.txt $OUT/04_CLUSTER/centroids/final_centroids.fasta > plant_novel.fasta
# python virome_submission_pipeline/virome_submission.py novel --fasta plant_novel.fasta ...

# ═══ 已知病毒提交 (known) — 需要先运行 full 阶段 ═══
# (当前全长组装未完成, 需先执行):
# python scripts/auto_known_virus.py --stage full \
#     --reads_dir $OUT/00b_HostDepletion/ \
#     --output_dir $OUT/known_viruses/ \
#     --ref_info /db/ref_info.tsv --reference /db/ref.fasta \
#     -t 40 -j 4

# 完成后运行:
# python virome_submission_pipeline/virome_submission.py known \
#     --fasta $OUT/known_viruses/3_Virus_assemblies_final/ \
#     --summary $OUT/known_viruses/1_FastViromeExplorer/summary/best.summary.tsv \
#     --suvtk-db $SUVTK_DB \
#     --output ./genbank_submission/known/ \
#     -t 40
```

### 7.3 植物病毒科参考 / Plant Virus Families Reference

| Category | Families |
|----------|----------|
| +ssRNA | Virgaviridae, Bromoviridae, Tombusviridae, Potyviridae, Closteroviridae, Luteoviridae, Tymoviridae, Secoviridae, Alphaflexiviridae, Betaflexiviridae, Gammaflexiviridae, Benyviridae, Sobemoviridae, Mayoviridae, Kitaviridae |
| -ssRNA | Rhabdoviridae, Phenuiviridae, Tospoviridae, Fimoviridae, Ophioviridae, Aspiviridae |
| dsRNA | Partitiviridae, Chrysoviridae, Totiviridae, Endornaviridae, Amalgaviridae, Megatotiviridae, Orthototiviridae |
| ssDNA | Geminiviridae, Nanoviridae |
| Unclassified | Botourmiaviridae, Mitoviridae, Narnaviridae (fungal/plant-associated) |

---

## 八、常见问题 / FAQ

**Q: taxonomy.tsv 中大量序列标记为 "Viruses;unclassified" 怎么办？**
A: 这些序列可能是真正新颖的病毒 (novel family+)。在 source.src 中将 Organism 设置为 "Viruses;unclassified RNA virus" 或类似的临时分类。

**Q: .tbl 特征表中的 CDS 有内部终止密码子？**
A: 检查序列是否被正确定向。如果 taxonomy 预测了负链 RNA 病毒但序列是正向的，features 步骤会自动反转。

**Q: table2asn 验证报告中有 ERROR？**
A: 常见原因: source.src 格式错误 (缺列/多列)、Organism 名称不符合 ICTV 规范、日期格式错误。逐条检查 .val 文件。

**Q: 如何获得 BioProject / BioSample 登录号？**
A: 在 NCBI Submission Portal (https://submit.ncbi.nlm.nih.gov/) 先注册 BioProject 和 BioSample，获得登录号后再填写到 source.src。

---

*参考: https://landerdc.github.io/suvtk/index.html | suvtk v0.1.1*
