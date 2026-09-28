# 数据库下载与配置指南

> 全部数据库默认存放路径: `~/database/virus-db/`

---

## 1. 同源比对依赖数据库

### 1.1 NCBI Virus RefSeq (蛋白)

```bash
# 下载病毒蛋白参考序列
wget https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.1.protein.faa.gz

# 建立 Diamond 索引 (含分类信息)
diamond makedb --in viral.1.protein.faa.gz \
    -d virus.pep.dmnd \
    --taxonmap ~/database/taxonomy/prot.accession2taxid \
    --taxonnodes ~/database/taxonomy/nodes.dmp \
    --threads 60
```

### 1.2 NCBI NR (蛋白)

```bash
# 下载 NCBI nr 数据库
wget https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz

# 建立 Diamond 索引
diamond makedb --in nr.gz \
    -d nr.dmnd \
    --taxonmap ~/database/taxonomy/prot.accession2taxid \
    --taxonnodes ~/database/taxonomy/nodes.dmp \
    --threads 60
```

> **与 §1.4 UniRef90 二选一（默认 UniRef90）**：管线入口分别是 pipeline_config.yaml 的
> `nr_db` 与 `uniprot_db`，默认走 `uniref90.dmnd`；NR 缺失时识别段自动降级（跳过 NR 抢救验证）。
> 仅当需要 NR 级全覆盖去假阳性时才部署 nr（~250GB）。

### 1.3 ClusteredNR (加速替代)

```
来源: http://bioinfo.bti.cornell.edu/ftp/program/VirusDetect/virus_database/v248/U100/
ClusteredNR 将 nr 蛋白聚类: 簇内 ≥90% 相同且长度 ≥ 最长成员 90%
每个簇挑选代表序列, 加速 BLAST 搜索
```

### 1.4 UniProt / UniRef90

```bash
# UniRef90 (diamond 格式)
wget https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/uniref90.fasta.gz
diamond makedb --in uniref90.fasta.gz -d uniref90.dmnd --threads 60
```

### 1.5 RVDB (参考病毒数据库)

```bash
# RVDB - Reference Viral DataBase
# https://fzer.github.io/rvdbtools/
wget https://github.com/fzer/rvdbtools/raw/main/db/RVDB_v31.fasta.gz

# 包含 viroid 子数据库
# 建立 BLAST 索引
makeblastdb -in RVDB_v31.fasta -dbtype nucl -out RVDB_v31.blast
```

### 1.6 Viroids 数据库

```bash
# 类病毒参考序列 (来自 VirusDetect)
wget http://bioinfo.bti.cornell.edu/ftp/program/VirusDetect/virus_database/v248/U100/viroids.fasta
makeblastdb -in viroids.fasta -dbtype nucl -out viroids.blast
```

### 1.7 ICTV/NCBI 完整植物病毒基因组

```bash
# 参考用途: CD-HIT 预聚类 + BLASTN 抢救
# 从 NCBI virus 下载完整植物病毒基因组
# 路径: ~/database/virus-db/ncbi-virus_ref/

# EVE Stage1 正向发现库 (同目录下蛋白版, endogenous_virus_pipeline)
wget https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.1.protein.faa.gz
diamond makedb --in viral.1.protein.faa -d ~/database/virus-db/ncbi-virus_ref/ncbi-virus_ref.pep.dmnd --threads 60
```

### 1.8 植物/病毒拆分库 (hybrid) + id2div — EVE 双侧判定专用

```bash
# 用途: endogenous_virus_pipeline Stage2 双侧判定
#   把 SwissProt 拆成 viral + plant 两个子集合进同一个 diamond 库,
#   id2div_final.tsv 记录每条 sseqid 归属 (viral|plant),
#   判定逻辑: viral 侧 bitscore ≥50 且 ≥ plant 侧 → viral_supported
# 路径: ~/database/uniport_db/hybrid/ (见 pipeline_config.yaml databases.eve_pv_db)

# 1) 下载 SwissProt (flat 文件含 OX/OC 分类行; fasta 用于抽子集)
wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz
wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
gunzip uniprot_sprot.dat.gz uniprot_sprot.fasta.gz

# 2) 按分类拆分 ( Viruses / Viridiplantae ) 并生成 id2div
python - <<'EOF'
from pathlib import Path

viral, plant, div = [], [], []
acc = None
lineage_viral = lineage_plant = False

def flush():
    if not acc:
        return
    for line, org in (("viral", lineage_viral), ("plant", lineage_plant)):
        if org:
            (viral if line == "viral" else plant).append(acc)
            div.append(f"{acc.split()[0]}\t{line}")
            break

with open("uniprot_sprot.dat", encoding="utf-8") as fh:
    for line in fh:
        if line.startswith("AC   "):
            acc = line[5:].rstrip().rstrip(";").split(";")[0].strip()
        elif line.startswith("OC   "):
            terms = line[5:].split(";")
            if any(t.strip().startswith("Viruses") for t in terms):
                lineage_viral = True
            if any(t.strip().startswith("Viridiplantae") for t in terms):
                lineage_plant = True
        elif line.startswith("//"):
            flush()
            acc, lineage_viral, lineage_plant = None, False, False

Path("id2div_final.tsv").write_text("\n".join(div) + "\n", encoding="utf-8")
print(f"viral={len(viral)} plant={len(plant)} -> id2div_final.tsv")
EOF

# 3) 拆分 fasta 并合并建库 (seqkit 按 accession 列表抽子集)
seqkit grep -f <(awk -F'\t' '$2=="viral"{print $1}' id2div_final.tsv) uniprot_sprot.fasta -o viral.fasta
seqkit grep -f <(awk -F'\t' '$2=="plant"{print $1}' id2div_final.tsv) uniprot_sprot.fasta -o plant.fasta
cat viral.fasta plant.fasta > plant_virus_final.fasta
diamond makedb --in plant_virus_final.fasta -d plant_virus_final.dmnd -p 16

mkdir -p ~/database/uniport_db/hybrid
mv plant_virus_final.dmnd id2div_final.tsv ~/database/uniport_db/hybrid/
```

> 说明: 2026 年 hi-fever 试点用的是已建好的 `plant_virus_final.dmnd`（SwissProt 全量拆
> viral + Viridiplantae，viral/plant 条目比例约 1:2）。若服务器上已有现成库，确认
> id2div 两列格式为 `sseqid<TAB>viral|plant` 即可直接用；`eve_screen.py` 启动时会
> 校验四个库文件存在，缺失即报错退出，不会静默跑错。

---

## 2. 同源比对非依赖数据库

### 2.1 geNomad

```bash
# 安装
pixi global install -c conda-forge -c bioconda genomad

# 下载数据库 (~14GB)
genomad download-database ~/database/virus-db/genomad_db

# 来源: https://zenodo.org/records/14886553
```

### 2.2 Cenote-Taker 3

```bash
# 安装
mamba install cenote-taker3=3.4.3

# 下载全部数据库
get_ct3_dbs -o ~/database/virus-db/ct3_DBs \
    --hmm T --hallmark_tax T --refseq_tax T \
    --mmseqs_cdd T --domain_list T \
    --hhCDD T --hhPFAM T --hhPDB T
```

### 2.3 VirSorter2

```bash
# 安装 + 数据库初始化
virsorter config --init-source --db-dir ~/database/virus-db/virsorter2_db/
```

### 2.4 ViralVerify

```bash
# https://github.com/ablab/viralVerify
# 需要 HMM 模型数据库
# 默认路径: ~/database/virus-db/viralverify_db/
```

### 2.5 CheckV

```bash
# https://bitbucket.org/berkeleylab/checkv/
checkv download_database ~/database/virus-db/checkv-db-v1.7/
```

### 2.6 ViraLM (DNABERT-2 模型)

```bash
# https://github.com/ChengPENG-wolf/ViraLM
git clone https://github.com/ChengPENG-wolf/ViraLM.git
conda env create -f envs/viralm.yaml -n viralm

# 下载模型 (~1.5GB)
gdown --id 1EQVPmFbpLGrBLU0xCtZBpwvXrtrRxic1
tar -xzvf model.tar.gz -C ~/database/virus-db/viralm_db/
```

### 2.7 VirHunter (深度学习权重)

```bash
# https://github.com/cbib/virhunter
# 权重下载: https://www.dropbox.com/scl/fi/vuln5dgqpfh5n73quya1r/
# 已整合到 biosoft/virhunter/weights/generalistic/
```

### 2.8 VirBot

```bash
# https://github.com/GreyGuoweiChen/VirBot
git clone https://github.com/GreyGuoweiChen/VirBot.git
# 参考数据集从 OneDrive 下载
```

### 2.9 RdRpCatch 数据库 (Stage 2 鉴定, RdRp 检测)

```bash
# https://github.com/dimitris-karapliafis/RdRpCATCH
# 新版命令为 `rdrpcatch databases` (旧版 `rdrpcatch download`), 支持追加自定义 pHMM 库
conda run -n rdrpcatch rdrpcatch databases -o ~/database/virus-db/rdrp-db/
# 产出目录即管线引用的 rdrpcatch_dbs:
#   ~/database/virus-db/rdrp-db/rdrpcatch_dbs/
# 完整离线包 (含数据库+脚本) 也可从官网下载: https://rdrpcatch.bioinformatics.nl/
```

---

## 3. 分类数据库

### 3.1 NCBI Taxonomy

```bash
# 下载 NCBI 分类数据库
wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
tar -xzf taxdump.tar.gz -C ~/database/taxonomy/

# 提取病毒 taxid (10239 = Viruses)
taxonkit list --ids 10239 --indent "" > viral_taxIDs.txt
```

### 3.2 MMseqs2 分类数据库

```bash
# RVDB 转 MMseqs2 格式
mmseqs createdb RVDB_v31.fasta RVDB.mmseqs
mmseqs createtaxdb RVDB.mmseqs tmp --tax-mapping-file taxon.map

# ICTV MMseqs2 蛋白数据库
# https://github.com/apcamargo/ictv-mmseqs2-protein-database
```

### 3.3 Accession2Taxid 映射

```bash
# 核酸 accession → taxid
wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz

# 蛋白 accession → taxid
wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz
```

### 3.4 MEGAN 映射文件

```bash
# MEGAN NCBI-nr 分类映射
wget https://software-ab.cs.uni-tuebingen.de/download/megan6/megan-map-Feb2022.db.zip
```

### 3.5 VITAP 数据库

```bash
# VMR-MSL40 分类数据库
# https://github.com/xishenyuzhou/VITAP
```

### 3.6 vConTACT3 参考数据库

```bash
# https://bitbucket.org/MAVERICLab/vcontact3/
```

### 3.7 ACVirus 数据库

```bash
# https://github.com/icelu/ACVirus
# 构建 (VMR xlsx 在仓库 database/ 下):
ACVirus create_db --vmr VMR_MSL41.v1.20260320.xlsx --dbpath ~/database/virus-db/acvirus_db
```

### 3.8 suvtk 数据库 (提交管线, taxonomy/features 注释)

```bash
# https://github.com/LanderDC/suvtk  (Submission of Uncultivated Viral genomes toolkit)
suvtk download-database -o ~/database/virus-db/suvtk_db/
# 备用: 手动从 Zenodo 下载解压 (https://doi.org/10.5281/zenodo.15374439)
# 管线调用: submission_pipeline.py --suvtk-db ~/database/virus-db/suvtk_db/
```

---

## 4. 宿主预测数据库

### 4.1 PhaBOX2

```bash
# https://github.com/KennthShang/PhaBOX
# 数据库路径: ~/database/virus-db/phabox_db_v2_2/
```

### 4.2 ICTV 宿主概率表

```bash
# 由 classify_contigs.py 使用的 cross_analysis/ 目录
```

### 4.3 NCBI Host 信息

```bash
# host.dmp 文件 (NCBI taxonomy 中提取)
```

---

## 5. 物种基因组长度参考

```bash
# NCBI Assembly Reports - 物种基因组大小
wget https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/species_genome_size.txt

# 提取病毒属的平均长度
# (脚本已收编至 virome_discovery_pipeline/utils/db_build/make_genus-length.py, 服务器源 ~/bin/)
python virome_discovery_pipeline/utils/db_build/make_genus-length.py --taxid 10239 \
    -o virus_genus_lens_stats.tsv \
    --genus-lens-output ~/database/virus-db/genus_lens
```

---

## 6. 宿主去除数据库 (本管线构建)

```bash
# 由 public_metadata_pipeline/build_host_pipeline.py 自动构建
python public_metadata_pipeline/build_host_pipeline.py \
    --species "<species_name>" --taxid <taxid> \
    --stage all --threads 30
```

产出:
```
hostdb/
├── kraken2/       # Kraken2 分类库
├── bowtie2/       # Bowtie2 比对索引
├── hisat2/        # HISAT2 比对索引
└── minimap2/      # Minimap2 比对索引
```

**[可选]** 标准 Kraken2 预建通用索引 (host_depletion.py 去宿主第一步, config 引用 `k2_pluspfp_20260226`)。默认流程是自建宿主库 (见 §6 开头的 build_host_pipeline.py, 有基因组数据即可建 kraken2/bowtie2/hisat2/minimap2 四索引); 无参考基因组或需通用去宿主时才用下面的大库:

```bash
# Ben Langmead 预建索引: https://benlangmead.github.io/aws-indexes/k2
mkdir -p ~/database/kraken2/k2_pluspfp_20260226
wget -c https://genome-idx.s3.amazonaws.com/kraken/k2_pluspfp_20260226.tar.gz \
    -O - | tar -xz -C ~/database/kraken2/k2_pluspfp_20260226
# 注: 新版体积远超早期 ~50GB(实测 0226 tar 已 ~180GB); 建议用 aria2c 多线程下载
# 若该日期版本已下架, 从索引页取最新 k2_pluspfp_YYYYMMDD 并同步改 pipeline_config.yaml
```

---

## 7. 数据库版本清单

| 数据库 | 版本 | 大小 | 用途 |
|--------|------|------|------|
| NCBI nr | latest | ~250GB | Diamond BLASTX 去假阳性 |
| RVDB | v31 | ~5GB | 病毒参考序列 |
| geNomad DB | v1.7 | ~14GB | 病毒鉴定 |
| CheckV DB | v1.7 | ~3GB | 完整性评估 |
| VirSorter2 DB | latest | ~2GB | 病毒分类 |
| ViraLM model | DNABERT-2 | ~1.5GB | DL 病毒鉴定 |
| VirHunter weights | generalistic | ~6MB | DL 病毒鉴定 |
| NCBI taxonomy | latest | ~200MB | 分类信息 |
| MEGAN map | Feb2022 | ~20GB | 分类注释 |
| Kraken2 host | custom | ~50GB | 宿主去除 |
| Bowtie2/HISAT2 host | custom | ~3GB | 宿主比对 |
| PhaBOX2 | v2.2 | ~2GB | 宿主预测 |

---

## 8. pixi global install 工具安装清单 (2026-08-26 更新)

> 全部走 `conda-forge + bioconda`；**flye 用 biosoft 源码版（2.9-b1779，已装，不走 pixi）；penguin / RNAVirHost 是 conda 官方版（mambaforge）**（详见 SOFTWARE_VERSIONS.txt）。

```bash
# A. 质控 / 去宿主
pixi global install -c conda-forge -c bioconda fastp seqkit clumpify bbmap seqtk pigz
pixi global install -c conda-forge -c bioconda kraken2 bowtie2 hisat2 minimap2 ribodetector
pixi global install -c conda-forge -c bioconda hocort hostil ribodetector fastp

# B. 组装
pixi global install -c conda-forge -c bioconda megahit spades rnaviralspades

# C. 病毒鉴定
pixi global install -c conda-forge -c bioconda diamond ncbi-blast+ genomad metabuli virsorter2

# D. 聚类 / 分类
pixi global install -c conda-forge -c bioconda cd-hit mmseqs2 taxonkit vclust coverm

# E. 比对 / 质量 / 建树 / 变异
pixi global install -c conda-forge -c bioconda samtools bwa-mem2 checkv
pixi global install -c conda-forge -c bioconda mafft iqtree
pixi global install -c conda-forge -c bioconda bcftools bgzip tabix ivar freebayes lofreq snpsift

# F. 数据获取
pixi global install -c conda-forge -c bioconda sra-tools aria2c ncbi-datasets-cli entrez-direct

# G. Python 包 (pixi 全局)
pixi global install -c conda-forge numpy scipy pandas polars biopython matplotlib seaborn plotly kaleido pysam scikit-learn tqdm psutil pyhmmer pyyaml requests openpyxl
```

> 不走 pixi 的：`CAT / VITAP / ACVirus / vConTACT3`(git 仓库, ACVirus 在 biosoft/ACVirus)、`viralverify / VirBot / VirHunter / ViraLM / RdRpCatch`(独立仓库或 conda env viralm/rdrpcatch/virhunter)。

## 9. 数据库核对状态 (2026-08-26)

| 依赖 | 路径 | 状态 |
|---|---|---|
| CDD 主库 | `~/database/cdd/cdd-db/cdd_db` (mmseqs 全套) | ✅ |
| CDD 病毒域白名单 | `~/MMPV-RNA/database/cdd/cdd_virus_final_v4.txt` (278KB) | ✅ (~/database/cdd 同款) |
| CDD 分类 taxid | `~/MMPV-RNA/database/cdd/cdd_classified_taxid.tsv` = **v2 (47.4MB, 08-06)** | ✅ 已统一（~/database/cdd 也已换 v2，旧 v1 备份为 cdd_classified_taxid_v1.tsv.bak） |
| CheckV DB | `~/database/virus-db/checkv-db-v1.7/` (6.5G, genome_db+hmm_db) | ✅ config 用 v1.7 |
| genus_lens (病毒属长度) | `~/database/virus-db/db/genus_lens` (63KB) + genus_lens_no_prefix.json | ✅ |
| ACVirus DB | `~/database/virus-db/acvirus_db/` (taxa.txt 2.8MB + all_virus.fasta 586MB) | ✅ 构建命令 `ACVirus create_db --vmr VMR_MSL40...xlsx --dbpath acvirus_db` |
| UniRef90 | `~/database/uniport_db/uniref90/uniref90.dmnd` | ✅ Diamond 索引 |
| ncbi-virus 参考 | `~/database/virus-db/ncbi-virus/` (nt 48GB) | ✅ |
| 建树/SDT/拯救脚本 | `~/bin/acvirus_tree_pro.py`、`virome_analysis_pipeline/sdt_genus_matrix.py`、`virome_discovery_pipeline/Virseqimprover.py` | ✅ |
| vConTACT3 DB | `~/database/virus-db/vConTACT3_db/` (14G, v232) | ✅ 库在，但**vcontact3 为 opt-in 工具，默认 `--tools all` 不包含、流程默认不跑**（需显式 `--tools vcontact3`） |
| VITAP DB | `~/database/virus-db/vitap-db/` (2.6G, VMR-MSL40) | ✅ 库在；VITAP 在默认 7 工具内，正常参与 05 分类 |
| CDD 构建脚本 | `~/database/cdd/build_cdd_databases.sh` (生成 cdd_db + 白名单 + 分类 taxid 表) | 📦 收编至 `virome_discovery_pipeline/utils/db_build/` |

---

## 10. 自写建库脚本收编 (utils/db_build/)

以下自写脚本此前只在服务器（`~/bin`、`~/database/cdd/`），仓库无版本化。已建收编目录
`virome_discovery_pipeline/utils/db_build/`（含 README 与拷贝命令），
`build_virus_db.py` 会自动把该目录插到 PATH 最前，拷入即用：

| 脚本 | 服务器源 | 用途 |
|---|---|---|
| `SearchAccessionIdToTaxId.py` | `~/bin/` | build_virus_db.py accession→taxid |
| `db_seqid2taxid_add_legth.py` | `~/bin/` | build_virus_db.py seqid2taxid 补长度 |
| `make_ktaxonomy.py` | `~/bin/` | build_virus_db.py ktaxonomy 映射 |
| `make_genus-length.py` | `~/bin/` | genus_lens 生成（§5） |
| `build_cdd_databases.sh` | `~/database/cdd/` | CDD 主库 + `cdd_virus_final_v4.txt` + `cdd_classified_taxid.tsv`（filter_virus.py / gen_final_judgement.py 硬依赖） |

> ⚠️ `database/cdd/` 两个产物文件被 `.gitignore` 的 `database/*` 忽略且本仓库工作副本中没有，
> 需从服务器拷回（命令见 db_build/README.md），否则本地跑 CDD 过滤/裁决阶段会直接失败。
