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
```

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
python ~/bin/make_genus-length.py --taxid 10239 \
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
