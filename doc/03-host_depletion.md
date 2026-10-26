# host_depletion.py — 去宿主 + 去rRNA / Host Depletion & rRNA Removal（`data_preprocessing_pipeline/`）

> Kraken2 / Kun-peng + Bowtie2/HISAT2/Minimap2 + Ribodetector/SILVA. Ablation experiments via `--steps`, resource monitoring per sample.

## 流程 / Workflow

```
分类初筛 (Kraken2 或 Kun-peng) → Bowtie2/HISAT2 去宿主 → rRNA 去除 (Ribodetector 或 SILVA)
```

两种分类后端输出完全一致（per-read 5 列 + kreport2）。Kun-peng 不写分流 FASTQ，
脚本会解析 per-read 结果后用 `seqkit grep` 按 read ID 反选。

## 用法

```bash
# Ribodetector (默认, 仅 rna-short)
python data_preprocessing_pipeline/host_depletion.py \
    --tool bowtie2 --seq-type rna-short \
    --kraken2_index <dir> --step2_index <dir> \
    --input-dir <dir> --outdir <dir> \
    --rrna --jobs 10 --threads 40

# SILVA Bowtie2 (所有 seq_type)
python data_preprocessing_pipeline/host_depletion.py \
    --tool bowtie2 --seq-type rna-short \
    --kraken2_index <dir> --step2_index <dir> \
    --input-dir <dir> --outdir <dir> \
    --rrna --rrna_tool silva --silva_index <prefix> \
    --jobs 10 --threads 40

# Kun-peng 后端 (同一个宿主库; 首次自动 hashshard 转换)
python data_preprocessing_pipeline/host_depletion.py \
    --classifier kunpeng --kunpeng-bin /opt/bin/kun_peng \
    --tool bowtie2 --seq-type rna-short \
    --kraken2_index <kraken2_or_kunpeng_db> --step2_index <dir> \
    --input-dir <dir> --outdir <dir> \
    --rrna --jobs 10 --threads 40
```

## 参数

| 参数 | 默认 | 说明 |
|------|------|------|
| `--tool` | 必需 | 比对工具: bowtie2/hisat2/minimap2 |
| `--seq-type` | 必需 | 测序类型: dna-short/rna-short/nanopore/pacbio |
| `--kraken2_index` | 必需 | 宿主分类库 (Kraken2 `hash.k2d` 或 Kun-peng `hash_*.k2d`) |
| `--step2_index` | 必需 | 比对索引前缀 |
| `--classifier` | kraken2 | 初筛分类器: kraken2 / kunpeng |
| `--kunpeng-bin` | `kun_peng` | Kun-peng 可执行文件 (或 `MMPV_KUNPENG_BIN`) |
| `--kunpeng-hash-capacity` | 1G | 自动 hashshard 转换时的分片容量 |
| `--min-hit-groups` | 2 | Kun-peng 最小 hit group 数 (kraken2 忽略) |
| `--keep-classifier-output` | 关闭 | 保留 Kun-peng per-read 输出/chunk (体积大) |
| `--input-dir` / `--input` | 必需 | 输入目录/文件 |
| `--outdir` / `--output` | 必需 | 输出目录 |
| `--jobs` | 1 | 并行样本数 |
| `--threads` | 4 | 单样本线程 |
| `--rrna` | 关闭 | 开启 rRNA 剔除 |
| `--rrna_tool` | ribodetector | rRNA 工具: ribodetector/silva |
| `--silva_index` | — | SILVA Bowtie2 索引 (silva 模式必需) |
| `--keep_rrna` | 关闭 | 保留分离出的 rRNA reads |
| `--chunk_size` | 256 | Ribodetector chunk_size |
| `--confidence` | 0.4 | 分类置信度 (两种后端通用) |
| `--keep-taxids` | 10239 | 从 classified 池救回的 taxID (默认病毒) |
| `--filter` | true | true=去除宿主, false=提取宿主 |
| `--steps` | kraken2,align,rrna | 执行步骤 (消融实验) |
| `--force` | — | 强制重跑 |

## Kun-peng 后端说明

- **库复用**：`kun_peng hashshard --db <kraken2库>` 原地追加 `hash_*.k2d` +
  `hash_config.k2d`（保留 `hash.k2d`），所以同一个宿主库可同时被
  `kraken2` 与 `kun_peng` 使用。`--classifier kunpeng` 首次运行会自动转换。
  ⚠ 分片是 hash 表的一份副本，首次转换会使库目录额外占用约等于原 `hash.k2d`
  的磁盘（k2_pluspfp 大库尤其明显）；`--kunpeng-hash-capacity` 只控制单片大小。
- **read 分流**：Kun-peng 的 `classify` 只产出 `output_*.txt`（Kraken 格式）
  与 `*.kreport2`，不提供 `--classified-out/--unclassified-out`。脚本解析
  per-read 结果得到宿主 ID，再用 `seqkit grep -v -f` 反选；`--keep-taxids`
  的救回逻辑不变（taxid 在保留集合中的 reads 不算宿主）。
  ⚠ `seqkit grep` 需要把宿主 ID 集合装进内存（约 80 B/ID）；数千万宿主 reads
  的单样本会额外占数 GB，并行跑多样本时请相应调低 `--jobs`。
- **双端**：`kun_peng classify -P R1 R2` 每个 read 对输出一行（ID 取 R1 去
  `/1` 后的 base ID）；脚本按输入 FASTQ 真实后缀（MGI `/1` `/2` vs
  Illumina 无后缀）写 seqkit 模式文件，保证两端都能精确命中。
- **对称性**：`--min-hit-groups` 默认 2（上游默认），比 kraken2 更严格；
  需要与 kraken2 完全对齐时可设 `--min-hit-groups 1`。

## SILVA 数据库准备

```bash
# 下载 SILVA
wget https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/SILVA_138.1_LSUParc_SSUParc_tax_silva.fasta.gz
gunzip SILVA_138.1_LSUParc_SSUParc_tax_silva.fasta.gz

# 建 Bowtie2 索引
bowtie2-build SILVA_138.1_LSUParc_SSUParc_tax_silva.fasta silva_index
```

## 输出

```
{outdir}/
├── {sample}/
│   ├── {sample}_clean_1.fastq.gz    清洁 R1
│   └── {sample}_clean_2.fastq.gz    清洁 R2
├── rrna/                            分离出的 rRNA (--keep_rrna)
├── logs/
│   ├── kraken2/
│   ├── bowtie2/
│   ├── ribodetector/ 或 silva/
│   └── seqkit/
└── summary.tsv                      各步骤统计
```

> 路径归属：`{outdir}` 就是去宿主阶段自己的输出目录 —— standard 布局下为
> `<项目>/02_Preprocessing/02_HostDepleted/`（默认），legacy 为
> `<共享根>/00b_HostDepletion/`。`rrna/` 恒为 `{outdir}` 的直接子目录
> （`--keep_rrna` 才产生，默认丢弃）；文件名按工具区分：
> ribodetector → `{sample}.rrna.[1|2].fq.gz`，silva → `{sample}.silva_rrna.[1|2].fq.gz`。
> rRNA 统计报告由 `--rrna_report` 指定（默认 `ribodetector.report.txt`），落在 `--logs_dir`。

## 断点续传

- 自动维护 `.checkpoints` 隐藏文件
- `--force` 清除 checkpoint
- 粒度: 每个样本独立追踪

## 宿主库版本与 EVE 掩码库

去宿主结果只取决于**用的是哪个宿主库**，所以库版本必须记录在样本表里（当前 `summary.tsv` 不自动写库版本，需人工登记或读 `db_provenance.json`）。

| 库目录 | 怎么来 | 什么时候用 |
|---|---|---|
| `<work>/hostdb/` | `build_host_pipeline.py --stage hostdb`（默认路径，行为与 v3.1 一致）| 常规；与已发表结果对齐 |
| `<work>/hostdb_evemask/` | `build_host_pipeline.py --stage eve-mask,hostdb --eve-verdict ...` | 敏感性分析 / 确认 EVE 造成假阴时；见 `doc/16-host-eve-masking.md` §6 |

两个库的目录布局相同（`kraken2/`、`bowtie2/host.*`、`hisat2/host.*`、`minimap2/host_<type>.mmi`），可直接互换。

⚠ **衔接坑**：`data_preprocessing.py` 的 `find_host_db()` 只会在 `--host_db` 给的目录里找 `hostdb/`（或目录本身含 `hash.k2d`），**不会自动发现同级的 `hostdb_evemask/`**。要掩码库必须显式写：

```bash
python data_preprocessing_pipeline/data_preprocessing.py \
    --host_db <work>/hostdb_evemask \
    --input <样本表> --output 02_Preprocessing --stage deplete \
    --tool bowtie2 --seq-type rna-short --threads 40
```

（`host_depletion.py` 自身的 `-k/--kraken2_index` + `-x/--step2_index` 则是逐目录显式传入，不存在这个问题。）

**掩码的代价**（为什么要单独当一个臂）：掩码区的全部宿主 reads 都会变成「非宿主」进入组装/鉴定，假阳性、嵌合与低覆盖 contig 会上升；且 `--keep-taxids 10239` 的病毒救援发生在 Kraken2 之后、比对之前，**掩码/未掩码都不影响 Stage 4 的「比对到宿主即删除」逻辑**——掩码是唯一能让 EVE 区外源病毒 reads 不被当作宿主删掉的手段。验收标准见 `doc/16-host-eve-masking.md` §5 的 A/B 表。

⚠ **版本前提**：上面这句以本仓 `HEAD` 版 `host_depletion.py` 为准。服务器上未入库的 7 月版在 `--keep-taxids` 非空时（默认 `10239`）**无条件跳过 align**，此时根本没有 Stage 4；跑掩码臂前先按 `doc/17-server-vs-local-audit-20261022.md` §3 确认样本是用哪一版跑的。
