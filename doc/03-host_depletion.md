# host_depletion.py — 去宿主 + 去 rRNA

## 流程

```
Kraken2 分类 → Bowtie2/HISAT2 去宿主 → rRNA 去除 (Ribodetector 或 SILVA)
```

## 用法

```bash
# Ribodetector (默认, 仅 rna-short)
python host_depletion.py \
    --tool bowtie2 --seq-type rna-short \
    --kraken2_index <dir> --step2_index <dir> \
    --input-dir <dir> --outdir <dir> \
    --rrna --jobs 10 --threads 40

# SILVA Bowtie2 (所有 seq_type)
python host_depletion.py \
    --tool bowtie2 --seq-type rna-short \
    --kraken2_index <dir> --step2_index <dir> \
    --input-dir <dir> --outdir <dir> \
    --rrna --rrna_tool silva --silva_index <prefix> \
    --jobs 10 --threads 40
```

## 参数

| 参数 | 默认 | 说明 |
|------|------|------|
| `--tool` | 必需 | 比对工具: bowtie2/hisat2/minimap2 |
| `--seq-type` | 必需 | 测序类型: dna-short/rna-short/nanopore/pacbio |
| `--kraken2_index` | 必需 | Kraken2 宿主库 |
| `--step2_index` | 必需 | 比对索引前缀 |
| `--input-dir` / `--input` | 必需 | 输入目录/文件 |
| `--outdir` / `--output` | 必需 | 输出目录 |
| `--jobs` | 1 | 并行样本数 |
| `--threads` | 4 | 单样本线程 |
| `--rrna` | 关闭 | 开启 rRNA 剔除 |
| `--rrna_tool` | ribodetector | rRNA 工具: ribodetector/silva |
| `--silva_index` | — | SILVA Bowtie2 索引 (silva 模式必需) |
| `--keep_rrna` | 关闭 | 保留分离出的 rRNA reads |
| `--chunk_size` | 256 | Ribodetector chunk_size |
| `--confidence` | 0.4 | Kraken2 置信度 |
| `--filter` | true | true=去除宿主, false=提取宿主 |
| `--steps` | kraken2,align,rrna | 执行步骤 (消融实验) |
| `--force` | — | 强制重跑 |

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

## 断点续传

- 自动维护 `.checkpoints` 隐藏文件
- `--force` 清除 checkpoint
- 粒度: 每个样本独立追踪
