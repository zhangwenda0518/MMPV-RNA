# clean-data.py — 数据清洗 / Data Cleaning（`data_preprocessing_pipeline/`）

> Fastp QC → Seqkit FASTQ-to-FASTA → Clumpify clustering/reordering for better compression. Checkpoint-resume with `.clean_checkpoints`.

## 流程 / Workflow

```
Fastp 质控 → Seqkit 统计 → Clumpify 聚类重排（提升压缩率）
```

## 用法

```bash
python data_preprocessing_pipeline/clean-data.py --input <dir> --output <dir> [参数]
```

## 参数

| 参数 | 默认 | 说明 |
|------|------|------|
| `--input` | 必需 | 输入 FASTQ 目录 |
| `--output` | 必需 | 输出目录 |
| `--fastp-threads` | 4 | Fastp 线程数 |
| `--jobs` | 1 | 并行任务数 |
| `--skip-clumpify` | — | 跳过 Clumpify 聚类重排（提升压缩率） |
| `--force` | — | 强制重跑 (清除 checkpoint) |
| `--dedup` | — | fastp 自带去重 (默认未开启) |
| `--clumpify-memory` | 10g | Clumpify 内存 |

## 输出

```
{output}/
├── 1.fastp/          Fastp 质控报告 (HTML/JSON)
├── 2.fasta/          清洗后 FASTA
├── 3.clumpify/       Clumpify 重排 reads
└── logs/
```

## 断点续传

- 自动维护 `.clean_checkpoints` 隐藏文件
- `--force` 清除 checkpoint 重新运行
- 粒度: 每个样本独立追踪
