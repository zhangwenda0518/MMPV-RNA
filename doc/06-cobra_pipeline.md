# cobra_pipeline.py — COBRA 批量延伸 / COBRA Batch Extension

> BWA-MEM2 → CoverM coverage → COBRA overlap extension. Auto-matches reads+contig+virus triples, JSON checkpoint.

## 流程 / Workflow

```
自动匹配 reads + contig + virus 三元组
BWA-MEM2 比对 → CoverM 覆盖度 → COBRA 重叠延伸 → CheckV 评估
```

## 用法

```bash
python cobra_pipeline.py --mode mix \
    --reads-dir <dir> --contigs-dir <dir> --virsorter-dir <dir> \
    --output-dir <dir> [参数]
```

## 参数

| 参数 | 默认 | 说明 |
|------|------|------|
| `--mode` | mix | mix / single |
| `--reads-dir` | 必需 | 清洁 reads 目录 |
| `--contigs-dir` | 必需 | 组装 contig 目录 |
| `--virsorter-dir` | 必需 | 病毒鉴定输出目录 |
| `--output-dir` | 必需 | 输出目录 |
| `--assembly-tools` | — | 逗号分隔: megahit,rnaviralspades,penguin |
| `--jobs` | 1 | 并行任务数 |
| `--threads` | 4 | 线程数 |
| `--resume` | — | 断点续传 |

## 输出

```
{output}/{sample}/
├── cobra_{tool}_result/
│   ├── *.cobra.fa              延伸后序列
│   └── checkv/quality_summary.tsv
└── logs/
```

## 断点续传

- JSON checkpoint 文件 (`checkpoint_status.json`)
- `--resume` 跳过已完成任务
- 粒度: 每个 {sample}.{mode}.{tool} 任务独立追踪
