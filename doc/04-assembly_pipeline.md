# assembly_pipeline.py — 宏转录组组装

## 流程

```
选择工具: penguin (默认) / megahit / rnaviralspades / all (三种并行)
```

## 用法

```bash
python assembly_pipeline.py --tool all --input <dir> --output-dir <dir> [参数]
```

## 参数

| 参数 | 默认 | 说明 |
|------|------|------|
| `--tool` | penguin | penguin / megahit / rnaviralspades / all |
| `--input` | 必需 | 输入 reads 目录 |
| `--output-dir` | 必需 | 输出目录 |
| `--threads` / `-n` | 4 | 线程数 |
| `--memory` / `-m` | 32 | 内存 GB |
| `--jobs` / `-j` | 1 | 并行样本数 |
| `--log_dirs` | — | 日志目录 |
| `--force` | — | 强制重跑 |

## 输出

```
{output}/
├── {sample}/
│   ├── {sample}_penguin.contig.fasta
│   ├── {sample}_megahit.contig.fasta
│   └── {sample}_rnaviralspades.contig.fasta
└── logs/
```

## 断点续传

- 检查最终 contig FASTA 是否存在
- `--force` 跳过检查, 重新组装
