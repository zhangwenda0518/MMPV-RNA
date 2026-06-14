# virus_identification16.py — 6工具病毒鉴定 / 6-Tool Virus Identification

> Genomad + Diamond BLASTX + VirSorter2 + ViralVerify + VirHunter + Metabuli. Post-filtering (UniProt/NR), Venn/Upset visualization.

## 6 个鉴定工具 / 6 Identification Tools

| 工具 | 原理 |
|------|------|
| Genomad | 深度学习病毒/质粒/前病毒分类 |
| Diamond BLASTX | RefSeq 病毒蛋白 + NR 库比对 |
| VirSorter2 | 隐马尔可夫模型病毒检测 |
| ViralVerify | 病毒蛋白验证 |
| VirHunter | 机器学习病毒鉴定 |
| Metabuli | 基于 k-mer 的分类 |

## 用法

```bash
python virus_identification16.py --input <fasta> --output <dir> --db_dir <dir> [参数]
```

## 参数

| 参数 | 默认 | 说明 |
|------|------|------|
| `--input` / `-i` | 必需 | 输入 contig FASTA (文件或目录) |
| `--output` / `-o` | 必需 | 输出目录 |
| `--db_dir` | 必需 | 病毒鉴定数据库根目录 |
| `--identify_tools` | all | 逗号分隔: genomad,diamond,virsorter2,... |
| `--threads` / `-t` | 4 | 线程数 |
| `--jobs` / `-j` | 1 | 并行样本数 |
| `--force` / `-f` | — | 强制重跑 |

## 输出

```
{output}/{sample}/
├── {sample}_genomad_taxonomy.tsv
├── {sample}_diamond_blastx.tsv
├── {sample}_virsorter2_result.tsv
├── {sample}_viralverify_result.tsv
├── {sample}_virhunter_result.tsv
├── {sample}_metabuli_result.tsv
├── {sample}_virus.all.candidate.fasta  候选病毒序列
└── {sample}.processing.done            完成标记
```

## 断点续传

- 检查 `.processing.done` 完成标记
- 各工具独立检查结果文件
- `--force` 重跑全部工具
