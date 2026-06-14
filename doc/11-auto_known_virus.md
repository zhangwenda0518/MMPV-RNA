# auto_known_virus.py — 已知病毒分析主控 / Known Virus Orchestrator

> 3-stage analysis: detect (rapid quantification) → variants (SnpEff+SnpGenie) → full (de novo assembly). Independent of the main discovery pipeline.

## 三个独立阶段 / Three Independent Stages

```
--stage detect    batch_virus_depth40.py   Salmon/Bowtie2 快速检测
--stage variants  batch_virus_variants.py  FreeBayes + SnpEff + SnpGenie
--stage full      batch_virus_full.py      virus-full 全长组装
```

## 用法

```bash
# 全流程
python auto_known_virus.py --stage all \
    --reads_dir 00b_HostDepletion/ --output_dir known_viruses/ \
    --ref_info ref_info.tsv --reference ref.fasta \
    --tool salmon --snpeff --snpgenie --threads 40 --jobs 4

# 分步
python auto_known_virus.py --stage detect    [参数]
python auto_known_virus.py --stage variants  [参数]
python auto_known_virus.py --stage full      [参数]
```

## 参数

### 通用

| 参数 | 默认 | 说明 |
|------|------|------|
| `--reads_dir` | 必需 | 清洁 reads 目录 |
| `--output_dir` | 必需 | 输出根目录 |
| `--ref_info` | 必需 | 参考元数据 TSV |
| `--reference` | 必需 | 参考序列 FASTA |
| `--stage` | all | detect/variants/full/all |
| `--resume` | — | 断点续传 |
| `--threads` | 40 | 线程数 |
| `--jobs` | 4 | 并行数 |

### detect (batch_virus_depth40.py)

| 参数 | 默认 | 说明 |
|------|------|------|
| `--tool` | bowtie2 | salmon/kallisto/bowtie2/bwa/... |
| `--batch_size` | 20 | 批大小 |
| `--coverage` | 10.0 | 全长覆盖度 % |
| `--ratio` | 0.3 | 泊松覆盖度比值 |
| `--meandepth` | 0.0 | 最小平均深度 |
| `--sp_thresh` | 95.0 | 物种 ANI 阈值 % |
| `--genes_cov` | — | 转录覆盖率文件 (双轨过滤) |

### variants (batch_virus_variants.py)

| 参数 | 默认 | 说明 |
|------|------|------|
| `--variant_caller` | freebayes | freebayes/ivar/lofreq |
| `--snpeff` | — | 启用 SnpEff 注释 |
| `--snpgenie` | — | 启用 SnpGenie 群体遗传 |

### full (batch_virus_full.py)

| 参数 | 默认 | 说明 |
|------|------|------|
| `--assembly_tools` | all | 组装工具 |
| `--min_covered` | 10.0 | 最小覆盖度 % |
| `--extra_args` | — | virus-full 额外参数 |

## 输出

```
known_viruses/
├── 1_FastViromeExplorer/       detect
│   └── summary/all_viruses.best.summary.tsv
├── 2_Virus_variants_Results/   variants
│   └── {virus}/{sample}/
└── 3_Virus_assemblies_final/   full
    └── {virus}/{sample}/*_final.fasta
```
