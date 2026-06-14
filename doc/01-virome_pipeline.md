# virome_pipeline.py — 主控流水线 / Master Orchestrator

> End-to-end automated virome pipeline controller. 11 independent stages, checkpoint-resume, sub-script progress passthrough.

宏病毒组端到端全自动主控，11 个独立阶段，断点续传，子脚本进度透传。

## 用法

```bash
python virome_pipeline.py --stage <stage> [参数]
```

## 11 个阶段

```
--stage clean              清洗 (Fastp + Seqkit + Clumpify)
--stage deplete            去宿主 (Kraken2 + Bowtie2/HISAT2 + rRNA)
--stage assembly           三工具组装 (Penguin + SPAdes + MEGAHIT)
--stage identification     六工具病毒鉴定
--stage cobra              COBRA 批量延伸
--stage cluster            聚类 (CD-HIT + vclust, 仅聚类不拯救)
--stage taxonomy           五工具分类 + R 共识
--stage host               宿主预测 (ICTV > RNAVirHost > PhaBOX2)
--stage checkv             按宿主 CheckV 预评估 (新增)
--stage rescue             宿主过滤 + 三支路级联拯救
--stage all                全流程串行
```

## 参数

| 参数 | 默认 | 说明 |
|------|------|------|
| `--input_reads` | — | 原始 FASTQ 目录 |
| `--output_dir` | **必需** | 项目输出根目录 |
| `--stage` | `all` | 运行阶段 |
| `--host-filter` | `Plant` | rescue 阶段目标宿主 |
| `--skip_clean` | — | 跳过清洗 |
| `--skip_depletion` | — | 跳过去宿主 |
| `--skip_clumpify` | — | 跳过 Clumpify |
| `--force` | — | 强制重跑 |
| `--kraken2_db` | — | Kraken2 宿主库 |
| `--host_align_db` | — | 宿主比对索引 |
| `--virus_db` | — | 病毒数据库根目录 |
| `--checkv_db` | — | CheckV 数据库 |
| `--phabox-db` | — | PhaBOX2 数据库 |
| `--prob-dir` | — | ICTV 宿主概率表 |
| `--ref-genomes` | — | ICTV/NCBI 参考基因组 |
| `--rrna` | 关闭 | 开启 rRNA 剔除 |
| `--rrna_tool` | `ribodetector` | `ribodetector` / `silva` |
| `--silva_index` | — | SILVA Bowtie2 索引 |
| `--assembler` | `penguin` | `penguin/megahit/rnaviralspades/all` |
| `--aligner` | `bowtie2` | `bowtie2/hisat2/minimap2` |
| `--seq_type` | `rna-short` | `dna-short/rna-short/nanopore/pacbio` |
| `--cluster_input` | — | 直接输入已合并病毒 FASTA |
| `-t, --threads` | `20` | 单任务线程数 |
| `-m, --memory` | `64` | 内存 GB |
| `-j, --jobs` | `2` | 并行任务数 |
| `--min-length` | `500` | 病毒最小长度 bp |
| `--ani` | `0.95` | 聚类 ANI 阈值 |
| `--qcov` | `0.85` | 聚类 QCOV 阈值 |

## 断点续传

默认跳过已完成步骤。`--force` 强制重跑。

```bash
# 中断后继续 — 直接重新运行相同命令
python virome_pipeline.py --stage all ...

# 只重跑某个阶段
python virome_pipeline.py --stage rescue --force ...
```

## 输出

```
out/
├── 00a_CleanData/         clean
├── 00b_HostDepletion/     deplete
├── 01_Assembly/           assembly
├── 02_Identification/     identification
├── 03_COBRA/              cobra
├── 04_CLUSTER/            cluster
├── 05_Taxonomy/           taxonomy
├── 06_HostPrediction/     host
├── 07_Checkv/             checkv
├── 08_Rescue/             rescue
└── orchestrator.log       全流程日志
```
