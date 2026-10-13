# rescue_pipeline.py — 两支路级联拯救 / Two-Branch Cascade Rescue

> Branch A (CheckV ≥90% 直通) → Branch B (Virseqimprover VSI reads 延伸). Multi-sample read aggregation, final vclust dedup.
>
> **v3.2 (2026-10-11)**: 原分支 C (BLASTN+ragtag 参考引导重建) 与分支 D (genus-length 兜底)
> 已移除归档至 [`archive/20261011_rescue_branch_cd/`](../archive/20261011_rescue_branch_cd/README.md)。
> 原因: C 的重建路径因自动建库缺 `-parse_seqids` 全部失败, 实际产出退化为"片段被参考包含即直通"
> (OneKP 9,212/9,212 直通、<1kb 占 11%; 枸杞三物种全部直通、<1kb 占 19–23%), 与高质量序列目标相悖;
> D 在全部实测数据集零产出且与 B 内置的属长回退重复。移除后 A/B 均失败的 contig 直接记为未拯救。

## 流程 / Workflow

```
接收: centroids + clusters.tsv + split_fastas (由 cluster_pipeline.py 产出)
    │
    ▼
分支 A: CheckV 并行评估 → completeness ≥ 90% → pass (HQ)
    │ fail
    ▼
分支 B: Virseqimprover reads 迭代延伸
    Salmon 定量 → BBMap 提取 → SPAdes 组装 → CheckV 复评
    (自动聚合 cluster 内所有样本的 reads, 最多 10 样本)
    CheckV=NA/不达标但长度 ≥85% 属均长 → 属长回退拯救 (节段属按范围判定)
    │ fail
    ▼
记为未拯救 (仅 branchB_fail.fasta 存档, 不入 HQ 目录)
    │
    ▼
合并 (A+B) + vclust 最终去重 (ANI 0.95 / qcov 0.85) → HQ vOTU
```

## 用法

```bash
python rescue_pipeline.py \
    -c centroids.fasta \
    --clusters-tsv vclust_clusters.tsv \
    --split-dir split_fastas/ \
    -o out/ -fq reads/ -cv checkv_db/ \
    -t 64 -j 20 [--resume]
```

## 参数

| 参数 | 默认 | 说明 |
|------|------|------|
| `-c, --centroids` | 必需 | centroids FASTA |
| `--clusters-tsv` | 必需 | vclust 聚类结果 TSV |
| `--split-dir` | 必需 | per-cluster 拆分目录 |
| `-o, --output-dir` | 必需 | 输出目录 |
| `-fq, --fastq-dir` | 必需 | reads 目录 |
| `-cv, --checkv-db` | 必需 | CheckV 数据库 |
| `--taxonomy-tsv` | — | 分类 TSV (B 属长回退判定用) |
| `--genus-len` | — | 属平均长度表 (B 属长回退判定用) |
| `--virseqimprover-path` | — | Virseqimprover.py 路径 |
| `--salmon-bin` | salmon | Salmon 路径 |
| `--max-vsi-samples` | 10 | VSI 最大合并样本数 |
| `--min-vsi-len` | 2000 | VSI 最小 contig 长度 bp (短于此记为未拯救) |
| `--checkv-threshold` | 90.0 | CheckV 通过阈值 |
| `-t, --threads` | 64 | 线程数 |
| `-j, --jobs` | 4 | VSI 并行数 |
| `--ani` / `--qcov` | 0.95 / 0.85 | 最终 vclust ANI / QCOV |
| `--resume` | 开 | 断点续传 |
| `-db, --cdhit-fa, --genus-tolerance` | — | **已废弃** (随分支 C/D 归档, 接受但忽略) |

## 分支 B 多样本 reads 聚合

Virseqimprover 原本只接受单样本 reads。rescue_pipeline.py 自动从 cluster 内所有成员 contig 名提取样本前缀, 合并所有样本的 reads, 显著提升低丰度病毒的延伸成功率。

```
cluster_0 成员:
  GS-3_clean_0   → GS-3
  NX-4_clean_12  → NX-4
  QH_clean_5     → QH

→ 合并 GS-3 + NX-4 + QH 三个样本的 reads → Virseqimprover
```

## 输出

```
{out}/
├── branch_a/
│   ├── branchA_pass.fasta     CheckV ≥90% 直通
│   ├── branchA_fail.fasta     → 分支 B
│   └── branch_a_prepass.txt   编排器预通过计数 (免拯救, 候选表外)
├── branch_b/
│   ├── merged_reads/          多样本合并 reads (跑完自动清理)
│   ├── branchB_pass.fasta     VSI 延伸 + CheckV/属长判定通过
│   └── branchB_fail.fasta     未拯救存档
├── merged/all_HQ.fasta        A+B 合并
├── centroids/final_centroids.fasta  ★ 最终 HQ vOTU (vclust 去重)
├── rescue_report.tsv          逐条追踪 (contig_id,input_len,final_len,branch,...)
└── rescue_summary.md          分支统计 + 未拯救原因分布
```

## 断点续传

- 每个分支检查输出 pass-FASTA 是否存在
- 分支 B 内部检查每个 contig 的 scaffold.fasta / run.ok / .vsi_error
- `--resume` 跳过已完成分支 (默认开启)

## 历史: 被归档的分支 C / D

| 分支 | 策略 | 归档原因 |
|------|------|----------|
| C | dc-megablast → VirASCA 三档 → RGA/ragtag 参考引导重建 | 重建路径 0 产出 (blastdbcmd 提取参考全失败), 退化为 qcov≥98% 直通筛选, 放行大量 sub-1kb 片段 |
| D | genus_len 属水平长度兜底 (≥85% 属均长) | 全部实测数据集 0 产出; 功能与 B 内置属长回退重复; 节段属增强版未接线 |

恢复方法与完整裁定证据见 [archive/20261011_rescue_branch_cd/README.md](../archive/20261011_rescue_branch_cd/README.md)。
