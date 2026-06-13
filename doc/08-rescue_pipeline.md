# rescue_pipeline.py — 三支路级联拯救

## 流程

```
接收: centroids + clusters.tsv + split_fastas (由 cluster_pipeline.py 产出)
    │
    ▼
分支 A: CheckV 并行评估 → completeness ≥ 90% → pass
    │ fail
    ▼
分支 C: Virseqimprover reads 迭代延伸
    Salmon 定量 → BBMap 提取 → SPAdes 组装 → 串联重复检测 → CheckV
    (自动聚合 cluster 内所有样本的 reads)
    │ fail
    ▼
分支 D: BLASTN 参考搜索 + CheckV + VSI 最后拯救
    │
    ▼
合并 (A+C+D) + vclust 最终去重 → HQ vOTU
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
| `-db, --blast-db` | — | BLAST 数据库 (分支 D) |
| `--virseqimprover-path` | — | Virseqimprover.py 路径 |
| `--salmon-bin` | salmon | Salmon 路径 |
| `-t, --threads` | 64 | 线程数 |
| `-j, --jobs` | 4 | VSI 并行数 |
| `--ani` | 0.95 | 最终 vclust ANI |
| `--qcov` | 0.85 | 最终 vclust QCOV |
| `--resume` | — | 断点续传 |

## 分支 C 多样本 reads 聚合

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
│   ├── branchA_pass.fasta     CheckV ≥90%
│   └── branchA_fail.fasta     → 分支 C
├── branch_c/
│   ├── merged_reads/          多样本合并 reads
│   ├── branchC_pass.fasta     VSI 成功
│   └── branchC_fail.fasta     → 分支 D
├── branch_d/
│   └── branchD_pass.fasta     BLASTN+VSI 成功
├── merged/all_HQ.fasta        A+C+D 合并
└── centroids/final_centroids.fasta  ★ 最终 HQ vOTU
```

## 断点续传

- 每个分支检查输出 pass-FASTA 是否存在
- 分支 C 内部检查每个 contig 的 scaffold.fasta
- `--resume` 跳过已完成分支
