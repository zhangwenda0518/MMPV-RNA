# cluster_pipeline.py — 聚类管道

## 流程

```
seqkit 长度过滤 (≥500bp)
    │
    ▼ (可选) CD-HIT 参考引导预聚类
vclust deduplicate 去重参考 → ref|/our| 标记 → vclust cd-hit 聚类
    ├─ Known: 代表含 ref|, 参考锚定
    └─ Novel: 代表含 our|, 无参考命中
    │
    ▼ vclust Leiden (仅 novel 部分)
→ centroids + clusters.tsv + split_fastas
```

## 用法

```bash
# 基础聚类
python cluster_pipeline.py -i input.fasta -o out/ -t 64 \
    --min-length 500 --ani 0.95 --qcov 0.85 --stop-after-vclust

# CD-HIT 参考引导
python cluster_pipeline.py -i input.fasta -o out/ -t 64 \
    --ref-genomes ICTV_plant_viruses.fasta \
    --cdhit-ani 0.95 --cdhit-qcov 0.85 --stop-after-vclust
```

## 参数

| 参数 | 默认 | 说明 |
|------|------|------|
| `-i, --input-fasta` | 必需 | 输入 FASTA |
| `-o, --output-dir` | 必需 | 输出目录 |
| `-t, --threads` | 64 | 线程数 |
| `--min-length` | 500 | 最小长度 bp |
| `--ani` | 0.95 | vclust ANI 阈值 |
| `--qcov` | 0.85 | vclust QCOV 阈值 |
| `--ref-genomes` | — | ICTV/NCBI 参考基因组 |
| `--cdhit-ani` | 0.95 | CD-HIT ANI 阈值 |
| `--cdhit-qcov` | 0.85 | CD-HIT QCOV 阈值 |
| `--skip-vclust` | — | 跳过 vclust (复用已有) |
| `--vclust-cluster-file` | — | 复用已有聚类 TSV |
| `--stop-after-vclust` | — | 聚类后停止 (编排器使用) |
| `--resume` | — | 断点续传 |

## 输出

```
{out}/
├── 1_seqkit/virus.candidate.fasta
├── 2_cdhit/                         (有 --ref-genomes 时)
│   ├── known_clusters/
│   ├── known_centroids.fasta
│   ├── known_association.tsv
│   └── novel_contigs.fasta
├── 3_vclust/
│   ├── vclust_clusters.tsv
│   └── split_fastas/
└── centroids/
    ├── final_centroids.fasta
    ├── known_ids.txt
    └── known_association.tsv
```

## CD-HIT 原理

vclust cd-hit 算法: 按长度降序排列, 最长序列为代表, 后续序列与代表比较。参考基因组通常更长, 自然成为代表, 碎片化 contig 只要与参考达 95% ANI, 85% QCOV 即归入同簇。

## 断点续传

- 检查 seqkit 和 vclust 输出文件
- CD-HIT 步骤始终重跑 (快速)
- `--resume` 跳过已完成步骤
