# Virseqimprover.py — 病毒contig迭代延伸引擎 / Iterative Extension Engine

> Salmon quantification → BBMap extraction → SPAdes assembly → loop until CheckV > 90% or no length gain. Circularity detection.

## 流程 / Workflow

```
输入 scaffold + reads
    │
    ▼ Salmon 定量 → 提取末端低覆盖区域
    ▼ BBMap 提取映射 reads
    ▼ SPAdes --trusted-contigs 迭代组装
    ▼ 串联重复检测 (连续 +Δ×3 拦截)
    ▼ CheckV 完整性检查
    │
    ▼ scaffold.fasta (延伸后序列)
```

## 用法

```bash
python Virseqimprover.py \
    -1 reads_1.fq.gz [-2 reads_2.fq.gz] \
    -scaffold ref.fasta -o out/ [参数]
```

## 参数

| 参数 | 默认 | 说明 |
|------|------|------|
| `-1` | 必需 | 正向 reads (FASTQ/FASTA, .gz) |
| `-2` | `""` | 反向 reads (空 = 单端模式) |
| `-scaffold` | 必需 | 输入 scaffold FASTA |
| `-o` | 必需 | 输出目录 |
| `-t, --threads` | 16 | 线程数 |
| `-salmon` | salmon | Salmon 二进制路径 |
| `-checkv_db` | `""` | CheckV 数据库 |
| `-spadeskmer` | default | SPAdes k-mer |
| `-minOverlapCircular` | 5000 | 环状检测最小重叠 |
| `-minIdentityCircular` | 95 | 环状检测最小 identity |
| `-readFrac` | 0 | Salmon read fraction |
| `-minSuspiciousLen` | 1000 | 可疑区域最小长度 |
| `-h, --help` | — | 帮助 |

## 输出

```
{out}/
├── scaffold.fasta              最终延伸产物
├── salmon-res/                 Salmon 定量结果
├── salmon-mapped.sam           映射 reads
├── spades-res/scaffolds.fasta  SPAdes 组装
├── checkv_tmp/                 CheckV 结果
└── *.log                       运行日志
```

## 串联重复检测

延伸过程中检测 scaffold 是否无限增长 (串联重复):
- 连续 3 次 +Δ → 拦截, 回退到备份
- 截断 300bp → 再延伸 → 截断 500bp → 优雅退出
- 防止 SPAdes 在重复区域循环组装
