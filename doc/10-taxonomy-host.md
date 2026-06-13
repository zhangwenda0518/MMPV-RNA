# taxonomy + host — 分类注释与宿主预测

## virus_classifier2.py — 5 工具分类

### 5 个分类工具

| 工具 | 原理 |
|------|------|
| genomad | 深度学习病毒分类 |
| mmseqs | 蛋白序列比对分类 |
| VITAP | 病毒蛋白分类 |
| ACVirus | 古菌病毒分类 |
| vcontact3 | 蛋白簇网络分类 |

### 用法

```bash
python virus_classifier2.py \
    -g centroids.fasta -s sample_name \
    -t genomad,mmseqs,VITAP,ACVirus,vcontact3 \
    -o out/ -p 64 --db-dir /db/virus_db/ [-f]
```

### 参数

| 参数 | 默认 | 说明 |
|------|------|------|
| `-g, --genomes` | 必需 | 输入 FASTA |
| `-s, --sample` | 必需 | 样本名 |
| `-t, --tools` | all | 逗号分隔工具列表 |
| `-o, --output-dir` | ./classify_output | 输出目录 |
| `-p, --threads` | 20 | 线程数 |
| `--db-dir` | ~/database/virus-db | 数据库目录 |
| `-f, --force` | — | 强制重跑 |

### R 共识整合

```bash
Rscript virus_classifier_analysis14.R \
    --combined combined_taxonomy.tsv --output out/
```

优先级: vcontact3 > vitap > acvirus > mmseqs > genomad
共识: 最多非 NA 层级的工具获胜

### 输出

```
{out}/
├── {sample}.virus_classed/
│   ├── {sample}_{tool}_taxonomy.tsv    各工具单独输出
│   └── {sample}_combined_taxonomy.tsv  5工具合并
└── integrated/
    └── final_integrated_classification.tsv ★ 最终分类
        列: contig_id, Realm, Kingdom, Phylum, Class,
             Order, Family, Genus, Species, Determination_Method
```

---

## run_host_prediction.py — 宿主预测决策树

### 三工具集成

| 工具 | 原理 | 覆盖范围 |
|------|------|----------|
| ICTV (C9) | 官方分类库查找 | 已知病毒宿主 |
| RNAVirHost | 全生态位模型 | plant/animal/fungi/bacteria |
| PhaBOX2 | CRISPR+AAI 网络 | 噬菌体 |

### 决策树

```
1. Class ∈ 噬菌体纲 → Bacteria (硬规则)
2. ICTV == RNAVirHost ≠ Unknown → 直接采用
3. 分歧时 PhaBOX2 决胜 → 支持哪方用哪方
4. 全部分歧/缺失 → ICTV > RNAVirHost > PhaBOX2
```

### 用法

```bash
python run_host_prediction.py \
    -i centroids.fasta --tax final_integrated_classification.tsv \
    -o out/ -t 64 --mode all \
    [--phabox-db /db/phabox/] [--prob-dir /db/cross_analysis/] [-f]
```

### 输出

```
{out}/
├── ensemble_host_summary.tsv         ★ 核心输出
│   列: contig_id, Class, Host_ICTV, pred|L1,
│        Final_Host, Decision_Method
└── host_classified_fasta/
    ├── Plant.classified.fasta
    ├── Bacteria.classified.fasta
    ├── Fungi.classified.fasta
    ├── Animal.classified.fasta
    └── Unknown.classified.fasta
```
