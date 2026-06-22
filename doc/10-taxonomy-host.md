# taxonomy + host — 分类注释与宿主预测 / Taxonomy & Host Prediction

> 8-tool parallel classification + weighted-vote consensus → 8-rank taxonomy with per-rank agreement tracking.

## virus_classifier.py — 8工具并行分类

### 8 个分类工具

| 工具 | 原理 |
|------|------|
| genomad | 深度学习病毒分类 |
| metabuli | k-mer 分类 |
| CAT | 蛋白比对分类 |
| diamond_lca | LCA 分类 |
| VITAP | 病毒蛋白分类 |
| mmseqs | 蛋白比对分类 |
| ACVirus | 古菌病毒分类 |
| vcontact3 | 蛋白簇网络分类 |

### 用法

```bash
python virus_classifier.py \
    -g centroids.fasta -s sample_name \
    -t all \
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
| `--validate-only` | — | 仅验证已有结果 (对比 combined 与重算) |

### 共识整合 (Python 或 R)

```bash
# Python
python virus_classifier_analysis.py --combined combined_taxonomy.tsv -o out/

# R
Rscript virus_classifier_analysis.R --combined combined_taxonomy.tsv -o out/
```

共识算法: 逐 rank 加权投票 (深度权重 Species×128 > Genus×64 > ... > Realm×1)
+ ICTV 后缀校验 + subrank 清除 + genus-species 一致性修复

### 输出

```
{out}/
├── {sample}.classed/
│   ├── {sample}_{tool}_taxonomy.tsv         各工具标准化输出
│   └── {sample}_combined_taxonomy.tsv       8工具合并 (virus_classifier.py)
└── {sample}.integrated/
    ├── final_integrated_classification.tsv  ★ 最终分类 + 工具一致信息
    │   列: contig_id, primary_tool, completeness, confidence,
    │        Realm..Species, Realm_agree..Species_agree
    │   agree 格式: "5/7: tool1,tool2,..." (同意数/总数: 工具列表)
    ├── intersection_upset.pdf                UpSet 交集图
    ├── agreement_rates.pdf                   各工具各 rank 一致率
    ├── consensus_summary.pdf                 完备度分布 + rank 填充率
    ├── standardized_{tool}.tsv              各工具标准化结果
    └── comparison_{level}.tsv               逐 contig 逐 rank 详细对比
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
