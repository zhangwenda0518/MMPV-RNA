# validate_novel_viruses.py — 交叉验证

基于 taxonomy 分类层级判断病毒新颖性，无需 BLASTN。

## 分类规则

```
Species ≠ NA              → ★ known          已知病毒
Genus ≠ NA, Species = NA  → ★★ novel_species 新种
Family ≠ NA, Genus = NA   → ★★ novel_genus   新属
Order ≠ NA, Family = NA   → ★★★ novel_family 新科
Class ≠ NA, Order = NA    → ★★★ novel_order  新目
全是 NA                    → ★★★ truly_novel  全新
```

CD-HIT 已知标记的 centroids 直接归为 ★ known。

## 用法

```bash
python validate_novel_viruses.py \
    -i rescue_Plant/centroids/final_centroids.fasta \
    --taxonomy final_integrated_classification.tsv \
    --cdhit-known known_association.tsv \
    --clusters-tsv vclust_clusters.tsv \
    --host ensemble_host_summary.tsv \
    -o 09_Validation/
```

## 参数

| 参数 | 默认 | 说明 |
|------|------|------|
| `-i, --input` | 必需 | vOTU centroids FASTA |
| `--taxonomy` | 必需 | 分类注释 TSV |
| `--cdhit-known` | — | CD-HIT known_association.tsv |
| `--clusters-tsv` | — | vclust_clusters.tsv (频率计算) |
| `--host` | — | 宿主预测 TSV |
| `--known-summary` | — | 已知病毒 summary (丰度补充, 可选) |
| `--ref-info` | — | 参考元数据 (物种名补充, 可选) |
| `-o, --output-dir` | 09_Validation | 输出目录 |

## 病毒频率

从 vclust_clusters.tsv 提取每个 cluster 的成员 contig, 取样本前缀去重:

```
cluster_0 成员:
  GS-3_clean_0   → GS-3
  NX-4_clean_12  → NX-4
  QH_clean_5     → QH

frequency = 3, samples = "GS-3,NX-4,QH"
```

## 输出

```
{out}/
├── novel_viruses.annotated.tsv       ★ 核心: 全部 vOTU 分类+频率
│   列: contig_id, category, frequency, samples,
│        Final_Host, Species, Genus, Family, ...
├── final_virus_catalog.fasta         按 category 标记的合并目录
└── validation_report.html            Plotly.js 交互报告
    ├── 饼图: ★已知/★★新种/★★新属/★★★新科+/★★★全新
    ├── 分类层级柱状图
    ├── 宿主分布柱状图
    └── 频率分布直方图 (singleton/2-5/6+)
```
