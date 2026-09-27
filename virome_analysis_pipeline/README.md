# Known Virus Analysis Pipeline

9-stage automated pipeline for known virus detection, variant analysis, full-length assembly, and post-hoc characterization.（正选择分析已于 v3.0 迁移至 `virome_phylo_pipeline/`，见文末注。）

## Quick Start

```bash
python auto_known_virus.py --stage all --filter \
  --reads_dir clean_reads/ \
  --output_dir ./virus_analysis/ \
  --ref_info ref_info.tsv \
  --reference ref.fasta \
  --genes_cov virus_genes_cov.tsv \
  --tool salmon --threads 60 --align_threads 8 \
  --variant_caller freebayes --snpeff --snpgenie \
  --resume
```

## Pipeline Flow

```
Input: clean_reads/ + ref_info.tsv + ref.fasta

[1/9] batch_virus_depth.py        → detect      Salmon/Kallisto/Bowtie2 定量 + Poisson 打假
[2/9] filter_summary.py           → filter      高置信度过滤
[3/9] batch_virus_variants.py     → variants    FreeBayes/iVar/LoFreq 变异检出 + 共识序列
[4/9] virus_vcf_pipeline.py       → post        VCF 可视化 + SnpEff 宏 + MAF 瀑布 + SnpGenie dN/dS
[5/9] batch_virus_full.py         → full        12 步全长基因组构建 (OmniVirusAssembler V9.0)
       └── virus-full.py (engine)
[6/9] extract_full_fasta.py       → extract     最长 contig 提取 + 参考 N-fill
[7/9] virus_auto_pipeline.py      → similarity  SDT 相似性热图 + 层次聚类
[8/9] batch_virema_dvg.py         → dvg         ViReMa DVG/重组检测 + Circos
       └── virema_summary_report.R   Circos 4-track
[9/9] generate_pipeline_report.py → report      交互式 HTML 综合报告 + AI 解读提示词
```

> **注**: 正选择分析 (`capheine` / HyPhy FEL·MEME·BUSTED·PRIME) 已于 v3.0 迁移至
> `../virome_phylo_pipeline/`（STAGE_ORDER 的 `capheine` 阶段：cawlign 密码子比对 →
> IQ-TREE 基因树 → HyPhy 系列检测 → DRHIP；入口 `capheine_pipeline.py` /
> `selection_suite.py`）。本管线自 v3.0 起为 9 阶段，阶段-脚本映射以
> `auto_known_virus.py` 的 `STAGES` 注册表为准。

## Stage Control

```bash
# Single stage
--stage detect|filter|variants|full|extract|post|similarity|dvg|report

# Full pipeline with filtering
--stage all --filter

# Resume from checkpoint
--resume

# Preview only
--dry_run
```

## Script Index

| Script | Stage | Lines | Function |
|--------|:----:|:----:|----------|
| `auto_known_virus.py` | — | 755 | 9-stage orchestrator |
| `batch_virus_depth.py` | 1 | 832 | Salmon pseudo-alignment + Poisson filtering |
| `filter_summary.py` | 2 | 236 | Coverage/depth/reads threshold filtering |
| `batch_virus_variants.py` | 3 | 1229 | FreeBayes/iVar/LoFreq 变异检出 + 共识序列 |
| `virus_vcf_pipeline.py` | 4 | 186 | post 阶段入口：bcftools merge + VCF2Dis + VCF2PCA |
| `batch_virus_full.py` | 5 | 240 | Assembly task scheduler |
| `virus-full.py` | 5 | 1120 | 12-step de novo assembly engine |
| `extract_full_fasta.py` | 6 | 123 | Extract longest contigs |
| `batch_plot_virus_depth.py` | 4 组件 | 784 | 3-in-1 viz (depth/freq/meta) |
| `virus_variants_analyzer.py` | 4 组件 | 477 | VCF landscape + PCA |
| `snpeff_analysis.py` | 4 组件 | 310 | SnpEff macro stats + OncoPrint |
| `snpeff2maf.py` | 4 组件 | 333 | VCF to MAF conversion |
| `viral_maftools.R` | 4 组件 | 233 | maftools waterfall + lollipop |
| `snpgenie_master.py` | 4 组件 | 671 | dN/dS + Auto-K 3D PCA + Kruskal-Wallis |
| `virus_auto_pipeline.py` | 7 | 1265 | Pairwise SDT similarity panorama |
| `batch_virema_dvg.py` | 8 | 544 | ViReMa DVG detection |
| `virema_summary_report.R` | 8 | 335 | Circos 4-track recombination plot |
| `generate_pipeline_report.py` | 9 | 374 | HTML summary + AI prompts |
| `consensus_extract.py` | util | 231 | Consensus sequence QA + N-filling |
| `panvirome_prevalence.py` | 汇总 | — | 序列级发生率总表 (发生率 × 分类 × CheckV × 属平均长度) |
| `panvirome_novel_known.py` | 汇总 | — | 新病毒 vs 已知病毒综合判定表 (核酸/蛋白/CDD + 02blast + 鉴定/分类工具) |

> 已迁出脚本：`capheine_pipeline.py` / `gbk_extractor.py` / `visual_codon_miner.py`
> 现位于 `../virome_phylo_pipeline/`（正选择分析）。

## Pan-Virome 汇总

`panvirome_prevalence.py` 是 detect 阶段之后的下游汇总工具，横跨 discovery（rescue/taxonomy/CheckV）与 analysis（detect）两条管线。它以 contig_id（vOTU 序列名）为最小单元统计序列级发生率，而不是按分类工具的物种标签聚合，避免物种名不确定带来的误导。

```bash
python panvirome_prevalence.py \
    --base ~/virus/data-2026/data-test \
    --species barbarum,ruthenicum,chinense,amarum \
    --detect-tag bowtie2 \
    --out prevalence_full_table.tsv
```

`panvirome_novel_known.py` 对每个物种的 rescue vOTU 做新病毒 vs 已知病毒综合判定，三层证据（核酸 blastn identity、蛋白 mmseqs identity、CDD 结构域）+ 流程溯源（02 鉴定阶段 blast 结果与工具、05 分类 primary_tool）。CDD 结构域用于区分真病毒（RdRp/RT）与宿主基因/逆转座子，是戳穿"分类工具误标"的硬证据。中间产物（blastn/mmseqs/cdd 结果）存在时自动复用，可用 `--rerun-*` 强制重跑。

```bash
python panvirome_novel_known.py \
    --base ~/virus/data-2026/data-test \
    --species barbarum,ruthenicum,chinense,amarum \
    --out novel_vs_known_final.tsv
```

`panvirome_novel_known_report.py` 读 `novel_vs_known_final.tsv`，生成自包含的深色 SCI 风格 HTML 报告（判定分布图 + 新病毒候选清单 + 方法说明），可直接在浏览器打开。

```bash
python panvirome_novel_known_report.py \
    --tsv novel_vs_known_final.tsv \
    --out novel_vs_known_report.html
```

## Dependencies

**Python packages**: polars, pysam, biopython, pandas, matplotlib, seaborn, scipy, scikit-learn, tqdm, colorlog

**External tools**: bowtie2, samtools, freebayes, bcftools, snpEff, salmon/kallisto, megahit/spades, mafft, iqtree, pandepth, ViReMa, BBMerge（正选择相关 hyphy/cawlign/drhip/multiqc 随 capheine 迁至 `../virome_phylo_pipeline/`）

**R packages**: circlize, ComplexHeatmap, maftools, ggplot2, dplyr, tidyr, viridis

## Archives

- `archive/` — 52 historical versions of batch_virus_depth*.py
- `utils/` — 6 standalone tools (snpeff_build, consensus_to_proteins, etc.)
