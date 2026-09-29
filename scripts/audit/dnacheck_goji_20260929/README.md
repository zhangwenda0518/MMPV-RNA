# 枸杞 DNA 病毒真伪判定 · 2026-09-29

问题：病毒组（RNA-seq）挖掘出的 DNA 病毒记录里，有多少是"真 DNA 病毒感染"而非 EVE/内源转录。

## 方法

1. 取 `plant_virus_DNA_RNA_final.tsv`（ICTV VMR 属优先裁定版）中 Nucleic_acid=DNA 且 Category ∈ {植物病毒, 植物病毒?} 的 **500 条记录**（转座子 185 / 非植物 / 未定已剔除）。
2. 从各队列 `03b_MergeSamples/all_sample_virus.fasta` 抽序列（495/500 成功）。
3. 与三个自组装枸杞参考基因组（zhonghua=中华 L. chinense、ningxia=宁夏 L. barbarum、heiguo=黑果 L. ruthenicum，`Lycium_EVE/00_edta/genomes/`）做三层比对：
   - minimap2 `-x asm5`（近期整合拷贝）
   - blastn dc-megablast e1e-10（中距同源）
   - blastn 经典任务 e1e-5（远缘同源，检出限 ~55% nt）
   联合库 goji3db（3.9 Gb，`~/goji_dnacheck_20260929/`）。
4. 判档：≥90% 且 ≥200nt = EVE 强同源（近期整合/同源拷贝）；无任何 ≥55%/150nt 证据 = 无基因组同源。

## 结果（500 条）

| 判定 | 条数 | 说明 |
|---|---|---|
| 基因组强同源（EVE 转录/同源拷贝） | **62** | ≥90% 匹配参考基因组 → 非真感染 |
| 无自身参考基因组（amarum + 3 个非枸杞队列） | **46** | 无法用基因组排除 |
| 无基因组同源（真病毒候选） | **392** | 其中 n≥2 的 66 条为优先复核池 |

- 候选 392 条中：220 条对参考基因组完全零命中、172 条只有 <55% 或 <100nt 的碎片命中 → 检测限以下的远缘分歧。
- 候选的种名宿主词几乎全为其他植物（maculaglycinis=大豆、venapetuniae=矮牵牛、venabougainvilleae=叶子花、dioscoreae=薯蓣…），而枸杞三个参考基因组本身带 213+ 个 pararetrovirus EVE 位点（Lycium_EVE 四象限终版 72/72/69）→ 大部分更可能是未测株系内源拷贝的远缘转录。

## 可确认的真 DNA 病毒

**1 种（2 条个体）**：Soymovirus（归 *Soymovirus maculaglycinis*），近全长 8420bp / 8312bp，均带 `_self_circ` 环化证据，覆盖 1411×（黑果 SRR36966021）与 1187×（中华 SRR16071801），与三参考基因组仅 58-62nt 碎片命中 → 全数据集唯一具备活性感染分子特征的 DNA 病毒。

## 待验证候选（优先池）

| contig | 队列 | n_samples | 覆盖 | 长度 |
|---|---|---|---|---|
| SRR31836787_NODE_19 *Begomovirus solanumflavusvietnamense* | 黑果 | **190（54.9% 检出率）** | 16× | 2569bp |
| SRR8293184_NODE_3 *Caulimovirus* sp. | 黑果 | 47 | 54.6× | 4995bp |
| CRR1440131_NODE_5 *Soymovirus maculaglycinis* | 宁夏 | 15 | 167.9× | 5023bp |
| SRR24320198_NODE_10 *Soymovirus maculaglycinis* | 中华 | 12 | 23.5× | 3610bp |

n=190 的 begomovirus avg_ANI=99.56，若经 read 级/RT-PCR 验证属实，将是继 cytorhabdovirus 与 PSTVd 之后枸杞第三大病毒，且是唯一 DNA 病毒——建议优先验证。

## 已知病毒侧（对照）

reference-based 47 种检出里唯一 DNA 病毒 = Banana streak IM virus（n=1，Badnavirus，著名内源化属，单样本不足以定性）。

## 产物

- 服务器：`~/goji_dnacheck_20260929/`（query.fasta、三轮比对输出、final_dna_verdict.tsv）
- 本地：本目录 `final_dna_verdict.tsv`（500 行终裁）、`genome_match_round1.tsv`（minimap2 中间表）
- 脚本：`scripts/audit/dnacheck_extract_queries.py`、`dnacheck_round1.py`、`dnacheck_final_verdict.py`、`dnacheck_fasubset.py`
