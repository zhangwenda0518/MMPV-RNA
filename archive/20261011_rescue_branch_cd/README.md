# 归档：rescue 分支 C/D（2026-10-11 移除）

## 归档内容

| 文件 | 说明 |
|---|---|
| `rescue_pipeline_4branch.py` | 移除前的完整四支路版 rescue_pipeline.py（v3.0–v3.1，原样备份） |
| `doc_08-rescue_pipeline_4branch.md` | 对应的四支路文档快照 |

现行版本只保留 **分支 A（CheckV ≥90% 直通）→ 分支 B（Virseqimprover reads 延伸）**，见 `virome_discovery_pipeline/rescue_pipeline.py` 与 `doc/08-rescue_pipeline.md`。

## 移除原因（2026-10-11 服务器实测裁定）

### 分支 C（BLASTN + ragtag 参考引导重建）：名不副实，实为"直通筛选"

- 参考引导重建路径（RGA / ragtag）在全部实测数据集上 **0 产出**：81 个属组（OneKP）、
  12/16/16 个属组（枸杞三物种）全部在 `blastdbcmd` 提取参考一步失败
  （`DB contains no accession info`——自动建库用内部 contig ID、`makeblastdb` 未加
  `-parse_seqids`；2026-06 的修复 161f074 只覆盖了外部库重建路径，自动建库路径至今未修）。
- 因此 C 的全部产出均为 **Completo 直通**（dc-megablast query 覆盖率 ≥98% 即原样放行，
  不做任何延伸）：
  - OneKP：9,212/9,212 条直通，0 条带 `method=` 标注；<1kb 1,013 条（11.0%）。
  - 枸杞 barbarum/ruthenicum/chinense：273/274/117 条全部直通，<1kb 占 23%/23%/19%。
- 直通判据本身有缺陷：qcov 是 query 覆盖率（"我被参考包含"），不是参考覆盖率
  （"我覆盖了完整基因组"），短片段天然满足；且最佳命中 95% 是库内其他 contig
  （聚类冗余），不构成完整性证据。
- 结论：C 实际语义是"片段被参考包含 → 放行入目录"，**不产出高质量序列**，
  与拯救目标相悖，移除。

### 分支 D（genus_len 属水平兜底）：零产出，与 B 重复

- OneKP 与枸杞全部数据集 D 均救回 0 条（branchD_pass.fasta 从未写出）。
- 属长判定在分支 B 内部已有同款回退（同 0.85 容忍、同长度判定），
  轮到 D 的只剩 VSI 彻底失败的原始长度 contig，几乎不可能达标。
- 节段属感知增强版（v3.1 `_genus_len_ok`）因编排器只递旧 2 列 genus_lens 表
  （CODE_AUDIT_2026-10-09 B3）从未生效。

## 影响与注意

- 移除后，A/B 均失败的 contig 直接记为未拯救（不再有片段直通入目录）。
- `rescue_report.tsv` 表头与 `_write_rescue_report` 四分支签名保持不变，
  `scripts/patch/regen_rescue_reports*.py` 对历史输出的重生成不受影响。
- 下游（report/integrate/panvirome/submission）只消费
  `centroids/final_centroids.fasta` 与 `rescue_report.tsv`，接口未变。
- 论文 Methods 中"four-branch cascade"表述需随重跑结果同步改写。

## 恢复方法

```bash
cp archive/20261011_rescue_branch_cd/rescue_pipeline_4branch.py \
   virome_discovery_pipeline/rescue_pipeline.py
```

（若需真正启用重建，须先修自动建库 ID 问题：`makeblastdb -parse_seqids` + 清洗 `|`，
或改为按 ID 直接从源 FASTA 提取参考；并把 Completo 判据改为参考覆盖率。）
