# 06 · EVE 内源病毒（集成于 ① 与 ③）

> **状态**：2026-10-25 重构。EVE 不再单独算一条管线，拆成两段接到宿主库构建（①）与
> 病毒发现的 DNA 候选分析（③）。旧的 `endogenous_virus_pipeline/` + `eve_distinguish/`
> 已归档至 `archive/20261011_endogenous_virus_pipeline/`，仅用于历史结果复现。

## 两段式接入

| 段 | 所在管线 | 模块 | 触发 | 产物 | 进哪个 HTML |
|---|---|---|---|---|---|
| ① 基因组侧 | `public_metadata_pipeline` | `eve_genome/` | `--eve-scan`（`eve` 阶段） | `<host_reference>/eve/04_Summary/eve_report.html`、`05_Elements/` | `Pipeline_Summary_Report.html` → **S7: EVE (genome)** |
| ③ 转录组侧 | `virome_discovery_pipeline` | `eve_discriminate/` | 09b 默认开（`--no-eve-discriminate` 关） | `09b_Analysis_Verify/virus_validation/eve/{eve_calls.tsv,eve_report.html}` | `10_Reports/pipeline_report.html` → 09b 卡片 **EVE discrimination (before KEEP)** |

**为什么 ① 的掩码在 hostdb 之前、分析汇总在之后**：掩码必须发生在建索引之前，宿主库才不会
包含 EVE 序列（否则 EVE 来源的 reads 会被当宿主丢掉）；而"分析汇总"不改变已建好的库，放在
`hostdb` 之后，失败也不影响建库产物。

**为什么 ③ 在 KEEP 之前**：KEEP 由 `integrate_rescue_evidence.py` 判定，`class_KEEP.fasta`
直接从它抽取。放在 KEEP 之后只能给已进目录的 EVE 贴标签；放在之前，EVE 证据才能参与裁决
（`EVE_like` → REVIEW）。只降不删，保留人工复核。

## ① 基因组侧（`public_metadata_pipeline/eve_genome/`）

```
eve_screen.py          三阶段扫描: discover → verdict → annotate → merge
eve_genome_scan.py     单基因组 worker
eve_scan_core.py       坐标还原/位点合并/双侧判定
eve_element_db.py      元件层: 聚类 → regions → consensus → RM/CDD/phylo
eve_scan_report.py     阶段 HTML → 04_Summary/eve_report.html
eve_prot_lib.py        consensus 核苷酸库 → 六框蛋白库 (.sixframe.faa)
```

- **阶段位置**：`public_data_pipeline.py` 的 `eve` 阶段（`hostref → hostdb → eve → report`）；
  `build_host_pipeline.py` 的 `eve-report` 阶段（`genome-down → eve-mask → hostdb → eve-report`）。
- **默认关**：基因组扫描（滑窗 blastx）耗时，只有显式 `--eve-scan` 才跑；不给时阶段跳过并记 done。
- **元件层**：`--eve-elements` 追加 `collect,cluster,consensus,cdd`，并自动产出六框蛋白库，
  供 ③ 的 Δbit-score 使用。
- **报告章节**：`Pipeline_Summary_Report.html` 的 S7 卡片显示基因组数、viral-supported、
  undetermined、host-like、元件数，并链接 `eve_report.html` 与蛋白库路径。

```bash
# 九阶段公共数据管线
python public_metadata_pipeline/public_data_pipeline.py --species <拉丁名> --taxid <ID> \
    --stage hostref,hostdb,eve,report --eve-scan --eve-elements \
    --work-dir <项目>/01_PublicData

# 独立宿主构建管线
python public_metadata_pipeline/build_host_pipeline.py --species <拉丁名> --taxid <ID> \
    --stage all --eve-scan --eve-elements --work-dir out/host
```

## ③ 转录组侧（`virome_discovery_pipeline/eve_discriminate/`）

```
eve_discriminate/
├── no_ref_genome/      序列判别: ML(16科) / Δbit-score / Nanoviridae 无监督 / 双库分类 / action 映射
├── with_ref_genome/    整合证据: 参考基因组回比 / 位点架构 / 转录本嵌合
├── eve_classify_report.py   阶段 HTML 汇总
├── eve_final_score.py       多路证据 → final_tiers.tsv
├── eve_final_report.py      自适应最终 HTML
├── eve_bam_coverage.py      [可选] BAM 覆盖注释
├── eve_reconstruct.py       位点蛋白重建 (miniprot + GeneWise)
└── eve_html.py              报告共用渲染层
```

- **位置**：`run_analysis_verify()` 内，`validate_rescue_cdd.py` 之后、`integrate_rescue_evidence.py`
  之前；候选为 `HQ_plant_viruses.fasta`，科标签来自 `final_integrated_classification.tsv`。
- **默认通道**：16 科 ML（模型在 `database/eve_ml/`）；可选 Δbit-score 需要
  `--eve-faa`（① 的 `.sixframe.faa`）+ `--eve-exo-dmnd`。
- **否决规则**：`integrate_rescue_evidence.py --eve-calls`，`EVE_like` 且 `verdict=KEEP` →
  `REVIEW`，输出列新增 `eve_call / eve_score / eve_flag`。
- **报告区块**：09b 卡片追加 EVE discrimination 区块（打分候选数 / EVE_like / KEEP→REVIEW 数），
  并复制 `eve_report.html` 为 `EVE_discrimination.html` 供点击查看。

```bash
python virome_discovery_pipeline/virome_pipeline.py --stage analysis_verify \
    --output_dir <发现项目根> \
    --eve-faa <host_reference>/eve/05_Elements/eve_consensus_library.sixframe.faa \
    --eve-exo-dmnd <ICTV 蛋白库>.dmnd
```

## 模型与数据

| 资产 | 位置 |
|---|---|
| 16 科 ML 模型（权威） | `database/eve_ml/models/` |
| OOF AUC / Youden J 阈值 | `database/eve_ml/index.json` |
| 1,030 维特征规格 | `database/eve_ml/FEATURE_SPEC.md` |
| 训练/构建脚本 | `database/eve_ml/build/` |
| 自检样本（1,920 条） | `virome_discovery_pipeline/eve_discriminate/no_ref_genome/sample/` |

16 科模型 OOF AUC 0.838–0.976（Caulimoviridae 0.976、Geminiviridae 0.974）。

## 使用边界

1. ML/Δbit 需要科标签（③ 的 taxonomy/DIAMOND 分类）
2. 分数是排序分，不是判决；最终靠基因组整合证据
3. 新近整合 EVE 会漏检（设计固有边界）
4. PPV 未知（阳性标签非金标准）
5. Nanoviridae 走无监督（位点太少，训不了监督模型）
6. `EVE_like` 只把 KEEP 降为 REVIEW，不直接删除：序列模型有假阳，删数据需要更强证据
