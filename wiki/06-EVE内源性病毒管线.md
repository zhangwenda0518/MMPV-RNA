# 06 · ⑥ EVE 内源性病毒管线 (endogenous_virus_pipeline)

**定位**: 宿主基因组侧的内源性病毒元件 (Endogenous Viral Elements) 筛查 — **与 reads 主链相互独立** (输入是宿主参考基因组 FASTA)。移植自 hi-fever eve_kingdom.py, 246 服务器千种植物实测。
**主编排器**: `eve_screen.py` → 单基因组 worker `eve_genome_scan.py` → 纯逻辑核心 `eve_scan_core.py` (坐标还原/位点合并/判定/汇总, 142 个单测)

## 三阶段筛查 (+汇总)

| 阶段 | 做什么 | 产物 (standard 根: `06_EVE/`) |
|---|---|---|
| discover | 基因组切块 (滑窗 50kb, `--overlap` 减半步长修跨界截断) → diamond blastx vs 病毒参考蛋白 → 位点合并 (MERGE_D=300nt) → samtools faidx 抽序列 (≥1.11, 开跑前体检) | `01_Loci/<NAME>/` (fna/chunks/s1_hits/loci.bed/loci.fa) |
| verdict | vs 植物/病毒拆分库双侧 bitscore 判定 (HOST_BS=VIRAL_BS=50) | `02_Verdict/<NAME>/` → **viral_supported / host_like / undetermined** |
| annotate | viral_supported + undetermined → RVDB 深度归属 | `03_RVDB/<NAME>/` |
| merge | 跨基因组汇总 | `04_Summary/`: `<NAME>_eve_summary.tsv` (12 列) + `kingdom_summary.tsv` + `family_by_genome.tsv` + ★ `eve_report.html` |

## 关键参数 (eve_screen.py)

| 参数 | 默认 | 说明 | 出处 |
|---|---|---|---|
| `-g <genome.fa>` / `-B batch.tsv` | 必填二选一 | 单基因组 / 批量 (NAME<TAB>路径, 重名开跑前报错; .tar.gz 包自动解) | eve_screen.py |
| `-t / -J` | — | 每基因组线程 × 并行基因组数 (256 核推荐 -J5 -t40) | 同上 |
| `--fast` | off | Stage1 快速模式 (约 6×, 位点 -54%; 千种全扫推荐) | 同上 |
| `--overlap` | off | 切块步长减半 (25kb 重叠) 修跨界截断 | 同上 |
| `--skip-s3 / --blastn-viroid` | off | 跳过 RVDB 层 / 附加类病毒 blastn 层 | 同上 |
| `--cmd-timeout` | — | 单条外部命令超时秒数 (防 diamond/seqkit 挂死占核) | 同上 |
| `--cleanup/--force/--merge/--list` | — | 中间清理 / 全量重跑 / 仅汇总 / 干跑 | 同上 |
| 数据库路径 | CLI > pipeline_config.yaml > env | `MMPV_EVE_REF_DB / MMPV_EVE_PV_DB / MMPV_EVE_ID2DIV / MMPV_EVE_RVDB_DB` | 同上 |

## 判别器 (后运行, B7: 吃 ③ 的发现输出)

`eve_distinguish/run_all.sh -D <③根> -A <1kp_assemblies>` 一键 6 步:
0 参考面板 (panel.fasta/baits.fa 随模块交付) → 1 建库 → 2 blastx → 3 结构/域
(`s1_decay_scan` 终止子富集 + `s2_domain_scan` MP/CP/AP/RT/RH 域记功) → 4 寄主
(`s2b_locus_scan` 候选回寄主) → 5/6 `s3_verdict` 结构退化×域架构二维判定 + 寄主否决
→ `s4_filter` 产出 ★ `dna_vs_eve_filter.tsv` (每候选一个 action: REMOVE / MOVE_EVE / KEEP_virus / REVIEW)

## 配套工具

`eve_report.py` (merge 产物 → 自包含 HTML 报告, 已过视觉验收) · `eve_qc_contam.py`
(污染/占比核查) · `eve_element_db.py` (以 04_Summary 为遍历键的元件库构建 —
summary 是唯一三表 join 好的单表) · `tests/` 142 个单测 (核心逻辑/判别逻辑/产物不变量)

## 上下游

- **上游 B6←①**: 宿主基因组 FASTA; **B7↔③**: 判别器消费 ③ 发现输出, filter 表回流指导候选取舍
