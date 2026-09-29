# eve_screen 自研管线 · 枸杞三基因组重分析 · 2026-09-29

用 MMPV-RNA 自研 EVE 管线 `endogenous_virus_pipeline/eve_screen.py`（与 eve_kingdom 250+ 基因组同一套方法与数据库）对枸杞三个自组装参考基因组做正向 EVE 筛查，替代/交叉验证此前 detectEVE + hi-fever + CAULIFINDER（Lycium_EVE 项目）的结果。

## 运行

- 服务器 246：`~/eve_screen_goji3_20260929/`（01_Loci / 02_Verdict / 03_RVDB / 04_Summary）
- 命令：`python3 eve_screen.py -B goji3_batch.tsv -o ~/eve_screen_goji3_20260929 -t 64`（eve_screen.py 与 eve_scan_core.py 服务器/本地 md5 一致，数据库取 pipeline_config.yaml 默认四库：ncbi-virus_ref.pep / plant_virus_final / id2div_final / U-RVDBv31）
- 基因组：zhonghua=L. chinense 中华、ningxia=L. barbarum 宁夏（clean）、heiguo=L. ruthenicum 黑果（clean）
- 类病毒层：管线自带 word 7 全基因组版代价 ~17h/基因组，改用 Lycium_EVE 协议快筛（word 11 + dust no，按染色体分片并行，过滤 aln≥80nt 且 qcov≥50%）；word 7 精查只对阳性 ±10kb 后补。

## 结果总览（04_Summary/kingdom_summary.tsv）

| 基因组 | viral_supported | host_like | undetermined |
|---|---|---|---|
| zhonghua | 262 | 62,693 | 49 |
| ningxia | 493 | 79,559 | 646 |
| heiguo | 466 | 96,331 | 1,275 |

viral_supported 位点科级构成（family_by_genome.tsv）：

| family | heiguo | ningxia | zhonghua |
|---|---|---|---|
| Caulimoviridae | 205 | 361 | 251 |
| Phage | 192 | 96 | 0 |
| Other_viral | 59 | 32 | 10 |
| NCLDV | 3 | 2 | 0 |
| Partitiviridae | 2 | 1 | 0 |
| Chrysoviridae | 2 | 1 | 0 |
| Amalgaviridae | 1 | 0 | 1 |
| Trichoviridae / Geminiviridae | 1 / 1 | 0 / 0 | 0 / 0 |

要点：
1. **Caulimoviridae 是绝对主力**（817/1,221 = 67%），与 Lycium_EVE 四象限结论一致。
2. **新发现：Phage 源位点**（黑果 192、宁夏 96、中华 0）——枸杞基因组装里的噬菌体/细菌源序列（可能来自 endophyte 污染或已整合的噬菌体 DNA），三个基因组不均衡，值得抽查。
3. **RNA 病毒 EVE 罕见但存在**：Partitiviridae / Chrysoviridae / Amalgaviridae / Trichoviridae 各 1-2 个位点（分体病毒/金色病毒/ amalgavirus 的内源化在植物基因组里罕见，每个都是个案例外）。
4. heiguo 有 1 个 Geminiviridae 位点，可对上 Lycium_EVE 的 17 个 Burewala 残迹位点（蛋白层可见的子集）。

## 与 Lycium_EVE 四象限（hi-fever+CAULIFINDER）的交叉

四象限全部位点（v2.bed，含 1/2a/2b/x_other）对扫描 verdict 的落点：

| 基因组 | 四象限位点 | 扫描=viral_supported | 扫描=host_like | 扫描完全未见 |
|---|---|---|---|---|
| zhonghua | 326 | 27 | 134 | 163 |
| ningxia | 413 | 28 | 143 | 240 |
| heiguo | 504 | 35 | 157 | 309 |

- 两边一致的高置信 EVE：27-35 个/基因组。
- 134-157 个是扫描找到了蛋白同源、但因**植物/TE 侧 bitscore 更高**被判 host_like——典型即 StubV 位点（zhonghua Chr01:7534541-7535836，植物侧 556 vs 病毒侧 181）。这是方法口径差异：eve_screen 的侧翼比较对"古降解 + TE 包裹"的位点保守，hi-fever/CAULIFINDER 的结构+属级域判据对这类敏感。
- 163-309 个扫描蛋白层完全不可见（CAULIFINDER 的 900bp placement 片段、LTR/gag-only 残迹等），属方法边界而非矛盾。
- **结论：两套方法正交互补**；合并清单 = eve_screen viral_supported（近期/完整插入）∪ 四象限位点（古降解位点）。

## 与病毒组 DNA contig（前一轮 500 条）的交叉

- 病毒组 62 条 EVE_strong（≥90% 匹配基因组）contig 的落点：**0 条**在扫描 viral_supported 位点；18 条落在 host_like 位点；44 条扫描未见 → 病毒组转录匹配到的全是**古降解内源拷贝**，进一步支持此前"392 条无同源候选大概率也是内源来源"的判断。
- 扫描的 262-493 个 viral_supported（近期完整插入）没有出现在病毒组 contig 里（转录本未进发现管线候选或被寄主去除步骤滤掉）。

## 类病毒层（新发现 · word 7 精查已完成 2026-09-29）

三个基因组都携带 **Citrus exocortis Yucatan viroid（CEVd-Y）样**序列（21/21 个 isolate ≥90% 覆盖）。word 7 精查（位点 ±10kb，`cevd_fine/`）结果：

| 基因组 | CEVd 样位点 | **全长 371nt 连续拷贝** | 最长单拷贝 |
|---|---|---|---|
| zhonghua | 10 | **0**（最高 84%，run 287nt） | — |
| ningxia | 14 | **3**（全部在 nSX00015） | 371/371 @98.7% |
| heiguo | 50 | **11**（Chr01/04/05×4/06/08 + 未挂载 hS00065/hS00104/hS00471） | 371/371 @93.5–98.7% |

### 结构判定：这是一个携带 CEVd-Y 样单体的重复家族，不是独立整合事件

- 16 个全长拷贝的接合处只有**两种模式**（同一元件的正反两个方向）：模式 A `CTCACGAAGGACTGCCAGTGATATACT | CEVd样371nt | CGTAACAAGGTA-GGTAGCCGTAGGGGAACCTGTGGCTGG`，模式 B 为 A 的反向互补。同一接合序列在两个物种、9+ 条染色体上逐碱基相同 → 侧翼属于重复元件本体，不是寄主靶位重复（无 TSD）。
- 宁夏 nSX00015 上三个全长拷贝与另 3 份侧翼拷贝呈 **~12.6kb 等距阵列**（127204/253672/432537 ± 同一上游 84nt 片段出现在 114573/241041/633490）；nSX00002 上还有一个 **1,305nt 的更长同源拷贝**（96.3%）。
- 外侧翼对 Dfam3.9 植物库零命中（1e-3）→ 不是已知植物 TE；元件本体（CEVd 核心 + 保守侧翼）的性质待 nt/Rfam 比对定名。

### 与 RNA-seq CEVd 检出的交叉（关键）

已知病毒侧 26 个 CEVd-Y 阳性样本的宿主归属：**L. ruthenicum 18 / L. barbarum 7 / L. chinense 1** —— 与基因组全长拷贝的分布（heiguo 11 全长 / ningxia 3 全长 / zhonghua 0）**完全吻合**；且检出载量极低（每样本 11–196 条 reads）。最简洁的解释是：**这 26 个"CEVd 检出"主要是内源 CEVd-Y 位点的转录**，而不是流行感染。
定案还需 read 级证据：把 CEVd 阳性样本的 reads 回比全长插入位点，看 (a) 371nt 是否全覆盖、(b) 有无跨接合处 read pair（证明转录自基因组位点）。

## 产物

- 服务器：`~/eve_screen_goji3_20260929/`（全部分层产物；`cevd_fine/` 为 word 7 精查）；本地本目录已回传 3 × eve_summary、kingdom_summary、family_by_genome、3 × viroid_word11、cevd_fine/
- 脚本：`scripts/audit/dnacheck_viroid_word11.sh`、`dnacheck_viroid_blastn_only.sh`、`dnacheck_cross_analysis.py`、`dnacheck_cevd_fine.py`
- 待办：全长 CEVd 位点的 read 级转录验证（回比 + 接合 read pair）；元件侧翼 nt/Rfam 定名；Phage 位点抽查；viral_supported 完整位点的 LTR 定年（接 EDTA 表）
