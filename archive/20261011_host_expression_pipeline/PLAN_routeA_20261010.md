# ⑧ host_expression_pipeline — route A 重规划（2026-10-10）

> 起因：枸杞（*Lycium barbarum*）route A 单网实测推翻了本模块 05_Coexp 现行的
> **「按病毒状态各建网 → 比 ΔzScore」** 设计。本文给出证据、新架构、文件级改动
> 清单与迁移验证方案。**旧设计不删除，降级为 `--coexp-mode dual` 的 legacy 路径。**

---

## 0. 一句话结论

| 项 | 旧（现行代码） | 新（route A） |
|---|---|---|
| 建网 | 阳性网 + 阴性网，各一张 | **全部样本一张网** |
| 条件注入 | 靠"两网之差" Δz | **靠 DEG**（网不动，事后过滤） |
| 主产物 | `diff_coexp_*.tsv`（Δz 显著对） | **DEG 限定共表达模块**（Louvain） |
| 分区 | 两网各分一次 k-means | 单网分一次；分区是**注释轴**，不是分组轴 |
| 分区秩 | `Get_Partition_Rankings.py`（有缺陷，且未 vendor） | **修正版 `partition_rankings.py`** |

---

## 1. 为什么要改（实测证据，全部可复现）

### 1.1 双网络 Δz 会把「样本量差」误读成「生物学差」

`zScore(Co-exp_Str_MR)` 是**网内全局** z 标准化（不是逐基因）。两张网的样本集不同
⇒ 两个 z 的基准不可比，Δz 被分布截断主导。

- 枸杞实测：`routeA_all410` 单网 max z = **1.81068**、min = **−2.0484**，
  z≥1.70 的基因对只占 **2.02%**（分布右尾很薄）。
- 分区基线共表达均值与「哪一分区拿到 argmax」的 Spearman = **+0.988**；
  样本数与 fold 的 Spearman = **−0.906** ⇒ **样本数几乎单定了结果**。

### 1.2 P6 专项判据被证伪（最直接的证据）

`ningxia00107--ningxia01221` 一类边在 P6 上的富集：

| 口径 | P6 fold | 判读 |
|---|---|---|
| 全量 136,331 边（照官方脚本字面读） | 1.1427（z=+17.6） | "P6 最强" |
| complete-case（10 个分区全有值的 92,672 边） | **0.9858（z=−1.44）** | **随机** |

机制：P5（16 样本、6.2% 阳性、29.5% nan）是最强霸榜者（fold 4.65），P6 只是在
P5 缺失的边上**补位**。照字面读全量会得出「**病毒阴性分区 P5 驱动最强共表达**」的荒谬结论。

分区内 z 标准化后排序直接反转：P4（113 样本 / 66.4% 阳）1.453 变第一，
P6 掉到 0.816（第 8）。

### 1.3 官方 `Get_Partition_Rankings.py:136-139` 的秩公式有两处缺陷

```python
scores = np.array([float(score) if score != "nan" else 0
                   for score in line_contents[1].split(",")])      # ① nan→字面量 0
ranks = stats.rankdata(np.array(scores), method="min", nan_policy="omit")
ranks_centered = ranks - int(np.nanmedian(ranks))                  # ② 截断取整
```

① **nan 被替换成字面量 0 参与排序**（不是剔除）⇒ m 个 nan 并列拿最低秩，
有效分区的秩被整体抬高 m 档。
实测 `ningxia00107--ningxia01221`：只有 P2 有效，却得秩 `0,0,9,0,...`；
`ningxia00057--ningxia26119`：4 个有效值得 `6,7,8,9`。

② `int(nanmedian)` 在 nan 占多数时恒为 1 ⇒ 中心化常数不稳定。

逐行核实：complete-case 13,590/13,592 一致、nan 行 6,408/6,408 一致
（即差异**全部**来自这两处公式缺陷，不是数据问题）。

### 1.4 分区主要分的是**组织**，不是病毒状态

`Partition_annotation/partition_annotation.tsv`（超几何 + BH，官方
`Annotate_partitions.py` 同思路）：

```
P0  (tissue:leaf of wolfberry : 3.7e-06)
P1  (tissue:fruit : 1.09e-07)
P2  —
P3  (tissue:leaf : 2.57e-49)
P4  (tissue:fruit : 2.08e-41)
P5  (tissue:fruit : 0.000215)
P6  (GCVA+ : 0.00105),(Virus_any+ : 0.00156),(PSTVd+ : 0.00178)
P7  (tissue:flower : 2.69e-25),(tissue:anther : 2.8e-06),(PSTVd+ : 0.0178)
P8  (tissue:fruit : 0.000439)
P9  (tissue:root : 1.09e-07)
```

**只有 P6 和 P7 有病毒标签富集，且 p 值（1e-3 量级）比组织 p 值（1e-41、1e-49）
弱 10 个数量级以上。** ⇒ 「P3 = 阴性分区」这种说法不成立 —— P3 是**叶簇**
（51.2% 阳性，并不显著阴性）。此前按阳性率贴的标签（阴性/阳性/PSTVd）只能当
**描述**，不能当**分区身份**。

### 1.5 k 选择不稳健

`silhouette_coefficients_across_k.tsv`：

```
k=5  0.3047   k=6  0.3154   k=7  0.3147   k=8  0.2812   k=9  0.2829
k=10 0.3160 ← 被选中     k=11 0.3026   k=12 0.2672  ...  k=20 0.3003
```

k=10（0.31597）与 k=6（0.31539）相差 **0.0006** —— 在 k=5…30 全区间
（0.235–0.316）里几乎是平的。**k=10 是"平顶上的任意一个点"**，不是明确的最优。
分区数不稳定 ⇒ 所有"分区 X 富集了什么"的结论都要带这个不确定性。

### 1.6 富集口径需要固化（见 §3）

### 1.7 缺失组件

| 缺什么 | 现状 | 影响 |
|---|---|---|
| 分区注释 | 本模块**没有**（goji 侧有 `annotate_partitions_routeA.py`） | 分区无法解读 |
| 分区秩 | 官方脚本**未 vendor**，且公式有缺陷 | 分区轴整条断链 |
| GO/KEGG 富集 | README §7 明写「**未纳入**」（需 eggNOG/Mercator） | 模块/分区无功能解读 |
| 分区签名基因集 | 无（分区是样本簇，直接富集无对象） | 见 §3 |

### 1.8 ★ routeA 的 410 样本集里混了 19 个近空文库样本（2026-10 新发现）

**这是比 1.1–1.5 更靠前的根因 —— 它污染的不只是分区，是全部下游。**

**事实**

- routeA 用的是 `05_expression/quant/RNA-Lycium_barbarum/matrix/TPM_matrix.tsv`
  （`routeA_all410/inputs/TPM_matrix.tsv` 是它的**符号链接**）
- 该矩阵 410 列，而 salmon 一共定量了 445 个样本 ⇒ 有一步 445→410 的过滤
- **该过滤只是「去掉 lib == 0」的零值判断**，不是量级判断：
  被去掉的 35 个样本恰好全部是 `lib == 0`；剩下 19 个「有值但极小」的样本原样进入
- 这 19 个：`library_size` 1,065 ~ 15,860（中位正常样本 10.8M），
  **全部 SE（R2 缺失按单端跑）、全部 `Virus_any = 阴性`**
  - 9 个 `CRR1295698–706`：`metadata/missing_read_samples.tsv` 记录 R2 缺失
  - 10 个 `SRR2163069–72` / `SRR36147586–91`

**后果（实测）**

| 指标 | 含这 19 个 | 剔除后 |
|---|---|---|
| size factor 范围 | 0.001 ~ 51.8 | 0.178 ~ 3.735 |
| leaf 内 M6（光合模块）log2FC | **−0.933** | **−0.158** |
| leaf 内显著 M6 中上调比例 | 13% | 28% |
| 全队列显著 M6 数 | 71 | 320 |

⇒ 污染**制造了一个假的强下调信号**。所有基于该矩阵的结论（共表达网、Louvain 模块、
k-means 分区、富集）都需要在剔除后重跑。

**门槛取值的实测依据**：`lib < 1e5` 的 54 个样本与 `lib ≥ 2.4e6` 的正常样本之间
存在 **约 100 倍空档**（最差正常样本 2,438,517 vs 最好问题样本 15,860）
⇒ 阈值在 [1e5, 2.4e6] 内取任何值，结果完全一致。

**修法**：把闸门前移到 salmon 阶段（§2b），让 `matrix/TPM_matrix.tsv` 出厂即干净。

**两级串联实测：`445 → 391 → 378`**

| 阶段 | 通过 | 剔除 | 原因 |
|---|---|---|---|
| 起点 | 445 | — | — |
| Gate 1（02_Quant 源头） | 391 | 54 | 35 个 `lib=0` + 19 个近空文库 |
| Gate 2（03_QC，**mapped ≥ 30%**） | **378** | 13 | 全部因 `mapped < 30%`；文库量 2.7M–5.3M、检出 ~24–26k，属「文库没问题但比对率低」 |

**阈值决定（2026-10-10，用户裁定）**：Gate 2 的 `--min-mapped` 默认由 **50 → 30**。
理由：枸杞公共多组织队列比对率中位仅 60.2%，按 50% 会多剔 50 个（掉到 341）；
按 30% 只剔 13 个。且 30% 与 Gate 1 的 WARN 线对齐（源头标出、策略剔除）。

**在 routeA 的 410 样本集里**：Gate 1 剔 19 → **391**；再叠 Gate 2(30%) 剔 13 → **378**。

---

## 2. 阳性/阴性到底怎么来的（回答"分区后怎么划分阳/阴"）

**关键：阳性/阴性不是从分区推出来的，两者是独立的两条线。**

```
④ 病毒检出表 ──derive_virus_labels.py──▶ 样本级阳性/阴性 ─┐
   (Covered_pct ≥ 90%)                                    │
                                                          ├─▶ 06_DE 分组
② clean reads ─salmon─▶ TPM 矩阵 ─Generate_partitions─▶ 分区 ┤
                              (k-means, 与病毒标签无关)      │
                                                          └─▶ 分区注释
                                                              (超几何: 组织 + 病毒)
```

1. **样本级阳性/阴性**：`derive_virus_labels.py` 读 ④ 检出表，
   逐样本逐病毒 `Rep_Coverage(%) ≥ --min-coverage`（默认 **90%**，与 goji 口径一致）
   判阳；`Virus_any` = 任一病毒阳性。同时保留连续载量列（`<病毒>_load`）供剂量扫描。
2. **分区**：`Generate_partitions.py` 对**标准化后的表达矩阵**做 PCA → k-means，
   k 由轮廓系数选。**整个过程不看病毒标签** —— 是无监督的样本分群。
3. **分区标签（阴性/阳性/PSTVd）**：事后由 `Annotate_partitions` 式超几何检验
   （分区 × 标签取值，背景 = 全部样本，BH 校正）得到。这才是"分区代表什么条件"的
   正式答案 —— 而 §1.4 显示它主要给的是**组织**。

> **对本模块的直接后果**：分区轴必须与分组轴**分开**。不能拿"分区 = 阳性"当分组
> 去建网（那正是双网 Δz 的错），也不能拿分区注释去替代样本级病毒标签。

---

## 3. 富集分析怎么做（回答"怎么进行富集"）

分区是**样本簇**，不是基因集 ⇒ 必须先把分区转成基因集。三条口径，按可靠性排序：

### 口径 A（推荐）：分区表达特征基因

```
X = log2(TPM + 1)                      # 基因 × 样本
分区 P:  A = X[:, inP]   B = X[:, ~inP]
        Welch t = (mA − mB) / sqrt(vA/nA + vB/nB)
        基因集 = |t| 最大的 top-N（N = 500）
```

优点：分区大小悬殊时比 fold 稳；不依赖网络。枸杞实测 8 个分区各 500 基因，
两两 Jaccard 全 <0.1（分区确实是 8 套独立表达程序）。

> **反例（不要用）**：按"分区内共表达最强的基因"定义基因集 —— 枸杞实测
> 8 个分区候选基因**全部为 0**，因为 top-5 边表每基因平均只有 ~7 条边，
> 「≥100 条有效边」的门槛根本够不着。

### 口径 B：分区内 z 标准化后的分区秩（修正版，替代官方脚本）

```python
# 修正 Get_Partition_Rankings.py 的两处缺陷
valid = ~np.isnan(scores)                      # ① 剔除, 不是置 0
r = stats.rankdata(scores[valid], method="min")
center = np.median(r)                          # ② 不取整
centered = r - center
# 未参与的分区输出空（不是 0，也不是最低秩）
```

### 口径 C：模块 / DEG（KingdomStress 主线）

`z ≥ 1.7 ∩ 两端均为 DEG` 的边 → Louvain(res=1.0) → 模块基因集。

### 统计与背景（三口径共用）

| 项 | 取值 | 理由 |
|---|---|---|
| 检验 | 超几何 ORA，`P(X≥k)=hypergeom.sf(k-1, N_BG, N_Term_BG, N_Set)` | 与 `clusterProfiler::enricher` 同口径 |
| N_BG | **背景 ∩ 注释**（可注释空间） | 无注释基因不可能命中，不该进抽样空间 |
| 检验族 | **背景内全部过规模门槛的 term**（不是"被命中的 term"） | 否则 BH 族由结果决定 → 假阳膨胀 |
| 条目过滤 | **5 ≤ term 背景 ≤ 300**（全转录组背景） | 防超级条目（"代谢过程"）与超小条目 |
| 校正 | BH-FDR，输出同时给 `N_BG` 与 `N_BG_Raw` | 跨工具可核对 |
| **背景选择** | **全转录组**（推荐）；DEG 背景仅作敏感性 | 枸杞实测：M6 光合 fold 15–30（全转录组）vs 5–6（DEG 背景），DEG 背景把信号稀释 |

### 解读纪律（用户裁定，必须写进产物）

1. **通路（条目）可以跨分区重复，基因不会** ⇒ 同一条目在多个分区出现，
   **不能当成"多次独立验证"**加权重。
2. **必须看基因级重叠，不看条目名**。枸杞实测：P3 与 P6 都富集
   `photosynthesis`，但 P3 命中 75 基因 / P6 命中 24 基因 → **重叠 0**；
   `ko00195 Photosynthesis` P3 43 / P6 45 → **重叠 0**。
3. 48 条跨分区共享条目中 **30 条（62.5%）基因重叠为 0**；只有
   P3&P8（花器官/RNA pol II 转录调控，最多 7 共同基因）与
   P7&P8（微管—细胞壁，最多 3）是真共用基因。
4. **注意：基因交集为 0 不是"分区"自动保证的** —— 同一基因在一个分区上调、
   在另一分区下调时，两次 |t| 都可能进 top-N。枸杞实测 4,000 个签名位
   → 去重 3,641 基因，359 个跨分区重复（`ningxia00622`→P4/P7/P8、
   `ningxia00703`→P3/P7/P8）。
5. **组织混杂必须声明**：分区组织构成会进入签名（枸杞 P3/P6 的光合信号
   高度可疑与叶样本比例相关）。

---

## 4. 新架构（route A）

```
02_Quant TPM ──┐
03_QC 名单 ────┼─▶ [matrices]  单网输入: all_samples.tsv + sample_manifest.tsv
               │                (旧: positive.tsv + negative.tsv  ← legacy)
               ▼
   [partitions] Generate_partitions ×1  (k 扫 + silhouette + 稳定性重采样)
               ▼
   [networks]   Run_TEA-GCN_mp ×1       (全部样本一张网)
               ▼
   [ranks]      Rank_transform ×1       → zScore(Co-exp_Str_MR)
               ▼
   [qc]         coexp_network_qc ×1     (max_z ≥ 1.7 ? 否则子网必空, Fig S8)
               ▼
   ┌───────────┴────────────────────────────┬──────────────────────┐
   ▼                                        ▼                      ▼
[modules]  DEG 注入                     [pannot] 分区注释        [prank] 分区秩
   z≥1.7 ∩ DEG×DEG 边                    分区 × (组织+病毒)        修正版秩公式
   → Louvain res=1.0                      超几何 + BH              (nan 剔除)
   → 模块表 / hub                                                       │
   ▼                                                                     ▼
[enrich]  模块 / DEG / 分区签名 基因集 ──▶ 超几何 + BH ──▶ 条目 + 基因级重叠核对
```

**条件注入的正确姿势**（KingdomStress 三步法）：一张网（全部样本）→ 按条件取
DEG → 在网里过滤到 `z ≥ 1.7 且两端均 DEG` 的边 → Louvain 得模块。
**从不按条件各建网、不比 Δz。**

---

## 5. 文件级改动清单

| 文件 | 改动 | 优先级 |
|---|---|---|
| `build_coexp_matrices.py` | 默认改为**不切分**，输出 `all_samples.tsv` + `sample_manifest.tsv`；切分逻辑保留在 `--coexp-mode dual` | P0 |
| `coexp_modules.py` | `--gcn-dir` 从「组名=目录（可重复）」改为**单目录**；`--group` 参数删除（模块只有一个来源网） | P0 |
| `coexp_network_qc.py` | 单网；闸门语义从「Δz 前体检」改为「**子网可行性体检**」：`max_z < 1.7` ⇒ 模块阶段必空，直接 skip | P0 |
| **新增** `annotate_partitions.py` | 分区 ×（组织 + 病毒标签）超几何 + BH；输出 `partition_annotation.tsv` / `partition_composition.tsv` / `partition_sample_assignments.tsv` | P0 |
| **新增** `partition_rankings.py` | 修正版官方 `Get_Partition_Rankings.py`（nan 剔除、不取整、未参与分区输出空） | P0 |
| **新增** `partition_signature.py` | 分区表达特征基因（Welch t top-N）→ 供富集消费 | P1 |
| `enrich_analysis.py` | ① 新增分区签名基因集自动收集；② 条目过滤默认 `5 ≤ bg ≤ 300`；③ 新增**基因级重叠核对**输出（跨集合共享 term 的基因交集） | P1 |
| `host_expression.py` | `COEXP_SUBSTEPS` 改为 `("matrices","partitions","networks","ranks","qc","pannot","prank","modules")`；新增 `--coexp-mode {routeA,dual}`（默认 routeA）；dry-run 文案同步 | P0 |
| `host_expression_report.py` | 新增分区注释段、分区秩段、基因级重叠段 | P2 |
| `tests/test_host_expression.py` | 新增：单网切分、分区注释超几何已知值、修正秩公式（nan 行）、签名基因集、重叠核对 | P0 |
| `README.md` | §3 标注 route A 为默认；新增 §3b 分区轴、§2b 定量 QC 闸门、§9b 富集口径；§7 对应表补缺失组件 | P0 |

**§1.8 的修法（已落地）**

| 文件 | 改动 | 优先级 |
|---|---|---|
| `salmon_auto_pipeline.py` | **新增定量后 QC 闸门**（默认开启）：`collect_quant_metrics` / `evaluate_quant_qc` / `write_quant_qc_report`；`merge_quant_results` 支持 `sample_filter`，主矩阵只含 PASS+WARN，全部样本另存 `*_matrix_all_samples.tsv`；新增 CLI `--min-library-size` / `--warn-library-size` / `--min-detected-genes` / `--warn-detected-genes` / `--min-mapping-rate` / `--warn-mapping-rate` / `--no-quant-qc` / `--qc-report-only` / **`--qc-only`**（从 quant 目录自举，给历史结果补做质控） | **P0（已实现）** |
| 产物 | `02_Quant/metadata/{quant_qc.tsv, qc_failed_samples.tsv, qc_pass_samples.txt}`；`run_summary.json` 增 `quant_qc_pass/fail/warn`、`matrix_samples`、`quant_qc_thresholds` | P0（已实现） |
| routeA 重跑 | 从 410 样本集剔除这 19 个 → **391 样本**，重跑网/模块/分区/富集；若同时启用 Gate 2（`mapped ≥ 30%`）→ **378 样本** | **P0（待执行）** |

**磁盘**：单网 RAW+Completed ≈ `70B × 基因数²`（39k 基因 ≈ **100 GB**，旧双网
200 GB）。单网后 `--coexp-workers` 与 `rank_transform_mp` 内存需求不变。

---

## 6. 迁移与验证

### 6.1 向后兼容
- `--coexp-mode dual`（默认关）保留旧全链路，产物目录不变，老项目可复现。
- `05_Coexp/` 产物**新增**而非覆盖：`partitions_single/`、`Partition_annotation/`、
  `Partition_rankings/`；旧 `matrices/{positive,negative}.tsv` 仅在 dual 模式生成。

### 6.2 验证（在枸杞 routeA 上对拍）
| 检查 | 期望 |
|---|---|
| 单网 z 分布 | `max_z ≈ 1.81068`、`min ≈ −2.0484`、z≥1.7 占比 **≈ 2.02%** |
| 子网 | DEG=10,641；`z≥1.7 ∩ DEG×DEG` ≈ **1,705,597 无向边 / 10,641 节点**（100% 覆盖，无孤立点，平均度 ≈320） |
| Louvain | Q ≈ **0.5316**（res=1.0, seed=42） |
| 模块数 | **7** 个（枸杞实测），方向极纯（5/7 单向，全局 DEG 上调仅 39.0%） |
| 分区注释 | P3=leaf 2.6e-49、P4=fruit 2.1e-41、P6=100% 阳性 1e-3 量级 |
| 修正秩公式 | complete-case 13,590/13,592、nan 行 6,408/6,408 与旧实现一致；**差异只出现在混合行** |
| 富集 | M6 `photosynthesis` fold ≈ 14.97、BH ≈ 1.1e-79（全转录组背景） |

### 6.3 尚未做（需用户确认后才做）
- **置换零模型**：打乱 DEG 标签 100–1000 次重跑 Louvain，检验 Q=0.53 与方向纯度
  是否超随机。这是审稿人一定会问的一条。
- **跨条件**：GCVA+ / PSTVd+ 各跑一遍模块，用 OC（`|OG₁∩OG₂|/min`）比三套模块。
- **模块-条目热图 / 网络可视化**（ggNetView 未纳入，需另接）。

---

## 7. 一页速查

- **阳性/阴性**：④ 检出表 → `derive_virus_labels.py`（Covered ≥ 90%）。**与分区无关**。
- **分区**：`Generate_partitions.py` 无监督 k-means，k 由 silhouette 选（枸杞 k=10，但曲线是平顶）。
- **分区身份**：事后超几何注释 —— 枸杞主要是**组织**（P3 叶 / P4 果实 / P7 花），
  只有 P6、P7 有病毒富集。
- **富集**：分区 → Welch t 签名基因集（top-500）→ 超几何 + BH，背景 = 全转录组，
  条目过滤 5≤bg≤300。
- **解读**：条目可重复，基因不重复；**看基因级重叠，不看条目名**。
- **建网**：全部样本**一张网**，条件靠 DEG 注入。**不比 Δz。**
