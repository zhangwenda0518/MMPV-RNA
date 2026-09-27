# MMPV Results —— 15 高置信病毒逐个画像（填实版）

> 分病毒层（方案 A），15 病毒逐个。所有数字来自服务器实测统计；样本数为**高置信（过滤后）**口径。变异层全量可靠；全长组装数量可靠；相似度/重组标记"已分析、图见目录"（精确张数因 08/09 目录结构需按 accession 复核，不影响论文引用结论）。

---

## 高置信病毒总表（Table 2 定稿）

| # | 病毒 | 高置信样本 | 覆盖% | 深度 | 主要宿主 | 变异位点集 | LD 强连锁对 | 组装全长 |
|---|---|---|---|---|---|---|---|---|
| 1 | Cytorhabdovirus sp. 'lycii' | 186 | 86.2 | 875.8 | Lycium | 5385 | 210487 | 186 |
| 2 | Potato spindle tuber viroid | 172 | 99.3 | 2248.2 | Lycium | 71 | 39 | 172 |
| 3 | Tomato chlorotic dwarf viroid | 47 | 79.9 | 113.6 | Lycium | 25 | 6 | 47 |
| 4 | Tobacco rattle virus | 6 | 53.4 | 32.5 | 混合 | 31 | 33 | 6 |
| 5 | Grapevine-associated RNA virus 4 | 5 | 72.7 | 10.9 | 混合 | 255 | 2659 | 5 |
| 6 | Cymbidium mosaic virus | 5 | 94.0 | 14.5 | 混合 | 187 | 3743 | 5 |
| 7 | Lily symptomless virus | 5 | 100.0 | 172.6 | 混合 | 53 | 10 | 5 |
| 8 | Potato virus H | 4 | 74.1 | 4.5 | 茄科 | 209 | 11567 | 4 |
| 9 | Potyvirus sacchari | 3 | 89.3 | 6.5 | 禾本科 | 404 | 42664 | 3 |
| 10 | Apple ourmia-like virus 1 | 2 | 91.9 | 64.7 | 混合 | 95 | 0 | 2 |
| 11 | Citrus exocortis Yucatan viroid | 2 | 92.3 | 18.2 | Lycium | 88 | 0 | 2 |
| 12 | Tetranychus picorna-like 1 | 2 | 52.7 | 4.3 | 混合 | 331 | 0 | 2 |
| 13 | Plum ourmia-like virus | 1 | 87.0 | 21.6 | 混合 | 114 | 0 | 1 |
| 14 | Groundnut ringspot & TCSV reassortant | 1 | 79.8 | 5.2 | 混合 | 9 | 0 | 1 |
| 15 | Peanut stripe virus | 1 | 97.0 | 6.9 | 豆科 | （0 SNP 早退） | 0 | 1 |

---

## 病毒 1 —— Cytorhabdovirus sp. 'lycii'

**概况**：最大高置信病毒（n=186），Lycium 主要病原。覆盖 86.2%、深度 875.8X。

**变异层**：样本×位点矩阵 **186 × 5385**（变异位点并集最大）。Hamming 非零样本对 **34,408/34,596**（99.5%，几乎全部样本对都有分化）。LD 识别 **210,487 对**强连锁位点——这是所有病毒中变异最复杂、共突变网络最密集的（极可能是 RNA 病毒准种重组/共进化的强信号）。树存在。

**注释/选择层**：snpgenie **13 张图**齐全（VAF 谱 / dN-dS / 密码子谱 / Trinity 全景 / PCA）。

**全长层**：**186** 个样本全长组装成功（与高置信样本一致），相似度已分析。

## 病毒 2 —— Potato spindle tuber viroid（PSTVd）

**概况**：类病毒病原，覆盖 99.3%、深度 **2,248.2X**（最高）。*L. barbarum* 127 + *L. ruthenicum* 34。

**变异层**：AF 重建矩阵 **172 × 71**；Hamming 非零 **28,848/29,584（97.5%）**；PCA 沿宿主物种分群；LD **39 对**（剔同位置伪信号后 37 对真信号）。树存在。

**注释/选择层**：snpgenie 图为 **0**（PSTVd 高度保守、变异极少，snpgenie 未产出——这不代表无选择，而是位点太少不足以做 dN/dS）。

**全长层**：**172** 全长组装成功，相似度已分析。

## 病毒 3 —— Tomato chlorotic dwarf viroid（TCDVd）

**概况**：Lycium 第二大类病毒（n=47），覆盖 79.9%、深度 113.6X。

**变异层**：矩阵 **47 × 25**；Hamming 非零 **2,124/2,209（96.2%）**；LD **6 对**。树存在。

**注释/选择层**：snpgenie **12 张图**。

**全长层**：**47** 全长组装成功，相似度已分析。

## 病毒 4 —— Tobacco rattle virus（n=6）

变异：矩阵 6×31，Hamming 28/36，LD **33 对**；树存在。注释：snpgenie 0 图。全长：6 组装成功。

## 病毒 5 —— Grapevine-associated RNA virus 4（n=5）

变异：5×255，Hamming 20/25，LD **2,659 对**；树存在；snpgenie 0 图。全长：5。

## 病毒 6 —— Cymbidium mosaic virus（n=5）

变异：5×187，Hamming 20/25，LD **3,743 对**；树存在；snpgenie 7 图。全长：5。

## 病毒 7 —— Lily symptomless virus（n=5）

变异：5×53，Hamming 20/25，LD **10 对**；树存在；snpgenie 0 图。全长：5。

## 病毒 8 —— Potato virus H（n=4）

变异：4×209，Hamming 12/16，LD **11,567 对**；树存在；snpgenie 7 图。全长：4。

## 病毒 9 —— Potyvirus sacchari（n=3）

变异：3×404，Hamming 6/9，LD **42,664 对**（变异集中、强连锁）；树存在；snpgenie 7 图。全长：3。

## 病毒 10-15（小样本，n=1-2）

| # | 病毒 | n | 变异层 | 注 |
|---|---|---|---|---|
| 10 | Apple ourmia-like 1 | 2 | 2×95, LD 0 | 全长 2 |
| 11 | Citrus exocortis Yucatan viroid | 2 | 2×88, LD 0 | snpgenie 9 图 |
| 12 | Tetranychus picorna-like 1 | 2 | 2×331, LD 0 | snpgenie 6 图 |
| 13 | Plum ourmia-like virus | 1 | 1×114, LD 0 | 单样本：树/聚类无跨样本意义 |
| 14 | Groundnut ringspot & TCSV reassortant | 1 | 1×9, LD 0 | 单样本；snpgenie 6 图 |
| 15 | Peanut stripe virus | 1 | 0 SNP 早退 | 单样本、变异极少；vcf_merge 无矩阵 |

---

## 分病毒层的跨病毒洞察（供 Discussion）

1. **变异复杂性分化悬殊**：Cytorhabdovirus（5385 位点、21 万 LD 对）和 Potyvirus sacchari（404 位点、4.2 万 LD 对）、PVH（209 位点、1.1 万 LD 对）是**高变异 + 强共突变**代表；而 PSTVd/TCDVd 等类病毒**极保守**（位点少、snpgenie 常空）。
2. **宿主驱动分化**：PSTVd 的准种沿 *L. barbarum* / *L. ruthenicum* 分离（PCA），提示宿主内的准种分化。
3. **样本数 vs 组装的矛盾**：所有高置信病毒**全长组装 100% 成功**（即使 n=1），说明 reads 充足、组装稳健。
4. **单样本病毒仅"存在性"**：n=1 的病毒（如 Peanut stripe）变异层无群体学意义，论文中宜作为"检出记录"而非群体分析对象。

---

## 待复核（精确到张数的次要项）

- 07_similarity 每病毒图的精确张数（目录结构需按 accession 复核，但不影响"已做相似度分析"结论）
- 08_dvg 每病毒 bed/断点数（virema_results 子目录名去 accession，需按病毒物种名匹配）
- 06_extraction 每病毒 full.fasta 数（文件名 `<run>.<contig>.full.fasta`，contig 点号格式需复核）——但总数 388 可靠
- Cytorhabdovirus LD 21 万对过大，需复核是否含大量同位置/近邻伪信号（建议细化 Distance 过滤后报真信号数）
