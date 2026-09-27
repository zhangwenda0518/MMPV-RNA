# 共识 + 逐级淘汰异议工具 规格（草案 v0.4）

时间：2026-08-30　｜　范围：`virus_classifier_analysis.R`（`build_consensus` :511）
状态：**只读预演完成，未改源码、未写产物**
预演对象：`goji-virome/02_novel_virus/RNA-Lycium_barbarum_out/05_Taxonomy/Votus.integrated/`

---

## 一、机制（大王 2026-08-30 定调）

不是拿分类谱系去过滤分类单元，**在工具之间做逐级淘汰**：

```
active = {全部工具}                          # 活跃工具集
for rank in [Realm, Kingdom, …, Species]:    # 从粗到细
    recs = active 中在该阶元报了有效值的工具
    if recs 为空:
        该阶元留空，active 不变               # NA 不是反对票
    elif recs 只有一个不同取值:
        采用它，不淘汰任何人                  # 只有 ACVirus 有分类就按 ACVirus
    else:
        共识 = recs 的加权票首
        报了别的值的工具 -> 移出 active        # 后面更细的阶元不再看它们
        该阶元 = 共识
```

三条硬边界（实测支撑，见第三节）：
1. **NA 不算异议**：工具在该阶元没值，不进淘汰名单。
2. **淘汰门槛用「票首占比 > 50%」**：无明确多数时不淘汰，直接取票首，避免误杀。
3. **被淘汰的工具不复活**：即便后续阶元活跃工具全体无值，也不回头取被淘汰工具的值。

另：删除「最完整工具行」整行补值兜底（`:537-546`，实测生效 0 格）。

## 二、实测结果（barbarum，21,630 contig，7 工具）

### 1. 科属/科种不匹配（VMR 谱系校验，仅算 VMR 内可查的行）

| 方案 | Species 侧科不一致 | Genus 侧科属不一致 |
|---|---|---|
| 现行（每阶元独立取加权票首） | 1,664 / 12,726（13.1%） | 1,597 / 14,402（11.1%） |
| **P1 逐级淘汰（保守门槛）** | **37 / 11,136（0.3%）** | **39 / 13,016（0.3%）** |
| P2 逐级淘汰（有分歧即淘汰） | 0 / 11,113（0.0%） | 2 / 12,994（0.0%） |
| P1b 逐级淘汰 + 允许复活 | 534 / 11,828（4.5%） | 1,162 / 14,343（8.1%） |

结论：**逐级淘汰本身够了，不需要再叠 VMR 谱系约束**。P2 只比 P1 多清掉 37 格，
代价是多留空 60 格，不值。P1b 证明复活会让不匹配回到 8%。

### 2. NA 原则保护的对象（只有 1 个工具报值的格数）

| 阶元 | 单工具独有 | 占有值格 |
|---|---|---|
| Realm | 7,786 | 40.9% |
| Family | 8,329 | 46.5% |
| Genus | 10,055 | **62.8%** |
| Species | 10,254 | **56.8%** |

过半的低阶元格子只有单一工具给值。若把"无值"当成"不一致"，这些格子会被整片清掉。

### 3. 淘汰事件分布（首次被淘汰的阶元）

| 阶元 | 主要被淘汰者 |
|---|---|
| Realm | VITAP 8,904、ACVirus 1,577、metabuli 832 |
| Order | ACVirus 2,555、metabuli 1,083、mmseqs 990 |
| Family | ACVirus 1,187、metabuli 490 |
| Genus | metabuli 935、ACVirus 446 |
| Species | metabuli 2,778、CAT 2,416、diamond_lca 1,481 |

VITAP 在 Realm 层就被大面积淘汰（与它在低阶元自举一致率 0.02 到 0.03 一致）。
大王举的例子真实存在：**门层 ACVirus 与多数不同 = 588 条 contig**。

### 4. 代价

| 方案 | Family 改变/变空 | Genus 改变/变空 | Species 改变/变空 | Order 改变/变空 |
|---|---|---|---|---|
| P1 | 892 / 350 | 1,835 / 1,429 | 2,035 / 837 | 778 / 254 |
| P2 | 940 / 362 | 1,866 / 1,451 | 2,067 / 851 | 835 / 264 |

淘汰后活跃工具数：剩 7 个 13,296 contig、6 个 4,553、5 个 2,597、4 个 1,003、3 个 155、2 个 26。
**没有出现只剩 1 个工具唱独角戏的情况。**

## 三、机制固有后果（需大王确认接受）

粗阶元一定音，整条链跟着走。实测样例：

```
CRR1126135_clean_NODE_138   Phylum  Nucleocytoviricota -> Uroviricota
                            Class   Megaviricetes      -> Caudoviricetes
                            Family  Mimiviridae        -> Herelleviridae
CRR1126135_clean_NODE_216   Order   Pantevenvirales    -> Pimascovirales
                            Family  Kyanoviridae       -> Marseilleviridae
                            Species Palaemonvirus pssm7 -> Marseillevirus LCMAC201
```

这些翻转是机制的性质而非缺陷：粗阶元的多数票赢了，细阶元服从。也不排除某些 contig
粗阶元本身判错。要判"某 contig 到底属哪边"只能看序列本体标志基因，需点名才跑。

## 四、顺带实锤的一个真 bug（占位值漏网）

`:622` 用 `tolower(get(col)) %in% PLACEHOLDER_TAXA` 做**整串精确匹配**，
而 `PLACEHOLDER_TAXA` 里只有 `"environmental"` 和 `"samples"` 两个独立整串，
多词形式全部漏网。盘上产物实查：

| 阶元 | 漏网取值 | 条数 |
|---|---|---|
| Species | `environmental samples` | 101 |
| Species | `Megaviridae environmental sample*` | 50 |
| Genus | `uncultured partitivirus` | 1 |

原始输入侧：mmseqs 168、metabuli 187、CAT 198、diamond_lca 194 条 Species 带
`environmental samples`（CAT 另有 163 条 `Megaviridae environmental sample*`）。

修法建议（待点名）：保留通用词的整串匹配，另加一条子串规则清掉
`environmental` / `uncultured` / `unclassified` / `unidentified` / `unassigned` 开头的多词串。

## 五、R 源码真实路径（更正）

| 文件 | md5 | 时间 | 备注 |
|---|---|---|---|
| `~/MMPV-RNA/virome_discovery_pipeline/virus_classifier_analysis.R` | 3e6c0cd8 | 2026-09-15 01:41 | **现行 v6.2（含闸门）** |
| `~/bin/virus_classifier_analysis.R` | b37f17b3 | | 另一份，内容不同 |

本地快照 `scripts/audit/snapshot/virus_classifier_analysis.R` 与现行一致（同 md5）。

## 六、待定

1. **落地方式**：改 R 源码 `build_consensus`（根治，05/09/10 需重跑 8 项目）
   还是 后置调和层（不动 R，只处理现有产物）。
2. 占位值 bug 是否一并修。
3. 第三节的整链翻转是否接受。

## 七、脚本与日志

| 脚本 | 日志 |
|---|---|
| `/tmp/rehearse_cascade.py` | `/tmp/cascade.log` |
| `/tmp/probe_placeholder.py` | 直接 stdout |
| `/tmp/rehearse_toptobottom.py` | `/tmp/t2b.log` |
| `/tmp/check_replica_fidelity.py` | `/tmp/fid.log` |

本地副本在 `MMPV-RNA/scripts/audit/`。基线保真度：复刻 vs 盘上 08-12 产物逐格一致率
Realm 0.990 / Family 0.964 / Species 0.906，差异集中在多值格且产物早于闸门，
改格数为量级参考。
