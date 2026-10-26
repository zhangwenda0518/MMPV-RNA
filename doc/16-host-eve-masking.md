# 16 · 宿主参考基因组的「内源病毒（EVE）屏蔽」：要不要做、做到哪一层

> 状态：决策文档（2026-10-21）；§6 的 `eve-mask` 阶段已落地（2026-10-21，默认关闭）
> 触发问题：`public_metadata_pipeline` / `data_preprocessing_pipeline` 是否需要补充内源病毒屏蔽。
> 配套工具：`public_metadata_pipeline/mask_host_genome.py`（本文件 §5 使用，含单测）、
> `public_metadata_pipeline/build_host_pipeline.py` 的 `eve-mask` 阶段（§6，含单测）
> 参考：EuPathDB 去污染（PLOS Comput Biol 10.1371/journal.pcbi.1006277）、
> GENEX（GigaScience 10.1093/gigascience/giaf108）、mask-reference-genome（GitHub megjohnson1999）

---

## 0. 结论（先看这一节）

**需要补的是「掩码能力 + 量化诊断」，不是「默认开启掩码」。**

| 管线 | 要不要改 | 改什么 | 默认行为 |
|---|---|---|---|
| `public_metadata_pipeline` | **要（已落地）** | 新增可选阶段 `eve-mask`（宿主基因组 EVE 判定表/BED → 掩码 FASTA + 掩码率报告 + 建库溯源 `db_provenance.json`），建掩码库到 `hostdb_evemask/`；P0 诊断仍用归档 `eve_screen.py` | 保持现状：**不掩码**（不给掩码源时 `--stage all` 跳过它，建库路径/产物与 v3.1 一致）|
| `data_preprocessing_pipeline` | **代码不改** | 掩码发生在宿主 FASTA 层，Kraken2/比对库同源重建即可；只需在文档与 `summary` 里记录库版本（见 `doc/03-host_depletion.md` 末节）| 保持现状 |

三条理由（按重要性）：

1. **主力信号不受影响，收益低。** 本研究 15 个高置信病毒里，`Cytorhabdovirus`（RNA）、PSTVd/TCDVd（类病毒，不整合）都不依赖 DNA 型 EPRV 的屏蔽收益；屏蔽只对「与宿主内源 DNA 病毒元件同源的外源病毒」有意义，而这类病毒在本数据集的判定表里落点是「新种候选 / 远缘 / 宿主污染 / 未定」四类（L. barbarum 284 条中 24+77+40+125），不是核心结论。
2. **屏蔽有明确的反向代价。** 掩码区不再参与建库 → 该区**所有**宿主 reads 都变成「非宿主」进入组装，假阳性、嵌合、低覆盖 contig 都会上升（GENEX 自己也只把 masking 当整条流程的一环，而不是单独解法）。没有量化前就改生产库，是用已知的假阴换未知的假阳。
3. **改生产库不可逆且昂贵。** 1,348 样本已按现库跑完；换库 = 全流程重跑 + 与已发表口径不一致。正确顺序是「先量化损失 → 再决定是否换库」，而不是先换库。

一句话：**屏蔽应当是「有证据支撑的可选臂」，不是「默认步骤」。**

---

## 1. 问题定义：两个方向，别混在一起

宿主参考库与病毒参考库都可能「夹带」对方的序列，方向相反、后果相反：

| 方向 | 库中夹带 | reads 后果 | 统计后果 | 本研究对应问题 |
|---|---|---|---|---|
| **A** | 宿主库里含病毒样序列（EVE/EPRV、污染片段） | 外源病毒 reads 被当宿主删除 | **假阴性**，近缘病毒检测灵敏度下降 | 本文件主题 |
| **B** | 病毒参考库里含宿主/污染序列 | 宿主 reads 被当成病毒 | **假阳性**，丰度虚高 | `database/final.cluster.ref.fasta`（15,202 条）的卫生问题 |

三篇材料分别落在哪：

- **EuPathDB（PLOS pcbi.1006277）**：把「基因组草图→清洗」做成多步流程（100 bp 重叠伪读段 → Kraken 对污染库 → Bowtie2 对人基因组 → Dustmasker 低复杂度）。它解决的是「参考序列本身脏」，方向是 **B 的治法**（清洗参考库），可迁移到「清洗病毒参考库」，**不能**直接回答「宿主 EVE 要不要掩码」——EVE 不是污染，是宿主基因组的真实组成部分。
- **mask-reference-genome**：把病毒库切 60–80 bp 伪片段比对参考基因组、屏蔽「病毒样区域」。与本文件的 A 方向**完全同构**，可直接用（见 §4 的三个前提）。
- **GENEX（giaf108）**：把 GTDB 微生物基因组切 60 bp 伪读长比对真核参考基因组，屏蔽「微生物样区域」，目的是防止细菌/古菌 reads 被真核参考抢走——方向是「eukaryote 参考抢 prokaryote reads」，与我们「宿主参考抢病毒 reads」机理相同、对象不同。它的价值在于**证明这类参考库污染是系统性、可量化的**（最高 70%，非哺乳脊椎动物 <0.2%），并提供可复用的伪读长流程。

---

## 2. 机理：现有管线里，EVE 造成假阴的**两个**入口

去宿主主流程在 `data_preprocessing_pipeline/host_depletion.py:652-991`（`run_single`）：

```
[Stage 2] Kraken2 初筛 ──classified──> 展开子 taxID，seqtk 抽 keep_taxids(10239) ──┐
            │                                                                        ├─> 合并池
            └──unclassified──────────────────────────────────────────────────────────┘
[Stage 4] 比对精筛: bowtie2/hisat2/minimap2 → SAM.view_filter(mfilter=True) → 抽回 fastq
[Stage 5] Name-Sort        [Stage 7] rRNA 去除        [Stage 8] gzip
```

两个入口：

**入口 1（Stage 4，全量删除）**：`SAM.view_filter`（`host_depletion.py:491-502`）在 PE 下用
`-f 12 -F 256`，只判「双端都没比对上」= 非宿主；**只要比对到宿主基因组任何位置，整对 reads 删除**，
不看落在哪、不看是否与病毒同源。EVE 区 → 外源病毒 reads 一律判宿主。

**入口 2（Stage 2，救援被覆盖）**：`--keep-taxids 10239` 的救援（`:702-807`）发生在 Kraken2 之后、
比对之前；救援出来的 reads 随后仍进 Stage 4 的同一个合并池，**照样被比对删除**。
稿件第 52 行「preserving viral taxid 10239」只在 Stage 2 成立——这是一个 Reviewer 可以直接问的点。

**入口 3（隐含，Stage 2 自身）**：Kraken2 的 `--confidence` 越低越激进（更少 k-mer 支持就判宿主）。
`host_depletion.py:1017` 默认 `0.4`，而稿件 Methods 写的是 `0.2`。若实际跑的是 0.2，
那么与 EVE 完全一致的 reads 在 Stage 2 就已经被删掉，连 Stage 4 都到不了。
**待核**：从历史日志 `[Stage 2] Kraken2 初筛 (Conf: ...)`（`:706`/`:809`）或 `summary.tsv` 确认实际值。

> 结论：EVE 假阴不是「一个 bug」，是「库的属性」。屏蔽（N 化）是唯一能**同时关掉三个入口**的做法——
> 因为掩码后该位置不再产生 k-mer、也不再产生比对种子。

> **版本前提（2026-10-22 补）**：以上三个入口以本仓 `HEAD` 版 `host_depletion.py` 为准。服务器上那份
> 未入库的 7 月版在 `--keep-taxids` 非空时**无条件跳过 align**（默认 `10239` → 默认只跑 Kraken2 一层），
> 入口 1/2 不成立，只剩入口 3。两版差异、确认命令与部署建议见
> [`doc/17-server-vs-local-audit-20261022.md`](17-server-vs-local-audit-20261022.md) §3。

---

## 3. 在本数据集上的影响面

| 事实 | 数值 / 位置 | 对「要不要屏蔽」的含义 |
|---|---|---|
| 样本规模 | 1,348 宏转录组（L. barbarum 为主） | 换库成本高，不能默认改 |
| 主力检出 | `Cytorhabdovirus` sp. 'lycii' 17.0%、PSTVd 16.4% | 均为 RNA 病毒/类病毒 → DNA 型 EPRV 屏蔽无收益 |
| 高置信病毒 | 15 个 / 285 样本 / 442 关联 | 结论层不需要屏蔽才能成立 |
| 受影响判定 | L. barbarum 284 条 = 24 新种候选 + 38 已知 + 77 远缘 + 40 宿主污染 + 125 未定 | 屏蔽/证据层影响的是这四类的边界 |
| 稿件口径 | Limitations 已声明 "RNA-only mining misses DNA viruses" | 屏蔽不改变该限制；但「宿主 EVE 会吃掉近缘病毒 reads」目前**没有**写进去 |
| 枸杞 EVE 现状 | 仓库内无 Lycium 基因组的 EVE 鉴定记录（`archive/20261026_artifacts/outputs/MMPV_EVE_ML_20261020/`、`reports/` 均无 Lycium 字样） | **必须先扫**，不能凭「茄科植物普遍有 EPRV」直接掩码 |
| 宿主库组成 | RefSeq 基因组 + GFF3 + 细胞器 + Kraken2 污染库（archaea/bacteria/plasmid/fungi/protozoa/UniVec） | 细胞器序列（叶绿体/线粒体）里也可能有病毒样片段，掩码范围要显式决定 |

**一个必须单独排查的特例：PSTVd（16.4%）。** 类病毒不整合，所以它不该出现在宿主库里；
但如果历史上有人把 viroid 库加进 `--add-library`/`--k2-libs`，或 RefSeq 宿主组装里混入了类病毒载体片段，
那 16.4% 的主力信号会被「合法地」当成宿主删掉。**这属于方向 B 的污染排查，优先级高于 EVE 屏蔽。**

---

## 4. 三份材料怎么用 / 哪些不能照搬

### 4.1 mask-reference-genome（最贴近，但有三个前提，都必须先验伪）

用户给的命令（伪读长层）：

```bash
python ~/bin/reference_masker.py \
  --viral-database ~/database/virus-db/RVDB-v31/RVDB_viroids.fasta \
  --reference-genomes ningxia.genome.fasta heiguo.genome.fasta zhonghua.genome.fasta \
  --output-dir mask-genome --threads 60
cat mask-genome/*masked.fasta > mask-genome/all.masked.fasta
```

**前提 1：`RVDB_viroids.fasta` 到底是全 RVDB 还是「只有类病毒」？**

- 若是 **viroid-only**：类病毒不整合 → 拿它去掩码宿主基因组，几乎什么都掩不到，这次运行在科学上是**空的**。
- 若是 **RVDB 核苷酸全库 + viroids 合并**：才有 Caulimoviridae / Geminiviridae 等植物 DNA 病毒，才有意义。
- 判定（两条，5 分钟）：

```bash
seqkit stats -a ~/database/virus-db/RVDB-v31/RVDB_viroids.fasta
# 全库: 百万条量级 / GB 级;  viroid-only: 数千条 / 数 Mb —— 差 2-3 个数量级
```

再用 2–3 条 Caulimoviridae/Badnavirus 参考做探针 blastn，确认能命中：

```bash
makeblastdb -in RVDB_viroids.fasta -dbtype nucl -out /tmp/rvdb_nt
blastn -query database/<3 条 Caulimoviridae 参考>.fa -db /tmp/rvdb_nt \
       -outfmt '6 qseqid sseqid pident length' | head
```

**前提 2：90% identity 的核酸层会漏掉老 EVE，而漏掉的 EVE 仍然吃 reads。**
`--min-identity 90 --min-length 60 --max-indel 2` 意味着 60–80 bp 上的同源度要 ≥90%。
但 bowtie2 `--end-to-end` 默认打分（`--score-min L,-0.6,-0.6`、mismatch 6）对 100 bp reads
允许约 10 个错配 ≈ **90% 一致性**；也就是说 80–90% 一致性的 EVE 可能躲过伪读长检测，
却仍然能把外源病毒 reads 拉走。→ 掩码必须与**蛋白层**（`eve_screen.py` 的 DIAMOND 流程）取并集，
不能只跑核酸层。

**前提 3：该脚本没有「病毒库自身污染」的排除逻辑。** 病毒库里的宿主片段、载体片段、低复杂度区
也会被切成伪片段并命中宿主基因组，掩出来的不全是 EVE。→ 必须看 report，按「命中片段来自哪些病毒分类单元」
人工过一遍（这正是 EuPathDB 流程在做的事）。

### 4.2 EuPathDB：用来洗**病毒参考库**（方向 B），不是用来洗宿主

它的流程（100 bp 伪读段 → Kraken 污染库 → Bowtie2 人基因组 → Dustmasker）可以直接改造为：
把 `database/final.cluster.ref.fasta`（15,202 条）比对到 3 个枸杞基因组 + Kraken2 污染库，
>`覆盖度 ≥50% 且一致性 ≥95%` 的参考序列列为「宿主污染候选」。这条与 EVE 屏蔽同批做，成本低、收益直接
（对应「Toxoplasma 被误判为几乎所有角膜样本主病原」那类假阳性）。

### 4.3 GENEX：证明「系统性」并提供流程，但方向相反

它解决的是「微生物 reads 被真核参考抢走」，我们的是「病毒 reads 被宿主参考抢走」——机理同、对象不同。
可复用的部分：60 bp 伪读长、Bowtie2 保留多映射（最多 10）以提高灵敏度、按基因组出 BED。
**不要**照搬它的结论阈值（「最高 70%」是 PhyloNorway 的植物基因组；枸杞未必如此）。

---

## 5. 方案：三档，独立验收，互不阻塞

```
P0 诊断（必须）            P1 contig 级证据（推荐）        P2 掩码宿主库（可选臂）
宿主基因组 EVE 扫描   ──>  候选 contig 与 EVE 位点相交 ──>  masked FASTA → 重建四索引
位点 BED + 掩码占比        新增"内源性"证据字段            3–5 样本 A/B 对比
（不改任何生产库）          （不改去宿主）                  （通过才推广）
```

### P0 · 诊断：先量出来「宿主基因组里有多少病毒样序列」

**输入**：`01_PublicData/host_reference/genome/all.genome.uniq.fasta`（或分物种的 3 个 FASTA）
**工具**：归档的 `archive/20261011_endogenous_virus_pipeline/eve_screen.py`

```bash
cd archive/20261011_endogenous_virus_pipeline

# 首扫建议 --fast（约快 6 倍，位点数约 -54%），拿到量级再决定要不要全量
python eve_screen.py -g ningxia.genome.fasta -n ningxia \
  --outdir 06_EVE --stage 1,2,3 --fast -t 60
python eve_screen.py -g heiguo.genome.fasta  -n heiguo  --outdir 06_EVE --fast -t 60
python eve_screen.py -g zhonghua.genome.fasta -n zhonghua --outdir 06_EVE --fast -t 60
python eve_screen.py --merge -o 06_EVE          # kingdom_summary.tsv / kingdom_loci_all.tsv.gz
```

**产物**：`06_EVE/01_Loci/<name>/<name>.loci.bed`、`06_EVE/02_Verdict/<name>/<name>.s2_verdict.tsv`、
`04_Summary/`、`kingdom_summary.tsv`。

**生成掩码 FASTA（本文件配套工具，只做坐标换算 + N 化，不依赖仓库其它模块）**：

```bash
python public_metadata_pipeline/mask_host_genome.py \
  -g ningxia.genome.fasta heiguo.genome.fasta zhonghua.genome.fasta \
  -v 06_EVE/kingdom_loci_all.tsv.gz \
  --pad 50 --merge-gap 50 \
  --out-dir 06_EVE/masked/
# -> 06_EVE/masked/<name>.masked.fasta + <name>.masked.report.tsv
```

**读法（这就是决策依据）**：

| `masked_frac`（合计） | 解释 | 决策 |
|---|---|---|
| `< 0.1%` 且无 Caulimoviridae/Geminiviridae | 枸杞基因组基本没有近期 EVE | **不掩码**，写进 Limitations 即可 |
| `0.1%–1%` | 有少量近期 EVE | 先做 P1（证据层）；P2 作为敏感性分析 |
| `> 1%`，或含大段 Caulimoviridae/Badnavirus | 明显 EPRV 负载 | 做 P2，并在 Methods 声明 |

> **口径提醒**：归档说明（`ARCHIVED_20261011.md`）写明这套 EVE 流程「方法不成熟」——判定阈值、
> 污染/退化区分、验证体系都不足。所以 P0 的数字只能用于**量级判断与掩码源生成**，
> 不能当作「枸杞有几个 EVE」的结论写进论文；P1 的证据层同理，只作辅助证据，单独不构成判定。

**同时做的三件排查**（都属于「参考库卫生」，与 EVE 并列）：

1. `--confidence` 实际值：从 `summary.tsv`/日志确认是 0.2 还是 0.4（稿件口径不一致，见 §2）。
2. 宿主库里有没有类病毒序列：用 PSTVd/TCDVd 参考 blastn 对 `all.genome.uniq.fasta` + 对少量样本的
   Kraken2 `classified` 池，若有高一致性命中 → 立即处理（这比 EVE 更要紧，因为占 16.4%）。
3. 病毒参考库污染（§4.2）：对 `final.cluster.ref.fasta` 做宿主污染扫描。

### P1 · contig 级「内源性」证据层（推荐先做，收益最高、风险为零）

不改去宿主、不换库，只给下游判定加一列证据：

```
候选 contig (final.fasta / rescue candidates)
        │  minimap2 -ax asm5  →  宿主基因组
        ▼
宿主侧翼 / 覆盖比例  ──┐
                       ├──> EVE 位点相交 (BED) ──> eve_locus_frac / host_flank / evidence_tier
EVE BED（P0 产物） ────┘
```

- 直接对应 `final_judgement_table.tsv` 里的「宿主污染 40」「新种候选 24/21」两类判定；
- 复用归档的 `eve_distinguish/s2c_refgenome_scan.py`、`eve_final_score.py` 的 tier 逻辑
  （`A_genomic_flank > B_locus_architecture > C_decay_only > D_none`）；
- 论文收益：把「内源 vs 外源」从 k-mer/同源层提升到**基因组侧证据**，这是 GENEX 类审稿意见的标准回应；
- 风险：零（只增证据，不动数据）。

### P2 · 掩码宿主库（只有 P0 显示高掩码率时才做）

```bash
# 1)+2) 一条命令完成「掩码 + 重建四索引」，与现库并行（新目录 hostdb_evemask，不覆盖旧库）
python public_metadata_pipeline/build_host_pipeline.py \
  --species "Lycium barbarum" --taxid 112863 --stage eve-mask,hostdb \
  --genome-fasta 01_PublicData/host_reference/genome/all.genome.uniq.fasta \
  --eve-verdict 06_EVE/kingdom_loci_all.tsv.gz --eve-pad 50 --eve-merge-gap 50 \
  --hostdb-tools kraken2,bowtie2,hisat2,minimap2 -t 60 \
  --seq-type dna-short,rna-short \
  --k2-libs archaea,bacteria,plasmid,fungi,protozoa,UniVec,UniVec_Core \
  --work-dir 01_PublicData/host_reference
# 产物: <work>/genome/all.genome.masked.fasta{,.report.tsv,.meta.json}
#       <work>/hostdb_evemask/{kraken2,bowtie2,hisat2,minimap2} + db_provenance.json

# 等价的低层写法（不想用 stage 时）：
#   python public_metadata_pipeline/mask_host_genome.py -g all.genome.uniq.fasta \
#     -v 06_EVE/kingdom_loci_all.tsv.gz --pad 50 --merge-gap 50 -o all.genome.masked.fasta
#   python public_metadata_pipeline/build_hostbase.py -i all.genome.masked.fasta \
#     -o hostdb_evemask_v1 --tool kraken2,bowtie2,hisat2,minimap2 -t 60 --taxid 112863 ...
```

**注意**：Kraken2 库与比对库必须**同源同批**重建，否则两层掩码不一致，Stage 2/4 行为会互相矛盾。

**A/B 验证（只跑 3–5 个已知阳性样本，成本可控）**：

| 指标 | 期望 | 通过阈值 |
|---|---|---|
| EVE **以外**区域的宿主去除率 | 与对照一致 | 差值 < 0.2%（证明没误伤正常去宿主） |
| EVE 区 reads 保留量 | 显著上升 | ≥ 对照的 2×，或「原被删 reads 中 ≥50% 被保留」 |
| 组装 contig 数 / N50 | 不劣化 | contig 数变化 < 10%，N50 不降 |
| 新增 contig 的病毒样比例 | 新增的主要是真病毒 | 人工看 top 20 新增 contig 的分类 |
| 15 个高置信病毒检出 | 不丢不增 | 100% 保持，且不新增假阳 |

**只有上表全过**才把掩码库用于新样本；**同一批样本必须用同一个库**（否则库间批次效应会污染多样本比较）。

---

## 6. 接线方式（已实现：`build_host_pipeline.py` 的 `eve-mask` 阶段，默认关闭）

`public_metadata_pipeline/build_host_pipeline.py` 的三阶段：

```
genome-down  ->  [eve-mask]  ->  hostdb
```

`eve-mask` 卡在中间：输入是 genome-down 产出的 `all.genome.uniq.fasta`（或 `--genome-fasta` 指定的 FASTA），输出 `all.genome.masked.fasta`（+ `.report.tsv` + `.meta.json`），hostdb 阶段再消费掩码后的 FASTA。**不启用时行为与改造前完全一致**（同样的 FASTA、同样的 `<work>/hostdb`、同样的索引布局）。

### 6.1 开关：给任意一个掩码源即启用

不显式打开就整段跳过——`--stage all` 也会把它从计划里剔除（不会留下「skipped」checkpoint 噪声）。

| 参数 | 含义 |
|---|---|
| `--eve-mask` | 显式开关（配合 `--eve-bed` / `--eve-verdict` 用）|
| `--eve-bed FILE...` | 直接给 BED（0 基半开）|
| `--eve-verdict FILE...` | 给 EVE 判定表（`<name>.s2_verdict.tsv` 或 `kingdom_loci_all.tsv.gz`），按 `--eve-verdict-values` 过滤，默认只掩 `viral_supported` |
| `--eve-scan` | 现场跑归档管线的 `eve_screen.py`（stage 1,2,3 `--fast` + `--merge`），产物落 `<work>/eve/`，再拿其结果当掩码源 |
| `--masked-genome FILE` | 用别人已掩码好的 FASTA：只登记 + 写 meta，不重算 |

其余参数：`--eve-pad`（默认 0）、`--eve-merge-gap`（默认 0）、`--eve-mask-mode {N,lower}`（默认 `N`）、`--eve-masked-fasta`（自定义输出路径）、`--masked-hostdb-dir`（自定义掩码库目录，默认 `<work>/hostdb_evemask`）、`--eve-scan-script`（默认指向归档脚本）。

### 6.2 三条命令

（1）已有 P0 判定表 → 建掩码库（推荐；已有基因组 FASTA，跳过下载）：

```bash
python public_metadata_pipeline/build_host_pipeline.py \
  --species "Lycium barbarum" --taxid 112863 --stage eve-mask,hostdb \
  --genome-fasta 01_PublicData/host_reference/genome/all.genome.uniq.fasta \
  --eve-verdict 06_EVE/kingdom_loci_all.tsv.gz --eve-pad 50 --eve-merge-gap 50 \
  --hostdb-tools kraken2,bowtie2,hisat2,minimap2 --threads 60 \
  --work-dir 01_PublicData/host_reference
```

（2）一步到位：让 stage 自己跑 EVE 扫描再建库：

```bash
python public_metadata_pipeline/build_host_pipeline.py \
  --species "Lycium barbarum" --taxid 112863 --stage eve-mask,hostdb --eve-scan \
  --genome-fasta 01_PublicData/host_reference/genome/all.genome.uniq.fasta \
  --eve-pad 50 --threads 60 --work-dir 01_PublicData/host_reference
```

（需要服务器上 `archive/20261011_endogenous_virus_pipeline/eve_screen.py` 与病毒库路径可用；脚本缺失直接报错，不静默跳过。）

（3）只建库、不掩码（等价于现状）：

```bash
python public_metadata_pipeline/build_host_pipeline.py \
  --species "Lycium barbarum" --taxid 112863 --stage hostdb \
  --genome-fasta 01_PublicData/host_reference/genome/all.genome.uniq.fasta \
  --threads 60 --work-dir 01_PublicData/host_reference
```

### 6.3 产物与自证

`<work>/genome/all.genome.masked.fasta` 旁边两个旁车文件：

- `.report.tsv`：逐 contig `genome  contig  length  n_loci  masked_bp  masked_frac`
- `.meta.json`（schema `mmpv-eve-mask/1`）：源/掩码 FASTA 的 md5、掩码源路径、`pad/merge_gap/mode`、`totals`、`advice`（掩码率判读建议）

`<work>/hostdb_evemask/db_provenance.json`（schema `mmpv-hostdb-provenance/1`）：species/taxid、构建所用基因组（路径 + md5）、`masked` 布尔、内嵌完整的 `eve_mask_meta`、工具与 seq-types、k2 libs、线程数。**这就是「这批样本用的是哪个库」的追溯凭证**；未掩码的库也会写（`masked: false`）。

**安全闸**：若判定表的 contig 名与 FASTA 头对不上（一条区间都没落到 contig 上），阶段会删掉刚生成的 masked FASTA 并失败（`no interval matched any contig`）——绝不留下「其实没掩码」的库冒充掩码库。

### 6.4 下游衔接

新库目录布局与 `build_hostbase.py` 完全一致（`kraken2/`、`bowtie2/host.*`、`hisat2/host.*`、`minimap2/host_<seqtype>.mmi`），所以 `data_preprocessing_pipeline` **不用改代码**，显式指过去即可：

```bash
python data_preprocessing_pipeline/data_preprocessing.py \
  --host_db 01_PublicData/host_reference/hostdb_evemask \
  --input <样本表> --output 02_Preprocessing --stage deplete \
  --tool bowtie2 --seq-type rna-short --threads 60
```

⚠ **坑**：`data_preprocessing.py` 的 `ensure_hostdb()` 只会自动找 `<host_reference>/hostdb`，**不会**发现 `hostdb_evemask/`。所以必须显式 `--host_db`，否则它会静默回落到旧库（或重新下载建库）。建议在新样本表里加一列「库版本」，并在 `summary.tsv` / 报告里带 `host_db_version`（当前管线不会自动写，需人工登记）。

`preprocess_unified.py` 的 `hostref` 阶段**不经过本编排器**（它直接调 `download_host_genome.py` + `build_hostbase.py`，路径写死 `hostdb/`），所以它永远不会自动用掩码库；要用就显式指索引：

```bash
python public_metadata_pipeline/preprocess_unified.py \
  --kraken2-db <work>/hostdb_evemask/kraken2 \
  --step2-index <work>/hostdb_evemask/bowtie2/host ...
```

---

## 7. 与稿件的关系

### 7.1 如果最终**不掩码**（当前最可能）

在 Limitations 增一句（现有 §Limitations 已声明 RNA-only，紧接其后）：

> **Host endogenous viral elements were not masked.** The Lycium reference genomes were used in full;
> reads mapping to endogenous viral elements (EVEs) are therefore classified as host and removed,
> which may reduce sensitivity for exogenous viruses closely related to integrated elements.
> A genome-wide EVE scan (`eve_screen`) and an optional masked-library arm were run as sensitivity
> analyses (<P0 结论一句话：X% of the genome masked, N loci, no Caulimoviridae/Badnavirus / or …>).

### 7.2 如果最终**掩码**（P2 通过）

Methods 的 host-depletion 段加一句：

> Endogenous viral elements (EVEs) annotated in the host genomes (n = N loci, X% of genome) were
> soft-masked (replaced by N) prior to host-database construction; all four indices (Kraken2, Bowtie2,
> HISAT2, Minimap2) were rebuilt from the masked FASTA so that EVE regions neither contribute k-mers nor
> alignment seeds. The full-length host database was retained for the main analysis, and the masked
> database was used as a sensitivity arm on the same samples.

两段都可以直接落到 `doc/METHODS_TEMPLATE.md` 的对应小节。

---

## 8. 待验伪清单（做之前先把这些钉死）

| # | 假设 | 验证方法 | 不成立的后果 |
|---|---|---|---|
| 1 | `RVDB_viroids.fasta` 含 DNA 病毒（非 viroid-only） | `seqkit stats -a` + blastn 探针 | 伪读长掩码为空操作 |
| 2 | 枸杞三个基因组有大段近期 EVE | P0 `eve_screen --fast` + `masked_frac` | 掩码收益 ≈ 0 |
| 3 | 现库 `--confidence` 实际值 | 日志/`summary.tsv` | 稿件参数不一致；假阴入口大小估计错 |
| 4 | 宿主库不含类病毒序列 | PSTVd/TCDVd blastn 对宿主 FASTA | **16.4% 主力信号被合法删除** |
| 5 | 病毒参考库无宿主污染（方向 B） | 参考序列 → 宿主基因组比对，覆盖≥50%/一致性≥95% | 假阳性检出 |
| 6 | 掩码后 EVE 区 reads 真的被保留 | P2 的 A/B 表 | 改了库但没解决问题 |

---

## 9. 参考文献

1. Lu J, Salzberg SL. Removing contaminants from databases of draft genomes.
   *PLOS Comput Biol* 2018;14(6):e1006277. https://doi.org/10.1371/journal.pcbi.1006277
   —— EuPathDB 245 个真核病原体基因组的去污染流程；`Plasmodium yoelii` 三株 ~60% 污染、
   `Entamoeba histolytica` Rahman ~50%；20 例人角膜样本验证：原始库把 *Toxoplasma gondii* 判为
   几乎所有样本的主病原，清洗库才正确识别 *Acanthamoeba* / *Aspergillus* / *Candida*。
2. Oskolkov N, et al. Improving taxonomic inference in ancient environmental metagenomes by masking
   microbe-like regions in eukaryotic reference genomes. *GigaScience* 2025;14:giaf108.
   https://doi.org/10.1093/gigascience/giaf108 — GENEX：GTDB 微生物基因组切 60 bp 伪读长、
   Bowtie2（保留 ≤10 多映射）比对约 4,300 个真核参考基因组；真核基因组中微生物样序列最高达 70%
   （*Claytonia eschscholtzii*），非哺乳脊椎动物 <0.2%；北极/格陵兰古 DNA 中原本判为
   *Hippuris vulgaris* 的 reads 分别有 97% / 74% 落在其微生物样区域。
3. Johnson M. *mask-reference-genome*（`reference_masker.py`）.
   https://github.com/megjohnson1999/mask-reference-genome —— 病毒库切 60–80 bp 伪片段（overlap 40、
   median 80、variance 5、min-length 60、min-identity 90、max-indel 2、entropy 0.5）比对参考基因组并屏蔽。
