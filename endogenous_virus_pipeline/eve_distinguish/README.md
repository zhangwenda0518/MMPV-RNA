# eve_distinguish — DNA 病毒候选的 EVE 判别 (v6.3)

对 OneKP 新候选（`<DISCOVERY>/eve_screen/evescan_query.fasta`，1,027 条 = 843 远缘 + 184 新种）
做"结构退化 × 基因座架构"判别，与 `eve_screen.py` 的 EVE 库同源检索**相互独立、互为补充**。

> **版本标记约定**：模块版本只有这里和 `run_all.sh` 头两处权威声明（当前 **v6.3**）。各 stage 文件里写的 `[判定逻辑最后变更: 模块 vN]` 指该文件自己的逻辑最后改于模块哪个版本，不是并列的另一套版本号 —— s2b 停在 v5.1 是正常的，它此后未再改动判定逻辑；s1/s2c 的判定逻辑最后变更于 v6.3。

> **v6.2 → v6.3 改了四件事（其中①②会改产线数字，见下方"v6.2 → v6.3"一节）**：
> ① s1 修 `stop_profile` 的 3nt 系统性高估；② 面板升 v3（bare polyprotein 兜底 +
> `--supplement` 属级补充，Badnavirus 13→39 条）；③ `verify/rerun_diff.py` 派生表
> 重算差分工具；④ 新增 s2c 参考基因组通道（可选 `-G`，详见 `run_all.sh` 头部与
> `s2c_refgenome_scan.py`）。

> **v5.1 修了四个会改产线数字的 bug**：s2b 把无寄主命中的候选并成伪位点、
> s2 的组件并集按字符迭代、s3 的寄主否决分不清自身物种与跨物种、s2b 输出漏表头。
> 详见"v4 → v5.1 修的问题"一节，`EVE_STRONG_provirus` 40→19、`MOVE_EVE` 228→68 等数字
> 都是这些 bug 的直接影响。**引用过 40/228/55/22 这组数要先换成 19/68/48/8；但 v6.1
> 又改了五处（见"v5.1 → v6.1 修的问题"），现行数字是 `REMOVE` 2 / `MOVE_EVE` 174 /
> `KEEP_virus` 2 / `REVIEW` 849，见"交接"一节——两轮数字不要混用。**

**v5 的定位变化：本模块不再是 `virome_discovery_pipeline` 的一个 stage，而是跑完发现
管道之后单独调用的后运行脚本。** v3 之前它挂在 `virome_pipeline.py --stage eve` 上，
输入输出全靠主管道用环境变量注入，忘开就整类假阳性静默漏掉；v5 起它搬进
`endogenous_virus_pipeline/`，自带 CLI 和完整参数校验，`virome_pipeline.py` 里的 eve
阶段接入已整体移除。判别逻辑（verdict→action 映射）也随之从主管道搬回本模块的
`s4_filter.py`——以前改一头忘另一头，现在两处在一起。

v4 要解决的是一个具体问题：**怎么把新 DNA 病毒和 EVE 分开**。两个方法各自的边界：

| 方法 | 有效 | 无效 / 陷阱 |
|---|---|---|
| 方案一 框内提前终止 | 古老整合（92 个终止的 6.8kb contig 是铁证） | 近期整合（还没积累终止）；被组装 indel 污染 |
| 方案二 基因功能完整性 | 判定"是不是整合的前病毒" | 不能单独用：RNA-seq 里真 DNA 病毒本来就拼不全，"不全"≠EVE |

v4 的两条核心修正都针对上表的右列：方案一改用**同批组装自身的错误率做基线**；
方案二把方向**倒置**过来——不再问"完不完整"，只问"单一位点上是否 MP+CP+AP+RT+RH
齐全且密码子退化"。

## 运行

```bash
# 最常用: 指向一次发现管道的输出目录, 候选/证据表按目录约定自动定位
#   候选  <DISCOVERY>/eve_screen/evescan_query.fasta  (备选 10_Reports/eve_candidates.fasta)
#   证据  <DISCOVERY>/10_Reports/rescue_evidence_scored.tsv
#         (只跑到 analysis_verify 的树里退回 09b_Analysis_Verify/virus_validation/)
bash run_all.sh -D /path/to/onekp-virus -A /path/to/1kp/assemblies -t 32

# 显式给候选 fasta 与证据表, 产物落自己指定的目录
bash run_all.sh -q candidates.fasta -e 10_Reports/rescue_evidence_scored.tsv \
                -A /path/to/1kp/assemblies -o out/ -t 32

# 只验证 S1/S2 通道 (不跑寄主比对): host_contamination_likely 检不出, 仅供调试
bash run_all.sh -q candidates.fasta -e evidence.tsv -o out/ --no-host

# 加开参考基因组通道 (v6.3, 可选): 有现成基因组的物种做整合信号判定,
# 清单里没有的物种自动跳过 (格式见 genome_manifest.template.tsv)
bash run_all.sh -q candidates.fasta -e evidence.tsv -A /asm -t 32 \
                -G genome_manifest.tsv
```

| 参数 | 默认 | 说明 |
|---|---|---|
| `-q/--query FILE` | 无 | 待判别候选 fasta（必填，除非用 `-D` 自动定位） |
| `-e/--evidence FILE` | 无 | `rescue_evidence_scored.tsv`，s3 靠它取 tax_family/CheckV（必填） |
| `-D/--discovery DIR` | 无 | 发现管道一次输出目录；自动找 query/evidence，OUT 默认 `<DISCOVERY>/08b_EVE_Distinguish` |
| `-o/--outdir DIR` | `-D` 时如上，否则 `./eve_distinguish_out` | 产物目录 |
| `-A/--asm-root DIR` | 无 | OneKP 转录组汇编根目录（映射表第三列相对它）；**寄主否决必需** |
| `-m/--host-map FILE` | 模块自带 `host_sample_map.tsv` | 三列 TSV：样本号 `<TAB>` 4位代码 `<TAB>` 物种目录名 |
| `--baits FILE` | 模块自带 `baits.fa` | Caulifinder banks（TRUE/FALSE 编码头） |
| `--locus-query FILE` | 同 `-q` | 只对子集做寄主 blastn（省时） |
| `-G/--genome-manifest FILE` | 无（通道关闭） | 参考基因组清单，三列 TSV：4位代码 `<TAB>` 注释 `<TAB>` fasta 路径（模板 `genome_manifest.template.tsv`）；给开才跑 s2c，清单外物种记 `no_ref_genome` 跳过；需要 `-m` 做样本归属 |
| `-t/--threads N` | 16 | DIAMOND 线程数 |
| `--subj-cap-gb N` | 40 | 解压后寄主组装超过这个体积就跳过，防 OOM |
| `--no-host` | 关 | 不跑寄主比对（S1/S2 仍可用；但 `host_contamination_likely` 检不出） |
| `-h/--help` | — | 打印完整用法 |

三条硬校验（缺一条直接退出 1，不静默降级）：

- 缺候选 fasta / 证据表 → 退出，提示用 `-q`/`-e` 或 `-D`；
- **缺 `-A` 且没给 `--no-host` → 退出**。寄主通道是唯一能识别"候选其实就是寄主基因"
  的手段，缺了 `host_contamination_likely` 一条都检不出来。第一版就因为 `ASM_ROOT`
  配错、23 个物种全记 `#NOGZ`，寄主污染只报出 16 条而实际 55 条；
- 找不到 `diamond` → 退出；`--no-host` 时不校验 `blastn`（根本不会调它）。

模块目录（`run_all.sh` 所在）只放**脚本与固定参考**：`panel.fasta`、`baits.fa`、
`host_sample_map.tsv`、`build_panel.py`；所有中间表与最终判别表全落在 `-o` 指定目录，
模块目录永远只读。`panel.fasta` 是随模块提交的固定产物，脚本直接复用（软链，Windows
上没有符号链接权限时自动退回复制），不重新联网抓取，保证结果可复现；只有模块里没有它
时才调 `build_panel.py`。

`host_sample_map.tsv` 三列 `样本号 <TAB> 4位代码 <TAB> ASM_ROOT 下的物种目录名`，
ASM_ROOT 即 `-A` 给的值（OneKP `CODE-SOAPdenovo-Trans-assembly.fa.gz`，注意
`denovo` 有两种拼写，脚本两种都试）。**这一份是可选的**：有它时自身寄主物种直接给定；
没有时 s2b 退化为按"最佳单条命中的 bitscore 所在物种"推断自身寄主，实测 1,027 条里
只有 1 条结论会变（见下）。

## 阶段

0. 参考面板 — 复用模块自带 `panel.fasta`（无则 `build_panel.py` 联网构建）
1. DIAMOND blastx（`--sensitive -e 1e-5`）vs panel + baits.fa
   （baits 来自 Caulifinder banks：`FALSE_*` = LTR 逆转录转座子域，命名条目 = Caulimoviridae）
2. `s1_decay_scan.py` 方案一：结构通道 **v4 — 同批组装错误率基线**
3. `s2_domain_scan.py` 方案二：功能通道 **v4 — 单基因座（不是 contig）架构**
4. `s2b_locus_scan.py`：候选 vs 自身寄主转录组 blastn → 位点归并；**无寄主命中的候选记
   `no_host_locus`**，不再参与"位点齐全"判定（`--no-host` 时整表只有表头，含义相同）
5. `s3_verdict.py` 合成判别 **v4 — 二维表 + 寄主同源否决** → `eve_distinguish_verdict.tsv`
6. `s4_filter.py` **v5 新增**：verdict→action 映射 → `dna_vs_eve_filter.tsv`（下游按列过滤）
7. `s2c_refgenome_scan.py` 参考基因组通道 **v6.3 新增**（可选 `-G`）：候选 vs 自身物种
   染色体级参考基因组 blastn → 整合信号判定 → `refgenome_evidence.tsv`；
   s3 以可选第 6 参单向并入（见下方"参考基因组通道"一节）

### 方案一 (s1_decay_scan.py v4)

v3 用绝对计数 `stops>=3` 判退化。问题是**一个 1nt indel 就足以让它之后的全部密码子错读**，
在远缘同源区（25-45% pid）上造出十几个"提前终止"，与真正的化石退化无法区分。

v4 的做法是拿 contig **自己**当基线：同源区同一段序列、同一条链、**另外两个读框**的
终止密度。同一批组装、同一读段、同一碱基组成，组装 indel 对三个读框的污染是同等的，
所以另外两框的终止密度就是这套组装的本地噪声底噪。

- `stop_enrichment = (病毒读框终止数+1) / (另外两框均数+1)`，**越大越脏**（= 病毒读框比
  另两框的噪声底噪脏多少倍）。<1 → `compositional_noise`
  （该读框并不比随机框更脏，终止来自组成/组装，无退化证据）
- ≥1 时再用**分布剖面**区分两种情形：
  - `max_stopfree_frac >= 0.25`：存在一段 ≥1/4 同源区长度的完整无终止区 → `assembly_breakpoint`
    （断裂点两侧各是一段完整可读区，最长无终止段必然占去大半）
  - `head_stops >= 2` 且无长完整段 → `distributed_decay`（化石签名）

  **不能用 "head==0 就推断裂"**：实测那条 92 终止的 6.8kb contig 前端恰好无终止
  （`head_stops=0`、`max_stopfree_frac=0.059`、`enrich=2.16`），显然不是单点断裂。
  **`head_stops` 要到 2 才算分布式退化**：单个头部终止同样可能只是那里恰好有个 indel，
  v4 的 `head_stops >= 1` 会把这类单 indel 化石误判成化石签名。

`compositional_noise` 里再分一档"编码其实完好"：`stop_enrichment <= CLEAN_ENRICH(0.5)`
且 `stops <= CLEAN_MAX_STOPS(5)` 时归 `coding_intact`——病毒读框比另两框干净一个数量级，
这不是"无退化证据"，而是"编码完好"，s3 据此在单一位点齐全时给 `virus_candidate`。
实测这一档从 `compositional_noise` 挪走了 10 条。

### 方案二 (s2_domain_scan.py + s2b_locus_scan.py)

v3 把"完整"当 EVE 证据（`loci>=2 或 span>=1500 → review_extend`），方向是反的：
RNA-seq 里真 DNA 病毒本来就拼不全（环状 7.5kb 基因组靠短读段只能拿到局部）。
v4 只问一件事——**单一位点上 MP+CP+AP+RT+RH 是否齐全**（`provirus_scale`：
单基因座 span≥1500nt 且 ≥4/5 组件），且必须与 S1 的密码子退化**联合**才有意义：

```
                 单一位点架构齐全        架构不全
密码子分布式退化   EVE_STRONG_provirus    EVE_suspect / ancient_EVE
                  (整合的前病毒)          (contig 级无法定案)
密码子完好        virus_candidate         virus_fragment_review
                  (感染中的真病毒, 要留)   (病毒碎片, 不判 EVE)
```

"位点"是**寄主基因组位点**，不是 contig——感染中的病毒会被打成若干个 contig 各带一两个
基因，而整合的前病毒在位点上是完整的、只是被组装切成碎片。所以 s2b 把候选 blastn 回
各物种 OneKP 转录组组装，按 `(样本, 物种, scaffold)` 聚合格座，输出 `locus_arch` /
`locus_ncomp` / `locus_hint` / `locus_key`。"s2 数组件时把融合蛋白（如 `AP+RT`）的两个
组件都算上，且**按整词合并**——同一位点内多个 HSP 的组件集取并集时按组件名并，
不是按字符并（否则 `MP` 会被拆成 `M`+`P`，第二个组件永远并不进来）。但**算进来不等于
有证据**：v6 起 `credit_components()` 只把"有独立 query 区间、且比对长度够放一个蛋白"的
组件计入 `best_locus_ncomp` / `provirus_scale`，每条 HSP 还按 `aln // COMP_AA_MIN` 封顶
——条目名里的四个组件是**参考序列**的事实，不是这段比对的事实。真正贡献了证据的 HSP
条数记在 `credit_hsps` 列。

**没有寄主命中的候选不归入任何位点**：旧版按 `(样本, host_scaffold)` 分组，而未命中
contig 的 `host_scaffold` 一律是 `-`，于是同一份样本里所有无命中 contig 被并成一个
巨型伪位点，各自的组件集取并集——实测 787 条无命中候选里有 **70 条因此被判成
`locus_full`**（例如 ERR2040363 样本 7 条互不相干的 contig 并出 `MP+CP+AP+RT+RH` 五件
套）。现在这些候选一律记 `no_host_locus`，s3 对它们最高只能给 `EVE_suspect`。

> OneKP 转录组组装本身是碎的，一个寄主基因常被打断在多个 scaffold 上（实测 HSP70 候选
> 最高单条比对只有 750nt/94.5%）。所以 s2b 的同源强度用**物种级 query 区间并集覆盖**，
> 不是单条比对覆盖。但"自身物种"与"跨物种旁证"要分开看：样本在 `host_sample_map.tsv`
> 里有条目才可能是 `host_scope=own_species`（实测只有 73/1027 条），其余 167 条命中来自
> 跨物种最佳比对，只是旁证（`host_scope=cross_species_guess`）。口径不同，采信强度不同。

#### 面板覆盖的现实（读数前先知道）

**v6.3 起面板为 v3（291 条）**：v2 的 222 条里 Badnavirus（科内最大属、已报道植物 EVE
的主力类群，Banana streak 系列即属此）只有 1 条 VAP——它不在 s2 的 CANON 五组件里，
等于该属的结构基因通道整属失明。根因是 Badnavirus 的 GenBank 蛋白大量只标
"polyprotein"，不含任何组件关键词，被 v2 的"无关键词即丢弃"规则整条扔掉。v3 修法：

- `build_panel.py` 标注规则加兜底：干净标题的 bare "polyprotein" 按科的 pol 区规范标
  **AP+RT+RH**（MP/CP 有独立 ORF/条目，不会因此虚标）；
- `build_panel.py --supplement` 属级补充模式：不重建面板（NCBI 内容漂移会让全部下游
  数字搬家），只为面板没有的物种按属级配额补条目（Badnavirus 25 种、其他属 8 种；
  每物种每组件组至多 1 条、取最长），追加写入，既有 222 条一条不动；
- 补充实测（2026-09-27）：222 → 291 条，Badnavirus 13 → 39（新增 25 个物种的
  MP/CP/AP+RT+RH），Ruflodivirus / Petuvirus / Vaccinivirus 三个零覆盖属补齐；
  Banana streak 系列 8 → 17 条。属级归属见 `panel_genus_coverage.tsv`。

v2 时代的口径记录（保留作历史对照）：222 条按组件展开 MP 62 / CP 62 / RT 60 / VAP 44 /
AP 15 / RH 6，其中 AP 和 RH 绝大多数只以融合体入池（`AP+RT` 13 条、`RT+RH` 3 条、
`MP+CP+RT+RH` 2 条；纯 AP 只有 2 条、纯 RH 只有 1 条）。后果：`locus_ncomp` 顶不满
4/5 可能只是**面板覆盖问题**而不是生物学事实，别把"组件少"直接读成"没整合完整"。

**注意既有 OUT 目录的复用语义**：`run_all.sh` 对已存在的 `panel.fasta` 会原样沿用
（同目录可复现），检测到模块面板更新时只打警告不自动换——要用 v3 面板请删掉 OUT 里的
`panel.fasta`/`panel.dmnd` 重跑或换新 `-o` 目录；面板换代后组件命中会变，新旧数字不可
直接对比。

### 寄主同源否决 (s3)

HSP70 假阳性的来源：`cd24029` = ASKHA_NBD_HSP70_DnaK_HscA_HscC（伴侣蛋白域，不是病毒域），
`pfam00012`/`cd10233` 同样。管线的 CDD 报告自己带 `hsp70cap=1` 标记，但两条 27-29% 流行率
的 DNA 候选仍凭这个域 PASS_VIRAL 混进新种候选——其一与**自身寄主**的转录组 scaffold
99.8%/2476nt 相同。污染源是 Phycodnaviridae（*Micromonas* 病毒，~197kb NCLDV 基因组）。

否定条件：**寄主聚合 ≥90% 一致且区间覆盖 ≥80%，且无 MP/CP/AP** → `host_contamination_likely`。
带 MP/CP/AP 的走 EVE 通道（那正是 EVE 问题，不是寄主污染问题）。
`host_xeno` 列另记最强的跨物种命中，用来发现样本串号；`host_scope` 列记这次命中的口径
（`own_species` / `cross_species_guess` / `none`），verdict_reason 里也会写明比较对象是谁
——旧版一律写成"自身寄主物种"，跨物种旁证被当成了自身物种证据。

## v3 → v4 踩坑记录

| # | v3 行为 | v4 修正 |
|---|---|---|
| 1 | 面板"基因组多聚蛋白"条目带全部组件标签，单条 HSP 即记 MP+CP+AP+RT 全分 → v2 已改为按 query 坐标合并基因座 | v4 进一步：只认**单一基因座**的 4/5 组件，跨基因座并集不算 |
| 2 | DIAMOND 负链 HSP `qstart>qend` 未归一化导致 locus 拆分虚高、终止计数为 0（v1 把 92 终止的 6.8kb contig 评成"五件套完整"） | 统一 (min,max) 归一化；v4 该条正确落到 `EVE_STRONG_provirus` |
| 3 | Caulimoviridae RT 与转座子 RT 同源，TE bitscore 对比误杀真病毒 | TE 否决仅适用于**无 MP/CP/AP**的 contig |
| 4 | `stops>=3` 绝对计数，组装 indel 污染 | 同链另两框基线 + 停止分布剖面 |
| 5 | "完整"当 EVE 证据 | 方向倒置：齐全+退化=前病毒；齐全+完好=真病毒 |
| 6 | 不比对寄主，无法识别寄主基因假阳性 | s2b 寄主 locus blastn + 物种级区间并集 + ≥90%/≥80% 否决 |
| 7 | `blastn -subject <(gunzip -c asm.gz)` 是**管道**，blastn 无法 seek，会把整个转录组汇编缓存进内存 —— dmesg 里被 OOM killer 干掉十几次（最大一次 rss 266 GB），寄主比对循环静默中断，`locus_blastn.tsv` 只有几十 KB | 先解压到临时文件再 `-subject <文件>`（blastn 流式读），并加 `SUBJ_CAP_GB`（默认 40 GB）体积闸门；解压失败/超限记进 `locus_errors.log`，跑完汇总"命中行数 / 成功 scaffold 数 / 跳过行数" |

## v4 → v5 踩坑记录

| # | v4 行为 | v5 修正 |
|---|---|---|
| 1 | 挂在 `virome_pipeline.py --stage eve` 上，`QUERY`/`EVIDENCE`/`OUT_DIR`/`ASM_ROOT` 全由主管道注入；忘开这一阶段就静默漏掉整类假阳性 | 独立后运行脚本 + 自有 CLI；缺 `-A` 直接报错，`--no-host` 才允许不带 |
| 2 | verdict→action 映射写在 `virome_pipeline.py` 的 `run_eve_distinguish` 里，与 s3 的 verdict 定义分居两处 | 搬到 `s4_filter.py`，与 verdict 生产者在同一模块；表外新 verdict 返回非零（已按 REVIEW 兜底写全），调用方不会不知道 |
| 3 | `OUT_DIR` 不给就落调用方 cwd，产物散落 | `-D` 时默认 `<DISCOVERY>/08b_EVE_Distinguish`，`-o` 显式覆盖 |
| 4 | Windows/Git Bash 上 `ln -sf` 失败（无符号链接权限），`set -e` 让整个运行当场终止 | `link_or_copy()`：软链不成就复制（面板/dmnd 都是只读复用，复制同样可复现） |
| 5 | 寄主 blastn 跑完不知道有没有真的比上（全 `#NOGZ` 也看不出来） | 循环后打印命中行数 / 成功 scaffold 数 / 跳过行数，零命中时明确警告"寄主同源否决整条失效" |
| 6 | `$DIAMOND blastx ... 2>/dev/null` 且不接 rc：diamond 被 OOM killer 干掉或段错误时 `set -e` 只让它静默退出，上游看不到任何原因 | `run_or_die`：失败即打印 `FAILED(rc)` + 完整命令行并中断，不留半成品表 |
| 7 | 0 命中时 `panel_hits.tsv` 缺失直接让下一步 `cat` 报看不懂的 "No such file or directory" | 比对后逐个补空表并打印"未生成, 按空表继续（疑似磁盘满/权限问题）" |
| 8 | Windows 上只装了 `python`（没有 `python3` 入口），`PYTHON=python3` 要等跑到第三个阶段才报 "command not found" | 启动时探测：无 `python3` 就降级 `python` 并打印一行；两个都没有则立刻报错 |
| 9 | 无 `baits.dmnd` 时 `cat panel_hits.tsv baits_hits.tsv` 直接中断 | 空表分支：缺库只让 TE 否决不可用，主流程跑完 |
| 10 | `--no-host` / 缺 `HOST_MAP` 时 `locus_architecture.tsv` 只有 3 列（重定向只写在三行 `printf` 的最后一行，前几列漏进 stdout） | `{ ...; } > file` 包住整段表头，两种分支都是完整 16 列；v5.1 起改成 `--emit-header` 由模块自己打表头，列数再变也不会和脚本脱节 |

## 分布对比（1,027 条实测，宿主 23 物种全量比对）

| verdict | v3 | v4 | v5 修 bug 后 | v6.1 | 含义 |
|---|---|---|---|---|---|
| review | 640 | 530 | 530 | 530 | 证据不足（多为 RNA 病毒科，Cauli 面板不适用） |
| EVE_LTR_TE | 221 | 171 | 48 | 155 | 转座子科 / 纯 pol 区且 TE 领先 |
| virus_fragment_review | — | 69 | 127 | 95 | 架构不全 + 编码完好：**真 DNA 病毒本就拼不全，不判 EVE** |
| EVE_suspect | 3 | 71 | 162 | 113 | 架构不全 + 分布式退化，contig 级无法定案 |
| host_contamination_likely | — | 55 | 48 | **2** | **寄主序列，非病毒**；v6 起要求三条上游同源通道全无记录，其余 51 条转 `host_conflict_review` 33 / `host_homology_cross_species` 18 |
| host_conflict_review | — | — | 0 | 33 | 自身物种命中 + 上游有病毒证据：两通道矛盾，不仲裁、不删 |
| host_homology_cross_species | — | — | 0 | 18 | 跨物种命中，同源强度再高也只是旁证 |
| EVE_STRONG_provirus | — | 40 | 19 | 4 | 单一位点 4/5 组件齐全 + 密码子退化 = 整合前病毒；v6 起另要求 `locus_full` |
| assembly_breakpoint_review | — | 28 | 55 | 35 | 终止集中在单个断裂点，救长后复判 |
| virus_candidate | 23 | 22 | 8 | 2 | 单一位点齐全 + 密码子完好 = 感染中的真病毒（**要留的**） |
| compositional_noise_review | 20 | 20 | 25 | 20 | 病毒读框不比其他读框更脏，无退化证据 |
| ancient_EVE | 142 | 17 | 1 | 15 | 纯 pol 区分布式退化 |
| structure_intact_review | 7 | 4 | 4 | 5 | 非 Cauli 结构完整 |

对应的 action：`REVIEW` 849 / `MOVE_EVE` 174 / `REMOVE_host_contamination` 2 /
`KEEP_virus` 2（v6.1）。**`REMOVE` 从 53 条掉到 2 条是本轮最大的数字变化**，但根因不是
"寄主污染变少了"，而是两个口径都改严了：33 条是自身物种命中却同时带病毒侧信号（旧版
一刀切删，现在转 `host_conflict_review` 人工定案），18 条只是跨物种旁证、样本压根没有
自身物种映射（旧版据此删，现在改记 `host_homology_cross_species`，不删）。

结构退化类的重新归因（`assembly_breakpoint_review` 55→28、`EVE_STRONG_provirus` 0→40、
`EVE_suspect` 3→71）来自修正 #4/#6：v3 的绝对阈值把组装断裂和化石退化混在一档。

> **v5 修 bug 后那一列的 123 条 `EVE_LTR_TE` 离开是修正，不是退化**：这些 contig 的
> 融合面板条目（如 `AP+RT|sp|Q02964.1|POL_CAMVE`）本来就同时覆盖两个组件，旧版把组件并集
> 按字符迭代（`for c in h["comps"]`），`AP+RT` 并进已有的 `RT` 变成 `RT+A+P`，于是
> `best_locus_ncomp` 只剩 1、`components` 只剩 `RT`，s3 据此判定"纯 pol 区" → `EVE_LTR_TE`。
> 修好后它们拿到真实组件（`RT+AP` / `RT+CP` / `RT+CP+AP+RH`…），按退化与架构重新判：
> 50 条 `EVE_suspect`、37 条 `virus_fragment_review`、27 条 `assembly_breakpoint_review`、
> 8 条 `compositional_noise_review`、1 条升级为 `EVE_STRONG_provirus`。同一个 bug 也是
> `best_locus_ncomp` 从"只有 1/2/4"变成 1–5 全覆盖的原因（274→31 条卡在 ncomp=1）。

> **宿主比对必须跑完**：`host_contamination_likely` 只在**自身物种**出现在 `HOST_MAP`
> 且该物种的 blastn 已完成时才可能触发。第一版判别只跑完一部分物种（blastn 被 OOM killer
> 打断，见坑 #7），只报出 16 条；23 个物种全量跑完后是 **55 条**（`review`→host 32、
> `structure_intact_review`→3、`virus_fragment_review`→2、`EVE_LTR_TE`→2，另有个别位点完备度
> 变化使 2 条 `EVE_STRONG_provirus` 降为 `EVE_suspect`）。新增 39 条的自身物种聚合是中位
> 99.4% / 0.98 覆盖，而它们的病毒证据全是远缘弱命中（aa 一致度中位 57%、最低 26%，多为
> KOG/pfam 通用域），正是寄主基因被误注释成病毒的样子。跑完 `locus_blastn.tsv` 后应确认
> `locus_errors.log` 里没有 `#NOGZ` / `#TOOBIG` / `#GUNZIP_FAIL`，并核对 v5 加的那行
> "命中行数 / 成功 scaffold 数"不为 0。

## v4 → v5.1 修的问题（服务器 1,027 条实测，前四个会改产线数字）

s1/s2/s2b/s3 各改了一两处，都不是加功能而是**修掉已经在悄悄改数字的错误**。前四条都在
`~/eve_verify/` 用真实输入重跑 s1→s4 复现过（不碰产线 `server-wip-20260905` 分支），
第五条是同一轮改出来的表头回归。

| # | 症状 | 根因 | 修后实测 |
|---|---|---|---|
| 1 | 787 条无寄主命中的候选中 70 条被判 `locus_full`，进而出 40 条 `EVE_STRONG_provirus` | s2b 按 `(样本, host_scaffold)` 分组，未命中者的 `host_scaffold` 全是 `-`，同样本未命中 contig 被并成一个伪位点，组件集跨 contig 取并集 | 未命中候选一律 `no_host_locus`；`EVE_STRONG_provirus` 40→19，多出的转到 `EVE_suspect`（+91）。不变量：`locus_full` 24 条全部挂在真实位点上、组件数均 ≥4 |
| 2 | `best_locus_comps` 出现 `RT+A+P`、`MP+++C`、`RT+H+A+P` 这种垃圾值（334 条） | `merge_loci` 里 `for c in h["comps"]` 迭代的是**字符串**，`AP+RT` 被拆成字符；只有第一个 HSP 的首个组件能完整留下 | 按整词并集：`RT+A+P`→`RT+AP`。非 CANON 片段 334→0；`best_locus_ncomp` 从只有 {1,2,4} 变成 1–5 全覆盖 |
| 3 | 123 条带真实病毒组件（RT+CP/AP/RH）的候被判 `EVE_LTR_TE` 丢掉 | 同 #2：`components` 少算组件 → s3 的 `pol_only` 误真 → "纯 pol 区 TE" | 这 123 条改按退化/架构判：50 `EVE_suspect`、37 `virus_fragment_review`、27 `assembly_breakpoint_review`、8 `compositional_noise_review`、1 升级 |
| 4 | 只跑一部分寄主物种时 `host_contamination_likely` 的 reason 一律写"自身寄主物种" | `host_scope` 不存在，跨物种旁证与自身物种命中混在一句话里，下游无法按口径过滤 | 新增 `host_scope` 列（实测 `own_species` 73 / `cross_species_guess` 167 / `none` 787），reason 分开写；s3/s4 表头随之多一列 |
| 5 | `s2b` 写出的 `locus_architecture.tsv` 没有表头行；`--no-host` 分支的表头只进了 stdout | 上一版把 `out.write(HDR)` 放进了 `--emit-header` 分支，正常写文件时漏掉 | s3 是按首行当列名读这份表的，没有表头就会把第一条候选当表头吃掉、后面整行列错位；修后两种分支共用同一个 `HDR` 常量，写文件必带表头 |

`locus_architecture.tsv` 从 16 列变 20 列（多出 `sample_flag` / `host_scope` /
`n_host_scaffolds` / `locus_key`），`eve_distinguish_verdict.tsv` 19→20 列、
`dna_vs_eve_filter.tsv` 12→13 列（都只多 `host_scope`）。四个阶段的表头都抽成了模块级
常量（`S1_HDR` / `S2_HDR` / `VERDICT_HDR` / `s2b.HDR`）并由测试直接 import，加列时不会
只有一边改。

## v5.1 → v6.1 修的问题（同一份 1,027 条输入逐条对比）

上一节那些改完之后，还能用同一套输入再改出四个**已经在悄悄改数字**的问题。这一轮的共同
前提：组件标签来自面板条目名，而**条目名是参考序列的事实，不是这段比对的事实**——一条
104 aa 的比对打中写着 `MP+CP+AP+RT` 的融合条目，旧版就记 4 个组件。

| # | 症状 | 根因 | 修后实测 |
|---|---|---|---|
| 1 | 一条 104 aa 的比对记 4 个组件，`best_locus_ncomp=4` 直接送进 `EVE_STRONG_provirus` | `merge_loci` 的组件并集按条目名取，不问这段比对撑不撑得起 | 新增 `credit_components()`：每条组件要有**独立 query 区间** + 至少 `COMP_AA_MIN=100` aa（=300 nt）的长度支撑，并按 `aln // comp_aa_min` 封顶。`best_locus_ncomp>=4` 48→18，`provirus_scale` 12→11，新增 `credit_hsps` 列 |
| 2 | `locus_full` 与 `provirus_scale` 取或 | s2b 是"这个基因座在寄主基因组里成立吗"的唯一判据；s2 的 contig 级规模只是旁证 | 纵向集成只认 `locus_full`；`locus_ncomp` 改取 s2b 值。6 条 `no_host_locus` 的候选退出 `EVE_STRONG_provirus` |
| 3 | 53 条 `host_contamination_likely`（会删数据）里逐条查上游**全都带病毒科注释** | 上游给了病毒科就是病毒侧的独立证据，与"这就是寄主基因"直接矛盾，旧版不理会照删 | 寄主否决拆三分支（见下），只有 `host_contamination_likely` 允许 REMOVE：53→2 |
| 4 | 有 2 条候选只带了 blastn 命中（71.1%/72.1%）就被判"无病毒信号" | "上游病毒信号"只数蛋白通道 + 科注释，漏了核苷酸通道 | 补上 `nt_pident>0 且 nt_species 非空`（v6.1）；`host_contamination_likely` 4→2 |
| 5 | `loci_detail` 印成 `L1:1677-4094({'MP', 'CP'},bit560)` | 直接把 `credit_components()` 返回的 set 抖进 f-string，既带 repr 又依赖 set 迭代顺序，**同一份输入两次跑出的这一列可能不同** | 新增 `fmt_comps()`，按 CANON 序输出 `MP+CP+RT+RH`。不改任何数字，但 s2b 按这列分组，列不可复现就没了（实测两次运行 `s2_domains.tsv` md5 相同） |

`credit_hsps` 是 `s2_domains.tsv` 的第 18 列（表头 `S2_HDR` 从模块 import，测试直接引用）。

### `multi_locus_arch` 现状：这一列现在是死的（20→0）

必须如实说：v6.1 在 1,027 条上 `multi_locus_arch` **全为 FALSE**（旧版 20 条 TRUE）。
这 20→0 是**两个改动叠出来的**，不是一次跌到底：

1. v5 把定义从永真式（"组件≥1 的基因座 ≥2"——每条基因座必有组件，恒真）改成
   "基因座 ≥2 **且**合计 ≥3 个不同组件"。这一步就把 20 条砍到 6 条。
2. v6 的组件计数修复把剩下 6 条的合计组件数全部打到 3 以下 → 0 条。

逐条查过第二步动的那 6 条，旧的 TRUE 全部来自超发计数（`credit_hsps` 是新列，记录真正
贡献了证据的 HSP 条数）：

| contig | 旧版 components / best_ncomp | 证据（新表 `evidence_detail`） | v6.1 components / best_ncomp / credit_hsps |
|---|---|---|---|
| `ERR2040603_..._NODE_3783` | `MP+CP+RT+RH` / 4 | 一条 104 aa / 29% 的一致度打中四组分融合条目，另三条 88 aa / 38% 的 RT | `MP` / 1 / 1 |
| `ERR2040665_..._NODE_15791` | `RT+MP+CP+RH` / 4 | 三条 52 aa / 46%，一条 90 aa / 33% 打中 `MP+CP+RT+RH` | `-` / 0 / 0 |
| `ERR2040332_..._NODE_791` | `RT+RH+CP` / 2 | 168 aa / 36% 打中 `AP+RT`，另两段 ~498 nt 与它几乎完全重叠 | `CP+AP` / 1 / 2 |
| `ERR2040573_..._NODE_8849` | `MP+AP+RT` / 2 | 87/71/77/78 aa、42–54% 一致度，全部够不到 100 aa 闸门 | `-` / 0 / 0 |
| `ERR2040628_..._NODE_162` | `AP+RT+CP` / 2 | 最长三段 MP/MP+CP 分别 321/339/342 aa，但彼此重叠在同一段 query 上 | `MP+AP` / 2 / 2 |
| `ERR2040363_..._NODE_6522` | `AP+RT+RH` / 2 | 318 nt / 116 aa 打中 `AP`；另两段 150 nt/80 aa、228 nt/80 aa 在闸门下 | `AP` / 1 / 1 |

所以这不是阈值要调，而是**先前的 TRUE 本来就不成立**。这一列保留（定义是对的，测试
锁着"非永真式"），但在当前数据集上没有判别力——下游不要按它筛东西，`best_locus_ncomp` /
`credit_hsps` 才有信息量。若将来要救活它，方向是把 `MULTI_MIN_COMPS=3` 与
`COMP_AA_MIN=100` 联动下调并重做敏感性分析，而不是把闸门撤回旧值。

### 寄主否决三分支（v6 起）

| verdict | 条件 | action |
|---|---|---|
| `host_homology_cross_species` | 跨物种命中（样本无自身物种映射） | REVIEW |
| `host_conflict_review` | 自身物种命中 **且** 上游三条同源通道任一有记录 | REVIEW |
| `host_contamination_likely` | 自身物种命中 **且** 三条通道全无记录 | `REMOVE_host_contamination` |

"上游病毒信号"只认三条通道：`tax_family` 非空、`aa_pident>=95%` 且 `aa_species` 非空、
`nt_pident>0` 且 `nt_species` 非空。**故意不认 `cdd_top` 与低一致性蛋白命中**：这类命中
里最多的是 dUTPase / RdP / RNase H 这类通用结构域（`KOG0947` / `KOG1098` / `pfam07780` /
`COG0293`），寄主基因照样带，把它们算成"病毒侧证据"会让矛盾分支无限膨胀。

## 关键个例（v6.1 实测；数字与上一版不同，见下）

| contig | v4 实测 | v6.1 实测 | 依据 |
|---|---|---|---|
| `ERR2040732_..._NODE_11` 6.8kb | `EVE_STRONG_provirus` | `EVE_suspect` | 92 终止 / enrich 2.16 / msf 0.059 + MP+CP+RT+RH 单位点 2417nt 都在，但 `host_scope=none`、`locus_arch=no_host_locus`（0 个寄主 scaffold） |
| `ERR2040803_..._NODE_1874` 2.5kb | `host_contamination_likely` | `host_conflict_review` | 自身寄主 IAJW 聚合 91.2%/0.985，无病毒结构基因；但核苷酸命中 71.1% 一致到 `PP728250.1` |
| `ERR3487378_..._NODE_1805` 2.1kb | `host_contamination_likely` | `host_conflict_review` | 自身寄主 BLAJ 聚合 94.4%/0.915，无病毒结构基因；但核苷酸命中 72.1% 一致到 `gb|HQ633072.1|` |
| `ERR2040641_..._NODE_491` 2.8kb | — | `host_contamination_likely`（当前仅 2 条 REMOVE 之一） | 自身寄主 EITK 聚合 98.4%/0.99，三条上游同源通道全无记录 |

**更正一条此前的错误陈述**：v4/v5.1 记录的"`ERR2040732_..._NODE_11` 的寄主位点是真实
存在的，不是伪位点"是错的。它在 `locus_architecture.tsv` 里是 `no_host_locus`、
`host_scope=none`、`n_host_scaffolds=0`——**这个候选一条寄主命中都没有**。v6 之前
`locus_full or provirus_scale` 的取或逻辑让"contig 级单基因座齐全"冒充了"寄主基因座
齐全"，把它一路送进 `EVE_STRONG_provirus`；改掉之后它退回 `EVE_suspect`（退化与
架构证据都在，缺的是基因座级证据）。它的 `credit_hsps=1`：2417 nt 的单个比对
（L1:1677-4094, bit560）支撑 MP+CP+RT+RH 四个组件——单条长比对打中四组分融合条目时，
四个组件是这段比对撑得起来的，不属于超发计数。

## `--no-host` 实测（同一份 1,027 条输入）

不跑寄主 blastn 时 `locus_architecture.tsv` 只有表头（0 条数据），所有候选记
`no_host_locus`，`host_scope` 全为 `none`。后果：`REMOVE_host_contamination` 归零
（检不出寄主基因），`EVE_STRONG_provirus` 与 `virus_candidate` **同时归零**（两者都
要求 `locus_full`，没有寄主比对就无从判定），只剩 `EVE_suspect` /
`virus_fragment_review` 这类 contig 级结论。verdict / action 分布：

```
review 574 / EVE_LTR_TE 156 / EVE_suspect 117 / virus_fragment_review 100
assembly_breakpoint_review 35 / compositional_noise_review 20
ancient_EVE 17 / structure_intact_review 8        (REVIEW 854 / MOVE_EVE 173)
```

（`host_conflict_review` / `host_homology_cross_species` / `host_contamination_likely` /
`EVE_STRONG_provirus` / `virus_candidate` 在无寄主比对时全部归零——它们都要读 s2b 的
位点结论。）

自校验口径要改一句：v5.1 记录的"`EVE_LTR_TE` 只取决于 s1/s2，与跑不跑寄主完全一致"
现在**不再严格成立**——156（无寄主）对 155（有寄主），差 1 条。差的是
`ERR2040626_clean_NODE_11631`：它有 `te_flag=TRUE` + 纯 RT 区，本身是转座子，但同时
命中自身物种 DLJZ（100%/0.84）且上游归了 Caulimoviridae，寄主否决分支先一步接管，
判 `host_conflict_review`。这是有意的优先级（寄主矛盾证据强于 TE 判据），但说明
`EVE_LTR_TE` 不再是可用来交叉验证的"纯 s1/s2 列"；要校验表头-only 输入有没有带偏
s3，应改用 `assembly_breakpoint_review`（35 = 35，完全一致）。

## 判据被校准到什么程度（v6.2 新增，`verify/` 实测）

上面所有小节说的都是"代码按定义执行"。**判据本身站不站得住**是另一个问题，四轮实测
见 [`verify/README.md`](verify/README.md)。结论直接改变对外口径，摘要如下：

**1. 整条流水线对随机序列零假阳性（强证据）。** 把 1,027 条候选整条单核苷酸洗牌
（严格保碱基组成）后重跑全流程：diamond **0 命中**、s2 组件 **0 条**、s3 实质判定
**0 条**，全部 1,027 条落 `review`。所以"它会开火"本身就是有效的信号。

**2. 但火力的来源不均：编码保留是强信号，化石退化是弱信号（强证据）。**
只洗同源区、保留序列其余部分，比较判定分布：

| 判定 | 真实 | 空模型 | 富集 |
|---|---|---|---|
| `coding_intact`（编码受选择保留） | 34.6% | 3.1% | **11.1×** |
| `distributed_decay`（化石退化） | 45.4% | 37.6% | **1.21×** |
| `assembly_breakpoint` | 13.4% | 13.7% | 0.98× |

即：**"病毒读框明显比同链另两框干净"是与组成无关的真信号；而"分布式退化"里只有约
17% 超出随机组成**（0.454 − 0.376 = 0.078，占真实触发的 17%）。因此凡**依赖退化通道**
的结论（`EVE_suspect` 113 / `ancient_EVE` 15 及其派生的 `MOVE_EVE`）证据强度明显低于
架构与编码保留一侧（`EVE_STRONG_provirus` 4 / `virus_candidate` 2）。对外表述时不要把
这两类混为一谈。

**3. 12 个阈值里 4 个对结论零影响（强证据）。** 逐个阈值扫描 1,027 条输入：
`PROVIRUS_SPAN`（500→4000）、`PROVIRUS_NCOMP`（3/4/5）、`MIN_HEAD_STOPS`（1→4）、
`AA_VIRAL_PID`（80→100）**变更 0 条**；真正决定结论的只有 `MIN_PID`（11.5%）与
`COMP_AA_MIN`（8.2%）。**这意味着"前病毒规模"（`PROVIRUS_*`，v4 那次方向倒置的核心
判据）在当前数据上已不影响任何判定**，与 `multi_locus_arch` 同性质 —— 是信息列，不是
判据。`provirus_scale` 列保留（定义对、测试锁着），但下游不要按它筛。

**4. 组件判定的独立复核：目前未定论（如实记录）。** 用 hmmscan + Pfam-A 独立复核
RT/RH/AP，先做阳性对照（面板参考蛋白标签已知）：RT 检出 **60/60 = 100%**（可用），
AP 3/15 = 20%、RH 0/6 = 0%（**不可用**，属 oracle 自身假阴性）。即可解读的只有 RT：
s2 声称 244 条、oracle 证实 65 条（**26.6%**），且 oracle 另在 **98 条**候选上看到 RT
而 s2 未记。26.6% 是**下界不是估计** —— oracle 只在全长保守参考上验过灵敏度，在 25–45%
一致度的碎片上没测过，所以**不能**据此说"七成 RT 声称是错的"。那 98 条是提高组件召回
率的具体线索。

> **复用派生表的坑**：`EVE_VERIFY_DIR` 指向的产物目录必须由当前代码生成。实测复用旧
> `s2_domains.tsv`（早于 `credit_hsps`）与旧 `locus_architecture.tsv` 会得到
> `EVE_STRONG_provirus` 11 条 / `virus_candidate` 5 条；从原始命中重算是 4 / 2，与本文档
> 记载的 v6.1 数字逐条一致。判别错了方向就是"结论变了"，务必重算。
> 这套手工重算命令现在固化成了工具：`verify/rerun_diff.py`（见下节）。

## v6.2 → v6.3：两处判据修复 + 一条新通道（2026-09-27）

### ① s1 `stop_profile` 的 3nt 系统性高估（已修，本地差分验证）

旧版"最长无终止段"按相邻终止的**起点差**计（`b - a`），把上游终止密码子自身占的 3nt
也算成了无终止，每个间隔高估 3nt（末段同理：`span - pos[-1]` 没扣最后那个终止的 3nt）。
后果集中在 `msf >= 0.25` 这道闸门上：间隔恰好落在 0.25×span 到 0.25×span+3nt 之间的
条目会被错放进 `assembly_breakpoint`（判词是"组装断裂"，实际该走"分布式退化"一侧）。

修法：内部间隔扣 3nt、末段从最后一个终止的**结尾**起算。语义不变、数值系统性变小
（最多 3nt/span）。**影响面用两层验证锁定**：

- 单元测试：新增边界用例（终止紧贴两端、相邻终止间隔恰 3nt），旧断言按正确算术更新；
- A/B 差分（`verify/selftest_fixture.py`，38 条合成 contigs：8 条设计 + 30 条随机，
  旧实现产基线、新实现经 `rerun_diff.py` 重算）：**恰好 1 条按设计翻转**
  （精心构造在 0.25 阈值上的 `assembly_breakpoint_review → ancient_EVE`），
  其余 37 条 decay/verdict 逐条不动——翻转需要某段间隔正落在 150–153nt/600nt
  这个 3nt 宽的窗口里，随机序列撞上的概率极低（30 条随机 contigs 零翻转）。

**真实 1,027 条的影响面要在服务器上量**（本地没有原始命中）：

```bash
python3 verify/rerun_diff.py --run-dir <既有run目录> --evidence <rescue_evidence_scored.tsv>
# 输出 verdict 迁移矩阵 + verdict_diff.tsv; 全对角 = 该 run 与旧版逐条一致
```

### ② 面板升 v3：Badnavirus 从 1 条到 39 条（`build_panel.py`）

v2 面板 222 条里 Badnavirus 只有 1 条 VAP（不在 CANON 五组件里）——该属是科内最大属、
已报道植物 EVE 的主力类群，结构基因通道对它整属失明。根因与修法见上文
"面板覆盖的现实"一节（bare polyprotein 兜底标 AP+RT+RH + `--supplement` 属级配额补充，
既有条目一条不动）。实测 222 → 291 条。

**面板换代与 s1 修复都会改组件命中与 decay 归类**——文档里引用的 v6.1 数字
（`REMOVE` 2 / `MOVE_EVE` 174 / `KEEP_virus` 2 / `REVIEW` 849）是**v2 面板 + 旧 s1**
口径，服务器用 v3 面板复跑后必须整体重算，不能只换局部数字。

### ③ `verify/rerun_diff.py` — 派生表重算差分工具（新增）

把"复用派生表的坑"一节的手工重算命令固化成工具：从 run 目录的**原始命中**
（panel_hits/baits_hits/locus_blastn）用当前代码重算 s1→s4，与目录里的旧 verdict
逐条 diff，输出 `decay_class` 与 `verdict` 两张迁移矩阵 + `verdict_diff.tsv`。
适用场景：判定逻辑任何改动之后对既有 run 量化影响面；迁移矩阵全对角 = 该 run 与
旧版逐条一致。`verify/selftest_fixture.py` 是它的合成数据驱动（同时充当 s1 修复的
A/B 回归：旧实现自动从当前源码生成，基线对比不需要手工维护两份代码）。

### ④ s2c 参考基因组通道（可选 `-G`，2026-09-27）

**解决什么**：s2b 比对的是 OneKP **转录组组装**（碎、无染色体上下文、映射表只覆盖
23 物种），回答"候选是不是寄主基因 / 位点上组件齐不齐"；s2c 比对的是**现成的染色体级
参考基因组**，回答"这条候选在寄主基因组里有没有整合证据"。`EVE_suspect` 一类
contig 级定不了案的候选，缺的正是这层基因座级证据。

**输入与跳过**：清单 `genome_manifest.template.tsv`（三列：4位代码 / 注释 / fasta
路径），有什么基因组填什么；**没有参考基因组的物种不写行**，通道对它们记
`no_ref_genome`（`rg_scope=genome_absent` 或 `no_sample_map`）——跳过不产生任何结论，
**也绝不因此判"真病毒"**。run_all.sh 逐基因组 blastn（megablast，.gz 解压到临时文件、
同 `--subj-cap-gb` 体积闸门），每个基因组比对完追加 `#BLASTED` 标记行，s2c 靠它把
"比过零命中（`no_hit`）"与"根本没比成（`genome_absent`）"分开——两者都不定案。

**两种整合判据**（qualifying 命中：`pident>=90` 且 `aln>=100nt`）：

| rg_call | 条件 | 含义 |
|---|---|---|
| `EVE_flank_confirmed` | contig **两端**各一段宿主命中，同一条染色体同链、中间缺失 `<=500kb`（`flank_gap` 列） | 基因组片段跨着整合位点被组装出来 = 前病毒铁证 |
| `flank_locus_conflict` | 两端都有宿主段但染色体/链/距离对不上 | 更像组装嵌合，只记列不定案 |
| `EVE_splice_chimera` | 内部宿主段呈剪接样式（同染色体相邻段 subject 间隔 `>=500nt`）或与任一端侧翼并存 | 宿主外显子嵌合的转座录 EVE；**contig 级代理，read 级 split-splice 才是金标准** |
| `host_mix_weak` | 内部有宿主段但不满足任何样式 | 仅旁证 |
| `host_dominant` | 自身参考基因组聚合覆盖 `>=80%` 且 `>=90%` 一致、无 MP/CP/AP | 是寄主序列，移交 s3 寄主否决（与 s2b 同一闸门） |

**s3 的单向并入**（可选第 6 参 `refgenome_evidence.tsv`；verdict 表新增
`rg_call` / `rg_scope` 两列，s4 的 filter 表照搬）：证据是**单向**的——只把 review 档
往 EVE 方向升（`EVE_flank_confirmed → EVE_STRONG_provirus`；`EVE_splice_chimera →
EVE_suspect`，不冒充铁证），或把矛盾转 `host_conflict_review` 人工定案（整合证据撞
`virus_candidate` / `host_contamination_likely`；`host_dominant` 撞上游病毒信号）。
`no_hit` / `no_ref_genome` / 无结构 **不改任何 verdict**——参考基因组缺失是常态，
"没比对上"不是真病毒的证据。


## 交接

- **调用方**：跑完 `virome_discovery_pipeline` 之后单独执行本脚本（见上文"运行"），
  产物默认落在 `<DISCOVERY>/08b_EVE_Distinguish/`。
- **下游清除顺序**（数字为 v6.1，1,027 条实测）：先看 `host_contamination_likely`
  （**2**，不再是 48 —— v6 起这一列要求自身物种命中 + 上游三条同源通道全无记录，
  逐个查证后大部分转去了 `host_conflict_review` 33 条，需人工定案而不是直接删）；
  `EVE_LTR_TE` 155 + `EVE_STRONG_provirus` 4 + `ancient_EVE` 15 移出/单列（EVE 合计
  174）；`virus_candidate`（2）是**要保留的真 DNA 病毒**；`virus_fragment_review`（95）
  进 rescue 延伸通道后复判，**不要按 EVE 清掉**。
- **`dna_vs_eve_filter.tsv`**（v5 起由本模块自己的 `s4_filter.py` 产出，不再依赖主管道）：
  每条候选一行，前 11 列原样搬运 s3 的关键列，末两列 `action` / `action_reason` 直接可执行：
  `REMOVE_host_contamination` **2** / `MOVE_EVE` **174** / `KEEP_virus` **2** / `REVIEW` **849**。
  下游按 `action` 列过滤即可，不必再理解 13 种 verdict。
- `EVE_suspect` 113 条的瓶颈仍是 contig 级证据不足（单 ~1kb 保守区无法区分近期 EVE 与
  病毒碎片），升级路径同 v3：cluster 级准基因组组装，或 **v6.3 起直接用 s2c 参考基因组
  通道**（`-G` 给上现成基因组的物种，两端侧翼/外显子嵌合能把其中一批钉死成
  `EVE_STRONG_provirus`；没有参考基因组的物种仍走 rescue 延伸）。
  其中多数**没有寄主位点**（787/1027 条无任何寄主命中），这是当前最主要的证据缺口——
  `host_sample_map.tsv` 只覆盖 23/1000+ 个样本，扩样本映射比调阈值更能提升判别力。
  反过来，能在 `host_sample_map` 里对上的只有 73 条，所以 `host_contamination_likely`
  的天花板本来就很低（≤73 条），别指望它能替代人工复核。
- 依赖：`diamond 2.2.5` 的 `--outfmt` 必须双横线；寄主比对另需 `blastn`+`gunzip`；
  Python 3.8+（PATH 里叫 `python3` 或 `python` 都行，也可用 `PYTHON=/path/to/python` 指定）；
  `host_sample_map.tsv` 是可选的（见上）。
- 测试：四套共 160 个用例，全部可离线跑（端到端用例用的是假 diamond shim，不需要真实外部工具）：
  `tests/test_eve_core.py` 68 个（`eve_scan_core.py` 的纯函数）+ `tests/test_eve_distinguish.py` 19 个
  （`s4_filter.py` 的 verdict→action 映射；`run_all.sh` 的 CLI 校验、外部工具崩掉的
  可见性、五个 python 阶段的文件交接）+ `tests/test_eve_distinguish_logic.py` 45 个
  （s1/s2/s2b/s3 的判定逻辑与边界：读框换算、位点合并、组件并集、位点归组、
  verdict 二维表、不变量）+ `tests/test_eve_refgenome_scan.py` 28 个（v6.3：s2c 的
  侧翼/嵌合/host_dominant 结构判定、no_hit 与 genome_absent 口径分离、s3 rg_adjust
  单向性与矛盾转复核、s4 搬运 rg 列）。表头一律从被测模块 import，不两头各抄一份。
  运行：`cd endogenous_virus_pipeline && python -m pytest tests/ -q`
  （注：不能用 `python -m unittest discover -s tests -t .` 从仓库根跑，tests 目录没有
  `__init__.py` 之外的包结构，根目录 discover 会报 "Start directory is not importable"）。
