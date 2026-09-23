# eve_distinguish — DNA 病毒候选的 EVE 判别 (v5)

对 OneKP 新候选（发现管道的 `10_Reports/eve_candidates.fasta`，1,027 条 = 843 远缘 + 184 新种）
做"结构退化 × 基因座架构"判别，与 `eve_screen.py` 的 EVE 库同源检索**相互独立、互为补充**。

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
# 最常用: 指向一次发现管道的输出目录, 候选/证据表按阶段目录约定自动定位
bash run_all.sh -D /path/to/onekp-virus -A /path/to/1kp/assemblies -t 32

# 显式给候选 fasta 与证据表, 产物落自己指定的目录
bash run_all.sh -q candidates.fasta -e 10_Reports/rescue_evidence_scored.tsv \
                -A /path/to/1kp/assemblies -o out/ -t 32

# 只验证 S1/S2 通道 (不跑寄主比对): host_contamination_likely 检不出, 仅供调试
bash run_all.sh -q candidates.fasta -e evidence.tsv -o out/ --no-host
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
4. `s2b_locus_scan.py`：候选 vs 自身寄主转录组 blastn → 位点归并（`--no-host` 时整表记 unresolved）
5. `s3_verdict.py` 合成判别 **v4 — 二维表 + 寄主同源否决** → `eve_distinguish_verdict.tsv`
6. `s4_filter.py` **v5 新增**：verdict→action 映射 → `dna_vs_eve_filter.tsv`（下游按列过滤）

### 方案一 (s1_decay_scan.py v4)

v3 用绝对计数 `stops>=3` 判退化。问题是**一个 1nt indel 就足以让它之后的全部密码子错读**，
在远缘同源区（25-45% pid）上造出十几个"提前终止"，与真正的化石退化无法区分。

v4 的做法是拿 contig **自己**当基线：同源区同一段序列、同一条链、**另外两个读框**的
终止密度。同一批组装、同一读段、同一碱基组成，组装 indel 对三个读框的污染是同等的，
所以另外两框的终止密度就是这套组装的本地噪声底噪。

- `stop_enrichment = (病毒读框终止数+1) / (另外两框均数+1)`，<1 → `compositional_noise`
  （该读框并不比随机框更脏，终止来自组成/组装，无退化证据）
- ≥1 时再用**分布剖面**区分两种情形：
  - `max_stopfree_frac >= 0.25`：存在一段 ≥1/4 同源区长度的完整无终止区 → `assembly_breakpoint`
    （断裂点两侧各是一段完整可读区，最长无终止段必然占去大半）
  - `head_stops >= 1 或 stops >= 3` 且无长完整段 → `distributed_decay`（化石签名）

  **不能用 "head==0 就推断裂"**：实测那条 92 终止的 6.8kb contig 前端恰好无终止
  （`head_stops=0`、`max_stopfree_frac=0.059`、`enrich=2.16`），显然不是单点断裂。

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
各物种 OneKP 转录组组装，按 (样本, scaffold) 聚合格座，输出 `locus_arch` /
`locus_ncomp` / `locus_hint`。s2 数组件时把融合蛋白（如 `AP+RT`）的两个组件都算上。

> OneKP 转录组组装本身是碎的，一个寄主基因常被打断在多个 scaffold 上（实测 HSP70 候选
> 最高单条比对只有 750nt/94.5%）。所以 s2b 的同源强度用**物种级 query 区间并集覆盖**，
> 不是单条比对覆盖；且只统计**自身物种**（`host_sample_map.tsv` 给定；缺表时按最佳
> bitscore 推断，旁系特别强的样本会判错）。

#### 面板覆盖的现实（读数前先知道）

`panel.fasta` 222 条蛋白，按组件展开计数：MP 62 / CP 62 / RT 60 / VAP 44 / AP 15 / RH 6。
其中 AP 和 RH **绝大多数只以融合体入池**（`AP+RT` 13 条、`RT+RH` 3 条、
`MP+CP+RT+RH` 2 条；纯 AP 只有 2 条、纯 RH 只有 1 条）。后果：`locus_ncomp` 顶不满 4/5
可能只是**面板覆盖问题**而不是生物学事实，别把"组件少"直接读成"没整合完整"。
`build_panel.py` 在组件参考数 <5 时本来就会告警（RH 正好压线）。

### 寄主同源否决 (s3)

HSP70 假阳性的来源：`cd24029` = ASKHA_NBD_HSP70_DnaK_HscA_HscC（伴侣蛋白域，不是病毒域），
`pfam00012`/`cd10233` 同样。管线的 CDD 报告自己带 `hsp70cap=1` 标记，但两条 27-29% 流行率
的 DNA 候选仍凭这个域 PASS_VIRAL 混进新种候选——其一与**自身寄主**的转录组 scaffold
99.8%/2476nt 相同。污染源是 Phycodnaviridae（*Micromonas* 病毒，~197kb NCLDV 基因组）。

否定条件：**自身寄主物种聚合 ≥90% 一致且区间覆盖 ≥80%，且无 MP/CP/AP** → `host_contamination_likely`。
带 MP/CP/AP 的走 EVE 通道（那正是 EVE 问题，不是寄主污染问题）。
`host_xeno` 列另记最强的跨物种命中，用来发现样本串号。

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
| 10 | `--no-host` / 缺 `HOST_MAP` 时 `locus_architecture.tsv` 只有 3 列（重定向只写在三行 `printf` 的最后一行，前几列漏进 stdout） | `{ ...; } > file` 包住整段表头，两种分支都是完整 16 列 |

## 分布对比（1,027 条实测，宿主 23 物种全量比对）

| verdict | v3 | v4 | 含义 |
|---|---|---|---|
| review | 640 | 530 | 证据不足（多为 RNA 病毒科，Cauli 面板不适用） |
| EVE_LTR_TE | 221 | 171 | 转座子科 / 纯 pol 区且 TE 领先 |
| virus_fragment_review | — | 69 | 架构不全 + 编码完好：**真 DNA 病毒本就拼不全，不判 EVE** |
| EVE_suspect | 3 | 71 | 架构不全 + 分布式退化，contig 级无法定案 |
| host_contamination_likely | — | **55** | **自身寄主序列，非病毒**（含两条 HSP70 假阳性） |
| EVE_STRONG_provirus | — | 40 | 单一位点 4/5 组件齐全 + 密码子退化 = 整合前病毒 |
| assembly_breakpoint_review | — | 28 | 终止集中在单个断裂点，救长后复判 |
| virus_candidate | 23 | **22** | 单一位点齐全 + 密码子完好 = 感染中的真病毒（**要留的**） |
| compositional_noise_review | 20 | 20 | 病毒读框不比其他读框更脏，无退化证据 |
| ancient_EVE | 142 | 17 | 纯 pol 区分布式退化 |
| structure_intact_review | 7 | 4 | 非 Cauli 结构完整 |

结构退化类的重新归因（`assembly_breakpoint_review` 55→28、`EVE_STRONG_provirus` 0→40、
`EVE_suspect` 3→71）来自修正 #4/#6：v3 的绝对阈值把组装断裂和化石退化混在一档。

> **宿主比对必须跑完**：`host_contamination_likely` 只在**自身物种**出现在 `HOST_MAP`
> 且该物种的 blastn 已完成时才可能触发。第一版判别只跑完一部分物种（blastn 被 OOM killer
> 打断，见坑 #7），只报出 16 条；23 个物种全量跑完后是 **55 条**（`review`→host 32、
> `structure_intact_review`→3、`virus_fragment_review`→2、`EVE_LTR_TE`→2，另有个别位点完备度
> 变化使 2 条 `EVE_STRONG_provirus` 降为 `EVE_suspect`）。新增 39 条的自身物种聚合是中位
> 99.4% / 0.98 覆盖，而它们的病毒证据全是远缘弱命中（aa 一致度中位 57%、最低 26%，多为
> KOG/pfam 通用域），正是寄主基因被误注释成病毒的样子。跑完 `locus_blastn.tsv` 后应确认
> `locus_errors.log` 里没有 `#NOGZ` / `#TOOBIG` / `#GUNZIP_FAIL`，并核对 v5 加的那行
> "命中行数 / 成功 scaffold 数"不为 0。

## 关键个例（v4 实测）

| contig | v3 | v4 | 依据 |
|---|---|---|---|
| `ERR2040732_..._NODE_11` 6.8kb | ancient_EVE（v1 曾误评"五件套完整"） | `EVE_STRONG_provirus`，hint=`provirus_integrated_signature` | 92 终止 / enrich 2.16 / msf 0.059 + MP+CP+RT+RH 单位点 2417nt |
| `ERR2040803_..._NODE_1874` 2.5kb 流行率 28.56% | 新种候选 PASS_VIRAL | `host_contamination_likely` | 自身寄主 IAJW 聚合 91.2%/0.985，无病毒结构基因（HSP70） |
| `ERR3487378_..._NODE_1805` 2.1kb 流行率 26.92% | 新种候选 PASS_VIRAL | `host_contamination_likely` | 自身寄主 BLAJ 聚合 94.4%/0.915（HSP70，碎片化分布在多个 scaffold） |

## 交接

- **调用方**：跑完 `virome_discovery_pipeline` 之后单独执行本脚本（见上文"运行"），
  产物默认落在 `<DISCOVERY>/08b_EVE_Distinguish/`。
- **下游清除顺序**：先剔 `host_contamination_likely`（55，寄主基因）；
  `EVE_LTR_TE` + `EVE_STRONG_provirus` + `ancient_EVE` 移出/单列（228，EVE）；
  `virus_candidate`（22）是**要保留的真 DNA 病毒**；`virus_fragment_review`（69）进
  rescue 延伸通道后复判，**不要按 EVE 清掉**。
- **`dna_vs_eve_filter.tsv`**（v5 起由本模块自己的 `s4_filter.py` 产出，不再依赖主管道）：
  每条候选一行，前 10 列原样搬运 s3 的关键列，末两列 `action` / `action_reason` 直接可执行：
  `REMOVE_host_contamination` 55 / `MOVE_EVE` 228 / `KEEP_virus` 22 / `REVIEW` 722。
  下游按 `action` 列过滤即可，不必再理解 12 种 verdict。
- `EVE_suspect` 71 条的瓶颈仍是 contig 级证据不足（单 ~1kb 保守区无法区分近期 EVE 与
  病毒碎片），升级路径同 v3：cluster 级准基因组组装或宿主基因组比对。
- 依赖：`diamond 2.2.5` 的 `--outfmt` 必须双横线；寄主比对另需 `blastn`+`gunzip`；
  Python 3.8+（PATH 里叫 `python3` 或 `python` 都行，也可用 `PYTHON=/path/to/python` 指定）；
  `host_sample_map.tsv` 是可选的（见上）。
- 测试：两套共 54 个用例，全部可离线跑（`test_eve_distinguish.py` 里的端到端用例用的是
  假 diamond shim，不需要真实外部工具）：
  `tests/test_eve_core.py` 35 个（`eve_scan_core.py` 的纯函数：清洗、位点合并、
  滑动窗口命令、坐标还原、工具检查、batch 任务）+ `tests/test_eve_distinguish.py` 19 个
  （`s4_filter.py` 的 verdict→action 映射；`run_all.sh` 的 CLI 校验、外部工具崩掉的
  可见性、五个 python 阶段的文件交接）。
  运行：`cd endogenous_virus_pipeline && python -m unittest discover -s tests`
  （注：不能用 `python -m unittest discover -s tests -t .` 从仓库根跑，tests 目录没有
  `__init__.py` 之外的包结构，根目录 discover 会报 "Start directory is not importable"）。
