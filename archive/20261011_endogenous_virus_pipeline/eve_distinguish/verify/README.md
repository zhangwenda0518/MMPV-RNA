# verify/ — 判别器的检验工具（空模型 / 敏感度 / oracle）

所属模块：`eve_distinguish` **v6.2**（v6.2 的产出与 v6.1 逐字节相同，新增的就是本目录）。

这个目录回答一个此前从没被回答过的问题：**这套判据说的是对的吗？**

前面几轮改的都是"代码没做到我们说的"（列错位、断点冻结、清理失效）。但从来没人量过
**"我们说的"本身站不站得住**。四个工具就是补这一步，跑法和实测结论都在下面。

原则三条：

1. **每个工具都要先有自己的阳性对照**。一个没被对照过的检验器给出的"未通过"，
   可能只是检验器自己瞎 —— 本目录就踩过这个坑（见 `oracle_pfam.py`）。
2. **结论按"证据强度"分层写**：实测数、下界、还是未定论，分清楚。
3. 工具本身进仓库、可重跑；不留"我跑过、数字在文档里"的一次性脚本。

---

## 1. `null_model.py` — 空模型特异性（洗牌对照）

把候选序列单核苷酸洗牌（严格保碱基组成、破坏密码子结构与读框连续性），拿洗牌序列
重跑同一套脚本同一套阈值。`P(判据|洗牌)` 就是"只有组成、没有编码结构"时该判据的
触发概率；真实/洗牌 = 富集倍数，接近 1 说明该判据只是组成噪声的代理。

### 1a. 端到端空模型（整条 contig 洗牌，重跑 diamond）

`bash verify/run_null_e2e.sh --mod MOD --run-dir REAL --work WORK --evidence E.tsv`

实测（1,027 条候选，2026-09-24）：

| stage | 真实 | 洗牌 |
|---|---|---|
| diamond panel 命中 | 7,862 行 | **0 行** |
| s2 任一组件 (ncomp≥1) | 409 条 | **0 条** |
| s2 前病毒规模 (provirus_scale) | 14 条 | **0 条** |
| s3 任何实质 verdict | 453 条 | **0 条** |
| s3 `review`（证据不足） | 574 条 | **1,027 条（全部）** |

**结论：整条流水线对组成匹配的随机序列零假阳性。** 洗牌后 blastx 一条命中都没有，
所有下游判据全部落在 `review`。这是这条流水线最硬的一个性质 ——
"它会开火"这件事本身是有意义的信号，不是噪声。

### 1b. 同源区细粒度空模型（只洗同源区，不需 diamond）

`python3 verify/null_model.py region --run-dir REAL --n-shuffle 20 --json OUT`

实测（454 条有可用同源区的候选，各洗牌 20 次）：

| decay_class | 真实 | 真实占比 | 空模型 | 空模型占比 | 富集 |
|---|---|---|---|---|---|
| `coding_intact`（编码受保留） | 157 | 34.6% | 283 | 3.1% | **11.08×** |
| `distributed_decay`（化石退化） | 206 | 45.4% | 3,410 | 37.6% | **1.21×** |
| `assembly_breakpoint` | 61 | 13.4% | 1,241 | 13.7% | 0.98× |
| `compositional_noise` | 30 | 6.6% | 4,146 | 45.7% | 0.14× |

**结论（本次最重要的方法学结论）**：

- `coding_intact` 有 **11 倍**富集 —— "病毒读框明显比同链另两框干净"是一个真实的、
  与组成无关的信号。这一路判据是可信的。
- `distributed_decay` 只有 **1.21 倍**富集 —— 化石退化签名里扣除组成噪声之后只剩
  **约 17%**（0.454 − 0.376 = 0.078，占真实触发的 17%）是超出随机的部分。
  换句话说：**当初级判据把一条候选标成"分布式退化"时，它更可能是在描述这段序列的
  碱基组成，而不是在描述"这是病毒化石"。**
- `assembly_breakpoint` 富集 0.98 —— 与随机不可区分，符合预期（洗牌会一并破坏
  组装 indel 结构，所以这个类别本身不该有富集）。

这条结论直接改变了对外口径：`EVE_suspect` / `ancient_EVE` / `MOVE_EVE` 这些**依赖退化
通道**的结论，其证据强度远低于 `coding_intact` 一侧。**"有病毒同源（架构/编码保留）"
是强信号；"像化石一样退化"是弱信号。**

## 2. `sensitivity.py` — 阈值敏感度扫描

对每个 (脚本, 常数, 取值) 组合把常数注入后重跑该 stage 及下游，统计相对基线的
verdict 变更条数。**变更 0 说明这个阈值在当前数据上不是判别参数，而是装饰。**

`python3 verify/sensitivity.py --run-dir REAL --evidence E.tsv --work WORK`

实测（1,027 条，基线 `review 574 / EVE_LTR_TE 156 / EVE_suspect 117 / ...`）：

| 阈值 | 扫描范围 | 最大变更 | 占 1,027 |
|---|---|---|---|
| `MIN_PID`（blastx 一致度门槛） | 20 → 45 | **118 条** | **11.5%** |
| `COMP_AA_MIN`（组件记功长度闸） | 50 → 300 | **84 条** | **8.2%** |
| `MIN_SPAN`（s1 基线最小同源区） | 60 → 300 | 42 条 | 4.1% |
| `MIN_ALN`（比对长度门槛） | 30 → 150 | 28 条 | 2.7% |
| `CLEAN_MAX_STOPS` | 2 → 20 | 5 条 | 0.5% |
| `CLEAN_ENRICH` | 0.3 → 1.0 | 4 条 | 0.4% |
| `MIN_BITS` | 30 → 80 | 2 条 | 0.2% |
| `MIN_DECAY_STOPS` | 2 → 8 | 2 条 | 0.2% |
| **`PROVIRUS_NCOMP`** | 3 / 4 / 5 | **0 条** | 0% |
| **`PROVIRUS_SPAN`** | 500 → 4,000 | **0 条** | 0% |
| **`MIN_HEAD_STOPS`** | 1 → 4 | **0 条** | 0% |
| **`AA_VIRAL_PID`** | 80 → 100 | **0 条** | 0% |

**结论**：

- 12 个阈值里 **4 个完全不动结论**。其中两个是**整套"前病毒规模"定义**
  （`PROVIRUS_SPAN` / `PROVIRUS_NCOMP`）—— 那是 v4 那次"方向倒置"的核心判据。
  v6 去掉了 `locus_full or provirus_scale` 的取或之后，`provirus_scale` 已经不再
  影响任何 verdict，**这一列现在是信息列，不是判据**（与 `multi_locus_arch` 同性质）。
- `MIN_HEAD_STOPS = 0` 影响，说明"分布式退化"的定性实际上由 `enrich ≥ 1` 单独决定，
  头部分布那套逻辑没起作用 —— 与 1b 的结论一致（退化通道整体偏弱）。
- 真正在决定结论的只有 **`MIN_PID`** 和 **`COMP_AA_MIN`**（后者是 v6 加的，8.2%）。
  这两个是需要优先校准的参数；其余阈值可以放宽或收紧而结论不变。

## 3. `oracle_pfam.py` — 组件判定的独立 oracle（**目前未定论**）

想法：`components` 列的 MP/CP/AP/RT/RH 来自**面板条目名**，整条链路没有任何一步看过
query 序列本身有没有对应结构域。所以用 hmmscan + Pfam-A 独立复核 RT/RH/AP
（MP/CP 无可靠单家族映射，不判定）。

```bash
python3 verify/oracle_pfam.py emit --run-dir REAL -o pep.fa
hmmscan --domtblout pfam.domtbl -E 1e-3 --cpu 32 Pfam-A.hmm pep.fa > /dev/null
python3 verify/oracle_pfam.py control --panel panel.fa --domtbl panel_pfam.domtbl  # 先看这个
python3 verify/oracle_pfam.py compare --s2 REAL/s2_domains.tsv --domtbl pfam.domtbl
```

### 阳性对照（必做，否则 compare 的结论是假的）

拿面板参考蛋白做阳性样本 —— 它们条目名里就写着 RT/RH/AP：

| 组件 | oracle 在**有该标注**的面板蛋白上的检出率 | 可用性 |
|---|---|---|
| RT | **60 / 60 = 100%** | 可用 |
| AP | 3 / 15 = 20% | **不可用** |
| RH | 0 / 6 = 0% | **不可用** |

**所以这个 oracle 只有 RT 那一列能读**，AP/RH 的"未检出"是 oracle 自己的假阴性。

### 候选侧结果（只有 RT 可解读）

| 组件 | s2 声称 | oracle 证实 | 一致率 | oracle 另发现 |
|---|---|---|---|---|
| RT | 244 条 | 65 条 | **26.6%** | 98 条 |
| AP | 74 条 | 4 条 | 5.4% | 15 条（不可解读） |
| RH | 21 条 | 0 条 | 0% | 0 条（不可解读） |

**结论：定不了案，但要记下来两件事。**

1. RT 的阳性一致率 26.6%。它是个**下界不是估计** —— oracle 对照只在"全长、保守的
   参考蛋白"上验到 100%，在"25–45% 一致度的碎片 query"上的灵敏度没测过，必然偏低。
   所以**不能**据此说"73% 的 RT 声称是错的"（我差一点就这么写了）。
2. 反过来 oracle 在 **98 条候选**上看到 RT 而 s2 没记 —— 这是**可执行的具体线索**：
   要么这些 HSP 被 `MIN_ALN/MIN_PID` 挡了，要么打中的面板条目名没带 RT。要提升组件
   通道的召回率，从这里查。

要把它变成能裁决的工具，需要换 Caulimoviridae 专属家族、用 Pfam gathering threshold
（`--cut_ga`）、或改用 CDD（服务器上有 `rpsblast`）。

---

## 4. 不变量检查 → 已进 `tests/test_eve_invariants.py`

原先散在仓库根目录的一次性核查脚本，现在是正式测试：

```bash
python -m pytest tests/test_eve_invariants.py -q                       # 谓词层, 离线
EVE_VERIFY_DIR=/path/to/run python -m pytest tests/test_eve_invariants.py -q   # 真实产物层
```

真实产物层实测 **21/21 通过**（1,027 条候选的跑完产物）。锁住的不变量：
`locus_full` 必须挂在真实寄主位点上、不许出现伪位点并组、`locus_ncomp` 与
`locus_comps` 必须一致、`components` 只许出现 CANON 整词、`loci_detail` 不许漏出
set repr、`provirus_scale`/`multi_locus_arch` 必须满足自身定义、
只有 `host_contamination_likely` 允许 REMOVE、action 表与 verdict 表逐行对齐。

### ⚠️ 复用派生表的坑（本次真踩到）

`EVE_VERIFY_DIR` 指向的目录**必须是用当前代码生成的**。实测：复用上一次运行留下的
`s2_domains.tsv`（17 列，早于 `credit_hsps` 那次修改）和 `locus_architecture.tsv`
（旧 s2b 的 `locus_ncomp` 偏高），得到的 verdict 是 `EVE_STRONG_provirus` 11 条、
`virus_candidate` 5 条；而**从原始命中重算**（`q.fa` + `panel_hits.tsv` + `baits_hits.tsv`
+ `locus_blastn.tsv`）得到 4 条 / 2 条，与文档记录的 v6.1 数字**逐条一致**。
判别错了方向就是"结论变了"，所以宁可多花几分钟重算：

```bash
python3 s1_decay_scan.py q.fa all_hits.tsv s1_decay.tsv
python3 s2_domain_scan.py q.fa panel_hits.tsv baits_hits.tsv s2_domains.tsv
python3 s2b_locus_scan.py locus_blastn.tsv s2_domains.tsv s1_decay.tsv locus_architecture.tsv host_sample_map.tsv
python3 s3_verdict.py s1_decay.tsv s2_domains.tsv locus_architecture.tsv rescue_evidence_scored.tsv verdict.tsv
python3 s4_filter.py verdict.tsv -o dna_vs_eve_filter.tsv
```

（`s1`/`s2` 是纯计算，`s2b` 复用缓存 `locus_blastn.tsv` 也是纯计算 —— 只有 blastx 与
host blastn 贵，那两份原始命中不用重跑。）

### rerun_diff.py — 上面那套手工命令的工具化（v6.3 新增）

```bash
python3 verify/rerun_diff.py --run-dir <既有run目录> --evidence <rescue_evidence_scored.tsv>
```

从 run 目录的**原始命中**用当前代码重算 s1→s4（产物落 `<run-dir>/rerun_currentcode/`，
不覆盖原表），与目录里的旧 `eve_distinguish_verdict.tsv` 逐条 diff：stdout 打
`decay_class` 与 `verdict` 两张迁移矩阵（旧行 → 新列），`verdict_diff.tsv` 只列变了的行，
退出码 1 = 有 verdict 变更。**适用场景：判定逻辑任何改动之后，对既有 run 一条命令量化
影响面**；迁移矩阵全对角 = 该 run 与旧版逐条一致。v6.3 的 s1 `stop_profile` 修复就是
先用它做了合成数据 A/B（见下），真实 1,027 条的复算也走这条命令。

### selftest_fixture.py — 合成数据 A/B 驱动（v6.3 新增）

```bash
python3 verify/selftest_fixture.py [输出目录]     # 默认 verify/_selftest_fixture/
```

38 条合成 contigs（8 条设计：coding_intact / distributed_decay / assembly_breakpoint /
0.25 阈值边界翻转 / compositional_noise / 无命中 / TE，+30 条随机）自动构造完整 run
目录，基线 verdict 用**旧版 stop_profile**（从当前 s1 源码文本替换生成，不需要手工维护
第二份代码）产出，再交给 `rerun_diff.py` 用当前代码重算对比，最后按设计意图逐条断言
decay/verdict。s1 判定算术再被改动时，这个脚本就是第一道回归。


---

## 一句话总结四个工具各自证明了什么

| 工具 | 结论 | 强度 |
|---|---|---|
| 端到端空模型 | 对组成匹配的随机序列零假阳性（0/1027 触发任何实质判定） | 强，可直接写进方法 |
| 同源区空模型 | `coding_intact` 11× 真信号；`distributed_decay` 仅 1.21× | 强，且**改变了对外口径** |
| 敏感度扫描 | 12 个阈值里 4 个不动结论；真正决定结论的只有 `MIN_PID` 与 `COMP_AA_MIN` | 强 |
| Pfam oracle | RT 一致率 26.6%（下界）、98 条候选 oracle 有而 s2 无；AP/RH 不可裁决 | **未定论**，已记录线索 |

---

## 校准层自身的验证（谁看住这些工具？）

上面四条结论都出自这套工具。**工具自己错了，结论就全错** —— 而这已经发生过三次（见下）。
所以每个工具都必须先过自己的对照，再谈它的结论。测试在 `tests/test_eve_verify_tools.py`
（19 个用例：17 离线 + 2 需 `EVE_VERIFY_RUN_DIR`）。

### 已经真实发生过的三类失效

| 失效 | 表现 | 只有什么能抓住 |
|---|---|---|
| **仪器瞎** | oracle 初版报"88% 组件声称未获证实"，看着像重大发现；真因是它自己键提取错 | **阳性对照**：拿条目名明写 RT 的面板参考蛋白去测 → RT 检出 0/60（不可能） |
| **干预没生效** | `COMP_AA_MIN` 初版报"影响 0%"，真因是默认参数在 def 时求值，改全局到不了代码路径 | **下游诊断量**：原自检断言"全局被改了"——它**通过**了，而结论依然错 |
| **静默退化** | `di` 洗牌器欧拉路径写错，走不完边就悄悄退回 mono（长度对、性质错） | **逐次试验断言性质**：60 次试验里 25 次不保二核苷酸 |

### ① 已知答案对照

```
python3 verify/null_model.py control --work /tmp/ctl          # 离线, 秒级
python3 verify/sensitivity.py ... --control                   # 需真实 run-dir
```

**null_model**：自己造两类答案由构造保证的输入（frame 1 无终止的编码序列 / 纯随机），
跑同一套判据看是否分得开。实测：

| 输入 | `coding_intact` 真实 | 空模型 | 富集 | 要求 |
|---|---|---|---|---|
| 已知编码 (orf) | 1.000 | 0.007 | **133×** | 必须 ≥3× 且真实占比 ≥0.5 |
| 已知随机 (random) | 0.000 | 0.005 | 0.00× | 必须 ≤3× |

**sensitivity**：用**不可能值**证明注入真的到达代码路径 —— 这是把"阈值在宽范围内不敏感"
（真发现）和"我的注入静默失败了"（假发现）分开的唯一办法：

| 注入 | 期望的下游必然结果 | 实测 |
|---|---|---|
| `MIN_PID = 101`（一致度上限 100） | 有组件的候选必须归零 | 归零 ✅ |
| `PROVIRUS_NCOMP = 6`（CANON 只有 5 个） | `provirus_scale` 必须全 FALSE | 11 → 0 ✅ |

`PROVIRUS_NCOMP=6` 那一针**顺手把一个可疑结论变成了已验证结论**：provirus_scale 归零
而 verdict 一条不变 → 独立印证"该列已不参与判定"。

### ② 干预生效性：看下游量，不看变量本身

`DIAGNOSTICS` 给每个参数配一个**必须随之移动的下游量**（`COMP_AA_MIN`→ncomp 分布、
`PROVIRUS_*`→provirus TRUE 数、`MIN_BITS`→TE flag 数、s1 系→decay_class 分布…）。
参数动而下游量一动不动时，再打一针不可能值，按结果三态判定：

- 诊断量在网格内**动** → 阈值在此范围有效或不敏感（结论可用）
- 网格内不动、**打针会动** → 该阈值在范围内不敏感（合法结论）
- **连不可能值都不动** → 该阈值**逻辑冗余**，或仪器坏（由 ① 的实验排除后者）

三处**诊断映射本身写错**也在这里被抓出来：`MIN_BITS` 作用在 baits 路径却被拿 panel 组件
当诊断量；`MIN_SPAN` 的诊断量固定写死 90（是常数，永远不动）。修正后才得到下表。

### ③ 零模型稳健性（决定论文口径）

三种零模型，保住的东西依次变多：

| 零模型 | 保住 | 破坏 |
|---|---|---|
| `mono` | 碱基总组成 | 密码子结构、位置特异组成 |
| `frame` | 总组成 **+ 位置特异组成** | 密码子顺序 |
| `di` | **全部二核苷酸频率**（欧拉路径法） | 更高阶结构 |

实测（1,027 条，454 条有可用同源区，各洗牌 20 次；真实占比与零模型无关，只有空模型率变）：

| 判定 | mono 富集 | frame 富集 | di 富集 | 结论 |
|---|---|---|---|---|
| `coding_intact` | 11.08× | 5.91× | 8.94× | **方向稳健**，幅度须报区间（≈6–11×） |
| `distributed_decay` | **1.21×** | **1.18×** | **1.26×** | **1.18–1.26×，零模型无关** |
| `assembly_breakpoint` | 0.98× | 1.22× | 0.91× | ~1，无信号 |

**这是 ③ 最重要的一条**：`distributed_decay` 的"只比随机高一点"在三种零模型下都成立
（1.18–1.26×），扣除噪声后超出随机的比例稳定在 **16–20%**。而 `coding_intact` 的**方向**
稳健（5.9–11.1×，始终强富集），但**幅度依赖零模型**，所以对外只能报区间并注明零模型。

### ④ 独立实现差分

线上实现 vs 暴力实现，用**完全不同的机制**重算同一件事：

| 函数 | 线上机制 | 暴力机制 | 结果 |
|---|---|---|---|
| `merge_loci` | 排序后单趟扩展 | 反复合并到不动点 | 300 组随机输入 + **全部真实候选** 一致 |
| `credit_components` | 区间减法 `subtract_ivs` | 占用位置**集合**取差集 | 300 组随机输入 + **全部真实候选** 一致 |
| `classify_decay` | 代码 | docstring 当规格重写 + 枚举网格 | 1,400 个网格点逐点一致 |

这类差分才是"我读过代码、看着对"的替代品 —— 同一次验证里最强的证据一直是这个形态
（坐标链当初靠纯 Python 切片 + `seqkit subseq` 两条独立参照逐条吻合）。

### 修正一条此前写错的结论

原来报的"`AA_VIRAL_PID` 80→100 变更 0 条 → 该阈值是装饰"**是错的**。它只作用于 s3 的
寄主否决分支，而那次扫描是 `--no-host`（locus 表仅表头），该分支**结构上不可能触发** ——
测不了不等于没用。现在它被显式标成"不可测"（`BLIND_NO_HOST`），不再计入"零影响"名单。
同理，`MIN_HEAD_STOPS` 的"0 条变更"经不可能值探针确认是**逻辑冗余**（同一 OR 分支里
`stops >= MIN_DECAY_STOPS` 恒先满足），这是关于**代码结构**的事实，不是关于数据的。

### 校准层的边界（不能被上述任何一项消除）

以上全部只能**界定误差率，不能确立正确性**：

| 能校准 | 不能校准 |
|---|---|
| 对组成匹配随机序列的特异性 ✅ | 整条流水线的灵敏度/召回 ❌ |
| 阈值敏感度、哪些阈值逻辑冗余 ✅ | 准确率 ❌ |
| 与正交方法的一致率（仅 RT）✅ | 其余组件的一致率 ❌ |

端到端空模型破坏了全部信号，**只能测特异性**。要测灵敏度需要一份**策展阳性集**
（文献已报道的植物 EPRV + 策展现存病毒 + 带病毒样结构域的寄主基因作阴性）——
这是唯一还缺的一环，需要外部策展数据。
