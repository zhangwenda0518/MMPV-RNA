# MMPV AI 驱动化设计与论证 — 四层递进架构

> 2026-10-08 成稿。回答一个问题: **MMPV 如何从"人驱动、AI 点缀"升级为"AI 驱动"**,
> 且不牺牲生信平台的生命线——确定性、可复现、可审计。
>
> 本文是论证 + 交付清单 + 上手指南。代码已全部落地并有单测守门
> (mmpv_common/tests/, 94 passed, 全离线可跑)。
> **落地顺序四阶段已全部完成** (2026-10-08): ① ai.py+digest+五报告接线 →
> ② triage+advisor → ③ copilot → ④ stage_tools+MCP+受控 agent。

---

## 0. 结论先行

MMPV 已有的资产——`--stage` 确定性阶段 API、检查点断点续传、退出码/blocked 语义
(`phylo_pipeline.py:325`)、五边界交接契约 (`io_layout.py`)、profile 化配置——
恰好是 AI Agent 最难自建的东西: **可安全调用的工具面**。

因此路线不是"用 LLM 重写分析逻辑", 而是:

> **确定性内核不动, 在外面套四层 AI 外壳: 状态数字化 → 解读 → 决策辅助 → 编排。
> AI 只读摘要、只写标注、只产校验过的命令草稿; 执行权永远在人的终端。**

四层每层独立可用、独立交付价值, 随时可以停在任意深度。

---

## 1. 目标与非目标

**目标**
1. 人从"记路径、翻报告、逐参数拼命令"中解放出来: 问一句就知道跑到哪、下一步跑什么;
2. 大规模候选 triage (工具投票冲突的 contig 等) 有 AI 辅助标注, 人只裁决 REVIEW;
3. 任意 AI 客户端 (ZCode / Claude / Cursor / MCP 生态) 能安全地"看见"MMPV;
4. 全部 AI 参与可溯源、可复现、可关掉 (离线时平台完整可用)。

**非目标 (红线, 见 §4)**
- 不用 LLM 替代任何统计算法 (分类投票 / CheckV / DESeq2 …);
- 不给 AI 执行管线命令的能力 (copilot 与 MCP 均无执行工具);
- 不让 AI 直接修改任何数据文件的数据列。

---

## 2. 为什么是"四层递进" — 方案论证

### 2.1 三个候选方案的否决理由

| 方案 | 内容 | 否决理由 |
|---|---|---|
| A. 全 Agent 化 | LLM 直接编排各 bioconda 工具 (类似把 MMPV 降级为提示词集) | 不可复现: 同一目标两次跑出不同命令链; 成本爆炸: 每步工具调用都是 LLM 决策; 无法过审: 论文方法章没法写"AI 当时决定用 rnaviralSPAdes" |
| B. AI 内嵌进管线 | 把 LLM 调用写进各编排器的关键判定点 (如 filter 阈值由 AI 定) | 侵入式改动 7 条管线; AI 故障 = 管线故障; 阈值类决策一旦 AI 化, 历史结果失去可比性 |
| C. **外挂分层 (采纳)** | 管线零改动; 新增 mmpv_common AI 层 + 统一入口; AI 消费机器可读状态, 产出建议/标注/校验草稿 | 渐进、可回退、离线完整可用; AI 故障不传染; 论文流程可用纯规则参数重放 |

### 2.2 分层的依据: "AI 能可靠做什么"决定它站哪层

| 层 | AI 任务 | 错误代价 | 对策 |
|---|---|---|---|
| L1 状态 | 不用 AI (纯规则抽取) | — | 确定性代码 |
| L2 解读 | 总结/解释 (读 digest 写人话) | 低 (说错话, 不改数据) | 提示词明令"数字只准引用 digest"; 全量溯源日志 |
| L3 决策辅助 | 标注 KEEP/REVIEW + 串建议成路径 | 中 (漏标) | 动作集合封闭 {KEEP, REVIEW}; 拿不准一律 REVIEW; 只写 ai_* 列 |
| L4 编排 | 自然语言→命令草稿 | 高 (跑错命令) | 脚本/参数/元字符三道白名单; 不执行; 危险 flag 强制确认 |

错误代价越高, AI 自由度越小——这是整个架构的第一性原理。

### 2.3 为什么必须先做 L1

现状下 AI (和人) 只能"翻 HTML 报告"了解运行状态——非结构化、贵、且每管线格式
不同。`run_digest` 把它变成一份 200 行内的 JSON 后: L2 提示词有了稳定输入,
L3 规则有了判定依据, L4 意图规划有了项目上下文。**没有 L1, 上层全是空中楼阁;
有了 L1, 上层每层都是薄封装。**

---

## 3. 分层设计 (已交付)

```
                         ┌─────────────────────────────────────┐
  人 / AI 客户端          │  mmpv_ai.py (统一 CLI)  mcp_server   │  L4 编排
                         │  chat / validate │ MCP: 5 tools      │  (只校验不执行)
                         └──────────────┬──────────────────────┘
                                        │
                 ┌──────────────────────┼──────────────────────┐
                 ▼                      ▼                      ▼
        ┌────────────────┐    ┌──────────────────┐    ┌──────────────────┐
        │ copilot.py     │    │ advisor.py       │    │ triage.py        │
        │ 白名单三道防线   │    │ 规则骨架+AI 叙事  │    │ 只加 ai_* 三列    │  L3 决策辅助
        └───────┬────────┘    └────────┬─────────┘    └────────┬─────────┘
                │                      │ digest                 │ 证据 TSV
                │                      ▼                        │
                │            ┌──────────────────┐              │
                │            │ run_digest.py    │              │                L1 状态
                │            │ 90_Handoff/state/│              │
                │            │  run_digest.json │              │
                │            └──────────────────┘              │
                │                      │
        ┌───────┴──────────────────────┴───────────────────────┘
        ▼
┌───────────────────────────────────────────────────────────────┐
│ ai.py — provider 抽象 (OpenAI 兼容) · 版本化提示词 PROMPTS      │  L2 解读
│ 溯源日志 ai_log.jsonl · 离线降级 (无 key 一切照常)              │
└───────────────────────────────────────────────────────────────┘
        ▲
┌────────┴──────────────────────────────────────────────────────┐
│ stage_cards.py — 8 管线 × stage 卡片 (单一真相源, 喂文档/AI/人) │  L1 知识
└───────────────────────────────────────────────────────────────┘
```

### 3.1 L1 状态层 `run_digest.py`

- 扫 standard/legacy 双布局项目根 → `90_Handoff/state/run_digest.json`
  (schema_version=1, 字段只加不改);
- 每管线: 各阶段目录存在性/文件数/字节 + 关键计数 (Core14 行数、clean reads 数、
  组装样本数、kingdom_summary 行数、报告 HTML 清单…);
- 单文件/单管线解析失败记 `error` 不拖垮整体; 目录缺失记 `absent` (项目是渐进跑的);
- `compact_for_ai()` 出瘦身版供提示词 (剥子目录名清单, AI 永不见原始数据路径以外的
  大对象)。

### 3.2 L1 知识层 `stage_cards.py`

8 条管线 × 每个 stage 的 `name/desc/inputs/outputs/risks/requires` 卡片,
内容取自各编排器 STAGE_ORDER 与 doc/*.md。`self_check()` 断言结构完整;
`summary_text()` 供 AI 上下文。一份真相源同时喂 `mmpv_ai.py cards`、AI 提示词
与后续文档生成——防"文档说 A 代码做 B"。

### 3.3 L2 解读层 `ai.py`

- **配置**: `MMPV_AI_API_KEY` > `DEEPSEEK_API_KEY` (与 gsa_sra.info.py 同口径),
  `MMPV_AI_BASE_URL`/`MMPV_AI_MODEL` 任意 OpenAI 兼容端点 (DeepSeek/Ollama/vLLM);
- **提示词注册表**: `register_prompt(name, version, template)`; 同名同版本改内容
  直接报错 (必须升版本); 渲染体带 sha256 入日志——同版本同输入必得同 prompt;
- **溯源**: 每次调用追加 `90_Handoff/state/ai_log.jsonl` (ts/op/model/prompt 元
  数据/usage/latency/输入哈希)——审稿人问"AI 参与了什么", 拿日志回答;
- **离线降级**: 所有高层接口 (`summarize_run` / `ai_interpretation_section` /
  `ask`) 无 key 时返回 `status="unavailable"` 结构, 绝不抛异常——**AI 是加速器,
  不是依赖**。

**报告接线 (已落地)**: 五份报告 (discovery/preprocess/eve/host_expression/
submission；其中 **eve / host_expression 报告随 ⑥⑧ 于 2026-10-11 归档**, 接线代码保留但不再产出) 的 `main()` 均已在生成 HTML 后调用统一接线点, 形如:

```python
try:
    from mmpv_common.ai import append_ai_section_to_html
    from mmpv_common.run_digest import guess_state_dir
    append_ai_section_to_html(html_path, stats_dict,
                              state_dir=guess_state_dir(outdir))
except Exception:
    pass
```

契约: 离线/失败时**不改动文件** (status != ok → 不写回), 报告本体零风险;
AI 有话时在 </body> 前插入一段带模型名与 prompt 版本的解读卡片。
discovery 报告既有的 `--ai-summary` 旁路总结保留, 两者并存互不影响。

### 3.4 L3 决策层 `advisor.py` + `triage.py`

**advisor**: 9 条确定性规则 (空项目警告 / 有原始数据没质控 / clean 没 deplete /
深度极差 ≥10× 提示 bbnorm / 样本极少提示跳过共组装 / 鉴定完没跑下游 / ④ 差报告 /
宿主基因组就绪提示 EVE【已归档 2026-10-11, 该规则已停用 return None】/ rescue 就绪提示提交), 每条建议自带 `command` 草稿——
草稿仍要过 L4 白名单才算数。`advise_with_ai()` 让 LLM 只写 ≤200 字"执行路径叙事",
提示词明令**不得新增建议、数字只准引用给定 JSON**。

**triage**: 把 `eve_distinguish/s4_filter.py` 的 verdict→action 模式推广为通用
标注器。宪法 (tests/test_ai_driven.py 有不变量断言):

1. 只追加 `ai_action / ai_reason / ai_source` 三列, 原有列逐字节原样回写;
2. 动作集合封闭 {KEEP, REVIEW}——REMOVE 类动作只属于 eve_distinguish 的规则链;
3. 输入已含 ai_* 列 → 拒绝写入 (防静默覆盖);
4. AI 模式: batch + JSON schema 强制, 漏标/坏 schema/不可解析 → 保守回退 REVIEW。

### 3.5 L4 编排层 `copilot.py` + `mcp_server.py`

**copilot 三道防线** (全部测试覆盖):
1. 脚本白名单: 只认 COMMAND_SPECS 里 8 个编排器 (flag 集合取自各 argparse, 2026-10 核对);
2. 参数白名单: flag 必须在允许集; `--stage` 词元必须在真实 choices 里; 条件必填
   (如 stage 含 clean → 必须 `--input_reads`) 与参数组 (eve: `-g` 或 `-B`) 校验;
3. 元字符零容忍: `; | & \` > < $( )` 换行 → 直接拒绝, 不进入解析。

另有: 危险 flag (`--force`) → `requires_confirmation` + 提醒; 内联 API key →
提醒改环境变量。**没有 execute()——这是特性, 不是缺陷。**

**过程报告层 `process_report.py` (2026-10-08 增补)**: 让最终报告承载**整个过程**
而不只读数据。数据源 = `digest_history.jsonl` (每次 digest 落盘追加的时间线快照) +
`ai_log.jsonl` 聚合 (AI 参与清单: 调用数/操作分布/模型/提示词版本/token) +
`agent_runs/` 留痕。三件交付:

1. `build_process_panel_html()` — **确定性**过程面板 (管线 chips ①→⑧、时间线表、
   AI 参与统计、治理声明), 零 AI 依赖, 离线报告也能展示全过程骨架;
2. `append_process_panel_to_html()` — 五份管线报告在生成后自动嵌入该面板
   (幂等: PANEL_ID 查重; 项目外运行优雅放弃);
3. `generate_master_report()` — 项目级主报告 `90_Handoff/MMPV_master_report.html`
   (`mmpv_ai.py report`): 管线全景 + 过程时间线 + 过程面板 + 规则建议 +
   AI 过程叙事 (新提示词 process_narrative: 推进程度/AI 参与/下一步三段,
   无 AI 记录时如实写"纯规则运行")。

**triage 实战化 (P2, 2026-10-08 完成)**: `mmpv_ai.py triage --preset ident --project <根>`
从 ②鉴定目录现场合并各工具 `*.result.id` → contig×工具 1/0 证据矩阵 (全票分母
= 实际在跑的工具, 零命中 contig 不进矩阵), 再走同一套 triage 宪法; `--preset
eve_distinguish` 定位 `dna_vs_eve_filter.tsv` 自动标注。REVIEW 计数由过程层扫描
`*.triage.tsv` 汇入过程面板与主报告 (KPI 卡 + 专区分区)。

**报告互链与演化图 (P5, 完成)**: 主报告顶部导航直达五管线报告 (只链真实存在
的文件); 过程面板反向链接主报告; digest 演化图 = 内联 SVG 阶梯折线 (在线管线数
× 快照序列, 零依赖)。**AI 用量面板**: ai_log 按 op / 按月聚合的调用与 token 表。

**provenance 一键导出 (P3, 完成)**: `mmpv_ai.py provenance <根>` →
`90_Handoff/state/provenance_export.json` (机器可读全量: ai_log 事件逐条 +
agent 轮次 + triage + 时间线 + 治理声明) 与 `90_Handoff/AI_PROVENANCE.md`
(人读方法学文档, 含**可直接粘贴的中英文方法章段落** —— 有 AI 参与时如实写
用途/模型/次数/留痕位置, 纯规则运行时如实写"未使用 LLM")。导出本身零 AI 调用。

**事后补 AI (`mmpv_ai.py enrich`, 2026-10-08 增补)**: 无 key 时管线照跑、报告照出
(AI 段缺席但文件完好); 之后配好 key 跑 `enrich` 一次性补齐——扫项目内已生成的
管线报告 HTML, 补插 AI 解读段 (id=mmpv-ai-section, 幂等可反复跑; 段内标注
"事后补解读 · 数据来源: run digest") + 过程面板 + 重新生成带叙事的主报告。
解读数据源 = 项目 run digest (事后唯一可信的机器状态); `compact_for_ai` 已剔除
generated_at 等易变字段, 同项目重复补叙事命中缓存零计费。残缺 HTML (无
</body>) 如实报 error 不改动。另: ④ 分析报告内建的浏览器端 AI 面板 (key 填在
页面里) 是查看时补充的老路径, 与 enrich 互补。

**AI 调用缓存** (同期增补): `chat(cache_path=...)` 以 (模型/端点/温度/json_mode/
messages) 哈希为键落 `state/ai_cache.json` —— 报告重渲染零计费且文案稳定;
`MMPV_AI_CACHE=0` 关闭; 命中在 ai_log 记 `cached: true`。

**MCP server**: 手写 stdio JSON-RPC (零第三方依赖, 与 io_layout 同级"可安全
import"), 7 个工具: digest / advise / cards / validate_command / stage_tools /
plan_chain / ask。刻意不提供 run_pipeline 类执行工具 (§4 红线 4)。接入见 §7。

**stages 工具化 `stage_tools.py`**: 8 个编排器 → 带 JSON Schema 的工具描述
(参数面取自 copilot 白名单; 前置条件 = digest 谓词, 即交接契约的机器可查形态;
后置条件 = 跑完后 digest 应出现的产物键; 退出码语义 0/1/2 沿 phylo 约定)。
`plan_chain()` 把 advisor 建议升格为按主链 ①→⑧ 排序的计划, `<项目>` 自动
替换, 其余 `<占位符>` 显式列为 missing_inputs —— 路径不臆造。

**agent 全链路 `agent.py` (受控循环)**: plan → (确认门) → execute →
read_digest → iterate。与"全自动 Agent"的边界写在实现里:

- 计划只来自规则 (advisor + plan_chain), AI 不决定跑什么;
- 每步命令先过 copilot 三道防线, 以 token 列表直 subprocess (不经 shell);
- 危险命令 (`--force` 类) 即使 `--yes` 也必须交互确认, 非交互一律跳过;
- 缺输入的步骤显式 skipped (`--set "占位符=值"` 补齐), 预算 `--max-steps` 封顶;
- 退出码 1=failed 停链, 2=blocked 跳过; 全程留痕 `90_Handoff/state/agent_runs/`。

---

## 4. 红线 (AI 平台宪法)

1. **确定性内核不可替代**: 分类=8 工具加权投票, 完整性=CheckV, 差异=DESeq2/TEA-GCN。
   AI 可以解读、质疑、标 REVIEW, 但不能改写任何算法结论;
2. **AI 永远不执行破坏性操作**: 只产出建议+理由, `--force` 类操作必须人确认;
3. **可复现性**: AI 决策入 `ai_log.jsonl` (prompt 版本+模型+输入哈希);
   发表流程必须能用纯规则参数重放;
4. **执行权在人**: copilot 与 MCP 均无执行工具; 命令草稿须过白名单并复制到终端;
5. **AI 写入只碰 ai_* 列**: 沿袭并收紧 gsa_sra.info.py 的"只回填空白"纪律;
6. **成本与隐私**: LLM 只见 digest/证据表/卡片, 永不见 fastq 与全量 TSV;
   key 走环境变量, copilot 见内联 key 即警告。

## 5. 威胁模型与对策

| 威胁 | 场景 | 对策 |
|---|---|---|
| 幻觉 | AI 编造"检测到新病毒" | 提示词明令数字只出自 digest; 解读与数据分离 (ai_* 列/独立 section) |
| 提示词注入 | 恶意构造的 contig id / 元数据进入 prompt 后指令劫持 | AI 输出只进封闭动作集与白名单校验器, 无既执行通道; validate 三道防线在 AI 下游 |
| 提权执行 | AI 客户端要求 MCP 跑管线 | 协议层不存在执行工具; -32602 |
| 密钥泄露 | key 进 shell 历史/报告 | 环境变量优先; copilot 检测内联 key 并警告; 日志只记 masked key |
| 离线断裂 | 无网/无额度 | 全部功能有规则回退: advisor 规则完整可用, triage 有 rules 模式, chat 有关键词匹配 |
| 成本失控 | 大 TSV 直接进 prompt | digest 瘦身上限 30k 字符; triage 分批 (默认 20 行/批); 无流式失控面 |

## 6. 与现有资产的对账 (复用而非重建)

| 现有机制 | 在本架构中的角色 |
|---|---|
| `--stage` + 检查点 + 退出码 0/1/2 | copilot 校验器的语义基础; 未来 agent 循环的重试依据 |
| `io_layout.py` 双布局 | run_digest 的目录解析底座 (standard/legacy 通吃) |
| `report_style.py` 调色板 | ai_interpretation_section 的样式来源 (视觉零漂移) |
| gsa_sra.info.py AI 清洗 | 提示词纪律与 env 口径 (`DEEPSEEK_API_KEY`) 的出处 |
| eve_distinguish s4_filter | triage 动作模式的模板 |
| `_record_provenance` | ai_log.jsonl 的对接点 (后续可并入统一 provenance) |

## 7. 快速上手

```bash
# AI 能力配置 (可选; 不配则全部功能走规则回退)
export MMPV_AI_API_KEY=sk-xxx            # 或 DEEPSEEK_API_KEY (gsa_sra 同口径)
# export MMPV_AI_BASE_URL=http://localhost:11434/v1   # Ollama 例
# export MMPV_AI_MODEL=deepseek-chat

python mmpv_ai.py digest  <项目根> --write     # L1: 状态落盘 90_Handoff/state/
python mmpv_ai.py cards   discovery rescue     # L1: 阶段卡片
python mmpv_ai.py advise  <项目根> --ai        # L3: 建议 + AI 叙事
python mmpv_ai.py triage  candidates.tsv --evidence genomad,virsorter2 --mode rules
python mmpv_ai.py validate "python virome_discovery_pipeline/virome_pipeline.py ..."
python mmpv_ai.py chat                         # L4: 自然语言副驾驶 (不执行)
python mmpv_ai.py agent <项目根> --dry-run     # L4: 受控 agent — 先看计划
python mmpv_ai.py agent <项目根> --yes \
    --set "项目名=demo" --set "原始FASTQ目录=D:/raw"
python mmpv_ai.py report <项目根>              # 项目级主报告: 结果×过程×AI 参与
python mmpv_ai.py enrich <项目根>              # 事后补 AI: 已生成报告补插解读+面板 (幂等)
python mmpv_ai.py triage --preset ident --project <项目根>   # 鉴定证据矩阵实战标注
python mmpv_ai.py provenance <项目根>          # AI 参与溯源 → AI_PROVENANCE.md (方法章)
python mmpv_ai.py mcp                          # L4: MCP server
```

MCP 接入 (ZCode 例, `~/.zcode/workspace/<ws>/.zcode/config.json`):

```json
"mcp": { "servers": { "mmpv": { "type": "stdio", "command": "python",
         "args": ["<仓库根>/mmpv_common/mcp_server.py"] } } }
```

测试: `python -m pytest mmpv_common/tests/ -q` (149 passed, 全离线, LLM/executor/confirm 全桩注入)。
配 key 后自验真实端到端: `python mmpv_ai.py advise <项目根> --ai` 看叙事是否只引用 digest 数字。

## 8. 交付清单 (2026-10-08, 落地顺序四阶段全完成)

| 文件 | 层 | 内容 |
|---|---|---|
| `mmpv_common/ai.py` | L2 | provider/提示词注册表/溯源日志/离线降级/6 组内置提示词/`append_ai_section_to_html` 统一接线点 |
| `mmpv_common/stage_cards.py` | L1 | 8 管线 ~60 stage 卡片 + self_check |
| `mmpv_common/run_digest.py` | L1 | digest 构建/落盘/瘦身/guess_state_dir + CLI |
| `mmpv_common/advisor.py` | L3 | 9 规则顾问 + AI 叙事 |
| `mmpv_common/triage.py` | L3 | 只加 ai_* 列的通用标注器 (rules/ai 双模式) |
| `mmpv_common/copilot.py` | L4 | 8 编排器命令白名单 + 三道防线 + dry-run + AI 意图映射 |
| `mmpv_common/stage_tools.py` | L4 | stages 工具化: JSON Schema/前置/后置/退出码语义 + plan_chain |
| `mmpv_common/agent.py` | L4 | 受控 agent 循环: 确认门/预算/退出码语义/agent_runs 留痕 |
| `mmpv_common/mcp_server.py` | L4 | 零依赖 stdio MCP (7 工具, 无执行) |
| `mmpv_common/process_report.py` | L2 | 过程报告层: 时间线/AI 参与聚合 → 过程面板 (嵌五报告) + 项目主报告 (含互链/演化图/用量/triage 分区); v0.2.1 起标记渲染走 Jinja2 模板 |
| `mmpv_common/report_template.py` | L2 | 报告层 Jinja2 统一出口: autoescape + StrictUndefined + `templates/` (base 骨架/macros 组件宏/过程面板/主报告/⑦ 提交报告); `report_style.py` 组件函数改为宏薄包装 (单一实现源) |
| `mmpv_common/provenance_export.py` | L3 | AI 参与溯源一键导出 (JSON + 方法学 MD, 中英文段落) |
| `mmpv_ai.py` | 入口 | digest/cards/advise/triage(+preset)/ask/chat/validate/agent/report/**provenance**/mcp |
| 五报告接线 | L2 | preprocess/eve/host_expression/submission/discovery 报告 main 里 try/except 接 `append_ai_section_to_html` (纯增量, 离线零改动; **eve/host_expression 报告随 ⑥⑧ 归档 2026-10-11, 接线代码保留**) |
| `mmpv_common/tests/` | 守门 | 149 tests: 溯源/降级/digest 不变量/宪法/白名单/MCP/agent 语义/接线离线不破坏/模板层 (存在性·StrictUndefined·autoescape·离线无外链·回归契约子串) |

## 9. 后续路线

- ~~P1 报告接线~~ ✅ 已完成 (五报告统一接 AI, 离线零影响; 真实文案待配 key 自验);
- ~~P1.5 过程报告~~ ✅ 已完成 (过程面板嵌入五报告 + 项目主报告 `mmpv_ai.py report`);
- ~~P2 triage 实战化~~ ✅ 已完成 (`--preset ident/eve_distinguish` + REVIEW 计数入面板/主报告);
- ~~P3 provenance 合流~~ ✅ 已完成 (`mmpv_ai.py provenance` → JSON + 方法学 MD;
  编排器 `_record_provenance` 侧的物理合并刻意不做 —— 平台留痕已单点汇聚于
  90_Handoff/state, 动编排器反而破坏"管线零侵入"原则);
- ~~P5 报告间导航~~ ✅ 已完成 (主报告↔管线报告互链 + digest 演化 SVG + token 用量面板);
- ~~P4 agent 循环~~ ✅ 已完成 (受控版: 计划规则化/确认门/预算/留痕; 全自动
  agent 需求驱动时再评估放宽, 红线 §4 不动)。

---

*MMPV = Massive Meta-mining of Plant Viruses · AI 层版本 mmpv_common 0.2.1*
