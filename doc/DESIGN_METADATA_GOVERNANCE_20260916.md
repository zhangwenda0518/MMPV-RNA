# Metadata 入口治理模块 — 设计方案与实施报告

**日期**：2026-09-16
**触发**：老师问「BioAider 里的时间解析脚本能不能吸收进来；时间和地理的检查和矫正，补充到我们的流程中对 metadata 进行格式化，以免下游错误。怎么设计合适？」
**状态**：**已实施并全量验证通过**（老师 2026-09-16 拍板三项决策后执行）

---

## 老师拍板结论（2026-09-16）

| 问题 | 老师决定 | 实施结果 |
|---|---|---|
| 模块范围 | **体检 + 格式化另存新文件** | 产出 `metadata_std.csv` / `metadata_report.tsv` / `metadata_summary.txt`，原件不动 ✅ |
| 地理 `Unknown` 层级 | **删除** | `parse_location(drop_unknown=True)` 删层并留痕于 `dropped_levels` ✅ |
| `virphy_bridge` 月初口径 | **接入统一入口 `mode='start'`** | `virphy_bridge._to_decimal` 改为委托，**30/30 逐位一致** ✅ |

### 实施产物

| 项 | 内容 |
|---|---|
| 新模块 | `virome_phylo_pipeline/utils/metadata_governance.py`（**711 行**，含独立 CLI 13 个参数） |
| 新 stage | `metadata`（`STAGE_ORDER` 第 2 位，`prep` 之后 / `online` 之前） |
| 改动文件 | `phylo_pipeline.py`（+stage/参数）、`virphy_bridge.py`（-内联，+委托） |
| 新参数 | 7 个 `--metadata_*`（第二轮追加 2 个：`--metadata_no_blank_placeholders` / `--metadata_drop_placeholder_rows`） |
| 验证脚本 | `verify_metadata_governance.py`(26) / `verify_metadata_stage.py`(17) / `verify_virphy_bridge_governed.py`(8) / `verify_metadata_dataquality.py`(20) |
| 回归 | **全绿，零 FAIL**（原 31 PASS 基线 + 新增 71 PASS） |

---

## 0. 先说结论（三句话）

1. **BioAider 的 `date_to_float` 不能直接吸收**——实测它有 3 处比我方更弱（只给年不转换 / 整百年闰年判错 / 静默回显脏值不报错）。**但它的「分隔符归一化」思路值得吸收**，正是我方 `decimal_year` 目前缺的一环。**已吸收这一点。**
2. 真正值得建的不是「再写一个日期解析器」，而是一层 **入口 metadata 治理 + 体检报告**：解析、归一化、标注精度、加地理校验，**并在越界/歧义时拦住或显式标注**，而不是静默放行。
3. 关键价值在于 **「体检报告」**：老师当前两批真实 metadata 里已经查出 **450 条坐标全是占位符 `XX.XX N XXX.XX E`、42 条 `collection_date` 是字面量 `YYYY-MM-DD`、36 条 `DD-Mon-YYYY` 被静默丢弃、Unknown 地理层级会被静默降级成国家级**。这些问题在实施前全都不报错，现在全部进报告。

---

## 1. BioAider `date_to_float` 实体取证

### 1.1 源码位置（两个副本，实现体逐行一致）

| 位置 | 行号 |
|---|---|
| `_bioaider_re/pyz_src/custom_libraries/public_functions.py` | `:8709` |
| `docs/bioaider_ref_src/public_functions.py` | `:497` |
| GUI 调用方 `_bioaider_re/pyz_src/convert_tools/DateToFolat_gui.py` | `from custom_libraries.public_functions import date_to_float` |

反编译自 `pyz_dec/custom_libraries/public_functions.pyc`。

### 1.2 它的实现要点

```python
date_list = re.split("[./_-]", date_str)     # ← 一个正则吃掉 4 种分隔符
if len(date_list) == 3:   # 年 月 日
    days = cum_month_days[0..month-2] + day
    results = str(float(year) + days / 366.00000001  or  365.00000001)
elif len(date_list) == 2:  # 年 月 → 补 1 日（月初）
    ...
else:
    results = date_str        # ← 原样回显（含只给年 "2013"）
except:
    results = date_str        # ← 静默吞掉一切
```

### 1.3 与我们 `utils/decimal_year.py` 的逐用例对拍（24 用例，**实测**）

脚本：`_consistency_check_20260916/verify_bioaider_date.py`

| 输入 | BioAider | 我方 | 差(天) | 判定 |
|---|---|---|---|---|
| `2013` | `'2013'`（原样回显） | `2013.455556` | 166.28 | BioAider 不转换 |
| `2014.05` | `2014.3315068` | `None`(拒) | — | **我方拒，BioAider 强解** |
| `2014_05` | `2014.3315068` | `None`(拒) | — | **我方拒，BioAider 强解** |
| `2022_11-11` | `2022.8630137` | `None`(拒) | — | **我方拒，BioAider 强解** |
| `2015/10/01` | `2015.7506849` | `2015.750000` | 0.25 | 一致（差 6 小时） |
| `2014-05` | `2014.3315068` | `2014.370968` | 14.40 | **口径不同** |
| `20-Jun-2019` | `'20-Jun-2019'`（回显） | `None`(拒) | — | 两者都拒 |
| `2019-13-01`（月越界） | `2020.0027397` ← **强解成 2020 年** | `None`(拒) | — | **BioAider 危险** |
| `12019` | `'12019'`（回显） | `None`(拒) | — | 两者都拒 |
| `NA` / `not a date` | 原样回显 | `None` | — | 两者都拒 |

**闰年判据**（实测）：`1900` → 真闰年 `False`，BioAider 判 `True`；`2100` 同。BioAider 用 `year % 4 == 0`，缺「整百年须 %400」规则。

### 1.4 结论：**不可直接吸收，可选择性吸收一条**

| 维度 | BioAider | 我方 | 谁更好 |
|---|---|---|---|
| 分隔符 `.` `/` `_` `-` 混用 | ✅ `re.split("[./_-]")` | ❌ 只 `replace('/','-')` | **BioAider** ← 唯一值得吸收 |
| 只给年 `YYYY` | 回显字符串（下游易崩） | 取 6-15（年中） | 我方（且老师已定「保留年初」，见 §4.3） |
| 闰年 | `%4`（1900/2100 错） | `calendar.monthrange` | 我方 |
| 脏值处理 | 静默回显 | `None` / fail-loudly | 我方 |
| 越界日期 | 静默强解（`2019-13-01`→2020） | 拒 | 我方 |
| 月→日补全 | 补 **1 日**（月初） | 补 **15 日**（月中） | 取决于口径，见 §4.3 |

> **一句话**：BioAider 这个脚本是「能用就行」的 GUI 工具逻辑，**不适合直接进科研管线**。它唯一优于我们的是分隔符容忍度——我们把这个思路拿过来，实现照自己的标准写。

---

## 2. 当前 metadata 的真实病灶（全量扫描，实测）

脚本：`_consistency_check_20260916/diagnose_metadata_20260916.py`

### 2.1 会议纪要

| 文件 | 行数 | 问题 |
|---|---|---|
| `virome_submission_pipeline/submission/unified_metadata.csv` | 450 | **坐标 450/450 全占位符** `XX.XX N XXX.XX E`；**42 条日期是字面量 `YYYY-MM-DD`**；42 条地理是字面量 `Country:Region` |
| `submission_gui/samples/unified_metadata_barbarum_real.csv` | 75 | **39 条 `YYYY`**；**36 条 `DD-Mon-YYYY` 被 `_to_decimal` 静默丢弃**；1 条 `Unknown:Unknown` |
| `submission_gui/samples/_test_roundtrip_out.csv` | 75 | 同上（roundtrip 复刻） |
| `submission_gui/samples/sample_selfseq.csv` | 4 | 4 条日期无法解析 |
| `public_metadata_pipeline/.../Global_Unified_Metadata_Core14.csv` | 148 | 日期 100% 规范，但 **`Location` 列 148/148 全空** |
| `biosoft/VirPhyKit/Example/TreeDater-LTT/*.csv` | 476 | 日期列是**已转好的小数年** `2007.569473`；地理 `VietNam`/`HongKong`/`New_Zealand` |
| `biosoft/VirPhyKit/Example/TreeTime-RTT/H3N2/*.csv` | 476 | 同上 |

### 2.2 四类「下游会错但不报错」的具体情形

**A. 坐标占位符被当成真坐标**
`XX.XX N XXX.XX E` 出现在 450 行。若下游有某条路径直接 `float()` 或正则抽数字，会得到垃圾；若走 `geo_resolver` 则被忽略，两种行为都不一致。

**B. 日期占位符被当成真日期**
`collection_date == "YYYY-MM-DD"` 42 条。它与真实 `2024-04-24` **视觉上同类**，任何「检查格式是否合法」的校验都会放行。

**C. 地理层级含 `Unknown` 却被静默降级**
```
China, Unknown, Yinchuan_AI      × 9
China, Unknown, Unknown_AI       × 2
```
`geo_resolver._resolve_offline` 在 `:193`/`:201` 会 `continue` 跳过 `Unknown` 段，于是这两串会被解析成 **China 国家质心 (35.00, 103.00)**。下游以为拿到了省级/市级坐标，实际拿到的是全国中心——**误差上千公里，且无任何告警**。

**D. 两种地理编码并存**
- `China:Ningxia`（冒号，GenBank `geo_loc_name` 标准）
- `China, Ningxia, Yinchuan_AI`（逗号+`_AI` 后缀）

`geo_resolver._norm` 只按 `[,;]` 拆段，**冒号不拆**。`'china:ningxia'` 整串进查表 → 不命中 → 走 L3 段内反查 → 侥幸命中 `ningxia`。**碰巧对，但脆弱**。

**E. 日期列可能已经是小数年**
VirPhyKit 的两份 Example 里 `date` 列是 `2007.569473`。如果同一列在不同数据集里有时是字符串、有时是小数年，**下游无法区分**，重复转换会得到荒谬结果（`float(2007.569473)` 再当年算 → 2007.57 年）。**目前没有任何地方检查这一点。**

---

## 3. 设计：`utils/metadata_governance.py`（新模块，不动老代码）

### 3.1 定位

放在**入口**（`prep` stage 之前 / `seq_clean` 之前），职责单一：**读原始 metadata → 输出「标准化 metadata + 体检报告」**。下游一律消费标准化产物。

```
原始 metadata.csv
        │
        ▼
┌───────────────────────────────┐
│ utils/metadata_governance.py  │
│  ① 列名归一（别名表）          │
│  ② 日期：解析→精度分档→小数年  │  ← 唯一权威
│  ③ 地理：占位符→层级拆分→校验  │
│  ④ 校验：越界/占位符/冲突      │
│  ⑤ 产出标准化 CSV + 体检报告   │
└───────────────────────────────┘
        │                    │
        ▼                    ▼
  metadata_std.csv      metadata_report.tsv
  （下游唯一输入）        （人看，含逐行问题）
```

### 3.2 日期侧怎么设计（核心）

**统一入口函数**，返回**结构化结果**而不是裸 float：

```python
@dataclass
class DateInfo:
    raw: str              # 原始字符串
    iso: Optional[str]    # 归一化后的 ISO（YYYY-MM-DD / YYYY-MM / YYYY）
    decimal: Optional[float]   # 小数年
    precision: str        # 'day' | 'month' | 'year' | 'none'
    confidence: float     # 1.0 精确 / 低值表示需人工确认
    issue: Optional[str]  # 越界/歧义/占位符 的说明
    source_format: str    # 'ISO' | 'DD-Mon-YYYY' | 'DD/MM/YYYY' | 'decimal_year' | ...
```

要点：

1. **分隔符归一化**（吸收 BioAider 的唯一优点）
   先把 `.` `_` `/` `-` 四种分隔符统一成一种再解析。这样 `2014_05`、`2014.05`、`2022_11-11` 都能被正确识别，而不是像现在直接失败。

2. **格式登记表**（显式白名单，不靠猜）

   | 格式 | 例 | 归一化 | precision |
   |---|---|---|---|
   | ISO 全 | `2024-04-24` | 原样 | day |
   | ISO 到月 | `2019-06` | 补 `-15` | month |
   | 仅年 | `2024` | 补 `-06-15` | year |
   | GenBank | `20-Jun-2019` | → `2019-06-20` | day |
   | 欧式 | `15/07/2016`、`15.07.2016` | → `2016-07-15` | day |
   | 已是小数年 | `2007.569473` | 反解为 `2007-07-27` 附近 | day(推断) |

3. **小数年的两套口径显式标注**（回应当前 `virphy_bridge:726` 与 `decimal_year` 的差异）
   不再让「年初 vs 年中」藏在实现里。`DateInfo.iso` 是唯一真相，小数年由**调用方显式指定口径**产出：
   - `mode='mid'` → 月中/年中（`decimal_year` 现口径）
   - `mode='start'` → 月初/年初（`virphy_bridge` 现口径，老师已定保留）
   两条路径各自的产物在报告里写明口径，**将来若要统一，只需改一个开关**。

4. **占位符拦截**
   字面量 `YYYY-MM-DD`、`DD-Mon-YYYY`（作为整列默认值出现时）、`Unknown`、`NA` 一律判为「占位符」→ `precision='none'` + `issue='placeholder'`，**不再进入下游**。

### 3.3 地理侧怎么设计

在 `geo_resolver` **之上**包一层校验，不改它：

1. **占位符识别**：`Country:Region`、`XX.XX N XXX.XX E`、`Unknown`、空
2. **层级一致性检查**（新增，最重要）：
   ```
   China, Unknown, Yinchuan_AI
   → 记录 components = ['China','Unknown','Yinchuan_AI']
   → 标 issue='unknown_level'（第 2 级未知）
   → 解析结果标 fallback_level='country'，精度=国家
   → 下游若要求省级，可据此拒绝
   ```
3. **编码归一**：`China:Ningxia` → 统一成 `China, Ningxia`；`_AI` 后缀剥离（但记录原值）
4. **冲突检查**：坐标列若真有值，与地名解析结果比对，距离超阈值 → 标 `issue='geo_conflict'`

### 3.4 关键设计原则（三条，请老师重点看）

| 原则 | 说明 |
|---|---|
| **不猜** | 解析不出就标 `none` + 记 issue，绝不套默认值（这是 BioAider 最大的问题） |
| **不静默** | 任何降精度（省级→国家级）、任何丢弃（无日期/无地理）都进报告，可统计条数 |
| **不改语义** | 模块只做「归一 + 标注」，**不改变**任何现有判定阈值；`virphy_bridge` 的年初口径保持不动 |

---

## 4. 三个问题的决策与实施

### 4.1 模块范围 → **A：体检 + 格式化，另存新文件** ✅

原件完全不动。产出：

| 文件 | 内容 |
|---|---|
| `data/metadata/metadata_std.csv` | 标准化 metadata（下游唯一输入）；新增 `_date_decimal` / `_date_precision` / `_date_issue` / `_geo_unresolved` 四列 |
| `data/metadata/metadata_report.tsv` | 逐行体检（17 列：日期原值/ISO/小数年/精度/格式/置信度/问题，地理原值/归一/层级/删除/未解析/问题） |
| `data/metadata/metadata_summary.txt` | 人类可读汇总 |

编排层把 `args.metadata` 改写为 `metadata_std.csv`，原值存 `prep['_metadata_original']`。
**失败降级**：治理失败只告警，继续用原 metadata（不阻断流程）。
**幂等**：输入已是 `metadata_std.csv` 时直接复用，不二次处理。

### 4.2 地理 `Unknown` 层级 → **删除** ✅

```python
parse_location("China, Unknown, Yinchuan_AI")
→ components=['China', 'Yinchuan_AI']      # 保留
  dropped_levels=['Unknown']               # 删除并留痕
  issue='unknown_level', unresolved=True   # 显式标注
```

`_is_unknown_level()` 还覆盖 `Unknown_AI` 这类「占位词 + 来源标记」写法（前缀匹配 `unknown_` / `unknown-` / `unknown `）。

> **与 `geo_resolver` 的区别**：`geo_resolver._resolve_offline` 会 `continue` 跳过 `Unknown` 段后**仍返回国家质心并标记为 `'province'`**（下游无从分辨）；本模块**删除该层 + 显式标 `unknown_level`**，下游必须看到才能决定是否使用。

### 4.3 `virphy_bridge` 接入 → **A：`mode='start'`，数值逐位不变** ✅

`virphy_bridge.run_treedater_ltt` 里的内联 `_to_decimal` 已改为：

```python
from utils.metadata_governance import decimal_year as _gov_decimal_year

def _to_decimal(d):
    try:
        return _gov_decimal_year(d, mode='start')
    except (ValueError, IndexError, TypeError):
        return None
```

**验证结果（30/30 逐位一致，真差异 0）**：

| 分组 | 结果 |
|---|---|
| 只给年（`2013`/`2019`/`2020`/`1999`/`2100`） | 一致 |
| 只给年月（`2014-05`/`2019-06`/`2014/05`/`2021-12`） | 一致 |
| 完整日期（9 个） | 一致 |
| 脏值（`''`/`NA`/`nan`/`not a date`） | 一致（都 None） |
| 越界（`2019-13-01`/`2019-99-99`） | 一致（都 None） |
| 占位符（`YYYY-MM-DD`/`YYYY`/`XX.XX N XXX.XX E`） | 一致（都 None） |

**3 处「有意收紧」**（新版拒绝、旧版返回荒谬值）：
`0000`(旧→0.0) / `9999`(旧→9999.0) / `2019-02-30`(旧→2019.1637)
→ **已全量扫描 883 个项目 CSV，零条真实数据落在此集**，无实际风险。

`verify_fix_decimal_year.py` 的 F5 组独立佐证：`2020` / `2020-06` / `2021` / `1999-11` 与旧实现**零差异**。

---

## 5. 验证计划（已按此执行完毕）

1. **对拍**：新模块的日期解析，对 24 个用例与 `decimal_year` / `virphy_bridge` 三路对拍，逐位比对（`verify_bioaider_date.py` 已备）✅
2. **零漂移**：`mode='mid'` / `mode='start'` 两条路径，对现有 metadata 全量算小数年，**与旧实现逐位一致**（浮点容差 0）✅
3. **端到端**：拿 `unified_metadata_barbarum_real.csv` 跑一遍，报告应精确报出 39 `YYYY` / 36 `DD-Mon-YYYY` / 1 `Unknown:Unknown`✅
4. **回归**：现有 31 PASS 全绿，不新增 FAIL ✅（叠加本轮新增 71 PASS，全绿）

---

## 6. 本文档的证据清单

| 脚本 | 用途 | 结果 |
|---|---|---|
| `verify_bioaider_date.py` | BioAider vs 我方 24 用例对拍 | 一致 0 / 不一致 12 / BioAider 回显 0（详见 §1.3） |
| `diagnose_metadata_20260916.py` | 全项目 metadata 质量扫描 | 见 §2.1 |
| `verify_metadata_governance.py` | 治理模块零漂移 | **26 PASS / 0 FAIL** |
| `verify_metadata_stage.py` | stage 端到端 + 幂等 + 降级 | **17 PASS / 0 FAIL** |
| `verify_virphy_bridge_governed.py` | virphy_bridge 接入前后逐位比对 | **8 PASS / 0 FAIL**（30/30 一致） |
| `verify_metadata_dataquality.py` | 数据质量处置（清空/剔除/坐标不触发/降级） | **20 PASS / 0 FAIL** |

---

## 7. 实施中发现并修复的两个真 bug

### 7.1 占位符正则**误伤真实地理**（`China:Ningxia` 72 条全被打成 placeholder）

**现象**：首版写了一条泛化正则 `^[a-z]+\s*[:：]\s*[a-z]+$`，想兜住「任意 `Word:Word` 占位符」。
**后果**：把**真实地理** `China:Ningxia` / `China:Beijing` / `China:Gansu` 全部判为 placeholder，`unified_metadata_barbarum_real.csv` 里 72 条地理列被清空。
**修复**：删掉泛化正则，只精确匹配已知占位字面量 `^country\s*[:：]\s*region$`。
**教训**：**这是本会话第 4 次踩「正则/字样误伤真实数据」**。判据必须先用真实数据回归，不能只靠构造用例。

> 若未发现：`metadata_std.csv` 的地理列会被大面积清空，且静默——正是本模块要防的那类错误。**说明判据必须跑真实数据，不能只跑构造用例。**

### 7.2 stage 会改写 `args.metadata`，重复调用会读到自己的产物

**现象**：`run_stage_metadata` 成功后把 `args.metadata` 指向 `metadata_std.csv`（为让下游消费产物）。
**后果**：同一 `args` 被重复调用时，第二次读到的输入已是产物 → 二次治理 / 永久改名。
**修复**：加**幂等检查**——输入 basename 为 `metadata_std.csv` 且父目录为 `metadata` 时直接复用，不重复处理。

---

## 8. 明确未做的事（防止误读）

| 项 | 状态 |
|---|---|
| 改动 `utils/geo_resolver.py` | **未动**（治理层在其之上包装，不改它） |
| 改动 `utils/decimal_year.py` | **未动**（`mode='mid'` 复刻其公式，逐位一致） |
| 改动 `utils/import_export.py` | **未动**（其 `_to_decimal_year` 早已委托 `decimal_year`） |
| 改动 `virphy_bridge` 的**数值口径** | **未动**（`mode='start'`，30/30 逐位一致） |
| 修改任何原始 metadata 文件 | **未动**（一律另存新文件） |
| 改任何判定阈值 / 门控方向 | **未动** |
| 修复 `DD-Mon-YYYY` 36 条**被静默丢弃**的下游影响 | **已根除**——新解析器能解析，`virphy_bridge` 接入后不再丢弃 |
| 修复 `XX.XX N XXX.XX E` / 42 条字面量日期 | **仅标注，未改数据**（属上游数据质量问题，需老师决定是否回补） |

---

## 9. 追加实施：数据质量问题的处置口径（2026-09-16 第二轮，老师拍板）

### 9.1 老师的三条口径澄清

老师先说「**数据质量问题 剔除**」，我随即用 AskUserQuestion 澄清三个必须定死的细节，老师答复：

| 追问 | 老师答复 | 实施方式 |
|---|---|---|
| 剔除**粒度**？ | **仅清空该字段** | 默认只把该单元格置空，**行保留** |
| 剔除**后果**？ | **照常剔除，报告警告**（全被剔光则降级用原文件） | `--metadata_drop_placeholder_rows` 可整行剔除；`kept=0` 时降级 |
| 坐标列怎么算「无地理」？ | **坐标列为空 → 判该行无地理**（由**地名列**决定） | 坐标占位符**只清空、不触发剔除** |

**为什么必须问这一条**：如果按"整行剔除"理解「剔除」，而坐标列 `XX.XX N XXX.XX E` 在 `unified_metadata.csv` 里是 **450/450 全占位符**，整表会被剔光——与老师"照常剔除"的意图南辕北辙。问清后实现才正确。

### 9.2 新增开关（3 个）

| 参数 | 默认 | 行为 |
|---|---|---|
| `--metadata_no_blank_placeholders` | False（即**默认清空**） | 不清空占位符字段，保留原值 |
| `--metadata_drop_placeholder_rows` | False | 含占位符字段的行**整行剔除** |
| `--metadata_drop_no_location` / `--metadata_drop_no_date` | False | 无地理 / 无日期的行剔除（前一轮已有） |

**关键逻辑（坐标占位符不触发剔除）**：

```python
blanked_for_drop = [b for b in blanked if b not in co_cs]   # co_cs = 坐标列集合
if drop_placeholder_rows and blanked_for_drop:
    drop_why = "placeholder"
```

即：**只有非坐标列的占位符**才构成整行剔除的理由。

### 9.3 实测口径（真实数据，可复现）

| 文件 | 行数 | 默认（仅清空） | `--drop-placeholder-rows` | `--drop-no-location` |
|---|---|---|---|---|
| `unified_metadata.csv` | 450 | **保留 450**；清空 Lat_Lon 450 / date 42 / geo 42 | **剔 42 → 保留 408** | **剔 42 → 保留 408** |
| `unified_metadata_barbarum_real.csv` | 75 | **保留 75**；清空坐标 9 | 保留 75（坐标占位符不触发） | 保留 75 |

> 坐标列 450/450 全占位符，若触发剔除会 `kept=0`——**修正前正是这个行为，是 bug**。

### 9.4 新增验证脚本

| 脚本 | 内容 | 结果 |
|---|---|---|
| `verify_metadata_dataquality.py` | Q1–Q8：清空 / 剔除 / 坐标不触发 / 报告结构 / 降级 | **20 PASS / 0 FAIL** |

同时 `STAGE_REFERENCE.md` 新增 **§2.3 占位符处置** 与 §2.4 实测表，把上述口径写进流程文档。

### 9.5 本节实施中修正的一个真 bug

**`drop_placeholder_rows` 把坐标占位符也算作触发条件** → `unified_metadata.csv` 保留 0 行（坐标 450/450 全占位符）。
这与老师「坐标列为空 → 判该行无地理」的答复**直接矛盾**。
**修复**：`blanked_for_drop` 排除坐标列（见 §9.2 代码）。修正后剔 42 保留 408，与预期一致。

> 本条与 §7.1 同源：**判据必须在真实数据上回归**。构造用例（坐标只有几条占位符）永远发现不了 450/450 全占位符的畸形态。

---

## 10. 追加实施：坐标「先补算，再清空」（2026-09-16 第三轮，老师拍板）

### 10.1 触发：老师质疑「之前可以解析，怎么到你就不行了」

老师原话：
> 「src-Lat_Lon 的 `XX.XX N XXX.XX E` 为什么不对呢，之前不是可以解析吗，怎么到你，就不行了呢。
>  collection_date 之前也可以的呀，而且我还让你吸收 BioAider 的处理方法，怎么不行呀。
>  src-geo_loc_name 之前也正常的呀」

**取证结论：老师是对的，本报告 §9 的表述有误。**

### 10.2 我上一版写错了什么

§9 把 450 条坐标占位符说成"数据缺失"倾向。**更正后的表述**：坐标**没被填上，但 408 条可以算出来**。

| 证据 | 内容 |
|---|---|
| `geo_to_latlon` 实测可工作 | `unified_metadata.py:216`；`China, Beijing, Beijing_AI` → `39.90 N 116.40 E` |
| 另一文件就是这么填的 | `unified_metadata_barbarum_real.csv` 真坐标 `37.48 N 105.68 E`(30) / `39.90 N 116.40 E`(19) / `38.47 N 106.27 E`(15)，**与 `geo_to_latlon` 输出逐字一致** |
| 补算覆盖率 | **408 / 450 = 90.7%**（408 行有真实地名） |

### 10.3 但原有解析能力**未被破坏**（附对拍证据）

- **时间线**：`unified_metadata.csv` 生成于 **2026-06-23**，本模块 **2026-09-16** —— 数据早 3 个月
- **新旧实现对拍（tol=0）**：

| 值 | 性质 | 旧实现 | 新实现 |
|---|---|---|---|
| `2024-04-24` / `2019-06` / `2024` / `18-Jan-2023` | **真数据** | ✓ | ✓ **完全一致** |
| `XX.XX N XXX.XX E` / `YYYY-MM-DD` / `Country:Region` | **占位符** | None | None |

→ 真数据一字未变；占位符新旧都返回 None（**解析能力未变**）。

### 10.4 真正的缺陷：只清空、没补算

**本报告 §9 的实现把 450 条坐标全部清空 → 覆盖率 0%，把 408 行本可挽救的坐标清掉了。**
老师拍板改为「**先补算，再清空**」。

### 10.5 实施

| 项 | 内容 |
|---|---|
| 新函数 | `metadata_governance.coord_from_location(loc) -> (coord, source)` |
| 逻辑来源 | 移植自 `virome_submission_pipeline/unified_metadata.py:216 geo_to_latlon`（30 城市字典 + 中国中心点兜底），**逐字一致**；**不 import** 该模块（避免跨 pipeline 耦合） |
| 唯一增强 | 额外返回 `source`：`city` / `country_centroid`，供下游分辨误差量级 |
| 新参数 | `govern_table(..., derive_coords=True)`；CLI `--no-derive-coords`；编排层 `--metadata_no_derive_coords` |
| stats 新增 | `coords_derived` / `coords_derived_by_source` / `coords_unrecoverable` |
| 报告新增列 | `coord_source`（如 `src-Lat_Lon:city`），报告共 **19 列** |
| 处理顺序 | 坐标占位符/空 → 先 `coord_from_location(地名)` → 成功填入并标来源；失败才清空 |

**关键：坐标补算值不触发整行剔除**（沿用 §9.2 约束），且**不覆盖已有真坐标**。

### 10.6 实测

| 文件 | 坐标占位符 | 补算成功 | 推不出 | 覆盖率 |
|---|---|---|---|---|
| `unified_metadata.csv` | 450 | **408**（city 406 / country_centroid 2） | 42 | **0% → 90.7%** |
| `unified_metadata_barbarum_real.csv` | 9 | **5**（均 `China:Ningxia`） | 4 | — |

推不出的 42 行：地名恰好也是 `Country:Region`（**与 42 条日期占位符是完全同一批行**）。
barbarum 推不出的 4 行：`United States:Maryland`×2 / `United States:Massachusetts`×1 / `Unknown:Unknown`×1。

### 10.7 能力边界（如实记录，非缺陷）

`geo_to_latlon` / `coord_from_location` **只支持中国地名**（30 城市 + 中国中心点兜底）。
非中国地名（如 `United States:Maryland`）**推不出**，会被清空。这是原实现的固有局限，未扩展。

### 10.8 验证

| 脚本 | 内容 | 结果 |
|---|---|---|
| `verify_metadata_coord_derive.py`（新） | C1–C12：补算/来源/清空/幂等/真值不覆盖/**与 `geo_to_latlon` 逐字对拍**/格式/能力边界 | **18 PASS / 0 FAIL** |
| `verify_metadata_dataquality.py`（更新判据） | Q2c/Q2d/Q2f/Q7b/Q8b 由「全清空」改为「先补算」；新增 Q2x 回归"不补算"旧行为、Q8a2/Q8b2 测新列 | **26 PASS / 0 FAIL**（原 20） |

**全量回归零 FAIL**：governance 26 / stage 17 / virphy_bridge_governed 8 /
dataquality 26 / coord_derive 18 / fix_decimal_year 32 / task12 17 /
fix_complete_deletion 14 / D11 8 / seq_clean 10 / clean_stage 9 / date_precision 5。

### 10.9 判据过时的说明（不是新 bug）

`verify_metadata_dataquality.py` 首跑 **5 FAIL**（Q2c/Q2d/Q2f/Q7b/Q8b）——
**全部是判据假设了"全清空"旧行为**，数值本身都对（42 / 4）。
已按新行为更新，并**另加 Q2x 段落**用 `derive_coords=False` 继续守住"清空"功能，避免旧能力失去覆盖。

---

## 11. 追加实施：坐标补算扩展为**全球地理**（2026-09-16 第四轮，老师拍板）

### 11.1 触发：老师指出能力边界不可接受

老师原话：
> 「修复，不能只是中国的省市呀，全世界的地理怎么做，
>  看看 D:\桌面\植物病毒分析平台\git-repo 的项目，有可以借鉴参考的处理方法嘛」

§10 交付时明确报告过一个边界：`geo_to_latlon` / 首版 `coord_from_location`
**只支持中国地名**（30 城市 + 中国中心点兜底），非中国地名推不出。
老师判断这不可接受 —— 正确。

### 11.2 参考项目的调研结论

`git-repo` 下 27 个仓库，逐个查了地理处理方式：

| 仓库 | 做法 | 能否借鉴 |
|---|---|---|
| **spreadgl2.github.io** | `src/lib/format/gazetteer.ts`：多级 gazetteer（country/admin1/region）+ ISO + 别名 + 转写归一 + 国家消歧；数据来自 **Natural Earth** 构建 | ✅ **主参考**（本方案即采用其架构） |
| **MAPLE** | `backend/app/services/location_resolver.py`：geopandas + Natural Earth，`representative_point` 取点 | ✅ 取点方法参考（但方向是坐标→地名，反了） |
| polio-wpv1-phylodynamics | `R/map.R:94`：投影 EPSG:3857 后取 centroid 再转回 | ✅ 佐证「不能直接对经纬度取 centroid」 |
| ggphylogeo | `maps::map.where()` 坐标→国家；`world.cities` 城市点 | ⚠ 反向，且依赖 R `maps` 包 |
| usat_snp | `centroid_sf.r`：对同 locality 的样本点求均值 | ⚠ 仅聚合已有坐标，非推算 |
| PTA / PhyloGeoPlot / 其余 | 无地名→坐标推算，仅消费已有坐标 | — |

**关键发现**：`spreadgl2` 的构建流程**完全公开可复现**
（`github.com/spreadgl2/spreadgl2-gazetteer`，MIT），文档写明：
> 「Use GeoPandas/Shapely `representative_point()` to compute one point
>  guaranteed to lie within each Natural Earth geometry.」
> 「`representative_point()` is used instead of polygon centroid because
>  centroids can fall outside concave or multipart polygons.」

**未复用其 `dist/gazetteer.json`**（那是其项目产物、且捆绑 US Census 覆盖），
改为从 **Natural Earth 原始数据自己重建** —— 保证可复跑、可审计、来源清晰。

### 11.3 数据源与许可（已核实）

- **Natural Earth 1:10m** cultural vectors（admin0 countries / admin1 states-provinces）
- 许可：**public domain**。官方 terms-of-use 原文：
  > "All versions of Natural Earth raster + vector map data found on this website
  >  are in the public domain. You may use the maps in any manner, including
  >  modifying the content and design, electronic dissemination, and offset
  >  printing. ... No permission is needed to use Natural Earth."
- 下载通道（均实测 200）：
  - `raw.githubusercontent.com/nvkelso/natural-earth-vector/master/geojson/`
  - `naciscdn.org/naturalearth/10m/cultural/`
- 采用 **GeoJSON** 而非 shapefile：无需 GDAL/Fiona/pyogrio 等重依赖，
  只装 `shapely` 即可（本环境 shapely 2.1.2）

### 11.4 实现

| 项 | 内容 |
|---|---|
| 构建脚本 | `_geo_build_20260916/scripts/build_gazetteer.py`（可复跑，`--check` 可单独校验） |
| 产物 | `_geo_build_20260916/dist/gazetteer.json`（1852 KB）+ `dist/sources.json`（含源文件 sha256） |
| 条目 | **4665** = country 255 + admin1 4410；覆盖 **236** 个国家 |
| 取点 | admin0 优先用 Natural Earth 人工校订的 `LABEL_X/LABEL_Y`；缺失时 `representative_point()`。admin1 全部用 `representative_point()`（该层 `label_x/label_y` 为空） |
| 别名 | `NAME_ALT`/`name_alt` 拆分 + `name_local` + `gn_name`/`woe_name` + **1010** 个 ASCII 折叠别名 + **2689** 个中文简称别名 |
| 中文支持 | 直接用 Natural Earth 的 `NAME_ZH`/`name_zh`（admin0 258/258、admin1 4589/4596 有值），**无需外部翻译表** |
| 解析层 | `metadata_governance.coord_from_location` 改为**多级查找**；新增 `_load_gazetteer` / `_lookup_gazetteer` / `_split_loc` / `_coord_from_gazetteer` / `_norm_key` / `_fmt_coord` |
| 消歧 | 国家前缀（`United States:Maryland`）→ 优先取该国同名地名。Natural Earth 共 **95 组**重名 admin1，故消歧必需 |
| 降级 | gazetteer 缺失时自动降级为「仅中国城市字典」，不报错、不阻断 |
| 依赖 | **不 import** submission pipeline；**不 import** 参考项目；数据以 JSON 落地 |

### 11.5 中国主权相关规范化（合规处理，重要）

上游 Natural Earth 的标注**与中华人民共和国官方立场不一致**，实测原始属性：

| ADMIN | TYPE | NAME_ZH | SOVEREIGNT |
|---|---|---|---|
| `Taiwan` | `Sovereign country` | **`中华民国`** | `Taiwan` |

构建脚本在生成时**予以规范化**（原始属性保留在 `normalizedFrom` 字段便于核对）：

| 上游 | 本项目输出 | 处理 |
|---|---|---|
| `Taiwan`（Sovereign country） | `Taiwan, China` / `中国台湾` / `CN-TW` | **kind 由 country 降为 admin1**，`countryIso2` 归 `CN` |
| `Hong Kong S.A.R.` | `Hong Kong, China` / `中国香港` / `CN-HK` | 同上 |
| `Macao S.A.R` | `Macao, China` / `中国澳门` / `CN-MO` | 同上 |
| admin1 中 `adm0_a3 ∈ {TWN, HKG, MAC}`（共 40 个） | 全部 `countryIso2='CN'` | 省市一并归入中国 |

解析层另设 `_CN_SUBREGION_KEYS`：港澳台被判为「比国家更细的层级」，
因此 `Taiwan, China` / `Macao, China` 会解析到**真实港澳台坐标**，
而**不会**输出任何国家层级的表述，也不会降级成中国中心点。
验证见 C11d。

### 11.6 过程中实测发现并修掉的 3 个 bug

| # | bug | 症状 | 修法 |
|---|---|---|---|
| 1 | `_GAZETTEER_PATH` 用 `parents[2]` + `"..",".."` 双重上跳 | 路径解析到 `MMPV-RNA\..\..\_geo_build...`，gazetteer 永远加载失败，静默降级 | 改为 `parents[2]` 直接拼 `_geo_build_20260916` |
| 2 | **噪声后缀正则 `[_\-]?(ai\|v1\|v2\|v3)$` 的 `?` 使分隔符可选** | `China:Qinghai` → `Qingh`、`Shanghai` → `Shangh`；**任何以 ai 结尾的地名被截断** | 改为**必须**有分隔符：`[_\-](?:ai\|v[1-9]\|...)$` |
| 3 | `_split_loc` 盲目反转顺序 | `Macao, China` → 反转后首段是 `China`，被解析成国家质心，港澳台降级 | 不反转；改为「由细到粗」再「由粗到细」**两轮**尝试 |

bug 2 尤其危险：它不报错、不抛异常，只是静默给出错误地区的坐标。

### 11.7 行为变更（老师逐条确认）

| 地名 | 旧 `geo_to_latlon` | 新实现 | 老师裁定 |
|---|---|---|---|
| 30 个中国城市 | ✓ | **逐字一致** | 回归硬约束 |
| `China`（无更细层级） | `35.86 N 104.19 E` | **逐字一致** | 零漂移（不引入 Natural Earth 的 `32.50 N 106.337 E`） |
| `China:Gansu` / `China:Qinghai` | 中国中心点（误差上千公里） | **真实省坐标** | ✅「用真实省坐标」 |
| `United States:Maryland` / `:Massachusetts` | **推不出** | 能算 | ✅「采用」 |

**全库对拍**：原实现可解析 **479** → 新实现 **482**，**一条未丢**；
差异**仅限**上述白名单 4 个地名（C6b 硬性判据）。

### 11.8 验证

| 脚本 | 内容 | 结果 |
|---|---|---|
| `verify_metadata_coord_derive.py`（改写判据） | **C6 拆为分层**：C6a 中国城市逐字一致 / C6b 差异仅限白名单 / C6c 一条未丢 / C6d 省级精度提升 / C6e 美国州新增 / C6f 可解析数不减少；**C11 新增全球覆盖**（16 个各洲代表地名 + 坐标合法范围 + 同名消歧 + 中国主权 C11d）；C12 能力边界改写 | **29 PASS / 0 FAIL**（原 18） |
| `verify_metadata_dataquality.py`（更新判据） | Q7b 由「补算 5 条」改为「8 条」（全球查找 +3） | **26 PASS / 0 FAIL** |
| 全量回归 | governance 26 / stage 17 / virphy_bridge_governed 8 / dataquality 26 / coord_derive 29 / fix_decimal_year 32 / task12 17 / fix_complete_deletion 14 / D11 8 / seq_clean 10 / clean_stage 9 / date_precision 4 | **零 FAIL** |

### 11.9 判据变更的说明（改的是判据，不是降低标准）

原 C6/C11 要求与 `geo_to_latlon` **逐字一致**。本轮行为**有意**改变
（老师批准），故判据相应分层 —— 但**同时新增了更强的约束**：

- 原判据只覆盖「450 行 + 30 城市」
- 新判据额外覆盖「全球 16 个代表性地名 + 同名消歧 + 中国主权 + 白名单外的任何差异一律 FAIL」

即：**对旧行为的约束收紧了（差异必须落在白名单内），对新能力加了覆盖**。

### 11.10 能力边界（如实记录）

- **层级精度**：gazetteer 只到 **admin1（省/州）**，没有城市级全球数据。
  城市级精度仍只在中国 30 城（手校字典）。非中国城市名（如 `Chennai`）
  若不在别名内则推不出。
- **覆盖范围**：Natural Earth admin1 为 4410 条（有标准 ISO 3166-2 编码者）。
  无标准编码的条目（如 `CN-X01~` 西沙群岛）未收录。
- **取点语义**：坐标是**几何代表点**，不是人口中心或采样点。
  admin0 多为标注点，admin1 为 `representative_point()`。
- **依赖体积**：gazetteer.json **1852 KB** 随仓库分发；缺失时自动降级。
