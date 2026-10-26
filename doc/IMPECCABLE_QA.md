# Impeccable 设计检测与报告层 QA 记录

- 日期：2026-10-17
- 范围：`reports/html_preview_20261016/` 下 6 份管线报告（① 公共数据、② 预处理、③ 发现、④ 分析、⑤ 系统发育、⑦ 提交）
- 本轮收口：⑤ `virome_phylo_pipeline/utils/report_builder.py`、⑦ `virome_submission_pipeline/report_html.py` + `mmpv_common/templates/submission_report.html.j2`，以及两份预览成品
- 验收口径：impeccable 结构性规则 0 findings；⑤ 52 项测试、mmpv_common 24 项测试全绿

## 1. impeccable 是什么

impeccable 是一套面向 UI 的确定性设计检测规则集与设计技能包（上游：`pbakaus/impeccable`，本机安装在 `%USERPROFILE%\.agents\skills\impeccable`）。它不依赖视觉模型，直接扫描 HTML/CSS/JSX 文本，按规则命中“AI 生成界面”的高频反模式：单侧彩条（`side-tab`）、圆角卡顶部彩条（`border-accent-on-rounded`）、布局属性动画（`layout-transition`）、空 `src` 图片（`broken-image`）、单字体页面（`single-font`）、装饰性编号（`numbered-section-markers`）、破折号滥用（`em-dash-overuse`）等。每条命中都有稳定 rule id，可用于 CI 门禁、项目级忽略配置或单文件内联豁免；退出码 `0` 表示无主要命中，`2` 表示有命中，`1` 表示扫描目标不可读。规则集会随版本演进，因此 QA 记录必须同时写明检测器版本。

## 2. 如何对 6 份管线报告跑检测

### 2.1 任务书里的 Node 检测器命令

```bat
node "%USERPROFILE%\.agents\skills\impeccable\scripts\detect.mjs" reports\html_preview_<date>\*.html
```

本轮基线日期为 `20261016`：

```bat
node "%USERPROFILE%\.agents\skills\impeccable\scripts\detect.mjs" reports\html_preview_20261016\*.html
```

### 2.2 npx / 当前 CLI 用法

```bat
npx impeccable detect reports\html_preview_20261016\*.html
npx impeccable detect --json reports\html_preview_20261016\*.html
npx impeccable detect --quiet reports\html_preview_20261016\*.html
```

本机 skill 已升级到 0.1.15（`scripts/impeccable.cmd` + `bin/windows-x64/impeccable.exe`），**不再随包提供 `scripts/detect.mjs`**。因此当前等价命令是：

```bat
"%USERPROFILE%\.agents\skills\impeccable\scripts\impeccable.cmd" detect reports\html_preview_20261016\*.html
```

或直接 `npx impeccable detect ...`。若必须复现本轮 24 处命中的原始规则集（skill v3.9.1 的 Node 检测器），从上游标签取同一脚本：

```bat
curl -L https://github.com/pbakaus/impeccable/archive/refs/tags/skill-v3.9.1.tar.gz -o skill-v3.9.1.tar.gz
tar -xzf skill-v3.9.1.tar.gz
node impeccable-skill-v3.9.1\skill\scripts\detect.mjs reports\html_preview_20261016\*.html
```

说明：

- `--no-config` 会同时关闭项目配置和内联 `impeccable-disable*` 豁免，只用于审计“豁免前”的原始命中；日常验收不要加。
- `--json` 输出机器可读数组，空数组 `[]` 即 0 findings。
- 新版检测器（当前 CLI / skill v4 系列）新增了 `low-contrast`、`tiny-text`、`undersized-ui-text` 等可访问性规则。当前 CLI 扫描 ⑤⑦ 会报这批配色/字号问题；它们不属于本轮 24 条结构性命中，建议作为下一轮独立配色审计处理，不要用整文件豁免掩盖。

## 3. 本轮 24 处命中与分类处置

基线：skill-v3.9.1 `detect.mjs` 对 6 份报告扫描，共 24 处。处置分为两类：

- **修复**：改结构或动画实现，命中消失，不需要豁免。
- **豁免**：内容本身是真实数据/科学文案，改写会损失信息，使用单文件内联注释精确豁免。

| # | 报告 | 命中（规则 / 位置） | 处置 | 说明 |
|---:|---|---|---|---|
| 1 | ① | `layout-transition` line 21 `transition: max-height` | 修复 | 删除 `max-height` 过渡，折叠改为瞬时开合（保留 `.nav-sub.collapsed` 状态类） |
| 2 | ① | `em-dash-overuse` 14 处 | 豁免 | CJK 科学文案中的破折号是语义标点，不改写 |
| 3 | ① | `numbered-section-markers` 序列 10,11,12 | 豁免 | 编号是真实阶段 ID（S10/S11/S12），不是装饰性章节序号 |
| 4 | ② | `side-tab` line 33 stage-pass/fail/skip | 修复 | 左侧 4px 彩条改为 1px 全边框 + 状态色浅底 |
| 5 | ② | `side-tab` line 53 tr-fail/tr-skip 首列 | 修复 | 行首列彩条移除，改为整行浅色背景（`.tr-pass/.tr-fail/.tr-skip`） |
| 6 | ② | `broken-image` line 83/84 lightbox | 修复 | 空 `src` 改为 1x1 透明 GIF data URI，真实图片仍由 `openLightbox()` 注入 |
| 7 | ③ | `side-tab` line 122 stage-pass/fail/skip | 修复 | 同 ②：1px 全边框 + 状态色浅底 |
| 8 | ③ | `side-tab` line 148 tr-fail/tr-skip 首列 | 修复 | 同 ②：整行浅色背景 |
| 9 | ③ | `side-tab` line 167 可见左边条 | 修复 | 移除该 3px 可见左边条；状态由圆点/徽章承载 |
| 10 | ③ | `border-accent-on-rounded` line 157 | 修复 | 移除圆角卡顶部 4px 彩条，保留 1px 边框 |
| 11 | ③ | `broken-image` line 225/226 lightbox | 修复 | 1x1 透明 GIF data URI |
| 12 | ③ | `em-dash-overuse` 7 处 | 豁免 | CJK 科学文案 |
| 13 | ③ | `numbered-section-markers` 序列 01–06 | 豁免 | 编号是管线阶段标识（S01–S06），不是装饰序号 |
| 14 | ④ | `side-tab` line 116 skip-chip | 修复 | `border-left:3px solid var(--warning)` 改为 1px 全边框 `#fde68a` + `#fffbeb` 浅底 |
| 15 | ④ | `layout-transition` line 229 `max-height` | 修复 | 删除过渡，折叠瞬时完成 |
| 16 | ④ | `broken-image` line 499/501 lightbox | 豁免 | `src` 由 `openLightbox()` 运行时注入；用 `impeccable-disable-next-line broken-image` 精确豁免占位 `<img>` |
| 17 | ④ | `em-dash-overuse` 15 处 | 豁免 | CJK 科学文案 |
| 18 | ④ | `numbered-section-markers` 序列 09–12 | 豁免 | 真实阶段编号 |
| 19 | ⑤ | `side-tab` line 47 `.asset-note` | 修复 | 左侧 3px 橙条改为 1px 全边框 `#f0d9a8` + `#fdf6e3` 浅底；同时移除侧栏 `.sb-item` 的 pass/fail/skip 左边条，状态只保留 `sb-dot` 圆点 |
| 20 | ⑤ | `single-font` line 9 全页 Arial | 修复 | 标题/hero 改衬线 `'EB Garamond', Georgia, 'Times New Roman', serif`，正文保持系统 sans；零外部字体资源 |
| 21 | ⑤ | `numbered-section-markers` 序列 10,11,12 | 豁免 | 模块顺序与病毒编号是真实序列（Data→…→Report），不是装饰序号 |
| 22 | ⑦ | `layout-transition` line 163 `.progress-bar .fill` | 修复 | 该条形在 ⑦ 不使用，但共享 `base_css` 的死规则仍会被扫描；预览原位改为 `transition: transform`，`report_html.py` 在渲染后做同字符串替换。`report_style.py` 是共享样式源，不在本轮文件白名单，后续可把同一改动上收 |
| 23 | ⑦ | `layout-transition` line 184 `#completeBar` | 修复 | 内层宽 100% + `transform: scaleX(0)` + `transform-origin:left`；`updateStats()` 改为 `style.transform='scaleX('+(pct/100)+')'`，不再写 `width` |
| 24 | ⑦ | `em-dash-overuse` 7 处 | 豁免 | CJK 科学文案；生成物表格里的 09/10/11 是数据值，同注释一并豁免 `numbered-section-markers` |

补充：⑦ 的共享 `base_css` 里还有一条 `.sidebar a { border-left: 3px solid transparent }` 及 active 态 `border-left-color`。⑦ 不渲染侧栏，这是死规则，但它会在静态扫描里留下 side-tab 形状；`report_html.py` 在渲染后把这两处字符串一并移除，预览 07 同步。`report_style.py` 是共享样式源，不在本轮白名单，后续可把 `.progress-bar .fill` 和 `.sidebar a` 两处改动上收到该文件。

### 3.1 ⑦ 进度条比例核对

`updateStats()` 的计算保持原口径：

```js
const pct = total > 0 ? Math.round((total - missing) / total * 100) : 100;
document.getElementById('completeness').textContent = pct + '%';
document.getElementById('completeBar').style.transform = 'scaleX(' + (pct / 100) + ')';
```

- `total` 是 `#srcTable` 中带 `required-col` 的单元格数；`missing` 是其中 `isMissing()` 为真的数量。
- 填充比例 `scaleX = (total - missing) / total`（四舍五入到整数百分比后再除以 100），与“Complete”卡片显示的 `pct` 完全一致。
- `missing=0` 时 `scaleX(1)` 满格；`missing=total` 时 `scaleX(0)` 空格；`total=0` 时按 100% 处理，保持原行为。
- 只改 `transform`，`background` 仍按 90/70 阈值切绿/黄/橙，不触发 layout。

## 4. 豁免注释机制

impeccable 支持三种内联豁免，注释语法不限（HTML/CSS/JS/`#` 均可识别）：

```html
<!-- 整个文件：规则 id 逗号分隔，冒号或 -- 后写理由 -->
<!-- impeccable-disable numbered-section-markers,em-dash-overuse: 编号是真实阶段 ID，正文为 CJK 科学文案 -->

<!-- 仅下一行 -->
<!-- impeccable-disable-next-line broken-image -- src is assigned by openLightbox JS -->

<!-- 仅同一行 -->
<img src="..." /> <!-- impeccable-disable-line broken-image -->
```

规则：

- 规则 id 只含单连字符；多个 id 用逗号或空格分隔。
- 省略规则 id 或写 `*` 等于豁免全部规则；不要用它来“过验收”，那会把真实缺陷一起藏掉。
- 理由写在 `:` 或 `--` 之后，检测器解析时丢弃，只用于代码审阅。
- 默认 `impeccable detect` 会读取并执行这些豁免；`--no-config` 与 `--no-inline-ignores` 会绕过它们，适合做“豁免前”审计。
- 仓库级、长期存在的例外优先写 `.impeccable/config.json` 的 `detector.ignoreRules/ignoreFiles/ignoreValues`；内联注释只用于“必须随单文件分发”的导出报告。本项目的报告是自包含 HTML，会离开仓库单独发送，所以采用内联豁免。

本轮落点：

- ⑤：`report_builder.py` 生成页 `<head>` 内直接输出 `<!-- impeccable-disable numbered-section-markers,em-dash-overuse: module order & CJK scientific copy -->`；预览 05 第 6 行同字符串。
- ⑦：`submission_report.html.j2` extends 的 `base.html.j2` 没有 head block，故由 `report_html.py` 在 `render()` 后把 `<!-- impeccable-disable numbered-section-markers,em-dash-overuse: CJK scientific copy & sequence IDs are data values -->` 插到 `<head>` 后；预览 07 第 4 行同字符串。
- ④：占位 `<img>` 上方已有 `impeccable-disable-next-line broken-image`，属于“动态 src”场景的精确豁免。

## 5. 后续新增报告模板的设计纪律

新增或修改报告模板时，先过下面清单，再跑检测器：

1. **禁单侧彩条**：卡片/提示条/列表项不要用 `border-left: 3px/4px solid <状态色>`。要状态色就用 1px 全边框（状态色浅阶）+ 浅底，或只留圆点/徽章。
2. **禁圆角卡顶部彩条**：`border-top: 4px solid` 与 `border-radius` 不要同时出现；要么去彩条，要么去圆角。
3. **禁空 `src` 图**：不要输出 `src=""`、无 `src`、`src="#"` 的 `<img>`。灯箱/懒加载占位用 1x1 透明 data URI，或先不渲染 `<img>`。
4. **动画只用 transform/opacity**：不要 `transition: width/height/padding/margin/max-height`。进度条用 `transform: scaleX()`；折叠面板用 `grid-template-rows` 或不加动画。动画时长控制在 150–300ms。
5. **字体成对出现**：标题/hero 用 display 字体（本项目用 `'EB Garamond', Georgia, 'Times New Roman', serif`，零外部资源），正文用系统 sans。不要全页只有一个 `Arial`；如果确实单字体，必须靠字重/字号拉开层级并在评审里说明。
6. **编号只给真实序列**：`01/02/03` 只用于真实阶段、模块、样本或病毒 ID；纯装饰性章节序号不要用。真实编号需要保留时，用内联 `impeccable-disable numbered-section-markers` 并写清理由。
7. **CJK 科学文案的破折号**：技术 prose 里的 `—` 是语义标点，不要为了检测器改成逗号/括号；按文件豁免 `em-dash-overuse`。
8. **报告必须自包含**：不引外部字体、CSS、JS；字体走系统栈，图片走 data URI 或相对链接。
9. **提交前跑检测**：`node ...detect.mjs <报告>.html` 或 `npx impeccable detect <报告>.html`，结构性规则应为 0；再用 `--no-config` 复核豁免是否最小化。
10. **升级检测器要重新基线**：规则集会变（例如新版新增 `low-contrast`/`tiny-text`）。升级后先记录新命中，再决定修复还是新增豁免，不要沿用旧版本的“0 findings”结论。

## 6. 本轮验证结果

```text
python -m py_compile virome_phylo_pipeline/utils/report_builder.py virome_submission_pipeline/report_html.py
→ py_compile OK

python -m pytest virome_phylo_pipeline/tests -q
→ 52 passed, 5 warnings in 50.59s

python -m pytest mmpv_common/tests/test_report_template.py mmpv_common/tests/test_process_report.py -q
→ 24 passed in 1.28s

node impeccable-skill-v3.9.1/skill/scripts/detect.mjs reports/html_preview_20261016/05_phylo_report.html reports/html_preview_20261016/07_submission_report.html
→ 无输出（0 findings），exit 0

node impeccable-skill-v3.9.1/skill/scripts/detect.mjs --json reports/html_preview_20261016/05_phylo_report.html reports/html_preview_20261016/07_submission_report.html
→ []
```

补充：本轮还用同一检测器扫描了由修改后源码重新生成的 ⑤/⑦ 临时报告，均为 0 findings，确认修复不只存在于预览成品，生成器本身也已干净。
