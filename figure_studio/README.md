# Figure Studio · 本地组图工作台

把已做好的单面板图（PNG / JPG / PDF / **SVG**）**拖拽自由组合**成一张大图，
所见即所得地调整位置、尺寸、面板字母、题注，然后**一键导出印刷级 PNG + PDF**。
布局存成 JSON —— 出图之后想再改，改的是**数据**，不是代码。

```
┌──────────┬────────────────────────────┬──────────────┐
│  面板库   │          画布               │   属性面板    │
│ 缩略图网格 │  （拖拽 / 缩放 / 吸附对齐）  │ 画布/面板/样式 │
└──────────┴────────────────────────────┴──────────────┘
```

**两级编辑，能力不同**（重要，决定了「什么图能改什么」）：

| 源格式 | 编辑层级 | 能改什么 |
|---|---|---|
| PNG / JPG / PDF | **图像级**（`composer.apply_adjust`） | 亮度/对比度/饱和度/灰度/反相/阈值/色调分离/底色/着色/**精确颜色替换（带容差）**/裁边 |
| **SVG** | **矢量级**（`svg_edit.apply_svg_edit`，改 XML DOM）+ 图像级 | 上面全部，外加：**语义色替换（精确匹配、无抗锯齿偏色、无容差）** · **图例隐藏/移动/缩放** · **文本内容替换** · **全局字号缩放 / 字体族 / 文字颜色 / 线宽缩放** · 按选择器隐藏任意元素 |

SVG 的算子顺序是**契约**：先矢量编辑（改 XML）→ 再栅格化 → 最后图像级调色。
PDF 导出时 SVG 面板**以矢量嵌入**（Form XObject），文字仍是可检索的 PDF 文本，
不会变成位图 —— 导出提示里会写 `PDF SVG vector embed: N panel(s) as vector`。

---

## 快速开始

```bash
C:/Python312/python.exe D:/桌面/延伸基因组/MMPV-RNA/figure_studio/app.py --port 8770
```

浏览器打开 <http://127.0.0.1:8770>

也可以双击 `start.bat`（等价命令，会自动 `cd` 到本目录）。

依赖：Python 3.12 + `flask` / `matplotlib` / `Pillow` / `numpy` / `PyMuPDF`（`C:\Python312` 已具备）。
不需要 npm、不需要联网。

---

## 界面

| 区域 | 能做什么 |
|---|---|
| **左：面板库** | 填目录 → 扫描。递归找 PNG/JPG/PDF/**SVG**，出缩略图。搜索框按**文件名 / 目录**过滤。单击加入画布，或**直接拖到画布上的目标位置** |
| **中：画布** | 拖动移动、右下角方块缩放、方向键微调（`Shift` 加速）、`Delete` 删除。拖动时按 **0.05 in 网格 + 其他面板边缘/中心 + 画布边距**吸附，并显示粉色对齐参考线。**`Shift`/`Ctrl` 点击加选**（多选后拖动 = 整组移动，方向键 = 整组微调，`Delete` = 批量删除） |
| **右：属性面板** | 画布尺寸/DPI、选中块的精确 X/Y/宽/高、面板字母、题注、层级；**调色**（滑块/开关/色板/颜色替换表）；全局字母/题注样式与颜色；总标题（含颜色）；导出与产物列表 |
| **工具栏** | **↶撤销 / ↷重做** · 横排一行 / 纵排一列 / 网格 2 列 / 自动字母 · 左对齐 / 水平居中 / 右对齐 / 顶对齐 / 垂直居中 / 底对齐 / 水平均分 / 垂直均分 · 统一高 / 统一宽 / 1:1 尺寸 / **裁到内容** |

### 多选与撤销（作用域 + 后悔药）

**多选**：`Shift` 或 `Ctrl` 点击面板 = 加选/取消。集合成员画**虚线**蓝框，
**主选中**（属性面板、调色、精修的作用对象）画**实线粗框**并带缩放手柄 ——
手柄只给主选中，因为多选状态下 `resize` 只作用于一个面板（整组缩放未实现，画全了会误导）。

多选决定**排版工具的作用域**：选中 ≥2 个时，对齐/均分/统一尺寸/自动字母只作用于这 2 个；
没选或只选 1 个时作用于全部。

> ⚠ 这条语义在代码注释里一直写着，但早期 `state.sel` 是单值，
> `filter(i => i.id === state.sel)` 最多命中 1 个 → 那个分支**永远进不去**，
> 等于死代码：所有工具只能作用于全部面板。`selSet` 才让它真正成立。

**撤销/重做**：整状态快照栈（上限 80 步）。粒度按「一次操作」——拖动/缩放松手记一步，
输入类（属性框、调色滑块、方向键连按）防抖 350ms 合并成一步。

### 快捷键

| 键 | 作用 |
|---|---|
| `Shift` / `Ctrl` + 点击 | 加选 / 取消某个面板 |
| `Ctrl+A` | 全选面板（焦点在输入框里时是选文字，不抢） |
| `Ctrl+Z` / `Ctrl+Shift+Z` / `Ctrl+Y` | 撤销 / 重做 |
| `Ctrl+D` | 复制（多选时整组复制，保持相对位置） |
| `Delete` / `Backspace` | 删除（多选时批量删） |
| 方向键 | 微调 0.01 in（`Shift` = 0.1 in；多选时整组移动） |
| `Ctrl+E` | 用 svguest 精修选中的 SVG 面板 |

⚠ 带修饰键的快捷键**先判**、裸键后判，且输入框内一律早退 —— 否则在属性框里
按 `Ctrl+A` 会把画布全选、按 `Delete` 会删面板。

### 关键：排版是按「资产比」反解的，不会形变

「横排一行」的算法与项目 v23 主图同源：

```
可用宽 AW = 画布宽 − 2×留边 − 间距×(n−1)
行高   h  = AW / Σ(1/aspect_i)          # 全行等高
列宽   w_i = h / aspect_i
```

⇒ 每个面板的**面板比恒等于源图比**，`imshow(aspect="auto")` 不会静默拉伸，
也不会留下拼接空白。纵排、网格同理。

### 「显示字母」开关

⚠ 本项目 `make_discovery_sci_figs.py::build_panels()` 出的单面板图**已经烘焙了面板字母**
（`panel_L.png` 左上角自带一个 "L"）。这种源图再叠平台字母会**双字母**
→ 把右侧「显示字母」取消勾选即可（或清空单个面板的字母字段）。

### 图型无关 —— 它只认像素，不认内容

合成器是**纯贴图**：读图 → 按宽高比放进你给的位置 → 叠字母/题注。
**不管面板里画的是条形图、环形树、网络图、桑基、热图还是折线，走的是同一条路。**

已用本项目档案库的**真实成品图**实测过 8 种图型（`demo_types.py`）：

| 字母 | 图型 | 源图 | 像素 |
|---|---|---|---|
| A | 环形 ML 系统发育树 | `03_MLtree_circular/成品图/7_MLtree_circular_host.png` | 1374×1374 |
| B | NeighborNet 网络图 | `04_NeighborNet/成品图/8_NeighborNet_host.png` | 2796×2787 |
| C | 桑基图 | `07_Fig_GojiVirome_NoveltySankey/成品图/…` | 2796×2520 |
| D | UpSet 交集图 | `09_UpSet_工具分类宿主/成品图/…` | 2796×1800 |
| E | 柱状（**R 出的调色板 PNG**） | `13_HostPrediction_宿主预测/成品图/…` | 2796×1800 |
| F | SDT 相似度热图 | `S10_SCI_Discovery/Fig_GojiVirome_Potexvirus_SDT.png` | 2327×2425 |
| G | 分布折线 | `…_SDTdist.png` | 3025×1711 |
| H | 覆盖度折线 | `…_Coverage.png` | 2379×617 |

```bash
C:/Python312/python.exe figure_studio/demo_types.py
# → outputs/demo_types.png|pdf（9.32 in 通栏）+ layouts/demo_types.json
```

结果：**零形变告警**，八块各按源图比铺满，题注完整。
（其中面板 E 是 `mode='P'` 的调色板 PNG —— 合成器在 `load_panel()` 里一律
`convert("RGB")`，正是本项目踩过的「白底变亮黄」那个坑。）

想直接看这张组合图：启动平台 → 「加载布局」→ 选 `demo_types`。

---

## 调色（图像级）

面板是**已渲染的栅格图**，所以调色只能是**图像级**（滤镜 / 颜色替换 / 底色处理），
**不是**语义级重绘 —— 平台不会去改「某根柱子」的颜色，只会改**像素**。
（如果面板源是 **SVG**，则还有一套**矢量级**编辑可用，见下一节。）

> ⚠⚠ **先选中面板，调色才生效**：所有调色控件只在「画布上点选了一个面板」
> （出现蓝框）时有反应；没选中时它们是 disabled、点击零反馈。所以：
> 点画布空白处取消选中后，右侧会出现**空态提示 + 遮罩**，点遮罩会告诉你去点面板；
> 「加载布局」也会自动选中第一个面板，就是为了避开这个坑。

每块面板各带自己的 `adjust`，**预览与成品走同一条 `apply_adjust()`**，所见即所得。

| 类别 | 控件 |
|---|---|
| 基础 | 亮度 · 对比度 · 饱和度 · 锐度 · 伽马 · 色相旋转 |
| 转换 | 灰度 · 反相 · 自动对比度 · 二值化（阈值） · 色调分离（位数） |
| 底色 | 透明底色（alpha 合成） · 近白替换 + 容差（去灰底） |
| 着色 | 单色着色（目标色 + 强度） |
| **颜色替换** | 把源图某个色精确换成另一个色（带容差） |
| 几何 | 裁白边（⚠ 会改变面板比例，自动重算宽高） |

**算子顺序是契约**（`apply_adjust()` 里固定）：颜色替换 → 反相/灰度 → 自动对比度 →
二值化/色调分离 → 亮度/对比度/饱和度/锐度 → 伽马 → 色相 → 着色 → 近白替换 → 裁边。
颜色替换**必须最先**：要在原图**精确色值**上匹配，被亮度/伽马改过就对不上了。

### 颜色替换：本项目最实用的那个

1. 选中面板 → 点 **「取主要颜色」** → 平台按像素占比列出该图主色（含占比）。
2. 把右侧色块改成目标色 —— 建议直接用上方**项目调色板**点选。
3. `from == to` 的行**不生效**（是「占位未填」状态），可以先全取出来再逐个改。

⚠⚠ **背景白会被自动跳过**（不生成替换行）。实测 matplotlib 图白底占比 **96%**，
它是「最大主色」，若生成 `white→white` 行，用户一把目标色改成彩色**整张图瞬间变色**，
看起来像改崩了。要处理背景请用上面的**「透明底色」/「近白替换」**；
真要把白色换成别色，点「+」手加一行即可。

**为什么它对你们特别有用**：R 的 `ggsave` 出的图常带**抗锯齿偏色**
（实测某面板主色是 `#0172b2`，比项目 `#0072B2` 差 1 个色阶）——
容差 6–8 的替换就能把它们**拉回项目色板**，不用回去改 R 脚本重出图。

⚠ 但**α 混合出来的颜色替换不了**：桑基节点的 `#add1e6` 是 `#56B4E9` 与底色的混合色，
不是「某个色值」，容差再大也匹配不准 → 这种图要用**整体着色（tint）**，别用 remap。

### 项目调色板

色块条从 `virome_phylo_pipeline/utils/sci_palette.py` **动态加载**（按显式路径
`spec_from_file_location`，不走 `import` —— 该模块名在多个管线里重名）：
`OKABE_ITO` 9 色（固定序）+ `OKABE_BRANCH` 5 色 + 6 个**语义色**（圆形标记）。

点色块 → 填入**最近点过的那个颜色框**。⚠ 语义色是**固定含义**的
（`C_UP`=上升/显著…），**别当普通色板随手换色**。

### 预览 vs 成品

预览是**降采样后**再调色，`remap` 的容差在抗锯齿边缘会漏掉一部分过渡像素，
所以**预览边缘看着不如成品干净**。成品在**原始分辨率**上调色 —— 以成品为准。

### 题注自动折行

题注（`fig.text`）**按面板宽自动折行**。原因：`fig.text` 溢出既不报错、也不被
`tight_layout` 搬动，长题注会**直接压在隔壁面板的题注上**（这个 bug 真实发生过）。
折行宽度可用 `caption_max_w_in` 覆盖。

### 调色对照图

`demo_adjust.py` 用**与 `demo_types.py` 完全相同的版式**（同样 8 块、同样的行列尺寸），
只改每块的 `adjust` —— 可以直接和 `demo_types.png` 对照着看每个算子的效果：

```bash
C:/Python312/python.exe figure_studio/demo_adjust.py
# → outputs/demo_adjust.png|pdf + layouts/demo_adjust.json
```

算子按**各面板的真实问题**选，不是随机套：A 容差替换拉回偏色 · B 灰度 ·
C 整体着色（α 混合色替换不了） · D 色相 +150° · E 近白替换 + 自动对比度 ·
F 亮度/对比度加深 · G 红色 → 项目 vermillion · H 青色 → 项目 sky blue。

---

## SVG 矢量编辑（颜色 / 图例 / 文字，改的是 XML 不是像素）

选中一个 `.svg` 面板，属性面板会出现**「SVG 矢量编辑」**块。点
**「读 SVG 结构」**后，平台用 lxml 解析这份 SVG，列出它的语义色 / 文本 / 组 / 图例，
之后所有改动都发生在 **XML DOM 层**：

| 能力 | 控件 | 说明 |
|---|---|---|
| 语义色替换 | 「读语义色」→ 每行一个色 → 改右侧目标色 + scope（all/fill/stroke） | **精确色值匹配**，无抗锯齿偏色、无容差；白/近白底色行默认跳过 |
| 图例 | 动作（隐藏 / 移动 / 缩放）+ dx/dy/字号倍率 + 选择器 | 作用于**所有命中的图例组**（genome tracks 有 `legend_1`+`legend_2`，早年只改一个 → 看着像没生效） |
| 文字 | 「读文本清单」→ 每行一条文本 → 改成新内容 | 精确匹配整段文本后替换；长度不同不报错也不重排（矢量文本没有自动布局） |
| 全局文字 | 字号倍率 · 字体族 · 文字颜色 · 线宽倍率 | 字号倍率对**所有**文本节点生效（含刻度、题注内的文字） |
| 隐藏元素 | 选择器（`#id` / `.class` / `tag` / `tag.class`）+ 「隐藏」 | 图例、刻度、某个标注都能藏；chips 可点 × 撤销 |

### 三个必须知道的坑（都是实测踩出来的，零报错）

1. **PyMuPDF 的 SVG 渲染器不解析 CSS `font` 简写里的字号** ——
   `font: 7.5px 'DejaVu Sans'` 改字号后栅格结果 **0 差异像素**。
   所以 `apply_svg_edit` 无条件先把 `font` 简写拆成 `font-size` + `font-family`。
   这同时是**基线保真度**问题：不拆时所有文字按默认字号渲染（实测 49k 差异像素）。
2. **`display: none` 被同一个渲染器忽略** —— 隐藏元素必须**从树里移除**
   （`parent.remove(el)`），实测移除 = 28803 差异像素、`display:none` = 0。
   浏览器里两种都生效，所以对浏览器交付的 SVG 没有副作用。
3. **SVG 尺寸单位不统一** —— matplotlib 写 `width="509.76pt"`，
   R/svglib 写裸 `width="920"`（px）。fitz 把裸 px 当成 pt → 物理尺寸虚胖 4/3。
   唯一真相源是 `svg_edit.svg_size_pt()`（px→pt 按 0.75 折算），
   `panels._meta_of` / `composer.load_panel` / `inspect_svg` 三处都从它取，
   selftest 有逐文件断言兜住（`SVG size 口径一致`）。

### 预览与栅格化

预览按目标像素长边**反解 dpi**（`svg_edit.fit_dpi_for`）。早年直接传
`dpi=max_px`：一张 12.8 in 宽的 SVG 在 1600 dpi 下要 20k×18k px ≈ 1 GB →
MemoryError → `/api/img` 静默 404（面板库卡片裂图）。selftest 现在对 4 种
SVG 逐个打预览端点。

### 演示

```bash
C:/Python312/python.exe figure_studio/demo_svg.py
# → outputs/demo_svg_mixed.png|pdf + layouts/demo_svg_mixed.json
```

4 个真实 SVG（matplotlib ×2、R/svglib ×1、0 文本的桑科图 ×1）+ 1 个真实 PNG
拼 10×7.909 in 通栏图：A 行改色 + 文字放大 1.35× · B 行两图例全隐藏 + 文字 1.2× ·
C 行 R 出的图改色 + 文字 1.25× · D 行 PNG 走图像级灰度。PDF 里 3 个 SVG 面板
是**矢量**（1361 条路径、SVG 文字成为可提取的 PDF 文本），PNG 面板仍是栅格。

---

## svguest 精修（交互式矢量编辑器，第三方组件）

上面的「矢量编辑」是**批量语义编辑**：按色值/图例/字号成批改 DOM，适合统一风格。
要**手动精修单张图**（点选、拖锚点、画连接器、吸附对齐、逐字改文本），
用右侧检查器的 **「精修（svguest 交互式）」** 区块 —— 它把 SVG 面板交给一个完整的
交互式矢量编辑器（[svguest](https://github.com/JieqiangHe/svguest)，GPL-3.0）。

```
选中 .svg 面板 → 「精修选中面板」（或 Ctrl+E）→ 平台内全屏浮层打开编辑器（图已载入）
  → 手动改 → Ctrl+S /「保存回工作台」→ 写回原文件 → 工作台预览自动刷新
  → 顶栏「关闭 ✕」收起（焦点在平台上时 Esc 亦可）
```

### 内嵌，不是弹窗

编辑器是**同源 iframe**（`#sgLayer > iframe#sgFrame`），不是 `window.open` 出来的独立窗口。
弹窗的老问题：可能被浏览器拦、保存后要手动切回工作台、关掉窗口还得回来刷新。
改成内嵌后一次点击就进出，且关浮层时把 iframe 卸回 `about:blank`，上游的全局监听、
IndexedDB 连接和定时器随之释放（不然每开一次就多留一份）。

⚠⚠ **载入后必须把键盘焦点推进 iframe**：点的是工作台上的按钮，焦点本来留在平台，
`Ctrl+S`、`m`/`v` 这些工具快捷键一个都到不了编辑器 —— 表现是「画不出东西、保存没反应」，
而且**零报错**。`openSvguest()` 在 iframe 的 `load` 事件里 `contentWindow.focus()` 兜住。

### 三档布局：全屏 / 右侧分屏 / 浮动窗口

顶栏「布局」按钮循环切换（双击标题栏 = 全屏 ↔ 分屏）：

| 档 | 用途 |
|---|---|
| **全屏** | 专心改图，工作台完全让位 |
| **右侧分屏**（58vw） | **左边工作台仍可见可点** —— 一边看画布/面板库，一边改图 |
| **浮动窗口** | 拖标题栏移动、右下角可缩放（CSS `resize:both`），四边不贴屏 |

选择记在 `localStorage`（`fs.sgMode`），下次打开沿用。

⚠⚠ **切档时必须把 float 写下的 inline 定位全清掉**（`left/top/right/bottom/width/height`
**六个**）。漏掉 `right/bottom` 时，`inset:0` 被 inline 的 `auto` 覆盖，元素从
「撑满屏幕」退化成 shrink-to-fit（实测 1680 px → 444 px），**零报错**，只是"看着不对"。
`uitest.py` 的 8b 段就是量这个：切回全屏后 `left==0` 且宽度回到原值。

⚠⚠ **焦点在 iframe 里时，Esc 归上游**，平台不抢。实测上游把每一次 `Escape` 都
`preventDefault`（取消当前工具 / 取消选择 / 退出文本编辑）。早先做过「Esc 关闭浮层」，
结果是用户想「取消选择」却把整个编辑器关掉，已撤除。退出只走顶栏按钮
（平台侧 Esc 仍可用 —— 那时焦点不在 iframe 里，不存在歧义）。

### 能做什么

| 入口 | 行为 |
|---|---|
| **精修选中面板**（Ctrl+E） | 把当前选中的 SVG 面板交给 svguest |
| **全部 SVG 面板（N）** | 把画布上所有 SVG 面板**一次性作为标签页**打开；浮条按钮自动变成「全部保存」，一键写回各自的原文件（svguest 自身上限 20 个标签） |
| 浮条 **PNG 300 / 600** | 用 svguest 的渲染器出图，**直接落到原图同目录**（`<名>_edit_300dpi.png`），不经过浏览器下载；若该目录就是当前面板库根，会自动重扫，导出的图立刻能当面板拖进画布 |
| 浮条 **JPG 300** | 同上，出 `<名>_edit_300dpi.jpg`（dpi 写进 JFIF 头） |
| 浮条 **PDF** | 出**矢量** PDF（`<名>_edit.pdf`）—— 由服务端用 PyMuPDF 生成，文字**可提取**、尺寸与 SVG 声明一致。⚠ 上游原生的 `Ctrl+P` 是「iframe + `window.print()`」，会弹打印对话框、还得手动「另存为 PDF」，这里已接管 |
| 浮条 **SVG 副本** | 另存一份 `<名>_edit.svg` 到同目录（保留原文件不动） |
| 浮条 **查隐藏信息** | 调起 svguest 自带的水印检查器（查隐藏的 AI/来源元数据） |
| 浮条 **另存本地** | svguest 原生保存（`showSaveFilePicker`），行为不变 |
| 检查器 **备份历史** | 列出该面板每次写回/回滚前的备份，点「回滚」即可还原（⚠ 回滚前会先把当前版本再备份一次，所以回滚本身也能撤销） |

### 为什么接管 `A` 表，而不是拦键盘

上游源码里的注释原话是「Buttons, menus, keys and the palette all call
`A[data-a]`」，主键盘处理器也确实写着：

```js
s: e.shiftKey ? A.svgAs : A.svg,     // Ctrl+S / Ctrl+Shift+S
p: e.shiftKey ? A.imp   : A.pdf,     // Ctrl+P
```

也就是说**菜单点击、快捷键、`Ctrl+K` 命令面板三个入口都经 `A` 表**。
早先只在 `keydown` 的 capture 阶段拦 Ctrl+S，只能挡住键盘：
从 File 菜单点「Export SVG」会去调 `showSaveFilePicker` 弹出文件选择器，
与「保存回工作台」的语义打架。现在改 `A.svg` / `A.pdf` 两个属性即可一次覆盖全部入口。

### 怎么接的（零改动移植 + 可重复构建）

上游 `app.html` **一行未改**。`static/vendor/svguest.html` 由
`tools/build_vendor.py` 自动合成：

```
vendor_src/app.html  +  vendor_src/bridge.css  +  vendor_src/bridge.js
```

svguest 的主 script 是**经典脚本**（既非模块也非 IIFE），顶层 `function` 声明挂在
`window` 上、顶层 `const/let` 落在共享的全局词法环境里，所以桥接可以直接调用
`readFile` / `svgText` / `A`（动作表）/ `TABS` / `setTheme`，无需上游暴露 API。

| 环节 | 实现 |
|---|---|
| 打开 | 工作台 `POST /api/svg/session` 拿到短 token → 编辑器在 `#sgLayer` 的**同源 iframe** 里打开，hash 只带 `#s=<token>&theme=light` → 桥接取清单后逐个 `readFile()`（它内部会自动开新标签）。token 每次都是新的，所以赋 `iframe.src` 必然是一次全新导航，不复用上一次的脏状态 |
| 保存 | 接管 `A.svg`：活动标签走 `svgText()`，后台标签直接序列化快照（等价于切过去再存，但不用真切换）→ `POST /api/svg/writeback(_batch)`；成功后同步上游的 `dirty` 状态，否则关标签会一直问要不要保存 |
| 备份 | 写回前把旧版存到 `outputs/_svguest_backup/`（⚠ **不放在面板目录**，否则 `.svg` 备份会被扫描当成面板） |
| 导出直通 | **常驻 `dl` 分流器**按 MIME 认：`image/png`→`png`、`image/jpeg`→`jpg`、`image/svg`→`svg` 走 `POST /api/svg/export`；`application/json`（FigKit 模板）、`application/zip`（Kept 下载）**原样交给浏览器**——按调用点分流会把这些也吞掉 |
| 矢量 PDF | 接管 `A.pdf` → `POST /api/svg/export {kind:'pdf', svg}` → 服务端 PyMuPDF 转 PDF 落盘 |
| 刷新 | 保存成功 `postMessage` 通知工作台（内嵌时发给 `window.parent`，独立窗口时发给 `window.opener`，两条都通） → 重取预览（缓存键含 mtime，自动失效）+ 刷新备份历史 |

**升级上游**：换掉 `vendor_src/app.html`，重跑构建即可（脚本会先剥掉旧桥接块）：

```bash
C:/Python312/python.exe figure_studio/tools/build_vendor.py --check   # 只校验
C:/Python312/python.exe figure_studio/tools/build_vendor.py           # 写出
```

### 五个实测踩出来的坑

1. ⚠⚠ **`Content-Type: image/svg+xml` 的 GET 在本机 Chrome 里会变成合成 204**
   （无 Content-Type、无 body、**零报错**）→ 桥接那边表现为「提示已载入但画布是空的」，
   再保存就把整张图覆盖成空白。`image/svg` / `text/plain` 正常（`/api/file` 一直能用
   就是因为 Windows 注册表把 `.svg` 猜成 `image/svg`）；POST 的 `/api/svg/source`
   用 `image/svg+xml` 反而没事 → 是「GET + 该 MIME」的组合问题。
   所以 `/api/svg/raw` 必须写 `mimetype="image/svg"`。
2. ⚠⚠ **`<style>` 里只能有 CSS 注释**。曾经在 `<style>` 开头混进一行 Markdown 标题
   （`# Figure Studio ⇄ svguest 桥接样式…`），浏览器把 `#` 当成 ID 选择器，
   紧接着的整条 `#fsBar{…}` 被当作无效选择器**静默丢弃** —— 浮条失去 `position:fixed`、
   挤进底部状态栏，全程零报错。现在构建脚本做 CSS 自检，
   `svguest_e2e.py` 也直接断言浮条的 `getComputedStyle`。
3. ⚠ **保存是覆盖写**：svguest 会用它自己的序列化器重写整个文件（注释/属性顺序可能变），
   所以每次写回前都备份；要退回旧版用检查器里的「回滚」，或去 `outputs/_svguest_backup/` 找。
4. ⚠⚠ **`svg_size_pt()` 传 str 且带 XML 声明会静默返回 `(None, None)`**。
   lxml 拒绝解析「带 `encoding=` 声明的 unicode 字符串」，异常被 `parse_svg` 吞掉 → 返回 `None`。
   而浏览器 `svgText()` 输出的 SVG **一定**带 `<?xml version="1.0" encoding="UTF-8"?>` ——
   于是尺寸校正永远不生效、PDF 页尺寸悄悄错 4/3，**零报错**。
   判据：`svg_size_pt(svg)` → `(None, None)` 而 `svg_size_pt(svg.encode())` → 正常值。
   **凡是把浏览器来的 SVG 交给 svg_edit，一律先 `.encode()`。**
5. ⚠⚠ **PyMuPDF 把无单位 `width="920"` 当成 920 pt**（正确是 690 pt）→ px 系 SVG 转 PDF 虚胖 4/3。
   `_svg_to_pdf()` 用 `svg_size_pt()` 复核，偏差 >0.5% 时 `show_pdf_page` 重排到正确尺寸；
   e2e 直接断言 PDF 页宽 == SVG 声明宽度。

⚠ **许可**：svguest 是 GPL-3.0（强 Copyleft）。目前仅本地自用不触发义务；
**若要公开发布本平台，须先评估**（详见 `static/vendor/NOTICE.md`）。

```bash
C:/Python312/python.exe figure_studio/svguest_e2e.py   # 端到端 75 项：CSS→载入→画形→写回→备份→回滚→导出→多标签→AI 代理→自动快照
```

---

## AI 绘图（SVG Chat，需自备 key）

svguest 自带一个**对话式改图**面板：描述你要的画法或改法（也能贴参考图），
它直接改画布。检查器底部 **「AI 绘图」** 区块负责配置与入口。

### 为什么经本机服务端转发

```
GLM      浏览器 ──✗ 预检无 ACAO──> open.bigmodel.cn     ← 直连不可用
GLM      浏览器 ──✓──> 127.0.0.1:8770/api/ai/chat ──✓──> 上游 API
DeepSeek 浏览器 ──✓ 直连本来也行──> api.deepseek.com     ← 走代理是为了密钥托管
```

实测两个 API 的 CORS 行为**并不一样**（别想当然）：

| | OPTIONS 预检 | POST 响应 | 浏览器直连 |
|---|---|---|---|
| `api.deepseek.com` | 200，`ACAO` **回显任意 origin** | 带 ACAO | ✅ 可以 |
| `open.bigmodel.cn` | 200，但**没有 ACAO** | 带 ACAO | ❌ 预检就被拦 |

所以：
- **GLM 必须走代理**。上游也只给 GLM 配了 Cloudflare Worker 中继、DeepSeek 直接
  `send()` —— 正好印证这个差异。但那个 Worker 只放行 `jieqianghe.github.io` /
  `svguest.online` 等三个域名，**`127.0.0.1` 不在其中**，本平台用不了它。
- **DeepSeek 走代理是为了密钥托管**：key 只存本机磁盘，不落浏览器 localStorage，
  并且两条链路共用一套错误呈现。

⚠ **流式必须逐块透传**。上游用 `r.body` 读 SSE 做打字机效果；若服务端攒成一坨，
功能不报错，只是「等 N 秒然后全文一次性蹦出来」。服务端用
`stream_with_context` + 不设 `Content-Length`，实测 5 个事件按 0.4 s 间隔
逐块到达（`selftest.py` 里固化成了回归项）。

### 配置

| 方式 | 位置 | 说明 |
|---|---|---|
| 工作台右侧「AI 绘图」 | `figure_studio/.ai_keys.json` | 推荐。填一次即可，前端只显示掩码（`sk-****abcd`） |
| svguest 自带「Models and keys」 | 浏览器 localStorage | 上游原方式，也支持（请求同样走代理，不会 CORS） |

支持 **DeepSeek**（`deepseek-v4-pro` / `deepseek-flash`）与
**GLM**（`glm-5.3` / `glm-5.3-flash`）。

```bash
# 也可以直接用接口写
curl -X POST http://127.0.0.1:8770/api/ai/config \
  -H "Content-Type: application/json" \
  -d '{"provider":"Ds","key":"sk-..."}'
```

### 怎么接的

上游把「发请求」收敛在一个函数 `apiPost(p, url, key, body, o)` 上、**返回标准
`Response`**，调用方用 `r.body` 读 SSE 做打字机效果。桥接只做一件事：

```js
window.apiPost = (p, u, key, body, o) => fetch('/api/ai/chat', {
  method: 'POST', headers: {'Content-Type': 'application/json'},
  body: JSON.stringify({ provider: p, url: u, key: key || '', body }),
  signal: o && o.signal,
});
```

于是流式解析、重试、错误提示、`Models and keys` 里的「测试连接」全部**无感复用**。
服务端逐块透传（`stream_with_context` + `direct_passthrough`，不设 `Content-Length`），
状态码与 `Content-Type` 原样回传。

### 两个坑

1. ⚠⚠ **错误体必须是 `{"error":{"message":…}}`，不能是 `{"error":"字符串"}`**。
   上游 `apiMsg()` 解析的是 `j.error.message`，给它字符串会把整段 JSON
   原样糊到用户脸上。写成对象后它才会渲染成 `[no_key] 还没有配置…`。
2. ⚠⚠ **`apiPost` 能替换，是因为它是顶层 `async function` 声明**（全局对象上的
   可写属性）。换成 `let/const` 就是不可重新赋值的词法绑定，替换会抛 `TypeError`。
   升级上游时若改了声明方式，这条会静默失效 —— e2e 里断言了 `apiPost.__fsProxied`。

⚠ **SSRF 防护**：`/api/ai/chat` 只接受 `https` + 主机精确匹配
`api.deepseek.com` / `open.bigmodel.cn`，其余一律 403。不要把本机变成
任意内网请求的跳板。

---

## 自动保存（两层，互补）

| 层 | 存在哪 | 何时写 | 换机器/清缓存 |
|---|---|---|---|
| **上游 Auto-save** | 浏览器 IndexedDB（`sgDocs`） | 编辑后 2–10 s，默认开启，上限 10 份 / 200 MB | ❌ 丢 |
| **磁盘自动快照**（本平台加的） | `outputs/_svguest_autosave/<名>__<hash>/` | 编辑停顿 **15 s**，每源保留 12 份 | ✅ 在 |

两层都不覆盖原文件 —— **写回仍然只由 Ctrl+S / 「保存回工作台」触发**。
这条边界是有意的：自动保存负责「别丢」，手动保存负责「确认落盘」。

检查器里 **「磁盘自动快照」** 列出每份快照的时间与大小，点「恢复」写回原文件
（写回前仍会进 `_svguest_backup/` 备份一次，所以恢复本身也能撤销）。

### 怎么实现的

触发信号用上游的 `dv`（document version，每次编辑自增）而不是监听 DOM：
比 `MutationObserver` 省，也比定时比对全文准。

```js
setInterval(() => {            // 2 s 轮询
  if (dv === lastDv) return;   // 没变化
  lastDv = dv;
  if (!dirty || !itemOf(TI)) return;   // 没改动 / 非工作台标签 → 不存
  clearTimeout(t); t = setTimeout(doSnap, 15000);   // 停顿 15 s 才写
}, 2000);
```

序列化用上游的 `svgRaw()`（注释原话：*the document as text without ending an edit
(autosave)*），不打断正在进行的文本编辑。

⚠ 空画布不存（复用写回那道 `collapsed()` 判定），否则会存一堆空白。
⚠ `#…&autoMs=1500` 可覆盖节流间隔 —— 给 e2e 真实验证节流链路用的，正常使用不用管。

---

## 三种典型用法

**A. 纯贴图自由组合**（最常用）
1. 左侧填面板目录 → 扫描 → 按顺序单击加入
2. 工具栏「横排一行」/「网格 2 列」先摆个大概
3. 拖拽微调、属性面板改精确坐标、加题注
4. 「裁到内容」收掉空白 → 「导出 PNG + PDF」

**B. 改版式（面板图不重出）**
打开已保存的布局 → 换顺序 / 调行列 / 改间距 / 改字母 → 重新导出。
**全程不碰任何绘图脚本。**

**C. 改内容（面板图必须重出）**
面板数据变了 → 用原脚本重出该面板（覆盖同名文件）→
回平台点「导出」，平台会自动读到新图。
⚠ 项目纪律：改源图必须**先写 `08_report/` 再重出主图**，顺序反了主图会嵌旧面板。

**D. 只改配色（面板图不重出）**
选中面板 → 「取主要颜色」→ 把偏色改成项目色 → 导出。
**R 脚本一行不动。** 这类改动同样只落在 layout JSON 的 `adjust` 里，随时可回退。

---

## 导出

* **PNG**：`canvas(in) × dpi` 像素，`int()` 截断（非整数英寸可能差 1 px）；
  自动剥掉 alpha 通道（RGB），避免投稿系统拒收。
* **PDF**：矢量，`pdf.fonttype = 42`（嵌入 TrueType，可在 Illustrator / Inkscape 里编辑文字）。
* ⚠ **不使用 `bbox_inches="tight"`** —— 画布尺寸严格等于你设定的尺寸，不被内容撑大/裁小。
* 导出前会逐块比对**面板比 vs 源图比**，偏差 >1% 会在弹窗里列出（这是防形变的最后一道闸）。
* 还会检查**静默失效**（matplotlib 对越界内容一律零报错）：
  面板块超出画布 / 题注超出画布下缘 / 题注右溢出 / 题注折行后仍超宽。命中即列入告警。

---

## 布局 JSON

存在 `figure_studio/layouts/`，可直接手改：

```json
{
  "canvas": { "width_in": 10.0, "height_in": 6.19, "dpi": 300 },
  "style":  { "frame": true, "letter_fs": 11, "letter_dy": 1.02,
              "letter_color": "#111111", "show_letters": true,
              "caption_fs": 7, "caption_gap_in": 0.06, "caption_color": "#222222",
              "margin_in": 0.08, "pack_margin_in": 0.08, "gap_in": 0.12 },
  "title":  { "text": "", "subtitle": "", "fs": 13, "y_in": 0.30,
              "color": "#000000", "sub_color": "#444444" },
  "items": [
    { "id": "i1", "src": "D:/.../panel_A.png",
      "x_in": 0.1, "y_in": 0.9, "w_in": 2.4, "h_in": 2.1,
      "letter": "A", "caption": "…", "z": 0,
      "adjust": { "grayscale": true,
                  "remap": [ { "from": "#0172b2", "to": "#0072B2", "tol": 8 } ] } }
  ]
}
```

* 单位一律**英寸**，原点**左上角**（与前端一致；后端渲染时转 matplotlib 左下原点）。
* `x/y/w/h` 与源图宽高比不符 → 会形变，导出时告警。
  ⚠ 带 `trim` 时比例以**调色后**的图为准（前后端都取处理后的尺寸，不会各算一套）。
* `adjust` 可省；只写要改的键，其余取默认。`caption_max_w_in` 可覆盖题注折行宽度。

命令行也能直接渲染（批量/无界面）：

```bash
C:/Python312/python.exe figure_studio/composer.py layouts/xxx.json -o outputs --stem Fig1
```

---

## HTTP 接口

| 方法 | 路径 | 作用 |
|---|---|---|
| GET | `/api/panels?root=&q=&limit=` | 扫描面板库 |
| GET | `/api/thumb` `/api/img?src=&max=&adjust=` | 缩略图 / 画布预览（`adjust` 为 URL 编码 JSON；SVG 的 `svg` 子树在这里先做矢量编辑再栅格化） |
| POST | `/api/probe` | `{src, dpi?, adjust?}` → 尺寸与宽高比（调色后重算用） |
| POST | `/api/colors` | `{src, n?}` → 主要颜色 + 色桶数 + 源图模式（SVG 额外给语义色/图例/文本/字体） |
| POST | `/api/svg/inspect` | `{src}` → SVG 尺寸/语义色/文本清单/组结构/图例候选 |
| POST | `/api/svg/source` | `{src, svg?}` → 施加矢量编辑后的 SVG 字节（浏览器可直接存） |
| GET | `/api/svg/raw?src=` | 原样输出 SVG 字节（⚠ MIME 必须是 `image/svg`，见「svguest 精修」节） |
| POST | `/api/svg/writeback` | `{src, svg}` → 写回原路径（先备份到 `outputs/_svguest_backup/`） |
| POST | `/api/svg/writeback_batch` | `{items:[{src,svg}]}` → 批量写回，逐项独立成败（一个失败不影响其余） |
| POST | `/api/svg/export` | 位图/SVG：`{src, kind:'png'\|'jpg'\|'svg', dpi, data(base64)}` → 落到**原图同目录**（`<名>_edit_300dpi.png` / `_edit.svg`）；矢量 PDF：`{src, kind:'pdf', svg(文本)}` → `<名>_edit.pdf`。重名自动加 `_2` |
| GET | `/api/svg/history?src=` | 该 SVG 的备份列表（新→旧，最多 50 条） |
| POST | `/api/svg/restore` | `{src, backup}` → 回滚到某个备份（⚠ 回滚前先把当前版本再备份一次） |
| POST | `/api/svg/session` | `{items:[{src,name}]}` → `{token}`，把「要精修哪些面板」存服务端，hash 只带短 token（避免中文路径反复编码 / URL 超长） |
| GET | `/api/svg/session/<token>` | 取该会话的面板清单 + 各接口地址（30 分钟过期） |
| POST | `/api/svg/autosave` | `{src, svg}` → 存一份磁盘快照到 `outputs/_svguest_autosave/<名>__<hash>/`（**不覆盖原文件**），每源保留 12 份 |
| GET | `/api/svg/autosave?src=` | 该源的快照列表（新的在前） |
| POST | `/api/svg/autosave_restore` | `{src, file}` → 把快照写回原文件（写回前仍进 `_svguest_backup/` 备份） |
| POST | `/api/ai/chat` | `{provider, url, key?, body}` → 流式转发到 DeepSeek / GLM（白名单 + https，防 SSRF） |
| GET/POST | `/api/ai/config` | 读（只回掩码）/ 写 AI 密钥到本机 `.ai_keys.json` |
| GET | `/api/palette` | 项目调色板（从 `sci_palette.py` 动态加载） |
| POST | `/api/render` | `{layout, stem, formats, dpi}` → 渲染 PNG/PDF（SVG 面板在 PDF 里矢量嵌入） |
| GET/POST/DELETE | `/api/layouts[/<name>]` | 布局列表 / 读取 / 保存 / 删除 |

---

## 自检

```bash
C:/Python312/python.exe figure_studio/selftest.py   # 后端：扫描/排版/渲染/尺寸/往返/调色/SVG 矢量编辑
C:/Python312/python.exe figure_studio/uitest.py     # 前端：真浏览器点一遍（需先启动服务）
```

`uitest.py` 判据不是「页面能打开」，而是面板库真的出卡片、点卡片画布真的出块、
「横排」真的改了几何且同行等高、勾灰度预览 URL 真的带 `adjust`、
「取主要颜色」真的生成替换行、SVG 面板的语义色/图例/字号改动真的进入 URL 且重置后消失、
导出真的拿到产物 URL，且**全程零 JS 错误 / 零 4xx**。

`selftest.py` 的调色判据是**裁出该块的像素区**再验 R==G==B，并确认**其他块仍是彩色**
（证明逐块生效，不是整图被改）。SVG 段的判据更硬：每个算子（改色/图例隐藏/字号/
字体族/文本替换/图例缩放/图例移动）都必须产生**差异像素**（>2000 px），
因为「改了 XML 但栅格 0 差异」正是渲染器缺陷的静默形态；另有物理尺寸三处口径
一致、4 种 SVG 预览端点 200、混排 PDF 的 SVG 面板是矢量（1361 条路径 + 可提取文本）。

---

## 目录

```
figure_studio/
├── app.py           Flask 服务与 API
├── composer.py      渲染核心（图像级调色 + layout JSON → matplotlib → PNG/PDF，SVG 矢量嵌入）
├── svg_edit.py      SVG 矢量编辑（解析/探测/DOM 改色/图例/文字/字号/线宽/隐藏 + fitz 栅格化）
├── panels.py        面板库扫描 + 缩略图/预览图缓存
├── demo_types.py    异构图型实测（8 种真实成品图拼通栏图）
├── demo_adjust.py   调色算子对照图（同版式，只改 adjust）
├── demo_svg.py      SVG 矢量编辑演示（4 SVG + 1 PNG 拼通栏图，PDF 矢量嵌入）
├── svguest_e2e.py   svguest 深度整合端到端测试（75 项：CSS→载入→画形→写回→备份→回滚→导出→多标签→AI 代理→自动快照）
├── selftest.py      后端自检
├── uitest.py        前端真浏览器测试
├── shot_flex.py     交付截图：多选 / 三档布局 / 工具栏
├── tools/
│   └── build_vendor.py  合成 vendor_src → static/vendor/svguest.html（含 CSS 自检 + JS 语法校验）
├── vendor_src/      svguest 上游原件 + 桥接源码（改这里，再跑 build_vendor.py）
│   ├── app.html     上游 v1.0.1 原样（未改动，sha256 记在 NOTICE.md）
│   ├── bridge.css · bridge.js   桥接（浮条 / 会话载入 / 写回 / 导出直通 / AI 代理 / 自动快照）
│   ├── NOTICE.md · LICENSE · README.upstream.md
├── static/
│   ├── index.html · app.js · style.css
│   └── vendor/      构建产物：svguest.html（勿手改）+ LICENSE + NOTICE.md
├── layouts/         保存的布局 JSON
├── .ai_keys.json    AI 密钥（本机，勿提交；前端只显示掩码）
├── outputs/         导出产物（`_svguest_backup/` = 写回/回滚前旧版；`_svguest_export/` = 无源文件时导出落点；`_svguest_autosave/` = 编辑过程磁盘快照）
└── cache/           缩略图 / 预览图 / 调色派生预览 / 元数据缓存
```

## 与项目现有管线的关系

| | 位置 | 形态 |
|---|---|---|
| 本平台 | `figure_studio/` | **通用**贴图合成器，任意目录的面板图都能组合 |
| phylo 主图 | `virome_phylo_pipeline/utils/core_results_fig.py` | 固定 `ROWS`/`PANELS` 的专用合成器 |
| discovery 主图 | `reports/virome_analysis_pipeline/make_discovery_sci_figs.py` | 绘图函数 + 贴图混合，`_v23_geometry()` 反解 |

本平台把前两者的**反解逻辑 + 面板字母偏移 + 零形变纪律**抽了出来，
并补上交互层。现有管线脚本不受影响，可继续用。
