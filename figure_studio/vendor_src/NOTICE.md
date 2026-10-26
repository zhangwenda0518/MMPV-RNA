# vendor/ · 第三方组件与改动说明

## svguest

- 上游：https://github.com/JieqiangHe/svguest
- 版本：**v1.0.1**（`vendor_src/app.html`，538624 字节）
- sha256：`54dae0ac49d4e3d64aca1fdb80eaa1e06c22aac45c8e881bc6a40a220c07839e`
- 许可：**GPL-3.0-or-later**（全文见同目录 `LICENSE`）
- 版权：© 2026 Jieqiang He

### 我们改了什么

**上游 `app.html` 一行未改。** `static/vendor/svguest.html` 由
`tools/build_vendor.py` 自动生成，等于：

```
vendor_src/app.html  +  vendor_src/bridge.css  +  vendor_src/bridge.js
```

桥接块插在 `</body>` 之前，用 `<!-- ▼▼▼ Figure Studio 桥接 …` 标记包裹，
升级上游时只换 `vendor_src/app.html` 再重跑构建即可（构建脚本会先剥掉旧块）。

桥接做的事：

| 动作 | 实现 |
|---|---|
| 从工作台载入面板 | `GET /api/svg/session/<token>` 拿文件清单，逐个 `readFile()`；`readFile` 内部会自动开新标签，因此多面板 = 多标签 |
| 保存回工作台 | **接管 `A.svg`**（Ctrl+S / File▸Export SVG / Ctrl+K 三入口）→ `svgText()` 或后台标签序列化快照 → `POST /api/svg/writeback(_batch)` |
| 导出 PNG / JPG | **常驻 `dl` 分流器**按 MIME 识别（`image/png`、`image/jpeg`）→ `POST /api/svg/export`，不经过浏览器下载 |
| 导出矢量 PDF | **接管 `A.pdf`**（Ctrl+P / File▸Export PDF）→ `POST /api/svg/export {kind:'pdf'}`，服务端 PyMuPDF 生成，不弹打印对话框 |
| 导出 SVG 副本 | 直接 `svgText()` → `POST /api/svg/export` |
| 查隐藏信息 | 调用上游自带的水印检查器 `A.watermark()` |
| **AI 绘图（SVG Chat）** | **替换 `apiPost(p,url,key,body,o)`** → `POST /api/ai/chat` 服务端转发。上游的流式解析、重试、错误提示、`Models and keys` 的测试连接全部无感复用 |
| **磁盘自动快照** | 轮询上游的 `dv`（document version）+ `dirty` → 停顿 15 s 调 `svgRaw()` → `POST /api/svg/autosave`（**不覆盖原文件**） |
| 主题跟随 | `setTheme()` |

⚠ 上游自带的 **Auto-save**（浏览器 IndexedDB `sgDocs`，默认开启）**保持原样**，
桥接不接管也不关闭 —— 它和上面的磁盘快照是互补的两层：
前者防浏览器崩溃，后者防换机器/清缓存。

### 为什么能「零改动」

上游主脚本是**经典脚本**（既非 `type=module` 也非 IIFE），所以顶层
`function` 声明会挂到 `window` 上、顶层 `const/let` 落在共享的全局词法环境里。
另一个 `<script>` 因此可以直接调用 `readFile` / `svgText` / `A` / `TABS` /
`setTheme` 等，无需上游暴露任何 API。

`A` 是 `const A = {...}`（动作表）—— `const` 只锁绑定、不锁属性，
所以桥接可以安全地 `A.svg = …` / `A.pdf = …` 来换掉单个动作，
上游其余 200+ 个动作不受影响。

### ⚠ 踩过的坑（别再犯）

1. `<style>` 里**只能有 CSS 注释**。曾经在 `<style>` 开头混进一行 Markdown 标题
   （`# Figure Studio ⇄ svguest 桥接样式…`），浏览器把 `#` 当成 ID 选择器，紧接着的
   整条 `#fsBar{…}` 规则被当作无效选择器**静默丢弃** —— 浮条失去 `position:fixed`、
   挤进底部状态栏，全程零报错。现在 `tools/build_vendor.py` 会做 CSS 自检，
   `svguest_e2e.py` 也会断言浮条的 computed style。
2. **别只在 `keydown` 上拦 Ctrl+S**。菜单点击走 `A[b.dataset.a]()`、快捷键走
   `A.svg`/`A.pdf`、命令面板也走 `A` —— 只有改 `A` 表才能一次覆盖三个入口；
   只拦键盘会漏掉 File 菜单，那边会去调 `showSaveFilePicker` 弹文件选择器。
3. **`dl()` 不是导出专用的**：FigKit 模板导出（`application/json`）和
   「Kept in this browser」（`application/zip`）也走它。分流要按 MIME 判断，
   按调用点分流会把这两类也吞掉，用户拿不到文件。
4. **`apiPost` 能替换，靠的是它是顶层 `async function` 声明**（全局对象上的可写
   属性）。若上游哪天改成 `let`/`const`，替换会抛 `TypeError` 而桥接会静默失效 ——
   e2e 里断言了 `apiPost.__fsProxied` 来兜这条。
5. **AI 错误体必须是 `{"error":{"message":…}}`**。上游 `apiMsg()` 读的是
   `j.error.message`；给它 `{"error":"字符串"}` 会把整段 JSON 原样显示给用户。
6. **内嵌到 iframe 之后，键盘焦点必须显式推进去**。点的是工作台上的按钮，
   焦点留在平台，`Ctrl+S` 和工具快捷键全部失效，且**零报错**。在 iframe 的
   `load` 里 `contentWindow.focus()` 兜住；e2e 断言平台侧 `activeElement.id ==
   'sgFrame'`。
7. **别抢 `Escape`**。上游把每一次 `Escape` 都 `preventDefault`（取消工具 /
   取消选择 / 退出文本编辑），冒泡阶段收得到但 `defaultPrevented` 恒为 true。
   按它关浮层会导致「想取消选择 → 编辑器没了」。退出只走顶栏按钮；平台侧
   Esc 仅当焦点不在 iframe 里时才生效。

### 许可提示

GPL-3.0 是**强 copyleft**。本仓库内部自用不触发分发义务；若要把
`figure_studio` 对外发布（尤其是以二进制/打包形式），需要评估是否构成
「分发衍生作品」，通常意味着整体也需以 GPL-3.0 兼容方式提供源码。
