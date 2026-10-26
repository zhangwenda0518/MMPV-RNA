# 图件四格式导出规范 / Figure Export Convention

> 2026-10-20 起，仓库内**所有**出图脚本默认同时产出四种文件，主名相同、仅扩展名不同。
> 这一口径由 `mmpv_common/figure_export.py`（Python）与 `mmpv_common/figure_export.R`（R）
> 统一收口，脚本本身不需要重复实现。

## 1. 默认产物

| 文件 | 用途 | 关键参数 |
|------|------|----------|
| `<stem>.svg` | 矢量，文本可编辑（Illustrator / Inkscape 改字） | `svg.fonttype='none'` |
| `<stem>.png` | 报告内嵌、快速预览 | 300 dpi、白底 |
| `<stem>.pdf` | 出版投稿、矢量排版 | 沿用调用方原有 dpi/参数 |
| `_hires_1200dpi/<stem>.png` | 印刷、审稿人索图 | 1200 dpi、白底 |

**为什么高清 PNG 放独立子目录**：报告层（`utils/report_builder.collect_plots`）
按 `(目录, 文件名主干)` 去重并只保留优先级最高的一种格式。若把 1200 dpi 副本
写成 `<stem>_1200dpi.png` 放在同目录，报告会把它当成**第二张图**重复内嵌。
放进 `_hires_1200dpi/` 后，报告的 glob 不会命中，版面不变。

## 2. 三种使用方式

### 2.1 零改动（既有脚本，已自动生效）

49 个 Python 出图脚本与 12 个 R 脚本已在文件头注入引导块，导入时安装
`Figure.savefig` 钩子（Python）或 `ggsave` 钩子（R）。原有的

```python
fig.savefig("out/fig1.png", dpi=300, bbox_inches="tight")
fig.savefig("out/fig1.pdf", dpi=300, bbox_inches="tight")
```

会自动补齐其余格式；第二次调用发现兄弟文件已存在则跳过，**不会重复渲染**。

### 2.2 显式调用（推荐新代码）

```python
from mmpv_common.figure_export import export_figure

export_figure(fig, "out/fig1.pdf", bbox_inches="tight")
# → {'svg': ..., 'png': ..., 'pdf': ..., 'hires': ...}
```

```r
source("mmpv_common/figure_export.R")

save_figure(p, "out/fig1.pdf", width = 7, height = 4)
```

### 2.3 直写设备的老封装

`pdf()` / `png()` 直接开设备的老代码不走 ggsave 钩子，用 `draw_figure_formats()`
包一层（`viral_maftools.R`、`panvirome_suite.R` 的 `safe_plot()` 已按此改造）：

```r
draw <- function() { print(p) }        # 或 base 绘图语句
draw_figure_formats(draw, "out/fig1.pdf", width = 8, height = 5)
```

## 3. 环境变量

| 变量 | 默认 | 说明 |
|------|------|------|
| `MMPV_FIG_EXPORT` | `1` | `0` 关闭**自动补齐**（钩子）；显式 `export_figure` 不受影响 |
| `MMPV_FIG_FORMATS` | `svg,png,pdf,hires` | 裁剪格式，如 `png,pdf` |
| `MMPV_FIG_PNG_DPI` | `300` | 普通 PNG 的 dpi |
| `MMPV_FIG_HIRES_DPI` | `1200` | 高清 PNG 的 dpi |
| `MMPV_FIG_HIRES_DIR` | `_hires_1200dpi` | 高清副本子目录名 |

`MMPV_FIG_EXPORT=0` 与 `formats=()` 的区别：前者只关自动补齐，后者让显式调用
也什么都不写。测试环境用前者（见仓库根 `conftest.py`），避免钩子污染
`pdf_to_png` 这类「先造 PDF 再断言兜底渲染出 PNG」的用例。

## 4. 安全边界

- **多页 PDF 不补格式**：`PdfPages` 与文件对象（`BytesIO`）目标原样透传。
- **单格式失败不中断管线**：某一种格式写失败只告警，其余格式照常落盘。
- **超大画布自动降级**：预估像素超过 4 亿时，高清 PNG 从 1200 dpi 降到 600 dpi
  并告警，避免 OOM。
- **中间文件豁免**（R）：路径落在 `tempdir()` 下或以 `tmp`/`_tmp` 开头的文件
  只写调用方指定的那一种格式 —— inset / 预览图不刷四份。

## 5. 未纳入范围

| 对象 | 原因 |
|------|------|
| `figure_studio/composer.py` | 交互式渲染核心，`formats` 由前端显式传入；每次预览都写四份会拖慢 UI |
| `archive/`、`reports/_archive/`、`_scratch_tests/` | 归档与一次性脚本，不在生产路径 |
| `biosoft/pypopart/` | 第三方工具源码副本，随上游维护 |
| `_refs/` | 外部参考仓库副本 |
| 第三方工具直写的 PDF（TreeTime / SpreaD3 / BEAST） | 改不了其内部；仍由 `utils/pdf_to_png.py` 兜底渲染 PNG |

## 6. 验证

```bash
# Python 侧
python -m pytest mmpv_common/tests/test_figure_export.py virome_phylo_pipeline/tests/test_plot_io.py -q

# 端到端：任一脚本跑完检查四格式
python doc/figures/make_module_overview.py
ls doc/figures/MMPV_module_overview.{svg,png,pdf} doc/figures/_hires_1200dpi/
```

实测参考值：`MMPV_module_overview`（190×128 mm）产出
PNG 2244×1511、高清 PNG 8976×6047、SVG/PDF 为矢量。
