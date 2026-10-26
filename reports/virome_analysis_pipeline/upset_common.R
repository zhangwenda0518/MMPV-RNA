# ── 图件四格式导出 (SVG / PNG / PDF / PNG@1200dpi) — 默认开启 ──────────
# 2026-10-20: 所有出图脚本统一口径; MMPV_FIG_EXPORT=0 可整体关闭。
# 钩子只拦截 ggsave(); pdf()/png() 直写的老代码请改用 save_figure()。
local({
  .find_repo <- function(start) {
    d <- normalizePath(start, mustWork = FALSE)
    repeat {
      if (dir.exists(file.path(d, "mmpv_common"))) return(d)
      p <- dirname(d)
      if (identical(p, d)) return(NA_character_)
      d <- p
    }
  }
  .here <- tryCatch({
    a <- commandArgs(FALSE)
    f <- sub("^--file=", "", a[grep("^--file=", a)])
    if (length(f)) dirname(normalizePath(f[1], mustWork = FALSE)) else getwd()
  }, error = function(e) getwd())
  .root <- .find_repo(.here)
  if (is.na(.root)) .root <- .find_repo(getwd())
  if (!is.na(.root)) {
    source(file.path(.root, "mmpv_common", "figure_export.R"), local = FALSE)
    # upset_common.R 是被 source() 的公共库 —— 只加载函数, 安装由调用方脚本负责
  } else {
    message("[figure_export] 未找到 mmpv_common/, 跳过四格式导出")
  }
})
# ──────────────────────────────────────────────────────────────────────
# upset_common.R — 枸杞 discovery 侧 UpSet 系列：共用配色与 ggplot 主题
# 配色与 virome_discovery_pipeline/utils/sci_palette.py 的 OKABE_ITO 同源（顺序固定，勿改）。

OKABE_ITO <- c(
  "#E69F00", # orange
  "#56B4E9", # sky blue
  "#009E73", # bluish green
  "#F0E442", # yellow
  "#0072B2", # blue
  "#D55E00", # vermillion
  "#CC79A7", # reddish purple
  "#000000", # black
  "#999999"  # grey
)

# 7 个分类工具的固定配色（供三张图统一，按 OKABE_ITO 顺序取前 7）
TOOLS7 <- c("diamond_lca", "mmseqs", "CAT", "metabuli", "genomad", "VITAP", "ACVirus")
TOOL_COLORS <- setNames(OKABE_ITO[1:7], TOOLS7)

# degree（重叠工具数）1..7 的连续配色：浅→深取自同色系，末档用灰示意"全票共识"
DEGREE_COLORS <- c(
  "1" = "#F0E442", "2" = "#B8DF8E", "3" = "#7CCD7C", "4" = "#56B4E9",
  "5" = "#449B99", "6" = "#0072B2", "7" = "#D55E00"
)

# 宿主预测器配色（图③）
HOSTPRED_COLORS <- c(
  "ICTV" = "#0072B2",
  "RVH" = "#D55E00",
  "PB2" = "#009E73",
  "Rule" = "#E69F00"
)

# ══════════════════════════════════════════════════════════════════════
# v22.2（2026-10-14）：**类别色 vs 数量色分离** —— 三张 UpSet 统一配色
# ══════════════════════════════════════════════════════════════════════
# 用户诉求：「三个 upset 的颜色可以调整下嘛？尤其是左侧的条形图」+「也一并协调顶部交集柱」
# 经确认采纳：① 左侧 set-size 条 = **身份色**（工具/预测器固定色）；
#             ② 顶部交集柱 = **三张统一同一套数量色阶**。
#
# ⚠⚠⚠ **为什么必须分离**（本次实测的撞色清单，勿凭直觉合并）：
#   Q 左侧若用 TOOL_COLORS(7 工具)，顶部原用 DEGREE_COLORS(共识度 1..7) →
#     **4 个完全同色**（#56B4E9/#F0E442/#0072B2/#D55E00）+ 2 个近色 → 图例串味；
#   P 左侧若用 HOSTPRED_COLORS(4 预测器)，顶部按命中数 → **2 个完全同色**
#     （#0072B2/#D55E00）。
#   ⇒ 正解 = 「**类别用彩色、数量用单色系渐变**」（UpSet 经典规范）。
#
# ⚠⚠ 数量色阶必须**跨图同值同色**：三张图分别出现 1..7 / 1..4 / 4 档，
#   若各自 `viridis(n)` 插值，同一个「2」在三张图里颜色会不同（实测常见误用）。
#   ⇒ 做法 = **固定取 viridis 的 7 档**，`QUANT_COLORS[k]` 直接按**档位**索引，
#      无论图里只有几档，档位 k 的颜色恒定。
# ⚠ viridis 的 1..7 顺序是「深紫 → 黄」，语义 = **数值越大越亮**。
#   用户既有口径（`DEGREE_COLORS`）是「1=浅黄 … 7=橙红」= 数值越大越深 →
#   ⚠ 两种方向**相反**。本次按用户「三张统一」的要求重定为 viridis 方向
#   （**大 = 亮**），并在各图题注里写清 `colour = number of agreeing tools`，
#   避免读者按旧习惯反读。改方向属**用户可见变更**，已在 README 记录。
QUANT_COLORS <- stats::setNames(
  viridisLite::viridis(7),
  as.character(1:7)
)
# 便捷取色：给一组整数档位，返回对应颜色（顺序与输入一致）
quant_colors <- function(levels) {
  lv <- as.character(levels)
  miss <- setdiff(lv, names(QUANT_COLORS))
  if (length(miss))
    stop("quant_colors(): 档位超出 1..7 范围 -> ", paste(miss, collapse = ", "))
  unname(QUANT_COLORS[lv])
}
# 单色系（蓝）数量渐变 —— 备选方案 B，保留供切换（`--quant-scheme blues`）。
QUANT_BLUES <- c("1" = "#C6DBEF", "2" = "#9ECAE1", "3" = "#6BAED6", "4" = "#4292C6",
                 "5" = "#2171B5", "6" = "#08519C", "7" = "#08306B")

# ⚠⚠⚠ **set-size 面板上色只能用「颜色向量」**（2026-10-14 实测三种写法全废）：
#   `upset_set_size(mapping = aes(fill = x))` → `object 'x' not found`
#     （该 mapping 被当**图层级之外**的表达式求值，见 ComplexUpset 1.3.3 源码）；
#   `geom_bar(aes(fill = set))` → `object 'set' not found`；
#   `geom_bar(aes(fill = x))` + `scale_fill_manual` → `Cannot add <ggproto> objects together`。
#   ✅ **正解 = `geom_bar(fill = <颜色向量>)`**（纯常量向量，不经 aes）。
#   ⚠⚠ 向量按**传入 `sets` 的顺序**逐条对应，**与 `sort_sets` 造成的渲染顺序无关**
#     （实测：默认 `sort_sets="descending"` 下，`sets[i]` 的条仍拿到 `vec[i]`）。
#     ⇒ 构造向量时**只按 sets 顺序排**，不要去跟渲染后的上下位置对齐。
#   实测证据（真实数据，逐个上洋红）：Q 传入
#     `c(diamond_lca, mmseqs, CAT, metabuli, genomad, VITAP, ACVirus)` →
#     每根条的颜色恒为同序号工具色，条长与该工具真实 set size 吻合。
set_size_colors <- function(sets, pal) {
  cols <- unname(pal[sets])
  if (anyNA(cols))
    stop("set_size_colors(): 集合缺色 -> ",
         paste(sets[is.na(cols)], collapse = ", "))
  cols
}

# ⚠⚠⚠ **颜色向量长度必须等于 ComplexUpset 内部实际的集合数**（2026-10-14 实测踩坑）：
#   `min_size` 会把「在**所有保留交集**里都不出现」的集合**静默裁掉**（它内部只保留
#   出现在 `filter_intersections` 后数据里的 set）。此时若仍按**传入的完整 sets**
#   给颜色向量 → `geom_bar` 报
#     `Aesthetics must be either length 1 or the same as the data (6). Fix: fill`
#   （实测：Species 的 `genomad` 仅 13 条，而 `min_size=68` → 被裁 → 7 色 vs 6 行）。
#   ⚠ 该错误**只在该 rank 触发**：Family 的 min_size=141 下 7 个工具全在保留交集里
#     → 不报错。**别以为 Family 跑通就万事大吉**（必须两个 rank 都跑）。
#   ✅ 正解 = 先按「该集合是否出现在保留交集中」过滤，再取色。
#      判定口径与 ComplexUpset 一致：集合 t 保留 ⟺ 保留交集的签名里该位出现过 "1"。
#   ⚠ 过滤后集合的**渲染顺序不变**（仍按传入 sets 的先后），故颜色对位仍然正确。
set_size_colors_filtered <- function(sets, pal, kept_signatures) {
  n <- nchar(kept_signatures[1])
  if (n != length(sets))
    stop("set_size_colors_filtered(): 签名长度 ", n, " != 集合数 ", length(sets))
  pos <- seq_along(sets)
  present <- vapply(pos, function(i)
    any(substr(kept_signatures, i, i) == "1"), logical(1))
  cols <- unname(pal[sets[present]])
  if (anyNA(cols))
    stop("set_size_colors_filtered(): 集合缺色 -> ",
         paste(sets[present][is.na(cols)], collapse = ", "))
  cols
}

theme_upset <- function(base_size = 8) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      axis.title = ggplot2::element_text(size = base_size),
      axis.text = ggplot2::element_text(size = base_size - 0.5, color = "black"),
      plot.title = ggplot2::element_text(size = base_size + 1.5, face = "bold", hjust = 0),
      plot.subtitle = ggplot2::element_text(size = base_size - 0.5, color = "#444444", hjust = 0),
      legend.title = ggplot2::element_text(size = base_size - 0.5),
      legend.text = ggplot2::element_text(size = base_size - 1),
      plot.margin = ggplot2::margin(6, 8, 6, 8)
    )
}

# 保存 PNG + PDF（源画布 = 排版宽度，避免字号被压小）
save_upset <- function(p, outdir, name, width, height) {
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  ggplot2::ggsave(file.path(outdir, paste0(name, ".png")), p,
    width = width, height = height, dpi = 300, bg = "white")
  ggplot2::ggsave(file.path(outdir, paste0(name, ".pdf")), p,
    width = width, height = height, bg = "white")
  invisible(file.path(outdir, paste0(name, ".png")))
}
