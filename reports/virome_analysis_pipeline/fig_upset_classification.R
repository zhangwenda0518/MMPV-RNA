#!/usr/bin/env Rscript
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
    install_ggsave_hook()
  } else {
    message("[figure_export] 未找到 mmpv_common/, 跳过四格式导出")
  }
})
# ──────────────────────────────────────────────────────────────────────
# fig_upset_classification.R — 图② 分类证据（工具赞成票）UpSet
#
# 数据：10_Reports/final_integrated_classification.tsv 的 `<Rank>_agree` 列
#   口径：agree 串 "n_agree/n_vote: toolA,toolB" → 括号内工具对该 contig
#         在该 rank 上"投了赞成票"（taxon 一致）。
#   本图取 **Species rank**（最下游、最严）：17,282 条有物种级赞成票的 vOTU
#   × 7 个分类工具 = "哪些工具在物种级互相印证"。
#
# 图形：ComplexUpset —— 交集大小柱（按 degree 着色）+ 成员点阵。
#
# 来源：数据 = 项目 05_Taxonomy 管线产物；绘图代码为本次独立编写（PlotCase 简洁风格）。
# 说明：PlotCase MCP（127.0.0.1:17890）连接检查未通过，故未经案例库改编。

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(ComplexUpset)
  library(grid)
  library(png)
  # ⭐ 2026-10-14 第 8 次修订：右上角 Venn inset 需要 ggrepel（防标签重叠）
  #   与 ggVennDiagram（与 ③ HostPrediction BPanel 的 Venn 同族，风格一致）。
  library(ggrepel)
  library(ggVennDiagram)
})
script_dir <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) dirname(f) else "."
}
source(file.path(script_dir(), "upset_common.R"))

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(k, d) {
  i <- which(args == k)
  if (length(i)) args[i + 1] else d
}
# ⚠ 阶元参数化（2026-10-13）：`--rank Family` 出科级版，输出名带 `_Family` 后缀。
#   默认仍是 Species（与既有 Fig_GojiVirome_UpsetClassification 的零回归对拍口径一致）。
rank <- get_arg("--rank", "Species")
mtx <- get_arg("--matrix", file.path(script_dir(), "upset_data",
                                     sprintf("agree_matrix_%s.tsv", rank)))
outdir <- get_arg("--outdir", "reports/_discovery_sci_20261012/_lt/10_Reports/S10_SCI_Discovery")
# ⚠⚠ `--top` **不能随便设**：自绘柱顶数字要求「前 N 个交集大小互不相同」
#   （有并列时 ComplexUpset 的先后与本表顺序可能不一致 → 数字整体错位，**零报错**）。
#   实测无并列的安全值：**Species = 16**（前 17 起出现并列 68）、**Family = 18**。
#   ⚠ 历史：Species 曾用 18（旧版靠 ComplexUpset 内置 `counts=TRUE` 标注，它用自己的
#     内部排序 → 不会错位）；改为自绘后必须降到 16。脚本末尾有自检，并列即 `stop()`。
.DEFAULT_TOP <- list(Species = 16L, Family = 18L)
top_n <- as.integer(get_arg("--top",
  as.character(if (!is.null(.DEFAULT_TOP[[rank]])) .DEFAULT_TOP[[rank]] else 18L)))
suffix <- get_arg("--suffix", if (rank == "Species") "" else paste0("_", rank))
# ⚠ 布尔开关：共享 dest 时**先注册者的 default 生效** → 用「只查存在性」的写法更稳。
want_timing <- !("--no-timing" %in% args)
# ⭐ v22.2（2026-10-14）配色开关：左侧 set-size 条改用 TOOL_COLORS（7 工具固定色）。
#   ⚠ 默认**开**（用户 2026-10-14 确认「方案A 身份色」）；`--plain-setsize` 回到旧深灰，
#     供零回归对拍（旧行为可一键复现，是本次回归验证的手段之一）。
want_tool_color <- !("--plain-setsize" %in% args)
# 顶部交集柱的数量色阶：viridis（默认，三张统一）| blues（备选单色系）| degree（旧 DEGREE_COLORS）
quant_scheme <- get_arg("--quant-scheme", "viridis")

m <- fread(mtx, sep = "\t", header = TRUE)
tools <- TOOLS7[TOOLS7 %in% names(m)]
m[, (tools) := lapply(.SD, as.logical), .SDcols = tools]
m[, count := rowSums(.SD), .SDcols = tools]

cat(sprintf("总计 %s 条 vOTU（%s rank 有赞成票）x %d 个分类工具\n",
  format(nrow(m), big.mark = ","), rank, length(tools)))
print(table(m$count))

sig <- m[, do.call(paste0, lapply(.SD, as.integer)), .SDcols = tools]
sizes <- sort(table(sig), decreasing = TRUE)
n_all <- length(sizes)
n_show <- min(top_n, n_all)
min_size <- as.integer(sizes[n_show])
cat(sprintf("非空交集共 %d 个；显示前 %d 个（size >= %d）\n", n_all, n_show, min_size))

# ── 柱顶标注用：前 n_show 个交集（按大小降序，与 ComplexUpset 的
#    `sort_intersections_by="cardinality"` 同序 → x = 1..n_show 逐位对应）──
# ⚠ 自检：前 n_show 个交集大小必须**互不相同** —— 有并列时 ComplexUpset 的先后
#   与本表顺序可能不一致 → 柱顶数字会整体错位（零报错）。并列时直接停。
sig_top <- data.table(sig = names(sizes)[seq_len(n_show)],
                      n = as.integer(sizes[seq_len(n_show)]))
if (anyDuplicated(sig_top$n)) {
  stop(sprintf("前 %d 个交集存在并列大小（%s）→ 柱顶标注会错位，请调 --top",
               n_show, paste(unique(sig_top$n[duplicated(sig_top$n)]), collapse = ", ")))
}

p <- upset(
  as.data.frame(m), tools,
  name = "Classification tools",
  # ⭐⭐ user 2026-10-14（第 7 次修订）「中间空白太多了」→ **加宽 set-size 列吃掉空白**。
  #   ⚠⚠ 机制（`upset()` L261-264）：`plot_layout(widths = c(width_ratio, 1 - width_ratio))`
  #     → 这个数是**两列「面板」宽之比**（不是像素，也不含坐标轴区）。
  #     实测：`width_ratio=0.26` → 面板 542 : 1543 px（= 0.26 : 0.74 ✓）。
  #   ⚠⚠⚠ **「中间空白」由两块拼成，且都跟 `width_ratio` 无关**（画布 2796×1800，逐行实测）：
  #     ① 列 2 的**左侧坐标轴区** —— 顶部交集柱面板的 y 轴刻度文字（`5000`…`0`）住在那里，
  #        ~80 px + 两侧 margin ≈ **150 px，压不掉**（它是柱图的 y 轴，必须留）。
  #     ② 矩阵那串**行标签槽** —— 只有把 `axis.text.y` 真正 `element_blank()` 掉才会归零
  #        （`colour="transparent"` 仍占 205 px 宽）。
  #     → 实测四种组合（同源对拍，Family，`--no-timing`）：
  #         透明标签 + wr 0.26 → 空白 **297 px**（矩阵面板 1054..2488）
  #         透明标签 + wr 0.45 → 空白 **297 px**（**不变**：间隙是结构性的，只把整体右移）
  #         blank 标签 + wr 0.26 → 空白 **149 px**（矩阵面板被拉宽并左移到 945..2488）
  #         **blank 标签 + wr 0.312 → 空白 150 px（矩阵面板精确回到 1054..2488）** ← 采纳
  #     ⭐ 正解 = **两件事一起做**：`element_blank()` 消掉 205 px 的透明标签槽，
  #       再把 `width_ratio` 从 0.26 提到 **0.312** —— 后者把矩阵面板**推回原位**
  #       （补偿 blank 带来的"面板变宽+左移"），同时 set-size 面板净加宽 **147 px**
  #       （254..757 → 254..904，+29 %）。⇒ **点阵逐像素不动，空白从 297 → 150 px（-50 %）**。
  #     ⚠⚠ 两个参数**互相绑定**，不能只改一个：只 blank → 点阵左移 109 px（user 已否决）；
  #       只调 wr → 空白一点不减。改任一个都必须重测矩阵面板位置。
  #     ⚠ 试过且**无效**：把 set-size 面板右 margin / 矩阵左 margin 归零（`margin(5,0,5,5)`）
  #       → 只省 **9 px**（150 → 141）—— 因为列 2 的左轴区由**顶部柱面板**的 y 轴决定，
  #       矩阵面板自己的 margin 不是瓶颈。故不留这段无用复杂度。
  width_ratio = as.numeric(get_arg("--width-ratio", "0.312")),
  min_size = min_size,
  sort_intersections = "descending",
  sort_intersections_by = "cardinality",
  # ⚠⚠ 必须 `encode_sets = FALSE`（README 硬规则 18）：工具名含下划线（diamond_lca）
  #   不是合法 R 名 → 默认 `encode_sets = TRUE` 会把**集合编码成整数 id**，
  #   于是左侧 set-size 面板的离散轴标签变成 `1..7` 的乱码数字
  #   （实测首版显示 "6 1 4 2 7 5 3"），且**零报错**。
  #   关掉后名称原样保留，由**最左侧 set-size 面板**显示一次
  #   （矩阵面板那份重复的已由下面 `themes[["intersections_matrix"]]` 关掉）。
  encode_sets = FALSE,
  # ⭐⭐ v22.2（2026-10-14 用户「左侧条形图配色调整」）：set-size 条按**工具身份**上色。
  #   ⚠⚠⚠ 只能用「颜色向量」，**不能**写 `aes(fill = ...)`（三种写法全废，实测）：
  #     `upset_set_size(mapping = aes(fill = x))` → `object 'x' not found`；
  #     `geom_bar(aes(fill = set))`             → `object 'set' not found`；
  #     `geom_bar(aes(fill = x)) + scale_...`   → `Cannot add <ggproto> objects together`。
  #   ✅ `geom_bar(fill = <向量>)` 生效；向量按**传入 `sets`（= tools）的顺序**逐条对应，
  #      **与 `sort_sets` 造成的渲染顺序无关**（默认 descending 下实测确认）。
  #   ⚠ 工具顺序 = `tools`（= TOOLS7 交集，固定序），故直接用 `TOOL_COLORS[tools]`。
  set_sizes = upset_set_size(
      geom = geom_bar(
        # ⚠⚠⚠ 必须用 `set_size_colors_filtered()`（不是 `set_size_colors()`）：
        #   `min_size` 会把「在保留交集里从不出现」的集合**静默裁掉**
        #   （Species 的 `genomad` 仅 13 条 vs min_size=68 → 被裁，7 色 vs 6 行
        #    → `Aesthetics must be either length 1 or the same as the data (6)`）。
        #   判定签名 = `sig_top$sig`（= 本图实际显示的 n_show 个交集）。
        fill = if (want_tool_color)
                 set_size_colors_filtered(tools, TOOL_COLORS, sig_top$sig)
               else "grey35",
        # ⚠ 宽度必须保持 0.6（= `upset_set_size()` 默认 geom_bar 的宽度）：
        #   写成 0.62 会让条宽变化 → 左侧 3,855 px 差异（实测）。
        #   本行只为**换填充色**，几何必须逐像素不动。
        width = 0.6)) +
    # ⭐⭐ user 2026-10-14（第 7 次修订）「红色区域翻转」→ **柱条左右镜像**：
    #   由「右端对齐刻度 0、柱向左生长」改为「**左端对齐刻度 0、柱向右生长**」。
    #   ⚠⚠ 机制：`upset()` L220-222 对**左侧** set-size 面板默认给 `scale_y_reverse()`
    #     （`default_scale`），再 `coord_flip()` → 横轴反向 → 柱从右缘往左长。
    #     L228 是 `scale_if_missing(set_sizes, axis = "y", scale = default_scale)`，
    #     而 `scale_if_missing()` = 「`get_scale(annotation,"y")` 非空就返回 `list()`」
    #     → **只要我们自己挂一个 y 标度，默认的 `scale_y_reverse()` 就完全不加**。
    #   ⚠ `expansion(mult = c(0, 0.05))`：**非 reverse 标度下按常规方向**（0 端不留白、
    #     高端留 5 %）。⚠ 若写反，柱起点会离开刻度 0（README 硬规则 ③ 的反面场景）。
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
    # ⚠ 坐标轴标题默认是 `group`（ComplexUpset 的通用措辞）→ 改为可读的语义名。
    #   集合名已在矩阵面板显示 → 这里把 y 轴标题置空、只留数值轴标题。
    xlab("Set size") +
    theme_upset(8) +
    theme(axis.title.y = element_blank()),
  base_annotations = list(
    # ⚠⚠⚠ **绝不能用 `counts = TRUE` 的内置柱顶文字**（用户 2026-10-13 报
    #   「第一条柱子的数字标在了柱子里面」）。
    #   实测（`ggplot_build()` 读回，Family 版）：**第一根柱（x=1，高 5,443）
    #   的文字被放在 y = 4,626.55 = 柱高的 85 % 处 —— 即柱内**，
    #   而其余柱（2,711 / 1,446 / …）都用**真实柱高**。
    #   → 根因 = ComplexUpset 对「全图最大值那根柱」特殊处理（防止文字顶出 panel），
    #     把它压到柱内；白字压深色柱，**零报错**，只有肉眼才看得出。
    #   ⚠ 试过且**无效**的修法：把 `scale_y_continuous` 的 `expand` 从 0.16 提到 0.26
    #     （顶部净空 82 px → 210 px）——文字仍在柱内，因为偏移是**柱高比例**不是像素。
    #   → 正解 = **关掉内置文字，自己按真实柱高 `annotate()`**（对齐 x = 1..N）。
    "Agreed vOTUs" = intersection_size(
      counts = FALSE,
      mapping = aes(fill = as.factor(count))
    ) +
    annotate("text",
             x = seq_len(nrow(sig_top)),
             y = sig_top$n,
             label = format(sig_top$n, big.mark = ","),
             vjust = -0.5, size = 2.4, colour = "#222222") +
    # ⭐ v22.2：顶部交集柱改用**三张图统一的数量色阶**（见 upset_common.R 的
    #   `QUANT_COLORS` / `quant_colors()`）。⚠ 旧 `DEGREE_COLORS`（1=浅黄…7=橙红）
    #   与左侧工具身份色撞了 4 个完全同色 → 故按用户「三张统一」的要求换掉。
    #   `--quant-scheme degree` 可一键回到旧口径（零回归对拍用）。
    scale_fill_manual(values = switch(quant_scheme,
                               "degree" = DEGREE_COLORS,
                               "blues"  = QUANT_BLUES,
                               stats::setNames(quant_colors(1:7), as.character(1:7))),
                      name = "Tools agreed", drop = FALSE) +
    # 顶部余量给足（文字在柱顶之上，需容纳 2.4 pt 字 + 偏移）
    scale_y_continuous(expand = expansion(mult = c(0, 0.10))) +
    # ⚠ user 2026-10-13「去除 upset 图 Y 轴的 Agreed vOTUs」。
    #   ⚠⚠⚠ **踩了 3 次坑才定案**（正解 = 就写在这里！）：
    #     ① 删 `ylab()` → 无效（该标题不由这里产出，是 `upset()` L176 内部补写的）
    #     ② `themes` 里写 key `"intersections"` → **静默失效**（key 名不匹配任何面板）
    #     ③ `themes` 里写 key `"Agreed vOTUs"` → **仍无效**！虽然探针实测该 theme
    #        确实被赋到 plot[[2]]（`axis.title.y = element_blank` 已在该对象上），
    #        但 `upset()` L169-170 的应用顺序是：
    #            `annotation_plot + selected_theme + do.call(theme, user_theme)`
    #        `user_theme` = **本 list 里跟的 `theme_upset(8)`**（即下面这行），
    #        **排在 `selected_theme` 之后** → 它把 `axis.title.y` 覆盖回 `element_text`。
    #     ✅ 正解 = 在本 theme 里 `axis.title.y = element_blank()`（最后应用，必胜）。
    #   ℹ 「Agreed vOTUs」作为 list 名仍被用作 patchwork 面板名，保留不动。
    theme_upset(8) +
    theme(
      axis.text.x = element_blank(), axis.title.x = element_blank(),
      panel.grid.major.x = element_blank(),
      # ⭐ 唯一生效的关闭点（见上）。
      axis.title.y = element_blank()
    )
  ),
  stripes = upset_stripes(mapping = aes(color = "black", fill = "grey92")),
  # ⭐⭐ user 2026-10-14「（点阵矩阵的行标签）放到左侧可以吗」→ 行标签只保留**最左一处**。
  #   实测（本图逐像素看）：set-size 面板（最左）与矩阵面板**各显示一遍完全相同的集合名**
  #   → 本就是重复，去掉矩阵那份后 = 「行标签移到最左」。
  #   ⚠⚠ user 2026-10-14 第二轮纠正：「点阵整体[平移]保持不对，只把行标签[移到]最左处」
  #     → **点阵的 x 几何必须原封不动**；第 7 次修订要"吃掉中间空白"时，改用
  #     `element_blank()` + `width_ratio = 0.312` 的组合，**仍然把点阵钉回逐像素原位**
  #     （见下面 `themes[["intersections_matrix"]]` 与上面 `width_ratio` 两处的实测表）。
  #   ⚠⚠⚠ **不能写在这个图层链里**（`matrix = ... + theme(...)`，实测**无效**）：
  #     该 theme 会先被 `upset()` L98-102 的 `+ themes$intersections_matrix` 追加，
  #     而内置 `upset_themes$intersections_matrix` 是**两段**：
  #       [[1]] = 144 元素的**完整主题**（显式含 `axis.text.y = element_text(...)`）
  #       [[2]] = 3 元素微调（`axis.title.y` / `axis.text.x` / `axis.ticks.x` = blank）
  #     `[[1]]` 排在最后 → 把图层链里的写法**原样覆盖回去** → 毫无效果，**零报错**。
  #   ⚠⚠⚠ 也不能简单写 `themes[["intersections_matrix"]] = theme_upset(8) + theme(...)`：
  #     `themes[[key]]` 是**整体替换**语义 → 丢掉 `[[2]]` 的内置 blank
  #     → 底部冒出**挤成一团的 x 轴集合名** + 竖排 `group` 标题（**零报错**，实测踩到）。
  #   ✅ 正解 = 见下面 `themes` 里的 `"intersections_matrix"` 项：**拼回内置两段**再叠加。
  matrix = intersection_matrix(geom = geom_point(shape = 21, size = 2.6, stroke = 0.8)) +
    scale_color_manual(values = c("TRUE" = "#D55E00", "FALSE" = "grey88"), guide = "none"),
  themes = upset_modify_themes(list(
    # ⚠⚠⚠ **key 必须 === `base_annotations` 的 list 元素名**（这里是 `"Agreed vOTUs"`）。
    #   写成 `"intersections"` 会**静默失效**：`upset_modify_themes()` 内部
    #   `c(..., upset_themes[setdiff(names(upset_themes), names(to_update))])`
    #   对不认识的 key 只是「附加」，不做任何报错 → theme 永不生效（实测踩坑 2 次）。
    #   ℹ 内置 key 全集 = `intersections_matrix` / `Intersection size` / `overall_sizes` / `default`
    #     其中 `"Intersection size"` 是**默认** list 名（我们自定义了名字，故不适用）。
    "Agreed vOTUs" = theme_upset(8) +
      theme(
        axis.text.y = element_text(size = 7.5, color = "black"),
        # ⚠ 这里的 `axis.title.y = element_blank()` **不生效**（会被 base_annotations
        #   里后应用的 `theme_upset(8)` 覆盖）→ 真正生效的那份写在上面
        #   `base_annotations` 的 theme 里。此处保留仅为「表达意图 + 万一上游改序」。
        axis.title.y = element_blank(),
        plot.margin = margin(4, 4, 2, 4)
      ),
    "overall_sizes" = theme_upset(8) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 7, color = "black"),
        panel.grid.major.x = element_blank()
      ),
    # ⭐⭐ user 2026-10-14「（点阵矩阵的行标签）放到左侧可以吗」→ 行标签只保留**最左一处**。
    #   ── 第 6 次修订（当时的目标：只藏起来、不动几何）──
    #     用过 `colour = "transparent"`：**文字仍占位**（grob 宽度照常计算）→ 所有 x 几何
    #     与原来**逐像素一致**，只是看不见 → 集合名只剩最左 set-size 面板那一处。
    #     ℹ 用 `"transparent"` 而不是 `"white"`：白字在非白底上会现形。
    #   ── 第 7 次修订（新目标：吃掉中间空白）──
    #     ⚠⚠ 透明标签**仍占 205 px 宽** → 它就是「中间空白」里除柱图 y 轴之外的另一半
    #       （实测：透明 = 空白 297 px；blank = 149 px）。
    #     ✅ 故改为 **`element_blank()`**（槽宽归还面板），并**同时**把 `width_ratio`
    #        从 0.26 提到 **0.312** 来补偿「面板变宽 + 左移」→ 矩阵面板**回到 1054..2488**，
    #        与第 6 次修订的产物**逐像素同位**（同源对拍实测），而 set-size 面板净宽 +147 px。
    #     ⚠⚠ 两个参数**互相绑定**：只 `element_blank()` → 点阵左移 109 px（user 已否决）；
    #        只调 `width_ratio` → 空白一点不减。改任一个都必须重测矩阵面板位置。
    #   ℹ 若将来又只需"藏起标签、不动几何"（例如把 set-size 面板调回窄版），退回
    #     `element_text(colour = "transparent")` 即可 —— 但那时空白会回到 297 px。
    "intersections_matrix" = Reduce(`+`, ComplexUpset::upset_themes$intersections_matrix) +
      theme(axis.text.y = element_blank())
  ))
)

p <- p +
  plot_annotation(
    title = sprintf("%s-level classification agreement across tools", rank),
    subtitle = sprintf(
      "Lycium barbarum vOTUs with a %s-level consensus (n = %s) - showing top %d of %d tool-agreement patterns; colour = number of agreeing tools",
      rank, format(nrow(m), big.mark = ","), n_show, n_all
    ),
    theme = theme(
      plot.title = element_text(size = 10.5, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = 7.5, color = "#444444", hjust = 0),
      plot.margin = margin(6, 8, 4, 8)
    )
  )

# ══════════════════════════════════════════════════════════════════
# 左上角：7 分类工具「全库耗时」条形图（用户 2026-10-13「时间的图呢，左上角没看到呀」）
# ══════════════════════════════════════════════════════════════════
# ⭐ 位置来源 = **本版图（干净基线）逐像素实测**，不是沿用任何旧坐标：
#   画布 2796 × 1800 px（9.32 × 6.0 in @300 dpi）。
#   标题 y 26..66、副标题 y 91..118（x0..880 列内）；
#   ⚠⚠ README 硬规则：「空白带必须用**本次改动后**的干净基线重测，禁沿用上一版坐标」。
#
# ⭐⭐ user 2026-10-13「时间的图 位置往下些，往右边些」→ 相对首版（x 0..880 / y 137..513）平移。
#   实测（**排除耗时图自身后**）`x 0..1010 / y 119..560` 内的障碍只有两处，都在 **x 815..1009**：
#     · y 284..305 = 交集柱面板的 y 轴刻度 `6000`
#     · y 527..541 = 同轴刻度 `4000`
#   → **右缘必须 ≤ 810**（否则撞刻度文字）、**下缘 ≤ 525**（否则撞 `4000`）。
#   故：**右移 +100 px、下移 +25 px**，框收窄为 **x 100..810（710 px）、y 162..525（363 px）**
#   = **2.3667 × 1.2100 in**。
#   ⚠⚠ 尺寸变了 → **预渲 PNG 必须按新尺寸重跑**（预渲资产不随参数更新，README 硬规则）。
#
# ⭐⭐⭐ user 2026-10-14（第 2 条）「把耗时条形图 再次调整红色区域的大小，使其与预期效果一致」
#   → **整块放大到用户红框尺寸**（本图第 4 次改耗时条框）。
#   红框映射（截图 1136×695 → 画布 2796×1800，用**柱面板 4 个 y 轴刻度**标定）：
#     `x_canvas = 46.7 + 2.424 * x_shot`、`y_canvas = 164.8 + 2.424 * y_shot`
#     （`6000` 刻度 shot y53.5 ↔ 画布 294.5；`0` 刻度 shot y389 ↔ 画布 1106.5；
#      四刻度间距 shot 111.8 / 画布 271 → **s = 2.424**，逐刻度残差 ≤ 1 px）
#   红框外缘 shot `x19.5..371.5 / y41.5..376.5` → 画布 **`x94..946 / y265..1076`**
#     （852 × 811 px = **2.8400 × 2.7033 in**，比例 1.05 —— 由"扁"(1.96) 变"近正方"）。
#   ⚠⚠⚠ **尺寸变了 → 预渲 PNG 必须重渲**（`--w 2.8400 --h 2.7033`），否则 `rasterGrob`
#     按框硬拉伸、**零报错**（README 硬规则 22）。
#   净空（`--no-timing` 底图逐像素实测）：**红框全区 0 墨迹**；
#     右缘 946 距柱面板 y 轴刻度（`x971..1009`）**25 px**；底缘 1076 距 set-size 面板（`y1177`）**101 px**。
#   ⚠ 沿革：第 1 次 `x0..880/y137..513` → 第 2 次 `x100..810/y162..525` → 第 3 次下移 `y222..585` → 本次放大。
TIM_L <- 94              / 2796; TIM_R <- 946  / 2796
TIM_T <- 1 - 265 / 1800;         TIM_B <- 1 - 1076 / 1800
tim_w_in <- (TIM_R - TIM_L) * 9.32
tim_h_in <- (TIM_T - TIM_B) * 6.0

timing_raster <- NULL
if (want_timing) {
  # ⚠⚠⚠ **默认名必须随 `--rank` 推导**（2026-10-14 踩到，零报错）：
  #   原先默认恒为 `cls_tool_timing_fig.png`（Species 那份）→ `--rank Family` 若**没显式传**
  #   `--timing-png`，会把 **Species 行序**的耗时图静默嵌进 Family 版
  #   （Family 正确降序 = … > genomad(654) > VITAP(607)；Species = … > VITAP(604) > genomad(13)
  #     → 末两行互换）。⚠ 肉眼只看「耗时条有没有画出来」完全看不出，**必须量蓝条长度**。
  #   实测：Family 正确资产长度序(上→下) = 447/1200/169/346/334/**874/1027**；
  #         错嵌 Species 资产 = …/**1027/874**（末两位互换）。
  #   ⭐ 与 `make_cls_tool_timing_fig.R` 的 `--out` 推导规则**必须完全一致**（同族硬规则）。
  tim_png <- get_arg("--timing-png",
                     file.path(script_dir(), "upset_data",
                               if (identical(rank, "Species")) "cls_tool_timing_fig.png"
                               else sprintf("cls_tool_timing_fig_%s.png", rank)))
  if (!file.exists(tim_png)) {
    stop(sprintf("缺耗时图预渲文件 %s，请先跑 make_cls_tool_timing_fig.R", tim_png))
  }
  timing_raster <- grid::rasterGrob(png::readPNG(tim_png), interpolate = TRUE,
    width = unit(tim_w_in, "in"), height = unit(tim_h_in, "in"))
  cat(sprintf("耗时图框：x %.3f..%.3f y %.3f..%.3f npc｜%.2f × %.2f in（预渲 PNG）\n",
              TIM_L, TIM_R, TIM_B, TIM_T, tim_w_in, tim_h_in))
} else {
  cat("已跳过耗时图（--no-timing）\n")
}

# ══════════════════════════════════════════════════════════════════
# 右上角：7 分类工具 Venn 图（user 2026-10-14 第 8 次修订「Q 也绘制 venn
#   补充到右上角的空白处」；参考 https://zhuanlan.zhihu.com/p/434441906）
# ══════════════════════════════════════════════════════════════════
# ⭐ 实现选型（实测对比）：
#   · 参考文章用 `venn` 包（2–7 组）——但 7 组时**只给区域数字、无集合轮廓的清晰分层**，
#     且 7 集共 127 个区域、77 个非空 → 数字全画必糊成一团（实测两版皆不可读）。
#   · `ggVennDiagram`（`shape_id="701"`，与 ③ BPanel 的 Venn **同族同风格**）可
#     按区域取质心坐标 → 才能做「阈值过滤 + 白底标签 + 防重叠」。
#   ⇒ 采纳 ggVennDiagram，与 P 面板风格一致（灰平填充 #F2F2F2 + 黑轮廓）。
#
# ⭐⭐ 可读性定案（用户需求「补充到空白处」⇒ 必须能读，不能只是贴上去）：
#   7 集 Venn 的区域标签有三个坑，逐个实测解决：
#     ① 全标（127/77 个）→ 完全糊死（实测 `venn_gg7.png`）。
#     ② 只留质心 + 阈值过滤 → 仍**压线**（数字与花瓣轮廓交叠，实测 `c200.png`）。
#     ③ ggrepel 强 repel（box.padding 0.4）→ 标签被推离质心、**撞上集合名**（实测）。
#   ✅ 正解 = **geom_label_repel（白底无边框 linewidth=0）**：白底遮住穿过的轮廓线、
#      `force=1.2` 轻度位移（不飞出去）+ `segment` 细引线回指质心。
#     ⚠ `stroke`/`stroke.colour` 不是 `geom_text` 参数（实测被静默忽略，只有 warning）；
#       `label.size` 在 ggplot2 ≥3.5 已弃用 → 必须用 `linewidth`。
#   ⭐ 阈值取「**前 N 大**」（N=15）而非固定数值 —— Species（49 非空，top12=165）与
#     Family（77 非空，top12=302）分布差异大，固定阈值会让两版标注数差一倍。
#     实测 N=15 → Species 阈值 94 / Family 阈值 211，两版都是 15 个标签。
#   ⚠ 中心区域（1/2/3/4/5/6/7 全票）实测 count = **0**（Family）→ 不标注，无误导。
VENN_TOP_N <- as.integer(get_arg("--venn-top", "15"))
VENN_COL_LABEL <- get_arg("--venn-color", "#111111")
want_venn <- !("--no-venn" %in% args)

# ── 框位置（px @300dpi，画布 2796×1800）：**两版（Species/Family）同时为空**，
#   含 14 px 净空；且在 panel Q 的标题裁切线（y=212）**之下** ⇒ 组合主图里不会被裁掉。
#   ⚠⚠ 判据（实测）：框内 `g<215` 的墨迹 = 0（`g<245` 的 13k 像素全是 4 条
#     极浅网格线，`g<220` 即为 0，不是文字/柱体）。
#   ⚠ 右缘 2500 距图例色块（x2538 起）38 px；左缘 1390 在最高柱（x1066..1136）
#     与次高柱（x1144..1214，Species 版 1147..1226）右侧之外。
VEN_L <- 1390 / 2796; VEN_R <- 2500 / 2796
VEN_T <- 1 - 220 / 1800; VEN_B <- 1 - 920 / 1800
ven_w_in <- (VEN_R - VEN_L) * 9.32
ven_h_in <- (VEN_T - VEN_B) * 6.0

venn_raster <- NULL
if (want_venn) {
  # ⭐ 集合成员：直接用同一张 agree 矩阵的 TRUE/FALSE 列（与 UpSet 主体**同源**，
  #   保证 Venn 与 UpSet 逐区域一致，不引入第二套口径）。
  vsets <- lapply(tools, function(t) which(m[[t]]))
  names(vsets) <- tools
  vdat <- as.data.frame(
    ggVennDiagram:::process_data(ggVennDiagram::Venn(vsets), shape_id = "701")$regionLabel)
  vdat <- vdat[order(-vdat$count), ]
  n_lab <- min(VENN_TOP_N, sum(vdat$count > 0))
  vlab <- vdat[seq_len(n_lab), c("X", "Y", "count")]
  vlab$lab <- format(vlab$count, big.mark = ",")
  # ⚠ 自检：阈值必须是"前 N 大"的自然断点（并列时 ggrepel 仍能分开，但记录实际阈值）
  vthr <- vlab$count[n_lab]
  cat(sprintf("Venn：%d 集｜非空区域 %d｜标注前 %d 个（阈值 %s）｜中心区域 count=%s\n",
              length(tools), sum(vdat$count > 0), n_lab, format(vthr, big.mark = ","),
              vdat$count[vdat$id == paste(seq_along(tools), collapse = "/")]))

  pv <- ggVennDiagram::ggVennDiagram(
        vsets, shape_id = "701", label = "none",
        set_size = 1.9, set_color = "#333333", edge_size = 0.35) +
    # ⚠ 平填充（README 配色纪律：数值渐变禁碰语义色）→ 统一 #F2F2F2，与 P 面板一致。
    scale_fill_gradient(low = "#F2F2F2", high = "#F2F2F2", guide = "none") +
    geom_label_repel(
        data = vlab, aes(x = X, y = Y, label = lab),
        size = 1.45, colour = VENN_COL_LABEL, fontface = "bold",
        fill = "white", linewidth = 0, label.padding = unit(0.02, "lines"),
        min.segment.length = 0.3, segment.size = 0.15, segment.colour = "#999999",
        box.padding = 0.16, point.padding = 0.03, force = 1.2,
        max.overlaps = Inf, seed = 42) +
    labs(title = "Tool-agreement overlap (Venn)") +
    theme_void(base_size = 8) +
    theme(
      plot.title = element_text(size = 9.0, face = "bold", hjust = 0,
                                margin = margin(2, 0, 2, 4)),
      plot.margin = margin(1, 3, 1, 3)
    )

  tmpv <- tempfile(fileext = ".png")
  ggsave(tmpv, pv, width = ven_w_in, height = ven_h_in, dpi = 300, bg = "white")
  venn_raster <- grid::rasterGrob(png::readPNG(tmpv), interpolate = TRUE,
    width = unit(ven_w_in, "in"), height = unit(ven_h_in, "in"))
  cat(sprintf("Venn 图框：x %.3f..%.3f y %.3f..%.3f npc｜%.2f × %.2f in（预渲 PNG）\n",
              VEN_L, VEN_R, VEN_B, VEN_T, ven_w_in, ven_h_in))
} else {
  cat("已跳过 Venn（--no-venn）\n")
}

# ── 保存：gtable + grid（overlays 逐层叠 rasterGrob）──
# ⭐ 为什么不用 `inset_element()`：作用于 patchwork 会**静默误解析框 npc**
#   （IdentFilter 第 6 次修订实测，框被缩了一半且零报错）→ 唯一零 reflow 的画法是
#   「预渲 PNG → rasterGrob → 居中 vp 叠在 `patchwork:::patchworkGrob(p)` 上」。
save_cls_upset <- function(p, outdir, name, width, height, overlays = list()) {
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  # 只负责往**当前已打开的设备**上画; 设备开关交给 draw_figure_formats
  draw <- function() {
    grid::grid.newpage()
    grid::grid.draw(patchwork:::patchworkGrob(p))
    for (ov in overlays) {
      if (is.null(ov$grob)) next
      grid::pushViewport(grid::viewport(x = ov$cx, y = ov$cy,
                                        width = ov$w, height = ov$h, just = "center"))
      grid::grid.draw(ov$grob)
      grid::popViewport()
    }
    invisible(NULL)
  }
  # ⚠ 尺寸用**逻辑英寸**（与 save_upset 的 width/height 同值）→ 两图版式一致。
  # 2026-10-20: 四格式统一 (SVG/PNG/PDF/PNG@1200dpi)
  draw_figure_formats(draw, file.path(outdir, paste0(name, ".pdf")),
                      width = width, height = height, bg = "white")
  invisible(file.path(outdir, paste0(name, ".png")))
}

overlays <- list()
if (!is.null(timing_raster))
  overlays[[length(overlays) + 1]] <- list(grob = timing_raster,
    cx = (TIM_L + TIM_R) / 2, cy = (TIM_B + TIM_T) / 2,
    w = TIM_R - TIM_L, h = TIM_T - TIM_B)
if (!is.null(venn_raster))
  overlays[[length(overlays) + 1]] <- list(grob = venn_raster,
    cx = (VEN_L + VEN_R) / 2, cy = (VEN_B + VEN_T) / 2,
    w = VEN_R - VEN_L, h = VEN_T - VEN_B)

save_cls_upset(p, outdir, paste0("Fig_GojiVirome_UpsetClassification", suffix),
               width = 9.32, height = 6.0, overlays = overlays)
cat("Written:", file.path(outdir, paste0("Fig_GojiVirome_UpsetClassification", suffix, ".png")), "\n")
