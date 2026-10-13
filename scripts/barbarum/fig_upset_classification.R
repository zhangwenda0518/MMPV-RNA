#!/usr/bin/env Rscript
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
mtx <- get_arg("--matrix", sprintf("scripts/barbarum/upset_data/agree_matrix_%s.tsv", rank))
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
  width_ratio = 0.26,
  min_size = min_size,
  sort_intersections = "descending",
  sort_intersections_by = "cardinality",
  # ⚠⚠ 必须 `encode_sets = FALSE`（README 硬规则 18）：工具名含下划线（diamond_lca）
  #   不是合法 R 名 → 默认 `encode_sets = TRUE` 会把**集合编码成整数 id**，
  #   于是左侧 set-size 面板的离散轴标签变成 `1..7` 的乱码数字
  #   （实测首版显示 "6 1 4 2 7 5 3"），且**零报错**。
  #   关掉后名称原样保留，配合下面 overall_sizes 的 axis.text.y = element_blank()
  #   （集合名由矩阵面板统一显示一次）。
  encode_sets = FALSE,
  set_sizes = upset_set_size() +
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
    scale_fill_manual(values = DEGREE_COLORS, name = "Tools agreed", drop = FALSE) +
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
      )
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
TIM_L <- 100             / 2796; TIM_R <- 810  / 2796
TIM_T <- 1 - 162 / 1800;         TIM_B <- 1 - 525  / 1800
tim_w_in <- (TIM_R - TIM_L) * 9.32
tim_h_in <- (TIM_T - TIM_B) * 6.0

timing_raster <- NULL
if (want_timing) {
  tim_png <- get_arg("--timing-png",
                     file.path("scripts", "barbarum", "upset_data", "cls_tool_timing_fig.png"))
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

# ── 保存：gtable + grid（overlays 逐层叠 rasterGrob）──
# ⭐ 为什么不用 `inset_element()`：作用于 patchwork 会**静默误解析框 npc**
#   （IdentFilter 第 6 次修订实测，框被缩了一半且零报错）→ 唯一零 reflow 的画法是
#   「预渲 PNG → rasterGrob → 居中 vp 叠在 `patchwork:::patchworkGrob(p)` 上」。
save_cls_upset <- function(p, outdir, name, width, height, overlays = list()) {
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  one <- function(path, dev_open) {
    dev_open(path)
    grid::grid.newpage()
    grid::grid.draw(patchwork:::patchworkGrob(p))
    for (ov in overlays) {
      if (is.null(ov$grob)) next
      grid::pushViewport(grid::viewport(x = ov$cx, y = ov$cy,
                                        width = ov$w, height = ov$h, just = "center"))
      grid::grid.draw(ov$grob)
      grid::popViewport()
    }
    grDevices::dev.off()
    invisible(path)
  }
  # ⚠ 尺寸用**逻辑英寸**（与 save_upset 的 width/height 同值）→ 两图版式一致。
  one(file.path(outdir, paste0(name, ".png")),
      function(f) grDevices::png(f, width = width, height = height,
                                 units = "in", res = 300, bg = "white"))
  one(file.path(outdir, paste0(name, ".pdf")),
      function(f) grDevices::pdf(f, width = width, height = height, bg = "white"))
  invisible(file.path(outdir, paste0(name, ".png")))
}

overlays <- list()
if (!is.null(timing_raster))
  overlays[[length(overlays) + 1]] <- list(grob = timing_raster,
    cx = (TIM_L + TIM_R) / 2, cy = (TIM_B + TIM_T) / 2,
    w = TIM_R - TIM_L, h = TIM_T - TIM_B)

save_cls_upset(p, outdir, paste0("Fig_GojiVirome_UpsetClassification", suffix),
               width = 9.32, height = 6.0, overlays = overlays)
cat("Written:", file.path(outdir, paste0("Fig_GojiVirome_UpsetClassification", suffix, ".png")), "\n")
