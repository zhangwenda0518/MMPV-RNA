#!/usr/bin/env Rscript
# fig_ident_filter_combined.R — 图⑦「鉴定 x 过滤」双联图（**同图**）
#
#   左：7 鉴定轨道 x 过滤结局的**堆叠 UpSet**（= fig_ident_filter_upset.R 的图体）
#   右：**按工具**的过滤结局 100 % 堆叠柱（= fig_ident_filter_by_tool.R 的图体）
#
# 用户 2026-10-13 指令：「不应该和堆叠的 upset 图一起展示放到右侧嘛」
#   → 单张图、UpSet 在左、按工具堆叠柱在**右侧**，共用一套配色与一个图例。
#
# 数据：identfilter_data/ident_filter_upset_signature.tsv（与上述两图**同一份**，不引入新口径）。
#
# ⚠ 本脚本**自带**数据装载与两个面板的构建（有意不去改那两个已验证的单图脚本，
#   以免破坏它们的零回归；若口径变动，三处需同步）。

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(ComplexUpset)
  library(grid)
})

script_dir <- function() {
  f <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  # ⚠ 含中文路径下 normalizePath() 会把路径编坏 → 只用 dirname()
  if (length(f)) return(dirname(sub("^--file=", "", f[1])))
  getwd()
}
SD <- script_dir()
source(file.path(SD, "upset_common.R"))

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(k, d) { i <- which(args == k); if (length(i)) args[i + 1] else d }
datadir <- get_arg("--data", file.path(SD, "identfilter_data"))
outdir  <- get_arg("--outdir",
  file.path(SD, "..", "..", "reports", "_discovery_sci_20261012", "_lt",
            "10_Reports", "S10_SCI_Discovery"))
top_n   <- as.integer(get_arg("--top", "16"))

TRACKS  <- c("viralverify", "metabuli", "blast", "genomad", "virbot", "rdrpcatch",
             "viroid_blastn")
VIROID  <- "viroid_blastn"
STATUS_LEVELS <- c("Retained", "Removed by UniProt", "Removed by CDD", "Viroid track")
STATUS_COLORS <- c(
  "Retained"           = "#009E73",
  "Removed by UniProt" = "#BFBFBF",
  "Removed by CDD"     = "#7A7A7A",
  "Viroid track"       = "#E69F00"
)

# ⚠ signature 是 7 位 0/1 串 → fread 会当整数读、**吃掉前导零** → 必须按字符读。
sig <- fread(file.path(datadir, "ident_filter_upset_signature.tsv"), sep = "\t",
             colClasses = list(character = "signature"))
stopifnot(all(nchar(sig$signature) == length(TRACKS)))

n_all      <- nrow(sig)
n_nofilter <- sum(sig$n_nofilter)
n_plotted  <- sum(sig$n_total) - n_nofilter
n_ret      <- sum(sig$n_retained)
n_union    <- sum(sig$n_total) - sum(sig$n_viroid) - n_nofilter
n_cdd      <- n_ret + sum(sig$n_cdd_only)

# ── 逐工具求和（集合重叠，故列向求和会重复计数）──
per_tool <- rbindlist(lapply(seq_along(TRACKS), function(i) {
  sel <- substr(sig$signature, i, i) == "1"
  data.table(track = TRACKS[i],
             Retained             = sum(sig$n_retained[sel]),
             `Removed by UniProt` = sum(sig$n_cdd_only[sel]),
             `Removed by CDD`     = sum(sig$n_removed[sel]),
             `Viroid track`       = sum(sig$n_viroid[sel]),
             nofilter             = sum(sig$n_nofilter[sel]))
}))
per_tool[, total   := Retained + `Removed by UniProt` + `Removed by CDD` + `Viroid track`]
per_tool[, n_hits  := total + nofilter]
per_tool[, pct_ret := 100 * Retained / total]

# ═══════════════ 面板 A：堆叠 UpSet（无标题、无图例） ═══════════════
exp <- sig[, {
  st <- c(rep("Retained", n_retained), rep("Removed by UniProt", n_cdd_only),
          rep("Removed by CDD", n_removed), rep("Viroid track", n_viroid))
  .(status = st)
}, by = signature]
for (i in seq_along(TRACKS)) exp[[TRACKS[i]]] <- substr(exp$signature, i, i) == "1"
exp[, status := factor(status, levels = STATUS_LEVELS)]
stopifnot(nrow(exp) == n_plotted)

sizes    <- sort(sig$n_total, decreasing = TRUE)
n_show   <- min(top_n, n_all)
min_size <- as.integer(sizes[n_show])
cum_shown <- sum(sizes[seq_len(n_show)])

sig_top <- head(sig, n_show)
# ⚠ 柱顶标注按序号定位的前提：前 N 个交集大小互不相同（并列会整体错位）
if (anyDuplicated(sig_top$n_total)) {
  stop(sprintf("前 %d 个交集存在并列大小 → 柱顶标注会错位，请调 --top", n_show))
}
sig_top[, bar_h   := n_total - n_nofilter]
sig_top[, pct_lab := ifelse(n_viroid > 0, "viroid",
                            sprintf("%.0f%%", 100 * n_retained / pmax(bar_h, 1)))]

df <- as.data.frame(exp[, c(TRACKS, "status"), with = FALSE])
set_sizes_n <- vapply(seq_along(TRACKS), function(i) {
  sum(sig$n_total[substr(sig$signature, i, i) == "1"])
}, numeric(1))
names(set_sizes_n) <- TRACKS
# ⚠ 集合标签必须**单行**（两行总高超过点阵行距 → 与上一行压字）
LABELS <- setNames(sprintf("%s (%s)", TRACKS, format(set_sizes_n, big.mark = ",")), TRACKS)
setnames(df, TRACKS, unname(LABELS[TRACKS]))
# ⚠ ComplexUpset 会把 intersect 顺序**上下翻转** → 传 rev() 才能让最大集在最上面
SET_ORDER <- rev(TRACKS)

pA <- upset(
  df, unname(LABELS[SET_ORDER]),
  name = "Identification track",
  width_ratio = 0.22,
  height_ratio = 0.36,
  min_size = min_size,
  sort_sets = FALSE,
  set_sizes = FALSE,
  sort_intersections = "descending",
  sort_intersections_by = "cardinality",
  base_annotations = list(
    "Sequences" = intersection_size(counts = FALSE, mapping = aes(fill = status)) +
      annotate("text", x = seq_len(nrow(sig_top)), y = sig_top$bar_h,
               label = sig_top$pct_lab, vjust = -0.45, size = 1.9, colour = "#1A5C3A") +
      scale_fill_manual(values = STATUS_COLORS, name = NULL, drop = FALSE) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
      ylab("Sequences") +
      theme_upset(7.5) +
      # ⚠ 图例只在右侧面板出（见文末「图例只出一份」说明）→ 左侧必须显式 none，
      #   否则会得到两份一模一样的图例（patchwork 的 guides="collect" 进不到
      #   upset() 这个**嵌套 patchwork** 里面）。
      theme(axis.text.x = element_blank(), axis.title.x = element_blank(),
            panel.grid.major.x = element_blank(),
            legend.position = "none")
  ),
  stripes = upset_stripes(mapping = aes(color = "black", fill = "grey92")),
  matrix = intersection_matrix(geom = geom_point(shape = 21, size = 2.3, stroke = 0.7)) +
    scale_color_manual(values = c("TRUE" = "#D55E00", "FALSE" = "grey88"), guide = "none"),
  themes = upset_modify_themes(list(
    "intersections" = theme_upset(7.5) +
      theme(axis.text.y = element_text(size = 6.8, color = "black"),
            axis.title.y = element_text(size = 7.2),
            plot.margin = margin(4, 3, 2, 3)),
    "overall_sizes" = theme_upset(7.5) +
      theme(axis.text.x = element_blank(), axis.title.x = element_blank())
  ))
)

# ═══════════════ 面板 B：按工具 100 % 堆叠柱（无标题、无图例） ═══════════════
long <- melt(per_tool[, c("track", STATUS_LEVELS), with = FALSE], id.vars = "track",
             variable.name = "status", value.name = "n")
long[, status := factor(status, levels = STATUS_LEVELS)]
long[, frac := n / sum(n), by = track]

ord  <- per_tool[track != VIROID][order(-n_hits)]$track
lev  <- c(ord, VIROID)
xpos <- setNames(c(seq_along(ord), length(ord) + 2L), lev)
long[, x := xpos[track]]
per_tool[, x := xpos[track]]

xlab <- setNames(sprintf("%s\n%s", lev,
                         format(per_tool$n_hits[match(lev, per_tool$track)],
                                big.mark = ",")), lev)
xlab[VIROID] <- sub("^viroid_blastn", "viroid_blastn*", xlab[VIROID])

pB <- ggplot(long, aes(x = x, y = frac, fill = status)) +
  # ⚠ 因子首档默认在顶 → 要 Retained 在**柱底**必须 reverse
  geom_col(width = 0.70, position = position_stack(reverse = TRUE),
           colour = "white", linewidth = 0.25) +
  annotate("rect", xmin = xpos[VIROID] - 0.36, xmax = xpos[VIROID] + 0.36,
           ymin = 0, ymax = 1, fill = NA, colour = "#7A4F00",
           linetype = "dashed", linewidth = 0.45) +
  # ⚠ 不能加 limits = c(0,1)：limits 是数据范围过滤器 → 会把 y=1.035 的标注变 NA 丢掉
  geom_text(data = per_tool, aes(x = x, y = 1.035,
                                 label = sprintf("%.1f%%", pct_ret)),
            inherit.aes = FALSE, size = 2.1, colour = "#1A5C3A") +
  scale_fill_manual(values = STATUS_COLORS, name = NULL, drop = FALSE) +
  scale_y_continuous(breaks = seq(0, 1, 0.25),
                     labels = function(v) sprintf("%d%%", round(100 * v)),
                     expand = expansion(mult = c(0, 0.085))) +
  scale_x_continuous(breaks = unname(xpos[lev]), labels = unname(xlab[lev]),
                     limits = c(0.4, xpos[VIROID] + 0.6), expand = expansion(mult = 0)) +
  labs(x = NULL, y = "share of that tool's candidate contigs") +
  theme_upset(7.5) +
  theme(axis.text.x = element_text(size = 6.2, colour = "black", lineheight = 1.05),
        panel.grid.major.x = element_blank(),
        legend.position = "bottom", legend.key.size = unit(0.24, "cm"),
        legend.text = element_text(size = 6.6))

# ═══════════════ 合成 ═══════════════
title_grob <- function(title, subtitle) {
  gTree(children = gList(
    textGrob(title, x = 0, hjust = 0, y = 0.72,
             gp = gpar(fontsize = 10.5, fontface = "bold")),
    textGrob(subtitle, x = 0, hjust = 0, y = 0.20,
             gp = gpar(fontsize = 6.4, col = "#444444", lineheight = 1.15))
  ))
}

sub_txt <- paste0(
  "Lycium barbarum virome, 395 samples. Left: which tool(s) found each contig and what survived the filter chain. ",
  "Right: the same filter outcome as 100% composition per tool.\n",
  "union ", format(n_union, big.mark = ","), " viral candidates -> CDD ", format(n_cdd, big.mark = ","),
  " kept -> UniProt ", format(n_ret, big.mark = ","), " kept (", sprintf("%.1f", 100 * n_ret / n_union), "% of union). ",
  "UpSet shows top ", n_show, " of ", n_all, " intersections (", sprintf("%.2f", 100 * cum_shown / sum(sig$n_total)), "% of sequences); ",
  "right panel is per-tool membership (sets overlap, so columns sum to more than the union).\n",
  "viroid_blastn* is a separate viroid track: 88 of its 604 hits lie inside the viral union and enter the filter chain, the other 516 do not. ",
  "27 contigs from 5 samples without 02b output are excluded."
)

final <- wrap_elements(title_grob(
  "Identification x filtering of L. barbarum virome candidates", sub_txt
)) / (pA | pB) / guide_area() +
  plot_layout(heights = c(0.15, 1, 0.05), widths = c(2.3, 1), guides = "collect")
# ⚠ 不要用 `& theme(legend.position=...)` 收尾：`&` 会把主题应用到**所有**子图，
#   包括 pA → 会把 pA 里刚设的 legend.position="none" 覆盖掉 → 图例又变两份。

save_upset(final, outdir, "Fig_GojiVirome_IdentFilterCombined", width = 12.6, height = 8.4)

cat(sprintf("签名 %d｜入图 %s｜retained %s（%.2f%% of union）\n", n_all,
            format(n_plotted, big.mark = ","), format(n_ret, big.mark = ","),
            100 * n_ret / n_union))
cat(sprintf("UpSet 显示前 %d/%d 个交集（size >= %d），覆盖 %.2f%%\n",
            n_show, n_all, min_size, 100 * cum_shown / sum(sig$n_total)))
print(per_tool[, .(track, n_hits, Retained, `Removed by UniProt`, `Removed by CDD`,
                   `Viroid track`, pct_ret = round(pct_ret, 2))], row.names = FALSE)
cat("Written:", file.path(outdir, "Fig_GojiVirome_IdentFilterCombined.png"), "\n")
