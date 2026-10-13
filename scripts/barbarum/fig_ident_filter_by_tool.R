#!/usr/bin/env Rscript
# fig_ident_filter_by_tool.R — 图⑥「按鉴定工具」过滤结局 100% 堆叠柱（枸杞 discovery 侧）
#
# 与 fig_ident_filter_upset.R 同源、同色、同数据；只把 x 轴由「交集」换成「鉴定工具」：
#   每根柱 = 该工具**全部命中序列**的过滤结局构成，归一化到 100%。
#
# 为什么要归一化（而不是画绝对堆叠柱）：
#   7 条轨道的量级差 **750 倍**（viralverify 452,904 vs viroid_blastn 604）。
#   绝对堆叠柱放在同一线性轴上，后 4 根（genomad 11,526 / virbot 3,219 /
#   rdrpcatch 2,229 / viroid_blastn 604）只有 viralverify 的 0.1%–2.5% 高 →
#   肉眼不可见。构成比才是这张图的论点：**命中越少、越"专一"的工具，保留率越高**。
#   绝对数不丢：写在每个工具名下方 + 柱顶标保留率。
#
# 数据：identfilter_data/ident_filter_upset_signature.tsv
#   （7 位 0/1 签名 × 结局计数；逐工具取 signature[i]=='1' 的行求和）
#   与 UpSet 图**同一份数据**，不引入新口径。
# 口径：union 556,950 ⊇ CDD 473,523 ⊇ Uniprot 94,826（17.03%）；viroid∩union = 88。
#   5 个样本无 02b 产物（27 条）→ 过滤结局未知，**不计入柱**，副标题声明。

script_dir <- function() {
  f <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  # ⚠ 含中文路径下 normalizePath() 会把路径编坏 → 只用 dirname()
  if (length(f)) return(dirname(sub("^--file=", "", f[1])))
  getwd()
}
SD <- script_dir()
source(file.path(SD, "upset_common.R"))
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(grid) })

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(k, d) { i <- which(args == k); if (length(i)) args[i + 1] else d }
datadir <- get_arg("--data", file.path(SD, "identfilter_data"))
outdir  <- get_arg("--outdir",
  file.path(SD, "..", "..", "reports", "_discovery_sci_20261012", "_lt",
            "10_Reports", "S10_SCI_Discovery"))

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
stopifnot(all(c("signature", "n_total", "n_retained", "n_cdd_only", "n_removed",
                "n_viroid", "n_nofilter") %in% names(sig)))

# ── 逐工具求和（同一签名可同时计入多个工具 —— 集合本就重叠）──
per_tool <- rbindlist(lapply(seq_along(TRACKS), function(i) {
  sel <- substr(sig$signature, i, i) == "1"
  data.table(track = TRACKS[i],
             Retained             = sum(sig$n_retained[sel]),
             `Removed by UniProt` = sum(sig$n_cdd_only[sel]),
             `Removed by CDD`     = sum(sig$n_removed[sel]),
             `Viroid track`       = sum(sig$n_viroid[sel]),
             nofilter             = sum(sig$n_nofilter[sel]))
}))
# total = 过滤结局**已知**的命中数（= 柱内四段之和，不含 nofilter）
# n_hits = 该工具**原始命中总数**（含结局未知者）—— 与报告口径一致，用于 x 轴标签
per_tool[, total   := Retained + `Removed by UniProt` + `Removed by CDD` + `Viroid track`]
per_tool[, n_hits  := total + nofilter]
per_tool[, pct_ret := 100 * Retained / total]
setorder(per_tool, -n_hits)

cat("--- 逐工具过滤结局 ---\n")
print(per_tool[, .(track, n_hits, nofilter, Retained, `Removed by UniProt`,
                   `Removed by CDD`, `Viroid track`, pct_ret = round(pct_ret, 2))],
      row.names = FALSE)
cat(sprintf("Σ(工具命中) = %s（重叠计数，非并集）｜nofilter 逐工具求和 = %d（去重后 27 条），不入柱\n",
            format(sum(per_tool$n_hits), big.mark = ","), sum(per_tool$nofilter)))

# ── 长表 + 归一化 ──
long <- melt(per_tool[, c("track", STATUS_LEVELS), with = FALSE], id.vars = "track",
             variable.name = "status", value.name = "n")
long[, status := factor(status, levels = STATUS_LEVELS)]
long[, frac := n / sum(n), by = track]
stopifnot(all(abs(long[, sum(frac), by = track]$V1 - 1) < 1e-9))

# ── x 位置：6 个病毒工具按总量降序 1..6，留 1 位间隔，类病毒第 8 位 ──
ord    <- per_tool[track != VIROID]$track
lev    <- c(ord, VIROID)
xpos   <- setNames(c(seq_along(ord), length(ord) + 2L), lev)
long[, x := xpos[track]]
per_tool[, x := xpos[track]]
GAPX   <- length(ord) + 1L

xlab <- setNames(sprintf("%s\n%s", lev,
                         format(per_tool$n_hits[match(lev, per_tool$track)],
                                big.mark = ",")), lev)
xlab[VIROID] <- sub("^viroid_blastn", "viroid_blastn*", xlab[VIROID])

p <- ggplot(long, aes(x = x, y = frac, fill = status)) +
  # ⚠ 想「Retained 在柱底」→ 因子首档默认在顶 → 必须 position_stack(reverse = TRUE)
  geom_col(width = 0.70, position = position_stack(reverse = TRUE),
           colour = "white", linewidth = 0.25) +
  # 类病毒轨道：橙色虚线框 + 间隔位，强调「独立轨道，不在病毒并集内」
  annotate("rect", xmin = xpos[VIROID] - 0.36, xmax = xpos[VIROID] + 0.36,
           ymin = 0, ymax = 1, fill = NA, colour = "#7A4F00",
           linetype = "dashed", linewidth = 0.45) +
  # 柱顶：最终保留率（全图的定量结论）
  geom_text(data = per_tool, aes(x = x, y = 1.035,
                                 label = sprintf("%.1f%%", pct_ret)),
            inherit.aes = FALSE, size = 2.3, colour = "#1A5C3A") +
  scale_fill_manual(values = STATUS_COLORS, name = NULL, drop = FALSE) +
  scale_y_continuous(breaks = seq(0, 1, 0.25),
                     labels = function(v) sprintf("%d%%", round(100 * v)),
                     expand = expansion(mult = c(0, 0.085))) +
  scale_x_continuous(breaks = unname(xpos[lev]), labels = unname(xlab[lev]),
                     limits = c(0.4, xpos[VIROID] + 0.6), expand = expansion(mult = 0)) +
  labs(
    x = NULL,
    y = "share of that tool's candidate contigs",
    title = "Filter outcome by identification tool",
    subtitle = paste0(
      "Lycium barbarum virome (395 samples). Each bar = 100% of one identification tool's candidate contigs;\n",
      "number under each tool = its total hits; value above each bar = share retained after CDD + UniProt.\n",
      "Retention is much higher for selective tools: viralverify 15.4% (452,904 hits) vs virbot 83.3% (3,219).\n",
      "viroid_blastn* = separate viroid track; only its 88 hits inside the viral union enter the filter chain."
    )
  ) +
  theme_upset(8.5) +
  theme(
    axis.text.x = element_text(size = 6.8, colour = "black", lineheight = 1.05),
    panel.grid.major.x = element_blank(),
    legend.position = "bottom",
    legend.key.size = unit(0.26, "cm"),
    legend.text = element_text(size = 7)
  )

save_upset(p, outdir, "Fig_GojiVirome_IdentFilterByTool", width = 9.32, height = 5.2)
cat("Written:", file.path(outdir, "Fig_GojiVirome_IdentFilterByTool.png"), "\n")
