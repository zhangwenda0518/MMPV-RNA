#!/usr/bin/env Rscript
# fig_ident_filter_stacked_all.R — 图⑧「一体式堆叠 UpSet」：原生 ComplexUpset 版式，
#   **左侧 set-size 面板 与 右上交集柱 两侧都按过滤结局堆叠**（用户 2026-10-13 选定）。
#
# 与图④（fig_ident_filter_upset.R）的差别：
#   图④ set_sizes = FALSE（左侧没有 set-size 面板）；
#   本图 set_sizes = upset_set_size(geom = geom_bar(aes(fill = status), ...))
#     → 左侧面板变成「每个鉴定工具命中序列的过滤结局构成」横向堆叠条。
#   ⚠ 这是把图⑥（按工具 100% 堆叠柱）**并进 UpSet 本体**的做法：
#     与图⑦（fig_ident_filter_combined.R，patchwork 左右拼两块）不同 —— 本图**只有一套坐标系**，
#     矩阵行、set-size 条、交集柱三者共享同一组 16 个交集位置，不是两张图并排。
#
# 集合（7）= 6 个有命中的病毒鉴定工具 + 1 条独立类病毒轨道：
#   viralverify / metabuli / blast / genomad / virbot / rdrpcatch / viroid_blastn
# 段 = 过滤结局：
#   Retained            通过 CDD + UniProt 两道过滤，最终保留
#   Removed by UniProt  通过 CDD、被 UniProt 滤掉
#   Removed by CDD      被 CDD 滤掉（未进入 cdd_kept）
#   Viroid track        viroid_blastn 命中且**不在**病毒过滤链内（独立轨道）
#
# 数据：identfilter_data/ident_filter_upset_signature.tsv（服务器从 02a/02b 逐样本汇总，395 样本）
# 口径：union 556,923（Σ557,466 剔除 27 条无 02b 产物）⊇ CDD 473,487 ⊇ UniProt 94,818 = 17.03%
#
# ⚠⚠ set-size 面板的 x 轴尺度（--set-size）：
#   stack = 绝对条长（忠实于原生 set-size 面板语义）；但本数据 7 条轨道相差 **750 倍**
#           （452,904 vs 604）→ 后 3 根在 2 in 宽面板里只有 0–4 px，段完全不可见。
#   fill  = 归一化 100%（每根等宽，构成完全可读）；绝对数由**行标签**（如 `virbot (3,219)`）承载。
#   两种都产出，默认 stack。
#
# 来源：数据 = 项目 discovery 管线产物；绘图代码为本次独立编写。

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
get_arg <- function(k, d) {
  i <- which(args == k)
  if (length(i)) args[i + 1] else d
}
datadir  <- get_arg("--data", file.path(SD, "identfilter_data"))
outdir   <- get_arg("--outdir",
  file.path(SD, "..", "..", "reports", "_discovery_sci_20261012", "_lt",
            "10_Reports", "S10_SCI_Discovery"))
top_n    <- as.integer(get_arg("--top", "16"))
ss_mode  <- get_arg("--set-size", "stack")          # stack | fill
if (!ss_mode %in% c("stack", "fill")) stop("--set-size 只能是 stack 或 fill")
# 右上角 6 工具 Venn 开关（用户 2026-10-13 追加）。默认开；`--no-venn` 关闭便于零回归对拍。
# ⚠ 布尔开关：argparse 共享 dest 时**先注册者的 default 生效** → 用「只查存在性」的写法更稳。
want_venn <- !("--no-venn" %in% args)
# 左上格 7 轨道耗时条形图开关（用户 2026-10-13 第 7 次修订追加）。默认开。
want_timing <- !("--no-timing" %in% args)
# 右上角图例 overlay 开关（用户 2026-10-13 第 8 次修订追加）。默认开；`--no-legend-ov` 关。
want_legend_ov <- !("--no-legend-ov" %in% args)
suffix   <- if (ss_mode == "fill") "_norm" else ""

TRACKS <- c("viralverify", "metabuli", "blast", "genomad", "virbot", "rdrpcatch",
            "viroid_blastn")
# 右上角 Venn 的 6 个集合 = **6 个病毒工具**（不含 viroid_blastn；用户 2026-10-13 选定）。
VENN_TRACKS <- c("viralverify", "metabuli", "blast", "genomad", "virbot", "rdrpcatch")
STATUS_LEVELS <- c("Retained", "Removed by UniProt", "Removed by CDD", "Viroid track")
STATUS_COLORS <- c(
  "Retained"           = "#009E73",  # 绿：通过两道过滤
  "Removed by UniProt" = "#BFBFBF",  # 浅灰：CDD 过、UniProt 滤掉
  "Removed by CDD"     = "#7A7A7A",  # 深灰：CDD 即滤掉
  "Viroid track"       = "#E69F00"   # 橙：独立类病毒轨道（不在病毒过滤链内）
)

# ⚠ signature 是 7 位 0/1 串（如 "0100000"）→ fread 会当整数读、**吃掉前导零**，
#   必须强制按字符读，否则 nchar 校验失败 / 位序整体错位。
sig <- fread(file.path(datadir, "ident_filter_upset_signature.tsv"), sep = "\t",
             colClasses = list(character = "signature"))
stopifnot(all(nchar(sig$signature) == length(TRACKS)))
stopifnot(all(c("signature", "n_total", "n_retained", "n_cdd_only", "n_removed",
                "n_viroid", "n_nofilter") %in% names(sig)))

n_all      <- nrow(sig)
n_nofilter <- sum(sig$n_nofilter)
n_plotted  <- sum(sig$n_total) - n_nofilter
n_ret      <- sum(sig$n_retained)
n_union    <- sum(sig$n_total) - sum(sig$n_viroid) - n_nofilter
n_cdd      <- n_ret + sum(sig$n_cdd_only)

cat(sprintf("签名数 %d｜序列（入图）%s｜其中 retained %s（%.2f%% of union）\n",
            n_all, format(n_plotted, big.mark = ","), format(n_ret, big.mark = ","),
            100 * n_ret / n_union))
cat(sprintf("剔除 nofilter = %d 条（5 个样本无 02b 产物）\n", n_nofilter))

# ── 展开：签名 -> 逐序列行（ComplexUpset 需要逐元素表）──
exp <- sig[, {
  st <- c(
    rep("Retained", n_retained),
    rep("Removed by UniProt", n_cdd_only),
    rep("Removed by CDD", n_removed),
    rep("Viroid track", n_viroid)
  )
  .(status = st)
}, by = signature]
for (i in seq_along(TRACKS)) {
  exp[[TRACKS[i]]] <- substr(exp$signature, i, i) == "1"
}
exp[, status := factor(status, levels = STATUS_LEVELS)]
stopifnot(nrow(exp) == n_plotted)
cat(sprintf("展开行数 %s\n", format(nrow(exp), big.mark = ",")))

# ── 显示阈值：前 TOP_N 个交集 ──
sizes <- sort(sig$n_total, decreasing = TRUE)
n_show <- min(top_n, n_all)
min_size <- as.integer(sizes[n_show])
cum_shown <- sum(sizes[seq_len(n_show)])
cat(sprintf("显示前 %d/%d 个交集（size >= %d），覆盖 %s 条 = %.2f%%\n",
            n_show, n_all, min_size, format(cum_shown, big.mark = ","),
            100 * cum_shown / sum(sig$n_total)))

# ── 柱顶标注用：前 TOP_N 个交集 ──
# ⚠ 自检：前 TOP_N 个交集大小必须互不相同 —— 有并列时 ComplexUpset 的先后与本表顺序
#   可能不一致，柱顶 % 会整体错位。
sig_top <- head(sig, n_show)
if (anyDuplicated(sig_top$n_total)) {
  stop(sprintf("前 %d 个交集存在并列大小（%s）→ 柱顶标注会错位，请调 --top 或改查表实现",
               n_show, paste(unique(sig_top$n_total[duplicated(sig_top$n_total)]), collapse = ", ")))
}
sig_top[, bar_h := n_total - n_nofilter]
sig_top[, pct := 100 * n_retained / pmax(bar_h, 1)]
sig_top[, pct_lab := ifelse(n_viroid > 0, "viroid",
                            sprintf("%.0f%%", pct))]

df <- as.data.frame(exp[, c(TRACKS, "status"), with = FALSE])

# ── 集合大小（记录口径：Σ n_total over 命中该工具的签名）──
set_sizes_n <- vapply(seq_along(TRACKS), function(i) {
  sum(sig$n_total[substr(sig$signature, i, i) == "1"])
}, numeric(1))
names(set_sizes_n) <- TRACKS

# ⚠⚠ 自检：set-size 面板实际统计的是**展开表**里命中该工具的行数，
#   而记录口径是 Σ n_total（含 n_nofilter 那 27 条）。两者必须一致，
#   否则面板上的条长与行标签里的数字对不上。
set_sizes_actual <- vapply(TRACKS, function(t) sum(df[[t]]), numeric(1))
cat("集合大小（记录口径 Σ n_total / 展开表行数 / 差）：\n")
print(data.frame(recorded = set_sizes_n, expanded = set_sizes_actual,
                 diff = set_sizes_n - set_sizes_actual))
# ⚠⚠ 两者必然不等：记录口径 Σ n_total 含 n_nofilter（5 个样本无 02b 产物），
#   而展开表**不含**那部分（无过滤结局）→ 逐工具差 19/5/8/0/0/1/0，合计 33。
#   ⚠ 33 ≠ 27：27 是**去重后**的 nofilter 序列数，33 是**按工具重复计数**（一条序列可被多工具命中）。
#   → set-size 面板画的是展开表行数，**行标签必须用展开口径**，否则标签与条长对不上。
LABELS <- setNames(
  sprintf("%s (%s)", TRACKS, format(set_sizes_actual, big.mark = ",")), TRACKS)
setnames(df, TRACKS, unname(LABELS[TRACKS]))

# ── 左侧 set-size 条：逐工具标保留率（用户 2026-10-13 追加「左侧也标保留率」）──
#   与右侧交集柱顶的 % 同一套语义：**保留率 = 绿段 / 该工具整根条**
#   （⚠ viroid_blastn 的 604 命中里 516 是类病毒轨道 → 它的 1.8 % 是「绿段占整根条」，
#     不是「占进入过滤链的 88 条」；后者是 12.5 %。副标题已声明 88/516 的拆分。）
#   ⚠⚠ 此处 df 的列名**已经**被换成带计数的标签（`setnames` 在上面）→ 必须用
#      `LABELS[[t]]` 取列，直接写 `df[[t]]` 会拿到 NULL、`sum()` 静默得 0（首版实测全 0）。
per_tool_ret <- vapply(TRACKS, function(t) sum(df[[LABELS[[t]]]] & df$status == "Retained"),
                       numeric(1))
names(per_tool_ret) <- TRACKS
ss_pct <- 100 * per_tool_ret / set_sizes_actual
cat("逐工具保留率（绿段 / 整根条）：\n")
print(data.frame(retained = per_tool_ret, hits = set_sizes_actual,
                 pct = round(ss_pct, 2)))
ss_lab <- data.frame(
  set = unname(LABELS[TRACKS]),
  # ⚠ 标注位置 = 该条**外端**（与右侧柱顶「标在柱顶」一致）：
  #   fill 模式每根等宽 → 外端都在同一处，标签自然排成一列；
  #   stack 模式外端 = 该工具命中数，最大那根贴面板左缘 → 故 y 轴要留出外端余量。
  y   = if (ss_mode == "fill") 1 else unname(set_sizes_actual),
  lab = sprintf("%.1f%%", ss_pct),
  stringsAsFactors = FALSE
)

# ⚠ ComplexUpset 会把 intersect 给的顺序**上下翻转**再画 → 传 rev() 才能让
#   viralverify（最大集）落在最上面。
SET_ORDER <- rev(TRACKS)

# ═══════════ 左侧 set-size 面板：横向堆叠 ═══════════
# ⚠⚠ 两个平台/版本无关的硬事实（实测，勿凭直觉改）：
#   ① `upset_set_size()` 的 `mapping=` 形参**不可用** —— 传 `mapping = aes(fill = x)`
#      会在 compute_aesthetics 阶段报 `object 'x' not found`（该 mapping 被当作**图层级**
#      之外的表达式求值）。**必须把 aes 写进 geom 里**：`geom_bar(aes(fill = x), ...)`。
#   ② set-size 面板的底层数据是 `upset_data()$presence`，**不随 min_size 裁剪**
#      （`filter_intersections` 默认 FALSE）→ 绝对条长是真集合大小，可以放心用。
# ⚠⚠⚠ `scale_y_reverse` 下 `expansion(mult = c(a, b))` 的**方向是反的**（实测，务必记住）：
#   两个数按「**变换后**的区间」分配 —— `trans = "reverse"` 把数据 [0, 1] 映成 [-1, 0]，
#   于是 `a` 落在变换后的下界 = 数据**上界 y = max**（即 set-size 条的**外端 / 左侧**），
#   `b` 落在变换后的上界 = 数据 0（**内侧 / 右侧**）。
#   → 保留率标签贴在 y = 外端，所以**余量必须写在第一个位置**。
#   ⚠ 写成 `c(0.004, 0.13)`（直觉上的"上界留白"）会把余量全部留到**右侧空处**，
#     左侧外端只剩 0.4 % → 标签被面板左缘裁掉，**且不报任何错、也不进 build 数据**
#     （`ggplot_build()` 里 x/y/label 三列都在、位置也对，只是渲染时被裁）——
#     本图首版 7 个标签全丢就是这么来的。实测对照：`c(0.004,0.13)` → 标签色 0 px；
#     `c(0.13,0.004)` → 97 px；再加到 `c(0.30,0.004)` 仍 97 px（说明方向对了之后余量已够）。
# ⚠ 余量取 16 %（不是"够用就行"）：标签宽度是**物理尺寸**（固定 px），余量是**面板宽的百分比**，
#   而服务器字形比本机宽约 12 %（实测）→ 必须留富余。
#   本机实测（fill 模式，3420 px 宽画布）：条长 737 px ↔ 数据 1 单位；面板宽 ≈ 857 px；
#   余量区 = 0.16 × 737 = 118 px；最长标签 "15.4%" 实宽 **64 px** → 左端富余 47 px；
#   服务器按 ×1.12 ≈ 72 px → 富余仍有 ~46 px。（曾误估标签仅 ~33 px，实为 64 px，勿凭估。）
# ⚠ fill 模式显式钉 breaks = 0/25/50/75/100：加了余量后 ggplot 自动 breaks 会漂成
#   0/30/60/90（实测），读数不如四分位直观。
ss_scale <- if (ss_mode == "fill") {
  scale_y_reverse(labels = function(v) sprintf("%.0f%%", 100 * v),
                  breaks = seq(0, 1, by = 0.25),
                  expand = expansion(mult = c(0.16, 0.004)))
} else {
  scale_y_reverse(labels = function(v) format(v, big.mark = ",", scientific = FALSE),
                  expand = expansion(mult = c(0.16, 0.004)))
}
# ⚠⚠ 堆叠顺序（用户 2026-10-13 追加：「保留的放里面」）：
#   两个面板都改成 **Retained 贴坐标轴（内侧）**，即
#     竖向交集柱 → 贴底边；横向 set-size 条 → 贴右端（该面板 x 轴是 reverse，0 在右）。
#   理由：① 两类柱子的"好结局"都从同一条基线起算 → **绿色段长度可直接横向比长短**
#            （原先绿色贴外端、各自终点对齐，长度反而不好比）；
#         ② Viroid track 作为独立轨道被推到最外端，不与病毒三态混在一起。
#   ⚠ 必须用 `position_*_reverse(TRUE)`：默认 stack 把**第 1 个因子水平**放在最外端。
#   ⚠ 竖向面板的 position 走 `intersection_size(position = )`（该形参存在，勿自己拼 geom）。
ss_pos <- if (ss_mode == "fill") position_fill(reverse = TRUE) else position_stack(reverse = TRUE)
set_sizes_ann <- upset_set_size(
    geom = geom_bar(aes(fill = status), width = 0.62, position = ss_pos),
    position = "left"
  ) +
  # ⚠⚠ 必须给 set-size 面板**单独**补 scale_fill_manual —— 它和交集柱面板是两个独立 ggplot，
  #   不共享 scale。漏了就会用 ggplot 默认色轮（红/橄榄/青/紫），与图例完全对不上（实测首版）。
  scale_fill_manual(values = STATUS_COLORS, drop = FALSE) +
  # 逐工具保留率标签（贴在条的外端；与右侧柱顶 % 同色同字号）
  # ⚠ 图层自带 data → 不受 `set_sizes %+% set_sizes_data` 换数据影响。
  geom_text(data = ss_lab, aes(x = set, y = y, label = lab),
            hjust = 1.08, vjust = 0.5, size = 2.05, colour = "#1A5C3A") +
  ss_scale +
  # ⚠ 轴标题写在 flip 之前：`ylab` 落在翻转后的**下边**（数值轴），
  #   `xlab` 落在翻转后的左边（集合名轴，已由 overall_sizes 主题置空）。
  ylab(if (ss_mode == "fill") "Share of hits per identification track"
       else "Hits per identification track") +
  xlab(NULL) +
  # ⚠ 图例只在右上交集柱面板出 → 左侧必须显式 none，否则同一份图例会出两次。
  theme(legend.position = "none")

# ═══════════ 本体 ═══════════
p_body <- upset(
  df, unname(LABELS[SET_ORDER]),
  name = "Identification track",
  width_ratio = 0.30,
  min_size = min_size,
  sort_sets = FALSE,
  # ⚠⚠ 必须 encode_sets = FALSE：集合名里含空格 / 括号 / 逗号（"viralverify (452,885)"）
  #   不是合法 R 名 → 默认 encode_sets = TRUE 会把集合**编码成整数 id**，
  #   于是 set-size 面板自己的离散轴标签变成 1..7 的乱码数字（实测首版左侧 "5 3 1 2 6 4 7"）。
  #   关掉后名称原样保留，配合下面 overall_sizes 的 axis.text.y = element_blank()
  #   （集合名由矩阵面板统一显示一次，与参考版式一致）。
  encode_sets = FALSE,
  set_sizes = set_sizes_ann,
  sort_intersections = "descending",
  sort_intersections_by = "cardinality",
  base_annotations = list(
    "Sequences" = intersection_size(
      counts = FALSE,
      mapping = aes(fill = status),
      # ⚠⚠ 让 Retained 贴底边（内侧），与左侧 set-size 条贴轴一致 —— 见上文 ss_pos 注释。
      position = position_stack(reverse = TRUE)
    ) +
      # 每根柱顶标「最终保留率」—— 全图的定量结论：工具共识越高，保留率越高
      # ⚠ x 位置 = 交集在 descending-cardinality 序下的 1..N；实测前 N 个交集
      #   大小互不相同（见上自检），ComplexUpset 的排序与 sig 表顺序一致。
      annotate("text",
               x = seq_len(nrow(sig_top)), y = sig_top$bar_h,
               label = sig_top$pct_lab,
               vjust = -0.45, size = 2.05, colour = "#1A5C3A") +
      scale_fill_manual(values = STATUS_COLORS, name = NULL, drop = FALSE) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.10))) +
      # ⚠⚠⚠ **第 10 次修订（用户 2026-10-13 21:03：「upset 的 Y 轴 Sequences 不显示」）**：
      #   原生 `ylab("Sequences")` **确实渲染了，但位置错**：ComplexUpset 的 intersections 面板
      #   在本 patchwork 里带着**巨大的左侧 margin**（set-size 列 + spacer 列），
      #   axis.title.y 按「整个面板含左 margin」定位 → 题注飘到 **x≈952..984**，
      #   而 y 轴刻度文字 `3e+05/2e+05/1e+05/0e+00` 在 **x≈1257..1499** →
      #   **题注离轴 386 px（≈1.3 in）**，孤零零浮在耗时图下方空白里 →
      #   用户视觉上判为「Y 轴题不显示」（其实在，只是脱轴 + 与耗时图 x 轴题挤在一起）。
      #   ⭐ 实测干净空档（本版渲染图逐像素扫）：**x1190..1250 在 y249..1582 全空** →
      #     正好是「紧贴刻度文字左缘、纵向覆盖整根柱面板」的题注位。
      #   → 正解 = **关掉原生 ylab，改由 overlay 在实测位画题注**（与本图 Venn/耗时/图例同一套
      #     「预渲 → rasterGrob」零 reflow 机制；见文件末尾 `seq_title_raster`）。
      #   ⚠⚠⚠ **踩坑：`labs(y = NULL)` 与 `axis.title.y = element_blank()` 都关不掉它！**
      #     ComplexUpset 源码（`upset()` 第 168–171 行）在拼 patchwork 时会**条件性补写**：
      #       `if (is.null(labels$y) || default_y == labels$y) annotation_plot + ylab(name)`
      #     `labs(y = NULL)` → `labels$y` 为 NULL → **条件成立 → 它又把 `ylab("Sequences")` 加回来**
      #     （theme 的 element_blank 只控制"怎么画"，控制不了"标签内容存在" → 实测两个题注并存）。
      #   ⭐ 正解 = **给一个非空、且 ≠ 默认 y 映射的哨兵串**（这里用单个空格）：
      #     `labels$y` 非 NULL 且 ≠ default_y → ComplexUpset **跳过** `ylab(name)`；
      #     再由下面 themes 里的 `axis.title.y = element_blank()` 把这个空格题注**隐藏**。
      labs(y = " ") +
      theme_upset(8) +
      theme(
        axis.text.x = element_blank(), axis.title.x = element_blank(),
        panel.grid.major.x = element_blank(),
        # ⚠ 图例**不在这里**：见下面「图例 → 左上角空白格」。
        #   这里必须显式 none，否则一份图例会出两次（本面板 + 左上格）。
        legend.position = "none"
      )
  ),
  stripes = upset_stripes(mapping = aes(color = "black", fill = "grey92")),
  matrix = intersection_matrix(geom = geom_point(shape = 21, size = 2.8, stroke = 0.8)) +
    scale_color_manual(values = c("TRUE" = "#D55E00", "FALSE" = "grey88"), guide = "none"),
  themes = upset_modify_themes(list(
    "intersections" = theme_upset(8) +
      theme(
        axis.text.y = element_text(size = 7.5, color = "black"),
        # ⚠⚠⚠ **第 10 次修订**：这里用 `element_blank()` **隐藏**上一步 `labs(y = " ")` 设的哨兵题注。
        #   ⚠ 但**单靠它关不掉原生题注** —— `element_blank()` 只控制"怎么画"，
        #   控制不了 ComplexUpset 是否**写入** `ylab(name)`（那是 `labels$y` 的空否决定的）。
        #   → 两道缺一不可：① `labs(y = " ")` 让源码跳过补写；② 这里把哨兵隐藏。
        axis.title.y = element_blank(),
        plot.margin = margin(4, 4, 2, 4)
      ),
    "overall_sizes" = theme_upset(8) +
      theme(
        # ⚠⚠ 翻转后主题元素名指的是**渲染后**的轴，不是逻辑轴：
        #   axis.text.y / axis.title.y = 垂直轴（集合名，由矩阵面板统一显示 → 这里置空）
        #   axis.text.x / axis.title.x = 水平轴（数值，即 set-size 的条长）
        axis.text.y  = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x  = element_text(size = 6.6, colour = "black"),
        plot.margin  = margin(4, 2, 4, 6)
      )
  ))
)

title_grob <- function(title, subtitle, y_title, y_sub) {
  # ⚠⚠ 必须用 **vjust = 1（顶对齐）**，不能沿用系列里「y 当中心点」的写法：
  #   多行副标题以中心锚定时，块高的一半会落到锚点下方 → 行数一多就被标题带**下边缘裁掉**
  #   （实测 4 行版第 4 行只剩 11 px 高，正是被裁）。
  gTree(children = gList(
    textGrob(title, x = 0, hjust = 0, y = y_title, vjust = 1,
             gp = gpar(fontsize = 10.5, fontface = "bold")),
    textGrob(subtitle, x = 0, hjust = 0, y = y_sub, vjust = 1,
             gp = gpar(fontsize = 6.4, col = "#444444", lineheight = 1.15))
  ))
}

# ⚠⚠ 题注必须**量宽 + 量高**折行，且必须按「最宽的一端」定稿：
#   本机与服务器字形宽度不同（实测服务器约宽 12%）→ 同一行在本机占画布 96%、
#   在服务器**直接被右边缘裁掉**（本图首版就是如此：第 2 行 268 字符，服务器上断在
#   "the absolute hit cou"）。**光看本机渲图"没超"不算通过**。
#   硬上限：每行 ≤ 170 字符；并在下面用 grobWidth/grobHeight 实测自检。
ss_note <- if (ss_mode == "fill") {
  paste0("Left: the same filter outcome per identification track, normalised to 100% ",
         "(hit totals differ 750-fold); the absolute hit count is in each set label.")
} else {
  paste0("Left: absolute hit counts per identification track ",
         "(viralverify 452,885 vs viroid_blastn 604, a 750-fold span); ",
         "right: the same filter outcome per intersection.")
}
sub_txt <- paste0(
  "Lycium barbarum virome: which identification tool(s) found each contig, and what survived the filter chain. ",
  "Both bar panels are stacked by filter outcome.\n",
  "395 samples; union ", format(n_union, big.mark = ","),
  " viral candidates -> CDD ", format(n_cdd, big.mark = ","),
  " kept -> UniProt ", format(n_ret, big.mark = ","),
  " kept (", sprintf("%.1f", 100 * n_ret / n_union), "% of union); top ", n_show,
  " of ", n_all, " intersections shown (",
  sprintf("%.2f", 100 * cum_shown / sum(sig$n_total)), "% of sequences).\n",
  ss_note, "\n",
  "viroid_blastn is a separate viroid track: 88 of its 604 hits enter the filter chain, the other 516 do not. ",
  "27 contigs from 5 samples lacking 02b output are excluded."
)

W_FIG   <- 11.4
H_FIG   <- 7.4
H_TITLE <- 0.22                      # 标题带权重（plot_layout 的 heights）
Y_TITLE <- 0.97                      # 标题顶对齐锚点（npc of 标题带）
Y_SUB   <- 0.80                      # 副标题顶对齐锚点
# ⚠⚠ plot_layout 的 heights 是**相对权重** → 标题带实际高度 = H_FIG × H_TITLE/(H_TITLE+1)，
#   不是 H_FIG × H_TITLE。按后者算自检会**偏大 18%** → 假通过（首版就是这么漏掉的）。
BAND_IN <- H_FIG * H_TITLE / (H_TITLE + 1)
cap_gp  <- gpar(fontsize = 6.4, lineheight = 1.15)
cap_w   <- convertWidth(grobWidth(textGrob(sub_txt, gp = cap_gp)), "in", valueOnly = TRUE)
cap_h   <- convertHeight(grobHeight(textGrob(sub_txt, gp = cap_gp)), "in", valueOnly = TRUE)
avail_h <- BAND_IN * Y_SUB - 0.04
cat(sprintf("标题带实测 %.2f in｜题注宽 %.2f in（可用 %.2f in）｜题注高 %.2f in（可用 %.2f in）\n",
            BAND_IN, cap_w, W_FIG - 0.30, cap_h, avail_h))
if (cap_w > W_FIG - 0.30) {
  stop(sprintf("题注过宽（%.2f in > %.2f in）→ 会被右边缘裁掉，请再折行",
               cap_w, W_FIG - 0.30))
}
if (cap_h > avail_h) {
  stop(sprintf("题注过高（%.2f in > %.2f in）→ 会被标题带下边缘裁掉，请减少行数或加高标题带",
               cap_h, avail_h))
}

# ═══════════ 图例 → 左上角空白格 ═══════════
# 用户 2026-10-13 追加：「图例放到左上角呢」。
# ⚠ `upset()` 返回的是 **patchwork**，grid 布局 = 2 列 × 2 行：
#     widths  = c(0.3, 0.7)  → [set-size 列 | 主列]
#     heights = c(1, 0.5)    → [交集柱行 | 矩阵行]
#   左上格 (行1, 列1) 是一个 **`spacer`** —— 整块空白（就是主图上那片大白区）→ 正好放图例。
#   ⚠⚠ `widths` 是**显式**值 → 该列宽固定，放进去的对象**超宽会被裁掉且不报错** → 必须量宽自检。
spacer_idx <- which(vapply(p_body$patches$plots, function(x) inherits(x, "spacer"), logical(1)))
if (length(spacer_idx) != 1) {
  stop(sprintf("预期左上角恰好 1 个 spacer 格子，实际 %d 个 → ComplexUpset 版式变了，请复查",
               length(spacer_idx)))
}
# ═══════════ 图例 → Venn 右侧一列（用户 2026-10-13：先「放右上角」，再改为「一列右侧」）═══
# ⚠⚠ 历史：图例原先在**左上角**，塞进 ComplexUpset patchwork 的 (行1,列1) `spacer` 格。
#   第 7 次修订把**耗时条形图**放进了 left panel 下方的空白带（见下），该带的上半部分
#   （y 280..470 附近）与图例（原占 x 67..376 / y 474..592）**同一水平区** →
#   两者挤在左上半区；第 8 次修订遂把图例挪到"右上角"（x2560..3320 / y435..640）。
#   ⚠ 但该带**仍与 Venn 顶部交叠**（Venn 墨迹 x 到 2856、y 从 406 起）→ 用户指出
#   「图例挡到 venn 了，一列右侧」→ 第 9 次修订改为**竖排一列、放 Venn 右侧空档**。
#
# ⭐ 实现选择（为什么不换格子、而用「预渲 PNG + rasterGrob」）：
#   ⚠ 右上角那个"空白区"**不是** patchwork 的格子 —— 它是交集柱面板（x1490..3362 的
#     ggplot 面板）内部的无人区，`p_body$patches$plots[...]` 换不到它（那格是柱面板本体，
#     非 spacer）。而可见的 patchwork 格子只有 4 个（2×2），全部已被占用。
#   → 故沿用文件内既有的、已验证「零 reflow」的通道：**右上角空白区 = 第三个 overlay**，
#     与 Venn/耗时图同一套 `save_stacked_all(overlays=)` 机制。
#
# ⚠⚠⚠ **第 9 次修订（用户 2026-10-13：「图例挡到 venn 了，一列右侧」）**：
#   第 8 次修订把图例放在 Venn **上方**（x2560..3320 / y435..640），但实测该带
#   与 Venn 顶部**仍有交叠** —— Venn 墨迹 y 从 406 起、x 到 2856，
#   而图例框从 x2560 起 → 图例左侧 300 px 横跨在 Venn 上（靠"图例在上"的叠加顺序
#   才勉强可见，视觉上仍是**压住 Venn 顶部两条标签**）。
#   → 用户要求：**挪到 Venn 右侧，排成一列**。
#   ⭐ 实测干净基线（`--no-legend-ov --no-timing`，3420×2220）右侧空档：
#     Venn 墨迹右缘 x=2858；面板右缘 x≈3362；**x2890..3350 从 y0 一直到 y1470 全空**
#     （该列唯一的墨迹从 y1478 起，是下方 x 轴条，与本框无关）。
#     → 可用列 = **x 2890..3350（460 px）、y 430..1290（860 px）= 1.53 × 2.87 in**，
#       与 Venn 纵向同高，视觉上左右并排。
#   ⚠ 图例本体只有 4 项 → 竖排自然尺寸 ~1.06 × 0.47 in，远小于该框；
#     **不把框缩小**（否则 overlay 框小、留白落在框外仍可能压到东西），
#     而是**框取满该列、图例在其中垂直居中** → 视觉上"贴在 Venn 右侧一列"。
LEG_BOX <- list(l = 2890, r = 3350, t = 430, b = 1290)   # 画布 px（3420×2220）
LEG_L <- LEG_BOX$l / 3420; LEG_R <- LEG_BOX$r / 3420
LEG_T <- 1 - LEG_BOX$t / 2220; LEG_B <- 1 - LEG_BOX$b / 2220
leg_w_in <- (LEG_R - LEG_L) * W_FIG
leg_h_in <- (LEG_T - LEG_B) * H_FIG

# ⚠ 图例本体用一个"无图层"的 ggplot + **竖排一列（4 行）**：
#   第 9 次修订起改为 **nrow = 4 / direction = "vertical"**（用户「一列」）。
#   ⚠⚠ 实测竖排自然尺寸（gtable guide-box 量法，见下）：
#        key=0.30cm text=7.0 → **1.06 × 0.47 in**（4 项 1 列的实测值）。
#   ⚠ 为什么不用 nrow=1 默认：`guide_legend` 默认 `nrow = NULL` → 按 direction 自动；
#     显式给 `nrow = 4` 才能保证**一列 4 行**（不受 legend.box 排列影响）。
legend_plot <- ggplot(
    data.frame(status = factor(STATUS_LEVELS, levels = STATUS_LEVELS)),
    aes(x = 1, y = 1, fill = status)) +
  # ⚠ 必须留一个 layer：ggplot 无图层会直接报 "no layers in plot"。
  #   用 alpha = 0 让它画不出东西，再在 guide 里 override.aes 把**图例色块**改回不透明
  #   （不 override 的话图例也跟着透明 → 图例整块看不见，实测过）。
  geom_col(alpha = 0) +
  scale_fill_manual(values = STATUS_COLORS, name = NULL, drop = FALSE,
                    guide = guide_legend(override.aes = list(alpha = 1), nrow = 4)) +
  theme_void(base_size = 8) +
  theme(
    legend.position      = c(0.5, 0.5),   # 居中（框本身就是量好的空区）
    legend.justification = c(0.5, 0.5),
    legend.direction     = "vertical",
    legend.key.size      = unit(0.30, "cm"),
    legend.text          = element_text(size = 7.0),
    legend.margin        = margin(0, 0, 0, 0),
    legend.box.spacing   = unit(0, "cm")
  )
# ⚠⚠ 量宽不能直接 `grobWidth(ggplotGrob(legend_plot))` —— ggplotGrob 的宽高是 **null 单位**，
#   转成 in 恒得 **0.00** → 自检**恒假通过**（首版就是这么写的，实测打印 `0.00 × 0.00 in`）。
#   正确量法：把图例改成**固定位置**（`right`）后从 gtable 里把 `guide-box` 抠出来量；
#   尺寸与数值定位版**完全一致**（同一套 key/text/direction 参数，只是摆放方式不同）。
#   ⚠⚠ 若改用**像素法**（渲 PNG 数墨迹）也可，但注意：`theme_void` 下 `legend.position="right"`
#     会把图例压在画布右中，量出的 h 会含额外留白 → 两种量法务必**同口径**。
lg_grob <- ggplotGrob(legend_plot + theme(legend.position = "right"))
gi <- which(vapply(lg_grob$grobs, function(x) x$name, character(1)) == "guide-box")
if (length(gi) != 1) {
  stop(sprintf("取不到图例 guide-box（命中 %d 个）→ ggplot2 版式变了，请复查", length(gi)))
}
leg_w   <- convertWidth(grobWidth(lg_grob$grobs[[gi]]), "in", valueOnly = TRUE)
leg_h   <- convertHeight(grobHeight(lg_grob$grobs[[gi]]), "in", valueOnly = TRUE)
# ⚠⚠ 判据口径：该框是**独立 overlay**，不受 patchwork `widths` 裁剪 → 可用宽 = 框宽。
#   仍留 10% 内缩余量（服务器字形比本机宽 ~12%）→ 换算到"框宽"上判。
#   ⭐ 第 9 次修订：框变成**竖长**（1.53 × 2.87 in）→ 判据仍是**宽与高都不超框的 90%**。
cat(sprintf("图例实测 %.2f × %.2f in｜右侧列框 %.2f × %.2f in\n",
            leg_w, leg_h, leg_w_in, leg_h_in))
if (leg_w > leg_w_in * 0.90 || leg_h > leg_h_in * 0.90) {
  stop(sprintf("图例 %.2f×%.2f in 超出右侧列框 %.2f×%.2f in（>90%%）→ 会被裁掉，请缩字号/减条目",
               leg_w, leg_h, leg_w_in, leg_h_in))
}
# ⚠⚠ 左侧 spacer 格**仍留占位**：若把该格内容移除，ComplexUpset 的 set-size 列宽可能变化
#   → 引发 reflow（会破坏"零回归"）。故仍放一个 `geom_blank` 的极简占位图。
#   真正的图例改由 overlay 画在 **Venn 右侧一列**（见文件末尾 legend_raster）。
spacer_placeholder <- ggplot() +
  geom_blank(aes(x = 1, y = 1)) +
  theme_void()
p_body$patches$plots[[spacer_idx]] <- spacer_placeholder

# ═══════════ 左侧：7 轨道单样本平均耗时条形图（用户 2026-10-13 第 7 次修订追加）═══════════
# 用户原话：「这一部分空白的，可以补充每个软件的时间消耗的条形图嘛」+「单样本平均耗时含 viroid_blastn」。
#
# ⭐⭐ **第 11 次修订（用户 2026-10-13 21:36，附红框截图「这个位置不行嘛」）**：
#   用户圈出**左上空白块的下方**（原耗时框 y280..870 的**正下方**那片更大的空区）并画箭头指向它
#   → 要求耗时图**下移**到那里。核实两个好处：① 原位置紧邻交集面板 y 轴题带（y893..1057），
#   下移后**彻底分开**；② 新块更宽更高。
#   ⭐ 实测新空块（本版渲染图逐像素扫，`ink=(R<235)|(G<235)|(B<235)`）：
#     **x 150..1250、y 1080..1560 墨迹 = 0 像素**（列带扫描 0..1240 全 0；
#      x1240 起才有 Sequences 题注，x1360 起才是柱面板）。
#     ⚠ 上界必须 ≥1080：Sequences 题注在本 x 区间内占 y893..1057 → 压上去会重叠。
#     ⚠ 下界必须 ≤1560：左侧 set-size 条从 y1572 起。
#   → 取框 **x 180..1250 / y 1080..1560** = 1070 × 480 px = **3.57 × 1.60 in**（比例 2.23，适合 7 根横条）。
#   ⚠ 原框（x400..1245 / y280..870）作废 —— 用户明确否定该位。
# ⚠⚠⚠ **如何把 Venn 叠到「已组装好的 patchwork」上** —— 本轮第 6 次修订的全部踩坑在此：
#   ① `layout + inset_element(...)` ✗ —— `layout` 是 patchwork，`align_to`
#      （`"full"`/`"plot"`/`"panel"` 三值都试过）**不以整画布为基准** → 框 npc 被误解析。
#      标记块实测：目标中心 x=2426 y=722，实际落 **x=2056 y=598**，且"框"报成 2.82×4.16 in
#      （真值 2.95×2.95）→ **零报错、零警告**。真图上表现为 Venn 只有 **1.40 × 1.37 in**
#      （应 2.60 × 2.68 in）→ 上一版「右上角」看着"放得下"其实是被静默缩了一半。
#   ② `wrap_elements(full = layout) + inset_element(...)` ✗ —— **patchwork 二次 reflow**：
#      整图重新分配（左侧 set-size / 右下点阵全跟着变），墨迹从 x182..3215 → **x203..3194**
#      → 不是刚体平移而是**重缩放**，对齐搜索到 ±6 px 边界仍有 9.88 残差。**不可接受**。
#   ③ `wrap_elements(full = p_body) + inset_element(...)`（p_body 层 inset、再拼标题）
#      ✗ —— 玩具探针看着"近乎中性"，真图上左半（x<1490）mean diff 仍有 **12.937**。
#   ④ `patchwork:::patchworkGrob(layout)` + `grid.draw()` + 自己 `pushViewport()`
#      ✓ **零 reflow**（基线/加图后墨迹逐字节一致），**但** child 直接 `grid.draw(ggplotGrob(p))`
#      会锚在 vp 的**左下角**且 `grobWidth(ggplotGrob(p))` 恒 **0.00 in**（null 单位）→
#      `just=` 给什么都无效，实测红块落 x2058..2270 之间，永远到不了 2426。
#   ⑤ ✓✓ **最终采用**：把 Venn **预渲成 PNG（尺寸 = 目标框）**，再 `grid.draw(rasterGrob(...))`
#      塞进居中的 vp。实测红块中心 **2424,721**（目标 2426,723，差 2 px），
#      且全图 `diff > 12` 的像素数 = **0**（最大 diff 3.0 = 抗锯齿噪声）→ 基底逐像素不变。
#      ⚠ 代价 = PDF 里 Venn 是位图 → 预渲按 **600 dpi** 出图（2× 过采样）保证印刷清晰。
#   → 故 `layout` **照常在 Venn 块之前**赋值为「标题带 / p_body + heights」的 patchwork；
#     Venn 不再用 `inset_element()`，而是在 save 阶段用 gtable+grid 叠上去（见文件末尾 draw_*）。
# ═══════════ 中部：6 个病毒工具的 Venn ═══════════
# 用户 2026-10-13 追加：「绘制6个软件的venn图」（先问过「饼图」，用户改口要 Venn）。
#   ⚠ 用户明确两点：①「过滤后的」→ 只画 **Retained**（通过 CDD+UniProt 两道过滤）的序列，
#     不是总命中；②「只算 6 个病毒工具」→ **不含** viroid_blastn（它是独立轨道）。
#   过滤后口径（实测）：并集 94,807；六工具分别命中 69,817 / 30,414 / 26,811 / 6,918 / 2,682 / 1,514；
#   被 k 个工具命中 = 65,760 / 18,343 / 8,578 / 1,179 / 422 / 525。
#   （对比总命中口径：并集 556,862，k=1 占 78.5%，六交集仅 533 → 6 椭圆 Venn 几乎全空、
#     数字放不下。用户选「过滤后」正是为了避开这个退化形状。）
#
# ⚠⚠ 位置：**不能用 `p_body$patches$plots[...]` 换格子** —— 那块格子里是交集柱面板本身（非空）。
#   ⚠⚠ 用户 2026-10-13 更正：「右上角」是口误，实际要 **放中间**（见下方第 6 次修订测绘）。
#   历史记录（第 5 次修订的「右上角」空白区）：x 1495..3362、y 约 400..660 —— 横长竖短，
#   近正方 Venn 塞进去会被按高缩到 ~0.87 in 宽，数字完全不可读（首版实测即此）。
#
# ⚠⚠ 六椭圆真 Venn 有 **63 个区域**，本数据大量区域为 0/个位数 → 数字必然挤。
#   实测两个致命点（都已踩）：
#     ① `label_alpha` 是**文字本身的透明度**，不是"标签背景"！写成 `label_alpha = 0`
#        → 全部数字**完全透明、一个都看不见**，但 `ggplot_build()` 里 x/y/label **三列齐全、
#        位置也对** → 又是「零报错的静默丢失」。**必须 `label_alpha = 1`**。
#     ② 不给 `scale_fill_gradient()` 时默认连续色阶把整片染成近黑 → 数字（默认黑）与底色糊在一起。
#       必须显式给白色→绿色的 fill 渐变。
#   `shape_id`：6 集合只有两个形状（`get_shapes()`）—— `601`（多边形/椭圆，经典 Venn 观感）
#   与 `601t`（三角放射）。用户要"椭圆 Venn" → **用 `601`**。
#   ⚠ `ggVennDiagram` **本机原本没有**（服务器 4.4.3 有）→ 已在本机 CRAN 安装（含依赖 venn/admisc/forcats），
#     保证双端同一份代码可跑，不许只在服务器装。
make_tool_venn <- function() {
  if (!requireNamespace("ggVennDiagram", quietly = TRUE)) {
    stop("缺少 ggVennDiagram（本机需 install.packages('ggVennDiagram')）")
  }
  # ⚠ 只取 **Retained** 且 **viroid_blastn 未命中**（signature 第 7 位 = 0）的签名
  vsig <- sig[substr(signature, 7, 7) == "0" & n_retained > 0, .(n = n_retained), by = signature]
  if (!nrow(vsig)) stop("过滤后没有可用于 Venn 的签名")
  # ⚠ 展开成逐序列 id：`mapply` 按 n 复制签名串，再按位取子集
  vids <- unlist(mapply(function(s, n) rep(s, n), vsig$signature, vsig$n, SIMPLIFY = FALSE),
                 use.names = FALSE)
  vsets <- lapply(seq_along(VENN_TRACKS), function(i) which(substr(vids, i, i) == "1"))
  names(vsets) <- VENN_TRACKS
  venn_n <- vapply(vsets, length, numeric(1))
  venn_union <- length(unique(unlist(vsets)))
  cat(sprintf("Venn（过滤后 Retained，6 病毒工具）：并集 %s｜各集合 %s\n",
              format(venn_union, big.mark = ","),
              paste(sprintf("%s %s", VENN_TRACKS, format(venn_n, big.mark = ",")), collapse = " / ")))
  # 自检：Venn 并集必须等于「6 工具任一命中」的 Retained 总数（口径一致性）
  ref_union <- sum(vsig$n) - sum(vsig$n[vsig$signature == "000000"])
  if (venn_union != ref_union) {
    stop(sprintf("Venn 并集(%d) ≠ 6 工具 Retained 并集(%d) → 展开逻辑有误", venn_union, ref_union))
  }
  # ⚠⚠ 区域标签**自己画**，不用 ggVennDiagram 的自动标注：
  #   ① 自动标注会给 **63 个区域全标**（实测其中 5 个 count=0）→ 一堆 "0" 白占空间还挤；
  #   ② 六椭圆有大量**区域小到放不下数字**，自动标注会与椭圆弧线**重叠**（实测 "44687" 被弧压住）。
  #   做法 = `label = "none"` 关掉自动标注，再从 `process_data()` 的 `regionLabel` 里
  #   取出每个区域的**质心 (X,Y) 与 count**，只标 `count > 0` 的 58 个，并做**白描边双绘**：
  #   先用白色粗体画一遍（halo），再用主题绿画一遍 → 压在线上的数字仍然清晰。
  #   ⚠ `ggVennDiagram:::process_data()` 是 **S4 泛型**，必须传 `Venn(x)` 对象（传裸 list 会
  #     "unable to find an inherited method"）。
  rdat <- as.data.frame(
    ggVennDiagram:::process_data(ggVennDiagram::Venn(vsets), shape_id = "601")$regionLabel)
  rlab <- rdat[rdat$count > 0, c("X", "Y", "count")]
  cat(sprintf("Venn 区域：%d 个（其中 %d 个 count>0 → 只标这些）\n", nrow(rdat), nrow(rlab)))
  vfmt <- function(v) format(v, big.mark = ",", scientific = FALSE)
  ggVennDiagram::ggVennDiagram(
      vsets, shape_id = "601",
      label = "none",
      set_color    = "#333333",
      set_size     = 1.9,
      edge_size    = 0.45) +
    # ⚠ 必须显式给渐变：默认色阶会把整片染黑（见上文）。
    scale_fill_gradient(low = "white", high = "#009E73", guide = "none") +
    # ① 白描边（halo）：先画一遍白色粗体 —— 数字压在椭圆弧上也能读。
    geom_text(data = rlab, aes(x = X, y = Y, label = vfmt(count)),
              size = 1.75, colour = "white", fontface = "bold") +
    # ② 正色：与右侧柱顶 % / 左侧保留率同色
    geom_text(data = rlab, aes(x = X, y = Y, label = vfmt(count)),
              size = 1.75, colour = "#1A5C3A") +
    theme_void(base_size = 8) +
    theme(plot.margin = margin(0, 0, 0, 0))
}
# ── 组装 layout（Venn 一律在此之后用 gtable+grid 叠，不再用 inset_element）──
layout <- wrap_elements(title_grob(
  "Identification x filtering of L. barbarum virome candidates", sub_txt, Y_TITLE, Y_SUB
)) / p_body +
  plot_layout(heights = c(H_TITLE, 1))

# ═══════════ Venn 叠加：预渲 PNG + rasterGrob（唯一零 reflow 且定位正确的方案）═══════════
venn_raster <- NULL            # 非 NULL 时才叠
if (want_venn) {
  venn_plot <- make_tool_venn()

  # ⚠⚠ 位置：**不能用 `p_body$patches$plots[...]` 换格子** —— 那块格子里是交集柱面板本身（非空）。
  #
  #   ⭐ 用户 2026-10-13 第 6 次修订：「Venn 图放大些 **放中间**」
  #     （随后更正第 5 次的「右上角」为口误，真意就是中间）。
  #     用 `--no-venn` 出的干净图逐像素扫描中部空区：
  #       · 左界 = 最高的第一根柱（x 1370..1481；x≥1490 起全空）；
  #       · 右界 = 面板右缘 x≈3362；
  #       · 上界 = 面板顶 y≈430（再上是 y 轴刻度「3e+05」标签区）；
  #       · **下界 = 柱顶 % 标签带，从 y≈1330 起**（"13%"/"21%"… 实测 y1330..1490 有内容）
  #         → 下缘必须 ≤ y 1320 才不压标签；
  #       · 内部唯一障碍 = y=670 一条浅灰网格线（235,235,235），Venn 内部为白 → 可**忽略**。
  #     → 可用方块 = x 1490..3362、y 430..1320（**1872 × 890 px = 6.24 × 2.97 in**）→ 高度是瓶颈。
  #     目标：Venn 尽量放大（受 y 向 890 px 限制）**且水平居中于交集柱面板**。
  #     ⚠ 先实测 Venn 的墨迹占比（90% × 93%，见下），再反推框尺寸：
  #       墨迹高 ≤ 890 px → **框高 = 890/0.93 = 957 px = 3.19 in**（墨迹 ≈ 861 px = 2.87 in）；
  #       框宽同取 957 px（保持正方）→ 墨迹宽 ≈ 861 px，左右各留 (1872−861)/2 ≈ 505 px。
  #     ⚠ 框中心：水平 = 面板视觉中心 **x 2426**（左缘 1448 离最左柱 x<1490 尚差 42 px，
  #       但**墨迹**左缘 1996 距柱 506 px，安全）；垂直 = 可用带中心上移 25 px → **y 850**
  #       （墨迹下缘 ~1293，离柱顶 % 标签带 y≥1330 留 **37 px** 富余；上缘 ~406 → 会略高于
  #        面板顶 430 —— 但上方是 y 轴刻度标签的**空白**区，实测无内容，可容）。
  VENN_L <- (2426 - 478.5) / 3420; VENN_R <- (2426 + 478.5) / 3420
  VENN_T <- 1 - (850 - 478.5) / 2220; VENN_B <- 1 - (850 + 478.5) / 2220
  venn_w_in <- (VENN_R - VENN_L) * W_FIG
  venn_h_in <- (VENN_T - VENN_B) * H_FIG
  venn_cx_px <- 2426
  venn_cy_px <- 850
  cat(sprintf("Venn 框：x %.3f..%.3f  y %.3f..%.3f npc｜%.2f × %.2f in｜中心 px (%.0f, %.0f)\n",
              VENN_L, VENN_R, VENN_B, VENN_T, venn_w_in, venn_h_in, venn_cx_px, venn_cy_px))

  # ⚠⚠ 量尺寸**不能**用 `grobWidth(ggplotGrob(p))` —— 全是 null 单位，回读恒 **0.00 in**
  #   （首版实测；与「图例量宽」同一个坑）。**可靠做法 = 渲染到临时 png 数像素**。
  tmpf <- tempfile(fileext = ".png")
  ggsave(tmpf, venn_plot, width = 3, height = 3, dpi = 300, bg = "white")
  vb <- png::readPNG(tmpf)
  vink <- (vb[, , 1] < 0.98) | (vb[, , 2] < 0.98) | (vb[, , 3] < 0.98)
  vr <- range(which(rowSums(vink) > 0)); vc <- range(which(colSums(vink) > 0))
  frac_w <- (diff(vc) + 1) / (3 * 300); frac_h <- (diff(vr) + 1) / (3 * 300)
  # 预渲到目标框尺寸后，墨迹实际占位 = 框 × frac（等比、无额外缩放，因框给的就是 in 尺寸）
  ink_w <- venn_w_in * frac_w; ink_h <- venn_h_in * frac_h
  cat(sprintf("Venn 自然墨迹占比 %.0f%% × %.0f%% → 框内墨迹 %.2f × %.2f in\n",
              100 * frac_w, 100 * frac_h, ink_w, ink_h))
  # ⚠ 自检：墨迹**不得**压住左界柱（x≥1490）、右界（x≤3362）、柱顶 % 标签带（y≤1320）。
  #   上界放宽到 y≥395（面板顶 430 之上是 y 轴刻度标签的空白区，实测无内容可容）。
  ink_l <- venn_cx_px - ink_w * 300 / 2; ink_r <- venn_cx_px + ink_w * 300 / 2
  ink_t <- venn_cy_px - ink_h * 300 / 2; ink_b <- venn_cy_px + ink_h * 300 / 2
  cat(sprintf("Venn 墨迹预估 px：x %.0f..%.0f  y %.0f..%.0f（可用区 x1490..3362 y395..1320）\n",
              ink_l, ink_r, ink_t, ink_b))
  if (ink_l < 1490 || ink_r > 3362 || ink_t < 395 || ink_b > 1320) {
    stop(sprintf("Venn 墨迹越界（x %.0f..%.0f y %.0f..%.0f）→ 会压住最左柱或柱顶 %% 标签带",
                 ink_l, ink_r, ink_t, ink_b))
  }

  # ★★ 最终方案：**预渲成 PNG（尺寸 = 框）→ rasterGrob → grid.draw 进居中 vp**。
  #   ⚠ 为什么不用 `inset_element()`：见文件上文「① ② ③」三条否决理由（静默误解析 / reflow）。
  #   ⚠ 为什么不用 `grid.draw(ggplotGrob(p))` 直接画：child 会锚在 vp **左下角**且
  #     `grobWidth` 恒 0 → `just=` 无效（实测红块落 2058~2270，永远到不了 2426）。
  #   ⚠ 过采样 600 dpi（框 2.95 in → 1772 px），保证 PDF 里位图印刷清晰。
  venn_pre <- tempfile(fileext = ".png")
  ggsave(venn_pre, venn_plot, width = venn_w_in, height = venn_h_in,
         dpi = 600, bg = "white")
  venn_raster <- grid::rasterGrob(
    png::readPNG(venn_pre), interpolate = TRUE,
    width = unit(venn_w_in, "in"), height = unit(venn_h_in, "in"))
} else {
  cat("已跳过 Venn（--no-venn）\n")
}

# ── 保存：PNG/PDF 都走 gtable + grid（overlays 逐层叠 rasterGrob）──
# ⭐ 为什么不用 `inset_element()`：作用于 patchwork 会**静默误解析框 npc**（第 6 次修订实测，
#   见上方硬规则）→ 唯一零 reflow 的画法是「预渲 PNG → rasterGrob → 居中 vp 叠在
#   `patchwork:::patchworkGrob(p)` 的 grid 底座上」。overlays = list(list(grob, cx, cy, w, h), ...)
save_stacked_all <- function(p, outdir, name, width, height, overlays = list()) {
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
  one(file.path(outdir, paste0(name, ".png")),
      function(f) grDevices::png(f, width = width, height = height,
                                 units = "in", res = 300, bg = "white"))
  one(file.path(outdir, paste0(name, ".pdf")),
      function(f) grDevices::pdf(f, width = width, height = height, bg = "white"))
  invisible(file.path(outdir, paste0(name, ".png")))
}

# ── 左侧耗时图（第 13 次修订：修正 Y 偏移）：预渲 PNG（尺寸 = 目标框 px）→ rasterGrob ──
# ⭐⭐ **第 13 次修订（用户 2026-10-13 21:55 红框「再次调整大小，符合红色区域的大小」）**：
#   ⚠⚠⚠ **第 12 次修订的红框→画布映射 Y 偏移错了 ≈110 px**（X 也偏 8..20 px），
#   导致耗时图整体**偏高 121 px**（用户看到"还是不在红框里"）。
#   正解 = 用**矩阵 7 行棕点**做锚点（唯一可靠基准），框取
#   **`x 51..1156 / y 545..1494`（3.6833 × 3.1633 in，比例 1.164）**。
#
#   ⭐⭐ 映射（截图 1161×612；列锚点 = 16 个矩阵点列、行锚点 = 7 个矩阵点行）：
#     `act_x = -41.1 + 2.9618 * shot_x`、`act_y = 462.5 + 2.9545 * shot_y`
#     红框外缘 shot x31..404 / y28..349 → 画布 **x51..1156 / y545..1494**。
#     标定残差（shot→画布，实测）：矩阵 7 行 **0 px**、16 列 **0 px**；
#     左面板首条 +25 px；「Identification track」+1 px。→ 亚像素级可信。
#     ⚠ 验证法：矩阵行**棕点数逐行递减**（12/11/13/8/4/1/1）是 UpSet 点阵的指纹，
#       两侧同现即证明这 7 行是同一组 → 可作锚点。
#   ⚠⚠ **第 12 次修订用的 `act_y = (shot_y + 114.4) / 0.3355` 在矩阵行上偏 −111 px**
#     （shot 383 → 预测 1483 / 实际 1594），**X 也偏 +8..+20 px** → 已弃用，勿复用。
#   ⚠ 截图**裁掉了顶部 341 px**（主标题+副标题未入镜）——第 12 次修订的该判断正确，保留。
#
#   ⭐ 净空实测（画布 3420×2220 逐像素扫）：
#     · 新占带 `x51..1156 / y1382..1494`（下移后新增）= **0 墨迹** → 可安全移入
#     · 让出带 `x51..1156 / y424..545` = 30,661 px（现图顶部，移走后腾空）
#     · 下方左面板上缘 `y≈1555`（`y1495..1555` = 0 墨迹）→ 新框底 1494 与面板留 **61 px**
#     · 右缘 1156 <「Sequences」题注左缘 1198（`SEQ_CX_PX`）→ 不重叠 ✓
#
#   ⚠⚠⚠ **框比 1.162 → 1.164（宽 −7 px、高 −8 px）→ 预渲 PNG 必须重渲**
#     （`--w 3.6833 --h 3.1633`），否则 `rasterGrob` 按框硬拉伸，**零报错**（README 硬规则 22）。
TIM_L <- 51 / 3420;   TIM_R <- 1156 / 3420
TIM_T <- 1 - 545 / 2220; TIM_B <- 1 - 1494 / 2220
tim_w_in <- (TIM_R - TIM_L) * W_FIG
tim_h_in <- (TIM_T - TIM_B) * H_FIG
timing_raster <- NULL
if (want_timing) {
  tim_png <- file.path(datadir, "tool_timing_fig.png")
  if (!file.exists(tim_png)) stop("缺耗时图预渲文件，请先跑 make_tool_timing_fig.R")
  timing_raster <- grid::rasterGrob(png::readPNG(tim_png), interpolate = TRUE,
    width = unit(tim_w_in, "in"), height = unit(tim_h_in, "in"))
  cat(sprintf("耗时图框：x %.3f..%.3f y %.3f..%.3f npc｜%.2f × %.2f in（预渲 PNG）\n",
              TIM_L, TIM_R, TIM_B, TIM_T, tim_w_in, tim_h_in))
}

# ── 图例（Venn 右侧一列）：预渲 PNG → rasterGrob（与耗时图 / Venn 同一套零 reflow 机制）──
# ⚠⚠ 为什么不直接把 `legend_plot` 画上去：`grid.draw(ggplotGrob(p))` 会锚在 vp 左下角
#   且 `grobWidth` 恒 0（null 单位）→ `just=` 无效（第 6 次修订实测，见上文 ④）。
#   → 一律走「预渲 PNG（尺寸 = 框）→ rasterGrob」这条已验证通道。
# ⚠ 预渲尺寸用**框尺寸**（不是内容自然尺寸）→ 图例在框内居中，四周留白由框吸收。
legend_raster <- NULL
if (want_legend_ov) {
  lg_pre <- file.path(datadir, "legend_ov.png")
  # ⚠⚠ 与其他预渲资产同理：**参数变了必须重渲**（预渲 PNG 不随参数更新）→
  #   这里每跑必渲（图例很便宜），避免"改了字号却看到旧图"的静默过期。
  ggsave(lg_pre, legend_plot, width = leg_w_in, height = leg_h_in,
         dpi = 600, bg = "white")
  legend_raster <- grid::rasterGrob(png::readPNG(lg_pre), interpolate = TRUE,
    width = unit(leg_w_in, "in"), height = unit(leg_h_in, "in"))
  cat(sprintf("图例框：x %.3f..%.3f y %.3f..%.3f npc｜%.2f × %.2f in（预渲 PNG）\n",
              LEG_L, LEG_R, LEG_B, LEG_T, leg_w_in, leg_h_in))
}

# ── 交集面板 y 轴题「Sequences」：预渲 PNG → rasterGrob ──
# ⚠⚠⚠ **第 10 次修订（用户 2026-10-13 21:03：「upset 的 Y 轴 Sequences 不显示」）**：
#   原生 ylab 在本 patchwork 里脱轴 386 px（详见上文 `labs(y = NULL)` 处注释）→ 改由本 overlay 画。
#   ⭐ 目标位实测（本版渲染图逐像素扫）：y 轴刻度文字左缘 x≈1257；
#      紧贴其左侧的 **x1190..1250 竖带在 y249..1582 全空** → 题注位取 **x≈1220**；
#      纵向中心 = 交集柱面板主区中心（柱顶 ≈y480、基线 ≈y1470）→ **y≈975**。
#   ⚠ 旋转 90°（自下而上）；字号沿用轴题规格 8 pt（与刻度 7.5 pt 协调）。
SEQ_TXT    <- "Sequences"
SEQ_CX_PX  <- 1220                     # 题注中心 x（画布 px）
SEQ_CY_PX  <- 975                      # 题注中心 y（画布 px）
SEQ_W_PX   <- 44                       # 框宽（px）：比 8 pt 字高（≈33 px @300dpi）留 1.33× 余量
SEQ_H_PX   <- 300                      # 框高（px）：纵向足够容纳旋转后整串字
SEQ_L <- (SEQ_CX_PX - SEQ_W_PX / 2) / 3420; SEQ_R <- (SEQ_CX_PX + SEQ_W_PX / 2) / 3420
SEQ_T <- 1 - (SEQ_CY_PX - SEQ_H_PX / 2) / 2220; SEQ_B <- 1 - (SEQ_CY_PX + SEQ_H_PX / 2) / 2220
seq_w_in <- (SEQ_R - SEQ_L) * W_FIG
seq_h_in <- (SEQ_T - SEQ_B) * H_FIG
# ⚠ 预渲：白底 → 与 Venn/图例同为不透明 raster。**顺序无重叠风险**（题注位在 x1190..1250，
#   耗时图右缘 x1245、Venn 左缘 x1947、图例左缘 x2890，三者均不与之横向重叠）+ 纵向 y249..1582 全空。
seq_title_raster <- NULL
{
  seq_pre <- file.path(datadir, "seq_title_ov.png")
  # 每跑必渲（同图例）：避免"改了字号却看到旧图"的静默过期
  ggsave(seq_pre,
         ggplot2::ggplot() +
           ggplot2::annotate("text", x = 0, y = 0, label = SEQ_TXT,
                             angle = 90, size = 8 / .pt, fontface = "plain") +
           ggplot2::xlim(-1, 1) + ggplot2::ylim(-1, 1) +
           ggplot2::theme_void(),
         width = seq_w_in, height = seq_h_in, dpi = 600, bg = "white")
  seq_title_raster <- grid::rasterGrob(png::readPNG(seq_pre), interpolate = TRUE,
    width = unit(seq_w_in, "in"), height = unit(seq_h_in, "in"))
  cat(sprintf("Y 轴题框：x %.3f..%.3f y %.3f..%.3f npc｜%.2f × %.2f in（px %d,%d）\n",
              SEQ_L, SEQ_R, SEQ_B, SEQ_T, seq_w_in, seq_h_in, SEQ_CX_PX, SEQ_CY_PX))
}

overlays <- list()
if (!is.null(timing_raster))
  overlays[[length(overlays) + 1]] <- list(grob = timing_raster,
    cx = (TIM_L + TIM_R) / 2, cy = (TIM_B + TIM_T) / 2, w = TIM_R - TIM_L, h = TIM_T - TIM_B)
if (!is.null(seq_title_raster))
  overlays[[length(overlays) + 1]] <- list(grob = seq_title_raster,
    cx = (SEQ_L + SEQ_R) / 2, cy = (SEQ_B + SEQ_T) / 2, w = SEQ_R - SEQ_L, h = SEQ_T - SEQ_B)
if (!is.null(venn_raster))
  overlays[[length(overlays) + 1]] <- list(grob = venn_raster,
    cx = (VENN_L + VENN_R) / 2, cy = (VENN_B + VENN_T) / 2,
    w = VENN_R - VENN_L, h = VENN_T - VENN_B)
# ⚠⚠⚠ **图例必须排在 Venn 之后**（第 8 次修订实测踩坑，第 9 次修订后仍保留此顺序）：
#   第 8 次修订时图例框 (x2560..3320) 与 Venn 框 (x1947..2905) **横向重叠** x2560..2905（245 px）。
#   Venn 的预渲 PNG 是**白底**（`bg = "white"`）→ 若先画图例、后画 Venn，
#   Venn 的白色矩形会把图例左半（"Retained"/"Removed by UniProt" 两行）**整块盖掉**，
#   实测图面只剩右半两条 —— 且**零报错**（看起来像"图例被裁"，实际是被覆盖）。
#   ⭐ 第 9 次修订把图例挪到 **x2890..3350**（Venn 右缘 2858 之右，**不再横向重叠**），
#     但**顺序仍保持"Venn 在下、图例在上"** —— 这是零成本保险：Venn 若将来再放大，
#     顺序对了才不会静默盖掉图例。⚠ 另：两框在纵向上有重叠（y430..1290 vs Venn y406..1294），
#     所以顺序仍是语义的一部分，不可随意调换。
#   ⚠ 副产品：若将来 Venn 改成**透明底** PNG，则顺序不再敏感 —— 但现在**不可依赖**。
if (!is.null(legend_raster))
  overlays[[length(overlays) + 1]] <- list(grob = legend_raster,
    cx = (LEG_L + LEG_R) / 2, cy = (LEG_B + LEG_T) / 2, w = LEG_R - LEG_L, h = LEG_T - LEG_B)
save_stacked_all(layout, outdir, paste0("Fig_GojiVirome_IdentFilterStackedAll", suffix),
                 W_FIG, H_FIG, overlays)
cat(sprintf("Written: %s\n",
            file.path(outdir, paste0("Fig_GojiVirome_IdentFilterStackedAll", suffix, ".png"))))
