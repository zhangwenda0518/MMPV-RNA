# figure_export.R — 图件四格式统一导出 (SVG / PNG / PDF / PNG@1200dpi)
#
# 口径与 Python 侧 `mmpv_common/figure_export.py` **完全一致**:
#
#     <stem>.svg                    矢量, 文本可编辑
#     <stem>.png                    300 dpi, 白底 (报告内嵌用)
#     <stem>.pdf                    矢量, 出版投稿用
#     _hires_1200dpi/<stem>.png     1200 dpi, 白底 (印刷/审稿人索图用)
#
# `_hires_1200dpi/` 是独立子目录: 报告层按 (目录, 主名) 去重, 同目录多一张
# PNG 会被当成第二张图重复内嵌。
#
# 用法 (推荐, ggplot 对象):
#     source("mmpv_common/figure_export.R")     # 或 source(<repo>/mmpv_common/figure_export.R)
#     save_figure(p, "out/fig1.pdf", width = 7, height = 4)
#     save_figure(p, "out/fig1.png", width = 7, height = 4, dpi = 300)
#     # 两者都会产出四种格式, 主名相同
#
# 用法 (base graphics / lattice, 传绘图函数):
#     save_figure(function() plot(x, y), "out/base.pdf", width = 7, height = 4)
#
# 环境变量开关 (与 Python 侧同名):
#     MMPV_FIG_EXPORT=0            关闭 (回退为只写调用方指定格式)
#     MMPV_FIG_FORMATS=svg,png,pdf,hires
#     MMPV_FIG_PNG_DPI=300
#     MMPV_FIG_HIRES_DPI=1200
#     MMPV_FIG_HIRES_DIR=_hires_1200dpi

.fig_env <- function(name, default) {
  v <- Sys.getenv(name, unset = "")
  if (identical(v, "")) default else v
}

.fig_enabled <- function() {
  !tolower(.fig_env("MMPV_FIG_EXPORT", "1")) %in% c("0", "false", "no", "off")
}

.fig_formats <- function() {
  raw <- .fig_env("MMPV_FIG_FORMATS", "svg,png,pdf,hires")
  want <- trimws(tolower(strsplit(raw, ",", fixed = TRUE)[[1]]))
  intersect(want, c("svg", "png", "pdf", "hires"))
}

.fig_int <- function(name, default) {
  v <- suppressWarnings(as.integer(.fig_env(name, as.character(default))))
  if (is.na(v) || v < 1L) default else v
}

# 任意目标路径 → 指定格式的落盘路径 (与 Python target_path 同规则)
fig_target_path <- function(path, fmt) {
  d <- dirname(path)
  stem <- sub("\\.[^./\\\\]*$", "", basename(path))
  if (!nzchar(stem)) stop("figure_export: 无法从 ", path, " 推断主名")
  if (identical(fmt, "hires")) {
    return(file.path(d, .fig_env("MMPV_FIG_HIRES_DIR", "_hires_1200dpi"),
                     paste0(stem, ".png")))
  }
  file.path(d, paste0(stem, ".", fmt))
}

# 按 figure 尺寸估算像素, 超限则降级 (与 Python 侧 MAX_HIRES_PIXELS 同口径)
.fig_effective_dpi <- function(width, height, dpi) {
  px <- width * dpi * height * dpi
  max_px <- 400e6
  if (is.finite(px) && px > max_px) {
    fallback <- max(1L, as.integer(dpi * sqrt(max_px / px)))
    fallback <- min(fallback, 600L)
    message(sprintf(
      "[figure_export] 警告: %.1fx%.1f in @ %d dpi ≈ %.0f Mpx 超限, 高清 PNG 降级到 %d dpi",
      width, height, dpi, px / 1e6, fallback))
    return(fallback)
  }
  dpi
}

# 把 "写一种格式" 收敛成单点: ggplot 用 print(), 函数用调用, 其他用 print()
.fig_draw <- function(x, dev_open) {
  dev_open()
  on.exit(grDevices::dev.off(), add = TRUE)
  if (is.function(x)) {
    x()
  } else {
    print(x)
  }
}

.fig_write <- function(x, dst, fmt, width, height, dpi, bg) {
  dir.create(dirname(dst), recursive = TRUE, showWarnings = FALSE)

  if (identical(fmt, "svg")) {
    if (!isTRUE(capabilities("cairo"))) {
      stop("SVG 需要 cairo 支持 (capabilities('cairo') 为 FALSE)")
    }
    .fig_draw(x, function() grDevices::svg(dst, width = width, height = height,
                                           bg = bg, onefile = TRUE))
    return(TRUE)
  }
  if (identical(fmt, "pdf")) {
    .fig_draw(x, function() grDevices::pdf(dst, width = width, height = height,
                                           bg = bg, onefile = TRUE))
    return(TRUE)
  }
  if (fmt %in% c("png", "hires")) {
    use_dpi <- if (identical(fmt, "hires")) {
      .fig_effective_dpi(width, height, dpi)
    } else {
      dpi
    }
    .fig_draw(x, function() grDevices::png(dst, width = width, height = height,
                                           units = "in", res = use_dpi,
                                           bg = bg, type = "cairo"))
    return(TRUE)
  }
  stop("未知格式: ", fmt)
}

#' 四格式导出 (ggplot / 绘图函数 / base 图)
#'
#' @param x ggplot 对象、lattice 对象, 或**绘图函数** (base graphics)
#' @param path 任意目标路径; 只取目录 + 主名, 扩展名由 formats 决定
#' @param width,height 英寸
#' @param dpi 普通 PNG 的 dpi (默认 300); 高清 PNG 用 MMPV_FIG_HIRES_DPI (默认 1200)
#' @param bg 底色 (默认 white)
#' @param formats 默认 svg/png/pdf/hires, 可用 MMPV_FIG_FORMATS 覆盖
#' @param overwrite FALSE 时已存在的文件跳过
#' @return 命名列表: 各格式的落盘路径 (失败为 NA)
save_figure <- function(x, path, width = 7, height = 5, dpi = NULL, bg = "white",
                        formats = NULL, overwrite = TRUE) {
  if (is.null(dpi)) dpi <- .fig_int("MMPV_FIG_PNG_DPI", 300L)
  hires_dpi <- .fig_int("MMPV_FIG_HIRES_DPI", 1200L)
  if (!.fig_enabled()) return(invisible(list()))

  fmts <- if (is.null(formats)) .fig_formats() else intersect(formats,
                                                              c("svg", "png", "pdf", "hires"))
  out <- stats::setNames(rep(NA_character_, 4L), c("svg", "png", "pdf", "hires"))

  for (fmt in fmts) {
    dst <- fig_target_path(path, fmt)
    out[[fmt]] <- dst
    if (!overwrite && file.exists(dst)) next
    use_dpi <- if (identical(fmt, "hires")) hires_dpi else dpi
    ok <- tryCatch({
      .fig_write(x, dst, fmt, width, height, use_dpi, bg)
      TRUE
    }, error = function(e) {
      message("[figure_export] 警告: ", fmt, " 写出失败 (", dst, "): ",
              conditionMessage(e))
      FALSE
    })
    if (!ok) out[[fmt]] <- NA_character_
  }
  invisible(out)
}

#' 只补写缺失的兄弟格式 (给已有 ggsave/pdf() 写法的老脚本用)
#'
#' 老脚本常常已经写了 png + pdf, 只需补 svg 与 hires。
#' @param x ggplot 对象或绘图函数
#' @param path 已知的任一路径 (通常是脚本刚写的那个)
#' @param formats 要补的格式, 默认 c("svg","hires")
ensure_figure_formats <- function(x, path, width = 7, height = 5, dpi = NULL,
                                  bg = "white", formats = c("svg", "hires")) {
  if (is.null(dpi)) dpi <- .fig_int("MMPV_FIG_PNG_DPI", 300L)
  if (!.fig_enabled()) return(invisible(list()))
  fmts <- intersect(formats, .fig_formats())
  out <- list()
  for (fmt in fmts) {
    dst <- fig_target_path(path, fmt)
    if (file.exists(dst)) next
    use_dpi <- if (identical(fmt, "hires")) {
      .fig_int("MMPV_FIG_HIRES_DPI", 1200L)
    } else {
      dpi
    }
    tryCatch({
      .fig_write(x, dst, fmt, width, height, use_dpi, bg)
      out[[fmt]] <- dst
    }, error = function(e) {
      message("[figure_export] 警告: ", fmt, " 写出失败 (", dst, "): ",
              conditionMessage(e))
    })
  }
  invisible(out)
}

#' 按固定尺寸重画一遍, 补齐缺失的格式 (给 pdf()/png() 直写的老封装用)
#'
#' 老封装常常这样写: 开 pdf 设备画一次 → 开 png 设备再画一次。改成:
#'
#'     draw <- function() { print(p) }          # 或 base 绘图语句包成的函数
#'     draw_figure_formats(draw, "out/fig.pdf", width = 8, height = 5)
#'
#' 函数会按 formats 逐个开设备、各画一次 (与老写法的重画次数同阶),
#' 已存在的格式直接跳过 —— 所以「脚本自己已经写了 pdf + png」时只补 svg/hires。
#'
#' @param draw 无参函数, 负责把图真正画到**当前设备**上
#' @param path 任意目标路径 (通常传脚本原本写的那一个)
#' @param width,height 英寸
#' @param dpi 普通 PNG 的 dpi; 高清用 MMPV_FIG_HIRES_DPI
#' @param bg 底色
#' @param formats 默认四格式, 可用 MMPV_FIG_FORMATS 覆盖
#' @return 命名列表: 实际写出/已存在的路径
draw_figure_formats <- function(draw, path, width = 7, height = 5, dpi = NULL,
                                bg = "white", formats = NULL) {
  if (is.null(dpi)) dpi <- .fig_int("MMPV_FIG_PNG_DPI", 300L)
  if (!.fig_enabled()) return(invisible(list()))
  fmts <- if (is.null(formats)) .fig_formats() else intersect(
    formats, c("svg", "png", "pdf", "hires"))
  hires_dpi <- .fig_int("MMPV_FIG_HIRES_DPI", 1200L)

  out <- list()
  for (fmt in fmts) {
    dst <- fig_target_path(path, fmt)
    if (file.exists(dst)) {
      out[[fmt]] <- dst
      next
    }
    use_dpi <- if (identical(fmt, "hires")) hires_dpi else dpi
    ok <- tryCatch({
      .fig_write(draw, dst, fmt, width, height, use_dpi, bg)
      TRUE
    }, error = function(e) {
      message("[figure_export] 警告: ", fmt, " 写出失败 (", dst, "): ",
              conditionMessage(e))
      FALSE
    })
    if (ok) out[[fmt]] <- dst
  }
  invisible(out)
}

# ══════════════════════════════════════════════════════════════════════
# ggsave 钩子 —— 零改动覆盖脚本里既有的 ggplot2::ggsave(...) 调用
# ══════════════════════════════════════════════════════════════════════

.fig_state <- new.env(parent = emptyenv())
.fig_state$installed <- FALSE
.fig_state$orig <- NULL
.fig_state$suppress <- 0L

#' 临时目录豁免: inset / 预览等中间文件不四格式化
#'
#' 判据 = 路径落在 tempdir() 下, 或文件名以 "tmp" / "_tmp" 开头。
#' 中间产物只写调用方要的那一种格式, 避免刷出一堆无用文件。
.fig_is_ephemeral <- function(path) {
  p <- tryCatch(normalizePath(path, mustWork = FALSE), error = function(e) path)
  td <- tryCatch(normalizePath(tempdir(), mustWork = FALSE),
                 error = function(e) tempdir())
  if (startsWith(p, td)) return(TRUE)
  grepl("^_?tmp", basename(p), ignore.case = TRUE)
}

#' 安装 ggsave 钩子 (幂等)
#'
#' 在 .GlobalEnv 定义 `ggsave` 覆盖 ggplot2 的版本: 每次调用先按原参数写目标
#' 格式, 再补齐其余三种。已存在的兄弟文件跳过 (脚本里 png + pdf 双写不会重复渲染)。
#' 中间文件 (tempdir 下 / tmp* 开头) 原样透传, 不补格式。
install_ggsave_hook <- function() {
  if (!.fig_enabled()) return(invisible(FALSE))
  if (isTRUE(.fig_state$installed)) return(invisible(TRUE))
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    message("[figure_export] ggplot2 不可用, 跳过 ggsave 钩子")
    return(invisible(FALSE))
  }

  .fig_state$orig <- ggplot2::ggsave
  orig <- .fig_state$orig

  hooked <- function(filename, plot = ggplot2::last_plot(), device = NULL,
                     path = NULL, scale = 1, width = NA, height = NA,
                     units = c("in", "cm", "mm", "px"), dpi = 300, limitsize = TRUE,
                     bg = NULL, ...) {
    # 计算真实落盘路径 (ggsave 的 path + filename 组合)
    full <- if (is.null(path)) filename else file.path(path, filename)

    # ① 中间文件 / 钩子自身补写 → 原样透传
    if (.fig_state$suppress > 0L || .fig_is_ephemeral(full)) {
      return(orig(filename = filename, plot = plot, device = device,
                  path = path, scale = scale, width = width, height = height,
                  units = units, dpi = dpi, limitsize = limitsize, bg = bg, ...))
    }

    res <- orig(filename = filename, plot = plot, device = device,
                path = path, scale = scale, width = width, height = height,
                units = units, dpi = dpi, limitsize = limitsize, bg = bg, ...)

    # ② 补齐兄弟格式
    tryCatch({
      u <- match.arg(units)
      w_in <- if (is.na(width)) 7 else switch(u, `in` = width, cm = width / 2.54,
                                              mm = width / 25.4, px = width / dpi)
      h_in <- if (is.na(height)) 5 else switch(u, `in` = height, cm = height / 2.54,
                                               mm = height / 25.4, px = height / dpi)
      fmts <- .fig_formats()
      target_fmt <- tolower(tools::file_ext(filename))
      for (fmt in fmts) {
        if (identical(fmt, target_fmt) && !identical(fmt, "hires")) next
        dst <- fig_target_path(full, fmt)
        if (file.exists(dst)) next
        use_dpi <- if (identical(fmt, "hires")) {
          .fig_int("MMPV_FIG_HIRES_DPI", 1200L)
        } else {
          dpi
        }
        .fig_state$suppress <- .fig_state$suppress + 1L
        tryCatch(
          .fig_write(plot, dst, fmt, w_in, h_in, use_dpi,
                     if (is.null(bg)) "white" else bg),
          finally = .fig_state$suppress <- .fig_state$suppress - 1L
        )
      }
    }, error = function(e) {
      message("[figure_export] 警告: 补齐兄弟格式失败: ", conditionMessage(e))
    })
    invisible(res)
  }

  assign("ggsave", hooked, envir = .GlobalEnv)
  .fig_state$installed <- TRUE
  invisible(TRUE)
}

#' 卸载 ggsave 钩子 (测试/调试用)
uninstall_ggsave_hook <- function() {
  if (isTRUE(.fig_state$installed)) {
    if (exists("ggsave", envir = .GlobalEnv, inherits = FALSE)) {
      rm("ggsave", envir = .GlobalEnv)
    }
    .fig_state$installed <- FALSE
  }
  invisible(TRUE)
}
