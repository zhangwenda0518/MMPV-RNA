# 实测 9/4 之前的 harmonize_genus_species 写法为何失效 (tstrsplit 未取 [[1]])
suppressMessages(library(data.table))

sp <- c("Guapo partitivirus", "NA", "Biavirus raunefjordenense", "Solo")

cat("=== OLD 写法 (Aug4~Sep4): ifelse + tstrsplit(keep=1L) 不取 [[1]] ===\n")
dt <- data.table(Species = sp)
r <- tryCatch({
  dt[, g := ifelse(grepl(" ", sp, fixed = TRUE), tstrsplit(sp, " ", keep = 1L), NA_character_)]
  cat("  未报错\n")
  cat("  列类型 : "); cat(class(dt$g), "\n")
  cat("  内容   : "); print(dt$g)
  cat("  nrow   : "); cat(nrow(dt), "\n")
  g2 <- dt$g
  cat("  grepl('virus$', g) 结果: "); print(tryCatch(grepl("virus$", g2, ignore.case = TRUE),
                                                     error = function(e) paste("报错:", conditionMessage(e))))
  cat("  is.na(g) 结果        : "); print(tryCatch(is.na(g2),
                                                   error = function(e) paste("报错:", conditionMessage(e))))
  TRUE
}, error = function(e) { cat("  报错 ->", conditionMessage(e), "\n"); FALSE })

cat("\n=== NEW 写法 (现行): fifelse + [[1]] + type.convert=FALSE ===\n")
dt2 <- data.table(Species = sp)
dt2[, g := fifelse(grepl(" ", sp, fixed = TRUE),
                   tstrsplit(sp, " ", keep = 1L, type.convert = FALSE)[[1]], NA_character_)]
cat("  列类型 : "); cat(class(dt2$g), "\n")
cat("  内容   : "); print(dt2$g)
