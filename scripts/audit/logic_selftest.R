#!/usr/bin/env Rscript
# ==============================================================================
# logic_selftest.R —— v6.4 分类共识规则层自测（人造数据 + 手算期望值）
#
# 目的：不依赖真实数据集，把 build_cascade_winners / species_conflict_idx /
#       enforce_rank_containment / compute_tool_weights / build_consensus 拆到
#       单函数层，用人造数据逐一验证「文档声称的规则」是否就是「代码实际行为」。
#       期望值全部由人手工推算写下，不从代码里反推。
#
# 用法：Rscript logic_selftest.R [补丁文件路径]
# 输出：PASS/FAIL 表 + 观测项（INFO）+ 退出码（有 FAIL 则非 0）
# ==============================================================================
suppressMessages(library(data.table))

SRC <- if (length(commandArgs(TRUE)) >= 1) commandArgs(TRUE)[1] else
  Sys.getenv("SELFTEST_R", unset = "/tmp/virus_classifier_analysis.R.cascade_v64")

# ---------- 0. 从补丁里只提取「函数定义」与「常量定义」，不执行主流程 ----------
WANT_CONST <- c("TAX_LEVELS", "TAX_V2_COLS", "KNOWN_REALMS", "SUBRANK_SUFFIX", "PLACEHOLDER_TAXA",
                "RANK_DEPTH_WEIGHTS", "TOOL_BIAS", "GENUS_FAMILY_REF", "TAX_GATE_VERSION",
                "CONSENSUS_MODE", "CASCADE_MIN_SHARE", "RATE_FLOOR", "TIE_BREAK_ORDER")
exprs <- parse(SRC)
n_fun <- 0L; got_const <- character(0)
for (e in exprs) {
  if (!is.call(e)) next
  op <- as.character(e[[1]])[1]
  if (!op %in% c("<-", "=")) next
  lhs <- e[[2]]; if (!is.name(lhs)) next
  nm <- as.character(lhs); rhs <- e[[3]]
  is_fun <- is.call(rhs) && identical(as.character(rhs[[1]])[1], "function")
  if (is_fun) { eval(e, envir = globalenv()); n_fun <- n_fun + 1L }
  else if (nm %in% WANT_CONST) { eval(e, envir = globalenv()); got_const <- c(got_const, nm) }
}
# v6.5/v6.6 补：环境变量派生的常量不在提取范围（Sys.getenv 不是 function 定义）
if (!exists("VOTE_WEIGHT_MODE")) {
  .vw <- Sys.getenv("MMPV_VOTE_WEIGHT", unset = "count")
  VOTE_WEIGHT_MODE <- if (.vw %in% c("count", "weighted")) .vw else "count"
}
if (!exists("TIE_BREAK_ORDER")) {
  TIE_BREAK_ORDER <- c("ACVirus", "CAT", "VITAP", "diamond_lca", "genomad", "metabuli", "mmseqs")
}
cat(sprintf("[环境] VOTE_WEIGHT_MODE=%s | TIE_BREAK_ORDER=%s | TAX_GATE_VERSION=%s\n",
            VOTE_WEIGHT_MODE, paste(TIE_BREAK_ORDER, collapse = ">"), TAX_GATE_VERSION))

cat(sprintf("[提取] 补丁=%s\n", SRC))
cat(sprintf("[提取] 函数 %d 个; 常量 %d/%d 缺失=%s\n", n_fun, length(got_const), length(WANT_CONST),
            paste(setdiff(WANT_CONST, got_const), collapse = ",")))
cat(sprintf("[环境] TAX_LEVELS=%s | CASCADE_MIN_SHARE=%s | CONSENSUS_MODE=%s | RATE_FLOOR=%s\n\n",
            paste(TAX_LEVELS, collapse = ">"), CASCADE_MIN_SHARE, CONSENSUS_MODE, RATE_FLOOR))

# ---------- 1. 断言框架 ----------
RES <- list()
chk <- function(case, got, want) {
  ok <- isTRUE(all.equal(got, want, check.attributes = FALSE))
  RES[[length(RES) + 1L]] <<- list(case = case, status = if (ok) "PASS" else "FAIL",
                                   got = paste(format(got), collapse = " | "),
                                   want = paste(format(want), collapse = " | "))
  invisible(ok)
}
inf <- function(case, got) {
  RES[[length(RES) + 1L]] <<- list(case = case, status = "INFO",
                                   got = paste(format(got), collapse = " | "), want = "-")
}

mk <- function(...) rbindlist(lapply(list(...), function(v)
  data.table(contig_id = v[1], Tool = v[2], Rank = v[3], Taxon = v[4], Weight = as.numeric(v[5]))))
win <- function(long, tools = c("A", "B", "C"), share = 0.5) {
  # v6.6：模拟管线调用点，先按 TIE_BREAK_ORDER 显式排定工具序（未登记名按原序追尾）
  tool_order <- order_tools_by_tiebreak(tools)
  w <- build_cascade_winners(long, tool_order, share)
  if (is.null(w) || nrow(w) == 0) return(data.table())
  dcast(w, contig_id ~ Rank, value.var = "Taxon")
}
gv <- function(dt, cid, rank) {
  if (is.null(dt) || nrow(dt) == 0 || !rank %in% names(dt)) return(NA_character_)
  v <- dt[contig_id == cid][[rank]]
  if (length(v) == 0) NA_character_ else as.character(v)
}
tdt <- function(cid, realm, kingdom, phylum, klass, order, family, genus, species)
  data.table(contig_id = cid, Realm = realm, Kingdom = kingdom, Phylum = phylum, Class = klass,
             Order = order, Family = family, Genus = genus, Species = species)
row8 <- function(dt, i = 1) unlist(as.list(dt[i, TAX_LEVELS, with = FALSE]))
w1 <- function(...) {
  l <- lapply(list(...), function(x) setNames(rep(1, length(TAX_LEVELS)), TAX_LEVELS))
  names(l) <- c(...); l
}
mkref <- function(path, tethys, mimi) {
  g <- c("simplexvirus", "tethysvirus", "mimivirus", "orthoflavivirus", "betavirus")
  ncbi <- c("herpesviridae", tethys, mimi, NA, "flaviviridae")
  vmr  <- c("herpesviridae", tethys, mimi, "flaviviridae", "poxviridae")
  fwrite(data.table(Genus = g, NCBI_Family = ncbi, VMR_Family = vmr,
                    NCBI_n = c(3, 2, 2, NA, 2), VMR_n = c(5, 4, 4, 4, 4), Domain = "Viruses"),
         path, sep = "\t")
  path
}
# 参照表 A：tethysvirus 的定型科 = herpesviridae（供 E1/E3/E4/E5/E6）
REFF_A <- mkref(file.path(tempdir(), "refA.tsv"), "herpesviridae", "mimiviridae")
# 参照表 B：tethysvirus 的定型科 = mimiviridae（供 E2 科-种闸门单独触发 / E7）
REFF_B <- mkref(file.path(tempdir(), "refB.tsv"), "mimiviridae", "mimiviridae")
inf("参照表 A/B 已写入", paste(REFF_A, REFF_B))

# ==============================================================================
# A. build_cascade_winners —— 逐级淘汰投票引擎
# ==============================================================================
cat("== A. build_cascade_winners ==\n")

L <- mk(c("c1", "A", "Realm", "Riboviria", 1), c("c1", "B", "Realm", "Riboviria", 1),
        c("c1", "A", "Kingdom", "Orthornavirae", 1), c("c1", "C", "Kingdom", "Shotokuvirae", 1),
        c("c1", "C", "Phylum", "Uroviricota", 1))
W <- win(L)
chk("A1 Realm 2/2 全票", gv(W, "c1", "Realm"), "Riboviria")
# v6.6：平票不再取行序先到者，改为先看下一级谐调度（C 报了 Phylum，A 没有 -> C 胜）
chk("A1 Kingdom 1/2=0.5 平票 -> 下一级一致性定（C 有 Phylum）", gv(W, "c1", "Kingdom"), "Shotokuvirae")
chk("A1 未被淘汰的 C 仍可在细阶元定音", gv(W, "c1", "Phylum"), "Uroviricota")

L <- mk(c("c2", "A", "Genus", "Genusalpha", 1), c("c2", "B", "Genus", "Genusbeta", 1),
        c("c2", "B", "Species", "Genusbeta spx", 1))
W <- win(L)
# v6.6：Genus 平票 1:1，B 的支持工具在 Species 层有值、A 没有 -> B 胜（旧版取行序先到的 A）
chk("A2 0.5 平票 -> 下一级一致性定（B 报了 Species）", gv(W, "c2", "Genus"), "Genusbeta")
chk("A2 0.5 未淘汰 B(其种仍可胜出)", gv(W, "c2", "Species"), "Genusbeta spx")

L <- mk(c("c3", "A", "Genus", "Genusalpha", 3), c("c3", "B", "Genus", "Genusbeta", 1),
        c("c3", "B", "Species", "Genusbeta spx", 1))
W <- win(L)
chk("A3 0.75>0.5 票首", gv(W, "c3", "Genus"), "Genusalpha")
chk("A3 被淘汰工具在细阶元不能复活", gv(W, "c3", "Species"), NA_character_)

L <- mk(c("c4", "A", "Genus", "Genusalpha", 0.9), c("c4", "B", "Genus", "Genusbeta", 0.05),
        c("c4", "C", "Genus", "Genusgamma", 0.05), c("c4", "B", "Species", "Genusbeta spx", 1))
W <- win(L)
chk("A4 1 票 0.9 击败 2 票 0.05+0.05", gv(W, "c4", "Genus"), "Genusalpha")
chk("A4 两个被淘汰者都不复活", gv(W, "c4", "Species"), NA_character_)

L1 <- mk(c("c5", "A", "Genus", "Genusalpha", 1), c("c5", "B", "Genus", "Genusbeta", 1))
L2 <- mk(c("c5", "B", "Genus", "Genusbeta", 1), c("c5", "A", "Genus", "Genusalpha", 1))
chk("A5 平票且下一级无信息 -> 按工具顺序取 A", gv(win(L1), "c5", "Genus"), "Genusalpha")
# v6.6 关鍵性质：平票裁决不再看输入行序，L2 必须与 L1 同结果（旧版会翻转成 Genusbeta）
chk("A5 换行序不翻转（顺序由 TIE_BREAK_ORDER 固定）", gv(win(L2), "c5", "Genus"), "Genusalpha")
# A5b：平票兜底按 TIE_BREAK_ORDER，而非输入行序（mmseqs 行写在前面也赢不了 ACVirus）
L3 <- mk(c("c5b", "mmseqs", "Genus", "Genusm", 1), c("c5b", "ACVirus", "Genus", "Genusa", 1))
chk("A5b 平票按 TIE_BREAK_ORDER（ACVirus 先于 mmseqs）", gv(win(L3, c("mmseqs", "ACVirus")), "c5b", "Genus"), "Genusa")
chk("A5e order_tools_by_tiebreak：登记名按 TIE_BREAK_ORDER 排，未登记名按原序追尾",
    order_tools_by_tiebreak(c("mmseqs", "vcontact3", "ACVirus", "contigtax")),
    c("ACVirus", "mmseqs", "vcontact3", "contigtax"))
chk("A5e 已是正序时幂等", order_tools_by_tiebreak(c("ACVirus", "mmseqs")), c("ACVirus", "mmseqs"))
# A5c：第一级优先于第二级——ACVirus 排序靠前但它的支持工具在下一级无言，mmseqs 有 -> mmseqs 胜
L4 <- mk(c("c5c", "ACVirus", "Genus", "Genusa", 1), c("c5c", "mmseqs", "Genus", "Genusm", 1),
         c("c5c", "mmseqs", "Species", "Genusm spx", 1))
chk("A5c 下一级一致性优先于工具顺序（mmseqs 报了种）", gv(win(L4, c("ACVirus", "mmseqs")), "c5c", "Genus"), "Genusm")
# A5d：两候选下一级谐调度相同（各 1/1）-> 回到第二级工具顺序，ACVirus 胜
L5 <- mk(c("c5d", "mmseqs", "Genus", "Genusm", 1), c("c5d", "ACVirus", "Genus", "Genusa", 1),
         c("c5d", "ACVirus", "Species", "Genusa spx", 1), c("c5d", "mmseqs", "Species", "Genusm spx", 1))
chk("A5d 两级谐调度相同 -> 回到工具顺序，ACVirus 胜", gv(win(L5, c("mmseqs", "ACVirus")), "c5d", "Genus"), "Genusa")

L <- mk(c("c6", "A", "Order", "Alphaovirales", 3), c("c6", "B", "Order", "Betaovirales", 1),
        c("c6", "C", "Order", "Betaovirales", 1), c("c6", "B", "Genus", "Genusbeta", 1),
        c("c6", "C", "Genus", "Genusgamma", 1))
chk("A6 share=0.6 门槛 0.5 时淘汰", gv(win(L, share = 0.5), "c6", "Genus"), NA_character_)
chk("A6 share=0.6 门槛 0.6 时不淘汰", gv(win(L, share = 0.6), "c6", "Genus"), "Genusbeta")

L <- mk(c("c7", "A", "Family", "Alphaviridae", 3), c("c7", "B", "Family", "Betaviridae", 1),
        c("c7", "B", "Genus", "Genusbeta", 1))
W <- win(L)
chk("A7 有值阶元照常淘汰", gv(W, "c7", "Family"), "Alphaviridae")
chk("A7 无值阶元留空", gv(W, "c7", "Order"), NA_character_)
chk("A7 无值阶元不重置活跃集", gv(W, "c7", "Genus"), NA_character_)

chk("A8 空 long -> 0 行", nrow(win(mk())), 0L)

L <- mk(c("c9", "A", "Genus", "Genusalpha", 0), c("c9", "B", "Genus", "Genusbeta", 1))
chk("A9 零权重工具不得胜", gv(win(L), "c9", "Genus"), "Genusbeta")

r <- tryCatch({ build_cascade_winners(mk(c("c10", "T1", "Genus", "Gx", 1)), paste0("T", 1:31)); "no-error" },
              error = function(e) "error")
chk("A10 >30 工具应报错而非静默", r, "error")

# ==============================================================================
# B. 相容性判据 —— 科-种闸门（species_conflict_idx）与科-属闸门镜像
# ==============================================================================
cat("== B. species_conflict_idx / enforce_rank_containment ==\n")
res_dt <- data.table(
  contig_id = paste0("r", 1:9),
  Family    = c("Mimiviridae", "Poxviridae", "Mimiviridae", "Flaviviridae", "Mimiviridae",
                NA_character_, "MIMIVIRIDAE", "Mimi viridae", "Poxviridae"),
  Genus     = c("Simplexvirus", "Simplexvirus", "Orthoflavivirus", "Orthoflavivirus", "Spumavirus",
                "Simplexvirus", "Simplexvirus", "Simplexvirus", "Betavirus"),
  Species   = c("Simplexvirus alpha", "Simplexvirus beta", "Orthoflavivirus gamma",
                "Orthoflavivirus delta", "Spumavirus", "Simplexvirus epsilon",
                "Simplexvirus zeta", "Simplexvirus eta", "Betavirus theta"))
ref_dt <- data.table(
  g_l      = c("simplexvirus", "orthoflavivirus", "betavirus"),
  ncbi_fam = c("poxviridae", NA_character_, "flaviviridae"),
  vmr_pad  = c(";herpesviridae;poxviridae;", ";flaviviridae;", ";poxviridae;"),
  ncbi_n   = c(3L, NA_integer_, 2L))
sc <- species_conflict_idx(res_dt, ref_dt)
chk("B1 互斥行下标", sc$idx, c(1L, 7L, 8L))
chk("B1 属名可查行数", sc$n_known, 7L)
chk("B1 单侧判据行数(XOR)", sc$n_side, 2L)

g1 <- enforce_rank_containment(copy(res_dt), REFF_B, blank_species = TRUE)
chk("B2 blank_species=TRUE 互斥行 Genus 置空", g1[1, Genus], NA_character_)
chk("B2 blank_species=TRUE 互斥行 Species 一并置空", g1[1, Species], NA_character_)
chk("B2 单侧行不动(r3)", c(g1[3, Family], g1[3, Genus]), c("Mimiviridae", "Orthoflavivirus"))
chk("B2 单侧行不动(r9)", g1[9, Genus], "Betavirus")
g2 <- enforce_rank_containment(copy(res_dt), REFF_B, blank_species = FALSE)
chk("B3 blank_species=FALSE 只清 Genus",
    c(g2[1, Genus], g2[1, Species]), c(NA_character_, "Simplexvirus alpha"))
inf("B4 科-属闸门处置行", paste(g1[is.na(Genus), contig_id], collapse = ","))
inf("B5 科-种闸门处置行", paste(res_dt[sc$idx, contig_id], collapse = ","))
inf("B6 r2(两侧相容)/r4(VMR 相容) 未被判定",
    paste0("r2.Genus=", g1[2, Genus], " r4.Genus=", g1[4, Genus]))

# ==============================================================================
# C. compute_tool_weights —— rate 下限（v6.4 修的「rate>0 倒挂」）
# ==============================================================================
cat("== C. compute_tool_weights ==\n")
# 本节断言全是 v6.4 的加权口径（RANK_DEPTH_WEIGHTS × bias × max(rate, RATE_FLOOR)），
# v6.5 起默认 count 模式会把这些全部压成 1，因此显式切回 weighted 才有意义；
# count 口径的断言放在本节末尾 C7。
.old_vw <- VOTE_WEIGHT_MODE; VOTE_WEIGHT_MODE <- "weighted"
dl <- list(A = data.table(contig_id = "c1"), B = data.table(contig_id = "c1"), C = data.table(contig_id = "c1"))
st <- list(A = list(Genus = 0, Species = NA_real_, Family = 0.4), B = list(Order = 0))
tw <- compute_tool_weights(dl, st, out_dir = NULL)
# 注意：工具名 A/B/C 都未登记在 TOOL_BIAS，实际取默认 bias 0.8；
# 首轮 5 个 FAIL 的根因就是期望值误按已登记工具的 1.2 推算（got 恰为 want/1.5），代码行为正确。
BIA_DEF <- 0.8
chk("C1 rate=0 罚到下限 0.25", unname(tw$A[["Genus"]]), 64 * BIA_DEF * 0.25)
chk("C1 rate=NA 不罚", unname(tw$A[["Species"]]), 128 * BIA_DEF)
chk("C1 rate=0.4 按 0.4 罚", unname(tw$A[["Family"]]), 32 * BIA_DEF * 0.4)
chk("C2 无 rate 条目保持原权重", unname(tw$A[["Realm"]]), 1 * BIA_DEF)
chk("C3 未登记 rate 的工具不罚", unname(tw$C[["Genus"]]), 64 * BIA_DEF)
tw2 <- compute_tool_weights(list(ZZZ = data.table(contig_id = "c1")), NULL, out_dir = NULL)
chk("C4 未登记工具名默认 bias 0.8", unname(tw2$ZZZ[["Genus"]]), 64 * 0.8)
# 已登记工具名的 bias 必须查得到 TOOL_BIAS 表（若查表失效会全部退化成 0.8）
tw4 <- compute_tool_weights(list(ACVirus = data.table(contig_id = "c1"), VITAP = data.table(contig_id = "c1"),
                                 vcontact3 = data.table(contig_id = "c1"), contigtax = data.table(contig_id = "c1")),
                            NULL, out_dir = NULL)
chk("C6 已登记工具 bias: ACVirus 1.2 / VITAP 1.1 / vcontact3 0.7 / contigtax 0.6",
    c(unname(tw4$ACVirus[["Genus"]]), unname(tw4$VITAP[["Genus"]]),
      unname(tw4$vcontact3[["Genus"]]), unname(tw4$contigtax[["Genus"]])),
    64 * c(1.2, 1.1, 0.7, 0.6))
old_floor <- RATE_FLOOR; RATE_FLOOR <- 0.1
tw3 <- compute_tool_weights(list(A = data.table(contig_id = "c1")), list(A = list(Genus = 0)), NULL)
RATE_FLOOR <- old_floor
chk("C5 RATE_FLOOR 旋钮生效(0.1)", unname(tw3$A[["Genus"]]), 64 * BIA_DEF * 0.1)
# v6.5 新增：count 模式（默认）下权重与阶元、bias、rate、RATE_FLOOR 全部无关，恒为 1
VOTE_WEIGHT_MODE <- "count"
twc <- compute_tool_weights(list(A = data.table(contig_id = "c1"), ACVirus = data.table(contig_id = "c1")),
                            list(A = list(Genus = 0, Species = NA_real_, Family = 0.4), ACVirus = list(Species = 0.1)),
                            out_dir = NULL)
chk("C7 count 模式权重恒 1（生物 bias / rate 均停用）",
    c(unname(twc$A[["Genus"]]), unname(twc$A[["Species"]]), unname(twc$A[["Family"]]), unname(twc$ACVirus[["Species"]])),
    c(1, 1, 1, 1))
VOTE_WEIGHT_MODE <- .old_vw

# ==============================================================================
# D. 值合法性判据
# ==============================================================================
cat("== D. is_valid_value_vec / species_quality_score_vec ==\n")
chk("D1 占位字面量全判无效",
    is_valid_value_vec(c("-", "NA", "na", "N/A", "no rank", "undefined", "unknown", "null",
                         "default", "Unclassified", "", NA)),
    rep(FALSE, 12))
chk("D1 正常值判有效", is_valid_value_vec(c("Riboviria", "Mimiviridae", "Simplexvirus alpha")),
    rep(TRUE, 3))
chk("D2 阶元后缀误置判 0 分", species_quality_score_vec(c("Alphaviridae", "Betavirinae")), c(0, 0))
chk("D2 双名法升权 1.5", species_quality_score_vec("Simplexvirus alpha"), 1.5)
chk("D2 单名降权 0.5", species_quality_score_vec("Simplexvirus"), 0.5)
inf("D3 sp. 写法实际得分（限定词降权与双名法升权叠加）",
    species_quality_score_vec(c("Simplexvirus sp.", "Simplexvirus cf. alpha", "unclassified Simplexvirus")))
chk("D3 unclassified 前缀降权 0.1", species_quality_score_vec("unclassified Simplexvirus"), 0.1)

# ==============================================================================
# E. build_consensus 端到端合成用例
# ==============================================================================
cat("== E. build_consensus 端到端 ==\n")
run_bc <- function(dl, tw, mode, reffile) {
  assign("CONSENSUS_MODE", mode, envir = globalenv())
  assign("GENUS_FAMILY_REF", reffile, envir = globalenv())
  od <- file.path(tempdir(), paste0("bc_", mode)); dir.create(od, showWarnings = FALSE)
  r <- build_consensus(dl, tw, od)
  list(res = r, outdir = od)
}
E1 <- tdt("e1", "Riboviria", "Heunggongvirae", "Peploviricota", "Herviviricetes", "Herpesvirales",
          "Herpesviridae", "Simplexvirus", "Simplexvirus alpha")
r1 <- run_bc(list(A = E1), w1("A"), "cascade", REFF_A)
chk("E1 全自洽 -> 8 阶元原样保留", row8(r1$res), row8(E1))
chk("E1 completeness=8", r1$res[1, completeness], 8L)
chk("E1 confidence=1", r1$res[1, confidence], 1)
chk("E1 Species_agree 前缀 1/1", substr(r1$res[1, Species_agree], 1, 3), "1/1")

E2v <- tdt("e2", "Riboviria", "Heunggongvirae", "Peploviricota", "Herviviricetes", "Herpesvirales",
           "Mimiviridae", "Tethysvirus", "Simplexvirus alpha")
r2 <- run_bc(list(A = E2v), w1("A"), "cascade", REFF_B)
chk("E2 科-种闸门只清 Species", r2$res[1, Species], NA_character_)
chk("E2 Genus 不受科-种闸门影响", r2$res[1, Genus], "Tethysvirus")
chk("E2 Family 不动", r2$res[1, Family], "Mimiviridae")
chk("E2 completeness=7", r2$res[1, completeness], 7L)
chk("E2 Species_agree 前缀 0/1", substr(r2$res[1, Species_agree], 1, 3), "0/1")
inf("E2 科-属闸门放行(属与科相容)但科-种闸门拦下并记录",
    if (file.exists(file.path(r2$outdir, "species_containment_blanked.tsv")))
      nrow(fread(file.path(r2$outdir, "species_containment_blanked.tsv"))) else "no-file")

E3 <- tdt("e3", "Riboviria", "Riboviria", "Peploviricota", "Herviviricetes", "Herpesvirales",
          "Herpesviridae", "Simplexvirus", "Simplexvirus alpha")
r3 <- run_bc(list(A = E3), w1("A"), "cascade", REFF_A)
chk("E3 非 -virae 王国置空", r3$res[1, Kingdom], NA_character_)
chk("E3 其余不动", r3$res[1, Realm], "Riboviria")
chk("E3 completeness=7", r3$res[1, completeness], 7L)

E4 <- tdt("e4", "Viruses", "Heunggongvirae", "Peploviricota", "Herviviricetes", "Herpesvirales",
          "Herpesviridae", "uncultured virus", "environmental samples")
r4 <- run_bc(list(A = E4), w1("A"), "cascade", REFF_A)
chk("E4 Realm 占位(整串)", r4$res[1, Realm], NA_character_)
chk("E4 Genus 占位(限定词子串)", r4$res[1, Genus], NA_character_)
chk("E4 Species 占位(限定词子串)", r4$res[1, Species], NA_character_)
chk("E4 completeness=5", r4$res[1, completeness], 5L)

# E5 设计：粗阶元(Realm..Class)三家报同一值，保证 A 不在粗层被淘汰；淘汰精确发生在 Order 层
#   A 的 Order 权重 3 vs B/C 各 1 -> 3/5=0.6>0.5 淘汰 B/C。
#   B/C 在 Family 的权重故意给 5（A 只有 1）：若实现允许「权重大者无视淘汰」，Family 就会落回
#   Mimiviridae；实测落 Herpesviridae，即证明「淘汰优先于权重」。
#   （首轮 E5 的 3 个 FAIL 是我把 A 的 Phylum 报成 Peploviricota、B/C 报 Uroviricota 而权重同为 1，
#     A 在 Phylum 层 1:2 就已出局，根本走不到 Order 层的权重对比；代码行为正确，期望值错。）
E5a <- tdt("e5", "Riboviria", "Heunggongvirae", "Peploviricota", "Herviviricetes", "Herpesvirales",
           "Herpesviridae", "Tethysvirus", "Tethysvirus alpha")
E5b <- tdt("e5", "Riboviria", "Heunggongvirae", "Peploviricota", "Herviviricetes", "Caudovirales",
           "Mimiviridae", "Mimivirus", "Mimivirus beta")
tw5 <- list(A = setNames(c(1, 1, 1, 1, 3, 1, 1, 1), TAX_LEVELS),
            B = setNames(c(1, 1, 1, 1, 1, 5, 5, 5), TAX_LEVELS),
            C = setNames(c(1, 1, 1, 1, 1, 5, 5, 5), TAX_LEVELS))
r5c <- run_bc(list(A = E5a, B = E5b, C = copy(E5b)), tw5, "cascade", REFF_A)
r5l <- run_bc(list(A = E5a, B = E5b, C = copy(E5b)), tw5, "legacy", REFF_A)
chk("E5 cascade Order 3/5=0.6 淘汰 B/C", r5c$res[1, Order], "Herpesvirales")
chk("E5 cascade Family 不复活 B/C(其权重 5 vs A 1)", r5c$res[1, Family], "Herpesviridae")
chk("E5 cascade Genus/Species 随科一致",
    c(r5c$res[1, Genus], r5c$res[1, Species]), c("Tethysvirus", "Tethysvirus alpha"))
chk("E5 legacy 同输入 Order 同样 A 胜(两家一致)", r5l$res[1, Order], "Herpesvirales")
chk("E5 legacy Family 被 B/C 权重 5+5 压倒", r5l$res[1, Family], "Mimiviridae")
chk("E5 legacy Genus/Species 随权重",
    c(r5l$res[1, Genus], r5l$res[1, Species]), c("Mimivirus", "Mimivirus beta"))
inf("E5 两模式粗阶元",
    paste0("cascade=", r5c$res[1, Phylum], "/", r5c$res[1, Class],
           " legacy=", r5l$res[1, Phylum], "/", r5l$res[1, Class]))

E6b <- copy(E1); E6b[, Species := NA_character_]
r6 <- run_bc(list(A = copy(E1), B = E6b), w1("A", "B"), "cascade", REFF_A)
chk("E6 缺值工具不计入该阶元分母(Genus)", substr(r6$res[1, Genus_agree], 1, 3), "2/2")
chk("E6 Species 分母只算有值工具", substr(r6$res[1, Species_agree], 1, 3), "1/1")
chk("E6 Species 取唯一有值工具", r6$res[1, Species], "Simplexvirus alpha")

E7 <- tdt("e7", "Riboviria", "Heunggongvirae", "Peploviricota", "Herviviricetes", "Herpesvirales",
          "Mimiviridae", "Mimivirus", "Tethysvirus alpha")
r7 <- run_bc(list(A = E7), w1("A"), "cascade", REFF_B)
inf("E7 最终 Genus（由种提属 Tethysvirus 被科-属校准改回 Mimivirus）", r7$res[1, Genus])
inf("E7 最终 Species（未与最终 Genus 复核）", r7$res[1, Species])
inf("E7 Genus_agree ;; Species_agree",
    paste0(r7$res[1, Genus_agree], " ;; ", r7$res[1, Species_agree]))
inf("E7 completeness", r7$res[1, completeness])

# ==============================================================================
# F. 淘汰优先级与门槛边界（单元层，绕过闸门，只测投票引擎）
# ==============================================================================
cat("== F. 淘汰优先级 / 门槛边界 ==\n")
# F1 淘汰发生在 Order，B/C 在 Family 的权重 5+5 远大于 A 的 1，检验「淘汰不复活」是否压过权重
Lf1 <- mk(c("f1", "A", "Order", "Herpesvirales", 3), c("f1", "B", "Order", "Caudovirales", 1),
          c("f1", "C", "Order", "Caudovirales", 1),
          c("f1", "A", "Family", "Herpesviridae", 1), c("f1", "B", "Family", "Mimiviridae", 5),
          c("f1", "C", "Family", "Mimiviridae", 5))
Wf1 <- win(Lf1)
chk("F1 淘汰优先于权重(细阶元 5 不复活)", gv(Wf1, "f1", "Family"), "Herpesviridae")
chk("F1 粗阶元 0.6 定音", gv(Wf1, "f1", "Order"), "Herpesvirales")
# F2 恰好 0.5 不淘汰（门槛是严格大于），细阶元回到权重多数
Lf2 <- mk(c("f2", "A", "Order", "Herpesvirales", 1), c("f2", "B", "Order", "Caudovirales", 1),
          c("f2", "A", "Family", "Herpesviridae", 1), c("f2", "B", "Family", "Mimiviridae", 5))
chk("F2 share 恰 0.5 不淘汰 -> 细阶元权重定音", gv(win(Lf2), "f2", "Family"), "Mimiviridae")
chk("F2 share 恰 0.5 平票 -> 下一级拼谐度亦相同，回到工具顺序取 A", gv(win(Lf2), "f2", "Order"), "Herpesvirales")
# F3 A 在 Order 层以 1:2 被淘汰，B/C 定音；未被淘汰的 C 在细阶元照常定音
Lf3 <- mk(c("f3", "A", "Order", "Alphaovirales", 1), c("f3", "B", "Order", "Betaovirales", 1),
          c("f3", "C", "Order", "Betaovirales", 1), c("f3", "C", "Family", "Gammaviridae", 1))
Wf3 <- win(Lf3)
chk("F3 1:2 时少数 A 被淘汰", gv(Wf3, "f3", "Order"), "Betaovirales")
chk("F3 幸存者 C 细阶元照常定音", gv(Wf3, "f3", "Family"), "Gammaviridae")

# ==============================================================================
# 汇总
# ==============================================================================
cat("\n================ 自测结果 ================\n")
for (x in RES) {
  if (x$status == "INFO") cat(sprintf("[INFO] %-44s %s\n", x$case, x$got))
  else if (x$status == "PASS") cat(sprintf("[PASS] %-44s %s\n", x$case, x$got))
  else cat(sprintf("[FAIL] %-44s got=%s  want=%s\n", x$case, x$got, x$want))
}
n_fail <- sum(vapply(RES, function(x) x$status == "FAIL", logical(1)))
n_pass <- sum(vapply(RES, function(x) x$status == "PASS", logical(1)))
n_info <- sum(vapply(RES, function(x) x$status == "INFO", logical(1)))
cat(sprintf("\n合计 PASS=%d FAIL=%d INFO=%d\n", n_pass, n_fail, n_info))
if (n_fail > 0) quit(status = 1)
