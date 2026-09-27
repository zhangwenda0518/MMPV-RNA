#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""一次性改造脚本：把本地 R 快照里的 repair_family_from_genus 换成 verify_family_genus_consistency。

按标记行边界整段替换，避免手写长文本匹配。
"""
import io
import sys

PATH = "D:/桌面/延伸基因组/MMPV-RNA/scripts/audit/snapshot/virus_classifier_analysis.R"

START_MARK = "# 科-属校准（以属为准修正科）：行内 Genus 的定型科"
END_MARK = "# 逐级相容性约束（科-属）：行内 Genus 主张必须与行内 Family 相容。"

NEW_BLOCK = '''# 参照表载入：属 -> 双参照定型科。返回 data.table(g_l, fam_raw, fam_l, dual, agree, ncbi_raw, vmr_raw)
load_genus_family_ref <- function(ref_file = GENUS_FAMILY_REF) {
  if (is.null(ref_file) || length(ref_file) == 0 || is.na(ref_file[1]) || !file.exists(ref_file[1])) {
    log_msg("INFO", "科-属参照表不可用: %s",
            if (is.null(ref_file) || length(ref_file) == 0) "NULL" else as.character(ref_file[1]))
    return(NULL)
  }
  ref <- tryCatch(fread(ref_file[1], sep = "\\t", quote = "", colClasses = "character", showProgress = FALSE),
                  error = function(e) NULL)
  if (is.null(ref) || nrow(ref) == 0) {
    log_msg("INFO", "科-属参照表为空或读取失败: %s", ref_file[1])
    return(NULL)
  }
  if (!all(c("Genus", "NCBI_Family", "VMR_Family") %in% names(ref))) {
    log_msg("INFO", "科-属参照表缺 Genus/NCBI_Family/VMR_Family 列: %s", ref_file[1])
    return(NULL)
  }
  norm_tax <- function(x) {
    x <- trimws(gsub("[\\"*]", "", as.character(x)))
    x[!nzchar(x)] <- NA_character_
    x[toupper(x) %in% c("NA", "N/A", "-")] <- NA_character_
    x
  }
  ref_dt <- data.table(
    g_l      = tolower(norm_tax(ref$Genus)),
    ncbi_raw = norm_tax(ref$NCBI_Family),
    vmr_raw  = norm_tax(ref$VMR_Family)
  )
  ref_dt[, ncbi_fam := tolower(ncbi_raw)]
  ref_dt[, vmr_fam  := tolower(vmr_raw)]
  ref_dt <- unique(ref_dt[!is.na(g_l) & (!is.na(ncbi_fam) | !is.na(vmr_fam))], by = "g_l")
  ref_dt[, dual := !is.na(ncbi_fam) & !is.na(vmr_fam)]
  ref_dt[, agree := dual & ncbi_fam == vmr_fam]
  ref_dt[, fam_raw := fcase(
    agree, ncbi_raw,
    is.na(ncbi_fam) & !is.na(vmr_fam), vmr_raw,
    !is.na(ncbi_fam) & is.na(vmr_fam), ncbi_raw,
    default = NA_character_
  )]
  ref_dt[, fam_l := tolower(fam_raw)]
  ref_dt[, .(g_l, fam_raw, fam_l, dual, agree, ncbi_raw, vmr_raw)]
}

# 跑完自检：最终表里不允许出现「属的双参照定型科 与 行内 Family 相符的参照一致却互斥」的行。
# 这类行就是"科属打架"的成品形态。期望值 0；大于 0 说明闸门失效，必须 ERROR 暴露而不是静默出货。
# 台账：out_dir/taxonomy_gate_check.tsv（check / n / expected / status）
verify_family_genus_consistency <- function(result_dt, ref_file = GENUS_FAMILY_REF, out_dir = NULL) {
  checks <- data.table(check = character(), n = integer(), expected = character(), status = character())
  add <- function(check, n, expected, status) {
    checks <<- rbind(checks, data.table(check = check, n = as.integer(n),
                                        expected = expected, status = status))
  }
  n_total <- if (is.null(result_dt)) 0L else nrow(result_dt)
  add("rows_total", n_total, "-", "INFO")
  if (n_total == 0 || !all(c("Family", "Genus") %in% names(result_dt))) {
    add("dual_ref_conflict", 0L, "0", "PASS")
    add("single_ref_conflict", 0L, "-", "INFO")
    return(invisible(checks))
  }
  norm_tax <- function(x) {
    x <- trimws(gsub("[\\"*]", "", as.character(x)))
    x[!nzchar(x)] <- NA_character_
    x[toupper(x) %in% c("NA", "N/A", "-")] <- NA_character_
    x
  }
  ref_dt <- load_genus_family_ref(ref_file)
  n_fam <- sum(!is.na(norm_tax(result_dt$Family)))
  n_gen <- sum(!is.na(norm_tax(result_dt$Genus)))
  add("family_filled", n_fam, "-", "INFO")
  add("genus_filled", n_gen, "-", "INFO")
  if (is.null(ref_dt)) {
    add("dual_ref_conflict", 0L, "0", "SKIP_REF_MISSING")
    add("single_ref_conflict", 0L, "-", "SKIP_REF_MISSING")
    if (!is.null(out_dir)) fwrite(checks, file.path(out_dir, "taxonomy_gate_check.tsv"), sep = "\\t")
    log_msg("WARN", "科-属自检: 参照表不可用, 自检降级为 SKIP")
    return(invisible(checks))
  }
  cur <- data.table(fam_l = tolower(norm_tax(result_dt$Family)),
                    g_l   = tolower(norm_tax(result_dt$Genus)))
  cur <- cur[!is.na(fam_l) & !is.na(g_l)]
  cur <- merge(cur, ref_dt[, .(g_l, fam_l_ref = fam_l, dual, agree)], by = "g_l", all.x = TRUE)
  cur <- cur[!is.na(fam_l_ref)]
  n_dual_conf <- cur[dual == TRUE & fam_l != fam_l_ref, .N]
  n_single_conf <- cur[dual == FALSE & fam_l != fam_l_ref, .N]
  n_ok <- cur[dual == TRUE & fam_l == fam_l_ref, .N]
  add("genus_ref_known", nrow(cur), "-", "INFO")
  add("dual_ref_consistent", n_ok, "-", "INFO")
  add("dual_ref_conflict", n_dual_conf, "0", if (n_dual_conf == 0) "PASS" else "FAIL")
  add("single_ref_conflict", n_single_conf, "unfixable_by_design", "INFO")
  if (!is.null(out_dir)) {
    fwrite(checks, file.path(out_dir, "taxonomy_gate_check.tsv"), sep = "\\t")
  }
  if (n_dual_conf > 0) {
    log_msg("ERROR", "科-属自检失败: 双参照一致却与行内 Family 互斥的行 %d 条（应为 0）", n_dual_conf)
  } else {
    log_msg("INFO", "科-属自检通过: 双参照一致冲突 0 条; 单侧冲突 %d 条（设计内不动）; 双参照一致相容 %d 条",
            n_single_conf, n_ok)
  }
  invisible(checks)
}

'''

with io.open(PATH, "r", encoding="utf-8") as fh:
    lines = fh.readlines()

si = ei = None
for i, ln in enumerate(lines):
    if si is None and ln.startswith(START_MARK):
        si = i
    if si is not None and ln.startswith(END_MARK):
        ei = i
        break
if si is None or ei is None:
    print("ERROR: 未找到边界标记 si=%s ei=%s" % (si, ei))
    sys.exit(1)

print("替换行区间: %d..%d (共 %d 行)" % (si + 1, ei, ei - si))
out = lines[:si] + [NEW_BLOCK] + lines[ei:]
with io.open(PATH, "w", encoding="utf-8", newline="") as fh:
    fh.writelines(out)
print("已写回, 原 %d 行 -> 新 %d 行" % (len(lines), len(out)))
