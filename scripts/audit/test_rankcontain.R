## Independent test for enforce_rank_containment() in virus_classifier_analysis.R
## Strategy: parse() the whole pipeline script, then eval ONLY the selected top-level
## assignments into a private env, so the main flow never runs. Then:
##   (a) 11 synthetic boundary cases against a synthetic reference table
##   (b) fail-safe cases (missing ref / wrong columns / blank_species=FALSE / empty table)
##   (c) the real 20892-row table against the real 16 MB reference table,
##       required to reproduce exactly the same blank set as apply_calib_A.py wrote.
## Any failed stopifnot() aborts with non-zero status.
suppressPackageStartupMessages(library(data.table))

SRC <- "/home/zhangwenda/MMPV-RNA/virome_discovery_pipeline/virus_classifier_analysis.R"
REF <- "/home/zhangwenda/database/taxonomy/genus_family_ref.tsv"
TBL <- "/home/zhangwenda/MMPV-paper/goji-virome/02_novel_virus/RNA-Lycium_barbarum_out/05_Taxonomy/Votus.integrated/final_integrated_classification.tsv"
CAL <- "/home/zhangwenda/MMPV-paper/goji-virome/02_novel_virus/RNA-Lycium_barbarum_out/05_Taxonomy/Votus.integrated/calibration_20260914/final_integrated_classification.calibrated.tsv"

strip <- function(x) gsub('"', "", x)

## ---------------------------------------------------------------- load function
exprs <- parse(file = SRC, encoding = "UTF-8")
cat(sprintf("[1] parse() OK: %d top-level expressions\n", length(exprs)))
want <- c("GENUS_FAMILY_REF", "log_msg", "is_valid_value_vec", "enforce_rank_containment")
env <- new.env(parent = globalenv())
for (e in exprs) {
  if (!is.call(e)) next
  if (!identical(as.character(e[[1]])[1], "<-")) next
  nm <- as.character(e[[2]])[1]
  if (nm %in% want) eval(e, env)
}
stopifnot(is.function(env$enforce_rank_containment))
cat(sprintf("[2] loaded into private env: %s\n", paste(ls(env), collapse = ", ")))
cat(sprintf("[2] GENUS_FAMILY_REF default = %s\n", env$GENUS_FAMILY_REF))

## ------------------------------------------------------- synthetic reference
TMP <- tempfile(fileext = ".tsv")
writeLines(c(
  "Genus\tNCBI_Family\tNCBI_n\tVMR_Family\tVMR_n\tDomain",
  "gblank\tFamOther\t1\tFamOther\t1\tRiboviria",
  "gok\tFamRow\t1\tFamRow\t1\tRiboviria",
  "gncbionly\tFamOther\t1\t\t0\tRiboviria",
  "gvmronly\t\t0\tFamOther\t1\tRiboviria",
  "gmulti\t\t0\tFamRow;FamOther\t2\tRiboviria",
  "gncbiok\tFamRow\t1\tFamOther\t1\tRiboviria",
  "gambig\tFamOther\t2\tFamOther\t1\tRiboviria"
), TMP)

## c01 double reject                 -> blank
## c02 quoted Family, genus in ref   -> untouched
## c03 NCBI says other, VMR absent   -> single side, untouched
## c04 VMR says other, NCBI absent   -> single side, untouched
## c05 VMR multi-family, one matches -> untouched
## c06 NCBI agrees, VMR rejects      -> conflict, untouched
## c07 genus not in ref              -> untouched
## c08 double reject, NCBI_n = 2     -> blank (ambiguity must not suppress action)
## c09 Family with surrounding space -> blank
## c10 Family NA                     -> untouched
## c11 Family empty string           -> untouched
tt <- data.table(
  contig_id = sprintf("c%02d", 1:11),
  Family = c("FamRow", '"FamRow"', "FamRow", "FamRow", "FamRow", "FamRow",
             "FamRow", "FamRow", " FamRow ", NA, ""),
  Genus = c("gblank", "gok", "gncbionly", "gvmronly", "gmulti", "gncbiok",
            "gabsent", "gambig", "gblank", "gblank", "gblank"),
  Species = sprintf("sp%02d", 1:11)
)
EXP_BLANK <- c("c01", "c08", "c09")

res <- env$enforce_rank_containment(copy(tt), TMP)
got <- res[is.na(Genus), contig_id]
stopifnot(setequal(got, EXP_BLANK))
stopifnot(all(is.na(res[contig_id %in% EXP_BLANK, Species])))
stopifnot(all(!is.na(res[!(contig_id %in% EXP_BLANK), Genus])))
stopifnot(identical(res$Family, tt$Family))
cat(sprintf("[3] synthetic cases PASS: blank = %s (of 11)\n", paste(got, collapse = ",")))

## ------------------------------------------------------------- fail-safe cases
res2 <- env$enforce_rank_containment(copy(tt), "/tmp/__no_such_ref_file__.tsv")
stopifnot(identical(res2$Genus, tt$Genus), identical(res2$Species, tt$Species))
BAD <- tempfile(fileext = ".tsv")
writeLines(c("Genus\tFoo", "gblank\tx"), BAD)
res3 <- env$enforce_rank_containment(copy(tt), BAD)
stopifnot(identical(res3$Genus, tt$Genus))
res4 <- env$enforce_rank_containment(copy(tt), TMP, blank_species = FALSE)
stopifnot(all(is.na(res4[contig_id %in% EXP_BLANK, Genus])))
stopifnot(all(!is.na(res4[contig_id %in% EXP_BLANK, Species])))
res5 <- env$enforce_rank_containment(tt[0], TMP)
stopifnot(nrow(res5) == 0)
res6 <- env$enforce_rank_containment(copy(tt), NULL)
stopifnot(identical(res6$Genus, tt$Genus))
cat("[4] fail-safe PASS: missing ref / wrong cols / blank_species=FALSE / empty table / ref=NULL\n")

## ------------------------------------------------------------------ real data
tbl <- fread(TBL, sep = "\t", quote = "", colClasses = "character", showProgress = FALSE)
setnames(tbl, strip(names(tbl)))
stopifnot(all(c("contig_id", "Family", "Genus", "Species") %in% names(tbl)))
cat(sprintf("[5] real table: %d rows x %d cols\n", nrow(tbl), ncol(tbl)))

res <- env$enforce_rank_containment(copy(tbl), REF)
b0 <- strip(tbl$Genus)
a0 <- strip(res$Genus)
bl <- strip(res$contig_id)[!is.na(b0) & b0 != "" & (is.na(a0) | a0 == "")]

cal <- fread(CAL, sep = "\t", quote = "", colClasses = "character", showProgress = FALSE)
setnames(cal, strip(names(cal)))
exp <- strip(cal$contig_id)[strip(cal$calib_action) == "blank"]
cat(sprintf("[5] R blanks = %d ; calib_action==blank = %d\n", length(bl), length(exp)))
cat(sprintf("[5] blank sets identical: %s\n", setequal(bl, exp)))
if (!setequal(bl, exp)) {
  cat(sprintf("      only in R: %s\n", paste(head(setdiff(bl, exp), 5), collapse = ",")))
  cat(sprintf("      only in py: %s\n", paste(head(setdiff(exp, bl), 5), collapse = ",")))
}
stopifnot(all(is.na(strip(res$Species)[strip(res$contig_id) %in% bl])))
stopifnot(identical(strip(res$Family), strip(tbl$Family)))
stopifnot(identical(strip(res$Order), strip(tbl$Order)))
stopifnot(identical(strip(res$Realm), strip(tbl$Realm)))
stopifnot(identical(strip(res$Genus)[!(strip(res$contig_id) %in% bl)], b0[!(strip(res$contig_id) %in% bl)]))
cat("[6] Species blanked with Genus; Family/Order/Realm/other-Genus untouched\n")
stopifnot(setequal(bl, exp))
cat("ALL R-SIDE TESTS PASS\n")
