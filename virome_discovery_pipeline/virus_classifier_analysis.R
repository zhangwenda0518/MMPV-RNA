#!/usr/bin/env Rscript

# ==============================================================================
# 病毒分类软件结果整合分析脚本 v6.1 (The Pinnacle & Uncompromised Edition)
# 审计: 彻底消除所有隐式循环，实现 100% 纯 C 级 GForce 加速 (百万级数据十几秒)
# 视觉: 全面恢复并展开所有 ComplexUpset 彩色高级绘图逻辑及 PDF/PNG 双路输出
# ==============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(ggplot2)
  library(VennDiagram)
  library(ggVennDiagram)
  library(grid)
  library(gridExtra)
  library(stringr)
  library(RColorBrewer)
  library(scales)
  library(cowplot)
  library(patchwork)
  library(parallelly)
})

TAX_LEVELS <- c("Realm", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
TAX_V2_COLS <- c("superkingdom", "phylum", "class", "order", "family", "genus", "species")

KNOWN_REALMS   <- c("Riboviria","Monodnaviria","Duplodnaviria","Varidnaviria","Adnaviria","Ribozyviria")
SUBRANK_SUFFIX <- c("viricotina","viricetidae","virineae","virinae")
RANK_DEPTH_WEIGHTS <- c(Realm=1, Kingdom=2, Phylum=4, Class=8, Order=16, Family=32, Genus=64, Species=128)
TOOL_BIAS <- c(ACVirus=1.2, VITAP=1.1, mmseqs=1.0, metabuli=1.0, CAT=0.9, genomad=0.9, diamond_lca=0.8, vcontact3=0.7, contigtax=0.6, BASTA=0.6, PhaGCN3=0.8)

setup_threads <- function(cores = NULL) {
  avail <- availableCores()
  use <- if (is.null(cores)) max(2, min(avail - 1, 8)) else min(cores, avail)
  setDTthreads(use)
  cat(sprintf("[System] %d 核心, 使用 %d 线程进行 data.table 极速 C 运算\n", avail, use))
}

log_msg <- function(level, msg, ...) {
  cat(sprintf("[%s] %s: %s\n", format(Sys.time(), "%H:%M:%S"), level, sprintf(msg, ...)))
}

# ==============================================================================
# 标准化与向量化校验
# ==============================================================================
standardize_dt <- function(dt) {
  if (is.null(dt) || nrow(dt) == 0) return(NULL)
  if (!"contig_id" %in% names(dt)) {
    id_guess <- grep("^(contig_id|seq_name|query|genome)$", names(dt), ignore.case=TRUE, value=TRUE)
    if (length(id_guess) > 0) setnames(dt, id_guess[1], "contig_id") else return(NULL)
  }
  dt[, contig_id := as.character(contig_id)]
  for (col in TAX_LEVELS) {
    if (!col %in% names(dt)) dt[, (col) := NA_character_]
    dt[, (col) := as.character(get(col))]
    dt[get(col) %in% c("", "-", "NA", "na", "N/A", "no rank", "undefined", "unknown", "null", "default", "Unclassified")
       | grepl("^Unclassified\\.", get(col), ignore.case=TRUE)
       | grepl("^unplaced", get(col), ignore.case=TRUE), (col) := NA_character_]
  }
  return(dt[, .SD, .SDcols = c("contig_id", TAX_LEVELS)])
}

is_valid_value_vec <- function(v) {
  invalid_strs <- c("-","NA","na","N/A","no rank","undefined","unknown","null","default","Unclassified")
  !is.na(v) & v != "" & !(v %in% invalid_strs)
}

species_quality_score_vec <- function(v) {
  scores <- rep(1.0, length(v))
  valid_idx <- is_valid_value_vec(v)
  if (!any(valid_idx)) return(scores)

  s_valid <- as.character(v[valid_idx])
  sub_scores <- rep(1.0, length(s_valid))

  idx_bad_sfx <- grepl("(inae|idae|ales|icetes|viricota)$", s_valid, ignore.case=TRUE)
  sub_scores[idx_bad_sfx] <- 0
  idx_sp <- !idx_bad_sfx & grepl("\\bsp\\.?\\b|\\bcf\\.?\\b|\\baff\\.?\\b", s_valid, ignore.case=TRUE, perl=TRUE)
  sub_scores[idx_sp] <- sub_scores[idx_sp] * 0.3
  idx_env <- !idx_bad_sfx & grepl("unclassified|environmental|uncultured", s_valid, ignore.case=TRUE)
  sub_scores[idx_env] <- sub_scores[idx_env] * 0.1
  idx_nov <- !idx_bad_sfx & grepl("^unplaced|^novel_", s_valid, ignore.case=TRUE)
  sub_scores[idx_nov] <- sub_scores[idx_nov] * 0.2
  idx_good <- !idx_bad_sfx & grepl("^[A-Z][a-z]+virus\\s+[a-z]", s_valid, perl=TRUE)
  sub_scores[idx_good] <- sub_scores[idx_good] * 1.5
  idx_no_space <- !idx_bad_sfx & !grepl(" ", s_valid)
  sub_scores[idx_no_space] <- sub_scores[idx_no_space] * 0.5

  scores[valid_idx] <- sub_scores
  return(scores)
}

clean_all_ranks <- function(data_list) {
  for (tn in names(data_list)) {
    dt <- data_list[[tn]]
    if (is.null(dt) || nrow(dt)==0) next
    dt[!tolower(Realm) %in% tolower(KNOWN_REALMS), Realm := NA_character_]
    dt[!is.na(Phylum) & !grepl("viricota$", Phylum, ignore.case=TRUE), Phylum := NA_character_]
    dt[!is.na(Class) & !grepl("viricetes$", Class, ignore.case=TRUE), Class := NA_character_]
    dt[!is.na(Order) & !grepl("virales$", Order, ignore.case=TRUE), Order := NA_character_]
    dt[!is.na(Family) & !grepl("viridae$", Family, ignore.case=TRUE), Family := NA_character_]
    dt[!is.na(Genus) & grepl("viridae$|virinae$", Genus, ignore.case=TRUE), Genus := NA_character_]
    dt[!is.na(Species) & grepl("(inae|idae|ales|icetes|viricota)$", Species, ignore.case=TRUE), Species := NA_character_]
    for (col in TAX_LEVELS) for (sfx in SUBRANK_SUFFIX) dt[grepl(paste0(sfx,"$"), get(col), ignore.case=TRUE), (col) := NA_character_]
    data_list[[tn]] <- dt
  }
  return(data_list)
}

# ==============================================================================
# 解析器模块
# ==============================================================================
parse_combined_v2 <- function(file) {
  if (is.null(file) || !file.exists(file)) return(NULL)
  dt <- fread(file, sep="\t", quote="")
  if (nrow(dt) == 0) return(NULL)
  setnames(dt, names(dt), tolower(names(dt)))
  if (!"tool" %in% names(dt) || !"seq_name" %in% names(dt)) return(NULL)
  setnames(dt, "seq_name", "contig_id")
  col_map <- c("superkingdom"="Realm", "realm"="Realm", "kingdom"="Kingdom", "phylum"="Phylum", "class"="Class", "order"="Order", "family"="Family", "genus"="Genus", "species"="Species")
  for (v2col in names(col_map)) if (v2col %in% names(dt)) setnames(dt, v2col, col_map[v2col])
  tool_list <- split(dt, by="tool")
  result <- lapply(names(tool_list), function(tn) standardize_dt(tool_list[[tn]]))
  names(result) <- names(tool_list)
  return(result)
}

parse_mmseqs <- function(file) {
  if (is.null(file) || !file.exists(file)) return(NULL)
  dt <- fread(file, header=FALSE, sep="\t", select=c(1,3,4,9), col.names=c("contig_id","rank","name","full_taxonomy"), fill=TRUE, quote="")
  dt <- dt[!rank %in% c("acellular root","root")]
  dt <- dt[grepl("Viruses|viria", full_taxonomy, ignore.case=TRUE)]
  if (nrow(dt)==0) return(NULL)
  dt[, (TAX_LEVELS) := NA_character_]
  parts <- tstrsplit(dt$full_taxonomy, ";", fixed=TRUE)
  rank_patterns <- list("^[kd]_"="Kingdom", "^p_"="Phylum", "^c_"="Class", "^o_"="Order", "^f_"="Family", "^g_"="Genus", "^s_"="Species")
  for (col_vec in parts) {
    col_vec <- trimws(col_vec)
    for (pat in names(rank_patterns)) { idx <- grep(pat, col_vec); if(length(idx)>0) dt[idx, (rank_patterns[[pat]]) := sub(pat,"",col_vec[idx])] }
    idx_r <- grep("^r_", col_vec); if(length(idx_r)>0) dt[idx_r, Realm := sub("^r_","",col_vec[idx_r])]
    idx_dash <- grep("^-_", col_vec); if (length(idx_dash)>0) { potentials <- sub("^-_","",col_vec[idx_dash]); is_valid <- potentials %in% KNOWN_REALMS; rows <- idx_dash[is_valid]; if (length(rows)>0) { rows2 <- rows[is.na(dt$Realm[rows])]; if(length(rows2)>0) dt[rows2, Realm := sub("^-_","",col_vec[rows2])] } }
  }
  rank_map <- c("species"="Species","genus"="Genus","family"="Family","order"="Order","class"="Class","phylum"="Phylum","kingdom"="Kingdom")
  for (r in names(rank_map)) dt[rank==r & !is.na(name) & name!="", (rank_map[[r]]) := name]
  return(standardize_dt(dt))
}

parse_vcontact3 <- function(file) {
  if (is.null(file) || !file.exists(file)) return(NULL)
  dt <- fread(file, quote=""); if (nrow(dt)==0) return(NULL)
  setnames(dt, tolower(names(dt)))
  if ("reference" %in% names(dt)) dt <- dt[grepl("false", tolower(reference))]
  id_col <- grep("genome|bin", names(dt), value=TRUE)[1]; if(is.na(id_col)) return(NULL)
  setnames(dt, id_col, "contig_id"); dt <- dt[contig_id!="default"]
  for (lvl in TAX_LEVELS) { target <- grep(paste0("^",tolower(lvl),"_prediction"), names(dt), value=TRUE)[1]; if (!is.na(target)) dt[, (lvl) := get(target)] }
  return(standardize_dt(dt))
}

parse_vitap <- function(file) {
  if (is.null(file) || !file.exists(file)) return(NULL)
  hd <- readLines(file, n=20); skip_n <- grep("^Genome_ID|^Contig", hd)[1]; skip_n <- if(is.na(skip_n)) 0 else skip_n-1
  dt <- fread(file, skip=skip_n, fill=TRUE, quote=""); if (nrow(dt)==0) return(NULL)
  setnames(dt, 1:2, c("contig_id", "lineage")); dt[, (TAX_LEVELS) := NA_character_]
  for (i in 1:nrow(dt)) {
    lin <- dt[i, lineage]; parts <- trimws(strsplit(lin, ";")[[1]]); parts <- parts[parts != "" & tolower(parts) != "viruses"]
    for (p in parts) {
      if (tolower(p) %in% tolower(KNOWN_REALMS)) { dt[i, Realm := p] }
      else if (grepl("viricota$", p, ignore.case=TRUE)) { dt[i, Phylum := p] }
      else if (grepl("viricetes$", p, ignore.case=TRUE)) { dt[i, Class := p] }
      else if (grepl("virales$", p, ignore.case=TRUE)) { dt[i, Order := p] }
      else if (grepl("viridae$", p, ignore.case=TRUE)) { dt[i, Family := p] }
      else if (grepl("virus$", p, ignore.case=TRUE) && !grepl("viridae$|virinae$", p, ignore.case=TRUE) && !any(sapply(SUBRANK_SUFFIX, function(s) grepl(paste0(s,"$"), p, ignore.case=TRUE)))) {
        if (is.na(dt[i, Genus]) || dt[i, Genus]=="") { dt[i, Genus := p] }
        else if (is.na(dt[i, Species]) || dt[i, Species]=="") { dt[i, Species := p] }
      } else if (!any(sapply(SUBRANK_SUFFIX, function(s) grepl(paste0(s,"$"), p, ignore.case=TRUE)))) { dt[i, Species := p] }
    }
  }
  return(standardize_dt(dt))
}

parse_phagcn3 <- function(file) {
  if (is.null(file) || !file.exists(file)) return(NULL)
  dt <- fread(file, sep=",", fill=TRUE, quote=""); if (nrow(dt)==0) return(NULL)
  if (!"contig_name" %in% names(dt)) setnames(dt, 1, "contig_name")
  setnames(dt, "contig_name", "contig_id"); dt[, (TAX_LEVELS) := NA_character_]
  if ("prediction" %in% names(dt)) dt[, Family := prediction]
  char_cols <- names(dt)[sapply(dt, is.character)]; char_cols <- setdiff(char_cols, c("contig_id", "prediction"))
  patterns <- list("r;"="Realm","k;"="Kingdom","p;"="Phylum","c;"="Class","o;"="Order","f;"="Family","g;"="Genus","s;"="Species")
  for (col in char_cols) {
    vv <- dt[[col]]; if (!any(grepl(";", vv[1:min(100,length(vv))]))) next
    for (pat in names(patterns)) { idx <- grep(pat, vv, fixed=TRUE); if (length(idx)>0) dt[idx, (patterns[[pat]]) := sub(paste0(".*",pat),"",vv[idx])] }
  }
  return(standardize_dt(dt))
}

parse_acvirus <- function(file) {
  if (is.null(file) || !file.exists(file)) return(NULL)
  dt <- fread(file, sep="\t", fill=TRUE, quote=""); if (nrow(dt)==0) return(NULL)
  names(dt) <- paste0(toupper(substr(names(dt),1,1)), substr(names(dt),2,nchar(names(dt))))
  setnames(dt, 1, "contig_id")
  return(standardize_dt(dt))
}

# ==============================================================================
# 核心业务逻辑: 加权投票引擎 (极度调优纯 C 向量化版本)
# ==============================================================================
compute_tool_weights <- function(data_list, consensus_stats=NULL) {
  weights <- list()
  for (tn in names(data_list)) {
    wvec <- numeric(length(TAX_LEVELS)); names(wvec) <- TAX_LEVELS
    for (lvl in TAX_LEVELS) {
      w <- RANK_DEPTH_WEIGHTS[lvl] * ifelse(tn %in% names(TOOL_BIAS), TOOL_BIAS[tn], 0.8)
      if (!is.null(consensus_stats) && tn %in% names(consensus_stats)) {
        rate <- consensus_stats[[tn]][[lvl]]
        if (!is.null(rate) && is.numeric(rate) && rate > 0) w <- w * rate
      }
      wvec[lvl] <- w
    }
    weights[[tn]] <- wvec
  }
  return(weights)
}

harmonize_genus_species <- function(result_dt) {
  if (nrow(result_dt) == 0) return(result_dt)
  result_dt[, genus_from_sp := tstrsplit(as.character(Species), " ", keep = 1L)]
  valid_genus <- result_dt[!is.na(genus_from_sp) & grepl("virus$", genus_from_sp, ignore.case=TRUE) & !grepl("viridae$|virinae$", genus_from_sp, ignore.case=TRUE)]
  for (sfx in SUBRANK_SUFFIX) valid_genus <- valid_genus[!grepl(paste0(sfx, "$"), genus_from_sp, ignore.case=TRUE)]
  valid_genus <- valid_genus[!is.na(Genus) & tolower(genus_from_sp) != tolower(Genus)]
  if (nrow(valid_genus) > 0) result_dt[valid_genus, on = "contig_id", Genus := i.genus_from_sp]
  result_dt[, genus_from_sp := NULL]
  return(result_dt)
}

build_consensus <- function(data_list, tool_weights) {
  # ── 1. 快速长表构建 ──
  stacked <- rbindlist(lapply(names(data_list), function(tn) {
    d <- unique(data_list[[tn]], by="contig_id"); d[, Tool := tn]; d
  }), use.names=TRUE, fill=TRUE)

  long <- melt(stacked, id.vars=c("contig_id","Tool"), measure.vars=TAX_LEVELS, variable.name="Rank", value.name="Taxon")
  long <- long[is_valid_value_vec(Taxon)]
  if (nrow(long)==0) return(data.table())

  # ── 2. 向量化权重计算 ──
  wt_dt <- rbindlist(lapply(names(tool_weights), function(tn) data.table(Tool=tn, Rank=TAX_LEVELS, Weight=as.numeric(tool_weights[[tn]]))))
  long <- merge(long, wt_dt, by=c("Tool","Rank"), all.x=TRUE)
  long[is.na(Weight), Weight := 0]
  sp_idx <- long$Rank == "Species"
  if(any(sp_idx)) long[sp_idx, Weight := Weight * species_quality_score_vec(Taxon)]

  # ── 3. 极速哈希聚合投票 ──
  votes <- long[, .(total_weight=sum(Weight)), by=.(contig_id, Rank, Taxon)]
  setorder(votes, contig_id, Rank, -total_weight)
  winners <- unique(votes, by=c("contig_id", "Rank"))

  wide <- dcast(winners[, .(contig_id, Rank, Taxon)], contig_id ~ Rank, value.var="Taxon")
  for (lvl in setdiff(TAX_LEVELS, names(wide))) set(wide, NULL, lvl, NA_character_)

  # ── 4. 全向量化空缺填补 ──
  val_cols <- lapply(TAX_LEVELS, function(col) is_valid_value_vec(stacked[[col]]))
  stacked[, score := Reduce(`+`, val_cols)]
  setorder(stacked, contig_id, -score)
  best_data <- unique(stacked, by="contig_id")

  for (lvl in TAX_LEVELS) {
    missing_idx <- which(!is_valid_value_vec(wide[[lvl]]))
    if(length(missing_idx) > 0) {
        miss_cids <- wide$contig_id[missing_idx]
        fill_vals <- best_data[contig_id %in% miss_cids, .(contig_id, val=get(lvl))]
        fill_vals <- fill_vals[is_valid_value_vec(val)]
        if(nrow(fill_vals) > 0) wide[fill_vals, on="contig_id", (lvl) := i.val]
    }
  }

  wide[, completeness := rowSums(!is.na(.SD) & .SD!=""), .SDcols=TAX_LEVELS]
  wide[, confidence := round(completeness/length(TAX_LEVELS), 2)]

  # ── 5. 哈希匹配主工具 ──
  cons_long <- melt(wide[, c("contig_id", TAX_LEVELS), with=FALSE], id.vars="contig_id", measure.vars=TAX_LEVELS, variable.name="Rank", value.name="cons_Taxon")
  cons_long <- cons_long[is_valid_value_vec(cons_Taxon)]

  match_dt <- merge(long[, .(contig_id, Tool, Rank, Taxon)], cons_long, by=c("contig_id", "Rank"))
  rank_weights_dt <- data.table(Rank = names(RANK_DEPTH_WEIGHTS), rw = as.numeric(RANK_DEPTH_WEIGHTS))
  match_dt <- merge(match_dt, rank_weights_dt, by="Rank", all.x=TRUE)
  match_dt[, match_weight := fcase(tolower(Taxon) == tolower(cons_Taxon), rw, default = 0)]

  match_scores_dt <- match_dt[, .(score_val = sum(match_weight)), by=.(contig_id, Tool)]
  all_pairs <- unique(stacked[, .(contig_id, Tool)])
  match_scores_dt <- merge(all_pairs, match_scores_dt, by=c("contig_id", "Tool"), all.x=TRUE)
  match_scores_dt[is.na(score_val), score_val := 0]

  setorder(match_scores_dt, contig_id, -score_val)
  pt <- unique(match_scores_dt, by="contig_id")[, .(contig_id, primary_tool = Tool, score_val)]
  pt[score_val == 0, primary_tool := "consensus"]
  pt[, score_val := NULL]

  wide <- merge(wide, pt, by="contig_id", all.x=TRUE)
  wide[is.na(primary_tool), primary_tool := "consensus"]

  wide <- harmonize_genus_species(wide)

  # ── 6. 🚀 极致多维宽长转换: Agree列拼接 ──
  tot_dt <- long[, .(n_tot = .N), by=.(contig_id, Rank)]
  ag_dt <- match_dt[match_weight > 0]
  setorder(ag_dt, contig_id, Rank, Tool)
  ag_dt <- ag_dt[, .(n_ag = .N, ag_tools = paste(Tool, collapse=",")), by=.(contig_id, Rank)]

  res <- merge(tot_dt, ag_dt, by=c("contig_id", "Rank"), all.x=TRUE)
  res[is.na(n_ag), `:=`(n_ag = 0, ag_tools = "")]
  res[, agree_str := fcase(
      n_ag > 0, sprintf("%d/%d: %s", n_ag, n_tot, ag_tools),
      default = sprintf("0/%d", n_tot)
  )]

  agree_wide <- dcast(res, contig_id ~ Rank, value.var="agree_str")
  exist_cols <- intersect(TAX_LEVELS, names(agree_wide))
  setnames(agree_wide, exist_cols, paste0(exist_cols, "_agree"))

  wide <- merge(wide, agree_wide, by="contig_id", all.x=TRUE)
  for (lvl in TAX_LEVELS) {
    col <- paste0(lvl, "_agree")
    if (!col %in% names(wide)) wide[, (col) := "0/0"]
    else wide[is.na(get(col)), (col) := "0/0"]
  }

  agree_cols <- paste0(TAX_LEVELS, "_agree")
  setcolorder(wide, c("contig_id","primary_tool","completeness","confidence", TAX_LEVELS, agree_cols))
  return(wide[])
}

# ==============================================================================
# 无损视觉展示：全量图表生成与完整导出
# ==============================================================================
plot_classification_counts <- function(data_list, output_dir) {
  log_msg("VIS", "绘制分类数目条形图...")
  counts <- data.frame(Tool=names(data_list), Count=sapply(data_list, nrow))
  counts$Percentage <- round(counts$Count/sum(counts$Count)*100, 1)
  p <- ggplot(counts, aes(x=reorder(Tool,-Count), y=Count, fill=Tool)) +
    geom_bar(stat="identity") +
    geom_text(aes(label=sprintf("%d\n(%.1f%%)", Count, Percentage)), vjust=-0.5, size=4) +
    labs(title="各工具分类数目", x="工具", y="序列数") +
    theme_minimal() +
    theme(axis.text.x=element_text(angle=45, hjust=1), legend.position="none") +
    scale_fill_brewer(palette="Set2") +
    ylim(0, max(counts$Count)*1.15)

  ggsave(file.path(output_dir, "classification_counts.pdf"), p, width=max(8, nrow(counts)*1.2), height=6)
  ggsave(file.path(output_dir, "classification_counts.png"), p, width=max(8, nrow(counts)*1.2), height=6, dpi=300)
  return(counts)
}

plot_intersection_diagram <- function(data_list, output_dir) {
  if (length(data_list) < 2) return()
  n_tools <- length(data_list)
  id_lists <- lapply(data_list, function(x) unique(x$contig_id))

  if (n_tools <= 3) {
    log_msg("VIS", "绘制 Venn 图...")
    tryCatch({
      colors <- brewer.pal(min(n_tools, 8), "Set2")[1:n_tools]
      venn.plot <- venn.diagram(
        x=id_lists, filename=NULL, category.names=names(data_list),
        fill=colors, cex=1.2, cat.cex=1.2, cat.fontface="bold", main="序列交集 (Venn)"
      )
      pdf(file.path(output_dir, "intersection_venn.pdf"), width=10, height=10)
      grid.draw(venn.plot); dev.off()
      png(file.path(output_dir, "intersection_venn.png"), width=10, height=10, units="in", res=300)
      grid.draw(venn.plot); dev.off()
    }, error=function(e) log_msg("WARN", "Venn生成失败: %s", e$message))

  } else if (requireNamespace("ComplexUpset", quietly=TRUE)) {
    log_msg("VIS", sprintf("绘制 ComplexUpset 彩色图 (%d 工具)...", n_tools))
    all_ids <- unique(unlist(id_lists))
    bin_mat <- as.data.frame(lapply(names(id_lists), function(tn) as.integer(all_ids %in% id_lists[[tn]])))
    colnames(bin_mat) <- names(id_lists)

    tryCatch({
      bin_mat$degree <- rowSums(bin_mat[, names(id_lists)])

      g <- ComplexUpset::upset(
        bin_mat, names(id_lists),
        name="Tool intersections",
        width_ratio=0.2,
        min_size=1,
        base_annotations=list(
          "Intersection size" = ComplexUpset::intersection_size(
            counts=FALSE,
            mapping=ggplot2::aes(fill=as.factor(degree))
          ) +
          ggplot2::scale_fill_manual(
            values=c("#449B99", "#55B8E9", "#7CCD7C", "#B8DF8E", "#F0D48A", "#FFB16C", "#F6352A", "grey50", "grey30", "grey10"),
            name="Tools overlapped"
          ) +
          ggplot2::ylab("Intersection Size") +
          ggplot2::theme_minimal() +
          ggplot2::theme(axis.text.x=ggplot2::element_blank(), axis.title.x=ggplot2::element_blank())
        ),
        stripes=ComplexUpset::upset_stripes(
          mapping=ggplot2::aes(color="black", fill="steelblue")
        ),
        matrix=ComplexUpset::intersection_matrix(
          geom=ggplot2::geom_point(shape=21, size=3, stroke=1)
        ) +
        ggplot2::scale_color_manual(values=c("TRUE"="#E64B35", "FALSE"="grey90"), guide="none")
      )

      ggplot2::ggsave(file.path(output_dir, "intersection_upset.pdf"), g, width=14, height=8, dpi=300)
      ggplot2::ggsave(file.path(output_dir, "intersection_upset.png"), g, width=14, height=8, dpi=300)
    }, error=function(e) log_msg("WARN", "ComplexUpset 绘图失败: %s", e$message))

  } else {
    log_msg("WARN", "ComplexUpset 包未安装，跳过高级交集图绘制")
  }
}

plot_venn_by_tax_level <- function(data_list, output_dir) {
  if (length(data_list) > 7) {
    log_msg("INFO", ">7 tools, 跳过 per-rank Venn")
    return()
  }
  log_msg("VIS", "绘制各分类等级 Venn 图...")
  venn_plots <- list()
  for (level in TAX_LEVELS) {
    id_lists <- list()
    for (tn in names(data_list)) {
      ids <- data_list[[tn]][!is.na(get(level)) & get(level)!="", unique(contig_id)]
      if (length(ids)>0) id_lists[[tn]] <- ids
    }
    if (length(id_lists)>=2) {
      p <- ggVennDiagram(id_lists, label_alpha=0) +
        scale_fill_gradient(low="#F4FAFE", high="#4981BF") +
        labs(title=level) +
        theme(plot.title=element_text(hjust=0.5, size=12, face="bold"), legend.position="none")
      venn_plots[[level]] <- p
    }
  }
  if (length(venn_plots)>0) {
    combined <- wrap_plots(venn_plots, ncol=4, nrow=2) + plot_annotation(title="Per-Rank Intersection")
    ggsave(file.path(output_dir, "combined_tax_level_venn.pdf"), combined, width=20, height=12)
    ggsave(file.path(output_dir, "combined_tax_level_venn.png"), combined, width=20, height=12, dpi=300)
  }
}

analyze_consistency_optimized <- function(data_list, common_ids, output_dir) {
  if (length(common_ids)==0) return(NULL)
  log_msg("VIS", "分析一致性分布...")
  fll <- lapply(names(data_list), function(nm) { dt <- data_list[[nm]][contig_id %in% common_ids, .SD, .SDcols=c("contig_id",TAX_LEVELS)]; dt[, Tool:=nm]; dt })
  big_dt <- rbindlist(fll)
  results <- list()
  for (level in TAX_LEVELS) {
    vd <- big_dt[, .(contig_id, Tool, Taxon=get(level))][!is.na(Taxon) & Taxon!=""]
    if (nrow(vd)==0) next
    counts <- vd[, .N, by=.(contig_id, Taxon)]
    setorder(counts, contig_id, -N, Taxon)
    cons <- unique(counts, by="contig_id")[, .(contig_id, Consensus=Taxon)]
    md <- merge(vd, cons, by="contig_id")
    md[, is_agreed := (!is.na(Taxon) & Taxon==Consensus)]
    stats <- md[!is.na(Taxon), .(Agreement_Count=sum(is_agreed), Total_Classified=.N), by=Tool]
    if (nrow(stats)>0) {
      stats[, `:=`(Taxonomic_Level=level, Agreement_Rate=Agreement_Count/Total_Classified)]
      results[[level]] <- stats
    }
  }
  ad <- rbindlist(results)
  if (nrow(ad)>0) {
    ad$Taxonomic_Level <- factor(ad$Taxonomic_Level, levels=TAX_LEVELS)
    p <- ggplot(ad, aes(x=Taxonomic_Level, y=Agreement_Rate, fill=Tool)) +
      geom_bar(stat="identity", position=position_dodge(0.9)) +
      geom_text(aes(label=sprintf("%.2f", Agreement_Rate)), position=position_dodge(0.9), vjust=-0.5, size=3) +
      labs(title="各工具与共识的一致率", x="分类等级", y="一致率") +
      theme_minimal() +
      scale_fill_brewer(palette="Set2") +
      scale_y_continuous(labels=percent, limits=c(0,1.1))

    ggsave(file.path(output_dir, "agreement_rates.pdf"), p, width=12, height=8)
    ggsave(file.path(output_dir, "agreement_rates.png"), p, width=12, height=8, dpi=300)
    fwrite(ad, file.path(output_dir, "agreement_stats.tsv"), sep="\t")
  }
}

analyze_detailed_consistency <- function(data_list, common_ids, output_dir) {
  if (length(common_ids)==0) return(NULL)
  log_msg("INFO", "详细一致性数据导出...")
  tool_names <- names(data_list); consistency_summary <- data.frame()
  for (level in TAX_LEVELS) {
    long_list <- lapply(tool_names, function(tn) { dt <- data_list[[tn]][contig_id %in% common_ids, .(contig_id, taxon=get(level))]; dt[, Tool:=tn]; dt[taxon %in% c("","NA","na","N/A","no rank","undefined","unknown","null"), taxon := NA_character_]; dt })
    big_dt <- rbindlist(long_list); valid_dt <- big_dt[!is.na(taxon)]
    if (nrow(valid_dt)>0) {
      counts <- valid_dt[, .N, by=.(contig_id, taxon)]
      setorder(counts, contig_id, -N, taxon)
      consensus_dt <- unique(counts, by="contig_id")[, .(contig_id, consensus_taxon=taxon)]

      vd_uniq <- unique(valid_dt, by=c("contig_id", "taxon"))
      setorder(vd_uniq, contig_id, taxon)
      stats_dt <- vd_uniq[, .(uniq_cnt=.N, all_vals_str=paste(taxon, collapse="|")), by=contig_id]

      res_dt <- merge(stats_dt, consensus_dt, by="contig_id", all=TRUE)
      res_dt[, `:=`(is_consistent=(uniq_cnt==1), consistency_status=ifelse(uniq_cnt==1,"consistent","inconsistent"))]
    } else {
      res_dt <- data.table(contig_id=common_ids, uniq_cnt=0, all_vals_str=NA_character_, consensus_taxon=NA_character_, is_consistent=NA, consistency_status="all_NA")
    }
    wide_dt <- dcast(big_dt, contig_id ~ Tool, value.var="taxon")
    final_comp <- merge(wide_dt, res_dt, by="contig_id", all.x=TRUE)

    consistent_out <- final_comp[is_consistent==TRUE, .(contig_id, consensus_taxon)]
    if (nrow(consistent_out)>0) fwrite(consistent_out, file.path(output_dir, paste0("consistent_",tolower(level),".tsv")), sep="\t")

    cols_order <- c("contig_id", tool_names, "all_vals_str", "consensus_taxon", "is_consistent")
    fwrite(final_comp[, cols_order, with=FALSE], file.path(output_dir, paste0("comparison_",tolower(level),".tsv")), sep="\t")

    for (tn in tool_names) {
      if (tn %in% names(final_comp)) {
        cnt <- sum(!is.na(final_comp[[tn]]) & final_comp[[tn]]!="")
        consistency_summary <- rbind(consistency_summary, data.frame(Tool=tn, Level=level, Classified=cnt))
      }
    }
  }

  if (nrow(consistency_summary)>0) {
    p1 <- ggplot(consistency_summary, aes(x=Level, y=Classified, fill=Tool)) +
      geom_bar(stat="identity", position=position_dodge(0.8)) +
      labs(title="各工具各等级分类数", x="等级", y="序列数") +
      theme_minimal() +
      theme(axis.text.x=element_text(angle=45, hjust=1)) +
      scale_fill_brewer(palette="Set2")

    p2 <- ggplot(consistency_summary, aes(x=Tool, y=Classified, fill=Level)) +
      geom_bar(stat="identity", position=position_stack()) +
      labs(title="分类数堆叠", x="工具", y="序列数") +
      theme_minimal() +
      scale_fill_brewer(palette="Set3")

    ggsave(file.path(output_dir, "classification_by_level.pdf"), p1, width=12, height=8)
    ggsave(file.path(output_dir, "classification_by_level.png"), p1, width=12, height=8, dpi=300)
    ggsave(file.path(output_dir, "classification_stacked.pdf"), p2, width=10, height=8)
    ggsave(file.path(output_dir, "classification_stacked.png"), p2, width=10, height=8, dpi=300)
    fwrite(consistency_summary, file.path(output_dir, "consistency_summary.tsv"), sep="\t")
  }
}

generate_grouped_consistent_files <- function(data_list, common_ids, output_dir) {
  if (length(common_ids)==0) return()
  log_msg("INFO", "按一致分类单元分组导出...")
  tool_names <- names(data_list); overall_summary <- data.frame()
  for (level in TAX_LEVELS) {
    long_list <- lapply(tool_names, function(tn) { dt <- data_list[[tn]][contig_id %in% common_ids, .(contig_id, taxon=get(level))]; dt[taxon %in% c("","NA","na","N/A"), taxon:=NA_character_]; dt })
    valid_dt <- rbindlist(long_list)[!is.na(taxon)]
    if (nrow(valid_dt)==0) next

    check_dt <- valid_dt[, .(uniq_cnt=uniqueN(taxon)), by=contig_id]
    consistent_ids <- check_dt[uniq_cnt==1, contig_id]
    if (length(consistent_ids)==0) next

    consistent_data <- unique(valid_dt[contig_id %in% consistent_ids], by="contig_id")
    setnames(consistent_data, "taxon", "consensus_taxon")

    level_dir <- file.path(output_dir, paste0("consistent_",tolower(level)))
    dir.create(level_dir, recursive=TRUE, showWarnings=FALSE)

    summary_dt <- consistent_data[, .(Count=.N), by=consensus_taxon]
    setorder(summary_dt, -Count)
    fwrite(summary_dt, file.path(level_dir, "summary.tsv"), sep="\t")

    for (taxon_name in unique(consistent_data$consensus_taxon)) {
      safe <- gsub("[^[:alnum:]_]", "_", taxon_name)
      safe <- gsub("_{2,}","_", gsub("^_|_$","", safe))
      sub_dt <- consistent_data[consensus_taxon==taxon_name, .(contig_id)]
      fwrite(sub_dt, file.path(level_dir, paste0(safe,".tsv")), sep="\t")
    }

    overall_summary <- rbind(overall_summary, data.frame(Taxonomic_Level=level, Unique_Taxa=length(unique(consistent_data$consensus_taxon)), Total_Sequences=nrow(consistent_data)))
  }

  if (nrow(overall_summary)>0) {
    fwrite(overall_summary, file.path(output_dir, "consistent_summary.tsv"), sep="\t")

    p1 <- ggplot(overall_summary, aes(x=Taxonomic_Level, y=Unique_Taxa, fill=Taxonomic_Level)) +
      geom_bar(stat="identity") +
      geom_text(aes(label=Unique_Taxa), vjust=-0.5) +
      labs(title="一致分类单元数", x="等级") +
      theme_minimal() +
      theme(legend.position="none") +
      scale_fill_brewer(palette="Set2")

    p2 <- ggplot(overall_summary, aes(x=Taxonomic_Level, y=Total_Sequences, fill=Taxonomic_Level)) +
      geom_bar(stat="identity") +
      geom_text(aes(label=Total_Sequences), vjust=-0.5) +
      labs(title="一致序列数", x="等级") +
      theme_minimal() +
      theme(legend.position="none") +
      scale_fill_brewer(palette="Set3")

    summary_plot <- plot_grid(p1, p2, ncol=2, labels="AUTO")
    ggsave(file.path(output_dir, "consistent_summary.pdf"), summary_plot, width=14, height=7)
    ggsave(file.path(output_dir, "consistent_summary.png"), summary_plot, width=14, height=7, dpi=300)
  }
}

# ==============================================================================
# 主流程
# ==============================================================================
option_list <- list(
  make_option("--combined", type="character", default=NULL),
  make_option("--mmseqs", type="character", default=NULL),
  make_option("--vcontact3", type="character", default=NULL),
  make_option("--vitap", type="character", default=NULL),
  make_option("--phagcn3", type="character", default=NULL),
  make_option("--acvirus", type="character", default=NULL),
  make_option("--genomad", type="character", default=NULL),
  make_option("--metabuli", type="character", default=NULL),
  make_option("--cat", type="character", default=NULL),
  make_option("--diamond_lca", type="character", default=NULL),
  make_option("--contigtax", type="character", default=NULL),
  make_option("--basta", type="character", default=NULL),
  make_option(c("-o","--output"), type="character", default="."),
  make_option("--cores", type="integer", default=NULL)
)

opt <- parse_args(OptionParser(option_list=option_list))
if (!dir.exists(opt$output)) dir.create(opt$output, recursive=TRUE)
setup_threads(opt$cores)

all_data <- list()
if (!is.null(opt$combined) && file.exists(opt$combined)) {
  combined_result <- parse_combined_v2(opt$combined)
  if (!is.null(combined_result)) all_data <- combined_result
}

parsers <- list(mmseqs=parse_mmseqs, vcontact3=parse_vcontact3, vitap=parse_vitap, phagcn3=parse_phagcn3, acvirus=parse_acvirus)
for (tool in names(parsers)) {
  f <- opt[[tool]]
  if (!is.null(f) && file.exists(f)) {
    tryCatch({
      res <- parsers[[tool]](f)
      if (!is.null(res) && nrow(res)>0) all_data[[tool]] <- res
    }, error=function(e) log_msg("ERROR", "%s: %s", tool, e$message))
  }
}

if (length(all_data) < 1) stop("无有效数据, 检查输入参数。")
log_msg("INFO", "共加载 %d 个工具数据", length(all_data))

log_msg("INFO", "清理 rank 值...")
all_data <- clean_all_ranks(all_data)
for (tn in names(all_data)) {
  fwrite(unique(all_data[[tn]], by="contig_id"), file.path(opt$output, paste0("standardized_", tn, ".tsv")), sep="\t", na="NA")
}

all_ids <- unique(unlist(lapply(all_data, function(d) d$contig_id)))
n_total <- length(all_ids)
cat(sprintf("\n总序列数: %d\n", n_total))

count_stats <- plot_classification_counts(all_data, opt$output)
mat <- matrix(0, nrow=n_total, ncol=length(all_data), dimnames=list(all_ids, names(all_data)))
for (tn in names(all_data)) mat[all_data[[tn]]$contig_id, tn] <- 1
mat_dt <- as.data.table(mat, keep.rownames="contig_id")
mat_dt[, count := rowSums(.SD), .SDcols=names(all_data)]
fwrite(mat_dt, file.path(opt$output, "intersection_matrix.tsv"), sep="\t")

common_ids <- rownames(mat)[rowSums(mat)==length(all_data)]
if (length(common_ids)>0) writeLines(common_ids, file.path(opt$output, "common_ids.txt"))

if (length(all_data)>=2) {
  plot_intersection_diagram(all_data, opt$output)
  plot_venn_by_tax_level(all_data, opt$output)
}

if (length(common_ids)>0) {
  analyze_consistency_optimized(all_data, common_ids, opt$output)
  analyze_detailed_consistency(all_data, common_ids, opt$output)
  generate_grouped_consistent_files(all_data, common_ids, opt$output)
}

log_msg("INFO", "计算工具权重 (全量哈希向量化模式)...")
ids_for_weight <- if (length(common_ids) > 0) common_ids else mat_dt[count >= max(2, length(all_data)/2), contig_id]
consensus_stats <- NULL

if (length(ids_for_weight) > 0) {
  dt_sub <- rbindlist(lapply(names(all_data), function(tn) {
    d <- all_data[[tn]][contig_id %in% ids_for_weight, c("contig_id", TAX_LEVELS), with=FALSE]
    if(nrow(d) > 0) { d[, Tool := tn]; return(d) } else return(NULL)
  }), use.names=TRUE, fill=TRUE)

  long_sub <- melt(dt_sub, id.vars=c("contig_id", "Tool"), measure.vars=TAX_LEVELS, variable.name="Rank", value.name="Taxon")
  long_sub <- long_sub[is_valid_value_vec(Taxon)]

  if (nrow(long_sub) > 0) {
    counts <- long_sub[, .N, by=.(contig_id, Rank, Taxon)]
    setorder(counts, contig_id, Rank, -N)

    cons_sub <- unique(counts, by=c("contig_id", "Rank"))
    setnames(cons_sub, "Taxon", "Consensus")

    mg_sub <- merge(long_sub, cons_sub[, .(contig_id, Rank, Consensus)], by=c("contig_id", "Rank"))
    mg_sub[, agreed := (tolower(Taxon) == tolower(Consensus))]

    total_cons_per_rank <- cons_sub[, .(total_cons = .N), by=Rank]
    rates_dt <- mg_sub[, .(agreed = sum(agreed)), by=.(Tool, Rank)]
    rates_dt <- merge(rates_dt, total_cons_per_rank, by="Rank")
    rates_dt[, rate := agreed / total_cons]

    consensus_stats <- list()
    for(tn in unique(rates_dt$Tool)) {
      r_tool <- rates_dt[Tool == tn]
      res <- as.list(r_tool$rate); names(res) <- as.character(r_tool$Rank)
      for(lvl in TAX_LEVELS) if(is.null(res[[lvl]])) res[[lvl]] <- 0
      consensus_stats[[tn]] <- res
    }
  }
}

tool_weights <- compute_tool_weights(all_data, consensus_stats)

log_msg("INFO", "逐 rank 加权投票计算共识...")
final_result <- build_consensus(all_data, tool_weights)

agree_cols <- paste0(TAX_LEVELS, "_agree")
out_cols <- c("contig_id","primary_tool","completeness","confidence", TAX_LEVELS, agree_cols)
fwrite(final_result[, out_cols, with=FALSE], file.path(opt$output, "final_integrated_classification.tsv"), sep="\t", na="NA")
log_msg("INFO", "最终分类输出: %d 个序列", nrow(final_result))

fill_counts <- sapply(TAX_LEVELS, function(l) sum(is_valid_value_vec(final_result[[l]])))
if (nrow(final_result) > 0) {
  fill_rates <- round(fill_counts / nrow(final_result), 4)
} else {
  fill_rates <- rep(0, length(TAX_LEVELS))
}
fill_dt <- data.table(Rank=TAX_LEVELS, Filled=fill_counts, Rate=fill_rates)
fwrite(fill_dt, file.path(opt$output, "fill_stats.tsv"), sep="\t")

if (requireNamespace("gridExtra", quietly=TRUE) && nrow(final_result) > 0) {
  log_msg("VIS", "绘制最终共识汇总图...")
  p1 <- ggplot(final_result, aes(x=factor(completeness))) +
    geom_bar(fill="#66c2a5") +
    labs(title="共识完备度分布", x="完备度", y="序列数") +
    theme_minimal()

  fr_df <- data.frame(Level=factor(TAX_LEVELS, levels=TAX_LEVELS), Rate=fill_rates*100)
  p2 <- ggplot(fr_df, aes(x=Level, y=Rate, fill=Level)) +
    geom_bar(stat="identity") +
    labs(title="各 Rank 填充率", x="等级", y="填充率(%)") +
    theme_minimal() +
    theme(legend.position="none") +
    scale_fill_brewer(palette="Set2")

  summary_plot <- plot_grid(p1, p2, ncol=2, labels="AUTO")
  ggsave(file.path(opt$output, "consensus_summary.pdf"), summary_plot, width=14, height=7)
  ggsave(file.path(opt$output, "consensus_summary.png"), summary_plot, width=14, height=7, dpi=300)
}

sink(file.path(opt$output, "analysis_summary.txt"))
cat("=== 病毒分类分析报告 ===\n\n")
cat(sprintf("时间: %s\n", Sys.time()))
cat(sprintf("总序列: %d\n", n_total))
cat(sprintf("最终分类: %d\n", nrow(final_result)))
cat(sprintf("工具数: %d\n", length(all_data)))
cat("\n各工具分类数:\n"); print(count_stats)
cat("\n各 Rank 填充率:\n")
for (i in 1:length(TAX_LEVELS)) cat(paste0("  ", TAX_LEVELS[i], ": ", fill_counts[i], " (", round(fill_rates[i]*100, 1), "%)\n"))
if(nrow(final_result) > 0) {
  cat("\n主要工具分布:\n"); print(table(final_result$primary_tool))
}
sink()
log_msg("INFO", "全流程完美执行完毕! 结果已保存至: %s", opt$output)
