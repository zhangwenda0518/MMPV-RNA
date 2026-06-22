#!/usr/bin/env Rscript

# ==============================================================================
# 病毒分类软件结果整合分析脚本 v4.0
# 兼容: virus_classifier.py 的 combined_taxonomy.tsv + 各工具独立输出
# 新增: genomad, metabuli, CAT, diamond_lca, contigtax, BASTA 解析器
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
# virus_classifier.py 的输出列名 (小写)
TAX_V2_COLS <- c("superkingdom", "phylum", "class", "order", "family", "genus", "species")

setup_threads <- function(cores = NULL) {
  avail <- availableCores()
  use <- if (is.null(cores)) max(2, min(avail - 1, 8)) else min(cores, avail)
  setDTthreads(use)
  cat(sprintf("[System] %d 核心, 使用 %d\n", avail, use))
}

log_msg <- function(level, msg, ...) {
  cat(sprintf("[%s] %s: %s\n", format(Sys.time(), "%H:%M:%S"), level, sprintf(msg, ...)))
}

# ==============================================================================
# 标准化
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
       | grepl("^Unclassified\\.", get(col))
       | grepl("^unplaced", get(col)), (col) := NA]
  }
  return(dt[, .SD, .SDcols = c("contig_id", TAX_LEVELS)])
}

# ==============================================================================
# 新增解析器: virus_classifier.py combined_taxonomy.tsv
# ==============================================================================
parse_combined_v2 <- function(file) {
  if (is.null(file) || !file.exists(file)) return(NULL)
  log_msg("INFO", "解析 combined_taxonomy: %s", basename(file))
  dt <- fread(file, sep="\t")
  if (nrow(dt) == 0) return(NULL)

  # 按 tool 列拆分
  setnames(dt, names(dt), tolower(names(dt)))
  if (!"tool" %in% names(dt) || !"seq_name" %in% names(dt)) {
    log_msg("WARN", "combined_taxonomy 缺少 tool/seq_name 列")
    return(NULL)
  }
  setnames(dt, "seq_name", "contig_id")

  # v2 列名映射到标准 TAX_LEVELS (已用 tolower, 用小写匹配)
  col_map <- c(
    "superkingdom"="Realm", "realm"="Realm", "kingdom"="Kingdom",
    "phylum"="Phylum", "class"="Class", "order"="Order",
    "family"="Family", "genus"="Genus", "species"="Species"
  )
  for (v2col in names(col_map)) {
    if (v2col %in% names(dt)) setnames(dt, v2col, col_map[v2col])
  }

  # 按 tool 拆分为列表
  tool_list <- split(dt, by="tool")
  result <- lapply(names(tool_list), function(tn) {
    standardize_dt(tool_list[[tn]])
  })
  names(result) <- names(tool_list)
  return(result)
}

# ==============================================================================
# 各工具独立解析器 (兼容旧版)
# ==============================================================================
parse_mmseqs <- function(file) {
  if (is.null(file) || !file.exists(file)) return(NULL)
  log_msg("INFO", "解析 MMseqs2: %s", basename(file))
  dt <- fread(file, header=FALSE, sep="\t", select=c(1,3,4,9),
              col.names=c("contig_id","rank","name","full_taxonomy"), fill=TRUE)
  dt <- dt[!rank %in% c("acellular root","root")]
  dt <- dt[grepl("Viruses|viria", full_taxonomy, ignore.case=TRUE)]
  if (nrow(dt)==0) return(NULL)

  dt[, (TAX_LEVELS) := NA_character_]
  parts <- tstrsplit(dt$full_taxonomy, ";", fixed=TRUE)
  rank_patterns <- list("^[kd]_"="Kingdom", "^p_"="Phylum", "^c_"="Class",
                        "^o_"="Order", "^f_"="Family", "^g_"="Genus", "^s_"="Species")
  valid_realms <- c("Adnaviria","Duplodnaviria","Monodnaviria","Riboviria","Ribozyviria","Varidnaviria")

  for (col_vec in parts) {
    col_vec <- trimws(col_vec)
    for (pat in names(rank_patterns)) {
      idx <- grep(pat, col_vec); if(length(idx)>0) dt[idx, (rank_patterns[[pat]]) := sub(pat,"",col_vec[idx])]
    }
    idx_r <- grep("^r_", col_vec); if(length(idx_r)>0) dt[idx_r, Realm := sub("^r_","",col_vec[idx_r])]
    idx_dash <- grep("^-_", col_vec)
    if (length(idx_dash)>0) {
      potentials <- sub("^-_","",col_vec[idx_dash])
      is_valid <- potentials %in% valid_realms
      rows <- idx_dash[is_valid]
      if (length(rows)>0) { rows2 <- rows[is.na(dt$Realm[rows])]; if(length(rows2)>0) dt[rows2, Realm := sub("^-_","",col_vec[rows2])] }
    }
  }
  rank_map <- c("species"="Species","genus"="Genus","family"="Family","order"="Order","class"="Class","phylum"="Phylum","kingdom"="Kingdom")
  for (r in names(rank_map)) dt[rank==r & !is.na(name) & name!="", (rank_map[[r]]) := name]
  return(standardize_dt(dt))
}

parse_vcontact3 <- function(file) {
  if (is.null(file) || !file.exists(file)) return(NULL)
  log_msg("INFO", "解析 vContact3: %s", basename(file))
  dt <- fread(file); if (nrow(dt)==0) return(NULL)
  setnames(dt, tolower(names(dt)))
  if ("reference" %in% names(dt)) dt <- dt[grepl("false", tolower(reference))]
  id_col <- grep("genome|bin", names(dt), value=TRUE)[1]; if(is.na(id_col)) return(NULL)
  setnames(dt, id_col, "contig_id"); dt <- dt[contig_id!="default"]
  for (lvl in TAX_LEVELS) {
    target <- grep(paste0("^",tolower(lvl),"_prediction"), names(dt), value=TRUE)[1]
    if (!is.na(target)) dt[, (lvl) := get(target)]
  }
  return(standardize_dt(dt))
}

parse_vitap <- function(file) {
  if (is.null(file) || !file.exists(file)) return(NULL)
  log_msg("INFO", "解析 VITAP: %s", basename(file))
  hd <- readLines(file, n=20); skip_n <- grep("^Genome_ID|^Contig", hd)[1]
  skip_n <- if(is.na(skip_n)) 0 else skip_n-1
  dt <- fread(file, skip=skip_n, fill=TRUE); if (nrow(dt)==0) return(NULL)
  setnames(dt, 1:2, c("contig_id", "lineage"))
  dt[, (TAX_LEVELS) := NA_character_]
  known_realms <- c("Adnaviria","Duplodnaviria","Monodnaviria","Riboviria","Ribozyviria","Varidnaviria")
  subrank_sfx <- c("viricotina","viricetidae","virineae","virinae")

  for (i in 1:nrow(dt)) {
    lin <- dt[i, lineage]
    parts <- trimws(strsplit(lin, ";")[[1]])
    parts <- parts[parts != "" & tolower(parts) != "viruses"]
    for (p in parts) {
      if (p %in% known_realms) { dt[i, Realm := p] }
      else if (grepl("viricota$", p)) { dt[i, Phylum := p] }
      else if (grepl("viricetes$", p)) { dt[i, Class := p] }
      else if (grepl("virales$", p)) { dt[i, Order := p] }
      else if (grepl("viridae$", p)) { dt[i, Family := p] }
      else if (grepl("virus$", p) && !grepl("viridae$|virinae$", p) &&
               !any(sapply(subrank_sfx, function(s) grepl(paste0(s,"$"), p)))) {
        if (is.na(dt[i, Genus]) || dt[i, Genus]=="") { dt[i, Genus := p] }
        else if (is.na(dt[i, Species]) || dt[i, Species]=="") { dt[i, Species := p] }
      } else if (!any(sapply(subrank_sfx, function(s) grepl(paste0(s,"$"), p)))) {
        dt[i, Species := p]
      }
    }
  }
  return(standardize_dt(dt))
}

parse_phagcn3 <- function(file) {
  if (is.null(file) || !file.exists(file)) return(NULL)
  log_msg("INFO", "解析 PhaGCN3: %s", basename(file))
  dt <- fread(file, sep=",", fill=TRUE); if (nrow(dt)==0) return(NULL)
  if (!"contig_name" %in% names(dt)) setnames(dt, 1, "contig_name")
  setnames(dt, "contig_name", "contig_id")
  dt[, (TAX_LEVELS) := NA_character_]
  if ("prediction" %in% names(dt)) dt[, Family := prediction]
  char_cols <- names(dt)[sapply(dt, is.character)]
  char_cols <- setdiff(char_cols, c("contig_id", "prediction"))
  patterns <- list("r;"="Realm","k;"="Kingdom","p;"="Phylum","c;"="Class","o;"="Order","f;"="Family","g;"="Genus","s;"="Species")
  for (col in char_cols) {
    vv <- dt[[col]]; if (!any(grepl(";", vv[1:min(100,length(vv))]))) next
    for (pat in names(patterns)) {
      idx <- grep(pat, vv, fixed=TRUE)
      if (length(idx)>0) dt[idx, (patterns[[pat]]) := sub(paste0(".*",pat),"",vv[idx])]
    }
  }
  return(standardize_dt(dt))
}

parse_acvirus <- function(file) {
  if (is.null(file) || !file.exists(file)) return(NULL)
  log_msg("INFO", "解析 ACVirus: %s", basename(file))
  dt <- fread(file, sep="\t", fill=TRUE); if (nrow(dt)==0) return(NULL)
  names(dt) <- paste0(toupper(substr(names(dt),1,1)), substr(names(dt),2,nchar(names(dt))))
  setnames(dt, 1, "contig_id")
  return(standardize_dt(dt))
}

# ==============================================================================
# 可视化 (保留原有全部函数)
# ==============================================================================
plot_classification_counts <- function(data_list, output_dir) {
  log_msg("VIS", "绘制分类数目条形图...")
  counts <- data.frame(Tool=names(data_list), Count=sapply(data_list, nrow))
  counts$Percentage <- round(counts$Count/sum(counts$Count)*100, 1)
  p <- ggplot(counts, aes(x=reorder(Tool,-Count), y=Count, fill=Tool)) +
    geom_bar(stat="identity") +
    geom_text(aes(label=sprintf("%d\n(%.1f%%)", Count, Percentage)), vjust=-0.5, size=4) +
    labs(title="各工具分类数目", x="工具", y="序列数") +
    theme_minimal() + theme(axis.text.x=element_text(angle=45, hjust=1), legend.position="none") +
    scale_fill_brewer(palette="Set2") + ylim(0, max(counts$Count)*1.15)
  ggsave(file.path(output_dir, "classification_counts.pdf"), p, width=max(8, nrow(counts)*1.2), height=6)
  return(counts)
}

plot_intersection_diagram <- function(data_list, output_dir) {
  if (length(data_list) < 2) return()
  n_tools <- length(data_list)
  id_lists <- lapply(data_list, function(x) unique(x$contig_id))

  # ── ≤3 工具: Venn ──
  if (n_tools <= 3) {
    log_msg("VIS", "绘制 Venn 图...")
    colors <- brewer.pal(min(n_tools, 8), "Set2")[1:n_tools]
    tryCatch({
      venn.plot <- venn.diagram(x=id_lists, filename=NULL, category.names=names(data_list),
        fill=colors, cex=1.2, cat.cex=1.2, cat.fontface="bold",
        main="序列交集 (Venn)")
      pdf(file.path(output_dir, "intersection_venn.pdf"), width=10, height=10)
      grid.draw(venn.plot); dev.off()
    }, error=function(e) log_msg("WARN", "Venn: %s", e$message))
    return()
  }

  # ── >3 工具: UpSet ──
  if (!requireNamespace("UpSetR", quietly=TRUE)) {
    log_msg("WARN", "UpSetR 未安装, 跳过 UpSet 图 (install.packages('UpSetR'))")
    return()
  }
  log_msg("VIS", sprintf("绘制 UpSet 图 (%d 工具)...", n_tools))

  # 构建二元矩阵
  all_ids <- unique(unlist(id_lists))
  bin_mat <- as.data.frame(lapply(names(id_lists), function(tn) {
    as.integer(all_ids %in% id_lists[[tn]])
  }))
  colnames(bin_mat) <- names(id_lists)

  tryCatch({
    pdf(file.path(output_dir, "intersection_upset.pdf"), width=14, height=8)
    suppressWarnings(UpSetR::upset(bin_mat, sets=names(id_lists), order.by="freq",
          main.bar.color="#4981BF", sets.bar.color="#66c2a5",
          nintersects=min(40, 2^n_tools), nsets=n_tools,
          text.scale=c(1.2,1.2,1,1,1.2,1.2)))
    dev.off()
  }, error=function(e) log_msg("WARN", "UpSet: %s", e$message))
}

plot_venn_by_tax_level <- function(data_list, output_dir) {
  if (length(data_list) > 7) {
    log_msg("INFO", ">7 tools, 跳过 per-rank Venn (使用 UpSet 替代)")
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
        labs(title=level) + theme(plot.title=element_text(hjust=0.5, size=12, face="bold"), legend.position="none")
      venn_plots[[level]] <- p
    }
  }
  if (length(venn_plots)>0) {
    combined <- wrap_plots(venn_plots, ncol=4, nrow=2) + plot_annotation(title="Per-Rank Intersection")
    ggsave(file.path(output_dir, "combined_tax_level_venn.pdf"), combined, width=20, height=12)
  }
}

analyze_consistency_optimized <- function(data_list, common_ids, output_dir) {
  if (length(common_ids)==0) return(NULL)
  log_msg("VIS", "分析一致性...")
  fll <- lapply(names(data_list), function(nm) {
    dt <- data_list[[nm]][contig_id %in% common_ids, .SD, .SDcols=c("contig_id",TAX_LEVELS)]; dt[, Tool:=nm]; dt
  })
  big_dt <- rbindlist(fll)
  results <- list()
  for (level in TAX_LEVELS) {
    vd <- big_dt[, .(contig_id, Tool, Taxon=get(level))][!is.na(Taxon) & Taxon!=""]
    if (nrow(vd)==0) next
    counts <- vd[, .N, by=.(contig_id, Taxon)]; setorder(counts, contig_id, -N, Taxon)
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
      theme_minimal() + scale_fill_brewer(palette="Set2") + scale_y_continuous(labels=percent, limits=c(0,1.1))
    ggsave(file.path(output_dir, "agreement_rates.pdf"), p, width=12, height=8)
    fwrite(ad, file.path(output_dir, "agreement_stats.tsv"), sep="\t")
  }
  return(ad)
}

analyze_detailed_consistency <- function(data_list, common_ids, output_dir) {
  if (length(common_ids)==0) return(NULL)
  log_msg("INFO", "详细一致性分析...")
  tool_names <- names(data_list)
  consistency_summary <- data.frame()
  for (level in TAX_LEVELS) {
    long_list <- lapply(tool_names, function(tn) {
      dt <- data_list[[tn]][contig_id %in% common_ids, .(contig_id, taxon=get(level))]; dt[, Tool:=tn]
      dt[taxon %in% c("","NA","na","N/A","no rank","undefined","unknown","null"), taxon := NA_character_]; dt
    })
    big_dt <- rbindlist(long_list); valid_dt <- big_dt[!is.na(taxon)]
    if (nrow(valid_dt)>0) {
      counts <- valid_dt[, .N, by=.(contig_id, taxon)]; setorder(counts, contig_id, -N, taxon)
      consensus_dt <- counts[, .SD[1], by=contig_id, .SDcols="taxon"]; setnames(consensus_dt, "taxon", "consensus_taxon")
      stats_dt <- valid_dt[, .(uniq_cnt=uniqueN(taxon), all_vals_str=paste(sort(unique(taxon)), collapse="|")), by=contig_id]
      res_dt <- merge(stats_dt, consensus_dt, by="contig_id", all=TRUE)
      res_dt[, `:=`(is_consistent=(uniq_cnt==1), consistency_status=ifelse(uniq_cnt==1,"consistent","inconsistent"))]
    } else {
      res_dt <- data.table(contig_id=common_ids, uniq_cnt=0, all_vals_str=NA_character_,
                           consensus_taxon=NA_character_, is_consistent=NA, consistency_status="all_NA")
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
      labs(title="各工具各等级分类数", x="等级", y="序列数") + theme_minimal() +
      theme(axis.text.x=element_text(angle=45, hjust=1)) + scale_fill_brewer(palette="Set2")
    p2 <- ggplot(consistency_summary, aes(x=Tool, y=Classified, fill=Level)) +
      geom_bar(stat="identity", position=position_stack()) +
      labs(title="分类数堆叠", x="工具", y="序列数") + theme_minimal() + scale_fill_brewer(palette="Set3")
    ggsave(file.path(output_dir, "classification_by_level.pdf"), p1, width=12, height=8)
    ggsave(file.path(output_dir, "classification_stacked.pdf"), p2, width=10, height=8)
    fwrite(consistency_summary, file.path(output_dir, "consistency_summary.tsv"), sep="\t")
  }
}

generate_grouped_consistent_files <- function(data_list, common_ids, output_dir) {
  if (length(common_ids)==0) return()
  log_msg("INFO", "按分类单元分组...")
  tool_names <- names(data_list)
  overall_summary <- data.frame()
  for (level in TAX_LEVELS) {
    long_list <- lapply(tool_names, function(tn) {
      dt <- data_list[[tn]][contig_id %in% common_ids, .(contig_id, taxon=get(level))]
      dt[taxon %in% c("","NA","na","N/A"), taxon:=NA_character_]; dt
    })
    valid_dt <- rbindlist(long_list)[!is.na(taxon)]
    if (nrow(valid_dt)==0) next
    check_dt <- valid_dt[, .(uniq_cnt=uniqueN(taxon)), by=contig_id]
    consistent_ids <- check_dt[uniq_cnt==1, contig_id]
    if (length(consistent_ids)==0) next
    consistent_data <- valid_dt[contig_id %in% consistent_ids, .SD[1], by=contig_id]
    setnames(consistent_data, "taxon", "consensus_taxon")
    level_dir <- file.path(output_dir, paste0("consistent_",tolower(level)))
    dir.create(level_dir, recursive=TRUE, showWarnings=FALSE)
    summary_dt <- consistent_data[, .(Count=.N), by=consensus_taxon]; setorder(summary_dt, -Count)
    fwrite(summary_dt, file.path(level_dir, "summary.tsv"), sep="\t")
    for (taxon_name in unique(consistent_data$consensus_taxon)) {
      safe <- gsub("[^[:alnum:]_]", "_", taxon_name); safe <- gsub("_{2,}","_", gsub("^_|_$","", safe))
      sub_dt <- consistent_data[consensus_taxon==taxon_name, .(contig_id)]
      fwrite(sub_dt, file.path(level_dir, paste0(safe,".tsv")), sep="\t")
    }
    overall_summary <- rbind(overall_summary, data.frame(
      Taxonomic_Level=level, Unique_Taxa=length(unique(consistent_data$consensus_taxon)),
      Total_Sequences=nrow(consistent_data)))
  }
  if (nrow(overall_summary)>0) {
    fwrite(overall_summary, file.path(output_dir, "consistent_summary.tsv"), sep="\t")
    p1 <- ggplot(overall_summary, aes(x=Taxonomic_Level, y=Unique_Taxa, fill=Taxonomic_Level)) +
      geom_bar(stat="identity") + geom_text(aes(label=Unique_Taxa), vjust=-0.5) +
      labs(title="一致分类单元数", x="等级") + theme_minimal() + theme(legend.position="none") + scale_fill_brewer(palette="Set2")
    p2 <- ggplot(overall_summary, aes(x=Taxonomic_Level, y=Total_Sequences, fill=Taxonomic_Level)) +
      geom_bar(stat="identity") + geom_text(aes(label=Total_Sequences), vjust=-0.5) +
      labs(title="一致序列数", x="等级") + theme_minimal() + theme(legend.position="none") + scale_fill_brewer(palette="Set3")
    ggsave(file.path(output_dir, "consistent_summary.pdf"), plot_grid(p1, p2, ncol=2, labels="AUTO"), width=14, height=7)
  }
}

# ==============================================================================
# 新增: Rank 清理 / 加权投票引擎
# ==============================================================================
KNOWN_REALMS   <- c("Riboviria","Monodnaviria","Duplodnaviria","Varidnaviria","Adnaviria","Ribozyviria")
SUBRANK_SUFFIX <- c("viricotina","viricetidae","virineae","virinae")

RANK_DEPTH_WEIGHTS <- c(Realm=1, Kingdom=2, Phylum=4, Class=8, Order=16, Family=32, Genus=64, Species=128)

TOOL_BIAS <- c(ACVirus=1.2, VITAP=1.1, mmseqs=1.0, metabuli=1.0,
               CAT=0.9, genomad=0.9, diamond_lca=0.8, vcontact3=0.7,
               contigtax=0.6, BASTA=0.6, PhaGCN3=0.8)

is_valid_value <- function(v) {
  if (length(v) != 1) return(FALSE)
  !is.na(v) && v != "" && !(v %in% c("-","NA","na","N/A","no rank","undefined","unknown","null","default","Unclassified"))
}

species_quality_score <- function(v) {
  if (!is_valid_value(v)) return(0)
  s <- as.character(v)
  if (grepl("(inae|idae|ales|icetes|viricota)$", s, ignore.case=TRUE)) return(0)
  score <- 1.0
  if (grepl("\\bsp\\.?\\b|\\bcf\\.?\\b|\\baff\\.?\\b", s, ignore.case=TRUE, perl=TRUE)) score <- score * 0.3
  if (grepl("unclassified|environmental|uncultured", s, ignore.case=TRUE)) score <- score * 0.1
  if (grepl("^unplaced|^novel_", s, ignore.case=TRUE)) score <- score * 0.2
  if (grepl("^[A-Z][a-z]+virus\\s+[a-z]", s, perl=TRUE)) score <- score * 1.5
  if (!grepl(" ", s)) score <- score * 0.5
  return(score)
}

clean_all_ranks <- function(data_list) {
  for (tn in names(data_list)) {
    dt <- data_list[[tn]]
    if (is.null(dt) || nrow(dt)==0) next

    # Realm: 必须在 KNOWN_REALMS 中
    dt[!Realm %in% KNOWN_REALMS, Realm := NA_character_]
    # Phylum: 必须以 -viricota 结尾
    dt[!is.na(Phylum) & !grepl("viricota$", Phylum), Phylum := NA_character_]
    # Class: 必须以 -viricetes 结尾
    dt[!is.na(Class) & !grepl("viricetes$", Class), Class := NA_character_]
    # Order: 必须以 -virales 结尾
    dt[!is.na(Order) & !grepl("virales$", Order), Order := NA_character_]
    # Family: 必须以 -viridae 结尾
    dt[!is.na(Family) & !grepl("viridae$", Family), Family := NA_character_]
    # Genus: 必须以 -virus 结尾, 但不能是 -viridae/-virinae
    dt[!is.na(Genus) & (!grepl("virus$", Genus) | grepl("viridae$|virinae$", Genus)), Genus := NA_character_]
    # Species: 不能以 subrank/高阶后缀结尾
    dt[!is.na(Species) & grepl("(inae|idae|ales|icetes|viricota)$", Species), Species := NA_character_]

    # 所有 rank: 清除 subrank 后缀值
    for (col in TAX_LEVELS) {
      for (sfx in SUBRANK_SUFFIX) {
        dt[grepl(paste0(sfx,"$"), get(col)), (col) := NA_character_]
      }
    }
    data_list[[tn]] <- dt
  }
  return(data_list)
}

compute_tool_weights <- function(data_list, consensus_stats=NULL) {
  weights <- list()
  for (tn in names(data_list)) {
    wvec <- numeric(length(TAX_LEVELS))
    names(wvec) <- TAX_LEVELS
    for (lvl in TAX_LEVELS) {
      w <- RANK_DEPTH_WEIGHTS[lvl] * ifelse(tn %in% names(TOOL_BIAS), TOOL_BIAS[tn], 0.8)
      if (!is.null(consensus_stats) && tn %in% names(consensus_stats)) {
        row <- consensus_stats[[tn]]
        rate <- row[[lvl]]
        if (!is.null(rate) && is.numeric(rate) && rate > 0) w <- w * rate
      }
      wvec[lvl] <- w
    }
    weights[[tn]] <- wvec
  }
  return(weights)
}

weighted_vote_rank <- function(cid, data_list, tool_weights, rank) {
  votes <- list()
  for (tn in names(data_list)) {
    dt <- data_list[[tn]]
    if (is.null(dt) || nrow(dt)==0) next
    row <- dt[contig_id==cid]
    if (nrow(row)==0) next
    v <- row[[rank]][1]
    if (!is_valid_value(v)) next
    w <- tool_weights[[tn]][rank]
    if (rank=="Species") w <- w * species_quality_score(v)
    if (is.null(votes[[v]])) votes[[v]] <- 0
    votes[[v]] <- votes[[v]] + w
  }
  if (length(votes)==0) return(NA_character_)
  names(votes)[which.max(unlist(votes))]
}

harmonize_genus_species <- function(result_dt) {
  sp_vec <- as.character(result_dt$Species)
  ge_vec <- as.character(result_dt$Genus)
  for (i in seq_len(nrow(result_dt))) {
    sp <- sp_vec[i]; ge <- ge_vec[i]
    if (is.na(sp) || is.na(ge) || sp=="" || ge=="") next
    parts <- strsplit(sp, " ")[[1]]
    if (length(parts) < 2) next
    genus_from_sp <- parts[1]
    if (!grepl("virus$", genus_from_sp) || grepl("viridae$|virinae$", genus_from_sp)) next
    if (any(sapply(SUBRANK_SUFFIX, function(s) grepl(paste0(s,"$"), genus_from_sp)))) next
    if (tolower(genus_from_sp) != tolower(ge)) {
      set(result_dt, i, "Genus", genus_from_sp)
    }
  }
  return(result_dt)
}

build_consensus <- function(data_list, tool_weights) {
  # ── 1. 去重 + 堆叠为 long format ──
  for (tn in names(data_list)) {
    data_list[[tn]] <- unique(data_list[[tn]], by="contig_id")
  }
  stacked <- rbindlist(lapply(names(data_list), function(tn) {
    d <- copy(data_list[[tn]]); d[, Tool := tn]
  }), use.names=TRUE, fill=TRUE)

  # 长格式: contig_id, Tool, Rank, Taxon
  long <- melt(stacked, id.vars=c("contig_id","Tool"),
               measure.vars=TAX_LEVELS, variable.name="Rank", value.name="Taxon")
  long <- long[!is.na(Taxon) & nchar(Taxon) > 0]

  if (nrow(long)==0) return(data.table())

  # ── 2. 权重表 (Tool × Rank) ──
  wt_dt <- rbindlist(lapply(names(tool_weights), function(tn) {
    data.table(Tool=tn, Rank=TAX_LEVELS, Weight=as.numeric(tool_weights[[tn]]))
  }))

  # ── 3. 合并权重 + Species 质量调整 ──
  long <- merge(long, wt_dt, by=c("Tool","Rank"), all.x=TRUE)
  long[is.na(Weight), Weight := 0]
  long[Rank=="Species", Weight := Weight * sapply(Taxon, species_quality_score)]

  # ── 4. 聚合投票: sum(Weight) per (contig, Rank, Taxon) ──
  votes <- long[, .(total_weight=sum(Weight)), by=.(contig_id, Rank, Taxon)]
  setorder(votes, contig_id, Rank, -total_weight)

  # ── 5. 每 (contig, Rank) 取票数最高者 ──
  winners <- votes[, .SD[1], by=.(contig_id, Rank)]

  # ── 6. Cast 回宽表 ──
  wide <- dcast(winners[, .(contig_id, Rank, Taxon)], contig_id ~ Rank, value.var="Taxon")
  for (lvl in setdiff(TAX_LEVELS, names(wide))) set(wide, NULL, lvl, NA_character_)

  # ── 7. 空缺填补: 从完备性最高的工具补 ──
  # 计算每个 tool-contig 的完备度
  best_per_contig <- stacked[, .(score=sum(!is.na(.SD) & .SD!=""), .SDcols=TAX_LEVELS),
                             by=.(contig_id, Tool)]
  setorder(best_per_contig, contig_id, -score)
  best <- best_per_contig[, .SD[1], by=contig_id]

  for (lvl in TAX_LEVELS) {
    missing <- wide[is.na(get(lvl)) | get(lvl)=="", which=TRUE]
    if (length(missing)==0) next
    cids_miss <- wide[missing, contig_id]
    fill_dt <- stacked[contig_id %in% cids_miss, .(contig_id, Tool, val=get(lvl))]
    fill_dt <- fill_dt[!is.na(val) & val!=""][, .SD[1], by=contig_id]
    wide[fill_dt, on="contig_id", (lvl) := fill_dt$val]
  }

  # ── 8. completeness + confidence ──
  wide[, completeness := rowSums(!is.na(.SD) & .SD!=""), .SDcols=TAX_LEVELS]
  wide[, confidence := round(completeness/length(TAX_LEVELS), 2)]

  # ── 9. 找 primary_tool (匹配度最高) ──
  best_tool <- stacked[wide[, .(contig_id)], on="contig_id", nomatch=NULL]
  match_score <- best_tool[, {
    sc <- 0
    wr <- wide[contig_id==.BY$contig_id]
    for (lvl in TAX_LEVELS) {
      tv <- get(lvl); fv <- wr[[lvl]]
      if (is_valid_value(tv) && is_valid_value(fv) &&
          tolower(as.character(tv))==tolower(as.character(fv))) {
        sc <- sc + RANK_DEPTH_WEIGHTS[lvl]
      }
    }
    .(match_score=sc)
  }, by=.(contig_id, Tool)]
  setorder(match_score, contig_id, -match_score)
  pt <- match_score[, .(primary_tool=Tool[1]), by=contig_id]
  wide <- merge(wide, pt, by="contig_id", all.x=TRUE)
  wide[is.na(primary_tool), primary_tool := "consensus"]

  # ── 10. genus-species 一致性 ──
  wide <- harmonize_genus_species(wide)

  # ── 11. 计算 agree 列 ──
  for (lvl in TAX_LEVELS) wide[, (paste0(lvl,"_agree")) := ""]
  agree_long <- long[, .(total=.N), by=.(contig_id, Rank)]
  for (i in 1:nrow(wide)) {
    cid <- wide[i, contig_id]
    for (lvl in TAX_LEVELS) {
      fv <- wide[i, get(lvl)]
      if (!is_valid_value(fv)) {
        tt <- agree_long[contig_id==cid & Rank==lvl, total]
        set(wide, i, paste0(lvl,"_agree"), sprintf("0/%d", if (length(tt)) tt[1] else 0L))
        next
      }
      agreed <- stacked[contig_id==cid & !is.na(get(lvl)) & get(lvl)!="" &
                        tolower(get(lvl))==tolower(fv), unique(Tool)]
      total <- stacked[contig_id==cid & !is.na(get(lvl)) & get(lvl)!="", uniqueN(Tool)]
      set(wide, i, paste0(lvl,"_agree"),
          sprintf("%d/%d: %s", length(agreed), total, paste(sort(agreed), collapse=",")))
    }
  }

  # ── 12. 列顺序 ──
  agree_cols <- paste0(TAX_LEVELS, "_agree")
  setcolorder(wide, c("contig_id","primary_tool","completeness","confidence",
                      TAX_LEVELS, agree_cols))
  return(wide[])
}

# ==============================================================================
# 主流程
# ==============================================================================
option_list <- list(
  make_option("--combined", type="character", default=NULL, help="virus_classifier.py 的 combined_taxonomy.tsv"),
  make_option("--mmseqs", type="character", default=NULL),
  make_option("--vcontact3", type="character", default=NULL),
  make_option("--vitap", type="character", default=NULL),
  make_option("--phagcn3", type="character", default=NULL),
  make_option("--acvirus", type="character", default=NULL),
  make_option("--genomad", type="character", default=NULL, help="genomad taxonomy.tsv"),
  make_option("--metabuli", type="character", default=NULL, help="metabuli taxonomy.tsv"),
  make_option("--cat", type="character", default=NULL, help="CAT taxonomy.tsv"),
  make_option("--diamond_lca", type="character", default=NULL),
  make_option("--contigtax", type="character", default=NULL),
  make_option("--basta", type="character", default=NULL),
  make_option(c("-o","--output"), type="character", default="."),
  make_option("--cores", type="integer", default=NULL)
)

opt <- parse_args(OptionParser(option_list=option_list))
if (!dir.exists(opt$output)) dir.create(opt$output, recursive=TRUE)
setup_threads(opt$cores)

# 优先从 combined_taxonomy.tsv 解析所有工具
all_data <- list()
if (!is.null(opt$combined) && file.exists(opt$combined)) {
  combined_result <- parse_combined_v2(opt$combined)
  if (!is.null(combined_result)) all_data <- combined_result
}

# 补充独立工具文件 (覆盖/追加)
parsers <- list(
  mmseqs=parse_mmseqs, vcontact3=parse_vcontact3, vitap=parse_vitap,
  phagcn3=parse_phagcn3, acvirus=parse_acvirus
)
for (tool in names(parsers)) {
  f <- opt[[tool]]
  if (!is.null(f) && file.exists(f)) {
    tryCatch({
      res <- parsers[[tool]](f)
      if (!is.null(res) && nrow(res)>0) all_data[[tool]] <- res
    }, error=function(e) log_msg("ERROR", "%s: %s", tool, e$message))
  }
}

if (length(all_data) < 1) stop("无有效数据, 检查 --combined 或 --mmseqs/--vitap 等参数。")
log_msg("INFO", "共加载 %d 个工具数据: %s", length(all_data), paste(names(all_data), collapse=", "))

# Phase 2: 清理 rank 值 (subrank/不规范后缀)
log_msg("INFO", "清理 rank 值 (subrank/后缀校验)...")
all_data <- clean_all_ranks(all_data)

# 输出标准化表
for (tn in names(all_data)) {
  fwrite(unique(all_data[[tn]], by="contig_id"), file.path(opt$output, paste0("standardized_", tn, ".tsv")), sep="\t", na="NA")
}

# 交叉分析
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

# Phase 5: 计算工具权重 (有一致率数据时用数据驱动)
log_msg("INFO", "计算工具权重...")
ids_for_weight <- if (length(common_ids) > 0) common_ids else
  mat_dt[count >= max(2, length(all_data)/2), contig_id]
consensus_stats <- NULL
if (length(ids_for_weight) > 0) {
  stats_list <- lapply(names(all_data), function(tn) {
    dt <- all_data[[tn]][contig_id %in% ids_for_weight]
    row_stats <- list()
    for (lvl in TAX_LEVELS) {
      valid_cnt <- sum(sapply(dt[[lvl]], is_valid_value))
      if (valid_cnt == 0) { row_stats[[lvl]] <- 0; next }
      # 计算共识值 (所有工具该 rank 的众数)
      agreed <- 0; total <- 0
      for (cid in dt$contig_id) {
        vals <- c()
        for (t2 in names(all_data)) {
          r2 <- all_data[[t2]][contig_id==cid]
          if (nrow(r2)>0 && is_valid_value(r2[[lvl]])) vals <- c(vals, r2[[lvl]])
        }
        if (length(vals)==0) next
        cons_val <- names(sort(table(vals), decreasing=TRUE))[1]
        tool_val <- all_data[[tn]][contig_id==cid][[lvl]]
        total <- total + 1
        if (!is.na(tool_val) && tool_val==cons_val) agreed <- agreed + 1
      }
      row_stats[[lvl]] <- if (total>0) agreed/total else 0
    }
    row_stats
  })
  names(stats_list) <- names(all_data)
  consensus_stats <- stats_list
  log_msg("INFO", "基于 %d 个序列计算一致率权重", length(ids_for_weight))
} else {
  log_msg("INFO", "无足够交集序列, 使用默认权重")
}

tool_weights <- compute_tool_weights(all_data, consensus_stats)

# Phase 6: 逐 rank 加权投票构建共识
log_msg("INFO", "逐 rank 加权投票...")
final_result <- build_consensus(all_data, tool_weights)

agree_cols <- paste0(TAX_LEVELS, "_agree")
out_cols <- c("contig_id","primary_tool","completeness","confidence", TAX_LEVELS, agree_cols)
fwrite(final_result[, out_cols, with=FALSE], file.path(opt$output, "final_integrated_classification.tsv"), sep="\t", na="NA")
log_msg("INFO", "最终分类: %d 个序列", nrow(final_result))

# 填充统计
fill_counts <- sapply(TAX_LEVELS, function(l) sum(sapply(final_result[[l]], is_valid_value)))
fill_rates  <- round(fill_counts / nrow(final_result), 4)
fill_dt <- data.table(Rank=TAX_LEVELS, Filled=fill_counts, Rate=fill_rates)
fwrite(fill_dt, file.path(opt$output, "fill_stats.tsv"), sep="\t")

# 共识汇总图
if (requireNamespace("gridExtra", quietly=TRUE)) {
  log_msg("VIS", "绘制共识汇总...")
  p1 <- ggplot(final_result, aes(x=factor(completeness))) +
    geom_bar(fill="#66c2a5") + labs(title="共识完备度分布", x="完备度", y="序列数") + theme_minimal()
  fr_df <- data.frame(Level=factor(TAX_LEVELS, levels=TAX_LEVELS), Rate=fill_rates*100)
  p2 <- ggplot(fr_df, aes(x=Level, y=Rate, fill=Level)) +
    geom_bar(stat="identity") + labs(title="各 Rank 填充率", x="等级", y="填充率(%)") +
    theme_minimal() + theme(legend.position="none") + scale_fill_brewer(palette="Set2")
  ggsave(file.path(opt$output, "consensus_summary.pdf"),
         plot_grid(p1, p2, ncol=2, labels="AUTO"), width=14, height=7)
}

# 报告
sink(file.path(opt$output, "analysis_summary.txt"))
cat("=== 病毒分类分析报告 ===\n\n")
cat(sprintf("时间: %s\n", Sys.time()))
cat(sprintf("总序列: %d\n", n_total))
cat(sprintf("最终分类: %d\n", nrow(final_result)))
cat(sprintf("工具数: %d\n", length(all_data)))
cat("\n各工具分类数:\n"); print(count_stats)
cat("\n各 Rank 填充率:\n")
for (i in 1:length(TAX_LEVELS))
  cat(paste0("  ", TAX_LEVELS[i], ": ", fill_counts[i],
             " (", round(fill_rates[i]*100, 1), "%)\n"))
cat("\n主要工具分布:\n"); print(table(final_result$primary_tool))
sink()
log_msg("INFO", "完成! 结果: %s", opt$output)
