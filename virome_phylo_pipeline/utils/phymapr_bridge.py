#!/usr/bin/env python3
"""
phymapr_bridge.py — phymapr R 桥 (TreeTime 路线系统地理)
========================================================

2026-10-12 新增, 参考 ADHS-Taylor/phymap-workflow (scripts/run_phymapr.R)。
与既有 BEAST 路线 (phylogeo→geo_paths) **并列**的第二条时间-地理路线:

    treetime stage (utils/treetime_bridge.py) → timetree.nexus
        └→ 本模块 → phymapr::generate_phylo_transmission()
              ├─ 自研 BM-ASR (lat/lon 连续祖先重建 → 吸附最近采样点)
              ├─ tree_plot.png              (ggtree 时间树按地点着色)
              ├─ transmission_map.gif       (gganimate 时间揭示动画)
              ├─ transmission_map_static.png(workflow 漏存, 本桥补上)
              └─ episodes.csv / pathways.csv(与 geo_paths 产物同构, 可双方法对照)

离线设计 (服务器无外网友好):
  · 元数据桥 build_phymapr_metadata() 预先解析坐标写入 lat/long 列 ——
    phymapr 检测到这两列即跳过 OSM 在线地理编码 (data_prep.R L80-82);
  · PHYMAPR_SRC 环境变量指向本地 phymapr 源码目录 (含 R/ 子目录) 时直接
    source() 加载, 不联网; 否则已装则 library, 最后才尝试 GitHub 安装。
    (workflow 原版硬编码 C:/Code/phymapr 等开发机路径的坑不搬。)

R 调用范式与 utils/virphy_bridge.py::run_treedater_ltt 一致:
模块级 R 脚本常量 + NamedTemporaryFile + subprocess + JSON 哨兵 +
永不抛异常 (失败一律 {"success": False, "error": ...})。
"""

import json
import os
import subprocess
import tempfile
from datetime import datetime
from typing import Dict, Optional

from utils.virphy_bridge import LogCollector

_JSON_SENTINEL = "___PHYMAPR_JSON_END___"

# phymapr 元数据桥的固定列 (generate_phylo_transmission 以参数指定列名, 不必叫
# accessionVersion/sampleCollectionDate —— 直接用管线三通道表的 name/date/location)
_PHYMAPR_META_COLS = ["name", "date", "location", "lat", "long"]


def _mid_date(raw: str) -> str:
    """日期精度中点规范化 (phymapr 侧 lubridate::as_date 只认完整 ISO 日期)。

    YYYY → YYYY-07-02; YYYY-MM → 该月中日; YYYY-MM-DD → 原样;
    其余原样返回 (R 侧 as.Date 解析失败 → NA → 该 tip 被剔除, 不静默造日期)。
    """
    import calendar

    s = str(raw).strip().rstrip("Z")
    parts = s.split("-")
    if len(parts) == 1 and parts[0].isdigit() and len(parts[0]) == 4:
        return f"{parts[0]}-07-02"
    if len(parts) == 2 and parts[0].isdigit() and parts[1].isdigit():
        y, mo = int(parts[0]), int(parts[1])
        if 1 <= mo <= 12:
            return f"{y:04d}-{mo:02d}-{calendar.monthrange(y, mo)[1] // 2 + 1:02d}"
    return s


def _resolve_loc(loc: str) -> Optional[tuple]:
    """coord_lookup 前置变体梯子。

    coord_lookup (geo_resolver L2 段拆分) 只认逗号/空格分隔的多级地名,
    对中文省市连写 ('宁夏银川' / '内蒙古自治区呼和浩特') 会整段未命中。
    这里补两级变体: ①按已知省名最长前缀拆出城市段重试 (城市级坐标优于省级);
    ②"省 + 空格 + 城" 拼接重试。命中语义与原串一致 (city > province 优先,
    由 lookup 内部表序保证); 全部落空返回 None, 绝不伪造坐标。
    """
    from utils.coord_lookup import lookup as _lk

    hit = _lk(loc)
    if hit:
        return hit
    try:
        from utils.geo_resolver import PROVINCE_LOOKUP
    except ImportError:  # noqa: BLE001 — geo_resolver 不可用时退化为单次 lookup
        return None
    prov_fallback = None
    for prov in sorted(PROVINCE_LOOKUP, key=len, reverse=True):
        if loc.startswith(prov) and len(loc) > len(prov):
            rest = loc[len(prov):].strip()
            for cand in (rest, f"{prov} {rest}"):
                hit = _lk(cand)
                if hit:
                    return hit
            # 城市段全落空 → 记省级坐标作最后兜底 (先走完下面的通用梯子,
            # 全称截尾/末段变体仍可能解出城市级坐标, 城市级优于省级)
            prov_fallback = _lk(prov)
            break
    # 全称前缀兜底 ('内蒙古自治区呼和浩特'): 省表是简称, 逐位截尾试城市段 ——
    # 只会在 city 精确命中时返回, 不会误配。
    if len(loc) >= 5:
        for i in range(2, len(loc) - 1):
            hit = _lk(loc[i:])
            if hit:
                return hit
    # 末段/去"市"后缀兜底 ('黑龙江省哈尔滨市' / 'USA California Los Angeles'):
    # 取逗号/空格分隔的末段, 以及对其去尾部行政后缀, 再试。
    import re as _re

    tail = _re.split(r"[,，\s]+", loc.strip())[-1] if loc.strip() else ""
    for cand in (tail, tail.rstrip("市县区")):
        if cand and cand != loc:
            hit = _lk(cand)
            if hit:
                return hit
    return prov_fallback


def build_phymapr_metadata(
    meta_csv: str,
    out_tsv: str,
    log: Optional[LogCollector] = None,
) -> Dict:
    """把管线 sample_metadata (name,date,location) 桥成 phymapr 输入 TSV。

    坐标走 utils/coord_lookup (管线唯一坐标源, 含 geo_cache 持久缓存) ——
    解析出的 lat/long 直接写进 TSV, phymapr 据此完全离线。
    未解析坐标/无日期的行**保留但留空** —— 由 R 侧统一 drop.tip,
    这样树-元数据名称对不齐的样本在 summary 里可见而不是静默消失。

    Returns
    -------
    dict: {success, out_tsv, n_rows, n_coord_ok, n_date_ok,
           missing_locations[...], error}
    """
    log = log or LogCollector()
    result: Dict = {"success": False, "out_tsv": out_tsv,
                    "n_rows": 0, "n_coord_ok": 0, "n_date_ok": 0,
                    "missing_locations": [], "error": None}
    try:
        import pandas as pd

        df = pd.read_csv(meta_csv, encoding="utf-8-sig")
        cols = {c.lower().strip(): c for c in df.columns}
        for need in ("name", "date", "location"):
            if need not in cols:
                result["error"] = f"元数据缺列 {need}: {meta_csv}"
                log.error(result["error"])
                return result
        meta = df[[cols["name"], cols["date"], cols["location"]]].copy()
        meta.columns = ["name", "date", "location"]
        meta["name"] = meta["name"].astype(str).str.strip()
        meta["location"] = meta["location"].fillna("").astype(str).str.strip()

        # 坐标按地点字符串去重解析 (_resolve_loc 含省市连写拆分; coord_lookup 自带
        # 磁盘缓存, 批量后 flush 落盘)
        coord_cache: Dict[str, Optional[tuple]] = {}
        for loc in sorted(set(meta["location"]) - {""}):
            try:
                hit = _resolve_loc(loc)
            except Exception:  # noqa: BLE001
                hit = None
            coord_cache[loc] = hit
        try:
            from utils.coord_lookup import flush as _coord_flush
            _coord_flush()
        except Exception:  # noqa: BLE001
            pass

        lats, lons = [], []
        for loc in meta["location"]:
            hit = coord_cache.get(loc)
            lats.append(hit[0] if hit else "")
            lons.append(hit[1] if hit else "")
        meta["lat"] = lats
        meta["long"] = lons
        meta["date"] = meta["date"].apply(_mid_date)

        meta = meta[_PHYMAPR_META_COLS]
        os.makedirs(os.path.dirname(os.path.abspath(out_tsv)), exist_ok=True)
        meta.to_csv(out_tsv, sep="\t", index=False)

        result.update(
            success=True,
            n_rows=len(meta),
            n_coord_ok=int((meta["lat"] != "").sum()),
            n_date_ok=int(meta["date"].str.fullmatch(r"\d{4}-\d{2}-\d{2}").sum()),
            missing_locations=sorted(l for l, h in coord_cache.items() if h is None),
        )
        log.emit(f"phymapr 元数据桥: {result['n_rows']} 行 "
                 f"(坐标 {result['n_coord_ok']}, 日期 {result['n_date_ok']}); "
                 f"未解析地点 {len(result['missing_locations'])}")
        return result
    except Exception as e:  # noqa: BLE001
        result["error"] = f"phymapr 元数据桥失败: {e}"
        log.error(result["error"])
        return result


# ═══════════════════════════════════════════════════════════════════
# R 侧脚本 (位置参数, 供 NamedTemporaryFile + Rscript --vanilla 调用)
# ═══════════════════════════════════════════════════════════════════
PHYMAPR_R_SCRIPT = r"""
## phymapr_bridge R 侧 — 依赖探测 → 树-元数据对齐 → generate_phylo_transmission → 落盘
## 输出指标经 ___PHYMAPR_JSON_END___ 哨兵回传 Python; 失败 quit(status>0) 且仍发 JSON。
args <- commandArgs(trailingOnly=TRUE)
timetree_file <- args[1]
meta_file     <- args[2]
out_dir       <- args[3]
map_style     <- args[4]
n_prune_early <- as.integer(args[5])
thr_direct    <- as.numeric(args[6])
thr_indirect  <- as.numeric(args[7])
distant_style <- args[8]
pace          <- as.numeric(args[9])

out <- list(success = FALSE)

emit_json <- function(status) {
  cat(jsonlite::toJSON(out, auto_unbox = TRUE, digits = 10), "\n")
  cat("___PHYMAPR_JSON_END___\n")
  quit(status = status)
}

## ── 0) jsonlite 单独先探 (发 JSON 全靠它) ──
if (!requireNamespace("jsonlite", quietly = TRUE)) {
  cat('{"success":false,"error":"jsonlite not installed (install.packages(jsonlite))"}\n')
  cat("___PHYMAPR_JSON_END___\n")
  quit(status = 3)
}

## ── 1) 依赖探测 ──
## hard: polygon 模式 + phymapr 内部代码路径必需 (tidygeocoder 在 DESCRIPTION
##       Imports 里, library(phymapr) 会强检; 本桥预写 lat/long 后运行时其实不联网)
## soft: 仅 library/install 模式需要 (source 模式不 library, tile 瓦片样式才用 sf/ggspatial)
hard <- c("ape","treeio","ggtree","ggplot2","dplyr","tidyr","tibble",
          "gganimate","gifski","lubridate","maps","tidygeocoder")
soft <- c("sf","ggspatial","phytools")
src_env <- Sys.getenv("PHYMAPR_SRC", "")
src_mode <- nzchar(src_env) && dir.exists(file.path(src_env, "R"))
missing <- hard[!vapply(hard, function(p) requireNamespace(p, quietly = TRUE), logical(1))]
if (!src_mode)
  missing <- union(missing, soft[!vapply(soft, function(p) requireNamespace(p, quietly = TRUE), logical(1))])
if (length(missing) > 0) {
  out$error <- paste0("Missing R packages: ", paste(missing, collapse = ", "))
  out$missing_packages <- missing
  out$install_hint <- c(
    "CRAN: install.packages(c('ape','ggplot2','dplyr','tidyr','tibble','gganimate','gifski','lubridate','maps','jsonlite','tidygeocoder','remotes','sf','ggspatial','phytools'))",
    "Bioconductor: BiocManager::install(c('ggtree','treeio'))",
    "phymapr: remotes::install_github('ADHS-Taylor/phymapr')  |  Offline: set env PHYMAPR_SRC=<local phymapr source dir containing R/>")
  emit_json(3)
}

## ── 2) 加载 phymapr (source > library > GitHub 安装) ──
phymapr_mode <- ""
if (src_mode) {
  for (f in list.files(file.path(src_env, "R"), pattern = "\\.[Rr]$", full.names = TRUE)) source(f)
  phymapr_mode <- paste0("source:", src_env)
} else if (requireNamespace("phymapr", quietly = TRUE)) {
  suppressMessages(library(phymapr))
  phymapr_mode <- "library"
} else {
  inst <- try({
    if (!requireNamespace("remotes", quietly = TRUE))
      install.packages("remotes", repos = "https://cloud.r-project.org")
    remotes::install_github("ADHS-Taylor/phymapr", upgrade = "never")
    suppressMessages(library(phymapr))
  }, silent = TRUE)
  if (inherits(inst, "try-error")) {
    out$error <- paste0("phymapr unavailable: GitHub install failed (offline?) and PHYMAPR_SRC not set. ",
                        conditionMessage(attr(inst, "condition")))
    emit_json(3)
  }
  phymapr_mode <- "github-install"
}
out$phymapr_mode <- phymapr_mode

## ── 3) 读元数据, 校验列 ──
meta <- read.delim(meta_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE,
                   colClasses = "character")
req_cols <- c("name","date","location","lat","long")
miss_cols <- setdiff(req_cols, colnames(meta))
if (length(miss_cols) > 0) {
  out$error <- paste0("metadata missing columns: ", paste(miss_cols, collapse = ", "))
  emit_json(2)
}
meta$lat  <- suppressWarnings(as.numeric(meta$lat))
meta$long <- suppressWarnings(as.numeric(meta$long))
meta$date_ok  <- !is.na(as.Date(meta$date, format = "%Y-%m-%d"))
meta$coord_ok <- !is.na(meta$lat) & !is.na(meta$long)

## ── 4) 树-元数据对齐: 无元数据/无日期/无坐标 tip 一并 drop ──
td  <- treeio::read.beast(timetree_file)
ptree <- ape::as.phylo(td)
tip_lab <- ptree$tip.label
m <- match(tip_lab, meta$name)
keep_tip <- !is.na(m) & meta$date_ok[m] & meta$coord_ok[m]
out$n_tips_total  <- length(tip_lab)
out$n_tips_used   <- sum(keep_tip)
out$n_tips_dropped <- length(tip_lab) - sum(keep_tip)
if (sum(keep_tip) < 5) {
  out$error <- sprintf("Usable tips = %d (<5): need name-matched tips with date+coords.", sum(keep_tip))
  emit_json(2)
}
if (any(!keep_tip)) td <- treeio::drop.tip(td, tip_lab[!keep_tip])

td_path   <- file.path(out_dir, "timetree_matched.nexus")
treeio::write.beast(td, td_path)
ptree2 <- ape::as.phylo(td)
meta2 <- meta[match(ptree2$tip.label, meta$name), req_cols]
meta_path <- file.path(out_dir, "phymapr_metadata_used.tsv")
write.table(meta2, meta_path, sep = "\t", quote = FALSE, row.names = FALSE)

## ── 5) generate_phylo_transmission ──
res <- try(generate_phylo_transmission(
  timetree_file = td_path,
  metadata_file = meta_path,
  col_specimen  = "name",
  col_date      = "date",
  col_location  = "location",
  map_style               = map_style,
  n_prune_early           = n_prune_early,
  threshold_direct_snp    = thr_direct,
  threshold_indirect_snp  = thr_indirect,
  distant_rendering_style = distant_style,
  animation_pace_balance  = pace), silent = TRUE)
if (inherits(res, "try-error")) {
  out$error <- paste0("generate_phylo_transmission failed: ",
                      conditionMessage(attr(res, "condition")))
  emit_json(4)
}

## ── 6) 落盘 (单件失败只降级, 不整体失败) ──
save_one <- function(expr) inherits(try(expr, silent = TRUE), "try-error")
ggsave(file.path(out_dir, "tree_plot.png"), res$tree_plot, width = 14, height = 9, dpi = 150)
if (!is.null(res$static_map))
  ggsave(file.path(out_dir, "transmission_map_static.png"), res$static_map,
         width = 12, height = 9, dpi = 150)
if (!is.null(res$map_animation)) {
  g <- try(animate(res$map_animation, width = 1200, height = 900, res = 120,
                   fps = 10, renderer = gifski_renderer()), silent = TRUE)
  if (!inherits(g, "try-error"))
    anim_save(file.path(out_dir, "transmission_map.gif"), animation = g)
}
if (!is.null(res$episodes))
  write.csv(res$episodes, file.path(out_dir, "episodes.csv"), row.names = FALSE)
if (!is.null(res$pathways))
  write.csv(res$pathways, file.path(out_dir, "pathways.csv"), row.names = FALSE)

## ── 7) 指标回传 ──
out$success    <- TRUE
out$n_pathways <- if (!is.null(res$pathways)) nrow(res$pathways) else 0
out$n_episodes <- if (!is.null(res$episodes)) nrow(res$episodes) else 0
out$n_locations <- length(unique(meta2$location))
if (!is.null(res$pathways) && nrow(res$pathways) > 0 && "transmission_type" %in% names(res$pathways))
  out$transmission_counts <- as.list(table(res$pathways$transmission_type))
emit_json(0)
"""


def run_phymapr(
    timetree_nexus: str,
    meta_tsv: str,
    output_dir: str,
    rscript_path: str = "Rscript",
    map_style: str = "polygon",
    n_prune_early: int = 4,
    direct_snp: float = 2.0,
    indirect_snp: float = 5.0,
    distant_style: str = "import_arrow",
    pace: float = 0.5,
    timeout: int = 3600,
    log: Optional[LogCollector] = None,
) -> Dict:
    """调 R 跑 phymapr 祖先重建 + 传播地图, 落盘全套产物。

    Returns
    -------
    dict: {success, outputs{tree_plot,gif,static_map,episodes_csv,pathways_csv,
          metadata_used,timetree_matched}, metrics{n_tips_*,n_pathways,...},
          summary_json, error}
          永不抛异常; 成功时写 phymapr_summary.json (stage checkpoint 标志物)。
    """
    log = log or LogCollector()
    os.makedirs(output_dir, exist_ok=True)
    result: Dict = {"success": False, "outputs": {}, "metrics": {},
                    "summary_json": None, "error": None}

    for label, path in (("时间树", timetree_nexus), ("元数据", meta_tsv)):
        if not path or not os.path.exists(path):
            result["error"] = f"输入缺失: {label} {path}"
            log.error(result["error"])
            return result

    with tempfile.NamedTemporaryFile(mode="w", suffix=".R", delete=False,
                                     encoding="utf-8") as f:
        f.write(PHYMAPR_R_SCRIPT)
        r_script_path = f.name

    cmd = [
        rscript_path, "--vanilla", r_script_path,
        timetree_nexus, meta_tsv, output_dir,
        str(map_style), str(int(n_prune_early)),
        str(float(direct_snp)), str(float(indirect_snp)),
        str(distant_style), str(float(pace)),
    ]
    log.emit(f"Running: {' '.join(cmd)}")
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True,
                              encoding="utf-8", errors="replace", timeout=timeout)
    except subprocess.TimeoutExpired:
        result["error"] = f"phymapr R 超时 (>{timeout}s; GIF 渲染大样本时偏慢, 可重跑续传)"
        log.error(result["error"])
        return result
    except Exception as e:  # noqa: BLE001
        result["error"] = f"Rscript 启动失败: {e}"
        log.error(result["error"])
        return result
    finally:
        try:
            os.unlink(r_script_path)
        except OSError:
            pass

    # ── 解析 JSON 哨兵 (与 TreeDater 同款: 取哨兵前最后一个 {..} 行) ──
    metrics: Dict = {}
    stdout = proc.stdout or ""
    if _JSON_SENTINEL in stdout:
        blob = stdout.split(_JSON_SENTINEL)[0].strip()
        for line in reversed(blob.splitlines()):
            line = line.strip()
            if line.startswith("{") and line.endswith("}"):
                try:
                    metrics = json.loads(line)
                except json.JSONDecodeError as e:
                    log.warning(f"phymapr 指标 JSON 解析失败 (非阻断): {e}")
                break

    if metrics.get("missing_packages"):
        hint = "; ".join(metrics.get("install_hint") or [])
        suffix = ("" if "PHYMAPR_SRC" in hint
                  else " | 或设 PHYMAPR_SRC 环境变量指向本地 phymapr 源码目录(含 R/)离线加载")
        result["error"] = f"R 依赖缺失: {', '.join(metrics['missing_packages'])}。安装: {hint}{suffix}"
        log.error(result["error"])
        result["metrics"] = metrics
        return result
    if proc.returncode != 0 or not metrics.get("success"):
        tail = " | ".join((proc.stderr or "").strip().splitlines()[-8:])
        result["error"] = (metrics.get("error")
                           or f"phymapr R 失败 (exit {proc.returncode}): {tail}")
        result["metrics"] = metrics
        log.error(result["error"])
        return result

    outputs = {
        "tree_plot": os.path.join(output_dir, "tree_plot.png"),
        "static_map": os.path.join(output_dir, "transmission_map_static.png"),
        "gif": os.path.join(output_dir, "transmission_map.gif"),
        "episodes_csv": os.path.join(output_dir, "episodes.csv"),
        "pathways_csv": os.path.join(output_dir, "pathways.csv"),
        "metadata_used": os.path.join(output_dir, "phymapr_metadata_used.tsv"),
        "timetree_matched": os.path.join(output_dir, "timetree_matched.nexus"),
    }
    missing_out = [k for k, v in outputs.items() if not os.path.exists(v)]
    if missing_out:
        log.warning(f"phymapr 产物缺件 (降级继续): {missing_out}")

    result["success"] = True
    result["outputs"] = {k: v for k, v in outputs.items() if os.path.exists(v)}
    result["metrics"] = metrics

    # checkpoint 标志物: phymapr_summary.json
    summary_path = os.path.join(output_dir, "phymapr_summary.json")
    summary = {"stage": "phymapr", "success": True,
               "generated_at": datetime.now().isoformat(timespec="seconds"),
               "metrics": metrics, "outputs": result["outputs"]}
    try:
        with open(summary_path, "w", encoding="utf-8") as f:
            json.dump(summary, f, ensure_ascii=False, indent=2)
        result["summary_json"] = summary_path
    except OSError as e:
        log.warning(f"phymapr_summary.json 写入失败: {e}")

    log.emit(f"phymapr 完成: pathways={metrics.get('n_pathways')} "
             f"episodes={metrics.get('n_episodes')} "
             f"locations={metrics.get('n_locations')} "
             f"(tips {metrics.get('n_tips_used')}/{metrics.get('n_tips_total')})")
    return result
