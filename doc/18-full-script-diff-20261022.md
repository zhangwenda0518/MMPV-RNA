# 18 · 全仓脚本比对：服务器 vs 本仓（2026-10-22）

> 缘由：用户担心服务器上正在跑的生产脚本与本地 `MMPV-RNA`（含本轮 eve-mask 改动）冲突，要求「全部比较，所有的脚本」（m00755）。
> 本文是逐文件比对结果、谁更新的判定与处置建议。
> 前一份 `doc/17-server-vs-local-audit-20261022.md` 只比了 3 个文件，且把「服务器 = 旧版」的判断写成了「服务器独有」，
> 已在本文 §4.4 与 doc/17 的修订中更正。

## 0. 结论摘要

- **规模**：服务器 `~/MMPV-RNA` 有 **1,491** 个脚本/配置文件，本仓 **2,058** 个（扩展名见 §1）。
- **不一致构成**（同路径 367 个）：同内容 **217**、服务器停在本仓历史旧版 **107**、两边都改过 **43**；另有仅服务器 **1,124**、仅本地 **1,691**。
- **与 eve-mask 改动不冲突**：同路径文件只有 `public_metadata_pipeline/build_host_pipeline.py`，服务器那份正是本仓 `24b32f7` 旧版且服务器侧无改动（`58eb7c7d`，16,029 B CRLF）；
  掩码相关新文件（`mask_host_genome.py`、两个测试）服务器都不存在；服务器全仓 `grep mask` 也没有任何既有掩码实现。
- **服务器不是「另一个分支」，是 `24b32f7`（2026-09-28）的快照**：本地 HEAD 之后又走了 37 个提交，服务器工作区则在 10-09～10-19 之间被单独改过。
  所以两边**各有对方没有的修改**，不是单向落后。
- **真正「服务器在跑、本仓没有」的代码只有两条已归档管线的旧版**：`endogenous_virus_pipeline`（v6.3，与本地 archive 版有 26 个文件内容不同）
  与 `host_expression_pipeline`（10 同 / 23 不同）。本地 `archive/20261011_*` 是更新版本（v7.x / 2026-10-11 归档），**不是缺失**。
- **已确认的唯一实质缺口**：发现管线 2026-10-19 的「丢弃审计」补丁（`run_host_prediction.py`、`virus_classifier.py` 输出 `dropped_ids.tsv`
  记录被静默丢弃的序列）本仓没有。其余 41 个 UNIQUE 文件是两条独立改线的差异，逐个在 §4 列出，需人工确认后决定是否回移。
  **→ 已于 2026-10-26 回移（含本文漏记的 R 侧第三处），见 §0.2。**
- 服务器独有 1,124 个脚本中 **643 个只是换了路径**（内容与本地某处完全相同，多为本地已归档的旧副本），真正独有内容的 481 个里 307 个是服务器自己的
  `archive/20261018_server_cleanup/` 清理包，只有 **165 个**属「活代码」（§5.2）。

### 0.1 处置建议（按紧急度）

| # | 动作 | 涉及 | 理由 |
|---:|---|---|---|
| 1 | **先冻结服务器工作区**：`git checkout -b server-wip-20261022 && git add -A && git commit` | 全仓 3,770 条未提交变更 | 不冻结就无法安全同步/回滚 |
| 2 | **确认 1,348 样本用的 `host_depletion.py` 版本** | `grep -rl "已启用，自动跳过 align" <输出根>/logs` | 服务器 7 月版会无条件跳过 align，直接影响稿件 Methods（doc/17 §3） |
| 3 | **回移 2026-10-19 丢弃审计** | `virome_discovery_pipeline/{run_host_prediction.py,virus_classifier.py}` | 本仓缺 `dropped_ids.tsv` 逻辑；曾有 74,943 条序列静默消失无据可查 |
| 4 | **逐个判定 43 个 UNIQUE 文件** | §4 | 两边各有对方没有的修改，不能整树覆盖 |
| 5 | **以本地为准覆盖 107 个 STALE** | §3 | 服务器内容就是本仓 `24b32f7` 的副本，无独立内容 |
| 6 | **忽略**归档、一次性调试、结果 JSON | §5.1 / §6 | 内容本地已有或本就无用 |

### 0.2 回移记录（2026-10-26，上表 #3 已完成）

服务器侧该补丁共 **3 处**，本文原只记了 2 个 py，R 侧 `virus_classifier_analysis.R` 的共识丢弃审计同样缺、当时未记。

| 文件 | 回移内容 | 相对服务器版的修正 |
|---|---|---|
| `virome_discovery_pipeline/run_host_prediction.py` | `run_ensemble` 内新增输入全集/输出全集差集审计 → `dropped_ids.tsv`（`lost_in_merge` / `absent_from_taxonomy`）+ `duplicate_tax_ids.tsv` | 原样（`-i` 已是 centroids FASTA，天然可用） |
| `virome_discovery_pipeline/virus_classifier.py` | `_tool_ids` 逐工具序列名收集 + `_write_classifier_drop_log()` → `{sample}_dropped_ids.tsv`（`no_valid_lineage` / `no_tool_hit`，含 `tools_that_saw_it`） | **服务器版在管线里是死代码**：它只在 `.classed/` 下 glob `{sample}*.fasta`，而输入 FASTA 在 `centroids/final_centroids.fasta`，glob 恒空 → 审计文件永不生成。改为显式传 `input_fasta`，保留 glob 兜底 |
| `virome_discovery_pipeline/virus_classifier_analysis.R` | `input_all_ids`（输入全集）+ 共识阶段 `dropped_ids.tsv`（`no_consensus` / `no_tool_hit`）；新增 `--fasta` | **服务器版读 header 用 `system("grep '^>' …")`**，Windows 下 cmd 不认单引号，pattern 变字面 `'^>'` → 恒不匹配，静默退化为「工具命中全集」。改为纯 R 分块 `readLines` + `grep("^>", …)` |
| `virome_discovery_pipeline/virome_pipeline.py` | R 调用补 `--fasta {centroids}` | 服务器版无此参数，R 侧只能靠同级探测（探不到就退化） |

另**顺带修回** `report_pipeline.py` 的 `filter_summary.tsv`：本地在 10-20 后的重写回退了服务器 10-19 的修复（失效符号链接/0 候选时整表为空）。
恢复为：表头含 `Note` 列、`is_file() or is_symlink()` 判定、`BROKEN_LINK` / `ZERO_CANDIDATES` 标注、`Retained(%)` 允许 `NA`、**始终写行**；其余本地新功能（s11 章节表、`_tbl_id()` 锚点、配色）保留。

验证：3 个独立测试脚本从真实源码行抽取执行（`_sync_audit/test_droplog_py.py`、`test_droplog_r.R`、`test_filter_summary.py`），覆盖正常路径、无输入退化路径、失效符号链接、0 候选、表缺失行；全部通过。

## 1. 比对口径

1. 两侧各生成 manifest（`_scratch_cmp/scan_files.py`）：递归遍历仓库，收脚本/配置扩展名
   `.py .sh .bash .R .r .pl .pm .yaml .yml .toml .cfg .ini .conf .json .ipynb .nf .smk .Rmd .js .ts .qmd`，
   跳过 `.git/__pycache__/.pytest_cache/node_modules/venv/database/search_results/biosoft/analysis_results*` 及 conda 目录，
   **CRLF→LF 归一后**算 sha1，单文件上限 4 MB → `path<TAB>sha1<TAB>size`。
2. 分类（`compare_manifest.py`）：同路径同 sha1 = 同内容；同路径不同 sha1 = 差异；单侧 = 仅服务器/仅本地。
3. 差异脚本判性质（`analyze_diffs.py`）：服务器文件内容 `git hash-object` 后到本仓对象库 `git cat-file -e` 查——
   查得到 = **STALE**（服务器 = 本仓某历史版本）；查不到 = **UNIQUE**（本仓历史里从无该内容）。
4. 服务器侧 150 个差异脚本 + 1,124 个独有脚本已全部下载到 `_server_diff/live/<相对路径>`（tar 5.5 MB，1,274 个文件）；
   Δ 由 `git diff --no-index --ignore-cr-at-eol --numstat` 得出，方向 = **本地相对服务器**（`+` 为本地独有行）。

**局限**：只比内容，不比权限（服务器 88 处 `644→755`）；不比运行时行为；`database/`、`biosoft/`、结果 JSON（1,605 个）被排除；
服务器工作区另有 3,421 删除 + 153 修改（多为 `.md/.html/.tsv` 与整目录删除），不在本矩阵。

## 2. 总账

| 指标 | 服务器 | 本仓 | 说明 |
|---|---:|---:|---|
| 脚本/配置文件总数 | 1,491 | 2,058 | 扩展名过滤后 |
| 同路径同内容 | 217 | 217 | 忽略 |
| 同路径、服务器 = 历史旧版 | 107 | — | STALE：覆盖即可 |
| 同路径、两边都改过 | 43 | 43 | UNIQUE：逐个判定 |
| 仅服务器有 | 1,124 | — | 643 只是换路径 |
| 仅本地有 | — | 1,691 | 多数为归档/结果目录 |

按目录（只列脚本）：

| 目录 | 同内容 | 旧版 | 内容分歧 | 仅服务器 | 仅本地 |
|---|---:|---:|---:|---:|---:|
| `archive/20261018_server_cleanup` | 0 | 0 | 0 | 495 | 0 |
| `scripts` | 26 | 1 | 0 | 354 | 64 |
| `reports` | 21 | 10 | 5 | 0 | 310 |
| `archive/20261007_audit_cleanup` | 0 | 0 | 0 | 0 | 332 |
| `virome_phylo_pipeline` | 54 | 35 | 15 | 72 | 122 |
| `archive/20260928_phylo_results_diag` | 0 | 0 | 0 | 0 | 174 |
| `archive/20261011_endogenous_virus_pipeline` | 0 | 0 | 0 | 0 | 138 |
| `archive/20260917_root_cleanup` | 0 | 0 | 0 | 0 | 136 |
| `virome_analysis_pipeline` | 71 | 36 | 4 | 9 | 2 |
| `archive/20260917_phylo_cleanup` | 0 | 0 | 0 | 0 | 121 |
| `virome_discovery_pipeline` | 33 | 5 | 9 | 31 | 8 |
| `archive/20260928_dev_cleanup` | 0 | 0 | 0 | 7 | 41 |
| `archive/server_local_diff_20260922` | 0 | 0 | 0 | 0 | 44 |
| `(根目录)` | 0 | 2 | 0 | 36 | 4 |
| `endogenous_virus_pipeline` | 0 | 0 | 0 | 33 | 0 |
| `host_expression_pipeline` | 0 | 0 | 0 | 33 | 0 |
| `archive/20261011_host_expression_pipeline` | 0 | 0 | 0 | 0 | 33 |
| `mmpv_common` | 2 | 1 | 1 | 0 | 24 |
| `archive/20260927_oneoff_cleanup` | 0 | 0 | 0 | 0 | 24 |
| `public_metadata_pipeline` | 6 | 7 | 3 | 0 | 6 |
| `archive/debug_scripts_2026-08-21` | 0 | 0 | 0 | 20 | 0 |
| `archive/deprecated_2026-08-01` | 0 | 0 | 0 | 20 | 0 |
| `archive/scripts_20260905` | 0 | 0 | 0 | 0 | 20 |
| `archive/external_ref` | 0 | 0 | 0 | 0 | 18 |
| `文献自动下载` | 0 | 0 | 0 | 0 | 16 |
| `_refs` | 0 | 0 | 0 | 0 | 12 |
| `virome_submission_pipeline` | 4 | 2 | 5 | 0 | 0 |
| `figure_studio` / `outputs` | 0 | 0 | 0 | 0 | 各 9 |
| `data_preprocessing_pipeline` | 0 | 5 | 1 | 0 | 2 |
| `scripts_20260825` | 0 | 0 | 0 | 8 | 0 |
| `archive/20261013_root_cleanup` | 0 | 0 | 0 | 0 | 8 |
| `MMPV-RNA-paper` | 0 | 0 | 0 | 0 | 7 |
| `archive/scripts_20260825` | 0 | 0 | 0 | 5 | 0 |
| `doc` | 0 | 3 | 0 | 0 | 0 |

合计：**同内容 217 / 旧版 107 / 内容分歧 43 / 仅服务器 1,124 / 仅本地 1,691**。

## 3. 服务器停在本仓旧版的 107 个脚本（STALE）

判定依据：服务器文件内容在本仓对象库里能查到，`git log --all -- <path>` 找到的匹配提交绝大多数是
**`24b32f7`（2026-09-28，v3.1 统一 I/O 布局）**——即服务器 HEAD `b8a73f9` 就是那天从本地 24b32f7 同步出去的快照，
之后本地又改了 1–3 次。**这类文件不需要三方合并，直接用本地版覆盖即可。**

表头：`路径` — 服务器大小 → 本地大小（+增/−删行，落后提交数，匹配提交）。

### 3.1 `virome_analysis_pipeline`（36 个）

- `auto_known_virus.py` — 48,693 → 50,660 B（+40/−15，落后 3 提交，`24b32f7`）
- `batch_plot_virus_depth.py` — 46,934 → 47,065 B（+20/−0，落后 0 提交，`f46a54a`）
- `batch_virema_dvg.py` — 26,507 → 29,802 B（+61/−5，落后 1 提交，`24b32f7`）
- `batch_virus_depth.py` — 55,587 → 57,172 B（+27/−4，落后 1 提交，`24b32f7`）
- `batch_virus_full.py` — 11,194 → 13,995 B（+68/−15，落后 2 提交，`24b32f7`）
- `batch_virus_variants.py` — 60,839 → 64,608 B（+110/−37，落后 2 提交，`24b32f7`）
- `consensus_extract.py` — 10,449 → 12,557 B（+37/−7，落后 1 提交，`24b32f7`）
- `extract_full_fasta.py` — 14,097 → 15,226 B（+20/−0，落后 0 提交，`24b32f7`）
- `filter_summary.py` — 11,957 → 13,246 B（+24/−4，落后 1 提交，`24b32f7`）
- `panvirome_landscape.py` — 7,133 → 8,262 B（+20/−0，落后 0 提交，`24b32f7`）
- `panvirome_novel_known.py` — 16,059 → 16,436 B（+1/−1，落后 2 提交，`24b32f7`）
- `panvirome_oncoprint.R` — 3,373 → 4,762 B（+31/−3，落后 0 提交，`24b32f7`）
- `panvirome_prevalence.py` — 10,174 → 12,165 B（+52/−18，落后 2 提交，`24b32f7`）
- `panvirome_suite.R` — 16,327 → 18,227 B（+53/−12，落后 1 提交，`24b32f7`）
- `regenerate_coverage_plot.py` — 16,331 → 17,124 B（+20/−0，落后 0 提交，`24b32f7`）
- `snpgenie_master.py` — 41,938 → 43,067 B（+20/−0，落后 0 提交，`f46a54a`）
- `stage_funnel_summary.py` — 21,431 → 22,622 B（+26/−6，落后 1 提交，`fdb417a`）
- `summarize_assemblies.py` — 7,662 → 9,031 B（+38/−6，落后 1 提交，`24b32f7`）
- `tbl2gb.py` — 6,890 → 1,648 B（+28/−175，落后 1 提交，`24b32f7`）
- `utils/consensus_extract.py` — 10,449 → 12,557 B（+37/−7，落后 1 提交，`ec2b324`）
- `utils/extract_full_fasta.py` — 16,390 → 20,762 B（+68/−4，落后 1 提交，`24b32f7`）
- `utils/filter_summary.py` — 14,712 → 16,657 B（+43/−12，落后 2 提交，`24b32f7`）
- `utils/regenerate_coverage_plot.py` — 16,331 → 17,124 B（+20/−0，落后 0 提交，`3c0d3cb`）
- `utils/tbl2gb.py` — 6,890 → 7,499 B（+9/−1，落后 2 提交，`e1bc793`）
- `utils/variant_species_association.py` — 17,194 → 18,323 B（+20/−0，落后 0 提交，`24b32f7`）
- `utils/verify_plots_real.sh` — 1,434 → 1,759 B（+5/−2，落后 0 提交，`24b32f7`）
- `utils/virus_frequency_plot.R` — 8,646 → 10,198 B（+33/−1，落后 1 提交，`3c0d3cb`）
- `utils/virus_haplotype_qst.py` — 36,160 → 37,953 B（+37/−6，落后 1 提交，`24b32f7`）
- `utils/virus_metadata_plot.py` — 20,895 → 22,966 B（+33/−0，落后 1 提交，`24b32f7`）
- `viral_maftools.R` — 20,589 → 21,451 B（+46/−14，落后 1 提交，`24b32f7`）
- `virema_summary_report.R` — 24,844 → 26,110 B（+28/−0，落后 0 提交，`3c0d3cb`）
- `virus_auto_pipeline.py` — 63,691 → 64,820 B（+20/−0，落后 0 提交，`24b32f7`）
- `virus_detection_summary.py` — 10,303 → 10,691 B（+4/−2，落后 2 提交，`24b32f7`）
- `virus_metadata_plot.py` — 28,564 → 30,055 B（+37/−18，落后 1 提交，`24b32f7`）
- `virus_vcf_pipeline.py` — 107,455 → 108,584 B（+20/−0，落后 0 提交，`f46a54a`）
- `visual_codon_miner.py` — 10,039 → 10,933 B（+20/−0，落后 0 提交，`24b32f7`）

### 3.2 `virome_phylo_pipeline`（35 个）

- `build_haplo_outputs.py` — 12,555 → 13,732 B（+28/−3，落后 0 提交，`e851ce1`）
- `codon_usage/codon_usage_plots.py` — 9,888 → 12,148 B（+41/−6，落后 1 提交，`24b32f7`）
- `gene_align_generator.py` — 12,492 → 13,101 B（+20/−9，落后 1 提交，`24b32f7`）
- `gene_partition_dating.py` — 31,545 → 31,773 B（+3/−1，落后 1 提交，`24b32f7`）
- `geo_analysis.py` — 58,793 → 58,421 B（+18/−10，落后 2 提交，`24b32f7`）
- `phylogeo_figures.py` — 34,703 → 37,307 B（+89/−44，落后 2 提交，`24b32f7`）
- `popgen_figures.py` — 9,124 → 10,253 B（+20/−0，落后 0 提交，`24b32f7`）
- `publication_figures.py` — 29,793 → 29,919 B（+6/−4，落后 0 提交，`e851ce1`）
- `recombination_analysis.py` — 12,136 → 12,464 B（+15/−12，落后 1 提交，`24b32f7`）
- `selection_analysis.py` — 24,237 → 24,471 B（+5/−2，落后 1 提交，`24b32f7`）
- `utils/align_qc.py` — 7,347 → 7,381 B（+4/−4，落后 1 提交，`24b32f7`）
- `utils/data_collector.py` — 18,482 → 19,067 B（+17/−2，落后 0 提交，`24b32f7`）
- `utils/diffusion_stats.py` — 16,970 → 17,067 B（+4/−4，落后 1 提交，`24b32f7`）
- `utils/distance_matrix.py` — 7,758 → 7,979 B（+10/−2，落后 1 提交，`24b32f7`）
- `utils/dnasp.py` — 176,731 → 177,039 B（+11/−7，落后 0 提交，`e851ce1`）
- `utils/genetic_subsample.py` — 18,272 → 18,321 B（+1/−1，落后 1 提交，`24b32f7`）
- `utils/geo_resolver.py` — 13,923 → 16,332 B（+46/−4，落后 1 提交，`24b32f7`）
- `utils/mcmc_diag_plot.py` — 14,901 → 15,044 B（+7/−4，落后 0 提交，`e851ce1`）
- `utils/metadata_governance.py` — 63,871 → 65,858 B（+81/−38，落后 1 提交，`24b32f7`）
- `utils/path_sampling.py` — 77,866 → 80,906 B（+90/−40，落后 1 提交，`24b32f7`）
- `utils/phylogeo_bridge.py` — 50,445 → 51,207 B（+14/−0，落后 0 提交，`24b32f7`）
- `utils/phymapr_bridge.py` — 21,218 → 21,576 B（+7/−1，落后 0 提交，`e851ce1`）
- `utils/pipeline_overview_fig.py` — 36,067 → 36,537 B（+20/−0，落后 0 提交，`e851ce1`）
- `utils/popgen_r_backend.R` — 1,111 → 2,492 B（+11/−0，落后 1 提交，`24b32f7`）
- `utils/pub_plots.py` — 7,530 → 8,544 B（+22/−8，落后 1 提交，`24b32f7`）
- `utils/rdp5_bdp_plot.py` — 8,722 → 9,400 B（+11/−1，落后 0 提交，`e851ce1`）
- `utils/rdp5_validate.py` — 14,215 → 14,733 B（+12/−2，落后 0 提交，`24b32f7`）
- `utils/seq_clean.py` — 22,968 → 24,085 B（+18/−4，落后 1 提交，`24b32f7`）
- `utils/seqgrouper_bridge.py` — 6,359 → 8,186 B（+30/−3，落后 1 提交，`24b32f7`）
- `utils/splitstree_viz.py` — 31,504 → 31,676 B（+5/−3，落后 0 提交，`e851ce1`）
- `utils/tempmig_full.py` — 14,949 → 15,342 B（+9/−1，落后 0 提交，`24b32f7`）
- `utils/transmission_map.py` — 52,009 → 52,085 B（+20/−0，落后 0 提交，`e851ce1`）
- `utils/transmission_paths.py` — 34,373 → 34,436 B（+2/−2，落后 1 提交，`24b32f7`）
- `utils/treetime_bridge.py` — 10,336 → 11,154 B（+14/−0，落后 0 提交，`e851ce1`）
- `utils/world_basemap.py` — 3,321 → 4,450 B（+20/−0，落后 0 提交，`e851ce1`）

### 3.3 `reports`（10 个，均无法定位版本）

`reports/virome_analysis_pipeline/` 下 `fig_discovery_supp.py`（8,263 → 9,218 B）、`fig_ident_filter_stacked_all.R`（70,702 → 70,998 B）、
`fig_novelty_sankey.py`（42,981 → 43,408 B）、`fig_novelty_sankey_supp.py`（17,155 → 18,284 B）、
`make_cls_tool_timing_fig.R`（6,404 → 7,670 B）、`make_discovery_novelty_figs.py`（43,992 → 44,169 B）、
`make_goji_virome_phylo_figs.py`（21,704 → 22,318 B）、`make_sci_composite_figs.py`（203,876 → 201,623 B）、
`make_tool_timing_fig.R`（9,867 → 11,133 B）、`redraw_barbarum_detection_figs.py`（18,696 → 19,825 B）。
本地这些文件从未入库（path 历史为空），服务器内容是本地某个未提交版本的产物；只影响出图，不影响管线逻辑。

### 3.4 `public_metadata_pipeline`（7 个）

- `build_host_pipeline.py` — 16,029 → 37,740 B（+443/−17，落后 2 提交，`24b32f7`；本地含本轮 eve-mask 阶段）
- `download_host_genome.py` — 24,521 → 24,664 B（+6/−4，落后 1 提交，`24b32f7`）
- `gsa_sra.down.py` — 22,737 → 22,884 B（+4/−0，落后 1 提交，`24b32f7`）
- `gsa_sra.info.py` — 74,234 → 77,465 B（+68/−9，落后 2 提交，`24b32f7`）
- `preprocess_unified.py` — 19,293 → 22,643 B（+75/−12，落后 3 提交，`24b32f7`）
- `public_data_pipeline.py` — 20,824 → 35,537 B（+282/−29，落后 3 提交，`24b32f7`）
- `tests/test_pipeline_utils.py` — 10,399 → 17,230 B（+177/−0，落后 0 提交，`24b32f7`）

### 3.5 `data_preprocessing_pipeline`（5 个）

- `clean-data.py` — 22,471 → 25,746 B（+73/−13，落后 2 提交，`24b32f7`）
- `data_preprocessing.py` — 19,339 → 19,816 B（+14/−8，落后 3 提交，`24b32f7`）
- `host_depletion.py` — 63,007 → 65,236 B（+68/−20，落后 2 提交，`24b32f7`；`--keep-taxids` 语义差异见 doc/17 §3）
- `preprocess.py` — 7,339 → 9,078 B（+33/−5，落后 0 提交，`24b32f7`）
- `run_bbnorm.py` — 7,025 → 7,991 B（+19/−4，落后 1 提交，`24b32f7`）

### 3.6 `virome_discovery_pipeline`（5 个）

- `graph_propagation.py` — 42,251 → 45,174 B（+57/−7，落后 1 提交，`8fea067`）
- `tests/test_cobra_source_fix.py` — 7,469 → 10,752 B（+66/−2，落后 1 提交，`8fea067`）
- `tests/test_graph_propagation.py` — 23,073 → 26,180 B（+54/−4，落后 1 提交，`8fea067`）
- `utils/acvirus_tree_pro.py` — 56,098 → 57,227 B（+20/−0，落后 0 提交，`24b32f7`）
- `utils/discovery2analysis.py` — 12,387 → 15,062 B（+56/−10，落后 1 提交，`8fea067`）

### 3.7 其余（`doc` 3 / 根目录 2 / `virome_submission_pipeline` 2 / `mmpv_common` 1 / `scripts` 1）

- `doc/figures/make_figure1_capabilities.py`（+20/−0）、`make_figure1_full.py`（+20/−0）、`make_module_overview.py`（+24/−13）
- 根：`pipeline_config.yaml` 11,087 → 11,964 B（+16/−6，落后 1 提交，`ad275e3`）、`pixi.toml`（+16/−5，落后 3 提交，`167cd6b`）
- `virome_submission_pipeline/report_html.py`（+268/−230）、`submission_gui.py`（+22/−36）
- `mmpv_common/__init__.py` — 154 → 1,354 B（+19/−2）
- `scripts/audit/check_io_handoffs.py`（+7/−3，落后 3 提交，`89c7a77`）

## 4. 两边都改过的 43 个脚本（UNIQUE）

这些文件在服务器上的内容**本仓历史里从未出现过**：两边在 `24b32f7` 之后各自改了同一文件。

### 4.1 谁更新？——三种证据互相矛盾，只能靠人工抽查

| 证据 | 结果 |
|---|---|
| 本地最后一次提交时间 | 多数停在 2026-10-13（本地 HEAD），服务器 mtime 是 10-09～10-19 → 看似服务器更新 |
| 文件内容里写的日期（`2026-MM-DD`） | 26 个本地更新（10-16～10-20，多数是未提交批次）、5 个服务器更新（10-18/19 `benchfix`）、12 个无/同 |
| 人工抽查 4 组 | `io_layout.py`：本地只改了注释/归档说明（**低风险**）；`pipeline_config.yaml`：本地只加文档性说明（**低风险**）；`filter_virus.py`/`virus_identification.py`：本地**已有等价修复**（最长前缀归属、类病毒并入）；`run_host_prediction.py`/`virus_classifier.py`：本地**缺** `dropped_ids.tsv` 审计（**真缺口**） |

**结论**：日期信号不可信（本地 10-20 的 `figure_export` 未提交批次会抬高所有文件的「内容日期」）。处置按 §4.2 分组；
已确认可直接忽略的：`io_layout.py`、`pipeline_config.yaml`；已确认必须回移的：发现管线「丢弃审计」2 个文件；其余按表格逐个人工看。

### 4.2 `virome_phylo_pipeline`（15 个）

| 文件 | Δ（本地 vs 服务器） | 判定 |
|---|---:|---|
| `phylo_pipeline.py` | +11/−4 | 小改，以本地为准 |
| `saturation_analysis.py` | +17/−8 | 本地 +17/−8 |
| `utils/beast_bridge.py` | +15/−0 | 本地新增兼容层 |
| `utils/beast_parser.py` | +3/−1 | 本地小改 |
| `utils/clock_analysis.py` | +20/−0 | 本地新增 |
| `utils/core_results_fig.py` | +20/−0 | 本地新增 |
| `utils/ess_convergence_fig.py` | +20/−0 | 本地新增（本地未入库） |
| `utils/popgen_plots.py` | +13/−6 | 本地小改 |
| `utils/report_builder.py` | +256/−36 | 本地大扩写（报告 v2），以本地为准 |
| `utils/sample_plot.py` | +9/−2 | 本地小改 |
| `utils/splitstree_native_viz.py` | +6/−3 | 本地小改 |
| `utils/stage_plots.py` | +41/−12 | 本地扩写 |
| `utils/temporal_signal.py` | +3/−1 | 本地小改 |
| `utils/virphy_bridge.py` | +52/−18 | 本地路径/中文修复，以本地为准 |
| `visual_codon_miner.py` | +13/−3 | 本地小改 |

### 4.3 `virome_discovery_pipeline`（9 个）

| 文件 | Δ | 判定 |
|---|---:|---|
| ★ `run_host_prediction.py` | +0/−43 | **本地缺 2026-10-19 丢弃审计**（`dropped_ids.tsv`、`duplicate_tax_ids.tsv`）→ 建议回移 |
| ★ `virus_classifier.py` | +9/−48 | **本地缺同一个丢弃审计**（`_write_classifier_drop_log`）→ 建议回移；本地另有 `encoding='utf-8'` 修复需保留 |
| `filter_virus.py` | +57/−49 | 本地已有等价「最长匹配前缀」修复（:190），差异在其它细节 → 人工看 |
| `virome_pipeline.py` | +62/−51 | 服务器有 2026-10-18 benchfix（分阶段 `reads_dir`、样本名前缀）；本地是 15 阶段重构 → **人工合** |
| `report_pipeline.py` | +212/−57 | 服务器有 10-19 符号链接判定修复；本地报告重构 → 人工确认该分支是否覆盖 |
| `virus_identification.py` | +97/−13 | 本地已有等价「类病毒并入 all.candidate」（:1342）→ 以本地为准 |
| `cobra_pipeline.py` | +17/−3 | 服务器 10-18 `_clean` 命名变体；本地是否有等价需确认 |
| `virus_classifier_analysis.R` | +32/−47 | R 绘图重构，抽查输出列名 |
| `utils/sci_palette.py` | +6/−0 | 本地新增色板，低风险 |

### 4.4 `virome_submission_pipeline`（5 个）

| 文件 | Δ | 判定 |
|---|---:|---|
| `submission_pipeline.py` | +130/−25 | 本地重构；提交链风险高，覆盖前隔离空跑 |
| `virome_submission.py` | +41/−41 | 近似重写，需三方比对 |
| `unified_metadata.py` | +24/−22 | 字段映射改写，抽查输出列 |
| `known_virus_submit.py` | +12/−4 | 小改 |
| `sequin_builder.py` | +7/−7 | 参数改名，抽查 `.sbt` 模板字段 |

### 4.5 `reports`（5 个，均为科学图脚本）

`virome_analysis_pipeline/` 下 `fig_host_prediction_b_panel.R`（+28/−0）、`fig_upset_classification.R`（+34/−9）、
`make_discovery_sci_figs.py`（+158/−16）、`make_goji_virome_phylosynteny_fig.py`（+20/−0）、`upset_common.R`（+28/−0）。
只读产物，以本地为准；抽查输出文件名是否被下游引用。

### 4.6 `virome_analysis_pipeline`（4 个）

| 文件 | Δ | 判定 |
|---|---:|---|
| `generate_pipeline_report.py` | +182/−59 | 本地报告重构；`_find_virus_files` 被删，确认新实现覆盖 |
| `snpeff_analysis.py` | +39/−10 | 本地改了 R 出图链，确认仍完整 |
| `virus-full.py` | +82/−54 | 双向改写，抽查命令行参数兼容性 |
| `virus_variants_analyzer.py` | +46/−1 | 低风险 |

### 4.7 `public_metadata_pipeline`（3 个）

| 文件 | Δ | 判定 |
|---|---:|---|
| `generate_report.py` | +91/−27 | 本地报告修复；抽查 HTML 段落是否缺块 |
| `gsa_sra.search.py` | +282/−9 | 本地检索重试/字段归一；生产在跑服务器版，覆盖前小样本试跑 |
| `gsa_sra.plot.py` | +20/−0 | 本地纯新增绘图修复，低风险 |

### 4.8 `data_preprocessing_pipeline`（1 个）/ `mmpv_common`（1 个）

| 文件 | Δ | 判定 |
|---|---:|---|
| `preprocess_report.py` | +514/−80 | 本地 HTML 报告 v2；服务器 27 KB 版可能各有字段名修复，抽查 `build_rows` 输出列 |
| `io_layout.py` | +18/−6 | **已抽查：本地只改注释/归档说明（10-11 归档标注 + `resolve_layout` docstring）**，无目录名变动 → 低风险，以本地为准 |

## 5. 仅服务器有的 1,124 个脚本

### 5.1 只是换了路径（643 个，内容与本地某处完全一致）

| 服务器目录 | 本地对应目录（内容相同） | 文件数 |
|---|---|---:|
| `scripts` | `archive/20261007_audit_cleanup` | 314 |
| `archive/20261018_server_cleanup` | `archive/20260928_phylo_results_diag` | 171 |
| `virome_phylo_pipeline` | `archive/20260917_phylo_cleanup` | 41 |
| `(根目录)` | `reports` | 27 |
| `virome_discovery_pipeline` | `archive/20260928_dev_cleanup` | 18 |
| `scripts` | `reports` | 12 |
| `host_expression_pipeline` | `archive/20261011_host_expression_pipeline` | 10 |
| `archive/20261018_server_cleanup` | `public_metadata_pipeline` | 7 |
| `endogenous_virus_pipeline` | `archive/20261011_endogenous_virus_pipeline` | 7 |
| `archive/20261018_server_cleanup` | `archive/server_local_diff_20260922` | 6 |

这些是「服务器上还没清理的旧副本」，本地已归档到 `archive/` 或 `reports/_archive/`，**不需要回迁**。

### 5.2 真正独有的 481 个（内容本仓从未有过）

| 目录 | 脚本数 | 占比 | 性质 |
|---|---:|---:|---|
| `archive/20261018_server_cleanup` | 307 | 64% | 服务器自己的清理包，可忽略 |
| `endogenous_virus_pipeline` | 26 | 5% | EVE v6.3，本地 archive 有 v7.x |
| `scripts` | 26 | 5% | 审计脚本 |
| `virome_phylo_pipeline` | 26 | 5% | 含 `_phymap_workflow` 6、`_review20260916` 4 |
| `host_expression_pipeline` | 23 | 5% | 本地 archive 有更新版 |
| `archive/debug_scripts_2026-08-21` | 20 | 4% | 一次性调试 |
| `archive/deprecated_2026-08-01` | 20 | 4% | 一次性调试 |
| `(根目录)` | 8 | 2% | Fig3/SplitsTree 色板调试 |
| `scripts_20260825` | 8 | 2% | 一次性 |
| `virome_discovery_pipeline` | 8 | 2% | 试验脚本 |

剔除 `archive/` 与 `scripts/` 归档后，真正「服务器在用、本地没有」的活代码 **165 个**（`_scratch_cmp/one_sided.tsv`）：

| 服务器目录 | 脚本数 | 示例 |
|---|---:|---|
| `(根目录)` | 36 | `_ab_r4.py`、`_check_fig3_overlap_mod.py`、`_probe_native_nb.py` |
| `endogenous_virus_pipeline` | 33 | `eve_screen.py`、`eve_scan_core.py`、`eve_distinguish/*`（v6.3） |
| `host_expression_pipeline` | 33 | `host_expression.py`、`salmon_auto_pipeline.py`、`teagcn/**` |
| `virome_discovery_pipeline` | 31 | `build_famtax.py`、`batch_sdt.sh`、`linclust.sh` |
| `scripts` | 16 | `1_cross_kingdom_overlap.py`、`3_host_probability.py` |
| `virome_phylo_pipeline` | 12 | `_backup_phylo_pipeline_v1_20260827.py`、`verify_docs_vs_code.py` |
| `virome_analysis_pipeline` | 4 | `_debug_stage1.py`、`_report_charts.py` |

## 6. 仅本地有的 1,691 个脚本

| 目录 | 脚本数 | 说明 |
|---|---:|---|
| `archive/20261007_audit_cleanup` | 332 | 本地归档（服务器已清理） |
| `reports` | 310 | 出图/归档 |
| `archive/20260928_phylo_results_diag` | 174 | 归档 |
| `archive/20261011_endogenous_virus_pipeline` | 138 | **EVE v7.x（更新版本）** |
| `archive/20260917_root_cleanup` | 136 | 归档 |
| `virome_phylo_pipeline` | 122 | 本地新增模块/工作流 |
| `archive/20260917_phylo_cleanup` | 121 | 归档 |
| `scripts` | 64 | 本地审计脚本 |
| `archive/server_local_diff_20260922` | 44 | 归档 |
| `archive/20260928_dev_cleanup` | 41 | 归档 |
| `archive/20261011_host_expression_pipeline` | 33 | **宿主表达更新版** |
| `mmpv_common` | 24 | 本地新增模块 |
| 其余（`archive/*`、`_refs`、`文献自动下载`、`figure_studio` 等） | … | 归档/工具 |
| `public_metadata_pipeline` | 6 | **含本轮 `mask_host_genome.py` 与 2 个测试** |

真正的**新增代码**只有三处需要关心：`archive/20261011_endogenous_virus_pipeline`（EVE v7.x，138 个）、
`archive/20261011_host_expression_pipeline`（33 个）、以及本轮 eve-mask 三件套；其余是归档与结果 JSON（1,605 个）。

## 7. 两条「搬家」管线的版本对照

服务器放在顶层目录，本地已挪进 `archive/20261011_*` 并继续改过。**同名文件内容不同 = 两边版本不同，不是分支冲突。**

### 7.1 内源病毒：服务器 v6.3 vs 本地 archive v7.x

| 指标 | 数量 |
|---|---:|
| 同名同内容 | 7 |
| 同名不同内容 | 26 |
| 仅服务器 | 0 |
| 仅本地 | 236 |

差异最大的 10 个：`tests/test_eve_distinguish_logic.py` +309/−28、`eve_distinguish/s3_verdict.py` +236/−25、
`eve_scan_core.py` +187/−20、`tests/test_eve_core.py` +165/−1、`eve_element_db.py` +115/−29、
`tests/test_eve_distinguish.py` +113/−5、`eve_distinguish/run_all.sh` +97/−15、`eve_distinguish/s2b_locus_scan.py` +98/−13、
`eve_screen.py` +86/−7、`eve_report.py` +65/−17。

两版 `s2_verdict.tsv` 表头**一致**（`locus	verdict	best_viral_id	best_viral_bs	...`），kingdom 汇总表也一致
→ `mask_host_genome.parse_verdict_tsv` 对两版都能读。但**判定数值会不同**（v7.x 修了读框偏移、证据闸门、侧翼上下文），
同一基因组在两版上可能得到不同的 `viral_supported/host_like` 计数 → 文档里必须写明用哪版。

### 7.2 宿主表达：服务器 vs 本地 archive

| 指标 | 数量 |
|---|---:|
| 同名同内容 | 10 |
| 同名不同内容 | 23 |
| 仅服务器 | 0 |
| 仅本地 | 40 |

差异最大：`salmon_auto_pipeline.py` +718/−200、`tests/test_host_expression.py` +454/−15、
`host_expression_report.py` +56/−25、`diff_coexp.py` +52/−13、`teagcn/coexpression/TEA_mp.py` +33/−9、
`coexp_module_conservation.py` +34/−4。

## 8. 与本轮 eve-mask 改动的交集（原问题的直接回答）

| 检查项 | 结果 | 结论 |
|---|---|---|
| 同路径文件是否被双方修改 | 只有 `public_metadata_pipeline/build_host_pipeline.py`；服务器内容 = 本仓 `24b32f7` 旧版（16,029 B CRLF） | **无三方合并**，直接用本地版覆盖 |
| 掩码脚本/测试是否存在冲突版本 | 服务器没有 `mask_host_genome.py`、`test_mask_host_genome.py`、`test_eve_mask_stage.py` | **纯新增** |
| 服务器是否已有其它掩码实现 | `grep -l mask` 全仓只命中 API-key 脱敏与 `virome_phylo_pipeline/mask_recombination.py`（重组检测） | **无功能重叠** |
| `--eve-scan` 调用的 EVE 脚本 | `_run_eve_scan` 先找 `archive/20261011_.../eve_screen.py`，再找服务器顶层 `endogenous_virus_pipeline/eve_screen.py` | 两边都能跑；表头兼容；建议 `--eve-scan-script` 固定版本 |
| `data_preprocessing_pipeline` 是否要改 | 代码不改；只需 `--host_db` 指向新库目录 | 与服务器旧版无冲突 |

**部署建议（B′ 方案，最小面）**——只搬本次改动，不整树同步：

1. 服务器先冻结：`git checkout -b server-wip-20261022 && git add -A && git commit -m "wip: freeze server worktree before eve-mask sync"`；
2. 传 6 个文件：`public_metadata_pipeline/{build_host_pipeline.py,mask_host_genome.py}`、
   `public_metadata_pipeline/tests/{test_mask_host_genome.py,test_eve_mask_stage.py}`、
   `doc/{16-host-eve-masking.md,17-server-vs-local-audit-20261022.md}`（+ 本文）；
3. 服务器上 `python -m pytest public_metadata_pipeline/tests/test_eve_mask_stage.py public_metadata_pipeline/tests/test_mask_host_genome.py -q` 应 35 passed；
4. 确认服务器 `endogenous_virus_pipeline/eve_screen.py`（v6.3）能被 `--eve-scan` 自动发现；否则加 `--eve-scan-script`。

## 9. 待确认与风险

1. **丢弃审计缺口（已确认）**：服务器 2026-10-19 在 `run_host_prediction.py` / `virus_classifier.py` 加了 `dropped_ids.tsv`
   记录被静默丢弃的序列（注释写明「610425 → 535485 的 74,943 条缺口无据可查」），本仓没有 → 建议优先回移。
2. **43 个 UNIQUE 中 41 个仍未逐行判定**：日期证据（mtime / commit date / 内容日期）互相矛盾，本报告只人工抽查了 4 组；
   其中 `virome_pipeline.py`、`report_pipeline.py`、`cobra_pipeline.py`、`submission_pipeline.py`、`virome_submission.py`
   属高风险，必须在覆盖前逐行看。
3. **107 个 STALE 里 12 个属生产链路**（`data_preprocessing_pipeline` 5、`public_metadata_pipeline` 7）——服务器现在跑的就是这些旧版；
   若要复跑历史样本，结果与本地不可比。
4. **服务器工作区 3,770 条未提交变更**（3,421 删除 + 153 修改 + 196 未跟踪）未纳入矩阵；整树 `git pull` 会把删除变成
   「已删除的已跟踪文件」，必须先冻结提交。
5. **88 处 `644→755` 权限差异**未比较，脚本可执行位可能影响 `subprocess` 调用。
6. **生产样本输出位置未确认**：`/data`、`/alldata`、`/mnt` 均为空，`/dev/sda2` 500G 只有系统与 home；
   1,348 样本大概率在别的节点/挂载点，需用户提供路径后才能核实 `host_depletion` 版本（doc/17 §3）。
7. **两版 EVE 管线的判定数值不可混用**（§7.1）：`--eve-scan` 在本地跑与在服务器跑，同一基因组的 `masked_frac` 可能不同。

## 10. 附录：复核命令

```bash
# 1) 重新生成两侧清单（本地）
python _scratch_cmp/scan_files.py --root . --out _scratch_cmp/manifest_local.tsv
# 2) 服务器侧（在服务器上）
python scan_files.py --root ~/MMPV-RNA --out /tmp/manifest_server.tsv
scp success:/tmp/manifest_server.tsv _scratch_cmp/
# 3) 分类 + 差异明细
python _scratch_cmp/compare_manifest.py
python _scratch_cmp/analyze_diffs.py
# 4) 查看某个文件的三方内容
git diff --no-index --ignore-cr-at-eol _server_diff/live/<path> <path>
# 5) 确认服务器某脚本属于哪个历史版本
git hash-object _server_diff/live/<path>   # 再到本仓 git cat-file -e <blob>
```

---

生成时间：2026-10-22 ｜ 数据：`_scratch_cmp/{manifest_server,manifest_local,diff_summary,unique_summary,uniq_dates,relocated,one_sided}.tsv|json`，
服务器副本 `_server_diff/live/`（1,274 个文件）。
