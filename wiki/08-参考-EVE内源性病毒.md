# 08-参考 · EVE内源性病毒 — 脚本 `--help` 参数库 (自动提取)

> 提取自本仓库 6 个核心脚本的实际 `--help`/docstring (2026-09-28)。
 提取超时或依赖缺失者回退 docstring 并标注 `[docstring]`。以代码为准。

## `endogenous_virus_pipeline/eve_element_db.py`

```text
usage: eve_element_db.py [-h] [-o OUTDIR] [-B BATCH] [--stages STAGES]
                         [--force] [-t THREADS] [-J JOBS]
                         [--cmd-timeout CMD_TIMEOUT] [--verdicts VERDICTS]
                         [--min-bs MIN_BS] [--id ID] [--cov COV]
                         [--er-dist ER_DIST] [--mafft-consensus] [--rm]
                         [--rm-bin RM_BIN] [--cdd-db CDD_DB]
                         [--cdd-sens CDD_SENS] [--cdd-max-seqs CDD_MAX_SEQS]
                         [--rt-ref RT_REF] [--run-phylo]
                         [--min-clusters MIN_CLUSTERS] [--mmseqs MMSEQS]
                         [--mafft MAFFT]

元件级 (element-level) EVE 后处理层: 跨基因组聚类 → 元件区 → consensus → 拷贝数/结构域/系统发育

options:
  -h, --help            show this help message and exit

输入/输出:
  -o OUTDIR, --outdir OUTDIR
                        eve_screen.py 的输出根目录 (也是本层输出根目录) (default:
                        eve_results)
  -B BATCH, --batch BATCH
                        批量 TSV: 每行 NAME<TAB>基因组路径 (rm 阶段必需; 与 eve_screen.py -B
                        同一张表) (default: None)
  --stages STAGES       阶段: collect/cluster/regions/consensus/rm/cdd/phylo
                        逗号组合 (别名 1..7) 或 all; all 会连带 rm/cdd/phylo, 需要相应参数
                        (default: collect,cluster,regions,consensus)
  --force               忽略断点, 重跑请求的阶段 (rm 阶段连带忽略每基因组的 .out 标记) (default:
                        False)

运行:
  -t THREADS, --threads THREADS
  -J JOBS, --jobs JOBS  rm 阶段并行基因组数 (总核数 = jobs x threads) (default: 4)
  --cmd-timeout CMD_TIMEOUT
                        单条外部命令超时 (秒); 0=不限 (防工具挂死占核) (default: 0)

筛选/聚类口径 (与源脚本一致):
  --verdicts VERDICTS   collect/regions 收哪些 verdict (逗号分隔) (default:
                        viral_supported)
  --min-bs MIN_BS       evidence_bs 下限 (max(ref_bitscore, rvdb_bitscore))
                        (default: 50.0)
  --id ID               easy-linclust --min-seq-id (default: 0.8)
  --cov COV             easy-linclust -c (default: 0.8)
  --er-dist ER_DIST     元件区归并距离 (nt, ±同家族位点) (default: 5000)
  --mafft-consensus     consensus 走 mafft 比对共识 (默认取簇内最长成员) (default: False)

需要额外依赖的阶段:
  --rm                  把 rm 阶段加进本次运行 (需 -B 与 RepeatMasker) (default: False)
  --rm-bin RM_BIN       RepeatMasker 可执行文件 (默认取 RM_BIN 环境变量, 再默认 PATH)
                        (default: RepeatMasker)
  --cdd-db CDD_DB       mmseqs profile 格式的 CDD 库前缀 (给定即把 cdd 阶段加进本次运行)
                        (default: )
  --cdd-sens CDD_SENS   mmseqs search -s (default: 4.0)
  --cdd-max-seqs CDD_MAX_SEQS
                        mmseqs search --max-seqs (default: 300)
  --rt-ref RT_REF       RT 参考集 fasta (--run-phylo 必需) (default: )
  --run-phylo           本机执行 phylo_commands.sh (需 --rt-ref 与
                        mafft/trimal/iqtree) (default: False)
  --min-clusters MIN_CLUSTERS
                        phylo: 家族至少要有这么多个簇才出代表集 (default: 10)

工具:
  --mmseqs MMSEQS       mmseqs 可执行文件 (默认 PATH) (default: None)
  --mafft MAFFT         mafft 可执行文件 (默认 PATH) (default: None)
```

## `endogenous_virus_pipeline/eve_genome_scan.py`

```text
usage: eve_genome_scan.py [-h] -g GENOME [-n NAME] -o OUTDIR [-t THREADS]
                          [--stage STAGE] [--fast] [--overlap]
                          [--blastn-viroid] [--skip-s3] [--cleanup] [--force]
                          [--log LOG] [--ref-db REF_DB] [--pv-db PV_DB]
                          [--id2div ID2DIV] [--rvdb-db RVDB_DB]
                          [--viroids-db VIROIDS_DB] [--diamond DIAMOND]
                          [--samtools SAMTOOLS] [--window WINDOW]
                          [--merge-distance MERGE_DISTANCE] [--evalue EVALUE]
                          [--host-bs HOST_BS] [--viral-bs VIRAL_BS]
                          [--cmd-timeout CMD_TIMEOUT]

单基因组 EVE 三阶段筛查 worker

options:
  -h, --help            show this help message and exit
  -g GENOME, --genome GENOME
                        宿主基因组 FASTA (支持 .fa/.fasta/.fna/.gz/.tar.gz) (default:
                        None)
  -n NAME, --name NAME  基因组名 (默认: 文件名去后缀) (default: None)
  -o OUTDIR, --outdir OUTDIR
                        输出根目录 (default: None)
  -t THREADS, --threads THREADS
  --stage STAGE, --stages STAGE
                        阶段: 1/discover,2/verdict,3/annotate 逗号组合 或 all
                        (default: 1,2,3)
  --fast                Stage1 用 --fast (约快6倍, 位点数约-54%) (default: False)
  --overlap             切块步长减半 (25kb 重叠, 修跨界截断; 重点基因组复扫用) (default: False)
  --blastn-viroid       附加类病毒 blastn 层 (默认关) (default: False)
  --skip-s3             跳过 RVDB 层 (default: False)
  --cleanup             完成后删除已被下游产物取代的中间文件 (default: False)
  --force               忽略断点全量重跑 (default: False)
  --log LOG             日志文件路径 (默认仅 stderr) (default: None)

数据库 (单跑时用 CLI 指定; 批量跑由编排器按 pipeline_config.yaml 注入):
  --ref-db REF_DB       Stage1 病毒参考蛋白 diamond 库 (default: None)
  --pv-db PV_DB         Stage2 植物/病毒拆分库 (default: None)
  --id2div ID2DIV       Stage2 sseqid->viral/plant 表 (default: None)
  --rvdb-db RVDB_DB     Stage3 RVDB diamond 库 (default: None)
  --viroids-db VIROIDS_DB
                        类病毒 fasta (blastn 层) (default: None)

工具:
  --diamond DIAMOND     diamond 可执行文件 (默认 PATH) (default: None)
  --samtools SAMTOOLS   samtools 可执行文件 (默认 PATH) (default: None)

运行参数:
  --window WINDOW
  --merge-distance MERGE_DISTANCE
  --evalue EVALUE
  --host-bs HOST_BS
  --viral-bs VIRAL_BS
  --cmd-timeout CMD_TIMEOUT
                        单条外部命令超时 (秒); 0=不限 (防工具挂死占核) (default: 0)
```

## `endogenous_virus_pipeline/eve_qc_contam.py`

```text
usage: eve_qc_contam.py [-h] --outdir OUTDIR [--layout {auto,ours,kingdom}]
                        [--frac FRAC] [--json JSON] [--md MD]

对一次已跑完的 EVE 筛查结果事后做测序污染筛查

options:
  -h, --help            show this help message and exit
  --outdir OUTDIR       筛查输出根目录（含 01_Loci/）
  --layout {auto,ours,kingdom}
                        ours=本管线布局(01_Loci/<NAME>/<NAME>.s1_raw.tsv, 占比按该
                        contig 全部命中算); kingdom=上游布局(<NAME>/<NAME>.s1_best.tsv,
                        占比按该基因组位点数算); auto=自动识别
  --frac FRAC           污染命中占比达到多少算'建议剔除'（默认 0.5；只影响标记，明细表恒含全部有污染命中的条目）
  --json JSON           汇总 JSON 输出路径
  --md MD               质控记录输出路径（默认 <outdir>/QC_CONTAMINATION.md）
```

## `endogenous_virus_pipeline/eve_report.py`

```text
usage: eve_report.py [-h] -o OUTDIR [--out OUT] [--top TOP]

EVE 筛查汇总 → 自包含 HTML 报告 (在 eve_screen --merge 之后运行)

options:
  -h, --help            show this help message and exit
  -o OUTDIR, --outdir OUTDIR
                        EVE 筛查输出根目录 (含 kingdom_summary.tsv 等 --merge 产物)
  --out OUT             输出 HTML 路径 (默认: <outdir>/04_Summary/eve_report.html)
  --top TOP             家族/位点明细表展示条数 (默认: 20)
```

## `endogenous_virus_pipeline/eve_scan_core.py`

```text
[docstring] eve_scan_core.py — 内源性病毒元件 (EVE) 筛查核心逻辑
====================================================
对宿主参考基因组做三阶段 EVE 筛查 (内源性病毒化石位点):

  Stage 1  discover  正向发现: 50kb 滑窗 → diamond blastx vs 病毒参考蛋白库
                      → 坐标还原 → 位点合并 (300nt) → 抽序列 → 最佳参考命中 + 科级映射
  Stage 2  verdict   双侧判定: loci vs 植物/病毒拆分库 (id2div 映射)
                      → viral_supported / host_like / undetermined
                      (宿主蛋白误匹配排除: 两侧 bitscore 对比, 防植物基因假阳性)
  Stage 3  annotate  RVDB 深度归属: viral_supported + undetermined 位点 vs RVDB

输出目录结构 (以 --output_dir 为根):
  01_Loci/<NAME>/         Stage1 产物 (fna/chunks/s1_raw/s1_hits/loci.bed/loci.fa/s1_best)
  02_Verdict/<NAME>/      Stage2 产物 (s2_raw/s2_verdict)
  03_RVDB/<NAME>/         Stage3 产物 (cand3.bed/cand3.fa/s3_rvdb)
  04_Summary/             <NAME>_eve_summary.tsv (+ 可选 <NAME>_viroid.tsv)
  logs/                   每基因组日志 (编排器写入)
  --merge 汇总产物写在根目录: kingdom_summary.tsv / kingdom_loci_all.tsv.gz /
                            family_by_genome.tsv

断点续传: 每个阶段产物文件存在即视为完成 (空文件 = 已完成且零结果),
--force 删除旧产物重跑.

来源: 移植自 hi-fever/eve_kingdom.py (2026-09-19 定版, 246 服务器千种植物验证),
数据库路径/工具路径全部参数化 (EveConfig), 适配 MMPV-RNA 配置体系.
数据文件读写统一 utf-8 (服务器 C locale / 本机 GBK 默认编码下不炸)
```

## `endogenous_virus_pipeline/eve_screen.py`

```text
usage: eve_screen.py [-h] [-g GENOME] [-n NAME] [-B BATCH] [--list]
                     [-o OUTDIR] [-t THREADS] [-J JOBS] [--stage STAGE]
                     [--merge] [--no-merge] [--force] [--fast] [--overlap]
                     [--skip-s3] [--blastn-viroid] [--cleanup]
                     [--config CONFIG] [--profile PROFILE] [--ref-db REF_DB]
                     [--pv-db PV_DB] [--id2div ID2DIV] [--rvdb-db RVDB_DB]
                     [--viroids-db VIROIDS_DB] [--diamond DIAMOND]
                     [--samtools SAMTOOLS] [--window WINDOW]
                     [--merge-distance MERGE_DISTANCE] [--evalue EVALUE]
                     [--host-bs HOST_BS] [--viral-bs VIRAL_BS]
                     [--cmd-timeout CMD_TIMEOUT]

内源性病毒元件 (EVE) 三阶段筛查编排器

options:
  -h, --help            show this help message and exit

输入:
  -g GENOME, --genome GENOME
                        单基因组 FASTA (.fa/.fasta/.fna/.gz/.tar.gz) (default:
                        None)
  -n NAME, --name NAME  基因组名 (默认: 文件名去后缀) (default: None)
  -B BATCH, --batch BATCH
                        批量 TSV: 每行 NAME<TAB>基因组路径 (default: None)
  --list                仅列出任务 (name\tpath) 不执行 (default: False)

运行:
  -o OUTDIR, --outdir OUTDIR
                        输出根目录 (default: eve_results)
  -t THREADS, --threads THREADS
  -J JOBS, --jobs JOBS  并行基因组数 (总核数 = jobs x threads) (default: 1)
  --stage STAGE         阶段: discover/verdict/annotate 逗号组合 或 all (default:
                        all)
  --merge               只跑跨基因组汇总 (不筛查) (default: False)
  --no-merge            批量完成后不自动汇总 (default: False)
  --force               忽略断点全量重跑 (已完成基因组也重跑) (default: False)
  --fast                Stage1 用 --fast (约快6倍, 位点数约-54%) (default: False)
  --overlap             切块 25kb 重叠 (修跨界截断; 重点基因组复扫用) (default: False)
  --skip-s3             跳过 RVDB 层 (default: False)
  --blastn-viroid       附加类病毒 blastn 层 (默认关) (default: False)
  --cleanup             完成后删除已无用中间文件 (chunks/raw/hits) (default: False)
  --config CONFIG       YAML 配置文件 (默认: 自动查找 pipeline_config.yaml) (default:
                        None)
  --profile PROFILE     配置预设 (default/plant, 可选见 pipeline_config.yaml)
                        (default: default)

数据库 (默认: 从 pipeline_config.yaml 读取):
  --ref-db REF_DB       Stage1 病毒参考蛋白 diamond 库 (default: None)
  --pv-db PV_DB         Stage2 植物/病毒拆分库 (default: None)
  --id2div ID2DIV       Stage2 sseqid->viral/plant 表 (default: None)
  --rvdb-db RVDB_DB     Stage3 RVDB diamond 库 (default: None)
  --viroids-db VIROIDS_DB
                        类病毒 fasta (blastn 层) (default: None)

工具:
  --diamond DIAMOND     diamond 可执行文件 (默认 PATH) (default: None)
  --samtools SAMTOOLS   samtools 可执行文件 (默认 PATH) (default: None)

EVE 参数:
  --window WINDOW       滑窗宽度 (nt) (default: 50000)
  --merge-distance MERGE_DISTANCE
                        位点合并距离 (nt) (default: 300)
  --evalue EVALUE
  --host-bs HOST_BS     Stage2 植物侧判定阈值 (bitscore) (default: 50)
  --viral-bs VIRAL_BS   Stage2 病毒侧判定阈值 (bitscore) (default: 50)
  --cmd-timeout CMD_TIMEOUT
                        单条外部命令超时 (秒); 0=不限 (防工具挂死占核) (default: 0)
```
