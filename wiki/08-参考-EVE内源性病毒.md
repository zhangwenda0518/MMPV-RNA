# 08-参考 · EVE内源性病毒 — 脚本 `--help` 参数库 (自动提取)

> **现行结构**：EVE 拆成两段 —— `public_metadata_pipeline/eve_genome/`（基因组侧）与
> `virome_discovery_pipeline/eve_discriminate/`（转录组候选判别）。本页由各脚本实际
> `--help` 提取；依赖缺失或超时的条目回退 docstring 并标注 `[docstring]`。以代码为准。

> 覆盖 20 个脚本：eve_genome 6 + eve_discriminate 14。

## `public_metadata_pipeline/eve_genome/eve_screen.py`

```text
usage: eve_screen.py [-h] [-g GENOME] [-n NAME] [-B BATCH] [--list]
                     [-o OUTDIR] [-t THREADS] [-J JOBS] [--stage STAGE]
                     [--merge] [--no-merge] [--force] [--fast] [--overlap]
                     [--skip-s3] [--blastn-viroid] [--cleanup]
                     [--allow-duplicate-assembly] [--config CONFIG]
                     [--profile PROFILE] [--ref-db REF_DB] [--pv-db PV_DB]
                     [--id2div ID2DIV] [--rvdb-db RVDB_DB]
                     [--viroids-db VIROIDS_DB] [--diamond DIAMOND]
                     [--samtools SAMTOOLS] [--window WINDOW]
                     [--merge-distance MERGE_DISTANCE] [--evalue EVALUE]
                     [--host-bs HOST_BS] [--viral-bs VIRAL_BS]
                     [--cmd-timeout CMD_TIMEOUT] [--flank-nt FLANK_NT]
                     [--flank-merge-d FLANK_MERGE_D]

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
                        输出根目录 (默认 01_Scan; 编排器传 out/eve/01_Scan) (default:
                        01_Scan)
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
  --allow-duplicate-assembly
                        放行 assembly 指纹重复的批次登记 (默认中止; 放行也记录在案) (default: False)
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
  --flank-nt FLANK_NT   侧翼上下文长度 (nt, 两侧各扩; 0=不抽) (default: 3000)
  --flank-merge-d FLANK_MERGE_D
                        抽侧翼前先合并位点的距离 (nt) (default: 1000)
```

## `public_metadata_pipeline/eve_genome/eve_genome_scan.py`

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
                          [--cmd-timeout CMD_TIMEOUT] [--flank-nt FLANK_NT]
                          [--flank-merge-d FLANK_MERGE_D]

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
  --flank-nt FLANK_NT   侧翼上下文长度 (nt, 两侧各扩; 0=不抽) (default: 3000)
  --flank-merge-d FLANK_MERGE_D
                        抽侧翼前先合并位点的距离 (nt) (default: 1000)
```

## `public_metadata_pipeline/eve_genome/eve_scan_core.py`

```text
[docstring]
eve_scan_core.py — 内源性病毒元件 (EVE) 筛查核心逻辑
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
数据文件读写统一 utf-8 (服务器 C locale / 本机 GBK 默认编码下不炸);
外部命令支持超时 (cmd_timeout>0) 防工具挂死白占并行核.
```

## `public_metadata_pipeline/eve_genome/eve_element_db.py`

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
                         [--mafft MAFFT] [--diamond DIAMOND]
                         [--ictv-prot ICTV_PROT]

元件级 (element-level) EVE 后处理层: 跨基因组聚类 → 元件区 → consensus → 拷贝数/结构域/系统发育

options:
  -h, --help            show this help message and exit

输入/输出:
  -o OUTDIR, --outdir OUTDIR
                        eve_screen.py 的输出根目录 (也是本层输出根目录; 默认 01_Scan) (default:
                        01_Scan)
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
  --diamond DIAMOND     DIAMOND 可执行文件 (默认 PATH; 不指定则不做 reclass) (default:
                        None)
  --ictv-prot ICTV_PROT
                        ICTV 六框翻译蛋白库 (.dmnd 或 .faa); 指定后在 cluster 阶段对
                        Other_viral 簇 DIAMOND 定科, 从源头消除 Other_viral 标签
                        (default: None)
```

## `public_metadata_pipeline/eve_genome/eve_scan_report.py`

```text
usage: eve_scan_report.py [-h] -o OUTDIR [--out OUT] [--top TOP]

EVE 筛查汇总 → 自包含 HTML 报告 (在 eve_screen --merge 之后运行)

options:
  -h, --help            show this help message and exit
  -o OUTDIR, --outdir OUTDIR
                        EVE 筛查输出根目录 (含 kingdom_summary.tsv 等 --merge 产物)
  --out OUT             输出 HTML 路径 (默认: <outdir>/04_Summary/eve_report.html)
  --top TOP             家族/位点明细表展示条数 (默认: 20)
```

## `public_metadata_pipeline/eve_genome/eve_prot_lib.py`

```text
usage: eve_prot_lib.py [-h] --in FA_IN --out OUT [--min-aa MIN_AA]

eve_prot_lib.py — 把 01 阶段的核苷酸元件库翻成蛋白参考库（六框）。

为什么需要: `02`（发现管线的 EVE 判别）里 Δbit-score 的 `--eve-faa` 要的是
蛋白库（diamond makedb + blastx），而 `eve_element_db.py` 的
`eve_consensus_library.fa` 是核苷酸。这里按六框翻译、按终止符切 ORF，
只保留 >= min_aa 的片段，与 ICTV 六框库同一口径。

用法:
  python eve_prot_lib.py --in eve_consensus_library.fa       --out eve_consensus_library.sixframe.faa [--min-aa 30]

options:
  -h, --help       show this help message and exit
  --in FA_IN       核苷酸元件库 (eve_consensus_library.fa)
  --out OUT        输出蛋白 FASTA
  --min-aa MIN_AA  保留的最短 ORF 氨基酸数
```

## `virome_discovery_pipeline/eve_discriminate/eve_classify_report.py`

```text
usage: eve_classify_report.py [-h] [--table NAME=PATH] [--dir DIR] --out OUT
                              [--title TITLE] [--subtitle SUBTITLE]
                              [--notes NOTES] [--top TOP]

eve_classify_report.py — 02_classify 阶段 HTML 汇总生成器.

把这一阶段任意组合的产物汇总成一页自包含 HTML:
  * 无基因组路径: ML 打分 / Δbit-score / Nanoviridae 无监督 / 双库分类
  * 有基因组路径: 参考基因组整合证据 / 位点架构 / 转录本嵌合
  * 两张路径的产物可以同时给, 章节按传入顺序排列

为什么表用 NAME=PATH 传: 阶段产物常常改文件名或挪目录, 报告不绑定固定路径;
表型由表头自动识别 (见 eve_html.kind_of), 名字只用于章节标题。

用法:
  python eve_classify_report.py --out out/02_Classify/classify_report.html       --table "ML 打分=scores/ml_scores.tsv"       --table "Δbit-score=dbitscore.tsv"       --table "参考基因组证据=refgenome_evidence.tsv"

  # 只给一个目录时自动收 *.tsv (文件名作章节名)
  python eve_classify_report.py --dir out/02_Classify --out out/02_Classify/classify_report.html

options:
  -h, --help           show this help message and exit
  --table NAME=PATH    阶段产物表, 可重复
  --dir DIR            从目录自动收 *.tsv / *.tsv.gz (与 --table 可叠加)
  --out OUT            输出 HTML
  --title TITLE
  --subtitle SUBTITLE
  --notes NOTES        追加到页首的说明 (通常由编排器写入输入来源)
  --top TOP            每张表明细最多显示行数
```

## `virome_discovery_pipeline/eve_discriminate/no_ref_genome/eve_score_genome.py`

```text
usage: eve_score_genome.py [-h] --candidates CANDIDATES --meta META --models
                           MODELS --index INDEX --out OUT

options:
  -h, --help            show this help message and exit
  --candidates CANDIDATES
  --meta META
  --models MODELS
  --index INDEX
  --out OUT
```

## `virome_discovery_pipeline/eve_discriminate/no_ref_genome/eve_dbitscore.py`

```text
usage: eve_dbitscore.py [-h] --queries QUERIES [--labels LABELS]
                        [--label-col LABEL_COL] [--id-col ID_COL]
                        [--fam-col FAM_COL] --exo-dmnd EXO_DMND --eve-faa
                        EVE_FAA --eve-dmnd EVE_DMND
                        [--exclude-clusters EXCLUDE_CLUSTERS] [--work WORK]
                        [--threads THREADS] [--diamond DIAMOND] --out OUT

options:
  -h, --help            show this help message and exit
  --queries QUERIES     候选核酸 FASTA
  --labels LABELS       TSV: seq_id/label/family（有标签时算 AUC）
  --label-col LABEL_COL
  --id-col ID_COL
  --fam-col FAM_COL
  --exo-dmnd EXO_DMND
  --eve-faa EVE_FAA
  --eve-dmnd EVE_DMND
  --exclude-clusters EXCLUDE_CLUSTERS
                        从 EVE 库剔除的簇 ID 文件（每行一个，如 C000001）
  --work WORK
  --threads THREADS
  --diamond DIAMOND
  --out OUT
```

## `virome_discovery_pipeline/eve_discriminate/no_ref_genome/eve_nanoviridae_unsup.py`

```text
usage: eve_nanoviridae_unsup.py [-h] --candidates CANDIDATES [--ids IDS]
                                [--ref REF] [--null-cache NULL_CACHE]
                                [--mc MC] [--fdr] --out OUT

options:
  -h, --help            show this help message and exit
  --candidates CANDIDATES
  --ids IDS
  --ref REF
  --null-cache NULL_CACHE
  --mc MC
  --fdr                 批内 BH-FDR 校正（候选多时用）
  --out OUT
```

## `virome_discovery_pipeline/eve_discriminate/no_ref_genome/eve_ref_classify.py`

```text
usage: eve_ref_classify.py [-h] --query-fa QUERY_FA --eve-hits EVE_HITS
                           --virus-hits VIRUS_HITS
                           [--eve-hits-prot EVE_HITS_PROT]
                           [--virus-hits-prot VIRUS_HITS_PROT]
                           [--delta-dominant DELTA_DOMINANT]
                           [--classify-on {nt,prot}] --out OUT
                           [--min-bits MIN_BITS] [--ratio RATIO]

eve_ref_classify.py — 候选 vs EVE 元件库 / 真实病毒库 的双侧判定

背景
----
机器学习实验证明: 序列内在特征**无法**区分内源与真实病毒
(硬负样本下 AUC 0.66)。可行的替代是**同源参照**:

  * 候选命中 **已知植物 EVE 元件库** (41,933 条, 有基因组坐标背书)
    -> 与"已经整合进宿主基因组的东西"同源 = 内源旁证
  * 候选命中 **真实病毒基因组库** -> 外源病毒证据
  * 两侧都命中 -> 该家族既有内源也有外源 (如 Caulimoviridae), 无法二分

判定规则 (双侧 bitscore 比较, 与 A 线 Stage2 同一哲学):
  EVE 侧最佳 bitscore >= 阈值 且 > 病毒侧  -> 内源旁证强
  病毒侧最佳 >= 阈值 且 > EVE 侧          -> 更像外源病毒
  两侧都 >= 阈值且接近                     -> 两者兼有 (需人工)
  都不过阈值                              -> 无法判定

输出: <out>/ref_classify.tsv
  contig_id  ref_class  eve_bits  eve_id  eve_pident  virus_bits  virus_id
  virus_pident  reason

options:
  -h, --help            show this help message and exit
  --query-fa QUERY_FA
  --eve-hits EVE_HITS   候选 vs EVE 元件 blastn
  --virus-hits VIRUS_HITS
                        候选 vs 真实病毒 blastn
  --eve-hits-prot EVE_HITS_PROT
                        可选: 候选 vs EVE 元件蛋白库 diamond blastx (v7.6 Δ 占优子标签)
  --virus-hits-prot VIRUS_HITS_PROT
                        可选: 候选 vs 真实病毒蛋白库 diamond blastx
  --delta-dominant DELTA_DOMINANT
                        both 类的 Δ 占优阈值: Δ=EXO−EVE (蛋白级), |Δ|>=该值给占优子标签
  --classify-on {nt,prot}
                        分类依据: nt=核苷级 bits (原始口径); prot=蛋白级 bits (需给 --*-hits-
                        prot, 敏感度更高: 1027 候选 nt undetermined 788 -> prot 170,
                        both 5 -> 287)。给 prot 时 eve_bits/virus_bits 列仍保持核苷级口径,
                        另加 classify_basis 列标注
  --out OUT
  --min-bits MIN_BITS   两侧判定阈值 (与 A 线 Stage2 的 HOST_BS/VIRAL_BS 一致)
  --ratio RATIO         两侧接近的判据: 较弱侧 >= 较强侧 * 该值 -> both_families
```

## `virome_discovery_pipeline/eve_discriminate/no_ref_genome/eve_action_map.py`

```text
usage: eve_action_map.py [-h] [-o OUT] [--strict-eve] [--require-flank]
                         verdict_tsv

eve_action_map.py — verdict -> 可执行动作 (filter 表) [判定逻辑最后变更: 模块 v7.1].

把逐候选 verdict 表翻译成下游能直接按列过滤的动作。这一步原本写在
virome_pipeline.py 的阶段后处理里; 判别模块独立后，表也得跟着搬过来，
否则输出的 verdict 没人翻译成 action, 下游过滤还是得手工做。

输入是通用契约: 只要表里有 verdict 列（可选 evidence_tier 列）就能跑，
不绑定某个上游脚本。原先由 s3_verdict.py 产出该表（已归档到
archive/20261011_endogenous_virus_pipeline/eve_distinguish/），该步退役后
本脚本仍可用于任何带 verdict 列的判别输出。

verdict -> action 的映射是一对多的, 共 17 种 verdict (下表列常用 9 种):

  host_contamination_likely    REMOVE_host_contamination  自身寄主物种高同源, 面板适用且上游无病毒信号
  host_homology_cross_species  REVIEW                     跨物种旁证口径不作自身寄主证据
  host_homology_panel_na       REVIEW                     自身寄主高同源但面板不适用, 结构基因阴性不可评估 (v7.1)
  host_conflict_review         REVIEW                     自身寄主高同源与上游病毒注释矛盾
  missing_evidence             REVIEW                     evidence 表缺行, 上游通道未评估 (v7 P0-2)
  EVE_STRONG_provirus          MOVE_EVE                   真位点齐全 + 密码子退化 = 整合前病毒
  EVE_LTR_TE                   MOVE_EVE                   转座子域/转座子科
  ancient_EVE                  MOVE_EVE                   纯 pol 区分布式退化
  virus_candidate              KEEP_virus                 单位点齐全 + 密码子完好 = 感染中的真病毒
  *review* / review            REVIEW                     证据不足, 不得据此下 EVE 结论

--strict-eve (v7.1): evidence_tier 为 C_decay_only / D_none 的 MOVE_EVE 降为 REVIEW
—— 退化通道对随机组成的富集只有 1.21x (修复前校准), 单独不足以判 EVE; 位点级
(B) 与基因组侧翼 (A) 证据才放行。下游按 evidence_tier 列可自行分层。

输入: 带 verdict 列的逐候选 TSV (可选 evidence_tier 列)
输出: dna_vs_eve_filter.tsv (在 verdict 表后面追加 action / action_reason 两列,
 并丢掉只对调试有用的中间列, 让下游 review 时一眼能看到判据)。

v6 注: 后两条 host_* 都是 REVIEW 而不是 REMOVE。`host_contamination_likely` 是唯一
能删数据的判词, 旧版把跨物种命中和上游有病毒科注释的候选也塞进它 —— 1027 条真实候选
上 48 次调用全部打错 (全是带病毒科注释的真病毒)。删数据需要最强的证据, 拿到弱证据时
只能升级人工复核。

positional arguments:
  verdict_tsv        逐候选 verdict 表 (需 verdict 列)

options:
  -h, --help         show this help message and exit
  -o OUT, --out OUT  输出路径
  --strict-eve       严格模式: 仅凭退化通道 (C 档) 或全无证据 (D 档) 的 MOVE_EVE 降为 REVIEW ——
                     退化对随机组成富集仅 1.21x, 不足以单独判 EVE
  --require-flank    更严: 只有基因组侧翼证据 (A 档, eve_host_integration_scan) 才允许
                     MOVE_EVE; B 档 (转录本 scaffold 位点) 也降为 REVIEW。投稿口径建议开
```

## `virome_discovery_pipeline/eve_discriminate/no_ref_genome/reproduce_check_genome.py`

```text
usage: reproduce_check_genome.py [-h] [--dir DIR]

genome bundle 自检: 用包内特征实现重算 sample/ 分数, 与 expected.tsv 逐条比对并复核 16 科 AUC

options:
  -h, --help  show this help message and exit
  --dir DIR   bundle 目录 (含 sample/ models/ index.json)
```

## `virome_discovery_pipeline/eve_discriminate/with_ref_genome/eve_host_integration_scan.py`

```text
usage: eve_host_integration_scan.py [-h] [--emit-header]
                                    [blast_tsv] [s2_tsv] [host_map_tsv]
                                    [out_tsv]

eve_host_integration_scan.py — 宿主参考基因组整合扫描 [判定逻辑最后变更: 模块 v6.3]:
候选 contig blastn 回**自身物种的染色体级参考基因组**, 用整合信号区分内源 (EVE) 与真病毒.

与 eve_locus_architecture.py 的分工: 后者比对的是 OneKP **转录组组装** (碎、无染色体
上下文), 回答"候选是不是寄主基因 / 位点上组件齐不齐"; 本脚本比对的是**现成的参考
基因组**, 回答"这条候选在寄主基因组里有没有整合证据". 两条证据正交, 本脚本只在给过
`-G/--genome-manifest` 时运行; 没有参考基因组的物种自动跳过 (记 no_ref_genome),
**跳过不是"真病毒"的证据** —— 参考基因组缺失是常态, 证据是单向的:

  比对上且结构成立  -> 内源 (EVE) 方向的证据 (下面两种 call)
  比对不上/没基因组 -> 不定案, 只记列, 不改任何 verdict 方向

两种整合判据 (最终判定采信强度不同, 见 03_report/eve_final_score.py 的 rg_adjust):

  判据一 两端宿主侧翼同一位点  contig 的两端各有一段高一致宿主比对, 两段落于**同一条
                              染色体的同一条链**, 间隔不超过 MAX_FLANK_GAP —— 这是
                              "基因组片段跨着整合位点被组装出来"的铁证 (前病毒).
                              两端都有宿主段但位点对不上 -> flank_locus_conflict
                              (更像组装嵌合, 只记列不定案).
  判据二 宿主外显子嵌合        contig 内部有宿主高一致段, 呈剪接样式 (同一染色体上
                              相邻段的 subject 间隔 >= SPLICE_GAP, 即"宿主内含子被
                              剪掉了"的转座录本) 或与任一端侧翼并存 -> 这是转座录的
                              EVE 与宿主外显子拼在一起的嵌合转录本. 注意这是 contig
                              级代理: read 级的 split-splice junction 才是金标准,
                              见 README "参考基因组通道" 一节.

另有一个移交信号: 自身参考基因组聚合覆盖 >=80% 且 >=90% 一致、又无 MP/CP/AP 病毒结构
基因 (host_dominant) —— 这就是寄主序列, 交给最终判定的寄主否决三分支 (与
eve_locus_architecture.py 同一闸门, 阈值一致). 全程遵守 v6 的教训: 带 MP/CP/AP 或上游有
病毒注释的候选, 不由单一同源通道单方面判死.

输入 refgenome_blastn.tsv: run_all.sh 生成, 每行 = 4位代码 + blastn outfmt 6 的 14 列
  (qseqid sseqid pident length mismatch gaps qstart qend sstart send evalue bitscore
   qlen slen), query = 候选 contig, subject = 该物种参考基因组; 每个物种比对完追加
  一行 `#BLASTED	<代码>` 标记 (零命中也留标记, 区分"比过没打中"与"根本没比").
  行首 `#` 的其他行一律当注释跳过.

用法: python eve_host_integration_scan.py <refgenome_blastn.tsv> <domains.tsv>           <host_map.tsv> <out.tsv>
      python eve_host_integration_scan.py --emit-header
host_map.tsv: 与 run_all.sh 同一份三列映射 (样本号 <TAB> 4位代码 <TAB> 物种目录名),
  用来把候选归属到"自身物种的参考基因组"; 样本不在映射表里 -> no_ref_genome
  (rg_scope=no_sample_map), 代码不在清单/比对失败 -> no_ref_genome (genome_absent).

输出 refgenome_evidence.tsv: 逐候选 rg_call 与结构特征, 最终判定以可选第 6 参读它.

positional arguments:
  blast_tsv      run_all.sh 生成的 refgenome_blastn.tsv
  s2_tsv         域注释表 (取病毒组件)
  host_map_tsv   host_map.tsv: 样本号 <TAB> 4位代码 <TAB> 物种目录名
  out_tsv        输出 refgenome_evidence.tsv

options:
  -h, --help     show this help message and exit
  --emit-header  只打印表头并退出
```

## `virome_discovery_pipeline/eve_discriminate/with_ref_genome/eve_locus_architecture.py`

```text
usage: eve_locus_architecture.py [-h] [--emit-header]
                                 [blast_tsv] [s2_tsv] [s1_tsv] [out_tsv]
                                 [host_map_tsv]

eve_locus_architecture.py — 位点归并 [判定逻辑最后变更: 模块 v5.1]: 把候选 contig 按寄主基因组 scaffold 聚合, 在位点层面算基因座架构.

为什么需要这一步: "单一位点上 MP+CP+AP+RT+RH 齐全" 里的"位点"是**寄主基因组位点**,
不是 contig. 一个感染中的病毒在 RNA-seq 里会被打成若干个 contig, 每个只带一两个基因;
而一个整合的前病毒在位点上是完整的, 只是被组装/注释切成碎片. 所以要把映射到同一个
寄主 scaffold 的候选 contig 合并起来再看架构.

v5 修订 (2026-09-23): 两个对判定可信度影响很大的问题, 都在这一版修掉:

  1. **未映射候选不得合并成伪位点** (旧版核心缺陷). 旧代码用
     `loci[(sample, host_scaffold)]` 归组, 而没有任何寄主命中的候选 host_scaffold
     是 "-", 于是同一样本下**所有**未映射候选被塌缩成一个位点 —— 若干各带一两个
     组件的 contig 并起来凑够 >=4 组件, 就伪造出"单位点齐全", 进而伪造 EVE_STRONG
     / virus_candidate. 实测一次 1027 候选的运行: 787 条无寄主命中, 484 条落进
     多人伪位点, 39 个 locus_full 里 27 个出自伪位点. 现在未映射候选每条**自成一位点**,
     arch 记为 `no_host_locus`, 不参与位点齐全判定 (下游 final_score 的 locus_full 只看真位点).
  2. **"自身物种"口径必须可见**. 旧版在样本没有映射表条目时静默退化成"最佳 bitscore
     所在物种"的跨物种启发式, 但同源强度和 reason 文本都还写着"自身寄主物种".
     现在逐候选记录 host_scope:
       own_species         样本在映射表内且该物种有比对 -> 可当自身寄主证据
       cross_species_guess 样本无映射条目(表里没有) -> 命中来自哪个物种就算哪个, 只是旁证
       none                无任何寄主命中 -> 位点通道无输入
     最终判定的寄主否决 reason 会带上这个口径, 下游可按 host_scope 过滤.

同时这一步输出"候选 vs 寄主"的同源强度, 用来剔除把寄主基因当病毒的假阳性
(实测: 两条 27-29% 流行率的 DNA 候选就是寄主 HSP70, 靠 CDD 的 hsp70 域 PASS_VIRAL 混进来的).

输入 ( blastn outfmt 6 ):
  qseqid sseqid pident length mismatch gaps qstart qend sstart send evalue bitscore qlen slen
  (query = 候选 contig, subject = 该样本物种的 OneKP 转录组组装)

用法: python eve_locus_architecture.py <locus_blastn.tsv> <domains.tsv> <decay.tsv> <out.tsv> [host_map.tsv]
      python eve_locus_architecture.py --emit-header   # 只打印表头 (无寄主比对的运行分支复用)
host_map.tsv: 三列 "样本号 <TAB> 4位代码 <TAB> 物种目录名" (run_all.sh 用的同一份).
  给了它, 寄主同源强度只统计**自身物种**的比对 —— 跨物种旁系同源会虚高, 不能算自身寄主.
输出: locus_architecture.tsv

positional arguments:
  blast_tsv      candidate vs OneKP blastn outfmt 6
  s2_tsv         域注释表 (需 components/best_locus_comps/te_flag 列)
  s1_tsv         退化表 (可选, 需 decay_class 列)
  out_tsv        输出 locus_architecture.tsv
  host_map_tsv   host_map.tsv: 样本号 <TAB> 4位代码 <TAB> 物种目录名

options:
  -h, --help     show this help message and exit
  --emit-header  只打印输出表头并退出 (run_all.sh 的 --no-host 分支复用)
```

## `virome_discovery_pipeline/eve_discriminate/with_ref_genome/eve_transcript_chimera.py`

```text
usage: eve_transcript_chimera.py [-h] --blastn BLASTN --query-fa QUERY_FA
                                 --out OUT [--end-frac END_FRAC]
                                 [--min-mid-gap MIN_MID_GAP]

eve_transcript_chimera.py — 转录组嵌合检测 (无宿主基因组时的整合旁证)

原理
----
整合进宿主基因组的病毒片段, 其转录本两端会是**宿主外显子**, 中间是病毒序列。
这种"宿主-病毒-宿主"嵌合结构, 在没有宿主基因组时是**最接近整合位点的证据**。

判定 (基于候选 vs 宿主转录组 blastn 的 query 坐标):
  * 候选的**前 X%** 有宿主命中 (高一致性)  -> 左侧宿主段
  * 候选的**后 X%** 有宿主命中            -> 右侧宿主段
  * 中间是病毒序列 (无宿主命中)
  * 且左右两段不重叠 (中间确实有间隔)

与 eve_host_integration_scan.py 的区别
-----------------------------------
后者要求两端落在**参考基因组同一染色体同一位点** (真正的整合位点);
本脚本只能要求"两端都是宿主", 证据强度低一档 (转录本嵌合也可能是
反式剪接/组装错误), 但远强于没有任何证据。

输出: <out>/chimera.tsv
  contig_id  chimera_call  left_cov  right_cov  mid_gap  host_identity
             n_host_segs  reason

options:
  -h, --help            show this help message and exit
  --blastn BLASTN       候选 vs 宿主转录组 blastn
  --query-fa QUERY_FA   候选 fasta (取全部 ID)
  --out OUT
  --end-frac END_FRAC
  --min-mid-gap MIN_MID_GAP
                        中间间隔下限 (nt); 小于此值视为比对空档而非嵌合
```

## `virome_discovery_pipeline/eve_discriminate/eve_final_score.py`

```text
usage: eve_final_score.py [-h] [--ref-classify REF_CLASSIFY]
                          [--chimera CHIMERA] [--locus-arch LOCUS_ARCH]
                          [--multicopy-clusters MULTICOPY_CLUSTERS]
                          [--bam-coverage BAM_COVERAGE] --query-fa QUERY_FA
                          --out OUT

eve_final_score.py — 综合多路证据，输出候选最终分层

整合的五路证据
--------------
1. **同源参照** (ref_classify): 与 41,933 个已知植物 EVE 元件的同源性
2. **转录组嵌合** (chimera): 两端宿主 + 中间病毒 = 整合旁证
3. **位点架构** (locus_architecture): 同一宿主 scaffold 上组件是否齐全
4. **多拷贝家族** (mmseqs): 是否有家族扩增
5. **BAM 覆盖** (bam_coverage): 是否有活跃病毒样覆盖

评分逻辑
--------
不给"是/否"的硬判定，而是给**证据计数 + 分层**:

  tier 1 (内源强证据): >=2 条独立的内源证据
  tier 2 (内源弱证据): 1 条内源证据
  tier 3 (无证据)    : 两侧都无
  tier 4 (偏外源)    : 有外源证据且无内源证据

**重要**: 这是"证据汇总"，不是"概率"。每个 tier 都附上具体是哪几条证据，
下游可自行按需过滤。

输入可缺省
----------
`--ref-classify / --chimera / --locus-arch` 都不是必填: 没给的通道按"无证据"计。
这样三种场景都能先出表 —— 有参考基因组 (三路齐全)、只有转录组 (三路全空,
结论请看 02_classify 的序列判别)、或只跑到一半。HTML 版见
`03_report/eve_final_report.py`。

options:
  -h, --help            show this help message and exit
  --ref-classify REF_CLASSIFY
                        双库分类表; 不给则该通道记为空 (无参考基因组时常见)
  --chimera CHIMERA     转录本嵌合表; 不给则该通道记为空
  --locus-arch LOCUS_ARCH
                        位点架构表; 不给则该通道记为空
  --multicopy-clusters MULTICOPY_CLUSTERS
  --bam-coverage BAM_COVERAGE
  --query-fa QUERY_FA   候选 fasta (定义报告要列哪些候选, 不参与打分)
  --out OUT
```

## `virome_discovery_pipeline/eve_discriminate/eve_final_report.py`

```text
usage: eve_final_report.py [-h] [--final-tiers FINAL_TIERS]
                           [--channel NAME=PATH] --out OUT [--title TITLE]
                           [--notes NOTES] [--top TOP]

eve_final_report.py — 03_report 自适应最终 HTML.

为什么"自适应": 这条管线的输入是两条互斥路径 ——
  * 有宿主参考基因组: 有基因组侧翼/剪接等**整合证据**, 走 with_ref_genome
  * 只有转录组: 只有序列判别 (ML / Δbit-score / 无监督), 没有整合证据
上游还可能只跑到一半 (例如只有 ML 分, 没算 Δbit)。所以最终报告不能写死章节,
要先看**哪些通道真的产出了表**, 再决定输出什么:

  有 final_tiers.tsv  -> 证据合并报告 (tier 分布 + 证据组合 + 明细)
  没有 final_tiers    -> 明确写出缺哪几路、当前结论只基于哪些通道, 再给通道汇总

用法:
  python eve_final_report.py --out out/03_Report/eve_final_report.html       --final-tiers out/03_Report/final_tiers.tsv       --channel "ML=out/02_Classify/ML/ml_scores.tsv"       --channel "Δbit=out/02_Classify/DeltaBit/dbitscore.tsv"

options:
  -h, --help            show this help message and exit
  --final-tiers FINAL_TIERS
                        03_report/eve_final_score.py 产出的 final_tiers.tsv
  --channel NAME=PATH   某一路证据表 (与 02_classify 产物同格式), 可重复
  --out OUT
  --title TITLE
  --notes NOTES
  --top TOP
```

## `virome_discovery_pipeline/eve_discriminate/eve_bam_coverage.py`

```text
usage: eve_bam_coverage.py [-h] --bam-dir BAM_DIR --sample-map SAMPLE_MAP
                           --out OUT [--samtools SAMTOOLS] [--limit LIMIT]

eve_bam_coverage.py — 用 BAM 覆盖特征区分"整合"与"活跃病毒"

背景
----
MMPV 项目的 BAM 是**reads 比对到 5,679 条参考病毒基因组**的结果 (不是比对到
候选 contig)。所以不能直接给候选 contig 算覆盖。

但仍有可用路径: 如果一条候选能通过 blastn 映射到某个参考病毒, 就可以用
**该参考病毒的 BAM 覆盖特征**来判断这条病毒在样本里是"活跃感染"还是
"低水平/整合样":

  * 覆盖均匀性 (CV): 活跃病毒沿基因组均匀; 整合/低水平常有断崖
  * 覆盖深度: 活跃病毒高; 整合样极低
  * 链偏好: 整合位点常有链偏好 (单向转录); 活跃病毒双链均衡
  * 跨样本一致性: 同物种多样本都低覆盖 -> 更像内源

输入:
  --bam-dir    BAM 目录 (含 <sample>.sorted.bam)
  --sample-map 样本号 -> 参考病毒 accession 的映射 (来自 blastn 结果)
  --ref-fa     参考病毒 fasta (取长度)
输出: <out>/bam_coverage.tsv

options:
  -h, --help            show this help message and exit
  --bam-dir BAM_DIR
  --sample-map SAMPLE_MAP
                        两列 TSV: 样本号 <TAB> 参考病毒 accession
  --out OUT
  --samtools SAMTOOLS
  --limit LIMIT
```

## `virome_discovery_pipeline/eve_discriminate/eve_reconstruct.py`

```text
[docstring]
eve_reconstruct.py — 位点蛋白重建 (两级: miniprot 定位 → GeneWise 计数)

为什么需要这一层 (v7.2, 审稿意见): s1 的退化判定直接数"提前终止密码子", 但
终止可能来自**组装 indel** 而不是真实退化 —— 一个单碱基 indel 就让其后全部密码子
错读。hi-fever 用 GeneWise 在核苷酸层面重建蛋白、显式建模移码与内含子, 才能把
"真退化"和"组装错误"分开。本模块补上这一层。

**v7.3 改为两级架构** (借鉴 EasyPseudogene, MIT):
  1. **miniprot** 快速全基因组/位点定位 (PAF, CIGAR 的 F/N 事件) —— 便宜, 用于
     筛出"这个位点有没有同源蛋白、有没有移码/内含子";
  2. **GeneWise** 对 miniprot 给出的候选做精确重建, 在**重建后的蛋白比对串**上
     数终止 (X) 与移码 (!) —— 贵, 但计数语义可靠。
为什么必须两级: 在原始核苷酸序列上数 TAA/TAG/TGA 会被组装 indel 污染 (错读框),
而 GeneWise 的比对串里每个 X 都已核对过对应密码子 (见 `verify_gw_stops`)。
miniprot 单独用于"判定"会漏掉 GeneWise 才能看见的逐密码子事件。

后端选择:
  * `miniprot` 必有 (定位层); `genewise` 可选 (计数层)。
  * genewise 缺失 -> 用 miniprot 的 CIGAR 事件计数, 并在 status 列标
    `miniprot_only` (证据等级较低, 不伪造 GeneWise 语义)。
  * 两者都没有 -> `tool_missing`, 不伪造结果。

关于 miniprot 能否替代 GeneWise:
  * 能替代: 定位、移码/内含子事件的有无、比对质量。
  * 不能完全替代: 逐密码子的终止核验 (GeneWise 会把 X 与显示出的密码子对照)。
    因此**计数层优先 GeneWise**, miniprot 计数只作降级路径。

输出: <outdir>/05_Reconstruct/<NAME>.recon.tsv (逐位点)
  locus  backend  n_frameshift  n_stop  n_intron  aln_aa  identity  coverage
  earliest_event_frac  local_flank_identity  status
```
