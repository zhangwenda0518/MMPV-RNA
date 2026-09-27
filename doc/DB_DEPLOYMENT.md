# MMPV 数据库部署方案（可复现）

> 目标：把全管线依赖的数据库集中到一个**共享库根**，使 (1) 换机器/换用户可复现，
> (2) 旧路径 `~/database`、`~/plant_virus_db` 通过软链接继续可用——核心代码里那些
> 硬编码路径**一行不用改**。
>
> 配套工具：`scripts/db_setup/mmpv_db.py`（部署/迁移/下载/校验）
> + `mmpv_db_manifest.tsv`（50 项依赖清单）+ `mmpv_db_fetch.tsv`（25 项获取配方）。
> 每库的下载与建索引原文命令见 `DATABASE_SETUP.md` §1–§6。
> 本文与 `doc/MMPV_dependencies.md` 的 **D 节（数据库/参考文件）** 和 **「路径配置」节** 配套：
> 那里按用途分组讲"有哪些库"，这里按部署单元讲"怎么装到新机器上"。

---

## 0. 结论：可以在服务器单独建 database 目录吗？

**可以。** 但三个前提，缺一会踩坑：

| # | 前提 | 为什么 |
|---|---|---|
| 1 | 旧路径必须留软链接 | 配置层（`pipeline_config.yaml`）虽然支持 `${MMPV_DB_ROOT:-...}`，但**核心代码有 6 处绕过配置硬编码 `~/database/...`**，其中 `filter_virus.py` 的 CDD 三件套是**当前实际生效**的路径（编排器根本不传 CLI 覆盖）。详见 §5。 |
| 2 | 软链接指向**绝对路径** | 相对软链接会随链接所在目录失效。`mmpv_db.py link` 统一建绝对链接。 |
| 3 | 共享盘不能是网络文件系统（至少 Kraken2 库不能） | `k2_pluspfp` 约 50GB，靠 page cache 加速；NFS 上会让去宿主阶段慢到不可用。其余库放网络盘可以。 |

反向来说**不需要**的前提：不需要改 `pipeline_config.yaml`、不需要改任何 Python 代码、
不需要重新下载数据库。数据搬到新位置后旧路径变成软链接，所有消费者照常工作。

---

## 1. 依赖全景

共 **50 项**（必需 25 / 可选 25），全量约 **710GB**；最小必需集约 **250GB**
（不含 nr 库 250GB、CAT 库 40GB、ncbi-virus_ref 48GB、现成 vConTACT3 14GB 等可选大件）。
其中 **25 项有自动化获取配方**（`down` 能直接下），24 项需共享池或人工，1 项随 git 自带；分类见 §3.5。

```bash
python scripts/db_setup/mmpv_db.py list    # 按类别打印 50 项
```

| 类别 | 数量 | 含义 | 复现方式 |
|---|---|---|---|
| `download` | 7 | 一条命令下预建库 | `down` 可自动 |
| `build` | 22 | 下原始数据后本地建索引 | 部分 `down` 可自动；CDD/genus_lens 等需共享池 |
| `local` | 20 | **无法重下**，只能从已有机器/共享池拷贝（含 5 个被 gitignore 的 `repo_*` 小文件） | `rsync` / `adopt --strategy link` / `link` |
| `repo` | 1 | 随 git 仓库自带，无需部署 | 无需处理（列出以保证清单完整） |

> ⚠️ `local` 类最容易被"复现方案"漏掉，多为第三方网页下载（CAT_pack、VITAP、PhaBOX2、
> vConTACT3、ncbi-virus_ref、ViraLM 模型、VirBot 参考）与被 gitignore 的仓库侧小文件。
> **共享池必须包含它们**，否则新用户永远配不齐。

---

## 2. 目标布局

```
${MMPV_DATA_SHARE}/                      # 共享库根，例: /data/mmpv
├── database/                            # ← ${MMPV_DB_ROOT}
│   ├── virus-db/
│   │   ├── RVDB-v31/                    # ← RVDB 派生的整套库都在这（见 §3.2）
│   │   │   ├── RVDB_viroids.diamond_db/U-RVDBv31.0-prot_unique.dmnd
│   │   │   ├── RVDB.mmseqs_db/RVDB.mmseqs
│   │   │   ├── RVDB_viroids.metabuli_db/
│   │   │   ├── CAT-db/{db,tax}          # CAT_pack 产出（无前缀）
│   │   │   └── U-RVDBv31.0-prot.{hmm,EX.acc.fasta,info.tab}
│   │   ├── RVDB-30/RVDB.mmseqs          # 与 v31 并存（YAML 与代码指向不同版本）
│   │   ├── genomad_db/  checkv-db-v1.7/  virsorter2_db/
│   │   ├── viralverify_db/nbc_hmms.hmm  ViralVerify/nbc_hmms.h3m   # 两个都留，见 §5
│   │   ├── viroids-db/  vitap-db/  acvirus_db/  vConTACT3_db/  ct3_DBs/
│   │   ├── ictv-db/  ncbi-virus_ref/  ncbi-virus/  Diamond_VirusProtein_db/
│   │   ├── db/genus_lens + genus_lens_no_prefix.json
│   │   ├── hmm/{viral/combined.hmm,pfam/Pfam-A-*.hmm}   # pyhmmer
│   │   ├── viralm_db/  VirBot/ref/  phabox_db_v2_2/  suvtk_db/  rdrp-db/
│   │   └── VMR_MSL41.v1.20260320.xlsx
│   ├── taxonomy/                        # taxdump + new_taxdump 产物 + accession2taxid
│   ├── cdd/cdd-db/                      # mmseqs 库（cdd_db 是文件前缀，单元是它所在目录）
│   ├── kraken2/k2_pluspfp_20260226/     # ~50GB，别放网络盘
│   ├── uniport_db/uniref90/             # 目录名拼写就是 uniport_db（历史笔误，勿改）
│   ├── nr_db/nr.dmnd
│   ├── nt-db/nt_viruses/
│   └── host_db/{kraken2,bowtie2,hisat2,minimap2}
├── plant_virus_db/                      # ← ${MMPV_PLANT_VIRUS_DB}
│   └── 3.final-ref-virus.db/
├── src/                                 # ← ${MMPV_SRC} 上游原始数据（RVDB fasta 等），建完索引可删
│   └── RVDB/RVDB_v31.fasta.gz
├── repo-overlay/                        # 仓库侧小文件单一副本（CDD 白名单、virus.taxid.txt、
│                                        # genus_lens、ICTV 表、cross_analysis 宿主概率表）
├── tmp/                                 # 大件下载临时区（支持断点续传）
└── logs/down/                           # down 的逐库日志
```

---

## 3. 部署路线

### 3.1 三条路线总览

| 路线 | 适用 | 命令序列 |
|---|---|---|
| **A. 本机已有数据库 → 迁移** | 当前服务器，想集中到独立大盘 | `check` → `layout --apply` → `adopt --apply` → `link --apply` → `verify` |
| **B. 全新机器 → 下载** | 换服务器、首次复现 | `layout --apply` → `down --apply` → `link --apply` → `verify` |
| **C. 已有共享池 → 挂载** | 多人共用 | `link --apply` → `env --write` → `verify` |

### 3.2 RVDB 派生链（路线 B 的主力）

`virome_discovery_pipeline/build_virus_db.py` 是仓库自带的**统一建库器**，接受上游 fasta，
产出 18 种软件的库（metabuli centrifuger kraken2 bracken krakenuniq kmcp ganon sylph
kun_peng blast kma salmon_kallisto lexicmap kaiju kraken2x diamond mmseqs CAT）。

```
RVDB_v31.fasta.gz（上游基底，条目 src_rvdb31）
      │
      ▼  build_virus_db.py  --work-dir $VIRUS_DB/RVDB-v31 --tax-dir $DB/taxonomy
      ├── --db-prefix RVDB_viroids  --databases diamond    →  RVDB_viroids.diamond_db/
      ├── --db-prefix RVDB_viroids  --databases metabuli   →  RVDB_viroids.metabuli_db/
      ├── --db-prefix RVDB_viroids  --databases CAT        →  RVDB_viroids.CAT-db/
      └── --db-prefix RVDB          --databases mmseqs     →  RVDB.mmseqs_db/RVDB.mmseqs
```

**产出名有两处与管线期望不一致**，配方里已补 `ln -sf` 别名（这是实测踩到的坑）：

| 管线期望 | 建库器产出 | 处理 |
|---|---|---|
| `RVDB_viroids.diamond_db/U-RVDBv31.0-prot_unique.dmnd` | `RVDB_viroids.diamond_db/RVDB_viroids.dmnd` | `ln -sf RVDB_viroids.dmnd U-RVDBv31.0-prot_unique.dmnd` |
| `RVDB-v31/CAT-db/` | `RVDB-v31/RVDB_viroids.CAT-db/` | `ln -sfn RVDB_viroids.CAT-db CAT-db` |
| `RVDB.mmseqs_db/RVDB.mmseqs` | 同（前缀取 `RVDB`） | 无需处理 |
| `RVDB_viroids.metabuli_db/` | 同（前缀取 `RVDB_viroids`） | 无需处理 |

> 运行前提：`build_virus_db.py` 的 `check_dependencies()` 强制要求 3 个自写辅助脚本
> （`SearchAccessionIdToTaxId.py` / `db_seqid2taxid_add_legth.py` / `make_ktaxonomy.py`）
> 和 `aria2c rapidgzip seqkit awk sed`。辅助脚本**目前不在仓库**——见 §3.5。

### 3.3 路线 A：既有机器迁移

```bash
cd ~/MMPV-RNA

python scripts/db_setup/mmpv_db.py --share /data/mmpv check                 # ① 体检
python scripts/db_setup/mmpv_db.py --share /data/mmpv layout --apply        # ② 建骨架

python scripts/db_setup/mmpv_db.py --share /data/mmpv adopt \               # ③ 演练后加 --apply
    --db-from ~/database --plant-from ~/plant_virus_db
python scripts/db_setup/mmpv_db.py --share /data/mmpv adopt \
    --db-from ~/database --plant-from ~/plant_virus_db --apply

python scripts/db_setup/mmpv_db.py --share /data/mmpv link --apply          # ④ 旧路径软链接
python scripts/db_setup/mmpv_db.py --share /data/mmpv publish --apply       # ⑤ 仓库小文件入池
python scripts/db_setup/mmpv_db.py --share /data/mmpv env --write /data/mmpv/mmpv-db.env
echo 'source /data/mmpv/mmpv-db.env' >> ~/.bashrc
python scripts/db_setup/mmpv_db.py --share /data/mmpv verify                # ⑥ 校验
```

`adopt` 行为要点：同文件系统用 `mv`（秒级），跨文件系统自动切 `rsync -a`（可续传）；
**递归合并**（先 `layout` 建的空骨架会就地让位）；**绝不覆盖已有非空数据**；
迁移后旧目录若已空，自动 `rmdir` 并换成软链接。

### 3.4 路线 B：全新机器下载

```bash
python scripts/db_setup/mmpv_db.py --share /data/mmpv layout --apply
python scripts/db_setup/mmpv_db.py --share /data/mmpv down                # 演练：看会下哪些
python scripts/db_setup/mmpv_db.py --share /data/mmpv down --apply        # 实际下载/构建
python scripts/db_setup/mmpv_db.py --share /data/mmpv verify              # 看还缺什么
```

`down` 的特点：
- **默认演练**，必须 `--apply` 才落盘；已就位的条目自动跳过（`--force` 可强制重跑）。
- 逐条检查**工具前置**（缺 `genomad`/`diamond` 就跳过并说明）与**先决条目**（`nr_db` 等
  `taxid_map_prot` 就位）。
- 跑完按清单哨兵**复验**：命令返回 0 但产物不对会报 `FAILED`，不会假成功。
- 大件落 `<share>/tmp`（在共享盘上，空间够），逐库日志在 `<share>/logs/down/`。
- `--only ID...` 单跑、`--kind download|build` 分类跑、`--jobs N` 并行、`--threads N` 建索引线程。
- `--fetch-file` 可换自定义配方（站点自维护 / 测试用）。

### 3.5 ⚠️ 当前复现阻塞点（必须知道）

**24 项没有配方**，`down` 会逐条列出并给出替代路径。但「需共享池」这个说法太粗——
按**真实原因**分四类，其中只有 C 类是真的拿不到：

| 类 | 数量 | 条目 | 性质与处置 |
|---|---|---|---|
| **A. 自家数据 / 自家服务器文件** | 7 | `plant_ref_db`(你们的植物病毒数据集)、`repo_cdd_whitelist`、`repo_cdd_taxid`、`repo_virus_taxid`、`repo_genus_lens`、`repo_ictv_msl41`、`repo_cross_analysis` | **对当前服务器根本不算缺**，只是别人 clone 不到。处置：`git add -f` 直接入库（都是小文件，最干脆），或 `publish` 入池 |
| **B. 被缺失脚本卡住** | 4 | `cdd_db`、`virus_db_genus_lens`、`virus_db_genus_lens_json`、`genus_family_ref` | 建库脚本只在服务器（`utils/db_build/` 待收编）。**收编 5 个脚本即全部解除**，同时解锁 §3.2 整条 RVDB 派生链 |
| **C. 第三方网页下载，无稳定命令行** | 11 | `virus_db_vitap`、`virus_db_vcontact3`、`virus_db_ncbi_ref`、`virus_db_ncbi_virus`、`virus_db_phabox2`、`virus_db_virbot_ref`、`virus_db_viralverify`、`virus_db_viralverify_alt`、`virus_db_rvdb31_hmm`、`virus_db_pyhmmer_hmm`、`virus_db_ictv` | **真要共享池**。多为 bitbucket / OneDrive / Google Drive / 作者自建库，无命令行获取方式 |
| **D. 需 per-run 参数，非「下不到」** | 2 | `host_db`（`build_host_pipeline.py --species X --taxid Y`，每个宿主物种一份）、`nt_viruses`（从 nt 抽病毒子集，或 `build_virus_db.py --databases blast`） | 工具都在仓库里，只是参数随项目变，不适合写死成配方 |

> 由此修正一个更准的结论：**换新服务器时，真正"必须从共享池拿"的只有 C 类 11 项**；
> A 类是自产文件（建议直接入库）、B 类收编脚本即可、D 类是参数化构建。
>
> 与 `down` 输出的关系：`down` 按 **kind（获取方式）** 分组打印，上表按 **缺失原因** 分组，
> 两者互补——同一个条目在 `down` 里可能显示为 `[local]`，但在上表里属于 A 类（自产）
> 或 C 类（外部网页）。
>
> 底线结论：**全新机器目前无法只靠仓库把 50 项配齐**——即便收编脚本、入库自产文件之后，
> C 类 11 项仍需共享池先由已有机器 `publish` + `rsync` 一次。这不是方案缺陷，
> 是上游依赖本身的形态；如实标注优于假装能下载。

### 3.6 路线 C：挂共享池

```bash
python scripts/db_setup/mmpv_db.py --share /mnt/mmpv-share link --apply    # 建 ~/database 软链接
python scripts/db_setup/mmpv_db.py --share /mnt/mmpv-share env --write ~/mmpv-db.env
source ~/mmpv-db.env
python scripts/db_setup/mmpv_db.py --share /mnt/mmpv-share verify
```

只读挂载也没问题：`link` 建的是**指向**挂载点的链接，不写挂载点本身。

---

## 4. 环境变量

`mmpv_db.py env` 会生成数据库相关的全部变量。完整清单与语义（含 phylo / submission 管线的
变量）见 `doc/MMPV_dependencies.md` 的「环境变量清单」节，这里只列数据库相关的：

**第一类：`pipeline_config.yaml` 读取（决定配置层解析结果）**

| 变量 | 默认 | 说明 |
|---|---|---|
| `MMPV_DB_ROOT` | `/home/zhangwenda/database` | **总开关**，覆盖它即迁移整批数据库 |
| `MMPV_VIRUS_DB` | `{DB_ROOT}/virus-db` | |
| `MMPV_HOST_DB` | `{DB_ROOT}/host_db` | |
| `MMPV_KRAKEN2_DB` | — | Kraken2 预建索引 |
| `MMPV_CHECKV_DB` | — | 仅 `plant` profile 用 |
| `MMPV_PLANT_VIRUS_DB` | `/home/zhangwenda/plant_virus_db` | 独立根 |
| `MMPV_SRC` | `{SHARE}/src` | 上游原始数据（本方案新增） |

**第二类：核心代码绕过配置直接读（不设就回落到 `~/database` 硬编码）**

| 变量 | 对应硬编码点 |
|---|---|
| `MMPV_VMR` | `utils/annotate_nucleic_acid.py:116` |
| `MMPV_GENUS_FAMILY_REF` | `virome_pipeline.py:1843`、`virus_classifier_analysis.R:36` |
| `MMPV_NT_DB` / `MMPV_AA_DB` | `gen_final_judgement.py:300/303` |
| `MMPV_CDD_TAXID_TSV` | `gen_final_judgement.py:44` |
| `UNIPROT_DB` / `VIRUS_TAXID` | `filter_virus.py:24-25`（`UNIPROT_DB` 要**不带 `.dmnd` 后缀**，代码自己拼） |
| `MMPV_T2ASN` | `virome_submission_pipeline/sync_sqn_from_csv.py:43` |

**第三类：只影响本工具**：`MMPV_DATA_SHARE`（共享库根）、`MMPV_BASH`（显式指定 bash；
Windows 上 PATH 里的 `bash.exe` 可能只是 WSL 安装提示桩）。

---

## 5. 兼容软链接层（为什么必须，覆盖了哪些）

`pipeline_config.yaml` 的路径全部包在 `${...:-...}` 里，理论上设个 `MMPV_DB_ROOT` 就迁移完了。
**但实际上有 6 处代码绕过了配置**，所以必须靠软链接兜底：

| # | 位置 | 硬编码内容 | 严重度 |
|---|---|---|---|
| 1 | `filter_virus.py:21-25` | `DEFAULTS` 里的 `CDD_DB` / `CDD_WHITELIST` / `CDD_CLASSIFIED` / `UNIPROT_DB` / `VIRUS_TAXID` | **最高**：`virome_pipeline.py` 传参时**没传** `--cdd-db/--cdd-whitelist/--cdd-classified`，所以这就是 02b_Filter 实际生效的路径 |
| 2 | `virus_classifier.py:757-775, 863` | 一整套并行 DB 树：`RVDB-v31/CAT-db`、`vitap-db/VMR-MSL40_DB`、`RVDB-v31/RVDB.mmseqs_db`、`acvirus_db`、`vConTACT3_db`、`genomad_db`、`metabuli` | 高：`db_paths` 的键只有 CLI 显式传才生效，编排器只传 `--db-dir` |
| 3 | `virus_identification.py:70, 327-390` | `DEFAULT_DB_DIR = ~/database/virus-db` **及全部子目录布局** | 高：单脚本独立运行时生效 |
| 4 | `validate_rescue_cdd.py:47-51` | CDD 三件套 + `plant_virus_db` 的 genera/family 表 | 高 |
| 5 | `virome_pipeline.py:52, 1489, 3019, 3426` | `taxonomy/genus_family_ref.tsv`、`virus-db/db/genus_lens`、`viroids.fasta.blast.db`、`VMR_MSL41*.xlsx` | 中 |
| 6 | 其他 | `hmm_ct3_evidence.py:83`、`run_host_prediction.py:814`、`virome_analysis.py:205-206`、`sdt_genus_matrix.py:218`、`integrate_rescue_evidence.py:206-207`、`panvirome_novel_known.py:49-58` | 中 |

另外 `virome_analysis_pipeline/auto_known_virus.py:268-313` 用**裸 `yaml.safe_load`，不做 `${}` 展开**——
所以分析管线**完全不受 `MMPV_DB_ROOT` 影响**，只能靠软链接。

**结论：软链接不是可选优化，是这套方案能成立的前提。**

两处 YAML 与代码的路径分歧也已固化在清单里，两个位置都保留（避免任一方代码改动后断链）：

| 键 | YAML 指向 | 代码指向 |
|---|---|---|
| `viralverify_hmm` | `virus-db/viralverify_db/nbc_hmms.hmm` | `virus-db/ViralVerify/nbc_hmms.h3m` |
| `mmseqs_db` | `virus-db/RVDB-30/RVDB.mmseqs` | `virus-db/RVDB-v31/RVDB.mmseqs_db/RVDB.mmseqs` |

---

## 6. 校验与排障

```bash
python scripts/db_setup/mmpv_db.py verify                         # 逐项表格，缺必需项退出码 1
python scripts/db_setup/mmpv_db.py verify --only cdd_db kraken2_pluspfp
python scripts/db_setup/mmpv_db.py which virus_db_genomad         # 单查
python scripts/db_setup/mmpv_db.py fetch virus_db_rvdb31_diamond  # 看某个库的配方
python scripts/db_setup/mmpv_db.py check                          # 含旧路径软链接健康度
```

`verify` 的状态值不只有「缺失」，区分这几种——**迁移后最常出的就是断链**：

| 状态 | 含义 | 处置 |
|---|---|---|
| `OK` | 存在且非空 | — |
| `缺失` | 路径不存在 | `down --only <id> --apply` 或从共享池拷 |
| `断链` | 软链接指向的目标不存在 | 共享盘没挂载 / 被挪走 → 检查挂载点 |
| `空目录` | 目录存在但没有内容 | `layout` 建的空壳，数据还没灌进去 |
| `空文件` | 文件存在但 0 字节 | 下载中断，重下 |
| `类型不符` | 声明 dir 实为 file（或反之） | 清单写错了，或路径被错误地建成了目录 |

**典型故障**：

- **`~/database/...` 全部报不存在，但共享根里明明有** → 旧路径软链接悬空。`check` 会明确报出
  「软链接 → 目标缺失（断链！）」。多为共享盘未挂载。
- **`down` 全部 SKIP，说缺 `genomad`/`diamond`** → 先装工具，见 `DATABASE_SETUP.md` §8。
- **`down` 说 `utils/db_build 缺 …`** → §3.5 的硬阻塞，先把 5 个脚本收编进仓库。
- **CDD 过滤阶段报找不到 `cdd_virus_final_v4.txt`** → 该文件被 `.gitignore` 忽略，`git clone`
  拿不到。从共享池 `link --apply` 挂进来。
- **`genus_lens` 被当成目录** → 它是 TSV **文件**（`integrate_rescue_evidence.py:212` 用 `isfile()`）。
  清单里已声明为 `file`；若手工建成目录会静默跳过长度门槛。
- **`cdd_db` 找不到** → `cdd_db` 是 mmseqs 的**文件前缀**，真正要保证的是它的**父目录**
  `cdd/cdd-db/` 里整套 `cdd_db*` 文件齐全。清单用「目录 + 哨兵 `cdd_db`」两级校验覆盖。
- **`fullnamelineage.dmp` 找不到** → 它**不在** `taxdump.tar.gz` 里，来自 `new_taxdump.tar.gz`。
  配方已处理（同时取 `rankedlineage.dmp` 供 `build_genus_family_ref.py` 用）。

---

## 7. 新用户复现清单（copy-paste）

```bash
# ── 0. 前置：拿到共享池路径（向维护者索取）
export MMPV_DATA_SHARE=/data/mmpv

# ── 1. 克隆仓库
git clone <repo> ~/MMPV-RNA && cd ~/MMPV-RNA

# ── 2. 环境变量
python scripts/db_setup/mmpv_db.py --share $MMPV_DATA_SHARE env --write ~/mmpv-db.env
source ~/mmpv-db.env && echo 'source ~/mmpv-db.env' >> ~/.bashrc

# ── 3. 旧路径软链接（核心代码的硬编码路径靠这个活）
python scripts/db_setup/mmpv_db.py --share $MMPV_DATA_SHARE link --apply

# ── 4. 补齐本机缺的库（能自动下的一批）
python scripts/db_setup/mmpv_db.py --share $MMPV_DATA_SHARE down
python scripts/db_setup/mmpv_db.py --share $MMPV_DATA_SHARE down --apply

# ── 5. 校验；剩余缺项按 §3.5 从共享池补
python scripts/db_setup/mmpv_db.py --share $MMPV_DATA_SHARE verify

# ── 6. 工具环境（与数据库无关，独立一步）
pixi install          # 或按 DATABASE_SETUP.md §8 逐个 pixi global install
```

---

## 8. 维护约定与已知陷阱

**改版本时三处要同步**，否则会漂移：
1. `scripts/db_setup/mmpv_db_manifest.tsv` 的 `path`/`probe`
2. `pipeline_config.yaml` 的对应键
3. 本文 §2 的目录树、§1 的类别表

**陷阱清单**：

- **别把 `uniport_db` 改成 `uniport`/`uniprot_db`**。历史笔误，代码和 YAML 都用 `uniport_db`，改了必然找不到库。
- **`~` 是按用户展开的**。`~/database` 软链接只对建它的那个用户有效。多人共用或提交 HPC 作业时，
  必须让作业环境 `source` 共享的环境变量文件，不能指望 `~`。
- **软链接可以是跨文件系统的**（硬链接不行）。数据放独立大盘、链接留在 `~`，完全没问题。
- **Kraken2 库别放网络盘**（见 §0 前提 3）。
- **`taxid_map_prot` / `taxid_map_nucl` 若保留 `.gz`**，`probe` 会报缺失——`diamond --taxonmap` 需要解压后的文件。
- **blast 库是文件组**（`xxx.blast.db.{nhr,nin,nsq}`），没有单一哨兵文件，所以清单里只做目录非空校验。
- **`database/*` 被 `.gitignore` 忽略**（只有 `final.cluster.ref.fasta` 例外）。所以仓库里的
  `database/cdd/*.txt`、`database/uniprot/virus.taxid.txt`、`database/cross_analysis/`、`database/genus_lens`
  在别人 clone 出来是**空的**，这是 `repo-overlay` 机制存在的原因。
- **`deplete:` / `db_dir` / `blast_db` 等配置键是死键**（`virome_pipeline.py` 读了但不传给子进程），
  迁移时别指望这些键，靠软链接。
- **YAML 浅合并陷阱**：给 profile 加锚点继承时嵌套段要各自加锚点，否则子键会整体丢失
  （详见 `doc/MMPV_dependencies.md`）。

---

## 附：工具子命令速查

| 命令 | 作用 | 默认是否落盘 |
|---|---|---|
| `list` | 打印 50 项依赖清单 | — |
| `check` | 部署前体检（可写性、旧路径现状、缺项统计） | 只读 |
| `verify` | 逐项校验，缺必需项退出码 1 | 只读 |
| `which <id>...` | 单查若干 id 的解析路径与状态 | 只读 |
| `fetch <id>...` | 看某个库的获取方式与配方命令 | 只读 |
| `env [--write F]` | 生成环境变量块 | 只读 |
| `down` | 按配方下载/构建（`--only/--kind/--jobs/--threads/--force`） | 需 `--apply` |
| `layout` | 建规范目录骨架 | 需 `--apply` |
| `adopt --db-from .. [--plant-from ..]` | 迁移/链接既有数据库树（`--strategy move\|link`） | 需 `--apply` |
| `link` | 建旧路径 + 仓库侧兼容软链接 | 需 `--apply` |
| `publish` | 仓库侧小文件发布到共享池 | 需 `--apply` |

所有写操作默认**演练**，必须显式 `--apply`；`adopt` 永不覆盖已有非空数据。
