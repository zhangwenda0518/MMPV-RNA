# 241 服务器从零部署与端到端运行测试报告

> 日期: 2026-09-28 凌晨 — 目标机 `10.202.2.241`（CentOS 8, 256 核 / 755G RAM / /home 146T）
> 路径: Windows 本机 → 246 (202.119.189.246, ProxyJump) → 241 (10.202.2.241)
> 部署用户: zhangwenda ｜ 项目: `~/MMPV-RNA` ｜ 测试输出: `~/mmpv_test/out_e2e`

---

## 1. 结论

**部署成功，七管线中的发现管线（discovery）已端到端跑通，且结果生物学正确。**
`mmpv_db.py verify --share ~`：**50 项依赖 OK 40，必需项仅剩 kraken2_pluspfp（解压中，见 §5）**。

合成数据端到端验证（TMV + Obuda pepper virus 两个全长病毒基因组 wgsim 模拟 14 万对
PE150 reads + 12 万条拟南芥宿主 reads，宿主库用 `build_host_pipeline.py` 现建）：

| 阶段 | 结果 |
|---|---|
| 00a Clean (fastp→seqkit→clumpify) | ✓ 560,004 reads → 279,956 |
| 00b HostDepletion (kraken2+bowtie2+ribodetector) | ✓ 宿主 reads 去除，病毒 reads 全保留 |
| 01 Assembly (megahit) | ✓ 3,482 contigs / 1.5 Mb，**含全长 TMV contig 6397/6397 bp (99.9% id)** |
| 02 Identification (10 工具) | ✓ 9 条病毒候选（viralverify=9, blast=5, metabuli=3） |
| 02b Filter (CDD+UniRef 分层) | ✓ CDD 保留 7/9（viral 6 + unclassified 1） |
| 03 COBRA / 03b Merge | ✓ merge=4（COBRA 两个已知问题见 §4） |
| 04 Cluster (CD-HIT+vclust) | ✓ 4 novel centroids + 18,629 known（plant_virus_db 正常参与） |
| 05 Taxonomy (7 工具+R 共识) | ✓ 7 条：Known=5 / NewFa=2；**Virgaviridae=2（正是 TMV/Obuda 所在科）** |
| 06 HostPrediction | ✓ Plant=4（TMV 为植物病毒，判定正确）、Bacteria=1、Animal=1、Fungi=1 |
| 07 CheckV | ✓ **Complete=2**（两个全长病毒基因组）、HQ=2/4 |
| 08 Rescue | ✓ 0 rescued（known=2，符合预期） |
| 10 Report | ✓ pipeline_report.html + Sankey 生成 |

EVE 管线单元测试（`pixi run eve-test`）：阳性/阴性对照全部通过。

## 2. 环境组成（241 实际部署形态）

- **仓库**: `git clone` 自 GitHub（246 中转 rsync；241 直连 GitHub 被间歇重置）。
  部署中发现的修复已推回仓库：`167cd6b`、`f6cc13b`、`cb20e69`。
- **主环境**: `pixi install`（`~/.pixi/bin/pixi` 0.67.2 + TUNA/USTC 双镜像）→
  `.pixi/envs/default` 11GB，48/49 个命令行工具就位。
- **独立 conda 环境**（Miniforge3，env 名与代码 `conda run -n X` 严格对齐）:
  `rdrpcatch` / `vs2`(virsorter 2.2.4) / `snpgenie` / `virhunter` / `viralm` / `wgsim`。
- **服务器侧 biosoft**（从 246 rsync，GitHub 仓库不含）:
  `biosoft/{Flye,salmon-2.5.1,ACVirus,CAPHEINE,VirSorter2-pyhmmer,rdp5}`、
  `biosoft/bin/{penguin,penguin2}`、`biosoft/virus/{CAT_pack,VCF2Dis-1.55,VCF2PCACluster}`。
- **PATH shim**（链入 `.pixi/envs/default/bin`）: `conda`、`virsorter`、`snpgenie.pl`、
  `penguin(+2)`、`VCF2Dis`、`VCF2PCACluster`；
  config 兼容路径软链: `~/mambaforge/bin/ragtag.py`、`~/biosoft/binary/diamond`、
  `~/biosoft/virus/VirBot/VirBot.py`。
- **数据库**: `~/database`（= 默认 `MMPV_DB_ROOT`，本机 146T 本地盘，非 NFS，满足 Kraken2 约束）+
  `~/plant_virus_db`。来源: 246 内网 rsync（P0-P3 优先级脚本 + 并行流），
  repo-overlay 小文件（CDD 白名单/taxid 表、cross_analysis、ICTV MSL41、VMR、genus_lens）
  只在 Windows 工作副本存在，已从本机推送到 241。

## 3. 本次实测发现并修复的 bug（均已提交 GitHub main）

| # | 位置 | 问题 | 修复 |
|---|---|---|---|
| 1 | `pixi.toml:66` | `r-data.table` TOML 裸键含点 → 解析成嵌套表，**`pixi install` 直接失败** | 键加引号（`167cd6b`） |
| 2 | `pixi.toml` | `java-openjdk` 包不存在 | 改 `openjdk` |
| 3 | `pixi.toml` | `R` 别名包不存在 | 删除（r-base 已有） |
| 4 | `pixi.toml` | `bgzip` 包不存在、`tabix` 停在 1.11 | `htslib`，删 tabix 键 |
| 5 | `pixi.toml` | `pandepth` conda 无包 | 移出手动安装 |
| 6 | `pixi.toml` | `cawlign>=1.0`/`virema>=1.0`/`drhip>=1.0` 约束超镜像上限 | 放宽；virema 包锁 bowtie≤1.0 不可解，移出（用 biosoft ViReMa 0.29） |
| 7 | `pixi.toml` | `snpgenie`(锁 perl5.22)/`rdrpcatch`(锁 py≥3.12)/`virsorter`(锁旧 ruamel) 与主环境互斥 | 移出独立 env，python 钉 `>=3.10,<3.11`（ribodetector 上限） |
| 8 | `pixi.toml` | `[feature.*]` 未挂载 → 依赖不解析不安装 | 新增 `[environments]`（default 并入 discovery/analysis/publicdata） |
| 9 | `build_host_pipeline.py` | datasets 拉该物种**全部**组装（拟南芥 374 个 >4GB） | 加 `--reference`（`f6cc13b`） |
| 10 | `build_host_pipeline.py` + `download_host_genome.py` | 合并 FASTA 写入 `#` 注释行 → **bowtie2-build 拒绝**（"Reference file does not seem to be a FASTA file"） | 元数据改旁车 `.meta.txt`（`cb20e69`） |

另发现（未改代码，见 §4）：checkpoint 失效逻辑、`mmpv_db.py` 默认 share 路径、
246 源端 `uniref90.dmnd.gz`/`hash.k2d.gz` 压缩存放（文档 §6 陷阱清单已有此条，本次实际踩中）。

## 4. 已知问题 / 行为备注（按优先级）

1. **COBRA 断点续跑把失败任务当完成**：run 中 bwa-mem2 失败的任务，下轮 resume 直接 0 任务跳过，
   无告警。需要清理 `03a_COBRA` 下任务状态后重跑。建议 resume 状态区分 failed。
2. **bwa-mem2 分派器在 241 上损坏**（"ERROR: prefix is too long!" 连 `version` 子命令都报），
   各架构子二进制（avx2/sse41/avx512bw 均 2.2.1）正常。已将 env bin 内 `bwa-mem2`
   软链到 `bwa-mem2.avx2`（CPU 支持 avx2）。⚠️ 重新 `pixi install` 会还原，建议改用
   非 dispatcher 的包 build 或上游反馈。生产 246 未复现（CPU 不同）。
3. **sample 命名约定**：编排器 clean 段接受 `{sample}_R1.fq.gz`，但 COBRA 只认
   `{sample}_1.fq.gz`。本次重命名解决。建议统一（COBRA 补 R1/R2 变体）。
4. **`mmpv_db.py verify` 默认按 `~/mmpv-db` share 布局解析**，直接放 `~/database` 时需
   显式 `--share /home/zhangwenda`。文档 §7 复现清单已有 `--share`，但默认值易误导。
5. **`filter_summary.tsv` 计数与实际不符**（显示 0/0，实际 CDD 保留 7/9）——报告聚合读的键
   与 UniProt 降级/正常两种输出的路径没对齐。 cosmetic。
6. **UniRef90 库在 246 上是 `.dmnd.gz` 压缩存放**，清单 probe 要求解压版；241 已 gunzip
   （84G）。k2 通用库同款压缩存放问题：20260226 的 `hash.k2d.gz` 解压方案**废弃**，
   k2_pluspfp 已降级为**可选**通用库（manifest required=no）：默认流程是自建宿主库
   （`build_host_pipeline.py` 四索引，E2E 已验证）。**最终决策：不部署 k2_pluspfp**
   （本次为部署测试，无需通用大库）；曾启动的 20260626 官方包下载(173GB)已中止，
   残留已清理，config/文档引用维持 0226；需要时按 `DATABASE_SETUP.md` §6 单点下载最新版。
7. **gpu env**：`pixi install` 只装 default；gpu 需 `pixi install -e gpu`（可选）。
8. **virome_phylo_pipeline（A12）**：BEAST/RDP5/PAML/TempMig 等按文档属独立部署，
   本次未部署（mafft/iqtree/treetime/R 等 pixi 侧已就位）。

## 5. 数据库清单核对（`mmpv_db.py --share ~ verify`，截至发稿）

- **必需项**: 24/24 全部就位（kraken2_pluspfp 已按决策从必需降为可选项，见 §4.6）。
- **可选项**: 就位 16/25（genomad/checkv/virsorter2/viralverify/viroids/RVDB 全家/
  vitap/acvirus/suvtk/ct3/CAT 外层/ncbi-virus/ncbi-virus_ref/plant_ref 等）。
  **nr_db 按决策不部署**（UniRef90 与 NR 二选一，默认 UniRef90，见 `DATABASE_SETUP.md` §1.2
  与 manifest 注记；半成品 152G 已清除，P0-P3 传输脚本已全部收尾）。
  其余 9 项可选缺失（`src_rvdb31`、`nt_viruses`、`RVDB-30`、`CAT-db` 内容、`ictv_nr_db`、
  `Diamond_VirusProtein_db`、`VirBot/ref`、`pyhmmer hmm` 等）经与 246 源端核对，
  多数在 246 上同样不存在——属清单"愿景项"，两端均自动降级，不阻塞任何默认管线阶段；
  需要时按 `DATABASE_SETUP.md` 对应小节单点补建。

## 6. 复现命令速查（241）

```bash
# 主链
cd ~/MMPV-RNA
pixi run python virome_discovery_pipeline/virome_pipeline.py \
    --input_reads ~/mmpv_test/reads --output_dir ~/mmpv_test/out_e2e \
    --stage all --host_db ~/mmpv_test/hostdb_build/hostdb \
    --checkv_db ~/database/virus-db/checkv-db-v1.7 --threads 32
# 校验
pixi run python scripts/db_setup/mmpv_db.py --share /home/zhangwenda verify
# 宿主库重建（新物种）
pixi run python public_metadata_pipeline/build_host_pipeline.py \
    --species "<sp>" --taxid <taxid> --stage genome-down hostdb \
    --taxonomy-dir ~/database/taxonomy --k2-libs none --work-dir ~/mmpv_test/hostdb_build
```

> 密码类凭据未写入本报告；SSH 免密链路已配（Windows ↔246 ↔241；241→246 单向不通，
> 246→241 与双向跳转可用）。
