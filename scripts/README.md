# scripts/ — 开发期一次性脚本归档

> 本目录收纳开发/调试/审计过程中产生的一次性脚本（2026-09-14 从仓库根目录整理归入）。
> 这些脚本**不属于核心管线**，多数硬编码了服务器路径（`/home/zhangwenda/...`）或本地结果目录，
> 仅供追溯分析过程时参考。核心管线见根目录四个 `*_pipeline/` 包及 `README.md`。
>
> ⚠️ **例外**：`db_setup/` 是**维护中的基础设施**，不是一次性脚本——部署数据库时要用它。

## 目录说明

| 子目录 | 内容 | 代表脚本 |
|---|---|---|
| `db_setup/` | **数据库部署/迁移/下载/校验（维护中，非一次性）** | `mmpv_db.py`、`mmpv_db_manifest.tsv`、`mmpv_db_fetch.tsv` |
| `blacklist/` | 宿主非植物黑名单的构建、同步与验证流水（layer1~4） | `bl_layer1_delete.py`、`bl_layer4_sdt.py`、`bl_backup.sh` |
| `db_query/` | 针对结果数据库/TSV 的临时查询脚本 | `q_host.py`、`q_mimi*.py`、`q_blacklist_candidates.py` |
| `patch/` | 已合入管线代码的一次性补丁留档 | `patch_host_blacklist*.py`、`patch_resume.py`、`rerun_sdt.py` |
| `audit/` | 结果审计、交叉核对与黑名单影响评估 | `audit_plant_fam_genera*.py`、`verify_blacklist*.py`、`scan_record_level_499.py` |
| `pilot_rvdb/` | RVDB 数据库试点分析（HMM 信号、shuffle 零假设等） | `pilot_rvdb_v3.py`、`analyze_rvdb_pilot.py`、`rvdb_hmm_per_contig.tsv` |
| `mapper_bench/` | 比对器基准测试（bowtie2/minibwa/ViReMa 对比） | `bt2_ab_test.sh`、`three_way_bench.sh`、`repro_virema.sh` |
| `utils/` | 通用小工具与服务器端启动器 | `md2docx.py`（Markdown→Word）、`run_onekp_resume.sh` |

> `db_setup/mmpv_db.py` 用法与完整部署方案见 [`../doc/DB_DEPLOYMENT.md`](../doc/DB_DEPLOYMENT.md)；
> 它是全管线数据库依赖的唯一权威清单（50 项），改库版本时需与 `pipeline_config.yaml`、
> `../doc/MMPV_dependencies.md` 的 D 节三处同步。

## 仓库完整性检查

- `../scripts/check_script_refs.py`：扫描核心管线中字符串形式引用的同仓库脚本，
  验证目标文件存在。整理/移动文件后建议运行一次。
