# MMPV — Massive Meta-mining of Plant Viruses · 平台 Wiki

> 大规模宏植物病毒挖掘与分析平台 v3.1。从公共测序档案 (NCBI SRA / CNCB GSA) 出发,
> 经 **病毒发现 → 已知病毒深度分析 → 系统发育/进化 → GenBank 提交** 的完整闭环,
> 外加 **EVE 内源性病毒筛查** 独立支线。Python 编排 + bioconda 工具 + 少量 R,
> 全部阶段带 checkpoint 与断点续跑。

## 导航

| 页面 | 内容 |
|---|---|
| [[00-架构总览]] | 七管线数据流、上下游边界契约 B1–B7、双布局 (legacy/standard)、全局机制 |
| [[01-公共数据管线]] | ① public_metadata_pipeline — SRA/GSA 检索下载 + 宿主参考建库 (8 阶段) |
| [[02-数据预处理管线]] | ② data_preprocessing_pipeline — 质控 + 去宿主去 rRNA |
| [[03-病毒发现管线]] | ③ virome_discovery_pipeline — 15 阶段 de novo 发现 (核心) |
| [[04-已知病毒分析管线]] | ④ virome_analysis_pipeline — 定量/变异/全长组装/相似性/DVG (9 阶段) |
| [[05-系统发育进化管线]] | ⑤ virome_phylo_pipeline — 18 stage / 8 大模块 (BEAST+pypopart) |
| [[06-EVE内源性病毒管线]] | ⑥ endogenous_virus_pipeline — 三阶段筛查 + 判别 (独立支线) |
| [[07-数据提交管线]] | ⑦ virome_submission_pipeline — taxonomy→.sqn (GenBank/CNCB) |
| [[08-脚本参数参考]] | 全部 467 个脚本的 `--help` 自动提取库 (按管线分组) |
| [[09-质量审查与审计]] | 审计器/单测/引用完整性检查/验收记录 |
| [[10-部署与运维]] | pixi 部署、50 项数据库、环境变量、双布局、三端同步 |

## 30 秒上手

```bash
git clone https://github.com/zhangwenda0518/MMPV-RNA.git && cd MMPV-RNA
pixi install                      # 100 个 conda 包 (手动工具见 doc/MMPV_dependencies.md)
python scripts/db_setup/mmpv_db.py check   # 数据库体检 (50 项, 必需 25)
python public_metadata_pipeline/public_data_pipeline.py --species "Lycium barbarum" \
    --taxid 112863 --stage all    # 检索→下载→转换→宿主建库→报告
python data_preprocessing_pipeline/data_preprocessing.py --stage all \
    --input_reads raw/ --host_db host_reference/hostdb/ -o out/
python virome_discovery_pipeline/virome_pipeline.py --stage all \
    --input_reads out/00b_HostDepletion --output_dir out/
```

## 平台事实卡 (以代码为准)

- **管线**: 7 条 (①–⑦) + 2 个桌面 GUI (metadata_gui / submission_gui)
- **脚本规模**: 467 个活动 Python 脚本 (本次 wiki 参数库的实际提取量)
- **阶段**: ③ 15 阶段 / ④ 9 阶段 / ⑤ 18 stage / ① 8 阶段 / ⑥ 3+判别 / ⑦ 7 项
- **断点续传**: 全部编排器 `--resume/--force`; `.ok` 标记 + 产物核验双条件
- **布局**: `legacy` (v3.0 目录名, 默认) / `standard` (管线独立根 01_PublicData…06_EVE)
- **回滚锚点**: 任何编排器均可 `--stage <子集>` 重跑, 不动上游
