# teagcn/ — vendored TEA-GCN 引擎（来源与改动记录）

上游: https://github.com/pengkenlim/TEA-GCN （MIT License, 见 LICENSE）
**上游锁定 commit**: `a953ac36cd6af9ccace82760343089e7dce70583`（= tag `v1.0.1`，亦为 main HEAD；提交时间 2026-05-11T03:57:50Z）。

抽取范围核对结论（**2026-10-09 重新复核**：GitHub trees API 取全树 + 逐文件 git-blob SHA-1 比对）：

1. 与上游**共有的 13 个文件**确实逐字节一致（`LICENSE`、`coexpression/{bicor,ensemble,pearson,spearman,TEA_fast}.py`、`data_processing/{expression_matrix,network,read_write}.py`、`main/{Generate_partitions,Run_TEA-GCN,Run_TEA-GCN_fast,Run_TEA-GCN_mp}.py`、`setup/requirements.txt`、`sklearn_wrappers/{kmeans,PCA}.py`）；
2. 其中 **3 处带已声明补丁（涉 4 个文件）**：`coexpression/{TEA.py, TEA_mp.py, TEA_fast.py}`（补丁 1、3）与 `main/Rank_transform.py`（补丁 2），详见下文；
3. 另有 **4 个本地新增文件在上游全历史中不存在**（`commits?path=…` API 均返回空数组）：
   `coexpression/TEA_mp.py`、`coexpression/TEA_fast.py`、`main/Run_TEA-GCN_mp.py`、`main/Run_TEA-GCN_fast.py`。
   其中 `TEA_mp.py` / `Run_TEA-GCN_mp.py` 是本项目**默认**并行入口，`TEA_fast.py` / `Run_TEA-GCN_fast.py` 为线程池版。

> ⚠️ **勘误（2026-10-09）**：本文件早先写作「未从 goji 工作副本带入任何上游没有的改动」，与上面第 3 条直接矛盾（`README.md §7` 也明说 mp/fast 是「你们的入口」）。**以本节为准。**
>
> 未在「除外项」中点名、但同样未 vendor 的上游文件：`README.md`、`data_processing/metadata.py`、`data_processing/nlp.py`、`__pycache__/*.pyc`。
>
> 本次复核同时确认：goji 工作副本（`D:\桌面\差异表达分析\TEA-GCN-main`）对上游 main 为 `same=46 diff=0 missing=0`，即上游 main 的忠实快照。

## 抽取范围

仅保留 coexp 阶段所需的引擎代码；**未**包含 colab_data/（~58MB 演示数据）与
下游 Colab 模块（Download_metadata / Annotate_partitions / Prepare_neighbourhoods /
Get_Partition_Rankings —— 这些走 ENA API 与 Colab notebook，与 Δz 分析无关）。

## 在上游基础上应用的两个补丁（均验证输出不变）

1. **死代码跳过**（`coexpression/TEA.py`、`coexpression/TEA_mp.py` 各 1 行）：
   `full=False`（默认）时，逐对基因的分区相关系数字符串拼接结果**从不写入输出
   文件**（上游 full 分支中带 cor 列的写入行本身是注释），但该拼接实测占单基因
   耗时 ~80%（39,224 基因 × k=19 微基准：每基因 311ms 中 249ms）。补丁后单网络
   串行 3.4h → 约 40min。**输出逐字节不变**（300 基因沙箱 diff 验证）。
2. **Windows 路径分割修复**（`main/Rank_transform.py` 2 处）：
   `network_path.split("/")[-1]` → `os.path.basename(...)`。上游在 Windows 下会把
   Completed_GCN 写到错误嵌套路径；Linux 上两者等价，行为不变。
3. **`--resume` 断点续跑会漏算基因**（`coexpression/TEA.py`、`TEA_mp.py`、`TEA_fast.py`；2026-10-09 新增）：
   入口脚本用 `offset = len(os.listdir(network_path))`（已完成文件**个数**）作游标，
   `calc_untargeted` 按 `gene_idx >= offset` 跳过——这假定「已完成 = 基因顺序前缀」。
   对串行版成立；但 `TEA_mp.py`（`pool.imap_unordered`）与 `TEA_fast.py` **乱序写盘**，
   文件个数 ≠ 前缀长度 → 已算完的基因可能被跳过、**永久缺失且不报错**。
   补丁把跳过条件改为**按文件名集合**判断（`done = set(os.listdir(network_path))`；
   `if gene not in done`），`offset` 仅作为「是否续跑」的开关；非续跑（offset=0）行为不变。
   **单元测试**：预置 A、C 已完成、B 缺失且 offset=2 —— 旧逻辑只处理 C（B 永久丢失），
   新逻辑恰好只处理 B，A/C 不重算。
   注：未触发过已知事故（历次运行均完整产 39,224 个文件）；属潜伏缺陷修复。

## 并行入口与结果等价性

| 入口 | 机制 | 平台 | 等价性证据 |
|---|---|---|---|
| `main/Run_TEA-GCN_mp.py`（默认） | fork 进程池，子进程 BLAS 单线程 | 仅 Linux | goji 真数据 39,224/39,224 文件 md5 与串行一致 |
| `main/Run_TEA-GCN_fast.py` | 线程池（受 GIL 限制，加速有限） | 跨平台 | 300 基因沙箱与串行逐字节一致 |
| `main/Run_TEA-GCN.py` | 串行（基准/兜底） | 跨平台 | 上游原版 |

mp 入口依赖 `multiprocessing.get_context("fork")`，Windows 不可用（编排器会前置
检查并提示换 `--teagcn-entry serial|fast`）。MMPV 目标平台为 linux-64，默认 mp。

## 环境依赖

`setup/requirements.txt` 为上游 pin（numpy==2.0.1 / pandas==2.2.2 / scipy==1.14.0 /
scikit-learn==1.5.1 / einsumt==0.9.4）。MMPV 主 pixi 环境的版本范围与其全部兼容，
仅缺 `einsumt`（conda 无包）——已加入 pixi.toml 的 `[pypi-dependencies]`。
注意：在非 pin 环境重跑可能与 goji 原结果存在浮点级差异（rank 并列翻转），
属环境固有差异；首次真机运行建议做一次串行 vs mp 的 md5 冒烟对比。
