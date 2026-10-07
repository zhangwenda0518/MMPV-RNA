# teagcn/ — vendored TEA-GCN 引擎（来源与改动记录）

上游: https://github.com/pengkenlim/TEA-GCN （MIT License, 见 LICENSE）
上游 main 最后提交: 2026-05-11。本目录抽取时与上游逐文件比对为**逐字节一致**
（2026-10-02 经 GitHub API 全文比对核实），未从 goji 工作副本带入任何上游没有的改动。

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
