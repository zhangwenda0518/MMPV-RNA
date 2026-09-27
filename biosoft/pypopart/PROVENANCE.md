# pypopart — 来源与版本

## 来源

- 上游仓库：https://github.com/Adamtaranto/pypopart
- 该副本导出日期：2026-09-22
- **上游 commit**：`4562cba8034981a600f2d97afb85c2962cfe81d4`
- commit 日期：2026-09-14
- commit 主题：Merge pull request #52 from Adamtaranto/feature/grid-overlaps-ticks-gui-jobs

本副本**不含 `.git` 目录**（导出时排除），仅供代码查阅与部署参考；
需要 git 操作请在上游克隆（服务器运行环境：`~/MMPV-RNA/biosoft/pypopart`，是完整 git 仓库）。

## 为什么固定在这个 commit

上游 v0.1.0（≤ `97d76b1`, 2025-11-19）存在一个会导致**单倍型网络分析必崩**的 bug：

- `core/haplotype.py: identify_haplotypes_from_alignment` 用 `remove_gaps()` 后的序列**同时**作分组 key
  和单倍型代表序列；
- 而下游 `core/distance.py: hamming_distance` 强制 `len(seq1) == len(seq2)`，否则
  `raise ValueError('Sequences must have same length: ...')`；
- 病毒比对必然含 indel → 各序列去掉 gap 后长度不等 → MJN / MSN / MST / TCS / TSW / parsimony
  等全部网络算法崩溃（都 import 该函数）。

上游 `4562cba` 已修：改为按**含 gap 的完整比对序列**分组（`key = seq.data`、
`sequence_map[key] = seq`），与 PopART 的 `condenseSeqs` 语义一致，保证所有单倍型共享比对长度。

## 与本仓库补丁的关系

`virome_phylo_pipeline/patch_pypopart_haplotype.py`（2026-08-28）是本仓库对此 bug 的**独立修复**，
比上游早 10 天发现。仍保留作为「无法升级时的权宜方案」，但**两者分组语义不同**：

| | 分组 key | gap 位置不同但碱基相同 | 适用 |
|---|---|---|---|
| 上游 `4562cba`（推荐） | 含 gap 序列 | 视为**不同**单倍型（PopART 语义） | 论文级结果，保证与 PopART 可比 |
| 本仓库补丁 | 去 gap 序列 | 归为**同一**单倍型 | 无法升级时的权宜 |

## 验证记录（2026-09-22）

在服务器上做新旧对比验证（同一测试脚本、同一输入）：

**含 indel 比对 → MJN 建网**（5 条序列，比对长 12，gap 数 [0,0,2,1,1]，去 gap 后长度 [12,12,10,11,11]）：

| 版本 | 单倍型代表序列长度 | MJN 建网 |
|---|---|---|
| v0.1.0（`97d76b1`） | `[11, 10, 11, 12]` 不等长 | ✗ `ValueError: Sequences must have same length: 11 vs 10` |
| `4562cba`（本副本） | `[12, 12, 12, 12]` 等长 | ✓ 通过（10 节点） |

**分组语义**（2 条序列，长度均 13，gap 分别在首/尾，去 gap 后碱基序列相同）：

| 版本 | 单倍型数 | 解读 |
|---|---|---|
| v0.1.0 | 1 | 合并（gap 位置差异被忽略） |
| `4562cba` | 2 | 区分（PopART `condenseSeqs` 语义） |

## 安装

上游是纯 Python 包，源码直用或 pip 安装：

```bash
# 方式 1：源码直用（项目现状，build_haplo_outputs.py 用 sys.path 插入 src/）
export PYPOPART_PATH=<此目录>/src

# 方式 2：pip 可编辑安装
pip install -e .
```

依赖：`numpy` / `numba` / `biopython` / `plotly` / `networkx` / `pandas`
（numba 需要 LLVM，Windows 原生环境可能无法运行；实际分析在 Linux 服务器执行）。
