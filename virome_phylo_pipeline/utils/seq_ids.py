"""seq_ids.py — 序列名/样品名 → 规范键 (单一口径)

## 为什么需要

流水线里"从序列名解析出样品名"这件事曾在至少三处各写一遍, 且规则不同:

| 位置 | 旧规则 | 对 `CRR1126135.OR489165.1` 的结果 |
|---|---|---|
| `utils/data_collector.py` | `split('.')[0]` | `CRR1126135` ✅ |
| `build_haplo_outputs.py` | `split('_')[0]` | `CRR1126135.OR489165.1` ❌ |
| `popgen_analysis._norm_seq_key` | 有 `/` 取 `split('/')[1]`, 否则 `split('_')[0]` | `CRR1126135.OR489165.1` ❌ |

后果是**静默的**: 元数据 CSV 的 `name` 列与比对里的序列名只要后缀形式不同,
匹配就落空 → 地点/宿主分组全变 `Unknown`、单倍型网络的 population 全丢、
Fst 分组为空 —— 但流程不报错, 只是结果悄悄变错。

## 现在的规则 (顺序很重要)

1. 去掉集合式前缀: `GB/CRR1126135` → `CRR1126135` (取最后一段)
2. 去掉 accession 版本后缀: `CRR1126135.OR489165.1` → `CRR1126135`
3. 去掉 `_` 后缀: `CRR1126135_2` → `CRR1126135`

与旧行为在**常见输入上完全一致**, 但两侧用同一个函数后, 混合写法也能对上。
"""
from __future__ import annotations


def base_sample_id(seq_id) -> str:
    """序列名 → 规范样品名。空/None 原样返回。"""
    if seq_id is None:
        return ""
    s = str(seq_id).strip()
    if not s:
        return s
    if "/" in s:              # ① 集合式前缀
        s = s.split("/")[-1]
    s = s.split(".")[0]       # ② accession 版本号
    s = s.split("_")[0]       # ③ 下划线后缀
    return s
