# fill 层根因报告：分类错误的上游源头（fill-fix1，2026-09-16）

数据范围：`~/MMPV-paper/goji-virome/02_novel_virus/RNA-*_out` 7 个 + `~/MMPV-paper/onekp-virome/onekp-virus`，共 8 数据集，成品 `<ds>/05_Taxonomy/Votus.integrated/final_integrated_classification.tsv`。

## 0. 一句话

`fill_taxonomy_na` 旧版不是「补 NA」，而是「用行内最深阶元名去 NCBI 查谱系、覆盖整行阶元」；名字→谱系映射按最后写入者胜，撞上 NCBI 同名歧义（`environmental samples` 有 7748 条同名行 / 47 个互不相干的病毒谱系）就把整块阶元写成无关谱系；`run_vc_analysis` 里又重复调用一次，第二遍把第一遍刚写入的通用名再当键查一遍，劫持被二次放大。这是「Phycodnaviridae / Mimiviridae 变成 Adenoviridae」的直接来源。

## 1. 代码证据

文件：`~/MMPV-RNA/virome_discovery_pipeline/virus_classifier.py`

| 版本 | 大小 | mtime | 说明 |
|---|---|---|---|
| `virus_classifier.py.bak_fillname_20260916` | 45923 B | 2026-07-25 12:07 | 修复前 |
| `virus_classifier.py` | 48433 B | 2026-09-16 05:39 | 修复后 |

调用点：

- L471 `fill_taxonomy_na(combined, combined)`：在 `merge_taxonomy_results` 内部，in-place 回填（保留）。
- 旧版 L826-830 `fill_taxonomy_na(merged, merged)`：在 `run_vc_analysis` 内部**又调一次**（已删除，位置留注释）。

旧实现的三处缺陷（`diff -u` 摘录）：

```python
# ① 映射构建：同名最后写入者胜
if any(r[rn] != "NA" for rn in RANK_NAMES):
    n2r[name.lower()] = r
    n2r[name] = r            # 原始大小写也存一份

# ② 行处理：拿行内最深非 NA 名去查谱系，然后覆盖整行
ictv = n2r.get(name.lower()) or n2r.get(name)
for rn in rank_cols:
    rv = ictv.get(rn, "NA")
    if rv != "NA":
        old = ps[idx].strip().strip('"')
        if old != rv and old != "":     # 只跳过空串
            ps[idx] = rv; corrected += 1   # 语义是「纠正/覆盖」，不是「回填 NA」
```

机制串联：

1. 行内最深非 NA 名恰好是通用串（`environmental samples`、`Bracoviriform`、`Polydnaviriformidae` 等）时，查表命中 7748 条同名行中的任意一条无关病毒谱系，把 Realm..Species 整块改成该谱系。
2. 旧版第二遍回填把第一遍刚写进去的通用串当作新的「最深名」再查一次，谱系继续漂移。onekp `Organic Lake phycodnavirus`（Genus 列已是物种级名）的 Family 被写成 `Adenoviridae` 就发生在这里。
3. 覆盖语义意味着即使工具原始输出正确（例如 metabuli 报的 `Floreoviria|Shotokuvirae|Cossaviricota|Quintoviricetes|Piccovirales|Parvoviridae`），也会被查表结果顶掉。

时间线判据：备份 mtime 2026-07-25 + 修复脚本 mtime 2026-09-16 05:39 + 8 数据集 R 步 05:44–06:17（R 脚本 md5 未变，无其他改动方）→ 归因确定。

## 2. 修复口径（fill-fix1，4 条）

1. 映射结构改为 `name -> {tuple(各阶元值): 出现次数}`，同名多谱系并存，不再互相覆盖。
2. 候选谱系必须与行内已有非 NA 阶元**零冲突**（`conf == 0`），有冲突直接跳过。
3. 在零冲突候选里取一致度（agreement）最高且**唯一**者（`len(top) == 1`）；并列即跳过，不猜、不按顺序兜底。
4. 只写 `NA`/`""`/`-` 位置，行内已有非 NA 值一律不动 → fill 变成幂等语义（重复调用不再改变结果），并删除重复调用。

新增自检输出：`[fill] 回填 N 个 NA / 命中映射 N 行 / 自洽性不足跳过 N 行`。

## 3. 实证收益

- family 级参照一致性闸门 `single_ref_conflict`：修前 8 数据集合计 **314**（onekp 300 / barbarum 8 / chinense 3 / ruthenicum 2 / Fusarium 1）→ 修后 **全 0**。
- 定点验证：onekp `Organic Lake phycodnavirus` 行，旧 Family `Adenoviridae`（300/300）→ 新 `Phycodnaviridae`（现该 Genus 下 674 行 Family=Phycodnaviridae），与参照表一致。
- 成品行数合计 588,953 → 589,441（+488）：旧版有的行因整块被写成同一谱系而与相邻行重复合并，修复后还原为独立行。

## 4. 残余（与 fill 缺陷无关，属下游规则/参照库/门槛）

本轮 878 个「值→NA」丢失格，按机制分五桶（`explain_lost_cells.py`）：

| 桶 | 含义 | goji | onekp | 合计 |
|---|---|---|---|---|
| A | 祖先链不相容（值在参照库里的祖先与行内粗阶元冲突）→ 由逐级相容性规则清空，设计使然 | 172 | 487 | 659 |
| A2 | 参照库对该名没有祖先记录 | 7 | 152 | 159 |
| B | 占位/伪名串票数最高（或并列取胜）后不回落 | 20 | 2 | 22 |
| C | 逐级淘汰（门槛 `share >= 0.5` 的直接后果） | 15 | 5 | 20 |
| D | 未解释 → 手工分解后全部归因（见下） | 1 | 17 | 18 |

D 桶分解：chinense 1 格为 Realm 白名单（`KNOWN_REALMS` 清 Floreoviria，见 §5）；onekp 17 格 = 6 格 Realm 白名单 + 4 格票面为空（旧值是 fill 构造）+ 7 格「真种 vs `environmental samples` 1:1 平票，占位串取胜后清空」（与 B 同机制，只是票首是 NA 不是占位串）。

独立第二轴（`probe_lost_backing.py`，判「旧值是否真有工具投票」）：

- goji 215 格：**213 格**旧值在新 combined 里仍被工具投过票（真丢失，由规则/门槛清空）；2 格为 fill 构造值（1 格新票面为空、1 格仅旧票面有值）。
- onekp 663 格：**640 格**真丢失；17 格新票面为空（该阶元没有任何工具投票，NA 是必然）；6 格仅旧票面有值（旧值是 fill 写入的伪造值）。
- 结论：fill 修复没有丢真数据，878 格里 25 格（2.8%）的「丢失」是撤销旧 fill 伪造值。

## 5. 顺带发现的独立缺陷：`KNOWN_REALMS` 白名单静默清空（12,797 行）

R 脚本 L33：

```r
KNOWN_REALMS <- c("Riboviria","Monodnaviria","Duplodnaviria","Varidnaviria","Adnaviria","Ribozyviria")
```

L215：`dt[!tolower(Realm) %in% tolower(KNOWN_REALMS), Realm := NA_character_]`。

参照库 `~/database/taxonomy/rankedlineage.dmp`（mtime 2026-05-23，375.9 MB）实际 realm 取值共 11 种，病毒侧为：`Riboviria 208937`、`Duplodnaviria 38070`、`Floreoviria 14881`、`Varidnaviria 3822`、`Ribozyviria 192`、`Adnaviria 157`、`Singelaviria 62`。**Floreoviria 与 Singelaviria 都在白名单外，而 Monodnaviria 已不在参照库里**（仍出现在 combined，来自工具自带库）。

后果：成品里 Realm=NA 而 combined 有越界 realm 的行共 **12,797**（Alternaria 13 / Aphis 27 / Fusarium 27 / amarum 14 / barbarum 391 / chinense 184 / ruthenicum 217 / onekp 11,924）。数据里实际出现的越界值只有 `Floreoviria`（onekp combined realm 清单：Varidnaviria 759436 / Duplodnaviria 392344 / Riboviria 168639 / Floreoviria 14576 / Adnaviria 5939 / Monodnaviria 436 / Ribozyviria 307）。

文献核实：JVI 2026 `10.1128/jvi.01019-26`（PubMed 42584058）载 2025 释放以 Efunaviria / Floreoviria / Pleomoviria / Volvereviria 四 realm 替代 Monodnaviria；Wikipedia Floreoviria 条目（2025 release，3 phyla / 18 orders / 209 families）一致。

修复建议见汇总文档待决项 ①（改 R 一处常量 + 统一重跑 R）。

## 6. 复核命令

```bash
# 根因 diff
cd ~/MMPV-RNA/virome_discovery_pipeline
diff -u virus_classifier.py.bak_fillname_20260916 virus_classifier.py | head -120

# 丢失格机制归因（需 --comb/--comb-old/--prod-old/--prod-new）
python3 /tmp/explain_lost_cells.py --label <ds> --comb <new combined> --prod-old <bak product> --prod-new <new product>

# 丢失格投票背书
python3 /tmp/probe_lost_backing.py --prod-old ... --prod-new ... --comb-old ... --comb-new ... --label <ds>

# realm 白名单影响面
python3 /tmp/probe_realm_whitelist.py ...
python3 /tmp/probe_realm_inventory.py   # 参照库 realm 全集
```

探针脚本本地副本：`D:\桌面\延伸基因组\MMPV-RNA\scripts\audit\fix_fill\`；服务器 `/tmp/` 同名。
