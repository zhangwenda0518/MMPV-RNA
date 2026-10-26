#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""make_novelty_sankey_data.py — 为「已知/新 → DNA/RNA/Viroid → 科 → 属 → 种」桑基图构建五层流量数据。

⚠ 口径常量集中在**本文件**（`TIER_OF_CATEGORY` / `EXPECTED` / `EXPECTED_LT` / 四个默认开关），
  图脚本 `fig_novelty_sankey.py` 不再另写一套 —— 改口径只改这里。

**当前默认口径 = 121**（四开关全开）：291 HQ 植物病毒（去 Unassigned 19、去 Other 87、去长尾 69 后
剩 116）+ 类病毒按物种去重 5 = 121。
  · 第 1 层 = 已知 / 新（「其他」已整体去除）
      已知 = 已知植物病毒(35) + 类病毒物种(5)                     = 40
      新   = 新种候选(21) + 远缘病毒候选 + 远缘(未定)（去长尾后）  = 81
  · 第 2 层 = DNA / RNA / Viroid（Unassigned 已去除）
  · 第 3/4/5 层 = 科 / 属 / 种（16 / 24 / 51 个唯一值，**无长尾聚合块**）
  · 新种命名（**仅 `category == 新种候选` 的 21 条**）= `<属>-goji`；
    同一属有多条时按 contig_id 升序加 `-1`/`-2`/…

**八组合口径**（`--viroid-unit` × `--drop-unassigned` × `--drop-other`）实测期望值见 `EXPECTED`：
190 / 277 / 193 / 296 / **789** / 876 / 792 / 895；
再叠加 `--drop-longtail`（默认开）后的期望值见 `EXPECTED_LT`：
121 / 163 / 122 / 175 / 720 / 762 / 721 / 774。

追加口径 1（用户 2026-10-13）：「类病毒太多了，其实种类没那么多，去冗余下」
  → `--viroid-unit species`（默认）：类病毒每个 species 只算 1 个单位（604 → 5）。
  依据：604 条命中只落在 129 参考 / 3 属 / 5 种（冗余 121×），病毒侧 291 条落在 152 种（冗余 1.9×）。

追加口径 2（用户 2026-10-13）：「去除unass结果」
  → `--drop-unassigned`（默认开）：剔除 `family == "(Unassigned)"` 的记录（19 条病毒），
    即「L2 基因组型未指派」10 条 + 「科未定」9 条；属/种层的 `(Unassigned)` 占位保留。
  ⚠⚠ 类病毒**强制豁免**：`Rubber viroid India/2009` 的科本身就是 `(Unassigned)`，
    不豁免会把 Viroid 从 5 砍到 4，破坏刚定的类病毒去冗余口径（规则 55）。

追加口径 3（用户 2026-10-13）：「去除 L1 的 other」
  → `--drop-other`（默认开）：剔除 tier == 「其他」的全部记录（87 条），L1 只剩 Known / Novel。

追加口径 4（用户 2026-10-13）：「为什么 L1 Novel 的结果，物种层面还是注释给了已知病毒的名字」
  再答「同样的属 goji 1 2 3 这样」
  → `--novel-naming all`（默认）：**全部 Novel 记录的物种层都改名为 `<属>-goji[-N]`**，
    不再显示源表的最近已知物种名。规则与 21 条新种候选一致：
      ① 属内先排「新种候选」（每条 contig 一个名字，按 contig_id 升序）
      ② 再排「远缘」（按源物种名去重，每个物种名一个名字，按 contig_id 首次出现排序）
      ③ 连续编号 1..N；**N == 1 时不加序号**（与旧规则一致）
      ④ `genus == "(Unassigned)"` 时用 `Unassigned-goji[-N]`
    ⚠ 依据：Novel = 新种候选 21 + 远缘病毒候选 75 + 远缘(未定) 51 = 147；
      后 126 条的源表 `Species` 填的是**最近的已知参考物种名**（如 `Potyvirus iristenuis` 23 条），
      图上会让读者误以为「新病毒 = 已知病毒」。改名后 Novel 分支不再出现任何已知物种名。
    `--novel-naming candidate` 可切回旧口径（只有 21 条新种候选改名，其余显示源表名）。
  ⚠ 原始物种名一律保留在 `raw_species` 列，可追溯。

追加口径 5（用户 2026-10-13）：「发生率可以补上嘛」→ 经二选一确认取
  **「L5 物种层 + 多条 contig 取 max」**：
  → 每个单位补 `n_samples` / `prevalence_pct` 两列（分母 **312**，见 `PREVALENCE_DENOM`）；
    图上在 **L5 物种层标签**后附 ` (xx.x%)`，节点含多条 contig 时取 **max**。
  ⚠⚠ 发生率是「**样本维度**」，桑基节点高是「**contig 维度**」→ 发生率**不参与任何几何计算**，
    只作标签附注（否则两种维度混算会让节点高失去「严格 ∝ 单位数」的性质）。
  病毒侧 = `frequency_table.tsv`（本图 185 条 **100% 命中**）；类病毒侧 = `species_info.tsv`（按 species）。

追加口径 6（用户 2026-10-13）：「所有的 other 去掉」→ 经选择确认取
  **「L3 科 / L4 属 / L5 种 的长尾聚合块都去除」**：
  → `--drop-longtail`（默认开）：**逐层迭代剔除** family / genus / species 中单位数 < `MIN_UNITS`(=2)
    的类群，直到收敛（实测 **3 轮**，190 → **121** 单位，实删 69 条 / 69 个类群）。
    ⚠ 必须**逐层迭代**而非只做一轮：科被剔 → 其下所有属/种一起走（否则留下「孤儿的种」）。
    ⚠ 两条豁免（与旧 `protect` 一致）：① 类病毒（按物种去重后每个 species 天然只有 1 条命中）；
      ② 21 条新种候选（`is_new_species`，本图「主角」）→ 剔除后仍剩 **26** 个 1 单位节点。
    ⚠⚠ 消费端**绝不能再聚合**：图脚本据 `manifest["drop_longtail"]` 把阈值置 0，
      否则这 26 个单例又会被聚成 `Other species`（图上判据 = 不出现任何 `Other ...` 标签）。
  → `--no-drop-longtail` 可切回旧口径（190，绘图端按 `THR` 聚合成 3 个 `Other ...` 块）。
  ⚠⚠ `is_new_species` 在**本文件内存里是 bool**（`cid in new_species`），写进 TSV 才是 `"True"` 字符串
    → 本文件内比较必须用真值判断；图脚本读 TSV 才用 `== "True"`。混用会**静默**让豁免集合只剩类病毒。

数据源（本地镜像 reports/_discovery_sci_20261012/data/）：
  · 10_Reports/HQ_plant_viruses_info.tsv         291 行；`Accession` == rescue 的 `contig_id`（已逐条验证 291/291）
  · 09b_Analysis_Verify/virus_validation/rescue_evidence_scored.tsv   291 行；`category` = 8 类中文
  · 10_Reports/Viroid.all_info.tsv               604 行；与 HQ **零重叠**
  · 09b_Analysis_Verify/keep_summary.tsv         101 行；**101 ⊂ 291**，本图不单独使用（仅做自检）
  · 10_Reports/frequency_table.tsv               291 行；`n_samples` / `prevalence_%`（分母 312）→ 发生率（病毒侧）
  · 10_Reports/Viroid.species_info.tsv           5 行；`n_samples` → 发生率（类病毒侧，按 species）
"""
import argparse
import csv
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "reports" / "_discovery_sci_20261012" / "data"
OUT = Path(__file__).resolve().parent / "novelty_sankey_data"

# ── 口径常量（改这里 = 改图，勿在图脚本里另写一套）──
TIER_OF_CATEGORY = {
    "已知植物病毒": "已知",
    "新种候选": "新", "远缘病毒候选": "新", "远缘(未定)": "新",
    "宿主基因污染": "其他", "昆虫病毒/污染": "其他", "真菌病毒/污染": "其他", "未定": "其他",
}
NEW_SPECIES_CATEGORY = "新种候选"
NEW_SPECIES_SUFFIX = "goji"
EXPECTED_NEW_SPECIES = 21

# ── 类病毒计数口径（用户 2026-10-13 追加指令：「类病毒太多了，其实种类没那么多，去冗余下」）──
#   contig  = 逐条 BLAST 命中计（604 条 → 全图 895，类病毒占 67.5 %）
#   species = **按物种去重**，每个物种算 1 个单位（5 个 → 全图 296，类病毒占 1.7 %）
# 病毒侧两种口径下都是 291 条 contig，不参与去重。
# 依据：604 条命中只落在 129 条参考序列 / 3 个属 / 5 个种上（冗余度 121×，
#       而病毒侧 291 条落在 152 个种上，冗余度仅 1.9×）。
VIROID_UNIT_DEFAULT = "species"
NOVEL_NAMING_DEFAULT = "all"
UNASSIGNED = "(Unassigned)"
# (viroid_unit, drop_unassigned, drop_other) -> 期望值（自检用；8 组合全部实测过）
# ⚠ 2026-10-13 修正：`(contig, True, True)` 原记 773 = 876-103-19，把两个过滤的**交集 16 条**
#   扣了两次（`drop_unassigned` 会顺带删掉 16 条 tier==「其他」的记录）。
#   正确 = 895 - |U ∪ O| = 895 - (19 + 103 - 16) = **789**。已实测复现。
EXPECTED = {
    ("contig",  False, False): {"total": 895, "tiers": {"已知": 642, "新": 150, "其他": 103}, "viroid": 604},
    ("contig",  False, True):  {"total": 792, "tiers": {"已知": 642, "新": 150},                 "viroid": 604},
    ("contig",  True,  False): {"total": 876, "tiers": {"已知": 642, "新": 147, "其他": 87},     "viroid": 604},
    ("contig",  True,  True):  {"total": 789, "tiers": {"已知": 642, "新": 147},                 "viroid": 604},
    ("species", False, False): {"total": 296, "tiers": {"已知": 43,  "新": 150, "其他": 103},    "viroid": 5},
    ("species", False, True):  {"total": 193, "tiers": {"已知": 43,  "新": 150},                 "viroid": 5},
    ("species", True,  False): {"total": 277, "tiers": {"已知": 43,  "新": 147, "其他": 87},     "viroid": 5},
    ("species", True,  True):  {"total": 190, "tiers": {"已知": 43,  "新": 147},                 "viroid": 5},
}

# ── 发生率（prevalence）口径（用户 2026-10-13 追加：「发生率可以补上嘛」）──
#   分母 = **312 个样本**；由 `frequency_table.tsv` 的 `n_samples / prevalence_%` 反推
#   （291/291 条吻合，假设 311 / 313 分别只有 225 / 200 条吻合 → 唯一解）。
#   病毒侧 = 逐 contig（`frequency_table.tsv`，本图用到的 185 条 **100% 命中**）；
#   类病毒侧 = 按 species（`Viroid.species_info.tsv` 的 `n_samples`，与「按物种去重」的单位口径天然对齐）。
#   ⚠⚠ 发生率是「**样本维度**」，而桑基节点高度是「**contig 维度**」——两者不可混算。
#      图上只在 **L5 物种层标签**里作为附注展示，**不参与任何几何计算**（节点高仍 ∝ 单位数）。
#   ⚠ 一个 L5 节点含多条 contig 时取 **max**（用户 2026-10-13 选定），含义 =
#      「该物种至少有一个分离株在 X% 的样本中检出」。
PREVALENCE_DENOM = 312
FREQ_TSV = "10_Reports/frequency_table.tsv"
VIROID_SPECIES_TSV = "10_Reports/Viroid.species_info.tsv"

# ── 长尾剔除口径（用户 2026-10-13 追加：「所有的 other 去掉」）──
#   用户选定：**L3 科 / L4 属 / L5 种 的长尾聚合块全部去除**（不是「改个名」）。
#   做法：不再把「单位数不足」的类群聚合成 `Other families/genera/species`，
#         而是**逐层迭代剔除**单位数 < MIN_UNITS 的类群，直到收敛。
#   ⚠ 两条豁免（`protect`）：① 类病毒（按物种去重后每个物种天然只有 1 条命中）；
#      ② 21 条新种候选（`is_new_species == "True"`，是本图的「主角」，必须保留）。
#      这两类即使单位数 = 1 也**不剔除、不聚合**（图上以细线单节点呈现）。
#   实测：190 → **121** 单位（3 轮收敛）；剔除后 family 19 / genus 25 / species 52 个唯一值，
#         三层均无「长尾聚合块」。
#   ⚠ 剔除后**消费端不得再聚合**：图上 THR 必须置 0（否则 26 个被豁免的单例又会被聚成
#     `Other species`，等于白做）。判据 = 图上不再出现任何 `Other ...` 标签。
MIN_UNITS = 2
# (viroid_unit, drop_unassigned, drop_other) -> 剔除长尾后的期望总单位数（自检用；8 组合实测）
EXPECTED_LT = {
    ("contig",  False, False): 774,
    ("contig",  False, True):  721,
    ("contig",  True,  False): 762,
    ("contig",  True,  True):  720,
    ("species", False, False): 175,
    ("species", False, True):  122,
    ("species", True,  False): 163,
    ("species", True,  True):  121,
}


def rd(p: Path):
    with p.open(encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def norm(x, blank=UNASSIGNED):
    x = (x or "").strip()
    return x if x and x.lower() not in ("na", "nan", "-") else blank


def main():
    ap = argparse.ArgumentParser(description="09b 桑基图数据层")
    ap.add_argument("--viroid-unit", choices=["contig", "species"], default=VIROID_UNIT_DEFAULT,
                    help="类病毒计数口径：contig=逐条命中（604）；species=按物种去重（5）。默认 species。")
    ap.add_argument("--drop-unassigned", action=argparse.BooleanOptionalAction, default=True,
                    help="剔除 `family == \"(Unassigned)\"` 的记录（类病毒强制豁免）。默认开；"
                         "用 --no-drop-unassigned 关掉。")
    ap.add_argument("--drop-other", action=argparse.BooleanOptionalAction, default=True,
                    help="剔除 tier == 「其他」的全部记录，L1 只剩 Known / Novel。默认开；"
                         "用 --no-drop-other 关掉。")
    ap.add_argument("--novel-naming", choices=["all", "candidate"], default=NOVEL_NAMING_DEFAULT,
                    help="Novel 物种层命名：all=全部 Novel 改名为 <属>-goji[-N]（默认）；"
                         "candidate=只有 21 条新种候选改名，其余沿用源表物种名。")
    ap.add_argument("--drop-longtail", action=argparse.BooleanOptionalAction, default=True,
                    help="剔除 L3/L4/L5 中单位数 < MIN_UNITS 的类群（即不再聚合成 Other 长尾块）。"
                         "类病毒与新种候选豁免。默认开；用 --no-drop-longtail 关掉（回到聚合口径）。")
    ap.add_argument("--novel-dedup-per-genus", type=int, default=1,
                    help="【用户 2026-10-13】同属 Novel 去冗余：每属「同一物种名(goji-N)只留 1 条"
                         "（取发生率最高）」→ 若去重后仍 > N 条则按发生率砍到 N 条。"
                         "用户最终定 **N=1（每属只展示 1 条代表）**，避免新病毒过多挤爆桑基；"
                         "完整明细由附图（全部 Novel contig 清单）承载。Known 完全不动。"
                         "设为 0 关闭。默认 1。")
    ap.add_argument("--novel-detail-out", default=None,
                    help="把「去冗余前的全部 Novel 明细」另存一份 TSV 作附图数据（默认写入 "
                         "novelty_sankey_data/novel_dedup_full.tsv）。")
    args = ap.parse_args()
    unit = args.viroid_unit
    exp = EXPECTED[(unit, args.drop_unassigned, args.drop_other)]

    hq = {r["Accession"]: r for r in rd(DATA / "10_Reports" / "HQ_plant_viruses_info.tsv")}
    res = {r["contig_id"]: r for r in rd(DATA / "09b_Analysis_Verify" / "virus_validation" /
                                         "rescue_evidence_scored.tsv")}
    vir = rd(DATA / "10_Reports" / "Viroid.all_info.tsv")
    # 发生率源（分母 312）：病毒逐 contig / 类病毒按 species
    freq = {r["contig_id"]: r for r in rd(DATA / FREQ_TSV)}
    vinf = {r["species"]: r for r in rd(DATA / VIROID_SPECIES_TSV)}

    # ── 自检 1：HQ 与 rescue 必须同集（否则 join 会静默丢行）──
    if set(hq) != set(res):
        only_hq, only_res = set(hq) - set(res), set(res) - set(hq)
        sys.exit(f"❌ HQ 与 rescue 的 contig_id 不同集：HQ 独有 {len(only_hq)} / rescue 独有 {len(only_res)}")

    # ── 自检 2：类病毒与病毒必须零重叠 ──
    vids = {r["contig_id"] for r in vir}
    if vids & set(hq):
        sys.exit(f"❌ 类病毒与 HQ 病毒重叠 {len(vids & set(hq))} 条")

    # ── 新种命名：同属多条 → -1/-2/…（按 contig_id 升序，确定性）──
    ns_ids = sorted(cid for cid, r in res.items() if r["category"] == NEW_SPECIES_CATEGORY)
    by_genus = defaultdict(list)
    for cid in ns_ids:
        by_genus[norm(hq[cid]["Genus"], "(无属)")].append(cid)
    new_species = {}
    for genus, cids in by_genus.items():
        for i, cid in enumerate(cids, 1):
            suffix = f"-{i}" if len(cids) > 1 else ""
            new_species[cid] = f"{genus}-{NEW_SPECIES_SUFFIX}{suffix}"

    # ── 逐 contig 打平（`raw_species` = 源表原物种名，改名前先留底，可追溯）──
    rows = []
    for cid, h in hq.items():
        tier = TIER_OF_CATEGORY.get(res[cid]["category"])
        if tier is None:
            sys.exit(f"❌ 未登记 category: {res[cid]['category']!r}（contig {cid}）")
        gt = norm(h["Genome_Type"], "Unassigned")
        if gt not in ("DNA", "RNA"):
            gt = "Unassigned"
        rows.append(dict(
            contig_id=cid, tier=tier, category=res[cid]["category"], genome_type=gt,
            family=norm(h["Family"]), genus=norm(h["Genus"]),
            species=new_species.get(cid, norm(h["Species"])), raw_species=norm(h["Species"]),
            is_new_species=cid in new_species,
            length=int(h["Length"]) if h["Length"].isdigit() else None,
            nt_pident=res[cid]["nt_pident"], nt_qcovs=res[cid]["nt_qcovs"],
        ))
    # ── 类病毒：按口径展开 ──
    if unit == "contig":
        for r in vir:
            rows.append(dict(
                contig_id=r["contig_id"], tier="已知", category="(类病毒)", genome_type="Viroid",
                family=norm(r["family"]), genus=norm(r["genus"]),
                species=norm(r["species"]), raw_species=norm(r["species"]), is_new_species=False,
                length=int(r["ref_len"]) if (r["ref_len"] or "").isdigit() else None,
                nt_pident=r["pident"], nt_qcovs="",
            ))
    else:
        # **按物种去重**：每个 species 只算 1 个单位；family/genus 取该物种的众数归属。
        by_sp = defaultdict(list)
        for r in vir:
            by_sp[norm(r["species"])].append(r)
        for sp in sorted(by_sp):
            grp = by_sp[sp]
            fam = Counter(norm(r["family"]) for r in grp).most_common(1)[0][0]
            gen = Counter(norm(r["genus"]) for r in grp).most_common(1)[0][0]
            pids = [float(r["pident"]) for r in grp if (r["pident"] or "").strip()]
            rows.append(dict(
                contig_id=f"viroid::{sp}", tier="已知", category="(类病毒)", genome_type="Viroid",
                family=fam, genus=gen, species=sp, raw_species=sp, is_new_species=False,
                length=max(int(r["ref_len"]) for r in grp if (r["ref_len"] or "").isdigit()),
                nt_pident=f"{max(pids):.1f}" if pids else "",
                nt_qcovs=f"n_hits={len(grp)}",
            ))

    # ── Unassigned 过滤（用户 2026-10-13 追加：「去除unass结果」→ 取「去 L2 + 科层未指派」）──
    #   规则 = 剔除 `family == "(Unassigned)"` 的记录：L2「基因组型未指派」10 条全在其中，
    #   加上「仅科未定」9 条病毒 = 19 条；属 / 种层的 `(Unassigned)` 占位保留。
    #   ⚠⚠ 类病毒**强制豁免**：`Rubber viroid India/2009` 的科本身就是 `(Unassigned)`，
    #      不豁免会把它删掉 → Viroid 5 → 4，破坏刚定的类病毒去冗余口径（规则 53 精神）。
    n_dropped = 0
    if args.drop_unassigned:
        kept = [r for r in rows
                if not (r["family"] == UNASSIGNED and r["genome_type"] != "Viroid")]
        n_dropped = len(rows) - len(kept)
        rows = kept

    # ── Other 过滤（用户 2026-10-13 追加：「去除 L1 的 other」）──
    #   L1 只剩 Known / Novel；「其他」（宿主基因污染 / 昆虫·真菌病毒污染 / 未定）全部剔除。
    n_dropped_other = 0
    if args.drop_other:
        kept = [r for r in rows if r["tier"] != "其他"]
        n_dropped_other = len(rows) - len(kept)
        rows = kept

    # ── Novel 物种层改名（用户 2026-10-13 追加：「为什么 L1 Novel 的结果，物种层面还是注释给了
    #    已知病毒的名字 → 同样的属 goji 1 2 3 这样」）──
    #   全部 Novel 记录的物种层改为 `<属>-goji[-N]`；Known 分支保持源表真名。
    #   规则：属内先排「新种候选」（每条 contig 一个名，按 contig_id 升序），
    #         再排「远缘」（按 `raw_species` 去重，每个物种名一个名，按 contig_id 首次出现排序）；
    #         连续编号 1..N；**N == 1 时不加序号**（与 21 条新种候选的旧规则一致）。
    #   ⚠ `genus == "(Unassigned)"` 时前缀用 `Unassigned`（避免 `(Unassigned)-goji` 这种怪名）。
    n_renamed = 0
    if args.novel_naming == "all":
        novel = [r for r in rows if r["tier"] == "新"]
        by_g = defaultdict(list)
        for r in novel:
            by_g[r["genus"]].append(r)
        for g, rs in by_g.items():
            prefix = "Unassigned" if g == UNASSIGNED else g
            cand = sorted([r for r in rs if r["category"] == NEW_SPECIES_CATEGORY],
                          key=lambda r: r["contig_id"])
            cand_ids = {r["contig_id"] for r in cand}
            dist_first = {}
            for r in sorted([r for r in rs if r["category"] != NEW_SPECIES_CATEGORY],
                            key=lambda r: r["contig_id"]):
                dist_first.setdefault(r["raw_species"], r)          # 同一物种名只取首次出现
            order = cand + [dist_first[k] for k in dist_first]
            n = len(order)
            name_of_unit = {}
            for i, r0 in enumerate(order, 1):
                nm = f"{prefix}-{NEW_SPECIES_SUFFIX}" + (f"-{i}" if n > 1 else "")
                key = r0["contig_id"] if r0["contig_id"] in cand_ids else ("sp", r0["raw_species"])
                name_of_unit[key] = nm
            for r in rs:
                key = r["contig_id"] if r["category"] == NEW_SPECIES_CATEGORY else ("sp", r["raw_species"])
                r["species"] = name_of_unit[key]
                n_renamed += 1

    # ── 自检 3：层 1 与总数必须与口径一致 ──
    tc = Counter(r["tier"] for r in rows)
    n_vir = sum(1 for r in rows if r["genome_type"] == "Viroid")
    if len(rows) != exp["total"] or dict(tc) != exp["tiers"] or n_vir != exp["viroid"]:
        sys.exit(f"❌ 口径对不上（unit={unit}, drop_unassigned={args.drop_unassigned}, "
                 f"drop_other={args.drop_other}）：总 {len(rows)}（期望 {exp['total']}），"
                 f"分层 {dict(tc)}（期望 {exp['tiers']}），类病毒 {n_vir}（期望 {exp['viroid']}）")
    if len(new_species) != EXPECTED_NEW_SPECIES:
        sys.exit(f"❌ 新种数 {len(new_species)}（期望 {EXPECTED_NEW_SPECIES}）")
    if args.drop_unassigned and any(r["family"] == UNASSIGNED and r["genome_type"] != "Viroid"
                                    for r in rows):
        sys.exit("❌ 过滤后仍存在非类病毒的 `(Unassigned)` 科记录")
    if args.drop_other and any(r["tier"] == "其他" for r in rows):
        sys.exit("❌ 过滤后仍存在 tier == 「其他」的记录")
    # ⚠ Novel 分支不得残留任何「已知物种名」：改名后每条 Novel 的 species 必须以 -goji 结尾
    if args.novel_naming == "all":
        bad = [r for r in rows if r["tier"] == "新"
               and f"-{NEW_SPECIES_SUFFIX}" not in r["species"]]
        if bad:
            sys.exit(f"❌ Novel 分支仍有 {len(bad)} 条未改名（例：{bad[0]['species']!r}）")

    # ── 长尾剔除（用户 2026-10-13 追加：「所有的 other 去掉」）──
    #   ⚠⚠ 2026-10-13 晚更新（用户「只对物种层去长尾」）：**科/属层不再剔除**。
    #      原因：同属 Novel 去冗余（cap 3/属）会把 `Closteroviridae`/`Alphaflexiviridae`
    #      等科压到 1 个单位 → 若仍按「各层 <MIN_UNITS 全剔」会把**真实新病毒整个删掉**
    #      （`Closterovirus-goji` 等），与「展示较完整」矛盾。
    #      → 现在**只在 species 层剔除**；科/属层保留全部（保证每个病毒都能追到完整链）。
    #   ⚠ 剔除仍**逐轮迭代**（species 被剔 → 其所属科/属计数下降，但科/属不再因此被剔）。
    n_dropped_lt = 0
    lt_dropped_log = []
    if args.drop_longtail:
        protect = {("species", r["species"]) for r in rows
                   if r["genome_type"] == "Viroid" or r["is_new_species"]}
        for _round in range(1, 51):
            changed = False
            for lv in ("species",):
                cnt = Counter(r[lv] for r in rows)
                for v, c in sorted(cnt.items()):
                    if c >= MIN_UNITS:
                        continue
                    keep = any(r[lv] == v and ("species", r["species"]) in protect for r in rows)
                    if keep:
                        continue
                    rows = [r for r in rows if r[lv] != v]
                    n_dropped_lt += c
                    lt_dropped_log.append((_round, lv, v, c))
                    changed = True
            if not changed:
                break
        else:
            sys.exit("❌ 长尾剔除 50 轮未收敛（疑似剔除规则自相矛盾）")

    # ── 自检 3b：剔除后物种层不得残留「未被豁免的」长尾类群 ──
    #   ⚠⚠ 2026-10-13 晚：**科/属层不再剔除**（见上）→ 本自检也只看 species 层。
    #     科/属层的单例是**合法保留**（保证科–属–种链完整），不再判为残留。
    if args.drop_longtail:
        if n_dropped_lt == 0:
            sys.exit("❌ 开了 --drop-longtail 但一条都没剔除（口径或输入异常）")
        # ⚠ 期望值：现在只剔 species 层 → EXPECTED_LT 表按「只种层」口径重算。
        #   若 --novel-dedup-per-genus 开启，则总数在去冗余后才定 → 此处跳过精确比对，
        #   改由「去冗余后」的自检（下方）兜底。
        _exp_lt = EXPECTED_LT[(unit, args.drop_unassigned, args.drop_other)]
        if not (args.novel_dedup_per_genus and args.novel_dedup_per_genus > 0):
            if len(rows) != _exp_lt:
                sys.exit(f"❌ 剔除长尾后总 {len(rows)}，期望 {_exp_lt}"
                         f"（unit={unit}, drop_unassigned={args.drop_unassigned}, "
                         f"drop_other={args.drop_other}）")
        prot = {("species", r["species"]) for r in rows
                if r["genome_type"] == "Viroid" or r["is_new_species"]}
        for lv in ("species",):
            bad = {v: c for v, c in Counter(r[lv] for r in rows).items()
                   if c < MIN_UNITS and not any(r[lv] == v and ("species", r["species"]) in prot
                                                for r in rows)}
            if bad:
                sys.exit(f"❌ {lv} 层仍残留未豁免的长尾类群：{list(bad.items())[:5]}")
        # ⚠ 剔除后必须重算分层计数（`tc` 是剔除前的快照；否则 meta.tiers 与 meta.total 打架）
        tc = Counter(r["tier"] for r in rows)
        n_vir = sum(1 for r in rows if r["genome_type"] == "Viroid")

    # ── 发生率 attach（用户 2026-10-13 追加：「发生率可以补上嘛」）──
    #   病毒按 contig_id 查 frequency_table；类病毒按 species 查 species_info。
    n_prev_missing = 0
    for r in rows:
        if r["genome_type"] == "Viroid":
            rec = vinf.get(r["species"])
        else:
            rec = freq.get(r["contig_id"])
        ns = rec["n_samples"].strip() if rec and (rec.get("n_samples") or "").strip() else ""
        r["n_samples"] = int(ns) if ns.isdigit() else ""
        r["prevalence_pct"] = round(r["n_samples"] / PREVALENCE_DENOM * 100, 2) if r["n_samples"] != "" else ""
        if r["n_samples"] == "":
            n_prev_missing += 1

    # ── 自检 8：分母 312 必须能复现 frequency_table 自带的 prevalence_% ──
    bad_denom = [cid for cid, f in freq.items()
                 if (f["n_samples"] or "").strip() and (f["prevalence_%"] or "").strip()
                 and abs(round(int(f["n_samples"]) / PREVALENCE_DENOM * 100, 2)
                         - float(f["prevalence_%"])) > 0.011]
    if bad_denom:
        sys.exit(f"❌ 发生率分母 {PREVALENCE_DENOM} 与 frequency_table 自带的 prevalence_% 不符"
                 f"（{len(bad_denom)}/{len(freq)} 条，例：{bad_denom[0]}）")
    # ── 自检 9：本图用到的单位必须**全部**拿到发生率（静默漏 join 会让图上出现空百分比）──
    if n_prev_missing:
        miss = [r["contig_id"] for r in rows if r["n_samples"] == ""][:3]
        sys.exit(f"❌ {n_prev_missing} 个单位缺发生率（例：{miss}）")

    # ══════════════════════════════════════════════════════════════════════
    # 同属 Novel 去冗余（用户 2026-10-13：「同属 物种去冗余 新病毒 每属保留前3个差异很大的」）
    # ══════════════════════════════════════════════════════════════════════
    #   动机：Novel 81 条把图挤爆（Potyvirus 24 / Betacytorhabdovirus 14）；
    #         但其中大量是「同一新种的多个 contig」（物种名同为 <属>-goji-N）→ 视觉重复。
    #   规则（用户 2026-10-13 最终定案：「同属的新病毒只展示一个」「详细结果放在附图里」）：
    #     ① 同属内按 `species` 分组（= 同一个 goji-N），每组只留 **1** 条代表：
    #        取 `prevalence_pct` 最高者（并列时取 contig_id 字典序最小 → 确定性）。
    #     ② 若属内去重后仍 > N 条（**默认 N=1**：每属只展示 1 条代表），
    #        按 `prevalence_pct` 降序砍到 N 条。
    #   ⭐ **为什么 N=1**：用户「新病毒太多了，同属 物种去冗余」→ 桑基图每属只留 1 条
    #      代表，避免 24 条 Potyvirus / 14 条 Betacytorhabdovirus 挤爆图。
    #      **完整明细另存 `novel_dedup_full.tsv`（附图数据）** → 图简洁 + 信息不丢。
    #   ⚠ **Known 完全不动** —— 用户明确「只 Novel，已知病毒不动」。
    #   ⚠ **命名不变**（用户「只下采样，重新命名不变」）→ 不重跑命名/去长尾链，
    #      本步只**减少展示条数**；因此 L4/L5 节点数可能随之变小（如 Potyvirus 24→1）。
    #   ⚠ 必须有**日志 + 断言**：静默下采样会让「total」这种数字对不上时无人察觉。
    #   ⚠ 本步在 `--drop-longtail` **之后**、发生率 attach **之后** —— 代表按发生率选，
    #      所以顺序不能前移（前移时 `prevalence_pct` 尚未回填）。
    n_novel_dedup_dropped = 0
    novel_dedup_log = []
    # ⭐ 附图数据：**去冗余前**的全部 Novel 明细（供「详细结果」附图消费，信息零丢失）
    _novel_full_rows = [r for r in rows if r["tier"] == "新"]
    if args.novel_dedup_per_genus and args.novel_dedup_per_genus > 0:
        _cap = args.novel_dedup_per_genus
        _known = [r for r in rows if r["tier"] != "新"]
        _novel = [r for r in rows if r["tier"] == "新"]
        _by_genus = defaultdict(list)
        for r in _novel:
            _by_genus[r["genus"]].append(r)
        _kept = []
        for _gen, _rs in sorted(_by_genus.items()):
            # ① 同物种名去重（每 goji-N 留 1 条，取发生率最高）
            _by_sp = defaultdict(list)
            for _r in _rs:
                _by_sp[_r["species"]].append(_r)
            _dedup = []
            for _sp, _cands in _by_sp.items():
                # 取发生率最高；并列时取 contig_id 字典序最小（确定性）
                _best = sorted(_cands, key=lambda r: (-(r["prevalence_pct"] or 0.0),
                                                      r["contig_id"]))[0]
                _dedup.append(_best)
                if len(_cands) > 1:
                    novel_dedup_log.append(("species", _gen, _sp, len(_cands)))
            # ② 仍 > cap → 按发生率砍
            _dedup.sort(key=lambda r: (-(r["prevalence_pct"] or 0.0), r["contig_id"]))
            if len(_dedup) > _cap:
                novel_dedup_log.append(("genus_cap", _gen, "", len(_dedup) - _cap))
            _kept.extend(_dedup[:_cap])
        n_novel_dedup_dropped = len(_novel) - len(_kept)
        # ⚠ 断言：去冗余后 Novel 条数必须**严格减少**（否则旗标形同虚设，静默无效）
        if n_novel_dedup_dropped <= 0:
            sys.exit("❌ 开了 --novel-dedup-per-genus 但一条都没删（口径或输入异常）")
        # ⚠ 断言：保留集必须是原集的子集（防拼接 bug 引入幻影行）
        _orig_ids = {r["contig_id"] for r in _novel}
        _kept_ids = {r["contig_id"] for r in _kept}
        if not _kept_ids <= _orig_ids:
            sys.exit(f"❌ 去冗余保留了不存在的 contig：{sorted(_kept_ids - _orig_ids)[:3]}")
        if len(_kept_ids) != len(_kept):
            sys.exit("❌ 去冗余保留集有重复 contig_id（同一 contig 被两条路径留下）")
        rows = _known + _kept
        # ⚠ 去冗余后必须重算分层计数（`tc` 是剔除前的快照）
        tc = Counter(r["tier"] for r in rows)
        n_vir = sum(1 for r in rows if r["genome_type"] == "Viroid")
        # ⚠ 新种候选属分布也要重算（`by_genus` 是去冗余前的快照，会与 total 打架）
        by_genus = defaultdict(list)
        _ns_kept = sorted(r["contig_id"] for r in rows if r["is_new_species"])
        for _cid in _ns_kept:
            by_genus[norm(hq[_cid]["Genus"], "(无属)")].append(_cid)

        # ⚠⚠ 去冗余后**不再做二次去长尾**（2026-10-13 决策）：
        #   去冗余会让很多属只剩 1 个物种（`Soymovirus`/`Vitivirus`/`Umbravirus`/`Closterovirus`…），
        #   若再按「物种单位数 < 2 剔除」，会把这 **26 个真实新病毒属整个删掉** ——
        #   与用户「展示较完整」的核心诉求直接矛盾（实测：58 vs 125 条）。
        #   → **去冗余本身就是「减少视觉冗余」的手段**，长尾剔除不再叠加。
        #   （第一次去长尾仍保留：它剔的是**去冗余前**的物种单例，即真正的长尾。）

    OUT.mkdir(parents=True, exist_ok=True)

    # ⭐ 附图数据：去冗余前的**全部 Novel 明细**（用户「详细结果放在附图里」）
    #   本表**不参与桑基几何**，只作「完整结果」附图/附表的数据源 → 信息零丢失。
    #   ⚠ 必须与桑基用同一份 `rows`（同一次运行、同一发生率口径），否则两图数字打架。
    _cols_full = ["contig_id", "tier", "category", "genome_type", "family", "genus",
                  "species", "raw_species", "is_new_species", "length", "nt_pident",
                  "nt_qcovs", "n_samples", "prevalence_pct", "shown_in_sankey"]
    _kept_ids_final = {r["contig_id"] for r in rows}
    _full_sorted = sorted(_novel_full_rows,
                          key=lambda r: (r["genus"], r["species"], -(r["prevalence_pct"] or 0.0),
                                         r["contig_id"]))
    _detail_path = (Path(args.novel_detail_out) if args.novel_detail_out
                    else OUT / "novel_dedup_full.tsv")
    _detail_path.parent.mkdir(parents=True, exist_ok=True)
    with _detail_path.open("w", encoding="utf-8", newline="\n") as fh:
        _w = csv.DictWriter(fh, fieldnames=_cols_full, delimiter="\t", extrasaction="ignore")
        _w.writeheader()
        for _r in _full_sorted:
            _row = dict(_r)
            _row["shown_in_sankey"] = "True" if _r["contig_id"] in _kept_ids_final else "False"
            _w.writerow(_row)
    _n_shown = sum(1 for r in _full_sorted if r["contig_id"] in _kept_ids_final)
    print(f"   📄 附图数据 {_detail_path.name}：Novel 全量 {len(_full_sorted)} 条"
          f"（桑基展示 {_n_shown} / 折叠 {len(_full_sorted) - _n_shown}）")

    # ① 逐 contig 明细（可追溯）
    cols = ["contig_id", "tier", "category", "genome_type", "family", "genus",
            "species", "raw_species", "is_new_species", "length", "nt_pident", "nt_qcovs",
            "n_samples", "prevalence_pct"]
    with (OUT / "contigs.tsv").open("w", encoding="utf-8", newline="\n") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader()
        for r in sorted(rows, key=lambda x: (x["tier"], x["genome_type"], x["family"], x["genus"], x["species"])):
            w.writerow(r)

    # ② 五层节点 / 边
    LEVELS = [("tier", "tier"), ("genome_type", "genome_type"),
              ("family", "family"), ("genus", "genus"), ("species", "species")]
    nodes, links = [], []
    for depth, (key, _) in enumerate(LEVELS):
        cnt = Counter(r[key] for r in rows)
        for name, n in cnt.items():
            nodes.append(dict(level=depth, key=key, name=name, n=n))
    for depth in range(len(LEVELS) - 1):
        a, b = LEVELS[depth][0], LEVELS[depth + 1][0]
        cnt = Counter((r[a], r[b]) for r in rows)
        for (sa, sb), n in cnt.items():
            links.append(dict(level=depth, source=sa, target=sb, n=n))
    with (OUT / "nodes.tsv").open("w", encoding="utf-8", newline="\n") as fh:
        w = csv.DictWriter(fh, fieldnames=["level", "key", "name", "n"], delimiter="\t")
        w.writeheader(); w.writerows(nodes)
    with (OUT / "links.tsv").open("w", encoding="utf-8", newline="\n") as fh:
        w = csv.DictWriter(fh, fieldnames=["level", "source", "target", "n"], delimiter="\t")
        w.writeheader(); w.writerows(links)

    # ③ 新种命名清单
    with (OUT / "new_species.tsv").open("w", encoding="utf-8", newline="\n") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["new_species", "genus", "old_species_name", "contig_id", "family",
                    "length", "nt_pident", "nt_qcovs", "genus_n_new"])
        for cid in sorted(new_species, key=lambda c: (new_species[c], c)):
            h = hq[cid]
            w.writerow([new_species[cid], norm(h["Genus"]), h["Species"], cid, norm(h["Family"]),
                        h["Length"], res[cid]["nt_pident"], res[cid]["nt_qcovs"],
                        len(by_genus[norm(h["Genus"])])])

    meta = dict(
        total=len(rows), tiers=dict(tc),
        viroid_unit=unit,
        viroid_units=n_vir,
        viroid_note=("类病毒按物种去重，每个 species 算 1 个单位（604 条 BLAST 命中 → 5 个种）"
                     if unit == "species" else "类病毒逐条 BLAST 命中计数（604 条）"),
        drop_unassigned=bool(args.drop_unassigned),
        n_dropped=n_dropped,
        drop_rule=("剔除 `family == \"(Unassigned)\"` 的非类病毒记录（L2 基因组型未指派 + 仅科未定）；"
                   "类病毒强制豁免，属/种层 (Unassigned) 占位保留"
                   if args.drop_unassigned else "未过滤，保留全部记录"),
        drop_other=bool(args.drop_other),
        n_dropped_other=n_dropped_other,
        drop_other_rule=("剔除 tier == 「其他」的全部记录（宿主基因污染 / 昆虫·真菌病毒污染 / 未定）；"
                         "L1 只剩 Known / Novel"
                         if args.drop_other else "未过滤，保留「其他」tier"),
        drop_longtail=bool(args.drop_longtail),
        min_units=MIN_UNITS,
        # ⚠ 2026-10-13 晚：用户定「只对物种层去长尾」→ 科/属层单例合法保留。
        longtail_scope=(["species"] if args.drop_longtail else []),
        n_dropped_longtail=n_dropped_lt,
        longtail_rule=(f"**仅在 species 层**逐轮剔除单位数 < {MIN_UNITS} 的物种，直到收敛；"
                       f"类病毒与 21 条新种候选豁免。⚠ 科/属层**不剔除**（同属 Novel 去冗余名"
                       f"会把科压到 1 单位，若也剔会删掉真实新病毒）"
                       if args.drop_longtail else
                       "未剔除，L3/L4/L5 的单位数不足类群由绘图端聚合成 `Other ...` 长尾块"),
        longtail_dropped_sample=lt_dropped_log[:12],
        novel_naming=args.novel_naming,
        n_novel_renamed=n_renamed,
        novel_naming_rule=("全部 Novel 记录的物种层改为 `<属>-goji[-N]`（不再显示源表最近已知物种名）；"
                           "属内先排新种候选（每条 contig 一名），再排远缘（按源物种名去重），"
                           "编号 1..N，N==1 不加序号"
                           if args.novel_naming == "all" else
                           "仅 21 条新种候选改名，其余沿用源表物种名"),
        novel_dedup_per_genus=(args.novel_dedup_per_genus
                               if args.novel_dedup_per_genus and args.novel_dedup_per_genus > 0
                               else None),
        n_novel_dedup_dropped=n_novel_dedup_dropped,
        novel_dedup_rule=("同属 Novel 去冗余：① 同一物种名(goji-N)只留 1 条（取发生率最高，"
                          "并列取 contig_id 字典序最小）；② 若去重后仍 > N 条则按发生率砍到 N 条。"
                          "Known 完全不动；命名不变（只下采样）。"
                          if args.novel_dedup_per_genus and args.novel_dedup_per_genus > 0 else
                          "未去冗余，Novel 全部 contig 均展示"),
        novel_dedup_dropped_sample=novel_dedup_log[:12],
        # ⭐ 附图数据（去冗余前全部 Novel 明细）：信息零丢失，供「详细结果」附图消费
        novel_detail_file=_detail_path.name,
        novel_detail_total=len(_full_sorted),
        novel_detail_shown=_n_shown,
        genome_types=dict(Counter(r["genome_type"] for r in rows)),
        n_nodes={LEVELS[d][0]: len({r[LEVELS[d][0]] for r in rows}) for d in range(len(LEVELS))},
        n_links=len(links),
        new_species=sum(1 for r in rows if r["is_new_species"]),
        new_species_genera={g: len(c) for g, c in sorted(by_genus.items())},
        species_singletons=sum(1 for _, n in Counter(r["species"] for r in rows).items() if n == 1),
        naming_rule=f"<Genus>-{NEW_SPECIES_SUFFIX}" + ("-N（同属多条时按 contig_id 升序）"),
        sources=["10_Reports/HQ_plant_viruses_info.tsv",
                 "09b_Analysis_Verify/virus_validation/rescue_evidence_scored.tsv",
                 "10_Reports/Viroid.all_info.tsv"],
        prevalence_denom=PREVALENCE_DENOM,
        prevalence_source={"virus": FREQ_TSV, "viroid": VIROID_SPECIES_TSV},
        prevalence_rule=("发生率 = n_samples / 312（**样本维度**）。病毒按 contig（frequency_table）、"
                         "类病毒按 species（species_info）；L5 节点含多条 contig 时图上取 **max**。"
                         "不参与任何几何计算（节点高仍 ∝ 单位数）"),
        prevalence_missing=n_prev_missing,
        tier_rule=TIER_OF_CATEGORY,
    )
    with (OUT / "manifest.json").open("w", encoding="utf-8") as fh:
        json.dump(meta, fh, ensure_ascii=False, indent=2)

    print(f"✅ {OUT}")
    print(f"   总计 {len(rows)}  分层 {dict(tc)}   类病毒口径 = {unit}（{n_vir} 个单位）")
    print(f"   去 Unassigned = {args.drop_unassigned}（实删 {n_dropped} 条）")
    print(f"   类型 {dict(Counter(r['genome_type'] for r in rows))}")
    print(f"   去 Other = {args.drop_other}（实删 {n_dropped_other} 条）")
    print(f"   Novel 命名 = {args.novel_naming}（改名 {n_renamed} 条）")
    print(f"   去长尾 = {args.drop_longtail}（实删 {n_dropped_lt} 条 / "
          f"{len(lt_dropped_log)} 个类群；MIN_UNITS={MIN_UNITS}）")
    if args.novel_dedup_per_genus and args.novel_dedup_per_genus > 0:
        print(f"   同属 Novel 去冗余 = cap {args.novel_dedup_per_genus}/属"
              f"（实删 {n_novel_dedup_dropped} 条；Known 不动）")
        for _k, _g, _s, _n in novel_dedup_log[:8]:
            if _k == "species":
                print(f"      · 物种去重 {_g} / {_s}：{_n} → 1")
            else:
                print(f"      · 属内超帽 {_g}：再删 {_n} 条")
    print(f"   节点数 {meta['n_nodes']}  边数 {len(links)}")
    print(f"   新种 {meta['new_species']} 条；多条的属: "
          f"{ {g: len(c) for g, c in sorted(by_genus.items()) if len(c) > 1} }")
    print(f"   单条种（长尾）{meta['species_singletons']} 个")
    for lv in ("family", "genus", "species"):
        print(f"   {lv:8s} 层 (Unassigned) 剩 {sum(1 for r in rows if r[lv] == UNASSIGNED)} 条")
    print("   Novel 的物种层（前 12 个，确认无已知物种名）：")
    seen, shown = set(), []
    for r in rows:
        if r["tier"] == "新" and r["species"] not in seen:
            seen.add(r["species"])
            shown.append(f"{r['species']} (n={sum(1 for x in rows if x['species'] == r['species'])})")
    print("     " + " | ".join(sorted(shown)[:12]))
    # ── 发生率摘要（分母 312；L5 节点取 max）──
    from collections import defaultdict as _dd
    _g = _dd(list)
    for r in rows:
        _g[r["species"]].append(r["prevalence_pct"])
    _l5 = {k: max(v) for k, v in _g.items() if v}
    print(f"   发生率（分母 {PREVALENCE_DENOM}）：L5 节点 {len(_l5)} 个，"
          f"0% {sum(1 for v in _l5.values() if v == 0)} 个 / <1% {sum(1 for v in _l5.values() if v < 1)} 个 / "
          f">=50% {sum(1 for v in _l5.values() if v >= 50)} 个")
    _top = sorted(_l5.items(), key=lambda kv: -kv[1])[:6]
    print("   L5 发生率 Top6: " + " | ".join(f"{k}={v:.2f}%" for k, v in _top))


if __name__ == "__main__":
    main()
