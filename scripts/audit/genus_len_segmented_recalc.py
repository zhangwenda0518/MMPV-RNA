# -*- coding: utf-8 -*-
"""属平均长度(genus_lens)重算: 补充分段病毒后验证 ±15%(>=85%) 规则.

背景:
    分支 D 拯救规则 = contig 长度 >= 属平均长度 * 0.85 (允许短 15%).
    database/genus_lens 与 VIGA db/genus_len 的属平均长度疑似按 accession
    逐条平均 -- 分段病毒(分体型, 如 begomovirus DNA-A/B, 分体 RNA 病毒)
    的单个 segment 远短于全基因组, 会把属平均长度严重拉低.

方法:
    1. 读 database/final.cluster.ref_info.tsv, 只取 Nuc_Completeness == complete.
    2. NonSegmented_*: 每个 accession 就是一个基因组 (length 原值).
       Segmented_*: 按 Taxid 聚合, 同 Taxid 内按 Segment 标签去重(取最长),
       各 segment 求和 -> 单基因组总长.
    3. 属 = VMR_Genus (缺则 Species_ICTV 首词), 与 dna_genus_length_analysis 一致.
    4. 参考长度规则 (用户指定): 节段(分段)病毒为主的属 (segmented>=50%)
       取该属**最短完整基因组**作为保留基准 -- 段组完整以"段数=属内众数"界定
       (排除残缺 segment 集), 属内单/双组分双峰时以短的完整型为准; 其余属用均值.
    5. 对比旧值 (database/genus_lens, VIGA db/genus_len), 判定旧口径.
    6. 15% 规则验证: 每属计算参考基因组落在 参考长度*0.85 以上的比例
       (旧值 / 新均值 / 新最短口径), 以及 双侧 ±15% 命中率.

输出: reports/genome_length_analysis/genus_len_recalc/
    genus_lens_recalc.tsv      属级新旧对照 + 统计
    validation_15pct.tsv       每属新旧口径的 >=85% / ±15% 命中率
    genus_lens_recalc_dropin.tsv  genus\ttotal (g__ 前缀, 可直接替换 database/genus_lens,
                               分段属=最短完整基因组, 其余=均值)
    REPORT.md
"""
from pathlib import Path

import numpy as np
import pandas as pd

pd.set_option("display.width", 200)

ROOT = Path(r"D:\桌面\延伸基因组\MMPV-RNA")
OUT = ROOT / "reports" / "genome_length_analysis" / "genus_len_recalc"
OUT.mkdir(parents=True, exist_ok=True)

GENUS_TOL = 0.85          # 分支 D: >= 85% 属平均长度即拯救 (允许短 15%)
TWO_SIDED = 0.15          # 双侧 ±15% 命中率
MIN_N = 3                 # 属内基因组数 < MIN_N 不参与规则验证(单条均值无意义)

print(f"pandas {pd.__version__}, numpy {np.__version__}")

# ------------------------------------------------------------------ load
ref = pd.read_csv(ROOT / "database" / "final.cluster.ref_info.tsv", sep="\t", dtype=str)
ref["Length"] = ref["Length"].astype(float)
ref["Segment"] = ref["Segment"].fillna("").str.strip()
print(f"ref_info 记录: {len(ref)}")

comp = ref[ref["Nuc_Completeness"] == "complete"].drop_duplicates(subset=["Accession"]).copy()
print(f"complete 记录: {len(comp)} "
      f"(Segmented_Complete {len(comp[comp['Category']=='Segmented_Complete'])}, "
      f"NonSegmented_Complete {len(comp[comp['Category']=='NonSegmented_Complete'])})")

comp["genus"] = comp["VMR_Genus"].fillna(
    comp["Species_ICTV"].astype(str).str.split().str[0].replace({"nan": None})
)

# ------------------------------------------------------- isolate grouping
# Segmented_* 按 Taxid 聚合; 同 Taxid 同 Segment 标签多条时保留最长 (冗余/多 isolate 兜底)
seg = comp[comp["Category"] == "Segmented_Complete"].copy()
non = comp[comp["Category"] == "NonSegmented_Complete"].copy()

seg["seg_key"] = seg["Segment"].str.replace(r"\s+", "", regex=True).str.upper()
dup_mask = seg.duplicated(subset=["Taxid", "seg_key"], keep=False)
n_dup_tax = seg.loc[dup_mask, "Taxid"].nunique()
seg = seg.sort_values("Length").drop_duplicates(subset=["Taxid", "seg_key"], keep="last")
iso_seg = (seg.groupby("Taxid")
              .agg(isolate_len=("Length", "sum"),
                   n_segments=("Length", "size"),
                   genus=("genus", "first"),
                   species=("Species_NCBI", "first"))
              .reset_index())
print(f"分段病毒: {len(comp[comp['Category']=='Segmented_Complete'])} accessions -> "
      f"{len(iso_seg)} 个 isolate (Taxid); 覆盖 {n_dup_tax} 个 Taxid 存在同段多条记录已去重")

iso_non = (non.rename(columns={"Length": "isolate_len"})
              .assign(n_segments=1, species=non["Species_NCBI"])
              [["Taxid", "isolate_len", "n_segments", "genus", "species"]])
iso = pd.concat([iso_seg, iso_non], ignore_index=True)
iso["is_segmented"] = ~iso["Taxid"].isin(iso_non["Taxid"])
# 完整段组启发式: 属内期望段数 = n_segments 众数; 段数=众数的 isolate 视为段组完整
mode_cnt = iso.groupby("genus")["n_segments"].agg(lambda x: x.mode().iat[0]).rename("mode_cnt")
iso = iso.merge(mode_cnt, on="genus")
print(f"完整基因组 isolate 合计: {len(iso)}")

# -------------------------------------------------- genus stats (new vs old)
g = (iso.groupby("genus")
        .agg(n_isolates=("isolate_len", "size"),
             n_segmented=("is_segmented", "sum"),
             new_mean=("isolate_len", "mean"),
             new_median=("isolate_len", "median"),
             new_min=("isolate_len", "min"),
             new_max=("isolate_len", "max"),
             cv=("isolate_len", lambda x: x.std() / x.mean() if len(x) > 1 else 0.0))
        .reset_index())
comp_pool = iso[iso["n_segments"] == iso["mode_cnt"]]
g = g.merge(comp_pool.groupby("genus")["isolate_len"].min().rename("min_complete_raw"),
            on="genus", how="left")

# min 护栏: 完整段组池内先剔除 < 0.6×池中位数的记录 (betasatellite/DNA-beta ~1.35kb
# 常挂在 begomovirus 等属名下, 会把 min 拖到假值), 再取最短
def _min_guarded(x: pd.Series) -> float:
    med = x.median()
    low = x[x >= 0.6 * med]
    return float(low.min()) if len(low) else float(x.min())

g = g.merge(comp_pool.groupby("genus")["isolate_len"].agg(_min_guarded).rename("min_complete"),
            on="genus", how="left")

old = pd.read_csv(ROOT / "database" / "genus_lens", sep="\t", dtype=str)
old.columns = ["genus", "old_len"]
old["old_len"] = old["old_len"].astype(float)
old["plain"] = old["genus"].str.replace(r"^[gG]__", "", regex=True)

# 用户规则: 分段为主属 (segmented>=50%) 参考 = 最短**完整段组**基因组
# (只取段数=属众数的 isolate 再取 min, 排除残缺 segment 集与离谱小记录); 其余属 = 均值
g["seg_frac"] = g["n_segmented"] / g["n_isolates"]
g["ref_final"] = np.where(g["seg_frac"] >= 0.5, g["min_complete"], g["new_mean"])
g["is_segmented"] = np.where(g["seg_frac"] >= 0.5, "yes", "no")
g["min_over_median"] = g["min_complete"] / g["new_median"]
suspect = g[(g["seg_frac"] >= 0.5) & (g["min_over_median"] < 0.5)]  # 完整段组内最短仍可疑偏小

# 节段属接受范围 [ref_min, ref_max]: 完整单元 = 完整 segment (护栏后) ∪ 完整基因组总长
seg_units = (seg.groupby("genus")["Length"]
                .apply(lambda s: s[s >= 0.6 * s.median()])
                .reset_index(name="L"))
seg_units = (seg_units.groupby("genus")["L"].agg(seg_min="min", seg_max="max")
             .reset_index())
iso_max = iso.groupby("genus")["isolate_len"].max().rename("iso_max").reset_index()
iso_med = iso.groupby("genus")["isolate_len"].median().rename("iso_med").reset_index()
units = (g[["genus", "min_complete"]].rename(columns={"min_complete": "pool_min"})
         .merge(seg_units, on="genus", how="left")
         .merge(iso_med, on="genus", how="left")
         .merge(iso_max, on="genus", how="left"))
# 上界封顶: 异常 taxid (如挂 satellite 的 4 段求和 10.9kb begomovirus) 不可信,
# 总长 max 封顶在 中位×2.5; 双组分 begomovirus (~1.9×中位) 仍在界内
units["iso_max_cap"] = np.minimum(units["iso_max"], units["iso_med"] * 2.5)
units["ref_min"] = units[["seg_min", "pool_min"]].min(axis=1)
units["ref_max"] = units[["seg_max", "iso_max_cap"]].max(axis=1)
g = g.merge(units[["genus", "ref_min", "ref_max"]], on="genus", how="left")

viga = pd.read_csv(r"D:\桌面\延伸基因组\VIGA\db\genus_len", sep="\t", dtype=str)
viga.columns = ["genus", "viga_len"]
viga["viga_len"] = viga["viga_len"].astype(float)
viga["plain"] = viga["genus"].str.replace(r"^[gG]__", "", regex=True)
viga = viga[viga["plain"] != ""]

# 旧口径判定: 旧值接近 "逐 accession 平均" 还是 "逐 isolate 求和平均"?
acc_mean = (comp.groupby("genus")["Length"].mean().rename("acc_mean").reset_index())
# 不完整记录拖累: 全部记录(complete+partial)逐 accession 均值 vs 仅 complete
ref["genus"] = ref["VMR_Genus"].fillna(
    ref["Species_ICTV"].astype(str).str.split().str[0].replace({"nan": None})
)
ref_len_ok = ref.dropna(subset=["genus", "Length"])
acc_mean_all = (ref_len_ok.groupby("genus")["Length"].mean().rename("acc_mean_all").reset_index())
n_partial = int((ref_len_ok["Nuc_Completeness"] != "complete").sum())
partial_only = sorted(set(ref_len_ok["genus"]) - set(iso["genus"]))
old_plain_set = set(old["plain"])
viga_plain_set = set(viga["plain"])

stat = (g.merge(old[["plain", "old_len"]], left_on="genus", right_on="plain", how="left")
         .merge(viga[["plain", "viga_len"]], on="plain", how="left")
         .merge(acc_mean, on="genus", how="left")
         .merge(acc_mean_all, on="genus", how="left"))
stat["ratio_old_mean"] = stat["old_len"] / stat["new_mean"]
stat["ratio_acc_mean"] = stat["acc_mean"] / stat["new_mean"]
stat["drag_partial"] = stat["acc_mean_all"] / stat["acc_mean"]
both = stat.dropna(subset=["old_len"])
print(f"\n与旧 database/genus_lens 可对照的属: {len(both)} (旧表 {len(old)} 属, ref_info 覆盖 {both.notna().shape[0]})")
for col in ["ratio_old_mean", "ratio_acc_mean"]:
    r = both[col].dropna()
    print(f"  {col}: 中位 {r.median():.3f}, 10-90 分位 "
          f"{r.quantile(0.1):.3f}-{r.quantile(0.9):.3f}, 与1偏差>10%的属 {int(((r-1).abs()>0.1).sum())}")
seg_heavy = both[both["n_segmented"] / both["n_isolates"] >= 0.5]
print(f"  分段为主的属 (segmented>=50%): {len(seg_heavy)}, 其 ratio_old_mean 中位 "
      f"{seg_heavy['ratio_old_mean'].median():.3f}  <-- 旧值系统性偏低即为 per-accession 口径")

# ---- 不完整记录适应情况 ----
dp = stat["drag_partial"].dropna()
print(f"\n不完整记录 (Nuc_Completeness != complete): {n_partial} 条 / {len(ref_len_ok)} 条; "
      f"只完整记录缺失的属 {len(partial_only)} 个 (其中旧表有值 {sum(g_ in old_plain_set for g_ in partial_only)}, "
      f"VIGA 有值 {sum(g_ in viga_plain_set for g_ in partial_only)})")
print(f"  若把 partial 混回逐 accession 均值: 属均值再被拖低 {dp.median():.3f} (10-90 分位 "
      f"{dp.quantile(0.1):.3f}-{dp.quantile(0.9):.3f})  <-- 本重算已全部剔除 partial, 未混入")
print(f"  分段集不全适应: 段数=属众数才进 min 池; min 护栏 0.6×池中位")

# ---- 三方表格交叉对比 (old vs VIGA vs 重算) ----
c3 = stat.dropna(subset=["old_len", "viga_len", "new_mean"])
print(f"\n三方均有值的属: {len(c3)}; 两两覆盖: old∩VIGA "
      f"{len(stat.dropna(subset=['old_len','viga_len']))}, old∩ref {len(both)}, VIGA∩ref "
      f"{len(stat.dropna(subset=['viga_len','new_mean']))}")
rv = (c3["old_len"] / c3["viga_len"]).dropna()
print(f"  old/VIGA: 中位 {rv.median():.3f}, 偏差>10% 的属 {int(((rv-1).abs()>0.1).sum())}/{len(rv)}, "
      f">2倍的属 {int((rv>2).sum())}, <0.5的属 {int((rv<0.5).sum())}")
for tag, sub_ in [("非分段属", c3[c3["n_segmented"]/c3["n_isolates"] < 0.5]),
                  ("分段为主属", c3[c3["n_segmented"]/c3["n_isolates"] >= 0.5])]:
    r1 = (sub_["old_len"]/sub_["new_mean"]).dropna()
    r2 = (sub_["viga_len"]/sub_["new_mean"]).dropna()
    print(f"  {tag} ({len(sub_)}): old/new 中位 {r1.median():.3f}, VIGA/new 中位 {r2.median():.3f}")
top_div = c3.assign(rv=abs(c3["old_len"]/c3["viga_len"]-1)).sort_values("rv", ascending=False).head(10)
print("  old 与 VIGA 分歧最大的属:")
print(top_div[["genus", "old_len", "viga_len", "new_mean", "n_segmented"]].to_string(index=False, float_format=lambda x: f"{x:.0f}"))

# ------------------------------------------------- 15% rule validation
def pass_rate(rows, ref_len):
    if ref_len is None or pd.isna(ref_len) or ref_len <= 0:
        return np.nan
    return float((rows >= ref_len * GENUS_TOL).mean()) * 100

def hit_rate(rows, ref_len):
    if ref_len is None or pd.isna(ref_len) or ref_len <= 0:
        return np.nan
    m = float(((rows >= ref_len * (1 - TWO_SIDED)) & (rows <= ref_len * (1 + TWO_SIDED))).mean())
    return m * 100

val = []
for _, row in stat[stat["n_isolates"] >= MIN_N].iterrows():
    sub = iso.loc[iso["genus"] == row["genus"], "isolate_len"]
    val.append({
        "genus": row["genus"], "n_isolates": int(row["n_isolates"]),
        "n_segmented": int(row["n_segmented"]),
        "old_len": row["old_len"], "new_mean": row["new_mean"],
        "ref_final": row["ref_final"],
        "ge85_old_pct": pass_rate(sub, row["old_len"]),
        "ge85_new_pct": pass_rate(sub, row["new_mean"]),
        "ge85_final_pct": pass_rate(sub, row["ref_final"]),
        "within15_old_pct": hit_rate(sub, row["old_len"]),
        "within15_new_pct": hit_rate(sub, row["new_mean"]),
        "within15_final_pct": hit_rate(sub, row["ref_final"]),
    })
val = pd.DataFrame(val)

print(f"\n规则验证 (n>={MIN_N} 的 {len(val)} 属):")
for tag, c_old, c_new in [("旧值", "ge85_old_pct", "within15_old_pct"),
                          ("新均值", "ge85_new_pct", "within15_new_pct"),
                          ("新规则(分段属取最短)", "ge85_final_pct", "within15_final_pct")]:
    s = val[c_old].dropna(); s2 = val[c_new].dropna()
    print(f"  {tag}: >=85% 中位 {s.median():.1f}%, ±15% 中位 {s2.median():.1f}%")

seg_val = val[val["n_segmented"] / val["n_isolates"] >= 0.5]
print(f"  分段为主属 {len(seg_val)} 个: ge85 旧中位 {seg_val['ge85_old_pct'].median():.1f}% "
      f"-> 新均值 {seg_val['ge85_new_pct'].median():.1f}% -> 新规则(最短) {seg_val['ge85_final_pct'].median():.1f}%; "
      f"±15% 旧 {seg_val['within15_old_pct'].median():.1f}% -> 均值 {seg_val['within15_new_pct'].median():.1f}% -> 最短 {seg_val['within15_final_pct'].median():.1f}%")
print(f"  最短值可疑属 (min/median<0.5, 疑似残缺集拉低): {len(suspect)} 个")

# ------------------------------------------------- architecture / outliers
seg_dist_rows = []
for genus_, sub in seg.groupby("genus"):
    cnt = sub.groupby("Taxid")["Length"].size()
    if len(cnt) >= 10:
        seg_dist_rows.append({"genus": genus_, "n_taxid": len(cnt),
                              "n1": int((cnt == 1).sum()), "n2": int((cnt == 2).sum()),
                              "n3p": int((cnt >= 3).sum())})
seg_dist = pd.DataFrame(seg_dist_rows).sort_values("n_taxid", ascending=False)

# 旧表离谱值: 与新 isolate 均值差 3 倍以上 (污染/求和错误)
outliers = both[(both["ratio_old_mean"] > 3) | (both["ratio_old_mean"] < 0.25)].copy()
outliers = outliers.sort_values("ratio_old_mean", ascending=False)

# ------------------------------------------------------------------ output
stat_out = stat[["genus", "n_isolates", "n_segmented", "is_segmented", "ref_final", "ref_min", "ref_max",
                 "new_mean", "new_min", "min_complete_raw",
                 "min_complete", "new_median", "min_over_median", "new_max", "cv",
                 "old_len", "ratio_old_mean", "acc_mean", "ratio_acc_mean", "viga_len",
                 "acc_mean_all", "drag_partial"]].copy()
stat_out = stat_out.sort_values("n_isolates", ascending=False)
stat_out.to_csv(OUT / "genus_lens_recalc.tsv", sep="\t", index=False, float_format="%.1f")
val.to_csv(OUT / "validation_15pct.tsv", sep="\t", index=False, float_format="%.2f")

# v2 表: 下游分支 D/B 节段范围判定用; 旧 loader 只读前 2 列, 天然向后兼容
# (total: 非节段=均值, 节段=最短完整段组; 节段属另给 [ref_min, ref_max])
v2 = stat_out[["genus", "ref_final", "is_segmented", "ref_min", "ref_max",
               "n_isolates", "n_segmented"]].rename(columns={"ref_final": "total"}).copy()
v2["genus"] = "g__" + v2["genus"].astype(str)
v2 = v2.dropna(subset=["total"])
v2.to_csv(ROOT / "database" / "genus_lens_v2.tsv", sep="\t", index=False, float_format="%.1f")
v2.to_csv(OUT / "genus_lens_v2.tsv", sep="\t", index=False, float_format="%.1f")

dropin = stat_out.dropna(subset=["ref_final"]).copy()
dropin["genus"] = "g__" + dropin["genus"].astype(str)
dropin[["genus", "ref_final"]].rename(columns={"ref_final": "total"}) \
    .to_csv(OUT / "genus_lens_recalc_dropin.tsv", sep="\t", index=False, float_format="%.1f")

changed = val.dropna(subset=["ge85_old_pct"]).copy()
changed["d_ge85"] = changed["ge85_new_pct"] - changed["ge85_old_pct"]
top_changed = changed.reindex(changed["d_ge85"].abs().sort_values(ascending=False).index).head(15)

with open(OUT / "REPORT.md", "w", encoding="utf-8") as f:
    f.write("# genus_lens 重算 (补充分段病毒) 与 ±15% 规则验证\n\n")
    f.write(f"- 完整基因组 isolate: {len(iso)} (分段 {len(iso_seg)} / 非分段 {len(iso_non)}), 属 {len(g)} 个\n")
    f.write(f"- 分段病毒聚合规则: 同 Taxid 各 segment 求和 (Segment 标签去重取最长)\n")
    f.write(f"- 规则: 拯救线 = 属均值 × {GENUS_TOL} (允许短 {int((1-GENUS_TOL)*100)}%); 双侧 ±{int(TWO_SIDED*100)}% 命中率\n\n")
    f.write("## 旧口径判定\n\n")
    f.write(f"- ratio_old_mean = 旧 genus_lens / 新 isolate 均值, 中位 **{both['ratio_old_mean'].median():.3f}**\n")
    f.write(f"- ratio_acc_mean = 逐 accession 均值 / 新 isolate 均值, 中位 **{both['ratio_acc_mean'].median():.3f}**\n")
    f.write(f"- 分段为主属 (segmented≥50%, {len(seg_heavy)} 个) 的 ratio_old_mean 中位 {seg_heavy['ratio_old_mean'].median():.3f}"
            f" —— 明显 <1 即说明旧值按 segment 逐条平均、被分段拉低\n\n")
    f.write("## 15% 规则命中率 (n≥3 的属)\n\n")
    f.write("| 口径 | ≥85% 命中率中位 | ±15% 命中率中位 |\n|---|---|---|\n")
    f.write(f"| 旧 genus_lens | {val['ge85_old_pct'].median():.1f}% | {val['within15_old_pct'].median():.1f}% |\n")
    f.write(f"| 新均值 | {val['ge85_new_pct'].median():.1f}% | {val['within15_new_pct'].median():.1f}% |\n")
    f.write(f"| 新规则(分段属取最短完整基因组) | {val['ge85_final_pct'].median():.1f}% | {val['within15_final_pct'].median():.1f}% |\n")
    f.write(f"\n分段为主属 ({len(seg_val)} 个): ≥85% 旧 {seg_val['ge85_old_pct'].median():.1f}% → 均值 {seg_val['ge85_new_pct'].median():.1f}% "
            f"→ 取最短 {seg_val['ge85_final_pct'].median():.1f}%; "
            f"±15% 旧 {seg_val['within15_old_pct'].median():.1f}% → 均值 {seg_val['within15_new_pct'].median():.1f}% → 取最短 {seg_val['within15_final_pct'].median():.1f}%\n\n")
    f.write("## 参考长度规则 (本次采用)\n\n")
    f.write("- 分段为主属 (segmented≥50%): **ref = 最短完整段组基因组** —— 先以属内段数众数定义\"完整段组\""
            "(排除只录了部分 segment 的 taxid), 池内再剔除 <0.6×池中位数的注释污染记录"
            "(如挂在 Begomovirus 属名下的 DNA-beta/betasatellite 1357nt, EU384593.1), 最后取 **min**。"
            "属内单/双组分双峰 (如 Begomovirus) 时以短的完整型为保留基准, 长的 (A+B 双组分) 因分支 D 为单侧 ≥85% 天然通过。\n")
    f.write(f"- 其余属: ref = 均值。\n")
    f.write(f"- 完整段组内最短仍 <中位数 0.5 倍的属 ({len(suspect)} 个, 建议人工复核): "
            + (", ".join(suspect["genus"].tolist()) or "无") + "\n\n")
    ex = stat[stat["genus"].isin(["Begomovirus", "Crinivirus", "Fijivirus", "Ilarvirus", "Orthotospovirus",
                                  "Ophiovirus", "Potyvirus", "Polerovirus"])]
    f.write("代表属对照 (old / 均值 / 全局min / 最短完整段组 → 最终 ref):\n\n")
    f.write(ex[["genus", "n_isolates", "n_segmented", "old_len", "new_mean", "new_min", "min_complete", "ref_final"]]
            .to_markdown(index=False, floatfmt=".0f"))
    f.write("\n\n")
    f.write("## 不完整记录的适应\n\n")
    f.write(f"- ref_info 共 {len(ref_len_ok)} 条有效长度记录, 其中 partial/fragment {n_partial} 条 "
            f"({n_partial/len(ref_len_ok)*100:.0f}%, 长度中位 ~1.1-1.7kb)。**本重算只取 Nuc_Completeness==complete 的 "
            f"{len(comp)} 条**, partial 一律不进属均值。\n")
    f.write(f"- 若把 partial 混回逐 accession 均值, 属均值中位再被拖低至 **{dp.median():.2f}×** "
            f"(如 Potyvirus 9038→6519, -28%)。\n")
    f.write(f"- 分段集不全的 \"complete\" 记录 (属内只录了部分 segment 的 taxid) 由段数众数规则排除在 min 之外。\n")
    f.write(f"- 仅有不完整记录、故无重算值的属 {len(partial_only)} 个 (旧表有值 "
            f"{sum(g_ in old_plain_set for g_ in partial_only)} 个, VIGA 有值 {sum(g_ in viga_plain_set for g_ in partial_only)} 个) "
            f"-- 这些属不在 drop-in 中, 分支 D 将落入 KEEP_no_genus_len 兜底。\n\n")
    f.write("## 三方表格交叉对比 (MMPV 旧表 vs VIGA vs 重算)\n\n")
    f.write(f"- 三方均有值 {len(c3)} 属; old∩VIGA {len(stat.dropna(subset=['old_len','viga_len']))} 属, "
            f"VIGA∩重算 {len(stat.dropna(subset=['viga_len','new_mean']))} 属。\n")
    f.write(f"- old/VIGA 直接对比: 中位比 {rv.median():.3f}, 偏差>10% 的属 {int(((rv-1).abs()>0.1).sum())}/{len(rv)} "
            f"(两表参考库不同, 属均值本身就有出入)。\n")
    for tag, sub_ in [("非分段属", c3[c3["n_segmented"]/c3["n_isolates"] < 0.5]),
                      ("分段为主属", c3[c3["n_segmented"]/c3["n_isolates"] >= 0.5])]:
        r1 = (sub_["old_len"]/sub_["new_mean"]).dropna()
        r2 = (sub_["viga_len"]/sub_["new_mean"]).dropna()
        f.write(f"- {tag} ({len(sub_)} 属): old/重算均值 中位 **{r1.median():.3f}**, VIGA/重算均值 中位 **{r2.median():.3f}**\n")
    f.write("\nold 与 VIGA 分歧最大的 10 个属:\n\n")
    f.write(top_div[["genus", "old_len", "viga_len", "new_mean", "n_segmented"]].to_markdown(index=False, floatfmt=".0f"))
    f.write("\n\n## 变化最大的 15 个属\n\n")
    f.write(top_changed[["genus", "n_isolates", "n_segmented", "old_len", "new_mean", "ge85_old_pct", "ge85_new_pct"]]
            .to_markdown(index=False, floatfmt=".0f"))
    f.write("\n\n## 旧表离谱值 (ratio_old_mean >3 或 <0.25)\n\n")
    f.write(outliers[["genus", "n_isolates", "old_len", "new_mean", "ratio_old_mean"]]
            .to_markdown(index=False, floatfmt=".1f"))
    f.write("\n\n## 属内分段架构 (n_taxid≥10)\n\n")
    f.write("单一段数 ≠ 双段数 的属存在双峰基因组架构 (如 Begomovirus 单组分 vs A+B 双组分), "
            "单一属均值对另一架构偏严, 分支 D 建议对这类属改用中位数或按架构分组:\n\n")
    f.write(seg_dist.to_markdown(index=False))
    f.write("\n\n## 下游利用与判定修改 (节段属感知)\n\n")
    f.write("判定规则: **先查属是否节段 (is_segmented 列, seg_frac>=0.5); 节段属要求 ref_min*0.85 <= 长度 <= ref_max*1.15"
            " (上下界各容忍 15%); 非节段属或表中无该属 → 沿用旧规则 长度 >= 属均值*0.85 (单侧)**。\n\n")
    f.write("| 脚本 | 用法 | 本次修改 |\n|---|---|---|\n")
    f.write("| rescue_pipeline.py 分支 D | len >= 属均值×0.85 拯救 | ✅ 改为节段范围判定 (branch_d 内自动加载 v2 列; 老表自动回退) |\n")
    f.write("| rescue_pipeline.py 分支 B | CheckV NA 但 len ≥ 均值×0.85 → genus_rescued | ✅ 同上 (seg_meta 参数) |\n")
    f.write("| virome_pipeline.py [9a] 建树过滤 | len_ratio≥0.7 才进属树 | ✅ 节段属 ratio 以 ref_min 为基准 (单段 contig 不再被误滤) |\n")
    f.write("| virome_pipeline.py [9b] KEEP 过滤表 | ratio<0.7 记 FILTERED | ✅ 同上 |\n")
    f.write("| integrate_rescue_evidence.py | core domain 校验长度门槛 (≥均值×20%) | ✅ 默认路径加 v2; 修复多列行 float 解析 |\n")
    f.write("| report_pipeline.py | 仅分支标签显示 | 无需改 |\n")
    f.write("| scripts/blacklist/bl_layer2c_trees.py | len_ratio≥0.7 建树 (黑名单侧脚本) | 建议照 [9b] 同样改, 未动 |\n\n")
    f.write("v2 表 `database/genus_lens_v2.tsv` 列: genus / total (旧口径兼容: 非节段=均值, 节段=最短完整段组) / "
            "is_segmented (yes/no) / ref_min / ref_max (节段属完整单元范围: 最短完整段到最长基因组, 上界封顶中位×2.5) / "
            "n_isolates / n_segmented。所有 loader 按列索引取值对多列天然兼容; 路径回退链 v2 → 旧表。\n\n")
    f.write("示例: Begomovirus [2491, 6903] → 接受 2118–7938bp (单段 DNA-A/B 2.5–2.9k ✓, 双组分全基因组 5.2k ✓, "
            "碎片 1.5k ✗, 混装异常 12k ✗); Ilarvirus [1956, 8807] → RNA3 单段 2.2k 从 ratio 0.28 (滤) 变 1.12 (留)。\n\n")
    f.write("\n\n## 文件\n\n- genus_lens_recalc.tsv 全量对照; validation_15pct.tsv 命中率; "
            "genus_lens_recalc_dropin.tsv (g__ 前缀 2 列, 兼容旧 loader); "
            "**database/genus_lens_v2.tsv** (genus/total/is_segmented/ref_min/ref_max/n_isolates/n_segmented, "
            "旧 loader 只读前 2 列向后兼容; 下游分支 D/B 据此做节段范围判定)\n")

print(f"\n输出: {OUT}")
print("Top 15 变化属:")
print(top_changed[["genus", "n_isolates", "n_segmented", "old_len", "new_mean", "ge85_old_pct", "ge85_new_pct"]].to_string(index=False, float_format=lambda x: f"{x:.0f}"))
