#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
补丁: 为 run_host_prediction.py 的 Family 级宿主判定加入非植物黑名单。

背景:
  C9 (ICTV) 按 Determination_Level=Family(科级) 判宿主, 科级注释无区分力,
  导致跨宿主科 (如 Partitiviridae) 的非植物属 (Biavirus 等) 被误判为 Plant。

数据来源:
  ~/plant_virus_db/1.virus-host_db/C-host_classify/classified_clean/Plant.tsv
  (PlantVirusDB 权威植物病毒宿主库, 51 合法科 / 2769 合法属)

规则:
  1. 整科黑名单: C9 判 Plant 但该科在权威库 Plant.tsv 完全不存在 → 否决
  2. 属黑名单:   C9 判 Plant 但该属在权威库 Plant.tsv 完全不存在 → 否决
  均只在 C9 科级/目级判定 (Determination_Level != Genus/Species) 时生效,
  属/种级判定保留 C9 原判, 避免误杀。
"""
import argparse
import shutil
import sys
from pathlib import Path

TARGET = Path("/home/zhangwenda/MMPV-RNA/virome_discovery_pipeline/run_host_prediction.py")

# ---------- 1) 整科黑名单 (权威库 Plant.tsv 中完全没有的科) ----------
NON_PLANT_FAMILIES = {
    # 巨型 DNA 病毒 (原生生物/藻类宿主)
    "Mimiviridae", "Phycodnaviridae", "Marseilleviridae",
    "Pithoviridae", "Mesomimiviridae", "Allomimiviridae",
    "Mamonoviridae", "Schizomimiviridae", "Maviroviridae",
    "Orpheoviridae", "Pandoraviridae",
    # 动物 / 昆虫病毒
    "Orthoherpesviridae", "Nimaviridae", "Astroviridae", "Baculoviridae",
    "Adenoviridae", "Poxviridae", "Iridoviridae", "Alloherpesviridae",
    "Papillomaviridae", "Picornaviridae", "Hepaciviridae", "Nudiviridae",
    "Dicistroviridae", "Arenaviridae", "Herelleviridae", "Picobirnaviridae",
    "Peduoviridae", "Lymphocystiviridae",
    # 真菌 / 真菌-植物跨界病毒
    "Tobaliviridae", "Ambiguiviridae", "Ourmiaviridae", "Deltaflexiviridae",
    "Gammaflexiviridae", "Hypoviridae", "Mitoviridae", "Pseudoviridae",
    "Genomoviridae", "Gandrviridae",
}

# ---------- 2) 属黑名单 (权威库 Plant.tsv 中完全没有的属) ----------
NON_PLANT_GENERA = {
    # Partitiviridae 跨宿主科内的非植物属
    "Biavirus", "Criusvirus", "Tethysvirus", "Klosneuvirus", "Hokovirus",
    "Catovirus", "Indivirus", "Delepquintavirus", "Clandestinovirus",
    "Bimevirus", "Chaphamavirus", "Prymnesiovirus",
    # Caulimoviridae / Tombusviridae / Geminiviridae 等跨宿主科
    "Sylvanvirus", "Rimosavirus", "Miraophiovirus", "Punuivirus",
    "Oleurovirus", "Jingmenvirus", "Hatfieldvirus", "Alphaorpheovirus",
    "Hubsclerovirus", "Geminivirus", "Rheavirus", "Patrovirus",
    "Zetanudivirus", "Solivirus",
    # 巨型 / 藻类病毒属
    "Chlorovirus", "Prasinovirus", "Coccolithovirus", "Theiavirus",
    "Varicellovirus", "Alphapithovirus", "Jouyvirus", "Tupanvirus",
    "Megavirus", "Kratosvirus", "Mavirus",
    # 动物 / 昆虫病毒属
    "Cytorhabdovirus", "Gammaretrovirus", "Bandavirus", "Avastrovirus",
    "Alphabaculovirus", "Lymphocystivirus", "Arenavirus", "Cypovirus",
    "Orthoflavivirus", "Iltovirus", "Simplexvirus", "Orthoreovirus",
    "Varicellovirus", "Grisebachstrassevirus", "Inibicvirus",
    "Timquatrovirus", "Velavnirus", "Velanvirus",
    # 其他
    "Alphambiguivirus", "Betambiguivirus", "Gammambiguivirus",
    "Tobalivirus", "Crucivirus",
}

BLOCK = '''
# ----------------- 非植物科属黑名单 (Family 级判定门控) -----------------
# 数据源: PlantVirusDB 权威库 classified_clean/Plant.tsv
#   (51 合法植物科 / 2769 合法植物属)
# 规则: C9 判定级别 < Genus 时 (Determination_Level 为 Family/Order/None),
#       若科或属不在权威库中, 则否决 Plant 判定, 交由下游 RVH/PB2 或标 Unknown。
NON_PLANT_FAMILIES = {NON_PLANT_FAMILIES_REPR}

NON_PLANT_GENERA = {NON_PLANT_GENERA_REPR}

# 视为"属/种级"的判定级别 (这些级别的 C9 判定可信, 不触发黑名单)
TRUSTED_LEVELS = {{'genus', 'species', 'species*'}}


def is_trusted_level(det_level):
    """Determination_Level 是否为属/种级 (可信级别)。"""
    if det_level is None:
        return False
    d = str(det_level).split('(')[0].strip().lower()
    return d in TRUSTED_LEVELS


def is_blacklisted(family, genus, det_level):
    """科级判定时命中黑名单 → True (应否决 Plant)。"""
    if is_trusted_level(det_level):
        return False
    fam = str(family).strip() if family is not None else ''
    gen = str(genus).strip() if genus is not None else ''
    if fam and fam in NON_PLANT_FAMILIES:
        return True
    if gen and gen in NON_PLANT_GENERA:
        return True
    return False

'''


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--apply", action="store_true", help="真正写入 (默认 dry-run)")
    args = p.parse_args()

    if not TARGET.is_file():
        sys.exit(f"ERROR: 找不到 {TARGET}")

    src = TARGET.read_text(encoding="utf-8")

    # --- 幂等检查 ---
    if "NON_PLANT_FAMILIES" in src:
        sys.exit("已打过补丁 (发现 NON_PLANT_FAMILIES), 退出。")

    block = BLOCK.format(
        NON_PLANT_FAMILIES_REPR=repr(sorted(NON_PLANT_FAMILIES)),
        NON_PLANT_GENERA_REPR=repr(sorted(NON_PLANT_GENERA)),
    )

    # ---------- 改动 1: 在 PHAGE_CLASSES 定义后插入黑名单块 ----------
    anchor1 = "PHAGE_CLASSES = {'leviviricetes', 'caudoviricetes', 'vidaverviricetes', 'faserviricetes'}\n"
    if anchor1 not in src:
        sys.exit("ERROR: 找不到 anchor1 (PHAGE_CLASSES 定义)")
    src = src.replace(anchor1, anchor1 + block, 1)

    # ---------- 改动 2: run_ensemble 的 c9 merge 扩展列 ----------
    anchor2 = (
        "        c9 = pd.read_csv(c9_tsv, sep='\\t')[['contig_id', 'Predicted_Host']]"
        ".rename(columns={'Predicted_Host': 'Host_ICTV'})\n"
    )
    new2 = (
        "        _c9_cols = ['contig_id', 'Predicted_Host']\n"
        "        _c9_extra = [c for c in ('Family', 'Genus', 'Determination_Level')\n"
        "                     if c in pd.read_csv(c9_tsv, sep='\\t', nrows=0).columns]\n"
        "        _c9_cols += _c9_extra\n"
        "        c9 = pd.read_csv(c9_tsv, sep='\\t')[_c9_cols]"
        ".rename(columns={'Predicted_Host': 'Host_ICTV'})\n"
    )
    if anchor2 not in src:
        sys.exit("ERROR: 找不到 anchor2 (c9 merge)")
    src = src.replace(anchor2, new2, 1)

    # ---------- 改动 3: decision_tree_cascade 加黑名单门控 ----------
    anchor3 = (
        "    h_ictv = normalize_c9(row.get('Host_ICTV', 'Unknown'))\n"
        "    if h_ictv != 'Unknown':\n"
        "        return h_ictv, 'ICTV_Preferred'\n"
    )
    new3 = (
        "    h_ictv = normalize_c9(row.get('Host_ICTV', 'Unknown'))\n"
        "    if h_ictv != 'Unknown':\n"
        "        # 黑名单门控: C9 科级判定 + 命中最植物黑名单 → 否决 Plant\n"
        "        if h_ictv == 'Plant' and is_blacklisted(\n"
        "                row.get('Family'), row.get('Genus'),\n"
        "                row.get('Determination_Level')):\n"
        "            h_ictv = 'Unknown'\n"
        "        else:\n"
        "            return h_ictv, 'ICTV_Preferred'\n"
    )
    if anchor3 not in src:
        sys.exit("ERROR: 找不到 anchor3 (cascade ICTV 分支)")
    src = src.replace(anchor3, new3, 1)

    # ---------- 改动 4: decision_tree(旧版) 同步 ----------
    anchor4 = (
        "    h_ictv = normalize_c9(row.get('Host_ICTV', 'Unknown'))\n"
        "    h_rvh  = parse_rvh(row)\n"
    )
    new4 = (
        "    h_ictv = normalize_c9(row.get('Host_ICTV', 'Unknown'))\n"
        "    if h_ictv == 'Plant' and is_blacklisted(\n"
        "            row.get('Family'), row.get('Genus'),\n"
        "            row.get('Determination_Level')):\n"
        "        h_ictv = 'Unknown'\n"
        "    h_rvh  = parse_rvh(row)\n"
    )
    if anchor4 not in src:
        sys.exit("ERROR: 找不到 anchor4 (decision_tree)")
    src = src.replace(anchor4, new4, 1)

    if not args.apply:
        print("=== DRY RUN ===")
        print(f"整科黑名单: {len(NON_PLANT_FAMILIES)} 个")
        print(f"属黑名单:   {len(NON_PLANT_GENERA)} 个")
        print("4 处改动全部 anchor 命中。加 --apply 写入。")
        return

    bak = TARGET.with_suffix(".py.bak_blacklist2_20260902")
    shutil.copy2(TARGET, bak)
    TARGET.write_text(src, encoding="utf-8")
    print(f"备份: {bak}")
    print(f"写入: {TARGET}")
    print(f"整科黑名单 {len(NON_PLANT_FAMILIES)} 个 / 属黑名单 {len(NON_PLANT_GENERA)} 个")


if __name__ == "__main__":
    main()
