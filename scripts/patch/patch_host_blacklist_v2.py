#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
补丁 v2: run_host_prediction.py 加 Family 级宿主判定的【属黑名单】门控。

核心口径 (用户确认):
  科级判定没有区分力 (同一科里植物属与非植物属混杂),
  所以当 C9 只报到科级时, 判据落在【属】上:
    - 属在非植物属黑名单 → 否决 Plant
    - 属不在黑名单 (或在权威库中) → 放行
  属/种级判定 (C9 有区分力的层) 不触发黑名单, 完全保留原判。

数据源:
  ~/plant_virus_db/1.virus-host_db/C-host_classify/classified_clean/Plant.tsv
"""
import argparse
import shutil
import sys
from pathlib import Path

TARGET = Path("/home/zhangwenda/MMPV-RNA/virome_discovery_pipeline/run_host_prediction.py")

# ---------- 属黑名单: 权威库 Plant.tsv 中完全没有的属 ----------
NON_PLANT_GENERA = {
    # Partitiviridae 跨宿主科内的非植物属 (Biavirus 485 条为最大项)
    "Biavirus", "Criusvirus", "Tethysvirus", "Klosneuvirus", "Hokovirus",
    "Catovirus", "Indivirus", "Delepquintavirus", "Clandestinovirus",
    "Bimevirus", "Chaphamavirus", "Prymnesiovirus",
    # Caulimoviridae / Tombusviridae / Geminiviridae 等跨宿主科的属级误判
    "Sylvanvirus", "Rimosavirus", "Miraophiovirus", "Punuivirus",
    "Oleurovirus", "Jingmenvirus", "Hatfieldvirus", "Alphaorpheovirus",
    "Hubsclerovirus", "Geminivirus", "Rheavirus", "Patrovirus",
    "Zetanudivirus", "Solivirus",
    # 巨型 DNA / 藻类病毒属
    "Chlorovirus", "Prasinovirus", "Coccolithovirus", "Theiavirus",
    "Varicellovirus", "Alphapithovirus", "Jouyvirus", "Tupanvirus",
    "Megavirus", "Kratosvirus", "Mavirus",
    # 动物 / 昆虫病毒属
    "Gammaretrovirus", "Bandavirus", "Avastrovirus", "Alphabaculovirus",
    "Lymphocystivirus", "Arenavirus", "Cypovirus", "Orthoflavivirus",
    "Iltovirus", "Simplexvirus", "Orthoreovirus",
    "Grisebachstrassevirus", "Inibicvirus", "Timquatrovirus", "Velanvirus",
    # 真菌/跨界病毒属
    "Alphambiguivirus", "Betambiguivirus", "Gammambiguivirus",
    "Tobalivirus", "Crucivirus",
}

# ---------- 科黑名单: 仅用于「无属信息」时的兜底 ----------
# 注意: 单独用科黑名单会误杀 (Mimiviridae 列下有 Fabavirus 等真植物病毒),
# 因此仅在 Genus 缺失时使用。
NON_PLANT_FAMILIES_FALLBACK = {
    "Mimiviridae", "Phycodnaviridae", "Marseilleviridae", "Pithoviridae",
    "Mesomimiviridae", "Allomimiviridae", "Mamonoviridae",
    "Schizomimiviridae", "Maviroviridae", "Orpheoviridae", "Pandoraviridae",
    "Orthoherpesviridae", "Nimaviridae", "Astroviridae", "Baculoviridae",
    "Adenoviridae", "Poxviridae", "Iridoviridae", "Alloherpesviridae",
    "Papillomaviridae", "Picornaviridae", "Hepaciviridae", "Nudiviridae",
    "Dicistroviridae", "Arenaviridae", "Herelleviridae", "Picobirnaviridae",
    "Peduoviridae", "Tobaliviridae", "Ambiguiviridae", "Ourmiaviridae",
    "Deltaflexiviridae", "Gammaflexiviridae", "Hypoviridae", "Mitoviridae",
    "Pseudoviridae", "Genomoviridae", "Gandrviridae",
}

BLOCK = '''
# --------------- 非植物属黑名单 (Family 级判定门控) ---------------
# 数据源: PlantVirusDB 权威库 classified_clean/Plant.tsv
# 口径: 科级判定无区分力 (同科内植物属与非植物属混杂),
#       故 C9 只报科级时, 以【属】为判据:
#         属在黑名单 → 否决 Plant; 否则放行 (含属在权威库的情形)。
#       属/种级判定不触发黑名单。
NON_PLANT_GENERA = {NON_PLANT_GENERA_REPR}

# 无属信息时的整科兜底黑名单 (权威库中无此科的记录)
NON_PLANT_FAMILIES_FALLBACK = {NON_PLANT_FAMILIES_FALLBACK_REPR}

TRUSTED_LEVELS = {{'genus', 'species', 'species*'}}


def is_trusted_level(det_level):
    """Determination_Level 是否为属/种级 (C9 有区分力的层级)。"""
    if det_level is None:
        return False
    d = str(det_level).split('(')[0].strip().lower()
    return d in TRUSTED_LEVELS


def is_blacklisted(family, genus, det_level):
    """科级判定时按属黑名单否决 Plant; 无属信息时按科兜底。"""
    if is_trusted_level(det_level):
        return False
    gen = str(genus).strip() if genus is not None else ''
    fam = str(family).strip() if family is not None else ''
    if gen and gen not in ('NA', 'nan'):
        return gen in NON_PLANT_GENERA
    # 属信息缺失: 整科兜底
    return bool(fam and fam in NON_PLANT_FAMILIES_FALLBACK)

'''


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--apply", action="store_true")
    args = p.parse_args()

    if not TARGET.is_file():
        sys.exit(f"ERROR: 找不到 {TARGET}")

    src = TARGET.read_text(encoding="utf-8")
    if "NON_PLANT_GENERA" in src:
        sys.exit("已打过补丁, 退出。")

    block = BLOCK.format(
        NON_PLANT_GENERA_REPR=repr(sorted(NON_PLANT_GENERA)),
        NON_PLANT_FAMILIES_FALLBACK_REPR=repr(sorted(NON_PLANT_FAMILIES_FALLBACK)),
    )

    anchor1 = "PHAGE_CLASSES = {'leviviricetes', 'caudoviricetes', 'vidaverviricetes', 'faserviricetes'}\n"
    if anchor1 not in src:
        sys.exit("ERROR anchor1")
    src = src.replace(anchor1, anchor1 + block, 1)

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
        sys.exit("ERROR anchor2")
    src = src.replace(anchor2, new2, 1)

    anchor3 = (
        "    h_ictv = normalize_c9(row.get('Host_ICTV', 'Unknown'))\n"
        "    if h_ictv != 'Unknown':\n"
        "        return h_ictv, 'ICTV_Preferred'\n"
    )
    new3 = (
        "    h_ictv = normalize_c9(row.get('Host_ICTV', 'Unknown'))\n"
        "    if h_ictv != 'Unknown':\n"
        "        if h_ictv == 'Plant' and is_blacklisted(\n"
        "                row.get('Family'), row.get('Genus'),\n"
        "                row.get('Determination_Level')):\n"
        "            h_ictv = 'Unknown'\n"
        "        else:\n"
        "            return h_ictv, 'ICTV_Preferred'\n"
    )
    if anchor3 not in src:
        sys.exit("ERROR anchor3")
    src = src.replace(anchor3, new3, 1)

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
        sys.exit("ERROR anchor4")
    src = src.replace(anchor4, new4, 1)

    if not args.apply:
        print(f"DRY RUN: 属黑名单 {len(NON_PLANT_GENERA)} 个 / 科兜底 {len(NON_PLANT_FAMILIES_FALLBACK)} 个")
        print("4 处 anchor 命中。加 --apply 写入。")
        return

    bak = TARGET.with_suffix(".py.bak_blacklist3_20260902")
    shutil.copy2(TARGET, bak)
    TARGET.write_text(src, encoding="utf-8")
    print(f"备份: {bak}\n写入: {TARGET}")


if __name__ == "__main__":
    main()
