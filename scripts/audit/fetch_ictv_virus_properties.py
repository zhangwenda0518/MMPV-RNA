# -*- coding: utf-8 -*-
"""
Fetch ICTV "Virus Properties by Family" and materialise two audit TSVs.

Scope / anchoring
  * Taxonomy anchor: ICTV MSL41 (2025-2026 release, Family = 427).
  * Data source: https://ictv.global/report/information/virus-properties
    (paginated; we request items_per_page=100 and walk pages until empty).
  * Ground truth for "current family or not" cross-checked against
    ICTV_Master_Species_List_2025_MSL41.v1.xlsx (Family count = 427).

Outputs (both under scripts/audit/):
  1. ICTV_HOST_SWEEP_233.tsv      columns: Family, ICTV_Host, Plants_Included, Source_URL, Note
  2. ICTV_PLANT_HOST_FLAG.tsv     columns: Family, Category, ICTV_Host, Note
                                  Category in {plant_included, host_undetermined, not_current_family}

This script only fetches + writes the two TSVs. It does not touch any
pipeline code (e.g. run_host_prediction.py).
"""

import csv
import os
import sys

import requests
from bs4 import BeautifulSoup

HERE = os.path.dirname(os.path.abspath(__file__))
SOURCE_URL = "https://ictv.global/report/information/virus-properties?items_per_page=100"
BASE = "https://ictv.global/report/information/virus-properties"
HEADERS = {
    "User-Agent": (
        "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 "
        "(KHTML, like Gecko) Chrome/124.0 Safari/537.36"
    )
}

# ---------------------------------------------------------------------------
# Target family list (as supplied for the sweep; 237 distinct names).
# ---------------------------------------------------------------------------
TARGET_NAMES = """
Ackermannviridae Adamaviridae Adomaviridae Ahmunviridae Aliceevansviridae Aliusviridae Alloherpesviridae Allomimiviridae
Alphaormycoviridae Alphatetraviridae Alternaviridae Amnoonviridae Anamaviridae Andersonviridae Anelloviridae Aoguangviridae
Apasviridae Arenaviridae Arenbergviridae Arteriviridae Artiviridae Artoviridae Ascoviridae Asfarviridae
Assiduviridae Astroviridae Autographiviridae Autonotataviridae Autoscriptoviridae Autosignataviridae Autotranscriptaviridae Bacilladnaviridae
Barnaviridae Belpaoviridae Berryhillviridae Betaormycoviridae Bidnaviridae Birnaviridae Blumeviridae Bornaviridae
Botybirnaviridae Buchnerviridae Caliciviridae Carmotetraviridae Casidaviridae Casjensviridae Chaseviridae Chimalliviridae
Chrysoviridae Chuviridae Circoviridae Clermontviridae Colingsworthviridae Connertonviridae Coronaviridae Cremegaviridae
Crevaviridae Cruciviridae Cruliviridae Curvulaviridae Darmviridae Deltaflexiviridae Deltaormycoviridae Demerecviridae
Dicistroviridae Discoviridae Dishuiviroviridae Draupnirviridae Drexlerviridae Dumbiviridae Duneviridae Ehrlichviridae
Endolinaviridae Epsomviridae Eupolintoviridae Euroniviridae Fervensiviridae Fiersviridae Filamentoviridae Filoviridae
Fusariviridae Fuselloviridae Gammaflexiviridae Gammaormycoviridae Gandrviridae Genomoviridae Giardiaviridae Graaviviridae
Grimontviridae Hadakaviridae Hafunaviridae Haloferuviridae Halomagnusviridae Halspiviridae Hantaviridae Helgolandviridae
Hepadnaviridae Hepaciviridae Hepeviridae Herelleviridae Hodgkinviridae Hydriviridae Hypoviridae Hytrosaviridae
Iflaviridae Inoviridae Inseviridae Intestiviridae Iridoviridae Itzamnaviridae Jeanschmidtviridae Kanorauviridae
Kirkoviridae Kleczkowskaviridae Kolmioviridae Konodaiviridae Krittikaviridae Kyanoviridae Lakviridae Lebotiviridae
Leishbuviridae Lindbergviridae Lipothrixviridae Lispiviridae Ludisviridae Malacoherpesviridae Mamonoviridae Marnaviridae
Marseilleviridae Matonaviridae Maviroviridae Medioniviridae Megabirnaviridae Megatotiviridae Mesomimiviridae Mesoniviridae
Mesyanzhinovviridae Microviridae Mimiviridae Mitoviridae Monocitiviridae Mononiviridae Mtkvariviridae Mycoalphaviridae
Mymonaviridae Nairoviridae Nanghoshaviridae Nanhypoviridae Narnaviridae Naryaviridae Natareviridae Nenyaviridae
Nimaviridae Nipumfusiviridae Nodaviridae Noraviridae Nudiviridae Nyamiviridae Oomyviridae Ootiviridae
Orpheoviridae Orthoherpesviridae Orthomyxoviridae Ouroboviridae Ourmiaviridae Pachyviridae Pandoraviridae Papillomaviridae
Parahypoviridae Paramyxoviridae Parvoviridae Pecoviridae Peduoviridae Permutotetraviridae Pestiviridae Phasmaviridae
Phlegiviridae Phycodnaviridae Phypoliviridae Picobirnaviridae Picornaviridae Pithoviridae Pneumoviridae Polycipiviridae
Polymycoviridae Polyomaviridae Pootjesviridae Poxviridae Pseudototiviridae Pseudoviridae Pyrstoviridae Quadriviridae
Redondoviridae Rhizouliviridae Roniviridae Rountreeviridae Salasmaviridae Sarkviridae Schitoviridae Schizomimiviridae
Shortaselviridae Sinhaliviridae Smacoviridae Solinviviridae Spiciviridae Splipalmiviridae Sputniviroviridae Stackebrandtviridae
Stanwilliamsviridae Steigviridae Steitzviridae Straboviridae Sunviridae Suoliviridae Tectiviridae Tobaniviridae
Tobaliviridae Tomosaviridae Toyamaviridae Trimbiviridae Tulasviridae Vandenendeviridae Vertoviridae Vilmaviridae
Winoviridae Wupedeviridae Xinmoviridae Yadokariviridae Yangangviridae Yaraviridae Yueviridae Zierdtviridae
Zimmerviridae Zobellviridae Reoviridae Virgaviridae Solemoviridae
""".split()

# Spelling aliases: name as supplied in the sweep list -> official MSL41 name.
ALIASES = {
    "Fervensiviridae": "Fervensviridae",  # sweep list typo; official = Fervensviridae
}

# Reasons for families absent from the MSL41 properties table (ground-truthed
# against ICTV_Master_Species_List_2025_MSL41.v1.xlsx and ICTV taxon reports).
NOT_CURRENT_REASON = {
    "Adomaviridae": "MSL41 无此科；系文献提出的 proposed 科名（鱼类/软骨鱼相关 DNA 病毒），未获 ICTV 采纳",
    "Autographiviridae": "MSL41 无此科；2025 年由科提升为目 Autographivirales（提案 2024.045B，含 4 个新科），原科名不再使用",
    "Cruciviridae": "MSL41 无此科；系 proposed 科名（CRESS-DNA 病毒类群），未获 ICTV 采纳",
    "Microviridae": "MSL41 无此科；2025 年已废止并升格为纲 Microviricetes（提案 2025.043B），原亚科升为目 Bullavirales/Gokushovirales",
    "Pandoraviridae": "MSL41 无此科（全表亦无 Pandoravirus 属）；系 proposed 巨病毒科名，未获 ICTV 采纳",
    "Reoviridae": "MSL41 无此科；2021 年（MSL36）已拆分为 Sedoreoviridae 与 Spinareoviridae（二者均为现行科）",
}

UNCERTAIN_TOKENS = ("predicted", "uncertain", "not determined", "undetermined", "soil")


def fetch_table():
    """Return (ordered list of (family, host), raw row count) from the ICTV page."""
    rows = []
    page = 0
    while True:
        resp = requests.get(
            BASE,
            params={"items_per_page": 100, "page": page},
            headers=HEADERS,
            timeout=60,
        )
        resp.raise_for_status()
        soup = BeautifulSoup(resp.text, "html.parser")
        table = soup.find("table")
        if table is None:
            break
        trs = table.find_all("tr")[1:]
        if not trs:
            break
        for tr in trs:
            tds = [td.get_text(" ", strip=True) for td in tr.find_all("td")]
            if len(tds) < 11:
                continue
            rows.append((tds[1], tds[10]))  # col 1 = Family, col 10 = Host
        page += 1
    return rows


def classify_plants(host):
    """Return (Plants_Included, note_suffix) for a found family."""
    h = host.lower()
    if "plant" in h:
        return "是", ""
    if any(tok in h for tok in UNCERTAIN_TOKENS) or host.strip() == "":
        return "未确定", "官方 Host 为预测值/环境来源，非确定宿主，故 Plants_Included 记未确定"
    return "否", ""


def main():
    rows = fetch_table()
    tbl = dict(rows)
    print("fetched families:", len(rows))

    uniq = list(dict.fromkeys(TARGET_NAMES))
    print("target names:", len(TARGET_NAMES), "distinct:", len(uniq))

    sweep = []          # [1]
    flag = []           # [2]
    plant_list = []

    for name in uniq:
        official = ALIASES.get(name, name)
        if official in tbl:
            host = tbl[official]
            plants, suffix = classify_plants(host)
            if plants == "是":
                note = ""
                plant_list.append(name)
            elif plants == "未确定":
                note = suffix
            else:
                note = ""
            if name != official:
                note = (
                    f"官方 MSL41 科名为 {official}（现行科，宿主 {tbl[official]}）；"
                    f"名单拼写 {name} 疑为笔误，已按官方科名核对填值。" + (note or "")
                )
            sweep.append([name, host, plants, SOURCE_URL, note])

            if plants == "是":
                flag.append([name, "plant_included", host, "官方 Host 列原文含 plants（MSL41）"])
            elif plants == "未确定":
                flag.append([name, "host_undetermined", host,
                             f"官方 Host 为预测值/环境来源（{host}），非确定宿主"])
        else:
            note = NOT_CURRENT_REASON.get(name, "MSL41 427 科表中未收录该名称")
            sweep.append([name, "未收录", "未确定", SOURCE_URL, note])
            flag.append([name, "not_current_family", "未收录", note])

    # ---- write [1] ----
    sweep_path = os.path.join(HERE, "ICTV_HOST_SWEEP_233.tsv")
    with open(sweep_path, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["Family", "ICTV_Host", "Plants_Included", "Source_URL", "Note"])
        w.writerows(sweep)

    # ---- write [2] ----
    order = {"plant_included": 0, "host_undetermined": 1, "not_current_family": 2}
    flag.sort(key=lambda r: order[r[1]])
    flag_path = os.path.join(HERE, "ICTV_PLANT_HOST_FLAG.tsv")
    with open(flag_path, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["Family", "Category", "ICTV_Host", "Note"])
        w.writerows(flag)

    print("sweep rows:", len(sweep), "->", sweep_path)
    print("flag rows:", len(flag), "->", flag_path)
    print("plant_included:", len(plant_list))
    for n in plant_list:
        print("  ", n)
    return 0


if __name__ == "__main__":
    sys.exit(main())
