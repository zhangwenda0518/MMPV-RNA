# ── 2026-09-15 植物白名单 (防误否决, 只做一致性断言, 不参与否决决策) ──
# 构造 A: PlantVirusDB 权威植物库 Plant.tsv 按位置口径解析 (Virus_lineage 9 段,
#         index5=Family / index6=Genus, 与 3_host_probability.py 的 levels 一致):
#         59 科 / 202 属 (Plant.tsv md5 见 HOST_REGISTRY_AUDIT_20260915.md)
# 构造 B: ICTV 官方 Virus Properties by Family Host 列含 plants, 但 Plant.tsv 零覆盖的 9 科:
#         Artoviridae, Chrysoviridae, Genomoviridae, Kanorauviridae, Mitoviridae, Ourmiaviridae, Pestiviridae, Pseudoviridae, Spiciviridae
#         这 9 科即 2026-09-15 从 NON_PLANT_FAMILIES_FALLBACK 撒回的全部科名
#         (审计: scripts/audit/HOST_REGISTRY_AUDIT_20260915.md, 证据: ICTV_HOST_SWEEP_233.tsv)。
# 用途: (1) 机器化不变量 "白名单 交 黑名单 = 空", 在导入时 fail-loud;
#       (2) 供审计脚本/人复核时直接检索, 无需重跑 Plant.tsv 解析。
#       注意: 白名单【不覆盖科级否决】—— 科级否决是刻意的嵌合拦截
#       (例: Potyvirus 属的嵌合行被标 Mimiviridae 时按科否决),
#       若让白名单越过否决会反向制造假阳性, 故此处仅登记, 不改变否决行为。
PLANT_FAMILIES_WHITELIST = frozenset([
    'Alphaflexiviridae', 'Alphasatellitidae', 'Amalgaviridae', 'Amesuviridae', 'Anicreviridae',
    'Artoviridae', 'Aspiviridae', 'Atkinsviridae', 'Avsunviroidae', 'Benyviridae',
    'Betaflexiviridae', 'Botourmiaviridae', 'Bromoviridae', 'Caulimoviridae', 'Chrysoviridae',
    'Closteroviridae', 'Endornaviridae', 'Fimoviridae', 'Flaviviridae', 'Fusagraviridae',
    'Geminiviridae', 'Genomoviridae', 'Geplanaviridae', 'Kanorauviridae', 'Kitaviridae',
    'Konkoviridae', 'Mahapunaviridae', 'Mayoviridae', 'Metaxyviridae', 'Mitoviridae',
    'Nanoviridae', 'Orthototiviridae', 'Ourmiaviridae', 'Pamosaviridae', 'Partitiviridae',
    'Peribunyaviridae', 'Pestiviridae', 'Phenuiviridae', 'Pospiviroidae', 'Potyviridae',
    'Pseudoviridae', 'Qinviridae', 'Retroviridae', 'Rhabdoviridae', 'Secoviridae',
    'Sedoreoviridae', 'Solemoviridae', 'Solspiviridae', 'Spiciviridae', 'Spinareoviridae',
    'Togaviridae', 'Tolecusatellitidae', 'Tombusviridae', 'Tomosaviridae', 'Tonesaviridae',
    'Tospoviridae', 'Tymoviridae', 'Vilyaviridae', 'Virgaviridae',
])

PLANT_GENERA_WHITELIST = frozenset([
    'Ageyesisatellite', 'Albetovirus', 'Alfamovirus', 'Allexivirus', 'Alohovirus',
    'Alphacarmovirus', 'Alphacytorhabdovirus', 'Alphaendornavirus', 'Alphagymnorhavirus', 'Alphanecrovirus',
    'Alphanucleorhabdovirus', 'Alphapartitivirus', 'Amalgavirus', 'Ampelovirus', 'Anulavirus',
    'Apscaviroid', 'Aranruthvirus', 'Arepavirus', 'Aumaivirus', 'Aureusvirus',
    'Avenavirus', 'Avsunviroid', 'Babusatellite', 'Babuvirus', 'Badnavirus',
    'Banmivirus', 'Barchevirus', 'Becurtovirus', 'Begomovirus', 'Benyvirus',
    'Betacarmovirus', 'Betacytorhabdovirus', 'Betanecrovirus', 'Betanucleorhabdovirus', 'Betapartitivirus',
    'Betasatellite', 'Bevemovirus', 'Birfovirus', 'Blunervirus', 'Bluvavirus',
    'Brambyvirus', 'Bromovirus', 'Bymovirus', 'Capillovirus', 'Capulavirus',
    'Carlavirus', 'Caulimovirus', 'Cavemovirus', 'Celavirus', 'Cheravirus',
    'Chordovirus', 'Cilevirus', 'Citlodavirus', 'Citrivirus', 'Clecrusatellite',
    'Closterovirus', 'Clostunsatellite', 'Cocadviroid', 'Cocosatellite', 'Cofodevirus',
    'Coguvirus', 'Colecusatellite', 'Coleviroid', 'Comovirus', 'Coprasatellite',
    'Crinivirus', 'Cucumovirus', 'Curtovirus', 'Deltanucleorhabdovirus', 'Deltapartitivirus',
    'Deltasatellite', 'Dianthovirus', 'Dichorhavirus', 'Dioscovirus', 'Divavirus',
    'Elaviroid', 'Emaravirus', 'Enamovirus', 'Episkevirus', 'Eragrovirus',
    'Fabavirus', 'Fabenesatellite', 'Fijivirus', 'Foveavirus', 'Furovirus',
    'Fusagravirus', 'Gallantivirus', 'Gammacarmovirus', 'Gammacytorhabdovirus', 'Gammanucleorhabdovirus',
    'Goravirus', 'Gosmusatellite', 'Grablovirus', 'Higrevirus', 'Hordeivirus',
    'Hostuviroid', 'Idaeovirus', 'Ilarvirus', 'Ipomovirus', 'Janusivirus',
    'Ketkevirus', 'Kobbarisatellite', 'Laulavirus', 'Lolavirus', 'Lumovirus',
    'Luteovirus', 'Macanavirus', 'Machlomovirus', 'Macluravirus', 'Maculavirus',
    'Maldovirus', 'Marafivirus', 'Mastrevirus', 'Mechlorovirus', 'Menthavirus',
    'Mersevirus', 'Milvetsatellite', 'Mivedwarsatellite', 'Mulcrilevirus', 'Muscarsatellite',
    'Nanovirus', 'Nepovirus', 'Nucleorhabdovirus', 'Oleavirus', 'Olivavirus',
    'Olpivirus', 'Ophiovirus', 'Opunvirus', 'Orthotospovirus', 'Oryzavirus',
    'Ourmiavirus', 'Panicovirus', 'Papanivirus', 'Pecluvirus', 'Pelamoviroid',
    'Pelarspovirus', 'Penoulivirus', 'Petuvirus', 'Phragmivirus', 'Phytoreovirus',
    'Platypuvirus', 'Poacevirus', 'Polemovirus', 'Polerovirus', 'Pomovirus',
    'Pospiviroid', 'Potexvirus', 'Potyvirus', 'Prunevirus', 'Pteridovirus',
    'Ravavirus', 'Riddikuvirus', 'Robigovirus', 'Rosadnavirus', 'Roymovirus',
    'Rubodvirus', 'Ruflodivirus', 'Rymovirus', 'Sadwavirus', 'Sequivirus',
    'Sobemovirus', 'Solendovirus', 'Somasatellite', 'Sophoyesatellite', 'Soymovirus',
    'Stralarivirus', 'Stupevirus', 'Subclovsatellite', 'Sustrivirus', 'Temfrudevirus',
    'Tenuivirus', 'Tepovirus', 'Tobamovirus', 'Tobravirus', 'Tombusvirus',
    'Topilevirus', 'Topocuvirus', 'Torradovirus', 'Totivirus', 'Tralespevirus',
    'Trichovirus', 'Trirhavirus', 'Tritimovirus', 'Tungrovirus', 'Turncurtovirus',
    'Tymovirus', 'Umbravirus', 'Vaccinivirus', 'Varicosavirus', 'Velarivirus',
    'Virtovirus', 'Vitivirus', 'Waikavirus', 'Wamavirus', 'Wehlfuvirus',
    'Welwivirus', 'Weothlivirus', 'Whahdcavirus', 'Whiflysatellite', 'Witirovirus',
    'Yermavirus', 'Zeavirus',
])

# ICTV 现行宿主列含 plants 的 12 科 (科名 -> ICTV Host 列原文), 审计取数用
ICTV_PLANT_INCLUDING_FAMILIES = {
    'Artoviridae': 'plants, invertebrates',
    'Chrysoviridae': 'fungi, plants, invertebrates',
    'Genomoviridae': 'fungi, plants, invertebrates, vertebrates',
    'Kanorauviridae': 'plants, invertebrates, vertebrates',
    'Mitoviridae': 'fungi, plants',
    'Ourmiaviridae': 'invertebrates, plants',
    'Pestiviridae': 'vertebrates, invertebrates, plants',
    'Pseudoviridae': 'protists, fungi, plants, invertebrates',
    'Solemoviridae': 'plants',
    'Spiciviridae': 'plants, invertebrates',
    'Tomosaviridae': 'plants',
    'Virgaviridae': 'plants',
}

_pl_wl_bad_f = sorted(set(NON_PLANT_FAMILIES_FALLBACK) & PLANT_FAMILIES_WHITELIST)
_pl_wl_bad_g = sorted(set(NON_PLANT_GENERA) & PLANT_GENERA_WHITELIST)
if _pl_wl_bad_f or _pl_wl_bad_g:
    raise RuntimeError(
        '宿主名单自相矛盾: 植物白名单与黑名单相交 '
        '科=' + repr(_pl_wl_bad_f) + ' 属=' + repr(_pl_wl_bad_g) +
        ' (历史上出现过 Plant.tsv 零记录即判非植物科的误否决, 见 HOST_REGISTRY_AUDIT_20260915.md)')
