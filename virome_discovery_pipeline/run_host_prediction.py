#!/usr/bin/env python3
"""
run_host_prediction.py — 宿主预测级联工作流 (v2.1, 2026-09-19)

级联方案 (--mode all, 默认):
  L0 Rule:  噬菌体纲 (leviviricetes 等) → Bacteria (查 taxonomy 表, 零成本)
  L0.5 Rule: 非植物纲 (Nucleocytoviricota 各纲/Tokiviricetes/Naldaviricetes/
             Pharingeaviricetes) → 永不 Plant, 按 --sample-clade-map 映射
             Algae/Protist/Unknown (v2.1 新增)
  L1 ICTV:  C9 lookup 全量查官方宿主记录
  L2 RVH:   仅对 L1 未命中的子集跑 RNAVirHost (classify_order + predict)
  L3 PB2:   仅对 L1+L2 双未命中的子集跑 PhaBOX2 CHERRY
  融合:     Rule → ICTV → RVH → PB2 → Unknown; 否决硬返回 (v2.1)

v2.1 加固 (2026-09-19, OneKP 宿主审计后):
  ① Class 级非植物否决 — 纲在 NON_PLANT_CLASSES 即永不 Plant
  ② RVH pred_low_confidence / unclassified 不作为定案依据
  ③ 黑名单否决硬返回 Unknown (Blacklist_Veto), 不再落 RVH/PB2 被旧结果复活
  ④ RVH_result 新鲜度校验 — 早于 --tax 的残留默认拒绝合并 (--allow-stale-rvh 放行)

收益实测 (上次全量 535,103 条模拟对拍):
  RVH 输入 535,103 → ~18k (↓96.6%); CHERRY 输入 → ~14k (↓97.4%)
  host 段 8.75h → ~10min; Final_Host 一致率 99.913% (465 条分歧全部为
  ICTV 官方记录 vs RVH/PB2 模型预测的分歧, 级联取 ICTV 更准)

单工具模式 (--mode ICTV/RNAVirHost/PhaBOX2) 保持原全量行为, 用于补跑分析。

Usage:
  python run_host_prediction.py \
    -i vclust_centroids.fasta \
    --tax taxonomy_out/integrated/final_integrated_classification.tsv \
    -o host_out \
    -t 40 \
    --phabox-db ~/database/virus-db/phabox_db_v2_2
"""

import argparse
import os
import subprocess
import time
import sys
from datetime import datetime
from pathlib import Path
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent

# 噬菌体纲: 分类学规则直接定 Bacteria (决策树 L0)
PHAGE_CLASSES = {'leviviricetes', 'caudoviricetes', 'vidaverviricetes', 'faserviricetes'}

# ---- Class 级非植物否决 (2026-09-19, OneKP 宿主审计后新增) ----
# 这些纲在 ICTV 框架下没有任何植物宿主记录, 实测 (OneKP ensemble) 却被查表/RVH
# 误标为 Plant ~1.7k 条。命中即永不 Plant, 按 --sample-clade-map 映射 Algae/Protist/Unknown:
#   Nucleocytoviricota 各纲 (巨病毒/藻病毒/痘病毒 — 宿主为藻类、原生生物、动物):
#     Megaviricetes (Mimiviridae/Phycodnaviridae/Marseilleviridae/Pithoviridae 等),
#     Mriyaviricetes, Pokkesviricetes
#   Tokiviricetes (Adnaviria, 古菌病毒); Naldaviricetes (杆状病毒类, 节肢动物);
#   Pharingeaviricetes (腺病毒, 脊椎动物)
# 注: 科/属黑名单 (is_blacklisted) 对 Class 级未定论的行不生效, 故需此纲级规则兜底。
NON_PLANT_CLASSES = ['Megaviricetes', 'Mriyaviricetes', 'Pokkesviricetes',
                     'Tokiviricetes', 'Naldaviricetes', 'Pharingeaviricetes']

# 样本 clade 文本的关键词 (小写子串匹配): 藻类 → Algae; 陆生植物 → Protist (表面微生物组)
ALGAL_CLADE_MARKERS = ('chlorophyta', 'chloro', 'rhodophyta', 'rhodo', 'glaucophyta',
                       'prasinophy', 'bacillariophy', 'phaeophyce', 'xanthophy', 'eustigmat',
                       'cryptophy', 'haptophy', 'charophy', 'klebsormid', 'coleochaet',
                       'zygnem', 'dinoflagell', 'algae')
LAND_PLANT_CLADE_MARKERS = ('angiosperm', 'gymnosperm', 'bryophyt', 'lycopod', 'fern',
                            'monilo', 'euphyllo', 'streptophyt')

# --------------- 非植物属黑名单 (Family 级判定门控) ---------------
# 数据源: PlantVirusDB 权威库 classified_clean/Plant.tsv
# 口径: 科级判定无区分力 (同科内植物属与非植物属混杂),
#       故 C9 只报科级时, 以【属】为判据:
#         属在黑名单 → 否决 Plant; 否则放行 (含属在权威库的情形)。
#       属/种级判定不触发黑名单。
NON_PLANT_GENERA = ['Alphabaculovirus', 'Alphambiguivirus', 'Alphaorpheovirus', 'Alphapithovirus', 'Arenavirus', 'Avastrovirus', 'Bandavirus', 'Betambiguivirus', 'Biavirus', 'Bimevirus', 'Catovirus', 'Chaphamavirus', 'Chlorovirus', 'Clandestinovirus', 'Coccolithovirus', 'Criusvirus', 'Crucivirus', 'Cypovirus', 'Delepquintavirus', 'Gammambiguivirus', 'Gammaretrovirus', 'Geminivirus', 'Grisebachstrassevirus', 'Hatfieldvirus', 'Hokovirus', 'Hubsclerovirus', 'Iltovirus', 'Indivirus', 'Inibicvirus', 'Jingmenvirus', 'Jouyvirus', 'Klosneuvirus', 'Kratosvirus', 'Lymphocystivirus', 'Mavirus', 'Megavirus', 'Miraophiovirus', 'Oleurovirus', 'Orthoflavivirus', 'Orthoreovirus', 'Patrovirus', 'Prasinovirus', 'Prymnesiovirus', 'Punuivirus', 'Rheavirus', 'Rimosavirus', 'Simplexvirus', 'Solivirus', 'Sylvanvirus', 'Tethysvirus', 'Theiavirus', 'Timquatrovirus', 'Tobalivirus', 'Tupanvirus', 'Varicellovirus', 'Velanvirus', 'Zetanudivirus', 'Betapartitivirus', 'Botrexvirus', 'Sclerodarnavirus', 'Draflysatellite', 'Unirnavirus', 'Zybavirus', 'Betaendornavirus', 'Cryspovirus', 'Gammapartitivirus']

# 无属信息时的整科兜底黑名单 (权威库中无此科的记录)
NON_PLANT_FAMILIES_FALLBACK = ['Adenoviridae', 'Alloherpesviridae', 'Allomimiviridae', 'Ambiguiviridae', 'Arenaviridae', 'Astroviridae', 'Baculoviridae', 'Deltaflexiviridae', 'Dicistroviridae', 'Gammaflexiviridae', 'Gandrviridae', 'Genomoviridae', 'Hepaciviridae', 'Herelleviridae', 'Hypoviridae', 'Iridoviridae', 'Mamonoviridae', 'Marseilleviridae', 'Maviroviridae', 'Mesomimiviridae', 'Mimiviridae', 'Mitoviridae', 'Nimaviridae', 'Nudiviridae', 'Orpheoviridae', 'Orthoherpesviridae', 'Ourmiaviridae', 'Pandoraviridae', 'Papillomaviridae', 'Peduoviridae', 'Phycodnaviridae', 'Picobirnaviridae', 'Picornaviridae', 'Pithoviridae', 'Poxviridae', 'Pseudoviridae', 'Schizomimiviridae', 'Tobaliviridae']
# ── 2026-09-15 追加: 非植物科登记扩容 (199 科) ──
# 判据: ① 该科在 PlantVirusDB 权威植物库 classified_clean/Plant.tsv 中零记录
#         (Plant.tsv md5 6e3200b1b2a769610f0780979a8dec26; 201,488 条植物病毒记录 / 解析 46 科 / 821 属)
#       ② 该科在枸杞(goji)或 OneKP 树 05 分类表 / 下游植物表中出现
# 逐科行数与代表属证据: scripts/audit/non_plant_family_registry.tsv
# 宿主域逐个核实结论: scripts/audit/TAXONOMY_05_CALIB_STATE.md
# 作用面: 仅在 C9 只给到科级且属缺失时生效; 属/种级判定不经过本名单
#         (属级短路见 is_blacklisted 的 is_trusted_level 分支)
NON_PLANT_FAMILIES_FALLBACK += [
    # 出现在下游植物病毒结果表中的 8 科 (承重条目, 逐科核定)
    'Metaviridae',           # 05:g133/o5333 植物表:g6/o102 | Metavirus(5253), Errantivirus(143), Badnavirus(26)
    'Draupnirviridae',       # 05:g6/o12 植物表:g4/o2 | Malificivirus(8), Torentivirus(3), Volantivirus(2)
    'Iflaviridae',           # 05:g25/o990 植物表:g4/o0 | Iflavirus(884)
    'Barnaviridae',          # 05:g5/o45 植物表:g3/o0 | Barnavirus(8), Enamovirus(2), Sobemovirus(1)
    'Discoviridae',          # 05:g17/o2 植物表:g3/o0 | Orthodiscovirus(22)
    'Hydriviridae',          # 05:g236/o0 植物表:g1/o0 | Alphahydrivirus(237)
    'Kanorauviridae',        # 05:g6/o64 植物表:g1/o0 | Shabtayvirus(50), Ninurtavirus(8), Shanivirus(4)
    'Narnaviridae',          # 05:g67/o1069 植物表:g1/o0 | Narnavirus(793), Ourmiavirus(1)
    # 仅出现在 05 分类表中的 191 科 (按字母序, 不参与当前下游表清理)
    'Ackermannviridae',
    'Adamaviridae',
    'Adomaviridae',
    'Ahmunviridae',
    'Aliceevansviridae',
    'Aliusviridae',
    'Alphaormycoviridae',
    'Alphatetraviridae',
    'Alternaviridae',
    'Amnoonviridae',
    'Anamaviridae',
    'Andersonviridae',
    'Anelloviridae',
    'Aoguangviridae',
    'Apasviridae',
    'Arenbergviridae',
    'Arteriviridae',
    'Artiviridae',
    'Artoviridae',
    'Ascoviridae',
    'Asfarviridae',
    'Assiduviridae',
    'Autographiviridae',
    'Autonotataviridae',
    'Autoscriptoviridae',
    'Autosignataviridae',
    'Autotranscriptaviridae',
    'Bacilladnaviridae',
    'Belpaoviridae',
    'Berryhillviridae',
    'Betaormycoviridae',
    'Bidnaviridae',
    'Birnaviridae',
    'Blumeviridae',
    'Bornaviridae',
    'Botybirnaviridae',
    'Buchnerviridae',
    'Caliciviridae',
    'Carmotetraviridae',
    'Casidaviridae',
    'Casjensviridae',
    'Chaseviridae',
    'Chimalliviridae',
    'Chrysoviridae',
    'Chuviridae',
    'Circoviridae',
    'Clermontviridae',
    'Colingsworthviridae',
    'Connertonviridae',
    'Coronaviridae',
    'Cremegaviridae',
    'Crevaviridae',
    'Cruciviridae',
    'Cruliviridae',
    'Curvulaviridae',
    'Darmviridae',
    'Deltaormycoviridae',
    'Demerecviridae',
    'Dishuiviroviridae',
    'Drexlerviridae',
    'Dumbiviridae',
    'Duneviridae',
    'Ehrlichviridae',
    'Endolinaviridae',
    'Epsomviridae',
    'Eupolintoviridae',
    'Euroniviridae',
    'Fervensviridae',
    'Fiersviridae',
    'Filamentoviridae',
    'Filoviridae',
    'Fusariviridae',
    'Fuselloviridae',
    'Gammaormycoviridae',
    'Giardiaviridae',
    'Graaviviridae',
    'Grimontviridae',
    'Hadakaviridae',
    'Hafunaviridae',
    'Haloferuviridae',
    'Halomagnusviridae',
    'Halspiviridae',
    'Hantaviridae',
    'Helgolandviridae',
    'Hepadnaviridae',
    'Hepeviridae',
    'Hodgkinviridae',
    'Hytrosaviridae',
    'Inoviridae',
    'Inseviridae',
    'Intestiviridae',
    'Itzamnaviridae',
    'Jeanschmidtviridae',
    'Kirkoviridae',
    'Kleczkowskaviridae',
    'Kolmioviridae',
    'Konodaiviridae',
    'Krittikaviridae',
    'Kyanoviridae',
    'Lakviridae',
    'Lebotiviridae',
    'Leishbuviridae',
    'Lindbergviridae',
    'Lipothrixviridae',
    'Lispiviridae',
    'Ludisviridae',
    'Malacoherpesviridae',
    'Marnaviridae',
    'Matonaviridae',
    'Medioniviridae',
    'Megabirnaviridae',
    'Megatotiviridae',
    'Mesoniviridae',
    'Mesyanzhinovviridae',
    'Microviridae',
    'Monocitiviridae',
    'Mononiviridae',
    'Mtkvariviridae',
    'Mycoalphaviridae',
    'Mymonaviridae',
    'Nairoviridae',
    'Nanghoshaviridae',
    'Nanhypoviridae',
    'Naryaviridae',
    'Natareviridae',
    'Nenyaviridae',
    'Nipumfusiviridae',
    'Nodaviridae',
    'Noraviridae',
    'Nyamiviridae',
    'Oomyviridae',
    'Ootiviridae',
    'Orthomyxoviridae',
    'Ouroboviridae',
    'Pachyviridae',
    'Parahypoviridae',
    'Paramyxoviridae',
    'Parvoviridae',
    'Pecoviridae',
    'Permutotetraviridae',
    'Pestiviridae',
    'Phasmaviridae',
    'Phlegiviridae',
    'Phypoliviridae',
    'Pneumoviridae',
    'Polycipiviridae',
    'Polymycoviridae',
    'Polyomaviridae',
    'Pootjesviridae',
    'Pseudototiviridae',
    'Pyrstoviridae',
    'Quadriviridae',
    'Redondoviridae',
    'Rhizouliviridae',
    'Roniviridae',
    'Rountreeviridae',
    'Salasmaviridae',
    'Sarkviridae',
    'Schitoviridae',
    'Shortaselviridae',
    'Sinhaliviridae',
    'Smacoviridae',
    'Solinviviridae',
    'Spiciviridae',
    'Splipalmiviridae',
    'Sputniviroviridae',
    'Stackebrandtviridae',
    'Stanwilliamsviridae',
    'Steigviridae',
    'Steitzviridae',
    'Straboviridae',
    'Sunviridae',
    'Suoliviridae',
    'Tectiviridae',
    'Tobaniviridae',
    'Toyamaviridae',
    'Trimbiviridae',
    'Tulasviridae',
    'Vandenendeviridae',
    'Vertoviridae',
    'Vilmaviridae',
    'Winoviridae',
    'Wupedeviridae',
    'Xinmoviridae',
    'Yadokariviridae',
    'Yangangviridae',
    'Yaraviridae',
    'Yueviridae',
    'Zierdtviridae',
    'Zimmerviridae',
    'Zobellviridae',
]

# ── 2026-09-15 复核撒回: 判据 "Plant.tsv 零记录" 在改名/别名与宿主列含植物两种情形下失效 ──
#    判定源: ICTV VMR MSL41 (VMR_MSL41.v1.20260729) + 官方 "Virus Properties by Family"
#            (全表 427 科, Host 列原文逐科核对): scripts/audit/ICTV_HOST_SWEEP_233.tsv
#    A. ICTV 现行宿主列含 plants (不得标为非植物科):
#       Ourmiaviridae  植物 + 无脊椎动物 (含植物病毒属 Ourmiavirus; 2025.013F 新建科)
#       Mitoviridae    真菌 + 植物 (植物线粒体 mitovirus 有独立侵染证据)
#       Artoviridae    plants, invertebrates
#       Chrysoviridae  fungi, plants, invertebrates (成员为真菌病毒)
#       Genomoviridae  fungi, plants, invertebrates, vertebrates (CRESS-DNA)
#       Pestiviridae   vertebrates, invertebrates, plants
#       Spiciviridae   plants, invertebrates
#    B. 宿主列含 plants, 但成员本质为内源反转座子 (语境非植物侵染), 不宜用 "非植物科" 表述:
#       Pseudoviridae (Ty1/copia 型 LTR 反转录元件) / Metaviridae (Ty3/gypsy 型)
#       两科更适用 "宿主基因组内源元件" 单列处置, 见 HOST_REGISTRY_AUDIT_20260915.md 待办 4
#    C. ICTV 表列 plants, ViralZone 记宿主未知 (污水/粪便), 两处冲突: Kanorauviridae (撒回待核)
#    撒回方式=显式移除, 原始定义行与追加块均保持零改动; 幂等, 可重复导入。
NON_PLANT_FAMILIES_UNDER_REVIEW = [
    # A. 宿主列含植物
    'Ourmiaviridae', 'Mitoviridae', 'Artoviridae', 'Chrysoviridae', 'Genomoviridae',
    'Pestiviridae', 'Spiciviridae',
    # B. 内源反转座子语境
    'Pseudoviridae', 'Metaviridae',
    # C. 宿主记录冲突待核
    'Kanorauviridae',
]
for _fam in NON_PLANT_FAMILIES_UNDER_REVIEW:
    if _fam in NON_PLANT_FAMILIES_FALLBACK:
        NON_PLANT_FAMILIES_FALLBACK.remove(_fam)

# ── 2026-09-15 复核撒回: 植物库其实认得这些属, 属级"零记录"判据失效 ──
#    Betapartitivirus 植物库属级记录 79 条 + 物种名 8 条 (Betapartitivirus primulae / trifolii)
#    Geminivirus       植物库物种名 1 条 (Geminivirus isolate Euphua Iguala); 该名为
#                      旧统称/非 ICTV 属, Geminiviridae 在植物库有 30001 条记录
#    实测影响: goji 树 0 行受影响 (撤回后 Plant 仍为 1028), 纯粹降低误否决风险
NON_PLANT_GENERA_UNDER_REVIEW = ['Betapartitivirus', 'Geminivirus']
for _gen in NON_PLANT_GENERA_UNDER_REVIEW:
    if _gen in NON_PLANT_GENERA:
        NON_PLANT_GENERA.remove(_gen)

# 科级否决开关: True=科在非植物科名单即否决 Plant (科级优先);
#               False=仅当属缺失时才查科名单 (旧口径, 属/种级判定直接放行)
FAMILY_FIRST_VETO = True

TRUSTED_LEVELS = {'genus', 'species', 'species*'}


def is_trusted_level(det_level):
    """Determination_Level 是否为属/种级 (C9 有区分力的层级)。"""
    if det_level is None:
        return False
    d = str(det_level).split('(')[0].strip().lower()
    return d in TRUSTED_LEVELS


def is_blacklisted(family, genus, det_level):
    """科级判定时按属黑名单否决 Plant; 无属信息时按科兜底。

    FAMILY_FIRST_VETO=True 时科名单优先: 科已判定为非植物科即否决 Plant,
    不受属/种级判定层影响 (否则同一 (科,属) 组合在属缺失时被否决、属冲突时反而不否决,
    逻辑不自洽: 若 Barnaviridae 非植物, 则 Barnaviridae/Sobemovirus 也不是植物病毒)。
    """
    fam = str(family).strip() if family is not None else ''
    if FAMILY_FIRST_VETO and fam and fam in NON_PLANT_FAMILIES_FALLBACK:
        return True
    if is_trusted_level(det_level):
        return False
    gen = str(genus).strip() if genus is not None else ''
    if gen and gen not in ('NA', 'nan'):
        return gen in NON_PLANT_GENERA
    # 属信息缺失: 整科兜底
    return bool(fam and fam in NON_PLANT_FAMILIES_FALLBACK)


def run_step(cmd, step_name):
    print(f"\n{'='*75}\n[{datetime.now().strftime('%H:%M:%S')}] {step_name}")
    head = ' '.join(cmd) if isinstance(cmd, list) else cmd
    print(f"[CMD] {head[:200]}\n{'='*75}")
    t0 = time.time()
    try:
        result = subprocess.run(cmd, shell=True if isinstance(cmd, str) else False,
                                check=True, capture_output=True, text=True)
        print(f"[OK] {step_name} — {time.time() - t0:.0f}s")
        return True
    except subprocess.CalledProcessError as e:
        print(f"[FAIL] {step_name} — {time.time() - t0:.0f}s (rc={e.returncode})")
        if e.stderr:
            for l in e.stderr.strip().split('\n')[-6:]: print(f"  [stderr] {l}")
        return False

def check_file(path, min_size=10):
    return Path(path).is_file() and Path(path).stat().st_size > min_size

def read_fasta_ids(fa):
    ids = set()
    with open(fa) as f:
        for line in f:
            if line.startswith('>'):
                ids.add(line[1:].split()[0])
    return ids

def extract_fasta_subset(input_fa, ids, out_fa):
    """只保留 header first-token 在 ids 集合中的序列, 返回条数"""
    n = 0
    with open(input_fa) as fin, open(out_fa, 'w') as fout:
        keep = False
        for line in fin:
            if line.startswith('>'):
                keep = line[1:].split()[0] in ids
                if keep: n += 1
            if keep:
                fout.write(line)
    return n

# ----------------- ICTV / RVH 命中判定 (与融合树同一套 normalize 语义) -----------------

RVH_MAP = {'viridiplantae': 'Plant', 'fungi': 'Fungi',
           'chordata': 'Animal', 'invertebrate': 'Animal', 'metazoa': 'Animal',
           'animal': 'Animal', 'bacteria': 'Bacteria'}

def normalize_c9(h):
    h = str(h).strip() if not pd.isna(h) else 'nan'
    if h in ('', 'Unknown', 'None', 'NA', 'nan'): return 'Unknown'
    if h in ('Insecta', 'Arachnida', 'Aves', 'Human', 'Animal_other'): return 'Animal'
    if h == 'Oomycetes': return 'Protist'
    return h

def read_ictv_hits(c9_tsv):
    """C9 lookup 结果 → 命中 ID 集合 (与 normalize_c9 同语义)"""
    c9 = pd.read_csv(c9_tsv, sep='\t')
    return {str(cid) for cid, h in zip(c9['contig_id'], c9['Predicted_Host'])
            if normalize_c9(h) != 'Unknown'}

def read_rvh_hits(rvh_csv):
    """RNAVirHost result.csv → 命中 ID 集合 (与 parse_rvh 同语义: 映射表内才算命中)"""
    rvh = pd.read_csv(rvh_csv)
    id_col = rvh.columns[0]
    if 'Unnamed' in id_col or 'y|virus order' in rvh.columns:
        rvh = rvh.rename(columns={id_col: 'contig_id'})
    hits = set()
    for cid, h, e in zip(rvh['contig_id'],
                         rvh.get('pred|L1', pd.Series(dtype=str)),
                         rvh.get('evidence', pd.Series(dtype=str))):
        hh = str(h).strip().lower()
        if hh in RVH_MAP and str(e).strip().lower() != 'unclassified':
            hits.add(str(cid))
    return hits

# ----------------- 数据准备与工具执行 (级联) -----------------

def run_tools(args):
    outdir = args.output_dir
    rvh_csv = os.path.join(outdir, "RVH_result", "result.csv")
    pb2_tsv = os.path.join(outdir, "phabox2_output", "final_prediction", "cherry_prediction.tsv")
    c9_tsv = os.path.join(outdir, "C9_ICTV_result", "classification_result.tsv")
    rvh_taxa = os.path.join(outdir, "RVH_taxa.csv")
    sub1_fa = os.path.join(outdir, "cascade_sub1_ictv_miss.fasta")   # L1 未命中 → RVH 输入
    sub2_fa = os.path.join(outdir, "cascade_sub2_rvh_miss.fasta")   # L1+L2 双未命中 → CHERRY 输入

    # ── L1: ICTV C9 Lookup (全量, 查表, 先行) ──
    ictv_ready = check_file(c9_tsv, 50)
    if not args.skip_ictv:
        if not args.force and ictv_ready:
            print(f"[SKIP] ICTV Lookup (C9) — exists: {c9_tsv}")
        elif check_file(args.tax):
            host_dir = os.path.join(outdir, "C9_ICTV_result")
            os.makedirs(host_dir, exist_ok=True)
            cmd = f"python {SCRIPT_DIR}/utils/classify_contigs.py -i {args.tax} -f {args.input} --output_dir {host_dir} --prob_dir {args.prob_dir} --mode high"
            if run_step(cmd, "ICTV Taxonomy Lookup (C9)"):
                ictv_ready = True
        else:
            print("[WARN] Valid taxonomy TSV required for C9 ICTV lookup.")

    # 单工具模式: 保持原全量行为
    if args.mode != "all":
        if args.mode == "RNAVirHost" and not args.skip_rnavirhost:
            if not args.force and check_file(rvh_csv, 50):
                print(f"[SKIP] RNAVirHost — exists: {rvh_csv}")
            else:
                if os.path.isdir(os.path.join(outdir, "RVH_result")):
                    import shutil; shutil.rmtree(os.path.join(outdir, "RVH_result"))
                if not check_file(rvh_taxa, 50):
                    run_step(["rnavirhost", "classify_order", "-i", args.input, "-o", rvh_taxa],
                             "RNAVirHost Classify Order")
                cmd2 = ["rnavirhost", "predict", "-i", args.input, "--taxa", rvh_taxa, "-o", os.path.join(outdir, "RVH_result")]
                run_step(cmd2, "RNAVirHost Predict")
        elif args.mode == "PhaBOX2" and not args.skip_phabox:
            if not args.force and check_file(pb2_tsv, 50):
                print(f"[SKIP] PhaBOX2 — exists: {pb2_tsv}")
            else:
                cmd = ["phabox2", "--task", "cherry", "--dbdir", args.phabox_db,
                       "--outpth", os.path.join(outdir, "phabox2_output"),
                       "--contigs", args.input, "--threads", str(args.threads), "--len", "500"]
                run_step(cmd, "PhaBOX2 CHERRY Host Prediction")
                _cleanup_phabox(outdir)
        return rvh_csv, pb2_tsv, c9_tsv

    # ── 级联模式 (--mode all) ──
    if not ictv_ready:
        print("[FATAL] 级联模式需要 ICTV lookup 先行 (C9_ICTV_result 缺失且无法生成)")
        sys.exit(1)
    ictv_hits = read_ictv_hits(c9_tsv)
    print(f"[Cascade L1] ICTV 命中 {len(ictv_hits)} 条")

    # L0 Rule 命中 (噬菌体纲, 查 taxonomy 表) — 无需 RVH/CHERRY
    rule_ids = set()
    if check_file(args.tax):
        tdf = pd.read_csv(args.tax, sep='\t', usecols=lambda c: c in ('contig_id', 'Class'))
        if 'Class' in tdf.columns:
            rule_ids = {str(c) for c, cl in zip(tdf['contig_id'], tdf['Class'])
                        if str(cl).strip().lower() in PHAGE_CLASSES}
    print(f"[Cascade L0] Rule (噬菌体纲) 命中 {len(rule_ids)} 条")

    all_ids = read_fasta_ids(args.input)
    # 子集只取 tax 表内的 contigs: tax 表外 (~75k) 无 Class 无 ICTV 记录,
    # ensemble 本就不覆盖 (left join 以 tax 为主表), 与旧版行为一致
    tax_ids = set()
    if check_file(args.tax):
        tdf = pd.read_csv(args.tax, sep='\t', usecols=lambda c: c in ('contig_id', 'Class'))
        tax_ids = {str(c) for c in tdf['contig_id']}
        if 'Class' in tdf.columns:
            rule_ids = {str(c) for c, cl in zip(tdf['contig_id'], tdf['Class'])
                        if str(cl).strip().lower() in PHAGE_CLASSES}
    n_tax_miss = len(tax_ids) - len(ictv_hits & tax_ids) - len(rule_ids - ictv_hits)
    # L2 RVH 输入子集: tax 表内 - ICTV 命中 - Rule 命中
    sub1_ids = tax_ids - ictv_hits - rule_ids
    n1 = extract_fasta_subset(args.input, sub1_ids, sub1_fa)
    print(f"[Cascade L2] RVH 子集: {n1} 条 (tax 内 {n_tax_miss} 未命中, 全库 {len(all_ids)}) → {sub1_fa}")

    if not args.skip_rnavirhost:
        if not args.force and check_file(rvh_csv, 50):
            print(f"[SKIP] RNAVirHost — exists: {rvh_csv}")
        else:
            if os.path.isdir(os.path.join(outdir, "RVH_result")):
                import shutil; shutil.rmtree(os.path.join(outdir, "RVH_result"))
            if not check_file(rvh_taxa, 50):
                run_step(["rnavirhost", "classify_order", "-i", sub1_fa, "-o", rvh_taxa],
                         "RNAVirHost Classify Order (cascade sub1)")
            cmd2 = ["rnavirhost", "predict", "-i", sub1_fa, "--taxa", rvh_taxa, "-o", os.path.join(outdir, "RVH_result")]
            run_step(cmd2, "RNAVirHost Predict (cascade sub1)")

    rvh_hits = read_rvh_hits(rvh_csv) if check_file(rvh_csv, 50) else set()
    print(f"[Cascade L2] RVH 命中 {len(rvh_hits & sub1_ids)} 条 (子集内)")

    # L3 CHERRY 子集: sub1 中 RVH 也未命中
    sub2_ids = sub1_ids - rvh_hits
    n2 = extract_fasta_subset(args.input, sub2_ids, sub2_fa)
    print(f"[Cascade L3] CHERRY 子集: {n2} 条 (↓{100 - 100.0 * n2 / max(len(tax_ids) or 1, 1):.1f}%) → {sub2_fa}")

    if not args.skip_phabox and n2 > 0:
        if not args.force and check_file(pb2_tsv, 50):
            print(f"[SKIP] PhaBOX2 — exists: {pb2_tsv}")
        else:
            cmd = ["phabox2", "--task", "cherry", "--dbdir", args.phabox_db,
                   "--outpth", os.path.join(outdir, "phabox2_output"),
                   "--contigs", sub2_fa, "--threads", str(args.threads), "--len", "500"]
            run_step(cmd, "PhaBOX2 CHERRY Host Prediction (cascade sub2)")
            _cleanup_phabox(outdir)
    elif n2 == 0:
        print("[Cascade L3] 子集为空, 跳过 CHERRY")

    return rvh_csv, pb2_tsv, c9_tsv

def _cleanup_phabox(outdir):
    # 清理 PhaBOX2 中间产物 (全对全蛋白聚类 ~179G)
    for junk in ["midfolder", "filtered_contigs.fa"]:
        p = os.path.join(outdir, "phabox2_output", junk)
        if os.path.exists(p):
            try:
                import shutil
                if os.path.isdir(p): shutil.rmtree(p)
                else: os.remove(p)
                print(f"  清理 PhaBOX2/{junk}")
            except Exception as e:
                print(f"  清理 PhaBOX2/{junk} 失败: {e}")

# ----------------- 决策树与整合逻辑 -----------------

def normalize_c9(h):
    if pd.isna(h) or h in ['Unknown', 'None']: return 'Unknown'
    h = str(h)
    if h in ['Insecta', 'Arachnida', 'Aves', 'Human', 'Animal_other']: return 'Animal'
    if h in ['Oomycetes']: return 'Protist'
    return h

# RVH evidence 中不得作为定案依据的标签 (2026-09-19):
#   unclassified = 工具自认无判断; pred_low_confidence = k-mer 低置信
#   (实测 OneKP: Adnaviria 古菌病毒 30 条凭它被判成 Plant)。
# read_rvh_hits (子集划分) 不受此门槛影响 — 低置信仍算"已作答", 不再送 PB2。
RVH_NON_DECISIVE_EVIDENCE = {'pred_low_confidence', 'unclassified'}

def parse_rvh(row):
    h = str(row.get('pred|L1', 'Unknown')).strip().lower()
    ev = str(row.get('evidence', '')).strip().lower()
    if h == 'unknown' or ev in RVH_NON_DECISIVE_EVIDENCE: return 'Unknown'
    return RVH_MAP.get(h, 'Unknown')

# ---------- Class 级非植物否决的宿主映射 (2026-09-19) ----------
DIRECT_HOST_LABELS = {'algae': 'Algae', 'protist': 'Protist', 'unknown': 'Unknown',
                      'plant': 'Plant', 'animal': 'Animal', 'fungi': 'Fungi',
                      'bacteria': 'Bacteria'}

def load_sample_clade_map(path):
    """TSV 两列: 样本/run id → clade 文本 (或直接给 Algae/Protist/Unknown 等终值)"""
    m = {}
    with open(path, encoding='utf-8') as f:
        for line in f:
            c = line.rstrip('\n').split('\t')
            if len(c) >= 2 and c[0].strip():
                m[c[0].strip()] = c[1].strip()
    return m

def class_veto_host(row):
    """Class ∈ NON_PLANT_CLASSES 时的宿主去向:
    样本 clade 命中藻类关键词 → Algae; 陆生植物关键词 → Protist (表面微生物组背景);
    clade 列直接给了终值标签 (Algae/Protist/...) 则原样采用; 无样本信息 → Unknown (宁缺勿错)。"""
    clade = str(row.get('sample_clade', '') or '').strip()
    if not clade:
        return 'Unknown'
    direct = DIRECT_HOST_LABELS.get(clade.lower())
    if direct:
        return direct
    low = clade.lower()
    if any(k in low for k in ALGAL_CLADE_MARKERS):
        return 'Algae'
    if any(k in low for k in LAND_PLANT_CLADE_MARKERS):
        return 'Protist'
    return 'Unknown'

def _rvh_freshness(rvh_csv, tax_path):
    """返回 'missing' | 'fresh' | 'stale':
    result.csv 的 mtime 早于 --tax (05 分类产物) = 上一轮运行残留 (级联时代的旧全量结果)。"""
    if not os.path.exists(rvh_csv): return 'missing'
    try:
        return 'fresh' if os.path.getmtime(rvh_csv) >= os.path.getmtime(tax_path) else 'stale'
    except OSError:
        return 'fresh'

def parse_pb2(row):
    lineage = (str(row.get('Host_NCBI_lineage', '')) + "|" +
               str(row.get('Host_GTDB_lineage', '')) + "|" +
               str(row.get('Host', ''))).lower()
    if 'bacteria' in lineage: return 'Bacteria'
    if 'archaea' in lineage: return 'Archaea'
    if 'streptophyta' in lineage or 'viridiplantae' in lineage or 'plant' in lineage: return 'Plant'
    if 'fungi' in lineage: return 'Fungi'
    if any(k in lineage for k in ['metazoa', 'animal', 'chordata', 'arthropoda', 'insecta']): return 'Animal'
    return 'Unknown'

def decision_tree_cascade(row):
    """级联融合: Rule(噬菌体纲) → Rule(非植物纲) → ICTV → RVH → PB2 → Unknown
    子集外的 RVH/PB2 列为 NaN → Unknown → 级联自然跳过。
    规则/黑名单否决一律硬返回, 不再落到下游工具
    (否则 ICTV-Plant 被否决后会由旧全量 RVH 的 k-mer 预测复活成 Plant, OneKP 实测路径)。"""
    tax_class = str(row.get('Class', '')).strip()
    if tax_class.lower() in PHAGE_CLASSES:
        return 'Bacteria', 'Rule_Class_Taxonomy'
    if tax_class in NON_PLANT_CLASSES:
        return class_veto_host(row), 'Rule_Class_NonPlant'
    h_ictv = normalize_c9(row.get('Host_ICTV', 'Unknown'))
    if h_ictv != 'Unknown':
        if h_ictv == 'Plant' and is_blacklisted(
                row.get('Family'), row.get('Genus'),
                row.get('Determination_Level')):
            return 'Unknown', 'Blacklist_Veto'
        return h_ictv, 'ICTV_Preferred'
    h_rvh = parse_rvh(row)
    if h_rvh != 'Unknown':
        return h_rvh, 'RVH_Preferred'
    h_pb2 = parse_pb2(row)
    if h_pb2 != 'Unknown':
        return h_pb2, 'PB2_Preferred'
    return 'Unknown', 'Unassigned'

def decision_tree(row):
    """旧版并行融合树 (三源全量投票), 保留用于对照与回退。"""
    tax_class = str(row.get('Class', '')).strip()
    if tax_class.lower() in ['leviviricetes', 'caudoviricetes', 'vidaverviricetes', 'faserviricetes']:
        return 'Bacteria', 'Rule_Class_Taxonomy'
    if tax_class in NON_PLANT_CLASSES:
        return class_veto_host(row), 'Rule_Class_NonPlant'
    h_ictv = normalize_c9(row.get('Host_ICTV', 'Unknown'))
    if h_ictv == 'Plant' and is_blacklisted(
            row.get('Family'), row.get('Genus'),
            row.get('Determination_Level')):
        h_ictv = 'Unknown'
    h_rvh  = parse_rvh(row)
    h_pb2  = parse_pb2(row)
    if h_ictv == 'Unknown' and h_rvh == 'Unknown' and h_pb2 == 'Unknown':
        return 'Unknown', 'Unassigned'
    if h_ictv == h_rvh and h_ictv != 'Unknown':
        return h_ictv, 'ICTV_RVH_Agree'
    if h_pb2 != 'Unknown':
        if h_pb2 == h_ictv: return h_ictv, 'ICTV_PB2_Agree_Tiebreaker'
        if h_pb2 == h_rvh:  return h_rvh,  'RVH_PB2_Agree_Tiebreaker'
    if h_ictv != 'Unknown': return h_ictv, 'ICTV_Preferred'
    if h_rvh != 'Unknown': return h_rvh,  'RVH_Preferred'
    if h_pb2  != 'Unknown': return h_pb2,  'PB2_Preferred'
    return 'Unknown', 'Fallback'

def run_ensemble(args, rvh_csv, pb2_tsv, c9_tsv):
    print(f"\n{'='*75}\n[Ensemble] {'Cascade' if args.mode == 'all' else 'Parallel'} Decision Tree Integration\n{'='*75}")
    df = pd.read_csv(args.tax, sep='\t') if check_file(args.tax) else pd.DataFrame(columns=['contig_id', 'Class'])

    if check_file(c9_tsv):
        _c9_cols = ['contig_id', 'Predicted_Host']
        # 只取 Determination_Level: C9 的 Family/Genus 与 args.tax 同源 (classify_contigs.py
        # 以 args.tax 为输入直通), 实测 20892/20892 行完全一致, 无需重复携带;
        # 若一并携带会与 df 同名, merge 后变成 Family_x/Family_y,
        # is_blacklisted(row.get('Family'), ...) 取到 None, 名单恒不生效 (实测命中 0/20892)。
        _c9_extra = [c for c in ('Determination_Level',)
                     if c in pd.read_csv(c9_tsv, sep='\t', nrows=0).columns]
        _c9_cols += _c9_extra
        c9 = pd.read_csv(c9_tsv, sep='\t')[_c9_cols].rename(columns={'Predicted_Host': 'Host_ICTV'})
        df = df.merge(c9, on='contig_id', how='left')

    # RVH_result 新鲜度校验 (2026-09-19): 早于 --tax 的 result.csv 是上一轮运行的
    # 旧全量残留 — 级联下若照常合并, ICTV 未判出的子集会被旧 k-mer 预测接手
    # (OneKP 实测: 9月4日全量 result.csv 存活到级联时代, Adnaviria 30 条低置信 Plant 由此复活)。
    # --allow-stale-rvh 仅限对拍/审计时显式放行。
    rvh_state = _rvh_freshness(rvh_csv, args.tax)
    if rvh_state == 'stale':
        if getattr(args, 'allow_stale_rvh', False):
            print("  [WARN] --allow-stale-rvh: 合并旧 RVH 结果 (仅限对拍/审计)")
        else:
            print(f"  [WARN] RVH 结果早于 taxonomy 输出 (上一轮运行残留), 已拒绝合并: {rvh_csv}")
            print("         如确需对拍旧结果, 请加 --allow-stale-rvh")
    if check_file(rvh_csv) and rvh_state != 'stale':
        rvh = pd.read_csv(rvh_csv)
        id_col = rvh.columns[0]
        if 'Unnamed' in id_col or 'y|virus order' in rvh.columns:
            rvh.rename(columns={id_col: 'contig_id'}, inplace=True)
            use_cols = ['contig_id'] + [c for c in ['pred|L1', 'pred|L2', 'evidence'] if c in rvh.columns]
            df = df.merge(rvh[use_cols], on='contig_id', how='left')

    # Class 级否决的宿主映射: contig 内嵌 run/样本 id → 样本 clade (可选)
    if getattr(args, 'sample_clade_map', None):
        import re as _re
        clade_map = load_sample_clade_map(args.sample_clade_map)
        def _clade_of(cid):
            cid = str(cid)
            if cid in clade_map: return clade_map[cid]
            m = _re.match(r'([A-Za-z]{2}\d+)', cid)
            if m and m.group(1) in clade_map: return clade_map[m.group(1)]
            return ''
        df['sample_clade'] = df['contig_id'].map(_clade_of)

    if check_file(pb2_tsv):
        pb2 = pd.read_csv(pb2_tsv, sep='\t').rename(columns={'Accession': 'contig_id'})
        use_cols = [c for c in ['contig_id', 'Host', 'Host_NCBI_lineage', 'Host_GTDB_lineage'] if c in pb2.columns]
        df = df.merge(pb2[use_cols], on='contig_id', how='left')

    tree = decision_tree_cascade if args.mode == 'all' else decision_tree
    results = df.apply(tree, axis=1)
    df['Final_Host'], df['Decision_Method'] = zip(*results)

    out_tsv = os.path.join(args.output_dir, "ensemble_host_summary.tsv")
    df.to_csv(out_tsv, sep='\t', index=False)

    print("\n  [Summary]")
    for host, cnt in df['Final_Host'].value_counts().items():
        print(f"    {host:<15s}: {cnt:>6} contigs")
    print("\n  [Decision Methods]")
    for m, cnt in df['Decision_Method'].value_counts().items():
        print(f"    {m:<28s}: {cnt:>6}")

    print("\n  [Splitting FastA]")
    host_dict = dict(zip(df['contig_id'].astype(str), df['Final_Host']))
    import re as _re
    _acc_to_cid = {}
    for _cid in host_dict:
        _acc_match = _re.search(r'[A-Z]{2}\d{6}\.\d+|[A-Z]{1,2}_\d+\.\d+', str(_cid))
        if _acc_match:
            _acc = _acc_match.group(0)
            if _acc not in _acc_to_cid: _acc_to_cid[_acc] = _cid
        _src_match = _re.search(r'(?:source|src)=(\S+)', str(_cid))
        if _src_match:
            _src = _src_match.group(1).rstrip('|')
            if _src not in _acc_to_cid: _acc_to_cid[_src] = _cid

    def _resolve_host(fasta_id):
        if fasta_id in host_dict:
            return host_dict[fasta_id]
        for _pat in [r'[A-Z]{2}\d{6}\.\d+', r'[A-Z]{1,2}_\d+\.\d+',
                     r'(?:source|src)=(\S+)']:
            _m = _re.search(_pat, fasta_id)
            if _m:
                _key = _m.group(1) if '=' in _pat else _m.group(0)
                _key = _key.rstrip('|')
                if _key in _acc_to_cid:
                    return host_dict[_acc_to_cid[_key]]
                if _key in host_dict:
                    return host_dict[_key]
        return 'Unknown'

    out_fastas = {}
    with open(args.input, 'r') as f:
        curr_id, curr_seq = None, []
        for line in f:
            if line.startswith('>'):
                if curr_id:
                    host = _resolve_host(curr_id)
                    out_fastas.setdefault(host, []).append(f">{curr_id}\n{''.join(curr_seq)}")
                curr_id = line[1:].split()[0]
                curr_seq = []
            else: curr_seq.append(line)
        if curr_id:
            host = _resolve_host(curr_id)
            out_fastas.setdefault(host, []).append(f">{curr_id}\n{''.join(curr_seq)}")

    fasta_dir = os.path.join(args.output_dir, "host_classified_fasta")
    os.makedirs(fasta_dir, exist_ok=True)
    for host, seqs in out_fastas.items():
        with open(os.path.join(fasta_dir, f"{host}.classified.fasta"), 'w') as f: f.write("".join(seqs))
    print(f"    Saved into: {fasta_dir}/\n")

def main():
    p = argparse.ArgumentParser(description="Cascade Host Prediction Workflow (ICTV > RVH > PB2)")
    p.add_argument("-i", "--input", required=True, help="Input FASTA file")
    p.add_argument("--tax", required=True, help="Integrated taxonomy TSV (from 01_run_taxonomy.py)")
    p.add_argument("-o", "--output-dir", default="host_out", help="Output directory")
    p.add_argument("-t", "--threads", type=int, default=40, help="Threads to use")
    p.add_argument("--phabox-db", default=os.path.expanduser("~/database/virus-db/phabox_db_v2_2"), help="PhaBOX2 DB")
    p.add_argument("--prob-dir", default=str(SCRIPT_DIR.parent / "database" / "cross_analysis"),
                   help="ICTV 宿主概率表目录 (classify_contigs 使用)")
    p.add_argument("--sample-clade-map", default=None,
                   help="TSV 两列 (样本/run id → clade 文本或 Algae/Protist 等终值), "
                        "供 Class 级非植物否决映射宿主; 不给则此类判 Unknown")
    p.add_argument("--allow-stale-rvh", action="store_true",
                   help="允许合并早于 --tax 的旧 RVH_result (仅对拍/审计)")

    p.add_argument("--mode", default="all",
                   choices=["all", "ICTV", "RNAVirHost", "PhaBOX2"],
                   help="all=级联 (ICTV 先行, RVH/PB2 只跑未命中子集); 单工具=全量补跑")
    p.add_argument("--skip-rnavirhost", action="store_true", help=argparse.SUPPRESS)
    p.add_argument("--skip-phabox", action="store_true", help=argparse.SUPPRESS)
    p.add_argument("--skip-ictv", action="store_true", help=argparse.SUPPRESS)
    p.add_argument("-f", "--force", action="store_true", help="Force re-run all stages")
    args = p.parse_args()

    # 进程守护: 停止主命令时连带杀掉所有子孙进程, 避免后台残留
    try:
        import process_guard
        process_guard.install()
    except Exception:
        pass

    # mode 自动设置 skip 标志 (--skip-* 显式传入时优先)
    if args.mode != "all":
        args.skip_rnavirhost = args.skip_rnavirhost or (args.mode != "RNAVirHost")
        args.skip_phabox     = args.skip_phabox     or (args.mode != "PhaBOX2")
        args.skip_ictv       = args.skip_ictv       or (args.mode != "ICTV")

    if not os.path.isfile(args.input): sys.exit(f"ERROR: FastA file not found: {args.input}")
    os.makedirs(args.output_dir, exist_ok=True)

    rvh_csv, pb2_tsv, c9_tsv = run_tools(args)
    run_ensemble(args, rvh_csv, pb2_tsv, c9_tsv)

if __name__ == "__main__":
    main()
