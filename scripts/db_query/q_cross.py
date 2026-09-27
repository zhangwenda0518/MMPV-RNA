import csv
import collections

PLANT = "/home/zhangwenda/plant_virus_db/1.virus-host_db/C-host_classify/classified_clean/Plant.tsv"
FINAL = "/home/zhangwenda/plant_virus_db/1.virus-host_db/C-host_classify/host_extract/Final_Virus_Host_Lineage.tsv"

# 1) 权威库 Plant.tsv: 收集所有合法植物病毒的科 / 属
auth_fam = set()
auth_gen = set()
with open(PLANT, errors="replace") as f:
    r = csv.DictReader(f, delimiter="\t")
    for row in r:
        lin = (row.get("Virus_lineage", "") or "")
        for p in lin.split(";"):
            p = p.strip()
            if not p:
                continue
            if p.endswith("viridae"):
                auth_fam.add(p)
            # 属: 特征是以 virus 结尾且不是科
            if p.endswith("virus") and not p.endswith("viridae"):
                auth_gen.add(p)

print(f"权威库 Plant.tsv: 科 {len(auth_fam)} 个, 属 {len(auth_gen)} 个")

# 2) 权威库全量: 收集所有非植物科/属 (原生生物/昆虫/真菌/藻)
final_hosts = collections.defaultdict(collections.Counter)
kin = {"protist": 0, "fungi": 0, "metazoa": 0, "algae": 0, "bacteria": 0}
nonplant_fam = set()
nonplant_gen = set()
with open(FINAL, errors="replace") as f:
    r = csv.DictReader(f, delimiter="\t")
    for row in r:
        lin = (row.get("Host_lineage", "") or "").lower()
        vlin = (row.get("Virus_lineage", "") or "")
        fam = ""
        gen = ""
        for p in vlin.split(";"):
            p = p.strip()
            if p.endswith("viridae"):
                fam = p
            elif p.endswith("virus"):
                gen = p
        is_nonplant = False
        if "amoebozoa" in lin or "discosea" in lin or "haptophyta" in lin or "alveolata" in lin or "stramenopiles" in lin:
            is_nonplant = True; kin["protist"] += 1
        elif "fungi" in lin:
            is_nonplant = True; kin["fungi"] += 1
        elif "metazoa" in lin or "arthropoda" in lin or "insecta" in lin:
            is_nonplant = True; kin["metazoa"] += 1
        elif "viridiplantae" in lin or "streptophyta" in lin:
            pass
        if is_nonplant:
            if fam: nonplant_fam.add(fam)
            if gen: nonplant_gen.add(gen)

print(f"权威库全量: 非植物宿主的科 {len(nonplant_fam)} 个, 属 {len(nonplant_gen)} 个")
print(f"  宿主类型计数: {dict(kin)}")

# 3) 关键: 在 C9 里判为 Plant 的科, 有哪些在权威库 Plant.tsv 里不存在?
c9_plant_fams = ['Partitiviridae','Caulimoviridae','Potyviridae','Tombusviridae','Tymoviridae',
    'Secoviridae','Rhabdoviridae','Endornaviridae','Closteroviridae','Geminiviridae',
    'Bromoviridae','Virgaviridae','Betaflexiviridae','Solemoviridae','Alphaflexiviridae',
    'Mimiviridae','Tospoviridae','Phenuiviridae','Amalgaviridae','Kitaviridae',
    'Nanoviridae','Spinareoviridae','Aspiviridae','Phycodnaviridae','Atkinsviridae',
    'Benyviridae','Tobaliviridae','Ambiguiviridae','Geplanaviridae','Metaviridae']

print("\n=== C9 判 Plant 的科 vs 权威库 ===")
print(f"{'科':<22} {'权威库Plant':<12} {'权威库非植物':<12}")
for fam in c9_plant_fams:
    in_auth = "✓" if fam in auth_fam else "✗ 无"
    in_non = "★有" if fam in nonplant_fam else "-"
    print(f"{fam:<22} {in_auth:<12} {in_non:<12}")
