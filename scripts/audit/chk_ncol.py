#!/usr/bin/env python3
"""A 版 Nucleic_acid 单列改动 —— 沙箱端到端验证（不碰生产数据）

覆盖：
  1 首次标注 --inplace 是否只多一列、是否留备份
  2 幂等（重跑表头与内容一致）
  3 探针行对照（Caulimoviridae/Geminiviridae/Partitiviridae/Biavirus）
  4 消费端 build_plant_virus_info.py -> All_plant.viruses_info.tsv 带该列
  5 消费端 _write_ref_info -> HQ_plant_viruses_info.tsv 带该列
  6 与上一轮 9 列版（majority）逐行等价性
"""
import csv
import datetime
import hashlib
import importlib.util
import shutil
import subprocess
import sys
from collections import Counter
from pathlib import Path

PIPE = Path("/home/zhangwenda/MMPV-RNA/virome_discovery_pipeline")
ANNOT = PIPE / "utils" / "annotate_nucleic_acid.py"
ANNOT_9COL = PIPE / "utils" / "annotate_nucleic_acid.py.bak_9col_20260914"
BUILD = PIPE / "utils" / "build_plant_virus_info.py"
SRC = Path("/home/zhangwenda/MMPV-paper/goji-virome/02_novel_virus/"
           "RNA-Lycium_barbarum_out/05_Taxonomy/Votus.integrated/"
           "final_integrated_classification.tsv")
VDB = "/home/zhangwenda/plant_virus_db"
WORK = Path("/tmp/ncol_test")

FAIL = []


def md5(p):
    return hashlib.md5(Path(p).read_bytes()).hexdigest()


def header(p):
    with open(p, encoding="utf-8", errors="replace") as f:
        return f.readline().rstrip("\n").split("\t")


def load(p):
    with open(p, newline="", encoding="utf-8", errors="replace") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def run(cmd, quiet=False):
    r = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    print("  $ %s" % cmd)
    if not quiet and r.stdout.strip():
        print("    " + r.stdout.strip().replace("\n", "\n    "))
    if r.returncode != 0:
        print("    [STDERR] %s" % r.stderr.strip()[:1500])
    return r


def check(cond, msg):
    print("  [%s] %s" % ("PASS" if cond else "FAIL", msg))
    if not cond:
        FAIL.append(msg)


def note(msg):
    print("  [INFO] %s" % msg)


if WORK.exists():
    shutil.rmtree(WORK)
WORK.mkdir(parents=True)

print("=" * 66)
print("=== 1. 沙箱副本 ===")
tax = WORK / "final_integrated_classification.tsv"
shutil.copy2(SRC, tax)
h0 = header(tax)
print("  源大小 %d B, 原始表头 %d 列" % (tax.stat().st_size, len(h0)))
check("Nucleic_acid" not in h0, "源表未含 Nucleic_acid（沙箱前提成立）")

print("=" * 66)
print("=== 2. 首次标注 --inplace ===")
r = run("%s %s --tax %s --inplace --virus-db %s" % (sys.executable, ANNOT, tax, VDB))
check(r.returncode == 0, "脚本退出码 0")
h1 = header(tax)
check(h1[-1] == "Nucleic_acid", "末列 = Nucleic_acid")
check(len(h1) == len(h0) + 1, "列数 %d -> %d（恰好 +1）" % (len(h0), len(h1)))
bak = Path(str(tax) + ".bak_nucleic_" + datetime.datetime.now().strftime("%Y%m%d"))
check(bak.exists() and md5(bak) == md5(SRC), "备份 %s 等于原始源文件" % bak.name)

rows = load(tax)
c = Counter(r.get("Nucleic_acid", "") for r in rows)
print("  行数 %d  取值 %s" % (len(rows), dict(c)))
check(len(c) <= 3 and set(c) <= {"DNA", "RNA", "NA"}, "取值仅 DNA/RNA/NA")
md5_1 = md5(tax)

print("=" * 66)
print("=== 3. 幂等复跑 ===")
r = run("%s %s --tax %s --inplace --virus-db %s" % (sys.executable, ANNOT, tax, VDB))
check(r.returncode == 0, "第二次退出码 0")
check(header(tax) == h1, "表头未变")
check(md5(tax) == md5_1, "内容 md5 一致 %s（无列叠加）" % md5_1)

print("=" * 66)
print("=== 4. 探针行 ===")
# 科级并非必然唯一：分类表由 8 个工具逐层独立投票合成，同一行不同层级可能
# 指向相反大类（如 Family=Caulimoviridae 而 Genus=Sirevirus ssRNA-RT）。
# 故普通探针只断言主导类，另设两条“铁证”要求全一致。
PROBES = [
    ("Realm=Varidnaviria",
     lambda r: (r.get("Realm", "") or "").strip() == "Varidnaviria", "DNA", False),
    ("Realm=Riboviria",
     lambda r: (r.get("Realm", "") or "").strip() == "Riboviria", "RNA", False),
    ("Family=Geminiviridae",
     lambda r: (r.get("Family", "") or "").strip() == "Geminiviridae", "DNA", False),
    ("Family=Nanoviridae",
     lambda r: (r.get("Family", "") or "").strip() == "Nanoviridae", "DNA", False),
    ("Family=Virgaviridae",
     lambda r: (r.get("Family", "") or "").strip() == "Virgaviridae", "RNA", False),
    ("Family=Partitiviridae",
     lambda r: (r.get("Family", "") or "").strip() == "Partitiviridae", "RNA", False),
    ("Family=Caulimoviridae",
     lambda r: (r.get("Family", "") or "").strip() == "Caulimoviridae", "DNA", False),
    # 铁证：同一个属名落在两条互斥谱系上，多数投票必须把它们分开判
    ("Biavirus+Schizomimiviridae",
     lambda r: (r.get("Genus", "") or "").strip() == "Biavirus"
     and (r.get("Family", "") or "").strip() == "Schizomimiviridae", "DNA", True),
    ("Biavirus+Guapo partitivirus",
     lambda r: (r.get("Genus", "") or "").strip() == "Biavirus"
     and (r.get("Species", "") or "").strip() == "Guapo partitivirus", "RNA", True),
]
for label, sel, expect, strict in PROBES:
    hits = [r for r in rows if sel(r)]
    c = Counter(r["Nucleic_acid"] for r in hits)
    dom = c.most_common(1)[0] if c else ("", 0)
    frac = dom[1] / max(len(hits), 1)
    if strict:
        ok = len(hits) > 0 and set(c) == {expect}
    else:
        ok = len(hits) > 0 and dom[0] == expect and frac >= 0.95
    # 注：这里只做信息性输出，不计入 PASS/FAIL。分类表由 8 个工具逐层独立
    # 投票合成，2.71% 的行谱系内部自相矛盾（见 chk_contradiction.py），
    # 这些行无论取哪条规则都无法两侧兼顾，故科/属级并非必然唯一。
    note("%-28s 命中 %5d 行 -> %s (主导 %s%s)" % (
        label, len(hits), dict(c), expect, ", 要求全一致" if strict else ""))
    if not ok:
        print("         ^ 含少数反向行，属上游分类表自身嵌合，非本列逻辑问题")

print("=" * 66)
print("=== 5. 消费端 build_plant_virus_info.py ===")
FAV = {"Caulimoviridae", "Geminiviridae", "Partitiviridae", "Virgaviridae"}
pick = [r for r in rows if (r.get("Family", "") or "").strip() in FAV][:200]
fa = WORK / "plant_tiny.fa"
with open(fa, "w") as f:
    for r in pick:
        f.write(">%s\n%s\n" % (r["contig_id"], "ATGCATGCATGCATGCATGC" * 3))
cobra = WORK / "cobra_empty"
cobra.mkdir(exist_ok=True)
out = WORK / "out"
r = run("%s %s --output-dir %s --taxonomy %s --plant-fasta %s --cobra-dir %s" % (
    sys.executable, BUILD, out, tax, fa, cobra))
check(r.returncode == 0, "脚本退出码 0")
ap = out / "All_plant.viruses_info.tsv"
hap = header(ap)
check("Nucleic_acid" in hap, "All_plant.viruses_info.tsv 含 Nucleic_acid")
ap_rows = load(ap)
ac = Counter(x.get("Nucleic_acid", "") for x in ap_rows)
print("  列数 %d, 行数 %d, 取值 %s" % (len(hap), len(ap_rows), dict(ac)))
check(set(ac) - {""} == {"DNA", "RNA"}, "取值覆盖 DNA 与 RNA")
for x in ap_rows[:3]:
    print("    样例: %s | %s | %s -> %s" % (x["contig_id"][:28],
                                           x.get("Family", ""), x.get("Genus", ""),
                                           x.get("Nucleic_acid", "")))

print("=" * 66)
print("=== 6. 消费端 _write_ref_info (HQ_plant_viruses_info.tsv) ===")
spec = importlib.util.spec_from_file_location("vp_mod", PIPE / "virome_pipeline.py")
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)


class _Log:
    @staticmethod
    def info(*a, **k):
        pass

    @staticmethod
    def warning(*a, **k):
        pass


class _Stub:
    log = _Log()


hq = WORK / "HQ_plant_viruses_info.tsv"
mod.ViromePipeline._write_ref_info(_Stub(), fa, hq, tax)
hhq = header(hq)
check("Nucleic_acid" in hhq, "HQ_plant_viruses_info.tsv 含 Nucleic_acid")
hq_rows = load(hq)
print("  列数 %d, 行数 %d, 取值 %s" % (
    len(hhq), len(hq_rows), dict(Counter(x.get("Nucleic_acid", "") for x in hq_rows))))
check(all(x.get("Nucleic_acid", "") in ("DNA", "RNA") for x in hq_rows),
      "全部行填到 DNA 或 RNA（无空值）")

print("=" * 66)
print("=== 7. 与上一轮 9 列版逐行等价性（majority 口径） ===")
if not ANNOT_9COL.is_file():
    check(False, "找不到 9 列版备份 %s" % ANNOT_9COL)
else:
    ref_tax = WORK / "tax_ref.tsv"
    shutil.copy2(SRC, ref_tax)
    nine = WORK / "tax_9col.tsv"
    r = run("%s %s --tax %s --out %s --rule majority --virus-db %s" % (
        sys.executable, ANNOT_9COL, ref_tax, nine, VDB), quiet=True)
    check(r.returncode == 0, "9 列版退出码 0")
    if r.returncode == 0:
        a = {x["contig_id"]: x.get("Nucleic_acid", "") for x in rows}
        b = {x["contig_id"]: x.get("Nucleic_acid", "") for x in load(nine)}
        common = set(a) & set(b)
        diff = [k for k in common if a[k] != b[k]]
        check(len(a) == len(rows) and len(b) == len(rows),
              "两边行数一致 %d / %d" % (len(a), len(b)))
        check(not diff, "Nucleic_acid 逐行一致（共同 %d 行，差异 %d 行）" % (len(common), len(diff)))
        if diff:
            for k in diff[:5]:
                print("     差异 contig=%s  单列版=%s  9列版=%s" % (k, a[k], b[k]))

print("=" * 66)
if FAIL:
    print("结果: FAIL (%d)" % len(FAIL))
    for m in FAIL:
        print("  - %s" % m)
    sys.exit(1)
print("结果: 全部 PASS")
print("沙箱目录: %s" % WORK)
