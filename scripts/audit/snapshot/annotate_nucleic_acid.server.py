#!/usr/bin/env python3
"""按分类谱系给病毒 contig 补一列 Nucleic_acid (DNA / RNA)。

判定依据（VMR_MSL41 全库实测）
------------------------------------------------------------
Family / Genus / Species 三个层级，同一层级名下跨 DNA-RNA 大界的
数量均为 0（Family 0/427, Genus 0/4149, Species 0/17554）。
即拿到科名，核酸类型就唯一确定，不需要人工裁决。
唯一存在歧义的层级区间是 Realm -> Order，且 100% 集中在反转录那一条链：
    Riboviria > Pararnavirae > Artverviricota > Revtraviricetes > Ortervirales
到 Family 级自动消解（Caulimoviridae / Hepadnaviridae = dsDNA-RT 属 DNA；
Belpaoviridae / Metaviridae / Pseudoviridae / Retroviridae = ssRNA-RT 属 RNA）。

取值规则
------------------------------------------------------------
从最细层级向上逐层回查 VMR 的 Genome 列。某层级名下多值但同属 DNA 或
同属 RNA 大类时仍可判定；真跨大界则不在该层级判定，继续上溯。

谱系上**每个能唯一判定大类的层级**各算一份证据，取多数派，即少数服从
多数；多数派里最细的层级仅用于日志回报，不额外输出列。平票时回退到最细
层级。之所以按多数而非直接取最细层级：该表由 8 个工具逐层独立投票合成，
同一行不同层级可能来自不同工具的命中，实测 barbarum 队列 20,892 行中
567 行（2.71%）谱系内部指向相反大类。最细层级优先存在已查实的失效模式，
一个孤立的假属名能压过其余 5~7 个互相印证的层级，两规则分歧 290 行，
其中 199 行为 DNA->RNA，而这 199 行里 157 行的少数派就是孤零零一个属名：
    Family=Partitiviridae(dsRNA) + Genus=Biavirus + Species=Guapo partitivirus
Biavirus 实测为 Schizomimiviridae 的 dsDNA 巨病毒属，是外来命中；该行六个
层级加种名全指向 partitivirus。取多数可纠正此类。

不靠组装长度反推类型（长度实测不可分：全量候选集 dsDNA 中位 1,870 bp
反而长于 RNA 923 bp），避免与完整度评估互相循环。

取值语义
------------------------------------------------------------
  DNA = dsDNA / ssDNA / dsDNA-RT / ssDNA(+) / ssDNA(-) / ssDNA(+/-)
  RNA = ssRNA(+) / ssRNA(-) / ssRNA(+/-) / ssRNA-RT / dsRNA
  NA  = 谱系全空或全落在多义层级，无法判定
类病毒（Pospiviroidae / Avsunviroidae）属 Subviral Agents，植物转录组里
的 Pseudoviridae / Metaviridae 是转座子，两者都不适用此二分，勿混入。

用途（RNA 测序场景）
------------------------------------------------------------
RNA 病毒在 RNA 数据里可望组装到全长，DNA 病毒通常只能拼到单基因或多基因
级别。该列用于按核酸类型分层做完整度评估、从「全长基因组」类结论里剔除
DNA 候选、以及解释 DNA contig 普遍偏短这一现象。

用法
------------------------------------------------------------
  python annotate_nucleic_acid.py --tax final_integrated_classification.tsv --inplace
  python annotate_nucleic_acid.py --tax in.tsv --out out.tsv
  python annotate_nucleic_acid.py --tax in.tsv --inplace --vmr /path/VMR_MSL41.tsv
  python annotate_nucleic_acid.py --tax in.tsv --inplace --virus-db /path/plant_virus_db

幂等：重跑会先剔除旧列再追加，不会叠加。
"""
import argparse
import csv
import glob
import os
import shutil
import sys
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path

# ── 核酸大类归属 ──
DNA_GENOMES = {"dsDNA", "ssDNA", "dsDNA-RT", "ssDNA(+)", "ssDNA(-)", "ssDNA(+/-)"}
RNA_GENOMES = {"ssRNA(+)", "ssRNA(-)", "ssRNA(+/-)", "ssRNA-RT", "dsRNA"}

RANKS = ["Realm", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
LOOKUP_ORDER = ["Genus", "Species", "Family", "Order", "Class", "Phylum", "Kingdom", "Realm"]

NEW_COL = "Nucleic_acid"
EMPTY = "NA"


def normalize_genome(g):
    """VMR 的 Genome 原始值 -> 标准细分类型（并列表述取第一个）"""
    g = (g or "").strip()
    if not g:
        return ""
    if g in DNA_GENOMES or g in RNA_GENOMES:
        return g
    head = g.split(";")[0].strip()
    if head in DNA_GENOMES or head in RNA_GENOMES:
        return head
    return g


def to_broad(genome):
    """细分类型 -> DNA / RNA / NA"""
    g = normalize_genome(genome)
    if g in DNA_GENOMES:
        return "DNA"
    if g in RNA_GENOMES:
        return "RNA"
    if not g:
        return "NA"
    if "DNA" in g:
        return "DNA"
    if "RNA" in g:
        return "RNA"
    return "NA"


def is_blank(v):
    return (not v) or v.strip() in ("NA", "nan", "None", "-", "")


# ── VMR 定位 ──
def find_vmr(explicit=None, virus_db=None):
    """按优先级寻找 VMR_MSL*.tsv。返回 (Path|None, 说明文本)。"""
    cands = []
    if explicit:
        cands.append(Path(explicit))
    env = os.environ.get("MMPV_VMR")
    if env:
        cands.append(Path(env))
    roots = []
    if virus_db:
        roots.append(Path(virus_db))
    roots += [
        Path("/home/zhangwenda/plant_virus_db"),
        Path.home() / "plant_virus_db",
        Path(__file__).resolve().parent.parent / "plant_virus_db",
    ]
    for root in roots:
        if not root.is_dir():
            continue
        for pat in ("1.virus-host_db/B-ictv/VMR_MSL*.tsv",
                    "*/1.virus-host_db/B-ictv/VMR_MSL*.tsv",
                    "**/VMR_MSL*.tsv"):
            cands += [Path(p) for p in glob.glob(str(root / pat), recursive=True)]

    def vkey(p):
        """VMR_MSL41 -> 41，用于挑最新版本"""
        stem = p.stem.replace("VMR_MSL", "")
        return int(stem) if stem.isdigit() else -1

    seen, uniq = set(), []
    for c in cands:
        rp = str(c.resolve()) if c.exists() else str(c)
        if rp in seen:
            continue
        seen.add(rp)
        uniq.append(c)
    good = [c for c in uniq if c.is_file() and c.stat().st_size > 1000]
    if good:
        return max(good, key=vkey), "ok"
    return None, "已搜索: " + "; ".join(str(c) for c in uniq[:8])


def load_vmr(path):
    """读 VMR 的 Genome 列 -> {层级: {名字: Counter(细分类型)}}"""
    lv = {r: defaultdict(Counter) for r in RANKS}
    with open(path, newline="", encoding="utf-8", errors="replace") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            g = normalize_genome(row.get("Genome"))
            if not g:
                continue
            for r in RANKS:
                v = (row.get(r) or "").strip()
                if v and v != "NA":
                    lv[r][v][g] += 1
    return {"lv": lv}


def resolve_rank(value, rank, idx):
    """在某一层级查大类。返回 DNA / RNA，或 None（空白 / 该层未收录 / 真跨大界）"""
    if is_blank(value):
        return None
    gc = idx["lv"][rank].get(value.strip())
    if not gc:
        return None
    votes = defaultdict(int)
    for g, n in gc.items():
        votes[to_broad(g)] += n
    decided = {b: n for b, n in votes.items() if b in ("DNA", "RNA")}
    if len(decided) != 1:
        return None                      # 该层全落在 NA，或真跨大界
    return list(decided)[0]


def decide(row, idx):
    """取谱系上各层证据的多数派。返回 (broad, rank)。

    多义层级（名下同时含 DNA 与 RNA 谱系，如 Pararnavirae）在 resolve_rank()
    阶段就返回 None，不会进入投票，故反转录嵌套不会误判。
    """
    per = {}
    for rank in LOOKUP_ORDER:
        b = resolve_rank(row.get(rank, ""), rank, idx)
        if b:
            per[rank] = b
    if not per:
        return EMPTY, "none"
    top = Counter(per.values()).most_common()
    if len(top) > 1 and top[0][1] == top[1][1]:
        for r in LOOKUP_ORDER:           # 平票 -> 回退最细层级
            if r in per:
                return per[r], r
    broad = top[0][0]                    # 多数派里最细的层级（仅用于日志）
    rank = min((r for r, b in per.items() if b == broad), key=LOOKUP_ORDER.index)
    return broad, rank


def main():
    ap = argparse.ArgumentParser(
        description="按分类谱系给病毒 contig 补一列 Nucleic_acid (DNA/RNA)",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tax", required=True,
                    help="含分类谱系的 TSV（final_integrated_classification.tsv）")
    ap.add_argument("--out", default=None,
                    help="输出 TSV（默认 <tax>.nucleic.tsv；--inplace 时忽略）")
    ap.add_argument("--inplace", action="store_true",
                    help="就地写回 --tax（首次自动备份 .bak_nucleic_<日期>）")
    ap.add_argument("--vmr", default=None, help="VMR_MSL*.tsv 路径（默认自动查找）")
    ap.add_argument("--virus-db", default=None,
                    help="病毒库根目录，用于自动定位 VMR（<virus_db>/1.virus-host_db/B-ictv/）")
    args = ap.parse_args()

    tax_path = Path(args.tax)
    if not tax_path.is_file():
        print("[ERROR] 找不到输入: %s" % tax_path, file=sys.stderr)
        return 1

    vmr, why = find_vmr(args.vmr, args.virus_db)
    if vmr:
        print("[INFO] VMR: %s" % vmr)
        idx = load_vmr(vmr)
    else:
        print("[WARN] 未找到 VMR_MSL*.tsv，核酸列将全部写 NA。\n"
              "       用 --vmr 显式指定，或设环境变量 MMPV_VMR 后重跑。\n"
              "       %s" % why, file=sys.stderr)
        idx = {"lv": {r: {} for r in RANKS}}

    with open(tax_path, newline="", encoding="utf-8", errors="replace") as f:
        rd = csv.DictReader(f, delimiter="\t")
        cols = rd.fieldnames or []
        rows = list(rd)
    if not cols:
        print("[ERROR] 输入为空: %s" % tax_path, file=sys.stderr)
        return 1

    # 幂等：剔除旧列后重新追加
    out_cols = [c for c in cols if c != NEW_COL] + [NEW_COL]

    n_acid, n_rank = Counter(), Counter()
    out_rows = []
    for row in rows:
        broad, rank = decide(row, idx)
        r2 = dict(row)
        r2[NEW_COL] = broad
        out_rows.append(r2)
        n_acid[broad] += 1
        n_rank[rank] += 1

    if args.inplace:
        bak = tax_path.with_name(tax_path.name + ".bak_nucleic_" +
                                 datetime.now().strftime("%Y%m%d"))
        if not bak.exists():
            shutil.copy2(tax_path, bak)
            print("[INFO] 备份: %s" % bak)
        out_path = tax_path
    else:
        out_path = Path(args.out) if args.out else tax_path.with_suffix(".nucleic.tsv")

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=out_cols, delimiter="\t",
                           extrasaction="ignore", lineterminator="\n")
        w.writeheader()
        for r in out_rows:
            w.writerow({k: r.get(k, "") for k in out_cols})

    n = len(out_rows)
    ok = n_acid["DNA"] + n_acid["RNA"]
    print("Nucleic_acid 补全完成 -> %s" % out_path)
    print("  总行数    : %d" % n)
    print("  已判定    : %d (%.2f%%)" % (ok, 100.0 * ok / max(n, 1)))
    print("  DNA / RNA : %d / %d" % (n_acid["DNA"], n_acid["RNA"]))
    print("  未定(NA)  : %d (%.2f%%)" % (n_acid["NA"], 100.0 * n_acid["NA"] / max(n, 1)))
    print("  判定层级  : %s" % dict(n_rank.most_common()))
    return 0


if __name__ == "__main__":
    sys.exit(main())
