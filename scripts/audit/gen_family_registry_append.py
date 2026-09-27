# -*- coding: utf-8 -*-
"""从审计产物生成 NON_PLANT_FAMILIES_FALLBACK 追加块文本 (纯追加补丁的载荷)。

输入: scripts/audit/non_plant_family_registry.tsv  (逐科计数 + 代表属)
输出: scripts/patch/_family_registry_append.txt
用法: python scripts/audit/gen_family_registry_append.py
"""
import csv
import io
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
TSV = os.path.join(HERE, "non_plant_family_registry.tsv")
OUT = os.path.join(ROOT, "scripts", "patch", "_family_registry_append.txt")

PLANT_DB_MD5 = "6e3200b1b2a769610f0780979a8dec26"
DATE = "2026-09-15"


def main():
    rows = list(csv.DictReader(open(TSV, encoding="utf-8"), delimiter="\t"))
    def n(r, k):
        try:
            return int(r.get(k) or 0)
        except ValueError:
            return 0

    plant = [r for r in rows if n(r, "gojiPL") or n(r, "onekpPL")]
    plant.sort(key=lambda r: -(n(r, "gojiPL") + n(r, "onekpPL")))
    rest = sorted(r["Family"] for r in rows if not (n(r, "gojiPL") or n(r, "onekpPL")))
    assert len(plant) + len(rest) == len(rows), (len(plant), len(rest), len(rows))

    b = io.StringIO()
    w = b.write
    w("# ── %s 追加: 非植物科登记扩容 (%d 科) ──\n" % (DATE, len(rows)))
    w("# 判据: ① 该科在 PlantVirusDB 权威植物库 classified_clean/Plant.tsv 中零记录\n")
    w("#         (Plant.tsv md5 %s; 201,488 条植物病毒记录 / 解析 46 科 / 821 属)\n" % PLANT_DB_MD5)
    w("#       ② 该科在枸杞(goji)或 OneKP 树 05 分类表 / 下游植物表中出现\n")
    w("# 逐科行数与代表属证据: scripts/audit/non_plant_family_registry.tsv\n")
    w("# 宿主域逐个核实结论: scripts/audit/TAXONOMY_05_CALIB_STATE.md\n")
    w("# 作用面: 仅在 C9 只给到科级且属缺失时生效; 属/种级判定不经过本名单\n")
    w("#         (属级短路见 is_blacklisted 的 is_trusted_level 分支)\n")
    w("NON_PLANT_FAMILIES_FALLBACK += [\n")
    w("    # 出现在下游植物病毒结果表中的 %d 科 (承重条目, 逐科核定)\n" % len(plant))
    for r in plant:
        f = r["Family"]
        top = (r.get("代表属(Top3)") or "").strip()
        top = top[:60] if top else "无属记录"
        w("    '%s',%s# 05:g%s/o%s 植物表:g%s/o%s | %s\n" % (
            f, " " * max(1, 22 - len(f)),
            r["goji05"], r["onekp05"], r["gojiPL"], r["onekpPL"], top))
    w("    # 仅出现在 05 分类表中的 %d 科 (按字母序, 不参与当前下游表清理)\n" % len(rest))
    for f in rest:
        w("    '%s',\n" % f)
    w("]\n")
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    open(OUT, "w", encoding="utf-8", newline="\n").write(b.getvalue())
    print("写出 %s: 共 %d 科 (植物表出现 %d + 仅 05 表 %d)" % (OUT, len(rows), len(plant), len(rest)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
