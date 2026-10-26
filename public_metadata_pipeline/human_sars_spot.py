#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""human_sars_spot.py — 人类样本中 SARS-CoV-2 的 spot 分布"""

import glob
import os
import sys

import duckdb

# 服务器数据根目录 (元数据分片 + virome 汇总表都在其下)。
# 换机器/换账号时用环境变量覆盖, 不必改源码: 默认值与服务器现网路径一致, 行为不变。
_SERVER_ROOT = os.environ.get("MMPV_SERVER_ROOT", "/home/zhangwenda/plant_virus_mining")
MD = _SERVER_ROOT + "/00_data/01_raw/metadata"
VIR = _SERVER_ROOT + "/10_core/12_virome/all_runs_virus_rows.parquet"
TARGET = "Severe acute respiratory syndrome coronavirus 2"


def main() -> int:
    files = sorted(f for f in glob.glob(MD + "/*") if os.path.basename(f).startswith("20261015_"))
    print("metadata shards: %d" % len(files))

    con = duckdb.connect()
    con.execute("SET threads TO 32")
    con.execute("SET memory_limit='150GB'")

    flist = "[" + ",".join("'" + f + "'" for f in files) + "]"
    con.execute("CREATE VIEW md AS SELECT * FROM read_parquet(%s)" % flist)

    n_all = con.execute("SELECT COUNT(*) FROM md").fetchone()[0]
    print("metadata rows: %s" % format(n_all, ","))

    # ---- 人类样本 ----
    con.execute("""CREATE TABLE human AS
        SELECT acc, assay_type, librarysource, libraryselection, organism,
               releasedate, center_name, bioproject
        FROM md WHERE organism = 'Homo sapiens'""")
    n_h = con.execute("SELECT COUNT(*) FROM human").fetchone()[0]
    print("Homo sapiens runs: %s" % format(n_h, ","))

    print("\n=== human assay_type x librarysource ===")
    for r in con.execute("""SELECT assay_type, librarysource, COUNT(*) n
        FROM human GROUP BY 1,2 ORDER BY n DESC LIMIT 10""").fetchall():
        print("  %-26s %-18s %s" % (str(r[0])[:26], str(r[1])[:18], format(r[2], ",")))

    # ---- SARS-CoV-2 in human runs ----
    print("\nscanning SARS-CoV-2 in all_runs ...")
    con.execute("""CREATE TABLE sars AS
        SELECT acc, self_count, total_count, rank
        FROM read_parquet('%s') WHERE name = '%s'""" % (VIR, TARGET))
    n_s = con.execute("SELECT COUNT(*) FROM sars").fetchone()[0]
    print("all_runs SARS-CoV-2 rows: %s" % format(n_s, ","))

    con.execute("""CREATE TABLE hit AS
        SELECT h.acc, h.assay_type, h.librarysource, h.libraryselection,
               h.releasedate, h.center_name, h.bioproject,
               s.self_count, s.total_count
        FROM human h JOIN sars s ON h.acc = s.acc""")
    n_hit = con.execute("SELECT COUNT(*) FROM hit").fetchone()[0]
    print("Homo sapiens runs with SARS-CoV-2: %s" % format(n_hit, ","))

    if n_hit == 0:
        print("no hits")
        return 0

    # ---- 分布 ----
    print("\n=== spot quantiles (self_count) ===")
    for q in (0.0, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99, 1.0):
        v = con.execute(
            "SELECT CAST(quantile_cont(self_count, ?) AS BIGINT) FROM hit", [q]).fetchone()[0]
        print("  p%-5.0f %14s" % (q * 100, format(int(v), ",")))

    print("\n=== threshold passing ===")
    for t in (1, 2, 3, 5, 10, 20, 50, 100, 200, 500, 1000, 5000, 10000, 100000):
        n = con.execute("SELECT COUNT(*) FROM hit WHERE self_count >= ?", [t]).fetchone()[0]
        print("  >= %-8s %10s (%6.2f%%)" % (format(t, ","), format(n, ","),
                                            100.0 * n / n_hit))

    print("\n=== bins ===")
    for lo, hi in ((1, 1), (2, 4), (5, 9), (10, 49), (50, 99), (100, 499),
                   (500, 999), (1000, 4999), (5000, 99999), (100000, None)):
        if hi is None:
            n = con.execute("SELECT COUNT(*) FROM hit WHERE self_count >= ?", [lo]).fetchone()[0]
            lab = ">=%s" % format(lo, ",")
        elif hi == lo:
            n = con.execute("SELECT COUNT(*) FROM hit WHERE self_count = ?", [lo]).fetchone()[0]
            lab = str(lo)
        else:
            n = con.execute("SELECT COUNT(*) FROM hit WHERE self_count BETWEEN ? AND ?",
                            [lo, hi]).fetchone()[0]
            lab = "%s-%s" % (format(lo, ","), format(hi, ","))
        print("  %-16s %10s (%6.2f%%)" % (lab, format(n, ","), 100.0 * n / n_hit))

    # ---- 按 assay_type 分层 ----
    print("\n=== by assay_type ===")
    for r in con.execute("""SELECT assay_type, COUNT(*) n,
            CAST(median(self_count) AS BIGINT) med,
            CAST(quantile_cont(self_count,0.9) AS BIGINT) p90,
            MAX(self_count) mx
        FROM hit GROUP BY 1 ORDER BY n DESC LIMIT 12""").fetchall():
        print("  %-26s n=%-9s med=%-10s p90=%-10s max=%s" % (
            str(r[0])[:26], format(r[1], ","), format(int(r[2]), ","),
            format(int(r[3]), ","), format(int(r[4]), ",")))

    # ---- 按 librarysource 分层 ----
    print("\n=== by librarysource ===")
    for r in con.execute("""SELECT librarysource, COUNT(*) n,
            CAST(median(self_count) AS BIGINT) med,
            CAST(quantile_cont(self_count,0.9) AS BIGINT) p90,
            MAX(self_count) mx
        FROM hit GROUP BY 1 ORDER BY n DESC LIMIT 12""").fetchall():
        print("  %-20s n=%-9s med=%-10s p90=%-10s max=%s" % (
            str(r[0])[:20], format(r[1], ","), format(int(r[2]), ","),
            format(int(r[3]), ","), format(int(r[4]), ",")))

    # ---- 导出 ----
    out = os.path.expanduser("~/MMPV-paper/pv_spot_work/human_sars")
    os.makedirs(out, exist_ok=True)
    con.execute("COPY (SELECT * FROM hit ORDER BY self_count DESC) TO '%s/human_sars_hits.tsv' (FORMAT CSV, HEADER, DELIMITER '\t')" % out)
    con.execute("""COPY (
        SELECT assay_type, librarysource, COUNT(*) n,
               CAST(quantile_cont(self_count,0.25) AS BIGINT) p25,
               CAST(median(self_count) AS BIGINT) med,
               CAST(quantile_cont(self_count,0.75) AS BIGINT) p75,
               CAST(quantile_cont(self_count,0.9) AS BIGINT) p90,
               MAX(self_count) mx
        FROM hit GROUP BY 1,2 ORDER BY n DESC
    ) TO '%s/human_sars_by_strategy.tsv' (FORMAT CSV, HEADER, DELIMITER '\t')""" % out)
    print("\nwrote %s/human_sars_hits.tsv" % out)
    print("wrote %s/human_sars_by_strategy.tsv" % out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
