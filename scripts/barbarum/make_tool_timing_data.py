#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
make_tool_timing_data.py —— 生成图⑧ 左上格「7 轨道单样本平均耗时」数据表。

数据源（服务器结果树，只读）：
  <base>/02a_Identification/<sample>/<sample>_resource.tsv
      header: sample  tool  wall_sec  cpu_sec  mem_mb  status
      由 virus_identification.py::run_timed() 写入。

⚠⚠⚠ 三个必须踩对的口径（2026-10-13 第 7 次修订实测）：
  ① **`viroid_blastn` 不在 `_resource.tsv` 里**。
     源码 `virus_identification.py:1282`：
         if has_viroid and "viroid" in target_tools:
             t0 = time.time()                       # ← 算了 t0
             run_viroid_identification(...)
             log_step("类病毒")                      # ← 却**从不** log_resource()
     其余 6 个工具都走 `run_timed()` 包装 → 有行。故 viroid 行**永不存在**。
     → 用 **文件 mtime 差分**重建：`viroid_blastn.log`（= blastn 进入时创建、0 字节）
       到 `<sample>_viroids.result.id`（= 该步全部结束）的时间差 = 整步 wall time。
       ⚠ `viroid_blastn.log` 全树 **436/436 均为 0 字节**，本身不含耗时，只能靠 mtime。
       ⚠ 跨样本用 `stat -c %Y`（整秒）；单样本精度足够（该步 3–94 s）。
  ② **同一 (sample, tool) 会有多行** —— `_resource.tsv` 是**逐次追加**的：
     行序 = 执行序（genomad → blast → virbot → rdrpcatch → viralverify → metabuli），
     断点续跑 / 二次 pass 会在**末尾**再追加。故**最终口径 = 每 (sample, tool) 取最后一行**。
     实测行数：virbot 851（二 pass 普遍）、genomad 514、rdrpcatch 455、blast 399、
     viralverify 396、metabuli 396（后两者基本单行）。
     ⚠ 早期误用 `sort -u` 会**打乱追加序**；必须按文件原序读。
  ③ **样本分母 = 395**（有 `_resource.tsv` 的样本数）；viroid 覆盖 436（⊇395，多出的
     41 个是病毒支无产出、故未写 resource 的样本）。→ 画图统一口径以 **resource 的 395**
     为准时，viroid 取**其有 timing 的 436 个**的全样本均值（两者 n 不同，须在图上标注）。

用法：
  python make_tool_timing_data.py                 # 本地（需 ssh 免密到 246）
  python make_tool_timing_data.py --local-raw <已导出的 raw_resource_ordered.tsv>
输出：
  identfilter_data/tool_timing.tsv   （tool  n  mean_sec  median_sec  min_sec  max_sec）
  identfilter_data/tool_timing_meta.json
"""
import argparse, io, json, os, statistics as st, subprocess, sys, collections

HOST = "zhangwenda@202.119.189.246"
BASE = "~/MMPV-paper/goji-virome/02_novel_virus/RNA-Lycium_barbarum_out"
IDENT = BASE + "/02a_Identification"

VIRUS_TOOLS = ["genomad", "rdrpcatch", "virbot", "viralverify", "metabuli", "blast"]
PIPE_ORDER  = VIRUS_TOOLS + ["viroid_blastn"]
# 展示顺序（耗时降序，画图时再排亦可）
OUT_ORDER   = ["genomad", "rdrpcatch", "virbot", "viralverify", "metabuli", "blast", "viroid_blastn"]
LABEL = {"viroid_blastn": "viroid_blastn"}


def sh(cmd, timeout=600):
    return subprocess.check_output(["ssh", HOST, cmd], stderr=subprocess.DEVNULL,
                                   timeout=timeout).decode("utf-8", "replace")


def fetch_resource_raw():
    """按文件原序抓取全部 *_resource.tsv（保持追加序，禁 sort）"""
    cmd = ("cd %s && for d in */; do d=${d%%/}; [ \"$d\" = summary_plots ] && continue; "
           "f=$d/${d}_resource.tsv; [ -f \"$f\" ] && cat \"$f\"; done" % IDENT)
    return sh(cmd)


def fetch_viroid_timing():
    """viroid 步耗时 = mtime(result.id) − mtime(viroid_blastn.log)"""
    cmd = r"""cd %s && echo -e "sample\tviroid_wall_sec"; for d in */; do
  s=${d%%/}; [ "$s" = summary_plots ] && continue
  log="$s/viroid_output/viroid_blastn.log"; end="$s/${s}_viroids.result.id"
  [ -f "$log" ] || continue; [ -f "$end" ] || continue
  a=$(stat -c %%Y "$log" 2>/dev/null); b=$(stat -c %%Y "$end" 2>/dev/null)
  [ -z "$a" ] && continue; [ -z "$b" ] && continue
  d=$((b-a)); [ "$d" -lt 0 ] && continue
  echo -e "${s}\t${d}"
done""" % IDENT
    return sh(cmd)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default=os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                                     "identfilter_data"))
    ap.add_argument("--local-raw", help="已导出的 *_resource.tsv 原序拼接文件（离线用）")
    ap.add_argument("--local-viroid", help="已导出的 viroid timing TSV（离线用）")
    a = ap.parse_args()

    raw = io.open(a.local_raw, encoding="utf-8").read() if a.local_raw else fetch_resource_raw()
    vt  = io.open(a.local_viroid, encoding="utf-8").read() if a.local_viroid else fetch_viroid_timing()

    # ── 解析：每 (sample, tool) 取最后一行 = 最终跑 ──
    last = {}
    samples = set()
    nrow = 0
    for ln in raw.split("\n"):
        p = ln.rstrip("\r").split("\t")
        if len(p) != 6:
            continue
        s, tool, wall, cpu, mem, status = p
        if tool == "tool" or s == "sample":
            continue
        try:
            wall = float(wall)
        except ValueError:
            continue
        samples.add(s)
        last[(s, tool)] = wall          # 后写覆盖前写 → 末行胜
        nrow += 1

    # ── viroid ──
    viroid = {}
    for ln in vt.split("\n")[1:]:
        p = ln.rstrip("\n").split("\t")
        if len(p) != 2:
            continue
        try:
            viroid[p[0]] = float(p[1])
        except ValueError:
            pass

    # ── 汇总 ──
    table = []
    for t in OUT_ORDER:
        if t == "viroid_blastn":
            vals = list(viroid.values())
        else:
            vals = [w for (s, tt), w in last.items() if tt == t]
        if not vals:
            continue
        table.append(dict(tool=t, n=len(vals),
                          mean=st.mean(vals), median=st.median(vals),
                          min=min(vals), max=max(vals)))

    os.makedirs(a.outdir, exist_ok=True)
    tsv = os.path.join(a.outdir, "tool_timing.tsv")
    with io.open(tsv, "w", encoding="utf-8", newline="") as f:
        f.write("tool\tn\tmean_sec\tmedian_sec\tmin_sec\tmax_sec\n")
        for r in table:
            f.write("%s\t%d\t%.2f\t%.2f\t%.2f\t%.2f\n" %
                    (r["tool"], r["n"], r["mean"], r["median"], r["min"], r["max"]))

    meta = dict(
        source="02a_Identification/<sample>/<sample>_resource.tsv (+ viroid mtime diff)",
        n_samples_resource=len(samples),
        n_samples_viroid=len(viroid),
        resource_rows_total=nrow,
        dedup_rule="last row per (sample, tool) — file is append-only, order = run order",
        viroid_method="mtime(<sample>_viroids.result.id) - mtime(viroid_output/viroid_blastn.log)",
        viroid_absent_reason="virus_identification.py:1282 只算 t0 不调 log_resource() → 永无行",
        tools_expected=VIRUS_TOOLS + ["viroid_blastn"],
        tools_present=[r["tool"] for r in table],
    )
    with io.open(os.path.join(a.outdir, "tool_timing_meta.json"), "w", encoding="utf-8") as f:
        f.write(json.dumps(meta, ensure_ascii=False, indent=2))

    print("样本 resource=%d  viroid=%d  原始行=%d" % (len(samples), len(viroid), nrow))
    print("%-16s %6s %10s %10s" % ("tool", "n", "mean", "median"))
    for r in table:
        print("%-16s %6d %10.2f %10.2f" % (r["tool"], r["n"], r["mean"], r["median"]))
    print("→ %s" % tsv)


if __name__ == "__main__":
    main()
