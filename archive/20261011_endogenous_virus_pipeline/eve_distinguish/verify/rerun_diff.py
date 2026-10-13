#!/usr/bin/env python3
"""rerun_diff.py — 用当前代码从原始命中重算 s1→s4, 与 run 目录里的旧 verdict 逐条 diff.

背景: s1/s2/s2b/s3 都是纯计算, 只有 blastx 与 host blastn 贵。改了判定逻辑之后,
正确做法不是相信文档里的旧数字, 而是从原始命中重算一遍, 对出"哪些候选的 verdict
变了、往哪个方向变"。这个脚本把那套手工命令 (README"复用派生表的坑"一节) 固化成
一个工具 —— 判定逻辑任何改动之后, 对既有 run 目录发一条命令即可量化影响面。

用法:
  python3 verify/rerun_diff.py --run-dir <既有run目录> --evidence <rescue_evidence_scored.tsv> \
      [--host-map <tsv>] [--out <目录>]

输入约定 (与 run_all.sh 的产物同名, 原始命中不重跑 diamond):
  <run-dir>/q.fa (或 evescan_query.fasta)     候选 fasta
  <run-dir>/panel_hits.tsv / baits_hits.tsv   原始 blastx 命中
  <run-dir>/locus_blastn.tsv                  原始寄主 blastn (--no-host 的 run 可以没有:
                                              此时按 --emit-header 走空位点表, 与 run_all 一致)
  <run-dir>/eve_distinguish_verdict.tsv       旧 verdict (diff 基准; 没有则只报新分布)

输出 (默认 <run-dir>/rerun_currentcode/, --out 可覆盖):
  s1_decay.tsv / s2_domains.tsv / locus_architecture.tsv /
  eve_distinguish_verdict.tsv / dna_vs_eve_filter.tsv     当前代码的重算产物
  verdict_diff.tsv                                        仅列 verdict 变了的行
stdout: verdict 迁移矩阵 + decay_class 迁移矩阵 (旧行 -> 新列)。
迁移矩阵非对角即判定逻辑改动在真实数据上的影响面; 全对角 = 产出与旧版一致。
"""
import argparse
import collections
import subprocess
import sys
import tempfile
from pathlib import Path

MOD = Path(__file__).resolve().parent.parent     # eve_distinguish/
S1, S2, S2B, S3, S4 = (MOD / n for n in
                       ("s1_decay_scan.py", "s2_domain_scan.py", "s2b_locus_scan.py",
                        "s3_verdict.py", "s4_filter.py"))


def find_query(run_dir):
    for name in ("q.fa", "evescan_query.fasta", "query.fasta"):
        p = run_dir / name
        if p.is_file() and p.stat().st_size:
            return p
    sys.exit(f"[rerun] 在 {run_dir} 找不到候选 fasta (q.fa / evescan_query.fasta)")


def read_col(path, col):
    """按表头列名读一列: {首列: 该列值}; 文件不存在/无表头返回 {}."""
    if not path or not Path(path).is_file():
        return {}
    with open(path, encoding="utf-8") as f:
        hdr = f.readline().rstrip("\n").split("\t")
        if col not in hdr:
            return {}
        ci = hdr.index(col)
        out = {}
        for line in f:
            p = line.rstrip("\n").split("\t")
            if p:
                out[p[0]] = p[ci] if ci < len(p) else ""
        return out


def transition(old, new, label):
    """迁移矩阵 (旧行 -> 新列), 只打非对角线明细."""
    keys = sorted(set(old.values()) | set(new.values()))
    mat = collections.Counter((old.get(k, "<absent>"), new.get(k, "<absent>"))
                              for k in set(old) | set(new))
    offdiag = {t: n for t, n in mat.items() if t[0] != t[1]}
    print(f"\n[{label}] 迁移矩阵 (共 {sum(mat.values())} 条, 判定值 {len(keys)} 种)")
    for (o, n) in sorted(mat):
        mark = "  " if o == n else "->"
        print(f"  {mark} {o:32s} => {n:32s} {mat[(o, n)]}")
    return offdiag


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run-dir", required=True)
    ap.add_argument("--evidence", required=True, help="rescue_evidence_scored.tsv")
    ap.add_argument("--host-map", default=str(MOD / "host_sample_map.tsv"))
    ap.add_argument("--out", default="")
    a = ap.parse_args()

    run_dir = Path(a.run_dir)
    out = Path(a.out) if a.out else run_dir / "rerun_currentcode"
    out.mkdir(parents=True, exist_ok=True)
    q = find_query(run_dir)
    panel = run_dir / "panel_hits.tsv"
    baits = run_dir / "baits_hits.tsv"
    for p in (panel, baits, Path(a.evidence)):
        if not p.is_file():
            sys.exit(f"[rerun] 缺输入: {p}")
    locus_blastn = run_dir / "locus_blastn.tsv"
    old_verdict = run_dir / "eve_distinguish_verdict.tsv"

    def sh(*args):
        r = subprocess.run([sys.executable, *args], capture_output=True, text=True)
        if r.returncode != 0:
            sys.exit(f"[rerun] 阶段失败 ({args[0]} rc={r.returncode}):\n{r.stderr}")
        tail = (r.stdout or "").strip().splitlines()
        if tail:
            print(f"  {tail[-1]}")

    all_hits = out / "all_hits.tsv"
    with open(all_hits, "w", encoding="utf-8", newline="\n") as w:
        for src in (panel, baits):
            if src.is_file():
                w.write(src.read_text(encoding="utf-8", errors="replace"))

    print(f"[rerun] run_dir={run_dir}\n[rerun] out={out}\n[rerun] query={q.name}")
    s1_tsv = out / "s1_decay.tsv"
    s2_tsv = out / "s2_domains.tsv"
    arch_tsv = out / "locus_architecture.tsv"
    new_verdict = out / "eve_distinguish_verdict.tsv"

    sh(str(S1), str(q), str(all_hits), str(s1_tsv))
    sh(str(S2), str(q), str(panel), str(baits), str(s2_tsv))
    if locus_blastn.is_file():
        sh(str(S2B), str(locus_blastn), str(s2_tsv), str(s1_tsv), str(arch_tsv),
           a.host_map if Path(a.host_map).is_file() else "")
    else:
        print("  (无 locus_blastn.tsv: 按 --no-host 口径出空位点表)")
        hdr = subprocess.run([sys.executable, str(S2B), "--emit-header"],
                             capture_output=True, text=True, check=True).stdout
        arch_tsv.write_text(hdr, encoding="utf-8")
    sh(str(S3), str(s1_tsv), str(s2_tsv), str(arch_tsv), str(a.evidence), str(new_verdict))
    sh(str(S4), str(new_verdict), "-o", str(out / "dna_vs_eve_filter.tsv"))

    old_v = read_col(old_verdict, "verdict")
    new_v = read_col(new_verdict, "verdict")
    if not old_v:
        print(f"\n[rerun] {old_verdict.name} 不存在, 只报新分布 (无 diff 基准)")
        return 0
    # decay_class 变化是 verdict 变化的上游原因, 两张矩阵一起给
    old_d = read_col(old_verdict, "decay_class")
    new_d = read_col(new_verdict, "decay_class")

    diff_rows = [(cid, old_v.get(cid, "<absent>"), new_v.get(cid, "<absent>"))
                 for cid in sorted(set(old_v) | set(new_v))
                 if old_v.get(cid, "<absent>") != new_v.get(cid, "<absent>")]
    with open(out / "verdict_diff.tsv", "w", encoding="utf-8", newline="\n") as w:
        w.write("contig_id\told_verdict\tnew_verdict\n")
        for row in diff_rows:
            w.write("\t".join(row) + "\n")

    transition(old_d, new_d, "decay_class")
    offdiag = transition(old_v, new_v, "verdict")
    print(f"\n[rerun] verdict 变更 {len(diff_rows)} 条 -> {out / 'verdict_diff.tsv'}")
    if offdiag:
        for (o, n), k in sorted(offdiag.items(), key=lambda x: -x[1])[:10]:
            print(f"  {k:5d}  {o} -> {n}")
        return 1
    print("[rerun] 全部 verdict 与旧版一致 (迁移矩阵无非对角项)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
