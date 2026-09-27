#!/usr/bin/env python3
"""一次性驱动: 合成 run 目录 + 旧 stop_profile 基线 -> 交给 rerun_diff.py 对比新代码.

设计: 8 条确定性 contigs (coding_intact / distributed_decay / assembly_breakpoint /
边界翻转 / compositional_noise / 无命中 / TE) + 30 条随机 contigs 估自然翻转率。
基线 verdict 用旧版 stop_profile (3nt 高估) 产出, rerun_diff 用当前代码重算。
填充只用 GCT: 与 TAA 的全部拼接 3-mer (CTT/TTA/TAA/AAG/AGC) 里只有 TAA 是终止,
所以每个 TAA 恰好制造一个终止, 终止位置可精确断言。
"""
import random
import shutil
import subprocess
import sys
from pathlib import Path

MOD = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(MOD))
import s1_decay_scan as s1  # noqa: E402  (借用 scan_frame_stops 做设计验证)

FIX = Path(sys.argv[1] if len(sys.argv) > 1 else MOD / "verify/_selftest_fixture")
FIX.mkdir(parents=True, exist_ok=True)

SPAN = 600
FILLER = "GCT"


def build_region(spec):
    """spec: {position: 任一非零值}; 在绝对位置 p 覆写 TAA, 其余 GCT."""
    seq = [FILLER] * (SPAN // 3)
    seq = "".join(seq)
    for p in sorted(spec):
        assert 0 <= p <= SPAN - 3, p
        seq = seq[:p] + "TAA" + seq[p + 3:]
    # 设计验证: 用 s1 同一个扫描器数终止, 每帧位置必须与设计一致
    for off in range(3):
        got = scan_abs(seq, off)
        want = sorted(p for p in spec if p % 3 == off)
        assert got == want, f"frame+{off+1} 设计不符: got {got} want {want}"
    return seq


def scan_abs(seq, offset):
    """s1.scan_frame_stops 的直通版 (它返回的已是相对区域起点的绝对坐标)."""
    return s1.scan_frame_stops(seq, 0, len(seq), offset)


def contig(cid, region, flanks=()):
    """contig 就是设计区本身 (不带侧翼): HSP qstart=1..600 与设计坐标一一同格."""
    return f">{cid}\n{region}\n"


def hsp_row(cid, comps, pident=35.0, aln_aa=250, bits=200.0):
    return (f"{cid}\t{comps}|TEST01.1|pol_test\t{pident}\t{aln_aa}\t1\t{SPAN}\t"
            f"1\t{aln_aa}\t1e-20\t{bits}\t1\tpolyprotein [Test virus]\n")


rng = random.Random(20260927)
fa, panel, baits, ev = [], [], [], []
ev_hdr = "contig_id\ttax_family\tcheckv_completeness\taa_pident\taa_species\tnt_pident\tnt_species\n"
ev.append(ev_hdr)


def add(cid, region, comps=None, te_bait=False, panel_pident=35.0):
    fa.append(contig(cid, region))
    if comps:
        panel.append(hsp_row(cid, comps, pident=panel_pident))
    if te_bait:
        baits.append(hsp_row(cid, "FALSE_pol_TE1", pident=90.0, aln_aa=200, bits=100.0))
    ev.append(f"{cid}\t\t\t\t\t\t\n")


DISTRIB = [30, 120, 210, 300, 390, 480]
# 边界翻转: 内部间隔 180->330 = 150nt, 旧 msf=150/600=0.25 (>=0.25 -> breakpoint),
# 新 msf=147/600=0.245 (<0.25 -> 分布式退化)。其余间隔和首末段都 < 150。
FLIP = [30, 180, 186, 189, 198, 228, 258, 288, 318, 348, 378, 408, 438, 468,
        498, 528, 558]
# compositional_noise: +1 两处, +2 两处, +3 三处 -> enrich=(2+1)/(2.5+1)=0.857 <1 且 >0.5
NOISE = [300, 597, 101, 403, 202, 352, 502]

add("c_intact_mp", build_region({}), "MP+CP")                       # -> virus_fragment_review
add("c_decay_rt", build_region({p: 1 for p in DISTRIB}), "RT")      # -> ancient_EVE
add("c_decay_mp", build_region({p: 1 for p in DISTRIB}), "MP+CP")   # -> EVE_suspect
add("c_break_mp", build_region({9: 1, 30: 1}), "MP+CP")             # -> assembly_breakpoint_review
add("c_flip_rt", build_region({p: 1 for p in FLIP}), "RT")          # 旧 review -> 新 ancient_EVE
add("c_noise_mp", build_region({p: 1 for p in NOISE}), "MP+CP")     # -> compositional_noise_review
add("c_none", build_region({}), None)                               # -> review (无 HSP)
add("c_te", build_region({}), "RT", te_bait=True, panel_pident=22.0)  # -> EVE_LTR_TE
# 随机 contigs: 随机槽位终止 (都落在 +1 帧, 组成噪声主要靠 msf 分布)
for i in range(30):
    cid = f"r{i:02d}"
    spec = {rng.randrange(0, SPAN // 3) * 3: 1 for _ in range(rng.randint(0, 12))}
    fa.append(contig(cid, build_region(spec)))
    panel.append(hsp_row(cid, "RT"))
    ev.append(f"{cid}\t\t\t\t\t\t\n")

# 设计意图自检: 当前 (新) 代码下, 8 条设计 contig 的 decay_class 必须命中预期
EXPECT_DECAY = {
    "c_intact_mp": "coding_intact", "c_decay_rt": "distributed_decay",
    "c_decay_mp": "distributed_decay", "c_break_mp": "assembly_breakpoint",
    "c_flip_rt": "distributed_decay",   # 新代码: msf=0.245 < 0.25
    "c_noise_mp": "compositional_noise", "c_te": "coding_intact",
    "c_none": "no_hsp",                 # 无任何 HSP -> decay=no_hsp
}

(FIX / "q.fa").write_text("".join(fa), encoding="utf-8", newline="\n")
(FIX / "panel_hits.tsv").write_text("".join(panel), encoding="utf-8", newline="\n")
(FIX / "baits_hits.tsv").write_text("".join(baits), encoding="utf-8", newline="\n")
(FIX / "locus_blastn.tsv").write_text("", encoding="utf-8", newline="\n")
(FIX / "rescue_evidence_scored.tsv").write_text("".join(ev), encoding="utf-8", newline="\n")

# ── 旧版基线: 复制 s1 并把 runs 计算换回 3nt 高估版 ──
cur = (MOD / "s1_decay_scan.py").read_text(encoding="utf-8")
NEW_BLOCK = ("    runs = ([pos[0]]\n"
             "            + [b - a - 3 for a, b in zip(pos, pos[1:])]\n"
             "            + [span_nt - pos[-1] - 3])")
OLD_BLOCK = ("    runs = [pos[0]] + [b - a for a, b in zip(pos, pos[1:])] "
             "+ [span_nt - pos[-1]]")
assert NEW_BLOCK in cur, "新版 runs 计算块没找到, s1 可能又被改过"
(FIX / "_old_s1.py").write_text(cur.replace(NEW_BLOCK, OLD_BLOCK), encoding="utf-8")


def run(*args):
    r = subprocess.run([sys.executable, *args], capture_output=True, text=True)
    if r.returncode != 0:
        sys.exit(f"FAILED rc={r.returncode}: {args}\n{r.stderr}")
    return r.stdout


# 翻转用例的算术自检: 新旧 msf 必须正好卡在 0.25 两侧
region = build_region({p: 1 for p in FLIP})
pos_new = scan_abs(region, 0)
prof_new = s1.stop_profile(pos_new, SPAN)
runs_old = [pos_new[0]] + [b - a for a, b in zip(pos_new, pos_new[1:])] + [SPAN - pos_new[-1]]
msf_old = round(max(runs_old) / SPAN, 4)
print(f"[fixture] c_flip_rt msf: old={msf_old} (>=0.25 -> assembly_breakpoint) "
      f"new={prof_new['max_stopfree_frac']} (<0.25 -> distributed_decay)")
assert msf_old >= 0.25 > prof_new["max_stopfree_frac"]

# 旧代码全链: s1_old -> s2 -> s2b(空位点表) -> s3, 产物作为 rerun_diff 的"旧版基准"
all_hits = FIX / "all_hits.tsv"
all_hits.write_text((FIX / "panel_hits.tsv").read_text(encoding="utf-8")
                    + (FIX / "baits_hits.tsv").read_text(encoding="utf-8"),
                    encoding="utf-8", newline="\n")
run(str(FIX / "_old_s1.py"), str(FIX / "q.fa"), str(all_hits), str(FIX / "s1_decay.tsv"))
run(str(MOD / "s2_domain_scan.py"), str(FIX / "q.fa"), str(FIX / "panel_hits.tsv"),
    str(FIX / "baits_hits.tsv"), str(FIX / "s2_domains.tsv"))
hdr = run(str(MOD / "s2b_locus_scan.py"), "--emit-header")
(FIX / "locus_architecture.tsv").write_text(hdr, encoding="utf-8", newline="\n")
run(str(MOD / "s3_verdict.py"), str(FIX / "s1_decay.tsv"), str(FIX / "s2_domains.tsv"),
    str(FIX / "locus_architecture.tsv"), str(FIX / "rescue_evidence_scored.tsv"),
    str(FIX / "eve_distinguish_verdict.tsv"))
print(f"[fixture] 合成 run 目录就绪: {FIX}")
print("[fixture] 旧版基线 verdict 已就位, 交给 rerun_diff.py 用当前代码重算对比:\n")
r = subprocess.run([sys.executable, str(MOD / "verify" / "rerun_diff.py"),
                    "--run-dir", str(FIX),
                    "--evidence", str(FIX / "rescue_evidence_scored.tsv")],
                   capture_output=True, text=True)
sys.stdout.write(r.stdout)
sys.stderr.write(r.stderr)

# ── 产物级断言: 设计 contig 的 decay/verdict 必须命中预期 (新代码列) ──
EXPECT_VERDICT = {
    "c_intact_mp": "virus_fragment_review", "c_decay_rt": "ancient_EVE",
    "c_decay_mp": "EVE_suspect", "c_break_mp": "assembly_breakpoint_review",
    "c_flip_rt": "ancient_EVE", "c_noise_mp": "compositional_noise_review",
    "c_none": "review", "c_te": "EVE_LTR_TE",
}
EXPECT_DECAY_OLD = dict(EXPECT_DECAY, c_flip_rt="assembly_breakpoint")


def col_of(path, cid, name):
    for line in open(path, encoding="utf-8"):
        p = line.rstrip("\n").split("\t")
        if p[0] == cid:
            hdr = open(path, encoding="utf-8").readline().rstrip("\n").split("\t")
            return p[hdr.index(name)]
    return "<missing>"


fails = []
for cid in EXPECT_DECAY:
    for tag, path, exp in (
            ("old-decay", FIX / "s1_decay.tsv", EXPECT_DECAY_OLD[cid]),
            ("new-decay", FIX / "rerun_currentcode" / "s1_decay.tsv", EXPECT_DECAY[cid]),
            ("new-verdict", FIX / "rerun_currentcode" / "eve_distinguish_verdict.tsv",
             EXPECT_VERDICT[cid])):
        got = col_of(path, cid, "decay_class" if "decay" in tag else "verdict")
        if got != exp:
            fails.append(f"{cid} {tag}: got {got} want {exp}")
# 旧 verdict 基准里 c_flip_rt 必须是 assembly_breakpoint_review (正是被修复翻转的那个)
if col_of(FIX / "eve_distinguish_verdict.tsv", "c_flip_rt", "verdict") != \
        "assembly_breakpoint_review":
    fails.append("c_flip_rt 旧版 verdict 不是 assembly_breakpoint_review, 翻转用例失效")

if fails:
    print("\n[fixture] 设计断言失败:")
    for f in fails:
        print("  " + f)
    print(f"[fixture] 现场保留在 {FIX} 供排查 (stage 脚本 Windows 默认文本写盘会带 CRLF, "
          "排查完请手动删掉, 否则 test_eve_invariants 的 CRLF 检查会咬它)")
    sys.exit(1)
shutil.rmtree(FIX, ignore_errors=True)
print("\n[fixture] 设计断言全部通过: 8 条设计 contig 的 decay/verdict 逐条命中预期; "
      "c_flip_rt 按设计翻转 (assembly_breakpoint_review -> ancient_EVE), "
      "其余 37 条 verdict 不动")
print(f"[fixture] 产物目录 {FIX} 已清理 (重跑即再生)")
sys.exit(0)
