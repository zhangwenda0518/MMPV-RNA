#!/usr/bin/env python3
"""terminal_repeats.py 单元测试 — DTR/ITR 检测、修剪与起点归一。

覆盖 8 类用例:
  1. 精确 DTR (contig = P + M + P)          -> fetch_dtr 命中
  2. 容错 DTR (3' 端拷贝单碱基错配)          -> exact 落空 / fuzzy 命中
  3. ITR (I + M + RC(I))                    -> fetch_itr 命中
  4. 无重复随机序列                          -> 阴性
  5. 过短序列 / 超长序列跳过                 -> 边界
  6. 用户声明环状 (--circ-file)              -> 不搜 DTR、不修剪
  7. 起点归一 repeat / orf / none            -> 旋转语义
  8. CLI 端到端                              -> 三件产物 + 汇总字段

用法:
    python test_terminal_repeats.py             # 跑全部 (临时目录自动清理)
    python test_terminal_repeats.py --keep      # 保留临时产物便于人工核对
    python -m pytest test_terminal_repeats.py -v
"""

import os
import random
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import terminal_repeats as tr  # noqa: E402

PIPELINE_DIR = Path(__file__).resolve().parent.parent
RND = random.Random(20260914)

RESULTS = []


def _rand(n):
    return "".join(RND.choice("ACGT") for _ in range(n))


def _check(label, cond, detail=""):
    RESULTS.append((label, bool(cond), detail))
    print("[%s] %s%s" % ("PASS" if cond else "FAIL", label,
                         ("  " + detail) if detail and not cond else ""))


# ── 合成序列 ──────────────────────────────────────────────────────────────
# 无 ORF 填充: "TAATAA" 重复在六框里都没有起始密码子
FILLER = "TAATAA" * 30
CDS = "ATG" + "GCT" * 50 + "TAA"                      # 51 aa 完整 ORF
ORF_LEAD = len(FILLER)                                # ORF 起点 = 180
ORF_SEQ = FILLER + CDS + "TAATAA" * 20

P = _rand(30)
M = _rand(300)
DTR_EXACT_SEQ = P + M + P                             # 360 bp, 3' 端重复 30 bp
P_MUT = P[:15] + ("A" if P[15] != "A" else "C") + P[16:]
DTR_FUZZY_SEQ = P + M + P_MUT                         # 3' 端拷贝 1 处错配
DTR_ALT_M = "N" * 40                                  # min_dtr 边界用

I_ = _rand(40)
M_ITR = "A" + _rand(298) + "A"
ITR_SEQ = I_ + M_ITR + tr.reverse_complement(I_)      # ITR 40 bp
PLAIN_SEQ = _rand(500)


def test_detection():
    print("\n── 1/2. DTR 检测 ──")
    # 1. 精确 DTR
    got = tr.fetch_dtr(DTR_EXACT_SEQ)
    _check("fetch_dtr 精确命中 P (30 bp)", got == P, "got=%r" % got[:20])

    # CT3 形状约束: contig 5' 端与 3' 端拷贝必须等长且完全一致
    _check("fetch_dtr 阴性: 无重复随机序列", tr.fetch_dtr(PLAIN_SEQ) == "")
    _check("fetch_dtr 阴性: 序列短于 2*min_dtr",
           tr.fetch_dtr("ACGT" * 5) == "")
    _check("fetch_dtr 阴性: 3' 端拷贝含错配 (exact 口径)",
           tr.fetch_dtr(DTR_FUZZY_SEQ) == "")

    # 2. 容错 DTR
    fuzzy_seq, ident = tr.fetch_dtr_fuzzy(DTR_FUZZY_SEQ, min_identity=0.90)
    _check("fetch_dtr_fuzzy 命中 3' 端拷贝 (30 bp)", len(fuzzy_seq) == 30)
    _check("fetch_dtr_fuzzy 一致率 = 29/30",
           abs(ident - 29 / 30) < 1e-9, "ident=%s" % ident)
    _check("fetch_dtr_fuzzy 阴性: 无重复随机序列",
           tr.fetch_dtr_fuzzy(PLAIN_SEQ, min_identity=0.90) == ("", 0.0))
    _check("fetch_dtr_fuzzy 阈值收紧后落空 (min_identity=0.99)",
           tr.fetch_dtr_fuzzy(DTR_FUZZY_SEQ, min_identity=0.99)[0] == "")

    # 3. ITR
    print("\n── 3. ITR 检测 ──")
    itr = tr.fetch_itr(ITR_SEQ, min_len=20, max_len=1000)
    _check("fetch_itr 命中 40 bp", len(itr) == 40, "len=%d" % len(itr))
    _check("fetch_itr 返回片段 == 5' 端 I", itr == I_)
    _check("fetch_itr 阴性: 无反向重复随机序列",
           tr.fetch_itr(PLAIN_SEQ, min_len=20, max_len=1000) == "")
    _check("fetch_itr 阴性: 5' 端 20 bp 不构成反向重复",
           tr.fetch_itr(_rand(40) + M_ITR + tr.reverse_complement(_rand(40))) == "")


def test_orf_and_rotation():
    print("\n── 7. ORF 扫描与起点归一 ──")
    orfs = tr.find_complete_orfs(ORF_SEQ, min_aa=30)
    _check("find_complete_orfs 只认 1 个完整正向 ORF",
           len(orfs) == 1 and orfs[0][2] == 1,
           "orfs=%r" % (orfs,))
    _check("find_complete_orfs 起点 = %d" % ORF_LEAD,
           orfs and orfs[0][0] == ORF_LEAD)

    off, rc_flag, reason = tr.compute_orf_rotation(ORF_SEQ, min_aa=30)
    _check("compute_orf_rotation 正向锚定 (offset=%d, reason=orf_forward)"
           % ORF_LEAD, off == ORF_LEAD and rc_flag is False
           and reason == "orf_forward", "off=%d rc=%s reason=%s" % (off, rc_flag, reason))
    rotated = tr.rotate_sequence(ORF_SEQ, off)
    _check("旋转后第 1 位是起始密码子", rotated[:3] in tr.START_CODONS)

    # 反向链占优: 把整条序列反向互补, 唯一完整 ORF 落到反向链
    rev_seq = tr.reverse_complement(ORF_SEQ)
    off2, rc2, reason2 = tr.compute_orf_rotation(rev_seq, min_aa=30)
    _check("compute_orf_rotation 反向锚定 (reason=orf_rc)",
           rc2 is True and reason2 == "orf_rc", "rc=%s reason=%s" % (rc2, reason2))
    out2 = tr.rotate_sequence(tr.reverse_complement(rev_seq), off2)
    _check("反向链旋转后同样落在起始密码子", out2[:3] in tr.START_CODONS)

    # 无完整 ORF -> 不旋转
    off3, rc3, reason3 = tr.compute_orf_rotation(FILLER, min_aa=30)
    _check("无完整 ORF 时不旋转 (no_complete_orf)",
           off3 == 0 and rc3 is False and reason3 == "no_complete_orf")
    _check("rotate_sequence 位移 0 为恒等",
           tr.rotate_sequence(PLAIN_SEQ, 0) == PLAIN_SEQ)


def test_cli(tmpdir):
    print("\n── 8. CLI 端到端 ──")
    d = Path(tmpdir)
    fa = d / "in.fasta"
    recs = [
        ("dtr_exact", "len=360 dtr=30", DTR_EXACT_SEQ),
        ("dtr_fuzzy", "len=360 dtr=30 mismatch", DTR_FUZZY_SEQ),
        ("itr_only", "len=380 itr=40", ITR_SEQ),
        ("plain", "len=500 none", PLAIN_SEQ),
        ("user_circ", "len=400 user", _rand(400)),
        ("long_one", "len=6000 too long", _rand(6000)),
    ]
    with open(fa, "w") as fh:
        for cid, desc, seq in recs:
            fh.write(">%s %s\n%s\n" % (cid, desc, seq))
    with open(d / "circ_ids.txt", "w") as fh:
        fh.write("user_circ\n")

    out_dir = d / "out"
    cmd = [sys.executable, str(PIPELINE_DIR / "terminal_repeats.py"),
           "-i", str(fa), "-o", str(out_dir), "-t", "2",
           "--rotation", "orf", "--max-length", "5000",
           "--circ-file", str(d / "circ_ids.txt")]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    _check("CLI 退出码 0", proc.returncode == 0, proc.stderr[-500:])

    f_out = out_dir / "trimmed_contigs.fasta"
    t_out = out_dir / "terminal_repeats_summary.tsv"
    c_out = out_dir / "circular_contigs.txt"
    _check("产出 trimmed_contigs.fasta", f_out.is_file())
    _check("产出 terminal_repeats_summary.tsv", t_out.is_file())
    _check("产出 circular_contigs.txt", c_out.is_file())
    if not (f_out.is_file() and t_out.is_file() and c_out.is_file()):
        return

    rows = {}
    with open(t_out) as fh:
        head = fh.readline().rstrip("\n").split("\t")
        for ln in fh:
            parts = ln.rstrip("\n").split("\t")
            rows[parts[0]] = dict(zip(head, parts))
    _check("表头前 5 列与 CT3 同名同序",
           head[:5] == ["contig", "in_length_contig", "out_length_contig",
                        "dtr_seq", "itr_seq"], str(head[:5]))
    _check("汇总表 6 行", len(rows) == 6, "rows=%d" % len(rows))

    r = rows.get("dtr_exact", {})
    _check("dtr_exact: dtr_length=30 / method=exact / topology=circular_DTR",
           r.get("dtr_length") == "30" and r.get("dtr_method") == "exact"
           and r.get("topology") == "circular_DTR", str(r))
    _check("dtr_exact: 修剪 in-30 == out",
           int(r.get("in_length_contig", 0)) - 30 == int(r.get("out_length_contig", -1)))

    r = rows.get("dtr_fuzzy", {})
    _check("dtr_fuzzy: method=fuzzy / dtr_length=30 / identity 记录",
           r.get("dtr_method") == "fuzzy" and r.get("dtr_length") == "30"
           and r.get("dtr_identity", "").startswith("0.96"), str(r))

    r = rows.get("itr_only", {})
    _check("itr_only: itr_length=40 / 无 DTR 不修剪 / 非环状",
           r.get("itr_length") == "40" and r.get("dtr_seq") == "NA"
           and r.get("topology") == "undetermined"
           and r.get("in_length_contig") == r.get("out_length_contig"), str(r))

    r = rows.get("plain", {})
    _check("plain: 全阴性且长度不变",
           r.get("dtr_seq") == "NA" and r.get("itr_seq") == "NA"
           and r.get("in_length_contig") == r.get("out_length_contig")
           and r.get("rotated") == "False", str(r))

    r = rows.get("user_circ", {})
    _check("user_circ: --circ-file 口径 (不搜 DTR/不修剪)",
           r.get("dtr_seq") == "User-provided circular"
           and r.get("topology") == "circular_user"
           and r.get("in_length_contig") == r.get("out_length_contig"), str(r))

    r = rows.get("long_one", {})
    _check("long_one: 超 --max-length 跳过 DTR",
           r.get("dtr_method") == "too_long" and r.get("dtr_seq") == "NA", str(r))

    # 序列输出核对 (本轮为 --rotation orf)
    out_seqs = {rec.id: str(rec.seq) for rec in tr.SeqIO.parse(str(f_out), "fasta")}
    _check("FASTA 记录数 = 6", len(out_seqs) == 6, "n=%d" % len(out_seqs))
    _check("plain 序列原样透传", out_seqs.get("plain") == PLAIN_SEQ)
    _check("itr_only 序列原样(仅未旋转)",
           out_seqs.get("itr_only") == ITR_SEQ)
    # orf 旋转: 输出 = 对 P+M 做左移 rotation_offset
    off = int(rows["dtr_exact"]["rotation_offset"])
    _check("dtr_exact 修剪 = P+M 且按 rotation_offset 左移",
           out_seqs.get("dtr_exact") == tr.rotate_sequence(P + M, off),
           "offset=%d rotated=%s" % (off, rows["dtr_exact"]["rotated"]))
    _check("dtr_fuzzy ORF 旋转后起点是起始密码子",
           out_seqs.get("dtr_fuzzy", "")[:3] in tr.START_CODONS,
           out_seqs.get("dtr_fuzzy", "")[:3])
    _check("dtr_exact 与 dtr_fuzzy 旋转后互为同一骨架",
           sorted(out_seqs.get("dtr_exact", "")) == sorted(P + M))

    with open(c_out) as fh:
        cids = [ln.strip() for ln in fh if ln.strip()]
    _check("circular_contigs.txt = {dtr_exact, dtr_fuzzy, user_circ}",
           sorted(cids) == ["dtr_exact", "dtr_fuzzy", "user_circ"], str(cids))

    # 二次运行: repeat 模式下 DTR contig 位移为 0
    out2 = d / "out_repeat"
    proc2 = subprocess.run(
        [sys.executable, str(PIPELINE_DIR / "terminal_repeats.py"),
         "-i", str(fa), "-o", str(out2), "-t", "1", "--rotation", "repeat",
         "--max-length", "5000", "--quiet"],
        capture_output=True, text=True)
    _check("CLI --rotation repeat 退出码 0", proc2.returncode == 0,
           proc2.stderr[-300:])
    with open(out2 / "terminal_repeats_summary.tsv") as fh:
        head = fh.readline().rstrip("\n").split("\t")
        rr = {p[0]: dict(zip(head, p)) for p in
              (ln.rstrip("\n").split("\t") for ln in fh)}
    _check("repeat 模式: rotated=False / rotation_offset=0 / rotation_mode=repeat",
           rr["dtr_exact"]["rotated"] == "False"
           and rr["dtr_exact"]["rotation_offset"] == "0"
           and rr["dtr_exact"]["rotation_mode"] == "repeat", str(rr["dtr_exact"]))
    out_seqs2 = {rec.id: str(rec.seq) for rec in
                 tr.SeqIO.parse(str(out2 / "trimmed_contigs.fasta"), "fasta")}
    _check("repeat 模式: dtr_exact 修剪后 = P+M, 起点落在末端重复处",
           out_seqs2.get("dtr_exact") == P + M
           and out_seqs2.get("dtr_exact", "")[:30] == P)
    _check("repeat 模式: dtr_fuzzy 修剪后 = P+M (错配 3' 拷贝被丢弃)",
           out_seqs2.get("dtr_fuzzy") == P + M, out_seqs2.get("dtr_fuzzy", "")[:30])


def main():
    keep = "--keep" in sys.argv
    tmpdir = tempfile.mkdtemp(prefix="tr_selftest_")
    print("=" * 72)
    print("terminal_repeats.py 自测 (临时目录 %s)" % tmpdir)
    print("=" * 72)
    try:
        test_detection()
        test_orf_and_rotation()
        test_cli(tmpdir)
    finally:
        if keep:
            print("\n[keep] 临时产物保留在 %s" % tmpdir)
        else:
            shutil.rmtree(tmpdir, ignore_errors=True)

    passed = sum(1 for _, ok, _ in RESULTS if ok)
    failed = len(RESULTS) - passed
    print("\n" + "=" * 72)
    print("结果: %d passed, %d failed / %d" % (passed, failed, len(RESULTS)))
    print("=" * 72)
    if failed:
        for label, ok, detail in RESULTS:
            if not ok:
                print("  FAIL: %s  %s" % (label, detail))
    return failed == 0


if __name__ == "__main__":
    sys.exit(0 if main() else 1)
