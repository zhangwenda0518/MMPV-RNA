#!/usr/bin/env python3
"""seq_clean.py — 分析前序列清洗（长度 / N / 大小写 / 非法字符 / 元数据缺失）

背景 (2026-09-16, 老师要求「分析前，再次长度检查，去除长度异常，N 碱基异常等，
修正大小写，去除时间和地理都没记录的序列等」)
=====================================================================
现有 `utils/align_qc.py` 已覆盖一部分（长度 / gap% / N% / identity），但实测发现
它有**四个缺口**。实测数据与结论见 `_consistency_check_20260916/verify_seq_clean.py`
（该脚本自带脏数据、可独立重跑）：

  输入 10 条脏序列, align_qc 判定「保留 8 条 / 剔除 2 条(short, withN)」——
  放行的 8 条里, 有 **2 条**是真正有问题的, seq_clean 把它们拦下:

    · star_stop  含 '*' (终止符)      —— align_qc **完全不检查非 IUPAC 字符**
    · digit      含 '9' (数字)        —— 同上

  另外两项是**能力缺口**而非"放行错误", 一并补齐:

    · 大小写:   全小写/混合大小写 —— align_qc 只在内部 .upper() 算指标,
                但**不改写输出**, 小写序列原样进入下游 (本模块统一转大写)
    · 元数据:   align_qc **看不到 metadata**, 无法做时间/地理筛除
    · 长度上限: align_qc 只查"过短", 没有"异常长"的上限检查

  ⚠ 更正 (2026-09-16): 本模块初稿的 docstring 曾写「4 条被放行」, 把
    lower_mix / dot 也算了进去。**实测不符**: lower_mix 属"大小写该修正"
    (修完可用, 本就不是该剔的), dot 见下条注释的口径问题。以实测脚本为准。

  ⚠ 关于 '.' / '-': 本模块是**比对前**, 此时不该有 gap 字符。默认
    `strict_gap_chars=True` 把它们也算异常; 若上游确会输出半比对序列,
    显式设 False 走宽松口径 (见参数文档)。

本模块补齐上述缺口，并明确**两道过滤**的分工：

  第一道 (本模块, 比对前): 纯序列本体
      大小写统一 → 非法字符检测/剔除 → 长度异常 → N 比例
      ↑ 必须在 MAFFT 之前: 空格/星号/小写会让比对程序本身出错或产出错位比对

  第二道 (align_qc, 比对后): 需要比对坐标系的指标
      gap% / 比对坐标系下的 N% / 与参考 identity
      ↑ gap 只有在比对后才存在, 比对前算没有意义

元数据门槛 (两个**独立开关**, 默认都不开):
  --drop-no-date      剔除日期为空的序列
  --drop-no-location  剔除地理为空的序列
  为什么不合成一个「时间地理都缺才剔」的口径:
    下游**各自独立**在自己那一侧处理缺失 ——
      · virphy_bridge.py:754  `dropna(_decimal_year)` 只剔无日期的
      · temporal_signal.py:106 同, 无日期不进随机化池
      · phylogeo_bridge.py:166-172 缺地理的**不剔**, 而是并成 "Other" 组
    所以「系统发育/饱和分析」根本不需要日期, 若入口按"都缺才剔"仍会保留
    只缺日期的样本(无害); 而"缺地理"的样本进系统地理时会变 "Other"(与现状一致)。
    拆成两个开关 = 让使用者按**当前要跑的分析**决定, 而不是替下游统一决定。

用法
----
  # 只看不改 (报告模式, 默认)
  python -m utils.seq_clean --fasta seqs.fa --report clean_report.tsv

  # 输出清洗后的序列 (+ 被剔除的单独存)
  python -m utils.seq_clean --fasta seqs.fa --out clean.fa \
      --max-n 0.05 --min-length-ratio 0.9

  # 开启元数据筛除 (需 metadata)
  python -m utils.seq_clean --fasta seqs.fa --out clean.fa \
      --metadata meta.csv --drop-no-date --drop-no-location

输出
----
  --out 指定的 clean FASTA (仅在给了 --out 时写)
  <out 同目录>/removed.fasta        被剔除序列 (原样保留, 便于复核)
  --report 指定的 TSV (默认 clean_report.tsv) 逐序列判定

安全约定
--------
· **默认只报告, 不写文件**; 必须显式给 --out 才产出 clean.fasta
· 严格判据: 任何 **非 IUPAC 核酸字符** 一律视为异常 (可选 --strip-illegal 剔除而非整条丢弃)
· 长度判据默认用**中位数**做基准 (无参考序列时), 给了 --reference 则用它
· 大小写统一: 一律转**大写** (MAFFT/BEAST/RAxML 通用约定)
"""
from __future__ import annotations

import argparse
import csv
import os
import re
import statistics
import sys
from collections import Counter, OrderedDict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

# 允许 `python utils/seq_clean.py` 与 `python -m utils.seq_clean` 两种跑法
_HERE = os.path.dirname(os.path.abspath(__file__))
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)

# IUPAC 核酸简并码 (含 U: 某些 RNA 病毒数据会带 U)
IUPAC_DNA = set("ACGTRYSWKMBDHVN")
IUPAC_RNA = IUPAC_DNA | {"U"}

# alignment 允许出现的 gap 字符
GAP_CHARS = set("-.")

# 日期列名候选 (与 temporal_signal.py:100 的列名判定保持一致)
DATE_COLS = ("date", "collectiondate", "collection_date", "sampling_date", "time")
LOC_COLS = ("location", "geo_loc_name", "geolocation", "country", "region",
            "state", "locality", "place", "site")
NAME_COLS = ("name", "id", "accession", "seq", "taxon", "isolate", "strain")


# ═══════════════════════════════════════════════════════════════════
# 内部工具
# ═══════════════════════════════════════════════════════════════════

def _read_fasta(path: str) -> "OrderedDict[str, Tuple[str, str]]":
    """读 FASTA → OrderedDict {id: (header_full, seq_raw)}。

    自己解析而不用 SeqIO: 需要看见**原始字符**(空格/星号/小写),
    Biopython 的 SeqIO 对非 IUPAC 字符会给出警告甚至改写。
    也保留 header 全文, 便于下游按 description 匹配。
    """
    out: "OrderedDict[str, Tuple[str, str]]" = OrderedDict()
    header, buf = None, []
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.rstrip("\r\n")
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    out[header[1:].split()[0]] = (header, "".join(buf))
                header, buf = line, []
            else:
                buf.append(line)
    if header is not None:
        out[header[1:].split()[0]] = (header, "".join(buf))
    return out


def _norm_header(s: str) -> str:
    """把元数据里的名字规整成可比较的键 (去引号/空白/大小写)。"""
    return str(s).strip().strip('"').strip("'").strip()


def _load_metadata(path: str) -> Tuple[Dict[str, dict], str, str, str]:
    """读 metadata CSV/TSV → ({name: row_dict}, name_col, date_col, loc_col)。

    列名判定与 temporal_signal.py 一致: 先按候选列名找, 找不到退回首列/次列。

    ⚠ 2026-09-16 修: 原实现三个 `pick` 各自独立回退, 列数不足时会**撞到同一列**
      (实测 2 列表 'a,b' → date_col 与 loc_col 都回退成 'b', 把日期当地理)。
      现在回退时排除已被占用的列; 若真的没有可用列则返回空串, 由调用方报错。
    """
    delim = "\t" if path.lower().endswith((".tsv", ".tab")) else ","
    rows: List[dict] = []
    with open(path, "r", encoding="utf-8", errors="replace", newline="") as f:
        rd = csv.DictReader(f, delimiter=delim)
        if not rd.fieldnames:
            raise SystemExit(f"[seq_clean] metadata 无表头: {path}")
        for r in rd:
            rows.append({(k or "").strip(): (v or "") for k, v in r.items()})
    if not rows:
        return {}, "", "", ""
    cols = list(rows[0].keys())

    def pick(cands, fallback_idx, taken):
        low = {c.lower().replace("-", "_"): c for c in cols}
        for c in cands:
            if c in low and low[c] not in taken:
                return low[c]
        # 回退: 取第 fallback_idx 列, 但跳过已占用的列
        for i, c in enumerate(cols):
            if i >= fallback_idx and c not in taken:
                return c
        return ""

    name_col = pick(NAME_COLS, 0, set())
    date_col = pick(DATE_COLS, 1, {name_col})
    loc_col = pick(LOC_COLS, 2, {name_col, date_col})
    return {_norm_header(r.get(name_col, "")): r for r in rows if r.get(name_col) is not None}, \
        name_col, date_col, loc_col


def _is_present(v: str) -> bool:
    """元数据字段是否算"有值"。空串 / NA / N/A / unknown / missing / - 都算缺。"""
    s = str(v or "").strip()
    return bool(s) and s.lower() not in (
        "na", "n/a", "nan", "none", "unknown", "missing", "not_provided",
        "not provided", "-", "--", ".", "?")


# ═══════════════════════════════════════════════════════════════════
# 主逻辑
# ═══════════════════════════════════════════════════════════════════

def clean_sequences(
    fasta_path: str,
    *,
    reference: Optional[str] = None,
    min_length_ratio: float = 0.9,
    max_n: float = 0.05,
    max_length_ratio: float = 1.5,
    min_abs_length: int = 30,
    allow_rna: bool = False,
    strip_illegal: bool = False,
    metadata_path: Optional[str] = None,
    drop_no_date: bool = False,
    drop_no_location: bool = False,
    strict_gap_chars: bool = True,
    log=None,
) -> dict:
    """分析前序列清洗（比对前第一道）。

    Parameters
    ----------
    fasta_path : str
    reference : str, optional
        参考序列 ID（子串匹配）。给了就用它的**去 gap 长度**作基准；
        否则用全部序列长度中位数。中位数更稳健（少数超短/超长不拉动基准）。
    min_length_ratio : float, default 0.9
        长度 / 基准 的下限。低于则剔。
    max_length_ratio : float, default 1.5
        长度 / 基准 的上限，抓**异常长的嵌合/重复序列**（旧版没有上限检查）。
    max_n : float, default 0.05
        N（及简并码总数）占比上限。
    min_abs_length : int, default 30
        绝对长度下限。比值判据在"全体都短"时会全部放行，这道兜底。
    allow_rna : bool, default False
        True 则把 U 视为合法（RNA 病毒）。
    strip_illegal : bool, default False
        True → 删除非法字符后继续用该序列（并在报告标注 stripped 数）；
        False（默认）→ 整条剔除。**默认保守**：非法字符往往是数据问题的信号，
        静默删掉可能造出一条看似正常、实则错位的序列。
    strict_gap_chars : bool, default True
        **本模块是比对前**，此时不该存在 gap 字符。默认 True 表示:
          `-` 与 `.` 出现在**未比对**序列里 = 异常, 整条剔除。
        设 False 则沿用宽松口径（把它们当 gap 从长度里减掉、其余照常）。
        取舍: 真实未比对 FASTA 里出现 `-`/`.` 几乎都是数据损坏
        (手工编辑残留、跨库拼接、从比对结果误取),
        但也可能来自上游确实输出了半比对序列; 若要吃这类输入请显式设 False。
    metadata_path / drop_no_date / drop_no_location : 元数据门槛，见模块 docstring。
    log : callable, optional

    Returns
    -------
    dict: records(保留), removed(剔除), report(逐序列), stats(汇总),
          issues(理由计数), reference_len, basis
    """
    def emit(m):
        if log is None:
            print(m)
        elif callable(log):
            log(m)
        elif hasattr(log, "emit"):     # utils/*_bridge.py 的 LogCollector
            log.emit(m)
        else:
            print(m)

    recs = _read_fasta(fasta_path)
    if not recs:
        raise SystemExit(f"[seq_clean] 空 FASTA: {fasta_path}")
    emit(f"[seq_clean] 读入 {len(recs)} 条序列")

    # ── 基准长度 ──
    bare_lens = {k: len(re.sub(r"\s", "", s)) for k, (_, s) in recs.items()}
    ref_id = None
    if reference:
        ref_id = next((k for k in recs if reference in k), None)
    if ref_id:
        basis = bare_lens[ref_id]
        basis_src = f"参考 {ref_id}"
    else:
        if reference:
            emit(f"[seq_clean] ⚠ 未找到参考 '{reference}', 改用中位数长度")
        basis = int(statistics.median(bare_lens.values()))
        basis_src = "中位数"
    emit(f"[seq_clean] 长度基准: {basis} bp ({basis_src})")

    # ── 元数据 ──
    meta: Dict[str, dict] = {}
    m_name = m_date = m_loc = ""
    if metadata_path:
        if not os.path.exists(metadata_path):
            raise SystemExit(f"[seq_clean] metadata 不存在: {metadata_path}")
        meta, m_name, m_date, m_loc = _load_metadata(metadata_path)
        emit(f"[seq_clean] metadata: {len(meta)} 行, "
             f"name={m_name!r} date={m_date!r} location={m_loc!r}")
        if drop_no_date or drop_no_location:
            emit(f"[seq_clean] 元数据门槛: drop_no_date={drop_no_date} "
                 f"drop_no_location={drop_no_location}")
            # 2026-09-16 护栏: 开了门槛但对应列没识别出来时, 逐行判定会拿
            # row.get("") == "" → has_date/has_loc=False → **整批序列被判"无日期/无地理"
            # 而全部删除** (实测: 两列 'name,date' 的表配 --drop-no-location 必清空)。
            # 这属于"静默批量删数据", 必须 fail-loud 而不是照删。
            if drop_no_date and not m_date:
                raise SystemExit("[seq_clean] metadata 里没有可识别的**日期**列, "
                                 "拒绝执行 --drop-no-date (避免整批误删)")
            if drop_no_location and not m_loc:
                raise SystemExit("[seq_clean] metadata 里没有可识别的**地名**列, "
                                 "拒绝执行 --drop-no-location (避免整批误删)")
    elif drop_no_date or drop_no_location:
        raise SystemExit("[seq_clean] 需要 --metadata 才能用 --drop-no-date/-location")

    legal = IUPAC_RNA if allow_rna else IUPAC_DNA
    # gap 字符是否算合法: 比对前不该有 gap → 默认视为异常
    gap_is_illegal = strict_gap_chars
    report: List[dict] = []
    keep: List[dict] = []
    removed: List[dict] = []
    issues = Counter()
    n_stripped_total = 0
    n_case_fixed = 0

    for sid, (header, raw) in recs.items():
        issues_this: List[str] = []

        # ① 大小写: 统一转大写 (不改内容, 只改 case)
        seq = raw.upper()
        had_lower = any(c.islower() for c in raw)
        if had_lower:
            n_case_fixed += 1

        # ② 非法字符检测
        #    先去掉空白 (FASTA 行内空格/制表符属排版, 不算序列内容)
        no_ws = re.sub(r"\s", "", seq)
        #    再区分: gap 字符 vs 真正的非法字符
        #    strict_gap_chars=True(默认) 时, '-'/'.' 在**未比对**序列里也算异常
        allowed = legal if gap_is_illegal else (legal | GAP_CHARS)
        illegal = sorted({c for c in no_ws if c not in allowed})
        n_illegal = sum(no_ws.count(c) for c in illegal)
        if illegal:
            issues_this.append(f"非法字符 {'/'.join(illegal)} x{n_illegal}")

        work = no_ws
        stripped = 0
        if illegal and strip_illegal:
            for c in illegal:
                stripped += work.count(c)
                work = work.replace(c, "")
            n_stripped_total += stripped

        if not work:
            issues_this.append("序列为空")

        # ③ 长度 (去 gap 后的裸长度)
        bare = work.replace("-", "").replace(".", "")
        L = len(bare)
        len_ratio = L / basis if basis else 0.0

        # ④ N / 简并码比例
        n_cnt = sum(bare.count(c) for c in "N")
        degen_cnt = sum(bare.count(c) for c in "RYSWKMBDHV")
        n_ratio = (n_cnt + degen_cnt) / L if L else 1.0

        # ⑤ 元数据
        row = meta.get(_norm_header(sid))
        has_date = has_loc = None
        if metadata_path:
            if row is None:
                has_date = has_loc = False
                issues_this.append("metadata 无此行")
            else:
                has_date = _is_present(row.get(m_date, ""))
                has_loc = _is_present(row.get(m_loc, ""))

        # ── 判定 ──
        if not illegal and L:
            if L < min_abs_length:
                issues_this.append(f"绝对长度过短 ({L} < {min_abs_length})")
            if len_ratio < min_length_ratio:
                issues_this.append(f"长度不足 ({len_ratio:.2f} < {min_length_ratio})")
            elif len_ratio > max_length_ratio:
                issues_this.append(f"长度异常长 ({len_ratio:.2f} > {max_length_ratio})")
        if n_ratio > max_n:
            issues_this.append(f"N/简并过多 ({n_ratio:.2f} > {max_n})")
        if drop_no_date and has_date is False:
            issues_this.append("无日期")
        if drop_no_location and has_loc is False:
            issues_this.append("无地理")

        status = "REMOVE" if issues_this else "KEEP"
        for i in issues_this:
            issues[i.split(" ")[0]] += 1

        rec = {"id": sid, "header": header, "seq": work,
               "raw_len": len(no_ws), "bare_len": L}
        (keep if status == "KEEP" else removed).append(rec)

        report.append({
            "id": sid,
            "raw_len": len(no_ws),
            "bare_len": L,
            "len_ratio": round(len_ratio, 4),
            "n_degen_ratio": round(n_ratio, 4),
            "illegal_chars": "".join(illegal) if illegal else "",
            "n_illegal": n_illegal,
            "case_fixed": "Y" if had_lower else "",
            "has_date": "" if has_date is None else ("Y" if has_date else "N"),
            "has_location": "" if has_loc is None else ("Y" if has_loc else "N"),
            "status": status,
            "reason": "; ".join(issues_this),
        })

    stats = {
        "n_input": len(recs),
        "n_keep": len(keep),
        "n_removed": len(removed),
        "basis": basis,
        "basis_src": basis_src,
        "reference_id": ref_id,
        "n_case_fixed": n_case_fixed,
        "n_stripped_chars": n_stripped_total,
    }
    emit(f"[seq_clean] 保留 {len(keep)} / 剔除 {len(removed)} "
         f"(共 {len(recs)}); 大小写修正 {n_case_fixed} 条")
    if issues:
        emit(f"[seq_clean] 剔除理由: {dict(issues)}")

    return {"records": keep, "removed": removed, "report": report,
            "stats": stats, "issues": dict(issues)}


def write_outputs(result: dict, out_fasta: Optional[str], report_path: str,
                  wrap: int = 60) -> dict:
    """写 clean.fasta / removed.fasta / 报告 TSV。

    **只有给了 out_fasta 才写 FASTA** —— 没给则纯报告模式, 绝不产生文件。
    """
    produced = {}
    rep = result["report"]
    if rep:
        with open(report_path, "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=list(rep[0].keys()), delimiter="\t")
            w.writeheader()
            w.writerows(rep)
        produced["report"] = report_path

    if not out_fasta:
        return produced

    def _dump(recs, path):
        with open(path, "w", encoding="utf-8", newline="\n") as f:
            for r in recs:
                f.write(f">{r['id']}\n")
                s = r["seq"]
                for i in range(0, len(s), wrap):
                    f.write(s[i:i + wrap] + "\n")

    os.makedirs(os.path.dirname(os.path.abspath(out_fasta)) or ".", exist_ok=True)
    _dump(result["records"], out_fasta)
    produced["clean"] = out_fasta

    # removed 放在 clean 同目录, 便于一并复核
    rem_path = os.path.join(os.path.dirname(os.path.abspath(out_fasta)), "removed.fasta")
    _dump(result["removed"], rem_path)
    produced["removed"] = rem_path
    return produced


def main(argv=None):
    ap = argparse.ArgumentParser(
        description="分析前序列清洗: 长度 / N / 大小写 / 非法字符 / 元数据缺失",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="默认只出报告不改数据; 加 --out 才写 clean.fasta(与 removed.fasta)")
    ap.add_argument("--fasta", required=True, help="输入 FASTA")
    ap.add_argument("--out", default=None, help="输出 clean FASTA (不给则只报告)")
    ap.add_argument("--report", default="clean_report.tsv", help="判定报告 TSV")
    ap.add_argument("--reference", default=None, help="参考序列 ID (长度基准)")
    ap.add_argument("--min-length-ratio", type=float, default=0.9)
    ap.add_argument("--max-length-ratio", type=float, default=1.5)
    ap.add_argument("--min-abs-length", type=int, default=30)
    ap.add_argument("--max-n", type=float, default=0.05)
    ap.add_argument("--allow-rna", action="store_true", help="允许 U (RNA 病毒)")
    ap.add_argument("--strip-illegal", action="store_true",
                    help="删除非法字符而非剔除整条 (默认整条剔除, 更保守)")
    ap.add_argument("--metadata", default=None, help="metadata CSV/TSV")
    ap.add_argument("--drop-no-date", action="store_true", help="剔除无日期序列")
    ap.add_argument("--drop-no-location", action="store_true", help="剔除无地理序列")
    ap.add_argument("--allow-gap-chars", action="store_true",
                    help="宽松: 允许未比对序列里的 '-'/'.' (默认视为异常剔除)")
    args = ap.parse_args(argv)

    r = clean_sequences(
        args.fasta,
        reference=args.reference,
        min_length_ratio=args.min_length_ratio,
        max_length_ratio=args.max_length_ratio,
        min_abs_length=args.min_abs_length,
        max_n=args.max_n,
        allow_rna=args.allow_rna,
        strip_illegal=args.strip_illegal,
        metadata_path=args.metadata,
        drop_no_date=args.drop_no_date,
        drop_no_location=args.drop_no_location,
        strict_gap_chars=not args.allow_gap_chars,
    )
    produced = write_outputs(r, args.out, args.report)

    if r["removed"]:
        print("\n被剔除的序列:")
        for row in r["report"]:
            if row["status"] == "REMOVE":
                print(f"  x {row['id']}: {row['bare_len']}bp "
                      f"len={row['len_ratio']:.2f} N={row['n_degen_ratio']:.2f} "
                      f"→ {row['reason']}")
    print("\n输出:")
    for k, v in produced.items():
        print(f"  {k}: {v}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
