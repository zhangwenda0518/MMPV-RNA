#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""eve_final_report.py — 03_report 自适应最终 HTML.

为什么"自适应": 这条管线的输入是两条互斥路径 ——
  * 有宿主参考基因组: 有基因组侧翼/剪接等**整合证据**, 走 with_ref_genome
  * 只有转录组: 只有序列判别 (ML / Δbit-score / 无监督), 没有整合证据
上游还可能只跑到一半 (例如只有 ML 分, 没算 Δbit)。所以最终报告不能写死章节,
要先看**哪些通道真的产出了表**, 再决定输出什么:

  有 final_tiers.tsv  -> 证据合并报告 (tier 分布 + 证据组合 + 明细)
  没有 final_tiers    -> 明确写出缺哪几路、当前结论只基于哪些通道, 再给通道汇总

用法:
  python eve_final_report.py --out out/03_Report/eve_final_report.html \
      --final-tiers out/03_Report/final_tiers.tsv \
      --channel "ML=out/02_Classify/ML/ml_scores.tsv" \
      --channel "Δbit=out/02_Classify/DeltaBit/dbitscore.tsv"
"""

import argparse
import sys
from pathlib import Path

_HERE = Path(__file__).resolve().parent
_PIPE = _HERE.parent
for _p in (_PIPE, _HERE):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

import eve_html as EH   # noqa: E402

# (key, 展示名, 证据方向). key 与 eve_html.kind_of 的表型名一致,
# 用来判断某个通道有没有真的给表.
EXPECTED = [
    ("ref_classify", "双库分类 (EVE 元件库 vs 真病毒库)", "内源/外源双向"),
    ("chimera", "转录本嵌合 (宿主-病毒-宿主)", "内源"),
    ("locus", "位点架构 (组件齐全)", "内源"),
    ("refgenome", "参考基因组整合 (侧翼/剪接)", "内源·金标准"),
    ("ml", "ML 序列打分 (16 科)", "序列侧"),
    ("dbitscore", "Δbit-score 双库竞争", "内源/外源双向"),
    ("nanovir", "Nanoviridae 无监督降解", "内源"),
]

# 哪些通道属于"整合证据"; 缺它们时报告必须写明结论强度受限
INTEGRATION = {"chimera", "locus", "refgenome"}


def parse_kv(s):
    if "=" in s:
        name, path = s.split("=", 1)
        return name.strip(), path.strip()
    return Path(s).stem, s.strip()


def main(argv=None):
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--final-tiers", default="",
                    help="03_report/eve_final_score.py 产出的 final_tiers.tsv")
    ap.add_argument("--channel", action="append", metavar="NAME=PATH",
                    help="某一路证据表 (与 02_classify 产物同格式), 可重复")
    ap.add_argument("--out", required=True)
    ap.add_argument("--title", default="EVE 判别 · 03_report 最终汇总")
    ap.add_argument("--notes", default="")
    ap.add_argument("--top", type=int, default=30)
    a = ap.parse_args(argv)

    channels = [parse_kv(s) for s in (a.channel or [])]
    seen_kinds, channel_rows, sections = set(), {}, []
    used_rows, missing_files = [], []
    for name, path in channels:
        if not Path(path).is_file():
            missing_files.append(str(path))
            used_rows.append({"table": name, "type": "-", "rows": "-",
                              "path": str(path) + "  (文件缺失)"})
            continue
        headers, rows = EH.read_table(path)
        kind = EH.kind_of(headers)
        seen_kinds.add(kind)
        channel_rows[kind] = len(rows)
        used_rows.append({
            "table": name,
            "type": EH.KIND_INFO.get(kind, EH.KIND_INFO["generic"])[0],
            "rows": f"{len(rows):,}",
            "path": str(path),
        })

    # 有 final_tiers 才是"完整合并"; 否则按可用通道出结论
    ft_headers, ft_rows = EH.read_table(a.final_tiers) if a.final_tiers else ([], [])
    has_final = bool(ft_rows)

    status_rows = []
    for key, label, direction in EXPECTED:
        provided = key in seen_kinds
        status_rows.append({
            "channel": label,
            "direction": direction,
            "status": "已提供" if provided else "未提供",
            "rows": f"{channel_rows.get(key, 0):,}" if provided else "-",
        })

    if has_final:
        sections.append(EH.section(
            "final_tiers", "多路证据合并 (final_tiers)", ft_rows, top=a.top,
            extra_note="tier1=≥2 条独立内源证据; tier2=1 条; "
                       "tier3=无证据; tier4=有外源证据且无内源证据。"))
    else:
        missing = [label for key, label, _d in EXPECTED if key not in seen_kinds]
        integ_got = [label for key, label, _d in EXPECTED
                     if key in seen_kinds and key in INTEGRATION]
        if not channels:
            diag = ("本阶段没有收到任何证据表, 也没有 final_tiers。"
                    "请先运行 02_classify 或直接给出 --channel。")
        elif not integ_got:
            diag = ("当前**没有任何整合证据** (参考基因组侧翼 / 位点架构 / 嵌合), "
                    "结论只基于序列判别 —— 这对无参考基因组的物种是预期结果, "
                    "不等同于已确认整合。缺失通道: " + "、".join(missing) + "。")
        else:
            diag = ("已有整合证据通道: " + "、".join(integ_got)
                    + "; 但未生成 final_tiers (缺少 ref_classify 等合并输入或未运行 "
                      "eve_final_score.py), 因此不输出 tier 分层。缺失通道: "
                    + "、".join(missing) + "。")
        sections.append('<div class="note">' + diag + "</div>")
        # 没有合并结果时, 把已有通道的汇总直接放进最终报告, 免得用户再翻 02 的报告
        for name, path in channels:
            if not Path(path).is_file():
                continue
            headers, rows = EH.read_table(path)
            kind = EH.kind_of(headers)
            sections.append(EH.section(
                kind, f"{name} · {EH.KIND_INFO.get(kind, EH.KIND_INFO['generic'])[0]}",
                rows, top=a.top,
                extra_note=EH.KIND_INFO.get(kind, EH.KIND_INFO["generic"])[1]))

    body = ["<h2>证据通道可用性</h2>"
            + EH.table_html(["channel", "direction", "status", "rows"],
                            status_rows, limit=50, num_cols=("rows",)),
            "<h2>本次实际使用的产物</h2>"]
    body.append(EH.table_html(["table", "type", "rows", "path"], used_rows,
                              limit=100, num_cols=("rows",))
                if used_rows else '<div class="meta">无</div>')

    if missing_files:
        body.insert(0, '<div class="note">以下证据表路径不存在, 已跳过: '
                    + "、".join(missing_files) + "</div>")
    if a.notes:
        body.insert(0, f'<div class="note">{a.notes}</div>')

    html = EH.render_page(
        a.title,
        "按实际可用通道自适应输出: 有整合证据时给 tier 分层, 只有序列判别时明确标注结论边界",
        body + sections,
        footer="03_report/eve_final_report.py")
    out = EH.write_page(a.out, html)
    print(f"[done] final_report -> {out}  "
          f"({'含 final_tiers' if has_final else '无 final_tiers, 按可用通道汇总'})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
