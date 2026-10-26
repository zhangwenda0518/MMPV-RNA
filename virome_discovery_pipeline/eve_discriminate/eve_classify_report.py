#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""eve_classify_report.py — 02_classify 阶段 HTML 汇总生成器.

把这一阶段任意组合的产物汇总成一页自包含 HTML:
  * 无基因组路径: ML 打分 / Δbit-score / Nanoviridae 无监督 / 双库分类
  * 有基因组路径: 参考基因组整合证据 / 位点架构 / 转录本嵌合
  * 两张路径的产物可以同时给, 章节按传入顺序排列

为什么表用 NAME=PATH 传: 阶段产物常常改文件名或挪目录, 报告不绑定固定路径;
表型由表头自动识别 (见 eve_html.kind_of), 名字只用于章节标题。

用法:
  python eve_classify_report.py --out out/02_Classify/classify_report.html \
      --table "ML 打分=scores/ml_scores.tsv" \
      --table "Δbit-score=dbitscore.tsv" \
      --table "参考基因组证据=refgenome_evidence.tsv"

  # 只给一个目录时自动收 *.tsv (文件名作章节名)
  python eve_classify_report.py --dir out/02_Classify --out out/02_Classify/classify_report.html
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


def parse_table_arg(s):
    if "=" in s:
        name, path = s.split("=", 1)
        return name.strip(), path.strip()
    return Path(s).stem, s.strip()


def collect(args):
    entries = []
    for s in args.table or []:
        entries.append(parse_table_arg(s))
    if args.dir:
        for p in sorted(Path(args.dir).glob("*.tsv")):
            entries.append((p.stem, str(p)))
        for p in sorted(Path(args.dir).glob("*.tsv.gz")):
            entries.append((p.name, str(p)))
    # 去重 (同名同路径只留一次)
    seen, out = set(), []
    for name, path in entries:
        key = (name, str(Path(path).resolve()))
        if key in seen:
            continue
        seen.add(key)
        out.append((name, path))
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--table", action="append", metavar="NAME=PATH",
                    help="阶段产物表, 可重复")
    ap.add_argument("--dir", default="",
                    help="从目录自动收 *.tsv / *.tsv.gz (与 --table 可叠加)")
    ap.add_argument("--out", required=True, help="输出 HTML")
    ap.add_argument("--title", default="EVE 候选判别 · 02_classify 阶段汇总")
    ap.add_argument("--subtitle", default="")
    ap.add_argument("--notes", default="",
                    help="追加到页首的说明 (通常由编排器写入输入来源)")
    ap.add_argument("--top", type=int, default=25, help="每张表明细最多显示行数")
    a = ap.parse_args(argv)

    tables = collect(a)
    if not tables:
        print("[error] 没有可汇总的表: 用 --table 或 --dir 指定", file=sys.stderr)
        return 2

    sections, entries = [], []
    total_rows = 0
    for name, path in tables:
        headers, rows = EH.read_table(path)
        kind = EH.kind_of(headers)
        total_rows = max(total_rows, len(rows))
        entries.append((name, path, kind, len(rows)))
        title = f"{name} · {EH.KIND_INFO.get(kind, EH.KIND_INFO['generic'])[0]}"
        note = EH.KIND_INFO.get(kind, EH.KIND_INFO["generic"])[1]
        if not rows:
            note += " —— 表为空 (上游可能跳过了该通道)"
        sections.append(EH.section(kind, title, rows, top=a.top,
                                   extra_note=note))

    head = EH.cards([
        ("表数", len(tables)),
        ("最大行数", f"{total_rows:,}"),
        ("识别的表型", len({e[2] for e in entries})),
    ])
    sections.insert(0, "<h2>产物概览</h2>" + head
                    + EH.overview_table(entries))
    if a.notes:
        sections.insert(0, f'<div class="note">{a.notes}</div>')

    subtitle = a.subtitle or ("有参考基因组路径 (整合证据) 与无参考基因组路径 "
                              "(序列判别) 的产物可在本页并列查看")
    html = EH.render_page(a.title, subtitle, sections,
                          footer="02_classify/eve_classify_report.py")
    out = EH.write_page(a.out, html)
    print(f"[done] classify_report -> {out}  ({len(tables)} 张表)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
