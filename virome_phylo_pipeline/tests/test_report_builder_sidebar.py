# -*- coding: utf-8 -*-
"""⑤ HTML 报告侧栏回归测试 (2026-10-17, 用户两次投诉侧栏消失后的护栏)。

覆盖三件事:
  1. 侧栏存在, 且**任何窗口宽度都不隐藏** —— 禁止在 @media(max-width:…) 断点里
     写 `.sidebar{display:none}` (③ 曾因断点隐藏被用户投诉, ⑤ 不许复现);
     仅 @media print 允许隐藏。
  2. 顶部 sticky chips (#8 的 .topnav) 已移除, 不许回流。
  3. 图片内联上限 = 1200 KB (300dpi PNG 能进内联分支), PDF 上限仍 6000 KB。
"""
import re
import sys
from pathlib import Path

import pytest

_HERE = Path(__file__).resolve().parent
_PHYLO = _HERE.parent
_REPO = _PHYLO.parent
for p in (str(_PHYLO), str(_REPO)):
    if p not in sys.path:
        sys.path.insert(0, p)

# ⑤ 的 utils/ 与其他管线的同名包冲突: 全量 pytest 时清缓存再把 ⑤ 插到最前
for _m in [m for m in list(sys.modules) if m == "utils" or m.startswith("utils.")]:
    del sys.modules[_m]
sys.path.insert(0, str(_PHYLO))

from utils.report_builder import (  # noqa: E402
    MAX_IMG_KB, MAX_PDF_KB, build_html_report, embed_asset,
)


def _media_blocks(css: str):
    """切出 [(条件, 块体)] —— 手工扫花括号配对, 不引入额外依赖。

    先剥掉 /* … */ 注释: 注释里可能出现 `@media(max-width:…)` 之类的说明文字
    (侧栏裁定注释就写了), 不剥会把注释误解析成一个断点块。
    """
    css = re.sub(r"/\*.*?\*/", "", css, flags=re.S)
    blocks, i = [], 0
    while True:
        m = re.search(r"@media([^{]*)\{", css[i:])
        if not m:
            break
        brace = i + m.end() - 1
        depth = 0
        for j in range(brace, len(css)):
            if css[j] == "{":
                depth += 1
            elif css[j] == "}":
                depth -= 1
                if depth == 0:
                    blocks.append((m.group(1).strip(), css[brace + 1:j]))
                    i = j + 1
                    break
        else:
            break
    return blocks


@pytest.fixture(scope="module")
def report_html(tmp_path_factory):
    """最小夹具: 1 病毒 + 空 work_dir, 走真实入口产出一份报告 HTML。"""
    tmp = tmp_path_factory.mktemp("phylo_report_sb")
    vdir = tmp / "Virus_A"
    vdir.mkdir()
    (vdir / ".checkpoints").mkdir()
    out = tmp / "out"
    out.mkdir()
    results = [{
        "virus": "Virus A", "phase": "batch", "success": True, "metrics": {},
        "virus_dir": str(vdir), "stages_done": ["tree"], "stages_skipped": [],
        "stages_failed": [], "stages_blocked": [],
    }]
    path = build_html_report(results, str(out))
    return Path(path).read_text(encoding="utf-8")


def test_sidebar_present(report_html):
    """侧栏本体 + 8 模块条目 + Viruses 分组 (③ 同款 sb 体系)。"""
    assert '<div class="sidebar">' in report_html
    assert '<div class="sb-title">Pipeline Modules</div>' in report_html
    assert '<div class="sb-title" style="margin-top:10px;">Viruses</div>' in report_html
    assert report_html.count('class="sb-dot"') >= 8          # 8 模块各一 dot
    assert 'class="sb-item sb-' in report_html
    assert 'id="virus-0"' in report_html                     # 病毒锚点仍可用
    for key in ("data", "phylogeny", "popgen", "recomb",
                "select", "time", "geography", "report"):
        assert f'href="#mod-{key}"' in report_html


def test_sidebar_never_hidden_by_width(report_html):
    """⚠ 核心裁定: 侧栏常驻 —— 任何 max-width 断点都不得隐藏它。"""
    css = report_html.split("<style>", 1)[1].split("</style>", 1)[0]
    # 基础规则: fixed 定位 + 224px 宽 (与 ③ 同款)
    assert re.search(r"\.sidebar\s*\{[^}]*position:\s*fixed", css)
    assert re.search(r"\.sidebar\s*\{[^}]*width:\s*224px", css)
    assert re.search(r"\.container\s*\{[^}]*padding-left:\s*240px", css)
    # 宽度断点里不得出现 .sidebar (尤其 display:none)
    for cond, body in _media_blocks(css):
        if "max-width" in cond or "min-width" in cond:
            assert ".sidebar" not in body, (
                f"@media({cond}) 里出现了 .sidebar —— 侧栏必须常驻: {body!r}")
    # 打印时允许隐藏
    print_blocks = [b for c, b in _media_blocks(css) if "print" in c]
    assert print_blocks and any(".sidebar" in b and "display: none" in b
                                for b in print_blocks)


def test_topnav_removed(report_html):
    """顶部 sticky chips (#8 的 .topnav / .chip / .nav-select) 已整体移除。"""
    assert "topnav" not in report_html
    assert 'class="chip"' not in report_html
    assert "nav-select" not in report_html


def test_scrollspy_present(report_html):
    """scroll-spy 脚本在 (侧栏高亮当前 section)。"""
    assert "sb-active" in report_html
    assert "addEventListener('scroll'" in report_html
    assert ".virus-card[id], .mod[id]" in report_html


def test_image_inline_cap_1200kb(tmp_path):
    """图片内联上限 1200 KB: 900 KB 内联, 1300 KB 落排除清单; PDF 上限不变。"""
    assert MAX_IMG_KB == 1200
    assert MAX_PDF_KB == 6000
    ok = tmp_path / "ok.png"
    ok.write_bytes(b"\x89PNG\r\n\x1a\n" + b"0" * (900 * 1024))
    kind, uri, reason = embed_asset(str(ok))
    assert kind == "img" and uri.startswith("data:image/png;base64,") and reason == ""
    too_big = tmp_path / "too_big.png"
    too_big.write_bytes(b"\x89PNG\r\n\x1a\n" + b"0" * (1300 * 1024))
    kind2, uri2, reason2 = embed_asset(str(too_big))
    assert kind2 == "img" and uri2 == "" and "inline limit" in reason2
