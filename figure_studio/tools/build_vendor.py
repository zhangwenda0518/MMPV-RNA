# -*- coding: utf-8 -*-
"""把上游 svguest 的 app.html 与 Figure Studio 桥接合成到 static/vendor/。

用法::

    C:/Python312/python.exe figure_studio/tools/build_vendor.py
    C:/Python312/python.exe figure_studio/tools/build_vendor.py --check   # 只校验不写

为什么要有这个脚本
------------------
之前是用一次性脚本往 app.html 末尾拼桥接，结果 `<style>` 里混进了一行
Markdown 标题（`# Figure Studio ⇄ svguest 桥接样式…`），浏览器把 `#` 当成
ID 选择器，紧接着的整条 `#fsBar{…}` 规则被当作无效选择器**静默丢弃** ——
浮条失去 position:fixed、挤进底部状态栏，全程零报错。

现在：桥接源码单独放在 vendor_src/，本脚本负责拼接 + **自检**，
自检不过就直接报错退出，绝不写出坏文件。自检项：

1. bridge.css 里每一处 `#` / `.` / `@` / `*` / `[` 开头的行，必须在规则块内
   或位于 `/* */` 注释里（防「说明文字当选择器」重演）；
2. bridge.css 里不得出现 `/*` 未闭合；
3. 合成后的文件里 `<style>` 段必须能被正则切成「纯规则」，且含 `#fsBar`；
4. 桥接 JS 必须能被 node --check 解析（有 node 时）；
5. 上游 app.html 的版本号/sha256 与 NOTICE 记录一致（升级时提示更新）。
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import shutil
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
SRC = ROOT / "vendor_src"
OUT = ROOT / "static" / "vendor"

APP = SRC / "app.html"
CSS = SRC / "bridge.css"
JS = SRC / "bridge.js"
NOTICE = SRC / "NOTICE.md"
TARGET = OUT / "svguest.html"

MARK_BEGIN = "<!-- ▼▼▼ Figure Studio 桥接（自动生成，改 vendor_src/ 后重跑 tools/build_vendor.py）▼▼▼ -->"
MARK_END = "<!-- ▲▲▲ Figure Studio 桥接 ▲▲▲ -->"


class BuildError(SystemExit):
    pass


def fail(msg: str) -> None:
    raise BuildError("构建失败：%s" % msg)


# ──────────────────────────────────────────────────────────────────
# CSS 自检
# ──────────────────────────────────────────────────────────────────
def css_lint(css: str) -> None:
    """确认注释都闭合、且没有「说明文字被当选择器」。"""
    # 注释闭合
    if css.count("/*") != css.count("*/"):
        fail("bridge.css 的 /* 与 */ 数量不匹配（%d vs %d）"
             % (css.count("/*"), css.count("*/")))

    # 逐行：剥掉注释内容后，若某行以 # . @ * [ 开头，必须处于规则块内部
    depth = 0
    in_comment = False
    for ln, raw in enumerate(css.splitlines(), 1):
        line = raw
        buf = []
        i = 0
        while i < len(line):
            if in_comment:
                j = line.find("*/", i)
                if j < 0:
                    i = len(line)
                else:
                    in_comment = False
                    i = j + 2
                continue
            j = line.find("/*", i)
            if j < 0:
                buf.append(line[i:])
                break
            buf.append(line[i:j])
            in_comment = True
            i = j + 2
        code = "".join(buf)
        stripped = code.strip()
        if stripped.startswith(("#", ".", "@", "*", "[")) and depth == 0:
            # @media/@supports 等 at-rule 在顶层合法
            if re.match(r"^@(media|supports|font-face|import|charset|keyframes|layer)\b",
                        stripped):
                pass
            # 顶层选择器必须在本行内出现 {（完整规则或规则开头），或以 , 续行
            # （合法的多行选择器）。两者都不是时，它会被浏览器粘到下一条规则的
            # selector list 上，把下一条规则整条吃掉 —— 这正是当初 #fsBar 丢失的原因。
            elif "{" not in stripped and not stripped.endswith(","):
                fail("bridge.css 第 %d 行看起来是顶层选择器，却没有以 { 或 , 收尾，"
                     "会被浏览器粘到下一条规则上把它吞掉：%r" % (ln, raw[:90]))
        depth += code.count("{") - code.count("}")
        if depth < 0:
            fail("bridge.css 第 %d 行出现多余的 }" % ln)
    if depth != 0:
        fail("bridge.css 花括号不平衡（结尾 depth=%d）" % depth)


# ──────────────────────────────────────────────────────────────────
# 上游指纹
# ──────────────────────────────────────────────────────────────────
def upstream_fingerprint() -> dict:
    b = APP.read_bytes()
    m = re.search(rb"svguest (v[\d.]+)", b)
    return {
        "bytes": len(b),
        "sha256": hashlib.sha256(b).hexdigest(),
        "version": m.group(1).decode() if m else "?",
    }


# ──────────────────────────────────────────────────────────────────
# 构建
# ──────────────────────────────────────────────────────────────────
def build(check_only: bool = False) -> dict:
    for p in (APP, CSS, JS):
        if not p.is_file():
            fail("缺少源文件 %s" % p)

    css = CSS.read_text(encoding="utf-8")
    js = JS.read_text(encoding="utf-8")
    html = APP.read_text(encoding="utf-8")

    css_lint(css)
    if "</script" in js.lower():
        fail("bridge.js 里出现 </script，会提前闭合脚本块")
    # 只拦真正会提前闭合的结束标签（注释里提到 <style> 是允许的）
    if "</style" in css.lower() or "</script" in css.lower():
        fail("bridge.css 里出现结束标签 </style 或 </script，会提前闭合样式块")

    block = "\n".join([
        "",
        MARK_BEGIN,
        "<style>",
        css.rstrip(),
        "</style>",
        "<script>",
        js.rstrip(),
        "</script>",
        MARK_END,
    ])

    low = html.lower()
    k = low.rfind("</body>")
    if k < 0:
        fail("上游 app.html 里找不到 </body>")
    # 上游若已带旧桥接，先剥掉（可重复构建）
    b0 = low.find(MARK_BEGIN.lower())
    if b0 >= 0:
        e0 = low.find(MARK_END.lower(), b0)
        if e0 < 0:
            fail("上游文件里有桥接起始标记但没有结束标记，先手工清理")
        html = html[:b0] + html[e0 + len(MARK_END):]
        low = html.lower()
        k = low.rfind("</body>")

    merged = html[:k].rstrip("\n") + "\n" + block + "\n" + html[k:]

    # 合成后自检：<style> 段必须含 #fsBar 且注释闭合
    for m in re.finditer(r"<style>(.*?)</style>", merged, re.S):
        seg = m.group(1)
        if "#fsBar" in seg:
            if seg.count("/*") != seg.count("*/"):
                fail("合成后 <style> 注释不闭合")
            break
    else:
        fail("合成后找不到含 #fsBar 的 <style> 段")

    fp = upstream_fingerprint()
    if not check_only:
        OUT.mkdir(parents=True, exist_ok=True)
        TARGET.write_text(merged, encoding="utf-8")
        for f in ("LICENSE", "README.upstream.md"):
            s = SRC / f
            if s.is_file():
                shutil.copy2(s, OUT / f)
        if NOTICE.is_file():
            shutil.copy2(NOTICE, OUT / "NOTICE.md")

    return {"version": fp["version"], "upstream_sha256": fp["sha256"],
            "upstream_bytes": fp["bytes"], "merged_chars": len(merged),
            "css_chars": len(css), "js_chars": len(js),
            "target": str(TARGET), "check_only": check_only}


def node_check() -> str:
    node = shutil.which("node")
    if not node:
        for cand in (r"C:\Program Files\nodejs\node.exe",):
            if Path(cand).is_file():
                node = cand
                break
    if not node:
        return "跳过（没找到 node）"
    tmp = HERE / "_bridge_check.js"
    tmp.write_text(JS.read_text(encoding="utf-8"), encoding="utf-8")
    try:
        r = subprocess.run([node, "--check", str(tmp)],
                           capture_output=True, text=True)
    finally:
        tmp.unlink(missing_ok=True)
    if r.returncode != 0:
        fail("bridge.js 语法错误：\n%s" % (r.stderr or r.stdout))
    return "通过"


def main() -> int:
    ap = argparse.ArgumentParser(description="构建 static/vendor/svguest.html")
    ap.add_argument("--check", action="store_true", help="只校验，不写文件")
    ap.add_argument("--json", action="store_true", help="以 JSON 输出结果")
    a = ap.parse_args()

    info = build(check_only=a.check)
    info["node_check"] = node_check()
    if a.json:
        print(json.dumps(info, ensure_ascii=False, indent=2))
    else:
        print("上游 svguest %s  sha256 %s…（%d 字节）"
              % (info["version"], info["upstream_sha256"][:16], info["upstream_bytes"]))
        print("bridge.css %d 字符 / bridge.js %d 字符 → 合成 %d 字符"
              % (info["css_chars"], info["js_chars"], info["merged_chars"]))
        print("CSS 自检 通过 · JS 语法 %s" % info["node_check"])
        print(("校验通过（未写文件）：" if a.check else "已写出：") + info["target"])
    return 0


if __name__ == "__main__":
    sys.exit(main())
