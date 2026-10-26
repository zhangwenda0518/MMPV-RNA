# -*- coding: utf-8 -*-
"""SVG 矢量级编辑：颜色 / 文字 / 字号 / 图例 / 元素显隐 / 线宽。

与栅格面板的「图像级调色」（composer.apply_adjust）互补：
PNG/PDF 只能改像素；SVG 是 XML，因此能做**语义级**编辑 ——
颜色替换是**精确色值匹配**（不靠容差，不会漏抗锯齿边缘）、
图例标签文字可以直接改、字号整体缩放、图例整体隐藏/移动。

⚠ 栅格化走 PyMuPDF：``fitz.open(stream=..., filetype="svg")``。
  本机 cairosvg 不可用（缺 cairo DLL），svglib 的 renderPM 也无后端。
⚠ 解析用 lxml，并**显式关掉实体解析/网络**（防 XXE）。
"""
from __future__ import annotations

import copy
import re
from collections import Counter
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from lxml import etree

SVG_NS = "http://www.w3.org/2000/svg"
XLINK_NS = "http://www.w3.org/1999/xlink"
NSMAP = {"svg": SVG_NS, "xlink": XLINK_NS}
_Q = lambda tag: "{%s}%s" % (SVG_NS, tag)

# ── 颜色归一化 ────────────────────────────────────────────────
_NAMED = {
    "white": "#ffffff", "black": "#000000", "red": "#ff0000",
    "green": "#008000", "blue": "#0000ff", "none": "none",
    "gray": "#808080", "grey": "#808080", "silver": "#c0c0c0",
}
_COLOR_ATTRS = ("fill", "stroke", "stop-color", "color", "flood-color",
                "lighting-color")
_DRAW_PROPS = ("fill", "stroke", "stroke-width", "fill-opacity",
               "stroke-opacity", "opacity", "stop-color", "font",
               "font-size", "font-family", "font-weight", "font-style",
               "letter-spacing", "text-anchor", "display",
               "stroke-linecap", "stroke-linejoin", "stroke-dasharray",
               "stroke-dashoffset", "fill-rule", "clip-rule", "visibility")


def _expand_short(hexv: str) -> str:
    """#abc -> #aabbcc；#abcd/#abcdabcd 只取前 6 位。"""
    h = hexv[1:]
    if len(h) in (3, 4):
        return "#" + "".join(c * 2 for c in h[:3])
    return ("#" + h[:6]) if len(h) >= 6 else ("#" + h.ljust(6, "0"))


def norm_color(v: Optional[str]) -> str:
    """把任意颜色写法归一成 ``#rrggbb``（小写）；``none`` 原样返回。"""
    if not v:
        return ""
    s = v.strip().lower()
    if s in ("none", "transparent", "currentcolor", "inherit"):
        return s
    if s.startswith("#"):
        return _expand_short(s)
    m = re.match(r"rgba?\(([^)]*)\)", s)
    if m:
        parts = [p.strip() for p in re.split(r"[,\s/]+", m.group(1)) if p.strip()]
        if len(parts) >= 3:
            ch = []
            for p in parts[:3]:
                if p.endswith("%"):
                    ch.append(int(round(float(p[:-1]) * 255 / 100.0)))
                else:
                    ch.append(int(round(float(p))))
            return "#%02x%02x%02x" % tuple(max(0, min(255, c)) for c in ch)
    return _NAMED.get(s, s)


def is_bg(hexv: str) -> bool:
    """是否白/近白（用于主色过滤：背景不当作「可替换的语义色」）。"""
    if not hexv or hexv == "none" or not hexv.startswith("#"):
        return False
    try:
        r, g, b = (int(hexv[i:i + 2], 16) for i in (1, 3, 5))
    except ValueError:
        return False
    return r >= 245 and g >= 245 and b >= 245


def _lname(el) -> str:
    t = el.tag
    return t.split("}", 1)[1] if isinstance(t, str) and t.startswith("{") else str(t)


# ── style 属性读写 ───────────────────────────────────────────
def parse_style(s: Optional[str]) -> Dict[str, str]:
    out: Dict[str, str] = {}
    if not s:
        return out
    for part in s.split(";"):
        if ":" in part:
            k, v = part.split(":", 1)
            k = k.strip()
            if k:
                out[k] = v.strip()
    return out


def dump_style(d: Dict[str, str]) -> str:
    return "; ".join("%s: %s" % (k, v) for k, v in d.items())


# ── 尺寸/字号/字体的轻量解析 ─────────────────────────────────
_LEN_RE = re.compile(r"^\s*([\d.]+)\s*(px|pt|mm|cm|in|pc|%)?\s*$")
_UNIT_TO_PT = {"px": 0.75, "pt": 1.0, "mm": 72.0 / 25.4, "cm": 72.0 / 2.54,
               "in": 72.0, "pc": 12.0, "%": 1.0}


def len_to_pt(v: Optional[str], viewport: float = 0.0) -> Optional[float]:
    if not v:
        return None
    m = _LEN_RE.match(v)
    if not m:
        return None
    n = float(m.group(1))
    if m.group(2) == "%":
        return n / 100.0 * viewport if viewport else None
    return n * _UNIT_TO_PT.get(m.group(2) or "px", 0.75)


_FONT_SIZE_RE = re.compile(r"(?<![\w.#-])(\d+(?:\.\d+)?)(px|pt)\b")


_FONT_STYLE_KW = {"normal", "italic", "oblique", "small-caps"}
_FONT_WEIGHT_KW = {"bold", "bolder", "lighter"}
_FONT_LEN_RE = re.compile(r"^\d+(?:\.\d+)?(?:px|pt|em|rem|%)$")
_FONT_TOKEN_RE = re.compile(r'"[^"]*"|\'[^\']*\'|\S+')
_FONT_WEIGHT_NUM_RE = re.compile(r"^[1-9]00$")


def split_font_shorthand(v: str):
    """拆 CSS `font` 简写。

    `italic bold 10px 'DejaVu Sans', sans-serif`
    → (font-style, font-weight, font-size, font-family)

    ⚠ 必须拆：PyMuPDF 的 SVG 渲染器**不解析 `font` 简写里的字号**
    （改 `font: 7.5px` 后栅格 0 差异；拆成 `font-size` 后 28k 差异像素）。
    """
    size, family, style_v, weight = "", "", "", ""
    rest = []
    for tk in _FONT_TOKEN_RE.findall(v.strip()):
        low = tk.lower()
        if _FONT_LEN_RE.match(low) and not size:
            size = tk
        elif low in _FONT_STYLE_KW:
            style_v = tk
        elif low in _FONT_WEIGHT_KW or _FONT_WEIGHT_NUM_RE.match(low):
            weight = tk
        else:
            rest.append(tk)
    return style_v, weight, size, " ".join(rest)


def expand_font_shorthand(root) -> int:
    """把整棵树的 `font` 简写统一拆成 font-size/font-family/… 。

    在 `apply_svg_edit` 里无条件先跑一遍 —— 否则任何字号/字体编辑
    都会静默失效（改 XML 有笔记，栅格结果 0 差异）。
    """
    n = 0
    for el in root.iter():
        if not isinstance(el.tag, str):
            continue
        st = parse_style(el.get("style"))
        src = None
        if st.get("font"):
            src = st.pop("font")
        elif el.get("font"):
            src = el.get("font")
            el.attrib.pop("font", None)
        if not src:
            continue
        style_v, weight, size, family = split_font_shorthand(src)
        if size:
            st["font-size"] = size
        if family:
            st["font-family"] = family
        if style_v:
            st["font-style"] = style_v
        if weight:
            st["font-weight"] = weight
        el.set("style", dump_style(st))
        n += 1
    return n


def scale_font_value(v: str, factor: float) -> str:
    """缩放 CSS font 简写或 font-size 值里的字号数字（只动字号那一项）。"""
    if factor == 1.0:
        return v
    def repl(m):
        return "%.4g%s" % (float(m.group(1)) * factor, m.group(2))
    return _FONT_SIZE_RE.sub(repl, v)


def font_size_of(style: Dict[str, str]) -> Optional[float]:
    """从 style 里取字号（pt）。"""
    for key in ("font", "font-size"):
        v = style.get(key)
        if not v:
            continue
        m = _FONT_SIZE_RE.search(v)
        if m:
            return float(m.group(1)) * _UNIT_TO_PT[m.group(2)]
        n = len_to_pt(v)
        if n:
            return n
    for key in ("font-size",):
        if style.get(key):
            n = len_to_pt(style[key])
            if n:
                return n
    return None


_FONT_FAM_RE = re.compile(r"'([^']+)'|\"([^\"]+)\"|([A-Za-z][\w -]*)")


def font_family_of(style: Dict[str, str]) -> str:
    for key in ("font", "font-family"):
        v = style.get(key)
        if not v:
            continue
        if key == "font":                      # 简写：取字号之后的部分
            m = _FONT_SIZE_RE.search(v)
            v = v[m.end():] if m else v
            v = re.sub(r"^/\s*[\d.]+\s*", "", v).strip()
        elif v:
            pass
        if v:
            v = v.split(",")[0].strip()
            return v.strip("'\"") or v
    return ""


# ── 安全解析 ────────────────────────────────────────────────
def _parser() -> etree.XMLParser:
    return etree.XMLParser(resolve_entities=False, no_network=True,
                           load_dtd=False, dtd_validation=False,
                           huge_tree=False, recover=True)


def parse_svg(source) -> Optional[etree._Element]:
    """source 可以是路径 / bytes / str。失败返回 None。"""
    try:
        if isinstance(source, (bytes, bytearray)):
            return etree.fromstring(bytes(source), _parser())
        if isinstance(source, str):
            p = Path(source)
            if p.is_file():
                return etree.parse(_winlong(str(p)), _parser()).getroot()
            return etree.fromstring(source, _parser())
        return None
    except Exception:
        return None


_SEP = chr(92)                       # backslash, written via chr() to survive tooling
_LONGP = _SEP * 2 + "?" + _SEP      # Win32 extended-length path prefix

def _winlong(p):
    """加 Win32 长路径前缀（仅 Windows）：绕过 MAX_PATH(260) 静默截断。"""
    import os
    ap = os.path.abspath(p)
    if os.name == "nt" and not ap.startswith(_LONGP):
        return _LONGP + ap
    return ap


# 探测：把 SVG 摊平成可编辑的清单
def _root_size(root):
    """返回 (w_pt, h_pt, viewBox)。宽高缺省时回落到 viewBox。"""
    vb = None
    if root.get("viewBox"):
        try:
            vb = [float(x) for x in re.split(r"[,\s]+", root.get("viewBox").strip())]
        except ValueError:
            vb = None
    w = len_to_pt(root.get("width"), vb[2] if vb and len(vb) >= 3 else 0.0)
    h = len_to_pt(root.get("height"), vb[3] if vb and len(vb) >= 4 else 0.0)
    if (w is None or h is None) and vb and len(vb) >= 4:
        w = w if w is not None else vb[2]
        h = h if h is not None else vb[3]
    return w, h, vb


def _child_bbox(el):
    """粗估子元素 bbox（只读 x/y/cx/cy/d，够用于图例定位显示）。"""
    xs = []
    ys = []

    def num(v):
        try:
            return float(re.sub(r"[a-z%]+$", "", v.strip()))
        except (ValueError, AttributeError):
            return None

    for a in ("x", "cx", "x1"):
        n = num(el.get(a)) if el.get(a) else None
        if n is not None:
            xs.append(n)
    for a in ("y", "cy", "y1"):
        n = num(el.get(a)) if el.get(a) else None
        if n is not None:
            ys.append(n)
    d = el.get("d")
    if d:
        pts = re.findall(r"(-?\d+(?:\.\d+)?)[ ,]+(-?\d+(?:\.\d+)?)", d)
        if pts:
            xs.extend(float(p[0]) for p in pts[:400])
            ys.extend(float(p[1]) for p in pts[:400])
    if not xs or not ys:
        for c in el:
            bb = _child_bbox(c)
            if bb:
                xs += [bb[0], bb[2]]
                ys += [bb[1], bb[3]]
    if not xs or not ys:
        return None
    return [min(xs), min(ys), max(xs), max(ys)]


def _flavor(ids):
    pats = [re.sub(r"\d+$", "", i) for i in ids]
    if "figure_" in pats or "axes_" in pats or "patch_" in pats or "line2d_" in pats:
        return "matplotlib"
    return "generic"


def svg_size_pt(source):
    """SVG 的物理尺寸 (w_pt, h_pt)，单位口径统一为 **pt**。

    ⚠⚠ 这是「尺寸」的唯一真相源：
    - fitz 把 `width="920"`（无单位 = px）当成 **920 pt** → 物理尺寸虚胖 4/3
      （920pt=12.78in vs 正确 690pt=9.58in）。`panels._meta_of` 早年直接用
      fitz 的 page.rect → px 系 SVG 的 w_in 全错、且与 inspect_svg 不一致。
    - matplotlib SVG 用 `509.76pt`，R/svglib 用裸 px → 两种口径必须归一。
    返回 (None, None) 表示解析失败。
    """
    root = parse_svg(source)
    if root is None:
        return None, None
    w, h, _vb = _root_size(root)
    if not w or not h or w <= 0 or h <= 0:
        return None, None
    return float(w), float(h)


def inspect_svg(source):
    """把 SVG 摊平成：尺寸 / 颜色 / 文本 / 组结构 / 图例候选。"""
    root = parse_svg(source)
    if root is None:
        return {"ok": False, "error": "SVG not valid XML"}
    w_pt, h_pt, vb = _root_size(root)

    colors = Counter()
    colors_fill = Counter()
    colors_stroke = Counter()
    named_colors = Counter()
    fonts = Counter()
    sizes = Counter()

    for el in root.iter():
        if not isinstance(el.tag, str):
            continue
        for a in _COLOR_ATTRS:
            v = el.get(a)
            if v:
                named_colors[v.strip().lower()] += 1
        st = parse_style(el.get("style"))
        if not st:
            continue
        for k, v in st.items():
            if k in _COLOR_ATTRS:
                named_colors[v.strip().lower()] += 1
        fs = font_size_of(st)
        if fs:
            sizes[round(fs, 2)] += 1
        ff = font_family_of(st)
        if ff:
            fonts[ff] += 1

    for raw, cnt in named_colors.items():
        n = norm_color(raw)
        if not n.startswith("#"):
            continue
        colors[n] += cnt

    for el in root.iter():
        if not isinstance(el.tag, str):
            continue
        st = parse_style(el.get("style"))
        fill = st.get("fill") or el.get("fill") or ""
        stroke = st.get("stroke") or el.get("stroke") or ""
        for raw, slot in ((fill, colors_fill), (stroke, colors_stroke)):
            n = norm_color(raw)
            if n.startswith("#"):
                slot[n] += 1

    texts = []
    for el in root.iter(_Q("text")):
        st = parse_style(el.get("style"))
        txt = re.sub(r"\s+", " ", "".join(el.itertext())).strip()
        if not txt:
            continue
        try:
            fx = float(el.get("x") or 0)
            fy = float(el.get("y") or 0)
        except ValueError:
            fx, fy = 0.0, 0.0
        texts.append({
            "idx": len(texts),
            "text": txt,
            "size_pt": font_size_of(st),
            "fill": norm_color(st.get("fill") or el.get("fill") or ""),
            "family": font_family_of(st),
            "x": fx,
            "y": fy,
            "anchor": st.get("text-anchor") or el.get("text-anchor") or "start",
        })

    groups = []
    legend = []
    g_ids = []
    for gi, el in enumerate(root.iter(_Q("g"))):
        gid = el.get("id") or ""
        cls = el.get("class") or ""
        if gid:
            g_ids.append(gid)
        kids = [c for c in el if isinstance(c.tag, str)]
        n_text = sum(1 for _ in el.iter(_Q("text")))
        n_draw = sum(1 for t in ("path", "use", "rect", "circle", "line",
                                 "ellipse", "polygon", "polyline")
                     for _ in el.iter(_Q(t)))
        bb = _child_bbox(el)
        hay = (gid + " " + cls).lower()
        entry = {
            "idx": gi, "id": gid, "class": cls, "n_kids": len(kids),
            "n_text": n_text, "n_draw": n_draw,
            "bbox": [round(v, 2) for v in bb] if bb else None,
            "transform": el.get("transform") or "",
            "is_legend": "legend" in hay,
        }
        groups.append(entry)
        if entry["is_legend"]:
            legend.append(entry)

    if not legend:
        for g in groups:
            if g["n_text"] >= 2 and g["n_draw"] >= 2:
                legend.append(dict(g, is_legend=True, heuristic=True))

    return {
        "ok": True,
        "size_pt": [round(w_pt, 3), round(h_pt, 3)] if w_pt and h_pt else None,
        "size_in": [round(w_pt / 72.0, 4), round(h_pt / 72.0, 4)] if w_pt and h_pt else None,
        "viewBox": vb,
        "flavor": _flavor(g_ids),
        "colors": [
            {"hex": h, "count": c, "fill": colors_fill.get(h, 0),
             "stroke": colors_stroke.get(h, 0), "is_bg": is_bg(h)}
            for h, c in colors.most_common(40)
        ],
        "texts": texts,
        "n_text": len(texts),
        "groups": groups[:200],
        "n_groups": len(groups),
        "legend": legend,
        "fonts": {"families": [k for k, _ in fonts.most_common(8)],
                  "sizes": [k for k, _ in sizes.most_common(10)]},
    }


# 选择器：支持 #id / .class / tag / tag.class / #id tag
_TRANSLATE_RE = re.compile(r"translate\(\s*([-\d.]+)[ ,]+([-\d.]+)")


def _sel_matches(el, sel):
    sel = sel.strip()
    if not sel:
        return False
    if sel == "*":
        return True
    parts = sel.split()
    last = parts[-1]
    name = _lname(el)
    ok = True
    if last.startswith("#"):
        ok = (el.get("id") or "") == last[1:]
    elif last.startswith("."):
        ok = last[1:] in (el.get("class") or "").split()
    elif "." in last:
        tag, cls = last.split(".", 1)
        ok = name == tag and cls in (el.get("class") or "").split()
    else:
        ok = name == last
    if not ok:
        return False
    anc = el.getparent()
    for p in parts[:-1]:
        while anc is not None:
            if _sel_matches(anc, p):
                break
            anc = anc.getparent()
        if anc is None:
            return False
    return True


def _set_color(el, key, newval):
    """把某个颜色属性写回元素（同时覆盖 style 简写与直接属性两种写法）。"""
    st = parse_style(el.get("style"))
    if key in st or el.get("style"):
        if key in st:
            st[key] = newval
            el.set("style", dump_style(st))
            return True
    if el.get(key) is not None:
        el.set(key, newval)
        return True
    return False


def _color_hits(el, key, target):
    """元素上该颜色属性的现值是否命中 target。"""
    st = parse_style(el.get("style"))
    v = st.get(key) or el.get(key) or ""
    return norm_color(v) == target


def apply_svg_edit(source, edit):
    """在内存里改 SVG，返回 (svg_bytes, notes)。

    edit 键（全部可选，可组合）：
      colors:   [{from,to,scope:"all|fill|stroke"}]  精确色值替换（不靠容差）
      texts:    [{from,to}]                          文本内容替换（图例标签直接改）
      font_scale / font_family: 全局字号倍数 / 字体
      text_color: 统一所有文字颜色
      stroke_scale: 线宽倍数
      hide:     ["#legend_1", ".ticks"]               按选择器隐藏
      legend:   {selector, action:"hide|move|scale", dx, dy, font_scale}
    """
    root = parse_svg(source)
    if root is None:
        return None, ["SVG parse failed"]
    tree = root.getroottree()
    notes = []

    # ⚠⚠ 先拆 `font` 简写：PyMuPDF 渲染器不解析它里面的字号，导致
    #   font_scale / legend scale / font_family 全部静默失效（0 差异像素）。
    #   同时修基线保真度：未拆时所有文字按默认字号渲染（实测 49k 差异像素）。
    n_split = expand_font_shorthand(root)
    if n_split and any(k in edit for k in
                       ("font_scale", "font_family", "legend", "text_color",
                        "texts", "stroke_scale")):
        notes.append("font shorthand expanded on %d nodes" % n_split)

    for rule in (edit.get("colors") or []):
        frm = norm_color(rule.get("from") or "")
        to = norm_color(rule.get("to") or "")
        if not frm.startswith("#") or not to.startswith("#"):
            notes.append("color rule skipped: " + str(rule.get("from")))
            continue
        scope = (rule.get("scope") or "all").lower()
        keys = ("fill", "stroke") if scope == "all" else (scope,)
        n = 0
        for el in root.iter():
            if not isinstance(el.tag, str):
                continue
            for k in keys:
                if _color_hits(el, k, frm) and _set_color(el, k, to.upper()):
                    n += 1
        notes.append("color %s -> %s on %d slots" % (frm, to.upper(), n))

    if edit.get("font_scale"):
        f = float(edit["font_scale"])
        n = 0
        for el in root.iter():
            if not isinstance(el.tag, str):
                continue
            st = parse_style(el.get("style"))
            hit = False
            for k in ("font", "font-size"):
                if k in st:
                    st[k] = scale_font_value(st[k], f)
                    hit = True
            if hit:
                el.set("style", dump_style(st))
                n += 1
            if el.get("font-size"):
                el.set("font-size", scale_font_value(el.get("font-size"), f))
                n += 1
        notes.append("font scale x%.3g on %d nodes" % (f, n))

    if edit.get("font_family"):
        fam = str(edit["font_family"])
        n = 0
        for el in root.iter():
            if not isinstance(el.tag, str):
                continue
            st = parse_style(el.get("style"))
            hit = False
            if "font" in st:
                st["font"] = re.sub(
                    r"(\d+(?:\.\d+)?(?:px|pt))(.*)$",
                    lambda m: m.group(1) + " " + _q(fam), st["font"])
                hit = True
            if "font-family" in st:
                st["font-family"] = _q(fam)
                hit = True
            if hit:
                el.set("style", dump_style(st))
                n += 1
        notes.append("font family -> %s on %d nodes" % (fam, n))

    if edit.get("text_color"):
        c = norm_color(edit["text_color"]).upper()
        n = 0
        for el in root.iter(_Q("text")):
            st = parse_style(el.get("style"))
            if st.get("fill") or el.get("fill"):
                st["fill"] = c
                el.set("style", dump_style(st))
                n += 1
        notes.append("text color -> %s on %d text nodes" % (c, n))

    if edit.get("stroke_scale"):
        f = float(edit["stroke_scale"])
        n = 0
        for el in root.iter():
            if not isinstance(el.tag, str):
                continue
            st = parse_style(el.get("style"))
            if "stroke-width" in st:
                m = _LEN_RE.match(st["stroke-width"])
                if m:
                    st["stroke-width"] = "%.4g%s" % (
                        float(m.group(1)) * f, m.group(2) or "")
                    el.set("style", dump_style(st))
                    n += 1
        notes.append("stroke width x%.3g on %d nodes" % (f, n))

    for rule in (edit.get("texts") or []):
        frm = (rule.get("from") or "").strip()
        to = rule.get("to") or ""
        if not frm:
            continue
        n = 0
        for el in root.iter(_Q("text")):
            cur = re.sub(r"\s+", " ", "".join(el.itertext()))
            if cur.strip() == frm:
                for ch in list(el):
                    el.remove(ch)
                el.text = to
                n += 1
        notes.append("text [%s] -> [%s] on %d nodes" % (frm, to, n))

    # ⚠ 隐藏必须是「从树里移除节点」而不是 `display:none`
    #   （PyMuPDF 的 SVG 渲染器忽略 display:none：改 XML 后栅格 0 差异，
    #    实测移除节点 = 28803 差异像素）。浏览器里两种都生效，故无回归。
    for sel in (edit.get("hide") or []):
        n = 0
        for el in list(root.iter()):
            if el is root or not isinstance(el.tag, str):
                continue
            if _sel_matches(el, sel):
                parent = el.getparent()
                if parent is not None:
                    parent.remove(el)
                    n += 1
        notes.append("hide %r on %d nodes" % (sel, n))

    lg = edit.get("legend") or {}
    if lg:
        sel = lg.get("selector") or ""
        action = (lg.get("action") or "hide").lower()
        # ⚠ 可能命中**多个**图例（如 genome tracks 有 legend_1 + legend_2）。
        #   早年只取「节点最多的那一个」→ 另一个留着，看起来像没生效。
        cand = [el for el in root.iter(_Q("g"))
                if (sel and _sel_matches(el, sel)) or _legendish(el)]
        if not cand:
            notes.append("legend node not found (selector=%r)" % sel)
        for el in cand:
            gid = el.get("id") or "?"
            if action == "hide":
                # 同 hide：移除节点（display:none 被 fitz 忽略）。
                parent = el.getparent()
                if parent is not None:
                    parent.remove(el)
                notes.append("legend removed: " + gid)
            elif action == "move":
                dx = float(lg.get("dx", 0) or 0)
                dy = float(lg.get("dy", 0) or 0)
                cur = el.get("transform") or ""
                m = _TRANSLATE_RE.search(cur)
                if m:
                    nx = float(m.group(1)) + dx
                    ny = float(m.group(2)) + dy
                    el.set("transform",
                           cur[:m.start()] + "translate(%.4g, %.4g)" % (nx, ny)
                           + cur[m.end():])
                else:
                    el.set("transform",
                           (cur + " translate(%.4g, %.4g)" % (dx, dy)).strip())
                notes.append("legend %s moved by (%g, %g)" % (gid, dx, dy))
            elif action == "scale":
                f = float(lg.get("font_scale", 1.0) or 1.0)
                n = 0
                for sub in el.iter():
                    if not isinstance(sub.tag, str):
                        continue
                    st = parse_style(sub.get("style"))
                    hit = False
                    for k in ("font", "font-size"):
                        if k in st:
                            st[k] = scale_font_value(st[k], f)
                            hit = True
                    if hit:
                        sub.set("style", dump_style(st))
                        n += 1
                notes.append("legend %s text scaled x%.3g on %d nodes"
                             % (gid, f, n))

    try:
        out = etree.tostring(tree, xml_declaration=True, encoding="UTF-8")
    except Exception as e:                                     # noqa: BLE001
        return None, ["serialize failed: %s" % e]
    return out, notes


def _q(s):
    """给字体名加单引号（避免源码里直接写引号组合被工具链吞字符）。"""
    return "'" + str(s).replace("'", "") + "'"


def _legendish(el):
    hay = ((el.get("id") or "") + " " + (el.get("class") or "")).lower()
    return "legend" in hay


# ══════════════════════════════════════════════════════════════
# 栅格化（fitz）
# ══════════════════════════════════════════════════════════════
def rasterize_svg(svg_bytes, dpi=300):
    """把 SVG 字节流栅格化为 RGB PIL 图。失败返回 None。

    ⚠ cairosvg 本机不可用（缺 cairo DLL）、svglib 的 renderPM 无后端
      ⇒ 统一走 PyMuPDF：fitz.open(stream=..., filetype="svg")。
    ⚠ SVG 是 pt 坐标系：zoom = dpi / 72。
    ⚠⚠ 先无条件拆 `font` 简写：PyMuPDF 不解析它里面的字号，未拆时
      所有文字按**默认字号**渲染（实测 49k 差异像素）——这是基线
      保真度问题，不拆等于预览/缩略图都失真。
    """
    try:
        import fitz
    except ImportError:
        return None
    try:
        expanded, _n = apply_svg_edit(svg_bytes, {})
        if expanded:
            svg_bytes = expanded
        doc = fitz.open(stream=svg_bytes, filetype="svg")
        try:
            if doc.page_count < 1:
                return None
            page = doc.load_page(0)
            zoom = float(dpi) / 72.0
            pix = page.get_pixmap(matrix=fitz.Matrix(zoom, zoom), alpha=False)
            from PIL import Image
            return Image.frombytes("RGB", (pix.width, pix.height), pix.samples)
        finally:
            doc.close()
    except Exception:
        return None


def fit_dpi_for(svg_bytes, max_px: int, floor: float = 36.0) -> float:
    """按目标像素长边反解栅格化 dpi（避免 `dpi=max_px` 直接爆内存）。

    ⚠ 早年生 `render_preview` 传 `dpi=max_px`：一张 12.8 in 宽的 SVG 在
      1600 dpi 下要 20k×18k px ≈ 1 GB → MemoryError → `/api/img` 404。
      口径与 PDF 分支一致：zoom = max_px / max(w_pt, h_pt)，dpi = 72·zoom。
    """
    w_pt, h_pt = svg_size_pt(svg_bytes)
    if not w_pt or not h_pt:
        return 150.0
    m = max(w_pt, h_pt)
    return max(floor, min(float(max_px), 72.0 * float(max_px) / m))


def render_svg(source, dpi=300, edit=None, max_px=None):
    """便入口：读 SVG → 施加矢量编辑 → 栅格化。

    返回 {"img", "px_w", "px_h", "aspect", "svg_bytes", "notes", "w_in", "h_in"}；
    任一步失败返回 None。``max_px`` 给定时按长边反解 dpi（见 fit_dpi_for）。
    """
    if edit:
        data, notes = apply_svg_edit(source, edit)
        if data is None:
            return None
    else:
        data = open(_winlong(str(source)), "rb").read() \
            if not isinstance(source, (bytes, bytearray)) else bytes(source)
        notes = []
    if max_px:
        dpi = fit_dpi_for(data, int(max_px))
    img = rasterize_svg(data, dpi=dpi)
    if img is None:
        return None
    # ⚠ w_in/h_in 报 **CSS 正确物理尺寸**（svg_size_pt），不是
    #   img.width/dpi：fitz 把无单位 px 当成 pt → 裸 px 的 SVG（R/svglib）
    #   会虚胖 4/3。两侧口径必须一致，否则「1:1」与面板库信息对不上。
    cw, ch = svg_size_pt(data)
    w_in = (cw / 72.0) if cw else img.width / float(dpi)
    h_in = (ch / 72.0) if ch else img.height / float(dpi)
    return {
        "img": img,
        "px_w": img.width,
        "px_h": img.height,
        "aspect": img.height / float(img.width),
        "w_in": w_in,
        "h_in": h_in,
        "src_dpi": float(dpi),
        "svg_bytes": data,
        "notes": notes,
        "raw_mode": "SVG",
    }
