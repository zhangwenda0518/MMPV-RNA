# -*- coding: utf-8 -*-
"""Figure Studio · 渲染核心。

把一份 layout JSON（面板图 + 位置/尺寸/字母/题注）渲染成印刷级 PNG / PDF。

坐标约定
--------
* layout 里所有几何量单位是 **英寸**，原点在 **左上角**（x 向右、y 向下）——
  这是前端编辑器的直觉坐标系。
* 渲染时转成 matplotlib 的左下原点：``bottom = H - y - h``。
* 输出 **不使用** ``bbox_inches="tight"``：编辑器里用户指定的画布尺寸必须
  被严格保持（tight 会按内容撑大/裁小画布，让尺寸不可预测）。

零形变纪律
----------
面板块的 ``w_in`` / ``h_in`` 与源图宽高比不一致时，``imshow(aspect="auto")``
会**静默拉伸**。本模块渲染时会逐块比对，把偏差写进返回的 ``warnings``
（前端据此提示），并提供 ``natural_size()`` 给前端「锁定比例 / 1:1 尺寸」用。

调色板 PNG 陷阱
---------------
R 的 ``ggsave`` / ``grDevices::png``（Windows GDI+）常写**调色板 PNG
（``mode='P'``）**，``np.asarray`` 会返回 (H,W) 索引数组 → 被 ``imshow``
当标量场套 ``image.cmap``（白底变亮黄，零报错）。
⇒ 本模块一律先把 alpha **合成到白底**再转 RGB。

调色（adjust）
--------------
面板是**已渲染的栅格图**，所以调色只能是**图像级**（滤镜 / 颜色替换 / 底色处理），
不是语义级重绘。每块面板可带一个 ``adjust`` dict（见 ``ADJUST_DEFAULTS``），
渲染与预览走**同一条** ``apply_adjust()``，保证所见即所得。

⚠ 算子顺序是**契约**（``apply_adjust`` 里按固定顺序执行），改顺序会改变结果：

1. ``bg`` 把透明区合成到底色（在 load_panel 里做）
2. ``remap`` 颜色替换 ← **必须最先**：要在原图**精确色值**上匹配，
   否则被亮度/伽马改过之后就对不上了（调色板 PNG 尤其依赖这点）
3. ``invert`` / ``grayscale``
4. ``autocontrast``
5. ``threshold`` / ``posterize``
6. ``brightness`` / ``contrast`` / ``saturation`` / ``sharpness``
7. ``gamma``
8. ``hue`` 色相旋转
9. ``tint`` 单色着色
10. ``flatten`` 近白像素 → 指定色（去灰底/去白底）
11. ``trim`` 裁掉四周纯色边框
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from PIL import Image, ImageEnhance, ImageOps  # noqa: E402

DEFAULT_DPI = 300
IMAGE_EXT = {".png", ".jpg", ".jpeg", ".tif", ".tiff", ".bmp", ".webp"}
PDF_EXT = {".pdf"}
SVG_EXT = {".svg"}
ALL_PANEL_EXT = IMAGE_EXT | PDF_EXT | SVG_EXT
ALL_EXT = IMAGE_EXT | PDF_EXT

# 含 CJK 字体回退，保证中文题注不出现方块
_FONT_STACK = [
    "Arial", "Helvetica", "DejaVu Sans",
    "Microsoft YaHei", "SimHei", "Noto Sans CJK SC", "PingFang SC",
]


# ══════════════════════════════════════════════════════════════════
# 调色（图像级）
# ══════════════════════════════════════════════════════════════════
ADJUST_DEFAULTS: Dict[str, Any] = {
    "brightness": 1.0,      # 1 = 原样
    "contrast": 1.0,
    "saturation": 1.0,      # 0 = 灰度
    "sharpness": 1.0,
    "gamma": 1.0,           # >1 提亮暗部
    "hue": 0.0,             # 色相旋转（度，-180..180）
    "grayscale": False,
    "invert": False,
    "autocontrast": False,
    "threshold": 0,         # 0=关；1..255 二值化阈值
    "posterize": 0,         # 0=关；2..8 每通道位数
    "bg": None,             # 透明区合成底色（"#ffffff"）
    "flatten": None,        # 近白像素 → 该色（"#ffffff" 去灰底）
    "flatten_tol": 12,
    "tint": None,           # {"color": "#E69F00", "amount": 0.35}
    "remap": None,          # [{"from": "#d62728", "to": "#E69F00", "tol": 6}]
    "trim": False,          # 裁掉四周纯色边框
    "trim_tol": 8,
}

# 只认这些键，防止前端塞进奇怪字段
ADJUST_KEYS = tuple(ADJUST_DEFAULTS.keys())


def hex_to_rgb(s: str) -> Tuple[int, int, int]:
    s = (s or "").strip().lstrip("#")
    if len(s) == 3:
        s = "".join(c * 2 for c in s)
    if len(s) != 6:
        raise ValueError(f"非法颜色: {s!r}")
    return tuple(int(s[i:i + 2], 16) for i in (0, 2, 4))          # type: ignore[return-value]


def rgb_to_hex(rgb) -> str:
    return "#%02x%02x%02x" % (int(rgb[0]), int(rgb[1]), int(rgb[2]))


def normalize_adjust(adj: Optional[Dict[str, Any]]) -> Dict[str, Any]:
    """补齐默认值、丢掉未知键。返回的 dict 里**只有非默认项**才可能被用到。"""
    out = dict(ADJUST_DEFAULTS)
    if isinstance(adj, dict):
        for k in ADJUST_KEYS:
            if k in adj and adj[k] is not None:
                out[k] = adj[k]
        # 「svg」子树是**矢量级**编辑指令，不经 normalize，原样透传
        if adj.get("svg"):
            out["svg"] = adj["svg"]
    return out


def is_noop(adj: Optional[Dict[str, Any]]) -> bool:
    """判断 adjust 是否等于「什么都不做」（用于跳过整段处理）。"""
    if not adj:
        return True
    a = normalize_adjust(adj)
    for k, v in ADJUST_DEFAULTS.items():
        if a[k] != v:
            return False
    return True


def _hsv_hue_shift(img: Image.Image, degrees: float) -> Image.Image:
    """色相旋转（PIL HSV 的 H 通道是 0..255 循环）。"""
    if not degrees:
        return img
    hsv = img.convert("HSV")
    arr = np.asarray(hsv, dtype=np.int16)
    shift = int(round(float(degrees) / 360.0 * 255.0))
    arr[:, :, 0] = (arr[:, :, 0] + shift) % 256
    return Image.fromarray(arr.astype(np.uint8), "HSV").convert("RGB")


def _flatten_near(img: Image.Image, color: str, tol: float) -> Image.Image:
    """把「每通道都 ≥ 255-tol」的近白像素替换成 color（去灰底 / 统一白底）。"""
    arr = np.asarray(img, dtype=np.float32)
    target = np.array(hex_to_rgb(color), dtype=np.float32)
    mask = (arr >= (255.0 - float(tol))).all(axis=2)
    arr[mask] = target
    return Image.fromarray(arr.astype(np.uint8), "RGB")


def _trim(img: Image.Image, tol: float) -> Image.Image:
    """裁掉四周与「左上角像素同色」的边框。整图同色时不裁（避免 0 尺寸）。"""
    arr = np.asarray(img, dtype=np.int16)
    bg = arr[0, 0].astype(np.int16)
    diff = np.abs(arr - bg).max(axis=2) > float(tol)
    rows = np.where(diff.any(axis=1))[0]
    cols = np.where(diff.any(axis=0))[0]
    if rows.size == 0 or cols.size == 0:
        return img
    box = (int(cols[0]), int(rows[0]), int(cols[-1]) + 1, int(rows[-1]) + 1)
    if box == (0, 0, img.width, img.height):
        return img
    return img.crop(box)


def apply_adjust(img: Image.Image, adj: Optional[Dict[str, Any]]) -> Image.Image:
    """按固定顺序施加图像级调色。顺序见模块 docstring（**契约，勿改**）。"""
    if is_noop(adj):
        return img
    a = normalize_adjust(adj)

    # 2. 颜色替换 —— 必须在所有增强之前，否则精确色值已变
    remap = a.get("remap") or []
    if remap:
        arr = np.asarray(img, dtype=np.float32)
        for rule in remap:
            try:
                # from == to 是「占位未填」状态（前端取色后默认填成同色）→ 跳过
                if str(rule["from"]).lower() == str(rule["to"]).lower():
                    continue
                frm = np.array(hex_to_rgb(rule["from"]), dtype=np.float32)
                to = np.array(hex_to_rgb(rule["to"]), dtype=np.float32)
            except Exception:                                    # noqa: BLE001
                continue
            tol = float(rule.get("tol", 0) or 0)
            mask = np.abs(arr - frm).max(axis=2) <= tol
            arr[mask] = to
        img = Image.fromarray(arr.astype(np.uint8), "RGB")

    # 3. 反相 / 灰度
    if a["invert"]:
        img = ImageOps.invert(img.convert("RGB"))
    if a["grayscale"]:
        img = img.convert("L").convert("RGB")

    # 4. 自动对比度
    if a["autocontrast"]:
        img = ImageOps.autocontrast(img.convert("RGB"))

    # 5. 二值化 / 色调分离
    thr = int(a["threshold"] or 0)
    if thr > 0:
        luma = np.asarray(img.convert("L"), dtype=np.uint8)
        bw = np.where(luma >= thr, 255, 0).astype(np.uint8)
        img = Image.fromarray(bw, "L").convert("RGB")
    pos = int(a["posterize"] or 0)
    if 2 <= pos <= 8:
        img = ImageOps.posterize(img.convert("RGB"), pos)

    # 6. 亮度 / 对比度 / 饱和度 / 锐度
    if float(a["brightness"]) != 1.0:
        img = ImageEnhance.Brightness(img).enhance(max(0.0, float(a["brightness"])))
    if float(a["contrast"]) != 1.0:
        img = ImageEnhance.Contrast(img).enhance(max(0.0, float(a["contrast"])))
    if float(a["saturation"]) != 1.0:
        img = ImageEnhance.Color(img).enhance(max(0.0, float(a["saturation"])))
    if float(a["sharpness"]) != 1.0:
        img = ImageEnhance.Sharpness(img).enhance(max(0.0, float(a["sharpness"])))

    # 7. 伽马
    g = float(a["gamma"])
    if g != 1.0 and g > 0:
        lut = (np.power(np.arange(256, dtype=np.float64) / 255.0, 1.0 / g) * 255.0)
        img = img.point(lut.astype(np.uint8).tolist() * 3)

    # 8. 色相
    if float(a["hue"]):
        img = _hsv_hue_shift(img, float(a["hue"]))

    # 9. 单色着色（按 amount 线性混向目标色）
    tint = a.get("tint")
    if isinstance(tint, dict) and tint.get("color"):
        try:
            amount = min(1.0, max(0.0, float(tint.get("amount", 0.3))))
            col = np.array(hex_to_rgb(tint["color"]), dtype=np.float32)
            arr = np.asarray(img, dtype=np.float32)
            arr = arr * (1.0 - amount) + col * amount
            img = Image.fromarray(arr.astype(np.uint8), "RGB")
        except Exception:                                        # noqa: BLE001
            pass

    # 10. 近白 → 指定色
    if a.get("flatten"):
        try:
            img = _flatten_near(img, a["flatten"], a.get("flatten_tol", 12))
        except Exception:                                        # noqa: BLE001
            pass

    # 11. 裁边
    if a.get("trim"):
        img = _trim(img, a.get("trim_tol", 8))

    return img


def dominant_colors(src: str, n: int = 12, dpi: int = DEFAULT_DPI,
                    adjust: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
    """提取面板的主要颜色（按像素占比降序）。

    给前端「颜色替换」用：``remap`` 要在**原图色值**上精确匹配，
    所以这里默认**不施加 adjust**（``adjust=None``），返回的色值可直接填进 from。
    桶宽 3 bit/通道（8 级）先聚合、再取桶内真实均值，避免抗锯齿碎色淹没主色。
    """
    # ── SVG 走语义色统计：数 DOM 里的色值出现次数，而不是像素占比。
    #   像素占比会被大面积背景/白色淹没，DOM 计数才是「这个色是不是语义色」。
    if Path(src).suffix.lower() in SVG_EXT:
        import svg_edit

        info = svg_edit.inspect_svg(src)
        if info.get("ok"):
            total = max(1, sum(c["count"] for c in info["colors"]))
            cols = [c for c in info["colors"]][:max(1, int(n))]
            for c in cols:
                c["frac"] = round(float(c["count"]) / total, 5)
            return {"colors": cols, "n_buckets": len(info["colors"]),
                    "px": 0, "semantic": True, "is_svg": True,
                    "aspect": (info["size_pt"][1] / info["size_pt"][0])
                    if info.get("size_pt") else None,
                    "src": src,
                    "legend": info.get("legend", []),
                    "n_text": info.get("n_text", 0),
                    "fonts": info.get("fonts", {}),
                    "flavor": info.get("flavor", "generic")}
    info = load_panel(src, dpi=dpi, adjust=adjust)
    im = info["img"]
    if max(im.size) > 480:
        im = im.copy()
        im.thumbnail((480, 480), Image.NEAREST)
    arr = np.asarray(im.convert("RGB"), dtype=np.float64).reshape(-1, 3)
    q = (arr.astype(np.uint16) >> 3).astype(np.int64)
    key = (q[:, 0] << 10) | (q[:, 1] << 5) | q[:, 2]
    uniq, inv, cnt = np.unique(key, return_inverse=True, return_counts=True)
    sums = np.zeros((uniq.size, 3), dtype=np.float64)
    np.add.at(sums, inv, arr)
    means = sums / cnt[:, None]
    order = np.argsort(-cnt)[:max(1, int(n))]
    total = float(arr.shape[0])
    out = []
    for i in order:
        rgb = means[i].round().astype(int)
        out.append({"hex": rgb_to_hex(rgb), "rgb": [int(v) for v in rgb],
                    "count": int(cnt[i]), "frac": round(float(cnt[i]) / total, 5)})
    return {"colors": out, "n_buckets": int(uniq.size),
            "px": int(info["px_w"] * info["px_h"]),
            "aspect": info["aspect"], "src": src}


def adjust_signature(adjust: Optional[Dict[str, Any]]) -> str:
    """adjust 的稳定短哈希（用于缓存文件名）。"""
    payload = json.dumps(normalize_adjust(adjust), sort_keys=True,
                         ensure_ascii=False, default=str)
    return hashlib.sha1(payload.encode("utf-8")).hexdigest()[:12]


def adjusted_preview(base_png: str, cache_dir: str,
                     adjust: Optional[Dict[str, Any]]) -> Optional[str]:
    """给**已缓存的基础预览图**施加调色，返回派生缓存文件路径。

    基础预览图由 ``panels.render_preview()`` 产出、已按 ``src+mtime`` 命名；
    这里再按 adjust 哈希派生一个文件 —— 同一面板的不同调色互不覆盖，
    也避免为了预览去重新栅格化 PDF。

    ⚠ 预览是**降采样后**再调色的：``remap`` 的容差在抗锯齿边缘会漏掉一部分
    过渡像素，所以预览边缘看着不如成品干净。**成品在原始分辨率上调色**，
    以成品为准。
    """
    if not adjust or is_noop(adjust):
        return base_png
    try:
        sig = adjust_signature(adjust)
        st = os.stat(base_png)
        stem = Path(base_png).stem
        adir = Path(cache_dir) / "adjusted"
        out = adir / f"{stem}_{int(st.st_mtime)}_{sig}.png"
        if out.is_file():
            return str(out)
        adir.mkdir(parents=True, exist_ok=True)
        with Image.open(base_png) as im:
            img = apply_adjust(im.convert("RGB"), adjust)
        img.save(out, "PNG", optimize=True)
        _prune_dir(adir, keep=400)
        return str(out)
    except Exception:                                            # noqa: BLE001
        return base_png


def _prune_dir(d: Path, keep: int = 400) -> None:
    """缓存目录按修改时间保留最近 keep 个文件（防无限膨胀）。"""
    try:
        files = sorted(d.glob("*.png"), key=lambda p: p.stat().st_mtime,
                       reverse=True)
        for f in files[keep:]:
            try:
                f.unlink()
            except OSError:
                pass
    except Exception:                                            # noqa: BLE001
        pass


# ══════════════════════════════════════════════════════════════════
# 源图读取
# ══════════════════════════════════════════════════════════════════
def _rasterize_pdf(path: Path, dpi: int) -> Image.Image:
    """把 PDF 第 1 页按 dpi 渲染成 RGB 位图（用 PyMuPDF）。"""
    import fitz  # PyMuPDF

    doc = fitz.open(str(path))
    try:
        if doc.page_count < 1:
            raise ValueError("PDF 无页面")
        page = doc.load_page(0)
        zoom = dpi / 72.0
        pix = page.get_pixmap(matrix=fitz.Matrix(zoom, zoom), alpha=False)
        img = Image.frombytes("RGB", (pix.width, pix.height), pix.samples)
    finally:
        doc.close()
    return img


def load_panel(src: str, dpi: int = DEFAULT_DPI,
               assume_dpi: int = DEFAULT_DPI,
               adjust: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
    """读取面板源图并施加调色，返回
    ``{img, px_w, px_h, aspect, w_in, h_in, src_dpi, raw_mode}``。

    ``aspect = px_h / px_w``（高/宽），是唯一真正要紧的量：
    它决定面板块「锁定比例」时的形状。⚠ **裁边（trim）会改变它** ——
    所以 aspect 一律取**调色之后**的图像尺寸，前后端才不会各算一套。
    ``w_in / h_in`` 是该图按自身 DPI 换算的**物理尺寸**（信息用，供「1:1」按钮）。
    """
    p = Path(src)
    ext = p.suffix.lower()
    a = normalize_adjust(adjust)
    raw_mode = "?"
    if ext in SVG_EXT:
        # ── SVG：先做**矢量级**编辑（改 XML），再栅格化，最后才走图像级调色。
        #   顺序是契约：颜色/文字/字号/图例必须在 XML 上做，栅格后就只剩像素可改了。
        import svg_edit

        r = svg_edit.render_svg(str(p), dpi=dpi, edit=a.get("svg") or None)
        if r is None:
            raise ValueError("SVG 读取/栅格化失败（fitz）：%s" % p.name)
        src_dpi = float(dpi)
        w_in, h_in = r["w_in"], r["h_in"]
        raw_mode = "SVG"
        img = r["img"]
        if is_noop({k: v for k, v in a.items() if k != "svg"}):
            return {
                "img": img,
                "px_w": img.width,
                "px_h": img.height,
                "aspect": r["aspect"],
                "w_in": w_in,
                "h_in": h_in,
                "src_dpi": src_dpi,
                "raw_mode": raw_mode,
                "svg_notes": r["notes"],
            }
        img = apply_adjust(img, {k: v for k, v in a.items() if k != "svg"})
        # trim 等算子会改像素尺寸 → w_in/h_in 按同一比例缩放（保持 CSS 口径）
        _px0_w, _px0_h = r["px_w"], r["px_h"]
        _kx = img.width / float(_px0_w) if _px0_w else 1.0
        _ky = img.height / float(_px0_h) if _px0_h else 1.0
        return {
            "img": img,
            "px_w": img.width,
            "px_h": img.height,
            "aspect": img.height / float(img.width),
            "w_in": w_in * _kx,
            "h_in": h_in * _ky,
            "src_dpi": src_dpi,
            "raw_mode": raw_mode,
            "svg_notes": r["notes"],
        }
    if ext in PDF_EXT:
        img = _rasterize_pdf(p, dpi)
        w_in, h_in, src_dpi = img.width / dpi, img.height / dpi, dpi
        raw_mode = "PDF"
    elif ext in IMAGE_EXT:
        raw = Image.open(p)
        raw_mode = raw.mode
        src_dpi = assume_dpi
        mdpi = raw.info.get("dpi")
        if mdpi and mdpi[0] and 20 < float(mdpi[0]) < 4000:
            src_dpi = float(mdpi[0])
        # ⚠ 调色板 PNG 陷阱：直接 np.asarray(P 图) 会得到 (H,W) 索引数组。
        #   同时 alpha 不能「丢掉」，要**合成到底色**（默认白），否则透明区取到的是
        #   存储层的原始色值（常为黑），白底会变黑块。
        has_alpha = "A" in raw.getbands() or (
            raw.mode == "P" and "transparency" in raw.info)
        if has_alpha:
            base = Image.new("RGB", raw.size,
                             hex_to_rgb(a.get("bg") or "#ffffff"))
            rgba = raw.convert("RGBA")
            base.paste(rgba, mask=rgba.split()[-1])
            img = base
        else:
            img = raw.convert("RGB")
        w_in, h_in = img.width / src_dpi, img.height / src_dpi
    else:
        raise ValueError(f"不支持的面板格式: {ext}")

    img = apply_adjust(img, adjust)
    # 调色（尤其 trim）可能改变尺寸 → 物理尺寸与比例一律按**结果**重算
    w_in, h_in = img.width / src_dpi, img.height / src_dpi
    return {
        "img": img,
        "px_w": img.width,
        "px_h": img.height,
        "aspect": img.height / float(img.width),
        "w_in": w_in,
        "h_in": h_in,
        "src_dpi": src_dpi,
        "raw_mode": raw_mode,
    }


def probe(src: str, dpi: int = DEFAULT_DPI,
          adjust: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
    """轻量探测：只返回尺寸信息（不含像素），给前端定初始大小 / 调色后重算比例。"""
    info = load_panel(src, dpi=dpi, adjust=adjust)
    return {k: v for k, v in info.items() if k != "img"}


# ══════════════════════════════════════════════════════════════════
# 渲染
# ══════════════════════════════════════════════════════════════════
def _apply_rcparams(dpi: int, font_family: Optional[List[str]] = None) -> None:
    plt.rcParams.update({
        "pdf.fonttype": 42,          # 嵌入 TrueType → 矢量可编辑（项目硬规则）
        "ps.fonttype": 42,
        "font.family": "sans-serif",
        "font.sans-serif": font_family or _FONT_STACK,
        "axes.unicode_minus": False,
        "savefig.dpi": dpi,
        "figure.dpi": 150,
    })


def _text_w_in(fig, txt: str, fontsize: float) -> float:
    """用**真实渲染器**量一段文字在画布上的宽度（英寸）。

    ⚠ 不能按「字符数 × 字号」估：中英混排、字重、字距都会让估算失准，
    而 ``fig.text`` 溢出**不会报错**、只会压到隔壁面板上（本项目踩过）。
    """
    t = fig.text(0.0, 0.0, txt, fontsize=fontsize)
    try:
        bb = t.get_window_extent(fig.canvas.get_renderer())
        return bb.width / fig.dpi
    finally:
        t.remove()


def _wrap_text(fig, txt: str, fontsize: float, max_w_in: float) -> List[str]:
    """把文字折行到 ``max_w_in`` 以内。先按空格切词，单个词超宽再按字符切。"""
    if max_w_in <= 0:
        return [str(txt)]
    out: List[str] = []
    for para in str(txt).split("\n"):
        cur = ""
        for word in para.split(" "):
            cand = (cur + " " + word) if cur else word
            if _text_w_in(fig, cand, fontsize) <= max_w_in:
                cur = cand
                continue
            if cur:
                out.append(cur)
            if _text_w_in(fig, word, fontsize) > max_w_in:
                buf = ""
                for ch in word:
                    if buf and _text_w_in(fig, buf + ch, fontsize) > max_w_in:
                        out.append(buf)
                        buf = ch
                    else:
                        buf += ch
                cur = buf
            else:
                cur = word
        out.append(cur)
    return [ln for ln in out if ln] or [""]


def _overlay_svg_vectors(pdf_path: str, pending) -> int:
    """把 SVG 面板以**矢量**形式嵌进已生成的 PDF。

    做法：先用 matplotlib 出好画布（面板位置留白），再用 PyMuPDF 的
    ``page.show_pdf_page(rect, srcdoc)`` 把「SVG 转成的 PDF」作为 Form
    XObject 叠进对应矩形 —— 文字与路径都保持矢量，不是 300dpi 位图。

    坐标：layout 用**左上原点、英寸**；PDF 是**pt（1in = 72pt）**，
    fitz 的 Rect 同样是左上原点，所以直接 x*72 / y*72 换算即可。
    任一面板转换失败就跳过（该区域留白），不影响其他面板。
    """
    import fitz

    n = 0
    closed = False
    try:
        doc = fitz.open(pdf_path)
    except Exception:                                            # noqa: BLE001
        return 0
    try:
        if doc.page_count < 1:
            return 0
        page = doc.load_page(0)
        for (x, y, w, h, src, adjust) in pending:
            try:
                import svg_edit

                r = svg_edit.render_svg(src, dpi=150,
                                        edit=(adjust or {}).get("svg"))
                if r is None:
                    continue
                svg_pdf = fitz.open(stream=r["svg_bytes"],
                                    filetype="svg").convert_to_pdf()
                srcdoc = fitz.open("pdf", svg_pdf)
                rect = fitz.Rect(x * 72.0, y * 72.0,
                                 (x + w) * 72.0, (y + h) * 72.0)
                page.show_pdf_page(rect, srcdoc)
                srcdoc.close()
                n += 1
            except Exception:                                    # noqa: BLE001
                continue
        if n:
            # ⚠ PyMuPDF 不允许「另存到原路径且非 incremental」→ 先写临时文件再替换。
            #   ⚠ 还必须**先 close 再 replace**：fitz 仍持有原文件句柄时
            #     os.replace 会 PermissionError(WinError 5)。
            tmp = pdf_path + ".tmp"
            doc.save(tmp, deflate=True, garbage=3)
            closed = True
            doc.close()
            os.replace(tmp, pdf_path)
            return n
    finally:
        if not closed:
            doc.close()


def render_layout(layout: Dict[str, Any], out_dir: str, stem: str,
                  dpi: Optional[int] = None,
                  formats: Tuple[str, ...] = ("png", "pdf"),
                  panel_dpi: Optional[int] = None,
                  svg_vector: bool = True) -> Dict[str, Any]:
    """按 layout 渲染并返回 ``{png, pdf, size_px, warnings, canvas_in}``。

    layout 结构::

        {
          "canvas": {"width_in": 10.0, "height_in": 12.0, "dpi": 300},
          "style":  {"bg": "white", "frame": false, "frame_color": "#bbbbbb",
                     "frame_lw": 0.6, "letter_fs": 11, "letter_dx": 0.0,
                     "letter_dy": 1.02, "caption_fs": 7.0,
                     "caption_gap_in": 0.06, "caption_linespacing": 1.2},
          "title":  {"text": "", "subtitle": "", "fs": 13, "sub_fs": 9.5,
                     "x_in": 0.08, "y_in": 0.30, "color": "#000000"},
          "items":  [{"id": "p1", "src": "/abs/a.png", "x_in": 0.08, "y_in": 0.08,
                      "w_in": 4.5, "h_in": 2.0, "letter": "A",
                      "caption": "", "z": 0,
                      "adjust": {"brightness": 1.05, "grayscale": false,
                                 "remap": [{"from": "#d62728",
                                            "to": "#D55E00", "tol": 6}]}}]
        }
        }
    """
    canvas = layout.get("canvas") or {}
    W = float(canvas.get("width_in", 10.0))
    H = float(canvas.get("height_in", 12.0))
    dpi = int(dpi or canvas.get("dpi", DEFAULT_DPI))
    style = layout.get("style") or {}
    items = list(layout.get("items") or [])

    _apply_rcparams(dpi, style.get("font_family"))

    out_dir = os.path.abspath(os.path.expanduser(out_dir))
    os.makedirs(out_dir, exist_ok=True)

    fig = plt.figure(figsize=(W, H))
    fig.patch.set_facecolor(style.get("bg", "white"))

    warnings: List[str] = []
    pdpi = int(panel_dpi or dpi)

    # z 序：数值大的画在上层（list 顺序即绘制顺序）
    items = sorted(items, key=lambda it: float(it.get("z", 0)))

    drawn = 0
    vector_pending = []          # [(x, y, w, h, src, adjust)] 走矢量嵌入的面板
    raster_for_pdf = []          # PNG 存盘后要移除的 AxesImage（让位给矢量）
    want_pdf = "pdf" in formats
    for it in items:
        src = it.get("src")
        if not src:
            continue
        is_svg_src = Path(src).suffix.lower() in SVG_EXT
        try:
            info = load_panel(src, dpi=pdpi, adjust=it.get("adjust"))
            for note in (info.get("svg_notes") or []):
                warnings.append(f"SVG {it.get('letter') or it.get('id')}: {note}")
        except Exception as e:                                   # noqa: BLE001
            warnings.append(f"面板 {it.get('letter') or it.get('id')} 读取失败：{e}")
            continue

        x = float(it.get("x_in", 0.0))
        y = float(it.get("y_in", 0.0))
        w = float(it.get("w_in", 1.0))
        h = float(it.get("h_in", 1.0))
        if w <= 0 or h <= 0:
            warnings.append(f"面板 {it.get('letter') or it.get('id')} 尺寸非正，已跳过")
            continue

        # ── 零形变检查 ──
        nat = info["aspect"]
        got = h / w
        dev = abs(got - nat) / nat if nat else 0.0
        if dev > 0.01:
            warnings.append(
                "面板 %s 形变 %.1f%%（面板比 %.4f vs 源图比 %.4f）"
                % (it.get("letter") or it.get("id"), dev * 100, got, nat)
            )

        # ── 出界检查 ──
        # matplotlib 对越界内容**零报错**，只静默裁掉（本项目踩过多次）→ 主动告警。
        if x < -1e-9 or y < -1e-9 or x + w > W + 1e-9 or y + h > H + 1e-9:
            warnings.append(
                "面板 %s 超出画布：块 (%.2f, %.2f) %.2f x %.2f in，画布 %.2f x %.2f in"
                % (it.get("letter") or it.get("id"), x, y, w, h, W, H)
            )

        bottom = H - y - h
        rect = [x / W, bottom / H, w / W, h / H]
        ax = fig.add_axes(rect)
        # ⚠ 必须**先照常画栅格**：PNG 就靠它。
        #   SVG 面板的矢量嵌入只对 PDF 有意义 —— 所以 PNG 存盘**之后**再把
        #   这些 AxesImage remove 掉，然后存 PDF、叠加矢量。
        #   （反过来做会让 PNG 里 SVG 面板变空白，且零报错。）
        ax_im = ax.imshow(np.asarray(info["img"]), interpolation="antialiased",
                          aspect="auto")
        if want_pdf and svg_vector and is_svg_src:
            vector_pending.append((x, y, w, h, src,
                                   normalize_adjust(it.get("adjust"))))
            raster_for_pdf.append(ax_im)
        ax.set_xticks([])
        ax.set_yticks([])

        if style.get("frame"):
            for s in ax.spines.values():
                s.set_visible(True)
                s.set_edgecolor(style.get("frame_color", "#bbbbbb"))
                s.set_linewidth(float(style.get("frame_lw", 0.6)))
        else:
            for s in ax.spines.values():
                s.set_visible(False)

        # ── 面板字母（轴分数坐标；dx=0 贴轴左缘，dy>1 抬到标题上方）──
        # ⚠ 若单面板源图里**已烘焙**字母（本项目 `build_panels()` 出的 panel_X.png
        #   自带 "X" 标号），此处再加会**双字母** → 用 `style.show_letters=False` 关掉。
        letter = it.get("letter") or ""
        if letter and style.get("show_letters", True):
            ax.text(float(style.get("letter_dx", 0.0)),
                    float(style.get("letter_dy", 1.02)),
                    str(letter), transform=ax.transAxes,
                    fontsize=float(style.get("letter_fs", 11)),
                    fontweight="bold", ha="left", va="bottom",
                    color=style.get("letter_color") or "#000000")

        # ── 题注（画布坐标，挂在面板下方左对齐）──
        cap = (it.get("caption") or "").strip()
        if cap:
            gap = float(style.get("caption_gap_in", 0.06))
            cfs = float(style.get("caption_fs", 7.0))
            # ⚠ 默认折到**面板宽**：题注是 fig.text，压到隔壁面板上不会报错。
            #   要更宽可显式给 caption_max_w_in。
            max_w = float(it.get("caption_max_w_in") or w)
            lines = _wrap_text(fig, cap, cfs, max_w)
            txt = fig.text(x / W, (H - y - h - gap) / H, "\n".join(lines),
                           fontsize=cfs,
                           ha="left", va="top", color=style.get("caption_color", "#222222"),
                           linespacing=float(style.get("caption_linespacing", 1.2)))
            # 题注是 fig.text：tight_layout 不会搬它，越界即被裁且**零报错**
            try:
                bb = txt.get_window_extent(fig.canvas.get_renderer())
                if bb.y0 < 0:
                    warnings.append(
                        "面板 %s 的题注超出画布下缘 %.2f in，会被裁掉"
                        % (it.get("letter") or it.get("id"), -bb.y0 / fig.dpi)
                    )
                if bb.x1 / fig.dpi > W + 1e-6:
                    warnings.append(
                        "面板 %s 的题注右溢出画布 %.2f in"
                        % (it.get("letter") or it.get("id"), bb.x1 / fig.dpi - W)
                    )
                # 折行是按面板宽算的；若还是压到隔壁，说明面板太窄或字号太大
                if len(lines) > 1 and _text_w_in(fig, max(lines, key=len), cfs) > max_w * 1.02:
                    warnings.append(
                        "面板 %s 的题注折行后仍超宽，请调小题注字号"
                        % (it.get("letter") or it.get("id")))
            except Exception:                                    # noqa: BLE001
                pass

        drawn += 1

    # ── 全局标题 / 副标题 ──
    title = layout.get("title") or {}
    if (title.get("text") or "").strip():
        fig.text(float(title.get("x_in", 0.08)) / W,
                 (H - float(title.get("y_in", 0.30))) / H,
                 title["text"], fontsize=float(title.get("fs", 13)),
                 fontweight="bold", ha="left", va="center",
                 color=title.get("color") or "#000000")
    if (title.get("subtitle") or "").strip():
        fig.text(float(title.get("x_in", 0.08)) / W,
                 (H - float(title.get("y_in", 0.30)) - 0.26) / H,
                 title["subtitle"], fontsize=float(title.get("sub_fs", 9.5)),
                 color=title.get("sub_color") or "#444444", ha="left", va="center")

    if drawn == 0:
        plt.close(fig)
        return {"png": None, "pdf": None, "size_px": [0, 0],
                "warnings": warnings + ["没有任何面板被绘制"], "canvas_in": [W, H]}

    out: Dict[str, Any] = {"png": None, "pdf": None,
                           # matplotlib 的 Agg 画布尺寸是 ``int(W*dpi)``（截断），
                           # 非整数英寸会差 1 px —— 先给截断值，PNG 落盘后再用
                           # 实际文件尺寸覆盖（见下），保证前端显示的就是真值。
                           "size_px": [int(W * dpi), int(H * dpi)],
                           "warnings": warnings, "canvas_in": [W, H]}

    base = Path(out_dir) / stem
    if "png" in formats:
        p = base.with_suffix(".png")
        # ⚠ 不用 bbox_inches="tight"：画布尺寸必须严格 = W x H
        fig.savefig(p, dpi=dpi, facecolor=fig.get_facecolor())
        # matplotlib 的 PNG 恒为 RGBA；白底图剥掉 alpha 通道 —— 体积更小，
        # 也避免部分投稿系统拒收带透明通道的图。
        try:
            with Image.open(p) as _im:
                out["size_px"] = [int(_im.width), int(_im.height)]
                if _im.mode != "RGB":
                    _im.convert("RGB").save(p, "PNG", optimize=True)
        except Exception:                                        # noqa: BLE001
            pass
        out["png"] = str(p)
    if "pdf" in formats:
        # PNG 已经拿过栅格了 —— 现在把 SVG 面板的栅格撤掉，让位给矢量嵌入。
        # ⚠ 顺序不能反：先存 PNG 再 remove，否则 PNG 里这些面板是空白（零报错）。
        for _im in raster_for_pdf:
            try:
                _im.remove()
            except Exception:                                    # noqa: BLE001
                pass
        p = base.with_suffix(".pdf")
        fig.savefig(p, dpi=dpi, facecolor=fig.get_facecolor())
        if vector_pending:
            n = _overlay_svg_vectors(str(p), vector_pending)
            warnings.append(
                "PDF SVG vector embed: %d panel(s) as vector, rest raster"
                "（矢量嵌入 %d 个 SVG 面板，其余为栅格）" % (n, n))
        out["pdf"] = str(p)

    plt.close(fig)
    return out


# ══════════════════════════════════════════════════════════════════
# CLI
# ══════════════════════════════════════════════════════════════════
def main(argv: Optional[List[str]] = None) -> int:
    ap = argparse.ArgumentParser(description="Figure Studio 渲染核心（layout JSON → PNG/PDF）")
    ap.add_argument("layout", help="layout JSON 路径")
    ap.add_argument("-o", "--out", required=True, help="输出目录")
    ap.add_argument("--stem", default="figure", help="输出文件名（不含扩展名）")
    ap.add_argument("--dpi", type=int, default=None, help="覆盖 canvas.dpi")
    ap.add_argument("--formats", default="png,pdf", help="png,pdf")
    args = ap.parse_args(argv)

    with open(args.layout, "r", encoding="utf-8") as f:
        layout = json.load(f)
    res = render_layout(layout, args.out, args.stem, dpi=args.dpi,
                        formats=tuple(s.strip() for s in args.formats.split(",") if s.strip()))
    for w in res["warnings"]:
        print("[warn]", w, file=sys.stderr)
    print(json.dumps({k: v for k, v in res.items() if k != "warnings"},
                     ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
