# -*- coding: utf-8 -*-
"""面板库：扫描目录 → 面板清单 + 缩略图缓存。

⚠ Windows MAX_PATH(260) 静默截断
--------------------------------
``os.walk`` / ``os.scandir`` 在超长路径下**枚举会静默截断且零报错**。
⇒ 本模块所有遍历/读取一律走 ``_long()`` 加 ``\\\\?\\`` 前缀（先 ``abspath``），
返回给前端的路径再 ``_short()`` 去掉前缀。
"""
from __future__ import annotations

import hashlib
import json
import os
import time
from pathlib import Path
from typing import Any, Dict, List, Optional

from PIL import Image

IMAGE_EXT = {".png", ".jpg", ".jpeg", ".tif", ".tiff", ".bmp", ".webp"}
PDF_EXT = {".pdf"}
SVG_EXT = {".svg"}
ALL_EXT = IMAGE_EXT | PDF_EXT | SVG_EXT

SKIP_DIRS = {
    ".git", ".svn", "__pycache__", "node_modules", ".ipynb_checkpoints",
    ".workbuddy-ai", "cache", "thumbs", ".idea", ".vscode", "$RECYCLE.BIN",
    "System Volume Information",
}

DEFAULT_MAX_DEPTH = 6
DEFAULT_LIMIT = 600
META_CACHE_TTL = 3600.0


def _long(p: str) -> str:
    """加 ``\\\\?\\`` 前缀绕过 MAX_PATH（仅 Windows）。"""
    ap = os.path.abspath(p)
    if os.name == "nt" and not ap.startswith("\\\\?\\"):
        return "\\\\?\\" + ap
    return ap


def _short(p: str) -> str:
    return p[4:] if p.startswith("\\\\?\\") else p


# ══════════════════════════════════════════════════════════════════
# 元数据
# ══════════════════════════════════════════════════════════════════
def _meta_of(path: str, dpi: int = 300) -> Optional[Dict[str, Any]]:
    """读单个面板的尺寸元数据（不加载全部像素）。失败返回 None。"""
    ext = Path(path).suffix.lower()
    try:
        if ext in PDF_EXT:
            import fitz

            doc = fitz.open(_long(path))
            try:
                if doc.page_count < 1:
                    return None
                r = doc.load_page(0).rect
                w_pt, h_pt = float(r.width), float(r.height)
            finally:
                doc.close()
            if w_pt <= 0 or h_pt <= 0:
                return None
            w_in, h_in = w_pt / 72.0, h_pt / 72.0
            px_w, px_h = int(round(w_in * dpi)), int(round(h_in * dpi))
        elif ext in IMAGE_EXT:
            with Image.open(_long(path)) as im:
                px_w, px_h = im.width, im.height
                src_dpi = dpi
                mdpi = im.info.get("dpi")
                if mdpi and mdpi[0] and 20 < float(mdpi[0]) < 4000:
                    src_dpi = float(mdpi[0])
                w_in, h_in = px_w / src_dpi, px_h / src_dpi
        elif ext in SVG_EXT:
            # SVG：物理尺寸走 svg_edit.svg_size_pt（**不要**用 fitz 的
            # page.rect：它把无单位 width="920" 当成 920pt → px 系 SVG
            # 的 w_in 虚胖 4/3，且与 inspect_svg 口径互相矛盾）。
            import svg_edit

            w_pt, h_pt = svg_edit.svg_size_pt(_long(path))
            if not w_pt or not h_pt:
                return None
            w_in, h_in = w_pt / 72.0, h_pt / 72.0
            px_w, px_h = int(round(w_in * dpi)), int(round(h_in * dpi))
        else:
            return None
    except Exception:                                            # noqa: BLE001
        return None
    if px_w <= 0 or px_h <= 0:
        return None
    return {"px_w": px_w, "px_h": px_h, "aspect": px_h / float(px_w),
            "w_in": round(w_in, 4), "h_in": round(h_in, 4)}


class MetaCache:
    """按 (path, mtime, size) 缓存元数据，避免每次全量重扫。"""

    def __init__(self, cache_dir: str):
        self.dir = cache_dir
        os.makedirs(cache_dir, exist_ok=True)
        self.file = os.path.join(cache_dir, "panel_meta.json")
        self.data: Dict[str, Any] = {}
        self._dirty = False
        try:
            with open(self.file, "r", encoding="utf-8") as f:
                self.data = json.load(f)
        except Exception:                                        # noqa: BLE001
            self.data = {}

    def get(self, path: str, dpi: int = 300) -> Optional[Dict[str, Any]]:
        try:
            st = os.stat(_long(path))
        except OSError:
            return None
        key = f"{path}|{int(st.st_mtime)}|{st.st_size}|{dpi}"
        hit = self.data.get(key)
        if hit is not None:
            return hit if hit else None
        m = _meta_of(path, dpi=dpi)
        self.data[key] = m or {}
        self._dirty = True
        return m

    def save(self) -> None:
        if not self._dirty:
            return
        try:
            # 只留最近 4000 条，防无限膨胀
            if len(self.data) > 4000:
                self.data = dict(list(self.data.items())[-4000:])
            with open(self.file, "w", encoding="utf-8") as f:
                json.dump(self.data, f, ensure_ascii=False)
            self._dirty = False
        except Exception:                                        # noqa: BLE001
            pass


# ══════════════════════════════════════════════════════════════════
# 扫描
# ══════════════════════════════════════════════════════════════════
def scan(root: str, cache: MetaCache, *, query: str = "",
         limit: int = DEFAULT_LIMIT, max_depth: int = DEFAULT_MAX_DEPTH,
         dpi: int = 300, newest_first: bool = True) -> Dict[str, Any]:
    """扫描 root 下的面板图，返回 ``{root, total, truncated, panels: [...]}``。"""
    root_abs = os.path.abspath(os.path.expanduser(root))
    if not os.path.isdir(root_abs):
        return {"root": root_abs, "total": 0, "truncated": False,
                "panels": [], "error": "目录不存在"}

    found: List[Dict[str, Any]] = []
    root_long = _long(root_abs)
    base_depth = root_long.rstrip("\\/").count(os.sep)

    for dirpath, dirnames, filenames in os.walk(root_long):
        dirnames[:] = [d for d in dirnames if d not in SKIP_DIRS and not d.startswith(".")]
        if dirpath.rstrip("\\/").count(os.sep) - base_depth >= max_depth:
            dirnames[:] = []
        for fn in filenames:
            ext = os.path.splitext(fn)[1].lower()
            if ext not in ALL_EXT:
                continue
            full = os.path.join(dirpath, fn)
            try:
                st = os.stat(full)
            except OSError:
                continue
            found.append({"path": _short(full), "name": fn,
                          "mtime": st.st_mtime, "size": st.st_size})

    found.sort(key=lambda d: d["mtime"], reverse=newest_first)

    q = (query or "").strip().lower()
    if q:
        found = [d for d in found if q in d["name"].lower()
                 or q in os.path.dirname(d["path"]).lower()]

    total = len(found)
    truncated = total > limit
    found = found[:limit]

    panels: List[Dict[str, Any]] = []
    for d in found:
        m = cache.get(d["path"], dpi=dpi)
        if not m:
            continue
        rel = os.path.relpath(d["path"], root_abs)
        panels.append({
            "src": d["path"],
            "rel": rel.replace("\\", "/"),
            "name": d["name"],
            "ext": os.path.splitext(d["name"])[1].lower(),
            "mtime": d["mtime"],
            "size_kb": round(d["size"] / 1024.0, 1),
            **m,
        })
    cache.save()
    return {"root": root_abs, "total": total, "truncated": truncated,
            "panels": panels}


# ══════════════════════════════════════════════════════════════════
# 缩略图
# ══════════════════════════════════════════════════════════════════
def _thumb_key(src: str, box: int) -> str:
    try:
        st = os.stat(_long(src))
        sig = f"{src}|{int(st.st_mtime)}|{st.st_size}|{box}"
    except OSError:
        sig = f"{src}|{box}"
    return hashlib.sha1(sig.encode("utf-8", "replace")).hexdigest()[:20]


def thumbnail(src: str, cache_dir: str, box: int = 300) -> Optional[str]:
    """生成/取用缩略图，返回缓存文件绝对路径。失败返回 None。"""
    ext = Path(src).suffix.lower()
    if ext not in ALL_EXT:
        return None
    os.makedirs(cache_dir, exist_ok=True)
    out = os.path.join(cache_dir, _thumb_key(src, box) + ".png")
    if os.path.exists(out):
        return out
    try:
        if ext in SVG_EXT:
            import svg_edit

            # SVG 缩略图：栅格化到 ~160 dpi 再缩，避免小字号 text 变成糊块
            im = svg_edit.rasterize_svg(open(_long(src), "rb").read(),
                                        dpi=160)
            if im is None:
                return None
        elif ext in PDF_EXT:
            import fitz

            doc = fitz.open(_long(src))
            try:
                if doc.page_count < 1:
                    return None
                page = doc.load_page(0)
                zoom = 2.0
                pix = page.get_pixmap(matrix=fitz.Matrix(zoom, zoom), alpha=False)
                im = Image.frombytes("RGB", (pix.width, pix.height), pix.samples)
            finally:
                doc.close()
        else:
            im = Image.open(_long(src)).convert("RGB")
        im.thumbnail((box, box), Image.LANCZOS)
        im.save(out, "PNG", optimize=True)
        return out
    except Exception:                                            # noqa: BLE001
        return None


def render_preview(src: str, cache_dir: str, max_px: int = 1600,
                   adjust: Optional[Dict[str, Any]] = None) -> Optional[str]:
    """生成/取用「画布用预览图」：长边不超过 max_px 的 RGB PNG。

    与 ``thumbnail()`` 分开缓存：缩略图给面板库网格（小、快），
    预览图给画布（大、清晰）。PDF 会被栅格化。
    """
    ext = Path(src).suffix.lower()
    if ext not in ALL_EXT:
        return None
    pdir = os.path.join(cache_dir, "preview")
    os.makedirs(pdir, exist_ok=True)
    # 缓存键必须含 adjust 指纹，否则改调色后预览会用旧缓存（静默不同步）
    import hashlib

    afp = ""
    if adjust:
        try:
            afp = hashlib.md5(json.dumps(
                adjust, sort_keys=True, ensure_ascii=False).encode(
                "utf-8")).hexdigest()[:10]
        except Exception:                                    # noqa: BLE001
            afp = "raw"
    out = os.path.join(pdir, _thumb_key(src, max_px) + ("_a" + afp if afp else "") + ".png")
    if os.path.exists(out):
        return out
    try:
        if ext in SVG_EXT:
            # 顺序与 composer.load_panel 保持一致：先矢量编辑，再栅格化，再图像级调色
            import svg_edit

            r = svg_edit.render_svg(src, max_px=max_px,
                                    edit=(adjust or {}).get("svg") or None)
            if r is None:
                return None
            im = r["img"]
            img_adj = {k: v for k, v in (adjust or {}).items() if k != "svg"}
            if img_adj:
                import composer

                if not composer.is_noop(img_adj):
                    im = composer.apply_adjust(im, img_adj)
        elif ext in PDF_EXT:
            import fitz

            doc = fitz.open(_long(src))
            try:
                if doc.page_count < 1:
                    return None
                page = doc.load_page(0)
                r = page.rect
                zoom = max(1.0, float(max_px) / max(r.width, r.height))
                pix = page.get_pixmap(matrix=fitz.Matrix(zoom, zoom), alpha=False)
                im = Image.frombytes("RGB", (pix.width, pix.height), pix.samples)
            finally:
                doc.close()
        else:
            im = Image.open(_long(src)).convert("RGB")
        if max(im.size) > max_px:
            im.thumbnail((max_px, max_px), Image.LANCZOS)
        im.save(out, "PNG", optimize=True)
        return out
    except Exception:                                            # noqa: BLE001
        return None


def clear_thumbs(cache_dir: str) -> int:
    n = 0
    if not os.path.isdir(cache_dir):
        return 0
    for fn in os.listdir(cache_dir):
        if fn.endswith(".png"):
            try:
                os.remove(os.path.join(cache_dir, fn))
                n += 1
            except OSError:
                pass
    return n


def _now() -> float:
    return time.time()
