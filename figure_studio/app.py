# -*- coding: utf-8 -*-
"""Figure Studio · 本地可视化组图平台。

启动::

    C:/Python312/python.exe figure_studio/app.py --port 8770

然后浏览器打开 http://127.0.0.1:8770

能力
----
* 面板库：扫描任意目录下的 PNG/JPG/PDF，缩略图网格 + 搜索
* 画布：拖拽 / 缩放 / 吸附 / 层级 / 面板字母 / 题注
* 调色：图像级（亮度/对比度/饱和度/锐度/伽马/色相/灰度/反相/二值化/色调分离/
  底色填充/近白替换/单色着色/颜色替换/裁边），预览与成品走同一条 ``apply_adjust``
* 导出：后端按 layout 精确渲染印刷级 PNG + PDF（画布尺寸严格保持）
* 布局：JSON 存取，随时改、随时重出（这是「出图后还能改」的落地方式）
"""
from __future__ import annotations

import json
import os
import re
import shutil
import sys
import time
import traceback
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional

from flask import (Flask, Response, abort, jsonify, request, send_file,
                   send_from_directory, stream_with_context)

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import composer  # noqa: E402
import panels as panel_lib  # noqa: E402

CACHE_DIR = HERE / "cache"
LAYOUT_DIR = HERE / "layouts"
OUTPUT_DIR = HERE / "outputs"
CONFIG_FILE = HERE / "config.json"

for d in (CACHE_DIR, LAYOUT_DIR, OUTPUT_DIR):
    d.mkdir(parents=True, exist_ok=True)

app = Flask(__name__, static_folder=str(HERE / "static"), static_url_path="/static")
app.config["JSON_AS_ASCII"] = False

META_CACHE = panel_lib.MetaCache(str(CACHE_DIR))

DEFAULT_ROOTS = [
    r"D:\桌面\延伸基因组\MMPV-RNA\reports",
    r"D:\桌面\延伸基因组\MMPV-RNA\virome_phylo_pipeline",
    r"D:\桌面\MMPV主图档案库_20261014",
]


# ══════════════════════════════════════════════════════════════════
# 配置
# ══════════════════════════════════════════════════════════════════
def load_config() -> Dict[str, Any]:
    cfg = {"roots": list(DEFAULT_ROOTS), "last_root": DEFAULT_ROOTS[0]}
    try:
        with open(CONFIG_FILE, "r", encoding="utf-8") as f:
            cfg.update(json.load(f))
    except Exception:                                            # noqa: BLE001
        pass
    cfg["roots"] = [r for r in cfg.get("roots", []) if isinstance(r, str)]
    return cfg


def save_config(cfg: Dict[str, Any]) -> None:
    try:
        with open(CONFIG_FILE, "w", encoding="utf-8") as f:
            json.dump(cfg, f, ensure_ascii=False, indent=2)
    except Exception:                                            # noqa: BLE001
        pass


def _exists(p: str) -> bool:
    try:
        return os.path.isfile(panel_lib._long(p))
    except Exception:                                            # noqa: BLE001
        return False


def _safe_name(name: str) -> str:
    name = re.sub(r"[^0-9A-Za-z_\-\u4e00-\u9fff]+", "_", name or "").strip("_")
    return name[:80] or "untitled"


def _parse_adjust(raw: Any) -> Optional[Dict[str, Any]]:
    """adjust 可能来自 URL 查询串（JSON 文本）或请求体（已是 dict）。"""
    if not raw:
        return None
    if isinstance(raw, dict):
        return raw
    if isinstance(raw, str):
        try:
            d = json.loads(raw)
            return d if isinstance(d, dict) else None
        except Exception:                                        # noqa: BLE001
            return None
    return None


# ══════════════════════════════════════════════════════════════════
# 页面
# ══════════════════════════════════════════════════════════════════
@app.route("/")
def index():
    return send_from_directory(str(HERE / "static"), "index.html")


@app.route("/api/health")
def health():
    return jsonify({"ok": True, "time": time.time(),
                    "python": sys.version.split()[0],
                    "cwd": str(HERE)})


# ══════════════════════════════════════════════════════════════════
# 面板库
# ══════════════════════════════════════════════════════════════════
@app.route("/api/config", methods=["GET", "POST"])
def api_config():
    cfg = load_config()
    if request.method == "POST":
        body = request.get_json(force=True, silent=True) or {}
        if "roots" in body and isinstance(body["roots"], list):
            cfg["roots"] = [str(r) for r in body["roots"] if str(r).strip()]
        if "last_root" in body:
            cfg["last_root"] = str(body["last_root"])
        save_config(cfg)
    return jsonify(cfg)


@app.route("/api/panels")
def api_panels():
    root = request.args.get("root") or load_config().get("last_root") or ""
    q = request.args.get("q", "")
    limit = int(request.args.get("limit", panel_lib.DEFAULT_LIMIT))
    depth = int(request.args.get("depth", panel_lib.DEFAULT_MAX_DEPTH))
    try:
        res = panel_lib.scan(root, META_CACHE, query=q, limit=limit,
                             max_depth=depth)
    except Exception as e:                                       # noqa: BLE001
        traceback.print_exc()
        return jsonify({"root": root, "total": 0, "panels": [],
                        "error": f"{type(e).__name__}: {e}"}), 500
    cfg = load_config()
    if os.path.isdir(root):
        cfg["last_root"] = root
        if root not in cfg["roots"]:
            cfg["roots"].insert(0, root)
        save_config(cfg)
    return jsonify(res)


@app.route("/api/thumb")
def api_thumb():
    src = request.args.get("src", "")
    box = int(request.args.get("box", 300))
    if not src or not _exists(src):
        abort(404)
    p = panel_lib.thumbnail(src, str(CACHE_DIR), box=box)
    if not p:
        abort(404)
    return send_file(p, mimetype="image/png")


@app.route("/api/img")
def api_img():
    """画布用的预览图（长边 ≤ max）。可带 ``adjust``（URL 编码 JSON）实时调色。"""
    src = request.args.get("src", "")
    max_px = int(request.args.get("max", 1600))
    if not src or not _exists(src):
        abort(404)
    adj = _parse_adjust(request.args.get("adjust")) or {}
    # SVG 的矢量编辑必须在这里做（栅格化前），所以 adjust 要传进 render_preview
    p = panel_lib.render_preview(src, str(CACHE_DIR), max_px=max_px, adjust=adj)
    if not p:
        abort(404)
    img_adj = {k: v for k, v in adj.items() if k != "svg"}
    if img_adj and not composer.is_noop(img_adj):
        p = composer.adjusted_preview(p, str(CACHE_DIR), img_adj) or p
    return send_file(p, mimetype="image/png")


@app.route("/api/file")
def api_file():
    """原始文件下载（交付/查看用）。"""
    src = request.args.get("src", "")
    if not src or not _exists(src):
        abort(404)
    return send_file(panel_lib._long(src), as_attachment=True,
                     download_name=os.path.basename(src))


@app.route("/api/svg/raw")
def api_svg_raw():
    """原样输出 SVG 字节（inline）—— 给 svguest 精修窗口 fetch 用。

    与 ``/api/file`` 的区别：那个是 ``as_attachment`` 下载；这个要能被
    ``fetch()`` 直接读成 Blob 再喂给 svguest 的 ``readFile()``。

    ⚠⚠ Content-Type 必须写 ``image/svg`` 而**不是** ``image/svg+xml``：
       实测本机 Chrome 会把 ``image/svg+xml``（及 ``application/octet-stream``）
       的 **GET** 响应换成合成 204（无 Content-Type、无 body、零报错），
       fetch/XHR/<img> 全部中招 → 桥接那边表现为「载入成功但画布是空的」。
       （POST 的 ``/api/svg/source`` 用 image/svg+xml 反而正常，所以是
       「GET + 该 MIME」的组合问题；``image/svg`` 与本机 ``/api/file`` 一致，可用。）
    """
    src = request.args.get("src", "")
    if not src or not _exists(src):
        abort(404)
    if not src.lower().endswith(".svg"):
        abort(400)
    return send_file(panel_lib._long(src), mimetype="image/svg")


@app.route("/api/svg/inspect", methods=["POST"])
def api_svg_inspect():
    """SVG 结构探测：尺寸 / 语义颜色 / 文本清单 / 组结构 / 图例候选。

    这是「全方位调整」的入口：前端靠它生成颜色替换表与文本替换表。
    """
    body = request.get_json(force=True, silent=True) or {}
    src = body.get("src", "")
    if not src or not _exists(src):
        return jsonify({"error": "file missing"}), 404
    if not src.lower().endswith(".svg"):
        return jsonify({"error": "not an SVG"}), 400
    import svg_edit

    return jsonify(svg_edit.inspect_svg(src))


@app.route("/api/svg/source", methods=["POST"])
def api_svg_source():
    """返回施加矢量编辑后的 SVG 文本（UTF-8），前端可直接塞进 <img> 预览。

    body: {src, edit}。⚠ 只在内存里改，不写回磁盘。
    """
    body = request.get_json(force=True, silent=True) or {}
    src = body.get("src", "")
    edit = body.get("edit") or {}
    if not src or not _exists(src):
        return jsonify({"error": "file missing"}), 404
    if not src.lower().endswith(".svg"):
        return jsonify({"error": "not an SVG"}), 400
    import svg_edit

    data, notes = svg_edit.apply_svg_edit(src, edit)
    if data is None:
        return jsonify({"error": "svg edit failed", "notes": notes}), 500
    resp = app.response_class(data, mimetype="image/svg+xml; charset=utf-8")
    resp.headers["X-SVG-Notes"] = "; ".join(notes)[:2000]
    return resp


@app.route("/api/svg/writeback", methods=["POST"])
def api_svg_writeback():
    """svguest 精修后把 SVG **写回原路径**（先备份）。

    body: ``{src, svg}``
    ⚠ 备份放在 ``outputs/_svguest_backup/``，**不放在面板目录** ——
      否则 `.svg` 备份会被面板库扫描出来、当成一块面板。
    ⚠ 预览缓存键含 mtime/size（``panels._thumb_key``），写回后自动失效，
      无需手动清缓存。
    """
    body = request.get_json(force=True, silent=True) or {}
    res = _write_svg_back(body.get("src", ""), body.get("svg"))
    if "error" in res:
        code = res.pop("_code", 400)
        return jsonify(res), code
    return jsonify(res)


# ══════════════════════════════════════════════════════════════════
# svguest 深度整合：批量写回 / 导出落盘 / 备份历史 / 回滚 / 会话
# ══════════════════════════════════════════════════════════════════
BACKUP_DIR = OUTPUT_DIR / "_svguest_backup"
EXPORT_DIR = OUTPUT_DIR / "_svguest_export"


def _write_svg_back(src: str, svg: Any) -> Dict[str, Any]:
    """把 svg 文本写回 src（先备份）。返回 dict；出错时含 ``error`` 与 ``_code``。"""
    if not src or not _exists(src):
        return {"error": "file missing", "_code": 404}
    if not src.lower().endswith(".svg"):
        return {"error": "not an SVG", "_code": 400}
    if not isinstance(svg, str) or "<svg" not in svg:
        return {"error": "empty svg", "_code": 400}
    lp = panel_lib._long(src)
    try:
        bdir = str(BACKUP_DIR)
        os.makedirs(bdir, exist_ok=True)
        stamp = datetime.now().strftime("%Y%m%d_%H%M%S_%f")[:-3]
        base = os.path.basename(src)
        bak = os.path.join(bdir, "%s.bak_%s.svg" % (base, stamp))
        if os.path.isfile(lp):
            shutil.copy2(lp, bak)
        else:
            bak = None
        with open(lp, "wb") as f:
            f.write(svg.encode("utf-8"))
    except OSError as e:
        return {"error": "write failed: %s" % e, "_code": 500}
    return {"ok": True, "src": src, "backup": bak, "bytes": len(svg)}


@app.route("/api/svg/writeback_batch", methods=["POST"])
def api_svg_writeback_batch():
    """多标签批量写回。body: ``{items:[{src, svg}]}``

    逐项独立成败——一个失败不影响其余，前端按 ``results`` 汇总。
    """
    body = request.get_json(force=True, silent=True) or {}
    items = body.get("items")
    if not isinstance(items, list) or not items:
        return jsonify({"error": "items 为空"}), 400
    if len(items) > 40:
        return jsonify({"error": "一次最多 40 项"}), 400
    results = []
    for it in items:
        if not isinstance(it, dict):
            results.append({"ok": False, "error": "bad item"})
            continue
        r = _write_svg_back(it.get("src", ""), it.get("svg"))
        r.pop("_code", None)
        results.append(r)
    ok_n = sum(1 for r in results if r.get("ok"))
    return jsonify({"ok": ok_n == len(results), "saved": ok_n,
                    "total": len(results), "results": results})


def _human(n: int) -> str:
    for unit in ("B", "KB", "MB", "GB"):
        if n < 1024 or unit == "GB":
            return ("%d %s" % (n, unit)) if unit == "B" else ("%.1f %s" % (n, unit))
        n /= 1024.0
    return "%d B" % n


def _svg_to_pdf(svg_text: str) -> bytes:
    """SVG → **矢量** PDF（文字可提取，不是栅格图）。

    为什么要走服务端而不是浏览器：svguest 原生的 ``pdf()`` 是
    「iframe + ``window.print()``」，会弹打印对话框、依赖纸张设置、
    还得手动选「另存为 PDF」。服务端转 PDF 不弹窗、尺寸可控、能批量。

    ⚠⚠ 两个「零报错但结果错」的坑（都已在 svg_edit 里踩过）：
      ① **`font` 简写必须先拆** —— PyMuPDF 的 SVG 渲染器不解析简写里的
         字号，不拆则全部文字按默认字号排版（改 XML 有变化、栅格 0 差异）。
      ② **fitz 把无单位 px 当成 pt** → px 系 SVG 页面虚胖 4/3。尺寸真相源
         是 ``svg_edit.svg_size_pt()``，偏差 >0.5% 时重排到正确尺寸。
    """
    import fitz
    from lxml import etree

    import svg_edit

    raw = svg_text.encode("utf-8")
    # ⚠⚠ 必须传 **bytes** 而不是 str：带 `<?xml … encoding="UTF-8"?>` 声明的
    #   unicode 字符串会被 lxml 直接拒绝（"Unicode strings with encoding
    #   declaration are not supported"）→ parse_svg 返回 None → 尺寸校正
    #   静默失效。而浏览器 svgText() 输出的 SVG **一定**带这行声明。
    root = svg_edit.parse_svg(raw)
    if root is not None:
        try:
            if svg_edit.expand_font_shorthand(root):
                raw = etree.tostring(root, xml_declaration=True, encoding="UTF-8")
        except Exception:                                        # noqa: BLE001
            pass                     # 拆不了就让 fitz 原样试，别让导出整体失败

    doc = fitz.open(stream=raw, filetype="svg")
    try:
        pdf_bytes = doc.convert_to_pdf()
    finally:
        doc.close()

    w_pt, h_pt = svg_edit.svg_size_pt(raw)
    out = fitz.open(stream=pdf_bytes, filetype="pdf")
    try:
        if w_pt and h_pt and out.page_count:
            cur = out[0].rect
            if (abs(cur.width - w_pt) > max(1.0, w_pt * 0.005)
                    or abs(cur.height - h_pt) > max(1.0, h_pt * 0.005)):
                fixed = fitz.open()
                pg = fixed.new_page(width=w_pt, height=h_pt)
                pg.show_pdf_page(pg.rect, out, 0)
                out.close()
                out = fixed
        return out.tobytes()
    finally:
        try:
            out.close()
        except Exception:                                        # noqa: BLE001
            pass


@app.route("/api/svg/export", methods=["POST"])
def api_svg_export():
    """svguest 导出的图像**直接落到原图同目录**，不经过浏览器下载。

    body:
      - 位图/SVG：``{src, kind:'png'|'jpg'|'svg', dpi:int, data:<base64>}``
      - 矢量 PDF：``{src, kind:'pdf', svg:<文本>}``
        （PDF 由服务端用 PyMuPDF 生成，**不收 base64** —— 前端根本不该
         为了一份 PDF 去调打印对话框）

    ⚠ 落盘位置取 ``src`` 所在目录（没有 src 时落 ``outputs/_svguest_export/``）——
      这样面板库重新扫描同一根目录就能直接把导出的 PNG 当面板用，
      省掉「下载到 Downloads 再手动搬过去」。
    """
    body = request.get_json(force=True, silent=True) or {}
    kind = str(body.get("kind") or "png").lower()
    if kind not in ("png", "jpg", "jpeg", "svg", "pdf"):
        return jsonify({"error": "kind 只支持 png/jpg/svg/pdf"}), 400

    if kind == "pdf":
        svg_in = body.get("svg")
        if not isinstance(svg_in, str) or "<svg" not in svg_in:
            return jsonify({"error": "svg 为空"}), 400
        if len(svg_in) > 60 * 1024 * 1024:
            return jsonify({"error": "SVG 过大"}), 413
        try:
            raw = _svg_to_pdf(svg_in)
        except Exception as e:                                   # noqa: BLE001
            return jsonify({"error": "PDF 生成失败: %s" % e}), 500
        if len(raw) < 512:
            return jsonify({"error": "生成的 PDF 太小（%d 字节）" % len(raw)}), 500
    else:
        data = body.get("data")
        if not isinstance(data, str) or not data:
            return jsonify({"error": "data 为空"}), 400
        if len(data) > 200 * 1024 * 1024:
            return jsonify({"error": "数据过大"}), 413
        import base64

        try:
            raw = base64.b64decode(data, validate=False)
        except Exception as e:                                   # noqa: BLE001
            return jsonify({"error": "base64 解码失败: %s" % e}), 400
        if len(raw) < 64:
            return jsonify({"error": "导出数据太小（%d 字节）" % len(raw)}), 400

    src = body.get("src") or ""
    if src and _exists(src):
        dest_dir = os.path.dirname(panel_lib._long(src)) or str(EXPORT_DIR)
        stem = os.path.splitext(os.path.basename(src))[0]
    else:
        dest_dir = str(EXPORT_DIR)
        stem = _safe_name(body.get("name") or "figure")
        stem = os.path.splitext(stem)[0]

    dpi = body.get("dpi")
    ext = ".jpg" if kind in ("jpg", "jpeg") else "." + kind
    tag = "_edit"
    if kind in ("png", "jpg", "jpeg") and dpi:
        tag += "_%ddpi" % int(dpi)

    os.makedirs(dest_dir, exist_ok=True)
    target = os.path.join(dest_dir, stem + tag + ext)
    n = 2
    while os.path.exists(target):
        target = os.path.join(dest_dir, "%s%s_%d%s" % (stem, tag, n, ext))
        n += 1
    try:
        with open(target, "wb") as f:
            f.write(raw)
    except OSError as e:
        return jsonify({"error": "写盘失败: %s" % e}), 500
    return jsonify({"ok": True, "path": target, "bytes": len(raw),
                    "human": _human(len(raw)), "kind": kind})


@app.route("/api/svg/history")
def api_svg_history():
    """某个 SVG 的备份历史（新的在前）。"""
    src = request.args.get("src", "")
    if not src:
        return jsonify({"error": "缺少 src"}), 400
    base = os.path.basename(src)
    out: List[Dict[str, Any]] = []
    try:
        if BACKUP_DIR.is_dir():
            for f in BACKUP_DIR.iterdir():
                nm = f.name
                if not nm.startswith(base + ".bak_"):
                    continue
                try:
                    st = f.stat()
                except OSError:
                    continue
                out.append({"file": nm, "bytes": st.st_size,
                            "mtime": st.st_mtime,
                            "ts": datetime.fromtimestamp(st.st_mtime).strftime("%Y-%m-%d %H:%M:%S")})
    except OSError:
        pass
    out.sort(key=lambda x: x["mtime"], reverse=True)
    return jsonify({"ok": True, "src": src, "backups": out[:50]})


@app.route("/api/svg/restore", methods=["POST"])
def api_svg_restore():
    """回滚到某个备份。body: ``{src, backup}``

    ⚠ 回滚前**先把当前版本再备份一次**，所以「回滚」本身也是可撤销的。
    """
    body = request.get_json(force=True, silent=True) or {}
    src = body.get("src", "")
    bak_name = body.get("backup", "")
    if not src or not _exists(src):
        return jsonify({"error": "file missing"}), 404
    # 只允许 BACKUP_DIR 里的文件名，防路径逃逸
    bak_name = os.path.basename(str(bak_name))
    bak = BACKUP_DIR / bak_name
    if not bak_name or not bak.is_file():
        return jsonify({"error": "备份不存在"}), 404
    lp = panel_lib._long(src)
    try:
        stamp = datetime.now().strftime("%Y%m%d_%H%M%S_%f")[:-3]
        pre = BACKUP_DIR / ("%s.bak_%s.svg" % (os.path.basename(src), stamp))
        shutil.copy2(lp, pre)
        shutil.copy2(str(bak), lp)
        size = os.path.getsize(lp)
    except OSError as e:
        return jsonify({"error": "回滚失败: %s" % e}), 500
    return jsonify({"ok": True, "src": src, "from": bak_name,
                    "pre_backup": str(pre), "bytes": size})


# ── svguest 会话：把「要精修哪些面板」存在服务端，hash 里只带一个短 token ──
SESSIONS: Dict[str, Dict[str, Any]] = {}
SESSION_TTL = 1800.0          # 30 分钟


@app.route("/api/svg/session", methods=["POST"])
def api_svg_session_create():
    """body: ``{items:[{src,name}]}`` → ``{token}``

    ⚠ 放在服务端而不是 hash 里：面板多时 URL 会超长，且中文路径要反复编码。
    """
    body = request.get_json(force=True, silent=True) or {}
    items = body.get("items")
    if not isinstance(items, list) or not items:
        return jsonify({"error": "items 为空"}), 400
    clean = []
    for it in items[:20]:                     # svguest 自身上限 20 个标签
        if not isinstance(it, dict):
            continue
        src = it.get("src") or ""
        if not src or not _exists(src) or not src.lower().endswith(".svg"):
            continue
        clean.append({"src": src,
                      "name": it.get("name") or os.path.basename(src)})
    if not clean:
        return jsonify({"error": "没有可用的 SVG 面板"}), 400
    token = "%x-%s" % (int(time.time() * 1000), os.urandom(4).hex())
    now = time.time()
    for k in [k for k, v in SESSIONS.items() if now - v["t"] > SESSION_TTL]:
        SESSIONS.pop(k, None)
    SESSIONS[token] = {"t": now, "items": clean}
    return jsonify({"ok": True, "token": token, "count": len(clean)})


@app.route("/api/svg/session/<token>")
def api_svg_session_get(token: str):
    s = SESSIONS.get(token)
    if not s:
        return jsonify({"error": "会话不存在或已过期"}), 404
    import urllib.parse

    items = [{"src": it["src"], "name": it["name"],
              "url": "/api/svg/raw?src=" + urllib.parse.quote(it["src"])}
             for it in s["items"]]
    return jsonify({"ok": True, "items": items,
                    "back": "/api/svg/writeback",
                    "backBatch": "/api/svg/writeback_batch",
                    "export": "/api/svg/export"})



@app.route("/api/probe", methods=["POST"])
def api_probe():
    body = request.get_json(force=True, silent=True) or {}
    src = body.get("src", "")
    if not src or not _exists(src):
        return jsonify({"error": "文件不存在"}), 404
    try:
        return jsonify(composer.probe(src, dpi=int(body.get("dpi", 300)),
                                      adjust=_parse_adjust(body.get("adjust"))))
    except Exception as e:                                       # noqa: BLE001
        return jsonify({"error": f"{type(e).__name__}: {e}"}), 500


# ══════════════════════════════════════════════════════════════════
# 调色辅助：主要颜色 / 项目调色板
# ══════════════════════════════════════════════════════════════════
_PALETTE_FALLBACK = {
    "name": "内置 OKABE_ITO",
    "path": None,
    "groups": [{"name": "OKABE_ITO（项目分类色·固定序）", "colors": [
        {"name": n, "hex": h} for n, h in zip(
            ["orange", "sky blue", "bluish green", "yellow", "blue",
             "vermillion", "reddish purple", "black", "grey"],
            ["#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
             "#D55E00", "#CC79A7", "#000000", "#999999"])]}],
    "note": "未找到项目 sci_palette.py，用内置副本",
}

_SEMANTIC_ATTRS = [
    ("C_UP", "上升 / 显著 / 观测值"),
    ("C_DOWN", "下降 / 不显著 / 期望值"),
    ("C_NEUTRAL", "参考线 / 灰底"),
    ("C_REF", "阈值线 / 参考线"),
    ("C_OBS", "观测统计量"),
    ("C_EXP", "随机化 / 期望分布"),
]


def _load_project_palette() -> Dict[str, Any]:
    """从项目 ``utils/sci_palette.py`` 取分类色与语义色。

    ⚠ 必须按**显式路径**用 ``spec_from_file_location`` 加载（项目纪律）：
    该模块名 ``sci_palette`` 在多个管线里重名，走 ``import`` 会命中错的那份。
    """
    import importlib.util

    for p in (HERE.parent / "virome_phylo_pipeline" / "utils" / "sci_palette.py",
              HERE.parent / "virome_discovery_pipeline" / "utils" / "sci_palette.py"):
        if not p.is_file():
            continue
        try:
            spec = importlib.util.spec_from_file_location("_fs_sci_palette", str(p))
            mod = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(mod)                          # type: ignore[union-attr]
        except Exception:                                        # noqa: BLE001
            continue
        groups: List[Dict[str, Any]] = []
        oi = getattr(mod, "OKABE_ITO", None)
        if oi:
            groups.append({"name": "OKABE_ITO（项目分类色·固定序）", "colors": [
                {"name": f"OI{i + 1}", "hex": c} for i, c in enumerate(oi)]})
        ob = getattr(mod, "OKABE_BRANCH", None)
        if ob:
            groups.append({"name": "OKABE_BRANCH（树 / 分支）", "colors": [
                {"name": f"OB{i + 1}", "hex": c} for i, c in enumerate(ob)]})
        sem = [{"name": f"{attr} · {label}", "hex": getattr(mod, attr)}
               for attr, label in _SEMANTIC_ATTRS if getattr(mod, attr, None)]
        if sem:
            groups.append({"name": "语义色（固定含义 · 跨图一致）", "colors": sem})
        if not groups:
            continue
        return {"name": "项目 sci_palette", "path": str(p), "groups": groups,
                "seq_cmap": getattr(mod, "SEQ_CMAP", None),
                "div_cmap": getattr(mod, "DIV_CMAP", None),
                "note": "语义色未登记会静默退化为 #999999 —— 换色时别把语义色当普通色板用"}
    return dict(_PALETTE_FALLBACK)


@app.route("/api/palette")
def api_palette():
    """项目调色板（分类色 / 语义色），供前端做颜色替换与着色取色。"""
    try:
        return jsonify(_load_project_palette())
    except Exception as e:                                       # noqa: BLE001
        traceback.print_exc()
        return jsonify({"error": f"{type(e).__name__}: {e}"}), 500


@app.route("/api/colors", methods=["POST"])
def api_colors():
    """面板的主要颜色（按占比降序）—— 颜色替换表的候选 from 色。

    body: ``{src, n?, dpi?, adjust?}``。默认**不施加 adjust**：
    remap 要在原图色值上精确匹配，所以候选色取的是源图原始色。
    """
    body = request.get_json(force=True, silent=True) or {}
    src = body.get("src", "")
    if not src or not _exists(src):
        return jsonify({"error": "文件不存在"}), 404
    try:
        res = composer.dominant_colors(
            src, n=int(body.get("n", 12)), dpi=int(body.get("dpi", 300)),
            adjust=_parse_adjust(body.get("adjust")))
        res["mode"] = composer.probe(src, dpi=int(body.get("dpi", 300)))["raw_mode"]
        if res.get("is_svg"):
            try:
                import svg_edit

                ins = svg_edit.inspect_svg(src)
                res["legend"] = ins.get("legend", [])
                res["groups"] = ins.get("groups", [])[:60]
                res["texts"] = ins.get("texts", [])[:200]
                res["fonts"] = ins.get("fonts", {})
                res["flavor"] = ins.get("flavor", "generic")
            except Exception:                                    # noqa: BLE001
                pass
        return jsonify(res)
    except Exception as e:                                       # noqa: BLE001
        traceback.print_exc()
        return jsonify({"error": f"{type(e).__name__}: {e}"}), 500


# ══════════════════════════════════════════════════════════════════
# 渲染
# ══════════════════════════════════════════════════════════════════
@app.route("/api/render", methods=["POST"])
def api_render():
    body = request.get_json(force=True, silent=True) or {}
    layout = body.get("layout")
    if not isinstance(layout, dict):
        return jsonify({"error": "缺少 layout"}), 400
    stem = _safe_name(body.get("stem") or "figure")
    formats = tuple(body.get("formats") or ("png", "pdf"))
    dpi = body.get("dpi")
    try:
        res = composer.render_layout(
            layout, str(OUTPUT_DIR), stem, dpi=dpi, formats=formats)
    except Exception as e:                                       # noqa: BLE001
        traceback.print_exc()
        return jsonify({"error": f"{type(e).__name__}: {e}"}), 500

    def _url(p: Optional[str]) -> Optional[str]:
        if not p:
            return None
        rel = os.path.relpath(p, str(OUTPUT_DIR)).replace("\\", "/")
        return f"/api/output/{rel}"

    res["urls"] = {"png": _url(res.get("png")), "pdf": _url(res.get("pdf"))}
    res["stem"] = stem
    return jsonify(res)


@app.route("/api/output/<path:rel>")
def api_output(rel: str):
    target = (OUTPUT_DIR / rel).resolve()
    if not str(target).startswith(str(OUTPUT_DIR.resolve())) or not target.is_file():
        abort(404)
    inline = request.args.get("download") != "1"
    return send_file(str(target), as_attachment=not inline)


@app.route("/api/outputs")
def api_outputs():
    items = []
    for f in sorted(OUTPUT_DIR.glob("*"), key=lambda p: p.stat().st_mtime,
                    reverse=True):
        if f.is_file() and f.suffix.lower() in (".png", ".pdf"):
            st = f.stat()
            items.append({"name": f.name, "url": f"/api/output/{f.name}",
                          "size_kb": round(st.st_size / 1024.0, 1),
                          "mtime": st.st_mtime})
    return jsonify({"items": items[:200]})


# ══════════════════════════════════════════════════════════════════
# 布局存取
# ══════════════════════════════════════════════════════════════════
@app.route("/api/layouts", methods=["GET"])
def api_layouts_list():
    items = []
    for f in sorted(LAYOUT_DIR.glob("*.json"),
                    key=lambda p: p.stat().st_mtime, reverse=True):
        try:
            with open(f, "r", encoding="utf-8") as fh:
                d = json.load(fh)
            n = len(d.get("items") or [])
        except Exception:                                        # noqa: BLE001
            n = -1
        items.append({"name": f.stem, "items": n,
                      "mtime": f.stat().st_mtime})
    return jsonify({"items": items})


@app.route("/api/layouts/<name>", methods=["GET", "POST", "DELETE"])
def api_layouts_item(name: str):
    safe = _safe_name(name)
    f = LAYOUT_DIR / f"{safe}.json"
    if request.method == "GET":
        if not f.is_file():
            return jsonify({"error": "不存在"}), 404
        with open(f, "r", encoding="utf-8") as fh:
            return jsonify(json.load(fh))
    if request.method == "DELETE":
        if f.is_file():
            f.unlink()
        return jsonify({"ok": True})
    body = request.get_json(force=True, silent=True)
    if not isinstance(body, dict):
        return jsonify({"error": "无效 JSON"}), 400
    with open(f, "w", encoding="utf-8") as fh:
        json.dump(body, fh, ensure_ascii=False, indent=2)
    return jsonify({"ok": True, "name": safe, "path": str(f)})


@app.route("/api/cache/clear", methods=["POST"])
def api_cache_clear():
    n = panel_lib.clear_thumbs(str(CACHE_DIR))
    return jsonify({"ok": True, "removed": n})


# ══════════════════════════════════════════════════════════════════
# AI 绘图（svguest 的 SVG Chat）——服务端代理 + 密钥托管
# ══════════════════════════════════════════════════════════════════
# ⚠⚠ 为什么走服务端代理（**实测过 CORS 头，别想当然**）：
#   ① GLM（open.bigmodel.cn）：OPTIONS 预检**不回 Access-Control-Allow-Origin**
#      （虽然 POST 的 401 响应带了），浏览器预检阶段就会被拦 → 直连不可用。
#      这正是上游只给 GLM 配 Cloudflare 中继、而 DeepSeek 直接 send() 的原因；
#      但那个中继的 Worker 只放行 jieqianghe.github.io / svguest.online，
#      127.0.0.1 不在白名单里，所以本平台用不了它。
#   ② DeepSeek（api.deepseek.com）：预检 200 且 ACAO **回显任意 origin**
#      → 浏览器其实可以直连。走代理是为了**密钥不落浏览器**（只存本机磁盘）
#      和统一错误呈现，不是被 CORS 逼的。
#   ⇒ 代理对 GLM 是硬需求，对 DeepSeek 是可选但更安全的路径。
# 上游把「发请求」收敛在一个函数 apiPost(p, url, key, body, o) 上，返回标准
# Response，调用方用 r.body 读 SSE —— 所以桥接只要替换 window.apiPost，
# 上游的流式解析、重试、错误提示全都无感复用。
AI_KEY_FILE = HERE / ".ai_keys.json"
AI_KEY_FILE.parent.mkdir(parents=True, exist_ok=True)
# 只允许这两个 API 主机（防 SSRF：不能让本机变成任意内网请求的跳板）
AI_HOSTS = {"api.deepseek.com", "open.bigmodel.cn"}
AI_PROVIDERS = {"Ds": "DeepSeek", "Glm": "GLM"}


def _load_ai_keys() -> Dict[str, str]:
    try:
        with open(AI_KEY_FILE, "r", encoding="utf-8") as fh:
            d = json.load(fh)
        return {k: v for k, v in d.items() if isinstance(v, str) and v}
    except (OSError, ValueError):
        return {}


def _mask(k: str) -> str:
    """密钥只回显可辨识的部分，绝不明文返回给前端。"""
    if not k:
        return ""
    if len(k) <= 10:
        return "*" * len(k)
    return k[:4] + "*" * 6 + k[-4:]


def _check_ai_url(url: str) -> Optional[str]:
    """校验目标 URL，返回错误信息（None = 通过）。"""
    from urllib.parse import urlsplit

    try:
        u = urlsplit(url)
    except ValueError:
        return "URL 无法解析"
    if u.scheme != "https":
        return "只允许 https"
    if u.username or u.password:
        return "URL 里不能带账号密码"
    if (u.hostname or "").lower() not in AI_HOSTS:
        return "不允许的主机：%s（只支持 %s）" % (u.hostname, "、".join(sorted(AI_HOSTS)))
    return None


@app.route("/api/ai/config", methods=["GET", "POST"])
def api_ai_config():
    """服务端托管的 AI 密钥。GET 只返回掩码，POST 写入/清除。

    ⚠ 文件是本机 ``.ai_keys.json``（已在 .gitignore 之外独立存放，不要提交）。
    """
    keys = _load_ai_keys()
    if request.method == "GET":
        return jsonify({
            "ok": True,
            "providers": {
                p: {"name": n, "configured": bool(keys.get(p)), "hint": _mask(keys.get(p, ""))}
                for p, n in AI_PROVIDERS.items()
            },
            "hosts": sorted(AI_HOSTS),
            "file": str(AI_KEY_FILE),
        })
    body = request.get_json(force=True, silent=True) or {}
    p = str(body.get("provider") or "").strip()
    if p not in AI_PROVIDERS:
        return jsonify({"error": "provider 只支持 Ds / Glm"}), 400
    key = body.get("key")
    if key is None:
        return jsonify({"error": "缺少 key"}), 400
    key = str(key).strip()
    if key and len(key) < 8:
        return jsonify({"error": "key 太短，像是没填对"}), 400
    if key:
        keys[p] = key
    else:
        keys.pop(p, None)
    try:
        with open(AI_KEY_FILE, "w", encoding="utf-8") as fh:
            json.dump(keys, fh, ensure_ascii=False, indent=2)
        try:                                   # 尽力收紧权限（Windows 上无效果但不报错）
            os.chmod(AI_KEY_FILE, 0o600)
        except OSError:
            pass
    except OSError as e:
        return jsonify({"error": "写入失败: %s" % e}), 500
    return jsonify({"ok": True, "provider": p, "configured": bool(key),
                    "hint": _mask(keys.get(p, ""))})


def _ai_err(msg: str, code: str, status: int, **extra: Any):
    """AI 相关错误的统一格式。

    ⚠ 必须是 ``{"error": {"message": ...}}`` 而不是 ``{"error": "字符串"}``：
      svguest 的 ``apiMsg()`` 解析的是 ``j.error.message``，给它一个字符串
      会把整段 JSON 原样显示给用户（``{"error":"..."}``）。
    """
    body: Dict[str, Any] = {"error": {"message": msg, "code": code}}
    body.update(extra)
    return jsonify(body), status


@app.route("/api/ai/chat", methods=["POST"])
def api_ai_chat():
    """把 svguest 的 LLM 请求转发给上游 API，**流式透传**响应。

    body: ``{provider:'Ds'|'Glm', url, key?, body}``
    - ``key`` 省略时用服务端托管的密钥（浏览器里就不需要存 key 了）。
    - 返回的是上游原样的字节流：非流式是 JSON，流式是 ``text/event-stream``。
      ⚠ 绝不能缓冲成一次性 JSON —— 上游用 ``r.body`` 逐块读 SSE 做打字机效果。
    """
    import requests

    body = request.get_json(force=True, silent=True) or {}
    provider = str(body.get("provider") or "Ds")
    if provider not in AI_PROVIDERS:
        return _ai_err("provider 只支持 Ds / Glm", "bad_provider", 400)
    url = str(body.get("url") or "")
    err = _check_ai_url(url)
    if err:
        return _ai_err(err, "bad_url", 403)
    payload = body.get("body")
    if not isinstance(payload, dict):
        return _ai_err("body 必须是对象", "bad_body", 400)
    key = str(body.get("key") or "").strip() or _load_ai_keys().get(provider, "")
    if not key:
        return _ai_err(
            "还没有配置 %s 的 API key。可以在工作台右侧「AI 绘图」里填，"
            "或在 svguest 的 Models and keys 里填。" % AI_PROVIDERS[provider],
            "no_key", 401, no_key=True)

    try:
        up = requests.post(
            url,
            headers={"Content-Type": "application/json",
                     "Authorization": "Bearer " + key},
            json=payload,
            stream=True,
            timeout=(15, 300),          # 连接 15s，读 300s（长回复不能中途断）
        )
    except requests.RequestException as e:
        return _ai_err("连不上 %s：%s" % (url.split("/")[2], e), "upstream_unreachable", 502)

    ctype = up.headers.get("Content-Type", "application/json")

    def gen():
        try:
            for chunk in up.iter_content(chunk_size=8192):
                if chunk:
                    yield chunk
        finally:
            up.close()

    resp = Response(stream_with_context(gen()), status=up.status_code,
                    content_type=ctype)
    # 逐块下发：不要 Content-Length（否则浏览器会等全量），也不要任何缓存
    resp.headers["Cache-Control"] = "no-cache, no-store, must-revalidate"
    resp.headers["X-Accel-Buffering"] = "no"
    return resp


# ══════════════════════════════════════════════════════════════════
# 磁盘自动快照（Auto-save 的服务端版本）
# ══════════════════════════════════════════════════════════════════
# 上游自带 Auto-save 存的是**浏览器 IndexedDB**（sgDocs，默认开启，10 份/200 MB）。
# 那是崩溃安全网，但换机器、清缓存就没了。这里再补一层**落盘**快照：
#   - 只在编辑停顿时写，**绝不覆盖原文件**（写回仍然只由 Ctrl+S 触发）；
#   - 每源保留最近 AUTOSAVE_KEEP 份，自动剪枝；
#   - 放在 outputs/ 下而不是面板目录（否则 .svg 会被面板库扫成面板）。
AUTOSAVE_DIR = OUTPUT_DIR / "_svguest_autosave"
AUTOSAVE_KEEP = 12


def _autosave_slug(src: str) -> str:
    """给源文件一个稳定的短目录名（避开中文/超长路径，且同名不同目录不冲突）。"""
    import hashlib

    h = hashlib.sha1(os.path.abspath(src).encode("utf-8")).hexdigest()[:10]
    base = re.sub(r"[^0-9A-Za-z._-]+", "_", os.path.basename(src))[:60] or "figure.svg"
    return "%s__%s" % (base, h)


@app.route("/api/svg/autosave", methods=["POST"])
def api_svg_autosave():
    """body: ``{src, svg}`` → 存一份快照，返回 ``{ok, file, n, keep}``。"""
    body = request.get_json(force=True, silent=True) or {}
    src = str(body.get("src") or "")
    svg = body.get("svg")
    if not src:
        return jsonify({"error": "缺少 src"}), 400
    if not isinstance(svg, str) or "<svg" not in svg:
        return jsonify({"error": "empty svg"}), 400
    if len(svg) < 64:
        return jsonify({"error": "内容太小（%d 字节），不存" % len(svg)}), 400
    d = AUTOSAVE_DIR / _autosave_slug(src)
    try:
        os.makedirs(str(d), exist_ok=True)
        stamp = datetime.now().strftime("%Y%m%d_%H%M%S_%f")[:-3]
        f = d / ("%s.svg" % stamp)
        with open(f, "wb") as fh:
            fh.write(svg.encode("utf-8"))
        # 剪枝：只留最近 N 份
        olds = sorted((p for p in d.iterdir() if p.suffix == ".svg"),
                      key=lambda p: p.stat().st_mtime, reverse=True)
        for p in olds[AUTOSAVE_KEEP:]:
            try:
                p.unlink()
            except OSError:
                pass
        n = min(len(olds), AUTOSAVE_KEEP)
    except OSError as e:
        return jsonify({"error": "快照写入失败: %s" % e}), 500
    return jsonify({"ok": True, "src": src, "file": f.name, "dir": str(d),
                    "bytes": len(svg), "n": n, "keep": AUTOSAVE_KEEP,
                    "ts": datetime.now().strftime("%H:%M:%S")})


@app.route("/api/svg/autosave")
def api_svg_autosave_list():
    """``?src=`` → 该源的快照列表（新的在前）。"""
    src = request.args.get("src", "")
    if not src:
        return jsonify({"error": "缺少 src"}), 400
    d = AUTOSAVE_DIR / _autosave_slug(src)
    out: List[Dict[str, Any]] = []
    try:
        if d.is_dir():
            for f in d.iterdir():
                if f.suffix != ".svg":
                    continue
                try:
                    st = f.stat()
                except OSError:
                    continue
                out.append({"file": f.name, "bytes": st.st_size, "mtime": st.st_mtime,
                            "ts": datetime.fromtimestamp(st.st_mtime).strftime("%Y-%m-%d %H:%M:%S")})
    except OSError:
        pass
    out.sort(key=lambda x: x["mtime"], reverse=True)
    return jsonify({"ok": True, "src": src, "snapshots": out, "keep": AUTOSAVE_KEEP})


@app.route("/api/svg/autosave_restore", methods=["POST"])
def api_svg_autosave_restore():
    """body: ``{src, file}`` → 把快照写回原文件（写回前仍会进 BACKUP_DIR 备份）。"""
    body = request.get_json(force=True, silent=True) or {}
    src = str(body.get("src") or "")
    name = os.path.basename(str(body.get("file") or ""))
    if not src or not _exists(src):
        return jsonify({"error": "file missing"}), 404
    snap = AUTOSAVE_DIR / _autosave_slug(src) / name
    if not name or not snap.is_file():
        return jsonify({"error": "快照不存在"}), 404
    try:
        svg = snap.read_text(encoding="utf-8")
    except OSError as e:
        return jsonify({"error": "读快照失败: %s" % e}), 500
    res = _write_svg_back(src, svg)
    if "error" in res:
        return jsonify(res), res.pop("_code", 400)
    res["from_snapshot"] = name
    return jsonify(res)


def main(argv: Optional[List[str]] = None) -> int:
    import argparse

    ap = argparse.ArgumentParser(description="Figure Studio 本地组图平台")
    ap.add_argument("--host", default="127.0.0.1")
    ap.add_argument("--port", type=int, default=8770)
    ap.add_argument("--debug", action="store_true")
    ap.add_argument("--root", default=None, help="启动时默认扫描目录")
    args = ap.parse_args(argv)

    if args.root:
        cfg = load_config()
        cfg["last_root"] = os.path.abspath(args.root)
        if cfg["last_root"] not in cfg["roots"]:
            cfg["roots"].insert(0, cfg["last_root"])
        save_config(cfg)

    print(f"Figure Studio -> http://{args.host}:{args.port}")
    app.run(host=args.host, port=args.port, debug=args.debug,
            use_reloader=False, threaded=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
