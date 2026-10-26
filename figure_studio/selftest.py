# -*- coding: utf-8 -*-
"""Figure Studio 端到端自检。

真实面板 → 按资产比反解排版 → 走后端 API 渲染 → 校验产物。

跑法::

    C:/Python312/python.exe figure_studio/selftest.py

校验项
------
1. 面板库扫描能返回真实面板（含 PNG 与 PDF 两种源）
2. 排版后**每块面板比 == 源图比**（零形变）
3. 导出 PNG 像素尺寸 == canvas(in) × dpi
4. 后端返回的 ``warnings`` 为空
5. PDF 产物存在且非空
"""
from __future__ import annotations

import json
import sys
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

from PIL import Image
import numpy as np

HERE = Path(__file__).resolve().parent
SVG_PANEL = str(HERE / "outputs/_svgdemo/panel_net.svg")
sys.path.insert(0, str(HERE))

import composer  # noqa: E402
import panels as panel_lib  # noqa: E402

BASE = "http://127.0.0.1:8770"
PANEL_DIR = (r"D:\桌面\延伸基因组\MMPV-RNA\reports\_archive\review_dirs"
             r"\_discovery_sci_20261012\local_test\panels")

DPI = 300
MARGIN = 0.10
TOP = 0.90          # 顶部给总标题留出空间
GAP = 0.12
CAP_ALLOW = 0.24    # ⚠ 每行题注占的高度；不预留 → 末行题注被画布下缘裁掉（零报错）
ROWS = [["A", "B", "C"], ["D", "L"]]      # 第 2 行放一个 PDF 源，验证 fitz 路径

ok_all = True


def check(label: str, cond: bool, detail: str = "") -> None:
    global ok_all
    ok_all = ok_all and cond
    print(("  [PASS] " if cond else "  [FAIL] ") + label + ("  " + detail if detail else ""))


def post(path: str, payload: dict):
    req = urllib.request.Request(
        BASE + path, data=json.dumps(payload).encode("utf-8"),
        headers={"Content-Type": "application/json"}, method="POST")
    with urllib.request.urlopen(req, timeout=600) as r:
        return json.loads(r.read().decode("utf-8"))


def get(path: str):
    with urllib.request.urlopen(BASE + path, timeout=600) as r:
        return json.loads(r.read().decode("utf-8"))


def solve_row(items, avail_w, gap):
    """行内按资产比反解：等高 h，w_i = h / aspect_i，Σw = avail_w - gap*(n-1)。"""
    s = sum(1.0 / it["aspect"] for it in items)
    h = (avail_w - gap * (len(items) - 1)) / s
    out = []
    x = 0.0
    for it in items:
        w = h / it["aspect"]
        out.append((it, x, w, h))
        x += w + gap
    return out


def _ai_proxy_stream_test() -> None:
    """AI 代理必须**逐块透传** SSE。

    上游用 `r.body` 读 SSE 做打字机效果。若服务端攒成一坨，功能不报错，
    只是「等 N 秒然后全文一次性蹦出来」—— 零报错的功能退化。

    ⚠ 假上游必须用 `Transfer-Encoding: chunked`（HTTP/1.1）。真实 LLM API 都是
      chunked；若用 HTTP/1.0「靠关连接表示结束」，`read(amt)` 会阻塞到 EOF，
      测出来必然像「被缓冲」—— 那是**夹具造成的假阳性**（实测踩过）。
    """
    import threading
    import time
    from http.server import BaseHTTPRequestHandler, HTTPServer

    port = 18925

    class Fake(BaseHTTPRequestHandler):
        protocol_version = "HTTP/1.1"

        def do_POST(self):                                   # noqa: N802
            n = int(self.headers.get("Content-Length") or 0)
            self.rfile.read(n)
            self.send_response(200)
            self.send_header("Content-Type", "text/event-stream")
            self.send_header("Transfer-Encoding", "chunked")
            self.end_headers()

            def put(b):
                self.wfile.write(("%x\r\n" % len(b)).encode() + b + b"\r\n")
                self.wfile.flush()

            for i in range(5):
                put(("data: " + json.dumps(
                    {"choices": [{"delta": {"content": "tok%d" % i}}]}) + "\n\n").encode())
                time.sleep(0.3)
            put(b"data: [DONE]\n\n")
            self.wfile.write(b"0\r\n\r\n")
            self.wfile.flush()

        def log_message(self, *a):
            pass

    srv = HTTPServer(("127.0.0.1", port), Fake)
    threading.Thread(target=srv.serve_forever, daemon=True).start()
    time.sleep(0.3)
    try:
        import app as A                                      # noqa: E402

        A.AI_HOSTS.add("127.0.0.1")
        orig = A._check_ai_url
        A._check_ai_url = lambda u: None if "127.0.0.1" in u else orig(u)
        try:
            c = A.app.test_client()
            t0 = time.time()
            r = c.post("/api/ai/chat", json={
                "provider": "Ds", "url": "http://127.0.0.1:%d/v1" % port,
                "key": "sk-test",
                "body": {"model": "x", "messages": [], "stream": True},
            }, buffered=False)
            arrivals, body = [], b""
            for ch in r.response:
                arrivals.append(time.time() - t0)
                body += ch
        finally:
            A._check_ai_url = orig

        spread = (arrivals[-1] - arrivals[0]) if len(arrivals) > 1 else 0.0
        check("AI proxy returns text/event-stream",
              (r.headers.get("Content-Type") or "").startswith("text/event-stream"),
              r.headers.get("Content-Type") or "")
        check("AI proxy streams chunk by chunk (not buffered)",
              len(arrivals) >= 3 and spread > 0.6,
              "%d chunks, spread %.1fs" % (len(arrivals), spread))
        check("AI proxy passes [DONE] through", b"[DONE]" in body,
              "%d bytes" % len(body))

        # 白名单外必须拒（防 SSRF：别把本机变成内网跳板）
        r403 = c.post("/api/ai/chat", json={
            "provider": "Ds", "url": "https://evil.example.com/x", "body": {}})
        check("AI proxy rejects non-whitelisted host",
              r403.status_code == 403, "status=%d" % r403.status_code)
        # 非 https 也拒
        r400 = c.post("/api/ai/chat", json={
            "provider": "Ds", "url": "http://api.deepseek.com/x", "body": {}})
        check("AI proxy rejects non-https", r400.status_code == 403,
              "status=%d" % r400.status_code)
    finally:
        srv.shutdown()


def main() -> int:
    print("== 1. 扫描面板库 ==")
    try:
        res = get(f"/api/panels?root={urllib.parse.quote(PANEL_DIR)}&limit=50")
    except urllib.error.URLError as e:
        print("  服务未启动，请先跑 app.py：", e)
        return 2
    found = {(Path(p["src"]).stem, p["ext"]): p for p in res["panels"]}
    check("扫描到面板", len(found) >= 6, f"{len(found)} 个")
    check("含 PNG 源", any(k[1] == ".png" for k in found))
    check("含 PDF 源", any(k[1] == ".pdf" for k in found))

    # 每个字母优先用 PNG，行 2 的 L 用 PDF（测 fitz）
    srcs = {}
    for row in ROWS:
        for L in row:
            ext = ".pdf" if (L == "L") else ".png"
            p = found.get((f"panel_{L}", ext))
            if p:
                srcs[L] = p
    check("取齐所需面板", len(srcs) == sum(len(r) for r in ROWS),
          f"{sorted(srcs)}")

    print("== 2. 排版（按资产比反解）==")
    items = []
    y = TOP
    for row in ROWS:
        row_items = [dict(letter=L, **srcs[L]) for L in row if L in srcs]
        if not row_items:
            continue
        avail = 10.0 - 2 * MARGIN
        for it, x, w, h in solve_row(row_items, avail, GAP):
            items.append({
                "id": "t_" + it["letter"], "src": it["src"],
                "x_in": round(MARGIN + x, 4), "y_in": round(y, 4),
                "w_in": round(w, 4), "h_in": round(h, 4),
                "letter": it["letter"],
                "caption": f"panel {it['letter']} · {it['px_w']}×{it['px_h']} px",
                "z": 0, "_aspect": it["aspect"],
            })
        y += items[-1]["h_in"] + GAP + CAP_ALLOW
    canvas_h = round(y - GAP + MARGIN, 4)

    worst = 0.0
    for it in items:
        dev = abs(it["h_in"] / it["w_in"] - it["_aspect"]) / it["_aspect"]
        worst = max(worst, dev)
    # 阈值 2e-4：受「几何量保留 4 位小数」取整限制（0.00005 in / ~0.5 in ≈ 1e-4）。
    # 真正的形变防线是后端渲染时的 1% 告警（见 composer.render_layout）。
    check("排版零形变（逐块面板比 == 源图比）", worst < 2e-4,
          f"最大偏差 {worst:.2e}（取整限）")

    print("== 3. 后端渲染 ==")
    layout = {
        "canvas": {"width_in": 10.0, "height_in": canvas_h, "dpi": DPI},
        "style": {"bg": "white", "frame": True, "frame_color": "#bbbbbb",
                  "frame_lw": 0.6, "letter_fs": 11, "letter_dx": 0.0,
                  "letter_dy": 1.02, "caption_fs": 7, "caption_gap_in": 0.06},
        "title": {"text": "Figure Studio selftest", "subtitle":
                  "aspect-preserving composite", "fs": 13, "sub_fs": 9.5,
                  "x_in": MARGIN, "y_in": 0.34},
        "items": [{k: v for k, v in it.items() if not k.startswith("_")}
                  for it in items],
    }
    (HERE / "layouts").mkdir(exist_ok=True)
    (HERE / "layouts" / "_selftest.json").write_text(
        json.dumps(layout, ensure_ascii=False, indent=2), encoding="utf-8")

    res = post("/api/render", {"layout": layout, "stem": "_selftest_fig",
                               "formats": ["png", "pdf"], "dpi": DPI})
    if res.get("error"):
        check("渲染返回", False, res["error"])
        return 1
    check("无形变/读取告警", not res["warnings"], str(res["warnings"]))

    png = res.get("png")
    pdf = res.get("pdf")
    check("PNG 产物存在", bool(png) and Path(png).is_file())
    check("PDF 产物存在", bool(pdf) and Path(pdf).is_file()
          and Path(pdf).stat().st_size > 1000)

    if png and Path(png).is_file():
        with Image.open(png) as im:
            exp_w, exp_h = 10.0 * DPI, canvas_h * DPI
            check("PNG 尺寸 == canvas × dpi（±1 px 取整）",
                  abs(im.size[0] - exp_w) <= 1 and abs(im.size[1] - exp_h) <= 1,
                  f"实际 {im.size} / 期望 {exp_w:.1f}×{exp_h:.1f}")
            check("后端报告尺寸 == 文件实际尺寸",
                  tuple(res["size_px"]) == im.size,
                  f"报告 {res['size_px']} / 文件 {list(im.size)}")
            check("PNG 为 RGB（非调色板）", im.mode == "RGB", im.mode)

    print("== 4. 布局往返 ==")
    post("/api/layouts/_selftest_rt", layout)
    back = get("/api/layouts/_selftest_rt")
    check("布局存取一致", len(back.get("items", [])) == len(layout["items"]))

    print("== 5. 调色（adjust）端到端 ==")
    dc = composer.dominant_colors(srcs[ROWS[0][0]]["src"], n=5)
    check("主要颜色提取（前 5）",
          len(dc["colors"]) == 5
          and all(c["hex"].startswith("#") and len(c["hex"]) == 7 for c in dc["colors"]),
          " / ".join(c["hex"] for c in dc["colors"]))

    # 只给第 1 块加灰度 → 渲染后裁出该块像素区，验证 R==G==B，
    # 并确认其他块仍是彩色（证明是**逐块**生效，不是整图被改）
    base_items = [{k: v for k, v in it.items() if not k.startswith("_")} for it in items]
    adj_items = [dict(x) for x in base_items]
    adj_items[0]["adjust"] = {"grayscale": True}
    layout_adj = dict(layout)
    layout_adj["items"] = adj_items
    res2 = post("/api/render", {"layout": layout_adj, "stem": "_selftest_adj",
                                "formats": ["png"], "dpi": DPI})
    if res2.get("error"):
        check("带 adjust 渲染返回", False, res2["error"])
    else:
        check("带 adjust 渲染无告警", not res2.get("warnings"), str(res2.get("warnings")))
    png2 = res2.get("png")
    if png2 and Path(png2).is_file():

        def _crop(im, it):
            x0 = int(round(it["x_in"] * DPI)) + 2
            y0 = int(round(it["y_in"] * DPI)) + 2
            x1 = int(round((it["x_in"] + it["w_in"]) * DPI)) - 2
            y1 = int(round((it["y_in"] + it["h_in"]) * DPI)) - 2
            return np.asarray(im.convert("RGB").crop((x0, y0, x1, y1)))

        with Image.open(png2) as im2:
            a = _crop(im2, adj_items[0])
            check("灰度面板渲染后 R==G==B",
                  bool((a[:, :, 0] == a[:, :, 1]).all() and (a[:, :, 1] == a[:, :, 2]).all()))
            b = _crop(im2, adj_items[1])
            check("其他面板仍为彩色（逐块生效）",
                  bool(((b[:, :, 0] != b[:, :, 1]) | (b[:, :, 1] != b[:, :, 2])).any()))

    print("== 6. SVG vector editing end-to-end ==")
    svg_src = SVG_PANEL
    if not Path(svg_src).is_file():
        check("SVG test panel exists", False, svg_src)
    else:
        import svg_edit

        ins = svg_edit.inspect_svg(svg_src)
        check("SVG inspect ok", bool(ins.get("ok")))
        check("SVG size parsed (in)",
              bool(ins.get("size_in")) and ins["size_in"][0] > 0,
              str(ins.get("size_in")))
        check("SVG semantic colours non-empty", len(ins.get("colors") or []) >= 3,
              "%d colours" % len(ins.get("colors") or []))
        check("SVG text nodes non-empty", (ins.get("n_text") or 0) > 0,
              str(ins.get("n_text")))
        check("SVG legend group detected", len(ins.get("legend") or []) >= 1,
              str([g.get("id") for g in (ins.get("legend") or [])]))

        # ── 物理尺寸三处口径必须一致 ────────────────────────────────
        #   fitz 把无单位 width="920"（px）当成 920pt → px 系 SVG 虚胖 4/3。
        #   唯一真相源是 svg_edit.svg_size_pt；meta / load_panel / inspect
        #   任何一处漂移都会让「1:1」和面板库信息对不上（且零报错）。
        for _f in ("panel_net.svg", "panel_tracks.svg", "panel_pca.svg",
                   "panel_map.svg"):
            _p = HERE / "outputs" / "_svgdemo" / _f
            if not _p.is_file():
                continue
            _m = panel_lib._meta_of(str(_p), DPI)
            _lp = composer.load_panel(str(_p), dpi=DPI)
            _w, _h = svg_edit.svg_size_pt(str(_p))
            check("SVG size 口径一致（meta/load/pt） " + _f,
                  _m is not None and _w
                  and abs(_m["w_in"] - _w / 72.0) < 0.01
                  and abs(_lp["w_in"] - _w / 72.0) < 0.01
                  and abs(_lp["h_in"] - _h / 72.0) < 0.01,
                  "meta %.3f / load %.3f / pt %.3f in"
                  % (_m["w_in"] if _m else -1, _lp["w_in"], _w / 72.0 if _w else -1))

        before = composer.load_panel(svg_src, dpi=200)["img"]
        after = composer.load_panel(svg_src, dpi=200,
                                    adjust={"svg": {"colors": [
                                        {"from": "#0f4d92", "to": "#D55E00"}]}})["img"]
        a0 = np.asarray(before)
        a1 = np.asarray(after)

        def _cnt(arr, hexv):
            r, g, b = (int(hexv[i:i + 2], 16) for i in (1, 3, 5))
            return int(((arr[:, :, 0] == r) & (arr[:, :, 1] == g)
                        & (arr[:, :, 2] == b)).sum())

        n_before = _cnt(a0, "#D55E00")
        n_after = _cnt(a1, "#D55E00")
        check("SVG vector recolour hits target colour",
              n_after > n_before + 500,
              "target-colour px %d -> %d" % (n_before, n_after))
        check("SVG recolour is exact (no anti-alias residue)",
              _cnt(a1, "#0f4d92") == 0, "residue %d px" % _cnt(a1, "#0f4d92"))

        hid = np.asarray(svg_edit.render_svg(
            svg_src, dpi=200, edit={"legend": {"action": "hide"}})["img"])
        check("SVG legend hide changes the picture",
              int((hid != a0).sum()) > 2000,
              "diff px %d" % int((hid != a0).sum()))

        fsc = np.asarray(svg_edit.render_svg(
            svg_src, dpi=200, edit={"font_scale": 1.5})["img"])
        check("SVG font scale changes the picture",
              int((fsc != a0).sum()) > 2000,
              "diff px %d" % int((fsc != a0).sum()))

        _first = next((t["text"] for t in ins["texts"] if t.get("text")), None)
        _txt_ok, _xml_hit = False, False
        if _first:
            _out, _notes = svg_edit.apply_svg_edit(
                svg_src, {"texts": [{"from": _first, "to": "ZZQQXX"}]})
            _txt_ok = _out is not None
            # 断言落到 DOM（不是"渲染没崩"就算过）
            _xml_hit = bool(_out) and b"ZZQQXX" in _out
        check("SVG text replace lands in DOM", _xml_hit, repr(_first)[:40])

        # 字体族：换成 fitz 确有差异的 serif（Arial 与本机 sans 同字形，0 差异正常）
        _fam = np.asarray(svg_edit.render_svg(
            svg_src, dpi=200,
            edit={"font_family": "Times New Roman"})["img"])
        check("SVG font family change changes the picture",
              int((_fam != a0).sum()) > 2000,
              "diff px %d" % int((_fam != a0).sum()))

        # 图例缩放 + 图例移动（图例操作的多图例语义已由 tracks 面板覆盖）
        _ls = np.asarray(svg_edit.render_svg(
            svg_src, dpi=200,
            edit={"legend": {"selector": "#legend_1", "action": "scale",
                             "font_scale": 1.6}})["img"])
        check("SVG legend font scale changes the picture",
              int((_ls != a0).sum()) > 500,
              "diff px %d" % int((_ls != a0).sum()))
        _lm = np.asarray(svg_edit.render_svg(
            svg_src, dpi=200,
            edit={"legend": {"selector": "#legend_1", "action": "move",
                             "dx": 30, "dy": 12}})["img"])
        check("SVG legend move changes the picture",
              int((_lm != a0).sum()) > 500,
              "diff px %d" % int((_lm != a0).sum()))

        # 双图例面板（tracks）：两个图例都要被同一动作命中
        _tr = HERE / "outputs/_svgdemo/panel_tracks.svg"
        if _tr.is_file():
            _t_base = np.asarray(svg_edit.render_svg(str(_tr), dpi=200)["img"])
            _t_out, _t_notes = svg_edit.apply_svg_edit(
                str(_tr), {"legend": {"action": "hide"}})
            _n_hidden = sum(1 for n in _t_notes if "legend removed" in n)
            check("SVG multi-legend panel hides all legends", _n_hidden >= 2,
                  "hidden %d / notes %s" % (_n_hidden, _t_notes))
            if _t_out:
                _t_after = np.asarray(
                    svg_edit.rasterize_svg(_t_out, dpi=200))
                check("SVG multi-legend hide changes the picture",
                      int((_t_after != _t_base).sum()) > 2000,
                      "diff px %d" % int((_t_after != _t_base).sum()))

        # ── 预览端点：每种 SVG 都必须出图 ────────────────────────────
        #   ⚠ 早年 render_preview 传 dpi=max_px：12.8 in 宽的 SVG 在
        #     1600 dpi 下要 20k×18k px ≈ 1 GB → 静默异常 → /api/img 404
        #     （面板库卡片裂图，而 selftest 只测了 net/tracks 两个小 SVG）。
        for _f in ("panel_net.svg", "panel_tracks.svg", "panel_pca.svg",
                   "panel_map.svg"):
            _u = ("/api/img?src=" + urllib.parse.quote(
                str(HERE / "outputs" / "_svgdemo" / _f)) + "&max=1600")
            try:
                with urllib.request.urlopen(BASE + _u, timeout=300) as _r:
                    _n = len(_r.read())
                check("SVG 预览端点 200 " + _f, _n > 5000, "%d bytes" % _n)
            except Exception as _e:                               # noqa: BLE001
                check("SVG 预览端点 200 " + _f, False, repr(_e)[:80])

        # ── SVG + PNG 混排：几何按**资产比**反解（零形变）────────────────
        #   ⚠ 早期硬写 w/h，SVG 面板形变 5.3%、PNG 面板形变 91.6%，
        #     而 matplotlib/fitz 对形变零报错 → 形变只能靠几何断言兜住。
        a_svg = ins["size_in"][1] / ins["size_in"][0]        # aspect = h/w
        png_src = srcs[ROWS[0][0]]["src"]
        with Image.open(png_src) as _im:
            a_png = _im.size[1] / float(_im.size[0])
        MG, GP, TOPM, CAP = 0.08, 0.10, 0.90, 0.24
        avail = 7.5 - 2 * MG
        h_mix = (avail - GP) / (1.0 / a_svg + 1.0 / a_png)
        w_svg, w_png = h_mix / a_svg, h_mix / a_png
        y_mix = TOPM + h_mix + GP + CAP
        v_layout = {
            "canvas": {"width_in": 7.5,
                       "height_in": round(y_mix - GP + MG, 4), "dpi": DPI},
            "style": {"bg": "white", "frame": True, "caption_fs": 7,
                      "caption_gap_in": 0.06, "margin_in": MG,
                      "gap_in": GP, "pack_margin_in": MG, "letter_fs": 11},
            "title": {"text": "", "subtitle": ""},
            "items": [
                {"id": "s1", "src": svg_src, "x_in": MG, "y_in": TOPM,
                 "w_in": round(w_svg, 4), "h_in": round(h_mix, 4),
                 "letter": "A", "caption": "svg panel", "z": 0,
                 "adjust": {"svg": {"colors": [
                     {"from": "#0f4d92", "to": "#0072B2"}]}}},
                {"id": "s2", "src": png_src, "x_in": round(MG + w_svg + GP, 4),
                 "y_in": TOPM, "w_in": round(w_png, 4),
                 "h_in": round(h_mix, 4), "letter": "B", "caption": "",
                 "z": 0},
            ],
        }
        res3 = post("/api/render", {"layout": v_layout, "stem": "_selftest_svg",
                                    "formats": ["png", "pdf"], "dpi": DPI})
        if res3.get("error"):
            check("SVG + PNG mixed composite", False, res3["error"])
        else:
            check("SVG + PNG mixed composite", True)
            check("mixed composite reports vector embed",
                  any(("vector embed" in w) or ("SVG vector" in w)
                      for w in res3.get("warnings", [])),
                  str(res3.get("warnings")))
            check("mixed composite has no distortion warning",
                  not any("distort" in w.lower() or "形变" in w
                          for w in res3.get("warnings", [])),
                  str(res3.get("warnings")))
            pdf3 = res3.get("pdf")
            if pdf3 and Path(pdf3).is_file():
                import fitz

                _d = fitz.open(pdf3)
                _pg = _d[_d.page_count - 1]
                # 布局里本就有 1 个 PNG 面板 → 必然 1 张位图；
                # 关键是 SVG 面板**不是**位图：矢量路径 + 可提取文本。
                check("mixed PDF has exactly the PNG panel as raster",
                      len(_pg.get_images(full=True)) == 1,
                      "images=%d" % len(_pg.get_images(full=True)))
                check("mixed PDF keeps SVG as vector (paths)",
                      len(_pg.get_drawings()) > 100,
                      "paths=%d" % len(_pg.get_drawings()))
                check("mixed PDF text extractable (SVG text -> PDF text)",
                      len(_pg.get_text().strip()) > 50,
                      "text len=%d" % len(_pg.get_text().strip()))
                _d.close()

    print("\n== AI 代理（SVG Chat）==")
    try:
        _ai_proxy_stream_test()
    except Exception as e:                                   # noqa: BLE001
        check("AI proxy stream test", False, "%s: %s" % (type(e).__name__, e))

    print("\n== 结果 ==")
    print("  全部通过 ✔" if ok_all else "  有失败项 ✘")
    return 0 if ok_all else 1


if __name__ == "__main__":
    raise SystemExit(main())
