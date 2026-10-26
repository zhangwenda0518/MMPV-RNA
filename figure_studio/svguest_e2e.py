# -*- coding: utf-8 -*-
"""svguest 深度整合端到端测试。

覆盖：工作台 → 会话 → svguest 多标签 → 编辑 → 写回 → 备份 → 回滚 → 导出直通。

⚠ 有一条**回归守卫**必须留着：浮条的 computed style。
   曾经 <style> 里混进一行 Markdown 标题，浏览器把 `#` 当 ID 选择器，
   整条 `#fsBar{…}` 被静默丢弃 —— 浮条失去 position:fixed、挤进状态栏，
   全程零报错。所以这里直接断言 getComputedStyle 的结果。
"""
import atexit
import shutil
from pathlib import Path

from playwright.sync_api import sync_playwright

HERE = Path(__file__).resolve().parent
DEMO = HERE / "outputs" / "_svgdemo"
TGT = HERE / "outputs" / "_svgedit_test"
BACKUP = HERE / "outputs" / "_svguest_backup"

# 干净的测试目录
if TGT.exists():
    shutil.rmtree(TGT)
TGT.mkdir(parents=True, exist_ok=True)
A = TGT / "a_net.svg"
B = TGT / "b_map.svg"
shutil.copy2(DEMO / "panel_net.svg", A)
shutil.copy2(DEMO / "panel_map.svg", B)
A_S, B_S = str(A.resolve()), str(B.resolve())
a_before = A.read_text(encoding="utf-8")
b_before = B.read_text(encoding="utf-8")
a_rects0 = a_before.count("<rect")

# 清掉本测试可能产生的旧备份，避免「备份数」断言被历史干扰
if BACKUP.is_dir():
    for f in BACKUP.glob("a_net.svg.bak_*"):
        f.unlink(missing_ok=True)

# 快照目录也必须清空：否则 shot_svguest.py 预置的旧快照会混进来，
# 而 rglob 的顺序不确定，断言可能拿到旧的那份（踩过）。
AUTOSAVE = HERE / "outputs" / "_svguest_autosave"
shutil.rmtree(AUTOSAVE, ignore_errors=True)

errs, fails, results = [], [], []


def check(name, ok, detail=""):
    results.append((name, bool(ok), str(detail)[:110]))


class SgFrame:
    """把 svguest 所在的**同源 iframe** 包装成够用的 Page 形状。

    平台已把精修从 window.open 改成内嵌浮层（#sgLayer > iframe#sgFrame），
    所以 ctx.expect_page() 那条路没了。Frame 与 Page 只差几个方法：
    locator / evaluate / url 直接转发；wait_for_timeout / keyboard / mouse
    借主 page（键盘鼠标本来就是全局输入）；close 改成点关闭按钮。

    ⚠ 关键前提：`locator.bounding_box()` 返回的是**主框架视口**坐标
      （Playwright 的约定，iframe 里的元素也一样），所以 page.mouse 可以
      直接用，不需要自己加 iframe 偏移 —— 加了反而全歪。
    """

    def __init__(self, page, sel="#sgFrame", name="svguest"):
        self._p, self._sel, self._name = page, sel, name

    @property
    def frame(self):
        f = self._p.frame(name=self._name)
        if f is None:
            raise AssertionError("iframe #sgFrame 还没挂载")
        return f

    def locator(self, sel):
        return self.frame.locator(sel)

    def evaluate(self, js, arg=None):
        return self.frame.evaluate(js) if arg is None else self.frame.evaluate(js, arg)

    def wait_for_timeout(self, ms):
        self._p.wait_for_timeout(ms)

    def wait_ready(self, timeout=60000):
        """等桥接浮条出现 = 会话取到、面板载入完成。"""
        self.frame.wait_for_selector("#fsBar", timeout=timeout)

    def focus(self):
        """把键盘焦点挪进 iframe，返回平台侧 activeElement 的 id。

        ⚠⚠ 内嵌之后这一步是**必须的**：以前是独立窗口，窗口自带焦点；
        现在点的是平台上的按钮，焦点留在平台上，Ctrl+S 和 m/v 这些工具
        快捷键一个都到不了 svguest —— 表现为「画不出东西、保存没反应」，
        而且**零报错**。所以先 window.focus()，不行就真点一下浮条。
        """
        self.frame.evaluate("() => { try { window.focus(); } catch (e) {} }")
        name = self._sel.lstrip("#")
        if self._p.evaluate("() => (document.activeElement || {}).id || ''") != name:
            self.locator("#fsName").click()
        return self._p.evaluate("() => (document.activeElement || {}).id || ''")

    @property
    def url(self):
        return self.frame.url

    @property
    def mouse(self):
        return self._p.mouse

    @property
    def keyboard(self):
        return self._p.keyboard

    def close(self):
        self._p.locator("#btnSgClose").click()
        self._p.wait_for_timeout(300)


with sync_playwright() as pw:
    br = pw.chromium.launch(channel="chrome", headless=True)
    ctx = br.new_context(viewport={"width": 1600, "height": 950})
    page = ctx.new_page()

    # ⚠ 7e 会**故意**打两条失败请求（未配 key → 401、白名单外主机 → 403）。
    #   以前精修在独立窗口里，那些响应落在另一个 Page 上、这里监听不到；
    #   改成内嵌 iframe 之后它们也归 page 管了，所以必须**精确**放行这两条，
    #   否则「errs/fails 为空」这条总闸会永远红 —— 而它正是回归的最后一道网。
    def _expected_bad(r):
        return "/api/ai/chat" in r.url and r.status in (401, 403)

    def _on_console(m):
        if m.type != "error":
            return
        # 网络失败的 console 噪声不带 URL，而 fails 里已经有带 URL 的精确记录，
        # 这里丢掉以免同一条错误既算 fails 又算 errs（没法精确过滤）。
        if m.text.startswith("Failed to load resource"):
            return
        errs.append(m.type + ": " + m.text)

    page.on("console", _on_console)
    page.on("response", lambda r: fails.append("HTTP %d %s" % (r.status, r.url[-70:]))
            if r.status >= 400 and not _expected_bad(r) else None)
    page.on("dialog", lambda d: d.accept())      # 回滚的 confirm()

    page.goto("http://127.0.0.1:8770", wait_until="load")
    page.wait_for_timeout(400)
    page.fill("#rootInput", str(TGT))
    page.locator("#btnScan").click()
    page.wait_for_selector(".lib-card", timeout=60000)
    page.wait_for_timeout(500)

    # ── 1. 选中 a_net.svg ────────────────────────────────────────
    cards = page.locator(".lib-card")
    idx_a = 0
    for i in range(cards.count()):
        if "a_net" in cards.nth(i).inner_text():
            idx_a = i
            break
    cards.nth(idx_a).click()
    page.wait_for_timeout(1000)

    check("精修区块存在", page.locator("#sgBlock").count() == 1)
    check("精修按钮对 SVG 可用", page.locator("#btnSgOne").is_enabled())
    check("全部 SVG 面板按钮可用", page.locator("#btnSgAll").is_enabled(),
          page.locator("#btnSgAll").inner_text())

    # ── 2. 打开 svguest（单面板，内嵌浮层） ─────────────────────
    check("精修浮层初始是隐藏的", not page.locator("#sgLayer.on").count())
    page.locator("#btnSgOne").click()
    sg = SgFrame(page)
    sg.wait_ready()
    sg.wait_for_timeout(2800)

    check("精修浮层弹出（内嵌 iframe，不是新窗口）",
          page.locator("#sgLayer.on").count() == 1)
    check("没有开新窗口（页面数仍是 1）", len(ctx.pages) == 1, "%d 个" % len(ctx.pages))
    check("iframe 指向 svguest", "svguest" in sg.url, sg.url[-40:])
    check("hash 带 session token", "#s=" in sg.url, sg.url[-40:])
    check("浮层标题栏写明文件名", "a_net.svg" in (page.locator("#sgLayerSub").inner_text() or ""),
          page.locator("#sgLayerSub").inner_text()[:50])
    check("桥接浮条出现", sg.locator("#fsBar").count() == 1)
    check("键盘焦点自动落进 iframe（否则 Ctrl+S / 工具键全失效）",
          sg.focus() == "sgFrame")

    # ⚠ 自查：应用搬进 iframe 之后，page.on('console') 还收不收得到里面的报错？
    #   收不到就等于错误监听静默失明 —— 后面「errs 为空」的断言会变成假阴性。
    #   所以主动打一条进去验证，验证完再摘掉，免得污染最终判定。
    sg.frame.evaluate("() => console.error('fs-iframe-console-probe')")
    sg.wait_for_timeout(400)
    probe_seen = any("fs-iframe-console-probe" in e for e in errs)
    errs[:] = [e for e in errs if "fs-iframe-console-probe" not in e]
    check("iframe 内的 console.error 能被收集（错误监听没失明）", probe_seen)

    # ⚠⚠ 回归守卫：浮条必须真的浮起来
    cs = sg.evaluate("""() => {
        const b = document.getElementById('fsBar');
        if (!b) return null;
        const s = getComputedStyle(b);
        const r = b.getBoundingClientRect();
        return {pos: s.position, bottom: s.bottom, z: s.zIndex, bg: s.backgroundColor,
                top: Math.round(r.top), h: Math.round(r.height), w: Math.round(r.width)};
    }""")
    check("浮条 CSS 真的生效（position:fixed）",
          cs and cs["pos"] == "fixed", cs and cs["pos"])
    check("浮条浮在画布上方（不挤进状态栏）",
          cs and cs["top"] > 100 and cs["h"] > 20,
          cs and "top=%s h=%s" % (cs["top"], cs["h"]))
    check("浮条有背景色", cs and cs["bg"] not in ("", "rgba(0, 0, 0, 0)"),
          cs and cs["bg"])
    check("浮条按钮齐全",
          sg.locator("#fsBar button").count() >= 6,
          "%d 个按钮" % sg.locator("#fsBar button").count())

    check("svguest 画布有内容", sg.locator("#mStage svg *").count() > 50,
          "%d 个节点" % sg.locator("#mStage svg *").count())
    check("浮条显示文件名", "a_net.svg" in (sg.locator("#fsName").inner_text() or ""),
          sg.locator("#fsName").inner_text()[:40])

    # ── 3. 真实编辑：矩形工具 ───────────────────────────────────
    box = sg.locator("#mStage").bounding_box() or {"x": 100, "y": 100,
                                                   "width": 900, "height": 600}
    cx, cy = box["x"] + box["width"] * 0.75, box["y"] + box["height"] * 0.8
    sg.keyboard.press("m")
    sg.wait_for_timeout(250)
    sg.mouse.move(cx, cy)
    sg.mouse.down()
    sg.mouse.move(cx + 120, cy - 70, steps=8)
    sg.mouse.up()
    sg.wait_for_timeout(500)
    n_rects = sg.evaluate("() => document.querySelectorAll('#mStage svg rect').length")
    check("矩形工具画出了形状", n_rects > 0, "%d 个 rect" % n_rects)

    # ── 4. Ctrl+S 写回 ─────────────────────────────────────────
    sg.keyboard.press("Control+s")
    sg.wait_for_timeout(2400)
    toast = sg.locator("#fsToast").inner_text() or ""
    check("保存后出现成功提示", "已保存回工作台" in toast, toast[:60])

    a_after = A.read_text(encoding="utf-8")
    check("磁盘文件被改写", a_after != a_before,
          "%d → %d 字符" % (len(a_before), len(a_after)))
    check("新方块写进 SVG", a_after.count("<rect") > a_rects0,
          "rect %d → %d" % (a_rects0, a_after.count("<rect")))
    baks = sorted(BACKUP.glob("a_net.svg.bak_*"))
    check("旧版已备份", len(baks) >= 1, "%d 份" % len(baks))

    sg.wait_for_timeout(1200)
    st = page.locator("#statusText").inner_text()
    check("工作台收到写回通知", "svguest" in st, st[:50])

    # ── 5. 工作台侧备份历史 + 回滚 ─────────────────────────────
    # ⚠ 精修浮层是**全屏**的，会挡住工作台 —— 要操作检查器必须先收起来。
    #   这不是妥协，正是真实用法：改完关掉，再在工作台里看历史/回滚。
    #   （旧的弹窗方案里父窗口一直可点，所以这步以前不需要。）
    sg.close()
    page.wait_for_timeout(900)
    rows = page.locator("#sgHistory .sg-row")
    check("检查器列出备份历史", rows.count() >= 1, "%d 行" % rows.count())
    check("历史行有回滚按钮",
          page.locator("#sgHistory .sg-restore").count() >= 1)

    page.locator("#sgHistory .sg-restore").first.click()
    page.wait_for_timeout(2200)
    a_restored = A.read_text(encoding="utf-8")
    check("回滚把文件还原了", a_restored == a_before,
          "%d → %d 字符" % (len(a_after), len(a_restored)))
    check("回滚也有状态提示", "回滚" in (page.locator("#statusText").inner_text() or ""),
          page.locator("#statusText").inner_text()[:50])

    # ── 6. 导出直通（PNG 300 dpi 落到原图同目录） ───────────────
    # 回滚完再开一次（每次都是新会话，从磁盘读回当前内容）
    page.locator("#btnSgOne").click()
    sg = SgFrame(page)
    sg.wait_ready()
    sg.focus()
    sg.wait_for_timeout(1500)
    sg.keyboard.press("Control+s")          # 再存一次，让 a 有内容可导
    sg.wait_for_timeout(2000)
    sg.locator("#fsPng300").click()
    sg.wait_for_timeout(6000)
    pngs = sorted(TGT.glob("a_net_edit_300dpi*.png"))
    check("PNG 导出落到原图同目录", len(pngs) == 1,
          pngs[0].name if pngs else "没找到")
    check("导出的 PNG 不是空文件",
          bool(pngs) and pngs[0].stat().st_size > 20000,
          "%d 字节" % (pngs[0].stat().st_size if pngs else 0))
    etoast = sg.locator("#fsToast").inner_text() or ""
    check("导出有提示", "已导出" in etoast, etoast[:70])

    # ── 7. SVG 副本导出 ────────────────────────────────────────
    sg.locator("#fsSvg").click()
    sg.wait_for_timeout(2500)
    svgs = sorted(TGT.glob("a_net_edit.svg"))
    check("SVG 副本导出成功", len(svgs) == 1,
          svgs[0].name if svgs else "没找到")

    # ── 7b. JPG 300 导出（走 dl 分流器，按 MIME 认出来） ────────
    sg.locator("#fsJpg300").click()
    sg.wait_for_timeout(7000)
    jpgs = sorted(TGT.glob("a_net_edit_300dpi*.jpg"))
    check("JPG 导出落到原图同目录", len(jpgs) == 1,
          jpgs[0].name if jpgs else "没找到")
    check("导出的 JPG 不是空文件",
          bool(jpgs) and jpgs[0].stat().st_size > 10000,
          "%d 字节" % (jpgs[0].stat().st_size if jpgs else 0))

    # ── 7c. 矢量 PDF：Ctrl+P 不再弹打印对话框 ───────────────────
    # 判据不是「文件存在」，而是「文字可提取」—— 那是矢量的定义。
    sg.keyboard.press("Control+p")
    sg.wait_for_timeout(9000)
    pdfs = sorted(TGT.glob("a_net_edit*.pdf"))
    check("Ctrl+P 直接落盘 PDF（无打印对话框）", len(pdfs) == 1,
          pdfs[0].name if pdfs else "没找到")
    if pdfs:
        import fitz
        d = fitz.open(str(pdfs[0]))
        ptxt = d[0].get_text().strip()
        pw_pt = d[0].rect.width
        d.close()
        check("PDF 是矢量（文字可提取）", len(ptxt) > 20, "%d 字符" % len(ptxt))
        check("PDF 页宽与 SVG 声明一致（无 px→pt 虚胖）",
              abs(pw_pt - 509.76) < 1.0, "%.2f pt（期望 509.76）" % pw_pt)

    # ── 7d. A 表接管：菜单 / 快捷键 / 命令面板是同一个入口 ──────
    hooked = sg.evaluate("""() => {
      try {
        if (typeof A === 'undefined' || !A) return 'no-A';
        return (A.__fsHooked ? 'hooked' : 'raw') + '|' + typeof A.svg + '|' + typeof A.pdf;
      } catch (e) { return 'err:' + e.message; }
    }""")
    check("A 表已接管（菜单/快捷键/命令面板同一入口）",
          str(hooked).startswith("hooked"), str(hooked))

    # File ▸ Export SVG 走的就是 A.svg —— 这里直接调它，等价于点菜单
    a_pre_menu = A.read_text(encoding="utf-8")
    bx = sg.locator("#mStage").bounding_box() or {"x": 100, "y": 100,
                                                 "width": 900, "height": 600}
    sg.keyboard.press("m")
    sg.wait_for_timeout(200)
    sg.mouse.move(bx["x"] + bx["width"] * 0.32, bx["y"] + bx["height"] * 0.32)
    sg.mouse.down()
    sg.mouse.move(bx["x"] + bx["width"] * 0.32 + 90,
                  bx["y"] + bx["height"] * 0.32 + 50, steps=6)
    sg.mouse.up()
    sg.wait_for_timeout(600)
    sg.evaluate("() => A.svg()")
    sg.wait_for_timeout(4000)
    a_post_menu = A.read_text(encoding="utf-8")
    check("A.svg() 写回原文件（等价 File▸Export SVG，不弹文件选择器）",
          a_post_menu != a_pre_menu,
          "%d → %d 字符" % (len(a_pre_menu), len(a_post_menu)))

    # ── 7e. AI 绘图：apiPost 被平台代理接管 ────────────────────
    # 浏览器直连 api.deepseek.com 必然 CORS 失败（上游靠自家 Cloudflare 中继
    # 绕过，而那个 Worker 不放行 127.0.0.1），所以必须走本机服务端转发。
    aip = sg.evaluate("""() => {
      try { return typeof apiPost !== 'function' ? 'missing'
                  : (apiPost.__fsProxied ? 'proxied' : 'raw'); }
      catch (e) { return 'err:' + e.message; }
    }""")
    check("apiPost 已改走平台代理", aip == "proxied", str(aip))

    # 真打一发：没配 key 时应由服务端回 401，且错误体是
    # {"error":{"message":...}} —— 上游 apiMsg() 解析的就是 error.message，
    # 给它字符串会把整段 JSON 糊到用户脸上。
    air = sg.evaluate("""async () => {
      try {
        const r = await apiPost('Ds', 'https://api.deepseek.com/chat/completions', '',
          { model: 'deepseek-chat', messages: [{ role: 'user', content: 'hi' }], max_tokens: 4 });
        return r.status + '|' + (await r.text()).slice(0, 200);
      } catch (e) { return 'throw:' + e.message; }
    }""")
    check("代理链路可用（未配 key → 401 且文案可读）",
          str(air).startswith("401") and "API key" in str(air), str(air)[:120])

    # 防 SSRF：白名单外的主机必须被服务端挡掉
    evil = sg.evaluate("""async () => {
      try {
        const r = await apiPost('Ds', 'https://evil.example.com/x', '', {});
        return r.status + '|' + (await r.text()).slice(0, 90);
      } catch (e) { return 'throw:' + e.message; }
    }""")
    check("白名单外主机被拒（403）", str(evil).startswith("403"), str(evil)[:110])

    # 上游「Models and keys」的测试连接也走 apiPost，顺便验证它的错误呈现。
    # ⚠ probe() 在失败时是**抛错**（不是返回数组），错误文案由 apiMsg() 拼；
    #   我们给的是 {"error":{"message","code"}}，所以会渲染成 "[no_key] 还没有…"。
    #   若错误体写成字符串，这里就会变成整段裸 JSON。
    probe_res = sg.evaluate("""async () => {
      try { return 'ok:' + JSON.stringify(await probe('Ds')).slice(0, 200); }
      catch (e) { return 'throw:' + String(e.message || e); }
    }""")
    check("上游 probe() 的错误信息可读（[code] 文案，不是裸 JSON）",
          "API key" in str(probe_res) and '{"error"' not in str(probe_res)
          and "no_key" in str(probe_res),
          str(probe_res)[:120])

    sg.close()

    # ── 7g. 浮层退出：关闭按钮；以及「Esc 不许被平台抢走」的守卫 ──
    check("关闭按钮把浮层收起", not page.locator("#sgLayer.on").count())
    check("收起后 iframe 被卸载回 about:blank",
          "about:blank" in (page.locator("#sgFrame").get_attribute("src") or ""),
          (page.locator("#sgFrame").get_attribute("src") or "")[:30])

    page.locator("#btnSgOne").click()
    sg2 = SgFrame(page)
    sg2.wait_ready()
    sg2.wait_for_timeout(1500)
    check("重新打开是干净的新会话（不是复用旧 iframe）",
          page.locator("#sgLayer.on").count() == 1
          and "about:blank" not in (page.locator("#sgFrame").get_attribute("src") or ""))
    check("键盘焦点自动落进 iframe（否则 Ctrl+S / 工具键全失效）",
          sg2.focus() == "sgFrame")
    # ⚠⚠ 回归守卫：焦点在 iframe 里时 Esc 属于上游（取消工具/取消选择），
    #    平台**不许**抢来关浮层。曾这么做过 —— 探针测出上游把每一次
    #    Escape 都 preventDefault，照此关闭会导致「想取消选择 → 编辑器没了」。
    sg2.keyboard.press("Escape")
    sg2.wait_for_timeout(800)
    check("iframe 里的 Esc 不会被平台抢去关浮层（仍归上游）",
          page.locator("#sgLayer.on").count() == 1)
    # 焦点回到平台时 Esc 才关（这时没有歧义）
    page.locator("#sgLayerSub").click()
    page.wait_for_timeout(200)
    check("点顶栏后焦点回到平台",
          page.evaluate("() => (document.activeElement||{}).id||''") != "sgFrame")
    page.keyboard.press("Escape")
    page.wait_for_timeout(600)
    check("焦点在平台上时 Esc 能收起浮层",
          not page.locator("#sgLayer.on").count())

    # ── 7f. 磁盘自动快照（Auto-save 的服务端那一层）─────────────
    # 上游自带 Auto-save 存浏览器 IndexedDB（默认开），那是崩溃安全网但换机器
    # 就没了；这里再验证服务端落盘的那层。autoMs 调到 1.5 s 以便真跑节流链路。
    import json as _json
    import urllib.request as _ur

    _req = _ur.Request(
        "http://127.0.0.1:8770/api/svg/session",
        data=_json.dumps({"items": [{"src": A_S, "name": "a_net.svg"}]}).encode(),
        headers={"Content-Type": "application/json"})
    _tok = _json.loads(_ur.urlopen(_req).read())["token"]
    asg = ctx.new_page()
    asg.on("console", lambda m: errs.append("auto " + m.type + ": " + m.text) if m.type == "error" else None)
    asg.goto("http://127.0.0.1:8770/static/vendor/svguest.html#s=%s&autoMs=1500" % _tok)
    asg.wait_for_timeout(4500)

    check("上游 IndexedDB Auto-save 默认开启（opts.auto）",
          asg.evaluate("() => (typeof opts !== 'undefined' ? !!opts.auto : 'missing')") is True)

    st0 = asg.evaluate("() => window.__fsAutoState ? window.__fsAutoState() : null")
    check("桥接快照节流器已启动（autoMs 覆盖生效）",
          isinstance(st0, dict) and st0.get("ms") == 1500, str(st0))

    # 画一个矩形 → dv 自增、dirty 变真
    abx = asg.locator("#mStage").bounding_box() or {"x": 100, "y": 100,
                                                    "width": 900, "height": 600}
    asg.keyboard.press("m")
    asg.wait_for_timeout(200)
    asg.mouse.move(abx["x"] + abx["width"] * 0.5, abx["y"] + abx["height"] * 0.5)
    asg.mouse.down()
    asg.mouse.move(abx["x"] + abx["width"] * 0.5 + 100,
                   abx["y"] + abx["height"] * 0.5 + 60, steps=6)
    asg.mouse.up()
    asg.wait_for_timeout(3000)          # 轮询 2 s + 节流 1.5 s

    st1 = asg.evaluate("() => window.__fsAutoState ? window.__fsAutoState() : null")
    check("节流器捕捉到编辑（lastDv 追平 dv 且 dirty）",
          isinstance(st1, dict) and st1.get("lastDv") == st1.get("dv")
          and st1.get("dirty") is True, str(st1))

    asg.wait_for_timeout(2500)
    # 按 mtime 取最新那份（rglob 顺序不确定，直接 [0] 会拿到旧的）
    snaps = sorted(AUTOSAVE.rglob("*.svg"),
                   key=lambda p: p.stat().st_mtime, reverse=True) if AUTOSAVE.is_dir() else []
    check("磁盘快照真的落盘（不覆盖原文件）",
          len(snaps) >= 1, "%d 份：%s" % (len(snaps),
                                          snaps[0].name if snaps else "无"))
    if snaps:
        # 快照里应该有刚画的矩形，且原文件此刻**还没**被改（写回只由 Ctrl+S 触发）
        snap_txt = snaps[0].read_text(encoding="utf-8")
        check("快照内容含本次编辑，原文件仍未被自动覆盖",
              snap_txt.count("<rect") > a_post_menu.count("<rect")
              and A.read_text(encoding="utf-8") == a_post_menu,
              "snap rect=%d / 原文件 rect=%d" % (snap_txt.count("<rect"),
                                                 A.read_text(encoding="utf-8").count("<rect")))

    auto_lbl = asg.locator("#fsAuto").inner_text() or ""
    check("浮条显示自动快照状态", "自动快照" in auto_lbl, auto_lbl[:40])

    # 上游那层（IndexedDB sgDocs）也应当在写
    idb = asg.evaluate("""() => new Promise((ok) => {
      try {
        const r = indexedDB.open('sgDocs', 1);
        r.onupgradeneeded = () => { try { r.result.createObjectStore('d', { keyPath: 'id' }); } catch (e) {} };
        r.onsuccess = () => {
          const db = r.result;
          if (!db.objectStoreNames.contains('d')) return ok('no-store');
          const q = db.transaction('d', 'readonly').objectStore('d').count();
          q.onsuccess = () => ok('n=' + q.result);
          q.onerror = () => ok('count-err');
        };
        r.onerror = () => ok('open-err');
      } catch (e) { ok('throw:' + e.message); }
    })""")
    check("上游 IndexedDB 副本也在写（sgDocs 有条目）",
          str(idb).startswith("n=") and str(idb) != "n=0", str(idb))
    asg.close()

    # ── 8. 多标签：全部 SVG 面板 ───────────────────────────────
    # 先把第二个 SVG 也加进画布，这样「全部 SVG 面板」才有 2 个
    page.wait_for_timeout(400)
    for i in range(cards.count()):
        if "b_map" in cards.nth(i).inner_text():
            cards.nth(i).click()
            break
    page.wait_for_timeout(900)
    n_svg = page.evaluate("() => document.querySelectorAll('.item').length")
    check("画布上有 2 个面板", n_svg == 2, "%d 个" % n_svg)
    check("全部按钮显示数量", "2" in (page.locator("#btnSgAll").inner_text() or ""),
          page.locator("#btnSgAll").inner_text())

    page.locator("#btnSgAll").click()
    mg = SgFrame(page)
    mg.wait_ready()
    mg.wait_for_timeout(7000)               # 两张图，等两个 readFile
    mg.focus()
    n_tabs = mg.evaluate("() => (typeof TABS !== 'undefined' ? TABS.length : -1)")
    check("多标签真的开了 2 个", n_tabs == 2, "TABS=%s" % n_tabs)
    check("浮条显示标签数",
          "标签" in (mg.locator("#fsN").inner_text() or ""),
          mg.locator("#fsN").inner_text()[:30])
    check("多标签时按钮变成「全部保存」",
          "全部保存" in (mg.locator("#fsBack").inner_text() or ""),
          mg.locator("#fsBack").inner_text()[:30])
    check("第二个面板也载入了",
          mg.locator("#mStage svg *").count() > 50,
          "%d 个节点" % mg.locator("#mStage svg *").count())

    # 在第二个标签上画个方块，然后全部保存
    box2 = mg.locator("#mStage").bounding_box() or {"x": 100, "y": 100,
                                                    "width": 900, "height": 600}
    mx, my = box2["x"] + box2["width"] * 0.6, box2["y"] + box2["height"] * 0.7
    mg.keyboard.press("m")
    mg.wait_for_timeout(250)
    mg.mouse.move(mx, my)
    mg.mouse.down()
    mg.mouse.move(mx + 110, my - 60, steps=8)
    mg.mouse.up()
    mg.wait_for_timeout(600)

    mg.keyboard.press("Control+s")
    mg.wait_for_timeout(6000)
    mtoast = mg.locator("#fsToast").inner_text() or ""
    check("全部保存有结果提示", "已保存" in mtoast, mtoast[:70])
    b_after = B.read_text(encoding="utf-8")
    check("第二个面板被写回", b_after != b_before,
          "%d → %d 字符" % (len(b_before), len(b_after)))
    check("批量保存两份备份都在",
          len(list(BACKUP.glob("a_net.svg.bak_*"))) >= 2
          and len(list(BACKUP.glob("b_map.svg.bak_*"))) >= 1,
          "a=%d b=%d" % (len(list(BACKUP.glob("a_net.svg.bak_*"))),
                         len(list(BACKUP.glob("b_map.svg.bak_*")))))
    mg.close()

    # ── 9. 空布局 + PNG 面板：精修按钮禁用 ─────────────────────
    png_dir = HERE.parent / "archive" / "20260917_phylo_cleanup" / "_absorb_test_20260916" / "show"
    if png_dir.is_dir():
        page.reload(wait_until="load")        # 重新加载 = 空画布
        page.wait_for_timeout(600)
        page.fill("#rootInput", str(png_dir))
        page.locator("#btnScan").click()
        page.wait_for_selector(".lib-card", timeout=60000)
        page.wait_for_timeout(500)
        check("空布局时全部精修禁用", not page.locator("#btnSgAll").is_enabled())
        page.locator(".lib-card").first.click()
        page.wait_for_timeout(900)
        check("PNG 面板时精修按钮禁用", not page.locator("#btnSgOne").is_enabled())
        check("PNG 面板时全部精修仍禁用", not page.locator("#btnSgAll").is_enabled())
        page.screenshot(path=str(HERE / "outputs" / "_ui" / "18_sg_block.png"))

    # ── 10. 平台侧：AI 绘图 + 自动快照两个新块 ─────────────────
    page.fill("#rootInput", str(TGT))
    page.locator("#btnScan").click()
    page.wait_for_selector(".lib-card", timeout=60000)
    page.wait_for_timeout(600)
    for i in range(cards.count()):
        if "a_net" in cards.nth(i).inner_text():
            cards.nth(i).click()
            break
    page.wait_for_timeout(1500)

    check("检查器出现「AI 绘图」区块", page.locator("#aiBlock").count() == 1)
    ai_tag = (page.locator("#aiTag").inner_text() or "").strip()
    check("AI 状态标签有内容", ai_tag != "", ai_tag)
    check("AI 两个密钥输入框都在",
          page.locator("#aiKeyDs").count() == 1 and page.locator("#aiKeyGlm").count() == 1)
    ai_hint = page.locator("#aiHint").inner_text() or ""
    check("AI 提示里写明密钥存放位置（.ai_keys.json）",
          "ai_keys.json" in ai_hint, ai_hint[:90])

    # 7f 刚给 a_net 存过快照，这里应能列出来
    as_txt = page.locator("#asList").inner_text() or ""
    check("自动快照列表已列出（含时间 + 恢复按钮）",
          "KB" in as_txt and "恢复" in as_txt, as_txt[:70])
    as_tag = (page.locator("#asTag").inner_text() or "").strip()
    check("自动快照标签显示份数", "份" in as_tag, as_tag)
    page.screenshot(path=str(HERE / "outputs" / "_ui" / "20_ai_autosave.png"))
    br.close()

# 测试产物清掉（快照目录是本次测试造出来的）
shutil.rmtree(HERE / "outputs" / "_svguest_autosave", ignore_errors=True)

_reported = False


def report():
    """打印结果并返回是否全绿。

    ⚠ 挂在 atexit 上：中途某个 locator 超时抛异常时，已经跑过的断言不至于
      连同一句诊断一起丢掉（踩过 —— 只看到 traceback，不知道前面哪步先坏）。
    """
    global _reported
    if _reported:
        return True
    _reported = True
    print("== svguest 深度整合端到端 ==")
    for n, ok, d in results:
        print("  [%s] %-42s %s" % ("PASS" if ok else "FAIL", n, d))
    bad = [n for n, ok, _ in results if not ok]
    print("console errors:", errs[:4])
    print("4xx/5xx:", fails[:4])
    ok_all = not bad and not errs and not fails
    print("结果:", "全部通过 ✔ (%d 项)" % len(results) if ok_all
          else "有失败: %s" % bad)
    return ok_all


atexit.register(report)
# 中途抛异常时先把结果打出来 —— 否则 traceback 会把前面所有诊断挤掉
import sys as _sys
_orig_hook = _sys.excepthook


def _hook(t, v, tb):
    report()
    _orig_hook(t, v, tb)


_sys.excepthook = _hook
raise SystemExit(0 if report() else 1)
