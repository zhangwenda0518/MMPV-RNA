# -*- coding: utf-8 -*-
"""Figure Studio 前端 UI 冒烟测试（真浏览器 + 真后端）。

跑法（需先启动 app.py）::

    C:/Python312/python.exe figure_studio/uitest.py

判据不是「页面能打开」，而是：
* 面板库真的渲染出卡片（扫描链路通）
* 点卡片后画布真的出现面板块（添加链路通）
* 「横排一行」真的改了块的位置/尺寸（排版链路通）
* 点导出真的拿到产物 URL（渲染链路通）
* SVG 面板走**矢量级**编辑：语义色/图例/字号改动进入预览 URL，重置后消失
* 全程 **零 JS 运行时错误**
"""
from __future__ import annotations

import sys
from pathlib import Path

from playwright.sync_api import sync_playwright

BASE = "http://127.0.0.1:8770"
OUT = Path(__file__).resolve().parent / "outputs" / "_ui"

ok_all = True


def check(label: str, cond: bool, detail: str = "") -> None:
    global ok_all
    ok_all = ok_all and bool(cond)
    print(("  [PASS] " if cond else "  [FAIL] ") + label + ("  " + detail if detail else ""))


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    js_errors: list[str] = []
    http_bad: list[str] = []

    with sync_playwright() as pw:
        try:
            browser = pw.chromium.launch(channel="chrome", headless=True)
            print("  [info] 使用本地 Chrome")
        except Exception as e:                                   # noqa: BLE001
            print("  [info] Chrome channel 不可用，回退自带 Chromium:", e)
            browser = pw.chromium.launch(headless=True)

        page = browser.new_page(viewport={"width": 1680, "height": 950})
        page.on("pageerror", lambda e: js_errors.append("pageerror: " + str(e)))
        page.on("console", lambda m: js_errors.append(f"console.error: {m.text}")
                if m.type == "error" else None)
        page.on("response", lambda r: http_bad.append(f"{r.status} {r.url}")
                if r.status >= 400 else None)

        print("== 1. 打开页面 ==")
        page.goto(BASE, wait_until="load")
        page.wait_for_timeout(600)
        check("标题正确", "Figure Studio" in page.title(), page.title())
        check("三栏骨架存在",
              page.locator(".panel-lib").count() == 1
              and page.locator(".canvas-wrap").count() == 1
              and page.locator(".inspector").count() == 1)

        print("== 2. 面板库扫描 ==")
        try:
            page.wait_for_selector(".lib-card", timeout=20000)
        except Exception:                                        # noqa: BLE001
            pass
        n_cards = page.locator(".lib-card").count()
        check("面板库渲染出卡片", n_cards > 0, f"{n_cards} 张")
        check("缩略图有实际内容",
              page.locator(".lib-card img").first.evaluate("e => e.naturalWidth > 0"))
        page.screenshot(path=str(OUT / "01_library.png"))

        print("== 3. 添加面板到画布 ==")
        for i in range(min(3, n_cards)):
            page.locator(".lib-card").nth(i).click()
            page.wait_for_timeout(350)
        n_items = page.locator(".item").count()
        check("画布出现面板块", n_items == min(3, n_cards), f"{n_items} 个")
        check("面板字母已自动分配",
              page.locator(".item .letter").first.inner_text().strip() != "")
        check("预览图加载成功",
              page.locator(".item img").first.evaluate("e => e.complete && e.naturalWidth > 0"))

        # 记录首块几何，用于验证排版确实改变了它
        geo0 = page.locator(".item").first.evaluate(
            "e => ({l: e.style.left, t: e.style.top, w: e.style.width, h: e.style.height})")
        page.screenshot(path=str(OUT / "02_added.png"))

        print("== 4. 排版工具（横排一行，等比反解）==")
        page.locator('[data-tool="row"]').click()
        page.wait_for_timeout(500)
        geo1 = page.locator(".item").first.evaluate(
            "e => ({l: e.style.left, t: e.style.top, w: e.style.width, h: e.style.height})")
        check("横排改变了块几何", geo0 != geo1, f"{geo0} -> {geo1}")
        # 横排后同一行内各块高度应一致（等比反解的核心性质）
        hs = page.locator(".item").evaluate_all("els => els.map(e => e.style.height)")
        check("同行各块等高（零形变反解）", len(set(hs)) == 1, str(hs))

        print("== 5. 属性面板双向绑定 ==")
        page.locator('[data-tool="autoletter"]').click()
        page.wait_for_timeout(200)
        # 「源图已烘焙字母」场景：关掉开关应清空画布字母
        page.locator('[data-bind="style.show_letters"]').uncheck()
        page.wait_for_timeout(250)
        letters = page.locator(".item .letter").evaluate_all(
            "els => els.map(e => e.textContent)")
        check("关闭「显示字母」后画布无字母", all(not s for s in letters), str(letters))
        page.locator('[data-bind="style.show_letters"]').check()
        page.wait_for_timeout(250)
        letters2 = page.locator(".item .letter").evaluate_all(
            "els => els.map(e => e.textContent)")
        check("重新打开后字母恢复", any(s for s in letters2), str(letters2))

        page.locator('[data-bind="canvas.width_in"]').fill("12")
        page.locator('[data-bind="canvas.width_in"]').dispatch_event("change")
        page.wait_for_timeout(300)
        w_px = page.locator("#canvas").evaluate("e => e.style.width")
        check("改画布宽后画布重算", float(w_px.replace("px", "")) > 0, w_px)
        page.screenshot(path=str(OUT / "03_arranged.png"))

        print("== 6. 裁到内容 ==")
        h_before = page.locator("#canvas").evaluate("e => e.style.height")
        page.locator('[data-tool="trim"]').click()
        page.wait_for_timeout(400)
        h_after = page.locator("#canvas").evaluate("e => e.style.height")
        check("裁到内容后画布收拢", h_before != h_after, f"{h_before} -> {h_after}")

        print("== 7. 调色（图像级）==")
        page.locator(".item").first.click()
        page.wait_for_timeout(300)
        tag = page.locator("#adjTag").inner_text().strip()
        check("选中面板后调色区显示字母", tag not in ("", "未选中"), tag)
        n_sw = page.locator(".pal-swatch").count()
        check("项目调色板已加载色块", n_sw >= 9, f"{n_sw} 个")

        # 灰度：即时 CSS 反馈 + 后端预览 URL 带 adjust
        page.locator('[data-adjt="grayscale"]').check()
        page.wait_for_timeout(200)
        filt = page.locator(".item img").first.evaluate("e => e.style.filter")
        check("灰度有即时反馈（CSS filter）", "grayscale" in filt, filt)
        page.wait_for_timeout(1200)
        src_gray = page.locator(".item img").first.evaluate("e => e.src")
        check("灰度后预览 URL 带 adjust",
              "adjust=" in src_gray and "grayscale" in src_gray,
              src_gray.split("&adjust=")[-1][:70])

        # 取主要颜色 → 生成颜色替换行
        page.locator("#btnLoadColors").click()
        try:
            page.wait_for_selector(".remap-row", timeout=90000)
        except Exception:                                        # noqa: BLE001
            pass
        n_rm = page.locator(".remap-row").count()
        check("「取主要颜色」生成替换行", n_rm > 0, f"{n_rm} 行")

        # 改目标色 → 预览 URL 带 remap
        if n_rm:
            page.locator('.remap-row input[type="color"]').first.evaluate(
                "e => { e.value = '#d55e00';"
                " e.dispatchEvent(new Event('input', {bubbles: true})); }")
            page.wait_for_timeout(1200)
            src_rm = page.locator(".item img").first.evaluate("e => e.src")
            check("改目标色后预览 URL 带 remap", "remap" in src_rm,
                  src_rm.split("&adjust=")[-1][:90])

        # 应用到全部
        page.locator("#btnAdjApplyAll").click()
        page.wait_for_timeout(1200)
        srcs = page.locator(".item img").evaluate_all("els => els.map(e => e.src)")
        check("「应用到全部」覆盖所有面板",
              all("adjust=" in s for s in srcs),
              f"{sum('adjust=' in s for s in srcs)}/{len(srcs)}")
        page.screenshot(path=str(OUT / "05_adjust.png"))

        # 重置本面板
        page.locator(".item").first.click()
        page.wait_for_timeout(250)
        page.locator("#btnAdjReset").click()
        page.wait_for_timeout(1200)
        src_reset = page.locator(".item img").first.evaluate("e => e.src")
        check("重置后预览 URL 去掉 adjust", "adjust=" not in src_reset,
              src_reset.split("?")[-1][:70])

        # 再勾一次灰度，保证下面的导出真的带调色
        page.locator('[data-adjt="grayscale"]').check()
        page.wait_for_timeout(1200)

        # ── 7b. 未选中时的空态引导 ─────────────────────────────
        # ⚠ 历史 bug：未选中时调色控件全 disabled、点击零反馈，用户以为"没动静"。
        #   判据：取消选中后空态提示+遮罩必须出现；点遮罩必须给出文字指引；
        #         重新点选后面板必须恢复。
        print("== 7b. 未选中空态（调色控件禁用的可发现性）==")
        # ⚠ 「裁到内容」后画布与面板等大，画布内部没有空白点 → 必须点到
        #   画布**外侧**的滚动区（target = canvasScroll / canvas-stage 同样取消选中）。
        free = page.evaluate("""() => {
          const sc = document.getElementById('canvasScroll');
          const cv = document.getElementById('canvas');
          if (!sc || !cv) return null;
          const s = sc.getBoundingClientRect(), c = cv.getBoundingClientRect();
          const cands = [
            [c.right + Math.min(24, (s.right - c.right) / 2), s.top + s.height / 2],
            [s.left + 8, s.bottom - 8],
            [s.right - 8, s.top + 8],
          ];
          for (const [x, y] of cands) {
            if (x < s.left || x > s.right || y < s.top || y > s.bottom) continue;
            const el = document.elementFromPoint(x, y);
            if (el && !el.closest('.item')) return [x, y];
          }
          return null;
        }""")
        check("找到画布外空白点（用于取消选中）", free is not None, str(free))
        page.mouse.click(*free)
        page.wait_for_timeout(400)
        check("取消选中后标签变未选中",
              page.locator("#adjTag").inner_text().strip() == "未选中",
              page.locator("#adjTag").inner_text().strip())
        check("空态提示出现",
              page.locator("#adjEmpty").evaluate("e => getComputedStyle(e).display") == "block")
        check("遮罩接住点击（控件全 disabled 也不至于零反馈）",
              page.locator("#adjCover").evaluate("e => getComputedStyle(e).display") == "flex")
        page.locator("#adjCover").click()
        page.wait_for_timeout(300)
        check("点遮罩给出文字指引",
              "点选" in page.locator("#statusText").inner_text(),
              page.locator("#statusText").inner_text()[:50])
        page.screenshot(path=str(OUT / "11_empty_state.png"))
        page.locator(".item").first.click()
        page.wait_for_timeout(500)
        check("重新点选后控件恢复",
              page.locator("#adjTag").inner_text().strip() not in ("", "未选中")
              and page.locator("#adjEmpty").evaluate(
                  "e => getComputedStyle(e).display") == "none")

        # 背景白必须被跳过（96% 白底的图，第一条规则若是 white→white，
        # 用户一改目标色整图变色，看起来像"改崩了"）
        page.locator("#btnLoadColors").click()
        try:
            page.wait_for_selector(".remap-row", timeout=90000)
        except Exception:                                        # noqa: BLE001
            pass
        froms = page.locator(".remap-row .sw").evaluate_all(
            "els => els.map(e => e.title.split(' ')[0])")
        near_white = [c for c in froms
                      if len(c) == 7 and all(int(c[i:i + 2], 16) >= 232
                                             for i in (1, 3, 5))]
        check("自动替换行跳过背景白", not near_white, f"行数 {len(froms)}，近白 {near_white}")

        print("== 7c. 拖右下角手柄缩放（回归：曾把宽度吸成 0，图片凭空消失）==")
        # ⚠⚠ 根因：resize 分支曾写成 snapAxis(it.x_in, nw, …) —— 拿**固定不动的
        #   左边缘**去比对吸附候选，而候选里就有 0 和左边距 m；面板本来就摆在 m 上，
        #   于是必然命中、nw = |m - m| = 0 → 宽高归零、图片消失，**全程零 JS 错误**。
        #   判据必须是「拖完尺寸真的变大且 img 仍有自然尺寸」；只看"没报错"抓不到，
        #   只看"元素还在 DOM 里"也抓不到（宽高为 0 的元素照样在）。
        page.locator(".item").first.click()
        page.wait_for_timeout(600)

        def item_box():
            return page.evaluate("""() => {
              const d = document.querySelector('.item');
              if (!d) return null;
              const r = d.getBoundingClientRect();
              const img = d.querySelector('img');
              return { w: r.width, h: r.height,
                       nw: img.naturalWidth, complete: img.complete };
            }""")

        def drag_handle(dx, dy, steps=5):
            hb = page.locator(".item .handle").first.bounding_box()
            cx, cy = hb["x"] + hb["width"] / 2, hb["y"] + hb["height"] / 2
            page.mouse.move(cx, cy)
            page.mouse.down()
            for i in range(1, steps + 1):
                page.mouse.move(cx + dx * i / steps, cy + dy * i / steps)
                page.wait_for_timeout(70)
            page.mouse.up()
            page.wait_for_timeout(450)

        b0 = item_box()
        drag_handle(72, 46)
        b1 = item_box()
        check("拖手柄放大后尺寸真的变大（不是被吸成 0）",
              b0 and b1 and b1["w"] > b0["w"] * 1.1 and b1["h"] > b0["h"] * 1.1,
              "%.0f×%.0f → %.0f×%.0f" % (b0["w"], b0["h"], b1["w"], b1["h"]))
        check("放大后图片仍可见（naturalWidth>0、宽高非 0）",
              b1 and b1["nw"] > 0 and b1["complete"] and b1["w"] > 10 and b1["h"] > 10,
              "nw=%s w=%.0f h=%.0f" % (b1["nw"], b1["w"], b1["h"]))

        drag_handle(-56, -36)
        b2 = item_box()
        check("反向拖能缩回去（双向可用，且不会归零）",
              b2 and b2["w"] < b1["w"] * 0.95 and b2["w"] > 10,
              "%.0f → %.0f" % (b1["w"], b2["w"]))

        print("== 7d. 多选 / 对齐作用域 / 撤销重做 ==")
        # ⚠⚠ 背景：targets() 的注释一直写着「选中 ≥2 个 → 只对选中的生效」，
        #   但 state.sel 是单个 id，filter 最多命中 1 个 → 那个分支**永远进不去**，
        #   等于死代码：所有排版工具只能作用于全部面板。
        #   这里验证 selSet 之后那句话真的成立 —— 判据必须是「第三个面板没动」，
        #   只看「对齐了」是抓不到这个 bug 的（对全部生效时也会"对齐"）。
        page.locator('[data-tool="row"]').click()          # 先摆成整齐一行
        page.wait_for_timeout(500)
        n_all = page.locator(".item").count()
        check("画布上有 3 个面板可供多选", n_all >= 3, "%d 个" % n_all)

        def xs():
            return page.locator(".item").evaluate_all(
                "els => els.map(e => ({l: parseFloat(e.style.left), "
                "r: parseFloat(e.style.left) + parseFloat(e.style.width)}))")

        # Ctrl+A 全选（先 blur，确保焦点不在输入框 —— 在输入框里 Ctrl+A 是选文字）
        page.evaluate("() => document.activeElement && document.activeElement.blur()")
        page.keyboard.press("Control+a")
        page.wait_for_timeout(300)
        check("Ctrl+A 全选所有面板",
              page.locator(".item.sel").count() == n_all,
              "%d/%d" % (page.locator(".item.sel").count(), n_all))
        check("全选后仍只有一个主选中",
              page.locator(".item.primary").count() == 1)
        check("选中标签显示 +N",
              "+" in (page.locator("#selTag").inner_text() or ""),
              page.locator("#selTag").inner_text())

        # Shift 点掉一个 → 剩 n-1 个
        page.locator(".item").first.click(modifiers=["Shift"])
        page.wait_for_timeout(300)
        n_sel = page.locator(".item.sel").count()
        check("Shift 点选可以取消某个面板", n_sel == n_all - 1,
              "%d 个选中" % n_sel)

        # 右对齐：只作用于选中的那些
        before = xs()
        page.locator('[data-align="right"]').click()
        page.wait_for_timeout(500)
        after = xs()
        # 未选中的那个 = 被 Shift 点掉的那个（DOM 里第一个）
        unsel_moved = abs(after[0]["l"] - before[0]["l"]) > 0.5
        sel_rs = [after[i]["r"] for i in range(1, n_all)]
        check("右对齐只作用于选中的面板（未选中的不动）", not unsel_moved,
              "未选中 l: %.1f → %.1f" % (before[0]["l"], after[0]["l"]))
        check("选中的面板右边缘真的对齐了",
              max(sel_rs) - min(sel_rs) < 1.0,
              "右边缘 %s" % [round(v, 1) for v in sel_rs])

        # 撤销 / 重做
        # ⚠ 判据必须比**整组**左边缘向量，不能只看某一个面板：右对齐的基准是
        #   max(x+w)，本来就是最右的那个不会移动 —— 只盯它会得到永远成立的假阳性。
        def same(a_, b_, tol=0.6):
            return len(a_) == len(b_) and all(abs(a_[i]["l"] - b_[i]["l"]) < tol
                                              for i in range(len(a_)))

        page.keyboard.press("Control+z")
        page.wait_for_timeout(500)
        undo_x = xs()
        check("Ctrl+Z 撤销回到对齐前", same(undo_x, before),
              "l %s vs %s" % ([round(v["l"], 1) for v in undo_x],
                              [round(v["l"], 1) for v in before]))
        page.keyboard.press("Control+Shift+z")
        page.wait_for_timeout(500)
        redo_x = xs()
        check("Ctrl+Shift+Z 重做回到对齐后", same(redo_x, after),
              "l %s vs %s" % ([round(v["l"], 1) for v in redo_x],
                              [round(v["l"], 1) for v in after]))
        check("撤销/重做按钮状态跟着变",
              page.locator("#btnUndo").is_enabled() and not page.locator("#btnRedo").is_enabled())

        # 多选删除
        page.keyboard.press("Delete")
        page.wait_for_timeout(400)
        check("Delete 一次删掉所有选中的面板",
              page.locator(".item").count() == n_all - n_sel,
              "%d → %d" % (n_all, page.locator(".item").count()))
        page.keyboard.press("Control+z")
        page.wait_for_timeout(500)
        check("撤销能恢复被删掉的面板",
              page.locator(".item").count() == n_all,
              "%d 个" % page.locator(".item").count())

        print("== 8. SVG 矢量编辑（颜色 / 图例 / 字号，改的是 XML 不是像素）==")
        svg_dir = Path(__file__).resolve().parent / "outputs" / "_svgdemo"
        page.fill("#rootInput", str(svg_dir))
        page.locator("#btnScan").click()
        try:
            page.wait_for_selector(".lib-card", timeout=60000)
        except Exception:                                        # noqa: BLE001
            pass
        n_svg_cards = page.locator(".lib-card").count()
        check("SVG 目录扫描出卡片", n_svg_cards >= 3, f"{n_svg_cards} 张")

        page.locator(".lib-card").first.click()
        page.wait_for_timeout(700)
        svg_img = page.locator(".item img").last
        svg_src = svg_img.evaluate("e => e.src")
        check("SVG 面板加入画布", ".svg" in svg_src.lower(), svg_src[-58:])
        page.locator(".item").last.click()
        page.wait_for_timeout(400)
        check("SVG 编辑块对 SVG 面板显示", page.locator("#svgBlock").is_visible())

        # 语义色 → 颜色替换行
        page.locator("#btnSvgColors").click()
        try:
            page.wait_for_selector("#svgColorList .remap-row", timeout=90000)
        except Exception:                                        # noqa: BLE001
            pass
        n_svg_col = page.locator("#svgColorList .remap-row").count()
        check("SVG 语义色生成颜色行", n_svg_col >= 3, f"{n_svg_col} 行")
        if n_svg_col:
            page.locator('#svgColorList input[type="color"]').first.evaluate(
                "e => { e.value = '#d55e00';"
                " e.dispatchEvent(new Event('input', {bubbles: true})); }")
            page.wait_for_timeout(1600)
            s2 = page.locator(".item img").last.evaluate("e => e.src")
            check("SVG 改色后预览 URL 带 svg 子树",
                  "svg" in s2 and "%23d55e00" in s2.lower(),
                  s2.split("adjust=")[-1][:80])

        # 图例隐藏
        page.select_option("#svgLegendAction", "hide")
        page.locator("#btnSvgLegend").click()
        page.wait_for_timeout(1600)
        s3 = page.locator(".item img").last.evaluate("e => e.src")
        check("SVG 图例隐藏进入 URL", "legend" in s3,
              s3.split("adjust=")[-1][:80])

        # 字号缩放 + 字体族
        page.fill("#svgFontScale", "1.35")
        page.locator("#btnSvgFontApply").click()
        page.wait_for_timeout(1600)
        s4 = page.locator(".item img").last.evaluate("e => e.src")
        check("SVG 字号缩放进入 URL", "font_scale" in s4,
              s4.split("adjust=")[-1][:80])

        # 重置 → svg 子树消失
        # ⚠ 判据必须只看 adjust 段：URL 的 src 段本身就含 ".svg"（文件名）
        page.locator("#btnSvgReset").click()
        page.wait_for_timeout(1600)
        s5 = page.locator(".item img").last.evaluate("e => e.src")
        adj5 = s5.split("adjust=")[-1] if "adjust=" in s5 else ""
        check("SVG 重置后 svg 子树消失", "svg" not in adj5,
              (adj5 or "(no adjust)")[:70])
        page.screenshot(path=str(OUT / "06_svg_edit.png"))

        # 收拢画布，避免 SVG 面板（9.6×8.75 in 源图）撑出画布
        page.locator('[data-tool="trim"]').click()
        page.wait_for_timeout(500)

        print("== 8b. 精修浮层三档布局（全屏 / 右侧分屏 / 浮动窗口）==")
        # 判据要点：不能只看「点了按钮没报错」。要量**几何**：
        #   分屏 → 宽度显著变窄且 left>0；浮动 → 四边都不贴边；切回全屏 → left 归 0。
        # 最后一条尤其重要：float 模式把 left/top/width/height 写成 inline，
        # 切走时若不清掉，full 就会停在浮动位置（零报错、只是"看着不对"）。
        page.locator(".item").last.click()
        page.wait_for_timeout(600)
        check("选中 SVG 面板后精修按钮可用", page.locator("#btnSgOne").is_enabled())
        page.locator("#btnSgOne").click()
        page.frame(name="svguest").wait_for_selector("#fsBar", timeout=60000)
        page.wait_for_timeout(1200)

        def layer_geom():
            return page.evaluate("""() => {
              const l = document.getElementById('sgLayer');
              const r = l.getBoundingClientRect();
              const f = document.getElementById('sgFrame').getBoundingClientRect();
              return {w: Math.round(r.width), h: Math.round(r.height),
                      left: Math.round(r.left), top: Math.round(r.top),
                      fw: Math.round(f.width), fh: Math.round(f.height),
                      mode: l.dataset.mode};
            }""")

        g_full = layer_geom()
        check("浮层默认全屏", g_full["mode"] == "full" and g_full["left"] == 0,
              "mode=%s left=%d w=%d" % (g_full["mode"], g_full["left"], g_full["w"]))

        page.locator("#btnSgMode").click()
        page.wait_for_timeout(400)
        g_split = layer_geom()
        check("切到右侧分屏：只占右半边（left>0、宽度明显变小）",
              g_split["mode"] == "split" and g_split["left"] > 0
              and g_split["w"] < g_full["w"] * 0.75,
              "mode=%s left=%d w %d→%d" % (g_split["mode"], g_split["left"],
                                          g_full["w"], g_split["w"]))
        check("分屏时 iframe 跟着变窄（不是整块盖住）",
              g_split["fw"] < g_full["fw"] * 0.8,
              "iframe 宽 %d → %d" % (g_full["fw"], g_split["fw"]))
        # 真正的判据是 elementFromPoint —— is_visible() 不考虑遮挡
        check("分屏时左侧工作台没被遮住（能点到）",
              page.evaluate("""() => {
                const el = document.elementFromPoint(30, 60);
                return !!el && !el.closest('#sgLayer');
              }"""),
              page.evaluate("""() => {
                const el = document.elementFromPoint(30, 60);
                return el ? (el.id || el.className) : 'null';
              }"""))
        page.screenshot(path=str(OUT / "07_sg_split.png"))

        page.locator("#btnSgMode").click()
        page.wait_for_timeout(400)
        g_float = layer_geom()
        check("切到浮动窗口：四边都不贴屏幕边",
              g_float["mode"] == "float" and g_float["left"] > 10
              and g_float["top"] > 10 and g_float["w"] < g_full["w"] - 10,
              "mode=%s left=%d top=%d w=%d" % (g_float["mode"], g_float["left"],
                                               g_float["top"], g_float["w"]))

        hb = page.locator("#sgHead").bounding_box()
        page.mouse.move(hb["x"] + 120, hb["y"] + hb["height"] / 2)
        page.mouse.down()
        page.mouse.move(hb["x"] + 220, hb["y"] + hb["height"] / 2 + 50, steps=6)
        page.mouse.up()
        page.wait_for_timeout(300)
        g_moved = layer_geom()
        check("浮动窗口可以拖标题栏移动",
              abs(g_moved["left"] - g_float["left"]) > 40,
              "left %d → %d" % (g_float["left"], g_moved["left"]))

        page.locator("#btnSgMode").click()          # 再循环一次 → 回全屏
        page.wait_for_timeout(400)
        g_back = layer_geom()
        check("循环切回全屏，且清掉了浮动的 inline 定位",
              g_back["mode"] == "full" and g_back["left"] == 0
              and g_back["w"] == g_full["w"],
              "mode=%s left=%d w=%d（期望 left=0 w=%d）"
              % (g_back["mode"], g_back["left"], g_back["w"], g_full["w"]))

        page.locator("#btnSgClose").click()
        page.wait_for_timeout(400)
        check("关闭后浮层收起", not page.locator("#sgLayer.on").count())

        print("== 9. 导出（真渲染，含 SVG 面板的矢量嵌入）==")
        page.fill("#stemInput", "_uitest_fig")
        page.locator("#btnExport").click()
        try:
            page.wait_for_selector("#modalMask.on", timeout=120000)
            modal = page.locator("#modalBox").inner_text()
            check("导出弹窗出现", "导出完成" in modal, modal.replace("\n", " | ")[:120])
            check("PNG 链接存在", page.locator('#modalBox a[href*=".png"]').count() >= 1)
            check("PDF 链接存在", page.locator('#modalBox a[href*=".pdf"]').count() >= 1)
            check("导出提示 SVG 矢量嵌入",
                  ("vector" in modal.lower()) or ("矢量" in modal),
                  modal.replace("\n", " | ")[:140])
        except Exception as e:                                   # noqa: BLE001
            check("导出弹窗出现", False, str(e))
        page.screenshot(path=str(OUT / "04_export.png"))

        print("== 10. JS 运行时错误 / 失败请求 ==")
        check("零 JS 错误", not js_errors, " | ".join(js_errors[:3]))
        check("零 4xx/5xx 请求", not http_bad, " | ".join(http_bad[:3]))

        browser.close()

    print("\n== 结果 ==")
    print("  全部通过 ✔" if ok_all else "  有失败项 ✘")
    return 0 if ok_all else 1


if __name__ == "__main__":
    raise SystemExit(main())
