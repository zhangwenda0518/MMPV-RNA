# -*- coding: utf-8 -*-
"""交付截图：自由度/灵活度这批新能力。

  22_multi_select.png  多选（虚线=集合成员，实线粗框=主选中）
  23_sg_split.png      精修浮层「右侧分屏」—— 左边工作台仍可见可点
  24_sg_float.png      精修浮层「浮动窗口」—— 可拖可缩，四边不贴屏
  25_undo_redo.png     撤销/重做与对齐补全的工具栏

自包含：目录不存在就建，素材从 outputs/_svgdemo/ 复制。
"""
import shutil
from pathlib import Path

from playwright.sync_api import sync_playwright

ROOT = Path(r"D:\桌面\延伸基因组\MMPV-RNA\figure_studio")
SRC = ROOT / "outputs" / "_svgdemo"
TGT = ROOT / "outputs" / "_flex_shot"
UI = ROOT / "outputs" / "_ui"

TGT.mkdir(parents=True, exist_ok=True)
UI.mkdir(parents=True, exist_ok=True)
for _s, _d in (("panel_net.svg", "a_net.svg"), ("panel_map.svg", "b_map.svg"),
               ("panel_tracks.svg", "c_tracks.svg")):
    if (SRC / _s).exists() and not (TGT / _d).exists():
        shutil.copy2(SRC / _s, TGT / _d)

with sync_playwright() as pw:
    b = pw.chromium.launch(channel="chrome", headless=True)
    ctx = b.new_context(viewport={"width": 1680, "height": 950})
    page = ctx.new_page()
    page.goto("http://127.0.0.1:8770", wait_until="load")
    page.wait_for_timeout(500)
    page.fill("#rootInput", str(TGT))
    page.locator("#btnScan").click()
    page.wait_for_selector(".lib-card", timeout=60000)
    page.wait_for_timeout(400)

    for i in range(min(3, page.locator(".lib-card").count())):
        page.locator(".lib-card").nth(i).click()
        page.wait_for_timeout(400)
    page.locator('[data-tool="row"]').click()
    page.wait_for_timeout(700)
    page.locator('[data-tool="autoletter"]').click()
    page.wait_for_timeout(500)
    n = page.locator(".item").count()
    print("画布面板数:", n)

    # ── 多选：Ctrl+A 后 Shift 点掉一个，留下「集合成员 + 主选中」两种框 ──
    page.evaluate("() => document.activeElement && document.activeElement.blur()")
    page.keyboard.press("Control+a")
    page.wait_for_timeout(300)
    page.locator(".item").first.click(modifiers=["Shift"])
    page.wait_for_timeout(400)
    print("选中:", page.locator(".item.sel").count(),
          "主选中:", page.locator(".item.primary").count(),
          "标签:", page.locator("#selTag").inner_text())
    page.screenshot(path=str(UI / "22_multi_select.png"))

    # 工具栏（撤销/重做 + 六向对齐）
    page.locator(".canvas-toolbar").screenshot(path=str(UI / "25_undo_redo.png"))

    # ── 精修浮层：分屏 ──
    page.locator(".item").last.click()
    page.wait_for_timeout(600)
    page.locator("#btnSgOne").click()
    page.frame(name="svguest").wait_for_selector("#fsBar", timeout=60000)
    page.wait_for_timeout(3500)
    page.locator("#btnSgMode").click()          # full → split
    page.wait_for_timeout(900)
    print("模式:", page.locator("#sgLayer").get_attribute("data-mode"))
    page.screenshot(path=str(UI / "23_sg_split.png"))

    # ── 浮动窗口 ──
    page.locator("#btnSgMode").click()          # split → float
    page.wait_for_timeout(900)
    print("模式:", page.locator("#sgLayer").get_attribute("data-mode"))
    page.screenshot(path=str(UI / "24_sg_float.png"))

    page.locator("#btnSgClose").click()
    b.close()
print("截图完成")
