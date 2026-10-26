# -*- coding: utf-8 -*-
"""交付截图：svguest 精修窗（浮条 + 多标签）+ 工作台的精修区块。

自包含：目录不存在就自己建、素材自己从 `outputs/_svgdemo/` 复制 ——
否则跑过一次清理就再也截不出图（`_svgedit_test` 是测试产物，随时会被删）。
"""
import shutil
from pathlib import Path

from playwright.sync_api import sync_playwright

ROOT = Path(r"D:\桌面\延伸基因组\MMPV-RNA\figure_studio")
TGT = ROOT / "outputs" / "_svgedit_test"
DEMO = ROOT / "outputs" / "_svgdemo"

TGT.mkdir(parents=True, exist_ok=True)
for _src, _dst in (("panel_net.svg", "a_net.svg"), ("panel_map.svg", "b_map.svg")):
    if not (TGT / _dst).exists():
        shutil.copy2(DEMO / _src, TGT / _dst)

# 先造一份自动快照，再开页面。
# ⚠ 必须放在最前面：快照列表有同源去重（src 没变就不重拉），
#   页面选中面板之后再造，列表不会自己刷新，截出来还是空态。
import json
import urllib.request

_A = (TGT / "a_net.svg").resolve()
try:
    _req = urllib.request.Request(
        "http://127.0.0.1:8770/api/svg/autosave",
        data=json.dumps({"src": str(_A),
                         "svg": _A.read_text(encoding="utf-8")}).encode(),
        headers={"Content-Type": "application/json"})
    urllib.request.urlopen(_req, timeout=30).read()
    print("预置快照 OK")
except Exception as _e:                                    # noqa: BLE001
    print("预置快照失败（不影响其余截图）:", _e)

with sync_playwright() as pw:
    b = pw.chromium.launch(channel="chrome", headless=True)
    ctx = b.new_context(viewport={"width": 1600, "height": 950})
    page = ctx.new_page()
    page.goto("http://127.0.0.1:8770", wait_until="load")
    page.wait_for_timeout(500)
    page.fill("#rootInput", str(TGT))
    page.locator("#btnScan").click()
    page.wait_for_selector(".lib-card", timeout=60000)
    page.wait_for_timeout(500)

    # 两张都加进画布 → 演示多标签
    cards = page.locator(".lib-card")
    for i in range(min(2, cards.count())):
        cards.nth(i).click()
        page.wait_for_timeout(700)
    page.wait_for_timeout(900)
    page.screenshot(path=str(ROOT / "outputs" / "_ui" / "18_sg_block.png"))
    print("工作台 SVG 数:", page.evaluate("() => document.querySelectorAll('.item').length"))

    # 精修现在是**内嵌浮层**（同源 iframe），不再弹独立窗口
    page.locator("#btnSgAll").click()
    frame = page.frame(name="svguest")
    assert frame is not None, "iframe #sgFrame 没挂上"
    frame.wait_for_selector("#fsBar", timeout=60000)
    page.wait_for_timeout(7000)
    page.screenshot(path=str(ROOT / "outputs" / "_ui" / "17_svguest_open.png"))
    print("svguest 标签数:", frame.evaluate("() => TABS.length"))
    print("浮条文案:", frame.locator("#fsBack").inner_text(), "|",
          frame.locator("#fsN").inner_text())
    page.locator("#btnSgClose").click()
    page.wait_for_timeout(500)

    # 精修区块 + 备份历史（滚到该区块再截）
    page.wait_for_timeout(900)
    page.locator("#sgBlock").scroll_into_view_if_needed()
    page.wait_for_timeout(600)
    page.locator("aside.inspector").screenshot(path=str(ROOT / "outputs" / "_ui" / "19_sg_history.png"))
    print("历史行数:", page.locator("#sgHistory .sg-row").count())

    # AI 绘图区块（再往下滚一屏）。快照已在页面打开前预置好。
    page.locator("#aiBlock").scroll_into_view_if_needed()
    page.wait_for_timeout(1200)
    page.locator("aside.inspector").screenshot(path=str(ROOT / "outputs" / "_ui" / "20_ai_autosave.png"))
    print("AI 状态:", page.locator("#aiTag").inner_text())
    print("快照行数:", page.locator("#asList .sg-row").count())
    b.close()
