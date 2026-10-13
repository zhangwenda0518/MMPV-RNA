// 实测报告 HTML 中 UpSet 是否真的渲染出来
// 用法: node _probe_upset_render.js <html_path>
const path = require('path');
const { chromium } = require(path.join(
  'C:/Users/17711/.workbuddy-ai/plugins/cache/codebuddy-plugins-official/playwright-cli/0.1.0/node_modules/playwright'));

(async () => {
  const htmlPath = process.argv[2];
  if (!htmlPath) { console.error('usage: node _probe_upset_render.js <html>'); process.exit(1); }
  const url = 'file:///' + htmlPath.replace(/\\/g, '/');

  const EXE = process.env.PW_CHROME ||
    'C:/Users/17711/AppData/Local/ms-playwright/chromium-1234/chrome-win64/chrome.exe';
  const browser = await chromium.launch({ headless: true, executablePath: EXE });
  const page = await browser.newPage({ viewport: { width: 1440, height: 1000 } });

  const logs = [];
  page.on('console', m => logs.push(`[${m.type()}] ${m.text()}`));
  page.on('pageerror', e => logs.push(`[PAGEERROR] ${e.message}`));

  console.log('loading', url);
  await page.goto(url, { waitUntil: 'load', timeout: 180000 });
  await page.waitForTimeout(4000);

  const probe = await page.evaluate(() => {
    const out = {};
    const el = document.getElementById('upset-container');
    out.containerExists = !!el;
    if (el) {
      out.containerRect = el.getBoundingClientRect().toJSON();
      out.containerChildCount = el.children.length;
      out.containerHTMLlen = el.innerHTML.length;
      out.containerHTMLhead = el.innerHTML.slice(0, 400);
      out.svgCount = el.querySelectorAll('svg').length;
      out.canvasCount = el.querySelectorAll('canvas').length;
      // 量墨水: 容器内文本
      out.containerText = (el.textContent || '').slice(0, 300);
    }
    const d = document.getElementById('upset-data');
    out.dataExists = !!d;
    if (d) {
      out.dataLen = d.textContent.length;
      try {
        const arr = JSON.parse(d.textContent);
        out.dataN = arr.length;
        out.sample = arr[0];
        const sets = new Set();
        arr.forEach(x => (x.sets || []).forEach(s => sets.add(s)));
        out.allSets = Array.from(sets);
      } catch (e) { out.dataParseErr = String(e); }
    }
    out.hasUpSetJS = !!window.UpSetJS;
    out.UpSetJSkeys = window.UpSetJS ? Object.keys(window.UpSetJS).slice(0, 20) : null;
    return out;
  });

  console.log('=== PROBE ===');
  console.log(JSON.stringify(probe, null, 2));
  console.log('=== CONSOLE (' + logs.length + ') ===');
  logs.slice(0, 40).forEach(l => console.log(l));

  // 截图 UpSet 区域
  const el = await page.$('#upset-container');
  if (el) {
    await el.screenshot({ path: 'D:/桌面/延伸基因组/MMPV-RNA/scripts/barbarum/_probe_upset.png' });
    console.log('screenshot -> _probe_upset.png');
  }

  await browser.close();
})().catch(e => { console.error('FATAL', e); process.exit(2); });
