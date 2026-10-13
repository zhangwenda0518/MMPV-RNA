// 验证修复: 对现有报告注入「sets 规范化」后 UpSet 是否能渲染
const path = require('path');
const { chromium } = require(path.join(
  'C:/Users/17711/.workbuddy-ai/plugins/cache/codebuddy-plugins-official/playwright-cli/0.1.0/node_modules/playwright'));

(async () => {
  const htmlPath = process.argv[2];
  const url = 'file:///' + htmlPath.replace(/\\/g, '/');
  const browser = await chromium.launch({
    headless: true,
    executablePath: 'C:/Users/17711/AppData/Local/ms-playwright/chromium-1234/chrome-win64/chrome.exe'
  });
  const page = await browser.newPage({ viewport: { width: 1440, height: 1000 } });
  const logs = [];
  page.on('console', m => logs.push(`[${m.type()}] ${m.text()}`));
  page.on('pageerror', e => logs.push(`[PAGEERROR] ${e.message}`));

  await page.goto(url, { waitUntil: 'load', timeout: 180000 });
  await page.waitForTimeout(3000);

  // 注入修复: 规范化 sets 后重新渲染
  const res = await page.evaluate(() => {
    const d = document.getElementById('upset-data');
    const el = document.getElementById('upset-container');
    if (!d || !el || !window.UpSetJS) return { err: 'missing deps' };
    el.innerHTML = '';
    const elems = JSON.parse(d.textContent);
    elems.forEach(function (e) {
      if (!e) return;
      if (e.sets == null) { e.sets = []; }
      else if (typeof e.sets === 'string') { e.sets = [e.sets]; }
      else if (!Array.isArray(e.sets)) { e.sets = Array.from(e.sets); }
    });
    const comb = window.UpSetJS.extractCombinations(elems);
    const _pal = ['#E69F00', '#56B4E9', '#009E73', '#F0E442', '#0072B2', '#D55E00', '#CC79A7', '#000000', '#999999'];
    comb.sets.forEach(function (s, i) { s.color = _pal[i % _pal.length]; });
    window.UpSetJS.renderUpSet(el, {
      sets: comb.sets, combinations: comb.combinations,
      width: 1020, height: 480, padding: 16,
      widthRatios: [0.24, 0.12, 0.64], heightRatios: [0.55, 0.45],
      barLabelOffset: 3, setLabelAlignment: 'center', fontFamily: 'sans-serif',
      exportButtons: false
    });
    return {
      setsCount: comb.sets.length,
      setNames: comb.sets.map(s => s.name || s.set || JSON.stringify(s).slice(0, 40)),
      combCount: comb.combinations.length,
      childCount: el.children.length,
      svgCount: el.querySelectorAll('svg').length,
      innerLen: el.innerHTML.length
    };
  });

  console.log('=== FIX VERIFY ===');
  console.log(JSON.stringify(res, null, 2));
  const errs = logs.filter(l => l.startsWith('[error]') || l.startsWith('[PAGEERROR]'));
  console.log('=== ERRORS (' + errs.length + ') ===');
  errs.slice(0, 10).forEach(l => console.log(l));

  const el = await page.$('#upset-container');
  if (el) {
    await page.setViewportSize({ width: 1440, height: 1400 });
    await el.screenshot({ path: 'D:/桌面/延伸基因组/MMPV-RNA/scripts/barbarum/_probe_upset_FIXED.png' });
    console.log('screenshot -> _probe_upset_FIXED.png');
  }
  await browser.close();
})().catch(e => { console.error('FATAL', e); process.exit(2); });
