/* ═══════════════════════════════════════════════════════════════════
 * Figure Studio ⇄ svguest 桥接
 *
 * 为什么能「零改动」：svguest 的 app.html 只有一个**经典脚本**
 * （既非 type=module 也非 IIFE），所以它的顶层 `function` 声明会挂到
 * window 上、顶层 `const/let` 落在共享的全局词法环境里。另一个 <script>
 * 因此可以直接调用 readFile / svgText / save / A / TABS / setTheme …
 * 上游一行都不用改，升级只要换 vendor_src/app.html 重新构建。
 *
 * 入口：location.hash 里带一个 session token
 *   #s=<token>&theme=light|dark
 *   → GET /api/svg/session/<token> 拿到
 *     { items:[{src,name,url}], back, backBatch, export }
 *
 * 也兼容早期的单文件写法（#open=…&back=…&src=…&name=…）。
 *
 * 提供的动作（**全部通过接管上游 `A` 表实现**，见 hookActions）：
 *   保存回工作台   A.svg  → POST back   （Ctrl+S、File▸Export SVG、Ctrl+K 三入口）
 *   导出矢量 PDF   A.pdf  → POST export （Ctrl+P、File▸Export PDF）
 *   导出 PNG/JPG   A.p300/p600/j300/j600 原生渲染 → dl() 分流器接走字节
 *   导出 SVG 副本  直接 POST svgText()
 *   检查隐藏信息   A.watermark（svguest 自带的水印/元数据检查器）
 *   另存到本地     A.svgAs 保持原生（showSaveFilePicker）
 *
 * 为什么接管 A 而不是拦 keydown：上游主键盘处理器写的是
 *   `s: e.shiftKey ? A.svgAs : A.svg`、`p: e.shiftKey ? A.imp : A.pdf`，
 * 而菜单点击走 `A[b.dataset.a]()`、命令面板也走 A —— 只有 A 表能一次
 * 覆盖三个入口；只在 keydown 上拦 Ctrl+S 会漏掉菜单和命令面板，
 * 它们会去调 showSaveFilePicker 弹文件选择器，与「写回」语义打架。
 * ═══════════════════════════════════════════════════════════════════ */
(function () {
  'use strict';

  // ── 参数 ──────────────────────────────────────────────────────────
  var Q = {};
  (location.hash || '').replace(/^#/, '').split('&').forEach(function (kv) {
    if (!kv) return;
    var i = kv.indexOf('=');
    if (i < 0) return;
    try { Q[decodeURIComponent(kv.slice(0, i))] = decodeURIComponent(kv.slice(i + 1)); }
    catch (e) { Q[kv.slice(0, i)] = kv.slice(i + 1); }
  });

  var S = null;               // session 对象
  var ITEMS = [];             // [{src,name,url}]
  var BACK = Q.back || '/api/svg/writeback';
  var BACK_BATCH = Q.backBatch || '/api/svg/writeback_batch';
  var EXPORT = Q.export || '/api/svg/export';
  var AI_CHAT = Q.aiChat || '/api/ai/chat';
  var AUTOSAVE = Q.autosave || '/api/svg/autosave';
  // 编辑停顿多久后落一份磁盘快照。默认 15 s；`#…&autoMs=1500` 可覆盖，
  // 让端到端测试能真实验证节流链路（而不是只调手动入口）。
  var AUTO_MS = parseInt(Q.autoMs, 10) > 0 ? parseInt(Q.autoMs, 10) : 15000;
  var LOADED_N = {};          // name → 载入时的对象数（防「空画布覆盖原图」）
  var BUSY = false;

  // ── 小工具 ────────────────────────────────────────────────────────
  function toast(t, kind) {
    var el = document.getElementById('fsToast');
    if (!el) {
      el = document.createElement('div');
      el.id = 'fsToast';
      document.body.appendChild(el);
    }
    el.textContent = t;
    el.className = 'on' + (kind === 'bad' ? ' bad' : kind === 'good' ? ' good' : '');
    clearTimeout(el._t);
    el._t = setTimeout(function () { el.className = kind === 'bad' ? 'bad' : ''; }, 4200);
  }

  function post(url, obj) {
    return fetch(url, {
      method: 'POST',
      credentials: 'same-origin',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify(obj)
    }).then(function (r) {
      return r.json().catch(function () { return { error: 'HTTP ' + r.status }; });
    });
  }

  function tabName(i) {
    if (i === TI) return name;
    return (TABS[i] && TABS[i].name) || '';
  }

  function nodeCount(i) {
    try {
      var r = i === TI ? root : (TABS[i] && TABS[i].root);
      return r && r.getElementsByTagName ? r.getElementsByTagName('*').length : 0;
    } catch (e) { return 0; }
  }

  /* 序列化第 i 个标签。
   * 活动标签走 svguest 的 svgText()（会先结束文本编辑/颜色预览）；
   * 后台标签没有进行中的编辑，按 svgRaw() 的口径从快照直接拼 —— 
   * 等价于切过去再存，但不用真的切换，20 个标签也不会卡。 */
  function serializeTab(i) {
    if (i === TI) return svgText();
    var t = TABS[i];
    if (!t || !t.root) return null;
    return '<?xml version="1.0" encoding="UTF-8" standalone="no"?>\n'
      + (t.prolog || '')
      + new XMLSerializer().serializeToString(t.root) + '\n';
  }

  // ── 载入 ──────────────────────────────────────────────────────────
  function fetchBlob(url) {
    return fetch(url, { credentials: 'same-origin' }).then(function (r) {
      if (!r.ok) throw new Error('HTTP ' + r.status);
      return r.blob();
    }).then(function (b) {
      // ⚠⚠ resp.ok 对 204 也成立：实测本机 Chrome 会把 image/svg+xml 的 GET
      //   换成合成 204（无 body、零报错），不校验就会载入空文件，
      //   一保存就把原图覆盖成空白。
      if (!b || b.size < 64) {
        throw new Error('取回的 SVG 是空的（' + (b ? b.size : 0) + ' 字节）'
          + ' —— 检查 /api/svg/raw 的 Content-Type 是否为 image/svg');
      }
      return b;
    });
  }

  async function loadAll() {
    var ok = 0, bad = 0;
    for (var i = 0; i < ITEMS.length; i++) {
      var it = ITEMS[i];
      try {
        var b = await fetchBlob(it.url);
        // readFile 内部会调 tabFresh()：当前标签干净就复用，否则自动开新标签
        await readFile(new File([b], it.name, { type: 'image/svg+xml' }), false);
        LOADED_N[it.name] = nodeCount(TI);
        ok++;
      } catch (e) {
        bad++;
        toast('载入失败 ' + it.name + '：' + e.message, 'bad');
      }
    }
    if (ITEMS.length > 1) {
      toast('已载入 ' + ok + '/' + ITEMS.length + ' 个面板为标签'
        + (bad ? '（' + bad + ' 个失败）' : ''), bad ? 'bad' : 'good');
    } else if (ok) {
      toast('已载入 ' + ITEMS[0].name + '（' + (LOADED_N[ITEMS[0].name] || 0) + ' 个对象）');
    }
    // 载入过程本身会让 dv 自增，这里对齐基线，避免刚打开就触发一次快照
    try { lastDv = dv; } catch (e) { /* 上游没有 dv 就算了 */ }
    syncBar();
  }

  // ── 保存 ──────────────────────────────────────────────────────────
  function collapsed(i, nm) {
    var now = nodeCount(i);
    var was = LOADED_N[nm] || 0;
    return was > 20 && now < Math.max(5, was * 0.1);
  }

  async function saveOne(i, nm) {
    var svg = serializeTab(i);
    if (!svg || svg.indexOf('<svg') < 0) throw new Error('画布是空的');
    if (collapsed(i, nm)) {
      throw new Error('画布只剩 ' + nodeCount(i) + ' 个对象（载入时 '
        + (LOADED_N[nm] || 0) + ' 个），已阻止写回以免把原图覆盖成空白');
    }
    return svg;
  }

  async function saveActive() {
    if (BUSY) return;
    var it = itemOf(TI);
    if (!it) { toast('这个标签不是从工作台打开的，用「另存到本地」', 'bad'); return; }
    BUSY = true;
    try {
      var svg = await saveOne(TI, it.name);
      var j = await post(BACK, { src: it.src, svg: svg });
      if (j.error) throw new Error(j.error);
      var tail = j.backup ? '（旧版备份 ' + String(j.backup).split(/[\\/]/).pop() + '）' : '';
      toast('已保存回工作台' + tail, 'good');
      markClean();
      notify({ fs: 'svguest-saved', src: it.src });
    } catch (e) {
      toast('保存失败：' + e.message, 'bad');
    } finally {
      BUSY = false;
      syncBar();
    }
  }

  function itemOf(i) {
    var nm = tabName(i);
    for (var k = 0; k < ITEMS.length; k++) {
      if (ITEMS[k].name === nm) return ITEMS[k];
    }
    return null;
  }

  async function saveAll() {
    if (BUSY) return;
    BUSY = true;
    try {
      var taken = {}, out = [], skipped = [];
      for (var i = 0; i < TABS.length; i++) {
        var nm = tabName(i), hit = -1;
        for (var k = 0; k < ITEMS.length; k++) {
          if (!taken[k] && ITEMS[k].name === nm) { hit = k; break; }
        }
        if (hit < 0) { skipped.push(nm || ('标签' + (i + 1))); continue; }
        taken[hit] = 1;
        try {
          out.push({ src: ITEMS[hit].src, svg: await saveOne(i, nm) });
        } catch (e) {
          toast(nm + '：' + e.message, 'bad');
        }
      }
      if (!out.length) { toast('没有可保存的标签', 'bad'); return; }
      var j = await post(BACK_BATCH, { items: out });
      var okN = (j.results || []).filter(function (r) { return r.ok; }).length;
      var msg = '已保存 ' + okN + '/' + out.length + ' 个面板回工作台';
      if (skipped.length) msg += '（跳过 ' + skipped.length + ' 个非工作台标签）';
      toast(msg, okN === out.length ? 'good' : 'bad');
      if (okN) markClean();
      notify({ fs: 'svguest-saved-batch', srcs: (j.results || []).filter(function (r) { return r.ok; }).map(function (r) { return r.src; }) });
    } catch (e) {
      toast('批量保存失败：' + e.message, 'bad');
    } finally {
      BUSY = false;
      syncBar();
    }
  }

  function notify(m) {
    try {
      /* 平台已把精修改成**内嵌 iframe**（同源），此时 opener 是 null，
       * 消息要往 parent 发；独立窗口打开的老路径也保留，两条都通。 */
      if (window.parent && window.parent !== window) {
        window.parent.postMessage(m, location.origin);
      } else if (window.opener) {
        window.opener.postMessage(m, location.origin);
      }
    } catch (e) { /* 跨窗口通知失败不影响保存本身 */ }
  }

  /* ⚠ 曾经试过「Esc 收起内嵌浮层」，**已否决**。
   * 实测（探针）上游自己把每一次 Escape 都 preventDefault 掉 —— 它用 Esc
   * 取消当前工具/取消选择/退出文本编辑。冒泡阶段能收到，但
   * e.defaultPrevented 恒为 true。若照此关闭浮层，用户想「取消选择」就会
   * 连整个编辑器一起关掉。所以退出只走浮层顶栏的关闭按钮（平台侧 Esc
   * 仍可用，那是在焦点不在 iframe 里的时候，不存在歧义）。 */

  // ── 导出直通：把 svguest 的导出字节送回工作台，而不是下载 ─────────
  /* 上游所有位图/SVG 导出最后都调 dl(blob, name)。
   * ⚠ 但 dl 不是导出专用的：FigKit 模板（application/json）和
   *   「Kept in this browser」（application/zip）也走它 —— 那些是 svguest
   *   自己的数据导出，落进面板目录毫无意义。所以**按 MIME 分流**，
   *   而不是按调用点分流，也不需要「临时替换再还原」那种时序脆弱的写法。 */
  function installDlRouter() {
    if (window.__fsDlRouted) return;
    var orig = window.dl;
    if (typeof orig !== 'function') return;
    window.__fsDlRouted = true;
    window.dl = function (b, n) {
      var t = String((b && b.type) || '');
      var kind = '';
      if (/^image\/png$/i.test(t)) kind = 'png';
      else if (/^image\/jpe?g$/i.test(t)) kind = 'jpg';
      else if (/^image\/svg/i.test(t)) kind = 'svg';
      if (!kind) return orig(b, n);            // json / zip → 原生下载
      if (!itemOf(TI)) return orig(b, n);      // 非工作台标签 → 原生下载
      blobToB64(b, function (b64) {
        sendExport(b64, kind, guessDpi(n));
      }, function (e) { toast(e.message, 'bad'); });
    };
  }

  /* 从 svguest 的导出名里取 dpi（outName 会拼成 name_300dpi.png）。 */
  function guessDpi(n) {
    var m = /(\d{3,4})\s*dpi/i.exec(String(n || ''));
    return m ? parseInt(m[1], 10) : 0;
  }

  function blobToB64(b, done, fail) {
    var fr = new FileReader();
    fr.onload = function () {
      var s = String(fr.result || '');
      var i = s.indexOf(',');
      done(i < 0 ? '' : s.slice(i + 1));
    };
    fr.onerror = function () { fail(new Error('读取导出数据失败')); };
    fr.readAsDataURL(b);
  }

  function sendExport(b64, kind, dpi) {
    var it = itemOf(TI);
    var src = it ? it.src : '';
    return post(EXPORT, { src: src, kind: kind, dpi: dpi || 0, data: b64 }).then(function (j) {
      if (j.error) throw new Error(j.error);
      toast('已导出到 ' + j.path + '（' + j.human + '）', 'good');
      notify({ fs: 'svguest-exported', path: j.path, src: src, kind: kind });
    }).catch(function (e) { toast('导出失败：' + e.message, 'bad'); });
  }

  /* 位图导出：直接调上游入口，字节由 dl 分流器接走。 */
  function exportRaster(act) {
    if (typeof A === 'undefined' || !A || typeof A[act] !== 'function') {
      toast('导出入口不可用', 'bad');
      return;
    }
    try { A[act](); } catch (e) { toast('导出失败：' + e.message, 'bad'); }
  }

  function exportSvgCopy() {
    var svg = svgText();
    if (!svg || svg.indexOf('<svg') < 0) { toast('画布是空的', 'bad'); return; }
    var b = new Blob([svg], { type: 'image/svg+xml' });
    blobToB64(b, function (b64) {
      sendExport(b64, 'svg', 0);
    }, function (e) { toast(e.message, 'bad'); });
  }

  /* 矢量 PDF：交给服务端用 PyMuPDF 生成。
   * 上游原生 pdf() 是 iframe + window.print()，会弹打印对话框、还得手动
   * 「另存为 PDF」；这里换掉以后 Ctrl+P 直接落盘到原图同目录。 */
  async function exportPdf() {
    var it = itemOf(TI);
    if (!it) { toast('这个标签不是从工作台打开的，用「另存本地」', 'bad'); return; }
    var svg = svgText();
    if (!svg || svg.indexOf('<svg') < 0) { toast('画布是空的', 'bad'); return; }
    toast('正在生成矢量 PDF…');
    var j = await post(EXPORT, { src: it.src, kind: 'pdf', svg: svg });
    if (j.error) { toast('PDF 导出失败：' + j.error, 'bad'); return; }
    toast('已导出矢量 PDF 到 ' + j.path + '（' + j.human + '）', 'good');
    notify({ fs: 'svguest-exported', path: j.path, src: it.src, kind: 'pdf' });
  }

  /* ── 接管 A 表 ────────────────────────────────────────────────────
   * 上游注释原话：「Buttons, menus, keys and the palette all call
   * A[data-a]」—— 菜单点击、快捷键、Ctrl+K 命令面板**三个入口都走 A**。
   * 所以改 A 的属性是唯一能一次覆盖全部入口的挂载点；只在 keydown 上
   * 拦 Ctrl+S 会漏掉「File ▸ Export SVG」和命令面板（它们会去调
   * showSaveFilePicker，弹出文件选择器，与「保存回工作台」语义打架）。 */
  function hookActions() {
    try {
      if (typeof A === 'undefined' || !A || A.__fsHooked) return;
      A.__fsHooked = true;
      A.svg = function () {                       // Ctrl+S / File▸Export SVG
        return TABS.length > 1 ? saveAll() : saveActive();
      };
      A.pdf = function () { return exportPdf(); };// Ctrl+P / File▸Export PDF
    } catch (e) { /* 上游结构变了也不该让整个桥接失效 */ }
  }

  /* 保存后同步 svguest 的「未保存」状态，否则关标签会一直问要不要保存。
   * `dirty` 是上游顶层的 let（不在 window 上，但同一 realm 的脚本可直接赋值）。 */
  function markClean() {
    try { dirty = false; lastDv = dv; } catch (e) { /* 拿不到就算了，不影响写回 */ }
    setAutoLabel('');
  }

  // ── AI 绘图（SVG Chat）走平台代理 ────────────────────────────────
  /* ⚠⚠ 为什么要代理（实测过，别想当然）：
   *   ① **GLM** 的 OPTIONS 预检不回 Access-Control-Allow-Origin（POST 的 401
   *      响应倒是带了）→ 浏览器预检就被拦，直连不可用。上游因此只给 GLM 配
   *      Cloudflare 中继；但那个 Worker 只放行 svguest.online 等三个域名，
   *      127.0.0.1 不在其中，本平台用不了。
   *   ② **DeepSeek** 的预检返回 200 且 ACAO **回显任意 origin** → 浏览器本来
   *      就能直连。走代理是为了让 key 只存在本机磁盘、不落浏览器
   *      localStorage，并统一错误呈现。
   * 上游把发请求收敛在 apiPost(p,url,key,body,o) 一个函数上、返回标准
   * Response，调用方用 r.body 读 SSE —— 所以替换它即可，流式解析、
   * 重试、错误提示全部无感复用。
   * ⚠ apiPost 是**顶层 async function 声明**（不是 let/const），它是全局
   *   对象上的可写属性，所以 `window.apiPost = …` 能生效；换成 let/const
   *   就会变成不可重新赋值的词法绑定，替换会抛 TypeError。 */
  function installAiProxy() {
    try {
      if (typeof apiPost !== 'function' || apiPost.__fsProxied) return false;
      var orig = apiPost;
      function proxied(p, u, key, body, o) {
        return fetch(AI_CHAT, {
          method: 'POST',
          credentials: 'same-origin',
          headers: { 'Content-Type': 'application/json' },
          body: JSON.stringify({ provider: p, url: u, key: key || '', body: body }),
          signal: o && o.signal
        });
      }
      proxied.__fsProxied = 1;
      proxied.__fsOrig = orig;
      window.apiPost = proxied;
      return true;
    } catch (e) {
      return false;                 // 上游改了结构也不该让整个桥接失效
    }
  }

  // ── 磁盘自动快照（Auto-save 的服务端版本）────────────────────────
  /* 上游自带的 Auto-save 存浏览器 IndexedDB（sgDocs，默认开启，10 份/200 MB），
   * 那是崩溃安全网，但换机器 / 清缓存就没了。这里再补一层**落盘**快照：
   *   - 只在编辑停顿时写，**绝不覆盖原文件**（写回仍只由 Ctrl+S 触发）；
   *   - 服务端每源保留最近 12 份并自动剪枝；
   *   - 存在 outputs/_svguest_autosave/，不在面板目录（否则 .svg 会被扫成面板）。
   * 触发信号用上游的 `dv`（document version，每次编辑自增）——比监听 DOM
   * 更省，也比定时比对全文更准；`dirty` 为假时跳过（切换标签不该产生快照）。 */
  var lastDv = -1, snapTimer = 0, snapBusy = false, snapN = 0;

  function setAutoLabel(t) {
    var el = document.getElementById('fsAuto');
    if (el) { el.textContent = t || ''; el.style.display = t ? '' : 'none'; }
  }

  function startAutoSnap() {
    setInterval(function () {
      try {
        if (typeof dv === 'undefined' || typeof dirty === 'undefined') return;
        if (dv === lastDv) return;
        lastDv = dv;
        if (!dirty || !itemOf(TI)) return;   // 没改动 / 非工作台标签 → 不存
        clearTimeout(snapTimer);
        snapTimer = setTimeout(doSnap, AUTO_MS);
      } catch (e) { /* 上游变量名变了就安静退出 */ }
    }, 2000);
  }

  async function doSnap() {
    if (snapBusy) return;
    var it = itemOf(TI);
    if (!it) return;
    var txt;
    try {
      // svgRaw() 是上游专为 autosave 准备的「不结束编辑的序列化」；
      // svgText() 会先提交正在进行的文本编辑，作为兜底。
      txt = (typeof svgRaw === 'function') ? svgRaw() : svgText();
    } catch (e) { return; }
    if (!txt || txt.indexOf('<svg') < 0 || txt.length < 64) return;
    if (collapsed(TI, it.name)) return;      // 空画布不存（否则会存一堆空白）
    snapBusy = true;
    try {
      var j = await post(AUTOSAVE, { src: it.src, svg: txt });
      if (j && !j.error) {
        snapN = j.n || 0;
        setAutoLabel('自动快照 ' + (j.ts || '') + ' · ' + snapN + '/' + (j.keep || 12));
      }
    } catch (e) { /* 自动快照失败不打扰用户，Ctrl+S 才是权威路径 */ }
    finally { snapBusy = false; }
  }

  /* 手动触发一次快照 + 暴露内部状态。
   * 一是给用户「现在就存一份」的入口，二是让 e2e 能断言节流器真的
   * 捕捉到了编辑（lastDv === dv），而不只是测一个手动函数。 */
  window.__fsSnapNow = doSnap;
  window.__fsAutoState = function () {
    try { return { lastDv: lastDv, dv: dv, dirty: dirty, ms: AUTO_MS, n: snapN }; }
    catch (e) { return { err: String(e) }; }
  };

  // ── 浮条 ──────────────────────────────────────────────────────────
  function syncBar() {
    var bar = document.getElementById('fsBar');
    if (!bar) return;
    var n = document.getElementById('fsN');
    if (n) {
      n.textContent = TABS.length > 1 ? (TABS.length + ' 个标签') : '';
      n.style.display = TABS.length > 1 ? '' : 'none';
    }
    var nm = document.getElementById('fsName');
    if (nm) nm.textContent = name || '';
    var bAll = document.getElementById('fsBack');
    if (bAll) {
      bAll.textContent = TABS.length > 1 ? '全部保存' : '保存回工作台';
      bAll.title = TABS.length > 1
        ? '把每个标签写回各自的原文件（Ctrl+S）'
        : '把当前画布写回工作台里的原文件（Ctrl+S）';
    }
  }

  function buildBar() {
    if (document.getElementById('fsBar')) return;
    var d = document.createElement('div');
    d.id = 'fsBar';
    d.innerHTML =
      '<span class="fs-brand">Figure Studio</span>'
      + '<button id="fsBack" class="fs-p">保存回工作台</button>'
      + '<button class="fs-act" id="fsPng300" title="按 300 dpi 渲染 PNG，直接落到原图同目录（不经过下载）">PNG 300</button>'
      + '<button class="fs-act" id="fsPng600" title="按 600 dpi 渲染 PNG，直接落到原图同目录">PNG 600</button>'
      + '<button class="fs-act" id="fsJpg300" title="按 300 dpi 渲染 JPG，直接落到原图同目录">JPG 300</button>'
      + '<button class="fs-act" id="fsPdf" title="生成**矢量** PDF（文字可提取，不是栅格图），直接落到原图同目录">PDF</button>'
      + '<button class="fs-act" id="fsSvg" title="把当前画布另存一份 SVG 副本到原图同目录">SVG 副本</button>'
      + '<button class="fs-act" id="fsWm" title="用 svguest 自带的水印检查器看这张图有没有隐藏的 AI/来源信息">查隐藏信息</button>'
      + '<button class="fs-act" id="fsAi" title="打开 SVG Chat：用自然语言改图（需在右侧「AI 绘图」里配 DeepSeek / GLM 的 key）">AI 改图</button>'
      + '<button class="fs-act" id="fsDl" title="下载到本地（svguest 原生保存）">另存本地</button>'
      + '<span class="fs-auto" id="fsAuto" style="display:none"></span>'
      + '<span class="fs-n" id="fsN" style="display:none"></span>'
      + '<span class="fs-name" id="fsName"></span>'
      + '<button id="fsMin" title="收起/展开" style="padding:5px 9px">−</button>';
    document.body.appendChild(d);

    document.getElementById('fsBack').onclick = function () {
      if (TABS.length > 1) saveAll(); else saveActive();
    };
    document.getElementById('fsPng300').onclick = function () { exportRaster('p300'); };
    document.getElementById('fsPng600').onclick = function () { exportRaster('p600'); };
    document.getElementById('fsJpg300').onclick = function () { exportRaster('j300'); };
    document.getElementById('fsPdf').onclick = exportPdf;
    document.getElementById('fsSvg').onclick = exportSvgCopy;
    document.getElementById('fsWm').onclick = function () {
      try { A.watermark(); } catch (e) { toast('打不开检查器：' + e.message, 'bad'); }
    };
    document.getElementById('fsAi').onclick = function () {
      try { A.chat(); } catch (e) { toast('打不开 SVG Chat：' + e.message, 'bad'); }
    };
    document.getElementById('fsDl').onclick = function () {
      try { save(); } catch (e) { toast('另存失败：' + e.message, 'bad'); }
    };
    document.getElementById('fsMin').onclick = function () {
      d.classList.toggle('fs-min');
      this.textContent = d.classList.contains('fs-min') ? '+' : '−';
    };

    // svguest 换标签/换文件名时刷新浮条
    try {
      var mo = new MutationObserver(function () { syncBar(); });
      var tabsEl = document.getElementById('mTabs');
      if (tabsEl) mo.observe(tabsEl, { childList: true, subtree: true, characterData: true });
    } catch (e) { /* 观察不到就只在动作后刷新 */ }
    syncBar();
  }

  // ── 启动 ──────────────────────────────────────────────────────────
  async function boot() {
    if (Q.theme) {
      try { setTheme(Q.theme === 'dark' ? 1 : 0); } catch (e) { /* 主题失败不致命 */ }
    }
    if (Q.s) {
      try {
        var r = await fetch('/api/svg/session/' + encodeURIComponent(Q.s),
          { credentials: 'same-origin' });
        if (!r.ok) throw new Error('HTTP ' + r.status);
        S = await r.json();
        if (S.error) throw new Error(S.error);
        ITEMS = S.items || [];
        if (S.back) BACK = S.back;
        if (S.backBatch) BACK_BATCH = S.backBatch;
        if (S.export) EXPORT = S.export;
      } catch (e) {
        toast('取会话失败：' + e.message, 'bad');
        return;
      }
    } else if (Q.open) {
      ITEMS = [{ src: Q.src || '', name: Q.name || 'figure.svg', url: Q.open }];
    } else {
      return;                       // 没带参数：当普通 svguest 用，不插浮条
    }
    // ⚠ 顺序：先装路由和动作接管，再建浮条（浮条按钮直接调 A.xxx）
    installDlRouter();
    hookActions();
    installAiProxy();
    buildBar();
    setTimeout(loadAll, 300);
    startAutoSnap();
    // 从工作台「打开 SVG Chat」进来时，直接把 AI 侧栏拉开（等图载入完再开，
    // 否则上游会在空画布上初始化对话上下文）
    if (Q.chat) {
      setTimeout(function () {
        try {
          if (!document.body.classList.contains('chat')) A.chat();
        } catch (e) { toast('打不开 SVG Chat：' + e.message, 'bad'); }
      }, 1200);
    }
  }

  if (document.readyState === 'loading') {
    document.addEventListener('DOMContentLoaded', boot);
  } else {
    boot();
  }
})();
