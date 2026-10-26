/* Figure Studio · 前端工作台
 *
 * 坐标纪律：state 里一切几何量单位都是**英寸**，原点**左上**（与后端 composer 一致）。
 * 屏幕像素只在这里由 pxPerIn 换算 —— 改缩放/窗口大小不改变 layout 数据。
 * 字号 WYSIWYG：输出图 1 in = 72 pt，故屏幕 px = pt / 72 * pxPerIn。
 */
(function () {
  'use strict';

  // ════════════════════════════════════════════════════════════
  // 状态
  // ════════════════════════════════════════════════════════════
  const DEFAULTS = () => ({
    canvas: { width_in: 10, height_in: 12, dpi: 300 },
    style: {
      bg: 'white', frame: false, frame_color: '#bbbbbb', frame_lw: 0.6,
      letter_fs: 11, letter_dx: 0.0, letter_dy: 1.02, show_letters: true,
      letter_color: '#111111',
      caption_fs: 7, caption_gap_in: 0.06, caption_linespacing: 1.2,
      caption_color: '#222222',
      margin_in: 0.08, pack_margin_in: 0.08, gap_in: 0.12,
    },
    title: {
      text: '', subtitle: '', fs: 13, sub_fs: 9.5, x_in: 0.08, y_in: 0.30,
      color: '#000000', sub_color: '#444444',
    },
    items: [],
  });

  // 调色默认值（必须与后端 composer.ADJUST_DEFAULTS 对齐）
  const ADJ_DEF = {
    brightness: 1, contrast: 1, saturation: 1, sharpness: 1, gamma: 1, hue: 0,
    grayscale: false, invert: false, autocontrast: false,
    threshold: 0, posterize: 0,
    bg: null, flatten: null, flatten_tol: 12,
    tint: null, remap: null, trim: false, trim_tol: 8,
  };

  // sel      = 「主选中」：属性面板 / 调色 / 精修 的目标，永远是最新点中的那一个
  // selSet   = 「多选集合」：对齐 / 均分 / 统一尺寸 / 自动字母 的作用域
  // ⚠ 早先只有 sel 单值，而 targets() 里却写着「选中 ≥2 个 → 只对选中的生效」——
  //   那个分支永远进不去（filter 单 id 最多命中 1 个），等于**死代码**，
  //   所有排版工具只能作用于全部面板。selSet 就是为了让那句话真正成立。
  let state = Object.assign(DEFAULTS(), {
    sel: null, selSet: new Set(), zoom: 1, snap: true, showGrid: true,
  });

  let ppi = 60;                 // px per inch（屏幕）
  let uid = 1;

  const r4 = (v) => Math.round(Number(v) * 10000) / 10000;
  const num = (v, d) => (isFinite(Number(v)) ? Number(v) : d);

  // ════════════════════════════════════════════════════════════
  // DOM
  // ════════════════════════════════════════════════════════════
  const $ = (id) => document.getElementById(id);
  const canvasEl = $('canvas');
  const guidesEl = $('guides');
  const scrollEl = $('canvasScroll');
  const libGrid = $('libGrid');
  const libStat = $('libStat');
  const statusText = $('statusText');
  const statusCanvas = $('statusCanvas');

  function setStatus(msg, kind) {
    statusText.textContent = msg || '';
    statusText.className = kind || '';
  }
  function fmtSize(bytes) { return bytes > 1048576 ? (bytes / 1048576).toFixed(1) + ' MB' : (bytes / 1024).toFixed(0) + ' KB'; }
  function fmtTime(t) { const d = new Date(t * 1000); return d.toLocaleString('zh-CN', { hour12: false }); }

  // ════════════════════════════════════════════════════════════
  // 调色（图像级）
  // ════════════════════════════════════════════════════════════
  function adjOf(it) { return (it && it.adjust) || {}; }

  /** 去掉 UI 私有字段（_off / _frac），只留后端认识的键。空 → null。 */
  function cleanAdjust(a) {
    if (!a) return null;
    const out = {};
    for (const k of Object.keys(a)) {
      if (k === 'remap') {
        const rows = (a.remap || [])
          .filter((r) => r && !r._off && r.from && r.to
            && String(r.from).toLowerCase() !== String(r.to).toLowerCase())
          .map((r) => ({ from: r.from, to: r.to, tol: num(r.tol, 6) }));
        if (rows.length) out.remap = rows;
      } else if (a[k] !== null && a[k] !== undefined && a[k] !== ADJ_DEF[k]) {
        out[k] = a[k];
      }
    }
    return Object.keys(out).length ? out : null;
  }

  function hasAdjust(it) { return !!cleanAdjust(adjOf(it)); }

  function setAdj(it, key, val) {
    if (!it) return;
    const a = Object.assign({}, it.adjust || {});
    const isDefault = (val === null || val === undefined || val === ADJ_DEF[key])
      || (typeof val === 'object' && !Array.isArray(val) && !Object.keys(val).length)
      || (Array.isArray(val) && !val.length);
    if (isDefault) delete a[key]; else a[key] = val;
    if (Object.keys(a).length) it.adjust = a; else delete it.adjust;
    schedulePreview(it);
    markDirty();          // 滑块高频 → 交给防抖，一次连续调节算一步撤销
  }

  /** CSS 滤镜近似（即时反馈用）。gamma/锐度/着色/替换/裁边 CSS 表达不了，
   *  那些改动要等后端预览回来才看得见 —— 后端预览才是准的。 */
  function adjFilter(a) {
    if (!a) return '';
    const p = [];
    if (num(a.brightness, 1) !== 1) p.push('brightness(' + num(a.brightness, 1) + ')');
    if (num(a.contrast, 1) !== 1) p.push('contrast(' + num(a.contrast, 1) + ')');
    if (num(a.saturation, 1) !== 1) p.push('saturate(' + num(a.saturation, 1) + ')');
    if (a.grayscale) p.push('grayscale(1)');
    if (a.invert) p.push('invert(1)');
    if (num(a.hue, 0)) p.push('hue-rotate(' + num(a.hue, 0) + 'deg)');
    return p.join(' ');
  }

  function imgUrl(src, adjust) {
    const c = cleanAdjust(adjust);
    return '/api/img?src=' + encodeURIComponent(src) + '&max=1600' +
      (c ? '&adjust=' + encodeURIComponent(JSON.stringify(c)) : '');
  }

  /** 把调色与文字颜色刷到该面板的 DOM 上（不动几何）。 */
  function applyItemVisual(it) {
    const d = dom.get(it.id);
    if (!d) return;
    d.querySelector('img').style.filter = adjFilter(it.adjust);
    d.querySelector('.letter').style.color = state.style.letter_color || '#111111';
    d.querySelector('.caption').style.color = state.style.caption_color || '#222222';
  }

  /** 调色改动后：先给 CSS 近似反馈，防抖 220ms 再拉后端真预览。 */
  function schedulePreview(it) {
    applyItemVisual(it);
    clearTimeout(it._pt);
    it._pt = setTimeout(async () => {
      const d = dom.get(it.id);
      if (!d) return;
      // 不同 adjust → 不同 URL，浏览器天然区分缓存，不需要额外 cache-bust
      d.querySelector('img').src = imgUrl(it.src, it.adjust);
      if (adjOf(it).trim) await resyncAspect(it);   // 裁边会改比例
      setStatus('调色已更新', 'ok');
    }, 220);
  }

  /** 调色改变了图像尺寸（裁边）后，重算 aspect / 自然尺寸并按锁定比例回正。 */
  async function resyncAspect(it) {
    try {
      const r = await fetch('/api/probe', {
        method: 'POST', headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ src: it.src, dpi: state.canvas.dpi, adjust: cleanAdjust(it.adjust) }),
      });
      const j = await r.json();
      if (j.error) { setStatus('重算比例失败：' + j.error, 'err'); return; }
      it.aspect = j.aspect; it.px_w = j.px_w; it.px_h = j.px_h;
      it.w0_in = j.w_in; it.h0_in = j.h_in;
      if (it.lock !== false) it.h_in = r4(it.w_in * j.aspect);
      renderCanvas(); refreshInspector();
    } catch (e) { /* 忽略 */ }
  }

  // ════════════════════════════════════════════════════════════
  // 画布渲染
  // ════════════════════════════════════════════════════════════
  const dom = new Map();   // itemId -> element

  function computePpi() {
    const availW = Math.max(120, scrollEl.clientWidth - 64);
    const availH = Math.max(120, scrollEl.clientHeight - 64);
    const fit = Math.min(availW / state.canvas.width_in, availH / state.canvas.height_in);
    ppi = Math.max(3, fit) * state.zoom;
    return ppi;
  }

  function renderCanvas() {
    computePpi();
    const W = state.canvas.width_in * ppi;
    const H = state.canvas.height_in * ppi;
    canvasEl.style.width = W + 'px';
    canvasEl.style.height = H + 'px';
    canvasEl.classList.toggle('show-grid', state.showGrid);
    canvasEl.style.backgroundSize = (ppi * 0.25) + 'px ' + (ppi * 0.25) + 'px';
    guidesEl.style.width = W + 'px';
    guidesEl.style.height = H + 'px';
    syncItems();
    updateTitle();
    $('zoomLabel').textContent = Math.round(state.zoom * 100) + '%';
    statusCanvas.textContent =
      `画布 ${state.canvas.width_in} × ${state.canvas.height_in} in @ ${state.canvas.dpi} dpi` +
      ` → ${Math.round(state.canvas.width_in * state.canvas.dpi)} × ${Math.round(state.canvas.height_in * state.canvas.dpi)} px` +
      ` · 面板 ${state.items.length}`;
    $('canvasHint').textContent = state.items.length
      ? (state.sel
          ? '拖动移动 · 右下角方块缩放 · Delete 删除 · 方向键微调（Shift 加速）'
          : '点一个面板图选中它（出现蓝框）→ 右侧即可调色；SVG 面板还有「矢量编辑」')
      : '从左侧面板库单击或拖拽，把面板加进来。';
  }

  function syncItems() {
    const alive = new Set();
    for (const it of state.items) {
      alive.add(it.id);
      let d = dom.get(it.id);
      if (!d) {
        d = document.createElement('div');
        d.className = 'item';
        d.dataset.id = it.id;
        d.innerHTML = '<img alt=""><div class="frame"></div>' +
          '<span class="letter"></span><div class="caption"></div><div class="handle"></div>';
        canvasEl.appendChild(d);
        d.querySelector('img').src = imgUrl(it.src, it.adjust);
        attachItemEvents(d, it.id);
        dom.set(it.id, d);
      }
      const x = it.x_in * ppi, y = it.y_in * ppi, w = it.w_in * ppi, h = it.h_in * ppi;
      d.style.left = x + 'px';
      d.style.top = y + 'px';
      d.style.width = w + 'px';
      d.style.height = h + 'px';
      d.style.zIndex = String(10 + num(it.z, 0));
      // 多选：集合内全部高亮，主选中再多一层更实的框（.primary）
      d.classList.toggle('sel', state.selSet.has(it.id));
      d.classList.toggle('primary', state.sel === it.id);

      const lf = Math.max(6, state.style.letter_fs / 72 * ppi);
      const le = d.querySelector('.letter');
      // show_letters=false 用于「源图已烘焙字母」的场景，避免双字母
      le.textContent = (state.style.show_letters !== false) ? (it.letter || '') : '';
      le.style.fontSize = lf + 'px';
      le.style.left = (state.style.letter_dx * w) + 'px';
      le.style.bottom = 'calc(100% + ' + ((state.style.letter_dy - 1) * h) + 'px)';

      const cf = Math.max(6, state.style.caption_fs / 72 * ppi);
      const ce = d.querySelector('.caption');
      ce.textContent = it.caption || '';
      ce.style.fontSize = cf + 'px';
      ce.style.top = 'calc(100% + ' + (state.style.caption_gap_in * ppi) + 'px)';

      applyItemVisual(it);          // 调色滤镜 + 字母/题注颜色
    }
    for (const [id, el] of dom) {
      if (!alive.has(id)) { el.remove(); dom.delete(id); }
    }
  }

  function updateTitle() {
    let t = canvasEl.querySelector('.canvas-title');
    const txt = (state.title.text || '').trim();
    const sub = (state.title.subtitle || '').trim();
    if (!txt && !sub) { if (t) t.remove(); return; }
    if (!t) {
      t = document.createElement('div');
      t.className = 'canvas-title';
      canvasEl.insertBefore(t, guidesEl);
    }
    t.style.cssText = `position:absolute;left:${state.title.x_in * ppi}px;` +
      `top:${(state.title.y_in - 0.12) * ppi}px;font-weight:700;pointer-events:none;` +
      `color:${state.title.color || '#000000'};` +
      `font-size:${state.title.fs / 72 * ppi}px;line-height:1.25;`;
    t.innerHTML = txt ? escapeHtml(txt) : '';
    if (sub) {
      t.innerHTML += `<div style="font-weight:400;color:${state.title.sub_color || '#444444'};font-size:${state.title.sub_fs / 72 * ppi}px">${escapeHtml(sub)}</div>`;
    }
  }

  function escapeHtml(s) {
    return String(s).replace(/[&<>"']/g, (c) => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' }[c]));
  }

  // ════════════════════════════════════════════════════════════
  // 拖动 / 缩放 / 吸附
  // ════════════════════════════════════════════════════════════
  function attachItemEvents(el, id) {
    el.addEventListener('pointerdown', (e) => {
      const it = state.items.find((x) => x.id === id);
      if (!it) return;
      // Shift / Ctrl / ⌘ = 加选（切换）；普通点选 = 单选
      select(id, e.shiftKey || e.ctrlKey || e.metaKey);
      if (e.target.classList.contains('handle')) return;   // resize 由 handle 处理
      e.preventDefault();
      el.setPointerCapture(e.pointerId);
      dragItem(e, el, it, 'move');
    });
    const h = el.querySelector('.handle');
    h.addEventListener('pointerdown', (e) => {
      const it = state.items.find((x) => x.id === id);
      if (!it) return;
      e.preventDefault(); e.stopPropagation();
      h.setPointerCapture(e.pointerId);
      dragItem(e, el, it, 'resize');
    });
  }

  function dragItem(e0, el, it, mode) {
    const start = { x: e0.clientX, y: e0.clientY };
    const orig = { x: it.x_in, y: it.y_in, w: it.w_in, h: it.h_in };
    const aspect = it.aspect || (orig.h / orig.w);
    let moved = false;

    // 多选整体拖动：被拖的面板属于多选集合且集合 >1 → 整组跟随位移。
    // ⚠ 吸附只按**被拖的那一个**算（把整组当刚体，算完把同一个位移
    //   ddx/ddy 应用给所有成员），否则每个成员各自吸附会把相对位置打散。
    const group = (mode === 'move' && state.selSet.size > 1 && state.selSet.has(it.id))
      ? state.items.filter((i) => state.selSet.has(i.id))
      : [it];
    const origs = group.map((i) => ({ it: i, x: i.x_in, y: i.y_in }));

    function onMove(e) {
      const dx = (e.clientX - start.x) / ppi;
      const dy = (e.clientY - start.y) / ppi;
      if (Math.abs(dx) > 0.001 || Math.abs(dy) > 0.001) moved = true;
      clearGuides();
      if (mode === 'move') {
        let nx = orig.x + dx, ny = orig.y + dy;
        if (state.snap) {
          const sx = snapAxis(nx, it.w_in, it, 'x');
          const sy = snapAxis(ny, it.h_in, it, 'y');
          nx = sx.v; ny = sy.v;
          sx.guide && drawGuide('v', sx.guide);
          sy.guide && drawGuide('h', sy.guide);
        }
        const ddx = nx - orig.x, ddy = ny - orig.y;
        for (const o of origs) {
          o.it.x_in = r4(o.x + ddx);
          o.it.y_in = r4(o.y + ddy);
        }
      } else {
        let nw = Math.max(0.05, orig.w + dx);
        let nh = Math.max(0.05, orig.h + dy);
        if (it.lock !== false) nh = nw * aspect;          // 默认锁定比例
        if (state.snap) {
          /* ⚠⚠ 这里吸附的必须是**右边缘**（x + nw），不是左边缘。
           * snapAxis(v, size) 的语义是「左边缘在 v、宽 size 的盒子找吸附目标」。
           * 早先写成 snapAxis(it.x_in, nw, …) —— 拿**固定不动的左边缘**去比对，
           * 而候选里就有 0 和左边距 m；面板本来就摆在 m 上，于是必然命中，
           * nw = |m - m| = 0 → 面板宽高归零、图片凭空消失，**全程零报错**。
           * 传 (it.x_in + nw, 0) 才是「右边缘吸附到某个位置」；
           * size=0 时候选正好是 0 / 边距 / 画布右边 / 其他面板的两条边。
           * isEdge=false：允许右边缘也吸附到画布中线。 */
          const sx = snapAxis(it.x_in + nw, 0, it, 'x', false);
          if (sx.hit) nw = Math.max(0.05, sx.v - it.x_in);
          if (it.lock !== false) nh = nw * aspect;
        }
        it.w_in = r4(nw); it.h_in = r4(nh);
      }
      renderCanvas();
    }
    function onUp() {
      clearGuides();
      window.removeEventListener('pointermove', onMove);
      window.removeEventListener('pointerup', onUp);
      // 一次拖动/缩放 = 一步撤销（不在这里 push 的话每次 move 都会记一格）
      if (moved) { pushHist(); refreshInspector(); }
    }
    window.addEventListener('pointermove', onMove);
    window.addEventListener('pointerup', onUp);
  }

  function snapAxis(v, size, exclude, axis, isEdge) {
    const TH = 7 / ppi;
    const cands = [];
    const cw = state.canvas.width_in, ch = state.canvas.height_in;
    const m = state.style.margin_in;
    cands.push(0, m, (axis === 'x' ? cw : ch) - size);
    if (!isEdge) cands.push((axis === 'x' ? cw : ch) / 2 - size / 2);
    for (const o of state.items) {
      if (o.id === exclude.id) continue;
      if (axis === 'x') cands.push(o.x_in, o.x_in + o.w_in - size, o.x_in + o.w_in / 2 - size / 2);
      else cands.push(o.y_in, o.y_in + o.h_in - size, o.y_in + o.h_in / 2 - size / 2);
    }
    let best = null, bd = Infinity;
    for (const c of cands) {
      const d = Math.abs(c - v);
      if (d < bd) { bd = d; best = c; }
    }
    if (best !== null && bd <= TH) return { v: r4(best), guide: best, hit: true };
    const g = r4(Math.round(v / 0.05) * 0.05);          // 退化到 0.05 in 网格
    return { v: g, guide: null, hit: false };
  }

  function drawGuide(kind, inchPos) {
    const g = document.createElement('div');
    g.className = 'guide ' + kind;
    if (kind === 'v') g.style.left = (inchPos * ppi) + 'px';
    else g.style.top = (inchPos * ppi) + 'px';
    guidesEl.appendChild(g);
  }
  function clearGuides() { guidesEl.innerHTML = ''; }

  /**
   * 选中。
   * @param id       面板 id；传 null / undefined 表示清空
   * @param additive true = Shift/Ctrl 点选：切换该面板在集合里的去留
   */
  function select(id, additive) {
    if (!id) {
      state.sel = null;
      state.selSet = new Set();
    } else if (additive) {
      const s = new Set(state.selSet);
      if (s.has(id)) {
        s.delete(id);
        // 取消的正好是主选中 → 主选中退到集合里最后一个，没有就置空
        state.sel = s.size ? [...s][s.size - 1] : null;
      } else {
        s.add(id);
        state.sel = id;               // 新加入的成为主选中
      }
      state.selSet = s;
    } else {
      state.sel = id;
      state.selSet = new Set([id]);
    }
    refreshInspector();
    renderCanvas();
  }

  // ════════════════════════════════════════════════════════════
  // 撤销 / 重做（整状态快照栈）
  //   粒度按「一次操作」：拖动/缩放松手记一步，输入类（属性框、滑块）
  //   防抖 350ms 合并成一步 —— 否则拖一次滑块能塞满整个栈。
  // ════════════════════════════════════════════════════════════
  const HIST_MAX = 80;
  let hist = [];
  let histIdx = -1;
  let histLock = false;        // 应用快照期间禁止记录（否则自己触发自己）

  /** 只序列化布局相关字段：拖动时的 _pt 定时器等临时属性不该进快照。 */
  function snapState() {
    return JSON.stringify({
      canvas: state.canvas, style: state.style, title: state.title,
      items: state.items.map((i) => ({
        id: i.id, src: i.src, name: i.name, px_w: i.px_w, px_h: i.px_h,
        w0_in: i.w0_in, h0_in: i.h0_in, aspect: i.aspect, lock: i.lock,
        x_in: i.x_in, y_in: i.y_in, w_in: i.w_in, h_in: i.h_in,
        letter: i.letter, caption: i.caption, z: i.z, adjust: i.adjust,
      })),
    });
  }

  function updateHistUI() {
    const u = $('btnUndo'); const r = $('btnRedo');
    if (u) {
      u.disabled = histIdx <= 0;
      u.title = histIdx > 0 ? `撤销（Ctrl+Z）· 还能退 ${histIdx} 步` : '撤销（Ctrl+Z）· 没有更早的状态';
    }
    if (r) {
      r.disabled = histIdx >= hist.length - 1;
      r.title = histIdx < hist.length - 1 ? '重做（Ctrl+Shift+Z）' : '重做（Ctrl+Shift+Z）· 没有更新的状态';
    }
  }

  function pushHist() {
    if (histLock) return;
    const s = snapState();
    if (hist[histIdx] === s) return;           // 状态没变 → 不占一格
    hist = hist.slice(0, histIdx + 1);         // 撤销后又改 → 丢掉原来的「未来」
    hist.push(s);
    histIdx = hist.length - 1;
    if (hist.length > HIST_MAX) { hist.shift(); histIdx--; }
    updateHistUI();
  }

  let histTimer = 0;
  function markDirty() {
    clearTimeout(histTimer);
    histTimer = setTimeout(pushHist, 350);
  }

  function applyHist(s) {
    histLock = true;
    clearTimeout(histTimer);   // 撤/重做时掐掉挂起的防抖记录，免得事后又记一格
    try {
      const o = JSON.parse(s);
      state.canvas = o.canvas;
      state.style = o.style;
      state.title = o.title;
      state.items = o.items;
      // ⚠ 快照里的 items 是全新对象 → 旧 DOM 必须整体丢掉重建。否则
      //   syncItems 会按 id 复用旧元素，图片与调色停在旧状态（零报错）。
      dom.forEach((el) => el.remove());
      dom.clear();
      const alive = new Set(state.items.map((i) => i.id));
      if (!alive.has(state.sel)) state.sel = null;
      state.selSet = new Set([...state.selSet].filter((id) => alive.has(id)));
      if (state.sel && !state.selSet.size) state.selSet = new Set([state.sel]);
      if (!state.sel) state.selSet = new Set();
    } finally {
      histLock = false;
    }
    renderCanvas();
    refreshInspector();
  }

  function undo() {
    if (histIdx <= 0) { setStatus('没有可撤销的操作', 'err'); return; }
    histIdx--;
    applyHist(hist[histIdx]);
    setStatus(`已撤销（还能退 ${histIdx} 步）`, 'ok');
    updateHistUI();
  }

  function redo() {
    if (histIdx >= hist.length - 1) { setStatus('没有可重做的操作', 'err'); return; }
    histIdx++;
    applyHist(hist[histIdx]);
    setStatus('已重做', 'ok');
    updateHistUI();
  }

  // ════════════════════════════════════════════════════════════
  // 面板库
  // ════════════════════════════════════════════════════════════
  let libPanels = [];

  async function scan() {
    const root = $('rootInput').value.trim();
    if (!root) { setStatus('请先填写面板图目录', 'err'); return; }
    setStatus('扫描中…');
    libGrid.innerHTML = '<div class="lib-empty">扫描中…</div>';
    try {
      const r = await fetch('/api/panels?root=' + encodeURIComponent(root) +
        '&q=' + encodeURIComponent($('searchInput').value.trim()));
      const j = await r.json();
      if (j.error) { setStatus('扫描失败：' + j.error, 'err'); libGrid.innerHTML = '<div class="lib-empty">' + escapeHtml(j.error) + '</div>'; return; }
      libPanels = j.panels || [];
      renderLib();
      libStat.textContent = `${libPanels.length} / ${j.total} 个` + (j.truncated ? '（已截断）' : '');
      setStatus(`扫描完成：${j.total} 个候选，展示 ${libPanels.length} 个`, 'ok');
      loadRoots();
    } catch (e) {
      setStatus('扫描失败：' + e.message, 'err');
      libGrid.innerHTML = '<div class="lib-empty">扫描失败</div>';
    }
  }

  function renderLib() {
    if (!libPanels.length) {
      libGrid.innerHTML = '<div class="lib-empty">没有匹配的面板图</div>';
      return;
    }
    const frag = document.createDocumentFragment();
    for (const p of libPanels) {
      const c = document.createElement('div');
      c.className = 'lib-card';
      c.draggable = true;
      c.title = p.src + `\n${p.px_w} × ${p.px_h} px\n物理 ${p.w_in} × ${p.h_in} in`;
      c.innerHTML = `<img loading="lazy" src="/api/thumb?src=${encodeURIComponent(p.src)}&box=320" alt="">` +
        `<div class="dims">${p.px_w}×${p.px_h}</div>` +
        `<div class="cap">${escapeHtml(p.name)}</div>`;
      c.addEventListener('click', () => addPanel(p));
      c.addEventListener('dragstart', (e) => {
        e.dataTransfer.setData('text/plain', p.src);
        e.dataTransfer.effectAllowed = 'copy';
      });
      frag.appendChild(c);
    }
    libGrid.innerHTML = '';
    libGrid.appendChild(frag);
  }

  async function loadRoots() {
    try {
      const j = await (await fetch('/api/config')).json();
      const dl = $('rootList');
      dl.innerHTML = (j.roots || []).map((r) => `<option value="${escapeHtml(r)}">`).join('');
      if (!state._rootInit && j.last_root) { $('rootInput').value = j.last_root; state._rootInit = true; }
    } catch (e) { /* 忽略 */ }
  }

  function addPanel(p, at) {
    const m = state.style.margin_in;
    const maxW = (state.canvas.width_in - 2 * m) * 0.5;
    const natW = p.w_in || maxW;
    let w = Math.min(maxW, natW);
    if (!isFinite(w) || w <= 0) w = maxW;
    const aspect = p.aspect || 0.6;
    const h = w * aspect;
    const spot = at ? { x: at.x, y: at.y } : freeSpot(w, h);
    const it = {
      id: 'i' + (uid++), src: p.src, rel: p.rel, name: p.name, ext: p.ext,
      px_w: p.px_w, px_h: p.px_h, aspect,
      w0_in: p.w_in, h0_in: p.h_in,
      x_in: r4(spot.x), y_in: r4(spot.y),
      w_in: r4(w), h_in: r4(h),
      letter: nextLetter(), caption: '', z: state.items.length, lock: true,
    };
    state.items.push(it);
    select(it.id);
    pushHist();
    setStatus('已加入：' + p.name, 'ok');
  }

  /* 找一块空地放新面板。
   * 原来固定用 0.18 in 级联 —— 一次拖进/加入好几个面板时它们几乎完全重叠，
   * 用户得逐个挪开（「全部 SVG 面板」这种批量流程尤其明显）。
   * 现在按 0.25 in 网格扫第一个不与现有面板相交的位置；实在没空位才退回级联。 */
  function freeSpot(w, h) {
    const m = state.style.margin_in;
    const CW = state.canvas.width_in;
    const CH = state.canvas.height_in;
    const boxes = state.items.map((i) => ({ x: i.x_in, y: i.y_in, w: i.w_in, h: i.h_in }));
    const hit = (x, y) => boxes.some((b) => !(
      x + w <= b.x + 1e-9 || b.x + b.w <= x + 1e-9 ||
      y + h <= b.y + 1e-9 || b.y + b.h <= y + 1e-9));
    const step = 0.25;
    const xMax = Math.max(m, CW - m - w);
    const yMax = Math.max(m, CH - m - h);
    for (let y = m; y <= yMax + 1e-9; y += step) {
      for (let x = m; x <= xMax + 1e-9; x += step) {
        if (!hit(x, y)) return { x, y };
      }
    }
    const k = state.items.length % 6;
    return { x: m + 0.3 * k, y: m + 0.3 * k };
  }

  function nextLetter() {
    const used = new Set(state.items.map((i) => i.letter).filter(Boolean));
    for (let i = 0; i < 26; i++) {
      const L = String.fromCharCode(65 + i);
      if (!used.has(L)) return L;
    }
    return '';
  }

  // ════════════════════════════════════════════════════════════
  // 属性面板绑定
  // ════════════════════════════════════════════════════════════
  function getPath(obj, path) { return path.split('.').reduce((o, k) => (o == null ? o : o[k]), obj); }
  function setPath(obj, path, val) {
    const ks = path.split('.'); const last = ks.pop();
    let o = obj; for (const k of ks) o = o[k];
    o[last] = val;
  }
  function selItem() { return state.items.find((i) => i.id === state.sel) || null; }

  function bindInspector() {
    document.querySelectorAll('[data-bind]').forEach((el) => {
      const path = el.dataset.bind;
      const h = () => {
        let v = el.type === 'checkbox' ? el.checked : el.value;
        if (el.type === 'number') v = num(v, getPath(state, path));
        setPath(state, path, v);
        renderCanvas();
        markDirty();
      };
      el.addEventListener('input', h);
      el.addEventListener('change', h);
    });
    document.querySelectorAll('[data-sel]').forEach((el) => {
      const key = el.dataset.sel;
      const h = () => {
        const it = selItem(); if (!it) return;
        let v = el.type === 'checkbox' ? el.checked : el.value;
        if (el.type === 'number') v = r4(num(v, it[key]));
        it[key] = v;
        if (key === 'w_in' && it.lock !== false) it.h_in = r4(it.w_in * (it.aspect || 1));
        if (key === 'h_in' && it.lock !== false) it.w_in = r4(it.h_in / (it.aspect || 1));
        renderCanvas();
        markDirty();
      };
      el.addEventListener('input', h);
      el.addEventListener('change', h);
    });
  }

  function refreshInspector() {
    const it = selItem();
    document.querySelectorAll('[data-bind]').forEach((el) => {
      const v = getPath(state, el.dataset.bind);
      if (el.type === 'checkbox') el.checked = !!v; else el.value = v == null ? '' : v;
    });
    document.querySelectorAll('[data-sel]').forEach((el) => {
      const v = it ? it[el.dataset.sel] : '';
      el.disabled = !it;
      if (el.type === 'checkbox') el.checked = it ? v !== false : false;
      else el.value = v == null ? '' : v;
    });
    const nSel = state.selSet.size;
    $('selTag').textContent = it
      ? (it.letter || '—') + (nSel > 1 ? ` +${nSel - 1}` : '')
      : '未选中';
    $('selTag').title = it
      ? (nSel > 1
          ? `已多选 ${nSel} 个：下面编辑的是「${it.name}」；对齐/均分/统一尺寸/自动字母 只作用于这 ${nSel} 个`
          : it.name || '')
      : '在画布上点一个面板（Shift/Ctrl 点击可加选）';
    const bsg = $('btnSvguest');
    if (bsg) {
      bsg.disabled = !isSvg(it);
      bsg.title = isSvg(it)
        ? '用 svguest 交互式精修 ' + (it.name || '') + '，改完保存回工作台'
        : 'svguest 只能精修 SVG 面板（先选中一个 .svg 面板）';
    }
    const bOne = $('btnSgOne');
    if (bOne) {
      bOne.disabled = !isSvg(it);
      bOne.title = isSvg(it)
        ? '用 svguest 精修 ' + (it.name || '') + '（Ctrl+E）'
        : '先选中一个 SVG 面板';
    }
    const nSvg = sgSvgItems().length;
    const bAll = $('btnSgAll');
    if (bAll) {
      bAll.disabled = !nSvg;
      bAll.title = nSvg
        ? `把 ${nSvg} 个 SVG 面板一次性作为标签页打开`
        : '当前布局里没有 SVG 面板';
      bAll.textContent = nSvg > 1 ? `全部 SVG 面板（${nSvg}）` : '全部 SVG 面板';
    }
    const h = $('sgHint');
    if (h) {
      h.textContent = nSvg
        ? `当前布局有 ${nSvg} 个 SVG 面板` +
          (it && isSvg(it) ? `，已选中 ${it.name}` : '，先在画布上点选一个')
        : '当前布局里没有 SVG 面板（PNG/JPG 走上面的调色）';
    }
    $('selInfo').textContent = it
      ? `${it.name}\n源 ${it.px_w}×${it.px_h} px · 比例 ${(1 / (it.aspect || 1)).toFixed(3)}:1 · 面板比 ${(it.h_in / it.w_in).toFixed(4)}` +
        (Math.abs(it.h_in / it.w_in - it.aspect) / it.aspect > 0.01 ? '  ⚠ 与源比例不符（会形变）' : '') +
        (hasAdjust(it) ? `\n调色：${Object.keys(cleanAdjust(it.adjust)).join(' / ')}` : '') +
        (nSel > 1 ? `\n多选：共 ${nSel} 个（排版工具只作用这 ${nSel} 个）` : '')
      : '';
    $('selInfo').style.whiteSpace = 'pre-wrap';
    refreshAdjUI(it);
    if (typeof refreshSvgUI === 'function') refreshSvgUI(it);
    if (typeof refreshSgHistory === 'function') refreshSgHistory();
  }

  // ════════════════════════════════════════════════════════════
  // 排版工具（全部按资产比反解 → 零形变）
  // ════════════════════════════════════════════════════════════
  function targets() {
    // 多选 ≥2 个 → 只对选中的生效；否则（没选 / 只选 1 个）→ 对全部生效。
    // 理由：排版/对齐类操作对单个面板无意义，单选时用户的意图一定是「整幅图」。
    const sel = state.items.filter((i) => state.selSet.has(i.id));
    return sel.length >= 2 ? sel : state.items.slice();
  }

  function arrangeRow(items, y0) {
    if (!items.length) return y0;
    const g = state.style.gap_in, m = state.style.pack_margin_in;
    const AW = state.canvas.width_in - 2 * m - g * (items.length - 1);
    const sumInv = items.reduce((s, i) => s + 1 / (i.aspect || 1), 0);
    const h = AW / sumInv;                       // 行高（所有面板等高）
    let x = m;
    for (const i of items) {
      i.h_in = r4(h); i.w_in = r4(h / (i.aspect || 1));
      i.x_in = r4(x); i.y_in = r4(y0);
      x += i.w_in + g;
    }
    return y0 + h + g;
  }

  function toolArrange(kind) {
    const items = targets();
    if (!items.length) { setStatus('先选中面板（不选则对全部生效）', 'err'); return; }
    const order = items.slice().sort((a, b) => (a.y_in - b.y_in) || (a.x_in - b.x_in));
    const m = state.style.pack_margin_in, g = state.style.gap_in;
    if (kind === 'row') {
      arrangeRow(order, m);
      setStatus('已横排 ' + order.length + ' 个（等比反解，零形变）', 'ok');
    } else if (kind === 'col') {
      const W = state.canvas.width_in - 2 * m;
      let y = m;
      for (const i of order) {
        i.w_in = r4(W); i.h_in = r4(W * (i.aspect || 1));
        i.x_in = r4(m); i.y_in = r4(y);
        y += i.h_in + g;
      }
      setStatus('已纵排 ' + order.length + ' 个', 'ok');
    } else if (kind === 'grid') {
      let y = m, k = 0;
      while (k < order.length) {
        const chunk = order.slice(k, k + 2);
        y = arrangeRow(chunk, y);
        k += 2;
      }
      setStatus('已按 2 列网格排布 ' + order.length + ' 个', 'ok');
    }
    autoGrowCanvas();
    renderCanvas(); refreshInspector(); pushHist();
  }

  function autoGrowCanvas() {
    if (!state.items.length) return;
    const m = state.style.pack_margin_in;
    const maxB = Math.max(...state.items.map((i) => i.y_in + i.h_in));
    const maxR = Math.max(...state.items.map((i) => i.x_in + i.w_in));
    if (maxB + m > state.canvas.height_in) state.canvas.height_in = r4(maxB + m);
    if (maxR + m > state.canvas.width_in) state.canvas.width_in = r4(maxR + m);
  }

  function trimCanvas() {
    if (!state.items.length) { setStatus('画布是空的', 'err'); return; }
    const m = state.style.pack_margin_in;
    const minX = Math.min(...state.items.map((i) => i.x_in));
    const minY = Math.min(...state.items.map((i) => i.y_in));
    const dx = m - minX, dy = m - minY;
    state.items.forEach((i) => { i.x_in = r4(i.x_in + dx); i.y_in = r4(i.y_in + dy); });
    state.title.x_in = r4(state.title.x_in + dx);
    state.title.y_in = r4(state.title.y_in + dy);
    const maxR = Math.max(...state.items.map((i) => i.x_in + i.w_in));
    const maxB = Math.max(...state.items.map((i) => i.y_in + i.h_in));
    state.canvas.width_in = r4(maxR + m);
    state.canvas.height_in = r4(maxB + m);
    renderCanvas(); refreshInspector(); pushHist();
    setStatus(`画布已收到内容边界：${state.canvas.width_in} × ${state.canvas.height_in} in`, 'ok');
  }

  function autoLetter() {
    const items = targets();
    items.sort((a, b) => (a.y_in - b.y_in) || (a.x_in - b.x_in));
    items.forEach((i, k) => { i.letter = k < 26 ? String.fromCharCode(65 + k) : ''; });
    renderCanvas(); refreshInspector(); pushHist();
    setStatus('已按「先上下后左右」重排字母 A…' + (items.length ? items[items.length - 1].letter : ''), 'ok');
  }

  function align(kind) {
    const items = targets();
    if (items.length < 2) { setStatus('至少选中 2 个面板', 'err'); return; }
    if (kind === 'left') { const v = Math.min(...items.map((i) => i.x_in)); items.forEach((i) => i.x_in = r4(v)); }
    if (kind === 'right') {
      const v = Math.max(...items.map((i) => i.x_in + i.w_in));
      items.forEach((i) => i.x_in = r4(v - i.w_in));
    }
    if (kind === 'top') { const v = Math.min(...items.map((i) => i.y_in)); items.forEach((i) => i.y_in = r4(v)); }
    if (kind === 'bottom') {
      const v = Math.max(...items.map((i) => i.y_in + i.h_in));
      items.forEach((i) => i.y_in = r4(v - i.h_in));
    }
    // 居中按**画布**算（不是按选中组包围盒）：这才是「放到页面正中间」的直觉
    if (kind === 'hcenter') {
      const c = state.canvas.width_in / 2;
      items.forEach((i) => i.x_in = r4(c - i.w_in / 2));
    }
    if (kind === 'vcenter') {
      const c = state.canvas.height_in / 2;
      items.forEach((i) => i.y_in = r4(c - i.h_in / 2));
    }
    if (kind === 'hdist') {
      const s = items.slice().sort((a, b) => a.x_in - b.x_in);
      const L = s[0].x_in, R = s[s.length - 1].x_in + s[s.length - 1].w_in;
      const total = s.reduce((a, i) => a + i.w_in, 0);
      const gap = (R - L - total) / (s.length - 1);
      let x = L;
      for (const i of s) { i.x_in = r4(x); x += i.w_in + gap; }
    }
    if (kind === 'vdist') {
      const s = items.slice().sort((a, b) => a.y_in - b.y_in);
      const T = s[0].y_in, B = s[s.length - 1].y_in + s[s.length - 1].h_in;
      const total = s.reduce((a, i) => a + i.h_in, 0);
      const gap = (B - T - total) / (s.length - 1);
      let y = T;
      for (const i of s) { i.y_in = r4(y); y += i.h_in + gap; }
    }
    renderCanvas(); refreshInspector(); pushHist();
    setStatus(`对齐完成（作用域：${items.length} 个面板）`, 'ok');
  }

  function sameSize(kind) {
    const items = targets();
    if (items.length < 2) { setStatus('至少选中 2 个面板', 'err'); return; }
    if (kind === 'sameh') {
      const h = Math.max(...items.map((i) => i.h_in));
      items.forEach((i) => { i.h_in = r4(h); i.w_in = r4(h / (i.aspect || 1)); });
    } else {
      const w = Math.max(...items.map((i) => i.w_in));
      items.forEach((i) => { i.w_in = r4(w); i.h_in = r4(w * (i.aspect || 1)); });
    }
    renderCanvas(); refreshInspector(); pushHist();
  }

  function fit1to1() {
    const items = targets();
    items.forEach((i) => {
      if (i.w0_in && i.h0_in) { i.w_in = r4(i.w0_in); i.h_in = r4(i.h0_in); }
    });
    renderCanvas(); refreshInspector(); pushHist();
    setStatus('已恢复为源图物理尺寸', 'ok');
  }

  function itemAction(kind) {
    const it = selItem();
    // del 支持多选批量：有集合就删集合，否则删主选中
    if (kind === 'del') {
      const ids = state.selSet.size ? new Set(state.selSet) : new Set(it ? [it.id] : []);
      if (!ids.size) { setStatus('先选中面板', 'err'); return; }
      state.items = state.items.filter((i) => !ids.has(i.id));
      ids.forEach((id) => { const el = dom.get(id); if (el) { el.remove(); dom.delete(id); } });
      state.sel = null; state.selSet = new Set();
      renderCanvas(); refreshInspector(); pushHist();
      setStatus(`已删除 ${ids.size} 个面板`, 'ok');
      return;
    }
    if (!it) { setStatus('先选中一个面板', 'err'); return; }
    if (kind === 'top') it.z = Math.max(0, ...state.items.map((i) => num(i.z, 0))) + 1;
    if (kind === 'bottom') it.z = Math.min(0, ...state.items.map((i) => num(i.z, 0))) - 1;
    if (kind === 'dup') {
      // 多选时整组复制（保持相对位置），否则只复制主选中
      const src = (state.selSet.size > 1)
        ? state.items.filter((i) => state.selSet.has(i.id))
        : [it];
      const made = [];
      for (const s of src) {
        const c = Object.assign({}, s, {
          id: 'i' + (uid++), x_in: r4(s.x_in + 0.2), y_in: r4(s.y_in + 0.2),
          letter: nextLetter(),
          // ⚠ 必须深拷贝 adjust：浅拷贝会让两个面板共用同一个对象，改一个动两个
          adjust: s.adjust ? JSON.parse(JSON.stringify(s.adjust)) : undefined,
        });
        state.items.push(c); made.push(c);
      }
      state.sel = made[made.length - 1].id;
      state.selSet = new Set(made.map((c) => c.id));
      setStatus(`已复制 ${made.length} 个面板`, 'ok');
    }
    if (kind === 'nat') { it.h_in = r4(it.w_in * (it.aspect || 1)); }
    renderCanvas(); refreshInspector(); pushHist();
  }

  // ════════════════════════════════════════════════════════════
  // 调色面板
  // ════════════════════════════════════════════════════════════
  let lastColorInput = null;   // 最近聚焦的颜色框 —— 调色板点击时往这里填
  let palNote = '';

  // SVG 矢量编辑：状态读写（编辑发生在 XML 层，不是像素层）
  const isSvg = (it) => !!it && /\.svg$/i.test(it.src || '');
  function svgOf(it) { return (it && it.adjust && it.adjust.svg) || {}; }

  /** 合并进 it.adjust.svg；null / 空对象 / 空数组 表示删除该键。 */
  function setSvg(it, patch) {
    if (!it) return;
    const cur = Object.assign({}, svgOf(it));
    for (const k of Object.keys(patch || {})) {
      const v = patch[k];
      const empty = (v === null || v === undefined
        || (typeof v === 'object' && !Array.isArray(v) && !Object.keys(v).length)
        || (Array.isArray(v) && !v.length));
      if (empty) delete cur[k]; else cur[k] = v;
    }
    const a = Object.assign({}, it.adjust || {});
    if (Object.keys(cur).length) a.svg = cur; else delete a.svg;
    if (Object.keys(a).length) it.adjust = a; else delete it.adjust;
    schedulePreview(it);
    refreshSvgUI(it);
  }

  function selItem() { return state.items.find((i) => i.id === state.sel) || null; }

  function refreshSvgUI(it) {
    const blk = $('svgBlock');
    if (!blk) return;
    const on = isSvg(it);
    blk.style.display = on ? '' : 'none';
    if (!on) { $('svgTag').textContent = ''; return; }
    const e = svgOf(it);
    const n = Object.keys(e).length;
    $('svgTag').textContent = (it.letter || '--') + (n ? ' . edited' : '');
    $('svgFontScale').value = num(e.font_scale, 1);
    $('svgStrokeScale').value = num(e.stroke_scale, 1);
    $('svgFontFamily').value = e.font_family || '';
    $('svgTextColor').value = e.text_color || '#222222';
    $('svgLegendAction').value = (e.legend && e.legend.action) || 'hide';
    $('svgLegendDx').value = num(e.legend && e.legend.dx, 0);
    $('svgLegendDy').value = num(e.legend && e.legend.dy, 0);
    $('svgLegendFs').value = num(e.legend && e.legend.font_scale, 1);
    renderSvgRows('svgColorList', e.colors || [],
      'no colour rule yet -- click "Load SVG semantic colours".');
    renderSvgRows('svgTextList', e.texts || [],
      'no text rule yet -- click "Load text list".');
    const hides = e.hide || [];
    $('svgHideList').innerHTML = hides.length
      ? hides.map((h) => `<span class="chip">${escapeHtml(h)}` +
        `<b data-hide="${escapeHtml(h)}" title="unhide">x</b></span>`).join('')
      : '<div class="hint">nothing hidden</div>';
    $('svgHideList').querySelectorAll('b[data-hide]').forEach((b) => {
      b.addEventListener('click', () => {
        const t = b.getAttribute('data-hide');
        setSvg(it, { hide: (svgOf(it).hide || []).filter((x) => x !== t) });
      });
    });
    if ($('svgMeta').dataset.src !== it.src) {
      $('svgMeta').dataset.src = it.src;
      fetch('/api/svg/inspect', {
        method: 'POST', headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ src: it.src }),
      }).then((r) => r.json()).then((j) => {
        if (j.error) { $('svgMeta').textContent = j.error; return; }
        const ins = (j.size_in || []).map((v) => v.toFixed(2)).join(' x ');
        const leg = (j.legend || []).map((g) => g.id || '(no id)').join(', ');
        $('svgMeta').textContent =
          `canvas ${ins} in | from ${j.flavor} | text ${j.n_text} | ` +
          `groups ${j.n_groups} | legend: ${leg || 'none'}`;
      }).catch((err) => { $('svgMeta').textContent = 'inspect failed: ' + err; });
    }
  }

  function renderSvgRows(boxId, rows, emptyText) {
    const box = $(boxId);
    const isColor = boxId === 'svgColorList';
    if (!rows.length) {
      box.innerHTML = `<div class="remap-empty">${emptyText}</div>`;
      return;
    }
    box.innerHTML = '';
    rows.forEach((r, i) => {
      const d = document.createElement('div');
      d.className = 'remap-row';
      const sw = isColor
        ? `<span class="sw" style="background:${escapeHtml(r.from)}" title="${escapeHtml(r.from)}"></span>`
        : `<span class="sw sw-text" title="text">T</span>`;
      const scope = isColor
        ? `<select title="scope"><option value="all">all</option>` +
          `<option value="fill"${r.scope === 'fill' ? ' selected' : ''}>fill</option>` +
          `<option value="stroke"${r.scope === 'stroke' ? ' selected' : ''}>stroke</option></select>`
        : '';
      const editor = isColor
        ? `<input type="color" value="${escapeHtml(r.to || r.from)}" title="target colour">`
        : `<input type="text" value="${escapeHtml(r.to || '')}" placeholder="new text" style="flex:1">`;
      const frm = r.from.length > 20 ? r.from.slice(0, 20) + '...' : r.from;
      d.innerHTML =
        `<input type="checkbox" ${r._off ? '' : 'checked'} title="enable / disable">` +
        sw + `<span class="arrow">to</span>` + editor +
        `<span class="frm" title="${escapeHtml(r.from)}">${escapeHtml(frm)}</span>` +
        scope + `<button class="del" title="delete">x</button>`;
      const cb = d.querySelector('input[type=checkbox]');
      // ⚠ 不能用 d.children[2]：那是箭头 span，真正的输入框是 children[3]
      const ed = d.querySelector('input[type=color], input[type=text]');
      const sel = d.querySelector('select');
      const del = d.querySelector('.del');
      cb.addEventListener('change', (e2) => {
        r._off = !e2.target.checked; commitSvgRows(boxId, rows);
      });
      ed.addEventListener('input', (e2) => {
        r.to = e2.target.value; commitSvgRows(boxId, rows);
      });
      if (sel) sel.addEventListener('change', (e2) => {
        r.scope = e2.target.value; commitSvgRows(boxId, rows);
      });
      del.addEventListener('click', () => {
        rows.splice(i, 1); commitSvgRows(boxId, rows.length ? rows : null);
      });
      box.appendChild(d);
    });
  }

  function commitSvgRows(boxId, rows) {
    const it = selItem();
    if (!it) return;
    const key = boxId === 'svgColorList' ? 'colors' : 'texts';
    setSvg(it, { [key]: rows && rows.length ? rows : null });
  }

  async function svgLoad(kind) {
    const it = selItem();
    if (!it) { setStatus('select a panel first', 'err'); return; }
    if (!isSvg(it)) {
      setStatus('not an SVG -- vector editing only applies to .svg; ' +
        'use the image-level panel above for PNG/PDF', 'err');
      return;
    }
    setStatus('reading SVG structure...');
    try {
      const j = await fetch('/api/svg/inspect', {
        method: 'POST', headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ src: it.src }),
      }).then((r) => r.json());
      if (j.error) { setStatus('read failed: ' + j.error, 'err'); return; }
      if (kind === 'colors') {
        const cur = (svgOf(it).colors || []).slice();
        const used = new Set(cur.map((x) => String(x.from).toLowerCase()));
        (j.colors || []).slice(0, 14).forEach((c) => {
          if (used.has(c.hex.toLowerCase()) || c.is_bg) return;
          cur.push({ from: c.hex, to: c.hex, scope: 'all', _frac: c.frac });
        });
        setSvg(it, { colors: cur });
        setStatus(`generated ${cur.length} colour rules (near-white skipped)`, 'ok');
      } else {
        const cur = (svgOf(it).texts || []).slice();
        const used = new Set(cur.map((x) => x.from));
        (j.texts || []).slice(0, 40).forEach((t) => {
          if (used.has(t.text) || !t.text) return;
          cur.push({ from: t.text, to: t.text });
        });
        setSvg(it, { texts: cur });
        setStatus(`generated ${cur.length} text rules`, 'ok');
      }
    } catch (err) {
      setStatus('read failed: ' + err, 'err');
    }
  }

  function adjValText(k, v) {
    if (k === 'hue') return Math.round(v) + '°';
    if (k === 'flatten_tol' || k === 'trim_tol') return String(Math.round(v));
    return Number(v).toFixed(2);
  }

  function refreshAdjUI(it) {
    const a = adjOf(it);
    $('adjTag').textContent = it
      ? ((it.letter || '—') + (hasAdjust(it) ? ' · 已调色' : ''))
      : '未选中';
    // 未选中时：整个调色块进空态（控件全部 disabled，点了没反应也不提示）
    const blk = $('adjBlock');
    if (blk) blk.classList.toggle('empty', !it);
    const em = $('adjEmpty');
    if (em) em.style.display = it ? 'none' : '';

    document.querySelectorAll('[data-adj]').forEach((el) => {
      const k = el.dataset.adj;
      const v = (k === 'tint_amount')
        ? (a.tint ? num(a.tint.amount, 0.3) : 0)
        : num(a[k], ADJ_DEF[k]);
      el.value = v;
      el.disabled = !it;
      const val = el.parentElement.querySelector('.val');
      if (val) val.textContent = adjValText(k, v);
    });
    document.querySelectorAll('[data-adjt]').forEach((el) => {
      el.checked = !!(it && a[el.dataset.adjt]);
      el.disabled = !it;
    });
    document.querySelectorAll('[data-adjn]').forEach((el) => {
      const k = el.dataset.adjn;
      el.value = (it && a[k] !== undefined) ? a[k] : ADJ_DEF[k];
      el.disabled = !it;
    });
    document.querySelectorAll('[data-adjc]').forEach((el) => {
      const k = el.dataset.adjc;
      let v = '#ffffff';
      if (k === 'tint') v = (a.tint && a.tint.color) || '#0072B2';
      else if (a[k]) v = a[k];
      el.value = v;
      el.disabled = !it;
    });
    document.querySelectorAll('[data-adjc-on]').forEach((el) => {
      el.checked = !!(it && a[el.dataset.adjcOn]);
      el.disabled = !it;
    });

    renderRemap(it);

    const c = cleanAdjust(a) || {};
    const notes = [];
    if (c.trim) notes.push('「裁白边」改变面板比例，已自动重算宽高');
    if (c.threshold) notes.push('二值化会丢弃彩色信息');
    if (c.posterize) notes.push('色调分离减少色阶');
    $('adjHint').textContent = notes.join('；') || palNote;
  }

  function renderRemap(it) {
    const box = $('remapList');
    const rows = (it && adjOf(it).remap) || [];
    if (!rows.length) {
      box.innerHTML = '<div class="remap-empty">还没有替换规则 —— 点「取主要颜色」按面板主色自动生成'
        + '（背景白会自动跳过，避免一改目标色整图变色），再把右侧色块改成目标色即可'
        + '（from == to 的规则不会生效）。</div>';
      return;
    }
    box.innerHTML = '';
    rows.forEach((r, i) => {
      const d = document.createElement('div');
      d.className = 'remap-row';
      const frac = (r._frac != null) ? `  ${(r._frac * 100).toFixed(1)}%` : '';
      d.innerHTML =
        `<input type="checkbox" ${r._off ? '' : 'checked'} title="启用 / 停用">` +
        `<span class="sw" style="background:${escapeHtml(r.from)}" title="${escapeHtml(r.from)}${frac}"></span>` +
        `<span class="arrow">→</span>` +
        `<input type="color" value="${escapeHtml(r.to || r.from)}" title="目标色">` +
        `<input type="number" min="0" max="80" step="1" value="${num(r.tol, 6)}" title="容差（0 = 精确匹配）">` +
        `<button class="del" title="删除">×</button>`;
      d.querySelector('input[type=checkbox]').addEventListener('change', (e) => {
        r._off = !e.target.checked;
        setAdj(it, 'remap', rows.slice());
      });
      d.querySelector('input[type=color]').addEventListener('input', (e) => {
        r.to = e.target.value;
        setAdj(it, 'remap', rows.slice());
      });
      d.querySelector('input[type=number]').addEventListener('input', (e) => {
        r.tol = num(e.target.value, 6);
        setAdj(it, 'remap', rows.slice());
      });
      d.querySelector('.del').addEventListener('click', () => {
        rows.splice(i, 1);
        setAdj(it, 'remap', rows.length ? rows.slice() : null);
        renderRemap(it);
      });
      box.appendChild(d);
    });
  }

  async function loadColors(it) {
    setStatus('提取主要颜色…');
    try {
      const r = await fetch('/api/colors', {
        method: 'POST', headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ src: it.src, n: 14, dpi: state.canvas.dpi }),
      });
      const j = await r.json();
      if (j.error) { setStatus('提取失败：' + j.error, 'err'); return; }
      const cur = adjOf(it).remap || [];
      const used = new Set(cur.map((x) => String(x.from).toLowerCase()));
      // ⚠⚠ 背景白必须跳过：实测 matplotlib 图白底占比 96%，它是"最大主色"，
      //   自动生成的第一条规则就是 white→white；用户一把目标色改成彩色，
      //   整张图瞬间变彩色，看起来像"改崩了"。改背景请走「透明底色 / 近白替换」。
      const nearWhite = (hex) => {
        const m = /^#?([0-9a-f]{6})$/i.exec(String(hex || ''));
        if (!m) return false;
        const v = parseInt(m[1], 16);
        return ((v >> 16) & 255) >= 232 && ((v >> 8) & 255) >= 232 && (v & 255) >= 232;
      };
      const all = j.colors || [];
      const skipBg = all.filter((c) => nearWhite(c.hex)).length;
      const add = all
        .filter((c) => !used.has(c.hex.toLowerCase()) && !nearWhite(c.hex))
        .slice(0, 14)
        .map((c) => ({ from: c.hex, to: c.hex, tol: 6, _frac: c.frac }));
      if (!add.length) {
        setStatus('没有新的主要颜色了' + (skipBg ? '（背景白已跳过）' : ''));
        return;
      }
      setAdj(it, 'remap', cur.concat(add));
      renderRemap(it);
      setStatus(`提取到 ${add.length} 个主要颜色（源图 ${j.n_buckets} 个色桶 · ${j.mode} 模式`
        + (skipBg ? ' · 已跳过背景白，改背景请用「透明底色 / 近白替换」' : '') + '）', 'ok');
    } catch (e) {
      setStatus('提取失败：' + e.message, 'err');
    }
  }

  function applyAdjToAll() {
    const it = selItem();
    if (!it) { setStatus('先选中一个面板', 'err'); return; }
    const src = cleanAdjust(it.adjust);
    if (!src) { setStatus('当前面板没有调色可复制', 'err'); return; }
    let n = 0;
    state.items.forEach((o) => {
      if (o.id === it.id) return;
      o.adjust = JSON.parse(JSON.stringify(src));
      schedulePreview(o);
      n++;
    });
    renderCanvas();
    pushHist();
    setStatus(`已把调色复制到另外 ${n} 个面板`, 'ok');
  }

  async function loadPalette() {
    try {
      const j = await (await fetch('/api/palette')).json();
      const strip = $('palStrip');
      strip.innerHTML = '';
      palNote = j.note || '';
      strip.title = (j.name || '') + (j.path ? '\n' + j.path : '') + (j.note ? '\n' + j.note : '');
      (j.groups || []).forEach((g, gi) => {
        if (gi) {
          const s = document.createElement('span');
          s.className = 'sep';
          strip.appendChild(s);
        }
        const semantic = /语义/.test(g.name || '');
        g.colors.forEach((c) => {
          const b = document.createElement('span');
          b.className = 'pal-swatch' + (semantic ? ' sem' : '');
          b.style.background = c.hex;
          b.title = `${c.name}  ${c.hex}\n点击填入「最近点过的颜色框」`;
          b.addEventListener('click', () => {
            if (lastColorInput && !lastColorInput.disabled) {
              lastColorInput.value = c.hex;
              lastColorInput.dispatchEvent(new Event('input', { bubbles: true }));
              setStatus('已填入 ' + c.hex, 'ok');
            } else {
              setStatus('先点一个颜色框，再点色板', 'err');
            }
          });
          strip.appendChild(b);
        });
      });
    } catch (e) { /* 忽略 */ }
  }

  function bindAdj() {
    document.querySelectorAll('[data-adj]').forEach((el) => {
      const k = el.dataset.adj;
      el.addEventListener('input', () => {
        const it = selItem(); if (!it) return;
        const v = num(el.value, 1);
        if (k === 'tint_amount') {
          const t = Object.assign({ color: '#0072B2', amount: 0.3 }, adjOf(it).tint || {});
          t.amount = v;
          setAdj(it, 'tint', v > 0 ? t : null);
        } else {
          setAdj(it, k, v);
        }
        const val = el.parentElement.querySelector('.val');
        if (val) val.textContent = adjValText(k, v);
      });
    });
    document.querySelectorAll('[data-adjt]').forEach((el) => {
      el.addEventListener('change', () => {
        const it = selItem(); if (it) setAdj(it, el.dataset.adjt, el.checked);
      });
    });
    document.querySelectorAll('[data-adjn]').forEach((el) => {
      el.addEventListener('input', () => {
        const it = selItem();
        if (it) setAdj(it, el.dataset.adjn, num(el.value, ADJ_DEF[el.dataset.adjn]));
      });
    });
    document.querySelectorAll('[data-adjc]').forEach((el) => {
      el.addEventListener('input', () => {
        const it = selItem(); if (!it) return;
        const k = el.dataset.adjc;
        if (k === 'tint') {
          const t = Object.assign({ amount: 0.3 }, adjOf(it).tint || {});
          t.color = el.value;
          setAdj(it, 'tint', t);
        } else {
          setAdj(it, k, el.value);
        }
      });
    });
    document.querySelectorAll('[data-adjc-on]').forEach((el) => {
      el.addEventListener('change', () => {
        const it = selItem(); if (!it) return;
        const k = el.dataset.adjcOn;
        const row = el.closest('.row');
        const colorEl = row ? row.querySelector('[data-adjc]') : null;
        if (!el.checked) { setAdj(it, k, null); refreshAdjUI(it); return; }
        if (k === 'tint') {
          const t = Object.assign({ amount: 0.3 }, adjOf(it).tint || {});
          if (colorEl) t.color = colorEl.value;
          setAdj(it, 'tint', t);
        } else {
          setAdj(it, k, (colorEl && colorEl.value) || '#ffffff');
        }
        refreshAdjUI(it);
      });
    });

    $('btnLoadColors').addEventListener('click', () => {
      const it = selItem();
      if (!it) { setStatus('先选中一个面板', 'err'); return; }
      loadColors(it);
    });
    $('btnAddRemap').addEventListener('click', () => {
      const it = selItem();
      if (!it) { setStatus('先选中一个面板', 'err'); return; }
      const rows = (adjOf(it).remap || []).slice();
      rows.push({ from: '#ffffff', to: '#ffffff', tol: 6 });
      setAdj(it, 'remap', rows);
      renderRemap(it);
    });
    // ════════════════════════════════════════════════════════
    // SVG 矢量编辑（仅 .svg 面板可用；编辑发生在 XML 层，不是像素层）
    // ════════════════════════════════════════════════════════
    $('btnSvgColors').addEventListener('click', () => svgLoad('colors'));
    $('btnSvgTexts').addEventListener('click', () => svgLoad('texts'));
    $('btnSvgLegend').addEventListener('click', () => {
      const it = selItem();
      if (!it) { setStatus('先选中一个 SVG 面板', 'err'); return; }
      setSvg(it, {
        legend: {
          action: $('svgLegendAction').value,
          dx: num($('svgLegendDx').value, 0),
          dy: num($('svgLegendDy').value, 0),
          font_scale: num($('svgLegendFs').value, 1),
        },
      });
      setStatus('已应用到图例', 'ok');
    });
    $('btnSvgFontApply').addEventListener('click', () => {
      const it = selItem();
      if (!it) { setStatus('先选中一个 SVG 面板', 'err'); return; }
      const fs = num($('svgFontScale').value, 1);
      const sw = num($('svgStrokeScale').value, 1);
      const fam = $('svgFontFamily').value.trim();
      const tc = $('svgTextColor').value;
      setSvg(it, {
        font_scale: fs !== 1 ? fs : null,
        stroke_scale: sw !== 1 ? sw : null,
        font_family: fam || null,
        text_color: tc && tc !== '#222222' ? tc : null,
      });
      setStatus('已应用文字/线宽', 'ok');
    });
    $('btnSvgFontReset').addEventListener('click', () => {
      $('svgFontScale').value = 1;
      $('svgStrokeScale').value = 1;
      $('svgFontFamily').value = '';
      $('svgTextColor').value = '#222222';
      const it = selItem();
      if (it) setSvg(it, { font_scale: null, stroke_scale: null,
                           font_family: null, text_color: null });
    });
    $('btnSvgHide').addEventListener('click', () => {
      const it = selItem();
      if (!it) { setStatus('先选中一个 SVG 面板', 'err'); return; }
      const sel = $('svgHideSel').value.trim();
      if (!sel) { setStatus('填一个选择器，例如 #text_12', 'err'); return; }
      const cur = (svgOf(it).hide || []).slice();
      if (cur.indexOf(sel) < 0) cur.push(sel);
      setSvg(it, { hide: cur });
      $('svgHideSel').value = '';
      setStatus('已隐藏 ' + sel, 'ok');
    });
    $('btnSvgReset').addEventListener('click', () => {
      const it = selItem();
      if (!it) return;
      const a = Object.assign({}, it.adjust || {});
      delete a.svg;
      if (Object.keys(a).length) it.adjust = a; else delete it.adjust;
      schedulePreview(it);
      refreshSvgUI(it);
      setStatus('已重置该面板矢量编辑', 'ok');
    });

    $('btnAdjApplyAll').addEventListener('click', applyAdjToAll);
    $('btnAdjReset').addEventListener('click', () => {
      const it = selItem();
      if (!it) { setStatus('先选中一个面板', 'err'); return; }
      delete it.adjust;
      schedulePreview(it);
      refreshAdjUI(it);
      setStatus('已重置该面板调色', 'ok');
    });

    // 未选中时，调色控件全是 disabled（点了没反应也不提示）→ 用遮罩接住点击，
    // 给明确指引，并闪一下画布上的面板，告诉用户"点这里"。
    $('adjCover').addEventListener('click', () => {
      setStatus('请先在画布上点选一个面板图（出现蓝框即选中）', 'err');
      document.querySelectorAll('.item').forEach((el) => {
        el.classList.remove('flash');
        void el.offsetWidth;               // 重启动画
        el.classList.add('flash');
        setTimeout(() => el.classList.remove('flash'), 1100);
      });
    });

    // 记录最近聚焦的颜色框（调色板点击时填这里）
    document.addEventListener('focusin', (e) => {
      if (e.target && e.target.type === 'color') lastColorInput = e.target;
    });
  }

  // ════════════════════════════════════════════════════════════
  // 导出 / 布局存取
  // ════════════════════════════════════════════════════════════
  function toLayout() {
    return {
      canvas: Object.assign({}, state.canvas),
      style: Object.assign({}, state.style),
      title: Object.assign({}, state.title),
      items: state.items.map((i) => ({
        id: i.id, src: i.src, x_in: i.x_in, y_in: i.y_in,
        w_in: i.w_in, h_in: i.h_in, letter: i.letter,
        caption: i.caption || '', z: num(i.z, 0),
        adjust: cleanAdjust(i.adjust) || undefined,
      })),
    };
  }

  async function doExport() {
    if (!state.items.length) { setStatus('画布是空的', 'err'); return; }
    const stem = $('stemInput').value.trim() || 'figure';
    setStatus('渲染中…');
    try {
      const r = await fetch('/api/render', {
        method: 'POST', headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ layout: toLayout(), stem, formats: ['png', 'pdf'], dpi: state.canvas.dpi }),
      });
      const j = await r.json();
      if (j.error) { setStatus('导出失败：' + j.error, 'err'); return; }
      const warn = (j.warnings || []).length ? ` · ${j.warnings.length} 条提示` : '';
      setStatus(`导出成功：${j.stem}（${j.size_px[0]}×${j.size_px[1]} px）${warn}`, 'ok');
      if ((j.warnings || []).length) console.warn('[Figure Studio]', j.warnings);
      loadOutputs();
      openModal(`<h2>导出完成</h2>
        <div class="list">
          ${j.urls.png ? `<a href="${j.urls.png}" target="_blank">PNG · ${j.size_px[0]}×${j.size_px[1]} px</a>` : ''}
          ${j.urls.pdf ? `<a href="${j.urls.pdf}" target="_blank">PDF · 矢量（pdf.fonttype=42）</a>` : ''}
        </div>
        ${(j.warnings || []).length ? '<div class="warn">提示：<br>' + j.warnings.map(escapeHtml).join('<br>') + '</div>' : ''}
        <div class="foot"><button onclick="document.getElementById('modalMask').classList.remove('on')">关闭</button></div>`);
    } catch (e) {
      setStatus('导出失败：' + e.message, 'err');
    }
  }

  async function loadOutputs() {
    try {
      const j = await (await fetch('/api/outputs')).json();
      $('outList').innerHTML = (j.items || []).slice(0, 40).map((o) =>
        `<a href="${o.url}" target="_blank"><span>${escapeHtml(o.name)}</span><span class="meta">${fmtSize(o.size_kb * 1024)}</span></a>`
      ).join('') || '<div class="meta">还没有导出产物</div>';
    } catch (e) { /* 忽略 */ }
  }

  function openModal(html) { $('modalBox').innerHTML = html; $('modalMask').classList.add('on'); }
  function closeModal() { $('modalMask').classList.remove('on'); }

  // ════════════════════════════════════════════════════════════
  // svguest 精修（交互式矢量编辑器，零改动移植）
  //   单面板 / 全部 SVG 面板（多标签）/ 写回 / 导出直通 / 备份回滚
  // ════════════════════════════════════════════════════════════
  const sgSvgItems = () => state.items.filter((i) => isSvg(i));

  // 内嵌浮层：同源 iframe，不弹独立窗口
  //   ⚠ 早期用 window.open 是把它当独立应用 —— 它确实自带全局快捷键和
  //     菜单，但弹窗会被浏览器拦、存完还得切回原窗口。同源 iframe 同样
  //     能跑（快捷键、菜单、视口都在 iframe 内自成一体），而且关掉即销毁。
  // 三档布局：full 全屏 / split 右侧分屏（工作台仍可点）/ float 浮动窗口
  const SG_MODES = ['full', 'split', 'float'];
  const SG_MODE_CN = { full: '全屏', split: '右侧分屏', float: '浮动窗口' };

  function sgMode() {
    const l = $('sgLayer');
    const m = l && l.dataset.mode;
    return SG_MODES.indexOf(m) >= 0 ? m : 'full';
  }

  function setSgMode(m) {
    const l = $('sgLayer');
    if (!l) return;
    if (SG_MODES.indexOf(m) < 0) m = 'full';
    l.dataset.mode = m;
    /* ⚠ float 模式靠 inline 定位与缩放（拖标题栏写 left/top/right/bottom，
     *  resize:both 写 width/height）。切走时必须**六个全清**，否则残留的
     *  inline 值会盖住 CSS 的 inset —— 漏掉 right/bottom 时元素会从
     *  「撑满屏幕」退化成 shrink-to-fit（实测 1680px → 444px），
     *  而且**零报错**，只是"看着不对"。 */
    l.style.left = ''; l.style.top = '';
    l.style.right = ''; l.style.bottom = '';
    l.style.width = ''; l.style.height = '';
    const b = $('btnSgMode');
    if (b) b.textContent = '布局：' + SG_MODE_CN[m];
    try { localStorage.setItem('fs.sgMode', m); } catch (e) { /* 隐私模式忽略 */ }
  }

  function cycleSgMode() {
    const i = SG_MODES.indexOf(sgMode());
    setSgMode(SG_MODES[(i + 1) % SG_MODES.length]);
    setStatus('精修面板布局：' + SG_MODE_CN[sgMode()] +
      (sgMode() === 'split' ? '（左侧工作台仍可操作）' : ''), 'ok');
  }

  function bindSgLayer() {
    const l = $('sgLayer'); const head = $('sgHead'); const b = $('btnSgMode');
    if (b) b.addEventListener('click', cycleSgMode);
    if (!head || !l) return;

    // 双击标题栏：全屏 ↔ 右侧分屏（最常用的两档来回切）
    head.addEventListener('dblclick', (e) => {
      if (e.target.closest('button')) return;
      setSgMode(sgMode() === 'full' ? 'split' : 'full');
    });

    // float 模式：拖标题栏移动窗口
    head.addEventListener('pointerdown', (e) => {
      if (sgMode() !== 'float') return;
      if (e.target.closest('button')) return;
      e.preventDefault();
      const r = l.getBoundingClientRect();
      const sx = e.clientX; const sy = e.clientY;
      const ox = r.left; const oy = r.top;
      function mv(ev) {
        const nx = Math.max(0, Math.min(window.innerWidth - 140, ox + ev.clientX - sx));
        const ny = Math.max(0, Math.min(window.innerHeight - 44, oy + ev.clientY - sy));
        l.style.left = nx + 'px';
        l.style.top = ny + 'px';
        l.style.right = 'auto';
        l.style.bottom = 'auto';
      }
      function up() {
        window.removeEventListener('pointermove', mv);
        window.removeEventListener('pointerup', up);
      }
      window.addEventListener('pointermove', mv);
      window.addEventListener('pointerup', up);
    });
  }

  function closeSvguest() {
    const l = $('sgLayer');
    if (!l || !l.classList.contains('on')) return;
    l.classList.remove('on');
    const f = $('sgFrame');
    // 回 about:blank 而不是只隐藏：让上游的全局监听、IndexedDB 连接、
    // 定时器随文档一起卸载，不然每开一次就多留一份。
    if (f) f.src = 'about:blank';
    refreshSgHistory(true);
  }

  async function openSvguest(items, opts) {
    const list = (items || []).filter(isSvg);
    if (!list.length) {
      setStatus('没有 SVG 面板可精修（PNG/JPG 请用上面的调色）', 'err');
      return;
    }
    const layer = $('sgLayer');
    const frame = $('sgFrame');
    if (!layer || !frame) {
      setStatus('页面缺少精修浮层（static/index.html 没更新？）', 'err');
      return;
    }
    setStatus('正在准备精修会话…');
    let j;
    try {
      const r = await fetch('/api/svg/session', {
        method: 'POST', headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({
          items: list.map((i) => ({ src: i.src, name: i.name })),
        }),
      });
      j = await r.json();
    } catch (e) {
      setStatus('准备会话失败：' + e.message, 'err');
      return;
    }
    if (j.error) { setStatus('准备会话失败：' + j.error, 'err'); return; }
    // ⚠ 平台目前**没有**深色主题（style.css 里没有 [data-theme] 规则，也没有
    //   切换入口），所以这里实际恒为 'light'。保留读取是为了将来加主题时
    //   不用再改这一行。
    const theme = document.documentElement.dataset.theme === 'dark' ? 'dark' : 'light';
    const url = '/static/vendor/svguest.html#s=' + encodeURIComponent(j.token) +
      '&theme=' + theme + ((opts && opts.chat) ? '&chat=1' : '');
    const sub = $('sgLayerSub');
    if (sub) {
      sub.textContent = list.length > 1
        ? `${list.length} 个面板（多标签）· 改完按 Ctrl+S 全部写回`
        : `${list[0].name} · 改完按 Ctrl+S 写回`;
    }
    layer.classList.add('on');
    // 恢复上次用的布局（全屏/分屏/浮动）
    let want = 'full';
    try { want = localStorage.getItem('fs.sgMode') || 'full'; } catch (e) { /* 忽略 */ }
    setSgMode(want);
    // 每次会话 token 都不同 → 赋 src 必然触发一次全新导航，不需要额外的
    // cache-busting 参数，也不会复用上一次的脏状态。
    // ⚠ 载入后必须把焦点推进 iframe：点的是平台上的按钮，焦点本来留在平台，
    //   Ctrl+S 和 m/v 这些工具快捷键就一个都到不了 svguest（零报错，像坏了）。
    frame.addEventListener('load', function focusOnce() {
      frame.removeEventListener('load', focusOnce);
      try { frame.contentWindow.focus(); } catch (e) { /* 跨域或未就绪就算了 */ }
    });
    frame.src = url;
    setStatus(list.length > 1
      ? `已在内嵌编辑器中打开 ${list.length} 个面板（多标签），改完按 Ctrl+S 全部保存`
      : `已在内嵌编辑器中打开 ${list[0].name}，改完按 Ctrl+S 写回`, 'ok');
  }

  function openSelectedInSvguest() {
    const it = selItem();
    if (!it) { setStatus('先选中一个 SVG 面板（画布上点一下出现蓝框）', 'err'); return; }
    if (!isSvg(it)) {
      setStatus('svguest 只能精修 SVG 面板（当前是 ' +
        (it.src.split('.').pop() || '?').toUpperCase() + '）', 'err');
      return;
    }
    openSvguest([it]);
  }

  // ── 备份历史 ────────────────────────────────────────────────
  // refreshInspector 在拖动时会被高频调用，这里做防抖 + 同源去重，
  // 免得每挪一下鼠标就打一次 /api/svg/history
  let sgHistTimer = null;
  let sgHistKey = '';
  function refreshSgHistory(force) {
    if (force) sgHistKey = '';
    clearTimeout(sgHistTimer);
    sgHistTimer = setTimeout(() => {
      doRefreshSgHistory();
      // 快照列表跟着一起刷（它自己的请求更便宜，也做同源去重）
      if (typeof doRefreshAsList === 'function') doRefreshAsList(force);
    }, force ? 0 : 160);
  }

  async function doRefreshSgHistory() {
    const box = $('sgHistory');
    const tag = $('sgTag');
    if (!box) return;
    const it = selItem();
    if (!it || !isSvg(it)) {
      sgHistKey = '';
      if (tag) { tag.textContent = '未选中 SVG'; tag.className = 'tag'; }
      box.innerHTML = '<div class="hint">选中一个 SVG 面板后显示</div>';
      return;
    }
    const key = it.src + '|' + (it.adjust ? JSON.stringify(it.adjust) : '');
    if (key === sgHistKey && box.dataset.ready === '1') return;
    let j;
    try {
      const r = await fetch('/api/svg/history?src=' + encodeURIComponent(it.src));
      j = await r.json();
    } catch (e) {
      box.innerHTML = '<div class="hint">读历史失败：' + escapeHtml(e.message) + '</div>';
      return;
    }
    sgHistKey = key;
    box.dataset.ready = '1';
    const bs = (j && j.backups) || [];
    if (tag) {
      tag.textContent = bs.length ? bs.length + ' 份备份' : '无备份';
      tag.className = 'tag' + (bs.length ? ' ok' : '');
    }
    if (!bs.length) {
      box.innerHTML = '<div class="hint">这个面板还没被精修过</div>';
      return;
    }
    box.innerHTML = bs.slice(0, 8).map((b) => (
      '<div class="sg-row">' +
      '<span class="sg-t">' + escapeHtml(b.ts) + '</span>' +
      '<span class="sg-s">' + (b.bytes / 1024).toFixed(0) + ' KB</span>' +
      '<button class="ghost sg-restore" data-bak="' + escapeHtml(b.file) + '"' +
      ' title="把这个备份写回原文件（当前版本会先被再备份一次，所以回滚本身也能撤销）">回滚</button>' +
      '</div>'
    )).join('');
    box.querySelectorAll('.sg-restore').forEach((btn) => {
      btn.addEventListener('click', () => restoreBackup(it.src, btn.dataset.bak, btn));
    });
  }

  async function restoreBackup(src, backup, btn) {
    if (!confirm('把 ' + backup + ' 的内容写回原文件？\n\n当前版本会先自动备份，回滚本身也可撤销。')) return;
    btn.disabled = true;
    btn.textContent = '回滚中';
    try {
      const r = await fetch('/api/svg/restore', {
        method: 'POST', headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ src: src, backup: backup }),
      });
      const j = await r.json();
      if (j.error) throw new Error(j.error);
      const it = state.items.find((x) => x.src === src);
      if (it) {
        const el = document.querySelector('.item[data-id="' + it.id + '"] img');
        if (el) el.src = imgUrl(it.src, it.adjust) + '&_t=' + Date.now();
      }
      setStatus('已回滚到 ' + backup, 'ok');
      refreshSgHistory(true);
    } catch (e) {
      setStatus('回滚失败：' + e.message, 'err');
      btn.disabled = false;
      btn.textContent = '回滚';
    }
  }

  // ── 磁盘自动快照（svguest 编辑停顿 15s 后服务端落盘的那份）──────
  // ⚠ 和上面的「备份历史」是两回事：备份只在**写回原文件**时产生（旧版副本），
  //   自动快照在**编辑过程中**产生，所以它能救「改到一半崩了、还没保存」。
  let asKey = '';

  async function doRefreshAsList(force) {
    const box = $('asList');
    const tag = $('asTag');
    if (!box) return;
    const it = selItem();
    if (!it || !isSvg(it)) {
      asKey = '';
      if (tag) { tag.textContent = ''; tag.className = 'tag'; }
      box.innerHTML = '<div class="hint">选中一个 SVG 面板后显示</div>';
      return;
    }
    const key = it.src;
    if (!force && key === asKey && box.dataset.ready === '1') return;
    let j;
    try {
      const r = await fetch('/api/svg/autosave?src=' + encodeURIComponent(it.src));
      j = await r.json();
    } catch (e) {
      box.innerHTML = '<div class="hint">读快照失败：' + escapeHtml(e.message) + '</div>';
      return;
    }
    asKey = key;
    box.dataset.ready = '1';
    const ss = (j && j.snapshots) || [];
    if (tag) {
      tag.textContent = ss.length ? ss.length + ' 份' : '无';
      tag.className = 'tag' + (ss.length ? ' ok' : '');
    }
    if (!ss.length) {
      box.innerHTML = '<div class="hint">还没有自动快照 —— 在 svguest 里改几笔、停 15 秒就会自动生成</div>';
      return;
    }
    box.innerHTML = ss.slice(0, 8).map((b) => (
      '<div class="sg-row">' +
      '<span class="sg-t">' + escapeHtml(b.ts) + '</span>' +
      '<span class="sg-s">' + (b.bytes / 1024).toFixed(0) + ' KB</span>' +
      '<button class="ghost as-restore" data-snap="' + escapeHtml(b.file) + '"' +
      ' title="把这份快照写回原文件（写回前仍会自动备份当前版本）">恢复</button>' +
      '</div>'
    )).join('');
    box.querySelectorAll('.as-restore').forEach((btn) => {
      btn.addEventListener('click', () => restoreSnapshot(it.src, btn.dataset.snap, btn));
    });
  }

  async function restoreSnapshot(src, file, btn) {
    if (!confirm('把自动快照 ' + file + ' 写回原文件？\n\n当前版本会先自动备份，所以这一步也能撤销。')) return;
    btn.disabled = true;
    btn.textContent = '恢复中';
    try {
      const r = await fetch('/api/svg/autosave_restore', {
        method: 'POST', headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ src: src, file: file }),
      });
      const j = await r.json();
      if (j.error) throw new Error(j.error.message || j.error);
      refreshPanelPreview(src);
      setStatus('已从自动快照恢复：' + file, 'ok');
      refreshSgHistory(true);
    } catch (e) {
      setStatus('恢复失败：' + e.message, 'err');
      btn.disabled = false;
      btn.textContent = '恢复';
    }
  }

  // ── AI 绘图（SVG Chat）密钥 ─────────────────────────────────
  let aiCfg = null;

  async function loadAiConfig() {
    const tag = $('aiTag');
    const hint = $('aiHint');
    try {
      aiCfg = await (await fetch('/api/ai/config')).json();
    } catch (e) {
      aiCfg = null;
    }
    if (!aiCfg || !aiCfg.providers) {
      if (tag) { tag.textContent = '读不到'; tag.className = 'tag'; }
      if (hint) hint.textContent = '服务端没有响应 /api/ai/config';
      return;
    }
    const P = aiCfg.providers;
    const on = ['Ds', 'Glm'].filter((k) => P[k] && P[k].configured);
    if (tag) {
      tag.textContent = on.length ? '已配 ' + on.map((k) => P[k].name).join(' + ') : '未配置';
      tag.className = 'tag' + (on.length ? ' ok' : '');
    }
    if (hint) {
      hint.innerHTML = ['Ds', 'Glm'].map((k) => {
        const p = P[k] || {};
        return '<div>' + escapeHtml(p.name || k) + '：' +
          (p.configured ? '已存 <code>' + escapeHtml(p.hint) + '</code>' : '未配置') + '</div>';
      }).join('') + '<div style="opacity:.65">密钥文件：<code>' +
        escapeHtml(aiCfg.file || '') + '</code></div>';
    }
    // 输入框留空 = 不改；已配置的用占位符提示掩码
    const ds = $('aiKeyDs');
    const gl = $('aiKeyGlm');
    if (ds && !ds.value) ds.placeholder = P.Ds.configured ? P.Ds.hint + '（留空不改）' : 'sk-…';
    if (gl && !gl.value) gl.placeholder = P.Glm.configured ? P.Glm.hint + '（留空不改）' : '…';
  }

  async function saveAiKeys() {
    const ds = (($('aiKeyDs') || {}).value || '').trim();
    const gl = (($('aiKeyGlm') || {}).value || '').trim();
    if (!ds && !gl) { setStatus('两个输入框都空着，没有要保存的内容', 'err'); return; }
    const jobs = [];
    if (ds) jobs.push(['Ds', ds]);
    if (gl) jobs.push(['Glm', gl]);
    const btn = $('btnAiSave');
    if (btn) btn.disabled = true;
    try {
      for (const [p, k] of jobs) {
        const r = await fetch('/api/ai/config', {
          method: 'POST', headers: { 'Content-Type': 'application/json' },
          body: JSON.stringify({ provider: p, key: k }),
        });
        const j = await r.json();
        if (j.error) throw new Error(j.error.message || j.error);
      }
      if ($('aiKeyDs')) $('aiKeyDs').value = '';
      if ($('aiKeyGlm')) $('aiKeyGlm').value = '';
      await loadAiConfig();
      setStatus('AI 密钥已存到本机（不写进浏览器、不随布局文件走）', 'ok');
    } catch (e) {
      setStatus('保存密钥失败：' + e.message, 'err');
    } finally {
      if (btn) btn.disabled = false;
    }
  }

  // svguest 的消息：写回 / 批量写回 / 导出
  window.addEventListener('message', (e) => {
    if (e.origin !== location.origin) return;
    const d = e.data || {};
    if (d.fs === 'svguest-saved') {
      refreshPanelPreview(d.src);
      setStatus('svguest 的修改已写回并刷新：' + shortName(d.src), 'ok');
      refreshSgHistory(true);
    } else if (d.fs === 'svguest-saved-batch') {
      (d.srcs || []).forEach(refreshPanelPreview);
      setStatus('svguest 已写回 ' + ((d.srcs || []).length) + ' 个面板', 'ok');
      refreshSgHistory(true);
    } else if (d.fs === 'svguest-exported') {
      onSvguestExported(d);
    }
  });

  function shortName(p) {
    const s = String(p || '');
    return s.split(/[\\/]/).pop() || s;
  }

  function refreshPanelPreview(src) {
    const it = state.items.find((x) => x.src === src);
    if (!it) return;
    const el = document.querySelector('.item[data-id="' + it.id + '"] img');
    if (el) el.src = imgUrl(it.src, it.adjust) + '&_t=' + Date.now();
  }

  // 导出的图落在原图同目录 → 若在当前面板库根目录下，顺手重扫一次，
  // 让导出的 PNG 立刻能当面板拖进画布
  async function onSvguestExported(d) {
    const p = String(d.path || '');
    const root = ($('rootInput').value || '').trim();
    const under = root && p.toLowerCase().startsWith(
      root.replace(/[\\/]+$/, '').toLowerCase());
    setStatus('已导出：' + p + (under ? '（面板库已刷新）' : ''), 'ok');
    if (under) {
      try { await scan(); } catch (e) { /* 刷新失败不影响导出本身 */ }
    }
  }

  async function saveLayout() {

    const name = prompt('布局名称：', ($('stemInput').value || 'layout').replace(/[^\w\u4e00-\u9fff\-]/g, '_'));
    if (!name) return;
    const r = await fetch('/api/layouts/' + encodeURIComponent(name), {
      method: 'POST', headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify(toLayout()),
    });
    const j = await r.json();
    setStatus(j.ok ? '布局已保存：' + j.name : '保存失败', j.ok ? 'ok' : 'err');
  }

  async function openLayout() {
    const j = await (await fetch('/api/layouts')).json();
    const rows = (j.items || []).map((o) =>
      `<button data-name="${escapeHtml(o.name)}"><span>${escapeHtml(o.name)}</span><span class="meta">${o.items} 面板 · ${fmtTime(o.mtime)}</span></button>`
    ).join('') || '<div class="meta">还没有保存过布局</div>';
    openModal(`<h2>打开布局</h2><div class="list">${rows}</div>
      <div class="foot"><button onclick="document.getElementById('modalMask').classList.remove('on')">取消</button></div>`);
    $('modalBox').querySelectorAll('button[data-name]').forEach((b) => {
      b.addEventListener('click', async () => {
        const d = await (await fetch('/api/layouts/' + encodeURIComponent(b.dataset.name))).json();
        applyLayout(d);
        closeModal();
        setStatus('已加载布局：' + b.dataset.name, 'ok');
      });
    });
  }

  function applyLayout(d) {
    const base = DEFAULTS();
    state = Object.assign(base, {
      // ⚠ selSet 必须显式带上：这个 Object.assign 是**整体换掉 state**，
      //   漏了它后面所有 `state.selSet.has(...)` 都会抛 TypeError（零提示崩溃）。
      sel: null, selSet: new Set(),
      zoom: state.zoom, snap: state.snap, showGrid: state.showGrid,
    });
    if (d.canvas) Object.assign(state.canvas, d.canvas);
    if (d.style) Object.assign(state.style, d.style);
    if (d.title) Object.assign(state.title, d.title);
    state.items = (d.items || []).map((i) => Object.assign({ aspect: 0.6 }, i));
    dom.forEach((el) => el.remove()); dom.clear();
    uid = state.items.length + 1;
    // 补齐 aspect / 自然尺寸（老布局可能缺）
    state.items.forEach((i) => { if (!i.aspect && i.px_h && i.px_w) i.aspect = i.px_h / i.px_w; });
    /* ⚠ 旧版 resize 的吸附 bug 会把 w_in/h_in 写成 0（面板"凭空消失"）。
     * 这种布局一旦保存过，再打开就永远看不见 —— 加载时兜底修复，
     * 用加入时的自然宽度（w0_in），拿不到就给一个安全默认值。 */
    const _m = state.style.margin_in;
    state.items.forEach((i) => {
      if (!isFinite(i.w_in) || i.w_in <= 0.01) {
        i.w_in = r4((isFinite(i.w0_in) && i.w0_in > 0.01)
          ? i.w0_in
          : Math.min(4, Math.max(0.5, state.canvas.width_in - 2 * _m)));
      }
      if (!isFinite(i.h_in) || i.h_in <= 0.01) {
        i.h_in = r4(i.w_in * (i.aspect || 0.6));
      }
      if (!isFinite(i.x_in)) i.x_in = r4(_m);
      if (!isFinite(i.y_in)) i.y_in = r4(_m);
    });
    // ⚠ 历史上加载布局后 sel=null → 调色控件全 disabled、点击零反馈，用户以为"没动静"。
    //   加载布局的目的就是接着调 → 自动选中第一个面板。
    if (state.items.length && !state.items.some((i) => i.id === state.sel)) {
      state.sel = state.items[0].id;
    }
    if (state.sel) state.selSet = new Set([state.sel]);
    refreshInspector(); renderCanvas();
    pushHist();          // 打开布局本身也是一步（可撤销回打开前的状态）
  }

  // ════════════════════════════════════════════════════════════
  // 事件绑定
  // ════════════════════════════════════════════════════════════
  function bindUI() {
    bindAdj();
    $('btnScan').addEventListener('click', scan);
    $('rootInput').addEventListener('keydown', (e) => { if (e.key === 'Enter') scan(); });
    $('searchInput').addEventListener('keydown', (e) => { if (e.key === 'Enter') scan(); });
    $('btnClearSearch').addEventListener('click', () => { $('searchInput').value = ''; scan(); });
    $('chkOnlyNew').addEventListener('change', () => { /* 提示性开关，扫描时统一 600 上限 */ });

    $('btnExport').addEventListener('click', doExport);
    $('btnExport2').addEventListener('click', doExport);
    $('btnSaveLayout').addEventListener('click', saveLayout);
    $('btnSaveLayout2').addEventListener('click', saveLayout);
    $('btnSvguest').addEventListener('click', openSelectedInSvguest);
    $('btnSgOne').addEventListener('click', openSelectedInSvguest);
    $('btnSgAll').addEventListener('click', () => openSvguest(sgSvgItems()));
    $('btnAiSave').addEventListener('click', saveAiKeys);
    $('btnAiChat').addEventListener('click', () => {
      const it = selItem();
      if (!it || !isSvg(it)) {
        setStatus('先在画布上选一个 SVG 面板，AI 改图要有对象', 'err');
        return;
      }
      openSvguest([it], { chat: true });
    });
    $('btnOpenLayout').addEventListener('click', openLayout);

    // 精修浮层的退出：顶栏关闭按钮 + 平台侧 Esc。
    // ⚠ Esc 只在**焦点不在 iframe 里**时生效 —— 焦点在 iframe 里时 Esc 属于
    //   上游（取消工具/取消选择），抢过来会让用户「想取消选择」却关掉整个
    //   编辑器。这个监听自然只在焦点落在平台时才收到，正好是我们要的语义。
    $('btnSgClose').addEventListener('click', closeSvguest);
    bindSgLayer();
    document.addEventListener('keydown', (e) => {
      if (e.key === 'Escape' && $('sgLayer').classList.contains('on')) closeSvguest();
    });

    $('btnZoomIn').addEventListener('click', () => { state.zoom = Math.min(4, state.zoom * 1.15); renderCanvas(); });
    $('btnZoomOut').addEventListener('click', () => { state.zoom = Math.max(0.15, state.zoom / 1.15); renderCanvas(); });
    $('btnFit').addEventListener('click', () => { state.zoom = 1; renderCanvas(); });
    $('chkSnap').addEventListener('change', (e) => { state.snap = e.target.checked; });
    $('chkGrid').addEventListener('change', (e) => { state.showGrid = e.target.checked; renderCanvas(); });

    document.querySelectorAll('[data-tool]').forEach((b) => {
      b.addEventListener('click', () => {
        const t = b.dataset.tool;
        if (t === 'autoletter') autoLetter();
        else if (t === 'trim') trimCanvas();
        else toolArrange(t);
      });
    });
    document.querySelectorAll('[data-align]').forEach((b) => {
      b.addEventListener('click', () => align(b.dataset.align));
    });
    document.querySelectorAll('[data-item]').forEach((b) => {
      b.addEventListener('click', () => itemAction(b.dataset.item));
    });

    $('modalMask').addEventListener('click', (e) => { if (e.target.id === 'modalMask') closeModal(); });

    // 撤销 / 重做
    if ($('btnUndo')) $('btnUndo').addEventListener('click', undo);
    if ($('btnRedo')) $('btnRedo').addEventListener('click', redo);

    // 画布空白点击 → 取消选中（同时清空多选集合）
    scrollEl.addEventListener('pointerdown', (e) => {
      if (e.target === canvasEl || e.target === guidesEl || e.target === scrollEl || e.target.classList.contains('canvas-stage')) {
        select(null);
      }
    });

    // 拖放落位
    canvasEl.addEventListener('dragover', (e) => { e.preventDefault(); e.dataTransfer.dropEffect = 'copy'; });
    canvasEl.addEventListener('drop', async (e) => {
      e.preventDefault();
      const src = e.dataTransfer.getData('text/plain');
      const p = libPanels.find((x) => x.src === src);
      if (!p) return;
      const rect = canvasEl.getBoundingClientRect();
      addPanel(p, { x: Math.max(0, (e.clientX - rect.left) / ppi), y: Math.max(0, (e.clientY - rect.top) / ppi) });
    });

    // 键盘
    //   ⚠ 顺序有意义：撤销/重做/全选/复制 这些**带修饰键**的先判，再判
    //     方向键/Delete 这类裸键。否则在输入框里打字会被下面的裸键逻辑吃掉
    //     （下面统一有 input/textarea 早退，这里也一样）。
    window.addEventListener('keydown', (e) => {
      const tag = (e.target.tagName || '').toLowerCase();
      const typing = tag === 'input' || tag === 'textarea' || e.target.isContentEditable;
      const mod = e.ctrlKey || e.metaKey;

      if (mod && !e.altKey) {
        const k = e.key.toLowerCase();
        if (k === 'z') {
          if (typing) return;                 // 输入框内保留浏览器原生撤销
          e.preventDefault();
          e.shiftKey ? redo() : undo();
          return;
        }
        if (k === 'y') { if (typing) return; e.preventDefault(); redo(); return; }
        if (k === 'a') {                      // 全选面板（不是全选文字）
          if (typing) return;
          e.preventDefault();
          if (!state.items.length) { setStatus('画布是空的', 'err'); return; }
          state.selSet = new Set(state.items.map((i) => i.id));
          state.sel = state.items[state.items.length - 1].id;
          refreshInspector(); renderCanvas();
          setStatus(`已全选 ${state.selSet.size} 个面板`, 'ok');
          return;
        }
        if (k === 'd') {                      // 复制（多选时整组复制）
          if (typing) return;
          e.preventDefault();
          if (state.selSet.size || state.sel) itemAction('dup');
          return;
        }
      }

      if (typing) return;
      const it = selItem();
      if (!it) return;
      const step = e.shiftKey ? 0.1 : 0.01;
      if (e.key === 'Delete' || e.key === 'Backspace') { e.preventDefault(); itemAction('del'); return; }
      // 方向键：多选时整组一起挪
      let dx = 0; let dy = 0;
      if (e.key === 'ArrowLeft') dx = -step;
      else if (e.key === 'ArrowRight') dx = step;
      else if (e.key === 'ArrowUp') dy = -step;
      else if (e.key === 'ArrowDown') dy = step;
      else return;
      e.preventDefault();
      const mv = (state.selSet.size > 1)
        ? state.items.filter((i) => state.selSet.has(i.id))
        : [it];
      for (const i of mv) { i.x_in = r4(i.x_in + dx); i.y_in = r4(i.y_in + dy); }
      renderCanvas(); refreshInspector();
      markDirty();                            // 连按合并成一步撤销
    });

    // Ctrl/⌘+E：把选中的 SVG 面板交给 svguest 精修
    window.addEventListener('keydown', (e) => {
      if (!(e.key === 'e' || e.key === 'E')) return;
      if (!(e.ctrlKey || e.metaKey) || e.shiftKey || e.altKey) return;
      const tag = (e.target.tagName || '').toLowerCase();
      if (tag === 'input' || tag === 'textarea') return;
      if (!selItem()) return;
      e.preventDefault();
      openSelectedInSvguest();
    });

    let rt = null;
    window.addEventListener('resize', () => { clearTimeout(rt); rt = setTimeout(renderCanvas, 120); });
  }

  // ════════════════════════════════════════════════════════════
  // 启动
  // ════════════════════════════════════════════════════════════
  async function init() {
    bindInspector();
    bindUI();
    await loadRoots();
    loadPalette();
    refreshAdjUI(null);
    renderCanvas();
    // ⚠ 初次加载也要同步一次检查器：否则空画布时「精修/全部精修」等按钮
    //   保持 DOM 默认的 enabled，点下去才发现无事可做（零反馈的又一种）
    refreshInspector();
    loadOutputs();
    loadAiConfig();
    // 撤销栈的基线：空布局记一格，这样第一次改动也能退回来
    pushHist();
    if ($('rootInput').value) scan();
    setStatus('就绪');
  }

  document.addEventListener('DOMContentLoaded', init);
})();
