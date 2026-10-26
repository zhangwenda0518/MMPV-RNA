"""table_widgets.py — MMPV 全管线 HTML 报告表格统一控件 (零依赖, 纯内联 CSS/JS)

给任意静态报告里的数据表自动挂四件套:
  1. 筛选   — 工具条 Filter 输入框, 逐行匹配
  2. 排序   — 点击表头切换升/降序 (数值感知), 箭头指示
  3. 列可选 — Columns(N) ▾ 面板: All/None/Key 快捷键 + 逐列复选框;
              >15 列且关键列命中 >=4 的宽表默认精简视图
  4. 列固定 — 每个表头悬停出现 📌, 点击把该列 sticky 钉在左侧,
              多列按 DOM 顺序累进偏移, 斑马纹背景保留, 隐藏列自动解钉
  附: Export CSV 导出 (含隐藏列), 分页/排序后状态自动保持 (MutationObserver)。

适用选择器: 页面内所有 <table>, 需满足 — 有 <thead> (首行全 TH 的无头表会
自动规范化) 且数据行 >=3; 布局用小表自动跳过; 容器加 class="tw-skip" 可豁免。

用法 (各报告生成器):
    try:
        from mmpv_common.table_widgets import TW_BLOCK as _TW
    except Exception:
        _TW = ""          # 独立运行/包不可达时静默降级
    然后在 HTML 模板 </body> 前放 {_TW} (f-string) 即可。
    只要 JS 不要 CSS 时用 TABLE_WIDGET_JS; 只要 CSS 用 TABLE_WIDGET_CSS。
"""

TABLE_WIDGET_CSS = """.twbar{display:flex;gap:8px;align-items:center;flex-wrap:wrap;margin:6px 0}
.twbar input{width:200px;padding:4px 8px;border:1px solid #ccc;border-radius:4px;font-size:11px}
.twbar button{padding:3px 10px;font-size:11px;border:1px solid #bbb;border-radius:4px;background:#fff;cursor:pointer;font-weight:600}
.twbar button:hover{background:#f2f0ea}
.twpanel{display:none;margin:4px 0;padding:6px 8px;background:#faf9f5;border:1px solid #eee;border-radius:4px;font-size:11px;max-height:240px;overflow:auto}
.twpanel label{display:inline-block;margin:2px 12px 2px 0;white-space:nowrap}
.twpanel .twq{margin:0 6px 4px 0;padding:2px 8px;font-size:11px;border:1px solid #ccc;border-radius:3px;background:#fff;cursor:pointer}
table.tw-done thead th{cursor:pointer;user-select:none}
.twpin{opacity:.25;cursor:pointer;margin-left:5px;font-size:9px;vertical-align:middle;filter:grayscale(1)}
thead th:hover .twpin{opacity:1}
.twpin.on{opacity:1;filter:none}
.twpinned{box-shadow:2px 0 0 0 rgba(0,0,0,.06),5px 0 8px -4px rgba(0,0,0,.18)}
.twarr{font-size:9px;margin-left:3px;color:#666}"""

TABLE_WIDGET_JS = """(function(){
var KEY=['contig_id','dataset','sample','run','length','Family','Genus','Species','Nucleic_acid','Genome_Type','checkv_completeness','n_samples','prevalence_%','nt_pident','aa_pident','final_verdict','final_call','prediction','novel','category_9a','host_assign'];
function normHead(t){
  if(t.tHead) return true;
  var r0=t.rows[0]; if(!r0||r0.cells.length<2) return false;
  for(var i=0;i<r0.cells.length;i++){ if(r0.cells[i].tagName!=='TH') return false; }
  var th=document.createElement('thead'); th.appendChild(r0); t.insertBefore(th,t.firstChild);
  return true;
}
function vis(t,items){t.querySelectorAll('tr').forEach(function(tr){items.forEach(function(it){var c=tr.cells[it.ci];if(c)c.style.display=it.cb.checked?'':'none';});});}
function cellsOf(t,ci){var out=[];for(var r=0;r<t.rows.length;r++){var c=t.rows[r].cells[ci];if(c)out.push(c);}return out;}
function install(){
 document.querySelectorAll('table').forEach(function(t){
  if(t._tw||t.closest('.tw-skip'))return;
  if(!normHead(t))return;
  var tb=t.tBodies[0];if(!tb||tb.rows.length<2)return;
  t._tw=1;t.classList.add('tw-done');
  var ths=[].slice.call(t.tHead.rows[0].cells);
  var wrap=t.closest('div')||t.parentNode,host=wrap.parentNode||wrap;
  /* ── 工具条: 筛选 + Columns + Export (先于 📌 图标构建, 列名才干净) ── */
  var bar=document.createElement('div');bar.className='twbar';
  var inp=document.createElement('input');inp.type='text';inp.placeholder='Filter...';bar.appendChild(inp);
  var panel=null,items=null;
  if(ths.length>=6){
    var colBtn=document.createElement('button');colBtn.textContent='Columns ('+ths.length+') \\u25BE';bar.appendChild(colBtn);
    panel=document.createElement('div');panel.className='twpanel';
    function qb(x){var b=document.createElement('button');b.className='twq';b.textContent=x;panel.appendChild(b);return b;}
    var allB=qb('All'),noneB=qb('None'),keyB=qb('Key');
    items=[];ths.forEach(function(th,i){
      var lab=document.createElement('label');var cb=document.createElement('input');cb.type='checkbox';cb.checked=true;
      lab.appendChild(cb);lab.appendChild(document.createTextNode(' '+th.textContent.trim().slice(0,28)));panel.appendChild(lab);
      items.push({cb:cb,ci:i,name:th.textContent.trim()});
    });
    var keyIdx=items.filter(function(it){return KEY.indexOf(it.name)>=0;}).map(function(it){return it.ci;});
    function setKey(){items.forEach(function(it){it.cb.checked=keyIdx.indexOf(it.ci)>=0;});vis(t,items);}
    allB.onclick=function(){items.forEach(function(it){it.cb.checked=true;});vis(t,items);layout();};
    noneB.onclick=function(){items.forEach(function(it){it.cb.checked=false;});vis(t,items);};
    keyB.onclick=setKey;
    panel.querySelectorAll('input').forEach(function(cb){cb.addEventListener('change',function(){
      items.forEach(function(it){ if(!it.cb.checked){var th=ths[it.ci];if(pinned.indexOf(th)>=0)toggle(th);} });
      vis(t,items);});
    });
    colBtn.onclick=function(){panel.style.display=panel.style.display==='none'?'block':'none';};
    if(keyIdx.length>=4&&ths.length>15)setKey();
  }
  var ex=document.createElement('button');ex.textContent='Export CSV';bar.appendChild(ex);
  host.insertBefore(bar,wrap);if(panel)host.insertBefore(panel,wrap);
  /* ── 筛选 ── */
  function applyFilter(){var q=inp.value.toLowerCase();t.querySelectorAll('tbody tr').forEach(function(tr){tr.style.display=tr.textContent.toLowerCase().indexOf(q)>=0?'':'none';});}
  inp.addEventListener('input',applyFilter);
  /* ── 排序 ── */
  ths.forEach(function(th){
    th.title='Click to sort';
    th.addEventListener('click',function(e){
      if(e.target&&e.target.classList&&e.target.classList.contains('twpin'))return;
      var tb=t.tBodies[0];if(!tb)return;
      var ci=th.cellIndex,asc=th._asc=!th._asc;
      var rows=Array.prototype.slice.call(tb.rows);
      rows.sort(function(a,b){
        var av=a.cells[ci]?a.cells[ci].textContent.trim():'';var bv=b.cells[ci]?b.cells[ci].textContent.trim():'';
        var na=parseFloat(av.replace(/,/g,'')),nb=parseFloat(bv.replace(/,/g,''));
        if(!isNaN(na)&&!isNaN(nb))return (na-nb)*(asc?1:-1);
        return av.localeCompare(bv)*(asc?1:-1);
      });
      rows.forEach(function(r){tb.appendChild(r);});
      [].forEach.call(t.tHead.rows[0].cells,function(h){var a=h.querySelector('.twarr');if(a)a.remove();});
      var arr=document.createElement('span');arr.className='twarr';arr.textContent=asc?'\\u25B2':'\\u25BC';th.appendChild(arr);
      if(inp.value)applyFilter();
    });
  });
  /* ── 导出 ── */
  ex.onclick=function(){
    var out=[];t.querySelectorAll('tr').forEach(function(tr){var row=[];tr.querySelectorAll('th,td').forEach(function(c){row.push('"'+c.innerText.replace(/"/g,'""')+'"');});out.push(row.join(','));});
    var blob=new Blob([out.join('\\n')],{type:'text/csv'});var a=document.createElement('a');a.href=URL.createObjectURL(blob);a.download=(t.id||'table')+'.csv';a.click();};
  /* ── 列固定: 每个表头 📌, sticky 累进偏移 ── */
  var pinned=[];
  function layout(){var x=0;pinned.forEach(function(th){var ci=th.cellIndex;cellsOf(t,ci).forEach(function(c){
    if(!c._twbg){c._twbg=1;c.dataset.twbg=c.style.backgroundColor||'';}
    c.style.position='sticky';c.style.left=x+'px';
    var bg=getComputedStyle(c).backgroundColor;
    c.style.background=(bg&&bg!=='rgba(0, 0, 0, 0)'&&bg!=='transparent')?bg:'#FCFAF5';
    c.classList.add('twpinned');});x+=th.offsetWidth;});}
  function unpinCell(ci){cellsOf(t,ci).forEach(function(c){c.style.position='';c.style.left='';c.style.background=c.dataset.twbg||'';c.classList.remove('twpinned');});}
  function toggle(th){
    var ix=pinned.indexOf(th),pin=th.querySelector('.twpin');
    if(ix>=0){pinned.splice(ix,1);unpinCell(th.cellIndex);if(pin)pin.classList.remove('on');}
    else{pinned.push(th);if(pin)pin.classList.add('on');}
    pinned.sort(function(a,b){return a.cellIndex-b.cellIndex;});layout();
  }
  /* 样式表里已有 sticky 首列的表 (如 KEEP 表) 视为已钉第 0 列 */
  if(ths.length&&getComputedStyle(ths[0]).position==='sticky'){pinned.push(ths[0]);var p0=document.createElement('span');p0.className='twpin on';p0.textContent='\\uD83D\\uDCCC';ths[0].appendChild(p0);p0.addEventListener('click',function(e){e.stopPropagation();toggle(ths[0]);});}
  if(ths.length>=4){
    ths.forEach(function(th){
      if(th.querySelector('.twpin'))return;
      var pin=document.createElement('span');pin.className='twpin';pin.textContent='\\uD83D\\uDCCC';pin.title='Pin/unpin column';
      th.appendChild(pin);
      pin.addEventListener('click',function(e){e.stopPropagation();toggle(th);});
    });
  }
  /* ── 分页/排序后保持状态 ── */
  new MutationObserver(function(){if(items)vis(t,items);layout();}).observe(tb,{childList:true});
  window.addEventListener('resize',layout);
 });
}
if(document.readyState==='loading')document.addEventListener('DOMContentLoaded',install);else install();
})();"""

TW_BLOCK = ('<style>\n' + TABLE_WIDGET_CSS + '\n</style>\n<script>\n' + TABLE_WIDGET_JS + '\n</script>')
