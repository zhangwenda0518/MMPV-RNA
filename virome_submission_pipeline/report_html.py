#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
report_html.py — 交互式全表编辑报告 v3.0
========================================

从 unified_metadata.csv + miuvig.tsv + assembly.tsv 生成完整可编辑 HTML。
所有表格均可双击编辑、批量填充、下载导出。

v3.0: 页面标记迁移到 mmpv_common/templates/submission_report.html.j2
(Jinja2 渲染, 表格走共享 data_table 宏, 元数据值经 autoescape); 交互 JS 仍在本文件
(JS_EDITABLE)。顺带修复了历史无定义 class (工具栏/统计卡/分区/Toast 此前无样式)。

用法:
  python report_html.py \
      --csv submission_metadata/unified_metadata.csv \
      --miuvig submission_metadata/miuvig.tsv \
      --assembly submission_metadata/assembly.tsv \
      --run-title my_project \
      -o submission_metadata/report.html
"""

import argparse, os
from pathlib import Path
from datetime import datetime
import pandas as pd

try:
    from mmpv_common.report_template import render
except ImportError:
    import sys; sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent.parent))
    from mmpv_common.report_template import render  # type: ignore

try:
    from unified_metadata import REQUIRED_COLS  # 同目录: 必填列单一事实源
except ImportError:
    import sys
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    try:
        from unified_metadata import REQUIRED_COLS  # type: ignore
    except ImportError:  # pragma: no cover - 缺依赖时退化为无必填标记, 不崩
        REQUIRED_COLS = []

JS_EDITABLE = """
const PLACEHOLDER = ['XXXX','YYYY','DD-Mmm','Country:Region','PRJNAXXXX','SAMNXXXX','XX.XX N','Author','not_provided'];

function isMissing(v){
  if(!v||v.trim()==='')return true;
  for(const p of PLACEHOLDER)if(v.includes(p))return true;
  return false;
}

function makeEditableTable(tableId){
  const tbl=document.getElementById(tableId);
  if(!tbl)return;
  const ths=tbl.querySelectorAll('thead th');
  tbl.querySelectorAll('tbody td').forEach(td=>{
    const full=td.getAttribute('data-full');
    const disp=td.textContent.trim();
    const val=full||disp;   /* edit/export use the original value (data-full); display stays truncated */
    const row=td.parentElement.rowIndex-1;
    const col=td.cellIndex;
    td.innerHTML='';
    const span=document.createElement('span');span.className='cell-content';span.textContent=disp;
    const inp=document.createElement('input');inp.value=val;
    inp.setAttribute('data-row',row);inp.setAttribute('data-col',col);
    td.appendChild(span);td.appendChild(inp);
    span.ondblclick=()=>{td.classList.add('editing');inp.focus();};
    inp.onblur=()=>{td.classList.remove('editing');span.textContent=inp.value;
      td.classList.toggle('missing',isMissing(inp.value));
      td.classList.toggle('filled',!isMissing(inp.value)&&inp.value.trim()!=='');
      updateStats();
    };
    inp.onkeydown=e=>{if(e.key==='Enter')inp.blur();if(e.key==='Escape'){inp.value=span.textContent;inp.blur();}};
    td.classList.toggle('missing',isMissing(val));
    td.classList.toggle('filled',!isMissing(val)&&val.trim()!=='');
  });
}

function collectTable(tableId){
  const tbl=document.getElementById(tableId);
  if(!tbl)return{headers:[],rows:[]};
  const headers=[];tbl.querySelectorAll('thead th').forEach(th=>headers.push(th.textContent.trim()));
  const rows=[];
  tbl.querySelectorAll('tbody tr').forEach(tr=>{
    const cols=[];tr.querySelectorAll('td').forEach(td=>{
      const inp=td.querySelector('input');
      const full=td.getAttribute('data-full');  /* export prefers the original value, else fall back to displayed text */
      cols.push(inp?inp.value:(full||td.querySelector('.cell-content')?.textContent||td.textContent||''));
    });rows.push(cols);
  });
  return{headers,rows};
}

function downloadTSV(tableId,filename){
  const{headers,rows}=collectTable(tableId);
  const lines=[headers.join('\\t')];
  rows.forEach(r=>lines.push(r.join('\\t')));
  const blob=new Blob([lines.join('\\n')],{type:'text/tab-separated-values'});
  const a=document.createElement('a');a.href=URL.createObjectURL(blob);
  a.download=filename;a.click();
  showToast(filename+' downloaded!');
}

function downloadCSV(tableId,filename){
  const{headers,rows}=collectTable(tableId);
  const lines=[headers.join(',')];
  rows.forEach(r=>lines.push(r.map(c=>'"'+c.replace(/"/g,'""')+'"').join(',')));
  const blob=new Blob([lines.join('\\n')],{type:'text/csv'});
  const a=document.createElement('a');a.href=URL.createObjectURL(blob);
  a.download=filename;a.click();
  showToast(filename+' saved!');
}

function batchFill(tableId,colName,val){
  if(!val)val=prompt('Enter value for '+colName+':');
  if(!val)return;
  const tbl=document.getElementById(tableId);
  const headers=[];tbl.querySelectorAll('thead th').forEach(th=>headers.push(th.textContent.trim()));
  const ci=headers.indexOf(colName);
  if(ci<0){showToast('column not found: '+colName,true);return;}
  tbl.querySelectorAll('tbody tr').forEach(tr=>{
    const td=tr.querySelectorAll('td')[ci];
    if(!td)return;
    const inp=td.querySelector('input');
    const span=td.querySelector('.cell-content');
    if(inp)inp.value=val;
    if(span)span.textContent=val;
    td.classList.remove('missing');td.classList.add('filled');
  });
  updateStats();showToast('Filled '+colName+' = '+val);
}

function updateStats(){
  let total=0,missing=0;
  document.querySelectorAll('#srcTable tbody td').forEach(td=>{
    const inp=td.querySelector('input');
    const v=inp?inp.value:(td.querySelector('.cell-content')?.textContent||'');
    if(td.classList.contains('required-col')){total++;if(isMissing(v))missing++;}
  });
  const pct=total>0?Math.round((total-missing)/total*100):100;
  document.getElementById('completeness').textContent=pct+'%';
  document.getElementById('completeBar').style.width=pct+'%';
  document.getElementById('completeBar').style.background=pct>=90?'#4caf50':pct>=70?'#f9a825':'#e65100';
  document.getElementById('missingCount').textContent=missing;
}

function showToast(msg,isErr){
  const t=document.createElement('div');t.className='toast'+(isErr?' error':'');
  t.textContent=msg;document.body.appendChild(t);
  setTimeout(()=>t.classList.add('show'),50);
  setTimeout(()=>{t.classList.remove('show');setTimeout(()=>t.remove(),300);},2200);
}

function filterTable(q){
  q=q.toLowerCase();
  document.querySelectorAll('#srcTable tbody tr').forEach(r=>{
    r.style.display=!q||r.textContent.toLowerCase().includes(q)?'':'none';
  });
}

function toggleSection(id){
  document.getElementById(id).classList.toggle('collapsed');
}

window.addEventListener('DOMContentLoaded',()=>{
  makeEditableTable('srcTable');
  makeEditableTable('cmtTable');
  makeEditableTable('bsTable');
  makeEditableTable('miuvigTable');
  makeEditableTable('asmTable');
  makeEditableTable('topoTable');
  updateStats();
});
"""

def _cell(row, col):
    return str(row.get(col, '')) if pd.notna(row.get(col, '')) else ''


def safe_row(row, cols, max_len=80):
    """Generate a row, truncating long values (显示用; 原值见 full_row)"""
    cells = []
    for c in cols:
        v = _cell(row, c)
        if len(v) > max_len:
            v = v[:max_len-3] + '...'
        cells.append(v)
    return cells


def full_row(row, cols):
    """未截断的行 (td data-full → CSV/TSV 导出保真, 不丢长字段)"""
    return [_cell(row, c) for c in cols]


def load_table(tsv_path):
    """读 miuvig/assembly 参数表 → (列名, 行)。NaN 归一为 '' (此前会渲染成字面 nan)。"""
    if not tsv_path or not os.path.exists(tsv_path):
        return [], []
    df = pd.read_csv(tsv_path, sep='\t' if str(tsv_path).endswith('.tsv') else None,
                     engine='python')
    rows = [['' if pd.isna(v) else v for v in list(df.iloc[i])]
            for i in range(len(df))]
    return list(df.columns), rows


def _kv_rows(rows):
    """miuvig/assembly 行 → [参数, 值] 二列 (缺值补空)。"""
    return [[str(r[0]), '' if len(r) < 2 or pd.isna(r[1]) else str(r[1])]
            for r in rows]


def generate_html(csv_path, run_title, output_path, miuvig_path=None, asm_path=None,
                   suvtk_tax_path=None, suvtk_feat_path=None, hypo_path=None,
                   topo_path=None, gb_dir=None, log=None):
    """Generate the full interactive report.

    渲染: mmpv_common/templates/submission_report.html.j2 (表格走共享 data_table
    宏 + autoescape, 元数据值里的 <>& 不再破坏标记)。
    """
    df = pd.read_csv(csv_path)
    n_seqs = len(df)
    n_viruses = df['organism'].nunique() if 'organism' in df.columns else 0

    # ── suvtk summary stats ──
    suvtk_parts = []
    if suvtk_tax_path and os.path.isfile(suvtk_tax_path):
        # 读侧 encoding 必标: suvtk tsv 为 UTF-8, Windows cp936 下读中文即崩
        with open(suvtk_tax_path, encoding='utf-8', errors='replace') as f:
            n_tax = sum(1 for _ in f) - 1
        suvtk_parts.append(f'suvtk taxonomy: {n_tax} sequences classified')
    if suvtk_feat_path and os.path.isfile(suvtk_feat_path):
        # .tbl CDS 行格式: <start>\t<end>\tCDS (负链为反转区间; suvtk features.py 原样写出)
        with open(suvtk_feat_path, encoding='utf-8', errors='replace') as f:
            n_cds = sum(1 for l in f if l.rstrip().endswith('CDS'))
        suvtk_parts.append(f'suvtk features: {n_cds} CDS annotated')
    if hypo_path and os.path.isfile(hypo_path):
        n_hypo = sum(1 for l in open(hypo_path, encoding='utf-8', errors='replace')
                     if 'hypothetical' in l.lower())
        n_anno = sum(1 for l in open(hypo_path, encoding='utf-8', errors='replace')
                     if 'inference\talignment' in l)
        suvtk_parts.append(f'hypothetical proteins: {n_anno} annotated, {n_hypo} remaining')
    if topo_path and os.path.isfile(topo_path):
        rows = []
        with open(topo_path, encoding='utf-8', errors='replace') as f:
            hdr = f.readline().strip().split('\t')
            for line in f:
                cols = line.strip().split('\t')
                if len(cols) >= len(hdr):
                    rows.append(dict(zip(hdr, cols)))
        n_linear = sum(1 for r in rows if 'linear' in r.get('final_topology', '').lower())
        n_circular = len(rows) - n_linear
        suvtk_parts.append(f'topology: {n_linear} linear / {n_circular} circular')
    if gb_dir and os.path.isdir(gb_dir):
        n_gb = len([f for f in os.listdir(gb_dir) if f.endswith('.gb')])
        if n_gb:
            suvtk_parts.append(f'GenBank files: {n_gb} .gb ready')

    # --- source.src table ---
    src_cols = [c for c in df.columns if not c.startswith('cmt-') and not c.startswith('bs-')]
    src_rows = [safe_row(df.iloc[i], src_cols) for i in range(n_seqs)]
    src_full = [full_row(df.iloc[i], src_cols) for i in range(n_seqs)]

    # --- cmt table ---
    cmt_cols = [c for c in df.columns if c.startswith('cmt-')]
    if cmt_cols:
        cmt_headers = cmt_cols
        cmt_rows = [safe_row(df.iloc[i], cmt_cols) for i in range(n_seqs)]
        cmt_full = [full_row(df.iloc[i], cmt_cols) for i in range(n_seqs)]
    else:
        cmt_headers, cmt_rows = ['Parameter', 'Value'], [['Assembly_Method', 'MEGAHIT;1.2.9']]
        cmt_full = None

    # --- bs table ---
    bs_cols = [c for c in df.columns if c.startswith('bs-')]
    if bs_cols:
        bs_headers = bs_cols
        bs_rows = [safe_row(df.iloc[i], bs_cols) for i in range(n_seqs)]
        bs_full = [full_row(df.iloc[i], bs_cols) for i in range(n_seqs)]
    else:
        bs_headers, bs_rows = ['Parameter', 'Value'], []
        bs_full = None

    # --- miuvig / assembly tables ---
    _, miuvig_raw = load_table(miuvig_path)
    _, asm_raw = load_table(asm_path)
    miuvig_rows = _kv_rows(miuvig_raw)
    asm_rows = _kv_rows(asm_raw)

    # --- topology table ---
    topo_rows = []
    topo_summary = 'No data'
    if topo_path and os.path.isfile(topo_path):
        topo_raw = []
        with open(topo_path, encoding='utf-8', errors='replace') as f:
            topo_hdr = f.readline().strip().split('\t')
            for line in f:
                cols = line.strip().split('\t')
                if len(cols) >= len(topo_hdr):
                    topo_raw.append(dict(zip(topo_hdr, cols)))
        if topo_raw:
            n_linear = sum(1 for r in topo_raw if 'linear' in r.get('final_topology', '').lower())
            n_circular = sum(1 for r in topo_raw if 'circular' in r.get('final_topology', '').lower())
            topo_summary = f'{len(topo_raw)} seqs: {n_linear} linear / {n_circular} circular'
            for row in topo_raw:
                topo = row.get('final_topology', '').lower()
                cls = 'filled' if 'linear' in topo else ('missing' if 'circular' not in topo else '')
                topo_rows.append({
                    'contig': row.get('contig', ''),
                    'taxonomy': row.get('taxonomy', ''),
                    'rule_topology': row.get('rule_topology', ''),
                    'final_topology': row.get('final_topology', ''),
                    'cls': cls,
                    'seq_length': row.get('seq_length', ''),
                    'dtr': '✓' if row.get('dtr_detected', 'False').lower() == 'true' else '',
                    'evidence': row.get('evidence', ''),
                })

    # --- checklist (由实际输入推导, 不再恒 done) ---
    _has_tax = bool(suvtk_tax_path) and os.path.isfile(suvtk_tax_path)
    _has_feat = bool(suvtk_feat_path) and os.path.isfile(suvtk_feat_path)
    _has_hypo = bool(hypo_path) and os.path.isfile(hypo_path)
    checklist = [
        ('done' if _has_tax else 'todo',
         'suvtk taxonomy — ICTV classification complete' if _has_tax
         else 'suvtk taxonomy — no input provided (--suvtk-tax)'),
        ('done' if _has_feat else 'todo',
         'suvtk features — ORF prediction + BFVD annotation' if _has_feat
         else 'suvtk features — no input provided (--suvtk-feat)'),
        ('done' if _has_hypo else 'todo',
         'analyze_hypothetical — diamond + HMM annotation' if _has_hypo
         else 'analyze_hypothetical — no input provided (--hypo)'),
        ('todo', 'Edit source.src above → download → save CSV'),
        ('todo', 'Register BioProject at https://submit.ncbi.nlm.nih.gov/'),
        ('todo', 'Register BioSamples using biosample_template.tsv'),
        ('todo', 'Run: suvtk comments -t miuvig_taxonomy -f miuvig_features -m miuvig.tsv -a assembly.tsv'),
        ('todo', 'Run: suvtk table2asn → email submission.sqn to gb-sub@ncbi.nlm.nih.gov'),
    ]

    html = render(
        "submission_report.html.j2",
        run_title=run_title,
        timestamp=datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
        pipeline='MMPV-RNA v2.3 + suvtk v0.1.1',
        n_seqs=n_seqs, n_viruses=n_viruses,
        suvtk_parts=suvtk_parts,
        src_headers=src_cols, src_rows=src_rows, src_full=src_full,
        required_cols=REQUIRED_COLS,
        cmt_headers=cmt_headers, cmt_rows=cmt_rows, cmt_full=cmt_full,
        bs_headers=bs_headers, bs_rows=bs_rows, bs_full=bs_full,
        miuvig_rows=miuvig_rows, asm_rows=asm_rows,
        topo_rows=topo_rows, topo_summary=topo_summary,
        checklist=checklist,
        js=JS_EDITABLE,
    )

    with open(output_path, 'w', encoding='utf-8') as f:
        f.write(html)

    if log:
        log.info("  → %s (%d seqs)", output_path, n_seqs)
    print(f"\n  Report: {output_path} | {n_seqs} seqs, {n_viruses} viruses")
    print(f"  Double-click any cell to edit | All tables downloadable")

    return output_path


def main():
    parser = argparse.ArgumentParser(description="report_html.py v3.0 — Full interactive submission report")
    parser.add_argument('--csv', required=True, help='unified_metadata.csv')
    parser.add_argument('--miuvig', help='miuvig.tsv')
    parser.add_argument('--assembly', help='assembly.tsv')
    parser.add_argument('--suvtk-tax', help='suvtk taxonomy.tsv')
    parser.add_argument('--suvtk-feat', help='suvtk featuretable.tbl')
    parser.add_argument('--hypo', help='analyze_hypothetical featuretable_updated.tbl')
    parser.add_argument('--topology', help='topology.tsv')
    parser.add_argument('--gb-dir', help='virus-annotations directory')
    parser.add_argument('--run-title', default='viral_submission')
    parser.add_argument('-o', '--output', default='report.html')
    args = parser.parse_args()
    generate_html(args.csv, args.run_title, args.output,
                  args.miuvig, args.assembly,
                  args.suvtk_tax, args.suvtk_feat, args.hypo,
                  args.topology, args.gb_dir)

    # AI 解读 section (可选, 纯增量: 离线/未配 key 时自动跳过, 不改动报告)
    try:
        from mmpv_common.ai import append_ai_section_to_html
        n_seqs = 0
        with open(args.csv, encoding="utf-8") as _f:
            for _line in _f:
                if _line.strip():
                    n_seqs += 1
        append_ai_section_to_html(Path(args.output), {
            "pipeline": "submission", "run_title": args.run_title,
            "n_seqs": max(n_seqs - 1, 0),
        })
    except Exception:
        pass

    # 过程面板 (确定性, 无需 AI): 项目 90_Handoff/state 在时嵌入全过程总览
    try:
        from mmpv_common.process_report import append_process_panel_to_html
        append_process_panel_to_html(Path(args.output), Path(args.output))
    except Exception:
        pass


if __name__ == '__main__':
    main()
