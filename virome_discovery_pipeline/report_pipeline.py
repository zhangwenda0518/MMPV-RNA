#!/usr/bin/env python3
"""
report_pipeline.py — 独立报告生成器 (MMPV-RNA v3.1)

根据流水线输出目录生成:
  - 各阶段汇总 TSV (data/assembly/ident/filter/cobra/hostdep/checkv)
  - Sankey 分类图 (全部 + 植物病毒)
  - pipeline_report.html (期刊级交互式 HTML)

用法:
  python report_pipeline.py -o out/                    # 从流水线输出根目录生成
  python report_pipeline.py -o out/ --skip-sankey      # 跳过 Sankey 图
"""

import argparse, csv, json, os, re, sys, subprocess, time
from datetime import datetime
from pathlib import Path
from Bio import SeqIO

# 跨管线统一 I/O 布局 (mmpv_common/, 仓库根): 目录名随 MMPV_IO_LAYOUT 解析
# (编排器已 normalize 环境变量, 子进程导入本模块时快照即正确布局)
_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)
from mmpv_common.io_layout import layout_dirs as _layout_dirs, build_discovery_dirs as _build_discovery_dirs
from mmpv_common.table_widgets import TABLE_WIDGET_CSS as _TW_CSS, TABLE_WIDGET_JS as _TW_JS
_D = _layout_dirs(os.environ.get("MMPV_IO_LAYOUT"))

SCRIPT_DIR = Path(__file__).resolve().parent
VERSION = "v3.1"

# ═══════════════════════════════════════════════════════════════
# 统一 SCI 配色 —— 唯一真源 = virome_discovery_pipeline/utils/sci_palette.py
# ═══════════════════════════════════════════════════════════════
# 2026-10-12: 此前本文件 chart.js 的**数据色**是 ~60 处游离 hex
# (#8B2C1F / #537D96 / #4A6B4A / #ffa726 / #00838f / #90caf9 / #88D6B2 …)，
# 与静态 SCI 图 (Fig_GojiVirome_Discovery / Fig_KnownVirome_Composite) 不同色系。
# 现统一到 sci_palette：
#   分类(无序) → OKABE_ITO   ·  序数(有序) → SEQ_CMAP=viridis
#   语义(好坏/判定) → SEMANTIC
# ⚠ 加载必须按**显式文件路径**：服务器 site-packages 里有同名第三方 `utils`
#   常规包 (有 __init__.py)，优先级高于本仓库的命名空间包，`import utils.…`
#   与 sys.path 顺序无关地必然失败。
_PAL = None
for _rel in ("virome_discovery_pipeline/utils/sci_palette.py",
             "virome_phylo_pipeline/utils/sci_palette.py"):
    _cand = Path(_REPO_ROOT) / _rel
    if _cand.is_file():
        import importlib.util as _ilu
        _spec = _ilu.spec_from_file_location("_sci_palette", _cand)
        _mod = _ilu.module_from_spec(_spec)
        sys.modules.setdefault("_sci_palette", _mod)
        _spec.loader.exec_module(_mod)
        _PAL = _mod
        break


class _PaletteFallback:                       # 真源缺失也不让报告生成失败
    OKABE_ITO = ['#E69F00', '#56B4E9', '#009E73', '#F0E442', '#0072B2',
                 '#D55E00', '#CC79A7', '#000000', '#999999']
    OKABE_BRANCH = ['#0072B2', '#D55E00', '#009E73', '#E69F00', '#CC79A7']
    SEQ_CMAP = "viridis"
    DIV_CMAP = "RdBu_r"
    SEMANTIC = {}


if _PAL is None:
    _PAL = _PaletteFallback()

OI = list(_PAL.OKABE_ITO)                     # 9 色，顺序固定
_SEM = getattr(_PAL, "SEMANTIC", None) or {}


def _sem(group, key, default):
    """取 sci_palette.SEMANTIC[group][key]，缺则回退 default。"""
    return (_SEM.get(group) or {}).get(key, default)


# 语义色（全部取自 OKABE_ITO，无游离 hex）
_C_TOTAL = "#DDDDDD"                          # 总量 / 中性底
_C_MUTED = OI[8]                              # #999999 次要 / 缺测
_C_BAD = OI[5]                                # #D55E00 负向 / 去除
_C_GOOD = OI[2]                               # #009E73 正向 / 保留
_C_WARN = OI[0]                               # #E69F00 警戒 / 中等
# 三档阈值配色（越高越好 / 越低越好），供 N50·contigs·dup 等阈值条共用
_TIER_GOOD = [OI[8], OI[0], OI[2]]            # 低 → 中 → 高
_TIER_BAD = [OI[2], OI[0], OI[5]]             # 好 → 警戒 → 差
# 序数色带：viridis(lo + (hi-lo)·i/(n-1))，lo=0.12 hi=0.90，与 sci_palette.SEQ_CMAP 同源。
# 本模块只产出 HTML/JS、不依赖 matplotlib，故色值必须字面量固化。
_VIRIDIS_RAMPS = {
    3: ['#472a7a', '#20928c', '#bddf26'],
    4: ['#472a7a', '#2c738e', '#2cb17e', '#bddf26'],
    5: ['#472a7a', '#33638d', '#20928c', '#46c06f', '#bddf26'],
    6: ['#472a7a', '#38598c', '#27808e', '#21a585', '#5ac864', '#bddf26'],
    7: ['#472a7a', '#3b528b', '#2c738e', '#20928c', '#2cb17e', '#69cd5b', '#bddf26'],
    8: ['#472a7a', '#3d4d8a', '#306a8e', '#25858e', '#1fa088', '#3aba76', '#73d056', '#bddf26'],
}


def _ramp(n):
    """viridis 上取 n 档（3–8；其余按 0.12–0.90 现插值）。"""
    if n in _VIRIDIS_RAMPS:
        return list(_VIRIDIS_RAMPS[n])
    if n == 1:
        return ['#20928c']
    lo, hi = 0.12, 0.90
    return ['#472a7a' if i == 0 else '#bddf26' if i == n - 1 else '#20928c'
            for i in range(n)]


# ═══════════════════════════════════════════════════════════════
# 辅助函数
# ═══════════════════════════════════════════════════════════════



_EN_MAP = {
'未定':'Undetermined','新病毒候选':'Novel candidate','宿主':'Host','非病毒':'Non-viral',
'未通过任何分支':'No branch passed','远缘病毒候选':'Distantly-related virus candidate','远缘':'Distantly related',
'宿主基因污染':'Host gene contamination','仅结构域证据':'Domain-only evidence','短片段':'Short fragment',
'待核验':'Pending','点击查看原图':'Click to view original','无病毒域':'No viral domain','无属分类':'No genus',
'组处理失败':'Group processing failed','蛋白':'Protein','免拯救':'No-rescue','样本':'Samples',
'已知植物病毒':'Known plant virus','昆虫病毒':'Insect virus','污染':'Contamination','已知病毒':'Known virus',
'下载':'Download','矩阵':'Matrix','新种候选':'Novel species candidate','长度达标但未被拯救':'Supported length, not rescued',
'已知':'Known','新种':'Novel species','分析验证':'Analysis & Validation','综合判定':'Final verdict','冲突':'Conflict',
'序列':'Sequence','新病毒判定':'Novelty call','序列快照':'Sequence snapshot','点击复制':'Click to copy','可复制':'Copyable',
'分类学':'Taxonomy','分类与':'Taxonomy &','鉴定':'Identification','新病毒':'Novel virus','高置信':'High confidence',
'进化树':'Phylogeny','个科建树':'families built','新属':'Novel genus',
'已知却去除':'Known but removed','蛋白近似却去除':'Protein-approx removed','有检出却去除':'Detected but removed',
'病毒':'virus','检测值轴':'value axis','横向图':'horizontal chart','值轴是':'value axis is',
'纵向图值轴是':'vertical chart value axis','近缘':'Closely related','无任何证据':'No evidence','延伸':'extended','延伸(partial)':'extended(partial)','延伸(circular)':'extended(circular)',
'昆虫病毒/污染':'Insect virus/contamination','真菌病毒/污染':'Fungal virus/contamination',
'噬菌体污染':'Phage contamination','无病毒域+blast弱':'No viral domain + weak blast','已检出但证据弱':'Detected but weak evidence',
'已运行但无最终结果':'Ran but no final result','未运行':'Not run','已运行':'Ran',
'其他':'Other','新科':'New family','远缘(未定)':'Distantly related (undetermined)','近缘(未定)':'Closely related (undetermined)',
# 数据派生词 (gen_final_judgement.py category 列中残存的混合词) → 显示层英文
'噬菌体Contamination':'Phage contamination','真菌病毒/Contamination':'Fungal virus/contamination',
'短片段/低覆盖':'Short fragment / low coverage',
}
def _en(t):
    """把显示文本中的中文词映射为英文 (report 全英文化)。"""
    for zh in sorted(_EN_MAP, key=len, reverse=True):
        t = t.replace(zh, _EN_MAP[zh])
    return t

_esf_js = "(function(){function esf(){document.querySelectorAll('canvas').forEach(function(cv){if(cv.closest('.chart-box')&&!cv._ex){cv._ex=true;var b=document.createElement('button');b.textContent='PNG';b.style.cssText='position:absolute;top:6px;right:8px;padding:2px 8px;font-size:11px;border:1px solid #ccc;border-radius:4px;background:#fff;cursor:pointer;z-index:5;';var cb=cv.closest('.chart-box');cb.style.position='relative';cb.appendChild(b);b.onclick=function(){try{var u=cv.toDataURL('image/png');var a=document.createElement('a');a.href=u;a.download=(cv.id||'chart')+'.png';a.click();}catch(e){}};}});document.querySelectorAll('.stage').forEach(function(st){if(st.querySelector('.stage-header')&&!st.querySelector('.chart-box')&&!st.querySelector('.stage-table')&&!st.querySelector('.stage-charts')){var d=document.createElement('div');d.style.cssText='padding:12px;color:var(--muted);font-size:12px;';d.textContent='No data for this stage';st.appendChild(d);}});}if(document.readyState==='loading')document.addEventListener('DOMContentLoaded',esf);else esf();})();"
# 表格四件套 (筛选/排序/列可选/列固定 + CSV 导出) 统一来自 mmpv_common.table_widgets
_keep_js = _TW_JS

def _resolve_acv_dir(base):
    """返回 09b 目录: 优先新目录 09b_Analysis_Verify, 回退旧名 09b_ACVirus_Analysis"""
    nd = Path(base) / _D['d_verify']
    if nd.is_dir():
        return nd
    return Path(base) / "09b_ACVirus_Analysis"

def _resolve_analysis_dir(base):
    """返回 09 分析目录: 按 io_layout 的 _analysis_dirs 顺序解析
    (legacy: 09a_Virome_Analysis 优先 → 回退旧名 09_Virome_Analysis;
     standard: 11_ViromeAnalysis)。新旧输出树都能读
    （同 _cobra_dirs / _centroids_v1 的双名容错做法）。
    """
    for d in _build_discovery_dirs(Path(base), Path(base),
                                   os.environ.get("MMPV_IO_LAYOUT"))['_analysis_dirs']:
        if d.is_dir():
            return d
    return Path(base) / _D['d_analysis']

# rescue 目录下的基础设施子目录 (非宿主), 不参与宿主子目录解析
_NON_HOST_SUBDIRS = ("known", "no_rescue", "checkv")


def _resolve_host_subdirs(root):
    """rescue 目录下的宿主子目录候选名 (优先 "Plant", 其后为其余宿主子目录)

    多宿主运行 (--host-filter Plant,Fungi → rescue 目录名 "Plant_Fungi") 下,
    硬编码 "Plant" 会取不到产物而静默为空。现枚举 rescue 下含
    centroids/final_centroids.fasta 或 rescue_report.tsv 的宿主子目录
    (排除 known/no_rescue/checkv), 按目录名排序保证确定性。
    返回候选名列表; "Plant" 始终在首位 (存在性由调用方 is_file() 判定)。
    """
    names = ["Plant"]
    rescue = Path(root) / _D['d_rescue']
    if rescue.is_dir():
        for d in sorted(p for p in rescue.iterdir() if p.is_dir()):
            if d.name in names or d.name in _NON_HOST_SUBDIRS:
                continue
            if (d / "centroids" / "final_centroids.fasta").is_file() \
               or (d / "rescue_report.tsv").is_file():
                names.append(d.name)
    return names

def _count_fasta(path):
    if not path or not os.path.isfile(str(path)): return 0
    with open(str(path), encoding='utf-8') as f:
        return sum(1 for _ in f if _.startswith('>'))

def _count_lines(path):
    if not path or not os.path.isfile(str(path)): return 0
    with open(str(path), encoding='utf-8') as f:
        return sum(1 for _ in f)

def _count_dir(d):
    return sum(1 for _ in Path(d).iterdir() if _.is_dir()) if d and Path(d).is_dir() else 0

def _tool_id_file(d, tool):
    """样本目录 d 下某工具的结果 id 文件。常规为 {sample}_virus.{tool}.result.id;
    类病毒 BLASTN 产物名特殊: {sample}_viroids.result.id (virus_identification.py)。"""
    if tool == 'viroid_blastn':
        return d / f"{d.name}_viroids.result.id"
    return d / f"{d.name}_virus.{tool}.result.id"

def _read_tsv(path):
    rows = []
    p = Path(path)
    if not p.is_file(): return rows
    # newline='' 保留原始行尾 → 显式剥 CR (生产端既有产物可能是 CRLF/混合行尾)
    with open(p, encoding='utf-8', newline='') as f:
        hdr = f.readline().rstrip('\r\n').split('\t')
        for line in f:
            line = line.rstrip('\r\n')
            if not line.strip(): continue
            rows.append(dict(zip(hdr, line.split('\t'))))
    return rows

def _n50_n90(lens):
    if not lens: return (0, 0)
    lens_s = sorted(lens, reverse=True)
    total = sum(lens_s); cum = 0; n50 = n90 = 0
    half, n90t = total / 2, total * 0.9
    for la in lens_s:
        cum += la
        if n50 == 0 and cum >= half: n50 = la
        if n90 == 0 and cum >= n90t: n90 = la
        if n50 > 0 and n90 > 0: break
    return (n50, n90)

def _parse_fasta_lens(fa_path):
    lens = []; seq = ""
    with open(str(fa_path), encoding='utf-8') as f:
        for line in f:
            s = line.strip()
            if s.startswith('>'):
                if seq: lens.append(len(seq))
                seq = ""
            else: seq += s
    if seq: lens.append(len(seq))
    return lens


def _img_to_base64(path, max_kb=3000):
    """Convert image to base64 for embedding in HTML. Skip if > max_kb."""
    import base64
    if not path or not Path(path).exists():
        return None
    size = Path(path).stat().st_size
    if size == 0:                      # 0 字节 (写失败/占位文件) → 不内嵌
        return None
    if size / 1024 > max_kb:
        return None
    with open(str(path), 'rb') as f:
        return base64.b64encode(f.read()).decode('utf-8')


# ═══════════════════════════════════════════════════════════════
# 阶段数据收集 (各阶段独立函数)
# ═══════════════════════════════════════════════════════════════

def _collect_cleandata(root, report_dir, _add):
    """00a_CleanData 阶段"""
    clean = root / _D['d_clean']
    if not clean.is_dir():
        _add(_D['d_clean'], "○", details="Not run"); return
    ns = len(set(f.name.split("_")[0] for d in [clean / "3.clumpify", clean / "2.fasta"] if d.is_dir() for f in d.iterdir() if f.is_file())) or _count_dir(clean / "3.clumpify") or _count_dir(clean / "2.fasta") or _count_dir(clean)
    _add(_D['d_clean'], "✓", key_metric=f"{ns} samples", details=str(clean))
    fp = clean / "logs"
    if not fp.is_dir(): return
    jf_list = list(fp.glob("*_fastp_report.json"))
    if not jf_list: return
    n_total_before, n_total_after = 0, 0
    b_total_before, b_total_after = 0, 0  # bases
    with open(report_dir / "data_summary.tsv", "w", encoding='utf-8') as ds:
        ds.write("Sample\tRaw_Reads\tClean_Reads\tRetained(%)\tRaw_Bases\tClean_Bases\tRaw_Q20(%)\tClean_Q20(%)\tRaw_Q30(%)\tClean_Q30(%)\tLowQ_Reads\tTooShort_Reads\tDup_Rate(%)\n")
        for jf in sorted(jf_list):
            try:
                with open(jf, encoding='utf-8') as jfh:
                    js = json.load(jfh)
                sn = jf.name.replace("_fastp_report.json", "")
                bef = js.get("summary", {}).get("before_filtering", {})
                aft = js.get("summary", {}).get("after_filtering", {})
                fil = js.get("filtering_result", {})
                dup = js.get("duplication", {})
                nb = bef.get("total_reads", 0); na = aft.get("total_reads", 0)
                bb = bef.get("total_bases", 0); ba = aft.get("total_bases", 0)
                n_total_before += nb; n_total_after += na
                b_total_before += bb; b_total_after += ba
                ds.write(f"{sn}\t{nb}\t{na}\t{round(na/max(nb,1)*100,1)}\t"
                         f"{bb}\t{ba}\t"
                         f"{round(bef.get('q20_rate',0)*100,1)}\t{round(aft.get('q20_rate',0)*100,1)}\t"
                         f"{round(bef.get('q30_rate',0)*100,1)}\t{round(aft.get('q30_rate',0)*100,1)}\t"
                         f"{fil.get('low_quality_reads',0)}\t{fil.get('too_short_reads',0)}\t"
                         f"{round(dup.get('rate',0)*100,3)}\n")
            except Exception as e: print(f"  [WARN] 解析 fastp JSON 失败 ({jf.name}): {e}")
        ds.write(f"TOTAL\t{n_total_before}\t{n_total_after}\t{round(n_total_after/max(n_total_before,1)*100,1)}\t"
                 f"{b_total_before}\t{b_total_after}\n")
    _add("  └ data_summary", "✓", key_metric=f"reads: {n_total_before:,}→{n_total_after:,}  bases: {b_total_before:,}→{b_total_after:,}")


def _collect_hostdepletion(root, report_dir, _add):
    """00b_HostDepletion 阶段"""
    hostdep = root / _D['d_hostdep']
    if not hostdep.is_dir():
        _add(_D['d_hostdep'], "○", details="Not run"); return
    ns = _count_dir(hostdep)
    hd_rows = {}
    sq_tsv = hostdep / "logs" / "host_depletion_seqkit_summary.tsv"
    if not sq_tsv.is_file():
        sq_tsv = hostdep / "host_depletion_seqkit_summary.tsv"
    if sq_tsv.is_file():
        for r in _read_tsv(sq_tsv):
            sn = r.get("Sample",""); stage = r.get("Stage","")
            nseq = int(r.get("num_seqs",0))
            if sn not in hd_rows: hd_rows[sn] = {"Sample":sn,"Raw":0,"After_Kraken2":0,"After_Host":0}
            if "Raw" in stage or "1_" in stage: hd_rows[sn]["Raw"] = max(hd_rows[sn]["Raw"], nseq)
            elif "Kraken" in stage or "2_" in stage: hd_rows[sn]["After_Kraken2"] = max(hd_rows[sn]["After_Kraken2"], nseq)
            elif "Host" in stage or "3_" in stage: hd_rows[sn]["After_Host"] = max(hd_rows[sn]["After_Host"], nseq)
    rr_tsv = hostdep / "logs" / "ribodetector.report.txt"
    if not rr_tsv.is_file():
        rr_tsv = hostdep / "ribodetector.report.txt"
    if rr_tsv.is_file():
        for r in _read_tsv(rr_tsv):
            sn = r.get("Sample","")
            if sn not in hd_rows: hd_rows[sn] = {"Sample":sn,"Raw":0,"After_Kraken2":0,"After_Host":0}
            hd_rows[sn]["rRNA"] = int(r.get("rRNA",0))
            hd_rows[sn]["non_rRNA"] = int(r.get("non_rRNA",0))
            hd_rows[sn]["Total_rRNA"] = int(r.get("Total_sequences",0))
    if hd_rows:
        with open(report_dir / "hostdep_summary.tsv", "w", encoding='utf-8') as hf:
            cols = ["Sample","Raw","After_Kraken2","After_Host","Total_rRNA","non_rRNA","rRNA"]
            hf.write("\t".join(cols)+"\n")
            for sn in sorted(hd_rows):
                r = hd_rows[sn]
                hf.write("\t".join(str(r.get(c,0)) for c in cols)+"\n")
    _add(_D['d_hostdep'], "✓", key_metric=f"{ns} samples", details=str(hostdep))


def _collect_bbnorm(root, report_dir, _add):
    """00c_BBnorm — 预处理侧的可选阶段 (非发现管线阶段)

    由 data_preprocessing_pipeline 产出 (data_preprocessing.py --stage bbnorm 或
    public_metadata_pipeline/preprocess_unified.py --bbnorm)，默认不跑；
    仅当预处理与发现管线共用同一输出根目录、且该目录存在时才展示。
    """
    bb = root / _D['d_bbnorm']
    if not bb.is_dir():
        return
    ns = _count_dir(bb)
    _add("00c_BBNorm", "✓", key_metric=f"{ns} samples normalized", details=str(bb))


def _collect_assembly(root, report_dir, _add):
    """01_Assembly 阶段"""
    asm = root / _D['d_asm']
    if not asm.is_dir():
        _add(_D['d_asm'], "○", details="Not run"); return
    ns = _count_dir(asm)
    total_contigs, total_bp = 0, 0
    asm_data = []
    with open(report_dir / "assembly_summary.tsv", "w", encoding='utf-8') as af:
        af.write("Sample\tAssembler\tSize(Mb)\tContigs\tMax_Len\tN50\tN90\t>500bp\t>500bp(%)\t>1000bp\t>1000bp(%)\n")
        for d in sorted(asm.iterdir()):
            if not d.is_dir(): continue
            sample_contigs = 0
            for fi, f in enumerate(sorted(d.glob("*.contig.fasta"))):
                lens = _parse_fasta_lens(f)
                if not lens: continue
                n = len(lens); total = sum(lens); mx = max(lens)
                n50, n90 = _n50_n90(lens)
                c500 = sum(1 for l in lens if l > 500); r500 = round(c500/n*100,1) if n else 0
                c1000 = sum(1 for l in lens if l > 1000); r1000 = round(c1000/n*100,1) if n else 0
                at = f.stem.replace(f"{d.name}_", "").replace(".contig", "")
                af.write(f"{d.name}\t{at}\t{total/1e6:.1f}\t{n}\t{mx}\t{n50}\t{n90}\t{c500}\t{r500}\t{c1000}\t{r1000}\n")
                asm_data.append({'s': d.name, 'a': at, 'n': n, 'total': total, 'n50': n50})
                if fi == 0:  # main assembler only
                    total_contigs += n; total_bp += total
                    sample_contigs += n
            _add(f"  └ {d.name}", "✓", key_metric=f"{sample_contigs} contigs")
        if len(asm_data) > 1:
            t_n = sum(r['n'] for r in asm_data); t_bp = sum(r['total'] for r in asm_data)
            af.write(f"TOTAL\tall\t{t_bp/1e6:.1f}\t{t_n}\t-\t-\t-\t-\t-\t-\t-\n")
    _add(_D['d_asm'], "✓", key_metric=f"{ns} samples, {total_contigs:,} contigs, {total_bp/1e6:.1f} Mb", details=str(asm))
    as_script = SCRIPT_DIR / "utils" / "assembly_stats.py"
    if as_script.is_file():
        as_out = report_dir / "assembly_detail"
        as_out.mkdir(parents=True, exist_ok=True)
        try: subprocess.run([sys.executable, str(as_script), "-a", str(asm), "-o", str(as_out / "assembly_summary.tsv")], capture_output=True, timeout=60)
        except Exception as e: print(f"  [WARN] assembly_stats.py 失败: {e}")


def _collect_identification(root, report_dir, _add):
    """02_Identification 阶段"""
    ident = root / _D['d_ident']
    if not ident.is_dir():
        _add("02_Identification", "○", details="Not run"); return
    ns = _count_dir(ident)
    n_virus = 0
    for d in ident.iterdir():
        if not d.is_dir(): continue
        for f in d.glob("*virus.all.candidate.fasta"): n_virus += _count_fasta(f)
    _add("02_Identification", "✓", key_metric=f"{ns} samples, {n_virus:,} viral sequences", details=str(ident))
    tools_list = ['genomad','blast','metabuli','virsorter2','viralverify','virhunter','virbot','viralm','rdrpcatch','viroid_blastn']
    ident_data = []
    with open(report_dir / "ident_summary.tsv", "w", encoding='utf-8') as ids:
        ids.write("Sample\tAll_Candidate\t" + "\t".join(tools_list) + "\n")
        for d in sorted(ident.iterdir()):
            if not d.is_dir(): continue
            all_ids = _count_fasta(d / f"{d.name}_virus.all.candidate.fasta")
            tcounts = {}
            for tool in tools_list:
                idf = _tool_id_file(d, tool)
                tcounts[tool] = _count_lines(idf) if idf.is_file() else 0
            ids.write(f"{d.name}\t{all_ids}\t" + "\t".join(str(tcounts[t]) for t in tools_list) + "\n")
            ident_data.append({'Sample': d.name, 'All': all_ids, **tcounts})
        if len(ident_data) > 1:
            total_all = sum(r['All'] for r in ident_data)
            total_tools = {t: sum(r[t] for r in ident_data) for t in tools_list}
            ids.write(f"TOTAL\t{total_all}\t" + "\t".join(str(total_tools[t]) for t in tools_list) + "\n")
            top_tools = sorted(total_tools.items(), key=lambda x: -x[1])[:3]
            best = " | ".join(f"{t}={c}" for t,c in top_tools)
            _add("  └ multi-sample", "✓", key_metric=f"total={total_all}, top={best}")
    filter_dir = root / _D['d_filter']
    if filter_dir.is_dir():
        filter_data = []; modes_seen = set()
        with open(report_dir / "filter_summary.tsv", "w", encoding='utf-8') as fs:
            fs.write("Sample\tMode\tAll_Candidate\tPassed\tRetained(%)\tNote\n")
            for d in sorted(filter_dir.iterdir()):
                if not d.is_dir(): continue
                all_n = _count_fasta(d / f"{d.name}_virus.all.candidate.fasta")
                for fm, ff in [('filter', d / f"{d.name}.virus.candidate_uniprot_filtered.fasta"),
                               ('strict', d / f"{d.name}.virus.candidate_filtered.fasta")]:
                    # 2026-10-19: 原实现用 ff.is_file() 判定，失效符号链接会被判 False
                    # 导致整表写 0 行。改为 is_symlink() 也纳入，并对失效链接单独计数。
                    nf = _count_fasta(ff) if (ff.is_file() or ff.is_symlink()) else 0
                    note = ''
                    if ff.is_symlink() and not ff.exists():
                        note = 'BROKEN_LINK'
                    # 始终写行（0 也是有效信息），避免报告整表为空
                    pct = "NA" if all_n == 0 else round(nf / all_n * 100, 1)
                    if all_n == 0:
                        note = (note + ";" if note else "") + "ZERO_CANDIDATES"
                    fs.write(f"{d.name}\t{fm}\t{all_n}\t{nf}\t{pct}\t{note}\n")
                    filter_data.append({'Sample': d.name, 'Mode': fm, 'All': all_n, 'Passed': nf})
                    modes_seen.add(fm)
            if filter_data:
                for m in sorted(modes_seen):
                    mr = [r for r in filter_data if r['Mode'] == m]
                    t_all = sum(r['All'] for r in mr); t_pass = sum(r['Passed'] for r in mr)
                    fs.write(f"TOTAL\t{m}\t{t_all}\t{t_pass}\t{round(t_pass/max(t_all,1)*100,1)}\t\n")
    # Per-tool raw/filter/strict 统计 (raw 取 02a 工具结果, filter/strict 回退 02b)
    with open(report_dir / "tool_filter_summary.tsv", "w", encoding='utf-8') as tfs:
        tfs.write("Tool\tRaw\tFilter\tStrict\n")
        tool_filter = {t: [0, 0, 0] for t in tools_list}  # [raw, filter, strict]
        for d in sorted(ident.iterdir()):
            if not d.is_dir(): continue
            fd2 = root / _D['d_filter'] / d.name
            for tool in tools_list:
                tool_filter[tool][0] += _count_lines(_tool_id_file(d, tool))
                for fi, sub in [(1,'candidate_uniprot_filtered'),(2,'candidate_filtered')]:
                    idf = fd2 / f"{d.name}.virus.{tool}.{sub}.id"
                    tool_filter[tool][fi] += _count_lines(idf) if idf.is_file() else 0
        for tool in sorted(tool_filter, key=lambda t: -tool_filter[t][0]):
            raw, filt, strict = tool_filter[tool]
            if raw > 0:
                tfs.write(f"{tool}\t{raw}\t{filt}\t{strict}\n")
    is_script = SCRIPT_DIR / "utils" / "ident_stats.py"
    if is_script.is_file():
        is_out = report_dir / "ident_detail"
        is_out.mkdir(parents=True, exist_ok=True)
        try: subprocess.run([sys.executable, str(is_script), "-i", str(ident), "-o", str(is_out)], capture_output=True, timeout=120)
        except Exception as e: print(f"  [WARN] ident_stats.py 失败: {e}")


def _collect_graph(root, report_dir, _add):
    """graph 阶段 (组装图标签传播; 目录 d_graph = 02b/03_Filter/graph, 默认不跑)

    有产物才登记: 逐样本 {sample}.graph_verdict.tsv → verdict 计数表
    (doc.md §6b; 契约见 graph_propagation.py VERDICT_COLUMNS)。
    """
    graph = root / _D['d_graph']
    if not graph.is_dir():
        return
    verdicts = sorted(graph.rglob("*.graph_verdict.tsv"))    # 产物在 d_graph/<样本>/ 下
    if not verdicts:
        return
    order = ['viral_seed', 'nonviral_seed', 'graph_viral', 'graph_nonviral', 'unknown', 'conflicting']
    counted = {}                      # {Sample: {verdict: n}}
    for vf in verdicts:
        sn = vf.name.replace(".graph_verdict.tsv", "")
        cnt = {}
        for r in _read_tsv(vf):
            v = (r.get('verdict') or '').strip() or 'unknown'
            cnt[v] = cnt.get(v, 0) + 1
        counted[sn] = cnt
    seen = {v for cnt in counted.values() for v in cnt}
    labels = [v for v in order if v in seen] + sorted(seen - set(order))
    with open(report_dir / "graph_verdict_summary.tsv", "w", encoding='utf-8') as gf:
        gf.write("Sample\t" + "\t".join(labels) + "\tTotal\n")
        tot = {v: 0 for v in labels}
        for sn in sorted(counted):
            cnt = counted[sn]
            for v in labels: tot[v] += cnt.get(v, 0)
            gf.write(f"{sn}\t" + "\t".join(str(cnt.get(v, 0)) for v in labels)
                     + f"\t{sum(cnt.values())}\n")
        gf.write("TOTAL\t" + "\t".join(str(tot[v]) for v in labels) + f"\t{sum(tot.values())}\n")
    _add(_D['d_graph'], "✓",
         key_metric=f"{len(counted)} samples, graph_viral={tot.get('graph_viral',0)}, graph_nonviral={tot.get('graph_nonviral',0)}",
         details=str(graph))


def _collect_cobra(root, report_dir, _add):
    """03_COBRA 阶段"""
    cobra = root / _D['d_cobra']
    if not cobra.is_dir():
        _add("03_COBRA", "○", details="Not run"); return
    ns = _count_dir(cobra)
    n_ext, n_queries, n_orphan, n_total_gain = 0, 0, 0, 0
    with open(report_dir / "cobra_summary.tsv", "w", encoding='utf-8') as cs:
        cs.write("Sample\tTotal_Queries\tExtended_Circular\tExtended_Partial\tExtended_Failed\tOrphan_End\tExtension_Rate(%)\tOrphan_Rate(%)\tExtended_Contigs\tTotal_Gain(bp)\n")
        for sd in sorted(cobra.iterdir()):
            if not sd.is_dir(): continue
            sn = sd.name
            sn_ext, sn_gain = 0, 0
            for cf in sd.rglob("*.cobra.fa"):
                if cf.stat().st_size > 0:
                    cobra_lens = _parse_fasta_lens(cf)
                    sn_ext += len(cobra_lens)
                    sn_gain += sum(cobra_lens)
            n_ext += sn_ext; n_total_gain += sn_gain
            logs = list(sd.rglob("log"))
            cobra_logs = [f for f in logs if 'COBRA' in str(f.parent.name)]
            if not cobra_logs: cobra_logs = logs[:1]
            tq = ec = ep = ef = oe = 0
            if cobra_logs:
                try:
                    with open(str(cobra_logs[0]), encoding='utf-8') as lf:
                        for line in lf:
                            s = line.strip()
                            try:
                                if s.startswith('# Total queries:'): tq = int(s.split(':')[1].strip())
                                elif 'Extended_circular' in s: ec = int(s.split(':')[1].strip().split()[0])
                                elif 'Extended_partial' in s: ep = int(s.split(':')[1].strip().split()[0])
                                elif 'Extended_failed' in s: ef = int(s.split(':')[1].strip())
                                elif 'Orphan end' in s: oe = int(s.split(':')[1].strip())
                            except (IndexError, ValueError): pass
                except Exception as e: print(f"  [WARN] COBRA log 解析失败 ({sd.name}): {e}")
            n_queries += tq; n_orphan += oe
            er = round((ec+ep)/max(tq,1)*100, 1); or_ = round(oe/max(tq,1)*100, 1)
            cs.write(f"{sn}\t{tq}\t{ec}\t{ep}\t{ef}\t{oe}\t{er}\t{or_}\t{sn_ext}\t{sn_gain}\n")
    _add("03_COBRA", "✓", key_metric=f"{ns} samples, {n_ext} extended, {n_queries} query, {n_orphan} orphan, {n_total_gain:,} bp gain")
    cs_script = SCRIPT_DIR / "utils" / "cobra_stats.py"
    if cs_script.is_file():
        cs_out = report_dir / "cobra_detail"
        cs_out.mkdir(parents=True, exist_ok=True)
        try: subprocess.run([sys.executable, str(cs_script), "-c", str(cobra), "-o", str(cs_out)], capture_output=True, timeout=120)
        except Exception as e: print(f"  [WARN] cobra_stats.py 失败: {e}")


def _collect_merge(root, report_dir, _add):
    """03b_MergeSamples + Flye 阶段"""
    merge_dir = root / _D['d_merge']
    if not merge_dir.is_dir():
        return  # merge是可选的，不显示
    merged_fa = merge_dir / "all_sample_virus.fasta"
    combined_fa = merge_dir / "all_sample_virus_combined.fasta"
    n_merged = _count_fasta(merged_fa) if merged_fa.is_file() else 0
    n_combined = _count_fasta(combined_fa) if combined_fa.is_file() else 0
    key = f"merge={n_merged:,}"
    flye_dir = merge_dir / "flye_coassembly"
    if flye_dir.is_dir():
        flye_fa = flye_dir / "flye_output" / "assembly.fasta"
        n_flye = _count_fasta(flye_fa) if flye_fa.is_file() else 0
        if n_flye > 0:
            key += f" + Flye={n_flye:,}"
    if n_combined > 0:
        key += f" → combined={n_combined:,}"
    _add("03b_Merge", "✓", key_metric=key, details=str(merge_dir))


def _collect_cluster(root, report_dir, _add):
    """04_CLUSTER 阶段"""
    cluster = root / _D['d_cluster']
    if not cluster.is_dir():
        _add(_D['d_cluster'], "○", details="Not run"); return
    centroids_fa = None
    for cand in ["4_centroids", "04_centroids", "centroids"]:
        cf = cluster / cand / "final_centroids.fasta"
        if cf.is_file(): centroids_fa = cf; break
    if centroids_fa is None: centroids_fa = cluster / "centroids" / "final_centroids.fasta"
    known_fa = cluster / "2_cdhit" / "known_centroids.fasta"
    n_centroids = _count_fasta(centroids_fa) if centroids_fa.is_file() else 0
    n_known = _count_fasta(known_fa) if known_fa.is_file() else 0
    _add(_D['d_cluster'], "✓", key_metric=f"{n_centroids:,} novel + {n_known:,} known centroids", details=str(cluster))
    ctsv = cluster / "3_vclust" / "vclust_clusters.tsv"
    if ctsv.is_file():
        n_contigs_in = _count_lines(ctsv) - 1
        cluster_map = {}
        with open(ctsv, encoding='utf-8') as cf:
            cf.readline()
            for line in cf:
                parts = line.strip().split()
                if len(parts) >= 2:
                    cid, clid = parts[0], parts[1]
                    cluster_map.setdefault(clid, []).append(cid)
        n_clusters = len(cluster_map)
        sizes = [len(v) for v in cluster_map.values()]
        singletons = sum(1 for sz in sizes if sz <= 1)
        max_sz = max(sizes) if sizes else 0
        _add("  └ vclust", "✓", key_metric=f"{n_clusters:,} clusters, {singletons} singletons, max cluster={max_sz}")
        size_bins = {'1 (singleton)':0,'2-3':0,'4-5':0,'6-10':0,'11-20':0,'21+':0}
        for sz in sizes:
            if sz <= 1: size_bins['1 (singleton)'] += 1
            elif sz <= 3: size_bins['2-3'] += 1
            elif sz <= 5: size_bins['4-5'] += 1
            elif sz <= 10: size_bins['6-10'] += 1
            elif sz <= 20: size_bins['11-20'] += 1
            else: size_bins['21+'] += 1
        with open(report_dir / 'cluster_size_distribution.tsv', 'w', encoding='utf-8') as csf:
            csf.write('Size\tCount\n')
            for k, v in size_bins.items():
                csf.write(f'{k}\t{v}\n')
        with open(report_dir / 'cluster_pipeline_reduction.tsv', 'w', encoding='utf-8') as prf:
            prf.write('Stage\tCount\n')
            prf.write(f'Input contigs\t{sum(sizes)}\n')
            prf.write(f'vOTU clusters\t{n_clusters}\n')
            prf.write(f'Final centroids\t{n_centroids}\n')
        skip_counts = {}
        # centroids 子目录名随布局变化 (写方在 4_centroids): 与上方最终 centroids
        # 同一三候选列表解析, 否则 cluster_skip_by_host.tsv 恒为空
        centroids_d = None
        for cand in ["4_centroids", "04_centroids", "centroids"]:
            cd = cluster / cand
            if cd.is_dir(): centroids_d = cd; break
        if centroids_d is not None:
            for sf in centroids_d.glob('skipped_*.fasta'):
                host = sf.stem.replace('skipped_', '').title()
                skip_counts[host] = _count_fasta(sf)
            unk = _count_fasta(centroids_d / 'unknown_votus.fasta')
            if unk > 0: skip_counts['Unknown'] = unk
            if skip_counts:
                with open(report_dir / 'cluster_skip_by_host.tsv', 'w', encoding='utf-8') as hsf:
                    hsf.write('Host\tCount\n')
                    for h, c in sorted(skip_counts.items(), key=lambda x: -x[1]):
                        hsf.write(f'{h}\t{c}\n')
    known_linked_fa = cluster / "2_cdhit" / "known_linked_centroids.fasta"
    n_linked = _count_fasta(known_linked_fa) if known_linked_fa.is_file() else 0
    if n_linked > 0: _add("  └ CD-HIT linked", "✓", key_metric=f"{n_linked} known clusters with linked contigs")
    if n_known - n_linked > 0: _add("  └ CD-HIT pure", "○", key_metric=f"{n_known - n_linked} pure reference clusters (not downstream)")


def _collect_taxonomy(root, report_dir, _add):
    """05_Taxonomy 阶段"""
    tax = root / _D['d_taxonomy']
    tax_files = _find_final_taxonomy(tax)
    # 多份产物 (逐样本 {sample}.integrated/) 按 contig_id 合并去重 (先出现者为准)
    _, tax_rows = _merge_tax_rows(tax_files) if tax_files else ([], [])
    if tax_rows:
        n = len(tax_rows)
        counts = {"Known": 0, "Novel_Species": 0, "Novel_Genus": 0, "Novel_Family": 0}
        rank_fill = {"Realm":0,"Kingdom":0,"Phylum":0,"Class":0,"Order":0,"Family":0,"Genus":0,"Species":0}
        for row in tax_rows:
            sp = row.get("Species", row.get("species", ""))
            ge = row.get("Genus", row.get("genus", ""))
            fa = row.get("Family", row.get("family", ""))
            if sp and sp not in ("NA", "-"): counts["Known"] += 1
            elif ge and ge not in ("NA", "-"): counts["Novel_Species"] += 1
            elif fa and fa not in ("NA", "-"): counts["Novel_Genus"] += 1
            else: counts["Novel_Family"] += 1
            for rk in ["Realm","Kingdom","Phylum","Class","Order","Family","Genus","Species"]:
                v = row.get(rk, row.get(rk.lower(), ""))
                if v and v not in ("NA", "-"): rank_fill[rk] += 1
        _add(_D['d_taxonomy'], "✓", key_metric=f"{n} seqs, ★{counts['Known']} known, ★★{counts['Novel_Species']} novel species", details=str(tax))
        _add("  └ Novel Rank", "✓", key_metric=f"Known={counts['Known']} NewSp={counts['Novel_Species']} NewGe={counts['Novel_Genus']} NewFa={counts['Novel_Family']}")
        rk_info = " ".join(f"{r}={rank_fill[r]}" for r in ["Realm","Phylum","Class","Order","Family","Genus","Species"] if rank_fill.get(r,0) > 0)
        if rk_info: _add("  └ Rank Fill", "✓", key_metric=rk_info[:120])
    elif tax.is_dir():
        _add(_D['d_taxonomy'], "○", key_metric="Ran but no final result", details=str(tax))
    else:
        _add(_D['d_taxonomy'], "○", details="Not run")

    # 最终植物病毒分类汇总 (HQ_plant_viruses.fasta × taxonomy)
    all_plant_fa = root / _D['d_rescue'] / "HQ_plant_viruses.fasta"
    if all_plant_fa.is_file() and tax_rows:
        plant_ids = set()
        for line in open(all_plant_fa, encoding='utf-8'):
            if line.startswith('>'): plant_ids.add(line[1:].split()[0])
        counts_p = {"Known": 0, "Novel_Species": 0, "Novel_Genus": 0, "Novel_Family": 0}
        for row in tax_rows:
            cid = row.get("contig_id", row.get("Contig_ID", ""))
            if cid not in plant_ids: continue
            sp = row.get("Species", row.get("species", ""))
            ge = row.get("Genus", row.get("genus", ""))
            fa = row.get("Family", row.get("family", ""))
            if sp and sp not in ("NA", "-"): counts_p["Known"] += 1
            elif ge and ge not in ("NA", "-"): counts_p["Novel_Species"] += 1
            elif fa and fa not in ("NA", "-"): counts_p["Novel_Genus"] += 1
            else: counts_p["Novel_Family"] += 1
        n_plant = sum(counts_p.values())
        if n_plant > 0:
            _add("  └ Plant Virus Taxonomy", "✓",
                 key_metric=f"{n_plant} seqs, Known={counts_p['Known']} NewSp={counts_p['Novel_Species']} NewGe={counts_p['Novel_Genus']} NewFa={counts_p['Novel_Family']}")

    # 收集 R 生成的 TSV + taxonomy composition + UpSet PNG
    integrated_dir = None
    for d in (root / _D['d_taxonomy']).glob("*.integrated"):
        if d.is_dir(): integrated_dir = d; break
    if integrated_dir:
        for src, dst in [("fill_stats.tsv", "taxonomy_fill_stats.tsv"),
                          ("agreement_stats.tsv", "taxonomy_agreement_stats.tsv"),
                          ("consistency_summary.tsv", "taxonomy_consistency_summary.tsv")]:
            sp = integrated_dir / src
            if sp.is_file():
                import shutil
                shutil.copy2(str(sp), str(report_dir / dst))
        # UpSet 交互式数据 (JSON) 优先, 静态 PNG 作为 fallback
        upset_data = integrated_dir / "upset_data.json"
        if upset_data.is_file():
            import shutil as _sh
            _sh.copy2(str(upset_data), str(report_dir / "upset_data.json"))
        else:
            upset_png = integrated_dir / "intersection_upset.png"
            b64 = _img_to_base64(upset_png, max_kb=2000)
            if b64:
                import json as _json
                with open(report_dir / "taxonomy_images.json", "w", encoding='utf-8') as imf:
                    _json.dump({"Tool Intersection UpSet": b64}, imf)
    # Taxonomy composition (from final_integrated_classification.tsv)
    if tax_rows:
        tax_comp = {}
        for row in tax_rows:
            for rk in ["Realm","Kingdom","Phylum","Class","Order","Family","Genus","Species"]:
                v = row.get(rk, row.get(rk.lower(), ""))
                if v and str(v).strip() not in ("","NA","-","N/A","None"):
                    tax_comp.setdefault(rk, {})
                    tax_comp[rk][v] = tax_comp[rk].get(v, 0) + 1
        with open(report_dir / "taxonomy_composition.tsv", "w", encoding='utf-8') as tcf:
            tcf.write("Rank\tTaxon\tCount\n")
            for rk in ["Realm","Kingdom","Phylum","Class","Order","Family","Genus","Species"]:
                if rk not in tax_comp: continue
                for taxon, cnt in sorted(tax_comp[rk].items(), key=lambda x: -x[1]):
                    tcf.write(f"{rk}\t{taxon}\t{cnt}\n")
        if "Family" in tax_comp:
            top_fam = sorted(tax_comp["Family"].items(), key=lambda x: -x[1])[:5]
            _add("  \u2514 Family", "\u2713", key_metric=" | ".join(f"{k}={v}" for k,v in top_fam)[:200])

    # 复制 final_integrated_classification.tsv 到 report_dir 供表格嵌入
    # (多份逐样本产物 → 先按 contig_id 合并去重再嵌入)
    if tax_files:
        import shutil as _shutil
        if len(tax_files) == 1:
            _shutil.copy2(str(tax_files[0]), str(report_dir / "final_integrated_classification.tsv"))
        else:
            _write_merged_tax(tax_files, report_dir / "final_integrated_classification.tsv")


def _collect_hostprediction(root, report_dir, _add):
    """06_HostPrediction 阶段"""
    host = root / _D['d_host_pred']
    host_summary = host / "ensemble_host_summary.tsv"
    if host_summary.is_file():
        n = _count_lines(host_summary) - 1
        try:
            import polars as pl
            hdf = pl.read_csv(str(host_summary), separator="\t", null_values=["NA","N/A",""])
            # Final_Host 分布
            if "Final_Host" in hdf.columns:
                hcounts = hdf.group_by("Final_Host").agg(pl.len().alias("count")).sort("count", descending=True)
                all_hosts = " | ".join(f"{r[0]}={r[1]}" for r in hcounts.iter_rows())
                _add(_D['d_host_pred'], "✓", key_metric=f"{n} seqs, {all_hosts}", details=str(host))
                # 写入 host_distribution.tsv
                with open(report_dir / "host_distribution.tsv", "w", encoding='utf-8') as hdfile:
                    hdfile.write("Host\tCount\n")
                    for r in hcounts.iter_rows():
                        hdfile.write(f"{r[0]}\t{r[1]}\n")
            # Decision_Method 分布
            if "Decision_Method" in hdf.columns:
                dm_counts = hdf.group_by("Decision_Method").agg(pl.len().alias("count")).sort("count", descending=True)
                with open(report_dir / "host_decision_method.tsv", "w", encoding='utf-8') as dmf:
                    dmf.write("Method\tCount\n")
                    for r in dm_counts.iter_rows():
                        dmf.write(f"{r[0]}\t{r[1]}\n")
        except Exception as e:
            print(f"  [WARN] HostPrediction polars 解析失败: {e}")
            _add(_D['d_host_pred'], "✓", key_metric=f"{n} seqs", details=str(host))
    elif host.is_dir():
        _add(_D['d_host_pred'], "○", key_metric="Ran but no final result", details=str(host))
    else:
        _add(_D['d_host_pred'], "○", details="Not run")

    # 复制 C9 ICTV 分类汇总表
    c9_summary = host / "C9_ICTV_result" / "classification_summary.tsv"
    if c9_summary.is_file():
        import shutil
        shutil.copy2(str(c9_summary), str(report_dir / "host_ictv_classification_summary.tsv"))
    c9_conf = host / "C9_ICTV_result" / "confidence_report.tsv"
    if c9_conf.is_file():
        import shutil
        shutil.copy2(str(c9_conf), str(report_dir / "host_ictv_confidence.tsv"))


def _collect_checkv(root, report_dir, _add):
    """07_CheckV 阶段"""
    cv_dir = root / _D['d_checkv']
    QUALITY_ORDER = ["Complete","High-quality","Medium-quality","Low-quality","Not-determined"]
    if not cv_dir.is_dir():
        _add("07_CheckV", "○", details="Not run"); return
    cv_tsvs = list(cv_dir.rglob("completeness.tsv"))
    if not cv_tsvs:
        _add("07_CheckV", "✓", key_metric="Ran", details=str(cv_dir)); return
    host_qd = {}; host_conf = {}
    for ct in cv_tsvs:
        host_name = ct.parent.name
        if host_name not in host_qd:
            host_qd[host_name] = dict.fromkeys(QUALITY_ORDER, 0)
            host_conf[host_name] = {}
        try:
            import polars as pl
            cv = pl.read_csv(str(ct), separator="\t", null_values=["NA","N/A",""])
            if "aai_confidence" in cv.columns:
                for row in cv.iter_rows(named=True):
                    conf = (row.get("aai_confidence") or "")
                    conf = str(conf).strip()   # None → "" (不能 str(None) 变 "None")
                    if not conf or conf.upper() in ("NA","N/A","NONE","NAN"): conf = "Not-determined"
                    host_conf[host_name][conf] = host_conf[host_name].get(conf,0) + 1
            comp_col = next((c for c in ["aai_completeness","completeness"] if c in cv.columns), None)
            if comp_col:
                for row in cv.iter_rows(named=True):
                    val = row.get(comp_col)
                    try:
                        v = float(val) if val and val != "NA" else None
                        if v is None: key = "Not-determined"
                        elif v >= 90: key = "Complete"
                        elif v >= 50: key = "High-quality"
                        elif v >= 10: key = "Medium-quality"
                        else: key = "Low-quality"
                    except (ValueError, TypeError): key = "Not-determined"
                    host_qd[host_name][key] += 1
        except Exception as e: print(f"  [WARN] CheckV 解析失败 ({ct.name}): {e}")
    global_qd = dict.fromkeys(QUALITY_ORDER, 0)
    for hqd in host_qd.values():
        for q in QUALITY_ORDER: global_qd[q] += hqd[q]
    n_total = sum(global_qd.values())
    n_hq = global_qd["Complete"] + global_qd["High-quality"]
    with open(report_dir / "checkv_summary.tsv", "w", encoding='utf-8') as cvf:
        cvf.write("Host\t" + "\t".join(QUALITY_ORDER) + "\tTotal\tHQ\n")
        for h in sorted(host_qd):
            hqd = host_qd[h]; t = sum(hqd.values()); hq = hqd["Complete"]+hqd["High-quality"]
            cvf.write(f"{h}\t"+"\t".join(str(hqd[q]) for q in QUALITY_ORDER)+f"\t{t}\t{hq}\n")
        cvf.write(f"TOTAL\t"+"\t".join(str(global_qd[q]) for q in QUALITY_ORDER)+f"\t{n_total}\t{n_hq}\n")
    if host_conf:
        all_confs = set()
        for hc in host_conf.values(): all_confs.update(hc.keys())
        conf_order = sorted(all_confs)
        with open(report_dir / "checkv_confidence.tsv", "w", encoding='utf-8') as cff:
            cff.write("Host\t"+"\t".join(conf_order)+"\tTotal\n")
            for h in sorted(host_conf):
                hc = host_conf[h]; t = sum(hc.values())
                cff.write(f"{h}\t"+"\t".join(str(hc.get(c,0)) for c in conf_order)+f"\t{t}\n")
    _add("07_CheckV", "✓", key_metric=f"{n_hq} HQ / {n_total} total", details=str(cv_dir))
    for h, hqd in sorted(host_qd.items(), key=lambda x: -sum(x[1].values()))[:8]:
        tot = sum(hqd.values())
        if tot > 0: _add(f"  └ {h}", "✓", key_metric=f"HQ={hqd['Complete']+hqd['High-quality']} total={tot}")
    qparts = [f"{q}={global_qd[q]}" for q in QUALITY_ORDER if global_qd[q] > 0]
    _add("  └ Distribution", "✓", key_metric=" ".join(qparts[:4]))


def _collect_rescue(root, report_dir, _add):
    """08_Rescue stage"""
    rescue = root / _D['d_rescue']
    if not rescue.is_dir():
        _add(_D['d_rescue'], "○", details="not run"); return
    hq_fa = rescue / "HQ_plant_viruses.fasta"
    n_hq = _count_fasta(hq_fa) if hq_fa.is_file() else 0
    known_fa = rescue / "known" / "centroids" / "final_centroids.fasta"
    n_known = _count_fasta(known_fa) if known_fa.is_file() else 0
    n_rescued = n_hq - n_known if n_hq > n_known else 0
    _add(_D['d_rescue'], "✓", key_metric=f"{n_hq} viruses (known={n_known} + rescued={n_rescued})", details=str(rescue))
    import shutil
    for rp in rescue.rglob("rescue_report.tsv"):
        shutil.copy2(str(rp), str(report_dir / "rescue_report.tsv")); break

def _collect_analysis(root, report_dir, _add):
    """10 Virome Analysis: result summary"""
    analysis = _resolve_analysis_dir(root)
    if not analysis.is_dir():
        _add("10_Virome_Analysis", "○", details="not run"); return
    parts = []
    hq_info = analysis / "HQ_analysis" / "HQ_plant_viruses_info.tsv"
    if hq_info.is_file():
        try:
            rows = _read_tsv(hq_info)
            parts.append(f"HQ viruses: {len(rows)} seqs")
        except: pass
    plant_info = analysis / "all_plant_analysis" / "All_plant.viruses_info.tsv"
    if plant_info.is_file():
        parts.append(f"All plant: {_count_lines(plant_info)-1} seqs")
    for fn, label in [("Viroid.all_info.tsv","Viroid all"),("Viroid.species_info.tsv","Viroid species")]:
        vf = analysis / "viroid_analysis" / fn
        if vf.is_file():
            parts.append(f"{label}: {_count_lines(vf)-1} seqs")
    # rescue 验证统计 (Part 4 检测频次 + virus_validation CDD 验证 + 证据整合, 09b)
    scored = (_resolve_acv_dir(analysis.parent) / "virus_validation" / "rescue_evidence_scored.tsv")
    if not scored.is_file():
        scored = analysis / "rescue_validation" / "rescue_evidence_scored.tsv"
    if scored.is_file():
        try:
            srows = _read_tsv(scored)
            keep = sum(1 for r in srows if r.get('verdict') == 'KEEP')
            review = sum(1 for r in srows if r.get('verdict') == 'REVIEW')
            drop = sum(1 for r in srows if r.get('verdict') == 'DROP')
            parts.append(f"Rescue: {keep}K/{review}R/{drop}D")
            # 检出频次统计 (n_samples > 0, 来自 salmon 检测)
            detected = 0
            for r in srows:
                try:
                    if int(float(r.get('n_samples') or 0)) > 0:
                        detected += 1
                except Exception:
                    pass
            parts.append(f"Detected: {detected}/{len(srows)}")
        except Exception:
            pass
    _add("10_Virome_Analysis", "✓", key_metric=" | ".join(parts) if parts else "run", details=str(analysis))


def _collect_analysis_files(root, report_dir):
    """从 09(a)_Virome_Analysis 子目录收集 TSV 文件 → report_dir (供图表/表格引用)"""
    import shutil
    analysis = _resolve_analysis_dir(root)
    if not analysis.is_dir():
        return
    copy_map = {
        analysis / "HQ_analysis" / "HQ_plant_viruses_info.tsv": "HQ_plant_viruses_info.tsv",
        analysis / "all_plant_analysis" / "All_plant.viruses_info.tsv": "All_plant.viruses_info.tsv",
        analysis / "all_plant_analysis" / "all_plant_viruses_genus_summary.tsv": "all_plant_viruses_genus_summary.tsv",
        analysis / "all_plant_analysis" / "plant_virus_cluster_summary.tsv": "plant_virus_cluster_summary.tsv",
        analysis / "viroid_analysis" / "Viroid.all_info.tsv": "Viroid.all_info.tsv",
        analysis / "viroid_analysis" / "Viroid.species_info.tsv": "Viroid.species_info.tsv",
        analysis / "viroid_analysis" / "Viroid.per_virus.tsv": "Viroid.per_virus.tsv",
        analysis / "viroid_analysis" / "Viroid.per_sample.tsv": "Viroid.per_sample.tsv",
        _resolve_acv_dir(analysis.parent) / "virus_validation" / "rescue_evidence_scored.tsv": ("rescue_evidence_scored.tsv", analysis / "rescue_validation" / "rescue_evidence_scored.tsv"),
        _resolve_acv_dir(analysis.parent) / "virus_validation" / "cdd_evidence_report.tsv": ("cdd_evidence_report.tsv", analysis / "rescue_validation" / "cdd_evidence_report.tsv"),
        _resolve_acv_dir(analysis.parent) / "virus_validation" / "hmm_ct3_hits.tsv": "hmm_ct3_hits.tsv",
        _resolve_acv_dir(analysis.parent) / "virus_validation" / "hmm_ct3_rescue_candidates.tsv": "hmm_ct3_rescue_candidates.tsv",
        _resolve_acv_dir(analysis.parent) / "keep_summary.tsv": "keep_summary.tsv",
        # s11 (09b) 章节表格: ACVirus 重分类/新病毒鉴定结果 (通用渲染循环只认
        # report_dir 相对文件名, 原目录不在其路径基准内 → 一并复制)
        _resolve_acv_dir(analysis.parent) / "acvirus_classify" / "final_result_with_confidence.tsv": "final_result_with_confidence.tsv",
        _resolve_acv_dir(analysis.parent) / "acvirus_identify" / "identification_result_score.csv": "identification_result_score.csv",
        analysis / "rescue_detection" / "frequency_table.tsv": "frequency_table.tsv",
        analysis / "rescue_detection" / "summary" / "all_viruses.best.summary.tsv": "rescue_detection_summary.tsv",
    }
    for src, entry in copy_map.items():
        dst_name = entry
        fallback = None
        if isinstance(entry, tuple):
            dst_name, fallback = entry
        srcf = src if src.is_file() else (fallback if fallback and fallback.is_file() else None)
        if srcf:
            shutil.copy2(str(srcf), str(report_dir / dst_name))


def collect_data(output_dir, report_dir, blast_db=None):
    """收集所有阶段数据 → 生成 TSV + 返回 stage_stats"""
    root = Path(output_dir).resolve()
    stage_stats = []

    def _add(stage, status, key_metric="", details=""):
        stage_stats.append({"Stage": stage, "Status": status, "Key_Metric": key_metric, "Details": details})

    _collect_cleandata(root, report_dir, _add)
    _collect_hostdepletion(root, report_dir, _add)
    _collect_bbnorm(root, report_dir, _add)
    _collect_assembly(root, report_dir, _add)
    _collect_identification(root, report_dir, _add)
    _collect_graph(root, report_dir, _add)
    _collect_cobra(root, report_dir, _add)
    _collect_merge(root, report_dir, _add)
    _collect_cluster(root, report_dir, _add)
    _collect_taxonomy(root, report_dir, _add)
    _collect_hostprediction(root, report_dir, _add)
    _collect_checkv(root, report_dir, _add)
    _collect_rescue(root, report_dir, _add)
    _collect_analysis(root, report_dir, _add)
    _collect_analysis_files(root, report_dir)
    _collect_acvirus(root, report_dir, _add)

    # 最终植物病毒汇总表 + 旭日图
    _generate_plant_virus_summary(root, report_dir, _add, blast_db)

    return stage_stats


def _collect_acvirus(root, report_dir, _add):
    """09b ACVirus Analysis - 收集 KEEP 重分类/新病毒鉴定/进化树"""
    acv = _resolve_acv_dir(root)
    if not acv.is_dir():
        # Stage 标签不能与 legacy 目录名逐字相同 (审计器 FAIL 精确字面量判定),
        # 保持 "09b" 前缀 → _stage_key_of 的 startswith('09b') → s11 映射不变
        _add("09b Analysis Verify", "○", "Not run")
        return
    try:
        import csv as _c
        import glob as _g
        # classify
        cls = acv / "acvirus_classify" / "final_result_with_confidence.tsv"
        keep_fa = acv / "class_KEEP.fasta"
        n_keep = sum(1 for _ in open(keep_fa, encoding='utf-8') if _.startswith('>')) if keep_fa.exists() else 0
        if cls.is_file():
            rows = list(_c.DictReader(open(cls, encoding='utf-8'), delimiter='\t'))
            fams = {}
            for r in rows:
                f = r.get('Family','') or ''
                if f and f != 'NA': fams[f] = fams.get(f,0)+1
            top = ', '.join(f"{k}×{v}" for k,v in sorted(fams.items(),key=lambda x:-x[1])[:5])
            _add("09b Analysis Verify", "✓", key_metric=f"{n_keep} KEEP, classified to {len(fams)} families ({top})",
                 details=str(acv))
        # identify
        idt = acv / "acvirus_identify" / "identification_result_score.csv"
        if idt.is_file():
            ir = list(_c.DictReader(open(idt, encoding='utf-8')))
            novel = sum(1 for r in ir if (r.get('Final_Verdict','') in ('High_Conf_Novel','Putative_Novel')))
            hn = sum(1 for r in ir if r.get('Final_Verdict','')=='High_Conf_Novel')
            _add("09b_Identification", "✓", key_metric=f"{len(ir)} sequences, novel {novel}, high-confidence {hn}")
        # trees
        trees = [d for d in acv.joinpath('acvirus_trees').iterdir() if d.is_dir()] if (acv/'acvirus_trees').is_dir() else []
        if trees:
            _add("09b_Phylogeny", "✓", key_metric=f"{len(trees)} families with trees: {', '.join(d.name for d in trees[:6])}")
        # keep_summary 综合判定 (final_call 整合 9a×ACVirus)
        ksum = acv / "keep_summary.tsv"
        if ksum.is_file():
            ks = list(_c.DictReader(open(ksum, encoding='utf-8'), delimiter='\t'))
            from collections import Counter as _Cnt
            fc = _Cnt(r.get('final_call','') for r in ks)
            conf = sum(1 for r in ks if r.get('conflict')=='Y')
            novel_total = sum(1 for r in ks if str(r.get('final_call','')).startswith('Novel'))
            _add("09b_Novelty_Call", "✓",
                 key_metric=f"Final verdict {len(ks)} KEEP: Novel×{novel_total} (dual {fc.get('Novel(dual)',0)}/blast {fc.get('Novel(blast)',0)}/ACV {fc.get('Novel(ACV)',0)}), Known×{fc.get('Known',0)+fc.get('Known(ACV-conflict)',0)}, conflicts {conf}")
        # SDT 属级一致性矩阵统计
        sdt_dir = acv / "SDT_matrix"
        if sdt_dir.is_dir():
            # 产物为 SDT_{属}_heatmap.png/.pdf (sdt_genus_matrix.py); 旧名一并兼容
            sdt_done = [d for d in sdt_dir.iterdir() if d.is_dir() and any(
                (d / f).exists() for f in (f"SDT_{d.name}_heatmap.png", f"SDT_{d.name}_heatmap.pdf",
                                           f"SDT_{d.name}.png", f"SDT_{d.name}.pdf"))]
            if sdt_done:
                _add("09b_SDT_matrix", "✓", key_metric=f"SDT genus matrix for {len(sdt_done)} genera: {', '.join(d.name for d in sdt_done[:6])}",
                     details=str(sdt_dir))
    except Exception as e:
        _add("09b Analysis Verify", "✗", key_metric=str(e)[:60])



def _generate_plant_virus_summary(root, report_dir, _add, blast_db=None):
    """生成 plant_virus_summary.tsv + taxonomy_sunburst.html"""
    all_plant_fa = root / _D['d_rescue'] / "HQ_plant_viruses.fasta"
    if not all_plant_fa.is_file(): return

    # 1. 读取 contig IDs + lengths + 来源 (宿主子目录动态解析, 兼容多宿主)
    plant_data = {}  # {contig_id: {length, source, ...}}
    host_dirs = _resolve_host_subdirs(root)
    src_cands = [("no-rescue", root / _D['d_rescue'] / "known" / "centroids" / "final_centroids.fasta")]
    src_cands += [("rescued", root / _D['d_rescue'] / hd / "centroids" / "final_centroids.fasta")
                  for hd in host_dirs]
    for src_label, src_fa in src_cands:
        if not src_fa.is_file(): continue
        for rec in SeqIO.parse(str(src_fa), "fasta"):
            plant_data[rec.id] = {"contig_id": rec.id, "length": len(rec.seq), "source": src_label}

    if not plant_data: return

    # 2. CheckV 数据: 优先 post-rescue (最新评估), 回退 CheckV 阶段原始值
    cv_data = {}
    cv_cands = [root / _D['d_rescue'] / "checkv" / hd / "completeness.tsv" for hd in host_dirs]
    cv_cands.append(root / _D['d_rescue'] / "checkv" / "no_rescue" / "completeness.tsv")
    cv_cands += [root / _D['d_checkv'] / hd / "completeness.tsv" for hd in host_dirs]
    for cv_tsv in cv_cands:
        if not cv_tsv.is_file(): continue
        rows = _read_tsv(cv_tsv)
        for r in rows:
            cid = r.get("contig_id","")
            if cid in plant_data and cid not in cv_data:
                cv_data[cid] = {
                    "aai_completeness": r.get("aai_completeness","NA"),
                    "aai_confidence": r.get("aai_confidence","NA"),
                    "viral_length": r.get("viral_length","NA"),
                    "aai_expected_length": r.get("aai_expected_length","NA"),
                    "kmer_freq": r.get("kmer_freq","NA"),
                }

    # 2.5. Rescue 数据: 挽救后长度 + 属平均长度 (来自 rescue_report.tsv)
    rescue_data = {}
    rs_cands = [root / _D['d_rescue'] / hd / "rescue_report.tsv" for hd in host_dirs]
    rs_cands.append(root / _D['d_rescue'] / "known" / "rescue_report.tsv")
    for rescue_tsv in rs_cands:
        if not rescue_tsv.is_file(): continue
        for r in _read_tsv(rescue_tsv):
            cid = r.get("contig_id","")
            if cid and cid in plant_data and cid not in rescue_data:
                rescue_data[cid] = {
                    "final_len": r.get("final_len",""),
                    "genus_avg_len": r.get("genus_avg_len",""),
                    "pct_of_genus": r.get("pct_of_genus",""),
                    "branch": r.get("branch",""),
                    "method": r.get("method",""),
                }

    # 3. Taxonomy (多份逐样本产物按 contig_id 合并去重, 先出现者为准)
    tax_data = {}
    tax_files = _find_final_taxonomy(root / _D['d_taxonomy'])
    for row in (_merge_tax_rows(tax_files)[1] if tax_files else []):
        cid = row.get("contig_id","")
        if cid in plant_data:
            tax_data[cid] = {rk: row.get(rk, row.get(rk.lower(),"")) for rk in
                             ["Realm","Kingdom","Phylum","Class","Order","Family","Genus","Species"]}

    # 写入 plant_virus_summary.tsv (全部作为新病毒 vOTU)
    cols = ["contig_id","length","final_len","genus_avg_len","pct_of_genus","source","branch","method",
            "aai_completeness","aai_confidence",
            "viral_length","aai_expected_length","kmer_freq",
            "Realm","Kingdom","Phylum","Class","Order","Family","Genus","Species"]
    with open(report_dir / "plant_virus_summary.tsv", "w", encoding='utf-8') as pvf:
        pvf.write("\t".join(cols) + "\n")
        for cid in sorted(plant_data):
            d = plant_data[cid]; cv = cv_data.get(cid, {}); tx = tax_data.get(cid, {}); rc = rescue_data.get(cid, {})
            vals = [cid, d["length"], rc.get("final_len",""), rc.get("genus_avg_len",""), rc.get("pct_of_genus",""),
                    d["source"], rc.get("branch",""), rc.get("method",""),
                    cv.get("aai_completeness","NA"), cv.get("aai_confidence","NA"),
                    cv.get("viral_length","NA"), cv.get("aai_expected_length","NA"), cv.get("kmer_freq","NA")]
            vals += [tx.get(rk,"") for rk in ["Realm","Kingdom","Phylum","Class","Order","Family","Genus","Species"]]
            pvf.write("\t".join(str(v) for v in vals) + "\n")

    n = len(plant_data)
    n_no_rescue = sum(1 for v in plant_data.values() if v["source"]=="no-rescue")
    n_rescued = n - n_no_rescue
    _add("09_Plant_virus", "✓",
         key_metric=f"{n} seqs (no-rescue={n_no_rescue} + rescued={n_rescued})",
         details=f"All treated as novel-virus vOTUs | {report_dir / 'plant_virus_summary.tsv'}")

    # 5. 旭日图 (Plotly Sunburst)
    if not tax_data: return
    try:
        import importlib.util
        if not importlib.util.find_spec("plotly"):
            subprocess.run([sys.executable, "-m", "pip", "install", "plotly", "-q"], check=False)
        import plotly.express as px
        import pandas as pd
        tax_rows = []
        for cid, d in plant_data.items():
            tx = tax_data.get(cid, {})
            row = {
                "Realm": tx.get("Realm","") or "Unclassified",
                "Kingdom": tx.get("Kingdom","") or "",
                "Phylum": tx.get("Phylum","") or "",
                "Class": tx.get("Class","") or "",
                "Order": tx.get("Order","") or "",
                "Family": tx.get("Family","") or "",
                "Genus": tx.get("Genus","") or "",
            }
            tax_rows.append(row)
        if tax_rows:
            df = pd.DataFrame(tax_rows)
            # 填充空值
            for col in df.columns:
                df[col] = df[col].replace("","Unclassified")
            path = ["Realm","Kingdom","Phylum","Class","Order","Family","Genus"]
            fig = px.sunburst(df, path=[p for p in path if p in df.columns],
                             title="Plant Virus Taxonomy Hierarchy",
                             height=700, width=900)
            fig.update_traces(textinfo="label+percent entry")
            fig.write_html(str(report_dir / "taxonomy_sunburst.html"),
                          include_plotlyjs='inline', full_html=True)
            print(f"  Sunburst → {report_dir / 'taxonomy_sunburst.html'}")
    except Exception as e:
        print(f"  [WARN] 旭日图生成失败: {e}")


# ═══════════════════════════════════════════════════════════════
# Sankey 图生成
# ═══════════════════════════════════════════════════════════════

def _find_final_taxonomy(tax_dir):
    """在 05_Taxonomy 下动态查找全部 final_integrated_classification.tsv (子目录名不固定)

    taxonomy 产物在逐样本布局下每样本一份 ({sample}.integrated/), 旧实现只取
    首份匹配 → 多样本项目丢失其余样本注释。现返回**全部**命中路径, 顺序为
    标准路径 → *integrated*/ → 其余子目录 (与 build_ref_info.py 的合并顺序一致);
    单份时列表长度 1, 行为与旧实现相同。
    """
    tax_dir = Path(tax_dir)
    cands = [tax_dir / "integrated" / "final_integrated_classification.tsv"]
    cands += sorted(tax_dir.glob("*integrated*/final_integrated_classification.tsv"))
    cands += sorted(tax_dir.glob("*/final_integrated_classification.tsv"))
    out = []
    for p in cands:
        if p.is_file() and p not in out:
            out.append(p)
    return out


def _merge_tax_rows(tax_files):
    """读取多份 taxonomy TSV, 按 contig_id 合并去重 (先出现者为准)

    参照 build_ref_info.py 的按 id 合并逻辑; 表头取首份文件。返回 (header, rows),
    rows 为 dict 列表; 单文件时即该文件原样行序。
    """
    header, rows, seen = None, [], set()
    for fp in tax_files:
        with open(fp, encoding='utf-8') as tf:
            rd = csv.DictReader(tf, delimiter="\t")
            if header is None:
                header = list(rd.fieldnames or [])
            for row in rd:
                cid = (row.get("contig_id") or row.get("Contig_ID") or "").strip()
                if cid:
                    if cid in seen: continue
                    seen.add(cid)
                rows.append(row)
    return (header or []), rows


def _write_merged_tax(tax_files, out_path):
    """把多份 taxonomy TSV 按 contig_id 合并写出 (供只接受单文件的下游脚本);
    单份输入时直接返回原路径, 不做拷贝 (行为与旧实现一致)。"""
    if len(tax_files) == 1:
        return tax_files[0]
    header, rows = _merge_tax_rows(tax_files)
    cols = header or ["contig_id"]
    with open(out_path, "w", encoding='utf-8', newline='') as of:
        of.write("\t".join(cols) + "\n")
        for row in rows:
            of.write("\t".join(str(row.get(c, "")) for c in cols) + "\n")
    return out_path

def generate_sankey(output_dir, report_dir):
    """生成交互式 taxonomy Sankey HTML (全部 + 植物病毒)"""
    tax = Path(output_dir) / _D['d_taxonomy']
    tax_files = _find_final_taxonomy(tax)
    if not tax_files: return
    # 多份产物 (逐样本 {sample}.integrated/) → 先合并为单文件再喂 sankey 脚本
    final_tax = _write_merged_tax(tax_files, report_dir / "final_integrated_classification.merged.tsv")
    sankey_script = SCRIPT_DIR.parent / "stats" / "taxonomic_sankey.py"
    if not sankey_script.is_file():
        sankey_script = SCRIPT_DIR / "utils" / "taxonomic_sankey.py"
    if not sankey_script.is_file(): return
    try:
        import importlib.util
        if not importlib.util.find_spec("plotly"):
            print("  安装 plotly...")
            subprocess.run([sys.executable, "-m", "pip", "install", "plotly", "-q"], check=False)
        # 生成交互式 HTML (--format html, 不设 title 避免重叠)
        _sk_html = report_dir / "classification_sankey.html"
        _r1 = subprocess.run([sys.executable, str(sankey_script),
                        "-i", str(final_tax), "-o", str(_sk_html),
                        "--format", "html", "--min-flow", "1", "--min-genus-flow", "10",
                        "--palette", "okabe_ito", "--height", "600", "--width", "800", "--node-pad", "30",
                        "--label-truncate", "25", "--font-size", "9", "--title-font-size", "14",
                        "--title", ""],
                       capture_output=True, timeout=120)
        if _r1.returncode == 0 and _sk_html.is_file():   # 成功才报 "→", 失败打 error (不再假成功)
            print(f"  Sankey HTML → {_sk_html}")
        else:
            print(f"  [ERROR] Sankey 生成失败 (rc={_r1.returncode}): "
                  f"{(_r1.stderr or b'').decode('utf-8','ignore')[-200:]}")
        # Plant-only
        host_summary = Path(output_dir) / _D['d_host_pred'] / "ensemble_host_summary.tsv"
        if host_summary.is_file():
            try:
                import polars as pl
                hdf = pl.read_csv(str(host_summary), separator="\t", null_values=["NA","N/A",""])
                plant_ids = set(hdf.filter(pl.col("Final_Host")=="Plant")["contig_id"].to_list())
                if plant_ids:
                    plant_tax = report_dir / "plant_final_taxonomy.tsv"
                    with open(final_tax, encoding='utf-8') as tf, open(plant_tax, "w", encoding='utf-8') as pf:
                        pf.write(tf.readline())
                        for line in tf:
                            cid = line.split('\t')[0].strip('"')
                            if cid in plant_ids: pf.write(line)
                    _sk_plant = report_dir / "classification_sankey_plant.html"
                    _r2 = subprocess.run([sys.executable, str(sankey_script),
                                    "-i", str(plant_tax), "-o", str(_sk_plant),
                                    "--format", "html", "--min-flow", "1", "--min-genus-flow", "5",
                                    "--palette", "okabe_ito", "--height", "600", "--width", "800", "--node-pad", "30",
                                    "--label-truncate", "25", "--font-size", "9", "--title-font-size", "14",
                                    "--title", ""],
                                   capture_output=True, timeout=120)
                    if _r2.returncode == 0 and _sk_plant.is_file():
                        print(f"  Plant Sankey HTML → {_sk_plant}")
                    else:
                        print(f"  [ERROR] Plant Sankey 生成失败 (rc={_r2.returncode}): "
                              f"{(_r2.stderr or b'').decode('utf-8','ignore')[-200:]}")
            except Exception as e: print(f"  [WARN] Plant Sankey 生成失败: {e}")
    except Exception as e: print(f"  [WARN] Sankey 生成失败: {e}")


# ═══════════════════════════════════════════════════════════════
# HTML 报告生成
# ═══════════════════════════════════════════════════════════════

def write_html_report(report_dir, stage_stats):
    """生成期刊级 HTML 流水线报告"""
    import json as _json

    # 内嵌 Chart.js (避免 CDN 不可用导致图表空白)
    chart_js_path = SCRIPT_DIR / "utils/chart.min.js"
    chart_js_inline = ""
    if chart_js_path.is_file():
        with open(chart_js_path, "r", encoding="utf-8") as cf:
            chart_js_inline = cf.read()
        # 同时复制到报告目录供离线使用
        import shutil
        shutil.copy2(str(chart_js_path), str(report_dir / "chart.min.js"))

    # 内嵌 UpSet.js (交互式 UpSet 图, 替代静态 PNG)
    upset_js_path = SCRIPT_DIR / "utils/upsetjs.umd.production.min.js"
    upset_js_inline = ""
    if upset_js_path.is_file():
        with open(upset_js_path, "r", encoding="utf-8") as uf:
            upset_js_inline = uf.read()

    upset_react_inline = ""
    upset_react_path = SCRIPT_DIR / "utils/upsetreact.umd.js"
    if upset_react_path.is_file():
        with open(upset_react_path, "r", encoding="utf-8") as rf:
            upset_react_inline = rf.read()

    # 读 upset_data.json (交互式 UpSet 数据, 由 R 脚本导出)
    upset_data_json = ""
    upset_container_html = ""
    upset_render_js = ""
    upset_data_path = report_dir / "upset_data.json"
    if upset_data_path.is_file():
        try:
            with open(upset_data_path, "r", encoding="utf-8") as uf:
                upset_data_json = uf.read().strip()
            if upset_data_json:
                upset_container_html = ('<div class="chart-box"><div class="chart-title">Tool Intersection (Interactive UpSet)</div>'
                                        '<div id="upset-container" style="min-height:420px"></div></div>')
                upset_render_js = '''(function(){
  var d=document.getElementById('upset-data');
  var el=document.getElementById('upset-container');
  if(!d||!el||!window.UpSetJS)return;
  try{
    var elems=JSON.parse(d.textContent);
    // ⚠ 防御: R jsonlite auto_unbox 会把单元素 sets 拆成字符串标量
    //   ("ACVirus" 而非 ["ACVirus"]) → UpSetJS.extractCombinations() 调
    //   .forEach() 抛 TypeError → 整图渲染中断。此处统一规范化为数组。
    elems.forEach(function(e){
      if(!e) return;
      if(e.sets==null){ e.sets=[]; }
      else if(typeof e.sets==='string'){ e.sets=[e.sets]; }
      else if(!Array.isArray(e.sets)){ e.sets=Array.from(e.sets); }
    });
    var comb=window.UpSetJS.extractCombinations(elems);
    var _pal=['#E69F00','#56B4E9','#009E73','#F0E442','#0072B2','#D55E00','#CC79A7','#000000','#999999'];
    comb.sets.forEach(function(s,i){s.color=_pal[i%_pal.length];});
    function fmt(n){if(n===null||n===undefined)return"";if(n>=1000)return (n/1000).toFixed(n>=10000?0:1)+"k";return String(n);}
    var _sel=null;
    function _r(){
      window.UpSetJS.renderUpSet(el,{sets:comb.sets,combinations:comb.combinations,width:1020,height:480,padding:16,fontSizes:{setLabel:"12px",axisTick:"9px",chartLabel:"13px",barLabel:"9px",legend:"9px",description:"12px",title:"16px",valueLabel:"11px",exportLabel:"9px"},widthRatios:[0.24,0.12,0.64],heightRatios:[0.55,0.45],barLabelOffset:3,setLabelAlignment:"center",fontFamily:"sans-serif",valueFormat:fmt,exportButtons:false,selection:_sel,onHover:function(sel){_sel=sel||null;_r();},onClick:function(sel){_sel=sel||null;_r();}});
    }
    _r();
  }catch(e){console.error("UpSet render failed:",e);}
})();'''
        except Exception:
            upset_data_json = ""

    main_stages = [s for s in stage_stats if not s["Stage"].startswith("  ")]

    S = {"✓": "pass", "○": "skip", "✗": "fail"}
    def _esc(v): return str(v).replace("&","&amp;").replace("<","&lt;").replace(">","&gt;").replace("'","&#39;")
    def _extract_kv(text, pattern):
        result = {}
        for m in re.finditer(pattern, text):
            result[m.group(1)] = int(m.group(2))
        return result

    n_pass = sum(1 for s in main_stages if s["Status"] == "✓")
    n_total = sum(1 for s in main_stages if s["Status"] != "○")
    pct = round(n_pass/max(n_total,1)*100)

    # ── 数据提取 ──
    chart_scripts = ""
    stage_has_chart = {}

    def _chart(canvas_id, chart_type, data_obj, options_obj=None):
        """生成 new Chart() JS 调用, 用 json.dumps 避免花括号转义问题"""
        cfg = {"type": chart_type, "data": data_obj}
        if options_obj:
            cfg["options"] = options_obj
        cfg_js = _json.dumps(cfg, ensure_ascii=False)
        return f"new Chart(document.getElementById('{canvas_id}'), {cfg_js});\n"

    def _tsv_data_uri(tsv_path):
        """将 TSV 文件转为 data:text/tab-separated-values;base64,... URI, 用于离线下载"""
        import base64
        if not Path(tsv_path).is_file(): return ""
        with open(tsv_path, "rb") as f:
            return "data:text/tab-separated-values;base64," + base64.b64encode(f.read()).decode()

    # S00a — 数据质量
    dq_rows = _read_tsv(report_dir / "data_summary.tsv")
    if dq_rows:
        dq_samples = [r.get("Sample","")[:14] for r in dq_rows if r.get("Sample","")!="TOTAL"]
        dq_raw = [int(r.get("Raw_Reads",0)) for r in dq_rows if r.get("Sample","")!="TOTAL"]
        dq_clean = [int(r.get("Clean_Reads",0)) for r in dq_rows if r.get("Sample","")!="TOTAL"]
        dq_q20b = [float(r.get("Raw_Q20(%)",0)) for r in dq_rows if r.get("Sample","")!="TOTAL"]
        dq_q20a = [float(r.get("Clean_Q20(%)",0)) for r in dq_rows if r.get("Sample","")!="TOTAL"]
        dq_dup = [float(r.get("Dup_Rate(%)",0)) for r in dq_rows if r.get("Sample","")!="TOTAL"]
        if dq_samples:
            stage_has_chart['s00a'] = True
            chart_scripts += _chart('chart_s00a', 'bar', {
                "labels": dq_samples,
                "datasets": [
                    {"label":"Raw Reads (M)","data":[round(v/1e6,1) for v in dq_raw],"backgroundColor":_C_MUTED,"yAxisID":"y"},
                    {"label":"Clean Reads (M)","data":[round(v/1e6,1) for v in dq_clean],"backgroundColor":OI[1],"yAxisID":"y"},
                    {"label":"Q20 Raw (%)","data":[round(v,1) for v in dq_q20b],"backgroundColor":_C_WARN,"yAxisID":"y1"},
                    {"label":"Q20 Clean (%)","data":[round(v,1) for v in dq_q20a],"backgroundColor":_C_GOOD,"yAxisID":"y1"},
                ]}, {"responsive":True,"plugins":{"title":{"display":True,"text":"Read Quality & Filtering"}},
                     "scales":{"y":{"beginAtZero":True,"position":"left","title":{"text":"Reads (M)"}},
                              "y1":{"beginAtZero":True,"position":"right","max":100,"grid":{"drawOnChartArea":False},"title":{"text":"Q20 (%)"}}}})
    if dq_rows and any(v > 0 for v in dq_dup):
        chart_scripts += _chart('chart_s00a_dup', 'bar', {
            "labels": dq_samples,
            "datasets": [{"label":"Duplication Rate (%)","data":[round(v,2) for v in dq_dup],
                          "backgroundColor":[_C_BAD if v>30 else _C_WARN if v>15 else _C_GOOD for v in dq_dup]}]},
            {"responsive":True,"plugins":{"title":{"display":True,"text":"Duplication Rate per Sample"}},
             "scales":{"y":{"beginAtZero":True,"title":{"text":"Dup Rate (%)"}}}})

    # S00b — 宿主去除
    hd_rows = _read_tsv(report_dir / "hostdep_summary.tsv")
    if hd_rows:
        hd_samples = [r.get("Sample","")[:12] for r in hd_rows]
        hd_raw = [int(r.get("Raw",0)) for r in hd_rows]
        hd_retained = []
        for r in hd_rows:
            nr = int(r.get("non_rRNA",0)); ah = int(r.get("After_Host",0)); ak = int(r.get("After_Kraken2",0))
            hd_retained.append(nr if nr > 0 else (ah if ah > 0 else ak))
        hd_removed = [max(0, hd_raw[i] - hd_retained[i]) for i in range(len(hd_samples))]
        if hd_samples and any(v > 0 for v in hd_raw):
            stage_has_chart['s00b'] = True
            chart_scripts += _chart('chart_s00b', 'bar', {
                "labels": hd_samples,
                "datasets": [{"label":"Retained (non-host)","data":hd_retained,"backgroundColor":_C_GOOD},
                              {"label":"Removed (host+rRNA)","data":hd_removed,"backgroundColor":_C_BAD}]},
                {"responsive":True,"indexAxis":"y","plugins":{"title":{"display":True,"text":"Host Depletion: Reads Retained vs Removed"}},
                 "scales":{"x":{"stacked":True,"beginAtZero":True,"title":{"text":"Reads"}},"y":{"stacked":True}}})

    # S01 — 组装
    asm_rows = _read_tsv(report_dir / "assembly_summary.tsv")
    if asm_rows:
        asm_samples = [(r["Sample"][:10] + " · " + (r.get("Assembler","") or "asm")[:12]) for r in asm_rows if r.get("Sample","")!="TOTAL"]
        asm_n50 = [int(r.get("N50",0)) for r in asm_rows if r.get("Sample","")!="TOTAL"]
        asm_contigs = [int(r.get("Contigs",0)) for r in asm_rows if r.get("Sample","")!="TOTAL"]
        if asm_samples:
            stage_has_chart['s01a'] = True
            chart_scripts += _chart('chart_s01a', 'bar', {
                "labels": asm_samples,
                "datasets": [{"label":"N50 (kb)","data":[round(v/1000,1) for v in asm_n50],
                              "backgroundColor":[_TIER_GOOD[2] if v>5000 else _TIER_GOOD[1] if v>1000 else _TIER_GOOD[0] for v in asm_n50]}]},
                {"responsive":True,"plugins":{"title":{"display":True,"text":"Assembly N50 (kb)"}},
                 "scales":{"y":{"beginAtZero":True,"title":{"text":"N50 (kb)"}}}})
            stage_has_chart['s01b'] = True
            chart_scripts += _chart('chart_s01b', 'bar', {
                "labels": asm_samples,
                "datasets": [{"label":"Contigs","data":asm_contigs,
                              "backgroundColor":[_TIER_GOOD[2] if v>5000 else _TIER_GOOD[1] if v>1000 else _TIER_GOOD[0] for v in asm_contigs]}]},
                {"responsive":True,"plugins":{"title":{"display":True,"text":"Assembly Contig Count"}},
                 "scales":{"y":{"beginAtZero":True,"title":{"text":"Contigs"}}}})

    # S02 — 每个工具的 raw/filter/strict 分组条形图
    tf_rows = _read_tsv(report_dir / "tool_filter_summary.tsv")
    if tf_rows:
        tools = [r["Tool"] for r in tf_rows]
        raw_vals = [int(r.get("Raw",0)) for r in tf_rows]
        filt_vals = [int(r.get("Filter",0)) for r in tf_rows]
        strict_vals = [int(r.get("Strict",0)) for r in tf_rows]
        if tools:
            stage_has_chart['s02'] = True
            chart_scripts += _chart('chart_s02', 'bar', {
                "labels": tools,
                "datasets": [
                    {"label":"Raw","data":raw_vals,"backgroundColor":_ramp(3)[0]},
                    {"label":"Filter","data":filt_vals,"backgroundColor":_ramp(3)[1]},
                    {"label":"Strict","data":strict_vals,"backgroundColor":_ramp(3)[2]},
                ]},
                {"responsive":True,"plugins":{"title":{"display":True,"text":"Per-Tool: Raw vs Filter vs Strict"}},
                 "scales":{"y":{"beginAtZero":True,"title":{"text":"Sequences"}}}})

    # S03 — COBRA
    cobra_rows = _read_tsv(report_dir / "cobra_summary.tsv")
    if cobra_rows:
        csamples = [r.get("Sample","")[:12] for r in cobra_rows]
        try: er_key = next(k for k in cobra_rows[0] if 'Extension_Rate' in k or 'extension_rate' in k)
        except: er_key = None
        try: or_key = next(k for k in cobra_rows[0] if 'Orphan_Rate' in k or 'orphan_rate' in k)
        except: or_key = None
        if er_key and csamples:
            stage_has_chart['s03'] = True
            cobra_ext = [float(r.get(er_key,0)) for r in cobra_rows]
            cobra_orph = [float(r.get(or_key,0)) for r in cobra_rows] if or_key else []
            ds_c = [{"label":"Extension Rate (%)","data":cobra_ext,"backgroundColor":_C_GOOD,"yAxisID":"y"}]
            if cobra_orph:
                ds_c.append({"label":"Orphan Rate (%)","data":cobra_orph,"backgroundColor":_C_BAD,"yAxisID":"y1"})
            scales_c = {"y":{"beginAtZero":True,"position":"left","title":{"text":"Rate (%)"}}}
            if cobra_orph:
                scales_c["y1"] = {"beginAtZero":True,"position":"right","grid":{"drawOnChartArea":False},"title":{"text":"Orphan (%)"}}
            chart_scripts += _chart('chart_s03', 'bar', {"labels":csamples,"datasets":ds_c},
                {"responsive":True,"plugins":{"title":{"display":True,"text":"COBRA Extension & Orphan Rate"}},"scales":scales_c})

    # S04 cluster charts
    csd_rows = _read_tsv(report_dir / "cluster_size_distribution.tsv")
    if csd_rows:
        csd_labels = [r.get("Size","") for r in csd_rows]
        csd_counts = [int(r.get("Count",0)) for r in csd_rows]
        csd_colors = _ramp(6)
        chart_scripts += _chart('chart_s04a', 'bar', {"labels": csd_labels, "datasets": [{"label":"Clusters","data":csd_counts,"backgroundColor":csd_colors}]}, {"responsive":True,"plugins":{"title":{"display":True,"text":"Cluster Size Distribution"},"legend":{"display":False}}, "scales":{"y":{"beginAtZero":True,"title":{"text":"Number of Clusters"}}}})
        stage_has_chart['s04a'] = True
    prd_rows = _read_tsv(report_dir / "cluster_pipeline_reduction.tsv")
    if prd_rows:
        prd_labels = [r.get("Stage","") for r in prd_rows]
        prd_counts = [int(r.get("Count",0)) for r in prd_rows]
        prd_colors = _ramp(3)
        chart_scripts += _chart('chart_s04c', 'bar', {"labels": prd_labels, "datasets": [{"label":"Count","data":prd_counts,"backgroundColor":prd_colors}]}, {"responsive":True,"indexAxis":"y", "plugins":{"title":{"display":True,"text":"Pipeline Reduction: Contigs to vOTUs"},"legend":{"display":False}}, "scales":{"x":{"beginAtZero":True,"title":{"text":"Count"}}}})
        stage_has_chart['s04c'] = True

    # S05 taxonomy composition charts
    tax_comp_rows = _read_tsv(report_dir / "taxonomy_composition.tsv")
    if tax_comp_rows:
        family_rows = [r for r in tax_comp_rows if r.get("Rank","") == "Family"]
        if family_rows:
            f_labels = [r.get("Taxon","")[:25] for r in family_rows[:20]]
            f_counts = [int(r.get("Count",0)) for r in family_rows[:20]]
            chart_scripts += _chart('chart_s05b', 'bar', {"labels": f_labels, "datasets": [{"label":"vOTUs","data":f_counts,"backgroundColor":OI[4]}]}, {"responsive":True,"indexAxis":"y","plugins":{"title":{"display":True,"text":"Family Distribution"},"legend":{"display":False}}, "scales":{"x":{"beginAtZero":True,"title":{"text":"vOTU Count"}}}})
            stage_has_chart['s05b'] = True
        genus_rows = [r for r in tax_comp_rows if r.get("Rank","") == "Genus"]
        if genus_rows:
            g_labels = [r.get("Taxon","")[:25] for r in genus_rows[:20]]
            g_counts = [int(r.get("Count",0)) for r in genus_rows[:20]]
            chart_scripts += _chart('chart_s05c', 'bar', {"labels": g_labels, "datasets": [{"label":"vOTUs","data":g_counts,"backgroundColor":OI[6]}]}, {"responsive":True,"indexAxis":"y","plugins":{"title":{"display":True,"text":"Genus Distribution"},"legend":{"display":False}}, "scales":{"x":{"beginAtZero":True,"title":{"text":"vOTU Count"}}}})
            stage_has_chart['s05c'] = True
    # R data charts
    fs_rows = _read_tsv(report_dir / "taxonomy_fill_stats.tsv")
    if fs_rows:
        fs_labels = [r.get("Rank","") for r in fs_rows]
        fs_rates = [round(float(r.get("Rate",0))*100, 1) for r in fs_rows]
        fs_filled = [int(r.get("Filled",0)) for r in fs_rows]
        chart_scripts += _chart('chart_s05d', 'bar', {"labels": fs_labels, "datasets": [{"label":"Fill Rate (%)","data":fs_rates,"backgroundColor":OI[4],"yAxisID":"y"}, {"label":"Filled vOTUs","data":fs_filled,"backgroundColor":OI[1],"yAxisID":"y1"}]}, {"responsive":True,"indexAxis":"y", "plugins":{"title":{"display":True,"text":"Taxonomy Fill Rate by Rank"}}, "scales":{"x":{"beginAtZero":True,"position":"bottom","title":{"text":"Fill Rate (%)"}},"y":{"position":"left"},"y1":{"display":False}}})
        stage_has_chart['s05d'] = True
    ag_rows = _read_tsv(report_dir / "taxonomy_agreement_stats.tsv")
    if ag_rows:
        levels_order = ["Realm","Kingdom","Phylum","Class","Order","Family","Genus","Species"]
        tools_set = sorted(set(r.get("Tool","") for r in ag_rows))
        ag_by_level = {}
        for r in ag_rows:
            lv = r.get("Taxonomic_Level","")
            if lv not in ag_by_level: ag_by_level[lv] = {}
            ag_by_level[lv][r.get("Tool","")] = round(float(r.get("Agreement_Rate",0))*100, 1)
        ag_levels = [lv for lv in levels_order if lv in ag_by_level]
        ag_colors = list(OI[:8])
        ag_datasets = []
        for ti, tool in enumerate(tools_set):
            ag_datasets.append({"label": tool, "data": [ag_by_level.get(lv,{}).get(tool, 0) for lv in ag_levels], "backgroundColor": ag_colors[ti % 8]})
        chart_scripts += _chart('chart_s05e', 'bar', {"labels": ag_levels, "datasets": ag_datasets}, {"responsive":True, "plugins":{"title":{"display":True,"text":"Tool Agreement Rate by Taxonomic Level"}, "legend":{"position":"bottom","labels":{"boxWidth":12,"font":{"size":10}}}}, "scales":{"y":{"beginAtZero":True,"max":100,"title":{"text":"Agreement Rate (%)"}}}})
        stage_has_chart['s05e'] = True
    cs_rows = _read_tsv(report_dir / "taxonomy_consistency_summary.tsv")
    if cs_rows:
        cs_by_level = {}
        for r in cs_rows:
            lv = r.get("Level","")
            if lv not in cs_by_level: cs_by_level[lv] = {}
            cs_by_level[lv][r.get("Tool","")] = int(r.get("Classified",0))
        cs_levels = [lv for lv in levels_order if lv in cs_by_level]
        cs_tools = sorted(set(r.get("Tool","") for r in cs_rows))
        cs_colors = list(OI[:8])
        cs_datasets = []
        for ti, tool in enumerate(cs_tools):
            cs_datasets.append({"label": tool, "data": [cs_by_level.get(lv,{}).get(tool, 0) for lv in cs_levels], "backgroundColor": cs_colors[ti % 8]})
        chart_scripts += _chart('chart_s05f', 'bar', {"labels": cs_levels, "datasets": cs_datasets}, {"responsive":True, "plugins":{"title":{"display":True,"text":"Tool Classification Count by Level"}, "legend":{"position":"bottom","labels":{"boxWidth":12,"font":{"size":10}}}}, "scales":{"x":{"stacked":True},"y":{"stacked":True,"beginAtZero":True,"title":{"text":"vOTU Count"}}}})
        stage_has_chart['s05f'] = True

    # S05 — Taxonomy novelty
    tax_novelty_kv = {}
    for s in stage_stats:
        if 'Novel Rank' in s['Stage'] and 'Known=' in s.get('Key_Metric',''):
            tax_novelty_kv = _extract_kv(s['Key_Metric'], r'(Known|NewSp|NewGe|NewFa)=(\d+)')
    if tax_novelty_kv:
        stage_has_chart['s05'] = True
        bar_colors = [_sem("novelty", "Known", OI[8]), _sem("novelty", "NewSp", OI[0]),
                      _sem("novelty", "NewGe", OI[5]), _sem("novelty", "NewFa", OI[6])]
        _tax_label_map = {"Known":"Known","NewSp":"Novel species","NewGe":"Novel genus","NewFa":"Novel family"}
        chart_scripts += _chart('chart_s05a', 'bar', {
            "labels": [_tax_label_map.get(k,k) for k in tax_novelty_kv.keys()],
            "datasets": [{"label":"Sequences","data":list(tax_novelty_kv.values()),"backgroundColor":bar_colors}]},
            {"responsive":True,"plugins":{"title":{"display":True,"text":"Taxonomy Novelty"},"legend":{"display":False}},
             "scales":{"y":{"beginAtZero":True,"title":{"text":"Sequences"}}}})

    # S06 — Host distribution
    host_kv = {}
    for s in stage_stats:
        if _D['d_host_pred'] in s['Stage'] and '=' in s.get('Key_Metric',''):
            host_kv = _extract_kv(s['Key_Metric'], r'(\w+)=(\d+)')
    if host_kv:
        stage_has_chart['s06'] = True
        host_colors = list(OI[:8])
        chart_scripts += _chart('chart_s06a', 'bar', {
            "labels": list(host_kv.keys()),
            "datasets": [{"label":"Sequences","data":list(host_kv.values()),"backgroundColor":host_colors[:len(host_kv)]}]},
            {"responsive":True,"plugins":{"title":{"display":True,"text":"Host Prediction Distribution"},"legend":{"display":False}},
             "scales":{"y":{"beginAtZero":True,"title":{"text":"Sequences"}}}})

    # S07 — CheckV quality
    cv_rows = _read_tsv(report_dir / "checkv_summary.tsv")
    if cv_rows:
        qlabels = ["Complete","High-quality","Medium-quality","Low-quality","Not-determined"]
        cv_hosts = [r["Host"][:14] for r in cv_rows if r.get("Host","")!="TOTAL"]
        cv_qdata = {q: [] for q in qlabels}
        for r in cv_rows:
            if r.get("Host","") == "TOTAL": continue
            for q in qlabels: cv_qdata[q].append(int(r.get(q,0)))
        if cv_hosts:
            stage_has_chart['s07a'] = True
            colors = [_sem("quality", q, d) for q, d in
                      (("Complete", OI[4]), ("High-quality", OI[2]), ("Medium-quality", OI[0]),
                       ("Low-quality", OI[5]), ("Not-determined", "#CCCCCC"))]
            cv_datasets = [{"label":q,"data":cv_qdata[q],"backgroundColor":colors[i]} for i,q in enumerate(qlabels)]
            chart_scripts += _chart('chart_s07a', 'bar', {"labels":cv_hosts,"datasets":cv_datasets},
                {"responsive":True,"plugins":{"title":{"display":True,"text":"CheckV Quality by Host"}},
                 "scales":{"x":{"stacked":True},"y":{"stacked":True,"beginAtZero":True,"title":{"text":"Contig Count"}}}})
    cv_conf_rows = _read_tsv(report_dir / "checkv_confidence.tsv")
    if cv_conf_rows:
        cql = [k for k in cv_conf_rows[0] if k not in ("Host","Total")]
        cf_hosts = [r["Host"][:14] for r in cv_conf_rows]
        cf_data = {q: [] for q in cql}
        for r in cv_conf_rows:
            for q in cql: cf_data[q].append(int(r.get(q,0)))
        if cf_hosts and cql:
            stage_has_chart['s07b'] = True
            colors2 = _ramp(5)
            cf_datasets = [{"label":q,"data":cf_data[q],"backgroundColor":colors2[i%5]} for i,q in enumerate(cql)]
            chart_scripts += _chart('chart_s07b', 'bar', {"labels":cf_hosts,"datasets":cf_datasets},
                {"responsive":True,"plugins":{"title":{"display":True,"text":"CheckV aai_confidence by Host"}},
                 "scales":{"x":{"stacked":True},"y":{"stacked":True,"beginAtZero":True,"title":{"text":"Contig Count"}}}})

    # S08 — Rescue branches
    rescue_bc = {}
    rr = report_dir / "rescue_report.tsv"
    if rr.is_file():
        import csv as _rc
        try:
            for row in _rc.DictReader(open(rr, encoding='utf-8'), delimiter='\t'):
                b = (row.get('branch') or '').strip()
                m = (row.get('method') or '').strip()
                if b in ('fail','F','') or b == '-':
                    k = 'Failed'
                else:
                    k = {'A':'CheckV','B':'VSI','C':'BLASTN','D':'genus_len'}.get(b, b)
                    if m and 'vsi' in m.lower(): k = 'VSI'
                rescue_bc[k] = rescue_bc.get(k, 0) + 1
        except Exception:
            rescue_bc = {}
    if rescue_bc:
        stage_has_chart['s08'] = True
        _labels = list(rescue_bc.keys()); _data = list(rescue_bc.values())
        chart_scripts += _chart('chart_s08', 'pie', {
            "labels": _labels,
            "datasets": [{"data":_data,"backgroundColor":[OI[2], OI[1], OI[0], OI[6], OI[5]]}]},
            {"responsive":True,"plugins":{"title":{"display":True,"text":"Rescue Branch Contributions"},"legend":{"position":"bottom"}}})

    # S09 — Viroid species summary
    viroid_sp = report_dir / "Viroid.species_info.tsv"
    if viroid_sp.is_file():
        vr_rows = _read_tsv(viroid_sp)
        if vr_rows:
            top_sp = sorted(vr_rows, key=lambda r: int(r.get('n_contigs',0)), reverse=True)[:20]
            sp_labels = [r['species'][:30] for r in top_sp]
            sp_counts = [int(r['n_contigs']) for r in top_sp]
            stage_has_chart['s09'] = True
            chart_scripts += _chart('chart_s09', 'bar', {
                "labels": sp_labels,
                "datasets": [{"label":"Contigs","data":sp_counts,"backgroundColor":_sem("genome","Viroid", OI[2])}]},
                {"indexAxis":"y","responsive":True,
                 "plugins":{"title":{"display":True,"text":"Top Viroid Species by Contig Count"},"legend":{"display":False}},
                 "scales":{"x":{"title":{"text":"Contig Count"},"beginAtZero":True}}})

    # S09b — Plant virus genus distribution
    plant_genus = report_dir / "all_plant_viruses_genus_summary.tsv"
    if plant_genus.is_file():
        pg_rows = _read_tsv(plant_genus)
        if pg_rows:
            # Top 20 genera by n_contigs, "Unclassified" to end
            pg_sorted = sorted(pg_rows, key=lambda r: (r.get("Genus","") == "Unclassified", -int(r.get("n_contigs",0))))
            top_pg = pg_sorted[:20]
            pg_labels = [r['Genus'][:25] for r in top_pg]
            pg_contigs = [int(r['n_contigs']) for r in top_pg]
            pg_rescued = [int(r['n_rescued']) for r in top_pg]
            stage_has_chart['s09b'] = True
            chart_scripts += _chart('chart_s09b', 'bar', {
                "labels": pg_labels,
                "datasets": [
                    {"label":"Total","data":pg_contigs,"backgroundColor":_C_GOOD},
                    {"label":"Rescued HQ","data":pg_rescued,"backgroundColor":OI[1]}]},
                {"indexAxis":"y","responsive":True,
                 "plugins":{"title":{"display":True,"text":"Plant Virus Genera (All Predicted)"},"legend":{"position":"bottom"}},
                 "scales":{"x":{"title":{"text":"Contig Count"},"beginAtZero":True,"stacked":True}}})

    # S09c — Rescue validation (CDD + blast 五层证据链)
    scored = report_dir / "rescue_evidence_scored.tsv"
    if scored.is_file():
        sc_rows = _read_tsv(scored)
        if sc_rows:
            verdict_cnt = {}
            for r in sc_rows:
                v = r.get('verdict', '')
                verdict_cnt[v] = verdict_cnt.get(v, 0) + 1
            v_order = ['KEEP', 'REVIEW', 'DROP']
            v_labels = [v for v in v_order if v in verdict_cnt]
            v_counts = [verdict_cnt[v] for v in v_labels]
            v_colors = {'KEEP': _sem("verdict","KEEP", OI[2]),
                        'REVIEW': _sem("verdict","REVIEW", OI[0]),
                        'DROP': _sem("verdict","DROP", OI[5])}
            stage_has_chart['s09c'] = True
            chart_scripts += _chart('chart_s09c', 'pie', {
                "labels": v_labels,
                "datasets": [{"data": v_counts, "backgroundColor": [v_colors[v] for v in v_labels]}]},
                {"responsive": True,
                 "plugins": {"title": {"display": True, "text": "Rescue Validation (CDD + BLAST Evidence)"},
                             "legend": {"position": "bottom"}}})
            # novelty 分布 (第二个图 s09d)
            nov_cnt = {}
            for r in sc_rows:
                n = r.get('novelty', '')
                nov_cnt[n] = nov_cnt.get(n, 0) + 1
            nov_labels = list(nov_cnt.keys())
            nov_counts = [nov_cnt[n] for n in nov_labels]
            stage_has_chart['s09d'] = True
            chart_scripts += _chart('chart_s09d', 'bar', {
                "labels": nov_labels,
                "datasets": [{"label": "Contigs", "data": nov_counts, "backgroundColor": OI[4]}]},
                {"indexAxis": "y", "responsive": True,
                 "plugins": {"title": {"display": True, "text": "Novelty Annotation"},
                             "legend": {"display": False}},
                 "scales": {"x": {"title": {"text": "Contig Count"}, "beginAtZero": True}}})

    # 检出频次分布 (chart_s09e)
            freq_bins = {'0 sample': 0, '1 sample': 0, '2-4 samples': 0, '5+ samples': 0}
            for r in sc_rows:
                try:
                    ns = int(float(r.get('n_samples') or 0))
                except Exception:
                    ns = 0
                if ns <= 0: freq_bins['0 sample'] += 1
                elif ns == 1: freq_bins['1 sample'] += 1
                elif ns <= 4: freq_bins['2-4 samples'] += 1
                else: freq_bins['5+ samples'] += 1
            stage_has_chart['s09e'] = True
            chart_scripts += _chart('chart_s09e', 'bar', {
                "labels": list(freq_bins.keys()),
                "datasets": [{"label": "Candidates", "data": list(freq_bins.values()), "backgroundColor": _C_GOOD}]},
                {"responsive": True,
                 "plugins": {"title": {"display": True, "text": "Detection Frequency"},
                             "legend": {"display": False}},
                 "scales": {"x": {"title": {"text": "Samples Detected"}, "beginAtZero": True}}})

            # ── S09f: Total Score Distribution (加权 30/30/40 评分分箱) ──
            def _binscore(v):
                try:
                    return float(v)
                except Exception:
                    return 0.0
            sc_bins = {'0-20': 0, '20-40': 0, '40-60': 0, '60-80': 0, '80-100': 0}
            for r in sc_rows:
                ts = _binscore(r.get('total_score', ''))
                if ts < 20: sc_bins['0-20'] += 1
                elif ts < 40: sc_bins['20-40'] += 1
                elif ts < 60: sc_bins['40-60'] += 1
                elif ts < 80: sc_bins['60-80'] += 1
                else: sc_bins['80-100'] += 1
            stage_has_chart['s09f'] = True
            chart_scripts += _chart('chart_s09f', 'bar', {
                "labels": list(sc_bins.keys()),
                "datasets": [{"label": "Candidates", "data": list(sc_bins.values()), "backgroundColor": _C_WARN}]},
                {"responsive": True,
                 "plugins": {"title": {"display": True, "text": "Total Evidence Score (30% nt + 30% aa + 40% Domain: max(CDD, HMM))"},
                             "legend": {"display": False}},
                 "scales": {"x": {"title": {"text": "Weighted Score"}, "beginAtZero": True}}})

            # ── S09g: Category Distribution (final_judgement 病毒属性) ──
            cat_cnt = {}
            for r in sc_rows:
                c = r.get('category', '') or '未定'
                cat_cnt[c] = cat_cnt.get(c, 0) + 1
            cat_order = ['已知植物病毒', '新种候选', '远缘(未定)', '近缘(未定)', '宿主基因污染',
                         '噬菌体污染', '昆虫病毒/污染', '真菌病毒/污染', '短片段/低覆盖', '未定']
            cat_labels = [c for c in cat_order if cat_cnt.get(c, 0) > 0]
            cat_vals = [cat_cnt[c] for c in cat_labels]
            stage_has_chart['s09g'] = True
            chart_scripts += _chart('chart_s09g', 'bar', {
                "labels": cat_labels,
                "datasets": [{"label": "Contigs", "data": cat_vals, "backgroundColor": OI[4]}]},
                {"indexAxis": "y", "responsive": True,
                 "plugins": {"title": {"display": True, "text": "Virus Attribute (Category)"},
                             "legend": {"display": False}},
                 "scales": {"x": {"title": {"text": "Contig Count"}, "beginAtZero": True}}})

            # ── S09h: Score Composition (top KEEP 候选的 blastn/blastx/CDD 构成) ──
            keep_scored = [r for r in sc_rows if r.get('verdict') == 'KEEP']
            keep_scored.sort(key=lambda r: _binscore(r.get('total_score', '')), reverse=True)
            topk = keep_scored[:10]
            if topk:
                lbls = []
                bblast, bblastx, bcdd, bhmm = [], [], [], []
                for r in topk:
                    cid = r.get('contig_id', '')
                    short = cid.split('_NODE_')[-1][:18] if '_NODE_' in cid else cid[:18]
                    lbls.append(short)
                    bblast.append(round(_binscore(r.get('blastn_score', '')), 1))
                    bblastx.append(round(_binscore(r.get('blastx_score', '')), 1))
                    # 域证据层 = max(CDD, HMM): CDD 段 + HMM 补充段 (hmm>cdd 时的高出差值) 恰好拼出 domain
                    _cdd = _binscore(r.get('cdd_score', ''))
                    bcdd.append(round(_cdd, 1))
                    bhmm.append(round(max(_binscore(r.get('domain_score', '')) - _cdd, 0.0), 1))
                stage_has_chart['s09h'] = True
                chart_scripts += _chart('chart_s09h', 'bar', {
                    "labels": lbls,
                    "datasets": [
                        {"label": "blastn (30%)", "data": bblast, "backgroundColor": _ramp(4)[0]},
                        {"label": "blastx (30%)", "data": bblastx, "backgroundColor": _ramp(4)[1]},
                        {"label": "CDD (domain layer)", "data": bcdd, "backgroundColor": _ramp(4)[2]},
                        {"label": "HMM complement (max part)", "data": bhmm, "backgroundColor": _ramp(4)[3]}]},
                    {"responsive": True,
                     "plugins": {"title": {"display": True, "text": "Top KEEP: Score Composition"},
                                 "legend": {"position": "bottom"}},
                     "scales": {"x": {"stacked": True, "title": {"text": "Candidate"}, "ticks": {"maxRotation": 45}},
                                "y": {"stacked": True, "title": {"text": "Score"}, "beginAtZero": True}}})

            # ── S09i: Removal/Retention Rationale (去假病毒判定依据) ──
            def _sf(v):
                try: return float(v) if v not in ('', None) else 0.0
                except Exception: return 0.0
            reason_keep, reason_drop = {}, {}
            for r in sc_rows:
                v = r.get('verdict', ''); cdd = r.get('cdd_evidence', '')
                hr = r.get('hmm_rescue', '') or ''
                nt = _sf(r.get('nt_pident')); aa = _sf(r.get('aa_pident')); ns = _sf(r.get('n_samples'))
                if v == 'KEEP':
                    if hr.endswith('->KEEP'): lab = 'CT3-HMM full-length profile'
                    elif cdd in ('PASS_CORE', 'PASS_CORE_FAM', 'PASS_VIRAL'): lab = 'Viral CDD domain'
                    elif nt >= 85: lab = 'BLAST-known (nt≥85)'
                    elif nt >= 40 or aa >= 40: lab = 'BLAST/protein moderate hit (40%+)'
                    elif ns > 0: lab = 'Detected but weak evidence'
                    else: lab = 'Other'
                    reason_keep[lab] = reason_keep.get(lab, 0) + 1
                elif v == 'DROP':
                    if cdd in ('PASS_CORE', 'PASS_CORE_FAM', 'PASS_VIRAL'): lab = '⚠ Strong CDD but removed'
                    elif nt >= 85: lab = 'BLAST-known but removed'
                    elif aa >= 40: lab = 'Protein-approx but removed'
                    elif ns > 0: lab = 'Detected but removed'
                    elif cdd in ('AMBIGUOUS', 'REVIEW'): lab = 'No viral domain + weak BLAST'
                    else: lab = 'No evidence'
                    reason_drop[lab] = reason_drop.get(lab, 0) + 1
            # CT3-HMM 第三路探针补票 (无 CDD 病毒域但全长 profile 命中, DROP→REVIEW)
            reason_hmm = {}
            for r in sc_rows:
                _hr = (r.get('hmm_rescue', '') or '').strip()
                if _hr.endswith('->KEEP') or _hr == 'DROP->REVIEW':
                    reason_hmm['HMM rescue (CT3)'] = reason_hmm.get('HMM rescue (CT3)', 0) + 1
            all_labs = []
            for k in list(reason_drop.keys()) + list(reason_keep.keys()) + list(reason_hmm.keys()):
                if k not in all_labs: all_labs.append(k)
            # 按 DROP 理由优先排序(去假病毒视角), 无则该类放后续
            drop_order = [k for k in all_labs if k in reason_drop]
            drop_order.sort(key=lambda k: -reason_drop[k])
            tail = [k for k in all_labs if k not in drop_order]
            labels_oi = drop_order + tail
            stage_has_chart['s09i'] = True
            chart_scripts += _chart('chart_s09i', 'bar', {
                "labels": labels_oi,
                "datasets": [
                    {"label": "DROP", "data": [reason_drop.get(k, 0) for k in labels_oi], "backgroundColor": _C_BAD},
                    {"label": "KEEP", "data": [reason_keep.get(k, 0) for k in labels_oi], "backgroundColor": _C_GOOD},
                    {"label": "HMM rescue", "data": [reason_hmm.get(k, 0) for k in labels_oi], "backgroundColor": OI[1]}]},
                {"responsive": True,
                 "plugins": {"title": {"display": True, "text": "Removal / Retention Rationale"},
                             "legend": {"position": "bottom"}},
                 "scales": {"x": {"title": {"text": "Candidates"}, "stacked": True, "beginAtZero": True},
                             "y": {"stacked": True, "beginAtZero": True, "title": {"text": "Count"}}}})

    # S11a — 综合新病毒判定 (9a × ACVirus 整合) — keep_summary final_call
    # (原位置在 s09 链中间, 会把 s09e~s09i 的 12 空格块错误并入 `if ks_rows:`
    #  且 scored 表缺失时 sc_rows 未绑定 → 崩溃; 移到链后独立成块)
    ksum = _resolve_acv_dir(report_dir.parent) / "keep_summary.tsv"
    if ksum.is_file():
        ks_rows = _read_tsv(ksum)
        if ks_rows:
            fc_cnt = {}
            for r in ks_rows:
                k = r.get('final_call', '') or 'Unknown'
                fc_cnt[k] = fc_cnt.get(k, 0) + 1
            order = ['Novel(dual)', 'Novel(blast)', 'Novel(ACV)', 'Known(ACV-conflict)', 'Known']
            labels = [k for k in order if k in fc_cnt]
            counts = [fc_cnt[k] for k in labels]
            colors = {'Novel(dual)': OI[6], 'Novel(blast)': OI[5], 'Novel(ACV)': OI[0],
                      'Known(ACV-conflict)': OI[1], 'Known': OI[8]}
            stage_has_chart['s11a'] = True
            chart_scripts += _chart('chart_s11a', 'pie', {
                "labels": labels,
                "datasets": [{"data": counts, "backgroundColor": [colors.get(k,'#ccc') for k in labels]}]},
                {"responsive": True,
                 "plugins": {"title": {"display": True,
                                       "text": "Comprehensive Novelty Call (9a BLAST × ACVirus)"},
                             "legend": {"position": "bottom"}}})

    # S10 — 09 Virome Analysis: 核酸类型 (DNA/RNA) + 科属分布, 类病毒单列不并入 RNA
    _s10_note = ""
    _TYPE_COLORS = {"DNA": _sem("genome","DNA", OI[4]), "RNA": _sem("genome","RNA", OI[5]),
                    "Viroid": _sem("genome","Viroid", OI[2]),
                    "Unassigned": _sem("genome","Unassigned", OI[8])}
    _TYPE_ORDER = ["DNA", "RNA", "Viroid", "Unassigned"]
    _VIROID_FAMS = ("Pospiviroidae", "Avsunviroidae")

    def _na_bucket(fam, gen, sp, gt):
        blob = ("%s %s %s" % (fam, gen, sp)).lower()
        if fam in _VIROID_FAMS or "viroid" in blob:
            return "Viroid"
        return gt if gt in ("DNA", "RNA") else "Unassigned"

    def _lab(x):
        x = (x or '').strip()
        return x if x and x.upper() not in ('NA', 'NAN') else '(Unassigned)'

    hq_info = report_dir / "HQ_plant_viruses_info.tsv"
    if hq_info.is_file():
        hq_rows = _read_tsv(hq_info)
        if hq_rows:
            _type_cnt = {}      # {类型: contig 数}
            _fam_by_type = {}   # {类型: {科: n}}
            _gen_by_type = {}   # {类型: {属: n}}
            for r in hq_rows:
                fam = (r.get('Family','') or '').strip()
                gen = (r.get('Genus','') or '').strip()
                sp  = (r.get('Species','') or '').strip()
                # 生产端列表头为 Nucleic_acid (virome_pipeline.py HQ 表 /
                # build_plant_virus_info.py); 兼容旧 Genome_Type 列名回退
                gt  = (r.get('Nucleic_acid') or r.get('Genome_Type') or '').strip()
                b = _na_bucket(fam, gen, sp, gt)
                _type_cnt[b] = _type_cnt.get(b, 0) + 1
                _fam_by_type.setdefault(b, {})[_lab(fam)] = _fam_by_type.setdefault(b, {}).get(_lab(fam), 0) + 1
                _gen_by_type.setdefault(b, {})[_lab(gen)] = _gen_by_type.setdefault(b, {}).get(_lab(gen), 0) + 1

            # 类病毒: 独立数据源 (viroid_analysis/, 与 HQ 病毒表零重叠), 单列一类
            _viroid_seen, _viroid_fam, _viroid_gen = set(), {}, {}
            _vip = report_dir / "Viroid.all_info.tsv"
            if _vip.is_file():
                for r in _read_tsv(_vip):
                    cid = (r.get('contig_id','') or '').strip()
                    if not cid or cid in _viroid_seen:
                        continue
                    _viroid_seen.add(cid)
                    vf = (r.get('family','') or '').strip() or 'Unclassified viroid'
                    vg = (r.get('genus','') or '').strip() or '(Unassigned)'
                    _viroid_fam[vf] = _viroid_fam.get(vf, 0) + 1
                    _viroid_gen[vg] = _viroid_gen.get(vg, 0) + 1
            _n_viroid = len(_viroid_seen)
            if _n_viroid:
                _type_cnt['Viroid'] = _type_cnt.get('Viroid', 0) + _n_viroid
                for _k, _v in _viroid_fam.items():
                    _fam_by_type.setdefault('Viroid', {})[_k] = _fam_by_type.setdefault('Viroid', {}).get(_k, 0) + _v
                for _k, _v in _viroid_gen.items():
                    _gen_by_type.setdefault('Viroid', {})[_k] = _gen_by_type.setdefault('Viroid', {}).get(_k, 0) + _v

            _n_hq = len(hq_rows)
            _n_all = sum(_type_cnt.get(k, 0) for k in _TYPE_ORDER) or 1

            # ① 核酸类型分布 (环形图) — DNA / RNA / Viroid / Unassigned
            _np = [(k, _type_cnt[k]) for k in _TYPE_ORDER if _type_cnt.get(k)]
            if _np:
                stage_has_chart['s10na'] = True
                chart_scripts += _chart('chart_s10na', 'doughnut', {
                    "labels": [k for k, _ in _np],
                    "datasets": [{"data": [v for _, v in _np],
                                  "backgroundColor": [_TYPE_COLORS[k] for k, _ in _np]}]},
                    {"responsive": True, "cutout": "50%",
                     "plugins": {"title": {"display": True,
                                           "text": "Nucleic Acid Type Distribution (HQ n=%d + viroid n=%d)" % (_n_hq, _n_viroid)},
                                 "legend": {"position": "bottom"}}})

            # ② 科分布 (top 15, 按核酸类型堆叠)
            _fam_total = {}
            for _b, _d in _fam_by_type.items():
                for _k, _v in _d.items():
                    _fam_total[_k] = _fam_total.get(_k, 0) + _v
            _top_fam = sorted(_fam_total.items(), key=lambda x: -x[1])[:15]
            if _top_fam:
                _fkeys = [k for k, _ in _top_fam]
                _fds = []
                for _b in _TYPE_ORDER:
                    _d = _fam_by_type.get(_b, {})
                    if any(_d.get(k) for k in _fkeys):
                        _fds.append({"label": _b, "data": [_d.get(k, 0) for k in _fkeys],
                                     "backgroundColor": _TYPE_COLORS[_b]})
                stage_has_chart['s10'] = True
                chart_scripts += _chart('chart_s10a', 'bar', {
                    "labels": [k[:28] for k in _fkeys], "datasets": _fds},
                    {"indexAxis": "y", "responsive": True,
                     "plugins": {"title": {"display": True,
                                           "text": "HQ Virus Family Distribution by Nucleic Acid Type"},
                                 "legend": {"position": "bottom"}},
                     "scales": {"x": {"stacked": True, "beginAtZero": True, "title": {"text": "Contig Count"}},
                                "y": {"stacked": True}}})

            # ③ 属分布 (top 15, 按核酸类型堆叠)
            _gen_total = {}
            for _b, _d in _gen_by_type.items():
                for _k, _v in _d.items():
                    _gen_total[_k] = _gen_total.get(_k, 0) + _v
            _top_gen = sorted(_gen_total.items(), key=lambda x: -x[1])[:15]
            if _top_gen:
                _gkeys = [k for k, _ in _top_gen]
                _gds = []
                for _b in _TYPE_ORDER:
                    _d = _gen_by_type.get(_b, {})
                    if any(_d.get(k) for k in _gkeys):
                        _gds.append({"label": _b, "data": [_d.get(k, 0) for k in _gkeys],
                                     "backgroundColor": _TYPE_COLORS[_b]})
                stage_has_chart['s10b'] = True
                chart_scripts += _chart('chart_s10b', 'bar', {
                    "labels": [k[:28] for k in _gkeys], "datasets": _gds},
                    {"indexAxis": "y", "responsive": True,
                     "plugins": {"title": {"display": True,
                                           "text": "HQ Virus Genus Distribution by Nucleic Acid Type"},
                                 "legend": {"position": "bottom"}},
                     "scales": {"x": {"stacked": True, "beginAtZero": True, "title": {"text": "Contig Count"}},
                                "y": {"stacked": True}}})

            # ④ 落盘: 类型汇总 + 科×类型 + 属×类型 (报告内嵌 + 可下载)
            try:
                with open(report_dir / "hq_genome_type_summary.tsv", "w", encoding="utf-8", newline="") as _f:
                    _w = csv.writer(_f, delimiter="\t")
                    _w.writerow(["Category", "n_contigs", "pct"])
                    for _k in _TYPE_ORDER:
                        if _type_cnt.get(_k):
                            _w.writerow([_k, _type_cnt[_k], round(100.0 * _type_cnt[_k] / _n_all, 2)])
                    _w.writerow(["Total", _n_all, 100.0])
                with open(report_dir / "hq_genome_type_family.tsv", "w", encoding="utf-8", newline="") as _f:
                    _w = csv.writer(_f, delimiter="\t")
                    _w.writerow(["Family", "DNA", "RNA", "Viroid", "Unassigned", "Total"])
                    for _k, _t in sorted(_fam_total.items(), key=lambda x: -x[1]):
                        _w.writerow([_k] + [_fam_by_type.get(b, {}).get(_k, 0) for b in ("DNA", "RNA", "Viroid", "Unassigned")] + [_t])
                with open(report_dir / "hq_genome_type_genus.tsv", "w", encoding="utf-8", newline="") as _f:
                    _w = csv.writer(_f, delimiter="\t")
                    _w.writerow(["Genus", "DNA", "RNA", "Viroid", "Total"])
                    for _k, _t in sorted(_gen_total.items(), key=lambda x: -x[1]):
                        _w.writerow([_k] + [_gen_by_type.get(b, {}).get(_k, 0) for b in ("DNA", "RNA", "Viroid")] + [_t])
            except Exception:
                pass

            _conflict = sum(1 for _k in _fam_total
                            if _fam_by_type.get('DNA', {}).get(_k, 0) > 0 and _fam_by_type.get('RNA', {}).get(_k, 0) > 0)
            _s10_note = ('<div class="chart-box" style="grid-column:1/-1">'
                         '<div class="chart-title">Note on nucleic acid type</div>'
                         '<div style="font-size:12px;line-height:1.7;color:var(--muted,#667)">'
                         'Type is resolved from VMR by precedence species &gt; genus &gt; family &gt; Realm/Kingdom '
                         '(column <code>Nucleic_acid</code>). Viroids are counted from <code>viroid_analysis/</code> and reported '
                         'as a separate class, never merged into RNA. %d families here carry a family label that disagrees '
                         'with the resolved type (mixed-tool taxonomy integration); for those rows the type, being the more '
                         'specific call, is the reliable one.'
                         '</div></div>') % _conflict

    # Sankey 交互式嵌入 (用 Blob URL 动态注入, 避免 data URI 大小限制)
    # 两张 Sankey 与旭日图都是分类学图 → 全挂 05 Taxonomy 章节
    sankey_by_stage = {}  # {stage_key: html_string}
    # 堆叠图默认切换为百分比模式 (延迟执行, 等图表全部渲染完)
    chart_scripts += "setTimeout(function(){document.querySelectorAll('.pct-btn').forEach(function(b){toggleStackedPct(b)})},100);\n"
    sankey_inject_scripts = ""
    sankey_map = [("classification_sankey.html","Taxonomy Classification Sankey","s05"),
                  ("classification_sankey_plant.html","Plant Virus Taxonomy Sankey","s05")]
    for i, (sname, stitle, stage_key) in enumerate(sankey_map):
        spath = report_dir / sname
        if spath.is_file():
            import base64
            with open(spath, "rb") as sf:
                sankey_b64 = base64.b64encode(sf.read()).decode()
            card = f'''<div class="sankey-card">
<h3>{stitle}</h3>
<iframe id="sankey_iframe_{i}" style="width:100%;height:480px;border:none;border-radius:4px" loading="lazy"></iframe>
</div>\n'''
            sankey_by_stage.setdefault(stage_key, "")
            sankey_by_stage[stage_key] += card
            sankey_inject_scripts += f"(function(){{var b='{sankey_b64}';var d=decodeURIComponent(escape(atob(b)));var u=URL.createObjectURL(new Blob([d],{{type:'text/html;charset=utf-8'}}));document.getElementById('sankey_iframe_{i}').src=u;}})();\n"

    # 旭日图 (Plotly Sunburst) → s05 section
    sunburst_path = report_dir / "taxonomy_sunburst.html"
    if sunburst_path.is_file():
        import base64
        with open(sunburst_path, "rb") as sf:
            sunburst_b64 = base64.b64encode(sf.read()).decode()
        sunburst_card = f'''<div class="sankey-card">
<h3>Plant Virus Taxonomy Sunburst</h3>
<iframe id="sunburst_iframe" style="width:100%;height:750px;border:none;border-radius:4px" loading="lazy"></iframe>
</div>\n'''
        sankey_by_stage.setdefault("s05", "")
        sankey_by_stage["s05"] += sunburst_card
        sankey_inject_scripts += f"(function(){{var b='{sunburst_b64}';var d=decodeURIComponent(escape(atob(b)));var u=URL.createObjectURL(new Blob([d],{{type:'text/html;charset=utf-8'}}));document.getElementById('sunburst_iframe').src=u;}})();\n"

    # ── KPI ──
    kpis = {}
    for s in stage_stats:
        if s['Stage'] in (_D['d_clean'], '  └ data_summary'):
            for m in re.finditer(r'reads:\s*([\d,]+)→([\d,]+)', s.get('Key_Metric','')):
                kpis['raw_reads'] = m.group(1); kpis['clean_reads'] = m.group(2)
            for m in re.finditer(r'bases:\s*([\d,]+)→([\d,]+)', s.get('Key_Metric','')):
                kpis['raw_bases'] = m.group(1); kpis['clean_bases'] = m.group(2)
        if s['Stage'] == _D['d_clean']:
            for m in re.finditer(r'(\d+)\s*samples', s.get('Key_Metric','')): kpis['n_sample'] = m.group(1)
        if s['Stage'] == _D['d_asm']:
            for m in re.finditer(r'([\d,]+)\s*contigs', s.get('Key_Metric','')): kpis['total_contigs'] = m.group(1)
            for m in re.finditer(r'([\d.]+)\s*Mb', s.get('Key_Metric','')): kpis['total_mb'] = m.group(1)
        if s['Stage'] == '02_Identification':
            for m in re.finditer(r'([\d,]+)\s*viral sequences', s.get('Key_Metric','')): kpis['virus_seqs'] = m.group(1)
        if 'Novel Rank' in s['Stage']:
            kpis['novelty'] = s.get('Key_Metric','')
            for m in re.finditer(r'Known=(\d+)', s.get('Key_Metric','')): kpis['n_known'] = m.group(1)
            for m in re.finditer(r'NewSp=(\d+)', s.get('Key_Metric','')): kpis['n_newsp'] = m.group(1)
            for m in re.finditer(r'NewGe=(\d+)', s.get('Key_Metric','')): kpis['n_newge'] = m.group(1)
            for m in re.finditer(r'NewFa=(\d+)', s.get('Key_Metric','')): kpis['n_newfa'] = m.group(1)
        if '  └ vclust' in s['Stage']:
            for m in re.finditer(r'([\d,]+)\s*clusters', s.get('Key_Metric','')): kpis['n_clusters'] = m.group(1)
        if s['Stage'] == _D['d_host_pred']:      # 精确匹配 (旧 startswith('06') 在 standard 布局失效)
            for m in re.finditer(r'([\d,]+)\s*seqs', s.get('Key_Metric','')): kpis['host_total'] = m.group(1)
        if s['Stage'] == '07_CheckV':
            for m in re.finditer(r'(\d+)\s*HQ', s.get('Key_Metric','')): kpis['hq_votus'] = m.group(1)
            for m in re.finditer(r'(\d+)\s*total', s.get('Key_Metric','')): kpis['cv_total'] = m.group(1)
        if s['Stage'] == _D['d_rescue']:
            for m in re.finditer(r'([\d,]+)\s*viruses', s.get('Key_Metric','')): kpis['rescued'] = m.group(1)
    # 从 hostdep_summary.tsv 提取宿主去除统计, 用量化的数据量(Gb/Mb)
    hd_rows = _read_tsv(report_dir / "hostdep_summary.tsv")
    if hd_rows:
        hd_raw_total = sum(int(r.get("Raw",0)) for r in hd_rows)
        hd_after_total = sum(int(r.get("After_Host",0)) for r in hd_rows)
        if hd_raw_total > 0:
            hd_pct = hd_after_total / hd_raw_total * 100
            kpis['hd_retained'] = f"{hd_pct:.1f}" if hd_pct >= 1 else f"{hd_pct:.2f}"
            # 用 QC 平均 read 长度换算 bases
            try:
                avg_len = int(kpis.get('raw_bases','0').replace(',','')) / max(int(kpis.get('raw_reads','1').replace(',','')), 1)
            except: avg_len = 150
            kpis['hd_raw_bp'] = hd_raw_total * avg_len
            kpis['hd_after_bp'] = hd_after_total * avg_len

    # ── Stage sections ──
    stage_defs = [
        ("s00a","CleanData","00a Data Preprocessing","QC","#607d8b"),
        ("s00b","HostDep","00b Host Depletion","DEP","#546e7a"),
        ("s00c","BBNorm","00c BBNorm","BBN","#00838f"),
        ("s01","Assembly","01 Assembly","ASM","#0072B2"),
        ("s02","Ident","02 Identification","ID","#537D96"),
        ("s03","COBRA","03 COBRA Extension","COBRA","#007A59"),
        ("s04","Cluster","04 Clustering","CLU","#A84400"),
        ("s05","Taxonomy","05 Taxonomy Classification","TAX","#9B4B72"),
        ("s06","Host","06 Host Prediction","HOST","#8B2C1F"),
        ("s07","CheckV","07 CheckV Quality","CV","#4A6B4A"),
        ("s08","Rescue","08 Rescue","RESCUE","#1B365D"),
        # ("s09","PlantVirus") removed,
        ("s10","Analysis","09 Virome Analysis","AN","#9B4B72"),
        ("s11","Verify","09b Analysis & Validation","AV","#537D96"),
    ]

    # ── 阶段 → 章节 精确映射 ───────────────────────────────────────────
    # Stage 字段 = 收集函数写入的字面量 (多数是当前布局的目录名, 少数写死)。
    # 旧实现按"旧编号前缀"匹配且不 break → standard 布局大面积串章
    # (01_CleanData 被当成 s01 Assembly, 05/07/08 张冠李戴 等)。
    # 现改为: d_* 逻辑键 → 两布局真实目录名 + 写死字面量别名 → 章节 key 的
    # 精确表; 一次匹配一次赋值 (首条为准, 主阶段行先于子行写入)。
    _stage_section = {}
    for _sk, _dkeys in (
        ('s00a', ('d_clean',)),    ('s00b', ('d_hostdep',)), ('s00c', ('d_bbnorm',)),
        ('s01',  ('d_asm',)),      ('s02',  ('d_ident', 'd_graph')),
        ('s03',  ('d_cobra',)),    ('s04',  ('d_cluster',)),
        ('s05',  ('d_taxonomy',)), ('s06',  ('d_host_pred',)),
        ('s07',  ('d_checkv',)),   ('s08',  ('d_rescue',)),
        ('s10',  ('d_analysis',)), ('s11',  ('d_verify',)),
    ):
        for _dk in _dkeys:
            for _lay in ('standard', 'legacy'):        # 两布局目录名均为候选
                _stage_section[_layout_dirs(_lay)[_dk]] = _sk
    _stage_section.update({                            # 收集函数写死的 Stage 字面量
        '00c_BBNorm': 's00c', '02_Identification': 's02', '03_COBRA': 's03',
        '03b_Merge': 's03', '07_CheckV': 's07', '10_Virome_Analysis': 's10',
    })

    def _stage_key_of(sn):
        """Stage 字面量 → 章节 key (09b_* 归 s11; 未收录者按旧编号回退)"""
        if sn in _stage_section:
            return _stage_section[sn]
        if sn.startswith('09b'):
            return 's11'
        return 's' + (sn.split('_')[0] if '_' in sn else sn)

    def _tbl_id(sk, tsv_name):
        """表容器 id 规则 (唯一真源, 与通用渲染循环共用): '.' → '_'。
        大表 (>500 行) 实际 id 追加 _large, 由渲染循环补同名空 span 别名锚。"""
        return f"stagetbl_{sk}_{tsv_name.replace(chr(46),chr(95))}"

    _subnav = {
        # _subnav 硬编码锚点必须与实际生成的 id 一致 (旧写法缺 _tsv 后缀 → 死锚),
        # 现与渲染循环共用 _tbl_id 规则: 大表也由别名锚保证恒命中
        's02': [('#' + _tbl_id('s02','ident_summary.tsv'),'02a Identification'),
                ('#' + _tbl_id('s02','filter_summary.tsv'),'02b Filter'),
                ('#' + _tbl_id('s02','graph_verdict_summary.tsv'),'graph verdicts')],
        's03': [('#' + _tbl_id('s03','cobra_summary.tsv'),'03a COBRA'),('#merge-flye','03b Merge')],
        's05': [('#upset-container','Tool UpSet')],
        's11': [('#keep-table','KEEP Sequences'),('#s11-tree','Phylogeny'),('#s11-sdt','SDT Matrix'),('#s11-cv','Coverage')],
    }

    stage_status = {}; stage_metric = {}
    for s in stage_stats:
        sk = _stage_key_of(s['Stage'])
        if sk in stage_status: continue                # 一次匹配一次赋值 (首条为准)
        stage_status[sk] = S.get(s['Status'], 'skip')
        stage_metric[sk] = s.get('Key_Metric', '')

    chart_map = {
        's00a': [('chart_s00a','Read Quality'),('chart_s00a_dup','Duplication Rate')],
        's00b': [('chart_s00b','Host Depletion')],
        's01':  [('chart_s01b','Contig Count'),('chart_s01a','N50 (kb)')],
        's02':  [('chart_s02','UniProt Filter')],
        's03':  [('chart_s03','COBRA Rates')],
        's04':  [('chart_s04a','Cluster Size'),('chart_s04c','Pipeline Reduction')],
        's05':  [('chart_s05a','Taxonomy Novelty'),('chart_s05b','Family Distribution'),('chart_s05c','Genus Distribution'),
                  ('chart_s05d','Fill Rate'),('chart_s05e','Agreement Rate'),('chart_s05f','Tool Classification')],
        's06':  [('chart_s06a','Host Distribution')],
        's07':  [('chart_s07a','CheckV Quality'),('chart_s07b','CheckV Confidence')],
        's08':  [('chart_s08','Rescue Branches')],
        # 's09' 主键在 stage_defs 中不存在(已删除) → 09 阶段图全部挂在 s10
        's10':  [('chart_s10na','Nucleic Acid Type'),('chart_s10a','Family Distribution'),('chart_s10b','Genus Distribution'),
                  ('chart_s09','Viroid Species Summary'), ('chart_s09b','Plant Virus Genera'),
                  ('chart_s09c','Rescue Validation'), ('chart_s09d','Novelty Annotation'),
                  ('chart_s09e','Detection Frequency'), ('chart_s09f','Total Evidence Score'),
                  ('chart_s09g','Virus Attribute'), ('chart_s09h','Top KEEP Score Composition'),
                  ('chart_s09i','Removal Rationale')],
        's11':  [('chart_s11a','Comprehensive Novelty Call')],
    }

    # 每个阶段对应的 TSV 文件, 用于在卡片内嵌入数据表
    # ⚠ 键必须存在于 stage_defs: 旧 's09' 键在 s09 章节移除后成死键 ——
    #   hmm_ct3_hits.tsv / hmm_ct3_rescue_candidates.tsv 被复制进报告目录却
    #   永不渲染, 's09' 其余条目已由 s08/s10 覆盖。现移入 s11 (09b) 章节。
    stage_tsv_map = {
        's00a': ['data_summary.tsv'],
        's00b': ['hostdep_summary.tsv'],
        's01':  ['assembly_summary.tsv'],
        's02':  ['ident_summary.tsv', 'filter_summary.tsv', 'graph_verdict_summary.tsv'],
        's03':  ['cobra_summary.tsv'],
        's04':  ['cluster_size_distribution.tsv', 'cluster_pipeline_reduction.tsv'],
        's05':  ['taxonomy_composition.tsv', 'final_integrated_classification.tsv'],
        's06':  ['host_distribution.tsv', 'host_decision_method.tsv', 'host_ictv_classification_summary.tsv', 'host_ictv_confidence.tsv'],
        's07':  ['checkv_summary.tsv', 'checkv_confidence.tsv'],
        's08':  ['rescue_report.tsv', 'plant_virus_summary.tsv'],
        's10':  ['hq_genome_type_summary.tsv', 'hq_genome_type_family.tsv', 'hq_genome_type_genus.tsv',
                  'HQ_plant_viruses_info.tsv', 'rescue_evidence_scored.tsv',
                  'frequency_table.tsv', 'All_plant.viruses_info.tsv',
                  'Viroid.species_info.tsv', 'all_plant_viruses_genus_summary.tsv'],
        # 09b 验证章节: ACVirus 重分类/新病毒鉴定 + CT3-HMM 探针证据
        # (通用渲染循环只认 report_dir 内的相对文件名 → 前两份由
        #  _collect_analysis_files 的 copy_map 一并复制进报告目录)
        's11':  ['final_result_with_confidence.tsv', 'identification_result_score.csv',
                 'hmm_ct3_hits.tsv', 'hmm_ct3_rescue_candidates.tsv'],
    }

    _sd = {}
    _stage_key_map = {sk: short for sk, short, _, _, _ in stage_defs}
    for _s in stage_stats:
        _sn = _s.get("Stage", "")
        if not _sn.startswith("  "):
            _key = _stage_key_of(_sn)                  # 与章节状态同一精确映射 (旧按编号切分)
            _sd[_key] = {"stage": _sn, "status": _s.get("Status",""), "metrics": _s.get("Key_Metric",""), "short": _stage_key_map.get(_key, "")}
    import json as _json4
    stage_data_json = _json4.dumps(_sd, ensure_ascii=False)
    sections_html = ""
    table_pager_js = ""
    stage_extra_html = {}
    if _s10_note:
        stage_extra_html['s10'] = _s10_note

    # ── 分组查表: 科/属/contig → 核酸类型 (DNA/RNA/Viroid/Other) ──
    # 供 s11 进化树与 SDT 分组、s11b KEEP 表核酸列共用。
    # 查表源: All_plant.viruses_info.tsv (全量) → HQ_plant_viruses_info.tsv (补充),
    # 列名兼容新 Nucleic_acid / 旧 Genome_Type, 同名多行取多数派;
    # 反转录科按 VMR 规则兜底 (与 utils/annotate_nucleic_acid.py 的判定一致)。
    def _s11_na_lookup():
        vf, vg, cid_na = {}, {}, {}
        for name in ("All_plant.viruses_info.tsv", "HQ_plant_viruses_info.tsv"):
            p = report_dir / name
            if not p.is_file(): continue
            try: rows = _read_tsv(p)
            except Exception: continue
            for r in rows:
                f = (r.get('Family','') or '').strip()
                g = (r.get('Genus','') or '').strip()
                t = (r.get('Nucleic_acid') or r.get('Genome_Type') or '').strip().upper()
                if t not in ("DNA", "RNA"): continue
                cid = (r.get('contig_id','') or '').strip()
                if cid and cid not in cid_na:
                    cid_na[cid] = t
                if f:
                    vf.setdefault(f, {"DNA": 0, "RNA": 0})[t] += 1
                if g:
                    vg.setdefault(g, {"DNA": 0, "RNA": 0})[t] += 1
        fm = {k: max(v.items(), key=lambda kv: kv[1])[0] for k, v in vf.items()}
        gm = {k: max(v.items(), key=lambda kv: kv[1])[0] for k, v in vg.items()}
        for f in ("Metaviridae", "Pseudoviridae", "Belpaoviridae", "Retroviridae"):
            fm.setdefault(f, "RNA")
        for f in ("Caulimoviridae", "Hepadnaviridae"):
            fm.setdefault(f, "DNA")
        return fm, gm, cid_na

    _fam_na, _gen_na, _na_by_contig = _s11_na_lookup()

    # ── s03 Flye co-assembly 报告图 (03b_MergeSamples/flye_coassembly/flye_report.png) ──
    try:
        flye_png = report_dir.parent / _D['d_merge'] / "flye_coassembly" / "flye_report.png"
        if flye_png.is_file():
            _fb = _img_to_base64(flye_png, max_kb=10000)
            if _fb:
                stage_extra_html['s03'] = ('<div class="chart-box" id="merge-flye"><div class="chart-title">Flye Co-Assembly</div>'
                                           '<img src="data:image/png;base64,' + _fb + '" loading="lazy" '
                                           'onclick="openLightbox(this.src)" title="Click to view" '
                                           'style="width:100%;height:auto;max-height:560px;object-fit:contain;display:block;cursor:zoom-in"></div>')
    except Exception:
        pass

    # ── S01/S02/S03 明细目录文件清单 ──
    # assembly_detail/ident_detail/cobra_detail 由 utils/*_stats.py 落盘到报告目录,
    # 此前复制后无任何链接入口 (文件在目录里但报告里找不到) → 章节末尾加折叠块,
    # 列目录内文件的相对链接 (新窗口打开); 目录不存在则整块不渲染。
    def _detail_files_block(dirname):
        _dd = report_dir / dirname
        if not _dd.is_dir(): return ""
        _dfiles = sorted(p for p in _dd.rglob("*") if p.is_file())
        if not _dfiles: return ""
        _lis = ""
        for _p in _dfiles:
            _rel = _p.relative_to(report_dir).as_posix()
            _lis += ('<li style="margin:2px 0"><a href="' + _esc(_rel) + '" target="_blank" rel="noopener" '
                     'style="color:var(--blue);font-size:.78rem">' + _esc(_rel) + '</a></li>')
        return ('<details class="stage-table-detail"><summary>Detail files ('
                + str(len(_dfiles)) + ' · ' + _esc(dirname) + ')</summary>'
                '<ul style="margin:6px 0 0 18px;padding:0;list-style:disc">' + _lis + '</ul></details>')

    for _det_sk, _det_dir in (('s01', 'assembly_detail'), ('s02', 'ident_detail'), ('s03', 'cobra_detail')):
        _det_blk = _detail_files_block(_det_dir)
        if _det_blk:
            stage_extra_html[_det_sk] = stage_extra_html.get(_det_sk, '') + _det_blk

    # ── S02 识别阶段图 (02_Identification) ──
    # summary_plots/All_Samples_Comparison.png 汇总图直接内嵌;
    # 逐样本 comparison_plots/{sample}_{venn_diagram,upset_plot,bar_charts}.png
    # 按样本分 details 折叠组 (默认折叠, 每样本最多 3 张全收); 样本数 >12 时
    # 仅内嵌前 12 个, 其余只列文件名清单。图片大小上限沿用 _img_to_base64。
    try:
        _ident_root = report_dir.parent / _D['d_ident']
        _s02_extra = ""
        _sum_png = _ident_root / "summary_plots" / "All_Samples_Comparison.png"
        if _sum_png.is_file():
            _b64s = _img_to_base64(_sum_png)
            if _b64s:
                _s02_extra += ('<div class="stage-charts" style="grid-template-columns:1fr"><div class="chart-box">'
                               '<div class="chart-title">All Samples Comparison</div>'
                               '<img src="data:image/png;base64,' + _b64s + '" loading="lazy" '
                               'onclick="openLightbox(this.src)" title="Click to zoom" '
                               'style="width:100%;height:auto;max-height:620px;object-fit:contain;display:block;cursor:zoom-in">'
                               '</div></div>')
        _s02_samples = []
        if _ident_root.is_dir():
            for _sd in sorted(_ident_root.iterdir()):
                if not _sd.is_dir(): continue
                _cp = _sd / "comparison_plots"
                if not _cp.is_dir(): continue
                _imgs = [(k, _cp / f"{_sd.name}_{k}.png")
                         for k in ("venn_diagram", "upset_plot", "bar_charts")]
                _imgs = [(k, p) for k, p in _imgs if p.is_file()]
                if _imgs: _s02_samples.append((_sd.name, _imgs))
        _S02_EMBED_MAX = 12
        _s02_blocks = ""
        for _sn, _imgs in _s02_samples[:_S02_EMBED_MAX]:
            _inner = ""
            for _kind, _p in _imgs[:3]:
                _b64i = _img_to_base64(_p)
                if not _b64i: continue
                _inner += ('<div class="chart-box"><div class="chart-title">' + _esc(_sn + " · " + _kind) + '</div>'
                           '<img src="data:image/png;base64,' + _b64i + '" loading="lazy" '
                           'onclick="openLightbox(this.src)" title="Click to zoom" '
                           'style="width:100%;height:auto;max-height:480px;object-fit:contain;display:block;cursor:zoom-in"></div>')
            if _inner:
                _ncol = "1fr 1fr 1fr" if len(_imgs) >= 3 else ("1fr 1fr" if len(_imgs) == 2 else "1fr")
                _s02_blocks += ('<details class="stage-table-detail"><summary>' + _esc(_sn) + ' — '
                                + str(len(_imgs)) + ' plots</summary>'
                                '<div class="stage-charts" style="grid-template-columns:' + _ncol
                                + ';padding:8px 0">' + _inner + '</div></details>')
        if _s02_blocks:
            _s02_extra += '<div class="stage-tables">' + _s02_blocks + '</div>'
        if len(_s02_samples) > _S02_EMBED_MAX:
            _rest = ""
            for _sn, _imgs in _s02_samples[_S02_EMBED_MAX:]:
                _rest += ('<li style="margin:2px 0">' + _esc(_sn) + ': '
                          + _esc(", ".join(k + ".png" for k, _p in _imgs)) + '</li>')
            _s02_extra += ('<div class="stage-tables"><details class="stage-table-detail"><summary>'
                           + str(len(_s02_samples) - _S02_EMBED_MAX) + ' more samples (file names only)</summary>'
                           '<ul style="margin:6px 0 0 18px;padding:0;list-style:disc">' + _rest + '</ul></details></div>')
        if _s02_extra:
            stage_extra_html['s02'] = stage_extra_html.get('s02', '') + _s02_extra
    except Exception:
        pass

    # ── S11b KEEP 序列快照表 (点击 contig_id 复制序列) ──
    keep_table_html = ""
    try:
        import csv as _kcsv
        import json as _kjson
        ksum = _resolve_acv_dir(report_dir.parent) / "keep_summary.tsv"
        keep_fa = _resolve_acv_dir(report_dir.parent) / "class_KEEP.fasta"
        if ksum.is_file():
            seqs = {}
            if keep_fa.is_file():
                _kid = None; _kbuf = []
                for _line in open(keep_fa, encoding='utf-8', errors='ignore'):
                    if _line.startswith('>'):
                        if _kid: seqs[_kid] = ''.join(_kbuf)
                        _kid = _line[1:].strip().split()[0]; _kbuf = []
                    else:
                        _kbuf.append(_line.strip())
                if _kid: seqs[_kid] = ''.join(_kbuf)
            ksrows = list(_kcsv.DictReader(open(ksum, encoding='utf-8'), delimiter='\t'))
            # CT3-HMM 第三路探针注释: 从 scored.tsv 按 contig_id 左连 ct3_* / hmm_rescue 列
            _ct3_by_cid = {}
            _scored_p = report_dir / "rescue_evidence_scored.tsv"
            if _scored_p.is_file():
                try:
                    for _sr in _kcsv.DictReader(open(_scored_p, encoding='utf-8'), delimiter='\t'):
                        _cid0 = (_sr.get('contig_id','') or '').strip()
                        if _cid0: _ct3_by_cid[_cid0] = _sr
                except Exception:
                    _ct3_by_cid = {}
            cols = ['contig_id','length','length_status','genus_avg_len','len_ratio',
                    'Realm','Kingdom','Phylum','Class','Order','Family','Genus','Species','Nucleic_acid','cf_max_conf',
                    'cdd_evidence','n_t1','cdd_top','max_cov_to_db','n_samples','prevalence_%',
                    'checkv_completeness','checkv_confidence',
                    'nt_pident','nt_qcovs','nt_species','aa_pident','aa_species','blastn_score','blastx_score',
                    'hmm_score','domain_score','total_score',
                    'final_verdict','final_call','conflict','novelty_score','prediction','novel','category_9a',
                    'tax_source','identify_tools','classify_tool','host_assign',
                    'ct3_n_viral','ct3_n_family','ct3_best_viral_tier','ct3_best_viral_model',
                    'ct3_best_viral_evalue','hmm_rescue']
            trs = []; seqdata = {}
            for r in ksrows:
                cid = (r.get('contig_id','') or '').strip()
                if cid in _ct3_by_cid:
                    r = {**_ct3_by_cid[cid], **r}
                # 核酸类型: contig 级左连 All_plant/HQ 表 (_na_by_contig),
                # 缺失时回退科/属多数派映射 (_fam_na/_gen_na, 含反转录科兜底)
                _na = (_na_by_contig.get(cid)
                       or _fam_na.get((r.get('Family','') or '').strip())
                       or _gen_na.get((r.get('Genus','') or '').strip()) or '')
                lr = (r.get('len_ratio','') or '').strip()
                try: lrf = float(lr)
                except: lrf = None
                if lrf is None: st = 'KEEP_no_genus_len'
                elif lrf < 0.7: st = 'FILTERED'
                else: st = 'KEEP'
                tds = []
                for c in cols:
                    if c == 'contig_id':
                        tds.append('<td><a class="kseq" href="javascript:void(0)" data-cid="'+_esc(cid)+'">'+_esc(cid)+'</a></td>')
                    elif c == 'Nucleic_acid':
                        _nc = {"DNA":"#1565c0","RNA":"#2e7d32","Viroid":"#6a1b9a"}.get(_na)
                        if _nc:
                            tds.append('<td><span style="display:inline-block;padding:1px 9px;border-radius:9px;'
                                       'font-size:.72rem;font-weight:600;color:#fff;background:'+_nc+'">'+_esc(_na)+'</span></td>')
                        else:
                            tds.append('<td style="color:var(--muted,#99a)">-</td>')
                    else:
                        tds.append('<td>'+_esc(r.get(c,''))+'</td>')
                tds.append('<td><span class="kflag '+st+'">'+st+'</span></td>')
                trs.append('<tr>'+''.join(tds)+'</tr>')
                if cid in seqs: seqdata[cid] = seqs[cid]
            ths = ''.join('<th>'+_esc(c)+'</th>' for c in cols)+'<th>filter_status</th>'
            sjson = _kjson.dumps(seqdata, ensure_ascii=False)
            keep_table_html = (
                '<div class="stage-tables"><div class="chart-title">KEEP Sequences — click a contig to copy</div>'
                '<div style="overflow-x:auto"><table class="app-table" id="keep-table" style="width:100%;white-space:nowrap;font-size:.8rem;border-collapse:collapse">'
                '<thead><tr>'+ths+'</tr></thead><tbody>'+''.join(trs)+'</tbody></table></div></div>'
                '<script>window.__keepSeqs='+sjson+';'
                '(function(){function cp(id){var s=window.__keepSeqs[id];if(!s)return;'
                'navigator.clipboard.writeText(">"+id+"\\n"+s).then(function(){var el=document.querySelector("[data-cid=\\""+id+"\\"]");if(el){var o=el.textContent;el.textContent=o+" ✓copied";setTimeout(function(){el.textContent=o;},900);}}).catch(function(){alert("copy failed");});}'
                'Array.prototype.forEach.call(document.querySelectorAll(".kseq"),function(el){el.onclick=function(){cp(el.getAttribute("data-cid"));};});})();</script>')
    except Exception:
        keep_table_html = ""
    # ── S05 Taxonomy 额外图区 (组件先累积, 末尾一次性挂载) ──
    # 组件: taxonomy_images.json 静态 PNG + 交互 UpSet 容器 + 07_Taxonomy 四类图
    _s05_extra = ""
    tax_img_json = report_dir / "taxonomy_images.json"
    if tax_img_json.is_file():
        try:
            import json as _json3
            with open(tax_img_json, encoding='utf-8') as tif:
                tax_images = _json3.load(tif)
            img_html = ""
            for label, b64 in tax_images.items():
                img_html += '<div class="chart-box"><div class="chart-title">' + _esc(label) + '</div>'
                img_html += '<img src="data:image/png;base64,' + b64 + '" loading="lazy" onclick="openLightbox(this.src)" title="Click to zoom" style="width:100%;height:auto;max-height:650px;object-fit:contain;display:block;cursor:zoom-in"></div>'
            if img_html:
                _s05_extra += '<div class="stage-charts" style="grid-template-columns:1fr">' + img_html + '</div>'
        except Exception: pass
    if upset_container_html:
        _s05_extra += '<div class="stage-charts" style="grid-template-columns:1fr">' + upset_container_html + '</div>'

    # ── S05 补充: 07_Taxonomy/{sample}.integrated 四类图 ──
    # intersection_venn / consensus_summary / classification_counts / agreement_rates,
    # 逐样本 {Votus}.integrated 目录; 优先 PNG base64 内嵌, 缺 PNG 回退 PDF 用
    # <embed type=application/pdf> (同 SDT 通道既有写法, 大小上限沿用其 15000KB 阈值)。
    try:
        _tax_root = report_dir.parent / _D['d_taxonomy']
        _tax_boxes = []
        if _tax_root.is_dir():
            for _idir in sorted(_tax_root.glob("*.integrated")):
                if not _idir.is_dir(): continue
                for _base, _lbl in (("intersection_venn", "Intersection Venn"),
                                    ("consensus_summary", "Consensus Summary"),
                                    ("classification_counts", "Classification Counts"),
                                    ("agreement_rates", "Agreement Rates")):
                    _pngf = _idir / f"{_base}.png"
                    _pdff = _idir / f"{_base}.pdf"
                    _imgf = _pngf if _pngf.is_file() else (_pdff if _pdff.is_file() else None)
                    if not _imgf: continue
                    _b64t = _img_to_base64(_imgf, max_kb=15000)
                    if not _b64t: continue
                    if _imgf.suffix.lower() == '.pdf':
                        _body = ('<embed src="data:application/pdf;base64,' + _b64t + '" type="application/pdf" '
                                 'style="width:100%;height:600px;display:block">')
                    else:
                        _body = ('<img src="data:image/png;base64,' + _b64t + '" loading="lazy" '
                                 'onclick="openLightbox(this.src)" title="Click to zoom" '
                                 'style="width:100%;height:auto;max-height:620px;object-fit:contain;display:block;cursor:zoom-in">')
                    _tax_boxes.append('<div class="chart-box"><div class="chart-title">'
                                      + _esc(_idir.name.replace(".integrated", "")) + ' · ' + _lbl
                                      + '</div>' + _body + '</div>')
        if _tax_boxes:
            _s05_extra += ('<div class="stage-charts" style="grid-template-columns:1fr 1fr">'
                           + ''.join(_tax_boxes) + '</div>')
    except Exception:
        pass
    # subnav #upset-container 指向 UpSet 图区, 但交互容器仅在 upset_data.json 存在时
    # 渲染 → 无 upset 时成死锚。恒存在的稳定别名锚: 交互容器未渲染时补一个空 span
    # 承载该 id (与实际容器互斥, 不产生重复 id), 交互/静态回退两条路径都命中。
    if upset_container_html:
        stage_extra_html['s05'] = _s05_extra
    else:
        stage_extra_html['s05'] = '<span id="upset-container"></span>' + _s05_extra

    _s11_tree = ""; _s11_sdt = ""; _s11_cv = ""

    # ── S11 分组展示 helper (核酸类型查表 _s11_na_lookup 见本函数前部) ──
    def _s11_bucket(name, tmap):
        if (name or '') in _VIROID_FAMS or "viroid" in (name or '').lower():
            return "Viroid"
        return tmap.get(name) or "Other"

    _S11_GROUPS = (("DNA", "DNA virus"), ("RNA", "RNA virus"), ("Viroid", "Viroid"), ("Other", "Unassigned type"))
    _S11_GROUP_COLOR = {"DNA": "#1565c0", "RNA": "#2e7d32", "Viroid": "#6a1b9a", "Other": "#78909c"}
    _S11_FAM_LABEL = {"DNA": "DNA virus families", "RNA": "RNA virus families",
                      "Viroid": "Viroid families", "Other": "Unassigned-type families"}
    _S11_GEN_LABEL = {"DNA": "DNA virus genera", "RNA": "RNA virus genera",
                      "Viroid": "Viroid genera", "Other": "Unassigned-type genera"}

    def _s11_grouped(entries, tmap):
        """entries: [(名称, box_html)] → {组: [box_html]}, 组内按名称排序"""
        buckets = {}
        for name, box in entries:
            buckets.setdefault(_s11_bucket(name, tmap), []).append((name, box))
        return {g: [b for _n, b in sorted(buckets.get(g, []), key=lambda x: x[0].lower())]
                for g, _ in _S11_GROUPS}

    def _s11_group_html(prefix, glbl, gkey, boxes, cols):
        color = _S11_GROUP_COLOR[gkey]
        head = (f'<div style="margin:2px 0 0;padding:7px 12px;border:1px solid {color};'
                'background:rgba(0,0,0,0.035);border-radius:4px;font-size:.92rem;">'
                f'{prefix} · <b style="color:{color}">{glbl}</b> '
                f'<span style="font-weight:400;color:var(--muted,#667)">(n={len(boxes)})</span></div>')
        grid = ('<div class="stage-charts" style="grid-template-columns:'
                + ' '.join(['1fr']*max(cols, 1)) + '">' + ''.join(boxes) + '</div>')
        return head + grid

    # ── S11 SDT 属级序列一致性矩阵图区 (DNA 属 / RNA 属分组) ──
    try:
        sdt_root = (_resolve_acv_dir(report_dir.parent) / "SDT_matrix")
        if sdt_root.is_dir():
            entries = []
            for d in sorted(sdt_root.iterdir()):
                if not d.is_dir(): continue
                png = d / f"SDT_{d.name}_heatmap.png"
                pdf = d / f"SDT_{d.name}_heatmap.pdf"
                imgf = png if png.exists() else (pdf if pdf.exists() else None)
                if not imgf: continue
                b64 = _img_to_base64(imgf, max_kb=15000)
                if not b64: continue
                csvf = d / f"sdt_matrix_{d.name}.csv"
                csv_uri = _tsv_data_uri(csvf) if csvf.is_file() else ""
                dl = (f' <a href="{csv_uri}" download style="font-size:.75rem;color:#00646e" title="Download {d.name} matrix CSV">matrix.csv</a>'
                      if csv_uri else '')
                # PDF 回退不能塞进 data:image/png (浏览器不识别) → 用 <embed type=application/pdf>
                if imgf.suffix.lower() == '.pdf':
                    body = ('<embed src="data:application/pdf;base64,' + b64 + '" type="application/pdf" '
                            'style="width:100%;height:600px;display:block">')
                else:
                    body = ('<img src="data:image/png;base64,' + b64 + '" loading="lazy" '
                            'onclick="openLightbox(this.src)" title="Click to zoom" '
                            'style="width:100%;height:auto;max-height:600px;object-fit:contain;display:block;cursor:zoom-in">')
                entries.append((d.name,
                    '<div class="chart-box"><div class="chart-title">SDT · ' + _esc(d.name) + dl + '</div>'
                    + body + '</div>'))
            if entries:
                _sdt_grp = _s11_grouped(entries, _gen_na)
                _sdt_parts = []
                for gkey, _glbl in _S11_GROUPS:
                    boxes = _sdt_grp.get(gkey) or []
                    if not boxes: continue
                    _sdt_parts.append(_s11_group_html("SDT matrix", _S11_GEN_LABEL[gkey], gkey, boxes,
                                                      cols=(2 if len(boxes) >= 2 else 1)))
                _s11_sdt = '<div id="s11-sdt">' + ''.join(_sdt_parts) + '</div>'
    except Exception:
        pass

    # ── S11 KEEP reads 覆盖曲线图区 ──
    try:
        cv_root = (_resolve_acv_dir(report_dir.parent) / "coverage_curves")
        if cv_root.is_dir():
            cv_imgs = []
            for p in sorted(cv_root.glob("*.png")):
                b64 = _img_to_base64(p, max_kb=10000)
                if not b64: continue
                lbl = p.stem.replace("__", " | ")
                cv_imgs.append(
                    '<div class="chart-box"><div class="chart-title">Coverage · ' + _esc(lbl) + '</div>'
                    '<img src="data:image/png;base64,' + b64 + '" loading="lazy" '
                    'onclick="openLightbox(this.src)" title="Click to zoom" '
                    'style="width:100%;height:auto;max-height:420px;object-fit:contain;display:block;cursor:zoom-in"></div>')
            if cv_imgs:
                ncol2 = 2 if len(cv_imgs) >= 2 else 1
                _s11_cv = ('<div class="stage-charts" id="s11-cv" style="grid-template-columns:'
                           + ' '.join(['1fr']*min(ncol2,2)) + '">' + ''.join(cv_imgs) + '</div>')
    except Exception:
        pass

    # ── S11 进化树图区 (acvirus_trees/<Family>/Tree_Synteny_Composite.png) ──
    # DNA 科与 RNA 科分开展示 (组内按科名字母序)
    try:
        tr_root = (_resolve_acv_dir(report_dir.parent) / "acvirus_trees")
        if tr_root.is_dir():
            entries = []
            for d in sorted(tr_root.iterdir()):
                if not d.is_dir(): continue
                imgf = d / "Tree_Synteny_Composite.png"
                if not imgf.is_file(): continue
                b64 = _img_to_base64(imgf, max_kb=22000)
                if not b64: continue
                entries.append((d.name,
                    '<div class="chart-box"><div class="chart-title">Phylogeny · ' + _esc(d.name) + '</div>'
                    '<img src="data:image/png;base64,' + b64 + '" loading="lazy" '
                    'onclick="openLightbox(this.src)" title="Click to zoom" '
                    'style="width:100%;height:auto;max-height:560px;object-fit:contain;display:block;cursor:zoom-in"></div>'))
            if entries:
                _tr_grp = _s11_grouped(entries, _fam_na)
                _tr_parts = []
                for gkey, _glbl in _S11_GROUPS:
                    boxes = _tr_grp.get(gkey) or []
                    if not boxes: continue
                    _tr_parts.append(_s11_group_html("Phylogeny", _S11_FAM_LABEL[gkey], gkey, boxes,
                                                     cols=(2 if len(boxes) >= 2 else 1)))
                _s11_tree = '<div id="s11-tree">' + ''.join(_tr_parts) + '</div>'
    except Exception:
        pass
    # ── EVE 判别 (KEEP 前的否决层) ──
    _s11_eve = ''
    try:
        from collections import Counter as _Counter
        _acv = _resolve_acv_dir(Path(report_dir).parent)
        _scored = _acv / 'virus_validation' / 'rescue_evidence_scored.tsv'
        if _scored.is_file():
            _cnt = _Counter()
            _veto = 0
            with open(_scored, encoding='utf-8', errors='replace') as _f:
                for _r in csv.DictReader(_f, delimiter='	'):
                    _c = (_r.get('eve_call') or '').strip()
                    if _c:
                        _cnt[_c] += 1
                    if (_r.get('eve_flag') or '').strip() == 'KEEP->REVIEW':
                        _veto += 1
            if _cnt:
                _n = sum(_cnt.values())
                _cards = (f'<div class="kpi-card"><div class="kpi-value">{_n}</div>'
                          f'<div class="kpi-label">Candidates scored</div></div>'
                          f'<div class="kpi-card"><div class="kpi-value">{_cnt.get("EVE_like", 0)}</div>'
                          f'<div class="kpi-label">EVE-like</div></div>'
                          f'<div class="kpi-card"><div class="kpi-value">{_cnt.get("virus_like", 0)}</div>'
                          f'<div class="kpi-label">Virus-like</div></div>'
                          f'<div class="kpi-card"><div class="kpi-value">{_veto}</div>'
                          f'<div class="kpi-label">KEEP to REVIEW (EVE veto)</div></div>')
                _link = ''
                _eh = _acv / 'virus_validation' / 'eve' / 'eve_report.html'
                if _eh.is_file():
                    _dst = Path(report_dir) / 'EVE_discrimination.html'
                    try:
                        if (not _dst.is_file()
                                or _dst.stat().st_mtime < _eh.stat().st_mtime):
                            _dst.write_bytes(_eh.read_bytes())
                        _link = ('<p style="font-size:12px;margin-top:8px">'
                                 '<a href="EVE_discrimination.html" target="_blank">'
                                 'Open EVE discrimination report</a></p>')
                    except OSError:
                        pass
                _s11_eve = ('<div style="margin-top:18px">'
                            '<h3 style="margin:0 0 8px">EVE discrimination (before KEEP)</h3>'
                            '<div class="kpi-row">' + _cards + '</div>' + _link + '</div>')
    except Exception:
        pass
    stage_extra_html['s11'] = _s11_eve + (_s11_tree or '') + (_s11_sdt or '') + (_s11_cv or '')
    for sk, short, full, icon, color in stage_defs:
        st = stage_status.get(sk, 'skip'); metric = stage_metric.get(sk, '')
        if st == 'pass': badge_cls, badge_txt, border_cls = 's-pass','✓ PASS','stage-pass'
        elif st == 'fail': badge_cls, badge_txt, border_cls = 's-fail','✗ FAIL','stage-fail'
        else: badge_cls, badge_txt, border_cls = 's-skip','○ SKIP','stage-skip'

        chart_html = ""
        if sk in chart_map:
            chs = chart_map[sk]; active = []
            for cid, _ in chs:
                if sk == 's00a':
                    if cid == 'chart_s00a' and stage_has_chart.get('s00a'): active.append(cid)
                    if cid == 'chart_s00a_dup' and dq_rows and any(float(r.get('Dup_Rate(%)',0))>0 for r in dq_rows if r.get('Sample','')!='TOTAL'): active.append(cid)
                elif sk == 's00b':
                    if stage_has_chart.get('s00b'): active.append(cid)
                elif sk == 's01':
                    if cid == 'chart_s01a' and stage_has_chart.get('s01a'): active.append(cid)
                    if cid == 'chart_s01b' and stage_has_chart.get('s01b'): active.append(cid)
                elif sk == 's02':
                    if stage_has_chart.get('s02'): active.append(cid)
                elif sk == 's04':
                    if cid == 'chart_s04a' and stage_has_chart.get('s04a'): active.append(cid)
                    if cid == 'chart_s04c' and stage_has_chart.get('s04c'): active.append(cid)
                elif sk == 's05':
                    if cid == 'chart_s05a' and stage_has_chart.get('s05'): active.append(cid)
                    if cid == 'chart_s05b' and stage_has_chart.get('s05b'): active.append(cid)
                    if cid == 'chart_s05c' and stage_has_chart.get('s05c'): active.append(cid)
                    if cid == 'chart_s05d' and stage_has_chart.get('s05d'): active.append(cid)
                    if cid == 'chart_s05e' and stage_has_chart.get('s05e'): active.append(cid)
                    if cid == 'chart_s05f' and stage_has_chart.get('s05f'): active.append(cid)
                elif sk == 's07':
                    if cid == 'chart_s07a' and stage_has_chart.get('s07a'): active.append(cid)
                    if cid == 'chart_s07b' and stage_has_chart.get('s07b'): active.append(cid)
                elif sk == 's09':
                    if cid == 'chart_s09' and stage_has_chart.get('s09'): active.append(cid)
                    if cid == 'chart_s09b' and stage_has_chart.get('s09b'): active.append(cid)
                    if cid == 'chart_s09c' and stage_has_chart.get('s09c'): active.append(cid)
                    if cid == 'chart_s09d' and stage_has_chart.get('s09d'): active.append(cid)
                    if cid == 'chart_s09e' and stage_has_chart.get('s09e'): active.append(cid)
                    if cid == 'chart_s09f' and stage_has_chart.get('s09f'): active.append(cid)
                    if cid == 'chart_s09g' and stage_has_chart.get('s09g'): active.append(cid)
                    if cid == 'chart_s09h' and stage_has_chart.get('s09h'): active.append(cid)
                elif sk == 's10':
                    if cid == 'chart_s10na' and stage_has_chart.get('s10na'): active.append(cid)
                    if cid == 'chart_s10a' and stage_has_chart.get('s10'): active.append(cid)
                    if cid == 'chart_s10b' and stage_has_chart.get('s10b'): active.append(cid)
                    if cid == 'chart_s09' and stage_has_chart.get('s09'): active.append(cid)
                    if cid == 'chart_s09b' and stage_has_chart.get('s09b'): active.append(cid)
                    if cid == 'chart_s09c' and stage_has_chart.get('s09c'): active.append(cid)
                    if cid == 'chart_s09d' and stage_has_chart.get('s09d'): active.append(cid)
                    if cid == 'chart_s09e' and stage_has_chart.get('s09e'): active.append(cid)
                    if cid == 'chart_s09f' and stage_has_chart.get('s09f'): active.append(cid)
                    if cid == 'chart_s09g' and stage_has_chart.get('s09g'): active.append(cid)
                    if cid == 'chart_s09h' and stage_has_chart.get('s09h'): active.append(cid)
                    if cid == 'chart_s09i' and stage_has_chart.get('s09i'): active.append(cid)
                elif sk == 's11':
                    if cid == 'chart_s11a' and stage_has_chart.get('s11a'): active.append(cid)
                else:
                    if stage_has_chart.get(sk): active.append(cid)
            if active:
                cols = '1fr' if len(active) == 1 else '1fr 1fr'
                chart_html = f'<div class="stage-charts" style="grid-template-columns:{cols}">'
                for cid in active:
                    is_stacked = cid in ('chart_s00b', 'chart_s07a', 'chart_s07b')
                    is_per_sample = cid in ('chart_s00a','chart_s00a_dup','chart_s01a','chart_s01b')
                    chart_html += f'<div class="chart-box">'
                    if is_stacked:
                        chart_html += f'<div style="display:flex;justify-content:flex-end;margin-bottom:6px"><button class="pct-btn" data-chart="{cid}" data-mode="abs" onclick="toggleStackedPct(this)">Show %</button></div>'
                    if is_per_sample:
                        chart_html += '<div class="chart-scroll">'
                    chart_html += f'<canvas id="{cid}" style="max-height:320px"></canvas>'
                    if is_per_sample:
                        chart_html += '</div>'
                    chart_html += '</div>'
                chart_html += '</div>'

        if sk in sankey_by_stage:
            chart_html += f'<div class="sankey-section">{sankey_by_stage[sk]}</div>'

        # 在卡片内嵌入对应 TSV 数据表
        table_html = ""
        _anchor_ids = set()                      # 本 section 表格锚点 id (真实 id + 别名锚)
        for tsv_name in stage_tsv_map.get(sk, []):
            tsv_path = report_dir / tsv_name
            if not tsv_path.is_file(): continue
            tsv_rows = _read_tsv(tsv_path)
            if not tsv_rows: continue
            is_large = len(tsv_rows) > 500
            use_pagination = len(tsv_rows) > 10 and not is_large
            preview = tsv_rows[:10]
            cols = list(preview[0].keys()) if preview else []
            th_h = "".join(f"<th>{_esc(c)}</th>" for c in cols)
            tr_h = ""
            if is_large:
                import json as _json3
                _load_rows = tsv_rows[:500]
                _js_data = _json3.dumps([{c: str(r.get(c,''))[:100] for c in cols} for r in _load_rows], ensure_ascii=False)
                tbid = _tbl_id(sk, tsv_name) + "_large"
                tr_h = '<tr><td colspan="' + str(len(cols)) + '" style="text-align:center;color:var(--muted)">Loading...</td></tr>'
            elif use_pagination:
                for r in tsv_rows:
                    tr_h += "<tr class='pag-row' style='display:none'>" + "".join(f"<td>{_esc(str(r.get(c,'')))[:60]}</td>" for c in cols) + "</tr>"
            else:
                for r in preview:
                    tr_h += "<tr>" + "".join(f"<td>{_esc(str(r.get(c,'')))[:60]}</td>" for c in cols) + "</tr>"
            more = f'<span style="color:var(--muted);font-size:11px">({len(tsv_rows)} rows{" - showing first 500" if is_large else ""})</span>'
            if is_large:
                dl_uri = ""
            else:
                dl_uri = _tsv_data_uri(tsv_path)
            dl_link = f'<a href="{dl_uri}" download="{tsv_name}" style="font-size:11px;margin-left:6px;color:var(--blue)">[download]</a>' if dl_uri else ''
            _short_tbid = _tbl_id(sk, tsv_name)      # 与 _subnav 共用同一 id 规则
            if not is_large:
                tbid = _short_tbid
            pager_html = ""
            if is_large:
                pager_html = f'<div id="pager_{tbid}" style="display:flex;justify-content:center;align-items:center;gap:8px;padding:8px;font-size:11px;color:var(--muted)"><button onclick="tblPrev_{tbid}()" style="border:1px solid #ccc;background:#fff;padding:3px 10px;border-radius:4px;cursor:pointer;font-size:11px">Prev</button><span id="info_{tbid}">Page 1</span><button onclick="tblNext_{tbid}()" style="border:1px solid #ccc;background:#fff;padding:3px 10px;border-radius:4px;cursor:pointer;font-size:11px">Next</button></div>'
                table_pager_js += f"""
var _tblRows_{tbid} = {_js_data};
var _tblCols_{tbid} = {_json3.dumps(cols, ensure_ascii=False)};
var _tblPage_{tbid} = 0;
var _tblPerPage_{tbid} = 10;
function tblShow_{tbid}(p){{
  _tblPage_{tbid} = p;
  var pp=_tblPerPage_{tbid}, start=p*pp;
  var page=_tblRows_{tbid}.slice(start,start+pp);
  var h='';
  page.forEach(function(r){{ h+='<tr>'+_tblCols_{tbid}.map(function(c){{ var v=r[c]||''; return '<td>'+v.substring(0,60)+'</td>'; }}).join('')+'</tr>'; }});
  document.getElementById('{tbid}').querySelector('tbody').innerHTML=h;
  var total=Math.ceil(_tblRows_{tbid}.length/pp);
  document.getElementById('info_{tbid}').textContent='Page '+(p+1)+' of '+total;
}}
function tblPrev_{tbid}(){{ if(_tblPage_{tbid}>0) tblShow_{tbid}(_tblPage_{tbid}-1); }}
function tblNext_{tbid}(){{ var total=Math.ceil(_tblRows_{tbid}.length/_tblPerPage_{tbid}); if(_tblPage_{tbid}<total-1) tblShow_{tbid}(_tblPage_{tbid}+1); }}
tblShow_{tbid}(0);
"""
            elif use_pagination:
                pager_html = f'<div id="pager_{tbid}" style="display:flex;justify-content:center;gap:8px;padding:8px;font-size:11px;color:var(--muted)"></div>'
                table_pager_js += f"""
(function(){{
  var tbody=document.getElementById('{tbid}').querySelector('tbody');
  var rows=Array.from(tbody.querySelectorAll('.pag-row'));
  var perPage=10,totalPages=Math.ceil(rows.length/perPage);
  var pager=document.getElementById('pager_{tbid}');
  var page=0;
  function show(p){{
    page=Math.max(0,Math.min(p,totalPages-1));
    rows.forEach(function(r,i){{r.style.display=(i>=page*perPage&&i<(page+1)*perPage)?'':'none'}});
    pager.innerHTML='<button onclick="void(0)" '+(page===0?'disabled':'')+' style="border:1px solid #ccc;background:#fff;padding:3px 10px;border-radius:4px;cursor:pointer;font-size:11px">Prev</button>'+'<span>Page <b>'+(page+1)+'</b> of '+totalPages+'</span>'+'<button '+(page===totalPages-1?'disabled':'')+' style="border:1px solid #ccc;background:#fff;padding:3px 10px;border-radius:4px;cursor:pointer;font-size:11px">Next</button>';
    pager.querySelectorAll('button')[0].onclick=function(){{show(page-1)}};
    pager.querySelectorAll('button')[1].onclick=function(){{show(page+1)}};
  }}
  show(0);
}})();
"""
            # 大表 (>500 行) id 带 _large 后缀 → _subnav 写死的短名锚点 (_short_tbid)
            # 会变死锚; 大表恒补一个稳定别名锚承载短名 (小表 id 即短名, 不重复输出)
            alias_anchor = f'<span id="{_short_tbid}"></span>' if is_large else ""
            _anchor_ids.add(tbid)
            if alias_anchor: _anchor_ids.add(_short_tbid)
            table_html += f'''{alias_anchor}<details class="stage-table-detail" open>
<summary>{tsv_name} {more} {dl_link}</summary>
<div style="overflow-x:auto;margin-top:6px">
<table class="app-table" id="{tbid}"><thead><tr>{th_h}</tr></thead><tbody>{tr_h}</tbody></table>
{pager_html}
</div>
</details>'''
        if table_html:
            if sk == 's04':
                table_html = f'<div class="stage-charts" style="grid-template-columns:1fr 1fr">{table_html}</div>'
            else:
                table_html = f'<div class="stage-tables">{table_html}</div>'

        extra_html = stage_extra_html.get(sk, "")
        if sk == 's11' and keep_table_html:
            extra_html = keep_table_html + extra_html
        # subnav 二级锚兜底: 对应表/图未渲染 (阶段未跑/无产物) 时同样补一个
        # 稳定别名锚, 保证 sidebar 二级锚不失效 (聚合本 section 实际输出的全部 id)
        _section_ids = set(_anchor_ids)
        # 只认属性位上的 id= (前导空白), 避免误匹配 data-id= / data-cid= 之类
        _section_ids |= set(re.findall(r'\sid="([^"]+)"', table_html + chart_html + extra_html))
        _subnav_anchor_html = ""
        for _shref, _slbl in _subnav.get(sk, []):
            _aid = _shref.lstrip('#')
            if _aid and _aid not in _section_ids:
                _subnav_anchor_html += f'<span id="{_aid}"></span>'
        sections_html += f'''
<section class="stage {border_cls}" id="stage-{short}">
  {_subnav_anchor_html}
  <div class="stage-header">
    <div class="stage-icon" style="background:{color}">{icon}</div>
    <div class="stage-title"><h2>{full} <button onclick="runStageAI('{sk}')" class="ai-stage-btn" title="AI summarize">AI</button></h2><span class="stage-metric">{_esc(metric)}</span></div>
    <span class="stage-badge {badge_cls}">{badge_txt}</span>
  </div>
  {table_html}
  {chart_html}
  {extra_html}
</section>'''

    # ── Table ──
    table_rows = ""
    for i, s in enumerate(main_stages):
        cls = S.get(s["Status"], "skip")
        badge_label = {"pass":"PASS","skip":"SKIP","fail":"FAIL"}.get(cls,"SKIP")
        table_rows += f'<tr class="tr-{cls}"><td class="td-num">{i}</td><td><b>{_esc(s["Stage"])}</b></td><td><span class="tb-badge tb-{cls}">{badge_label}</span></td><td class="td-metric">{_esc(s.get("Key_Metric",""))}</td></tr>'

    # ── Sidebar nav ──
    sidebar_items = ""
    for sk, short, full, icon, color in stage_defs:
        st = stage_status.get(sk, 'skip')
        item_cls = 'sb-pass' if st=='pass' else ('sb-fail' if st=='fail' else 'sb-skip')
        metric_text = stage_metric.get(sk, '')
        sidebar_items += f'<a href="#stage-{short}" class="sb-item {item_cls}" title="{metric_text}"><span class="sb-dot" style="background:{color}"></span><span class="sb-label">{full}</span></a>'
        for sub_href, sub_label in _subnav.get(sk, []):
            sidebar_items += f'<a href="{sub_href}" class="sb-sub" title="{sub_label}"><span class="sb-dot sb-dot-sub"></span><span class="sb-label">{sub_label}</span></a>' 



    gen_time = datetime.now().strftime('%Y-%m-%d %H:%M:%S')

    # KPI cards
    def _fmt_bases(v):
        """格式化碱基数为人可读: 1234567890 → 1.23 Gb"""
        try: n = int(v.replace(',',''))
        except: return v
        if n >= 1e9: return f"{n/1e9:.1f} Gb"
        if n >= 1e6: return f"{n/1e6:.1f} Mb"
        if n >= 1e3: return f"{n/1e3:.1f} Kb"
        return str(n)
    kpi_cards = ""
    raw_b = kpis.get('raw_bases',''); clean_b = kpis.get('clean_bases','')
    n_sample = kpis.get('n_sample','—')
    sample_kpi = "—"
    if raw_b and clean_b:
        sample_kpi = f"{_fmt_bases(raw_b)} raw → {_fmt_bases(clean_b)} clean"
    # Novelty calculation
    n_known = int(kpis.get('n_known','0').replace(',','') or 0)
    n_newsp = int(kpis.get('n_newsp','0').replace(',','') or 0)
    n_newge = int(kpis.get('n_newge','0').replace(',','') or 0)
    n_newfa = int(kpis.get('n_newfa','0').replace(',','') or 0)
    n_total_tax = n_known + n_newsp + n_newge + n_newfa
    n_novel = n_newsp + n_newge + n_newfa
    novelty_pct = f"{n_novel/n_total_tax*100:.0f}%" if n_total_tax > 0 else "—"
    # CheckV
    hq = kpis.get('hq_votus','—')
    cv_t = kpis.get('cv_total','')
    checkv_kpi = f"{hq} Complete / {cv_t} total" if cv_t else str(hq)
    # Host depletion — 显示数据量 (Gb/Mb), 从 reads 按平均 read 长度换算
    hd_kpi = "—"
    if kpis.get('hd_retained'):
        hd_raw_str = _fmt_bases(str(int(kpis['hd_raw_bp'])))
        hd_after_str = _fmt_bases(str(int(kpis['hd_after_bp'])))
        hd_kpi = f"{hd_raw_str} → {hd_after_str}<br>({kpis['hd_retained']}% retained)"
    kpi_items = [
        ("Samples", f"{sample_kpi}<br>{n_sample} samples", "🧬"),
        ("Host Depletion", hd_kpi, "🧹"),
        ("Assembly", f"{kpis.get('total_contigs','—')} contigs<br>{kpis.get('total_mb','—')} Mb", "🔧"),
        ("Viruses", f"{kpis.get('virus_seqs','—')} identified", "🦠"),
        ("vOTU Clusters", f"{kpis.get('n_clusters','—')}", "📦"),
        ("Novelty", f"{n_novel}/{n_total_tax} novel<br>({novelty_pct})" if n_total_tax > 0 else "—", "🆕"),
        ("Hosts", f"{kpis.get('host_total','—')} classified", "🌐"),
        ("CheckV", checkv_kpi, "✅"),
    ]
    for title, value, icon in kpi_items:
        kpi_cards += f'<div class="kpi-card"><div class="kpi-icon">{icon}</div><div class="kpi-value">{value}</div><div class="kpi-label">{title}</div></div>'

    # ── AI 总结 (IMRaD 格式) ──
    ai_summary_html = ""
    if getattr(generate_ai_summary, '_api_key', ''):
        print("  生成 AI 总结 (IMRaD)...")
        ai_summary_html = generate_ai_summary(stage_stats, kpis, report_dir) or ""

    # ── Pipeline Flow 图: 展示从 raw reads → HQ vOTUs 的逐级筛选 ──
    def _intv(s):
        try: return int(s.replace(',',''))
        except: return None
    flow_stages = []
    r_raw = _intv(kpis.get('raw_reads',''))
    r_clean = _intv(kpis.get('clean_reads',''))
    r_contig = _intv(kpis.get('total_contigs',''))
    r_virus = _intv(kpis.get('virus_seqs',''))
    r_votu = _intv(kpis.get('n_clusters',''))
    r_hq = _intv(kpis.get('hq_votus',''))
    # host-free reads: 从 hostdep_summary 汇总 After_Host
    r_hostfree = None
    hds = _read_tsv(report_dir / "hostdep_summary.tsv")
    if hds:
        s = sum(int(r.get("After_Host",0)) for r in hds)
        if s > 0: r_hostfree = s
    # 序数色带（raw → HQ vOTUs 逐级筛选）：viridis 7 档，与静态图 SEQ_CMAP 同源
    _flow_c = _ramp(7)
    if r_raw: flow_stages.append(("Raw Reads",r_raw,_flow_c[0]))
    if r_clean: flow_stages.append(("Clean Reads",r_clean,_flow_c[1]))
    if r_hostfree: flow_stages.append(("Host-free Reads",r_hostfree,_flow_c[2]))
    if r_contig: flow_stages.append(("Contigs",r_contig,_flow_c[3]))
    if r_virus: flow_stages.append(("Viral Seqs",r_virus,_flow_c[4]))
    if r_votu: flow_stages.append(("vOTU Clusters",r_votu,_flow_c[5]))
    if r_hq: flow_stages.append(("HQ vOTUs",r_hq,_flow_c[6]))

    flow_html = ""
    ai_pager_js = """
function toggleChat(){var p=document.getElementById("chat-panel");p.style.display=p.style.display==="none"?"flex":"none";}
function toggleAIPanel(){var p=document.getElementById("ai-panel");p.style.display=p.style.display==="none"?"block":"none";}
function onProviderChange(){var p=document.getElementById("ai-provider").value;var m={openai:"gpt-4o-mini",deepseek:"deepseek-v4-flash",ollama:"qwen2.5:7b",moonshot:"moonshot-v1-8k",custom:""};document.getElementById("ai-model").value=m[p]||"";}
function getAIURL(){var p=document.getElementById("ai-provider").value;if(p==="openai")return"https://api.openai.com/v1/chat/completions";if(p==="deepseek")return"https://api.deepseek.com/v1/chat/completions";if(p==="moonshot")return"https://api.moonshot.cn/v1/chat/completions";if(p==="ollama")return"http://localhost:11434/v1/chat/completions";return prompt("Enter API base URL:")||"";}
function sendChat(){var inp=document.getElementById("chat-input");var q=inp.value.trim();if(!q)return;var key=document.getElementById("ai-api-key").value||sessionStorage.getItem("ai_key");if(!key){toggleAIPanel();alert("Please enter your API key");return;}sessionStorage.setItem("ai_key",key);var msgs=document.getElementById("chat-msgs");msgs.innerHTML+="<div class=chat-msg user>"+q.replace(/</g,"&lt;")+"</div>";inp.value="";msgs.scrollTop=msgs.scrollHeight;var th=document.createElement("div");th.className="chat-msg ai";th.textContent="...";msgs.appendChild(th);var ctx=document.getElementById("stage-data")?document.getElementById("stage-data").textContent:"{}";fetch(getAIURL(),{method:"POST",headers:{"Content-Type":"application/json",Authorization:"Bearer "+key},body:JSON.stringify({model:document.getElementById("ai-model").value,messages:[{role:"system",content:"You are a virology research assistant. Answer concisely in English using the report context: "+ctx},{role:"user",content:q}],temperature:0.3,max_tokens:500})}).then(function(r){return r.json()}).then(function(d){th.textContent=d.choices[0].message.content.trim();msgs.scrollTop=msgs.scrollHeight;}).catch(function(e){th.textContent="Error: "+e;});}
function runStageAI(sn){var key=document.getElementById("ai-api-key").value||sessionStorage.getItem("ai_key");if(!key){toggleAIPanel();alert("Please enter your API key");return;}sessionStorage.setItem("ai_key",key);var stageInfo=_stageData[sn]||{};var sel=stageInfo.short||sn;var btn=document.querySelector("#stage-"+sel+" .ai-stage-btn");if(!btn)return;var orig=btn.textContent;btn.textContent="...";btn.disabled=true;var ctx=document.getElementById("stage-data")?document.getElementById("stage-data").textContent:"{}";fetch(getAIURL(),{method:"POST",headers:{"Content-Type":"application/json",Authorization:"Bearer "+key},body:JSON.stringify({model:document.getElementById("ai-model").value,messages:[{role:"system",content:"You are a virology scientist. Summarize this pipeline stage concisely in English (2-3 sentences)."},{role:"user",content:"Stage: "+sn+"\\nData: "+ctx}],temperature:0.3,max_tokens:300})}).then(function(r){return r.json()}).then(function(d){btn.textContent=orig;btn.disabled=false;var t=document.querySelector("#stage-"+sel+" .ai-stage-summary");if(!t){t=document.createElement("p");t.className="ai-stage-summary";t.style.display="block";var h=document.querySelector("#stage-"+sel+" .stage-charts, #stage-"+sel+" .stage-tables");if(h)h.before(t);else{var s=document.querySelector("#stage-"+sel+" .stage-header");if(s)s.after(t);}}t.textContent=d.choices[0].message.content.trim();}).catch(function(e){btn.textContent=orig;btn.disabled=false;alert("Error: "+e);});}
var _stageData=JSON.parse(document.getElementById("stage-data").textContent);var _saved=sessionStorage.getItem("ai_key");if(_saved)document.getElementById("ai-api-key").value=_saved;
document.getElementById("ai-provider").value=sessionStorage.getItem("ai_provider")||"openai";
document.getElementById("ai-model").value=sessionStorage.getItem("ai_model")||"gpt-4o-mini";
"""
    if len(flow_stages) >= 3:
        f_labels = [s[0] for s in flow_stages]
        f_values = [s[1] for s in flow_stages]
        f_colors = [s[2] for s in flow_stages]
        flow_chart_id = "chart_pipeline_flow"
        # Chart.js v4 的 tooltip 回调必须是真函数: 函数体放进 python 字典会被
        # json.dumps 变成字符串 ("function..." 不执行)。先用唯一占位符走 _chart,
        # 出 JS 后再把引号包裹的占位替换为函数字面量。
        chart_scripts += _chart(flow_chart_id, 'bar', {
            "labels": f_labels,
            "datasets": [{"label":"Count","data":f_values,"backgroundColor":f_colors,
                          "borderColor":f_colors,"borderWidth":0}]},
            {"indexAxis":"y","responsive":True,
             "plugins":{"title":{"display":True,"text":"Pipeline Flow — Reads → HQ vOTUs","font":{"size":14}},
                        "legend":{"display":False},
                        "tooltip":{"callbacks":{"label":"__FLOW_LABEL_CB__"}}},
             "scales":{"x":{"type":"logarithmic","title":{"text":"Count (log scale)","display":True}}}}
        ).replace('"__FLOW_LABEL_CB__"',
                  "function(ctx){var v=ctx.raw;if(v>=1e6)return (v/1e6).toFixed(1)+' M';"
                  "if(v>=1e3)return (v/1e3).toFixed(1)+' K';return v}")
        flow_html = f'<div class="flow-section"><div class="chart-box" style="max-width:800px;margin:0 auto"><canvas id="{flow_chart_id}" style="max-height:380px"></canvas></div></div>'

    html = f'''<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width,initial-scale=1.0">
<!-- impeccable-disable numbered-section-markers,em-dash-overuse: stage IDs are pipeline identifiers; CJK scientific copy -->
<title>MMPV-RNA {VERSION} — Pipeline Report</title>
<script>/*! Chart.js | https://www.chartjs.org | MIT License */
{{CHART_JS_INLINE}}</script>
<script>/*! UpSet.js | https://upset.js.org */
{{UPSET_JS_INLINE}}</script>
<style>
:root{{
  --bg:#F8F4ED;--card-bg:#FCFAF5;--text:#3B3D3F;--muted:#6B6F73;
  --navy:#1B365D;--indigo:#1B365D;--blue:#537D96;--green:#4A6B4A;
  --red:#8B3A3A;--amber:#9D5F4D;--border:rgba(122,96,88,.18);
  --shadow:0 1px 3px rgba(60,48,40,.05);--shadow-lg:0 3px 10px rgba(60,48,40,.08);
  --radius:4px;--radius-sm:2px;
}}
*{{box-sizing:border-box;margin:0;padding:0}}
body{{font-family:'EB Garamond','Noto Serif SC','Source Han Serif SC','Songti SC','STSong',serif;background:var(--bg);color:var(--text);line-height:1.7}}
.container{{max-width:1200px;margin:0 auto;padding:0 20px 40px;padding-left:240px}}
.hero{{background:#1B365D;color:#FFFDF7;padding:40px 32px 32px;border-radius:0 0 8px 8px;margin-bottom:28px;position:relative;overflow:hidden}}
.hero h1{{font-size:26px;font-weight:500;margin-bottom:6px;position:relative;z-index:1}}
.hero .subtitle{{font-size:13px;opacity:.85;position:relative;z-index:1}}
.hero .gen-time{{font-size:12px;opacity:.7;margin-top:6px;position:relative;z-index:1}}
.kpi-row{{display:grid;grid-template-columns:repeat(auto-fit,minmax(180px,1fr));gap:14px;margin-bottom:20px}}
.flow-section{{background:var(--card-bg);border-radius:var(--radius);box-shadow:var(--shadow);padding:16px 22px;margin-bottom:28px}}
.kpi-card{{background:var(--card-bg);border-radius:var(--radius);padding:18px 20px;box-shadow:var(--shadow);text-align:center;border:1px solid var(--border)}}
.kpi-icon{{font-size:24px;margin-bottom:6px}}
.kpi-value{{font-size:13px;font-weight:600;color:var(--text);line-height:1.4}}
.kpi-label{{font-size:11px;color:var(--muted);margin-top:4px;text-transform:uppercase;letter-spacing:.5px}}
.sidebar{{position:fixed;top:0;left:0;width:224px;height:100vh;background:var(--card-bg);border-right:1px solid var(--border);box-shadow:2px 0 8px rgba(0,0,0,.04);z-index:200;display:flex;flex-direction:column;padding-top:12px;overflow-y:auto}}
.sb-title{{font-size:11px;font-weight:700;color:var(--muted);text-transform:uppercase;letter-spacing:.8px;padding:8px 16px 6px;margin:0}}
.sb-item{{display:flex;align-items:center;gap:10px;padding:9px 16px 9px 19px;text-decoration:none;color:var(--text);font-size:12.5px;font-weight:500;transition:all .15s}}
.sb-item:hover{{background:rgba(83,125,150,.08);color:#3F6179}}
.sb-dot{{width:10px;height:10px;border-radius:50%;flex-shrink:0}}
.sb-label{{overflow:hidden;text-overflow:ellipsis;white-space:nowrap}}
.sb-sub{{display:flex;align-items:center;gap:8px;padding:6px 16px 6px 37px;text-decoration:none;color:var(--muted);font-size:11.5px;font-weight:400;transition:all .15s}}
.sb-sub:hover{{color:var(--text);background:rgba(83,125,150,.06)}}
.sb-dot-sub{{width:6px;height:6px;border-radius:50%;flex-shrink:0;background:var(--blue)}}
#keep-table th:first-child,#keep-table td:first-child{{position:sticky;left:0;background:#FCFAF5;z-index:1;border-right:1px solid var(--border)}}
.sb-pass{{background:rgba(74,107,74,.05)}}.sb-fail{{background:rgba(139,58,58,.05)}}.sb-skip{{background:rgba(107,111,115,.035)}}
.sb-active{{background:rgba(83,125,150,.1);font-weight:500;color:#3F6179}}
.stage{{background:var(--card-bg);border-radius:var(--radius);box-shadow:var(--shadow);margin-bottom:20px;overflow:hidden}}
.stage-pass{{border:1px solid rgba(74,107,74,.38);background:rgba(74,107,74,.045)}}.stage-fail{{border:1px solid rgba(139,58,58,.38);background:rgba(139,58,58,.045)}}.stage-skip{{border:1px solid var(--border);background:rgba(107,111,115,.04)}}
.stage-header{{display:flex;align-items:center;gap:16px;padding:18px 22px;border-bottom:1px solid var(--border);background:var(--bg)}}
.stage-icon{{width:44px;height:44px;border-radius:var(--radius-sm);display:flex;align-items:center;justify-content:center;font-size:11px;font-weight:500;color:#FFFDF7;flex-shrink:0}}
.stage-title{{flex:1;min-width:0}}
.stage-title h2{{font-size:16px;font-weight:500;color:var(--text);margin-bottom:2px}}
.stage-metric{{font-size:12px;color:var(--muted);display:block;overflow:hidden;text-overflow:ellipsis;white-space:nowrap}}
.stage-badge{{font-size:11px;font-weight:500;padding:4px 12px;border-radius:var(--radius-sm);flex-shrink:0}}
.s-pass{{background:rgba(74,107,74,.12);color:var(--green)}}.s-fail{{background:rgba(139,58,58,.12);color:var(--red)}}.s-skip{{background:rgba(0,0,0,.05);color:var(--muted)}}
.stage-charts{{display:grid;gap:16px;padding:18px 22px;grid-template-columns:repeat(auto-fit,minmax(0,1fr))}}
.chart-box{{background:var(--bg);border-radius:var(--radius-sm);padding:12px;border:1px solid var(--border);position:relative;min-width:0}}
.chart-scroll{{overflow-x:auto;max-width:100%}}
.pct-btn{{font-size:11px;padding:4px 14px;border:1px solid var(--blue);border-radius:var(--radius-sm);background:var(--blue);cursor:pointer;color:#FFFDF7;font-weight:500;transition:all .15s}}
.pct-btn:hover{{background:var(--accent-hover,#3F6179);border-color:#3F6179}}
.sankey-section{{padding:18px 22px;display:flex;flex-direction:column;gap:16px}}
.sankey-card{{background:var(--bg);border-radius:var(--radius-sm);padding:14px;border:1px solid var(--border)}}
.sankey-card h3{{font-size:14px;color:var(--indigo);margin-bottom:8px;text-align:center;font-weight:500}}
.sankey-card iframe{{display:block;width:100%;border-radius:var(--radius-sm)}}
.stage-tables{{padding:0 22px 14px;display:flex;flex-direction:column;gap:8px}}
.stage-table-detail{{border-top:1px solid var(--border);padding:8px 0 4px}}
.stage-table-detail summary{{cursor:pointer;font-size:12px;padding:2px 0;color:var(--muted)}}
.stage-table-detail summary:hover{{color:var(--blue)}}
.table-wrap{{background:var(--card-bg);border-radius:var(--radius);box-shadow:var(--shadow);overflow:hidden;margin-bottom:20px}}
.table-wrap h3{{font-size:15px;padding:16px 22px;border-bottom:1px solid var(--border);color:var(--text);font-weight:500}}
table{{width:100%;border-collapse:collapse;font-size:13px}}
th{{background:var(--bg);text-align:left;padding:10px 16px;font-size:11px;color:var(--muted);text-transform:uppercase;letter-spacing:.5px;border-bottom:2px solid var(--border)}}
td{{padding:9px 16px;border-bottom:1px solid var(--border)}}
.tr-pass{{background:rgba(74,107,74,.07)}}.tr-fail{{background:rgba(139,58,58,.07)}}.tr-skip{{background:rgba(107,111,115,.05)}}
.td-num{{text-align:center;color:var(--muted);font-size:11px;width:36px}}
.td-sub{{padding-left:40px!important;color:var(--muted);font-size:12px}}
.td-metric{{font-size:12px;color:var(--muted)}}
.tb-badge{{display:inline-block;padding:1px 8px;border-radius:var(--radius-sm);font-size:11px;font-weight:500}}
.tb-pass{{background:rgba(74,107,74,.12);color:var(--green)}}.tb-fail{{background:rgba(139,58,58,.12);color:var(--red)}}.tb-skip{{background:rgba(0,0,0,.05);color:var(--muted)}}
.app-table{{width:100%;border-collapse:collapse;font-size:11px;margin-bottom:8px}}
.app-table th{{background:#f5f5f5;padding:6px 8px;font-size:11px;text-align:left;border:1px solid #e0e0e0;white-space:nowrap}}
.app-table td{{padding:4px 8px;border:1px solid #f0f0f0;max-width:200px;overflow:hidden;text-overflow:ellipsis;white-space:nowrap;font-size:11px}}
.ai-section{{background:var(--card-bg);border-radius:var(--radius);box-shadow:var(--shadow);padding:20px 24px;margin-bottom:28px;border:1px solid var(--border)}}
.ai-header{{display:flex;align-items:center;gap:10px;margin-bottom:14px;border-bottom:1px solid var(--border);padding-bottom:10px}}
.ai-title{{font-size:15px;font-weight:500;color:var(--blue)}}
.ai-model-tag{{font-size:11px;background:var(--blue);color:#FFFDF7;padding:2px 8px;border-radius:var(--radius-sm)}}
.ai-format-tag{{font-size:11px;background:rgba(83,125,150,.12);color:var(--blue);padding:2px 8px;border-radius:var(--radius-sm);font-weight:500}}
.ai-brief{{font-size:14px;line-height:1.8;color:var(--text);background:var(--bg);padding:14px 18px;border-radius:var(--radius-sm);margin-bottom:14px;text-align:justify;font-weight:400}}
.ai-detail-toggle{{margin-top:8px}}
.ai-detail-toggle summary{{cursor:pointer;font-size:13px;font-weight:500;color:var(--blue);padding:10px 14px;background:var(--bg);border-radius:var(--radius-sm);user-select:none;transition:all .15s}}
.ai-detail-toggle summary:hover{{background:rgba(83,125,150,.1)}}
.ai-detail-content{{margin-top:14px}}
.ai-block{{margin-bottom:14px;padding:12px 14px;border:1px solid var(--border);border-radius:var(--radius-sm);background:rgba(83,125,150,.04)}}
.ai-block:last-child{{margin-bottom:0}}
.ai-block-title{{font-size:12px;font-weight:500;color:var(--blue);margin-bottom:4px}}
.ai-block-text{{font-size:13px;line-height:1.9;color:var(--text);text-align:justify}}
.ai-bg{{border-color:rgba(83,125,150,.45)}}.ai-methods{{border-color:rgba(74,107,74,.45)}}.ai-results{{border-color:rgba(157,95,77,.45)}}.ai-discussion{{border-color:rgba(107,111,115,.45)}}
.ai-error{{color:var(--red);font-size:12px;padding:8px}}
.ai-stage-btn{{font-size:11px;font-weight:500;padding:1px 8px;margin-left:8px;background:#537D96;color:#fff;border:none;border-radius:var(--radius-sm);cursor:pointer;vertical-align:middle;opacity:.85}}
.ai-stage-btn:hover{{opacity:1}}
.ai-stage-summary{{font-size:12px;line-height:1.6;padding:8px 12px;margin:6px 22px;background:var(--bg);border-radius:var(--radius-sm);display:none}}
#chat-btn{{position:fixed;bottom:72px;right:24px;width:38px;height:38px;background:#537D96;color:#fff;border:none;border-radius:50%;font-size:16px;font-weight:500;cursor:pointer;z-index:200;box-shadow:0 2px 8px rgba(60,48,40,.15);transition:background .15s}}
#chat-btn:hover{{background:#3F6179}}
#chat-panel{{position:fixed;bottom:120px;right:24px;width:360px;max-height:500px;background:var(--card-bg);border:1px solid var(--border);border-radius:var(--radius);box-shadow:var(--shadow-lg);z-index:199;display:none;flex-direction:column}}
.chat-header{{display:flex;justify-content:space-between;align-items:center;padding:8px 12px;border-bottom:1px solid var(--border);font-size:13px;font-weight:500}}
.chat-messages{{flex:1;overflow-y:auto;padding:8px 10px;max-height:360px;display:flex;flex-direction:column;gap:6px}}
.chat-msg{{font-size:12px;line-height:1.5;padding:6px 10px;border-radius:var(--radius-sm);max-width:85%}}
.chat-msg.user{{align-self:flex-end;background:var(--blue);color:#fff}}
.chat-msg.ai{{align-self:flex-start;background:var(--bg);color:var(--text)}}
.chat-input-row{{display:flex;gap:6px;padding:8px 10px;border-top:1px solid var(--border)}}
.chat-input-row input{{flex:1;padding:6px 8px;border:1px solid var(--border);border-radius:var(--radius-sm);font-size:12px}}
.chat-send-btn{{padding:6px 12px;background:var(--blue);color:#fff;border:none;border-radius:var(--radius-sm);font-size:12px;cursor:pointer}}
#ai-btn{{position:fixed;bottom:72px;right:72px;width:38px;height:38px;background:#537D96;color:#fff;border:none;border-radius:50%;font-size:13px;font-weight:500;cursor:pointer;z-index:200;box-shadow:0 2px 8px rgba(60,48,40,.15);transition:background .15s}}
#ai-btn:hover{{background:#3F6179}}
#ai-panel{{position:fixed;bottom:120px;right:24px;width:300px;background:var(--card-bg);border:1px solid var(--border);border-radius:var(--radius);box-shadow:var(--shadow-lg);z-index:198;display:none}}
.footer{{text-align:center;padding:24px;color:var(--muted);font-size:12px;line-height:1.8}}
.footer a{{color:var(--blue);text-decoration:none}}
@media(prefers-reduced-motion:reduce){{.kpi-card,.sb-item,.pct-btn,#chat-btn,#ai-btn,.stage-table-detail summary{{transition:none!important}}
.chart-box img{{transition:none!important}}}}
/* 断点 560px（2026-10-17 用户裁定：侧栏常驻，不再随宽度隐藏） */
@media(max-width:560px){{
  .hero{{padding:24px 20px 20px}}.hero h1{{font-size:20px}}
  .kpi-row{{grid-template-columns:repeat(2,1fr)}}
  .stage-header{{flex-wrap:wrap;gap:10px}}
  .stage-charts{{grid-template-columns:1fr!important}}
  .footer{{font-size:12px}}}}
@media print{{
  body{{background:#fff;font-size:11px}}
  .hero{{background:#1B365D!important;-webkit-print-color-adjust:exact}}
  .stage,.table-wrap,.kpi-card{{box-shadow:none;border:1px solid var(--border);break-inside:avoid}}
  .sidebar{{display:none}}
  .container{{padding-left:20px}}
}}
{_TW_CSS}
</style>
</head>
<body>
<div id="lightbox" onclick="closeLightbox()" style="display:none;position:fixed;top:0;left:0;width:100%;height:100%;background:rgba(0,0,0,.85);z-index:999;cursor:zoom-out;align-items:center;justify-content:center">
  <img id="lightbox-img" src="data:image/gif;base64,R0lGODlhAQABAIAAAAAAAP///yH5BAEAAAAALAAAAAABAAEAAAIBRAA7" alt="" style="max-width:95%;max-height:95%;object-fit:contain" onclick="event.stopPropagation()">
</div>
<button id="chat-btn" onclick="toggleChat()" title="Ask questions about this report">?</button>
<div id="chat-panel" style="display:none">
    <div class="chat-header">
        <span>Report Q&amp;A</span>
        <button onclick="toggleChat()" style="background:none;border:none;color:var(--muted);font-size:16px;cursor:pointer">&times;</button>
    </div>
    <div class="chat-messages" id="chat-msgs"></div>
    <div class="chat-input-row">
        <input id="chat-input" type="text" placeholder="Ask about this report..." onkeydown="if(event.key==='Enter')sendChat()">
        <button onclick="sendChat()" class="chat-send-btn">Send</button>
    </div>
</div>
<button id="ai-btn" onclick="toggleAIPanel()" title="AI Settings">AI</button>
<div id="ai-panel" style="display:none">
    <div style="display:flex;justify-content:space-between;align-items:center;padding:10px 14px;border-bottom:1px solid var(--border)">
        <span style="font-weight:700;font-size:14px;color:#1a237e">AI Settings</span>
        <button onclick="toggleAIPanel()" style="background:none;border:none;font-size:18px;cursor:pointer;color:var(--muted)">&times;</button>
    </div>
    <div style="padding:10px 14px">
        <label style="font-size:11px;color:var(--muted)">API Key (saved in this browser)</label>
        <input id="ai-api-key" type="password" placeholder="sk-..." style="width:100%;padding:6px 8px;margin:4px 0 8px;border:1px solid var(--border);border-radius:4px;font-size:12px" autocomplete="off">
        <label style="font-size:11px;color:var(--muted)">Provider</label>
        <select id="ai-provider" onchange="onProviderChange()" style="width:100%;padding:6px 8px;margin:4px 0 8px;border:1px solid var(--border);border-radius:4px;font-size:12px">
            <option value="openai">OpenAI</option>
            <option value="deepseek">DeepSeek</option>
            <option value="moonshot">Kimi (Moonshot)</option>
            <option value="ollama">Ollama (localhost)</option>
            <option value="custom">Custom</option>
        </select>
        <label style="font-size:11px;color:var(--muted)">Model</label>
        <input id="ai-model" value="gpt-4o-mini" style="width:100%;padding:6px 8px;margin:4px 0 8px;border:1px solid var(--border);border-radius:4px;font-size:12px">
    </div>
</div>
<div class="sidebar">
  <div class="sb-title">Pipeline Modules</div>
  {sidebar_items}
</div>
<div class="container">
<div class="hero">
  <h1>MMPV-RNA {VERSION} — Pipeline Report</h1>
  <div class="subtitle">Metatranscriptomic Virus Discovery — End-to-End Analysis</div>
  <div class="gen-time">Generated: {gen_time} &nbsp;|&nbsp; {n_pass}/{n_total} stages completed ({pct}%)</div>
</div>
<div class="kpi-row">{kpi_cards}</div>
{flow_html}
{sections_html}
<div class="table-wrap">
  <div style="display:flex;justify-content:space-between;align-items:center;padding-right:16px">
    <h3>Pipeline Stage Summary</h3>
    <button onclick="exportTable()" style="background:var(--indigo);color:#fff;border:none;padding:6px 16px;border-radius:6px;cursor:pointer;font-size:12px">Export CSV</button>
  </div>
  <div id="table-scroll" style="max-height:520px;overflow-y:auto">
  <table id="summary-table"><thead><tr><th style="width:36px">#</th><th>Stage</th><th style="width:70px">Status</th><th>Key Metrics</th></tr></thead>
  <tbody>{table_rows}</tbody></table>
  </div>
  <div id="table-pager" style="display:flex;justify-content:center;align-items:center;gap:8px;padding:12px;border-top:1px solid var(--border);font-size:12px;color:var(--muted)"></div>
</div>
{ai_summary_html}
<div class="footer">
  <strong>MMPV-RNA {VERSION}</strong> — Metatranscriptomic Virus Discovery Pipeline<br>
  Generated by <code>report_pipeline.py</code> &nbsp;|&nbsp; {gen_time}
</div>
</div>
<script id="stage-data" type="application/json">{stage_data_json}</script>
<script id="upset-data" type="application/json">{upset_data_json}</script>
<script>
{chart_scripts}
{upset_render_js}
{sankey_inject_scripts}
(function(){{
  const tbody=document.querySelector('#summary-table tbody');
  if(!tbody)return;
  const rows=Array.from(tbody.querySelectorAll('tr'));
  const perPage=10;
  const totalPages=Math.ceil(rows.length/perPage);
  if(totalPages<=1)return;
  const pager=document.getElementById('table-pager');
  if(!pager)return;
  let page=0;
  function show(p){{
    page=Math.max(0,Math.min(p,totalPages-1));
    rows.forEach((r,i)=>{{r.style.display=(i>=page*perPage&&i<(page+1)*perPage)?'':'none'}});
    pager.innerHTML='<button onclick="window._tblPage('+(page-1)+')" '+(page===0?'disabled':'')+' style="border:1px solid #ccc;background:#fff;padding:4px 12px;border-radius:4px;cursor:pointer">← Prev</button>'+
      '<span>Page <b>'+(page+1)+'</b> of '+totalPages+'</span>'+
      '<button onclick="window._tblPage('+(page+1)+')" '+(page===totalPages-1?'disabled':'')+' style="border:1px solid #ccc;background:#fff;padding:4px 12px;border-radius:4px;cursor:pointer">Next →</button>';
  }}
  window._tblPage=function(p){{show(p)}};
  show(0);
}})();
{_esf_js}
{_keep_js}
function openLightbox(src){{document.getElementById('lightbox-img').src=src;document.getElementById('lightbox').style.display='flex';}}
function closeLightbox(){{document.getElementById('lightbox').style.display='none';}}
function exportTable(){{
  const tbl=document.getElementById('summary-table');
  if(!tbl)return;
  let csv='';
  tbl.querySelectorAll('tr').forEach(tr=>{{
    let row=[];
    tr.querySelectorAll('th,td').forEach(cell=>row.push('"'+cell.innerText.replace(/"/g,'""')+'"'));
    csv+=row.join(',')+'\\n';
  }});
  const blob=new Blob([csv],{{type:'text/csv'}});
  const a=document.createElement('a');
  a.href=URL.createObjectURL(blob);
  a.download='pipeline_summary.csv';
  a.click();
}}
{table_pager_js}
// Sidebar scroll-spy
(function(){{
  const items=document.querySelectorAll('.sb-item');
  const stages=document.querySelectorAll('.stage[id]');
  if(!items.length||!stages.length)return;
  function onScroll(){{
    let current='';
    stages.forEach(s=>{{if(s.getBoundingClientRect().top<=160)current=s.id}});
    items.forEach(a=>{{
      a.classList.toggle('sb-active',a.getAttribute('href')==='#'+current);
    }});
  }}
  window.addEventListener('scroll',onScroll,{{passive:true}});
  onScroll();
}})();
// Stacked bar chart % toggle
(function(){{
  var _orig={{}};
  window.toggleStackedPct=function(btn){{
    var cid=btn.getAttribute('data-chart');
    var chart=Chart.getChart(cid);
    if(!chart)return;
    var ds=chart.data.datasets;
    var n=chart.data.labels.length;
    // Value axis: horizontal charts (indexAxis='y') have the value axis on x; vertical charts on y
    var valAxis=chart.options.indexAxis==='y'?'x':'y';
    if(btn.getAttribute('data-mode')==='abs'){{
      if(!_orig[cid])_orig[cid]=ds.map(function(d){{return d.data.slice()}});
      for(var i=0;i<n;i++){{
        var tot=0;
        for(var j=0;j<ds.length;j++)tot+=Number(_orig[cid][j][i])||0;
        for(var j=0;j<ds.length;j++)ds[j].data[i]=tot>0?((Number(_orig[cid][j][i])||0)/tot*100):0;
      }}
      chart.options.scales[valAxis].title={{display:true,text:'%'}};
      chart.options.scales[valAxis].max=100;chart.options.scales[valAxis].min=0;
      btn.setAttribute('data-mode','pct');btn.textContent='Show Count';
    }}else{{
      if(_orig[cid])for(var j=0;j<ds.length;j++)ds[j].data=_orig[cid][j].slice();
      chart.options.scales[valAxis].title={{display:true,text:'Reads'}};
      chart.options.scales[valAxis].max=undefined;chart.options.scales[valAxis].min=0;
      btn.setAttribute('data-mode','abs');btn.textContent='Show %';
    }}
    chart.update();
  }};
}})();
</script>
<script>
{ai_pager_js}
</script>
</body></html>'''

    with open(report_dir / "pipeline_report.html", "w", encoding="utf-8") as hf:
        hf.write(_en(html.replace("{CHART_JS_INLINE}", chart_js_inline).replace("{UPSET_JS_INLINE}", upset_js_inline).replace("{UPSET_REACT_INLINE}", upset_react_inline)))


# ═══════════════════════════════════════════════════════════════
# Main
# ═══════════════════════════════════════════════════════════════

def generate_ai_summary(stage_stats, kpis, report_dir):
    """调用 LLM 生成双层 AI 总结: 简报 + IMRaD 详报, 返回 HTML 或空"""
    import json as _json
    import re

    # 收集管线数据
    main_stages = [s for s in stage_stats if not s["Stage"].startswith("  ")]
    stage_lines = []
    for s in main_stages:
        st = s["Status"]; badge = "✓" if st == "✓" else ("✗" if st == "✗" else "○")
        stage_lines.append(f"  {badge} {s['Stage']}: {s.get('Key_Metric','')}")

    # 宿主分布 + CheckV 详情
    host_detail = ""
    cv_detail = ""
    for s in stage_stats:
        if s["Stage"].startswith("  ") and s["Stage"].strip() in ("Bacteria","Protist","Animal","Plant","Algae","Mammalia","Fungi","Unknown"):
            host_detail += f"    {s['Stage'].strip()}: {s.get('Key_Metric','')}\n"
        if s["Stage"].startswith("  ") and "HQ=" in s.get("Key_Metric",""):
            cv_detail += f"    {s['Stage'].strip()}: {s.get('Key_Metric','')}\n"

    prompt = f"""你是一位病毒宏基因组学领域的研究科学家。请基于以下 MMPV-RNA {VERSION} 病毒发现管线的分析结果, 撰写双层研究总结: 先写一段精炼的研究简报, 再写一份完整的 IMRaD 详报。严格按标签输出。

## 管线数据
- 样本数: {kpis.get('n_sample','?')}
- 数据量: QC前 {kpis.get('raw_bases','?')} → QC后 {kpis.get('clean_bases','?')}
- 宿主去除总保留率: {kpis.get('hd_retained','?')}%
- 各阶段结果:
{chr(10).join(stage_lines)}
- 宿主分布 (各宿主HQ/total):
{host_detail if host_detail else '  见各阶段'}
- CheckV 各宿主质量:
{cv_detail if cv_detail else '  见各阶段'}
- 新颖性: 已知={kpis.get('n_known','?')} | 新种={kpis.get('n_newsp','?')} | 新属={kpis.get('n_newge','?')} | 新科={kpis.get('n_newfa','?')}
- HQ vOTUs: {kpis.get('hq_votus','?')} / {kpis.get('cv_total','?')}

## 输出格式 (严格按此结构)

[BRIEF]
A 150-200 word English research highlight. Cover: study purpose, sample scale, key findings (viruses identified, novelty ratio, HQ vOTU count), main conclusions. Highlight plant viruses if present. Style: top-journal Highlights.

[DETAILED]
Write the detailed report in the following IMRaD structure (400-600 words, English):

[BACKGROUND]
1-2句。病毒宏基因组学背景 + 本研究目标。

[METHODS]
4-6句。完整列出 MMPV-RNA {VERSION} 全部阶段:
- 00a: fastp 质控 (Q20/Q30过滤) + clumpify 聚类重排 → 高质量 clean reads
- 00b: Kraken2 去除宿主 reads + ribodetector 去除 rRNA → 非宿主非核糖体 reads
- 01: SPAdes/MEGAHIT 从头组装 → contigs
- 02: 多工具病毒鉴定 (Viralm, BLASTx/n, VirHunter, CAT, genomad, metabuli, mmseqs, vcontact3) → 病毒候选序列
- 03: COBRA 延伸 (BLAST-based contig extension)
- 04: CD-HIT 聚类 (去冗余) + vclust → vOTU 簇
- 05: 分类学注释 (CAT + VITAP + ACVirus 集成) → 科/属/种级别分类
- 06: 宿主预测 (ICTV + VITAP + CAT + ACVirus + BLAST 共识决策) → Final_Host
- 07: CheckV 完整性评估 (AAI + HMM) → HQ/MQ/LQ 分级
- 08: Rescue (从低质量序列中恢复 HQ vOTU)
每步说明目的, 用 → 连接。

[RESULTS]
4-6句。报告数值发现: 数据量与QC保留率、组装contigs数和总长、病毒序列数与vOTU簇数、新颖性分布 (已知/新种/新属/新科) 及比例、宿主分布 (如有植物病毒需重点描述: 数量、分类层级、完整性)、CheckV 质量分布 (Complete/HQ/MQ/LQ)、Rescue结果。

[DISCUSSION]
3-5句。解读意义: 新颖性比例的生物学含义、宿主分布的特征趋势 (噬菌体 vs 真核病毒)、植物病毒的发现意义和完整性。管线优势与局限。下一步: 功能注释、比较基因组、系统发育、宿主-病毒互作网络。

## 要求
- 专业学术语气, 精确引用数据中的数字
- Brief 后空一行再输出 DETAILED 部分
- DETAILED 内四段以 [BACKGROUND][METHODS][RESULTS][DISCUSSION] 开头
- 仅输出要求的格式, 不要额外说明"""

    try:
        import html as _html
        provider = getattr(generate_ai_summary, '_provider', 'openai')
        model = getattr(generate_ai_summary, '_model', 'gpt-4o-mini')
        api_key = getattr(generate_ai_summary, '_api_key', '')
        base_url = getattr(generate_ai_summary, '_base_url', '')

        if not api_key:
            return '<div class="ai-error">⚠ AI summary requires the --ai-key option</div>'

        import urllib.request
        url = base_url or "https://api.openai.com/v1/chat/completions"
        body = _json.dumps({
            "model": model,
            "messages": [
                {"role": "system", "content": "你是病毒宏基因组学研究科学家。严格按要求的双层格式输出: BRIEF(简洁亮点摘要)+DETAILED(IMRaD完整详报)。in English, objective and precise, cite actual numbers, do not fabricate."},
                {"role": "user", "content": prompt}
            ],
            "temperature": 0.3,
            "max_tokens": 4000
        }).encode()
        req = urllib.request.Request(url, data=body, headers={
            "Content-Type": "application/json",
            "Authorization": f"Bearer {api_key}"
        })
        with urllib.request.urlopen(req, timeout=90) as resp:
            result = _json.loads(resp.read())
        raw = result["choices"][0]["message"]["content"].strip()
        print(f"  AI Summary generated ({len(raw)} chars)")

        # 解析 Brief 和 Detailed
        brief_m = re.search(r'\[BRIEF\]\s*(.*?)(?=\[DETAILED\]|\Z)', raw, re.DOTALL | re.IGNORECASE)
        brief = brief_m.group(1).strip() if brief_m else ""

        detailed_raw = ""
        det_m = re.search(r'\[DETAILED\]\s*(.*)', raw, re.DOTALL | re.IGNORECASE)
        if det_m: detailed_raw = det_m.group(1).strip()

        def _extract_section(text, tag):
            m = re.search(
                rf'\[{tag}\]\s*(.*?)(?=\[(?:BACKGROUND|METHODS|RESULTS|DISCUSSION)\]|\Z)',
                text, re.DOTALL | re.IGNORECASE)
            return m.group(1).strip() if m else ""

        sections = {}
        if detailed_raw:
            sections = {
                "background": _extract_section(detailed_raw, "BACKGROUND"),
                "methods": _extract_section(detailed_raw, "METHODS"),
                "results": _extract_section(detailed_raw, "RESULTS"),
                "discussion": _extract_section(detailed_raw, "DISCUSSION"),
            }
        if not any(sections.values()):
            sections = {"background": detailed_raw or raw, "methods": "", "results": "", "discussion": ""}

        # 构建 HTML: Brief + 可折叠 Detail
        sec_labels = {
            "background": ("📖 Background", "ai-bg"),
            "methods": ("🔬 Methods", "ai-methods"),
            "results": ("📊 Results", "ai-results"),
            "discussion": ("💡 Discussion", "ai-discussion"),
        }
        detail_items = ""
        for key in ("background","methods","results","discussion"):
            text = sections.get(key, "")
            if not text: continue
            label, cls = sec_labels[key]
            detail_items += f'<div class="ai-block {cls}"><div class="ai-block-title">{label}</div><div class="ai-block-text">{_html.escape(text)}</div></div>'

        return f'''<div class="ai-section" id="ai-summary">
<div class="ai-header">
  <span class="ai-title">🤖 AI Research Brief</span>
  <span class="ai-model-tag">{model}</span>
  <span class="ai-format-tag">Brief + IMRaD</span>
</div>
<div class="ai-brief">{_html.escape(brief)}</div>
<details class="ai-detail-toggle">
  <summary>📋 Expand detailed report (IMRaD)</summary>
  <div class="ai-detail-content">{detail_items}</div>
</details>
</div>'''
    except Exception as e:
        print(f"  [WARN] AI Summary failed: {e}")
        return f'<div class="ai-error">⚠ AI summary failed: {e}</div>'


def main():
    p = argparse.ArgumentParser(description=f"MMPV-RNA {VERSION} — 独立报告生成器")
    p.add_argument("-o", "--output-dir", required=True, help="流水线输出根目录 (包含 00a_CleanData/ ... 10_Reports/)")
    p.add_argument("--skip-sankey", action="store_true", help="跳过 Sankey 图生成")
    p.add_argument("--skip-html", action="store_true", help="仅生成 TSV, 不生成 HTML")
    p.add_argument("--blast-db", help="BLAST 参考数据库路径 (用于序列相似度分类)")
    p.add_argument("--ai-summary", action="store_true", help="生成 AI 管线总结 (需 --ai-key)")
    p.add_argument("--ai-provider", default="openai", choices=["openai","ollama","deepseek","moonshot","custom"], help="AI 提供商 (default: openai)")
    p.add_argument("--ai-model", default="gpt-4o-mini", help="模型名 (default: gpt-4o-mini)")
    p.add_argument("--ai-key", default="", help="API Key (或 ollama 时留空)")
    p.add_argument("--ai-base-url", default="", help="自定义 API 地址 (如 http://localhost:11434/v1/chat/completions)")
    args = p.parse_args()
    # 进程守护: 停止主命令时连带杀掉所有子孙进程, 避免后台残留
    try:
        import process_guard
        process_guard.install()
    except Exception:
        pass

    root = Path(args.output_dir).resolve()
    if not root.is_dir():
        sys.exit(f"ERROR: 目录不存在: {root}")

    report_dir = root / _D['d_reports']
    report_dir.mkdir(parents=True, exist_ok=True)

    print(f"{'='*60}")
    print(f"Report Pipeline {VERSION}")
    print(f"  Output: {root}")
    print(f"  Reports: {report_dir}")
    print(f"{'='*60}")

    # 1. 收集数据 + 生成 TSV
    print("\n[1/4] Collecting stage data...")
    t0 = time.time()
    stage_stats = collect_data(root, report_dir, args.blast_db)

    # 2. 写入 stage_summary.tsv + 目录树
    with open(report_dir / "stage_summary.tsv", "w", newline="", encoding='utf-8') as sf:
        w = csv.DictWriter(sf, fieldnames=["Stage","Status","Key_Metric","Details"], delimiter="\t")
        w.writeheader()
        for s in stage_stats: w.writerow(s)

    tree_file = report_dir / "directory_tree.txt"
    with open(tree_file, "w", encoding='utf-8') as tf:
        for d in sorted(root.iterdir()):
            if not d.is_dir(): continue
            tf.write(f"{d.name}/\n")
            for sd in sorted(d.iterdir()):
                if sd.is_dir():
                    tf.write(f"  {sd.name}/\n")
                    for f in sorted(sd.iterdir())[:5]:
                        tf.write(f"    {f.name}\n")
                    rest = sum(1 for _ in sd.iterdir()) - 5
                    if rest > 0: tf.write(f"    ... +{rest} more\n")
    print(f"  stage_summary.tsv ({len(stage_stats)} rows), directory_tree.txt — {time.time()-t0:.0f}s")

    # 3. Sankey
    if not args.skip_sankey:
        print("\n[2/4] Generating Sankey diagrams...")
        t0 = time.time()
        generate_sankey(root, report_dir)
        print(f"  Done — {time.time()-t0:.0f}s")
    else:
        print("\n[2/4] Sankey: skipped")

    # 4. HTML
    if not args.skip_html:
        # 配置 AI 总结参数 (通过函数属性传递, 避免改 write_html_report 签名)
        if args.ai_summary:
            generate_ai_summary._provider = args.ai_provider
            generate_ai_summary._model = args.ai_model
            generate_ai_summary._api_key = args.ai_key
            # 自动设置 base_url
            if args.ai_base_url:
                generate_ai_summary._base_url = args.ai_base_url
            elif args.ai_provider == "ollama":
                generate_ai_summary._base_url = "http://localhost:11434/v1/chat/completions"
            elif args.ai_provider == "deepseek":
                generate_ai_summary._base_url = "https://api.deepseek.com/v1/chat/completions"
            elif args.ai_provider == "moonshot":
                generate_ai_summary._base_url = "https://api.moonshot.cn/v1/chat/completions"
            else:
                generate_ai_summary._base_url = ""
        else:
            generate_ai_summary._api_key = ""  # 禁用
        print("\n[3/4] Generating HTML report...")
        t0 = time.time()
        write_html_report(report_dir, stage_stats)
        print(f"  pipeline_report.html — {time.time()-t0:.0f}s")
        # AI 解读 section (可选, 纯增量: 离线/未配 key 时自动跳过, 不改动报告;
        # 与本文件既有 --ai-summary 旁路总结并存, 互不影响)
        try:
            from mmpv_common.ai import append_ai_section_to_html
            from mmpv_common.run_digest import guess_state_dir
            append_ai_section_to_html(report_dir / "pipeline_report.html",
                                      {"pipeline": "discovery", "stage_stats": stage_stats},
                                      state_dir=guess_state_dir(root))
        except Exception:
            pass
        # 过程面板 (确定性, 无需 AI): 项目 90_Handoff/state 在时嵌入全过程总览
        try:
            from mmpv_common.process_report import append_process_panel_to_html
            append_process_panel_to_html(report_dir / "pipeline_report.html", root)
        except Exception:
            pass
    else:
        print("\n[3/4] HTML: skipped")

    # Summary
    print(f"\n[4/4] {'='*50}")
    print(f"  Report complete!")
    print(f"    {report_dir}/pipeline_report.html")
    print(f"    {report_dir}/stage_summary.tsv")
    for tsv in ["data_summary","assembly_summary","ident_summary","filter_summary",
                "cobra_summary","hostdep_summary","checkv_summary","checkv_confidence"]:
        p = report_dir / f"{tsv}.tsv"
        if p.is_file(): print(f"    {report_dir}/{tsv}.tsv")
    print(f"{'='*50}")


if __name__ == "__main__":
    main()
