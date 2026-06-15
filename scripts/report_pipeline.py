#!/usr/bin/env python3
"""
report_pipeline.py — 独立报告生成器 (MMPV-RNA v2.3)

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

SCRIPT_DIR = Path(__file__).resolve().parent

# ═══════════════════════════════════════════════════════════════
# 辅助函数
# ═══════════════════════════════════════════════════════════════

def _count_fasta(path):
    if not path or not os.path.isfile(str(path)): return 0
    with open(str(path)) as f:
        return sum(1 for _ in f if _.startswith('>'))

def _count_lines(path):
    if not path or not os.path.isfile(str(path)): return 0
    with open(str(path)) as f:
        return sum(1 for _ in f)

def _count_dir(d):
    return sum(1 for _ in Path(d).iterdir() if _.is_dir()) if d and Path(d).is_dir() else 0

def _read_tsv(path):
    rows = []
    p = Path(path)
    if not p.is_file(): return rows
    with open(p) as f:
        hdr = f.readline().strip().split('\t')
        for line in f:
            if not line.strip(): continue
            rows.append(dict(zip(hdr, line.strip().split('\t'))))
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
    with open(str(fa_path)) as f:
        for line in f:
            s = line.strip()
            if s.startswith('>'):
                if seq: lens.append(len(seq))
                seq = ""
            else: seq += s
    if seq: lens.append(len(seq))
    return lens


# ═══════════════════════════════════════════════════════════════
# 阶段数据收集 (各阶段独立函数)
# ═══════════════════════════════════════════════════════════════

def _collect_cleandata(root, report_dir, _add):
    """00a_CleanData 阶段"""
    clean = root / "00a_CleanData"
    if not clean.is_dir():
        _add("00a_CleanData", "○", details="未运行"); return
    ns = _count_dir(clean / "3.clumpify") or _count_dir(clean / "2.fasta") or _count_dir(clean)
    _add("00a_CleanData", "✓", key_metric=f"{ns} 样本", details=str(clean))
    fp = clean / "logs"
    if not fp.is_dir(): return
    jf_list = list(fp.glob("*_fastp_report.json"))
    if not jf_list: return
    n_total_before, n_total_after = 0, 0
    b_total_before, b_total_after = 0, 0  # bases
    with open(report_dir / "data_summary.tsv", "w") as ds:
        ds.write("Sample\tRaw_Reads\tClean_Reads\tRetained(%)\tRaw_Bases\tClean_Bases\tRaw_Q20(%)\tClean_Q20(%)\tRaw_Q30(%)\tClean_Q30(%)\tLowQ_Reads\tTooShort_Reads\tDup_Rate(%)\n")
        for jf in sorted(jf_list):
            try:
                with open(jf) as jfh:
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
    hostdep = root / "00b_HostDepletion"
    if not hostdep.is_dir():
        _add("00b_HostDepletion", "○", details="未运行"); return
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
        with open(report_dir / "hostdep_summary.tsv", "w") as hf:
            cols = ["Sample","Raw","After_Kraken2","After_Host","Total_rRNA","non_rRNA","rRNA"]
            hf.write("\t".join(cols)+"\n")
            for sn in sorted(hd_rows):
                r = hd_rows[sn]
                hf.write("\t".join(str(r.get(c,0)) for c in cols)+"\n")
    _add("00b_HostDepletion", "✓", key_metric=f"{ns} 样本", details=str(hostdep))


def _collect_assembly(root, report_dir, _add):
    """01_Assembly 阶段"""
    asm = root / "01_Assembly"
    if not asm.is_dir():
        _add("01_Assembly", "○", details="未运行"); return
    ns = _count_dir(asm)
    total_contigs, total_bp = 0, 0
    asm_data = []
    with open(report_dir / "assembly_summary.tsv", "w") as af:
        af.write("Sample\tAssembler\tSize(Mb)\tContigs\tMax_Len\tN50\tN90\t>500bp\t>500bp(%)\t>1000bp\t>1000bp(%)\n")
        for d in sorted(asm.iterdir()):
            if not d.is_dir(): continue
            sample_contigs = 0
            for f in d.glob("*.contig.fasta"):
                lens = _parse_fasta_lens(f)
                if not lens: continue
                n = len(lens); total = sum(lens); mx = max(lens)
                n50, n90 = _n50_n90(lens)
                c500 = sum(1 for l in lens if l > 500); r500 = round(c500/n*100,1) if n else 0
                c1000 = sum(1 for l in lens if l > 1000); r1000 = round(c1000/n*100,1) if n else 0
                at = f.stem.replace(f"{d.name}_", "").replace(".contig", "")
                af.write(f"{d.name}\t{at}\t{total/1e6:.1f}\t{n}\t{mx}\t{n50}\t{n90}\t{c500}\t{r500}\t{c1000}\t{r1000}\n")
                asm_data.append({'s': d.name, 'n': n, 'total': total, 'n50': n50})
                total_contigs += n; total_bp += total
                sample_contigs += n
            _add(f"  └ {d.name}", "✓", key_metric=f"{sample_contigs} contigs")
        if len(asm_data) > 1:
            t_n = sum(r['n'] for r in asm_data); t_bp = sum(r['total'] for r in asm_data)
            af.write(f"TOTAL\tall\t{t_bp/1e6:.1f}\t{t_n}\t-\t-\t-\t-\t-\t-\t-\n")
    _add("01_Assembly", "✓", key_metric=f"{ns} 样本, {total_contigs:,} contigs, {total_bp/1e6:.1f} Mb", details=str(asm))
    as_script = SCRIPT_DIR.parent / "analysis" / "assembly_stats.py"
    if as_script.is_file():
        as_out = report_dir / "assembly_detail"
        as_out.mkdir(parents=True, exist_ok=True)
        try: subprocess.run([sys.executable, str(as_script), "-a", str(asm), "-o", str(as_out / "assembly_summary.tsv")], capture_output=True, timeout=60)
        except Exception as e: print(f"  [WARN] assembly_stats.py 失败: {e}")


def _collect_identification(root, report_dir, _add):
    """02_Identification 阶段"""
    ident = root / "02_Identification"
    if not ident.is_dir():
        _add("02_Identification", "○", details="未运行"); return
    ns = _count_dir(ident)
    n_virus = 0
    for d in ident.iterdir():
        if not d.is_dir(): continue
        for f in d.glob("*virus.all.candidate.fasta"): n_virus += _count_fasta(f)
    _add("02_Identification", "✓", key_metric=f"{ns} 样本, {n_virus:,} 病毒序列", details=str(ident))
    tools_list = ['genomad','blast','metabuli','virsorter2','viralverify','virhunter','virbot','viralm','rdrpcatch']
    ident_data = []
    with open(report_dir / "ident_summary.tsv", "w") as ids:
        ids.write("Sample\tAll_Candidate\t" + "\t".join(tools_list) + "\n")
        for d in sorted(ident.iterdir()):
            if not d.is_dir(): continue
            all_ids = _count_fasta(d / f"{d.name}_virus.all.candidate.fasta")
            tcounts = {}
            for tool in tools_list:
                idf = d / f"{d.name}_virus.{tool}.result.id"
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
    filter_data = []; modes_seen = set()
    with open(report_dir / "filter_summary.tsv", "w") as fs:
        fs.write("Sample\tMode\tAll_Candidate\tPassed\tRetained(%)\n")
        for d in sorted(ident.iterdir()):
            if not d.is_dir(): continue
            all_n = _count_fasta(d / f"{d.name}_virus.all.candidate.fasta")
            for fm, fd in [('filter','uniprot_filter_output_filter'),('strict','uniprot_filter_output_strict'),('comb','uniprot_filter_output')]:
                ff = d / fd / f"{d.name}_virus.uniprot_filtered.fasta"
                nf = _count_fasta(ff) if ff.is_file() else 0
                if nf > 0 or (nf == 0 and fd == 'uniprot_filter_output'):
                    fs.write(f"{d.name}\t{fm}\t{all_n}\t{nf}\t{round(nf/max(all_n,1)*100,1)}\n")
                    filter_data.append({'Sample': d.name, 'Mode': fm, 'All': all_n, 'Passed': nf})
                    modes_seen.add(fm)
        if filter_data:
            for m in sorted(modes_seen):
                mr = [r for r in filter_data if r['Mode'] == m]
                t_all = sum(r['All'] for r in mr); t_pass = sum(r['Passed'] for r in mr)
                fs.write(f"TOTAL\t{m}\t{t_all}\t{t_pass}\t{round(t_pass/max(t_all,1)*100,1)}\n")
    # Per-tool raw/filter/strict 统计
    with open(report_dir / "tool_filter_summary.tsv", "w") as tfs:
        tfs.write("Tool\tRaw\tFilter\tStrict\n")
        tool_filter = {t: [0, 0, 0] for t in tools_list}  # [raw, filter, strict]
        for d in sorted(ident.iterdir()):
            if not d.is_dir(): continue
            for tool in tools_list:
                tool_filter[tool][0] += _count_lines(d / f"{d.name}_virus.{tool}.result.id")
                tool_filter[tool][1] += _count_lines(d / "uniprot_filter_output_filter" / f"{d.name}_virus.{tool}.uniprot_filtered.id")
                tool_filter[tool][2] += _count_lines(d / "uniprot_filter_output_strict" / f"{d.name}_virus.{tool}.uniprot_filtered.id")
        for tool in sorted(tool_filter, key=lambda t: -tool_filter[t][0]):
            raw, filt, strict = tool_filter[tool]
            if raw > 0:
                tfs.write(f"{tool}\t{raw}\t{filt}\t{strict}\n")
    is_script = SCRIPT_DIR.parent / "analysis" / "ident_stats.py"
    if is_script.is_file():
        is_out = report_dir / "ident_detail"
        is_out.mkdir(parents=True, exist_ok=True)
        try: subprocess.run([sys.executable, str(is_script), "-i", str(ident), "-o", str(is_out)], capture_output=True, timeout=120)
        except Exception as e: print(f"  [WARN] ident_stats.py 失败: {e}")


def _collect_cobra(root, report_dir, _add):
    """03_COBRA 阶段"""
    cobra = root / "03_COBRA"
    if not cobra.is_dir():
        _add("03_COBRA", "○", details="未运行"); return
    ns = _count_dir(cobra)
    n_ext, n_queries, n_orphan, n_total_gain = 0, 0, 0, 0
    with open(report_dir / "cobra_summary.tsv", "w") as cs:
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
                    with open(str(cobra_logs[0])) as lf:
                        for line in lf:
                            s = line.strip()
                            if s.startswith('# Total queries:'): tq = int(s.split(':')[1].strip())
                            elif 'Extended_circular' in s: ec = int(s.split(':')[1].strip().split()[0])
                            elif 'Extended_partial' in s: ep = int(s.split(':')[1].strip().split()[0])
                            elif 'Extended_failed' in s: ef = int(s.split(':')[1].strip())
                            elif 'Orphan end' in s: oe = int(s.split(':')[1].strip())
                except Exception as e: print(f"  [WARN] COBRA log 解析失败 ({sd.name}): {e}")
            n_queries += tq; n_orphan += oe
            er = round((ec+ep)/max(tq,1)*100, 1); or_ = round(oe/max(tq,1)*100, 1)
            cs.write(f"{sn}\t{tq}\t{ec}\t{ep}\t{ef}\t{oe}\t{er}\t{or_}\t{sn_ext}\t{sn_gain}\n")
    _add("03_COBRA", "✓", key_metric=f"{ns} 样本, {n_ext} 延伸, {n_queries} query, {n_orphan} orphan, {n_total_gain:,} bp gain")
    cs_script = SCRIPT_DIR.parent / "analysis" / "cobra_stats.py"
    if cs_script.is_file():
        cs_out = report_dir / "cobra_detail"
        cs_out.mkdir(parents=True, exist_ok=True)
        try: subprocess.run([sys.executable, str(cs_script), "-c", str(cobra), "-o", str(cs_out)], capture_output=True, timeout=120)
        except Exception as e: print(f"  [WARN] cobra_stats.py 失败: {e}")


def _collect_cluster(root, report_dir, _add):
    """04_CLUSTER 阶段"""
    cluster = root / "04_CLUSTER"
    if not cluster.is_dir():
        _add("04_CLUSTER", "○", details="未运行"); return
    centroids_fa = cluster / "centroids" / "final_centroids.fasta"
    known_fa = cluster / "2_cdhit" / "known_centroids.fasta"
    n_centroids = _count_fasta(centroids_fa) if centroids_fa.is_file() else 0
    n_known = _count_fasta(known_fa) if known_fa.is_file() else 0
    _add("04_CLUSTER", "✓", key_metric=f"{n_centroids:,} novel + {n_known:,} known centroids", details=str(cluster))
    ctsv = cluster / "3_vclust" / "vclust_clusters.tsv"
    if ctsv.is_file():
        n_clusters = _count_lines(ctsv) - 1
        sizes = []
        with open(ctsv) as cf:
            cf.readline()
            for line in cf:
                members = line.strip().split('\t')[1] if '\t' in line else ""
                sizes.append(len(members.split(',')) if members else 0)
        singletons = sum(1 for sz in sizes if sz <= 1)
        max_sz = max(sizes) if sizes else 0
        _add("  └ vclust", "✓", key_metric=f"{n_clusters:,} 簇, {singletons} 单例, 最大簇={max_sz}")
    known_linked_fa = cluster / "2_cdhit" / "known_linked_centroids.fasta"
    n_linked = _count_fasta(known_linked_fa) if known_linked_fa.is_file() else 0
    if n_linked > 0: _add("  └ CD-HIT linked", "✓", key_metric=f"{n_linked} 有关联contig的已知簇")
    if n_known - n_linked > 0: _add("  └ CD-HIT pure", "○", key_metric=f"{n_known - n_linked} 纯参考簇 (不进下游)")


def _collect_taxonomy(root, report_dir, _add):
    """05_Taxonomy 阶段"""
    tax = root / "05_Taxonomy"
    final_tax = tax / "integrated" / "final_integrated_classification.tsv"
    if final_tax.is_file():
        n = _count_lines(final_tax) - 1
        counts = {"Known": 0, "Novel_Species": 0, "Novel_Genus": 0, "Novel_Family": 0}
        rank_fill = {"Realm":0,"Kingdom":0,"Phylum":0,"Class":0,"Order":0,"Family":0,"Genus":0,"Species":0}
        with open(final_tax) as tf:
            for row in csv.DictReader(tf, delimiter="\t"):
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
        _add("05_Taxonomy", "✓", key_metric=f"{n} 条, ★{counts['Known']} 已知, ★★{counts['Novel_Species']} 新种", details=str(tax))
        _add("  └ Novel Rank", "✓", key_metric=f"Known={counts['Known']} NewSp={counts['Novel_Species']} NewGe={counts['Novel_Genus']} NewFa={counts['Novel_Family']}")
        rk_info = " ".join(f"{r}={rank_fill[r]}" for r in ["Realm","Phylum","Class","Order","Family","Genus","Species"] if rank_fill.get(r,0) > 0)
        if rk_info: _add("  └ Rank Fill", "✓", key_metric=rk_info[:120])
    elif tax.is_dir():
        _add("05_Taxonomy", "○", key_metric="已运行但无最终结果", details=str(tax))
    else:
        _add("05_Taxonomy", "○", details="未运行")

    # 最终植物病毒分类汇总 (all_plant_viruses.fasta × taxonomy)
    all_plant_fa = root / "08_Rescue" / "all_plant_viruses.fasta"
    if all_plant_fa.is_file() and final_tax.is_file():
        plant_ids = set()
        for line in open(all_plant_fa):
            if line.startswith('>'): plant_ids.add(line[1:].split()[0])
        counts_p = {"Known": 0, "Novel_Species": 0, "Novel_Genus": 0, "Novel_Family": 0}
        with open(final_tax) as tf:
            for row in csv.DictReader(tf, delimiter="\t"):
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
                 key_metric=f"{n_plant} 条, Known={counts_p['Known']} NewSp={counts_p['Novel_Species']} NewGe={counts_p['Novel_Genus']} NewFa={counts_p['Novel_Family']}")
            # 写入 plant_virus_taxonomy.tsv 供 HTML 图表使用
            with open(report_dir / "plant_virus_taxonomy.tsv", "w") as pvf:
                pvf.write("Category\tCount\n")
                for k, v in counts_p.items(): pvf.write(f"{k}\t{v}\n")


def _collect_hostprediction(root, report_dir, _add):
    """06_HostPrediction 阶段"""
    host = root / "06_HostPrediction"
    host_summary = host / "ensemble_host_summary.tsv"
    if host_summary.is_file():
        n = _count_lines(host_summary) - 1
        try:
            import polars as pl
            hdf = pl.read_csv(str(host_summary), separator="\t", null_values=["NA","N/A",""])
            if "Final_Host" in hdf.columns:
                hcounts = hdf.group_by("Final_Host").agg(pl.len()).sort("len", descending=True)
                all_hosts = " | ".join(f"{r[0]}={r[1]}" for r in hcounts.iter_rows())
                _add("06_HostPrediction", "✓", key_metric=f"{n} 条, {all_hosts}", details=str(host))
        except Exception as e:
            print(f"  [WARN] HostPrediction polars 解析失败: {e}")
            _add("06_HostPrediction", "✓", key_metric=f"{n} 条", details=str(host))
    elif host.is_dir():
        _add("06_HostPrediction", "○", key_metric="已运行但无最终结果", details=str(host))
    else:
        _add("06_HostPrediction", "○", details="未运行")


def _collect_checkv(root, report_dir, _add):
    """07_CheckV 阶段"""
    cv_dir = root / "07_Checkv"
    QUALITY_ORDER = ["Complete","High-quality","Medium-quality","Low-quality","Not-determined"]
    if not cv_dir.is_dir():
        _add("07_CheckV", "○", details="未运行"); return
    cv_tsvs = list(cv_dir.rglob("completeness.tsv"))
    if not cv_tsvs:
        _add("07_CheckV", "✓", key_metric="已运行", details=str(cv_dir)); return
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
                    conf = str(row.get("aai_confidence","")).strip()
                    if not conf or conf in ("NA","N/A",""): conf = "Not-determined"
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
    with open(report_dir / "checkv_summary.tsv", "w") as cvf:
        cvf.write("Host\t" + "\t".join(QUALITY_ORDER) + "\tTotal\tHQ\n")
        for h in sorted(host_qd):
            hqd = host_qd[h]; t = sum(hqd.values()); hq = hqd["Complete"]+hqd["High-quality"]
            cvf.write(f"{h}\t"+"\t".join(str(hqd[q]) for q in QUALITY_ORDER)+f"\t{t}\t{hq}\n")
        cvf.write(f"TOTAL\t"+"\t".join(str(global_qd[q]) for q in QUALITY_ORDER)+f"\t{n_total}\t{n_hq}\n")
    if host_conf:
        all_confs = set()
        for hc in host_conf.values(): all_confs.update(hc.keys())
        conf_order = sorted(all_confs)
        with open(report_dir / "checkv_confidence.tsv", "w") as cff:
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
    """08_Rescue 阶段"""
    rescue = root / "08_Rescue"
    if not rescue.is_dir():
        _add("08_Rescue", "○", details="未运行"); return
    rescue_finals = list(rescue.rglob("final_centroids.fasta"))
    n_rescued = 0
    for rf in rescue_finals:
        if "branch" not in str(rf) and "known" not in str(rf): n_rescued += _count_fasta(rf)
    branch_info = []
    for bname, blabel in [("branch_a","CheckV"),("branch_b","VSI"),("branch_c","BLASTN")]:
        for bd in rescue.rglob(bname):
            if not bd.is_dir(): continue
            pass_fa = bd / f"{'branchA' if bname=='branch_a' else 'branchB' if bname=='branch_b' else 'branchC'}_pass.fasta"
            if not pass_fa.is_file():
                pass_fa = bd / f"{'branchB' if bname=='branch_b' else 'branchC'}_pass.fasta"
            if pass_fa.is_file():
                bp = _count_fasta(pass_fa)
                if bp > 0: branch_info.append(f"{blabel}={bp}")
    _add("08_Rescue", "✓", key_metric=f"{n_rescued:,} HQ vOTU ({' | '.join(branch_info) if branch_info else '0 pas'})", details=str(rescue))


def collect_data(output_dir, report_dir, blast_db=None):
    """收集所有阶段数据 → 生成 TSV + 返回 stage_stats"""
    root = Path(output_dir).resolve()
    stage_stats = []

    def _add(stage, status, key_metric="", details=""):
        stage_stats.append({"Stage": stage, "Status": status, "Key_Metric": key_metric, "Details": details})

    _collect_cleandata(root, report_dir, _add)
    _collect_hostdepletion(root, report_dir, _add)
    _collect_assembly(root, report_dir, _add)
    _collect_identification(root, report_dir, _add)
    _collect_cobra(root, report_dir, _add)
    _collect_cluster(root, report_dir, _add)
    _collect_taxonomy(root, report_dir, _add)
    _collect_hostprediction(root, report_dir, _add)
    _collect_checkv(root, report_dir, _add)
    _collect_rescue(root, report_dir, _add)

    # 最终植物病毒汇总表 + 旭日图
    _generate_plant_virus_summary(root, report_dir, _add, blast_db)

    return stage_stats


def _generate_plant_virus_summary(root, report_dir, _add, blast_db=None):
    """生成 plant_virus_summary.tsv + taxonomy_sunburst.html"""
    all_plant_fa = root / "08_Rescue" / "all_plant_viruses.fasta"
    if not all_plant_fa.is_file(): return

    # 1. 读取 contig IDs + lengths + 来源
    plant_data = {}  # {contig_id: {length, source, ...}}
    for src_label, src_fa in [("免拯救", root / "08_Rescue" / "known" / "centroids" / "final_centroids.fasta"),
                               ("rescued", root / "08_Rescue" / "Plant" / "centroids" / "final_centroids.fasta")]:
        if not src_fa.is_file(): continue
        for rec in SeqIO.parse(str(src_fa), "fasta"):
            plant_data[rec.id] = {"contig_id": rec.id, "length": len(rec.seq), "source": src_label}

    if not plant_data: return

    # 2. CheckV 数据: 优先 post-rescue (最新评估), 回退 CheckV 阶段原始值
    cv_data = {}
    for cv_tsv in [root / "08_Rescue" / "checkv" / "Plant" / "completeness.tsv",
                   root / "08_Rescue" / "checkv" / "no_rescue" / "completeness.tsv",
                   root / "07_Checkv" / "Plant" / "completeness.tsv"]:
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

    # 3. Taxonomy
    tax_tsv = root / "05_Taxonomy" / "integrated" / "final_integrated_classification.tsv"
    tax_data = {}
    if tax_tsv.is_file():
        with open(tax_tsv) as tf:
            for row in csv.DictReader(tf, delimiter="\t"):
                cid = row.get("contig_id","")
                if cid in plant_data:
                    tax_data[cid] = {rk: row.get(rk, row.get(rk.lower(),"")) for rk in
                                     ["Realm","Kingdom","Phylum","Class","Order","Family","Genus","Species"]}

    # 4. BLAST 相似度分类 (参考 VIGA 阈值: Known≥95% | NewSp≥78% | NewGe≥65% | NewFa<65%)
    plant_records = list(SeqIO.parse(str(all_plant_fa), "fasta"))
    novel_by_sim = {}
    _bd = blast_db
    if _bd:
        _bd = Path(_bd)
        # BLAST DB 判定: 父目录中有 {db_basename}.nin 文件
        db_dir = _bd.parent if _bd.parent.is_dir() else None
        db_basenames = set(os.listdir(str(db_dir))) if db_dir else set()
        bname = _bd.name
        if not any(bname+ext in db_basenames for ext in [".nin",".nal"]):
            print(f"  [WARN] BLAST 数据库无效: {_bd} (缺少 {bname}.nin)"); _bd = None
    if not _bd:
        for p in [Path("/home/zhangwenda/db/viral_nt"), Path("/home/zhangwenda/db/nt"),
                  root.parent / "database" / "viral_nt",
                  Path("/home/zhangwenda/database/virus-db/ncbi-virus_ref/ncbi-virus_ref.blast.db")]:
            dp = p.parent; bn = p.name
            bfs = set(os.listdir(str(dp))) if dp.is_dir() else set()
            if any(bn+ext in bfs for ext in [".nin",".nal"]):
                _bd = p; break
    if not _bd and os.environ.get("BLAST_DB"):
        _bd = Path(os.environ["BLAST_DB"])
    blast_db = _bd

    if blast_db:
        print(f"  BLAST 参考数据库: {blast_db}")
        cq_fa = report_dir / "tmp_plant_virus.fa"
        SeqIO.write(plant_records, str(cq_fa), "fasta")
        blast_out = report_dir / "tmp_blast.tsv"
        subprocess.run(["blastn", "-task", "dc-megablast", "-query", str(cq_fa),
                      "-db", str(blast_db), "-outfmt", "6 qseqid sseqid pident length",
                      "-max_target_seqs", "1", "-evalue", "1e-5", "-num_threads", "4",
                      "-out", str(blast_out)], capture_output=True, check=False)
        if blast_out.is_file():
            n_hits = 0
            for line in open(blast_out):
                parts = line.strip().split('\t')
                if len(parts) >= 3:
                    cid, pident = parts[0], float(parts[2])
                    for label, thresh in [("Known",95),("NewSp",78),("NewGe",65)]:
                        if pident >= thresh:
                            novel_by_sim[cid] = label; break
                    else: novel_by_sim[cid] = "NewFa"
                    n_hits += 1
            print(f"  BLAST 命中: {n_hits}/{len(plant_records)} 条")
        # 无命中用 plant_records 补齐
        for rec in plant_records:
            if rec.id not in novel_by_sim:
                novel_by_sim[rec.id] = "NewFa"
        for _tf in [cq_fa, blast_out]:
            try: _tf.unlink()
            except: pass
    else:
        print("  [SKIP] BLAST 相似度分类: 未找到参考数据库 (设置 --blast-db 或 BLAST_DB)")

    # 写入 plant_virus_summary.tsv
    cols = ["contig_id","length","source","aai_completeness","aai_confidence",
            "viral_length","aai_expected_length","kmer_freq",
            "Realm","Kingdom","Phylum","Class","Order","Family","Genus","Species",
            "novelty_tool","novelty_similarity"]
    with open(report_dir / "plant_virus_summary.tsv", "w") as pvf:
        pvf.write("\t".join(cols) + "\n")
        for cid in sorted(plant_data):
            d = plant_data[cid]; cv = cv_data.get(cid, {}); tx = tax_data.get(cid, {})
            # tool-based novelty
            sp, ge, fa = tx.get("Species",""), tx.get("Genus",""), tx.get("Family","")
            nt = "Known" if (sp and sp not in ("NA","-")) else "NewSp" if (ge and ge not in ("NA","-")) else "NewGe" if (fa and fa not in ("NA","-")) else "NewFa"
            ns = novel_by_sim.get(cid, "NA")
            vals = [cid, d["length"], d["source"],
                    cv.get("aai_completeness","NA"), cv.get("aai_confidence","NA"),
                    cv.get("viral_length","NA"), cv.get("aai_expected_length","NA"), cv.get("kmer_freq","NA")]
            vals += [tx.get(rk,"") for rk in ["Realm","Kingdom","Phylum","Class","Order","Family","Genus","Species"]]
            vals += [nt, ns]
            pvf.write("\t".join(str(v) for v in vals) + "\n")

    # 按 BLAST 相似度分类更新环形图数据
    if novel_by_sim:
        p_counts = {"Known":0,"Novel_Species":0,"Novel_Genus":0,"Novel_Family":0}
        for cid, label in novel_by_sim.items():
            if label in ("Known",): p_counts["Known"] += 1
            elif label in ("NewSp",): p_counts["Novel_Species"] += 1
            elif label in ("NewGe",): p_counts["Novel_Genus"] += 1
            else: p_counts["Novel_Family"] += 1
        with open(report_dir / "plant_virus_taxonomy.tsv", "w") as pvf:
            pvf.write("Category\tCount\n")
            for k, v in p_counts.items(): pvf.write(f"{k}\t{v}\n")

    n = len(plant_data)
    n_no_rescue = sum(1 for v in plant_data.values() if v["source"]=="免拯救")
    n_rescued = n - n_no_rescue
    _add("09_Plant_virus", "✓",
         key_metric=f"{n} 条 (免拯救={n_no_rescue} + rescued={n_rescued})",
         details=str(report_dir / "plant_virus_summary.tsv"))

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
                          include_plotlyjs='cdn', full_html=True)
            print(f"  Sunburst → {report_dir / 'taxonomy_sunburst.html'}")
    except Exception as e:
        print(f"  [WARN] 旭日图生成失败: {e}")


# ═══════════════════════════════════════════════════════════════
# Sankey 图生成
# ═══════════════════════════════════════════════════════════════

def generate_sankey(output_dir, report_dir):
    """生成交互式 taxonomy Sankey HTML (全部 + 植物病毒)"""
    tax = Path(output_dir) / "05_Taxonomy"
    final_tax = tax / "integrated" / "final_integrated_classification.tsv"
    if not final_tax.is_file(): return
    sankey_script = SCRIPT_DIR.parent / "analysis" / "taxonomic_sankey.py"
    if not sankey_script.is_file(): return
    try:
        import importlib.util
        if not importlib.util.find_spec("plotly"):
            print("  安装 plotly...")
            subprocess.run([sys.executable, "-m", "pip", "install", "plotly", "-q"], check=False)
        # 生成交互式 HTML (--format html, 不设 title 避免重叠)
        subprocess.run([sys.executable, str(sankey_script),
                        "-i", str(final_tax), "-o", str(report_dir / "classification_sankey.html"),
                        "--format", "html", "--min-flow", "1", "--min-genus-flow", "10",
                        "--palette", "set3", "--height", "860", "--width", "1000", "--node-pad", "30",
                        "--label-truncate", "25", "--font-size", "9", "--title-font-size", "14",
                        "--title", ""],
                       capture_output=True, timeout=120)
        print(f"  Sankey HTML → {report_dir / 'classification_sankey.html'}")
        # Plant-only
        host_summary = Path(output_dir) / "06_HostPrediction" / "ensemble_host_summary.tsv"
        if host_summary.is_file():
            try:
                import polars as pl
                hdf = pl.read_csv(str(host_summary), separator="\t", null_values=["NA","N/A",""])
                plant_ids = set(hdf.filter(pl.col("Final_Host")=="Plant")["contig_id"].to_list())
                if plant_ids:
                    plant_tax = report_dir / "plant_final_taxonomy.tsv"
                    with open(final_tax) as tf, open(plant_tax, "w") as pf:
                        pf.write(tf.readline())
                        for line in tf:
                            cid = line.split('\t')[0].strip('"')
                            if cid in plant_ids: pf.write(line)
                    subprocess.run([sys.executable, str(sankey_script),
                                    "-i", str(plant_tax), "-o", str(report_dir / "classification_sankey_plant.html"),
                                    "--format", "html", "--min-flow", "1", "--min-genus-flow", "5",
                                    "--palette", "set3", "--height", "860", "--width", "1000", "--node-pad", "30",
                                    "--label-truncate", "25", "--font-size", "9", "--title-font-size", "14",
                                    "--title", ""],
                                   capture_output=True, timeout=120)
                    print(f"  Plant Sankey HTML → {report_dir / 'classification_sankey_plant.html'}")
            except Exception as e: print(f"  [WARN] Plant Sankey 生成失败: {e}")
    except Exception as e: print(f"  [WARN] Sankey 生成失败: {e}")


# ═══════════════════════════════════════════════════════════════
# HTML 报告生成
# ═══════════════════════════════════════════════════════════════

def write_html_report(report_dir, stage_stats):
    """生成期刊级 HTML 流水线报告"""
    import json as _json

    # 内嵌 Chart.js (避免 CDN 不可用导致图表空白)
    chart_js_path = SCRIPT_DIR / "chart.min.js"
    chart_js_inline = ""
    if chart_js_path.is_file():
        with open(chart_js_path, "r", encoding="utf-8") as cf:
            chart_js_inline = cf.read()
        # 同时复制到报告目录供离线使用
        import shutil
        shutil.copy2(str(chart_js_path), str(report_dir / "chart.min.js"))

    main_stages = [s for s in stage_stats if not s["Stage"].startswith("  ")]

    S = {"✓": "pass", "○": "skip", "✗": "fail"}
    def _esc(v): return str(v).replace("&","&amp;").replace("<","&lt;").replace(">","&gt;")
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
                    {"label":"Raw Reads (M)","data":[round(v/1e6,1) for v in dq_raw],"backgroundColor":"#90a4ae","yAxisID":"y"},
                    {"label":"Clean Reads (M)","data":[round(v/1e6,1) for v in dq_clean],"backgroundColor":"#42a5f5","yAxisID":"y"},
                    {"label":"Q20 Raw (%)","data":[round(v,1) for v in dq_q20b],"backgroundColor":"#ffa726","yAxisID":"y1"},
                    {"label":"Q20 Clean (%)","data":[round(v,1) for v in dq_q20a],"backgroundColor":"#66bb6a","yAxisID":"y1"},
                ]}, {"responsive":True,"plugins":{"title":{"display":True,"text":"Read Quality & Filtering"}},
                     "scales":{"y":{"beginAtZero":True,"position":"left","title":{"text":"Reads (M)"}},
                              "y1":{"beginAtZero":True,"position":"right","max":100,"grid":{"drawOnChartArea":False},"title":{"text":"Q20 (%)"}}}})
    if dq_rows and any(v > 0 for v in dq_dup):
        chart_scripts += _chart('chart_s00a_dup', 'bar', {
            "labels": dq_samples,
            "datasets": [{"label":"Duplication Rate (%)","data":[round(v,2) for v in dq_dup],
                          "backgroundColor":['#ef5350' if v>30 else '#ffa726' if v>15 else '#66bb6a' for v in dq_dup]}]},
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
                "datasets": [{"label":"Retained (non-host)","data":hd_retained,"backgroundColor":"#66bb6a"},
                              {"label":"Removed (host+rRNA)","data":hd_removed,"backgroundColor":"#ef5350"}]},
                {"responsive":True,"indexAxis":"y","plugins":{"title":{"display":True,"text":"Host Depletion: Reads Retained vs Removed"}},
                 "scales":{"x":{"stacked":True,"beginAtZero":True,"title":{"text":"Reads"}},"y":{"stacked":True}}})

    # S01 — 组装
    asm_rows = _read_tsv(report_dir / "assembly_summary.tsv")
    if asm_rows:
        asm_samples = [r["Sample"][:12] for r in asm_rows if r.get("Sample","")!="TOTAL"]
        asm_n50 = [int(r.get("N50",0)) for r in asm_rows if r.get("Sample","")!="TOTAL"]
        asm_contigs = [int(r.get("Contigs",0)) for r in asm_rows if r.get("Sample","")!="TOTAL"]
        if asm_samples:
            stage_has_chart['s01a'] = True
            chart_scripts += _chart('chart_s01a', 'bar', {
                "labels": asm_samples,
                "datasets": [{"label":"N50 (kb)","data":[round(v/1000,1) for v in asm_n50],
                              "backgroundColor":['#1565c0' if v>5000 else '#42a5f5' if v>1000 else '#90caf9' for v in asm_n50]}]},
                {"responsive":True,"plugins":{"title":{"display":True,"text":"Assembly N50 (kb)"}},
                 "scales":{"y":{"beginAtZero":True,"title":{"text":"N50 (kb)"}}}})
            stage_has_chart['s01b'] = True
            chart_scripts += _chart('chart_s01b', 'bar', {
                "labels": asm_samples,
                "datasets": [{"label":"Contigs","data":asm_contigs,
                              "backgroundColor":['#66bb6a' if v>5000 else '#a5d6a7' if v>1000 else '#c8e6c9' for v in asm_contigs]}]},
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
                    {"label":"Raw","data":raw_vals,"backgroundColor":"#90caf9"},
                    {"label":"Filter","data":filt_vals,"backgroundColor":"#42a5f5"},
                    {"label":"Strict","data":strict_vals,"backgroundColor":"#1565c0"},
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
            ds_c = [{"label":"Extension Rate (%)","data":cobra_ext,"backgroundColor":"#66bb6a","yAxisID":"y"}]
            if cobra_orph:
                ds_c.append({"label":"Orphan Rate (%)","data":cobra_orph,"backgroundColor":"#ef5350","yAxisID":"y1"})
            scales_c = {"y":{"beginAtZero":True,"position":"left","title":{"text":"Rate (%)"}}}
            if cobra_orph:
                scales_c["y1"] = {"beginAtZero":True,"position":"right","grid":{"drawOnChartArea":False},"title":{"text":"Orphan (%)"}}
            chart_scripts += _chart('chart_s03', 'bar', {"labels":csamples,"datasets":ds_c},
                {"responsive":True,"plugins":{"title":{"display":True,"text":"COBRA Extension & Orphan Rate"}},"scales":scales_c})

    # S05 — Taxonomy novelty
    tax_novelty_kv = {}
    for s in stage_stats:
        if 'Novel Rank' in s['Stage'] and 'Known=' in s.get('Key_Metric',''):
            tax_novelty_kv = _extract_kv(s['Key_Metric'], r'(Known|NewSp|NewGe|NewFa)=(\d+)')
    if tax_novelty_kv:
        stage_has_chart['s05'] = True
        bar_colors = ["#2e7d32","#1565c0","#ef6c00","#c62828"]
        _tax_label_map = {"Known":"已知","NewSp":"新种","NewGe":"新属","NewFa":"新科"}
        chart_scripts += _chart('chart_s05a', 'bar', {
            "labels": [_tax_label_map.get(k,k) for k in tax_novelty_kv.keys()],
            "datasets": [{"label":"序列数","data":list(tax_novelty_kv.values()),"backgroundColor":bar_colors}]},
            {"responsive":True,"plugins":{"title":{"display":True,"text":"Taxonomy Novelty"},"legend":{"display":False}},
             "scales":{"y":{"beginAtZero":True,"title":{"text":"序列数"}}}})

    # S05b — 最终植物病毒分类 (all_plant_viruses.fasta × taxonomy)
    plant_tax_rows = _read_tsv(report_dir / "plant_virus_taxonomy.tsv")
    if plant_tax_rows:
        pt_kv = {r["Category"]: int(r["Count"]) for r in plant_tax_rows if int(r.get("Count",0)) > 0}
        if pt_kv:
            stage_has_chart['s05b'] = True
            pt_colors = ["#2e7d32","#1565c0","#ef6c00","#c62828"]
            _pt_label_map = {"Known":"已知","Novel_Species":"新种","Novel_Genus":"新属","Novel_Family":"新科"}
            chart_scripts += _chart('chart_s05b', 'doughnut', {
                "labels": [_pt_label_map.get(k,k) for k in pt_kv.keys()],
                "datasets": [{"data":list(pt_kv.values()),"backgroundColor":pt_colors[:len(pt_kv)]}]},
                {"responsive":True,"plugins":{"title":{"display":True,"text":"Plant Virus Taxonomy"},"legend":{"position":"bottom"}}})

    # S06 — Host distribution
    host_kv = {}
    for s in stage_stats:
        if '06_HostPrediction' in s['Stage'] and '=' in s.get('Key_Metric',''):
            host_kv = _extract_kv(s['Key_Metric'], r'(\w+)=(\d+)')
    if host_kv:
        stage_has_chart['s06'] = True
        host_colors = ["#42a5f5","#66bb6a","#ffa726","#ef5350","#ab47bc","#26c6da","#7e57c2","#78909c"]
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
            colors = ['#2e7d32','#1565c0','#ef6c00','#c62828','#9e9e9e']
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
            colors2 = ['#2e7d32','#ef6c00','#c62828','#1565c0','#9e9e9e']
            cf_datasets = [{"label":q,"data":cf_data[q],"backgroundColor":colors2[i%5]} for i,q in enumerate(cql)]
            chart_scripts += _chart('chart_s07b', 'bar', {"labels":cf_hosts,"datasets":cf_datasets},
                {"responsive":True,"plugins":{"title":{"display":True,"text":"CheckV aai_confidence by Host"}},
                 "scales":{"x":{"stacked":True},"y":{"stacked":True,"beginAtZero":True,"title":{"text":"Contig Count"}}}})

    # S08 — Rescue branches
    rescue_kv = {}
    for s in stage_stats:
        if '08_Rescue' in s['Stage']:
            rescue_kv = _extract_kv(s['Key_Metric'], r'(CheckV|VSI|BLASTN)=(\d+)')
    if rescue_kv:
        stage_has_chart['s08'] = True
        chart_scripts += _chart('chart_s08', 'pie', {
            "labels": list(rescue_kv.keys()),
            "datasets": [{"data":list(rescue_kv.values()),"backgroundColor":["#66bb6a","#42a5f5","#ffa726"]}]},
            {"responsive":True,"plugins":{"title":{"display":True,"text":"Rescue Branch Contributions"},"legend":{"position":"bottom"}}})

    # Sankey 交互式嵌入 (用 Blob URL 动态注入, 避免 data URI 大小限制)
    # s05=全部分类, s06=植物病毒
    sankey_by_stage = {}  # {stage_key: html_string}
    # 堆叠图默认切换为百分比模式 (延迟执行, 等图表全部渲染完)
    chart_scripts += "setTimeout(function(){document.querySelectorAll('.pct-btn').forEach(function(b){toggleStackedPct(b)})},100);\n"
    sankey_inject_scripts = ""
    sankey_map = [("classification_sankey.html","Taxonomy Classification Sankey","s05"),
                  ("classification_sankey_plant.html","Plant Virus Taxonomy Sankey","s06")]
    for i, (sname, stitle, stage_key) in enumerate(sankey_map):
        spath = report_dir / sname
        if spath.is_file():
            import base64
            with open(spath, "rb") as sf:
                sankey_b64 = base64.b64encode(sf.read()).decode()
            card = f'''<div class="sankey-card">
<h3>{stitle}</h3>
<iframe id="sankey_iframe_{i}" style="width:100%;height:700px;border:none;border-radius:4px" loading="lazy"></iframe>
</div>\n'''
            sankey_by_stage.setdefault(stage_key, "")
            sankey_by_stage[stage_key] += card
            sankey_inject_scripts += f"(function(){{var b='{sankey_b64}';var d=atob(b);var u=URL.createObjectURL(new Blob([d],{{type:'text/html'}}));document.getElementById('sankey_iframe_{i}').src=u;}})();\n"

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
        sankey_by_stage.setdefault("s09", "")
        sankey_by_stage["s09"] += sunburst_card
        sankey_inject_scripts += f"(function(){{var b='{sunburst_b64}';var d=atob(b);var u=URL.createObjectURL(new Blob([d],{{type:'text/html'}}));document.getElementById('sunburst_iframe').src=u;}})();\n"

    # ── KPI ──
    kpis = {}
    for s in stage_stats:
        if s['Stage'] in ('00a_CleanData', '  └ data_summary'):
            for m in re.finditer(r'reads:\s*([\d,]+)→([\d,]+)', s.get('Key_Metric','')):
                kpis['raw_reads'] = m.group(1); kpis['clean_reads'] = m.group(2)
            for m in re.finditer(r'bases:\s*([\d,]+)→([\d,]+)', s.get('Key_Metric','')):
                kpis['raw_bases'] = m.group(1); kpis['clean_bases'] = m.group(2)
        if s['Stage'] == '00a_CleanData':
            for m in re.finditer(r'(\d+)\s*样本', s.get('Key_Metric','')): kpis['n_sample'] = m.group(1)
        if s['Stage'] == '01_Assembly':
            for m in re.finditer(r'([\d,]+)\s*contigs', s.get('Key_Metric','')): kpis['total_contigs'] = m.group(1)
            for m in re.finditer(r'([\d.]+)\s*Mb', s.get('Key_Metric','')): kpis['total_mb'] = m.group(1)
        if s['Stage'] == '02_Identification':
            for m in re.finditer(r'([\d,]+)\s*病毒序列', s.get('Key_Metric','')): kpis['virus_seqs'] = m.group(1)
        if 'Novel Rank' in s['Stage']:
            kpis['novelty'] = s.get('Key_Metric','')
            for m in re.finditer(r'Known=(\d+)', s.get('Key_Metric','')): kpis['n_known'] = m.group(1)
            for m in re.finditer(r'NewSp=(\d+)', s.get('Key_Metric','')): kpis['n_newsp'] = m.group(1)
            for m in re.finditer(r'NewGe=(\d+)', s.get('Key_Metric','')): kpis['n_newge'] = m.group(1)
            for m in re.finditer(r'NewFa=(\d+)', s.get('Key_Metric','')): kpis['n_newfa'] = m.group(1)
        if '  └ vclust' in s['Stage']:
            for m in re.finditer(r'([\d,]+)\s*簇', s.get('Key_Metric','')): kpis['n_clusters'] = m.group(1)
        if '06_HostPrediction' in s['Stage'] and s['Stage'].startswith('06'):
            for m in re.finditer(r'([\d,]+)\s*条', s.get('Key_Metric','')): kpis['host_total'] = m.group(1)
        if '07_CheckV' in s['Stage'] and s['Stage'].startswith('07'):
            for m in re.finditer(r'(\d+)\s*HQ', s.get('Key_Metric','')): kpis['hq_votus'] = m.group(1)
            for m in re.finditer(r'(\d+)\s*total', s.get('Key_Metric','')): kpis['cv_total'] = m.group(1)
        if '08_Rescue' in s['Stage'] and s['Stage'].startswith('08'):
            for m in re.finditer(r'([\d,]+)\s*HQ vOTU', s.get('Key_Metric','')): kpis['rescued'] = m.group(1)
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
        ("s01","Assembly","01 Assembly","ASM","#1565c0"),
        ("s02","Ident","02 Identification","ID","#5c6bc0"),
        ("s03","COBRA","03 COBRA Extension","COBRA","#00897b"),
        ("s04","Cluster","04 Clustering","CLU","#ef6c00"),
        ("s05","Taxonomy","05 Taxonomy Classification","TAX","#6a1b9a"),
        ("s06","Host","06 Host Prediction","HOST","#c62828"),
        ("s07","CheckV","07 CheckV Quality","CV","#2e7d32"),
        ("s08","Rescue","08 Rescue","RESCUE","#37474f"),
        ("s09","PlantVirus","09 Plant Virus Collection","PV","#00838f"),
    ]

    _skey_to_num = {'s00a':'00a','s00b':'00b','s01':'01','s02':'02','s03':'03',
                    's04':'04','s05':'05','s06':'06','s07':'07','s08':'08','s09':'09'}
    stage_status = {}; stage_metric = {}
    for s in stage_stats:
        sn = s['Stage']
        for sk, _, _, _, _ in stage_defs:
            snum = _skey_to_num.get(sk,'')
            if sn.startswith(snum) or (snum+'_') in sn or sn == snum:
                stage_status[sk] = S.get(s['Status'],'skip')
                stage_metric[sk] = s.get('Key_Metric','')
                break

    chart_map = {
        's00a': [('chart_s00a','Read Quality'),('chart_s00a_dup','Duplication Rate')],
        's00b': [('chart_s00b','Host Depletion')],
        's01':  [('chart_s01b','Contig Count'),('chart_s01a','N50 (kb)')],
        's02':  [('chart_s02','UniProt Filter')],
        's03':  [('chart_s03','COBRA Rates')],
        's05':  [('chart_s05a','Taxonomy Novelty')],
        's09':  [('chart_s05b','Plant Virus Taxonomy')],
        's06':  [('chart_s06a','Host Distribution')],
        's07':  [('chart_s07a','CheckV Quality'),('chart_s07b','CheckV Confidence')],
        's08':  [('chart_s08','Rescue Branches')],
    }

    # 每个阶段对应的 TSV 文件, 用于在卡片内嵌入数据表
    stage_tsv_map = {
        's00a': ['data_summary.tsv'],
        's00b': ['hostdep_summary.tsv'],
        's01':  ['assembly_summary.tsv'],
        's02':  ['ident_summary.tsv', 'filter_summary.tsv'],
        's03':  ['cobra_summary.tsv'],
        's05':  [],
        's06':  [],
        's07':  ['checkv_summary.tsv', 'checkv_confidence.tsv'],
        's08':  [],
        's09':  ['plant_virus_summary.tsv'],
    }

    sections_html = ""
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
                elif sk == 's05':
                    if cid == 'chart_s05a' and stage_has_chart.get('s05'): active.append(cid)
                elif sk == 's09':
                    if cid == 'chart_s05b' and stage_has_chart.get('s05b'): active.append(cid)
                elif sk == 's07':
                    if cid == 'chart_s07a' and stage_has_chart.get('s07a'): active.append(cid)
                    if cid == 'chart_s07b' and stage_has_chart.get('s07b'): active.append(cid)
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
        for tsv_name in stage_tsv_map.get(sk, []):
            tsv_path = report_dir / tsv_name
            if not tsv_path.is_file(): continue
            tsv_rows = _read_tsv(tsv_path)
            if not tsv_rows: continue
            max_preview = 200 if tsv_name == "plant_virus_summary.tsv" else 8
            preview = tsv_rows[:max_preview]
            cols = list(preview[0].keys())
            th_h = "".join(f"<th>{_esc(c)}</th>" for c in cols)
            tr_h = ""
            for r in preview:
                tr_h += "<tr>" + "".join(f"<td>{_esc(str(r.get(c,'')))[:50]}</td>" for c in cols) + "</tr>"
            more = f' <span style="color:var(--muted);font-size:10px">(+{len(tsv_rows)-8} more)</span>' if len(tsv_rows) > 8 else ""
            dl_uri = _tsv_data_uri(tsv_path)
            dl_link = f'<a href="{dl_uri}" download="{tsv_name}" style="font-size:10px;margin-left:6px;color:var(--blue)">[download]</a>' if dl_uri else ""
            table_html += f'''<details class="stage-table-detail" open>
<summary>{tsv_name} — {len(tsv_rows)} rows{more} {dl_link}</summary>
<div style="overflow-x:auto;margin-top:6px"><table class="app-table"><thead><tr>{th_h}</tr></thead><tbody>{tr_h}</tbody></table></div>
</details>'''
        if table_html:
            table_html = f'<div class="stage-tables">{table_html}</div>'

        sections_html += f'''
<section class="stage {border_cls}" id="stage-{short}">
  <div class="stage-header">
    <div class="stage-icon" style="background:{color}">{icon}</div>
    <div class="stage-title"><h2>{full}</h2><span class="stage-metric">{_esc(metric)}</span></div>
    <span class="stage-badge {badge_cls}">{badge_txt}</span>
  </div>
  {table_html}
  {chart_html}
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
        ("Samples", f"{sample_kpi}<br>{n_sample} 样本", "🧬"),
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
    if r_raw: flow_stages.append(("Raw Reads",r_raw,"#0d1b3e"))
    if r_clean: flow_stages.append(("Clean Reads",r_clean,"#1565c0"))
    if r_hostfree: flow_stages.append(("Host-free Reads",r_hostfree,"#546e7a"))
    if r_contig: flow_stages.append(("Contigs",r_contig,"#00897b"))
    if r_virus: flow_stages.append(("Viral Seqs",r_virus,"#ef6c00"))
    if r_votu: flow_stages.append(("vOTU Clusters",r_votu,"#6a1b9a"))
    if r_hq: flow_stages.append(("HQ vOTUs",r_hq,"#2e7d32"))

    flow_html = ""
    if len(flow_stages) >= 3:
        f_labels = [s[0] for s in flow_stages]
        f_values = [s[1] for s in flow_stages]
        f_colors = [s[2] for s in flow_stages]
        flow_chart_id = "chart_pipeline_flow"
        chart_scripts += _chart(flow_chart_id, 'bar', {
            "labels": f_labels,
            "datasets": [{"label":"Count","data":f_values,"backgroundColor":f_colors,
                          "borderColor":f_colors,"borderWidth":0}]},
            {"indexAxis":"y","responsive":True,
             "plugins":{"title":{"display":True,"text":"Pipeline Flow — Reads → HQ vOTUs","font":{"size":14}},
                        "legend":{"display":False},
                        "tooltip":{"callbacks":{"label":"function(ctx){var v=ctx.raw;if(v>=1e6)return (v/1e6).toFixed(1)+' M';if(v>=1e3)return (v/1e3).toFixed(1)+' K';return v}"}}},
             "scales":{"x":{"type":"logarithmic","title":{"text":"Count (log scale)","display":True}}}}
        )
        flow_html = f'<div class="flow-section"><div class="chart-box" style="max-width:800px;margin:0 auto"><canvas id="{flow_chart_id}" style="max-height:380px"></canvas></div></div>'

    html = f'''<!DOCTYPE html>
<html lang="zh">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width,initial-scale=1.0">
<title>MMPV-RNA v2.3 — Pipeline Report</title>
<script>/*! Chart.js | https://www.chartjs.org | MIT License */
{{CHART_JS_INLINE}}</script>
<style>
:root{{
  --bg:#f4f6f9;--card-bg:#fff;--text:#263238;--muted:#78909c;
  --navy:#0d1b3e;--indigo:#1a237e;--blue:#1565c0;--green:#2e7d32;
  --red:#c62828;--amber:#ef6c00;--border:#e0e0e0;
  --shadow:0 2px 8px rgba(0,0,0,.06);--shadow-lg:0 4px 16px rgba(0,0,0,.1);
  --radius:10px;--radius-sm:6px;
}}
*{{box-sizing:border-box;margin:0;padding:0}}
body{{font-family:-apple-system,BlinkMacSystemFont,'Segoe UI','Noto Sans SC',sans-serif;background:var(--bg);color:var(--text);line-height:1.6}}
.container{{max-width:1200px;margin:0 auto;padding:0 20px 40px;padding-left:240px}}
.hero{{background:linear-gradient(135deg,var(--navy) 0%,var(--indigo) 60%,#283593 100%);color:#fff;padding:40px 32px 32px;border-radius:0 0 16px 16px;margin-bottom:28px;position:relative;overflow:hidden}}
.hero::after{{content:'';position:absolute;top:-50%;right:-20%;width:500px;height:500px;background:rgba(255,255,255,.03);border-radius:50%}}
.hero h1{{font-size:26px;font-weight:700;margin-bottom:6px;position:relative;z-index:1}}
.hero .subtitle{{font-size:13px;opacity:.8;position:relative;z-index:1}}
.hero .gen-time{{font-size:12px;opacity:.65;margin-top:6px;position:relative;z-index:1}}
.kpi-row{{display:grid;grid-template-columns:repeat(auto-fit,minmax(180px,1fr));gap:14px;margin-bottom:20px}}
.flow-section{{background:var(--card-bg);border-radius:var(--radius);box-shadow:var(--shadow);padding:16px 22px;margin-bottom:28px}}
.kpi-card{{background:var(--card-bg);border-radius:var(--radius);padding:18px 20px;box-shadow:var(--shadow);text-align:center;transition:transform .15s}}
.kpi-card:hover{{transform:translateY(-2px);box-shadow:var(--shadow-lg)}}
.kpi-icon{{font-size:24px;margin-bottom:6px}}
.kpi-value{{font-size:13px;font-weight:600;color:var(--text);line-height:1.4}}
.kpi-label{{font-size:11px;color:var(--muted);margin-top:4px;text-transform:uppercase;letter-spacing:.5px}}
.sidebar{{position:fixed;top:0;left:0;width:224px;height:100vh;background:var(--card-bg);border-right:1px solid var(--border);box-shadow:2px 0 8px rgba(0,0,0,.04);z-index:200;display:flex;flex-direction:column;padding-top:12px;overflow-y:auto}}
.sb-title{{font-size:11px;font-weight:700;color:var(--muted);text-transform:uppercase;letter-spacing:.8px;padding:8px 16px 6px;margin:0}}
.sb-item{{display:flex;align-items:center;gap:10px;padding:9px 16px;text-decoration:none;color:var(--text);font-size:12.5px;font-weight:500;border-left:3px solid transparent;transition:all .15s}}
.sb-item:hover{{background:#f0f4ff;border-left-color:var(--blue)}}
.sb-dot{{width:10px;height:10px;border-radius:50%;flex-shrink:0}}
.sb-label{{overflow:hidden;text-overflow:ellipsis;white-space:nowrap}}
.sb-pass{{border-left-color:var(--green)}}.sb-fail{{border-left-color:var(--red)}}.sb-skip{{border-left-color:#e0e0e0}}
.sb-pass .sb-dot{{box-shadow:0 0 0 2px #66bb6a}}.sb-fail .sb-dot{{box-shadow:0 0 0 2px #ef5350}}.sb-skip .sb-dot{{box-shadow:0 0 0 2px #bdbdbd}}
.sb-active{{background:#e8eaf6;font-weight:700;border-left-color:var(--indigo)!important}}
.stage{{background:var(--card-bg);border-radius:var(--radius);box-shadow:var(--shadow);margin-bottom:20px;overflow:hidden}}
.stage-pass{{border-left:4px solid var(--green)}}.stage-fail{{border-left:4px solid var(--red)}}.stage-skip{{border-left:4px solid #bdbdbd}}
.stage-header{{display:flex;align-items:center;gap:16px;padding:18px 22px;border-bottom:1px solid var(--border);background:#fafbfc}}
.stage-icon{{width:44px;height:44px;border-radius:var(--radius-sm);display:flex;align-items:center;justify-content:center;font-size:11px;font-weight:800;color:#fff;flex-shrink:0}}
.stage-title{{flex:1;min-width:0}}
.stage-title h2{{font-size:16px;font-weight:700;color:var(--text);margin-bottom:2px}}
.stage-metric{{font-size:12px;color:var(--muted);display:block;overflow:hidden;text-overflow:ellipsis;white-space:nowrap}}
.stage-badge{{font-size:11px;font-weight:700;padding:4px 12px;border-radius:12px;flex-shrink:0}}
.s-pass{{background:#e8f5e9;color:var(--green)}}.s-fail{{background:#fce4ec;color:var(--red)}}.s-skip{{background:#f5f5f5;color:#9e9e9e}}
.stage-charts{{display:grid;gap:16px;padding:18px 22px}}
.chart-box{{background:#fafbfc;border-radius:var(--radius-sm);padding:12px;border:1px solid var(--border);position:relative}}
.chart-scroll{{overflow-x:auto;max-width:100%}}
.pct-btn{{font-size:11px;padding:4px 14px;border:1px solid var(--indigo);border-radius:4px;background:var(--indigo);cursor:pointer;color:#fff;font-weight:600;transition:all .15s}}
.pct-btn:hover{{background:#283593;border-color:#283593}}
.sankey-section{{padding:18px 22px;display:flex;flex-direction:column;gap:16px}}
.sankey-card{{background:#fafbfc;border-radius:var(--radius-sm);padding:14px;border:1px solid var(--border)}}
.sankey-card h3{{font-size:14px;color:var(--indigo);margin-bottom:8px;text-align:center}}
.sankey-card iframe{{display:block;width:100%;border-radius:4px}}
.stage-tables{{padding:0 22px 14px;display:flex;flex-direction:column;gap:8px}}
.stage-table-detail{{border-top:1px solid var(--border);padding:8px 0 4px}}
.stage-table-detail summary{{cursor:pointer;font-size:12px;padding:2px 0;color:var(--muted)}}
.stage-table-detail summary:hover{{color:var(--blue)}}
.table-wrap{{background:var(--card-bg);border-radius:var(--radius);box-shadow:var(--shadow);overflow:hidden;margin-bottom:20px}}
.table-wrap h3{{font-size:15px;padding:16px 22px;border-bottom:1px solid var(--border);color:var(--text)}}
table{{width:100%;border-collapse:collapse;font-size:13px}}
th{{background:#fafbfc;text-align:left;padding:10px 16px;font-size:11px;color:var(--muted);text-transform:uppercase;letter-spacing:.5px;border-bottom:2px solid var(--border)}}
td{{padding:9px 16px;border-bottom:1px solid #f0f0f0}}
.tr-pass td:first-child{{border-left:3px solid var(--green)}}.tr-fail td:first-child{{border-left:3px solid var(--red)}}.tr-skip td:first-child{{border-left:3px solid #ccc}}
.td-num{{text-align:center;color:var(--muted);font-size:11px;width:36px}}
.td-sub{{padding-left:40px!important;color:#78909c;font-size:12px}}
.td-metric{{font-size:12px;color:#607d8b}}
.tb-badge{{display:inline-block;padding:1px 8px;border-radius:10px;font-size:10px;font-weight:700}}
.tb-pass{{background:#e8f5e9;color:var(--green)}}.tb-fail{{background:#fce4ec;color:var(--red)}}.tb-skip{{background:#f5f5f5;color:#9e9e9e}}
.app-table{{width:100%;border-collapse:collapse;font-size:11px;margin-bottom:8px}}
.app-table th{{background:#f5f5f5;padding:6px 8px;font-size:10px;text-align:left;border:1px solid #e0e0e0;white-space:nowrap}}
.app-table td{{padding:4px 8px;border:1px solid #f0f0f0;max-width:200px;overflow:hidden;text-overflow:ellipsis;white-space:nowrap;font-size:11px}}
.ai-section{{background:var(--card-bg);border-radius:var(--radius);box-shadow:var(--shadow);padding:20px 24px;margin-bottom:28px;border-top:4px solid var(--indigo)}}
.ai-header{{display:flex;align-items:center;gap:10px;margin-bottom:14px;border-bottom:1px solid #e8eaf6;padding-bottom:10px}}
.ai-title{{font-size:15px;font-weight:700;color:var(--indigo)}}
.ai-model-tag{{font-size:10px;background:var(--indigo);color:#fff;padding:2px 8px;border-radius:8px}}
.ai-format-tag{{font-size:10px;background:#e8eaf6;color:var(--indigo);padding:2px 8px;border-radius:8px;font-weight:600}}
.ai-brief{{font-size:14px;line-height:1.8;color:#1a237e;background:linear-gradient(135deg,#e8eaf6 0%,#f0f4ff 100%);padding:14px 18px;border-radius:var(--radius-sm);margin-bottom:14px;text-align:justify;font-weight:500}}
.ai-detail-toggle{{margin-top:8px}}
.ai-detail-toggle summary{{cursor:pointer;font-size:13px;font-weight:600;color:var(--indigo);padding:10px 14px;background:#f5f6fa;border-radius:var(--radius-sm);user-select:none;transition:all .15s}}
.ai-detail-toggle summary:hover{{background:#e8eaf6}}
.ai-detail-content{{margin-top:14px}}
.ai-block{{margin-bottom:14px;padding-left:16px;border-left:3px solid #e0e0e0}}
.ai-block:last-child{{margin-bottom:0}}
.ai-block-title{{font-size:12px;font-weight:700;color:var(--indigo);margin-bottom:4px}}
.ai-block-text{{font-size:13px;line-height:1.9;color:#37474f;text-align:justify}}
.ai-bg{{border-left-color:#1565c0}}.ai-methods{{border-left-color:#00897b}}.ai-results{{border-left-color:#ef6c00}}.ai-discussion{{border-left-color:#6a1b9a}}
.ai-error{{color:var(--red);font-size:12px;padding:8px}}
.footer{{text-align:center;padding:24px;color:var(--muted);font-size:11px;line-height:1.8}}
.footer a{{color:var(--blue);text-decoration:none}}
@media(max-width:768px){{
  .sidebar{{display:none}}
  .container{{padding-left:20px}}
  .hero{{padding:24px 20px 20px}}.hero h1{{font-size:20px}}
  .kpi-row{{grid-template-columns:repeat(2,1fr)}}
  .stage-header{{flex-wrap:wrap;gap:10px}}
  .stage-charts{{grid-template-columns:1fr!important}}
}}
@media print{{
  body{{background:#fff;font-size:11px}}
  .hero{{background:#1a237e!important;-webkit-print-color-adjust:exact}}
  .stage,.table-wrap,.kpi-card{{box-shadow:none;border:1px solid #ddd;break-inside:avoid}}
  .sidebar{{display:none}}
  .container{{padding-left:20px}}
}}
</style>
</head>
<body>
<div class="sidebar">
  <div class="sb-title">Pipeline Modules</div>
  {sidebar_items}
</div>
<div class="container">
<div class="hero">
  <h1>MMPV-RNA v2.3 — Pipeline Report</h1>
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
  <strong>MMPV-RNA v2.3</strong> — Metatranscriptomic Virus Discovery Pipeline<br>
  Generated by <code>report_pipeline.py</code> &nbsp;|&nbsp; {gen_time}
</div>
</div>
<script>
{chart_scripts}
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
    // 检测值轴: 横向图(indexAxis='y')值轴是x, 纵向图值轴是y
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
</body></html>'''

    with open(report_dir / "pipeline_report.html", "w", encoding="utf-8") as hf:
        hf.write(html.replace("{CHART_JS_INLINE}", chart_js_inline))


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

    prompt = f"""你是一位病毒宏基因组学领域的研究科学家。请基于以下 MMPV-RNA v2.3 病毒发现管线的分析结果, 撰写双层研究总结: 先写一段精炼的研究简报, 再写一份完整的 IMRaD 详报。严格按标签输出。

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
一段 150-200 字中文精炼摘要。涵盖: 研究目的、样本规模、关键发现(鉴定病毒数、新颖性比例、HQ vOTU数)、主要结论。如果数据中有植物病毒, 必须重点提及。风格类似顶刊 Highlights。

[DETAILED]
按以下 IMRaD 结构撰写详报 (400-600字):

[BACKGROUND]
1-2句。病毒宏基因组学背景 + 本研究目标。

[METHODS]
4-6句。完整列出 MMPV-RNA v2.3 全部阶段:
- 00a: fastp 质控 (Q20/Q30过滤) + clumpify 去重 → 高质量 clean reads
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
        provider = getattr(generate_ai_summary, '_provider', 'openai')
        model = getattr(generate_ai_summary, '_model', 'gpt-4o-mini')
        api_key = getattr(generate_ai_summary, '_api_key', '')
        base_url = getattr(generate_ai_summary, '_base_url', '')

        if not api_key:
            return '<div class="ai-error">⚠ AI 总结需要 --ai-key 参数</div>'

        import urllib.request
        url = base_url or "https://api.openai.com/v1/chat/completions"
        body = _json.dumps({
            "model": model,
            "messages": [
                {"role": "system", "content": "你是病毒宏基因组学研究科学家。严格按要求的双层格式输出: BRIEF(简洁亮点摘要)+DETAILED(IMRaD完整详报)。用中文,客观精确,引用数据数字,不编造。"},
                {"role": "user", "content": prompt}
            ],
            "temperature": 0.3,
            "max_tokens": 2000
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
            "background": ("📖 研究背景", "ai-bg"),
            "methods": ("🔬 完整方法", "ai-methods"),
            "results": ("📊 关键结果", "ai-results"),
            "discussion": ("💡 讨论展望", "ai-discussion"),
        }
        detail_items = ""
        for key in ("background","methods","results","discussion"):
            text = sections.get(key, "")
            if not text: continue
            label, cls = sec_labels[key]
            detail_items += f'<div class="ai-block {cls}"><div class="ai-block-title">{label}</div><div class="ai-block-text">{text}</div></div>'

        return f'''<div class="ai-section" id="ai-summary">
<div class="ai-header">
  <span class="ai-title">🤖 AI 研究简报</span>
  <span class="ai-model-tag">{model}</span>
  <span class="ai-format-tag">Brief + IMRaD</span>
</div>
<div class="ai-brief">{brief}</div>
<details class="ai-detail-toggle">
  <summary>📋 展开详细报告 (IMRaD)</summary>
  <div class="ai-detail-content">{detail_items}</div>
</details>
</div>'''
    except Exception as e:
        print(f"  [WARN] AI Summary failed: {e}")
        return f'<div class="ai-error">⚠ AI 总结生成失败: {e}</div>'


def main():
    p = argparse.ArgumentParser(description="MMPV-RNA v2.3 — 独立报告生成器")
    p.add_argument("-o", "--output-dir", required=True, help="流水线输出根目录 (包含 00a_CleanData/ ... 09_Reports/)")
    p.add_argument("--skip-sankey", action="store_true", help="跳过 Sankey 图生成")
    p.add_argument("--skip-html", action="store_true", help="仅生成 TSV, 不生成 HTML")
    p.add_argument("--blast-db", help="BLAST 参考数据库路径 (用于序列相似度分类)")
    p.add_argument("--ai-summary", action="store_true", help="生成 AI 管线总结 (需 --ai-key)")
    p.add_argument("--ai-provider", default="openai", choices=["openai","ollama","deepseek","custom"], help="AI 提供商 (default: openai)")
    p.add_argument("--ai-model", default="gpt-4o-mini", help="模型名 (default: gpt-4o-mini)")
    p.add_argument("--ai-key", default="", help="API Key (或 ollama 时留空)")
    p.add_argument("--ai-base-url", default="", help="自定义 API 地址 (如 http://localhost:11434/v1/chat/completions)")
    args = p.parse_args()

    root = Path(args.output_dir).resolve()
    if not root.is_dir():
        sys.exit(f"ERROR: 目录不存在: {root}")

    report_dir = root / "09_Reports"
    report_dir.mkdir(parents=True, exist_ok=True)

    print(f"{'='*60}")
    print(f"Report Pipeline v2.3")
    print(f"  Output: {root}")
    print(f"  Reports: {report_dir}")
    print(f"{'='*60}")

    # 1. 收集数据 + 生成 TSV
    print("\n[1/4] Collecting stage data...")
    t0 = time.time()
    stage_stats = collect_data(root, report_dir, args.blast_db)

    # 2. 写入 stage_summary.tsv + 目录树
    with open(report_dir / "stage_summary.tsv", "w", newline="") as sf:
        w = csv.DictWriter(sf, fieldnames=["Stage","Status","Key_Metric","Details"], delimiter="\t")
        w.writeheader()
        for s in stage_stats: w.writerow(s)

    tree_file = report_dir / "directory_tree.txt"
    with open(tree_file, "w") as tf:
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
            else:
                generate_ai_summary._base_url = ""
        else:
            generate_ai_summary._api_key = ""  # 禁用
        print("\n[3/4] Generating HTML report...")
        t0 = time.time()
        write_html_report(report_dir, stage_stats)
        print(f"  pipeline_report.html — {time.time()-t0:.0f}s")
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
