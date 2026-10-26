"""preprocess_report.py 解析/汇总回归测试。

覆盖 2026-10-12 审查修复:
  P1-1 去宿主输出带 R 标签 (_clean_R1/_clean_R2) 不得拆成幻影样本
  P1-2/P3-14 交接清单每样本只留一套完整 PE/SE, 优先 bbnorm > 去宿主 > clean
  P1-3 bowtie2 真实输出 '<pct>% overall alignment rate' (百分比在前) 可解析
  P1-5 seqkit PE 两端逐行 → 按端累加为合计
  P2-6 表格/JSON 转义; P2-10 中位数偶数取均值
  (2026-10-13) seqkit 重跑临时随机路径批次不双倍累加 (只取最后一批;
  FINAL 稳定路径行为不回归); 比对率样本键带点样本名不被截断
  (2026-10-16) 分节 HTML + 侧栏导航; 表尾 TOTAL 行; Fastp_Dup(%) 列;
  每样本 fastp / kraken2 报告与资源表的相对链接 (缺失不渲染)
  (2026-10-16 晚) 版式改回 ③ report_pipeline 老 stage 卡体系:
  #stage-CleanData/#stage-HostDep/[#stage-BBNorm] + sb-item 侧栏 + chart-box
  (chart_s00b + Show % 切换) + app-table/tfoot + 页尾 Handoff details
"""
import importlib.util
import json
from pathlib import Path

import pytest

_SPEC = importlib.util.spec_from_file_location(
    "preprocess_report", Path(__file__).resolve().parents[1] / "preprocess_report.py")
pr = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(pr)


def _touch(p: Path):
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_bytes(b"\x1f\x8b")          # 内容无关: 报告只解析文件名


# ── P1-1 / P1-2 / P3-14 ──────────────────────────────────────────

def test_r_tag_no_phantom_sample(tmp_path):
    """输入带 R 标签时 host_depletion 产物为 X_clean_R1/_R2, 不得拆成两个样本。"""
    dep = tmp_path / "dep"
    for n in ("ERR1_clean_R1.fa.gz", "ERR1_clean_R2.fa.gz"):
        _touch(dep / n)
    files = pr.collect_final_files(dep, None)
    assert set(files) == {"ERR1"}
    assert files["ERR1"]["source"] == "depleted"
    assert [p.name for p in files["ERR1"]["files"]] == ["ERR1_clean_R1.fa.gz", "ERR1_clean_R2.fa.gz"]


def test_one_set_per_sample_depleted_beats_clean(tmp_path):
    """同一样本在去宿主与 clean 目录都有产物时, 清单只留去宿主一套。"""
    dep, clean = tmp_path / "dep", tmp_path / "clean"
    _touch(dep / "S1_clean_1.fa.gz"); _touch(dep / "S1_clean_2.fa.gz")
    _touch(clean / "3.clumpify" / "S1_1.fa.gz"); _touch(clean / "3.clumpify" / "S1_2.fa.gz")
    _touch(clean / "2.fasta" / "S1_1.fa.gz"); _touch(clean / "2.fasta" / "S1_2.fa.gz")
    files = pr.collect_final_files(dep, clean)
    assert files["S1"]["source"] == "depleted"
    assert len(files["S1"]["files"]) == 2
    assert all(p.parent == dep for p in files["S1"]["files"])


def test_clean_fallback_when_depleted_missing(tmp_path):
    """去宿主无产物 (或 PE 不完整) 时回退 clean 的 3.clumpify。"""
    clean = tmp_path / "clean"
    _touch(clean / "3.clumpify" / "S2_1.fa.gz"); _touch(clean / "3.clumpify" / "S2_2.fa.gz")
    _touch(clean / "3.clumpify" / "S3_1.fa.gz")          # PE 只有 R1 → 不算完整
    files = pr.collect_final_files(tmp_path / "nodep", clean)
    assert files["S2"]["source"] == "clean" and len(files["S2"]["files"]) == 2
    assert "S3" not in files


def test_bbnorm_highest_priority_pe_and_se(tmp_path):
    """bbnorm 产物优先于去宿主; PE 名残留 _clean, SE 名残留扩展名 (run_bbnorm 命名)。"""
    dep, bb = tmp_path / "dep", tmp_path / "bb"
    _touch(dep / "S1_clean_R1.fa.gz"); _touch(dep / "S1_clean_R2.fa.gz")
    _touch(dep / "S2_clean.fa.gz")
    _touch(bb / "S1_clean_norm_R1.fq.gz"); _touch(bb / "S1_clean_norm_R2.fq.gz")
    _touch(bb / "S2_clean.fa_norm_SE.fq.gz")
    files = pr.collect_final_files(dep, None, bb)
    assert files["S1"]["source"] == "normalized"
    assert [p.name for p in files["S1"]["files"]] == ["S1_clean_norm_R1.fq.gz", "S1_clean_norm_R2.fq.gz"]
    assert files["S2"]["source"] == "normalized" and len(files["S2"]["files"]) == 1


def test_rrna_subdir_excluded(tmp_path):
    """--keep_rrna 的 rrna/ 子目录产物不是最终交付 reads。"""
    dep = tmp_path / "dep"
    _touch(dep / "rrna" / "S1.rrna.1.fq.gz")
    _touch(dep / "S1_clean_1.fa.gz"); _touch(dep / "S1_clean_2.fa.gz")
    files = pr.collect_final_files(dep, None)
    assert [p.name for p in files["S1"]["files"]] == ["S1_clean_1.fa.gz", "S1_clean_2.fa.gz"]


# ── P1-3 比对率 ──────────────────────────────────────────────────

def test_aligner_rate_bowtie2_real_format(tmp_path):
    """bowtie2 真实排版 '95.50% overall alignment rate' (百分比在前) + samtools 日志不覆盖。"""
    dep = tmp_path / "dep"; (dep / "logs").mkdir(parents=True)
    (dep / "logs" / "S1.step2_align_filter.0_bowtie2.log").write_text(
        "=== 命令: bowtie2 -x host ===\n95.50% overall alignment rate\n", encoding="utf-8")
    (dep / "logs" / "S1.step2_align_filter.1_samtools.log").write_text(
        "=== 命令: samtools view ===\n", encoding="utf-8")
    assert pr.collect_aligner(dep) == {"S1": 95.5}


def test_aligner_rate_legacy_format_and_unparsed(tmp_path):
    """旧排版 'overall alignment rate: 87.0%' 兼容; 日志无该行 → UNPARSED (不留空)。"""
    dep = tmp_path / "dep"; (dep / "logs").mkdir(parents=True)
    (dep / "logs" / "A.step2_align_filter.0_hisat2.log").write_text(
        "overall alignment rate: 87.0%\n", encoding="utf-8")
    (dep / "logs" / "B.step2_align_filter.0_minimap2.log").write_text("no stats here\n", encoding="utf-8")
    out = pr.collect_aligner(dep)
    assert out["A"] == 87.0 and out["B"] == "UNPARSED"


def test_aligner_rate_dotted_sample_name(tmp_path):
    """样本名自带点 (如 ERR1.2): 按 .step2_align_filter 标记 rsplit, 不截断成 ERR1。"""
    dep = tmp_path / "dep"; (dep / "logs").mkdir(parents=True)
    (dep / "logs" / "ERR1.2.step2_align_filter.0_bowtie2.log").write_text(
        "93.10% overall alignment rate\n", encoding="utf-8")
    (dep / "logs" / "ERR1.2.step2_align_filter.1_samtools.log").write_text(
        "samtools view: no rate line\n", encoding="utf-8")
    out = pr.collect_aligner(dep)
    assert out == {"ERR1.2": 93.1}


# ── P1-5 seqkit 两端 ─────────────────────────────────────────────

def test_seqkit_pe_accumulates_both_ends(tmp_path):
    """PE 每 stage 两行 (R1/R2) 累加为合计; 同 file 重跑 append 不重复累加。"""
    dep = tmp_path / "dep"; (dep / "logs").mkdir(parents=True)
    hdr = "Sample\tStage\tfile\tformat\ttype\tnum_seqs\tsum_len\tmin_len\tavg_len\tmax_len"
    rows = [
        f"{hdr}",
        "S1\t1_Raw\tS1_1.fa.gz\tFASTA\tDNA\t1000\t0\t0\t100\t100",
        "S1\t1_Raw\tS1_2.fa.gz\tFASTA\tDNA\t1000\t0\t0\t100\t100",
        "S1\t1_Raw\tS1_1.fa.gz\tFASTA\tDNA\t1000\t0\t0\t100\t100",   # 重跑重复行
        "S1\t3_Host_Filtered\tS1_h1.fa.gz\tFASTA\tDNA\t600\t0\t0\t100\t100",
        "S1\t3_Host_Filtered\tS1_h2.fa.gz\tFASTA\tDNA\t400\t0\t0\t100\t100",
    ]
    (dep / "logs" / "host_depletion_seqkit_summary.tsv").write_text("\n".join(rows) + "\n", encoding="utf-8")
    out = pr.collect_seqkit(dep)
    assert out["S1"]["input"] == 2000 and out["S1"]["host"] == 1000


def test_rows_retained_uses_both_ends(tmp_path):
    """build_rows 的留存率 = 两端合计的 host / 两端合计的 input。"""
    dep = tmp_path / "dep"; (dep / "logs").mkdir(parents=True)
    hdr = "Sample\tStage\tfile\tformat\ttype\tnum_seqs\tsum_len\tmin_len\tavg_len\tmax_len"
    (dep / "logs" / "host_depletion_seqkit_summary.tsv").write_text("\n".join([
        hdr,
        "S1\t1_Raw\ta\tFASTA\tDNA\t1000\t0\t0\t100\t100",
        "S1\t1_Raw\tb\tFASTA\tDNA\t1000\t0\t0\t100\t100",
        "S1\t3_Host_Filtered\tc\tFASTA\tDNA\t400\t0\t0\t100\t100",
        "S1\t3_Host_Filtered\td\tFASTA\tDNA\t200\t0\t0\t100\t100",
    ]) + "\n", encoding="utf-8")
    _touch(dep / "S1_clean_1.fa.gz"); _touch(dep / "S1_clean_2.fa.gz")
    row, = pr.build_rows(None, dep, 20.0)
    assert row["Depletion_Input"] == 2000 and row["Retained(%)"] == 30.0
    assert row["Final_Files"] == 2 and row["Final_Source"] == "depleted"


def test_seqkit_rerun_temp_batch_not_double_counted(tmp_path):
    """失败重跑: 中间产物落在新临时目录 → 同样本两批记录不双计, 只取最后一批。

    旧实现按 file 列去重, 临时随机路径两批文件名不同 → After_Kraken2/
    After_HostAlign 双倍累加。"""
    dep = tmp_path / "dep"; (dep / "logs").mkdir(parents=True)
    hdr = "Sample\tStage\tfile\tformat\ttype\tnum_seqs\tsum_len\tmin_len\tavg_len\tmax_len"
    rows = [
        hdr,
        # 第一批 (失败中断): 临时目录 tmp_a
        "S1\t2_Kraken2_Filtered\t/tmp/tmp_a/out_1.fasta\tFASTA\tDNA\t900\t0\t0\t100\t100",
        "S1\t2_Kraken2_Filtered\t/tmp/tmp_a/out_2.fasta\tFASTA\tDNA\t900\t0\t0\t100\t100",
        "S1\t3_Host_Filtered\t/tmp/tmp_a/S1_hostclean.1.fasta\tFASTA\tDNA\t500\t0\t0\t100\t100",
        "S1\t3_Host_Filtered\t/tmp/tmp_a/S1_hostclean.2.fasta\tFASTA\tDNA\t500\t0\t0\t100\t100",
        # 第二批 (重跑成功): 新临时目录 tmp_b → 只算这批
        "S1\t2_Kraken2_Filtered\t/tmp/tmp_b/out_1.fasta\tFASTA\tDNA\t950\t0\t0\t100\t100",
        "S1\t2_Kraken2_Filtered\t/tmp/tmp_b/out_2.fasta\tFASTA\tDNA\t950\t0\t0\t100\t100",
        "S1\t3_Host_Filtered\t/tmp/tmp_b/S1_hostclean.1.fasta\tFASTA\tDNA\t450\t0\t0\t100\t100",
        "S1\t3_Host_Filtered\t/tmp/tmp_b/S1_hostclean.2.fasta\tFASTA\tDNA\t450\t0\t0\t100\t100",
    ]
    (dep / "logs" / "host_depletion_seqkit_summary.tsv").write_text("\n".join(rows) + "\n", encoding="utf-8")
    out = pr.collect_seqkit(dep)
    assert out["S1"]["kraken2"] == 1900 and out["S1"]["host"] == 900


def test_seqkit_final_stage_stable_path_rerun_no_regression(tmp_path):
    """FINAL 阶段 (4_rRNA_Filtered, 稳定终产物路径): 重跑同名行只算最后一次, 不双倍。"""
    dep = tmp_path / "dep"; (dep / "logs").mkdir(parents=True)
    hdr = "Sample\tStage\tfile\tformat\ttype\tnum_seqs\tsum_len\tmin_len\tavg_len\tmax_len"
    rows = [
        hdr,
        "S1\t4_rRNA_Filtered\tS1_clean_1.fa.gz\tFASTA\tDNA\t700\t0\t0\t100\t100",
        "S1\t4_rRNA_Filtered\tS1_clean_2.fa.gz\tFASTA\tDNA\t700\t0\t0\t100\t100",
        "S1\t4_rRNA_Filtered\tS1_clean_1.fa.gz\tFASTA\tDNA\t800\t0\t0\t100\t100",
        "S1\t4_rRNA_Filtered\tS1_clean_2.fa.gz\tFASTA\tDNA\t800\t0\t0\t100\t100",
    ]
    (dep / "logs" / "host_depletion_seqkit_summary.tsv").write_text("\n".join(rows) + "\n", encoding="utf-8")
    out = pr.collect_seqkit(dep)
    assert out["S1"]["rrna"] == 1600


# ── P2-6 / P2-10 ─────────────────────────────────────────────────

def test_median_even_takes_mean():
    assert pr._median([60.0, 30.0]) == 45.0
    assert pr._median([10.0, 20.0, 60.0]) == 20.0
    assert pr._median([]) is None


def test_generate_html_escapes_cells_and_json(tmp_path):
    """样本名进表格要 HTML 转义; 进 <script> 的 JSON 要转义 '<' (不产生裸 script)。

    表尾带 TOTAL 行 (<tfoot>), Fastp_Dup(%) 列与 TSV 同步出现在表头。
    """
    evil = "A&B<script>alert(1)</script>"
    rows = [{"Sample": evil, "Raw_Reads": 100, "Fastp_Clean_Reads": 90, "Fastp_Q30(%)": 90.0,
             "Fastp_Dup(%)": 12.3, "Depletion_Input": 90, "After_Kraken2": None,
             "After_HostAlign": 80, "After_rRNA": None, "Retained(%)": 88.9,
             "Aligner_Rate(%)": "UNPARSED", "Final_Files": 2, "Final_Source": "depleted",
             "Status": "PASS"}]
    out = tmp_path / "r.html"
    pr.generate_html(rows, out, 20.0)
    html = out.read_text(encoding="utf-8")
    assert evil not in html
    assert "&lt;script&gt;" in html
    assert "\\u003cscript" in html
    assert json.dumps({"k": 1}) and "var CD=" in html
    assert "<th>Final_Source</th>" in html and "<th>Aligner_Rate(%)</th>" in html
    assert "<th>Fastp_Dup(%)</th>" in html
    assert "<tfoot>" in html and ">TOTAL<" in html


def test_summary_tsv_has_source_column(tmp_path):
    """TSV 含来源列, TOTAL 行来源留空 (求和列不含字符串列)。"""
    rows = [{"Sample": "S1", "Raw_Reads": 10, "Fastp_Clean_Reads": None, "Fastp_Q30(%)": None,
             "Fastp_Dup(%)": None, "Depletion_Input": 10, "After_Kraken2": None,
             "After_HostAlign": 5, "After_rRNA": None,
             "Retained(%)": 50.0, "Aligner_Rate(%)": 95.5, "Final_Files": 2,
             "Final_Source": "normalized", "Status": "PASS"}]
    out = tmp_path / "s.tsv"
    pr.write_summary_tsv(rows, out)
    lines = [l for l in out.read_text(encoding="utf-8").splitlines() if l]
    assert lines[0].split("\t") == pr.SUMMARY_COLS
    assert lines[1].split("\t")[pr.SUMMARY_COLS.index("Final_Source")] == "normalized"
    total = lines[-1].split("\t")
    assert total[0] == "TOTAL" and total[pr.SUMMARY_COLS.index("Final_Source")] == ""
    assert total[pr.SUMMARY_COLS.index("Fastp_Dup(%)")] == ""    # 比率列不求和


# ── fastp 重复率列 / 分节侧栏 / 明细链接 (2026-10-16 审计修复) ──────

def _mk_row(sample="S1", source="depleted", status="PASS", **kw):
    row = {"Sample": sample, "Raw_Reads": 1000, "Fastp_Clean_Reads": 900, "Fastp_Q30(%)": 93.0,
           "Fastp_Dup(%)": 7.4, "Depletion_Input": 900, "After_Kraken2": None,
           "After_HostAlign": 800, "After_rRNA": None, "Retained(%)": 88.9,
           "Aligner_Rate(%)": 95.5, "Final_Files": 2, "Final_Source": source,
           "Status": status}
    row.update(kw)
    return row


def test_fastp_dup_rate_into_rows_and_tsv(tmp_path):
    """fastp duplication.rate 不再丢弃: build_rows / SUMMARY_COLS / TSV 都带 Fastp_Dup(%)。"""
    clean = tmp_path / "clean"; (clean / "logs").mkdir(parents=True)
    (clean / "logs" / "S1_fastp_report.json").write_text(json.dumps({
        "summary": {"before_filtering": {"total_reads": 1000},
                    "after_filtering": {"total_reads": 900, "q30_rate": 0.93}},
        "duplication": {"rate": 0.074},
    }), encoding="utf-8")
    _touch(clean / "3.clumpify" / "S1_1.fa.gz"); _touch(clean / "3.clumpify" / "S1_2.fa.gz")
    row, = pr.build_rows(clean, None, 20.0)
    assert "Fastp_Dup(%)" in pr.SUMMARY_COLS
    assert row["Fastp_Q30(%)"] == 93.0 and row["Fastp_Dup(%)"] == 7.4
    out = tmp_path / "s.tsv"
    pr.write_summary_tsv([row], out)
    lines = [l for l in out.read_text(encoding="utf-8").splitlines() if l]
    i = pr.SUMMARY_COLS.index("Fastp_Dup(%)")
    assert lines[0].split("\t")[i] == "Fastp_Dup(%)"
    assert lines[1].split("\t")[i] == "7.4"
    assert lines[2].split("\t")[i] == ""                   # TOTAL 行比率列留空


def test_html_stage_cards_sidebar_and_bbnorm_visibility(tmp_path):
    """老版式锚点: #stage-CleanData/#stage-HostDep (+bbnorm 卡仅在产物出现时) 与侧栏 sb-item。"""
    out = tmp_path / "r.html"
    pr.generate_html([_mk_row()], out, 20.0)
    html = out.read_text(encoding="utf-8")
    assert 'class="sidebar"' in html and '<div class="sb-title">Pipeline Modules</div>' in html
    for sid in ("stage-CleanData", "stage-HostDep"):
        assert f'id="{sid}"' in html and f'href="#{sid}"' in html
    assert 'id="stage-BBNorm"' not in html and 'href="#stage-BBNorm"' not in html
    assert html.count('class="sb-item') == 2 and html.count('class="stage stage-') == 2
    # Handoff 不再是独立侧栏节, 只在页尾 details
    assert 'id="handoff"' in html and 'href="#handoff"' not in html
    # ① 布局相关产物目录存在 → 渲染 + 侧栏出现
    bb = tmp_path / pr._D["d_bbnorm"]
    bb.mkdir()
    pr.generate_html([_mk_row()], out, 20.0)
    html = out.read_text(encoding="utf-8")
    assert 'id="stage-BBNorm"' in html and 'href="#stage-BBNorm"' in html
    assert html.count('class="sb-item') == 3
    # ② 无目录但 Final_Source=normalized → 同样渲染
    bb.rmdir()
    pr.generate_html([_mk_row(source="normalized")], out, 20.0)
    html = out.read_text(encoding="utf-8")
    assert 'id="stage-BBNorm"' in html and "samples normalized" in html


def test_html_old_visual_language_fragments(tmp_path):
    """③ 老报告同款版式: stage 卡徽章/图标色、KPI 行、app-table+tfoot、chart-box
    + Show % 切换、表格分页 JS、scroll-spy、打印媒体查询。"""
    out = tmp_path / "r.html"
    pr.generate_html([_mk_row(status="PASS"), _mk_row(sample="S2", status="LOW_RETAINED")],
                     out, 20.0)
    html = out.read_text(encoding="utf-8")
    # stage 卡 (≥2 张 pass) + 图标缩写/色值严格按老报告
    assert html.count('class="stage stage-pass"') >= 2
    for icon, color in ((">CLN<", "#607d8b"), (">DEP<", "#546e7a")):
        assert icon in html and color in html
    assert 'class="stage-badge s-pass">✓ PASS<' in html
    assert html.count("kpi-card") >= 5 and 'class="kpi-row"' in html
    # 表格: app-table + tfoot TOTAL + 状态徽章 (老 tb-badge)
    assert '<table class="app-table" id="stagetbl_s00b_preprocessing_summary">' in html
    assert "<tfoot>" in html and ">TOTAL<" in html and 'class="tb-badge tb-pass"' in html
    # 图表: chart-box + chart_s00b + Show % 切换 + scroll-spy (③ 原样实现)
    assert 'class="chart-box"' in html and 'id="chart_s00b"' in html
    assert 'class="pct-btn" data-chart="chart_s00b" data-mode="abs"' in html
    assert "toggleStackedPct" in html and "sb-active" in html
    assert "window.addEventListener('scroll'" in html
    assert "new Chart(document.getElementById('chart_s00b')" in html
    # 打印/窄屏: 侧栏隐藏
    assert "@media print" in html and ".sidebar{display:none}" in html
    # 页尾 Handoff details + 侧栏无该项
    assert 'id="handoff"' in html and "assembly_ready.list" in html


def test_html_table_pagination_over_ten_rows(tmp_path):
    """>10 行时走 ③ 同款分页: .pag-row 隐藏行 + 空 pager div + Prev/Next 分页 JS。"""
    rows = [_mk_row(sample=f"S{i:02d}", source="depleted") for i in range(12)]
    out = tmp_path / "r.html"
    pr.generate_html(rows, out, 20.0)
    html = out.read_text(encoding="utf-8")
    assert html.count('class="tr-pass pag-row"') == 12
    assert 'id="pager_stagetbl_s00b_preprocessing_summary"' in html
    assert "document.getElementById('stagetbl_s00b_preprocessing_summary').querySelector('tbody')" in html
    assert "totalPages" in html and "Page <b>" in html
    assert "<tfoot>" in html and ">TOTAL<" in html          # TOTAL 行不受分页影响
    # 3 行 (默认夹具) 不分页
    out2 = tmp_path / "r2.html"
    pr.generate_html([_mk_row()], out2, 20.0)
    html2 = out2.read_text(encoding="utf-8")
    assert "pag-row" not in html2 and "pager_stagetbl_" not in html2


def test_html_detail_links_fastp_kraken2_resource(tmp_path):
    """Clean 段列每样本 fastp 报告; Host Depletion 段列 kraken2 报告 + 资源表
    (均为 os.path.relpath 相对链接); 文件/目录缺失时对应块不渲染。"""
    clean, dep = tmp_path / "clean", tmp_path / "dep"
    (clean / "logs").mkdir(parents=True)
    (dep / "logs" / "kraken2_report").mkdir(parents=True)
    (clean / "logs" / "S1_fastp_report.html").write_text("<html></html>", encoding="utf-8")
    (dep / "logs" / "kraken2_report" / "S1_kraken2_report.txt").write_text("x", encoding="utf-8")
    (dep / "logs" / "host_depletion_resource_usage.tsv").write_text("a\tb\n", encoding="utf-8")
    out = tmp_path / "r.html"
    pr.generate_html([_mk_row()], out, 20.0, None, clean_dir=clean, deplete_dir=dep)
    html = out.read_text(encoding="utf-8")
    assert "Per-sample fastp QC reports (1)" in html
    assert "clean/logs/S1_fastp_report.html" in html
    assert "1 per-sample classification reports" in html
    assert "dep/logs/kraken2_report/S1_kraken2_report.txt" in html
    assert "Resource usage table:" in html
    assert "dep/logs/host_depletion_resource_usage.tsv" in html
    # 未给目录 (旧 3 参调用) → 三个明细块都不渲染
    out2 = tmp_path / "r2.html"
    pr.generate_html([_mk_row()], out2, 20.0)
    html2 = out2.read_text(encoding="utf-8")
    assert "Per-sample fastp QC reports" not in html2
    assert "classification reports" not in html2 and "Resource usage table" not in html2
