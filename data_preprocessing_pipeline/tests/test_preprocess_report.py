"""preprocess_report.py 解析/汇总回归测试。

覆盖 2026-10-12 审查修复:
  P1-1 去宿主输出带 R 标签 (_clean_R1/_clean_R2) 不得拆成幻影样本
  P1-2/P3-14 交接清单每样本只留一套完整 PE/SE, 优先 bbnorm > 去宿主 > clean
  P1-3 bowtie2 真实输出 '<pct>% overall alignment rate' (百分比在前) 可解析
  P1-5 seqkit PE 两端逐行 → 按端累加为合计
  P2-6 表格/JSON 转义; P2-10 中位数偶数取均值
  (2026-10-13) seqkit 重跑临时随机路径批次不双倍累加 (只取最后一批;
  FINAL 稳定路径行为不回归); 比对率样本键带点样本名不被截断
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
    """样本名进表格要 HTML 转义; 进 <script> 的 JSON 要转义 '<' (不产生裸 script)。"""
    evil = "A&B<script>alert(1)</script>"
    rows = [{"Sample": evil, "Raw_Reads": 100, "Fastp_Clean_Reads": 90, "Fastp_Q30(%)": 90.0,
             "Depletion_Input": 90, "After_Kraken2": None, "After_HostAlign": 80,
             "After_rRNA": None, "Retained(%)": 88.9, "Aligner_Rate(%)": "UNPARSED",
             "Final_Files": 2, "Final_Source": "depleted", "Status": "PASS"}]
    out = tmp_path / "r.html"
    pr.generate_html(rows, out, 20.0)
    html = out.read_text(encoding="utf-8")
    assert evil not in html
    assert "&lt;script&gt;" in html
    assert "\\u003cscript" in html
    assert json.dumps({"k": 1}) and "var CD=" in html
    assert "<th>Final_Source</th>" in html and "<th>Aligner_Rate(%)</th>" in html


def test_summary_tsv_has_source_column(tmp_path):
    """TSV 含来源列, TOTAL 行来源留空 (求和列不含字符串列)。"""
    rows = [{"Sample": "S1", "Raw_Reads": 10, "Fastp_Clean_Reads": None, "Fastp_Q30(%)": None,
             "Depletion_Input": 10, "After_Kraken2": None, "After_HostAlign": 5, "After_rRNA": None,
             "Retained(%)": 50.0, "Aligner_Rate(%)": 95.5, "Final_Files": 2,
             "Final_Source": "normalized", "Status": "PASS"}]
    out = tmp_path / "s.tsv"
    pr.write_summary_tsv(rows, out)
    lines = [l for l in out.read_text(encoding="utf-8").splitlines() if l]
    assert lines[0].split("\t") == pr.SUMMARY_COLS
    assert lines[1].split("\t")[pr.SUMMARY_COLS.index("Final_Source")] == "normalized"
    total = lines[-1].split("\t")
    assert total[0] == "TOTAL" and total[pr.SUMMARY_COLS.index("Final_Source")] == ""
