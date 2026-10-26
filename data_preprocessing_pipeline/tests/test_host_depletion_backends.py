"""host_depletion.py 分类后端 (kraken2 / kunpeng) 回归测试。

覆盖:
  · Kun-peng per-read 输出解析 (C/U 行、taxid 0、/1 /2 后缀剥离)
  · MGI / Illumina 两种 header 的配对后缀探测
  · Kun-peng 库就绪判定 (hash_*.k2d + hash_config.k2d + taxo/opts)
  · keep-taxids 救回语义: 病毒 taxid 不算宿主, unclassified 也不算宿主
  · mfilter=False 提取宿主时只留 C 行宿主 reads
"""
import gzip
import importlib.util
from pathlib import Path

_SPEC = importlib.util.spec_from_file_location(
    "host_depletion", Path(__file__).resolve().parents[1] / "host_depletion.py")
hd = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(hd)


def _write(p: Path, text: str, gz: bool = False):
    p.parent.mkdir(parents=True, exist_ok=True)
    if gz:
        with gzip.open(p, "wt", encoding="utf-8") as f:
            f.write(text)
    else:
        p.write_text(text, encoding="utf-8")


def _fq(ids):
    return "".join(f"@{i}\nACGT\n+\nIIII\n" for i in ids)


def _read_ids(path: Path):
    out = []
    lines = path.read_text(encoding="utf-8").splitlines()
    for i in range(0, len(lines), 4):
        tok = lines[i].split(None, 1)[0]
        out.append(tok[1:])
    return out


def _ready_db(p: Path):
    p.mkdir(parents=True, exist_ok=True)
    for name in ("hash_config.k2d", "opts.k2d", "taxo.k2d"):
        (p / name).write_bytes(b"x")
    (p / "hash_1.k2d").write_bytes(b"x")


def test_parse_kunpeng_output_lines():
    assert hd.process_kraken_output_line(
        "C\tread1\t1000\t150\t1000:12 1000:9") == (1000, "read1")
    assert hd.process_kraken_output_line(
        "U\tread2\t0\t150\t0:21") == (0, "read2")
    # 列数不足 → -1, 空 ID
    assert hd.process_kraken_output_line("C\tread3\t1000") == (-1, "")


def test_strip_and_detect_read_suffix(tmp_path):
    assert hd.strip_read_suffix("LH00278:1:1/1") == "LH00278:1:1"
    assert hd.strip_read_suffix("LH00278:1:1") == "LH00278:1:1"

    mgi = tmp_path / "mgi_R1.fastq"
    _write(mgi, _fq(["a/1", "b/1", "c/1"]))
    assert hd.detect_read_suffix(mgi) == "/1"

    illumina = tmp_path / "illumina_R1.fastq"
    _write(illumina, "".join(f"@r{i} 1:N:0:0\nACGT\n+\nIIII\n" for i in range(3)))
    assert hd.detect_read_suffix(illumina) == ""

    gz = tmp_path / "mgi_R2.fastq.gz"
    _write(gz, _fq(["a/2", "b/2"]), gz=True)
    assert hd.detect_read_suffix(gz) == "/2"


def test_write_suffixed_patterns(tmp_path):
    src = tmp_path / "ids.txt"
    src.write_text("read1\nread2\n", encoding="utf-8")
    assert hd.write_suffixed_patterns(src, tmp_path / "x.txt", "") == src
    dst = hd.write_suffixed_patterns(src, tmp_path / "ids_1.txt", "/1")
    assert Path(dst).read_text(encoding="utf-8") == "read1/1\nread2/1\n"


def test_kunpeng_db_ready(tmp_path):
    db = tmp_path / "db"
    assert not hd.kunpeng_db_ready(db)
    _ready_db(db)
    assert hd.kunpeng_db_ready(db)
    (db / "hash_1.k2d").unlink()
    assert not hd.kunpeng_db_ready(db)


def test_prepare_kunpeng_db_auto_hashshard(tmp_path, monkeypatch):
    db = tmp_path / "k2db"
    db.mkdir()
    (db / "hash.k2d").write_bytes(b"x")
    (db / "opts.k2d").write_bytes(b"x")
    (db / "taxo.k2d").write_bytes(b"x")
    calls = []

    def fake_run(cmd, check=False):
        calls.append(cmd)
        (db / "hash_config.k2d").write_bytes(b"x")
        (db / "hash_1.k2d").write_bytes(b"x")

    monkeypatch.setattr(hd.subprocess, "run", fake_run)
    assert hd.prepare_kunpeng_db(str(db), "/opt/bin/kun_peng", "1G") == str(db)
    assert calls and calls[0][:3] == ["/opt/bin/kun_peng", "hashshard", "--db"]
    assert hd.kunpeng_db_ready(db)


def test_prepare_kunpeng_db_half_converted_exits(tmp_path):
    db = tmp_path / "k2db"
    db.mkdir()
    (db / "hash.k2d").write_bytes(b"x")
    (db / "hash_config.k2d").write_bytes(b"x")
    import pytest
    with pytest.raises(SystemExit):
        hd.prepare_kunpeng_db(str(db))


def _fake_seqkit(seqkit, src, dst, pattern_file, invert, threads=4):
    """用 Python 模拟 seqkit grep [-v] -f：按首 token 精确匹配。"""
    pats = set(Path(pattern_file).read_text(encoding="utf-8").splitlines())
    ids = _read_ids(Path(src))
    keep = [i for i in ids if ((i not in pats) if invert else (i in pats))]
    Path(dst).parent.mkdir(parents=True, exist_ok=True)
    Path(dst).write_text(_fq(keep), encoding="utf-8")


def test_kunpeng_split_rescues_virus_and_drops_host(tmp_path, monkeypatch):
    db = tmp_path / "db"
    _ready_db(db)
    r1 = tmp_path / "S1_R1.fastq"
    r2 = tmp_path / "S1_R2.fastq"
    _write(r1, _fq(["host1/1", "virus1/1", "unclass1/1", "host2/1"]))
    _write(r2, _fq(["host1/2", "virus1/2", "unclass1/2", "host2/2"]))
    report = tmp_path / "report.txt"

    def fake_execute(cmds, **kwargs):
        out_dir = Path(cmds[0][cmds[0].index("--output-dir") + 1])
        out_dir.mkdir(parents=True, exist_ok=True)
        (out_dir / "output_1.txt").write_text(
            "C\thost1\t1000\t150\t1000:12\n"
            "C\tvirus1\t10239\t150\t10239:12\n"
            "U\tunclass1\t0\t150\t0:21\n"
            "C\thost2\t1000\t150\t1000:12\n",
            encoding="utf-8")
        (out_dir / "output_1.kreport2").write_text(
            "100.00\t4\t0\tR\t1\troot\n"
            "50.00\t2\t2\tS\t1000\t  Host\n"
            "25.00\t1\t1\tS\t10239\t  Viruses\n",
            encoding="utf-8")
        return [0], 0.0, 0.0

    monkeypatch.setattr(hd, "execute", fake_execute)
    monkeypatch.setattr(hd, "seqkit_filter_by_ids", _fake_seqkit)

    temp = tmp_path / "tmp"
    temp.mkdir()
    t1, t2, stats = hd.kunpeng_classify_and_split(
        str(db), str(r1), str(r2), threads=2, confidence=0.2,
        keep_taxids_list=[10239], mfilter=True, temp_dir=str(temp),
        report_file=str(report), ext_tmp=".fastq", prefix="[t]",
        sample_name="S1", tool_logs_dir=str(tmp_path / "logs"))

    assert stats["records"] == 4 and stats["host"] == 2 and stats["rescued"] == 1
    # 病毒 + unclassified 保留, 两个宿主 reads 被剔除
    assert _read_ids(Path(t1)) == ["virus1/1", "unclass1/1"]
    assert _read_ids(Path(t2)) == ["virus1/2", "unclass1/2"]
    assert report.is_file()


def test_kunpeng_split_extract_host_only(tmp_path, monkeypatch):
    db = tmp_path / "db"
    _ready_db(db)
    r1 = tmp_path / "S1.fastq"
    _write(r1, _fq(["host1", "unclass1"]))
    report = tmp_path / "report.txt"

    def fake_execute(cmds, **kwargs):
        out_dir = Path(cmds[0][cmds[0].index("--output-dir") + 1])
        out_dir.mkdir(parents=True, exist_ok=True)
        (out_dir / "output_1.txt").write_text(
            "C\thost1\t1000\t150\t1000:12\nU\tunclass1\t0\t150\t0:21\n",
            encoding="utf-8")
        return [0], 0.0, 0.0

    monkeypatch.setattr(hd, "execute", fake_execute)
    monkeypatch.setattr(hd, "seqkit_filter_by_ids", _fake_seqkit)

    temp = tmp_path / "tmp"
    temp.mkdir()
    t1, t2, stats = hd.kunpeng_classify_and_split(
        str(db), str(r1), None, threads=2, confidence=0.2,
        keep_taxids_list=[], mfilter=False, temp_dir=str(temp),
        report_file=str(report), ext_tmp=".fastq", prefix="[t]",
        sample_name="S1", tool_logs_dir=str(tmp_path / "logs"))
    assert t2 is None and stats["host"] == 1
    assert _read_ids(Path(t1)) == ["host1"]
