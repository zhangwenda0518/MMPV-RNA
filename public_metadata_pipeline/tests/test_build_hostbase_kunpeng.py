"""build_hostbase.py Kun-peng 建库路径的回归测试 (不依赖真实 kun_peng 二进制)。"""
import importlib.util
from pathlib import Path

_SPEC = importlib.util.spec_from_file_location(
    "build_hostbase", Path(__file__).resolve().parents[1] / "build_hostbase.py")
bh = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(bh)


def _taxonomy(p: Path):
    p.mkdir(parents=True, exist_ok=True)
    (p / "nodes.dmp").write_text("1\t|\t1\t|\tno rank\t|\n", encoding="utf-8")
    (p / "names.dmp").write_text("1\t|\troot\t|\t\t|\tscientific name\t|\n", encoding="utf-8")


def _ready_k2_db(p: Path):
    p.mkdir(parents=True, exist_ok=True)
    for name in ("hash.k2d", "opts.k2d", "taxo.k2d"):
        (p / name).write_bytes(b"x")


def test_kunpeng_native_build_commands(tmp_path, monkeypatch):
    fasta = tmp_path / "host.fa"
    fasta.write_text(">chr1 desc\nACGT\n>chr2\nTTTT\n", encoding="utf-8")
    tax = tmp_path / "taxonomy"
    _taxonomy(tax)
    db = tmp_path / "hostdb" / "kunpeng"
    cmds = []

    def fake_exec(cmd, tool_name="SYS"):
        cmds.append(cmd)
        if "add-library" in cmd:
            d = Path(cmd[cmd.index("--db") + 1])
            (d / "seqid2taxid.map").write_text("chr1\t123\nchr2\t123\n", encoding="utf-8")
        elif "build-db" in cmd:
            d = Path(cmd[cmd.index("--db") + 1])
            for name in ("hash_config.k2d", "opts.k2d", "taxo.k2d", "hash_1.k2d"):
                (d / name).write_bytes(b"x")
        return 0

    monkeypatch.setattr(bh, "execute_realtime", fake_exec)
    rc = bh.IndexBuilders.kunpeng(
        str(fasta), str(db), threads=4, taxonomy_path=str(tax), taxid=123,
        hash_capacity="256M", force=False, kunpeng_bin="kun_peng")

    assert rc == 0
    add = [c for c in cmds if "add-library" in c]
    build = [c for c in cmds if "build-db" in c]
    assert add and add[0][0] == "kun_peng"
    assert build and "--hash-capacity" in build[0] and "256M" in build[0]
    # taxid 注入 + taxonomy 复制
    tagged = (db / "prep" / "host_with_taxid.fasta").read_text(encoding="utf-8")
    assert ">chr1|kraken:taxid|123 desc" in tagged
    assert ">chr2|kraken:taxid|123" in tagged
    assert (db / "taxonomy" / "nodes.dmp").is_file()
    assert bh._kunpeng_db_ready(str(db))


def test_kunpeng_from_kraken2_hashshard(tmp_path, monkeypatch):
    k2 = tmp_path / "kraken2"
    _ready_k2_db(k2)
    cmds = []

    def fake_exec(cmd, tool_name="SYS"):
        cmds.append(cmd)
        d = Path(cmd[cmd.index("--db") + 1])
        (d / "hash_config.k2d").write_bytes(b"x")
        (d / "hash_1.k2d").write_bytes(b"x")
        return 0

    monkeypatch.setattr(bh, "execute_realtime", fake_exec)
    rc = bh.IndexBuilders.kunpeng_from_kraken2(
        str(k2), hash_capacity="1G", force=False, kunpeng_bin="kun_peng")

    assert rc == 0
    assert cmds[0][:3] == ["kun_peng", "hashshard", "--db"]
    assert "1G" in cmds[0]
    assert bh._kunpeng_db_ready(str(k2))
