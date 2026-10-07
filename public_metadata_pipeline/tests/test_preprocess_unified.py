"""preprocess_unified 的 convert 归集逻辑单元测试 (审核 F13)。

sra2fastx 是**就地转换**（fasterq-dump -O 指向 .sra 所在目录, 产物落在原地,
成功后删 .sra），preprocess_unified 必须把它归集到 00_Converted/ —— 否则
clean 阶段必然"未找到任何 fastq"（README 的一条命令路径断在这里）。

这里不跑 sra2fastx（需要 sratoolkit），只用假 .sra 与假产物验证移动逻辑的
单元路径: 产物被搬进 outdir、同名不互相覆盖、outdir 里有 reads。

Run:  python -m pytest public_metadata_pipeline/tests/test_preprocess_unified.py -v
"""
import shutil
import sys
import tempfile
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from preprocess_unified import collect_converted_files, has_reads  # noqa: E402


class _Log:
    def __init__(self):
        self.lines = []

    def info(self, *a, **k):
        self.lines.append(("info", a))

    def warning(self, *a, **k):
        self.lines.append(("warning", a))


@pytest.fixture
def workdir():
    d = tempfile.mkdtemp(prefix="preprocess_test_")
    yield Path(d)
    shutil.rmtree(d, ignore_errors=True)


class TestCollectConvertedFiles:
    def test_moves_produced_fastq_into_outdir(self, workdir):
        """就地产物 (含子目录) 全部归集到 outdir, 且 outdir 成为有效 reads 目录。"""
        sra_dir = workdir / "sra"
        out = workdir / "00_Converted"
        (sra_dir / "sub").mkdir(parents=True)
        (sra_dir / "SRR100.sra").write_bytes(b"")
        (sra_dir / "sub" / "SRR200.sra").write_bytes(b"")
        # sra2fastx 成功后的世界: .sra 已删, fastq.gz 落在原地 (含子目录)
        (sra_dir / "SRR100_1.fastq.gz").write_bytes(b"a")
        (sra_dir / "SRR100_2.fastq.gz").write_bytes(b"b")
        (sra_dir / "sub" / "SRR200.fastq.gz").write_bytes(b"c")

        moved = collect_converted_files(
            [sra_dir / "SRR100.sra", sra_dir / "sub" / "SRR200.sra"], out, _Log())

        assert moved == 3
        assert sorted(p.name for p in out.glob("*.gz")) == [
            "SRR100_1.fastq.gz", "SRR100_2.fastq.gz", "SRR200.fastq.gz"]
        assert not list(sra_dir.rglob("*.fastq.gz")), "产物必须离开原目录"
        assert has_reads(out), "归集后 outdir 必须能被 has_reads 认定为有 reads"

    def test_same_name_from_subdirs_does_not_overwrite(self, workdir):
        """不同子目录的同名样本: 加父目录前缀, 两份产物都要保住。"""
        sra_dir = workdir / "sra"
        out = workdir / "00_Converted"
        for sub in ("runA", "runB"):
            (sra_dir / sub).mkdir(parents=True)
            (sra_dir / sub / "SRR300.sra").write_bytes(b"")
            (sra_dir / sub / "SRR300.fastq.gz").write_bytes(sub.encode())

        moved = collect_converted_files(
            [sra_dir / "runA" / "SRR300.sra", sra_dir / "runB" / "SRR300.sra"],
            out, _Log())

        assert moved == 2
        assert (out / "SRR300.fastq.gz").read_bytes() == b"runA"
        assert (out / "runB_SRR300.fastq.gz").read_bytes() == b"runB"

    def test_nothing_produced_is_reported(self, workdir):
        """没有产物时返回 0 (调用方据此把 convert 判失败, 不让 clean 空跑)。"""
        sra_dir = workdir / "sra"
        sra_dir.mkdir(parents=True)
        (sra_dir / "SRR400.sra").write_bytes(b"")
        assert collect_converted_files([sra_dir / "SRR400.sra"],
                                       workdir / "00_Converted", _Log()) == 0
        assert not has_reads(workdir / "00_Converted")
