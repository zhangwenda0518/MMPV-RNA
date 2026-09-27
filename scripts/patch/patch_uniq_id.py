#!/usr/bin/env python3
"""
修 FastTree "Non-unique name" 崩溃 (Geminiviridae rc=1):

根因: 自研 KEEP contig 的 ID 与 DB 参考序列 ID 相同 (如 BK061149.1),
      analysis_sequences.fasta 出现两条同名序列, FastTree 拒绝重名。
      全库扫描确认仅 Geminiviridae 1 例。

修法: 合并 contigs + db_extracted 时做唯一性保证。
      DB 记录优先保留原名 (参考序列必须可辨识),
      contig 记录默认保持原名; 仅当与已出现的 ID 冲突时, 才追加 "sample|" 前缀,
      仍冲突则再加 __2/__3。这样未重名的 contig ID 保持原样, 不破坏血缘追踪。
"""
import shutil
import sys
from pathlib import Path

TREE = Path("/home/zhangwenda/MMPV-RNA/virome_discovery_pipeline/utils/acvirus_tree_pro.py")
BAK = TREE.with_suffix(".py.bak_uniqid_20260910")

if not TREE.is_file():
    sys.exit(f"找不到: {TREE}")

src = TREE.read_text(encoding="utf-8")
orig = src
applied = []

old = """    with final_fasta.open("wb") as out_f:
        if args.contigs and args.contigs.exists():
            with args.contigs.open("rb") as inf: shutil.copyfileobj(inf,out_f)
        if db_extracted.exists():
            with db_extracted.open("rb") as inf: shutil.copyfileobj(inf,out_f)"""
new = """    # 合并 contigs + DB 序列: 保证全局唯一 ID
    # (自研 contig ID 可能与 DB 参考同名, 如 BK061149.1, FastTree 会因重名崩溃)
    _seen_ids = set()

    def _read_fasta_records(_path):
        _rid, _buf = None, []
        with _path.open() as _h:
            for _ln in _h:
                if _ln.startswith(">"):
                    if _rid is not None:
                        yield _rid, _buf
                    _rid = _ln[1:].split()[0].strip()
                    _buf = []
                else:
                    _buf.append(_ln)
            if _rid is not None:
                yield _rid, _buf

    def _uniq(_base, _mark_prefix=False):
        # 默认原名; 冲突时 (或显式要求标记) 加 sample| 前缀, 再冲突加 __2/__3
        _cand = f"sample|{_base}" if (_mark_prefix or _base in _seen_ids) else _base
        if _cand in _seen_ids:
            _i = 2
            while f"{_cand}__{_i}" in _seen_ids:
                _i += 1
            _cand = f"{_cand}__{_i}"
        _seen_ids.add(_cand)
        return _cand

    with final_fasta.open("w") as out_f:
        # DB 先写, 优先占住原始名 (参考序列保持可辨识)
        if db_extracted.exists():
            for _rid, _buf in _read_fasta_records(db_extracted):
                out_f.write(f">{_uniq(_rid)}\\n")
                out_f.writelines(_buf)
        if args.contigs and args.contigs.exists():
            for _rid, _buf in _read_fasta_records(args.contigs):
                out_f.write(f">{_uniq(_rid)}\\n")
                out_f.writelines(_buf)"""
if old in src:
    src = src.replace(old, new, 1)
    applied.append("合并时保证唯一 ID (仅冲突时加 sample| 前缀)")
elif "_uniq(" in src:
    applied.append("已修复")
else:
    sys.exit("锚点未匹配")

if src != orig:
    if not BAK.exists():
        shutil.copy2(TREE, BAK)
        print(f"备份: {BAK}")
    TREE.write_text(src, encoding="utf-8")
    print("[acvirus_tree_pro.py]")
    for a in applied:
        print("  -", a)
