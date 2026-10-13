#!/usr/bin/env python3
"""
tbl2gb.py — 顶层转发薄壳 (实现唯一维护在 utils/tbl2gb.py)
==========================================================
本文件曾是 utils/tbl2gb.py 的旧副本 (未含 P1-20 负链坐标修复), 无调用方但可被
import — 两份实现并存会漂移, 故改为转发壳: 所有函数与 CLI 行为等价于
utils/tbl2gb.py (含负链 FeatureLocation(min-1, max+1) 修复)。

用法 (与 utils/tbl2gb.py 完全一致):
  python tbl2gb.py \
    --tbl analyze_hypothetical/featuretable_updated.tbl \
    --fasta suvtk.features_output/reoriented_nucleotide_sequences.fna \
    -o virus-annotations/

注意: 用 importlib 按文件路径加载 utils/tbl2gb.py, 不做 `import utils.*`
      (本仓多个管线各有同名 `utils` 包, 按 sys.path 解析可能命中他人实现)。
"""

import importlib.util
import os
import sys

_IMPL_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "utils", "tbl2gb.py")


def _load_impl():
    spec = importlib.util.spec_from_file_location("tbl2gb_impl", _IMPL_PATH)
    if spec is None or spec.loader is None:
        raise ImportError(f"tbl2gb.py: 无法加载实现文件 {_IMPL_PATH}")
    mod = importlib.util.module_from_spec(spec)
    sys.modules["tbl2gb_impl"] = mod   # 供 dataclass/pickle 场景引用
    spec.loader.exec_module(mod)
    return mod


_impl = _load_impl()

# re-export 全部公开符号 (parse_feature_table / build_gb_record / load_taxonomy / main ...)
for _name in dir(_impl):
    if not _name.startswith("_"):
        globals()[_name] = getattr(_impl, _name)

if __name__ == "__main__":
    sys.exit(_impl.main())
