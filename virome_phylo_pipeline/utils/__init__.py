"""virome_phylo_pipeline.utils — 依赖模块包

让 utils/ 内模块之间可以 `from utils.X import ...` 互相引用，
且 utils/ 内的独立库（如 dnasp.py）可直接 `import dnasp`：
任何 utils 模块被导入时，把 pipeline 根目录与 utils 目录都加入 sys.path。
"""
import os
import sys

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_UTIL = os.path.dirname(os.path.abspath(__file__))
for _p in (_ROOT, _UTIL):
    if _p not in sys.path:
        sys.path.insert(0, _p)
