"""Pytest 路径引导: 让 tests/ 下的用例可以直接 import 管线目录中的模块。"""
import sys
from pathlib import Path

PIPELINE_DIR = Path(__file__).resolve().parent.parent
if str(PIPELINE_DIR) not in sys.path:
    sys.path.insert(0, str(PIPELINE_DIR))
