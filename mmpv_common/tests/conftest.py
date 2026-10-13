# -*- coding: utf-8 -*-
"""pytest 共享配置: 把仓库根塞进 sys.path, 使 `from mmpv_common import ...` 可用。"""

import sys
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parent.parent.parent
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))
