#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""batch_draw_pymol.py — 兼容垫片 (2026-09-15, 审查 P3)

背景
----
历史上根目录与本包 `utils/` 下各有一份**内容不同**的同名实现:

  · 根目录版 (旧, 22,459 B): `__main__` 守卫正确, 但正选择位点解析读的是
    **不存在的列** `meme_marker` (真实列是 `meme_qval`/`meme_pval`),
    且缺 `prime_marker` 分支 → MEME 检出位点静默丢失、prime 标记位点被忽略。
  · `utils/` 版 (新, 21,589 B): marker 解析已修, 但重执行块与末尾 `main()`
    都没有 `__main__` 守卫 → 被 import 时会换掉/杀掉当前进程。

`run_postprocess("batch_draw_pymol", ...)` 因为 `utils/` 下有同名文件, 实际解析到
`utils.batch_draw_pymol` —— 也就是说**跑的是新逻辑**, 而根目录那份是死代码但会
误导人 (改错文件不生效)。

本文件现统一为垫片, 单一实现保留在 `utils/batch_draw_pymol.py`。
两个入口现在完全等价:

    python -m batch_draw_pymol --input x.csv --output out --pdb_dir pdbs
    python -m utils.batch_draw_pymol --input x.csv --output out --pdb_dir pdbs

原根目录实现已备份到 _fix_backup_20260915/batch_draw_pymol.root.py。
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from utils.batch_draw_pymol import main  # noqa: E402


if __name__ == '__main__':
    main()
