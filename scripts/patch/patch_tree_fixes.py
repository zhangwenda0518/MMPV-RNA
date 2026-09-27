#!/usr/bin/env python3
"""
修 9b 建树链路的 3 个问题:

P1（巨型科全量兜底 → MAFFT 4h 超时）
  _build_one 里 `use_fa = fam_keep if fam_keep.exists() else keep_fa`:
  某科在 KEEP 里 0 条 contig 时 fallback 到全量 class_KEEP.fasta (2365 条),
  MAFFT --auto 处理 ~2400 条必然跑满 4h 超时。这类科本就没有自研序列要定位。
  → 无 keep_<Fam>.fasta 时跳过该科, 记 rc=3 (区别于真失败)。

P2（trimal 删空序列后 rc=0 但无树 → 假成功）
  acvirus_tree_pro.py 的 n_after<3 分支 return (退出码 0), 无 treefile 产出,
  管线看到 rc=0 当成功, 静默缺树。
  → 改 sys.exit(3) 显式标记; 管线侧识别 rc=3 为跳过。

P3（重启后已完成科会全部重跑）
  建树循环无 skip 机制, 重跑浪费数小时。
  → 已存在 *.treefile 的科直接跳过。
"""
import shutil
import sys
from pathlib import Path

PIPE = Path("/home/zhangwenda/MMPV-RNA/virome_discovery_pipeline/virome_pipeline.py")
TREE = Path("/home/zhangwenda/MMPV-RNA/virome_discovery_pipeline/utils/acvirus_tree_pro.py")
BAK_PIPE = PIPE.with_suffix(".py.bak_treefix_20260910")
BAK_TREE = TREE.with_suffix(".py.bak_treefix_20260910")

for p in (PIPE, TREE):
    if not p.is_file():
        sys.exit(f"找不到: {p}")

# ══════════ virome_pipeline.py ══════════
src = PIPE.read_text(encoding="utf-8")
orig = src
applied = []

# P1 + P3: 合并改 _build_one 与调用
old = """                def _build_one(fam):
                    fam_out = tree_root / fam
                    try:
                        # 科内 KEEP 子集 (无则用全量 KEEP 兼底)
                        fam_keep = fam_out / f'keep_{fam}.fasta'
                        use_fa = fam_keep if fam_keep.exists() else keep_fa
                        r = _sp.run(['python3', str(acv_script),"""
new = """                def _build_one(fam):
                    fam_out = tree_root / fam
                    try:
                        # P3: 已有树则跳过 (重启时复用, 免数小时重跑)
                        if list(fam_out.glob('*.treefile')):
                            return fam, 4, 'tree exists'
                        # P1: 无本科 KEEP 子集则跳过 (无自研 contig; 全量兼底会让 MAFFT 超时)
                        fam_keep = fam_out / f'keep_{fam}.fasta'
                        if not fam_keep.exists():
                            return fam, 3, 'no family-specific KEEP contig'
                        use_fa = fam_keep
                        r = _sp.run(['python3', str(acv_script),"""
if old in src:
    src = src.replace(old, new, 1)
    applied.append("P1+P3: 无 keep 跳过 / 已有树跳过")
elif "P3: 已有树则跳过" in src:
    applied.append("P1+P3: 已修复")
else:
    sys.exit("P1+P3 锚点未匹配")

# P1b/P2b: 管线侧日志区分
old2 = """                    for fu in as_completed(futs):
                        fam, rc, err = fu.result()
                        self.log.info('  [9b] tree %s rc=%d %s', fam, rc, err if rc else '')"""
new2 = """                    for fu in as_completed(futs):
                        fam, rc, err = fu.result()
                        if rc == 0:
                            self.log.info('  [9b] tree %s ✓', fam)
                        elif rc in (3, 4):
                            self.log.info('  [9b] tree %s ⊘ %s', fam, err)
                        else:
                            self.log.warning('  [9b] tree %s ✗ rc=%d %s', fam, rc, err)"""
if old2 in src:
    src = src.replace(old2, new2, 1)
    applied.append("日志区分 rc=0/3/4/其它")
elif "tree %s ✓" in src:
    applied.append("日志: 已修复")
else:
    sys.exit("日志锚点未匹配")

if src != orig:
    if not BAK_PIPE.exists():
        shutil.copy2(PIPE, BAK_PIPE)
        print(f"备份: {BAK_PIPE}")
    PIPE.write_text(src, encoding="utf-8")
    print("[virome_pipeline.py]")
    for a in applied:
        print("  -", a)

# ══════════ acvirus_tree_pro.py ══════════
tsrc = TREE.read_text(encoding="utf-8")
torig = tsrc
tapplied = []

old3 = """    n_after=_count_phylip(trimmed)
    if n_after<3:
        print(f"[SKIP] {outdir.name}: only {n_after} sequences after trimming, skip tree")
        return"""
new3 = """    n_after=_count_phylip(trimmed)
    if n_after<3:
        # 退出码 3 = 跳过 (trimal 删空后无产物; 需与真失败区分)
        print(f"[SKIP] {outdir.name}: only {n_after} sequences after trimming, skip tree")
        sys.exit(3)"""
if old3 in tsrc:
    tsrc = tsrc.replace(old3, new3, 1)
    tapplied.append("trimal 删空分支改 sys.exit(3)")
elif "sys.exit(3)" in tsrc:
    tapplied.append("已修复")
else:
    sys.exit("acvirus P2 锚点未匹配")

if tsrc != torig:
    if not BAK_TREE.exists():
        shutil.copy2(TREE, BAK_TREE)
        print(f"备份: {BAK_TREE}")
    TREE.write_text(tsrc, encoding="utf-8")
    print("[acvirus_tree_pro.py]")
    for a in tapplied:
        print("  -", a)

print("\n完成")
