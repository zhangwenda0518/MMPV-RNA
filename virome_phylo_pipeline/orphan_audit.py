"""孤儿审计 v4 (AST 级) — 修复 v3 盲区, 引用判定升级

v3 (正则) 的盲区:
  函数内缩进导入 `from utils import X as _x` —— v3 的 `from utils.X import` /
  `import utils.X` 两条正则只认点分形式, 漏掉包导入形式 → seq_clean /
  metadata_channels 被误报孤儿 (phylo_pipeline.py:1043/2825 实际在用)。
  v3 已覆盖的形式 (v4 保留): `import X`、`from X import`、`from utils.X import`、
  `"-m", "X"` 列表调用、`X.py` 字面量、`run_postprocess("X", ...)`。

v4 引用判据 (并集):
  · ast.Import / ast.ImportFrom 全形态 (含函数内缩进导入、`from utils import X`)
  · ast.List 中 `..., "-m", "<mod>", ...` 相邻元素对
  · ast.Call func 名为 run_postprocess 的首个字面量参数
  · 全文 `<stem>.py` 字面量 (含注释/文档字符串, 保守口径)
输出分区: 管线可达 / 预期孤儿(文档登记的手工CLI) / 真孤儿。
"""
import ast, os, re

# 脚本所在目录即管线根（避免写死本机绝对路径）
ROOT = os.path.dirname(os.path.abspath(__file__))

MAIN_ENTRIES = ("phylo_pipeline.py", "capheine_pipeline.py", "gene_partition_dating.py")
QA_EXEMPT = ("dependency_audit.py", "orphan_audit.py", "patch_pypopart_haplotype.py")
# 文档登记过的独立 CLI / 手工工具 (RUN_GUIDE.md:356, OPTIMIZATION_20260916.md:312, phylo_pipeline.py:29/3270)
DOCUMENTED_MANUAL = {
    "multigene_tools.py",
    "codon_usage/codon_usage.py", "codon_usage/codon_usage_plots.py",
    "batch_draw_pymol.py", "utils/batch_draw_pymol.py",
    "utils/glm_predictors.py", "utils/import_export.py",
    "utils/pub_plots.py", "utils/splitstree_enrich.py",
    "utils/splitstree_viz.py", "utils/tiger_sites.py",
}

# ── 收集活代码 (根 / utils / codon_usage; 排除 _*、*.bak*、__init__) ──
files = {}
for sub in ("", "utils", "codon_usage"):
    d = os.path.join(ROOT, sub) if sub else ROOT
    if not os.path.isdir(d):
        continue
    for f in sorted(os.listdir(d)):
        if not f.endswith(".py") or f == "__init__.py":
            continue
        if f.startswith("_") or ".bak" in f:
            continue
        files[(sub + "/" + f) if sub else f] = os.path.join(d, f)

stems = {}
for key in files:
    stems.setdefault(os.path.basename(key)[:-3], set()).add(key)

def local_modules(ref):
    """import/-m/run_postprocess 引用 → 本仓库文件 key

    裸 stem (run_postprocess 的实参、`-m` 的模块名) 运行时可能命中根目录
    或 utils/ (phylo_pipeline.run_postprocess 会自动加 utils. 前缀),
    因此凡同名一律算引用, 不按目录过滤。"""
    out = set()
    parts = ref.split(".")
    top = parts[0]
    if top in ("utils", "codon_usage"):
        if len(parts) >= 2:
            for k in stems.get(parts[1], ()):
                if k.startswith(top + "/"):
                    out.add(k)
        return out
    return set(stems.get(top, ()))

def refs_of(key):
    src = open(files[key], encoding="utf-8", errors="replace").read()
    refs = set()
    try:
        tree = ast.parse(src)
    except SyntaxError as e:
        print(f"  [SYNTAX ERR] {key}: {e}")
        return refs
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for a in node.names:
                refs |= local_modules(a.name)
        elif isinstance(node, ast.ImportFrom):
            if node.module:
                refs |= local_modules(node.module)
            for a in node.names:
                refs |= local_modules(f"{node.module}.{a.name}" if node.module else a.name)
        elif isinstance(node, ast.List):
            elts = node.elts
            for i, elt in enumerate(elts[:-1]):
                if (isinstance(elt, ast.Constant) and elt.value == "-m"
                        and isinstance(elts[i + 1], ast.Constant)
                        and isinstance(elts[i + 1].value, str)):
                    refs |= local_modules(elts[i + 1].value)
        elif isinstance(node, ast.Call):
            fn = node.func
            if isinstance(fn, ast.Name) and fn.id == "run_postprocess" and node.args:
                a0 = node.args[0]
                if isinstance(a0, ast.Constant) and isinstance(a0.value, str):
                    refs |= local_modules(a0.value)
    for stem in stems:  # 保守口径: `<stem>.py` 任意出现即算引用
        if re.search(rf"\b{re.escape(stem)}\.py\b", src):
            refs |= stems[stem]
    return refs

reached, stack = set(MAIN_ENTRIES), list(MAIN_ENTRIES)
while stack:
    n = stack.pop()
    for r in refs_of(n) - reached:
        reached.add(r)
        stack.append(r)

candidates = [k for k in sorted(files)
              if k not in MAIN_ENTRIES and k not in QA_EXEMPT]
orphans = [k for k in candidates if k not in reached]
manual = [k for k in orphans if k in DOCUMENTED_MANUAL]
true_orphans = [k for k in orphans if k not in DOCUMENTED_MANUAL]

print(f"活代码: {len(files)} 个 | 管线可达: {len(reached)} | "
      f"预期孤儿(手工CLI): {len(manual)} | 真孤儿: {len(true_orphans)}")

if manual:
    print(f"\n=== 预期孤儿 — 文档登记的独立 CLI/手工工具 ({len(manual)}) ===")
    for k in manual:
        print(f"  {k}")
if true_orphans:
    print(f"\n=== 真孤儿 ({len(true_orphans)}) ===")
    for k in true_orphans:
        stem = os.path.basename(k)[:-3]
        mention = [k2 for k2 in files if k2 != k
                   and re.search(rf"\b{re.escape(stem)}\b",
                                 open(files[k2], encoding="utf-8", errors="replace").read())]
        has_main = bool(re.search(r"if __name__\s*==\s*['\"]__main__['\"]",
                                  open(files[k], encoding="utf-8", errors="replace").read()))
        loc = "根" if "/" not in k else k.split("/")[0]
        print(f"  {os.path.basename(k):28s} [{loc}] {'有main' if has_main else '纯库'} "
              f"被提及: {', '.join(mention[:4]) if mention else '—'}")
if not true_orphans:
    print("\n✅ 无真孤儿 — 所有未接入主流程的脚本均有文档登记 (见文件头 DOCUMENTED_MANUAL)")
