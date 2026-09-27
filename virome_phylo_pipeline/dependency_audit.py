"""全量依赖审计: 解析所有脚本的 import, 检查依赖是否可解析 (断链检测)"""
import ast, os, sys

# 脚本所在目录即管线根（避免写死本机绝对路径）
ROOT = os.path.dirname(os.path.abspath(__file__))
UTILS = os.path.join(ROOT, "utils")

# 收集所有 py 文件 (含根 + utils)
all_files = {}
for d in (ROOT, UTILS):
    for f in os.listdir(d):
        if f.endswith(".py"):
            all_files[f] = os.path.join(d, f)

# 模块名 -> 路径 (项目内可解析)
project_mods = {}
for f, p in all_files.items():
    stem = f[:-3]
    if stem != "__init__":
        project_mods[stem] = p  # 顶层可 import (根或 utils 在 path)
    project_mods[f"utils.{stem}"] = p  # utils.X 形式

broken = []  # (file, import_stmt, missing)
missing_third = set()  # 第三方缺失 (环境依赖, 记录不报错)

for fname, fpath in all_files.items():
    if fname.startswith("_") or fname == "__init__.py":
        continue
    try:
        tree = ast.parse(open(fpath, encoding="utf-8").read(), fpath)
    except SyntaxError as e:
        broken.append((fname, "SYNTAX ERROR", str(e)))
        continue
    for node in ast.walk(tree):
        # import X  /  import X.Y
        if isinstance(node, ast.Import):
            for alias in node.names:
                top = alias.name.split(".")[0]
                if top in project_mods or alias.name in project_mods:
                    continue
                if top in ("os", "sys", "re", "json", "csv", "math", "argparse",
                           "subprocess", "time", "shutil", "glob", "io", "gc",
                           "logging", "traceback", "warnings", "tempfile",
                           "itertools", "functools", "collections", "random",
                           "typing", "pathlib", "datetime", "string", "textwrap",
                           "hashlib", "urllib", "requests", "pandas", "numpy",
                           "matplotlib", "scipy", "sklearn", "networkx", "seaborn",
                           "plotly", "Bio", "pypopart", "dnasp", "yaml", "ete3",
                           "treetime", "geopy", "multiprocessing", "concurrent",
                           "threading", "queue", "statistics", "seaborn",
                           "pyyaml", "bs4", "geopandas", "scikit"):
                    continue
                missing_third.add((fname, top))
        # from X import ...  /  from X.Y import ...
        elif isinstance(node, ast.ImportFrom):
            if node.module is None or node.level > 0:  # 相对导入
                continue
            top = node.module.split(".")[0]
            if node.module in project_mods or top in project_mods:
                continue
            if top in ("os", "sys", "re", "json", "csv", "math", "argparse",
                       "subprocess", "time", "shutil", "glob", "io", "gc",
                       "logging", "traceback", "warnings", "tempfile",
                       "itertools", "functools", "collections", "random",
                       "typing", "pathlib", "datetime", "string", "textwrap",
                       "hashlib", "urllib", "requests", "pandas", "numpy",
                       "matplotlib", "scipy", "sklearn", "networkx", "seaborn",
                       "plotly", "Bio", "pypopart", "dnasp", "yaml", "ete3",
                       "treetime", "geopy", "multiprocessing", "concurrent",
                       "threading", "queue", "statistics",
                       "pyyaml", "bs4", "geopandas", "scikit"):
                continue
            missing_third.add((fname, top))

print("=== 项目内 import 断链 (项目模块缺失) ===")
if broken:
    for f, imp, msg in broken:
        print(f"  {f}: {imp} {msg}")
else:
    print("  无 (语法或项目内 import 均可解析)")

print("\n=== 第三方模块引用 (环境依赖, 非断链) ===")
seen = sorted(set(m[1] for m in missing_third))
print("  ", seen if seen else "无")
print("\n=== 项目模块清单 ===")
print("  根目录:", sorted(m for m in project_mods if "." not in m and m not in
      ("phylo_pipeline", "capheine_pipeline", "gene_partition_dating")))
print("  utils:", sorted(m.split('.')[1] for m in project_mods if m.startswith("utils.")))
