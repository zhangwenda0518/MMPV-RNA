"""端到端验证: 用打过补丁的模块直接测 is_blacklisted + merge 列可用性。"""
import sys
import pandas as pd

RUN = "/home/zhangwenda/MMPV-RNA/virome_discovery_pipeline"
sys.path.insert(0, RUN)

# 直接 exec 模块里的函数定义 (绕过 main 的 argparse), 用 ast 抽函数
import ast
src = open(f"{RUN}/run_host_prediction.py").read()
tree = ast.parse(src)
funcs = [n for n in tree.body if isinstance(n, ast.FunctionDef)
         and n.name in ("is_trusted_level", "is_blacklisted")]
mod = ast.Module(body=funcs, type_ignores=[])
ns = {}
# 先注入依赖的集合
assigns = [n for n in tree.body if isinstance(n, ast.Assign)]
set_assigns = [n for n in assigns
               if any(isinstance(t, ast.Name) and
                      t.id in ("NON_PLANT_GENERA", "NON_PLANT_FAMILIES_FALLBACK",
                               "TRUSTED_LEVELS") for t in n.targets)]
exec(compile(ast.Module(body=set_assigns + funcs, type_ignores=[]), "<m>", "exec"), ns)

is_blacklisted = ns["is_blacklisted"]

# 真实案例测试
cases = [
    # (family, genus, level, 期望被否决)
    ("Partitiviridae", "Biavirus", "Family(via Family)", True),
    ("Partitiviridae", "Deltapartitivirus", "Family(via Family)", False),
    ("Mimiviridae", "Fabavirus", "Species(via Species)", False),
    ("Mimiviridae", "Nepovirus", "Genus(via Genus)", False),
    ("Mimiviridae", "Megavirus", "Family(via Family)", True),
    ("Tombusviridae", "Rimosavirus", "Family(via Family)", True),
    ("Potyviridae", "Zetanudivirus", "Family(via Family)", True),
    ("Caulimoviridae", "Sylvanvirus", "Family(via Family)", True),
    ("Partitiviridae", "NA", "Family(via Family)", False),   # 无属→科不在兜底→放行
    ("Mimiviridae", "NA", "Family(via Family)", True),        # 无属→科在兜底→否决
    ("Potyviridae", "Potyvirus", "Species(via Species)", False),
]
print(f"{'科':<18} {'属':<20} {'级别':<22} {'期望':<6} {'实际':<6} {'结果'}")
allok = True
for fam, gen, lvl, exp in cases:
    got = is_blacklisted(fam, gen, lvl)
    ok = "✓" if got == exp else "✗ FAIL"
    if got != exp: allok = False
    print(f"{fam:<18} {gen:<20} {lvl:<22} {str(exp):<6} {str(got):<6} {ok}")

print("\n" + ("全部通过 ✓" if allok else "存在失败 ✗"))
sys.exit(0 if allok else 1)
