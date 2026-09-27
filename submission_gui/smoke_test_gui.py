#!/usr/bin/env python3
# smoke_test_gui.py — submission_gui 导入/编辑/导出功能冒烟测试 (offscreen)
import os, sys
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pandas as pd

# 只测数据层 + 必要的 Qt model (不弹窗口)
from submission_gui import SubmissionStore, is_placeholder

CSV = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                   "samples", "unified_metadata_barbarum_real.csv")
OUT_CSV = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       "samples", "_test_roundtrip_out.csv")
OUT_XLSX = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "samples", "_test_roundtrip_out.xlsx")

fails = []
def check(name, cond):
    print(("PASS " if cond else "FAIL ") + name)
    if not cond:
        fails.append(name)

s = SubmissionStore()
ok = s.load(CSV)
check("1. 导入真实 barbarum CSV (75条)", ok and s.row_count == 75)

st = s.stats()
print(f"   统计: {st['total']} 条 | {st['filled']}/{st['cells']} 已填 ({st['pct']:.1f}%)")
check("2. 统计功能 (填充率>20%)", st["pct"] > 20)

issues = s.validate()
n_missing_cells = sum(i["missing_count"] for i in issues if i["missing_count"] > 0)
print(f"   验证: {len(issues)} 个必填列存在占位, 共 {n_missing_cells} 个占位单元格")
for i in issues:
    print(f"     [{i['missing_count']}] {i['column']}  e.g. {i['examples'][:2]}")
check("3. 验证能识别占位符 (bioproject=PRJNAXXXXXX)", len(issues) > 0)

# 编辑单格
col_bp = s.columns.index("bioproject")
s.set_cell(0, col_bp, "PRJNA999999")
check("4. 双击编辑单格", s.get_cell(0, col_bp) == "PRJNA999999")

# 批量填充占位符
n1 = s.batch_fill_placeholder("biosample", "SAMN11111111")
real_before = sum(1 for v in s._df["biosample"] if v == "SAMN11111111")
check(f"5. 批量填充 biosample 占位 ({n1} 格)", n1 > 0 and n1 == real_before)
# 真实值不被覆盖 (quick fill 只填占位)
s.batch_fill_placeholder("bioproject", "PRJNA888888")
untouched = sum(1 for v in s._df["bioproject"] if v == "PRJNA999999")
check("6. 已编辑值不被批量填充覆盖", untouched == 1)

# 替换式批量
n2 = s.batch_replace("bioproject", "PRJNA888888", "PRJNA777777")
check(f"7. 批量替换 old→new ({n2} 格)", n2 == n1 if n1 else True)

# 保存 round-trip
ok = s.save(OUT_CSV)
df_out = pd.read_csv(OUT_CSV, dtype=str, keep_default_na=False) if ok else None
check("8. 另存 CSV 后重读一致", ok and df_out is not None and
      len(df_out) == 75 and list(df_out.columns) == s.columns)

# Excel 导出再导入
ok = s.export_excel(OUT_XLSX)
df_x = pd.read_excel(OUT_XLSX, dtype=str, keep_default_na=False) if ok else None
check("9. 导出 Excel 并可重读", ok and df_x is not None and len(df_x) == 75)

# 占位高亮逻辑抽查: 对照旧 schema 检查新列名
from submission_gui import REQUIRED_COLS
missing_cols = [c for c in REQUIRED_COLS if c not in s.columns]
check("10. 新 CSV 覆盖全部必填列", not missing_cols)
if missing_cols:
    print("   缺失必填列:", missing_cols)

print()
if fails:
    print("结果:", len(fails), "项失败:", fails); sys.exit(1)
print("结果: 全部通过 ✓")
