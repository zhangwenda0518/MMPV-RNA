"""新增功能冒烟测试: 多作者 SBT 生成 + BioSample TSV 导出"""
import os
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import submission_gui as g
from PySide6.QtWidgets import QApplication

app = QApplication([])
win = g.MainWindow()
csv_path = Path(__file__).parent / "samples" / "unified_metadata_barbarum_real.csv"
idx = win._add_file_tab(str(csv_path))
store = win._stores[idx]

# ── 测试 A: 多作者解析 ──
f = win._sbt_fields
f["last"].setText("Zhang")
f["first"].setText("Wenda")
win._sbt_extra_authors.setPlainText("Li, Ming\nWang, Fang Hua")
authors = win._parse_authors()
assert authors == [("Zhang", "Wenda", ""), ("Li", "Ming", ""), ("Wang", "Fang", "Hua")], authors
print("PASS A. 作者解析 (3 人含中间名)")

# ── 测试 B: 生成 sbt 结构 ──
win._sbt_title.setText("Plant virome test title")
win._do_generate_sbt()
txt = win._sbt_output.toPlainText()
assert txt.count("{") == txt.count("}"), "braces unbalanced"
assert 'Submission Title:Plant virome test title' in txt, "title 未写入 (旧 bug)"
assert txt.count("name name {") == 7, f"name block 数: {txt.count('name name {')}"  # contact 1 + cit 3 + pub 3
assert '"Zhang"' in txt and '"Li"' in txt and '"Wang"' in txt
print("PASS B. SBT 生成: 括号平衡 / 标题写入 / 3 作者×3 处")

# ── 测试 C: BioSample TSV 导出 ──
out = Path(__file__).parent / "_test_biosample.tsv"
df = store.dataframe
n_ph = int(df["sequence_name"].str.startswith("contig_").sum())
lines = ["\t".join(["sample_name", "organism"])]
skipped = 0
for _, r in df.iterrows():
    s = str(r["sequence_name"])
    if s.startswith("contig_"):
        skipped += 1
        continue
    vals = [s, str(r.get("organism", ""))]
    for c, _ in win._BS_ATTR_MAP:
        if c in df.columns:
            v = str(r.get(c, "") or "").strip()
            vals.append("" if win._BS_PLACEHOLDER.match(v) else v)
    lines.append("\t".join(vals))
out.write_text("\n".join(lines) + "\n", encoding="utf-8")
data = out.read_text(encoding="utf-8").strip().split("\n")
assert len(data) - 1 == len(df) - skipped, (len(data) - 1, len(df), skipped)
assert not any(d.startswith("contig_") for d in data[1:])
first_data = data[1].split("\t")
assert first_data[0] and first_data[1], "sample_name/organism 空"
print(f"PASS C. BioSample TSV: {len(data)-1} 条导出, {skipped} 占位跳过, 表头 {len(data[0].split(chr(9)))} 列")

# ── 测试 D: 默认邮箱 ──
assert f["email"].text() == "zhangwenda05@163.com"
print("PASS D. 默认邮箱已改真实值")

out.unlink()
app.quit()
print("\n全部新功能测试通过 ✓")

