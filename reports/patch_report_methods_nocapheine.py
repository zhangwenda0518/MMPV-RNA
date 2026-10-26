# -*- coding: utf-8 -*-
"""
给 generate_pipeline_report.py 打第 4 批补丁：把「capheine 残留」也做成条件化。

背景：当前已知病毒管道是 9 阶段（07_capheine 已移出到 virome_phylo_pipeline）。
前 3 批补丁已处理 导航/流程条/段编号/目录名，但 Methods 段与 stageFocus JS
仍是硬编码 10 阶段，无 capheine 时报告里会残留：
  - 脚本映射表第 7 项 "Selection Analysis / capheine_pipeline.py"
  - <details> "Selection Analysis (Stage 7)"
  - "through 10 automated stages"
  - stageFocus[7] = "该病毒位于Stage 7 Capheine区域..."

本补丁：
  A. 新增 _stage_methods_for(has_capheine)，无 capheine 时删第 7 项、8/9/10 前移一位
  B. _build_methods_section 增加 has_capheine 参数并按裁剪后的表渲染
  C. "through 10 automated stages" -> 动态阶段数；Selection 要点改为条件输出
  D. 调用点传入 HAS_CAPHEINE
  E. stageFocus 后追加 JS 重映射：无 capheine 时 delete 7 并把 8/9 前移

写盘放最后 + 每处锚点计数断言，失败即整体回滚（原子）。
"""
import os
import sys

# 本文件支持的环境变量：MMPV_SERVER_ROOT、MMPV_PAPER_ROOT（未设时保持原服务器字面量，路径逐字不变）
_SERVER = os.environ.get("MMPV_SERVER_ROOT", "/home/zhangwenda")
_PAPER  = os.environ.get("MMPV_PAPER_ROOT",  f"{_SERVER}/MMPV-paper")
_REPO   = f"{_SERVER}/MMPV-RNA"


TARGET = f"{_REPO}/virome_analysis_pipeline/generate_pipeline_report.py"

with open(TARGET, "r", encoding="utf-8") as f:
    src = f.read()

orig = src
edits = []


def rep(old, new, n_expect, tag):
    """带计数断言的替换"""
    global src
    n = src.count(old)
    if n != n_expect:
        print(f"FAIL[{tag}] 期望 {n_expect} 处锚点，实际 {n} 处 -> 整体回滚，不写盘")
        print("---- 锚点片段 ----")
        print(old[:400])
        sys.exit(2)
    src = src.replace(old, new, n_expect)
    edits.append(tag)


# ---------------- A. 新增裁剪helper（挂在 STAGE_METHODS 列表之后） ----------------
OLD_A = '''    (10, "Report Generation", "generate_pipeline_report.py",
     "This interactive HTML report was auto-generated from pipeline outputs, embedding "
     "all figures, tables, and result files. Includes AI-assisted interpretation via "
     "LLM integration (OpenAI/DeepSeek/Ollama compatible)."),
]
'''
NEW_A = '''    (10, "Report Generation", "generate_pipeline_report.py",
     "This interactive HTML report was auto-generated from pipeline outputs, embedding "
     "all figures, tables, and result files. Includes AI-assisted interpretation via "
     "LLM integration (OpenAI/DeepSeek/Ollama compatible)."),
]


def _stage_methods_for(has_capheine=True):
    """按当前管道实际阶段裁剪 STAGE_METHODS。

    07_capheine 已迁出已知病毒管道(-> virome_phylo_pipeline)，当前为 9 阶段。
    无 capheine 产物时删掉编号 7 的「Selection Analysis」，并把 8/9/10 前移一位，
    避免出现 6 -> 8 的跳号，同时保证「阶段数」与实际目录一致。
    """
    if has_capheine:
        return list(STAGE_METHODS)
    out = []
    for num, name, script, desc in STAGE_METHODS:
        if num == 7:
            continue
        if num > 7:
            num -= 1
        out.append((num, name, script, desc))
    return out
'''
rep(OLD_A, NEW_A, 1, "A:helper")

# ---------------- B. 函数签名 + 渲染改用裁剪后的表 ----------------
OLD_B = '''def _build_methods_section(viruses, overview_metrics, total_records):
    """Generate dynamic methods section from STAGE_METHODS metadata."""
    n_viruses = len(viruses)

    # Pipeline summary table
    rows = ""
    for num, name, script, _ in STAGE_METHODS:
        rows += f"<tr><td>{num}. {name}</td><td>{_esc(script)}</td></tr>"

    # Stage details (collapsible)
    details = ""
    for num, name, script, desc in STAGE_METHODS:
'''
NEW_B = '''def _build_methods_section(viruses, overview_metrics, total_records, has_capheine=True):
    """Generate dynamic methods section from STAGE_METHODS metadata."""
    n_viruses = len(viruses)
    methods = _stage_methods_for(has_capheine)

    # Pipeline summary table
    rows = ""
    for num, name, script, _ in methods:
        rows += f"<tr><td>{num}. {name}</td><td>{_esc(script)}</td></tr>"

    # Stage details (collapsible)
    details = ""
    for num, name, script, desc in methods:
'''
rep(OLD_B, NEW_B, 1, "B:signature")

# ---------------- C. 阶段数动态化 + Selection 要点条件化 ----------------
OLD_C2 = '''{overview_metrics.get("n_samples","?")} samples through 10 automated stages.</p>'''
NEW_C2 = '''{overview_metrics.get("n_samples","?")} samples through {len(methods)} automated stages.</p>'''
rep(OLD_C2, NEW_C2, 1, "C:stage-count")

OLD_C3 = '''            <li><strong>Selection</strong>: HyPhy codon-level testing (FEL/MEME/BUSTED/PRIME) for coding viruses</li>
'''
NEW_C3 = '''            {sel_bullet}
'''
rep(OLD_C3, NEW_C3, 1, "C:sel-bullet")

# 在 return f''' 之前定义 sel_bullet
OLD_C4 = '''    return f\'\'\'
    <h3 id="paper-reference">Pipeline Methods</h3>'''
NEW_C4 = '''    sel_bullet = ("<li><strong>Selection</strong>: HyPhy codon-level testing "
                  "(FEL/MEME/BUSTED/PRIME) for coding viruses</li>") if has_capheine else ""

    return f\'\'\'
    <h3 id="paper-reference">Pipeline Methods</h3>'''
rep(OLD_C4, NEW_C4, 1, "C:sel-bullet-def")

# ---------------- D. 调用点传入 HAS_CAPHEINE ----------------
OLD_D = '''        {_build_methods_section(viruses, overview_metrics, total_records)}'''
NEW_D = '''        {_build_methods_section(viruses, overview_metrics, total_records, HAS_CAPHEINE)}'''
rep(OLD_D, NEW_D, 1, "D:call-site")

# ---------------- E0. 暴露 JS 布尔量 ----------------
OLD_E0 = '''    HAS_CAPHEINE = (out / "07_capheine").is_dir()'''
NEW_E0 = '''    HAS_CAPHEINE = (out / "07_capheine").is_dir()
    _has_cap_js = "true" if HAS_CAPHEINE else "false"'''
rep(OLD_E0, NEW_E0, 1, "E0:js-bool")

# ---------------- E. stageFocus 重映射 ----------------
OLD_E = """\u6709\u5173\u3002'}};\n"""   # 有关。'}};
NEW_E = (
    "\u6709\u5173\u3002'}};\n"
    "    if(!{_has_cap_js}){{delete stageFocus[7];var _sr={{}};"
    "for(var _k in stageFocus){{var _n=+_k;_sr[_n>7?_n-1:_n]=stageFocus[_k];}}"
    "stageFocus=_sr;}}\n"
)
rep(OLD_E, NEW_E, 1, "E:stageFocus-remap")

# ---------------- 写盘（放最后） ----------------
tmp = TARGET + ".tmp_patch4"
with open(tmp, "w", encoding="utf-8") as f:
    f.write(src)
os.replace(tmp, TARGET)

print("PATCHED OK ->", TARGET)
print("edits:", ", ".join(edits))
print("size:", len(orig), "->", len(src), f"(delta {len(src)-len(orig)})")
