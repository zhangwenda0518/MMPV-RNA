#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""管线脚本引用完整性检查:
扫描核心管线 .py 文件中以字符串形式引用的同仓库脚本 (SCRIPT_DIR/'x.py', Path(__file__).../'x.py',
PREPROC_DIR/'x.py', '../pipeline/x.py', subprocess 里的 x.py 等), 验证目标文件在磁盘上存在。

目录常量 (SCRIPT_DIR / PREPROC_DIR / script_dir 等) 会先按脚本内赋值解析成真实目录,
再用于拼接裸文件名引用, 避免把 `PREPROC_DIR / 'x.py'` 这类正确引用误报为断链。
"""
import re, sys
from pathlib import Path

R = Path(r"D:\桌面\延伸基因组\MMPV-RNA")
CORE_DIRS = ["virome_discovery_pipeline", "virome_analysis_pipeline",
             "public_metadata_pipeline", "virome_submission_pipeline",
             "virome_phylo_pipeline", "endogenous_virus_pipeline",
             "metadata_gui", "submission_gui", "biosoft"]

# 排除目录 (统一成 / 分隔再比较, 否则 Windows 反斜杠路径永远匹配不上)
EXCLUDE = ("/archive", "__pycache__", ".bak", "_backup", "_scratch", "_archive",
           "server_checks", "/build/", "/dist/", "node_modules")

DIR_CONST_RE = re.compile(
    r"^([A-Za-z_][A-Za-z0-9_]*)\s*=\s*"
    r"(?:Path\(__file__\)(?P<parents>(?:\.parent)*)"
    r"|(?P<base>[A-Za-z_][A-Za-z0-9_]*)(?P<bparents>(?:\.parent)*))\s*"
    r"(?:/\s*['\"](?P<seg>[\w\-.]+)['\"])?",
    re.M)


def dir_constants(py, text):
    """解析以 __file__ 为锚的目录常量 → {常量名: 目录 Path}

    支持 SCRIPT_DIR = Path(__file__).resolve().parent
         PREPROC_DIR = SCRIPT_DIR.parent / "data_preprocessing_pipeline"
    """
    consts = {}
    pending = list(DIR_CONST_RE.finditer(text))
    for _ in range(len(pending) + 1):          # 反复解析, 处理常量之间的依赖
        progressed = False
        for m in pending:
            name, base, seg = m.group(1), m.group("base"), m.group("seg")
            if name in consts:
                continue
            if base is None:                    # Path(__file__) 锚点
                d = py.parent
                for _p in re.findall(r"\.parent", m.group("parents") or ""):
                    d = d.parent
            else:
                if base not in consts:
                    continue
                d = consts[base]
                for _p in re.findall(r"\.parent", m.group("bparents") or ""):
                    d = d.parent
            if seg:
                d = d / seg
            consts[name] = d
            progressed = True
        if not progressed:
            break
    return consts


missing = []
checked = 0
for d in CORE_DIRS:
    base = R / d
    if not base.exists():
        print(f"DIR MISSING: {d}"); missing.append((str(base), "directory"))
        continue
    for py in base.rglob("*.py"):
        s = str(py)
        if any(x in s.replace("\\", "/") for x in EXCLUDE):
            continue
        try:
            text = py.read_text(encoding="utf-8", errors="ignore")
        except Exception:
            continue
        consts = dir_constants(py, text)
        # 找出所有 "*.py" 字符串字面量
        for m in re.finditer(r"""['"]([\w\-./\\]+\.(?:py|R|sh|jar))['"]""", text):
            ref = m.group(1).replace("\\", "/")
            if ref.startswith(("biosoft/", "virome_", "public_", "metadata_", "submission_", "scripts/", "doc/")):
                cand = [R / ref]
            elif "/" in ref or ref.endswith(".py"):
                # 相对引用: 相对脚本目录 / 脚本目录的父级 / 仓库根, 以及显式拼接的目录常量
                cands = [py.parent / ref, py.parent.parent / ref, R / ref]
                for cname, cdir in consts.items():
                    if re.search(rf"{re.escape(cname)}\s*/\s*['\"]{re.escape(m.group(1))}['\"]", text):
                        cands.append(cdir / ref)
                cand = [c for c in cands if c.exists()]
                checked += 1
                if not cand and not any(c.exists() for c in cands):
                    # 只报告指向仓库内的引用
                    probe = (py.parent / ref).resolve()
                    if str(probe).startswith(str(R)):
                        missing.append((s.replace(str(R) + "\\", ""), ref))
                else:
                    continue
            else:
                continue
            if cand and not cand[0].exists():
                missing.append((s.replace(str(R) + "\\", ""), ref))

print(f"references checked: {checked}")
print(f"=== BROKEN REFERENCES ({len(missing)}) ===")
for src, ref in sorted(set(missing)):
    print(f"  {src}  ->  {ref}")
