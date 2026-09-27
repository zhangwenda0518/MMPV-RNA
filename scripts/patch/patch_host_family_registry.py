#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""纯追加补丁: 扩容 run_host_prediction.py 的 NON_PLANT_FAMILIES_FALLBACK。

形态: 原定义行一字不动, 紧随其后追加 `NON_PLANT_FAMILIES_FALLBACK += [...]`。
      回滚 = 删掉追加块 (或还原 .bak_famregistry_20260915)。
幂等: 检测到追加标记即跳过, 重复执行安全。

用法:
    python3 patch_host_family_registry.py            # dry-run, 只报告
    python3 patch_host_family_registry.py --apply    # 实际写入 (自动备份)
    python3 patch_host_family_registry.py --rollback # 从备份还原
"""
import argparse
import ast
import hashlib
import os
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
TARGET = os.environ.get(
    "HOST_PRED_TARGET",
    "/home/zhangwenda/MMPV-RNA/virome_discovery_pipeline/run_host_prediction.py")
PAYLOAD = os.path.join(HERE, "_family_registry_append.txt")
BAK = TARGET + ".bak_famregistry_20260915"
MARKER = "NON_PLANT_FAMILIES_FALLBACK += ["
DEF_MARK = "NON_PLANT_FAMILIES_FALLBACK = ["


def md5(path):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def newline_of(text):
    return "\r\n" if "\r\n" in text else "\n"


def parse_list(src, name):
    """按语句顺序累加取值: 兼容 `X = [...]` 与 `X += [...]` 两种形态。

    返回 None 表示文件里没有该名字的定义。
    """
    tree = ast.parse(src)
    found = False
    out = []
    for node in tree.body:
        if isinstance(node, ast.Assign) and any(
                isinstance(t, ast.Name) and t.id == name for t in node.targets):
            out = list(ast.literal_eval(node.value))
            found = True
        elif isinstance(node, ast.AugAssign) and isinstance(node.op, ast.Add) \
                and isinstance(node.target, ast.Name) and node.target.id == name:
            out = out + list(ast.literal_eval(node.value))
            found = True
    return out if found else None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--apply", action="store_true")
    ap.add_argument("--rollback", action="store_true")
    args = ap.parse_args()

    if args.rollback:
        if not os.path.exists(BAK):
            print("FAIL: 备份不存在 %s" % BAK)
            return 1
        shutil.copy2(BAK, TARGET)
        print("已回滚: %s <- %s (md5 %s)" % (TARGET, BAK, md5(TARGET)))
        print("当前条数:", len(parse_list(open(TARGET, encoding="utf-8").read(),
                                      "NON_PLANT_FAMILIES_FALLBACK") or []))
        return 0

    src = open(TARGET, encoding="utf-8", newline="").read()
    nl = newline_of(src)
    before = parse_list(src, "NON_PLANT_FAMILIES_FALLBACK")
    if before is None:
        print("FAIL: 未在目标文件中找到 NON_PLANT_FAMILIES_FALLBACK 定义")
        return 1

    if MARKER in src:
        print("SKIP: 追加块已存在 (幂等), 当前 %d 条" % len(before))
        return 0

    payload = open(PAYLOAD, encoding="utf-8").read()
    ns = {"NON_PLANT_FAMILIES_FALLBACK": []}
    exec(payload, ns)
    added = ns["NON_PLANT_FAMILIES_FALLBACK"]
    dup = sorted(set(added) & set(before))
    print("现有 %d 条; 追加 %d 条; 与现有重复 %d 条 %s" % (len(before), len(added), len(dup), dup[:8]))

    lines = src.split(nl)
    idx = [i for i, l in enumerate(lines) if l.startswith(DEF_MARK)]
    if len(idx) != 1:
        print("FAIL: 期望 1 处定义行, 实际 %d" % len(idx))
        return 1
    i = idx[0]
    new_lines = lines[:i + 1] + payload.rstrip(nl).split(nl) + lines[i + 1:]
    out = nl.join(new_lines)

    after = parse_list(out, "NON_PLANT_FAMILIES_FALLBACK")
    if after is None:
        print("FAIL: 改写后 AST 解析失败")
        return 1
    print("改写后 %d 条 (净增 %d)" % (len(after), len(after) - len(before)))
    if set(after) != set(before) | set(added):
        print("FAIL: 追加后集合与预期不符 (期望 %d 个唯一值, 实得 %d)"
              % (len(set(before) | set(added)), len(set(after))))
        return 1
    if lines[i] != new_lines[i]:
        print("FAIL: 原定义行被改动")
        return 1
    print("原定义行逐字节未变: OK")

    if not args.apply:
        print("DRY-RUN: 未写入。加 --apply 生效。")
        return 0

    if not os.path.exists(BAK):
        shutil.copy2(TARGET, BAK)
        print("备份: %s (md5 %s)" % (BAK, md5(BAK)))
    open(TARGET, "w", encoding="utf-8", newline="").write(out)
    print("写入: %s (md5 %s -> %s)" % (TARGET, md5(BAK), md5(TARGET)))

    r = subprocess.run([sys.executable, "-c",
                        "import importlib.util,sys;"
                        "spec=importlib.util.spec_from_file_location('rhp',r'%s');"
                        "m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);"
                        "print('import OK; 科兜底',len(m.NON_PLANT_FAMILIES_FALLBACK),"
                        "'属黑名单',len(m.NON_PLANT_GENERA));"
                        "print('样例', m.is_blacklisted('Mimiviridae','','Family'),"
                        "m.is_blacklisted('Metaviridae','','Family'),"
                        "m.is_blacklisted('Potyviridae','','Family'),"
                        "m.is_blacklisted('Mimiviridae','Potyvirus','Genus'))" % TARGET],
                       capture_output=True, text=True)
    print("[verify rc=%d] %s%s" % (r.returncode, r.stdout.strip(), r.stderr.strip()[-300:]))
    return 0 if r.returncode == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
