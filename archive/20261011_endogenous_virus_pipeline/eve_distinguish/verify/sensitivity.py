#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
sensitivity.py — 阈值敏感度扫描（谁在决定结论？）
=================================================
问题: 这条流水线的每个阈值都是"选"出来的常数, 从来没人量过它们各自能移动多少条判定。
`credit_components` 的 COMP_AA_MIN、`verdict_loci` 的 bitscore 门槛、s1 的终止切点 ——
如果某个阈值在很宽的范围内扫来扫去都几乎不动结论, 那它就不是一个"判别参数", 而是一份
装饰; 反过来, 结论对某个阈值极度敏感, 那它必须被校准过才敢用。

做法: 对每个 (脚本, 常数, 取值) 组合, 把常数**追加覆盖**到该脚本源码末尾 exec 成独立
模块（不碰磁盘上的源文件）, 用真 argv 调它的 main(), 重跑该 stage 及其下游, 然后统计
verdict 分布与"相对基线的变更条数"。

为什么可以追加覆盖: 这些常数都在函数体内被读取, 调用时查模块全局 —— 追加在末尾的赋值
就生效。唯一例外是派生量 (COMP_NT_MIN = COMP_AA_MIN*3), 扫 COMP_AA_MIN 时一并重算。

用法:
  python3 sensitivity.py --run-dir REAL --evidence E.tsv --work WORK [--only PARAM]

产物: WORK/sens.json + 屏幕表格
"""
import argparse
import collections
import json
import os
import re
import sys
import types

HERE = os.path.dirname(os.path.abspath(__file__))
MOD = os.path.join(HERE, "..")

WORK_RUN = WORK_EVID = WORK_CACHE = ""   # 由 main() 填, 供 run_probe 用
STAGES = {  # 脚本 -> 它的 argv 模板 (用 {out} 占位输出文件)
    "s1": ("s1_decay_scan.py", lambda q, d, out: [q, os.path.join(d, "all_hits.tsv"), out]),
    "s2": ("s2_domain_scan.py",
           lambda q, d, out: [q, os.path.join(d, "panel_hits.tsv"),
                              os.path.join(d, "baits_hits.tsv"), out]),
}

# 扫描计划: 参数 -> (脚本, [取值]). 取值含基线值, 便于同一张表里看趋势.
PLAN = {
    # s2 组件记功门槛: 一条比对要多少 aa 才证明一个组件
    "COMP_AA_MIN": ("s2", [50, 75, 100, 150, 200, 300]),
    # s2 HSP 入场门槛
    "MIN_ALN": ("s2", [30, 50, 75, 100, 150]),
    "MIN_PID": ("s2", [20, 25, 30, 35, 45]),
    "MIN_BITS": ("s2", [30, 40, 60, 80]),
    # s2 前病毒规模定义
    "PROVIRUS_NCOMP": ("s2", [3, 4, 5]),
    "PROVIRUS_SPAN": ("s2", [500, 1000, 1500, 2500, 4000]),
    # s1 退化切点
    "MIN_HEAD_STOPS": ("s1", [1, 2, 3, 4]),
    "MIN_DECAY_STOPS": ("s1", [2, 3, 4, 6, 8]),
    "CLEAN_ENRICH": ("s1", [0.3, 0.5, 0.7, 1.0]),
    "CLEAN_MAX_STOPS": ("s1", [2, 5, 10, 20]),
    "MIN_SPAN": ("s1", [60, 90, 150, 300]),
    # s3 上游病毒信号
    "AA_VIRAL_PID": ("s3", [80.0, 90.0, 95.0, 98.0, 100.0]),
}
DERIVED = {"COMP_AA_MIN": ["COMP_NT_MIN"]}   # 跟着重算的派生量



# ─────────── 干预生效性: 看下游效应, 不看变量本身 ───────────
# 教训: 我原来的自检断言的是"模块全局被改成了新值"—— 它**通过**了, 而结论依然错, 因为
# `credit_components(hsps, comp_aa_min=COMP_AA_MIN)` 的默认参数在 def 时求值, 改全局到不了
# 那条代码路径。所以必须检查**一个必须随之移动的下游量**: 参数动而下游量一动不动, 就
# 只能报"仪器可疑", 不能报"该阈值无影响"。
def _diag_s2(path, kind):
    """s2 产物的下游诊断量."""
    rows = _tsv(path)
    if kind == "ncomp_dist":
        c = collections.Counter(int(r.get("best_locus_ncomp") or 0) for r in rows)
        return tuple(sorted(c.items()))
    if kind == "provirus_n":
        return sum(1 for r in rows if r.get("provirus_scale") == "TRUE")
    if kind == "with_comps":
        return sum(1 for r in rows if (r.get("components") or "-") != "-")
    if kind == "credit_hsps":
        return sum(int(r.get("credit_hsps") or 0) for r in rows)
    if kind == "te_n":
        # MIN_BITS 只作用在 baits 路径(转座子域/panel 同源), 与 panel 组件无关 ——
        # 拿 with_comps 当它的诊断量, 会把它误判成"仪器可疑"。
        return sum(1 for r in rows if (r.get("te_flag") or "").strip().upper() == "TRUE")
    raise KeyError(kind)


def _diag_s1(path, kind):
    rows = _tsv(path)
    if kind == "decay_dist":
        return tuple(sorted(collections.Counter(
            r.get("decay_class") or "" for r in rows).items()))
    raise KeyError(kind)


def _diag_s3(path, kind):
    rows = _tsv(path)
    if kind == "host_branch":
        return sum(1 for r in rows if (r.get("verdict") or "").startswith("host_"))
    raise KeyError(kind)


def _tsv(path):
    out = []
    with open(path, encoding="utf-8") as f:
        hdr = f.readline().rstrip("\n").split("\t")
        for ln in f:
            if not ln.strip():
                continue
            p = ln.rstrip("\n").split("\t")
            p += [""] * (len(hdr) - len(p))
            out.append(dict(zip(hdr, p)))
    return out


# 参数 -> (诊断取自哪个 stage 产物, 诊断名). 诊断必须随该参数移动, 否则仪器可疑.
DIAGNOSTICS = {
    "COMP_AA_MIN": ("s2", "ncomp_dist"),
    "MIN_ALN": ("s2", "with_comps"),
    "MIN_PID": ("s2", "with_comps"),
    "MIN_BITS": ("s2", "te_n"),
    "PROVIRUS_NCOMP": ("s2", "provirus_n"),
    "PROVIRUS_SPAN": ("s2", "provirus_n"),
    "MIN_HEAD_STOPS": ("s1", "decay_dist"),
    "MIN_DECAY_STOPS": ("s1", "decay_dist"),
    "CLEAN_ENRICH": ("s1", "decay_dist"),
    "CLEAN_MAX_STOPS": ("s1", "decay_dist"),
    "MIN_SPAN": ("s1", "decay_dist"),
    # AA_VIRAL_PID 只作用于 s3 的寄主否决分支; 本次扫描是 --no-host (locus 表仅表头),
    # 该分支结构上不可能被触发 —— 所以它测不了, 不是"这个阈值没用"。
    "AA_VIRAL_PID": ("s3", "host_branch"),
}
BLIND_NO_HOST = {"AA_VIRAL_PID"}
# 每个阈值都有一个"按构造必然生效"的值。诊断量在**整个取值网格**里不动时, 再打这一针:
#   打了针动了  -> 阈值在范围内真的不敏感 (合法结论)
#   打针也不动  -> 该阈值在本分支里逻辑冗余(同一 OR 分支里另一个条件恒先满足), 或仪器坏
# 区分这两者看控制实验的结论: MIN_PID/PROVIRUS_NCOMP 已证明注入链路是通的。
IMPOSSIBLE = {
    "COMP_AA_MIN": 10 ** 6, "MIN_ALN": 10 ** 6, "MIN_PID": 101,
    "MIN_BITS": 10 ** 6, "PROVIRUS_NCOMP": 6, "PROVIRUS_SPAN": 10 ** 9,
    "MIN_HEAD_STOPS": 10 ** 6, "MIN_DECAY_STOPS": 10 ** 6,
    "CLEAN_ENRICH": -1.0, "CLEAN_MAX_STOPS": -1, "MIN_SPAN": 10 ** 6,
    "AA_VIRAL_PID": 101.0,
}


def load_module(script, overrides, name="patched"):
    """把 script 读进来, 替换/追加覆盖赋值, exec 成独立模块 (不动磁盘源码).

    需要两种覆盖手段, 缺一不可:
      * 追加赋值 —— 常数在函数体内被读取, 调用时查模块全局, 追加在末尾即生效;
      * 源码内替换 —— 形如 `def f(x=CONST)` 的**默认参数在 def 时求值**, 之后再改全局
        没用 (`credit_components(hsps, comp_aa_min=COMP_AA_MIN)` 就是这种)。所以先把
        源码里所有 `=CONST` 引用换掉 (负向后顾排除 `==CONST` 比较), 再追加。
      派生量 (COMP_NT_MIN = COMP_AA_MIN*3) 一并重算。
    """
    path = os.path.join(MOD, script)
    with open(path, encoding="utf-8") as f:
        src = f.read()
    for k, v in overrides.items():
        # 只替换赋值/默认参数形式的引用, 不碰 == 比较
        src, n = re.subn(r"(?<!=)=%s\b" % re.escape(k), "=%r" % (v,), src)
    extra = ["", "# ---- 敏感度扫描覆盖 (运行时注入, 不改源文件) ----"]
    for k, v in overrides.items():
        extra.append("%s = %r" % (k, v))
        for d in DERIVED.get(k, []):
            extra.append("%s = %s * 3" % (d, k))
    src = "\n".join([src] + extra) + "\n"
    mod = types.ModuleType(name)
    mod.__dict__["__file__"] = path
    mod.__dict__["__name__"] = name
    exec(compile(src, path, "exec"), mod.__dict__)
    # 自检一: 模块全局确实被改成了新值。
    # 但**这一条不够** —— 见 control(): 默认参数在 def 时求值, 全局改了也到不了代码路径。
    for k, v in overrides.items():
        got = mod.__dict__.get(k)
        if got != v:
            raise AssertionError("参数覆盖未生效: %s=%r (期望 %r)" % (k, got, v))
    return mod


def run_main(mod, argv):
    old = sys.argv
    sys.argv = [mod.__file__] + argv
    try:
        mod.main()
    finally:
        sys.argv = old


def run_probe(param, value, work):
    """打一针不可能值, 返回该参数对应的下游诊断量 (供汇总判定用)."""
    os.makedirs(work, exist_ok=True)
    q = os.path.join(WORK_RUN, "q.fa")
    locus = os.path.join(WORK_CACHE, "locus_nohost.tsv")
    dstage, dname = DIAGNOSTICS[param]
    stage = PLAN[param][0]
    if stage == "s1":
        f = os.path.join(work, "s1.tsv")
        run_main(load_module("s1_decay_scan.py", {param: value}),
                 STAGES["s1"][1](q, WORK_RUN, f))
        return _diag_s1(f, dname)
    if stage == "s2":
        f = os.path.join(work, "s2.tsv")
        run_main(load_module("s2_domain_scan.py", {param: value}),
                 STAGES["s2"][1](q, WORK_RUN, f))
        return _diag_s2(f, dname)
    v = os.path.join(work, "verdict.tsv")
    base_s1 = os.path.join(WORK_CACHE, "s1_base.tsv")
    base_s2 = os.path.join(WORK_CACHE, "s2_base.tsv")
    run_main(load_module("s3_verdict.py", {param: value}),
             [base_s1, base_s2, locus, WORK_EVID, v])
    return _diag_s3(v, dname)


def control(run_dir, evidence, work):
    """仪器自检: 用**不可能值**证明注入真的到达了代码路径.

    这是把两种情形分开的唯一办法 ——
      * "阈值在很宽的范围内都不影响结论" (真发现)
      * "我的注入静默失败了"            (假发现, 而且看起来一模一样)
    做法: 注入一个按构造必须生效的值, 断言**下游诊断量**随之移动:

      MIN_PID = 101   一致度上限是 100, 任何 HSP 都不该留下 -> 有组件的候选数必须归零
      PROVIRUS_NCOMP = 6  CANON 只有 5 个组件 -> provirus_scale 必须全 FALSE

    MIN_PID 那条还会带动 verdict 大改; PROVIRUS_NCOMP=6 那条则**不**该改 verdict
    (正好用来印证"provirus_scale 不是判据"), 但它必须改诊断量 —— 否则仪器就是坏的。
    """
    os.makedirs(work, exist_ok=True)
    q = os.path.join(run_dir, "q.fa")
    locus = os.path.join(work, "locus_nohost.tsv")
    with open(locus, "w", encoding="utf-8") as fo:
        old = sys.stdout
        sys.stdout = fo
        try:
            run_main(load_module("s2b_locus_scan.py", {}), ["--emit-header"])
        finally:
            sys.stdout = old

    def stage(param, val, tag):
        """注入一次, 返回 (诊断量, verdict 表路径)."""
        d = os.path.join(work, tag)
        os.makedirs(d, exist_ok=True)
        s1_f = os.path.join(d, "s1.tsv")
        s2_f = os.path.join(d, "s2.tsv")
        v_f = os.path.join(d, "verdict.tsv")
        run_main(load_module("s1_decay_scan.py", {}), STAGES["s1"][1](q, run_dir, s1_f))
        run_main(load_module("s2_domain_scan.py", {param: val}),
                 STAGES["s2"][1](q, run_dir, s2_f))
        run_main(load_module("s3_verdict.py", {}),
                 [s1_f, s2_f, locus, evidence, v_f])
        dstage, dname = DIAGNOSTICS[param]
        diag = _diag_s2(s2_f, dname) if dstage == "s2" else _diag_s3(v_f, dname)
        return diag, v_f

    print("=" * 74)
    print("sensitivity 仪器自检: 不可能值必须产生可观测的下游变化")
    print("=" * 74)
    fails = []

    base_comp, base_v = stage("MIN_PID", 25, "ctl_base_pid")
    kill_comp, kill_v = stage("MIN_PID", 101, "ctl_pid101")
    bv, kv = load_verdict(base_v), load_verdict(kill_v)
    changed = sum(1 for k in bv if bv.get(k) != kv.get(k))
    print("  MIN_PID=101 (不可能值): 有组件的候选 %s -> %s, verdict 变更 %d 条"
          % (base_comp, kill_comp, changed))
    # 断言要落在**结构性的必然结果**上, 不要用"改了百分之多少"这种拍出来的阈值:
    # 一致度上限 100, 所以一个 HSP 都不该留下 -> 有组件的候选必须为 0。
    if kill_comp != 0:
        fails.append("MIN_PID=101 仍留下 %s 个有组件的候选" % kill_comp)
    if changed == 0:
        fails.append("MIN_PID=101 一条 verdict 都没改, 注入没到代码路径")

    base_pro, base_pv = stage("PROVIRUS_NCOMP", 4, "ctl_base_pnc")
    kill_pro, kill_pv = stage("PROVIRUS_NCOMP", 6, "ctl_pnc6")
    pv, kv2 = load_verdict(base_pv), load_verdict(kill_pv)
    changed2 = sum(1 for k in pv if pv.get(k) != kv2.get(k))
    print("  PROVIRUS_NCOMP=6 (不可能值): provirus_scale TRUE %s -> %s, "
          "verdict 变更 %d 条" % (base_pro, kill_pro, changed2))
    if kill_pro != 0:
        fails.append("PROVIRUS_NCOMP=6 仍有 %s 条 provirus_scale=TRUE" % kill_pro)
    if changed2:
        print("     (verdict 也变了 %d 条 — 说明 provirus_scale 其实参与判定, "
              "与「它已不是判据」的结论冲突, 需要复核)" % changed2)

    print()
    if fails:
        for f in fails:
            print("  **失败**: %s" % f)
        print("\n!! 仪器自检未通过, 本次扫描的任何「无影响」结论都不可信。")
        return 1
    print("  自检通过: 注入确实到达代码路径, 扫描结果可读。")
    print("  附带结论: PROVIRUS_NCOMP=6 使 provirus_scale 全 FALSE 而 verdict 不变,")
    print("           独立印证了「该列已不参与判定」。")
    return 0


def load_verdict(path):
    """读 s3 verdict 表 -> {contig: verdict}."""
    out = {}
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return out
    with open(path, encoding="utf-8") as f:
        hdr = f.readline().rstrip("\n").split("\t")
        if "verdict" not in hdr:
            return out
        vi, ci = hdr.index("verdict"), hdr.index("contig_id")
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) > max(vi, ci):
                out[p[ci]] = p[vi]
    return out


def tabulate(v, order):
    c = collections.Counter(v.values())
    return {k: c.get(k, 0) for k in order}


def main():
    ap = argparse.ArgumentParser(description="EVE 判别阈值敏感度扫描")
    ap.add_argument("--run-dir", required=True)
    ap.add_argument("--evidence", required=True)
    ap.add_argument("--work", required=True)
    ap.add_argument("--only", default="", help="只扫某个参数")
    ap.add_argument("--json", default="", dest="json_out")
    ap.add_argument("--control", action="store_true",
                    help="只跑仪器自检 (不可能值必须产生可观测的下游变化), 不扫描")
    a = ap.parse_args()
    os.makedirs(a.work, exist_ok=True)
    global WORK_RUN, WORK_EVID, WORK_CACHE
    WORK_RUN, WORK_EVID, WORK_CACHE = a.run_dir, a.evidence, a.work
    if a.control:
        return control(a.run_dir, a.evidence, a.work)
    q = os.path.join(a.run_dir, "q.fa")
    # 关寄主通道: 本扫描只测 s1/s2/s3 的判定, 寄主侧要另跑
    locus = os.path.join(a.work, "locus_nohost.tsv")
    # 五个阶段的 main() 都是无参、直接读 sys.argv, 所以统一用 run_main 注入 argv.
    # s2b --emit-header 把表头打到 stdout, 这里重定向到文件.
    s2b = load_module("s2b_locus_scan.py", {})
    with open(locus, "w", encoding="utf-8") as fo:
        old = sys.stdout
        sys.stdout = fo
        try:
            run_main(s2b, ["--emit-header"])
        finally:
            sys.stdout = old
    print(f"[sens] 寄主通道关闭 (locus 表仅表头), 输出目录 {a.work}\n")

    base_s1 = os.path.join(a.work, "s1_base.tsv")
    base_s2 = os.path.join(a.work, "s2_base.tsv")
    run_main(load_module("s1_decay_scan.py", {}), STAGES["s1"][1](q, a.run_dir, base_s1))
    run_main(load_module("s2_domain_scan.py", {}), STAGES["s2"][1](q, a.run_dir, base_s2))
    base_v = os.path.join(a.work, "verdict_base.tsv")
    s3 = load_module("s3_verdict.py", {})
    run_main(s3, [base_s1, base_s2, locus, a.evidence, base_v])
    baseline = load_verdict(base_v)
    order = sorted(set(baseline.values()))
    print(f"基线 verdict 分布 ({len(baseline)} 条):")
    for k, n in sorted(tabulate(baseline, order).items(), key=lambda x: -x[1]):
        print(f"    {k:<30}{n:>6}")
    print()

    results = {}
    for param, (stage, values) in PLAN.items():
        if a.only and a.only != param:
            continue
        rows = []
        for val in values:
            d = os.path.join(a.work, "%s_%s" % (param, str(val).replace(".", "p")))
            os.makedirs(d, exist_ok=True)
            s1_f, s2_f = base_s1, base_s2
            if stage == "s1":
                s1_f = os.path.join(d, "s1.tsv")
                run_main(load_module("s1_decay_scan.py", {param: val}),
                         STAGES["s1"][1](q, a.run_dir, s1_f))
            elif stage == "s2":
                s2_f = os.path.join(d, "s2.tsv")
                run_main(load_module("s2_domain_scan.py", {param: val}),
                         STAGES["s2"][1](q, a.run_dir, s2_f))
            v_f = os.path.join(d, "verdict.tsv")
            ov = {param: val} if stage == "s3" else {}
            run_main(load_module("s3_verdict.py", ov),
                     [s1_f, s2_f, locus, a.evidence, v_f])
            cur = load_verdict(v_f)
            changed = sum(1 for k in baseline if baseline.get(k) != cur.get(k))
            tab = tabulate(cur, sorted(set(order) | set(cur.values())))
            dstage, dname = DIAGNOSTICS.get(param, (None, None))
            diag = None
            if dstage == "s2":
                diag = _diag_s2(s2_f, dname)
            elif dstage == "s1":
                diag = _diag_s1(s1_f, dname)
            elif dstage == "s3":
                diag = _diag_s3(v_f, dname)
            rows.append({"value": val, "changed": changed, "dist": tab,
                         "diag": diag})
            results.setdefault(param, rows)
            print("  %-18s = %-8s 变更 %4d / %d 条" %
                  (param, val, changed, len(baseline)))
        print()
    print("=" * 74)
    print("敏感度汇总 (变更条数 = 相对基线的 verdict 差异; 越大说明该阈值越关键)")
    print("=" * 74)
    print(f"{'参数':<20}{'最小变更':>10}{'最大变更':>10}{'跨度占比':>10}   取值->变更")
    for param, rows in results.items():
        ch = [r["changed"] for r in rows]
        span = (max(ch) - min(ch)) / max(len(baseline), 1)
        diags = [r.get("diag") for r in rows]
        moved = len(set(map(repr, diags))) > 1
        if param in BLIND_NO_HOST:
            verdict = "不可测 (本次关寄主通道, 该分支结构上不触发)"
        elif not moved:
            # 网格内不动 -> 补一针不可能值, 区分"范围外才生效"与"逻辑冗余"
            probe = IMPOSSIBLE.get(param)
            pd = None
            if probe is not None:
                pd = run_probe(param, probe, os.path.join(a.work, "probe_" + param))
            if pd is None:
                verdict = "**仪器可疑**: 参数变了而下游量一动不动"
            elif pd != diags[0]:
                verdict = ("范围外才生效: 网格内不动, 打不可能值(%s)诊断量会动 "
                           "-> 阈值在此范围内不敏感" % probe)
            else:
                verdict = ("**逻辑冗余**: 连不可能值(%s)都不改诊断量 "
                           "-> 同一分支里另有条件恒先满足" % probe)
        elif max(ch) == 0:
            verdict = "真无影响 (诊断量确实动了)"
        else:
            verdict = "有效"
        print(f"{param:<20}{min(ch):>10}{max(ch):>10}{span:>10.1%}   "
              + " ".join("%s:%d" % (r["value"], r["changed"]) for r in rows))
        print(f"{'':<20}诊断量移动: {'是' if moved else '否'}   -> {verdict}")
    if a.json_out:
        with open(a.json_out, "w", encoding="utf-8") as f:
            json.dump({"baseline": baseline, "results": results}, f,
                      ensure_ascii=False, indent=1)
        print(f"\n[sens] -> {a.json_out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
