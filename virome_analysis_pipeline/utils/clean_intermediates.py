#!/usr/bin/env python3
"""clean_intermediates.py — known-virus 管线中间文件清理 (默认 dry-run)

只清 **成功完成任务** 的中间文件，失败/未完成任务一律不动。

成功判据（逐任务）:
  <run>/11.Ultimate_Circular_Result.fasta 存在、非空、且含至少 1 条序列
  (仅 header 无序列的 header-only 残片判为失败, 不清)

默认清理层 (--layer steps):
  <05_Assembly>/<病毒>/<样本>/<run>/<编号目录>/   e.g. 1.DeNovo_Assembly …
  终产物 11.Ultimate / Fully-assembled.ok / Global_Evolution_Stats.tsv /
  logs/ 与 06_Extraction/ 全部保留;
  12.Coverage_Visualization/ 是报告 Stage 5 内嵌的图件, 同样保留

可选层 (显式 --layer 指定, 默认不删):
  virus-bam     03_Extract/virus-bam          (仅当不再跑 --stage variants)
  virus-reads   03_Extract/virus_reads        (仅当不再重跑 --stage full)
  detect-index  01_Detection/index, bam       (仅当不再重跑 --stage detect)

用法:
  python clean_intermediates.py -o <输出根目录>                    # dry-run, 只报告
  python clean_intermediates.py -o <输出根目录> --apply            # 真删 steps 层
  python clean_intermediates.py -o <输出根目录> --apply --layer steps,virus-bam
  python clean_intermediates.py -o <输出根目录> --apply --layer all   # 全部可选层
"""
import argparse
import re
import shutil
import sys
from pathlib import Path

NUM_DIR_RE = re.compile(r"^\d+\.")
PRODUCT = "11.Ultimate_Circular_Result.fasta"

# 编号目录里属于**交付物**而非中间件的: 12.Coverage_Visualization 是
# 逐任务覆盖度双轨图, 被报告 Stage 5 区块内嵌引用 —— 默认保留
KEEP_DIRS = {"12.Coverage_Visualization"}


def dir_size(p: Path) -> int:
    """递归累计字节数 (硬链接按一次算, 失败静默跳过)。"""
    total, seen = 0, set()
    try:
        for f in p.rglob("*"):
            try:
                if not f.is_file() or f.is_symlink():
                    continue
                st = f.stat()
                key = (st.st_dev, st.st_ino)
                if key in seen:
                    continue
                seen.add(key)
                total += st.st_size
            except OSError:
                continue
    except OSError:
        pass
    return total


def fmt(n: int) -> str:
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if n < 1024 or unit == "TB":
            return f"{n:.1f}{unit}"
        n /= 1024.0
    return f"{n:.1f}TB"


def product_ok(run_dir: Path) -> tuple:
    """返回 (ok, 说明)。ok=True 表示任务确实成功、可清中间件。"""
    hits = [run_dir / PRODUCT] if (run_dir / PRODUCT).is_file() else list(run_dir.rglob(PRODUCT))
    if not hits:
        return False, "无终产物"
    best = max(hits, key=lambda f: f.stat().st_size)
    try:
        if best.stat().st_size == 0:
            return False, "终产物 0 字节"
        seq_chars, has_head = 0, False
        with open(best, "r", errors="ignore") as fh:
            for line in fh:
                if line.startswith(">"):
                    has_head = True
                elif has_head:
                    seq_chars += len(line.strip())
                if seq_chars > 0:
                    return True, f"ok ({best.stat().st_size}B)"
    except OSError as e:
        return False, f"读取失败: {e}"
    return False, "header-only (无序列)"


def find_run_dirs(assembly: Path):
    """yield (virus_dir, sample_dir, run_dir) —— run_dir 内应含编号步骤目录。"""
    for virus_dir in sorted(d for d in assembly.iterdir() if d.is_dir()):
        for sample_dir in sorted(d for d in virus_dir.iterdir() if d.is_dir()):
            for run_dir in sorted(d for d in sample_dir.iterdir() if d.is_dir()):
                if any(NUM_DIR_RE.match(s.name) and s.is_dir() for s in run_dir.iterdir()):
                    yield virus_dir, sample_dir, run_dir


def plan_steps(assembly: Path, log):
    """返回 (待删目录列表, 成功任务数, 跳过任务数)。"""
    targets, n_ok, n_skip = [], 0, 0
    for virus_dir, sample_dir, run_dir in find_run_dirs(assembly):
        ok, why = product_ok(run_dir)
        tag = f"{virus_dir.name}/{sample_dir.name}/{run_dir.name}"
        if not ok:
            n_skip += 1
            log(f"  ⏭️  跳过 (未成功): {tag} — {why}")
            continue
        n_ok += 1
        steps = [s for s in run_dir.iterdir()
                 if s.is_dir() and NUM_DIR_RE.match(s.name) and s.name not in KEEP_DIRS]
        if not steps:
            continue
        sz = sum(dir_size(s) for s in steps)
        log(f"  ✅ {tag} — {len(steps)} 个步骤目录, {fmt(sz)}")
        targets.extend(steps)
    return targets, n_ok, n_skip


def plan_simple(paths, label, log):
    """可选层: 整目录删除 (只有存在才计入)。"""
    out = []
    for p in paths:
        if p.exists():
            sz = dir_size(p)
            log(f"  🗑️  [{label}] {p} — {fmt(sz)}")
            out.append(p)
    return out


def main():
    ap = argparse.ArgumentParser(
        description="known-virus 中间文件清理 (只清成功完成任务; 默认 dry-run)",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-o", "--output_dir", required=True, help="管线输出根目录 (含 05_Assembly)")
    ap.add_argument("--apply", action="store_true", help="真正删除 (默认只报告)")
    ap.add_argument("--layer", default="steps",
                    help="清理层, 逗号分隔: steps,virus-bam,virus-reads,detect-index,all (默认 steps)")
    ap.add_argument("--keep-bam-for-variants", action="store_true",
                    help="保护 virus-bam (计划以后再跑 --stage variants 时用)")
    args = ap.parse_args()

    out = Path(args.output_dir).resolve()
    assembly = out / "05_Assembly"
    if not assembly.is_dir():
        sys.exit(f"ERROR: 找不到 {assembly} (先跑 --stage full)")

    layers = {s.strip().lower() for s in args.layer.split(",") if s.strip()}
    if "all" in layers:
        layers |= {"steps", "virus-bam", "virus-reads", "detect-index"}

    mode = "APPLY (真删)" if args.apply else "DRY-RUN (只报告, 加 --apply 才删)"
    print(f"=== known-virus 中间文件清理 | {mode} ===")
    print(f"输出根: {out}")
    print(f"层: {','.join(sorted(layers))}\n")

    def log(msg):
        print(msg, flush=True)

    targets, n_ok, n_skip = [], 0, 0
    if "steps" in layers:
        print("[层 1] 05_Assembly 编号步骤目录 (只清成功任务)")
        targets, n_ok, n_skip = plan_steps(assembly, log)
    else:
        # 未清 steps 层时也统计成功/跳过, 保持报告口径一致
        for _v, _s, _r in find_run_dirs(assembly):
            ok, _why = product_ok(_r)
            n_ok += 1 if ok else 0
            n_skip += 0 if ok else 1

    if "virus-bam" in layers:
        print("\n[层 2] 03_Extract/virus-bam")
        if args.keep_bam_for_variants:
            log("  ⏭️  已用 --keep-bam-for-variants 保护, 跳过")
        else:
            targets += plan_simple([out / "03_Extract" / "virus-bam"], "virus-bam", log)
    if "virus-reads" in layers:
        print("\n[层 3] 03_Extract/virus_reads")
        targets += plan_simple([out / "03_Extract" / "virus_reads"], "virus-reads", log)
    if "detect-index" in layers:
        print("\n[层 4] 01_Detection 索引与 BAM")
        targets += plan_simple([out / "01_Detection" / "index", out / "01_Detection" / "bam"],
                               "detect-index", log)

    print("\n" + "=" * 60)
    if "steps" in layers:
        print(f"任务账本: 成功 {n_ok} | 跳过(未成功) {n_skip} | 待删目录 {len(targets)}")
    total = sum(dir_size(t) for t in targets)
    print(f"可回收空间: {fmt(total)}")

    if not args.apply:
        print("\n[DRY-RUN] 未删除任何文件。确认无误后加 --apply 执行。")
        return

    if not targets:
        print("无可清理对象。")
        return

    freed, failed = 0, 0
    for t in targets:
        sz = dir_size(t)
        try:
            shutil.rmtree(t)
            freed += sz
            log(f"  ✔ 已删 {t.relative_to(out)} ({fmt(sz)})")
        except OSError as e:
            failed += 1
            log(f"  ✘ 删除失败 {t}: {e}")
    print("\n" + "=" * 60)
    print(f"完成: 释放 {fmt(freed)} | 失败 {failed} 个")
    print("保留: 11.Ultimate_Circular_Result.fasta / Fully-assembled.ok / "
          "Global_Evolution_Stats.tsv / logs / 06_Extraction 全部产物")


if __name__ == "__main__":
    main()
