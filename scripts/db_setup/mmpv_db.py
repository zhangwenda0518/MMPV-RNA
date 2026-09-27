#!/usr/bin/env python3
"""
MMPV 数据库部署 / 迁移 / 校验工具

把 MMPV 全管线依赖的数据库集中到一个「共享库根」，使：(1) 换机器/换用户时可复现，
(2) 旧路径（~/database、~/plant_virus_db）通过软链接继续可用，核心代码里那些
硬编码路径无需改动。

用法速览
--------
    python scripts/db_setup/mmpv_db.py list                  # 打印依赖清单
    python scripts/db_setup/mmpv_db.py env                   # 打印环境变量块
    python scripts/db_setup/mmpv_db.py check                 # 部署前体检
    python scripts/db_setup/mmpv_db.py verify                # 校验现有库是否齐全

    # 换到新服务器：下载能自动获取的，其余靠共享池
    python scripts/db_setup/mmpv_db.py down                  # 演练：看会下哪些
    python scripts/db_setup/mmpv_db.py down --apply          # 实际下载/构建
    python scripts/db_setup/mmpv_db.py down --only virus_db_genomad --apply

    # 本机已有数据库：迁移到共享库根 + 建兼容软链接
    python scripts/db_setup/mmpv_db.py layout --apply
    python scripts/db_setup/mmpv_db.py adopt --db-from ~/database --apply
    python scripts/db_setup/mmpv_db.py link --apply
    python scripts/db_setup/mmpv_db.py publish --apply       # 仓库侧小文件入池

设计要点见 doc/DB_DEPLOYMENT.md。仅用标准库，无第三方依赖。
"""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
import unicodedata
from dataclasses import dataclass
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parent.parent
MANIFEST = HERE / "mmpv_db_manifest.tsv"
FETCH = HERE / "mmpv_db_fetch.tsv"

# 旧路径（核心代码硬编码引用），迁移后由 link 子命令指到共享根
LEGACY_DB = Path.home() / "database"
LEGACY_PLANT = Path.home() / "plant_virus_db"

REPO_OVERLAY_DIRNAME = "repo-overlay"


# ────────────────────────────── 清单 ──────────────────────────────

@dataclass
class Entry:
    id: str
    kind: str          # download | build | local | repo
    type: str          # dir | file
    root: str          # db | plant | repo
    path: str
    probe: str
    size: str
    required: bool
    note: str

    @property
    def required_mark(self) -> str:
        return "必需" if self.required else "可选"


FIELDS = 9


def load_manifest(path: Path = MANIFEST) -> list[Entry]:
    if not path.is_file():
        sys.exit(f"清单不存在: {path}")
    entries: list[Entry] = []
    seen: set[str] = set()
    with path.open(encoding="utf-8") as fh:
        for lineno, raw in enumerate(fh, 1):
            line = raw.rstrip("\n")
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            f = line.split("\t")
            if f[0].strip() == "id" and lineno < 5:
                continue          # 表头：留给 pandas 等外部读取，加载器跳过
            if len(f) < FIELDS:
                sys.exit(f"清单第 {lineno} 行字段不足（需 {FIELDS} 列，实得 {len(f)}）: {line!r}")
            e = Entry(*[x.strip() for x in f[:FIELDS]])
            e.required = e.required.strip().lower() in ("yes", "y", "true", "1")
            if e.id in seen:
                sys.exit(f"清单第 {lineno} 行 id 重复: {e.id}")
            if e.root not in ("db", "plant", "repo", "src"):
                sys.exit(f"清单第 {lineno} 行 root 非法: {e.root}（须为 db|plant|repo|src）")
            if e.type not in ("dir", "file"):
                sys.exit(f"清单第 {lineno} 行 type 非法: {e.type}（须为 dir|file）")
            if e.path.startswith("/") or ".." in Path(e.path).parts:
                sys.exit(f"清单第 {lineno} 行 path 必须是 root 下的相对路径: {e.path}")
            seen.add(e.id)
            entries.append(e)
    return entries


# ─────────────────────────── 根目录解析 ───────────────────────────

@dataclass
class Roots:
    share: Path
    repo: Path
    db: Path
    plant: Path
    src: Path

    def resolve(self, e: Entry) -> Path:
        if e.root == "db":
            return self.db / e.path
        if e.root == "plant":
            return self.plant / e.path
        if e.root == "src":
            return self.src / e.path
        return self.repo / e.path

    def overlay(self, e: Entry) -> Path:
        """仓库侧小文件在共享池中的位置（剥掉前导 database/）。"""
        rel = e.path
        if rel.startswith("database/"):
            rel = rel[len("database/"):]
        return self.share / REPO_OVERLAY_DIRNAME / rel


def resolve_roots(args: argparse.Namespace) -> Roots:
    if getattr(args, "share", None):
        share = Path(args.share).expanduser()
    elif os.environ.get("MMPV_DATA_SHARE"):
        share = Path(os.environ["MMPV_DATA_SHARE"]).expanduser()
    else:
        share = Path.home() / "mmpv-db"

    repo = Path(args.repo).expanduser() if getattr(args, "repo", None) else REPO_ROOT

    db = Path(args.db_root).expanduser() if getattr(args, "db_root", None) \
        else Path(os.environ.get("MMPV_DB_ROOT") or (share / "database"))
    plant = Path(args.plant_root).expanduser() if getattr(args, "plant_root", None) \
        else Path(os.environ.get("MMPV_PLANT_VIRUS_DB") or (share / "plant_virus_db"))

    src = Path(args.src_root).expanduser() if getattr(args, "src_root", None)         else Path(os.environ.get("MMPV_SRC") or (share / "src"))

    return Roots(share=share, repo=repo, db=db, plant=plant, src=src)


# ─────────────────────────── 状态探测 ───────────────────────────

MISSING, BROKEN, EMPTY_DIR, EMPTY_FILE, WRONG_TYPE, OK = \
    "缺失", "断链", "空目录", "空文件", "类型不符", "OK"


def probe(target: Path, want: str) -> str:
    """判定一个路径是否存在且非空。

    want = "dir"/"file" 时校验类型；"any" 时只查存在与非空。
    能区分断链——迁移/换盘后最常见的故障。
    """
    if target.is_symlink() and not target.exists():
        return BROKEN
    if not target.exists():
        return MISSING
    if target.is_dir():
        if want == "file":
            return WRONG_TYPE
        try:
            next(iter(target.iterdir()))
        except StopIteration:
            return EMPTY_DIR
        except OSError:
            return MISSING
        return OK
    if want == "dir":
        return WRONG_TYPE
    try:
        return OK if target.stat().st_size > 0 else EMPTY_FILE
    except OSError:
        return MISSING


def probe_entry(e: Entry, roots: Roots) -> tuple[str, Path]:
    """两级检查：单元自身须符合声明的 dir/file 类型；哨兵只要求存在且非空。

    哨兵常是单元内部的一个异类文件（如 taxonomy 目录里的 nodes.dmp），
    所以不能把单元的类型约束套到哨兵上。
    """
    base = roots.resolve(e)
    if not e.probe:
        return probe(base, e.type), base
    st = probe(base, e.type)
    if st != OK:
        return st, base
    return probe(base / e.probe, "any"), base


# ─────────────────────────── 输出工具 ───────────────────────────

def dwidth(s: str) -> int:
    """显示宽度：CJK 全角算 2 列，保证中文表格对齐。"""
    return sum(2 if unicodedata.east_asian_width(c) in "WF" else 1 for c in s)


def pad(s: str, n: int) -> str:
    return s + " " * max(0, n - dwidth(s))


def measure(rows: list[list[str]], header: list[str]) -> list[int]:
    widths = [dwidth(c) for c in header]
    for r in rows:
        for i, c in enumerate(r):
            widths[i] = max(widths[i], dwidth(c))
    return widths


def table(rows: list[list[str]], header: list[str], widths: list[int] | None = None) -> str:
    """渲染表格。widths 传入时复用外部列宽，使多个分组表格列对齐一致。"""
    widths = widths or measure(rows, header)
    out = ["  ".join(pad(c, widths[i]) for i, c in enumerate(header)).rstrip()]
    out.append("  ".join("-" * w for w in widths))
    for r in rows:
        out.append("  ".join(pad(c, widths[i]) for i, c in enumerate(r)).rstrip())
    return "\n".join(out)


def info(msg: str) -> None:
    print(msg, flush=True)


def warn(msg: str) -> None:
    print(f"⚠️  {msg}", flush=True)


def die(msg: str) -> None:
    sys.exit(f"❌ {msg}")


# ─────────────────────────── 获取配方 ───────────────────────────

FETCH_FIELDS = 5

# 配方里若有命令用到了 dl，先注入这个断点续传函数（curl 优先，回退 wget）
DL_PREAMBLE = r'''
dl() {
  local url="$1" out="$2"
  mkdir -p "$(dirname "$out")"
  if command -v curl >/dev/null 2>&1; then
    curl -fL --retry 5 --retry-delay 5 -C - -o "$out" "$url"
  elif command -v wget >/dev/null 2>&1; then
    wget -c -O "$out" "$url"
  else
    echo "ERROR: 需要 curl 或 wget 之一" >&2; return 1
  fi
}
'''


@dataclass
class Step:
    id: str
    step: int
    needs: list[str]
    deps: list[str]
    cmd: str


def load_fetch(path: Path = FETCH) -> dict[str, list[Step]]:
    """读获取配方。没有配方的 manifest id = 仓库内无可靠复现途径，down 会如实报告。"""
    recipes: dict[str, list[Step]] = {}
    if not path.is_file():
        return recipes
    with path.open(encoding="utf-8") as fh:
        for lineno, raw in enumerate(fh, 1):
            line = raw.rstrip("\n")
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            f = line.split("\t")
            if f[0].strip() == "id" and lineno < 5:
                continue
            if len(f) < FETCH_FIELDS:
                sys.exit(f"配方第 {lineno} 行字段不足（需 {FETCH_FIELDS} 列，实得 {len(f)}）: {line!r}")
            sid, step, needs, deps, cmd = [x.strip() for x in f[:FETCH_FIELDS]]
            try:
                istep = int(step)
            except ValueError:
                sys.exit(f"配方第 {lineno} 行 step 不是整数: {step!r}")
            recipes.setdefault(sid, []).append(Step(
                id=sid, step=istep,
                needs=[t for t in needs.split(",") if t],
                deps=[d for d in deps.split(",") if d],
                cmd=cmd,
            ))
    for steps in recipes.values():
        steps.sort(key=lambda s: s.step)
    return recipes


BUILD_DB = REPO_ROOT / "virome_discovery_pipeline" / "build_virus_db.py"
DB_BUILD_SCRIPTS = REPO_ROOT / "virome_discovery_pipeline" / "utils" / "db_build"
DB_BUILD_NEEDS = ["SearchAccessionIdToTaxId.py", "db_seqid2taxid_add_legth.py", "make_ktaxonomy.py"]


def missing_tools(needs: list[str]) -> list[str]:
    """返回缺失的依赖。

    'dl'       特指 curl/wget 至少有一个
    'db_build' 特指 build_virus_db.py 能用——即它的 3 个自写辅助脚本已收编进
               utils/db_build/（这是当前「其他用户复现」的硬阻塞点，见该目录 README）
    """
    miss = []
    for t in needs:
        if t == "dl":
            if not (shutil.which("curl") or shutil.which("wget")):
                miss.append("curl|wget")
        elif t == "db_build":
            absent = [s for s in DB_BUILD_NEEDS if not (DB_BUILD_SCRIPTS / s).is_file()]
            if absent:
                miss.append("utils/db_build 缺 " + ",".join(absent))
            elif not BUILD_DB.is_file():
                miss.append("virome_discovery_pipeline/build_virus_db.py")
        elif not shutil.which(t):
            miss.append(t)
    return miss


# ─────────────────────────── 子命令 ───────────────────────────

def cmd_list(args: argparse.Namespace) -> int:
    entries = load_manifest()
    groups = [
        ("download", "预建库（一条命令下载）"),
        ("build", "需本地建索引"),
        ("local", "无法重下，须从共享池拷贝"),
        ("repo", "随仓库的小文件"),
    ]
    header = ["id", "root", "相对路径", "大小", "要求", "说明"]
    all_rows = [
        [e.id, e.root, e.path, e.size, e.required_mark, e.note]
        for _, _ in groups for e in entries if e.kind == _
    ]
    widths = measure(all_rows, header)
    for kind, title in groups:
        group = [e for e in entries if e.kind == kind]
        if not group:
            continue
        info(f"\n【{title}】{len(group)} 项")
        rows = [[e.id, e.root, e.path, e.size, e.required_mark, e.note] for e in group]
        info(table(rows, header, widths))
    req = sum(1 for e in entries if e.required)
    info(f"\n共 {len(entries)} 项；必需 {req} 项，可选 {len(entries) - req} 项。")
    return 0


def cmd_which(args: argparse.Namespace) -> int:
    roots = resolve_roots(args)
    entries = {e.id: e for e in load_manifest()}
    rc = 0
    for eid in args.id:
        e = entries.get(eid)
        if e is None:
            warn(f"未知 id: {eid}")
            rc = 1
            continue
        st, base = probe_entry(e, roots)
        info(f"{pad(eid, 34)} {st:4}  {base}")
    return rc


def cmd_env(args: argparse.Namespace) -> int:
    roots = resolve_roots(args)
    v = roots.db / "virus-db"
    try:
        self_rel = Path(__file__).resolve().relative_to(REPO_ROOT.resolve())
    except ValueError:
        self_rel = Path(__file__).resolve()
    lines = [
        f"# MMPV 数据库环境变量 — 由 {self_rel} env 生成",
        f"# 共享库根: {roots.share}",
        "export MMPV_DATA_SHARE=" + str(roots.share),
        "",
        "# —— pipeline_config.yaml 读取的变量（决定配置层解析结果）——",
        "export MMPV_DB_ROOT=" + str(roots.db),
        "export MMPV_VIRUS_DB=" + str(v),
        "export MMPV_HOST_DB=" + str(roots.db / "host_db"),
        "export MMPV_KRAKEN2_DB=" + str(roots.db / "kraken2" / "k2_pluspfp_20260226"),
        "export MMPV_CHECKV_DB=" + str(v / "checkv-db-v1.7"),
        "export MMPV_PLANT_VIRUS_DB=" + str(roots.plant),
        "# 上游原始数据（RVDB fasta 等）；建完索引可删",
        "export MMPV_SRC=" + str(roots.src),
        "",
        "# —— 核心代码绕过配置、直接读的变量（不设就会回落到 ~/database 硬编码）——",
        "export MMPV_VMR=" + str(v / "VMR_MSL41.v1.20260320.xlsx"),
        "export MMPV_GENUS_FAMILY_REF=" + str(roots.db / "taxonomy" / "genus_family_ref.tsv"),
        "export MMPV_NT_DB=" + str(v / "ncbi-virus" / "ncbi_virus.nucl.fasta"),
        "export MMPV_AA_DB=" + str(v / "ncbi-virus" / "ncbi_virus.prot.fasta"),
        "export MMPV_CDD_TAXID_TSV=" + str(roots.repo / "database/cdd/cdd_classified_taxid.tsv"),
        "# filter_virus.py 用的是不带 .dmnd 后缀的前缀",
        "export UNIPROT_DB=" + str(roots.db / "uniport_db" / "uniref90" / "uniref90"),
        "export VIRUS_TAXID=" + str(roots.repo / "database/uniprot/virus.taxid.txt"),
        "",
    ]
    block = "\n".join(lines)

    if args.write:
        dest = Path(args.write).expanduser()
        dest.parent.mkdir(parents=True, exist_ok=True)
        dest.write_text(block, encoding="utf-8")
        info(f"已写入 {dest}")
        info(f"启用：  echo 'source {dest}' >> ~/.bashrc")
    else:
        print(block)
    return 0


def cmd_check(args: argparse.Namespace) -> int:
    roots = resolve_roots(args)
    entries = load_manifest()
    info("== MMPV 数据库体检 ==")
    info(f"清单      : {MANIFEST}  （{len(entries)} 项）")
    info(f"共享库根  : {roots.share}")
    info(f"  DB_ROOT : {roots.db}")
    info(f"  PLANT   : {roots.plant}")
    info(f"仓库      : {roots.repo}")

    # 共享根可写性
    probe_dir = roots.share if roots.share.exists() else roots.share.parent
    if os.access(probe_dir, os.W_OK):
        info(f"可写性    : OK（{probe_dir}）")
    else:
        warn(f"不可写: {probe_dir} —— layout/adopt/link 会失败；verify 不受影响")

    # 旧路径现状：决定 adopt 该做什么
    info("\n-- 旧路径现状（核心代码硬编码引用）--")
    for legacy in (LEGACY_DB, LEGACY_PLANT):
        if legacy.is_symlink():
            dest = os.readlink(legacy)
            ok = "→ 目标存在" if legacy.exists() else "→ 目标缺失（断链！）"
            info(f"  {legacy}  软链接  {dest}  {ok}")
        elif legacy.is_dir():
            n = sum(1 for _ in legacy.iterdir())
            info(f"  {legacy}  真实目录  含 {n} 个条目  → 可用 adopt 迁移")
        else:
            info(f"  {legacy}  不存在  → 新用户走 layout + fetch；老机器请确认库在哪")

    # 与共享根的包含关系
    for legacy in (LEGACY_DB, LEGACY_PLANT):
        if legacy.is_dir() and not legacy.is_symlink():
            try:
                legacy.relative_to(roots.share)
                warn(f"{legacy} 在共享根内部，adopt 会拒绝（避免自我搬迁）")
            except ValueError:
                pass

    # 逐项统计
    stats: dict[str, int] = {}
    missing_required: list[Entry] = []
    for e in entries:
        st, _ = probe_entry(e, roots)
        stats[st] = stats.get(st, 0) + 1
        if e.required and st != OK:
            missing_required.append(e)

    info("\n-- 逐项状态 --")
    info("  " + "  ".join(f"{k}:{v}" for k, v in sorted(stats.items())))
    if missing_required:
        info(f"\n缺必需项 {len(missing_required)} 个：")
        for e in missing_required[:20]:
            st, base = probe_entry(e, roots)
            info(f"  [{st}] {e.id}  ({e.kind})  {base}")
        if len(missing_required) > 20:
            info(f"  ... 另有 {len(missing_required) - 20} 个")
        info("\n跑 `verify` 看完整表格；跑 `layout` 建骨架，`adopt`/`fetch` 补数据。")
        return 1
    info("\n全部必需项就位。")
    return 0


def cmd_verify(args: argparse.Namespace) -> int:
    roots = resolve_roots(args)
    entries = load_manifest()
    wanted = set(args.only) if args.only else None

    # 旧路径悬空是迁移后最常见也最难查的故障（共享盘没挂上/被挪走），先于逐项校验报出来
    for legacy in (LEGACY_DB, LEGACY_PLANT):
        if legacy.is_symlink() and not legacy.exists():
            warn(f"旧路径软链接已悬空: {legacy} → {os.readlink(legacy)}"
                 f"（共享库根未挂载或被移动；核心代码的硬编码路径会全部失效）")

    rows, bad_req, bad_opt = [], [], []
    for e in entries:
        if wanted and e.id not in wanted:
            continue
        st, base = probe_entry(e, roots)
        rows.append([st, e.id, e.kind, e.required_mark, str(base)])
        if st != OK:
            (bad_req if e.required else bad_opt).append(e.id)

    info(table(rows, ["状态", "id", "类型", "要求", "解析路径"]))

    total = len(rows)
    ok = total - len(bad_req) - len(bad_opt)
    info(f"\n合计 {total} 项：OK {ok}，缺失/异常 必需 {len(bad_req)} + 可选 {len(bad_opt)}")
    if bad_req:
        info("缺必需项: " + ", ".join(bad_req))
        return 1
    if bad_opt:
        info("缺可选项（对应工具会自动降级或跳过）: " + ", ".join(bad_opt))
    return 0


def cmd_layout(args: argparse.Namespace) -> int:
    roots = resolve_roots(args)
    entries = load_manifest()
    plan: list[Path] = []
    for root in (roots.db, roots.plant, roots.src, roots.share / REPO_OVERLAY_DIRNAME):
        plan.append(root)
    for e in entries:
        if e.root == "repo":
            continue
        p = roots.resolve(e)
        # dir 单元：建出该目录本身；file 单元：只建父目录，文件留给 adopt/fetch
        plan.append(p.parent if e.type == "file" else p)

    created = 0
    for p in sorted(set(plan)):
        if p.exists():
            continue
        if args.apply:
            p.mkdir(parents=True, exist_ok=True)
        created += 1
        info(("  [建] " if args.apply else "  [待建] ") + str(p))
    info(f"\n{'已创建' if args.apply else '将创建'} {created} 个目录"
         f"（共 {len(set(plan))} 个）。" + ("" if args.apply else " 加 --apply 落盘。"))
    return 0


def _top_level_units(src: Path) -> list[Path]:
    return sorted((p for p in src.iterdir() if not p.name.startswith(".")), key=lambda p: p.name)


def _same_filesystem(a: Path, b: Path) -> bool:
    try:
        probe = b if b.exists() else b.parent
        while not probe.exists() and probe != probe.parent:
            probe = probe.parent
        return a.stat().st_dev == probe.stat().st_dev
    except OSError:
        return False


def _empty_dir(p: Path) -> bool:
    return p.is_dir() and not p.is_symlink() and not any(p.iterdir())


def _place_unit(src: Path, dst: Path, *, move: bool, apply: bool, cross: bool,
                indent: str, stats: dict[str, int]) -> None:
    """把一个单元安置到 dst；目标已存在时递归下潜合并，绝不覆盖已有数据。

    这样 `layout`（先建空骨架）与 `adopt`（再灌数据）可以组合使用：
    layout 建出的空目录会被就地替换成真实内容或软链接。
    """
    if dst.is_symlink():
        warn(f"{indent}{dst.name}: 已是软链接（→ {os.readlink(dst)}），跳过")
        stats["skip"] += 1
        return

    if not dst.exists():
        pass                                  # 直接安置
    elif _empty_dir(dst) and src.is_dir():
        dst.rmdir()                           # layout 留下的空壳，让位
    elif dst.is_dir() and src.is_dir():
        info(f"{indent}{dst.name}/ 已有内容，下潜合并")
        for child in sorted(src.iterdir(), key=lambda p: p.name):
            _place_unit(child, dst / child.name, move=move, apply=apply,
                        cross=cross, indent=indent + "  ", stats=stats)
        # 迁移模式下，子项搬空后要顺手删掉遗留的空壳，否则旧路径永远退不掉
        if apply and move and not any(src.iterdir()):
            src.rmdir()
        return
    else:
        warn(f"{indent}{dst.name}: 目标已存在且非空目录，跳过（不覆盖）")
        stats["skip"] += 1
        return

    if move:
        tail = "（跨文件系统，rsync 拷贝）" if cross else "（同文件系统，秒级）"
        info(f"{indent}[迁] {src.name}{tail}")
        stats["move"] += 1
        if apply:
            dst.parent.mkdir(parents=True, exist_ok=True)
            if cross and shutil.which("rsync"):
                subprocess.run(["rsync", "-a", str(src) + "/", str(dst) + "/"], check=True)
                shutil.rmtree(src)
            else:
                subprocess.run(["mv", str(src), str(dst)], check=True)
    else:
        info(f"{indent}[链] {dst.name} → {src}")
        stats["link"] += 1
        if apply:
            dst.parent.mkdir(parents=True, exist_ok=True)
            dst.symlink_to(src.resolve())


def cmd_adopt(args: argparse.Namespace) -> int:
    """把既有数据库树迁入共享根，或原地保留并把共享根做成软链接农场。"""
    roots = resolve_roots(args)
    jobs: list[tuple[Path, Path, str]] = []   # (src, dst_root, label)

    if args.db_from:
        jobs.append((Path(args.db_from).expanduser(), roots.db, "DB"))
    if args.plant_from:
        jobs.append((Path(args.plant_from).expanduser(), roots.plant, "PLANT"))
    if not jobs:
        die("至少给一个来源：--db-from ~/database 或 --plant-from ~/plant_virus_db")

    # 先做全部安全检查，任何一条不过就整体不执行
    for src, dst, label in jobs:
        if not src.is_dir():
            die(f"[{label}] 来源不是目录: {src}")
        if src.is_symlink():
            die(f"[{label}] 来源本身已是软链接: {src} → {os.readlink(src)}\n"
                f"     它已经指向别处了，无需再 adopt。")
        try:
            roots.share.resolve().relative_to(src.resolve())
            die(f"[{label}] 共享根 {roots.share} 位于来源 {src} 内部，会造成自我搬迁。")
        except ValueError:
            pass
        try:
            src.resolve().relative_to(roots.share.resolve())
            die(f"[{label}] 来源 {src} 位于共享根内部，无需迁移。")
        except ValueError:
            pass

    move = args.strategy == "move"
    mode = "迁移（移动真实数据，旧路径留软链接）" if move else "原地保留（共享根建软链接指向旧位置）"
    info(f"== adopt ==  策略: {mode}   {'【实际执行】' if args.apply else '【演练，未改动磁盘】'}\n")

    if move and not shutil.which("rsync"):
        warn("未找到 rsync；跨文件系统时将退化为 mv（慢，且中断不便续传）")

    for src, dst, label in jobs:
        units = _top_level_units(src)
        if not units:
            info(f"[{label}] {src} 为空，跳过。")
            continue
        info(f"[{label}] {src} → {dst}   共 {len(units)} 个顶层单元")
        cross = move and not _same_filesystem(src, dst)
        stats = {"move": 0, "link": 0, "skip": 0}
        for u in units:
            _place_unit(u, dst / u.name, move=move, apply=args.apply,
                        cross=cross, indent="  ", stats=stats)

        verb = "迁移" if move else "链接"
        info(f"  → {verb} {stats['move'] + stats['link']}，跳过 {stats['skip']}")

        if move and args.apply:
            leftovers = _top_level_units(src)
            if not leftovers:
                src.rmdir()
                src.symlink_to(dst.resolve())
                info(f"  [链] {src} → {dst}   （旧路径已改为软链接）")
            else:
                warn(f"  {src} 还剩 {len(leftovers)} 个条目未迁移，未改软链接: "
                     + ", ".join(p.name for p in leftovers[:8]))

    if not args.apply:
        info("\n以上为演练。确认无误后加 --apply 执行。")
    else:
        info("\n下一步：")
        info("  1) link --apply     建仓库侧/旧路径兼容软链接")
        info(f"  2) env --write {roots.share}/mmpv-db.env   固化环境变量")
        info("  3) verify           逐项校验")
    return 0


def cmd_link(args: argparse.Namespace) -> int:
    """建兼容软链接：旧路径 → 共享根；共享池 → 仓库侧小文件。"""
    roots = resolve_roots(args)
    entries = load_manifest()
    info(f"== link ==   {'【实际执行】' if args.apply else '【演练，未改动磁盘】'}\n")

    acts: list[tuple[str, Path, Path]] = []   # (说明, link路径, 指向)

    for legacy, dst in ((LEGACY_DB, roots.db), (LEGACY_PLANT, roots.plant)):
        if not dst.exists():
            acts.append((f"跳过（共享根尚无 {dst}）", legacy, dst))
            continue
        acts.append(("旧路径软链接", legacy, dst))

    for e in entries:
        if e.root != "repo":
            continue
        target = roots.resolve(e)
        src = roots.overlay(e)
        if not src.exists():
            acts.append((f"跳过（池中暂无 {src.name}）", target, src))
        elif target.is_symlink():
            acts.append(("已是软链接，跳过", target, src))
        elif target.exists():
            acts.append(("仓库已有真实文件，保留不覆盖", target, src))
        else:
            acts.append(("仓库侧软链接", target, src))

    for note, link_path, dest in acts:
        prefix = note.split("（")[0]
        info(f"  [{prefix}] {link_path}")
        info(f"        → {dest}")
        if not args.apply:
            continue
        if not dest.exists():
            continue
        # 只在 link 位置为空时创建，绝不覆盖真实数据
        if link_path.is_symlink():
            continue
        if link_path.exists():
            warn(f"        {link_path} 已存在且非软链接，保留原样")
            continue
        link_path.parent.mkdir(parents=True, exist_ok=True)
        link_path.symlink_to(dest.resolve())

    if not args.apply:
        info("\n以上为演练。确认无误后加 --apply 执行。")
    return 0


def cmd_publish(args: argparse.Namespace) -> int:
    """把仓库里被 .gitignore 忽略、别人又拿不到的小文件发布到共享池。"""
    roots = resolve_roots(args)
    entries = [e for e in load_manifest() if e.root == "repo"]
    info(f"== publish ==   仓库 → {roots.share / REPO_OVERLAY_DIRNAME}")
    info(f"{'【实际执行】' if args.apply else '【演练，未改动磁盘】'}\n")
    copied = skipped = missing = 0
    for e in entries:
        src = roots.resolve(e)
        dst = roots.overlay(e)
        if not src.exists():
            warn(f"  仓库中缺失: {src}")
            missing += 1
            continue
        if dst.exists():
            info(f"  池中已有，跳过: {dst.name}")
            skipped += 1
            continue
        info(f"  [发] {src.name}  →  {dst}")
        if args.apply:
            dst.parent.mkdir(parents=True, exist_ok=True)
            if src.is_dir():
                shutil.copytree(src, dst)
            else:
                shutil.copy2(src, dst)
        copied += 1
    info(f"\n{'已发布' if args.apply else '将发布'} {copied}，跳过 {skipped}，仓库缺失 {missing}。")
    if missing and not args.apply:
        pass
    info("发布后，其他用户用 `link --apply` 即可把这些文件以软链接形式挂进自己的仓库。")
    return 0


def _entry_env(e: Entry, roots: Roots, args: argparse.Namespace, tmp_dir: Path) -> tuple[dict, Path]:
    dst = roots.resolve(e)
    env = os.environ.copy()
    env.update({
        "DST": str(dst),
        "PARENT": str(dst.parent),
        "DB": str(roots.db),
        "VIRUS_DB": str(roots.db / "virus-db"),
        "PLANT": str(roots.plant),
        "SHARE": str(roots.share),
        "SRC": str(roots.src),
        "REPO": str(roots.repo),
        "TMP": str(tmp_dir),
        "THREADS": str(args.threads),
    })
    return env, dst


def _run_one(e: Entry, steps: list[Step], roots: Roots, args: argparse.Namespace,
             log_dir: Path, tmp_dir: Path, shell: str) -> tuple[str, str]:
    """跑一个条目的全部步骤。返回 (结局, 说明)；结局 ∈ OK|FAILED|SKIP|PLAN。"""
    env, dst = _entry_env(e, roots, args, tmp_dir)
    log = log_dir / f"{e.id}.log"

    for st in steps:
        miss = missing_tools(st.needs)
        if miss:
            return "SKIP", f"缺少工具 {', '.join(miss)}（装好后重跑）"
        if not args.apply:
            info(f"  [待跑] {e.id} step {st.step}: {st.cmd}")
            continue
        script = DL_PREAMBLE + "\nset -euo pipefail\n" + st.cmd
        log_dir.mkdir(parents=True, exist_ok=True)
        with log.open("a", encoding="utf-8") as fh:
            fh.write(f"\n===== {e.id} step {st.step} =====\n{st.cmd}\n")
            fh.flush()
            r = subprocess.run([shell, "-c", script], env=env, cwd=str(roots.share),
                               stdout=fh, stderr=subprocess.STDOUT)
        if r.returncode != 0:
            return "FAILED", f"step {st.step} 退出码 {r.returncode}，看日志 {log}"

    if not args.apply:
        return "PLAN", ""

    # 命令跑完不等于库就位——按清单的哨兵复验一次
    status, _ = probe_entry(e, roots)
    if status == OK:
        return "OK", ""
    return "FAILED", f"命令返回 0 但校验结果仍是「{status}」；日志 {log}"


def resolve_shell() -> str:
    """定位一个**真正可用**的 bash。

    优先 $MMPV_BASH，其次 PATH 上的 bash，最后 /bin/bash。必须实测一遍——
    Windows 上 PATH 里的 bash.exe 可能只是 WSL 安装提示桩，直接用会报一堆乱码。
    """
    cand = os.environ.get("MMPV_BASH") or shutil.which("bash") or "/bin/bash"
    try:
        r = subprocess.run([cand, "-c", 'printf %s "$BASH_VERSION"'],
                           capture_output=True, text=True, timeout=20)
    except (OSError, subprocess.SubprocessError) as exc:
        die(f"找不到可用的 bash（试过 {cand}）: {exc}\n"
            f"     设 MMPV_BASH 指向真实 bash，例如 MMPV_BASH=/bin/bash")
    if r.returncode != 0 or not r.stdout.strip():
        die(f"{cand} 不是可用的 bash（可能只是 WSL 提示桩）。\n"
            f"     设 MMPV_BASH 指向真实 bash，例如 MMPV_BASH=/bin/bash")
    return cand


def cmd_down(args: argparse.Namespace) -> int:
    """按 mmpv_db_fetch.tsv 的配方下载/构建数据库。"""
    roots = resolve_roots(args)
    entries = {e.id: e for e in load_manifest()}
    recipes = load_fetch(Path(args.fetch_file).expanduser() if args.fetch_file else FETCH)

    if args.only:
        unknown = [i for i in args.only if i not in entries]
        if unknown:
            die(f"未知 id: {', '.join(unknown)}（用 list 查看可用 id）")
        selected = [entries[i] for i in args.only]
    else:
        selected = list(entries.values())
    if args.kind:
        selected = [e for e in selected if e.kind in args.kind]

    runnable = [e for e in selected if e.id in recipes]
    manual = [e for e in selected if e.id not in recipes]

    info("== down ==")
    info(f"共享库根 : {roots.share}")
    info(f"下载临时区: {roots.share / 'tmp'}   （大件先落这里，支持断点续传）")
    info(f"日志      : {roots.share / 'logs' / 'down'}")
    info(f"模式      : {'【实际下载/构建】' if args.apply else '【演练，不改磁盘】'}"
         f"   线程 {args.threads}   并行 {args.jobs}")
    info(f"配方覆盖  : {len(runnable)}/{len(selected)} 项"
         f"（全清单 {len(entries)} 项，其中 {len(recipes)} 项有配方）")
    shell = resolve_shell() if args.apply else "（演练不解析）"
    info(f"shell     : {shell}\n")

    if not runnable:
        if manual:
            info("所选条目都没有获取配方，见下方清单。")
        return _down_report([], manual, roots, entries, recipes)

    if args.apply:
        (roots.share / "tmp").mkdir(parents=True, exist_ok=True)
        (roots.share / "logs" / "down").mkdir(parents=True, exist_ok=True)

    results: list[tuple[str, str, str]] = []   # (id, 结局, 说明)

    def work(e: Entry) -> tuple[str, str, str]:
        # 已就位就跳过（--force 可强制重跑）
        status, path = probe_entry(e, roots)
        if status == OK and not args.force:
            return e.id, "SKIP", f"已就位（{path}）"
        # 依赖先决条目
        missing_deps = [d for d in {d for s in recipes[e.id] for d in s.deps}
                        if probe_entry(entries[d], roots)[0] != OK]
        if missing_deps:
            return e.id, "SKIP", f"先决条目未就位: {', '.join(missing_deps)}"
        info(f"\n── {e.id}  ({e.kind}, {e.size}) → {path}")
        if status != OK:
            info(f"   当前状态: {status}")
        return (e.id, *_run_one(e, recipes[e.id], roots, args,
                                roots.share / "logs" / "down", roots.share / "tmp", shell))

    if args.jobs > 1 and args.apply:
        from concurrent.futures import ThreadPoolExecutor
        with ThreadPoolExecutor(max_workers=args.jobs) as ex:
            for res in ex.map(work, runnable):
                results.append(res)
    else:
        for e in runnable:
            results.append(work(e))

    return _down_report(results, manual, roots, entries, recipes, dry=not args.apply)


def _down_report(results: list[tuple[str, str, str]], manual: list[Entry], roots: Roots,
                 entries: dict[str, Entry], recipes: dict[str, list[Step]],
                 dry: bool = False) -> int:
    if results:
        info("\n== 结果 ==")
        rows = [[r[1], r[0], r[2]] for r in results]
        info(table(rows, ["结局", "id", "说明"]))
        tally: dict[str, int] = {}
        for _, verdict, _ in results:
            tally[verdict] = tally.get(verdict, 0) + 1
        info("  " + "  ".join(f"{k}:{v}" for k, v in sorted(tally.items())))

    # kind=repo 是随 git 自带的（clone 就有），不该混进「需要共享池」的清单里
    shipped = [e for e in manual if e.kind == "repo"]
    manual = [e for e in manual if e.kind != "repo"]
    if shipped:
        info("\n== 随 git 自带，无需部署 ==")
        for e in shipped:
            info(f"      {pad(e.id, 30)} {e.path}")

    if manual:
        info("\n== 无获取配方（仓库内无可靠复现途径）==")
        by_kind: dict[str, list[Entry]] = {}
        for e in manual:
            by_kind.setdefault(e.kind, []).append(e)
        advice = {
            "local": "无法重下 → 从已有机器/共享池 rsync，或 adopt --strategy link 原地引用",
            "build": "建库脚本不在仓库（见 utils/db_build/README.md 待收编清单）→ 从共享池拷贝",
            "download": "上游命令随版本变化 → 见 DATABASE_SETUP.md 对应小节",
        }
        for kind in ("local", "build", "download"):
            group = by_kind.get(kind, [])
            if not group:
                continue
            info(f"\n  [{kind}] {advice.get(kind, '')}")
            for e in group:
                info(f"      {pad(e.id, 30)} {e.path}")
        info(f"\n  逐库命令见 DATABASE_SETUP.md §1–§6；架构与替代获取方式见 doc/DB_DEPLOYMENT.md")

    if dry and results:
        info("\n以上为演练。确认无误后加 --apply 执行。")

    failed = [r[0] for r in results if r[1] == "FAILED"]
    if failed:
        info(f"\n失败 {len(failed)} 项: {', '.join(failed)}")
        return 1
    return 0


def cmd_fetch(args: argparse.Namespace) -> int:
    """查看若干 id 的获取方式：有配方就打印配方命令，没配方说明该走哪条路。"""
    roots = resolve_roots(args)
    entries = {e.id: e for e in load_manifest()}
    recipes = load_fetch(Path(args.fetch_file).expanduser() if args.fetch_file else FETCH)
    rc = 0
    for eid in args.id:
        e = entries.get(eid)
        if e is None:
            warn(f"未知 id: {eid}")
            rc = 1
            continue
        dst = roots.resolve(e)
        status, _ = probe_entry(e, roots)
        info(f"\n── {eid}  ({e.kind}, {e.size}, {e.required_mark})  当前: {status} ──")
        info(f"   目标: {dst}")
        info(f"   说明: {e.note}")
        if eid in recipes:
            info("   配方（可直接跑）:")
            for st in recipes[eid]:
                need = f"  需 {','.join(st.needs)}" if st.needs else ""
                dep = f"  先决 {','.join(st.deps)}" if st.deps else ""
                info(f"     step {st.step}{need}{dep}")
                info(f"       {st.cmd}")
            info(f"   执行: mmpv_db.py down --only {eid} --apply")
        elif e.kind == "local":
            info("   无配方: 该类库无法重下。从已有机器/共享池 rsync，或原地引用：")
            info("           mmpv_db.py adopt --strategy link --db-from <含此库的目录> --apply")
        elif e.kind == "repo":
            info("   无配方: 被 .gitignore 忽略的小文件，从共享池取用：")
            info("           mmpv_db.py link --apply   （维护者需先 publish --apply）")
        else:
            info("   无配方: 建库脚本不在仓库（见 virome_discovery_pipeline/utils/db_build/README.md）")
            info("           或上游命令随版本变化 → 见 DATABASE_SETUP.md 对应小节，或从共享池拷贝")
    return rc


# ─────────────────────────── 入口 ───────────────────────────

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="mmpv_db.py",
        description="MMPV 数据库部署 / 迁移 / 校验工具",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__.split("用法速览")[0].strip(),
    )
    p.add_argument("--share", help="共享库根（默认 $MMPV_DATA_SHARE 或 ~/mmpv-db）")
    p.add_argument("--repo", help="仓库根（默认由本脚本位置推导）")
    p.add_argument("--db-root", help="覆盖 ${MMPV_DB_ROOT}（默认 <share>/database）")
    p.add_argument("--plant-root", help="覆盖 ${MMPV_PLANT_VIRUS_DB}（默认 <share>/plant_virus_db）")
    p.add_argument("--src-root", help="上游原始数据目录（默认 $MMPV_SRC 或 <share>/src）")
    p.add_argument("--fetch-file", metavar="TSV",
                   help="用自定义获取配方替换 mmpv_db_fetch.tsv（站点自维护/测试用）")
    sub = p.add_subparsers(dest="cmd", required=True)

    sub.add_parser("list", help="打印依赖清单").set_defaults(func=cmd_list)

    sp = sub.add_parser("which", help="打印若干 id 的解析路径与状态")
    sp.add_argument("id", nargs="+")
    sp.set_defaults(func=cmd_which)

    sp = sub.add_parser("env", help="打印/写出环境变量块")
    sp.add_argument("--write", metavar="FILE", help="写到文件而不是打到 stdout")
    sp.set_defaults(func=cmd_env)

    sub.add_parser("check", help="部署前体检：可写性、旧路径现状、缺项统计").set_defaults(func=cmd_check)

    sp = sub.add_parser("verify", help="逐项校验数据库是否就位")
    sp.add_argument("--only", nargs="+", help="只校验这些 id")
    sp.set_defaults(func=cmd_verify)

    sp = sub.add_parser("layout", help="建规范目录骨架")
    sp.add_argument("--apply", action="store_true", help="真正落盘（默认演练）")
    sp.set_defaults(func=cmd_layout)

    sp = sub.add_parser("adopt", help="把既有数据库树迁入共享根")
    sp.add_argument("--db-from", metavar="DIR", help="既有数据库根，如 ~/database")
    sp.add_argument("--plant-from", metavar="DIR", help="既有植物库根，如 ~/plant_virus_db")
    sp.add_argument("--strategy", choices=["move", "link"], default="move",
                    help="move=搬真实数据并给旧路径留软链接（推荐）；link=数据不动，共享根建软链接")
    sp.add_argument("--apply", action="store_true", help="真正执行（默认演练）")
    sp.set_defaults(func=cmd_adopt)

    sp = sub.add_parser("link", help="建旧路径与仓库侧的兼容软链接")
    sp.add_argument("--apply", action="store_true", help="真正执行（默认演练）")
    sp.set_defaults(func=cmd_link)

    sp = sub.add_parser("publish", help="把仓库侧小文件发布到共享池")
    sp.add_argument("--apply", action="store_true", help="真正执行（默认演练）")
    sp.set_defaults(func=cmd_publish)

    sp = sub.add_parser("down", help="按配方下载/构建数据库（换服务器时用）")
    sp.add_argument("--only", nargs="+", help="只处理这些 id")
    sp.add_argument("--kind", nargs="+", choices=["download", "build", "local", "repo"],
                    help="只处理这些类别的条目")
    sp.add_argument("--threads", type=int, default=60, help="建索引线程数（默认 60）")
    sp.add_argument("--jobs", type=int, default=1, help="并行条目数（默认 1）")
    sp.add_argument("--force", action="store_true", help="已就位也重跑")
    sp.add_argument("--apply", action="store_true", help="真正下载/构建（默认演练）")
    sp.set_defaults(func=cmd_down)

    sp = sub.add_parser("fetch", help="查看某个库的获取方式与配方")
    sp.add_argument("id", nargs="+")
    sp.set_defaults(func=cmd_fetch)

    return p


def main() -> int:
    args = build_parser().parse_args()
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
