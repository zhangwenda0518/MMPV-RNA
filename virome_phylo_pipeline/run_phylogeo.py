#!/usr/bin/env python3
"""统一 BEAST 系统地理学执行入口 (MMPV-RNA virome_phylo_pipeline)

用法:
  python run_phylogeo.py --virus GCVA --prior skyline --chains 5 --threads 8
  python run_phylogeo.py --virus PSTVd --prior bdsky --chains 3 --chain-length 20000000

流程: 数据定位 → XML 生成 → 多链提交 → 状态记录(供监控) → 合并 → MCC → 参数解析
已知数据集 (GCVA/PSTVd) 自动定位; 其他用 --fasta/--metadata 指定

管线整合: submit_phylogeo_chains() 供 phylo_pipeline 的 phylogeo stage 复用
(管线生成 XML 后默认多链提交, 同一份提交逻辑, 不重复实现)。
"""
import argparse, csv, json, os, re, shlex, subprocess, sys, time
from pathlib import Path
from typing import Dict, List, Optional

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'utils'))

import logging
logging.basicConfig(level=logging.INFO, format='%(message)s')
logger = logging.getLogger()
from utils.virphy_bridge import LogCollector
from utils.beast1_bridge import generate_beast1_phylogeo_xml
from utils.dataset_config import (dataset as get_dataset, env as get_env,
                                  known_dataset_names)

# BDSKY 用独立工作目录 (避免与 skyline 结果混淆)
PRIOR_DIR = {'skyline': 'phylogeography_skyline', 'bdsky': 'phylogeography_bdsky'}
# 全链路统一 burnin 口径 (2026-08-27 定案): 10%。
# 历史坑: 合并用 1% / 诊断图用 10%, 同一条链两套口径。统一到保守值 10%,
# 与 Tracer 惯例一致; 下游 auto_parse_params 的 BURNIN=0 (吃已烧过的 merged.log) 语义不变。
MERGE_BURNIN_FRAC = 0.10
# 需要按链改名的产物扩展名 (log/trees/rates log/ops)
CHAIN_EXTS = ('.log', '.trees', '.ops.txt')


def detect_chain_tokens(xml_text: str):
    """从 XML 自动检测产物文件名 token (log/trees/ops), 支持 beast1/beast2 命名。
    返回基础名列表, 如 ['phylogeo_beast1.log', 'phylogeo_beast1.trees']
    或 ['phylogeo.log', 'phylogeo.trees']。"""
    fnames = re.findall(r'fileName="([^"]+)"', xml_text)
    # 历史坑: operatorAnalysis 用专属属性而非 fileName, 不加则所有链的 ops 写同一文件互相覆盖
    fnames += re.findall(r'operatorAnalysis="([^"]+)"', xml_text)
    toks = sorted(set(
        os.path.basename(f) for f in fnames if f.endswith(CHAIN_EXTS)
    ))
    return toks


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def _pid_alive(pid) -> bool:
    """进程是否存活 (POSIX os.kill(pid,0)); Windows / 权限异常时返回 True (不误判为死)。"""
    try:
        pid = int(pid)
    except (TypeError, ValueError):
        return True
    if pid <= 0:
        return True
    if os.name == 'nt':
        return True          # 目标机为 Linux; 本机不适用时不臆断
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    except OSError:
        return True
    return True


def detect_xml_log_every(xml_path: str) -> Optional[int]:
    """从 XML 提取 logEvery (state 记录间隔), 供 merge_results 精确判定链是否跑完。
    取所有 logEvery 中的最大值 (主 log/tree log 的间隔通常最大)。"""
    try:
        with open(xml_path, encoding='utf-8', errors='replace') as f:
            txt = f.read()
    except OSError:
        return None
    vals = [int(v) for v in re.findall(r'logEvery="(\d+)"', txt)]
    return max(vals) if vals else None


def resolve_paths(virus, fasta, metadata):
    """定位数据文件：优先已知数据集注册表；否则要求 --fasta/--metadata。"""
    known_key = virus.upper()
    k = get_dataset(known_key)
    if k and not (fasta and metadata):
        return (f'{k["dir"]}/{k["fasta"]}', f'{k["dir"]}/{k["metadata"]}',
                f'{k["dir"]}/{k["work"]}')
    if not (fasta and metadata):
        avail = ', '.join(known_dataset_names()) or '(无配置)'
        raise SystemExit(f"需提供 --fasta 和 --metadata (或使用已知数据集: {avail})")
    return (fasta, metadata, os.path.join(os.path.dirname(fasta), f'phylogeography_{virus}'))


def submit_phylogeo_chains(
    xml_path: str,
    work_dir: str,
    chains: int = 3,
    threads: int = 8,
    prefix: str = "phylogeo_sky",
    chain_length: int = 10_000_000,
    beast_jlp: Optional[str] = None,
    beast_cp: Optional[str] = None,
    logger=None,
    label: str = "phylogeo",
    seed_base: Optional[int] = None,
) -> List[Dict]:
    """多链 BEAST 提交 (管线 phylogeo/beast 全长 stage 复用)

    对每条链: 生成 {prefix}{i}.xml (替换 log/trees 文件名 token) →
    nohup java dr.app.beast.BeastMain -seed {seed} 后台运行。
    返回 [{prefix, seed}, ...]。label 仅用于日志前缀。

    可复现性 (2026-09-15 修复): 旧版 seed = i*1024 + int(time.time()*1000)%100000,
    种子随墙上时间变化 → 同一输入两次运行结果必然不同。现改为确定性种子
    seed_base + i*1024, seed_base 默认取 PHYLO_SEED_BASE 环境变量, 否则用固定常量。
    已落盘的 seed 会写进 run_status.json, 便于复现与追溯。
    """
    log = logger or globals()['log']
    # 统一日志接口: callable 函数 / logging.Logger(.info) / LogCollector(.emit)
    if callable(log):
        emit = log
    elif hasattr(log, 'info'):
        emit = log.info
    else:
        emit = log.emit
    beast_jlp = beast_jlp or get_env("beast_jlp") or os.environ.get("beast_jlp", "")
    beast_cp = beast_cp or get_env("beast_cp") or ""
    # 确定性种子基准 (2026-09-15 修复): 保证同输入可复现
    if seed_base is None:
        try:
            seed_base = int(os.environ.get("PHYLO_SEED_BASE", "") or 20260915)
        except ValueError:
            seed_base = 20260915
    if not os.path.exists(xml_path):
        emit(f"  [{label}] 无 XML: {xml_path}, 跳过提交")
        return []
    with open(xml_path, encoding='utf-8', errors='replace') as _f:
        base = _f.read()
    tokens = detect_chain_tokens(base)
    # BEAST1 vs BEAST2 自动检测: BEAST2 XML 是 <beast version="2.x"...> 语法,
    # BEAST1 运行器 dr.app.beast.BeastMain 无法解析; 反之 BEAST2 launcher 不吃 BEAST1 XML。
    is_beast2 = bool(re.search(r'<beast\s+version="2', base))
    emit(f"  [{label}] 检测到产物 token: {tokens}"
         f" (BEAST{'2' if is_beast2 else '1'} XML)")
    emit(f"  [{label}] 提交 {chains} 链 × {threads} 线程...")
    launched = []
    for i in range(1, chains + 1):
        p = f'{prefix}{i}' if i > 1 else f'{prefix}1'
        xml = base
        # 每个 token 的主干换成链前缀 (phylogeo_beast1.log → phylogeo_sky2.log)
        # 注意: 只替换主前缀 (第一个 '.' 之前), 保留子 token 后缀
        # (如 phylogeo_beast1.location.rates.log → phylogeo_sky2.location.rates.log)。
        # 历史坑: 用 rsplit('.',1)[0] 取 stem 会把 .location.rates 也吃掉,
        # 导致 rateMatrixLog 与主 log 撞名, 主 log 被 location 列覆盖 (缺核心参数)。
        for tok in tokens:
            stem = tok.rsplit('.', 1)[0]
            first_dot = stem.find('.')
            # 注意: 不能用 base 命名 (外层 base = 原始 XML 文本, 覆盖会导致后续链拿到残片)
            tok_base = stem[:first_dot] if first_dot > 0 else stem
            new_name = tok.replace(tok_base, p)
            xml = xml.replace(tok, new_name)
            xml = xml.replace(os.path.basename(tok), new_name)
        chain_xml = os.path.join(work_dir, f'{p}.xml')
        with open(chain_xml, 'w', encoding='utf-8') as f:
            f.write(xml)
        # 确定性种子 (2026-09-15 修复): 旧版混入 int(time.time()*1000)%100000 → 不可复现
        seed = int(seed_base) + i * 1024
        # shell 拼接的路径统一加引号 (2026-09-15 修复): 路径含空格即被拆成多个参数
        _wd = shlex.quote(work_dir)
        _jlp = shlex.quote(beast_jlp)
        _cp = shlex.quote(beast_cp)
        _xmlname = shlex.quote(f'{p}.xml')
        _runlog = shlex.quote(f'{p}_run.log')
        if is_beast2:
            # BEAST2: 用发行版 launcher (mambaforge/bin/beast 等), 自动带全套 package classpath。
            # 注意: 服务器默认 PATH 里的 beast 是 BEAST1.10.4, 不是 BEAST2!
            # 无 beast2_bin 显式配置时直接报错, 避免 BEAST1 跑 BEAST2 XML 静默崩 (SEVERE)。
            beast2_bin = get_env("beast2_bin") or os.environ.get("beast2_bin") or ""
            if not beast2_bin:
                emit(f"  [{label}] ⚠️ 检测到 BEAST2 XML ({p}.xml), 但未配置 beast2_bin: "
                     f"服务器 beast 命令是 BEAST1.10.4, 无法运行 BEAST2 XML。"
                     f"请改用 BEAST1 流程 (config beast.version: 1, 植物病毒文献标准) "
                     f"或在 datasets.yaml env 配置 beast2_bin。跳过该链。")
                continue
            cmd = (f'cd {_wd} && nohup {shlex.quote(beast2_bin)} -overwrite '
                   f'-threads {int(threads)} -seed {seed} {_xmlname} '
                   f'> {_runlog} 2>&1 & echo $!')
        else:
            cmd = (f'cd {_wd} && nohup java -Xms64m -Xmx2048m '
                   f'-Djava.library.path={_jlp} -cp {_cp} '
                   f'dr.app.beast.BeastMain -overwrite -threads {int(threads)} '
                   f'-seed {seed} {_xmlname} > {_runlog} 2>&1 & echo $!')
        # 2026-09-15 修复: 旧版 subprocess.run(cmd, shell=True) 不查 returncode,
        # 只 sleep 3 不验证链是否真的起来 → 提交失败被静默当成成功。
        try:
            _r = subprocess.run(cmd, shell=True, capture_output=True, text=True, encoding='utf-8', errors='replace', timeout=60)
        except subprocess.TimeoutExpired:
            emit(f"  [{label}] {p}: 提交命令超时 (60s)")
            continue
        _pid = None
        _out = (_r.stdout or '').strip().splitlines()
        if _out and _out[-1].strip().isdigit():
            _pid = int(_out[-1].strip())
        if _r.returncode != 0:
            emit(f"  [{label}] {p}: 提交失败 (exit={_r.returncode}): "
                 f"{(_r.stderr or '').strip()[-200:]}")
            continue
        if _pid is not None and not _pid_alive(_pid):
            emit(f"  [{label}] {p}: 已提交但进程 {_pid} 未存活, 请检查 {p}_run.log")
        launched.append({'prefix': p, 'seed': seed, 'pid': _pid})
        time.sleep(3)
    emit(f"  [{label}] 共 {len(launched)} 链已提交 → {work_dir}")

    # 写 run_status.json (merge_results 的驱动; 链完成后 logcombiner+treeannotator)
    if launched:
        _le = detect_xml_log_every(xml_path)
        status = {
            'label': label, 'work_dir': work_dir,
            'chains': chains, 'chain_length': chain_length,
            'launched_at': time.strftime('%Y-%m-%d %H:%M:%S'),
            'prefixes': [l['prefix'] for l in launched],
            'seeds': [l['seed'] for l in launched],
            # 2026-09-15 修复: 记录 PID 与 log_every。
            # · PID → merge_results.all_done 只需检查本任务的进程, 不再依赖
            #   "全机器 java beast 进程数为 0"(并行跑第二个数据集时永不 done)。
            # · log_every → 精确计算链的最后一个 state (chain_length 非 logEvery
            #   倍数时, 旧版按 5000 猜会算错 target)。
            'pids': [l.get('pid') for l in launched],
            'log_every': _le,
            'merge_burnin': int(chain_length * MERGE_BURNIN_FRAC),
        }
        try:
            with open(os.path.join(work_dir, 'run_status.json'), 'w', encoding='utf-8') as _sf:
                json.dump(status, _sf, indent=2)
        except Exception as e:
            emit(f"  [{label}] 写 run_status.json 失败: {e}")
    return launched


# 通用别名: beast 全长定年等多链场景复用同一份提交逻辑
def submit_mcmc_chains(*args, **kwargs):
    kwargs.setdefault("label", "beast")
    return submit_phylogeo_chains(*args, **kwargs)


def main():
    ap = argparse.ArgumentParser(description='统一 BEAST 系统地理学入口')
    ap.add_argument('--virus', default='GCVA', help='病毒名 (GCVA/PSTVd 或自定义)')
    ap.add_argument('--prior', default='skyline', choices=['skyline', 'constant', 'bdsky'],
                    help='树先验 (默认 skyline, 文献标准)')
    ap.add_argument('--chains', type=int, default=5, help='并行链数 (默认 5)')
    ap.add_argument('--threads', type=int, default=8, help='每链线程数 (默认 8)')
    ap.add_argument('--chain-length', type=int, default=50_000_000, help='每链 state 数')
    ap.add_argument('--fasta', default=None, help='自定义比对文件')
    ap.add_argument('--metadata', default=None, help='自定义元数据 CSV')
    ap.add_argument('--work-dir', default=None, help='自定义工作目录 (默认自动)')
    ap.add_argument('--submit-only', action='store_true', help='只提交, 不生成 XML')
    args = ap.parse_args()

    virus = args.virus
    fasta, metadata, work_dir = resolve_paths(virus, args.fasta, args.metadata)
    if args.work_dir:
        work_dir = args.work_dir
    if args.prior == 'bdsky':
        work_dir = work_dir.replace('phylogeography_skyline', 'phylogeography_bdsky')
    os.makedirs(work_dir, exist_ok=True)

    # ── 1. 杀残留 + 清理 (独立运行专用; 管线内不触发) ──
    prefix = f'{virus.lower()}_sky'
    if args.prior == 'bdsky':
        prefix = f'{virus.lower()}_bd'
    os.system(f'pkill -f "beast.*{prefix}" 2>/dev/null')
    time.sleep(3)
    for f in os.listdir(work_dir):
        if f.endswith('.xml') or f.endswith('.log') or f.endswith('.trees') or f.endswith('.ops.txt'):
            try:
                os.remove(f'{work_dir}/{f}')
            except OSError:
                pass

    # ── 2. 生成 XML ──
    xml_path = os.path.join(work_dir, 'phylogeo_beast1.xml')
    if not args.submit_only:
        log(f"生成 XML: {virus} / {args.prior} / {args.chain_length:,} states")
        xml_r = generate_beast1_phylogeo_xml(
            fasta_file=fasta, metadata_csv=metadata,
            output_dir=work_dir, clock_model='ucln',
            tree_prior=args.prior, discretize_locations=True,
            substitution_model='auto', chain_length=args.chain_length,
            log=LogCollector(logger))
        if not xml_r['success']:
            raise SystemExit(f"XML 生成失败: {xml_r.get('error')}")
        log(f"XML: {xml_path} ({os.path.getsize(xml_path)/1e6:.1f}MB), "
            f"{xml_r.get('n_taxa')} taxa, {xml_r.get('n_locations')} 地点")

    # ── 3. 提交多链 (复用 submit_phylogeo_chains) ──
    launched = submit_phylogeo_chains(
        xml_path=xml_path, work_dir=work_dir, chains=args.chains,
        threads=args.threads, prefix=prefix, chain_length=args.chain_length,
        logger=log)

    # ── 4. 状态文件 (供监控) ──
    status = {
        'virus': virus, 'prior': args.prior, 'work_dir': work_dir,
        'chains': args.chains, 'chain_length': args.chain_length,
        'launched_at': time.strftime('%Y-%m-%d %H:%M:%S'),
        'prefixes': [l['prefix'] for l in launched],
        'seeds': [l['seed'] for l in launched],
        'merge_burnin': int(args.chain_length * MERGE_BURNIN_FRAC),
    }
    status_path = os.path.join(work_dir, 'run_status.json')
    # 2026-09-15: 句柄泄漏 + 缺 encoding (与 merge_results 读侧口径对齐)
    with open(status_path, 'w', encoding='utf-8') as _sp:
        json.dump(status, _sp, indent=2)
    log(f"状态文件: {status_path}")

    # ── 5. 等待启动并验证 ──
    time.sleep(120)
    # 2026-09-15: ① prefix 未加引号, 含空格即被拆成多个参数;
    # ② 非 POSIX 环境 (Windows) 没有 ps -> 原来直接 FileNotFoundError 崩掉;
    # ③ 探测失败不应影响已经提交的链, 故只告警不抛。
    try:
        _ps_cmd = (f'ps aux | grep "[b]east" | grep {shlex.quote(prefix)} '
                   f'| grep java | wc -l')
        _r = subprocess.run(_ps_cmd, shell=True, capture_output=True, text=True,
                            encoding='utf-8', errors='replace', timeout=60)
        _n_proc = (_r.stdout or '').strip()
    except Exception as _e:
        _n_proc = '?'
        log(f"⚠️ 进程计数失败 ({type(_e).__name__}: {_e}); 以 log 文件判定为准")
    log(f"进程: {_n_proc}/{args.chains}")
    for p in launched:
        rl = f'{work_dir}/{p["prefix"]}_run.log'
        if os.path.exists(rl):
            c = open(rl, encoding='utf-8', errors='replace').read()
            n_sev = c.count('SEVERE')
            if n_sev:
                for line in c.split('\n'):
                    if 'SEVERE' in line:
                        log(f"  ⚠️ {p['prefix']}: {line.strip()[:150]}")
            else:
                log(f"  ✅ {p['prefix']}: 启动正常")


if __name__ == '__main__':
    main()
