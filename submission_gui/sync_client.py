#!/usr/bin/env python3
# sync_client.py v1.0 — 一键同步: 本地 CSV → 服务器 → 重建 .sqn → 拉回本地
# 用法: python sync_client.py <dataset> <local_csv>
# 例:   python sync_client.py barbarum ./edited/unified_metadata.csv
import subprocess
import sys
from pathlib import Path

HOST = 'zhangwenda@202.119.189.246'
REMOTE_DIR = '/home/zhangwenda/MMPV-RNA/virome_submission_pipeline'
HERE = Path(__file__).resolve().parent
LOCAL_SQN_DIR = HERE.parent / 'submission_gui' / 'sqn_sync' / '{name}'


def run(cmd, **kw):
    r = subprocess.run(cmd, shell=True, capture_output=True, text=True, **kw)
    if r.returncode != 0:
        print(f'[FAIL] {cmd[:120]}\n{r.stderr[-1500:]}')
        sys.exit(1)
    return r


def ssh_quote(s):
    return "'" + s.replace("'", "'\\''") + "'"


def main():
    if len(sys.argv) < 3:
        print(__doc__)
        sys.exit(1)
    name, local_csv = sys.argv[1], Path(sys.argv[2])
    if not local_csv.exists():
        print(f'[FAIL] 找不到 {local_csv}')
        sys.exit(1)

    # 部署远端脚本 (幂等)
    run(f'scp "{HERE / "sync_sqn_from_csv.py"}" {HOST}:/tmp/sync_sqn_from_csv.py')

    # 上传 CSV
    remote_csv = f'/tmp/submission_csv_{name}.csv'
    run(f'scp {ssh_quote(str(local_csv))} {HOST}:{remote_csv}')

    # 服务器执行
    r = run(f'ssh -o BatchMode=yes {HOST} '
            f'"python3 {REMOTE_DIR}/sync_sqn_from_csv.py '
            f'{name} {remote_csv} 2>&1 | tail -5"')
    print(r.stdout.strip())

    # 拉回 .sqn + source_filled.src (审计用)
    out_dir = Path(str(LOCAL_SQN_DIR).format(name=name))
    out_dir.mkdir(parents=True, exist_ok=True)
    sqn_remote = {
        'barbarum': '/home/zhangwenda/virus/data-2026/data-test/RNA-Lycium_barbarum_out',
        'ruthenicum': '/home/zhangwenda/virus/data-2026/data-test/RNA-Lycium_ruthenicum_out',
        'chinense': '/home/zhangwenda/virus/data-2026/data-test/RNA-Lycium_chinense_out',
        'amarum': '/home/zhangwenda/virus/data-2026/data-test/RNA-Lycium_amarum_out',
        'Fusarium': '/home/zhangwenda/virus/data-2026/data-test/RNA-Fusarium_nematophilum_out',
        'Alternaria': '/home/zhangwenda/virus/data-2026/data-test/RNA-Alternaria_alternata_out',
        'Aphis': '/home/zhangwenda/virus/data-2026/data-test/RNA-Aphis_gossypii_out',
        'mix': '/home/zhangwenda/virus/ningxiagouqi/11.merge_assembly/mix/out',
    }[name]
    sub = f'{sqn_remote}/submission_{name}_virome'
    run(f'scp {HOST}:{sub}/suvtk_submission/submission.sqn {ssh_quote(str(out_dir))}')
    run(f'scp {HOST}:{sub}/id_norm/source_filled.src {ssh_quote(str(out_dir))}')
    print(f'[OK] 已拉回: {out_dir}\\submission.sqn')


if __name__ == '__main__':
    main()
