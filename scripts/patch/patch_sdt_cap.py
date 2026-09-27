#!/usr/bin/env python3
"""
修 SDT 属级矩阵的规模爆炸问题:

sdt_genus_matrix.py → build_mat_sdt_exact 是 O(N^2) 全对全精确比对
(每对跑一次外部 MAFFT L-INS-i)。属内序列数一大就不可行:
  Biavirus 519 条 → 134,421 对; Metavirus 227 条 → 25,651 对
  两属合计占全部 178,798 对的 90%, 单属跑一天都跑不完。

实测分布: 120 个属 <=20 条 (总 1664 对), 只有 2 个属 >100 条。
→ 加 sdt_max_seqs 上限 (默认 80), 超限的属跳过并记日志。
   该阈值覆盖 72 个属 / 18,714 对, 排除只有 2 个的巨型属。
"""
import re
import shutil
import sys
from pathlib import Path

PIPE = Path("/home/zhangwenda/MMPV-RNA/virome_discovery_pipeline/virome_pipeline.py")
BAK = PIPE.with_suffix(".py.bak_sdtcap_20260910")

if not PIPE.is_file():
    sys.exit(f"找不到: {PIPE}")

src = PIPE.read_text(encoding="utf-8")
orig = src
applied = []

# ① 循环里加规模上限
old = """        sdt_th = getattr(self.args, 'sdt_threads', 8)
        n_run = 0
        for g, cids in sorted(gen_of.items(), key=lambda x: -len(x[1])):
            user_ok = [c for c in cids if ratio_map.get(c, 1.0) >= min_ratio]
            if len(user_ok) < 2:
                continue"""
new = """        sdt_th = getattr(self.args, 'sdt_threads', 8)
        # SDT 是 O(N^2) 全对全精确比对 (每对一次外部 MAFFT), 属内序列过多时
        # 组合数爆炸 (519 条 → 13.4 万对), 必须设上限。超限属跳过并记日志。
        sdt_max = getattr(self.args, 'sdt_max_seqs', 80)
        n_run = 0
        n_skip_big = 0
        for g, cids in sorted(gen_of.items(), key=lambda x: -len(x[1])):
            user_ok = [c for c in cids if ratio_map.get(c, 1.0) >= min_ratio]
            if len(user_ok) < 2:
                continue
            if len(user_ok) > sdt_max:
                n_skip_big += 1
                self.log.info('  [9b] SDT %s 跳过: %d 条 > 上限 %d (O(N^2) 组合 %d 对)',
                              g, len(user_ok), sdt_max, len(user_ok)*(len(user_ok)-1)//2)
                continue"""
if old in src:
    src = src.replace(old, new, 1)
    applied.append("(1) SDT 加 sdt_max_seqs 上限 (默认 80)")
elif "sdt_max_seqs" in src:
    applied.append("(1) 已修复")
else:
    sys.exit("锚点① 未匹配")

# ② 汇总日志补跳过数
old2 = """        self.log.info('  [9b] SDT 矩阵: %d 个属 → %s', n_run, sdt_out)"""
new2 = """        self.log.info('  [9b] SDT 矩阵: %d 个属 → %s (超上限跳过 %d 个)',
                      n_run, sdt_out, n_skip_big)"""
if old2 in src:
    src = src.replace(old2, new2, 1)
    applied.append("(2) 汇总日志补超限跳过计数")
elif "超上限跳过" in src:
    applied.append("(2) 已修复")
else:
    sys.exit("锚点② 未匹配")

# ③ 命令行参数 (锚点: g.add_argument('--sdt-threads' ...) 行)
old3 = "    g.add_argument('--sdt-threads', type=int, default=8, help='SDT 逐对比对线程数')"
new3 = ("    g.add_argument('--sdt-threads', type=int, default=8, help='SDT 逐对比对线程数')\n"
        "    g.add_argument('--sdt-max-seqs', type=int, default=80,\n"
        "                   help='SDT 属级矩阵单属序列上限 (超过则跳过, 防 O(N^2) 组合爆炸)')")
if old3 in src:
    src = src.replace(old3, new3, 1)
    applied.append("(3) 新增 --sdt-max-seqs 参数 (默认 80)")
elif "--sdt-max-seqs" in src:
    applied.append("(3) 已修复")
else:
    applied.append("(3) 未找到 sdt-threads 参数定义 (跳过)")

if src != orig:
    if not BAK.exists():
        shutil.copy2(PIPE, BAK)
        print(f"备份: {BAK}")
    PIPE.write_text(src, encoding="utf-8")
    print("[virome_pipeline.py]")
    for a in applied:
        print("  -", a)
else:
    print("无改动")
