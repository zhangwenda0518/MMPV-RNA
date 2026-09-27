#!/usr/bin/env python3
"""
修两个 stage 缺陷:

缺陷 1: no_reads_stages 缺少 'analysis_verify'
  -> `--stage analysis analysis_verify report` 触发 __init__ 里的 scan_samples_in_dir,
     对整个 output_dir 做 12 次 rglob (含 08_Rescue 8000+ SPAdes 临时目录), 卡死数小时。
     而 analysis_verify / report 根本不需要样本 reads。

缺陷 2: Part 4 rescue_detection 的 reads_dir 只判断目录存在, 不判断是否含序列文件
  -> 00a_CleanData 只剩日志 (序列已清理), 仍被优先选中, batch_virus_depth 报"未找到样本文件"。
     改为"目录存在 且 含序列文件"才选用, 否则回退 hostdep -> raw。
"""
import re
import shutil
import sys
from pathlib import Path

TARGET = Path("/home/zhangwenda/MMPV-RNA/virome_discovery_pipeline/virome_pipeline.py")
BAK = TARGET.with_suffix(".py.bak_stagefix_20260902")

if not TARGET.is_file():
    sys.exit(f"找不到目标文件: {TARGET}")

src = TARGET.read_text(encoding="utf-8")
orig = src
changes = []

# ── 缺陷 1: 白名单补 analysis_verify ──
old1 = ("no_reads_stages = {'identification', 'filter', 'cobra', 'merge', 'cluster', "
        "'taxonomy', 'host', 'checkv', 'report', 'rescue', 'analysis'}")
new1 = ("no_reads_stages = {'identification', 'filter', 'cobra', 'merge', 'cluster', "
        "'taxonomy', 'host', 'checkv', 'report', 'rescue', 'analysis', 'analysis_verify'}")
if old1 in src:
    src = src.replace(old1, new1, 1)
    changes.append("缺陷1: no_reads_stages 补入 analysis_verify")
elif new1 in src:
    changes.append("缺陷1: 已修复 (跳过)")
else:
    sys.exit("缺陷1 锚点未匹配")

# ── 缺陷 2: reads_dir 判断目录含序列文件 ──
old2 = ("                reads_dir = self.d['clean'] if self.d['clean'].is_dir() "
        "else (self.d['hostdep'] if self.d['hostdep'].is_dir() else self.d['raw'])")
new2 = ("                def _has_seq(d):\n"
        "                    # 目录存在且含序列文件才可用 (00a_CleanData 可能只剩日志)\n"
        "                    if not d.is_dir():\n"
        "                        return False\n"
        "                    for f in d.iterdir():\n"
        "                        if f.is_file() and any(\n"
        "                                f.name.lower().endswith(e) for e in SEQ_EXTS):\n"
        "                            return True\n"
        "                    return False\n"
        "                reads_dir = (self.d['clean'] if _has_seq(self.d['clean'])\n"
        "                             else (self.d['hostdep'] if _has_seq(self.d['hostdep'])\n"
        "                                   else self.d['raw']))")
if old2 in src:
    src = src.replace(old2, new2, 1)
    changes.append("缺陷2: reads_dir 改为需含序列文件")
elif "def _has_seq(d):" in src:
    changes.append("缺陷2: 已修复 (跳过)")
else:
    sys.exit("缺陷2 锚点未匹配")

if src == orig:
    print("无需修改 (可能已全部修复)")
    sys.exit(0)

if not BAK.exists():
    shutil.copy2(TARGET, BAK)
    print(f"备份: {BAK}")

TARGET.write_text(src, encoding="utf-8")
print("已应用:")
for c in changes:
    print("  -", c)
