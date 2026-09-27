#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""2026-08-30 补丁: (1) 撒回 NON_PLANT_GENERA 里两个真植物病毒属 Miraophiovirus /
Oleurovirus (ICTV VMR MSL41 宿主列 = plants, 非 (S) 型);
(2) 合并双 def normalize_c9 (497 行死副本 + 655 行生效版), 采用带 strip/NA 归一的版本。
带备份 + 断言 + 语法检查 + md5。用法: python3 /tmp/patch_genusretract_e1e2.py
"""
import hashlib
import shutil
import subprocess
import sys
import time

PIPE = "/home/zhangwenda/MMPV-RNA/virome_discovery_pipeline"
CUR = PIPE + "/run_host_prediction.py"
BAK = CUR + ".bak_genusretract_20260830"

old_comment = """#    实测影响: goji 树 0 行受影响 (撤回后 Plant 仍为 1028), 纯粹降低误否决风险
NON_PLANT_GENERA_UNDER_REVIEW = ['Betapartitivirus', 'Geminivirus']"""

new_comment = """#    实测影响: goji 树 0 行受影响 (撤回后 Plant 仍为 1028), 纯粹降低误否决风险
# ── 2026-08-30 二次复核撒回: 属级证据来自 ICTV VMR MSL41 "Host source" 列 (非植物库) ──
#    Miraophiovirus  VMR 13 条记录全为 plants (Aspiviridae); ICTV_MSL41_Genus.tsv 在册,
#                    管线参考库 17 条命中 → 属级证据明确为植物病毒属
#    Oleurovirus     VMR 1 条记录 host=plants (Geminiviridae); ICTV 属表在册
#    撒回判据: VMR 宿主列写 plants 且不带 (S) 标记 (=ICTV 确认宿主, 而非样本/环境来源)
#    实测影响 (内存撤回对拍, 两棵树): goji Plant 1042→1042 (+0);
#                                     onekp Plant 19967→19990 (+23, 明细见审计文档 §10.9)
#    取证脚本: scripts/audit/vmr_genus_host_audit.py / measure_genus_retract_e1.py
NON_PLANT_GENERA_UNDER_REVIEW = ['Betapartitivirus', 'Geminivirus',
                                 'Miraophiovirus', 'Oleurovirus']"""

old_def_dead = """def normalize_c9(h):
    h = str(h).strip() if not pd.isna(h) else 'nan'
    if h in ('', 'Unknown', 'None', 'NA', 'nan'): return 'Unknown'
    if h in ('Insecta', 'Arachnida', 'Aves', 'Human', 'Animal_other'): return 'Animal'
    if h == 'Oomycetes': return 'Protist'
    return h"""

new_def_dead = """# normalize_c9 单一定义在下方 "决策树与整合逻辑" 段 (2026-08-30 合并双定义,
# 原此处副本与下方定义重复; 实测两版对真实 C9 输入判定完全一致, 见审计文档 §10.9)"""

old_def_live = """def normalize_c9(h):
    if pd.isna(h) or h in ['Unknown', 'None']: return 'Unknown'
    h = str(h)
    if h in ['Insecta', 'Arachnida', 'Aves', 'Human', 'Animal_other']: return 'Animal'
    if h in ['Oomycetes']: return 'Protist'
    return h"""

new_def_live = """def normalize_c9(h):
    \"\"\"C9 的 Predicted_Host → 统一宿主类别 (与 read_ictv_hits 同语义)。

    2026-08-30 合并双定义: 文件里原有两处 def normalize_c9, 运行时只有本处生效。
    两版对真实 C9 输入 (14 类封闭取值, 两棵树 0 行含 NA/nan/空白) 判定 0 差异;
    此处保留带 strip + NA 归一的版本, 避免未来 C9 出现 'NA'/空白时把无效值
    当成有效宿主类别计入命中。
    \"\"\"
    h = str(h).strip() if not pd.isna(h) else 'nan'
    if h in ('', 'Unknown', 'None', 'NA', 'nan'): return 'Unknown'
    if h in ('Insecta', 'Arachnida', 'Aves', 'Human', 'Animal_other'): return 'Animal'
    if h == 'Oomycetes': return 'Protist'
    return h"""

src = open(CUR, encoding="utf-8").read()
before_md5 = hashlib.md5(src.encode("utf-8")).hexdigest()
print("改动前 md5:", before_md5)

for tag, old in (("A_属撒回注释+名单", old_comment),
                 ("B_死副本 normalize_c9", old_def_dead),
                 ("C_生效 normalize_c9", old_def_live)):
    n = src.count(old)
    print("  锚点 %-24s 命中 %d 次" % (tag, n))
    if n != 1:
        sys.exit("锚点不唯一/缺失: " + tag)

shutil.copy2(CUR, BAK)
print("备份:", BAK)

src = src.replace(old_comment, new_comment, 1)
src = src.replace(old_def_dead, new_def_dead, 1)
src = src.replace(old_def_live, new_def_live, 1)
open(CUR, "w", encoding="utf-8", newline="").write(src)
after_md5 = hashlib.md5(open(CUR, "rb").read()).hexdigest()
print("改动后 md5:", after_md5)

r = subprocess.run([sys.executable, "-m", "py_compile", CUR], capture_output=True, text=True)
print("py_compile rc=%d %s" % (r.returncode, r.stderr.strip()[:300]))
for line in open(CUR, encoding="utf-8").read().splitlines():
    if "genusretract" in line or "NON_PLANT_GENERA_UNDER_REVIEW" in line or "Miraophiovirus" in line:
        print("   >", line[:120])
print("时间:", time.strftime("%Y-%m-%d %H:%M:%S"))
