#!/bin/bash
# 只读：6 个争议属名的宿主证据盘点（ICTV 属表 / VMR / 管线参考库 / PlantVirusDB 植物库）
echo "=== [0] 关键参照文件定位 ==="
find ~/plant_virus_db -maxdepth 4 -name "VMR_MSL41.tsv" 2>/dev/null
find ~/plant_virus_db -maxdepth 4 -name "Plant.tsv" 2>/dev/null
find ~/MMPV-RNA ~/MMPV-paper -maxdepth 4 -name "ICTV_MSL41_Genus.tsv" 2>/dev/null
find ~/MMPV-RNA -maxdepth 4 -name "final.cluster.ref_info.tsv" 2>/dev/null
find ~/ -maxdepth 2 -name "*.tsv" -path "*ictv*" 2>/dev/null | head

echo
echo "=== [1] 参照文件表头（列号） ==="
python3 - <<'PY'
import os
cands = [
 "/home/zhangwenda/plant_virus_db/1.virus-host_db/B-ictv/VMR_MSL41.tsv",
 "/home/zhangwenda/database/virus-db/acvirus_db/taxa.txt",
]
for p in cands:
    print("---", p, "EXISTS" if os.path.isfile(p) else "MISSING")
    if os.path.isfile(p):
        with open(p, encoding="utf-8", errors="replace") as f:
            first = f.readline().rstrip("\n")
        sep = "\t" if "\t" in first else ("," if "," in first else None)
        cols = first.split(sep) if sep else [first]
        for i, c in enumerate(cols, 1):
            print("   %2d %s" % (i, c.strip()))
PY
