#!/bin/bash
set -u
cd /home/zhangwenda/MMPV-RNA/virome_discovery_pipeline || exit 1
echo "########## 1. 我引入的 4 个名字在文件里的出现位置 ##########"
for n in PLANT_FAMILIES_WHITELIST PLANT_GENERA_WHITELIST ICTV_PLANT_INCLUDING_FAMILIES WHITELIST_OVERRIDES_FAMILY_VETO; do
  echo "--- $n"
  grep -n "$n" run_host_prediction.py
done
echo
echo "########## 2. 谁 import 了 run_host_prediction ##########"
grep -rn --include=*.py -e "import run_host_prediction" -e "from run_host_prediction" . | grep -v "^./run_host_prediction.py" || echo "(无其他文件 import)"
echo
echo "########## 3. run_host_prediction.py 的模块级可执行语句 (import 时是否有副作用) ##########"
grep -n -E "^(print|os\.|subprocess|sys\.exit|raise|logging)" run_host_prediction.py || echo "(无模块级副作用语句)"
echo
echo "########## 4. 是否有其他脚本自己解析 Plant.tsv / 维护同类名单 ##########"
grep -rln "Plant.tsv" --include=*.py . | head -20
echo
echo "########## 5. 备份链 md5 ##########"
md5sum run_host_prediction.py run_host_prediction.py.bak_plantwl_20260915 2>/dev/null
