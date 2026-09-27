#!/bin/bash
cd ~/MMPV-RNA/virome_discovery_pipeline || exit 1
source ~/mambaforge/etc/profile.d/conda.sh
conda activate base
python - <<'EOF'
import sys
for m in list(sys.modules):
    if m.startswith("vp_"): del sys.modules[m]
sys.path.insert(0, "/tmp/coasm_retro")
import vp_new
results=[]
def check(n,ok): results.append(ok); print(("PASS" if ok else "FAIL"), n)
src = open("/tmp/coasm_retro/vp_new.py", encoding="utf-8").read()
i = src.index("    if getattr(args, 'coassembly', False):")
j = src.index("    stages_to_run =")
seg = src[i:j]
lines = seg.splitlines()
lines = [l for l in lines if not l.strip().startswith("#")]
# 整体去 4 格缩进
ded = "\n".join(l[4:] if l.startswith("    ") else l for l in lines)
ns = {"stage_order": ['clean','deplete','assembly','identification','filter','cobra','merge','cluster',
                      'taxonomy','host','checkv','rescue','analysis','analysis_verify','report']}
exec("stages={'all'}\n_all='all' in stages\ncoassembly=True\nargs=type('A',(),{'coassembly':True})()\nlogger=type('L',(),{'info':lambda *a,**k:None})()\n" + ded, ns)
check("Stage03(cobra+merge) 剔除", 'cobra' not in ns['stage_order'] and 'merge' not in ns['stage_order'])
check("cluster 及下游保留", all(s in ns['stage_order'] for s in ['cluster','taxonomy','host','checkv','rescue','analysis','analysis_verify','report']))
print("stage_order =", ns['stage_order'])
print("summary:", sum(results), "/", len(results))
sys.exit(0 if all(results) else 1)
EOF
