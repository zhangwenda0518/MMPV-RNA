#!/bin/bash
cd ~/MMPV-RNA/virome_discovery_pipeline || exit 1
echo "=== backup ==="
cp virome_pipeline.py virome_pipeline.py.bak_skip03_20260827
ls -la virome_pipeline.py.bak_skip03_20260827
echo "=== deploy new (already scp'd to /tmp/coasm_retro/vp_new.py) ==="
cp /tmp/coasm_retro/vp_new.py virome_pipeline.py
echo "=== py_compile ==="
source ~/mambaforge/etc/profile.d/conda.sh
conda activate base
python -m py_compile virome_pipeline.py && echo DEPLOY_OK
echo "=== md5 confirm new version live ==="
md5sum virome_pipeline.py /tmp/coasm_retro/vp_new.py
echo "=== grep stage-03 skip marker present ==="
grep -c "if s not in ('cobra', 'merge')" virome_pipeline.py
