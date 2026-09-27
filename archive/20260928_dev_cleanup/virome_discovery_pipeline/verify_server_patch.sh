#!/bin/bash
echo "--- stray files in ~ ---"
find ~ -maxdepth 1 -name "*D:*" 2>/dev/null
echo "--- stray files in pipeline dir ---"
find ~/MMPV-RNA/virome_discovery_pipeline -maxdepth 1 -name "*D:*" 2>/dev/null
echo "done"
