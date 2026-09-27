import sys
sys.path.insert(0, "/tmp")
# 读 v2 补丁里的黑名单
import importlib.util
spec = importlib.util.spec_from_file_location("pv2", "/tmp/patch_host_blacklist_v2.py")
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)

auth = set(l.strip() for l in open("/tmp/auth_plant_genera.txt") if l.strip())

print(f"黑名单科: {len(m.NON_PLANT_FAMILIES)}, 黑名单属: {len(m.NON_PLANT_GENERA)}")
print(f"权威库属: {len(auth)}")

# 检查黑名单属是否误入权威库 (冲突项)
conflict = sorted(g for g in m.NON_PLANT_GENERA if g in auth)
print(f"\n=== 黑名单属 ∩ 权威库 (冲突/需人工复核) : {len(conflict)} 个 ===")
for g in conflict:
    print(f"  {g}")

print(f"\n黑名单属中不在权威库的: {len(m.NON_PLANT_GENERA) - len(conflict)} 个")
