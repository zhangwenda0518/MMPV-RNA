#!/usr/bin/env python3
# 重建 <sample>_combined_taxonomy.tsv (merge_taxonomy_results + fill_taxonomy_na)
# 用法:
#   python3 /tmp/rerun_merge_fill.py --module <virus_classifier.py> --classed <dir> --outdir <dir> \
#       [--baseline <现有combined>] [--fill-extra N]
# 说明:
#   --module    指定 virus_classifier.py 路径 (指向备份版可做"重建保真"验证)
#   --baseline  用现有 combined 推断 tools_ran 顺序, 保证行序一致
#   --fill-extra merge 内部已回填 1 次; 额外再回填 N 次 (N=1 可复现线上旧行为=回填两遍)
import os, glob, argparse, importlib.util


def load_module(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def tool_order_from_baseline(path):
    order, cur = [], None
    with open(path) as f:
        f.readline()
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < 2:
                continue
            if p[1] != cur:
                cur = p[1]
                if cur not in order:
                    order.append(cur)
    return order


ap = argparse.ArgumentParser()
ap.add_argument("--module", required=True)
ap.add_argument("--classed", required=True)
ap.add_argument("--outdir", required=True)
ap.add_argument("--sample", default="Votus")
ap.add_argument("--baseline", default="")
ap.add_argument("--tools", default="")
ap.add_argument("--fill-extra", type=int, default=0)
args = ap.parse_args()

vc = load_module(args.module, "vc_mod")

if args.baseline and os.path.exists(args.baseline):
    tools = tool_order_from_baseline(args.baseline)
    src = "baseline"
elif args.tools:
    tools = [t for t in args.tools.split(",") if t]
    src = "arg"
else:
    tools = []
    for f in sorted(glob.glob(os.path.join(args.classed, args.sample + "_*_taxonomy.tsv"))):
        mid = os.path.basename(f)[len(args.sample) + 1: -len("_taxonomy.tsv")]
        tools.append(mid)
    src = "glob"

os.makedirs(args.outdir, exist_ok=True)
# 用 symlink 把 per-tool 文件摆到 outdir, merge 只认 outdir 下的固定文件名
for t in tools:
    s = os.path.join(args.classed, "%s_%s_taxonomy.tsv" % (args.sample, t))
    d = os.path.join(args.outdir, "%s_%s_taxonomy.tsv" % (args.sample, t))
    if not os.path.exists(d):
        os.symlink(os.path.abspath(s), d)

print("[drv] tools(%s) = %s" % (src, ",".join(tools)))
out = vc.merge_taxonomy_results(args.sample, args.outdir, tools)
for _ in range(args.fill_extra):
    print("[drv] fill extra pass")
    vc.fill_taxonomy_na(out, out)
print("[drv] wrote %s" % out)
