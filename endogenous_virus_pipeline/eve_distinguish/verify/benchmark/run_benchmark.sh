#!/bin/bash
# run_benchmark.sh — 把策展基准集喂进判别器, 出逐类判定矩阵
# ============================================================================
# 这是灵敏度/召回那一环: 与空模型相反, 空模型证明"随机不触发", 这里量"真的会不会被抓住"。
#
# 与产线同口径的取舍 (必须写明, 否则数字会被误读):
#   * **寄主通道关闭** (s2b --emit-header): 基准序列不在任何 OneKP 样本里, 没有寄主
#     转录组可比, 所以 host_* 三类结论结构上不出现。这也意味着本矩阵**不测**寄主否决。
#   * 证据表沿用产线那份: 基准集的 contig_id 不在其中 -> 上游 tax_family / CheckV /
#     各条同源通道全为空。所以测的是 **S1×S2 的判定逻辑**, 不含上游证据通道。
#   * diamond 参数与 run_all.sh 完全一致 (--outfmt 多值形式不能加引号)。
#
# 用法: bash run_benchmark.sh --bench BENCH --work WORK --mod MOD --evidence E.tsv
set -euo pipefail

BENCH=""; WORK=""; MOD=""; EVID=""; THREADS=32
while [ $# -gt 0 ]; do
    case "$1" in
        --bench) BENCH="$2"; shift 2;;
        --work) WORK="$2"; shift 2;;
        --mod) MOD="$2"; shift 2;;
        --evidence) EVID="$2"; shift 2;;
        --threads) THREADS="$2"; shift 2;;
        *) echo "未知参数: $1" >&2; exit 2;;
    esac
done
for v in BENCH WORK MOD EVID; do : "${!v:?缺 --$(echo $v | tr 'A-Z_' 'a-z-')}"; done

PY=${PYTHON:-python3}
DIAMOND=${DIAMOND:-diamond}
HERE="$(cd "$(dirname "$0")" && pwd)"
mkdir -p "$WORK"; cd "$WORK"

echo "=== 1/4 去重 + 类别映射 ==="
$PY "$HERE/prepare.py" --bench "$BENCH" --out "$WORK"

echo
echo "=== 2/4 diamond blastx (panel + baits, 参数同产线) ==="
OUTFMT='6 qseqid sseqid pident length qstart qend sstart send evalue bitscore qframe stitle'
for pair in "panel:panel_hits" "baits:baits_hits"; do
    d="${pair%%:*}"; stem="${pair##*:}"
    [ -f "$MOD/$d.dmnd" ] || { echo "缺 $MOD/$d.dmnd" >&2; exit 1; }
    echo "  blastx vs $d"
    $DIAMOND blastx -q bench.fa -d "$MOD/$d.dmnd" -o "bench_${stem}.tsv" \
        --outfmt $OUTFMT -e 1e-5 --max-target-seqs 20 --threads "$THREADS" --sensitive
done
cat bench_panel_hits.tsv bench_baits_hits.tsv > bench_all.tsv

echo
echo "=== 3/4 s1 / s2 / s3 / s4 ==="
$PY "$MOD/s1_decay_scan.py" bench.fa bench_all.tsv s1_decay.tsv
$PY "$MOD/s2_domain_scan.py" bench.fa bench_panel_hits.tsv bench_baits_hits.tsv s2_domains.tsv
$PY "$MOD/s2b_locus_scan.py" --emit-header > locus_architecture.tsv
$PY "$MOD/s3_verdict.py" s1_decay.tsv s2_domains.tsv locus_architecture.tsv \
    "$EVID" eve_distinguish_verdict.tsv
$PY "$MOD/s4_filter.py" eve_distinguish_verdict.tsv -o dna_vs_eve_filter.tsv || true

echo
echo "=== 4/4 判定矩阵 ==="
$PY "$HERE/tabulate.py" --verdict eve_distinguish_verdict.tsv \
    --filter dna_vs_eve_filter.tsv --map class_map.tsv \
    --json "$WORK/benchmark_result.json"
