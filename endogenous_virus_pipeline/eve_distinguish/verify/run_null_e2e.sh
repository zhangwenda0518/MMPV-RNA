#!/bin/bash
# run_null_e2e.sh — 端到端空模型: 把洗牌候选当真候选, 跑完整条判别流水线
# ============================================================================
#
# 这是"诱饵库"式阴性对照的判别器版本。把每条候选的核苷酸单核苷酸洗牌（严格保组成、
# 破密码子结构与读框连续性），重跑 diamond + s1 + s2 + s3 + s4，再和真实候选的产物
# 逐 stage 比。判据的富集倍数 = P(判据|真实) / P(判据|洗牌): 接近 1 说明该判据只是
# 组成噪声的代理。
#
# 寄主通道两条腿都关掉（--no-host 口径）: 洗牌序列去 blastn 寄主转录组没有意义，
# 真实侧也关掉才可比。s3 的寄主否决分支因此不参与，本对照只测 s1×s2 的判别力。
#
# 用法:
#   bash run_null_e2e.sh --mod DIR --run-dir REAL --work WORK \
#        --evidence evidence.tsv [--threads 32] [--seed 20260924]
#
# 产物 (全在 WORK/): shuf.fa, *_shuf.tsv, real_nohost/*.tsv
set -euo pipefail

MOD=""; REAL=""; WORK=""; EVID=""; THREADS=32; SEED=20260924
while [ $# -gt 0 ]; do
    case "$1" in
        --mod) MOD="$2"; shift 2;;
        --run-dir) REAL="$2"; shift 2;;
        --work) WORK="$2"; shift 2;;
        --evidence) EVID="$2"; shift 2;;
        --threads) THREADS="$2"; shift 2;;
        --seed) SEED="$2"; shift 2;;
        *) echo "未知参数: $1" >&2; exit 2;;
    esac
done
for v in MOD REAL WORK EVID; do
    [ -n "${!v}" ] || { echo "缺 --$(echo $v | tr 'A-Z_' 'a-z-')" >&2; exit 2; }
done

PY=${PYTHON:-python3}
DIAMOND=${DIAMOND:-diamond}
mkdir -p "$WORK"
cd "$WORK"
echo "[null] MOD=$MOD  REAL=$REAL  WORK=$WORK"

# ── 1) 洗牌候选 ────────────────────────────────────────────────
$PY "$MOD/verify/null_model.py" emit --run-dir "$REAL" --mode contig \
    -o shuf.fa --seed "$SEED"

# ── 2) 洗牌候选的 blastx (参数与 run_all.sh 完全一致) ──────────
# 注意 --outfmt 后面必须**不加引号**让 shell 拆词: diamond 的 --outfmt 是"吃掉后面
# 所有非选项参数"的多值形式, 整体作为一个 argv 传进去会被判 Invalid output format.
OUTFMT='6 qseqid sseqid pident length qstart qend sstart send evalue bitscore qframe stitle'
for pair in "panel:panel_hits" "baits:baits_hits"; do
    dmnd="${pair%%:*}"; stem="${pair##*:}"
    if [ ! -f "$MOD/$dmnd.dmnd" ]; then
        echo "[null] 缺 $MOD/$dmnd.dmnd, 先建库或从已跑过的目录软链过来" >&2; exit 1
    fi
    echo "[null] diamond blastx ($dmnd)"
    $DIAMOND blastx -q shuf.fa -d "$MOD/$dmnd.dmnd" -o "shuf_${stem}.tsv" \
        --outfmt $OUTFMT -e 1e-5 --max-target-seqs 20 --threads "$THREADS" --sensitive
done
cat shuf_panel_hits.tsv shuf_baits_hits.tsv > shuf_all.tsv

# ── 3) s1 / s2 ─────────────────────────────────────────────────
echo "[null] s1 + s2"
$PY "$MOD/s1_decay_scan.py" shuf.fa shuf_all.tsv s1_shuf.tsv
$PY "$MOD/s2_domain_scan.py" shuf.fa shuf_panel_hits.tsv shuf_baits_hits.tsv s2_shuf.tsv

# ── 4) s3 / s4: 两侧都关寄主通道, 才可比 ───────────────────────
mkdir -p real_nohost
$PY "$MOD/s2b_locus_scan.py" --emit-header > real_nohost/locus_architecture.tsv
echo "[null] s3 + s4 (洗牌侧)"
$PY "$MOD/s3_verdict.py" s1_shuf.tsv s2_shuf.tsv real_nohost/locus_architecture.tsv \
    "$EVID" verdict_shuf.tsv
$PY "$MOD/s4_filter.py" verdict_shuf.tsv -o filter_shuf.tsv || true
echo "[null] s3 + s4 (真实侧, 同口径)"
$PY "$MOD/s3_verdict.py" "$REAL/s1_decay.tsv" "$REAL/s2_domains.tsv" \
    real_nohost/locus_architecture.tsv "$EVID" real_nohost/verdict.tsv
$PY "$MOD/s4_filter.py" real_nohost/verdict.tsv -o real_nohost/filter.tsv || true

# ── 5) 比较 ────────────────────────────────────────────────────
echo
$PY "$MOD/verify/null_model.py" compare --stage s1 --real "$REAL/s1_decay.tsv" \
    --null s1_shuf.tsv --json cmp_s1.json
echo
$PY "$MOD/verify/null_model.py" compare --stage s2 --real "$REAL/s2_domains.tsv" \
    --null s2_shuf.tsv --json cmp_s2.json
echo
$PY "$MOD/verify/null_model.py" compare --stage s3 --real real_nohost/verdict.tsv \
    --null verdict_shuf.tsv --json cmp_s3.json
echo
echo "[null] 完成, 各 stage 对照 JSON 在 $WORK/"
