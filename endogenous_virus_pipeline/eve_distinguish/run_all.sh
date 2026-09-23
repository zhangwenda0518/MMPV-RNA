#!/bin/bash
# eve_distinguish — DNA 病毒候选 vs 内源病毒元件 (EVE) 判别, 一键后运行脚本 (v5)
# ============================================================================
#
# 两个方法 (详见 README.md):
#   方案一 结构通道 (s1_decay_scan.py)   同源区框内提前终止, 以同一条链另外两个读框的
#                                        终止密度当**同批组装错误率基线**; 再用分布剖面
#                                        区分"单个 indel 断裂"与"分布式化石退化".
#   方案二 功能通道 (s2/s2b_locus_scan)  只问"单一位点上 MP/CP/AP/RT/RH 是否齐全".
#                                        方向倒置: "不全" 不判 EVE (RNA-seq 里真 DNA
#                                        病毒本就拼不全); 有判别力的是"位点齐全 +
#                                        密码子退化" = 整合前病毒.
#   S2b                           候选 blastn 回自身寄主的 OneKP 转录组组装, 按物种级
#                                        区间并集识别"候选就是寄主基因"的假阳性
#                                        (实测两条 HSP70 候选就是这样混过 CDD 的).
#
# 用法 (独立运行; 这里**不是** virome_discovery_pipeline 的一个 stage, 跑完发现管道
# 再单独调用, 以免把判别逻辑焊死在主流程里):
#
#   # 最常用: 指向一次发现管道的输出目录, 候选/证据表按目录约定自动定位
#   bash run_all.sh -D /path/to/onekp-virus -A /path/to/1kp/assemblies -o out/ -t 32
#
#   # 显式给候选 fasta 与证据表
#   bash run_all.sh -q candidates.fasta -e 10_Reports/rescue_evidence_scored.tsv \
#                   -A /path/to/1kp/assemblies -o out/ -t 32
#
#   # 只验证 S1/S2 通道 (不跑寄主比对): host_contamination_likely 无法检测, 仅供调试
#   bash run_all.sh -q candidates.fasta -e evidence.tsv -o out/ --no-host
#
# 参数:
#   -q, --query FILE       待判别候选 fasta (必填, 除非用 -D 自动定位)
#   -e, --evidence FILE    rescue_evidence_scored.tsv, s3 靠它取 tax_family/CheckV (必填)
#   -D, --discovery DIR    发现管道一次输出目录; 按阶段目录约定自动找 query/evidence/out
#   -o, --outdir DIR       产物目录 (默认: -D 时为 <DISCOVERY>/08b_EVE_Distinguish, 否则 ./eve_distinguish_out)
#   -A, --asm-root DIR     OneKP 转录组汇编根目录 (映射表第三列相对它); 寄主否决必需
#   -m, --host-map FILE    三列 TSV: 样本号 <TAB> 4位代码 <TAB> 物种目录名 (默认: 模块自带)
#       --baits FILE       Caulifinder banks (TRUE/FALSE 编码头; 默认: 模块自带 baits.fa)
#       --locus-query FILE 只对 fasta 的子集做寄主 blastn (省时; 默认同 query)
#   -t, --threads N        DIAMOND 线程数 (默认 16)
#       --subj-cap-gb N    解压后寄主组装超过这个体积就跳过, 防 OOM (默认 40)
#       --no-host          不跑寄主比对 (S1/S2 仍可用; 但 host_contamination_likely 检不出)
#   -h, --help             打印本段用法
#
# 产物 (全部落在 -o 指定目录):
#   panel_hits.tsv / baits_hits.tsv / all_hits.tsv   DIAMOND blastx 原始命中
#   s1_decay.tsv / s2_domains.tsv                    两通道逐候选特征
#   locus_blastn.tsv / locus_errors.log              候选 vs 寄主转录组 blastn
#   locus_architecture.tsv                           位点归属 + 寄主同源强度
#   eve_distinguish_verdict.tsv                      逐候选 verdict (s3)
#   dna_vs_eve_filter.tsv                            每条候选一个 action (s4) <- 下游按列过滤
#
set -euo pipefail

MOD="$(cd "$(dirname "$0")" && pwd)"      # 模块目录: 脚本 / panel.fasta / baits.fa / host_sample_map.tsv
PYTHON=${PYTHON:-python3}
# Windows 上多半只装了 python (没有 python3 这个入口), 这里自动降级, 免得 six 个阶段
# 跑到第三个才报 "python3: command not found"
if ! command -v "$PYTHON" >/dev/null 2>&1; then
    if command -v python >/dev/null 2>&1; then
        echo "[run_all] 无 python3, 改用 python ($(command -v python))"
        PYTHON=python
    else
        echo "[run_all] 找不到 python3/python: 设 PYTHON=/path/to/python 再跑" >&2
        exit 1
    fi
fi

usage() {
    sed -n '2,/^set -euo pipefail$/p' "$0" | sed '$d'
}

# Windows/Git Bash 在没有符号链接权限时 ln 会失败, 而 set -e 会让整个运行当场终止,
# 所以软链不成就退回复制: 面板与 dmnd 库都是只读复用, 复制同样可复现。
link_or_copy() {
    if ln -sf "$1" "$2" 2>/dev/null; then
        echo "  软链 $2 -> $1"
    else
        cp -f "$1" "$2"
        echo "  无法软链 (无符号链接权限?), 已复制 $2"
    fi
}

# ─────────────────────────── 参数 ───────────────────────────
OUT=""
QUERY=""
EVIDENCE=""
DISCOVERY=""
ASM_ROOT=""
HOST_MAP="$MOD/host_sample_map.tsv"
BAITS="$MOD/baits.fa"
DIAMOND=${DIAMOND:-diamond}
THREADS=16
LOCUS_QUERY=""
SUBJ_CAP_GB=40
NO_HOST=0

while [ $# -gt 0 ]; do
    case "$1" in
        -q|--query)        QUERY="$2"; shift 2 ;;
        -e|--evidence)     EVIDENCE="$2"; shift 2 ;;
        -D|--discovery)    DISCOVERY="$2"; shift 2 ;;
        -o|--outdir)       OUT="$2"; shift 2 ;;
        -A|--asm-root)     ASM_ROOT="$2"; shift 2 ;;
        -m|--host-map)     HOST_MAP="$2"; shift 2 ;;
        --baits)           BAITS="$2"; shift 2 ;;
        --locus-query)     LOCUS_QUERY="$2"; shift 2 ;;
        -t|--threads)      THREADS="$2"; shift 2 ;;
        --subj-cap-gb)     SUBJ_CAP_GB="$2"; shift 2 ;;
        --no-host)         NO_HOST=1; shift ;;
        -h|--help)         usage; exit 0 ;;
        *) echo "[run_all] 未知参数: $1 (--help 看用法)" >&2; usage >&2; exit 1 ;;
    esac
done

# ─────────────────────── 输入定位与校验 ───────────────────────
# -D 给的是发现管道的一次输出目录, 按它的阶段目录约定找候选与证据表:
#   10_Reports/rescue_evidence_scored.tsv       科属/CheckV/CDD 证据 (s3 读它)
#   10_Reports/eve_candidates.fasta             新候选集合 (rescue 产出)
#   eve_screen/evescan_query.fasta              早期筛选集合
if [ -n "$DISCOVERY" ]; then
    if [ -z "$QUERY" ]; then
        for c in "$DISCOVERY/10_Reports/eve_candidates.fasta" \
                 "$DISCOVERY/eve_screen/evescan_query.fasta"; do
            if [ -s "$c" ]; then QUERY="$c"; break; fi
        done
    fi
    if [ -z "$EVIDENCE" ]; then
        for e in "$DISCOVERY/10_Reports/rescue_evidence_scored.tsv"; do
            if [ -s "$e" ]; then EVIDENCE="$e"; break; fi
        done
    fi
    [ -n "$OUT" ] || OUT="$DISCOVERY/08b_EVE_Distinguish"
fi
[ -n "$OUT" ] || OUT="$PWD/eve_distinguish_out"

if [ -z "$QUERY" ]; then
    echo "[run_all] 缺少候选 fasta. 用 -q/--query 指定, 或用 -D/--discovery 指向发现管道" >&2
    echo "          的一次输出目录 (候选默认可落在 <DISCOVERY>/10_Reports/eve_candidates.fasta" >&2
    echo "          或 <DISCOVERY>/eve_screen/evescan_query.fasta)" >&2
    exit 1
fi
if [ ! -s "$QUERY" ]; then
    echo "[run_all] 候选 fasta 不存在或为空: $QUERY" >&2; exit 1
fi
if [ -z "$EVIDENCE" ] || [ ! -s "$EVIDENCE" ]; then
    echo "[run_all] 缺少证据表 rescue_evidence_scored.tsv (s3 靠它取 tax_family/CheckV)." >&2
    echo "          用 -e/--evidence 指定, 或用 -D/--discovery 指向发现管道输出目录." >&2
    exit 1
fi
[ -n "$LOCUS_QUERY" ] || LOCUS_QUERY="$QUERY"

# 寄主通道是唯一能识别"候选其实就是寄主基因"的手段, 缺了会静默漏掉整类假阳性
# (早期版本就因为 ASM_ROOT 配错, 23 个物种全记 #NOGZ, 宿主污染 55 条只报出 16 条).
# 所以这里默认直接报错, 必须显式 --no-host 才允许不带.
if [ "$NO_HOST" -eq 0 ] && { [ -z "$ASM_ROOT" ] || [ ! -d "$ASM_ROOT" ]; }; then
    echo "[run_all] 缺少 OneKP 转录组汇编根目录 (-A/--asm-root): 寄主同源否决无法运行," >&2
    echo "          host_contamination_likely 将一条都检不出来. 确认不需要就加 --no-host." >&2
    exit 1
fi
# 只校验这一步真用得上的命令: --no-host 时 blastn 根本不会被调用, 不该因为它缺失而失败
command -v "$DIAMOND" >/dev/null 2>&1 || { echo "[run_all] 找不到命令: $DIAMOND" >&2; exit 1; }
if [ "$NO_HOST" -eq 0 ]; then
    command -v blastn >/dev/null 2>&1 || { echo "[run_all] 找不到命令: blastn (寄主比对需要)" >&2; exit 1; }
fi

mkdir -p "$OUT"
cd "$OUT"
echo "[run_all] MOD=$MOD"
echo "[run_all] OUT=$OUT"
echo "[run_all] QUERY=$QUERY ($(grep -c '^>' "$QUERY") 条)"
echo "[run_all] EVIDENCE=$EVIDENCE"
if [ "$NO_HOST" -eq 1 ]; then
    echo "[run_all] 寄主比对已关闭 (--no-host): 不检测 host_contamination_likely"
else
    echo "[run_all] ASM_ROOT=$ASM_ROOT  HOST_MAP=$HOST_MAP"
fi

# ─────────────────────── [0/6] 参考面板 ───────────────────────
# dmnd 复用: 模块里已建过库(且比源库新)就直接软链, 省一次 makedb
if [ ! -f panel.dmnd ] && [ -f "$MOD/panel.dmnd" ] && [ "$MOD/panel.dmnd" -nt "$MOD/panel.fasta" ]; then
    link_or_copy "$MOD/panel.dmnd" panel.dmnd
fi
if [ ! -f baits.dmnd ] && [ -f "$MOD/baits.dmnd" ] && [ -f "$BAITS" ] && [ "$MOD/baits.dmnd" -nt "$BAITS" ]; then
    link_or_copy "$MOD/baits.dmnd" baits.dmnd
fi

echo "=== [0/6] 参考面板 ==="
if [ ! -f panel.fasta ]; then
    if [ -f "$MOD/panel.fasta" ]; then
        # 面板是随模块提交的固定产物: 直接用, 不重新联网抓取 (保证可复现)
        link_or_copy "$MOD/panel.fasta" panel.fasta
        echo "复用模块面板 $MOD/panel.fasta"
    else
        echo "模块无 panel.fasta, 走联网构建 (需要 NCBI 可达)"
        $PYTHON "$MOD/build_panel.py" .
    fi
else
    echo "panel.fasta 已存在, 跳过"
fi

# ─────────────────────── [1/6] DIAMOND 建库 ───────────────────────
# stderr 被丢掉的命令必须自己接住 rc: 否则 diamond 被 OOM killer 干掉或段错误时,
# set -e 只让它静默退出, 上游看不到任何原因
run_or_die() {
    local what="$1"; shift
    if "$@"; then return 0; fi
    local rc=$?
    echo "[run_all] FAILED($rc): $what" >&2
    echo "[run_all]   命令: $*" >&2
    echo "[run_all]   $DIAMOND 是否可执行? 内存/磁盘是否够? 先把这条单独跑一遍看真实报错" >&2
    exit "$rc"
}

echo "=== [1/6] DIAMOND 建库 ==="
[ -f panel.dmnd ] || $DIAMOND makedb --in panel.fasta -d panel >/dev/null \
    || { echo "[run_all] FAILED: DIAMOND makedb (panel) rc=$?" >&2; exit 1; }
if [ -f "$BAITS" ]; then
    [ -f baits.dmnd ] || $DIAMOND makedb --in "$BAITS" -d baits >/dev/null \
        || { echo "[run_all] FAILED: DIAMOND makedb (baits) rc=$?" >&2; exit 1; }
else
    echo "!! 未找到 $BAITS (Caulifinder banks), 功能通道只跑 panel -> TE 否决不可用"
    : > baits_hits.tsv
fi

# ─────────────────────── [2/6] blastx ───────────────────────
OUTFMT="6 qseqid sseqid pident length qstart qend sstart send evalue bitscore qframe stitle"
echo "=== [2/6] blastx vs panel / baits ==="
run_or_die "DIAMOND blastx (panel)" $DIAMOND blastx -q "$QUERY" -d panel.dmnd -o panel_hits.tsv \
    --outfmt $OUTFMT -e 1e-5 --max-target-seqs 20 --threads "$THREADS" --sensitive
if [ -f baits.dmnd ]; then
    run_or_die "DIAMOND blastx (baits)" $DIAMOND blastx -q "$QUERY" -d baits.dmnd -o baits_hits.tsv \
        --outfmt $OUTFMT -e 1e-5 --max-target-seqs 20 --threads "$THREADS" --sensitive
else
    echo "!! 无 baits.dmnd, 功能通道只跑 panel -> TE 否决不可用"
    : > baits_hits.tsv
fi
# 0 命中也该有个空文件; 没有就说明上面哪一步没落盘, 补空表让下游能跑完并留下痕迹
for f in panel_hits.tsv baits_hits.tsv; do
    [ -f "$f" ] || { echo "!! $f 未生成, 按空表继续 (疑似磁盘满/权限问题)" >&2; : > "$f"; }
done

# ───────────────── [3/6] S1 结构通道 + S2 功能通道 ─────────────────
echo "=== [3/6] S1 结构通道 + S2 功能通道 ==="
cat panel_hits.tsv baits_hits.tsv > all_hits.tsv
$PYTHON "$MOD/s1_decay_scan.py" "$QUERY" all_hits.tsv s1_decay.tsv
$PYTHON "$MOD/s2_domain_scan.py" "$QUERY" panel_hits.tsv baits_hits.tsv s2_domains.tsv

# ───────── [4/6] 候选 vs 寄主转录组 blastn (位点归属) + S2b ─────────
# HOST_MAP: 三列 TSV  "样本号 <TAB> 4位代码 <TAB> ASM_ROOT 下的物种目录名"
#           只有具备 OneKP 汇编的行会被 blast; 无汇编的样本自动跳过并记 locus_errors.log
FMT='6 qseqid sseqid pident length mismatch gaps qstart qend sstart send evalue bitscore qlen slen'
echo "=== [4/6] 候选 vs 寄主转录组 blastn + S2b ==="
if [ "$NO_HOST" -eq 0 ] && [ -f "$HOST_MAP" ] && [ -s "$HOST_MAP" ]; then
    : > locus_blastn.tsv
    : > locus_errors.log
    SUBJ_CAP_GB=${SUBJ_CAP_GB:-40}                # 解压后超过这个体积就跳过, 防 OOM
    while IFS=$'\t' read -r sample code dir || [ -n "$sample" ]; do
        [ -z "$sample" ] && continue
        asm=""
        for f in "$ASM_ROOT/$dir/${code}-SOAPdenovo-Trans-assembly.fa.gz" \
                 "$ASM_ROOT/$dir/${code}-SOAPdennov-Trans-assembly.fa.gz"; do
            if [ -s "$f" ]; then asm="$f"; break; fi
        done
        if [ -z "$asm" ]; then
            echo "#NOGZ $sample $dir/${code}-SOAPdenovo-Trans-assembly.fa.gz" >> locus_errors.log
            continue
        fi
        # 必须落到临时文件再比对: `-subject <(gunzip -c ...)` 是管道, blastn 无法 seek,
        # 会把整个 assembly 缓存进内存 —— 实测被 OOM killer 干掉十几次 (最大一次 rss 266GB).
        # 用真实文件时 blastn 流式读 subject, 内存占用降一个量级.
        subj=$(mktemp "${TMPDIR:-/tmp}/evehost_${code}_XXXXXX.fa")
        if ! gunzip -c "$asm" > "$subj" 2>>locus_errors.log; then
            echo "#GUNZIP_FAIL $sample $asm" >> locus_errors.log
            rm -f "$subj"; continue
        fi
        sz_gb=$(( $(stat -c%s "$subj" 2>/dev/null || echo 0) / 1073741824 ))
        if [ "$sz_gb" -gt "$SUBJ_CAP_GB" ]; then
            echo "#TOOBIG $sample $asm ${sz_gb}GB > ${SUBJ_CAP_GB}GB" >> locus_errors.log
            rm -f "$subj"; continue
        fi
        blastn -query "$LOCUS_QUERY" -subject "$subj" \
               -task megablast -evalue 1e-5 -outfmt "$FMT" -num_threads 1 \
               -max_target_seqs 200 2>>locus_errors.log >> locus_blastn.tsv
        rm -f "$subj"
    done < "$HOST_MAP"
    n_host=$(cut -f2 locus_blastn.tsv | sort -u | grep -c . || true)
    n_err=$(grep -c . locus_errors.log || true)
    echo "寄主 blastn: $(wc -l < locus_blastn.tsv) 行命中 / $n_host 个寄主 scaffold; 跳过或失败 $n_err 行 (见 locus_errors.log)"
    if [ "$n_host" -eq 0 ]; then
        echo "[run_all] !! 没有任何 scaffold 比对成功, 寄主同源否决整条失效 ——" >&2
        echo "           检查 -A 指向的目录里是否真有 <物种目录>/<代码>-SOAPdenovo-Trans-assembly.fa.gz," >&2
        echo "           以及 HOST_MAP 第三列与目录名是否一致 (locus_errors.log 里的 #NOGZ 行即没找到的)." >&2
    fi
    $PYTHON "$MOD/s2b_locus_scan.py" locus_blastn.tsv s2_domains.tsv s1_decay.tsv \
        locus_architecture.tsv "$HOST_MAP"
elif [ "$NO_HOST" -eq 1 ]; then
    echo "--no-host: 跳过寄主比对, 位点级架构全部记为 unresolved"
    {
    printf 'contig_id\tsample\thost_scaffold\thost_code\thost_wpid\thost_cov\thost_xeno\t'
    printf 'comps_contig\tbest_locus_comps\tdecay_class\tte_flag\tlocus_arch\tlocus_ncomp\t'
    printf 'locus_comps\tlocus_members\tlocus_hint\n'
    } > locus_architecture.tsv
else
    echo "!! 未找到 $HOST_MAP, 跳过寄主比对: 位点级架构全部记为 unresolved" >&2
    echo "   (可用 -m/--host-map 指定三列文件: 样本号/4位代码/物种目录名)" >&2
    {
    printf 'contig_id\tsample\thost_scaffold\thost_code\thost_wpid\thost_cov\thost_xeno\t'
    printf 'comps_contig\tbest_locus_comps\tdecay_class\tte_flag\tlocus_arch\tlocus_ncomp\t'
    printf 'locus_comps\tlocus_members\tlocus_hint\n'
    } > locus_architecture.tsv
fi

# ─────────────────────── [5/6] S3 合成判别 ───────────────────────
echo "=== [5/6] S3 合成判别 ==="
$PYTHON "$MOD/s3_verdict.py" s1_decay.tsv s2_domains.tsv locus_architecture.tsv "$EVIDENCE" \
    eve_distinguish_verdict.tsv

# ─────────────────── [6/6] verdict -> action ───────────────────
echo "=== [6/6] verdict -> 可执行动作 ==="
$PYTHON "$MOD/s4_filter.py" eve_distinguish_verdict.tsv -o dna_vs_eve_filter.tsv || {
    rc=$?; echo "[run_all] s4_filter 返回 $rc: verdict 表里有未映射的新值, 已按 REVIEW 兜底" >&2; }

echo
echo "全部完成 -> $OUT"
echo "  逐候选判别  $OUT/eve_distinguish_verdict.tsv"
echo "  下游过滤表  $OUT/dna_vs_eve_filter.tsv"
