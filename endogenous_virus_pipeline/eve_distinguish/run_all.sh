#!/bin/bash
# eve_distinguish — DNA 病毒候选 vs 内源病毒元件 (EVE) 判别, 一键后运行脚本 (v6.3)
# v6.3 变更: ① s1 修 stop_profile 的 3nt 高估 (旧版把终止密码子自身 3nt 算进"无终止段");
#            ② 面板升 v3 (bare polyprotein 兜底标 AP+RT+RH + build_panel --supplement
#            属级补充, Badnavirus 13 -> 39 条, 三个零覆盖属补齐);
#            ③ 新增 s2c 参考基因组通道 (可选 -G, 见下方"方案三")。
#            对既有 run 目录量化影响面: python3 verify/rerun_diff.py --run-dir ... --evidence ...
# 版本标记约定: **模块版本**只有一处权威声明 —— 本文件头与 README 标题 (当前 v6.3)。
# 各 stage 文件里写的 `[判定逻辑最后变更: 模块 vN]` 指的是**该文件自己的逻辑**最后一次改动
# 发生在模块的哪个版本, 不是另一套并列的版本号。所以 s2b 停在 v5.1 是正常的:
# 它自那以后没再改过判定逻辑, 而模块整体继续往前走 (s1/s2c 的逻辑变更于 v6.3)。
# ============================================================================
#
# 两个方法 (详见 README.md):
#   方案一 结构通道 (s1_decay_scan.py)   同源区框内提前终止, 以同一条链另外两个读框的
#                                        终止密度当**本地读框基线**(扣掉"这段序列本来
#                                        就容易读出终止"的组成本底; 注意它扣不掉组装
#                                        indel); 再用分布剖面区分"单个 indel 断裂"与
#                                        "分布式化石退化".
#   方案二 功能通道 (s2/s2b_locus_scan)  只问"单一位点上 MP/CP/AP/RT/RH 是否**各有独立
#                                        证据**". 方向倒置: "不全"不判 EVE (RNA-seq 里
#                                        真 DNA 病毒本就拼不全); 有判别力的是"位点齐全 +
#                                        密码子退化" = 整合前病毒.
#   S2b                           候选 blastn 回自身寄主的 OneKP 转录组组装, 按物种级
#                                        区间并集识别"候选就是寄主基因"的假阳性
#                                        (实测两条 HSP70 候选就是这样混过 CDD 的).
#   方案三 参考基因组通道 (s2c, v6.3, 可选 -G)
#                                 候选 blastn 回自身物种的**染色体级参考基因组**: 两端
#                                        宿主侧翼同一位点 = 前病毒铁证; 宿主外显子嵌合
#                                        = 转座录 EVE 的代理证据. 没有现成基因组的物种
#                                        自动跳过 (no_ref_genome), 跳过绝不是"真病毒"
#                                        的证据 —— 这条通道的证据是单向的.
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
#   # 加开参考基因组通道: 有现成基因组的物种做整合信号判定 (清单格式见模板)
#   bash run_all.sh -q candidates.fasta -e evidence.tsv -A /asm -t 32 \
#                   -G genome_manifest.tsv
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
#   -G, --genome-manifest FILE  参考基因组清单 (三列: 4位代码 <TAB> 注释 <TAB> fasta 路径;
#                          模板 genome_manifest.template.tsv)。给开才跑 s2c 参考基因组
#                          通道; 清单里没有的物种自动跳过。需要 -m 映射表做样本归属.
#   -t, --threads N        DIAMOND 线程数 (默认 16)
#       --subj-cap-gb N    解压后寄主组装/参考基因组超过这个体积就跳过, 防 OOM (默认 40)
#       --no-host          不跑寄主比对 (S1/S2 仍可用; 但 host_contamination_likely 检不出)
#   -h, --help             打印本段用法
#
# 产物 (全部落在 -o 指定目录):
#   panel_hits.tsv / baits_hits.tsv / all_hits.tsv   DIAMOND blastx 原始命中
#   s1_decay.tsv / s2_domains.tsv                    两通道逐候选特征
#   locus_blastn.tsv / locus_errors.log              候选 vs 寄主转录组 blastn
#   locus_architecture.tsv                           位点归属 (未映射候选各自成组, 不伪合并)
#                                                   + 寄主同源强度与口径 (host_scope)
#   refgenome_blastn.tsv / refgenome_errors.log      候选 vs 自身参考基因组 blastn (仅 -G)
#   refgenome_evidence.tsv                           逐候选整合证据 rg_call (仅 -G)
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
GENOME_MANIFEST=""

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
        -G|--genome-manifest) GENOME_MANIFEST="$2"; shift 2 ;;
        -t|--threads)      THREADS="$2"; shift 2 ;;
        --subj-cap-gb)     SUBJ_CAP_GB="$2"; shift 2 ;;
        --no-host)         NO_HOST=1; shift ;;
        -h|--help)         usage; exit 0 ;;
        *) echo "[run_all] 未知参数: $1 (--help 看用法)" >&2; usage >&2; exit 1 ;;
    esac
done

# ─────────────────────── 输入定位与校验 ───────────────────────
# -D 给的是发现管道的一次输出目录, 按它的目录约定找候选与证据表:
#   eve_screen/evescan_query.fasta              新候选集合 (EVE 筛查产出的查询集)
#   10_Reports/eve_candidates.fasta             历史命名, 只有早期目录里才有
#   10_Reports/rescue_evidence_scored.tsv       科属/CheckV/CDD 证据 (s3 读它);
#                                               report 阶段未跑时退回 09b 原位
if [ -n "$DISCOVERY" ]; then
    if [ -z "$QUERY" ]; then
        for c in "$DISCOVERY/eve_screen/evescan_query.fasta" \
                 "$DISCOVERY/10_Reports/eve_candidates.fasta"; do
            if [ -s "$c" ]; then QUERY="$c"; break; fi
        done
    fi
    if [ -z "$EVIDENCE" ]; then
        # 10_Reports 里的是 report 阶段拷过去的副本; 只跑到 analysis_verify 的树里
        # 原件在 09b_Analysis_Verify/virus_validation/ 下 (与 report_pipeline 同一回退)
        for e in "$DISCOVERY/10_Reports/rescue_evidence_scored.tsv" \
                 "$DISCOVERY/09b_Analysis_Verify/virus_validation/rescue_evidence_scored.tsv"; do
            if [ -s "$e" ]; then EVIDENCE="$e"; break; fi
        done
    fi
    [ -n "$OUT" ] || OUT="$DISCOVERY/08b_EVE_Distinguish"
fi
[ -n "$OUT" ] || OUT="$PWD/eve_distinguish_out"

# ── 相对路径统一转绝对, 且必须在这里做 ──
# 脚本稍后会 `cd "$OUT"`, 之后所有相对路径都相对 OUT 解释 —— README 里写的
# `-q candidates.fasta -e 10_Reports/rescue_evidence_scored.tsv` 会在 cd 之后
# 变成在 OUT 下面找 candidates.fasta, 直接 FileNotFoundError 而整批挂掉。
# 与其要求调用者处处写全路径, 不如在改工作目录之前解析掉。
abspath() {
    case "$1" in
        /*|[A-Za-z]:[/\\]*) printf '%s' "$1" ;;   # 已是绝对 (POSIX 或 Windows 盘符)
        *)                  printf '%s/%s' "$PWD" "$1" ;;
    esac
}
for _v in QUERY EVIDENCE BAITS HOST_MAP ASM_ROOT LOCUS_QUERY GENOME_MANIFEST OUT; do
    eval "_cur=\${$_v:-}"
    [ -n "$_cur" ] && eval "$_v=\$(abspath \"\$_cur\")"
done
unset _v _cur

if [ -z "$QUERY" ]; then
    echo "[run_all] 缺少候选 fasta. 用 -q/--query 指定, 或用 -D/--discovery 指向发现管道" >&2
    echo "          的一次输出目录 (候选默认可落在 <DISCOVERY>/eve_screen/evescan_query.fasta" >&2
    echo "          或 <DISCOVERY>/10_Reports/eve_candidates.fasta)" >&2
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
# 只校验这一步真用得上的命令: --no-host 且没开 -G 时 blastn 根本不会被调用,
# 不该因为它缺失而失败
command -v "$DIAMOND" >/dev/null 2>&1 || { echo "[run_all] 找不到命令: $DIAMOND" >&2; exit 1; }
if [ "$NO_HOST" -eq 0 ] || [ -n "$GENOME_MANIFEST" ]; then
    command -v blastn >/dev/null 2>&1 || { echo "[run_all] 找不到命令: blastn (寄主/参考基因组比对需要)" >&2; exit 1; }
fi
# s2c 按样本号 -> 4位代码 -> 清单 归属"自身参考基因组", 没有映射表就一条都归属不了,
# 与其跑完发现全表 no_sample_map, 不如在这里直接拦下
if [ -n "$GENOME_MANIFEST" ]; then
    if [ ! -s "$GENOME_MANIFEST" ]; then
        echo "[run_all] 参考基因组清单不存在或为空: $GENOME_MANIFEST" >&2
        echo "          格式见模块内 genome_manifest.template.tsv; 没有现成基因组就" >&2
        echo "          不要给 -G (通道自动整体跳过, 不影响其余阶段)" >&2
        exit 1
    fi
    if [ ! -s "$HOST_MAP" ]; then
        echo "[run_all] 参考基因组通道需要 -m/--host-map 映射表做样本归属: $HOST_MAP 不存在" >&2
        exit 1
    fi
fi

mkdir -p "$OUT"
cd "$OUT"
echo "[run_all] eve_distinguish v6.3"
echo "[run_all]   (各 stage 的 [done] 横幅写的是该文件逻辑的版本, 不是模块版本)"
echo "[run_all] MOD=$MOD"
echo "[run_all] OUT=$OUT"
echo "[run_all] QUERY=$QUERY ($(grep -c '^>' "$QUERY") 条)"
echo "[run_all] EVIDENCE=$EVIDENCE"
if [ "$NO_HOST" -eq 1 ]; then
    echo "[run_all] 寄主比对已关闭 (--no-host): 不检测 host_contamination_likely"
else
    echo "[run_all] ASM_ROOT=$ASM_ROOT  HOST_MAP=$HOST_MAP"
fi
if [ -n "$GENOME_MANIFEST" ]; then
    echo "[run_all] 参考基因组清单 GENOME_MANIFEST=$GENOME_MANIFEST (s2c 通道开启)"
else
    echo "[run_all] 未给 -G/--genome-manifest: 参考基因组通道整体跳过 (属正常用法)"
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
elif [ -f "$MOD/panel.fasta" ] && [ "$MOD/panel.fasta" -nt panel.fasta ]; then
    # 既有 OUT 目录沿用旧面板是刻意设计 (同目录可复现), 但必须让人知道面板已经换代,
    # 否则新旧 run 混在一个 OUT 里、组件命中数悄悄不同还没人知道为什么。
    echo "!! OUT 里的 panel.fasta 比模块自带的旧 (模块面板已升 v3, Badnavirus 已补):" >&2
    echo "   本次沿用旧面板继续。要用新面板请删掉本目录的 panel.fasta / panel.dmnd 重跑," >&2
    echo "   或换一个 -o 输出目录; 面板换代后组件命中会变, 旧数字不可直接对比。" >&2
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
    n_q=$(grep -c '^>' "$LOCUS_QUERY" || true)
    n_hit=$(cut -f1 locus_blastn.tsv 2>/dev/null | sort -u | grep -c . || true)
    echo "寄主 blastn: $(wc -l < locus_blastn.tsv) 行命中 / $n_host 个寄主 scaffold / $n_q 条查询中 $n_hit 条有命中; 跳过或失败 $n_err 行 (见 locus_errors.log)"
    if [ "$n_host" -eq 0 ]; then
        echo "[run_all] !! 没有任何 scaffold 比对成功, 寄主同源否决整条失效 ——" >&2
        echo "           检查 -A 指向的目录里是否真有 <物种目录>/<代码>-SOAPdenovo-Trans-assembly.fa.gz," >&2
        echo "           以及 HOST_MAP 第三列与目录名是否一致 (locus_errors.log 里的 #NOGZ 行即没找到的)." >&2
    elif [ "$n_q" -gt 0 ] && [ "$n_hit" -lt $((n_q / 2)) ]; then
        # 部分失败比分文不值更危险: 有命中的那部分看着一切正常, 没命中的候选会被
        # 当成"无寄主映射"静默漏掉, 于是寄主否决和位点归并都只覆盖了一小部分候选.
        echo "[run_all] !! 只有 $n_hit/$n_q 条候选有寄主命中: 寄主通道部分失效. 没命中的候选" >&2
        echo "           s2b 会记为 no_host_locus (不参与位点齐全判定), host_contamination_likely" >&2
        echo "           也检不出来. 常见原因: 部分物种目录名与 HOST_MAP 第三列不符(#NOGZ), 或" >&2
        echo "           gunzip/ blastn 只在一部分样本上失败(#GUNZIP_FAIL/#TOOBIG) —— 看 locus_errors.log." >&2
    fi
    $PYTHON "$MOD/s2b_locus_scan.py" locus_blastn.tsv s2_domains.tsv s1_decay.tsv \
        locus_architecture.tsv "$HOST_MAP"
elif [ "$NO_HOST" -eq 1 ]; then
    echo "--no-host: 跳过寄主比对, 位点级架构全部记为 unresolved"
    # 表头由 s2b 自己吐, 免得 bash 与 python 两头各写一份、列数不一致
    # (旧版这里 printf 16 列, 与 s2b 的实际输出对不上, 靠列名取值会静默变空)
    $PYTHON "$MOD/s2b_locus_scan.py" --emit-header > locus_architecture.tsv
else
    echo "!! 未找到 $HOST_MAP, 跳过寄主比对: 位点级架构全部记为 unresolved" >&2
    echo "   (可用 -m/--host-map 指定三列文件: 样本号/4位代码/物种目录名)" >&2
    $PYTHON "$MOD/s2b_locus_scan.py" --emit-header > locus_architecture.tsv
fi

# ──────────── [4b/6] 候选 vs 自身参考基因组 blastn + S2c (可选, -G) ────────────
# 清单三列 "4位代码 <TAB> 注释 <TAB> fasta 路径" (模板 genome_manifest.template.tsv)。
# 没有现成基因组的物种**不写行**: s2c 对它们记 no_ref_genome (本通道跳过, 不产生结论)。
# 比对零命中与"根本没比"要分得开: 每个物种比对完追加 #BLASTED 标记行, s2c 靠它区分
# no_hit (比过零命中, 不定案) 与 genome_absent (没比成, 同样不定案)。
RG_ARG=""
if [ -n "$GENOME_MANIFEST" ]; then
    echo "=== [4b/6] 候选 vs 自身参考基因组 blastn + S2c ==="
    : > refgenome_blastn.tsv
    : > refgenome_errors.log
    n_rg=0
    while IFS=$'\t' read -r gcode gnote gpath || [ -n "$gcode" ]; do
        case "$gcode" in ''|\#*) continue ;; esac
        # 同一代物代码只取第一行: 后面的行多半是笔误或旧版残留, 记日志不静默吞
        if grep -q "^#BLASTED$(printf '\t')${gcode}$" refgenome_blastn.tsv; then
            echo "#DUPLICATE $gcode 清单中重复, 只取第一行 ($gpath)" >> refgenome_errors.log
            continue
        fi
        if [ ! -s "$gpath" ]; then
            echo "#NOGENOME $gcode $gpath" >> refgenome_errors.log
            continue
        fi
        # .gz 解压到临时文件 (与 [4/6] 同一个 OOM 教训: blastn 无法 seek 管道);
        # 纯 fasta 直接用, 不做无谓拷贝
        subj=""
        case "$gpath" in
            *.gz|*.bgz)
                subj=$(mktemp "${TMPDIR:-/tmp}/everg_${gcode}_XXXXXX.fa")
                if ! gunzip -c "$gpath" > "$subj" 2>>refgenome_errors.log; then
                    echo "#GUNZIP_FAIL $gcode $gpath" >> refgenome_errors.log
                    rm -f "$subj"; continue
                fi
                sz_gb=$(( $(stat -c%s "$subj" 2>/dev/null || echo 0) / 1073741824 ))
                if [ "$sz_gb" -gt "$SUBJ_CAP_GB" ]; then
                    echo "#TOOBIG $gcode $gpath ${sz_gb}GB > ${SUBJ_CAP_GB}GB" >> refgenome_errors.log
                    rm -f "$subj"; continue
                fi
                ;;
            *) subj="$gpath" ;;
        esac
        tmp_hits=$(mktemp "${TMPDIR:-/tmp}/everg_hits_${gcode}_XXXXXX.tsv")
        if ! blastn -query "$LOCUS_QUERY" -subject "$subj" \
               -task megablast -evalue 1e-5 -outfmt "$FMT" -num_threads 1 \
               -max_target_seqs 200 -out "$tmp_hits" 2>>refgenome_errors.log; then
            echo "#BLASTN_FAIL $gcode $gpath" >> refgenome_errors.log
            rm -f "$tmp_hits"; [ "$subj" != "$gpath" ] && rm -f "$subj"
            continue
        fi
        awk -v c="$gcode" 'BEGIN{OFS="\t"}{print c, $0}' "$tmp_hits" >> refgenome_blastn.tsv
        rm -f "$tmp_hits"; [ "$subj" != "$gpath" ] && rm -f "$subj"
        printf '#BLASTED\t%s\n' "$gcode" >> refgenome_blastn.tsv
        n_rg=$((n_rg + 1))
    done < "$GENOME_MANIFEST"
    n_rghit=$(grep -cv '^#' refgenome_blastn.tsv || true)
    echo "参考基因组 blastn: $n_rg 个基因组比对完成, $n_rghit 行命中" \
         "(失败/跳过见 refgenome_errors.log)"
    if [ "$n_rg" -eq 0 ]; then
        echo "[run_all] !! 清单里没有一个基因组可比 (#NOGENOME/#TOOBIG/#GUNZIP_FAIL)," >&2
        echo "           s2c 将把所有候选记 no_ref_genome, 本通道整体空转" >&2
    fi
    $PYTHON "$MOD/s2c_refgenome_scan.py" refgenome_blastn.tsv s2_domains.tsv \
        "$HOST_MAP" refgenome_evidence.tsv
    RG_ARG="refgenome_evidence.tsv"
fi

# ─────────────────────── [5/6] S3 合成判别 ───────────────────────
echo "=== [5/6] S3 合成判别 ==="
# $RG_ARG 不带引号: 为空时整体不出现 (未开 -G 的运行与 v6.2 逐字节同参)
$PYTHON "$MOD/s3_verdict.py" s1_decay.tsv s2_domains.tsv locus_architecture.tsv "$EVIDENCE" \
    eve_distinguish_verdict.tsv $RG_ARG

# ─────────────────── [6/6] verdict -> action ───────────────────
echo "=== [6/6] verdict -> 可执行动作 ==="
$PYTHON "$MOD/s4_filter.py" eve_distinguish_verdict.tsv -o dna_vs_eve_filter.tsv || {
    rc=$?; echo "[run_all] s4_filter 返回 $rc: verdict 表里有未映射的新值, 已按 REVIEW 兜底" >&2; }

echo
echo "全部完成 -> $OUT"
echo "  逐候选判别  $OUT/eve_distinguish_verdict.tsv"
echo "  下游过滤表  $OUT/dna_vs_eve_filter.tsv"
if [ -n "$GENOME_MANIFEST" ]; then
    echo "  参考基因组证据 $OUT/refgenome_evidence.tsv (rg_call 列已并入上面两表的 verdict)"
fi
