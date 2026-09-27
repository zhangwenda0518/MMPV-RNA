#!/bin/bash
# build_benchmark.sh — 判别器的策展阳性/阴性集 (判据灵敏度唯一缺的一环)
# ============================================================================
#
# 为什么需要它: `verify/` 那四个工具只能测**特异性**。端到端空模型把信号全破坏了,
# 它证明"随机序列不会触发", 证明不了"真病毒会被抓住"。要测灵敏度/召回, 必须有
# **标签已知**的序列。本脚本用 NCBI 上的策展资源凑出三类:
#
#   A 现存病毒 (extant_virus)   完整基因组、编码完好、**未整合** —— 生物上不是 EVE。
#                               期望: 不进 ancient_EVE / EVE_STRONG / EVE_LTR_TE;
#                               应落 coding_intact 一侧 (virus_candidate / fragment)。
#   B 文献 EPRV (eprv)          已被文献报道**整合进植物基因组**的元件。分两档:
#                               B1 重排/退化型 (非侵染性等位) -> 这才是干净的 EVE 阳性;
#                               B2 可激活型 (可侵染) -> **标签有歧义**: 它既是整合元件
#                                  又编码完好, 判定为"像病毒"在生物学上并不算错。
#                                  所以 B2 单列, 不混进"必须判成 EVE"的分母。
#   C 带病毒样结构域的寄主基因 (host_viral_domain)  端粒酶(TERT, 含 RT 域)、
#                                RNase H / dUTPase 类寄主基因 —— 生物上不是病毒也
#                                不是 EVE。期望: 不进任何 EVE/virus 结论。
#
# 第四类 TE(转座子) 不在这里下载: 服务器上已有策展库 (Dfam3.9_plants 16,304 个植物
# 重复家族 + REXdb/GyDB2.hmm 的 LTR 反转座子域模型), 而且 pipeline 自己的 baits
# 就是 Caulifinder banks —— 直接用它们更贴近实际运行条件。
#
# 用法:
#   bash build_benchmark.sh --out DIR [--dry-run]
# 依赖: efetch / esearch (entrez-direct) 与网络
set -euo pipefail

OUT=""; DRY=0
while [ $# -gt 0 ]; do
    case "$1" in
        --out) OUT="$2"; shift 2;;
        --dry-run) DRY=1; shift;;
        *) echo "未知参数: $1" >&2; exit 2;;
    esac
done
[ -n "$OUT" ] || { echo "需要 --out DIR" >&2; exit 2; }
mkdir -p "$OUT"

fetch() {   # fetch <输出文件> <nuccore 查询串>
    local f="$1"; shift
    local q="$1"
    if [ "$DRY" -eq 1 ]; then
        local n
        n=$(esearch -db nuccore -query "$q" 2>/dev/null | grep -oE '<Count>[0-9]+</Count>' \
            | head -1 | grep -oE '[0-9]+')
        printf "  %-34s %4s 条\t%s\n" "$(basename "$f")" "${n:-0}" "$q"
        return 0
    fi
    esearch -db nuccore -query "$q" 2>/dev/null \
        | efetch -format fasta 2>/dev/null > "$f" || true
    printf "  %-34s %5s 条\n" "$(basename "$f")" "$(grep -c '^>' "$f" || echo 0)"
}

fetch_acc() {   # fetch_acc <输出文件> <登录号列表>
    local f="$1"; shift
    if [ "$DRY" -eq 1 ]; then printf "  %-34s 按登录号: %s\n" "$(basename "$f")" "$*"; return 0; fi
    efetch -db nuccore -id "$*" -format fasta 2>/dev/null > "$f" || true
    printf "  %-34s %5s 条\n" "$(basename "$f")" "$(grep -c '^>' "$f" || echo 0)"
}

echo "=== A 现存病毒: Caulimoviridae 完整基因组 (RefSeq, 按属) ==="
for g in Caulimovirus Soymovirus Cavemovirus Petuvirus Tungrovirus Badnavirus Solendovirus; do
    fetch "$OUT/A_extant_$g.fna" "$g[Organism] AND \"complete genome\"[Title] AND srcdb_refseq[PROP]"
done
# 同属外的 DNA 病毒也在候选范围内 —— 它们同样"编码完好 + 未整合"
fetch "$OUT/A_extant_Geminiviridae.fna" "Geminiviridae[Organism] AND \"complete genome\"[Title] AND srcdb_refseq[PROP]"

echo
echo "=== B1 文献 EPRV: 重排/退化型整合元件 (干净的 EVE 阳性) ==="
# NCBI 上没有 Florendovirus 这个 organism 名 (ICTV 属名未被索引), 所以走"标题/全字段
# 里带 endogenous 的 Caulimoviridae"这几条路子 —— 实测能凑到 ~120 条。
# Musa balbisiana BAC 克隆里带着 BSGFV 的整合等位 (EPRV-7/9), 是基因组上下文完整的
# 整合元件, 比单独一条病毒参考更能代表"真的整合进去了"。
fetch "$OUT/B1_eprv_endogenous_title.fna" "endogenous[Title] AND Caulimoviridae[Organism]"
fetch "$OUT/B1_eprv_endogenous_any.fna" "endogenous pararetrovirus"
fetch "$OUT/B1_eprv_musa_bac.fna" "Musa balbisiana[Organism] AND bac[Title]"
fetch_acc "$OUT/B1_eprv_musa_integrant.fna" "AP009325 AP009326"

echo
echo "=== B1b 文献 EPRV: 糖甜菜个体基因组拷贝 (退化型, 带侧翼上下文) ==="
# 为什么单列这一类: B1 里从 NCBI 取到的 EPRV 大多是**完整/大片段**的整合元件
# (RefSeq complete genome 提交、BAC 克隆), 也就是"完好/可激活"那一档 —— 而本方法的
# 设计目标是**退化型**元件 ("不全"不算 EVE 证据, 是 v4 那个倒置的核心)。把两种难度
# 混在一个分母里, 灵敏度就没有意义。这一类是已发表研究在一个真实植物基因组里
# 全基因组筛出的个体拷贝, 带染色体位置与方向, 才是"目标类别"。
# 来源: 糖甜菜 (Beta vulgaris) EPRV 研究, Zenodo 3888270, CC 授权, 直接文件 URL。
ZEN="https://zenodo.org/api/records/3888270/files"
beet_fetch() {   # beet_fetch <输出> <URL 编码后的文件名>
    local out="$1" name="$2"
    if [ "$DRY" -eq 1 ]; then printf "  %-34s (zenodo)\n" "$(basename "$out")"; return 0; fi
    if curl -sSL --max-time 120 "$ZEN/$name/content" -o "$out.raw"; then
        # 这几个文件是**比对结果** (带 gap)。带 gap 的序列喂 blastx 会读成移码,
        # 必须先去掉 gap —— 这一步不能省, 否则测的是"方法对带 gap 输入的反应"。
        seqkit seq -g -w 0 "$out.raw" > "$out" 2>/dev/null || cp "$out.raw" "$out"
        rm -f "$out.raw"
        printf "  %-34s %5s 条  (去 gap 后)\n" "$(basename "$out")" "$(grep -c '^>' "$out" || echo 0)"
    else
        printf "  %-34s 下载失败\n" "$(basename "$out")"
    fi
}
beet_fetch "$OUT/B1b_eprv_decayed_beet1.fna" \
    "S2_Multiple%20sequence%20alignment%20of%2027%20beetEPRV1%20sequences.fasta"
beet_fetch "$OUT/B1b_eprv_decayed_beet2a.fna" \
    "S3_Multiple%20sequence%20alignment%20of%2042%20beetEPRV2%20component%20A%20sequences.fasta"
beet_fetch "$OUT/B1b_eprv_decayed_beet2b.fna" \
    "S4_Multiple%20sequence%20alignment%20of%2023%20beetEPRV2%20component%20B%20sequences.fasta"
beet_fetch "$OUT/B1b_eprv_decayed_beet3.fna" \
    "S5_Multiple%20sequence%20alignment%20of%2011%20beetEPRV3%20sequences.fasta"

echo
echo "=== B2 可激活型 EPRV (标签有歧义, 单列) ==="
fetch "$OUT/B2_eprv_activatable.fna" "Tobacco vein clearing virus[Organism] AND \"complete genome\"[Title]"
fetch "$OUT/B2_eprv_activatable2.fna" "Petunia vein clearing virus[Organism] AND \"complete genome\"[Title]"
fetch "$OUT/B2_eprv_activatable3.fna" "Banana streak virus[Organism] AND \"complete genome\"[Title]"

echo
echo "=== C 带病毒样结构域的寄主基因 (阴性) ==="
# 端粒酶含逆转录酶域, 是"RT 域但完全不是病毒"的经典陷阱
fetch "$OUT/C_host_telomerase.fna" \
      "telomerase reverse transcriptase[Title] AND plants[Organism] AND srcdb_refseq[PROP]"
# 寄主自己的 RNase H 域蛋白 / dUTPase (v6.1 已确认这三类是寄主基因也带的通用域)
fetch "$OUT/C_host_rnaseh.fna" "ribonuclease H[Title] AND Viridiplantae[Organism] AND srcdb_refseq[PROP]"
fetch "$OUT/C_host_dutpase.fna" "dUTPase[Title] AND Viridiplantae[Organism] AND srcdb_refseq[PROP]"

if [ "$DRY" -eq 0 ]; then
    echo
    echo "=== 汇总 ==="
    for f in "$OUT"/*.fna; do printf "  %-40s %6s 条\n" "$(basename "$f")" "$(grep -c '^>' "$f" || echo 0)"; done
    cat "$OUT"/A_*.fna > "$OUT/A_extant_virus.fna" 2>/dev/null || true
    cat "$OUT"/B1_eprv*.fna > "$OUT/B1_eprv_clean.fna" 2>/dev/null || true
    cat "$OUT"/B1b_*.fna > "$OUT/B1b_eprv_decayed.fna" 2>/dev/null || true
    cat "$OUT"/B2_*.fna > "$OUT/B2_eprv_activatable_all.fna" 2>/dev/null || true
    cat "$OUT"/C_*.fna > "$OUT/C_host_viral_domain.fna" 2>/dev/null || true
    echo "  ── 合并后 ──"
    for c in A_extant_virus B1_eprv_clean B1b_eprv_decayed B2_eprv_activatable_all C_host_viral_domain; do
        printf "  %-40s %6s 条\n" "$c.fna" "$(grep -c '^>' "$OUT/$c.fna" || echo 0)"
    done
fi
