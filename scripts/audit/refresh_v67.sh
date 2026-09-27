#!/bin/bash
# refresh_v67.sh <label>
# 门槛新语义（>=0.5，半数票也淘汰，gate 6.7）下重跑单个数据集的共识分类，并刷新其派生报告。
# 原则：只动本数据集目录；覆盖前一律 cp -a 备份；日志落 /tmp/refresh_v67/<label>/；只读脚本不动。
# 用法: bash /tmp/refresh_v67.sh barbarum|Alternaria|Aphis|Fusarium|amarum|chinense|ruthenicum|onekp
set -u
LABEL="${1:-}"
ROOT=/home/zhangwenda/MMPV-paper
case "$LABEL" in
  Alternaria) D=$ROOT/goji-virome/02_novel_virus/RNA-Alternaria_alternata_out ;;
  Aphis)      D=$ROOT/goji-virome/02_novel_virus/RNA-Aphis_gossypii_out ;;
  Fusarium)   D=$ROOT/goji-virome/02_novel_virus/RNA-Fusarium_nematophilum_out ;;
  amarum)     D=$ROOT/goji-virome/02_novel_virus/RNA-Lycium_amarum_out ;;
  barbarum)   D=$ROOT/goji-virome/02_novel_virus/RNA-Lycium_barbarum_out ;;
  chinense)   D=$ROOT/goji-virome/02_novel_virus/RNA-Lycium_chinense_out ;;
  ruthenicum) D=$ROOT/goji-virome/02_novel_virus/RNA-Lycium_ruthenicum_out ;;
  onekp)      D=$ROOT/onekp-virome/onekp-virus ;;
  *) echo "unknown label: $LABEL"; exit 2 ;;
esac
DS=$(basename "$D")
LOG=/tmp/refresh_v67/$LABEL
mkdir -p "$LOG"
STAMP=20260916
INT="$D/05_Taxonomy/Votus.integrated"
PROD="$INT/final_integrated_classification.tsv"
BAK="$D/05_Taxonomy/Votus.integrated.bak_v66f_$STAMP"
SCRIPT=/home/zhangwenda/MMPV-RNA/virome_discovery_pipeline/virus_classifier_analysis.R
REF=/home/zhangwenda/database/taxonomy/genus_family_ref.tsv
COMB="$D/05_Taxonomy/Votus.classed/Votus_combined_taxonomy.tsv"

echo "=== $LABEL ($DS) START $(date '+%F %T') ==="
echo "根目录: $D"

# ---- 0) 前置检查 ----
[ -f "$COMB" ] || { echo "FATAL 缺 combined: $COMB"; exit 3; }
[ -f "$PROD" ] || { echo "FATAL 缺成品: $PROD"; exit 3; }
if [ -e "$BAK" ]; then echo "FATAL 备份已存在(不覆盖): $BAK"; exit 3; fi

# ---- 1) 备份 + 基线 ----
OLD=$LOG/old.tsv
cp -p "$PROD" "$OLD"
cp -a "$INT" "$BAK" || { echo "FATAL 备份失败"; exit 3; }
echo "备份 -> $BAK"
echo "旧: 行数=$(($(wc -l < "$OLD") - 1)) md5=$(md5sum "$OLD" | cut -c1-32)"

# ---- 2) R 重跑（就地覆盖）----
( cd "$INT" && MMPV_GENUS_FAMILY_REF="$REF" Rscript "$SCRIPT" --combined "$COMB" --output "$INT" ) > "$LOG/r.log" 2>&1
RC=$?
NEWLINE=$(($(wc -l < "$PROD") - 1))
NEWMD5=$(md5sum "$PROD" | cut -c1-32)
echo "R rc=$RC 新: 行数=$NEWLINE md5=$NEWMD5"
grep -e 逐级淘汰 "$LOG/r.log" | tail -2
if [ $RC -ne 0 ]; then echo "FATAL R 失败, 见 $LOG/r.log; 备份在 $BAK"; exit 4; fi

# ---- 3) 自检 ----
echo "--- gate stamp ---"
grep -e gate_version -e script_md5 -e script_sha256 "$INT/taxonomy_gate_stamp.tsv" 2>/dev/null
echo "--- gate check ---"
grep -e dual_ref_conflict "$INT/taxonomy_gate_check.tsv" 2>/dev/null

# ---- 4) 一致率对照（旧 vs 新；chimera_multi 的 online=新 / v66f=旧）----
MAN=$LOG/manifest.tsv
printf '%s\t%s\t%s\n' "$LABEL" "$PROD" "$OLD" > "$MAN"
echo "--- chimera_multi: online=新产物, v66f=旧备份 ---"
python3 /tmp/chimera_multi.py "$MAN" 2>&1 | tee "$LOG/chimera.txt" | tail -6
echo "--- cmp_thr_vs_v66f (old -> new) ---"
python3 /tmp/cmp_thr_vs_v66f.py "$OLD" "$PROD" 2>&1 | tee "$LOG/cmp.txt"

# ---- 5) Stage 10 报告刷新 ----
if [ -d "$D/10_Reports" ] && [ ! -e "$D/10_Reports.bak_v66f_$STAMP" ]; then
  cp -a "$D/10_Reports" "$D/10_Reports.bak_v66f_$STAMP" && echo "备份 -> $D/10_Reports.bak_v66f_$STAMP"
fi
python3 /home/zhangwenda/MMPV-RNA/virome_discovery_pipeline/report_pipeline.py -o "$D" > "$LOG/report.log" 2>&1
echo "report rc=$? 10_Reports=$(du -sh "$D/10_Reports" | cut -f1)"
tail -4 "$LOG/report.log"

# ---- 6) rescue 报告刷新 ----
if [ -d "$D/08_Rescue/Plant" ]; then
  python3 /tmp/regen_rescue_reports.py "$DS" > "$LOG/rescue.log" 2>&1
  echo "rescue rc=$?"
  tail -4 "$LOG/rescue.log"
else
  echo "rescue 跳过 (无 08_Rescue/Plant)"
fi

echo "=== $LABEL DONE $(date '+%F %T') 日志目录 $LOG ==="
