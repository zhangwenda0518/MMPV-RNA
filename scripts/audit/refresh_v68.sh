#!/bin/bash
# refresh_v68.sh <label|all>
# 上游 fill 修正 (fill-fix1) + gate 6.7 重跑: 重建 combined -> 重跑 R -> 校验 gate -> 刷新报告/rescue
# 原则：只动本数据集目录；覆盖前一律 cp -a 备份；日志落 /tmp/refresh_v68/<label>/；只读脚本不动。
# 用法: bash /tmp/refresh_v68.sh all|Alternaria|Aphis|Fusarium|amarum|barbarum|chinense|ruthenicum|onekp
set -u
ROOT=/home/zhangwenda/MMPV-paper
PIPE=/home/zhangwenda/MMPV-RNA/virome_discovery_pipeline
STAMP=20260916

run_one() {
LABEL="$1"
case "$LABEL" in
  Alternaria) D=$ROOT/goji-virome/02_novel_virus/RNA-Alternaria_alternata_out ;;
  Aphis)      D=$ROOT/goji-virome/02_novel_virus/RNA-Aphis_gossypii_out ;;
  Fusarium)   D=$ROOT/goji-virome/02_novel_virus/RNA-Fusarium_nematophilum_out ;;
  amarum)     D=$ROOT/goji-virome/02_novel_virus/RNA-Lycium_amarum_out ;;
  barbarum)   D=$ROOT/goji-virome/02_novel_virus/RNA-Lycium_barbarum_out ;;
  chinense)   D=$ROOT/goji-virome/02_novel_virus/RNA-Lycium_chinense_out ;;
  ruthenicum) D=$ROOT/goji-virome/02_novel_virus/RNA-Lycium_ruthenicum_out ;;
  onekp)      D=$ROOT/onekp-virome/onekp-virus ;;
  *) echo "unknown label: $LABEL"; return 2 ;;
esac
DS=$(basename "$D")
LOG=/tmp/refresh_v68/$LABEL
mkdir -p "$LOG"
INT="$D/05_Taxonomy/Votus.integrated"
PROD="$INT/final_integrated_classification.tsv"
BAK="$D/05_Taxonomy/Votus.integrated.bak_fillfix_$STAMP"
SCRIPT=$PIPE/virus_classifier_analysis.R
REF=/home/zhangwenda/database/taxonomy/genus_family_ref.tsv
CLASSED="$D/05_Taxonomy/Votus.classed"
COMB="$CLASSED/Votus_combined_taxonomy.tsv"
CBBAK="$CLASSED/Votus_combined_taxonomy.tsv.bak_fillfix_$STAMP"

echo "=== $LABEL ($DS) START $(date '+%F %T') ==="
[ -f "$COMB" ] || { echo "FATAL 缺 combined: $COMB"; return 3; }
[ -f "$PROD" ] || { echo "FATAL 缺成品: $PROD"; return 3; }
if [ -e "$BAK" ] || [ -e "$CBBAK" ]; then echo "FATAL 备份已存在(不覆盖)"; return 3; fi

# ---- 1) 备份 + 基线 ----
OLD=$LOG/old.tsv
cp -p "$PROD" "$OLD"
cp -a "$INT" "$BAK" || { echo "FATAL 备份失败"; return 3; }
cp -p "$COMB" "$CBBAK"
echo "备份 -> $BAK ; $CBBAK"
echo "旧: 成品行数=$(($(wc -l < "$OLD") - 1)) md5=$(md5sum "$OLD" | cut -c1-32) combined行数=$(($(wc -l < "$CBBAK") - 1)) md5=$(md5sum "$CBBAK" | cut -c1-32)"

# ---- 2) combined 重建 (merge + fill-fix1) ----
python3 /tmp/rerun_merge_fill.py --module $PIPE/virus_classifier.py --classed "$CLASSED" --outdir "$CLASSED" --baseline "$CBBAK" > "$LOG/merge.log" 2>&1
MRC=$?
echo "merge rc=$MRC 新 combined行数=$(($(wc -l < "$COMB") - 1)) md5=$(md5sum "$COMB" | cut -c1-32)"
grep -e '\[fill\]' "$LOG/merge.log" | tail -2
if [ $MRC -ne 0 ]; then echo "FATAL merge 失败, 见 $LOG/merge.log; 备份在 $CBBAK"; return 4; fi
python3 /tmp/diff_combined.py "$CBBAK" "$COMB" --top 14 --quiet > "$LOG/cmp_combined.txt" 2>&1
head -8 "$LOG/cmp_combined.txt"

# ---- 3) R 重跑（就地覆盖成品）----
( cd "$INT" && MMPV_GENUS_FAMILY_REF="$REF" Rscript "$SCRIPT" --combined "$COMB" --output "$INT" ) > "$LOG/r.log" 2>&1
RC=$?
NEWLINE=$(($(wc -l < "$PROD") - 1))
NEWMD5=$(md5sum "$PROD" | cut -c1-32)
echo "R rc=$RC 新成品: 行数=$NEWLINE md5=$NEWMD5"
if [ $RC -ne 0 ]; then echo "FATAL R 失败, 见 $LOG/r.log; 备份在 $BAK"; return 5; fi

# ---- 4) gate 自检 ----
echo "--- gate stamp ---"
grep -e gate_version -e script_md5 -e script_sha256 "$INT/taxonomy_gate_stamp.tsv" 2>/dev/null
echo "--- gate check ---"
grep -e dual_ref_conflict -e species_dual_ref_conflict "$INT/taxonomy_gate_check.tsv" 2>/dev/null
python3 /tmp/diff_combined.py "$OLD" "$PROD" --top 10 --quiet > "$LOG/cmp_prod.txt" 2>&1
echo "--- 成品 old->new ---"
head -5 "$LOG/cmp_prod.txt"

# ---- 5) 一致率对照 ----
MAN=$LOG/manifest.tsv
printf '%s\t%s\t%s\n' "$LABEL" "$PROD" "$OLD" > "$MAN"
python3 /tmp/chimera_multi.py "$MAN" 2>&1 | tee "$LOG/chimera.txt" | tail -5

# ---- 6) Stage 10 报告刷新 ----
if [ -d "$D/10_Reports" ] && [ ! -e "$D/10_Reports.bak_fillfix_$STAMP" ]; then
  cp -a "$D/10_Reports" "$D/10_Reports.bak_fillfix_$STAMP" && echo "备份 -> $D/10_Reports.bak_fillfix_$STAMP"
fi
python3 $PIPE/report_pipeline.py -o "$D" > "$LOG/report.log" 2>&1
echo "report rc=$? 10_Reports=$(du -sh "$D/10_Reports" 2>/dev/null | cut -f1)"

# ---- 7) rescue 报告刷新 ----
if [ -d "$D/08_Rescue/Plant" ]; then
  python3 /tmp/regen_rescue_reports.py "$DS" > "$LOG/rescue.log" 2>&1
  echo "rescue rc=$?"
else
  echo "rescue 跳过 (无 08_Rescue/Plant)"
fi

echo "=== $LABEL DONE $(date '+%F %T') 日志 $LOG ==="
return 0
}

if [ "${1:-}" = "all" ]; then
  for L in Alternaria Aphis Fusarium amarum barbarum chinense ruthenicum onekp; do
    run_one "$L"
  done
else
  run_one "${1:-}"
fi
