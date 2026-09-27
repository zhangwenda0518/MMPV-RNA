#!/bin/bash
# Read-only schema / contig-set comparison, new vs old product. Temp only under /tmp/indep_verify.
set -u
OUT=/tmp/indep_verify/schema.txt
TMP=/tmp/indep_verify
: > "$OUT"
for d in Alternaria amarum Aphis barbarum chinense Fusarium onekp ruthenicum; do
  p=$(cut -f2 /tmp/refresh_v67/$d/manifest.tsv)
  o=/tmp/refresh_v67/$d/old.tsv
  {
    echo "==== $d"
    hn=$(head -1 "$p"); ho=$(head -1 "$o")
    echo "ncol_new=$(printf '%s' "$hn" | awk -F'\t' '{print NF}') ncol_old=$(printf '%s' "$ho" | awk -F'\t' '{print NF}')"
    if [ "$hn" = "$ho" ]; then echo "header_same=YES"; else echo "header_same=NO"; fi
    echo "has_Nucleic_acid_new=$(printf '%s' "$hn" | grep -c Nucleic_acid || true)"
    nn=$(tail -n +2 "$p" | wc -l); no=$(tail -n +2 "$o" | wc -l)
    echo "rows_new=$nn rows_old=$no"
    tail -n +2 "$p" | cut -f1 | sort > "$TMP/$d.new.ids"
    tail -n +2 "$o" | cut -f1 | sort > "$TMP/$d.old.ids"
    echo "ids_new_lines=$(wc -l < "$TMP/$d.new.ids") ids_new_uniq=$(sort -u "$TMP/$d.new.ids" | wc -l) ids_old_lines=$(wc -l < "$TMP/$d.old.ids") ids_old_uniq=$(sort -u "$TMP/$d.old.ids" | wc -l)"
    if diff -q "$TMP/$d.new.ids" "$TMP/$d.old.ids" > /dev/null; then
      echo "contig_set_equal=YES"
    else
      echo "contig_set_equal=NO new_only=$(comm -23 "$TMP/$d.new.ids" "$TMP/$d.old.ids" | wc -l) old_only=$(comm -13 "$TMP/$d.new.ids" "$TMP/$d.old.ids" | wc -l)"
      echo "--- first diffs"
      diff "$TMP/$d.new.ids" "$TMP/$d.old.ids" | head -20
    fi
  } >> "$OUT"
done
cat "$OUT"
