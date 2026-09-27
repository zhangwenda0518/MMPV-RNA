#!/bin/bash
# Read-only dmp structure inspector
F=/home/zhangwenda/database/taxonomy/rankedlineage.dmp
echo "== NF_first (fields when split on [TAB|TAB])"
head -1 "$F" | awk -F'\t\|\t' '{print NF}'
echo "== LINE1 fields"
head -1 "$F" | awk -F'\t\|\t' '{for(i=1;i<=NF;i++) printf("%d:[%s] ",i,$i); print ""}'
echo "== ARCHAEA line fields"
sed -n 3p "$F" | awk -F'\t\|\t' '{for(i=1;i<=NF;i++) printf("%d:[%s] ",i,$i); print ""}'
echo "== NARNAVIRIDAE first line fields"
grep -m1 Narnaviridae "$F" | awk -F'\t\|\t' '{for(i=1;i<=NF;i++) printf("%d:[%s] ",i,$i); print ""}'
echo "== count of lines containing Narnaviridae"
grep -c Narnaviridae "$F"
echo "== distinct nonempty species-col approx (first 200000 lines)"
head -200000 "$F" | awk -F'\t\|\t' '$3!=""{print $3}' | sort -u | wc -l
echo "== distinct nonempty superkingdom-col (first 200000 lines)"
head -200000 "$F" | awk -F'\t\|\t' '$10!=""{c[$10]++} END{for(k in c) print k, c[k]}'
echo "== tail test: last line NF"
tail -1 "$F" | awk -F'\t\|\t' '{print NF}'
