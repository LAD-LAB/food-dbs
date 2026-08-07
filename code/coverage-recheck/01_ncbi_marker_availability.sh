#!/bin/bash
# Query NCBI nucleotide for marker availability, one species per line.
# Usage: ncbi_recheck.sh <species_list.txt> <marker> <out.csv>
# Resumable: skips species already present in <out.csv>.
set -u
IN="$1"; MARKER="$2"; OUT="$3"
EMAIL="sarahbrown2648@gmail.com"
TOOL="foodseq-coverage-recheck"
EUTILS="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"

[ -f "$OUT" ] || echo "scientific_name,marker,ncbi_hits" > "$OUT"

total=$(grep -cve '^\s*$' "$IN")
n=0
while IFS= read -r sp; do
  [ -z "$sp" ] && continue
  n=$((n+1))
  # resume: skip if already recorded with a numeric hit count
  if grep -qF "\"$sp\",$MARKER," "$OUT" 2>/dev/null; then continue; fi

  term="\"$sp\"[Organism] AND $MARKER"
  enc=$(printf '%s' "$term" | sed 's/"/%22/g; s/\[/%5B/g; s/\]/%5D/g; s/ /+/g')
  url="$EUTILS?db=nucleotide&retmax=0&email=$EMAIL&tool=$TOOL&term=$enc"

  count=""
  for try in 1 2 3; do
    resp=$(curl -s --max-time 30 "$url")
    count=$(printf '%s' "$resp" | grep -o '<Count>[0-9]*</Count>' | head -1 | grep -o '[0-9]*')
    if [ -n "$count" ]; then break; fi
    sleep 2
  done

  if [ -z "$count" ]; then
    echo "\"$sp\",$MARKER,ERR" >> "$OUT.errors"
  else
    echo "\"$sp\",$MARKER,$count" >> "$OUT"
  fi
  sleep 0.4
done < "$IN"

echo "DONE: $MARKER — processed $n/$total, results in $OUT" >> "$OUT.log"
