#!/bin/bash
# For each species, count 12S records restricted to Vertebrata.
# 0 => not a vertebrate (12SV5 primers will not amplify it).
# Usage: vert_check.sh <species_list.txt> <out.csv>   (resumable)
set -u
IN="$1"; OUT="$2"
EMAIL="sarahbrown2648@gmail.com"
TOOL="foodseq-vert-check"
EUTILS="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"

[ -f "$OUT" ] || echo "scientific_name,vert_12S_hits" > "$OUT"

while IFS= read -r sp; do
  [ -z "$sp" ] && continue
  if grep -qF "\"$sp\"," "$OUT" 2>/dev/null; then continue; fi

  term="\"$sp\"[Organism] AND 12S AND Vertebrata[Organism]"
  enc=$(printf '%s' "$term" | sed 's/"/%22/g; s/\[/%5B/g; s/\]/%5D/g; s/ /+/g')
  url="$EUTILS?db=nucleotide&retmax=0&email=$EMAIL&tool=$TOOL&term=$enc"

  count=""
  for try in 1 2 3; do
    resp=$(curl -s --max-time 30 "$url")
    count=$(printf '%s' "$resp" | grep -o '<Count>[0-9]*</Count>' | head -1 | grep -o '[0-9]*')
    [ -n "$count" ] && break
    sleep 2
  done

  if [ -z "$count" ]; then
    echo "\"$sp\"" >> "$OUT.errors"
  else
    echo "\"$sp\",$count" >> "$OUT"
  fi
  sleep 0.4
done < "$IN"

echo "DONE vert check" >> "$OUT.log"
