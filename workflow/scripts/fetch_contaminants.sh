#!/usr/bin/env bash
# Download the spike-in contaminant genomes listed in config/contaminants.tsv from NCBI nuccore.
#
#   workflow/scripts/fetch_contaminants.sh [--force] [tsv]   (default: config/contaminants.tsv)
#
# The TSV is the single source of truth (columns: taxon, group, reference, accessions).
# Each row's comma-separated nuccore accessions are fetched into its `reference` path; an optional
# ':<copies>' suffix (copy number, used by the read simulators) is ignored here.
# Rows whose reference already exists as a non-empty FASTA are skipped (--force re-downloads them);
# rows with no accessions must point to a local FASTA, which is never overwritten. Run from the repository root. Set NCBI_API_KEY to raise the rate limit.
set -euo pipefail

force=0
if [[ ${1:-} == --force || ${1:-} == -f ]]; then force=1; shift; fi
tsv=${1:-config/contaminants.tsv}
EUTILS=https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi
key=${NCBI_API_KEY:+&api_key=$NCBI_API_KEY}

col() { head -1 "$tsv" | tr '\t' '\n' | grep -nx "$1" | cut -d: -f1; }
c_taxon=$(col taxon) c_ref=$(col reference) c_acc=$(col accessions)
[[ -n $c_taxon && -n $c_ref && -n $c_acc ]] || { echo "ERROR: $tsv needs columns taxon, reference, accessions" >&2; exit 1; }

# cut keeps empty fields (bash `read` would merge consecutive tabs and shift the columns)
field() { cut -f"$1" <<< "$2"; }

tail -n +2 "$tsv" | while IFS= read -r line; do
    taxon=$(field "$c_taxon" "$line") fa=$(field "$c_ref" "$line")
    ids=$(field "$c_acc" "$line" | sed -E 's/:[0-9.]+//g; s/ //g')
    [[ -n $taxon ]] || continue
    if [[ $force -eq 0 && -s $fa && $(head -c1 "$fa") == ">" ]]; then
        echo "skip  $taxon ($fa exists)"; continue
    fi
    if [[ -z $ids ]]; then
        [[ -s $fa ]] && { echo "keep  $taxon (local file, no accessions)"; continue; }
        echo "ERROR: $taxon: $fa missing and no accessions given" >&2; exit 1
    fi
    echo "fetch $taxon: $ids"
    mkdir -p "$(dirname "$fa")"
    curl -sSf --retry 5 --retry-delay 3 -o "$fa.tmp" "$EUTILS?db=nuccore&id=$ids&rettype=fasta&retmode=text$key"
    n=$(grep -c '^>' "$fa.tmp" || true)
    expected=$(tr ',' '\n' <<< "$ids" | grep -c .)
    if [[ $n -ne $expected ]]; then
        echo "ERROR: $taxon: got $n sequences, expected $expected" >&2
        rm -f "$fa.tmp"; exit 1
    fi
    mv "$fa.tmp" "$fa"
    sleep 0.4    # stay under NCBI's 3 requests/s without an API key
done

echo "contaminant references:"
tail -n +2 "$tsv" | while IFS= read -r line; do
    printf "  %-34s %s bp\n" "$(field "$c_taxon" "$line")" "$(grep -v '^>' "$(field "$c_ref" "$line")" | tr -d '\n' | wc -c)"
done
