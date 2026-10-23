#!/usr/bin/env bash
# Download the fungal reference genomes listed in config/sim_samples.tsv from NCBI.
#
#   workflow/scripts/fetch_samples.sh [--force] [tsv]   (default: config/sim_samples.tsv)
#
# Each row's `accessions` (comma-separated) are downloaded into its `reference` path, in order:
#   GCF_/GCA_ assembly accession  -> whole genome FASTA via the NCBI Datasets API
#   any other accession           -> one sequence via E-utilities (nuccore), e.g. a mitogenome;
#                                    skipped if the assembly already contains it
# An optional ':<copies>' suffix (copies per nuclear genome, used by the read simulators, e.g. for
# the mitochondrion) is ignored here. Rows whose reference already exists as a non-empty FASTA are
# skipped (--force re-downloads them); rows without accessions must point to a local FASTA.
# Rows starting with '#' are commented out and skipped. Run from the repository root. Set NCBI_API_KEY to raise the rate limit.
set -euo pipefail

force=0
if [[ ${1:-} == --force || ${1:-} == -f ]]; then force=1; shift; fi
tsv=${1:-config/sim_samples.tsv}
EUTILS=https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi
DATASETS=https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession
key=${NCBI_API_KEY:+&api_key=$NCBI_API_KEY}
tmp=$(mktemp -d); trap 'rm -rf "$tmp"' EXIT

col() { head -1 "$tsv" | tr '\t' '\n' | grep -nx "$1" | cut -d: -f1; }
c_name=$(col organism) c_ref=$(col reference) c_acc=$(col accessions)
[[ -n $c_name && -n $c_ref && -n $c_acc ]] || { echo "ERROR: $tsv needs columns organism, reference, accessions" >&2; exit 1; }

fetch_assembly() {   # $1 accession, $2 output fasta (appended)
    curl -sSf --retry 5 --retry-delay 3 -o "$tmp/asm.zip" \
        "$DATASETS/$1/download?include_annotation_type=GENOME_FASTA${NCBI_API_KEY:+&api_key=$NCBI_API_KEY}"
    rm -rf "$tmp/asm"; unzip -q -o "$tmp/asm.zip" -d "$tmp/asm"
    local fna; fna=$(find "$tmp/asm" -name '*.fna' | sort)
    [[ -n $fna ]] || { echo "ERROR: $1: no genome FASTA in the Datasets download" >&2; return 1; }
    cat $fna >> "$2"
}

fetch_nuccore() {    # $1 accession, $2 output fasta (appended)
    curl -sSf --retry 5 --retry-delay 3 "$EUTILS?db=nuccore&id=$1&rettype=fasta&retmode=text$key" > "$tmp/seq.fa"
    [[ $(grep -c '^>' "$tmp/seq.fa") -eq 1 ]] || { echo "ERROR: $1: expected one sequence" >&2; return 1; }
    cat "$tmp/seq.fa" >> "$2"
}

# cut keeps empty fields (bash `read` would merge consecutive tabs and shift the columns)
field() { cut -f"$1" <<< "$2"; }

tail -n +2 "$tsv" | awk '!/^[[:space:]]*#/' | while IFS= read -r line; do
    name=$(field "$c_name" "$line") fa=$(field "$c_ref" "$line")
    accs=$(field "$c_acc" "$line" | sed -E 's/:[0-9.]+//g; s/ //g')
    [[ -n $name ]] || continue
    if [[ $force -eq 0 && -s $fa && $(head -c1 "$fa") == ">" ]]; then
        echo "skip  $name ($fa exists)"; continue
    fi
    if [[ -z $accs ]]; then
        [[ -s $fa ]] && { echo "keep  $name (local file, no accessions)"; continue; }
        echo "ERROR: $name: $fa missing and no accessions given" >&2; exit 1
    fi
    echo "fetch $name: $accs"
    mkdir -p "$(dirname "$fa")"; : > "$fa.tmp"
    for acc in ${accs//,/ }; do
        if [[ $acc =~ ^GC[AF]_[0-9]+\.[0-9]+$ ]]; then
            fetch_assembly "$acc" "$fa.tmp"
        elif grep -q "^>$acc[[:space:]]" "$fa.tmp"; then
            echo "      $acc already in the assembly"; continue
        else
            fetch_nuccore "$acc" "$fa.tmp"
        fi
        sleep 0.4    # stay under NCBI's request rate without an API key
    done
    mv "$fa.tmp" "$fa"
done

echo "sample references (compare with the genome_size column):"
tail -n +2 "$tsv" | awk '!/^[[:space:]]*#/' | while IFS= read -r line; do
    fa=$(field "$c_ref" "$line")
    printf "  %-32s %5d sequences  %12s bp\n" "$(field "$c_name" "$line")" "$(grep -c '^>' "$fa")" "$(grep -v '^>' "$fa" | tr -d '\n' | wc -c)"
done
