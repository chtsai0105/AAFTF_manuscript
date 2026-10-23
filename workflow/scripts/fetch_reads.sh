#!/usr/bin/env bash
# Download the real sequencing reads listed in config/real_samples.tsv from SRA (sra-tools).
#
#   workflow/scripts/fetch_reads.sh [--force] [tsv]   (default: config/real_samples.tsv)
#
# For each row:
#   sra_short (comma-separated paired-end runs, e.g. SRR1234567) -> read1 / read2 paths
#   sra_long  (comma-separated ONT/PacBio runs)                    -> longreads path
# Several runs are concatenated. Rows whose target files already exist (non-empty) are skipped
# (--force re-downloads them); rows without SRA runs must point to local FASTQ files.
# The same-strain references are downloaded with: workflow/scripts/fetch_samples.sh config/real_samples.tsv
# Rows starting with '#' are commented out and skipped. Run from the repository root. Set THREADS (default 4) for fasterq-dump / pigz.
set -euo pipefail

force=0
if [[ ${1:-} == --force || ${1:-} == -f ]]; then force=1; shift; fi
tsv=${1:-config/real_samples.tsv}
threads=${THREADS:-4}
tmp=$(mktemp -d "${TMPDIR:-/tmp}/fetch_reads.XXXXXX"); trap 'rm -rf "$tmp"' EXIT

col() { head -1 "$tsv" | tr '\t' '\n' | grep -nx "$1" | cut -d: -f1; }
# cut keeps empty fields (bash `read` would merge consecutive tabs and shift the columns)
field() { local c; c=$(col "$1"); [[ -n $c ]] && cut -f"$c" <<< "$2" || true; }

have() { [[ $force -eq 0 && -s $1 ]]; }

# $1 runs (comma-separated), $2 mode (paired|single), $3.. output paths
dump() {
    local runs=$1 mode=$2; shift 2
    local outs=("$@") run
    for run in ${runs//,/ }; do
        echo "      fasterq-dump $run"
        rm -rf "$tmp/$run"
        fasterq-dump --split-files --skip-technical -e "$threads" -t "$tmp" -O "$tmp/$run" "$run"
    done
    for o in "${outs[@]}"; do mkdir -p "$(dirname "$o")"; : > "$o.tmp"; done
    for run in ${runs//,/ }; do
        if [[ $mode == paired ]]; then
            [[ -s $tmp/$run/${run}_1.fastq && -s $tmp/$run/${run}_2.fastq ]] \
                || { echo "ERROR: $run is not paired-end (no _1/_2 files)" >&2; exit 1; }
            pigz -p "$threads" -c "$tmp/$run/${run}_1.fastq" >> "${outs[0]}.tmp"
            pigz -p "$threads" -c "$tmp/$run/${run}_2.fastq" >> "${outs[1]}.tmp"
        else
            cat "$tmp/$run"/*.fastq | pigz -p "$threads" -c >> "${outs[0]}.tmp"
        fi
        rm -rf "${tmp:?}/$run"
    done
    for o in "${outs[@]}"; do mv "$o.tmp" "$o"; done
}

tail -n +2 "$tsv" | awk '!/^[[:space:]]*#/' | while IFS= read -r line; do
    name=$(field organism "$line"); [[ -n $name ]] || continue
    r1=$(field read1 "$line") r2=$(field read2 "$line") lr=$(field longreads "$line")
    s_short=$(field sra_short "$line" | tr -d ' ') s_long=$(field sra_long "$line" | tr -d ' ')

    if [[ -n $r1 ]]; then
        if have "$r1" && { [[ -z $r2 ]] || have "$r2"; }; then echo "skip  $name short reads ($r1 exists)"
        elif [[ -n $s_short ]]; then
            [[ -n $r2 ]] || { echo "ERROR: $name: sra_short needs both read1 and read2 paths" >&2; exit 1; }
            echo "fetch $name short: $s_short"; dump "$s_short" paired "$r1" "$r2"
        else echo "ERROR: $name: $r1 missing and no sra_short given" >&2; exit 1; fi
    fi
    if [[ -n $lr ]]; then
        if have "$lr"; then echo "skip  $name long reads ($lr exists)"
        elif [[ -n $s_long ]]; then echo "fetch $name long: $s_long"; dump "$s_long" single "$lr"
        else echo "ERROR: $name: $lr missing and no sra_long given" >&2; exit 1; fi
    fi
done
