#!/usr/bin/env bash
# Build a writable compleasm library (-L) whose lineage folders are symlinks into a BUSCO download
# folder, so BUSCO and compleasm use the same lineage files without a second download.
#
#   workflow/scripts/setup_compleasm_library.sh [busco_downloads] [compleasm_library]
#   (defaults: busco_downloads / compleasm_library in config/params.yaml)
#
# Lineages are taken from the `lineage` column of config/samples.tsv and config/real_samples.tsv,
# plus eukaryota_<odb>, which compleasm always loads. For each lineage L: <library>/L ->
# <busco_downloads>/lineages/L and the L.done marker compleasm checks. placement_files/ (+ .done) is
# created empty: compleasm only uses it for --autolineage, and otherwise downloads it on first run.
# file_versions.tsv is copied so compleasm does not re-download it. Existing entries are kept, so
# parallel pipeline jobs never start a download. Run from the repository root.
set -euo pipefail

param() { awk -v k="$1" '$1 == k":" {print $2}' config/params.yaml; }
src=${1:-$(param busco_downloads)}
lib=${2:-$(param compleasm_library)}
[[ -n $src && $src != null && -d $src/lineages ]] || { echo "ERROR: BUSCO download folder '$src' has no lineages/ folder" >&2; exit 1; }
[[ -n $lib && $lib != null ]] || { echo "ERROR: set compleasm_library in config/params.yaml" >&2; exit 1; }
mkdir -p "$lib"

lineages=$(for t in config/samples.tsv config/real_samples.tsv; do
    [[ -f $t ]] || continue
    c=$(head -1 "$t" | tr '\t' '\n' | grep -nx lineage | cut -d: -f1)
    [[ -n $c ]] && tail -n +2 "$t" | cut -f"$c"
done | grep . | sort -u)
odbs=$(sed -E 's/.*_(odb[0-9.]+)$/\1/' <<< "$lineages" | sort -u)
lineages=$(printf "%s\n" $lineages $(sed 's/^/eukaryota_/' <<< "$odbs") | sort -u)

for l in $lineages; do
    [[ -d $src/lineages/$l ]] || { echo "ERROR: $l not in $src/lineages" >&2; exit 1; }
    if [[ -e $lib/$l ]]; then echo "keep   $lib/$l"
    else ln -s "$src/lineages/$l" "$lib/$l"; echo "linked $lib/$l -> $src/lineages/$l"; fi
    touch "$lib/$l.done"
done
mkdir -p "$lib/placement_files" && touch "$lib/placement_files.done"
[[ -f $src/file_versions.tsv ]] && cp -f "$src/file_versions.tsv" "$lib/file_versions.tsv"
echo "compleasm library ready: $lib"
