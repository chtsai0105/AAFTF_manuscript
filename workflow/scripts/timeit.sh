#!/usr/bin/env bash
# timeit.sh <dataset> <checkpoint> <step> <command...>
# Runs the command under GNU time and appends one row to <step>.time.tsv:
# dataset checkpoint step wall_s max_rss_mb cpu_pct
set -euo pipefail
ds=$1 cp=$2 step=$3; shift 3
rc=0
/usr/bin/time -v -o "$step.time.raw" "$@" || rc=$?
awk -v ds="$ds" -v cp="$cp" -v st="$step" -F': ' '
  /Elapsed \(wall clock\)/ { n=split($2,t,":"); w=(n==3)? t[1]*3600+t[2]*60+t[3] : t[1]*60+t[2] }
  /Maximum resident set size/ { rss=$2/1024 }
  /Percent of CPU/ { cpu=$2; sub("%","",cpu) }
  /Exit status/ { ex=$2 }
  END { printf "dataset\tcheckpoint\tstep\twall_s\tmax_rss_mb\tcpu_pct\texit_status\n%s\t%s\t%s\t%.1f\t%.1f\t%s\t%s\n", ds, cp, st, w, rss, cpu, ex }
' "$step.time.raw" > "$step.time.tsv"
exit $rc
