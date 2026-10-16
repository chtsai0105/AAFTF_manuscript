#!/usr/bin/env python3
"""Combine QUAST, compleasm, BUSCO, decontamination and timing outputs into a tidy final_report.tsv.

Output columns: dataset, sample, spike_pct, track, checkpoint, tool, metric, value
Evaluation directories are named <dataset>__<checkpoint>; datasets are named <sample>_spike<pct>_<track> (track: short, long, hybrid).
"""
import argparse, glob, os, re
import pandas as pd

QUAST_KEEP = ['Genome fraction (%)', 'N50', 'NGA50', 'auN', 'auNGA', '# contigs', '# contigs (>= 0 bp)',
              '# contigs (>= 1000 bp)', 'Total length', 'Largest contig', 'GC (%)',
              '# misassemblies', 'Unaligned length']
MISASM_KEEP = ['# c. relocations', '# c. translocations', '# c. inversions']

ap = argparse.ArgumentParser(description=__doc__)
ap.add_argument('--quast', nargs='*', default=[], help='QUAST output dirs')
ap.add_argument('--compleasm', nargs='*', default=[], help='compleasm output dirs')
ap.add_argument('--busco', nargs='*', default=[], help='BUSCO output dirs')
ap.add_argument('--decon', nargs='*', default=[], help='calculate_decon_metrics.py TSVs')
ap.add_argument('--times', nargs='*', default=[], help='timeit.sh TSVs')
ap.add_argument('-o', '--out', default='final_report.tsv')
a = ap.parse_args()

rows = []
def add(ds, cp, tool, metric, value):
    rows.append(dict(dataset=ds, checkpoint=cp, tool=tool, metric=metric, value=value))

def ds_cp(path):
    return os.path.basename(os.path.normpath(path)).split('__', 1)

def transposed(fn, keep):
    if not os.path.exists(fn): return {}
    df = pd.read_csv(fn, sep='\t')
    return {k: df[k].iloc[0] for k in keep if k in df.columns}

for d in a.quast:
    ds, cp = ds_cp(d)
    vals = transposed(f'{d}/transposed_report.tsv', QUAST_KEEP)
    vals.update(transposed(f'{d}/contigs_reports/transposed_report_misassemblies.tsv', MISASM_KEEP))
    for k, v in vals.items(): add(ds, cp, 'quast', k, v)

for d in a.compleasm:
    ds, cp = ds_cp(d)
    for fn in glob.glob(f'{d}/summary.txt'):
        for key, pct in re.findall(r'^([SDFIM]):([\d.]+)%', open(fn).read(), re.M):
            add(ds, cp, 'compleasm', f'{key}%', float(pct))

for d in a.busco:
    ds, cp = ds_cp(d)
    for fn in glob.glob(f'{d}/**/short_summary*.txt', recursive=True)[:1]:
        m = re.search(r'C:([\d.]+)%\[S:([\d.]+)%,D:([\d.]+)%\],F:([\d.]+)%,M:([\d.]+)%', open(fn).read())
        if m:
            for key, v in zip('CSDFM', m.groups()): add(ds, cp, 'busco', f'{key}%', float(v))

for fn in a.decon:
    df = pd.read_csv(fn, sep='\t')
    for _, r in df.iterrows():
        for k in df.columns.drop(['dataset', 'checkpoint']):
            add(r.dataset, r.checkpoint, 'decon', k, r[k])

for fn in a.times:
    df = pd.read_csv(fn, sep='\t')
    for _, r in df.iterrows():
        for k in ['wall_s', 'max_rss_mb', 'cpu_pct']:
            add(r.dataset, r.checkpoint, 'time', f'{r.step}:{k}', r[k])

out = pd.DataFrame(rows, columns=['dataset', 'checkpoint', 'tool', 'metric', 'value'])
parts = out.dataset.str.extract(r'^(?P<sample>.+)_spike(?P<spike_pct>[\d.]+)_(?P<track>short|long|hybrid)$')
out.insert(1, 'sample', parts['sample'])
out.insert(2, 'spike_pct', parts['spike_pct'])
out.insert(3, 'track', parts['track'])
out.sort_values(['dataset', 'checkpoint', 'tool', 'metric']).to_csv(a.out, sep='\t', index=False)
