#!/usr/bin/env python3
"""Badread --quantity (bases) per contaminant so contaminants make up --percent of all long-read bases.

Contaminant bases = fungal_bases * p / (100 - p), split by group weight, then equally among the
taxa of a group. Prints: <contaminant fasta>\t<bases>
"""
import argparse, collections, gzip, os

ap = argparse.ArgumentParser(description=__doc__)
ap.add_argument('--fungus', required=True)
ap.add_argument('--coverage', type=float, required=True)
ap.add_argument('--percent', type=float, required=True)
ap.add_argument('--contam', nargs='*', default=[])
ap.add_argument('--groups', nargs='*', default=[], help='taxon:group for each contaminant')
ap.add_argument('--weights', default='bacteria:0.65,host:0.15,vector:0.1,fungal:0.1')
a = ap.parse_args()

def genome_size(fn):
    with (gzip.open if fn.endswith('.gz') else open)(fn, 'rt') as fh:
        return sum(len(l.strip()) for l in fh if l[0] != '>')

def taxon_of(fn):
    """Contaminant FASTAs are named CONTAM_<taxon>.fa"""
    return os.path.basename(fn).removeprefix('CONTAM_').rsplit('.', 1)[0]

if a.percent > 0:
    fungal_bases = genome_size(a.fungus) * a.coverage
    contam_bases = fungal_bases * a.percent / (100 - a.percent)
    weights = {k: float(v) for k, v in (x.split(':') for x in a.weights.split(','))}
    group_of = dict(x.split(':') for x in a.groups)
    per_group = collections.Counter(group_of.values())
    wsum = sum(weights[g] for g in per_group)
    for fn in a.contam:
        g = group_of[taxon_of(fn)]
        print(f'{fn}\t{int(contam_bases * weights[g] / wsum / per_group[g])}')
