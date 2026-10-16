#!/usr/bin/env python3
"""Write an InSilicoSeq --readcount_file with exact read counts per record.

Fungal reads: genome_size * coverage / read_length, spread over records by length.
Contaminant reads: N_fungus * p / (100 - p), split by group weight, then equally among taxa
of the group, then by record length x copy number (depth= in the header, default 1, e.g. an
organelle at depth=300 gets 300x the coverage of that organism's nuclear chromosomes).
Counts are rounded to even numbers (read pairs).
"""
import argparse, collections, gzip, os, re

ap = argparse.ArgumentParser(description=__doc__)
ap.add_argument('--fungus', required=True)
ap.add_argument('--contam', nargs='*', default=[], help='tagged contaminant FASTAs named CONTAM_<taxon>.fa')
ap.add_argument('--groups', nargs='*', default=[], help='taxon:group for each contaminant')
ap.add_argument('--weights', default='bacteria:0.65,host:0.15,vector:0.1,fungal:0.1')
ap.add_argument('--coverage', type=float, required=True)
ap.add_argument('--read_length', type=int, default=150)
ap.add_argument('--percent', type=float, required=True)
ap.add_argument('--out_fungus', required=True)
ap.add_argument('--out_contam', required=True)
a = ap.parse_args()

def lengths(fn):
    """Return {record: length x depth}; depth= from the header (default 1)."""
    d, depth, name = {}, {}, None
    with (gzip.open if fn.endswith('.gz') else open)(fn, 'rt') as fh:
        for l in fh:
            if l[0] == '>':
                name = l[1:].split()[0]; d[name] = 0
                m = re.search(r'depth=([\d.]+)', l)
                depth[name] = float(m.group(1)) if m else 1.0
            else:
                d[name] += len(l.strip())
    return {k: v * depth[k] for k, v in d.items()}

even = lambda x: int(round(x / 2)) * 2

def write(fn, lens, total):
    size = sum(lens.values())
    with open(fn, 'w') as out:
        for rec, L in lens.items():
            n = even(total * L / size)
            if n: out.write(f'{rec}\t{n}\n')

fl = lengths(a.fungus)
n_fungus = sum(fl.values()) * a.coverage / a.read_length
write(a.out_fungus, fl, n_fungus)

open(a.out_contam, 'w').close()
if a.percent > 0 and a.contam:
    n_contam = n_fungus * a.percent / (100 - a.percent)
    weights = {k: float(v) for k, v in (x.split(':') for x in a.weights.split(','))}
    group_of = dict(x.split(':') for x in a.groups)
    per_group = collections.Counter(group_of.values())
    # renormalise in case a group has no taxa
    wsum = sum(weights[g] for g in per_group)
    with open(a.out_contam, 'w') as out:
        for fn in a.contam:
            lens = lengths(fn)
            taxon = os.path.basename(fn).removeprefix('CONTAM_').rsplit('.', 1)[0]
            g = group_of[taxon]
            taxon_reads = n_contam * weights[g] / wsum / per_group[g]
            size = sum(lens.values())
            for rec, L in lens.items():
                n = even(taxon_reads * L / size)
                if n: out.write(f'{rec}\t{n}\n')
