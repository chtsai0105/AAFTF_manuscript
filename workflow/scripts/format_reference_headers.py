#!/usr/bin/env python3
"""Rename every FASTA record to {prefix}_{seqid}, e.g. FUNGUS_chr1 or CONTAM_Ralstonia_NC_010682.1.

Optionally tags records with a copy number for the read simulators, taken from the `accessions`
column of contaminants.tsv (seqid[:copies], comma-separated; default 1 copy per nuclear genome):
  --accessions NC_003075.7,NC_000932.1:300,NC_037304.1:20  ->  depth=300 / depth=20 in the headers
Badread reads depth= directly; make_readcounts.py uses it for InSilicoSeq read counts. Only the first
word of a header is the ID, so QUAST, minimap2 and the decontamination scoring ignore the tag.
"""
import argparse, gzip, re

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument('fasta')
ap.add_argument('--prefix', required=True, help='FUNGUS or CONTAM_<TAXON>')
ap.add_argument('--accessions', default='', help='comma list of seqid[:copies]')
a = ap.parse_args()

copies = {}
for x in filter(None, (x.strip() for x in a.accessions.split(','))):
    acc, _, n = x.partition(':')
    copies[acc] = float(n) if n else 1.0
seen = set()

opener = gzip.open if a.fasta.endswith('.gz') else open
with opener(a.fasta, 'rt') as fh:
    for line in fh:
        if line.startswith('>'):
            orig = line[1:].split()[0]
            seqid = re.sub(r'[^A-Za-z0-9_.-]', '_', orig)
            seen.add(orig)
            tag = f' depth={copies[orig]:g}' if copies.get(orig, 1.0) != 1.0 else ''
            print(f'>{a.prefix}_{seqid}{tag}')
        else:
            print(line.rstrip())

missing = {k for k, v in copies.items() if v != 1.0} - seen
if missing:
    raise SystemExit(f'ERROR: {a.fasta}: copy numbers given for sequences not in the file: {", ".join(sorted(missing))}')
