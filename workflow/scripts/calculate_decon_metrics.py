#!/usr/bin/env python3
"""Contig taxonomy and decontamination confusion matrix from minimap2 PAF files.

Each contig is assigned the class (Fungal / Contaminant) with the most aligned query bases
in its PAF hits (targets ^FUNGUS_ / ^CONTAM_); contigs with no hit are Unassigned.

Always reports the composition (bp per class) of --fasta.
With --before_fasta/--before_paf, also reports the confusion matrix of contigs removed between
the two assemblies (matched by contig name): TP = contaminant bp removed, FN = contaminant bp kept,
FP = fungal bp removed, TN = fungal bp kept. Unassigned contigs are excluded from the matrix.
"""
import argparse, gzip
import pandas as pd

PAF_COLS = ['qname', 'qlen', 'qstart', 'qend', 'strand', 'tname']

def fasta_lengths(fn):
    d, name = {}, None
    with (gzip.open if fn.endswith('.gz') else open)(fn, 'rt') as fh:
        for line in fh:
            if line[0] == '>':
                name = line[1:].split()[0]; d[name] = 0
            else:
                d[name] += len(line.strip())
    return pd.DataFrame({'contig': list(d), 'length': list(d.values())})

def classify(fasta, paf):
    contigs = fasta_lengths(fasta)
    try:
        hits = pd.read_csv(paf, sep='\t', header=None, usecols=range(6), names=PAF_COLS)
    except pd.errors.EmptyDataError:
        hits = pd.DataFrame(columns=PAF_COLS)
    hits['cls'] = hits.tname.str.extract(r'^(FUNGUS|CONTAM)_', expand=False).map(
        {'FUNGUS': 'Fungal', 'CONTAM': 'Contaminant'})
    hits['aln'] = hits.qend - hits.qstart
    best = (hits.dropna(subset=['cls']).groupby(['qname', 'cls']).aln.sum().reset_index()
                .sort_values('aln', ascending=False).drop_duplicates('qname'))
    contigs = contigs.merge(best[['qname', 'cls']], left_on='contig', right_on='qname', how='left')
    contigs['cls'] = contigs.cls.fillna('Unassigned')
    return contigs[['contig', 'length', 'cls']]

ap = argparse.ArgumentParser(description=__doc__)
ap.add_argument('--dataset', required=True)
ap.add_argument('--checkpoint', required=True)
ap.add_argument('--fasta', required=True)
ap.add_argument('--paf', required=True)
ap.add_argument('--before_fasta')
ap.add_argument('--before_paf')
ap.add_argument('-o', '--out', default='precision_recall_f1.tsv')
a = ap.parse_args()

after = classify(a.fasta, a.paf)
bp = after.groupby('cls').length.sum()
row = {'dataset': a.dataset, 'checkpoint': a.checkpoint,
       'fungal_bp': int(bp.get('Fungal', 0)), 'contam_bp': int(bp.get('Contaminant', 0)),
       'unassigned_bp': int(bp.get('Unassigned', 0)),
       'contam_contigs': int((after.cls == 'Contaminant').sum())}

if a.before_fasta and a.before_paf:
    before = classify(a.before_fasta, a.before_paf)
    before['removed'] = ~before.contig.isin(after.contig)
    sel = lambda cls, removed: int(before.loc[(before.cls == cls) & (before.removed == removed), 'length'].sum())
    tp, fn = sel('Contaminant', True), sel('Contaminant', False)
    fp, tn = sel('Fungal', True), sel('Fungal', False)
    prec = tp / (tp + fp) if tp + fp else float('nan')
    rec = tp / (tp + fn) if tp + fn else float('nan')
    f1 = 2 * prec * rec / (prec + rec) if prec + rec > 0 else float('nan')
    row.update(TP_bp=tp, FN_bp=fn, FP_bp=fp, TN_bp=tn, precision=prec, recall=rec, F1=f1)

pd.DataFrame([row]).to_csv(a.out, sep='\t', index=False)
