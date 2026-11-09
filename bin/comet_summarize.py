#!/usr/bin/env python3
"""Merge chunk results: per-read table, read lists, composition and genome-wide methylation.

Inputs
  prob     concatenated assignment output (read_id, n_cpg, label, P(1)..P(K), lineage)
  gw       concatenated per-read call counts (read_id, np, n_meth, n_unmeth)
Outputs (under --out-dir, prefixed by --sample)
  01.cell_type_assignment/cell_type_prob/<s>.cell_type_prob.tsv  per-read posteriors, cell type, lineage
  01.cell_type_assignment/read_list/cell_type/<s>.<type>.reads.txt
  01.cell_type_assignment/read_list/brain_lineage_rescued/<s>.<lineage>.reads.txt
  02.cell_type_quantification/<s>.composition.tsv  reads per cell type / lineage, three proportions:
                                    frac_single         among single-cell-type assigned reads
                                    frac_single_lineage among single-type + lineage-rescued reads
                                    frac_all            among all scored reads (incl. unassigned)
  02.cell_type_quantification/<s>.gw_methylation.tsv  genome-wide CpG methylation per read set
"""
import os
import argparse
import pandas as pd


def frac(n, d):
    return '%.6f' % (n / d) if d else 'NA'


def level(m, u):
    return '%.6f' % (m / (m + u)) if m + u > 0 else 'NA'


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('prob')
    p.add_argument('gw')
    p.add_argument('--out-dir', required=True)
    p.add_argument('--sample', required=True)
    p.add_argument('--names', required=True, help='comma-separated cell-type names, atlas order')
    p.add_argument('--lineages', default='', help='comma-separated lineage names, if rescue was on')
    a = p.parse_args()

    names = a.names.split(',')
    lineages = [x for x in a.lineages.split(',') if x]
    k = len(names)
    adir = os.path.join(a.out_dir, '01.cell_type_assignment')
    qdir = os.path.join(a.out_dir, '02.cell_type_quantification')
    os.makedirs(adir, exist_ok=True)
    os.makedirs(qdir, exist_ok=True)

    cols = ['read_id', 'n_cpg', 'label'] + ['P_' + n for n in names] + ['rescued_lineage']
    if os.path.getsize(a.prob) > 0:
        df = pd.read_csv(a.prob, sep='\t', header=None, names=cols, dtype={'read_id': str},
                         keep_default_na=False)
    else:
        df = pd.DataFrame(columns=cols)
    df.insert(3, 'cell_type', [names[l - 1] if l > 0 else 'Unassigned' for l in df['label']])
    pdir = os.path.join(adir, 'cell_type_prob')
    os.makedirs(pdir, exist_ok=True)
    df.to_csv(os.path.join(pdir, a.sample + '.cell_type_prob.tsv'), sep='\t', index=False)

    groups = {n: df.loc[df['label'] == i, 'read_id'].tolist() for i, n in enumerate(names, 1)}
    for g in lineages:
        groups[g] = df.loc[df['rescued_lineage'] == g, 'read_id'].tolist()
    unassigned = df.loc[(df['label'] == 0) & (df['rescued_lineage'] == '.'), 'read_id'].tolist()

    for g, reads in groups.items():
        sub = 'brain_lineage_rescued' if g in lineages else 'cell_type'
        rdir = os.path.join(adir, 'read_list', sub)
        os.makedirs(rdir, exist_ok=True)
        with open(os.path.join(rdir, '%s.%s.reads.txt' % (a.sample, g)), 'w') as f:
            f.writelines(r + '\n' for r in reads)

    n_single = sum(len(groups[n]) for n in names)
    n_lin = sum(len(groups[g]) for g in lineages)
    n_all = len(df)
    with open(os.path.join(qdir, a.sample + '.composition.tsv'), 'w') as f:
        f.write('group\tn_reads\tfrac_single\tfrac_single_lineage\tfrac_all\n')
        for n in names:
            c = len(groups[n])
            f.write('%s\t%d\t%s\t%s\t%s\n' % (n, c, frac(c, n_single), frac(c, n_single + n_lin),
                                              frac(c, n_all)))
        for g in lineages:
            c = len(groups[g])
            f.write('%s\t%d\tNA\t%s\t%s\n' % (g, c, frac(c, n_single + n_lin), frac(c, n_all)))
        f.write('Unassigned\t%d\tNA\tNA\t%s\n' % (len(unassigned), frac(len(unassigned), n_all)))

    # genome-wide methylation: Bulk = every read with calls, then each read set
    read_group = {}
    for g, reads in groups.items():
        for r in reads:
            read_group[r] = g
    for r in unassigned:
        read_group[r] = 'Unassigned'
    order = ['Bulk'] + list(groups) + ['Unassigned']
    tot = {g: [0, 0, 0] for g in order}           # reads, meth, unmeth
    with open(a.gw) as f:
        for line in f:
            rid, _npc, m, u = line.rstrip('\n').split('\t')
            m, u = int(m), int(u)
            gs = ['Bulk', read_group[rid]] if rid in read_group else ['Bulk']
            for g in gs:
                t = tot[g]
                t[0] += 1
                t[1] += m
                t[2] += u
    with open(os.path.join(qdir, a.sample + '.gw_methylation.tsv'), 'w') as f:
        f.write('group\tn_reads\tn_meth_calls\tn_unmeth_calls\tmethylation\n')
        for g in order:
            r, m, u = tot[g]
            f.write('%s\t%d\t%d\t%d\t%s\n' % (g, r, m, u, level(m, u)))
    print('Scored reads: %d | single-type assigned: %d | lineage rescued: %d | unassigned: %d' % (
        n_all, n_single, n_lin, len(unassigned)))


if __name__ == '__main__':
    main()
