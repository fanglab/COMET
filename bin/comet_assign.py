#!/usr/bin/env python3
"""Read-level cell-type assignment (Methods 4.1-4.4).

Input : calls matched to the atlas (bedtools intersect -wo), columns
        0-6  chrom, pos, pos+1, ML, np, strand, read_id
        7-9  atlas chrom, start, end
        10.. K atlas beta values, then the overlap length
Output: read_id, n_cpg, label (0 = unassigned), P(1) .. P(K), lineage ('.' if none)

Per read: binarise calls (ML >= tau -> 1, ML <= 255 - tau -> 0, else dropped),
Bernoulli log-likelihood per cell type with beta clipped to [1e-10, 1 - 1e-10],
uniform prior, softmax posterior. Assign the top type if P1 >= theta and P1 >= 2 * P2.

Lineage rescue (Methods 4.5, --lineage NAME=i,j,..): an unassigned read goes to lineage G
when, for G or one of its subsets S (|S| >= 2, largest first), sum(P[S]) >= --lineage-thres
and min(P[S]) > max(P[not S]) + --lineage-margin.
"""
import argparse
from itertools import combinations

import numpy as np
import pandas as pd
from scipy.special import softmax

EPS = 1e-10


def log_likelihood_binary(x, p):
    p = np.clip(p, EPS, 1 - EPS)
    return np.sum(x * np.log(p) + (1 - x) * np.log(1 - p))


def classify(x, betas):
    k = len(betas)
    log_scores = np.array([log_likelihood_binary(x, b) + np.log(1.0 / k) for b in betas])
    probs = softmax(log_scores)
    return int(np.argmax(log_scores)) + 1, probs


def check_assign(label, probs, theta):
    s = np.sort(probs)[::-1]
    if s[0] >= theta and s[0] >= 2 * s[1]:
        return label
    return 0


def parse_lineages(specs, k):
    out = []
    for spec in specs:
        name, members = spec.split('=', 1)
        idx = sorted({int(m) for m in members.split(',')})
        if len(idx) < 2 or idx[0] < 1 or idx[-1] > k:
            raise SystemExit('ERROR: --lineage %s needs >= 2 cell types in 1..%d' % (spec, k))
        out.append((name, idx))
    return out


def rescue(probs, lineages, thres, margin):
    k = len(probs)
    for name, members in lineages:
        for size in range(len(members), 1, -1):
            for sub in combinations(members, size):
                inside = [probs[i - 1] for i in sub]
                outside = [probs[i] for i in range(k) if i + 1 not in sub]
                if sum(inside) >= thres and min(inside) > max(outside) + margin:
                    return name
    return '.'


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('matched')
    p.add_argument('out')
    p.add_argument('--n-ct', type=int, required=True)
    p.add_argument('--cutoff', type=int, default=180)
    p.add_argument('--assign-thres', type=float, required=True)
    p.add_argument('--min-cpg', type=int, default=5)
    p.add_argument('--pn-table', default=None,
                   help='pass-number table (PacBio); omit for a fixed cutoff')
    p.add_argument('--lineage', action='append', default=[])
    p.add_argument('--lineage-thres', type=float, default=0.65)
    p.add_argument('--lineage-margin', type=float, default=0.10)
    a = p.parse_args()

    k = a.n_ct
    lineages = parse_lineages(a.lineage, k)
    beta_cols = list(range(10, 10 + k))
    try:
        df = pd.read_csv(a.matched, sep='\t', header=None, usecols=[3, 4, 6] + beta_cols,
                         dtype={6: str})
    except pd.errors.EmptyDataError:
        open(a.out, 'w').close()
        return

    ml = df[3].to_numpy()
    if a.pn_table:
        pn_df = pd.read_csv(a.pn_table, sep='\t', header=None, index_col=0)[[4]]
        tau = df[4].map(pn_df.to_dict()[4]).fillna(256).to_numpy()
    else:
        tau = np.full(len(df), a.cutoff)
    state = np.where(ml >= tau, 1.0, np.where(ml <= 255 - tau, 0.0, np.nan))
    keep = ~np.isnan(state)
    df = df.loc[keep]
    state = state[keep]

    reads = df[6].to_numpy()
    betas = [np.ascontiguousarray(df[c].to_numpy(dtype=np.float64)) for c in beta_cols]

    # group rows by read, keeping each read's original row order
    order = np.argsort(reads, kind='stable')
    reads_sorted = reads[order]
    bounds = np.flatnonzero(reads_sorted[1:] != reads_sorted[:-1]) + 1
    starts = np.concatenate(([0], bounds))
    ends = np.concatenate((bounds, [len(order)]))

    with open(a.out, 'w') as out:
        for s, e in zip(starts, ends):
            n = e - s
            if n < a.min_cpg:
                continue
            idx = order[s:e]
            x = state[idx]
            label, probs = classify(x, [b[idx] for b in betas])
            label = check_assign(label, probs, a.assign_thres)
            lin = rescue(probs, lineages, a.lineage_thres, a.lineage_margin) \
                if label == 0 and lineages else '.'
            out.write(reads_sorted[s] + '\t' + str(n) + '\t' + str(label) + '\t'
                      + '\t'.join(str(v) for v in probs) + '\t' + lin + '\n')


if __name__ == '__main__':
    main()
