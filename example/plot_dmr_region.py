#!/usr/bin/env python3
"""Methylation tracks of one region from per-site BED files.

usage: plot_dmr_region.py --region chr:start-end --track LABEL=FILE[,FILE...] [--track ...] --out PREFIX
                          [--min-depth 2]

One track per --track (top to bottom): the per-CpG methylation profile (n_mod / n_total of CpGs
with >= --min-depth reads), drawn as a line with the area below filled. Several comma-separated
FILEs in one track are pooled (n_mod and n_total summed per CpG).
FILE: bgzipped, tabix-indexed per-site BED (chrom, start, end, n_mod, n_total), as written by
comet assign. Writes PREFIX.pdf and PREFIX.png.
"""
import os
import sys
import glob
import argparse

import pysam
import matplotlib
matplotlib.use('Agg')
from matplotlib import font_manager
import matplotlib.pyplot as plt

COLOR = '#4C72B0'   # matplotlib / seaborn default blue


def use_arial():
    for d in (os.environ.get('COMET_FONT_DIR', ''), os.path.join(sys.prefix, 'fonts'),
              os.path.join(os.environ.get('CONDA_PREFIX', ''), 'fonts')):
        if d and os.path.isdir(d):
            for f in glob.glob(os.path.join(d, '[Aa]rial*.ttf')):
                font_manager.fontManager.addfont(f)
    if 'Arial' in {f.name for f in font_manager.fontManager.ttflist}:
        matplotlib.rcParams['font.family'] = 'Arial'
    else:
        print('WARNING: Arial not found (install mscorefonts or set COMET_FONT_DIR)', file=sys.stderr)
    matplotlib.rcParams.update({'font.size': 7, 'pdf.fonttype': 42, 'ps.fonttype': 42,
                                'savefig.dpi': 300, 'savefig.bbox': 'tight'})


def profile(paths, chrom, start, end, min_depth):
    counts = {}
    for path in paths:
        with pysam.TabixFile(path) as tb:
            for line in tb.fetch(chrom, start, end):
                r = line.split('\t')
                c = counts.setdefault(int(r[1]), [0, 0])
                c[0] += int(r[3])
                c[1] += int(r[4])
    x, y = [], []
    for pos in sorted(counts):
        n_mod, n_tot = counts[pos]
        if n_tot >= min_depth:
            x.append(pos)
            y.append(n_mod / n_tot)
    return x, y


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--region', required=True)
    p.add_argument('--track', action='append', required=True, help='LABEL=FILE[,FILE...]')
    p.add_argument('--out', required=True)
    p.add_argument('--min-depth', type=int, default=2)
    a = p.parse_args()

    chrom, rng = a.region.replace(',', '').split(':')
    start, end = (int(v) for v in rng.split('-'))
    use_arial()
    n = len(a.track)
    fig, axes = plt.subplots(n, 1, figsize=(3.6, 0.42 * n + 0.4), sharex=True, squeeze=False,
                             gridspec_kw={'hspace': 0.25})
    axes = axes[:, 0]
    for ax, t in zip(axes, a.track):
        label, path = t.split('=', 1)
        x, y = profile(path.split(','), chrom, start, end, a.min_depth)
        ax.plot(x, y, color=COLOR, lw=0.8)
        ax.fill_between(x, y, 0, color=COLOR, alpha=0.6, lw=0)
        ax.set_ylim(0, 1.2)
        ax.set_yticks([0, 1])
        ax.set_ylabel(label, rotation=0, ha='right', va='center', fontsize=7)
        for side in ('top', 'right'):
            ax.spines[side].set_visible(False)
        ax.tick_params(axis='x', length=0)
        ax.tick_params(axis='y', length=2, labelsize=6)
    last = axes[-1]
    last.set_xlim(start, end)
    ticks = [start, (start + end) // 2, end]
    last.set_xticks(ticks)
    last.set_xticklabels(['{:,}'.format(v) for v in ticks], fontsize=6)
    last.tick_params(axis='x', length=2)
    last.set_xlabel(chrom, fontsize=7)
    for ext in ('pdf', 'png'):
        fig.savefig('%s.%s' % (a.out, ext))
    print('%s.pdf / .png: %d tracks, %s:%d-%d' % (a.out, n, chrom, start, end))


if __name__ == '__main__':
    main()
