#!/usr/bin/env python3
"""Plots (and combined tables) of comet assign results for one or more samples.

usage: comet_plot_assign.py --sample NAME=QUANT_DIR [--sample ...] --out-prefix PREFIX [--cell-types A,B,C]
                            [--meth-samples S1,S2] [--tables] [--no-plot]

QUANT_DIR is a sample's 02.cell_type_quantification folder (NAME.composition.tsv,
NAME.gw_methylation.tsv). Cell types default to the single-type groups of the first
composition table. Writes
  PREFIX.composition.{pdf,png}            stacked horizontal bars, one row per sample:
                                          proportions of the reads assigned to a single cell type
  PREFIX.cell_type_methylation.{pdf,png}  dot plot, one row per sample: genome-wide CpG methylation
                                          of each cell type's reads (--meth-samples: rows to show)
  with --tables: PREFIX.composition.tsv and PREFIX.gw_methylation.tsv (all samples, long format)
Samples are drawn top to bottom in the order given. Text is set in Arial when it is installed
(conda package mscorefonts, or a folder of .ttf files in $COMET_FONT_DIR).
"""
import os
import sys
import glob
import argparse
import warnings

# one colour per cell type, in panel column order (colour-blind-safe, Okabe-Ito based)
PALETTE = ['#0072B2', '#E69F00', '#009E73', '#CC79A7', '#56B4E9', '#D55E00', '#F0E442', '#999999',
           '#000000', '#882255', '#44AA99', '#117733']


def setup_matplotlib():
    import matplotlib
    matplotlib.use('Agg')
    from matplotlib import font_manager
    dirs = [os.environ.get('COMET_FONT_DIR', ''), os.path.join(sys.prefix, 'fonts'),
            os.path.join(os.environ.get('CONDA_PREFIX', ''), 'fonts'), os.path.expanduser('~/.fonts')]
    for d in dirs:
        if d and os.path.isdir(d):
            for f in glob.glob(os.path.join(d, '[Aa]rial*.ttf')):
                font_manager.fontManager.addfont(f)
    names = {f.name for f in font_manager.fontManager.ttflist}
    family = 'Arial'
    if 'Arial' not in names:
        family = next((n for n in ('Liberation Sans', 'Helvetica', 'DejaVu Sans') if n in names), 'sans-serif')
        print('WARNING: Arial not found (install mscorefonts or set COMET_FONT_DIR); using %s' % family,
              file=sys.stderr)
    matplotlib.rcParams.update({
        'font.family': family, 'font.size': 9, 'pdf.fonttype': 42, 'ps.fonttype': 42,
        'axes.spines.top': False, 'axes.spines.right': False, 'axes.linewidth': 0.8,
        'xtick.major.width': 0.8, 'ytick.major.width': 0.8, 'legend.frameon': False,
        'savefig.dpi': 300, 'savefig.bbox': 'tight'})
    import matplotlib.pyplot as plt
    return plt


def read_tsv(path):
    with open(path) as f:
        rows = [l.rstrip('\n').split('\t') for l in f if l.strip()]
    return [dict(zip(rows[0], r)) for r in rows[1:]]


def save(fig, prefix, name):
    for ext in ('pdf', 'png'):
        fig.savefig('%s.%s.%s' % (prefix, name, ext))


def plot_composition(plt, samples, comp, cell_types, prefix):
    n = len(samples)
    fig, ax = plt.subplots(figsize=(6.0, 0.45 * n + 1.1))
    y = list(range(n))[::-1]                 # first sample on top
    left = [0.0] * n
    for k, ct in enumerate(cell_types):
        w = [comp[s].get(ct, 0.0) for s in samples]
        ax.barh(y, w, left=left, height=0.7, color=PALETTE[k % len(PALETTE)], edgecolor='white',
                linewidth=0.6, label=ct)
        for yi, l, wi in zip(y, left, w):
            if wi >= 0.06:
                ax.text(l + wi / 2, yi, '%.0f%%' % (100 * wi), ha='center', va='center', fontsize=7.5,
                        color='white')
        left = [l + wi for l, wi in zip(left, w)]
    ax.set_yticks(y)
    ax.set_yticklabels(samples)
    ax.set_xlim(0, 1)
    ax.set_xticks([0, 0.25, 0.5, 0.75, 1])
    ax.set_xticklabels(['0', '25', '50', '75', '100'])
    ax.set_xlabel('Cell-type composition (% of assigned reads)')
    ax.spines['left'].set_visible(False)
    ax.tick_params(axis='y', length=0)
    ax.legend(loc='lower left', bbox_to_anchor=(0, 1.0), ncol=min(len(cell_types), 4), fontsize=8,
              handlelength=1.0, columnspacing=1.2)
    save(fig, prefix, 'composition')
    plt.close(fig)


def plot_methylation(plt, samples, meth, cell_types, prefix):
    n = len(samples)
    fig, ax = plt.subplots(figsize=(5.0, 0.45 * n + 1.1))
    y = list(range(n))[::-1]
    vals = [meth[s][ct] for s in samples for ct in cell_types if ct in meth[s]]
    for yi, s in zip(y, samples):
        v = [meth[s][ct] for ct in cell_types if ct in meth[s]]
        if v:
            ax.plot([min(v), max(v)], [yi, yi], color='#BBBBBB', lw=1.2, zorder=1)
    for k, ct in enumerate(cell_types):
        pts = [(meth[s][ct], yi) for yi, s in zip(y, samples) if ct in meth[s]]
        ax.scatter([p[0] for p in pts], [p[1] for p in pts], s=46, color=PALETTE[k % len(PALETTE)],
                   edgecolor='white', linewidth=0.6, zorder=2, label=ct)
    lo, hi = min(vals), max(vals)
    pad = max(0.02, 0.08 * (hi - lo))
    ax.set_xlim(max(0, lo - pad), min(1, hi + pad))
    ax.set_ylim(-0.6, n - 0.4)
    ax.set_yticks(y)
    ax.set_yticklabels(samples)
    ax.set_xlabel('Genome-wide CpG methylation')
    ax.spines['left'].set_visible(False)
    ax.tick_params(axis='y', length=0)
    ax.grid(axis='x', color='#EEEEEE', lw=0.6)
    ax.set_axisbelow(True)
    ax.legend(loc='lower left', bbox_to_anchor=(0, 1.0), ncol=min(len(cell_types), 4), fontsize=8,
              handletextpad=0.2, columnspacing=1.2)
    save(fig, prefix, 'cell_type_methylation')
    plt.close(fig)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--sample', action='append', required=True, help='NAME=QUANT_DIR')
    p.add_argument('--cell-types', default='')
    p.add_argument('--out-prefix', required=True)
    p.add_argument('--meth-samples', default='', help='samples shown in the methylation plot [all]')
    p.add_argument('--tables', action='store_true', help='also write combined long-format tables')
    p.add_argument('--no-plot', action='store_true', help='tables only')
    a = p.parse_args()

    cell_types = [c for c in a.cell_types.split(',') if c]
    samples, comp, meth, comp_rows, meth_rows = [], {}, {}, [], []
    for x in a.sample:
        name, d = x.split('=', 1)
        samples.append(name)
        c = read_tsv(os.path.join(d, name + '.composition.tsv'))
        g = read_tsv(os.path.join(d, name + '.gw_methylation.tsv'))
        if not cell_types:
            cell_types = [r['group'] for r in c if r['frac_single'] != 'NA' and r['group'] != 'Unassigned']
        comp[name] = {r['group']: float(r['frac_single']) for r in c
                      if r['group'] in cell_types and r['frac_single'] != 'NA'}
        meth[name] = {r['group']: float(r['methylation']) for r in g
                      if r['group'] in cell_types and r['methylation'] not in ('NA', '')}
        comp_rows += [[name] + [r[k] for k in c[0].keys()] for r in c] if c else []
        meth_rows += [[name] + [r[k] for k in g[0].keys()] for r in g] if g else []
        comp_hdr, meth_hdr = list(c[0].keys()), list(g[0].keys())

    if a.tables:
        for rows, hdr, name in ((comp_rows, comp_hdr, 'composition'), (meth_rows, meth_hdr, 'gw_methylation')):
            with open('%s.%s.tsv' % (a.out_prefix, name), 'w') as f:
                f.write('\t'.join(['sample'] + hdr) + '\n')
                for r in rows:
                    f.write('\t'.join(r) + '\n')

    if a.no_plot:
        return
    try:
        plt = setup_matplotlib()
    except ImportError:
        print('WARNING: matplotlib not installed; plots skipped', file=sys.stderr)
        return
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        plot_composition(plt, samples, comp, cell_types, a.out_prefix)
        ms = [s for s in a.meth_samples.split(',') if s] or samples
        unknown = [s for s in ms if s not in meth]
        if unknown:
            sys.exit('ERROR: --meth-samples not among the samples: ' + ','.join(unknown))
        if any(meth[s] for s in ms):
            plot_methylation(plt, ms, meth, cell_types, a.out_prefix)


if __name__ == '__main__':
    main()
